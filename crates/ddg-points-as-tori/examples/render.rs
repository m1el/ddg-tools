//! Sphere-traced rendering of a Points-as-Tori field (Fig. 12 analogues, SPEC.md §8.5).
//!
//! ```text
//! cargo run --release -p ddg-points-as-tori --example render -- <cloud> <out.png> [key=value ...]
//! ```
//!
//! `<cloud>` is a binary little-endian PLY with x y z nx ny nz floats, an OBJ with `v`/`vn` lines,
//! an OBJ/triangle mesh without normals (sampled at `points=2048`), or `analytic:<seed>` (the
//! analytic set of SPEC §11, 1024 points). The cloud is centred and scaled into `[-1, 1]³`.
//!
//! Keys (defaults in brackets): `res` [512] image side (`resx`/`resy` override), `threads` [8],
//! `level` [0], `offset` (an offset applied to the field), `op` (`union|intersect|subtract`) with
//! `other=<cloud>` and `other_shift=x,y,z`, `cut=<axis>` with `cut_offset` [0], `cut_rot_x`,
//! `cut_rot_y` (degrees), `az` and `el` (degrees; orbit camera instead of the demo camera),
//! `up=y|z` [y], `demo=1` (the demo's constants), `sky=forest|white` [forest], `alpha=1`
//! (transparent background PNG), `safety`, `steps`, `device=cpu|cuda` [cpu], `points` [2048],
//! `seed` [7] (mesh sampling).

use ddg_points_as_tori::csg::Csg;
use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::generate::sampling::sample_mesh_uniform;
use ddg_points_as_tori::render::{Camera, CutPlane, Environment, Tracer, TraceOptions};
use ddg_points_as_tori::{BuildOptions, Device, Model, NormalPolicy, PointCloud, PointsAsTori, TorusField};
use std::collections::HashMap;

type Oriented = (Vec<[f64; 3]>, Vec<[f64; 3]>);

fn read_ply(bytes: &[u8]) -> Result<Oriented, String> {
    let end = bytes.windows(10).position(|w| w == b"end_header").ok_or("no end_header")?;
    let header = std::str::from_utf8(&bytes[..end]).map_err(|e| e.to_string())?;
    if !header.contains("binary_little_endian") {
        return Err("only binary_little_endian PLY is supported".into());
    }
    let n: usize = header.lines().find_map(|l| l.strip_prefix("element vertex ")).ok_or("no vertex element")?.trim().parse().map_err(|e| format!("{e}"))?;
    if header.lines().filter(|l| l.starts_with("property float ")).count() != 6 {
        return Err("expected 6 float properties (x y z nx ny nz)".into());
    }
    let start = bytes[end..].iter().position(|b| *b == b'\n').ok_or("bad header")? + end + 1;
    let d = &bytes[start..];
    if d.len() < n * 24 {
        return Err("truncated vertex data".into());
    }
    let f = |i: usize| f32::from_le_bytes([d[4 * i], d[4 * i + 1], d[4 * i + 2], d[4 * i + 3]]) as f64;
    Ok(((0..n).map(|i| [f(6 * i), f(6 * i + 1), f(6 * i + 2)]).collect(), (0..n).map(|i| [f(6 * i + 3), f(6 * i + 4), f(6 * i + 5)]).collect()))
}

fn read_obj_normals(text: &str) -> Option<Oriented> {
    let (mut p, mut n) = (Vec::new(), Vec::new());
    for l in text.lines() {
        let mut it = l.split_whitespace();
        let tag = it.next();
        if tag == Some("v") || tag == Some("vn") {
            let v: Vec<f64> = it.take(3).filter_map(|s| s.parse().ok()).collect();
            if v.len() < 3 {
                return None;
            }
            if tag == Some("v") { p.push([v[0], v[1], v[2]]) } else { n.push([v[0], v[1], v[2]]) }
        }
    }
    (p.len() == n.len() && !p.is_empty()).then_some((p, n))
}

fn load(spec: &str, points: usize, seed: u64) -> Result<Oriented, Box<dyn std::error::Error>> {
    let (p, n) = if let Some(s) = spec.strip_prefix("analytic:") {
        ddg_points_as_tori::analytic::Instance::from_seed(s.parse()?)?.solid.sample_surface(1024, seed)?
    } else if spec.ends_with(".ply") {
        read_ply(&std::fs::read(spec)?)?
    } else {
        let text = std::fs::read_to_string(spec)?;
        match read_obj_normals(&text) {
            Some(x) => x,
            None => {
                let m = ddg_core::io::read_obj(spec)?;
                let mut faces = Vec::new();
                for (f, &s) in m.face_starts.iter().enumerate() {
                    let e = m.face_starts.get(f + 1).copied().unwrap_or(m.corner_verts.len() as u32);
                    for k in s + 1..e.saturating_sub(1) {
                        faces.push([m.corner_verts[s as usize], m.corner_verts[k as usize], m.corner_verts[k as usize + 1]]);
                    }
                }
                sample_mesh_uniform(&m.positions, &faces, points, &mut Rng::new(seed))?
            }
        }
    };
    let (mean, scale) = {
        let k = p.len() as f64;
        let mean = [0, 1, 2].map(|a| p.iter().map(|x| x[a]).sum::<f64>() / k);
        let scale = p.iter().flat_map(|x| (0..3).map(move |a| (x[a] - mean[a]).abs())).fold(0.0, f64::max);
        (mean, scale)
    };
    let p = p.iter().map(|x| [0, 1, 2].map(|a| (x[a] - mean[a]) / scale)).collect();
    Ok((p, n))
}

fn build(spec: &str, shift: [f64; 3], points: usize, seed: u64, threads: usize, device: Device, model: &Model) -> Result<TorusField, Box<dyn std::error::Error>> {
    let (p, n) = load(spec, points, seed)?;
    let p = p.iter().map(|x| [x[0] + shift[0], x[1] + shift[1], x[2] + shift[2]]).collect();
    let (cloud, _) = PointCloud::new(p, n, NormalPolicy::Normalize)?;
    let mut bo = BuildOptions::default().with_threads(Some(threads));
    bo.eval.device = device;
    Ok(PointsAsTori::build(&cloud, model, &bo)?.0)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        return Err("usage: render <cloud.ply|cloud.obj|analytic:SEED> <out.png> [key=value ...]".into());
    }
    let kv: HashMap<String, String> = args[3..].iter().filter_map(|a| a.split_once('=')).map(|(k, v)| (k.to_string(), v.to_string())).collect();
    let get = |k: &str| kv.get(k).map(|s| s.as_str());
    let num = |k: &str, d: f64| -> Result<f64, Box<dyn std::error::Error>> { Ok(get(k).map(|s| s.parse()).transpose()?.unwrap_or(d)) };
    let threads = num("threads", 8.0)? as usize;
    let points = num("points", 2048.0)? as usize;
    let seed = num("seed", 7.0)? as u64;
    let res = num("res", 512.0)? as usize;
    let (w, h) = (num("resx", res as f64)? as usize, num("resy", res as f64)? as usize);
    let device = match get("device").unwrap_or("cpu") {
        "cuda" => Device::Cuda { ordinal: 0 },
        _ => Device::Cpu,
    };
    let model = Model::reference()?;
    let t = std::time::Instant::now();
    let a = build(&args[1], [0.0; 3], points, seed, threads, device, &model)?;
    let torus_only = get("op").is_none() && get("offset").is_none();
    let (single, mut scene) = if torus_only { (Some(a), None) } else { (None, Some(Csg::field(a))) };
    if let Some(op) = get("op") {
        let other = get("other").ok_or("op needs other=<cloud>")?;
        let sh: Vec<f64> = get("other_shift").unwrap_or("0,0,0").split(',').map(|s| s.parse()).collect::<Result<_, _>>()?;
        if sh.len() != 3 {
            return Err("other_shift needs x,y,z".into());
        }
        let b = Csg::field(build(other, [sh[0], sh[1], sh[2]], points, seed, threads, device, &model)?);
        let a = scene.take().unwrap();
        scene = Some(match op {
            "union" => a.union(b),
            "intersect" => a.intersect(b),
            "subtract" => a.subtract(b),
            o => return Err(format!("op = {o} (union, intersect or subtract)").into()),
        });
    }
    if let Some(c) = get("offset") {
        let c: f64 = c.parse()?;
        scene = scene.map(|s| s.offset(c));
    }
    eprintln!("built in {:.2} s", t.elapsed().as_secs_f64());

    let mut o = if get("demo") == Some("1") { TraceOptions::demo() } else { TraceOptions::default() };
    o.level = num("level", 0.0)?;
    o.safety = num("safety", o.safety)?;
    o.max_steps = num("steps", o.max_steps as f64)? as u32;
    o.threads = Some(threads);
    o.device = device;
    if get("sky") == Some("white") {
        o.environment = Environment::Constant([1.0, 1.0, 1.0]);
    }
    let up_z = get("up") == Some("z");
    if up_z {
        o.ground.as_mut().unwrap().up = [0.0, 0.0, 1.0];
    }
    if let Some(ax) = get("cut") {
        o.cut_plane = Some(CutPlane {
            axis: ax.parse()?,
            offset: num("cut_offset", 0.0)?,
            rotation_x: num("cut_rot_x", 0.0)?.to_radians(),
            rotation_y: num("cut_rot_y", 0.0)?.to_radians(),
            ..Default::default()
        });
    }
    let camera = match (get("az"), get("el")) {
        (None, None) => Camera::demo([0.0; 3], 1.0, w, h)?,
        (az, el) => {
            let (az, el) = (az.map(|s| s.parse::<f64>()).transpose()?.unwrap_or(0.0), el.map(|s| s.parse::<f64>()).transpose()?.unwrap_or(15.0));
            Camera::orbit([0.0; 3], 10.0, az.to_radians(), el.to_radians(), if up_z { 2 } else { 1 }, 25.1, w, h)?
        }
    };
    let tracer = Tracer::new(o)?;
    // a plain TorusField can run on the GPU (device=cuda); CSG scenes render on the CPU
    let (img, rep) = match (&single, &scene) {
        (Some(f), _) => tracer.render_torus(f, &camera)?,
        (None, Some(sc)) => tracer.render(sc, &camera)?,
        _ => unreachable!(),
    };
    img.write_png(std::path::Path::new(&args[2]), get("alpha") == Some("1"))?;
    println!(
        "{}: {w}x{h}, {:.2} s on {} threads | rays {} (clipped {}), hits {}, misses {}, step limit {}, non-finite {} | steps mean {:.1} max {} | shadow rays {} hits {} limit {} | cut-plane px {}",
        args[2], rep.seconds, rep.threads, rep.rays, rep.clipped, rep.hits, rep.misses, rep.step_limit, rep.non_finite, rep.steps_mean, rep.steps_max, rep.shadow_rays, rep.shadow_hits, rep.shadow_step_limit, rep.cut_plane_pixels
    );
    Ok(())
}
