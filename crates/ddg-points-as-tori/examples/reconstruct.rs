//! Oriented point cloud (binary little-endian PLY with x y z nx ny nz floats, or OBJ with `v` and
//! `vn` lines) -> Points as Tori -> marching cubes -> OBJ.
//!
//! ```text
//! cargo run --release -p ddg-points-as-tori --example reconstruct -- <cloud> <out-dir> [res,res,..] [threads] [refine] [bound]
//! ```
//!
//! The cloud is centred and scaled into `[-1, 1]³` like the reference's demo. The tori are fitted
//! once; marching cubes then runs at every resolution (default 64,128,256) in the box `[-1.1, 1.1]³`.
//! Prints timings, the narrow-band statistics and the topology check of every mesh.

use ddg_points_as_tori::mc::{Band, McOptions, MeshCheck, Resolution, marching_cubes};
use ddg_points_as_tori::{BuildOptions, Model, NormalPolicy, PointCloud, PointsAsTori, TorusField};

/// Positions and normals.
type Oriented = (Vec<[f64; 3]>, Vec<[f64; 3]>);

fn read_ply(bytes: &[u8]) -> Result<Oriented, String> {
    let end = bytes.windows(10).position(|w| w == b"end_header").ok_or("no end_header")?;
    let header = std::str::from_utf8(&bytes[..end]).map_err(|e| e.to_string())?;
    if !header.contains("binary_little_endian") {
        return Err("only binary_little_endian PLY is supported".into());
    }
    let n: usize = header
        .lines()
        .find_map(|l| l.strip_prefix("element vertex "))
        .ok_or("no vertex element")?
        .trim()
        .parse()
        .map_err(|e| format!("{e}"))?;
    let props = header.lines().filter(|l| l.starts_with("property float ")).count();
    if props != 6 {
        return Err(format!("expected 6 float properties (x y z nx ny nz), found {props}"));
    }
    let start = bytes[end..].iter().position(|b| *b == b'\n').ok_or("bad header")? + end + 1;
    let data = &bytes[start..];
    if data.len() < n * 24 {
        return Err("truncated vertex data".into());
    }
    let f = |i: usize| f32::from_le_bytes([data[4 * i], data[4 * i + 1], data[4 * i + 2], data[4 * i + 3]]) as f64;
    let p = (0..n).map(|i| [f(6 * i), f(6 * i + 1), f(6 * i + 2)]).collect();
    let nr = (0..n).map(|i| [f(6 * i + 3), f(6 * i + 4), f(6 * i + 5)]).collect();
    Ok((p, nr))
}

fn read_obj(text: &str) -> Result<Oriented, String> {
    let (mut p, mut n) = (Vec::new(), Vec::new());
    for l in text.lines() {
        let mut it = l.split_whitespace();
        let tag = it.next();
        if tag == Some("v") || tag == Some("vn") {
            let v: Vec<f64> = it.take(3).map(|s| s.parse::<f64>().map_err(|e| format!("{l}: {e}"))).collect::<Result<_, _>>()?;
            if v.len() < 3 {
                return Err(format!("short line: {l}"));
            }
            if tag == Some("v") { p.push([v[0], v[1], v[2]]) } else { n.push([v[0], v[1], v[2]]) }
        }
    }
    if p.len() != n.len() || p.is_empty() {
        return Err(format!("{} v lines but {} vn lines (this example needs one normal per vertex)", p.len(), n.len()));
    }
    Ok((p, n))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        return Err("usage: reconstruct <cloud.ply|cloud.obj> <out-dir> [res,res,..] [threads] [refine_steps] [bound|0=dense]".into());
    }
    let path = &args[1];
    let out_dir = std::path::Path::new(&args[2]);
    let resolutions: Vec<usize> = match args.get(3) {
        Some(s) => s.split(',').map(|x| x.parse()).collect::<Result<_, _>>()?,
        None => vec![64, 128, 256],
    };
    let threads: Option<usize> = args.get(4).map(|s| s.parse()).transpose()?;
    let refine: usize = args.get(5).map(|s| s.parse()).transpose()?.unwrap_or(0);
    // Lipschitz bound of the narrow band (0: dense grid)
    let bound: f64 = args.get(6).map(|s| s.parse()).transpose()?.unwrap_or(2.0);
    std::fs::create_dir_all(out_dir)?;
    let name = std::path::Path::new(path).file_stem().and_then(|s| s.to_str()).unwrap_or("cloud").to_string();

    let (mut p, n) = if path.ends_with(".ply") { read_ply(&std::fs::read(path)?)? } else { read_obj(&std::fs::read_to_string(path)?)? };
    let k = p.len() as f64;
    let mean = [0, 1, 2].map(|a| p.iter().map(|x| x[a]).sum::<f64>() / k);
    let scale = p.iter().flat_map(|x| (0..3).map(move |a| (x[a] - mean[a]).abs())).fold(0.0, f64::max);
    for x in p.iter_mut() {
        for a in 0..3 {
            x[a] = (x[a] - mean[a]) / scale;
        }
    }
    let (cloud, crep) = PointCloud::new(p, n, NormalPolicy::Normalize)?;
    println!("{name}: {} points, {} normals renormalized", crep.n, crep.normalized);

    let model = Model::reference()?;
    let mut bo = BuildOptions::default().with_threads(threads);
    bo.eval = TorusField::paper_options(&cloud, threads)?.0;
    let (field, brep) = PointsAsTori::build(&cloud, &model, &bo)?;
    println!("build: {:.2} s (knn {:.2}, network {:.2}, fit {:.2})", brep.seconds, brep.knn_seconds, brep.predict.forward_seconds, brep.fit.seconds);

    let load = std::fs::read_to_string("/proc/loadavg").unwrap_or_default();
    println!("loadavg before marching cubes: {}", load.trim());
    for &res in &resolutions {
        let opts = McOptions {
            bbox: Some([[-1.1; 3], [1.1; 3]]),
            resolution: Resolution::Longest(res),
            refine_steps: refine,
            band: if bound == 0.0 { Band::Dense } else { Band::Lipschitz { bound, block: 8 } },
            threads,
            ..Default::default()
        };
        let t = std::time::Instant::now();
        let (mesh, rep) = marching_cubes(&field, &opts)?;
        let secs = t.elapsed().as_secs_f64();
        let topo = match &rep.mesh {
            MeshCheck::Manifold { boundary_edges, euler_characteristic, components } => {
                format!("manifold, boundary edges {boundary_edges}, chi {euler_characteristic}, components {components}")
            }
            other => format!("{other:?}"),
        };
        println!(
            "res {res:4}: {secs:7.3} s (grid evals {:.3} s, cells {:.3} s, vertices {:.3} s) | {} V {} F | grid values {} of {} ({:.1}%), blocks {}/{} | max slope {:.3} (bound {:?}) | ambiguous faces {} | {topo}",
            rep.seconds[1],
            rep.seconds[2],
            rep.seconds[3],
            rep.vertices,
            rep.triangles,
            rep.grid_values,
            rep.grid_points,
            100.0 * rep.grid_values as f64 / rep.grid_points as f64,
            rep.blocks_active,
            rep.blocks,
            rep.max_edge_slope,
            rep.lipschitz_bound,
            rep.ambiguous_cell_faces,
        );
        let corners: Vec<u32> = mesh.faces.iter().flatten().copied().collect();
        let starts: Vec<u32> = (0..mesh.faces.len() as u32).map(|i| 3 * i).collect();
        let file = out_dir.join(format!("{name}-{res}.obj"));
        ddg_core::io::write_obj(&file, &mesh.positions, &corners, &starts, None)?;
    }
    Ok(())
}
