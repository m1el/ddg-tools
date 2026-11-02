//! Black-box comparison with the authors' program (`SPEC.md` §14.7, `COMPARISON.md`).
//!
//! ```text
//! compare_reference run SCENE --mode compat|default|noedge|compat-noedge --steps N --out PREFIX
//!                     [--threads T] [--obj-every K] [--check-every K] [--no-remesh]
//! compare_reference stats [--alpha A --beta B] [--exact-max F] OBJ...
//! compare_reference swept A.obj B.obj
//! ```
//!
//! * `run` steps our flow on a scene file, like `run_scene`, with the mode
//!   `compat` = `reference_compat` and no CCD (the reference has none), `default` = our defaults
//!   (CCD on unless the exponents are subcritical), `noedge` = defaults with edge splits and
//!   collapses off (the non-orientable variant of the reference does not split or collapse),
//!   `compat-noedge` = both. It writes `PREFIX.csv` (one row per step: cumulative seconds, step ms, `F`,
//!   |V|, |F|, area, volume, whether the step end self-intersects) and `PREFIX.obj` (final mesh,
//!   plus `PREFIX_NNNNN.obj` every `--obj-every` steps).
//! * `stats` prints one CSV line per OBJ: |V|, |F|, area, volume, bounding box, principal
//!   standard deviations of the area-weighted point cloud (a rotation-invariant shape measure),
//!   the minimal distance of non-adjacent elements, whether it self-intersects, and the tangent-point
//!   energy (exact up to `--exact-max` faces, Barnes–Hut with θ = 0.25 above).
//! * `swept` runs the exact-time CCD (gap 1e-6 ℓ̄) on the linear motion from mesh A to mesh B (same
//!   connectivity) and says whether any pair of elements meets in between (tunnelling).

use std::io::Write;
use std::time::Instant;

use ddg_repulsive_surfaces::ccd::{self, CcdOptions, first_contact};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::{face_quantities, mean_edge_length};
use ddg_repulsive_surfaces::step::scene_options;
use ddg_repulsive_surfaces::{
    BhOptions, EnergyEval, Exponents, SceneSpec, Session, Status, Surface, energy, read_obj, write_obj,
};

fn area_volume(s: &Surface, p: &[[f64; 3]]) -> (f64, f64) {
    let (mut a, mut v) = (0.0, 0.0);
    for f in s.faces() {
        let t = f.map(|i| p[i as usize]);
        a += face_quantities(t).0;
        let c = [
            t[1][1] * t[2][2] - t[1][2] * t[2][1],
            t[1][2] * t[2][0] - t[1][0] * t[2][2],
            t[1][0] * t[2][1] - t[1][1] * t[2][0],
        ];
        v += (t[0][0] * c[0] + t[0][1] * c[1] + t[0][2] * c[2]) / 6.0;
    }
    (a, v)
}

/// Principal standard deviations of the area-weighted surface point cloud (face centroids).
fn principal_sd(s: &Surface, p: &[[f64; 3]]) -> [f64; 3] {
    let (mut w, mut m) = (0.0, [0.0; 3]);
    let mut items = Vec::new();
    for f in s.faces() {
        let (a, _, x) = face_quantities(f.map(|i| p[i as usize]));
        items.push((a, x));
        w += a;
        for k in 0..3 {
            m[k] += a * x[k];
        }
    }
    for k in 0..3 {
        m[k] /= w;
    }
    let mut c = [[0.0; 3]; 3];
    for (a, x) in items {
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] += a * (x[i] - m[i]) * (x[j] - m[j]) / w;
            }
        }
    }
    // Jacobi eigenvalues of the symmetric 3x3 matrix.
    for _ in 0..60 {
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            if c[i][j].abs() < 1e-300 {
                continue;
            }
            let th = 0.5 * (c[j][j] - c[i][i]) / c[i][j];
            let t = th.signum() / (th.abs() + (th * th + 1.0).sqrt());
            let (cs, sn) = (1.0 / (t * t + 1.0).sqrt(), t / (t * t + 1.0).sqrt());
            for k in 0..3 {
                let (a, b) = (c[k][i], c[k][j]);
                c[k][i] = cs * a - sn * b;
                c[k][j] = sn * a + cs * b;
            }
            for k in 0..3 {
                let (a, b) = (c[i][k], c[j][k]);
                c[i][k] = cs * a - sn * b;
                c[j][k] = sn * a + cs * b;
            }
        }
    }
    let mut e = [c[0][0].max(0.0).sqrt(), c[1][1].max(0.0).sqrt(), c[2][2].max(0.0).sqrt()];
    e.sort_by(|a, b| b.total_cmp(a));
    e
}

fn bbox(p: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for q in p {
        for k in 0..3 {
            lo[k] = lo[k].min(q[k]);
            hi[k] = hi[k].max(q[k]);
        }
    }
    (lo, hi)
}

fn arg<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).map(String::as_str)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let r = match args.first().map(String::as_str) {
        Some("run") => run(&args[1..]),
        Some("stats") => stats(&args[1..]),
        Some("swept") => swept(&args[1..]),
        _ => Err("usage: compare_reference run|stats|swept ... (see the file header)".into()),
    };
    if let Err(e) = r {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run(args: &[String]) -> Result<(), String> {
    let scene = args.first().ok_or("scene file")?;
    let mode = arg(args, "--mode").unwrap_or("default");
    let steps: usize = arg(args, "--steps").map(|s| s.parse().unwrap()).ok_or("--steps")?;
    let prefix = arg(args, "--out").ok_or("--out PREFIX")?;
    let threads: Option<usize> = arg(args, "--threads").map(|s| s.parse().unwrap());
    let obj_every: usize = arg(args, "--obj-every").map_or(0, |s| s.parse().unwrap());
    let check_every: usize = arg(args, "--check-every").map_or(1, |s| s.parse().unwrap());
    let spec = SceneSpec::load(scene).map_err(|e| e.to_string())?;
    let subcritical = spec.exponents.subcritical();
    let mut o = scene_options(&spec).map_err(|e| e.to_string())?;
    o.max_iterations = steps;
    o.threads = threads;
    o.stop = ddg_repulsive_surfaces::StopRule::Off;
    let (compat, noedge) = match mode {
        "compat" => (true, false),
        "default" => (false, false),
        "noedge" => (false, true),
        "compat-noedge" => (true, true),
        m => return Err(format!("unknown mode {m}")),
    };
    if compat {
        o.reference_compat = true;
        o.line_search.ccd = None;
    }
    if arg(args, "--ccd") == Some("off") {
        o.line_search.ccd = None;
    }
    if args.iter().any(|a| a == "--no-remesh") {
        o.remesh = None;
    }
    if noedge && let Some(r) = o.remesh.as_mut() {
        r.enabled_ops.split = false;
        r.enabled_ops.collapse = false;
    }
    let mut session = Session::new(spec.surface, spec.positions, o).map_err(|e| e.to_string())?;
    let mut csv = std::io::BufWriter::new(std::fs::File::create(format!("{prefix}.csv")).map_err(|e| e.to_string())?);
    writeln!(csv, "step,cum_s,step_ms,gmres,F,V,Fc,area,volume,intersects,tau,tau_ccd,halvings,splits,collapses,flips").unwrap();
    let (a0, v0) = area_volume(session.surface(), session.positions());
    let e0 = session.energy();
    writeln!(csv, "0,0,0,0,{e0:.9e},{},{},{a0:.9e},{v0:.9e},0,0,,0,0,0,0", session.surface().n_vertices(), session.surface().n_faces()).unwrap();
    if obj_every > 0 {
        write_obj(format!("{prefix}_{:05}.obj", 0), session.surface(), session.positions()).map_err(|e| e.to_string())?;
    }
    let t0 = Instant::now();
    let (mut inter_steps, mut stuck, mut ms_sum, mut gm_sum) = (0usize, 0usize, 0.0, 0.0);
    let mut last_status = Status::Running;
    while session.iteration() < steps {
        let rep = session.step().map_err(|e| e.to_string())?;
        let (s, p) = (session.surface(), session.positions());
        let (a, v) = area_volume(s, p);
        let inter = if !subcritical && check_every > 0 && rep.iteration % check_every == 0 {
            let x = ccd::intersects(s, p).map_err(|e| e.to_string())?;
            inter_steps += x as usize;
            x as i32 as f64
        } else {
            f64::NAN
        };
        stuck += (rep.status == ddg_repulsive_surfaces::StepStatus::Stuck) as usize;
        let ms = rep.seconds.total * 1e3;
        ms_sum += ms;
        let gm = rep.linear.total_iterations() as f64;
        gm_sum += gm;
        let (sp, co, fl) = rep.remesh.as_ref().map_or((0, 0, 0), |m| (m.splits, m.collapses, m.flips));
        writeln!(
            csv,
            "{},{:.3},{ms:.1},{gm},{:.9e},{},{},{a:.9e},{v:.9e},{inter},{:.4e},{},{},{sp},{co},{fl}",
            rep.iteration,
            t0.elapsed().as_secs_f64(),
            rep.energy_after,
            rep.vertices,
            rep.faces,
            rep.tau,
            rep.tau_ccd.map_or(String::new(), |t| format!("{t:.4e}")),
            rep.halvings
        )
        .unwrap();
        csv.flush().unwrap();
        if obj_every > 0 && rep.iteration % obj_every == 0 {
            write_obj(format!("{prefix}_{:05}.obj", rep.iteration), s, p).map_err(|e| e.to_string())?;
        }
        last_status = session.status();
        if matches!(last_status, Status::Error(_)) {
            break;
        }
    }
    write_obj(format!("{prefix}.obj"), session.surface(), session.positions()).map_err(|e| e.to_string())?;
    let n = session.iteration().max(1) as f64;
    println!(
        "summary scene={scene} mode={mode} steps={} status={last_status:?} F0={e0:.6e} F={:.6e} V={} Fc={} ms_per_step={:.1} gmres_per_step={:.1} \
         intersecting_step_ends={inter_steps} stuck={stuck} seconds={:.1}",
        session.iteration(),
        session.energy(),
        session.surface().n_vertices(),
        session.surface().n_faces(),
        ms_sum / n,
        gm_sum / n,
        t0.elapsed().as_secs_f64()
    );
    let c = session.constraints();
    println!("final constraints (targets at the end): {c:?}");
    Ok(())
}

fn stats(args: &[String]) -> Result<(), String> {
    let alpha: f64 = arg(args, "--alpha").map_or(6.0, |s| s.parse().unwrap());
    let beta: f64 = arg(args, "--beta").map_or(12.0, |s| s.parse().unwrap());
    let exact_max: usize = arg(args, "--exact-max").map_or(16000, |s| s.parse().unwrap());
    let skip_inter = args.iter().any(|a| a == "--no-intersect");
    let e = Exponents::new(alpha, beta).map_err(|e| e.to_string())?;
    let mut files = Vec::new();
    let mut i = 0;
    while i < args.len() {
        if args[i].starts_with("--") && args[i] != "--no-intersect" {
            i += 2;
        } else if args[i] == "--no-intersect" {
            i += 1;
        } else {
            files.push(args[i].clone());
            i += 1;
        }
    }
    println!("file,V,F,area,volume,bbox_x,bbox_y,bbox_z,sd1,sd2,sd3,min_dist,intersects,energy,energy_kind");
    for f in files {
        let (s, p) = read_obj(&f).map_err(|e| format!("{f}: {e}"))?;
        let (a, v) = area_volume(&s, &p);
        let (lo, hi) = bbox(&p);
        let sd = principal_sd(&s, &p);
        let (md, inter) = if skip_inter {
            (f64::NAN, f64::NAN)
        } else {
            (
                ccd::min_distance(&s, &p).map_err(|e| e.to_string())? / mean_edge_length(&s, &p),
                ccd::intersects(&s, &p).map_err(|e| e.to_string())? as i32 as f64,
            )
        };
        let (eval, kind) = if s.n_faces() <= exact_max {
            (EnergyEval::Exact, "exact")
        } else {
            (EnergyEval::BarnesHut(BhOptions { theta: 0.25, leaf_size: 8 }), "bh0.25")
        };
        let en = energy(&s, &p, e, &eval, &[], None).map_err(|e| e.to_string())?.energy;
        println!(
            "{f},{},{},{a:.8e},{v:.8e},{:.5},{:.5},{:.5},{:.5},{:.5},{:.5},{md:.4},{inter},{en:.8e},{kind}",
            s.n_vertices(),
            s.n_faces(),
            hi[0] - lo[0],
            hi[1] - lo[1],
            hi[2] - lo[2],
            sd[0],
            sd[1],
            sd[2]
        );
    }
    Ok(())
}

fn swept(args: &[String]) -> Result<(), String> {
    let (s, a) = read_obj(args.first().ok_or("A.obj")?).map_err(|e| e.to_string())?;
    let (s2, b) = read_obj(args.get(1).ok_or("B.obj")?).map_err(|e| e.to_string())?;
    if s.faces() != s2.faces() {
        return Err("different connectivity".into());
    }
    let vel: Vec<[f64; 3]> = a.iter().zip(&b).map(|(x, y)| [y[0] - x[0], y[1] - x[1], y[2] - x[2]]).collect();
    let pool = make_pool(None).map_err(|e| e.to_string())?;
    let opts = CcdOptions { gap: 1e-6, reference_length: Some(mean_edge_length(&s, &a)), close_factor: 0.0, max_iterations: 1_000_000, ..CcdOptions::default() };
    let r = first_contact(&pool, &s, &a, &vel, 1.0, &[], &opts).map_err(|e| e.to_string())?;
    let maxmove = vel.iter().map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()).fold(0.0, f64::max) / mean_edge_length(&s, &a);
    println!("tau_max={:.6} contact={:?} max_vertex_move_in_mean_edges={maxmove:.3}", r.tau_max, r.contact);
    Ok(())
}
