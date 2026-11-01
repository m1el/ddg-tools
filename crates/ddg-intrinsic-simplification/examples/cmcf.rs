//! Conformalized mean curvature flow (Kazhdan et al. 2012) accelerated by the intrinsic
//! multigrid hierarchy (paper Sec. 8.7, Fig. 30; SPEC §9.3), against a direct solve per step.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example cmcf -- MESH.obj \
//!     [--steps N] [--h H] [--out DIR] [--write-every K] [--sweeps S] [--tol T] [--no-compare] [--no-rescale] [--subdivide K]
//! ```
//!
//! After every step the surface is re-centred and rescaled to the input's area (`--no-rescale`
//! to disable: the flow then shrinks to a point in finite time). The input is translated to its centroid and scaled to unit surface area, so `H` (default
//! `1e-3`) is an area fraction. The hierarchy of the input is built once (reported). Each step
//! runs (a) the direct solver on the direct trajectory and (b) PCG-multigrid on its own
//! trajectory; also (c) one multigrid step from the direct state (the per-step difference). The
//! flowed meshes (multigrid trajectory) are written as `DIR/step_XXXX.obj` (default
//! `runs_tmp/cmcf_<mesh>/`), every `K`-th step and the last. Timings are single-thread; run on a
//! quiet machine.

use std::time::Instant;

use ddg_core::TriMesh;
use ddg_core::io::{read_obj, write_obj};
use ddg_intrinsic_simplification::cmcf::{CmcfFlow, CmcfMultigridOptions};
use ddg_intrinsic_simplification::hierarchy::{HierarchyOptions, build_hierarchy};
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions};

fn arg<T: std::str::FromStr>(args: &[String], name: &str) -> Option<T> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|x| x.parse().ok())
}

fn die(msg: String) -> ! {
    eprintln!("{msg}");
    std::process::exit(1);
}

fn rel_diff(x: &[[f64; 3]], y: &[[f64; 3]]) -> f64 {
    let mut d = 0.0;
    let mut s = 0.0;
    for (a, b) in x.iter().zip(y) {
        for c in 0..3 {
            d += (a[c] - b[c]).powi(2);
            s += b[c].powi(2);
        }
    }
    (d / s).sqrt()
}

fn area(faces: &[[u32; 3]], p: &[[f64; 3]]) -> f64 {
    faces
        .iter()
        .map(|f| {
            let [a, b, c] = f.map(|i| p[i as usize]);
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let w = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt()
        })
        .sum()
}

/// One 4-to-1 midpoint split (the numbering of `ddg_intrinsic_apps::fem::regular_subdivision`).
fn subdivide(faces: &[[u32; 3]], pos: &[[f64; 3]]) -> (Vec<[u32; 3]>, Vec<[f64; 3]>) {
    let mesh = TriMesh::new(pos.len(), faces.to_vec()).unwrap_or_else(|e| die(e.to_string()));
    let n = pos.len() as u32;
    let mut p = pos.to_vec();
    for &[a, b] in mesh.edges() {
        let (pa, pb) = (pos[a as usize], pos[b as usize]);
        p.push([0, 1, 2].map(|k| 0.5 * (pa[k] + pb[k])));
    }
    let mut out = Vec::new();
    for (f, &[a, b, c]) in faces.iter().enumerate() {
        let [mab, mbc, mca] = [0, 1, 2].map(|k| n + mesh.edge(3 * f + k) as u32);
        out.extend([[a, mab, mca], [mab, b, mbc], [mca, mbc, c], [mab, mbc, mca]]);
    }
    (out, p)
}

fn write(dir: &str, step: usize, pos: &[[f64; 3]], faces: &[[u32; 3]]) {
    let corners: Vec<u32> = faces.iter().flatten().copied().collect();
    let starts: Vec<u32> = (0..faces.len() as u32).map(|f| 3 * f).collect();
    let path = format!("{dir}/step_{step:04}.obj");
    write_obj(&path, pos, &corners, &starts, None).unwrap_or_else(|e| die(format!("{path}: {e}")));
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        die("usage: cmcf MESH.obj [--steps N] [--h H] [--out DIR] [--write-every K]".into());
    };
    let steps: usize = arg(&args, "--steps").unwrap_or(10);
    let h: f64 = arg(&args, "--h").unwrap_or(1e-3);
    let every: usize = arg(&args, "--write-every").unwrap_or(1);
    let sweeps: usize = arg(&args, "--sweeps").unwrap_or(2);
    let tol: f64 = arg(&args, "--tol").unwrap_or(1e-10);
    let rescale = !args.iter().any(|a| a == "--no-rescale");
    let compare = !args.iter().any(|a| a == "--no-compare");
    let name = std::path::Path::new(path)
        .file_stem()
        .map_or("mesh".to_string(), |s| s.to_string_lossy().to_string());
    let out: String = arg(&args, "--out").unwrap_or(format!("runs_tmp/cmcf_{name}"));
    std::fs::create_dir_all(&out).unwrap_or_else(|e| die(format!("{out}: {e}")));

    let obj = read_obj(path).unwrap_or_else(|e| die(format!("{path}: {e}")));
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let mut pos = obj.positions.clone();
    for _ in 0..arg::<usize>(&args, "--subdivide").unwrap_or(0) {
        (faces, pos) = subdivide(&faces, &pos);
    }
    let n = pos.len() as f64;
    let centroid = [0, 1, 2].map(|c| pos.iter().map(|p| p[c]).sum::<f64>() / n);
    let scale = 1.0 / area(&faces, &pos).sqrt();
    for p in &mut pos {
        for c in 0..3 {
            p[c] = (p[c] - centroid[c]) * scale;
        }
    }
    let mesh =
        TriMesh::new(pos.len(), faces.clone()).unwrap_or_else(|e| die(format!("{path}: {e}")));

    let t = Instant::now();
    let coarse = CoarseMesh::from_positions(&pos, &faces, &[], CoarseMeshOptions::default())
        .unwrap_or_else(|e| die(format!("{path}: {e}")));
    let hierarchy = build_hierarchy(&coarse, &HierarchyOptions::default())
        .unwrap_or_else(|e| die(format!("{path}: {e}")));
    let build_seconds = t.elapsed().as_secs_f64();
    println!(
        "{name}: {} vertices, hierarchy {:?} built once in {build_seconds:.3} s, h = {h}",
        pos.len(),
        hierarchy.report().counts
    );
    let mut opts = CmcfMultigridOptions::default();
    opts.multigrid.pre_sweeps = sweeps;
    opts.multigrid.post_sweeps = sweeps;
    opts.multigrid.tolerance = tol;

    let mut direct = CmcfFlow::new(&mesh, &pos, h).unwrap_or_else(|e| die(format!("{name}: {e}")));
    let mut mg = direct.clone();
    write(&out, 0, &pos, &faces);
    println!(
        "step | direct s (factor+solve) | mg s (assemble+galerkin+pcg) | speedup | speedup w/o assembly | pcg iters | one-step diff | trajectory diff"
    );
    let (mut sum_d, mut sum_m) = (0.0, 0.0);
    for s in 1..=steps {
        let mut probe = direct.clone();
        let rd = direct
            .step_direct()
            .unwrap_or_else(|e| die(format!("step {s}: {e}")));
        if rescale {
            direct
                .normalize()
                .unwrap_or_else(|e| die(format!("step {s}: {e}")));
        }
        let rm = mg
            .step_multigrid(&hierarchy, &opts)
            .unwrap_or_else(|e| die(format!("step {s}: {e}")));
        if rescale {
            mg.normalize()
                .unwrap_or_else(|e| die(format!("step {s}: {e}")));
        }
        let one_step = if compare {
            probe
                .step_multigrid(&hierarchy, &opts)
                .unwrap_or_else(|e| die(format!("step {s}: {e}")));
            if rescale {
                probe
                    .normalize()
                    .unwrap_or_else(|e| die(format!("step {s}: {e}")));
            }
            rel_diff(probe.positions(), direct.positions())
        } else {
            f64::NAN
        };
        let traj = rel_diff(mg.positions(), direct.positions());
        sum_d += rd.total_seconds();
        sum_m += rm.total_seconds();
        println!(
            "{s:4} | {:.4} ({:.4}+{:.4}) | {:.4} ({:.4}+{:.4}+{:.4}) | {:.2}x | {:.2}x | {:?}{} | {one_step:.2e} | {traj:.2e}",
            rd.total_seconds(),
            rd.setup_seconds,
            rd.solve_seconds,
            rm.total_seconds(),
            rm.assemble_seconds,
            rm.setup_seconds,
            rm.solve_seconds,
            rd.total_seconds() / rm.total_seconds(),
            (rd.setup_seconds + rd.solve_seconds) / (rm.setup_seconds + rm.solve_seconds),
            rm.iterations,
            if rm.converged { "" } else { " NOT CONVERGED" },
        );
        if s % every == 0 || s == steps {
            write(&out, s, mg.positions(), &faces);
        }
    }
    println!(
        "total over {steps} steps: direct {sum_d:.3} s, multigrid {sum_m:.3} s (hierarchy {build_seconds:.3} s once), speedup {:.2}x; area of the last mesh {:.4} of the first",
        sum_d / sum_m,
        area(&faces, mg.positions())
    );
}
