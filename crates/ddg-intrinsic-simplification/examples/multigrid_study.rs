//! Multigrid on the intrinsic hierarchy (SPEC §9.1-9.2): builds the hierarchy of a mesh (one run
//! with snapshots, S24, and the reference's fresh levels, R33), then solves `L + εM` and the
//! reference's Dirichlet example (`L z = M 1`, vertices 0, 1, 2 fixed to 0, random start) with
//! V-cycles and with PCG for several sweep counts, against a direct solve.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example multigrid_study -- MESH.obj \
//!     [--eps E] [--ratio r] [--n-min N] [--tol t]
//! ```
//!
//! Output: one line per (hierarchy variant, problem, sweeps, solver). Single thread.

use std::time::Instant;

use ddg_core::generators::SplitMix64;
use ddg_core::io::read_obj;
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::SparseMatrix;
use ddg_intrinsic_simplification::hierarchy::{HierarchyOptions, build_hierarchy};
use ddg_intrinsic_simplification::multigrid::MultigridOptions;
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions, SimplifyOptions, Target};

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

fn direct(a: &SparseMatrix, b: &[f64], known: &[u32], vals: &[f64]) -> (Vec<f64>, f64) {
    let t = Instant::now();
    let n = a.nrows();
    let mut is_k = vec![false; n];
    for &k in known {
        is_k[k as usize] = true;
    }
    let unk: Vec<u32> = (0..n as u32).filter(|&i| !is_k[i as usize]).collect();
    let auu = a.submatrix(&unk, &unk).expect("submatrix");
    let auk = a.submatrix(&unk, known).expect("submatrix");
    let mut rhs: Vec<f64> = unk.iter().map(|&i| b[i as usize]).collect();
    let mut tmp = vec![0.0; unk.len()];
    auk.mul_vec(vals, &mut tmp).expect("mul");
    for (r, ti) in rhs.iter_mut().zip(&tmp) {
        *r -= ti;
    }
    let f = Symbolic::analyze(&auu, Method::Llt)
        .expect("analyze")
        .factor(&auu, Parallelism::Sequential)
        .expect("factor");
    f.solve_in_place(&mut rhs, Parallelism::Sequential)
        .expect("solve");
    let mut x = vec![0.0; n];
    for (&i, &v) in unk.iter().zip(&rhs) {
        x[i as usize] = v;
    }
    for (&i, &v) in known.iter().zip(vals) {
        x[i as usize] = v;
    }
    (x, t.elapsed().as_secs_f64())
}

fn rel_diff(x: &[f64], y: &[f64]) -> f64 {
    let d: f64 = x.iter().zip(y).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();
    d / y.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        die("usage: multigrid_study MESH.obj [--eps E] [--ratio r] [--n-min N] [--tol t]".into());
    };
    let eps: f64 = arg(&args, "--eps").unwrap_or(1.0);
    let tol: f64 = arg(&args, "--tol").unwrap_or(1e-10);
    let defaults = HierarchyOptions::default();
    let ratio: f64 = arg(&args, "--ratio").unwrap_or(defaults.ratio);
    let n_min: usize = arg(&args, "--n-min").unwrap_or(defaults.n_min);
    let obj = read_obj(path).unwrap_or_else(|e| die(format!("{path}: {e}")));
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let name = std::path::Path::new(path)
        .file_stem()
        .map_or("mesh".to_string(), |s| s.to_string_lossy().to_string());
    let mesh = CoarseMesh::from_positions(&obj.positions, &faces, &[], CoarseMeshOptions::default())
        .unwrap_or_else(|e| die(format!("{path}: {e}")));
    for fresh in [false, true] {
        let variant = if fresh { "fresh" } else { "snapshots" };
        let options = HierarchyOptions {
            ratio,
            n_min,
            fresh_levels: fresh,
            simplify: if args.iter().any(|a| a == "--compat") {
                SimplifyOptions::reference_compat(0.0, Target::AllRemovable)
            } else {
                SimplifyOptions::default()
            },
            ..HierarchyOptions::default()
        };
        let t = Instant::now();
        let h = build_hierarchy(&mesh, &options).unwrap_or_else(|e| die(format!("{name}: {e}")));
        let rep = h.report();
        println!(
            "{name} {variant}: levels {:?} (targets {:?}), build {:.3} s (simplify {:.3}, operators {:.3}), P nnz {:?}",
            rep.counts,
            rep.targets,
            t.elapsed().as_secs_f64(),
            rep.simplify_seconds,
            rep.operators_seconds,
            rep.prolongation_nnz
        );
        let l0 = &h.levels()[0];
        let n = l0.n_vertices();
        let total: f64 = l0.mass.iter().sum();
        // Problem A: L + eps M (M normalized to total area 1), random rhs.
        let mut a_eps = l0.laplacian.clone();
        a_eps
            .add_to_diagonal(&l0.mass.iter().map(|m| eps * m / total).collect::<Vec<_>>())
            .expect("diag");
        let mut rng = SplitMix64::new(7);
        let b_eps: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
        // Problem B: the reference's example, L z = M 1, z = 0 at 0, 1, 2.
        let a_dir = l0.laplacian.clone();
        let mut b_dir = l0.mass.clone();
        let known = [0u32, 1, 2];
        for &k in &known {
            b_dir[k as usize] = 0.0;
        }
        let z0: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let problems: [(&str, &SparseMatrix, &[f64], &[u32]); 2] = [
            ("L+eps*M", &a_eps, &b_eps, &[]),
            ("dirichlet", &a_dir, &b_dir, &known),
        ];
        for (pname, a, b, kn) in problems {
            let vals = vec![0.0; kn.len()];
            let (xd, tdirect) = direct(a, b, kn, &vals);
            println!("  {pname}: direct (submatrix + factor + solve) {tdirect:.4} s");
            for sweeps in [1usize, 2, 3] {
                let opts = MultigridOptions {
                    pre_sweeps: sweeps,
                    post_sweeps: sweeps,
                    tolerance: tol,
                    max_cycles: 200,
                    ..MultigridOptions::default()
                };
                let t = Instant::now();
                let mg = h
                    .multigrid(a, kn, opts)
                    .unwrap_or_else(|e| die(format!("{name}: {e}")));
                let setup = t.elapsed().as_secs_f64();
                let x0 = if kn.is_empty() { None } else { Some(&z0[..]) };
                let (x, r) = mg.solve(b, &vals, x0).expect("solve");
                let abs0 = r.residuals[0];
                let abs_n = *r.residuals.last().unwrap();
                let abs_cycles = r.residuals.iter().position(|&v| v < 4.2e-11);
                println!(
                    "    sweeps {sweeps}+{sweeps} vcycle: setup {setup:.4} s, {} cycles, factor {:.3}, residual {abs0:.3e} -> {abs_n:.3e} (rel {:.2e}), solve {:.4} s, diff vs direct {:.2e}, converged {}, cycles to abs 4.2e-11: {abs_cycles:?}",
                    r.cycles,
                    r.convergence_factor,
                    r.relative_residual,
                    r.solve_seconds,
                    rel_diff(&x, &xd),
                    r.converged
                );
                let (x, r) = mg.pcg(b, &vals, x0, 200).expect("pcg");
                println!(
                    "    sweeps {sweeps}+{sweeps} pcg:    {} iterations, factor {:.3}, rel {:.2e}, solve {:.4} s, diff vs direct {:.2e}",
                    r.cycles,
                    r.convergence_factor,
                    r.relative_residual,
                    r.solve_seconds,
                    rel_diff(&x, &xd)
                );
            }
        }
    }
}
