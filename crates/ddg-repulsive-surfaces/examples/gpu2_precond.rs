//! GPU-2 measurement of one preconditioner application (`NOTES-GPU2.md`): the fine `P⁻¹` on the
//! CPU (chunked and faer `L′` solves, the given threads) and on the GPU (GPU-1's split: fractional
//! product only; GPU-2: everything), for batches of 1, 2, 4 vectors; median of `--reps`.
//!
//! ```text
//! gpu2_precond <mesh> [--threads N] [--reps N] [--pins K]     <mesh>: tunnel8:FACES | fig12:I
//! ```

use std::sync::Arc;
use std::time::Instant;

use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::generate::shapes::{Knot, fig12, knotted_tunnel};
use ddg_repulsive_surfaces::gpu::{CudaOptions, Determinism, GpuContext};
use ddg_repulsive_surfaces::solve::precond::{PreconditionerOptions, SparsePreconditioner, TriangularSolve};
use ddg_repulsive_surfaces::{Exponents, Surface};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str, d: usize| -> usize {
        args.iter().position(|a| a == k).map(|i| args[i + 1].parse().unwrap()).unwrap_or(d)
    };
    let spec = args.first().cloned().unwrap_or_else(|| "fig12:3".into());
    let (threads, reps, pins) = (get("--threads", 4), get("--reps", 21), get("--pins", 0));
    let parts: Vec<&str> = spec.split(':').collect();
    let m = match parts[0] {
        "tunnel8" => knotted_tunnel(Knot::FigureEight, parts[1].parse().unwrap()).unwrap().mesh,
        _ => fig12(parts[1].parse().unwrap()).unwrap().mesh,
    };
    let (s, p) = (Surface::new(m.positions.len(), m.faces).unwrap(), m.positions);
    let pool = Arc::new(make_pool(Some(threads)).unwrap());
    let mut cons = vec![Constraint::Barycenter { component: None, target: None }];
    if let Some(step) = s.n_vertices().checked_div(pins) {
        cons.push(Constraint::Pins { vertices: (0..pins).map(|k| (k * step) as u32).collect(), targets: None });
    }
    let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
    let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
    let e = Exponents::default();
    let build = |o: &PreconditionerOptions| {
        SparsePreconditioner::new(pool.clone(), &s, &p, &e, &cq, o, None, None).unwrap()
    };
    let chunked = PreconditionerOptions::default();
    let faer = PreconditionerOptions { triangular: TriangularSolve::Faer, ..Default::default() };
    let split = PreconditionerOptions { gpu_solves: false, ..Default::default() };
    let g1 = GpuContext::new(&CudaOptions { ordinal: 0, determinism: Determinism::RunToRun }).unwrap();
    let g2 = GpuContext::new(&CudaOptions { ordinal: 0, determinism: Determinism::MatchCpu }).unwrap();
    let variants: Vec<(&str, SparsePreconditioner)> = vec![
        ("cpu chunked", build(&chunked)),
        ("cpu faer", build(&faer)),
        ("gpu-1 split L1", build(&split).with_gpu(&g1).unwrap()),
        ("gpu-2 L1", build(&chunked).with_gpu(&g1).unwrap()),
        ("gpu-2 L2", build(&chunked).with_gpu(&g2).unwrap()),
    ];
    let n = cq.ncols() + cq.nrows();
    println!("mesh {spec}: {} faces, dim {n}, threads {threads}", s.n_faces());
    if let ddg_repulsive_surfaces::solve::precond::FractionalFactor::Vertices(t) = variants[0].1.fractional() {
        println!("vertex-pair tree: {:?}", t.stats());
        let gv = t.upload(&g1).unwrap();
        let nv = s.n_vertices();
        for nb in [1usize, 2, 4, 8] {
            let xs: Vec<Vec<f64>> = (0..nb).map(|q| (0..3 * nv).map(|i| ((i + q) as f64).sin()).collect()).collect();
            let xr: Vec<&[f64]> = xs.iter().map(Vec::as_slice).collect();
            let mut ys = vec![vec![0.0; 3 * nv]; nb];
            let mut tt: Vec<f64> = (0..reps)
                .map(|_| {
                    let mut yr: Vec<&mut [f64]> = ys.iter_mut().map(Vec::as_mut_slice).collect();
                    let t0 = Instant::now();
                    gv.apply_batch(&xr, &mut yr).unwrap();
                    t0.elapsed().as_secs_f64()
                })
                .collect();
            tt.sort_by(f64::total_cmp);
            println!("  vertex-pair product on the GPU (L1, with PCIe), {nb} vectors: {:.3} ms", 1e3 * tt[reps / 2]);
        }
    }
    for nb in [1usize, 2, 4] {
        let rs: Vec<Vec<f64>> = (0..nb)
            .map(|q| (0..n).map(|i| ((i * (q + 3)) as f64 * 0.37).sin()).collect())
            .collect();
        let rr: Vec<&[f64]> = rs.iter().map(Vec::as_slice).collect();
        for (name, pc) in &variants {
            pc.apply_many(&rr).unwrap();
            let mut t: Vec<f64> = (0..reps)
                .map(|_| {
                    let t0 = Instant::now();
                    pc.apply_many(&rr).unwrap();
                    t0.elapsed().as_secs_f64()
                })
                .collect();
            t.sort_by(f64::total_cmp);
            println!("  {nb} vectors, {name:>16}: {:8.3} ms", 1e3 * t[reps / 2]);
            if let Some(g) = pc.gpu_preconditioner() {
                let (parts, whole) = g.profile_solve(&rr).unwrap();
                println!(
                    "      one L' solve: {:.3} ms; synchronized per launch: gather/scatter {:.3}, partials {:.3}, fwd chains {:.3}, bwd chains {:.3} ms",
                    1e3 * whole,
                    1e3 * parts[0],
                    1e3 * parts[1],
                    1e3 * parts[2],
                    1e3 * parts[3]
                );
            }
        }
    }
}
