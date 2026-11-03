//! The GPU back end (`NOTES-GPU.md`) against the CPU: products of `B̃`, `B̃⁰` (with the obstacle
//! diagonal) and `L̃^V_{2−s}` bitwise equal at determinism level 2 (`MatchCpu`), within round-off
//! and bitwise reproducible at level 1 (`RunToRun`); the pairwise-tree dot product; device
//! selection errors. Needs `--features cuda` and a CUDA device; without a device the tests print
//! a note and pass unless `DDG_REQUIRE_GPU` is set.

#![cfg(feature = "cuda")]

use std::sync::Arc;

use ddg_core::generators::{perturb, torus};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::mean_edge_length;
use ddg_repulsive_surfaces::gpu::{CudaOptions, Determinism, Device, GpuContext, pairwise_dot, pairwise_sum};
use ddg_repulsive_surfaces::metric::FaceDiagonal;
use ddg_repulsive_surfaces::metric::hierarchical::{ClusterTreeOptions, FaceMetricTree, VertexFractionalTree};
use ddg_repulsive_surfaces::{Exponents, Surface};

fn gpu(determinism: Determinism) -> Option<Arc<GpuContext>> {
    match GpuContext::new(&CudaOptions { ordinal: 0, determinism }) {
        Ok(g) => Some(g),
        Err(e) => {
            assert!(std::env::var_os("DDG_REQUIRE_GPU").is_none(), "no GPU: {e}");
            eprintln!("skipped (no usable GPU): {e}");
            None
        }
    }
}

fn mesh(faces: usize, seed: u64) -> (Surface, Vec<[f64; 3]>) {
    let nm = ((faces as f64 / 5.0).sqrt().round() as usize).max(3);
    let nmaj = ((faces as f64 / (2.0 * nm as f64)).round() as usize).max(3);
    let m = torus(nmaj, nm, 1.5, 0.6).unwrap();
    let s = Surface::new(m.positions.len(), m.faces).unwrap();
    let mut p = m.positions;
    let l = mean_edge_length(&s, &p);
    perturb(&mut p, [0.1 * l; 3], seed);
    (s, p)
}

fn vectors(n: usize, k: usize, seed: usize) -> Vec<Vec<f64>> {
    (0..k)
        .map(|q| (0..n).map(|i| ((i * (q + 3) + seed) as f64 * 0.37).sin() + 0.1 * (q as f64)).collect())
        .collect()
}

fn max_rel(a: &[f64], b: &[f64]) -> f64 {
    let scale = b.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
    a.iter().zip(b).fold(0.0f64, |m, (x, y)| m.max((x - y).abs())) / scale
}

struct Case {
    w_hi: f64,
    w_lo: f64,
    diag: Option<FaceDiagonal>,
}

fn cases(nf: usize) -> Vec<Case> {
    let d = |k: f64| (0..nf).map(|s| 0.5 + ((s as f64) * k).cos().abs()).collect::<Vec<f64>>();
    vec![
        Case { w_hi: 1.0, w_lo: 0.7, diag: None },
        Case { w_hi: 1.3, w_lo: 0.0, diag: None },
        Case { w_hi: 0.0, w_lo: 1.0, diag: None },
        Case { w_hi: 2.0, w_lo: 0.5, diag: Some(FaceDiagonal { hi: d(0.3), lo: d(0.7) }) },
        Case { w_hi: 1.0, w_lo: 0.0, diag: Some(FaceDiagonal { hi: Vec::new(), lo: d(1.1) }) },
    ]
}

/// Face metric products: level 2 bitwise equal to the CPU (several meshes, χ, leaf sizes, batch
/// sizes 1–6, with and without diagonals); level 1 within 1e-12 and bitwise reproducible.
#[test]
fn face_metric_products_match_cpu() {
    let Some(g2) = gpu(Determinism::MatchCpu) else { return };
    let Some(g1) = gpu(Determinism::RunToRun) else { return };
    let pool = make_pool(Some(3)).unwrap();
    let ex = Exponents::default();
    for (faces, seed, chi, leaf) in [(2000, 1, 0.5, 8), (5000, 2, 0.5, 8), (3000, 3, 0.25, 4), (1500, 4, 0.0, 8), (4000, 5, 1.0, 16)] {
        let (s, p) = mesh(faces, seed);
        let opts = ClusterTreeOptions { chi, leaf_size: leaf, leaf_pairs: leaf * leaf };
        let tree = FaceMetricTree::build(&pool, &s, &p, &ex, &opts).unwrap();
        let nd = 3 * s.n_vertices();
        for case in cases(s.n_faces()) {
            let d2 = tree.upload(&g2, case.diag.as_ref()).unwrap();
            let d1 = tree.upload(&g1, case.diag.as_ref()).unwrap();
            for nb in [1usize, 2, 3, 6] {
                let xs = vectors(nd, nb, faces);
                let xr: Vec<&[f64]> = xs.iter().map(Vec::as_slice).collect();
                let mut cpu = vec![vec![0.0; nd]; nb];
                {
                    let mut yr: Vec<&mut [f64]> = cpu.iter_mut().map(Vec::as_mut_slice).collect();
                    tree.apply_batch(&pool, case.w_hi, case.w_lo, case.diag.as_ref(), &xr, &mut yr).unwrap();
                }
                let run = |d: &ddg_repulsive_surfaces::gpu::GpuFaceMetric| {
                    let mut ys = vec![vec![0.0; nd]; nb];
                    let mut yr: Vec<&mut [f64]> = ys.iter_mut().map(Vec::as_mut_slice).collect();
                    d.apply_batch(case.w_hi, case.w_lo, &xr, &mut yr).unwrap();
                    ys
                };
                let y2 = run(&d2);
                for q in 0..nb {
                    let same = y2[q].iter().zip(&cpu[q]).all(|(a, b)| a.to_bits() == b.to_bits());
                    assert!(
                        same,
                        "level 2 not bitwise: faces {faces} chi {chi} leaf {leaf} w ({}, {}) nb {nb} q {q}: max rel {:e}",
                        case.w_hi,
                        case.w_lo,
                        max_rel(&y2[q], &cpu[q])
                    );
                }
                let y1 = run(&d1);
                let y1b = run(&d1);
                for q in 0..nb {
                    let e = max_rel(&y1[q], &cpu[q]);
                    assert!(e < 1e-12, "level 1 error {e:e} (faces {faces}, nb {nb})");
                    assert!(y1[q].iter().zip(&y1b[q]).all(|(a, b)| a.to_bits() == b.to_bits()), "level 1 not reproducible");
                }
            }
        }
    }
}

/// The vertex-pair fractional operator of the preconditioner: as above.
#[test]
fn vertex_fractional_products_match_cpu() {
    let Some(g2) = gpu(Determinism::MatchCpu) else { return };
    let Some(g1) = gpu(Determinism::RunToRun) else { return };
    let pool = make_pool(Some(2)).unwrap();
    let ex = Exponents::default();
    for (faces, seed, chi) in [(2000, 7, 0.5), (6000, 8, 0.5), (1000, 9, 0.0)] {
        let (s, p) = mesh(faces, seed);
        let opts = ClusterTreeOptions { chi, ..Default::default() };
        let t = VertexFractionalTree::for_surface(&pool, &s, &p, 2.0 - ex.s(), &opts).unwrap();
        let (d2, d1) = (t.upload(&g2).unwrap(), t.upload(&g1).unwrap());
        let nd = 3 * s.n_vertices();
        for nb in [1usize, 2, 5] {
            let xs = vectors(nd, nb, seed as usize);
            let xr: Vec<&[f64]> = xs.iter().map(Vec::as_slice).collect();
            let mut cpu = vec![vec![0.0; nd]; nb];
            {
                let mut yr: Vec<&mut [f64]> = cpu.iter_mut().map(Vec::as_mut_slice).collect();
                t.apply_batch(&pool, &xr, &mut yr).unwrap();
            }
            for (d, level) in [(&d2, 2), (&d1, 1)] {
                let mut ys = vec![vec![0.0; nd]; nb];
                {
                    let mut yr: Vec<&mut [f64]> = ys.iter_mut().map(Vec::as_mut_slice).collect();
                    d.apply_batch(&xr, &mut yr).unwrap();
                }
                for q in 0..nb {
                    if level == 2 {
                        assert!(
                            ys[q].iter().zip(&cpu[q]).all(|(a, b)| a.to_bits() == b.to_bits()),
                            "level 2 not bitwise (faces {faces}, nb {nb}): {:e}",
                            max_rel(&ys[q], &cpu[q])
                        );
                    } else {
                        assert!(max_rel(&ys[q], &cpu[q]) < 1e-12);
                    }
                }
            }
        }
    }
}

/// The pairwise tree: GPU = CPU bitwise for lengths around powers of two.
#[test]
fn pairwise_dot_matches_cpu() {
    let Some(g) = gpu(Determinism::RunToRun) else { return };
    for n in [0usize, 1, 2, 3, 255, 256, 257, 1000, 4097, 100_003] {
        let a: Vec<f64> = (0..n).map(|i| ((i as f64) * 0.71).sin() * 1e3_f64.powi((i % 7) as i32 - 3)).collect();
        let b: Vec<f64> = (0..n).map(|i| ((i as f64) * 1.3).cos()).collect();
        let c = pairwise_dot(&a, &b);
        let d = g.dot(&a, &b).unwrap();
        assert_eq!(c.to_bits(), d.to_bits(), "n = {n}: {c} vs {d}");
    }
    assert_eq!(pairwise_sum(&[1.0, 2.0, 3.0]), 6.0);
}

/// `Device` resolution: `Cpu` gives no context, `Cuda` with a missing ordinal is an error,
/// `Auto` with a missing ordinal falls back to the CPU.
#[test]
fn device_selection() {
    assert!(Device::Cpu.resolve().unwrap().is_none());
    let bad = CudaOptions { ordinal: 99, determinism: Determinism::RunToRun };
    assert!(matches!(Device::Cuda(bad).resolve(), Err(ddg_repulsive_surfaces::RsError::Gpu { .. })));
    assert!(Device::Auto(bad).resolve().unwrap().is_none());
}
