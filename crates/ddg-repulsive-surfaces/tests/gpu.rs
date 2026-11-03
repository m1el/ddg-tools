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

/// Level 1's device-built near fields (`build_on_gpu`): `H̃a` and the products within round-off
/// of the CPU-built trees (the device `pow` differs by ulps), bitwise reproducible; the
/// structure-only CPU trees refuse CPU products.
#[test]
fn gpu_built_near_fields_match_cpu() {
    use ddg_repulsive_surfaces::geometry::{FaceGeometry, vertex_areas};
    use ddg_repulsive_surfaces::metric::FaceOperators;
    use ddg_repulsive_surfaces::spatial::FaceTree;
    let Some(g1) = gpu(Determinism::RunToRun) else { return };
    let pool = make_pool(Some(3)).unwrap();
    let ex = Exponents::default();
    for (faces, seed, chi) in [(3000, 21, 0.5), (6000, 22, 0.5), (1200, 23, 0.0)] {
        let (s, p) = mesh(faces, seed);
        let opts = ClusterTreeOptions { chi, ..Default::default() };
        let cpu = FaceMetricTree::build(&pool, &s, &p, &ex, &opts).unwrap();
        let geom = FaceGeometry::new(&s, &p).unwrap();
        let ft = FaceTree::build(&s, &p, &geom, opts.leaf_size).unwrap();
        let ops = FaceOperators::from_geometry(&s, &p, geom);
        let diag = cases(s.n_faces()).pop().unwrap().diag;
        let (gt, dev) = FaceMetricTree::build_on_gpu(&pool, ft, ops, &ex, &opts, &g1, diag.as_ref()).unwrap();
        assert!(!gt.has_cpu_near_field());
        for (a, b) in gt.kernel_times_areas().iter().zip(cpu.kernel_times_areas()) {
            for k in 0..2 {
                assert!((a[k] - b[k]).abs() <= 1e-13 * b[k].abs(), "H~a {a:?} vs {b:?}");
            }
        }
        let nd = 3 * s.n_vertices();
        let xs = vectors(nd, 2, seed as usize);
        let xr: Vec<&[f64]> = xs.iter().map(Vec::as_slice).collect();
        let mut want = vec![vec![0.0; nd]; 2];
        {
            let mut yr: Vec<&mut [f64]> = want.iter_mut().map(Vec::as_mut_slice).collect();
            cpu.apply_batch(&pool, 1.0, 0.7, diag.as_ref(), &xr, &mut yr).unwrap();
        }
        let run = || {
            let mut ys = vec![vec![0.0; nd]; 2];
            let mut yr: Vec<&mut [f64]> = ys.iter_mut().map(Vec::as_mut_slice).collect();
            dev.apply_batch(1.0, 0.7, &xr, &mut yr).unwrap();
            ys
        };
        let (y, y2) = (run(), run());
        for q in 0..2 {
            let e = max_rel(&y[q], &want[q]);
            assert!(e < 1e-11, "faces {faces}: device-built product error {e:e}");
            assert!(y[q].iter().zip(&y2[q]).all(|(a, b)| a.to_bits() == b.to_bits()));
        }
        let mut out = vec![vec![0.0; nd]; 1];
        let mut yr: Vec<&mut [f64]> = out.iter_mut().map(Vec::as_mut_slice).collect();
        assert!(gt.apply_batch(&pool, 1.0, 0.7, None, &xr[..1], &mut yr).is_err());
        // The vertex-pair operator.
        let mass = vertex_areas(&s, &FaceGeometry::new(&s, &p).unwrap().areas);
        let tau = 2.0 - ex.s();
        let vc = VertexFractionalTree::build(&pool, &p, &mass, tau, &opts).unwrap();
        let (_, vd) = VertexFractionalTree::build_on_gpu(&pool, &p, &mass, tau, &opts, &g1).unwrap();
        let mut vw = vec![vec![0.0; nd]; 2];
        {
            let mut yr: Vec<&mut [f64]> = vw.iter_mut().map(Vec::as_mut_slice).collect();
            vc.apply_batch(&pool, &xr, &mut yr).unwrap();
        }
        let mut vy = vec![vec![0.0; nd]; 2];
        {
            let mut yr: Vec<&mut [f64]> = vy.iter_mut().map(Vec::as_mut_slice).collect();
            vd.apply_batch(&xr, &mut yr).unwrap();
        }
        for q in 0..2 {
            let e = max_rel(&vy[q], &vw[q]);
            assert!(e < 1e-11, "faces {faces}: device-built vertex product error {e:e}");
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

struct Trajectory {
    energies: Vec<u64>,
    iterations: Vec<usize>,
    faces: Vec<usize>,
    positions: Vec<[f64; 3]>,
}

fn flow(device: Device, steps: usize, faces: usize) -> Trajectory {
    use ddg_repulsive_surfaces::{Backend, Constraint, FlowOptions, HierarchicalOptions, Session, StopRule};
    let (s, p) = mesh(faces, 11);
    let o = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        backend: Backend::Hierarchical(HierarchicalOptions { device, ..HierarchicalOptions::default() }),
        threads: Some(3),
        stop: StopRule::Off,
        ..FlowOptions::default()
    };
    let mut session = Session::new(s, p, o).unwrap();
    let mut t = Trajectory { energies: Vec::new(), iterations: Vec::new(), faces: Vec::new(), positions: Vec::new() };
    for _ in 0..steps {
        let r = session.step().unwrap();
        t.energies.push(r.energy_after.to_bits());
        t.iterations.push(r.linear.total_iterations());
        t.faces.push(session.surface().n_faces());
    }
    t.positions = session.positions().to_vec();
    t
}

/// Whole flows (remeshing on): level 2 gives the CPU trajectory bit for bit (energies, GMRES
/// iterations, face counts, positions); level 1 is bitwise reproducible and stays within
/// round-off effects of the CPU (iterations within ±2 per step, energies within 1e-6 relative
/// over a few steps).
#[test]
fn flow_trajectories_match_cpu() {
    if gpu(Determinism::MatchCpu).is_none() {
        return;
    }
    let (steps, faces) = (6, 4000);
    let cpu = flow(Device::Cpu, steps, faces);
    let l2 = flow(Device::Cuda(CudaOptions { ordinal: 0, determinism: Determinism::MatchCpu }), steps, faces);
    assert_eq!(cpu.energies, l2.energies, "level 2 energies");
    assert_eq!(cpu.iterations, l2.iterations, "level 2 GMRES iterations");
    assert_eq!(cpu.faces, l2.faces);
    assert!(cpu.positions.iter().flatten().zip(l2.positions.iter().flatten()).all(|(a, b)| a.to_bits() == b.to_bits()));
    let l1 = flow(Device::Cuda(CudaOptions { ordinal: 0, determinism: Determinism::RunToRun }), steps, faces);
    let l1b = flow(Device::Cuda(CudaOptions { ordinal: 0, determinism: Determinism::RunToRun }), steps, faces);
    assert_eq!(l1.energies, l1b.energies, "level 1 reproducible");
    assert_eq!(l1.iterations, l1b.iterations);
    assert!(l1.positions.iter().flatten().zip(l1b.positions.iter().flatten()).all(|(a, b)| a.to_bits() == b.to_bits()));
    for k in 0..steps {
        let (a, b) = (f64::from_bits(cpu.energies[k]), f64::from_bits(l1.energies[k]));
        assert!((a - b).abs() <= 1e-6 * a.abs(), "step {k}: energy {a} vs {b}");
        assert!(cpu.iterations[k].abs_diff(l1.iterations[k]) <= 2, "step {k}: iterations {:?} vs {:?}", cpu.iterations, l1.iterations);
    }
    eprintln!(
        "cpu iterations {:?}, level 1 {:?}; final energy cpu {} level 1 {}",
        cpu.iterations,
        l1.iterations,
        f64::from_bits(cpu.energies[steps - 1]),
        f64::from_bits(l1.energies[steps - 1])
    );
}
