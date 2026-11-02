//! Hierarchical products, the sparse preconditioner and the iterative back end (`SPEC.md` §8.3,
//! §6.3, §13.3): `χ = 0` exactness against WP2's dense operators, the compressed null space,
//! symmetry, the error in `χ`, GMRES + Schur against the dense direct solve, determinism across
//! thread counts, `Backend::Auto` switching.

#![allow(clippy::needless_range_loop, clippy::type_complexity)]

use std::sync::Arc;

use ddg_core::SparseMatrix;
use ddg_core::generators::{GeneratedMesh, SplitMix64, bumpy_grid_disk, icosphere, perturb, torus};
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_surfaces::geometry::{FaceGeometry, mean_edge_length, vertex_areas};
use ddg_repulsive_surfaces::metric::dense::{
    fractional_kernel, high_order_kernel, kernel_matrix_apply, low_order_kernel,
    vertex_fractional_operator,
};
use ddg_repulsive_surfaces::metric::hierarchical::{
    ClusterTreeOptions, FaceFractionalTree, FaceMetricTree, VertexFractionalTree,
};
use ddg_repulsive_surfaces::metric::{
    Channel, DenseMatrix, FaceDiagonal, FaceOperators, FractionalForm, SobolevWeights,
};
use ddg_repulsive_surfaces::solve::iterative::{
    Backend, HierarchicalOptions, MetricExtras, dense_metric_with_extras, hierarchical_saddle,
    saddle_solver,
};
use ddg_repulsive_surfaces::solve::precond::{
    LaplacianSaddle, PreconditionerForm, PreconditionerOptions, SparsePreconditioner,
};
use ddg_repulsive_surfaces::solve::{
    DensePreconditioner, DenseSaddle, DenseSaddleOperator, GmresOptions, LinearOperator,
    Preconditioner, SaddleSolver, SchurComplement,
};
use ddg_repulsive_surfaces::{Exponents, Surface};

fn perturbed(m: GeneratedMesh, amp: f64, seed: u64) -> (Surface, Vec<[f64; 3]>) {
    let mut p = m.positions;
    let s = Surface::new(p.len(), m.faces).unwrap();
    let l = mean_edge_length(&s, &p);
    perturb(&mut p, [amp * l; 3], seed);
    (s, p)
}

fn two_spheres() -> (Surface, Vec<[f64; 3]>) {
    let m = icosphere(1).unwrap();
    let mut p = m.positions.clone();
    let mut f = m.faces.clone();
    let off = p.len() as u32;
    p.extend(m.positions.iter().map(|q| [0.7 * q[0] + 2.2, 0.7 * q[1] + 0.3, 0.7 * q[2]]));
    f.extend(m.faces.iter().map(|t| [t[0] + off, t[1] + off, t[2] + off]));
    perturbed(GeneratedMesh { positions: p, faces: f }, 0.1, 7)
}

/// Small test meshes: (name, surface, positions).
fn meshes() -> Vec<(&'static str, Surface, Vec<[f64; 3]>)> {
    let (s1, p1) = perturbed(icosphere(2).unwrap(), 0.2, 1);
    let (s2, p2) = perturbed(torus(14, 8, 1.5, 0.6).unwrap(), 0.1, 2);
    let (s3, p3) = perturbed(bumpy_grid_disk(8, 8, 0.3).unwrap(), 0.1, 3);
    let (s4, p4) = two_spheres();
    vec![
        ("icosphere 320", s1, p1),
        ("torus 224", s2, p2),
        ("bumpy disk", s3, p3),
        ("two spheres", s4, p4),
    ]
}

fn random_vec(n: usize, seed: u64) -> Vec<f64> {
    let mut r = SplitMix64::new(seed);
    (0..n).map(|_| r.uniform(-1.0, 1.0)).collect()
}

fn to3(x: &[f64]) -> Vec<[f64; 3]> {
    x.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()
}

fn flat(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn rel(a: &[f64], b: &[f64]) -> f64 {
    let d: Vec<f64> = a.iter().zip(b).map(|(x, y)| x - y).collect();
    norm(&d) / norm(b).max(1e-300)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn dense_apply3(m: &DenseMatrix, x: &[f64]) -> Vec<f64> {
    let n = m.n();
    let mut y = vec![0.0; 3 * n];
    for i in 0..n {
        let row = m.row(i);
        for j in 0..n {
            for k in 0..3 {
                y[3 * i + k] += row[j] * x[3 * j + k];
            }
        }
    }
    y
}

fn opts(chi: f64, leaf: usize) -> ClusterTreeOptions {
    ClusterTreeOptions {
        chi,
        leaf_size: leaf,
        ..Default::default()
    }
}

/// The three hierarchical products and their dense counterparts (WP2) at `x`.
struct Products {
    b: (Vec<f64>, Vec<f64>),
    b0: (Vec<f64>, Vec<f64>),
    lf: (Vec<f64>, Vec<f64>),
    lv: (Vec<f64>, Vec<f64>),
}

fn products(
    pool: &rayon::ThreadPool,
    s: &Surface,
    p: &[[f64; 3]],
    e: &Exponents,
    o: &ClusterTreeOptions,
    x: &[f64],
) -> Products {
    let n = x.len();
    let ops = FaceOperators::new(s, p).unwrap();
    let tree = FaceMetricTree::build(pool, s, p, e, o).unwrap();
    let mut b = vec![0.0; n];
    tree.apply_high(pool, x, &mut b).unwrap();
    let mut b0 = vec![0.0; n];
    tree.apply_low(pool, x, &mut b0).unwrap();
    let tau = 2.0 - e.s();
    let ft = FaceFractionalTree::build(pool, tree.face_tree(), &ops, tau, o).unwrap();
    let mut lf = vec![0.0; n];
    ft.apply(pool, x, &mut lf).unwrap();
    let vt = VertexFractionalTree::for_surface(pool, s, p, tau, o).unwrap();
    let mut lv = vec![0.0; n];
    vt.apply(pool, x, &mut lv).unwrap();

    let x3 = to3(x);
    let hb = high_order_kernel(pool, &ops, e).unwrap();
    let db = flat(&kernel_matrix_apply(pool, &ops, &hb, Channel::Gradient, &x3).unwrap());
    let hl = low_order_kernel(pool, &ops, e).unwrap();
    let dl = flat(&kernel_matrix_apply(pool, &ops, &hl, Channel::Average, &x3).unwrap());
    let hf = fractional_kernel(pool, &ops, tau).unwrap();
    let df = flat(&kernel_matrix_apply(pool, &ops, &hf, Channel::Average, &x3).unwrap());
    let geom = FaceGeometry::new(s, p).unwrap();
    let mass = vertex_areas(s, &geom.areas);
    let lvd = vertex_fractional_operator(pool, p, &mass, tau).unwrap();
    let dv = dense_apply3(&lvd, x);
    Products {
        b: (b, db),
        b0: (b0, dl),
        lf: (lf, df),
        lv: (lv, dv),
    }
}

#[test]
fn chi_zero_is_exact() {
    let pool = make_pool(Some(4)).unwrap();
    for e in [Exponents::default(), Exponents::new(4.5, 9.5).unwrap()] {
        for (name, s, p) in meshes() {
            for leaf in [1, 3, 8] {
                let x = random_vec(3 * s.n_vertices(), 11);
                let pr = products(&pool, &s, &p, &e, &opts(0.0, leaf), &x);
                let errs = [
                    rel(&pr.b.0, &pr.b.1),
                    rel(&pr.b0.0, &pr.b0.1),
                    rel(&pr.lf.0, &pr.lf.1),
                    rel(&pr.lv.0, &pr.lv.1),
                ];
                eprintln!("{name} leaf {leaf} s {:.3}: B {:.1e} B0 {:.1e} Lf {:.1e} Lv {:.1e}", e.s(), errs[0], errs[1], errs[2], errs[3]);
                for err in errs {
                    assert!(err < 1e-12, "{name} leaf {leaf}: {errs:?}");
                }
            }
        }
    }
}

#[test]
fn error_decays_with_chi_and_is_small_at_default() {
    let pool = make_pool(Some(4)).unwrap();
    let e = Exponents::default();
    let (s, p) = perturbed(torus(if cfg!(debug_assertions) { 30 } else { 40 }, if cfg!(debug_assertions) { 16 } else { 20 }, 1.5, 0.6).unwrap(), 0.1, 5);
    let x = random_vec(3 * s.n_vertices(), 3);
    // A smooth input as well (the far field is accurate on smooth functions).
    let smooth: Vec<f64> = p
        .iter()
        .flat_map(|q| [q[0].sin() + q[2], q[1] * q[0], (2.0 * q[2]).cos()])
        .collect();
    let mut prev = [f64::INFINITY; 4];
    for chi in [1.0, 0.5, 0.25] {
        for (kind, input) in [("random", &x), ("smooth", &smooth)] {
            let pr = products(&pool, &s, &p, &e, &opts(chi, 8), input);
            let errs = [
                rel(&pr.b.0, &pr.b.1),
                rel(&pr.b0.0, &pr.b0.1),
                rel(&pr.lf.0, &pr.lf.1),
                rel(&pr.lv.0, &pr.lv.1),
            ];
            eprintln!("torus, chi {chi}, {kind}: B {:.2e} B0 {:.2e} Lf {:.2e} Lv {:.2e}", errs[0], errs[1], errs[2], errs[3]);
            // SPEC §13.3: the metric (B, B⁰) within 5 % at χ = 0.5. The preconditioner's
            // fractional factors are only reported (on smooth inputs `L u` is a small difference).
            if chi == 0.5 {
                for err in &errs[..2] {
                    assert!(*err < 5e-2, "chi 0.5 {kind}: {errs:?}");
                }
            }
            if kind == "random" {
                for k in 0..4 {
                    assert!(errs[k] <= prev[k] * 1.01 + 1e-15, "not decreasing: {errs:?} after {prev:?}");
                }
                prev = errs;
            }
        }
    }
}

#[test]
fn compressed_null_space_and_symmetry() {
    let pool = make_pool(Some(4)).unwrap();
    let e = Exponents::default();
    for (name, s, p) in meshes() {
        let o = opts(0.5, 4);
        let n = 3 * s.n_vertices();
        let ones = vec![1.0; n];
        let ops = FaceOperators::new(&s, &p).unwrap();
        let tree = FaceMetricTree::build(&pool, &s, &p, &e, &o).unwrap();
        assert!(tree.stats().admissible_blocks > 0, "{name}: no far field");
        let tau = 2.0 - e.s();
        let ft = FaceFractionalTree::build(&pool, tree.face_tree(), &ops, tau, &o).unwrap();
        let vt = VertexFractionalTree::for_surface(&pool, &s, &p, tau, &o).unwrap();
        let x = random_vec(n, 5);
        let mut y = vec![0.0; n];
        // Low-order and fractional operators: exactly zero on constants.
        tree.apply_low(&pool, &ones, &mut y).unwrap();
        assert!(y.iter().all(|v| *v == 0.0), "{name}: B0 1 != 0");
        ft.apply(&pool, &ones, &mut y).unwrap();
        assert!(y.iter().all(|v| *v == 0.0), "{name}: Lf 1 != 0");
        vt.apply(&pool, &ones, &mut y).unwrap();
        assert!(y.iter().all(|v| *v == 0.0), "{name}: Lv 1 != 0");
        // B 1 = 0 up to D_f 1 = 0 (round-off).
        let mut bx = vec![0.0; n];
        tree.apply_high(&pool, &x, &mut bx).unwrap();
        tree.apply_high(&pool, &ones, &mut y).unwrap();
        assert!(norm(&y) <= 1e-12 * norm(&bx) * norm(&ones) / norm(&x), "{name}: B 1 = {:e}", norm(&y));
        // Symmetry of the compressed operators.
        let u = random_vec(n, 6);
        let v = random_vec(n, 7);
        type Apply<'a> = Box<dyn Fn(&[f64], &mut [f64]) + 'a>;
        let cases: Vec<(&str, Apply)> = vec![
            ("B + B0", Box::new(|a: &[f64], b: &mut [f64]| tree.apply(&pool, 1.0, 1.0, None, a, b).unwrap())),
            ("Lf", Box::new(|a: &[f64], b: &mut [f64]| ft.apply(&pool, a, b).unwrap())),
            ("Lv", Box::new(|a: &[f64], b: &mut [f64]| vt.apply(&pool, a, b).unwrap())),
        ];
        for (op, f) in cases {
            let mut au = vec![0.0; n];
            let mut av = vec![0.0; n];
            f(&u, &mut au);
            f(&v, &mut av);
            let (uav, vau) = (dot(&u, &av), dot(&v, &au));
            let asym = (uav - vau).abs() / (norm(&u) * norm(&av)).max(norm(&v) * norm(&au));
            eprintln!("{name} {op}: asymmetry {asym:.1e}, uᵀAu {:.3e}", dot(&u, &au));
            assert!(asym < 1e-13, "{name} {op}: {asym:e}");
            assert!(dot(&u, &au) > 0.0, "{name} {op}: not positive");
        }
    }
}

#[test]
fn extras_and_weights_match_the_dense_metric() {
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let e = Exponents::default();
    for (name, s, p) in meshes() {
        let nf = s.n_faces();
        let mut r = SplitMix64::new(9);
        let obstacle = FaceDiagonal {
            hi: (0..nf).map(|_| r.uniform(0.0, 2.0)).collect(),
            lo: (0..nf).map(|_| r.uniform(0.0, 2.0)).collect(),
        };
        let extras = MetricExtras { willmore: None, obstacle: Some(obstacle) };
        for weights in [
            SobolevWeights { tp_weight: 2.5, low_order: true, willmore_weight: 0.3 },
            SobolevWeights { tp_weight: 1.0, low_order: false, willmore_weight: 0.0 },
        ] {
            let set = ConstraintSet::new(&s, &p, &[Constraint::Barycenter { component: None, target: None }], BarycenterMode::PerComponent).unwrap();
            let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
            let mut ho = HierarchicalOptions::default();
            ho.cluster.chi = 0.0;
            let solver = hierarchical_saddle(pool.clone(), &s, &p, &e, &weights, &extras, &cq, &ho, None).unwrap();
            let a = dense_metric_with_extras(&pool, &s, &p, &e, &weights, &extras).unwrap();
            let dop = DenseSaddleOperator::new(&a, &cq).unwrap();
            let n = dop.dim();
            let x = random_vec(n, 4);
            let mut y1 = vec![0.0; n];
            let mut y2 = vec![0.0; n];
            solver.operator().apply(&x, &mut y1).unwrap();
            dop.apply(&x, &mut y2).unwrap();
            let err = rel(&y1, &y2);
            eprintln!("{name} {weights:?}: |A'x - dense| rel {err:.1e}");
            assert!(err < 1e-12, "{name}: {err:e}");
        }
    }
}

#[test]
fn sparse_preconditioner_matches_the_dense_one() {
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let e = Exponents::default();
    for (name, s, p) in meshes() {
        let cons = vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Pins { vertices: vec![1], targets: None },
        ];
        let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
        let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
        for pairs in [FractionalForm::VertexPairs, FractionalForm::FacePairs] {
            let mut po = PreconditionerOptions {
                fractional: pairs,
                cluster: opts(0.0, 8),
                form: PreconditionerForm::Congruence,
                // The exact preconditioner (refined L' solves), compared to round-off.
                refine: true,
                ..Default::default()
            };
            let sp = SparsePreconditioner::new(pool.clone(), &s, &p, &e, &cq, &po, None, None).unwrap();
            assert!(sp.laplacian_saddle().is_per_axis());
            assert!(!sp.epsilon().1);
            let dp = DensePreconditioner::new(&pool, &s, &p, &e, &cq, po.eps_rel, pairs).unwrap();
            let n = sp.dim();
            let nd = 3 * s.n_vertices();
            let r = random_vec(n, 8);
            let z1 = sp.apply_full(&r).unwrap();
            let z2 = dp.apply_full(&r).unwrap();
            let err = rel(&z1, &z2);
            // C_Q z = 0.
            let mut cz = vec![0.0; cq.nrows()];
            cq.mul_vec(&z1[..nd], &mut cz).unwrap();
            eprintln!("{name} {pairs:?}: |P r - dense| rel {err:.1e}, |C z| {:.1e}", norm(&cz) / norm(&z1));
            assert!(err < 1e-9, "{name} {pairs:?}: {err:e}");
            assert!(norm(&cz) < 1e-10 * norm(&z1), "{name}");
            // Symmetric: uᵀ P v = vᵀ P u.
            let u = random_vec(n, 21);
            let pu = sp.apply_full(&u).unwrap();
            let asym = (dot(&u, &z1) - dot(&r, &pu)).abs() / (norm(&u) * norm(&z1));
            assert!(asym < 1e-10, "{name}: preconditioner asymmetry {asym:e}");
            // The augmented form: symmetric, C_Q P⁻¹ [r; b] = b, and its position block on
            // force-free right sides with b = 0 is the congruence form's plus L′⁻¹'s.
            po.form = PreconditionerForm::Augmented;
            let sa = SparsePreconditioner::new(pool.clone(), &s, &p, &e, &cq, &po, None, None).unwrap();
            let za = sa.apply_full(&r).unwrap();
            let pa = sa.apply_full(&u).unwrap();
            let asym = (dot(&u, &za) - dot(&r, &pa)).abs() / (norm(&u) * norm(&za));
            cq.mul_vec(&za[..nd], &mut cz).unwrap();
            let cerr = rel(&cz, &r[nd..]);
            eprintln!("{name} {pairs:?} augmented: asymmetry {asym:.1e}, |C z - b|/|b| {cerr:.1e}");
            assert!(asym < 1e-10 && cerr < 1e-10, "{name}: {asym:e} {cerr:e}");
        }
    }
}

#[test]
fn laplacian_saddle_solves_accurately() {
    // The unpivoted LDLᵀ of [[W + εM, Cᵀ], [C, 0]] with ε = 1e-10 relative, with refinement.
    for (name, s, p) in meshes() {
        let set = ConstraintSet::new(&s, &p, &[Constraint::Barycenter { component: None, target: None }], BarycenterMode::PerComponent).unwrap();
        let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
        let w = ddg_repulsive_surfaces::geometry::cotan_stiffness(&s, &p).unwrap();
        let geom = FaceGeometry::new(&s, &p).unwrap();
        let mass = vertex_areas(&s, &geom.areas);
        let eps = 1e-10 * w.diagonal().iter().sum::<f64>() / mass.iter().sum::<f64>();
        let mut trip: Vec<(u32, u32, f64)> = w.triplets().collect();
        for (v, m) in mass.iter().enumerate() {
            trip.push((v as u32, v as u32, eps * m));
        }
        let k = SparseMatrix::from_triplets(s.n_vertices(), s.n_vertices(), &trip).unwrap();
        let ls = LaplacianSaddle::new(&k, &cq).unwrap();
        let mut kd = DenseMatrix::zeros(s.n_vertices());
        for (i, j, x) in k.triplets() {
            kd.as_mut_slice()[i as usize * s.n_vertices() + j as usize] += x;
        }
        let op = DenseSaddleOperator::new(&kd, &cq).unwrap();
        let r = random_vec(ls.dim(), 2);
        let y = ls.solve(&r).unwrap();
        let mut ay = vec![0.0; r.len()];
        op.apply(&y, &mut ay).unwrap();
        let res = rel(&ay, &r);
        eprintln!("{name}: L' solve residual {res:.1e}");
        assert!(res < 1e-10, "{name}: {res:e}");
    }
}

#[test]
fn global_barycenter_on_several_components_uses_eps_multi() {
    let pool = Arc::new(make_pool(Some(2)).unwrap());
    let (s, p) = two_spheres();
    let e = Exponents::default();
    for (mode, multi) in [(BarycenterMode::PerComponent, false), (BarycenterMode::Global, true)] {
        let set = ConstraintSet::new(&s, &p, &[Constraint::Barycenter { component: None, target: None }], mode).unwrap();
        let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
        let sp = SparsePreconditioner::new(pool.clone(), &s, &p, &e, &cq, &PreconditionerOptions::default(), None, None).unwrap();
        assert_eq!(sp.epsilon().1, multi, "{mode:?}");
    }
}

struct Case {
    name: &'static str,
    s: Surface,
    p: Vec<[f64; 3]>,
    cons: Vec<Constraint>,
}

fn solve_cases() -> Vec<Case> {
    let (s1, p1) = perturbed(icosphere(2).unwrap(), 0.2, 1);
    let (s2, p2) = perturbed(torus(14, 8, 1.5, 0.6).unwrap(), 0.1, 2);
    let (s3, p3) = perturbed(bumpy_grid_disk(8, 8, 0.3).unwrap(), 0.1, 3);
    let (s4, p4) = two_spheres();
    vec![
        Case {
            name: "icosphere: barycenter, area, volume",
            s: s1,
            p: p1,
            cons: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: Some(13.0) },
                Constraint::Volume { component: None, target: Some(4.0) },
            ],
        },
        Case {
            name: "torus: pins, volume",
            s: s2,
            p: p2,
            cons: vec![
                Constraint::Pins { vertices: vec![0, 17, 33], targets: None },
                Constraint::Volume { component: None, target: None },
            ],
        },
        Case {
            name: "disk: barycenter, area",
            s: s3,
            p: p3,
            cons: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: Some(1.1) },
            ],
        },
        Case {
            name: "two spheres: barycenters, area of 1, volume of 0",
            s: s4,
            p: p4,
            cons: vec![
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: Some(1), target: None },
                Constraint::Volume { component: Some(0), target: Some(3.0) },
            ],
        },
    ]
}

#[test]
fn gmres_schur_equals_the_dense_direct_solve() {
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let e = Exponents::default();
    let weights = SobolevWeights::default();
    for case in solve_cases() {
        let (s, p) = (&case.s, &case.p);
        let set = ConstraintSet::new(s, p, &case.cons, BarycenterMode::PerComponent).unwrap();
        let q = set.positional_rows();
        let nl = set.nonlinear_rows();
        let cq = set.jacobian(p, &q).unwrap();
        let cn = set.jacobian(p, &nl).unwrap();
        let geom = FaceGeometry::new(s, p).unwrap();
        let (_, g) = energy_and_differential(&pool, s, p, &geom, &e).unwrap();
        let g = flat(&g);
        let phi = set.values(p).unwrap();
        let phin: Vec<f64> = nl.iter().map(|&r| phi[r]).collect();

        let a = dense_metric_with_extras(&pool, s, p, &e, &weights, &MetricExtras::default()).unwrap();
        let mut dense = DenseSaddle::new(&pool, &a, &cq).unwrap();
        let dschur = SchurComplement::new(&mut dense, nl.clone(), cn.clone()).unwrap();
        let ddir = dschur.direction(&mut dense, &g).unwrap();
        let dh = dschur.correction(&phin);

        for (chi, tol, bound) in [(0.0, 1e-12, 1e-8), (0.5, 1e-4, 5e-2)] {
            for pairs in [FractionalForm::VertexPairs, FractionalForm::FacePairs] {
                if pairs == FractionalForm::FacePairs && case.name.starts_with("disk") {
                    // The face form stalls on the bumpy disk (D32; WP2 measured a stall at 4e-2).
                    continue;
                }
                let mut ho = HierarchicalOptions::default();
                ho.cluster.chi = chi;
                ho.preconditioner.fractional = pairs;
                ho.gmres = GmresOptions { tolerance: tol, max_iterations: 500, ..Default::default() };
                let mut it = hierarchical_saddle(pool.clone(), s, p, &e, &weights, &MetricExtras::default(), &cq, &ho, None).unwrap();
                let schur = SchurComplement::new(&mut it, nl.clone(), cn.clone()).unwrap();
                let dir = schur.direction(&mut it, &g).unwrap();
                let h = schur.correction(&phin);
                let ex = rel(&dir.x, &ddir.x);
                let eh = rel(&h, &dh);
                let iters: Vec<usize> = schur.records().iter().map(|r| r.iterations).chain([dir.record.iterations]).collect();
                eprintln!("{} chi {chi} {pairs:?}: |x - x_dense| {ex:.1e}, |h - h_dense| {eh:.1e}, iterations {iters:?}", case.name);
                assert!(ex < bound && eh < bound, "{}: {ex:e} {eh:e}", case.name);
                // C x = 0 on all rows (up to the GMRES tolerance).
                let all: Vec<usize> = q.iter().chain(&nl).copied().collect();
                let c_all = set.jacobian(p, &all).unwrap();
                let mut cx = vec![0.0; all.len()];
                c_all.mul_vec(&dir.x, &mut cx).unwrap();
                let scale = c_all.values().iter().fold(0.0f64, |m, v| m.max(v.abs())) * norm(&dir.x);
                assert!(norm(&cx) < 10.0 * tol.max(1e-10) * scale, "{}: |Cx| {:e}", case.name, norm(&cx));
            }
        }
    }
}

fn bits(x: &[f64]) -> Vec<u64> {
    x.iter().map(|v| v.to_bits()).collect()
}

#[test]
fn bitwise_independent_of_threads() {
    let e = Exponents::default();
    let (s, p) = perturbed(torus(30, 16, 1.5, 0.6).unwrap(), 0.1, 4);
    let set = ConstraintSet::new(&s, &p, &[Constraint::Barycenter { component: None, target: None }, Constraint::Area { component: None, target: None }], BarycenterMode::PerComponent).unwrap();
    let q = set.positional_rows();
    let nl = set.nonlinear_rows();
    let cq = set.jacobian(&p, &q).unwrap();
    let cn = set.jacobian(&p, &nl).unwrap();
    let x = random_vec(3 * s.n_vertices(), 1);
    let mut reference: Option<Vec<Vec<u64>>> = None;
    for threads in [1, 3, 8] {
        let pool = Arc::new(make_pool(Some(threads)).unwrap());
        let geom = FaceGeometry::new(&s, &p).unwrap();
        let (_, g) = energy_and_differential(&pool, &s, &p, &geom, &e).unwrap();
        let g = flat(&g);
        let mut out = Vec::new();
        for pairs in [FractionalForm::VertexPairs, FractionalForm::FacePairs] {
            let mut ho = HierarchicalOptions::default();
            ho.preconditioner.fractional = pairs;
            let mut it = hierarchical_saddle(pool.clone(), &s, &p, &e, &SobolevWeights::default(), &MetricExtras::default(), &cq, &ho, None).unwrap();
            let mut y = vec![0.0; x.len()];
            it.operator().apply_metric(&x, &mut y).unwrap();
            out.push(bits(&y));
            let schur = SchurComplement::new(&mut it, nl.clone(), cn.clone()).unwrap();
            let dir = schur.direction(&mut it, &g).unwrap();
            out.push(bits(&dir.x));
        }
        match &reference {
            None => reference = Some(out),
            Some(r) => assert!(r == &out, "results differ at {threads} threads"),
        }
    }
}

#[test]
fn auto_backend_switches_on_the_face_count() {
    let pool = Arc::new(make_pool(Some(2)).unwrap());
    let e = Exponents::default();
    let (s, p) = perturbed(icosphere(2).unwrap(), 0.2, 1);
    let set = ConstraintSet::new(&s, &p, &[Constraint::Barycenter { component: None, target: None }], BarycenterMode::PerComponent).unwrap();
    let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
    let w = SobolevWeights::default();
    let x = ConstraintSet::new(&s, &p, &[], BarycenterMode::PerComponent);
    assert!(x.is_err(), "translations must be fixed");
    let g = random_vec(3 * s.n_vertices(), 3);
    let mut sols = Vec::new();
    for (backend, label) in [
        (Backend::Auto { dense_max_faces: 320, hierarchical: HierarchicalOptions::default() }, "dense LU (per axis)"),
        (Backend::Auto { dense_max_faces: 319, hierarchical: HierarchicalOptions::default() }, "gmres (hierarchical)"),
        (Backend::Dense, "dense LU (per axis)"),
        (Backend::Hierarchical(HierarchicalOptions::default()), "gmres (hierarchical)"),
        (Backend::default(), "dense LU (per axis)"),
    ] {
        let mut solver = saddle_solver(pool.clone(), &s, &p, &e, &w, &MetricExtras::default(), &cq, &backend).unwrap();
        assert_eq!(solver.base_report().solver, label, "{backend:?}");
        let sol = solver.solve(&g, &vec![0.0; cq.nrows()]).unwrap();
        sols.push(sol.x);
    }
    let err = rel(&sols[1], &sols[0]);
    eprintln!("hierarchical (chi 0.5, tol 1e-4) vs dense direction: {err:.1e}");
    assert!(err < 5e-2);
}

/// The preconditioner on the dense operator gives few iterations (vertex pairs) on a mesh where
/// the face form needs many (WP2's measurement, here with the sparse preconditioner).
#[test]
fn vertex_pairs_precondition_better_than_face_pairs() {
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let e = Exponents::default();
    let (s, p) = perturbed(icosphere(3).unwrap(), 0.1, 2);
    let set = ConstraintSet::new(&s, &p, &[Constraint::Barycenter { component: None, target: None }], BarycenterMode::PerComponent).unwrap();
    let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
    let geom = FaceGeometry::new(&s, &p).unwrap();
    let (_, g) = energy_and_differential(&pool, &s, &p, &geom, &e).unwrap();
    let g = flat(&g);
    let mut iters = Vec::new();
    for pairs in [FractionalForm::VertexPairs, FractionalForm::FacePairs] {
        let mut ho = HierarchicalOptions::default();
        ho.preconditioner.fractional = pairs;
        let mut it = hierarchical_saddle(pool.clone(), &s, &p, &e, &SobolevWeights::default(), &MetricExtras::default(), &cq, &ho, None).unwrap();
        let sol = it.solve(&g, &vec![0.0; cq.nrows()]).unwrap();
        eprintln!("icosphere 1280 {pairs:?}: {} iterations, residual {:.1e}", sol.record.iterations, sol.record.relative_residual);
        iters.push(sol.record.iterations);
    }
    assert!(iters[0] <= 15, "{iters:?}");
    assert!(iters[0] < iters[1], "{iters:?}");
}

/// The SPEC's congruence form `L′⁻¹ blockdiag(L, 0) L′⁻¹` has no multiplier component in its
/// range: right GMRES stalls on a right side with a net force (here random) and with pins (their
/// reaction forces), left GMRES "converges" with a large true residual. The augmented form
/// (default) converges in both (proposed D33).
#[test]
fn augmented_form_reaches_the_multipliers() {
    use ddg_repulsive_surfaces::solve::{GmresSide, gmres};
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let e = Exponents::default();
    let (s1, p1) = perturbed(icosphere(2).unwrap(), 0.2, 1);
    let (s2, p2) = perturbed(torus(14, 8, 1.5, 0.6).unwrap(), 0.1, 2);
    let cases = [
        ("icosphere 320, barycenter, random force", s1, p1, vec![Constraint::Barycenter { component: None, target: None }], true),
        ("torus 224, 3 pins, energy gradient", s2, p2, vec![Constraint::Pins { vertices: vec![0, 17, 33], targets: None }], false),
    ];
    for (name, s, p, cons, random) in cases {
        let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
        let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
        let mut rhs = if random {
            random_vec(3 * s.n_vertices(), 3)
        } else {
            let geom = FaceGeometry::new(&s, &p).unwrap();
            flat(&energy_and_differential(&pool, &s, &p, &geom, &e).unwrap().1)
        };
        rhs.extend(vec![0.0; cq.nrows()]);
        for form in [PreconditionerForm::Augmented, PreconditionerForm::Congruence] {
            for side in [GmresSide::Right, GmresSide::Left] {
                // The solver defaults this test was written for (tolerance 1e-4, restart 30).
                let mut ho = HierarchicalOptions {
                    gmres: ddg_repulsive_surfaces::solve::GmresOptions::default(),
                    ..HierarchicalOptions::default()
                };
                ho.preconditioner.form = form;
                // The exact preconditioner (refined L' solves): the D33 property of the forms.
                ho.preconditioner.refine = true;
                ho.gmres.side = side;
                ho.gmres.max_iterations = 300;
                let it = hierarchical_saddle(pool.clone(), &s, &p, &e, &SobolevWeights::default(), &MetricExtras::default(), &cq, &ho, None).unwrap();
                let (res, its) = match gmres(it.operator(), it.preconditioner(), &rhs, &ho.gmres) {
                    Ok((_, rep)) => (rep.relative_residual, rep.iterations),
                    Err(ddg_repulsive_surfaces::RsError::SolverDiverged { residual }) => (residual, 300),
                    Err(e) => panic!("{e}"),
                };
                eprintln!("{name}, {form:?} {side:?}: {its} iterations, true residual {res:.1e}");
                match form {
                    // Left GMRES stops on the preconditioned residual, whose weighting of the
                    // position and constraint parts follows the fractional scale (SOLVERFIX:
                    // 6.5e-3 with `FractionalScale::Symbol`, right GMRES 5.8e-5).
                    PreconditionerForm::Augmented => {
                        let tol = if side == GmresSide::Left { 1e-2 } else { 5e-3 };
                        assert!(its <= 30 && res < tol, "{name} {side:?}: {its} {res:e}")
                    }
                    PreconditionerForm::Congruence => assert!(res > 1e-2, "{name} {side:?}: {res:e}"),
                }
            }
        }
    }
}
