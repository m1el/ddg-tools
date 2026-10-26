//! The Schur-complement dense saddle solver against the LU (`SPEC.md` §17l): same solutions to
//! 1e-10 relative on the test curves (closed, open with pins and an end tangent, a network with
//! junctions, planar, few constraint rows), bitwise identical across thread counts, and the LU
//! fallback on singular saddles.

#![allow(clippy::needless_range_loop)]

use std::sync::Arc;

use ddg_core::sparse::SparseMatrix;
use ddg_repulsive_curves::constraints::{Constraint, ConstraintSet, TangentSite};
use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::metric::{GramQuadrature, sobolev_gram};
use ddg_repulsive_curves::solve::{DenseMethod, DenseSaddle, SaddleSolver};
use ddg_repulsive_curves::{Axis, CurveNetwork, Exponents, RcError, read_curve_obj};

type P = [f64; 3];

fn pool(n: usize) -> Arc<rayon::ThreadPool> {
    Arc::new(make_pool(Some(n)).unwrap())
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn rel(a: &[f64], b: &[f64]) -> f64 {
    let d: Vec<f64> = a.iter().zip(b).map(|(p, q)| p - q).collect();
    norm(&d) / norm(b).max(f64::MIN_POSITIVE)
}

fn pseudo(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        })
        .collect()
}

fn closed(pos: Vec<P>) -> (CurveNetwork, Vec<P>) {
    let n = pos.len();
    let edges = (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    closed(
        (0..n)
            .map(|k| {
                let t = std::f64::consts::TAU * k as f64 / n as f64;
                [
                    t.sin() + 2.0 * (2.0 * t).sin(),
                    t.cos() - 2.0 * (2.0 * t).cos(),
                    -(3.0 * t).sin(),
                ]
            })
            .collect(),
    )
}

/// A planar closed curve in z = 0 (a lobed loop).
fn planar_loop(n: usize) -> (CurveNetwork, Vec<P>) {
    closed(
        (0..n)
            .map(|k| {
                let t = std::f64::consts::TAU * k as f64 / n as f64;
                let r = 1.0 + 0.3 * (3.0 * t).cos();
                [r * t.cos(), r * t.sin(), 0.0]
            })
            .collect(),
    )
}

/// An open helix.
fn helix(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos: Vec<P> = (0..n)
        .map(|k| {
            let t = 4.0 * std::f64::consts::PI * k as f64 / (n - 1) as f64;
            [t.cos(), t.sin(), 0.15 * t]
        })
        .collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

/// A theta graph: two junctions joined by three arcs.
fn theta(m: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos: Vec<P> = vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let mut edges = Vec::new();
    for (arc, (y, z)) in [(0.8, 0.1), (0.0, -0.3), (-0.8, 0.2)].iter().enumerate() {
        let mut prev = 0u32;
        for k in 1..m {
            let t = k as f64 / m as f64;
            let s = (std::f64::consts::PI * t).sin();
            pos.push([-1.0 + 2.0 * t, y * s, z * s + 0.05 * arc as f64 * s * s]);
            let v = pos.len() as u32 - 1;
            edges.push([prev, v]);
            prev = v;
        }
        edges.push([prev, 1]);
    }
    (CurveNetwork::new(pos.len(), edges).unwrap(), pos)
}

struct Case {
    name: &'static str,
    net: CurveNetwork,
    pos: Vec<P>,
    ex: Exponents,
    set: ConstraintSet,
    c: SparseMatrix,
    diff: Vec<f64>,
}

fn case(
    name: &'static str,
    (net, pos): (CurveNetwork, Vec<P>),
    ex: Exponents,
    cons: Vec<Constraint>,
    planar: Option<Axis>,
) -> Case {
    let set = ConstraintSet::new(&net, &pos, &cons, planar).unwrap();
    let c = set.jacobian(&net, &pos).unwrap();
    let (_, d) = energy_and_differential(&pool(8), &net, &pos, ex);
    let diff = set.layout().gather(&d);
    Case {
        name,
        net,
        pos,
        ex,
        set,
        c,
        diff,
    }
}

fn cases() -> Vec<Case> {
    let ex = Exponents::new(3.0, 6.0).unwrap();
    let ex2 = Exponents::new(2.0, 4.5).unwrap();
    let default = ddg_repulsive_curves::default_constraints;
    let mut v = vec![
        case("trefoil 120", trefoil(120), ex, default(), None),
        case(
            "trefoil 90, barycenter + total length",
            trefoil(90),
            ex2,
            vec![
                Constraint::Barycenter { target: None },
                Constraint::TotalLength { target: None },
            ],
            None,
        ),
        case(
            "open helix 80, pinned ends, end tangent, edge lengths",
            helix(80),
            ex,
            vec![
                Constraint::Point {
                    vertex: 0,
                    target: None,
                },
                Constraint::Point {
                    vertex: 79,
                    target: None,
                },
                Constraint::Tangent {
                    site: TangentSite::Vertex(0),
                    direction: None,
                },
                Constraint::EdgeLengths {
                    edges: None,
                    targets: None,
                },
            ],
            None,
        ),
        case(
            "theta network 3 x 30",
            theta(30),
            ex,
            vec![
                Constraint::Point {
                    vertex: 0,
                    target: None,
                },
                Constraint::EdgeLengths {
                    edges: None,
                    targets: None,
                },
            ],
            None,
        ),
        case(
            "planar loop 100",
            planar_loop(100),
            ex2,
            default(),
            Some(Axis::Z),
        ),
    ];
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/repulsive-curves/FreedmanHeWang/freedman-he-wang-600.obj");
    if let Ok(c) = read_curve_obj(&path) {
        v.push(case("Freedman-He-Wang 592", c, ex2, default(), None));
    } else {
        eprintln!(
            "skipping Freedman-He-Wang 592: {} not found",
            path.display()
        );
    }
    v
}

fn solver(c: &Case, threads: usize, method: DenseMethod) -> DenseSaddle {
    let p = pool(threads);
    let g = sobolev_gram(&p, &c.net, &c.pos, c.ex, GramQuadrature::Midpoint);
    DenseSaddle::with_method(p, &g.a, &c.c, c.set.layout(), method).unwrap()
}

#[test]
fn schur_matches_lu_to_1e_10() {
    for c in cases() {
        let mut lu = solver(&c, 4, DenseMethod::Lu);
        let mut sc = solver(&c, 4, DenseMethod::Schur);
        assert_eq!(lu.method(), DenseMethod::Lu);
        assert_eq!(
            sc.method(),
            DenseMethod::Schur,
            "{}: Schur fell back",
            c.name
        );
        let k = c.c.nrows();
        let nd = sc.n_dofs();
        // Eq. 22 (a = dF, b = 0) and a projection right side (a = 0, b random), and both.
        let rhs = [
            (c.diff.clone(), vec![0.0; k]),
            (
                vec![0.0; nd],
                pseudo(k, 3).iter().map(|v| 1e-3 * v).collect(),
            ),
            (pseudo(nd, 5), pseudo(k, 7)),
        ];
        for (j, (a, b)) in rhs.iter().enumerate() {
            let x_lu = lu.solve(a, b).unwrap();
            let x_sc = sc.solve(a, b).unwrap();
            let ex = rel(&x_sc.x, &x_lu.x);
            let el = rel(&x_sc.lambda, &x_lu.lambda);
            eprintln!(
                "{} rhs {j}: |x - x_LU| {ex:.1e}, |λ - λ_LU| {el:.1e}, residual Schur {:.1e} LU {:.1e}, \
                 pivot ratio Schur {:.1e} LU {:.1e}",
                c.name,
                x_sc.report.relative_residual,
                x_lu.report.relative_residual,
                sc.pivot_ratio(),
                lu.pivot_ratio()
            );
            assert!(ex <= 1e-10, "{} rhs {j}: x differs by {ex:.2e}", c.name);
            assert!(el <= 1e-10, "{} rhs {j}: λ differs by {el:.2e}", c.name);
            assert!(
                x_sc.report.relative_residual <= 1e-12f64.max(x_lu.report.relative_residual),
                "{} rhs {j}: residual {:.2e} (LU {:.2e})",
                c.name,
                x_sc.report.relative_residual,
                x_lu.report.relative_residual
            );
            assert_eq!(x_sc.report.solver, "dense Schur");
        }
    }
}

#[test]
fn schur_is_bitwise_identical_across_threads() {
    for c in cases() {
        let a = c.diff.clone();
        let b = pseudo(c.c.nrows(), 11);
        let mut first: Option<(Vec<f64>, Vec<f64>)> = None;
        for t in [1, 3, 8, 16] {
            let mut s = solver(&c, t, DenseMethod::Schur);
            let r = s.solve(&a, &b).unwrap();
            match &first {
                None => first = Some((r.x, r.lambda)),
                Some((x, l)) => {
                    assert!(
                        r.x == *x && r.lambda == *l,
                        "{}: {t} threads differ from 1",
                        c.name
                    );
                }
            }
        }
    }
}

#[test]
fn singular_saddles_fall_back_to_the_lu_and_fail_as_before() {
    // Edge lengths only: the translations are not fixed (W singular).
    let (net, pos) = trefoil(60);
    let ex = Exponents::new(3.0, 6.0).unwrap();
    let cons = vec![Constraint::EdgeLengths {
        edges: None,
        targets: None,
    }];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    let c = set.jacobian(&net, &pos).unwrap();
    let p = pool(2);
    let g = sobolev_gram(&p, &net, &pos, ex, GramQuadrature::Midpoint);
    let lu = DenseSaddle::with_method(p.clone(), &g.a, &c, set.layout(), DenseMethod::Lu);
    let sc = DenseSaddle::with_method(p, &g.a, &c, set.layout(), DenseMethod::Schur);
    match (lu, sc) {
        (Err(e1), Err(e2)) => {
            assert!(matches!(e1, RcError::SingularSaddle { .. }));
            assert_eq!(e1, e2);
        }
        (Ok(_), Ok(s)) => {
            // Not singular to the LU's tolerance: then the Schur method must have fallen back
            // (W is singular in exact arithmetic).
            assert_eq!(s.method(), DenseMethod::Lu);
        }
        (a, b) => panic!("LU and Schur disagree: {:?} vs {:?}", a.err(), b.err()),
    }
}
