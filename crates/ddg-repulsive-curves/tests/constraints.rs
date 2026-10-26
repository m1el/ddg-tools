//! Tests of the constraints and the saddle solves (`SPEC.md` §5.1–5.2, §6.1–6.3, §13.4–13.5):
//! Jacobians by central differences, the D9 redundancy patterns, rank and translation errors,
//! the projected gradient (`C g̃ = 0`, the `Ā`-closest such direction), the chord projection,
//! the projector identities with Eq. 30 and 33, determinism and scaling.

#![allow(clippy::needless_range_loop)]

use std::sync::Arc;

use ddg_core::sparse::SparseMatrix;
use ddg_repulsive_curves::constraints::{
    Constraint, ConstraintSet, ImplicitSurface, RowKind, TangentSite, default_constraints,
    tangent_frame,
};
use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::metric::{DofLayout, GramQuadrature, Metric, sobolev_gram};
use ddg_repulsive_curves::solve::{
    DenseSaddle, ProjectionOptions, ProjectionStatus, Projector, SaddleSolver, Tolerance,
    classic_solver, dense_sobolev_solver, project_constraints, projected_gradient,
};
use ddg_repulsive_curves::{Axis, CurveNetwork, Exponents, RcError};

type P = [f64; 3];

// ---------------------------------------------------------------- fixtures

struct Sphere {
    c: P,
    r: f64,
}

impl ImplicitSurface for Sphere {
    /// `|x − c|² − r²` (deliberately not a signed distance).
    fn value(&self, x: P) -> f64 {
        let d = [x[0] - self.c[0], x[1] - self.c[1], x[2] - self.c[2]];
        d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - self.r * self.r
    }
    fn gradient(&self, x: P) -> P {
        [
            2.0 * (x[0] - self.c[0]),
            2.0 * (x[1] - self.c[1]),
            2.0 * (x[2] - self.c[2]),
        ]
    }
}

fn closed_edges(n: usize, offset: u32) -> Vec<[u32; 2]> {
    (0..n as u32)
        .map(|k| [offset + k, offset + (k + 1) % n as u32])
        .collect()
}

fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [
                t.sin() + 2.0 * (2.0 * t).sin(),
                t.cos() - 2.0 * (2.0 * t).cos(),
                -(3.0 * t).sin(),
            ]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n, 0)).unwrap(), pos)
}

/// Open helix-like curve; the last edge is reversed so a degree-1 vertex is an `i₂`.
fn open_curve(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = 3.0 * k as f64 / n as f64;
            [t.cos(), t.sin(), 0.4 * t + 0.1 * (2.0 * t).sin()]
        })
        .collect();
    let mut edges: Vec<[u32; 2]> = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    let last = edges.len() - 1;
    edges[last] = [n as u32 - 1, n as u32 - 2];
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

/// A planar open curve in the plane z = 0.25.
fn planar_open(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = 4.0 * k as f64 / n as f64;
            [t.cos() * (1.0 + 0.2 * t), t.sin() * (1.0 + 0.2 * t), 0.25]
        })
        .collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

/// A planar closed flower in z = 0.5.
fn planar_flower(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            let r = 1.0 + 0.3 * (3.0 * t).cos();
            [r * t.cos(), r * t.sin(), 0.5]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n, 0)).unwrap(), pos)
}

fn tripod(m: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos = vec![[0.0, 0.0, 0.0]];
    let mut edges = Vec::new();
    let dirs = [[1.0, 0.0, 0.2], [-0.5, 0.9, -0.1], [-0.5, -0.8, 0.3]];
    for (a, d) in dirs.iter().enumerate() {
        let mut prev = 0u32;
        for k in 1..=m {
            let s = k as f64 / m as f64;
            pos.push([
                d[0] * s + 0.3 * s * s * (a as f64 + 1.0).cos(),
                d[1] * s + 0.3 * s * s * (a as f64 + 1.0).sin(),
                d[2] * s + 0.2 * (3.0 * s + a as f64).sin(),
            ]);
            let id = pos.len() as u32 - 1;
            edges.push([prev, id]);
            prev = id;
        }
    }
    (CurveNetwork::new(pos.len(), edges).unwrap(), pos)
}

fn hopf(n: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos = Vec::new();
    for k in 0..n {
        let t = std::f64::consts::TAU * k as f64 / n as f64;
        pos.push([t.cos(), t.sin(), 0.05 * (2.0 * t).sin()]);
    }
    for k in 0..n {
        let t = std::f64::consts::TAU * k as f64 / n as f64;
        pos.push([1.0 + t.cos(), 0.03 * (3.0 * t).cos(), t.sin()]);
    }
    let mut edges = closed_edges(n, 0);
    edges.extend(closed_edges(n, n as u32));
    (CurveNetwork::new(2 * n, edges).unwrap(), pos)
}

fn rand_vec(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        })
        .collect()
}

fn dense(c: &SparseMatrix) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; c.ncols()]; c.nrows()];
    for (r, col, v) in c.triplets() {
        m[r as usize][col as usize] += v;
    }
    m
}

fn mat_vec(c: &SparseMatrix, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; c.nrows()];
    c.mul_vec(x, &mut y).unwrap();
    y
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn dotv(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn mean_len(net: &CurveNetwork, pos: &[P]) -> f64 {
    net.mean_edge_length(pos)
}

/// `Ā x` for the dense metric of `pos`.
fn apply_metric(a: &ddg_repulsive_curves::metric::DenseMatrix, d: usize, x: &[f64]) -> Vec<f64> {
    let n = a.n();
    let mut y = vec![0.0; n * d];
    for v in 0..n {
        for w in 0..n {
            for s in 0..d {
                y[v * d + s] += a.get(v, w) * x[w * d + s];
            }
        }
    }
    y
}

// ---------------------------------------------------------------- Jacobians

fn all_constraints(net: &CurveNetwork, pos: &[P], planar: bool) -> Vec<Constraint> {
    let n = net.n_vertices() as u32;
    let sphere = Arc::new(Sphere {
        c: [0.3, -0.2, if planar { 0.25 } else { 0.1 }],
        r: 1.3,
    });
    let _ = pos;
    vec![
        Constraint::Barycenter { target: None },
        Constraint::TotalLength { target: None },
        Constraint::EdgeLengths {
            edges: Some(vec![0, 3, 5]),
            targets: None,
        },
        Constraint::Point {
            vertex: 4,
            target: None,
        },
        Constraint::Surface {
            vertices: Some(vec![2, 7]),
            surface: sphere,
        },
        Constraint::Tangent {
            site: TangentSite::Edge(6),
            direction: None,
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(0),
            direction: None,
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(9),
            direction: None,
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(n - 1),
            direction: None,
        },
    ]
}

fn check_jacobian(net: &CurveNetwork, pos: &[P], planar: Option<Axis>) {
    let cons = all_constraints(net, pos, planar.is_some());
    let set = ConstraintSet::new(net, pos, &cons, planar).unwrap();
    let layout = set.layout();
    let c = dense(&set.jacobian(net, pos).unwrap());
    let k = set.n_active_rows();
    assert_eq!(c.len(), k);
    let h = 1e-6 * mean_len(net, pos);
    let mut fd = vec![vec![0.0; layout.n_dofs(net.n_vertices())]; k];
    let x0 = layout.gather(pos);
    for j in 0..x0.len() {
        let mut p = pos.to_vec();
        let mut e = vec![0.0; x0.len()];
        e[j] = h;
        layout.add_scaled(&mut p, 1.0, &e);
        let fp = set.values(net, &p).unwrap();
        let mut m = pos.to_vec();
        layout.add_scaled(&mut m, -1.0, &e);
        let fm = set.values(net, &m).unwrap();
        for r in 0..k {
            fd[r][j] = (fp[r] - fm[r]) / (2.0 * h);
        }
    }
    for r in 0..k {
        let scale = c[r].iter().fold(0.0f64, |a, x| a.max(x.abs()));
        let err = c[r]
            .iter()
            .zip(&fd[r])
            .fold(0.0f64, |a, (x, y)| a.max((x - y).abs()));
        let row = set.rows()[set.active_rows()[r]];
        assert!(
            err <= 1e-6 * scale,
            "row {r} ({:?}): error {err:e}, scale {scale:e}",
            row.kind
        );
    }
}

#[test]
fn jacobians_match_central_differences() {
    let (net, pos) = open_curve(14);
    check_jacobian(&net, &pos, None);
    let (net, pos) = planar_open(14);
    check_jacobian(&net, &pos, Some(Axis::Z));
}

#[test]
fn row_layout_and_planar_rows() {
    let (net, pos) = planar_open(14);
    let cons = all_constraints(&net, &pos, true);
    let s3 = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    let s2 = ConstraintSet::new(&net, &pos, &cons, Some(Axis::Z)).unwrap();
    // 3 + 1 + 3 + 3 + 2 + 4 · 2 rows in 3-D; planar: 2 + 1 + 3 + 2 + 2 + 4 · 1.
    assert_eq!(s3.n_active_rows(), 20);
    assert_eq!(s2.n_active_rows(), 14);
    assert_eq!(s2.jacobian(&net, &pos).unwrap().ncols(), 2 * 14);
    // Rows come in list order.
    let kinds: Vec<RowKind> = s2.rows().iter().map(|r| r.kind).collect();
    assert_eq!(kinds[0], RowKind::Barycenter { coord: 0 });
    assert_eq!(kinds[1], RowKind::Barycenter { coord: 1 });
    assert_eq!(kinds[2], RowKind::TotalLength);
    // Every constraint holds at the start (targets are the current values; the sphere is not
    // through the curve, so only the surface rows are violated).
    let tol = s3
        .row_tolerances(
            &net,
            &pos,
            Tolerance::MeanEdgeLength(1e-4),
            mean_len(&net, &pos),
        )
        .unwrap();
    let phi = s3.values(&net, &pos).unwrap();
    for (r, (&p, t)) in phi.iter().zip(&tol).enumerate() {
        let kind = s3.rows()[s3.active_rows()[r]].kind;
        if !matches!(kind, RowKind::Surface { .. }) {
            assert!(p.abs() <= *t, "row {r} {kind:?}: {p}");
        }
    }
}

#[test]
fn tangent_frame_is_orthonormal() {
    for (i, x) in [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 2.0, 3.0],
        [-0.3, 0.1, -2.0],
        [1e-9, 1.0, 0.0],
    ]
    .iter()
    .enumerate()
    {
        let f = tangent_frame(*x);
        let n = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        let u = [x[0] / n, x[1] / n, x[2] / n];
        let d = |a: P, b: P| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        assert!((d(f[0], f[0]) - 1.0).abs() < 1e-14, "case {i}");
        assert!((d(f[1], f[1]) - 1.0).abs() < 1e-14);
        assert!(d(f[0], f[1]).abs() < 1e-14);
        assert!(d(f[0], u).abs() < 1e-14);
        assert!(d(f[1], u).abs() < 1e-14);
    }
}

// ---------------------------------------------------------------- rank and redundancy

#[test]
fn d9_redundancies_are_dropped_and_reported() {
    // An edge with both endpoints pinned: its length row is dropped.
    let (net, pos) = open_curve(12);
    let cons = vec![
        Constraint::Barycenter { target: None },
        Constraint::EdgeLengths {
            edges: None,
            targets: None,
        },
        Constraint::Point {
            vertex: 3,
            target: None,
        },
        Constraint::Point {
            vertex: 4,
            target: None,
        },
    ];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    assert_eq!(set.dropped_rows(), &[3 + 3]); // barycenter rows 0..3, edge 3 is row 6
    assert_eq!(set.rows()[6].kind, RowKind::EdgeLength { edge: 3, slot: 3 });
    let report = set.analyze(&net, &pos).unwrap();
    assert_eq!(report.dropped_rows, vec![6]);
    assert_eq!(report.n_active_rows, 3 + 11 - 1 + 6);
    assert!(report.min_scaled_eigenvalue > 1e-6);
    let l0 = mean_len(&net, &pos);
    set.check_feasible(&net, &pos, Tolerance::MeanEdgeLength(1e-4), l0)
        .unwrap();
    // The dropped row is still checked for consistency with the point targets.
    let mut targets: Vec<f64> = net.edge_geometry(&pos).length;
    targets[3] *= 1.01;
    let mut cons2 = cons.clone();
    cons2[1] = Constraint::EdgeLengths {
        edges: None,
        targets: Some(targets),
    };
    let set2 = ConstraintSet::new(&net, &pos, &cons2, None).unwrap();
    assert!(matches!(
        set2.check_feasible(&net, &pos, Tolerance::MeanEdgeLength(1e-4), l0),
        Err(RcError::InfeasibleConstraint { constraint: 1, .. })
    ));

    // Every vertex pinned: the barycenter rows (and all edge rows) are dropped.
    let (net, pos) = trefoil(7);
    let mut cons = default_constraints();
    for v in 0..7 {
        cons.push(Constraint::Point {
            vertex: v,
            target: None,
        });
    }
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    let dropped: Vec<usize> = (0..3 + 7).collect();
    assert_eq!(set.dropped_rows(), dropped.as_slice());
    let report = set.analyze(&net, &pos).unwrap();
    assert_eq!(report.n_active_rows, 21);
}

fn redundant_rows(r: Result<ddg_repulsive_curves::constraints::RankReport, RcError>) -> Vec<usize> {
    match r {
        Err(RcError::RedundantConstraints { rows }) => rows,
        other => panic!("expected RedundantConstraints, got {other:?}"),
    }
}

#[test]
fn other_rank_deficiencies_are_errors() {
    let (net, pos) = trefoil(10);
    // Total length = Σ edge lengths.
    let mut cons = default_constraints();
    cons.push(Constraint::TotalLength { target: None });
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    let rows = redundant_rows(set.analyze(&net, &pos));
    let expected: Vec<usize> = (3..14).collect(); // the 10 edge rows and the total-length row
    assert_eq!(rows, expected);

    // A surface row at a pinned vertex.
    let sphere = Arc::new(Sphere {
        c: [0.0, 0.0, 0.0],
        r: 1.0,
    });
    let cons = vec![
        Constraint::Barycenter { target: None },
        Constraint::Point {
            vertex: 2,
            target: None,
        },
        Constraint::Surface {
            vertices: Some(vec![2]),
            surface: sphere,
        },
    ];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    let rows = redundant_rows(set.analyze(&net, &pos));
    assert!(rows.contains(&6), "{rows:?}");
    assert!(rows.iter().all(|&r| (3..7).contains(&r)), "{rows:?}");

    // A vertex pinned twice.
    let cons = vec![
        Constraint::Point {
            vertex: 1,
            target: None,
        },
        Constraint::Point {
            vertex: 1,
            target: None,
        },
    ];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    assert_eq!(
        redundant_rows(set.analyze(&net, &pos)),
        vec![0, 1, 2, 3, 4, 5]
    );

    // A tangent edge with both endpoints pinned.
    let cons = vec![
        Constraint::Point {
            vertex: 4,
            target: None,
        },
        Constraint::Point {
            vertex: 5,
            target: None,
        },
        Constraint::Tangent {
            site: TangentSite::Edge(4),
            direction: None,
        },
    ];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    let rows = redundant_rows(set.analyze(&net, &pos));
    assert!(rows.contains(&6) || rows.contains(&7), "{rows:?}");
}

#[test]
fn translations_must_be_fixed() {
    let (net, pos) = trefoil(10);
    let only_lengths = vec![Constraint::EdgeLengths {
        edges: None,
        targets: None,
    }];
    let set = ConstraintSet::new(&net, &pos, &only_lengths, None).unwrap();
    assert_eq!(set.analyze(&net, &pos), Err(RcError::TranslationNotFixed));
    let none: Vec<Constraint> = Vec::new();
    let set = ConstraintSet::new(&net, &pos, &none, None).unwrap();
    assert_eq!(set.analyze(&net, &pos), Err(RcError::TranslationNotFixed));
    // One surface row fixes one direction only.
    let sphere = Arc::new(Sphere {
        c: [0.0, 0.0, 0.0],
        r: 2.0,
    });
    let cons = vec![
        Constraint::TotalLength { target: None },
        Constraint::Surface {
            vertices: Some(vec![0]),
            surface: sphere.clone(),
        },
    ];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    assert_eq!(set.analyze(&net, &pos), Err(RcError::TranslationNotFixed));
    // Three surface rows with independent normals fix all three.
    let cons = vec![Constraint::Surface {
        vertices: Some(vec![0, 3, 6]),
        surface: sphere,
    }];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    set.analyze(&net, &pos).unwrap();
    // A single point fixes them.
    let cons = vec![Constraint::Point {
        vertex: 0,
        target: None,
    }];
    let set = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
    set.analyze(&net, &pos).unwrap();
    // Planar: two translations.
    let (net, pos) = planar_flower(11);
    let set = ConstraintSet::new(&net, &pos, &default_constraints(), Some(Axis::Z)).unwrap();
    set.analyze(&net, &pos).unwrap();
}

#[test]
fn validation_errors() {
    let (net, pos) = tripod(4);
    let bad = |c: Constraint| ConstraintSet::new(&net, &pos, &[c], None).unwrap_err();
    assert_eq!(
        bad(Constraint::Tangent {
            site: TangentSite::Vertex(0),
            direction: None
        }),
        RcError::InvalidTangentSite {
            vertex: 0,
            degree: 3
        }
    );
    assert_eq!(
        bad(Constraint::EdgeLengths {
            edges: Some(vec![1, 99]),
            targets: None
        }),
        RcError::ConstraintOutOfRange {
            constraint: 0,
            index: 99
        }
    );
    assert_eq!(
        bad(Constraint::Point {
            vertex: 13,
            target: None
        }),
        RcError::ConstraintOutOfRange {
            constraint: 0,
            index: 13
        }
    );
    assert!(matches!(
        bad(Constraint::EdgeLengths {
            edges: Some(vec![1, 2]),
            targets: Some(vec![1.0])
        }),
        RcError::InvalidOption { .. }
    ));
    assert!(matches!(
        bad(Constraint::Tangent {
            site: TangentSite::Edge(0),
            direction: Some([0.0, 0.0, 0.0])
        }),
        RcError::InvalidOption { .. }
    ));
    // Antiparallel direction: infeasible (s·X ≤ 0).
    let t = net.edge_geometry(&pos).tangent[2];
    assert!(matches!(
        bad(Constraint::Tangent {
            site: TangentSite::Edge(2),
            direction: Some([-t[0], -t[1], -t[2]])
        }),
        RcError::InfeasibleConstraint { constraint: 0, .. }
    ));
    // A point target off its vertex is caught by the feasibility check.
    let l0 = mean_len(&net, &pos);
    let mut p = pos[5];
    p[1] += 1e-2 * l0;
    let set = ConstraintSet::new(
        &net,
        &pos,
        &[
            Constraint::Barycenter { target: None },
            Constraint::Point {
                vertex: 5,
                target: Some(p),
            },
        ],
        None,
    )
    .unwrap();
    assert!(matches!(
        set.check_feasible(&net, &pos, Tolerance::MeanEdgeLength(1e-4), l0),
        Err(RcError::InfeasibleConstraint { constraint: 1, .. })
    ));
    // Planar: a tangent direction off the plane.
    let (net, pos) = planar_open(10);
    assert!(matches!(
        ConstraintSet::new(
            &net,
            &pos,
            &[Constraint::Tangent {
                site: TangentSite::Edge(1),
                direction: Some([1.0, 0.0, 0.1])
            }],
            Some(Axis::Z)
        ),
        Err(RcError::InvalidOption { .. })
    ));
    // Setters.
    let mut set = ConstraintSet::new(&net, &pos, &default_constraints(), None).unwrap();
    assert!(set.set_total_length_target(1, 3.0).is_err());
    assert!(set.set_edge_length_targets(1, &[1.0]).is_err());
    assert!(set.set_edge_length_targets(1, &[0.5; 9]).is_ok());
    assert!(matches!(
        set.set_edge_length_targets(7, &[1.0]),
        Err(RcError::ConstraintOutOfRange { .. })
    ));
}

// ---------------------------------------------------------------- saddle solves

struct Problem {
    net: CurveNetwork,
    pos: Vec<P>,
    set: ConstraintSet,
    diff: Vec<P>,
}

fn problem(net: CurveNetwork, pos: Vec<P>, cons: &[Constraint], planar: Option<Axis>) -> Problem {
    let pool = make_pool(Some(4)).unwrap();
    let (_, diff) = energy_and_differential(&pool, &net, &pos, Exponents::default());
    let set = ConstraintSet::new(&net, &pos, cons, planar).unwrap();
    set.analyze(&net, &pos).unwrap();
    Problem {
        net,
        pos,
        set,
        diff,
    }
}

fn dense_solver(p: &Problem, threads: usize) -> DenseSaddle {
    dense_sobolev_solver(
        Arc::new(make_pool(Some(threads)).unwrap()),
        &p.net,
        &p.pos,
        Exponents::default(),
        GramQuadrature::Midpoint,
        &p.set,
    )
    .unwrap()
}

#[test]
fn projected_gradient_is_tangent_and_closest() {
    let cases = [
        problem(trefoil(40).0, trefoil(40).1, &default_constraints(), None),
        problem(
            open_curve(30).0,
            open_curve(30).1,
            &[
                Constraint::Barycenter { target: None },
                Constraint::EdgeLengths {
                    edges: None,
                    targets: None,
                },
                Constraint::Tangent {
                    site: TangentSite::Vertex(0),
                    direction: None,
                },
            ],
            None,
        ),
        problem(hopf(15).0, hopf(15).1, &default_constraints(), None),
        problem(
            planar_flower(31).0,
            planar_flower(31).1,
            &default_constraints(),
            Some(Axis::Z),
        ),
    ];
    for (ci, p) in cases.iter().enumerate() {
        let mut s = dense_solver(p, 4);
        let layout = p.set.layout();
        let g = projected_gradient(&mut s, layout, &p.diff).unwrap();
        assert!(
            g.linear.relative_residual < 1e-10,
            "case {ci}: {:?}",
            g.linear
        );
        let x = layout.gather(&g.direction);
        let df = layout.gather(&p.diff);
        let c = p.set.jacobian(&p.net, &p.pos).unwrap();
        // C g̃ = 0.
        let cx = mat_vec(&c, &x);
        let cnorm = dense(&c).iter().map(|r| norm(r)).fold(0.0f64, f64::max);
        assert!(
            norm(&cx) <= 1e-10 * cnorm * norm(&x),
            "case {ci}: |C g| = {:e}",
            norm(&cx)
        );
        // Planar: no component along the held axis.
        if layout.planar().is_some() {
            assert!(g.direction.iter().all(|v| v[2] == 0.0));
        }
        // Descent: dF·g̃ = g̃ᵀĀg̃ > 0.
        let a = s.metric().clone();
        let ax = apply_metric(&a, layout.dim(), &x);
        let quad = dotv(&x, &ax);
        assert!(g.sobolev_norm_sq > 0.0);
        assert!(
            (g.sobolev_norm_sq - quad).abs() <= 1e-9 * quad,
            "case {ci}: {} vs {quad}",
            g.sobolev_norm_sq
        );
        assert!((g.gradient_norm - norm(&x)).abs() <= 1e-15 * norm(&x));
        // Ā-closest: Āg̃ − dF ⟂ ker C, and J(h) = ½hᵀĀh − dF·h is minimal at g̃ on ker C.
        let proj = Projector::new(&c, None).unwrap();
        let resid: Vec<f64> = ax.iter().zip(&df).map(|(p, q)| p - q).collect();
        let jg = 0.5 * quad - dotv(&df, &x);
        for seed in 0..4 {
            let w = proj.project(&rand_vec(x.len(), 7 + seed)).unwrap();
            assert!(norm(&mat_vec(&c, &w)) <= 1e-12 * cnorm * norm(&w));
            let o = dotv(&w, &resid);
            assert!(
                o.abs() <= 1e-9 * norm(&w) * norm(&df),
                "case {ci}: orthogonality {o:e}"
            );
            let wscaled: Vec<f64> = w.iter().map(|v| v * 1e-2 * norm(&x) / norm(&w)).collect();
            let h: Vec<f64> = x.iter().zip(&wscaled).map(|(a, b)| a + b).collect();
            let jh = 0.5 * dotv(&h, &apply_metric(&a, layout.dim(), &h)) - dotv(&df, &h);
            assert!(jh > jg, "case {ci}: J(g̃ + w) = {jh} ≤ J(g̃) = {jg}");
        }
    }
}

#[test]
fn classic_metric_gradients() {
    let p = problem(trefoil(40).0, trefoil(40).1, &default_constraints(), None);
    let c = p.set.jacobian(&p.net, &p.pos).unwrap();
    for metric in [Metric::L2, Metric::H1, Metric::H2] {
        let mut s = classic_solver(&p.net, &p.pos, &metric, &p.set).unwrap();
        let g = projected_gradient(&mut s, p.set.layout(), &p.diff).unwrap();
        assert!(
            g.linear.relative_residual < 1e-10,
            "{metric:?}: {:?}",
            g.linear
        );
        let x = p.set.layout().gather(&g.direction);
        assert!(norm(&mat_vec(&c, &x)) <= 1e-10 * norm(&x), "{metric:?}");
        assert!(g.sobolev_norm_sq > 0.0);
    }
    // Two components, translations fixed only globally: L-type metrics are singular, L² and the
    // Sobolev metric are not (D5).
    let p = problem(hopf(12).0, hopf(12).1, &default_constraints(), None);
    for metric in [Metric::H1, Metric::H2] {
        assert_eq!(
            classic_solver(&p.net, &p.pos, &metric, &p.set).unwrap_err(),
            RcError::SingularMetric { components: 2 }
        );
    }
    classic_solver(&p.net, &p.pos, &Metric::L2, &p.set).unwrap();
    dense_solver(&p, 2);
    // Pinning a vertex of the second component fixes its translation.
    let mut cons = default_constraints();
    cons.push(Constraint::Point {
        vertex: 15,
        target: None,
    });
    let p = problem(hopf(12).0, hopf(12).1, &cons, None);
    classic_solver(&p.net, &p.pos, &Metric::H2, &p.set).unwrap();
    assert!(classic_solver(&p.net, &p.pos, &Metric::default(), &p.set).is_err());
}

#[test]
fn collinear_network_is_a_singular_metric() {
    let n = 8;
    let pos: Vec<P> = (0..n).map(|k| [k as f64 * 0.7, 0.0, 0.0]).collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    let net = CurveNetwork::new(n, edges).unwrap();
    let set = ConstraintSet::new(&net, &pos, &default_constraints(), None).unwrap();
    set.analyze(&net, &pos).unwrap();
    let err = dense_sobolev_solver(
        Arc::new(make_pool(Some(1)).unwrap()),
        &net,
        &pos,
        Exponents::default(),
        GramQuadrature::Midpoint,
        &set,
    )
    .unwrap_err();
    assert_eq!(err, RcError::SingularMetric { components: 1 });
}

#[test]
fn dense_solve_is_bitwise_independent_of_thread_count() {
    let p = problem(trefoil(400).0, trefoil(400).1, &default_constraints(), None);
    let layout = p.set.layout();
    let reference = {
        let mut s = dense_solver(&p, 1);
        projected_gradient(&mut s, layout, &p.diff).unwrap()
    };
    for t in [7, 64] {
        let mut s = dense_solver(&p, t);
        let g = projected_gradient(&mut s, layout, &p.diff).unwrap();
        assert!(
            g.direction == reference.direction && g.lambda == reference.lambda,
            "{t} threads"
        );
    }
}

#[test]
fn projected_direction_scales_with_the_curve() {
    let e = Exponents::default();
    let (net, pos) = trefoil(30);
    let pool = make_pool(Some(2)).unwrap();
    let dir = |pos: &[P]| {
        let (_, diff) = energy_and_differential(&pool, &net, pos, e);
        let set = ConstraintSet::new(&net, pos, &default_constraints(), None).unwrap();
        let mut s = dense_sobolev_solver(
            Arc::new(make_pool(Some(2)).unwrap()),
            &net,
            pos,
            e,
            GramQuadrature::Midpoint,
            &set,
        )
        .unwrap();
        projected_gradient(&mut s, set.layout(), &diff)
            .unwrap()
            .direction
    };
    let g = dir(&pos);
    let sigma = e.sigma();
    for c in [1e-3, 0.7, 1e3] {
        let scaled: Vec<P> = pos.iter().map(|p| [c * p[0], c * p[1], c * p[2]]).collect();
        let gs = dir(&scaled);
        let f = c.powf(2.0 * sigma + 2.0 + e.alpha() - e.beta());
        let gmax = g.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()));
        for (a, b) in g.iter().flatten().zip(gs.iter().flatten()) {
            assert!((a * f - b).abs() <= 1e-9 * gmax * f, "c = {c}");
        }
    }
}

// ---------------------------------------------------------------- projection

fn perturbed_projection(p: &Problem, trial: Vec<P>, what: &str) {
    let mut s = dense_solver(p, 4);
    let l0 = mean_len(&p.net, &p.pos);
    let opts = ProjectionOptions::default();
    let r = project_constraints(&mut s, &p.set, &p.net, &trial, &opts, l0).unwrap();
    assert_eq!(r.status, ProjectionStatus::Converged, "{what}: {r:?}");
    assert!(
        r.iterations >= 1 && r.iterations <= 3,
        "{what}: {} iterations",
        r.iterations
    );
    assert!(*r.history.last().unwrap() <= 1.0);
    let tol = p
        .set
        .row_tolerances(&p.net, &r.positions, opts.tolerance, l0)
        .unwrap();
    let phi = p.set.values(&p.net, &r.positions).unwrap();
    for (x, t) in phi.iter().zip(&tol) {
        assert!(x.abs() <= *t);
    }
}

#[test]
fn projection_converges_in_at_most_three_chord_iterations() {
    // A descent step on a trefoil, barycenter + edge lengths.
    let p = problem(trefoil(64).0, trefoil(64).1, &default_constraints(), None);
    let l0 = mean_len(&p.net, &p.pos);
    let mut s = dense_solver(&p, 4);
    let g = projected_gradient(&mut s, p.set.layout(), &p.diff).unwrap();
    let gmax = g
        .direction
        .iter()
        .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
        .fold(0.0f64, f64::max);
    for frac in [0.05, 0.2] {
        let tau = frac * l0 / gmax;
        let trial: Vec<P> = p
            .pos
            .iter()
            .zip(&g.direction)
            .map(|(x, d)| [x[0] - tau * d[0], x[1] - tau * d[1], x[2] - tau * d[2]])
            .collect();
        perturbed_projection(&p, trial, &format!("descent step {frac} l0"));
    }
    // A random perturbation of 1 % of ℓ̄₀.
    let noise = rand_vec(3 * 64, 99);
    let trial: Vec<P> = p
        .pos
        .iter()
        .enumerate()
        .map(|(v, x)| {
            [
                x[0] + 1e-2 * l0 * noise[3 * v],
                x[1] + 1e-2 * l0 * noise[3 * v + 1],
                x[2] + 1e-2 * l0 * noise[3 * v + 2],
            ]
        })
        .collect();
    perturbed_projection(&p, trial, "noise");

    // An open curve with every constraint type (the sphere passes through the surface vertices).
    let (net, pos) = open_curve(24);
    let centre = [0.1, 0.2, 0.3];
    let r = |v: usize| {
        let d = [
            pos[v][0] - centre[0],
            pos[v][1] - centre[1],
            pos[v][2] - centre[2],
        ];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    let cons = vec![
        Constraint::Barycenter { target: None },
        Constraint::EdgeLengths {
            edges: None,
            targets: None,
        },
        Constraint::Point {
            vertex: 5,
            target: None,
        },
        Constraint::Surface {
            vertices: Some(vec![12]),
            surface: Arc::new(Sphere {
                c: centre,
                r: r(12),
            }),
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(0),
            direction: None,
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(17),
            direction: None,
        },
    ];
    let p = problem(net, pos, &cons, None);
    let l0 = mean_len(&p.net, &p.pos);
    p.set
        .check_feasible(&p.net, &p.pos, Tolerance::MeanEdgeLength(1e-4), l0)
        .unwrap();
    let noise = rand_vec(3 * 24, 5);
    let trial: Vec<P> = p
        .pos
        .iter()
        .enumerate()
        .map(|(v, x)| {
            [
                x[0] + 1e-2 * l0 * noise[3 * v],
                x[1] + 1e-2 * l0 * noise[3 * v + 1],
                x[2] + 1e-2 * l0 * noise[3 * v + 2],
            ]
        })
        .collect();
    perturbed_projection(&p, trial, "open curve, all types");
}

#[test]
fn projection_reports_failures() {
    let p = problem(trefoil(30).0, trefoil(30).1, &default_constraints(), None);
    let l0 = mean_len(&p.net, &p.pos);
    let mut s = dense_solver(&p, 2);
    // Already feasible: zero iterations.
    let r = project_constraints(
        &mut s,
        &p.set,
        &p.net,
        &p.pos,
        &ProjectionOptions::default(),
        l0,
    )
    .unwrap();
    assert_eq!((r.status, r.iterations), (ProjectionStatus::Converged, 0));
    // A huge perturbation with a single iteration allowed cannot converge.
    let noise = rand_vec(90, 3);
    let trial: Vec<P> = p
        .pos
        .iter()
        .enumerate()
        .map(|(v, x)| {
            [
                x[0] + 0.3 * l0 * noise[3 * v],
                x[1] + 0.3 * l0 * noise[3 * v + 1],
                x[2] + 0.3 * l0 * noise[3 * v + 2],
            ]
        })
        .collect();
    let opts = ProjectionOptions {
        tolerance: Tolerance::Absolute(1e-12),
        max_iterations: 1,
    };
    let r = project_constraints(&mut s, &p.set, &p.net, &trial, &opts, l0).unwrap();
    assert_ne!(r.status, ProjectionStatus::Converged);
    assert_eq!(r.history.len(), 2);
    let bad = ProjectionOptions {
        tolerance: Tolerance::Absolute(0.0),
        max_iterations: 3,
    };
    assert!(project_constraints(&mut s, &p.set, &p.net, &trial, &bad, l0).is_err());
}

// ---------------------------------------------------------------- projector

#[test]
fn projector_identities_and_eq_30_33() {
    let p = problem(trefoil(30).0, trefoil(30).1, &default_constraints(), None);
    let c = p.set.jacobian(&p.net, &p.pos).unwrap();
    let proj = Projector::new(&c, None).unwrap();
    let nd = c.ncols();
    let k = c.nrows();
    let x = rand_vec(nd, 1);
    let y = rand_vec(nd, 2);
    let px = proj.project(&x).unwrap();
    let py = proj.project(&y).unwrap();
    // C P = 0.
    assert!(norm(&mat_vec(&c, &px)) <= 1e-12 * norm(&x));
    // P² = P.
    let ppx = proj.project(&px).unwrap();
    let d: Vec<f64> = ppx.iter().zip(&px).map(|(a, b)| a - b).collect();
    assert!(norm(&d) <= 1e-12 * norm(&x));
    // Pᵀ = P.
    assert!((dotv(&px, &y) - dotv(&x, &py)).abs() <= 1e-12 * norm(&x) * norm(&y));
    // C C† b = b, and C† b ⟂ ker C.
    let b = rand_vec(k, 3);
    let z = proj.pseudo_inverse(&b).unwrap();
    let cz = mat_vec(&c, &z);
    let d: Vec<f64> = cz.iter().zip(&b).map(|(a, q)| a - q).collect();
    assert!(norm(&d) <= 1e-12 * norm(&b));
    assert!(dotv(&z, &px).abs() <= 1e-12 * norm(&z) * norm(&px));
    // (C Cᵀ) (C Cᵀ)⁻¹ b = b.
    let s = proj.solve_gram(&b).unwrap();
    let mut cts = vec![0.0; nd];
    c.transpose().mul_vec(&s, &mut cts).unwrap();
    let back = mat_vec(&c, &cts);
    let d: Vec<f64> = back.iter().zip(&b).map(|(a, q)| a - q).collect();
    assert!(norm(&d) <= 1e-12 * norm(&b));

    // Eq. 30: x = P y with Pᵀ Ā P y = Pᵀ a solves the saddle system with b = 0. CG on the
    // projected system against the dense saddle solve.
    let mut solver = dense_solver(&p, 2);
    let a_mat = solver.metric().clone();
    let a = p.set.layout().gather(&p.diff);
    let exact = solver.solve(&a, &vec![0.0; k]).unwrap().x;
    let op = |v: &[f64]| {
        proj.projected_operator(|w| apply_metric(&a_mat, 3, w), v)
            .unwrap()
    };
    let rhs = proj.project(&a).unwrap();
    let mut yv = vec![0.0; nd];
    let mut r = rhs.clone();
    let mut dir = r.clone();
    let mut rr = dotv(&r, &r);
    for _ in 0..5000 {
        if rr.sqrt() <= 1e-14 * norm(&rhs) {
            break;
        }
        let q = op(&dir);
        let alpha = rr / dotv(&dir, &q);
        for i in 0..nd {
            yv[i] += alpha * dir[i];
            r[i] -= alpha * q[i];
        }
        let rr2 = dotv(&r, &r);
        for i in 0..nd {
            dir[i] = r[i] + rr2 / rr * dir[i];
        }
        rr = rr2;
    }
    let xcg = proj.project(&yv).unwrap();
    let d: Vec<f64> = xcg.iter().zip(&exact).map(|(a, b)| a - b).collect();
    assert!(
        norm(&d) <= 1e-7 * norm(&exact),
        "Eq. 30: {:e}",
        norm(&d) / norm(&exact)
    );

    // Eq. 33: z = C† b, y solves the b = 0 system with right side Ā z, x = z − y solves
    // [Ā Cᵀ; C 0][x; µ] = [0; b]; and the general form with a ≠ 0.
    let bb = rand_vec(k, 11);
    let z = proj.pseudo_inverse(&bb).unwrap();
    let az = apply_metric(&a_mat, 3, &z);
    let yz = solver.solve(&az, &vec![0.0; k]).unwrap().x;
    let x33: Vec<f64> = z.iter().zip(&yz).map(|(a, b)| a - b).collect();
    let direct = solver.solve(&vec![0.0; nd], &bb).unwrap().x;
    let d: Vec<f64> = x33.iter().zip(&direct).map(|(a, b)| a - b).collect();
    assert!(
        norm(&d) <= 1e-9 * norm(&direct),
        "Eq. 33: {:e}",
        norm(&d) / norm(&direct)
    );
    let mut kernel_solver = dense_solver(&p, 2);
    let general = proj
        .solve_with_kernel_solver(
            |w| apply_metric(&a_mat, 3, w),
            |r| Ok(kernel_solver.solve(r, &vec![0.0; k])?.x),
            &a,
            &bb,
        )
        .unwrap();
    let direct = solver.solve(&a, &bb).unwrap().x;
    let d: Vec<f64> = general.iter().zip(&direct).map(|(a, b)| a - b).collect();
    assert!(
        norm(&d) <= 1e-9 * norm(&direct),
        "Eq. 33 (a ≠ 0): {:e}",
        norm(&d) / norm(&direct)
    );

    // No rows: P = I.
    let empty = SparseMatrix::from_triplets(0, nd, &[]).unwrap();
    let id = Projector::new(&empty, None).unwrap();
    assert_eq!(id.project(&x).unwrap(), x);
}

#[test]
fn solver_rejects_wrong_lengths() {
    let p = problem(trefoil(12).0, trefoil(12).1, &default_constraints(), None);
    let mut s = dense_solver(&p, 1);
    assert!(s.solve(&[1.0, 2.0], &[]).is_err());
    let g = sobolev_gram(
        &make_pool(Some(1)).unwrap(),
        &p.net,
        &p.pos,
        Exponents::default(),
        GramQuadrature::Midpoint,
    );
    let c = p.set.jacobian(&p.net, &p.pos).unwrap();
    assert!(
        DenseSaddle::new(
            Arc::new(make_pool(Some(1)).unwrap()),
            &g.a,
            &c,
            DofLayout::new(Some(Axis::Y))
        )
        .is_err()
    );
}

/// Timings of the dense assembly and one saddle solve on the Freedman–He–Wang curves at 1 and
/// 64 threads, with bitwise comparison of the results. Run alone:
/// `cargo test --release -p ddg-repulsive-curves --test constraints freedman_timings -- --ignored --nocapture`.
#[test]
#[ignore]
fn freedman_timings() {
    use std::time::Instant;
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/repulsive-curves/FreedmanHeWang"
    );
    for size in [600, 1200, 2400] {
        let path = format!("{dir}/freedman-he-wang-{size}.obj");
        let Ok((net, pos)) = ddg_repulsive_curves::read_curve_obj(&path) else {
            eprintln!("TIMING skip {path} (not found)");
            continue;
        };
        let e = Exponents::default();
        let (_, diff) = energy_and_differential(&make_pool(None).unwrap(), &net, &pos, e);
        let set = ConstraintSet::new(&net, &pos, &default_constraints(), None).unwrap();
        set.analyze(&net, &pos).unwrap();
        let c = set.jacobian(&net, &pos).unwrap();
        let mut results = Vec::new();
        let threads_list: Vec<usize> = std::env::var("RC_TIMING_THREADS")
            .unwrap_or_else(|_| "1,64".into())
            .split(',')
            .filter_map(|t| t.trim().parse().ok())
            .collect();
        for &threads in &threads_list {
            let pool = Arc::new(make_pool(Some(threads)).unwrap());
            let t0 = Instant::now();
            let g = sobolev_gram(&pool, &net, &pos, e, GramQuadrature::Midpoint);
            let t_gram = t0.elapsed().as_secs_f64();
            let t1 = Instant::now();
            let mut s = DenseSaddle::new(pool.clone(), &g.a, &c, set.layout()).unwrap();
            let t_factor = t1.elapsed().as_secs_f64();
            let t2 = Instant::now();
            let pg = projected_gradient(&mut s, set.layout(), &diff).unwrap();
            let t_solve = t2.elapsed().as_secs_f64();
            eprintln!(
                "TIMING |V| = {} |E| = {} k = {} dim = {} threads = {threads}: gram {t_gram:.3} s, \
                 factor {t_factor:.3} s, solve {t_solve:.3} s (residual {:.1e}, pivot ratio {:.1e})",
                net.n_vertices(),
                net.n_edges(),
                c.nrows(),
                pg.linear.dimension,
                pg.linear.relative_residual,
                s.pivot_ratio()
            );
            results.push((g, pg.direction));
        }
        eprintln!(
            "TIMING |V| = {}: gram bitwise equal across threads: {}, direction bitwise equal: {}",
            net.n_vertices(),
            results.iter().all(|r| r.0 == results[0].0),
            results.iter().all(|r| r.1 == results[0].1)
        );
    }
}
