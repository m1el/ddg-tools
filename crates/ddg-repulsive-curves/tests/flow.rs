//! The flow (`SPEC.md` §7, §13.5): energy decrease, the isotopy guarantee of the
//! collision-guarded search (every accepted step audited by an independent CCD sweep, a knot
//! invariant), tunnelling of the normalized search, schedules and subdivision, planar mode,
//! every stop status, the soft methods, the session state, determinism, bad options; the
//! Freedman–He–Wang unknot (release, ignored) and the line-search comparison driver.

#![allow(clippy::needless_range_loop)]

use std::sync::Arc;
use std::sync::atomic::Ordering;
use std::time::{Duration, Instant};

use ddg_repulsive_curves::ccd::{CcdOptions, CcdReport, segment_segment_distance as seg_dist};
use ddg_repulsive_curves::optimize::{Hencky, Method};
use ddg_repulsive_curves::solve::{Projection, ProjectionStatus};
use ddg_repulsive_curves::step::line_search::{
    LineSearch, SearchEnv, SearchInput, StartRule, TrialOutcome, line_search,
};
use ddg_repulsive_curves::step::schedule::LengthSchedule;
use ddg_repulsive_curves::step::stop::{Status, StepStatus, StopRule};
use ddg_repulsive_curves::step::{
    Backend, BarnesHutEnergy, Control, EnergyModel, FlowOptions, PenaltyTerm, Session, Solver,
    StepReport, flow, soft_options,
};
use ddg_repulsive_curves::subdivide::SubdivisionMap;
use ddg_repulsive_curves::{
    Axis, BhOptions, BhTree, Constraint, CurveNetwork, EnergyEval, Exponents, RcError, Tolerance,
    read_curve_obj,
};

type P = [f64; 3];

// ---------------------------------------------------------------------------------------------
// Curves
// ---------------------------------------------------------------------------------------------

fn loop_edges(n: usize) -> Vec<[u32; 2]> {
    (0..n as u32).map(|i| [i, (i + 1) % n as u32]).collect()
}

/// The trefoil `(sin t + 2 sin 2t, cos t − 2 cos 2t, −sin 3t)`, slightly perturbed.
fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            let w = 0.05 * (5.0 * t).sin();
            [
                t.sin() + 2.0 * (2.0 * t).sin() + w,
                t.cos() - 2.0 * (2.0 * t).cos(),
                -(3.0 * t).sin() + 0.5 * w,
            ]
        })
        .collect();
    (CurveNetwork::new(n, loop_edges(n)).unwrap(), pos)
}

/// A strongly perturbed trefoil (same knot type), far from its energy minimum.
fn rough_trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            let w = 0.35 * (4.0 * t).sin();
            [
                t.sin() + 2.0 * (2.0 * t).sin() + w,
                t.cos() - 2.0 * (2.0 * t).cos() + 0.5 * w,
                -(3.0 * t).sin() + 0.3 * (7.0 * t).cos(),
            ]
        })
        .collect();
    (CurveNetwork::new(n, loop_edges(n)).unwrap(), pos)
}

/// The vertices of a single closed loop in curve order, starting at vertex 0.
fn loop_order(net: &CurveNetwork) -> Vec<usize> {
    let n = net.n_vertices();
    let mut order = Vec::with_capacity(n);
    let (mut prev, mut v) = (usize::MAX, 0usize);
    for _ in 0..n {
        order.push(v);
        let next = net
            .incident_edges(v as u32)
            .iter()
            .map(|&e| {
                let [a, b] = net.edges()[e as usize];
                if a as usize == v {
                    b as usize
                } else {
                    a as usize
                }
            })
            .find(|&w| w != prev)
            .unwrap();
        prev = v;
        v = next;
    }
    order
}

/// The knot determinant of a single loop given by its network.
fn knot_det(net: &CurveNetwork, pos: &[P]) -> i64 {
    let ordered: Vec<P> = loop_order(net).iter().map(|&v| pos[v]).collect();
    knot_determinant(&ordered)
}

/// A wobbly planar flower in the plane z = 0.3.
fn flower(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            let r = 1.0 + 0.3 * (3.0 * t).cos();
            [r * t.cos(), r * t.sin(), 0.3]
        })
        .collect();
    (CurveNetwork::new(n, loop_edges(n)).unwrap(), pos)
}

fn small_options() -> FlowOptions {
    FlowOptions {
        threads: Some(4),
        dense_threads: Some(2),
        ..FlowOptions::default()
    }
}

fn mean_len(net: &CurveNetwork, pos: &[P]) -> f64 {
    net.mean_edge_length(pos)
}

// ---------------------------------------------------------------------------------------------
// Independent audit: brute force over all pairs, plain Lipschitz constant, sampling
// ---------------------------------------------------------------------------------------------

fn norm3(a: P) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Whether any two edges without a common vertex touch along the linear motion `a → b`
/// (all pairs; conservative advancement with the plain constant `max|u| + max|w|` down to
/// `1e-9 d₀`; when the advances run out, 400 samples decide).
fn piece_crosses(net: &CurveNetwork, a: &[P], b: &[P]) -> Option<(usize, usize)> {
    let e = net.edges();
    let at = |v: usize, t: f64| -> P {
        [
            a[v][0] + t * (b[v][0] - a[v][0]),
            a[v][1] + t * (b[v][1] - a[v][1]),
            a[v][2] + t * (b[v][2] - a[v][2]),
        ]
    };
    let speed = |v: usize| norm3([b[v][0] - a[v][0], b[v][1] - a[v][1], b[v][2] - a[v][2]]);
    for i in 0..e.len() {
        for j in i + 1..e.len() {
            if net.edges_adjacent(i, j) {
                continue;
            }
            let [i1, i2] = [e[i][0] as usize, e[i][1] as usize];
            let [j1, j2] = [e[j][0] as usize, e[j][1] as usize];
            let dist = |t: f64| seg_dist(at(i1, t), at(i2, t), at(j1, t), at(j2, t));
            let l = speed(i1).max(speed(i2)) + speed(j1).max(speed(j2));
            let d0 = dist(0.0);
            if d0 <= 0.0 {
                return Some((i, j));
            }
            let mut t = 0.0;
            let mut done = false;
            for _ in 0..10_000 {
                let d = dist(t);
                if d <= 1e-9 * d0 {
                    return Some((i, j));
                }
                if l == 0.0 {
                    done = true;
                    break;
                }
                t += d / l;
                if t >= 1.0 {
                    done = true;
                    break;
                }
            }
            if !done && (0..=400).any(|k| dist(k as f64 / 400.0) <= 0.0) {
                return Some((i, j));
            }
        }
    }
    None
}

fn audit_path(net: &CurveNetwork, path: &[Vec<P>]) -> Option<(usize, usize, usize)> {
    for (k, w) in path.windows(2).enumerate() {
        if let Some((i, j)) = piece_crosses(net, &w[0], &w[1]) {
            return Some((k, i, j));
        }
    }
    None
}

// ---------------------------------------------------------------------------------------------
// Knot determinant |Δ(−1)| from a generic projection (Fox colouring matrix)
// ---------------------------------------------------------------------------------------------

/// `|det|` of any `(n−1)` minor of the colouring matrix of the diagram seen along a fixed
/// generic direction: 1 for the unknot, 3 for the trefoil. Single closed loop `0 → 1 → … → 0`.
fn knot_determinant(pos: &[P]) -> i64 {
    let n = pos.len();
    // A fixed generic rotation.
    let (a, b, c): (f64, f64, f64) = (0.3712, 1.1093, 0.5521);
    let rot = |p: P| -> P {
        let (ca, sa, cb, sb, cc, sc) = (a.cos(), a.sin(), b.cos(), b.sin(), c.cos(), c.sin());
        let x = [p[0], ca * p[1] - sa * p[2], sa * p[1] + ca * p[2]];
        let y = [cb * x[0] + sb * x[2], x[1], -sb * x[0] + cb * x[2]];
        [cc * y[0] - sc * y[1], sc * y[0] + cc * y[1], y[2]]
    };
    let q: Vec<P> = pos.iter().map(|&p| rot(p)).collect();
    // Crossings: (over position, under position) as curve parameters edge + s.
    let mut crossings: Vec<(f64, f64)> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            let (p1, p2) = (q[i], q[(i + 1) % n]);
            let (r1, r2) = (q[j], q[(j + 1) % n]);
            let d1 = [p2[0] - p1[0], p2[1] - p1[1]];
            let d2 = [r2[0] - r1[0], r2[1] - r1[1]];
            let den = d1[0] * d2[1] - d1[1] * d2[0];
            if den == 0.0 {
                continue;
            }
            let w = [r1[0] - p1[0], r1[1] - p1[1]];
            let s = (w[0] * d2[1] - w[1] * d2[0]) / den;
            let t = (w[0] * d1[1] - w[1] * d1[0]) / den;
            if (0.0..1.0).contains(&s) && (0.0..1.0).contains(&t) {
                let zi = p1[2] + s * (p2[2] - p1[2]);
                let zj = r1[2] + t * (r2[2] - r1[2]);
                let (pi, pj) = (i as f64 + s, j as f64 + t);
                crossings.push(if zi > zj { (pi, pj) } else { (pj, pi) });
            }
        }
    }
    let m = crossings.len();
    if m < 3 {
        return 1;
    }
    let mut unders: Vec<f64> = crossings.iter().map(|c| c.1).collect();
    unders.sort_by(f64::total_cmp);
    // Arc k runs from under-crossing k to under-crossing k + 1 (cyclically).
    let arc_of = |p: f64| -> usize {
        match unders.iter().rposition(|&u| u < p) {
            Some(k) => k,
            None => m - 1,
        }
    };
    let mut mat = vec![vec![0.0f64; m]; m];
    for (row, &(over, under)) in crossings.iter().enumerate() {
        let k = unders.iter().position(|&u| u == under).unwrap();
        let incoming = (k + m - 1) % m;
        let outgoing = k;
        mat[row][arc_of(over)] += 2.0;
        mat[row][incoming] -= 1.0;
        mat[row][outgoing] -= 1.0;
    }
    // Delete the last row and column; determinant by Gaussian elimination.
    let k = m - 1;
    let mut a: Vec<Vec<f64>> = mat[..k].iter().map(|r| r[..k].to_vec()).collect();
    let mut det = 1.0;
    for col in 0..k {
        let piv = (col..k)
            .max_by(|&x, &y| a[x][col].abs().total_cmp(&a[y][col].abs()))
            .unwrap();
        if a[piv][col].abs() < 1e-12 {
            return 0;
        }
        if piv != col {
            a.swap(piv, col);
            det = -det;
        }
        det *= a[col][col];
        for r in col + 1..k {
            let f = a[r][col] / a[col][col];
            for cc in col..k {
                a[r][cc] -= f * a[col][cc];
            }
        }
    }
    det.abs().round() as i64
}

#[test]
fn knot_determinant_of_known_curves() {
    let (_, t) = trefoil(60);
    assert_eq!(knot_determinant(&t), 3);
    let (_, f) = flower(40);
    assert_eq!(knot_determinant(&f), 1);
    // A figure-eight knot: determinant 5.
    let n = 120;
    let e8: Vec<P> = (0..n)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            [
                (2.0 + (2.0 * t).cos()) * (3.0 * t).cos(),
                (2.0 + (2.0 * t).cos()) * (3.0 * t).sin(),
                (4.0 * t).sin(),
            ]
        })
        .collect();
    assert_eq!(knot_determinant(&e8), 5);
}

// ---------------------------------------------------------------------------------------------
// Flow tests
// ---------------------------------------------------------------------------------------------

fn check_constraints(s: &Session) {
    let phi = s.constraints().values(s.network(), s.positions()).unwrap();
    let tol = s
        .constraints()
        .row_tolerances(
            s.network(),
            s.positions(),
            s.options().projection.tolerance,
            s.mean_edge_length0(),
        )
        .unwrap();
    for (k, (p, t)) in phi.iter().zip(&tol).enumerate() {
        assert!(p.abs() <= *t, "row {k}: violation {p} > {t}");
    }
}

#[test]
fn guarded_flow_decreases_energy_and_keeps_the_knot() {
    let (net, pos) = rough_trefoil(48);
    assert_eq!(knot_determinant(&pos), 3);
    let opts = FlowOptions {
        record_path: true,
        max_iterations: 25,
        ..small_options()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let e0 = s.energy();
    let mut accepted = 0;
    for _ in 0..25 {
        let r = s.step().unwrap();
        assert!(r.energy_after <= r.energy_before, "{r:?}");
        if r.status == StepStatus::Accepted {
            accepted += 1;
            assert!(r.tau > 0.0 && r.tau <= r.tau_start);
            assert!(r.tau_max.unwrap() > 0.0);
            // Armijo on the projected curve.
            assert!(r.energy_after <= r.energy_before - 1e-4 * r.tau * r.sobolev_norm_sq);
            let path = r.path.as_ref().unwrap();
            assert!(path.len() >= 2);
            assert_eq!(audit_path(s.network(), path), None, "step {}", r.iteration);
            assert_eq!(
                knot_det(s.network(), s.positions()),
                3,
                "step {}",
                r.iteration
            );
        }
        check_constraints(&s);
        if s.status().is_terminal() {
            break;
        }
    }
    assert!(accepted >= 10, "{accepted} accepted steps");
    eprintln!(
        "guarded trefoil: {e0} -> {} in {accepted} steps",
        s.energy()
    );
    assert!(s.energy() < 0.9 * e0, "{} -> {}", e0, s.energy());
}

/// A test environment: a linear objective that decreases along the motion, no constraints,
/// the library's CCD.
struct LinearEnv<'a> {
    net: &'a CurveNetwork,
    weights: Vec<P>,
    pool: rayon::ThreadPool,
    ccd: CcdOptions,
}

impl SearchEnv for LinearEnv<'_> {
    fn objective(&mut self, p: &[P]) -> Result<f64, RcError> {
        Ok(p.iter()
            .zip(&self.weights)
            .map(|(x, w)| x[0] * w[0] + x[1] * w[1] + x[2] * w[2])
            .sum())
    }
    fn project(&mut self, trial: Vec<P>) -> Result<(Projection, Vec<Vec<P>>), RcError> {
        Ok((
            Projection {
                positions: trial.clone(),
                status: ProjectionStatus::Converged,
                iterations: 0,
                max_violation: 0.0,
                history: vec![],
                linear: vec![],
            },
            vec![trial],
        ))
    }
    fn ccd(&mut self, s: &[P], v: &[P], t: f64) -> Result<CcdReport, RcError> {
        ddg_repulsive_curves::ccd::first_contact(&self.pool, self.net, s, v, t, &[], &self.ccd)
    }
    fn cancelled(&self) -> bool {
        false
    }
}

#[test]
fn adversarial_step_tunnels_only_without_the_guard() {
    // A large square in z = 0 and a small square just above one of its edges, pushed down
    // through it: a thin gap (0.05) and a long step.
    let mut pos: Vec<P> = vec![[-2., -2., 0.], [2., -2., 0.], [2., 2., 0.], [-2., 2., 0.]];
    pos.extend([
        [1.5, -0.5, 0.05],
        [2.5, -0.5, 0.05],
        [2.5, 0.5, 0.05],
        [1.5, 0.5, 0.05],
    ]);
    let edges = vec![
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 4],
    ];
    let net = CurveNetwork::new(8, edges).unwrap();
    // Objective Σ z of the small square: the direction p = (0, 0, 1) on it (trials γ − τ p).
    let dir: Vec<P> = (0..8)
        .map(|v| if v >= 4 { [0., 0., 1.] } else { [0.; 3] })
        .collect();
    let weights = dir.clone();
    let slope: f64 = 4.0;
    let mut env = LinearEnv {
        net: &net,
        weights,
        pool: rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap(),
        ccd: CcdOptions::default(),
    };
    let f0 = env.objective(&pos).unwrap();
    let input = SearchInput {
        base: &pos,
        direction: &dir,
        slope,
        f0,
        previous_tau: Some(1.0),
        max_start: None,
    };
    // Normalized: ‖p‖ = 2, τ₀ = min(1.5 τ_prev, 4/2) = 1.5 → the small square moves down by
    // 1.5 and passes through the edge x = 2 of the large one.
    let r = line_search(&LineSearch::normalized_armijo(), input, &mut env).unwrap();
    assert!(r.accepted);
    assert_eq!(r.tau, 1.5);
    assert!(
        audit_path(&net, &r.path).is_some(),
        "the audit must see the crossing"
    );
    // Guarded: τ_max from the CCD, the step stays above the edge.
    let r = line_search(&LineSearch::default(), input, &mut env).unwrap();
    assert!(r.accepted);
    let tm = r.tau_max.unwrap();
    assert!(tm > 0.0 && tm <= 0.05 * 0.9 + 1e-12, "{tm}");
    assert!(r.tau <= tm);
    assert!(audit_path(&net, &r.path).is_none());
    assert_eq!(r.trials[0].outcome, TrialOutcome::Accepted);
}

#[test]
fn growing_length_with_subdivision_keeps_constraints() {
    // Edge lengths grow by a common factor to 2.5×; one subdivision when the mean edge length
    // passes 2 ℓ̄₀.
    let (net, pos) = trefoil(40);
    let l0 = mean_len(&net, &pos);
    let opts = FlowOptions {
        schedule: Some(LengthSchedule {
            factor: 2.5,
            rate: 0.1,
        }),
        max_subdivisions: 1,
        max_iterations: 40,
        ..small_options()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let mut subdivided_at = None;
    for k in 0..40 {
        let r = s.step().unwrap();
        assert_eq!(s.positions().len(), s.network().n_vertices());
        assert_eq!(s.constraints().n_vertices(), s.network().n_vertices());
        check_constraints(&s);
        if r.subdivided {
            subdivided_at = Some(k);
            assert_eq!(s.network().n_vertices(), 80);
            assert_eq!(s.network().n_edges(), 80);
        }
        if s.schedule_progress() == Some(1.5) && subdivided_at.is_some() {
            break;
        }
    }
    let k = subdivided_at.expect("a subdivision");
    assert!(k > 5);
    assert_eq!(s.schedule_progress(), Some(1.5));
    // Every edge target is 2.5 × (its base) / 2 after the split.
    let ml = mean_len(s.network(), s.positions());
    assert!((ml / l0 - 1.25).abs() < 1e-3, "{}", ml / l0);
    assert_eq!(knot_det(s.network(), s.positions()), 3);
    // More steps on the refined curve keep everything valid.
    for _ in 0..3 {
        let r = s.step().unwrap();
        assert_eq!(r.n_vertices, 80);
        check_constraints(&s);
    }
}

#[test]
fn total_length_schedule_and_subdivision_session_state() {
    // Barycenter + total length, shrink-free growth to 2.2× with subdivision twice allowed.
    let (net, pos) = trefoil(36);
    let total0: f64 = net.edge_geometry(&pos).length.iter().sum();
    let opts = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { target: None },
            Constraint::TotalLength { target: None },
        ],
        schedule: Some(LengthSchedule {
            factor: 2.2,
            rate: 0.1,
        }),
        max_subdivisions: 2,
        max_iterations: 60,
        record_path: true,
        ..small_options()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let mut n_sub = 0;
    for _ in 0..30 {
        let before = s.network().clone();
        let r = s.step().unwrap();
        if let Some(p) = &r.path {
            assert_eq!(audit_path(&before, p), None);
        }
        if r.subdivided {
            n_sub += 1;
        }
        check_constraints(&s);
        if s.schedule_progress() == Some(1.2) && r.status == StepStatus::Accepted {
            break;
        }
    }
    assert!(n_sub >= 1);
    let total: f64 = s.network().edge_geometry(s.positions()).length.iter().sum();
    assert!((total / total0 - 2.2).abs() < 1e-3, "{}", total / total0);
    assert_eq!(knot_det(s.network(), s.positions()), 3);
}

#[test]
fn planar_mode_stays_planar() {
    let (net, pos) = flower(40);
    let opts = FlowOptions {
        planar: Some(Axis::Z),
        max_iterations: 15,
        ..small_options()
    };
    let res = flow(net, pos, opts, None).unwrap();
    assert!(res.report.steps.len() >= 5);
    assert!(res.positions.iter().all(|p| p[2] == 0.3));
    let e: Vec<f64> = res.report.steps.iter().map(|r| r.energy_after).collect();
    assert!(e.windows(2).all(|w| w[1] <= w[0]));
    // Not planar input is rejected.
    let (net, mut pos) = flower(40);
    pos[7][2] = 0.30001;
    let opts = FlowOptions {
        planar: Some(Axis::Z),
        ..small_options()
    };
    assert!(matches!(
        Session::new(net, pos, opts),
        Err(RcError::NotPlanar { vertex: 7, .. })
    ));
}

/// A penalty that fails at its `fail_at`-th evaluation with a gradient.
#[derive(Debug)]
struct Failing {
    count: std::sync::atomic::AtomicUsize,
    fail_at: usize,
}

impl PenaltyTerm for Failing {
    fn name(&self) -> &str {
        "failing"
    }
    fn evaluate(&self, _: &CurveNetwork, _: &[P], g: Option<&mut [P]>) -> Result<f64, RcError> {
        if g.is_some() && self.count.fetch_add(1, Ordering::SeqCst) + 1 >= self.fail_at {
            return Err(RcError::InvalidOption {
                name: "failing",
                value: "test".into(),
                requirement: "none",
            });
        }
        Ok(0.0)
    }
    fn subdivided(&self, _: &SubdivisionMap) -> Result<Arc<dyn PenaltyTerm>, RcError> {
        Err(RcError::Unsupported { what: "test" })
    }
}

#[test]
fn every_stop_status_is_reached() {
    let (net, pos) = trefoil(30);
    let run = |opts: FlowOptions| {
        let mut s = Session::new(net.clone(), pos.clone(), opts).unwrap();
        s.run(None).unwrap()
    };
    // Converged: a loose threshold.
    let r = run(FlowOptions {
        stop: StopRule::GradientNorm {
            norm: Default::default(),
            eps: 1e3,
        },
        ..small_options()
    });
    assert_eq!(r.status, Status::Converged);
    assert_eq!(r.steps.last().unwrap().status, StepStatus::Converged);
    // Converged: relative gradient and energy stall.
    let r = run(FlowOptions {
        stop: StopRule::EnergyStall {
            relative: 0.5,
            window: 2,
        },
        ..small_options()
    });
    assert_eq!(r.status, Status::Converged);
    assert_eq!(r.steps.len(), 2);
    // MaxIterations.
    let r = run(FlowOptions {
        max_iterations: 2,
        ..small_options()
    });
    assert_eq!(r.status, Status::MaxIterations);
    assert_eq!(r.steps.len(), 2);
    // TimeLimit.
    let r = run(FlowOptions {
        time_limit: Some(Duration::ZERO),
        ..small_options()
    });
    assert_eq!(r.status, Status::TimeLimit);
    assert!(r.steps.is_empty());
    // Stuck: no projection iterations and a tiny absolute tolerance.
    let r = run(FlowOptions {
        projection: ddg_repulsive_curves::ProjectionOptions {
            tolerance: Tolerance::Absolute(1e-15),
            max_iterations: 0,
        },
        line_search: LineSearch::CollisionGuarded {
            armijo: 1e-4,
            shrink: 0.5,
            start: StartRule::TwoThirdsTauMax,
            gap: 0.1,
            max_halvings: 2,
        },
        ..small_options()
    });
    assert_eq!(r.status, Status::Stuck);
    let last = r.steps.last().unwrap();
    assert_eq!(last.status, StepStatus::Stuck);
    assert_eq!(last.trials.len(), 3);
    assert!(
        last.trials
            .iter()
            .all(|t| matches!(t.outcome, TrialOutcome::ProjectionFailed(_)))
    );
    // Cancelled by the flag, and by the callback.
    let mut s = Session::new(net.clone(), pos.clone(), small_options()).unwrap();
    s.cancel_flag().store(true, Ordering::Relaxed);
    assert_eq!(s.run(None).unwrap().status, Status::Cancelled);
    assert!(matches!(s.step(), Err(RcError::Cancelled)));
    let mut s = Session::new(net.clone(), pos.clone(), small_options()).unwrap();
    let mut cb = |r: &StepReport| {
        if r.iteration >= 2 {
            Control::Stop
        } else {
            Control::Continue
        }
    };
    let r = s.run(Some(&mut cb)).unwrap();
    assert_eq!(r.status, Status::Cancelled);
    assert_eq!(r.steps.len(), 2);
    // Error: a penalty that fails at the second step.
    let mut s = Session::new(
        net.clone(),
        pos.clone(),
        FlowOptions {
            penalties: vec![Arc::new(Failing {
                count: Default::default(),
                fail_at: 2,
            })],
            ..small_options()
        },
    )
    .unwrap();
    let before = s.positions().to_vec();
    let r = s.run(None).unwrap();
    assert!(matches!(
        r.status,
        Status::Error(RcError::InvalidOption { .. })
    ));
    assert_eq!(r.steps.len(), 1);
    assert_ne!(s.positions(), &before[..]);
    // RetargetStalled: no chord iterations allowed, so neither the re-target nor the descent
    // can move (as in the Stuck case), while the schedule still wants to move the targets.
    let mut s = Session::new(
        net.clone(),
        pos.clone(),
        FlowOptions {
            schedule: Some(LengthSchedule {
                factor: 2.0,
                rate: 0.1,
            }),
            projection: ddg_repulsive_curves::ProjectionOptions {
                tolerance: Tolerance::Absolute(1e-15),
                max_iterations: 0,
            },
            line_search: LineSearch::CollisionGuarded {
                armijo: 1e-4,
                shrink: 0.5,
                start: StartRule::TwoThirdsTauMax,
                gap: 0.1,
                max_halvings: 2,
            },
            ..small_options()
        },
    )
    .unwrap();
    let r = s.step().unwrap();
    let rt = r.retarget.as_ref().unwrap();
    assert_eq!(rt.attempts, 6);
    assert_eq!(
        rt.status,
        ddg_repulsive_curves::step::RetargetStatus::Stalled
    );
    assert_eq!(r.status, StepStatus::RetargetStalled);
    let r = s.run(None).unwrap();
    assert_eq!(r.status, Status::Stuck);
}

#[test]
fn soft_methods_decrease_energy() {
    let (net, pos) = rough_trefoil(40);
    for method in [Method::Lbfgs { memory: 30 }, Method::Ncg] {
        for ls in [LineSearch::default(), LineSearch::normalized_armijo()] {
            let mut opts = soft_options(method, 10.0, &net, &pos).unwrap();
            opts.threads = Some(3);
            opts.dense_threads = Some(1);
            opts.line_search = ls;
            opts.max_iterations = 20;
            opts.record_path = ls.is_guarded();
            let mut s = Session::new(net.clone(), pos.clone(), opts).unwrap();
            let e0 = s.energy();
            let mut prev = e0;
            let mut acc = 0;
            for _ in 0..20 {
                let r = s.step().unwrap();
                assert!(r.energy_after <= prev, "{method:?} {r:?}");
                assert!(r.soft.is_some());
                prev = r.energy_after;
                if r.status == StepStatus::Accepted {
                    acc += 1;
                    if let Some(p) = &r.path {
                        assert_eq!(audit_path(s.network(), p), None);
                    }
                }
                if s.status().is_terminal() {
                    break;
                }
            }
            assert!(acc >= 10, "{method:?} {ls:?}: {acc}");
            assert!(prev < 0.9 * e0, "{method:?} {ls:?}: {e0} -> {prev}");
            let h = s.penalties()[0]
                .evaluate(s.network(), s.positions(), None)
                .unwrap();
            assert!(h.is_finite());
            eprintln!(
                "{method:?} guarded={}: {e0:.4} -> {prev:.4}",
                ls.is_guarded()
            );
        }
    }
    // The Hencky penalty follows a subdivision.
    let h = Hencky::at_rest(2.0, &net, &pos).unwrap();
    let map = SubdivisionMap {
        old_vertices: 40,
        old_edges: 40,
    };
    let h2 = h.subdivided(&map).unwrap();
    let (net2, pos2, _) = ddg_repulsive_curves::subdivide::subdivide(&net, &pos).unwrap();
    assert!(h2.evaluate(&net2, &pos2, None).unwrap().abs() < 1e-24);
}

#[test]
fn classic_metric_flows_decrease_energy() {
    let (net, pos) = trefoil(40);
    for metric in [
        ddg_repulsive_curves::Metric::L2,
        ddg_repulsive_curves::Metric::H1,
        ddg_repulsive_curves::Metric::H2,
    ] {
        let opts = FlowOptions {
            metric,
            max_iterations: 3,
            ..small_options()
        };
        let res = flow(net.clone(), pos.clone(), opts, None).unwrap();
        let s = &res.report.steps;
        assert!(
            s.iter().any(|r| r.status == StepStatus::Accepted),
            "{metric:?}"
        );
        assert!(s.iter().all(|r| r.energy_after <= r.energy_before));
        assert_eq!(s[0].solver_backend, "direct");
    }
}

#[test]
fn results_are_bitwise_independent_of_thread_counts() {
    let (net, pos) = trefoil(36);
    let run = |threads: usize, dense: usize| {
        let opts = FlowOptions {
            threads: Some(threads),
            dense_threads: Some(dense),
            max_iterations: 4,
            ..FlowOptions::default()
        };
        flow(net.clone(), pos.clone(), opts, None).unwrap()
    };
    let a = run(1, 1);
    for (t, d) in [(7, 3), (64, 16)] {
        let b = run(t, d);
        assert_eq!(a.positions, b.positions, "threads {t}, dense {d}");
        assert_eq!(a.report.final_energy, b.report.final_energy);
        let ta: Vec<f64> = a.report.steps.iter().map(|r| r.tau).collect();
        let tb: Vec<f64> = b.report.steps.iter().map(|r| r.tau).collect();
        assert_eq!(ta, tb);
    }
}

#[test]
fn bad_options_are_errors() {
    let (net, pos) = trefoil(30);
    let new = |o: FlowOptions| Session::new(net.clone(), pos.clone(), o);
    let bad_ls = FlowOptions {
        line_search: LineSearch::CollisionGuarded {
            armijo: 1.5,
            shrink: 0.5,
            start: StartRule::TwoThirdsTauMax,
            gap: 0.1,
            max_halvings: 4,
        },
        ..small_options()
    };
    assert!(matches!(new(bad_ls), Err(RcError::InvalidOption { .. })));
    let bad_gap = FlowOptions {
        line_search: LineSearch::CollisionGuarded {
            armijo: 1e-4,
            shrink: 0.5,
            start: StartRule::TwoThirdsTauMax,
            gap: 1.0,
            max_halvings: 4,
        },
        ..small_options()
    };
    assert!(new(bad_gap).is_err());
    let soft_with_constraints = FlowOptions {
        method: Method::Ncg,
        ..small_options()
    };
    assert!(matches!(
        new(soft_with_constraints),
        Err(RcError::MethodNeedsSoftProblem)
    ));
    let no_len = FlowOptions {
        constraints: vec![Constraint::Barycenter { target: None }],
        schedule: Some(LengthSchedule::default()),
        ..small_options()
    };
    assert!(matches!(
        new(no_len),
        Err(RcError::ScheduleWithoutLengthConstraint)
    ));
    let zero_threads = FlowOptions {
        threads: Some(0),
        ..small_options()
    };
    assert!(new(zero_threads).is_err());
    let zero_dense = FlowOptions {
        dense_threads: Some(0),
        ..small_options()
    };
    assert!(new(zero_dense).is_err());
    let bh = FlowOptions {
        backend: ddg_repulsive_curves::step::Backend::Explicit {
            energy: EnergyEval::BarnesHut(BhOptions {
                theta_x: 0.0,
                ..BhOptions::default()
            }),
            solver: ddg_repulsive_curves::step::Solver::Dense,
        },
        ..small_options()
    };
    assert!(matches!(new(bh), Err(RcError::InvalidOption { .. })));
    let no_translation = FlowOptions {
        constraints: vec![Constraint::TotalLength { target: None }],
        ..small_options()
    };
    assert!(matches!(
        new(no_translation),
        Err(RcError::TranslationNotFixed)
    ));
    let bad_schedule = FlowOptions {
        schedule: Some(LengthSchedule {
            factor: -1.0,
            rate: 0.01,
        }),
        ..small_options()
    };
    assert!(new(bad_schedule).is_err());
    let bad_stop = FlowOptions {
        stop: StopRule::EnergyStall {
            relative: 0.1,
            window: 0,
        },
        ..small_options()
    };
    assert!(new(bad_stop).is_err());
    let infeasible = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { target: None },
            Constraint::TotalLength { target: Some(1.0) },
        ],
        ..small_options()
    };
    assert!(matches!(
        new(infeasible),
        Err(RcError::InfeasibleConstraint { .. })
    ));
    let mut wrong = pos.clone();
    wrong.pop();
    assert!(Session::new(net.clone(), wrong, small_options()).is_err());
    let lbfgs0 = FlowOptions {
        method: Method::Lbfgs { memory: 0 },
        constraints: vec![],
        ..small_options()
    };
    assert!(new(lbfgs0).is_err());
    assert!(
        Exponents::new(3.0, 6.0).is_ok() && Exponents::new(1.0, 6.0).is_err(),
        "exponents are validated by construction"
    );
}

// ---------------------------------------------------------------------------------------------
// Freedman–He–Wang 600 (release; ignored)
// ---------------------------------------------------------------------------------------------

/// Roundness of a closed loop (vertices in curve order): `(max r − min r)/mean r` of the
/// distances to the length-weighted barycenter,
/// the planarity `(smallest singular value of the centred positions)/mean r`, and the total
/// turning angle divided by 2π (1 for a planar convex curve, Fenchel).
fn roundness(pos: &[P]) -> (f64, f64, f64) {
    let n = pos.len() as f64;
    // The length-weighted barycenter (the centre of a circle whatever its vertex spacing).
    let mut c = [0.0; 3];
    let mut total = 0.0;
    for i in 0..pos.len() {
        let (a, b) = (pos[i], pos[(i + 1) % pos.len()]);
        let l = norm3([b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
        total += l;
        for k in 0..3 {
            c[k] += l * 0.5 * (a[k] + b[k]);
        }
    }
    for k in 0..3 {
        c[k] /= total;
    }
    let r: Vec<f64> = pos
        .iter()
        .map(|p| norm3([p[0] - c[0], p[1] - c[1], p[2] - c[2]]))
        .collect();
    let mean = r.iter().sum::<f64>() / n;
    let spread =
        (r.iter().cloned().fold(0.0, f64::max) - r.iter().cloned().fold(f64::MAX, f64::min)) / mean;
    // Covariance and its smallest eigenvalue by power iteration on (tr·I − C).
    let mut cov = [[0.0; 3]; 3];
    for p in pos {
        let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j] / n;
            }
        }
    }
    let tr = cov[0][0] + cov[1][1] + cov[2][2];
    let mut v = [0.3, 0.5, 0.8];
    let mut lam = 0.0;
    for _ in 0..500 {
        let w: Vec<f64> = (0..3)
            .map(|i| tr * v[i] - (0..3).map(|j| cov[i][j] * v[j]).sum::<f64>())
            .collect();
        let nw = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        lam = nw;
        v = [w[0] / nw, w[1] / nw, w[2] / nw];
    }
    let smallest = (tr - lam).max(0.0).sqrt();
    let m = pos.len();
    let mut turning = 0.0;
    for i in 0..m {
        let a = pos[(i + m - 1) % m];
        let b = pos[i];
        let cc = pos[(i + 1) % m];
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [cc[0] - b[0], cc[1] - b[1], cc[2] - b[2]];
        let cosang = (u[0] * w[0] + u[1] * w[1] + u[2] * w[2]) / (norm3(u) * norm3(w));
        turning += cosang.clamp(-1.0, 1.0).acos();
    }
    (spread, smallest / mean, turning / std::f64::consts::TAU)
}

fn data_dir() -> String {
    format!("{}/../../data/repulsive-curves", env!("CARGO_MANIFEST_DIR"))
}

/// `cargo test --release -p ddg-repulsive-curves --test flow freedman -- --ignored --nocapture`.
#[test]
#[ignore]
fn freedman_600_flows_to_a_circle() {
    freedman_flow(EnergyEval::Exact);
}

/// The Freedman flow on the Barnes–Hut energy (dense solver): `RC_ALPHA`, `RC_BETA`, `RC_STEPS`
/// as above. Prints the exact energy of the final curve and compares it with the exact flow's
/// (`6.0816` at `|g| < 1e-4`, `runs_tmp/rc-wp3-notes.md`).
#[test]
#[ignore]
fn freedman_600_barnes_hut_flows_to_a_circle() {
    let s = freedman_flow(EnergyEval::BarnesHut(BhOptions::default()));
    let exact = ddg_repulsive_curves::energy(
        s.network(),
        s.positions(),
        Exponents::new(2.0, 4.5).unwrap(),
        &EnergyEval::Exact,
        Some(8),
    )
    .unwrap()
    .energy;
    eprintln!(
        "RESULT bh final: Barnes-Hut energy {:.5}, exact energy of the same curve {exact:.5}",
        s.energy()
    );
    // Since the differential is the exact derivative of the evaluated energy (SPEC 17i) the flow
    // converges (step 711, |g| 9.86e-5, exact energy 6.08164) instead of getting stuck at the
    // noise floor of the old differential.
    assert_eq!(s.status(), Status::Converged);
    // The same near-circle: the exact energy is within 0.5 % of the exact flow's 6.0816.
    assert!((exact - 6.0816).abs() < 0.005 * 6.0816, "exact {exact}");
    assert!((s.energy() - exact).abs() < 0.02 * exact);
}

fn freedman_flow(eval: EnergyEval) -> Session {
    let path = format!("{}/FreedmanHeWang/freedman-he-wang-600.obj", data_dir());
    let (net, pos) = read_curve_obj(&path).unwrap();
    let alpha: f64 = std::env::var("RC_ALPHA")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2.0);
    let beta: f64 = std::env::var("RC_BETA")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4.5);
    let budget: usize = std::env::var("RC_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1500);
    let opts = FlowOptions {
        exponents: Exponents::new(alpha, beta).unwrap(),
        max_iterations: budget,
        record_path: false,
        backend: Backend::Explicit {
            energy: eval,
            solver: Solver::Dense,
        },
        ..FlowOptions::default()
    };
    let t0 = Instant::now();
    let mut s = Session::new(net, pos, opts).unwrap();
    let e0 = s.energy();
    let (sp0, pl0, tu0) = roundness(s.positions());
    eprintln!(
        "start: |V| {} energy {e0:.4} spread {sp0:.3} planarity {pl0:.3} turning/2pi {tu0:.3}",
        s.network().n_vertices()
    );
    let mut round_at = None;
    let mut split = [0.0f64; 6];
    let mut min_dist = f64::INFINITY;
    loop {
        let r = s.step().unwrap();
        split[0] += r.seconds.energy;
        split[1] += r.seconds.metric;
        split[2] += r.seconds.solve;
        split[3] += r.seconds.line_search_energy;
        split[4] += r.seconds.projection;
        split[5] += r.seconds.ccd;
        let (sp, pl, tu) = roundness(s.positions());
        if r.iteration.is_multiple_of(25) || s.status().is_terminal() {
            min_dist = ddg_repulsive_curves::ccd::min_pair_distance(s.network(), s.positions());
            eprintln!(
                "step {:4} {:?} E {:.5} |g| {:.2e} tau {:.2e} tau_max {:?} halv {} proj {} \
                 spread {sp:.4} planarity {pl:.4} turning/2pi {tu:.4} min dist {min_dist:.3e} t {:.1}s",
                r.iteration,
                r.status,
                r.energy_after,
                r.gradient_norm,
                r.tau,
                r.tau_max,
                r.halvings,
                r.projection_iterations,
                t0.elapsed().as_secs_f64()
            );
        }
        if round_at.is_none() && sp < 0.02 && pl < 0.02 && tu < 1.01 {
            round_at = Some((r.iteration, t0.elapsed().as_secs_f64()));
            eprintln!(
                "near-circle at step {} ({:.1} s)",
                r.iteration,
                t0.elapsed().as_secs_f64()
            );
        }
        if s.status().is_terminal() {
            break;
        }
    }
    eprintln!(
        "end: status {:?} steps {} energy {:.5} time {:.1}s; split energy {:.1} metric {:.1} \
         solve {:.1} ls-energy {:.1} projection {:.1} ccd {:.1}",
        s.status(),
        s.iteration(),
        s.energy(),
        t0.elapsed().as_secs_f64(),
        split[0],
        split[1],
        split[2],
        split[3],
        split[4],
        split[5]
    );
    check_constraints(&s);
    assert!(min_dist > 0.0);
    let (it, secs) = round_at.expect("the unknot becomes a near-circle within the budget");
    eprintln!("RESULT near-circle iterations {it} seconds {secs:.1}");
    if matches!(eval, EnergyEval::BarnesHut(_)) {
        bh_differential_consistency(&s, eval);
    }
    s
}

/// Diagnostic for the Barnes-Hut stall: along random directions `v`, compares the central
/// difference of the evaluated Barnes-Hut energy (tree refitted, topology fixed, as in the line
/// search), the Barnes-Hut differential `<dE_bh, v>` and the exact differential `<dE, v>`.
fn bh_differential_consistency(s: &Session, eval: EnergyEval) {
    let EnergyEval::BarnesHut(opts) = eval else {
        return;
    };
    let e = Exponents::new(2.0, 4.5).unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();
    let (net, pos) = (s.network(), s.positions());
    let tree = BhTree::build(net, pos, opts).unwrap();
    let (_, d_bh, _) = tree.energy_and_differential(&pool, e);
    let lists = tree.interaction_lists(&pool);
    let exact = ddg_repulsive_curves::energy(net, pos, e, &EnergyEval::Exact, Some(8)).unwrap();
    let mut seed = 12345u64;
    let mut rnd = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    };
    for trial in 0..4 {
        let v: Vec<P> = pos.iter().map(|_| [rnd(), rnd(), rnd()]).collect();
        let h = 1e-5;
        let at_frozen = |sign: f64| {
            let q: Vec<P> = pos
                .iter()
                .zip(&v)
                .map(|(p, d)| {
                    [
                        p[0] + sign * h * d[0],
                        p[1] + sign * h * d[1],
                        p[2] + sign * h * d[2],
                    ]
                })
                .collect();
            let mut t = tree.clone();
            t.refit(&q).unwrap();
            t.energy_with_lists(&lists, &pool, e)
        };
        let fd_frozen = (at_frozen(1.0) - at_frozen(-1.0)) / (2.0 * h);
        let at = |sign: f64| {
            let q: Vec<P> = pos
                .iter()
                .zip(&v)
                .map(|(p, d)| {
                    [
                        p[0] + sign * h * d[0],
                        p[1] + sign * h * d[1],
                        p[2] + sign * h * d[2],
                    ]
                })
                .collect();
            let mut t = tree.clone();
            t.refit(&q).unwrap();
            t.energy(&pool, e).0
        };
        let fd = (at(1.0) - at(-1.0)) / (2.0 * h);
        let dot = |d: &[P]| -> f64 {
            d.iter()
                .zip(&v)
                .map(|(a, b)| a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
                .sum()
        };
        let (bh, ex) = (dot(&d_bh), dot(&exact.differential));
        eprintln!(
            "DIFF-CHECK {trial}: FD(E_bh) {fd:.6e}  <dE_bh,v> {bh:.6e}  <dE_exact,v> {ex:.6e}  \
             FD(frozen lists) {fd_frozen:.6e}  |FD-dE_bh| {:.2e}  |dE_bh-dE_exact| {:.2e}",
            (fd - bh).abs(),
            (bh - ex).abs()
        );
    }
}

// ---------------------------------------------------------------------------------------------
// Line-search comparison driver (one knot per process; ignored)
// ---------------------------------------------------------------------------------------------

/// `RC_KNOT=path RC_LS=guarded|grow|armijo cargo test --release --test flow ls_compare -- --ignored
/// --nocapture`: flows one knot with barycenter + edge lengths and audits every accepted step
/// with the CCD (tiny gap, many advances). Prints one `RESULT` line.
#[test]
#[ignore]
fn ls_compare_one_knot() {
    let Ok(path) = std::env::var("RC_KNOT") else {
        eprintln!("RC_KNOT not set; nothing to do");
        return;
    };
    let ls = match std::env::var("RC_LS").as_deref() {
        Ok("armijo") => LineSearch::normalized_armijo(),
        Ok("grow") => LineSearch::guarded_with_start(StartRule::grow_previous()),
        _ => LineSearch::guarded_with_start(StartRule::TwoThirdsTauMax),
    };
    let ls_name = std::env::var("RC_LS").unwrap_or_else(|_| "guarded".into());
    let alpha: f64 = std::env::var("RC_ALPHA")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2.0);
    let beta: f64 = std::env::var("RC_BETA")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4.5);
    let budget: usize = std::env::var("RC_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(300);
    let threads: usize = std::env::var("RC_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);
    let (net, pos) = read_curve_obj(&path).unwrap();
    let opts = FlowOptions {
        exponents: Exponents::new(alpha, beta).unwrap(),
        line_search: ls,
        projection: ddg_repulsive_curves::ProjectionOptions {
            tolerance: Tolerance::Absolute(1e-4),
            max_iterations: 3,
        },
        max_iterations: budget,
        record_path: true,
        threads: Some(threads),
        dense_threads: Some(threads),
        ..FlowOptions::default()
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let audit_opts = CcdOptions {
        gap: 1e-6,
        max_advances: 100_000,
    };
    let t0 = Instant::now();
    let mut s = Session::new(net, pos, opts).unwrap();
    let e0 = s.energy();
    let mut violations = 0;
    let mut first_violation = None;
    let (mut audit_s, mut trials, mut stuck_steps) = (0.0, 0, 0);
    loop {
        let r = s.step().unwrap();
        trials += r.trials.len();
        if r.status == StepStatus::Stuck {
            stuck_steps += 1;
        }
        println!(
            "TRAJ {} {:.4} {:.10} {:.3e} {}",
            r.iteration,
            t0.elapsed().as_secs_f64() - audit_s,
            r.energy_after,
            r.gradient_norm,
            r.trials.len()
        );
        if let Some(p) = &r.path {
            let ta = Instant::now();
            let mut rep = CcdReport::default();
            let hit = ddg_repulsive_curves::ccd::path_contact(
                &pool,
                s.network(),
                p,
                &[],
                &audit_opts,
                &mut rep,
            )
            .unwrap();
            audit_s += ta.elapsed().as_secs_f64();
            if let Some((piece, c)) = hit {
                violations += 1;
                first_violation.get_or_insert((r.iteration, piece, c.primitive));
            }
        }
        if s.status().is_terminal() {
            break;
        }
    }
    let secs = t0.elapsed().as_secs_f64() - audit_s;
    let last_g = s.iteration();
    println!(
        "RESULT knot {} ls {} status {:?} iterations {} seconds {:.2} energy0 {:.6} energy {:.6} \
         trials {} stuck_steps {} violations {} first_violation {:?} audit_seconds {:.1}",
        std::path::Path::new(&path)
            .file_name()
            .unwrap()
            .to_string_lossy(),
        ls_name,
        s.status(),
        last_g,
        secs,
        e0,
        s.energy(),
        trials,
        stuck_steps,
        violations,
        first_violation,
        audit_s
    );
}

#[test]
#[ignore]
fn debug_trefoil_steps() {
    let (net, pos) = trefoil(48);
    let ls = match std::env::var("RC_LS").as_deref() {
        Ok("armijo") => LineSearch::normalized_armijo(),
        _ => LineSearch::default(),
    };
    let opts = FlowOptions {
        line_search: ls,
        max_iterations: 40,
        ..small_options()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    for _ in 0..40 {
        let r = s.step().unwrap();
        eprintln!(
            "{} {:?} E {:.5} -> {:.5} |g| {:.3e} slope {:.3e} tau0 {:.3e} tau {:.3e} tau_hi {:?} tau_max {:?} halv {} trials {:?}",
            r.iteration,
            r.status,
            r.energy_before,
            r.energy_after,
            r.gradient_norm,
            r.sobolev_norm_sq,
            r.tau_start,
            r.tau,
            r.tau_hi,
            r.tau_max,
            r.halvings,
            r.trials
                .iter()
                .map(|t| format!("{:?}", t.outcome))
                .collect::<Vec<_>>()
        );
    }
}

/// Factorization time of the dense saddle of Freedman–He–Wang 600/1200 against the size of
/// the dense pool (`cargo test --release --test flow dense_pool -- --ignored --nocapture`).
#[test]
#[ignore]
fn dense_pool_scaling() {
    use ddg_repulsive_curves::constraints::ConstraintSet;
    use ddg_repulsive_curves::metric::{DofLayout, GramQuadrature, sobolev_gram};
    use ddg_repulsive_curves::solve::DenseSaddle;
    for size in [600, 1200] {
        let path = format!("{}/FreedmanHeWang/freedman-he-wang-{size}.obj", data_dir());
        let (net, pos) = read_curve_obj(&path).unwrap();
        let cons = ddg_repulsive_curves::default_constraints();
        let cs = ConstraintSet::new(&net, &pos, &cons, None).unwrap();
        let c = cs.jacobian(&net, &pos).unwrap();
        let main = rayon::ThreadPoolBuilder::new()
            .num_threads(64)
            .build()
            .unwrap();
        let g = sobolev_gram(
            &main,
            &net,
            &pos,
            Exponents::default(),
            GramQuadrature::Midpoint,
        );
        for threads in [1, 4, 8, 16, 24, 32, 48, 64] {
            let pool = Arc::new(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap(),
            );
            let mut best = f64::INFINITY;
            for _ in 0..3 {
                let s = DenseSaddle::new(pool.clone(), &g.a, &c, DofLayout::new(None)).unwrap();
                best = best.min(s.factor_seconds());
            }
            eprintln!(
                "RESULT |V| {} threads {threads:2} factor {best:.3} s",
                net.n_vertices()
            );
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Barnes–Hut in the session (D26, SPEC §17d)
// ---------------------------------------------------------------------------------------------

fn bh_flow_options(threads: usize, dense: usize) -> FlowOptions {
    FlowOptions {
        backend: Backend::Explicit {
            energy: EnergyEval::BarnesHut(BhOptions::default()),
            solver: Solver::Dense,
        },
        max_iterations: 8,
        threads: Some(threads),
        dense_threads: Some(dense),
        ..FlowOptions::default()
    }
}

#[test]
fn barnes_hut_refit_matches_rebuild() {
    let (net, pos) = rough_trefoil(80);
    let e = Exponents::new(3.0, 6.0).unwrap();
    let pool = Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap(),
    );
    let opts = BhOptions::default();
    // Unchanged positions: the refit is the tree that was built.
    let mut model = BarnesHutEnergy::new(pool.clone(), e, opts);
    model.prepare(&net, &pos).unwrap();
    let built = BhTree::build(&net, &pos, opts).unwrap().energy(&pool, e).0;
    assert_eq!(model.energy(&net, &pos).unwrap().to_bits(), built.to_bits());
    // A small move, as the line search makes: the refit tree is a valid (slightly different)
    // tree, so the energy agrees with a rebuild and with the exact energy to the approximation.
    let moved: Vec<P> = pos
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let t = 0.01 * (i as f64 * 0.7).sin();
            [p[0] + t, p[1] - 0.5 * t, p[2] + 0.3 * t]
        })
        .collect();
    let refit = model.energy(&net, &moved).unwrap();
    let rebuilt = BhTree::build(&net, &moved, opts)
        .unwrap()
        .energy(&pool, e)
        .0;
    let exact = ddg_repulsive_curves::energy(&net, &moved, e, &EnergyEval::Exact, Some(2))
        .unwrap()
        .energy;
    assert!(
        (refit - rebuilt).abs() < 1e-3 * rebuilt,
        "{refit} {rebuilt}"
    );
    assert!((refit - exact).abs() < 1e-2 * exact, "{refit} {exact}");
    // A new `prepare` (next step, or after a subdivision) rebuilds and so equals a rebuild.
    let (net2, pos2, _) = ddg_repulsive_curves::subdivide(&net, &moved).unwrap();
    model.prepare(&net2, &pos2).unwrap();
    let ev = model.energy_and_differential(&net2, &pos2).unwrap();
    let want = BhTree::build(&net2, &pos2, opts)
        .unwrap()
        .energy(&pool, e)
        .0;
    assert_eq!(ev.energy.to_bits(), want.to_bits());
    assert_eq!(ev.differential.len(), net2.n_vertices());
}

#[test]
fn barnes_hut_flow_follows_the_exact_flow() {
    let (net, pos) = rough_trefoil(60);
    let mut bh = Session::new(net.clone(), pos.clone(), bh_flow_options(2, 2)).unwrap();
    let mut ex = Session::new(
        net.clone(),
        pos,
        FlowOptions {
            backend: Backend::Explicit {
                energy: EnergyEval::Exact,
                solver: Solver::Dense,
            },
            ..bh_flow_options(2, 2)
        },
    )
    .unwrap();
    let e0 = bh.energy();
    let (r_bh, r_ex) = (bh.run(None).unwrap(), ex.run(None).unwrap());
    assert!(r_bh.steps.iter().all(|s| s.status != StepStatus::Stuck));
    assert!(bh.energy() < 0.9 * e0, "the flow decreases the energy");
    // Same energy within the approximation (the tree is accurate to about 1 %).
    assert!(
        (bh.energy() - ex.energy()).abs() < 0.03 * ex.energy(),
        "{} {}",
        bh.energy(),
        ex.energy()
    );
    // Every accepted step decreased the Barnes-Hut energy (both sides of Armijo, §17d).
    let mut prev = e0;
    for s in &r_bh.steps {
        assert!(s.energy_after <= prev, "{} {prev}", s.energy_after);
        prev = s.energy_after;
    }
    drop(r_ex);
    check_constraints(&bh);
}

#[test]
fn barnes_hut_flow_is_deterministic_across_thread_counts() {
    let (net, pos) = rough_trefoil(60);
    let run = |threads: usize| {
        let mut s = Session::new(net.clone(), pos.clone(), bh_flow_options(threads, 2)).unwrap();
        s.run(None).unwrap();
        (
            s.energy().to_bits(),
            s.positions()
                .iter()
                .flat_map(|p| p.map(f64::to_bits))
                .collect::<Vec<_>>(),
        )
    };
    let (e1, p1) = run(1);
    for threads in [3, 8] {
        let (e, p) = run(threads);
        assert_eq!(e, e1, "energy differs with {threads} threads");
        assert_eq!(p, p1, "positions differ with {threads} threads");
    }
}

#[test]
fn grow_previous_start_keeps_the_guard_and_needs_fewer_trials() {
    let (net, pos) = rough_trefoil(60);
    let run = |start: StartRule| {
        let mut s = Session::new(
            net.clone(),
            pos.clone(),
            FlowOptions {
                line_search: LineSearch::guarded_with_start(start),
                max_iterations: 12,
                record_path: true,
                ..small_options()
            },
        )
        .unwrap();
        let mut trials = 0;
        for _ in 0..12 {
            let r = s.step().unwrap();
            trials += r.trials.len();
            // Every trial is bounded by tau_max, whatever the start rule.
            if let Some(tm) = r.tau_max {
                assert!(r.trials.iter().all(|t| t.tau <= tm * (1.0 + 1e-12)));
            }
        }
        (trials, s.energy())
    };
    let (paper, e_paper) = run(StartRule::TwoThirdsTauMax);
    let (grow, e_grow) = run(StartRule::grow_previous());
    assert!(grow <= paper, "grow {grow} paper {paper}");
    assert!(e_grow.is_finite() && e_paper.is_finite());
    assert!(
        LineSearch::guarded_with_start(StartRule::GrowPrevious { factor: 0.5 })
            .validate()
            .is_err()
    );
}

#[test]
fn barnes_hut_flow_survives_subdivision() {
    // Regression: after a live subdivision the Barnes-Hut tree has the old vertex count
    // (PositionCountMismatch at the first energy evaluation of the new network).
    let (net, pos) = trefoil(40);
    let opts = FlowOptions {
        schedule: Some(LengthSchedule {
            factor: 2.5,
            rate: 0.1,
        }),
        max_subdivisions: 1,
        max_iterations: 40,
        ..bh_flow_options(2, 2)
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let mut subdivided = false;
    for _ in 0..40 {
        let r = s.step().unwrap();
        assert_eq!(s.positions().len(), s.network().n_vertices());
        subdivided |= r.subdivided;
    }
    assert!(subdivided, "the schedule subdivides");
    assert_eq!(s.network().n_vertices(), 80);
    assert!(s.energy().is_finite());
}
