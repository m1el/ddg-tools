//! The back ends of the session (WP3c, `SPEC.md` §8.5, §17h): the thread pools reach the
//! solvers, `Solver::Multigrid` with live subdivision and constraints, the dense fallback,
//! `Backend::Auto`.

use std::sync::Arc;

use ddg_repulsive_curves::{
    Backend, BhOptions, CurveNetwork, EnergyEval, FlowOptions, MgOptions, Session, Solver,
};

type P = [f64; 3];

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
    let edges = (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

fn small_mg() -> MgOptions {
    MgOptions {
        coarsest_vertices: 40,
        ..MgOptions::default()
    }
}

fn options(energy: EnergyEval, solver: Solver, threads: usize) -> FlowOptions {
    FlowOptions {
        backend: Backend::Explicit { energy, solver },
        threads: Some(threads),
        ..FlowOptions::default()
    }
}

/// Three steps with every energy / solver pairing give bitwise identical curves at 1, 3 and 8
/// threads: the solvers run on the session's pools only.
#[test]
fn every_backend_is_bitwise_independent_of_the_thread_count() {
    let (net, pos) = trefoil(200);
    let energies = [
        ("exact", EnergyEval::Exact),
        ("bh", EnergyEval::BarnesHut(BhOptions::default())),
    ];
    for (en, energy) in energies {
        for sn in ["dense", "multigrid"] {
            let mut curves = Vec::new();
            for t in [1, 3, 8] {
                let solver = match sn {
                    "dense" => Solver::Dense,
                    _ => Solver::Multigrid(small_mg()),
                };
                let mut s =
                    Session::new(net.clone(), pos.clone(), options(energy, solver, t)).unwrap();
                for _ in 0..3 {
                    s.step().unwrap();
                }
                curves.push(s.positions().to_vec());
            }
            assert!(
                curves[1] == curves[0] && curves[2] == curves[0],
                "{en} + {sn}: curves differ between 1, 3 and 8 threads"
            );
        }
    }
}

// ------------------------------------------------------------------ Solver::Multigrid

use ddg_repulsive_curves::constraints::{Constraint, TangentSite};
use ddg_repulsive_curves::{DirectSolvers, SolverRequest};
use ddg_repulsive_curves::{LengthSchedule, RcError, SaddleSolution, SaddleSolver, SolverFactory};

fn mg_options(energy: EnergyEval, threads: usize) -> FlowOptions {
    options(energy, Solver::Multigrid(small_mg()), threads)
}

fn max_violation(s: &Session) -> f64 {
    let phi = s.constraints().values(s.network(), s.positions()).unwrap();
    phi.iter().fold(0.0, |m, p| m.max(p.abs()))
}

/// `Solver::Multigrid` follows the dense solver closely with both energy models, reports the
/// multigrid solve and never falls back on a healthy problem.
#[test]
fn multigrid_session_follows_the_dense_session() {
    let (net, pos) = trefoil(160);
    for energy in [
        EnergyEval::Exact,
        EnergyEval::BarnesHut(BhOptions::default()),
    ] {
        let mut mg = Session::new(net.clone(), pos.clone(), mg_options(energy, 4)).unwrap();
        let mut dense =
            Session::new(net.clone(), pos.clone(), options(energy, Solver::Dense, 4)).unwrap();
        for _ in 0..6 {
            let rm = mg.step().unwrap();
            let rd = dense.step().unwrap();
            assert_eq!(rm.solver_backend, "multigrid");
            assert_eq!(rm.linear.solver, "multigrid");
            assert!(rm.linear.levels >= 3, "{:?}", rm.linear);
            assert!(rm.solver_fallback.is_none());
            assert!(rd.solver_fallback.is_none());
            let rel = (rm.energy_after - rd.energy_after).abs() / rd.energy_after.abs();
            assert!(
                rel < 1e-3,
                "energies {} {}",
                rm.energy_after,
                rd.energy_after
            );
        }
        assert!(max_violation(&mg) < 1e-3);
    }
}

/// Live subdivision: the hierarchy is rebuilt for the refined network (more levels), the
/// constraints follow, and no step needs the fallback. Both energy models.
#[test]
fn multigrid_session_survives_live_subdivision() {
    let (net, pos) = trefoil(60);
    for energy in [
        EnergyEval::Exact,
        EnergyEval::BarnesHut(BhOptions::default()),
    ] {
        let opts = FlowOptions {
            schedule: Some(LengthSchedule {
                factor: 2.5,
                rate: 0.1,
            }),
            max_subdivisions: 1,
            max_iterations: 40,
            ..mg_options(energy, 4)
        };
        let mut s = Session::new(net.clone(), pos.clone(), opts).unwrap();
        let (mut before, mut after) = (0, 0);
        for _ in 0..40 {
            let r = s.step().unwrap();
            assert!(r.solver_fallback.is_none(), "{:?}", r.solver_fallback);
            assert_eq!(s.positions().len(), s.network().n_vertices());
            assert_eq!(r.solver_backend, "multigrid");
            if r.n_vertices == 60 {
                before = before.max(r.linear.levels);
            } else {
                after = after.max(r.linear.levels);
            }
        }
        assert_eq!(s.network().n_vertices(), 120);
        assert!(before >= 2 && after > before, "levels {before} -> {after}");
        assert!(max_violation(&s) < 1e-3);
    }
}

/// Constraints on an open curve: both ends pinned, one end tangent fixed, edge lengths fixed.
/// The coarse levels keep the pins and the tangent site (SPEC §8.4) and the flow agrees with the
/// dense one.
#[test]
fn multigrid_session_handles_pins_and_tangents() {
    let n = 120usize;
    let pos: Vec<P> = (0..n)
        .map(|k| {
            let t = k as f64 / (n - 1) as f64;
            [
                3.0 * t,
                0.6 * (6.0 * t).sin(),
                0.4 * (9.0 * t).cos() * (1.0 - t),
            ]
        })
        .collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    let net = CurveNetwork::new(n, edges).unwrap();
    let constraints = vec![
        Constraint::EdgeLengths {
            edges: None,
            targets: None,
        },
        Constraint::Point {
            vertex: 0,
            target: None,
        },
        Constraint::Point {
            vertex: n as u32 - 1,
            target: None,
        },
        Constraint::Tangent {
            site: TangentSite::Vertex(0),
            direction: None,
        },
    ];
    let mk = |solver| FlowOptions {
        constraints: constraints.clone(),
        ..options(EnergyEval::Exact, solver, 4)
    };
    let mut mg = Session::new(net.clone(), pos.clone(), mk(Solver::Multigrid(small_mg()))).unwrap();
    let mut dense = Session::new(net, pos.clone(), mk(Solver::Dense)).unwrap();
    // The block cluster tree's product is an approximation of the dense metric (θ = 0.5): the
    // projected gradients of the first step agree to 2 %, later steps are different flows.
    let (mut previous, e0) = (mg.energy(), mg.energy());
    for k in 0..8 {
        let rm = mg.step().unwrap();
        let rd = dense.step().unwrap();
        assert!(rm.solver_fallback.is_none(), "{:?}", rm.solver_fallback);
        assert!(rm.linear.levels >= 2, "{:?}", rm.linear);
        if k == 0 {
            let rel = (rm.gradient_norm - rd.gradient_norm).abs() / rd.gradient_norm;
            assert!(
                rel < 0.02,
                "gradient norms {} {}",
                rm.gradient_norm,
                rd.gradient_norm
            );
        }
        assert!(rm.energy_after <= previous);
        previous = rm.energy_after;
    }
    assert!(previous < 0.97 * e0, "{previous} {e0}");
    assert!(max_violation(&mg) < 1e-3);
    assert!(
        mg.positions()[0] == pos[0] || {
            let d = mg.positions()[0]
                .iter()
                .zip(&pos[0])
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            d < 1e-3
        }
    );
}

// ------------------------------------------------------------------ fallback

/// Multigrid options that cannot converge in one cycle at a 1e-14 tolerance.
fn failing_mg() -> MgOptions {
    MgOptions {
        tolerance: 1e-14,
        failure_residual: 1e-13,
        max_cycles: 1,
        ..small_mg()
    }
}

/// A multigrid solve that stagnates above `failure_residual` ends in the dense solver of the
/// same step: the step's result equals the dense session's bit for bit, and the report says so.
#[test]
fn stagnating_multigrid_falls_back_to_the_dense_solver() {
    let (net, pos) = trefoil(100);
    let mut mg = Session::new(
        net.clone(),
        pos.clone(),
        options(EnergyEval::Exact, Solver::Multigrid(failing_mg()), 4),
    )
    .unwrap();
    let mut dense = Session::new(net, pos, options(EnergyEval::Exact, Solver::Dense, 4)).unwrap();
    for _ in 0..3 {
        let rm = mg.step().unwrap();
        let rd = dense.step().unwrap();
        let f = rm.solver_fallback.as_ref().expect("the step fell back");
        assert_eq!((f.from, f.to), ("multigrid", "direct"));
        assert!(
            matches!(f.cause, RcError::MultigridNotConverged { .. }),
            "{:?}",
            f.cause
        );
        assert!(f.dense_solves >= 1);
        assert!(f.cause.to_string().contains("did not converge"));
        assert_eq!(rm.linear.solver, rd.linear.solver);
        assert!(rd.solver_fallback.is_none());
    }
    assert!(mg.positions() == dense.positions());
}

/// Above `dense_fallback_max_vertices` the error goes up and the status is `Error`.
#[test]
fn fallback_is_refused_above_the_size_limit() {
    let (net, pos) = trefoil(100);
    let opts = FlowOptions {
        dense_fallback_max_vertices: 99,
        ..options(EnergyEval::Exact, Solver::Multigrid(failing_mg()), 4)
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let e = s.step().unwrap_err();
    assert!(matches!(e, RcError::MultigridNotConverged { .. }), "{e:?}");
    assert!(matches!(s.status(), ddg_repulsive_curves::Status::Error(_)));
}

/// A saddle solver that fails with `MultigridDiverged` on the first solve.
struct Diverging(Box<dyn SaddleSolver>, bool);

impl SaddleSolver for Diverging {
    fn n_dofs(&self) -> usize {
        self.0.n_dofs()
    }
    fn n_constraints(&self) -> usize {
        self.0.n_constraints()
    }
    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError> {
        if self.1 {
            return Err(RcError::MultigridDiverged {
                residuals: vec![1e-2, 1.0, 40.0],
            });
        }
        self.0.solve(a, b)
    }
}

struct DivergingFactory;

impl SolverFactory for DivergingFactory {
    fn name(&self) -> &'static str {
        "diverging"
    }
    fn build(&self, r: &SolverRequest<'_>) -> Result<Box<dyn SaddleSolver>, RcError> {
        Ok(Box::new(Diverging(DirectSolvers::new().build(r)?, true)))
    }
}

/// `MultigridDiverged` from any iterative factory also falls back, and the report names the
/// factory.
#[test]
fn diverging_solver_falls_back_and_names_the_cause() {
    let (net, pos) = trefoil(60);
    let mut s = Session::new(
        net,
        pos,
        options(
            EnergyEval::Exact,
            Solver::Custom(Arc::new(DivergingFactory)),
            2,
        ),
    )
    .unwrap();
    let r = s.step().unwrap();
    let f = r.solver_fallback.expect("fell back");
    assert_eq!(f.from, "diverging");
    assert!(matches!(f.cause, RcError::MultigridDiverged { .. }));
    assert!(f.cause.to_string().contains("diverged"));
    assert!(f.dense_solves >= 1);
}

// ------------------------------------------------------------------ crossover (SPEC §17h)

/// Time per step of one back end on Freedman–He–Wang resampled to `RC_SIZE` vertices
/// (300: every other vertex of the 600 file). Environment: `RC_SIZE`, `RC_ENERGY`
/// (`exact`|`bh`), `RC_SOLVER` (`dense`|`mg`), `RC_STEPS` (5), `RC_THREADS` (8). Prints one
/// line. `cargo test --release --test session_backends crossover -- --ignored --nocapture`.
#[test]
#[ignore = "measurement; reads data/repulsive-curves"]
fn crossover_one_configuration() {
    let var = |k: &str, d: &str| std::env::var(k).unwrap_or_else(|_| d.to_string());
    let size: usize = var("RC_SIZE", "600").parse().unwrap();
    let steps: usize = var("RC_STEPS", "5").parse().unwrap();
    let threads: usize = var("RC_THREADS", "8").parse().unwrap();
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/repulsive-curves");
    let file = match size {
        300 => 600,
        600 | 1200 | 2400 => size,
        _ => 1200,
    };
    let Ok((net, pos)) = ddg_repulsive_curves::read_curve_obj(
        dir.join(format!("FreedmanHeWang/freedman-he-wang-{file}.obj")),
    ) else {
        eprintln!("data/repulsive-curves not present; skipping");
        return;
    };
    let (net, pos) = if file == size {
        (net, pos)
    } else if size == 300 {
        let p: Vec<P> = pos.iter().step_by(2).copied().collect();
        let n = p.len();
        let e = (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect();
        (CurveNetwork::new(n, e).unwrap(), p)
    } else {
        // Between the files: split evenly spread edges at their midpoints of the 1200 curve
        // (the shape is unchanged, the edge lengths differ by a factor 2).
        let n = pos.len();
        assert!(size > n, "sizes between 1200 and 2400 only");
        let extra = size - n;
        let mut split = vec![false; n];
        for j in 0..extra {
            split[j * n / extra] = true;
        }
        let mut p: Vec<P> = Vec::with_capacity(size);
        for i in 0..n {
            let (a, b) = (pos[i], pos[(i + 1) % n]);
            p.push(a);
            if split[i] {
                p.push([0, 1, 2].map(|c| 0.5 * (a[c] + b[c])));
            }
        }
        let e = (0..size as u32)
            .map(|k| [k, (k + 1) % size as u32])
            .collect();
        (CurveNetwork::new(size, e).unwrap(), p)
    };
    {
        let g = net.edge_geometry(&pos);
        let (mn, mx) = g
            .length
            .iter()
            .fold((f64::MAX, 0.0f64), |(a, b), &l| (a.min(l), b.max(l)));
        eprintln!(
            "resampled: {} vertices, edge length {mn:.3e} .. {mx:.3e}",
            pos.len()
        );
    }
    let energy_name = var("RC_ENERGY", "exact");
    let solver_name = var("RC_SOLVER", "dense");
    let energy = match energy_name.as_str() {
        "bh" => EnergyEval::BarnesHut(BhOptions::default()),
        _ => EnergyEval::Exact,
    };
    let solver = match solver_name.as_str() {
        "mg" => Solver::Multigrid(MgOptions::default()),
        _ => Solver::Dense,
    };
    let opts = FlowOptions {
        exponents: ddg_repulsive_curves::Exponents::new(2.0, 4.5).unwrap(),
        dense_threads: Some(threads),
        ..options(energy, solver, threads)
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let (mut tot, mut en, mut me, mut so, mut ls) = (vec![], 0.0, 0.0, 0.0, 0.0);
    let mut fallbacks = 0;
    for _ in 0..steps {
        let r = s.step().unwrap();
        tot.push(r.seconds.total);
        en += r.seconds.energy;
        me += r.seconds.metric;
        so += r.seconds.solve;
        ls += r.seconds.line_search;
        fallbacks += usize::from(r.solver_fallback.is_some());
    }
    let n = steps as f64;
    let mut sorted = tot.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eprintln!(
        "CROSSOVER size {size} energy {energy_name} solver {solver_name} threads {threads} steps {steps} \
         mean_step {:.3} median_step {:.3} energy {:.3} metric {:.3} solve {:.3} line_search {:.3} fallbacks {fallbacks} E {:.6}",
        tot.iter().sum::<f64>() / n,
        sorted[steps / 2],
        en / n,
        me / n,
        so / n,
        ls / n,
        s.energy()
    );
}

// ------------------------------------------------------------------ Backend::Auto

/// The rule of `Backend::Auto`: exact + dense up to the crossover, Barnes–Hut + multigrid above.
#[test]
fn auto_resolves_by_vertex_count() {
    let auto = Backend::Auto {
        dense_max_vertices: 100,
    };
    for (n, accelerated) in [(1, false), (100, false), (101, true), (100_000, true)] {
        let (e, s) = auto.resolve(n);
        match (accelerated, e, s) {
            (false, EnergyEval::Exact, Solver::Dense) => {}
            (true, EnergyEval::BarnesHut(o), Solver::Multigrid(m)) => {
                assert_eq!(o, BhOptions::default());
                assert_eq!(m, MgOptions::default());
            }
            (_, e, s) => panic!("{n} vertices: {e:?} + {s:?}"),
        }
    }
    assert_eq!(
        Backend::default().resolve(2500).1.to_string_for_test(),
        "Dense"
    );
    assert_eq!(
        Backend::default().resolve(2501).1.to_string_for_test(),
        "Multigrid"
    );
}

trait SolverName {
    fn to_string_for_test(&self) -> &'static str;
}

impl SolverName for Solver {
    fn to_string_for_test(&self) -> &'static str {
        match self {
            Solver::Dense => "Dense",
            Solver::Multigrid(_) => "Multigrid",
            Solver::Custom(_) => "Custom",
        }
    }
}

/// A session on `Auto` uses exact + dense on a small curve and switches to Barnes–Hut +
/// multigrid when a live subdivision takes the curve above the crossover; the flow goes on.
#[test]
fn auto_session_switches_back_ends_after_a_subdivision() {
    let (net, pos) = trefoil(60);
    let opts = FlowOptions {
        backend: Backend::Auto {
            dense_max_vertices: 100,
        },
        threads: Some(4),
        schedule: Some(LengthSchedule {
            factor: 2.5,
            rate: 0.1,
        }),
        max_subdivisions: 1,
        max_iterations: 40,
        ..FlowOptions::default()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let mut seen = std::collections::BTreeSet::new();
    for _ in 0..40 {
        let r = s.step().unwrap();
        assert!(r.solver_fallback.is_none());
        seen.insert((r.n_vertices, r.energy_backend, r.solver_backend));
    }
    assert_eq!(s.network().n_vertices(), 120);
    let want: std::collections::BTreeSet<_> =
        [(60, "exact", "direct"), (120, "barnes-hut", "multigrid")]
            .into_iter()
            .collect();
    assert_eq!(seen, want);
}

// ------------------------------------------------------------------ run vs step on errors

struct FailingAfter(std::sync::atomic::AtomicUsize, usize);

impl SolverFactory for FailingAfter {
    fn name(&self) -> &'static str {
        "failing-after"
    }
    fn build(&self, r: &SolverRequest<'_>) -> Result<Box<dyn SaddleSolver>, RcError> {
        let k = self.0.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        if k >= self.1 {
            return Err(RcError::SingularMetric { components: 1 });
        }
        DirectSolvers::new().build(r)
    }
}

/// A failing step: `step` returns the error and sets `Status::Error`; `run` returns `Ok` with
/// `Status::Error` carrying the same error, the steps completed before it, and the curve of the
/// last completed step.
#[test]
fn run_ends_with_status_error_while_step_returns_the_error() {
    let (net, pos) = trefoil(40);
    let mk = || {
        Session::new(
            net.clone(),
            pos.clone(),
            options(
                EnergyEval::Exact,
                Solver::Custom(Arc::new(FailingAfter(Default::default(), 3))),
                2,
            ),
        )
        .unwrap()
    };
    // run
    let mut s = mk();
    let report = s.run(None).unwrap();
    assert_eq!(report.steps.len(), 3);
    assert_eq!(
        report.status,
        ddg_repulsive_curves::Status::Error(RcError::SingularMetric { components: 1 })
    );
    assert_eq!(
        report.error(),
        Some(&RcError::SingularMetric { components: 1 })
    );
    assert_eq!(s.iteration(), 3);
    let after_three = s.positions().to_vec();
    // step, same failure
    let mut t = mk();
    for _ in 0..3 {
        t.step().unwrap();
    }
    assert!(t.positions() == after_three.as_slice());
    let e = t.step().unwrap_err();
    assert_eq!(e, RcError::SingularMetric { components: 1 });
    assert_eq!(t.status(), ddg_repulsive_curves::Status::Error(e));
    assert!(t.positions() == after_three.as_slice());
    // Cancellation is its own status, in both.
    let mut c = Session::new(
        net.clone(),
        pos.clone(),
        options(EnergyEval::Exact, Solver::Dense, 2),
    )
    .unwrap();
    c.cancel_flag()
        .store(true, std::sync::atomic::Ordering::Relaxed);
    assert_eq!(c.step().unwrap_err(), RcError::Cancelled);
    assert_eq!(c.status(), ddg_repulsive_curves::Status::Cancelled);
    assert_eq!(
        c.run(None).unwrap().status,
        ddg_repulsive_curves::Status::Cancelled
    );
}

/// The report of a length schedule: the factor asked for, the factor imposed and the measured one.
#[test]
fn step_report_exposes_the_schedule_progress() {
    let (net, pos) = trefoil(40);
    let opts = FlowOptions {
        schedule: Some(LengthSchedule {
            factor: 1.5,
            rate: 0.1,
        }),
        max_iterations: 20,
        ..options(EnergyEval::Exact, Solver::Dense, 2)
    };
    let mut s = Session::new(net.clone(), pos.clone(), opts).unwrap();
    let mut last = None;
    for k in 0..8 {
        let r = s.step().unwrap();
        let p = r.schedule.expect("a schedule is set");
        assert_eq!(p.factor_target, 1.5);
        assert!((p.span - 0.5).abs() < 1e-15);
        assert!((p.factor_imposed - (1.0 + p.progress)).abs() < 1e-12);
        assert!((p.fraction - p.progress / 0.5).abs() < 1e-12);
        assert!(
            (p.factor_measured - p.factor_imposed).abs() < 1e-3,
            "step {k}: measured {} imposed {}",
            p.factor_measured,
            p.factor_imposed
        );
        assert_eq!(p.finished, p.progress >= 0.5);
        last = Some(p);
    }
    let p = last.unwrap();
    assert!(p.finished && p.fraction == 1.0 && (p.factor_measured - 1.5).abs() < 1e-3);
    // No schedule: no field.
    let mut t = Session::new(net, pos, options(EnergyEval::Exact, Solver::Dense, 2)).unwrap();
    assert!(t.step().unwrap().schedule.is_none());
}
