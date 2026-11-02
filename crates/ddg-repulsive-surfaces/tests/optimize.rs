//! The comparison optimizers (WP9): decrease, two-loop = BFGS, AQP momentum, soft problems,
//! reference variants, determinism.

use ddg_repulsive_curves::optimize::Lbfgs;
use ddg_repulsive_surfaces::constraints::Constraint;
use ddg_repulsive_surfaces::generate::Rng;
use ddg_repulsive_surfaces::generate::primitives as generate;
use ddg_repulsive_surfaces::optimize::{aqp_theta, soft_problem};
use ddg_repulsive_surfaces::{
    BcqnVariant, FlowOptions, LbfgsVariant, Method, Metric, Penalty, RemeshOptions, RsError, Session,
    StopRule, Surface,
};

fn perturbed(mesh: ddg_repulsive_surfaces::generate::Mesh, a: f64, seed: u64) -> (Surface, Vec<[f64; 3]>) {
    let mut rng = Rng::new(seed);
    let p: Vec<[f64; 3]> = mesh
        .positions
        .iter()
        .map(|q| {
            let s = 1.0 + a * rng.uniform(-1.0, 1.0);
            [q[0] * s, q[1] * s * 0.9, q[2] * s * 1.2]
        })
        .collect();
    (Surface::new(p.len(), mesh.faces.clone()).unwrap(), p)
}

fn sphere() -> (Surface, Vec<[f64; 3]>) {
    perturbed(generate::sphere(2, 1.0, [0.0; 3]).unwrap(), 0.08, 3)
}

fn torus() -> (Surface, Vec<[f64; 3]>) {
    perturbed(generate::torus(16, 8, 1.0, 0.4).unwrap(), 0.05, 4)
}

fn hard() -> FlowOptions {
    FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        stop: StopRule::Off,
        threads: Some(2),
        ..FlowOptions::default()
    }
}

fn soft(s: &Surface, p: &[[f64; 3]], method: Method, threads: usize) -> FlowOptions {
    let mut o = soft_problem(s, p, &hard()).unwrap();
    o.method = method;
    o.threads = Some(threads);
    o
}

fn methods() -> Vec<Method> {
    vec![
        Method::Aqp { kappa: 100.0 },
        Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Standard },
        Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Reference },
        Method::Bcqn { memory: 20, variant: BcqnVariant::Standard },
        Method::Bcqn { memory: 20, variant: BcqnVariant::Reference },
    ]
}

fn run(s: &Surface, p: &[[f64; 3]], o: FlowOptions, steps: usize) -> Session {
    let mut ses = Session::new(s.clone(), p.to_vec(), o).unwrap();
    for _ in 0..steps {
        ses.step().unwrap();
    }
    ses
}

#[test]
fn soft_problem_replaces_area() {
    let (s, p) = sphere();
    let o = soft_problem(&s, &p, &hard()).unwrap();
    assert_eq!(o.constraints.len(), 1);
    assert!(matches!(o.penalties.as_slice(), [Penalty::AreaDeviation { weight, .. }] if *weight == 100.0));
}

#[test]
fn needs_soft_problem_names_constraint() {
    let (s, p) = sphere();
    let o = FlowOptions { method: Method::Aqp { kappa: 100.0 }, ..hard() };
    let e = Session::new(s.clone(), p.clone(), o).unwrap_err();
    assert!(matches!(&e, RsError::MethodNeedsSoftProblem { constraint, .. } if constraint.contains("Area")));
    assert!(e.to_string().contains("Area"));
    let o = FlowOptions {
        method: Method::Bcqn { memory: 5, variant: BcqnVariant::Standard },
        constraints: vec![Constraint::Volume { component: None, target: None }],
        ..hard()
    };
    assert!(matches!(Session::new(s, p, o).unwrap_err(), RsError::MethodNeedsSoftProblem { .. }));
}

#[test]
fn every_method_decreases_energy() {
    for (name, (s, p)) in [("sphere", sphere()), ("torus", torus())] {
        for m in methods() {
            let o = soft(&s, &p, m, 2);
            let mut ses = Session::new(s.clone(), p.clone(), o).unwrap();
            let e0 = ses.energy();
            for _ in 0..12 {
                ses.step().unwrap();
            }
            let e1 = ses.energy();
            assert!(e1 < e0, "{name} {m:?}: {e0} -> {e1}");
        }
    }
}

#[test]
fn first_lbfgs_step_is_the_h1_step() {
    let (s, p) = sphere();
    let m = Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Standard };
    let a = run(&s, &p, soft(&s, &p, m, 2), 1);
    let mut o = soft(&s, &p, Method::ProjectedGradient, 2);
    o.metric = Metric::H1;
    let b = run(&s, &p, o, 1);
    assert_eq!(a.positions(), b.positions());
}

#[test]
fn aqp_momentum_and_resets() {
    assert!((aqp_theta(100.0) - 9.0 / 11.0).abs() < 1e-15);
    let (s, p) = sphere();
    let o = soft(&s, &p, Method::Aqp { kappa: 100.0 }, 2);
    let ses = run(&s, &p, o.clone(), 4);
    let d = ses.method_diagnostics().unwrap();
    assert_eq!(d.theta, Some(aqp_theta(100.0)));
    let mut h1 = soft(&s, &p, Method::ProjectedGradient, 2);
    h1.metric = Metric::H1;
    let plain = run(&s, &p, h1, 4);
    assert_ne!(ses.positions(), plain.positions());
    // Remeshing resets the momentum.
    let (s2, p2) = perturbed(generate::sphere(2, 1.0, [0.0; 3]).unwrap(), 0.3, 9);
    let mut ro = soft(&s2, &p2, Method::Aqp { kappa: 100.0 }, 2);
    ro.remesh = Some(RemeshOptions::default());
    let ses = run(&s2, &p2, ro, 5);
    let d = ses.method_diagnostics().unwrap();
    assert!(d.remeshes > 0, "{d:?}");
}

#[test]
fn remeshing_resets_lbfgs_memory() {
    let (s, p) = perturbed(generate::sphere(2, 1.0, [0.0; 3]).unwrap(), 0.3, 9);
    let m = Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Standard };
    let mut o = soft(&s, &p, m, 2);
    o.remesh = Some(RemeshOptions::default());
    let ses = run(&s, &p, o, 5);
    let d = ses.method_diagnostics().unwrap();
    assert!(d.remeshes > 0, "{d:?}");
}

/// The two-loop recursion equals the explicit BFGS inverse update of `H₀ = γ M⁻¹`.
#[test]
fn two_loop_equals_explicit_bfgs() {
    let n = 6;
    let mut rng = Rng::new(11);
    let b: Vec<f64> = (0..n * n).map(|_| rng.uniform(-1.0, 1.0)).collect();
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] =
                (0..n).map(|k| b[i * n + k] * b[j * n + k]).sum::<f64>() + if i == j { 1.0 } else { 0.0 };
        }
    }
    let mulv = |v: &[f64]| -> Vec<f64> { (0..n).map(|i| (0..n).map(|j| a[i * n + j] * v[j]).sum()).collect() };
    let dot = |x: &[f64], y: &[f64]| -> f64 { x.iter().zip(y).map(|(p, q)| p * q).sum() };
    let mut lb = Lbfgs::new(10);
    let mut pairs = Vec::new();
    for _ in 0..4 {
        let s: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let y = mulv(&s);
        assert!(lb.push(s.clone(), y.clone()));
        pairs.push((s, y));
    }
    let g: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
    let m_inv = |q: &[f64]| -> Vec<f64> { q.iter().map(|v| v / 2.0).collect() };
    let z = lb.direction(&g, &mut |q| Ok(m_inv(q))).unwrap();
    let (sl, yl) = pairs.last().unwrap();
    let myl = m_inv(yl);
    let gamma = dot(sl, yl) / dot(yl, &myl);
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        h[i * n + i] = gamma / 2.0;
    }
    for (s, y) in &pairs {
        let rho = 1.0 / dot(s, y);
        let mut v = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                v[i * n + j] = if i == j { 1.0 } else { 0.0 } - rho * y[i] * s[j];
            }
        }
        let mut hv = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                hv[i * n + j] = (0..n).map(|k| h[i * n + k] * v[k * n + j]).sum();
            }
        }
        for i in 0..n {
            for j in 0..n {
                h[i * n + j] = (0..n).map(|k| v[k * n + i] * hv[k * n + j]).sum::<f64>() + rho * s[i] * s[j];
            }
        }
    }
    let hg: Vec<f64> = (0..n).map(|i| (0..n).map(|j| h[i * n + j] * g[j]).sum()).collect();
    for i in 0..n {
        assert!((hg[i] - z[i]).abs() < 1e-12 * (1.0 + z[i].abs()), "{i}: {} vs {}", hg[i], z[i]);
    }
}

#[test]
fn reference_variants_differ() {
    let (s, p) = sphere();
    for (a, b) in [
        (
            Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Standard },
            Method::H1Lbfgs { memory: 20, variant: LbfgsVariant::Reference },
        ),
        (
            Method::Bcqn { memory: 20, variant: BcqnVariant::Standard },
            Method::Bcqn { memory: 20, variant: BcqnVariant::Reference },
        ),
    ] {
        let x = run(&s, &p, soft(&s, &p, a, 2), 6);
        let y = run(&s, &p, soft(&s, &p, b, 2), 6);
        assert_ne!(x.positions(), y.positions(), "{a:?}");
    }
    let x = run(&s, &p, soft(&s, &p, Method::Bcqn { memory: 20, variant: BcqnVariant::Standard }, 2), 4);
    let y = run(&s, &p, soft(&s, &p, Method::Bcqn { memory: 20, variant: BcqnVariant::Reference }, 2), 4);
    let (dx, dy) = (x.method_diagnostics().unwrap(), y.method_diagnostics().unwrap());
    assert!(dx.last_norm.unwrap() < dy.last_norm.unwrap(), "{dx:?} {dy:?}");
}

#[test]
fn deterministic_across_threads() {
    let (s, p) = torus();
    for m in methods() {
        let a = run(&s, &p, soft(&s, &p, m, 1), 5);
        let b = run(&s, &p, soft(&s, &p, m, 4), 5);
        assert_eq!(a.positions(), b.positions(), "{m:?}");
    }
}
