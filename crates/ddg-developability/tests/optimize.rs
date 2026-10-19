//! Tests of the generic optimizer (SPEC §4, §10.5) on classic smooth and nonsmooth functions.

use ddg_developability::optimize::{
    ClosureObjective, Control, LineSearch, Method, Minimized, OptimizeError, OptimizeOptions,
    OptimizeStatus, Problem, minimize,
};

fn run<F: FnMut(&[f64], &mut [f64]) -> Option<f64>>(
    f: F,
    x0: &[f64],
    group: usize,
    opts: &OptimizeOptions,
) -> Minimized {
    let mut obj = ClosureObjective(f);
    minimize(&mut obj, x0, &Problem { length_scale: 1.0, group }, opts, None).unwrap()
}

fn rosenbrock(x: &[f64], g: &mut [f64]) -> Option<f64> {
    let (a, b) = (x[0], x[1]);
    g[0] = -2.0 * (1.0 - a) - 400.0 * a * (b - a * a);
    g[1] = 200.0 * (b - a * a);
    Some((1.0 - a).powi(2) + 100.0 * (b - a * a).powi(2))
}

fn quadratic(x: &[f64], g: &mut [f64]) -> Option<f64> {
    // diag(1, 10, ..., 10^3) spread, minimum at (1, 2, ..., n)
    let n = x.len();
    let mut f = 0.0;
    for i in 0..n {
        let h = 10f64.powf(3.0 * i as f64 / (n - 1) as f64);
        let r = x[i] - (i as f64 + 1.0);
        f += 0.5 * h * r * r;
        g[i] = h * r;
    }
    Some(f)
}

fn report_line(name: &str, m: &Minimized) {
    eprintln!(
        "{name}: status {:?}, iterations {}, evaluations {}, f {:.3e}, resets {}, skipped {}, weak {}, capped {}",
        m.report.status,
        m.report.iterations,
        m.report.evaluations,
        m.report.energy_final,
        m.report.memory_resets,
        m.report.skipped_pairs,
        m.report.weak_steps,
        m.report.capped_steps
    );
}

fn monotone(m: &Minimized) {
    for w in m.report.energy_history.windows(2) {
        assert!(w[1] <= w[0], "energy increased: {} -> {}", w[0], w[1]);
    }
}

#[test]
fn lbfgs_quadratic_converges() {
    let x0 = vec![0.0; 10];
    let m = run(quadratic, &x0, 1, &OptimizeOptions { max_iterations: 500, ..Default::default() });
    report_line("lbfgs quadratic n=10", &m);
    monotone(&m);
    assert!(m.report.energy_final < 1e-14, "{:?}", m.report.status);
    for (i, v) in m.x.iter().enumerate() {
        assert!((v - (i as f64 + 1.0)).abs() < 1e-6);
    }
    // superlinear: the last decades of gradient norm take few iterations
    let gh = &m.report.gradient_norm_history;
    assert!(gh.len() > 5 && *gh.last().unwrap() < 1e-6 * gh[0]);
}

#[test]
fn lbfgs_rosenbrock_converges() {
    let m = run(rosenbrock, &[-1.2, 1.0], 1, &OptimizeOptions::default());
    report_line("lbfgs rosenbrock", &m);
    monotone(&m);
    assert!((m.x[0] - 1.0).abs() < 1e-6 && (m.x[1] - 1.0).abs() < 1e-6);
}

#[test]
fn gradient_descent_both_line_searches_decrease() {
    for (name, ls) in [
        ("gd weak wolfe", OptimizeOptions::default().line_search),
        (
            "gd backtracking",
            LineSearch::Backtracking { c1: 1e-4, shrink: 0.5, max_shrinks: 50 },
        ),
    ] {
        let opts = OptimizeOptions {
            method: Method::GradientDescent,
            line_search: ls,
            max_iterations: 300,
            ..Default::default()
        };
        let m = run(quadratic, &[0.0; 4], 1, &opts);
        report_line(name, &m);
        monotone(&m);
        assert!(m.report.energy_final < 0.01 * m.report.energy_initial);
    }
}

#[test]
fn fixed_step_counts_increases() {
    // 1D f = x^2 /2, step 0.5 contracts; step 2.5 diverges (energy increases every time).
    let f = |x: &[f64], g: &mut [f64]| {
        g[0] = x[0];
        Some(0.5 * x[0] * x[0])
    };
    let good = OptimizeOptions { method: Method::FixedStep { step: 0.5 }, max_iterations: 20, ..Default::default() };
    let m = run(f, &[1.0], 1, &good);
    assert_eq!(m.report.energy_increases, 0);
    assert!((m.x[0] - 0.5f64.powi(20)).abs() < 1e-15);
    let bad = OptimizeOptions { method: Method::FixedStep { step: 2.5 }, max_iterations: 5, ..Default::default() };
    let m = run(f, &[1.0], 1, &bad);
    assert_eq!(m.report.energy_increases, 5);
}

#[test]
fn abs_value_bfgs_alternates() {
    // L&O Sec. 4: BFGS on |x| with the weak Wolfe search; iterates alternate in sign, |x| decreases.
    let f = |x: &[f64], g: &mut [f64]| {
        g[0] = if x[0] >= 0.0 { 1.0 } else { -1.0 };
        Some(x[0].abs())
    };
    let opts = OptimizeOptions {
        method: Method::Lbfgs { memory: 20 },
        max_iterations: 40,
        gradient_tolerance: 0.0,
        energy_tolerance: 0.0,
        ..Default::default()
    };
    let mut xs = vec![];
    let mut obj = ClosureObjective(f);
    let mut cb = |p: &ddg_developability::optimize::FlatProgress| {
        xs.push(p.x[0]);
        Control::Continue
    };
    let m = minimize(&mut obj, &[3.3], &Problem { length_scale: 10.0, group: 1 }, &opts, Some(&mut cb)).unwrap();
    report_line("bfgs |x|", &m);
    monotone(&m);
    assert!(m.report.energy_final < 1e-3, "f = {}", m.report.energy_final);
    // After the first iterations the sign alternates.
    let alt = xs.windows(2).filter(|w| w[0] * w[1] < 0.0).count();
    assert!(alt >= xs.len() / 3, "sign changes {alt} of {}", xs.len());
}

#[test]
fn euclidean_norm_nonsmooth() {
    let f = |x: &[f64], g: &mut [f64]| {
        let n = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        for (gi, xi) in g.iter_mut().zip(x) {
            *gi = if n > 0.0 { xi / n } else { 0.0 };
        }
        Some(n)
    };
    let opts = OptimizeOptions { max_iterations: 300, gradient_tolerance: 0.0, energy_tolerance: 0.0, ..Default::default() };
    let m = run(f, &[1.0, -2.0, 0.5, 3.0, -1.0], 1, &opts);
    report_line("lbfgs ||x||", &m);
    monotone(&m);
    assert!(m.report.energy_final < 1e-8, "f = {}", m.report.energy_final);
}

/// Nesterov's Chebyshev-Rosenbrock functions (Gurbuzbalaban & Overton 2012). `squared` selects the
/// first variant `(x0-1)^2/4 + sum |x_{i+1} - 2 x_i^2 + 1|` (one Clarke stationary point, the
/// minimizer 0 at (1, ..., 1)); otherwise the second `|x0-1|/4 + sum |x_{i+1} - 2|x_i| + 1|`, which
/// has 2^(n-1) Clarke stationary points.
fn chebyshev_rosenbrock(squared: bool) -> impl FnMut(&[f64], &mut [f64]) -> Option<f64> {
    move |x: &[f64], g: &mut [f64]| {
        let n = x.len();
        g.fill(0.0);
        let sgn = |v: f64| if v >= 0.0 { 1.0 } else { -1.0 };
        let mut f;
        if squared {
            f = 0.25 * (x[0] - 1.0).powi(2);
            g[0] += 0.5 * (x[0] - 1.0);
        } else {
            f = 0.25 * (x[0] - 1.0).abs();
            g[0] += 0.25 * sgn(x[0] - 1.0);
        }
        for i in 0..n - 1 {
            let (r, dr) = if squared {
                (x[i + 1] - 2.0 * x[i] * x[i] + 1.0, -4.0 * x[i])
            } else {
                (x[i + 1] - 2.0 * x[i].abs() + 1.0, -2.0 * sgn(x[i]))
            };
            f += r.abs();
            g[i + 1] += sgn(r);
            g[i] += dr * sgn(r);
        }
        Some(f)
    }
}

#[test]
fn chebyshev_rosenbrock_nonsmooth() {
    let opts = OptimizeOptions {
        method: Method::Lbfgs { memory: 30 },
        max_iterations: 2000,
        gradient_tolerance: 0.0,
        energy_tolerance: 0.0,
        ..Default::default()
    };
    let mut rng = ddg_core::generators::SplitMix64::new(7);
    let (mut ok, trials) = (0, 20);
    let (mut its, mut evs) = (0, 0);
    for _ in 0..trials {
        let x0: Vec<f64> = (0..3).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let m = run(chebyshev_rosenbrock(true), &x0, 1, &opts);
        monotone(&m);
        its += m.report.iterations;
        evs += m.report.evaluations;
        ok += usize::from(m.report.energy_final < 1e-4);
        // second variant: BFGS typically stalls at a non-minimizing Clarke stationary point
        // (Gurbuzbalaban & Overton); only the descent property is asserted.
        let m2 = run(chebyshev_rosenbrock(false), &x0, 1, &opts);
        monotone(&m2);
    }
    eprintln!("chebyshev-rosenbrock (squared) n=3: {ok}/{trials} reach f < 1e-4; mean iterations {}, evaluations {}", its / trials, evs / trials);
    assert!(ok >= trials * 19 / 20, "{ok} of {trials} reached the minimizer");
}

#[test]
fn undefined_region_is_infinite() {
    // f = x - 0.01 ln x, undefined for x <= 0, minimum at 0.01.
    let f = |x: &[f64], g: &mut [f64]| {
        if x[0] <= 0.0 {
            return None;
        }
        g[0] = 1.0 - 0.01 / x[0];
        Some(x[0] - 0.01 * x[0].ln())
    };
    for method in [Method::Lbfgs { memory: 5 }, Method::GradientDescent] {
        let opts = OptimizeOptions { method, max_iterations: 500, ..Default::default() };
        let m = run(f, &[1.0], 1, &opts);
        report_line(&format!("undefined region {method:?}"), &m);
        monotone(&m);
        assert!((m.x[0] - 0.01).abs() < 1e-5, "x = {}", m.x[0]);
    }
    // start inside the undefined region is an error
    let mut obj = ClosureObjective(f);
    let e = minimize(&mut obj, &[-1.0], &Problem { length_scale: 1.0, group: 1 }, &OptimizeOptions::default(), None);
    assert!(matches!(e, Err(OptimizeError::UndefinedAtStart(_))));
}

#[test]
fn accepted_weak_wolfe_steps_satisfy_s_and_c() {
    for (name, f) in [
        ("rosenbrock", rosenbrock as fn(&[f64], &mut [f64]) -> Option<f64>),
    ] {
        let mut xs = vec![vec![-1.2, 1.0]];
        let mut obj = ClosureObjective(f);
        let mut cb = |p: &ddg_developability::optimize::FlatProgress| {
            xs.push(p.x.to_vec());
            Control::Continue
        };
        let opts = OptimizeOptions { max_iterations: 20, ..Default::default() };
        // c2 = 0.5, c1 = 0 (defaults)
        let _ = minimize(&mut obj, &[-1.2, 1.0], &Problem { length_scale: 1.0, group: 1 }, &opts, Some(&mut cb)).unwrap();
        let mut weak = 0;
        for w in xs.windows(2) {
            let (mut g0, mut g1) = ([0.0; 2], [0.0; 2]);
            let f0 = f(&w[0], &mut g0).unwrap();
            let f1 = f(&w[1], &mut g1).unwrap();
            let d: Vec<f64> = w[1].iter().zip(&w[0]).map(|(a, b)| a - b).collect();
            let s: f64 = g0.iter().zip(&d).map(|(a, b)| a * b).sum();
            let c: f64 = g1.iter().zip(&d).map(|(a, b)| a * b).sum();
            assert!(f1 < f0, "{name}: S fails");
            if c > 0.5 * s {
                continue;
            }
            weak += 1; // a capped/weak step
        }
        assert!(weak <= 2, "{name}: {weak} steps without C");
    }
}

#[test]
fn zero_iterations_and_pinned_and_cancel() {
    let opts = OptimizeOptions { max_iterations: 0, ..Default::default() };
    let m = run(rosenbrock, &[-1.2, 1.0], 1, &opts);
    assert_eq!(m.x, vec![-1.2, 1.0]);
    assert_eq!(m.report.status, OptimizeStatus::MaxIterations);
    assert_eq!(m.report.iterations, 0);

    // pinned: groups of 2 (two "vertices"), pin vertex 1; quadratic pulls both toward the origin
    let q = |x: &[f64], g: &mut [f64]| {
        for (gi, xi) in g.iter_mut().zip(x) {
            *gi = *xi;
        }
        Some(0.5 * x.iter().map(|v| v * v).sum::<f64>())
    };
    let opts = OptimizeOptions { pinned: vec![1], ..Default::default() };
    let x0 = [1.0, 2.0, 3.0, 4.0];
    let m = run(q, &x0, 2, &opts);
    assert_eq!(&m.x[2..], &[3.0, 4.0]);
    assert!(m.x[0].abs() < 1e-6 && m.x[1].abs() < 1e-6);
    let bad = OptimizeOptions { pinned: vec![2], ..Default::default() };
    let mut obj = ClosureObjective(q);
    assert!(matches!(
        minimize(&mut obj, &x0, &Problem { length_scale: 1.0, group: 2 }, &bad, None),
        Err(OptimizeError::PinnedOutOfRange { .. })
    ));

    // cancel after 3 iterations returns the state of iteration 3
    let mut last = vec![];
    let mut obj = ClosureObjective(rosenbrock);
    let mut cb = |p: &ddg_developability::optimize::FlatProgress| {
        last = p.x.to_vec();
        if p.iteration == 3 { Control::Break } else { Control::Continue }
    };
    let m = minimize(&mut obj, &[-1.2, 1.0], &Problem { length_scale: 1.0, group: 1 }, &OptimizeOptions::default(), Some(&mut cb)).unwrap();
    assert_eq!(m.report.status, OptimizeStatus::Cancelled);
    assert_eq!(m.report.iterations, 3);
    assert_eq!(m.x, last);
}

#[test]
fn options_validate() {
    assert!(OptimizeOptions::default().validate().is_ok());
    let bad = |f: fn(&mut OptimizeOptions)| {
        let mut o = OptimizeOptions::default();
        f(&mut o);
        o.validate().is_err()
    };
    assert!(bad(|o| o.method = Method::Lbfgs { memory: 0 }));
    assert!(bad(|o| o.initial_step = 0.0));
    assert!(bad(|o| o.max_displacement = f64::NAN));
    assert!(bad(|o| o.gradient_tolerance = -1.0));
    assert!(bad(|o| o.line_search = LineSearch::WeakWolfe { c1: 0.6, c2: 0.5, max_doublings: 1, max_bisections: 1 }));
    assert!(bad(|o| o.line_search = LineSearch::Backtracking { c1: 1e-4, shrink: 1.0, max_shrinks: 3 }));
}

