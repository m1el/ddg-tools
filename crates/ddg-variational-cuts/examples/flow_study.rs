//! Quantitative experiments with the flow on cuts whose motion is known in closed form.
//!
//! ```text
//! cargo run --release -p ddg-variational-cuts --example flow_study -- circle
//! cargo run --release -p ddg-variational-cuts --example flow_study -- latitude
//! ```

use std::f64::consts::PI;

use ddg_core::generators::{SplitMix64, grid_disk, icosphere};
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CutProblem, Energy, EvalOptions, Extension, Flow, FlowOptions, LengthStep, PenaltyRule,
    Redistance, RegionWeights, SpeedScaling, StepSize,
};

struct Shape {
    mesh: TriMesh,
    geometry: Geometry,
    positions: Vec<[f64; 3]>,
}

fn shape(positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Shape {
    let mesh = TriMesh::new(positions.len(), faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
    Shape {
        mesh,
        geometry,
        positions,
    }
}

/// The unit square with `n × n` cells; interior vertices moved by up to `jitter` cells.
fn square(n: usize, jitter: f64, seed: u64) -> Shape {
    let g = grid_disk(n, n).unwrap();
    let mut positions = g.positions.clone();
    let mut rng = SplitMix64::new(seed);
    let h = 1.0 / n as f64;
    for p in &mut positions {
        let interior =
            p[0] > 0.5 * h && p[0] < 1.0 - 0.5 * h && p[1] > 0.5 * h && p[1] < 1.0 - 0.5 * h;
        if interior {
            p[0] += jitter * h * rng.uniform(-1.0, 1.0);
            p[1] += jitter * h * rng.uniform(-1.0, 1.0);
        }
    }
    shape(positions, g.faces)
}

fn sphere(level: u32) -> Shape {
    let g = icosphere(level).unwrap();
    shape(g.positions, g.faces)
}

/// Mean and standard deviation of `f` over the cut points of the current state.
fn cut_statistics(s: &Shape, flow: &mut Flow, f: impl Fn([f64; 3]) -> f64) -> (f64, f64, usize) {
    let e = flow
        .evaluate(
            &CutProblem::new(Energy::Dirichlet, 0.0),
            &EvalOptions::default(),
        )
        .unwrap();
    let points = e.cut.point_positions(&s.mesh, &s.positions).unwrap();
    let n = points.len();
    if n == 0 {
        return (f64::NAN, f64::NAN, 0);
    }
    let values: Vec<f64> = points.iter().map(|p| f(*p)).collect();
    let mean = values.iter().sum::<f64>() / n as f64;
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
    (mean, var.sqrt(), n)
}

/// A variant of the flow: where the length term is integrated, how the speed is scaled, which
/// matrix extends it.
type Variant = (&'static str, LengthStep, SpeedScaling, Extension);

/// The variants of the flow that are compared.
const VARIANTS: [Variant; 4] = [
    (
        "paper: length step on the mesh, harmonic extension",
        LengthStep::OnMesh,
        SpeedScaling::Unit,
        Extension::Harmonic,
    ),
    (
        "length step on the cut, harmonic extension",
        LengthStep::OnCut,
        SpeedScaling::Unit,
        Extension::Harmonic,
    ),
    (
        "length step on the cut, extension from the cut",
        LengthStep::OnCut,
        SpeedScaling::Unit,
        Extension::Local,
    ),
    (
        "length step on the cut, extension from the cut, speed times mean |grad phi|",
        LengthStep::OnCut,
        SpeedScaling::MeanGradient,
        Extension::Local,
    ),
];

fn options(tau: f64, variant: Variant) -> FlowOptions<'static> {
    // Knobs for experiments: STUDY_DELTA (minimal crossing offset), STUDY_C (smoothing length).
    let env = |name: &str| std::env::var(name).ok().map(|v| v.parse::<f64>().unwrap());
    let mut eval = EvalOptions::default();
    if let Some(delta) = env("STUDY_DELTA") {
        eval.min_crossing_offset = delta;
    }
    if let Some(c) = env("STUDY_C") {
        eval.flux_recovery = ddg_variational_cuts::FluxRecovery::Smoothed { length: c };
    }
    FlowOptions {
        eval,
        length_step: variant.1,
        speed_scaling: variant.2,
        extension: variant.3,
        redistance: Redistance::KeepAdjacent,
        penalty_rule: PenaltyRule::Gradient,
        ..FlowOptions::new(StepSize::Fixed(tau))
    }
}

/// Flat disk, circular cut: `ρ(t) = √(ρ₀² − 2 α_L t)`.
fn circle() {
    let (rho0, alpha, t_end) = (0.3, 1.0, 0.03);
    let center = [0.511, 0.493];
    let exact = |t: f64| (rho0 * rho0 - 2.0 * alpha * t).sqrt();
    println!(
        "flat unit square, circle of radius {rho0}, alpha_L = {alpha}: radius at t = {t_end}, exact {:.6}",
        exact(t_end)
    );
    println!(
        "columns: relative error of the mean radius / roundness (std over mean) / |grad phi| mean at the last step"
    );
    for variant in VARIANTS {
        println!("{}:", variant.0);
        for (n, jitter) in [(32, 0.3), (64, 0.3), (128, 0.3), (256, 0.3), (64, 0.0)] {
            let s = square(n, jitter, 5);
            let h = 1.0 / n as f64;
            let radius =
                |p: [f64; 3]| ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2)).sqrt();
            let phi0: Vec<f64> = s.positions.iter().map(|p| radius(*p) - rho0).collect();
            let mut line = format!("  grid {n:>3}, jitter {jitter}:");
            for steps in [6, 24, 96, 384] {
                if n == 256 && steps == 384 {
                    continue;
                }
                let tau = t_end / steps as f64;
                let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
                let problem = CutProblem::new(Energy::Dirichlet, alpha);
                let opts = options(tau, variant);
                let mut last = None;
                for _ in 0..steps {
                    last = Some(flow.step(&problem, &opts).unwrap());
                }
                let (mean, std, _) = cut_statistics(&s, &mut flow, radius);
                let g = last.unwrap().gradient_norm.map_or(f64::NAN, |g| g.mean);
                line += &format!(
                    "  [{steps} steps, tau = {:.1} h^2] {:+.2e} / {:.1e} / {:.3}",
                    tau / (h * h),
                    mean / exact(t_end) - 1.0,
                    std / mean,
                    g
                );
            }
            println!("{line}");
        }
    }
}

/// Sphere, latitude circle of polar angle `ψ`, `α_L = 0`: `F(ψ(t)) = F(ψ₀) + 4 t / R³` with
/// `F(ψ) = ln(sec ψ + tan ψ) − sin ψ`.
fn latitude() {
    let r = 1.0 / (4.0 * PI).sqrt();
    let f = |psi: f64| (1.0 / psi.cos() + psi.tan()).ln() - psi.sin();
    // Invert F by bisection on (0, π/2).
    let inverse = |value: f64| {
        let (mut lo, mut hi) = (1e-9, 0.5 * PI - 1e-12);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid) < value {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    let psi0 = 0.8;
    let t_end = (f(1.2) - f(psi0)) * r.powi(3) / 4.0; // the exact circle reaches psi = 1.2
    println!(
        "unit-area sphere, latitude circle from psi = {psi0}, alpha_L = 0, t_end = {t_end:.4e} (exact psi = 1.2)"
    );
    for variant in VARIANTS {
        println!("{}:", variant.0);
        for level in [3, 4, 5] {
            let s = sphere(level);
            let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
            let phi0: Vec<f64> = s.positions.iter().map(|p| polar(*p) - psi0).collect();
            for steps in [10, 40, 160, 640] {
                let tau = t_end / steps as f64;
                let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
                let problem = CutProblem::new(Energy::Dirichlet, 0.0);
                let opts = options(tau, variant);
                let mut last = None;
                let mut energies = Vec::new();
                for _ in 0..steps {
                    let report = flow.step(&problem, &opts).unwrap();
                    energies.push(report.energy.total);
                    last = Some(report);
                }
                let (mean, std, count) = cut_statistics(&s, &mut flow, polar);
                let g = last.unwrap().gradient_norm.unwrap();
                let increases = energies.windows(2).filter(|w| w[1] > w[0]).count();
                println!(
                    "  level {level}  steps {:>4}  psi {:.5}  error {:>+9.2e}  spread {:>8.2e}  cut points {:>4}  |grad| {:.3}/{:.3}/{:.3}  energy increases {increases}",
                    steps,
                    mean,
                    mean - inverse(f(psi0) + 4.0 * t_end / r.powi(3)),
                    std,
                    count,
                    g.min,
                    g.mean,
                    g.max
                );
            }
        }
    }
    // The equator is stationary.
    println!("equator, 100 steps of tau = 1e-4:");
    for level in [3, 4, 5] {
        let s = sphere(level);
        let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
        let phi0: Vec<f64> = s
            .positions
            .iter()
            .map(|p| polar(*p) - 0.5 * PI - 1e-3)
            .collect();
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
        let problem = CutProblem::new(Energy::Dirichlet, 0.0);
        let opts = options(1e-4, VARIANTS[2]);
        for _ in 0..100 {
            flow.step(&problem, &opts).unwrap();
        }
        let (mean, std, _) = cut_statistics(&s, &mut flow, polar);
        println!(
            "  level {level}  psi - pi/2 = {:+.3e}  spread {:.2e}",
            mean - 0.5 * PI,
            std
        );
    }
}

/// Sphere, latitude circle with `α_L R > 4`: the stable equilibrium is `sin ψ* = 4/(α_L R)`.
/// One-sided cap with an area reward: `tan(ψ*/2) = R √(−β₀)`.
fn balance() {
    let r = 1.0 / (4.0 * PI).sqrt();
    let alpha = 20.0;
    let exact = (4.0 / (alpha * r)).asin();
    println!("unit-area sphere, alpha_L = {alpha}: equilibrium latitude psi* = {exact:.5}");
    for variant in VARIANTS {
        let name = variant.0;
        for level in [3, 4, 5] {
            let s = sphere(level);
            let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
            for psi0 in [0.5, 1.1] {
                let phi0: Vec<f64> = s.positions.iter().map(|p| polar(*p) - psi0).collect();
                let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
                let problem = CutProblem::new(Energy::Dirichlet, alpha);
                let opts = options(2e-5, variant);
                let mut line = format!("  {name}, level {level} from {psi0}:");
                for block in 0..6 {
                    let mut last = None;
                    for _ in 0..200 {
                        last = Some(flow.step(&problem, &opts).unwrap());
                    }
                    let (mean, std, _) = cut_statistics(&s, &mut flow, polar);
                    let g = last.unwrap().gradient_norm.unwrap();
                    if block % 2 == 1 {
                        line += &format!(
                            "  psi {:.4} (spread {:.1e}, |grad| {:.3}/{:.3}/{:.3})",
                            mean, std, g.min, g.mean, g.max
                        );
                    }
                }
                println!("{line}");
            }
        }
    }
    // One-sided cap: region 0 (the cap) pays for its distortion and is paid for its area;
    // region 1 costs nothing. Growing the cap by a unit of area costs
    // tan²(ψ/2)/R² + β₀ + α_L cot ψ / R, which vanishes at the equilibria.
    let (alpha, target) = (2.0, 1.2f64);
    let beta = -((0.5 * target).tan().powi(2) / (r * r) + alpha / (r * target.tan()));
    println!(
        "one-sided cap, w = (1, 0), beta = ({beta:.4}, 0), alpha_L = {alpha}: stable equilibrium psi* = {target}"
    );
    for level in [3, 4, 5] {
        let s = sphere(level);
        let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
        for psi0 in [1.0, 1.4] {
            let phi0: Vec<f64> = s.positions.iter().map(|p| polar(*p) - psi0).collect();
            let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
            let problem = CutProblem {
                regions: vec![
                    RegionWeights {
                        distortion: 1.0,
                        area: beta,
                    },
                    RegionWeights {
                        distortion: 0.0,
                        area: 0.0,
                    },
                ],
                ..CutProblem::new(Energy::Dirichlet, alpha)
            };
            let opts = FlowOptions::new(StepSize::default());
            let mut line = format!("  level {level} from {psi0}:");
            for _ in 0..4 {
                let mut motion = 0.0;
                for _ in 0..100 {
                    motion = flow.step(&problem, &opts).unwrap().cut_motion;
                }
                let (mean, std, _) = cut_statistics(&s, &mut flow, polar);
                line += &format!("  psi {:.4} ({:.1e}, motion {:.1e})", mean, std, motion);
            }
            println!("{line}");
        }
    }
}

/// The history of one circle run, to see how an instability develops.
fn trace() {
    let n: usize = std::env::args().nth(3).map_or(64, |a| a.parse().unwrap());
    let (rho0, alpha, t_end, steps) = (0.3, 1.0, 0.03, 384);
    let center = [0.511, 0.493];
    let s = square(n, 0.3, 5);
    let radius = |p: [f64; 3]| ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2)).sqrt();
    let phi0: Vec<f64> = s.positions.iter().map(|p| radius(*p) - rho0).collect();
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    let which: usize = std::env::args().nth(2).map_or(2, |a| a.parse().unwrap());
    let opts = options(t_end / steps as f64, VARIANTS[which]);
    let dump: Vec<usize> = std::env::args().nth(4).map_or(Vec::new(), |a| {
        a.split(',').map(|x| x.parse().unwrap()).collect()
    });
    for k in 0..steps {
        let before = flow.level_set();
        let e = flow.evaluate(&problem, &opts.eval).unwrap();
        let (r, stages) = flow.step_with_stages(&problem, &opts).unwrap();
        let g = r.gradient_norm.unwrap();
        if dump.contains(&k) {
            // The face with the largest |grad phi| before this step, and the cut around it.
            let tri = &e.triangles;
            let mut worst = (0.0, 0usize);
            for i in 0..e.cut.segments.len() {
                let t = tri.segment_triangle[i] as usize;
                let a = tri.faces[t][0] as usize;
                let gf = before[a].abs() * tri.segment_length[i] / (2.0 * tri.area[t]);
                if gf > worst.0 {
                    worst = (gf, i);
                }
            }
            let seg = e.cut.segments[worst.1];
            let face = s.mesh.faces()[seg.face as usize];
            println!(
                "  before step {k}: largest |grad| {:.2} in face {} = {:?}",
                worst.0, seg.face, face
            );
            let exact_rho = (rho0 * rho0 - 2.0 * alpha * (flow.time() - r.tau)).sqrt();
            for &v in &face {
                let v = v as usize;
                let nb: Vec<String> = s
                    .mesh
                    .vertex_neighbors(v)
                    .map(|u| {
                        let edge = s.mesh.find_edge(v, u).unwrap();
                        let point = e.cut.edge_point[edge];
                        if point == u32::MAX {
                            format!("{u}:{:+.5}", before[u])
                        } else {
                            let pt = e.cut.points[point as usize];
                            format!(
                                "{u}:{:+.5}[t={:.4} lam={:+.1}]",
                                before[u], pt.t, stages.length_speed[point as usize]
                            )
                        }
                    })
                    .collect();
                println!(
                    "    vertex {v}: exact distance {:+.5}, phi {:+.6}, sigma {:+.2}; neighbours {}",
                    radius(s.positions[v]) - exact_rho,
                    before[v],
                    stages.extension[v],
                    nb.join(" ")
                );
            }
        }
        if k % 8 == 0 || g.max > 3.0 {
            let (mean, std, count) = cut_statistics(&s, &mut flow, radius);
            println!(
                "step {k:>3}  radius {mean:.5} (exact {:.5})  std {std:.1e}  points {count}  |grad| {:.3}/{:.3}/{:.3}  max lambda {:.2}  motion {:.2e}  clamped {}",
                (rho0 * rho0 - 2.0 * alpha * flow.time()).sqrt(),
                g.min,
                g.mean,
                g.max,
                r.max_length_speed,
                r.cut_motion,
                r.n_clamped_crossings
            );
        }
    }
}

/// The amplitude of `cos(m θ)` in the polar angle of the cut, and the rms of what is left after
/// removing the mean and that mode. Points are weighted by the gaps in `θ` to their neighbours.
fn mode_amplitude(s: &Shape, flow: &mut Flow, m: usize) -> (f64, f64, f64) {
    let e = flow
        .evaluate(
            &CutProblem::new(Energy::Dirichlet, 0.0),
            &EvalOptions::default(),
        )
        .unwrap();
    let points = e.cut.point_positions(&s.mesh, &s.positions).unwrap();
    let mut samples: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            (
                p[1].atan2(p[0]),
                (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]),
            )
        })
        .collect();
    samples.sort_by(|a, b| a.0.total_cmp(&b.0));
    let n = samples.len();
    if n < 3 {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    let weight = |k: usize| {
        let before = samples[(k + n - 1) % n].0;
        let after = samples[(k + 1) % n].0;
        let gap = (after - before).rem_euclid(2.0 * PI);
        0.5 * gap / (2.0 * PI)
    };
    let mean: f64 = (0..n).map(|k| weight(k) * samples[k].1).sum();
    let amplitude: f64 = (0..n)
        .map(|k| 2.0 * weight(k) * (samples[k].1 - mean) * (m as f64 * samples[k].0).cos())
        .sum();
    let rest: f64 = (0..n)
        .map(|k| {
            let r = samples[k].1 - mean - amplitude * (m as f64 * samples[k].0).cos();
            weight(k) * r * r
        })
        .sum();
    (mean, amplitude, rest.sqrt())
}

/// Sphere, equator with the perturbation `ε cos(m θ)` of the polar angle. Linearizing the flow
/// about the equator gives the growth rate `(m − 1) [4/R³ − α_L (m + 1)/R²]` of the amplitude:
/// the distortion term amplifies every mode `m ≥ 2` (rate `∝ m`), the length term damps it
/// (rate `∝ m²`). The equator is stable iff `α_L > 4/(3R)`.
fn modes() {
    let r = 1.0 / (4.0 * PI).sqrt();
    let rate = |m: usize, alpha: f64| {
        (m as f64 - 1.0) * (4.0 / r.powi(3) - alpha * (m as f64 + 1.0) / (r * r))
    };
    println!(
        "unit-area sphere, equator + eps cos(m theta), eps = 0.02: growth rate of the amplitude; alpha_L* = 4/(3R) = {:.4}",
        4.0 / (3.0 * r)
    );
    let variants = [VARIANTS[0], VARIANTS[3]];
    for variant in variants {
        println!("{}:", variant.0);
        for alpha in [0.0, 2.0, 4.0, 4.726, 6.0, 10.0] {
            for m in [2usize, 3, 4, 6] {
                let exact = rate(m, alpha);
                let horizon = (0.7 / exact.abs().max(30.0)).min(4e-3);
                let mut line = format!(
                    "  alpha_L {alpha:>5}  m {m}  exact rate {exact:>+8.1}  T {horizon:.2e}  measured (level 3, 4, 5; 40 and 160 steps):"
                );
                for level in [3u32, 4, 5] {
                    let s = sphere(level);
                    let eps = 0.02;
                    let phi0: Vec<f64> = s
                        .positions
                        .iter()
                        .map(|p| {
                            let theta = p[1].atan2(p[0]);
                            let psi = (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
                            psi - 0.5 * PI - eps * (m as f64 * theta).cos()
                        })
                        .collect();
                    for steps in [40usize, 160] {
                        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
                        flow.redistance().unwrap();
                        let (_, a0, _) = mode_amplitude(&s, &mut flow, m);
                        let problem = CutProblem::new(Energy::Dirichlet, alpha);
                        let opts = options(horizon / steps as f64, variant);
                        for _ in 0..steps {
                            flow.step(&problem, &opts).unwrap();
                        }
                        let (_, a1, rest) = mode_amplitude(&s, &mut flow, m);
                        line +=
                            &format!("  {:>+8.1} (rest {:.0e})", (a1 / a0).ln() / horizon, rest);
                    }
                }
                println!("{line}");
            }
        }
    }
}

/// Sphere, latitude circle with a length weight: `ψ′ = cos ψ/(R² sin ψ) · (4/(R sin ψ) − α_L)`.
fn latitude_ode() {
    let r = 1.0 / (4.0 * PI).sqrt();
    let rhs =
        |psi: f64, alpha: f64| psi.cos() / (r * r * psi.sin()) * (4.0 / (r * psi.sin()) - alpha);
    let exact = |psi0: f64, alpha: f64, t: f64| {
        let n = 20000;
        let h = t / n as f64;
        let mut psi = psi0;
        for _ in 0..n {
            let k1 = rhs(psi, alpha);
            let k2 = rhs(psi + 0.5 * h * k1, alpha);
            let k3 = rhs(psi + 0.5 * h * k2, alpha);
            let k4 = rhs(psi + h * k3, alpha);
            psi += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        psi
    };
    println!("unit-area sphere, latitude circle with a length weight, against the ODE (RK4):");
    let variants = [VARIANTS[0], VARIANTS[3]];
    for variant in variants {
        println!("{}:", variant.0);
        for (alpha, psi0, t_end) in [
            (8.0, 0.8, 4e-3),
            (8.0, 0.8, 1.2e-2),
            (20.0, 0.5, 4e-3),
            (20.0, 1.1, 4e-3),
            (20.0, 1.1, 2e-2),
            (0.0, 0.8, 1e-3),
        ] {
            let target = exact(psi0, alpha, t_end);
            let mut line = format!(
                "  alpha_L {alpha:>4}  psi0 {psi0}  t {t_end:.1e}  exact psi {target:.5}  error (spread) at level 3, 4, 5 with 50 and 200 steps:"
            );
            for level in [3u32, 4, 5] {
                let s = sphere(level);
                let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
                let phi0: Vec<f64> = s.positions.iter().map(|p| polar(*p) - psi0).collect();
                for steps in [50usize, 200] {
                    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
                    let problem = CutProblem::new(Energy::Dirichlet, alpha);
                    let opts = options(t_end / steps as f64, variant);
                    for _ in 0..steps {
                        flow.step(&problem, &opts).unwrap();
                    }
                    let (mean, std, _) = cut_statistics(&s, &mut flow, polar);
                    line += &format!("  {:>+8.1e} ({:.0e})", mean - target, std);
                }
            }
            println!("{line}");
        }
    }
}

/// `sphere-trace LEVEL ALPHA PSI0 TAU STEPS [VARIANT]`: one line per step.
fn sphere_trace() {
    let arg = |k: usize| std::env::args().nth(k).unwrap();
    let level: u32 = arg(2).parse().unwrap();
    let alpha: f64 = arg(3).parse().unwrap();
    let psi0: f64 = arg(4).parse().unwrap();
    let tau: f64 = arg(5).parse().unwrap();
    let steps: usize = arg(6).parse().unwrap();
    let variant = std::env::args()
        .nth(7)
        .map_or(3, |v| v.parse::<usize>().unwrap());
    let s = sphere(level);
    let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
    let phi0: Vec<f64> = s.positions.iter().map(|p| polar(*p) - psi0).collect();
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    let opts = options(tau, VARIANTS[variant]);
    println!("{}", VARIANTS[variant].0);
    for step in 0..steps {
        let r = flow.step(&problem, &opts).unwrap();
        let (mean, std, count) = cut_statistics(&s, &mut flow, polar);
        let g = r.gradient_norm.unwrap();
        println!(
            "step {step:>3}  psi {mean:.5} spread {std:.2e} points {count:>4}  E {:.5}  max|s| {:.2e} max|lambda| {:.2e} omega {:.2e}  |grad| {:.3}/{:.3}/{:.3}  motion {:.2e} clamped {}",
            r.energy.total,
            r.max_speed,
            r.max_length_speed,
            r.roughening_rate,
            g.min,
            g.mean,
            g.max,
            r.cut_motion,
            r.n_clamped_crossings
        );
    }
}

/// Number of steps after which the energy was higher than before, and the largest such rise
/// relative to the total drop of the run.
fn energy_rises(report: &ddg_variational_cuts::RunReport) -> (usize, f64) {
    let energies: Vec<f64> = report
        .steps
        .iter()
        .map(|r| r.energy.total)
        .chain([report.final_energy.total])
        .collect();
    let drop = energies[0] - energies[energies.len() - 1];
    let rises: Vec<f64> = energies
        .windows(2)
        .map(|w| w[1] - w[0])
        .filter(|d| *d > 0.0)
        .collect();
    (
        rises.len(),
        rises.iter().fold(0.0f64, |m, d| m.max(*d)) / drop,
    )
}

/// Number of times `values` changes sign along a closed sequence.
fn sign_changes(values: &[f64]) -> usize {
    let n = values.len();
    (0..n)
        .filter(|&k| (values[k] > 0.0) != (values[(k + 1) % n] > 0.0))
        .count()
}

/// Paper Fig. 1: the equator of a sphere with a small perturbation, for three length weights.
/// The linear theory says: the mode `m` grows iff `α_L < 4/((m + 1) R)`, i.e. below 4.73 for
/// `m = 2`, 3.54 for `m = 3`, 2.84 for `m = 4`.
fn fig1() {
    let level: u32 = std::env::args().nth(2).map_or(5, |v| v.parse().unwrap());
    let s = sphere(level);
    let r = 1.0 / (4.0 * PI).sqrt();
    let mut rng = SplitMix64::new(5);
    let phases: Vec<f64> = (0..9).map(|_| rng.uniform(0.0, 2.0 * PI)).collect();
    let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
    let phi0: Vec<f64> = s
        .positions
        .iter()
        .map(|p| {
            let theta = p[1].atan2(p[0]);
            let bump: f64 = (2..=8)
                .map(|m| 0.01 * (m as f64 * theta + phases[m]).cos())
                .sum();
            polar(*p) - 0.5 * PI - bump
        })
        .collect();
    println!(
        "Fig. 1: icosphere({level}), equator (length {:.4}) + 0.01 cos(m theta + phase), m = 2..8",
        2.0 * PI * r
    );
    for alpha in [6.0, 3.0, 1.5] {
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
        flow.redistance().unwrap();
        let problem = CutProblem::new(Energy::Dirichlet, alpha);
        let opts = FlowOptions::new(StepSize::default());
        let stop = ddg_variational_cuts::StopRule {
            max_steps: 3000,
            motion_tolerance: 1e-5,
            energy_tolerance: 0.0,
            energy_window: 1,
        };
        let start = std::time::Instant::now();
        let report = flow.run(&problem, &opts, &stop).unwrap();
        let seconds = start.elapsed().as_secs_f64();
        let first = &report.steps[0].energy;
        let last = &report.final_energy;
        let (increases, largest_rise) = energy_rises(&report);
        let e = flow.evaluate(&problem, &opts.eval).unwrap();
        let lines = e.cut.polylines(&s.mesh, &s.positions).unwrap();
        let heights: Vec<f64> = lines[0]
            .positions
            .iter()
            .map(|p| polar(*p) - 0.5 * PI)
            .collect();
        let largest = heights.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        println!(
            "  alpha_L {alpha}: {} steps ({:?}, {seconds:.1} s, flow time {:.3}), energy {:.4} -> {:.4}, distortion {:.4} -> {:.4}, length {:.4} -> {:.4}, {} piece(s), {} lobes, largest |psi - pi/2| {:.3}, energy went up in {increases} steps (at most {largest_rise:.1e} of the drop), areas {:.4}/{:.4}",
            report.steps.len(),
            report.stop,
            flow.time(),
            first.total,
            last.total,
            first.distortion[0] + first.distortion[1],
            last.distortion[0] + last.distortion[1],
            first.length,
            last.length,
            lines.len(),
            sign_changes(&heights) / 2,
            largest,
            last.area[0],
            last.area[1]
        );
    }
}

/// Paper Fig. 5: a flat disk with a bump; a straight cut beside the bump moves to a longer cut
/// through it.
fn fig5() {
    let n: usize = std::env::args().nth(2).map_or(64, |v| v.parse().unwrap());
    let g = grid_disk(n, n).unwrap();
    let height: f64 = std::env::args().nth(3).map_or(0.05, |v| v.parse().unwrap());
    let width = 0.1;
    let positions: Vec<[f64; 3]> = g
        .positions
        .iter()
        .map(|p| {
            let r2 = (p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2);
            [p[0], p[1], height * (-r2 / (2.0 * width * width)).exp()]
        })
        .collect();
    let s = shape(positions, g.faces);
    let scale = s.geometry.total_area().sqrt();
    println!(
        "Fig. 5: unit square, {n} x {n} cells, bump of height {height} and width {width} in the middle (area {:.4}); cut x = 0.2",
        s.geometry.total_area()
    );
    let x0: f64 = std::env::args().nth(4).map_or(0.2, |v| v.parse().unwrap());
    let motion: f64 = std::env::args().nth(5).map_or(0.25, |v| v.parse().unwrap());
    println!("  start x = {x0}, at most {motion} of an edge per step");
    let phi0: Vec<f64> = s.positions.iter().map(|p| p[0] - x0).collect();
    for (name, energy) in [("Dirichlet", Energy::Dirichlet), ("Hencky", Energy::Hencky)] {
        // The local roughening rate at the top of the bump (the largest one of a cut through
        // it) sets the scale of the length weight: with alpha_L = omega * l / (2 pi), ripples
        // longer than l grow there. The rate of that state (the 0.9-quantile, SPEC D22) is
        // printed next to it: most of the line is flat.
        let through: Vec<f64> = s.positions.iter().map(|p| p[0] - 0.5 - 1e-3).collect();
        let mut probe = Flow::new(&s.mesh, &s.geometry, &through).unwrap();
        let probed = probe
            .step(
                &CutProblem::new(energy.clone(), 0.0),
                &FlowOptions::new(StepSize::Fixed(1e-12)),
            )
            .unwrap();
        let omega_bump = probed.roughening_rate_max;
        println!(
            "  {name}: largest local roughening rate of the cut x = 0.5: {omega_bump:.4} (rate of the state {:.4})",
            probed.roughening_rate
        );
        let wavelengths: Vec<f64> = std::env::var("STUDY_WAVELENGTHS").map_or_else(
            |_| vec![0.8, 0.4, 0.2],
            |v| v.split(',').map(|x| x.parse().unwrap()).collect(),
        );
        for wavelength in wavelengths {
            let alpha = omega_bump * wavelength / (2.0 * PI);
            let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
            flow.redistance().unwrap();
            let problem = CutProblem::new(energy.clone(), alpha);
            // Far from the bump the speed is tiny: the step is limited by the motion only.
            let opts = FlowOptions::new(StepSize::Adaptive {
                time_scale: 0.05,
                max_motion: motion,
                tau_max: 1e3,
            });
            let stop = ddg_variational_cuts::StopRule {
                max_steps: 3000,
                motion_tolerance: 1e-5,
                energy_tolerance: 0.0,
                energy_window: 1,
            };
            if std::env::var("STUDY_TRACE").is_ok() {
                for step in 0..12 {
                    let (r, stages) = flow.step_with_stages(&problem, &opts).unwrap();
                    let (mean, std, _) = cut_statistics(&s, &mut flow, |p| p[0]);
                    let lambda = &stages.length_speed;
                    println!(
                        "    step {step}: tau {:.3e} max|s| {:.3e} max|lambda| {:.3e} lambda ends {:+.3e} {:+.3e} mid {:+.3e} rate min/max {:+.3e}/{:+.3e} x {mean:.5} +- {std:.1e} E {:.6}",
                        r.tau,
                        r.max_speed,
                        r.max_length_speed,
                        lambda[0],
                        lambda[lambda.len() - 1],
                        lambda[lambda.len() / 2],
                        stages.rate.iter().fold(f64::INFINITY, |m, x| m.min(*x)),
                        stages.rate.iter().fold(f64::NEG_INFINITY, |m, x| m.max(*x)),
                        r.energy.total
                    );
                }
            }
            let start = std::time::Instant::now();
            let report = flow.run(&problem, &opts, &stop).unwrap();
            let seconds = start.elapsed().as_secs_f64();
            let first = &report.steps[0].energy;
            let last = &report.final_energy;
            let (increases, largest_rise) = energy_rises(&report);
            let omega = report.steps.last().unwrap().roughening_rate;
            let e = flow.evaluate(&problem, &opts.eval).unwrap();
            let points = e.cut.point_positions(&s.mesh, &s.positions).unwrap();
            let nearest = points
                .iter()
                .map(|p| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt())
                .fold(f64::INFINITY, f64::min);
            let (lo, hi) = points.iter().fold((f64::INFINITY, 0.0f64), |(lo, hi), p| {
                (lo.min(p[0]), hi.max(p[0]))
            });
            println!(
                "  {name}, alpha_L {alpha:.5}: {} steps ({:?}, {seconds:.1} s), energy {:.5} -> {:.5}, distortion {:.5} -> {:.5}, length {:.4} -> {:.4} (file units {:.4} -> {:.4}), cut spans x in [{lo:.3}, {hi:.3}], comes within {nearest:.3} of the top of the bump, {} piece(s), final omega {omega:.3} (ripples longer than {:.3} grow), energy went up in {increases} steps (at most {largest_rise:.1e} of the drop)",
                report.steps.len(),
                report.stop,
                first.total,
                last.total,
                first.distortion[0] + first.distortion[1],
                last.distortion[0] + last.distortion[1],
                first.length,
                last.length,
                first.length * scale,
                last.length * scale,
                e.cut.polylines(&s.mesh, &s.positions).unwrap().len(),
                2.0 * PI * alpha / omega
            );
        }
    }
}

/// Per-stage wall-clock times of one step, on one thread: the median over the steps after the
/// first (the first one also allocates the factorizations).
fn timings() {
    let spot = std::env::args().nth(2);
    let mut shapes: Vec<(String, Shape, Vec<f64>, CutProblem)> = Vec::new();
    for level in [6u32, 7] {
        let s = sphere(level);
        let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
        let phi: Vec<f64> = s.positions.iter().map(|p| polar(*p) - 1.0).collect();
        shapes.push((
            format!("icosphere({level})"),
            s,
            phi,
            CutProblem::new(Energy::Dirichlet, 6.0),
        ));
    }
    if let Some(path) = spot {
        let obj = ddg_core::io::read_obj(&path).unwrap();
        let s = shape(obj.positions.clone(), obj.triangles().unwrap());
        let middle = s.positions.iter().map(|p| p[0]).sum::<f64>() / s.positions.len() as f64;
        let phi: Vec<f64> = s.positions.iter().map(|p| p[0] - middle).collect();
        shapes.push((path, s, phi, CutProblem::new(Energy::Hencky, 0.3)));
    }
    println!(
        "one step, one thread, seconds (median of 5 steps after the first); in brackets the first step"
    );
    for (name, s, phi, problem) in &shapes {
        println!(
            "{name}: {} vertices, {} faces",
            s.mesh.n_vertices(),
            s.mesh.n_faces()
        );
        // The stages of the paper are timed with the flux recovery of SPEC §3.5: as printed
        // (nodal normal derivative) the cut is thrown off the sphere within a few steps.
        let mut paper = FlowOptions::paper(StepSize::Fixed(1e-5));
        paper.eval = EvalOptions::default();
        for (label, options) in [
            ("default stages", FlowOptions::new(StepSize::Fixed(1e-5))),
            ("stages of the paper", paper),
        ] {
            let start = std::time::Instant::now();
            let mut flow = Flow::new(&s.mesh, &s.geometry, phi).unwrap();
            flow.redistance().unwrap();
            let setup = start.elapsed().as_secs_f64();
            let mut rows: Vec<[f64; 5]> = Vec::new();
            let mut points = 0;
            for _ in 0..6 {
                let start = std::time::Instant::now();
                let r = flow.step(problem, &options).unwrap();
                let total = start.elapsed().as_secs_f64();
                let t = r.timings;
                rows.push([t.evaluate, t.length, t.extension, t.redistance, total]);
                points = r.n_cut_points;
            }
            let first = rows[0];
            let median = |k: usize| {
                let mut column: Vec<f64> = rows[1..].iter().map(|r| r[k]).collect();
                column.sort_by(f64::total_cmp);
                column[column.len() / 2]
            };
            println!(
                "  {label} ({points} cut points): setup {setup:.4}; evaluate {:.4} [{:.4}]  length {:.4} [{:.4}]  extension {:.4} [{:.4}]  redistance {:.4} [{:.4}]  step {:.4} [{:.4}]",
                median(0),
                first[0],
                median(1),
                first[1],
                median(2),
                first[2],
                median(3),
                first[3],
                median(4),
                first[4]
            );
        }
    }
}

/// `speeds MESH.obj`: how evenly each energy moves a cut. For three planes through the centroid,
/// the largest and the median `|s|` over the cut points, and the roughening rate.
fn speeds() {
    let path = std::env::args().nth(2).unwrap();
    let obj = ddg_core::io::read_obj(&path).unwrap();
    let s = shape(obj.positions.clone(), obj.triangles().unwrap());
    let n = s.positions.len() as f64;
    let centroid = [0, 1, 2].map(|k| s.positions.iter().map(|p| p[k]).sum::<f64>() / n);
    let mut energies = vec![
        ("Dirichlet".to_string(), Energy::Dirichlet),
        ("Hencky".to_string(), Energy::Hencky),
    ];
    for r_eps in [0.01, 0.03, 0.1, 0.3, 1.0, 10.0] {
        energies.push((
            format!("scaled, r_eps {r_eps}"),
            Energy::ScaledDirichlet(ddg_variational_cuts::ScaleWeight::FeatureSize {
                r_eps,
                smoothing: ddg_variational_cuts::PAPER_FEATURE_SMOOTHING,
            }),
        ));
    }
    println!(
        "{path}: largest |s| / median |s| over the cut (largest |s|, omega), planes x, y, z through the centroid"
    );
    for (name, energy) in &energies {
        let mut line = format!("  {name:<22}");
        for axis in 0..3 {
            let phi: Vec<f64> = s
                .positions
                .iter()
                .map(|p| p[axis] - centroid[axis])
                .collect();
            let mut flow = Flow::new(&s.mesh, &s.geometry, &phi).unwrap();
            flow.redistance().unwrap();
            let problem = CutProblem::new(energy.clone(), 0.0);
            let e = flow.evaluate(&problem, &EvalOptions::default()).unwrap();
            let mut speeds: Vec<f64> = e.speed.iter().map(|x| x.abs()).collect();
            speeds.sort_by(f64::total_cmp);
            let (largest, median) = (speeds[speeds.len() - 1], speeds[speeds.len() / 2]);
            let omega = flow
                .step(&problem, &FlowOptions::new(StepSize::Fixed(1e-12)))
                .unwrap()
                .roughening_rate;
            line += &format!("  {:>8.1} ({largest:.2e}, {omega:.2e})", largest / median);
        }
        println!("{line}");
    }
}

/// `omega MESH.obj [ENERGY]`: how the local roughening rate `ω_c = 2 Σ_r Q_r,c` and the speed
/// are distributed along a planar cut (the start of `cut_flow`): quantiles weighted by the
/// length of the cut at each cut point.
fn omega() {
    let path = std::env::args().nth(2).unwrap();
    let energy = match std::env::args().nth(3).as_deref() {
        Some("dirichlet") => Energy::Dirichlet,
        _ => Energy::Hencky,
    };
    let s = if let Some(level) = path.strip_prefix("icosphere:") {
        sphere(level.parse().unwrap())
    } else {
        let obj = ddg_core::io::read_obj(&path).unwrap();
        shape(obj.positions.clone(), obj.triangles().unwrap())
    };
    let n = s.positions.len() as f64;
    let centroid = [0, 1, 2].map(|k| s.positions.iter().map(|p| p[k]).sum::<f64>() / n);
    let extent = [0, 1, 2].map(|k| {
        let (lo, hi) = s
            .positions
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), p| {
                (lo.min(p[k]), hi.max(p[k]))
            });
        hi - lo
    });
    let axis = (0..3)
        .max_by(|a, b| extent[*a].total_cmp(&extent[*b]))
        .unwrap();
    let phi: Vec<f64> = s
        .positions
        .iter()
        .map(|p| p[axis] - centroid[axis])
        .collect();
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi).unwrap();
    flow.redistance().unwrap();
    let problem = CutProblem::new(energy, 0.0);
    let e = flow.evaluate(&problem, &EvalOptions::default()).unwrap();
    let tri = &e.triangles;
    let n_points = e.cut.points.len();
    let mut dual = vec![0.0; n_points];
    for (segment, length) in e.cut.segments.iter().zip(&tri.segment_length) {
        for p in segment.points {
            dual[p as usize] += 0.5 * length;
        }
    }
    let local: Vec<f64> = (0..n_points)
        .map(|p| {
            2.0 * tri
                .point_vertices(p)
                .iter()
                .map(|&v| e.gradient[v])
                .sum::<f64>()
        })
        .collect();
    let quantiles = |values: &[f64]| {
        let mut order: Vec<usize> = (0..values.len()).collect();
        order.sort_by(|a, b| values[*a].total_cmp(&values[*b]));
        let total: f64 = dual.iter().sum();
        [0.1, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0].map(|q| {
            let mut seen = 0.0;
            for &i in &order {
                seen += dual[i];
                if seen >= q * total {
                    return values[i];
                }
            }
            values[order[order.len() - 1]]
        })
    };
    let total: f64 = dual.iter().sum();
    let mean = |values: &[f64]| values.iter().zip(&dual).map(|(x, l)| x * l).sum::<f64>() / total;
    let speed: Vec<f64> = e.speed.iter().map(|x| x.abs()).collect();
    println!(
        "{path}: {} vertices, {n_points} cut points, cut length {total:.4}; quantiles 10 50 75 90 95 99 100 %, weighted by length",
        s.positions.len()
    );
    let row = |values: &[f64]| {
        quantiles(values)
            .iter()
            .map(|x| format!("{x:>10.3e}"))
            .collect::<Vec<_>>()
            .join(" ")
    };
    println!("  omega_c: {}  mean {:.3e}", row(&local), mean(&local));
    println!("  |s_c|:   {}  mean {:.3e}", row(&speed), mean(&speed));
}

fn main() {
    let what = std::env::args().nth(1).unwrap_or_default();
    match what.as_str() {
        "trace" => trace(),
        "circle" => circle(),
        "latitude" => latitude(),
        "balance" => balance(),
        "modes" => modes(),
        "omega" => omega(),
        "speeds" => speeds(),
        "timings" => timings(),
        "fig1" => fig1(),
        "fig5" => fig5(),
        "sphere-trace" => sphere_trace(),
        "latitude-ode" => latitude_ode(),
        _ => {
            circle();
            latitude();
            balance();
        }
    }
}
