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
        "paper: length step on the mesh, cotan extension",
        LengthStep::OnMesh,
        SpeedScaling::Unit,
        Extension::Cotan,
    ),
    (
        "length step on the cut, cotan extension",
        LengthStep::OnCut,
        SpeedScaling::Unit,
        Extension::Cotan,
    ),
    (
        "length step on the cut, monotone extension",
        LengthStep::OnCut,
        SpeedScaling::Unit,
        Extension::Monotone,
    ),
    (
        "length step on the cut, monotone extension, speed times |grad phi|",
        LengthStep::OnCut,
        SpeedScaling::GradientNorm,
        Extension::Monotone,
    ),
];

fn options(tau: f64, variant: Variant) -> FlowOptions<'static> {
    FlowOptions {
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
    let beta = -4.0;
    let exact = 2.0 * (r * (-beta as f64).sqrt()).atan();
    println!(
        "one-sided cap, w = (1, 0), beta = ({beta}, 0), alpha_L = 0: equilibrium psi* = {exact:.5}"
    );
    for level in [3, 4, 5] {
        let s = sphere(level);
        let polar = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2]);
        for psi0 in [0.5, 1.4] {
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
                ..CutProblem::new(Energy::Dirichlet, 0.0)
            };
            let opts = FlowOptions::new(StepSize::Capped {
                tau_max: 1e-2,
                max_motion: 0.25,
            });
            let mut line = format!("  level {level} from {psi0}:");
            for _ in 0..4 {
                for _ in 0..100 {
                    flow.step(&problem, &opts).unwrap();
                }
                let (mean, std, _) = cut_statistics(&s, &mut flow, polar);
                line += &format!("  psi {:.4} ({:.1e})", mean, std);
            }
            println!("{line}");
        }
    }
}

/// The history of one circle run, to see how an instability develops.
fn trace() {
    let (rho0, alpha, t_end, n, steps) = (0.3, 1.0, 0.03, 64, 384);
    let center = [0.511, 0.493];
    let s = square(n, 0.3, 5);
    let radius = |p: [f64; 3]| ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2)).sqrt();
    let phi0: Vec<f64> = s.positions.iter().map(|p| radius(*p) - rho0).collect();
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    let which: usize = std::env::args().nth(2).map_or(2, |a| a.parse().unwrap());
    let opts = options(t_end / steps as f64, VARIANTS[which]);
    for k in 0..steps {
        let r = flow.step(&problem, &opts).unwrap();
        let g = r.gradient_norm.unwrap();
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

fn main() {
    let what = std::env::args().nth(1).unwrap_or_default();
    match what.as_str() {
        "trace" => trace(),
        "circle" => circle(),
        "latitude" => latitude(),
        "balance" => balance(),
        _ => {
            circle();
            latitude();
            balance();
        }
    }
}
