//! The flow (SPEC §3.8–§3.12, §4): the laws a cut must follow (SPEC §6 items 4 and 5), the
//! controls of SPEC §5, the stop rule, and the helpers for initial cuts.
//!
//! The tolerances are a small multiple of what `examples/flow_study.rs` measures on the same
//! meshes; the tables for finer meshes and smaller steps are in the README.

mod common;

use std::f64::consts::PI;

use common::{Shape, bumpy, polar_angle, sphere, unit_area_radius};
use ddg_core::generators::{SplitMix64, grid_disk};
use ddg_variational_cuts::{
    CutError, CutProblem, Energy, EvalOptions, Flow, FlowOptions, PAPER_FEATURE_SMOOTHING,
    PenaltyRule, PointwiseEnergy, RegionWeights, ScaleWeight, StepSize, StopReason, StopRule, init,
};

/// The unit square with `n × n` cells; interior vertices moved by up to `jitter` cells.
fn square(n: usize, jitter: f64) -> Shape {
    let g = grid_disk(n, n).unwrap();
    let mut positions = g.positions.clone();
    let mut rng = SplitMix64::new(7);
    let h = 1.0 / n as f64;
    for p in &mut positions {
        let inside =
            p[0] > 0.5 * h && p[0] < 1.0 - 0.5 * h && p[1] > 0.5 * h && p[1] < 1.0 - 0.5 * h;
        if inside {
            p[0] += jitter * h * rng.uniform(-1.0, 1.0);
            p[1] += jitter * h * rng.uniform(-1.0, 1.0);
        }
    }
    Shape::new(&format!("square {n}"), positions, g.faces)
}

/// Mean and standard deviation of `f` over the cut points of the current state, and their
/// number.
fn cut_statistics(s: &Shape, flow: &mut Flow, f: impl Fn([f64; 3]) -> f64) -> (f64, f64, usize) {
    let e = flow
        .evaluate(
            &CutProblem::new(Energy::Dirichlet, 0.0),
            &EvalOptions::default(),
        )
        .unwrap();
    let points = e.cut.point_positions(&s.mesh, &s.positions).unwrap();
    let n = points.len();
    let values: Vec<f64> = points.iter().map(|p| f(*p)).collect();
    let mean = values.iter().sum::<f64>() / n as f64;
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
    (mean, var.sqrt(), n)
}

fn fixed(tau: f64) -> FlowOptions<'static> {
    FlowOptions::new(StepSize::Fixed(tau))
}

fn advance(flow: &mut Flow, problem: &CutProblem, options: &FlowOptions<'_>, steps: usize) {
    for _ in 0..steps {
        flow.step(problem, options).unwrap();
    }
}

/// SPEC §6 item 4: on a flat surface `s = 0`, and a circular cut shrinks as
/// `ρ(t) = √(ρ₀² − 2 α_L t)`.
#[test]
fn circle_shrinks_by_the_curve_shortening_law() {
    let s = square(32, 0.3);
    let (rho0, alpha, t_end) = (0.3f64, 1.0, 0.03);
    let exact = (rho0 * rho0 - 2.0 * alpha * t_end).sqrt();
    let radius = |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt();
    let phi0 = s.sample(|p| radius(p) - rho0);
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    let mut errors = Vec::new();
    for steps in [24usize, 96] {
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
        let options = fixed(t_end / steps as f64);
        let mut speed = 0.0f64;
        for _ in 0..steps {
            speed = speed.max(flow.step(&problem, &options).unwrap().max_speed);
        }
        let (mean, std, _) = cut_statistics(&s, &mut flow, radius);
        println!(
            "circle, {steps} steps: radius {mean:.5} (exact {exact:.5}), roundness {:.1e}",
            std / mean
        );
        // The surface is flat: the distortion speed is rounding.
        assert!(speed < 1e-9, "speed {speed}");
        assert!(std / mean < 1e-2, "roundness {}", std / mean);
        errors.push(mean / exact - 1.0);
    }
    // Measured: +9.6e-3 and +2.0e-3. The error is that of the explicit Euler method for the
    // radius, `½ α τ ln(ρ₀/ρ) / ρ²` relative: first order in the step.
    assert!(errors[0].abs() < 1.5e-2, "{errors:?}");
    assert!(errors[1].abs() < 4e-3, "{errors:?}");

    // The stages of the paper do not follow the law: the radius is 17 % too small on this mesh,
    // whatever the step (SPEC D14).
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let mut options = FlowOptions::paper(StepSize::Fixed(t_end / 24.0));
    options.eval = EvalOptions::default();
    advance(&mut flow, &problem, &options, 24);
    let (mean, _, _) = cut_statistics(&s, &mut flow, radius);
    println!("circle, paper stages: radius {mean:.5} (exact {exact:.5})");
    assert!(mean / exact - 1.0 < -5e-2, "{}", mean / exact - 1.0);
}

/// The polar angle of a latitude circle on the sphere of unit area after the time `t`:
/// `ψ′ = cos ψ / (R² sin ψ) · (4 / (R sin ψ) − α_L)`, by the classical Runge–Kutta method.
fn latitude_after(psi0: f64, alpha: f64, t: f64) -> f64 {
    let r = unit_area_radius();
    let rhs = |psi: f64| psi.cos() / (r * r * psi.sin()) * (4.0 / (r * psi.sin()) - alpha);
    let n = 4000;
    let h = t / n as f64;
    let mut psi = psi0;
    for _ in 0..n {
        let k1 = rhs(psi);
        let k2 = rhs(psi + 0.5 * h * k1);
        let k3 = rhs(psi + 0.5 * h * k2);
        let k4 = rhs(psi + h * k3);
        psi += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    psi
}

/// SPEC §6 item 5: a latitude circle on the sphere moves at the speed the exact `s` and the
/// geodesic curvature predict. In the cap of polar angle `ψ`, `∂u/∂n = tan(ψ/2)/R`, so
/// `G_0 − G_1 = −4 cos ψ / (R² sin² ψ)`, and the curvature of the circle is `cot ψ / R`.
#[test]
fn latitude_circle_follows_its_ode() {
    let s = sphere(3);
    let r = unit_area_radius();
    let phi_of = |psi0: f64| s.sample(|p| polar_angle(p) - psi0);

    // With a length weight: toward the equator, slowed down by the curvature of the circle.
    let (alpha, psi0, t_end) = (8.0, 0.8, 4e-3);
    let exact = latitude_after(psi0, alpha, t_end);
    let mut errors = Vec::new();
    for steps in [25usize, 100] {
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi_of(psi0)).unwrap();
        let problem = CutProblem::new(Energy::Dirichlet, alpha);
        advance(&mut flow, &problem, &fixed(t_end / steps as f64), steps);
        let (mean, std, _) = cut_statistics(&s, &mut flow, polar_angle);
        println!(
            "latitude, alpha 8, {steps} steps: psi {mean:.5} (exact {exact:.5}), spread {std:.1e}"
        );
        assert!(std < 3e-3, "spread {std}");
        errors.push(mean - exact);
    }
    // Measured: +4.8e-3 and +1.4e-3 (of a motion by 0.315): first order in the step.
    assert!(errors[0].abs() < 1e-2, "{errors:?}");
    assert!(errors[1].abs() < 3e-3, "{errors:?}");

    // Without a length weight the closed form is F(ψ(t)) = F(ψ₀) + 4t/R³ with
    // F(ψ) = ln(sec ψ + tan ψ) − sin ψ. The flow is then ill-posed (every ripple grows, the
    // shorter the faster; SPEC D18), so the horizon is short and the mesh coarse.
    let f = |psi: f64| (1.0 / psi.cos() + psi.tan()).ln() - psi.sin();
    let t_end = 1e-3;
    let target = f(psi0) + 4.0 * t_end / r.powi(3);
    let (mut lo, mut hi) = (psi0, 0.5 * PI - 1e-9);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let exact = 0.5 * (lo + hi);
    assert!((exact - latitude_after(psi0, 0.0, t_end)).abs() < 1e-9);
    let problem = CutProblem::new(Energy::Dirichlet, 0.0);
    // The stages of the paper, once with the flux recovery of SPEC §3.5 and once as printed.
    let mut paper_stages = FlowOptions::paper(StepSize::Fixed(t_end / 25.0));
    paper_stages.eval = EvalOptions::default();
    let as_printed = FlowOptions::paper(StepSize::Fixed(t_end / 25.0));
    for (name, options, works) in [
        ("default stages", fixed(t_end / 25.0), true),
        ("stages of the paper, recovered flux", paper_stages, true),
        ("stages of the paper as printed", as_printed, false),
    ] {
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi_of(psi0)).unwrap();
        advance(&mut flow, &problem, &options, 25);
        let (mean, std, _) = cut_statistics(&s, &mut flow, polar_angle);
        println!("latitude, alpha 0, {name}: psi {mean:.5} (exact {exact:.5}), spread {std:.1e}");
        if works {
            // Measured: +8.6e-4 and +2.7e-4, of a motion by 0.183.
            assert!((mean - exact).abs() < 2e-3, "{name}: {}", mean - exact);
            assert!(std < 5e-3, "{name}: spread {std}");
        } else {
            // The nodal normal derivative is wrong by orders of magnitude where the cut passes
            // a vertex (SPEC D12), and without a length weight nothing damps the result:
            // measured ψ = 1.336 with a spread of 0.43.
            assert!(std > 0.05, "{name}: spread {std}");
        }
    }
}

/// The amplitude of `cos(m θ)` in the polar angle of the cut.
fn mode_amplitude(s: &Shape, flow: &mut Flow, m: usize) -> f64 {
    let e = flow
        .evaluate(
            &CutProblem::new(Energy::Dirichlet, 0.0),
            &EvalOptions::default(),
        )
        .unwrap();
    let points = e.cut.point_positions(&s.mesh, &s.positions).unwrap();
    let mut samples: Vec<(f64, f64)> = points
        .iter()
        .map(|p| (p[1].atan2(p[0]), polar_angle(*p)))
        .collect();
    samples.sort_by(|a, b| a.0.total_cmp(&b.0));
    let n = samples.len();
    let weight = |k: usize| {
        let gap = (samples[(k + 1) % n].0 - samples[(k + n - 1) % n].0).rem_euclid(2.0 * PI);
        0.5 * gap / (2.0 * PI)
    };
    let mean: f64 = (0..n).map(|k| weight(k) * samples[k].1).sum();
    (0..n)
        .map(|k| 2.0 * weight(k) * (samples[k].1 - mean) * (m as f64 * samples[k].0).cos())
        .sum()
}

/// SPEC §6 item 5 and D18: the equator of the sphere is stationary, and a ripple
/// `ε cos(m θ)` of it grows at the rate `(m − 1) [4/R³ − α_L (m + 1)/R²]`: the distortion term
/// amplifies it, the length term damps it, and the equator is stable iff `α_L > 4/(3R)`.
#[test]
fn ripples_of_the_equator_grow_and_decay_at_the_linear_rates() {
    let s = sphere(4);
    let r = unit_area_radius();
    let rate = |m: usize, alpha: f64| {
        (m as f64 - 1.0) * (4.0 / r.powi(3) - alpha * (m as f64 + 1.0) / (r * r))
    };
    // (α_L, m, relative tolerance, measured rate).
    for (alpha, m, tolerance, measured) in [
        (0.0, 4usize, 0.08, 515.6),
        (2.0, 3, 0.15, 140.0),
        (6.0, 4, 0.30, -474.9),
    ] {
        let exact = rate(m, alpha);
        let horizon = 0.7 / exact.abs();
        let phi0 =
            s.sample(|p| polar_angle(p) - 0.5 * PI - 0.02 * (m as f64 * p[1].atan2(p[0])).cos());
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
        flow.redistance().unwrap();
        let before = mode_amplitude(&s, &mut flow, m);
        let problem = CutProblem::new(Energy::Dirichlet, alpha);
        advance(&mut flow, &problem, &fixed(horizon / 32.0), 32);
        let after = mode_amplitude(&s, &mut flow, m);
        let got = (after / before).ln() / horizon;
        println!(
            "alpha {alpha}, m {m}: rate {got:+.1} (exact {exact:+.1}, measured {measured:+.1})"
        );
        assert!(
            (got - exact).abs() < tolerance * exact.abs(),
            "alpha {alpha}, m {m}: {got} against {exact}"
        );
    }

    // A great circle, with a length weight above 4/(3R) = 4.73: it stays. Its axis is in no
    // special position relative to the mesh (along the equator of the icosphere the cut runs
    // through vertices, where every crossing is moved by the minimal offset).
    let axis = [0.3f64, 0.5, 0.81];
    let norm = axis.iter().map(|x| x * x).sum::<f64>().sqrt();
    let polar_angle = |p: [f64; 3]| {
        ((p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2]) / norm)
            .clamp(-1.0, 1.0)
            .acos()
    };
    let phi0 = s.sample(|p| polar_angle(p) - 0.5 * PI);
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let problem = CutProblem::new(Energy::Dirichlet, 6.0);
    let mut fastest = 0.0f64;
    let mut strongest = 0.0f64;
    for _ in 0..20 {
        let report = flow.step(&problem, &fixed(2e-4)).unwrap();
        fastest = fastest.max(report.max_speed);
        strongest = strongest.max(report.roughening_rate);
    }
    let (mean, std, _) = cut_statistics(&s, &mut flow, polar_angle);
    println!(
        "equator: psi - pi/2 = {:+.2e}, spread {std:.1e}, largest |s| {fastest:.2}, omega {strongest:.1}",
        mean - 0.5 * PI
    );
    // The mesh and the circle are symmetric under the antipodal map, so the mean is π/2 by
    // symmetry; what the circle does is in the spread. Measured: 1.0e-3 (the mesh width is
    // 0.075 in this angle).
    assert!((mean - 0.5 * PI).abs() < 1e-9, "{}", mean - 0.5 * PI);
    assert!(std < 3e-3, "spread {std}");
    // |s| is a discretization error there: G_0 = G_1 = 1/R² = 12.57 each.
    assert!(fastest < 0.05 * 2.0 / (r * r), "largest speed {fastest}");
    // ω = 2 (G_0 + G_1) = 4/R².
    assert!(
        (strongest - 4.0 / (r * r)).abs() < 0.05 * 4.0 / (r * r),
        "omega {strongest}"
    );
}

/// SPEC §6 item 5: the energy decreases along the flow for small steps, in every energy, on a
/// closed mesh and on a mesh with boundary that the cut ends on.
#[test]
fn energy_decreases_for_small_steps_in_every_energy() {
    let ball = sphere(3);
    let disk = bumpy(24, 0.08);
    let ball_phi = ball.sample(|p| polar_angle(p) - 1.0 - 0.1 * (3.0 * p[1].atan2(p[0])).cos());
    let disk_phi = disk.sample(|p| p[0] - 0.42 - 0.1 * (2.0 * PI * p[1]).sin());
    let energies = [
        ("Dirichlet", Energy::Dirichlet),
        (
            "scaled Dirichlet",
            Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
                r_eps: 0.3,
                smoothing: PAPER_FEATURE_SMOOTHING,
            }),
        ),
        ("Hencky", Energy::Hencky),
        (
            "pointwise",
            Energy::Pointwise(PointwiseEnergy::new(
                |u| u * u + u.powi(4),
                |u| 2.0 * u + 4.0 * u.powi(3),
            )),
        ),
    ];
    for (shape, phi0) in [(&ball, &ball_phi), (&disk, &disk_phi)] {
        for (name, energy) in &energies {
            let mut flow = Flow::new(&shape.mesh, &shape.geometry, phi0).unwrap();
            // The length weight is set from the roughening rate of the start state so that the
            // length term controls ripples shorter than a tenth of the surface (SPEC D18).
            let omega = flow
                .roughening_rate(
                    &CutProblem::new(energy.clone(), 0.0),
                    &EvalOptions::default(),
                )
                .unwrap();
            let alpha = ddg_variational_cuts::length_weight_for_ripples(omega, 0.1);
            assert!((alpha - omega * 0.1 / (2.0 * PI)).abs() < 1e-15 * omega);
            let problem = CutProblem::new(energy.clone(), alpha);
            let options = FlowOptions::new(StepSize::Capped {
                tau_max: 0.02 / omega,
                max_motion: 0.1,
            });
            let report = flow
                .run(
                    &problem,
                    &options,
                    &StopRule {
                        max_steps: 20,
                        motion_tolerance: 0.0,
                        energy_tolerance: 0.0,
                        energy_window: 1,
                    },
                )
                .unwrap();
            let energy_of: Vec<f64> = report
                .steps
                .iter()
                .map(|r| r.energy.total)
                .chain([report.final_energy.total])
                .collect();
            let drop = energy_of[0] - energy_of[energy_of.len() - 1];
            let worst = energy_of
                .windows(2)
                .map(|w| w[1] - w[0])
                .fold(f64::NEG_INFINITY, f64::max);
            println!(
                "{}, {name}: alpha {alpha:.3e}, energy {:.6e} -> {:.6e}, largest change of a step {:+.2e} of the total drop",
                shape.name,
                energy_of[0],
                energy_of[energy_of.len() - 1],
                worst / drop
            );
            assert!(drop > 0.0, "{}, {name}: {energy_of:?}", shape.name);
            assert!(worst <= 0.0, "{}, {name}: {energy_of:?}", shape.name);
        }
    }
}

/// For every vertex of the sphere, the index of the antipodal vertex.
fn antipodes(s: &Shape) -> Vec<usize> {
    let key = |p: [f64; 3]| p.map(|x| (x * 1e6).round() as i64);
    let index: std::collections::HashMap<[i64; 3], usize> = s
        .positions
        .iter()
        .enumerate()
        .map(|(i, p)| (key(*p), i))
        .collect();
    s.positions
        .iter()
        .map(|p| index[&key([-p[0], -p[1], -p[2]])])
        .collect()
}

/// SPEC §6 item 5, §3.11: a symmetric start stays symmetric, to rounding without symmetry
/// classes and exactly with them; classes make an asymmetric start symmetric.
#[test]
fn symmetric_input_stays_symmetric() {
    let s = sphere(3);
    let opposite = antipodes(&s);
    let problem = CutProblem::new(Energy::Dirichlet, 6.0);
    let asymmetry = |flow: &Flow| {
        let phi = flow.level_set();
        (0..phi.len())
            .map(|i| (phi[i] - phi[opposite[i]]).abs())
            .fold(0.0, f64::max)
    };
    // Region 0 is a band around a great circle; the cut is two circles, images of each other
    // under the antipodal map, which is a symmetry of the mesh.
    let band = |axis: [f64; 3]| {
        let norm = axis.iter().map(|x| x * x).sum::<f64>().sqrt();
        s.sample(|p| {
            let cosine = (p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2]) / norm;
            (cosine.clamp(-1.0, 1.0).acos() - 0.5 * PI).abs() - 0.45
        })
    };

    let phi0 = band([0.3, 0.5, 0.81]);
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    advance(&mut flow, &problem, &fixed(1e-4), 10);
    let (_, _, n) = cut_statistics(&s, &mut flow, polar_angle);
    assert!(n > 0);
    println!(
        "generic axis, without classes: asymmetry {:.1e}",
        asymmetry(&flow)
    );
    assert!(asymmetry(&flow) < 1e-12);

    // Around the axis of the icosphere the circles are parallel to mesh edges, the two
    // diagonals of the quads they cut off are equally long, and the rule for the tie
    // (SPEC §3.3, D13c) refers to the orientation of the face, which the antipodal map
    // reverses: the two halves are triangulated differently, and the symmetry holds to the
    // discretization error only. Measured: 1.2e-5 after 10 steps.
    let phi0 = band([0.0, 0.0, 1.0]);
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let ties = flow
        .evaluate(&problem, &EvalOptions::default())
        .unwrap()
        .report
        .n_diagonal_near_ties;
    advance(&mut flow, &problem, &fixed(1e-4), 10);
    println!(
        "axis of the mesh, without classes: asymmetry {:.1e} ({ties} near ties of the diagonal)",
        asymmetry(&flow)
    );
    assert!(ties > 0);
    assert!(asymmetry(&flow) < 1e-4);

    // With classes the symmetry is exact, ties or not.
    let classes: Vec<u32> = (0..opposite.len())
        .map(|i| i.min(opposite[i]) as u32)
        .collect();
    let mut options = fixed(1e-4);
    options.symmetry = Some(classes);
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    advance(&mut flow, &problem, &options, 10);
    assert_eq!(asymmetry(&flow), 0.0);

    // An asymmetric start: after one step with classes the function is symmetric.
    let tilted: Vec<f64> = phi0
        .iter()
        .zip(&s.positions)
        .map(|(x, p)| x + 0.05 * p[2])
        .collect();
    let mut flow = Flow::new(&s.mesh, &s.geometry, &tilted).unwrap();
    assert!(asymmetry(&flow) > 1e-2);
    advance(&mut flow, &problem, &options, 1);
    assert_eq!(asymmetry(&flow), 0.0);
}

/// SPEC §6 item 5, §2: one step does not depend on the numbering of vertices and faces, and
/// scaling the mesh scales the level-set function.
#[test]
fn one_step_is_invariant_under_relabeling_and_covariant_under_scaling() {
    let mut g = ddg_core::generators::bumpy_grid_disk(14, 14, 0.1).unwrap();
    ddg_core::generators::perturb(&mut g.positions, [0.01, 0.01, 0.01], 3);
    let s = Shape::from_generated("bumpy", g);
    let phi0 = s.sample(|p| ((p[0] - 0.45).powi(2) + (p[1] - 0.55).powi(2)).sqrt() - 0.27);
    let problem = CutProblem {
        penalty: Some(s.sample(|p| 0.5 + p[0])),
        ..CutProblem::new(Energy::Hencky, 0.4)
    };
    for options in [fixed(2e-4), FlowOptions::paper(StepSize::Fixed(2e-4))] {
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
        advance(&mut flow, &problem, &options, 3);
        let reference = flow.level_set();
        let size = reference.iter().fold(0.0f64, |m, x| m.max(x.abs()));

        // A permutation of the vertices, the faces in another order and rotated.
        let n = s.positions.len();
        let mut rng = SplitMix64::new(11);
        let mut new_index: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            new_index.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let mut positions = vec![[0.0; 3]; n];
        let mut phi = vec![0.0; n];
        let mut penalty = vec![0.0; n];
        for i in 0..n {
            positions[new_index[i]] = s.positions[i];
            phi[new_index[i]] = phi0[i];
            penalty[new_index[i]] = problem.penalty.as_ref().unwrap()[i];
        }
        let mut faces: Vec<[u32; 3]> = s
            .mesh
            .faces()
            .iter()
            .enumerate()
            .map(|(f, face)| {
                let k = f % 3;
                [0, 1, 2].map(|c| new_index[face[(c + k) % 3] as usize] as u32)
            })
            .collect();
        faces.reverse();
        let relabeled = Shape::new("relabeled", positions, faces);
        let mut flow = Flow::new(&relabeled.mesh, &relabeled.geometry, &phi).unwrap();
        let permuted = CutProblem {
            penalty: Some(penalty),
            ..problem.clone()
        };
        advance(&mut flow, &permuted, &options, 3);
        let got = flow.level_set();
        let worst = (0..n)
            .map(|i| (got[new_index[i]] - reference[i]).abs())
            .fold(0.0, f64::max);
        println!("relabeling: largest difference {:.1e} of {size:.2}", worst);
        assert!(worst < 1e-9 * size, "{worst}");

        // The same surface three times as large: the parameters refer to the unit-area surface,
        // so the function is three times as large.
        let scaled = Shape::new(
            "scaled",
            s.positions.iter().map(|p| p.map(|x| 3.0 * x)).collect(),
            s.mesh.faces().to_vec(),
        );
        let phi: Vec<f64> = phi0.iter().map(|x| 3.0 * x).collect();
        let mut flow = Flow::new(&scaled.mesh, &scaled.geometry, &phi).unwrap();
        advance(&mut flow, &problem, &options, 3);
        let got = flow.level_set();
        let worst = (0..n)
            .map(|i| (got[i] - 3.0 * reference[i]).abs())
            .fold(0.0, f64::max);
        println!(
            "scaling: largest difference {:.1e} of {:.2}",
            worst,
            3.0 * size
        );
        assert!(worst < 1e-9 * size, "{worst}");
    }
}

/// SPEC §3.10, D8: a straight cut across a penalty that grows linearly, `η = A x`, on a flat
/// surface. The energy of the cut at `x` is `(α_L + A x)` times its length, so its gradient
/// flow moves the cut down the penalty at the speed `A`. The rule of the paper keeps
/// `(α_L + η) κ` only, and the curvature of a straight cut is zero: it does not move.
#[test]
fn penalty_moves_the_cut_with_the_gradient_rule_only() {
    let s = square(24, 0.3);
    let (slope, alpha, x0, t_end) = (2.0, 0.5, 0.6, 0.05);
    let phi0 = s.sample(|p| p[0] - x0);
    let problem = CutProblem {
        penalty: Some(s.sample(|p| slope * p[0])),
        ..CutProblem::new(Energy::Dirichlet, alpha)
    };
    let steps = 25;
    for (name, base) in [
        (
            "default stages",
            FlowOptions::new(StepSize::Fixed(t_end / steps as f64)),
        ),
        (
            "stages of the paper",
            FlowOptions::paper(StepSize::Fixed(t_end / steps as f64)),
        ),
    ] {
        for rule in [PenaltyRule::Gradient, PenaltyRule::PaperDiagonal] {
            let mut options = base.clone();
            options.penalty_rule = rule;
            let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
            let before = flow
                .evaluate(&problem, &options.eval)
                .unwrap()
                .energy
                .penalty;
            advance(&mut flow, &problem, &options, steps);
            let after = flow
                .evaluate(&problem, &options.eval)
                .unwrap()
                .energy
                .penalty;
            let (mean, std, _) = cut_statistics(&s, &mut flow, |p| p[0]);
            println!(
                "{name}, {rule:?}: cut at x = {mean:.4} (spread {std:.1e}), penalty {before:.4} -> {after:.4}"
            );
            match rule {
                PenaltyRule::Gradient => {
                    // Measured: 0.5000 (default stages), 0.5033 (stages of the paper).
                    assert!((mean - (x0 - slope * t_end)).abs() < 0.01, "{name}: {mean}");
                    assert!(std < 5e-3, "{name}: spread {std}");
                    assert!(after < before - 0.15, "{name}: {before} -> {after}");
                }
                PenaltyRule::PaperDiagonal => {
                    // Exactly at rest with the length terms on the cut. On the mesh, the
                    // implicit step flattens the function at the two sides of the square that
                    // the cut does not meet, and the cut creeps by 0.004.
                    let tolerance = if name == "default stages" { 1e-9 } else { 1e-2 };
                    assert!((mean - x0).abs() < tolerance, "{name}: {mean}");
                    assert!(
                        (after - before).abs() < 2.0 * tolerance,
                        "{name}: {before} -> {after}"
                    );
                }
            }
        }
    }
}

/// Paper Fig. 14, SPEC §5: distortion is counted on one side only, and that side is paid for
/// its area. On the sphere a cap of polar angle `ψ` then gains
/// `tan²(ψ/2)/R² + β₀ + α_L cot ψ / R` per unit of area it grows by; with `β₀ < 0` chosen so
/// that this vanishes at `ψ* = 1.2`, a cap grows from `ψ = 1` and stops there, and a larger
/// one shrinks to it.
#[test]
fn one_sided_patch_with_an_area_reward_grows_and_stops() {
    let s = sphere(3);
    let r = unit_area_radius();
    let (alpha, target) = (2.0, 1.2f64);
    let beta = -((0.5 * target).tan().powi(2) / (r * r) + alpha / (r * target.tan()));
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
    let stop = StopRule {
        max_steps: 400,
        motion_tolerance: 5e-5,
        energy_tolerance: 0.0,
        energy_window: 1,
    };
    for psi0 in [1.0, 1.4] {
        let phi0 = s.sample(|p| polar_angle(p) - psi0);
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
        let report = flow
            .run(&problem, &FlowOptions::new(StepSize::default()), &stop)
            .unwrap();
        let (mean, std, _) = cut_statistics(&s, &mut flow, polar_angle);
        let first = report.steps[0].energy.total;
        println!(
            "one-sided patch from {psi0}: psi {mean:.4} (spread {std:.1e}) after {} steps, {:?}, energy {first:.5} -> {:.5}, area of the patch {:.4}",
            report.steps.len(),
            report.stop,
            report.final_energy.total,
            report.final_energy.area[0]
        );
        // Measured: 1.165 from below, 1.222 from above; the mesh width is 0.15 in this angle.
        // On icosphere(4) and (5): 1.202 and 1.205, 1.195 and 1.201.
        assert!((mean - target).abs() < 0.05, "from {psi0}: {mean}");
        assert!(std < 1e-2, "from {psi0}: spread {std}");
        assert_eq!(report.stop, StopReason::Converged, "from {psi0}");
        assert!(report.final_energy.total < first);
        // The other side is not solved for nothing: it is reported, and costs nothing.
        assert!(report.final_energy.distortion[1] > 0.0);
    }
}

/// SPEC §4: the three ways a run ends, and the reports of a run.
#[test]
fn a_run_stops_for_the_reason_it_reports() {
    // A small circle on a flat surface vanishes at t = ρ₀²/(2 α_L) = 0.0098.
    let s = square(16, 0.3);
    let radius = |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt();
    let phi0 = s.sample(|p| radius(p) - 0.14);
    let problem = CutProblem::new(Energy::Dirichlet, 1.0);
    let options = fixed(5e-4);
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let stop = StopRule {
        max_steps: 100,
        motion_tolerance: 0.0,
        energy_tolerance: 0.0,
        energy_window: 1,
    };
    let report = flow.run(&problem, &options, &stop).unwrap();
    println!(
        "vanishing circle: {} steps, time {:.4}",
        report.steps.len(),
        flow.time()
    );
    assert_eq!(report.stop, StopReason::CutVanished);
    assert!(report.steps.last().unwrap().cut_vanished);
    assert!(report.steps.len() < 40 && report.steps.len() > 10);
    assert_eq!(report.final_energy.length, 0.0);
    assert_eq!(flow.n_steps(), report.steps.len());
    // Times add up, energies are those of the start states, and a further step does nothing.
    for pair in report.steps.windows(2) {
        assert!((pair[0].time + pair[0].tau - pair[1].time).abs() < 1e-15);
        assert!(pair[1].energy.length < pair[0].energy.length);
    }
    let after = flow.level_set();
    let idle = flow.step(&problem, &options).unwrap();
    assert!(idle.cut_vanished && idle.tau == 0.0 && idle.n_cut_points == 0);
    assert_eq!(flow.level_set(), after);

    // The limit on the number of steps.
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let report = flow
        .run(
            &problem,
            &options,
            &StopRule {
                max_steps: 3,
                motion_tolerance: 0.0,
                energy_tolerance: 0.0,
                energy_window: 1,
            },
        )
        .unwrap();
    assert_eq!((report.stop, report.steps.len()), (StopReason::MaxSteps, 3));

    // The energy criterion: with a tolerance of 1000 % over two steps it ends the run as soon
    // as it can compare.
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let report = flow
        .run(
            &problem,
            &options,
            &StopRule {
                max_steps: 100,
                motion_tolerance: 0.0,
                energy_tolerance: 10.0,
                energy_window: 2,
            },
        )
        .unwrap();
    assert_eq!(
        (report.stop, report.steps.len()),
        (StopReason::EnergyStalled, 3)
    );

    // Convergence: a straight cut across a flat square does not move.
    let phi0 = s.sample(|p| p[0] - 0.47);
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let report = flow.run(&problem, &options, &StopRule::default()).unwrap();
    assert_eq!(
        (report.stop, report.steps.len()),
        (StopReason::Converged, 1)
    );
    assert!(report.steps[0].cut_motion < 1e-12);
}

/// SPEC §4, §5: the state is the level-set function; a flow can be resumed from a stored
/// function, and every control can change between steps.
#[test]
fn a_flow_resumes_from_its_state_and_takes_new_controls_between_steps() {
    let s = sphere(2);
    let phi0 = s.sample(|p| polar_angle(p) - 1.0);
    let problem = CutProblem::new(Energy::Dirichlet, 6.0);
    let options = fixed(1e-4);
    let mut whole = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    advance(&mut whole, &problem, &options, 6);

    let mut first = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    advance(&mut first, &problem, &options, 3);
    let mut second = Flow::new(&s.mesh, &s.geometry, &first.level_set()).unwrap();
    advance(&mut second, &problem, &options, 3);
    let worst = whole
        .level_set()
        .iter()
        .zip(second.level_set())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(worst < 1e-12, "{worst}");
    assert_eq!((whole.n_steps(), second.n_steps()), (6, 3));
    assert!((whole.time() - 6e-4).abs() < 1e-15);

    // The schedule of a caller: every step with other controls.
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let n = s.positions.len();
    let controls = [
        CutProblem::new(Energy::Dirichlet, 6.0),
        CutProblem::new(Energy::Hencky, 0.5),
        CutProblem {
            penalty: Some(s.sample(|p| 1.0 + p[2])),
            ..CutProblem::new(
                Energy::ScaledDirichlet(ScaleWeight::PerVertex(vec![0.5; n])),
                3.0,
            )
        },
        CutProblem {
            curvature: ddg_variational_cuts::CurvatureSource::Uniform(4.0 * PI),
            regions: vec![
                RegionWeights {
                    distortion: 2.0,
                    area: -1.0,
                },
                RegionWeights {
                    distortion: 0.5,
                    area: 0.3,
                },
            ],
            ..CutProblem::new(Energy::Dirichlet, 6.0)
        },
    ];
    let steps = [
        StepSize::Fixed(1e-4),
        StepSize::Capped {
            tau_max: 1e-3,
            max_motion: 0.2,
        },
    ];
    for (k, problem) in controls.iter().enumerate() {
        let mut options = FlowOptions::new(steps[k % 2]);
        if k == 2 {
            options.penalty_rule = PenaltyRule::PaperDiagonal;
        }
        let report = flow.step(problem, &options).unwrap();
        assert!(report.tau > 0.0 && report.cut_motion > 0.0 && report.n_cut_points > 0);
        assert!(report.energy.total.is_finite());
    }

    // The matrix of the paper's implicit step is factorized once and reused while the step
    // and the weights stay the same.
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    let problem = CutProblem::new(Energy::Dirichlet, 6.0);
    let mut options = FlowOptions::paper(StepSize::Fixed(1e-4));
    let factorized = |flow: &mut Flow, problem: &CutProblem, options: &FlowOptions<'_>| {
        flow.step(problem, options).unwrap().length_factorized
    };
    assert!(factorized(&mut flow, &problem, &options));
    assert!(!factorized(&mut flow, &problem, &options));
    options.step_size = StepSize::Fixed(2e-4);
    assert!(factorized(&mut flow, &problem, &options));
    assert!(!factorized(&mut flow, &problem, &options));
    assert!(factorized(
        &mut flow,
        &CutProblem::new(Energy::Dirichlet, 5.0),
        &options
    ));
    options.penalty_rule = PenaltyRule::PaperDiagonal;
    assert!(factorized(
        &mut flow,
        &CutProblem::new(Energy::Dirichlet, 5.0),
        &options
    ));
}

/// Bad input to the flow and to the helpers is an error that says what is wrong.
#[test]
fn bad_input_is_an_error() {
    let s = sphere(1);
    let n = s.positions.len();
    let phi0 = s.sample(|p| polar_angle(p) - 1.0);
    assert!(matches!(
        Flow::new(&s.mesh, &s.geometry, &phi0[1..]),
        Err(CutError::WrongLength { .. })
    ));
    let mut bad = phi0.clone();
    bad[3] = f64::NAN;
    assert!(matches!(
        Flow::new(&s.mesh, &s.geometry, &bad),
        Err(CutError::NonFiniteLevelSet { vertex: 3, .. })
    ));
    let other = sphere(2);
    assert!(matches!(
        Flow::new(&s.mesh, &other.geometry, &phi0),
        Err(CutError::Geometry(_))
    ));

    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
    assert!(matches!(
        flow.set_level_set(&bad),
        Err(CutError::NonFiniteLevelSet { .. })
    ));
    let problem = CutProblem::new(Energy::Dirichlet, 1.0);
    for step_size in [
        StepSize::Fixed(0.0),
        StepSize::Fixed(-1.0),
        StepSize::Fixed(f64::NAN),
        StepSize::Capped {
            tau_max: 0.0,
            max_motion: 0.5,
        },
        StepSize::Capped {
            tau_max: 1e-3,
            max_motion: f64::INFINITY,
        },
        StepSize::Adaptive {
            time_scale: -0.05,
            max_motion: 0.25,
            tau_max: 1.0,
        },
    ] {
        let error = flow
            .step(&problem, &FlowOptions::new(step_size))
            .unwrap_err();
        assert!(matches!(error, CutError::InvalidOptions { .. }), "{error}");
        assert!(!error.to_string().is_empty());
    }
    let mut options = fixed(1e-4);
    options.symmetry = Some(vec![0; n - 1]);
    assert!(matches!(
        flow.step(&problem, &options),
        Err(CutError::WrongLength { .. })
    ));
    // A failed step leaves the state alone.
    assert_eq!(flow.n_steps(), 0);
    let before = flow.level_set();
    let _ = flow.step(&CutProblem::new(Energy::Dirichlet, -1.0), &fixed(1e-4));
    assert_eq!(flow.level_set(), before);
    // The rule of the paper divides by α_L + η.
    let zero = CutProblem {
        penalty: Some((0..n).map(|i| (i % 2) as f64).collect()),
        ..CutProblem::new(Energy::Dirichlet, 0.0)
    };
    let mut options = FlowOptions::paper(StepSize::Fixed(1e-4));
    options.penalty_rule = PenaltyRule::PaperDiagonal;
    let error = flow.step(&zero, &options).unwrap_err();
    assert!(matches!(error, CutError::InvalidOptions { .. }), "{error}");

    // Helpers.
    assert!(matches!(
        init::geodesic_circle(&s.mesh, &s.geometry, n, 0.1),
        Err(CutError::IndexOutOfRange { .. })
    ));
    for radius in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            init::geodesic_circle(&s.mesh, &s.geometry, 0, radius),
            Err(CutError::InvalidOptions { .. })
        ));
    }
    assert!(matches!(
        init::geodesic_circle(&s.mesh, &other.geometry, 0, 0.1),
        Err(CutError::Geometry(_))
    ));
    assert!(matches!(
        init::from_vertex_labels(&s.mesh, &s.geometry, &vec![0; n + 1]),
        Err(CutError::WrongLength { .. })
    ));
    let mut labels = vec![0u32; n];
    labels[5] = 2;
    let error = init::from_vertex_labels(&s.mesh, &s.geometry, &labels).unwrap_err();
    assert!(error.to_string().contains("vertex 5"), "{error}");
    assert!(matches!(
        init::from_face_labels(&s.mesh, &s.geometry, &[0, 1]),
        Err(CutError::WrongLength { .. })
    ));
    let mut labels = vec![0u32; s.mesh.n_faces()];
    labels[7] = 3;
    let error = init::from_face_labels(&s.mesh, &s.geometry, &labels).unwrap_err();
    assert!(error.to_string().contains("face 7"), "{error}");
    assert!(matches!(
        init::redistanced(&s.mesh, &s.geometry, &bad),
        Err(CutError::NonFiniteLevelSet { vertex: 3, .. })
    ));
    assert!(matches!(
        init::redistanced(&s.mesh, &s.geometry, &phi0[1..]),
        Err(CutError::WrongLength { .. })
    ));
    let mut visibility = vec![1.0; n];
    visibility[2] = -0.5;
    assert!(matches!(
        init::visibility_penalty(&s.mesh, &s.positions, [0.0, 0.0, 5.0], &visibility),
        Err(CutError::InvalidPenalty { vertex: 2, .. })
    ));
    assert!(matches!(
        init::visibility_penalty(&s.mesh, &s.positions, [0.0, f64::NAN, 5.0], &vec![1.0; n]),
        Err(CutError::InvalidOptions { .. })
    ));
    assert!(matches!(
        init::visibility_penalty(&s.mesh, &s.positions[1..], [0.0, 0.0, 5.0], &vec![1.0; n]),
        Err(CutError::WrongLength { .. })
    ));
}

/// SPEC §5: the helpers for initial cuts, with their rules for ties.
#[test]
fn helpers_for_initial_cuts() {
    // A geodesic disk: on a flat regular grid, around a vertex.
    let s = square(32, 0.0);
    let center = s
        .positions
        .iter()
        .position(|p| (p[0] - 0.5).abs() < 1e-12 && (p[1] - 0.5).abs() < 1e-12)
        .unwrap();
    let phi = init::geodesic_circle(&s.mesh, &s.geometry, center, 0.25).unwrap();
    assert_eq!(phi[center], -0.25);
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi).unwrap();
    let radius = |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt();
    let (mean, std, _) = cut_statistics(&s, &mut flow, radius);
    println!("geodesic disk of radius 0.25: cut at {mean:.4}, spread {std:.1e}");
    // Fast marching overestimates distances, so the disk is a little smaller.
    assert!(mean < 0.25 && mean > 0.235, "{mean}");
    assert!(std < 1e-2);
    // A vertex at exactly the radius is inside (along the grid axes distances are exact).
    let on_circle = s
        .positions
        .iter()
        .position(|p| (p[0] - 0.75).abs() < 1e-12 && (p[1] - 0.5).abs() < 1e-12)
        .unwrap();
    assert!(
        phi[on_circle].abs() < 1e-12
            && ddg_variational_cuts::region_of(phi[on_circle].min(0.0)) == 0
    );

    // Vertex labels: the cut crosses exactly the edges whose ends differ; on congruent
    // triangles at their midpoints.
    let labels: Vec<u32> = s.positions.iter().map(|p| u32::from(p[0] > 0.52)).collect();
    let phi = init::from_vertex_labels(&s.mesh, &s.geometry, &labels).unwrap();
    for (label, value) in labels.iter().zip(&phi) {
        assert_eq!(*label, ddg_variational_cuts::region_of(*value));
    }
    let cut = ddg_variational_cuts::Cut::from_level_set(&s.mesh, &phi, 0.0).unwrap();
    assert!(!cut.points.is_empty());
    for point in &cut.points {
        assert!((point.t - 0.5).abs() < 1e-12, "{}", point.t);
    }
    // The function is a distance: one cell further its value is one cell larger.
    let h = 1.0 / 32.0;
    for (p, value) in s.positions.iter().zip(&phi) {
        let expected = p[0] - (17.0 * h - 0.5 * h);
        assert!(
            (value - expected).abs() < 0.3 * h,
            "{p:?}: {value} against {expected}"
        );
    }
    // On an irregular mesh the crossed edges are the same; the crossings are near the middle.
    let rough = square(32, 0.3);
    let labels: Vec<u32> = rough
        .positions
        .iter()
        .map(|p| u32::from(radius(*p) > 0.3))
        .collect();
    let phi = init::from_vertex_labels(&rough.mesh, &rough.geometry, &labels).unwrap();
    let cut = ddg_variational_cuts::Cut::from_level_set(&rough.mesh, &phi, 0.0).unwrap();
    let crossed = rough
        .mesh
        .edges()
        .iter()
        .filter(|e| labels[e[0] as usize] != labels[e[1] as usize])
        .count();
    assert_eq!(cut.points.len(), crossed);
    assert!(cut.points.iter().all(|p| p.t > 0.2 && p.t < 0.8));
    // One label only: no cut, and the sign of the label.
    let phi = init::from_vertex_labels(&s.mesh, &s.geometry, &vec![1; s.positions.len()]).unwrap();
    assert!(phi.iter().all(|x| *x > 0.0));

    // Face labels: a vertex takes the label that covers more of its angle; a tie gives 0.
    let face_labels: Vec<u32> = s
        .mesh
        .faces()
        .iter()
        .map(|f| {
            let x = f.iter().map(|v| s.positions[*v as usize][0]).sum::<f64>() / 3.0;
            u32::from(x > 0.5)
        })
        .collect();
    let vertex_labels = init::vertex_labels_from_faces(&s.mesh, &s.geometry, &face_labels).unwrap();
    for (p, label) in s.positions.iter().zip(&vertex_labels) {
        // The vertices on the line x = 0.5 have half of their angle on each side: label 0.
        assert_eq!(*label, u32::from(p[0] > 0.5 + 1e-9), "{p:?}");
    }
    let phi = init::from_face_labels(&s.mesh, &s.geometry, &face_labels).unwrap();
    let mut flow = Flow::new(&s.mesh, &s.geometry, &phi).unwrap();
    let (mean, std, _) = cut_statistics(&s, &mut flow, |p| p[0]);
    // The cut runs through the faces next to the interface, half a cell into label 1.
    assert!(
        (mean - (0.5 + 0.5 * h)).abs() < 1e-12 && std < 1e-12,
        "{mean} {std}"
    );
    // A single face of label 1 has no vertex of its own: it vanishes.
    let mut single = vec![0u32; s.mesh.n_faces()];
    single[100] = 1;
    let labels = init::vertex_labels_from_faces(&s.mesh, &s.geometry, &single).unwrap();
    assert!(labels.iter().all(|l| *l == 0));

    // Redistancing an arbitrary function keeps its cut and makes it a distance elsewhere.
    let wild =
        rough.sample(|p| (radius(p) - 0.3) * (1.0 + 5.0 * p[0]) * (3.0 + (9.0 * p[1]).sin()));
    let (phi, unreached) = init::redistanced(&rough.mesh, &rough.geometry, &wild).unwrap();
    assert_eq!(unreached, 0);
    let before = ddg_variational_cuts::Cut::from_level_set(&rough.mesh, &wild, 0.0).unwrap();
    let after = ddg_variational_cuts::Cut::from_level_set(&rough.mesh, &phi, 0.0).unwrap();
    assert_eq!(before.edge_point, after.edge_point);
    assert_eq!(before.segments, after.segments);
    for (p, q) in before.points.iter().zip(&after.points) {
        assert!(p.edge == q.edge && (p.t - q.t).abs() < 1e-12);
    }
    let far = rough
        .positions
        .iter()
        .zip(&phi)
        .filter(|(p, _)| radius(**p) < 0.1)
        .map(|(p, value)| value / (radius(*p) - 0.3))
        .fold((f64::INFINITY, 0.0f64), |(lo, hi), x| {
            (lo.min(x), hi.max(x))
        });
    println!("redistanced: phi / (signed distance) in the middle of the disk: {far:?}");
    assert!(wild.iter().zip(&phi).any(|(a, b)| (a - b).abs() > 0.1));
    // The slope of the input varies eightfold along the cut; the mean is taken out.
    assert!(far.0 > 0.8 && far.1 < 1.25, "{far:?}");
    // A multiple of a signed distance becomes that distance.
    let exact = rough.sample(|p| radius(p) - 0.3);
    let multiple: Vec<f64> = exact.iter().map(|x| 7.0 * x).collect();
    let (phi, _) = init::redistanced(&rough.mesh, &rough.geometry, &multiple).unwrap();
    let adjacent = ddg_variational_cuts::Cut::from_level_set(&rough.mesh, &multiple, 0.0)
        .unwrap()
        .cut_adjacent_vertices(&rough.mesh);
    for v in (0..phi.len()).filter(|v| adjacent[*v]) {
        assert!(
            (phi[v] - exact[v]).abs() < 2e-3,
            "{} against {}",
            phi[v],
            exact[v]
        );
    }
    // No cut: unchanged.
    let positive = vec![2.0; phi.len()];
    assert_eq!(
        init::redistanced(&rough.mesh, &rough.geometry, &positive)
            .unwrap()
            .0,
        positive
    );

    // The visibility penalty: the cosine of the angle to the camera, zero on the far side.
    let ball = sphere(3);
    let n = ball.positions.len();
    let camera = [0.0, 0.0, 1e6];
    let eta = init::visibility_penalty(&ball.mesh, &ball.positions, camera, &vec![2.0; n]).unwrap();
    for (p, value) in ball.positions.iter().zip(&eta) {
        let expected = 2.0 * p[2].max(0.0);
        assert!(
            (value - expected).abs() < 0.03,
            "{p:?}: {value} against {expected}"
        );
        assert!(*value >= 0.0);
    }
    let hidden: Vec<f64> = ball
        .positions
        .iter()
        .map(|p| f64::from(p[0] > 0.0))
        .collect();
    let eta = init::visibility_penalty(&ball.mesh, &ball.positions, camera, &hidden).unwrap();
    assert!(
        ball.positions
            .iter()
            .zip(&eta)
            .all(|(p, e)| p[0] > 0.0 || *e == 0.0)
    );
}

/// SPEC §7: the cut as ordered polylines, closed or ending on the boundary.
#[test]
fn cut_as_polylines() {
    // Two circles on a sphere.
    let s = sphere(3);
    let phi = s.sample(|p| (polar_angle(p) - 0.5 * PI).abs() - 0.45);
    let cut = ddg_variational_cuts::Cut::from_level_set(&s.mesh, &phi, 1e-3).unwrap();
    let lines = cut.polylines(&s.mesh, &s.positions).unwrap();
    assert_eq!(lines.len(), 2);
    let total: usize = lines.iter().map(|l| l.points.len()).sum();
    assert_eq!(total, cut.points.len());
    let at = cut.point_positions(&s.mesh, &s.positions).unwrap();
    for line in &lines {
        assert!(line.closed);
        assert_eq!(line.points.len(), line.positions.len());
        // It starts at its lowest point and goes on to the lower of the two neighbours.
        assert_eq!(line.points[0], *line.points.iter().min().unwrap());
        assert!(line.points[1] < line.points[line.points.len() - 1]);
        // One side of the equator, and consecutive points share a face.
        let side = line.positions[0][2].signum();
        let n = line.points.len();
        for k in 0..n {
            assert_eq!(line.positions[k], at[line.points[k] as usize]);
            assert_eq!(line.positions[k][2].signum(), side);
            let (p, q) = (line.points[k], line.points[(k + 1) % n]);
            assert!(
                cut.segments
                    .iter()
                    .any(|s| s.points == [p, q] || s.points == [q, p]),
                "points {p} and {q} are not joined by a segment"
            );
        }
    }
    assert!(lines[0].points[0] < lines[1].points[0]);

    // A line across a square and a circle beside it: the open piece comes first.
    let s = square(24, 0.3);
    let phi = s.sample(|p| {
        let circle = ((p[0] - 0.7).powi(2) + (p[1] - 0.5).powi(2)).sqrt() - 0.15;
        circle.min(p[0] - 0.3)
    });
    let cut = ddg_variational_cuts::Cut::from_level_set(&s.mesh, &phi, 1e-3).unwrap();
    let lines = cut.polylines(&s.mesh, &s.positions).unwrap();
    assert_eq!(lines.len(), 2);
    assert!(!lines[0].closed && lines[1].closed);
    let ends = [lines[0].points[0], *lines[0].points.last().unwrap()];
    assert!(ends[0] < ends[1]);
    for end in ends {
        assert!(
            s.mesh
                .is_boundary_edge(cut.points[end as usize].edge as usize)
        );
    }
    for line in &lines {
        let inside = line
            .points
            .iter()
            .filter(|p| {
                !s.mesh
                    .is_boundary_edge(cut.points[**p as usize].edge as usize)
            })
            .count();
        assert_eq!(inside, line.points.len() - if line.closed { 0 } else { 2 });
    }
    assert!(matches!(
        cut.polylines(&s.mesh, &s.positions[1..]),
        Err(CutError::WrongLength { .. })
    ));
    // No cut, no polylines.
    let none =
        ddg_variational_cuts::Cut::from_level_set(&s.mesh, &vec![1.0; phi.len()], 0.0).unwrap();
    assert!(none.polylines(&s.mesh, &s.positions).unwrap().is_empty());
}

/// SPEC §4 (D22): the weighted quantile that the roughening rate and the step rule use.
#[test]
fn weighted_quantile_of_the_step_rule() {
    use ddg_variational_cuts::weighted_quantile;
    let values = [3.0, 1.0, 2.0, 1.0e6];
    let weights = [1.0, 1.0, 1.0, 0.01];
    // In ascending order the weights add up to 1, 2, 3, 3.01.
    assert_eq!(weighted_quantile(&values, &weights, 1.0), 1.0e6);
    assert_eq!(weighted_quantile(&values, &weights, 0.9), 3.0);
    assert_eq!(weighted_quantile(&values, &weights, 0.5), 2.0);
    assert_eq!(weighted_quantile(&values, &weights, 0.33), 1.0);
    // Equal values keep the order of their index; a single value is every quantile.
    assert_eq!(weighted_quantile(&[5.0, 5.0], &[1.0, 3.0], 0.2), 5.0);
    assert_eq!(weighted_quantile(&[7.0], &[2.0], 0.9), 7.0);
    assert_eq!(weighted_quantile(&[], &[], 0.9), 0.0);
    // Negative values are ordered like any others.
    assert_eq!(
        weighted_quantile(&[-2.0, 4.0, -1.0], &[1.0, 1.0, 1.0], 0.6),
        -1.0
    );
}

/// SPEC §1, §4, D22: one spike of curvature on the cut does not set the roughening rate, the
/// automatic length weight or the step. A sphere cut along a latitude circle, with the density
/// of the curvature raised in one face next to the cut, against the same sphere without the
/// spike.
#[test]
fn one_spike_does_not_set_the_length_weight_or_the_step() {
    use ddg_variational_cuts::{CurvatureSource, StepStatistics, length_weight_for_ripples};
    let s = sphere(3);
    // The vertex nearest to the circle of polar angle 1.
    let spike = (0..s.positions.len())
        .min_by(|a, b| {
            let d = |v: usize| (polar_angle(s.positions[v]) - 1.0).abs();
            d(*a).total_cmp(&d(*b))
        })
        .unwrap() as u32;
    let plain = vec![4.0 * PI; s.mesh.n_faces()];
    // One face at that vertex carries 300 times the density (a fifth of the total curvature).
    let face = s
        .mesh
        .faces()
        .iter()
        .position(|face| face.contains(&spike))
        .unwrap();
    let mut spiked = plain.clone();
    spiked[face] *= 300.0;
    let phi0 = s.sample(|p| polar_angle(p) - 1.0);
    let step = StepSize::default();
    let mut rows = Vec::new();
    for (name, density) in [("sphere", plain), ("sphere with a spike", spiked)] {
        let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
        flow.redistance().unwrap();
        let mut problem = CutProblem {
            curvature: CurvatureSource::PerFace(density),
            ..CutProblem::new(Energy::Dirichlet, 0.0)
        };
        let omega = flow
            .roughening_rate(&problem, &EvalOptions::default())
            .unwrap();
        problem.length_weight = length_weight_for_ripples(omega, 0.25);
        let start = flow.level_set();
        let robust = flow.step(&problem, &FlowOptions::new(step)).unwrap();
        assert_eq!(robust.roughening_rate, omega);
        flow.set_level_set(&start).unwrap();
        let extremes = flow
            .step(
                &problem,
                &FlowOptions {
                    step_statistics: StepStatistics::Extremes,
                    ..FlowOptions::new(step)
                },
            )
            .unwrap();
        assert_eq!(extremes.n_limited, 0);
        // Ten steps with the robust statistics: the energy falls.
        flow.set_level_set(&start).unwrap();
        let report = flow
            .run(
                &problem,
                &FlowOptions::new(step),
                &StopRule {
                    max_steps: 10,
                    motion_tolerance: 0.0,
                    energy_tolerance: 0.0,
                    energy_window: 1,
                },
            )
            .unwrap();
        println!(
            "{name}: omega {:.3} (largest local rate {:.3e}), alpha_L {:.4}; tau {:.3e} with the robust statistics ({} cut points limited), {:.3e} with the extremes; energy over ten steps {:.5} -> {:.5}, time {:.3e}",
            omega,
            robust.roughening_rate_max,
            problem.length_weight,
            robust.tau,
            robust.n_limited,
            extremes.tau,
            report.steps[0].energy.total,
            report.final_energy.total,
            flow.time()
        );
        assert!(report.final_energy.total < report.steps[0].energy.total);
        rows.push((
            omega,
            robust.roughening_rate_max,
            robust.tau,
            extremes.tau,
            robust.n_limited,
        ));
    }
    let (plain, with_spike) = (rows[0], rows[1]);
    // Without a spike the two statistics see the same cut: the quantile is the maximum up to
    // the variation of the rate along a latitude circle of the icosphere.
    assert!(plain.1 / plain.0 < 1.2, "{plain:?}");
    // With it the largest local rate is far above the rate of the state.
    // Measured: 1724 against 148 (and 92 without the spike).
    assert!(with_spike.1 > 8.0 * with_spike.0, "{with_spike:?}");
    // The step: the extremes shrink it by the spike, the robust statistics do not, and the cut
    // points at the spike take a smaller step of their own.
    // Measured: 2.9e-5 with the extremes, 2.8e-4 with the robust statistics (3.2e-4 without the
    // spike), 14 of 86 cut points limited.
    assert!(with_spike.3 < 0.2 * with_spike.2, "{rows:?}");
    assert!(with_spike.2 > 0.5 * plain.2, "{rows:?}");
    assert!(with_spike.4 >= 1 && with_spike.4 <= 20, "{rows:?}");
    // The rate of the state, and with it the automatic length weight, moves by less than a
    // factor of two (the spike adds a fifth to the total curvature).
    assert!(
        with_spike.0 < 2.0 * plain.0 && with_spike.0 > plain.0,
        "{rows:?}"
    );
}

/// SPEC §3.9, D30: a face with a moved crossing does not count in the mean slope. A circle on
/// the flat grid that passes a vertex at the distance `ε`: the crossings on the edges at that
/// vertex are moved away from it by `δ`, and the slope formula of such a face gives
/// `ε / (δ · h)`. Counted, it made the harmonic mean of the whole cut `0.5` for `ε = 10⁻⁸` and
/// `0.01` for `10⁻¹⁰`, and the first step moved the cut by that much less. `ε = 0` was the one
/// case that was left out before.
#[test]
fn a_vertex_almost_on_the_cut_does_not_scale_the_rate_of_the_whole_cut() {
    let g = grid_disk(32, 32).unwrap();
    let mesh = ddg_core::TriMesh::new(g.positions.len(), g.faces.clone()).unwrap();
    let geometry = ddg_core::Geometry::from_positions(&mesh, &g.positions).unwrap();
    let problem = CutProblem::new(Energy::Hencky, 1.0);
    let options = FlowOptions::new(StepSize::default());
    // The vertex (22, 16)/32 is at the distance 6/32 from the middle of the square.
    let first_step = |eps: f64| {
        let phi0: Vec<f64> = g
            .positions
            .iter()
            .map(|p| (p[0] - 0.5).hypot(p[1] - 0.5) - (6.0 / 32.0 + eps))
            .collect();
        let mut flow = Flow::new(&mesh, &geometry, &phi0).unwrap();
        let moved = flow
            .evaluate(&problem, &options.eval)
            .unwrap()
            .cut
            .moved
            .iter()
            .filter(|m| **m)
            .count();
        let (report, stages) = flow.step_with_stages(&problem, &options).unwrap();
        assert_eq!(report.n_clamped_crossings, moved);
        let scale = stages
            .gradient_scale
            .iter()
            .fold(f64::INFINITY, |m, g| m.min(*g));
        (report, scale, moved)
    };
    let (clean, scale, moved) = first_step(1e-2);
    assert_eq!(moved, 0);
    assert!((scale - 1.0).abs() < 0.02);
    for eps in [1e-6, 1e-8, 1e-10, 1e-13, 0.0, -1e-10] {
        let (report, scale, moved) = first_step(eps);
        assert!(moved >= 8, "{eps:e}: {moved} crossings moved");
        assert!((scale - 1.0).abs() < 0.02, "{eps:e}: mean slope {scale}");
        let slopes = report.gradient_norm.unwrap();
        assert!(slopes.min > 0.9 && slopes.max < 1.1, "{eps:e}: {slopes:?}");
        assert!(
            (report.cut_motion / clean.cut_motion - 1.0).abs() < 0.2,
            "{eps:e}: the first step moves the cut by {:e}",
            report.cut_motion
        );
        assert!((report.predicted_motion / clean.predicted_motion - 1.0).abs() < 0.2);
    }
}
