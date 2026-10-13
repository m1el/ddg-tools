//! SPEC §6 item 2 and the checks of SPEC §1: analytic solutions.
//!
//! * A cap of the unit-area sphere of polar angle `ψ₀`: `u = 2 ln(cos(ψ₀/2)/cos(ψ/2))`,
//!   `∂u/∂n = tan(ψ₀/2)/R`, `J = −8π ln cos(ψ₀/2) − 4π sin²(ψ₀/2)`; its complement is the cap of
//!   polar angle `π − ψ₀` about the other pole.
//! * A flat disk of radius `ρ` with the constant source `K`: `u = K (r² − ρ²)/4`,
//!   `∫|∇u|² = π K² ρ⁴/8`, `∫u² = π K² ρ⁶/48`, `∂u/∂n = K ρ/2`, `∂p/∂n = −K ρ³/8` (Hencky).

mod common;

use std::f64::consts::PI;

use common::{Shape, polar_angle, sphere, unit_area_radius};
use ddg_core::generators::grid_disk;
use ddg_variational_cuts::{
    CurvatureSource, CutProblem, Energy, EvalOptions, Evaluation, FluxRecovery, RegionWeights,
    ScaleWeight,
};

const PSI0: f64 = 1.0;

/// Largest and mean of `|value/exact − 1|`.
fn relative_errors(values: impl Iterator<Item = f64>, exact: f64) -> (f64, f64) {
    let (mut max, mut sum, mut n) = (0.0f64, 0.0, 0usize);
    for v in values {
        let e = (v / exact - 1.0).abs();
        max = max.max(e);
        sum += e;
        n += 1;
    }
    (max, sum / n as f64)
}

/// The errors of one evaluation of the cap on one icosphere.
struct CapRow {
    level: u32,
    n_vertices: usize,
    n_cut: usize,
    n_sliver: usize,
    smallest_offset: f64,
    /// `max |u − exact|` over the vertices of the input mesh.
    u: f64,
    /// Relative errors of `J_0`, `J_1`, the cut length and the area of region 0.
    j: [f64; 2],
    length: f64,
    area: f64,
    /// `(max, mean)` relative error of `∂u/∂n` on the side of region 0 and of region 1, for the
    /// nodal recovery, the nodal recovery with one round of sliver smoothing (square root of
    /// the smoothed gradient density), and the smoothed recovery.
    nodal: [(f64, f64); 2],
    nodal_smoothed: [(f64, f64); 2],
    recovered: [(f64, f64); 2],
    /// `(max, mean)` relative error of the speed `s` for the same three.
    speed: [(f64, f64); 3],
}

fn cap_problem(weights: [f64; 2]) -> CutProblem {
    CutProblem {
        regions: weights
            .iter()
            .map(|&w| RegionWeights {
                distortion: w,
                area: 0.0,
            })
            .collect(),
        ..CutProblem::new(Energy::Dirichlet, 0.0)
    }
}

fn cap_row(level: u32) -> CapRow {
    let shape = sphere(level);
    let mut evaluator = shape.evaluator();
    let phi = shape.sample(|p| polar_angle(p) - PSI0);
    let r = unit_area_radius();
    let exact_d = [(0.5 * PSI0).tan() / r, 1.0 / ((0.5 * PSI0).tan() * r)];
    let exact_s = exact_d[0] * exact_d[0] - exact_d[1] * exact_d[1];
    let both = cap_problem([1.0, 1.0]);

    let nodal_options = EvalOptions {
        flux_recovery: FluxRecovery::Nodal,
        ..EvalOptions::default()
    };
    let e = evaluator
        .evaluate(&phi, &both, &EvalOptions::default())
        .unwrap();
    let nodal = evaluator.evaluate(&phi, &both, &nodal_options).unwrap();
    let paper = evaluator
        .evaluate(&phi, &both, &EvalOptions::paper())
        .unwrap();
    assert!(e.report.yamabe_residual < 1e-9);

    // u against the exact solution of each cap.
    let mut u_error = 0.0f64;
    for (v, &p) in shape.positions.iter().enumerate() {
        let psi = polar_angle(p);
        let exact = if e.cut.vertex_region[v] == 0 {
            2.0 * ((0.5 * PSI0).cos() / (0.5 * psi).cos()).ln()
        } else {
            2.0 * ((0.5 * PSI0).sin() / (0.5 * psi).sin()).ln()
        };
        u_error = u_error.max((e.log_scale[v] - exact).abs());
    }
    let exact_j = [
        -8.0 * PI * (0.5 * PSI0).cos().ln() - 4.0 * PI * (0.5 * PSI0).sin().powi(2),
        -8.0 * PI * (0.5 * PSI0).sin().ln() - 4.0 * PI * (0.5 * PSI0).cos().powi(2),
    ];
    let side = |e: &Evaluation, r: usize| -> Vec<f64> {
        (0..e.cut.points.len())
            .map(|k| e.normal_derivative[e.triangles.point_vertices_by_region(k)[r]])
            .collect()
    };
    let per_side = |e: &Evaluation| -> [(f64, f64); 2] {
        [0, 1].map(|r| relative_errors(side(e, r).into_iter(), exact_d[r]))
    };
    // The paper smooths the speed, not the normal derivative. With only one region weighted
    // the speed is ±(∂u/∂n)² of that region, so its square root shows what smoothing does.
    let nodal_smoothed = [0, 1].map(|r| {
        let mut weights = [0.0, 0.0];
        weights[r] = 1.0;
        let one = evaluator
            .evaluate(&phi, &cap_problem(weights), &EvalOptions::paper())
            .unwrap();
        let sign = if r == 0 { 1.0 } else { -1.0 };
        relative_errors(
            one.speed.iter().map(|s| (sign * s).abs().sqrt()),
            exact_d[r],
        )
    });
    let speed = |e: &Evaluation| relative_errors(e.speed.iter().copied(), exact_s);
    CapRow {
        level,
        n_vertices: shape.positions.len(),
        n_cut: e.cut.points.len(),
        n_sliver: e.report.n_sliver_points,
        smallest_offset: e.report.smallest_crossing_offset,
        u: u_error,
        j: [0, 1].map(|r| (e.energy.distortion[r] / exact_j[r] - 1.0).abs()),
        length: (e.energy.length / (2.0 * PI * r * PSI0.sin()) - 1.0).abs(),
        area: (e.energy.area[0] / (0.5 * (1.0 - PSI0.cos())) - 1.0).abs(),
        nodal: per_side(&nodal),
        nodal_smoothed,
        recovered: per_side(&e),
        speed: [speed(&nodal), speed(&paper), speed(&e)],
    }
}

#[test]
fn sphere_cap_converges() {
    let rows: Vec<CapRow> = (3..=5).map(cap_row).collect();
    for pair in rows.windows(2) {
        let (coarse, fine) = (&pair[0], &pair[1]);
        // u, J_0, the length and the area are second order.
        assert!(fine.u < 0.35 * coarse.u, "u: {} -> {}", coarse.u, fine.u);
        assert!(fine.j[0] < 0.35 * coarse.j[0], "J_0");
        assert!(fine.length < 0.35 * coarse.length && fine.area < 0.35 * coarse.area);
    }
    // Measured (levels 3, 4, 5; level 6 is in the printed table):
    //   max |u − u*|            2.7e-3  6.7e-4  1.7e-4
    //   J_0                     5.5e-3  1.4e-3  3.8e-4      J_1  7.7e-5  1.2e-5  7.0e-6
    //   du/dn max, region 0     1.7e-2  5.9e-3  3.8e-3      region 1  7.3e-3  4.8e-3  4.1e-3
    //   s max                   1.7e-2  1.0e-2  9.5e-3
    let (first, last) = (&rows[0], &rows[2]);
    assert!(last.u < 3e-4, "{}", last.u);
    assert!(last.j[0] < 8e-4 && last.j[1] < 1e-4, "{:?}", last.j);
    for r in 0..2 {
        // The recovered normal derivative converges in the maximum norm and in the mean.
        assert!(last.recovered[r].0 < 0.01, "{:?}", last.recovered);
        assert!(last.recovered[r].1 < 3e-3, "{:?}", last.recovered);
        assert!(last.recovered[r].0 < 0.7 * first.recovered[r].0);
        assert!(last.recovered[r].1 < 0.5 * first.recovered[r].1);
    }
    assert!(last.speed[2].0 < 0.02, "{:?}", last.speed);
    assert!(last.speed[2].0 < 0.7 * first.speed[2].0 && last.speed[2].1 < 0.5 * first.speed[2].1);
}

/// Prints the convergence table of the report. Run with
/// `cargo test --release -p ddg-variational-cuts --test analytic -- --ignored --nocapture`.
#[test]
#[ignore = "prints the sphere-cap convergence table (levels 3 to 6)"]
fn print_sphere_cap_table() {
    let rows: Vec<CapRow> = (3..=6).map(cap_row).collect();
    println!("sphere cap, polar angle {PSI0}; relative errors");
    println!(
        "{:>5} {:>8} {:>6} {:>7} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}",
        "level",
        "vertices",
        "cut",
        "slivers",
        "min off",
        "max|u-u*|",
        "J_0",
        "J_1",
        "length",
        "area_0"
    );
    for r in &rows {
        println!(
            "{:>5} {:>8} {:>6} {:>7} {:>9.1e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e}",
            r.level,
            r.n_vertices,
            r.n_cut,
            r.n_sliver,
            r.smallest_offset,
            r.u,
            r.j[0],
            r.j[1],
            r.length,
            r.area
        );
    }
    println!();
    println!("du/dn at the cut vertices, relative error max / mean");
    println!(
        "{:>5} {:>6} | {:>19} {:>19} {:>19}",
        "level", "region", "nodal v/l", "nodal + 1 sliver pass", "smoothed (default)"
    );
    let f = |(max, mean): (f64, f64)| format!("{max:>9.2e} /{mean:>9.2e}");
    for r in &rows {
        for region in 0..2 {
            println!(
                "{:>5} {:>6} | {} {} {}",
                r.level,
                region,
                f(r.nodal[region]),
                f(r.nodal_smoothed[region]),
                f(r.recovered[region])
            );
        }
    }
    println!();
    println!("speed s = G_0 - G_1, relative error max / mean");
    println!(
        "{:>5} | {:>19} {:>19} {:>19}",
        "level", "nodal v/l", "nodal + 1 sliver pass", "smoothed (default)"
    );
    for r in &rows {
        println!(
            "{:>5} | {} {} {}",
            r.level,
            f(r.speed[0]),
            f(r.speed[1]),
            f(r.speed[2])
        );
    }
}

/// The unit square, `n × n` cells, flat.
fn square(n: usize) -> Shape {
    Shape::from_generated(&format!("grid {n}"), grid_disk(n, n).unwrap())
}

#[test]
fn flat_disk_with_a_uniform_source() {
    // SPEC §1, "Checks". The unit square has area 1, so there is no change of unit.
    let (k, rho) = (5.0, 0.3);
    let center = [0.513, 0.479];
    let shape = square(64);
    let mut evaluator = shape.evaluator();
    let radius = |p: [f64; 3]| ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2)).sqrt();
    let phi = shape.sample(|p| radius(p) - rho);
    let options = EvalOptions::default();
    let uniform = |energy: Energy| CutProblem {
        curvature: CurvatureSource::Uniform(k),
        ..CutProblem::new(energy, 0.0)
    };

    let d = evaluator
        .evaluate(&phi, &uniform(Energy::Dirichlet), &options)
        .unwrap();
    assert!(d.adjoint.is_none() && d.report.adjoint_residual.is_none());
    let mut worst = 0.0f64;
    for (v, &p) in shape.positions.iter().enumerate() {
        if d.cut.vertex_region[v] == 0 {
            let exact = 0.25 * k * (radius(p).powi(2) - rho * rho);
            worst = worst.max((d.log_scale[v] - exact).abs());
        }
    }
    assert!(worst < 2e-4, "max |u - exact| = {worst}");
    let exact_j = PI * k * k * rho.powi(4) / 8.0;
    assert!((d.energy.distortion[0] / exact_j - 1.0).abs() < 2e-3);
    let inner = |e: &Evaluation, x: &[f64]| -> Vec<f64> {
        (0..e.cut.points.len())
            .map(|c| x[e.triangles.point_vertices_by_region(c)[0]])
            .collect()
    };
    let (max, mean) = relative_errors(inner(&d, &d.normal_derivative).into_iter(), 0.5 * k * rho);
    assert!(max < 0.02 && mean < 0.004, "du/dn: max {max}, mean {mean}");
    // G_0 = (K rho / 2)^2 > 0: the disk shrinks if only it is weighted.
    let (max, _) = relative_errors(inner(&d, &d.gradient).into_iter(), (0.5 * k * rho).powi(2));
    assert!(max < 0.04);

    // Hencky: the adjoint and its normal derivative.
    let h = evaluator
        .evaluate(&phi, &uniform(Energy::Hencky), &options)
        .unwrap();
    assert_eq!(h.log_scale, d.log_scale);
    assert!(h.report.adjoint_residual.unwrap() < 1e-9);
    let exact_j = PI * k * k * rho.powi(6) / 48.0;
    assert!((h.energy.distortion[0] / exact_j - 1.0).abs() < 1e-2);
    let dp = h.adjoint_normal_derivative.as_ref().unwrap();
    let (max, mean) = relative_errors(inner(&h, dp).into_iter(), -k * rho.powi(3) / 8.0);
    assert!(max < 0.03 && mean < 0.006, "dp/dn: max {max}, mean {mean}");
    let (max, _) = relative_errors(
        inner(&h, &h.gradient).into_iter(),
        k * k * rho.powi(4) / 16.0,
    );
    assert!(max < 0.05, "G (Hencky): max {max}");

    // Scaled Dirichlet with a = 1 is the Dirichlet energy: p = −2u, q = −2v.
    let ones = Energy::ScaledDirichlet(ScaleWeight::PerVertex(vec![1.0; shape.positions.len()]));
    let s = evaluator.evaluate(&phi, &uniform(ones), &options).unwrap();
    let (p, q) = (
        s.adjoint.as_ref().unwrap(),
        s.adjoint_flux.as_ref().unwrap(),
    );
    for v in 0..s.triangles.n_vertices() {
        assert!((p[v] + 2.0 * d.log_scale[v]).abs() < 1e-11);
        assert!((q[v] + 2.0 * d.flux[v]).abs() < 1e-11);
    }
    for c in 0..d.speed.len() {
        assert!((s.speed[c] - d.speed[c]).abs() < 1e-9 * d.speed[c].abs().max(1.0));
    }
    assert!((s.energy.total - d.energy.total).abs() < 1e-12);
}

#[test]
fn flat_surface_has_no_distortion() {
    // SPEC §6 item 4, the static part: on a flat disk u = 0 and s = 0 for every cut.
    let shape = square(24);
    let mut evaluator = shape.evaluator();
    let phi = shape.sample(|p| ((p[0] - 0.45).powi(2) + (p[1] - 0.55).powi(2)).sqrt() - 0.31);
    for energy in [Energy::Dirichlet, Energy::Hencky] {
        let e = evaluator
            .evaluate(&phi, &CutProblem::new(energy, 1.0), &EvalOptions::default())
            .unwrap();
        assert!(e.log_scale.iter().all(|u| u.abs() < 1e-12));
        assert!(e.speed.iter().all(|s| s.abs() < 1e-18));
        assert!(e.energy.distortion.iter().all(|j| j.abs() < 1e-20));
        // Only the length is left: a circle of radius 0.31, as a polygon.
        let circle = 2.0 * PI * 0.31;
        assert!((e.energy.total / circle - 1.0).abs() < 2e-3);
        assert_eq!(e.energy.total, e.energy.length);
    }
}

#[test]
fn length_and_penalty_terms_on_the_sphere() {
    // E_L = ∫ (α_L + η) ds on the circle of polar angle ψ₀ with η = η(z): the kernel's value,
    // and its derivative when φ is shifted, against SPEC §1: dE = ∫ [∂w/∂n + w κ] σ ds.
    let shape = sphere(5);
    let mut evaluator = shape.evaluator();
    let r = unit_area_radius();
    let alpha = 0.7;
    let eta = |z: f64| 0.3 + 0.2 * z + 0.1 * z * z;
    let problem = CutProblem {
        penalty: Some(shape.sample(|p| eta(p[2]))),
        // No distortion: only the length terms are left in the total.
        regions: vec![
            RegionWeights {
                distortion: 0.0,
                area: 0.0
            };
            2
        ],
        ..CutProblem::new(Energy::Dirichlet, alpha)
    };
    let total = |evaluator: &mut ddg_variational_cuts::Evaluator, psi0: f64| -> f64 {
        let phi = shape.sample(|p| polar_angle(p) - psi0);
        evaluator
            .evaluate(&phi, &problem, &EvalOptions::default())
            .unwrap()
            .energy
            .total
    };
    let exact = |psi0: f64| (alpha + eta(psi0.cos())) * 2.0 * PI * r * psi0.sin();
    let e = total(&mut evaluator, PSI0);
    assert!(
        (e / exact(PSI0) - 1.0).abs() < 1e-3,
        "{e} vs {}",
        exact(PSI0)
    );

    // φ is the distance on the unit sphere of the caller, where the circle has the geodesic
    // curvature κ = cot ψ₀ and ∂η/∂n = −η'(z) sin ψ₀. Raising φ by ε moves the cut by
    // σ = −ε along n, and lengths on the unit-area sphere are r times those of the caller:
    // dE/dε = −r ∮ [∂η/∂n + (α + η) κ] ds_caller.
    let z = PSI0.cos();
    let d_eta = 0.2 + 0.2 * z;
    let circumference = 2.0 * PI * PSI0.sin();
    let predicted =
        -r * circumference * (-d_eta * PSI0.sin() + (alpha + eta(z)) * PSI0.cos() / PSI0.sin());
    let step = 0.02;
    let measured =
        (total(&mut evaluator, PSI0 - step) - total(&mut evaluator, PSI0 + step)) / (2.0 * step);
    assert!(
        (measured / predicted - 1.0).abs() < 5e-3,
        "measured {measured}, predicted {predicted}"
    );
}
