//! SPEC §6 item 3: finite differences of the discrete energy against the change predicted by
//! the speed on the cut. This decides every sign of SPEC §1 (in particular D2: which region
//! enters the speed with a plus).
//!
//! For a perturbation `φ + εψ` the cut point `c` sweeps the area `dA_c` into region 0, and the
//! shape derivative of SPEC §1 predicts `dE = Σ_c (G_0 − G_1)_c dA_c = Σ_c s_c dA_c`. The left
//! side is measured by a central difference of the discrete energy, the right side uses the
//! speed `s` of the evaluation and the exact swept areas of the cut mesh. The two agree up to
//! the discretization error of `s`, which must shrink under refinement; a wrong sign would give
//! a relative error of 2.

mod common;

use common::{Gradient, Shape, bumpy, check_gradient, shift_away_from_vertices, sphere};
use ddg_variational_cuts::{
    CutProblem, Energy, EvalOptions, PAPER_FEATURE_SMOOTHING, PointwiseEnergy, RegionWeights,
    ScaleWeight,
};

/// A wobbly, non-symmetric closed curve around the north pole of the unit sphere.
fn sphere_level_set(p: [f64; 3]) -> f64 {
    let [x, y, z] = p;
    0.35 - (0.9 * z - 0.3 * x - 0.2 * y - 0.25 * x * y - 0.15 * (2.5 * x + 0.5).sin() * y)
}

/// A curve from the bottom edge to the top edge of the unit square.
fn disk_level_set(p: [f64; 3]) -> f64 {
    let [x, y, _] = p;
    x - 0.46 - 0.13 * (2.7 * y + 0.4).sin() - 0.06 * y * y
}

/// A positive, non-constant weight for the scaled Dirichlet energy.
fn scale_weight(p: [f64; 3]) -> f64 {
    0.6 + 0.4 * p[0] * p[0] + 0.3 * (2.0 * p[1] + 1.0).sin().powi(2) + 0.2 * p[2] * p[2]
}

const DIRECTIONS: [&str; 2] = ["uniform", "wave"];

fn direction(name: &str, p: [f64; 3]) -> f64 {
    match name {
        "uniform" => 1.0,
        _ => 0.3 + (3.3 * p[0] + 0.7).sin() * (2.1 * p[1] - 0.4).cos() + 0.4 * p[2],
    }
}

fn problems(shape: &Shape) -> Vec<(&'static str, CutProblem)> {
    let weighted = vec![
        RegionWeights {
            distortion: 0.7,
            area: 0.3,
        },
        RegionWeights {
            distortion: 1.9,
            area: -0.2,
        },
    ];
    let only = |r: usize| -> Vec<RegionWeights> {
        (0..2)
            .map(|k| RegionWeights {
                distortion: if k == r { 1.0 } else { 0.0 },
                area: 0.0,
            })
            .collect()
    };
    let scaled = Energy::ScaledDirichlet(ScaleWeight::PerVertex(shape.sample(scale_weight)));
    vec![
        ("Dirichlet", CutProblem::new(Energy::Dirichlet, 0.0)),
        (
            "Dirichlet, region 0 only",
            CutProblem {
                regions: only(0),
                ..CutProblem::new(Energy::Dirichlet, 0.0)
            },
        ),
        (
            "Dirichlet, region 1 only",
            CutProblem {
                regions: only(1),
                ..CutProblem::new(Energy::Dirichlet, 0.0)
            },
        ),
        ("Hencky", CutProblem::new(Energy::Hencky, 0.0)),
        (
            "scaled Dirichlet, a(x)",
            CutProblem::new(scaled.clone(), 0.0),
        ),
        (
            "scaled Dirichlet, feature size",
            CutProblem::new(
                Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
                    r_eps: 0.3,
                    smoothing: PAPER_FEATURE_SMOOTHING,
                }),
                0.0,
            ),
        ),
        (
            "Dirichlet, w = (0.7, 1.9), beta = (0.3, -0.2)",
            CutProblem {
                regions: weighted.clone(),
                ..CutProblem::new(Energy::Dirichlet, 0.0)
            },
        ),
        (
            "Hencky, w = (0.7, 1.9), beta = (0.3, -0.2)",
            CutProblem {
                regions: weighted.clone(),
                ..CutProblem::new(Energy::Hencky, 0.0)
            },
        ),
        (
            "scaled Dirichlet, a(x), w, beta",
            CutProblem {
                regions: weighted.clone(),
                ..CutProblem::new(scaled, 0.0)
            },
        ),
        (
            // b(0) and b'(0) are not zero: the constant term of G and the boundary term of the
            // adjoint flux q are exercised.
            "pointwise b = 0.5 + 0.3u + u^2 + u^4, w, beta",
            CutProblem {
                regions: weighted,
                ..CutProblem::new(
                    Energy::Pointwise(PointwiseEnergy::new(
                        |u| 0.5 + 0.3 * u + u * u + u.powi(4),
                        |u| 0.3 + 2.0 * u + 4.0 * u.powi(3),
                    )),
                    0.0,
                )
            },
        ),
    ]
}

/// The level-set function of the shape, shifted so that no crossing is close to a vertex.
fn level_set(shape: &Shape, f: fn([f64; 3]) -> f64) -> Vec<f64> {
    let h = shape.geometry.edge_lengths().iter().sum::<f64>() / shape.mesh.n_edges() as f64;
    let shifts: Vec<f64> = (-30..=30).map(|j| 0.01 * h * f64::from(j)).collect();
    let (phi, offset) = shift_away_from_vertices(shape, &shape.sample(f), &shifts);
    assert!(offset > 1e-4, "{}: smallest offset {offset}", shape.name);
    phi
}

/// One check: `(problem, direction, with the default options, with the options of the paper)`.
type Row = (&'static str, &'static str, Gradient, Option<Gradient>);

/// How much of the matrix of problems and directions a run covers.
#[derive(Clone, Copy, PartialEq)]
enum Coverage {
    /// Every problem, both directions, default options and the options of the paper.
    Tables,
    /// Every problem, both directions, default options.
    Full,
    /// The first five problems (each energy, each region alone), the non-uniform direction.
    Short,
}

/// All checks on one shape.
fn run(shape: &Shape, f: fn([f64; 3]) -> f64, coverage: Coverage) -> Vec<Row> {
    let phi = level_set(shape, f);
    let mut evaluator = shape.evaluator();
    // Nothing may be clamped: a clamped crossing does not follow the level-set function.
    let default = EvalOptions {
        min_crossing_offset: 0.0,
        ..EvalOptions::default()
    };
    let paper = EvalOptions {
        min_crossing_offset: 0.0,
        ..EvalOptions::paper()
    };
    let mut out = Vec::new();
    let short = coverage == Coverage::Short;
    for (name, problem) in problems(shape).into_iter().take(if short { 5 } else { 99 }) {
        for dir in DIRECTIONS.into_iter().skip(usize::from(short)) {
            let psi = shape.sample(|p| direction(dir, p));
            let a = check_gradient(shape, &mut evaluator, &phi, &psi, &problem, &default);
            let b = (coverage == Coverage::Tables)
                .then(|| check_gradient(shape, &mut evaluator, &phi, &psi, &problem, &paper));
            out.push((name, dir, a, b));
        }
    }
    out
}

fn print(shape: &Shape, rows: &[Row]) {
    println!();
    println!(
        "{} ({} vertices, {} cut points)",
        shape.name,
        shape.positions.len(),
        rows[0].2.n_cut_points
    );
    println!(
        "{:<46} {:<8} {:>13} {:>13} {:>9} {:>9} | {:>13} {:>9}",
        "energy",
        "psi",
        "finite diff",
        "predicted",
        "err/scale",
        "rel err",
        "paper pred.",
        "err/scale"
    );
    for (name, dir, a, b) in rows {
        let b = b.expect("the tables include the options of the paper");
        println!(
            "{:<46} {:<8} {:>13.6e} {:>13.6e} {:>9.2e} {:>9.2e} | {:>13.6e} {:>9.2e}",
            name,
            dir,
            a.finite_difference,
            a.predicted,
            a.error(),
            a.relative_error(),
            b.predicted,
            b.error()
        );
    }
}

/// The largest error over all problems and directions, with the default options.
fn worst(rows: &[Row]) -> f64 {
    rows.iter().map(|r| r.2.error()).fold(0.0, f64::max)
}

/// Checks every problem and direction on one shape against `bound` and returns the largest
/// error.
fn check(shape: &Shape, f: fn([f64; 3]) -> f64, bound: f64, coverage: Coverage) -> f64 {
    let rows = run(shape, f, coverage);
    for (name, dir, a, _) in &rows {
        // A wrong sign of a term gives an error of the order of 1 (one region) or 2.
        assert!(
            a.error() < bound,
            "{} {name} {dir}: finite difference {} against predicted {}, error {}",
            shape.name,
            a.finite_difference,
            a.predicted,
            a.error()
        );
    }
    // With a single weighted region the sign of the change is the sign of D2: shifting phi up
    // shrinks region 0, which lowers its own energy and raises that of region 1.
    for (name, dir, a, _) in &rows {
        if *dir == "uniform" && name.ends_with("region 0 only") {
            assert!(a.finite_difference < 0.0 && a.predicted < 0.0);
        }
        if *dir == "uniform" && name.ends_with("region 1 only") {
            assert!(a.finite_difference > 0.0 && a.predicted > 0.0);
        }
    }
    assert!(coverage != Coverage::Tables);
    worst(&rows)
}

#[test]
fn sphere_with_a_non_symmetric_cut() {
    // Measured on levels 2 to 6: 5.5e-2, 1.8e-2, 4.4e-3, 1.1e-3, 2.7e-4: second order.
    let bounds = [0.09, 0.03, 0.008];
    let errors: Vec<f64> = [2, 3, 4]
        .iter()
        .zip(bounds)
        .map(|(&level, bound)| {
            // The finest level takes most of the time; it runs the short matrix.
            let coverage = if level == 4 {
                Coverage::Short
            } else {
                Coverage::Full
            };
            check(&sphere(level), sphere_level_set, bound, coverage)
        })
        .collect();
    assert!(
        errors[1] < 0.5 * errors[0] && errors[2] < 0.5 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn bumpy_disk_with_a_cut_that_meets_the_boundary() {
    // Measured: 2.0e-1, 5.0e-2, 1.7e-2 (and 3.9e-3 on the 128 grid).
    let bounds = [0.3, 0.08, 0.03];
    let errors: Vec<f64> = [16, 32, 64]
        .iter()
        .zip(bounds)
        .map(|(&n, bound)| {
            let shape = bumpy(n, 0.12);
            // The cut ends on the boundary of the mesh.
            let phi = level_set(&shape, disk_level_set);
            let cut = ddg_variational_cuts::Cut::from_level_set(&shape.mesh, &phi, 0.0).unwrap();
            let on_boundary = cut
                .points
                .iter()
                .filter(|p| shape.mesh.is_boundary_edge(p.edge as usize))
                .count();
            assert_eq!(on_boundary, 2);
            let coverage = if n == 64 {
                Coverage::Short
            } else {
                Coverage::Full
            };
            check(&shape, disk_level_set, bound, coverage)
        })
        .collect();
    assert!(
        errors[1] < 0.5 * errors[0] && errors[2] < 0.5 * errors[1],
        "{errors:?}"
    );
}

/// Prints the tables of the report. Run with
/// `cargo test --release -p ddg-variational-cuts --test finite_differences -- --ignored --nocapture`.
#[test]
#[ignore = "prints the finite-difference tables (levels up to 6, grids up to 128)"]
fn print_finite_difference_tables() {
    for level in [2, 3, 4, 5, 6] {
        let shape = sphere(level);
        let rows = run(&shape, sphere_level_set, Coverage::Tables);
        print(&shape, &rows);
    }
    for n in [16, 32, 64, 128] {
        let shape = bumpy(n, 0.12);
        let rows = run(&shape, disk_level_set, Coverage::Tables);
        print(&shape, &rows);
    }
}
