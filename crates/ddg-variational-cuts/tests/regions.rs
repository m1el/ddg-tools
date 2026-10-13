//! `N` regions (SPEC §8): the tests of §8.6.
//!
//! The tolerances are a small multiple of what the same runs measure; the tables are in the
//! README.

mod common;

use std::f64::consts::PI;

use common::{Shape, bumpy, polar_angle, sphere};
use ddg_core::generators::SplitMix64;
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CurveEnd, Cut, CutError, CutProblem, Energy, EvalOptions, Flow, FlowOptions,
    PAPER_FEATURE_SMOOTHING, PointwiseEnergy, RegionFlow, RegionWeights, ScaleWeight, StepSize,
    StopRule, init,
};

fn fixed(tau: f64) -> FlowOptions<'static> {
    FlowOptions::new(StepSize::Fixed(tau))
}

fn paper(tau: f64) -> FlowOptions<'static> {
    FlowOptions::paper(StepSize::Fixed(tau))
}

/// The stages of the paper with the flux recovery of SPEC §3.5 (what the oracle implements).
fn paper_stages(tau: f64) -> FlowOptions<'static> {
    let mut options = FlowOptions::paper(StepSize::Fixed(tau));
    options.eval = EvalOptions::default();
    options
}

fn largest_difference(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// SPEC §8.6, §6 item 7: the state `[φ, −φ]` takes the steps of the two-region flow, for the
/// default stages, the stages of the paper, and the paper as printed, on a closed mesh and on
/// one with boundary, with a penalty, weights and symmetry classes in play.
#[test]
fn two_regions_reproduce_the_two_region_flow() {
    let ball = sphere(3);
    let disk = bumpy(20, 0.08);
    let ball_phi = ball.sample(|p| polar_angle(p) - 1.0 - 0.1 * (3.0 * p[1].atan2(p[0])).cos());
    let disk_phi = disk.sample(|p| ((p[0] - 0.45).powi(2) + (p[1] - 0.55).powi(2)).sqrt() - 0.27);
    for (shape, phi) in [(&ball, &ball_phi), (&disk, &disk_phi)] {
        let n = shape.positions.len();
        let problems = [
            CutProblem::new(Energy::Dirichlet, 2.0),
            CutProblem {
                penalty: Some(shape.sample(|p| 1.5 + p[0])),
                regions: vec![
                    RegionWeights {
                        distortion: 1.5,
                        area: -0.5,
                    },
                    RegionWeights {
                        distortion: 0.5,
                        area: 0.25,
                    },
                ],
                ..CutProblem::new(Energy::Hencky, 0.3)
            },
            CutProblem::new(
                Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
                    r_eps: 0.3,
                    smoothing: PAPER_FEATURE_SMOOTHING,
                }),
                0.01,
            ),
        ];
        let mut with_classes = fixed(2e-4);
        with_classes.symmetry = Some((0..n as u32).map(|v| v / 2).collect());
        let stage_sets = [
            ("default stages", fixed(2e-4)),
            (
                "default stages, adaptive step",
                FlowOptions::new(StepSize::default()),
            ),
            ("default stages, symmetry classes", with_classes),
            ("stages of the paper", paper_stages(2e-4)),
            ("the paper as printed", paper(2e-5)),
        ];
        let mut worst = 0.0f64;
        for problem in &problems {
            for (name, options) in &stage_sets {
                let mut two = Flow::new(&shape.mesh, &shape.geometry, phi).unwrap();
                let minus: Vec<f64> = phi.iter().map(|x| -x).collect();
                let mut many =
                    RegionFlow::new(&shape.mesh, &shape.geometry, &[phi.clone(), minus]).unwrap();
                assert_eq!(many.n_regions(), 2);
                for step in 0..4 {
                    let a = two.step(problem, options).unwrap();
                    let b = many.step(problem, options).unwrap();
                    let size = two.level_set().iter().fold(0.0f64, |m, x| m.max(x.abs()));
                    let sets = many.level_sets();
                    let d0 = largest_difference(&sets[0], &two.level_set()) / size;
                    let opposite: Vec<f64> = sets[1].iter().map(|x| -x).collect();
                    let d1 = largest_difference(&opposite, &two.level_set()) / size;
                    worst = worst.max(d0).max(d1);
                    assert!(
                        d0 < 1e-12 && d1 < 1e-12,
                        "{}, {name}, step {step}: {d0:e} {d1:e}",
                        shape.name
                    );
                    // The reports agree too.
                    assert_eq!(a.n_cut_points, b.n_cut_points);
                    assert_eq!(a.n_unreached, b.n_unreached);
                    assert_eq!(a.cut_vanished, b.cut_vanished);
                    for (x, y) in [
                        (a.tau, b.tau),
                        (a.energy.total, b.energy.total),
                        (a.max_speed, b.max_speed),
                        (a.max_length_speed, b.max_length_speed),
                        (a.roughening_rate, b.roughening_rate),
                        (a.cut_motion, b.cut_motion),
                    ] {
                        assert!(
                            (x - y).abs() <= 1e-10 * x.abs().max(1.0),
                            "{}, {name}, step {step}: {x} against {y}",
                            shape.name
                        );
                    }
                }
                assert_eq!(
                    many.regions(),
                    two.evaluate(problem, &options.eval)
                        .unwrap()
                        .cut
                        .vertex_region
                );
            }
        }
        println!(
            "{}: two regions as [phi, -phi], largest relative difference of the functions over all runs {worst:.1e}",
            shape.name
        );
    }

    // The cut and the evaluation are the same bit for bit.
    let cut = Cut::from_level_set(&ball.mesh, &ball_phi, 1e-3).unwrap();
    let minus: Vec<f64> = ball_phi.iter().map(|x| -x).collect();
    let network = Cut::from_regions(&ball.mesh, &[ball_phi.clone(), minus.clone()], 1e-3).unwrap();
    assert_eq!(cut, network);
    let problem = CutProblem::new(Energy::Hencky, 0.3);
    let mut evaluator = ball.evaluator();
    let a = evaluator
        .evaluate(&ball_phi, &problem, &EvalOptions::default())
        .unwrap();
    let b = evaluator
        .evaluate_regions(
            &[ball_phi.clone(), minus],
            &problem,
            &EvalOptions::default(),
        )
        .unwrap();
    assert_eq!(a.speed, b.speed);
    assert_eq!(a.energy, b.energy);
    assert_eq!(a.gradient, b.gradient);
}

/// A flat torus of width `w` and height `h` as a mesh with edge lengths and nothing else:
/// `nx × ny` rectangles, each split by a diagonal. Returns the coordinates of the vertices in
/// the chart `[0, w) × [0, h)` too, for setting up functions and measuring.
fn flat_torus(nx: usize, ny: usize, w: f64, h: f64) -> (TriMesh, Geometry, Vec<[f64; 2]>) {
    let index = |i: usize, j: usize| ((i % nx) + nx * (j % ny)) as u32;
    let mut faces = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let (a, b, c, d) = (
                index(i, j),
                index(i + 1, j),
                index(i + 1, j + 1),
                index(i, j + 1),
            );
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    let mesh = TriMesh::new(nx * ny, faces).unwrap();
    let (dx, dy) = (w / nx as f64, h / ny as f64);
    let chart: Vec<[f64; 2]> = (0..nx * ny)
        .map(|v| [(v % nx) as f64 * dx, (v / nx) as f64 * dy])
        .collect();
    let lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            let d = wrap(
                [
                    chart[b as usize][0] - chart[a as usize][0],
                    chart[b as usize][1] - chart[a as usize][1],
                ],
                w,
                h,
            );
            (d[0] * d[0] + d[1] * d[1]).sqrt()
        })
        .collect();
    let geometry = Geometry::from_edge_lengths(&mesh, lengths).unwrap();
    (mesh, geometry, chart)
}

/// The shortest representative of a difference vector on the torus.
fn wrap(d: [f64; 2], w: f64, h: f64) -> [f64; 2] {
    [d[0] - w * (d[0] / w).round(), d[1] - h * (d[1] / h).round()]
}

/// The honeycomb with three colours on the torus `3d × √3 d` (six hexagons with centres `d`
/// apart, two of each colour, twelve triple points, total edge length `6 √3 d`): for every
/// colour a function that is negative on its hexagons, at the points `x`.
fn honeycomb(d: f64, x: &[[f64; 2]]) -> Vec<Vec<f64>> {
    let (w, h) = (3.0 * d, 3f64.sqrt() * d);
    let normals = [0.0f64, PI / 3.0, 2.0 * PI / 3.0].map(|a| [a.cos(), a.sin()]);
    let mut phis = vec![vec![f64::INFINITY; x.len()]; 3];
    for i in -4..=7i32 {
        for j in -4..=5i32 {
            let centre = [
                (i as f64 + 0.5 * j as f64) * d + 0.123 * d,
                0.5 * 3f64.sqrt() * j as f64 * d + 0.071 * d,
            ];
            let colour = (i - j).rem_euclid(3) as usize;
            for (v, p) in x.iter().enumerate() {
                let r = wrap([p[0] - centre[0], p[1] - centre[1]], w, h);
                let hexagonal = normals
                    .iter()
                    .map(|n| (r[0] * n[0] + r[1] * n[1]).abs())
                    .fold(0.0, f64::max);
                phis[colour][v] = phis[colour][v].min(hexagonal - 0.5 * d);
            }
        }
    }
    phis
}

/// The largest deviation from 120° of the angles under which the three curves leave the triple
/// points, in degrees, and the number of triple points. The direction of a curve is taken
/// between its cut points number `reach` and `2 · reach` from the junction, so the place of the
/// triple point inside its face (the centroid of three cut points, up to a third of a mesh
/// width off the junction of the three curves) does not enter.
fn junction_angles(
    mesh: &TriMesh,
    cut: &Cut,
    chart: &[[f64; 2]],
    w: f64,
    h: f64,
    reach: usize,
) -> (f64, usize) {
    let at_point = |p: usize| {
        let point = cut.points[p];
        let [a, b] = mesh.edge_vertices(point.edge as usize);
        let d = wrap([chart[b][0] - chart[a][0], chart[b][1] - chart[a][1]], w, h);
        [chart[a][0] + point.t * d[0], chart[a][1] + point.t * d[1]]
    };
    let at_triple = |m: usize| {
        let triple = cut.triple_points[m];
        let face = mesh.faces()[triple.face as usize];
        let origin = chart[face[0] as usize];
        let mut x = origin;
        for c in 1..3 {
            let d = wrap(
                [
                    chart[face[c] as usize][0] - origin[0],
                    chart[face[c] as usize][1] - origin[1],
                ],
                w,
                h,
            );
            x[0] += triple.barycentric[c] * d[0];
            x[1] += triple.barycentric[c] * d[1];
        }
        x
    };
    // The curves need positions only for their order; give them the chart.
    let positions: Vec<[f64; 3]> = chart.iter().map(|p| [p[0], p[1], 0.0]).collect();
    let curves = cut.curves(mesh, &positions).unwrap();
    let _ = &at_triple;
    let mut directions: Vec<Vec<f64>> = vec![Vec::new(); cut.triple_points.len()];
    for curve in &curves {
        let n = curve.points.len();
        let (near, far) = (reach.min(n - 1), (2 * reach).min(n - 1));
        for (end, from, to) in [(0usize, near, far), (1, n - 1 - near, n - 1 - far)] {
            if let CurveEnd::Triple(m) = curve.ends[end] {
                let origin = at_point(curve.points[from] as usize);
                let target = at_point(curve.points[to] as usize);
                let d = wrap([target[0] - origin[0], target[1] - origin[1]], w, h);
                directions[m as usize].push(d[1].atan2(d[0]));
            }
        }
    }
    let mut worst = 0.0f64;
    for angles in &mut directions {
        assert_eq!(angles.len(), 3);
        angles.sort_by(f64::total_cmp);
        for k in 0..3 {
            let gap = (angles[(k + 1) % 3] - angles[k]).rem_euclid(2.0 * PI);
            worst = worst.max((gap.to_degrees() - 120.0).abs());
        }
    }
    (worst, directions.len())
}

/// SPEC §8.6: pure length flow of a honeycomb on a flat torus that is given by its edge
/// lengths alone. The honeycomb is stationary; a perturbed one returns: its curves straighten
/// and its junctions go back toward 120°.
///
/// What is not met, by either stage set: the length of the discrete network is not monotone
/// (the triple point is the centroid of three cut points and hops from face to face; single
/// steps raise the length by up to 0.2 % while it falls by 2 % overall), and at this mesh
/// width the junctions settle 5° to 9° off 120° over the first few segments, with the network
/// 0.7 % to 0.8 % longer than the honeycomb. Both errors fall when the mesh width is halved
/// (`examples/network_study.rs honeycomb 84`; README).
#[test]
fn honeycomb_on_a_flat_torus() {
    // Unit area: 3 √3 d² = 1.
    let d = (1.0 / (3.0 * 3f64.sqrt())).sqrt();
    let (w, h) = (3.0 * d, 3f64.sqrt() * d);
    let exact = 6.0 * 3f64.sqrt() * d;
    let (nx, ny) = (42, 24);
    let (mesh, geometry, chart) = flat_torus(nx, ny, w, h);
    let width = w / nx as f64;
    let problem = CutProblem::for_regions(Energy::Dirichlet, 1.0, 3);
    let stage_sets = [
        ("default stages", fixed(5e-4)),
        ("stages of the paper", paper_stages(5e-4)),
    ];

    for (name, options) in &stage_sets {
        // Stationary.
        let mut flow = RegionFlow::new(&mesh, &geometry, &honeycomb(d, &chart)).unwrap();
        flow.redistance().unwrap();
        let start = flow.evaluate(&problem, &options.eval).unwrap();
        assert_eq!(start.report.n_triple_points, 12);
        assert!(start.energy.distortion.iter().all(|j| j.abs() < 1e-20));
        let (angle, _) = junction_angles(&mesh, &start.cut, &chart, w, h, 3);
        let mut moved = 0.0f64;
        for _ in 0..40 {
            let report = flow.step(&problem, options).unwrap();
            assert!(report.max_speed < 1e-9);
            moved = moved.max(report.cut_motion);
        }
        let end = flow.evaluate(&problem, &options.eval).unwrap();
        let (angle_end, _) = junction_angles(&mesh, &end.cut, &chart, w, h, 3);
        println!(
            "{name}: honeycomb at rest: length {:.5} -> {:.5} (exact {exact:.5}), junction angles off by {angle:.2} -> {angle_end:.2} degrees, largest motion of a step {:.1e} mesh widths, {} triple points",
            start.energy.length,
            end.energy.length,
            moved / width,
            end.report.n_triple_points
        );
        assert_eq!(end.report.n_triple_points, 12);
        // The triple point is the centroid of three cut points, not the junction of the three
        // curves: the polyline has a kink of the size of a face at every junction, and is longer
        // than the honeycomb by a multiple of the mesh width (1.0 % here, 0.5 % at rest).
        assert!((start.energy.length / exact - 1.0).abs() < 1.5e-2);
        assert!((end.energy.length / exact - 1.0).abs() < 1.5e-2, "{name}");
        assert!(angle_end < 6.0, "{name}: {angle_end}");

        // Perturbed: every point is displaced by up to 1.3 mesh widths.
        let displaced: Vec<[f64; 2]> = chart
            .iter()
            .map(|p| {
                [
                    p[0] + 0.04 * (2.0 * PI * p[1] / h).sin(),
                    p[1] + 0.03 * (4.0 * PI * p[0] / w).cos(),
                ]
            })
            .collect();
        let mut flow = RegionFlow::new(&mesh, &geometry, &honeycomb(d, &displaced)).unwrap();
        flow.redistance().unwrap();
        let start = flow.evaluate(&problem, &options.eval).unwrap();
        let (angle, _) = junction_angles(&mesh, &start.cut, &chart, w, h, 3);
        let mut lengths = vec![start.energy.length];
        for _ in 0..160 {
            let report = flow.step(&problem, options).unwrap();
            lengths.push(report.energy.length);
        }
        let end = flow.evaluate(&problem, &options.eval).unwrap();
        lengths.push(end.energy.length);
        let (angle_end, _) = junction_angles(&mesh, &end.cut, &chart, w, h, 3);
        let rises: Vec<f64> = lengths
            .windows(2)
            .map(|pair| pair[1] - pair[0])
            .filter(|d| *d > 0.0)
            .collect();
        let rise = rises.iter().fold(0.0f64, |m, d| m.max(*d));
        println!(
            "{name}: perturbed honeycomb: length {:.5} -> {:.5} (exact {exact:.5}), junction angles off by {angle:.1} -> {angle_end:.1} degrees, the length rose in {} of {} steps, by at most {:.2} % of the length, {} triple points",
            start.energy.length,
            end.energy.length,
            rises.len(),
            lengths.len() - 1,
            100.0 * rise / exact,
            end.report.n_triple_points
        );
        assert_eq!(end.report.n_triple_points, 12, "{name}");
        assert!(start.energy.length > 1.01 * exact);
        // Measured: +0.71 % and 8.7° (default stages), +0.78 % and 5.1° (stages of the paper),
        // from +2.74 % and 28°.
        assert!((end.energy.length / exact - 1.0).abs() < 1.2e-2, "{name}");
        assert!(
            angle_end < 12.0 && angle_end < 0.5 * angle,
            "{name}: {angle_end}"
        );
        assert!(rise < 3e-3 * exact, "{name}: the length rose by {rise}");
    }
}

/// The level-set functions of the partition of the sphere into the cells of the directions
/// `axes`: region `k` is where `x · axes[k]` is largest.
fn cells(shape: &Shape, axes: &[[f64; 3]]) -> Vec<Vec<f64>> {
    axes.iter()
        .map(|a| shape.sample(|p| -(p[0] * a[0] + p[1] * a[1] + p[2] * a[2])))
        .collect()
}

/// SPEC §8.6: the tetrahedral partition of the sphere (4 regions, 4 triple points) and the
/// cubical one (6 regions, 8 triple points) are equilibria by symmetry. The mesh does not have
/// their symmetry, so they are at rest up to the discretization error: over 30 steps the
/// regions keep their areas and the triple points stay where they are.
#[test]
fn symmetric_partitions_of_the_sphere_are_stationary() {
    let s = sphere(3);
    let r = 1.0 / 3f64.sqrt();
    let tetrahedron = vec![[r, r, r], [r, -r, -r], [-r, r, -r], [-r, -r, r]];
    let cube = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    // A rotation that puts the partitions in no special position relative to the mesh.
    let rotate = |a: [f64; 3]| {
        let (c1, s1) = (0.4f64.cos(), 0.4f64.sin());
        let (c2, s2) = (0.7f64.cos(), 0.7f64.sin());
        let b = [c1 * a[0] - s1 * a[1], s1 * a[0] + c1 * a[1], a[2]];
        [b[0], c2 * b[1] - s2 * b[2], s2 * b[1] + c2 * b[2]]
    };
    for (name, axes, triples, alpha) in [
        ("tetrahedral", tetrahedron, 4usize, 3.0),
        ("cubical", cube, 8, 2.0),
    ] {
        let axes: Vec<[f64; 3]> = axes.into_iter().map(rotate).collect();
        let n = axes.len();
        let problem = CutProblem::for_regions(Energy::Dirichlet, alpha, n);
        for (stages, options) in [
            ("default stages", fixed(1e-4)),
            ("stages of the paper", paper_stages(1e-4)),
        ] {
            let mut flow = RegionFlow::new(&s.mesh, &s.geometry, &cells(&s, &axes)).unwrap();
            flow.redistance().unwrap();
            let start = flow.evaluate(&problem, &options.eval).unwrap();
            assert_eq!(start.report.n_triple_points, triples);
            let before = start.cut.triple_positions(&s.mesh, &s.positions).unwrap();
            let mut speed = 0.0f64;
            for _ in 0..30 {
                speed = speed.max(flow.step(&problem, &options).unwrap().max_speed);
            }
            let end = flow.evaluate(&problem, &options.eval).unwrap();
            assert_eq!(end.report.n_triple_points, triples, "{name}, {stages}");
            let after = end.cut.triple_positions(&s.mesh, &s.positions).unwrap();
            let moved = before
                .iter()
                .map(|p| {
                    after
                        .iter()
                        .map(|q| (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt())
                        .fold(f64::INFINITY, f64::min)
                })
                .fold(0.0, f64::max);
            let area_error = end
                .energy
                .area
                .iter()
                .map(|a| (a * n as f64 - 1.0).abs())
                .fold(0.0, f64::max);
            println!(
                "{name}, {stages}: areas off by {:.1e} relative, triple points moved by {moved:.3} (mesh width 0.15), length {:.4} -> {:.4}, energy {:.4} -> {:.4}, largest |s| {speed:.2}",
                area_error,
                start.energy.length,
                end.energy.length,
                start.energy.total,
                end.energy.total
            );
            assert!(area_error < 0.05, "{name}, {stages}: {area_error}");
            assert!(moved < 0.15, "{name}, {stages}: {moved}");
            assert!((end.energy.length / start.energy.length - 1.0).abs() < 0.02);
        }
    }
}

/// SPEC §8.6: the energy decreases along the flow from a generic start with several regions,
/// for every energy and both stage sets.
#[test]
fn energy_decreases_from_generic_starts() {
    let s = sphere(3);
    // Five cells of directions in no special position, with a ripple.
    let mut rng = SplitMix64::new(3);
    let axes: Vec<[f64; 3]> = (0..5)
        .map(|_| {
            let v = [0, 1, 2].map(|_| rng.uniform(-1.0, 1.0));
            let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            v.map(|x| x / norm)
        })
        .collect();
    let mut phis = cells(&s, &axes);
    for (k, phi) in phis.iter_mut().enumerate() {
        for (x, p) in phi.iter_mut().zip(&s.positions) {
            *x += 0.05 * ((3 + k) as f64 * p[0]).sin() * (2.0 * p[2]).cos();
        }
    }
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
    for (name, energy) in &energies {
        for stages in ["default stages", "stages of the paper"] {
            let mut flow = RegionFlow::new(&s.mesh, &s.geometry, &phis).unwrap();
            flow.redistance().unwrap();
            let mut problem = CutProblem::for_regions(energy.clone(), 0.0, 5);
            let omega = flow
                .roughening_rate(&problem, &EvalOptions::default())
                .unwrap();
            problem.length_weight = ddg_variational_cuts::length_weight_for_ripples(omega, 0.15);
            problem.regions[2] = RegionWeights {
                distortion: 0.5,
                area: 0.1 * omega,
            };
            let step = StepSize::Capped {
                tau_max: 0.01 / omega,
                max_motion: 0.1,
            };
            let mut options = FlowOptions::new(step);
            if stages != "default stages" {
                options = FlowOptions::paper(step);
                options.eval = EvalOptions::default();
            }
            let report = flow
                .run(
                    &problem,
                    &options,
                    &StopRule {
                        max_steps: 15,
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
                "{name}, {stages}: energy {:.6e} -> {:.6e}, largest change of a step {:+.2e} of the total drop",
                energy_of[0],
                energy_of[energy_of.len() - 1],
                worst / drop
            );
            assert!(drop > 0.0, "{name}, {stages}: {energy_of:?}");
            assert!(worst <= 0.0, "{name}, {stages}: {energy_of:?}");
        }
    }
}

/// SPEC §8.6: a step does not depend on the numbering of the vertices and faces or on the
/// order of the regions, and scales with the mesh.
///
/// What limits this (the tie rules): a vertex at which the two smallest functions are equal
/// goes to the region with the lower index, so a permutation of the regions can move it; the
/// diagonal of a cut quad is chosen by the roles of its corners when both are equally long
/// (SPEC §3.3); fast marching accepts the vertex with the lower index among equal values
/// (§3.12). None of the three occurs on a mesh in general position, which this one is.
#[test]
fn a_step_is_invariant_under_renumbering_and_region_permutation() {
    let mut g = ddg_core::generators::icosphere(2).unwrap();
    ddg_core::generators::perturb(&mut g.positions, [0.01, 0.01, 0.01], 5);
    let s = Shape::from_generated("perturbed sphere", g);
    let n = s.positions.len();
    let axes = [
        [0.8, 0.5, 0.33],
        [-0.6, 0.7, -0.4],
        [0.1, -0.9, 0.42],
        [-0.3, -0.2, -0.93],
    ];
    let phis = cells(&s, &axes);
    let problem = CutProblem {
        penalty: Some(s.sample(|p| 0.5 + 0.3 * p[0])),
        regions: (0..4)
            .map(|k| RegionWeights {
                distortion: 1.0 + 0.25 * k as f64,
                area: 0.1 * k as f64,
            })
            .collect(),
        ..CutProblem::new(Energy::Hencky, 0.4)
    };
    for (name, options) in [
        ("default stages", fixed(2e-4)),
        ("stages of the paper", paper_stages(2e-4)),
    ] {
        let mut flow = RegionFlow::new(&s.mesh, &s.geometry, &phis).unwrap();
        flow.redistance().unwrap();
        let start = flow.level_sets();
        assert!(
            flow.evaluate(&problem, &options.eval)
                .unwrap()
                .report
                .n_triple_points
                >= 4
        );
        for _ in 0..3 {
            flow.step(&problem, &options).unwrap();
        }
        let reference = flow.level_sets();
        let size = reference[0].iter().fold(0.0f64, |m, x| m.max(x.abs()));

        // The regions in another order.
        let order = [2usize, 0, 3, 1];
        let permuted: Vec<Vec<f64>> = order.iter().map(|&k| start[k].clone()).collect();
        let permuted_problem = CutProblem {
            regions: order.iter().map(|&k| problem.regions[k]).collect(),
            ..problem.clone()
        };
        let mut flow = RegionFlow::new(&s.mesh, &s.geometry, &permuted).unwrap();
        for _ in 0..3 {
            flow.step(&permuted_problem, &options).unwrap();
        }
        let got = flow.level_sets();
        let worst = (0..4)
            .map(|k| largest_difference(&got[k], &reference[order[k]]))
            .fold(0.0, f64::max);
        println!(
            "{name}: regions permuted: largest difference {:.1e} of {size:.2}",
            worst
        );
        assert!(worst < 1e-9 * size, "{name}: {worst}");

        // The vertices and faces in another order.
        let mut rng = SplitMix64::new(11);
        let mut new_index: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            new_index.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let mut positions = vec![[0.0; 3]; n];
        let mut moved = vec![vec![0.0; n]; 4];
        let mut penalty = vec![0.0; n];
        for i in 0..n {
            positions[new_index[i]] = s.positions[i];
            penalty[new_index[i]] = problem.penalty.as_ref().unwrap()[i];
            for k in 0..4 {
                moved[k][new_index[i]] = start[k][i];
            }
        }
        let mut faces: Vec<[u32; 3]> = s
            .mesh
            .faces()
            .iter()
            .enumerate()
            .map(|(f, face)| [0, 1, 2].map(|c| new_index[face[(c + f % 3) % 3] as usize] as u32))
            .collect();
        faces.reverse();
        let relabeled = Shape::new("relabeled", positions, faces);
        let relabeled_problem = CutProblem {
            penalty: Some(penalty),
            ..problem.clone()
        };
        let mut flow = RegionFlow::new(&relabeled.mesh, &relabeled.geometry, &moved).unwrap();
        for _ in 0..3 {
            flow.step(&relabeled_problem, &options).unwrap();
        }
        let got = flow.level_sets();
        let worst = (0..4)
            .map(|k| {
                (0..n)
                    .map(|i| (got[k][new_index[i]] - reference[k][i]).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max);
        println!(
            "{name}: vertices renumbered: largest difference {:.1e} of {size:.2}",
            worst
        );
        assert!(worst < 1e-9 * size, "{name}: {worst}");

        // Three times as large.
        let scaled = Shape::new(
            "scaled",
            s.positions.iter().map(|p| p.map(|x| 3.0 * x)).collect(),
            s.mesh.faces().to_vec(),
        );
        let tripled: Vec<Vec<f64>> = start
            .iter()
            .map(|phi| phi.iter().map(|x| 3.0 * x).collect())
            .collect();
        let mut flow = RegionFlow::new(&scaled.mesh, &scaled.geometry, &tripled).unwrap();
        for _ in 0..3 {
            flow.step(&problem, &options).unwrap();
        }
        let got = flow.level_sets();
        let worst = (0..4)
            .map(|k| {
                (0..n)
                    .map(|i| (got[k][i] - 3.0 * reference[k][i]).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max);
        println!(
            "{name}: scaled: largest difference {:.1e} of {:.2}",
            worst,
            3.0 * size
        );
        assert!(worst < 1e-9 * size, "{name}: {worst}");
    }
}

/// SPEC §8.2, §8.3, §8.7: the network of a state with several regions. Counts, the cut mesh,
/// the flux through triple points, the curves, the helpers, and bad input.
#[test]
fn network_of_several_regions() {
    let s = sphere(3);
    let r = 1.0 / 3f64.sqrt();
    let axes = [[r, r, r], [r, -r, -r], [-r, r, -r], [-r, -r, r]].map(|a| {
        let (c, sn) = (0.3f64.cos(), 0.3f64.sin());
        [c * a[0] - sn * a[2], a[1], sn * a[0] + c * a[2]]
    });
    let phis = cells(&s, &axes);
    let cut = Cut::from_regions(&s.mesh, &phis, 1e-3).unwrap();
    assert_eq!(cut.triple_points.len(), 4);
    assert_eq!(
        cut.vertex_region,
        ddg_variational_cuts::regions_of(&s.mesh, &phis).unwrap()
    );
    // Every triple point lies in a face with three regions, is the centroid of its three cut
    // points, and has three segments.
    for (m, triple) in cut.triple_points.iter().enumerate() {
        let face = s.mesh.faces()[triple.face as usize];
        let mut regions: Vec<u32> = face
            .iter()
            .map(|&v| cut.vertex_region[v as usize])
            .collect();
        regions.sort_unstable();
        regions.dedup();
        assert_eq!(regions.len(), 3);
        assert!((triple.barycentric.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        assert!(triple.barycentric.iter().all(|x| *x > 0.0));
        let node = (cut.points.len() + m) as u32;
        let ends: Vec<u32> = cut
            .segments
            .iter()
            .filter(|segment| segment.points[1] == node)
            .map(|segment| segment.points[0])
            .collect();
        assert_eq!(ends, triple.points.to_vec());
    }
    for (i, segment) in cut.segments.iter().enumerate() {
        let [lower, higher] = cut.segment_regions(&s.mesh, i);
        assert!(lower < higher);
        assert!((segment.points[0] as usize) < cut.points.len());
    }

    // The evaluation: the cut mesh is a valid mesh with one piece per region; the areas add
    // up; the flux of every region is the curvature it contains.
    let mut evaluator = s.evaluator();
    let problem = CutProblem::for_regions(Energy::Hencky, 0.5, 4);
    let e = evaluator
        .evaluate_regions(&phis, &problem, &EvalOptions::default())
        .unwrap();
    assert_eq!(e.report.n_triple_points, 4);
    assert!((e.energy.area.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    for area in &e.energy.area {
        assert!((area - 0.25).abs() < 0.02);
    }
    let cut_mesh = e.cut_mesh().unwrap();
    assert_eq!(cut_mesh.mesh.n_components(), 4);
    assert_eq!(cut_mesh.mesh.euler_characteristic(), 4);
    let tri = &e.triangles;
    assert_eq!(
        tri.n_vertices(),
        s.positions.len() + 2 * cut.points.len() + 12
    );
    let at = cut_mesh.positions(&e.cut, &s.mesh, &s.positions).unwrap();
    let triples = e.cut.triple_positions(&s.mesh, &s.positions).unwrap();
    for m in 0..4 {
        for v in tri.triple_vertices(m) {
            assert_eq!(at[v], triples[m]);
        }
        // The three copies lie in the three regions of the face.
        let mut regions = tri
            .triple_vertices(m)
            .map(|v| tri.vertex_region[v])
            .to_vec();
        regions.sort_unstable();
        regions.dedup();
        assert_eq!(regions.len(), 3);
    }
    // Gauss–Bonnet per region: the fluxes of a region sum to minus its curvature (4π/4 each,
    // up to the mesh), and the recovered normal derivative keeps that total.
    let mut flux = vec![0.0; 4];
    let mut recovered = vec![0.0; 4];
    for v in s.positions.len()..tri.n_vertices() {
        flux[tri.vertex_region[v] as usize] += e.flux[v];
        recovered[tri.vertex_region[v] as usize] += e.dual_length[v] * e.normal_derivative[v];
    }
    for k in 0..4 {
        assert!((flux[k] - PI).abs() < 0.1, "{flux:?}");
        assert!(
            (flux[k] - recovered[k]).abs() < 1e-10,
            "{flux:?} {recovered:?}"
        );
    }
    assert!((flux.iter().sum::<f64>() - 4.0 * PI).abs() < 1e-9);
    // The length of the tetrahedral net on the sphere of unit area: six arcs of
    // arccos(−1/3) · R.
    let exact = 6.0 * (-1.0f64 / 3.0).acos() / (4.0 * PI).sqrt();
    assert!(
        (e.energy.length / exact - 1.0).abs() < 0.01,
        "{}",
        e.energy.length
    );

    // The curves: six, between the four triple points, each between two regions.
    let curves = e.cut.curves(&s.mesh, &s.positions).unwrap();
    assert_eq!(curves.len(), 6);
    let mut pairs: Vec<[u32; 2]> = curves.iter().map(|c| c.regions).collect();
    pairs.sort_unstable();
    assert_eq!(pairs, vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]);
    let mut seen = vec![0usize; e.cut.points.len()];
    for curve in &curves {
        assert!(
            curve
                .ends
                .iter()
                .all(|end| matches!(end, CurveEnd::Triple(_)))
        );
        assert_eq!(curve.positions.len(), curve.points.len() + 2);
        for &p in &curve.points {
            seen[p as usize] += 1;
            assert_eq!(e.cut.point_regions(&s.mesh, p as usize), curve.regions);
        }
        if curve.points.len() > 1 {
            assert!(curve.points[0] < curve.points[curve.points.len() - 1]);
        }
    }
    assert!(seen.iter().all(|&count| count == 1));

    // Helpers: states from labels.
    let labels = ddg_variational_cuts::regions_of(&s.mesh, &phis).unwrap();
    let from_labels = init::regions_from_vertex_labels(&s.mesh, &s.geometry, &labels, 4).unwrap();
    assert_eq!(
        ddg_variational_cuts::regions_of(&s.mesh, &from_labels).unwrap(),
        labels
    );
    let network = Cut::from_regions(&s.mesh, &from_labels, 0.0).unwrap();
    assert_eq!(network.triple_points.len(), 4);
    assert_eq!(network.edge_point, cut.edge_point);
    // An empty region is allowed.
    let five = init::regions_from_vertex_labels(&s.mesh, &s.geometry, &labels, 5).unwrap();
    assert_eq!(five.len(), 5);
    assert!(five[4].iter().all(|x| *x > 0.0));
    assert_eq!(
        ddg_variational_cuts::regions_of(&s.mesh, &five).unwrap(),
        labels
    );
    // Face labels: a vertex takes the label that covers most of its angle.
    let face_labels: Vec<u32> = s
        .mesh
        .faces()
        .iter()
        .map(|face| {
            let c = [0, 1, 2].map(|k| {
                face.iter()
                    .map(|&v| s.positions[v as usize][k])
                    .sum::<f64>()
            });
            (0..4)
                .max_by(|&a, &b| {
                    let dot = |k: usize| (0..3).map(|x| c[x] * axes[k][x]).sum::<f64>();
                    dot(a).total_cmp(&dot(b))
                })
                .unwrap() as u32
        })
        .collect();
    let vertex_labels =
        init::region_labels_from_faces(&s.mesh, &s.geometry, &face_labels, 4).unwrap();
    let differing = vertex_labels
        .iter()
        .zip(&labels)
        .filter(|(a, b)| a != b)
        .count();
    assert!(
        differing * 20 < labels.len(),
        "{differing} of {}",
        labels.len()
    );
    let from_faces = init::regions_from_face_labels(&s.mesh, &s.geometry, &face_labels, 4).unwrap();
    assert_eq!(
        ddg_variational_cuts::regions_of(&s.mesh, &from_faces).unwrap(),
        vertex_labels
    );
    // Clusters of normals: on the sphere they are caps around the cluster centres.
    let clusters = init::normal_clusters(&s.mesh, &s.positions, 6, 50).unwrap();
    let mut sizes = vec![0usize; 6];
    for &c in &clusters {
        sizes[c as usize] += 1;
    }
    println!("six clusters of normals on the sphere: sizes {sizes:?}");
    assert!(
        sizes
            .iter()
            .all(|&size| size * 12 > clusters.len() && size * 4 < clusters.len())
    );
    assert_eq!(
        clusters,
        init::normal_clusters(&s.mesh, &s.positions, 6, 50).unwrap()
    );

    // Bad input.
    assert!(matches!(
        Cut::from_regions(&s.mesh, &phis[..1], 1e-3),
        Err(CutError::UnsupportedRegionCount { got: 1 })
    ));
    let mut short = phis.to_vec();
    short[2].pop();
    assert!(matches!(
        RegionFlow::new(&s.mesh, &s.geometry, &short),
        Err(CutError::WrongLength { .. })
    ));
    let mut bad = phis.to_vec();
    bad[3][7] = f64::NAN;
    assert!(matches!(
        RegionFlow::new(&s.mesh, &s.geometry, &bad),
        Err(CutError::NonFiniteLevelSet { vertex: 7, .. })
    ));
    let mut flow = RegionFlow::new(&s.mesh, &s.geometry, &phis).unwrap();
    let error = flow
        .step(&CutProblem::new(Energy::Dirichlet, 1.0), &fixed(1e-4))
        .unwrap_err();
    assert!(
        matches!(
            error,
            CutError::WrongLength {
                expected: 4,
                got: 2,
                ..
            }
        ),
        "{error}"
    );
    assert!(matches!(
        flow.step(&problem, &fixed(-1.0)),
        Err(CutError::InvalidOptions { .. })
    ));
    assert_eq!(flow.n_steps(), 0);
    let mut wrong = labels.clone();
    wrong[9] = 4;
    let error = init::regions_from_vertex_labels(&s.mesh, &s.geometry, &wrong, 4).unwrap_err();
    assert!(error.to_string().contains("vertex 9"), "{error}");
    assert!(matches!(
        init::normal_clusters(&s.mesh, &s.positions, 0, 10),
        Err(CutError::InvalidOptions { .. })
    ));
    // The two-region evaluation still refuses three regions' worth of weights.
    assert!(matches!(
        s.evaluator()
            .evaluate(&phis[0], &problem, &EvalOptions::default()),
        Err(CutError::UnsupportedRegionCount { got: 4 })
    ));
}
