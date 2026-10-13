//! Experiments with the flow for several regions (SPEC §8): the honeycomb on a flat torus,
//! atlases from clusters of normals, different starts, timings.
//!
//! ```text
//! cargo run --release -p ddg-variational-cuts --example network_study -- honeycomb [NX] [STAGES] [TAU] [STEPS]
//! ```

use std::f64::consts::PI;

use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CurveEnd, Cut, CutProblem, Energy, EvalOptions, FlowOptions, RegionFlow, StepSize,
};

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

fn options_for(stages: &str, step: StepSize) -> FlowOptions<'static> {
    if stages == "paper" {
        let mut options = FlowOptions::paper(step);
        options.eval = EvalOptions::default();
        options
    } else {
        FlowOptions::new(step)
    }
}

/// `honeycomb [NX] [STAGES] [TAU] [STEPS] [AMPLITUDE]`: the pure length flow of a perturbed
/// honeycomb, one line every few steps.
fn honeycomb_run() {
    let arg = |k: usize, default: &str| std::env::args().nth(k).unwrap_or(default.to_string());
    let nx: usize = arg(2, "42").parse().unwrap();
    let stages = arg(3, "default");
    let tau: f64 = arg(4, "5e-4").parse().unwrap();
    let steps: usize = arg(5, "160").parse().unwrap();
    let amplitude: f64 = arg(6, "1").parse().unwrap();
    let d = (1.0 / (3.0 * 3f64.sqrt())).sqrt();
    let (w, h) = (3.0 * d, 3f64.sqrt() * d);
    let exact = 6.0 * 3f64.sqrt() * d;
    let ny = (nx as f64 * h / w).round() as usize;
    let (mesh, geometry, chart) = flat_torus(nx, ny, w, h);
    let displaced: Vec<[f64; 2]> = chart
        .iter()
        .map(|p| {
            [
                p[0] + amplitude * 0.04 * (2.0 * PI * p[1] / h).sin(),
                p[1] + amplitude * 0.03 * (4.0 * PI * p[0] / w).cos(),
            ]
        })
        .collect();
    let problem = CutProblem::for_regions(Energy::Dirichlet, 1.0, 3);
    let options = options_for(&stages, StepSize::Fixed(tau));
    let mut flow = RegionFlow::new(&mesh, &geometry, &honeycomb(d, &displaced)).unwrap();
    flow.redistance().unwrap();
    println!(
        "honeycomb on the torus {nx} x {ny}, mesh width {:.4}, exact length {exact:.5}, {stages} stages, tau {tau:e}, displaced by {amplitude} x (0.04, 0.03)",
        w / nx as f64
    );
    let every = (steps / 16).max(1);
    let mut lengths: Vec<f64> = Vec::new();
    for step in 0..=steps {
        if step % every == 0 {
            let e = flow.evaluate(&problem, &options.eval).unwrap();
            let (angle, triples) = junction_angles(&mesh, &e.cut, &chart, w, h, 2);
            let (angle_far, _) = junction_angles(&mesh, &e.cut, &chart, w, h, nx / 10);
            println!(
                "             junction angles off by {angle_far:.2} degrees between the cut points {} and {} from the junction",
                nx / 10,
                2 * (nx / 10)
            );
            let n_points = e.cut.points.len();
            let at_triples: f64 = e
                .cut
                .segments
                .iter()
                .zip(&e.triangles.segment_length)
                .filter(|(segment, _)| segment.points[1] as usize >= n_points)
                .map(|(_, length)| *length)
                .sum();
            let shortest = e
                .triangles
                .segment_length
                .iter()
                .fold(f64::INFINITY, |m, l| m.min(*l));
            println!(
                "             {n_points} cut points, length of the segments at triple points {at_triples:.4}, shortest segment {:.1e} mesh widths",
                shortest / (w / nx as f64)
            );
            println!(
                "  step {step:>4}  length {:.5} ({:+.2} %)  junction angles off by {angle:>5.2} degrees  {triples} triple points  areas {:.4} {:.4} {:.4}",
                e.energy.length,
                100.0 * (e.energy.length / exact - 1.0),
                e.energy.area[0],
                e.energy.area[1],
                e.energy.area[2]
            );
        }
        if step < steps {
            let r = flow.step(&problem, &options).unwrap();
            lengths.push(r.energy.length);
            if step % every == 0 {
                let g = r.gradient_norm.unwrap();
                println!(
                    "             max |lambda| {:.2e}  motion {:.2e}  |grad| {:.3}/{:.3}/{:.3}  clamped {}",
                    r.max_length_speed, r.cut_motion, g.min, g.mean, g.max, r.n_clamped_crossings
                );
            }
        }
    }
    let rises: Vec<f64> = lengths
        .windows(2)
        .map(|pair| pair[1] - pair[0])
        .filter(|d| *d > 0.0)
        .collect();
    println!(
        "  the length rose in {} of {} steps, by at most {:.1e} ({:.2e} of the length), in total {:.1e}",
        rises.len(),
        lengths.len().saturating_sub(1),
        rises.iter().fold(0.0f64, |m, d| m.max(*d)),
        rises.iter().fold(0.0f64, |m, d| m.max(*d)) / exact,
        rises.iter().sum::<f64>()
    );
}

/// `stripes [NX] [STAGES] [TAU] [STEPS]`: three bands around the torus with rippled borders and
/// no triple point. Every border is a curve of its own and must straighten.
fn stripes_run() {
    let arg = |k: usize, default: &str| std::env::args().nth(k).unwrap_or(default.to_string());
    let nx: usize = arg(2, "42").parse().unwrap();
    let stages = arg(3, "default");
    let tau: f64 = arg(4, "5e-4").parse().unwrap();
    let steps: usize = arg(5, "80").parse().unwrap();
    let d = (1.0 / (3.0 * 3f64.sqrt())).sqrt();
    let (w, h) = (3.0 * d, 3f64.sqrt() * d);
    let ny = (nx as f64 * h / w).round() as usize;
    let (mesh, geometry, chart) = flat_torus(nx, ny, w, h);
    // Band k lies between the borders k and k + 1; border k is x = (k + 0.37) w / 3 + ripple.
    let border =
        |k: usize, y: f64| (k as f64 + 0.37) * w / 3.0 + 0.03 * (2.0 * PI * y / h + k as f64).sin();
    let phis: Vec<Vec<f64>> = (0..3)
        .map(|k| {
            chart
                .iter()
                .map(|p| {
                    let left = wrap([p[0] - border(k, p[1]), 0.0], w, h)[0];
                    let right = wrap([p[0] - border((k + 1) % 3, p[1]), 0.0], w, h)[0];
                    // Inside the band: right of its left border and left of its right border.
                    (-left)
                        .max(right)
                        .max(-0.5 * w / 3.0 - left.abs().min(right.abs()) * 0.0)
                })
                .collect()
        })
        .collect();
    let problem = CutProblem::for_regions(Energy::Dirichlet, 1.0, 3);
    let options = options_for(&stages, StepSize::Fixed(tau));
    let mut flow = RegionFlow::new(&mesh, &geometry, &phis).unwrap();
    flow.redistance().unwrap();
    println!(
        "three bands on the torus {nx} x {ny}, exact length at rest {:.5}, {stages} stages, tau {tau:e}; a ripple decays like exp(-{:.1} t)",
        3.0 * h,
        (2.0 * PI / h).powi(2)
    );
    let every = (steps / 10).max(1);
    for step in 0..=steps {
        if step % every == 0 {
            let e = flow.evaluate(&problem, &options.eval).unwrap();
            println!(
                "  step {step:>4}  t {:.4}  length {:.5} ({:+.3} %)  {} triple points  areas {:.4} {:.4} {:.4}",
                flow.time(),
                e.energy.length,
                100.0 * (e.energy.length / (3.0 * h) - 1.0),
                e.report.n_triple_points,
                e.energy.area[0],
                e.energy.area[1],
                e.energy.area[2]
            );
        }
        if step < steps {
            flow.step(&problem, &options).unwrap();
        }
    }
}

struct Model {
    mesh: TriMesh,
    geometry: Geometry,
    positions: Vec<[f64; 3]>,
}

fn load(path: &str) -> Model {
    if let Some(level) = path.strip_prefix("icosphere:") {
        let g = ddg_core::generators::icosphere(level.parse().unwrap()).unwrap();
        let mesh = TriMesh::new(g.positions.len(), g.faces).unwrap();
        let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
        return Model {
            mesh,
            geometry,
            positions: g.positions,
        };
    }
    let obj = ddg_core::io::read_obj(path).unwrap();
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let corners = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..corners.len().saturating_sub(1) {
            faces.push([corners[0], corners[k], corners[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
    Model {
        mesh,
        geometry,
        positions: obj.positions,
    }
}

fn energy_named(name: &str) -> Energy {
    match name {
        "dirichlet" => Energy::Dirichlet,
        "hencky" => Energy::Hencky,
        _ => Energy::ScaledDirichlet(ddg_variational_cuts::ScaleWeight::FeatureSize {
            r_eps: 1.0,
            smoothing: ddg_variational_cuts::PAPER_FEATURE_SMOOTHING,
        }),
    }
}

/// Labels from seeds: every vertex goes to the nearest of `n` seed vertices (straight-line
/// distance). The first seed is `first`, every further one the vertex farthest from the seeds
/// so far.
fn voronoi_labels(positions: &[[f64; 3]], n: usize, first: usize) -> Vec<u32> {
    let dist = |a: [f64; 3], b: [f64; 3]| (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>();
    let mut nearest = vec![f64::INFINITY; positions.len()];
    let mut labels = vec![0u32; positions.len()];
    let mut seed = first;
    for k in 0..n {
        for (v, p) in positions.iter().enumerate() {
            let d = dist(*p, positions[seed]);
            if d < nearest[v] {
                nearest[v] = d;
                labels[v] = k as u32;
            }
        }
        seed = (0..positions.len())
            .max_by(|a, b| nearest[*a].total_cmp(&nearest[*b]))
            .unwrap();
    }
    labels
}

/// One run from vertex labels: what came out, and the cut-adjacent vertices of the result.
struct Outcome {
    line: String,
    adjacent: Vec<bool>,
}

fn run_from_labels(
    model: &Model,
    labels: &[u32],
    n: usize,
    energy: &Energy,
    ripple: f64,
    stages: &str,
    max_steps: usize,
) -> Outcome {
    let phis = ddg_variational_cuts::init::regions_from_vertex_labels(
        &model.mesh,
        &model.geometry,
        labels,
        n,
    )
    .unwrap();
    let mut flow = RegionFlow::new(&model.mesh, &model.geometry, &phis).unwrap();
    let mut problem = CutProblem::for_regions(energy.clone(), 0.0, n);
    // STUDY_MOTION: the largest motion per step as a fraction of an edge (default 0.25).
    let motion: f64 = std::env::var("STUDY_MOTION").map_or(0.25, |v| v.parse().unwrap());
    let options = options_for(
        stages,
        StepSize::Adaptive {
            time_scale: 0.05,
            max_motion: motion,
            tau_max: 1.0,
        },
    );
    let omega = flow.roughening_rate(&problem, &options.eval).unwrap();
    // STUDY_ALPHA: one length weight for every start, instead of the ripple rule.
    problem.length_weight = std::env::var("STUDY_ALPHA").map_or_else(
        |_| ddg_variational_cuts::length_weight_for_ripples(omega, ripple),
        |v| v.parse().unwrap(),
    );
    let start = flow.evaluate(&problem, &options.eval).unwrap();
    let patches =
        |e: &ddg_variational_cuts::Evaluation| e.cut_mesh().map_or(0, |m| m.mesh.n_components());
    let clock = std::time::Instant::now();
    let report = flow
        .run(
            &problem,
            &options,
            &ddg_variational_cuts::StopRule {
                max_steps,
                ..Default::default()
            },
        )
        .unwrap();
    let seconds = clock.elapsed().as_secs_f64();
    let end = flow.evaluate(&problem, &options.eval).unwrap();
    let total = |e: &ddg_variational_cuts::Evaluation| e.energy.distortion.iter().sum::<f64>();
    let rises = report
        .steps
        .windows(2)
        .filter(|pair| pair[1].energy.total > pair[0].energy.total)
        .count();
    let line = format!(
        "alpha_L {:.4e} (omega {:.3e}); {} steps, {:?}, {seconds:.1} s; patches {} -> {}; triple points {} -> {}; cut length {:.4} -> {:.4}; distortion {:.5e} -> {:.5e}; energy {:.5e} -> {:.5e} (rose in {rises} steps)",
        problem.length_weight,
        omega,
        report.steps.len(),
        report.stop,
        patches(&start),
        patches(&end),
        start.report.n_triple_points,
        end.report.n_triple_points,
        start.energy.length,
        end.energy.length,
        total(&start),
        total(&end),
        start.energy.total,
        end.energy.total
    );
    Outcome {
        line,
        adjacent: end.cut.cut_adjacent_vertices(&model.mesh),
    }
}

/// `atlas MESH [N] [RIPPLE] [ENERGY] [STEPS]`: an atlas in the manner of the paper's Fig. 16:
/// `N` clusters of normals as the start, a strong length weight, both stage sets.
fn atlas() {
    let arg = |k: usize, default: &str| std::env::args().nth(k).unwrap_or(default.to_string());
    let path = arg(2, "icosphere:4");
    let n: usize = arg(3, "8").parse().unwrap();
    let ripple: f64 = arg(4, "0.5").parse().unwrap();
    let energy = energy_named(&arg(5, "scaled"));
    let steps: usize = arg(6, "1000").parse().unwrap();
    let model = load(&path);
    let labels =
        ddg_variational_cuts::init::normal_clusters(&model.mesh, &model.positions, n, 100).unwrap();
    println!(
        "{path}: {} vertices; {n} clusters of normals, energy {}, ripples longer than {ripple} grow",
        model.mesh.n_vertices(),
        arg(5, "scaled")
    );
    for stages in ["default", "paper"] {
        let outcome = run_from_labels(&model, &labels, n, &energy, ripple, stages, steps);
        println!("  {stages} stages: {}", outcome.line);
    }
}

/// `starts MESH [N] [RIPPLE] [ENERGY] [STEPS]`: the paper's Fig. 10: different starts on one
/// shape. For every pair of results, the share of the cut-adjacent vertices of one that lie
/// within two edges of a cut-adjacent vertex of the other.
fn starts() {
    let arg = |k: usize, default: &str| std::env::args().nth(k).unwrap_or(default.to_string());
    let path = arg(2, "icosphere:4");
    let n: usize = arg(3, "6").parse().unwrap();
    let ripple: f64 = arg(4, "0.5").parse().unwrap();
    let energy = energy_named(&arg(5, "scaled"));
    let steps: usize = arg(6, "1000").parse().unwrap();
    let model = load(&path);
    let nv = model.mesh.n_vertices();
    let mut rng = ddg_core::generators::SplitMix64::new(17);
    let random: Vec<u32> = (0..nv)
        .map(|_| (rng.next_u64() % n as u64) as u32)
        .collect();
    let kinds: Vec<(&str, Vec<u32>)> = vec![
        (
            "clusters of normals",
            ddg_variational_cuts::init::normal_clusters(&model.mesh, &model.positions, n, 100)
                .unwrap(),
        ),
        (
            "cells of seed vertices (from vertex 0)",
            voronoi_labels(&model.positions, n, 0),
        ),
        (
            "cells of seed vertices (from the middle vertex)",
            voronoi_labels(&model.positions, n, nv / 2),
        ),
        ("random label per vertex", random),
    ];
    let stages = arg(7, "paper");
    println!(
        "{path}: {nv} vertices, {n} regions, energy {}, ripples longer than {ripple} grow, {stages} stages",
        arg(5, "scaled")
    );
    let mut results = Vec::new();
    for (name, labels) in &kinds {
        let outcome = run_from_labels(&model, labels, n, &energy, ripple, &stages, steps);
        println!("  {name}: {}", outcome.line);
        results.push(outcome.adjacent);
    }
    // Proximity: breadth-first rings from the cut-adjacent vertices.
    let near = |adjacent: &[bool]| {
        let mut ring: Vec<bool> = adjacent.to_vec();
        for _ in 0..2 {
            let mut next = ring.clone();
            for &[a, b] in model.mesh.edges() {
                if ring[a as usize] || ring[b as usize] {
                    next[a as usize] = true;
                    next[b as usize] = true;
                }
            }
            ring = next;
        }
        ring
    };
    println!("  share of the cut of the row within two edges of the cut of the column:");
    for (i, a) in results.iter().enumerate() {
        let mut line = format!("    {:<48}", kinds[i].0);
        for b in &results {
            let close = near(b);
            let (inside, all) = a
                .iter()
                .zip(&close)
                .filter(|(x, _)| **x)
                .fold((0usize, 0usize), |(i, n), (_, c)| {
                    (i + usize::from(*c), n + 1)
                });
            line += &format!(" {:>5.1} %", 100.0 * inside as f64 / all.max(1) as f64);
        }
        println!("{line}");
    }
}

/// `timings [MESH...]`: one step for 2, 4 and 8 regions on every mesh, one thread, seconds
/// (median of five steps after the first).
fn timings() {
    let mut paths: Vec<String> = std::env::args().skip(2).collect();
    if paths.is_empty() {
        paths.push("icosphere:6".to_string());
    }
    for path in &paths {
        let model = load(path);
        println!("{path}: {} vertices", model.mesh.n_vertices());
        for n in [2usize, 4, 8] {
            let labels =
                ddg_variational_cuts::init::normal_clusters(&model.mesh, &model.positions, n, 100)
                    .unwrap();
            let phis = ddg_variational_cuts::init::regions_from_vertex_labels(
                &model.mesh,
                &model.geometry,
                &labels,
                n,
            )
            .unwrap();
            for stages in ["default", "paper"] {
                let mut flow = RegionFlow::new(&model.mesh, &model.geometry, &phis).unwrap();
                let mut problem = CutProblem::for_regions(Energy::Hencky, 0.0, n);
                let options = options_for(stages, StepSize::Fixed(1e-7));
                let omega = flow.roughening_rate(&problem, &options.eval).unwrap();
                problem.length_weight =
                    ddg_variational_cuts::length_weight_for_ripples(omega, 0.25);
                let mut rows: Vec<[f64; 5]> = Vec::new();
                let mut points = 0;
                for _ in 0..6 {
                    let clock = std::time::Instant::now();
                    let r = flow.step(&problem, &options).unwrap();
                    let total = clock.elapsed().as_secs_f64();
                    let t = r.timings;
                    rows.push([t.evaluate, t.length, t.extension, t.redistance, total]);
                    points = r.n_cut_points;
                }
                let median = |k: usize| {
                    let mut column: Vec<f64> = rows[1..].iter().map(|r| r[k]).collect();
                    column.sort_by(f64::total_cmp);
                    column[column.len() / 2]
                };
                println!(
                    "  N = {n}, {stages} stages ({points} cut points): evaluate {:.4}  length {:.4}  extension {:.4}  redistance {:.4}  step {:.4}",
                    median(0),
                    median(1),
                    median(2),
                    median(3),
                    median(4)
                );
            }
        }
    }
}

fn main() {
    match std::env::args().nth(1).unwrap_or_default().as_str() {
        "honeycomb" => honeycomb_run(),
        "stripes" => stripes_run(),
        "atlas" => atlas(),
        "starts" => starts(),
        "timings" => timings(),
        other => eprintln!("unknown experiment {other:?}; see the top of the example"),
    }
}
