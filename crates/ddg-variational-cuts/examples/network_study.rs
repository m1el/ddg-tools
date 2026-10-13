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

fn main() {
    match std::env::args().nth(1).unwrap_or_default().as_str() {
        "honeycomb" => honeycomb_run(),
        "stripes" => stripes_run(),
        other => eprintln!("unknown experiment {other:?}; see the top of the example"),
    }
}
