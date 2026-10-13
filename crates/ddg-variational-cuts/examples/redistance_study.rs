//! Why the shrinking-circle law fails with the redistancing of the paper, and which variant
//! repairs it: the length step of the flow on the flat unit square, with the redistancing done
//! here in several ways from positions (the flow itself runs with `Redistance::Off`).
//!
//! ```text
//! cargo run --release -p ddg-variational-cuts --example redistance_study
//! ```

use ddg_core::distance::upwind_distance;
use ddg_core::generators::{SplitMix64, grid_disk};
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    Cut, CutProblem, Energy, Flow, FlowOptions, Redistance, StepSize, region_of,
};

struct Shape {
    mesh: TriMesh,
    geometry: Geometry,
    positions: Vec<[f64; 3]>,
}

fn square(n: usize, jitter: f64, seed: u64) -> Shape {
    let g = grid_disk(n, n).unwrap();
    let mut positions = g.positions.clone();
    let mut rng = SplitMix64::new(seed);
    let h = 1.0 / n as f64;
    for p in &mut positions {
        let inside = |x: f64| x > 0.5 * h && x < 1.0 - 0.5 * h;
        if inside(p[0]) && inside(p[1]) {
            p[0] += jitter * h * rng.uniform(-1.0, 1.0);
            p[1] += jitter * h * rng.uniform(-1.0, 1.0);
        }
    }
    let mesh = TriMesh::new(positions.len(), g.faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
    Shape {
        mesh,
        geometry,
        positions,
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Mode {
    /// Cut-adjacent vertices keep their values, the others by fast marching (the paper).
    KeepFmm,
    /// Cut-adjacent vertices keep their values, the others get the exact distance to the cut.
    KeepExact,
    /// Every vertex gets the exact distance to the cut polyline.
    AllExact,
    /// Cut-adjacent vertices get their distance to the nearest line of a cut segment in an
    /// incident crossed face (`φ/|∇φ|` of that face), the others by fast marching from them.
    LineFmm,
    /// As `LineFmm` with the distance to the segment instead of its line.
    SegmentFmm,
    /// As `LineFmm`, but the others get the exact distance to the cut.
    LineExact,
    /// The exact circle through the mean radius of the cut: the ideal redistancing.
    Circle,
}

fn dist_segment(x: [f64; 3], p: [f64; 3], q: [f64; 3], clamp: bool) -> f64 {
    let (dx, dy) = (q[0] - p[0], q[1] - p[1]);
    let len2 = dx * dx + dy * dy;
    let mut t = ((x[0] - p[0]) * dx + (x[1] - p[1]) * dy) / len2;
    if clamp {
        t = t.clamp(0.0, 1.0);
    }
    let (cx, cy) = (p[0] + t * dx - x[0], p[1] + t * dy - x[1]);
    (cx * cx + cy * cy).sqrt()
}

fn redistance(s: &Shape, phi: &[f64], mode: Mode, center: [f64; 2]) -> Vec<f64> {
    let cut = Cut::from_level_set(&s.mesh, phi, 0.0).unwrap();
    let points = cut.point_positions(&s.mesh, &s.positions).unwrap();
    let n = phi.len();
    let sign = |v: usize| if region_of(phi[v]) == 0 { -1.0 } else { 1.0 };
    let adjacent = cut.cut_adjacent_vertices(&s.mesh);
    let exact = |v: usize| -> f64 {
        cut.segments
            .iter()
            .map(|seg| {
                dist_segment(
                    s.positions[v],
                    points[seg.points[0] as usize],
                    points[seg.points[1] as usize],
                    true,
                )
            })
            .fold(f64::INFINITY, f64::min)
    };
    if mode == Mode::Circle {
        let radius = |p: [f64; 3]| ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2)).sqrt();
        let mean = points.iter().map(|p| radius(*p)).sum::<f64>() / points.len() as f64;
        return s.positions.iter().map(|p| radius(*p) - mean).collect();
    }
    if mode == Mode::AllExact {
        return (0..n).map(|v| sign(v) * exact(v)).collect();
    }
    // Values of the cut-adjacent vertices.
    let mut out = phi.to_vec();
    if matches!(mode, Mode::LineFmm | Mode::SegmentFmm | Mode::LineExact) {
        let mut best = vec![f64::INFINITY; n];
        for seg in &cut.segments {
            let (p, q) = (
                points[seg.points[0] as usize],
                points[seg.points[1] as usize],
            );
            for &v in &s.mesh.faces()[seg.face as usize] {
                let d = dist_segment(s.positions[v as usize], p, q, mode == Mode::SegmentFmm);
                best[v as usize] = best[v as usize].min(d);
            }
        }
        for v in 0..n {
            if adjacent[v] {
                out[v] = sign(v) * best[v];
            }
        }
    }
    match mode {
        Mode::KeepExact | Mode::LineExact => {
            for v in 0..n {
                if !adjacent[v] {
                    out[v] = sign(v) * exact(v);
                }
            }
        }
        _ => {
            let seeds: Vec<(u32, f64)> = (0..n)
                .filter(|&v| adjacent[v])
                .map(|v| (v as u32, out[v].abs()))
                .collect();
            let region: Vec<u32> = phi.iter().map(|x| region_of(*x)).collect();
            let d = upwind_distance(&s.mesh, &s.geometry, &seeds, Some(&region)).unwrap();
            for v in 0..n {
                if !adjacent[v] {
                    out[v] = sign(v) * d.values[v];
                }
            }
        }
    }
    out
}

fn main() {
    let (rho0, alpha, t_end) = (0.3, 1.0, 0.03);
    let center = [0.511, 0.493];
    let exact = (rho0 * rho0 - 2.0 * alpha * t_end as f64).sqrt();
    let modes = [
        Mode::Circle,
        Mode::KeepFmm,
        Mode::KeepExact,
        Mode::AllExact,
        Mode::LineFmm,
        Mode::SegmentFmm,
        Mode::LineExact,
    ];
    println!(
        "flat unit square, circle rho0 = {rho0}, alpha_L = {alpha}, t = {t_end}: exact radius {exact:.6}"
    );
    println!(
        "entries: relative error of the mean radius / roundness (std of the radius over its mean), in percent"
    );
    for (n, jitter) in [(32, 0.3), (64, 0.3), (128, 0.3), (64, 0.0)] {
        let s = square(n, jitter, 5);
        let radius = |p: [f64; 3]| ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2)).sqrt();
        let phi0: Vec<f64> = s.positions.iter().map(|p| radius(*p) - rho0).collect();
        println!("grid {n} x {n}, jitter {jitter}");
        print!("{:>6} {:>8}", "steps", "tau/h^2");
        for mode in modes {
            print!(" {:>15}", format!("{mode:?}"));
        }
        println!();
        for steps in [6, 24, 96, 384] {
            let tau = t_end / steps as f64;
            print!("{:>6} {:>8.2}", steps, tau * (n * n) as f64);
            for mode in modes {
                let mut flow = Flow::new(&s.mesh, &s.geometry, &phi0).unwrap();
                let problem = CutProblem::new(Energy::Dirichlet, alpha);
                let options = FlowOptions {
                    redistance: Redistance::Off,
                    ..FlowOptions::new(StepSize::Fixed(tau))
                };
                for _ in 0..steps {
                    flow.step(&problem, &options).unwrap();
                    let phi = redistance(&s, &flow.level_set(), mode, center);
                    flow.set_level_set(&phi).unwrap();
                }
                let cut = Cut::from_level_set(&s.mesh, &flow.level_set(), 0.0).unwrap();
                let points = cut.point_positions(&s.mesh, &s.positions).unwrap();
                let r: Vec<f64> = points.iter().map(|p| radius(*p)).collect();
                let mean = r.iter().sum::<f64>() / r.len() as f64;
                let std =
                    (r.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / r.len() as f64).sqrt();
                print!(
                    " {:>+7.2} /{:>6.2}",
                    100.0 * (mean / exact - 1.0),
                    100.0 * std / mean
                );
            }
            println!();
        }
    }
}
