//! Shared helpers of the integration tests: meshes with positions, level-set functions, and the
//! finite-difference check of the gradient (SPEC §6 item 3).
#![allow(dead_code)]

use std::f64::consts::PI;

use ddg_core::generators::{GeneratedMesh, bumpy_grid_disk, icosphere};
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{Cut, CutProblem, EvalOptions, Evaluation, Evaluator};

/// A mesh with positions and geometry.
pub struct Shape {
    pub name: String,
    pub positions: Vec<[f64; 3]>,
    pub mesh: TriMesh,
    pub geometry: Geometry,
}

impl Shape {
    pub fn new(name: &str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Shape {
        let mesh = TriMesh::new(positions.len(), faces).expect("valid mesh");
        let geometry = Geometry::from_positions(&mesh, &positions).expect("valid geometry");
        Shape {
            name: name.to_string(),
            positions,
            mesh,
            geometry,
        }
    }

    pub fn from_generated(name: &str, g: GeneratedMesh) -> Shape {
        Shape::new(name, g.positions, g.faces)
    }

    pub fn evaluator(&self) -> Evaluator {
        Evaluator::new(&self.mesh, &self.geometry).expect("evaluator")
    }

    /// A function of position at every vertex.
    pub fn sample(&self, f: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
        self.positions.iter().map(|&p| f(p)).collect()
    }
}

/// The unit sphere, `10 · 4^level + 2` vertices.
pub fn sphere(level: u32) -> Shape {
    Shape::from_generated(&format!("icosphere {level}"), icosphere(level).unwrap())
}

/// The unit square with height `amplitude · sin(3πx) cos(2πy)`, `n × n` cells.
pub fn bumpy(n: usize, amplitude: f64) -> Shape {
    Shape::from_generated(
        &format!("bumpy disk {n}"),
        bumpy_grid_disk(n, n, amplitude).unwrap(),
    )
}

/// Radius of the sphere of unit area.
pub fn unit_area_radius() -> f64 {
    1.0 / (4.0 * PI).sqrt()
}

/// Polar angle of a point of the unit sphere, from the `z` axis.
pub fn polar_angle(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2])
}

/// Shifts `phi` by one of the `candidates` so that no crossing is close to a vertex: returns the
/// shifted function and the smallest crossing offset it has. Finite differences of the energy
/// are only meaningful where no crossing is clamped and none is about to jump to another edge.
pub fn shift_away_from_vertices(shape: &Shape, phi: &[f64], candidates: &[f64]) -> (Vec<f64>, f64) {
    let mut best: Option<(Vec<f64>, f64)> = None;
    for &shift in candidates {
        let shifted: Vec<f64> = phi.iter().map(|x| x + shift).collect();
        let offset = Cut::from_level_set(&shape.mesh, &shifted, 0.0)
            .unwrap()
            .smallest_offset;
        if best.as_ref().is_none_or(|b| offset > b.1) {
            best = Some((shifted, offset));
        }
    }
    best.expect("at least one candidate")
}

/// For every cut point, the rate `dA/dε` at which the region with the lower id gains area
/// through the motion of that cut point when `φ` becomes `φ + ε ψ`. It is exact for the cut
/// mesh: in a face of area `A_f` whose lone vertex `a` has the crossings at the fractions
/// `s_p`, `s_q` of its two edges, the corner triangle has the area `s_p s_q A_f`.
pub fn swept_area_rates(shape: &Shape, e: &Evaluation, phi: &[f64], psi: &[f64]) -> Vec<f64> {
    let mesh = &shape.mesh;
    let cut = &e.cut;
    // dt/dε of every cut point: t = φ_a/(φ_a − φ_b) on the edge [a, b].
    let dt: Vec<f64> = cut
        .points
        .iter()
        .map(|point| {
            let [a, b] = mesh.edge_vertices(point.edge as usize);
            ((1.0 - point.t) * psi[a] + point.t * psi[b]) / (phi[a] - phi[b])
        })
        .collect();
    let mut rates = vec![0.0; cut.points.len()];
    for segment in &cut.segments {
        let f = segment.face as usize;
        let face = mesh.faces()[f];
        let regions = face.map(|v| cut.vertex_region[v as usize]);
        let k = (0..3)
            .find(|&k| regions[(k + 1) % 3] == regions[(k + 2) % 3])
            .unwrap();
        let a = face[k];
        let area = e.scale.area_to_unit_area(shape.geometry.face_areas()[f]);
        // The fraction of each point from `a`, and its derivative with respect to `t`.
        let from_a = |point: u32| -> (f64, f64) {
            let p = cut.points[point as usize];
            if mesh.edges()[p.edge as usize][0] == a {
                (p.t, 1.0)
            } else {
                (1.0 - p.t, -1.0)
            }
        };
        let [p, q] = segment.points;
        let ((sp, dp), (sq, dq)) = (from_a(p), from_a(q));
        // The corner triangle belongs to the region of `a`.
        let other = regions[(k + 1) % 3];
        let sign = if regions[k] < other { 1.0 } else { -1.0 };
        rates[p as usize] += sign * area * sq * dp * dt[p as usize];
        rates[q as usize] += sign * area * sp * dq * dt[q as usize];
    }
    rates
}

/// The outcome of one finite-difference check.
#[derive(Debug, Clone, Copy)]
pub struct Gradient {
    /// `(E(φ + εψ) − E(φ − εψ)) / 2ε`.
    pub finite_difference: f64,
    /// `Σ_c s_c · dA_c/dε`: the change predicted by the speed of the evaluation.
    pub predicted: f64,
    /// `Σ_c |s_c · dA_c/dε|`: the scale against which the difference is measured, so that
    /// cancellation between parts of the cut does not inflate the error.
    pub scale: f64,
    /// Number of cut points.
    pub n_cut_points: usize,
}

impl Gradient {
    /// `|finite difference − predicted| / Σ |s dA|`.
    pub fn error(&self) -> f64 {
        (self.finite_difference - self.predicted).abs() / self.scale
    }

    /// `|finite difference − predicted| / |finite difference|`.
    pub fn relative_error(&self) -> f64 {
        (self.finite_difference - self.predicted).abs() / self.finite_difference.abs()
    }
}

/// Compares the central finite difference of the energy in the direction `psi` with the change
/// predicted by the speed on the cut: with `dA_c` the area that region 0 gains at cut point `c`,
/// `dE = Σ_c (G_0 − G_1)_c dA_c = Σ_c s_c dA_c` (SPEC §1; for a signed distance function
/// `dA_c ≈ −ε ψ_c ℓ_c`, which gives the form of SPEC §6 item 3).
///
/// The length weight and the penalty must be zero: the kernel does not compute their gradient.
pub fn check_gradient(
    shape: &Shape,
    evaluator: &mut Evaluator,
    phi: &[f64],
    psi: &[f64],
    problem: &CutProblem,
    options: &EvalOptions<'_>,
) -> Gradient {
    assert!(problem.length_weight == 0.0 && problem.penalty.is_none());
    let e = evaluator.evaluate(phi, problem, options).unwrap();
    assert_eq!(
        e.report.n_clamped_crossings, 0,
        "a clamped crossing does not move"
    );
    let rates = swept_area_rates(shape, &e, phi, psi);
    let predicted: f64 = e.speed.iter().zip(&rates).map(|(s, r)| s * r).sum();
    let scale: f64 = e.speed.iter().zip(&rates).map(|(s, r)| (s * r).abs()).sum();
    // A step that moves no crossing by more than 1e-6 of its edge.
    let fastest = e
        .cut
        .points
        .iter()
        .map(|point| {
            let [a, b] = shape.mesh.edge_vertices(point.edge as usize);
            (((1.0 - point.t) * psi[a] + point.t * psi[b]) / (phi[a] - phi[b])).abs()
        })
        .fold(0.0, f64::max);
    let eps = 1e-6 / fastest;
    assert!(
        e.report.smallest_crossing_offset > 1e-4,
        "{}: a crossing is too close to a vertex ({})",
        shape.name,
        e.report.smallest_crossing_offset
    );
    let mut shifted = |sign: f64| -> f64 {
        let moved: Vec<f64> = phi
            .iter()
            .zip(psi)
            .map(|(x, d)| x + sign * eps * d)
            .collect();
        let m = evaluator.evaluate(&moved, problem, options).unwrap();
        assert_eq!(m.cut.vertex_region, e.cut.vertex_region);
        assert_eq!(m.report.n_clamped_crossings, 0);
        m.energy.total
    };
    Gradient {
        finite_difference: (shifted(1.0) - shifted(-1.0)) / (2.0 * eps),
        predicted,
        scale,
        n_cut_points: e.cut.points.len(),
    }
}
