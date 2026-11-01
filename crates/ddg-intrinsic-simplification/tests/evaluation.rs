//! The distortion measure and the local unfolding distance (SPEC §10.3), and the extrinsic maps
//! as applications (§10.4).

use ddg_core::TriMesh;
use ddg_core::generators::{grid_disk, icosphere};
use ddg_core::sparse::ThreadPool;
use ddg_intrinsic_simplification::applications::CoarseApplications;
use ddg_intrinsic_simplification::baseline::{QemOptions, qem_simplify};
use ddg_intrinsic_simplification::evaluation::{
    DistortionOptions, ExtrinsicMap, PointDistance, coarse_area, distortion, extrinsic_applications,
};
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions};

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
    fn bary(&mut self) -> [f64; 3] {
        let (a, b) = (self.next(), self.next());
        let (a, b) = if a + b > 1.0 {
            (1.0 - a, 1.0 - b)
        } else {
            (a, b)
        };
        [1.0 - a - b, a, b]
    }
}

fn pool() -> ThreadPool {
    ThreadPool::new(2).unwrap()
}

#[test]
fn unfolding_distance_on_a_plane_is_euclidean() {
    let g = grid_disk(8, 8).unwrap();
    assert!(g.positions.iter().all(|p| p[2] == 0.0));
    let mesh =
        CoarseMesh::from_positions(&g.positions, &g.faces, &[], CoarseMeshOptions::default())
            .unwrap();
    let app = CoarseApplications::new(&mesh).unwrap();
    let out = app.output();
    let deep = DistortionOptions {
        max_depth: 100,
        max_nodes: 100_000,
        ..DistortionOptions::default()
    };
    let pd = PointDistance::new(out, &deep).unwrap();
    let pos = |f: u32, b: [f64; 3]| {
        let c = out.faces[f as usize].map(|v| g.positions[out.vertices[v as usize] as usize]);
        [0, 1].map(|k| b[0] * c[0][k] + b[1] * c[1][k] + b[2] * c[2][k])
    };
    let mut rng = Lcg(7);
    let nf = out.faces.len() as f64;
    let mut worst: f64 = 0.0;
    for _ in 0..300 {
        let (fa, fb) = ((rng.next() * nf) as u32, (rng.next() * nf) as u32);
        let (ba, bb) = (rng.bary(), rng.bary());
        let (pa, pb) = (pos(fa, ba), pos(fb, bb));
        let expect = ((pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2)).sqrt();
        let d = pd.distance((fa, ba), (fb, bb)).expect("a path");
        worst = worst.max((d - expect).abs());
    }
    assert!(worst < 1e-7, "worst error {worst}");
}

/// A cone: an apex and a ring of `n` vertices, `n` faces of apex angle `gamma`, slant length one.
/// Points in polar coordinates of the development have the closed-form distance
/// `√(r₁² + r₂² − 2 r₁ r₂ cos δ)` with `δ = min(β, Θ − β)`, or `r₁ + r₂` when `δ ≥ π`.
fn cone_check(n: usize, gamma: f64, expect_through_apex: bool) {
    let faces: Vec<[u32; 3]> = (0..n)
        .map(|k| [0, 1 + k as u32, 1 + ((k + 1) % n) as u32])
        .collect();
    let tri = TriMesh::new(n + 1, faces).unwrap();
    let chord = 2.0 * (gamma / 2.0).sin();
    let lengths: Vec<f64> = tri
        .edges()
        .iter()
        .map(|e| if e[0] == 0 || e[1] == 0 { 1.0 } else { chord })
        .collect();
    let mesh = CoarseMesh::from_lengths(&tri, lengths, &[], CoarseMeshOptions::default()).unwrap();
    let app = CoarseApplications::new(&mesh).unwrap();
    let out = app.output();
    let apex = out.coarse_of_fine[0].unwrap();
    let theta_total = n as f64 * gamma;
    let ring_of = |v: u32| out.vertices[v as usize] as usize - 1;
    // The development of a point: polar coordinates (r, θ) with θ ∈ [0, Θ). Face k is spanned
    // by the ring vertices k and k + 1 (mod n), at angles kγ and (k + 1)γ.
    let dev = |f: u32, b: [f64; 3]| {
        let c = out.faces[f as usize];
        let mut ring: Vec<(usize, f64)> = (0..3)
            .filter(|&s| c[s] != apex)
            .map(|s| (ring_of(c[s]), b[s]))
            .collect();
        ring.sort_by_key(|x| x.0);
        // The face (n - 1, 0) has ring indices {0, n - 1}: its first corner is n - 1.
        let (k, w0, w1) = if ring[0].0 == 0 && ring[1].0 == n - 1 {
            (n - 1, ring[1].1, ring[0].1)
        } else {
            (ring[0].0, ring[0].1, ring[1].1)
        };
        let x = w0 + w1 * gamma.cos();
        let y = w1 * gamma.sin();
        ((x * x + y * y).sqrt(), k as f64 * gamma + y.atan2(x))
    };
    let pd = PointDistance::new(out, &DistortionOptions::default()).unwrap();
    let mut rng = Lcg(11);
    let nf = out.faces.len() as f64;
    let (mut worst, mut through) = (0.0f64, 0usize);
    for _ in 0..400 {
        let (fa, fb) = ((rng.next() * nf) as u32, (rng.next() * nf) as u32);
        let (ba, bb) = (rng.bary(), rng.bary());
        let ((ra, ta), (rb, tb)) = (dev(fa, ba), dev(fb, bb));
        let beta = (ta - tb).abs();
        let delta = beta.min(theta_total - beta);
        let expect = if delta >= std::f64::consts::PI {
            through += 1;
            ra + rb
        } else {
            (ra * ra + rb * rb - 2.0 * ra * rb * delta.cos()).sqrt()
        };
        let d = pd.distance((fa, ba), (fb, bb)).expect("a path");
        worst = worst.max((d - expect).abs());
    }
    if expect_through_apex {
        assert!(through > 10, "only {through} paths through the apex");
    }
    assert!(
        worst < 1e-6,
        "worst error {worst} (through apex: {through})"
    );
}

#[test]
fn unfolding_distance_on_cones_wraps_and_bends_at_the_apex() {
    // Cone angle below 2π: the short way round can be the other side of the apex.
    cone_check(12, 0.5, false);
    // Cone angle above 2π (negative curvature): shortest paths bend at the apex.
    cone_check(12, 0.55, true);
}

fn flat_pair(sx: f64, sy: f64) -> (Vec<[u32; 3]>, Vec<[f64; 3]>, CoarseApplications) {
    let g = grid_disk(6, 6).unwrap();
    let scaled: Vec<[f64; 3]> = g
        .positions
        .iter()
        .map(|p| [p[0] * sx, p[1] * sy, 0.0])
        .collect();
    let mesh =
        CoarseMesh::from_positions(&scaled, &g.faces, &[], CoarseMeshOptions::default()).unwrap();
    let app = CoarseApplications::new(&mesh).unwrap();
    (g.faces, g.positions, app)
}

#[test]
fn distortion_of_known_maps() {
    // Identity.
    let (faces, pos, app) = flat_pair(1.0, 1.0);
    let r = distortion(&pos, &faces, &app, &pool(), &DistortionOptions::default()).unwrap();
    assert_eq!(r.n_invalid_faces, 0);
    assert!(
        (r.aniso_mean - 1.0).abs() < 1e-6 && (r.aniso_max - 1.0).abs() < 1e-6,
        "{r:?}"
    );
    assert!(r.area_max_abs < 1e-6);
    // Uniform scale 2: isotropic, area distortion log(4).
    let (faces, pos, app) = flat_pair(2.0, 2.0);
    let r = distortion(&pos, &faces, &app, &pool(), &DistortionOptions::default()).unwrap();
    assert!((r.aniso_max - 1.0).abs() < 1e-6, "{r:?}");
    assert!((r.area_mean_signed - 4.0f64.ln()).abs() < 1e-6, "{r:?}");
    // Stretch by (2, 1): σ = (2, 1) on every face: γ_max/γ_min = 4, σ₁/σ₂ = 2, log σ₁σ₂ = ln 2.
    let (faces, pos, app) = flat_pair(2.0, 1.0);
    let r = distortion(&pos, &faces, &app, &pool(), &DistortionOptions::default()).unwrap();
    assert!(
        (r.aniso_mean - 4.0).abs() < 1e-6 && (r.aniso_max - 4.0).abs() < 1e-6,
        "{r:?}"
    );
    assert!((r.sigma_ratio_max - 2.0).abs() < 1e-6, "{r:?}");
    assert!((r.area_mean_abs - 2.0f64.ln()).abs() < 1e-6, "{r:?}");
    assert_eq!(r.per_face.len(), faces.len());
}

#[test]
fn distortion_of_the_qem_maps_runs_and_is_finite() {
    let s = icosphere(3).unwrap();
    let o = QemOptions {
        successive_map: true,
        ..QemOptions::new(80)
    };
    let q = qem_simplify(&s.positions, &s.faces, &o).unwrap();
    for map in [ExtrinsicMap::ClosestPoint, ExtrinsicMap::Successive] {
        let app = extrinsic_applications(&q, s.positions.len(), map).unwrap();
        let area = coarse_area(&app).unwrap();
        assert!((area - q.report.coarse_area).abs() < 1e-9 * area);
        let r = distortion(
            &s.positions,
            &s.faces,
            &app,
            &pool(),
            &DistortionOptions::default(),
        )
        .unwrap();
        assert_eq!(r.n_unresolved, 0, "{r:?}");
        assert!(
            r.aniso_mean.is_finite() && r.aniso_mean >= 1.0,
            "{map:?} {r:?}"
        );
        eprintln!(
            "{map:?}: aniso mean {:.3} max {:.1}, |area| mean {:.3}, invalid {} ({:.3}), edges same/unf/via/fb {}/{}/{}/{}",
            r.aniso_mean,
            r.aniso_max,
            r.area_mean_abs,
            r.n_invalid_faces,
            r.invalid_area_fraction,
            r.n_same_face,
            r.n_unfolded,
            r.n_via_vertex,
            r.n_fallback
        );
    }
}
