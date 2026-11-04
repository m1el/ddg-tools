//! Torus fitting on analytic cases (SPEC.md §15.1-1, §15.1-3): spheres, tori at their outer and
//! inner equators, planes, cylinders; clamps and errors; invariance of the features and the
//! pipeline under translation and scaling.

mod common;

use common::{ellipsoid_cloud, lcg};
use ddg_points_as_tori::features::{FeatureOrder, FeatureScratch, local_frame, neighborhood_features};
use ddg_points_as_tori::fit::fit_torus;
use ddg_points_as_tori::sdf::torus_sdf;
use ddg_points_as_tori::spatial::{KdTree, KdTreeOptions, fitting_neighborhoods};
use ddg_points_as_tori::{FitOptions, Model, NormalPolicy, PatError, PointCloud, PredictOptions, SignRule, fit_tori, predict_coefficients};

fn unit(v: [f64; 3]) -> [f64; 3] {
    let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / l, v[1] / l, v[2] / l]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Coefficients of the second-order Monge patch with principal curvatures `k1` along `e1` and
/// `k2` along `n × e1` (normal `n`), in the frame of `n`: `Q = ½ (II_ss s² + 2 II_st s t + II_tt t²)`.
fn monge(n: [f64; 3], e1: [f64; 3], k1: f64, k2: f64, a00: f64) -> [f64; 6] {
    let e2 = cross(n, e1);
    let (s, t) = local_frame(n);
    let ii = |a: [f64; 3], b: [f64; 3]| k1 * dot(a, e1) * dot(b, e1) + k2 * dot(a, e2) * dot(b, e2);
    // (a00, a01, a10, a11, a02, a20)
    [a00, 0.0, 0.0, ii(s, t), 0.5 * ii(t, t), 0.5 * ii(s, s)]
}

fn random_unit(s: &mut u64) -> [f64; 3] {
    loop {
        let v = [2.0 * lcg(s) - 1.0, 2.0 * lcg(s) - 1.0, 2.0 * lcg(s) - 1.0];
        let l = dot(v, v).sqrt();
        if (0.1..1.0).contains(&l) {
            return unit(v);
        }
    }
}

fn close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() <= tol)
}

#[test]
fn sphere_gives_the_tangent_sphere() {
    let mut s = 1;
    for _ in 0..50 {
        let n = random_unit(&mut s);
        let p = [lcg(&mut s), lcg(&mut s), lcg(&mut s)];
        let rho = 0.05 + lcg(&mut s);
        let a00 = 0.01 * (lcg(&mut s) - 0.5);
        let a = monge(n, unit(cross(n, [0.3, 0.2, 0.1])), -1.0 / rho, -1.0 / rho, a00);
        let (t, d) = fit_torus(p, n, a, &FitOptions::default());
        let want_c = [p[0] + (a00 - rho) * n[0], p[1] + (a00 - rho) * n[1], p[2] + (a00 - rho) * n[2]];
        assert!(t.major.abs() < 1e-9, "R = {}", t.major);
        assert!((t.minor - rho).abs() < 1e-9 * rho.max(1.0), "r = {} vs {rho}", t.minor);
        assert!(close(t.center, want_c, 1e-9), "{:?} vs {want_c:?}", t.center);
        assert!((d.kappa_min + 1.0 / rho).abs() < 1e-9 / rho);
        // the SDF is that of the sphere
        for _ in 0..10 {
            let x = [p[0] + lcg(&mut s) - 0.5, p[1] + lcg(&mut s) - 0.5, p[2] + lcg(&mut s) - 0.5];
            let r = ((x[0] - want_c[0]).powi(2) + (x[1] - want_c[1]).powi(2) + (x[2] - want_c[2]).powi(2)).sqrt();
            assert!((torus_sdf(&t, x) - (r - rho)).abs() < 1e-9);
        }
    }
}

#[test]
fn torus_patch_gives_the_torus_back() {
    let mut s = 2;
    for _ in 0..50 {
        let axis = random_unit(&mut s);
        let c = [lcg(&mut s), lcg(&mut s), lcg(&mut s)];
        let r0 = 0.05 + 0.2 * lcg(&mut s);
        let big = 2.2 * r0 + lcg(&mut s); // R > 2r so that the inner equator has |κ_parallel| < |κ_meridian|
        let radial = unit(cross(axis, random_unit(&mut s)));
        let tangent = cross(axis, radial);
        for outer in [true, false] {
            // outer equator: p = c + (R + r) e, n = e; inner: p = c + (R − r) e, n = −e
            let (p, n, k_par) = if outer {
                (std::array::from_fn(|i| c[i] + (big + r0) * radial[i]), radial, -1.0 / (big + r0))
            } else {
                (std::array::from_fn(|i| c[i] + (big - r0) * radial[i]), radial.map(|x| -x), 1.0 / (big - r0))
            };
            let a = monge(n, tangent, k_par, -1.0 / r0, 0.0);
            let geo = FitOptions { sign_rule: SignRule::Geometric, ..Default::default() };
            let (t, d) = fit_torus(p, n, a, &geo);
            let (tr, _) = fit_torus(p, n, a, &FitOptions::default());
            if outer {
                assert_eq!(t, tr, "the rules agree away from saddles");
            } else {
                // the reference's saddle torus is the mirror image: center on the far side, and
                // the SDF decreases along the outward normal at p
                let mirrored: [f64; 3] = std::array::from_fn(|i| c[i] + 2.0 * (big - r0) * radial[i]);
                assert!(close(tr.center, mirrored, 1e-9 * (1.0 + big)), "{:?} vs {mirrored:?}", tr.center);
                let h = 1e-5;
                let dn = (torus_sdf(&tr, std::array::from_fn(|i| p[i] + h * n[i])) - torus_sdf(&tr, std::array::from_fn(|i| p[i] - h * n[i]))) / (2.0 * h);
                assert!((dn + 1.0).abs() < 1e-6, "reference saddle ∂g/∂n = {dn}");
                let dn = (torus_sdf(&t, std::array::from_fn(|i| p[i] + h * n[i])) - torus_sdf(&t, std::array::from_fn(|i| p[i] - h * n[i]))) / (2.0 * h);
                assert!((dn - 1.0).abs() < 1e-6, "geometric saddle ∂g/∂n = {dn}");
            }
            let tol = 1e-9 * (1.0 + big);
            assert!((t.major - big).abs() < tol, "outer {outer}: R {} vs {big}", t.major);
            assert!((t.minor - r0).abs() < tol, "outer {outer}: r {} vs {r0}", t.minor);
            assert!(close(t.center, c, tol), "outer {outer}: center {:?} vs {c:?}", t.center);
            assert!((dot(t.axis, axis).abs() - 1.0).abs() < 1e-9, "axis {:?} vs {axis:?}", t.axis);
            assert!(dot(d.v_min, tangent).abs() > 1.0 - 1e-9);
            // contact: the torus passes through Q*(0,0) = p
            assert!(torus_sdf(&t, p).abs() < 1e-9);
        }
    }
}

#[test]
fn plane_gives_a_half_space_and_counts_clamps() {
    let cloud = PointCloud::new(vec![[0.1, 0.2, 0.3], [1.0, 0.0, 0.0]], vec![unit([0.2, -0.3, 1.0]), [1.0, 0.0, 0.0]], NormalPolicy::default())
        .unwrap()
        .0;
    let co = vec![[0.0f32; 6], [0.0, 0.0, 0.0, 0.0, -0.5, -0.5]];
    let (p, n) = (cloud.positions()[0], cloud.normals()[0]);
    // the JAX rule: R = r = 1e10, s = −1, c = p + 1e10 n*: the core circle passes through p
    let (tori, rep) = fit_tori(&cloud, &co, &FitOptions::default()).unwrap();
    assert_eq!((rep.clamped_curvature, rep.zero_curvature), (1, 1), "{rep:?}");
    let t = tori.tori[0];
    assert!((t.major - 1e10).abs() < 1.0 && (t.minor + 1e10).abs() < 1.0, "{t:?}");
    assert!((torus_sdf(&t, p) - 1e10).abs() < 1.0, "reference rule: g(p) = {}", torus_sdf(&t, p));
    // ours: a sphere of radius 1e10 tangent at p, i.e. the half-space
    let opts = FitOptions { sign_rule: SignRule::Geometric, ..Default::default() };
    let (tori, rep) = fit_tori(&cloud, &co, &opts).unwrap();
    assert_eq!(rep.clamped_curvature, 1, "{rep:?}");
    assert_eq!((rep.clamped_disc, rep.umbilic), (2, 2), "{rep:?}"); // the plane and the umbilic sphere
    let t = tori.tori[0];
    assert!(t.major.abs() < 1.0 && (t.minor + 1e10).abs() < 1.0, "{t:?}");
    let mut s = 3;
    for _ in 0..100 {
        let x = [p[0] + 0.2 * (lcg(&mut s) - 0.5), p[1] + 0.2 * (lcg(&mut s) - 0.5), p[2] + 0.2 * (lcg(&mut s) - 0.5)];
        let h = dot([x[0] - p[0], x[1] - p[1], x[2] - p[2]], n);
        assert!((torus_sdf(&t, x) - h).abs() < 1e-5, "{} vs {h}", torus_sdf(&t, x));
    }
}

#[test]
fn cylinder_gives_a_tube() {
    let mut s = 4;
    for _ in 0..20 {
        let n = random_unit(&mut s);
        let generator = unit(cross(n, random_unit(&mut s)));
        let rho = 0.1 + lcg(&mut s);
        let p = [lcg(&mut s), lcg(&mut s), lcg(&mut s)];
        // curvature 0 along the generator, −1/ρ across
        let a = monge(n, generator, 0.0, -1.0 / rho, 0.0);
        let a = a.map(|x| x as f32 as f64);
        let opts = FitOptions { sign_rule: SignRule::Geometric, ..Default::default() };
        let (t, d) = fit_torus(p, n, a, &opts);
        assert!(d.kappa_min.abs() < 1e-6);
        if d.kappa_min * d.kappa_max > 0.0 {
            // tiny but non-zero, not a saddle: both rules agree exactly (with κmin > 0 from
            // round-off the reference puts the torus on the wrong side, see SignRule)
            assert_eq!(fit_torus(p, n, a, &FitOptions::default()).0, t);
        }
        assert!((t.minor.abs() - rho).abs() < 1e-6 * rho);
        // near p the SDF is the distance to the cylinder of radius ρ about the line p − ρ n + s·generator
        let axis_p = [p[0] - rho * n[0], p[1] - rho * n[1], p[2] - rho * n[2]];
        for _ in 0..20 {
            let x = [p[0] + 0.1 * (lcg(&mut s) - 0.5), p[1] + 0.1 * (lcg(&mut s) - 0.5), p[2] + 0.1 * (lcg(&mut s) - 0.5)];
            let w = [x[0] - axis_p[0], x[1] - axis_p[1], x[2] - axis_p[2]];
            let along = dot(w, generator);
            let perp = (dot(w, w) - along * along).sqrt();
            assert!((torus_sdf(&t, x) - (perp - rho)).abs() < 1e-4, "{} vs {}", torus_sdf(&t, x), perp - rho);
        }
    }
}

#[test]
fn bad_inputs_are_errors() {
    let cloud = ellipsoid_cloud(10, 1, [1.0, 1.0, 1.0]);
    assert!(matches!(fit_tori(&cloud, &[[0.0; 6]; 9], &FitOptions::default()), Err(PatError::SizeMismatch { .. })));
    let mut co = vec![[0.0f32, 0.0, 0.0, 0.0, -0.5, -0.5]; 10];
    co[7][3] = f32::NAN;
    assert!(matches!(fit_tori(&cloud, &co, &FitOptions::default()), Err(PatError::NonFinite { index: 7, .. })));
    let bad = FitOptions { eig_eps: 0.0, ..Default::default() };
    assert!(matches!(fit_tori(&cloud, &co, &bad), Err(PatError::InvalidOption { .. })));
    // cloud validation
    assert!(matches!(PointCloud::new(vec![[0.0; 3]], vec![[1.0, 0.0, 0.0]], NormalPolicy::default()), Err(PatError::TooFewPoints { n: 1, need: 2 })));
    assert!(matches!(
        PointCloud::new(vec![[0.0; 3]; 3], vec![[1.0, 0.0, 0.0]; 2], NormalPolicy::default()),
        Err(PatError::LengthMismatch { positions: 3, normals: 2 })
    ));
    assert!(matches!(
        PointCloud::new(vec![[0.0; 3]; 3], vec![[1.0, 0.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0]], NormalPolicy::Normalize),
        Err(PatError::ZeroNormal { index: 1 })
    ));
    assert!(matches!(
        PointCloud::new(vec![[0.0; 3]; 3], vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]], NormalPolicy::default()),
        Err(PatError::NonUnitNormal { index: 1, .. })
    ));
    assert!(matches!(
        PointCloud::new(vec![[0.0; 3], [f64::INFINITY, 0.0, 0.0]], vec![[1.0, 0.0, 0.0]; 2], NormalPolicy::default()),
        Err(PatError::NonFinite { index: 1, what: "position" })
    ));
    let (c, rep) = PointCloud::new(vec![[0.0; 3], [0.0; 3], [1.0, 0.0, 0.0]], vec![[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 3.0]], NormalPolicy::Normalize).unwrap();
    assert_eq!((rep.normalized, rep.exact_duplicates), (2, 1));
    assert_eq!(c.normals()[0], [1.0, 0.0, 0.0]);
}

#[test]
fn features_and_coefficients_are_translation_and_scale_covariant() {
    let model = Model::reference().unwrap();
    let cloud = ellipsoid_cloud(400, 9, [0.5, 0.4, 0.3]);
    let (p, n) = (cloud.positions(), cloud.normals());
    let tree = KdTree::new(p, KdTreeOptions::default()).unwrap();
    let nb = fitting_neighborhoods(&tree, p, 64);
    let (c, shift) = (3.0, [10.0, -4.0, 2.5]);
    let moved: Vec<[f64; 3]> = p.iter().map(|q| std::array::from_fn(|i| c * q[i] + shift[i])).collect();
    let cloud2 = PointCloud::new(moved.clone(), n.to_vec(), NormalPolicy::default()).unwrap().0;
    // features (difference-first): equal to round-off
    let mut s = FeatureScratch::default();
    let (mut f1, mut f2) = (vec![0f32; 64 * 6], vec![0f32; 64 * 6]);
    for i in (0..400).step_by(37) {
        let s1 = neighborhood_features(p, n, nb.row(i), FeatureOrder::DifferenceFirst, &mut s, &mut f1);
        let s2 = neighborhood_features(&moved, n, nb.row(i), FeatureOrder::DifferenceFirst, &mut s, &mut f2);
        assert!((s2 / s1 - c).abs() < 1e-12);
        let d = f1.iter().zip(&f2).map(|(a, b)| (a - b).abs()).fold(0f32, f32::max);
        assert!(d <= 2.0 * f32::EPSILON, "features differ by {d}");
    }
    // the neighborhoods are the same sets, so the coefficients scale: a00·c, slopes, curvatures/c
    let tree2 = KdTree::new(&moved, KdTreeOptions::default()).unwrap();
    let nb2 = fitting_neighborhoods(&tree2, &moved, 64);
    let (a1, _, _) = predict_coefficients(&model, &cloud, &nb, &PredictOptions::default()).unwrap();
    let (a2, _, _) = predict_coefficients(&model, &cloud2, &nb2, &PredictOptions::default()).unwrap();
    let scale = [c, 1.0, 1.0, 1.0 / c, 1.0 / c, 1.0 / c];
    let mut worst: f64 = 0.0;
    for (x, y) in a1.iter().zip(&a2) {
        for k in 0..6 {
            worst = worst.max(((y[k] as f64) - scale[k] * x[k] as f64).abs() / (1.0 + (y[k] as f64).abs()));
        }
    }
    assert!(worst < 1e-5, "coefficients not covariant: {worst:e}");
}
