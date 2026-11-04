//! Torus SDF and blend: finite differences (SPEC.md §15.2), singular sets (D12), convexity, the
//! identical-tori blend, the WP3 resolutions (min-ρ shift for constant λ, `DegenerateQuery`,
//! Laplacian 0 at singular points).

mod common;

use common::{ellipsoid_cloud, lcg};
use ddg_points_as_tori::blend::{BlendWant, GradientMode, LambdaRule};
use ddg_points_as_tori::sdf::{Singular, torus_sdf, torus_sdf_gradient, torus_sdf_laplacian};
use ddg_points_as_tori::{EvalOptions, EvalOut, Neighborhood, NormalPolicy, PatError, PointCloud, Tori, Torus, TorusField};

fn unit(v: [f64; 3]) -> [f64; 3] {
    let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / l, v[1] / l, v[2] / l]
}

fn random_torus(s: &mut u64) -> Torus {
    let axis = unit([lcg(s) - 0.5, lcg(s) - 0.5, lcg(s) - 0.5]);
    let major = 0.05 + lcg(s);
    let minor = (0.02 + 1.5 * major * lcg(s)) * if lcg(s) < 0.5 { -1.0 } else { 1.0 }; // includes spindles
    Torus { center: [lcg(s) - 0.5, lcg(s) - 0.5, lcg(s) - 0.5], axis, major, minor }
}

fn shifted(x: [f64; 3], k: usize, h: f64) -> [f64; 3] {
    let mut y = x;
    y[k] += h;
    y
}

#[test]
fn torus_gradient_and_laplacian_match_differences() {
    let mut s = 21;
    let mut n = 0;
    while n < 400 {
        let t = random_torus(&mut s);
        let x = [t.center[0] + 2.0 * lcg(&mut s) - 1.0, t.center[1] + 2.0 * lcg(&mut s) - 1.0, t.center[2] + 2.0 * lcg(&mut s) - 1.0];
        let (g, gr, lap, sing) = torus_sdf_laplacian(&t, x);
        assert_eq!(sing, Singular::None);
        // stay away from the singular sets, where g is not smooth
        let rv = [x[0] - t.center[0], x[1] - t.center[1], x[2] - t.center[2]];
        let z = rv[0] * t.axis[0] + rv[1] * t.axis[1] + rv[2] * t.axis[2];
        let v = ((rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]) - z * z).max(0.0).sqrt();
        let nd = ((v - t.major).powi(2) + z * z).sqrt();
        if v < 0.05 || nd < 0.05 {
            continue;
        }
        n += 1;
        assert_eq!(g, torus_sdf(&t, x));
        let h = 1e-6 * (1.0f64).max((rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]).sqrt());
        let mut l2 = 0.0;
        for k in 0..3 {
            let (gp, gm) = (torus_sdf(&t, shifted(x, k, h)), torus_sdf(&t, shifted(x, k, -h)));
            let fd = (gp - gm) / (2.0 * h);
            assert!((fd - gr[k]).abs() < 1e-7 * (1.0 + gr[k].abs()), "∇g[{k}] {} vs FD {fd}", gr[k]);
            let h2 = 1e-4;
            l2 += (torus_sdf(&t, shifted(x, k, h2)) - 2.0 * g + torus_sdf(&t, shifted(x, k, -h2))) / (h2 * h2);
        }
        assert!((l2 - lap).abs() < 1e-5 * (1.0 + lap.abs()), "Δg {lap} vs FD {l2}");
    }
}

#[test]
fn singular_sets_follow_the_fixed_rules() {
    let t = Torus { center: [0.0; 3], axis: [0.0, 0.0, 1.0], major: 0.5, minor: 0.2 };
    // on the axis: ∇g = (z/‖d‖) u, Δg = 0
    let (g, gr, lap, sing) = torus_sdf_laplacian(&t, [0.0, 0.0, 0.3]);
    assert_eq!(sing, Singular::Axis);
    let nd = (0.25f64 + 0.09).sqrt();
    assert!((g - (nd - 0.2)).abs() < 1e-15);
    assert_eq!(gr[..2], [0.0, 0.0]);
    assert!((gr[2] - 0.3 / nd).abs() < 1e-15);
    assert_eq!(lap, 0.0);
    // on the core circle: ∇g = 0, Δg = 0
    let (g, gr, lap, sing) = torus_sdf_laplacian(&t, [0.5, 0.0, 0.0]);
    assert_eq!((sing, gr, lap), (Singular::Core, [0.0; 3], 0.0));
    assert!((g + 0.2).abs() < 1e-15);
    let (_, _, sing) = torus_sdf_gradient(&t, [0.0, 0.0, 0.0]);
    assert_eq!(sing, Singular::Axis);
}

fn field(cloud: &PointCloud, tori: Vec<Torus>, opts: EvalOptions) -> TorusField {
    let n = tori.len();
    TorusField::new(cloud, Tori { tori, coefficients: vec![[0.0; 6]; n], sigma: vec![f64::NAN; n] }, opts).unwrap()
}

fn cloud_and_tori(n: usize) -> (PointCloud, Vec<Torus>) {
    let cloud = ellipsoid_cloud(n, 3, [0.5, 0.4, 0.3]);
    // tangent spheres of radius 0.4 (inside) as stand-in tori, plus noise in the radii
    let mut s = 8;
    let tori = cloud
        .positions()
        .iter()
        .zip(cloud.normals())
        .map(|(p, nr)| {
            let r = 0.3 + 0.2 * lcg(&mut s);
            Torus {
                center: [p[0] - r * nr[0], p[1] - r * nr[1], p[2] - r * nr[2]],
                axis: unit([lcg(&mut s) - 0.5, lcg(&mut s) - 0.5, 1.0]),
                major: 0.05 * lcg(&mut s),
                minor: r,
            }
        })
        .collect();
    (cloud, tori)
}

#[test]
fn blend_gradients_match_differences() {
    let (cloud, tori) = cloud_and_tori(200);
    let cases = [
        // constant λ: frozen = exact = FD
        (LambdaRule::Constant { lambda: 40.0 }, GradientMode::Frozen),
        (LambdaRule::Constant { lambda: 40.0 }, GradientMode::Exact),
        // per-query λ over all points: only the exact gradient is the derivative
        (LambdaRule::PerQuery { c: 64.0 }, GradientMode::Exact),
    ];
    let mut s = 31;
    let qs: Vec<[f64; 3]> = (0..60).map(|_| [lcg(&mut s) - 0.5, lcg(&mut s) - 0.5, lcg(&mut s) - 0.5]).collect();
    for (rule, mode) in cases {
        let f = field(&cloud, tori.clone(), EvalOptions { neighborhood: Neighborhood::All, lambda: rule, gradient: mode, ..Default::default() });
        for &x in &qs {
            let (_, gr) = f.distance_gradient(x).unwrap();
            let h = 1e-6;
            for k in 0..3 {
                let fd = (f.distance(shifted(x, k, h)).unwrap() - f.distance(shifted(x, k, -h)).unwrap()) / (2.0 * h);
                assert!((fd - gr[k]).abs() < 1e-6 * (1.0 + gr[k].abs()), "{rule:?} {mode:?}: ∇φ[{k}] {} vs FD {fd}", gr[k]);
            }
        }
    }
    // with per-query λ the frozen gradient is not the derivative (the reference's App. D)
    let f = field(&cloud, tori.clone(), EvalOptions { neighborhood: Neighborhood::All, ..Default::default() });
    let x = qs[0];
    let (_, gr) = f.distance_gradient(x).unwrap();
    let fd = (f.distance(shifted(x, 0, 1e-6)).unwrap() - f.distance(shifted(x, 0, -1e-6)).unwrap()) / 2e-6;
    eprintln!("per-query λ, frozen ∇φ_x {} vs FD {fd}", gr[0]);
    // Laplacian against second differences (constant λ, all points)
    let f = field(&cloud, tori, EvalOptions { neighborhood: Neighborhood::All, lambda: LambdaRule::Constant { lambda: 40.0 }, ..Default::default() });
    for &x in &qs[..20] {
        let lap = f.laplacian(x).unwrap();
        // second differences with Richardson extrapolation (λ = 40 makes φ'''' ~ λ⁴ large)
        let phi = f.distance(x).unwrap();
        let d2 = |h: f64| -> f64 {
            (0..3).map(|k| (f.distance(shifted(x, k, h)).unwrap() - 2.0 * phi + f.distance(shifted(x, k, -h)).unwrap()) / (h * h)).sum()
        };
        let fd = (4.0 * d2(5e-5) - d2(1e-4)) / 3.0;
        assert!((fd - lap).abs() < 1e-5 * (1.0 + lap.abs()), "Δφ {lap} vs FD {fd}");
    }
}

#[test]
fn blend_is_convex_and_identical_tori_give_the_torus() {
    let (cloud, tori) = cloud_and_tori(150);
    let f = field(&cloud, tori.clone(), EvalOptions::default());
    let mut s = 41;
    for _ in 0..200 {
        let x = [lcg(&mut s) - 0.5, lcg(&mut s) - 0.5, lcg(&mut s) - 0.5];
        let (v, nb) = f.blend_at(x, BlendWant::Distance).unwrap();
        let gs: Vec<f64> = nb.iter().map(|&i| torus_sdf(&tori[i as usize], x)).collect();
        let (lo, hi) = gs.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &g| (a.0.min(g), a.1.max(g)));
        assert!(v.distance >= lo - 1e-12 && v.distance <= hi + 1e-12);
        assert_eq!(nb.len(), 32);
    }
    let one = tori[0];
    let f = field(&cloud, vec![one; cloud.len()], EvalOptions::default());
    for _ in 0..50 {
        let x = [lcg(&mut s) - 0.5, lcg(&mut s) - 0.5, lcg(&mut s) - 0.5];
        let (d, g) = f.distance_gradient(x).unwrap();
        let (g0, gr0, _) = torus_sdf_gradient(&one, x);
        assert!((d - g0).abs() < 1e-13);
        // frozen gradient: A + φB = ∇g (the λ r̂ terms cancel)
        for k in 0..3 {
            assert!((g[k] - gr0[k]).abs() < 1e-11, "{g:?} vs {gr0:?}");
        }
    }
}

#[test]
fn degenerate_queries_and_overflow() {
    // two points; a radius set holding only the point at the query gives σ_x = 0
    let cloud = PointCloud::new(vec![[0.0; 3], [1.0, 0.0, 0.0]], vec![[0.0, 0.0, 1.0]; 2], NormalPolicy::default()).unwrap().0;
    let t = Torus { center: [0.0, 0.0, -1.0], axis: [0.0, 0.0, 1.0], major: 0.0, minor: 1.0 };
    let f = field(&cloud, vec![t; 2], EvalOptions { neighborhood: Neighborhood::Radius { r: 0.1, fallback_k: 1 }, ..Default::default() });
    assert!(matches!(f.distance([0.0; 3]), Err(PatError::DegenerateQuery { index: 0, .. })));
    let mut d = vec![0.0; 3];
    let r = f.evaluate(&[[0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [1.0, 0.0, 0.0]], EvalOut { distance: Some(&mut d), ..Default::default() });
    assert!(matches!(r, Err(PatError::DegenerateQuery { index: 2, .. })), "{r:?}");
    // constant λ is fine there, and far queries do not overflow (shift = min ρ)
    let f = field(
        &cloud,
        vec![t; 2],
        EvalOptions { neighborhood: Neighborhood::Radius { r: 0.1, fallback_k: 2 }, lambda: LambdaRule::Constant { lambda: 1e4 }, ..Default::default() },
    );
    assert!(f.distance([0.0; 3]).unwrap().abs() < 1e-15);
    let far = f.distance_gradient([30.0, 40.0, 0.0]).unwrap();
    assert!(far.0.is_finite() && far.1.iter().all(|x| x.is_finite()), "{far:?}");
    let (v, _) = f.blend_at([30.0, 40.0, 0.0], BlendWant::Distance).unwrap();
    // the shift is min ρ (the nearest member has weight 1), not σ_x = ½ max ρ
    assert!((v.shift - (30.0f64.powi(2) + 40.0f64.powi(2)).sqrt().min((29.0f64.powi(2) + 1600.0).sqrt())).abs() < 1e-12);
    assert!(v.shift > v.sigma);
    // a query at a data point: r̂ := 0, counted
    let f = field(&cloud, vec![t; 2], EvalOptions { neighborhood: Neighborhood::All, ..Default::default() });
    let mut g = vec![[0.0; 3]; 1];
    let rep = f.evaluate(&[[0.0; 3]], EvalOut { distance: Some(&mut d[..1]), gradient: Some(&mut g), laplacian: None }).unwrap();
    assert_eq!(rep.at_data_point, 1);
    assert!(g[0].iter().all(|x| x.is_finite()));
    assert!(matches!(f.distance([f64::NAN, 0.0, 0.0]), Err(PatError::NonFinite { .. })));
}
