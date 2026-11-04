//! GPU back end (feature `cuda`; SPEC.md §5.4, §7.5, §18.2 D4). Every test skips with a message
//! when no CUDA device is usable, unless `DDG_REQUIRE_GPU` is set.

#![cfg(feature = "cuda")]

mod common;

use common::{bits, ellipsoid_cloud, lcg};
use ddg_points_as_tori::features::{FeatureOrder, FeatureScratch, neighborhood_features};
use ddg_points_as_tori::fit::{EigenRule, fit_tori_with};
use ddg_points_as_tori::gpu::{GpuContext, fit_tori as gfit, fitting_neighborhoods as gfit_nb, knn_batch as gknn};
use ddg_points_as_tori::network::GpuMatmul;
use ddg_points_as_tori::network::cuda::forward_features_cuda;
use ddg_points_as_tori::network::{PredictOptions, forward_features};
use ddg_points_as_tori::spatial::{KdTree, KdTreeOptions, fitting_neighborhoods, knn_batch};
use ddg_points_as_tori::{
    BuildOptions, Device, EvalOptions, EvalOut, EvalPrecision, EvalReport, FitOptions, FitReport, GradientMode, LambdaRule, Model,
    Neighborhood, PatError, PointCloud, PointsAsTori, SignRule, Tori, TorusField,
};

fn gpu() -> Option<usize> {
    match ddg_points_as_tori::gpu::probe(0) {
        Ok(info) => {
            eprintln!("GPU: {} (cc {:?}, {} MB)", info.name, info.compute_capability, info.total_memory >> 20);
            Some(0)
        }
        Err(e) => {
            assert!(std::env::var_os("DDG_REQUIRE_GPU").is_none(), "DDG_REQUIRE_GPU is set but no GPU: {e}");
            eprintln!("SKIP: no usable CUDA device ({e})");
            None
        }
    }
}

fn features(cloud: &PointCloud, k: usize) -> Vec<f32> {
    let tree = KdTree::new(cloud.positions(), KdTreeOptions::default()).unwrap();
    let nb = fitting_neighborhoods(&tree, cloud.positions(), k);
    let mut scr = FeatureScratch::default();
    let mut out = vec![0f32; cloud.len() * k * 6];
    for i in 0..cloud.len() {
        neighborhood_features(cloud.positions(), cloud.normals(), nb.row(i), FeatureOrder::DifferenceFirst, &mut scr, &mut out[i * k * 6..(i + 1) * k * 6]);
    }
    out
}

fn hash32(v: &[f32]) -> u64 {
    v.iter().fold(0xcbf29ce484222325u64, |h, x| (h ^ x.to_bits() as u64).wrapping_mul(0x100000001b3))
}

fn rel(a: &[f32], b: &[f32]) -> f32 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs() / (1.0 + x.abs())).fold(0f32, f32::max)
}

#[test]
fn network_run_to_run_and_batch_shapes() {
    let Some(ord) = gpu() else { return };
    let model = Model::reference().unwrap();
    let cloud = ellipsoid_cloud(3000, 7, [0.6, 0.45, 0.3]);
    let f = features(&cloud, 64);
    let cpu = forward_features(&model, &f, 64, &PredictOptions { threads: Some(16), ..Default::default() }).unwrap();
    for mode in [GpuMatmul::Tf32, GpuMatmul::Split3] {
        let base = forward_features_cuda(&model, &f, 64, ord, None, mode).unwrap();
        let again = forward_features_cuda(&model, &f, 64, ord, None, mode).unwrap();
        eprintln!("{mode:?}: hash {:016x} / {:016x}, cpu vs gpu raw max |d|/(1+|a|) {:e}", hash32(&base), hash32(&again), rel(&cpu, &base));
        for b in [1usize, 7, 64, 256, 1000, 4096] {
            let r = forward_features_cuda(&model, &f, 64, ord, Some(b), mode).unwrap();
            let nd = r.iter().zip(&base).filter(|(a, b)| a.to_bits() != b.to_bits()).count();
            eprintln!("  batch {b:>5}: hash {:016x}, differs from default in {nd} values", hash32(&r));
        }
        assert_eq!(hash32(&base), hash32(&again));
    }
}

#[test]
#[ignore]
fn network_throughput() {
    let Some(ord) = gpu() else { return };
    let model = Model::reference().unwrap();
    let n = 40_000;
    let cloud = ellipsoid_cloud(n, 9, [0.6, 0.45, 0.3]);
    let f = features(&cloud, 64);
    for mode in [GpuMatmul::Tf32, GpuMatmul::Split3] {
        for b in [256usize, 512, 1024, 2048] {
            let r = ddg_points_as_tori::network::cuda::forward_features_cuda(&model, &f[..b * 64 * 6], 64, ord, Some(b), mode).unwrap();
            let t = std::time::Instant::now();
            let r2 = forward_features_cuda(&model, &f, 64, ord, Some(b), mode).unwrap();
            let s = t.elapsed().as_secs_f64();
            eprintln!("{mode:?} batch {b}: {n} in {s:.2} s = {:.0}/s ({} {})", n as f64 / s, r.len(), r2.len());
        }
    }
}

// ------------------------------------------------------------------------------------------------
// Our kernels: kNN, fit, blend (level 1 and level 2)
// ------------------------------------------------------------------------------------------------

fn grid_cloud(n: usize, seed: u64) -> Vec<[f64; 3]> {
    let mut s = seed;
    (0..n).map(|_| [lcg(&mut s), lcg(&mut s), lcg(&mut s)].map(|c| (c * 6.0).floor() / 6.0)).collect()
}

fn rand_points(n: usize, seed: u64) -> Vec<[f64; 3]> {
    let mut s = seed;
    (0..n).map(|_| [lcg(&mut s), lcg(&mut s), lcg(&mut s)]).collect()
}

#[test]
fn knn_bitwise_equal_to_cpu() {
    let Some(ord) = gpu() else { return };
    let g = GpuContext::new(ord).unwrap();
    for (n, grid) in [(1usize, false), (17, true), (300, true), (5000, false), (20_000, true)] {
        let pts = if grid { grid_cloud(n, 3 + n as u64) } else { rand_points(n, 3 + n as u64) };
        let tree = KdTree::new(&pts, KdTreeOptions::default()).unwrap();
        let mut qs = rand_points(2000, 99);
        qs.extend(pts.iter().take(500).copied());
        for k in [1usize, 5, 32, 33, 64, 100, 128] {
            let cpu = knn_batch(&tree, &qs, k);
            let gpu = gknn(&g, &tree, &qs, k).unwrap();
            assert!(cpu.0 == gpu.0, "n {n} k {k}: indices differ");
            assert!(bits(&cpu.1) == bits(&gpu.1), "n {n} k {k}: distances differ");
        }
        for k in [1usize, 16, 64, 128] {
            let cpu = fitting_neighborhoods(&tree, &pts, k);
            let gpu = gfit_nb(&g, &tree, &pts, k).unwrap();
            assert_eq!(cpu, gpu, "n {n} k {k}: fitting neighborhoods differ");
        }
    }
    assert!(gknn(&g, &KdTree::new(&[[0.0; 3]], KdTreeOptions::default()).unwrap(), &[[f64::NAN, 0.0, 0.0]], 1).is_err());
    assert!(gknn(&g, &KdTree::new(&[[0.0; 3]], KdTreeOptions::default()).unwrap(), &[[0.0; 3]], 129).is_err());
}

fn random_coeffs(n: usize, seed: u64) -> Vec<[f32; 6]> {
    let mut s = seed;
    (0..n)
        .map(|i| {
            let mut a = [0f32; 6];
            for v in a.iter_mut() {
                *v = (4.0 * lcg(&mut s) - 2.0) as f32;
            }
            // edge cases: umbilics, planes, cylinders, v2 = 0, saddles, tiny curvature
            match i % 9 {
                0 => a = [0.0; 6],
                1 => a = [0.01, 0.0, 0.0, 0.0, 0.5, 0.5],
                2 => a = [0.0, 0.0, 0.0, 0.0, 0.0, 0.7],
                3 => a = [0.0, 0.0, 0.0, 0.0, -0.3, 0.6],
                4 => a = [0.0, 0.1, -0.2, 0.0, 1e-12, 2e-12],
                _ => {}
            }
            a
        })
        .collect()
}

#[test]
fn fit_bitwise_equal_to_cpu() {
    let Some(ord) = gpu() else { return };
    let g = GpuContext::new(ord).unwrap();
    let cloud = ellipsoid_cloud(4000, 21, [0.6, 0.45, 0.3]);
    let co = random_coeffs(cloud.len(), 5);
    let sigma: Vec<f64> = (0..cloud.len()).map(|i| 0.01 + i as f64 * 1e-6).collect();
    for sign_rule in [SignRule::Reference, SignRule::Geometric] {
        for eigen_rule in [EigenRule::Reference, EigenRule::Robust] {
            let o = FitOptions { sign_rule, eigen_rule, ..Default::default() };
            let (tc, rc) = fit_tori_with(&cloud, &co, &sigma, 0, &o).unwrap();
            let (tg, rg) = gfit(&g, &cloud, &co, &sigma, 0, &o).unwrap();
            let tb = |t: &Tori| -> Vec<u64> { t.tori.iter().flat_map(|t| bits(&[t.center[0], t.center[1], t.center[2], t.axis[0], t.axis[1], t.axis[2], t.major, t.minor])).collect() };
            assert!(tb(&tc) == tb(&tg), "{sign_rule:?} {eigen_rule:?}: tori differ");
            let strip = |mut r: FitReport| {
                r.seconds = 0.0;
                r
            };
            assert_eq!(strip(rc.clone()), strip(rg), "{sign_rule:?} {eigen_rule:?}: reports differ");
            eprintln!("fit {sign_rule:?}/{eigen_rule:?}: bitwise equal; saddle {} disc {} eigvec {} curvature {} zero {}", rc.saddle, rc.clamped_disc, rc.clamped_eigvec, rc.clamped_curvature, rc.zero_curvature);
        }
    }
    let mut bad = co.clone();
    bad[7][2] = f32::NAN;
    assert!(gfit(&g, &cloud, &bad, &sigma, 0, &FitOptions::default()).is_err());
    assert!(gfit(&g, &cloud, &co[..10], &sigma, 0, &FitOptions::default()).is_err());
}

#[test]
fn exp_against_host_libm() {
    let Some(ord) = gpu() else { return };
    let g = GpuContext::new(ord).unwrap();
    let mut s = 17;
    let x: Vec<f64> = (0..2_000_000).map(|_| 128.0 * lcg(&mut s) - 64.0).collect();
    let y = g.exp_values(&x).unwrap();
    let mut diff = 0usize;
    let mut maxulp = 0i64;
    for (a, b) in x.iter().zip(&y) {
        let c = a.exp();
        if c.to_bits() != b.to_bits() {
            diff += 1;
            maxulp = maxulp.max((c.to_bits() as i64 - b.to_bits() as i64).abs());
        }
    }
    eprintln!("exp on [-64, 64]: {diff} of {} differ from the host libm (max {maxulp} ulp)", x.len());
    assert!(maxulp <= 2);
}

/// A field on an ellipsoid cloud with the network's tori, on `device` with `opts`.
fn field_with(cloud: &PointCloud, tori: &Tori, mut opts: EvalOptions, device: Device) -> TorusField {
    opts.device = device;
    opts.threads = Some(16);
    TorusField::new(cloud, tori.clone(), opts).unwrap()
}

fn queries(nq: usize, seed: u64, cloud: &PointCloud) -> Vec<[f64; 3]> {
    let mut s = seed;
    let mut qs: Vec<[f64; 3]> = (0..nq).map(|_| [1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8]).collect();
    qs.extend(cloud.positions().iter().take(50).copied()); // at data points
    qs
}

struct Eval {
    d: Vec<f64>,
    g: Vec<[f64; 3]>,
    l: Vec<f64>,
    rep: Vec<EvalReport>,
}

fn eval_all(f: &TorusField, qs: &[[f64; 3]]) -> Eval {
    let n = qs.len();
    let mut d = vec![0.0; n];
    let mut g = vec![[0.0; 3]; n];
    let mut l = vec![0.0; n];
    let mut d2 = vec![0.0; n];
    let mut rep = vec![f.evaluate(qs, EvalOut { distance: Some(&mut d2), gradient: None, laplacian: None }).unwrap()];
    let mut d3 = vec![0.0; n];
    let mut g3 = vec![[0.0; 3]; n];
    rep.push(f.evaluate(qs, EvalOut { distance: Some(&mut d3), gradient: Some(&mut g3), laplacian: None }).unwrap());
    if f.options().gradient == GradientMode::Frozen {
        rep.push(f.evaluate(qs, EvalOut { distance: Some(&mut d), gradient: Some(&mut g), laplacian: Some(&mut l) }).unwrap());
    }
    d.extend(d2);
    d.extend(d3);
    g.extend(g3);
    Eval { d, g, l, rep }
}

fn rel_err(a: f64, b: f64) -> f64 {
    (a - b).abs() / (1.0 + a.abs())
}

/// Max |Δ|/(1 + |a|) of distance, gradient and Laplacian, and the bitwise-equal fraction.
fn compare(name: &str, c: &Eval, gq: &Eval, tol: [f64; 3]) -> [f64; 3] {
    let mut worst = [0f64; 3];
    let (mut same, mut total) = (0usize, 0usize);
    let mut acc = |k: usize, a: f64, b: f64| {
        worst[k] = worst[k].max(rel_err(a, b));
        same += (a.to_bits() == b.to_bits()) as usize;
        total += 1;
    };
    for (a, b) in c.d.iter().zip(&gq.d) {
        acc(0, *a, *b);
    }
    for (a, b) in c.g.iter().zip(&gq.g) {
        for t in 0..3 {
            acc(1, a[t], b[t]);
        }
    }
    for (a, b) in c.l.iter().zip(&gq.l) {
        acc(2, *a, *b);
    }
    for (rc, rg) in c.rep.iter().zip(&gq.rep) {
        assert_eq!(
            (rc.queries, rc.singular_axis, rc.singular_core, rc.at_data_point, rc.radius_fallbacks, rc.terms),
            (rg.queries, rg.singular_axis, rg.singular_core, rg.at_data_point, rg.radius_fallbacks, rg.terms),
            "{name}: counters differ"
        );
        assert_eq!(rc.device, Device::Cpu);
        assert!(matches!(rg.device, Device::Cuda { .. }), "{name}: the GPU report says {:?}", rg.device);
    }
    let frac = same as f64 / total as f64;
    eprintln!(
        "{name:<38} max |Δ|/(1+|a|): φ {:8.1e}, ∇φ {:8.1e}, Δφ {:8.1e}; bitwise equal {:6.2} %",
        worst[0],
        worst[1],
        worst[2],
        100.0 * frac
    );
    for k in 0..3 {
        assert!(worst[k] <= tol[k], "{name}: quantity {k}: {:e} > {:e}", worst[k], tol[k]);
    }
    worst
}

fn hash64(e: &Eval) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for x in e.d.iter().chain(&e.l).chain(e.g.iter().flatten()) {
        h = (h ^ x.to_bits()).wrapping_mul(0x100000001b3);
    }
    h
}

fn net_tori(cloud: &PointCloud) -> Tori {
    let model = Model::reference().unwrap();
    let (field, _) = PointsAsTori::build(cloud, &model, &BuildOptions::default().with_threads(Some(16))).unwrap();
    field.tori().clone()
}

#[test]
fn blend_matches_cpu_f64() {
    let Some(ord) = gpu() else { return };
    let cloud = ellipsoid_cloud(3000, 11, [0.6, 0.45, 0.3]);
    let tori = net_tori(&cloud);
    let qs = queries(3000, 5, &cloud);
    let mut s = 3;
    let area: Vec<f64> = (0..cloud.len()).map(|_| 0.5 + lcg(&mut s)).collect();
    let (paper, _) = TorusField::paper_options(&cloud, Some(16)).unwrap();
    let cases: Vec<(&str, EvalOptions)> = vec![
        ("knn32 per-query frozen", EvalOptions::default()),
        ("knn32 per-query exact", EvalOptions { gradient: GradientMode::Exact, ..Default::default() }),
        ("knn64 constant λ=40", EvalOptions { neighborhood: Neighborhood::Knn { k: 64 }, lambda: LambdaRule::Constant { lambda: 40.0 }, ..Default::default() }),
        ("knn100 per-query, areas, isovalue", EvalOptions { neighborhood: Neighborhood::Knn { k: 100 }, area_weights: Some(area.clone()), isovalue: 0.01, ..Default::default() }),
        ("knn k≥n (index order)", EvalOptions { neighborhood: Neighborhood::Knn { k: 5000 }, ..Default::default() }),
        ("all", EvalOptions { neighborhood: Neighborhood::All, gradient: GradientMode::Exact, ..Default::default() }),
        ("radius 0.1 fallback 32", EvalOptions { neighborhood: Neighborhood::Radius { r: 0.1, fallback_k: 32 }, ..Default::default() }),
        ("paper heuristic (radius, λ = 10³/D)", paper),
    ];
    for (name, o) in cases {
        let fc = field_with(&cloud, &tori, o.clone(), Device::Cpu);
        let fg = field_with(&cloud, &tori, o.clone(), Device::Cuda { ordinal: ord });
        assert!(matches!(fg.device(), Device::Cuda { .. }));
        let c = eval_all(&fc, &qs);
        let a = eval_all(&fg, &qs);
        let b = eval_all(&fg, &qs);
        assert_eq!(hash64(&a), hash64(&b), "{name}: two GPU runs differ (level 1)");
        compare(name, &c, &a, [1e-14, 1e-10, 1e-7]);
    }
}

#[test]
fn blend_f32_stable_tolerance() {
    let Some(ord) = gpu() else { return };
    for (cloud_name, cloud) in [
        ("ellipsoid", ellipsoid_cloud(3000, 11, [0.6, 0.45, 0.3])),
        ("flat ellipsoid 1:1:0.01", ellipsoid_cloud(3000, 12, [1.0, 1.0, 0.01])),
        ("plane, |κ| ≤ 2e-4 and 0 (R ~ 1e10)", plane_cloud(3000, 41)),
        ("ellipsoid at offset 1000", {
            let c = ellipsoid_cloud(3000, 13, [0.6, 0.45, 0.3]);
            let p: Vec<[f64; 3]> = c.positions().iter().map(|p| [p[0] + 1000.0, p[1] - 500.0, p[2] + 250.0]).collect();
            PointCloud::new(p, c.normals().to_vec(), ddg_points_as_tori::NormalPolicy::default()).unwrap().0
        }),
    ] {
        let tori = if cloud_name.starts_with("plane") { plane_tori(&cloud) } else { net_tori(&cloud) };
        let big = tori.tori.iter().map(|t| t.major.abs().max(t.minor.abs())).fold(0.0, f64::max);
        // random queries (moved with the cloud) and the data points themselves
        let mut s = 7;
        let mut qs: Vec<[f64; 3]> = (0..3000).map(|_| [1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8]).collect();
        if cloud_name.contains("offset") {
            qs = qs.iter().map(|q| [q[0] + 1000.0, q[1] - 500.0, q[2] + 250.0]).collect();
        }
        qs.extend(cloud.positions().iter().take(50).copied());
        for o in [EvalOptions::default(), EvalOptions { neighborhood: Neighborhood::All, ..Default::default() }] {
            let fc = field_with(&cloud, &tori, o.clone(), Device::Cpu);
            let fg = field_with(&cloud, &tori, EvalOptions { precision: EvalPrecision::F32Stable, ..o.clone() }, Device::Cuda { ordinal: ord });
            let n = qs.len();
            let (mut dc, mut dg) = (vec![0.0; n], vec![0.0; n]);
            let (mut gc, mut gg) = (vec![[0.0; 3]; n], vec![[0.0; 3]; n]);
            let (mut lc, mut lgp) = (vec![0.0; n], vec![0.0; n]);
            fc.evaluate(&qs, EvalOut { distance: Some(&mut dc), gradient: Some(&mut gc), laplacian: Some(&mut lc) }).unwrap();
            fg.evaluate(&qs, EvalOut { distance: Some(&mut dg), gradient: Some(&mut gg), laplacian: Some(&mut lgp) }).unwrap();
            let el = lc.iter().zip(&lgp).map(|(a, b)| (a - b).abs() / (1.0 + a.abs())).fold(0.0, f64::max);
            if std::env::var_os("PAT_DEBUG").is_some() {
                let wi = (0..n).max_by(|&i, &j| {
                    let e = |k: usize| ((gc[k][0] - gg[k][0]).powi(2) + (gc[k][1] - gg[k][1]).powi(2) + (gc[k][2] - gg[k][2]).powi(2)).sqrt();
                    e(i).partial_cmp(&e(j)).unwrap()
                }).unwrap();
                let li = (0..n).max_by(|&i, &j| ((lc[i] - lgp[i]).abs() / (1.0 + lc[i].abs())).partial_cmp(&((lc[j] - lgp[j]).abs() / (1.0 + lc[j].abs()))).unwrap()).unwrap();
                eprintln!("  worst grad at q{wi} {:?}: cpu {:?} gpu {:?} φ {} / {}", qs[wi], gc[wi], gg[wi], dc[wi], dg[wi]);
                eprintln!("  worst lap at q{li} {:?}: cpu {} gpu {} φ {} / {}", qs[li], lc[li], lgp[li], dc[li], dg[li]);
            }
            let mut dg2 = vec![0.0; n];
            fg.evaluate(&qs, EvalOut { distance: Some(&mut dg2), gradient: None, laplacian: None }).unwrap();
            assert!(bits(&dg) == bits(&dg2), "F32Stable: distance with and without the gradient differ");
            // the scale: the cloud's extent
            let ed = dc.iter().zip(&dg).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            let eg = gc.iter().zip(&gg).map(|(a, b)| ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()).fold(0.0, f64::max);
            eprintln!(
                "F32Stable {cloud_name:<26} {:?}: max |Δφ| {ed:9.2e}, |Δ∇φ| {eg:9.2e}, |ΔΔφ|/(1+|Δφ|) {el:9.2e} (largest radius {big:.1e})",
                o.neighborhood
            );
            assert!(ed < 3e-6 && eg < 5e-4 && el < 0.5, "{cloud_name}: f32 stable error too large");
        }
    }
}

#[test]
fn device_selection() {
    let cloud = ellipsoid_cloud(300, 3, [0.6, 0.45, 0.3]);
    let tori = Tori { tori: vec![ddg_points_as_tori::Torus { center: [0.0; 3], axis: [0.0, 0.0, 1.0], major: 1.0, minor: 0.5 }; 300], coefficients: vec![[0.0; 6]; 300], sigma: vec![1.0; 300] };
    // a bad ordinal is an error for Cuda, the CPU for Auto
    let bad = TorusField::new(&cloud, tori.clone(), EvalOptions { device: Device::Cuda { ordinal: 99 }, ..Default::default() });
    assert!(matches!(bad, Err(PatError::Gpu { .. })), "{bad:?}");
    assert!(ddg_points_as_tori::gpu::probe(99).is_err());
    let f32_cpu = TorusField::new(&cloud, tori.clone(), EvalOptions { precision: EvalPrecision::F32Stable, ..Default::default() });
    assert!(f32_cpu.is_err());
    let Some(ord) = gpu() else { return };
    let auto = TorusField::new(&cloud, tori.clone(), EvalOptions { device: Device::Auto, ..Default::default() }).unwrap();
    assert_eq!(auto.device(), Device::Cuda { ordinal: ord });
    let big_k = TorusField::new(&cloud, tori.clone(), EvalOptions { device: Device::Cuda { ordinal: ord }, neighborhood: Neighborhood::Knn { k: 200 }, ..Default::default() });
    assert!(matches!(big_k, Err(PatError::Gpu { .. })));
    let big_k_auto = TorusField::new(&cloud, tori, EvalOptions { device: Device::Auto, neighborhood: Neighborhood::Knn { k: 200 }, ..Default::default() }).unwrap();
    assert_eq!(big_k_auto.device(), Device::Cpu);
    let d = ddg_points_as_tori::gpu::devices().unwrap();
    assert!(!d.is_empty());
}

#[test]
fn degenerate_query_on_gpu() {
    let Some(ord) = gpu() else { return };
    // every neighbor at the query: PerQuery σ = 0
    let p = vec![[0.0, 0.0, 0.0]; 40];
    let nr = vec![[0.0, 0.0, 1.0]; 40];
    let (cloud, _) = PointCloud::new(p, nr, ddg_points_as_tori::NormalPolicy::default()).unwrap();
    let tori = Tori { tori: vec![ddg_points_as_tori::Torus { center: [0.0; 3], axis: [0.0, 0.0, 1.0], major: 1.0, minor: 0.5 }; 40], coefficients: vec![[0.0; 6]; 40], sigma: vec![1.0; 40] };
    let f = TorusField::new(&cloud, tori, EvalOptions { device: Device::Cuda { ordinal: ord }, ..Default::default() }).unwrap();
    let mut d = vec![0.0; 3];
    let r = f.evaluate(&[[1.0, 0.0, 0.0], [0.0; 3], [0.0, 1.0, 0.0]], EvalOut { distance: Some(&mut d), gradient: None, laplacian: None });
    assert!(matches!(r, Err(PatError::DegenerateQuery { index: 1, .. })), "{r:?}");
}

#[test]
fn pipeline_gpu_vs_cpu() {
    let Some(ord) = gpu() else { return };
    let model = Model::reference().unwrap();
    let cloud = ellipsoid_cloud(5000, 31, [0.6, 0.45, 0.3]);
    let cpu = BuildOptions::default().with_threads(Some(16));
    let mut gpuo = cpu.clone();
    gpuo.predict.device = Device::Cuda { ordinal: ord };
    gpuo.eval.device = Device::Cuda { ordinal: ord };
    let (fc, rc) = PointsAsTori::build(&cloud, &model, &cpu).unwrap();
    let (fg, rg) = PointsAsTori::build(&cloud, &model, &gpuo).unwrap();
    let (fg2, _) = PointsAsTori::build(&cloud, &model, &gpuo).unwrap();
    eprintln!("devices: knn {:?}, network {:?}, fit {:?}, eval {:?}", rg.knn_device, rg.predict.device, rg.fit_device, fg.device());
    assert_eq!(rc.predict.device, Device::Cpu);
    assert_eq!(rg.predict.device, Device::Cuda { ordinal: ord });
    // run to run on the GPU: bitwise
    assert!(fg.tori().coefficients == fg2.tori().coefficients, "two GPU builds differ");
    let tb = |t: &Tori| -> Vec<u64> { t.tori.iter().flat_map(|t| bits(&[t.center[0], t.center[1], t.center[2], t.axis[0], t.axis[1], t.axis[2], t.major, t.minor])).collect() };
    assert!(tb(fg.tori()) == tb(fg2.tori()));
    // against the CPU: coefficients within 1e-5 (1 + |a|)
    let mut worst = 0f64;
    for (a, b) in fc.tori().coefficients.iter().zip(&fg.tori().coefficients) {
        for t in 0..6 {
            worst = worst.max(((a[t] - b[t]).abs() / (1.0 + a[t].abs())) as f64);
        }
    }
    let qs = queries(4000, 9, &cloud);
    let (mut dc, mut dg) = (vec![0.0; qs.len()], vec![0.0; qs.len()]);
    fc.evaluate(&qs, EvalOut { distance: Some(&mut dc), gradient: None, laplacian: None }).unwrap();
    fg.evaluate(&qs, EvalOut { distance: Some(&mut dg), gradient: None, laplacian: None }).unwrap();
    let mut ds: Vec<f64> = dc.iter().zip(&dg).map(|(a, b)| (a - b).abs()).collect();
    ds.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eprintln!(
        "pipeline: coefficients max |Δa|/(1+|a|) {worst:.2e}; SDF |Δφ| median {:.1e}, p99.9 {:.1e}, max {:.1e}; saddles CPU {} GPU {}",
        ds[ds.len() / 2],
        ds[ds.len() * 999 / 1000],
        ds[ds.len() - 1],
        rc.fit.saddle,
        rg.fit.saddle
    );
    assert!(worst < 1e-5, "coefficients: {worst:e}");
}

/// Debug: the stable multi-pass blend of `kernels.cu` written in Rust, generic over f32/f64.
#[allow(dead_code)]
mod stable_ref {
    use ddg_points_as_tori::Torus;
    pub trait R: Copy + std::ops::Add<Output = Self> + std::ops::Sub<Output = Self> + std::ops::Mul<Output = Self> + std::ops::Div<Output = Self> + PartialOrd {
        fn f(x: f64) -> Self;
        fn sqrt(self) -> Self;
        fn exp(self) -> Self;
        fn to(self) -> f64;
    }
    impl R for f32 {
        fn f(x: f64) -> Self { x as f32 }
        fn sqrt(self) -> Self { f32::sqrt(self) }
        fn exp(self) -> Self { f32::exp(self) }
        fn to(self) -> f64 { self as f64 }
    }
    impl R for f64 {
        fn f(x: f64) -> Self { x }
        fn sqrt(self) -> Self { f64::sqrt(self) }
        fn exp(self) -> Self { f64::exp(self) }
        fn to(self) -> f64 { self }
    }
    fn sdf<T: R>(p: [f64; 3], t: &Torus, y: [T; 3]) -> (T, [T; 3], T) {
        let (w, vp, d1p, d2p, k2) = ddg_points_as_tori::gpu::stable_constants(p, t);
        let w = w.map(T::f);
        let u = t.axis.map(T::f);
        let (vp, d1p, d2p, k2) = (T::f(vp), T::f(d1p), T::f(d2p), T::f(k2));
        let two = T::f(2.0);
        let zero = T::f(0.0);
        let yu = y[0] * u[0] + y[1] * u[1] + y[2] * u[2];
        let a = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) + two * (y[0] * w[0] + y[1] * w[1] + y[2] * w[2]) - yu * yu - two * yu * d2p;
        let zu = yu + d2p;
        let pz = [y[0] + w[0] - zu * u[0], y[1] + w[1] - zu * u[1], y[2] + w[2] - zu * u[2]];
        let v = (pz[0] * pz[0] + pz[1] * pz[1] + pz[2] * pz[2]).sqrt();
        let den = v + vp;
        let e1 = if den > zero { a / den } else { zero };
        let e2 = yu;
        let (d1, d2) = (d1p + e1, d2p + e2);
        let nd = (d1 * d1 + d2 * d2).sqrt();
        let s = if t.minor >= 0.0 { T::f(1.0) } else { T::f(-1.0) };
        let r = T::f(t.minor.abs());
        let e = k2 + two * (d1p * e1 + d2p * e2) + (e1 * e1 + e2 * e2);
        let g = s * (e / (nd + r));
        let (w1, w2) = (d1 / nd, d2 / nd);
        let gr = [0, 1, 2].map(|c| s * (w1 * pz[c] / v + w2 * u[c]));
        let lg = s * (T::f(1.0) / nd + d1 / (nd * v));
        (g, gr, lg)
    }
    /// (φ, ∇φ, Δφ) at x over the neighbor list, PerQuery C.
    pub fn blend<T: R>(pos: &[[f64; 3]], tori: &[Torus], x: [f64; 3], nb: &[u32], c: f64) -> (f64, [f64; 3], f64) {
        let ys: Vec<[T; 3]> = nb.iter().map(|&i| { let p = pos[i as usize]; [T::f(x[0]) - T::f(p[0]), T::f(x[1]) - T::f(p[1]), T::f(x[2]) - T::f(p[2])] }).collect();
        let rho: Vec<T> = ys.iter().map(|y| (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]).sqrt()).collect();
        let mut m = T::f(0.0);
        for r in &rho { if *r > m { m = *r; } }
        let sigma = T::f(0.5) * m;
        let lambda = T::f(c) / sigma;
        let zero = T::f(0.0);
        let terms: Vec<(T, [T; 3], T, T)> = nb.iter().zip(&ys).zip(&rho).map(|((&i, y), &r)| {
            let (g, gr, lg) = sdf::<T>(pos[i as usize], &tori[i as usize], *y);
            let w = (zero - lambda * (r - sigma)).exp();
            (g, gr, lg, w)
        }).collect();
        let wsum = terms.iter().fold(zero, |a, t| a + t.3);
        let phi = terms.iter().fold(zero, |a, t| a + t.3 * t.0) / wsum;
        let mut gd = [zero; 3];
        for ((t, y), &r) in terms.iter().zip(&ys).zip(&rho) {
            let rh = if r > zero { [y[0] / r, y[1] / r, y[2] / r] } else { [zero; 3] };
            for c in 0..3 { gd[c] = gd[c] + t.3 * t.1[c] / wsum + lambda * (t.3 * rh[c] * (phi - t.0)) / wsum; }
        }
        let mut l = zero;
        for ((t, y), &r) in terms.iter().zip(&ys).zip(&rho) {
            let (rh, inv) = if r > zero { ([y[0] / r, y[1] / r, y[2] / r], T::f(2.0) / r) } else { ([zero; 3], zero) };
            let dt = rh[0] * (t.1[0] - gd[0]) + rh[1] * (t.1[1] - gd[1]) + rh[2] * (t.1[2] - gd[2]);
            l = l + t.3 * (t.2 - T::f(2.0) * lambda * dt + (t.0 - phi) * lambda * (lambda - inv));
        }
        (phi.to(), gd.map(|v| v.to()), (l / wsum).to())
    }
}

/// Points on the plane z = 0 in [−0.8, 0.8]² (unit normals +z).
fn plane_cloud(n: usize, seed: u64) -> PointCloud {
    let mut s = seed;
    let p: Vec<[f64; 3]> = (0..n).map(|_| [1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8, 0.0]).collect();
    PointCloud::new(p, vec![[0.0, 0.0, 1.0]; n], ddg_points_as_tori::NormalPolicy::default()).unwrap().0
}

/// Nearly flat tori on [`plane_cloud`]: curvatures up to 2e-4 (radii ≥ 5e3), every 5th exactly flat
/// (the 1e-10 floor: radii 1e10).
fn plane_tori(cloud: &PointCloud) -> Tori {
    let mut s = 77;
    let co: Vec<[f32; 6]> = (0..cloud.len())
        .map(|i| if i % 5 == 0 { [0.0; 6] } else { [0.0, 0.0, 0.0, (2e-4 * lcg(&mut s) - 1e-4) as f32, (2e-4 * lcg(&mut s) - 1e-4) as f32, (2e-4 * lcg(&mut s) - 1e-4) as f32] })
        .collect();
    let o = FitOptions { sign_rule: SignRule::Geometric, ..Default::default() };
    ddg_points_as_tori::fit_tori(cloud, &co, &o).unwrap().0
}

/// The stable form of `kernels.cu` (f32 path) written in f64 equals the CPU's Eq. 23 blend: the
/// rearrangement is exact, so the F32Stable error is f32 round-off only.
#[test]
fn stable_form_in_f64_equals_cpu() {
    let cloud = ellipsoid_cloud(2000, 11, [0.6, 0.45, 0.3]);
    let tori = net_tori(&cloud);
    let fc = field_with(&cloud, &tori, EvalOptions::default(), Device::Cpu);
    let mut s = 23;
    let mut worst = [0f64; 3];
    let mut qs: Vec<[f64; 3]> = (0..500).map(|_| [1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8]).collect();
    qs.extend(cloud.positions().iter().take(20).copied());
    for x in qs {
        let (v, nb) = fc.blend_at(x, ddg_points_as_tori::blend::BlendWant::Laplacian).unwrap();
        let a = stable_ref::blend::<f64>(cloud.positions(), &tori.tori, x, &nb, 64.0);
        worst[0] = worst[0].max(rel_err(v.distance, a.0));
        for t in 0..3 {
            worst[1] = worst[1].max(rel_err(v.gradient[t], a.1[t]));
        }
        worst[2] = worst[2].max(rel_err(v.laplacian, a.2));
    }
    eprintln!("stable form in f64 vs Eq. 23 (CPU): φ {:.1e}, ∇φ {:.1e}, Δφ {:.1e}", worst[0], worst[1], worst[2]);
    assert!(worst[0] < 1e-13 && worst[1] < 1e-11 && worst[2] < 1e-9, "{worst:?}");
}

// ------------------------------------------------------------------------------------------------
// Timings (ignored; run alone: `--ignored timing_ --test-threads 1 --nocapture`). PAT_SIZES sets the
// cloud sizes (default 100000,1000000), PAT_CPU=0 skips the CPU rows.
// ------------------------------------------------------------------------------------------------

fn machine(tag: &str) {
    let up = std::process::Command::new("uptime").output().map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string()).unwrap_or_default();
    let smi = std::process::Command::new("nvidia-smi")
        .args(["--query-gpu=memory.used,utilization.gpu,temperature.gpu,clocks.sm", "--format=csv,noheader"])
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default();
    eprintln!("[{tag}] uptime: {up} | nvidia-smi: {smi}");
}

/// `n` points sampled area-uniformly on Keenan Crane's spot (data/models), with interpolated vertex
/// normals; `None` when the mesh is absent.
fn spot_cloud(n: usize, seed: u64) -> Option<PointCloud> {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models/spot/spot/spot_triangulated.obj");
    let mesh = ddg_core::io::read_obj(&path).ok()?;
    let tris = mesh.triangles()?;
    let v = &mesh.positions;
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let mut vn = vec![[0.0f64; 3]; v.len()];
    let mut cum = Vec::with_capacity(tris.len());
    let mut tot = 0.0;
    for t in &tris {
        let nf = cross(sub(v[t[1] as usize], v[t[0] as usize]), sub(v[t[2] as usize], v[t[0] as usize]));
        for &i in t {
            for a in 0..3 {
                vn[i as usize][a] += nf[a];
            }
        }
        tot += 0.5 * (nf[0] * nf[0] + nf[1] * nf[1] + nf[2] * nf[2]).sqrt();
        cum.push(tot);
    }
    let mut s = seed;
    let mut p = Vec::with_capacity(n);
    let mut nr = Vec::with_capacity(n);
    for _ in 0..n {
        let f = cum.partition_point(|&c| c < lcg(&mut s) * tot).min(tris.len() - 1);
        let (r1, r2) = (lcg(&mut s).sqrt(), lcg(&mut s));
        let w = [1.0 - r1, r1 * (1.0 - r2), r1 * r2];
        let t = tris[f];
        let mut q = [0.0; 3];
        let mut m = [0.0; 3];
        for (k, &i) in t.iter().enumerate() {
            for a in 0..3 {
                q[a] += w[k] * v[i as usize][a];
                m[a] += w[k] * vn[i as usize][a];
            }
        }
        let l = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
        p.push(q);
        nr.push([m[0] / l, m[1] / l, m[2] / l]);
    }
    Some(PointCloud::new(p, nr, ddg_points_as_tori::NormalPolicy::default()).ok()?.0)
}

fn sizes() -> Vec<usize> {
    std::env::var("PAT_SIZES").ok().map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect()).unwrap_or(vec![100_000, 1_000_000])
}

fn build_row(tag: &str, n: usize, r: &ddg_points_as_tori::BuildReport) {
    eprintln!(
        "{tag:<22} n {n:>8}: kNN {:7.2} s ({:?}), features {:7.2} s (CPU{}), network {:7.2} s ({:?}), predict wall {:7.2} s, fit {:6.2} s ({:?}), total {:7.2} s = {:.1} s per 1M",
        r.knn_seconds,
        r.knn_device,
        r.predict.features_seconds,
        if matches!(r.predict.device, Device::Cpu) { ", summed over workers" } else { ", overlapped" },
        r.predict.forward_seconds,
        r.predict.device,
        r.predict.seconds,
        r.fit.seconds,
        r.fit_device,
        r.seconds,
        r.seconds * 1e6 / n as f64
    );
}

#[test]
#[ignore]
fn timing_precompute() {
    let Some(ord) = gpu() else { return };
    let model = Model::reference().unwrap();
    let cpu_rows = std::env::var("PAT_CPU").map(|v| v != "0").unwrap_or(true);
    for n in sizes() {
        let Some(cloud) = spot_cloud(n, 1) else {
            eprintln!("SKIP: data/models/spot missing");
            return;
        };
        // stage comparisons: fitting kNN and fit, CPU (64 threads) vs GPU
        let pool = ddg_points_as_tori::pool::make_pool(None).unwrap();
        let t = std::time::Instant::now();
        let tree = pool.install(|| KdTree::new(cloud.positions(), KdTreeOptions::default())).unwrap();
        let t_tree = t.elapsed().as_secs_f64();
        let t = std::time::Instant::now();
        let nb_cpu = pool.install(|| fitting_neighborhoods(&tree, cloud.positions(), 64));
        let t_knn_cpu = t.elapsed().as_secs_f64();
        let g = GpuContext::new(ord).unwrap();
        let _ = gfit_nb(&g, &tree, &cloud.positions()[..1000.min(n)], 64).ok();
        let t = std::time::Instant::now();
        let nb_gpu = gfit_nb(&g, &tree, cloud.positions(), 64).unwrap();
        let t_knn_gpu = t.elapsed().as_secs_f64();
        assert_eq!(nb_cpu, nb_gpu);
        machine("stages");
        eprintln!("stages n {n}: kd-tree (CPU) {t_tree:.3} s, 64-NN CPU 64 threads {t_knn_cpu:.3} s, GPU {t_knn_gpu:.3} s (incl. upload/download)");
        drop((nb_cpu, nb_gpu));
        for mode in [GpuMatmul::Split3, GpuMatmul::Tf32] {
            let mut o = BuildOptions::default();
            o.predict.device = Device::Cuda { ordinal: ord };
            o.predict.gpu_matmul = mode;
            // warm-up: the first forward pass JIT-compiles CubeCL's fused kernels (once per process)
            let small = spot_cloud(3000, 9).unwrap();
            let t = std::time::Instant::now();
            let (_, rw) = PointsAsTori::build(&small, &model, &o).unwrap();
            eprintln!("warm-up {mode:?}: 3000 points in {:.2} s (network {:.2} s, includes the JIT compile)", t.elapsed().as_secs_f64(), rw.predict.forward_seconds);
            machine("before");
            let (field, r) = PointsAsTori::build(&cloud, &model, &o).unwrap();
            machine("after");
            build_row(&format!("GPU {mode:?}"), n, &r);
            if mode == GpuMatmul::Split3 {
                let co = field.tori().coefficients.clone();
                let sg = field.tori().sigma.clone();
                let t = std::time::Instant::now();
                let _ = fit_tori_with(&cloud, &co, &sg, 0, &FitOptions::default()).unwrap();
                let t_fit_cpu = t.elapsed().as_secs_f64();
                let t = std::time::Instant::now();
                let _ = gfit(&g, &cloud, &co, &sg, 0, &FitOptions::default()).unwrap();
                eprintln!("stages n {n}: fit CPU 64 threads {t_fit_cpu:.3} s, GPU {:.3} s (incl. transfers)", t.elapsed().as_secs_f64());
            }
        }
        if cpu_rows && n <= 1_000_000 {
            machine("before");
            let (_, r) = PointsAsTori::build(&cloud, &model, &BuildOptions::default()).unwrap();
            machine("after");
            build_row("CPU 64 threads", n, &r);
        }
    }
}

#[test]
#[ignore]
fn timing_evaluation() {
    let Some(ord) = gpu() else { return };
    let model = Model::reference().unwrap();
    let n = std::env::var("PAT_EVAL_N").ok().and_then(|s| s.parse().ok()).unwrap_or(100_000usize);
    let Some(cloud) = spot_cloud(n, 2) else { return };
    let mut o = BuildOptions::default();
    o.predict.device = Device::Cuda { ordinal: ord };
    let (field, _) = PointsAsTori::build(&cloud, &model, &o).unwrap();
    let tori = field.tori().clone();
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in cloud.positions() {
        for a in 0..3 {
            lo[a] = lo[a].min(p[a]);
            hi[a] = hi[a].max(p[a]);
        }
    }
    let grid = |m: usize| -> Vec<[f64; 3]> {
        let mut g = Vec::with_capacity(m * m * m);
        for i in 0..m {
            for j in 0..m {
                for k in 0..m {
                    let f = |a: usize, t: usize| lo[a] - 0.05 + (hi[a] - lo[a] + 0.1) * t as f64 / (m - 1) as f64;
                    g.push([f(0, i), f(1, j), f(2, k)]);
                }
            }
        }
        g
    };
    let fields = [
        ("CPU f64 64 threads", field_with(&cloud, &tori, EvalOptions { threads: None, ..Default::default() }, Device::Cpu)),
        ("GPU f64", field_with(&cloud, &tori, EvalOptions::default(), Device::Cuda { ordinal: ord })),
        ("GPU F32Stable", field_with(&cloud, &tori, EvalOptions { precision: EvalPrecision::F32Stable, ..Default::default() }, Device::Cuda { ordinal: ord })),
    ];
    for m in [64usize, 256] {
        let qs = grid(m);
        for (name, f) in &fields {
            for (what, grad, lap) in [("φ", false, false), ("φ+∇φ", true, false), ("φ+∇φ+Δφ", true, true)] {
                if m == 256 && lap {
                    continue;
                }
                let mut d = vec![0.0; qs.len()];
                let mut gr = vec![[0.0; 3]; if grad { qs.len() } else { 0 }];
                let mut l = vec![0.0; if lap { qs.len() } else { 0 }];
                machine("before");
                let t = std::time::Instant::now();
                let rep = f
                    .evaluate(
                        &qs,
                        EvalOut { distance: Some(&mut d), gradient: if grad { Some(&mut gr) } else { None }, laplacian: if lap { Some(&mut l) } else { None } },
                    )
                    .unwrap();
                let s = t.elapsed().as_secs_f64();
                eprintln!("eval n {n} grid {m}³ {name:<20} {what:<10}: {s:8.3} s ({:.2e} s/query, device {:?})", s / qs.len() as f64, rep.device);
            }
        }
    }
    // per-query latency: single queries on the CPU (one thread, sequential), one-query batches on the GPU
    let mut s = 3;
    let qs: Vec<[f64; 3]> = (0..20_000).map(|_| [lo[0] + (hi[0] - lo[0]) * lcg(&mut s), lo[1] + (hi[1] - lo[1]) * lcg(&mut s), lo[2] + (hi[2] - lo[2]) * lcg(&mut s)]).collect();
    machine("latency");
    let t = std::time::Instant::now();
    let mut acc = 0.0;
    for q in &qs {
        acc += fields[0].1.distance(*q).unwrap();
    }
    let cpu_lat = t.elapsed().as_secs_f64() / qs.len() as f64;
    let t = std::time::Instant::now();
    let mut acc2 = 0.0;
    for q in &qs {
        let (_, gq) = fields[0].1.distance_gradient(*q).unwrap();
        acc2 += gq[0];
    }
    let cpu_lat_g = t.elapsed().as_secs_f64() / qs.len() as f64;
    let mut d = [0.0];
    let t = std::time::Instant::now();
    for q in qs.iter().take(2000) {
        fields[1].1.evaluate(std::slice::from_ref(q), EvalOut { distance: Some(&mut d), gradient: None, laplacian: None }).unwrap();
        acc += d[0];
    }
    let gpu_lat = t.elapsed().as_secs_f64() / 2000.0;
    eprintln!(
        "latency n {n}: CPU sequential 32-NN φ {cpu_lat:.2e} s, φ+∇φ {cpu_lat_g:.2e} s per query; GPU one-query batch {gpu_lat:.2e} s ({acc:.3} {acc2:.3})"
    );
}

/// Five identical builds per matmul mode: the network time per build and the coefficient hash
/// (bitwise equal in every repetition, NOTES-GPU.md).
#[test]
#[ignore]
fn timing_repeated_builds() {
    let Some(ord) = gpu() else { return };
    let model = Model::reference().unwrap();
    let n: usize = std::env::var("PAT_N").ok().and_then(|s| s.parse().ok()).unwrap_or(50_000);
    let c = spot_cloud(n, 1).unwrap();
    for mode in [GpuMatmul::Tf32, GpuMatmul::Split3] {
        let mut o = BuildOptions::default();
        o.predict.device = Device::Cuda { ordinal: ord };
        o.predict.gpu_matmul = mode;
        for rep in 0..5 {
            let (f, r) = PointsAsTori::build(&c, &model, &o).unwrap();
            let h = hash32(&f.tori().coefficients.iter().flatten().copied().collect::<Vec<f32>>());
            machine(&format!("{mode:?} rep {rep}: network {:.2} s = {:.0}/s, coefficients {h:016x}", r.predict.forward_seconds, n as f64 / r.predict.forward_seconds));
        }
    }
}

/// One network on the device, six passes over the same 50,000 neighborhoods.
#[test]
#[ignore]
fn timing_one_network() {
    let Some(ord) = gpu() else { return };
    let model = Model::reference().unwrap();
    let n = 50_000;
    let f = features(&spot_cloud(n, 1).unwrap(), 64);
    for mode in [GpuMatmul::Tf32, GpuMatmul::Split3] {
        let net = ddg_points_as_tori::network::cuda::CudaNetwork::new(&model, ord, mode).unwrap();
        for rep in 0..6 {
            let t = std::time::Instant::now();
            let y = net.forward(&f, 64).unwrap();
            let s = t.elapsed().as_secs_f64();
            machine(&format!("one runner {mode:?} rep {rep}: {:.0}/s {:016x}", n as f64 / s, hash32(&y)));
        }
    }
}
