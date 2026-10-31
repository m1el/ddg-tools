//! Random Dijkstra paths (paper Sec. 5): against the exact polyhedral distance (SPEC §7
//! test 4), determinism (§7 test 10), and the evaluation on `data/` meshes (an ignored test,
//! run in release mode alone for its timings: SPEC §9.1, §9.3).

mod common;

use std::time::Instant;

use common::{Shape, data_shapes, load_data, random_pairs, run};
use ddg_core::generators::{self, annulus, bumpy_grid_disk, icosphere, torus};
use ddg_flip_geodesics::{
    CertifyOptions, FlipError, GeodesicNetwork, ShortenOptions, ShortenStatus, Side, dijkstra_path,
};
use ddg_intrinsic::IntrinsicTriangulation;
use ddg_intrinsic_apps::exact::{ExactOptions, exact_distance};

/// Per shape: `n_sources` sources, `n_targets` targets each; the exact distances from each
/// source once.
fn against_exact(shape: &Shape, n_sources: usize, n_targets: usize) -> (usize, usize, usize, f64) {
    let tri = shape.triangulation();
    let n = shape.mesh.n_vertices();
    let (mut minimal, mut near, mut total, mut worst) = (0, 0, 0, 0.0f64);
    for (k, (source, _)) in random_pairs(n, n_sources, 17).into_iter().enumerate() {
        let d = exact_distance(
            &shape.mesh,
            &shape.geometry,
            &[source],
            &ExactOptions::default(),
        )
        .unwrap_or_else(|e| panic!("{}: {e}", shape.name))
        .distances;
        for (target, _) in random_pairs(n, n_targets, 1000 + k as u64) {
            if target == source || !d[target].is_finite() {
                continue;
            }
            let path = dijkstra_path(&tri, source, target).unwrap();
            let out = run(tri.clone(), &[path], &[], &ShortenOptions::default());
            assert!(
                out.certificate.geodesic,
                "{}: {source} -> {target}: {:?}",
                shape.name, out.certificate.failures
            );
            let length = out.network.length();
            assert!(
                length >= d[target] * (1.0 - 1e-9),
                "{}: {source} -> {target}: {length} below the exact distance {}",
                shape.name,
                d[target]
            );
            if length <= d[target] * (1.0 + 1e-6) {
                minimal += 1;
            }
            if length <= d[target] * (1.0 + 1e-3) {
                near += 1;
            }
            worst = worst.max(length / d[target] - 1.0);
            if std::env::var("FO_VERBOSE").is_ok() && length > d[target] * (1.0 + 1e-6) {
                eprintln!(
                    "  {}: {source} -> {target}: excess {:.3e}, flips {}, m0 -> m {}",
                    shape.name,
                    length / d[target] - 1.0,
                    out.report.flips,
                    out.network.network().n_segments()
                );
            }
            total += 1;
        }
    }
    (minimal, near, total, worst)
}

/// SPEC §7 test 4: `|γ| ≥ d_exact (1 − 1e-9)` always; the fraction with `|γ| ≤ d_exact (1 + 1e-6)`
/// is reported (paper Sec. 6.5: more than 95 % minimal).
#[test]
fn never_shorter_than_the_exact_distance() {
    let mut jittered = generators::icosphere(2).unwrap();
    generators::perturb(&mut jittered.positions, [0.05; 3], 7);
    let mut shapes = vec![
        Shape::generated("icosphere 3", icosphere(3).unwrap()),
        Shape::generated("jittered icosphere 2", jittered),
        Shape::generated("torus 24x10", torus(24, 10, 1.0, 0.35).unwrap()),
        Shape::generated("bumpy grid 14", bumpy_grid_disk(14, 14, 0.25).unwrap()),
        Shape::generated("annulus 24x4", annulus(24, 4, 0.4, 1.0).unwrap()),
    ];
    shapes.extend(load_data("models/spot/spot/spot_triangulated.obj"));
    let (mut all_minimal, mut all_near, mut all_total) = (0, 0, 0);
    for shape in &shapes {
        let (minimal, near, total, worst) = against_exact(shape, 6, 10);
        eprintln!(
            "{}: {minimal}/{total} minimal to 1e-6, {near}/{total} to 1e-3, worst excess {:.3} %",
            shape.name,
            100.0 * worst
        );
        all_minimal += minimal;
        all_near += near;
        all_total += total;
    }
    eprintln!(
        "all: {all_minimal}/{all_total} = {:.1} % minimal to 1e-6, {:.1} % to 1e-3 (the others \
         are other local geodesics: every result is certified)",
        100.0 * all_minimal as f64 / all_total as f64,
        100.0 * all_near as f64 / all_total as f64
    );
    assert!(all_total >= 250);
}

/// SPEC §7 test 10: reruns are bitwise identical, whichever thread runs them and however many
/// run at once (the algorithm is single-threaded, F24); refused inputs change nothing.
#[test]
fn deterministic_and_refusals_change_nothing() {
    let mut shapes = vec![
        Shape::generated("torus 16x8", torus(16, 8, 1.0, 0.35).unwrap()),
        Shape::generated("bumpy grid 10", bumpy_grid_disk(10, 10, 0.25).unwrap()),
    ];
    shapes.extend(data_shapes());
    type Result1 = (Vec<ddg_flip_geodesics::EdgePath>, Vec<u64>, u64, String);
    let one = |tri: &IntrinsicTriangulation, a: usize, b: usize| -> Result1 {
        let path = dijkstra_path(tri, a, b).unwrap();
        let mut net = GeodesicNetwork::new(tri.clone(), &[path], &[], None).unwrap();
        let report = net.shorten(&ShortenOptions::default()).unwrap();
        (
            net.paths(),
            net.triangulation()
                .edge_lengths()
                .iter()
                .map(|x| x.to_bits())
                .collect(),
            net.length().to_bits(),
            format!("{report:?}"),
        )
    };
    for shape in &shapes {
        let tri = shape.triangulation();
        let pairs = random_pairs(shape.mesh.n_vertices(), 16, 5);
        let serial: Vec<Result1> = pairs.iter().map(|&(a, b)| one(&tri, a, b)).collect();
        let parallel: Vec<Result1> = std::thread::scope(|s| {
            let handles: Vec<_> = pairs
                .iter()
                .map(|&(a, b)| {
                    let tri = &tri;
                    s.spawn(move || one(tri, a, b))
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        assert!(serial == parallel, "{}: reruns differ", shape.name);

        // Refused: a joint that is not shortenable, bad options.
        let path = dijkstra_path(&tri, pairs[0].0, pairs[0].1).unwrap();
        let net = GeodesicNetwork::new(tri.clone(), &[path], &[], None).unwrap();
        let mut copy = net.clone();
        for joint in net.network().joints() {
            for side in [Side::Left, Side::Right] {
                let g = net.joint_geometry(joint).unwrap();
                if g.angle(side) < std::f64::consts::PI && g.is_flexible(side) {
                    continue;
                }
                let r = copy.shorten_joint(joint, side);
                assert!(matches!(r, Err(FlipError::NotShortenable { .. })), "{r:?}");
            }
        }
        let bad = ShortenOptions {
            straight_tolerance: f64::NAN,
            ..Default::default()
        };
        assert!(copy.shorten(&bad).is_err());
        assert_eq!(copy.paths(), net.paths());
        assert_eq!(copy.length().to_bits(), net.length().to_bits());
        assert!(
            copy.triangulation()
                .edge_lengths()
                .iter()
                .zip(net.triangulation().edge_lengths())
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
}

/// One evaluation trial.
struct Trial {
    m: usize,
    geodesic: bool,
    status: ShortenStatus,
    flips: usize,
    iterations: usize,
    refused: usize,
    stalled: usize,
    sweeps: usize,
    stale: usize,
    dijkstra_s: f64,
    flipout_s: f64,
    ratio_exact: Option<f64>,
}

fn quantile(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(f64::total_cmp);
    v[((v.len() - 1) as f64 * q).round() as usize]
}

/// The evaluation of SPEC §9.3 on the meshes of `data/` (paper Sec. 5): per mesh `FO_TRIALS`
/// (default 200) random vertex pairs, Dijkstra path, Alg. 2; success is the certificate at the
/// paper's `1e-4` rad (Sec. 5.3); times of Dijkstra and of FlipOut (the clone of the
/// triangulation excluded, F16; no extraction, which is work package F3); `|γ| / d_exact` on
/// meshes up to 20 000 vertices. Run alone, in release mode, single-threaded:
/// `cargo test --release -p ddg-flip-geodesics --test random_paths -- --ignored --nocapture
/// --test-threads=1`.
#[test]
#[ignore]
fn evaluation_on_data_meshes() {
    let trials: usize = std::env::var("FO_TRIALS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    let files = [
        "models/spot/spot/spot_triangulated.obj",
        "models/bob/bob/bob_tri.obj",
        "models/blub/blub/blub_triangulated.obj",
        "models/crumpleddevelopable/CrumpledDevelopable.obj",
        "optimal-cones/optimal-cones-data/data/camelbad.obj",
        "optimal-cones/optimal-cones-data/data/bunny.obj",
        "optimal-cones/optimal-cones-data/data/cow1.obj",
        "optimal-cones/optimal-cones-data/data/rhino.obj",
        "optimal-cones/optimal-cones-data/data/tree.obj",
    ];
    eprintln!(
        "| mesh | V | trials | geodesic | median ms (FlipOut) | p99 ms | max ms | FlipOut/Dijkstra (median) | median flips | max flips | median m | refused | stalled | sweeps>0 | minimal 1e-6 |"
    );
    eprintln!("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|");
    let mut all: Vec<(f64, f64)> = Vec::new();
    for file in files {
        let Some(shape) = load_data(file) else {
            continue;
        };
        let tri = shape.triangulation();
        let n = shape.mesh.n_vertices();
        let pairs = random_pairs(n, trials, 2020);
        let mut results = Vec::new();
        // Exact distances per distinct source (one source per 10 trials, as the paper's 10
        // trials per mesh).
        let exact_ok = n <= 20_000;
        let mut exact_cache: Option<(usize, Vec<f64>)> = None;
        for (t, &(a0, b)) in pairs.iter().enumerate() {
            let a = if exact_ok { pairs[t - t % 10].0 } else { a0 };
            if a == b {
                continue;
            }
            let t0 = Instant::now();
            let Ok(path) = dijkstra_path(&tri, a, b) else {
                continue;
            };
            let dijkstra_s = t0.elapsed().as_secs_f64();
            let m = path.segments.len();
            let mut net = GeodesicNetwork::new(tri.clone(), &[path], &[], None).unwrap();
            let t1 = Instant::now();
            let report = net.shorten(&ShortenOptions::default()).unwrap();
            let flipout_s = t1.elapsed().as_secs_f64();
            let cert = net.certify(&CertifyOptions::default()).unwrap();
            let ratio_exact = exact_ok.then(|| {
                if exact_cache.as_ref().is_none_or(|(s, _)| *s != a) {
                    let d = exact_distance(
                        &shape.mesh,
                        &shape.geometry,
                        &[a],
                        &ExactOptions::default(),
                    )
                    .unwrap()
                    .distances;
                    exact_cache = Some((a, d));
                }
                net.length() / exact_cache.as_ref().unwrap().1[b]
            });
            results.push(Trial {
                m,
                geodesic: cert.geodesic,
                status: report.status,
                flips: report.flips,
                iterations: report.iterations,
                refused: report.refused_flips,
                stalled: report.stalled,
                sweeps: report.sweeps,
                stale: report.stale_entries,
                dijkstra_s,
                flipout_s,
                ratio_exact,
            });
        }
        let k = results.len();
        let geodesic = results.iter().filter(|r| r.geodesic).count();
        let mut times: Vec<f64> = results.iter().map(|r| 1e3 * r.flipout_s).collect();
        let mut ratio: Vec<f64> = results
            .iter()
            .map(|r| r.flipout_s / r.dijkstra_s.max(1e-9))
            .collect();
        let mut flips: Vec<f64> = results.iter().map(|r| r.flips as f64).collect();
        let mut ms: Vec<f64> = results.iter().map(|r| r.m as f64).collect();
        let refused: usize = results.iter().map(|r| r.refused).sum();
        let stalled: usize = results.iter().map(|r| r.stalled).sum();
        let swept = results.iter().filter(|r| r.sweeps > 0).count();
        let minimal = results
            .iter()
            .filter(|r| r.ratio_exact.is_some_and(|x| x <= 1.0 + 1e-6))
            .count();
        let with_exact = results.iter().filter(|r| r.ratio_exact.is_some()).count();
        let below = results
            .iter()
            .filter(|r| r.ratio_exact.is_some_and(|x| x < 1.0 - 1e-9))
            .count();
        assert_eq!(
            below, 0,
            "{}: a path shorter than the exact distance",
            shape.name
        );
        let not_geodesic: Vec<String> = results
            .iter()
            .filter(|r| !r.geodesic || r.status != ShortenStatus::Geodesic)
            .map(|r| {
                format!(
                    "{:?} m {} it {} stale {}",
                    r.status, r.m, r.iterations, r.stale
                )
            })
            .collect();
        eprintln!(
            "| {} | {n} | {k} | {geodesic}/{k} | {:.3} | {:.3} | {:.3} | {:.2} | {} | {} | {} | {refused} | {stalled} | {swept} | {} |",
            shape.name,
            quantile(&mut times, 0.5),
            quantile(&mut times, 0.99),
            quantile(&mut times, 1.0),
            quantile(&mut ratio, 0.5),
            quantile(&mut flips, 0.5),
            quantile(&mut flips, 1.0),
            quantile(&mut ms, 0.5),
            if with_exact > 0 {
                format!("{minimal}/{with_exact}")
            } else {
                "-".into()
            },
        );
        if !not_geodesic.is_empty() {
            eprintln!(
                "  {}: not geodesic or not status Geodesic: {not_geodesic:?}",
                shape.name
            );
        }
        all.extend(results.iter().map(|r| (r.m as f64, r.flips as f64)));
    }
    // Fit flips ∝ m^p over trials with flips > 0 (paper Fig. 16: p ≈ 1.5).
    let pts: Vec<(f64, f64)> = all
        .iter()
        .filter(|&&(m, f)| m > 1.0 && f > 0.0)
        .map(|&(m, f)| (m.ln(), f.ln()))
        .collect();
    if pts.len() > 2 {
        let n = pts.len() as f64;
        let (sx, sy) = pts
            .iter()
            .fold((0.0, 0.0), |(a, b), &(x, y)| (a + x, b + y));
        let (mx, my) = (sx / n, sy / n);
        let (sxy, sxx) = pts.iter().fold((0.0, 0.0), |(a, b), &(x, y)| {
            (a + (x - mx) * (y - my), b + (x - mx) * (x - mx))
        });
        eprintln!("flips ∝ m^{:.2} over {} trials", sxy / sxx, pts.len());
    }
}
