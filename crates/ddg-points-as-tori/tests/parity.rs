//! Parity with the authors' program (SPEC.md §15.1-4, §16.1) on the 12 golden inputs of
//! `runs_tmp/pat-ref/` (the reference's `data/` clouds, four training clouds, two synthetic
//! clouds). Skips with a message when the golden files are absent.
//!
//! Run with `--release --nocapture` for the full 64³ grids and the tables; debug builds use every
//! 4th grid point per axis.

#![allow(clippy::needless_range_loop)]

mod common;

use common::{Golden, INPUTS, diff_stats, golden_root};
use ddg_points_as_tori::features::{FeatureOrder, FeatureScratch, neighborhood_features};
use ddg_points_as_tori::network::{PredictOptions, forward_features, predict_coefficients, scale_outputs};
use ddg_points_as_tori::spatial::{KdTree, KdTreeOptions, Neighborhoods, fitting_neighborhoods};
use ddg_points_as_tori::{
    BuildOptions, EvalOptions, EvalOut, FitOptions, Model, Neighborhood, PointsAsTori, Tori, TorusField, fit_tori,
};

fn grid_subset(g: &Golden) -> (Vec<[f64; 3]>, Vec<usize>) {
    let grid = g.grid();
    let step = if cfg!(debug_assertions) { 4 } else { 1 };
    let mut ids = Vec::new();
    for ix in (0..64).step_by(step) {
        for iy in (0..64).step_by(step) {
            for iz in (0..64).step_by(step) {
                ids.push((ix * 64 + iy) * 64 + iz);
            }
        }
    }
    (ids.iter().map(|&i| grid[i]).collect(), ids)
}

fn golden_nbhds(g: &Golden) -> Neighborhoods {
    let a = g.arr("knn_indices");
    Neighborhoods { k: a.shape[1], indices: a.data.iter().map(|&x| x as u32).collect() }
}

#[test]
fn knn_matches_reference() {
    let Some(root) = golden_root() else { return };
    for name in INPUTS {
        let g = Golden::open(&root, name);
        let p = g.cloud.positions();
        let tree = KdTree::new(p, KdTreeOptions::default()).unwrap();
        let ours = fitting_neighborhoods(&tree, p, 64);
        let refn = golden_nbhds(&g);
        assert_eq!(ours.k, refn.k, "{name}: k");
        let (mut same_order, mut same_set, mut ties) = (0, 0, 0);
        for i in 0..p.len() {
            let (a, b) = (ours.row(i), refn.row(i));
            assert_eq!(a[0] as usize, i, "{name}: row {i} must start with itself");
            if a == b {
                same_order += 1;
                continue;
            }
            let (mut sa, mut sb) = (a.to_vec(), b.to_vec());
            sa.sort_unstable();
            sb.sort_unstable();
            if sa == sb {
                same_set += 1;
                continue;
            }
            // different sets: only allowed at a tie of the farthest distance
            let d = |j: u32| {
                let q = p[j as usize];
                let c = p[i];
                (q[0] - c[0]).powi(2) + (q[1] - c[1]).powi(2) + (q[2] - c[2]).powi(2)
            };
            let far = d(*a.last().unwrap());
            for &j in a.iter().filter(|j| !b.contains(j)) {
                assert_eq!(d(j), far, "{name}: row {i}: our neighbor {j} is not in the reference's row and not a tie");
            }
            ties += 1;
        }
        eprintln!("knn {name:>12}: n {:>6}, same order {same_order}, same set {same_set}, tie at k {ties}", p.len());
    }
}

#[test]
fn sigma_and_network_match_reference() {
    let Some(root) = golden_root() else { return };
    let model = Model::reference().unwrap();
    let opts = PredictOptions::default();
    let mut worst_raw: f64 = 0.0;
    let mut worst_coef: f64 = 0.0;
    for name in INPUTS {
        let g = Golden::open(&root, name);
        let nb = golden_nbhds(&g);
        let (p, n) = (g.cloud.positions(), g.cloud.normals());
        let t = nb.k;
        // features and σ in the reference's arithmetic
        let mut feat = vec![0f32; p.len() * t * 6];
        let mut sig = Vec::new();
        let mut scr = FeatureScratch::default();
        for i in 0..p.len() {
            sig.push(neighborhood_features(p, n, nb.row(i), FeatureOrder::Reference, &mut scr, &mut feat[i * t * 6..(i + 1) * t * 6]));
        }
        let gs = g.arr("net_sigma").data;
        let sig_rel = sig.iter().zip(&gs).map(|(a, b)| ((a - b) / b).abs()).fold(0.0, f64::max);
        assert!(sig_rel <= 2.0 * f32::EPSILON as f64, "{name}: σ relative difference {sig_rel:e}");
        let raw = forward_features(&model, &feat, t, &opts).unwrap();
        let graw = g.arr("net_raw_outputs").data;
        let raw_d = raw.iter().zip(&graw).map(|(a, b)| (*a as f64 - b).abs()).fold(0.0, f64::max);
        let gco = g.coeffs();
        let mut coef_d: f64 = 0.0;
        for i in 0..p.len() {
            let y = [raw[6 * i], raw[6 * i + 1], raw[6 * i + 2], raw[6 * i + 3], raw[6 * i + 4], raw[6 * i + 5]];
            let a = scale_outputs(y, sig[i], FeatureOrder::Reference);
            for c in 0..6 {
                let (x, r) = (a[c] as f64, gco[i][c] as f64);
                coef_d = coef_d.max((x - r).abs() / (1.0 + r.abs()));
            }
        }
        // our default arithmetic (difference-first in f64) with our own neighborhoods
        let tree = KdTree::new(p, KdTreeOptions::default()).unwrap();
        let ours = fitting_neighborhoods(&tree, p, 64);
        let (co, _, _) = predict_coefficients(&model, &g.cloud, &ours, &opts).unwrap();
        let mut coef_default: f64 = 0.0;
        for i in 0..p.len() {
            for c in 0..6 {
                let (x, r) = (co[i][c] as f64, gco[i][c] as f64);
                coef_default = coef_default.max((x - r).abs() / (1.0 + r.abs()));
            }
        }
        eprintln!(
            "network {name:>12}: σ rel {sig_rel:.1e}, raw max |Δ| {raw_d:.2e}, coeff max |Δ|/(1+|a|) {coef_d:.2e} (reference order), {coef_default:.2e} (default order, our kNN)"
        );
        worst_raw = worst_raw.max(raw_d);
        worst_coef = worst_coef.max(coef_d);
    }
    assert!(worst_raw <= 2e-5, "raw outputs differ by {worst_raw:e}");
    assert!(worst_coef <= 2e-5, "coefficients differ by {worst_coef:e}");
}

#[test]
fn fit_matches_reference() {
    let Some(root) = golden_root() else { return };
    for name in INPUTS {
        let g = Golden::open(&root, name);
        let co = g.coeffs();
        let (tori, rep) = fit_tori(&g.cloud, &co, &FitOptions::default()).unwrap();
        let refs = g.tori();
        let mut worst = [0f64; 4];
        for (a, b) in tori.tori.iter().zip(&refs) {
            let rel = |x: f64, y: f64| (x - y).abs() / y.abs().max(1.0);
            for c in 0..3 {
                worst[0] = worst[0].max(rel(a.center[c], b.center[c]));
                worst[1] = worst[1].max((a.axis[c] - b.axis[c]).abs());
            }
            worst[2] = worst[2].max(rel(a.major, b.major));
            worst[3] = worst[3].max(rel(a.minor, b.minor));
        }
        eprintln!(
            "fit {name:>12}: center {:.1e}, axis {:.1e}, R {:.1e}, r {:.1e} (relative to max(1,|x|)); clamps disc {} eig {} κ {}",
            worst[0], worst[1], worst[2], worst[3], rep.clamped_disc, rep.clamped_eigvec, rep.clamped_curvature
        );
        assert!(worst.iter().all(|&w| w <= 1e-6), "{name}: tori differ: {worst:?}");
    }
}

fn eval(field: &TorusField, xs: &[[f64; 3]]) -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut d = vec![0.0; xs.len()];
    let mut gr = vec![[0.0; 3]; xs.len()];
    field.evaluate(xs, EvalOut { distance: Some(&mut d), gradient: Some(&mut gr), laplacian: None }).unwrap();
    (d, gr)
}

fn flat3(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flat_map(|x| x.iter().copied()).collect()
}

#[test]
fn sdf_matches_reference_given_tori() {
    let Some(root) = golden_root() else { return };
    for name in INPUTS {
        let g = Golden::open(&root, name);
        let tori = Tori { tori: g.tori(), coefficients: g.coeffs(), sigma: vec![f64::NAN; g.cloud.len()] };
        let (grid, ids) = grid_subset(&g);
        let qs = g.queries();
        for (tag, nbh) in [("accel", Neighborhood::Knn { k: 32 }), ("noaccel", Neighborhood::All)] {
            let field = TorusField::new(&g.cloud, tori.clone(), EvalOptions { neighborhood: nbh, ..Default::default() }).unwrap();
            let (d, gr) = eval(&field, &grid);
            let gd = g.arr(&format!("grid_sdf_{tag}")).data;
            let gg = g.arr(&format!("grid_grad_{tag}")).data;
            let gd: Vec<f64> = ids.iter().map(|&i| gd[i]).collect();
            let gg: Vec<f64> = ids.iter().flat_map(|&i| [gg[3 * i], gg[3 * i + 1], gg[3 * i + 2]]).collect();
            let (sd_max, sd_med, _) = diff_stats(&d, &gd);
            let (gr_max, _, _) = diff_stats(&flat3(&gr), &gg);
            let (qd, qg) = eval(&field, &qs);
            let (qd_max, _, _) = diff_stats(&qd, &g.arr(&format!("query_sdf_{tag}")).data);
            let (qg_max, _, _) = diff_stats(&flat3(&qg), &g.arr(&format!("query_grad_{tag}")).data);
            eprintln!(
                "sdf {name:>12} {tag:>7}: grid |Δφ| max {sd_max:.1e} median {sd_med:.1e}, |Δ∇φ| max {gr_max:.1e}; queries |Δφ| {qd_max:.1e}, |Δ∇φ| {qg_max:.1e}"
            );
            assert!(sd_max <= 1e-9 && qd_max <= 1e-9, "{name} {tag}: SDF differs by {sd_max:e} / {qd_max:e}");
            assert!(gr_max <= 1e-7 && qg_max <= 1e-7, "{name} {tag}: gradient differs by {gr_max:e} / {qg_max:e}");
        }
    }
}

#[test]
fn end_to_end_matches_reference() {
    let Some(root) = golden_root() else { return };
    let model = Model::reference().unwrap();
    for order in [FeatureOrder::Reference, FeatureOrder::DifferenceFirst] {
        for name in INPUTS {
            let g = Golden::open(&root, name);
            let mut opts = BuildOptions::default();
            opts.predict.feature_order = order;
            let (field, _rep) = PointsAsTori::build(&g.cloud, &model, &opts).unwrap();
            let (grid, ids) = grid_subset(&g);
            let (d, gr) = eval(&field, &grid);
            let gd = g.arr("grid_sdf_accel").data;
            let gg = g.arr("grid_grad_accel").data;
            let gd: Vec<f64> = ids.iter().map(|&i| gd[i]).collect();
            let gg: Vec<f64> = ids.iter().flat_map(|&i| [gg[3 * i], gg[3 * i + 1], gg[3 * i + 2]]).collect();
            let (mx, med, p999) = diff_stats(&d, &gd);
            let (gmx, gmed, _) = diff_stats(&flat3(&gr), &gg);
            eprintln!(
                "e2e {order:?} {name:>12}: grid |Δφ| max {mx:.1e} median {med:.1e} p99.9 {p999:.1e}; |Δ∇φ| max {gmx:.1e} median {gmed:.1e}"
            );
            assert!(med <= 1e-6 && p999 <= 1e-4, "{name} ({order:?}): SDF median {med:e}, p99.9 {p999:e}");
        }
    }
}

/// Report only: what [`SignRule::Geometric`] changes against the reference's rule (saddles and
/// exact zeros), on the golden coefficients; error against the exact SDF of the two synthetic
/// clouds.
#[test]
fn geometric_sign_rule_effect() {
    use ddg_points_as_tori::SignRule;
    let Some(root) = golden_root() else { return };
    for name in INPUTS {
        let g = Golden::open(&root, name);
        let co = g.coeffs();
        let (grid, _) = grid_subset(&g);
        let mut res = Vec::new();
        let mut saddles = 0;
        for rule in [SignRule::Reference, SignRule::Geometric] {
            let (tori, rep) = fit_tori(&g.cloud, &co, &FitOptions { sign_rule: rule, ..Default::default() }).unwrap();
            saddles = rep.saddle;
            let field = TorusField::new(&g.cloud, tori, EvalOptions::default()).unwrap();
            res.push(eval(&field, &grid).0);
        }
        let (mx, med, p999) = diff_stats(&res[0], &res[1]);
        let exact: Option<Box<dyn Fn([f64; 3]) -> f64>> = match name {
            "synth_sphere" => Some(Box::new(|x: [f64; 3]| (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt() - 0.5)),
            "synth_torus" => Some(Box::new(|x: [f64; 3]| ((x[0].hypot(x[1]) - 0.35).hypot(x[2])) - 0.15)),
            _ => None,
        };
        let err = exact.map(|f| {
            let truth: Vec<f64> = grid.iter().map(|&x| f(x)).collect();
            let mean = |v: &[f64]| v.iter().zip(&truth).map(|(a, b)| (a - b).abs()).sum::<f64>() / v.len() as f64;
            format!(", mean |φ − exact| reference {:.3e}, geometric {:.3e}", mean(&res[0]), mean(&res[1]))
        });
        eprintln!(
            "sign rule {name:>12}: saddles {saddles:>3}, grid |φ_ref − φ_geo| max {mx:.1e} median {med:.1e} p99.9 {p999:.1e}{}",
            err.unwrap_or_default()
        );
    }
}
