//! Training core tests (WP7, SPEC.md §10, §15.2). Run with `--features train`.
//!
//! The ignored tests need the training data and write to `runs_tmp/pat-wp7/` (see NOTES-WP7.md).
#![cfg(feature = "train")]

use burn::tensor::DType;
use ddg_points_as_tori::fit::{EigenRule, FitOptions, SignRule, fit_torus};
use ddg_points_as_tori::network::burn::cpu_device;
use ddg_points_as_tori::train::data::TrainShape;
use ddg_points_as_tori::train::loss::{CloudFit, LossOptions, build_plan, host_blend, loss_and_grad, tensor_tori};
use ddg_points_as_tori::train::net::{NetParams, PT};
use ddg_points_as_tori::train::npz::{Npz, NpzArray, NpzData, write_npz};
use ddg_points_as_tori::train::optim::{AdamOptions, AdamTrainer, Schedule, learning_rate, warmup_cosine};
use ddg_points_as_tori::train::queries::QueryOptions;
use ddg_points_as_tori::train::step::{StepOptions, prepare_shape, step};
use ddg_points_as_tori::train::tfit::{Ctx, TensorFitOptions};
use ddg_points_as_tori::train::truth::{GwnPrecision, MeshTruth, point_triangle_d2};
use ddg_points_as_tori::{ModelConfig, Torus};
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::{RngCore, SeedableRng};

fn u(r: &mut ChaCha8Rng) -> f64 {
    (r.next_u64() >> 11) as f64 / (1u64 << 53) as f64
}

fn ctx64() -> Ctx {
    Ctx { dev: cpu_device(), dtype: DType::F64 }
}

fn ctx32() -> Ctx {
    Ctx { dev: cpu_device(), dtype: DType::F32 }
}

/// An icosphere of radius `rad` (outward faces), `sub` subdivisions.
fn icosphere(sub: usize, rad: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0], [1.0, t, 0.0], [-1.0, -t, 0.0], [1.0, -t, 0.0],
        [0.0, -1.0, t], [0.0, 1.0, t], [0.0, -1.0, -t], [0.0, 1.0, -t],
        [t, 0.0, -1.0], [t, 0.0, 1.0], [-t, 0.0, -1.0], [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[u32; 3]> = vec![
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ];
    for _ in 0..sub {
        let mut mid = std::collections::HashMap::new();
        let mut nf = Vec::new();
        for tri in &f {
            let mut m = [0u32; 3];
            for e in 0..3 {
                let (a, b) = (tri[e].min(tri[(e + 1) % 3]), tri[e].max(tri[(e + 1) % 3]));
                m[e] = *mid.entry((a, b)).or_insert_with(|| {
                    let (pa, pb) = (v[a as usize], v[b as usize]);
                    v.push([(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0]);
                    (v.len() - 1) as u32
                });
            }
            nf.extend([[tri[0], m[0], m[2]], [tri[1], m[1], m[0]], [tri[2], m[2], m[1]], [m[0], m[1], m[2]]]);
        }
        f = nf;
    }
    for p in v.iter_mut() {
        let l = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        *p = [rad * p[0] / l, rad * p[1] / l, rad * p[2] / l];
    }
    (v, f)
}

/// A training shape: `n` points on an ellipsoid-ish blob (the icosphere stretched), normals from the
/// stretched sphere, `k` nearest neighbors (self first), the stretched mesh.
fn synth_shape(n: usize, k: usize, seed: u64) -> TrainShape {
    let (mut v, f) = icosphere(3, 0.5);
    let st = [1.0, 0.8, 0.6];
    for p in v.iter_mut() {
        *p = [p[0] * st[0], p[1] * st[1], p[2] * st[2]];
    }
    let mut r = ChaCha8Rng::seed_from_u64(seed);
    let mut pos = Vec::new();
    let mut nrm = Vec::new();
    for _ in 0..n {
        let (a, b) = (u(&mut r), u(&mut r));
        let z = 2.0 * a - 1.0;
        let ph = 2.0 * std::f64::consts::PI * b;
        let s = (1.0 - z * z).sqrt();
        let d = [s * ph.cos(), s * ph.sin(), z];
        let p = [0.5 * d[0] * st[0], 0.5 * d[1] * st[1], 0.5 * d[2] * st[2]];
        let g = [p[0] / (st[0] * st[0]), p[1] / (st[1] * st[1]), p[2] / (st[2] * st[2])];
        let l = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        pos.push(p);
        nrm.push([g[0] / l, g[1] / l, g[2] / l]);
    }
    let mut nb = vec![0i32; n * k];
    for i in 0..n {
        let mut d: Vec<(f64, usize)> = (0..n)
            .map(|j| {
                let e = [pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], pos[j][2] - pos[i][2]];
                (e[0] * e[0] + e[1] * e[1] + e[2] * e[2], j)
            })
            .collect();
        d.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        d.retain(|x| x.1 != i);
        nb[i * k] = i as i32;
        for c in 1..k {
            nb[i * k + c] = d[c - 1].1 as i32;
        }
    }
    TrainShape { source: ("synth".into(), 0), rows: n, n, k, positions: pos, normals: nrm, neighbors: nb, vertices: v, faces: f }
}

fn small_cfg() -> ModelConfig {
    ModelConfig { embed: 16, layers: 2, heads: 2, mlp: 32, k: 16, n_freq: 0 }
}

/// Small random parameters (head kernel non-zero so that every parameter matters).
fn small_params(ctx: &Ctx, seed: u64) -> NetParams {
    let cfg = small_cfg();
    let p = NetParams::init(&cfg, seed, &ctx32());
    let mut w = p.to_weights();
    let mut r = ChaCha8Rng::seed_from_u64(seed + 1);
    for x in w.head.w.iter_mut() {
        *x = (0.05 * (2.0 * u(&mut r) - 1.0)) as f32;
    }
    for x in w.head.b.iter_mut() {
        *x += (0.05 * (2.0 * u(&mut r) - 1.0)) as f32;
    }
    for b in w.blocks.iter_mut() {
        for x in b.norm1.bias.iter_mut().chain(b.ffn1.b.iter_mut()) {
            *x = (0.1 * (2.0 * u(&mut r) - 1.0)) as f32;
        }
    }
    NetParams::from_weights(&cfg, &w, ctx)
}

#[test]
fn npz_round_trip_and_errors() {
    let dir = std::env::temp_dir().join(format!("pat-wp7-npz-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("a.npz");
    let a = NpzArray { shape: vec![2, 3], data: NpzData::F64(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.5]) };
    let b = NpzArray { shape: vec![4], data: NpzData::I32(vec![-1, 0, 7, 9]) };
    let c = NpzArray { shape: vec![3], data: NpzData::Bool(vec![true, false, true]) };
    write_npz(&p, &[("x", &a), ("y", &b), ("z", &c)]).unwrap();
    let mut z = Npz::open(&p).unwrap();
    assert_eq!(z.read("x").unwrap(), a);
    assert_eq!(z.read("y").unwrap(), b);
    assert_eq!(z.read_bool("z").unwrap().1, vec![true, false, true]);
    let e = z.read_i32("x").unwrap_err().to_string();
    assert!(e.contains("member x") && e.contains("expected <i4"), "{e}");
    assert!(z.read("nope").unwrap_err().to_string().contains("member nope"));
    assert!(Npz::open(&dir.join("missing.npz")).is_err());
    std::fs::write(dir.join("bad.npz"), b"not a zip").unwrap();
    assert!(Npz::open(&dir.join("bad.npz")).is_err());
    let _ = std::fs::remove_dir_all(&dir);
}

fn write_shape_npz(path: &std::path::Path, s: &TrainShape, normal_scale: f64, break_mask: bool, break_self: bool) {
    let p = s.rows;
    let mut masks = vec![true; p];
    if break_mask {
        masks[3] = false;
    }
    let mut nb = s.neighbors.clone();
    if break_self {
        nb[5 * s.k] = 0;
    }
    let f = s.faces.len();
    let flat = |v: &[[f64; 3]], sc: f64| v.iter().flat_map(|x| x.map(|y| y * sc)).collect::<Vec<f64>>();
    let pos = NpzArray { shape: vec![1, p, 3], data: NpzData::F64(flat(&s.positions, 1.0)) };
    let nrm = NpzArray { shape: vec![1, p, 3], data: NpzData::F64(flat(&s.normals, normal_scale)) };
    let nbr = NpzArray { shape: vec![1, p, s.k], data: NpzData::I32(nb) };
    let msk = NpzArray { shape: vec![1, p], data: NpzData::Bool(masks) };
    let ver = NpzArray { shape: vec![1, s.vertices.len(), 3], data: NpzData::F64(flat(&s.vertices, 1.0)) };
    let fac = NpzArray { shape: vec![1, f, 3], data: NpzData::I32(s.faces.iter().flat_map(|t| t.map(|j| j as i32)).collect()) };
    let fm = NpzArray { shape: vec![1, f], data: NpzData::Bool(vec![true; f]) };
    write_npz(path, &[("positions", &pos), ("normals", &nrm), ("neighbors", &nbr), ("masks", &msk), ("mesh_vertices", &ver), ("mesh_faces", &fac), ("face_masks", &fm)]).unwrap();
}

#[test]
fn data_validation_and_normals() {
    use ddg_points_as_tori::train::data::{DataOptions, load_file};
    let dir = std::env::temp_dir().join(format!("pat-wp7-data-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let s = synth_shape(80, 16, 3);
    let p = dir.join("ok.npz");
    write_shape_npz(&p, &s, 2.0, false, false);
    let (sh, rep) = load_file(&p, &DataOptions::default()).unwrap();
    assert_eq!(rep.nonunit_normals, 80);
    assert!((rep.normal_length_range.0 - 2.0).abs() < 1e-9);
    let l = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    assert!((l(sh[0].normals[7]) - 1.0).abs() < 1e-12);
    let (raw, _) = load_file(&p, &DataOptions { raw_normals: true, ..Default::default() }).unwrap();
    assert!((l(raw[0].normals[7]) - 2.0).abs() < 1e-9);
    assert_eq!(sh[0].positions, s.positions);
    let pm = dir.join("mask.npz");
    write_shape_npz(&pm, &s, 1.0, true, false);
    let e = load_file(&pm, &DataOptions::default()).unwrap_err().to_string();
    assert!(e.contains("not a prefix"), "{e}");
    let ps = dir.join("self.npz");
    write_shape_npz(&ps, &s, 1.0, false, true);
    let e = load_file(&ps, &DataOptions::default()).unwrap_err().to_string();
    assert!(e.contains("row 5 starts with 0"), "{e}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ground_truth_distance_and_sign() {
    let (v, f) = icosphere(4, 0.5);
    let m = MeshTruth::new(&v, &f);
    let mut r = ChaCha8Rng::seed_from_u64(9);
    for _ in 0..300 {
        let q = [2.0 * u(&mut r) - 1.0, 2.0 * u(&mut r) - 1.0, 2.0 * u(&mut r) - 1.0];
        let brute = f.iter().map(|t| point_triangle_d2(q, v[t[0] as usize], v[t[1] as usize], v[t[2] as usize])).fold(f64::INFINITY, f64::min).sqrt();
        let d = m.distance(q);
        assert!((d - brute).abs() < 1e-14, "{d} {brute}");
        let sd = m.signed_distance(q, GwnPrecision::F64);
        let rq = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        assert_eq!(sd < 0.0, rq < 0.5 - 1e-3 || (rq < 0.5 && sd < 0.0), "q {q:?} sd {sd}");
        assert!((sd - (rq - 0.5)).abs() < 5e-3, "sd {sd} vs sphere {}", rq - 0.5);
        if (rq - 0.5).abs() > 1e-2 {
            assert_eq!(m.signed_distance(q, GwnPrecision::ReferenceF32).signum(), sd.signum());
        }
    }
    // closest-point regions against dense sampling of one triangle
    let (a, b, c) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 0.7, 0.0]);
    for _ in 0..200 {
        let q = [3.0 * u(&mut r) - 1.0, 3.0 * u(&mut r) - 1.0, u(&mut r) - 0.5];
        let mut best = f64::INFINITY;
        for i in 0..=200 {
            for j in 0..=(200 - i) {
                let (s, t) = (i as f64 / 200.0, j as f64 / 200.0);
                let p = [a[0] + s * (b[0] - a[0]) + t * (c[0] - a[0]), a[1] + s * (b[1] - a[1]) + t * (c[1] - a[1]), 0.0];
                best = best.min((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2));
            }
        }
        let d = point_triangle_d2(q, a, b, c);
        assert!(d <= best + 1e-12 && best.sqrt() - d.sqrt() < 1e-2, "{d} {best}");
    }
}

/// Random coefficient rows in the trained range (curvatures of a few 1/σ, slopes below 1).
fn random_coefficients(r: &mut ChaCha8Rng, n: usize) -> Vec<f64> {
    (0..n)
        .flat_map(|_| {
            let s = |r: &mut ChaCha8Rng, w: f64| w * (2.0 * u(r) - 1.0);
            [s(r, 0.02), s(r, 0.5), s(r, 0.5), s(r, 3.0), s(r, 3.0), s(r, 3.0)]
        })
        .collect()
}

#[test]
fn tensor_fit_matches_the_f64_fit() {
    let s = synth_shape(200, 16, 5);
    let ctx = ctx64();
    let mut r = ChaCha8Rng::seed_from_u64(11);
    let a = random_coefficients(&mut r, 200);
    let cloud = CloudFit { p: &s.positions, n: &s.normals };
    for sign_rule in [SignRule::Reference, SignRule::Geometric] {
        for eigen_rule in [EigenRule::Reference, EigenRule::Robust] {
            let tt = tensor_tori(&ctx, &cloud, &a, &TensorFitOptions { sign_rule, eigen_rule, ..Default::default() });
            let fo = FitOptions { sign_rule, eigen_rule, ..Default::default() };
            for i in 0..200 {
                let ai = [a[6 * i], a[6 * i + 1], a[6 * i + 2], a[6 * i + 3], a[6 * i + 4], a[6 * i + 5]];
                let (t, _) = fit_torus(s.positions[i], s.normals[i], ai, &fo);
                let rel = |x: f64, y: f64| (x - y).abs() / (1.0 + y.abs());
                let e = (0..3).map(|k| rel(tt[i].center[k], t.center[k])).fold(0.0, f64::max).max(rel(tt[i].major, t.major)).max(rel(tt[i].minor, t.minor));
                let dot = (0..3).map(|k| tt[i].axis[k] * t.axis[k]).sum::<f64>().abs();
                assert!(e < 1e-11 && (1.0 - dot) < 1e-12, "{sign_rule:?} {eigen_rule:?} point {i}: {:?} vs {:?}", tt[i], t);
            }
        }
    }
}

#[test]
fn robust_eigen_rule() {
    // in f64 the two rules agree away from v2 = 0
    let s = synth_shape(50, 16, 6);
    let mut r = ChaCha8Rng::seed_from_u64(12);
    let a = random_coefficients(&mut r, 50);
    for i in 0..50 {
        let ai = [a[6 * i], a[6 * i + 1], a[6 * i + 2], a[6 * i + 3], a[6 * i + 4], a[6 * i + 5]];
        let (t0, _) = fit_torus(s.positions[i], s.normals[i], ai, &FitOptions::default());
        let (t1, d1) = fit_torus(s.positions[i], s.normals[i], ai, &FitOptions { eigen_rule: EigenRule::Robust, ..Default::default() });
        assert!((t0.major - t1.major).abs() < 1e-10 * (1.0 + t0.major.abs()) && (t0.minor - t1.minor).abs() < 1e-10 * (1.0 + t0.minor.abs()));
        assert!(!d1.clamps.eigen_switched);
    }
    // v2 = 0 exactly (a11 = 0, a01 = a10 = 0): the reference collapses to a sphere, robust does not
    let ai = [0.0, 0.0, 0.0, 0.0, -0.2, -0.5];
    let (_, dr) = fit_torus([0.0; 3], [0.0, 0.0, 1.0], ai, &FitOptions::default());
    let (_, db) = fit_torus([0.0; 3], [0.0, 0.0, 1.0], ai, &FitOptions { eigen_rule: EigenRule::Robust, ..Default::default() });
    assert!((dr.kappa_min - dr.kappa_max).abs() < 1e-9, "reference: {} {}", dr.kappa_min, dr.kappa_max);
    assert!((db.kappa_min + 0.4).abs() < 1e-12 && (db.kappa_max + 1.0).abs() < 1e-12, "robust: {} {}", db.kappa_min, db.kappa_max);
    assert!(db.clamps.eigen_switched);
    // f32 tensor fit: v2 = a11 → 0 with v1 = a20 − a02 = −0.3 (ev₊ cancels); the robust curvatures
    // match the f64 fit, the reference's degrade once |v2/v1| reaches the f32 round-off
    let es = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 0.0];
    let cloud_p = vec![[0.0; 3]; es.len()];
    let cloud_n = vec![[0.0, 0.0, 1.0]; es.len()];
    let rows: Vec<f64> = es.iter().flat_map(|&e| [0.0, 0.0, 0.0, e, -0.2, -0.5]).collect();
    let cloud = CloudFit { p: &cloud_p, n: &cloud_n };
    let t_ref = tensor_tori(&ctx32(), &cloud, &rows, &TensorFitOptions::default());
    let t_rob = tensor_tori(&ctx32(), &cloud, &rows, &TensorFitOptions { eigen_rule: EigenRule::Robust, ..Default::default() });
    let mut worst_rob = 0.0f64;
    for (i, e) in es.iter().enumerate() {
        let a6 = [0.0, 0.0, 0.0, *e, -0.2, -0.5];
        let (t64, _) = fit_torus(cloud_p[i], cloud_n[i], a6, &FitOptions { eigen_rule: EigenRule::Robust, ..Default::default() });
        let err = |t: &Torus| ((t.major - t64.major).abs() + (t.minor - t64.minor).abs()) / (1.0 + t64.minor.abs());
        eprintln!("a11 = {e:.0e}: f32 reference rule error {:.2e}, robust {:.2e}", err(&t_ref[i]), err(&t_rob[i]));
        worst_rob = worst_rob.max(err(&t_rob[i]));
        if *e == 0.0 {
            assert!(err(&t_ref[i]) > 0.1, "the reference rule collapses at v2 = 0");
        }
    }
    assert!(worst_rob < 1e-5, "{worst_rob}");
}

#[test]
fn robust_rule_is_continuous_across_the_switch() {
    // v1 = a20 − a02 (a01 = a10 = 0) crosses 0; derivative of the loss in a20 from autodiff on both
    // sides equals the central difference, and the tori are continuous
    let s = synth_shape(120, 16, 8);
    let ctx = ctx64();
    let mut r = ChaCha8Rng::seed_from_u64(21);
    let mut a = random_coefficients(&mut r, 120);
    let i0 = 0;
    a[6 * i0 + 1] = 0.0;
    a[6 * i0 + 2] = 0.0;
    a[6 * i0 + 3] = 0.3;
    a[6 * i0 + 4] = -1.0;
    let o = LossOptions { fit: TensorFitOptions { eigen_rule: EigenRule::Robust, ..Default::default() }, cutoff: f64::INFINITY, ..Default::default() };
    let cloud = CloudFit { p: &s.positions, n: &s.normals };
    let mesh = MeshTruth::new(&s.vertices, &s.faces);
    let mut q = Vec::new();
    for k in 0..60 {
        let p = s.positions[i0];
        let n = s.normals[i0];
        let t = 0.002 * (k as f64 - 30.0);
        q.push([p[0] + t * n[0] + 0.01 * (k % 7) as f64, p[1] + t * n[1], p[2] + t * n[2] - 0.004 * (k % 5) as f64]);
    }
    let truth: Vec<f64> = q.iter().map(|&x| mesh.signed_distance(x, GwnPrecision::F64)).collect();
    let plan = build_plan(&q, &truth, &s.positions, &o, None).unwrap();
    let mut grads = Vec::new();
    for a20 in [-1.0 - 1e-4, -1.0 + 1e-4] {
        let mut b = a.clone();
        b[6 * i0 + 5] = a20;
        let (_, g) = loss_and_grad(&ctx, &cloud, &b, &plan, &o, true).unwrap();
        let g = g.unwrap()[6 * i0 + 5];
        let h = 1e-7;
        let mut bp = b.clone();
        bp[6 * i0 + 5] += h;
        let mut bm = b.clone();
        bm[6 * i0 + 5] -= h;
        let lp = loss_and_grad(&ctx, &cloud, &bp, &plan, &o, false).unwrap().0.loss;
        let lm = loss_and_grad(&ctx, &cloud, &bm, &plan, &o, false).unwrap().0.loss;
        let fd = (lp - lm) / (2.0 * h);
        assert!((g - fd).abs() < 1e-5 * (1.0 + fd.abs()), "a20 {a20}: autodiff {g} fd {fd}");
        grads.push(g);
    }
    assert!((grads[0] - grads[1]).abs() < 1e-2 * (1.0 + grads[0].abs()), "{grads:?}");
}

/// A shape and a loss setting for the finite-difference tests (f64, no truncation).
fn fd_setup(n: usize, seed: u64) -> (TrainShape, StepOptions) {
    let s = synth_shape(n, 16, seed);
    let o = StepOptions {
        loss: LossOptions { cutoff: f64::INFINITY, pair_chunk: 1500, ..Default::default() },
        queries: QueryOptions { per_kind: 4, ..Default::default() },
        network_chunk: 3,
        groups: 5,
        ..Default::default()
    };
    (s, o)
}

#[test]
fn tensor_blend_matches_host_blend() {
    let s = synth_shape(150, 16, 9);
    let ctx = ctx64();
    let mut r = ChaCha8Rng::seed_from_u64(13);
    let a = random_coefficients(&mut r, 150);
    let o = LossOptions { cutoff: f64::INFINITY, pair_chunk: 2000, ..Default::default() };
    let mesh = MeshTruth::new(&s.vertices, &s.faces);
    let q: Vec<[f64; 3]> = (0..200).map(|_| [u(&mut r) - 0.5, 0.8 * u(&mut r) - 0.4, 0.6 * u(&mut r) - 0.3]).collect();
    let truth: Vec<f64> = q.iter().map(|&x| mesh.signed_distance(x, GwnPrecision::F64)).collect();
    let plan = build_plan(&q, &truth, &s.positions, &o, None).unwrap();
    let tori: Vec<Torus> = (0..150).map(|i| fit_torus(s.positions[i], s.normals[i], std::array::from_fn(|c| a[6 * i + c]), &FitOptions::default()).0).collect();
    let host = host_blend(&tori, &plan);
    let cloud = CloudFit { p: &s.positions, n: &s.normals };
    let (v, _) = loss_and_grad(&ctx, &cloud, &a, &plan, &o, false).unwrap();
    let (mut sd, mut se) = (0.0, 0.0);
    for (j, (phi, g)) in host.iter().enumerate() {
        sd += (phi - plan.truth[j]).abs();
        se += (1.0 - (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()).abs();
    }
    let n = host.len() as f64;
    assert!((v.distance - sd / n).abs() < 1e-11 && (v.eikonal - se / n).abs() < 1e-10, "{v:?} vs {} {}", sd / n, se / n);
    // truncation at T = 40 changes the loss by round-off only
    let plan40 = build_plan(&q, &truth, &s.positions, &LossOptions { cutoff: 40.0, ..o }, None).unwrap();
    let (v40, _) = loss_and_grad(&ctx, &cloud, &a, &plan40, &o, false).unwrap();
    eprintln!("pairs: all {}, T = 40: {} (max dropped mass {:.1e})", plan.pairs(), plan40.pairs(), plan40.max_dropped_mass);
    assert!((v40.loss - v.loss).abs() < 1e-12, "{} {}", v40.loss, v.loss);
    assert!(plan40.pairs() < plan.pairs() / 2);
}

#[test]
fn fd_loss_wrt_coefficients() {
    let (s, _) = fd_setup(150, 10);
    let ctx = ctx64();
    let mut r = ChaCha8Rng::seed_from_u64(14);
    let a = random_coefficients(&mut r, 150);
    for fit in [TensorFitOptions::default(), TensorFitOptions { sign_rule: SignRule::Geometric, eigen_rule: EigenRule::Robust, ..Default::default() }] {
        let o = LossOptions { fit, cutoff: f64::INFINITY, pair_chunk: 3000, ..Default::default() };
        let mesh = MeshTruth::new(&s.vertices, &s.faces);
        let q: Vec<[f64; 3]> = (0..240).map(|_| [u(&mut r) - 0.5, 0.8 * u(&mut r) - 0.4, 0.6 * u(&mut r) - 0.3]).collect();
        let truth: Vec<f64> = q.iter().map(|&x| mesh.signed_distance(x, GwnPrecision::F64)).collect();
        let plan = build_plan(&q, &truth, &s.positions, &o, None).unwrap();
        let cloud = CloudFit { p: &s.positions, n: &s.normals };
        let (_, g) = loss_and_grad(&ctx, &cloud, &a, &plan, &o, true).unwrap();
        let g = g.unwrap();
        let mut worst = 0.0f64;
        let mut checked = 0;
        for i in (0..150).step_by(7) {
            for c in 0..6 {
                let j = 6 * i + c;
                let h = 1e-6 * (1.0f64).max(a[j].abs());
                let mut ap = a.clone();
                ap[j] += h;
                let mut am = a.clone();
                am[j] -= h;
                let lp = loss_and_grad(&ctx, &cloud, &ap, &plan, &o, false).unwrap().0.loss;
                let lm = loss_and_grad(&ctx, &cloud, &am, &plan, &o, false).unwrap().0.loss;
                let fd = (lp - lm) / (2.0 * h);
                let scale = g.iter().map(|x| x.abs()).fold(0.0, f64::max);
                worst = worst.max((fd - g[j]).abs() / scale);
                checked += 1;
            }
        }
        eprintln!("∂L/∂a {fit:?}: {checked} entries, worst |fd − autodiff|/max|g| = {worst:.2e}");
        assert!(worst < 1e-5, "{worst}");
    }
}

#[test]
fn fd_blend_gradient_wrt_x() {
    // ∇φ of the framework path against central differences of φ with σ (hence λ) frozen
    let s = synth_shape(150, 16, 15);
    let ctx = ctx64();
    let mut r = ChaCha8Rng::seed_from_u64(16);
    let a = random_coefficients(&mut r, 150);
    let o = LossOptions { cutoff: f64::INFINITY, ..Default::default() };
    let cloud = CloudFit { p: &s.positions, n: &s.normals };
    let tori = tensor_tori(&ctx, &cloud, &a, &o.fit);
    let mut worst = 0.0f64;
    for _ in 0..60 {
        let q = [u(&mut r) - 0.5, 0.8 * u(&mut r) - 0.4, 0.6 * u(&mut r) - 0.3];
        let plan0 = build_plan(&[q], &[0.0], &s.positions, &o, None).unwrap();
        let sig = plan0.sigma[0];
        let g = host_blend(&tori, &plan0)[0].1;
        for k in 0..3 {
            let h = 1e-6;
            let mut qp = q;
            qp[k] += h;
            let mut qm = q;
            qm[k] -= h;
            let pp = build_plan(&[qp], &[0.0], &s.positions, &o, Some(&[sig])).unwrap();
            let pm = build_plan(&[qm], &[0.0], &s.positions, &o, Some(&[sig])).unwrap();
            let fd = (host_blend(&tori, &pp)[0].0 - host_blend(&tori, &pm)[0].0) / (2.0 * h);
            worst = worst.max((fd - g[k]).abs() / (1.0 + g[k].abs()));
        }
    }
    eprintln!("∇φ vs FD (frozen σ): worst relative {worst:.2e}");
    assert!(worst < 1e-6, "{worst}");
}

#[test]
fn fd_loss_wrt_network_parameters() {
    let (s, o) = fd_setup(120, 17);
    let ctx = ctx64();
    let params = small_params(&ctx, 3);
    let prep = prepare_shape(&s, &o, (1, 2, 3)).unwrap();
    let out = step(&params, &ctx, &prep, 0, 60, &o, true).unwrap();
    let grads = out.grads.unwrap();
    let gmax = grads.iter().flat_map(|g| g.to_f64()).map(f64::abs).fold(0.0, f64::max);
    let mut worst = 0.0f64;
    let mut checked = 0;
    let mut r = ChaCha8Rng::seed_from_u64(18);
    for (ti, item) in params.items.iter().enumerate() {
        let vals = item.value.to_f64();
        let g = grads[ti].to_f64();
        for _ in 0..2 {
            let e = (r.next_u64() % vals.len() as u64) as usize;
            let h = 1e-6 * (1.0f64).max(vals[e].abs());
            let lfun = |d: f64| {
                let mut v = vals.clone();
                v[e] += d;
                let mut items = params.items.iter().map(|p| p.value.clone()).collect::<Vec<PT>>();
                items[ti] = item.value.like(v, &ctx);
                let p2 = params.with_values(items);
                step(&p2, &ctx, &prep, 0, 60, &o, false).unwrap().value.loss
            };
            let fd = (lfun(h) - lfun(-h)) / (2.0 * h);
            worst = worst.max((fd - g[e]).abs() / gmax);
            checked += 1;
        }
    }
    eprintln!("∂L/∂θ (2 layers, 2 heads, 16 tokens, f64): {checked} entries, worst |fd − autodiff|/max|g| = {worst:.2e}, loss {:.4}", out.value.loss);
    assert!(worst < 1e-5, "{worst}");
}

#[test]
fn two_pass_is_independent_of_chunks_and_groups() {
    let (s, o) = fd_setup(100, 19);
    let ctx = ctx64();
    let params = small_params(&ctx, 4);
    let prep = prepare_shape(&s, &o, (1, 2, 3)).unwrap();
    let a = step(&params, &ctx, &prep, 0, 100, &o, true).unwrap();
    let b = step(&params, &ctx, &prep, 0, 100, &StepOptions { network_chunk: 17, groups: 1, ..o }, true).unwrap();
    assert!((a.value.loss - b.value.loss).abs() < 1e-14);
    for (x, y) in a.grads.unwrap().iter().zip(b.grads.unwrap().iter()) {
        for (p, q) in x.to_f64().iter().zip(y.to_f64()) {
            assert!((p - q).abs() < 1e-12 * (1.0 + q.abs()), "{p} {q}");
        }
    }
}

#[test]
fn f32_step_is_bitwise_independent_of_threads() {
    let (s, o) = fd_setup(100, 20);
    let o = StepOptions { loss: LossOptions { cutoff: 40.0, ..o.loss }, ..o };
    let ctx = ctx32();
    let params = small_params(&ctx, 5);
    let run = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        pool.install(|| {
            let prep = prepare_shape(&s, &o, (1, 2, 3)).unwrap();
            let out = step(&params, &ctx, &prep, 0, 100, &o, true).unwrap();
            let g: Vec<u32> = out.grads.unwrap().iter().flat_map(|g| g.to_f32()).map(f32::to_bits).collect();
            (out.value.loss.to_bits(), g)
        })
    };
    let a = run(1);
    assert_eq!(a, run(1));
    assert_eq!(a, run(4));
    assert_eq!(a, run(13));
}

#[test]
fn initial_output_is_the_tangent_sphere() {
    let ctx = ctx32();
    let cfg = ModelConfig::default();
    let p = NetParams::init(&cfg, 42, &ctx);
    assert_eq!(p.count(), 1_588_358);
    let s = synth_shape(300, 64, 22);
    let o = StepOptions::default();
    let prep = prepare_shape(&s, &o, (0, 0, 0)).unwrap();
    let f = &prep.feats;
    let (net, _) = p.build(&ctx, false).unwrap();
    let y = ddg_points_as_tori::train::net::forward(&net, &ctx, f, 0, 10);
    let y = ddg_points_as_tori::train::tfit::to_f64(y);
    for r in 0..10 {
        assert_eq!(&y[6 * r..6 * r + 6], &[0.0, 0.0, 0.0, 0.0, -0.5, -0.5]);
    }
    // a = (0,0,0,0,−½/σ,−½/σ): the sphere of radius σ tangent at p, center p − σ n
    let i = 3;
    let sg = f.sigma[i];
    let (t, _) = fit_torus(s.positions[i], s.normals[i], [0.0, 0.0, 0.0, 0.0, -0.5 / sg, -0.5 / sg], &FitOptions::default());
    assert!(t.major.abs() < 1e-12 && (t.minor - sg).abs() < 1e-12, "{t:?} σ {sg}");
    for k in 0..3 {
        assert!((t.center[k] - (s.positions[i][k] - sg * s.normals[i][k])).abs() < 1e-12);
    }
    // the safetensors layout round-trips through the crate's loader
    let b = p.to_safetensors(&Default::default()).unwrap();
    let m = ddg_points_as_tori::Model::from_safetensors(&b).unwrap();
    assert_eq!(m.weights(), &p.to_weights());
}

#[test]
fn adam_and_schedules() {
    // the optax schedule
    assert_eq!(warmup_cosine(0, 1e-3, 0, 100_000_000), 1e-3);
    assert_eq!(warmup_cosine(0, 1e-5, 100, 100_000_000), 0.0);
    assert!((warmup_cosine(50, 1e-5, 100, 100_000_000) - 0.5e-5).abs() < 1e-20);
    assert!((learning_rate(Schedule::Code, 96_000, 1e-5, 100, 96_000) - 1e-5).abs() < 1e-5 * 1e-5);
    assert!((learning_rate(Schedule::Cosine, 100 + 96_000, 1e-5, 100, 96_000) - 1e-7).abs() < 1e-15);
    assert!((learning_rate(Schedule::Cosine, 100 + 48_000, 1e-5, 100, 96_000) - 0.505e-5).abs() < 1e-12);
    // burn's Adam per tensor equals optax's update m̂/(√v̂ + ε)
    let ctx = ctx64();
    let p = small_params(&ctx, 6);
    let mut tr = AdamTrainer::new(&p, AdamOptions::default());
    let g1: Vec<PT> = p.items.iter().map(|it| it.value.scale(0.5).add(&it.value.like(vec![0.01; it.value.len()], &ctx))).collect();
    let g2: Vec<PT> = p.items.iter().map(|it| it.value.scale(-0.3)).collect();
    let p1 = tr.update(&p, &g1, 1e-3).unwrap();
    let p2 = tr.update(&p1, &g2, 2e-3).unwrap();
    for (i, it) in p.items.iter().enumerate() {
        let x0 = it.value.to_f64();
        let (a1, a2) = (g1[i].to_f64(), g2[i].to_f64());
        let got = p2.items[i].value.to_f64();
        for e in 0..x0.len().min(20) {
            let (b1, b2, eps) = (0.9f64, 0.999f64, 1e-8f64);
            let m1 = (1.0 - b1) * a1[e];
            let v1 = (1.0 - b2) * a1[e] * a1[e];
            let x1 = x0[e] - 1e-3 * (m1 / (1.0 - b1)) / ((v1 / (1.0 - b2)).sqrt() + eps);
            let m2 = b1 * m1 + (1.0 - b1) * a2[e];
            let v2 = b2 * v1 + (1.0 - b2) * a2[e] * a2[e];
            let x2 = x1 - 2e-3 * (m2 / (1.0 - b1 * b1)) / ((v2 / (1.0 - b2 * b2)).sqrt() + eps);
            // burn keeps β and ε in f32 (as optax does in f32 training)
            assert!((got[e] - x2).abs() < 2e-5 * 3e-3, "{} {e}: {} vs {x2}", it.name, got[e]);
        }
    }
}

#[test]
fn curriculum_file() {
    use ddg_points_as_tori::train::Curriculum;
    let c = Curriculum::reference().unwrap();
    let names: Vec<&str> = c.phases.iter().map(|p| p.name.as_str()).collect();
    assert_eq!(names, ["pretrain", "phase0", "phase1", "phase2", "phase3", "phase4"]);
    assert_eq!(c.phases[0].batch_rows, 32);
    assert!(c.phases[0].from_scratch && !c.phases[1].from_scratch && c.phases[2].interleave);
    assert_eq!(c.phases.iter().map(|p| p.epochs).collect::<Vec<_>>(), [64, 24, 8, 8, 5, 2]);
    assert!(Curriculum::from_toml("[[phase]]\nname = 1").is_err());
    let dir = std::path::Path::new(DATA);
    if dir.is_dir() {
        // the step counts of NOTES-REFERENCE.md §9.3
        let mut total = 0u64;
        let want = [4_096u64, 96_000, 64_000, 119_344, 124_590, 89_828];
        for (p, w) in c.phases.iter().zip(want) {
            let files = ddg_points_as_tori::train::run::phase_files(dir, p).unwrap();
            let mut shapes = 0u64;
            for f in &files {
                let mut z = Npz::open(f).unwrap();
                shapes += z.read_bool("masks").unwrap().0[0] as u64;
            }
            let steps = shapes * 2048u64.div_ceil(p.batch_rows as u64) * p.epochs;
            assert_eq!(steps, w, "{}", p.name);
            total += steps;
        }
        assert_eq!(total, 497_858);
    }
}

const DATA: &str = "/var/data/discrete-geometry/PointsAsTori/PAT_training_data/PAT_training_data";

fn out_dir(name: &str) -> std::path::PathBuf {
    let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../runs_tmp/pat-wp7").join(name);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn peak_rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("VmHWM")).and_then(|l| l.split_whitespace().nth(1)).and_then(|v| v.parse::<f64>().ok()))
        .map(|kb| kb / 1024.0)
        .unwrap_or(0.0)
}

/// Our loss with the authors' weights on `ABC17_validation_3.npz` (raw normals and the reference's
/// feature order, as the reference computes it), and an export of the exact inputs for
/// `runs_tmp/pat-wp7/refloss/ref_loss.py` (the reference's `SDF_loss` in JAX).
#[test]
#[ignore]
fn reference_loss_export() {
    use ddg_points_as_tori::features::FeatureOrder;
    use ddg_points_as_tori::train::data::{DataOptions, load_file};
    let dir = out_dir("refloss");
    let file = std::path::Path::new(DATA).join("ABC17_validation_3.npz");
    let model = ddg_points_as_tori::Model::reference().unwrap();
    let ctx = ctx32();
    let params = NetParams::from_model(&model, &ctx);
    let mut json = String::from("[\n");
    for raw_normals in [true, false] {
        let (shapes, rep) = load_file(&file, &DataOptions { raw_normals, ..Default::default() }).unwrap();
        eprintln!("{} shapes, {} non-unit normals, lengths {:?}", rep.shapes, rep.nonunit_normals, rep.normal_length_range);
        for (si, s) in shapes.iter().enumerate() {
            let base = StepOptions { feature_order: FeatureOrder::Reference, ..Default::default() };
            let prep = prepare_shape(s, &base, (42, 0, 10_000 + si as u64)).unwrap();
            let prep_df = prepare_shape(s, &StepOptions { feature_order: FeatureOrder::DifferenceFirst, ..base }, (42, 0, 10_000 + si as u64)).unwrap();
            let mut variants: Vec<(&str, StepOptions)> = vec![("T40", base)];
            if raw_normals && si < 3 && std::env::var("PAT_TINF").is_ok() {
                variants.push(("Tinf", StepOptions { loss: LossOptions { cutoff: f64::INFINITY, ..base.loss }, ..base }));
            }
            variants.extend([

                ("robust", StepOptions { loss: LossOptions { fit: TensorFitOptions { eigen_rule: EigenRule::Robust, ..Default::default() }, ..base.loss }, ..base }),
                ("geometric", StepOptions { loss: LossOptions { fit: TensorFitOptions { sign_rule: SignRule::Geometric, ..Default::default() }, ..base.loss }, ..base }),
                ("diff_first", StepOptions { feature_order: FeatureOrder::DifferenceFirst, ..base }),
            ]);
            // saddle tori of the authors' network on this cloud (they are where Geometric differs)
            let cloud = ddg_points_as_tori::PointCloud::new(s.positions[..s.n].to_vec(), s.normals[..s.n].to_vec(), ddg_points_as_tori::NormalPolicy::Normalize).unwrap().0;
            let nb = ddg_points_as_tori::spatial::Neighborhoods { k: s.k, indices: s.neighbors[..s.n * s.k].iter().map(|&j| j as u32).collect() };
            let po = ddg_points_as_tori::PredictOptions { feature_order: FeatureOrder::Reference, ..Default::default() };
            let (co, _, _) = ddg_points_as_tori::predict_coefficients(&model, &cloud, &nb, &po).unwrap();
            let (_, fr) = ddg_points_as_tori::fit_tori(&cloud, &co, &FitOptions::default()).unwrap();
            let (_, fr32) = ddg_points_as_tori::fit_tori(&cloud, &co, &FitOptions { eigen_rule: EigenRule::Robust, ..Default::default() }).unwrap();
            eprintln!("shape {si}: {} saddle tori, {} where the robust rule changes the f64 fit", fr.saddle, fr32.eigen_switched);
            for (name, o) in variants {
                let pr = if name == "diff_first" { &prep_df } else { &prep };
                for (b, lo) in [0usize, 1024].iter().enumerate() {
                    let t = std::time::Instant::now();
                    let v = step(&params, &ctx, pr, *lo, lo + 1024, &o, false).unwrap().value;
                    json += &format!(
                        "{{\"raw_normals\": {raw_normals}, \"shape\": {si}, \"batch\": {b}, \"variant\": \"{name}\", \"distance\": {:e}, \"eikonal\": {:e}, \"queries\": {}, \"seconds\": {:.3}}},\n",
                        v.distance, v.eikonal, v.queries, t.elapsed().as_secs_f64()
                    );
                }
            }
            if raw_normals {
                // the inputs for the JAX script: the cloud as stored, our queries and ground truth
                let p = s.rows;
                let m = 120;
                let flat = |v: &[[f64; 3]]| v.iter().flatten().copied().collect::<Vec<f64>>();
                let pos = NpzArray { shape: vec![p, 3], data: NpzData::F64(flat(&s.positions)) };
                let nrm = NpzArray { shape: vec![p, 3], data: NpzData::F64(flat(&s.normals)) };
                let nbr = NpzArray { shape: vec![p, s.k], data: NpzData::I32(s.neighbors.clone()) };
                let msk = NpzArray { shape: vec![p], data: NpzData::Bool((0..p).map(|i| i < s.n).collect()) };
                let q = NpzArray { shape: vec![s.n, m, 3], data: NpzData::F64(flat(&prep.queries.points)) };
                let tr = NpzArray { shape: vec![s.n, m], data: NpzData::F64(prep.queries.truth.clone()) };
                let ver = NpzArray { shape: vec![s.vertices.len(), 3], data: NpzData::F64(flat(&s.vertices)) };
                let fac = NpzArray { shape: vec![s.faces.len(), 3], data: NpzData::I32(s.faces.iter().flat_map(|t| t.map(|j| j as i32)).collect()) };
                write_npz(&dir.join(format!("shape{si}.npz")), &[("positions", &pos), ("normals", &nrm), ("neighbors", &nbr), ("mask", &msk), ("queries", &q), ("truth", &tr), ("vertices", &ver), ("faces", &fac)]).unwrap();
            }
        }
    }
    json.truncate(json.len() - 2);
    json += "\n]\n";
    std::fs::write(dir.join("ours.json"), json).unwrap();
}

/// Step time and peak memory on the CPU: a blob shape at batch 1024 (phases 0–4) and the
/// pretraining blob at batch 32.
#[test]
#[ignore]
fn step_timing() {
    use ddg_points_as_tori::train::data::{DataOptions, load_file};
    let ctx = ctx32();
    let model = ddg_points_as_tori::Model::reference().unwrap();
    let params = NetParams::from_model(&model, &ctx);
    let threads: usize = std::env::var("PAT_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or(64);
    let chunk: usize = std::env::var("PAT_CHUNK").ok().and_then(|v| v.parse().ok()).unwrap_or(4);
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    let o = StepOptions { network_chunk: chunk, ..Default::default() };
    for (file, batch) in [("pretraining_blob.npz", 32usize), ("pretraining_blob.npz", 1024)] {
        let (shapes, _) = load_file(&std::path::Path::new(DATA).join(file), &DataOptions::default()).unwrap();
        pool.install(|| {
            let t = std::time::Instant::now();
            let prep = prepare_shape(&shapes[0], &o, (1, 0, 0)).unwrap();
            let tp = t.elapsed().as_secs_f64();
            for rep in 0..3 {
                let t = std::time::Instant::now();
                let out = step(&params, &ctx, &prep, 0, batch, &o, true).unwrap();
                eprintln!(
                    "{file} batch {batch} threads {threads} chunk {chunk} rep {rep}: step {:.2} s (plan {:.2}, forward {:.2}, loss {:.2}, backward {:.2}); prepare {tp:.2} s; rows fwd/bwd {}/{}; pairs {}; loss {:.4}; peak RSS {:.0} MB",
                    t.elapsed().as_secs_f64(), out.seconds[0], out.seconds[1], out.seconds[2], out.seconds[3], out.rows_forward, out.rows_backward, out.pairs, out.value.loss, peak_rss_mb()
                );
            }
        });
    }
}

/// A short CPU run of the pretraining phase on the blob, from scratch (logged to
/// `runs_tmp/pat-wp7/short/log.jsonl`), with a resume check.
#[test]
#[ignore]
fn short_run() {
    use ddg_points_as_tori::train::{Curriculum, TrainOptions, train};
    let steps: u64 = std::env::var("PAT_STEPS").ok().and_then(|v| v.parse().ok()).unwrap_or(200);
    let dir = out_dir("short");
    for f in ["checkpoint.safetensors", "checkpoint-adam.safetensors", "log.jsonl"] {
        let _ = std::fs::remove_file(dir.join(f));
    }
    let mut c = Curriculum::reference().unwrap();
    c.phases.truncate(1);
    let mut o = TrainOptions::reference(dir.clone()).unwrap();
    o.curriculum = c;
    o.max_steps = Some(steps);
    o.checkpoint_every = 64;
    o.validate = false;
    let t = std::time::Instant::now();
    let mut losses = Vec::new();
    let (_, rep) = train(std::path::Path::new(DATA), None, &o, &mut |l| {
        losses.push(l.loss);
        if l.step % 10 == 0 {
            eprintln!("step {} epoch {} lr {:.1e} loss {:.5} (distance {:.5}, eikonal {:.5}) {:.2} s", l.step, l.epoch, l.lr, l.loss, l.distance, l.eikonal, l.seconds.iter().sum::<f64>());
        }
    })
    .unwrap();
    eprintln!("{rep:?}; {:.1} s, peak RSS {:.0} MB", t.elapsed().as_secs_f64(), peak_rss_mb());
    let first: f64 = losses[..10].iter().sum::<f64>() / 10.0;
    let last: f64 = losses[losses.len() - 10..].iter().sum::<f64>() / 10.0;
    // the per-step loss is over different rows and queries each step (noisy); `short_curve` tracks a
    // fixed evaluation set
    eprintln!("mean loss of the first 10 steps {first:.5}, of the last 10 {last:.5}");
    assert_eq!(rep.steps, steps);
    assert!(losses.iter().all(|l| l.is_finite()));
}

/// Resumption: 2 + 4 steps with a restart in between equal 6 steps straight (bitwise).
#[test]
#[ignore]
fn resume_is_exact() {
    use ddg_points_as_tori::train::{Curriculum, TrainOptions, train};
    let mut c = Curriculum::reference().unwrap();
    c.phases.truncate(1);
    c.phases[0].batch_rows = 1024; // 2 steps per epoch of the one-shape phase
    let run = |name: &str, plan: &[u64]| {
        let dir = out_dir(name);
        for f in ["checkpoint.safetensors", "checkpoint-adam.safetensors", "log.jsonl"] {
            let _ = std::fs::remove_file(dir.join(f));
        }
        let mut model = None;
        for &m in plan {
            let mut o = TrainOptions::reference(dir.clone()).unwrap();
            o.curriculum = c.clone();
            o.max_steps = Some(m);
            o.checkpoint_every = 1;
            o.validate = false;
            model = Some(train(std::path::Path::new(DATA), None, &o, &mut |_| {}).unwrap().0);
        }
        model.unwrap()
    };
    // stop at a shape boundary (the end of epoch 1), resume at epoch 2
    let a = run("resume-a", &[6]);
    let b = run("resume-b", &[2, 6]);
    assert_eq!(a.weights(), b.weights());
}

/// The loss alone on one thread (profiling aid): 2048 queries of the pretraining blob.
#[test]
#[ignore]
fn loss_profile() {
    use ddg_points_as_tori::train::data::{DataOptions, load_file};
    let (shapes, _) = load_file(&std::path::Path::new(DATA).join("pretraining_blob.npz"), &DataOptions::default()).unwrap();
    let s = &shapes[0];
    let o = StepOptions::default();
    let prep = prepare_shape(s, &o, (1, 0, 0)).unwrap();
    let q = &prep.queries.points[..2048];
    let t = &prep.queries.truth[..2048];
    let plan = build_plan(q, t, &s.positions, &o.loss, None).unwrap();
    let mut r = ChaCha8Rng::seed_from_u64(1);
    let a: Vec<f64> = (0..s.n).flat_map(|_| [0.001 * u(&mut r), 0.1 * u(&mut r), 0.1 * u(&mut r), u(&mut r), -20.0 - u(&mut r), -20.0 - u(&mut r)]).collect();
    let cloud = CloudFit { p: &s.positions, n: &s.normals };
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    pool.install(|| {
        for want in [false, true, false, true] {
            let tm = std::time::Instant::now();
            let (v, _) = loss_and_grad(&ctx32(), &cloud, &a, &plan, &o.loss, want).unwrap();
            eprintln!("grad {want}: {:.3} s, {} pairs, loss {:.4}", tm.elapsed().as_secs_f64(), plan.pairs(), v.loss);
        }
    });
}

/// Padded clouds (D19): `n < P` rows and neighborhoods shorter than `k`.
#[test]
fn padded_clouds() {
    let ctx = ctx32();
    let params = small_params(&ctx, 7);
    for n in [10usize, 60] {
        let full = synth_shape(n, 16.min(n), 23);
        // pad to P = 80 rows and k = 16 columns
        let (p, k) = (80, 16);
        let mut s = full.clone();
        s.rows = p;
        s.k = k;
        s.positions.resize(p, [0.0; 3]);
        s.normals.resize(p, [0.0; 3]);
        s.neighbors = vec![-1; p * k];
        for i in 0..n {
            for c in 0..full.k {
                s.neighbors[i * k + c] = full.neighbors[i * full.k + c];
            }
        }
        let o = StepOptions { queries: QueryOptions { per_kind: 4, ..Default::default() }, ..Default::default() };
        let prep = prepare_shape(&s, &o, (1, 2, 3)).unwrap();
        assert_eq!(prep.net_rows.len(), n);
        assert_eq!(prep.feats.pad.is_some(), n < k);
        let out = step(&params, &ctx, &prep, 0, 64, &o, true).unwrap();
        assert!(out.value.loss.is_finite() && out.nonfinite_grads == 0, "n {n}: {:?} {}", out.value, out.nonfinite_grads);
        // the pad mask makes the padded tokens irrelevant: garbage in them changes nothing
        if n < k {
            let mut prep2 = prep.clone();
            for r in 0..prep2.feats.b {
                for t in n..k {
                    prep2.feats.feats[(r * k + t) * 6] = 7.0;
                }
            }
            let out2 = step(&params, &ctx, &prep2, 0, 64, &o, false).unwrap();
            assert_eq!(out2.value.loss, out.value.loss);
        }
        // the reference's padding: all P rows through the network, −1 → row P − 1
        let or = StepOptions { reference_padding: true, ..o };
        let prep_r = prepare_shape(&s, &or, (1, 2, 3)).unwrap();
        assert_eq!(prep_r.net_rows.len(), p);
        let out_r = step(&params, &ctx, &prep_r, 0, 64, &or, true).unwrap();
        assert!(out_r.value.loss.is_finite());
        eprintln!("n = {n}: loss {:.5} (ours) vs {:.5} (reference padding), non-finite gradient entries zeroed {}", out.value.loss, out_r.value.loss, out_r.nonfinite_grads);
        // empty slices have no queries
        assert_eq!(step(&params, &ctx, &prep, 64, 80, &o, true).unwrap().value.queries, 0);
    }
}

/// The pretraining setting (blob, batch 32, Adam, from scratch) through the step API, with the loss
/// on a fixed evaluation set (all 2048 rows, fixed queries) every 16 steps; `PAT_LR`, `PAT_STEPS`,
/// `PAT_THREADS`. Writes `runs_tmp/pat-wp7/curve-<lr>.tsv`.
#[test]
#[ignore]
fn short_curve() {
    use ddg_points_as_tori::train::data::{DataOptions, load_file};
    use ddg_points_as_tori::train::step::evaluate_shape;
    let lr: f64 = std::env::var("PAT_LR").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-3);
    let steps: usize = std::env::var("PAT_STEPS").ok().and_then(|v| v.parse().ok()).unwrap_or(128);
    let threads: usize = std::env::var("PAT_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or(64);
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    let ctx = ctx32();
    let (shapes, _) = load_file(&std::path::Path::new(DATA).join("pretraining_blob.npz"), &DataOptions::default()).unwrap();
    let s = &shapes[0];
    let o = StepOptions::default();
    let mut params = NetParams::init(&ModelConfig::default(), 42, &ctx);
    let mut adam = AdamTrainer::new(&params, AdamOptions::default());
    let mut out = String::from("step\ttrain_loss\teval_loss\teval_distance\teval_eikonal\tseconds\n");
    pool.install(|| {
        let eval_prep = prepare_shape(s, &o, (7, 7, 7)).unwrap();
        let eval = |p: &NetParams| {
            let v = evaluate_shape(p, &ctx, &eval_prep, 1024, &o).unwrap();
            let n = v.len() as f64;
            (v.iter().map(|x| x.loss).sum::<f64>() / n, v.iter().map(|x| x.distance).sum::<f64>() / n, v.iter().map(|x| x.eikonal).sum::<f64>() / n)
        };
        let t0 = std::time::Instant::now();
        let e = eval(&params);
        out += &format!("0\t\t{:.6}\t{:.6}\t{:.6}\t0\n", e.0, e.1, e.2);
        eprintln!("lr {lr:e} step 0: eval loss {:.5} (distance {:.5}, eikonal {:.5})", e.0, e.1, e.2);
        let mut prep = None;
        for k in 0..steps {
            let (epoch, b) = (k / 64, k % 64);
            if b == 0 {
                prep = Some(prepare_shape(s, &o, (42, epoch as u64, 0)).unwrap());
            }
            let r = step(&params, &ctx, prep.as_ref().unwrap(), 32 * b, 32 * b + 32, &o, true).unwrap();
            params = adam.update(&params, r.grads.as_ref().unwrap(), lr).unwrap();
            if (k + 1) % 16 == 0 {
                let e = eval(&params);
                out += &format!("{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.1}\n", k + 1, r.value.loss, e.0, e.1, e.2, t0.elapsed().as_secs_f64());
                eprintln!("lr {lr:e} step {}: train {:.5} eval loss {:.5} (distance {:.5}, eikonal {:.5}) {:.0} s", k + 1, r.value.loss, e.0, e.1, e.2, t0.elapsed().as_secs_f64());
            }
        }
    });
    std::fs::write(out_dir("curves").join(format!("curve-{lr:e}.tsv")), out).unwrap();
}
