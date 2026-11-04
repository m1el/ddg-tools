//! Normal orientation (SPEC.md §8.4, Fig. 27): flip a fraction of clean normals at random and
//! run the iteration. The sweep (`#[ignore]`) writes `runs_tmp/pat-wp9a/orient-sweep.tsv`.

mod common;

use ddg_points_as_tori::analytic::Instance;
use ddg_points_as_tori::apps::orient::{Epsilon, OrientOptions, orient_normals, orientation_metric, pca_normals_random_sign};
use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::{Model, PatError};
use std::io::Write;

fn flip_fraction(normals: &[[f64; 3]], frac: f64, seed: u64) -> (Vec<[f64; 3]>, usize) {
    let mut rng = Rng::new(seed);
    let mut out = normals.to_vec();
    let mut c = 0;
    for v in out.iter_mut() {
        if rng.f64() < frac {
            *v = [-v[0], -v[1], -v[2]];
            c += 1;
        }
    }
    (out, c)
}

fn wrong(a: &[[f64; 3]], truth: &[[f64; 3]]) -> usize {
    a.iter().zip(truth).filter(|(x, t)| x[0] * t[0] + x[1] * t[1] + x[2] * t[2] < 0.0).count()
}

/// Clean clouds: (name, positions, true normals).
type Named = (String, Vec<[f64; 3]>, Vec<[f64; 3]>);

fn clouds(golden: bool) -> Vec<Named> {
    let mut v = Vec::new();
    for seed in [3u64, 12, 26, 40, 52] {
        let inst = Instance::from_seed(seed).expect("instance");
        let (p, n) = inst.solid.sample_surface(2048, seed).expect("sample");
        v.push((inst.name(), p, n));
    }
    if golden && let Some(root) = common::golden_root() {
        for name in ["synth_torus", "tooth", "knot", "block", "femur", "snail"] {
            let g = common::Golden::open(&root, name);
            v.push((name.to_string(), g.cloud.positions().to_vec(), g.cloud.normals().to_vec()));
        }
    }
    v
}

#[test]
fn orientation_fixes_thirty_percent_flips() {
    let model = Model::reference().unwrap();
    let mut worst: f64 = 0.0;
    for (name, p, n) in clouds(false).into_iter().take(3) {
        let (bad, c) = flip_fraction(&n, 0.3, 7);
        let opts = OrientOptions { threads: Some(16), ..Default::default() };
        let (out, rep) = orient_normals(&p, &bad, &model, &opts).unwrap();
        let w = wrong(&out, &n);
        let frac = w as f64 / n.len() as f64;
        eprintln!("{name}: flipped {c} -> wrong {w} ({:.4}), flips {:?}, negative_final {}", frac, rep.flips_per_iteration, rep.negative_final);
        worst = worst.max(frac);
    }
    assert!(worst < 0.05, "residual flipped fraction {worst}");
}

#[test]
fn clean_input_is_a_fixed_point_up_to_a_few_points() {
    let model = Model::reference().unwrap();
    let (name, p, n) = clouds(false).remove(0);
    let opts = OrientOptions { threads: Some(16), max_iterations: 3, ..Default::default() };
    let (out, rep) = orient_normals(&p, &n, &model, &opts).unwrap();
    let w = wrong(&out, &n);
    eprintln!("{name}: clean input, {w} wrong after, flips {:?}", rep.flips_per_iteration);
    assert!(w * 100 < n.len(), "{w} of {} flipped", n.len());
}

#[test]
fn bad_input_is_an_error_not_a_panic() {
    let model = Model::reference().unwrap();
    let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let n = vec![[0.0, 0.0, 1.0]];
    assert!(matches!(orient_normals(&p, &n, &model, &OrientOptions::default()), Err(PatError::LengthMismatch { .. })));
    let n2 = vec![[0.0, 0.0, 1.0]; 2];
    let o = OrientOptions { k: 0, ..Default::default() };
    assert!(matches!(orient_normals(&p, &n2, &model, &o), Err(PatError::InvalidOption { .. })));
    let o = OrientOptions { epsilon: Epsilon::Constant(f64::NAN), threads: Some(2), ..Default::default() };
    assert!(orient_normals(&p, &n2, &model, &o).is_err());
    assert!(pca_normals_random_sign(&p, 2, 0, Some(1)).is_err());
}

#[test]
fn metric_is_a_norm_for_nonzero_weights() {
    use ddg_points_as_tori::FitDiagnostics;
    use ddg_points_as_tori::fit::ClampFlags;
    let d = FitDiagnostics {
        kappa_min: 0.5,
        kappa_max: -2.0,
        v_min: [1.0, 0.0, 0.0],
        v_max: [0.0, 1.0, 0.0],
        n_star: [0.0, 0.0, 1.0],
        clamps: ClampFlags::default(),
    };
    let nrm = [0.0, 0.0, 1.0];
    let mut s = 5u64;
    for _ in 0..200 {
        let r = [common::lcg(&mut s) - 0.5, common::lcg(&mut s) - 0.5, common::lcg(&mut s) - 0.5];
        let q = [common::lcg(&mut s) - 0.5, common::lcg(&mut s) - 0.5, common::lcg(&mut s) - 0.5];
        let sum = [r[0] + q[0], r[1] + q[1], r[2] + q[2]];
        let m = |x| orientation_metric(x, &d, nrm, 3.0);
        assert!(m(sum) <= m(r) + m(q) + 1e-12);
        assert!((m([2.0 * r[0], 2.0 * r[1], 2.0 * r[2]]) - 2.0 * m(r)).abs() < 1e-12);
        assert!(m(r) > 0.0);
    }
}

/// The experiment of NOTES-WP9a: sweep eps and the flip fraction over many clouds, in parallel.
/// `cargo test --release --test orient -- --ignored --nocapture sweep`
#[test]
#[ignore]
fn sweep() {
    let model = Model::reference().unwrap();
    let cl = clouds(true);
    let eps: Vec<(&str, Epsilon)> = vec![
        ("0", Epsilon::Constant(0.0)),
        ("0.25/s", Epsilon::InverseSigma { scale: 0.25 }),
        ("0.5/s", Epsilon::InverseSigma { scale: 0.5 }),
        ("1/s", Epsilon::InverseSigma { scale: 1.0 }),
        ("2/s", Epsilon::InverseSigma { scale: 2.0 }),
        ("4/s", Epsilon::InverseSigma { scale: 4.0 }),
        ("16/s", Epsilon::InverseSigma { scale: 16.0 }),
    ];
    // jobs: (cloud, eps, frac)
    let mut jobs = Vec::new();
    for ci in 0..cl.len() {
        for (ei, _) in eps.iter().enumerate() {
            for frac in [0.3, 0.5] {
                jobs.push((ci, ei, frac));
            }
        }
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let out = std::sync::Mutex::new(Vec::new());
    std::thread::scope(|s| {
        for _ in 0..4 {
            s.spawn(|| loop {
                let j = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                if j >= jobs.len() {
                    break;
                }
                let (ci, ei, frac) = jobs[j];
                let (name, p, n) = &cl[ci];
                let (bad, c) = flip_fraction(n, frac, 11);
                let opts = OrientOptions { epsilon: eps[ei].1, threads: Some(4), ..Default::default() };
                let r = orient_normals(p, &bad, &model, &opts);
                let line = match r {
                    Ok((o, rep)) => format!(
                        "{name}\t{}\t{}\t{}\t{}\t{:.5}\t{}\t{:?}",
                        n.len(), eps[ei].0, frac, c, wrong(&o, n) as f64 / n.len() as f64, rep.negative_final, rep.flips_per_iteration
                    ),
                    Err(e) => format!("{name}\t{}\t{}\t{}\tERROR {e}", n.len(), eps[ei].0, frac),
                };
                eprintln!("{line}");
                out.lock().unwrap().push((j, line));
            });
        }
    });
    let mut v = out.into_inner().unwrap();
    v.sort();
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../runs_tmp/pat-wp9a/orient-sweep.tsv");
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, "cloud\tn\teps\tflipped_frac\tflipped\twrong_frac_after\tneg_final\tflips_per_iter").unwrap();
    for (_, l) in v {
        writeln!(f, "{l}").unwrap();
    }
}

/// PCA normals with random signs as input (the "unoriented estimate" question).
#[test]
#[ignore]
fn pca_random_signs() {
    let model = Model::reference().unwrap();
    for (name, p, n) in clouds(true) {
        let est = pca_normals_random_sign(&p, 16, 5, Some(4)).unwrap();
        let line_err = est.iter().zip(&n).filter(|(a, b)| (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]).abs() < 0.9).count();
        let w0 = wrong(&est, &n);
        let opts = OrientOptions { threads: Some(16), ..Default::default() };
        let (o, rep) = orient_normals(&p, &est, &model, &opts).unwrap();
        // up to a global sign
        let w = wrong(&o, &n).min(n.len() - wrong(&o, &n));
        eprintln!("{name}\tn={}\tline_off={line_err}\twrong_before={w0}\twrong_after(global sign free)={w}\tflips={:?}", n.len(), rep.flips_per_iteration);
    }
}
