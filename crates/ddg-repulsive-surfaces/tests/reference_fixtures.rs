//! Parity with the authors' program (`SPEC.md` §14.7, WP18): replays the fixtures of
//! `tests/fixtures/reference/`, dumped by the black-box build of the reference
//! (`/var/data/discrete-geometry/RepulsiveSurfaces-build`, `rsurfaces_headless --dump-initial
//! --dump-metric`, one thread), against our energy, differential and dense metric.
//!
//! Per mesh `NAME` (`NAME.obj` is the reference's input, read by both sides):
//! * `NAME.txt`: exponents, `θ`, sizes, the all-pairs energy and the Barnes–Hut energy at `θ = 0.5`;
//! * `NAME.dE_all.txt`, `NAME.dE_bh.txt`: the differentials (`|V| × 3`);
//! * `NAME.metric.txt`: the constrained dense `Hs` matrix of the reference
//!   (`GetHsMatrixConstrained`, Barnes–Hut energy, constraints barycenter + area + volume), which is
//!   `A ⊗ I₃` in its top-left block; we keep the `x`-channel block `A[i][j] = M[3i][3j]` for the
//!   first `rows` vertices and the 5 constraint rows in full (`3|V|` numbers each). The
//!   generation checked that the full matrix is exactly block-triplicated (defect 0) and symmetric
//!   (1e-13 absolute).
//!
//! The tests print the measured agreement (`cargo test --release --test reference_fixtures --
//! --nocapture`).

#![allow(clippy::needless_range_loop)]

use std::path::PathBuf;

use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::metric::SobolevWeights;
use ddg_repulsive_surfaces::metric::dense::dense_metric;
use ddg_repulsive_surfaces::{BhOptions, EnergyEval, Exponents, Surface, energy, read_obj};

struct Fixture {
    surface: Surface,
    positions: Vec<[f64; 3]>,
    exponents: Exponents,
    theta: f64,
    all_pairs: f64,
    barnes_hut: f64,
    d_all: Vec<[f64; 3]>,
    d_bh: Vec<[f64; 3]>,
}

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/reference")
}

fn read_text(file: &str) -> String {
    std::fs::read_to_string(dir().join(file)).unwrap_or_else(|e| panic!("{file}: {e}"))
}

fn numbers(line: &str) -> Vec<f64> {
    line.split_whitespace().map(|t| t.parse().expect("number")).collect()
}

fn triples(file: &str) -> Vec<[f64; 3]> {
    read_text(file)
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let v = numbers(l);
            [v[0], v[1], v[2]]
        })
        .collect()
}

fn load(name: &'static str) -> Fixture {
    let (surface, positions) = read_obj(dir().join(format!("{name}.obj"))).expect("obj");
    let mut kv = std::collections::HashMap::new();
    for l in read_text(&format!("{name}.txt")).lines() {
        let mut it = l.split_whitespace();
        kv.insert(it.next().unwrap().to_string(), it.next().unwrap().parse::<f64>().unwrap());
    }
    assert_eq!(kv["nV"] as usize, surface.n_vertices());
    assert_eq!(kv["nF"] as usize, surface.n_faces());
    Fixture {
        surface,
        positions,
        exponents: Exponents::new(kv["alpha"], kv["beta"]).unwrap(),
        theta: kv["theta"],
        all_pairs: kv["tpe_all_pairs"],
        barnes_hut: kv["tpe_barnes_hut"],
        d_all: triples(&format!("{name}.dE_all.txt")),
        d_bh: triples(&format!("{name}.dE_bh.txt")),
    }
}

/// `max |a − b|` over all components relative to `max |b|`.
fn rel_max(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()));
    let diff = a
        .iter()
        .zip(b)
        .flat_map(|(p, q)| (0..3).map(move |k| (p[k] - q[k]).abs()))
        .fold(0.0_f64, f64::max);
    diff / scale
}

/// Relative `2`-norm error.
fn rel_norm(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let num: f64 = a.iter().zip(b).flat_map(|(p, q)| (0..3).map(move |k| (p[k] - q[k]).powi(2))).sum();
    let den: f64 = b.iter().flatten().map(|x| x * x).sum();
    (num / den).sqrt()
}

const NAMES: [&str; 3] = ["tube_regular", "torus192", "torus504"];

/// Exact all-pairs energy and differential against the reference's `TPEnergyAllPairs`.
#[test]
fn all_pairs_energy_and_differential_match() {
    for name in NAMES {
        let f = load(name);
        let r = energy(&f.surface, &f.positions, f.exponents, &EnergyEval::Exact, &[], Some(1)).unwrap();
        let e_rel = (r.energy - f.all_pairs).abs() / f.all_pairs;
        let d_max = rel_max(&r.differential, &f.d_all);
        let d_nrm = rel_norm(&r.differential, &f.d_all);
        println!(
            "{name}: all-pairs E ours {:.15e} ref {:.15e} rel {e_rel:.2e}; dE max-rel {d_max:.2e} 2-norm-rel {d_nrm:.2e}",
            r.energy, f.all_pairs
        );
        assert!(e_rel < 1e-12, "{name}: energy relative error {e_rel:e}");
        assert!(d_max < 1e-11, "{name}: differential max-relative error {d_max:e}");
    }
}

/// Barnes–Hut at `θ = 0.5`: the two implementations split their trees differently (ours: median
/// split along the longest axis, leaves of up to `leaf_size` faces; the reference: its own
/// clustering), so the far-field approximation errors differ; both must be within the
/// approximation error of the exact value, and our value must agree with the reference's to about
/// that error.
#[test]
fn barnes_hut_agrees_within_the_approximation_error() {
    for name in NAMES {
        let f = load(name);
        for leaf_size in [1usize, 8] {
            let bh = EnergyEval::BarnesHut(BhOptions { theta: f.theta, leaf_size });
            let r = energy(&f.surface, &f.positions, f.exponents, &bh, &[], Some(1)).unwrap();
            let vs_ref = (r.energy - f.barnes_hut).abs() / f.barnes_hut;
            let ours_vs_exact = (r.energy - f.all_pairs).abs() / f.all_pairs;
            let ref_vs_exact = (f.barnes_hut - f.all_pairs).abs() / f.all_pairs;
            let d_vs_ref = rel_norm(&r.differential, &f.d_bh);
            let d_ours_exact = rel_norm(&r.differential, &f.d_all);
            let d_ref_exact = rel_norm(&f.d_bh, &f.d_all);
            println!(
                "{name} leaf {leaf_size}: BH E ours {:.12e} ref {:.12e}; |ours-ref|/ref {vs_ref:.2e}; ours vs exact {ours_vs_exact:.2e}, \
                 ref vs exact {ref_vs_exact:.2e}; dE 2-norm: ours vs ref {d_vs_ref:.2e}, ours vs exact {d_ours_exact:.2e}, ref vs exact {d_ref_exact:.2e}",
                r.energy, f.barnes_hut
            );
            assert!(ours_vs_exact < 5e-3 && ref_vs_exact < 5e-3, "{name}: BH energy too far from exact");
            assert!(vs_ref < 1e-2, "{name} leaf {leaf_size}: BH energies differ by {vs_ref:e}");
            assert!(d_vs_ref < 5e-2, "{name} leaf {leaf_size}: BH differentials differ by {d_vs_ref:e}");
        }
    }
}

/// `reference_compat`'s metric: the dense `A` (Eq. 12 + 13, `w_TP = 1`, `B⁰` on) against the
/// `x`-channel block of the reference's constrained matrix, and the constraint rows
/// (barycenter, area, volume) against its last five rows. The reference's metric uses its
/// Barnes–Hut-independent dense assembly (`FillMatrixHigh/Low`), so the entries agree to round-off;
/// the volume row is the reference's `3×` (D9).
#[test]
fn dense_constrained_metric_matches() {
    for name in NAMES {
        let f = load(name);
        let n = f.surface.n_vertices();
        let text = read_text(&format!("{name}.metric.txt"));
        let mut lines = text.lines();
        let nv: usize = lines.next().unwrap().split_whitespace().nth(1).unwrap().parse().unwrap();
        let rows: usize = lines.next().unwrap().split_whitespace().nth(1).unwrap().parse().unwrap();
        let crows: usize = lines.next().unwrap().split_whitespace().nth(1).unwrap().parse().unwrap();
        assert_eq!((nv, crows), (n, 5));
        let block: Vec<Vec<f64>> = (0..rows).map(|_| numbers(lines.next().unwrap())).collect();
        let cons: Vec<Vec<f64>> = (0..crows).map(|_| numbers(lines.next().unwrap())).collect();

        let pool = make_pool(Some(1)).unwrap();
        let w = SobolevWeights::default();
        let m = dense_metric(&pool, &f.surface, &f.positions, &f.exponents, &w, None).unwrap();
        let mut worst_abs = 0.0_f64;
        let mut worst_rel = 0.0_f64;
        let scale = block.iter().flatten().fold(0.0_f64, |s, x| s.max(x.abs()));
        for i in 0..rows {
            for j in 0..n {
                let (a, b) = (m.a.get(i, j), block[i][j]);
                worst_abs = worst_abs.max((a - b).abs());
                worst_rel = worst_rel.max((a - b).abs() / b.abs().max(1e-6 * scale));
            }
        }
        println!("{name}: metric A rows 0..{rows}: max |ours-ref| {worst_abs:.2e} (largest entry {scale:.3e}), entry-wise rel {worst_rel:.2e}");
        assert!(worst_abs / scale < 1e-10, "{name}: metric entries differ by {worst_abs:e}");

        let set = ConstraintSet::new(
            &f.surface,
            &f.positions,
            &[
                Constraint::Barycenter { component: None, target: None },
                Constraint::Area { component: None, target: None },
                Constraint::Volume { component: None, target: None },
            ],
            BarycenterMode::Global,
        )
        .unwrap();
        assert_eq!(set.n_rows(), 5);
        let jac = set.jacobian(&f.positions, &[0, 1, 2, 3, 4]).unwrap();
        for (r, label, factor) in [
            (0, "barycenter x", 1.0),
            (1, "barycenter y", 1.0),
            (2, "barycenter z", 1.0),
            (3, "area", 1.0),
            (4, "volume", 3.0),
        ] {
            let mut worst = 0.0_f64;
            let mut sc = 0.0_f64;
            for c in 0..3 * n {
                let (a, b) = (factor * jac.get(r, c), cons[r][c]);
                worst = worst.max((a - b).abs());
                sc = sc.max(b.abs());
            }
            println!("{name}: constraint row {label}: max |ours*{factor}-ref| {worst:.2e} (largest entry {sc:.3e})");
            assert!(worst <= 1e-12 * sc.max(1.0), "{name}: {label} row differs by {worst:e}");
        }
    }
}
