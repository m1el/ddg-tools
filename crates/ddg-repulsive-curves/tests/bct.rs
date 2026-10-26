//! Tests of work package WP5 (`SPEC.md` §8.3, §13.4, D16, D25): the block cluster tree, the
//! hierarchical `Ā` product against the dense one, and a projected CG solve on top of it.
//!
//! The timing test is `#[ignore]`d; run it alone, in release:
//! `cargo test -p ddg-repulsive-curves --release --test bct -- --ignored --nocapture timing`.

#![allow(clippy::needless_range_loop)]

use std::sync::Arc;

use ddg_repulsive_curves::constraints::{ConstraintSet, default_constraints};
use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::metric::hierarchical::{
    CgOptions, HierarchicalCgSaddle, HierarchicalMetric,
};
use ddg_repulsive_curves::metric::{
    DenseMatrix, DofLayout, GramGeometry, GramQuadrature, kernel_matrices, sobolev_gram,
};
use ddg_repulsive_curves::solve::{SaddleSolver, dense_sobolev_solver, projected_gradient};
use ddg_repulsive_curves::spatial::bct::{BctOptions, BlockClusterTree};
use ddg_repulsive_curves::{Axis, CurveNetwork, Exponents, RcError, read_curve_obj};

type P = [f64; 3];

// ---------------------------------------------------------------- helpers

fn data(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/repulsive-curves")
        .join(rel)
}

fn have_data() -> bool {
    let ok = data("FreedmanHeWang/freedman-he-wang-600.obj").exists();
    if !ok {
        eprintln!("data/repulsive-curves not present; skipping");
    }
    ok
}

fn closed_edges(n: usize) -> Vec<[u32; 2]> {
    (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect()
}

fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [
                t.sin() + 2.0 * (2.0 * t).sin(),
                t.cos() - 2.0 * (2.0 * t).cos(),
                -(3.0 * t).sin(),
            ]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n)).unwrap(), pos)
}

fn circle(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [t.cos(), t.sin(), 0.0]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n)).unwrap(), pos)
}

fn planar_flower(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            let r = 1.0 + 0.3 * (3.0 * t).cos();
            [r * t.cos(), r * t.sin(), 0.5]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n)).unwrap(), pos)
}

/// A 3-junction graph: three arms from a common centre plus a loop through their ends.
fn tripod(m: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos = vec![[0.0, 0.0, 0.0]];
    let mut edges = Vec::new();
    for arm in 0..3 {
        let phi = std::f64::consts::TAU * arm as f64 / 3.0;
        let mut prev = 0u32;
        for k in 1..=m {
            let r = k as f64 / m as f64;
            pos.push([
                r * phi.cos(),
                r * phi.sin(),
                0.2 * r * r * (arm as f64 - 1.0),
            ]);
            let v = pos.len() as u32 - 1;
            edges.push([prev, v]);
            prev = v;
        }
    }
    (CurveNetwork::new(pos.len(), edges).unwrap(), pos)
}

/// Deterministic pseudo-random numbers in `[-0.5, 0.5)`.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    }
}

fn rand_vec(n: usize, seed: u64) -> Vec<f64> {
    let mut r = Lcg(seed);
    (0..n).map(|_| r.next()).collect()
}

/// Linear midpoint subdivision of every edge.
fn subdivide(net: &CurveNetwork, pos: &[P]) -> (CurveNetwork, Vec<P>) {
    let mut p = pos.to_vec();
    let mut edges = Vec::new();
    for &[a, b] in net.edges() {
        let m = p.len() as u32;
        let (pa, pb) = (pos[a as usize], pos[b as usize]);
        p.push([
            0.5 * (pa[0] + pb[0]),
            0.5 * (pa[1] + pb[1]),
            0.5 * (pa[2] + pb[2]),
        ]);
        edges.push([a, m]);
        edges.push([m, b]);
    }
    (CurveNetwork::new(p.len(), edges).unwrap(), p)
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn dotv(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn rel_err(a: &[f64], b: &[f64]) -> f64 {
    let d: Vec<f64> = a.iter().zip(b).map(|(x, y)| x - y).collect();
    norm(&d) / norm(b)
}

/// `(A ⊗ I_d) x` in the interleaved layout, rows in parallel on `pool`.
fn dense_apply_par(p: &rayon::ThreadPool, a: &DenseMatrix, d: usize, x: &[f64]) -> Vec<f64> {
    use rayon::prelude::*;
    let n = a.n();
    let rows: Vec<Vec<f64>> = p.install(|| {
        (0..n)
            .into_par_iter()
            .map(|u| {
                let row = a.row(u);
                (0..d)
                    .map(|s| (0..n).map(|v| row[v] * x[v * d + s]).sum())
                    .collect()
            })
            .collect()
    });
    rows.concat()
}

/// `(A ⊗ I_d) x` in the interleaved layout.
fn dense_apply(a: &DenseMatrix, d: usize, x: &[f64]) -> Vec<f64> {
    let n = a.n();
    let mut out = vec![0.0; n * d];
    for u in 0..n {
        let row = a.row(u);
        for s in 0..d {
            out[u * d + s] = (0..n).map(|v| row[v] * x[v * d + s]).sum();
        }
    }
    out
}

fn pool(n: usize) -> Arc<rayon::ThreadPool> {
    Arc::new(make_pool(Some(n)).unwrap())
}

fn opts(theta: f64) -> BctOptions {
    BctOptions {
        theta_x: theta,
        theta_t: theta,
        ..BctOptions::default()
    }
}

fn build(
    p: &Arc<rayon::ThreadPool>,
    net: &CurveNetwork,
    pos: &[P],
    ex: Exponents,
    q: GramQuadrature,
    layout: DofLayout,
    o: BctOptions,
) -> HierarchicalMetric {
    HierarchicalMetric::build(p.clone(), net, pos, ex, q, layout, o).unwrap()
}

fn flat(pos: &[P]) -> Vec<f64> {
    pos.iter().flatten().copied().collect()
}

// ---------------------------------------------------------------- options and errors

#[test]
fn bad_options_and_inputs_are_errors_not_panics() {
    let (net, pos) = trefoil(40);
    let p = pool(2);
    let ex = Exponents::default();
    let l = DofLayout::new(None);
    let q = GramQuadrature::Midpoint;
    for o in [
        BctOptions {
            theta_x: -0.1,
            ..BctOptions::default()
        },
        BctOptions {
            theta_t: f64::NAN,
            ..BctOptions::default()
        },
        BctOptions {
            theta_x: f64::INFINITY,
            ..BctOptions::default()
        },
        BctOptions {
            leaf_size: 0,
            ..BctOptions::default()
        },
    ] {
        let r = HierarchicalMetric::build(p.clone(), &net, &pos, ex, q, l, o);
        assert!(matches!(r, Err(RcError::InvalidOption { .. })), "{o:?}");
    }
    let r = HierarchicalMetric::build(p.clone(), &net, &pos[..39], ex, q, l, opts(0.5));
    assert!(matches!(r, Err(RcError::PositionCountMismatch { .. })));
    let mut bad = pos.clone();
    bad[3][1] = f64::NAN;
    let r = HierarchicalMetric::build(p.clone(), &net, &bad, ex, q, l, opts(0.5));
    assert!(matches!(r, Err(RcError::NonFinitePosition { vertex: 3 })));
    let mut bad = pos.clone();
    bad[5] = bad[4];
    let r = HierarchicalMetric::build(p.clone(), &net, &bad, ex, q, l, opts(0.5));
    assert!(matches!(r, Err(RcError::ZeroLengthEdge { edge: 4, .. })));

    let h = build(&p, &net, &pos, ex, q, l, opts(0.5));
    assert!(matches!(
        h.apply(&[0.0; 5]),
        Err(RcError::InvalidOption { .. })
    ));
    assert!(matches!(
        h.tree().multiply(&p, &[0usize], &[[1.0]; 3]),
        Err(RcError::InvalidOption { .. })
    ));
    assert!(matches!(
        h.tree()
            .multiply(&p, &[2usize], &vec![[1.0]; net.n_edges()]),
        Err(RcError::InvalidOption { .. })
    ));
    // A Jacobian of the wrong width, a bad CG tolerance, wrong right-hand sides.
    let set = ConstraintSet::new(&net, &pos, &default_constraints(), None).unwrap();
    let c = set.jacobian(&net, &pos).unwrap();
    let h2 = build(
        &p,
        &net,
        &pos,
        ex,
        q,
        DofLayout::new(Some(Axis::Z)),
        opts(0.5),
    );
    assert!(matches!(
        HierarchicalCgSaddle::new(h2, &c, CgOptions::default()),
        Err(RcError::InvalidOption { .. })
    ));
    let bad_tol = CgOptions {
        tolerance: 0.0,
        ..CgOptions::default()
    };
    assert!(matches!(
        HierarchicalCgSaddle::new(h.clone(), &c, bad_tol),
        Err(RcError::InvalidOption { .. })
    ));
    let mut s = HierarchicalCgSaddle::new(h, &c, CgOptions::default()).unwrap();
    assert!(matches!(
        s.solve(&[0.0; 3], &[]),
        Err(RcError::InvalidOption { .. })
    ));
}

// ---------------------------------------------------------------- block structure

/// Every ordered pair `(I, J)` lies in exactly one block; admissible blocks contain neither
/// `I = J` nor adjacent pairs; the partition is symmetric.
fn check_partition(net: &CurveNetwork, h: &HierarchicalMetric) -> (usize, u64) {
    let m = net.n_edges();
    let t = h.tree();
    let mut cover = vec![0u8; m * m];
    let mut adm = std::collections::BTreeSet::new();
    let mut n_adm = 0;
    for (a, b, k) in t.admissible_blocks() {
        n_adm += 1;
        adm.insert((a, b));
        assert!(a != b);
        assert!(k[0] > 0.0 && k[1] >= 0.0, "kappa {k:?}");
        for &i in t.elements_of(a) {
            for &j in t.elements_of(b) {
                let (i, j) = (i as usize, j as usize);
                assert!(
                    !net.edges_adjacent(i, j),
                    "admissible block ({a}, {b}) contains the adjacent pair ({i}, {j})"
                );
                cover[i * m + j] += 1;
            }
        }
    }
    for (a, b) in t.inadmissible_blocks() {
        for &i in t.elements_of(a) {
            for &j in t.elements_of(b) {
                cover[i as usize * m + j as usize] += 1;
            }
        }
    }
    assert!(cover.iter().all(|&c| c == 1), "blocks do not partition");
    for &(a, b) in &adm {
        assert!(
            adm.contains(&(b, a)),
            "partition not symmetric at ({a}, {b})"
        );
    }
    assert_eq!(n_adm, h.stats().admissible_blocks);
    (n_adm, h.stats().admissible_pairs)
}

#[test]
fn blocks_partition_the_matrix_symmetrically() {
    let p = pool(4);
    let ex = Exponents::default();
    let mut cases = vec![trefoil(200), tripod(40), planar_flower(150)];
    if have_data() {
        cases.push(read_curve_obj(data("FreedmanHeWang/freedman-he-wang-600.obj")).unwrap());
    }
    for (net, pos) in &cases {
        for theta in [0.0, 0.25, 0.5, 1.0, 4.0] {
            let h = build(
                &p,
                net,
                pos,
                ex,
                GramQuadrature::Midpoint,
                DofLayout::new(None),
                opts(theta),
            );
            let (n_adm, pairs) = check_partition(net, &h);
            if theta == 0.0 {
                assert_eq!(n_adm, 0);
            }
            eprintln!(
                "{} edges, theta {theta}: {n_adm} admissible blocks, {:.1} % of pairs, {} near entries",
                net.n_edges(),
                100.0 * pairs as f64 / (net.n_edges() as f64).powi(2),
                h.stats().near_entries
            );
        }
    }
}

/// With huge θ the spatial and tangent criteria accept nearly every pair of disjoint
/// clusters, so the endpoint-box guard (D16) is what keeps neighbouring edges apart: no
/// admissible block contains an adjacent pair, and every adjacent pair is in an exact block,
/// where its entry is zero.
#[test]
fn guard_keeps_adjacent_pairs_out_of_admissible_blocks() {
    let p = pool(4);
    let ex = Exponents::default();
    let mut cases = vec![circle(64), trefoil(300), tripod(50)];
    if have_data() {
        cases.push(read_curve_obj(data("FreedmanHeWang/freedman-he-wang-1200.obj")).unwrap());
    }
    for (net, pos) in &cases {
        let o = BctOptions {
            theta_x: 1e3,
            theta_t: 1e3,
            leaf_pairs: 0,
            leaf_size: 1,
        };
        let h = build(
            &p,
            net,
            pos,
            ex,
            GramQuadrature::Midpoint,
            DofLayout::new(None),
            o,
        );
        let (n_adm, _) = check_partition(net, &h);
        assert!(n_adm > 0);
        // Every adjacent pair lies in an inadmissible block.
        let m = net.n_edges();
        let t = h.tree();
        let mut exact = vec![false; m * m];
        for (a, b) in t.inadmissible_blocks() {
            for &i in t.elements_of(a) {
                for &j in t.elements_of(b) {
                    exact[i as usize * m + j as usize] = true;
                }
            }
        }
        for i in 0..m {
            for j in 0..m {
                if net.edges_adjacent(i, j) {
                    assert!(exact[i * m + j], "adjacent pair ({i}, {j}) not exact");
                }
            }
        }
        // The product still has Ā 𝟙 = 0 (translations in the kernel) and K entries 0 on
        // adjacent pairs: compare K̃ψ with ψ supported on the two neighbours of one edge.
        let ones = vec![1.0; 3 * net.n_vertices()];
        let a1 = h.apply(&ones).unwrap();
        let x = flat(pos);
        let scale = norm(&h.apply(&x).unwrap()) / norm(&x);
        assert!(
            norm(&a1) <= 1e-12 * scale * norm(&ones),
            "Ā𝟙 = {:e}",
            norm(&a1)
        );
        let i0 = 0usize;
        let mut psi = vec![[0.0f64; 2]; m];
        for j in 0..m {
            if net.edges_adjacent(i0, j) {
                psi[j] = [1.0, 1.0];
            }
        }
        let phi = t.multiply(&p, &[0, 1], &psi).unwrap();
        assert_eq!(phi[i0], [0.0, 0.0], "adjacent entries must be exactly zero");
    }
}

// ---------------------------------------------------------------- low-rank blocks

/// The hierarchical `K̃` assembled explicitly from the blocks: rank-1 entries
/// `ℓ_I κ(𝒜, ℬ) ℓ_J` on admissible blocks, the exact entries elsewhere.
fn assemble(net: &CurveNetwork, pos: &[P], ex: Exponents, h: &HierarchicalMetric) -> [Vec<f64>; 2] {
    let m = net.n_edges();
    let pl = make_pool(Some(4)).unwrap();
    let geom = GramGeometry::new(net, pos, ex);
    let k = kernel_matrices(&pl, net, &geom, h.quadrature());
    let mut out = [vec![0.0; m * m], vec![0.0; m * m]];
    let t = h.tree();
    for (a, b, kap) in t.admissible_blocks() {
        for &i in t.elements_of(a) {
            for &j in t.elements_of(b) {
                let (i, j) = (i as usize, j as usize);
                let ll = geom.length[i] * geom.length[j];
                out[0][i * m + j] = ll * kap[0];
                out[1][i * m + j] = ll * kap[1];
            }
        }
    }
    for (a, b) in t.inadmissible_blocks() {
        for &i in t.elements_of(a) {
            for &j in t.elements_of(b) {
                let (i, j) = (i as usize, j as usize);
                out[0][i * m + j] = k.kb.get(i, j);
                out[1][i * m + j] = k.k0.get(i, j);
            }
        }
    }
    out
}

/// The rank-1 blocks against the explicit blocks: `κ` equals the lumped kernel at the cluster
/// centroids and normalized mean tangents (computed here from the edges), the fast-multipole
/// product equals the explicitly assembled `K̃` to rounding, `K̃` is exactly symmetric, and the
/// entrywise error of the rank-1 blocks shrinks with θ.
#[test]
fn lowrank_blocks_match_the_explicit_blocks() {
    let p = pool(4);
    let (net, pos, ex) = if have_data() {
        let (n, x) = read_curve_obj(data("FreedmanHeWang/freedman-he-wang-600.obj")).unwrap();
        (n, x, Exponents::new(2.0, 4.5).unwrap())
    } else {
        let (n, x) = trefoil(400);
        (n, x, Exponents::default())
    };
    let m = net.n_edges();
    let geom = GramGeometry::new(&net, &pos, ex);
    let pl = make_pool(Some(4)).unwrap();
    let kd = kernel_matrices(&pl, &net, &geom, GramQuadrature::Midpoint);
    let sigma = ex.sigma();
    let mut last_err = f64::INFINITY;
    for theta in [1.0, 0.5, 0.25, 0.1] {
        let h = build(
            &p,
            &net,
            &pos,
            ex,
            GramQuadrature::Midpoint,
            DofLayout::new(None),
            opts(theta),
        );
        let t = h.tree();
        // κ against the lumped kernel computed from scratch.
        let centroid = |a: usize| {
            let (mut l, mut x, mut tt) = (0.0, [0.0; 3], [0.0; 3]);
            for &e in t.elements_of(a) {
                let e = e as usize;
                l += geom.length[e];
                for k in 0..3 {
                    x[k] += geom.length[e] * geom.midpoint[e][k];
                    tt[k] += geom.length[e] * geom.tangent[e][k];
                }
            }
            let tn = norm(&tt);
            (
                [x[0] / l, x[1] / l, x[2] / l],
                [tt[0] / tn, tt[1] / tn, tt[2] / tn],
            )
        };
        let mut worst_kappa: f64 = 0.0;
        let (mut worst_block, mut sum_block, mut nb) = (0.0f64, 0.0, 0usize);
        for (a, b, kap) in t.admissible_blocks() {
            let ((xa, ta), (xb, tb)) = (centroid(a), centroid(b));
            let d = [xa[0] - xb[0], xa[1] - xb[1], xa[2] - xb[2]];
            let r2 = dotv(&d, &d);
            let r = r2.sqrt();
            let kb = 2.0 * r.powf(-(2.0 * sigma + 1.0));
            let (pa, pb) = (dotv(&ta, &d), dotv(&tb, &d));
            let k0 = ((r2 - pa * pa) + (r2 - pb * pb)) * r.powf(-(2.0 * sigma + 5.0));
            worst_kappa = worst_kappa
                .max(((kap[0] - kb) / kb).abs())
                .max((kap[1] - k0).abs() / (kb / r2));
            // Entrywise error of the block (K^B, relative to the block's largest entry).
            let (mut emax, mut kmax) = (0.0f64, 0.0f64);
            for &i in t.elements_of(a) {
                for &j in t.elements_of(b) {
                    let (i, j) = (i as usize, j as usize);
                    let exact = kd.kb.get(i, j);
                    let low = geom.length[i] * geom.length[j] * kap[0];
                    emax = emax.max((exact - low).abs());
                    kmax = kmax.max(exact.abs());
                }
            }
            worst_block = worst_block.max(emax / kmax);
            sum_block += emax / kmax;
            nb += 1;
        }
        assert!(
            worst_kappa <= 1e-12,
            "theta {theta}: kappa off by {worst_kappa:e}"
        );
        let mean_block = sum_block / nb.max(1) as f64;
        eprintln!(
            "theta {theta}: {nb} admissible blocks, K^B block error max {worst_block:.3e} mean {mean_block:.3e}"
        );
        if nb > 0 {
            assert!(
                mean_block < last_err,
                "block error does not shrink with theta"
            );
            last_err = mean_block;
        }
        // Explicit K̃: exactly symmetric, and the FMM product equals it.
        let kt = assemble(&net, &pos, ex, &h);
        for k in 0..2 {
            for i in 0..m {
                for j in 0..i {
                    assert_eq!(
                        kt[k][i * m + j].to_bits(),
                        kt[k][j * m + i].to_bits(),
                        "K~ not symmetric at ({i}, {j})"
                    );
                }
            }
        }
        let psi_v = rand_vec(2 * m, 11);
        let psi: Vec<[f64; 2]> = (0..m).map(|i| [psi_v[2 * i], psi_v[2 * i + 1]]).collect();
        let phi = t.multiply(&p, &[0, 1], &psi).unwrap();
        for k in 0..2 {
            let want: Vec<f64> = (0..m)
                .map(|i| (0..m).map(|j| kt[k][i * m + j] * psi[j][k]).sum())
                .collect();
            let got: Vec<f64> = phi.iter().map(|v| v[k]).collect();
            let e = rel_err(&got, &want);
            assert!(
                e <= 1e-13,
                "theta {theta}, kernel {k}: FMM vs explicit {e:e}"
            );
            // Row sums are K~ 𝟙.
            let rs: Vec<f64> = t.row_sums().iter().map(|v| v[k]).collect();
            let want1: Vec<f64> = (0..m)
                .map(|i| kt[k][i * m..(i + 1) * m].iter().sum())
                .collect();
            assert!(rel_err(&rs, &want1) <= 1e-13);
        }
    }
}

// ---------------------------------------------------------------- against the dense product

/// Relative errors of `Ā_h x` for random `x` and for the smooth `x = γ`, and of `K̃^B ψ`,
/// `K̃⁰ ψ` for random `ψ`. `Āγ` is small against `‖Ā‖ ‖γ‖` (the operator nearly annihilates
/// smooth fields), so the second error is amplified by that cancellation, rounding included
/// (≈ 1e-11 at θ = 0 on Freedman 2400).
struct Err4 {
    rand: f64,
    smooth: f64,
    kb: f64,
    k0: f64,
}

#[allow(clippy::too_many_arguments)]
fn errors(
    net: &CurveNetwork,
    pos: &[P],
    ex: Exponents,
    q: GramQuadrature,
    layout: DofLayout,
    a_dense: &DenseMatrix,
    kd: &ddg_repulsive_curves::metric::KernelMatrices,
    o: BctOptions,
    p: &Arc<rayon::ThreadPool>,
) -> (Err4, HierarchicalMetric) {
    let h = build(p, net, pos, ex, q, layout, o);
    let d = layout.dim();
    let x = rand_vec(d * net.n_vertices(), 5);
    let rand = rel_err(&h.apply(&x).unwrap(), &dense_apply(a_dense, d, &x));
    let xs = layout.gather(pos);
    let smooth = rel_err(&h.apply(&xs).unwrap(), &dense_apply(a_dense, d, &xs));
    let m = net.n_edges();
    let psi_v = rand_vec(2 * m, 9);
    let psi: Vec<[f64; 2]> = (0..m).map(|i| [psi_v[2 * i], psi_v[2 * i + 1]]).collect();
    let phi = h.tree().multiply(p, &[0, 1], &psi).unwrap();
    let ek = |k: usize, km: &DenseMatrix| {
        let got: Vec<f64> = phi.iter().map(|v| v[k]).collect();
        let col: Vec<f64> = psi.iter().map(|v| v[k]).collect();
        rel_err(&got, &km.mul_vec(&col))
    };
    let e = Err4 {
        rand,
        smooth,
        kb: ek(0, &kd.kb),
        k0: ek(1, &kd.k0),
    };
    (e, h)
}

/// With every block inadmissible (θ = 0, or a huge `leaf_pairs`) the hierarchical product is
/// the dense one to rounding, for both quadratures, in 3-D and planar layouts; as θ → 0 the
/// error goes to 0.
#[test]
fn converges_to_the_dense_product_and_equals_it_when_all_blocks_are_exact() {
    let p = pool(8);
    let pl = make_pool(Some(8)).unwrap();
    let ex = Exponents::default();
    let cases = [
        (trefoil(300), DofLayout::new(None)),
        (tripod(60), DofLayout::new(None)),
        (planar_flower(200), DofLayout::new(Some(Axis::Z))),
    ];
    for ((net, pos), layout) in &cases {
        for q in [GramQuadrature::Midpoint, GramQuadrature::Trapezoid] {
            let g = sobolev_gram(&pl, net, pos, ex, q);
            let geom = GramGeometry::new(net, pos, ex);
            let kd = kernel_matrices(&pl, net, &geom, q);
            let mut prev = f64::INFINITY;
            for theta in [1.0, 0.5, 0.25, 0.1, 0.03, 0.0] {
                let (e, h) = errors(net, pos, ex, q, *layout, &g.a, &kd, opts(theta), &p);
                eprintln!(
                    "{} edges {q:?} theta {theta:<5}: rand {:.2e} smooth {:.2e} K^B {:.2e} K0 {:.2e} ({} adm blocks)",
                    net.n_edges(),
                    e.rand,
                    e.smooth,
                    e.kb,
                    e.k0,
                    h.stats().admissible_blocks
                );
                assert!(
                    e.rand <= prev.max(1e-12) * 1.5 + 1e-15,
                    "error grows at {theta}"
                );
                prev = e.rand;
                if theta == 0.0 {
                    for v in [e.rand, e.smooth, e.kb, e.k0] {
                        assert!(v <= 1e-12, "theta 0: {v:e}");
                    }
                }
            }
            // All blocks inadmissible through leaf_pairs as well.
            let o = BctOptions {
                leaf_pairs: usize::MAX,
                ..BctOptions::default()
            };
            let (e, h) = errors(net, pos, ex, q, *layout, &g.a, &kd, o, &p);
            assert_eq!(h.stats().admissible_blocks, 0);
            assert!(e.rand <= 1e-12 && e.smooth <= 1e-12 && e.kb <= 1e-12 && e.k0 <= 1e-12);
        }
    }
}

/// Error of the hierarchical product against ε = θ_x = θ_T on the Freedman–He–Wang curves
/// ((2, 4.5)) and a perturbed Knot128 knot ((3, 6)); the table is printed (SPEC §13.4) and
/// the error at the default (0.5) is checked.
#[test]
fn error_against_eps_on_datasets() {
    if !have_data() {
        return;
    }
    let files = [
        (
            "freedman-600",
            "FreedmanHeWang/freedman-he-wang-600.obj",
            Exponents::new(2.0, 4.5).unwrap(),
        ),
        (
            "freedman-1200",
            "FreedmanHeWang/freedman-he-wang-1200.obj",
            Exponents::new(2.0, 4.5).unwrap(),
        ),
        (
            "freedman-2400",
            "FreedmanHeWang/freedman-he-wang-2400.obj",
            Exponents::new(2.0, 4.5).unwrap(),
        ),
        (
            "knot128-perturbed-0001",
            "KnotDatasets/Knot128/perturbed/knot0001.obj",
            Exponents::default(),
        ),
    ];
    let p = pool(16);
    let pl = make_pool(Some(16)).unwrap();
    let q = GramQuadrature::Midpoint;
    let layout = DofLayout::new(None);
    eprintln!(
        "| curve | edges | eps | Ā err (random x) | Ā err (x = γ) | K^B err | K⁰ err | adm. pairs % | near entries |"
    );
    for (name, rel, ex) in files {
        let (net, pos) = read_curve_obj(data(rel)).unwrap();
        let g = sobolev_gram(&pl, &net, &pos, ex, q);
        let geom = GramGeometry::new(&net, &pos, ex);
        let kd = kernel_matrices(&pl, &net, &geom, q);
        let m = net.n_edges() as f64;
        let mut at_default = None;
        for eps in [1.0, 0.5, 1.0 / 3.0, 0.25, 0.1, 0.03, 0.0] {
            let (e, h) = errors(&net, &pos, ex, q, layout, &g.a, &kd, opts(eps), &p);
            let s = h.stats();
            eprintln!(
                "| {name} | {} | {eps:.3} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.1} | {} |",
                net.n_edges(),
                e.rand,
                e.smooth,
                e.kb,
                e.k0,
                100.0 * s.admissible_pairs as f64 / (m * m),
                s.near_entries
            );
            if eps == 0.5 {
                at_default = Some(e.rand);
            }
            if eps == 0.25 {
                assert!(e.rand <= 1e-2, "{name}: error {:e} at eps 0.25", e.rand);
            }
            if eps == 0.0 {
                assert!(e.rand <= 1e-12 && e.smooth <= 1e-9, "{name}: eps 0");
            }
        }
        // The SPEC's target is 1e-2 at the defaults (§13.4); the default 0.5 of §10 misses it
        // (3-4e-2 for random x; reported, to be decided in §14.3). Recorded and bounded here.
        let e = at_default.unwrap();
        assert!(e <= 0.1, "{name}: error {e:e} at the default eps");
    }
}

// ---------------------------------------------------------------- symmetry and determinism

/// `uᵀ Ā_h v = vᵀ Ā_h u` to rounding: `K̃` is exactly symmetric, so the defect is only the
/// floating-point evaluation (≤ 1e-13 relative to `‖u‖ ‖Ā_h v‖ + ‖v‖ ‖Ā_h u‖`).
#[test]
fn hierarchical_operator_is_symmetric() {
    let p = pool(8);
    let mut cases = vec![(trefoil(400), Exponents::default())];
    if have_data() {
        cases.push((
            read_curve_obj(data("FreedmanHeWang/freedman-he-wang-1200.obj")).unwrap(),
            Exponents::new(2.0, 4.5).unwrap(),
        ));
    }
    for ((net, pos), ex) in &cases {
        for theta in [0.25, 0.5, 1.0] {
            let h = build(
                &p,
                net,
                pos,
                *ex,
                GramQuadrature::Midpoint,
                DofLayout::new(None),
                opts(theta),
            );
            let n = h.n_dofs();
            for seed in 0..3 {
                let u = rand_vec(n, 100 + seed);
                let v = rand_vec(n, 200 + seed);
                let (au, av) = (h.apply(&u).unwrap(), h.apply(&v).unwrap());
                let (uav, vau) = (dotv(&u, &av), dotv(&v, &au));
                let scale = norm(&u) * norm(&av) + norm(&v) * norm(&au);
                let defect = (uav - vau).abs() / scale;
                assert!(defect <= 1e-13, "theta {theta}: defect {defect:e}");
                // Positive semidefinite on these samples.
                assert!(dotv(&u, &au) >= -1e-13 * norm(&u) * norm(&au));
            }
        }
    }
}

#[test]
fn bitwise_identical_at_1_7_and_64_threads() {
    let mut cases = vec![(trefoil(500), Exponents::default())];
    if have_data() {
        cases.push((
            read_curve_obj(data("FreedmanHeWang/freedman-he-wang-600.obj")).unwrap(),
            Exponents::new(2.0, 4.5).unwrap(),
        ));
    }
    for ((net, pos), ex) in &cases {
        for theta in [0.25, 0.5] {
            let x = rand_vec(3 * net.n_vertices(), 3);
            let mut reference: Option<(Vec<u64>, Vec<u64>)> = None;
            for threads in [1, 7, 64] {
                let h = build(
                    &pool(threads),
                    net,
                    pos,
                    *ex,
                    GramQuadrature::Midpoint,
                    DofLayout::new(None),
                    opts(theta),
                );
                let y: Vec<u64> = h.apply(&x).unwrap().iter().map(|v| v.to_bits()).collect();
                let rs: Vec<u64> = h
                    .tree()
                    .row_sums()
                    .iter()
                    .flatten()
                    .map(|v| v.to_bits())
                    .collect();
                match &reference {
                    None => reference = Some((y, rs)),
                    Some((y0, rs0)) => {
                        assert!(y == *y0, "product differs at {threads} threads");
                        assert!(rs == *rs0, "row sums differ at {threads} threads");
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------- solves

/// CG on the projected hierarchical system (Eq. 30 through WP2's projector) reproduces the
/// dense projected gradient (Eq. 22) on Freedman 600 with barycenter + edge lengths: to the CG
/// tolerance when every block is exact, and to the approximation error of the operator at the
/// default θ. Also a right side with `b ≠ 0` (Eq. 33) on a small curve.
#[test]
fn projected_cg_reproduces_the_dense_projected_gradient() {
    let (net, pos, ex) = if have_data() {
        let (n, x) = read_curve_obj(data("FreedmanHeWang/freedman-he-wang-600.obj")).unwrap();
        (n, x, Exponents::new(2.0, 4.5).unwrap())
    } else {
        let (n, x) = trefoil(300);
        (n, x, Exponents::default())
    };
    let p = pool(16);
    let q = GramQuadrature::Midpoint;
    let set = ConstraintSet::new(&net, &pos, &default_constraints(), None).unwrap();
    set.analyze(&net, &pos).unwrap();
    let layout = set.layout();
    let c = set.jacobian(&net, &pos).unwrap();
    let (_, diff) = energy_and_differential(&p, &net, &pos, ex);
    let mut dense = dense_sobolev_solver(p.clone(), &net, &pos, ex, q, &set).unwrap();
    let gd = projected_gradient(&mut dense, layout, &diff).unwrap();
    let xd = layout.gather(&gd.direction);
    for (theta, bound) in [(0.0, 1e-7), (0.25, 2e-2), (0.5, 5e-2)] {
        let h = build(&p, &net, &pos, ex, q, layout, opts(theta));
        let o = CgOptions {
            tolerance: 1e-12,
            max_iterations: 100_000,
        };
        let mut s = HierarchicalCgSaddle::new(h, &c, o).unwrap();
        let gh = projected_gradient(&mut s, layout, &diff).unwrap();
        let xh = layout.gather(&gh.direction);
        let e = rel_err(&xh, &xd);
        let mut cx = vec![0.0; c.nrows()];
        c.mul_vec(&xh, &mut cx).unwrap();
        eprintln!(
            "{} edges, theta {theta}: {} CG iterations, saddle residual {:.2e}, |g_h - g_d|/|g_d| = {e:.2e}, |C g_h| = {:.2e}, sobolev norm² {:.6e} vs {:.6e}",
            net.n_edges(),
            gh.linear.iterations,
            gh.linear.relative_residual,
            norm(&cx),
            gh.sobolev_norm_sq,
            gd.sobolev_norm_sq
        );
        assert!(gh.linear.relative_residual <= 1e-8, "{:?}", gh.linear);
        assert!(norm(&cx) <= 1e-10 * norm(&xh));
        assert!(e <= bound, "theta {theta}: {e:e} > {bound:e}");
    }

    // Eq. 33 with b ≠ 0 on a small curve, against the dense saddle solve.
    let (net, pos) = trefoil(60);
    let ex = Exponents::default();
    let set = ConstraintSet::new(&net, &pos, &default_constraints(), None).unwrap();
    let c = set.jacobian(&net, &pos).unwrap();
    let mut dense = dense_sobolev_solver(p.clone(), &net, &pos, ex, q, &set).unwrap();
    let h = build(&p, &net, &pos, ex, q, set.layout(), opts(0.0));
    let o = CgOptions {
        tolerance: 1e-13,
        max_iterations: 100_000,
    };
    let mut s = HierarchicalCgSaddle::new(h, &c, o).unwrap();
    let a = rand_vec(s.n_dofs(), 21);
    let b = rand_vec(s.n_constraints(), 22);
    let sd = dense.solve(&a, &b).unwrap();
    let sh = s.solve(&a, &b).unwrap();
    assert!(
        rel_err(&sh.x, &sd.x) <= 1e-8,
        "x: {:e}",
        rel_err(&sh.x, &sd.x)
    );
    assert!(
        rel_err(&sh.lambda, &sd.lambda) <= 1e-7,
        "lambda: {:e}",
        rel_err(&sh.lambda, &sd.lambda)
    );
}

// ---------------------------------------------------------------- generic use

/// The tree is generic: a scalar kernel `1/|p − q|` on random points in `R²` (no curve type).
#[test]
fn generic_tree_with_a_point_kernel() {
    use ddg_repulsive_curves::spatial::aabb::Aabb;
    use ddg_repulsive_curves::spatial::bct::ClusterKernel;
    use ddg_repulsive_curves::spatial::bvh::{Bvh, Element};

    struct Pt([f64; 2]);
    impl Element<2> for Pt {
        fn point(&self) -> [f64; 2] {
            self.0
        }
        fn mass(&self) -> f64 {
            1.0
        }
    }
    struct Coulomb<'a> {
        pts: &'a [Pt],
        boxes: Vec<(Aabb<2>, [f64; 2])>,
        theta: f64,
    }
    fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
    }
    impl ClusterKernel<1> for Coulomb<'_> {
        fn admissible(&self, a: usize, b: usize) -> bool {
            if self.theta == 0.0 {
                return false;
            }
            let (ba, ca) = &self.boxes[a];
            let (bb, cb) = &self.boxes[b];
            ba.half_diagonal().max(bb.half_diagonal()) <= self.theta * dist(*ca, *cb)
        }
        fn far(&self, a: usize, b: usize) -> [f64; 1] {
            [1.0 / dist(self.boxes[a].1, self.boxes[b].1)]
        }
        fn near(&self, i: usize, j: usize) -> Option<[f64; 1]> {
            Some([1.0 / dist(self.pts[i].0, self.pts[j].0)])
        }
        fn weight(&self, _: usize) -> f64 {
            1.0
        }
    }
    let n = 500;
    let r = rand_vec(2 * n, 77);
    let pts: Vec<Pt> = (0..n).map(|i| Pt([r[2 * i], r[2 * i + 1]])).collect();
    let bvh = Bvh::build(&pts, 4).unwrap();
    let boxes: Vec<(Aabb<2>, [f64; 2])> = bvh
        .nodes()
        .iter()
        .map(|nd| (nd.bbox, nd.centroid()))
        .collect();
    let pl = make_pool(Some(4)).unwrap();
    let psi: Vec<[f64; 1]> = rand_vec(n, 78).into_iter().map(|v| [v]).collect();
    let exact: Vec<f64> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i)
                .map(|j| psi[j][0] / dist(pts[i].0, pts[j].0))
                .sum()
        })
        .collect();
    let mut prev = f64::INFINITY;
    for theta in [0.5, 0.25, 0.1, 0.0] {
        let k = Coulomb {
            pts: &pts,
            boxes: boxes.clone(),
            theta,
        };
        let t = BlockClusterTree::build(&pl, &bvh, &k, 16).unwrap();
        let got: Vec<f64> = t
            .multiply(&pl, &[0], &psi)
            .unwrap()
            .iter()
            .map(|v| v[0])
            .collect();
        let e = rel_err(&got, &exact);
        eprintln!("coulomb theta {theta}: err {e:.2e}, {:?}", t.stats());
        assert!(e < prev || e <= 1e-14);
        prev = e;
    }
    assert!(prev <= 1e-13);
}

// ---------------------------------------------------------------- timings

fn time_ms<T>(reps: usize, mut f: impl FnMut() -> T) -> f64 {
    let _ = f();
    let t = std::time::Instant::now();
    for _ in 0..reps {
        std::hint::black_box(f());
    }
    t.elapsed().as_secs_f64() * 1e3 / reps as f64
}

/// BCT build and one `Ā` product against the dense Gram assembly and the dense product, at the
/// default options, 1 and 64 threads. Freedman 600/1200/2400 and 2400 refined twice and four
/// times. Run alone; prints a table. `BCT_TIMING_SKIP_DENSE_ABOVE` (vertices, default 10000)
/// bounds the dense runs.
#[test]
#[ignore]
fn timing_bct_against_dense() {
    let ex = Exponents::new(2.0, 4.5).unwrap();
    let q = GramQuadrature::Midpoint;
    let layout = DofLayout::new(None);
    let skip: usize = std::env::var("BCT_TIMING_SKIP_DENSE_ABOVE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(10_000);
    let mut cases: Vec<(String, CurveNetwork, Vec<P>)> = Vec::new();
    for n in [600, 1200, 2400] {
        let (net, pos) =
            read_curve_obj(data(&format!("FreedmanHeWang/freedman-he-wang-{n}.obj"))).unwrap();
        cases.push((format!("freedman-{n}"), net, pos));
    }
    let (n1, p1) = (cases[2].1.clone(), cases[2].2.clone());
    let (n2, p2) = subdivide(&n1, &p1);
    let (n4, p4) = subdivide(&n2, &p2);
    cases.push(("freedman-2400 x2".into(), n2, p2));
    cases.push(("freedman-2400 x4".into(), n4, p4));
    eprintln!(
        "| curve | edges | threads | BCT build ms | Ā_h product ms | dense Gram ms | dense product ms | product speed-up | Ā err | adm. pairs % | near entries |"
    );
    for (name, net, pos) in &cases {
        for threads in [1, 64] {
            let p = pool(threads);
            let build_ms = time_ms(3, || {
                build(&p, net, pos, ex, q, layout, BctOptions::default())
            });
            let h = build(&p, net, pos, ex, q, layout, BctOptions::default());
            let x = rand_vec(h.n_dofs(), 1);
            let reps = if threads == 1 { 5 } else { 20 };
            let prod_ms = time_ms(reps, || h.apply(&x).unwrap());
            let (gram_ms, dprod_ms, err) = if net.n_vertices() <= skip {
                let t = std::time::Instant::now();
                let g = sobolev_gram(&p, net, pos, ex, q);
                let gram_ms = t.elapsed().as_secs_f64() * 1e3;
                let dprod_ms = time_ms(3, || dense_apply_par(&p, &g.a, 3, &x));
                let err = rel_err(&h.apply(&x).unwrap(), &dense_apply(&g.a, 3, &x));
                (gram_ms, dprod_ms, err)
            } else {
                (f64::NAN, f64::NAN, f64::NAN)
            };
            let s = h.stats();
            let m = net.n_edges() as f64;
            eprintln!(
                "| {name} | {} | {threads} | {build_ms:.2} | {prod_ms:.3} | {gram_ms:.1} | {dprod_ms:.2} | {:.1} | {err:.2e} | {:.1} | {} |",
                net.n_edges(),
                dprod_ms / prod_ms,
                100.0 * s.admissible_pairs as f64 / (m * m),
                s.near_entries
            );
        }
    }
}
