//! Tests of work package WP4 (`SPEC.md` §8.1–8.2, §13.2, §13.4, D15, D25, D26): the generic
//! BVH, and the Barnes–Hut energy and differential against the exact ones.
//!
//! The timing tests are `#[ignore]`d; run them alone, in release:
//! `cargo test -p ddg-repulsive-curves --release --test barnes_hut -- --ignored --nocapture timing`.

#![allow(clippy::needless_range_loop)]

use ddg_repulsive_curves::energy::barnes_hut::BhTree;
use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::spatial::aabb::Aabb;
use ddg_repulsive_curves::spatial::bvh::{Bvh, Element};
use ddg_repulsive_curves::{
    BhOptions, CurveNetwork, EnergyEval, Exponents, RcError, energy, read_curve_obj,
};

type P = [f64; 3];

// ---------------------------------------------------------------- helpers

fn data(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/repulsive-curves")
        .join(rel)
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

fn perturb(net: &CurveNetwork, pos: &[P], amount: f64, seed: u64) -> Vec<P> {
    let lbar = net.mean_edge_length(pos);
    let mut r = Lcg(seed);
    pos.iter()
        .map(|p| {
            [
                p[0] + amount * lbar * r.next(),
                p[1] + amount * lbar * r.next(),
                p[2] + amount * lbar * r.next(),
            ]
        })
        .collect()
}

/// Linear midpoint subdivision of every edge (the curve as a set is unchanged).
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

fn norm2(g: &[P]) -> f64 {
    g.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()
}

fn diff_err(a: &[P], b: &[P]) -> f64 {
    let d: Vec<P> = a
        .iter()
        .zip(b)
        .map(|(x, y)| [x[0] - y[0], x[1] - y[1], x[2] - y[2]])
        .collect();
    norm2(&d) / norm2(b)
}

fn pool(n: usize) -> rayon::ThreadPool {
    make_pool(Some(n)).unwrap()
}

/// `θ_T = ε`, `θ_x = ε/2` (the default is `ε = 0.25`).
fn opts_eps(eps: f64) -> BhOptions {
    BhOptions {
        theta_x: eps / 2.0,
        theta_t: eps,
        leaf_size: 8,
    }
}

struct Cmp {
    e_rel: f64,
    g_rel: f64,
    pairs: u64,
    lumped: u64,
}

fn compare(
    net: &CurveNetwork,
    pos: &[P],
    ex: Exponents,
    o: BhOptions,
    p: &rayon::ThreadPool,
) -> Cmp {
    let (e0, g0) = energy_and_differential(p, net, pos, ex);
    let tree = BhTree::build(net, pos, o).unwrap();
    let (e1, g1, st) = tree.energy_and_differential(p, ex);
    Cmp {
        e_rel: ((e1 - e0) / e0).abs(),
        g_rel: diff_err(&g1, &g0),
        pairs: st.pairs_exact,
        lumped: st.nodes_lumped,
    }
}

// ---------------------------------------------------------------- options

#[test]
fn bad_options_are_errors_not_panics() {
    let (net, pos) = trefoil(40);
    let ok = BhOptions::default();
    for bad in [
        BhOptions { theta_x: 0.0, ..ok },
        BhOptions {
            theta_x: -1.0,
            ..ok
        },
        BhOptions {
            theta_x: f64::NAN,
            ..ok
        },
        BhOptions { theta_t: 0.0, ..ok },
        BhOptions {
            theta_t: f64::INFINITY,
            ..ok
        },
        BhOptions { leaf_size: 0, ..ok },
    ] {
        assert!(
            matches!(bad.validate(), Err(RcError::InvalidOption { .. })),
            "{bad:?}"
        );
        assert!(matches!(
            BhTree::build(&net, &pos, bad),
            Err(RcError::InvalidOption { .. })
        ));
        assert!(matches!(
            energy(
                &net,
                &pos,
                Exponents::default(),
                &EnergyEval::BarnesHut(bad),
                Some(1)
            ),
            Err(RcError::InvalidOption { .. })
        ));
    }
    // The error message names the option.
    let e = BhOptions { leaf_size: 0, ..ok }.validate().unwrap_err();
    assert!(e.to_string().contains("leaf_size"), "{e}");
    // Generic BVH.
    let none: Vec<Pt<3>> = Vec::new();
    assert!(Bvh::<3>::build(&none, 4).is_err());
    let one = vec![Pt::<3> {
        p: [0.0; 3],
        r: 0.0,
        m: 1.0,
    }];
    assert!(matches!(
        Bvh::<3>::build(&one, 0),
        Err(RcError::InvalidOption { .. })
    ));
    // Wrong position count on refit.
    let mut tree = BhTree::build(&net, &pos, ok).unwrap();
    assert!(tree.refit(&pos[..3]).is_err());
}

// ---------------------------------------------------------------- convergence to exact

#[test]
fn converges_to_exact_as_eps_goes_to_zero() {
    let (net, pos) = trefoil(240);
    let pos = perturb(&net, &pos, 0.05, 7);
    let p = pool(4);
    for ex in [
        Exponents::new(3.0, 6.0).unwrap(),
        Exponents::new(2.0, 4.5).unwrap(),
    ] {
        let mut prev = (f64::INFINITY, f64::INFINITY);
        let mut errs = Vec::new();
        for eps in [4.0, 2.0, 1.0, 0.5, 0.25, 0.125] {
            let c = compare(&net, &pos, ex, opts_eps(eps), &p);
            eprintln!(
                "trefoil240 ({},{}) eps {eps:<7} E err {:.2e}  g err {:.2e}  exact pairs {}  lumped {}",
                ex.alpha(),
                ex.beta(),
                c.e_rel,
                c.g_rel,
                c.pairs,
                c.lumped
            );
            errs.push((eps, c.e_rel, c.g_rel));
            // The differential error is roughly monotone in eps (3x slack: the tree changes with
            // eps). The signed energy error may cross zero by luck (eps = 2 below is better than
            // eps = 1), so the energy is only checked at the ends.
            assert!(c.g_rel <= 3.0 * prev.1 + 1e-14, "g err not decreasing");
            prev = (c.e_rel, c.g_rel);
        }
        let first = errs[0];
        let last = *errs.last().unwrap();
        assert!(
            last.1 < 0.02 * first.1 && last.2 < 0.05 * first.2,
            "{errs:?}"
        );
        // Tiny eps: nothing is admissible; equal to the exact values up to the summation order.
        let c = compare(&net, &pos, ex, opts_eps(1e-12), &p);
        assert_eq!(c.lumped, 0);
        assert!(
            c.e_rel < 1e-13 && c.g_rel < 1e-13,
            "{} {}",
            c.e_rel,
            c.g_rel
        );
    }
}

#[test]
fn all_inadmissible_settings_equal_exact() {
    let (net, pos) = trefoil(150);
    let pos = perturb(&net, &pos, 0.05, 3);
    let ex = Exponents::default();
    let p = pool(3);
    let tiny = 1e-300;
    let cases = [
        // spatial criterion never holds
        BhOptions {
            theta_x: tiny,
            theta_t: 10.0,
            leaf_size: 4,
        },
        // tangent criterion never holds
        BhOptions {
            theta_x: 10.0,
            theta_t: tiny,
            leaf_size: 4,
        },
        BhOptions {
            theta_x: tiny,
            theta_t: tiny,
            leaf_size: 1,
        },
        // the root is a leaf: everything is exact whatever the thetas
        BhOptions {
            theta_x: 10.0,
            theta_t: 10.0,
            leaf_size: 150,
        },
    ];
    for o in cases {
        let c = compare(&net, &pos, ex, o, &p);
        assert_eq!(c.lumped, 0, "{o:?}");
        assert_eq!(c.pairs, 150 * 147, "{o:?}"); // all ordered pairs without a common vertex
        assert!(
            c.e_rel < 1e-13 && c.g_rel < 1e-13,
            "{o:?}: {} {}",
            c.e_rel,
            c.g_rel
        );
    }
}

/// Error against `ε` on Freedman–He–Wang 600/1200/2400 (exponents (2, 4.5)) and a perturbed
/// Knot128 curve ((3, 6)); the table is printed and the default is checked.
#[test]
fn error_against_eps_on_datasets() {
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
    if !data("").exists() {
        eprintln!("data/repulsive-curves not present; skipping");
        return;
    }
    let p = pool(16);
    for (name, rel, ex) in files {
        let (net, pos) = read_curve_obj(data(rel)).unwrap();
        let (e0, g0) = energy_and_differential(&p, &net, &pos, ex);
        eprintln!("{name}: {} edges, exact energy {e0:.6e}", net.n_edges());
        let mut at_default = None;
        for eps in [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03] {
            let tree = BhTree::build(&net, &pos, opts_eps(eps)).unwrap();
            let (e1, g1, st) = tree.energy_and_differential(&p, ex);
            let (er, gr) = (((e1 - e0) / e0).abs(), diff_err(&g1, &g0));
            eprintln!(
                "  eps {eps:<7} E err {er:.2e}   g err {gr:.2e}   exact pairs {:>9}  lumped {:>8}",
                st.pairs_exact, st.nodes_lumped
            );
            if eps == 0.25 {
                at_default = Some((er, gr));
            }
        }
        let (er, gr) = at_default.unwrap();
        assert!(
            er < 2e-2 && gr < 5e-2,
            "{name}: default eps errors {er} {gr}"
        );
    }
}

// ---------------------------------------------------------------- neighbour guard

/// For every edge and every lumped node, no edge of the node is the edge itself or shares a
/// vertex with it, and lumped nodes plus exact leaf edges cover each `J` with `(I, J) ∈ 𝒫`
/// exactly once. Huge thetas make the guard the only thing that keeps nodes from being lumped.
#[test]
fn neighbour_guard_excludes_adjacent_pairs_exactly() {
    let (net, pos) = trefoil(96);
    let pos = perturb(&net, &pos, 0.02, 11);
    let (net2, pos2) = subdivide(&net, &pos);
    for (net, pos) in [(&net, &pos), (&net2, &pos2)] {
        for leaf_size in [1usize, 2, 8] {
            let o = BhOptions {
                theta_x: 1e6,
                theta_t: 1e6,
                leaf_size,
            };
            let tree = BhTree::build(net, pos, o).unwrap();
            let m = net.n_edges();
            let mut n_lumped = 0;
            for i in 0..m {
                let (lumped, exact) = tree.classify(i);
                n_lumped += lumped.len();
                let mut covered = vec![0u32; m];
                for n in lumped {
                    for &j in tree.bvh().elements_of(n) {
                        let j = j as usize;
                        assert!(
                            !net.edges_adjacent(i, j),
                            "edge {i} lumps adjacent edge {j}"
                        );
                        covered[j] += 1;
                    }
                }
                for j in exact {
                    assert!(!net.edges_adjacent(i, j));
                    covered[j] += 1;
                }
                for j in 0..m {
                    let want = u32::from(!net.edges_adjacent(i, j));
                    assert_eq!(covered[j], want, "edge {i}, partner {j}, leaf {leaf_size}");
                }
            }
            assert!(
                n_lumped > 0,
                "the guard test lumps nothing; it tests nothing"
            );
        }
    }
}

/// On a short closed curve (6 edges, `P` has 6 * 3 ordered pairs) the Barnes–Hut energy with
/// every node open equals the exact one, and with lumping allowed it stays finite and close.
#[test]
fn short_curve_matches_exact() {
    let n = 6usize;
    let pos: Vec<P> = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [t.cos(), t.sin(), 0.3 * (2.0 * t).sin()]
        })
        .collect();
    let net = CurveNetwork::new(n, closed_edges(n)).unwrap();
    let ex = Exponents::default();
    let p = pool(2);
    for leaf_size in [1usize, 2, 6, 8] {
        let o = BhOptions {
            theta_x: 1e-300,
            theta_t: 1e-300,
            leaf_size,
        };
        let c = compare(&net, &pos, ex, o, &p);
        assert_eq!(c.lumped, 0);
        assert_eq!(c.pairs, 6 * 3);
        assert!(c.e_rel < 1e-14 && c.g_rel < 1e-14);
    }
    // Everything lumpable (but never an adjacent pair): the answer is finite.
    let o = BhOptions {
        theta_x: 1e6,
        theta_t: 1e6,
        leaf_size: 1,
    };
    let c = compare(&net, &pos, ex, o, &p);
    assert!(c.lumped > 0 && c.e_rel.is_finite() && c.g_rel.is_finite());
}

// ---------------------------------------------------------------- determinism, refit

fn bits(g: &[P]) -> Vec<u64> {
    g.iter().flatten().map(|x| x.to_bits()).collect()
}

#[test]
fn bitwise_identical_at_1_7_and_64_threads() {
    let mut cases: Vec<(CurveNetwork, Vec<P>)> = Vec::new();
    let (n, p) = trefoil(700);
    cases.push((n.clone(), perturb(&n, &p, 0.05, 5)));
    if data("").exists() {
        cases.push(read_curve_obj(data("FreedmanHeWang/freedman-he-wang-600.obj")).unwrap());
    }
    for (net, pos) in cases {
        for o in [BhOptions::default(), opts_eps(0.7)] {
            let mut results = Vec::new();
            for t in [1usize, 7, 64] {
                let r = energy(
                    &net,
                    &pos,
                    Exponents::default(),
                    &EnergyEval::BarnesHut(o),
                    Some(t),
                )
                .unwrap();
                assert!(r.nodes_lumped > 0);
                results.push((
                    r.energy.to_bits(),
                    bits(&r.differential),
                    r.pairs_exact,
                    r.nodes_lumped,
                ));
            }
            assert_eq!(results[0], results[1], "1 vs 7 threads");
            assert_eq!(results[0], results[2], "1 vs 64 threads");
            // energy-only path
            let tree = BhTree::build(&net, &pos, o).unwrap();
            let es: Vec<u64> = [1usize, 7, 64]
                .iter()
                .map(|&t| tree.energy(&pool(t), Exponents::default()).0.to_bits())
                .collect();
            assert!(es[0] == es[1] && es[0] == es[2]);
            // the energy-only value equals the energy of the full evaluation
            assert_eq!(es[0], results[0].0);
        }
    }
}

#[test]
fn refit_equals_rebuild_for_unchanged_topology() {
    let (net, pos) = trefoil(300);
    let pos0 = perturb(&net, &pos, 0.05, 21);
    let ex = Exponents::default();
    let p = pool(4);
    let o = BhOptions::default();
    // (a) Refit on the same positions reproduces the build bit for bit.
    let fresh = BhTree::build(&net, &pos0, o).unwrap();
    let mut refit = BhTree::build(&net, &pos0, o).unwrap();
    refit.refit(&pos0).unwrap();
    assert_eq!(fresh.bvh(), refit.bvh());
    let (ea, ga, sa) = fresh.energy_and_differential(&p, ex);
    let (eb, gb, sb) = refit.energy_and_differential(&p, ex);
    assert_eq!(ea.to_bits(), eb.to_bits());
    assert_eq!(bits(&ga), bits(&gb));
    assert_eq!(sa, sb);
    // (b) Move the curve slightly: the refitted tree keeps the topology; its node boxes are those
    // a build over the new points would give on the same topology (checked by bounds), and the
    // result is an equally good approximation. With tiny eps it equals the exact value and the
    // rebuilt one to roundoff.
    let pos1 = perturb(&net, &pos0, 0.02, 22);
    let order_before = refit.bvh().order().to_vec();
    refit.refit(&pos1).unwrap();
    assert_eq!(refit.bvh().order(), &order_before[..]);
    let (e0, g0) = energy_and_differential(&p, &net, &pos1, ex);
    let rebuilt = BhTree::build(&net, &pos1, o).unwrap();
    let (e_rf, g_rf, _) = refit.energy_and_differential(&p, ex);
    let (e_rb, g_rb, _) = rebuilt.energy_and_differential(&p, ex);
    eprintln!(
        "after a move: refit E err {:.2e} g err {:.2e}; rebuild E err {:.2e} g err {:.2e}",
        ((e_rf - e0) / e0).abs(),
        diff_err(&g_rf, &g0),
        ((e_rb - e0) / e0).abs(),
        diff_err(&g_rb, &g0)
    );
    assert!(((e_rf - e0) / e0).abs() < 2e-2 && diff_err(&g_rf, &g0) < 5e-2);
    // Boxes of the refitted tree contain the new tangent-points of their edges.
    let geom = net.edge_geometry(&pos1);
    for (n, nd) in refit.bvh().nodes().iter().enumerate() {
        for &e in refit.bvh().elements_of(n) {
            let e = e as usize;
            let (x, t) = (geom.midpoint[e], geom.tangent[e]);
            let pt = [x[0], x[1], x[2], t[0], t[1], t[2]];
            assert!(nd.bbox.contains_point(pt));
        }
    }
    // Tiny eps: refit == rebuild == exact (to the summation order).
    let o = opts_eps(1e-12);
    let mut t = BhTree::build(&net, &pos0, o).unwrap();
    t.refit(&pos1).unwrap();
    let (e_t, g_t, _) = t.energy_and_differential(&p, ex);
    assert!(((e_t - e0) / e0).abs() < 1e-13 && diff_err(&g_t, &g0) < 1e-13);
}

// ---------------------------------------------------------------- finite differences

/// Central differences of the Barnes–Hut energy (the tree refitted, topology frozen: D26)
/// against (a) the Barnes–Hut differential and (b) the exact differential.
///
/// They are different approximations: the evaluated energy lumps only the pairs `(I, J)` with
/// `T_I` (Eq. 26) and its derivative also moves the lumped centres `x̄_N` and masses with the
/// vertices; the differential of §8.2 is the derivative of the lumped `(I, J) + (J, I)` pair with
/// `L_N, x̄_N, T̂_N` held fixed. Both converge to the exact differential with eps; the test
/// asserts the relative difference `|FD - g_BH|/|g_BH|` (2-norm of the whole vector) stays below
/// the tolerances listed in the loop and shrinks with eps. Measured on a perturbed trefoil with
/// 300 edges, leaf size 4, `eps = θ_T = 2 θ_x`: 0.19, 0.16, 0.041, 0.0035, 0.00087 for eps = 4, 2,
/// 1, 0.5, 0.25 (the BH differential differs from the exact one by 0.29, 0.22, 0.064,
/// 0.0024, 0.00033); the tolerances are about 1.5 times these (the first two 1.2 times). At eps -> 0 nothing is lumped and the
/// finite differences match the differential to the accuracy of the central differences.
#[test]
fn finite_differences_of_bh_energy_are_close_to_bh_differential() {
    let (net, pos) = trefoil(300);
    let pos = perturb(&net, &pos, 0.05, 9);
    let ex = Exponents::default();
    let p = pool(8);
    let h = 1e-6 * net.mean_edge_length(&pos);
    let (_, g_exact) = energy_and_differential(&p, &net, &pos, ex);
    let mut last = f64::INFINITY;
    for (eps, tol) in [
        (4.0, 0.29),
        (2.0, 0.24),
        (1.0, 0.062),
        (0.5, 5.3e-3),
        (0.25, 1.3e-3),
    ] {
        let o = BhOptions {
            leaf_size: 4,
            ..opts_eps(eps)
        };
        let mut tree = BhTree::build(&net, &pos, o).unwrap();
        let (_, g_bh, st) = tree.energy_and_differential(&p, ex);
        let mut fd = vec![[0.0; 3]; pos.len()];
        let mut q = pos.clone();
        for v in 0..pos.len() {
            for c in 0..3 {
                q[v][c] = pos[v][c] + h;
                tree.refit(&q).unwrap();
                let ep = tree.energy(&p, ex).0;
                q[v][c] = pos[v][c] - h;
                tree.refit(&q).unwrap();
                let em = tree.energy(&p, ex).0;
                q[v][c] = pos[v][c];
                fd[v][c] = (ep - em) / (2.0 * h);
            }
        }
        let (a, b) = (diff_err(&fd, &g_bh), diff_err(&fd, &g_exact));
        eprintln!(
            "eps {eps:<5} lumped {:>6}: |FD(BH E) - BH g|/|BH g| {a:.2e}; |FD(BH E) - exact g|/|g| {b:.2e}; |BH g - exact g|/|g| {:.2e}",
            st.nodes_lumped,
            diff_err(&g_bh, &g_exact)
        );
        assert!(a < tol, "eps {eps}: {a} >= {tol}");
        assert!(a <= last * 1.5 + 1e-9, "not shrinking with eps");
        last = a;
    }
    // With nothing lumped the two are the same function: FD against the differential to 1e-6.
    let o = opts_eps(1e-12);
    let mut tree = BhTree::build(&net, &pos, o).unwrap();
    let (_, g_bh, _) = tree.energy_and_differential(&p, ex);
    let mut worst: f64 = 0.0;
    let mut q = pos.clone();
    for v in (0..pos.len()).step_by(7) {
        for c in 0..3 {
            q[v][c] = pos[v][c] + h;
            tree.refit(&q).unwrap();
            let ep = tree.energy(&p, ex).0;
            q[v][c] = pos[v][c] - h;
            tree.refit(&q).unwrap();
            let em = tree.energy(&p, ex).0;
            q[v][c] = pos[v][c];
            worst = worst.max(((ep - em) / (2.0 * h) - g_bh[v][c]).abs() / norm2(&g_bh));
        }
    }
    assert!(worst < 1e-6, "{worst}");
}

// ---------------------------------------------------------------- the generic BVH

#[derive(Clone)]
struct Pt<const D: usize> {
    p: [f64; D],
    r: f64,
    m: f64,
}

impl<const D: usize> Element<D> for Pt<D> {
    fn point(&self) -> [f64; D] {
        self.p
    }
    fn bounds(&self) -> Aabb<D> {
        let mut b = Aabb::empty();
        b.grow(std::array::from_fn(|k| self.p[k] - self.r));
        b.grow(std::array::from_fn(|k| self.p[k] + self.r));
        b
    }
    fn mass(&self) -> f64 {
        self.m
    }
}

fn random_points<const D: usize>(n: usize, seed: u64, clustered: bool) -> Vec<Pt<D>> {
    let mut r = Lcg(seed);
    (0..n)
        .map(|i| {
            let c = if clustered && i % 3 == 0 { 0.01 } else { 1.0 };
            Pt {
                p: std::array::from_fn(|_| c * r.next() * 10.0),
                r: if i % 5 == 0 {
                    0.0
                } else {
                    (r.next() + 0.5) * 0.3
                },
                m: 0.1 + r.next() + 0.5,
            }
        })
        .collect()
}

fn check_bvh<const D: usize>(bvh: &Bvh<D>, el: &[Pt<D>], leaf_size: usize, min_fraction: f64) {
    let n = el.len();
    // order is a permutation
    let mut seen = vec![0u32; n];
    for &e in bvh.order() {
        seen[e as usize] += 1;
    }
    assert!(seen.iter().all(|&c| c == 1));
    // every element is in exactly one leaf; leaves respect the leaf size; boxes contain elements
    let mut in_leaf = vec![0u32; n];
    for (i, nd) in bvh.nodes().iter().enumerate() {
        assert!(nd.end > nd.start);
        if nd.is_leaf() {
            assert!(
                nd.len() <= leaf_size,
                "leaf with {} > {leaf_size}",
                nd.len()
            );
            for &e in bvh.elements_of(i) {
                in_leaf[e as usize] += 1;
            }
        } else {
            let (l, r) = (&bvh.nodes()[i + 1], &bvh.nodes()[nd.right as usize]);
            assert!(nd.len() > leaf_size);
            let lo = ((min_fraction * nd.len() as f64).ceil() as usize).clamp(1, nd.len() / 2);
            assert!(
                l.len() >= lo && r.len() >= lo,
                "unbalanced split {} + {}",
                l.len(),
                r.len()
            );
            assert_eq!(l.start, nd.start);
            assert_eq!(l.end, r.start);
            assert_eq!(r.end, nd.end);
            assert!(nd.bbox.contains_box(&l.bbox) && nd.bbox.contains_box(&r.bbox));
            assert!((nd.mass - (l.mass + r.mass)).abs() <= 1e-12 * nd.mass);
            assert!(nd.right as usize > i + 1);
        }
        let mut mass = 0.0;
        for &e in bvh.elements_of(i) {
            assert!(
                nd.bbox.contains_box(&el[e as usize].bounds()),
                "node {i} misses element {e}"
            );
            mass += el[e as usize].m;
        }
        assert!((nd.mass - mass).abs() <= 1e-12 * mass);
        let c = nd.centroid();
        for k in 0..D {
            let want: f64 = bvh
                .elements_of(i)
                .iter()
                .map(|&e| el[e as usize].m * el[e as usize].p[k])
                .sum::<f64>()
                / mass;
            assert!((c[k] - want).abs() <= 1e-10 * (1.0 + want.abs()));
        }
    }
    assert!(
        in_leaf.iter().all(|&c| c == 1),
        "an element is in zero or several leaves"
    );
    assert_eq!(bvh.nodes()[0].len(), n);
}

fn property<const D: usize>() {
    for (seed, n, leaf, clustered) in [
        (1u64, 1usize, 1usize, false),
        (2, 2, 1, false),
        (3, 7, 8, false),
        (4, 100, 1, false),
        (5, 100, 3, true),
        (6, 1000, 8, true),
        (7, 1000, 8, false),
        (8, 513, 5, true),
    ] {
        for frac in [0.0, 0.1, 0.25, 0.5] {
            let el: Vec<Pt<D>> = random_points(n, seed * 100 + D as u64, clustered);
            let bvh = Bvh::build_with(&el, leaf, frac).unwrap();
            check_bvh(&bvh, &el, leaf, frac);
            // Refit after moving every element: same topology, new valid boxes.
            let moved: Vec<Pt<D>> = random_points(n, seed * 977 + D as u64, clustered);
            let mut b2 = bvh.clone();
            b2.refit(&moved).unwrap();
            assert_eq!(b2.order(), bvh.order());
            check_bvh(&b2, &moved, leaf, frac);
            // Refit on unchanged elements is the identity.
            let mut b3 = bvh.clone();
            b3.refit(&el).unwrap();
            assert_eq!(b3, bvh);
            // fold_up: element counts agree with the node ranges.
            let counts = bvh.fold_up(|idx| idx.len(), |a, b| a + b);
            for (c, nd) in counts.iter().zip(bvh.nodes()) {
                assert_eq!(*c, nd.len());
            }
            // Determinism: building twice gives the same tree.
            assert_eq!(Bvh::build_with(&el, leaf, frac).unwrap(), bvh);
        }
    }
}

#[test]
fn bvh_property_test_d2() {
    property::<2>();
}
#[test]
fn bvh_property_test_d3() {
    property::<3>();
}
#[test]
fn bvh_property_test_d6() {
    property::<6>();
}

#[test]
fn bvh_with_identical_points_terminates() {
    let el: Vec<Pt<6>> = (0..200)
        .map(|_| Pt {
            p: [1.0; 6],
            r: 0.0,
            m: 1.0,
        })
        .collect();
    for frac in [0.0, 0.25] {
        let bvh = Bvh::build_with(&el, 8, frac).unwrap();
        check_bvh(&bvh, &el, 8, frac);
    }
}

fn depth_of<const D: usize>(bvh: &Bvh<D>) -> usize {
    let nodes = bvh.nodes();
    let mut depth = vec![0usize; nodes.len()];
    let mut max = 0;
    for i in 0..nodes.len() {
        if nodes[i].is_leaf() {
            max = max.max(depth[i]);
        } else {
            depth[i + 1] = depth[i] + 1;
            depth[nodes[i].right as usize] = depth[i] + 1;
        }
    }
    max
}

/// SPEC D15 as printed (split at the minimum of the sum of squared diameters over all positions)
/// builds a degenerate tree on a tangled curve: it peels single tangent-points off the root. The
/// default bounds the imbalance and gives depth `<= log_{4/3}(n / leaf) + 1`.
#[test]
fn unrestricted_split_rule_is_degenerate_on_curves() {
    if !data("").exists() {
        return;
    }
    let (net, pos) = read_curve_obj(data("FreedmanHeWang/freedman-he-wang-1200.obj")).unwrap();
    let o = BhOptions::default();
    let literal = BhTree::build_with(&net, &pos, o, 0.0).unwrap();
    let balanced = BhTree::build(&net, &pos, o).unwrap();
    let (dl, db) = (depth_of(literal.bvh()), depth_of(balanced.bvh()));
    let bound = ((net.n_edges() as f64 / 8.0).ln() / (4.0f64 / 3.0).ln()).ceil() as usize + 1;
    eprintln!("freedman-1200: depth {dl} with the printed rule, {db} balanced (bound {bound})");
    assert!(db <= bound, "{db} > {bound}");
    assert!(
        dl > 4 * db,
        "the printed rule was expected to degenerate: {dl} vs {db}"
    );
}

/// The experimental lumping of admissible leaves stays an approximation of the same quantities.
#[test]
fn lumping_leaves_is_a_looser_approximation_of_the_same_thing() {
    let (net, pos) = trefoil(400);
    let pos = perturb(&net, &pos, 0.05, 13);
    let ex = Exponents::default();
    let p = pool(4);
    let (e0, g0) = energy_and_differential(&p, &net, &pos, ex);
    let mut tree = BhTree::build(&net, &pos, opts_eps(0.25)).unwrap();
    let (_, _, a) = tree.energy_and_differential(&p, ex);
    tree.set_lump_leaves(true);
    let (e1, g1, b) = tree.energy_and_differential(&p, ex);
    assert!(b.pairs_exact < a.pairs_exact);
    assert!(((e1 - e0) / e0).abs() < 1e-2 && diff_err(&g1, &g0) < 5e-2);
}

// ---------------------------------------------------------------- timings (ignored)

fn time_ms<T>(reps: usize, mut f: impl FnMut() -> T) -> f64 {
    let _ = f();
    let t = std::time::Instant::now();
    for _ in 0..reps {
        std::hint::black_box(f());
    }
    t.elapsed().as_secs_f64() * 1e3 / reps as f64
}

/// Barnes-Hut (build + energy + differential) against exact, default options, 1 and 64 threads.
/// Run alone; prints a table. `BH_TIMING_SKIP_EXACT_ABOVE` (edges, default 12000) bounds the
/// exact runs.
#[test]
#[ignore]
fn timing_bh_against_exact() {
    let ex = Exponents::new(2.0, 4.5).unwrap();
    let base = read_curve_obj(data("FreedmanHeWang/freedman-he-wang-2400.obj")).unwrap();
    let mut cases: Vec<(String, CurveNetwork, Vec<P>)> = Vec::new();
    for n in [600, 1200, 2400] {
        let (net, pos) =
            read_curve_obj(data(&format!("FreedmanHeWang/freedman-he-wang-{n}.obj"))).unwrap();
        cases.push((format!("freedman-{n}"), net, pos));
    }
    let (n2, p2) = subdivide(&base.0, &base.1);
    let (n4, p4) = subdivide(&n2, &p2);
    cases.push(("freedman-2400 x2 (4800)".into(), n2, p2));
    cases.push(("freedman-2400 x4 (9600)".into(), n4, p4));
    eprintln!(
        "case, edges, threads, exact ms, BH ms (build+eval), speed-up, E err, g err, exact pairs, lumped"
    );
    for (name, net, pos) in cases {
        for threads in [1usize, 64] {
            let p = pool(threads);
            let reps = if net.n_edges() <= 2400 { 3 } else { 1 };
            let t_exact = time_ms(reps, || energy_and_differential(&p, &net, &pos, ex));
            let o = BhOptions::default();
            let t_bh = time_ms(reps.max(3), || {
                let tree = BhTree::build(&net, &pos, o).unwrap();
                tree.energy_and_differential(&p, ex)
            });
            let tree = BhTree::build(&net, &pos, o).unwrap();
            let t_eval = time_ms(reps.max(3), || tree.energy_and_differential(&p, ex));
            let t_build = time_ms(reps.max(3), || BhTree::build(&net, &pos, o).unwrap());
            let (e0, g0) = energy_and_differential(&p, &net, &pos, ex);
            let (e1, g1, st) = tree.energy_and_differential(&p, ex);
            eprintln!(
                "TIMING {name}, {} edges, {threads} threads: exact {t_exact:.1} ms; BH {t_bh:.1} ms (build {t_build:.1}, eval {t_eval:.1}); x{:.1}; E err {:.1e}, g err {:.1e}; exact pairs {}, lumped {}",
                net.n_edges(),
                t_exact / t_bh,
                ((e1 - e0) / e0).abs(),
                diff_err(&g1, &g0),
                st.pairs_exact,
                st.nodes_lumped
            );
        }
    }
}
