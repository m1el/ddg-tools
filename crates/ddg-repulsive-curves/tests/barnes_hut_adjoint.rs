//! Tests of work package WP3d (`SPEC.md` §17i, D25, D26): the Barnes–Hut differential is the
//! exact derivative (discrete adjoint) of the evaluated Barnes–Hut energy on fixed interaction
//! lists ([`BhTree::energy_with_lists`]).
//!
//! * central differences (Richardson-extrapolated) of `energy_with_lists` against the differential:
//!   coordinate by coordinate on small curves, along random directions and along the differential
//!   itself on larger ones; closed curve, open curve, a network with junctions of degree 3 and 4;
//!   `ε = θ_T = 2θ_x ∈ {0.25, 0.5}` and a coarser `θ_x = θ_T = 0.5`; at the base curve and at a
//!   moved curve with the base curve's lists (the line search's situation);
//! * `θ → 0` (nothing lumped) gives the exact differential;
//! * bitwise determinism across 1, 3 and 8 threads;
//! * mismatched lists are errors, not panics;
//! * (ignored, release) the Freedman–He–Wang 600 Barnes–Hut flow converges to `‖g‖ ≤ 1e-4`:
//!   `cargo test --release -p ddg-repulsive-curves --test barnes_hut_adjoint -- --ignored --nocapture`.

#![allow(clippy::needless_range_loop)]

use std::time::Instant;

use ddg_repulsive_curves::energy::barnes_hut::{BhTree, InteractionLists};
use ddg_repulsive_curves::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_curves::step::stop::Status;
use ddg_repulsive_curves::step::{Backend, FlowOptions, Session, Solver};
use ddg_repulsive_curves::{
    BhOptions, CurveNetwork, EnergyEval, Exponents, RcError, read_curve_obj,
};

type P = [f64; 3];

// ---------------------------------------------------------------- curves

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

/// An open helix of `n` vertices, 1.5 turns (endpoints of degree 1).
fn helix(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = 1.5 * std::f64::consts::TAU * k as f64 / (n - 1) as f64;
            [t.cos(), t.sin(), 0.4 * t]
        })
        .collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

/// Four meridian arcs between the poles `A = (0, 0, 1)` and `B = (0, 0, −1)` (junctions of
/// degree 4), each with `per_arc` interior vertices, plus an open radial branch of `branch`
/// edges from the middle vertex of the first arc (a junction of degree 3).
fn network(per_arc: usize, branch: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos: Vec<P> = vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
    let mut edges = Vec::new();
    let mut mid_of_first = 0u32;
    for arc in 0..4 {
        let phi = std::f64::consts::FRAC_PI_2 * arc as f64 + 0.1;
        let mut prev = 0u32;
        for k in 1..=per_arc {
            let t = std::f64::consts::PI * k as f64 / (per_arc + 1) as f64;
            let r = t.sin();
            pos.push([r * phi.cos(), r * phi.sin(), t.cos()]);
            let v = pos.len() as u32 - 1;
            edges.push([prev, v]);
            prev = v;
            if arc == 0 && k == per_arc.div_ceil(2) {
                mid_of_first = v;
            }
        }
        edges.push([prev, 1]);
    }
    let m = pos[mid_of_first as usize];
    let h = std::f64::consts::PI / (per_arc + 1) as f64;
    let mut prev = mid_of_first;
    for k in 1..=branch {
        let s = 1.0 + h * k as f64;
        pos.push([m[0] * s, m[1] * s, m[2] + 0.3 * h * k as f64]);
        let v = pos.len() as u32 - 1;
        edges.push([prev, v]);
        prev = v;
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

fn random_dir(n: usize, seed: u64) -> Vec<P> {
    let mut r = Lcg(seed);
    (0..n).map(|_| [r.next(), r.next(), r.next()]).collect()
}

fn norm2(g: &[P]) -> f64 {
    g.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()
}

fn dotv(a: &[P], b: &[P]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| x[0] * y[0] + x[1] * y[1] + x[2] * y[2])
        .sum()
}

fn diff_err(a: &[P], b: &[P]) -> f64 {
    let d: Vec<P> = a
        .iter()
        .zip(b)
        .map(|(x, y)| [x[0] - y[0], x[1] - y[1], x[2] - y[2]])
        .collect();
    norm2(&d) / norm2(b)
}

fn axpy_pos(pos: &[P], s: f64, v: &[P]) -> Vec<P> {
    pos.iter()
        .zip(v)
        .map(|(p, d)| [p[0] + s * d[0], p[1] + s * d[1], p[2] + s * d[2]])
        .collect()
}

fn pool(n: usize) -> rayon::ThreadPool {
    make_pool(Some(n)).unwrap()
}

/// `θ_T = ε`, `θ_x = ε/2` (the default is `ε = 0.25`).
fn opts_eps(eps: f64, leaf_size: usize) -> BhOptions {
    BhOptions {
        theta_x: eps / 2.0,
        theta_t: eps,
        leaf_size,
    }
}

fn bits(g: &[P]) -> Vec<u64> {
    g.iter().flatten().map(|x| x.to_bits()).collect()
}

// ---------------------------------------------------------------- finite differences

/// `Ê_bh` on fixed lists at `pos` (the tree refitted, D26).
fn e_lists(
    tree: &BhTree,
    lists: &InteractionLists,
    pos: &[P],
    p: &rayon::ThreadPool,
    ex: Exponents,
) -> f64 {
    let mut t = tree.clone();
    t.refit(pos).unwrap();
    t.energy_with_lists(lists, p, ex)
}

/// Directional derivative of `Ê_bh` on fixed lists along `v` at `pos`: central differences at
/// `h` and `h/2`, Richardson-extrapolated (error `O(h⁴)`).
fn fd_dir(
    tree: &BhTree,
    lists: &InteractionLists,
    pos: &[P],
    v: &[P],
    h: f64,
    p: &rayon::ThreadPool,
    ex: Exponents,
) -> f64 {
    let d = |h: f64| {
        (e_lists(tree, lists, &axpy_pos(pos, h, v), p, ex)
            - e_lists(tree, lists, &axpy_pos(pos, -h, v), p, ex))
            / (2.0 * h)
    };
    (4.0 * d(0.5 * h) - d(h)) / 3.0
}

struct FdResult {
    /// Worst `|FD − ⟨g, v⟩| / (‖g‖ ‖v‖)` over random `v`.
    random: f64,
    /// `|FD − ⟨g, g⟩| / ‖g‖²` along `v = g`.
    along_g: f64,
    /// The same two numbers for the old (lumped-pairs) differential.
    old_random: f64,
    old_along_g: f64,
    lumped: u64,
    exact: u64,
}

/// FD checks of the differential at `at` (the tree built at `base`, its lists frozen).
fn fd_check(
    net: &CurveNetwork,
    base: &[P],
    at: &[P],
    o: BhOptions,
    ex: Exponents,
    p: &rayon::ThreadPool,
    n_dirs: usize,
) -> FdResult {
    let tree0 = BhTree::build(net, base, o).unwrap();
    let lists = tree0.interaction_lists(p);
    let mut tree = tree0.clone();
    tree.refit(at).unwrap();
    let (e, g, st) = tree
        .energy_and_differential_with_lists(&lists, p, ex)
        .unwrap();
    // The energy is energy_with_lists, bit for bit.
    assert_eq!(e.to_bits(), tree.energy_with_lists(&lists, p, ex).to_bits());
    let (_, g_old, _) = tree.energy_and_differential_lumped_pairs(p, ex);
    let h = 1e-3 * net.mean_edge_length(at);
    let mut res = FdResult {
        random: 0.0,
        along_g: 0.0,
        old_random: 0.0,
        old_along_g: 0.0,
        lumped: st.nodes_lumped,
        exact: st.pairs_exact,
    };
    for k in 0..n_dirs {
        let v = random_dir(at.len(), 1000 + k as u64);
        let vmax = v.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()));
        let fd = fd_dir(&tree0, &lists, at, &v, h / vmax, p, ex);
        let scale = norm2(&g) * norm2(&v);
        res.random = res.random.max((fd - dotv(&g, &v)).abs() / scale);
        res.old_random = res.old_random.max((fd - dotv(&g_old, &v)).abs() / scale);
    }
    let gmax = g.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()));
    let fd = fd_dir(&tree0, &lists, at, &g, h / gmax, p, ex);
    let gg = dotv(&g, &g);
    res.along_g = (fd - gg).abs() / gg;
    res.old_along_g = (fd - dotv(&g_old, &g)).abs() / gg;
    res
}

fn report(name: &str, r: &FdResult) {
    eprintln!(
        "{name:<44} lumped {:>6} exact {:>7} | adjoint: rand {:.1e} along g {:.1e} | old: rand {:.1e} along g {:.1e}",
        r.lumped, r.exact, r.random, r.along_g, r.old_random, r.old_along_g
    );
}

/// The tolerance of the task: ~1e-7 relative.
const TOL: f64 = 1e-7;

#[test]
fn finite_differences_match_the_adjoint_on_curves_and_networks() {
    let p = pool(8);
    let ex36 = Exponents::new(3.0, 6.0).unwrap();
    let ex245 = Exponents::new(2.0, 4.5).unwrap();
    let curves: Vec<(&str, CurveNetwork, Vec<P>)> = {
        let (n1, p1) = trefoil(240);
        let p1 = perturb(&n1, &p1, 0.1, 3);
        let (n2, p2) = helix(300);
        let p2 = perturb(&n2, &p2, 0.1, 4);
        let (n3, p3) = network(40, 12);
        let p3 = perturb(&n3, &p3, 0.1, 5);
        vec![
            ("trefoil 240", n1, p1),
            ("open helix 300", n2, p2),
            ("network 4 arcs + branch", n3, p3),
        ]
    };
    let settings = [
        ("eps 0.25", opts_eps(0.25, 4)),
        ("eps 0.5", opts_eps(0.5, 4)),
        (
            "theta_x = theta_t = 0.5",
            BhOptions {
                theta_x: 0.5,
                theta_t: 0.5,
                leaf_size: 4,
            },
        ),
    ];
    for (cname, net, pos) in &curves {
        for (oname, o) in &settings {
            for ex in [ex36, ex245] {
                let name = format!("{cname} {oname} ({},{})", ex.alpha(), ex.beta());
                let r = fd_check(net, pos, pos, *o, ex, &p, 4);
                report(&name, &r);
                assert!(r.lumped > 0, "{name}: nothing lumped");
                assert!(r.random < TOL && r.along_g < TOL, "{name}");
            }
        }
    }
}

/// The line search's situation: the lists of the base curve, the energy and differential at a
/// moved curve.
#[test]
fn finite_differences_match_with_frozen_lists_at_a_moved_curve() {
    let p = pool(8);
    let ex = Exponents::default();
    let (net, pos) = network(40, 12);
    let base = perturb(&net, &pos, 0.1, 6);
    let moved = perturb(&net, &base, 0.3, 7);
    for o in [opts_eps(0.25, 4), opts_eps(0.5, 8)] {
        let r = fd_check(&net, &base, &moved, o, ex, &p, 4);
        report("network, moved curve, base lists", &r);
        assert!(r.lumped > 0);
        assert!(r.random < TOL && r.along_g < TOL);
    }
}

/// Coordinate by coordinate on a small curve: `‖FD − g‖/‖g‖`.
#[test]
fn finite_differences_coordinatewise_on_a_small_network() {
    let p = pool(8);
    let ex = Exponents::new(2.0, 4.5).unwrap();
    let (net, pos) = network(16, 6);
    let pos = perturb(&net, &pos, 0.1, 8);
    let o = BhOptions {
        theta_x: 0.5,
        theta_t: 0.5,
        leaf_size: 2,
    };
    let tree = BhTree::build(&net, &pos, o).unwrap();
    let lists = tree.interaction_lists(&p);
    let (_, g, st) = tree.energy_and_differential(&p, ex);
    let (_, g_old, _) = tree.energy_and_differential_lumped_pairs(&p, ex);
    assert!(st.nodes_lumped > 0);
    let h = 1e-3 * net.mean_edge_length(&pos);
    let mut fd = vec![[0.0; 3]; pos.len()];
    for v in 0..pos.len() {
        for c in 0..3 {
            let mut e = vec![[0.0; 3]; pos.len()];
            e[v][c] = 1.0;
            fd[v][c] = fd_dir(&tree, &lists, &pos, &e, h, &p, ex);
        }
    }
    let (a, b) = (diff_err(&g, &fd), diff_err(&g_old, &fd));
    eprintln!(
        "network {} vertices, lumped {}: |g - FD|/|FD| adjoint {a:.2e}, old {b:.2e}",
        pos.len(),
        st.nodes_lumped
    );
    assert!(a < TOL, "{a}");
}

// ---------------------------------------------------------------- limits and identities

/// Nothing admissible: the differential is the exact one (up to the summation order).
#[test]
fn tiny_theta_gives_the_exact_differential() {
    let p = pool(4);
    for (net, pos) in [trefoil(150), network(20, 5)] {
        let pos = perturb(&net, &pos, 0.05, 9);
        for ex in [Exponents::default(), Exponents::new(2.0, 4.5).unwrap()] {
            let (e0, g0) = energy_and_differential(&p, &net, &pos, ex);
            let tree = BhTree::build(&net, &pos, opts_eps(1e-12, 8)).unwrap();
            let (e1, g1, st) = tree.energy_and_differential(&p, ex);
            assert_eq!(st.nodes_lumped, 0);
            let (er, gr) = (((e1 - e0) / e0).abs(), diff_err(&g1, &g0));
            eprintln!("theta 1e-12: E err {er:.1e}, g err {gr:.1e}");
            assert!(er < 1e-13 && gr < 1e-13, "{er} {gr}");
        }
    }
}

/// Translation invariance of the Barnes–Hut energy on fixed lists implies `Σ_v g_v = 0`; the
/// downward pass is re-centred at every node, so this holds far from the origin too.
#[test]
fn differential_sums_to_zero_far_from_the_origin() {
    let p = pool(4);
    let (net, pos) = trefoil(300);
    let pos = perturb(&net, &pos, 0.1, 10);
    let far: Vec<P> = pos
        .iter()
        .map(|q| [q[0] + 1e4, q[1] - 3e3, q[2] + 7e3])
        .collect();
    for at in [&pos, &far] {
        let tree = BhTree::build(&net, at, opts_eps(0.5, 4)).unwrap();
        let (_, g, st) = tree.energy_and_differential(&p, Exponents::default());
        assert!(st.nodes_lumped > 0);
        let mut s = [0.0; 3];
        for v in &g {
            for c in 0..3 {
                s[c] += v[c];
            }
        }
        let rel = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt() / norm2(&g);
        eprintln!("|sum g| / |g| = {rel:.1e}");
        assert!(rel < 1e-12, "{rel}");
    }
}

/// The adjoint at the base equals the adjoint with the tree's own lists, and the old
/// differential and the adjoint approximate the same exact differential.
#[test]
fn adjoint_is_close_to_the_exact_differential() {
    let p = pool(8);
    let (net, pos) = trefoil(400);
    let pos = perturb(&net, &pos, 0.05, 11);
    let ex = Exponents::default();
    let (_, g0) = energy_and_differential(&p, &net, &pos, ex);
    for eps in [1.0, 0.5, 0.25, 0.125] {
        let tree = BhTree::build(&net, &pos, opts_eps(eps, 8)).unwrap();
        let (e, g, _) = tree.energy_and_differential(&p, ex);
        let lists = tree.interaction_lists(&p);
        let (e2, g2, _) = tree
            .energy_and_differential_with_lists(&lists, &p, ex)
            .unwrap();
        assert_eq!(e.to_bits(), e2.to_bits());
        assert_eq!(bits(&g), bits(&g2));
        let (_, g_old, _) = tree.energy_and_differential_lumped_pairs(&p, ex);
        let (a, b) = (diff_err(&g, &g0), diff_err(&g_old, &g0));
        eprintln!("eps {eps:<6} |g_adj - g|/|g| {a:.2e}   |g_old - g|/|g| {b:.2e}");
        if eps <= 0.25 {
            assert!(a < 5e-3, "{a}");
        }
    }
}

// ---------------------------------------------------------------- determinism and errors

#[test]
fn bitwise_identical_at_1_3_and_8_threads() {
    let ex = Exponents::default();
    let curves = [
        {
            let (n, q) = trefoil(600);
            let q = perturb(&n, &q, 0.1, 12);
            (n, q)
        },
        {
            let (n, q) = network(60, 20);
            let q = perturb(&n, &q, 0.1, 13);
            (n, q)
        },
    ];
    for (net, pos) in &curves {
        let mut results = Vec::new();
        for t in [1, 3, 8] {
            let p = pool(t);
            let tree = BhTree::build(net, pos, opts_eps(0.5, 4)).unwrap();
            let (e, g, st) = tree.energy_and_differential(&p, ex);
            assert!(st.nodes_lumped > 0);
            let lists = tree.interaction_lists(&p);
            let moved = perturb(net, pos, 0.2, 14);
            let mut t2 = tree.clone();
            t2.refit(&moved).unwrap();
            let (e2, g2, _) = t2
                .energy_and_differential_with_lists(&lists, &p, ex)
                .unwrap();
            results.push((e.to_bits(), bits(&g), e2.to_bits(), bits(&g2)));
        }
        assert_eq!(results[0], results[1], "1 vs 3 threads");
        assert_eq!(results[0], results[2], "1 vs 8 threads");
    }
}

#[test]
fn mismatched_lists_are_errors_not_panics() {
    let p = pool(2);
    let ex = Exponents::default();
    let (net, pos) = trefoil(120);
    let tree = BhTree::build(&net, &pos, opts_eps(0.5, 4)).unwrap();
    let (net2, pos2) = trefoil(90);
    let other = BhTree::build(&net2, &pos2, opts_eps(0.5, 4)).unwrap();
    let bad = other.interaction_lists(&p);
    for lists in [bad, InteractionLists::default()] {
        let r = tree.energy_and_differential_with_lists(&lists, &p, ex);
        assert!(
            matches!(r, Err(RcError::InvalidOption { name: "lists", .. })),
            "{r:?}"
        );
    }
}

// ---------------------------------------------------------------- the flow (release, ignored)

fn data_dir() -> String {
    format!("{}/../../data/repulsive-curves", env!("CARGO_MANIFEST_DIR"))
}

/// Freedman–He–Wang 600, `(α, β) = (2, 4.5)`, Barnes–Hut energy (default options) with the
/// dense solver and otherwise the default flow options (8 threads): with the exact differential
/// of the Barnes–Hut energy the flow is a descent method on the function it evaluates and
/// reaches the stop rule `‖g̃‖ < 1e-4` (`Converged`) instead of ending `Stuck` at `‖g‖ ≈ 1e-3`
/// (§17f). The energy of the final curve is within 0.5 % of the exact flow's `6.0816`.
#[test]
#[ignore]
fn freedman_600_barnes_hut_flow_converges() {
    let path = format!("{}/FreedmanHeWang/freedman-he-wang-600.obj", data_dir());
    let (net, pos) = read_curve_obj(&path).unwrap();
    let ex = Exponents::new(2.0, 4.5).unwrap();
    let opts = FlowOptions {
        exponents: ex,
        max_iterations: 2000,
        record_path: false,
        threads: Some(8),
        backend: Backend::Explicit {
            energy: EnergyEval::BarnesHut(BhOptions::default()),
            solver: Solver::Dense,
        },
        ..FlowOptions::default()
    };
    let t0 = Instant::now();
    let mut s = Session::new(net, pos, opts).unwrap();
    let last_g = loop {
        let r = s.step().unwrap();
        if r.iteration.is_multiple_of(50) || s.status().is_terminal() {
            eprintln!(
                "step {:4} {:?} E {:.6} |g| {:.2e} tau {:.2e} halvings {} t {:.1}s",
                r.iteration,
                r.status,
                r.energy_after,
                r.gradient_norm,
                r.tau,
                r.halvings,
                t0.elapsed().as_secs_f64()
            );
        }
        if s.status().is_terminal() {
            break r.gradient_norm;
        }
    };
    let exact =
        ddg_repulsive_curves::energy(s.network(), s.positions(), ex, &EnergyEval::Exact, Some(8))
            .unwrap()
            .energy;
    eprintln!(
        "RESULT status {:?} steps {} |g| {last_g:.2e} Barnes-Hut energy {:.5} exact energy {exact:.5} \
         time {:.1}s",
        s.status(),
        s.iteration(),
        s.energy(),
        t0.elapsed().as_secs_f64()
    );
    assert!(matches!(s.status(), Status::Converged), "{:?}", s.status());
    assert!(last_g <= 1e-4, "{last_g}");
    assert!((exact - 6.0816).abs() < 0.005 * 6.0816, "exact {exact}");
    assert!(
        (s.energy() - 6.0816).abs() < 0.005 * 6.0816,
        "bh {}",
        s.energy()
    );
}

// ---------------------------------------------------------------- cost (release, ignored)

/// Time of the old and the new differential at Freedman 600/2400 and a 9472-vertex curve (the
/// 2368-edge Freedman curve subdivided twice), default options, 8 threads; minimum of 5 runs.
/// `cargo test --release -p ddg-repulsive-curves --test barnes_hut_adjoint -- --ignored --nocapture timing`.
#[test]
#[ignore]
fn timing_old_and_new_differential() {
    let p = pool(8);
    let ex = Exponents::new(2.0, 4.5).unwrap();
    let subdivide = |net: &CurveNetwork, pos: &[P]| {
        let mut q = pos.to_vec();
        let mut edges = Vec::new();
        for &[a, b] in net.edges() {
            let m = q.len() as u32;
            let (pa, pb) = (pos[a as usize], pos[b as usize]);
            q.push([
                0.5 * (pa[0] + pb[0]),
                0.5 * (pa[1] + pb[1]),
                0.5 * (pa[2] + pb[2]),
            ]);
            edges.push([a, m]);
            edges.push([m, b]);
        }
        (CurveNetwork::new(q.len(), edges).unwrap(), q)
    };
    let mut cases = Vec::new();
    for size in [600, 2400] {
        let path = format!("{}/FreedmanHeWang/freedman-he-wang-{size}.obj", data_dir());
        cases.push(read_curve_obj(&path).unwrap());
    }
    let (n, q) = subdivide(&cases[1].0, &cases[1].1);
    cases.push(subdivide(&n, &q));
    for (net, pos) in &cases {
        let tree = BhTree::build(net, pos, BhOptions::default()).unwrap();
        let time = |f: &dyn Fn()| {
            (0..5)
                .map(|_| {
                    let t = Instant::now();
                    f();
                    t.elapsed().as_secs_f64()
                })
                .fold(f64::INFINITY, f64::min)
        };
        let lists = tree.interaction_lists(&p);
        let t_old = time(&|| {
            tree.energy_and_differential_lumped_pairs(&p, ex);
        });
        let t_new = time(&|| {
            tree.energy_and_differential(&p, ex);
        });
        let t_new_lists = time(&|| {
            tree.energy_and_differential_with_lists(&lists, &p, ex)
                .unwrap();
        });
        let t_lists = time(&|| {
            tree.interaction_lists(&p);
        });
        let t_energy = time(&|| {
            tree.energy_with_lists(&lists, &p, ex);
        });
        eprintln!(
            "RESULT {} vertices: old {:.1} ms, adjoint {:.1} ms (with given lists {:.1} ms; lists {:.1} ms; \
             energy on lists {:.1} ms); lumped {} exact {}",
            pos.len(),
            1e3 * t_old,
            1e3 * t_new,
            1e3 * t_new_lists,
            1e3 * t_lists,
            1e3 * t_energy,
            lists.n_lumped(),
            lists.n_exact()
        );
    }
}
