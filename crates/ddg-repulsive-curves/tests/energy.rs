//! Tests of work package WP1 (`SPEC.md` §13.2–13.3, §2.3): finite differences of the energy
//! against the differential, invariances and scaling, validation errors, dataset loading and a
//! closed-form check on circles.

use ddg_repulsive_curves::energy::exact::{energy_and_differential, energy_only, make_pool};
use ddg_repulsive_curves::kernel::{kernel, kernel_with_gradient};
use ddg_repulsive_curves::scene::read_scene_curve;
use ddg_repulsive_curves::{
    Axis, CurveNetwork, EnergyEval, Exponents, RcError, check_planar, energy, parse_curve_obj,
    read_curve_obj, write_curve_obj,
};
use rayon::prelude::*;

type P = [f64; 3];

// ---------------------------------------------------------------- curve generators

fn closed_edges(n: usize) -> Vec<[u32; 2]> {
    (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect()
}

fn circle(n: usize, r: f64) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [r * t.cos(), r * t.sin(), 0.0]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n)).unwrap(), pos)
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

fn open_curve(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = 3.0 * k as f64 / n as f64;
            [t.cos(), t.sin(), 0.4 * t + 0.1 * (2.0 * t).sin()]
        })
        .collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

/// Three arms of `m` vertices each, joined at vertex 0 (a degree-3 junction), bent in 3-D.
fn tripod(m: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos = vec![[0.0, 0.0, 0.0]];
    let mut edges = Vec::new();
    let dirs = [[1.0, 0.0, 0.2], [-0.5, 0.9, -0.1], [-0.5, -0.8, 0.3]];
    for (a, d) in dirs.iter().enumerate() {
        let mut prev = 0u32;
        for k in 1..=m {
            let s = k as f64 / m as f64;
            let bend = 0.3 * s * s;
            let p = [
                d[0] * s + bend * (a as f64 + 1.0).cos(),
                d[1] * s + bend * (a as f64 + 1.0).sin(),
                d[2] * s + 0.2 * (3.0 * s + a as f64).sin(),
            ];
            pos.push(p);
            let id = pos.len() as u32 - 1;
            edges.push([prev, id]);
            prev = id;
        }
    }
    (CurveNetwork::new(pos.len(), edges).unwrap(), pos)
}

/// Two linked circles (Hopf link), two components.
fn hopf(n: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos = Vec::new();
    for k in 0..n {
        let t = std::f64::consts::TAU * k as f64 / n as f64;
        pos.push([t.cos(), t.sin(), 0.0]);
    }
    for k in 0..n {
        let t = std::f64::consts::TAU * k as f64 / n as f64;
        pos.push([1.0 + t.cos(), 0.0, t.sin()]);
    }
    let mut edges = Vec::new();
    for c in 0..2u32 {
        let o = c * n as u32;
        for k in 0..n as u32 {
            edges.push([o + k, o + (k + 1) % n as u32]);
        }
    }
    (CurveNetwork::new(2 * n, edges).unwrap(), pos)
}

/// A wiggly closed curve in the plane `z = 0.3`.
fn planar_curve(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            let r = 1.0 + 0.25 * (3.0 * t).cos();
            [1.3 * r * t.cos(), r * t.sin(), 0.3]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n)).unwrap(), pos)
}

/// Deterministic pseudo-random perturbation (LCG), relative to the mean edge length.
fn perturb(net: &CurveNetwork, pos: &[P], amount: f64, seed: u64) -> Vec<P> {
    let lbar = net.mean_edge_length(pos);
    let mut s = seed;
    let mut next = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    };
    pos.iter()
        .map(|p| {
            [
                p[0] + amount * lbar * next(),
                p[1] + amount * lbar * next(),
                p[2] + amount * lbar * next(),
            ]
        })
        .collect()
}

fn exponent_pairs() -> Vec<Exponents> {
    vec![
        Exponents::new(3.0, 6.0).unwrap(),
        Exponents::new(2.0, 4.5).unwrap(),
        Exponents::new(4.0, 8.0).unwrap(),
    ]
}

fn pool() -> rayon::ThreadPool {
    make_pool(Some(4)).unwrap()
}

fn all_curves() -> Vec<(&'static str, CurveNetwork, Vec<P>)> {
    let mut out = Vec::new();
    let (n, p) = circle(24, 1.3);
    out.push(("circle", n.clone(), perturb(&n, &p, 0.05, 1)));
    let (n, p) = trefoil(30);
    out.push(("trefoil", n.clone(), perturb(&n, &p, 0.05, 2)));
    let (n, p) = open_curve(20);
    out.push(("open", n.clone(), perturb(&n, &p, 0.05, 3)));
    let (n, p) = tripod(7);
    out.push(("junction", n.clone(), perturb(&n, &p, 0.05, 4)));
    let (n, p) = hopf(16);
    out.push(("two components", n.clone(), perturb(&n, &p, 0.05, 5)));
    let (n, p) = planar_curve(26);
    out.push(("planar", n, p));
    out
}

fn max_abs(g: &[P]) -> f64 {
    g.iter().flatten().fold(0.0, |m, x| m.max(x.abs()))
}

// ---------------------------------------------------------------- finite differences (13.2)

#[test]
fn differential_matches_central_differences() {
    let pool = pool();
    for e in exponent_pairs() {
        for (name, net, pos) in all_curves() {
            let (_, g) = energy_and_differential(&pool, &net, &pos, e);
            let h = 1e-6 * net.mean_edge_length(&pos);
            let mut fd = vec![[0.0; 3]; pos.len()];
            for v in 0..pos.len() {
                for c in 0..3 {
                    let mut p = pos.clone();
                    p[v][c] += h;
                    let ep = energy_only(&pool, &net, &p, e);
                    p[v][c] -= 2.0 * h;
                    let em = energy_only(&pool, &net, &p, e);
                    fd[v][c] = (ep - em) / (2.0 * h);
                }
            }
            let err = g
                .iter()
                .zip(&fd)
                .flat_map(|(a, b)| (0..3).map(move |c| (a[c] - b[c]).abs()))
                .fold(0.0, f64::max);
            let rel = err / max_abs(&g);
            assert!(
                rel < 1e-6,
                "{name} ({}, {}): relative FD error {rel:e}",
                e.alpha(),
                e.beta()
            );
        }
    }
}

#[test]
fn planar_curve_has_no_gradient_along_axis_and_passes_planar_check() {
    let (net, pos) = planar_curve(26);
    check_planar(&pos, Axis::Z).unwrap();
    let r = energy(
        &net,
        &pos,
        Exponents::default(),
        &EnergyEval::Exact,
        Some(2),
    )
    .unwrap();
    let scale = max_abs(&r.differential);
    assert!(r.differential.iter().all(|g| g[2].abs() < 1e-14 * scale));
    assert!(matches!(
        check_planar(&pos, Axis::X),
        Err(RcError::NotPlanar {
            vertex: 1,
            axis: 0,
            ..
        })
    ));
}

#[test]
fn kernel_gradient_is_continuous_at_nu_zero() {
    for e in [
        Exponents::new(3.0, 6.0).unwrap(),
        Exponents::new(1.5, 3.6).unwrap(),
        Exponents::new(2.0, 4.5).unwrap(),
    ] {
        let t = [1.0, 0.0, 0.0];
        let q = [0.0; 3];
        // p exactly on the tangent line: nu = 0.
        let on = kernel_with_gradient([0.7, 0.0, 0.0], q, t, e);
        assert!(
            on.value == 0.0
                && on.grad_p.iter().all(|x| x.is_finite())
                && on.grad_t.iter().all(|x| x.is_finite())
        );
        // one-sided limits agree with the value at nu = 0.
        for eps in [1e-3f64, 1e-6, 1e-9] {
            let bound = 50.0 * eps.powf(e.alpha() - 1.0) + 50.0 * eps * eps;
            for p in [[0.7, eps, 0.0], [0.7, 0.0, -eps], [0.7, -eps, eps]] {
                let near = kernel_with_gradient(p, q, t, e);
                for c in 0..3 {
                    assert!(
                        (near.grad_p[c] - on.grad_p[c]).abs() < bound,
                        "alpha {} eps {eps}",
                        e.alpha()
                    );
                }
            }
        }
        // Tiny nu must not give NaN even for alpha close to 1.
        let tiny = kernel_with_gradient([0.7, 1e-300, 0.0], q, t, e);
        assert!(tiny.grad_p.iter().all(|x| x.is_finite()));
    }
}

#[test]
fn kernel_gradient_matches_finite_differences_pointwise() {
    let e = Exponents::new(2.0, 4.5).unwrap();
    let t = {
        let v: [f64; 3] = [0.3, -0.8, 0.5];
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / n, v[1] / n, v[2] / n]
    };
    let (p, q) = ([0.4, 0.9, -0.2], [-0.3, 0.1, 0.35]);
    let ke = kernel_with_gradient(p, q, t, e);
    assert_eq!(ke.value, kernel(p, q, t, e));
    let h = 1e-6;
    for c in 0..3 {
        let (mut a, mut b) = (p, p);
        a[c] += h;
        b[c] -= h;
        let fd = (kernel(a, q, t, e) - kernel(b, q, t, e)) / (2.0 * h);
        assert!((fd - ke.grad_p[c]).abs() < 1e-7 * (1.0 + fd.abs()));
    }
}

// ---------------------------------------------------------------- invariances (13.3)

fn rot(p: P, axis: P, angle: f64) -> P {
    let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let k = [axis[0] / n, axis[1] / n, axis[2] / n];
    let (s, c) = angle.sin_cos();
    let kxp = [
        k[1] * p[2] - k[2] * p[1],
        k[2] * p[0] - k[0] * p[2],
        k[0] * p[1] - k[1] * p[0],
    ];
    let kp = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
    [0, 1, 2].map(|i| p[i] * c + kxp[i] * s + k[i] * kp * (1.0 - c))
}

#[test]
fn gradient_identities_and_euler() {
    let pool = pool();
    for e in exponent_pairs() {
        for (name, net, pos) in all_curves() {
            let (en, g) = energy_and_differential(&pool, &net, &pos, e);
            let scale = max_abs(&g) * pos.len() as f64;
            let mut sum = [0.0; 3];
            let mut torque = [0.0; 3];
            let mut euler = 0.0;
            for (p, gv) in pos.iter().zip(&g) {
                for c in 0..3 {
                    sum[c] += gv[c];
                }
                torque[0] += p[1] * gv[2] - p[2] * gv[1];
                torque[1] += p[2] * gv[0] - p[0] * gv[2];
                torque[2] += p[0] * gv[1] - p[1] * gv[0];
                euler += p[0] * gv[0] + p[1] * gv[1] + p[2] * gv[2];
            }
            let tol = 1e-12 * scale;
            assert!(sum.iter().all(|x| x.abs() < tol), "{name}: sum {sum:?}");
            assert!(
                torque.iter().all(|x| x.abs() < tol * 3.0),
                "{name}: torque {torque:?}"
            );
            let target = e.energy_scaling() * en;
            assert!(
                (euler - target).abs() < 1e-12 * en.abs().max(scale),
                "{name}: euler {euler} vs {target}"
            );
        }
    }
}

#[test]
fn rigid_motions_reflection_reorientation_renumbering() {
    let pool = pool();
    let e = Exponents::default();
    for (name, net, pos) in all_curves() {
        let (e0, g0) = energy_and_differential(&pool, &net, &pos, e);
        let tol = 1e-12;
        let gmax = max_abs(&g0);
        // translation
        let tr: Vec<P> = pos
            .iter()
            .map(|p| [p[0] + 3.0, p[1] - 2.0, p[2] + 0.5])
            .collect();
        let (e1, g1) = energy_and_differential(&pool, &net, &tr, e);
        assert!((e1 - e0).abs() < tol * e0, "{name} translation");
        for (a, b) in g0.iter().zip(&g1) {
            for c in 0..3 {
                assert!((a[c] - b[c]).abs() < tol * gmax);
            }
        }
        // rotation: energy invariant, gradient rotates
        let (ax, an) = ([0.3, -1.0, 0.7], 1.1);
        let rp: Vec<P> = pos.iter().map(|&p| rot(p, ax, an)).collect();
        let (e2, g2) = energy_and_differential(&pool, &net, &rp, e);
        assert!((e2 - e0).abs() < tol * e0, "{name} rotation");
        for (a, b) in g0.iter().zip(&g2) {
            let ra = rot(*a, ax, an);
            for c in 0..3 {
                assert!((ra[c] - b[c]).abs() < tol * gmax);
            }
        }
        // reflection
        let rf: Vec<P> = pos.iter().map(|p| [-p[0], p[1], p[2]]).collect();
        let (e3, _) = energy_and_differential(&pool, &net, &rf, e);
        assert!((e3 - e0).abs() < tol * e0, "{name} reflection");
        // edge re-orientation: flip every other edge
        let edges: Vec<[u32; 2]> = net
            .edges()
            .iter()
            .enumerate()
            .map(|(i, &[a, b])| if i % 2 == 0 { [b, a] } else { [a, b] })
            .collect();
        let net2 = CurveNetwork::new(net.n_vertices(), edges).unwrap();
        let (e4, g4) = energy_and_differential(&pool, &net2, &pos, e);
        assert!((e4 - e0).abs() < tol * e0, "{name} reorientation");
        for (a, b) in g0.iter().zip(&g4) {
            for c in 0..3 {
                assert!((a[c] - b[c]).abs() < tol * gmax);
            }
        }
        // renumbering: reverse the vertex order
        let n = pos.len();
        let pos5: Vec<P> = (0..n).map(|v| pos[n - 1 - v]).collect();
        let edges5: Vec<[u32; 2]> = net
            .edges()
            .iter()
            .map(|&[a, b]| [n as u32 - 1 - a, n as u32 - 1 - b])
            .collect();
        let net5 = CurveNetwork::new(n, edges5).unwrap();
        let (e5, g5) = energy_and_differential(&pool, &net5, &pos5, e);
        assert!((e5 - e0).abs() < tol * e0, "{name} renumbering");
        for v in 0..n {
            for c in 0..3 {
                assert!((g0[v][c] - g5[n - 1 - v][c]).abs() < tol * gmax);
            }
        }
    }
}

#[test]
fn energy_scales_with_the_exponent_2_plus_alpha_minus_beta() {
    let pool = pool();
    for e in exponent_pairs() {
        for (name, net, pos) in all_curves() {
            let (e0, g0) = energy_and_differential(&pool, &net, &pos, e);
            for c in [1e-3, 0.7, 1e3] {
                let sp: Vec<P> = pos.iter().map(|p| [c * p[0], c * p[1], c * p[2]]).collect();
                let (e1, g1) = energy_and_differential(&pool, &net, &sp, e);
                let want = c.powf(e.energy_scaling()) * e0;
                assert!(
                    ((e1 - want) / want).abs() < 1e-12,
                    "{name} c={c} ({}, {}): {e1} vs {want}",
                    e.alpha(),
                    e.beta()
                );
                // the differential scales with c^{1+α−β}
                let gs = c.powf(1.0 + e.alpha() - e.beta());
                for (a, b) in g0.iter().zip(&g1) {
                    for k in 0..3 {
                        assert!(
                            (b[k] - gs * a[k]).abs() < 1e-11 * gs * max_abs(&g0),
                            "{name} c={c}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn result_is_bitwise_independent_of_the_thread_count() {
    let (net, pos) = trefoil(150);
    let pos = perturb(&net, &pos, 0.03, 9);
    let e = Exponents::default();
    let r1 = energy(&net, &pos, e, &EnergyEval::Exact, Some(1)).unwrap();
    for t in [2, 3, 7, 16] {
        let r = energy(&net, &pos, e, &EnergyEval::Exact, Some(t)).unwrap();
        assert_eq!(r1.energy.to_bits(), r.energy.to_bits(), "threads {t}");
        for (a, b) in r1.differential.iter().zip(&r.differential) {
            for c in 0..3 {
                assert_eq!(a[c].to_bits(), b[c].to_bits(), "threads {t}");
            }
        }
    }
    assert_eq!(r1.pairs_exact, 150 * 149 - 150 * 2);
    assert!(matches!(
        energy(&net, &pos, e, &EnergyEval::Exact, Some(0)),
        Err(RcError::InvalidOption {
            name: "threads",
            ..
        })
    ));
    // The Barnes-Hut energy and differential are deterministic across thread counts too.
    let bh = EnergyEval::BarnesHut(Default::default());
    let (net, pos) = trefoil(600);
    let pos = perturb(&net, &pos, 0.03, 9);
    let b1 = energy(&net, &pos, e, &bh, Some(1)).unwrap();
    assert!(b1.nodes_lumped > 0, "the test must exercise the far field");
    for t in [2, 3, 7, 16] {
        let r = energy(&net, &pos, e, &bh, Some(t)).unwrap();
        assert_eq!(b1.energy.to_bits(), r.energy.to_bits(), "bh threads {t}");
        assert_eq!(b1.pairs_exact, r.pairs_exact, "bh threads {t}");
        assert_eq!(b1.nodes_lumped, r.nodes_lumped, "bh threads {t}");
        for (a, b) in b1.differential.iter().zip(&r.differential) {
            for c in 0..3 {
                assert_eq!(a[c].to_bits(), b[c].to_bits(), "bh threads {t}");
            }
        }
    }
    assert!(matches!(
        energy(&net, &pos, e, &bh, Some(0)),
        Err(RcError::InvalidOption {
            name: "threads",
            ..
        })
    ));
}

// ---------------------------------------------------------------- closed forms

/// For a round circle of radius r, |T×d| = |d|²/(2r), so the continuous energy is
/// 2π r² (2r)^{α−β} ∫_0^{2π} sin(θ/2)^{2α−β} dθ. For (3,6) this is π²/(2r).
#[test]
fn circle_energy_converges_to_closed_form() {
    let pool = pool();
    let e = Exponents::new(3.0, 6.0).unwrap();
    for r in [1.0, 2.5] {
        let exact = std::f64::consts::PI.powi(2) / (2.0 * r);
        let mut prev = f64::INFINITY;
        for n in [50usize, 100, 200, 400, 800] {
            let (net, pos) = circle(n, r);
            let en = energy_only(&pool, &net, &pos, e);
            let rel = (en - exact).abs() / exact;
            // Omitted pairs: 3 of n edges per row, so the deficit is about 3/n (Sec. 5.1).
            assert!(rel < 5.0 / n as f64, "n={n} r={r}: rel {rel:e}");
            assert!(
                rel < 0.6 * prev || n == 50,
                "no convergence: {rel:e} after {prev:e}"
            );
            prev = rel;
        }
    }
    // (2, 4.5): ∫_0^{2π} sin(θ/2)^{-1/2} dθ = 2 √π Γ(1/4)/Γ(3/4); singular, slower convergence.
    let e = Exponents::new(2.0, 4.5).unwrap();
    let b = 5.244115108584239; // √π Γ(1/4) / Γ(3/4)
    let exact = 2.0 * std::f64::consts::PI * 2f64.powf(-2.5) * 2.0 * b;
    let mut prev = f64::INFINITY;
    for n in [100usize, 200, 400, 800, 1600] {
        let (net, pos) = circle(n, 1.0);
        let en = energy_only(&pool, &net, &pos, e);
        let rel = (en - exact).abs() / exact;
        assert!(rel < prev, "no convergence: {rel:e} after {prev:e}");
        prev = rel;
    }
    assert!(prev < 0.02, "final relative error {prev:e}");
}

// ---------------------------------------------------------------- validation (2.3)

#[test]
fn topology_errors() {
    let pe = |e: Vec<[u32; 2]>, n| CurveNetwork::new(n, e).unwrap_err();
    let err = pe(vec![[0, 1], [1, 5], [2, 3]], 4);
    assert_eq!(
        err,
        RcError::VertexOutOfRange {
            edge: 1,
            vertex: 5,
            n_vertices: 4
        }
    );
    assert!(err.to_string().contains("edge 1"));
    let err = pe(vec![[0, 1], [3, 3], [2, 3]], 4);
    assert_eq!(err, RcError::SelfLoop { edge: 1, vertex: 3 });
    assert!(err.to_string().contains("edge 1 (3, 3) is degenerate"));
    let err = pe(vec![[0, 1], [2, 3], [1, 0]], 4);
    assert_eq!(
        err,
        RcError::DuplicateEdge {
            a: 0,
            b: 1,
            edges: [0, 2]
        }
    );
    let err = pe(vec![[0, 1], [2, 3]], 5);
    assert_eq!(err, RcError::IsolatedVertex { vertex: 4 });
    // single edge, triangle, star: P empty
    for (edges, n) in [
        (vec![[0, 1]], 2usize),
        (vec![[0, 1], [1, 2], [2, 0]], 3),
        (vec![[0, 1], [0, 2], [0, 3]], 4),
    ] {
        assert!(matches!(
            CurveNetwork::new(n, edges),
            Err(RcError::NoInteractingPairs { .. })
        ));
    }
    // a square is fine
    CurveNetwork::new(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]]).unwrap();
}

fn two_edges() -> (CurveNetwork, Vec<P>) {
    let net = CurveNetwork::new(4, vec![[0, 1], [2, 3]]).unwrap();
    (
        net,
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
    )
}

#[test]
fn geometry_errors() {
    let (net, pos) = two_edges();
    net.validate_positions(&pos).unwrap();
    assert_eq!(
        net.validate_positions(&pos[..3]).unwrap_err(),
        RcError::PositionCountMismatch {
            expected: 4,
            found: 3
        }
    );
    let mut p = pos.clone();
    p[2][1] = f64::NAN;
    assert_eq!(
        net.validate_positions(&p).unwrap_err(),
        RcError::NonFinitePosition { vertex: 2 }
    );
    p[2][1] = f64::INFINITY;
    assert_eq!(
        net.validate_positions(&p).unwrap_err(),
        RcError::NonFinitePosition { vertex: 2 }
    );
    let mut p = pos.clone();
    p[3] = p[2];
    let err = net.validate_positions(&p).unwrap_err();
    assert_eq!(
        err,
        RcError::ZeroLengthEdge {
            edge: 1,
            a: 2,
            b: 3
        }
    );
    assert!(err.to_string().contains("edge 1"));
    // coincident vertices of different edges
    let mut p = pos.clone();
    p[2] = p[1];
    p[3] = [1.0, 1.0, 0.0];
    let err = net.validate_positions(&p).unwrap_err();
    assert_eq!(err, RcError::CoincidentVertices { a: 1, b: 2 });
    assert!(err.to_string().contains("vertices 1 and 2 coincide"));
    // crossing edges
    let p = vec![
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    match net.validate_positions(&p).unwrap_err() {
        RcError::EdgesIntersect {
            a: 0,
            b: 1,
            distance,
        } => assert_eq!(distance, 0.0),
        other => panic!("{other:?}"),
    }
    // touching T-junction without a shared vertex
    let p = vec![
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    assert!(matches!(
        net.validate_positions(&p),
        Err(RcError::EdgesIntersect { a: 0, b: 1, .. })
    ));
    // close but separate is fine
    let p = vec![
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1e-6, 0.0],
        [0.0, 1.0, 0.0],
    ];
    net.validate_positions(&p).unwrap();
    // exponents
    for (a, b) in [
        (1.0, 3.0),
        (0.5, 3.0),
        (3.0, 4.9),
        (3.0, 7.0),
        (3.0, 7.5),
        (f64::NAN, 6.0),
        (3.0, f64::INFINITY),
    ] {
        assert!(
            matches!(Exponents::new(a, b), Err(RcError::InvalidExponents { .. })),
            "({a}, {b})"
        );
    }
    Exponents::new(3.0, 5.0).unwrap();
    Exponents::new(3.0, 6.999).unwrap();
    assert_eq!(Exponents::default(), Exponents::new(3.0, 6.0).unwrap());
    // the energy entry point validates too
    assert!(matches!(
        energy(
            &net,
            &p[..2],
            Exponents::default(),
            &EnergyEval::Exact,
            Some(1)
        ),
        Err(RcError::PositionCountMismatch { .. })
    ));
}

#[test]
fn components_and_degrees() {
    let (net, _) = hopf(8);
    assert_eq!(net.components().len(), 2);
    assert!(
        net.components()
            .iter()
            .all(|c| c.closed && c.vertices.len() == 8 && c.edges.len() == 8)
    );
    assert_eq!(net.component_of(9), 1);
    let (net, pos) = tripod(5);
    assert_eq!(net.degree(0), 3);
    assert_eq!(net.components().len(), 1);
    assert!(!net.components()[0].closed);
    let m = net.dual_lengths(&pos);
    let total: f64 = net.edge_geometry(&pos).length.iter().sum();
    assert!((m.iter().sum::<f64>() - total).abs() < 1e-13);
}

// ---------------------------------------------------------------- OBJ I/O and datasets

#[test]
fn obj_parser() {
    let txt = "# c\nv 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nl 1 2 3 4 1\n";
    let (net, pos) = parse_curve_obj(txt).unwrap();
    assert_eq!((net.n_vertices(), net.n_edges(), pos.len()), (4, 4, 4));
    assert!(net.components()[0].closed);
    // two-index lines, slash forms, negative indices
    let txt = "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nl 1/1 2/2\nl 2 3\nl -2 -1\nl 4 1\n";
    let (net, _) = parse_curve_obj(txt).unwrap();
    assert_eq!(net.edges(), &[[0, 1], [1, 2], [2, 3], [3, 0]]);
    // errors with line numbers
    let e = parse_curve_obj("v 0 0\n").unwrap_err();
    assert!(matches!(e, RcError::Parse { line: 1, .. }), "{e}");
    let e = parse_curve_obj("v 0 0 0\nv 1 0 0\nl 1 7\n").unwrap_err();
    assert!(matches!(e, RcError::Parse { line: 3, .. }), "{e}");
    let e = parse_curve_obj("v 0 0 0\nv 1 0 0\nl 1\n").unwrap_err();
    assert!(matches!(e, RcError::Parse { line: 3, .. }), "{e}");
    let e = parse_curve_obj("v 0 0 x\n").unwrap_err();
    assert!(matches!(e, RcError::Parse { line: 1, .. }), "{e}");
    let e = parse_curve_obj("v 0 0 0\nv 1 0 0\nl 1 0\n").unwrap_err();
    assert!(matches!(e, RcError::Parse { line: 3, .. }), "{e}");
    // vertices without lines: isolated
    assert!(matches!(
        parse_curve_obj("v 0 0 0\nv 1 0 0\n"),
        Err(RcError::IsolatedVertex { vertex: 0 })
    ));
    assert!(matches!(
        read_curve_obj("/nonexistent/x.obj"),
        Err(RcError::Io { .. })
    ));
}

#[test]
fn obj_roundtrip() {
    let (net, pos) = trefoil(40);
    let pos = perturb(&net, &pos, 0.01, 3);
    let dir = std::env::temp_dir().join(format!("rc-wp1-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("t.obj");
    write_curve_obj(&path, &net, &pos).unwrap();
    let (net2, pos2) = read_curve_obj(&path).unwrap();
    assert_eq!(net, net2);
    assert_eq!(pos, pos2);
    std::fs::remove_dir_all(&dir).ok();
    assert!(matches!(
        write_curve_obj("/nonexistent/dir/x.obj", &net, &pos),
        Err(RcError::Io { .. })
    ));
    assert!(matches!(
        write_curve_obj("/nonexistent/x.obj", &net, &pos[..3]),
        Err(RcError::PositionCountMismatch { .. })
    ));
}

fn collect_objs(dir: &std::path::Path, out: &mut Vec<std::path::PathBuf>) {
    let mut entries: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .collect();
    entries.sort();
    for p in entries {
        if p.is_dir() {
            collect_objs(&p, out);
        } else if p.extension().is_some_and(|x| x == "obj") {
            out.push(p);
        }
    }
}

/// Every dataset file loads and has a finite positive energy. The three Freedman-He-Wang files
/// (up to 2400 edges) dominate the time; the whole test is parallel over files.
#[test]
fn all_dataset_files_load_with_finite_energy() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/repulsive-curves");
    if !root.exists() {
        eprintln!("data/repulsive-curves not present; skipping");
        return;
    }
    let mut files = Vec::new();
    collect_objs(&root, &mut files);
    assert!(files.len() >= 229, "only {} files found", files.len());
    // Surface and obstacle meshes (several faces, no `l` line) are not curves.
    let is_surface = |f: &std::path::PathBuf| {
        let text = std::fs::read_to_string(f).unwrap_or_default();
        let tag = |t: &str| {
            text.lines()
                .filter(|l| l.split_whitespace().next() == Some(t))
                .count()
        };
        tag("l") == 0 && tag("f") > 1
    };
    let (surfaces, files): (Vec<_>, Vec<_>) = files.into_iter().partition(is_surface);
    assert!(
        surfaces.len() <= 13 && files.len() >= 216,
        "{} surface meshes, {} curve files",
        surfaces.len(),
        files.len()
    );
    let outer = make_pool(None).unwrap();
    let inner = make_pool(Some(1)).unwrap();
    // Small files in parallel (one thread each); the big ones afterwards with the full pool.
    let (small, big): (Vec<_>, Vec<_>) = files
        .iter()
        .partition(|f| std::fs::metadata(f).map_or(0, |m| m.len()) < 40_000);
    let run = |f: &std::path::PathBuf, pool: &rayon::ThreadPool| -> Result<(String, f64), String> {
        let name = f.display().to_string();
        let (net, pos) = read_scene_curve(f).map_err(|e| format!("{name}: {e}"))?;
        let (en, g) = energy_and_differential(pool, &net, &pos, Exponents::default());
        if !(en.is_finite() && en > 0.0 && g.iter().flatten().all(|x| x.is_finite())) {
            return Err(format!("{name}: energy {en}"));
        }
        Ok((name, en))
    };
    let mut results: Vec<Result<(String, f64), String>> =
        outer.install(|| small.par_iter().map(|f| run(f, &inner)).collect());
    results.extend(big.iter().map(|f| run(f, &outer)));
    let errors: Vec<_> = results.iter().filter_map(|r| r.as_ref().err()).collect();
    assert!(
        errors.is_empty(),
        "{} files failed: {:?}",
        errors.len(),
        &errors[..errors.len().min(5)]
    );
    eprintln!(
        "{} curve files loaded, {} surface meshes skipped",
        files.len(),
        surfaces.len()
    );
}

/// The authors' two single-polygon curve files are closed loops through the face's vertices.
#[test]
fn scene_curves_given_as_one_polygon_face_are_closed_loops() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/repulsive-curves");
    for rel in [
        "scenes/DoubleTorus/InitialLoop.obj",
        "scenes/FreedmanUnknot/FreedmanKnotPlot.obj",
    ] {
        let path = root.join(rel);
        if !path.exists() {
            eprintln!("{rel} not present; skipping");
            continue;
        }
        let (net, pos) = read_scene_curve(&path).unwrap();
        assert_eq!(net.n_edges(), net.n_vertices(), "{rel}");
        assert!(
            (0..net.n_vertices() as u32).all(|v| net.degree(v) == 2),
            "{rel}"
        );
        assert_eq!(net.components().len(), 1, "{rel}");
        assert!(
            read_curve_obj(&path).is_err(),
            "{rel}: the plain reader rejects the file"
        );
        assert!(pos.iter().flatten().all(|x| x.is_finite()));
    }
    // A face with a repeated first vertex and `v/vt/vn` indices still reads.
    let text = "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1/1/1 2/2/2 3/3/3 4/4/4\n";
    let (net, _) = ddg_repulsive_curves::scene::parse_scene_curve(text).unwrap();
    assert_eq!((net.n_vertices(), net.n_edges()), (4, 4));
    // Two faces: not a curve.
    let two = format!("{text}f 1 2 3\n");
    assert!(ddg_repulsive_curves::scene::parse_scene_curve(&two).is_err());
}

/// Timing of energy + differential on Freedman-He-Wang 2400 (run alone, in release:
/// `cargo test -p ddg-repulsive-curves --release --test energy -- --ignored --nocapture timing`).
#[test]
#[ignore]
fn timing_freedman_2400() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/repulsive-curves/FreedmanHeWang/freedman-he-wang-2400.obj");
    let (net, pos) = read_curve_obj(path).unwrap();
    for threads in [1usize, 64] {
        let pool = make_pool(Some(threads)).unwrap();
        let e = Exponents::default();
        let _ = energy_and_differential(&pool, &net, &pos, e);
        let t = std::time::Instant::now();
        let reps = 5;
        for _ in 0..reps {
            let _ = energy_and_differential(&pool, &net, &pos, e);
        }
        eprintln!(
            "freedman-2400, {threads} threads: {:.1} ms per energy+differential",
            t.elapsed().as_secs_f64() * 1e3 / reps as f64
        );
    }
}
