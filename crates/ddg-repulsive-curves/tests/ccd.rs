//! Continuous collision detection (`SPEC.md` §7.4, §13.5): constructed cases with known first
//! contact times, randomized pairs against a sampled-and-bisected ground truth, networks,
//! obstacles, determinism and bad input.

use std::sync::Arc;

use ddg_repulsive_curves::CurveNetwork;
use ddg_repulsive_curves::RcError;
use ddg_repulsive_curves::ccd::{
    CcdOptions, CollisionObstacle, ContactPrimitive, PreparedObstacle, first_contact,
    min_pair_distance, point_triangle_advance, point_triangle_distance, segment_segment_advance,
    segment_segment_distance,
};

type P = [f64; 3];

fn pool(n: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .unwrap()
}

fn at(p: P, u: P, t: f64) -> P {
    [p[0] + t * u[0], p[1] + t * u[1], p[2] + t * u[2]]
}

/// A small deterministic generator (splitmix64).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    }
    fn sym(&mut self, s: f64) -> f64 {
        s * (2.0 * self.next() - 1.0)
    }
    fn vec(&mut self, s: f64) -> P {
        [self.sym(s), self.sym(s), self.sym(s)]
    }
}

fn pair_dist(p: [P; 2], u: [P; 2], q: [P; 2], w: [P; 2], t: f64) -> f64 {
    segment_segment_distance(
        at(p[0], u[0], t),
        at(p[1], u[1], t),
        at(q[0], w[0], t),
        at(q[1], w[1], t),
    )
}

/// First time in `[0, 1]` where `f(t) ≤ 0`, by sampling `n` points and bisecting to 1e-14.
fn first_root(f: impl Fn(f64) -> f64, n: usize) -> Option<f64> {
    if f(0.0) <= 0.0 {
        return Some(0.0);
    }
    let mut prev = 0.0;
    for k in 1..=n {
        let t = k as f64 / n as f64;
        if f(t) <= 0.0 {
            let (mut lo, mut hi) = (prev, t);
            while hi - lo > 1e-14 {
                let m = 0.5 * (lo + hi);
                if f(m) <= 0.0 {
                    hi = m;
                } else {
                    lo = m;
                }
            }
            return Some(lo);
        }
        prev = t;
    }
    None
}

const OPT: CcdOptions = CcdOptions {
    gap: 0.1,
    max_advances: 100,
};

fn ss(p: [P; 2], u: [P; 2], q: [P; 2], w: [P; 2]) -> Option<f64> {
    let d0 = segment_segment_distance(p[0], p[1], q[0], q[1]);
    segment_segment_advance(p, u, q, w, d0, 0.0, 10.0, &OPT).0
}

#[test]
fn distances_of_constructed_configurations() {
    // Parallel overlapping, parallel disjoint, crossing, skew, touching endpoints.
    let d = segment_segment_distance([0., 0., 0.], [1., 0., 0.], [0.5, 1., 0.], [2., 1., 0.]);
    assert!((d - 1.0).abs() < 1e-15);
    let d = segment_segment_distance([0., 0., 0.], [1., 0., 0.], [3., 0., 0.], [4., 0., 0.]);
    assert!((d - 2.0).abs() < 1e-15);
    let d = segment_segment_distance([-1., 0., 0.], [1., 0., 0.], [0., -1., 0.5], [0., 1., 0.5]);
    assert!((d - 0.5).abs() < 1e-15);
    let d = segment_segment_distance([-1., 0., 0.], [1., 0., 0.], [0., -1., 0.], [0., 1., 0.]);
    assert_eq!(d, 0.0);
    let d = segment_segment_distance([0., 0., 0.], [1., 0., 0.], [1., 0., 0.], [1., 1., 1.]);
    assert_eq!(d, 0.0);
    // Skew, closest points at the ends of the second segment.
    let d = segment_segment_distance([-1., 0., 0.], [1., 0., 0.], [0., 2., 1.], [0., 3., 1.]);
    assert!((d - 5f64.sqrt()).abs() < 1e-15);
    // Point–triangle: above the face, beyond an edge, beyond a corner, degenerate triangle.
    let t = [[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]];
    assert!((point_triangle_distance([0.5, 0.5, 0.7], t[0], t[1], t[2]) - 0.7).abs() < 1e-15);
    assert!((point_triangle_distance([1., -1., 0.], t[0], t[1], t[2]) - 1.0).abs() < 1e-15);
    assert!(
        (point_triangle_distance([-1., -1., 0.], t[0], t[1], t[2]) - 2f64.sqrt()).abs() < 1e-15
    );
    let dt = point_triangle_distance([1., 1., 0.], [0., 0., 0.], [1., 0., 0.], [2., 0., 0.]);
    assert!((dt - 1.0).abs() < 1e-15);
}

#[test]
fn exact_first_contact_times() {
    let z = [0.0; 3];
    // Head-on, parallel: d = 1 − τ, gap 0.1 → 0.9.
    let t = ss(
        [[0., 0., 0.], [1., 0., 0.]],
        [z, z],
        [[0., 1., 0.], [1., 1., 0.]],
        [[0., -1., 0.], [0., -1., 0.]],
    )
    .unwrap();
    assert!((t - 0.9).abs() < 1e-12, "{t}");
    // Crossing (perpendicular), the upper one falls at speed 2: d = 1 − 2τ → 0.45.
    let t = ss(
        [[-1., 0., 0.], [1., 0., 0.]],
        [z, z],
        [[0., -1., 1.], [0., 1., 1.]],
        [[0., 0., -2.], [0., 0., -2.]],
    )
    .unwrap();
    assert!((t - 0.45).abs() < 1e-12, "{t}");
    // Both move, towards each other: d = 2 − 2τ → 0.9.
    let t = ss(
        [[-1., 0., 0.], [1., 0., 0.]],
        [[0., 0., 1.], [0., 0., 1.]],
        [[0., -1., 2.], [0., 1., 2.]],
        [[0., 0., -1.], [0., 0., -1.]],
    )
    .unwrap();
    assert!((t - 0.9).abs() < 1e-12, "{t}");
    // Skew, oblique approach (several advances): d = 1 − τ → 0.9, never later.
    let t = ss(
        [[-5., 0., 0.], [5., 0., 0.]],
        [z, z],
        [[0., -1., 1.], [0., 1., 1.]],
        [[1., 0., -1.], [1., 0., -1.]],
    )
    .unwrap();
    assert!(t <= 0.9 && t > 0.9 - 1e-8, "{t}");
    // Touching at endpoints, collinear: d = 1 − τ → 0.9.
    let t = ss(
        [[0., 0., 0.], [1., 0., 0.]],
        [z, z],
        [[2., 0., 0.], [3., 0., 0.]],
        [[-1., 0., 0.], [-1., 0., 0.]],
    )
    .unwrap();
    assert!((t - 0.9).abs() < 1e-12, "{t}");
    // Planar: an endpoint hits the interior of the other segment in the plane z = 0.
    let t = ss(
        [[0., 0., 0.], [1., 0., 0.]],
        [z, z],
        [[0.5, 1., 0.], [0.5, 2., 0.]],
        [[0., -1., 0.], [0., -1., 0.]],
    )
    .unwrap();
    assert!((t - 0.9).abs() < 1e-12, "{t}");
    // Rigid translation of the pair together: never.
    assert!(
        ss(
            [[0., 0., 0.], [1., 0., 0.]],
            [[1., 2., 3.], [1., 2., 3.]],
            [[0., 1., 0.], [1., 1., 0.]],
            [[1., 2., 3.], [1., 2., 3.]],
        )
        .is_none()
    );
    // Parallel sliding: the distance stays 1.
    assert!(
        ss(
            [[0., 0., 0.], [1., 0., 0.]],
            [z, z],
            [[0., 1., 0.], [1., 1., 0.]],
            [[5., 0., 0.], [5., 0., 0.]],
        )
        .is_none()
    );
    // Rotating segment (endpoints with different velocities): against the bisected truth.
    let (p, u) = ([[0., 0., 0.], [4., 0., 0.]], [z, z]);
    let (q, w) = (
        [[1., 0.5, 1.], [1., 1.5, 3.]],
        [[0.5, -0.5, -1.2], [0.0, -1.0, -2.0]],
    );
    let d0 = segment_segment_distance(p[0], p[1], q[0], q[1]);
    let truth = first_root(|t| pair_dist(p, u, q, w, t) - 0.1 * d0, 100_000).unwrap();
    let t = segment_segment_advance(p, u, q, w, d0, 0.0, 1.0, &OPT)
        .0
        .unwrap();
    assert!(t <= truth + 1e-12 && t > truth - 1e-6, "{t} vs {truth}");
    // Point into a triangle face: d = 1 − 2τ → 0.45.
    let tri = [[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]];
    let (t, _) = point_triangle_advance([0.5, 0.5, 1.], [0., 0., -2.], tri, 1.0, 0.0, 1.0, &OPT);
    assert!((t.unwrap() - 0.45).abs() < 1e-12);
}

/// Random moving segment pairs: the CCD time is never after the first time the distance
/// reaches `s d₀`, and the distance stays above `s d₀` on `[0, τ_ccd]`.
fn random_pairs(n_pairs: usize, samples: usize, seed: u64) {
    let mut rng = Rng(seed);
    let (mut contacts, mut free, mut capped) = (0, 0, 0);
    for k in 0..n_pairs {
        let c = rng.vec(1.0);
        let p = [rng.vec(1.0), rng.vec(1.0)];
        let q = [
            [
                c[0] + rng.sym(1.0),
                c[1] + rng.sym(1.0),
                c[2] + rng.sym(1.0),
            ],
            [
                c[0] + rng.sym(1.0),
                c[1] + rng.sym(1.0),
                c[2] + rng.sym(1.0),
            ],
        ];
        let speed = if k % 3 == 0 { 0.3 } else { 3.0 };
        let u = [rng.vec(speed), rng.vec(speed)];
        let w = [rng.vec(speed), rng.vec(speed)];
        let d0 = segment_segment_distance(p[0], p[1], q[0], q[1]);
        if d0 < 1e-6 {
            continue;
        }
        let target = OPT.gap * d0;
        let (t, evals) = segment_segment_advance(p, u, q, w, d0, 0.0, 1.0, &OPT);
        let truth = first_root(|t| pair_dist(p, u, q, w, t) - target, samples);
        let end = t.unwrap_or(1.0);
        for s in 0..=samples {
            let tt = end * s as f64 / samples as f64;
            let d = pair_dist(p, u, q, w, tt);
            assert!(
                d >= target * (1.0 - 1e-8),
                "pair {k}: distance {d} < gap {target} at {tt} (ccd {t:?})"
            );
        }
        match (t, truth) {
            (Some(t), Some(r)) => {
                assert!(t <= r + 1e-12, "pair {k}: ccd {t} after the truth {r}");
                contacts += 1;
                if evals > OPT.max_advances {
                    capped += 1;
                }
            }
            (None, Some(r)) => panic!("pair {k}: ccd missed the contact at {r}"),
            (Some(_), None) => capped += 1,
            (None, None) => free += 1,
        }
    }
    eprintln!("random pairs: {contacts} contacts, {free} free, {capped} capped/conservative");
    assert!(contacts > n_pairs / 10 && free > n_pairs / 10);
}

#[test]
fn random_pairs_are_conservative() {
    random_pairs(3000, 400, 1);
}

/// The 10⁵-pair version of §13.5 (`cargo test --release --test ccd -- --ignored`).
#[test]
#[ignore]
fn random_pairs_are_conservative_1e5() {
    random_pairs(100_000, 2000, 7);
}

fn ring(n: usize, r: f64, z: f64, offset: u32) -> (Vec<[u32; 2]>, Vec<P>) {
    let edges = (0..n as u32)
        .map(|i| [offset + i, offset + (i + 1) % n as u32])
        .collect();
    let pos = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            [r * a.cos(), r * a.sin(), z]
        })
        .collect();
    (edges, pos)
}

#[test]
fn network_contact_matches_brute_force() {
    // Two stacked rings; the upper one falls and tilts.
    let (mut e, mut pos) = ring(16, 1.0, 0.0, 0);
    let (e2, p2) = ring(12, 0.8, 1.0, 16);
    e.extend(e2);
    pos.extend(p2);
    let net = CurveNetwork::new(pos.len(), e).unwrap();
    let mut rng = Rng(3);
    let vel: Vec<P> = pos
        .iter()
        .enumerate()
        .map(|(v, _)| {
            let mut x = rng.vec(0.2);
            if v >= 16 {
                x[2] -= 2.0;
            }
            x
        })
        .collect();
    for threads in [1, 7] {
        let pl = pool(threads);
        let r = first_contact(&pl, &net, &pos, &vel, 1.0, &[], &OPT).unwrap();
        let c = r.contact.expect("the rings collide");
        assert!(matches!(c.primitive, ContactPrimitive::EdgeEdge { .. }));
        // Brute force over all pairs in 𝒫 with the sampled truth.
        let edges = net.edges();
        let mut truth = f64::INFINITY;
        for i in 0..edges.len() {
            for j in i + 1..edges.len() {
                if net.edges_adjacent(i, j) {
                    continue;
                }
                let [a, b] = edges[i];
                let [cc, d] = edges[j];
                let p = [pos[a as usize], pos[b as usize]];
                let u = [vel[a as usize], vel[b as usize]];
                let q = [pos[cc as usize], pos[d as usize]];
                let w = [vel[cc as usize], vel[d as usize]];
                let d0 = segment_segment_distance(p[0], p[1], q[0], q[1]);
                if let Some(t) = first_root(|t| pair_dist(p, u, q, w, t) - 0.1 * d0, 2000) {
                    truth = truth.min(t);
                }
                // Nothing touches before the reported time.
                for s in 0..=50 {
                    let t = c.tau * s as f64 / 50.0;
                    assert!(pair_dist(p, u, q, w, t) > 0.0);
                }
            }
        }
        assert!(c.tau <= truth + 1e-12, "{} vs {}", c.tau, truth);
        assert!(
            c.tau > 0.5 * truth,
            "too conservative: {} vs {}",
            c.tau,
            truth
        );
        assert!(r.windows >= 1 && r.candidates > 0);
    }
    // Bitwise independent of the thread count.
    let a = first_contact(&pool(1), &net, &pos, &vel, 1.0, &[], &OPT).unwrap();
    let b = first_contact(&pool(7), &net, &pos, &vel, 1.0, &[], &OPT).unwrap();
    assert_eq!(a.contact, b.contact);
    assert_eq!(a.candidates, b.candidates);
    assert_eq!(a.advances, b.advances);
}

#[test]
fn broad_phase_skips_far_pairs() {
    // A large circle moving rigidly but slowly: few candidates per window, no contact.
    let (e, pos) = ring(2000, 100.0, 0.0, 0);
    let net = CurveNetwork::new(pos.len(), e).unwrap();
    let vel = vec![[1.0, 0.5, 0.0]; pos.len()];
    let r = first_contact(&pool(4), &net, &pos, &vel, 0.01, &[], &OPT).unwrap();
    assert!(r.contact.is_none());
    let m = net.n_edges() as u64;
    assert!(
        r.candidates < 20 * m,
        "{} candidates for {} edges",
        r.candidates,
        m
    );
}

#[test]
fn obstacles_stop_points_and_segments() {
    let tri: Vec<[[f64; 3]; 3]> = vec![[[-1., -1., 0.], [1., -1., 0.], [0., 1., 0.]]];
    let ob: Arc<dyn CollisionObstacle> = Arc::new(tri);
    let prepared = [PreparedObstacle::new(ob).unwrap()];
    // A square loop above the triangle falling onto its face: the vertex at the centre hits.
    let pos: Vec<P> = vec![
        [0., 0., 1.],
        [3., 0., 3.],
        [3., 3., 3.],
        [0., 3., 3.],
        [-3., 3., 3.],
    ];
    let e: Vec<[u32; 2]> = (0..5).map(|i| [i, (i + 1) % 5]).collect();
    let net = CurveNetwork::new(5, e).unwrap();
    let vel = vec![[0., 0., -2.]; 5];
    let r = first_contact(&pool(2), &net, &pos, &vel, 1.0, &prepared, &OPT).unwrap();
    let c = r.contact.unwrap();
    assert!((c.tau - 0.45).abs() < 1e-12, "{c:?}");
    assert!(matches!(
        c.primitive,
        ContactPrimitive::VertexTriangle { vertex: 0, .. }
    ));
    // A long segment sweeping horizontally through the triangle (its endpoints never come close
    // to the face): stopped by an edge–edge contact before it pierces.
    let pos: Vec<P> = vec![[-5., 0., -3.], [-5., 0., 3.], [-8., 9., 3.], [-9., 9., -3.]];
    let net = CurveNetwork::new(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]]).unwrap();
    let vel = vec![[10., 0., 0.]; 4];
    let r = first_contact(&pool(2), &net, &pos, &vel, 1.0, &prepared, &OPT).unwrap();
    let c = r.contact.unwrap();
    assert!(matches!(
        c.primitive,
        ContactPrimitive::EdgeTriangleEdge { edge: 0, .. }
    ));
    // The segment's point (x, 0, 0), x = −5 + 10τ, has distance |2x + 1|/√5 to the triangle's
    // left edge (on the edge's interior); the initial distance is √17 (to the corner (−1, −1)),
    // so the gap 0.1 √17 is reached at 2x + 1 = −0.1 √17 √5.
    let x = -(1.0 + 0.1 * 17f64.sqrt() * 5f64.sqrt()) / 2.0;
    let expected = (x + 5.0) / 10.0;
    assert!(
        c.tau <= expected && c.tau > expected - 1e-8,
        "{c:?} vs {expected}"
    );
}

#[test]
fn bad_input_is_an_error() {
    let (e, pos) = ring(8, 1.0, 0.0, 0);
    let net = CurveNetwork::new(8, e).unwrap();
    let pl = pool(1);
    let vel = vec![[0.0; 3]; 8];
    assert!(matches!(
        first_contact(&pl, &net, &pos[..7], &vel, 1.0, &[], &OPT),
        Err(RcError::PositionCountMismatch { .. })
    ));
    let mut bad = vel.clone();
    bad[3][1] = f64::NAN;
    assert!(matches!(
        first_contact(&pl, &net, &pos, &bad, 1.0, &[], &OPT),
        Err(RcError::InvalidOption { .. })
    ));
    assert!(first_contact(&pl, &net, &pos, &vel, -1.0, &[], &OPT).is_err());
    assert!(first_contact(&pl, &net, &pos, &vel, f64::INFINITY, &[], &OPT).is_err());
    let o = CcdOptions { gap: 0.0, ..OPT };
    assert!(first_contact(&pl, &net, &pos, &vel, 1.0, &[], &o).is_err());
    let nan_tri: Vec<[[f64; 3]; 3]> = vec![[[0., 0., f64::NAN], [1., 0., 0.], [0., 1., 0.]]];
    assert!(PreparedObstacle::new(Arc::new(nan_tri)).is_err());
    // No motion: no contact.
    let r = first_contact(&pl, &net, &pos, &vel, 1.0, &[], &OPT).unwrap();
    assert!(r.contact.is_none());
    // Minimal pair distance of a regular octagon: the distance of edges two apart.
    let d = min_pair_distance(&net, &pos);
    let side = 2.0 * (std::f64::consts::PI / 8.0).sin();
    assert!(d > 0.5 * side && d < 1.5 * side, "{d}");
}

/// A rotation of `R³` from a seed (unit quaternion; deterministic).
fn rotation(rng: &mut Rng) -> [[f64; 3]; 3] {
    let q = [rng.sym(1.0), rng.sym(1.0), rng.sym(1.0), rng.sym(1.0)];
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    let [w, x, y, z] = [q[0] / n, q[1] / n, q[2] / n, q[3] / n];
    [
        [
            1. - 2. * (y * y + z * z),
            2. * (x * y - z * w),
            2. * (x * z + y * w),
        ],
        [
            2. * (x * y + z * w),
            1. - 2. * (x * x + z * z),
            2. * (y * z - x * w),
        ],
        [
            2. * (x * z - y * w),
            2. * (y * z + x * w),
            1. - 2. * (x * x + y * y),
        ],
    ]
}

fn rot(r: &[[f64; 3]; 3], v: P) -> P {
    [
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    ]
}

/// Survey finding 1 (`lean/proofs/FORMALIZATION-SURVEY.md`): two nearly parallel segments whose
/// lines cross at an interior point of both, at height `h`, so the distance is exactly `h`. The
/// former code skipped the interior pair for `denom ≤ 1e-14·a·e` (`θ ≲ 1e-7`) and evaluated the
/// closed form `s = (bf − ce)/denom` badly above it; it overestimated by up to `≈ 3e-6·L`
/// (measured against 80-digit arithmetic), the unsafe direction for the CCD.
#[test]
fn near_parallel_crossing_segments_have_the_exact_distance() {
    let mut rng = Rng(7);
    for &theta in &[
        1e-12, 1e-10, 1e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2,
    ] {
        for &h in &[0.0, 1e-12, 1e-9, 1e-6] {
            for k in 0..8 {
                let r = if k == 0 {
                    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
                } else {
                    rotation(&mut rng)
                };
                let off = rng.vec(1.0);
                let (s0, t0) = (0.2 + 0.6 * rng.next(), 0.2 + 0.6 * rng.next());
                let (c, s) = (f64::cos(theta), f64::sin(theta));
                let m = |v: P| {
                    let x = rot(&r, v);
                    [x[0] + off[0], x[1] + off[1], x[2] + off[2]]
                };
                let p1 = m([-2.0 * s0, 0., 0.]);
                let q1 = m([2.0 * (1.0 - s0), 0., 0.]);
                let p2 = m([-3.0 * t0 * c, -3.0 * t0 * s, h]);
                let q2 = m([3.0 * (1.0 - t0) * c, 3.0 * (1.0 - t0) * s, h]);
                for d in [
                    segment_segment_distance(p1, q1, p2, q2),
                    segment_segment_distance(q1, p1, p2, q2),
                    segment_segment_distance(p2, q2, p1, q1),
                ] {
                    // The inputs are rounded to f64 after the rotation: allow a few ulps of the
                    // coordinates.
                    assert!(
                        (d - h).abs() <= 2e-14,
                        "theta {theta:e}, h {h:e}, rotation {k}: {d:e}"
                    );
                }
            }
        }
    }
}

/// The end-to-end consequence of finding 1: a nearly parallel segment moving through another
/// one. The true distance reaches `0` at `τ = ½`; the CCD must report a contact no later.
/// (With the former distance the pair was measured at `≈ θ` instead of `h₀`, the advance jumped
/// over the crossing and no contact was reported.)
#[test]
fn near_parallel_segment_passing_through_is_caught() {
    let mut rng = Rng(11);
    let z = [0.0; 3];
    for &theta in &[1e-9, 1e-8, 1e-7, 1e-6, 1e-5] {
        for &h0 in &[1e-10, 1e-9, 1e-7] {
            for k in 0..6 {
                let r = if k == 0 {
                    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
                } else {
                    rotation(&mut rng)
                };
                let (c, s) = (f64::cos(theta), f64::sin(theta));
                let p = [rot(&r, [-1., 0., 0.]), rot(&r, [1., 0., 0.])];
                let q = [rot(&r, [-c, -s, h0]), rot(&r, [c, s, h0])];
                // Down with speed 2 h0: the crossing is at τ = ½.
                let v = rot(&r, [0., 0., -2.0 * h0]);
                let d0 = segment_segment_distance(p[0], p[1], q[0], q[1]);
                assert!(
                    (d0 - h0).abs() <= 1e-15,
                    "theta {theta:e}, h0 {h0:e}: d0 {d0:e}"
                );
                let (t, _) = segment_segment_advance(p, [z, z], q, [v, v], d0, 0.0, 1.0, &OPT);
                let t = t.unwrap_or_else(|| panic!("theta {theta:e}, h0 {h0:e}: no contact"));
                assert!(t <= 0.5, "theta {theta:e}, h0 {h0:e}: contact at {t}");
            }
        }
    }
}
