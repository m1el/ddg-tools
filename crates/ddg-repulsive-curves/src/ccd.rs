//! Continuous collision detection (`SPEC.md` §7.4, D11).
//!
//! A *motion* moves every vertex linearly, `γ_v(τ) = γ_v + τ u_v` for `τ ∈ [0, τ_end]`.
//! [`first_contact`] returns a conservative first contact time `τ_max ≤` the first time any two
//! primitives touch:
//!
//! * segment–segment for every pair of curve edges without a common vertex (the pairs `𝒫`);
//! * with obstacles ([`CollisionObstacle`], rigid triangle soups): point–triangle (curve vertex,
//!   obstacle triangle) and segment–segment (curve edge, triangle edge). A moving segment can
//!   only start to pierce a triangle by touching one of its edges or by an endpoint touching its
//!   face, so these primitives see every first contact.
//!
//! **Narrow phase: conservative advancement** (Additive CCD, Li, Kaufman, Jiang, *Codimensional
//! Incremental Potential Contact*, ACM TOG 40(4), 2021). The distance of two linearly moving
//! primitives is Lipschitz in `τ` with constant `l` = the sum of the largest corner speeds on
//! each side, so `dist(τ) ≥ dist(τ_k) − (τ − τ_k) l` and the advance
//! `τ_{k+1} = τ_k + (dist(τ_k) − s d₀)/l` keeps `dist ≥ s d₀` on `[τ_k, τ_{k+1}]` (`s` the gap,
//! `d₀` the pair's distance at `τ = 0`). As in Additive CCD, the corner velocities are taken
//! relative to the mean velocity of the pair's corners: a common translation does not change the
//! distance, and a pair that translates rigidly gets `l = 0`. Each pair advances until its
//! distance reaches `s d₀` (contact reported at that time), the end of the time window
//! (no contact), or `max_advances` steps (contact reported at the current, still safe, time).
//!
//! **Broad phase.** The time range is cut into windows `[0, t₁], [t₁, 2t₁], [2t₁, 4t₁], …` with
//! `t₁` the time in which the fastest vertex moves one mean edge length. In each window a
//! bounding-box tree over the swept boxes of the edges (a linearly moving segment stays in the
//! convex hull of its four corner positions) gives the candidate pairs; the first window with a
//! contact gives `τ_max`. Pairs whose swept boxes are disjoint cannot touch in the window. Cost
//! per window: `O(|E| log |E| + candidates)`, so short motions on large curves do not visit all
//! `|E|²` pairs. (A simple box tree; the generic BVH of `spatial` may replace it.)
//!
//! Distances are exact closest-point distances (Ericson, *Real-Time Collision Detection*, 2004,
//! §5.1.9 segment–segment and §5.1.5 point–triangle), computed as the minimum over the
//! candidate closest-point pairs; they use no coplanarity polynomials and work unchanged in
//! planar mode.
//!
//! **Determinism.** Every pair's result depends only on its own data; per-element results are
//! reduced in index order, so the result is bitwise independent of the thread count (D25).
//!
//! [formula: rc-ccd-advancement]
//! [formula: rc-accd]
//! [formula: rc-closest-points]

use std::sync::Arc;
use std::time::Instant;

use rayon::prelude::*;

use crate::error::RcError;
use crate::kernel::v3::{add, cross, dot, norm, scale, sub};
use crate::network::CurveNetwork;

/// A rigid obstacle for the CCD: a triangle soup (no manifold requirement). Implemented by the
/// obstacles of work package WP7 (`obstacle::Obstacle`); a plain `Vec<[[f64; 3]; 3]>`
/// implements it too.
pub trait CollisionObstacle: Send + Sync {
    /// The triangles, corner coordinates.
    fn triangles(&self) -> &[[[f64; 3]; 3]];
}

impl CollisionObstacle for Vec<[[f64; 3]; 3]> {
    fn triangles(&self) -> &[[[f64; 3]; 3]] {
        self
    }
}

/// Options of the CCD.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CcdOptions {
    /// The gap `s ∈ (0, 1)`: contact is reported when a pair's distance reaches `s d₀`.
    pub gap: f64,
    /// Maximal number of advances per pair and window.
    pub max_advances: u32,
}

impl Default for CcdOptions {
    /// Gap 0.1, 100 advances (`SPEC.md` §7.4).
    fn default() -> Self {
        CcdOptions {
            gap: 0.1,
            max_advances: 100,
        }
    }
}

impl CcdOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a gap outside `(0, 1)` or `max_advances = 0`.
    pub fn validate(&self) -> Result<(), RcError> {
        if !(self.gap > 0.0 && self.gap < 1.0) {
            return Err(RcError::InvalidOption {
                name: "ccd.gap",
                value: self.gap.to_string(),
                requirement: "0 < gap < 1",
            });
        }
        if self.max_advances == 0 {
            return Err(RcError::InvalidOption {
                name: "ccd.max_advances",
                value: "0".into(),
                requirement: "at least one advance",
            });
        }
        Ok(())
    }
}

/// The two primitives of a contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactPrimitive {
    /// Curve edges `a < b` (no common vertex).
    EdgeEdge {
        /// First edge.
        a: u32,
        /// Second edge.
        b: u32,
    },
    /// A curve vertex and an obstacle triangle.
    VertexTriangle {
        /// The curve vertex.
        vertex: u32,
        /// Index of the obstacle in the list passed to [`first_contact`].
        obstacle: usize,
        /// The triangle.
        triangle: u32,
    },
    /// A curve edge and an edge of an obstacle triangle.
    EdgeTriangleEdge {
        /// The curve edge.
        edge: u32,
        /// Index of the obstacle.
        obstacle: usize,
        /// The triangle.
        triangle: u32,
        /// The triangle edge: corners `(side, side + 1 mod 3)`.
        side: u8,
    },
}

/// A (conservative) contact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    /// The time `τ_max`: no two primitives touch on `[0, τ_max]`, and the reporting pair's
    /// distance has come down to (about) `s d₀` there, or its advancement ran out of steps.
    pub tau: f64,
    /// The pair.
    pub primitive: ContactPrimitive,
}

/// Result of a CCD query.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CcdReport {
    /// The first contact, `None` if none before the end of the motion.
    pub contact: Option<Contact>,
    /// Number of candidate pairs given to the narrow phase (all windows).
    pub candidates: u64,
    /// Number of advances (distance evaluations) of the narrow phase.
    pub advances: u64,
    /// Number of time windows visited.
    pub windows: u32,
    /// Wall time.
    pub seconds: f64,
}

impl CcdReport {
    /// Adds the counters of `other` (contact not merged).
    pub fn accumulate(&mut self, other: &CcdReport) {
        self.candidates += other.candidates;
        self.advances += other.advances;
        self.windows += other.windows;
        self.seconds += other.seconds;
    }
}

/// An obstacle with its broad-phase trees (built once, the obstacle is rigid).
pub struct PreparedObstacle {
    obstacle: Arc<dyn CollisionObstacle>,
    triangle_tree: BoxTree,
    edge_tree: BoxTree,
}

impl std::fmt::Debug for PreparedObstacle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PreparedObstacle")
            .field("triangles", &self.obstacle.triangles().len())
            .finish_non_exhaustive()
    }
}

impl PreparedObstacle {
    /// Builds the trees over the triangles and over their edges.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-finite coordinate (naming the triangle).
    pub fn new(obstacle: Arc<dyn CollisionObstacle>) -> Result<PreparedObstacle, RcError> {
        let tris = obstacle.triangles();
        if let Some(t) = tris
            .iter()
            .position(|t| !t.iter().flatten().all(|c| c.is_finite()))
        {
            return Err(RcError::InvalidOption {
                name: "obstacle triangle",
                value: t.to_string(),
                requirement: "finite corner coordinates",
            });
        }
        let tri_boxes: Vec<Aabb> = tris
            .iter()
            .map(|t| Aabb::of_points(&[t[0], t[1], t[2]]))
            .collect();
        let edge_boxes: Vec<Aabb> = tris
            .iter()
            .flat_map(|t| (0..3).map(move |s| Aabb::of_points(&[t[s], t[(s + 1) % 3]])))
            .collect();
        Ok(PreparedObstacle {
            triangle_tree: BoxTree::build(&tri_boxes),
            edge_tree: BoxTree::build(&edge_boxes),
            obstacle,
        })
    }

    /// The obstacle.
    pub fn obstacle(&self) -> &Arc<dyn CollisionObstacle> {
        &self.obstacle
    }
}

// ---------------------------------------------------------------------------------------------
// Distances (Ericson 2004, §5.1.2, §5.1.5, §5.1.9)
// ---------------------------------------------------------------------------------------------

/// Distance from `p` to the segment `a b`.
///
/// [formula: rc-closest-points]
pub fn point_segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let ap = sub(p, a);
    let e = dot(ab, ab);
    let t = if e > 0.0 {
        (dot(ap, ab) / e).clamp(0.0, 1.0)
    } else {
        0.0
    };
    norm(sub(ap, scale(ab, t)))
}

/// Distance between the segments `p1 q1` and `p2 q2`: the minimum over the four
/// endpoint–segment distances and the pair `(s, t)` found from the closest points of the lines:
/// `s` clamped to `[0, 1]`, then `t` the clamped minimizer for that `s`, then `s` the clamped
/// minimizer for that `t` (Ericson's re-projection). Every candidate is the distance of an actual
/// pair of points. When the lines are not parallel and their closest points lie on both
/// segments, the pair is that interior pair, so in exact arithmetic the result is the distance
/// (`DDG.RC.segSegCode_eq`, `segSegDist_le`, `segSegDist_attained` in `lean/proofs`,
/// `RepulsiveCurves/ClosestPoints.lean`; the cross-product form: `cross_formulas_segSeg`).
///
/// The line parameter is `s = ((p₂ − p₁) × d₂)·n / |n|²` with `n = d₁ × d₂`, which equals the
/// closed form `(bf − ce)/(ae − b²)` (Lagrange's identity) but is well conditioned for nearly
/// parallel lines. The former code skipped the interior pair for `ae − b² ≤ 1e-14·ae` and used
/// the closed form above it; for nearly parallel crossing segments it overestimated the
/// distance (in exact arithmetic `d̃² ≤ d² + ¼·1e-14·min(a, e)`, tight; in f64 up to `≈ 3e-6·L`
/// for angles near `1e-6`), the unsafe direction for the conservative advancement.
///
/// [formula: rc-closest-points]
pub fn segment_segment_distance(p1: [f64; 3], q1: [f64; 3], p2: [f64; 3], q2: [f64; 3]) -> f64 {
    let mut best = point_segment_distance(p1, p2, q2)
        .min(point_segment_distance(q1, p2, q2))
        .min(point_segment_distance(p2, p1, q1))
        .min(point_segment_distance(q2, p1, q1));
    let d1 = sub(q1, p1);
    let d2 = sub(q2, p2);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let n = cross(d1, d2);
    let nn = dot(n, n);
    if a > 0.0 && e > 0.0 && nn > 0.0 {
        let r = sub(p1, p2);
        let b = dot(d1, d2);
        let c = dot(d1, r);
        let f = dot(d2, r);
        // `clamp` returns NaN for NaN (overflowing inputs); `min` below then ignores the pair.
        let s = (dot(cross(sub(p2, p1), d2), n) / nn).clamp(0.0, 1.0);
        let t = ((b * s + f) / e).clamp(0.0, 1.0);
        let s = ((b * t - c) / a).clamp(0.0, 1.0);
        let c1 = add(p1, scale(d1, s));
        let c2 = add(p2, scale(d2, t));
        best = best.min(norm(sub(c1, c2)));
    }
    best
}

/// Distance from `p` to the triangle `a b c` (degenerate triangles: distance to their edges). In
/// exact arithmetic it is the distance for every triangle (`DDG.RC.ptTriDist_le`,
/// `ptTriDist_attained`; the cross-product forms: `cross_formulas_ptTri`).
///
/// [formula: rc-closest-points]
pub fn point_triangle_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let edges = point_segment_distance(p, a, b)
        .min(point_segment_distance(p, b, c))
        .min(point_segment_distance(p, c, a));
    let ab = sub(b, a);
    let ac = sub(c, a);
    let n = cross(ab, ac);
    let nn = dot(n, n);
    if nn <= 0.0 {
        return edges;
    }
    // Barycentric coordinates of the projection of p onto the plane.
    let ap = sub(p, a);
    let v = dot(cross(ap, ac), n) / nn;
    let w = dot(cross(ab, ap), n) / nn;
    if v >= 0.0 && w >= 0.0 && v + w <= 1.0 {
        let h = dot(ap, n).abs() / nn.sqrt();
        edges.min(h)
    } else {
        edges
    }
}

// ---------------------------------------------------------------------------------------------
// Narrow phase
// ---------------------------------------------------------------------------------------------

/// Conservative advancement of one pair on `[a, b]` from the distance function `dist`, with
/// Lipschitz constant `l`, target distance `target = s d₀`. Returns the contact time (if the
/// distance reaches the target within the window or the advances run out) and the number of
/// distance evaluations.
///
/// [formula: rc-ccd-advancement]
fn advance(
    dist: impl Fn(f64) -> f64,
    l: f64,
    target: f64,
    d0: f64,
    a: f64,
    b: f64,
    max_advances: u32,
) -> (Option<f64>, u32) {
    let mut tau = a;
    let mut d = dist(tau);
    let mut evals = 1;
    // Reaching the target up to rounding counts as reaching it.
    let slack = 1e-9 * d0;
    for _ in 0..max_advances {
        if d.is_nan() || d - target <= slack {
            return (Some(tau), evals);
        }
        if l <= 0.0 {
            return (None, evals);
        }
        let next = tau + (d - target) / l;
        if next >= b {
            return (None, evals);
        }
        tau = next;
        d = dist(tau);
        evals += 1;
    }
    // Reached the target, or out of advances: the current time is still safe.
    (Some(tau), evals)
}

#[inline]
fn at(p: [f64; 3], u: [f64; 3], t: f64) -> [f64; 3] {
    [p[0] + t * u[0], p[1] + t * u[1], p[2] + t * u[2]]
}

/// Largest speed of the corners relative to their mean velocity, split by side.
fn relative_lipschitz(side_a: &[[f64; 3]], side_b: &[[f64; 3]]) -> f64 {
    let n = (side_a.len() + side_b.len()) as f64;
    let mut mean = [0.0; 3];
    for v in side_a.iter().chain(side_b) {
        mean = add(mean, *v);
    }
    mean = scale(mean, 1.0 / n);
    let la = side_a
        .iter()
        .map(|v| norm(sub(*v, mean)))
        .fold(0.0f64, f64::max);
    let lb = side_b
        .iter()
        .map(|v| norm(sub(*v, mean)))
        .fold(0.0f64, f64::max);
    // A relative margin against rounding in the speeds.
    (la + lb) * (1.0 + 1e-12)
}

/// Conservative first contact time on `[t_start, t_end]` of the segments `p0 p1` (velocities
/// `u0, u1`) and `q0 q1` (velocities `w0, w1`), positions at `τ` being `p + τ u`. `d0` is the
/// pair's distance at `τ = 0` (the reference of the gap). Returns `None` if the distance stays
/// above `gap · d0` on the whole interval, and the number of distance evaluations.
///
/// [formula: rc-ccd-advancement]
#[allow(clippy::too_many_arguments)]
pub fn segment_segment_advance(
    p: [[f64; 3]; 2],
    u: [[f64; 3]; 2],
    q: [[f64; 3]; 2],
    w: [[f64; 3]; 2],
    d0: f64,
    t_start: f64,
    t_end: f64,
    options: &CcdOptions,
) -> (Option<f64>, u32) {
    let l = relative_lipschitz(&u, &w);
    advance(
        |t| {
            segment_segment_distance(
                at(p[0], u[0], t),
                at(p[1], u[1], t),
                at(q[0], w[0], t),
                at(q[1], w[1], t),
            )
        },
        l,
        options.gap * d0,
        d0,
        t_start,
        t_end,
        options.max_advances,
    )
}

/// As [`segment_segment_advance`] for a moving point `p + τ u` and a fixed triangle.
///
/// [formula: rc-ccd-advancement]
#[allow(clippy::too_many_arguments)]
pub fn point_triangle_advance(
    p: [f64; 3],
    u: [f64; 3],
    tri: [[f64; 3]; 3],
    d0: f64,
    t_start: f64,
    t_end: f64,
    options: &CcdOptions,
) -> (Option<f64>, u32) {
    let l = norm(u) * (1.0 + 1e-12);
    advance(
        |t| point_triangle_distance(at(p, u, t), tri[0], tri[1], tri[2]),
        l,
        options.gap * d0,
        d0,
        t_start,
        t_end,
        options.max_advances,
    )
}

// ---------------------------------------------------------------------------------------------
// Broad phase
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
struct Aabb {
    lo: [f64; 3],
    hi: [f64; 3],
}

impl Aabb {
    fn of_points(ps: &[[f64; 3]]) -> Aabb {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in ps {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        Aabb { lo, hi }
    }

    fn overlaps(&self, o: &Aabb) -> bool {
        (0..3).all(|k| self.lo[k] <= o.hi[k] && o.lo[k] <= self.hi[k])
    }
}

/// A median-split box tree (leaves of at most 4 boxes).
struct BoxTree {
    nodes: Vec<BoxNode>,
    perm: Vec<u32>,
}

struct BoxNode {
    bounds: Aabb,
    /// Leaf: range of `perm`; inner: the two children.
    a: u32,
    b: u32,
    leaf: bool,
}

impl BoxTree {
    const LEAF: usize = 4;

    fn build(boxes: &[Aabb]) -> BoxTree {
        let mut perm: Vec<u32> = (0..boxes.len() as u32).collect();
        let mut nodes = Vec::with_capacity(2 * boxes.len() / Self::LEAF + 1);
        if !boxes.is_empty() {
            Self::build_rec(boxes, &mut perm, 0, boxes.len(), &mut nodes);
        }
        BoxTree { nodes, perm }
    }

    fn build_rec(
        boxes: &[Aabb],
        perm: &mut [u32],
        start: usize,
        end: usize,
        nodes: &mut Vec<BoxNode>,
    ) -> u32 {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &i in &perm[start..end] {
            let b = boxes[i as usize];
            for k in 0..3 {
                lo[k] = lo[k].min(b.lo[k]);
                hi[k] = hi[k].max(b.hi[k]);
            }
        }
        let id = nodes.len() as u32;
        nodes.push(BoxNode {
            bounds: Aabb { lo, hi },
            a: start as u32,
            b: end as u32,
            leaf: true,
        });
        if end - start > Self::LEAF {
            let ext = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
            let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] {
                0
            } else if ext[1] >= ext[2] {
                1
            } else {
                2
            };
            let centre = |i: u32| boxes[i as usize].lo[axis] + boxes[i as usize].hi[axis];
            perm[start..end]
                .sort_unstable_by(|&x, &y| centre(x).total_cmp(&centre(y)).then(x.cmp(&y)));
            let mid = (start + end) / 2;
            let l = Self::build_rec(boxes, perm, start, mid, nodes);
            let r = Self::build_rec(boxes, perm, mid, end, nodes);
            let n = &mut nodes[id as usize];
            n.leaf = false;
            n.a = l;
            n.b = r;
        }
        id
    }

    /// Indices of the boxes overlapping `q`, appended to `out` (order: deterministic).
    fn query(&self, q: &Aabb, stack: &mut Vec<u32>, out: &mut Vec<u32>) {
        if self.nodes.is_empty() {
            return;
        }
        stack.clear();
        stack.push(0);
        while let Some(id) = stack.pop() {
            let n = &self.nodes[id as usize];
            if !n.bounds.overlaps(q) {
                continue;
            }
            if n.leaf {
                out.extend_from_slice(&self.perm[n.a as usize..n.b as usize]);
            } else {
                stack.push(n.a);
                stack.push(n.b);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Queries
// ---------------------------------------------------------------------------------------------

fn bad(name: &'static str, value: String, requirement: &'static str) -> RcError {
    RcError::InvalidOption {
        name,
        value,
        requirement,
    }
}

/// Earlier of two contacts; ties keep `a` (the earlier one in index order).
fn earlier(a: Option<Contact>, b: Option<Contact>) -> Option<Contact> {
    match (a, b) {
        (Some(x), Some(y)) => Some(if y.tau < x.tau { y } else { x }),
        (x, None) => x,
        (None, y) => y,
    }
}

/// The first contact of the motion `γ_v(τ) = start_v + τ velocity_v`, `τ ∈ [0, tau_end]`, of the
/// curve against itself (pairs of edges without a common vertex) and against `obstacles`
/// (rigid). See the module documentation for the method and its guarantee: no two primitives
/// touch on `[0, contact.tau]`, and if `contact` is `None`, on `[0, tau_end]`.
///
/// Runs on `pool`; the result is bitwise independent of its size.
///
/// # Errors
/// [`RcError::PositionCountMismatch`] for arrays of the wrong length; [`RcError::InvalidOption`]
/// for invalid options, a negative or non-finite `tau_end`, or non-finite positions or
/// velocities.
#[allow(clippy::too_many_arguments)]
pub fn first_contact(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    start: &[[f64; 3]],
    velocity: &[[f64; 3]],
    tau_end: f64,
    obstacles: &[PreparedObstacle],
    options: &CcdOptions,
) -> Result<CcdReport, RcError> {
    let t0 = Instant::now();
    options.validate()?;
    let n = network.n_vertices();
    for arr in [start, velocity] {
        if arr.len() != n {
            return Err(RcError::PositionCountMismatch {
                expected: n,
                found: arr.len(),
            });
        }
    }
    if !(tau_end.is_finite() && tau_end >= 0.0) {
        return Err(bad(
            "ccd time range",
            tau_end.to_string(),
            "a finite non-negative end time",
        ));
    }
    if let Some(v) = start
        .iter()
        .chain(velocity)
        .position(|p| !p.iter().all(|c| c.is_finite()))
    {
        return Err(bad(
            "ccd motion",
            format!("vertex {}", v % n),
            "finite positions and velocities",
        ));
    }
    let mut report = CcdReport::default();
    let vmax = velocity.iter().map(|v| norm(*v)).fold(0.0f64, f64::max);
    if vmax == 0.0 || tau_end == 0.0 {
        report.seconds = t0.elapsed().as_secs_f64();
        return Ok(report);
    }
    let mean_len = network.mean_edge_length(start);
    let t1 = (mean_len / vmax).min(tau_end);
    let edges = network.edges();
    let mut a = 0.0;
    let mut b = if t1 > 0.0 { t1 } else { tau_end };
    loop {
        report.windows += 1;
        let (contact, cand, adv) = window_contact(
            pool, network, edges, start, velocity, a, b, obstacles, options,
        );
        report.candidates += cand;
        report.advances += adv;
        if contact.is_some() {
            report.contact = contact;
            break;
        }
        if b >= tau_end {
            break;
        }
        a = b;
        b = (2.0 * b).min(tau_end);
    }
    report.seconds = t0.elapsed().as_secs_f64();
    Ok(report)
}

/// One time window of [`first_contact`]: (earliest contact, candidates, advances).
#[allow(clippy::too_many_arguments)]
fn window_contact(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    edges: &[[u32; 2]],
    start: &[[f64; 3]],
    velocity: &[[f64; 3]],
    a: f64,
    b: f64,
    obstacles: &[PreparedObstacle],
    options: &CcdOptions,
) -> (Option<Contact>, u64, u64) {
    let swept: Vec<Aabb> = edges
        .iter()
        .map(|&[i, j]| {
            let (i, j) = (i as usize, j as usize);
            Aabb::of_points(&[
                at(start[i], velocity[i], a),
                at(start[j], velocity[j], a),
                at(start[i], velocity[i], b),
                at(start[j], velocity[j], b),
            ])
        })
        .collect();
    let tree = BoxTree::build(&swept);
    // Edge pairs: per edge I, the pairs (I, J) with J > I.
    let per_edge: Vec<(Option<Contact>, u64, u64)> = pool.install(|| {
        (0..edges.len())
            .into_par_iter()
            .map_init(
                || (Vec::new(), Vec::new()),
                |(stack, hits), i| {
                    hits.clear();
                    tree.query(&swept[i], stack, hits);
                    hits.sort_unstable();
                    let [i1, i2] = edges[i];
                    let p = [start[i1 as usize], start[i2 as usize]];
                    let u = [velocity[i1 as usize], velocity[i2 as usize]];
                    let mut best: Option<Contact> = None;
                    let (mut cand, mut adv) = (0u64, 0u64);
                    for &j in hits.iter() {
                        let j = j as usize;
                        if j <= i || network.edges_adjacent(i, j) {
                            continue;
                        }
                        let [j1, j2] = edges[j];
                        let q = [start[j1 as usize], start[j2 as usize]];
                        let w = [velocity[j1 as usize], velocity[j2 as usize]];
                        let d0 = segment_segment_distance(p[0], p[1], q[0], q[1]);
                        cand += 1;
                        let (t, e) = segment_segment_advance(p, u, q, w, d0, a, b, options);
                        adv += e as u64;
                        if let Some(tau) = t {
                            best = earlier(
                                best,
                                Some(Contact {
                                    tau,
                                    primitive: ContactPrimitive::EdgeEdge {
                                        a: i as u32,
                                        b: j as u32,
                                    },
                                }),
                            );
                        }
                    }
                    (best, cand, adv)
                },
            )
            .collect()
    });
    let mut best = None;
    let (mut cand, mut adv) = (0u64, 0u64);
    for (c, k, e) in per_edge {
        best = earlier(best, c);
        cand += k;
        adv += e;
    }
    for (oi, ob) in obstacles.iter().enumerate() {
        let (c, k, e) = obstacle_window(pool, edges, start, velocity, a, b, oi, ob, options);
        best = earlier(best, c);
        cand += k;
        adv += e;
    }
    (best, cand, adv)
}

/// Curve against one rigid obstacle in a window.
#[allow(clippy::too_many_arguments)]
fn obstacle_window(
    pool: &rayon::ThreadPool,
    edges: &[[u32; 2]],
    start: &[[f64; 3]],
    velocity: &[[f64; 3]],
    a: f64,
    b: f64,
    oi: usize,
    ob: &PreparedObstacle,
    options: &CcdOptions,
) -> (Option<Contact>, u64, u64) {
    let tris = ob.obstacle.triangles();
    let n = start.len();
    let zero = [[0.0; 3]; 2];
    // Vertices against triangles.
    let per_vertex: Vec<(Option<Contact>, u64, u64)> = pool.install(|| {
        (0..n)
            .into_par_iter()
            .map_init(
                || (Vec::new(), Vec::new()),
                |(stack, hits), v| {
                    hits.clear();
                    let bx = Aabb::of_points(&[
                        at(start[v], velocity[v], a),
                        at(start[v], velocity[v], b),
                    ]);
                    ob.triangle_tree.query(&bx, stack, hits);
                    hits.sort_unstable();
                    let mut best = None;
                    let (mut cand, mut adv) = (0u64, 0u64);
                    for &t in hits.iter() {
                        let tri = tris[t as usize];
                        let d0 = point_triangle_distance(start[v], tri[0], tri[1], tri[2]);
                        cand += 1;
                        let (c, e) =
                            point_triangle_advance(start[v], velocity[v], tri, d0, a, b, options);
                        adv += e as u64;
                        if let Some(tau) = c {
                            best = earlier(
                                best,
                                Some(Contact {
                                    tau,
                                    primitive: ContactPrimitive::VertexTriangle {
                                        vertex: v as u32,
                                        obstacle: oi,
                                        triangle: t,
                                    },
                                }),
                            );
                        }
                    }
                    (best, cand, adv)
                },
            )
            .collect()
    });
    // Curve edges against triangle edges.
    let per_edge: Vec<(Option<Contact>, u64, u64)> = pool.install(|| {
        (0..edges.len())
            .into_par_iter()
            .map_init(
                || (Vec::new(), Vec::new()),
                |(stack, hits), i| {
                    hits.clear();
                    let [i1, i2] = edges[i];
                    let (i1, i2) = (i1 as usize, i2 as usize);
                    let bx = Aabb::of_points(&[
                        at(start[i1], velocity[i1], a),
                        at(start[i2], velocity[i2], a),
                        at(start[i1], velocity[i1], b),
                        at(start[i2], velocity[i2], b),
                    ]);
                    ob.edge_tree.query(&bx, stack, hits);
                    hits.sort_unstable();
                    let p = [start[i1], start[i2]];
                    let u = [velocity[i1], velocity[i2]];
                    let mut best = None;
                    let (mut cand, mut adv) = (0u64, 0u64);
                    for &k in hits.iter() {
                        let (t, side) = (k / 3, (k % 3) as usize);
                        let tri = tris[t as usize];
                        let q = [tri[side], tri[(side + 1) % 3]];
                        let d0 = segment_segment_distance(p[0], p[1], q[0], q[1]);
                        cand += 1;
                        let (c, e) = segment_segment_advance(p, u, q, zero, d0, a, b, options);
                        adv += e as u64;
                        if let Some(tau) = c {
                            best = earlier(
                                best,
                                Some(Contact {
                                    tau,
                                    primitive: ContactPrimitive::EdgeTriangleEdge {
                                        edge: i as u32,
                                        obstacle: oi,
                                        triangle: t,
                                        side: side as u8,
                                    },
                                }),
                            );
                        }
                    }
                    (best, cand, adv)
                },
            )
            .collect()
    });
    let mut best = None;
    let (mut cand, mut adv) = (0u64, 0u64);
    for (c, k, e) in per_vertex.into_iter().chain(per_edge) {
        best = earlier(best, c);
        cand += k;
        adv += e;
    }
    (best, cand, adv)
}

/// Checks a piecewise linear path of configurations `path[0] → path[1] → …` (each piece a
/// linear motion over `τ ∈ [0, 1]`). Returns the first piece with a contact and the contact, or
/// `None` if the whole path is free; the reports of all pieces checked are accumulated into
/// `report`.
///
/// # Errors
/// As [`first_contact`].
pub fn path_contact(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    path: &[Vec<[f64; 3]>],
    obstacles: &[PreparedObstacle],
    options: &CcdOptions,
    report: &mut CcdReport,
) -> Result<Option<(usize, Contact)>, RcError> {
    for (k, w) in path.windows(2).enumerate() {
        let vel: Vec<[f64; 3]> = w[0].iter().zip(&w[1]).map(|(p, q)| sub(*q, *p)).collect();
        let r = first_contact(pool, network, &w[0], &vel, 1.0, obstacles, options)?;
        report.accumulate(&r);
        if let Some(c) = r.contact {
            return Ok(Some((k, c)));
        }
    }
    Ok(None)
}

/// The smallest distance over the pairs of edges without a common vertex (broad phase by a box
/// tree with a growing search radius; exact segment distances). For audits and tests.
pub fn min_pair_distance(network: &CurveNetwork, positions: &[[f64; 3]]) -> f64 {
    let edges = network.edges();
    let boxes: Vec<Aabb> = edges
        .iter()
        .map(|&[i, j]| Aabb::of_points(&[positions[i as usize], positions[j as usize]]))
        .collect();
    let tree = BoxTree::build(&boxes);
    let mut radius = network.mean_edge_length(positions);
    let (mut stack, mut hits) = (Vec::new(), Vec::new());
    loop {
        let mut best = f64::INFINITY;
        for (i, bx) in boxes.iter().enumerate() {
            let q = Aabb {
                lo: [bx.lo[0] - radius, bx.lo[1] - radius, bx.lo[2] - radius],
                hi: [bx.hi[0] + radius, bx.hi[1] + radius, bx.hi[2] + radius],
            };
            hits.clear();
            tree.query(&q, &mut stack, &mut hits);
            for &j in &hits {
                let j = j as usize;
                if j <= i || network.edges_adjacent(i, j) {
                    continue;
                }
                let [a, b] = edges[i];
                let [c, d] = edges[j];
                best = best.min(segment_segment_distance(
                    positions[a as usize],
                    positions[b as usize],
                    positions[c as usize],
                    positions[d as usize],
                ));
            }
        }
        if best <= radius || (!best.is_finite() && radius > 1e300) {
            return best;
        }
        if best.is_finite() {
            // Every pair closer than `best` is found with radius `best`.
            radius = best;
        } else if radius > 0.0 {
            radius *= 4.0;
        } else {
            radius = 1.0;
        }
    }
}
