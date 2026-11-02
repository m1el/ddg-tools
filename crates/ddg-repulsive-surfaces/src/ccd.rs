//! Continuous collision detection for moving triangle meshes (`SPEC.md` §7.3, D20).
//!
//! A *motion* moves every vertex linearly, `f_v(t) = start_v + t · velocity_v` for
//! `t ∈ [0, tau_end]` (the flow passes `velocity = −x`). [`first_contact`] returns the largest
//! `τ_max ≤ tau_end` such that on `[0, τ_max]`
//!
//! * every vertex–triangle pair `(v, f)` with `v ∉ f` and every edge–edge pair without a common
//!   vertex of the mesh, and
//! * every pair of a mesh primitive and a static [`CcdObstacle`] (mesh vertex against obstacle
//!   triangle, obstacle vertex against mesh triangle, edge against edge)
//!
//! stays at distance `≥ gap` (see below for the pairs that start closer). This is the
//! conservative version of "no two non-adjacent triangles touch": at the first touching time of
//! two triangles, a vertex lies on a triangle or two edges meet, and elements that share a vertex
//! are not tested.
//!
//! # Gap and initially close pairs
//!
//! The gap is `opts.gap · ℓ̄` with `ℓ̄ = opts.reference_length` (the session passes the mean edge
//! length `ℓ̄₀` of the initial mesh; `None` takes the mean edge length of `start`). A pair that
//! starts closer than the gap must not freeze the flow, so, as in the curves crate, its target is
//! measured against its *current* distance `d₀`: such a pair is only stopped when its distance
//! has come down to `close_factor · d₀` (default one half). A pair at distance `0` at `t = 0`
//! (the mesh already touches itself or the obstacle) cannot be measured; it is not constrained
//! and counted in [`CcdReport::touching_pairs`]. Both counts are in the report.
//!
//! # Inputs that already intersect: the exempt band
//!
//! A mesh that intersects itself at `t = 0` (an immersed Klein bottle, the subcritical runs of the
//! paper's Fig. 26) cannot be held to the gap rule: the intersection curve slides along the mesh
//! and every pair near it closes. [`IntersectionBand::new`] finds the faces that intersect
//! (exactly, [`intersecting_pairs`], BVH-accelerated), takes their vertices widened by one ring
//! (`I`), and [`first_contact_with_band`] leaves a pair unconstrained when both of its elements
//! have a vertex in `I` (counted in [`CcdReport::exempt_pairs`]); every other pair, in particular
//! one between two regions that were disjoint at the start, keeps the gap rule. A band of an
//! embedded mesh is empty and changes nothing (bitwise). [`CcdCollapseGuard::set_band`] applies
//! the same rule to collapses. Mesh-obstacle pairs are never exempt.
//!
//! # Narrow phase: conservative advancement with moving triangles
//!
//! The distance `d(t)` of two linearly moving primitives (convex hulls of their corners) is
//! Lipschitz with constant `l`, the sum of the largest corner speeds of the two sides measured
//! in a common frame (`max |u_p| + max |u_tri|`; a common translation does not change `d`, so
//! the better of the frame at rest and the frame of the mean corner velocity is used). Hence
//! `d(t) ≥ d(t_k) − (t − t_k) l`, and `t_{k+1} = t_k + (d(t_k) − target)/l` keeps
//! `d ≥ target` on `[t_k, t_{k+1}]`; the advance stops when `d − target` is zero up to rounding,
//! at the end of the window, or after `max_iterations` advances (the current time, still safe,
//! is reported).
//!
//! # Broad phase
//!
//! The time range is cut into windows `[0, t₁], [t₁, 2t₁], [2t₁, 4t₁], …` where `t₁` is the time
//! in which the fastest vertex moves one mean edge length. In each window two box trees (Morton-sorted blocks merged into a binary tree, built in parallel) are built
//! over the swept boxes (hull of the corner positions at both ends of the window, inflated by
//! `gap/2`) of the triangles and of the edges; vertices and edges query them, obstacle trees are
//! static. The first window that contains a contact gives `τ_max` (the minimum over its pairs).
//!
//! # Determinism
//!
//! Every pair's result depends only on its own data; results are combined with the total order
//! `(τ, primitive)` and with integer sums, so the report is bitwise independent of the thread
//! count.
//!
//! [formula: rs-ccd-advancement]

use std::collections::HashMap;
use std::time::Instant;

use ddg_repulsive_curves::Obstacle;
use ddg_repulsive_curves::ccd::{point_triangle_distance, segment_segment_distance};
use ddg_repulsive_curves::spatial::Aabb;
use rayon::prelude::*;

use crate::error::RsError;
use crate::geometry::mean_edge_length;
use crate::geometry::v3::{add, cross, dot, norm, scale, sub};
use crate::remesh::{CollapseGuard, CollapseQuery};
use crate::surface::Surface;

// ---------------------------------------------------------------------------------------------
// Options and reports
// ---------------------------------------------------------------------------------------------

/// Options of the CCD.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CcdOptions {
    /// The gap relative to the reference length, `gap > 0` (default `0.05`).
    pub gap: f64,
    /// The reference length `ℓ̄₀` (the mean edge length of the initial mesh); `None`: the mean
    /// edge length of the positions passed to [`first_contact`].
    pub reference_length: Option<f64>,
    /// Maximal number of advances per pair and window (default 100).
    pub max_iterations: u32,
    /// A pair closer than the gap at `t = 0` may close to this fraction of its current distance,
    /// `0 ≤ close_factor < 1` (default `0.5`).
    pub close_factor: f64,
    /// Leaf size of the broad-phase trees (default 8).
    pub leaf_size: usize,
}

impl Default for CcdOptions {
    fn default() -> Self {
        CcdOptions {
            gap: 0.05,
            reference_length: None,
            max_iterations: 100,
            close_factor: 0.5,
            leaf_size: 8,
        }
    }
}

fn bad(name: &'static str, value: String, requirement: &'static str) -> RsError {
    RsError::InvalidOption {
        name,
        value,
        requirement,
    }
}

impl CcdOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] naming the option.
    pub fn validate(&self) -> Result<(), RsError> {
        if !(self.gap.is_finite() && self.gap > 0.0) {
            return Err(bad("ccd.gap", self.gap.to_string(), "finite and positive"));
        }
        if let Some(l) = self.reference_length
            && !(l.is_finite() && l > 0.0)
        {
            return Err(bad(
                "ccd.reference_length",
                l.to_string(),
                "finite and positive",
            ));
        }
        if self.max_iterations == 0 {
            return Err(bad("ccd.max_iterations", "0".into(), "at least 1"));
        }
        if !(self.close_factor >= 0.0 && self.close_factor < 1.0) {
            return Err(bad(
                "ccd.close_factor",
                self.close_factor.to_string(),
                "in [0, 1)",
            ));
        }
        if self.leaf_size == 0 {
            return Err(bad("ccd.leaf_size", "0".into(), "at least 1"));
        }
        Ok(())
    }
}

/// The two primitives of a contact. Ids of the mesh are those of the [`Surface`]; edges are
/// edge ids of [`Surface::edges`]; `obstacle` indexes the slice passed to [`first_contact`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Primitive {
    /// Mesh vertex and mesh triangle.
    VertexTriangle {
        /// The vertex.
        vertex: u32,
        /// The triangle.
        triangle: u32,
    },
    /// Two mesh edges, `a < b`.
    EdgeEdge {
        /// First edge.
        a: u32,
        /// Second edge.
        b: u32,
    },
    /// Mesh vertex and obstacle triangle.
    VertexObstacleTriangle {
        /// The mesh vertex.
        vertex: u32,
        /// The obstacle.
        obstacle: usize,
        /// The obstacle triangle.
        triangle: u32,
    },
    /// Mesh triangle and obstacle vertex.
    TriangleObstacleVertex {
        /// The mesh triangle.
        triangle: u32,
        /// The obstacle.
        obstacle: usize,
        /// The obstacle vertex.
        vertex: u32,
    },
    /// Mesh edge and obstacle edge (index in [`CcdObstacle::edges`]).
    EdgeObstacleEdge {
        /// The mesh edge.
        edge: u32,
        /// The obstacle.
        obstacle: usize,
        /// The obstacle edge.
        obstacle_edge: u32,
    },
}

/// A contact: the earliest safe time of the reporting pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    /// The time `τ_max`.
    pub tau: f64,
    /// The pair whose advancement ended there.
    pub primitive: Primitive,
}

/// Result of [`first_contact`].
#[derive(Debug, Clone, PartialEq)]
pub struct CcdReport {
    /// The largest safe time `≤ tau_end` (`tau_end` itself without contact).
    pub tau_max: f64,
    /// The pair that limits `tau_max`, `None` if `tau_max = tau_end`.
    pub contact: Option<Contact>,
    /// The gap in absolute units, `gap · ℓ̄`.
    pub gap_abs: f64,
    /// Number of distance evaluations of the narrow phase.
    pub iterations: u64,
    /// Number of candidate pairs given to the narrow phase (all windows).
    pub pairs_tested: u64,
    /// Number of time windows visited.
    pub windows: u32,
    /// Number of pairs closer than the gap at `t = 0` (they are measured against their own
    /// distance).
    pub close_pairs: u64,
    /// Number of pairs at distance 0 at `t = 0` (not constrained).
    pub touching_pairs: u64,
    /// Number of broad-phase candidate pairs left unconstrained by the [`IntersectionBand`]
    /// (both elements have a vertex in the band; counted on the first window).
    pub exempt_pairs: u64,
    /// Wall time.
    pub seconds: f64,
}

// ---------------------------------------------------------------------------------------------
// Obstacles
// ---------------------------------------------------------------------------------------------

fn center(b: &Aabb<3>) -> [f64; 3] {
    [
        0.5 * (b.lo[0] + b.hi[0]),
        0.5 * (b.lo[1] + b.hi[1]),
        0.5 * (b.lo[2] + b.hi[2]),
    ]
}

fn inflate(b: &Aabb<3>, h: f64) -> Aabb<3> {
    Aabb {
        lo: [b.lo[0] - h, b.lo[1] - h, b.lo[2] - h],
        hi: [b.hi[0] + h, b.hi[1] + h, b.hi[2] + h],
    }
}

#[inline]
fn overlaps(a: &Aabb<3>, b: &Aabb<3>) -> bool {
    (0..3).all(|k| a.lo[k] <= b.hi[k] && b.lo[k] <= a.hi[k])
}

/// A static tree over boxes: the boxes are sorted along a Morton curve of their centres (ties by
/// index), cut into blocks of `leaf` boxes, and the blocks are merged pairwise into a complete
/// binary tree (`levels[0]` the blocks, the last level the root). The build is parallel
/// (sort, box unions) and deterministic.
#[derive(Debug, Clone)]
struct BoxTree {
    order: Vec<u32>,
    sboxes: Vec<Aabb<3>>,
    levels: Vec<Vec<Aabb<3>>>,
    leaf: usize,
}

fn split3(mut x: u64) -> u64 {
    x &= 0x1f_ffff;
    x = (x | x << 32) & 0x001f_0000_0000_ffff;
    x = (x | x << 16) & 0x001f_0000_ff00_00ff;
    x = (x | x << 8) & 0x100f_00f0_0f00_f00f;
    x = (x | x << 4) & 0x10c3_0c30_c30c_30c3;
    x = (x | x << 2) & 0x1249_2492_4924_9249;
    x
}

fn build_tree(boxes: &[Aabb<3>], leaf: usize) -> Result<Option<BoxTree>, RsError> {
    if boxes.is_empty() {
        return Ok(None);
    }
    let centers: Vec<[f64; 3]> = boxes.par_iter().map(center).collect();
    let (lo, hi) = centers
        .par_iter()
        .fold(
            || ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |(mut lo, mut hi), c| {
                for k in 0..3 {
                    lo[k] = lo[k].min(c[k]);
                    hi[k] = hi[k].max(c[k]);
                }
                (lo, hi)
            },
        )
        .reduce(
            || ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |a, b| {
                let mut lo = a.0;
                let mut hi = a.1;
                for k in 0..3 {
                    lo[k] = lo[k].min(b.0[k]);
                    hi[k] = hi[k].max(b.1[k]);
                }
                (lo, hi)
            },
        );
    let q = |c: &[f64; 3], k: usize| -> u64 {
        let w = hi[k] - lo[k];
        if w > 0.0 {
            (((c[k] - lo[k]) / w) * 2097151.0) as u64
        } else {
            0
        }
    };
    let mut keys: Vec<(u64, u32)> = centers
        .par_iter()
        .enumerate()
        .map(|(i, c)| {
            (
                split3(q(c, 0)) | split3(q(c, 1)) << 1 | split3(q(c, 2)) << 2,
                i as u32,
            )
        })
        .collect();
    keys.par_sort_unstable();
    let order: Vec<u32> = keys.iter().map(|k| k.1).collect();
    let sboxes: Vec<Aabb<3>> = order.par_iter().map(|&i| boxes[i as usize]).collect();
    let mut levels = vec![
        sboxes
            .par_chunks(leaf)
            .map(|c| c.iter().fold(Aabb::<3>::empty(), |a, b| a.union(b)))
            .collect::<Vec<_>>(),
    ];
    while levels.last().expect("a level").len() > 1 {
        let next: Vec<Aabb<3>> = levels
            .last()
            .expect("a level")
            .par_chunks(2)
            .map(|c| c.iter().fold(Aabb::<3>::empty(), |a, b| a.union(b)))
            .collect();
        levels.push(next);
    }
    Ok(Some(BoxTree {
        order,
        sboxes,
        levels,
        leaf,
    }))
}

/// Calls `f(e)` for every element whose box overlaps `q`.
fn query(tree: &BoxTree, q: &Aabb<3>, mut f: impl FnMut(u32)) {
    let mut stack = [(0u8, 0u32); 96];
    stack[0] = ((tree.levels.len() - 1) as u8, 0);
    let mut sp = 1;
    while sp > 0 {
        sp -= 1;
        let (k, i) = stack[sp];
        if !overlaps(&tree.levels[k as usize][i as usize], q) {
            continue;
        }
        if k == 0 {
            let a = i as usize * tree.leaf;
            let b = (a + tree.leaf).min(tree.order.len());
            for pos in a..b {
                if overlaps(&tree.sboxes[pos], q) {
                    f(tree.order[pos]);
                }
            }
        } else {
            let below = &tree.levels[k as usize - 1];
            stack[sp] = (k - 1, 2 * i);
            sp += 1;
            if ((2 * i + 1) as usize) < below.len() {
                stack[sp] = (k - 1, 2 * i + 1);
                sp += 1;
            }
        }
    }
}

/// A static triangle obstacle with its broad-phase trees (built once).
#[derive(Debug, Clone)]
pub struct CcdObstacle {
    vertices: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    edges: Vec<[u32; 2]>,
    tri_tree: BoxTree,
    edge_tree: BoxTree,
    vert_tree: BoxTree,
}

impl CcdObstacle {
    /// Builds the trees over the triangles, the (unique) edges and the vertices.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for `leaf_size = 0`, [`RsError::Curves`] if a tree fails.
    pub fn new(obstacle: &Obstacle, leaf_size: usize) -> Result<CcdObstacle, RsError> {
        if leaf_size == 0 {
            return Err(bad("ccd.leaf_size", "0".into(), "at least 1"));
        }
        let vertices = obstacle.vertices().to_vec();
        let triangles = obstacle.triangle_indices().to_vec();
        let mut edges: Vec<[u32; 2]> = Vec::new();
        for t in &triangles {
            for s in 0..3 {
                let (a, b) = (t[s], t[(s + 1) % 3]);
                edges.push([a.min(b), a.max(b)]);
            }
        }
        edges.sort_unstable();
        edges.dedup();
        let tb: Vec<Aabb<3>> = triangles
            .iter()
            .map(|t| hull(&[vertices[t[0] as usize], vertices[t[1] as usize], vertices[t[2] as usize]]))
            .collect();
        let eb: Vec<Aabb<3>> = edges
            .iter()
            .map(|e| hull(&[vertices[e[0] as usize], vertices[e[1] as usize]]))
            .collect();
        let vb: Vec<Aabb<3>> = vertices.iter().map(|&v| Aabb::from_point(v)).collect();
        Ok(CcdObstacle {
            tri_tree: build_tree(&tb, leaf_size)?.ok_or_else(|| {
                bad("obstacle", "no triangles".into(), "at least one triangle")
            })?,
            edge_tree: build_tree(&eb, leaf_size)?.expect("edges exist with triangles"),
            vert_tree: build_tree(&vb, leaf_size)?.expect("vertices exist"),
            vertices,
            triangles,
            edges,
        })
    }

    /// The vertices.
    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.vertices
    }

    /// The triangles.
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// The unique edges `[a, b]`, `a < b`, sorted.
    pub fn edges(&self) -> &[[u32; 2]] {
        &self.edges
    }

    fn tri(&self, t: u32) -> [[f64; 3]; 3] {
        self.triangles[t as usize].map(|i| self.vertices[i as usize])
    }
}

fn hull(ps: &[[f64; 3]]) -> Aabb<3> {
    let mut b = Aabb::<3>::empty();
    for &p in ps {
        b.grow(p);
    }
    b
}

// ---------------------------------------------------------------------------------------------
// Conservative advancement
// ---------------------------------------------------------------------------------------------

#[inline]
fn at(p: [f64; 3], u: [f64; 3], t: f64) -> [f64; 3] {
    [p[0] + t * u[0], p[1] + t * u[1], p[2] + t * u[2]]
}

/// The Lipschitz constant of the distance of two primitives with corner velocities `ua`, `ub`:
/// `max |u_a − c| + max |u_b − c|` for the better of `c = 0` and `c` = the mean corner velocity.
fn lipschitz(ua: &[[f64; 3]], ub: &[[f64; 3]]) -> f64 {
    let n = (ua.len() + ub.len()) as f64;
    let mut mean = [0.0; 3];
    for v in ua.iter().chain(ub) {
        mean = add(mean, *v);
    }
    mean = scale(mean, 1.0 / n);
    let l_of = |c: [f64; 3]| {
        let la = ua.iter().map(|v| norm(sub(*v, c))).fold(0.0f64, f64::max);
        let lb = ub.iter().map(|v| norm(sub(*v, c))).fold(0.0f64, f64::max);
        la + lb
    };
    // A relative margin against rounding in the speeds.
    l_of([0.0; 3]).min(l_of(mean)) * (1.0 + 1e-12)
}

/// Relative margin against rounding in the speeds.
const MARGIN: f64 = 1.0 + 1e-12;

/// Counters of the narrow phase.
#[derive(Debug, Clone, Copy, Default)]
struct Tally {
    iterations: u64,
    pairs: u64,
    close: u64,
    touching: u64,
    exempt: u64,
}

impl Tally {
    fn merge(mut self, o: Tally) -> Tally {
        self.iterations += o.iterations;
        self.pairs += o.pairs;
        self.close += o.close;
        self.touching += o.touching;
        self.exempt += o.exempt;
        self
    }
}

/// Conservative advancement of one pair on `[lo, hi]` with distance function `dist(t)` (valid
/// from `t = 0`), Lipschitz constant `l` (`SPEC.md` §7.3). Returns the earliest safe time at
/// which the distance has reached the pair's target, if before `hi`.
///
/// The target is `gap_abs`, or `close_factor · d₀` for a pair with `d₀ < gap_abs`; a pair with
/// `d₀ = 0` is unconstrained.
///
/// [formula: rs-ccd-advancement]
#[allow(clippy::too_many_arguments)]
fn advance(
    dist: impl Fn(f64) -> f64,
    l_cheap: f64,
    l_refined: impl Fn() -> f64,
    gap_abs: f64,
    opts: &CcdOptions,
    lo: f64,
    hi: f64,
    tally: &mut Tally,
) -> Option<f64> {
    let d0 = dist(0.0);
    tally.pairs += 1;
    tally.iterations += 1;
    if lo == 0.0 {
        if d0 <= 1e-12 * gap_abs {
            tally.touching += 1;
            return None;
        }
        if d0 < gap_abs {
            tally.close += 1;
        }
    } else if d0 <= 1e-12 * gap_abs {
        return None;
    }
    let target = if d0 > gap_abs {
        gap_abs
    } else {
        opts.close_factor * d0
    };
    let slack = 1e-9 * d0;
    let mut l = l_cheap;
    let mut refined = false;
    let mut tau = lo;
    let mut d = if lo == 0.0 { d0 } else { dist(lo) };
    if lo != 0.0 {
        tally.iterations += 1;
    }
    for _ in 0..opts.max_iterations {
        if d.is_nan() || d - target <= slack {
            return Some(tau);
        }
        let mut next = tau + (d - target) / l;
        if next >= hi && !refined {
            // The cheap bound (largest speeds, no common frame) did not clear the pair: try
            // the better one.
            l = l.min(l_refined());
            refined = true;
            next = tau + (d - target) / l;
        }
        if next >= hi || next.is_nan() {
            return None;
        }
        tau = next;
        d = dist(tau);
        tally.iterations += 1;
    }
    Some(tau)
}

/// Conservative advancement of a moving vertex against a moving triangle on `[t_start, t_end]`
/// (the gap `target` is absolute; positions are `p + t u`). Returns the earliest time at which
/// the distance is within `rounding` of `target`, and the number of distance evaluations. For
/// tests and other callers; [`first_contact`] uses the same step with per-pair targets.
///
/// [formula: rs-ccd-advancement]
#[allow(clippy::too_many_arguments)]
pub fn vertex_triangle_advance(
    p: [f64; 3],
    u: [f64; 3],
    tri: [[f64; 3]; 3],
    w: [[f64; 3]; 3],
    target: f64,
    t_start: f64,
    t_end: f64,
    max_iterations: u32,
) -> (Option<f64>, u64) {
    let l = lipschitz(&[u], &w);
    let dist = |t: f64| {
        point_triangle_distance(
            at(p, u, t),
            at(tri[0], w[0], t),
            at(tri[1], w[1], t),
            at(tri[2], w[2], t),
        )
    };
    raw_advance(dist, l, target, t_start, t_end, max_iterations)
}

/// As [`vertex_triangle_advance`] for two moving segments.
///
/// [formula: rs-ccd-advancement]
#[allow(clippy::too_many_arguments)]
pub fn edge_edge_advance(
    p: [[f64; 3]; 2],
    u: [[f64; 3]; 2],
    q: [[f64; 3]; 2],
    w: [[f64; 3]; 2],
    target: f64,
    t_start: f64,
    t_end: f64,
    max_iterations: u32,
) -> (Option<f64>, u64) {
    let l = lipschitz(&u, &w);
    let dist = |t: f64| {
        segment_segment_distance(
            at(p[0], u[0], t),
            at(p[1], u[1], t),
            at(q[0], w[0], t),
            at(q[1], w[1], t),
        )
    };
    raw_advance(dist, l, target, t_start, t_end, max_iterations)
}

fn raw_advance(
    dist: impl Fn(f64) -> f64,
    l: f64,
    target: f64,
    t_start: f64,
    t_end: f64,
    max_iterations: u32,
) -> (Option<f64>, u64) {
    let mut tau = t_start;
    let mut d = dist(tau);
    let mut evals = 1;
    for _ in 0..max_iterations {
        if d.is_nan() || d - target <= 1e-9 * target.max(f64::MIN_POSITIVE) {
            return (Some(tau), evals);
        }
        if l <= 0.0 {
            return (None, evals);
        }
        let next = tau + (d - target) / l;
        if next >= t_end {
            return (None, evals);
        }
        tau = next;
        d = dist(tau);
        evals += 1;
    }
    (Some(tau), evals)
}

// ---------------------------------------------------------------------------------------------
// First contact
// ---------------------------------------------------------------------------------------------

/// What one query (a vertex, an edge or a triangle) found in one window.
#[derive(Debug, Clone, Copy, Default)]
struct Found {
    best: Option<Contact>,
    tally: Tally,
}

fn earlier(a: Option<Contact>, b: Option<Contact>) -> Option<Contact> {
    match (a, b) {
        (None, x) | (x, None) => x,
        (Some(x), Some(y)) => {
            let o = x.tau.total_cmp(&y.tau).then(x.primitive.cmp(&y.primitive));
            Some(if o.is_le() { x } else { y })
        }
    }
}

impl Found {
    fn merge(self, o: Found) -> Found {
        Found {
            best: earlier(self.best, o.best),
            tally: self.tally.merge(o.tally),
        }
    }

    /// Records a candidate time of `primitive`, keeping the earliest.
    fn offer(&mut self, tau: Option<f64>, primitive: Primitive) {
        if let Some(tau) = tau {
            self.best = earlier(self.best, Some(Contact { tau, primitive }));
        }
    }

    /// The window end for the next pair of this query: nothing later than the best so far
    /// can change the result.
    fn horizon(&self, hi: f64) -> f64 {
        self.best.map_or(hi, |c| c.tau.min(hi))
    }
}

struct Motion<'a> {
    surface: &'a Surface,
    start: &'a [[f64; 3]],
    vel: &'a [[f64; 3]],
    speed: &'a [f64],
    gap_abs: f64,
    opts: &'a CcdOptions,
    band: Option<&'a [bool]>,
}

impl Motion<'_> {
    #[inline]
    fn pos(&self, v: u32, t: f64) -> [f64; 3] {
        at(self.start[v as usize], self.vel[v as usize], t)
    }
    #[inline]
    fn tri_at(&self, f: u32, t: f64) -> [[f64; 3]; 3] {
        self.surface.faces()[f as usize].map(|v| self.pos(v, t))
    }
    /// Is any of `vs` in the exempt band?
    #[inline]
    fn in_band(&self, vs: &[u32]) -> bool {
        self.band.is_some_and(|b| vs.iter().any(|&v| b[v as usize]))
    }
    /// The largest corner speed of face `f` (with the rounding margin).
    #[inline]
    fn tri_speed(&self, f: u32) -> f64 {
        let t = self.surface.faces()[f as usize];
        self.speed[t[0] as usize]
            .max(self.speed[t[1] as usize])
            .max(self.speed[t[2] as usize])
    }
    fn tri_vel(&self, f: u32) -> [[f64; 3]; 3] {
        self.surface.faces()[f as usize].map(|v| self.vel[v as usize])
    }
}

/// The largest `τ_max ≤ tau_end` such that `start + t · velocity` keeps all pairs of non-adjacent
/// elements of `surface` (and of the mesh and the static `obstacles`) at distance `≥ gap` for
/// `t ∈ [0, τ_max]`. See the module documentation. `pool` is the thread pool of the session.
///
/// [formula: rs-ccd-advancement]
///
/// # Errors
/// The errors of [`CcdOptions::validate`] and [`Surface::validate_positions`] (for `start`);
/// [`RsError::InvalidOption`] for a velocity of the wrong length or not finite, or a `tau_end`
/// negative or not finite.
pub fn first_contact(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    start: &[[f64; 3]],
    velocity: &[[f64; 3]],
    tau_end: f64,
    obstacles: &[CcdObstacle],
    opts: &CcdOptions,
) -> Result<CcdReport, RsError> {
    first_contact_with_band(pool, surface, start, velocity, tau_end, obstacles, opts, None)
}

/// [`first_contact`] with an exempt band: a mesh-mesh pair whose two elements both have a vertex
/// in `band` is unconstrained (see the module documentation). `None` (or an empty band) is
/// bitwise [`first_contact`].
///
/// # Errors
/// As [`first_contact`]; [`RsError::InvalidOption`] for a band of another surface.
#[allow(clippy::too_many_arguments)]
pub fn first_contact_with_band(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    start: &[[f64; 3]],
    velocity: &[[f64; 3]],
    tau_end: f64,
    obstacles: &[CcdObstacle],
    opts: &CcdOptions,
    band: Option<&IntersectionBand>,
) -> Result<CcdReport, RsError> {
    let band = match band {
        Some(b) if b.is_empty() => None,
        Some(b) => {
            if b.in_band.len() != surface.n_vertices() {
                return Err(bad("band", format!("{} flags", b.in_band.len()), "one flag per vertex"));
            }
            Some(&b.in_band[..])
        }
        None => None,
    };
    let t_begin = Instant::now();
    opts.validate()?;
    surface.validate_positions(start)?;
    if velocity.len() != surface.n_vertices() {
        return Err(bad(
            "velocity",
            format!("{} vectors", velocity.len()),
            "one vector per vertex",
        ));
    }
    if !velocity.iter().flatten().all(|c| c.is_finite()) {
        return Err(bad("velocity", "non-finite".into(), "finite components"));
    }
    if !(tau_end.is_finite() && tau_end >= 0.0) {
        return Err(bad("tau_end", tau_end.to_string(), "finite and non-negative"));
    }
    let ell = opts
        .reference_length
        .unwrap_or_else(|| mean_edge_length(surface, start));
    let gap_abs = opts.gap * ell;
    let speed: Vec<f64> = velocity.iter().map(|&v| norm(v)).collect();
    let motion = Motion {
        surface,
        start,
        vel: velocity,
        speed: &speed,
        gap_abs,
        opts,
        band,
    };
    let vmax = velocity.iter().map(|&v| norm(v)).fold(0.0f64, f64::max);
    let mut report = CcdReport {
        tau_max: tau_end,
        contact: None,
        gap_abs,
        iterations: 0,
        pairs_tested: 0,
        windows: 0,
        close_pairs: 0,
        touching_pairs: 0,
        exempt_pairs: 0,
        seconds: 0.0,
    };
    // Time in which the fastest vertex moves one mean edge length.
    let t1 = if vmax > 0.0 {
        mean_edge_length(surface, start) / vmax
    } else {
        f64::INFINITY
    };
    // Even without motion (vmax = 0) window 0 runs once so the close/touching counts are filled.
    let mut lo = 0.0;
    let mut hi = t1.min(tau_end);
    loop {
        if vmax == 0.0 {
            hi = 0.0;
        }
        report.windows += 1;
        let found = window(pool, &motion, obstacles, lo, hi)?;
        report.iterations += found.tally.iterations;
        report.pairs_tested += found.tally.pairs;
        if lo == 0.0 {
            report.close_pairs = found.tally.close;
            report.touching_pairs = found.tally.touching;
            report.exempt_pairs = found.tally.exempt;
        }
        if let Some(c) = found.best {
            report.tau_max = c.tau;
            report.contact = Some(c);
            break;
        }
        if hi >= tau_end || vmax == 0.0 {
            break;
        }
        lo = hi;
        hi = (2.0 * hi).max(hi + t1).min(tau_end);
        // Window k covers [2^{k-1} t1, 2^k t1]: lo doubles each time.
        if hi <= lo {
            break;
        }
    }
    report.seconds = t_begin.elapsed().as_secs_f64();
    Ok(report)
}

/// All pairs on the window `[lo, hi]`.
fn window(
    pool: &rayon::ThreadPool,
    m: &Motion<'_>,
    obstacles: &[CcdObstacle],
    lo: f64,
    hi: f64,
) -> Result<Found, RsError> {
    let s = m.surface;
    let faces = s.faces();
    let edges = s.edges();
    let h = 0.5 * m.gap_abs * (1.0 + 1e-6);
    let leaf = m.opts.leaf_size;
    // Swept boxes (uninflated hulls of the corner positions at both ends of the window).
    let vbox: Vec<Aabb<3>> = pool.install(|| {
        (0..s.n_vertices() as u32)
            .into_par_iter()
            .map(|v| hull(&[m.pos(v, lo), m.pos(v, hi)]))
            .collect()
    });
    let union3 = |a: &Aabb<3>, b: &Aabb<3>, c: &Aabb<3>| a.union(b).union(c);
    let fbox: Vec<Aabb<3>> = pool.install(|| {
        faces
            .par_iter()
            .map(|f| union3(&vbox[f[0] as usize], &vbox[f[1] as usize], &vbox[f[2] as usize]))
            .collect()
    });
    let ebox: Vec<Aabb<3>> = pool.install(|| {
        edges
            .par_iter()
            .map(|e| vbox[e[0] as usize].union(&vbox[e[1] as usize]))
            .collect()
    });
    let (ftree, etree) = pool.install(|| {
        rayon::join(
            || {
                let b: Vec<Aabb<3>> = fbox.iter().map(|b| inflate(b, h)).collect();
                build_tree(&b, leaf)
            },
            || {
                let b: Vec<Aabb<3>> = ebox.iter().map(|b| inflate(b, h)).collect();
                build_tree(&b, leaf)
            },
        )
    });
    let (ftree, etree) = (ftree?, etree?);
    let (ftree, etree) = (
        ftree.expect("a surface has faces"),
        etree.expect("a surface has edges"),
    );

    let mut total = Found::default();

    // Vertex against mesh triangle.
    let vt = pool.install(|| {
        (0..s.n_vertices() as u32)
            .into_par_iter()
            .map(|v| {
                let mut fnd = Found::default();
                let q = inflate(&vbox[v as usize], h);
                query(&ftree, &q, |f| {
                    let tri = faces[f as usize];
                    if tri.contains(&v) {
                        return;
                    }
                    if m.in_band(&[v]) && m.in_band(&tri) {
                        fnd.tally.exempt += 1;
                        return;
                    }
                    let (p, u) = (m.start[v as usize], m.vel[v as usize]);
                    let (tp, tu) = (m.tri_at(f, 0.0), m.tri_vel(f));
                    let lc = (m.speed[v as usize] + m.tri_speed(f)) * MARGIN;
                    let dist = |t: f64| {
                        let c = [at(tp[0], tu[0], t), at(tp[1], tu[1], t), at(tp[2], tu[2], t)];
                        point_triangle_distance(at(p, u, t), c[0], c[1], c[2])
                    };
                    let hz = fnd.horizon(hi);
                    let r = advance(
                        dist,
                        lc,
                        || lipschitz(&[u], &tu),
                        m.gap_abs,
                        m.opts,
                        lo,
                        hz,
                        &mut fnd.tally,
                    );
                    fnd.offer(
                        r,
                        Primitive::VertexTriangle {
                            vertex: v,
                            triangle: f,
                        },
                    );
                });
                fnd
            })
            .reduce(Found::default, Found::merge)
    });
    total = total.merge(vt);

    // Edge against mesh edge.
    let ee = pool.install(|| {
        (0..s.n_edges() as u32)
            .into_par_iter()
            .map(|e| {
                let mut fnd = Found::default();
                let [a, b] = edges[e as usize];
                let q = inflate(&ebox[e as usize], h);
                query(&etree, &q, |e2| {
                    if e2 <= e {
                        return;
                    }
                    let [c, d] = edges[e2 as usize];
                    if c == a || c == b || d == a || d == b {
                        return;
                    }
                    if m.in_band(&[a, b]) && m.in_band(&[c, d]) {
                        fnd.tally.exempt += 1;
                        return;
                    }
                    let pp = [m.start[a as usize], m.start[b as usize]];
                    let uu = [m.vel[a as usize], m.vel[b as usize]];
                    let qq = [m.start[c as usize], m.start[d as usize]];
                    let ww = [m.vel[c as usize], m.vel[d as usize]];
                    let lc = (m.speed[a as usize].max(m.speed[b as usize])
                        + m.speed[c as usize].max(m.speed[d as usize]))
                        * MARGIN;
                    let dist = |t: f64| {
                        segment_segment_distance(
                            at(pp[0], uu[0], t),
                            at(pp[1], uu[1], t),
                            at(qq[0], ww[0], t),
                            at(qq[1], ww[1], t),
                        )
                    };
                    let hz = fnd.horizon(hi);
                    let r = advance(
                        dist,
                        lc,
                        || lipschitz(&uu, &ww),
                        m.gap_abs,
                        m.opts,
                        lo,
                        hz,
                        &mut fnd.tally,
                    );
                    fnd.offer(r, Primitive::EdgeEdge { a: e, b: e2 });
                });
                fnd
            })
            .reduce(Found::default, Found::merge)
    });
    total = total.merge(ee);

    for (oi, ob) in obstacles.iter().enumerate() {
        total = total.merge(obstacle_window(pool, m, oi, ob, &vbox, &fbox, &ebox, lo, hi));
    }
    Ok(total)
}

#[allow(clippy::too_many_arguments)]
fn obstacle_window(
    pool: &rayon::ThreadPool,
    m: &Motion<'_>,
    oi: usize,
    ob: &CcdObstacle,
    vbox: &[Aabb<3>],
    fbox: &[Aabb<3>],
    ebox: &[Aabb<3>],
    lo: f64,
    hi: f64,
) -> Found {
    let s = m.surface;
    let g = m.gap_abs * (1.0 + 1e-6);
    // Mesh vertex against obstacle triangle.
    let a = pool.install(|| {
        (0..s.n_vertices() as u32)
            .into_par_iter()
            .map(|v| {
                let mut fnd = Found::default();
                query(&ob.tri_tree, &inflate(&vbox[v as usize], g), |t| {
                    let tri = ob.tri(t);
                    let (p, u) = (m.start[v as usize], m.vel[v as usize]);
                    let l = m.speed[v as usize] * MARGIN;
                    let dist = |tm: f64| point_triangle_distance(at(p, u, tm), tri[0], tri[1], tri[2]);
                    let hz = fnd.horizon(hi);
                    let r = advance(dist, l, || l, m.gap_abs, m.opts, lo, hz, &mut fnd.tally);
                    fnd.offer(
                        r,
                        Primitive::VertexObstacleTriangle {
                            vertex: v,
                            obstacle: oi,
                            triangle: t,
                        },
                    );
                });
                fnd
            })
            .reduce(Found::default, Found::merge)
    });
    // Mesh triangle against obstacle vertex.
    let b = pool.install(|| {
        (0..s.n_faces() as u32)
            .into_par_iter()
            .map(|f| {
                let mut fnd = Found::default();
                query(&ob.vert_tree, &inflate(&fbox[f as usize], g), |ov| {
                    let p = ob.vertices[ov as usize];
                    let (tp, tu) = (m.tri_at(f, 0.0), m.tri_vel(f));
                    let l = m.tri_speed(f) * MARGIN;
                    let dist = |t: f64| {
                        let c = [at(tp[0], tu[0], t), at(tp[1], tu[1], t), at(tp[2], tu[2], t)];
                        point_triangle_distance(p, c[0], c[1], c[2])
                    };
                    let hz = fnd.horizon(hi);
                    let r = advance(dist, l, || l, m.gap_abs, m.opts, lo, hz, &mut fnd.tally);
                    fnd.offer(
                        r,
                        Primitive::TriangleObstacleVertex {
                            triangle: f,
                            obstacle: oi,
                            vertex: ov,
                        },
                    );
                });
                fnd
            })
            .reduce(Found::default, Found::merge)
    });
    // Mesh edge against obstacle edge.
    let c = pool.install(|| {
        (0..s.n_edges() as u32)
            .into_par_iter()
            .map(|e| {
                let mut fnd = Found::default();
                let [ea, eb] = s.edges()[e as usize];
                query(&ob.edge_tree, &inflate(&ebox[e as usize], g), |oe| {
                    let [oa, obb] = ob.edges[oe as usize];
                    let q = [ob.vertices[oa as usize], ob.vertices[obb as usize]];
                    let pp = [m.start[ea as usize], m.start[eb as usize]];
                    let uu = [m.vel[ea as usize], m.vel[eb as usize]];
                    let l = m.speed[ea as usize].max(m.speed[eb as usize]) * MARGIN;
                    let dist = |t: f64| {
                        segment_segment_distance(at(pp[0], uu[0], t), at(pp[1], uu[1], t), q[0], q[1])
                    };
                    let hz = fnd.horizon(hi);
                    let r = advance(dist, l, || l, m.gap_abs, m.opts, lo, hz, &mut fnd.tally);
                    fnd.offer(
                        r,
                        Primitive::EdgeObstacleEdge {
                            edge: e,
                            obstacle: oi,
                            obstacle_edge: oe,
                        },
                    );
                });
                fnd
            })
            .reduce(Found::default, Found::merge)
    });
    a.merge(b).merge(c)
}

// ---------------------------------------------------------------------------------------------
// Static queries
// ---------------------------------------------------------------------------------------------

fn orient(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    dot(cross(sub(b, a), sub(c, a)), sub(d, a))
}

/// Does the closed segment `p q` meet the closed triangle `a b c`?
fn segment_meets_triangle(p: [f64; 3], q: [f64; 3], tri: [[f64; 3]; 3]) -> bool {
    let [a, b, c] = tri;
    let d1 = orient(a, b, c, p);
    let d2 = orient(a, b, c, q);
    if (d1 > 0.0 && d2 > 0.0) || (d1 < 0.0 && d2 < 0.0) {
        return false;
    }
    // Coplanar up to rounding: both endpoints within 1e-12 of the longest length involved from
    // the plane (`|orient| = |n| · distance`). Exact zeros are not enough: points in the plane
    // of a generic triangle give `orient` = ±1e-20, and the sign test below is then arbitrary
    // (it reported coplanar neighbours with one common vertex as intersecting).
    let n = norm(cross(sub(b, a), sub(c, a)));
    let len = norm(sub(q, p)).max(norm(sub(b, a))).max(norm(sub(c, a))).max(norm(sub(c, b)));
    let flat = 1e-12 * n * len;
    if d1.abs() <= flat && d2.abs() <= flat {
        // Coplanar: touching up to rounding.
        let eps = 1e-14 * (norm(sub(b, a)) + norm(sub(c, a)));
        return point_triangle_distance(p, a, b, c) <= eps
            || point_triangle_distance(q, a, b, c) <= eps
            || segment_segment_distance(p, q, a, b) <= eps
            || segment_segment_distance(p, q, b, c) <= eps
            || segment_segment_distance(p, q, c, a) <= eps;
    }
    let s1 = orient(p, q, a, b);
    let s2 = orient(p, q, b, c);
    let s3 = orient(p, q, c, a);
    (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) || (s1 <= 0.0 && s2 <= 0.0 && s3 <= 0.0)
}

/// Do two triangles with no common vertex intersect (closed)?
pub fn triangles_intersect(t1: [[f64; 3]; 3], t2: [[f64; 3]; 3]) -> bool {
    (0..3).any(|k| segment_meets_triangle(t1[k], t1[(k + 1) % 3], t2))
        || (0..3).any(|k| segment_meets_triangle(t2[k], t2[(k + 1) % 3], t1))
}

/// The distance of two disjoint-vertex triangles: `0` if they intersect, else the minimum over
/// the vertex–triangle and edge–edge distances.
pub fn triangle_distance(t1: [[f64; 3]; 3], t2: [[f64; 3]; 3]) -> f64 {
    if triangles_intersect(t1, t2) {
        return 0.0;
    }
    let mut d = f64::INFINITY;
    for k in 0..3 {
        d = d.min(point_triangle_distance(t1[k], t2[0], t2[1], t2[2]));
        d = d.min(point_triangle_distance(t2[k], t1[0], t1[1], t1[2]));
        for j in 0..3 {
            d = d.min(segment_segment_distance(
                t1[k],
                t1[(k + 1) % 3],
                t2[j],
                t2[(j + 1) % 3],
            ));
        }
    }
    d
}

fn face_boxes(surface: &Surface, positions: &[[f64; 3]]) -> Vec<Aabb<3>> {
    surface
        .faces()
        .iter()
        .map(|f| hull(&[positions[f[0] as usize], positions[f[1] as usize], positions[f[2] as usize]]))
        .collect()
}

/// The pairs of faces `[f, g]`, `f < g`, that intersect: faces without a common vertex that
/// touch each other, and faces with exactly one common vertex whose opposite edge meets the
/// other face (faces that share an edge are not tested). Sorted.
///
/// # Errors
/// The errors of [`Surface::validate_positions`].
pub fn intersecting_pairs(
    surface: &Surface,
    positions: &[[f64; 3]],
) -> Result<Vec<[u32; 2]>, RsError> {
    surface.validate_positions(positions)?;
    let faces = surface.faces();
    let boxes = face_boxes(surface, positions);
    let tree = build_tree(&boxes, 8)?.expect("a surface has faces");
    let tri = |f: u32| faces[f as usize].map(|v| positions[v as usize]);
    let mut out: Vec<[u32; 2]> = (0..faces.len() as u32)
        .into_par_iter()
        .flat_map_iter(|f| {
            let mut hits = Vec::new();
            query(&tree, &boxes[f as usize], |g| {
                if g <= f {
                    return;
                }
                let (a, b) = (faces[f as usize], faces[g as usize]);
                let mut common = [0usize; 3];
                let mut n_common = 0;
                for i in 0..3 {
                    if b.contains(&a[i]) {
                        common[n_common] = i;
                        n_common += 1;
                    }
                }
                let hit = match n_common {
                    0 => triangles_intersect(tri(f), tri(g)),
                    1 => {
                        let i = common[0];
                        let (p, q) = (positions[a[(i + 1) % 3] as usize], positions[a[(i + 2) % 3] as usize]);
                        let j = (0..3).find(|&j| b[j] == a[i]).expect("common vertex");
                        let (r, s) = (positions[b[(j + 1) % 3] as usize], positions[b[(j + 2) % 3] as usize]);
                        segment_meets_triangle(p, q, tri(g)) || segment_meets_triangle(r, s, tri(f))
                    }
                    _ => false,
                };
                if hit {
                    hits.push([f, g]);
                }
            });
            hits
        })
        .collect();
    out.sort_unstable();
    Ok(out)
}

/// The vertices `I` that are exempt from the gap rule because the mesh intersects itself near
/// them: the vertices of the faces in [`intersecting_pairs`], widened by one ring (every vertex
/// that shares an edge with one of them). Empty for an embedded mesh.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntersectionBand {
    in_band: Vec<bool>,
    n_pairs: usize,
    n_vertices: usize,
}

impl IntersectionBand {
    /// The band of `positions`; the intersection test runs on `pool`.
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`].
    pub fn new(
        pool: &rayon::ThreadPool,
        surface: &Surface,
        positions: &[[f64; 3]],
    ) -> Result<IntersectionBand, RsError> {
        let pairs = pool.install(|| intersecting_pairs(surface, positions))?;
        let nv = surface.n_vertices();
        let mut core = vec![false; nv];
        for &[f, g] in &pairs {
            for &fid in &[f, g] {
                for &v in &surface.faces()[fid as usize] {
                    core[v as usize] = true;
                }
            }
        }
        let mut in_band = core.clone();
        if !pairs.is_empty() {
            for &[a, b] in surface.edges() {
                if core[a as usize] || core[b as usize] {
                    in_band[a as usize] = true;
                    in_band[b as usize] = true;
                }
            }
        }
        let n_vertices = in_band.iter().filter(|&&x| x).count();
        Ok(IntersectionBand { in_band, n_pairs: pairs.len(), n_vertices })
    }

    /// No vertex in the band (the mesh is embedded).
    pub fn is_empty(&self) -> bool {
        self.n_vertices == 0
    }

    /// Number of intersecting face pairs at the start.
    pub fn intersecting_pairs(&self) -> usize {
        self.n_pairs
    }

    /// Number of vertices in the band.
    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    /// Per-vertex flags.
    pub fn flags(&self) -> &[bool] {
        &self.in_band
    }
}

/// Does the mesh intersect itself? (see [`intersecting_pairs`]).
///
/// # Errors
/// The errors of [`Surface::validate_positions`].
pub fn intersects(surface: &Surface, positions: &[[f64; 3]]) -> Result<bool, RsError> {
    Ok(!intersecting_pairs(surface, positions)?.is_empty())
}

/// The smallest distance over all pairs of non-adjacent elements (vertex–triangle with the
/// vertex not in the triangle, edge–edge without a common vertex); `+∞` if there are none.
///
/// # Errors
/// The errors of [`Surface::validate_positions`].
pub fn min_distance(surface: &Surface, positions: &[[f64; 3]]) -> Result<f64, RsError> {
    surface.validate_positions(positions)?;
    let faces = surface.faces();
    let edges = surface.edges();
    let fb = face_boxes(surface, positions);
    let eb: Vec<Aabb<3>> = edges
        .iter()
        .map(|e| hull(&[positions[e[0] as usize], positions[e[1] as usize]]))
        .collect();
    let ft = build_tree(&fb, 8)?.expect("faces");
    let et = build_tree(&eb, 8)?.expect("edges");
    let mut r = mean_edge_length(surface, positions);
    loop {
        let best = (0..surface.n_vertices() as u32)
            .into_par_iter()
            .map(|v| {
                let p = positions[v as usize];
                let mut best = f64::INFINITY;
                query(&ft, &inflate(&Aabb::from_point(p), r), |f| {
                    let t = faces[f as usize];
                    if !t.contains(&v) {
                        best = best.min(point_triangle_distance(
                            p,
                            positions[t[0] as usize],
                            positions[t[1] as usize],
                            positions[t[2] as usize],
                        ));
                    }
                });
                best
            })
            .reduce(|| f64::INFINITY, f64::min)
            .min(
                (0..edges.len() as u32)
                    .into_par_iter()
                    .map(|e| {
                        let [a, b] = edges[e as usize];
                        let mut best = f64::INFINITY;
                        query(&et, &inflate(&eb[e as usize], r), |e2| {
                            if e2 <= e {
                                return;
                            }
                            let [c, d] = edges[e2 as usize];
                            if c == a || c == b || d == a || d == b {
                                return;
                            }
                            best = best.min(segment_segment_distance(
                                positions[a as usize],
                                positions[b as usize],
                                positions[c as usize],
                                positions[d as usize],
                            ));
                        });
                        best
                    })
                    .reduce(|| f64::INFINITY, f64::min),
            );
        // Every pair within `r` was examined, so a value `≤ r` is the exact minimum.
        if best <= r {
            return Ok(best);
        }
        let all = fb
            .iter()
            .fold(Aabb::<3>::empty(), |acc, b| acc.union(b));
        let diam = norm(sub(all.hi, all.lo));
        if r > 2.0 * diam {
            return Ok(best);
        }
        r *= 2.0;
    }
}

// ---------------------------------------------------------------------------------------------
// Collapse guard
// ---------------------------------------------------------------------------------------------

type Cell = [i64; 3];

/// The CCD gap test of the remesher (`SPEC.md` §10): refuses a collapse whose new star
/// triangles intersect, or come within the gap of, a triangle that shares no vertex with them.
///
/// Use one guard per [`remesh`](crate::remesh::remesh) call: it builds a vertex grid on the
/// first query and keeps it up to date by watching which collapses were applied.
#[derive(Debug)]
pub struct CcdCollapseGuard {
    gap_abs: f64,
    grid: Option<HashMap<Cell, Vec<u32>>>,
    cell: f64,
    max_edge: f64,
    last: Option<(u32, u32, [f64; 3])>,
    refused: u64,
    tested: u64,
    band: Option<Vec<bool>>,
    exempt: u64,
}

impl CcdCollapseGuard {
    /// A guard with the absolute gap `opts.gap · l0` (`l0` = `opts.reference_length` if set).
    ///
    /// # Errors
    /// The errors of [`CcdOptions::validate`]; `l0` must be finite and positive.
    pub fn new(opts: &CcdOptions, l0: f64) -> Result<CcdCollapseGuard, RsError> {
        opts.validate()?;
        let l = opts.reference_length.unwrap_or(l0);
        if !(l.is_finite() && l > 0.0) {
            return Err(bad("l0", l.to_string(), "finite and positive"));
        }
        Ok(CcdCollapseGuard {
            gap_abs: opts.gap * l,
            grid: None,
            cell: 0.0,
            max_edge: 0.0,
            last: None,
            refused: 0,
            tested: 0,
            band: None,
            exempt: 0,
        })
    }

    /// Exempts triangle pairs that both have a vertex in `band` (indices of the surface passed
    /// to the remesher; vertices created by the remesher are not in the band, a collapse's
    /// survivor is when either end of its edge is). An empty band is ignored.
    pub fn set_band(&mut self, band: &IntersectionBand) {
        self.band = if band.is_empty() { None } else { Some(band.in_band.clone()) };
    }

    /// Number of triangle pairs skipped because of the band.
    pub fn exempt_pairs(&self) -> u64 {
        self.exempt
    }

    /// Number of collapses refused so far.
    pub fn refused(&self) -> u64 {
        self.refused
    }

    /// Number of triangle pairs examined so far.
    pub fn pairs_tested(&self) -> u64 {
        self.tested
    }

    fn cell_of(&self, p: [f64; 3]) -> Cell {
        [
            (p[0] / self.cell).floor() as i64,
            (p[1] / self.cell).floor() as i64,
            (p[2] / self.cell).floor() as i64,
        ]
    }

    fn build(&mut self, mesh: &crate::remesh::EditMesh) {
        let pos = mesh.positions();
        let mut max_edge = 0.0f64;
        for f in 0..mesh.n_face_slots() as u32 {
            if mesh.is_face_alive(f) {
                let t = mesh.face(f);
                for k in 0..3 {
                    max_edge = max_edge.max(norm(sub(pos[t[k] as usize], pos[t[(k + 1) % 3] as usize])));
                }
            }
        }
        self.max_edge = max_edge.max(self.gap_abs);
        self.cell = self.max_edge;
        let mut grid: HashMap<Cell, Vec<u32>> = HashMap::new();
        for v in 0..mesh.n_vertex_slots() as u32 {
            if mesh.is_vertex_alive(v) {
                grid.entry(self.cell_of(pos[v as usize])).or_default().push(v);
            }
        }
        self.grid = Some(grid);
    }
}

impl CollapseGuard for CcdCollapseGuard {
    fn allow(&mut self, q: &CollapseQuery<'_>) -> bool {
        let mesh = q.mesh;
        if self.grid.is_none() {
            self.build(mesh);
        }
        // Was the previous query applied? Then its survivor moved.
        if let Some((a, b, x)) = self.last.take() {
            let (survivor, removed) = if q.mesh.is_vertex_alive(a) { (a, b) } else { (b, a) };
            if !mesh.is_vertex_alive(removed) && mesh.is_vertex_alive(survivor) {
                let c = self.cell_of(x);
                self.grid.as_mut().expect("built").entry(c).or_default().push(survivor);
            }
        }
        let [ea, eb] = q.edge;
        self.last = Some((ea, eb, q.position));
        let pos = mesh.positions();
        let at_new = |v: u32| {
            if v == q.survivor {
                q.position
            } else {
                pos[v as usize]
            }
        };
        let mut star_max = 0.0f64;
        for t in q.new_star {
            for k in 0..3 {
                star_max = star_max.max(norm(sub(at_new(t[k]), at_new(t[(k + 1) % 3]))));
            }
        }
        self.max_edge = self.max_edge.max(star_max);
        let reach = self.gap_abs + self.max_edge;
        // The star's box and the cells it reaches.
        let mut bx = Aabb::<3>::empty();
        for t in q.new_star {
            for &v in t {
                bx.grow(at_new(v));
            }
        }
        let bx = inflate(&bx, reach);
        let (c0, c1) = (self.cell_of(bx.lo), self.cell_of(bx.hi));
        let grid = self.grid.as_ref().expect("built");
        let mut faces: Vec<u32> = Vec::new();
        for i in c0[0]..=c1[0] {
            for j in c0[1]..=c1[1] {
                for k in c0[2]..=c1[2] {
                    if let Some(vs) = grid.get(&[i, j, k]) {
                        for &w in vs {
                            if w != ea && w != eb && mesh.is_vertex_alive(w) && bx_contains(&bx, pos[w as usize]) {
                                faces.extend_from_slice(mesh.vertex_faces(w));
                            }
                        }
                    }
                }
            }
        }
        faces.sort_unstable();
        faces.dedup();
        let mut tested = 0;
        let mut exempt = 0;
        let mut refuse = false;
        let flag = |v: u32| self.band.as_ref().is_some_and(|b| b.get(v as usize).copied().unwrap_or(false));
        let star_in: Vec<bool> = q
            .new_star
            .iter()
            .map(|t| t.iter().any(|&v| flag(v) || (v == q.survivor && (flag(ea) || flag(eb)))))
            .collect();
        'outer: for &g in &faces {
            let gv = mesh.face(g);
            if gv.contains(&ea) || gv.contains(&eb) {
                continue;
            }
            let gt = gv.map(|v| pos[v as usize]);
            let gbox = inflate(&hull(&gt), self.gap_abs);
            let g_in = self.band.is_some() && gv.iter().any(|&v| flag(v));
            for (ti, t) in q.new_star.iter().enumerate() {
                if t.iter().any(|v| gv.contains(v)) {
                    continue;
                }
                if g_in && star_in[ti] {
                    exempt += 1;
                    continue;
                }
                let tt = t.map(at_new);
                if !overlaps(&gbox, &hull(&tt)) {
                    continue;
                }
                tested += 1;
                if triangle_distance(tt, gt) < self.gap_abs {
                    refuse = true;
                    break 'outer;
                }
            }
        }
        self.tested += tested;
        self.exempt += exempt;
        if refuse {
            self.refused += 1;
        }
        !refuse
    }
}

fn bx_contains(b: &Aabb<3>, p: [f64; 3]) -> bool {
    b.contains_point(p)
}
