//! Delaunay flipping (paper Sec. 4.1).
//!
//! An interior edge with opposite angles `α`, `β` is *Delaunay* when `α + β ≤ π`, which is the
//! same as a non-negative cotan weight `½ (cot α + cot β)`. The algorithm keeps a queue of
//! edges, flips the first one if it is not Delaunay, and puts its four neighbours back on the
//! queue. In exact arithmetic every edge that is not Delaunay can be flipped, the algorithm
//! terminates, and the result is unique up to the choice of diagonal in quadrilaterals whose
//! four vertices lie on a circle (`α + β = π`: a *tie*). The paper gives no bound on the number
//! of flips.
//!
//! In floating point arithmetic the test needs a tolerance: without one, rounding makes both
//! diagonals of a cocircular quadrilateral look non-Delaunay in turn. [`DelaunayCriterion`]
//! says which quantity the tolerance is applied to.

use std::collections::VecDeque;

use ddg_core::topology::NONE;

use crate::complex::{Complex, next, prev};
use crate::error::{FlipRefusal, IntrinsicError};
use crate::flip::{FlipRequest, FlipShapeTest};
use crate::formulas::{delaunay_angle_excess, opposite_cotan_sum};
use crate::interface::{Triangulation, opposite_cotan};
use crate::triangulation::IntrinsicTriangulation;

/// When an interior edge counts as not Delaunay, and is flipped.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DelaunayCriterion {
    /// The sum of the opposite angles exceeds `π` by more than `tolerance` radians (paper
    /// Eq. 4 with a tolerance). After flipping, every interior edge has
    /// `α + β ≤ π + tolerance`, hence a cotan weight of at least
    /// `−tolerance / (2 sin α sin β)` to first order.
    AngleSum {
        /// The tolerance in radians; finite and not negative.
        tolerance: f64,
    },
    /// The cotan weight `½ (cot α + cot β)` is below `−tolerance`. This is the test of the
    /// authors' implementation, which uses `1e-6`. After flipping, every interior edge has a
    /// cotan weight of at least `−tolerance`; in terms of angles the tolerance is
    /// `2 · tolerance · sin α sin β`, tighter the thinner the triangles.
    CotanWeight {
        /// The tolerance on the weight (dimensionless); finite and not negative.
        tolerance: f64,
    },
}

impl DelaunayCriterion {
    /// The test of the authors' implementation: [`DelaunayCriterion::CotanWeight`] with
    /// tolerance `1e-6`.
    pub const REFERENCE: DelaunayCriterion = DelaunayCriterion::CotanWeight { tolerance: 1e-6 };

    fn tolerance(&self) -> f64 {
        match *self {
            DelaunayCriterion::AngleSum { tolerance }
            | DelaunayCriterion::CotanWeight { tolerance } => tolerance,
        }
    }
}

impl Default for DelaunayCriterion {
    /// [`DelaunayCriterion::AngleSum`] with a tolerance of `1e-10` radians: far above the
    /// rounding error of an angle sum (about `1e-15`), so exact ties are left alone, and small
    /// enough that the cotan weights that remain negative are of the order of `1e-10` on
    /// triangles that are not needles.
    fn default() -> Self {
        DelaunayCriterion::AngleSum { tolerance: 1e-10 }
    }
}

/// How many flips [`IntrinsicTriangulation::flip_to_delaunay`] may make.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlipLimit {
    /// At most this many flips for every edge of the triangulation (the limit is the product).
    PerEdge(usize),
    /// At most this many flips.
    Flips(usize),
    /// No limit: the paper's algorithm. Termination is then guaranteed in exact arithmetic
    /// only.
    Unlimited,
}

impl Default for FlipLimit {
    /// `PerEdge(64)`. On the paper's data (Fig. 10) the number of flips rarely exceeds the
    /// number of edges.
    fn default() -> Self {
        FlipLimit::PerEdge(64)
    }
}

/// Options of [`IntrinsicTriangulation::flip_to_delaunay`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct DelaunayOptions {
    /// When an edge is flipped.
    pub criterion: DelaunayCriterion,
    /// The limit on the number of flips.
    pub max_flips: FlipLimit,
}

/// How Delaunay flipping ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DelaunayStatus {
    /// Every interior edge satisfies the criterion.
    Delaunay,
    /// The limit on the number of flips was reached; some edges do not satisfy the criterion.
    FlipLimit,
    /// The queue ran empty, but some edges that do not satisfy the criterion could not be
    /// flipped: in floating point arithmetic their two triangles are not a strictly convex
    /// quadrilateral, or the flip would create a triangle too thin to be described by its
    /// edge lengths ([`FlipRefusal`](crate::FlipRefusal)). This happens only with nearly
    /// degenerate triangles.
    Blocked,
}

/// What [`IntrinsicTriangulation::flip_to_delaunay`] did and what it left.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelaunayReport {
    /// How it ended.
    pub status: DelaunayStatus,
    /// Number of flips made.
    pub n_flips: usize,
    /// Number of edges taken from the queue and tested.
    pub n_tested: usize,
    /// Number of times an edge that failed the criterion was refused by the flip.
    pub n_refused: usize,
    /// Number of interior edges that do not satisfy the criterion at the end; zero exactly
    /// when the status is [`DelaunayStatus::Delaunay`].
    pub n_not_delaunay: usize,
    /// Number of interior edges at the end whose distance from a tie is within the tolerance
    /// of the criterion (`|α + β − π| ≤ tolerance`, or `|weight| ≤ tolerance`). The result is
    /// the Delaunay triangulation up to flips of these edges.
    pub n_ties: usize,
    /// The largest `α + β − π` over the interior edges at the end, in radians (`−∞` if there
    /// is no interior edge).
    pub max_angle_excess: f64,
    /// The smallest cotan weight `½ (cot α + cot β)` over the interior edges at the end (`+∞`
    /// if there is no interior edge).
    pub min_cotan_weight: f64,
    /// Number of constrained edges ([`Complex::is_constrained`]). They are left alone and are
    /// not counted in the numbers above, which are over the interior edges that are not
    /// constrained.
    pub n_constrained: usize,
}

/// A first-in first-out queue of edges that holds every edge at most once.
#[derive(Debug, Clone)]
pub(crate) struct EdgeQueue {
    queue: VecDeque<u32>,
    queued: Vec<bool>,
}

impl EdgeQueue {
    /// A queue of all edges, in ascending order.
    pub(crate) fn all(n_edges: usize) -> EdgeQueue {
        EdgeQueue {
            queue: (0..n_edges as u32).collect(),
            queued: vec![true; n_edges],
        }
    }

    /// An empty queue for a triangulation with `n_edges` edges.
    pub(crate) fn empty(n_edges: usize) -> EdgeQueue {
        EdgeQueue {
            queue: VecDeque::new(),
            queued: vec![false; n_edges],
        }
    }

    /// Puts an edge on the queue unless it is on it; the queue grows with the triangulation.
    pub(crate) fn push(&mut self, e: u32) {
        if e as usize >= self.queued.len() {
            self.queued.resize(e as usize + 1, false);
        }
        let slot = &mut self.queued[e as usize];
        if !*slot {
            *slot = true;
            self.queue.push_back(e);
        }
    }

    /// Follows a removal of a vertex: its removed edges leave the queue, the edges it
    /// renumbered keep their places under their new ids.
    pub(crate) fn apply_removal(&mut self, removal: &crate::remove::RemovalStep) {
        let n_edges = removal.counts()[1] as usize;
        let removed = removal.removed_edges();
        let mut queued = vec![false; n_edges];
        self.queue.retain_mut(|e| {
            if removed.contains(e) {
                return false;
            }
            if let Some(&(_, new)) = removal.edge_moves().iter().find(|m| m.0 == *e) {
                *e = new;
            }
            queued[*e as usize] = true;
            true
        });
        self.queued = queued;
    }

    fn pop(&mut self) -> Option<u32> {
        let e = self.queue.pop_front()?;
        self.queued[e as usize] = false;
        Some(e)
    }
}

/// Counters of one run of the queue.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct QueueRun {
    pub(crate) n_flips: usize,
    pub(crate) n_tested: usize,
    pub(crate) n_refused: usize,
    /// Flips that the shape limit of SPEC N5 refused and that were made under the second test
    /// given to [`run_flip_queue`].
    pub(crate) n_thin: usize,
    pub(crate) limit_reached: bool,
}

/// The five lengths `[ℓ_ij, ℓ_jk, ℓ_ki, ℓ_il, ℓ_lj]` of the two triangles `ijk`, `jil` on
/// interior edge `e = ij`, or `None` for a boundary edge. For an edge whose two sides are the
/// same triangle the triangle is listed twice, which is what the angle sums need.
#[inline]
pub(crate) fn diamond_lengths(complex: &Complex, lengths: &[f64], e: usize) -> Option<[f64; 5]> {
    let h0 = complex.edge_halfedge[e] as usize;
    let h1 = complex.twin[h0];
    if h1 == NONE {
        return None;
    }
    let h1 = h1 as usize;
    let len = |h: usize| lengths[complex.halfedge_edge[h] as usize];
    Some([
        lengths[e],
        len(next(h0)),
        len(prev(h0)),
        len(next(h1)),
        len(prev(h1)),
    ])
}

/// See [`IntrinsicTriangulation::edge_angle_excess`].
pub(crate) fn edge_angle_excess(complex: &Complex, lengths: &[f64], e: usize) -> Option<f64> {
    let [l_ij, l_jk, l_ki, l_il, l_lj] = diamond_lengths(complex, lengths, e)?;
    Some(
        delaunay_angle_excess(l_ij, l_jk, l_ki, l_il, l_lj)
            .expect("every face keeps the triangle inequalities"),
    )
}

/// See [`IntrinsicTriangulation::edge_cotan_weight`].
pub(crate) fn edge_cotan_weight(complex: &Complex, lengths: &[f64], e: usize) -> f64 {
    let h = complex.edge_halfedge[e] as usize;
    let mut sum = opposite_cotan(complex, lengths, h);
    if let Some(t) = complex.twin(h) {
        sum += opposite_cotan(complex, lengths, t);
    }
    0.5 * sum
}

/// Whether interior edge `e` fails the criterion. `false` for a boundary edge and for a
/// constrained edge.
#[inline]
pub(crate) fn fails_delaunay(
    complex: &Complex,
    lengths: &[f64],
    e: usize,
    criterion: &DelaunayCriterion,
) -> bool {
    if complex.constrained[e] {
        return false;
    }
    let Some([l_ij, l_jk, l_ki, l_il, l_lj]) = diamond_lengths(complex, lengths, e) else {
        return false;
    };
    let Some(cotans) = opposite_cotan_sum(l_ij, l_jk, l_ki, l_il, l_lj) else {
        return false;
    };
    match *criterion {
        DelaunayCriterion::CotanWeight { tolerance } => 0.5 * cotans < -tolerance,
        DelaunayCriterion::AngleSum { tolerance } => {
            // The cotangent sum has the sign of π − (α + β); only when it is negative are
            // the angles themselves needed.
            cotans < 0.0
                && delaunay_angle_excess(l_ij, l_jk, l_ki, l_il, l_lj)
                    .is_some_and(|excess| excess > tolerance)
        }
    }
}

/// Flips edges from the queue until it is empty or `budget` flips were made, with the flips
/// of the representation. The two faces of every flipped edge are appended to `touched`, if
/// given. A refused flip ([`IntrinsicError::FlipRefused`]) is counted and skipped; any other
/// error of the representation (e.g. an overflow of integer coordinates) ends the run and is
/// returned, so it is never mistaken for a refusal.
///
/// `thin`: a flip that the shape limit of SPEC N5 refuses ([`FlipRefusal::TooThin`]) is tried
/// again under this test (refinement, IC-SPEC I30); `None` leaves it refused.
///
/// [`FlipRefusal::TooThin`]: crate::FlipRefusal::TooThin
pub(crate) fn run_flip_queue<T: Triangulation>(
    tri: &mut T,
    queue: &mut EdgeQueue,
    criterion: &DelaunayCriterion,
    budget: usize,
    mut touched: Option<&mut Vec<u32>>,
    thin: Option<FlipShapeTest>,
) -> Result<QueueRun, IntrinsicError> {
    let mut run = QueueRun::default();
    while let Some(e) = queue.pop() {
        run.n_tested += 1;
        if !fails_delaunay(tri.complex(), tri.edge_lengths(), e as usize, criterion) {
            continue;
        }
        if run.n_flips >= budget {
            queue.push(e);
            run.limit_reached = true;
            break;
        }
        let mut flipped = tri.flip_edge(e as usize);
        if let (
            Some(shape),
            Err(IntrinsicError::FlipRefused {
                reason: FlipRefusal::TooThin { .. },
                ..
            }),
        ) = (thin, &flipped)
        {
            flipped = tri.flip_edge_with(
                e as usize,
                FlipRequest {
                    shape,
                    corner_angles: false,
                },
            );
            run.n_thin += usize::from(flipped.is_ok());
        }
        match flipped {
            Ok(_) => {}
            Err(IntrinsicError::FlipRefused { .. }) => {
                run.n_refused += 1;
                continue;
            }
            Err(other) => return Err(other),
        }
        run.n_flips += 1;
        let complex = tri.complex();
        let h = complex.edge_halfedge[e as usize] as usize;
        let t = complex.twin[h] as usize;
        for g in [next(h), prev(h), next(t), prev(t)] {
            queue.push(complex.halfedge_edge[g]);
        }
        if let Some(touched) = touched.as_deref_mut() {
            touched.extend([(h / 3) as u32, (t / 3) as u32]);
        }
    }
    Ok(run)
}

/// Delaunay flipping over any representation; see
/// [`IntrinsicTriangulation::flip_to_delaunay`].
pub(crate) fn flip_to_delaunay<T: Triangulation>(
    tri: &mut T,
    options: &DelaunayOptions,
) -> Result<DelaunayReport, IntrinsicError> {
    let tolerance = options.criterion.tolerance();
    if !(tolerance.is_finite() && tolerance >= 0.0) {
        return Err(IntrinsicError::InvalidOptions {
            reason: format!(
                "the tolerance of the Delaunay criterion must be finite and not negative, got \
                 {tolerance}"
            ),
        });
    }
    let n_edges = tri.complex().n_edges();
    let budget = match options.max_flips {
        FlipLimit::PerEdge(k) => k.saturating_mul(n_edges),
        FlipLimit::Flips(n) => n,
        FlipLimit::Unlimited => usize::MAX,
    };
    let mut queue = EdgeQueue::all(n_edges);
    let run = run_flip_queue(tri, &mut queue, &options.criterion, budget, None, None)?;

    let (complex, lengths) = (tri.complex(), tri.edge_lengths());
    let mut n_not_delaunay = 0;
    let mut n_ties = 0;
    let mut n_constrained = 0;
    let mut max_angle_excess = f64::NEG_INFINITY;
    let mut min_cotan_weight = f64::INFINITY;
    for e in 0..complex.n_edges() {
        if complex.constrained[e] {
            n_constrained += 1;
            continue;
        }
        let Some(excess) = edge_angle_excess(complex, lengths, e) else {
            continue;
        };
        let weight = edge_cotan_weight(complex, lengths, e);
        max_angle_excess = max_angle_excess.max(excess);
        min_cotan_weight = min_cotan_weight.min(weight);
        if fails_delaunay(complex, lengths, e, &options.criterion) {
            n_not_delaunay += 1;
        }
        let tie = match options.criterion {
            DelaunayCriterion::AngleSum { tolerance } => excess.abs() <= tolerance,
            DelaunayCriterion::CotanWeight { tolerance } => weight.abs() <= tolerance,
        };
        if tie {
            n_ties += 1;
        }
    }
    let status = if n_not_delaunay == 0 {
        DelaunayStatus::Delaunay
    } else if run.limit_reached {
        DelaunayStatus::FlipLimit
    } else {
        DelaunayStatus::Blocked
    };
    Ok(DelaunayReport {
        status,
        n_flips: run.n_flips,
        n_tested: run.n_tested,
        n_refused: run.n_refused,
        n_not_delaunay,
        n_ties,
        max_angle_excess,
        min_cotan_weight,
        n_constrained,
    })
}

impl IntrinsicTriangulation {
    /// The sum of the two angles opposite edge `e` minus `π`, in radians; `None` for a
    /// boundary edge. The edge is Delaunay exactly when this is not positive (paper Eq. 4).
    /// For an edge with an end point of degree one both angles belong to the same triangle.
    ///
    /// # Panics
    /// If `e` is out of range.
    pub fn edge_angle_excess(&self, e: usize) -> Option<f64> {
        edge_angle_excess(&self.complex, &self.edge_lengths, e)
    }

    /// The cotan weight of edge `e`: half the sum of the cotangents of the one or two angles
    /// opposite it. `−L_ij` of the [cotan Laplacian](Self::cotan_laplacian) is the sum of the
    /// weights of the edges between `i` and `j`.
    ///
    /// # Panics
    /// If `e` is out of range.
    pub fn edge_cotan_weight(&self, e: usize) -> f64 {
        edge_cotan_weight(&self.complex, &self.edge_lengths, e)
    }

    /// Whether interior edge `e` fails the criterion. `false` for a boundary edge and for a
    /// constrained edge.
    pub(crate) fn fails_delaunay(&self, e: usize, criterion: &DelaunayCriterion) -> bool {
        fails_delaunay(&self.complex, &self.edge_lengths, e, criterion)
    }

    /// Whether every interior edge that is not constrained satisfies the criterion.
    pub fn is_delaunay(&self, criterion: &DelaunayCriterion) -> bool {
        (0..self.n_edges()).all(|e| !self.fails_delaunay(e, criterion))
    }

    /// Flips edges until the triangulation is an intrinsic Delaunay triangulation (paper
    /// Sec. 4.1): all edges are put on a first-in first-out queue in ascending order; the
    /// first edge is taken off and, if it fails the criterion, flipped, and the four other
    /// edges of its two triangles are put back on the queue unless they are on it.
    ///
    /// Boundary edges and constrained edges ([`Complex::is_constrained`]) are never flipped
    /// and count as Delaunay. An edge that fails the criterion but cannot be flipped
    /// ([`FlipRefusal`](crate::FlipRefusal)) is left and counted; it is tested again when one
    /// of its neighbours is flipped.
    ///
    /// The same code runs for every representation ([`Triangulation::flip_to_delaunay`]).
    /// The report says how it ended; reaching the flip limit is not an error.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidOptions`] for a tolerance that is negative or not finite.
    pub fn flip_to_delaunay(
        &mut self,
        options: &DelaunayOptions,
    ) -> Result<DelaunayReport, IntrinsicError> {
        flip_to_delaunay(self, options)
    }
}
