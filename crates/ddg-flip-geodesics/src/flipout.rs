//! Local shortening, Alg. 1 `FlipOut(T, γ_abc)` (SPEC §4.1, paper Sec. 4.1 and App. A), and
//! the App. B move of a single-segment loop (SPEC §4.3).
//!
//! # The wedge and its outer arc
//!
//! The wedge on side `s` of a joint `γ_abc` at `b` is a fan of corners at `b`
//! ([`JointGeometry::corners_from_incoming`]), ordered from `ba` to `bc`: corner `i` is the
//! corner at `b` of the triangle `b n_i n_{i+1}`, with `n_0 = a`, `n_k = c`, `k = deg W`. A
//! corner is named by the halfedge `h_i` leaving `b` in its face, so the face's halfedge
//! opposite `b` is `next(h_i)`: the outer-arc edge `n_i n_{i+1}`. Going counter-clockwise
//! (the right wedge) `h_i` runs `b → n_i`; going clockwise (the left wedge) `h_i` runs
//! `b → n_{i+1}`. Hence, for `0 < j < k`,
//!
//! | | right wedge | left wedge |
//! |---|---|---|
//! | edge `b n_j` | `edge(h_j)` | `edge(h_{j−1})` |
//! | corner of face `j − 1` at `n_j` | `prev(h_{j−1})` | `next(h_{j−1})` |
//! | corner of face `j` at `n_j` | `next(h_j)` | `prev(h_j)` |
//! | outer-arc segment `n_i → n_{i+1}` | `next(h_i)` | twin of `next(h_i)` |
//!
//! and `β_j` is the sum of the two corners at `n_j` (SPEC id `fo-outer-arc-angle`, not
//! tagged until work package F5 adds it to Lean).
//!
//! Nothing here unfolds anything: the wedge is recomputed from the connectivity after every
//! flip (a flip renumbers the halfedges of its two faces, IT N9). That is what makes the
//! Δ-complex cases of App. A work without special code: a self-edge at `b` with both ends in
//! the wedge (twice incident) is the edge `b n_j` for two values of `j`, and flipping it
//! changes both; an outer-arc vertex of degree one has `face(h_{j−1}) = face(h_j)`.
//!
//! # The flip order
//!
//! As SPEC §4.1: scan `j = 1, 2, …`; flip `b n_j` when `β_j < π`; after a flip resume at
//! `j − 1` (only `β_{j−1}` can newly be `< π`), else go on with `j + 1`. Stop at `j = k`. A flip
//! the triangulation refuses (`TooThin`, or `NotConvex` by rounding although `β_j < π`, or a
//! topological refusal) is skipped and counted (F6); the scan goes on past it, as the
//! reference's `flipEdgeIfPossible` loop does (`gc: fg.cpp:785-819`). In exact arithmetic the
//! scan ends with every `β_j ≥ π` after at most `k₀²` flips (Thm A.7, A.8; SPEC ids
//! `fo-wedge-flip-valid`, `fo-flipout-terminates`); more than `k₀² + 2k₀ + 8` flips is reported
//! as [`FlipError::Internal`] rather than looping.

use std::f64::consts::PI;

use ddg_intrinsic::{Complex, Triangulation};

use crate::error::FlipError;
use crate::network::{Joint, Network, Replaced, Segment};
use crate::report::Shortening;
use crate::wedge::{JointGeometry, Side};

/// Counts of one Alg. 1 call.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct FlipOutStats {
    pub(crate) flips: usize,
    pub(crate) refused: usize,
    pub(crate) twice_incident: bool,
    pub(crate) degree_one: bool,
}

/// The result of one Alg. 1 call (or loop move): what happened, the counts, and the network
/// bookkeeping of the replacement (`None` when stalled).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FlipOutResult {
    pub(crate) shortening: Shortening,
    pub(crate) stats: FlipOutStats,
    pub(crate) replaced: Option<Replaced>,
}

/// The edge `b n_j` and `β_j` for `0 < j < k` (see the module documentation). SPEC id
/// `fo-outer-arc-angle` (not tagged until work package F5).
fn interior<T: Triangulation>(tri: &T, side: Side, corners: &[u32], j: usize) -> (usize, f64) {
    let c = tri.complex();
    let (hp, hj) = (corners[j - 1] as usize, corners[j] as usize);
    let (edge, at_prev, at_this) = match side {
        Side::Right => (c.edge(hj), c.prev(hp), c.next(hj)),
        Side::Left => (c.edge(hp), c.next(hp), c.prev(hj)),
    };
    (edge, tri.corner_angle(at_prev) + tri.corner_angle(at_this))
}

/// Whether the wedge has an edge incident on it twice (an interior edge `b n_j` that occurs
/// for two values of `j`), and whether an interior outer-arc vertex has degree one.
fn delta_cases(complex: &Complex, side: Side, corners: &[u32]) -> (bool, bool) {
    let k = corners.len();
    let mut edges: Vec<usize> = (1..k)
        .map(|j| {
            let h = match side {
                Side::Right => corners[j],
                Side::Left => corners[j - 1],
            };
            complex.edge(h as usize)
        })
        .collect();
    let degree_one =
        (1..k).any(|j| complex.face(corners[j - 1] as usize) == complex.face(corners[j] as usize));
    edges.sort_unstable();
    let twice = edges.windows(2).any(|w| w[0] == w[1]);
    (twice, degree_one)
}

/// The outer arc of the wedge as segments in path order, each with the halfedge of its edge
/// in the wedge's face (the side the segment is pushed on, SPEC §4.1 step 2).
fn outer_arc(complex: &Complex, side: Side, corners: &[u32]) -> Vec<(Segment, usize)> {
    corners
        .iter()
        .map(|&h| {
            let g = complex.next(h as usize);
            let along = complex.is_forward(g);
            let forward = match side {
                Side::Right => along,
                Side::Left => !along,
            };
            (Segment::new(complex.edge(g), forward), g)
        })
        .collect()
}

/// Alg. 1 on the wedge on `side` of `joint`, which the caller has checked: not at a pinned
/// vertex, flexible, not across the boundary, angle below `π`.
///
/// * `deg W = 0` (out and back along one edge, F26): the two segments are removed (F12).
/// * A single-segment loop (`ab = bc`): the App. B move ([`move_loop`]).
/// * Otherwise: flip until every `β_j ≥ π` (module documentation), then replace the joint by
///   the outer arc if it is strictly shorter than `ℓ_ab + ℓ_bc` (the length guard of
///   `gc: fg.cpp:836-838`, F6); else nothing is replaced and the flips are kept (stalled).
///
/// # Errors
/// [`FlipError::Internal`] when an edge to be flipped carries a segment (F7) or the flips do
/// not terminate within the bound of Thm A.8; errors of the triangulation's flips.
pub(crate) fn flip_out<T: Triangulation>(
    tri: &mut T,
    network: &mut Network,
    joint: Joint,
    side: Side,
) -> Result<FlipOutResult, FlipError> {
    if joint.incoming == joint.outgoing {
        return move_loop(tri, network, joint, side);
    }
    let mut stats = FlipOutStats::default();
    let length_of = |tri: &T, network: &Network, id: u32| {
        tri.edge_length(network.node(id).segment.edge as usize)
    };
    let old_length =
        length_of(tri, network, joint.incoming.0) + length_of(tri, network, joint.outgoing.0);
    let run = [joint.incoming.0, joint.outgoing.0];

    let mut geometry = JointGeometry::new(tri, network, joint)?;
    let k0 = geometry.wedge(side).degree();
    if k0 == 0 {
        let replaced = network.replace_run(tri.complex(), &run, &[])?;
        return Ok(FlipOutResult {
            shortening: Shortening::Contracted {
                old_length,
                path_contracted: replaced.contracted,
            },
            stats,
            replaced: Some(replaced),
        });
    }

    let max_flips = k0 * k0 + 2 * k0 + 8;
    let mut j = 1;
    let corners = loop {
        let corners = geometry.corners_from_incoming(side);
        let k = corners.len();
        let (twice, degree_one) = delta_cases(tri.complex(), side, &corners);
        stats.twice_incident |= twice;
        stats.degree_one |= degree_one;
        if j >= k {
            break corners;
        }
        let (edge, beta) = interior(tri, side, &corners, j);
        if beta < PI {
            if network.carries_segment(edge) {
                return Err(FlipError::Internal {
                    what: format!(
                        "Alg. 1 at vertex {}: edge {edge} inside a flexible wedge carries a \
                         segment (F7)",
                        geometry.vertex()
                    ),
                });
            }
            if tri.can_flip(edge)?.is_ok() {
                tri.flip_edge(edge)?;
                stats.flips += 1;
                if stats.flips > max_flips {
                    return Err(FlipError::Internal {
                        what: format!(
                            "Alg. 1 at vertex {} made more than {max_flips} flips on a wedge of \
                             degree {k0} (Thm A.8 bounds them by {})",
                            geometry.vertex(),
                            k0 * k0
                        ),
                    });
                }
                j = (j - 1).max(1);
                geometry = JointGeometry::new(tri, network, joint)?;
                continue;
            }
            stats.refused += 1;
        }
        j += 1;
    };

    let new = outer_arc(tri.complex(), side, &corners);
    let new_length: f64 = new
        .iter()
        .map(|(s, _)| tri.edge_length(s.edge as usize))
        .sum();
    if !(new_length < old_length) {
        return Ok(FlipOutResult {
            shortening: Shortening::Stalled {
                old_length,
                new_length,
            },
            stats,
            replaced: None,
        });
    }
    let replaced = network.replace_run(tri.complex(), &run, &new)?;
    Ok(FlipOutResult {
        shortening: Shortening::Replaced {
            new_segments: replaced
                .new_ids
                .iter()
                .map(|&id| crate::SegmentId(id))
                .collect(),
            old_length,
            new_length,
        },
        stats,
        replaced: Some(replaced),
    })
}

/// The App. B move of a single-segment loop `(aa)` (SPEC §4.3; `gc: fg.cpp:855-917`): the
/// loop is replaced by the two other edges of the triangle on the side of the wedge, an
/// isotopy across that face. On the left the face is that of the halfedge running like the
/// segment, `(a, a, x)` with `next = a → x` and `prev = x → a`, and the new loop is
/// `a → x` along `prev`'s edge, then `x → a` along `next`'s; on the right the face is that of
/// the opposite halfedge and the new loop runs along its `next`, then its `prev`. The new
/// segments are pushed on the side of that face. The length grows; App. B shows the next
/// step does not undo it (SPEC id `fo-loop-perturbation`, not tagged until F5).
fn move_loop<T: Triangulation>(
    tri: &mut T,
    network: &mut Network,
    joint: Joint,
    side: Side,
) -> Result<FlipOutResult, FlipError> {
    let complex = tri.complex();
    let segment = network.node(joint.outgoing.0).segment;
    let old_length = tri.edge_length(segment.edge as usize);
    let g = match side {
        Side::Left => segment.halfedge(complex),
        Side::Right => segment.reversed().halfedge(complex),
    };
    let Some(g) = g else {
        return Err(FlipError::Internal {
            what: format!(
                "single-segment loop on boundary edge {}: no face on its {side:?} side",
                segment.edge
            ),
        });
    };
    let (n, p) = (complex.next(g), complex.prev(g));
    let along = |h: usize, forward: bool| {
        (
            Segment::new(complex.edge(h), complex.is_forward(h) == forward),
            h,
        )
    };
    let new = match side {
        Side::Left => [along(p, false), along(n, false)],
        Side::Right => [along(n, true), along(p, true)],
    };
    let new_length =
        tri.edge_length(new[0].0.edge as usize) + tri.edge_length(new[1].0.edge as usize);
    let replaced = network.replace_run(complex, &[joint.outgoing.0], &new)?;
    Ok(FlipOutResult {
        shortening: Shortening::LoopMoved {
            new_segments: replaced
                .new_ids
                .iter()
                .map(|&id| crate::SegmentId(id))
                .collect(),
            old_length,
            new_length,
        },
        stats: FlipOutStats::default(),
        replaced: Some(replaced),
    })
}
