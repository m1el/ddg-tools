//! Local shortening, Alg. 1 `FlipOut(T, γ_abc)` (SPEC §4.1, paper Sec. 4.1 and App. A), and
//! the App. B move of a single-segment loop (SPEC §4.3).
//!
//! # The wedge and its outer arc
//!
//! The wedge on side `s` of a joint `γ_abc` at `b` is a fan of corners at `b`
//! ([`JointGeometry::corners_from_incoming`](crate::JointGeometry::corners_from_incoming)),
//! ordered from `ba` to `bc`: corner `i` is the corner at `b` of the triangle `b n_i n_{i+1}`,
//! with `n_0 = a`, `n_k = c`, `k = deg W`. A
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

use ddg_intrinsic::{Complex, FlipRequest, FlipShapeTest, IntrinsicError, Triangulation};

use crate::error::FlipError;
use crate::network::{Joint, Network, Segment};
use crate::report::Shortening;
use crate::scratch::{AngleCache, Located, Scratch};
use crate::wedge::Side;

/// Counts of one Alg. 1 call.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct FlipOutStats {
    pub(crate) flips: usize,
    pub(crate) refused: usize,
    pub(crate) twice_incident: bool,
    pub(crate) degree_one: bool,
}

/// What one Alg. 1 call (or loop move) did; the new segments of a replacement are
/// `scratch.replaced.new_ids`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Outcome {
    /// [`Shortening::Replaced`].
    Replaced { old_length: f64, new_length: f64 },
    /// [`Shortening::Contracted`].
    Contracted {
        old_length: f64,
        path_contracted: bool,
    },
    /// [`Shortening::LoopMoved`].
    LoopMoved { old_length: f64, new_length: f64 },
    /// [`Shortening::Stalled`].
    Stalled { old_length: f64, new_length: f64 },
}

impl Outcome {
    /// Whether the network changed (`scratch.replaced` is then this call's).
    pub(crate) fn replaced(self) -> bool {
        !matches!(self, Outcome::Stalled { .. })
    }

    /// [`Shortening::length_change`], bitwise.
    pub(crate) fn length_change(self) -> f64 {
        match self {
            Outcome::Replaced {
                old_length,
                new_length,
            }
            | Outcome::LoopMoved {
                old_length,
                new_length,
            } => new_length - old_length,
            Outcome::Contracted { old_length, .. } => -old_length,
            Outcome::Stalled { .. } => 0.0,
        }
    }

    /// The public [`Shortening`], with the new segments `new_ids`.
    pub(crate) fn shortening(self, new_ids: &[u32]) -> Shortening {
        let ids = || new_ids.iter().map(|&id| crate::SegmentId(id)).collect();
        match self {
            Outcome::Replaced {
                old_length,
                new_length,
            } => Shortening::Replaced {
                new_segments: ids(),
                old_length,
                new_length,
            },
            Outcome::Contracted {
                old_length,
                path_contracted,
            } => Shortening::Contracted {
                old_length,
                path_contracted,
            },
            Outcome::LoopMoved {
                old_length,
                new_length,
            } => Shortening::LoopMoved {
                new_segments: ids(),
                old_length,
                new_length,
            },
            Outcome::Stalled {
                old_length,
                new_length,
            } => Shortening::Stalled {
                old_length,
                new_length,
            },
        }
    }
}

/// The result of one Alg. 1 call (or loop move): what happened and the counts. The network
/// bookkeeping of a replacement is in `scratch.replaced`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FlipOutResult {
    pub(crate) outcome: Outcome,
    pub(crate) stats: FlipOutStats,
}

/// The edge `b n_j` and whether `β_j < π` for `0 < j < k` (see the module documentation):
/// bitwise the test `θ_prev + θ_this < π` on the two corner angles at `n_j`
/// ([`Triangulation::corner_angle`]), decided by [`AngleCache::sum_below_pi`]. SPEC id
/// `fo-outer-arc-angle` (not tagged until work package F5).
fn interior<T: Triangulation>(
    tri: &T,
    angles: &mut AngleCache,
    side: Side,
    corners: &[u32],
    j: usize,
) -> (usize, bool) {
    let c = tri.complex();
    let (hp, hj) = (corners[j - 1] as usize, corners[j] as usize);
    let (edge, at_prev, at_this) = match side {
        Side::Right => (c.edge(hj), c.prev(hp), c.next(hj)),
        Side::Left => (c.edge(hp), c.next(hp), c.prev(hj)),
    };
    (edge, angles.sum_below_pi(tri, at_prev, at_this))
}

/// Whether the wedge has an edge incident on it twice (an interior edge `b n_j` that occurs
/// for two values of `j`), and whether an interior outer-arc vertex has degree one. `edges`
/// is working memory.
fn delta_cases(
    complex: &Complex,
    side: Side,
    corners: &[u32],
    edges: &mut Vec<usize>,
) -> (bool, bool) {
    let k = corners.len();
    let interior = |j: usize| match side {
        Side::Right => corners[j] as usize,
        Side::Left => corners[j - 1] as usize,
    };
    let degree_one =
        (1..k).any(|j| complex.face(corners[j - 1] as usize) == complex.face(corners[j] as usize));
    // An edge met twice has both its halfedges leaving `b`: it is a loop at `b`. Only then
    // are the edges compared.
    let twice = (1..k).any(|j| {
        let h = interior(j);
        complex.from_vertex(h) == complex.to_vertex(h)
    }) && {
        edges.clear();
        edges.extend((1..k).map(|j| complex.edge(interior(j))));
        edges.sort_unstable();
        edges.windows(2).any(|w| w[0] == w[1])
    };
    (twice, degree_one)
}

/// The name of the corner that replaces corners `j − 1` and `j` of the wedge when `b n_j` is
/// flipped: the halfedge of `b n_j` leaving `b` (in its face it keeps the corner at `b`, now of
/// the triangle `b n_{j−1} n_{j+1}`). `None` when that edge is a loop or another corner of the
/// wedge lies in one of its two faces; the corners are then recomputed from the star.
fn merged_corner(complex: &Complex, side: Side, corners: &[u32], j: usize) -> Option<u32> {
    let h = match side {
        Side::Right => corners[j],
        Side::Left => corners[j - 1],
    };
    if complex.from_vertex(h as usize) == complex.to_vertex(h as usize) {
        return None;
    }
    // No other corner of the wedge lies in the two faces when `b` is a corner of each only
    // once.
    let b = complex.from_vertex(h as usize);
    let once = |c: u32| {
        let f = complex.face(c as usize);
        (0..3)
            .filter(|&i| complex.from_vertex(3 * f + i) == b)
            .count()
            == 1
    };
    (once(corners[j - 1]) && once(corners[j])).then_some(h)
}

/// What changed in the corners of the wedge since the Δ-complex cases were last evaluated.
#[derive(Clone, Copy)]
enum Changed {
    No,
    All,
    /// Two corners were merged into the one at this index.
    Merged(usize),
}

/// The outer arc of the wedge as segments in path order, each with the halfedge of its edge
/// in the wedge's face (the side the segment is pushed on, SPEC §4.1 step 2), into `out`.
fn outer_arc(complex: &Complex, side: Side, corners: &[u32], out: &mut Vec<(Segment, usize)>) {
    out.clear();
    out.extend(corners.iter().map(|&h| {
        let g = complex.next(h as usize);
        let along = complex.is_forward(g);
        let forward = match side {
            Side::Right => along,
            Side::Left => !along,
        };
        (Segment::new(complex.edge(g), forward), g)
    }));
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
/// `at` is the joint, located last in `scratch` (its star is current). After a flip only the
/// corners of the wedge are updated, not its angle: the flip of `b n_j` merges the corners
/// `j − 1` and `j` into the one named by the halfedge of `b n_j` that left `b` (the faces keep
/// their ids, `Complex::flip_edge`), the other corners keep their names; when the flipped edge
/// is a loop or another corner lies in one of its two faces (Δ-complex cases) the corners are
/// recomputed from the star of `b` instead, as [`JointGeometry`](crate::JointGeometry) does.
/// A flip is asked of [`Triangulation::flip_edge_with`] with the shape test `test`
/// ([`ShortenOptions::flip_shape`](crate::ShortenOptions::flip_shape)) and the request for its
/// two corner angles (IT-SPEC N52; the only caller that asks) directly, which refuses
/// exactly when [`Triangulation::can_flip_with`] would and then changes nothing.
///
/// # Errors
/// [`FlipError::Internal`] when an edge to be flipped carries a segment (F7) or the flips do
/// not terminate within the bound of Thm A.8; errors of the triangulation's flips.
pub(crate) fn flip_out<T: Triangulation>(
    tri: &mut T,
    network: &mut Network,
    at: Located,
    side: Side,
    scratch: &mut Scratch,
    test: FlipShapeTest,
) -> Result<FlipOutResult, FlipError> {
    let joint = at.joint;
    if joint.incoming == joint.outgoing {
        return move_loop(tri, network, joint, side, scratch);
    }
    let mut stats = FlipOutStats::default();
    let length_of = |tri: &T, network: &Network, id: u32| {
        tri.edge_length(network.node(id).segment.edge as usize)
    };
    let old_length =
        length_of(tri, network, joint.incoming.0) + length_of(tri, network, joint.outgoing.0);
    let run = [joint.incoming.0, joint.outgoing.0];

    let vertex = at.vertex;
    scratch.corners_from_incoming(&*tri, at, side);
    let k0 = scratch.corners.len();
    if k0 == 0 {
        network.replace_run(tri.complex(), &run, &[], &mut scratch.replaced)?;
        return Ok(FlipOutResult {
            outcome: Outcome::Contracted {
                old_length,
                path_contracted: scratch.replaced.contracted,
            },
            stats,
        });
    }

    let request = FlipRequest {
        shape: test,
        corner_angles: true,
    };
    let max_flips = k0 * k0 + 2 * k0 + 8;
    let mut j = 1;
    // The Δ-complex cases are evaluated whenever the corners change (they only accumulate).
    // After a merge at `m` the wedge's interior edges have lost one and gained none, and only
    // the pairs of consecutive corners at `m` are new: only those are tested.
    let mut changed = Changed::All;
    loop {
        let k = scratch.corners.len();
        match changed {
            Changed::All => {
                let (twice, degree_one) =
                    delta_cases(tri.complex(), side, &scratch.corners, &mut scratch.edges);
                stats.twice_incident |= twice;
                stats.degree_one |= degree_one;
            }
            Changed::Merged(m) => {
                let c = tri.complex();
                let face = |i: usize| c.face(scratch.corners[i] as usize);
                stats.degree_one |=
                    (m > 0 && face(m - 1) == face(m)) || (m + 1 < k && face(m) == face(m + 1));
                #[cfg(debug_assertions)]
                {
                    let (twice, degree_one) =
                        delta_cases(c, side, &scratch.corners, &mut scratch.edges);
                    assert!(
                        (!twice || stats.twice_incident) && (!degree_one || stats.degree_one),
                        "Alg. 1 at vertex {vertex}: Δ-complex cases missed after a merge"
                    );
                }
            }
            Changed::No => {}
        }
        changed = Changed::No;
        if j >= k {
            break;
        }
        let (edge, beta_below_pi) = interior(tri, &mut scratch.angles, side, &scratch.corners, j);
        if beta_below_pi {
            if network.carries_segment(edge) {
                return Err(FlipError::Internal {
                    what: format!(
                        "Alg. 1 at vertex {vertex}: edge {edge} inside a flexible wedge carries \
                         a segment (F7)"
                    ),
                });
            }
            let merge = merged_corner(tri.complex(), side, &scratch.corners, j);
            let faces = {
                let c = tri.complex();
                let h = c.edge_halfedge(edge);
                [c.face(h), c.twin(h).map_or(c.face(h), |t| c.face(t))]
            };
            match tri.flip_edge_with(edge, request) {
                Ok(flip) => {
                    // The two faces keep their ids and change their lengths; the flip computed
                    // a corner angle in each (IT-SPEC N52).
                    for f in faces {
                        scratch.angles.invalidate(f);
                    }
                    if let Some(angles) = flip.new_corner_angles {
                        for (h, angle) in flip.new_corners.into_iter().zip(angles) {
                            scratch.angles.insert(&*tri, h as usize, angle);
                        }
                    }
                    stats.flips += 1;
                    if stats.flips > max_flips {
                        return Err(FlipError::Internal {
                            what: format!(
                                "Alg. 1 at vertex {vertex} made more than {max_flips} flips on a \
                                 wedge of degree {k0} (Thm A.8 bounds them by {})",
                                k0 * k0
                            ),
                        });
                    }
                    changed = Changed::All;
                    match merge {
                        Some(h) => {
                            changed = Changed::Merged(j - 1);
                            scratch.corners[j - 1] = h;
                            scratch.corners.remove(j);
                            #[cfg(debug_assertions)]
                            {
                                let kept = scratch.corners.clone();
                                let at = scratch.locate(tri, network, joint)?;
                                scratch.corners_from_incoming(&*tri, at, side);
                                assert_eq!(
                                    kept, scratch.corners,
                                    "Alg. 1 at vertex {vertex}: corners merged after the flip of \
                                     edge {edge} differ from the star's"
                                );
                            }
                        }
                        None => {
                            let at = scratch.locate(tri, network, joint)?;
                            scratch.corners_from_incoming(&*tri, at, side);
                        }
                    }
                    j = (j - 1).max(1);
                    continue;
                }
                Err(IntrinsicError::FlipRefused { .. }) => stats.refused += 1,
                Err(e) => return Err(e.into()),
            }
        }
        j += 1;
    }

    outer_arc(tri.complex(), side, &scratch.corners, &mut scratch.new);
    let new_length: f64 = scratch
        .new
        .iter()
        .map(|(s, _)| tri.edge_length(s.edge as usize))
        .sum();
    if new_length.partial_cmp(&old_length) != Some(std::cmp::Ordering::Less) {
        return Ok(FlipOutResult {
            outcome: Outcome::Stalled {
                old_length,
                new_length,
            },
            stats,
        });
    }
    network.replace_run(tri.complex(), &run, &scratch.new, &mut scratch.replaced)?;
    Ok(FlipOutResult {
        outcome: Outcome::Replaced {
            old_length,
            new_length,
        },
        stats,
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
    scratch: &mut Scratch,
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
    network.replace_run(complex, &[joint.outgoing.0], &new, &mut scratch.replaced)?;
    Ok(FlipOutResult {
        outcome: Outcome::LoopMoved {
            old_length,
            new_length,
        },
        stats: FlipOutStats::default(),
    })
}
