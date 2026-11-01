//! Merging the two halves of a split edge back into one edge (SPEC N53; IC-SPEC §4.6, I27).
//!
//! A vertex made by [`split_edge`](crate::Triangulation::split_edge) has a *straight* star: two
//! of its edges (the halves) lie on one geodesic, so the corner angles between them sum to `π`
//! on each side. Removal by flips (Gillespie, Sharp, Crane 2021, Alg. 8) cannot always remove
//! such a vertex: the edges to the opposite corners have an angle of `π` at the vertex and are
//! refused as too thin, the halves flip only where their quadrilateral is convex at the far
//! end. Such a vertex is removed by undoing the split instead: [`Complex::merge_edge`] (the
//! topology) after [`plan_edge_merge`] (the decision and the length, from the connectivity and
//! the lengths alone, so every representation merges the same vertices with bitwise the same
//! length).
//!
//! **Straight.** The two halves are straight when the corner angles at the vertex between them
//! sum to `π` within [`MERGE_STRAIGHT_TOLERANCE`] on each side (on the boundary: the one side,
//! whose angle sum is that of the vertex). The merged edge then has the length `ℓ₁ + ℓ₂`; for
//! halves whose angle sums are `π ± δ` the true distance between the far ends is shorter by
//! about `ℓ₁ℓ₂ δ² / (2(ℓ₁ + ℓ₂))`, at most `ℓ δ²/8 ≈ 1.3e-19 ℓ` at the tolerance, below the
//! rounding of the sum. A merged triangle must satisfy the strict triangle inequality
//! ([`RemovalRefusal::Degenerate`]).

use std::f64::consts::PI;

use ddg_core::geometry::corner_angle;
use ddg_core::topology::NONE;

use crate::complex::{Complex, EdgeMerge, prev};
use crate::error::{IntrinsicError, RemovalRefusal};
use crate::flip::triangle_is_valid;
use crate::interface::face_lengths;
use crate::point::SurfacePoint;
use crate::surface::Surface;
use crate::triangulation::IntrinsicTriangulation;

/// How far from `π` (radians) the corner angles between two halves may sum on each side for
/// [`plan_edge_merge`] to call them straight. The corner angles of a split vertex come from the
/// lengths of [`plan_edge_split`](crate::plan_edge_split) and sum to `π` within rounding
/// (measured: within `5.0e-11` over 3 397 splits of the test shapes after random flips, needles
/// included, `tests/merge.rs`; SPEC N53); the edges to the opposite
/// corners of a star this close to straight are refused by every flip (their flip would make a
/// triangle of shape below `1e-9`, far under [`FLIP_SHAPE_LIMIT`](crate::formulas::FLIP_SHAPE_LIMIT)).
pub const MERGE_STRAIGHT_TOLERANCE: f64 = 1e-9;

/// What [`plan_edge_merge`] decided: the merge is possible, with these numbers.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeMergePlan {
    /// The two halves (edge ids): the edge given to [`plan_edge_merge`] first on an interior
    /// vertex; on the boundary the boundary edge leaving the vertex first.
    pub halves: [u32; 2],
    /// The length of the merged edge, `ℓ₁ + ℓ₂`.
    pub length: f64,
    /// The corner angles between the halves, summed on each side, minus `π` (the second `0`
    /// on the boundary).
    pub excess: [f64; 2],
    /// The star for [`Complex::merge_edge`].
    pub(crate) star: Vec<usize>,
}

/// The corner angle at the tail of halfedge `h` from the lengths, or the face as degenerate.
fn angle_at(complex: &Complex, lengths: &[f64], h: usize) -> Result<f64, IntrinsicError> {
    let l = face_lengths(complex, lengths, h / 3);
    let k = h % 3;
    corner_angle(l[(k + 1) % 3], l[k], l[(k + 2) % 3]).ok_or(IntrinsicError::DegenerateFace {
        surface: IntrinsicTriangulation::NAME,
        face: (h / 3) as u32,
        lengths: l,
    })
}

/// Checks that the two halves at vertex `v` of which `e` is one can be merged
/// ([`Complex::merge_edge`]) in a triangulation with this connectivity and these lengths, and
/// returns the merged length: the topological pattern of [`Complex::merge_edge`], the halves
/// straight within [`MERGE_STRAIGHT_TOLERANCE`], each merged triangle valid. Whether `v` may be
/// removed at all (an input vertex may not) is the representation's business.
///
/// # Errors
/// [`IntrinsicError::IndexOutOfRange`], [`IntrinsicError::DegenerateFace`], and
/// [`IntrinsicError::RemovalRefused`] with the refusals of [`Complex::merge_edge`],
/// [`RemovalRefusal::NotStraight`] and [`RemovalRefusal::Degenerate`] (the lengths of the first
/// merged triangle that is not valid).
pub fn plan_edge_merge(
    complex: &Complex,
    lengths: &[f64],
    v: usize,
    e: usize,
) -> Result<EdgeMergePlan, IntrinsicError> {
    for (what, index, bound) in [
        ("vertex", v, complex.n_vertices()),
        ("edge", e, complex.n_edges()),
    ] {
        if index >= bound {
            return Err(IntrinsicError::IndexOutOfRange { what, index, bound });
        }
    }
    let refused = |reason| IntrinsicError::RemovalRefused {
        vertex: v as u32,
        reason,
    };
    let star = complex.merge_star(v, e).map_err(refused)?;
    let halves = complex.merge_halves(&star);
    let n_sides = star.len() / 2;
    let mut excess = [0.0; 2];
    for (s, x) in excess.iter_mut().enumerate().take(n_sides) {
        let (a, b) = (star[2 * s], star[2 * s + 1]);
        *x = angle_at(complex, lengths, a)? + angle_at(complex, lengths, b)? - PI;
    }
    if !excess.iter().all(|x| x.abs() <= MERGE_STRAIGHT_TOLERANCE) {
        return Err(refused(RemovalRefusal::NotStraight {
            edges: halves.map(|x| x as u32),
            excess,
        }));
    }
    let length = lengths[halves[0]] + lengths[halves[1]];
    for s in 0..n_sides {
        let (a, b) = (star[2 * s], star[2 * s + 1]);
        let tri = [
            length,
            lengths[complex.edge(complex.next(a))],
            lengths[complex.edge(complex.next(b))],
        ];
        if !triangle_is_valid(tri[0], tri[1], tri[2]) {
            return Err(refused(RemovalRefusal::Degenerate { lengths: tri }));
        }
    }
    Ok(EdgeMergePlan {
        halves: halves.map(|x| x as u32),
        length,
        excess,
        star,
    })
}

/// Two edges at vertex `v` that lie on one straight line through it, if there are such: the
/// pair whose corner angles between them sum to `π` on both sides (interior) or the two
/// boundary edges of a boundary vertex whose angle sum is `π`, within
/// [`MERGE_STRAIGHT_TOLERANCE`], at any degree. Of several pairs the one with the smallest
/// largest excess wins, ties by the order of the outgoing halfedges from
/// [`Complex::outgoing_halfedge`]. Self-edges at `v` and degenerate faces give `None`. Returns
/// the two edges, the first one the half of the earlier outgoing halfedge.
pub fn straight_pair(complex: &Complex, lengths: &[f64], v: usize) -> Option<[u32; 2]> {
    let out: Vec<usize> = complex.outgoing_halfedges(v).collect();
    if out.iter().any(|&h| complex.to_vertex(h) == v) {
        return None;
    }
    let angles: Option<Vec<f64>> = out
        .iter()
        .map(|&h| angle_at(complex, lengths, h).ok())
        .collect();
    let angles = angles?;
    if complex.is_boundary_vertex(v) {
        let sum: f64 = angles.iter().sum();
        if (sum - PI).abs() > MERGE_STRAIGHT_TOLERANCE || out.len() < 2 {
            return None;
        }
        let last = *out.last()?;
        return Some([complex.edge(out[0]) as u32, complex.edge(prev(last)) as u32]);
    }
    let d = out.len();
    let mut best: Option<(f64, usize, usize)> = None;
    for i in 0..d {
        for j in i + 2..d {
            if i + d - j < 2 {
                continue;
            }
            let first: f64 = angles[i..j].iter().sum();
            let second: f64 = angles[j..].iter().chain(&angles[..i]).sum();
            let worst = (first - PI).abs().max((second - PI).abs());
            if worst <= MERGE_STRAIGHT_TOLERANCE && best.is_none_or(|b| worst < b.0) {
                best = Some((worst, i, j));
            }
        }
    }
    let (_, i, j) = best?;
    let pair = [complex.edge(out[i]) as u32, complex.edge(out[j]) as u32];
    (pair[0] != pair[1] && pair.iter().all(|&x| x != NONE)).then_some(pair)
}

impl IntrinsicTriangulation {
    /// Merges the two halves at the inserted vertex `v` of which `e` is one back into one edge
    /// and removes `v` (SPEC N53): [`plan_edge_merge`], then [`Complex::merge_edge`], whose
    /// [`EdgeMerge`] says how ids are renumbered. The merged edge gets the length `ℓ₁ + ℓ₂`
    /// and lies along the input edge both halves lie along (if they do; or one half does, the
    /// other is not known to lie along any, and `v` lies on that input edge); the surviving
    /// halfedges keep their signposts (they leave their tails in the same directions);
    /// nothing else changes. Merging right after [`split_edge`](Self::split_edge) gives back
    /// the triangulation, bitwise except the merged length (`t ℓ + (1 − t) ℓ`, within an ulp
    /// or two of `ℓ`).
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`]; [`IntrinsicError::RemovalRefused`] with
    /// [`RemovalRefusal::InputVertex`], the refusals of [`plan_edge_merge`], and
    /// [`RemovalRefusal::HalvesDisagree`] when the halves lie along different input edges, or
    /// one does and `v` is not on that input edge. The triangulation is unchanged then.
    pub fn merge_edge(&mut self, v: usize, e: usize) -> Result<EdgeMerge, IntrinsicError> {
        if v >= self.n_vertices() {
            return Err(IntrinsicError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: self.n_vertices(),
            });
        }
        let refused = |reason| IntrinsicError::RemovalRefused {
            vertex: v as u32,
            reason,
        };
        if v < self.n_input_vertices() {
            return Err(refused(RemovalRefusal::InputVertex));
        }
        let plan = plan_edge_merge(&self.complex, &self.edge_lengths, v, e)?;
        let [h0, h1] = plan.halves.map(|x| x as usize);
        // The input edge the merged edge lies along: that of both halves; or that of one half
        // when the other is not known to lie along an input edge and `v` lies on that input
        // edge (then the straight continuation through `v` lies along it too).
        let (x0, x1) = (self.edge_on_input[h0], self.edge_on_input[h1]);
        let on_v = match self.locations[v] {
            SurfacePoint::Edge { edge, .. } => edge,
            _ => NONE,
        };
        let on_input = match (x0, x1) {
            (x, y) if x == y => x,
            (x, NONE) | (NONE, x) if x == on_v => x,
            _ => {
                return Err(refused(RemovalRefusal::HalvesDisagree {
                    edges: plan.halves,
                }));
            }
        };
        let merge = self.complex.merge_edge_unchecked(v, &plan.star);
        let keep = merge.halves[0] as usize;
        self.edge_lengths[keep] = plan.length;
        self.edge_on_input[keep] = on_input;
        merge.apply_halfedge_data(&mut self.signposts);
        merge.apply_edge_data(&mut self.edge_lengths);
        merge.apply_edge_data(&mut self.edge_on_input);
        merge.apply_vertex_data(&mut self.angle_sums);
        merge.apply_vertex_data(&mut self.locations);
        Ok(merge)
    }
}
