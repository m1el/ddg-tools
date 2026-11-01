//! Merging the two halves of a split edge back into one edge (SPEC §4.6, I27;
//! `ddg_intrinsic::merge`, IT-SPEC N53).
//!
//! An inserted vertex `m` with a straight star (two edges, the *halves*, on one geodesic
//! through it, as after an edge split) is removed by merging the halves: the decision and the
//! length `ℓ₁ + ℓ₂` come from [`ddg_intrinsic::plan_edge_merge`] (bitwise as the signposts),
//! the connectivity from [`Complex::merge_edge`](ddg_intrinsic::Complex::merge_edge) (the
//! inverse of `Complex::split_edge`, ids renumbered as a removal renumbers them).
//!
//! **The integers are exact.** The merged edge is the union of the two halves as curves, and
//! every other edge of the merged faces is an edge of the old ones, so the realization of `T⁰`
//! over the new `T¹` is the old one and only the merged edge has a new normal coordinate:
//! [`merge_edge_normal`]. Its crossings are those of the halves, plus one at `m` when `m` lies
//! on an input edge `E` and the halves do not run along it: `E` passes through `m`, where (by
//! SPEC §2.4 item 4) it ends two arcs of the curve of `E` at the vertex and is counted by
//! neither half; merged, it crosses the edge transversally there. If the halves run along `E`
//! (`n = −1`) the merged edge runs along `E` (`n = −1`). A merge never changes a roundabout:
//! the halfedge from the far end of each half towards `m` becomes the merged edge's halfedge
//! from that end and leaves in the same direction, and every other surviving halfedge is an
//! outer halfedge of the star. This inverts Eq. 6 (`ic-split-edge-normal`: halves `−1`, `−1`
//! along the split input edge) and `ic-split-crossed-edge` (halves `before`, `n − before`, the
//! split point inside an input face, never on a crossing): split then merge restores `n`
//! exactly, and the roundabouts, `q⁰` and the connectivity with their ids.
//!
//! **Combinatorially straight.** "Straight" is decided from the lengths, within
//! [`MERGE_STRAIGHT_TOLERANCE`](ddg_intrinsic::MERGE_STRAIGHT_TOLERANCE). A vertex inserted
//! within rounding of an input edge `E` (a `Face` point at a barycentric distance of `1e-10`
//! from `E`, say) can have a pair that is straight by its angles while the polyline of the two
//! halves passes on the other side of `E` than the geodesic between the far ends: then the
//! sum is wrong (measured: the merged edge was `E` itself with `n = 0` instead of `−1`, and the
//! certificate failed). So the merge also checks the integers (`check_straight_arcs`): on each
//! side, with faces `(m, X0, X1)` and `(m, X1, X2)`, no arc may run through the spoke `m X1`
//! next to `m` from the half `m X0` or the corner `X0` to the half `m X2` or the corner `X2`
//! (their number is `min(c_m + e_X0, c'_m + e'_X2)`); and a vertex on an input edge must have
//! one end of that input edge on each side. Otherwise the merge is refused with
//! [`RemovalRefusal::CurvesBend`] (the signposts, which cannot see it, merge there: such a
//! removal is `RemovalStuck` in this representation, as before I27 was resolved).

use ddg_intrinsic::complex::EdgeMerge;
use ddg_intrinsic::{IntrinsicError, RemovalRefusal, SurfacePoint, plan_edge_merge};

use crate::coords::IntegerTriangulation;
use crate::error::{IntegerError, Place, ValidityItem};

/// The normal coordinate of the edge made by merging two straight halves with normal
/// coordinates `n₁`, `n₂` at an inserted vertex `m` (see the [module documentation](self)):
///
/// ```text
/// n = −1               if n₁ = n₂ = −1 (both halves along the input edge through m)
/// n = n₁ + n₂ + 1      if n₁, n₂ ≥ 0 and m lies on an input edge
/// n = n₁ + n₂          if n₁, n₂ ≥ 0 and m lies inside an input face
/// ```
///
/// `None` where no straight edge has such halves (exactly one half along an input edge, or both
/// along one with `m` inside a face) and when the sum leaves the `i64` range.
///
/// [formula: ic-merge-edge-normal]
#[inline]
pub fn merge_edge_normal(n1: i64, n2: i64, m_on_input_edge: bool) -> Option<i64> {
    match (n1 == -1, n2 == -1) {
        (true, true) => m_on_input_edge.then_some(-1),
        (false, false) if n1 >= 0 && n2 >= 0 => {
            n1.checked_add(n2)?.checked_add(i64::from(m_on_input_edge))
        }
        _ => None,
    }
}

impl IntegerTriangulation {
    /// Merges the two halves at the inserted vertex `v` of which `e` is one back into one edge
    /// and removes `v` (SPEC §4.6, I27): [`plan_edge_merge`] (degree four, or three on the
    /// boundary; straight within
    /// [`MERGE_STRAIGHT_TOLERANCE`](ddg_intrinsic::MERGE_STRAIGHT_TOLERANCE); merged triangles
    /// valid; the length `ℓ₁ + ℓ₂`), then
    /// [`Complex::merge_edge`](ddg_intrinsic::Complex::merge_edge), whose [`EdgeMerge`] says how
    /// ids are renumbered (it moves `n`, `r`, the lengths and `q⁰`). The merged edge's normal
    /// coordinate is [`merge_edge_normal`]; nothing else changes. Merging right after
    /// [`split_edge`](Self::split_edge) gives back the state, bitwise except the merged length
    /// (`t ℓ + (1 − t) ℓ`, within an ulp of `ℓ`).
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::Intrinsic`] with
    /// [`IntrinsicError::RemovalRefused`] for an input vertex, the refusals of
    /// [`plan_edge_merge`], [`RemovalRefusal::HalvesDisagree`] when exactly one half runs
    /// along an input edge, and [`RemovalRefusal::CurvesBend`] when the input edges do not cross
    /// the halves as one straight edge (module documentation); [`IntegerError::InvalidCoordinates`]
    /// when both halves run along an input edge and `v` is not on one; [`IntegerError::Overflow`].
    /// The state is unchanged then.
    pub fn merge_edge(&mut self, v: usize, e: usize) -> Result<EdgeMerge, IntegerError> {
        let nv = self.complex.n_vertices();
        if v >= nv {
            return Err(IntegerError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: nv,
            });
        }
        if self.is_input_vertex(v) {
            return Err(IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                vertex: v as u32,
                reason: RemovalRefusal::InputVertex,
            }));
        }
        let plan = plan_edge_merge(&self.complex, &self.lengths, v, e)?;
        let [h0, h1] = plan.halves.map(|x| x as usize);
        let (n0, n1) = (self.normal[h0], self.normal[h1]);
        let on_edge = matches!(self.locations[v], SurfacePoint::Edge { .. });
        let Some(n) = merge_edge_normal(n0, n1, on_edge) else {
            return Err(if (n0 == -1) != (n1 == -1) {
                IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                    vertex: v as u32,
                    reason: RemovalRefusal::HalvesDisagree { edges: plan.halves },
                })
            } else if n0 == -1 {
                IntegerError::invalid(
                    ValidityItem::AlongInputEdge,
                    Place::Vertex(v as u32),
                    format!(
                        "the edges {h0} and {h1} at vertex {v} run along input edges (n = -1), \
                         but the vertex lies inside an input face ({})",
                        self.locations[v]
                    ),
                )
            } else {
                IntegerError::Overflow { edge: h0 as u32 }
            });
        };
        self.check_straight_arcs(v, plan.halves, n0 == -1)?;
        let merge = self.complex.merge_edge(v, e)?;
        let keep = merge.halves[0] as usize;
        self.lengths[keep] = plan.length;
        self.normal[keep] = n;
        merge.apply_halfedge_data(&mut self.roundabouts);
        merge.apply_edge_data(&mut self.lengths);
        merge.apply_edge_data(&mut self.normal);
        merge.apply_vertex_data(&mut self.locations);
        Ok(merge)
    }
}

impl IntegerTriangulation {
    /// Whether the input edges cross the star of `v` as they cross one straight edge through
    /// it (see the [module documentation](self), "Combinatorially straight"). On each side of
    /// the halves, with faces `(v, X0, X1)` and `(v, X1, X2)` and the spoke `v X1` between
    /// them, no arc may pass from the first face to the second through the spoke next to `v`
    /// while starting at the half `v X0` or at `X0` and ending at the half `v X2` or at `X2`:
    /// those are the arcs the polyline `X0 v X2` meets and a straight edge `X0 X2` would not
    /// (a bigon with the merged edge, an input edge from `X0` into the sliver, or the input
    /// edge `X0 X2` itself). Their number is `min(c_v + e_X0, c'_v + e'_X2)` (the arcs on the
    /// spoke counted from `v` in both faces, Eq. 1, 2). And when `v` lies on an input edge
    /// that the halves do not run along, each side holds exactly one of its two ends at `v`.
    fn check_straight_arcs(
        &self,
        v: usize,
        halves: [u32; 2],
        along: bool,
    ) -> Result<(), IntegerError> {
        let c = &self.complex;
        let mut out: Vec<usize> = c.outgoing_halfedges(v).collect();
        if !c.is_boundary_vertex(v)
            && let Some(i) = out.iter().position(|&h| c.edge(h) as u32 == halves[0])
        {
            out.rotate_left(i);
        }
        let on_edge = matches!(self.locations[v], SurfacePoint::Edge { .. });
        let bend = || {
            IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                vertex: v as u32,
                reason: RemovalRefusal::CurvesBend { edges: halves },
            })
        };
        for s in 0..out.len() / 2 {
            let (a, b) = (out[2 * s], out[2 * s + 1]);
            let (fa, fb) = (self.face_coordinates(a / 3)?, self.face_coordinates(b / 3)?);
            let (k, m) = (a % 3, b % 3);
            let near_first = fa.c[k] + fa.e[(k + 1) % 3];
            let near_second = fb.c[m] + fb.e[(m + 2) % 3];
            if near_first.min(near_second) > 0 {
                return Err(bend());
            }
            if on_edge && !along {
                let spoke_along = i64::from(self.normal[c.edge(b)] == -1);
                if fa.e[k] + fb.e[m] + spoke_along != 1 {
                    return Err(bend());
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::merge_edge_normal;

    #[test]
    fn merge_inverts_the_splits() {
        // Eq. 6 (n = −1): the halves keep −1; the merge gives −1 back.
        assert_eq!(merge_edge_normal(-1, -1, true), Some(-1));
        // A crossed edge split before `b` of its `n` crossings, inside an input face.
        for n in 0..20 {
            for before in 0..=n {
                assert_eq!(merge_edge_normal(before, n - before, false), Some(n));
            }
        }
        // Across an input edge at the vertex.
        assert_eq!(merge_edge_normal(0, 0, true), Some(1));
        assert_eq!(merge_edge_normal(3, 4, true), Some(8));
        // No straight edge.
        assert_eq!(merge_edge_normal(-1, 0, true), None);
        assert_eq!(merge_edge_normal(2, -1, false), None);
        assert_eq!(merge_edge_normal(-1, -1, false), None);
        assert_eq!(merge_edge_normal(i64::MAX, 1, false), None);
    }
}
