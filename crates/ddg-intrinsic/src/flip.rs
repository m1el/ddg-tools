//! The edge flip (paper Sec. 3.3.1, Alg. 5).

use ddg_core::geometry::{corner_angle, corner_cotan, triangle_area};
use ddg_core::topology::NONE;

use crate::error::{FlipRefusal, IntrinsicError};
use crate::formulas::{FLIP_SHAPE_LIMIT, signpost_update, unfold_diamond};
use crate::triangulation::IntrinsicTriangulation;

/// What a flip did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Flip {
    /// The edge; it keeps its id.
    pub edge: u32,
    /// Its end points before the flip (tail and head of its halfedge).
    pub old_vertices: [u32; 2],
    /// Its end points after the flip: the vertex that was opposite on the right of the old
    /// halfedge, then the one that was opposite on its left.
    pub new_vertices: [u32; 2],
    /// Its length before the flip.
    pub old_length: f64,
    /// Its length after the flip.
    pub new_length: f64,
}

/// Whether three lengths form a triangle that every formula of this crate can be evaluated on:
/// positive area, and positive finite angles and finite cotangents at all corners. This is the
/// standard of `ddg_core::Geometry`. For lengths between `1e-100` and `1e100` a positive area
/// implies the rest (no product of two factors can under- or overflow), so only the area is
/// computed.
pub(crate) fn triangle_is_valid(a: f64, b: f64, c: f64) -> bool {
    const LO: f64 = 1e-100;
    const HI: f64 = 1e100;
    if triangle_area(a, b, c).is_none() {
        return false;
    }
    let sane = |x: f64| (LO..=HI).contains(&x);
    if sane(a) && sane(b) && sane(c) {
        return true;
    }
    [(a, b, c), (b, c, a), (c, a, b)]
        .iter()
        .all(|&(o, s, t)| corner_angle(o, s, t).is_some() && corner_cotan(o, s, t).is_some())
}

/// The six halfedges of the two triangles on an interior edge and its five lengths.
struct Quad {
    /// `i → j`, the edge's halfedge, in face `f0`.
    h0: usize,
    /// `j → i`, in face `f1`.
    h1: usize,
    /// `[ℓ_ij, ℓ_jk, ℓ_ki, ℓ_il, ℓ_lj]`.
    lengths: [f64; 5],
}

impl IntrinsicTriangulation {
    /// The two triangles on edge `e`, or why there are no two (or why the edge is left
    /// alone).
    fn quad(&self, e: usize) -> Result<Quad, FlipRefusal> {
        if let Some(reason) = self.complex.flip_refusal(e) {
            return Err(reason);
        }
        let h0 = self.complex.edge_halfedge[e] as usize;
        let h1 = self.complex.twin[h0] as usize;
        let len = |h: usize| self.edge_lengths[self.complex.halfedge_edge[h] as usize];
        Ok(Quad {
            h0,
            h1,
            lengths: [
                self.edge_lengths[e],
                len(self.next(h0)),
                len(self.prev(h0)),
                len(self.next(h1)),
                len(self.prev(h1)),
            ],
        })
    }

    /// Flips edge `e`: its two triangles `ijk`, `jil` are replaced by `ilk`, `jkl` (paper
    /// Alg. 5). The new length is the other diagonal of the two triangles unfolded into the
    /// plane ([`Diamond::flipped_length`](crate::formulas::Diamond::flipped_length)), and the
    /// signposts of the two new halfedges follow from those of their clockwise neighbours
    /// ([`signpost_update`]). Nothing else changes: the edge keeps its id, the two faces keep
    /// theirs, no vertex moves.
    ///
    /// Flipping the same edge again restores its end points, its length and its signposts up
    /// to rounding (with its direction reversed; four flips restore that too).
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], and [`IntrinsicError::FlipRefused`] with the
    /// reason ([`FlipRefusal`]): a boundary edge, an end point of degree one, two triangles
    /// that are not a strictly convex quadrilateral, or a quadrilateral so barely convex that
    /// the flip would create a triangle too thin to be described by its edge lengths
    /// ([`FLIP_SHAPE_LIMIT`]). The triangulation is unchanged when the flip is refused.
    pub fn flip_edge(&mut self, e: usize) -> Result<Flip, IntrinsicError> {
        self.check_edge(e)?;
        self.try_flip(e).map_err(|reason| {
            let [a, b] = self.edge_vertices(e);
            IntrinsicError::FlipRefused {
                edge: e as u32,
                vertices: [a as u32, b as u32],
                reason,
            }
        })
    }

    /// Whether [`flip_edge`](Self::flip_edge) would flip edge `e`, without changing anything.
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`]; a refusal is the `Err` of the inner result.
    pub fn can_flip(&self, e: usize) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        self.check_edge(e)?;
        Ok(self.plan_flip(e).map(|_| ()))
    }

    /// Checks that edge `e` can be flipped and returns its halfedges and its new length.
    fn plan_flip(&self, e: usize) -> Result<(Quad, f64), FlipRefusal> {
        let quad = self.quad(e)?;
        let [l_ij, l_jk, l_ki, l_il, l_lj] = quad.lengths;
        // The lengths of valid faces always unfold; `None` cannot happen.
        let diamond = unfold_diamond(l_ij, l_jk, l_ki, l_il, l_lj).ok_or(FlipRefusal::TooThin {
            new_length: 0.0,
            shape: 0.0,
        })?;
        let [area_at_tail, area_at_head] = diamond.flipped_areas();
        if !(area_at_tail > 0.0 && area_at_head > 0.0) {
            return Err(FlipRefusal::NotConvex {
                area_at_tail,
                area_at_head,
            });
        }
        let new_length = diamond.flipped_length();
        // The triangles i l k and j k l that the flip creates (SPEC N5).
        let shape = (area_at_tail / (new_length * new_length.max(l_il).max(l_ki)))
            .min(area_at_head / (new_length * new_length.max(l_jk).max(l_lj)));
        if !(shape >= FLIP_SHAPE_LIMIT
            && triangle_is_valid(l_il, new_length, l_ki)
            && triangle_is_valid(l_jk, new_length, l_lj))
        {
            return Err(FlipRefusal::TooThin { new_length, shape });
        }
        Ok((quad, new_length))
    }

    /// [`flip_edge`](Self::flip_edge) for an edge id that is known to be in range.
    pub(crate) fn try_flip(&mut self, e: usize) -> Result<Flip, FlipRefusal> {
        let (quad, new_length) = self.plan_flip(e)?;
        let [l_ij, l_jk, l_ki, l_il, l_lj] = quad.lengths;
        let (h0, h1) = (quad.h0, quad.h1);
        // The topology (see `ComplexFlip`): f0 = (i, j, k), f1 = (j, i, l) become
        // f0 = (i, l, k), f1 = (j, k, l); i→l and j→k change faces, keeping their signposts.
        let flip = self.complex.flip_unchecked(e);
        flip.move_halfedge_data(&mut self.signposts);
        let [i, j] = flip.old_vertices;
        let [l, k] = flip.new_vertices;
        let [h0n, h1n] = flip.new_halfedges.map(|h| h as usize);
        let (h0p, h1p) = (self.prev(h0), self.prev(h1));

        self.edge_lengths[e] = new_length;
        self.edge_on_input[e] = NONE;

        // Signposts: l→k follows l→j, k→l follows k→i.
        let (k, l) = (k as usize, l as usize);
        // Both angles exist: the two new triangles were checked by `plan_flip`.
        let at_l = corner_angle(l_jk, l_lj, new_length).unwrap_or(0.0);
        let at_k = corner_angle(l_il, l_ki, new_length).unwrap_or(0.0);
        self.signposts[h0n] = signpost_update(
            self.signposts[h1p],
            at_l,
            self.angle_sums[l],
            self.is_boundary_vertex(l),
        );
        self.signposts[h1n] = signpost_update(
            self.signposts[h0p],
            at_k,
            self.angle_sums[k],
            self.is_boundary_vertex(k),
        );

        Ok(Flip {
            edge: e as u32,
            old_vertices: [i, j],
            new_vertices: [l as u32, k as u32],
            old_length: l_ij,
            new_length,
        })
    }
}
