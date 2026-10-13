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
    /// The two triangles on edge `e`, or why there are no two.
    fn quad(&self, e: usize) -> Result<Quad, FlipRefusal> {
        let h0 = self.edge_halfedge[e] as usize;
        let h1 = self.twin[h0];
        if h1 == NONE {
            return Err(FlipRefusal::BoundaryEdge);
        }
        let h1 = h1 as usize;
        if h0 / 3 == h1 / 3 {
            // The face is glued to itself along the edge. Its third halfedge is a loop at one
            // end point; the other end point has this one corner only.
            let loop_halfedge = 3 * (h0 / 3) + (3 - h0 % 3 - h1 % 3);
            let at_loop = self.from_vertex(loop_halfedge);
            let [a, b] = [self.from_vertex(h0), self.to_vertex(h0)];
            return Err(FlipRefusal::DegreeOneEndpoint {
                vertex: (if a == at_loop { b } else { a }) as u32,
            });
        }
        let len = |h: usize| self.edge_lengths[self.halfedge_edge[h] as usize];
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

    /// The five lengths `[ℓ_ij, ℓ_jk, ℓ_ki, ℓ_il, ℓ_lj]` of the two triangles `ijk`, `jil` on
    /// interior edge `e = ij`, or `None` for a boundary edge. For an edge whose two sides are
    /// the same triangle the triangle is listed twice, which is what the angle sums need.
    pub(crate) fn diamond_lengths(&self, e: usize) -> Option<[f64; 5]> {
        let h0 = self.edge_halfedge[e] as usize;
        let h1 = self.twin[h0];
        if h1 == NONE {
            return None;
        }
        let h1 = h1 as usize;
        let len = |h: usize| self.edge_lengths[self.halfedge_edge[h] as usize];
        Some([
            self.edge_lengths[e],
            len(self.next(h0)),
            len(self.prev(h0)),
            len(self.next(h1)),
            len(self.prev(h1)),
        ])
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
        let (f0, a) = (h0 / 3, h0 % 3);
        let (f1, b) = (h1 / 3, h1 % 3);
        // Before:  f0 = (i, j, k) with halfedges h0 = i→j, h0n = j→k, h0p = k→i,
        //          f1 = (j, i, l) with halfedges h1 = j→i, h1n = i→l, h1p = l→j.
        // After:   f0 = (i, l, k) with halfedges h0 = i→l, h0n = l→k, h0p = k→i,
        //          f1 = (j, k, l) with halfedges h1 = j→k, h1n = k→l, h1p = l→j.
        // The halfedges k→i and l→j stay where they are; i→l and j→k change faces.
        let (h0n, h0p) = (3 * f0 + (a + 1) % 3, 3 * f0 + (a + 2) % 3);
        let (h1n, h1p) = (3 * f1 + (b + 1) % 3, 3 * f1 + (b + 2) % 3);
        let (i, j) = (self.faces[f0][a], self.faces[f1][b]);
        let (k, l) = (self.faces[f0][(a + 2) % 3], self.faces[f1][(b + 2) % 3]);

        // Twins, edges and signposts of the two halfedges that move.
        let moved = |h: u32| {
            if h == h0n as u32 {
                h1 as u32
            } else if h == h1n as u32 {
                h0 as u32
            } else {
                h
            }
        };
        let twin_il = moved(self.twin[h1n]);
        let twin_jk = moved(self.twin[h0n]);
        let (edge_il, edge_jk) = (self.halfedge_edge[h1n], self.halfedge_edge[h0n]);
        let (phi_il, phi_jk) = (self.signposts[h1n], self.signposts[h0n]);

        self.faces[f0][(a + 1) % 3] = l;
        self.faces[f1][(b + 1) % 3] = k;

        self.twin[h0] = twin_il;
        if twin_il != NONE {
            self.twin[twin_il as usize] = h0 as u32;
        }
        self.twin[h1] = twin_jk;
        if twin_jk != NONE {
            self.twin[twin_jk as usize] = h1 as u32;
        }
        self.twin[h0n] = h1n as u32;
        self.twin[h1n] = h0n as u32;

        self.halfedge_edge[h0] = edge_il;
        self.halfedge_edge[h1] = edge_jk;
        self.halfedge_edge[h0n] = e as u32;
        self.halfedge_edge[h1n] = e as u32;
        if self.edge_halfedge[edge_il as usize] == h1n as u32 {
            self.edge_halfedge[edge_il as usize] = h0 as u32;
        }
        if self.edge_halfedge[edge_jk as usize] == h0n as u32 {
            self.edge_halfedge[edge_jk as usize] = h1 as u32;
        }
        self.edge_halfedge[e] = h0n as u32;

        if self.vertex_halfedge[i as usize] == h1n as u32 {
            self.vertex_halfedge[i as usize] = h0 as u32;
        }
        if self.vertex_halfedge[j as usize] == h0n as u32 {
            self.vertex_halfedge[j as usize] = h1 as u32;
        }

        self.edge_lengths[e] = new_length;
        self.edge_on_input[e] = NONE;

        // Signposts: i→l and j→k keep theirs; l→k follows l→j, k→l follows k→i.
        self.signposts[h0] = phi_il;
        self.signposts[h1] = phi_jk;
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
