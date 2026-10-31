//! The edge flip (paper Sec. 3.3.1, Alg. 5).

use ddg_core::geometry::{corner_angle, corner_cotan, triangle_area};
use ddg_core::topology::NONE;

use crate::complex::{Complex, next, prev};
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
    /// The two corner angles the flip computes at the new edge's end points (SPEC N52), in the
    /// order of [`new_vertices`](Self::new_vertices): the angle at `l` in the new face
    /// `(j, k, l)` (the corner of the halfedge `l → j`, [`new_corners`](Self::new_corners)`[0]`)
    /// and the angle at `k` in the new face `(i, l, k)` (the corner of the halfedge `k → i`,
    /// `new_corners[1]`). Each is bitwise what
    /// [`Triangulation::corner_angle`](crate::Triangulation::corner_angle) gives for that
    /// corner after the flip. The other two corners at `l` and `k` (of the halfedges `l → k`
    /// and `k → l`) are not computed.
    pub new_corner_angles: [f64; 2],
    /// The halfedges whose corners [`new_corner_angles`](Self::new_corner_angles) are: `l → j`
    /// in face `(j, k, l)` and `k → i` in face `(i, l, k)`, ids after the flip (they are the
    /// halfedges that stayed in place, [`ComplexFlip`](crate::complex::ComplexFlip)).
    pub new_corners: [u32; 2],
}

/// The test a flip's two new triangles must pass besides strict convexity and the triangle
/// inequalities (SPEC N5, N52). The caller chooses it per flip
/// ([`Triangulation::flip_edge_with`](crate::Triangulation::flip_edge_with),
/// [`plan_flip_with`]); every operation of this crate, and
/// [`flip_edge`](crate::Triangulation::flip_edge), use the default, [`FlipShapeTest::N5`].
///
/// Both tests are invariant under scaling. With `A_tail`, `A_head` twice the areas of the new
/// triangles `ilk` (at the old tail `i`) and `jkl` (at the old head `j`), `ℓ_new` the new length
/// and `ℓ_max(T)` the longest side of triangle `T`:
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FlipShapeTest {
    /// SPEC N5: refused unless each new triangle `T` has `A_T ≥ limit · ℓ_new · ℓ_max(T)`
    /// (its shape number against a fixed limit; the default is
    /// [`FLIP_SHAPE_LIMIT`](crate::formulas::FLIP_SHAPE_LIMIT)).
    Absolute {
        /// The smallest shape number accepted.
        limit: f64,
    },
    /// The reference's test (geometry-central `flipEdgeIfPossible`,
    /// `signpost_intrinsic_triangulation.cpp:333-338`): refused unless
    /// `min(A_tail, A_head) ≥ ratio · (A_tail + A_head)`, each new triangle's area against the
    /// two together. It accepts a needle whose short side is not the new edge as long as the
    /// two areas are balanced, which the absolute test refuses.
    Relative {
        /// The smallest share of the two areas together that each must have.
        ratio: f64,
    },
}

impl FlipShapeTest {
    /// The default: SPEC N5, `Absolute { limit: FLIP_SHAPE_LIMIT }`.
    pub const N5: FlipShapeTest = FlipShapeTest::Absolute {
        limit: FLIP_SHAPE_LIMIT,
    };
    /// The reference's relative test with its constant, `Relative { ratio: 1e-6 }`
    /// (`triangleTestEPS`).
    pub const REFERENCE: FlipShapeTest = FlipShapeTest::Relative { ratio: 1e-6 };

    /// Whether the two new triangles pass, from twice their areas (`A_tail`, `A_head`) and
    /// their shape numbers' smaller value `shape` (SPEC N5). `false` for NaN.
    #[inline]
    fn accepts(self, area_at_tail: f64, area_at_head: f64, shape: f64) -> bool {
        match self {
            FlipShapeTest::Absolute { limit } => shape >= limit,
            FlipShapeTest::Relative { ratio } => {
                area_at_tail.min(area_at_head) >= ratio * (area_at_tail + area_at_head)
            }
        }
    }
}

impl Default for FlipShapeTest {
    fn default() -> Self {
        FlipShapeTest::N5
    }
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

/// What a flip of an edge would do to the lengths, decided before anything changes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlipPlan {
    /// `[ℓ_ij, ℓ_jk, ℓ_ki, ℓ_il, ℓ_lj]`: the edge `ij` and the sides of its two triangles
    /// `ijk` (on the side of its halfedge) and `jil`.
    pub lengths: [f64; 5],
    /// The length of the flipped edge `lk` (`it-flip-length`).
    pub new_length: f64,
}

impl FlipPlan {
    /// The corner angles at the new edge's end points that a flip by this plan creates, as
    /// [`Flip::new_corner_angles`]: at `l` in `(j, k, l)`, at `k` in `(i, l, k)`.
    #[inline]
    pub fn new_corner_angles(&self) -> [f64; 2] {
        let [_, l_jk, l_ki, l_il, l_lj] = self.lengths;
        // Both angles exist: the two new triangles were checked by `plan_flip_with`.
        [
            corner_angle(l_jk, l_lj, self.new_length).unwrap_or(0.0),
            corner_angle(l_il, l_ki, self.new_length).unwrap_or(0.0),
        ]
    }
}

/// Checks that edge `e` (in range) of a triangulation with this connectivity and these
/// lengths can be flipped, and returns its new length: every refusal of
/// [`IntrinsicTriangulation::flip_edge`] (the topological ones of
/// [`Complex::flip_refusal`], then strict convexity and the shape limit of SPEC N5), decided
/// from the lengths alone, for every representation (SPEC N50). It is [`plan_flip_with`] with
/// [`FlipShapeTest::N5`].
///
/// # Errors
/// The [`FlipRefusal`].
#[inline]
pub fn plan_flip(complex: &Complex, lengths: &[f64], e: usize) -> Result<FlipPlan, FlipRefusal> {
    plan_flip_with(complex, lengths, e, FlipShapeTest::N5)
}

/// [`plan_flip`] with the shape test chosen by the caller (SPEC N52): the topological
/// refusals, strict convexity, `test`, and the triangle inequalities of both new triangles.
///
/// # Errors
/// The [`FlipRefusal`]; a flip that `test` refuses is [`FlipRefusal::TooThin`] with the shape
/// number of N5 whichever test refused it.
pub fn plan_flip_with(
    complex: &Complex,
    lengths: &[f64],
    e: usize,
    test: FlipShapeTest,
) -> Result<FlipPlan, FlipRefusal> {
    if let Some(reason) = complex.flip_refusal(e) {
        return Err(reason);
    }
    let h0 = complex.edge_halfedge[e] as usize;
    let h1 = complex.twin[h0] as usize;
    let len = |h: usize| lengths[complex.halfedge_edge[h] as usize];
    let quad = [
        lengths[e],
        len(next(h0)),
        len(prev(h0)),
        len(next(h1)),
        len(prev(h1)),
    ];
    let [l_ij, l_jk, l_ki, l_il, l_lj] = quad;
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
    if !(test.accepts(area_at_tail, area_at_head, shape)
        && triangle_is_valid(l_il, new_length, l_ki)
        && triangle_is_valid(l_jk, new_length, l_lj))
    {
        return Err(FlipRefusal::TooThin { new_length, shape });
    }
    Ok(FlipPlan {
        lengths: quad,
        new_length,
    })
}

impl IntrinsicTriangulation {
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
    #[inline]
    pub fn flip_edge(&mut self, e: usize) -> Result<Flip, IntrinsicError> {
        self.flip_edge_with(e, FlipShapeTest::N5)
    }

    /// [`flip_edge`](Self::flip_edge) with the shape test chosen by the caller (SPEC N52)
    /// instead of that of N5.
    ///
    /// # Errors
    /// As [`flip_edge`](Self::flip_edge).
    pub fn flip_edge_with(
        &mut self,
        e: usize,
        test: FlipShapeTest,
    ) -> Result<Flip, IntrinsicError> {
        self.check_edge(e)?;
        self.try_flip(e, test).map_err(|reason| {
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
    #[inline]
    pub fn can_flip(&self, e: usize) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        self.can_flip_with(e, FlipShapeTest::N5)
    }

    /// Whether [`flip_edge_with`](Self::flip_edge_with) would flip edge `e` under `test`,
    /// without changing anything.
    ///
    /// # Errors
    /// As [`can_flip`](Self::can_flip).
    pub fn can_flip_with(
        &self,
        e: usize,
        test: FlipShapeTest,
    ) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        self.check_edge(e)?;
        Ok(self.plan_flip(e, test).map(|_| ()))
    }

    /// [`plan_flip_with`] on this triangulation.
    fn plan_flip(&self, e: usize, test: FlipShapeTest) -> Result<FlipPlan, FlipRefusal> {
        plan_flip_with(&self.complex, &self.edge_lengths, e, test)
    }

    /// [`flip_edge_with`](Self::flip_edge_with) for an edge id that is known to be in range.
    pub(crate) fn try_flip(&mut self, e: usize, test: FlipShapeTest) -> Result<Flip, FlipRefusal> {
        let FlipPlan {
            lengths: [l_ij, l_jk, l_ki, l_il, l_lj],
            new_length,
        } = self.plan_flip(e, test)?;
        let h0 = self.complex.edge_halfedge[e] as usize;
        let h1 = self.complex.twin[h0] as usize;
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
            new_corner_angles: [at_l, at_k],
            new_corners: [h1p as u32, h0p as u32],
        })
    }
}
