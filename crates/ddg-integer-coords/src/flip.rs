//! The edge flip (paper Sec. 3.2, Eq. 3, 4; SPEC §4.2).
//!
//! The refusals and the new length come from [`ddg_intrinsic::plan_flip`], so the integer and
//! the signpost representation refuse the same flips and give the flipped edge bitwise the same
//! length (SPEC §2.1, §4.3). The integers of the new edge are computed from the state before the
//! flip; only then is the connectivity changed ([`Complex::flip_edge`]), so a flip that fails
//! changes nothing.
//!
//! With the flipped edge `ij`, its faces `ijk` (on the side of its halfedge) and `jil`, and the
//! new edge `kl` (SPEC §4.2; [`ComplexFlip`] for the ids):
//!
//! * `n_kl` by Eq. 3 ([`formulas::flip_normal`]) with the corrected first terms
//!   `c^{ji}_l + c^{ij}_k`; a value below `−1` (two or more input edges from `k` to `l` inside the
//!   quadrilateral, which straight edges cannot do) is [`IntegerError::InvalidCoordinates`]
//!   naming the edge; a value outside `i64` is [`IntegerError::Overflow`];
//! * `r_kl` and `r_lk` by Eq. 4 ([`formulas::roundabout_update`], the update of the roundabouts,
//!   `ic-roundabout-update`, SPEC §9.2 tier C) for the ends in `V⁰`, from the roundabouts of
//!   `k → i` and `l → j`, the emanating counts of the corners `k` of `kil` and `l` of `ljk`, and
//!   `n⁻_ki`, `n⁻_lj`; [`NO_ROUNDABOUT`] for ends in `V*`;
//! * the halfedges `i → l` and `j → k`, which change faces, keep their roundabouts
//!   ([`ComplexFlip::move_halfedge_data`]); nothing else changes.
//!
//! [`Complex::flip_edge`]: ddg_intrinsic::Complex::flip_edge
//! [`ComplexFlip`]: ddg_intrinsic::complex::ComplexFlip
//! [`ComplexFlip::move_halfedge_data`]: ddg_intrinsic::complex::ComplexFlip::move_halfedge_data

use ddg_intrinsic::{
    DelaunayOptions, DelaunayReport, Flip, FlipRefusal, IntrinsicError, Triangulation, plan_flip,
};

use crate::certify::CertifyReport;
use crate::coords::{IntegerTriangulation, NO_ROUNDABOUT};
use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::{self, FlipNormalError, emanating, n_minus};

/// The integers a flip gives the new edge, computed before anything changes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FlippedIntegers {
    /// `n_kl`.
    pub(crate) normal: i64,
    /// The roundabout of `l → k` (it goes to `ComplexFlip::new_halfedges[0]`).
    pub(crate) r_lk: u32,
    /// The roundabout of `k → l` (`ComplexFlip::new_halfedges[1]`).
    pub(crate) r_kl: u32,
}

impl IntegerTriangulation {
    /// Whether [`flip_edge`](Self::flip_edge) would flip edge `e`, without changing anything:
    /// the refusals of [`plan_flip`] (boundary, constrained, degree-one end point, not strictly
    /// convex, too thin), decided from the lengths alone as in the signpost representation.
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; a refusal is the `Err` of the inner result.
    pub fn can_flip(&self, e: usize) -> Result<Result<(), FlipRefusal>, IntegerError> {
        self.check_edge(e)?;
        Ok(plan_flip(&self.complex, &self.lengths, e).map(|_| ()))
    }

    /// Flips edge `e` (paper Sec. 3.2): its faces `ijk`, `jil` become `ilk`, `jkl` (ids as
    /// [`Complex::flip_edge`](ddg_intrinsic::Complex::flip_edge)), the new edge `kl` gets the
    /// length of `it-flip-length` (from [`plan_flip`], bitwise as the signpost
    /// representation), the normal coordinate of Eq. 3 and, at its ends in `V⁰`, the
    /// roundabouts of Eq. 4 (see the [module documentation](crate::flip)). Every other edge
    /// keeps its id, length, direction and normal coordinate, and every other halfedge its
    /// roundabout.
    ///
    /// Flipping the new edge again restores the integers exactly (`ic-flip-involution`) and the
    /// length up to rounding.
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::Intrinsic`] with
    /// [`IntrinsicError::FlipRefused`] and the [`FlipRefusal`]; [`IntegerError::Overflow`]
    /// naming the edge when `n_kl` leaves `i64`; [`IntegerError::InvalidCoordinates`] for a
    /// face of the quadrilateral that is not valid, an `n_kl < −1`, or a roundabout outside its
    /// range. Nothing is changed in every case.
    ///
    /// [formula: ic-flip-normal]
    pub fn flip_edge(&mut self, e: usize) -> Result<Flip, IntegerError> {
        self.check_edge(e)?;
        let plan = plan_flip(&self.complex, &self.lengths, e).map_err(|reason| {
            let [a, b] = self.complex.edge_vertices(e);
            IntegerError::Intrinsic(IntrinsicError::FlipRefused {
                edge: e as u32,
                vertices: [a as u32, b as u32],
                reason,
            })
        })?;
        let integers = self.flipped_integers(e)?;
        let flip = self.complex.flip_edge(e)?;
        flip.move_halfedge_data(&mut self.roundabouts);
        let [h_lk, h_kl] = flip.new_halfedges.map(|h| h as usize);
        self.roundabouts[h_lk] = integers.r_lk;
        self.roundabouts[h_kl] = integers.r_kl;
        self.normal[e] = integers.normal;
        self.lengths[e] = plan.new_length;
        Ok(Flip {
            edge: e as u32,
            old_vertices: flip.old_vertices,
            new_vertices: flip.new_vertices,
            old_length: plan.lengths[0],
            new_length: plan.new_length,
        })
    }

    /// Delaunay flipping (paper Sec. 4.1, SPEC §4.3) followed by the certificate (SPEC §2.4, I4:
    /// certify at the end of every operation that does a variable amount of work). The
    /// flipping is the shared [`Triangulation::flip_to_delaunay`], so it makes the same flips in
    /// the same order as the signpost representation on the same input; a flip that fails with
    /// an error of this crate (an overflow, say) is counted as refused there
    /// (`DelaunayReport::n_refused`) and changes nothing.
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] for invalid options (nothing changed); the error of
    /// [`certify`](Self::certify) if the result does not certify (the flips are then done:
    /// the state is not that of any realization, which would be a bug of this crate).
    pub fn flip_to_delaunay_certified(
        &mut self,
        options: &DelaunayOptions,
    ) -> Result<(DelaunayReport, CertifyReport), IntegerError> {
        let report = Triangulation::flip_to_delaunay(self, options)?;
        let certificate = self.certify()?;
        Ok((report, certificate))
    }

    /// Eq. 3 and 4 for flipping the interior edge `e` (whose two faces differ), from the current
    /// state, without changing it.
    pub(crate) fn flipped_integers(&self, e: usize) -> Result<FlippedIntegers, IntegerError> {
        let c = &self.complex;
        let h0 = c.edge_halfedge(e);
        let h1 = c
            .twin(h0)
            .expect("plan_flip refuses boundary edges before this is called");
        // f0 = (i, j, k): h0 = i→j, h0n = j→k, h0p = k→i; f1 = (j, i, l): h1 = j→i,
        // h1n = i→l, h1p = l→j.
        let (h0n, h0p, h1n, h1p) = (c.next(h0), c.prev(h0), c.next(h1), c.prev(h1));
        let n = |h: usize| self.normal[c.edge(h)];
        let (n_ij, n_jk, n_ki, n_il, n_lj) = (self.normal[e], n(h0n), n(h0p), n(h1n), n(h1p));
        let invalid =
            |item, detail: String| IntegerError::invalid(item, Place::Edge(e as u32), detail);
        let n_kl = match formulas::flip_normal(n_ij, n_jk, n_ki, n_il, n_lj) {
            Ok(x) => x,
            Err(FlipNormalError::Overflow) => {
                return Err(IntegerError::Overflow { edge: e as u32 });
            }
            Err(FlipNormalError::InvalidFace) => {
                return Err(invalid(
                    ValidityItem::FaceParity,
                    format!(
                        "a face of the flip's quadrilateral is not valid (n_ij {n_ij}, n_jk \
                         {n_jk}, n_ki {n_ki}, n_il {n_il}, n_lj {n_lj})"
                    ),
                ));
            }
            Err(FlipNormalError::OddDoubled(d)) => {
                return Err(invalid(
                    ValidityItem::FaceParity,
                    format!("twice Eq. 3 is odd ({d}); it is not halved (I14)"),
                ));
            }
        };
        if n_kl < -1 {
            return Err(invalid(
                ValidityItem::Ranges,
                format!(
                    "Eq. 3 gives n_kl = {n_kl} < −1 (n_ij {n_ij}, n_jk {n_jk}, n_ki {n_ki}, n_il \
                     {n_il}, n_lj {n_lj}): two or more input edges from k to l inside the \
                     quadrilateral, which straight edges cannot do"
                ),
            ));
        }
        // New faces: kil (corner k emanates to il) and ljk (corner l emanates to jk).
        let e_k = emanating(n_il, n_ki, n_kl);
        let e_l = emanating(n_jk, n_lj, n_kl);
        let k = c.from_vertex(h0p);
        let l = c.from_vertex(h1p);
        let r_kl = self.updated_roundabout(e, k, h0p, e_k, n_minus(n_ki))?;
        let r_lk = self.updated_roundabout(e, l, h1p, e_l, n_minus(n_lj))?;
        Ok(FlippedIntegers {
            normal: n_kl,
            r_lk,
            r_kl,
        })
    }

    /// Eq. 4 at vertex `v` from the roundabout of the halfedge `before` (just clockwise of the
    /// new one), for the new edge `e`: [`NO_ROUNDABOUT`] if `v ∈ V*`.
    fn updated_roundabout(
        &self,
        e: usize,
        v: usize,
        before: usize,
        e_new: i64,
        n_minus_before: i64,
    ) -> Result<u32, IntegerError> {
        if !self.is_input_vertex(v) {
            return Ok(NO_ROUNDABOUT);
        }
        let deg = self.local.degree(v);
        let bad = |detail: String| {
            IntegerError::invalid(ValidityItem::Ranges, Place::Halfedge(before as u32), detail)
        };
        let r = self.roundabouts[before];
        if r >= deg {
            return Err(bad(format!(
                "the roundabout {r} of halfedge {before} at input vertex {v} is not below \
                 deg₀ = {deg}"
            )));
        }
        let x = formulas::roundabout_update(r as i64, e_new, n_minus_before, self.local.modulus(v))
            .ok_or(IntegerError::Overflow { edge: e as u32 })?;
        if !(0..deg as i64).contains(&x) {
            return Err(bad(format!(
                "Eq. 4 gives the new edge {e} the roundabout {x} at boundary vertex {v}, not \
                 below deg₀ = {deg} (r {r}, e {e_new}, n⁻ {n_minus_before})"
            )));
        }
        Ok(x as u32)
    }

    /// [`IntegerError::IndexOutOfRange`] unless `e` is an edge.
    pub(crate) fn check_edge(&self, e: usize) -> Result<(), IntegerError> {
        let bound = self.complex.n_edges();
        if e < bound {
            Ok(())
        } else {
            Err(IntegerError::IndexOutOfRange {
                what: "edge",
                index: e,
                bound,
            })
        }
    }
}
