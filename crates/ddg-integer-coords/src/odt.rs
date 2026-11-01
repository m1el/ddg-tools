//! The intrinsic optimal Delaunay triangulation on integer coordinates (SPEC §4.11, I29;
//! IT-SPEC N21–N25: Sharp, Soliman, Crane 2019, Sec. 4.3, after Chen and Holst 2011).
//!
//! The iteration is `ddg_intrinsic`'s, written once over [`OdtMoves`]
//! ([`ddg_intrinsic::odt::optimal_delaunay`]): the same splits (longest first, at the
//! midpoint), the same Delaunay flips, the same Gauss–Seidel order of the moves and the same
//! halving of a step. What this representation supplies:
//!
//! * **The frame of a vertex** is [`AngleFrames`] (SPEC I9: direction `0` along
//!   `Complex::outgoing_halfedge`, the others from the corner angles); the update (the
//!   area-weighted mean of the vectors to the circumcentres, barycentres at the boundary) is
//!   [`odt_target_of`] in that frame, the signposts' formula.
//! * **The test of a step** is the signposts' ([`check_star_move`](ddg_intrinsic::reposition::check_star_move)): the moved star must keep
//!   every triangle's area positive and its shape number above the bound, or the step is
//!   halved. The integers could move a vertex farther (Sec. 3.6 is not limited to the star);
//!   the same test keeps the iteration the signposts' one.
//! * **The move** is Sec. 3.6 ([`IntegerTriangulation::move_vertex`]): trace the new place over
//!   `T¹`, insert a vertex there, remove the old one. Its connectivity afterwards differs from
//!   the signposts' (which keep the star and change only its lengths) until the iteration's
//!   Delaunay flips, which make both the Delaunay triangulation of the same points. A removal
//!   that gets stuck (I27: mostly the straight star of a vertex made by an edge split, which
//!   is what the iteration's splits make) refuses the move with
//!   [`IntrinsicError::RemovalRefused`], and the iteration halves the step as for a move that
//!   leaves the star. The signposts never refuse a move for that reason, so the two
//!   representations can end with different vertices (measured in SPEC §4.11).
//!
//! Inserted vertices on the boundary are not moved (I7): with
//! [`BoundaryVertices::AlongBoundary`](ddg_intrinsic::BoundaryVertices::AlongBoundary) the
//! iteration fails with [`IntrinsicError::Unsupported`]; the default keeps them fixed.

use ddg_intrinsic::odt::{OdtMoves, optimal_delaunay};
use ddg_intrinsic::reposition::{VertexMove, odt_target_of, star_min_shape_of};
use ddg_intrinsic::{
    AngleFrames, Exportable, IntrinsicError, OptimalDelaunayOptions, OptimalDelaunayReport,
};

use crate::certify::CertifyReport;
use crate::coords::IntegerTriangulation;
use crate::error::IntegerError;
use crate::interface::to_intrinsic_error;

/// Options of [`IntegerTriangulation::optimal_delaunay`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntegerOdtOptions {
    /// The options of the shared iteration (target length, iterations, step, halvings, shape
    /// bound, boundary vertices, centres, insertion limit, Delaunay flips).
    pub odt: OptimalDelaunayOptions,
    /// Whether to run the certificate (SPEC §2.4) on the result (I4: default yes).
    pub certify: bool,
}

impl Default for IntegerOdtOptions {
    fn default() -> Self {
        IntegerOdtOptions {
            odt: OptimalDelaunayOptions::default(),
            certify: true,
        }
    }
}

/// What [`IntegerTriangulation::optimal_delaunay`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct IntegerOdtReport {
    /// The shared iteration's report.
    pub odt: OptimalDelaunayReport,
    /// The certificate of the result (`None` if [`IntegerOdtOptions::certify`] is off).
    pub certificate: Option<CertifyReport>,
}

impl OdtMoves for IntegerTriangulation {
    /// The export's split ([`Exportable::split_for_export`]): the midpoint, snapped onto the
    /// crossing when it is within rounding of one (SPEC I28), as the signposts do.
    fn split_for_odt(&mut self, e: usize) -> Result<(), IntrinsicError> {
        Exportable::split_for_export(self, e).map(|_| ())
    }

    fn odt_target(
        &self,
        v: usize,
        barycenters_at_boundary: bool,
    ) -> Result<[f64; 2], IntrinsicError> {
        let frames = AngleFrames::of(self);
        odt_target_of(self, v, barycenters_at_boundary, &|h| frames.signpost(h))
    }

    fn star_min_shape(&self, v: usize) -> Option<f64> {
        star_min_shape_of(self, v)
    }

    /// [`IntegerTriangulation::move_vertex_in_star`]: the signposts' local update (the star
    /// keeps its triangles) on integer coordinates. [`VertexMove::back_direction`] is not
    /// computed (`NaN`) and [`VertexMove::snapped`] is `None` (a new place within rounding of an
    /// input edge is refused, the step halved).
    fn move_vertex_with_shape(
        &mut self,
        v: usize,
        displacement: [f64; 2],
        min_shape: f64,
    ) -> Result<VertexMove, IntrinsicError> {
        self.move_vertex_in_star(v, displacement, min_shape)
            .map_err(|e| to_intrinsic_error(e, "move_vertex_in_star"))
    }
}

impl IntegerTriangulation {
    /// The vector from the inserted vertex `v` to the target of its optimal Delaunay update,
    /// in the [`AngleFrames`] frame of `v` that [`move_vertex`](Self::move_vertex) takes
    /// (the signposts' `odt_target` in this frame).
    ///
    /// # Errors
    /// As `ddg_intrinsic::IntrinsicTriangulation::odt_target`.
    pub fn odt_target(
        &self,
        v: usize,
        barycenters_at_boundary: bool,
    ) -> Result<[f64; 2], IntrinsicError> {
        OdtMoves::odt_target(self, v, barycenters_at_boundary)
    }

    /// The intrinsic optimal Delaunay triangulation (see the [module documentation](self)):
    /// the shared iteration on this representation, then the certificate.
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] with the errors of
    /// [`ddg_intrinsic::odt::optimal_delaunay`] (invalid options, constrained edges, boundary
    /// moves, this crate's errors passed as [`IntrinsicError::Representation`]), the state as
    /// the iteration left it (valid); the certificate's error if the result does not certify.
    pub fn optimal_delaunay(
        &mut self,
        options: &IntegerOdtOptions,
    ) -> Result<IntegerOdtReport, IntegerError> {
        let odt = optimal_delaunay(self, &options.odt).map_err(IntegerError::Intrinsic)?;
        let certificate = if options.certify {
            Some(self.certify()?)
        } else {
            None
        };
        Ok(IntegerOdtReport { odt, certificate })
    }
}
