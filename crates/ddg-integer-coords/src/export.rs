//! The triangulation as a simplicial mesh, and the way back (SPEC §4.10, I29): `ddg_intrinsic`'s
//! export ([`ddg_intrinsic::export::export_simplicial`], written once over
//! [`Exportable`]) on integer coordinates.
//!
//! * **The export** splits and flips by the shared interface, so it makes the same splits and
//!   flips as on signposts wherever the lengths agree. The places of the vertices are the
//!   stored `q⁰` ([`Triangulation::vertex_location`](ddg_intrinsic::Triangulation::vertex_location)):
//!   for a vertex an edge split made, interpolated between the crossings around it (§4.5); for
//!   one a refinement inserted in a face, App. A (§4.4).
//! * **Refinement between rounds of splits** ([`ExportOptions::refine`]) is
//!   [`IntegerTriangulation::refine`] with the signposts' rules
//!   ([`RefineRules::Ruppert`], IT-SPEC N1, N4) and the certificate, so that the export is the
//!   same function of the triangulation in both representations.
//! * **The fingerprint** ([`SimplicialExport::fingerprint`]) is
//!   [`IntegerTriangulation::fingerprint`]: FNV-1a over the faces, the length of every
//!   halfedge's edge, the normal coordinates, the roundabouts and the places of the vertices.
//!   It differs from a signpost export's fingerprint of the same mesh (which hashes no integers),
//!   so a state made on one representation is refused by the other.
//! * **The way back** is [`IntegerTriangulation::segments_on_input`] ([`crate::segments`]).

use ddg_intrinsic::{
    EdgeTraceOptions, ExportOptions, Exportable, IntrinsicError, PathsOnInput, RefineOptions,
    RefineReport, RefineRules, Segment, SimplicialExport, SurfacePoint,
};

use crate::coords::IntegerTriangulation;
use crate::error::IntegerError;
use crate::interface::to_intrinsic_error;
use crate::refine::IntegerRefineOptions;

/// The parameters tried, in order, when the midpoint of an edge is still refused as
/// [`NotLocated`](ddg_intrinsic::InsertRefusal::NotLocated) (SPEC I28): `½ ± 2⁻¹⁶`, then
/// `½ ± 2⁻¹²`. Since I28 a midpoint within rounding of a crossing is not refused (the split
/// snaps onto the crossing, as the signposts do), so this is left for a midpoint within rounding
/// of an input edge that the edge does not cross there, or of a corner of the input face: none
/// in the evaluation.
pub const OFF_MIDPOINT: [f64; 4] = [
    0.5 + 1.0 / 65536.0,
    0.5 - 1.0 / 65536.0,
    0.5 + 1.0 / 4096.0,
    0.5 - 1.0 / 4096.0,
];

impl Exportable for IntegerTriangulation {
    /// The midpoint (SPEC I28): a midpoint within rounding of a crossing of the edge with an
    /// input edge, as on symmetric meshes where an edge crosses an input edge exactly at its
    /// middle (torus 5×3), is split at that crossing (`ic-split-at-crossing`), the vertex on the
    /// input edge, which is what the signposts do. Only if the midpoint is still refused with
    /// `InsertRefusal::NotLocated`, the first of [`OFF_MIDPOINT`] that is not (`Ok(true)`,
    /// counted in `SimplicialExport::n_off_midpoint_splits`; 0 in the evaluation).
    fn split_for_export(&mut self, e: usize) -> Result<bool, IntrinsicError> {
        let near_crossing = |err: &IntegerError| {
            matches!(
                err,
                IntegerError::Intrinsic(IntrinsicError::InsertionRefused {
                    reason: ddg_intrinsic::InsertRefusal::NotLocated { .. },
                    ..
                })
            )
        };
        let mut last = match self.split_edge(e, 0.5) {
            Ok(_) => return Ok(false),
            Err(err) if near_crossing(&err) && self.normal[e] >= 0 => err,
            Err(err) => return Err(to_intrinsic_error(err, "split_edge")),
        };
        for t in OFF_MIDPOINT {
            match self.split_edge(e, t) {
                Ok(_) => return Ok(true),
                Err(err) if near_crossing(&err) => last = err,
                Err(err) => return Err(to_intrinsic_error(err, "split_edge")),
            }
        }
        Err(to_intrinsic_error(last, "split_edge"))
    }

    fn refine_for_export(
        &mut self,
        options: &RefineOptions,
    ) -> Result<RefineReport, IntrinsicError> {
        let options = IntegerRefineOptions {
            refine: *options,
            rules: RefineRules::Ruppert,
            ..IntegerRefineOptions::default()
        };
        self.refine(&options)
            .map(|report| report.refine)
            .map_err(|e| to_intrinsic_error(e, "refine"))
    }

    fn state_fingerprint(&self) -> u64 {
        self.fingerprint()
    }

    fn segments_on_input(
        &self,
        nodes: &[SurfacePoint],
        segments: &[Segment],
        options: &EdgeTraceOptions,
    ) -> Result<PathsOnInput, IntrinsicError> {
        IntegerTriangulation::segments_on_input(self, nodes, segments, options)
            .map_err(|e| to_intrinsic_error(e, "segments_on_input"))
    }
}

impl IntegerTriangulation {
    /// Makes the triangulation a simplicial complex by splitting edges and hands it out as a
    /// `TriMesh` with a `Geometry`: `ddg_intrinsic`'s
    /// [`export_simplicial`](ddg_intrinsic::export::export_simplicial) on this representation
    /// (see the [module documentation](self) for what differs). The triangulation is changed
    /// by this; the way back ([`segments_on_input`](Self::segments_on_input)) is asked of it in
    /// that state.
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] with the export's errors
    /// (`ddg_intrinsic::IntrinsicError::ExportFailed`, …); this crate's errors of the splits,
    /// flips and refinement passed as `IntrinsicError::Representation`.
    pub fn export_simplicial(
        &mut self,
        options: &ExportOptions,
    ) -> Result<SimplicialExport, IntegerError> {
        ddg_intrinsic::export::export_simplicial(self, options).map_err(IntegerError::Intrinsic)
    }
}
