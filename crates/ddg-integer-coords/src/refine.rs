//! Intrinsic Delaunay refinement on integer coordinates (paper Sec. 4.2, Supp. Alg. 10,
//! App. C; SPEC §5, I11, I12, I23).
//!
//! The driver is `ddg_intrinsic`'s ([`ddg_intrinsic::refine::refine`]), written once over the
//! shared interface [`Triangulation`]: the same queue, the same circumcentres, the same
//! boundary rule and the same flips as the signpost representation (I11), so that both insert
//! the same points on the same input wherever their floating-point decisions agree. What this
//! crate supplies is the representation: Delaunay flips by Eq. 3–4, the circumcentre traced over
//! `T¹` with angle frames (I9; from the barycentre the path depends on the lengths alone), face
//! splits by Supp. Alg. 5–6 with App. A, edge splits by Alg. 7 / I20, and removals by Alg. 8
//! with the flips of I24.
//!
//! **Rules.** [`IntegerRefineOptions::rules`] defaults to [`RefineRules::Chew`], the paper's:
//! a boundary edge is split at its midpoint only when a circumcentre path leaves the surface
//! through it, then the inserted interior vertices within Dijkstra distance `ℓ` of the new
//! vertex are removed (Xia 2013, Cor. 1; [formula: ic-dijkstra-ball]) and the triangulation is
//! flipped to Delaunay; the narrow-cone rule of I23. [`RefineRules::Ruppert`] runs the signpost
//! crate's default rules (IT-SPEC N1, N4) on integer coordinates. Termination is the paper's
//! Thm. C.1 for closed surfaces with angle sums of at least 60° and bounds up to 30°
//! ([formula: ic-refinement-termination]); with boundary it is not proved, and the insertion
//! limit holds in every case.
//!
//! **What floating point decides here** (SPEC §3.2): the circumcentre's face and barycentric
//! coordinates, the region of a face split (I6: the distance to the nearest arc is reported,
//! [`IntegerRefineReport::min_curve_distance`]), the crossings before an edge split point, and
//! the places `q⁰` (App. A). A face split whose App. A place is outside its input face by more
//! than rounding ([`IntegerError::RecoverFailed`]) is refused like a refused insertion: the face
//! is left and counted ([`IntegerRefineReport::n_recover_failed`]). Crossing parameters clamped
//! into `[0, 1]` are counted (I5). The certificate (§2.4) runs at the end (I4).

use ddg_intrinsic::refine::refine;
use ddg_intrinsic::{
    Complex, EdgeTrace, EdgeTraceOptions, Flip, FlipRefusal, InputSurface, InsertRefusal,
    InsertedVertexRemoval, Insertion, IntrinsicError, RefineOptions, RefineReport, RefineRules,
    SurfacePoint, Triangulation,
};

use crate::certify::CertifyReport;
use crate::coords::IntegerTriangulation;
use crate::error::IntegerError;
use crate::insert::FaceSplitOptions;
use crate::interface::to_intrinsic_error;

/// The distance to the nearest arc (in shortest sides of the face, SPEC I6) below which a face
/// split counts as near an arc in [`IntegerRefineReport::n_near_arc`]: where the signposts may
/// put the point on the other side of the arc by rounding (WP4 measured such disagreements up
/// to `1.2e-10`).
pub const NEAR_ARC: f64 = 1e-9;

/// Options of [`IntegerTriangulation::refine`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntegerRefineOptions {
    /// The options of the shared driver: angle bound (default 25°, the paper's evaluation),
    /// circumradius bound, insertion limit (default 10 per vertex, at least 1000), Delaunay
    /// test.
    pub refine: RefineOptions,
    /// The rules (default [`RefineRules::Chew`], the paper's; see the
    /// [module documentation](self)).
    pub rules: RefineRules,
    /// The options of the face splits (SPEC I6). Default: never refused for being near an arc;
    /// the integers are valid whichever side of the arc rounding puts the point.
    pub face_split: FaceSplitOptions,
    /// Whether to run the certificate (SPEC §2.4) on the result (I4: default yes).
    pub certify: bool,
}

impl Default for IntegerRefineOptions {
    fn default() -> Self {
        IntegerRefineOptions {
            refine: RefineOptions::default(),
            rules: RefineRules::Chew,
            face_split: FaceSplitOptions::default(),
            certify: true,
        }
    }
}

/// What [`IntegerTriangulation::refine`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct IntegerRefineReport {
    /// The shared driver's report: insertions, splits, flips, removals and refused removals,
    /// trace failures, the angles and what is left below the bound, the narrow-cone faces.
    pub refine: RefineReport,
    /// Face splits made (circumcentres inserted in a face).
    pub n_face_splits: usize,
    /// Face splits within [`NEAR_ARC`] of an arc (SPEC I6).
    pub n_near_arc: usize,
    /// The smallest distance of a face split to the nearest arc, in shortest sides of its face
    /// (`+∞` if no face split had an arc).
    pub min_curve_distance: f64,
    /// Face splits refused because App. A placed the point outside its input face by more than
    /// rounding ([`IntegerError::RecoverFailed`]) or, with a positive
    /// [`FaceSplitOptions::too_close`], because it was too near an arc; the face was left.
    pub n_recover_failed: usize,
    /// Places `q⁰` within rounding outside their input face that were moved inside (I16).
    pub n_places_moved: usize,
    /// Crossings whose parameters were clamped into `[0, 1]` in the face and edge splits
    /// (SPEC I5).
    pub clamped_crossings: usize,
    /// The largest such clamp.
    pub max_clamp: f64,
    /// The certificate of the result (`None` if [`IntegerRefineOptions::certify`] is off).
    pub certificate: Option<CertifyReport>,
}

/// An [`IntegerTriangulation`] seen by the shared driver: every operation is the
/// triangulation's own; face and edge splits record what floating point decided.
struct Recording<'a> {
    tri: &'a mut IntegerTriangulation,
    face_split: FaceSplitOptions,
    stats: IntegerRefineReport,
    /// The last error of this crate passed to the driver as
    /// [`IntrinsicError::Representation`], to return it as itself.
    last_error: Option<IntegerError>,
}

impl Recording<'_> {
    fn pass(&mut self, error: IntegerError, operation: &'static str) -> IntrinsicError {
        self.last_error = Some(error.clone());
        to_intrinsic_error(error, operation)
    }

    fn clamps(&mut self, n: usize, max: f64) {
        self.stats.clamped_crossings += n;
        self.stats.max_clamp = self.stats.max_clamp.max(max);
    }
}

impl Triangulation for Recording<'_> {
    #[inline]
    fn complex(&self) -> &Complex {
        &self.tri.complex
    }
    #[inline]
    fn input(&self) -> &InputSurface {
        &self.tri.input
    }
    #[inline]
    fn edge_lengths(&self) -> &[f64] {
        &self.tri.lengths
    }
    #[inline]
    fn vertex_angle_sum(&self, v: usize) -> f64 {
        self.tri.vertex_angle_sum(v)
    }
    #[inline]
    fn vertex_location(&self, v: usize) -> SurfacePoint {
        self.tri.locations[v]
    }
    fn can_flip(&self, e: usize) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        Triangulation::can_flip(&*self.tri, e)
    }
    fn flip_edge(&mut self, e: usize) -> Result<Flip, IntrinsicError> {
        match self.tri.flip_edge(e) {
            Ok(flip) => Ok(flip),
            Err(err) => Err(self.pass(err, "flip_edge")),
        }
    }
    fn insert_vertex(&mut self, f: usize, bary: [f64; 3]) -> Result<Insertion, IntrinsicError> {
        match self.tri.split_face(f, bary, &self.face_split) {
            Ok(report) => {
                let location = &report.location;
                self.stats.n_face_splits += 1;
                if let Some(d) = location.choice.curve_distance {
                    self.stats.min_curve_distance = self.stats.min_curve_distance.min(d);
                    self.stats.n_near_arc += usize::from(d < NEAR_ARC);
                }
                self.stats.n_places_moved += usize::from(report.recover.moved > 0.0);
                let (n, max) = (
                    location.sides.clamped_crossings + location.parameter_clamps,
                    location.sides.max_clamp.max(location.max_parameter_clamp),
                );
                self.clamps(n, max);
                Ok(report.insertion)
            }
            Err(
                err @ (IntegerError::RecoverFailed { .. } | IntegerError::TooCloseToCurve { .. }),
            ) => {
                self.stats.n_recover_failed += 1;
                Err(IntrinsicError::InsertionRefused {
                    at: SurfacePoint::Face {
                        face: f as u32,
                        bary,
                    },
                    reason: InsertRefusal::NotLocated {
                        detail: err.to_string(),
                    },
                })
            }
            Err(err) => Err(self.pass(err, "insert_vertex")),
        }
    }
    fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntrinsicError> {
        if e < self.tri.complex.n_edges() && self.tri.normal[e] > 0 {
            // The crossings the split divides, for the clamps (I5); the split extracts them
            // again.
            match self
                .tri
                .extract_intrinsic_edge_with_report(e, &EdgeTraceOptions::default())
            {
                Ok((_, strip)) => self.clamps(strip.clamped_crossings, strip.max_clamp),
                Err(err) => return Err(self.pass(err, "split_edge")),
            }
        }
        match self.tri.split_edge(e, t) {
            Ok(insertion) => Ok(insertion),
            Err(err) => Err(self.pass(err, "split_edge")),
        }
    }
    fn set_constrained(&mut self, e: usize, constrained: bool) -> Result<(), IntrinsicError> {
        Triangulation::set_constrained(&mut *self.tri, e, constrained)
    }
    fn trace_edge_on_input(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, IntrinsicError> {
        Triangulation::trace_edge_on_input(&*self.tri, e, options)
    }
    fn remove_inserted_vertex(
        &mut self,
        v: usize,
    ) -> Result<InsertedVertexRemoval, IntrinsicError> {
        Triangulation::remove_inserted_vertex(&mut *self.tri, v)
    }
}

impl IntegerTriangulation {
    /// Intrinsic Delaunay refinement (paper Sec. 4.2, Supp. Alg. 10; SPEC §5): flips to
    /// Delaunay, then inserts circumcentres (and splits boundary edges, removing the inserted
    /// vertices near them, by the rules chosen) until no triangle needs refinement, then
    /// certifies the result. See the [module documentation](crate::refine).
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] with the errors of
    /// [`ddg_intrinsic::refine::refine`] (invalid options, constrained edges); the errors of
    /// this crate's operations other than refusals (e.g. [`IntegerError::Overflow`]), with the
    /// state as that operation left it (valid, partly refined); the certificate's error if
    /// the result does not certify (a bug). Reaching the insertion limit is not an error.
    pub fn refine(
        &mut self,
        options: &IntegerRefineOptions,
    ) -> Result<IntegerRefineReport, IntegerError> {
        self.refine_with(options, &|_, _| false)
    }

    /// [`refine`](Self::refine) with an extra criterion: a face `f` for which
    /// `extra(triangulation, f)` returns true is refined as well (the hook of adaptive
    /// refinement; see `ddg_intrinsic`'s `IntrinsicTriangulation::refine_with` for the
    /// conditions it must meet).
    ///
    /// # Errors
    /// As [`refine`](Self::refine).
    pub fn refine_with(
        &mut self,
        options: &IntegerRefineOptions,
        extra: &dyn Fn(&IntegerTriangulation, usize) -> bool,
    ) -> Result<IntegerRefineReport, IntegerError> {
        let mut recording = Recording {
            tri: self,
            face_split: options.face_split,
            stats: IntegerRefineReport {
                refine: placeholder_report(options.rules),
                n_face_splits: 0,
                n_near_arc: 0,
                min_curve_distance: f64::INFINITY,
                n_recover_failed: 0,
                n_places_moved: 0,
                clamped_crossings: 0,
                max_clamp: 0.0,
                certificate: None,
            },
            last_error: None,
        };
        let result = refine(
            &mut recording,
            &options.refine,
            options.rules,
            &|r: &Recording, f| extra(r.tri, f),
        );
        let Recording {
            mut stats,
            last_error,
            ..
        } = recording;
        stats.refine = match (result, last_error) {
            (Ok(report), _) => report,
            (Err(IntrinsicError::Representation { .. }), Some(original)) => return Err(original),
            (Err(err), _) => return Err(IntegerError::Intrinsic(err)),
        };
        if options.certify {
            stats.certificate = Some(self.certify()?);
        }
        Ok(stats)
    }
}

/// A report to fill in.
fn placeholder_report(rules: RefineRules) -> RefineReport {
    RefineReport {
        status: ddg_intrinsic::RefineStatus::Done,
        n_inserted: 0,
        n_face_insertions: 0,
        n_edge_splits: 0,
        n_boundary_splits: 0,
        n_flips: 0,
        n_snapped: 0,
        min_angle: f64::INFINITY,
        max_angle: 0.0,
        min_angle_reachable: f64::INFINITY,
        below_bound: Vec::new(),
        too_large: Vec::new(),
        n_not_delaunay: 0,
        n_encroached: 0,
        rules,
        n_removals: 0,
        n_removal_flips: 0,
        n_removals_refused: 0,
        n_trace_failures: 0,
        n_narrow_cone_faces: 0,
        min_angle_valid: f64::INFINITY,
    }
}
