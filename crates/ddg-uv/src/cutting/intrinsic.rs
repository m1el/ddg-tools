//! The intrinsic triangulation the cut flow can run on, and the state that goes with it.
//!
//! "Variational Surface Cutting" is sensitive to badly shaped triangles. With
//! [`Triangulation::Intrinsic`] the session builds, from the seam-cut triangle mesh it already
//! makes (so the seams stay boundary), an intrinsic triangulation of the same surface
//! (Sharp, Soliman, Crane, *Navigating Intrinsic Triangulations*, 2019; crate `ddg-intrinsic`):
//! it is flipped to Delaunay, optionally refined (vertices are inserted, also on boundary
//! edges) and exported as a simplicial mesh with a geometry. The flow runs on that export. Its
//! first vertices are the vertices of the seam-cut mesh; the others lie inside faces or on
//! edges of it.
//!
//! **The state** of such a session has one value per vertex of the export, which is more than
//! the polygon mesh has. [`SessionState`] holds it together with the fingerprint of the export
//! and the options that reproduce the triangulation (refinement and export are
//! deterministic). It can be written to bytes and read back, and is refused by a session whose
//! export has another fingerprint.

use std::time::Instant;

use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    BelowBound, DelaunayCriterion, DelaunayOptions, EdgeTraceOptions, ExportOptions, FlipLimit,
    InsertionLimit, IntrinsicTriangulation, RefineOptions, RefineReport, SimplicialExport,
    SurfacePoint,
};

use super::error::CuttingError;

/// The triangulation the cut flow runs on.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(clippy::large_enum_variant)]
pub enum Triangulation {
    /// The triangle mesh of the caller's triangulation, cut open at the seams. The behaviour of
    /// every earlier version, bit for bit ([`Session::new`](super::Session::new) and
    /// [`Session::new_with`](super::Session::new_with) use it).
    Caller,
    /// An intrinsic triangulation of that mesh, flipped to Delaunay and exported; see
    /// [`IntrinsicOptions`].
    Intrinsic(IntrinsicOptions),
}

/// Options of [`Triangulation::Intrinsic`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntrinsicOptions {
    /// The Delaunay flips made first.
    pub delaunay: DelaunayOptions,
    /// Refinement to an angle bound after the flips (`None`: only flips).
    pub refine: Option<RefineOptions>,
    /// The export to a simplicial mesh: splits of what is not simplicial, flips again, and
    /// optionally refinement between rounds of splits (see `ddg_intrinsic::ExportOptions`).
    pub export: ExportOptions,
    /// How accurately a cut is carried back to the caller's triangles (the tolerance and the
    /// crossing limit of `ddg_intrinsic::IntrinsicTriangulation::segments_on_input`).
    pub trace: EdgeTraceOptions,
    /// The smallest distance, as a fraction of a side or of the barycentric coordinates, that a
    /// point of the cut keeps from the vertices and sides of the caller's triangle it lies in.
    /// A point of an edge is moved into `[min_offset, 1 − min_offset]`; a point inside a
    /// triangle has every barycentric coordinate raised to `min_offset` and is renormalized.
    /// Moved points are counted in [`CutOnPolygons::n_clamped_crossings`](super::CutOnPolygons).
    /// This is the rule that keeps the parts of a polygon from degenerating when the cut passes
    /// arbitrarily close to a vertex or an edge.
    pub min_offset: f64,
}

impl IntrinsicOptions {
    /// Flips to Delaunay, refines to 25 degrees (the default of `ddg-intrinsic`), exports with
    /// the refinement kept between the rounds of splits. Refinement is an explicit choice.
    pub fn refined() -> IntrinsicOptions {
        IntrinsicOptions {
            delaunay: DelaunayOptions::default(),
            refine: Some(RefineOptions::default()),
            export: ExportOptions {
                refine: Some(RefineOptions::default()),
                ..ExportOptions::default()
            },
            trace: EdgeTraceOptions::default(),
            min_offset: 1e-6,
        }
    }

    /// Delaunay flips only: no vertex is inserted except what the export needs to make the
    /// triangulation a simplicial complex. The default. When the report says no vertex was
    /// added ([`IntrinsicReport::n_added`] is 0, [`Session::state_is_per_corner`](super::Session::state_is_per_corner)),
    /// the state has one value per vertex of the polygon mesh (per seam copy) and the per-corner
    /// view is the whole state.
    pub fn flips_only() -> IntrinsicOptions {
        IntrinsicOptions {
            refine: None,
            export: ExportOptions::default(),
            ..IntrinsicOptions::refined()
        }
    }

    pub(super) fn validate(&self) -> Result<(), CuttingError> {
        if !(self.min_offset.is_finite() && (0.0..0.5).contains(&self.min_offset)) {
            return Err(CuttingError::InvalidOption {
                name: "min_offset",
                reason: format!("must be in [0, 0.5), got {}", self.min_offset),
            });
        }
        Ok(())
    }
}

impl Default for IntrinsicOptions {
    /// [`IntrinsicOptions::flips_only`]: the measurements of `ddg-intrinsic` and of the cut
    /// flow on it show that flipping to Delaunay is what makes badly shaped meshes behave, and
    /// that refinement adds 30 % to 140 % vertices for little gain on most meshes.
    fn default() -> Self {
        IntrinsicOptions::flips_only()
    }
}

impl Default for Triangulation {
    /// Intrinsic, Delaunay flips only ([`IntrinsicOptions::default`]). Note that
    /// [`Session::new`](super::Session::new) keeps the caller's triangulation; ask for this one
    /// with [`Session::with_triangulation`](super::Session::with_triangulation).
    fn default() -> Self {
        Triangulation::Intrinsic(IntrinsicOptions::default())
    }
}

/// What the construction of the intrinsic triangulation did.
#[derive(Debug, Clone, PartialEq)]
pub struct IntrinsicReport {
    /// Vertices of the seam-cut mesh.
    pub n_input_vertices: usize,
    /// Vertices of the triangulation the flow runs on.
    pub n_vertices: usize,
    /// Vertices added (refinement and the export's splits): the difference.
    pub n_added: usize,
    /// Flips of the first Delaunay pass.
    pub n_initial_flips: usize,
    /// Vertices inserted by the refinement before the export.
    pub n_refined: usize,
    /// Edges the export split to make the triangulation a simplicial complex.
    pub n_splits: usize,
    /// Flips made by the export (its refinements included).
    pub n_export_flips: usize,
    /// Refinements the export made between rounds of splits.
    pub n_refine_rounds: usize,
    /// Whether the export is Delaunay.
    pub is_delaunay: bool,
    /// The smallest corner angle of the export.
    pub min_angle: f64,
    /// The largest corner angle of the export.
    pub max_angle: f64,
    /// The angle bound that was asked for (`None`: no refinement).
    pub bound: Option<f64>,
    /// Corners of the export below the bound, counted on the export itself.
    pub n_below_bound: usize,
    /// The corners the last refinement left below the bound, with the reason: those of the
    /// export's own refinement if the export ended with one (then they are the corners of the
    /// export), otherwise those of the refinement before the export, which splits may have
    /// changed (see `ddg_intrinsic::SimplicialExport::min_angle`).
    pub below_bound: Vec<BelowBound>,
    /// Seconds spent building the triangulation and the export.
    pub seconds: f64,
    /// The fingerprint of the export.
    pub fingerprint: u64,
}

/// The intrinsic triangulation of a session, in the state of its export.
#[derive(Debug, Clone)]
pub(super) struct IntrinsicState {
    pub options: IntrinsicOptions,
    pub triangulation: IntrinsicTriangulation,
    pub export: SimplicialExport,
    pub report: IntrinsicReport,
}

impl IntrinsicState {
    /// Builds the triangulation of `mesh` with `geometry`.
    pub fn build(
        mesh: &TriMesh,
        geometry: &Geometry,
        options: &IntrinsicOptions,
    ) -> Result<IntrinsicState, CuttingError> {
        options.validate()?;
        let start = Instant::now();
        let mut tri = IntrinsicTriangulation::new(mesh, geometry)?;
        let initial = tri.flip_to_delaunay(&options.delaunay)?;
        let refinement: Option<RefineReport> = match &options.refine {
            Some(r) => Some(tri.refine(r)?),
            None => None,
        };
        let n_refined = refinement.as_ref().map_or(0, |r| r.n_inserted);
        let export = tri.export_simplicial(&options.export)?;
        let bound = options
            .refine
            .map(|r| r.min_angle)
            .or(options.export.refine.map(|r| r.min_angle));
        let n_below_bound = bound.map_or(0, |b| {
            export
                .geometry
                .corner_angles()
                .iter()
                .flatten()
                .filter(|&&a| a < b)
                .count()
        });
        let below_bound = match (&export.refinement, &refinement) {
            (Some(r), _) => r.below_bound.clone(),
            (None, Some(r)) => r.below_bound.clone(),
            _ => Vec::new(),
        };
        let n_input_vertices = mesh.n_vertices();
        let report = IntrinsicReport {
            n_input_vertices,
            n_vertices: export.mesh.n_vertices(),
            n_added: export.mesh.n_vertices() - n_input_vertices,
            n_initial_flips: initial.n_flips,
            n_refined,
            n_splits: export.n_splits,
            n_export_flips: export.n_flips,
            n_refine_rounds: export.n_refine_rounds,
            is_delaunay: export.is_delaunay,
            min_angle: export.min_angle,
            max_angle: export.max_angle,
            bound,
            n_below_bound,
            below_bound,
            seconds: start.elapsed().as_secs_f64(),
            fingerprint: export.fingerprint,
        };
        Ok(IntrinsicState {
            options: *options,
            triangulation: tri,
            export,
            report,
        })
    }

    /// Values on the vertices of the export from values on the vertices of the seam-cut mesh
    /// `mesh`: the value of an inserted vertex is the linear interpolation of the values of the
    /// mesh at the vertex's place on it (the piecewise linear function of the mesh evaluated
    /// there).
    pub fn extend(&self, mesh: &TriMesh, values: &[f64]) -> Vec<f64> {
        self.export
            .locations
            .iter()
            .map(|location| match *location {
                SurfacePoint::Vertex(v) => values[v as usize],
                SurfacePoint::Edge { edge, t } => {
                    let [a, b] = mesh.edge_vertices(edge as usize);
                    (1.0 - t) * values[a] + t * values[b]
                }
                SurfacePoint::Face { face, bary } => {
                    let f = mesh.faces()[face as usize];
                    (0..3).map(|k| bary[k] * values[f[k] as usize]).sum()
                }
            })
            .collect()
    }
}

/// A hash of a caller-triangulation session's mesh: the faces and the edge lengths of the
/// seam-cut mesh (the identity of the state when no intrinsic triangulation is used).
pub(super) fn mesh_fingerprint(mesh: &TriMesh, geometry: &Geometry) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    let mut word = |w: u64| {
        for byte in w.to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    word(mesh.n_vertices() as u64);
    for face in mesh.faces() {
        for &v in face {
            word(u64::from(v));
        }
    }
    for l in geometry.edge_lengths() {
        word(l.to_bits());
    }
    h
}

/// The state of a session: the level-set functions on the vertices of the triangulation the
/// flow runs on, the fingerprint of that triangulation, and the options that reproduce it.
///
/// Opaque on purpose: a caller stores it (`to_bytes`) and hands it back
/// ([`Session::restore`](super::Session::restore),
/// [`Session::new_from_state`](super::Session::new_from_state)); a session whose triangulation
/// differs refuses it with [`CuttingError::StaleState`]. The per-corner view of the state
/// ([`Session::state`](super::Session::state)) is for display and for starting; it does not
/// contain the values at inserted vertices.
#[derive(Debug, Clone, PartialEq)]
pub struct SessionState {
    pub(super) fingerprint: u64,
    pub(super) n_regions: usize,
    pub(super) region_flow: bool,
    pub(super) values: Vec<Vec<f64>>,
    pub(super) triangulation: Triangulation,
}

const MAGIC: &[u8; 4] = b"DDGS";
const VERSION: u8 = 1;

impl SessionState {
    /// The fingerprint of the triangulation the state was made for.
    pub fn fingerprint(&self) -> u64 {
        self.fingerprint
    }

    /// The number of regions.
    pub fn n_regions(&self) -> usize {
        self.n_regions
    }

    /// The number of values per function (vertices of the triangulation the flow runs on).
    pub fn n_values(&self) -> usize {
        self.values.first().map_or(0, Vec::len)
    }

    /// The triangulation option that reproduces the triangulation the state belongs to: build
    /// a session with it (the same mesh and seams) and the fingerprints agree.
    pub fn triangulation(&self) -> Triangulation {
        self.triangulation
    }

    /// The same state with other values (one function per entry of `values`, in the same
    /// units): for editing a stored state, and for tests.
    ///
    /// # Errors
    /// [`CuttingError::BadState`] if the number of functions or of values per function differs
    /// from this state's, or a value is not finite.
    pub fn with_values(&self, values: Vec<Vec<f64>>) -> Result<SessionState, CuttingError> {
        if values.len() != self.values.len() || values.iter().any(|f| f.len() != self.n_values()) {
            return Err(bad(
                "the number of functions or of values differs from the state's",
            ));
        }
        if values.iter().flatten().any(|x| !x.is_finite()) {
            return Err(bad("a value is not finite"));
        }
        Ok(SessionState {
            values,
            ..self.clone()
        })
    }

    /// The functions: for the two-region flow one level-set function `φ`, otherwise one per
    /// region (negative inside the region), each with one value per vertex of the triangulation
    /// the flow runs on, in the caller's length unit.
    pub fn values(&self) -> &[Vec<f64>] {
        &self.values
    }

    /// The state as bytes (little endian, versioned).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut w = Writer(Vec::new());
        w.0.extend_from_slice(MAGIC);
        w.0.push(VERSION);
        w.u64(self.fingerprint);
        w.u64(self.n_regions as u64);
        w.u8(u8::from(self.region_flow));
        w.u64(self.values.len() as u64);
        w.u64(self.n_values() as u64);
        for f in &self.values {
            for &x in f {
                w.f64(x);
            }
        }
        match &self.triangulation {
            Triangulation::Caller => w.u8(0),
            Triangulation::Intrinsic(o) => {
                w.u8(1);
                w.intrinsic(o);
            }
        }
        w.0
    }

    /// Reads a state written by [`to_bytes`](Self::to_bytes).
    ///
    /// # Errors
    /// [`CuttingError::BadState`] for bytes that are not a state (wrong magic or version,
    /// truncated, trailing bytes, an unknown tag).
    pub fn from_bytes(bytes: &[u8]) -> Result<SessionState, CuttingError> {
        let mut r = Reader { bytes, at: 0 };
        if r.take(4)? != MAGIC {
            return Err(bad(
                "the bytes do not start with the magic of a session state",
            ));
        }
        let version = r.u8()?;
        if version != VERSION {
            return Err(bad(&format!(
                "version {version} of the format is not supported"
            )));
        }
        let fingerprint = r.u64()?;
        let n_regions = r.usize()?;
        let region_flow = r.u8()? != 0;
        let n_functions = r.usize()?;
        let n_values = r.usize()?;
        let total = n_functions
            .checked_mul(n_values)
            .filter(|t| t.checked_mul(8).is_some_and(|b| b <= bytes.len()))
            .ok_or_else(|| bad("the number of values does not fit in the bytes"))?;
        let mut values = Vec::with_capacity(n_functions);
        for _ in 0..n_functions {
            let mut f = Vec::with_capacity(n_values);
            for _ in 0..n_values {
                f.push(r.f64()?);
            }
            values.push(f);
        }
        let _ = total;
        let triangulation = match r.u8()? {
            0 => Triangulation::Caller,
            1 => Triangulation::Intrinsic(r.intrinsic()?),
            t => return Err(bad(&format!("unknown triangulation tag {t}"))),
        };
        if r.at != bytes.len() {
            return Err(bad("trailing bytes after the state"));
        }
        Ok(SessionState {
            fingerprint,
            n_regions,
            region_flow,
            values,
            triangulation,
        })
    }
}

fn bad(reason: &str) -> CuttingError {
    CuttingError::BadState {
        reason: reason.to_string(),
    }
}

struct Writer(Vec<u8>);

impl Writer {
    fn u8(&mut self, x: u8) {
        self.0.push(x);
    }
    fn u64(&mut self, x: u64) {
        self.0.extend_from_slice(&x.to_le_bytes());
    }
    fn f64(&mut self, x: f64) {
        self.u64(x.to_bits());
    }
    fn opt_f64(&mut self, x: Option<f64>) {
        match x {
            Some(v) => {
                self.u8(1);
                self.f64(v);
            }
            None => self.u8(0),
        }
    }
    fn criterion(&mut self, c: &DelaunayCriterion) {
        match *c {
            DelaunayCriterion::AngleSum { tolerance } => {
                self.u8(0);
                self.f64(tolerance);
            }
            DelaunayCriterion::CotanWeight { tolerance } => {
                self.u8(1);
                self.f64(tolerance);
            }
        }
    }
    fn delaunay(&mut self, d: &DelaunayOptions) {
        self.criterion(&d.criterion);
        let (tag, n) = match d.max_flips {
            FlipLimit::PerEdge(n) => (0, n),
            FlipLimit::Flips(n) => (1, n),
            FlipLimit::Unlimited => (2, 0),
        };
        self.u8(tag);
        self.u64(n as u64);
    }
    fn refine(&mut self, r: &RefineOptions) {
        self.f64(r.min_angle);
        self.opt_f64(r.max_circumradius);
        let (tag, n) = match r.max_insertions {
            InsertionLimit::PerVertex(n) => (0, n),
            InsertionLimit::Vertices(n) => (1, n),
            InsertionLimit::Unlimited => (2, 0),
        };
        self.u8(tag);
        self.u64(n as u64);
        self.criterion(&r.criterion);
    }
    fn opt_refine(&mut self, r: &Option<RefineOptions>) {
        match r {
            Some(r) => {
                self.u8(1);
                self.refine(r);
            }
            None => self.u8(0),
        }
    }
    fn intrinsic(&mut self, o: &IntrinsicOptions) {
        self.delaunay(&o.delaunay);
        self.opt_refine(&o.refine);
        match &o.export.delaunay {
            Some(d) => {
                self.u8(1);
                self.delaunay(d);
            }
            None => self.u8(0),
        }
        self.opt_refine(&o.export.refine);
        self.u64(o.export.max_refine_rounds as u64);
        self.u64(o.export.max_splits as u64);
        self.f64(o.trace.tolerance);
        match o.trace.max_crossings {
            Some(n) => {
                self.u8(1);
                self.u64(n as u64);
            }
            None => self.u8(0),
        }
        self.f64(o.min_offset);
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Reader<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8], CuttingError> {
        let end = self
            .at
            .checked_add(n)
            .filter(|&e| e <= self.bytes.len())
            .ok_or_else(|| bad("the bytes end inside the state"))?;
        let out = &self.bytes[self.at..end];
        self.at = end;
        Ok(out)
    }
    fn u8(&mut self) -> Result<u8, CuttingError> {
        Ok(self.take(1)?[0])
    }
    fn u64(&mut self) -> Result<u64, CuttingError> {
        let b: [u8; 8] = self.take(8)?.try_into().expect("8 bytes");
        Ok(u64::from_le_bytes(b))
    }
    fn usize(&mut self) -> Result<usize, CuttingError> {
        usize::try_from(self.u64()?).map_err(|_| bad("a count does not fit in usize"))
    }
    fn f64(&mut self) -> Result<f64, CuttingError> {
        Ok(f64::from_bits(self.u64()?))
    }
    fn opt_f64(&mut self) -> Result<Option<f64>, CuttingError> {
        Ok(if self.u8()? == 1 {
            Some(self.f64()?)
        } else {
            None
        })
    }
    fn criterion(&mut self) -> Result<DelaunayCriterion, CuttingError> {
        let tag = self.u8()?;
        let tolerance = self.f64()?;
        match tag {
            0 => Ok(DelaunayCriterion::AngleSum { tolerance }),
            1 => Ok(DelaunayCriterion::CotanWeight { tolerance }),
            t => Err(bad(&format!("unknown Delaunay criterion tag {t}"))),
        }
    }
    fn delaunay(&mut self) -> Result<DelaunayOptions, CuttingError> {
        let criterion = self.criterion()?;
        let tag = self.u8()?;
        let n = self.usize()?;
        let max_flips = match tag {
            0 => FlipLimit::PerEdge(n),
            1 => FlipLimit::Flips(n),
            2 => FlipLimit::Unlimited,
            t => return Err(bad(&format!("unknown flip limit tag {t}"))),
        };
        Ok(DelaunayOptions {
            criterion,
            max_flips,
        })
    }
    fn refine(&mut self) -> Result<RefineOptions, CuttingError> {
        let min_angle = self.f64()?;
        let max_circumradius = self.opt_f64()?;
        let tag = self.u8()?;
        let n = self.usize()?;
        let max_insertions = match tag {
            0 => InsertionLimit::PerVertex(n),
            1 => InsertionLimit::Vertices(n),
            2 => InsertionLimit::Unlimited,
            t => return Err(bad(&format!("unknown insertion limit tag {t}"))),
        };
        let criterion = self.criterion()?;
        Ok(RefineOptions {
            min_angle,
            max_circumradius,
            max_insertions,
            criterion,
        })
    }
    fn opt_refine(&mut self) -> Result<Option<RefineOptions>, CuttingError> {
        Ok(if self.u8()? == 1 {
            Some(self.refine()?)
        } else {
            None
        })
    }
    fn intrinsic(&mut self) -> Result<IntrinsicOptions, CuttingError> {
        let delaunay = self.delaunay()?;
        let refine = self.opt_refine()?;
        let export_delaunay = if self.u8()? == 1 {
            Some(self.delaunay()?)
        } else {
            None
        };
        let export_refine = self.opt_refine()?;
        let max_refine_rounds = self.usize()?;
        let max_splits = self.usize()?;
        let tolerance = self.f64()?;
        let max_crossings = if self.u8()? == 1 {
            Some(self.usize()?)
        } else {
            None
        };
        let min_offset = self.f64()?;
        Ok(IntrinsicOptions {
            delaunay,
            refine,
            export: ExportOptions {
                delaunay: export_delaunay,
                refine: export_refine,
                max_refine_rounds,
                max_splits,
            },
            trace: EdgeTraceOptions {
                tolerance,
                max_crossings,
            },
            min_offset,
        })
    }
}
