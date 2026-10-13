//! The cut flow for the Python module: a session on a polygon mesh, in plain Rust.
//!
//! This is the binding of [`ddg_uv::cutting`], the polygon-mesh front end of the cut flow of
//! "Variational Surface Cutting" (Sharp, Crane, ACM TOG 2018; `ddg-variational-cuts`, implemented
//! clean-room from the text of the paper). Like the rest of the crate it contains no Python:
//! [`crate::python`] turns buffers into the types here and the results into arrays and
//! dictionaries.
//!
//! * [`Controls`] holds every control of the flow as a plain value with the library's default
//!   (energy, length weight or ripple wavelength, penalty, region weights, symmetry classes,
//!   curvature override, stage set and individual stages, step-size rule, stop rule) and
//!   turns them into the library's [`CutProblem`], [`FlowOptions`] and [`StopRule`].
//! * [`CutSession`] wraps [`Session`]: the state as one flat row-major array, `step`, `run`
//!   (reports as [`Dict`]s), the preview, the cut on polygons ([`DescriptionArrays`]), the
//!   insertion ([`InsertedArrays`]) and the end-to-end unwrap ([`UnwrapArrays`]).
//! * [`BindError`] sorts every failure into input, mesh, "not available yet" and solve errors;
//!   the message is the library's, unchanged.
//!
//! **Units.** As in the library, the weights, the penalty, the step sizes and the ripple
//! wavelength refer to the surface scaled to unit total area ([`CutSession::length_scale`] is
//! the factor from a length of the mesh to a length there); the state and the positions are
//! in the unit of the mesh. A pointwise energy given as a callable is not offered.

use std::time::Instant;

use ddg_core::TriMesh;
use ddg_uv::cutting::{
    ConePlacement, CutOnPolygons, CutUnwrapOptions, CuttingError, DEFAULT_CROSSING_OFFSET,
    Inserted, PathEnd, PreviewPolyline, Session, Start, VertexSource, insert,
};
use ddg_uv::{BoundaryValues, ConeSigns, ErrorKind, MeshInput, PolygonMesh};
use ddg_variational_cuts as vc;
use ddg_variational_cuts::{
    CurvatureSource, CutError, CutProblem, Energy, FlowOptions, FluxRecovery, RegionWeights,
    ScaleWeight, StepReport, StepSize, StopReason, StopRule,
};

use crate::contract::UnwrapInput;

// ----------------------------------------------------------------------------------------------
// Values and errors
// ----------------------------------------------------------------------------------------------

/// A value of a report. The Python side makes `Bool`, `Int`, `Float`, `Text` plain scalars,
/// `TextList` and `List` lists and `Dict` a dictionary.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// A flag.
    Bool(bool),
    /// An integer.
    Int(i64),
    /// A real number.
    Float(f64),
    /// A string.
    Text(String),
    /// A list of strings.
    TextList(Vec<String>),
    /// A list of values.
    List(Vec<Value>),
    /// A dictionary (insertion order is kept).
    Dict(Dict),
}

/// A dictionary of a report: keys in insertion order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Dict(pub Vec<(String, Value)>);

impl Dict {
    /// Adds an entry.
    pub fn set(&mut self, key: &str, value: impl Into<Value>) {
        self.0.push((key.to_owned(), value.into()));
    }

    /// The value of a key.
    pub fn get(&self, key: &str) -> Option<&Value> {
        self.0.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    /// The real number under a key.
    pub fn float(&self, key: &str) -> Option<f64> {
        match self.get(key) {
            Some(Value::Float(x)) => Some(*x),
            _ => None,
        }
    }
}

impl From<bool> for Value {
    fn from(x: bool) -> Value {
        Value::Bool(x)
    }
}
impl From<usize> for Value {
    fn from(x: usize) -> Value {
        Value::Int(x as i64)
    }
}
impl From<i64> for Value {
    fn from(x: i64) -> Value {
        Value::Int(x)
    }
}
impl From<f64> for Value {
    fn from(x: f64) -> Value {
        Value::Float(x)
    }
}
impl From<&str> for Value {
    fn from(x: &str) -> Value {
        Value::Text(x.to_owned())
    }
}
impl From<String> for Value {
    fn from(x: String) -> Value {
        Value::Text(x)
    }
}
impl From<Dict> for Value {
    fn from(x: Dict) -> Value {
        Value::Dict(x)
    }
}

/// Why a call failed. The message is the library's (it names vertices, corners and polygons in
/// the caller's numbering) or, for an option of this binding, names the keyword.
#[derive(Debug, Clone, PartialEq)]
pub enum BindError {
    /// An argument is malformed or outside its range (`InputError`).
    Input(String),
    /// The mesh cannot be used as it is (`MeshError`).
    Mesh(String),
    /// The request is valid but not implemented yet, for example more than two regions
    /// (`NotAvailableError`).
    NotAvailable(String),
    /// A numerical or internal failure (`SolveError`).
    Solve(String),
}

impl BindError {
    /// The message.
    pub fn message(&self) -> &str {
        match self {
            BindError::Input(m)
            | BindError::Mesh(m)
            | BindError::NotAvailable(m)
            | BindError::Solve(m) => m,
        }
    }
}

impl std::fmt::Display for BindError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.message())
    }
}

impl std::error::Error for BindError {}

fn input_error(message: impl Into<String>) -> BindError {
    BindError::Input(message.into())
}

impl From<CutError> for BindError {
    fn from(e: CutError) -> BindError {
        let message = e.to_string();
        match e {
            CutError::UnsupportedRegionCount { .. } | CutError::UnsupportedTriplePoint { .. } => {
                BindError::NotAvailable(message)
            }
            CutError::WrongLength { .. }
            | CutError::NonFiniteLevelSet { .. }
            | CutError::InvalidRegionWeight { .. }
            | CutError::InvalidLengthWeight { .. }
            | CutError::InvalidPenalty { .. }
            | CutError::InvalidScaleWeight { .. }
            | CutError::InvalidCurvature { .. }
            | CutError::InvalidFeatureSize { .. }
            | CutError::InvalidOptions { .. }
            | CutError::IndexOutOfRange { .. } => BindError::Input(message),
            CutError::Mesh(_) | CutError::InvalidTotalArea { .. } => BindError::Mesh(message),
            _ => BindError::Solve(message),
        }
    }
}

impl From<CuttingError> for BindError {
    fn from(e: CuttingError) -> BindError {
        let message = e.to_string();
        match e {
            CuttingError::Mesh(inner) => match inner.kind() {
                ErrorKind::Input => BindError::Input(message),
                ErrorKind::Mesh => BindError::Mesh(message),
                _ => BindError::Solve(message),
            },
            CuttingError::Flow(inner) => {
                // Keep the class of the flow's error, but the message with its prefix.
                match BindError::from(inner) {
                    BindError::Input(_) => BindError::Input(message),
                    BindError::Mesh(_) => BindError::Mesh(message),
                    BindError::NotAvailable(_) => BindError::NotAvailable(message),
                    BindError::Solve(_) => BindError::Solve(message),
                }
            }
            CuttingError::RegionsNotAvailable { .. } | CuttingError::Unsupported { .. } => {
                BindError::NotAvailable(message)
            }
            CuttingError::Internal { .. } => BindError::Solve(message),
            _ => BindError::Input(message),
        }
    }
}

// ----------------------------------------------------------------------------------------------
// Controls
// ----------------------------------------------------------------------------------------------

/// The energy of the distortion term (`pointwise` energies are not offered).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyKind {
    /// `∫ |∇u|²`.
    Dirichlet,
    /// `∫ a |∇u|²` with the curvature-based weight, or a per-vertex weight.
    ScaledDirichlet,
    /// `∫ u²`.
    Hencky,
}

/// The step-size rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepRule {
    /// [`StepSize::Adaptive`].
    Adaptive,
    /// [`StepSize::Capped`].
    Capped,
    /// [`StepSize::Fixed`].
    Fixed,
}

/// Every control of one `step` or `run`, as plain values. `None` means "the library's default".
///
/// The defaults of this struct are the library's: [`FlowOptions::new`] (the default stages),
/// [`StepSize::default`], [`StopRule::default`], [`EvalOptions::default`],
/// [`vc::PAPER_FEATURE_SMOOTHING`], `PenaltyRule::Gradient`. Three values are chosen by the
/// example program of the library and by nothing else: the energy ([`EnergyKind::Hencky`]), the
/// cap `r_ε = 1` of the curvature-based weight, and the ripple wavelength `0.25` that sets the
/// length weight when none is given.
#[derive(Debug, Clone)]
pub struct Controls {
    /// The distortion energy.
    pub energy: EnergyKind,
    /// A weight `a > 0` per vertex for [`EnergyKind::ScaledDirichlet`]; `None`: curvature-based.
    pub scale_weights: Option<Vec<f64>>,
    /// The cap `r_ε` of the curvature-based weight.
    pub feature_r_eps: f64,
    /// The smoothing constant `c` of the curvature-based weight.
    pub feature_smoothing: f64,
    /// The length weight `α_L`, explicit.
    pub length_weight: Option<f64>,
    /// The ripple wavelength: sets `α_L` from the roughening rate of the state.
    pub ripple_wavelength: Option<f64>,
    /// The penalty `η` per vertex.
    pub penalty: Option<Vec<f64>>,
    /// How the penalty enters the implicit step.
    pub penalty_rule: Option<vc::PenaltyRule>,
    /// The weight `w_k` of the distortion of every region.
    pub distortion_weights: Option<Vec<f64>>,
    /// The weight `β_k` of the area of every region.
    pub area_weights: Option<Vec<f64>>,
    /// Symmetry classes, one per vertex.
    pub symmetry: Option<Vec<u32>>,
    /// A uniform curvature density instead of the Gaussian curvature.
    pub curvature_uniform: Option<f64>,
    /// A curvature density per triangle of the triangulation.
    pub curvature_per_triangle: Option<Vec<f64>>,
    /// The stages of the paper ([`FlowOptions::paper`]) instead of the defaults.
    pub paper_stages: bool,
    /// How the normal derivative on the cut is recovered.
    pub flux_recovery: Option<FluxChoice>,
    /// The smoothing length of the smoothed recovery.
    pub flux_smoothing: Option<f64>,
    /// `ρ₀`, the aspect ratio above which a cut vertex is a sliver vertex.
    pub sliver_aspect_ratio: Option<f64>,
    /// Rounds of sliver smoothing of the speed.
    pub smoothing_passes: Option<usize>,
    /// Where the length and penalty terms are integrated.
    pub length_step: Option<vc::LengthStep>,
    /// How the rate on the cut is extended to the vertices.
    pub extension: Option<vc::Extension>,
    /// How the speed becomes a rate of the level-set function.
    pub speed_scaling: Option<vc::SpeedScaling>,
    /// What the redistancing does.
    pub redistance: Option<vc::Redistance>,
    /// The step-size rule.
    pub step_rule: StepRule,
    /// `time_scale` of [`StepSize::Adaptive`].
    pub time_scale: Option<f64>,
    /// `max_motion` of [`StepSize::Adaptive`] and [`StepSize::Capped`].
    pub max_motion: Option<f64>,
    /// `tau_max` of [`StepSize::Adaptive`] and [`StepSize::Capped`].
    pub tau_max: Option<f64>,
    /// `τ` of [`StepSize::Fixed`] (required for that rule).
    pub tau: Option<f64>,
    /// `δ`, the smallest distance of a crossing from a vertex, as a fraction of its edge.
    pub min_crossing_offset: Option<f64>,
    /// The stop rule: at most this many steps.
    pub max_steps: Option<usize>,
    /// The stop rule: cut motion tolerance (`0` never stops for it).
    pub motion_tolerance: Option<f64>,
    /// The stop rule: relative energy tolerance (`0` never stops for it).
    pub energy_tolerance: Option<f64>,
    /// The stop rule: window of the energy criterion.
    pub energy_window: Option<usize>,
}

/// The recovery of the normal derivative on the cut.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FluxChoice {
    /// [`FluxRecovery::Nodal`].
    Nodal,
    /// [`FluxRecovery::Smoothed`].
    Smoothed,
}

/// The ripple wavelength used when neither a length weight nor a wavelength is given (the
/// default of the library's example program, on the unit-area surface).
pub const DEFAULT_RIPPLE_WAVELENGTH: f64 = 0.25;

/// The cap `r_ε` of the curvature-based weight when none is given (the default of the library's
/// example program, on the unit-area surface).
pub const DEFAULT_FEATURE_R_EPS: f64 = 1.0;

impl Default for Controls {
    fn default() -> Controls {
        Controls {
            energy: EnergyKind::Hencky,
            scale_weights: None,
            feature_r_eps: DEFAULT_FEATURE_R_EPS,
            feature_smoothing: vc::PAPER_FEATURE_SMOOTHING,
            length_weight: None,
            ripple_wavelength: None,
            penalty: None,
            penalty_rule: None,
            distortion_weights: None,
            area_weights: None,
            symmetry: None,
            curvature_uniform: None,
            curvature_per_triangle: None,
            paper_stages: false,
            flux_recovery: None,
            flux_smoothing: None,
            sliver_aspect_ratio: None,
            smoothing_passes: None,
            length_step: None,
            extension: None,
            speed_scaling: None,
            redistance: None,
            step_rule: StepRule::Adaptive,
            time_scale: None,
            max_motion: None,
            tau_max: None,
            tau: None,
            min_crossing_offset: None,
            max_steps: None,
            motion_tolerance: None,
            energy_tolerance: None,
            energy_window: None,
        }
    }
}

impl Controls {
    /// The step size of the rule, with the library's defaults for what is not given.
    fn step_size(&self) -> Result<StepSize, BindError> {
        let StepSize::Adaptive {
            time_scale: d_time,
            max_motion: d_motion,
            tau_max: d_tau_max,
        } = StepSize::default()
        else {
            unreachable!("the default step size is adaptive");
        };
        let unused = |name: &str, given: bool, rule: &str| {
            if given {
                Err(input_error(format!(
                    "{name} does not apply to step_rule={rule:?}"
                )))
            } else {
                Ok(())
            }
        };
        Ok(match self.step_rule {
            StepRule::Adaptive => {
                unused("tau", self.tau.is_some(), "adaptive")?;
                StepSize::Adaptive {
                    time_scale: self.time_scale.unwrap_or(d_time),
                    max_motion: self.max_motion.unwrap_or(d_motion),
                    tau_max: self.tau_max.unwrap_or(d_tau_max),
                }
            }
            StepRule::Capped => {
                unused("tau", self.tau.is_some(), "capped")?;
                unused("time_scale", self.time_scale.is_some(), "capped")?;
                StepSize::Capped {
                    tau_max: self.tau_max.unwrap_or(d_tau_max),
                    max_motion: self.max_motion.unwrap_or(d_motion),
                }
            }
            StepRule::Fixed => {
                unused("time_scale", self.time_scale.is_some(), "fixed")?;
                unused("max_motion", self.max_motion.is_some(), "fixed")?;
                unused("tau_max", self.tau_max.is_some(), "fixed")?;
                StepSize::Fixed(
                    self.tau
                        .ok_or_else(|| input_error("tau is required for step_rule=\"fixed\""))?,
                )
            }
        })
    }

    /// The options of the flow.
    ///
    /// # Errors
    /// [`BindError::Input`] for contradictory options.
    pub fn flow_options(&self) -> Result<FlowOptions<'static>, BindError> {
        let step = self.step_size()?;
        let mut o = if self.paper_stages {
            FlowOptions::paper(step)
        } else {
            FlowOptions::new(step)
        };
        if let Some(delta) = self.min_crossing_offset {
            o.eval.min_crossing_offset = delta;
        }
        match (self.flux_recovery, self.flux_smoothing) {
            (Some(FluxChoice::Nodal), Some(_)) => {
                return Err(input_error(
                    "flux_smoothing does not apply to flux_recovery=\"nodal\"",
                ));
            }
            (Some(FluxChoice::Nodal), None) => o.eval.flux_recovery = FluxRecovery::Nodal,
            (Some(FluxChoice::Smoothed), length) | (None, length @ Some(_)) => {
                let default = match FluxRecovery::default() {
                    FluxRecovery::Smoothed { length } => length,
                    FluxRecovery::Nodal => 0.0,
                };
                o.eval.flux_recovery = FluxRecovery::Smoothed {
                    length: length.unwrap_or(default),
                };
            }
            (None, None) => {}
        }
        if let Some(r) = self.sliver_aspect_ratio {
            o.eval.sliver_aspect_ratio = r;
        }
        if let Some(n) = self.smoothing_passes {
            o.eval.smoothing_passes = n;
        }
        if let Some(x) = self.length_step {
            o.length_step = x;
        }
        if let Some(x) = self.extension {
            o.extension = x;
        }
        if let Some(x) = self.speed_scaling {
            o.speed_scaling = x;
        }
        if let Some(x) = self.redistance {
            o.redistance = x;
        }
        if let Some(x) = self.penalty_rule {
            o.penalty_rule = x;
        }
        o.symmetry = self.symmetry.clone();
        Ok(o)
    }

    /// The stop rule of a run.
    pub fn stop_rule(&self) -> StopRule {
        let d = StopRule::default();
        StopRule {
            max_steps: self.max_steps.unwrap_or(d.max_steps),
            motion_tolerance: self.motion_tolerance.unwrap_or(d.motion_tolerance),
            energy_tolerance: self.energy_tolerance.unwrap_or(d.energy_tolerance),
            energy_window: self.energy_window.unwrap_or(d.energy_window),
        }
    }

    /// The problem for a session with `n_regions` regions and the length weight `alpha`.
    ///
    /// # Errors
    /// [`BindError::Input`] for contradictory options or arrays of the wrong length.
    pub fn problem(&self, n_regions: usize, alpha: f64) -> Result<CutProblem, BindError> {
        let energy = match (self.energy, &self.scale_weights) {
            (EnergyKind::Dirichlet, None) => Energy::Dirichlet,
            (EnergyKind::Hencky, None) => Energy::Hencky,
            (EnergyKind::ScaledDirichlet, Some(a)) => {
                Energy::ScaledDirichlet(ScaleWeight::PerVertex(a.clone()))
            }
            (EnergyKind::ScaledDirichlet, None) => {
                Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
                    r_eps: self.feature_r_eps,
                    smoothing: self.feature_smoothing,
                })
            }
            (_, Some(_)) => {
                return Err(input_error(
                    "scale_weights is only used by energy=\"scaled_dirichlet\"",
                ));
            }
        };
        let curvature = match (&self.curvature_uniform, &self.curvature_per_triangle) {
            (Some(_), Some(_)) => {
                return Err(input_error(
                    "curvature_uniform and curvature_per_triangle exclude each other",
                ));
            }
            (Some(k), None) => CurvatureSource::Uniform(*k),
            (None, Some(k)) => CurvatureSource::PerFace(k.clone()),
            (None, None) => CurvatureSource::Gaussian,
        };
        let weights = |name: &str, values: &Option<Vec<f64>>, default: f64| match values {
            None => Ok(vec![default; n_regions]),
            Some(v) if v.len() == n_regions => Ok(v.clone()),
            Some(v) => Err(input_error(format!(
                "{name} has {} entries but needs one per region ({n_regions})",
                v.len()
            ))),
        };
        let distortion = weights("distortion_weights", &self.distortion_weights, 1.0)?;
        let area = weights("area_weights", &self.area_weights, 0.0)?;
        Ok(CutProblem {
            energy,
            curvature,
            regions: distortion
                .into_iter()
                .zip(area)
                .map(|(distortion, area)| RegionWeights { distortion, area })
                .collect(),
            length_weight: alpha,
            penalty: self.penalty.clone(),
        })
    }

    /// Whether the length weight comes from a ripple wavelength (no explicit weight).
    fn wavelength(&self) -> Result<Option<f64>, BindError> {
        match (self.length_weight, self.ripple_wavelength) {
            (Some(_), Some(_)) => Err(input_error(
                "length_weight and ripple_wavelength exclude each other",
            )),
            (Some(_), None) => Ok(None),
            (None, Some(w)) => Ok(Some(w)),
            (None, None) => Ok(Some(DEFAULT_RIPPLE_WAVELENGTH)),
        }
    }
}

// ----------------------------------------------------------------------------------------------
// The session
// ----------------------------------------------------------------------------------------------

/// How a session starts or restarts.
#[derive(Debug, Clone)]
pub enum StartSpec {
    /// The geodesic disk of a radius around a vertex: region 0.
    Circle {
        /// The vertex.
        vertex: u32,
        /// The radius, in the unit of the mesh.
        radius: f64,
    },
    /// One label per polygon.
    Labels(Vec<u32>),
    /// Region functions per corner, row-major `(n_regions, n_corners)`.
    Corners(Vec<f64>),
}

/// The options of the end-to-end unwrap.
#[derive(Debug, Clone)]
pub struct UnwrapSettings {
    /// `λ` of the cone solver.
    pub lambda: f64,
    /// Only cone tips.
    pub positive_cones_only: bool,
    /// Let the solver choose the scale along boundaries.
    pub optimize_boundary: bool,
    /// Whether the unwrap may place cones besides the cut.
    pub cones: bool,
}

impl Default for UnwrapSettings {
    fn default() -> Self {
        let d = ddg_uv::UnwrapOptions::default();
        UnwrapSettings {
            lambda: d.lambda,
            positive_cones_only: true,
            optimize_boundary: false,
            cones: false,
        }
    }
}

/// A cut session on a polygon mesh: [`Session`] with the sizes the binding checks against.
#[derive(Debug)]
pub struct CutSession {
    inner: Session,
    n_regions: usize,
    n_vertices: usize,
    n_corners: usize,
    n_polygons: usize,
    n_triangles: usize,
}

fn start_of(spec: &StartSpec, n_regions: usize, n_corners: usize) -> Result<Start, BindError> {
    Ok(match spec {
        StartSpec::Circle { vertex, radius } => Start::GeodesicCircle {
            vertex: *vertex,
            radius: *radius,
        },
        StartSpec::Labels(labels) => Start::PolygonLabels(labels.clone()),
        StartSpec::Corners(flat) => {
            if flat.len() != n_regions * n_corners {
                return Err(input_error(format!(
                    "corner_values has {} entries but needs n_regions * n_corners = {} * {} = {}",
                    flat.len(),
                    n_regions,
                    n_corners,
                    n_regions * n_corners
                )));
            }
            Start::Corners(flat.chunks(n_corners.max(1)).map(<[f64]>::to_vec).collect())
        }
    })
}

fn prescribed_cones(input: &UnwrapInput) -> Vec<ddg_uv::Cone> {
    input
        .prescribed_cones()
        .iter()
        .map(|c| ddg_uv::Cone {
            vertex: c.vertex,
            angle: c.angle,
        })
        .collect()
}

fn mesh_input<'a>(input: &'a UnwrapInput, prescribed: &'a [ddg_uv::Cone]) -> MeshInput<'a> {
    MeshInput {
        positions: input.positions(),
        corner_vertices: input.corner_vertices(),
        polygon_offsets: input.polygon_offsets(),
        triangles: input.triangles(),
        seam_edges: input.seam_edges(),
        pinned_vertices: input.pinned_vertices(),
        prescribed_cones: prescribed,
        forbidden_vertices: input.forbidden_vertices(),
    }
}

impl CutSession {
    /// Builds a session on the validated mesh `input` with `n_regions` regions.
    ///
    /// # Errors
    /// As [`Session::new`], sorted by [`BindError`].
    pub fn new(
        input: &UnwrapInput,
        n_regions: usize,
        start: &StartSpec,
    ) -> Result<CutSession, BindError> {
        let prescribed = prescribed_cones(input);
        let mesh = mesh_input(input, &prescribed);
        let start = start_of(start, n_regions, input.n_corners())?;
        let inner = Session::new(&mesh, n_regions, &start)?;
        Ok(CutSession {
            inner,
            n_regions,
            n_vertices: input.n_vertices(),
            n_corners: input.n_corners(),
            n_polygons: input.n_polygons(),
            n_triangles: input.n_triangles(),
        })
    }

    /// The number of regions.
    pub fn n_regions(&self) -> usize {
        self.n_regions
    }
    /// The number of vertices of the polygon mesh.
    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }
    /// The number of polygon corners.
    pub fn n_corners(&self) -> usize {
        self.n_corners
    }
    /// The number of polygons.
    pub fn n_polygons(&self) -> usize {
        self.n_polygons
    }
    /// The number of triangles of the supplied triangulation.
    pub fn n_triangles(&self) -> usize {
        self.n_triangles
    }
    /// The factor from a length of the mesh to a length on the unit-area surface.
    pub fn length_scale(&self) -> f64 {
        self.inner.flow().scale().length_scale
    }
    /// The total area of the mesh (of the triangulation), in the unit of the mesh squared.
    pub fn total_area(&self) -> f64 {
        self.inner.flow().scale().total_area
    }
    /// The flow time since the last start or `set_state`.
    pub fn time(&self) -> f64 {
        self.inner.flow().time()
    }
    /// The number of steps since the last start or `set_state`.
    pub fn n_steps(&self) -> usize {
        self.inner.flow().n_steps()
    }
    /// Warnings of the validation of the mesh.
    pub fn warnings(&self) -> &[String] {
        self.inner.warnings()
    }
    /// The seams that were cut open before the flow, as vertex pairs `[min, max]`.
    pub fn seams_cut(&self) -> &[[u32; 2]] {
        self.inner.seams_cut()
    }

    /// Replaces the state by a new start.
    ///
    /// # Errors
    /// As [`Session::restart`].
    pub fn restart(&mut self, start: &StartSpec) -> Result<(), BindError> {
        let start = start_of(start, self.n_regions, self.n_corners)?;
        Ok(self.inner.restart(&start)?)
    }

    /// The state, row-major `(n_regions, n_corners)`.
    pub fn state(&self) -> Vec<f64> {
        self.inner.state().concat()
    }

    /// Replaces the state; `flat` is row-major `(n_regions, n_corners)`.
    ///
    /// # Errors
    /// [`BindError::Input`] for a wrong length, as [`Session::set_state`] otherwise.
    pub fn set_state(&mut self, flat: &[f64]) -> Result<(), BindError> {
        if flat.len() != self.n_regions * self.n_corners {
            return Err(input_error(format!(
                "state has {} entries but needs n_regions * n_corners = {} * {} = {}",
                flat.len(),
                self.n_regions,
                self.n_corners,
                self.n_regions * self.n_corners
            )));
        }
        let rows: Vec<&[f64]> = flat.chunks(self.n_corners.max(1)).collect();
        Ok(self.inner.set_state(&rows)?)
    }

    /// The region of every polygon corner.
    pub fn corner_regions(&self) -> Vec<u32> {
        self.inner.corner_regions()
    }

    /// The length weight of this call: the explicit one, or the one the ripple wavelength gives
    /// at the current state. Returns `(α_L, ω)`; `ω` is `None` for an explicit weight.
    fn resolve_length_weight(
        &mut self,
        c: &Controls,
        options: &FlowOptions<'_>,
    ) -> Result<(f64, Option<f64>), BindError> {
        match c.wavelength()? {
            None => Ok((c.length_weight.unwrap_or(0.0), None)),
            Some(wavelength) => {
                if !(wavelength.is_finite() && wavelength > 0.0) {
                    return Err(input_error(format!(
                        "ripple_wavelength must be positive and finite, got {wavelength}"
                    )));
                }
                let problem = c.problem(self.n_regions, 0.0)?;
                let omega = self.inner.roughening_rate(&problem, &options.eval)?;
                Ok((
                    vc::length_weight_for_ripples(omega, wavelength),
                    Some(omega),
                ))
            }
        }
    }

    /// The roughening rate `ω` of the current state (the scale of the length weight).
    ///
    /// # Errors
    /// As [`Session::roughening_rate`].
    pub fn roughening_rate(&mut self, c: &Controls) -> Result<f64, BindError> {
        let options = c.flow_options()?;
        let problem = c.problem(self.n_regions, 0.0)?;
        Ok(self.inner.roughening_rate(&problem, &options.eval)?)
    }

    /// One step of the flow; the report of the step plus `length_weight`.
    ///
    /// # Errors
    /// As [`Session::step`].
    pub fn step(&mut self, c: &Controls) -> Result<Dict, BindError> {
        let started = Instant::now();
        let options = c.flow_options()?;
        let (alpha, omega) = self.resolve_length_weight(c, &options)?;
        let problem = c.problem(self.n_regions, alpha)?;
        let report = self.inner.step(&problem, &options)?;
        let mut d = step_dict(&report, self.length_scale());
        d.set("length_weight", alpha);
        if let Some(omega) = omega {
            d.set("length_weight_roughening_rate", omega);
        }
        d.set("seconds_total", started.elapsed().as_secs_f64());
        Ok(d)
    }

    /// Steps until the stop rule says so. The length weight, if it comes from a wavelength, is
    /// set once at the start of the run.
    ///
    /// # Errors
    /// As [`Session::run`].
    pub fn run(&mut self, c: &Controls) -> Result<Dict, BindError> {
        let started = Instant::now();
        let options = c.flow_options()?;
        let (alpha, omega) = self.resolve_length_weight(c, &options)?;
        let problem = c.problem(self.n_regions, alpha)?;
        let report = self.inner.run(&problem, &options, &c.stop_rule())?;
        let scale = self.length_scale();
        let mut d = Dict::default();
        d.set("steps_taken", report.steps.len());
        d.set("stop_reason", stop_name(report.stop));
        let before = report.steps.first().map(|s| &s.energy);
        let after = &report.final_energy;
        let before = before.unwrap_or(after);
        d.set("energy_before", before.total);
        d.set("energy_after", after.total);
        d.set("cut_length_before", before.length / scale);
        d.set("cut_length_after", after.length / scale);
        d.set("cut_length_before_unit_area", before.length);
        d.set("cut_length_after_unit_area", after.length);
        for (k, x) in after.distortion.iter().enumerate() {
            d.set(&format!("distortion_{k}"), *x);
        }
        for (k, x) in after.area.iter().enumerate() {
            d.set(&format!("area_{k}"), *x);
        }
        d.set("penalty_after", after.penalty);
        d.set("length_weight", alpha);
        if let Some(omega) = omega {
            d.set("length_weight_roughening_rate", omega);
        }
        d.set("time", self.time());
        d.set("n_steps_total", self.n_steps());
        let history: Vec<Value> = report
            .steps
            .iter()
            .map(|s| Value::Dict(step_dict(s, scale)))
            .collect();
        d.set(
            "seconds_steps",
            report
                .steps
                .iter()
                .map(|s| {
                    s.timings.evaluate
                        + s.timings.extension
                        + s.timings.length
                        + s.timings.redistance
                })
                .sum::<f64>(),
        );
        d.0.push(("steps".to_owned(), Value::List(history)));
        d.set("seconds_total", started.elapsed().as_secs_f64());
        Ok(d)
    }

    /// The cut as polylines for display.
    ///
    /// # Errors
    /// As [`Session::preview`].
    pub fn preview(&self, delta: Option<f64>) -> Result<PreviewArrays, BindError> {
        let lines = self
            .inner
            .preview(delta.unwrap_or(DEFAULT_CROSSING_OFFSET))?;
        Ok(PreviewArrays::from_polylines(&lines))
    }

    /// The cut on the polygons.
    ///
    /// # Errors
    /// As [`Session::describe_cut`].
    pub fn describe_cut(&self, delta: Option<f64>) -> Result<CutOnPolygons, BindError> {
        Ok(self
            .inner
            .describe_cut(delta.unwrap_or(DEFAULT_CROSSING_OFFSET))?)
    }

    /// Inserts a cut description (or the current cut) into the session's polygon mesh.
    ///
    /// # Errors
    /// As [`insert`].
    pub fn insert(
        &self,
        cut: Option<&CutOnPolygons>,
        delta: Option<f64>,
    ) -> Result<Inserted, BindError> {
        let own;
        let cut = match cut {
            Some(c) => c,
            None => {
                own = self.describe_cut(delta)?;
                &own
            }
        };
        Ok(insert(self.inner.mesh(), cut)?)
    }

    /// Describes the current cut, inserts it and unwraps the result.
    ///
    /// # Errors
    /// As [`Session::unwrap`].
    pub fn unwrap(
        &self,
        delta: Option<f64>,
        settings: &UnwrapSettings,
    ) -> Result<(CutOnPolygons, ddg_uv::cutting::CutUnwrap), BindError> {
        let options = CutUnwrapOptions {
            unwrap: ddg_uv::UnwrapOptions {
                lambda: settings.lambda,
                cone_signs: if settings.positive_cones_only {
                    ConeSigns::PositiveWhereFeasible
                } else {
                    ConeSigns::Any
                },
                boundary: if settings.optimize_boundary {
                    BoundaryValues::Optimized
                } else {
                    BoundaryValues::Zero
                },
                // `ddg_tools.pack_charts` packs.
                packing: None,
                ..ddg_uv::UnwrapOptions::default()
            },
            cones: if settings.cones {
                ConePlacement::Optimal
            } else {
                ConePlacement::Off
            },
        };
        Ok(self
            .inner
            .unwrap(delta.unwrap_or(DEFAULT_CROSSING_OFFSET), &options)?)
    }
}

fn stop_name(r: StopReason) -> &'static str {
    match r {
        StopReason::MaxSteps => "max_steps",
        StopReason::Converged => "converged",
        StopReason::EnergyStalled => "energy_stalled",
        StopReason::CutVanished => "cut_vanished",
    }
}

/// The report of one step as a dictionary. Lengths are given in the unit of the mesh
/// (`cut_length`) and on the unit-area surface (`cut_length_unit_area`); everything else is as
/// the library reports it.
fn step_dict(r: &StepReport, length_scale: f64) -> Dict {
    let mut d = Dict::default();
    d.set("time", r.time);
    d.set("tau", r.tau);
    d.set("energy_total", r.energy.total);
    for (k, x) in r.energy.distortion.iter().enumerate() {
        d.set(&format!("distortion_{k}"), *x);
    }
    for (k, x) in r.energy.area.iter().enumerate() {
        d.set(&format!("area_{k}"), *x);
    }
    d.set("cut_length", r.energy.length / length_scale);
    d.set("cut_length_unit_area", r.energy.length);
    d.set("penalty", r.energy.penalty);
    d.set("max_speed", r.max_speed);
    d.set("max_length_speed", r.max_length_speed);
    d.set("roughening_rate", r.roughening_rate);
    d.set("cut_motion", r.cut_motion);
    if let Some(g) = r.gradient_norm {
        d.set("gradient_norm_min", g.min);
        d.set("gradient_norm_mean", g.mean);
        d.set("gradient_norm_max", g.max);
    }
    d.set("n_cut_points", r.n_cut_points);
    d.set("n_clamped_crossings", r.n_clamped_crossings);
    d.set("n_uncut_closed_components", r.n_uncut_closed_components);
    d.set("n_unreached", r.n_unreached);
    d.set("extension_residual", r.extension_residual);
    if let Some(x) = r.length_residual {
        d.set("length_residual", x);
    }
    d.set("length_factorized", r.length_factorized);
    d.set("cut_vanished", r.cut_vanished);
    d.set("seconds_evaluate", r.timings.evaluate);
    d.set("seconds_extension", r.timings.extension);
    d.set("seconds_length", r.timings.length);
    d.set("seconds_redistance", r.timings.redistance);
    d
}

// ----------------------------------------------------------------------------------------------
// Arrays that leave the session
// ----------------------------------------------------------------------------------------------

/// The preview polylines as flat arrays.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PreviewArrays {
    /// 3 per point, the polylines one after the other.
    pub positions: Vec<f64>,
    /// `n_polylines + 1` offsets into the points (not into `positions`).
    pub starts: Vec<u32>,
    /// Whether the polyline is closed.
    pub closed: Vec<bool>,
}

impl PreviewArrays {
    fn from_polylines(lines: &[PreviewPolyline]) -> PreviewArrays {
        let mut out = PreviewArrays {
            starts: vec![0],
            ..PreviewArrays::default()
        };
        for line in lines {
            out.positions.extend(line.positions.iter().flatten());
            out.starts
                .push(out.starts.last().unwrap() + line.positions.len() as u32);
            out.closed.push(line.closed);
        }
        out
    }
}

/// The cut on polygons as flat arrays (see the README for the meaning of every array).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DescriptionArrays {
    /// 2 per edge point: the vertex pair `[min, max]` of the polygon edge.
    pub edge_vertices: Vec<u32>,
    /// The parameter from the first vertex.
    pub edge_t: Vec<f64>,
    /// 3 per edge point.
    pub edge_positions: Vec<f64>,
    /// The polygon whose side of a seam the point is on, or `-1`.
    pub edge_side: Vec<i64>,
    /// The polygon of every path.
    pub path_polygon: Vec<u32>,
    /// The start of every path: an edge point `>= 0`, or a junction `-(j + 1)`.
    pub path_start: Vec<i64>,
    /// The end of every path, encoded like the start.
    pub path_end: Vec<i64>,
    /// 2 per path: the regions on its left and right.
    pub path_regions: Vec<u32>,
    /// `n_paths + 1` offsets into the interior points.
    pub path_interior_starts: Vec<u32>,
    /// 2 per interior point: the polygon corners of the diagonal it lies on.
    pub interior_corners: Vec<u32>,
    /// The parameter from the first corner.
    pub interior_t: Vec<f64>,
    /// 3 per interior point.
    pub interior_positions: Vec<f64>,
    /// 3 per junction.
    pub junction_positions: Vec<f64>,
    /// The number of crossings the flow moved away from a vertex.
    pub n_clamped_crossings: usize,
}

fn end_code(e: PathEnd) -> i64 {
    match e {
        PathEnd::EdgePoint(i) => i64::from(i),
        PathEnd::Junction(j) => -(i64::from(j) + 1),
    }
}

impl DescriptionArrays {
    /// The arrays of a description.
    pub fn from_cut(cut: &CutOnPolygons) -> DescriptionArrays {
        let mut a = DescriptionArrays {
            path_interior_starts: vec![0],
            n_clamped_crossings: cut.n_clamped_crossings,
            ..DescriptionArrays::default()
        };
        for p in &cut.edge_points {
            a.edge_vertices.extend(p.edge);
            a.edge_t.push(p.t);
            a.edge_positions.extend(p.position);
            a.edge_side.push(p.side.map_or(-1, i64::from));
        }
        for p in &cut.paths {
            a.path_polygon.push(p.polygon);
            a.path_start.push(end_code(p.start));
            a.path_end.push(end_code(p.end));
            a.path_regions.extend(p.regions);
            for q in &p.interior {
                a.interior_corners.extend(q.corners);
                a.interior_t.push(q.t);
                a.interior_positions.extend(q.position);
            }
            a.path_interior_starts.push(a.interior_t.len() as u32);
        }
        for j in &cut.junctions {
            a.junction_positions.extend(j.position);
        }
        a
    }
}

/// The inserted mesh as flat arrays.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct InsertedArrays {
    /// 3 per vertex.
    pub positions: Vec<f64>,
    /// The vertex of every corner.
    pub corner_vertices: Vec<u32>,
    /// The first corner of every polygon.
    pub polygon_starts: Vec<u32>,
    /// 3 corner indices per triangle.
    pub triangles: Vec<u32>,
    /// 2 per seam edge of the cut, `[min, max]`.
    pub new_seams: Vec<u32>,
    /// The old polygon of every new polygon.
    pub polygon_parent: Vec<u32>,
    /// 3 per new vertex: kind (0 old vertex, 1 edge point, 2 interior point), then the index
    /// (old vertex, edge point) or the path and the position in it.
    pub vertex_source: Vec<u32>,
    /// The vertex of every edge point.
    pub edge_point_vertex: Vec<u32>,
    /// The number of vertices of the old mesh.
    pub n_old_vertices: usize,
}

impl InsertedArrays {
    /// The arrays of an insertion.
    pub fn from_inserted(i: &Inserted) -> InsertedArrays {
        let m: &PolygonMesh = &i.mesh;
        InsertedArrays {
            positions: m.positions.iter().flatten().copied().collect(),
            corner_vertices: m.corner_vertices.clone(),
            polygon_starts: m.polygon_offsets[..m.polygon_offsets.len() - 1].to_vec(),
            triangles: m.triangles.iter().flatten().copied().collect(),
            new_seams: i.new_seams.iter().flatten().copied().collect(),
            polygon_parent: i.polygon_parent.clone(),
            vertex_source: i
                .vertex_source
                .iter()
                .flat_map(|s| match *s {
                    VertexSource::Old(v) => [0, v, 0],
                    VertexSource::EdgePoint(e) => [1, e, 0],
                    VertexSource::Interior { path, index } => [2, path, index],
                })
                .collect(),
            edge_point_vertex: i.edge_point_vertex.clone(),
            n_old_vertices: i.n_old_vertices,
        }
    }
}

/// The result of the end-to-end unwrap as flat arrays and a report.
#[derive(Debug, Clone, PartialEq)]
pub struct UnwrapArrays {
    /// The insertion.
    pub inserted: InsertedArrays,
    /// 2 per seam of the new mesh (old seams and the cut).
    pub seams: Vec<u32>,
    /// 2 per corner of the new mesh.
    pub uvs: Vec<[f64; 2]>,
    /// The chart of every polygon of the new mesh.
    pub polygon_charts: Vec<u32>,
    /// The vertices with a cone.
    pub cone_vertices: Vec<u32>,
    /// The cone angles.
    pub cone_angles: Vec<f64>,
    /// The edges the unwrap cut along.
    pub cut_edges: Vec<u32>,
    /// The report of the unwrap.
    pub report: Dict,
    /// The warnings of the unwrap.
    pub warnings: Vec<String>,
}

impl UnwrapArrays {
    /// The arrays of an unwrap.
    pub fn from_unwrap(u: &ddg_uv::cutting::CutUnwrap) -> UnwrapArrays {
        let mut report = Dict::default();
        for (key, value) in u.map.report.summary() {
            match value {
                ddg_uv::Value::Bool(v) => report.set(key, v),
                ddg_uv::Value::Int(v) => report.set(key, v),
                ddg_uv::Value::Float(v) => report.set(key, v),
                ddg_uv::Value::Text(v) => report.set(key, v),
            }
        }
        UnwrapArrays {
            inserted: InsertedArrays::from_inserted(&u.inserted),
            seams: u.seams.iter().flatten().copied().collect(),
            uvs: u.map.uvs.clone(),
            polygon_charts: u.map.polygon_charts.clone(),
            cone_vertices: u.map.cones.iter().map(|c| c.vertex).collect(),
            cone_angles: u.map.cones.iter().map(|c| c.angle).collect(),
            cut_edges: u.map.cut_edges.iter().flatten().copied().collect(),
            report,
            warnings: u.map.report.warnings.clone(),
        }
    }
}

/// The visibility penalty of the paper (Sec. 7.1), `η_i = ν_i · max(n_i · e_i, 0)`, for the
/// triangulation of a validated polygon mesh. `visibility` is `ν` per vertex (1: the camera
/// sees the vertex), or `None` for 1 everywhere.
///
/// # Errors
/// [`BindError::Input`] for a visibility of the wrong length or a camera that is not finite,
/// [`BindError::Mesh`] if the triangulation is not a manifold triangle mesh.
pub fn visibility_penalty(
    input: &UnwrapInput,
    camera: [f64; 3],
    visibility: Option<&[f64]>,
) -> Result<Vec<f64>, BindError> {
    let faces: Vec<[u32; 3]> = input
        .triangles()
        .iter()
        .map(|t| t.map(|c| input.corner_vertices()[c as usize]))
        .collect();
    let mesh = TriMesh::new(input.n_vertices(), faces)
        .map_err(|e| BindError::Mesh(format!("the triangulation is not a manifold mesh: {e}")))?;
    let ones;
    let visibility = match visibility {
        Some(v) => v,
        None => {
            ones = vec![1.0; input.n_vertices()];
            &ones
        }
    };
    Ok(vc::init::visibility_penalty(
        &mesh,
        input.positions(),
        camera,
        visibility,
    )?)
}

/// `α_L = ω ℓ / (2π)`: the length weight at which ripples longer than `wavelength` grow.
pub fn length_weight_for_ripples(roughening_rate: f64, wavelength: f64) -> f64 {
    vc::length_weight_for_ripples(roughening_rate, wavelength)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::MeshBuffers;
    use ddg_uv::primitives::uv_sphere;

    fn sphere() -> UnwrapInput {
        let mesh = uv_sphere(24, 12, 1.0).expect("valid parameters");
        let buffers = MeshBuffers {
            positions: mesh.positions.iter().flatten().copied().collect(),
            corner_vertices: mesh.corner_vertices.clone(),
            polygon_starts: mesh.polygon_offsets[..mesh.polygon_offsets.len() - 1].to_vec(),
            triangles: mesh.triangles.iter().flatten().copied().collect(),
            ..MeshBuffers::default()
        };
        UnwrapInput::from_buffers(buffers).expect("valid mesh")
    }

    fn circle() -> StartSpec {
        StartSpec::Circle {
            vertex: 5,
            radius: 0.9,
        }
    }

    fn flow_controls() -> Controls {
        Controls {
            length_weight: Some(0.05),
            ..Controls::default()
        }
    }

    #[test]
    fn step_changes_the_state_and_reports_energy() {
        let input = sphere();
        let mut s = CutSession::new(&input, 2, &circle()).unwrap();
        assert_eq!(s.state().len(), 2 * s.n_corners());
        let before = s.state();
        let report = s.step(&flow_controls()).unwrap();
        assert!(report.float("energy_total").unwrap().is_finite());
        assert_ne!(before, s.state());
        assert_eq!(s.n_steps(), 1);
    }

    #[test]
    fn k_steps_twice_equal_2k_steps_once() {
        let input = sphere();
        let controls = Controls {
            max_steps: Some(4),
            motion_tolerance: Some(0.0),
            energy_tolerance: Some(0.0),
            ..flow_controls()
        };
        let mut once = CutSession::new(&input, 2, &circle()).unwrap();
        once.run(&Controls {
            max_steps: Some(8),
            ..controls.clone()
        })
        .unwrap();
        let mut twice = CutSession::new(&input, 2, &circle()).unwrap();
        twice.run(&controls).unwrap();
        let mid = twice.state();
        twice.set_state(&mid).unwrap();
        twice.run(&controls).unwrap();
        let (a, b) = (once.state(), twice.state());
        let worst = a
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        assert!(worst < 1e-12, "{worst}");
    }

    #[test]
    fn run_reports_energy_before_and_after() {
        let input = sphere();
        let mut s = CutSession::new(&input, 2, &circle()).unwrap();
        let report = s
            .run(&Controls {
                max_steps: Some(5),
                motion_tolerance: Some(0.0),
                energy_tolerance: Some(0.0),
                ..flow_controls()
            })
            .unwrap();
        assert_eq!(report.get("steps_taken"), Some(&Value::Int(5)));
        assert_eq!(
            report.get("stop_reason"),
            Some(&Value::Text("max_steps".into()))
        );
        let (b, a) = (
            report.float("energy_before").unwrap(),
            report.float("energy_after").unwrap(),
        );
        assert!(a < b, "energy {b} -> {a}");
    }

    #[test]
    fn ripple_wavelength_sets_the_length_weight() {
        let input = sphere();
        let mut s = CutSession::new(&input, 2, &circle()).unwrap();
        let c = Controls {
            ripple_wavelength: Some(0.25),
            ..Controls::default()
        };
        let omega = s.roughening_rate(&c).unwrap();
        let report = s.step(&c).unwrap();
        let alpha = report.float("length_weight").unwrap();
        assert!((alpha - vc::length_weight_for_ripples(omega, 0.25)).abs() <= 1e-12 * alpha.abs());
    }

    #[test]
    fn contradictory_options_are_input_errors() {
        let input = sphere();
        let mut s = CutSession::new(&input, 2, &circle()).unwrap();
        for c in [
            Controls {
                length_weight: Some(0.1),
                ripple_wavelength: Some(0.2),
                ..Controls::default()
            },
            Controls {
                step_rule: StepRule::Fixed,
                length_weight: Some(0.1),
                ..Controls::default()
            },
            Controls {
                scale_weights: Some(vec![1.0; 3]),
                length_weight: Some(0.1),
                ..Controls::default()
            },
            Controls {
                distortion_weights: Some(vec![1.0; 3]),
                length_weight: Some(0.1),
                ..Controls::default()
            },
            Controls {
                penalty: Some(vec![0.0; 3]),
                length_weight: Some(0.1),
                ..Controls::default()
            },
        ] {
            assert!(matches!(s.step(&c), Err(BindError::Input(_))), "{c:?}");
        }
    }

    #[test]
    fn more_than_two_regions_is_not_available() {
        let input = sphere();
        let e = CutSession::new(&input, 3, &circle()).unwrap_err();
        assert!(matches!(e, BindError::NotAvailable(_)), "{e}");
        assert!(e.message().contains("not available"), "{e}");
        let e = CutSession::new(&input, 1, &circle()).unwrap_err();
        assert!(matches!(e, BindError::Input(_)), "{e}");
    }

    #[test]
    fn bad_state_is_an_input_error_naming_the_problem() {
        let input = sphere();
        let mut s = CutSession::new(&input, 2, &circle()).unwrap();
        let e = s.set_state(&[0.0; 5]).unwrap_err();
        assert!(
            matches!(e, BindError::Input(_)) && e.message().contains("n_regions"),
            "{e}"
        );
        let mut state = s.state();
        let n = s.n_corners();
        state[0] = f64::NAN;
        let e = s.set_state(&state).unwrap_err();
        assert!(matches!(e, BindError::Input(_)), "{e}");
        assert_eq!(state.len(), 2 * n);
    }

    #[test]
    fn cut_description_inserts_and_unwraps() {
        let input = sphere();
        let mut s = CutSession::new(&input, 2, &circle()).unwrap();
        s.run(&Controls {
            max_steps: Some(3),
            motion_tolerance: Some(0.0),
            energy_tolerance: Some(0.0),
            ..flow_controls()
        })
        .unwrap();
        let cut = s.describe_cut(None).unwrap();
        let arrays = DescriptionArrays::from_cut(&cut);
        assert_eq!(arrays.edge_t.len(), cut.edge_points.len());
        assert_eq!(arrays.path_interior_starts.len(), cut.paths.len() + 1);
        let inserted = s.insert(Some(&cut), None).unwrap();
        let flat = InsertedArrays::from_inserted(&inserted);
        assert!(flat.polygon_starts.len() > input.n_polygons());
        assert_eq!(flat.vertex_source.len(), 3 * (flat.positions.len() / 3));
        let (_, unwrapped) = s.unwrap(None, &UnwrapSettings::default()).unwrap();
        let u = UnwrapArrays::from_unwrap(&unwrapped);
        assert_eq!(u.uvs.len(), u.inserted.corner_vertices.len());
        let preview = s.preview(None).unwrap();
        assert_eq!(preview.starts.len(), preview.closed.len() + 1);
        assert_eq!(
            preview.positions.len(),
            3 * *preview.starts.last().unwrap() as usize
        );
    }

    #[test]
    fn visibility_penalty_is_zero_facing_away() {
        let input = sphere();
        let eta = visibility_penalty(&input, [0.0, 0.0, 10.0], None).unwrap();
        assert_eq!(eta.len(), input.n_vertices());
        let top = eta.iter().cloned().fold(0.0, f64::max);
        assert!(top > 0.99 && eta.contains(&0.0), "{top}");
        let e = visibility_penalty(&input, [0.0, 0.0, 10.0], Some(&[1.0])).unwrap_err();
        assert!(matches!(e, BindError::Input(_)), "{e}");
    }
}
