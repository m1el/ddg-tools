//! Developability of triangle meshes: the binding of `ddg-developability`.
//!
//! Three functions, all on a polygon mesh with Blender's triangulation (the arrays of
//! [`UnwrapInput`]; the triangulation of the caller is used, never one chosen here):
//!
//! * [`developability`]: the per-vertex developability energy (a preview);
//! * [`make_developable`]: the developability flow of Stein, Grinspun, Crane (ACM TOG 2018,
//!   Sec. 4.2 and 5.2: rounds of L-BFGS with 4-1 subdivision between them and small-angle
//!   remeshing during them), SPEC §7.2 and decision D24: the result carries positions, and the
//!   triangle mesh only when the topology changed (or when asked for);
//! * [`developable_unwrap`]: cut through the vertices of high energy and flatten the charts with
//!   Boundary First Flattening (Sec. 4.4); the cut follows polygon edges only, so every polygon
//!   lies in one chart and has a UV per corner. Charts are packed into the unit square.
//!
//! The Python glue decodes the buffers; everything here is plain Rust and is exercised by
//! `cargo test -p ddg-py`. Mesh errors name vertices and polygons of the caller's mesh.
//!
//! # Face provenance
//!
//! The whole flow runs in one library call. The library reports the input triangle every output
//! triangle descends from (`FlowResult::face_parent`: exact through subdivision, flips and
//! collapses), and `face_polygon` is the polygon of that triangle.

use std::time::Instant;

use ddg_core::topology::{MeshError, NONE, TriMesh, compact_vertices};
use ddg_developability as dev;
use ddg_developability::{
    CutOptions, Energy, EnergyOptions, FlowOptions, LineSearch, Method, OptimizeStatus, RemeshOptions,
};

use crate::contract::UnwrapInput;
use crate::error::{InputError, UnwrapError};
use crate::pack::pack_charts;

/// Default free space around every chart when packing, as a fraction of the unit square (the
/// same as the unwrap of the add-on).
pub const DEFAULT_PACK_MARGIN: f64 = 1.0 / 256.0;

/// The names of the six energies, in the order of [`Energy::ALL`].
pub const ENERGY_NAMES: [&str; 6] =
    ["combinatorial", "combinatorial_max", "covariance", "covariance_max", "intrinsic", "intrinsic_max"];

/// The name of an energy.
pub fn energy_name(e: Energy) -> &'static str {
    ENERGY_NAMES[Energy::ALL.iter().position(|&x| x == e).unwrap_or(2)]
}

/// Parses an energy name.
///
/// # Errors
/// [`InputError::Rejected`] listing the valid names.
pub fn parse_energy(name: &'static str, value: &str) -> Result<Energy, UnwrapError> {
    match ENERGY_NAMES.iter().position(|&n| n == value) {
        Some(i) => Ok(Energy::ALL[i]),
        None => Err(InputError::Rejected {
            name,
            reason: format!("'{value}' is not an energy; use one of {}", ENERGY_NAMES.join(", ")),
        }
        .into()),
    }
}

// ---------------------------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------------------------

/// A mesh validation error in the caller's numbering.
fn mesh_message(e: &MeshError, new_to_old: &[u32], triangle_polygon: &[u32]) -> String {
    let v = |x: u32| new_to_old.get(x as usize).copied().unwrap_or(x);
    let p = |f: u32| triangle_polygon.get(f as usize).copied().unwrap_or(f);
    let edge = |e: &[u32; 2]| {
        let (a, b) = (v(e[0]), v(e[1]));
        (a.min(b), a.max(b))
    };
    let list = |polys: Vec<u32>| polys.iter().map(u32::to_string).collect::<Vec<_>>().join(", ");
    match e {
        MeshError::NoFaces => "the mesh has no faces".to_owned(),
        MeshError::RepeatedVertexInFace { face, vertices } => {
            let repeated = if vertices[0] == vertices[1] || vertices[0] == vertices[2] { vertices[0] } else { vertices[1] };
            format!("polygon {} uses vertex {} more than once in a triangle", p(*face), v(repeated))
        }
        MeshError::DuplicateFace { first, second, .. } => {
            format!("polygons {} and {} have the same vertices", p(*first), p(*second))
        }
        MeshError::NonManifoldEdge { edge: ed, faces } => {
            let mut polys: Vec<u32> = faces.iter().map(|&f| p(f)).collect();
            polys.sort_unstable();
            polys.dedup();
            let (a, b) = edge(ed);
            if polys.len() == faces.len() {
                format!(
                    "the mesh is not manifold: edge ({a}, {b}) is shared by {} faces: polygons {}",
                    faces.len(),
                    list(polys)
                )
            } else {
                format!(
                    "the mesh is not manifold: edge ({a}, {b}) is used by {} triangles of polygons {} \
                     (a polygon that touches itself)",
                    faces.len(),
                    list(polys)
                )
            }
        }
        MeshError::InconsistentOrientation { edge: ed, faces } => {
            let (a, b) = edge(ed);
            let (p0, p1) = (p(faces[0]), p(faces[1]));
            format!(
                "the mesh is not consistently oriented: polygons {} and {} traverse their shared edge \
                 ({a}, {b}) in the same direction (flip the normals so that they agree)",
                p0.min(p1),
                p0.max(p1)
            )
        }
        MeshError::NonManifoldVertex { vertex, fans } => format!(
            "the mesh is not manifold: vertex {} is shared by {fans} separate fans of polygons",
            v(*vertex)
        ),
        other => format!("the triangle mesh was refused: {other}"),
    }
}

/// What the library errors mean to the caller. `input` says that vertex and face numbers refer
/// to the prepared input mesh and can be translated.
fn dev_error(e: dev::DevError, prep: &Prepared, input: bool) -> UnwrapError {
    use dev::DevError as E;
    match e {
        E::Mesh(m) => UnwrapError::Mesh(mesh_message(&m, &prep.new_to_old, &prep.triangle_polygon)),
        E::Geometry(g) => UnwrapError::Mesh(format!("positions: {g}")),
        E::DegenerateFace { face } if input => UnwrapError::Mesh(format!(
            "polygon {} has a degenerate triangle (zero area): the energy needs the normal of every face",
            prep.triangle_polygon.get(face as usize).copied().unwrap_or(face)
        )),
        E::DegenerateFace { face } => UnwrapError::Mesh(format!(
            "a degenerate triangle (zero area) appeared during the computation (face {face} of the \
             current mesh)"
        )),
        E::NearDegenerateFace { face, min_angle, limit } if input => UnwrapError::Mesh(format!(
            "polygon {} has a nearly degenerate triangle (corner angle {:.3e} degrees, the limit is \
             {:.3e} degrees): repair the sliver or lower min_face_angle (0 switches the check off)",
            prep.triangle_polygon.get(face as usize).copied().unwrap_or(face),
            min_angle.to_degrees(),
            limit.to_degrees()
        )),
        E::NearDegenerateFace { face, min_angle, limit } => UnwrapError::Solve(format!(
            "a nearly degenerate triangle appeared (face {face} of the current mesh, corner angle \
             {:.3e} degrees, limit {:.3e} degrees)",
            min_angle.to_degrees(),
            limit.to_degrees()
        )),
        E::UndefinedEnergy { vertex, reason } if input => UnwrapError::Mesh(format!(
            "the energy is undefined at vertex {}: {reason}",
            prep.new_to_old.get(vertex as usize).copied().unwrap_or(vertex)
        )),
        E::UndefinedEnergy { vertex, reason } => UnwrapError::Mesh(format!(
            "the energy became undefined at vertex {vertex} of the current mesh: {reason}"
        )),
        E::NonFinitePosition { vertex } => {
            UnwrapError::Mesh(format!("the position of vertex {vertex} is not finite"))
        }
        E::InvalidOption { name, value, requirement } => UnwrapError::Input(InputError::BadOption {
            name,
            value: value.to_string(),
            requirement,
        }),
        E::VertexOutOfRange { what, vertex, n_vertices } => UnwrapError::Input(InputError::Rejected {
            name: what,
            reason: format!("vertex {vertex} is not a vertex of the mesh ({n_vertices} vertices)"),
        }),
        E::LengthMismatch { what, expected, got } => {
            UnwrapError::Solve(format!("internal error: {what} has {got} entries, expected {expected}"))
        }
        E::Bff(b) => UnwrapError::Solve(format!("boundary first flattening failed: {b}")),
        E::Eig(g) => UnwrapError::Solve(format!("eigensolver failed: {g}")),
        E::Internal(m) => UnwrapError::Solve(format!("internal error: {m}")),
        // Variants added by the library later: its message is meant for the user.
        #[allow(unreachable_patterns)]
        other => UnwrapError::Solve(other.to_string()),
    }
}

// ---------------------------------------------------------------------------------------------
// The triangle mesh of the loop triangles
// ---------------------------------------------------------------------------------------------

/// The triangle mesh of the caller's loop triangles, without unused vertices.
#[derive(Debug, Clone)]
pub struct Prepared {
    /// The mesh; vertex `i` is vertex `new_to_old[i]` of the caller.
    pub mesh: TriMesh,
    /// Its positions.
    pub positions: Vec<[f64; 3]>,
    /// For each vertex, the caller's vertex.
    pub new_to_old: Vec<u32>,
    /// For each caller's vertex, its vertex here, or [`NONE`] if no polygon uses it.
    pub old_to_new: Vec<u32>,
    /// The polygon of every triangle (the triangles are those of the input, in order).
    pub triangle_polygon: Vec<u32>,
    /// The vertices of the caller's mesh.
    pub n_input_vertices: usize,
}

/// Builds the triangle mesh of `input` and validates it (manifold, consistently oriented).
///
/// # Errors
/// [`UnwrapError::Mesh`] naming vertices and polygons of the caller.
pub fn prepare(input: &UnwrapInput) -> Result<Prepared, UnwrapError> {
    let cv = input.corner_vertices();
    let corner_polygon = input.corner_polygons();
    let faces: Vec<[u32; 3]> = input.triangles().iter().map(|t| t.map(|c| cv[c as usize])).collect();
    let triangle_polygon: Vec<u32> = input.triangles().iter().map(|t| corner_polygon[t[0] as usize]).collect();
    let n = input.n_vertices();
    let compact = compact_vertices(n, &faces).map_err(|e| UnwrapError::Mesh(e.to_string()))?;
    let mesh = TriMesh::new(compact.new_to_old.len(), compact.faces)
        .map_err(|e| UnwrapError::Mesh(mesh_message(&e, &compact.new_to_old, &triangle_polygon)))?;
    let positions = compact.new_to_old.iter().map(|&o| input.positions()[o as usize]).collect();
    Ok(Prepared {
        mesh,
        positions,
        new_to_old: compact.new_to_old,
        old_to_new: compact.old_to_new,
        triangle_polygon,
        n_input_vertices: n,
    })
}

fn check_vertices(name: &'static str, vertices: &[u32], n: usize) -> Result<(), UnwrapError> {
    match vertices.iter().position(|&v| v as usize >= n) {
        Some(position) => Err(InputError::IndexOutOfRange {
            name,
            position,
            value: vertices[position],
            count: n,
            of: "vertices",
        }
        .into()),
        None => Ok(()),
    }
}

/// Caller's vertices to vertices of the prepared mesh (vertices no polygon uses are dropped).
fn to_prepared(prep: &Prepared, vertices: &[u32]) -> Vec<u32> {
    let mut out: Vec<u32> = vertices.iter().map(|&v| prep.old_to_new[v as usize]).filter(|&v| v != NONE).collect();
    out.sort_unstable();
    out.dedup();
    out
}

// ---------------------------------------------------------------------------------------------
// The energy
// ---------------------------------------------------------------------------------------------

/// Options of the energy (SPEC §3).
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyParams {
    /// Which energy.
    pub energy: Energy,
    /// Count boundary vertices.
    pub include_boundary: bool,
    /// Leave out interior vertices with exactly three faces.
    pub omit_valence_three: bool,
    /// Worker threads of the evaluation (`None`: all).
    pub threads: Option<usize>,
}

impl Default for EnergyParams {
    fn default() -> Self {
        let d = EnergyOptions::default();
        EnergyParams {
            energy: d.energy,
            include_boundary: d.include_boundary,
            omit_valence_three: d.omit_valence_three,
            threads: d.threads,
        }
    }
}

impl EnergyParams {
    fn options(&self) -> EnergyOptions {
        EnergyOptions {
            energy: self.energy,
            include_boundary: self.include_boundary,
            omit_valence_three: self.omit_valence_three,
            threads: self.threads,
        }
    }
}

/// The result of [`developability`].
#[derive(Debug, Clone, PartialEq)]
pub struct DevelopabilityOutput {
    /// The energy at every vertex of the caller's mesh (0 for a vertex no polygon uses, and for
    /// a vertex that does not count, see `counted`).
    pub vertex_energy: Vec<f64>,
    /// Whether the vertex counts toward the total.
    pub counted: Vec<bool>,
    /// The total energy.
    pub energy: f64,
    /// Number of vertices whose active piece of the energy was a near tie.
    pub near_ties: usize,
    /// Wall time of the evaluation, seconds.
    pub seconds: f64,
}

/// The per-vertex developability energy of the mesh as it is (nothing is modified).
///
/// # Errors
/// [`UnwrapError::Mesh`] for a mesh that is not a valid manifold triangle mesh, a degenerate
/// triangle, or a vertex where the energy is undefined; [`UnwrapError::Input`] for a bad option.
pub fn developability(input: &UnwrapInput, params: &EnergyParams) -> Result<DevelopabilityOutput, UnwrapError> {
    let prep = prepare(input)?;
    let start = Instant::now();
    let ev = dev::evaluate(&prep.mesh, &prep.positions, &params.options()).map_err(|e| dev_error(e, &prep, true))?;
    let seconds = start.elapsed().as_secs_f64();
    let mut vertex_energy = vec![0.0; prep.n_input_vertices];
    let mut counted = vec![false; prep.n_input_vertices];
    for (i, &o) in prep.new_to_old.iter().enumerate() {
        vertex_energy[o as usize] = ev.vertex_energy[i];
        counted[o as usize] = ev.counted[i];
    }
    Ok(DevelopabilityOutput { vertex_energy, counted, energy: ev.energy, near_ties: ev.near_ties, seconds })
}

// ---------------------------------------------------------------------------------------------
// The flow
// ---------------------------------------------------------------------------------------------

/// Options of [`make_developable`].
#[derive(Debug, Clone, PartialEq)]
pub struct MakeOptions {
    /// The energy.
    pub energy: EnergyParams,
    /// L-BFGS, gradient descent or fixed step.
    pub method: Method,
    /// Line search of L-BFGS and gradient descent.
    pub line_search: LineSearch,
    /// Iterations per round.
    pub iterations: usize,
    /// Number of 4-1 subdivisions (rounds after the first).
    pub subdivisions: usize,
    /// The energy of round 0 when `subdivisions > 0` (the paper's coarse phase, Sec. 5.2, D36);
    /// the rounds after a subdivision use `energy`. `None`: `energy` in every round (the
    /// `FlowOptions::uniform` flow). Ignored when `subdivisions` is 0 (the one round uses
    /// `energy`).
    pub coarse_energy: Option<Energy>,
    /// Small-angle remeshing; `None` switches it off.
    pub remesh: Option<RemeshOptions>,
    /// Step cap in mean edge lengths.
    pub max_displacement: f64,
    /// Stop when `ℓ̄ ‖g‖_∞` is below this.
    pub gradient_tolerance: f64,
    /// Relative energy decrease that counts as a stall.
    pub energy_tolerance: f64,
    /// Consecutive stalled iterations that end a round.
    pub stall_iterations: usize,
    /// First trial step in mean edge lengths.
    pub initial_step: f64,
    /// Admissibility floor on the smallest corner angle of a face, radians (D31); 0 is off.
    pub min_face_angle: f64,
    /// Vertices of the caller that do not move.
    pub pinned: Vec<u32>,
    /// Return the triangle mesh even if the topology did not change.
    pub return_faces: bool,
    /// Call the progress function at most this often (seconds); `0` calls it after every
    /// accepted iteration.
    pub progress_interval: f64,
}

impl Default for MakeOptions {
    fn default() -> Self {
        let o = dev::OptimizeOptions::default();
        MakeOptions {
            energy: EnergyParams::default(),
            method: o.method,
            line_search: o.line_search,
            iterations: 500,
            subdivisions: 2,
            coarse_energy: Some(Energy::Combinatorial),
            remesh: Some(RemeshOptions::default()),
            max_displacement: o.max_displacement,
            gradient_tolerance: o.gradient_tolerance,
            energy_tolerance: o.energy_tolerance,
            stall_iterations: o.stall_iterations,
            initial_step: o.initial_step,
            min_face_angle: o.min_face_angle,
            pinned: Vec::new(),
            return_faces: false,
            progress_interval: 0.1,
        }
    }
}

/// What the progress function is told.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProgressInfo {
    /// Round (0 is the first).
    pub round: usize,
    /// Accepted iterations in this round.
    pub iteration: usize,
    /// Current energy.
    pub energy: f64,
    /// `ℓ̄ ‖g‖_∞`.
    pub gradient_norm: f64,
}

/// One round of the flow, summarized.
#[derive(Debug, Clone, PartialEq)]
pub struct RoundSummary {
    /// Vertices at the end of the round.
    pub n_vertices: usize,
    /// Triangles at the end of the round.
    pub n_faces: usize,
    /// Energy at the start (after the subdivision).
    pub energy_initial: f64,
    /// Energy at the end.
    pub energy_final: f64,
    /// Accepted iterations (over restarts after remeshing).
    pub iterations: usize,
    /// Why the optimizer ended (`gradient_tolerance`, `energy_stalled`, `max_iterations`,
    /// `line_search_failed`, `stationary`, `cancelled`, `non_descent`).
    pub status: &'static str,
    /// Edge flips of all remeshing runs of the round.
    pub flips: usize,
    /// Edge collapses of all remeshing runs of the round.
    pub collapses: usize,
    /// The energy minimized in this round (its name, see [`ENERGY_NAMES`]).
    pub energy: &'static str,
    /// Of `flips`, those applied as a fallback (D35).
    pub fallback_flips: usize,
    /// Of `collapses`, those applied as a fallback (D35).
    pub fallback_collapses: usize,
    /// Remeshing runs of the round that the library rejected as a whole (it keeps the mesh).
    pub remesh_rejected: usize,
    /// Trial steps of the round that the admissibility floor refused (D31).
    pub inadmissible_trials: usize,
    /// Wall time of the round.
    pub seconds: f64,
}

/// Why the flow ended before a round (D32): the energy of that round is undefined at its
/// starting state.
#[derive(Debug, Clone, PartialEq)]
pub struct StopInfo {
    /// The round that was not run (0 is the first).
    pub round: usize,
    /// Whether the undefined state is the result of that round's subdivision.
    pub after_subdivision: bool,
    /// The vertex where the energy is undefined, in the numbering of the mesh at that moment
    /// (after subdivisions and remeshing; not the caller's).
    pub vertex: u32,
    /// Why it is undefined.
    pub reason: &'static str,
}

/// The report of [`make_developable`].
#[derive(Debug, Clone, PartialEq)]
pub struct MakeReport {
    /// `Some` if the flow ended before a round because the energy is undefined there; the result
    /// is the state after the previous round.
    pub stopped: Option<StopInfo>,
    /// Triangles of the result with a corner angle below the small angle of the remeshing
    /// (1 degree if remeshing is off) that nothing could repair.
    pub bad_faces: usize,
    /// Trial steps of all rounds refused by the admissibility floor.
    pub inadmissible_trials: usize,
    /// Energy before the flow.
    pub energy_initial: f64,
    /// Energy after the flow.
    pub energy_final: f64,
    /// The progress function stopped the flow; the state at that moment is returned.
    pub cancelled: bool,
    /// Whether the faces of the result differ from the input's triangles (subdivision or an
    /// applied remeshing operation).
    pub topology_changed: bool,
    /// Vertices and triangles of the input (vertices no polygon uses are not counted among the
    /// triangles' vertices but are in `n_input_vertices`).
    pub n_input_vertices: usize,
    /// Triangles of the input.
    pub n_input_triangles: usize,
    /// Vertices of the result.
    pub n_vertices: usize,
    /// Triangles of the result.
    pub n_triangles: usize,
    /// Total accepted iterations.
    pub iterations: usize,
    /// Largest displacement of a vertex of the input from its start, in mean edge lengths, over
    /// the last round (the library reports it per round).
    pub max_displacement: f64,
    /// Mean edge length of the input (the unit of the steps).
    pub mean_edge_length: f64,
    /// Smallest corner angle of the result, radians.
    pub min_angle: f64,
    /// Largest non-planarity of a polygon with more than three corners, before the flow: the
    /// largest distance of a corner from the polygon's mean plane (centroid and Newell normal).
    /// `None` if all polygons are triangles.
    pub nonplanarity_before: Option<f64>,
    /// The same after the flow, `None` if all polygons are triangles or the topology changed.
    pub nonplanarity_after: Option<f64>,
    /// `nonplanarity_after` relative to the mean side length of that polygon.
    pub nonplanarity_after_relative: Option<f64>,
    /// The polygon with the largest non-planarity after the flow.
    pub nonplanarity_polygon: Option<u32>,
    /// One entry per round that ran.
    pub rounds: Vec<RoundSummary>,
    /// Warnings.
    pub warnings: Vec<String>,
    /// Wall time, seconds.
    pub seconds: f64,
}

/// The result of [`make_developable`].
#[derive(Debug, Clone, PartialEq)]
pub struct MakeOutput {
    /// 3 reals per vertex. If the topology did not change: one per vertex of the caller's mesh
    /// (a vertex no polygon uses keeps its position). Otherwise one per vertex of the new mesh.
    pub positions: Vec<f64>,
    /// 3 vertex indices per triangle of the new mesh; `None` if the topology did not change
    /// (and `return_faces` was off).
    pub faces: Option<Vec<u32>>,
    /// The polygon of the caller's mesh every triangle of `faces` descends from; present with
    /// `faces`.
    pub face_polygon: Option<Vec<u32>>,
    /// For every vertex of the result, the vertex of the caller's mesh it descends from, or
    /// [`NONE`] for a vertex created by subdivision.
    pub vertex_source: Vec<u32>,
    /// What happened.
    pub report: MakeReport,
}

fn status_name(s: OptimizeStatus) -> &'static str {
    match s {
        OptimizeStatus::GradientTolerance => "gradient_tolerance",
        OptimizeStatus::EnergyStalled => "energy_stalled",
        OptimizeStatus::MaxIterations => "max_iterations",
        OptimizeStatus::LineSearchFailed => "line_search_failed",
        OptimizeStatus::Cancelled => "cancelled",
        OptimizeStatus::NonDescent => "non_descent",
        OptimizeStatus::Stationary => "stationary",
    }
}

fn newell(points: &[[f64; 3]]) -> [f64; 3] {
    let mut n = [0.0; 3];
    for i in 0..points.len() {
        let (a, b) = (points[i], points[(i + 1) % points.len()]);
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    n
}

/// The largest distance of a corner from the mean plane over the polygons with more than three
/// corners, with the polygon and its mean side length; `None` if there are none.
pub fn max_nonplanarity(input: &UnwrapInput, positions: &[[f64; 3]]) -> Option<(f64, u32, f64)> {
    let mut best: Option<(f64, u32, f64)> = None;
    let cv = input.corner_vertices();
    for p in 0..input.n_polygons() {
        let r = input.polygon_corners(p);
        if r.len() <= 3 {
            continue;
        }
        let pts: Vec<[f64; 3]> = cv[r].iter().map(|&v| positions[v as usize]).collect();
        let k = pts.len() as f64;
        let c = [0, 1, 2].map(|a| pts.iter().map(|q| q[a]).sum::<f64>() / k);
        let n = newell(&pts);
        let nl = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if nl.is_nan() || nl <= 0.0 {
            continue;
        }
        let d = pts
            .iter()
            .map(|q| ((q[0] - c[0]) * n[0] + (q[1] - c[1]) * n[1] + (q[2] - c[2]) * n[2]).abs() / nl)
            .fold(0.0, f64::max);
        let side = (0..pts.len())
            .map(|i| {
                let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
                ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
            })
            .sum::<f64>()
            / k;
        if best.is_none_or(|b| d > b.0) {
            best = Some((d, p as u32, side));
        }
    }
    best
}

/// Runs the developability flow on the triangle mesh of the input. See the
/// [module documentation](self).
///
/// `progress` is called (at most every `options.progress_interval` seconds, from the calling
/// thread) with the state of the optimizer; returning `false` stops the flow, and the state at
/// that moment is returned with `cancelled` set.
///
/// # Errors
/// [`UnwrapError::Mesh`] for an invalid mesh, a degenerate triangle or an energy that is
/// undefined at the input; [`UnwrapError::Input`] for options outside their domain or pinned
/// vertices that are not vertices; [`UnwrapError::Solve`] for a failure of the computation.
pub fn make_developable(
    input: &UnwrapInput,
    options: &MakeOptions,
    mut progress: Option<&mut dyn FnMut(&ProgressInfo) -> bool>,
) -> Result<MakeOutput, UnwrapError> {
    let start = Instant::now();
    let prep = prepare(input)?;
    check_vertices("pinned_vertices", &options.pinned, prep.n_input_vertices)?;
    let bad = |name, value: f64, requirement| -> UnwrapError {
        InputError::BadOption { name, value: value.to_string(), requirement }.into()
    };
    if options.iterations == 0 {
        return Err(bad("iterations", 0.0, ">= 1"));
    }
    let energy = options.energy.options();
    energy.validate().map_err(|e| dev_error(e, &prep, true))?;
    let oo = dev::OptimizeOptions {
        method: options.method,
        line_search: options.line_search,
        max_iterations: options.iterations,
        gradient_tolerance: options.gradient_tolerance,
        energy_tolerance: options.energy_tolerance,
        stall_iterations: options.stall_iterations,
        initial_step: options.initial_step,
        max_displacement: options.max_displacement,
        min_face_angle: options.min_face_angle,
        pinned: Vec::new(),
    };
    oo.validate().map_err(|e| {
        UnwrapError::Input(InputError::BadOption { name: e.name, value: e.value.to_string(), requirement: e.requirement })
    })?;
    if let Some(r) = &options.remesh {
        r.validate().map_err(|e| dev_error(e, &prep, true))?;
    }
    if !(options.progress_interval.is_finite() && options.progress_interval >= 0.0) {
        return Err(bad("progress_interval", options.progress_interval, "finite and >= 0"));
    }

    // The input state; this also finds degenerate triangles and undefined energies in the
    // caller's numbering.
    let energy_initial =
        dev::evaluate_energy(&prep.mesh, &prep.positions, &energy).map_err(|e| dev_error(e, &prep, true))?.0;
    if options.remesh.is_none() && options.min_face_angle > 0.0 {
        // Without remeshing the flow starts on the input itself: report in the caller's numbering.
        for (f, tri) in prep.mesh.faces().iter().enumerate() {
            let a = dev::geometry::corner_angles(tri.map(|v| prep.positions[v as usize]));
            let min_angle = a[0].min(a[1]).min(a[2]);
            if min_angle < options.min_face_angle {
                let e = dev::DevError::NearDegenerateFace { face: f as u32, min_angle, limit: options.min_face_angle };
                return Err(dev_error(e, &prep, true));
            }
        }
    }
    let nonplanarity_before = max_nonplanarity(input, input.positions()).map(|b| b.0);

    let fo = FlowOptions {
        optimize: oo,
        remesh: options.remesh.clone(),
        pinned: to_prepared(&prep, &options.pinned),
        threads: options.energy.threads,
        ..FlowOptions::uniform(energy.clone(), options.subdivisions, options.iterations)
    };
    let mut fo = fo;
    if let (Some(coarse), true) = (options.coarse_energy, options.subdivisions > 0) {
        fo.rounds[0].energy = EnergyOptions { energy: coarse, ..energy.clone() };
    }
    let round_energies: Vec<&'static str> = fo.rounds.iter().map(|r| energy_name(r.energy.energy)).collect();
    let interval = options.progress_interval;
    let mut last_call: Option<Instant> = None;
    let mut cb = |p: &dev::Progress| -> dev::Control {
        let Some(user) = progress.as_deref_mut() else { return dev::Control::Continue };
        let now = Instant::now();
        if let Some(t) = last_call
            && now.duration_since(t).as_secs_f64() < interval
        {
            return dev::Control::Continue;
        }
        last_call = Some(now);
        let info = ProgressInfo { round: p.round, iteration: p.iteration, energy: p.energy, gradient_norm: p.gradient_norm };
        if user(&info) { dev::Control::Continue } else { dev::Control::Break }
    };
    let res = dev::flow(&prep.mesh, &prep.positions, &fo, Some(&mut cb)).map_err(|e| dev_error(e, &prep, false))?;
    let rounds: Vec<RoundSummary> = res
        .report
        .rounds
        .iter()
        .enumerate()
        .map(|(i, rr)| RoundSummary {
            energy: round_energies[i],
            fallback_flips: rr.remesh.iter().map(|x| x.fallback_flips).sum(),
            fallback_collapses: rr.remesh.iter().map(|x| x.fallback_collapses).sum(),
            n_vertices: rr.n_vertices,
            n_faces: rr.n_faces,
            energy_initial: rr.energy_after_subdivision,
            energy_final: rr.optimize.energy_final,
            iterations: rr.optimize.iterations,
            status: status_name(rr.optimize.status),
            flips: rr.remesh.iter().map(|x| x.flips).sum(),
            collapses: rr.remesh.iter().map(|x| x.collapses).sum(),
            remesh_rejected: rr.remesh.iter().filter(|x| x.rejected.is_some()).count(),
            inadmissible_trials: rr.optimize.inadmissible_trials,
            seconds: rr.seconds_subdivide + rr.seconds_remesh + rr.seconds_optimize,
        })
        .collect();
    let (max_disp, mean_edge) = res
        .report
        .rounds
        .last()
        .map_or((0.0, 0.0), |rr| (rr.optimize.max_vertex_displacement, rr.optimize.mean_edge_length));
    let cancelled = res.report.cancelled;
    let stopped = res.report.stopped.as_ref().map(|s| StopInfo {
        round: s.round,
        after_subdivision: s.after_subdivision,
        vertex: s.vertex,
        reason: s.reason,
    });
    let bad_faces = res.report.bad_faces;
    let face_polygon: Vec<u32> = res.face_parent.iter().map(|&t| prep.triangle_polygon[t as usize]).collect();
    let source = res.vertex_source;
    let (mesh, pos) = (res.mesh, res.positions);
    let changed = mesh.n_vertices() != prep.mesh.n_vertices()
        || mesh.faces() != prep.mesh.faces()
        || source.iter().enumerate().any(|(i, &s)| s as usize != i);
    let energy_final = dev::evaluate_energy(&mesh, &pos, &energy).map_err(|e| dev_error(e, &prep, false))?.0;

    let mut warnings = Vec::new();
    if cancelled {
        warnings.push("the flow was stopped before it finished".to_owned());
    }
    if let Some(st) = &stopped {
        warnings.push(format!(
            "the flow stopped before round {}{}: the energy is undefined at vertex {} of the mesh at that \
             moment ({}); the result is the state after round {}",
            st.round + 1,
            if st.after_subdivision { " (right after its subdivision)" } else { "" },
            st.vertex,
            st.reason,
            st.round
        ));
    }
    if bad_faces > 0 {
        warnings.push(format!(
            "{bad_faces} triangles have a corner angle below the small angle and could not be repaired"
        ));
    }
    let (positions, faces, face_poly, vertex_source);
    let mut np_after = None;
    if !changed {
        // Same triangles on the same vertices: write the moved positions back by caller's id.
        let mut all: Vec<[f64; 3]> = input.positions().to_vec();
        for (i, &o) in prep.new_to_old.iter().enumerate() {
            all[o as usize] = pos[i];
        }
        np_after = max_nonplanarity(input, &all);
        positions = all.iter().flatten().copied().collect();
        vertex_source = (0..prep.n_input_vertices as u32).collect();
        if options.return_faces {
            faces = Some(prep.mesh.faces().iter().flatten().copied().map(|v| prep.new_to_old[v as usize]).collect());
            face_poly = Some(prep.triangle_polygon.clone());
        } else {
            faces = None;
            face_poly = None;
        }
    } else {
        positions = pos.iter().flatten().copied().collect();
        vertex_source = source.iter().map(|&s| if s == NONE { NONE } else { prep.new_to_old[s as usize] }).collect();
        faces = Some(mesh.faces().iter().flatten().copied().collect());
        face_poly = Some(face_polygon);
    }
    let report = MakeReport {
        stopped,
        bad_faces,
        inadmissible_trials: rounds.iter().map(|x| x.inadmissible_trials).sum(),
        energy_initial,
        energy_final,
        cancelled,
        topology_changed: changed,
        n_input_vertices: prep.n_input_vertices,
        n_input_triangles: prep.mesh.n_faces(),
        n_vertices: mesh.n_vertices(),
        n_triangles: mesh.n_faces(),
        iterations: rounds.iter().map(|x| x.iterations).sum(),
        max_displacement: max_disp,
        mean_edge_length: mean_edge,
        min_angle: mesh
            .faces()
            .iter()
            .flat_map(|f| dev::geometry::corner_angles(f.map(|v| pos[v as usize])))
            .fold(f64::INFINITY, f64::min),
        nonplanarity_before,
        nonplanarity_after: np_after.map(|b| b.0),
        nonplanarity_after_relative: np_after.map(|b| if b.2 > 0.0 { b.0 / b.2 } else { 0.0 }),
        nonplanarity_polygon: np_after.map(|b| b.1),
        rounds,
        warnings,
        seconds: start.elapsed().as_secs_f64(),
    };
    Ok(MakeOutput { positions, faces, face_polygon: face_poly, vertex_source, report })
}

// ---------------------------------------------------------------------------------------------
// Cut and flatten
// ---------------------------------------------------------------------------------------------

/// Options of [`developable_unwrap`] (`CutOptions` of the library, with the pack margin).
#[derive(Debug, Clone, PartialEq)]
pub struct UnwrapDevOptions {
    /// The cut options; `may_cut` and `must_cut` are filled in from the polygons and the seams.
    pub cut: CutOptions,
    /// Free space around every chart when packing, `0 <= margin < 0.5`.
    pub margin: f64,
}

impl Default for UnwrapDevOptions {
    fn default() -> Self {
        UnwrapDevOptions { cut: CutOptions::default(), margin: DEFAULT_PACK_MARGIN }
    }
}

/// One chart of [`developable_unwrap`].
#[derive(Debug, Clone, PartialEq)]
pub struct ChartInfo {
    /// Lowest-numbered polygon of the chart.
    pub first_polygon: u32,
    /// Number of polygons.
    pub n_polygons: usize,
    /// Triangles whose image is flipped.
    pub flipped_triangles: usize,
    /// Triangles whose image has zero area.
    pub degenerate_triangles: usize,
    /// Boundary vertices around which the layout covers more than a full turn.
    pub wedges_over_2pi: usize,
    /// Crossings of the boundary of the layout with itself.
    pub boundary_crossings: usize,
    /// Sampled fraction of the chart's footprint that is covered more than once.
    pub overlapped_fraction: Option<f64>,
    /// No flips, no degenerate triangles, no wedges over a full turn, no crossings.
    pub is_injective: bool,
    /// Quasi-conformal error, mean (1 is conformal).
    pub quasi_conformal_mean: f64,
    /// Quasi-conformal error, largest.
    pub quasi_conformal_max: f64,
    /// Smallest and largest log scale factor over the triangles.
    pub log_scale_range: [f64; 2],
}

/// The report of [`developable_unwrap`].
#[derive(Debug, Clone, PartialEq)]
pub struct UnwrapDevReport {
    /// The energy that chose the terminals.
    pub energy: &'static str,
    /// Number of charts.
    pub n_charts: usize,
    /// Number of edges cut.
    pub n_seam_edges: usize,
    /// Terminals from the thresholds and the extra terminals.
    pub n_terminals: usize,
    /// Terminals the library added on closed components that had fewer than two.
    pub n_added_terminals: usize,
    /// The largest per-vertex energy.
    pub max_vertex_energy: f64,
    /// Flipped triangles, all charts.
    pub flipped_triangles: usize,
    /// Zero-area triangles, all charts.
    pub degenerate_triangles: usize,
    /// No chart has a flip, a degenerate triangle, an over-full wedge or a boundary crossing.
    pub is_injective: bool,
    /// Area-weighted mean of `|log scale|` over the triangles: 0 for an isometric development.
    pub mean_abs_log_scale: f64,
    /// The charts.
    pub charts: Vec<ChartInfo>,
    /// Scale applied by the packing.
    pub pack_scale: f64,
    /// Margin used by the packing.
    pub pack_margin: f64,
    /// Fraction of the unit square covered by the chart bounding boxes.
    pub pack_coverage: f64,
    /// User seams that lie on the boundary of the mesh (nothing to cut).
    pub seams_on_boundary: usize,
    /// Warnings.
    pub warnings: Vec<String>,
    /// Wall time, seconds.
    pub seconds: f64,
}

/// The result of [`developable_unwrap`].
#[derive(Debug, Clone, PartialEq)]
pub struct UnwrapDevOutput {
    /// A UV per polygon corner, packed into the unit square.
    pub corner_uvs: Vec<[f64; 2]>,
    /// The chart of every polygon.
    pub polygon_chart: Vec<u32>,
    /// The log scale factor `½ ln(A_uv / A)` of the triangle that owns the corner (`NaN` for a
    /// triangle with zero area in the plane).
    pub corner_log_scale: Vec<f64>,
    /// The edges that were cut, as vertex pairs `[min, max]` of the caller's mesh, ascending.
    pub seam_edges: Vec<[u32; 2]>,
    /// The terminals, vertices of the caller's mesh, ascending.
    pub terminals: Vec<u32>,
    /// The per-vertex energy that selected the terminals (every vertex of the caller's mesh; 0
    /// for an unused vertex).
    pub vertex_energy: Vec<f64>,
    /// The angle defect `Ω_v` of every interior vertex (0 on the boundary and for unused vertices).
    pub vertex_defect: Vec<f64>,
    /// The star embeddedness diagnostic of the library (`NaN` for unused vertices).
    pub star_embeddedness: Vec<f64>,
    /// What happened.
    pub report: UnwrapDevReport,
}

/// Cuts the mesh through the vertices of high energy (polygon edges only) and flattens the
/// charts; `input.seam_edges()` are cut as well. See the [module documentation](self).
///
/// `extra_terminals` are vertices of the caller's mesh that become terminals as well.
///
/// # Errors
/// As for [`make_developable`]; also a seam that is not a polygon edge of the mesh, an extra
/// terminal that is not a vertex, and an invalid margin.
pub fn developable_unwrap(
    input: &UnwrapInput,
    extra_terminals: &[u32],
    options: &UnwrapDevOptions,
) -> Result<UnwrapDevOutput, UnwrapError> {
    let start = Instant::now();
    let prep = prepare(input)?;
    check_vertices("extra_terminals", extra_terminals, prep.n_input_vertices)?;
    let mesh = &prep.mesh;

    // Polygon edges may be cut, diagonals of the triangulation may not.
    let cv = input.corner_vertices();
    let mut may_cut = vec![false; mesh.n_edges()];
    for p in 0..input.n_polygons() {
        let r = input.polygon_corners(p);
        let n = r.len();
        for i in 0..n {
            let (a, b) = (cv[r.start + i], cv[r.start + (i + 1) % n]);
            let (a, b) = (prep.old_to_new[a as usize], prep.old_to_new[b as usize]);
            if let Some(e) = mesh.find_edge(a as usize, b as usize) {
                may_cut[e] = true;
            }
        }
    }
    let mut must_cut = Vec::new();
    let mut seams_on_boundary = 0;
    for &[a, b] in input.seam_edges() {
        let (na, nb) = (prep.old_to_new[a as usize], prep.old_to_new[b as usize]);
        let edge = if na == NONE || nb == NONE { None } else { mesh.find_edge(na as usize, nb as usize) };
        match edge {
            Some(e) if may_cut[e] => {
                if mesh.is_boundary_edge(e) {
                    seams_on_boundary += 1;
                } else {
                    must_cut.push(e as u32);
                }
            }
            _ => {
                return Err(InputError::Rejected {
                    name: "seam_edges",
                    reason: format!("({a}, {b}) is not an edge of any polygon of the mesh"),
                }
                .into());
            }
        }
    }
    let mut cut = options.cut.clone();
    cut.may_cut = Some(may_cut);
    cut.must_cut = must_cut;
    cut.extra_terminals = to_prepared(&prep, extra_terminals);
    cut.validate().map_err(|e| dev_error(e, &prep, true))?;

    let flat = dev::cut_and_flatten(mesh, &prep.positions, &cut).map_err(|e| dev_error(e, &prep, true))?;

    // Triangle corners to polygon corners. All triangles of a polygon share one chart, so a
    // corner has one UV whichever triangle supplies it.
    let n_corners = input.n_corners();
    let mut uvs = vec![[f64::NAN; 2]; n_corners];
    let mut log_scale = vec![f64::NAN; n_corners];
    let mut written = vec![false; n_corners];
    for (t, tri) in input.triangles().iter().enumerate() {
        for (k, &c) in tri.iter().enumerate() {
            let c = c as usize;
            if !written[c] {
                written[c] = true;
                uvs[c] = flat.corner_uvs[3 * t + k];
                log_scale[c] = flat.face_log_scale[t];
            }
        }
    }
    let mut polygon_chart = vec![u32::MAX; input.n_polygons()];
    for (t, &p) in prep.triangle_polygon.iter().enumerate() {
        let chart = flat.face_chart[t];
        let slot = &mut polygon_chart[p as usize];
        if *slot == u32::MAX {
            *slot = chart;
        } else if *slot != chart {
            return Err(UnwrapError::Solve(format!(
                "internal error: polygon {p} was split between charts {} and {chart}",
                *slot
            )));
        }
    }
    let pack = pack_charts(&mut uvs, input.polygon_offsets(), &polygon_chart, options.margin)?;

    let mut charts: Vec<ChartInfo> = Vec::with_capacity(flat.charts.len());
    let mut polygons_per_chart = vec![0usize; flat.charts.len()];
    let mut first_polygon = vec![u32::MAX; flat.charts.len()];
    for (p, &c) in polygon_chart.iter().enumerate() {
        polygons_per_chart[c as usize] += 1;
        first_polygon[c as usize] = first_polygon[c as usize].min(p as u32);
    }
    for (i, c) in flat.charts.iter().enumerate() {
        let d = &c.distortion;
        charts.push(ChartInfo {
            first_polygon: first_polygon[i],
            n_polygons: polygons_per_chart[i],
            flipped_triangles: d.flipped_triangles,
            degenerate_triangles: d.degenerate_triangles,
            wedges_over_2pi: c.overlap.wedges_over_2pi,
            boundary_crossings: c.overlap.boundary_crossings,
            overlapped_fraction: c.overlap.overlapped_fraction(),
            is_injective: c.is_injective(),
            quasi_conformal_mean: d.quasi_conformal.mean,
            quasi_conformal_max: d.quasi_conformal.max,
            log_scale_range: [d.log_scale.min, d.log_scale.max],
        });
    }
    let old = |v: u32| prep.new_to_old[v as usize];
    let seam_edges: Vec<[u32; 2]> = flat
        .cut_vertex_pairs
        .iter()
        .map(|&[a, b]| {
            let (a, b) = (old(a), old(b));
            [a.min(b), a.max(b)]
        })
        .collect();
    let mut terminals: Vec<u32> = flat.terminals.iter().chain(&flat.added_terminals).map(|&v| old(v)).collect();
    terminals.sort_unstable();
    terminals.dedup();
    let mut warnings = Vec::new();
    if !flat.is_injective() {
        warnings.push(format!(
            "the layout is not injective: {} flipped triangles, {} charts with overlap",
            flat.flipped_triangles(),
            charts.iter().filter(|c| !c.is_injective).count()
        ));
    }
    let report = UnwrapDevReport {
        energy: energy_name(cut.energy),
        n_charts: charts.len(),
        n_seam_edges: seam_edges.len(),
        n_terminals: flat.terminals.len(),
        n_added_terminals: flat.added_terminals.len(),
        max_vertex_energy: flat.vertex_energy.iter().copied().filter(|x| x.is_finite()).fold(0.0, f64::max),
        flipped_triangles: charts.iter().map(|c| c.flipped_triangles).sum(),
        degenerate_triangles: charts.iter().map(|c| c.degenerate_triangles).sum(),
        is_injective: flat.is_injective(),
        mean_abs_log_scale: flat.mean_abs_log_scale,
        charts,
        pack_scale: pack.scale,
        pack_margin: pack.margin,
        pack_coverage: pack.coverage,
        seams_on_boundary,
        warnings,
        seconds: start.elapsed().as_secs_f64(),
    };
    let spread = |values: &[f64], fill: f64| {
        let mut out = vec![fill; prep.n_input_vertices];
        for (i, &o) in prep.new_to_old.iter().enumerate() {
            out[o as usize] = values[i];
        }
        out
    };
    Ok(UnwrapDevOutput {
        corner_uvs: uvs,
        polygon_chart,
        corner_log_scale: log_scale,
        seam_edges,
        terminals,
        vertex_energy: spread(&flat.vertex_energy, 0.0),
        vertex_defect: spread(&flat.vertex_defect, 0.0),
        star_embeddedness: spread(&flat.star_embeddedness, f64::NAN),
        report,
    })
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::contract::MeshBuffers;
    use crate::contract::test_meshes::{buffers, cube, grid};

    fn input(b: MeshBuffers) -> UnwrapInput {
        UnwrapInput::from_buffers(b).unwrap()
    }

    /// A `n x n` grid of quads with a smooth bump: not developable.
    fn bump(n: u32) -> MeshBuffers {
        let mut b = grid(n, n);
        for (i, p) in b.positions.chunks_exact_mut(3).enumerate() {
            let (x, y) = ((i as u32 % (n + 1)) as f64 / n as f64 - 0.5, (i as u32 / (n + 1)) as f64 / n as f64 - 0.5);
            p[2] = 0.2 * (-(x * x + y * y) * 8.0).exp();
        }
        b
    }

    /// An open cylinder of `m` quads around, `h` rows, radius 1.
    fn cylinder(m: u32, h: u32) -> MeshBuffers {
        let mut pos = Vec::new();
        for j in 0..=h {
            for i in 0..m {
                let a = 2.0 * PI * i as f64 / m as f64;
                pos.push([a.cos(), a.sin(), j as f64 * 0.3]);
            }
        }
        let quads: Vec<[u32; 4]> = (0..h)
            .flat_map(|j| (0..m).map(move |i| (i, j)))
            .map(|(i, j)| [j * m + i, j * m + (i + 1) % m, (j + 1) * m + (i + 1) % m, (j + 1) * m + i])
            .collect();
        let polys: Vec<&[u32]> = quads.iter().map(|q| q.as_slice()).collect();
        buffers(&pos, &polys)
    }

    #[test]
    fn python_defaults_match_the_library() {
        // The literals in the Python signatures (src/python/developability.rs).
        let c = CutOptions::default();
        assert_eq!((c.epsilon, c.length_weight, c.alignment_weight), (1e-3, 1.0, 0.0));
        assert_eq!(c.defect_tolerance, Some(ddg_developability::DEFAULT_DEFECT_TOLERANCE));
        let o = MakeOptions::default();
        assert_eq!((o.iterations, o.subdivisions, o.max_displacement), (500, 2, 0.5));
        assert!(matches!(o.method, Method::Lbfgs { memory: 10 }));
        assert!(matches!(o.line_search, LineSearch::WeakWolfe { .. }));
        let r = o.remesh.unwrap();
        assert_eq!((r.small_angle, r.collapse_boundary_edges), (PI / 180.0, true));
        assert_eq!(EnergyParams::default().energy, Energy::Covariance);
        assert_eq!(o.coarse_energy, Some(Energy::Combinatorial));
        assert!(r.fallback_operations);
    }

    #[test]
    fn the_coarse_energy_is_the_energy_of_round_zero_and_none_is_the_uniform_flow() {
        let i = input(bump(4));
        let mut o = MakeOptions::default();
        o.iterations = 12;
        o.remesh = None;
        let names: Vec<&str> = make_developable(&i, &o, None).unwrap().report.rounds.iter().map(|r| r.energy).collect();
        assert_eq!(names, ["combinatorial", "covariance", "covariance"]);
        o.coarse_energy = None;
        let uniform = make_developable(&i, &o, None).unwrap();
        assert!(uniform.report.rounds.iter().all(|r| r.energy == "covariance"));
        // `None` is bitwise the library's uniform flow; the coarse energy changes the result.
        let fo = FlowOptions { remesh: None, ..FlowOptions::uniform(EnergyOptions::new(Energy::Covariance), 2, 12) };
        let prep = prepare(&i).unwrap();
        let lib = dev::flow(&prep.mesh, &prep.positions, &fo, None).unwrap();
        let flat: Vec<f64> = lib.positions.iter().flatten().copied().collect();
        assert_eq!(uniform.positions, flat);
        o.coarse_energy = Some(Energy::Combinatorial);
        assert_ne!(make_developable(&i, &o, None).unwrap().positions, flat);
        // With no subdivision the one round uses `energy`.
        o.subdivisions = 0;
        assert_eq!(make_developable(&i, &o, None).unwrap().report.rounds[0].energy, "covariance");
    }

    #[test]
    fn energy_names_round_trip() {
        for (i, e) in Energy::ALL.iter().enumerate() {
            assert_eq!(energy_name(*e), ENERGY_NAMES[i]);
            assert_eq!(parse_energy("energy", ENERGY_NAMES[i]).unwrap(), *e);
        }
        assert!(parse_energy("energy", "bogus").unwrap_err().to_string().contains("covariance_max"));
    }

    #[test]
    fn a_flat_grid_has_no_energy_and_a_bump_has() {
        let flat = developability(&input(grid(6, 6)), &EnergyParams::default()).unwrap();
        assert_eq!(flat.vertex_energy.len(), 49);
        assert!(flat.energy.abs() < 1e-12, "{}", flat.energy);
        let bumpy = developability(&input(bump(8)), &EnergyParams::default()).unwrap();
        assert!(bumpy.energy > 1e-4);
        assert!(bumpy.counted.iter().any(|&c| c));
        assert!(bumpy.vertex_energy.iter().all(|e| e.is_finite() && *e >= 0.0));
    }

    #[test]
    fn unused_vertices_are_allowed_and_get_zero() {
        let mut b = grid(2, 2);
        b.positions.extend([5.0, 5.0, 5.0]);
        let out = developability(&input(b), &EnergyParams::default()).unwrap();
        assert_eq!(out.vertex_energy.len(), 10);
        assert_eq!(out.vertex_energy[9], 0.0);
        assert!(!out.counted[9]);
    }

    #[test]
    fn non_manifold_and_unoriented_meshes_name_polygons() {
        let fins = crate::contract::test_meshes::three_fins();
        let e = developability(&input(fins), &EnergyParams::default()).unwrap_err();
        assert!(matches!(&e, UnwrapError::Mesh(m) if m.contains("not manifold") && m.contains("edge (0, 1)")), "{e}");
        let mut b = grid(2, 1);
        // Reverse the second quad.
        let n = b.corner_vertices.len();
        b.corner_vertices[n - 4..].reverse();
        b.triangles = crate::contract::test_meshes::fan_triangles(&b.polygon_starts, n as u32);
        let e = developability(&input(b), &EnergyParams::default()).unwrap_err();
        assert!(matches!(&e, UnwrapError::Mesh(m) if m.contains("polygons 0 and 1")), "{e}");
    }

    #[test]
    fn bad_options_are_input_errors() {
        let i = input(grid(3, 3));
        let mut o = MakeOptions::default();
        o.iterations = 0;
        assert!(matches!(make_developable(&i, &o, None), Err(UnwrapError::Input(_))));
        let mut o = MakeOptions::default();
        o.max_displacement = -1.0;
        let e = make_developable(&i, &o, None).unwrap_err();
        assert!(e.to_string().contains("max_displacement"), "{e}");
        let mut o = MakeOptions::default();
        o.pinned = vec![99];
        let e = make_developable(&i, &o, None).unwrap_err();
        assert!(e.to_string().starts_with("pinned_vertices"), "{e}");
        let mut o = MakeOptions::default();
        o.remesh = Some(RemeshOptions { small_angle: 2.0, ..RemeshOptions::default() });
        let e = make_developable(&i, &o, None).unwrap_err();
        assert!(e.to_string().contains("small_angle"), "{e}");
    }

    #[test]
    fn the_flow_lowers_the_energy_and_keeps_the_quads_when_the_topology_is_unchanged() {
        let i = input(bump(8));
        let mut o = MakeOptions::default();
        o.subdivisions = 0;
        o.iterations = 60;
        o.remesh = None;
        let out = make_developable(&i, &o, None).unwrap();
        assert!(out.faces.is_none() && out.face_polygon.is_none());
        assert!(!out.report.topology_changed);
        assert_eq!(out.positions.len(), 81 * 3);
        assert_eq!(out.vertex_source, (0..81).collect::<Vec<u32>>());
        assert!(out.report.energy_final < 0.5 * out.report.energy_initial, "{:?}", out.report);
        assert!(out.report.nonplanarity_after.is_some());
        assert!(out.report.nonplanarity_before.unwrap() >= 0.0);
        let mut o2 = o.clone();
        o2.return_faces = true;
        let out2 = make_developable(&i, &o2, None).unwrap();
        assert_eq!(out2.faces.as_ref().unwrap().len(), 3 * 128);
        assert_eq!(out2.face_polygon.as_ref().unwrap().len(), 128);
        assert_eq!(out2.positions, out.positions);
    }

    #[test]
    fn subdivision_returns_a_triangle_mesh_with_polygons_and_sources() {
        let i = input(bump(4));
        let mut o = MakeOptions::default();
        o.subdivisions = 1;
        o.iterations = 20;
        o.remesh = None;
        let out = make_developable(&i, &o, None).unwrap();
        let faces = out.faces.as_ref().unwrap();
        assert_eq!(faces.len(), 3 * 4 * 32);
        let fp = out.face_polygon.as_ref().unwrap();
        // Every child of a triangle belongs to the polygon of the triangle.
        assert!(fp.chunks_exact(4).all(|c| c.iter().all(|&p| p == c[0])));
        assert_eq!(*fp.iter().max().unwrap(), 15);
        let nv = out.positions.len() / 3;
        assert_eq!(out.vertex_source.len(), nv);
        assert!(out.vertex_source[..25].iter().enumerate().all(|(k, &s)| s as usize == k));
        assert!(out.vertex_source[25..].iter().all(|&s| s == NONE));
        assert!(out.report.topology_changed && out.report.rounds.len() == 2);
        assert_eq!(out.report.nonplanarity_after, None);
    }

    #[test]
    fn pinned_vertices_do_not_move() {
        let i = input(bump(6));
        let mut o = MakeOptions::default();
        o.subdivisions = 0;
        o.iterations = 30;
        o.remesh = None;
        o.pinned = vec![24, 25];
        let out = make_developable(&i, &o, None).unwrap();
        for v in [24usize, 25] {
            assert_eq!(out.positions[3 * v..3 * v + 3], i.positions()[v]);
        }
    }

    #[test]
    fn progress_can_stop_the_flow_and_the_flow_is_deterministic() {
        let i = input(bump(8));
        let mut o = MakeOptions::default();
        o.subdivisions = 1;
        o.iterations = 200;
        o.progress_interval = 0.0;
        let mut calls = 0;
        let mut stop = |p: &ProgressInfo| {
            calls += 1;
            p.iteration < 5
        };
        let out = make_developable(&i, &o, Some(&mut stop)).unwrap();
        assert!(out.report.cancelled);
        assert_eq!(out.report.rounds.len(), 1);
        assert_eq!(calls, 5);
        assert_eq!(out.report.rounds[0].status, "cancelled");
        o.iterations = 25;
        let a = make_developable(&i, &o, None).unwrap();
        let b = make_developable(&i, &o, None).unwrap();
        assert_eq!(a.positions, b.positions);
        assert_eq!(a.faces, b.faces);
    }

    #[test]
    fn a_cylinder_unwraps_with_one_seam_and_no_flips() {
        let i = input(cylinder(12, 4));
        let out = developable_unwrap(&i, &[], &UnwrapDevOptions::default()).unwrap();
        let r = &out.report;
        assert_eq!(out.corner_uvs.len(), 12 * 4 * 4);
        assert_eq!(r.flipped_triangles, 0, "{r:?}");
        assert!(r.is_injective, "{r:?}");
        assert_eq!(r.n_charts, 1);
        assert!(r.mean_abs_log_scale < 1e-6, "{}", r.mean_abs_log_scale);
        assert!(out.corner_uvs.iter().all(|uv| uv[0].is_finite() && (-1e-9..=1.0 + 1e-9).contains(&uv[0])));
        // The open cylinder is already a disk? No: an annulus, so the cut is needed.
        assert!(!out.seam_edges.is_empty());
        assert_eq!(out.polygon_chart, vec![0; 48]);
    }

    #[test]
    fn the_cube_is_cut_into_charts_along_polygon_edges() {
        let i = input(cube());
        let out = developable_unwrap(&i, &[], &UnwrapDevOptions::default()).unwrap();
        assert_eq!(out.report.flipped_triangles, 0, "{:?}", out.report);
        // Cube corners have defect pi/2 > 0.1 rad: all eight are terminals.
        assert_eq!(out.report.n_terminals, 8);
        assert!(out.seam_edges.len() >= 7);
        for e in &out.seam_edges {
            // Seams are cube edges (no quad diagonals).
            assert_eq!((e[0] ^ e[1]).count_ones(), 1, "{e:?}");
        }
        let mut opts = UnwrapDevOptions::default();
        opts.cut.defect_tolerance = None;
        let o2 = developable_unwrap(&i, &[], &opts).unwrap();
        assert_eq!(o2.report.flipped_triangles, 0);
    }

    #[test]
    fn seams_must_be_polygon_edges() {
        let mut b = cube();
        b.seam_edges = vec![0, 3]; // the diagonal of the z = 0 face
        let e = developable_unwrap(&input(b), &[], &UnwrapDevOptions::default()).unwrap_err();
        assert!(e.to_string().starts_with("seam_edges"), "{e}");
        let mut b = cube();
        b.seam_edges = vec![0, 1];
        let out = developable_unwrap(&input(b), &[], &UnwrapDevOptions::default()).unwrap();
        assert!(out.seam_edges.contains(&[0, 1]));
        let e = developable_unwrap(&input(cube()), &[77], &UnwrapDevOptions::default()).unwrap_err();
        assert!(e.to_string().starts_with("extra_terminals"), "{e}");
    }

    #[test]
    fn unwrap_is_deterministic_and_rejects_a_bad_margin() {
        let i = input(cylinder(8, 3));
        let a = developable_unwrap(&i, &[], &UnwrapDevOptions::default()).unwrap();
        let b = developable_unwrap(&i, &[], &UnwrapDevOptions::default()).unwrap();
        assert_eq!(a.corner_uvs, b.corner_uvs);
        let mut o = UnwrapDevOptions::default();
        o.margin = 0.7;
        assert!(developable_unwrap(&i, &[], &o).is_err());
    }

    #[test]
    fn face_polygon_is_exact_over_two_rounds_and_valid_after_remeshing() {
        let i = input(bump(4));
        let mut o = MakeOptions::default();
        o.subdivisions = 2;
        o.iterations = 15;
        o.remesh = None;
        let out = make_developable(&i, &o, None).unwrap();
        let fp = out.face_polygon.as_ref().unwrap();
        assert_eq!(fp.len(), 16 * 2 * 16);
        // 16 polygons, 2 triangles each, 16 descendants per triangle: 32 per polygon, in order.
        assert!(fp.chunks_exact(16).all(|c| c.iter().all(|&p| p == c[0])));
        let mut count = vec![0usize; 16];
        fp.iter().for_each(|&p| count[p as usize] += 1);
        assert!(count.iter().all(|&c| c == 32), "{count:?}");
        assert!(out.report.stopped.is_none() && out.report.rounds.len() == 3);
        // With an aggressive small angle the remeshing runs; the polygons stay valid.
        let mut o = MakeOptions::default();
        o.subdivisions = 1;
        o.iterations = 30;
        o.remesh = Some(RemeshOptions { small_angle: 0.8, min_face_angle: 0.0, ..RemeshOptions::default() });
        o.min_face_angle = 0.0;
        let out = make_developable(&i, &o, None).unwrap();
        let fp = out.face_polygon.as_ref().unwrap();
        assert_eq!(fp.len() * 3, out.faces.as_ref().unwrap().len());
        assert!(fp.iter().all(|&p| p < 16));
        assert_eq!(out.report.n_triangles, fp.len());
    }

    #[test]
    fn a_nearly_degenerate_input_is_a_mesh_error_naming_the_polygon() {
        let mut b = grid(2, 2);
        // The centre vertex (1, 1) moves next to the edge (0,0)-(1,0): needle triangles.
        b.positions[3 * 4] = 1.0;
        b.positions[3 * 4 + 1] = 1e-9;
        let mut o = MakeOptions::default();
        o.remesh = None;
        let m = make_developable(&input(b.clone()), &o, None).unwrap_err().to_string();
        assert!(m.contains("polygon") && m.contains("nearly degenerate") && m.contains("min_face_angle"), "{m}");
        o.min_face_angle = 0.0;
        assert!(make_developable(&input(b), &o, None).is_ok());
    }
}
