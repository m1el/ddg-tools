//! The PyO3 glue of `ddg_tools.developability`, `make_developable` and `developable_unwrap`.
//! As in the parent module no logic lives here: buffers are decoded, [`crate::developability`]
//! does the work with the GIL released (taking it back only to call the progress function and to
//! check for Ctrl-C), results go back as `array.array` objects and plain dictionaries.

use std::time::Instant;

use ddg_developability::{CutOptions, DEFAULT_DEFECT_TOLERANCE, FirstTrial, LineSearch, Method, RemeshOptions};
use ddg_bff::Extension;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use super::{InputError, index_array, indices, new_array, real_array, reals, unwrap_error, uv_array};
use crate::contract::{MeshBuffers, UnwrapInput};
use crate::developability as dv;

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(developability, m)?)?;
    m.add_function(wrap_pyfunction!(make_developable, m)?)?;
    m.add_function(wrap_pyfunction!(developable_unwrap, m)?)?;
    m.add_class::<DevelopabilityResult>()?;
    m.add_class::<MakeDevelopableResult>()?;
    m.add_class::<DevelopableUnwrapResult>()?;
    m.add_function(wrap_pyfunction!(developable_seams, m)?)?;
    m.add_class::<DevelopableSeamsResult>()?;
    m.add("ENERGIES", PyList::new(m.py(), dv::ENERGY_NAMES)?)?;
    m.add("DEFAULT_DEFECT_TOLERANCE", DEFAULT_DEFECT_TOLERANCE)?;
    Ok(())
}

/// `array.array('i')` from signed values; [`ddg_core::topology::NONE`] becomes -1.
fn signed_array<'py>(py: Python<'py>, name: &str, values: &[u32]) -> PyResult<Bound<'py, PyAny>> {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for (i, &v) in values.iter().enumerate() {
        let x = if v == ddg_core::topology::NONE {
            -1
        } else {
            i32::try_from(v).map_err(|_| {
                super::SolveError::new_err(format!("internal error: {name}[{i}] = {v} does not fit int32"))
            })?
        };
        bytes.extend_from_slice(&x.to_ne_bytes());
    }
    new_array(py, "i", &bytes)
}

fn bool_array<'py>(py: Python<'py>, values: &[bool]) -> PyResult<Bound<'py, PyAny>> {
    let bytes: Vec<u8> = values.iter().map(|&b| u8::from(b)).collect();
    new_array(py, "B", &bytes)
}

/// The mesh arguments, decoded and validated. Steering arguments are not part of these
/// functions; `seam_edges` is optional.
fn mesh<'py>(
    positions: &Bound<'py, PyAny>,
    corner_vertices: &Bound<'py, PyAny>,
    polygon_starts: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    seam_edges: Option<&Bound<'py, PyAny>>,
) -> PyResult<MeshBuffers> {
    Ok(MeshBuffers {
        positions: reals("positions", positions, 3, "vertex")?,
        corner_vertices: indices("corner_vertices", corner_vertices, 1, "corner")?,
        polygon_starts: indices("polygon_starts", polygon_starts, 1, "polygon")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
        seam_edges: match seam_edges {
            Some(obj) => indices("seam_edges", obj, 2, "edge")?,
            None => Vec::new(),
        },
        ..MeshBuffers::default()
    })
}

fn optional_indices<'py>(
    name: &'static str,
    obj: Option<&Bound<'py, PyAny>>,
) -> PyResult<Vec<u32>> {
    match obj {
        Some(obj) => indices(name, obj, 1, "vertex"),
        None => Ok(Vec::new()),
    }
}

fn threads(n: Option<i64>) -> PyResult<Option<usize>> {
    match n {
        None => Ok(None),
        Some(n) if n >= 1 => Ok(Some(n as usize)),
        Some(n) => Err(InputError::new_err(format!("threads: {n} must be >= 1 (or None)"))),
    }
}

fn count(name: &str, n: i64, min: i64) -> PyResult<usize> {
    if n < min {
        return Err(InputError::new_err(format!("{name}: {n} must be >= {min}")));
    }
    Ok(n as usize)
}

// -------------------------------------------------------------------------------------------
// developability
// -------------------------------------------------------------------------------------------

/// The result of `developability`.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct DevelopabilityResult {
    /// `array('d')`, one per vertex: the energy of the vertex (0 where it does not count).
    vertex_energy: Py<PyAny>,
    /// `array('B')`, one per vertex: 1 where the vertex counts toward the total.
    counted: Py<PyAny>,
    /// The total energy.
    energy: f64,
    /// `dict` of plain values.
    report: Py<PyDict>,
}

#[pymethods]
impl DevelopabilityResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "DevelopabilityResult(vertices={}, energy={:.6e})",
            self.vertex_energy.bind(py).len()?,
            self.energy
        ))
    }
}

/// The developability energy of a polygon mesh, per vertex. Returns a DevelopabilityResult.
///
/// This is the preview of the developability flow (Stein, Grinspun, Crane, "Developability of
/// Triangle Meshes", 2018): a vertex has energy 0 exactly when the faces around it are
/// developable there (its normals lie on a great circle, in the sense of the chosen energy). Nothing
/// is modified.
///
/// positions         3 reals per vertex.
/// corner_vertices   1 index per polygon corner ("loop"): its vertex.
/// polygon_starts    1 index per polygon: its first corner (Blender's loop_start).
/// triangles         3 corner indices per triangle: Blender's mesh.loop_triangles[...].loops.
///                   The triangulation of the caller is used.
/// energy            "covariance" (default; E_lambda, the smallest eigenvalue of the normal
///                   covariance), "combinatorial" (E_P), "combinatorial_max", "covariance_max",
///                   "intrinsic" or "intrinsic_max" (the module constant ENERGIES lists them).
/// include_boundary  count the vertices on the boundary (default True).
/// omit_valence_three  leave out interior vertices with exactly 3 faces (default True; their
///                   energy is 0 and not counted).
/// threads           worker threads (default: all). The result does not depend on it.
///
/// The result has vertex_energy (array('d'), one per vertex of the mesh; a vertex no polygon uses
/// gets 0), counted (array('B')), energy (the total over the counted vertices) and a report with
/// near_ties, n_counted, max_vertex_energy, seconds, seconds_total.
///
/// Raises InputError for malformed arguments and unknown energies, MeshError for a mesh that is
/// not a manifold, consistently oriented triangle mesh, a degenerate triangle, or a vertex where
/// an intrinsic energy is undefined (the message names polygons and vertices of the caller's
/// mesh), TypeError for an argument without the buffer protocol.
#[pyfunction]
#[pyo3(signature = (
    positions,
    corner_vertices,
    polygon_starts,
    triangles,
    *,
    energy = "covariance",
    include_boundary = true,
    omit_valence_three = true,
    threads = None,
))]
#[allow(clippy::too_many_arguments)]
fn developability<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    corner_vertices: &Bound<'py, PyAny>,
    polygon_starts: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    energy: &str,
    include_boundary: bool,
    omit_valence_three: bool,
    threads: Option<i64>,
) -> PyResult<DevelopabilityResult> {
    let start = Instant::now();
    let buffers = mesh(positions, corner_vertices, polygon_starts, triangles, None)?;
    let params = dv::EnergyParams {
        energy: dv::parse_energy("energy", energy).map_err(unwrap_error)?,
        include_boundary,
        omit_valence_three,
        threads: self::threads(threads)?,
    };
    let out = py
        .detach(move || {
            let input = UnwrapInput::from_buffers(buffers)?;
            dv::developability(&input, &params)
        })
        .map_err(unwrap_error)?;
    let report = PyDict::new(py);
    report.set_item("near_ties", out.near_ties)?;
    report.set_item("n_counted", out.counted.iter().filter(|&&c| c).count())?;
    report.set_item("max_vertex_energy", out.vertex_energy.iter().copied().fold(0.0, f64::max))?;
    report.set_item("seconds", out.seconds)?;
    report.set_item("seconds_total", start.elapsed().as_secs_f64())?;
    Ok(DevelopabilityResult {
        vertex_energy: real_array(py, &out.vertex_energy)?.unbind(),
        counted: bool_array(py, &out.counted)?.unbind(),
        energy: out.energy,
        report: report.unbind(),
    })
}

// -------------------------------------------------------------------------------------------
// make_developable
// -------------------------------------------------------------------------------------------

/// The result of `make_developable`.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct MakeDevelopableResult {
    /// `array('d')`, 3 per vertex. With `faces` None: one per vertex of the input mesh. Otherwise
    /// one per vertex of the new mesh.
    positions: Py<PyAny>,
    /// `array('i')`, 3 vertex indices per triangle of the new mesh, or None if the topology did
    /// not change (and `return_faces` was off).
    faces: Option<Py<PyAny>>,
    /// `array('i')`, one per triangle of `faces`: the polygon of the input mesh it descends from
    /// (present exactly when `faces` is).
    face_polygon: Option<Py<PyAny>>,
    /// `array('i')`, one per vertex of the result: the vertex of the input mesh it descends from,
    /// -1 for a vertex created by subdivision.
    vertex_source: Py<PyAny>,
    /// `dict` of plain values; `rounds` is a list of dicts and `warnings` a list of str.
    report: Py<PyDict>,
}

#[pymethods]
impl MakeDevelopableResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let get = |key: &str| -> PyResult<String> {
            Ok(match self.report.bind(py).get_item(key)? {
                Some(value) => value.to_string(),
                None => String::new(),
            })
        };
        Ok(format!(
            "MakeDevelopableResult(vertices={}, topology_changed={}, energy {} -> {})",
            self.vertex_source.bind(py).len()?,
            get("topology_changed")?,
            get("energy_initial")?,
            get("energy_final")?
        ))
    }
}

/// Run the developability flow on a polygon mesh: move the vertices so that the surface becomes
/// developable (flat faces around a vertex fold along straight lines). Returns a
/// MakeDevelopableResult. The input is not modified.
///
/// Method (Stein, Grinspun, Crane 2018, Sec. 4.2 and 5.2): rounds of quasi-Newton minimization of
/// the developability energy on the triangle mesh of the loop triangles, with a 4-1 subdivision
/// between rounds and, during a round, the removal of triangles with a corner angle below
/// small_angle by edge flips and collapses. Every triangle step is capped at max_displacement
/// mean edge lengths.
///
/// Arrays and energy options are as for `developability`.
///
/// Optimization
/// method            "lbfgs" (default), "gradient_descent" or "fixed_step".
/// memory            L-BFGS memory (default 10).
/// step              fixed step, in units of the squared mean edge length (default 1e-3).
/// line_search       "wolfe" (default; weak Wolfe, Lewis-Overton) or "backtracking" (Armijo);
///                   ignored by "fixed_step".
/// first_trial       first trial step of the "wolfe" line search: "unit" (default; 1 at every
///                   iteration) or "grow" (min(first_trial_growth x the previous accepted step,
///                   first_trial_max); SPEC D37). "grow" ends lower on most energies and meshes
///                   (28 of 36 study runs) but can get stuck; with a larger max_displacement it
///                   goes further, by larger moves. Ignored by "backtracking" and "fixed_step".
/// first_trial_growth, first_trial_max
///                   the "grow" rule's factor (>= 1, default 2) and bound (> 0, default 100).
/// iterations        iterations per round (default 500), over restarts after remeshing.
/// subdivisions      number of 4-1 subdivisions; the flow has subdivisions + 1 rounds (default 2).
/// coarse_energy     energy of round 0 when subdivisions > 0 (default "combinatorial"; one of the
///                   energy names, or None for `energy` in every round, the flow of earlier
///                   versions). The paper minimizes E_P on the coarse input and E_lambda after
///                   each subdivision (Sec. 5.2, Figs. 23-24; SPEC D36): E_lambda alone keeps
///                   flattening the shape. Rounds after a subdivision use `energy`. Ignored when
///                   subdivisions is 0. Meant for a coarse, uniform input (about 1k triangles).
/// max_displacement  step cap in mean edge lengths (default 0.5).
/// gradient_tolerance, energy_tolerance, stall_iterations, initial_step
///                   stopping rules and first trial step; None keeps the library's value.
/// min_face_angle    admissibility floor (radians; None keeps the library's value 1e-4, 0
///                   switches the check off): a step that would make a corner angle smaller than
///                   this is refused, and so is a collapse that would. Must be below small_angle.
///                   An input that violates it is a MeshError ("nearly degenerate triangle").
/// pinned_vertices   vertices of the input that do not move (a midpoint created by subdivision
///                   is pinned iff both ends are).
/// Remeshing
/// remesh            small-angle remeshing on (default True).
/// small_angle       a triangle with a corner angle below this (radians, 0 < x < pi/3) is
///                   repaired (default 1 degree).
/// fallback_operations  when the paper's operation for a small-angle triangle is refused, try the
///                   other flips/collapses of that triangle (default True; SPEC D35).
/// collapse_boundary_edges  boundary edges may be collapsed (default True).
/// Output
/// return_faces      return `faces` and `face_polygon` even if the topology did not change (the
///                   triangulation then is the input's, with polygons as given).
/// threads           worker threads (default: all). Results do not depend on it.
/// progress          a callable progress(round, iteration, energy, gradient_norm), called from
///                   the calling thread at most every progress_interval seconds (default 0.1;
///                   0 calls it after every iteration). Returning False stops the flow; the state
///                   at that moment is returned and report["cancelled"] is True. An exception
///                   (and Ctrl-C) stops the flow and propagates.
///
/// faces is None exactly when no subdivision ran and no remeshing operation was applied; then
/// only the positions of the input's vertices change and the caller keeps its polygons. The report
/// then has the largest non-planarity of a polygon with more than three corners before and after
/// (nonplanarity_before, nonplanarity_after, nonplanarity_after_relative, nonplanarity_polygon;
/// None if there is none): the flow makes the triangulation developable, which polygons with
/// more than three corners follow only as far as they are planar.
///
/// The report has energy_initial, energy_final, cancelled, topology_changed, n_input_vertices,
/// n_input_triangles, n_vertices, n_triangles, iterations, max_displacement (last round, in
/// mean edge lengths), mean_edge_length, nonplanarity_* (above), rounds (a list of dicts with
/// energy (name of the energy minimized in that round; energy_initial and energy_final are in it),
/// n_vertices, n_faces, energy_initial, energy_final, iterations, status (gradient_tolerance,
/// energy_stalled, max_iterations, line_search_failed, stationary, cancelled, non_descent),
/// flips, collapses, fallback_flips, fallback_collapses (of the flips and collapses, those applied
/// as a fallback), remesh_rejected, inadmissible_trials, seconds), min_angle (smallest corner angle of the result,
/// radians), bad_faces (triangles of the result below small_angle that nothing could repair),
/// inadmissible_trials (trial steps refused by min_face_angle, all rounds), stopped (None, or a
/// dict round/after_subdivision/vertex/reason: the flow ended before that round because the
/// energy is undefined at its starting state, and the result is the state before it; also named
/// in warnings), warnings, seconds, seconds_total.
///
/// face_polygon is exact: the library tracks the input triangle of every output triangle through
/// subdivision, flips and collapses, and its polygon is reported.
///
/// Raises InputError for malformed arguments or options outside their domain, MeshError for a
/// mesh that is not a manifold, consistently oriented triangle mesh, or has a degenerate
/// triangle (named by polygon), SolveError if the computation fails, KeyboardInterrupt on
/// Ctrl-C, whatever the progress function raises, TypeError for an argument without the buffer
/// protocol.
#[pyfunction]
#[pyo3(signature = (
    positions,
    corner_vertices,
    polygon_starts,
    triangles,
    *,
    energy = "covariance",
    include_boundary = true,
    omit_valence_three = true,
    method = "lbfgs",
    memory = 10,
    step = 1e-3,
    line_search = "wolfe",
    first_trial = "unit",
    first_trial_growth = 2.0,
    first_trial_max = 100.0,
    iterations = 500,
    subdivisions = 2,
    coarse_energy = Some("combinatorial"),
    max_displacement = 0.5,
    gradient_tolerance = None,
    energy_tolerance = None,
    stall_iterations = None,
    initial_step = None,
    min_face_angle = None,
    pinned_vertices = None,
    remesh = true,
    small_angle = std::f64::consts::PI / 180.0,
    collapse_boundary_edges = true,
    fallback_operations = true,
    return_faces = false,
    threads = None,
    progress = None,
    progress_interval = 0.1,
))]
#[allow(clippy::too_many_arguments)]
fn make_developable<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    corner_vertices: &Bound<'py, PyAny>,
    polygon_starts: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    energy: &str,
    include_boundary: bool,
    omit_valence_three: bool,
    method: &str,
    memory: i64,
    step: f64,
    line_search: &str,
    first_trial: &str,
    first_trial_growth: f64,
    first_trial_max: f64,
    iterations: i64,
    subdivisions: i64,
    coarse_energy: Option<&str>,
    max_displacement: f64,
    gradient_tolerance: Option<f64>,
    energy_tolerance: Option<f64>,
    stall_iterations: Option<i64>,
    initial_step: Option<f64>,
    min_face_angle: Option<f64>,
    pinned_vertices: Option<&Bound<'py, PyAny>>,
    remesh: bool,
    small_angle: f64,
    collapse_boundary_edges: bool,
    fallback_operations: bool,
    return_faces: bool,
    threads: Option<i64>,
    progress: Option<Bound<'py, PyAny>>,
    progress_interval: f64,
) -> PyResult<MakeDevelopableResult> {
    let start = Instant::now();
    let buffers = mesh(positions, corner_vertices, polygon_starts, triangles, None)?;
    let reject = |name: &'static str, value: &str, allowed: &str| {
        InputError::new_err(format!("{name}: '{value}' is not valid; use {allowed}"))
    };
    let mut o = dv::MakeOptions::default();
    o.energy = dv::EnergyParams {
        energy: dv::parse_energy("energy", energy).map_err(unwrap_error)?,
        include_boundary,
        omit_valence_three,
        threads: self::threads(threads)?,
    };
    o.method = match method {
        "lbfgs" => Method::Lbfgs { memory: count("memory", memory, 1)? },
        "gradient_descent" => Method::GradientDescent,
        "fixed_step" => Method::FixedStep { step },
        other => return Err(reject("method", other, "lbfgs, gradient_descent or fixed_step")),
    };
    o.line_search = match line_search {
        "wolfe" => o.line_search,
        "backtracking" => LineSearch::Backtracking { c1: 1e-4, shrink: 0.5, max_shrinks: 60 },
        other => return Err(reject("line_search", other, "wolfe or backtracking")),
    };
    o.first_trial = match first_trial {
        "unit" => FirstTrial::Unit,
        "grow" => FirstTrial::GrowPrevious { growth: first_trial_growth, max: first_trial_max },
        other => return Err(reject("first_trial", other, "unit or grow")),
    };
    o.iterations = count("iterations", iterations, 1)?;
    o.subdivisions = count("subdivisions", subdivisions, 0)?;
    o.coarse_energy = coarse_energy
        .map(|name| dv::parse_energy("coarse_energy", name).map_err(unwrap_error))
        .transpose()?;
    o.max_displacement = max_displacement;
    if let Some(x) = gradient_tolerance {
        o.gradient_tolerance = x;
    }
    if let Some(x) = energy_tolerance {
        o.energy_tolerance = x;
    }
    if let Some(x) = stall_iterations {
        o.stall_iterations = count("stall_iterations", x, 1)?;
    }
    if let Some(x) = initial_step {
        o.initial_step = x;
    }
    if let Some(x) = min_face_angle {
        o.min_face_angle = x;
    }
    o.pinned = optional_indices("pinned_vertices", pinned_vertices)?;
    // The remeshing's collapse floor is the same admissibility floor (it must be below small_angle).
    o.remesh = remesh.then(|| RemeshOptions {
        small_angle,
        collapse_boundary_edges,
        fallback_operations,
        min_face_angle: o.min_face_angle,
        ..RemeshOptions::default()
    });
    o.return_faces = return_faces;
    o.progress_interval = progress_interval;
    let progress: Option<Py<PyAny>> = progress.map(Bound::unbind);

    let mut raised: Option<PyErr> = None;
    let out = {
        let raised = &mut raised;
        py.detach(move || {
            let input = UnwrapInput::from_buffers(buffers)?;
            // Takes the GIL back for the callback: Ctrl-C is checked and the user's function
            // called. A False return stops the flow; an exception stops it and is raised after.
            let mut callback = |info: &dv::ProgressInfo| -> bool {
                Python::attach(|py| {
                    if let Err(e) = py.check_signals() {
                        *raised = Some(e);
                        return false;
                    }
                    let Some(f) = &progress else { return true };
                    match f.call1(py, (info.round, info.iteration, info.energy, info.gradient_norm)) {
                        Ok(value) => !matches!(value.extract::<bool>(py), Ok(false)),
                        Err(e) => {
                            *raised = Some(e);
                            false
                        }
                    }
                })
            };
            dv::make_developable(&input, &o, Some(&mut callback))
        })
    };
    if let Some(e) = raised {
        return Err(e);
    }
    let out = out.map_err(unwrap_error)?;

    let r = &out.report;
    let report = PyDict::new(py);
    report.set_item("energy_initial", r.energy_initial)?;
    report.set_item("energy_final", r.energy_final)?;
    report.set_item("cancelled", r.cancelled)?;
    report.set_item("topology_changed", r.topology_changed)?;
    report.set_item("n_input_vertices", r.n_input_vertices)?;
    report.set_item("n_input_triangles", r.n_input_triangles)?;
    report.set_item("n_vertices", r.n_vertices)?;
    report.set_item("n_triangles", r.n_triangles)?;
    report.set_item("iterations", r.iterations)?;
    report.set_item("max_displacement", r.max_displacement)?;
    report.set_item("mean_edge_length", r.mean_edge_length)?;
    report.set_item("min_angle", r.min_angle)?;
    report.set_item("bad_faces", r.bad_faces)?;
    report.set_item("inadmissible_trials", r.inadmissible_trials)?;
    match &r.stopped {
        Some(s) => {
            let d = PyDict::new(py);
            d.set_item("round", s.round)?;
            d.set_item("after_subdivision", s.after_subdivision)?;
            d.set_item("vertex", s.vertex)?;
            d.set_item("reason", s.reason)?;
            report.set_item("stopped", d)?;
        }
        None => report.set_item("stopped", py.None())?,
    }
    report.set_item("nonplanarity_before", r.nonplanarity_before)?;
    report.set_item("nonplanarity_after", r.nonplanarity_after)?;
    report.set_item("nonplanarity_after_relative", r.nonplanarity_after_relative)?;
    report.set_item("nonplanarity_polygon", r.nonplanarity_polygon)?;
    let rounds = PyList::empty(py);
    for round in &r.rounds {
        let d = PyDict::new(py);
        d.set_item("energy", round.energy)?;
        d.set_item("n_vertices", round.n_vertices)?;
        d.set_item("n_faces", round.n_faces)?;
        d.set_item("energy_initial", round.energy_initial)?;
        d.set_item("energy_final", round.energy_final)?;
        d.set_item("iterations", round.iterations)?;
        d.set_item("status", round.status)?;
        d.set_item("flips", round.flips)?;
        d.set_item("collapses", round.collapses)?;
        d.set_item("fallback_flips", round.fallback_flips)?;
        d.set_item("fallback_collapses", round.fallback_collapses)?;
        d.set_item("remesh_rejected", round.remesh_rejected)?;
        d.set_item("inadmissible_trials", round.inadmissible_trials)?;
        d.set_item("seconds", round.seconds)?;
        rounds.append(d)?;
    }
    report.set_item("rounds", rounds)?;
    report.set_item("warnings", PyList::new(py, &r.warnings)?)?;
    report.set_item("seconds", r.seconds)?;
    report.set_item("seconds_total", start.elapsed().as_secs_f64())?;
    let faces = match &out.faces {
        Some(f) => Some(index_array(py, "faces", f.iter().copied())?.unbind()),
        None => None,
    };
    let face_polygon = match &out.face_polygon {
        Some(f) => Some(index_array(py, "face_polygon", f.iter().copied())?.unbind()),
        None => None,
    };
    Ok(MakeDevelopableResult {
        positions: real_array(py, &out.positions)?.unbind(),
        faces,
        face_polygon,
        vertex_source: signed_array(py, "vertex_source", &out.vertex_source)?.unbind(),
        report: report.unbind(),
    })
}


#[allow(clippy::too_many_arguments)]
fn cut_options(
    energy: Option<&str>,
    epsilon: f64,
    defect_tolerance: Option<f64>,
    length_weight: f64,
    alignment_weight: f64,
    extension: &str,
    overlap_samples_per_triangle: Option<f64>,
    threads: Option<i64>,
) -> PyResult<CutOptions> {
    let mut cut = CutOptions::default();
    if let Some(name) = energy {
        cut.energy = dv::parse_energy("energy", name).map_err(unwrap_error)?;
    }
    cut.epsilon = epsilon;
    cut.defect_tolerance = defect_tolerance;
    cut.length_weight = length_weight;
    cut.alignment_weight = alignment_weight;
    cut.extension = match extension {
        "harmonic" => Extension::Harmonic,
        "conjugate" => Extension::Conjugate,
        other => {
            return Err(InputError::new_err(format!(
                "extension: '{other}' is not valid; use harmonic or conjugate"
            )));
        }
    };
    if let Some(x) = overlap_samples_per_triangle {
        cut.overlap_samples_per_triangle = x;
    }
    cut.threads = self::threads(threads)?;
    Ok(cut)
}

// -------------------------------------------------------------------------------------------
// developable_unwrap
// -------------------------------------------------------------------------------------------

/// The result of `developable_unwrap`.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct DevelopableUnwrapResult {
    /// `array('f')`, 2 per polygon corner, packed into the unit square.
    corner_uvs: Py<PyAny>,
    /// `array('i')`, 2 per edge: the edges the surface was cut along, as vertex pairs
    /// `(min, max)`, ascending.
    seam_edges: Py<PyAny>,
    /// `array('i')`, one per polygon: its chart.
    polygon_chart: Py<PyAny>,
    /// `array('d')`, one per polygon corner: the log scale factor of the UV map at the corner
    /// (`0.5 * ln(uv area / area)` of the triangle that owns it; 0 for an isometry).
    corner_log_scale: Py<PyAny>,
    /// `array('i')`: the vertices the cut goes through, ascending.
    terminals: Py<PyAny>,
    /// `array('d')`, one per vertex: the energy that selected the terminals.
    vertex_energy: Py<PyAny>,
    /// `array('d')`, one per vertex: the angle defect (radians; 0 on the boundary).
    vertex_defect: Py<PyAny>,
    /// `array('d')`, one per vertex: the star embeddedness diagnostic (see the README).
    star_embeddedness: Py<PyAny>,
    /// `dict` of plain values; `charts` is a list of dicts and `warnings` a list of str.
    report: Py<PyDict>,
}

#[pymethods]
impl DevelopableUnwrapResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let get = |key: &str| -> PyResult<String> {
            Ok(match self.report.bind(py).get_item(key)? {
                Some(value) => value.to_string(),
                None => String::new(),
            })
        };
        Ok(format!(
            "DevelopableUnwrapResult(corners={}, charts={}, seam_edges={}, flipped_triangles={})",
            self.corner_uvs.bind(py).len()? / 2,
            get("n_charts")?,
            self.seam_edges.bind(py).len()? / 2,
            get("flipped_triangles")?
        ))
    }
}

/// Cut a polygon mesh open through the vertices where it is not developable and flatten the
/// pieces: a UV map. Returns a DevelopableUnwrapResult. The mesh is not changed.
///
/// Method (Stein, Grinspun, Crane 2018, Sec. 4.4): the vertices whose developability energy
/// exceeds epsilon (and, as an addition of this implementation, interior vertices whose angle
/// defect exceeds defect_tolerance) are connected by a cut graph (Steiner tree approximation plus
/// the loops a handle needs), the surface is cut along it, and each chart is flattened with
/// Boundary First Flattening with the lengths along its boundary kept: where the surface is
/// developable and the cut passes through its singular vertices, the charts are isometric. The
/// cut follows polygon edges only (diagonals of the triangulation are never cut), so every
/// polygon lies in one chart. The charts are packed into the unit square by one common scale,
/// without rotation (pack_charts).
///
/// Arrays and energy options are as for `developability`.
///
/// seam_edges        2 vertex indices per edge that must be cut (user seams): edges of polygons.
/// extra_terminals   vertices the cut must go through as well.
/// energy            the energy that selects the terminals; None (default) is the library's
///                   default (currently "covariance"). Evaluated with boundary vertices included
///                   and valence-3 vertices counted.
/// epsilon           energy above which a vertex is a terminal (default 1e-3; the paper
///                   uses 0.01 to 0.001).
/// defect_tolerance  radians; an interior vertex with |angle defect| above this is a terminal as
///                   well (default DEFAULT_DEFECT_TOLERANCE, 0.1). None switches the rule off
///                   (the paper's rule: energy only).
/// length_weight, alignment_weight
///                   the cut edge weight c_len * length + c_align * (1 - alignment of the edge
///                   with the rulings at its ends); defaults 1 and 0 (length only).
/// extension         "harmonic" (default; the boundary curve is reproduced exactly) or
///                   "conjugate" (BFF's automatic mode).
/// overlap_samples_per_triangle  sample density of the overlap estimate (None: library default,
///                   0 switches it off).
/// margin            free space around every chart in the packing, 0 <= margin < 0.5 (default
///                   1/256).
/// threads           worker threads of the energy (default: all).
///
/// The report has energy, n_charts, n_seam_edges, n_terminals, n_added_terminals (terminals the
/// library added on closed pieces with fewer than two), max_vertex_energy, flipped_triangles,
/// degenerate_triangles, is_injective, mean_abs_log_scale (0 for an isometric development),
/// charts (a list of dicts: first_polygon, n_polygons, flipped_triangles, degenerate_triangles,
/// wedges_over_2pi, boundary_crossings, overlapped_fraction, is_injective, quasi_conformal_mean,
/// quasi_conformal_max, log_scale_min, log_scale_max), pack_scale, pack_margin, pack_coverage,
/// seams_on_boundary, warnings, seconds, seconds_total.
///
/// Raises InputError for malformed arguments, options outside their domain, a seam that is not
/// a polygon edge, or a terminal that is not a vertex; MeshError for a mesh that is not a
/// manifold, consistently oriented triangle mesh or has a degenerate triangle; SolveError if the
/// cut or the flattening fails; TypeError for an argument without the buffer protocol.
#[pyfunction]
#[pyo3(signature = (
    positions,
    corner_vertices,
    polygon_starts,
    triangles,
    *,
    seam_edges = None,
    extra_terminals = None,
    energy = None,
    epsilon = 1e-3,
    defect_tolerance = Some(DEFAULT_DEFECT_TOLERANCE),
    length_weight = 1.0,
    alignment_weight = 0.0,
    extension = "harmonic",
    overlap_samples_per_triangle = None,
    margin = dv::DEFAULT_PACK_MARGIN,
    threads = None,
))]
#[allow(clippy::too_many_arguments)]
fn developable_unwrap<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    corner_vertices: &Bound<'py, PyAny>,
    polygon_starts: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    seam_edges: Option<&Bound<'py, PyAny>>,
    extra_terminals: Option<&Bound<'py, PyAny>>,
    energy: Option<&str>,
    epsilon: f64,
    defect_tolerance: Option<f64>,
    length_weight: f64,
    alignment_weight: f64,
    extension: &str,
    overlap_samples_per_triangle: Option<f64>,
    margin: f64,
    threads: Option<i64>,
) -> PyResult<DevelopableUnwrapResult> {
    let start = Instant::now();
    let buffers = mesh(positions, corner_vertices, polygon_starts, triangles, seam_edges)?;
    let extra = optional_indices("extra_terminals", extra_terminals)?;
    let cut = cut_options(
        energy, epsilon, defect_tolerance, length_weight, alignment_weight, extension,
        overlap_samples_per_triangle, threads,
    )?;
    let options = dv::UnwrapDevOptions { cut, margin };
    let out = py
        .detach(move || {
            let input = UnwrapInput::from_buffers(buffers)?;
            dv::developable_unwrap(&input, &extra, &options)
        })
        .map_err(unwrap_error)?;

    let r = &out.report;
    let report = PyDict::new(py);
    report.set_item("energy", r.energy)?;
    report.set_item("n_charts", r.n_charts)?;
    report.set_item("n_seam_edges", r.n_seam_edges)?;
    report.set_item("n_terminals", r.n_terminals)?;
    report.set_item("n_added_terminals", r.n_added_terminals)?;
    report.set_item("max_vertex_energy", r.max_vertex_energy)?;
    report.set_item("flipped_triangles", r.flipped_triangles)?;
    report.set_item("degenerate_triangles", r.degenerate_triangles)?;
    report.set_item("is_injective", r.is_injective)?;
    report.set_item("mean_abs_log_scale", r.mean_abs_log_scale)?;
    let charts = PyList::empty(py);
    for c in &r.charts {
        let d = PyDict::new(py);
        d.set_item("first_polygon", c.first_polygon)?;
        d.set_item("n_polygons", c.n_polygons)?;
        d.set_item("flipped_triangles", c.flipped_triangles)?;
        d.set_item("degenerate_triangles", c.degenerate_triangles)?;
        d.set_item("wedges_over_2pi", c.wedges_over_2pi)?;
        d.set_item("boundary_crossings", c.boundary_crossings)?;
        d.set_item("overlapped_fraction", c.overlapped_fraction)?;
        d.set_item("is_injective", c.is_injective)?;
        d.set_item("quasi_conformal_mean", c.quasi_conformal_mean)?;
        d.set_item("quasi_conformal_max", c.quasi_conformal_max)?;
        d.set_item("log_scale_min", c.log_scale_range[0])?;
        d.set_item("log_scale_max", c.log_scale_range[1])?;
        charts.append(d)?;
    }
    report.set_item("charts", charts)?;
    report.set_item("pack_scale", r.pack_scale)?;
    report.set_item("pack_margin", r.pack_margin)?;
    report.set_item("pack_coverage", r.pack_coverage)?;
    report.set_item("seams_on_boundary", r.seams_on_boundary)?;
    report.set_item("warnings", PyList::new(py, &r.warnings)?)?;
    report.set_item("seconds", r.seconds)?;
    report.set_item("seconds_total", start.elapsed().as_secs_f64())?;
    Ok(DevelopableUnwrapResult {
        corner_uvs: uv_array(py, &out.corner_uvs)?.unbind(),
        seam_edges: index_array(py, "seam_edges", out.seam_edges.as_flattened().iter().copied())?.unbind(),
        polygon_chart: index_array(py, "polygon_chart", out.polygon_chart.iter().copied())?.unbind(),
        corner_log_scale: real_array(py, &out.corner_log_scale)?.unbind(),
        terminals: index_array(py, "terminals", out.terminals.iter().copied())?.unbind(),
        vertex_energy: real_array(py, &out.vertex_energy)?.unbind(),
        vertex_defect: real_array(py, &out.vertex_defect)?.unbind(),
        star_embeddedness: real_array(py, &out.star_embeddedness)?.unbind(),
        report: report.unbind(),
    })
}

/// The result of `developable_seams`.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct DevelopableSeamsResult {
    /// `array('i')`, 2 per edge: the edges to cut, as vertex pairs `(min, max)`, ascending: the
    /// user seams and the cut through the terminals. Polygon edges only.
    seam_edges: Py<PyAny>,
    /// `array('i')`: the vertices the cut goes through, ascending.
    terminals: Py<PyAny>,
    /// `array('d')`, one per vertex: the energy that selected the terminals.
    vertex_energy: Py<PyAny>,
    /// `array('d')`, one per vertex: the angle defect (radians; 0 on the boundary).
    vertex_defect: Py<PyAny>,
    /// `array('d')`, one per vertex: the star embeddedness diagnostic.
    star_embeddedness: Py<PyAny>,
    /// `dict`: energy, n_seam_edges, n_terminals, n_added_terminals, max_vertex_energy,
    /// seams_on_boundary, seconds, seconds_total.
    report: Py<PyDict>,
}

#[pymethods]
impl DevelopableSeamsResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "DevelopableSeamsResult(seam_edges={}, terminals={})",
            self.seam_edges.bind(py).len()? / 2,
            self.terminals.bind(py).len()?
        ))
    }
}

/// The seams of a developable surface: where to cut a polygon mesh so that it can be flattened
/// with little distortion. Returns a DevelopableSeamsResult. The mesh is not changed.
///
/// This is the cut half of `developable_unwrap` (same arguments, same cut rule: terminals from
/// the energy threshold epsilon and the angle-defect rule, a Steiner tree approximation through
/// them, handle loops, polygon edges only), without packing. The UVs are then best made by the
/// general unwrap (`unwrap` with these edges as seam_edges): the plain per-chart flattening of
/// `developable_unwrap` keeps the lengths along the cut but a chart may overlap itself.
///
/// Arguments and errors as for `developable_unwrap` (margin does not exist).
///
/// The result has seam_edges (array('i'), 2 per edge, user seams included), terminals,
/// vertex_energy, vertex_defect, star_embeddedness (array('d'), one per vertex) and a report.
#[pyfunction]
#[pyo3(signature = (
    positions,
    corner_vertices,
    polygon_starts,
    triangles,
    *,
    seam_edges = None,
    extra_terminals = None,
    energy = None,
    epsilon = 1e-3,
    defect_tolerance = Some(DEFAULT_DEFECT_TOLERANCE),
    length_weight = 1.0,
    alignment_weight = 0.0,
    extension = "harmonic",
    overlap_samples_per_triangle = Some(0.0),
    threads = None,
))]
#[allow(clippy::too_many_arguments)]
fn developable_seams<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    corner_vertices: &Bound<'py, PyAny>,
    polygon_starts: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    seam_edges: Option<&Bound<'py, PyAny>>,
    extra_terminals: Option<&Bound<'py, PyAny>>,
    energy: Option<&str>,
    epsilon: f64,
    defect_tolerance: Option<f64>,
    length_weight: f64,
    alignment_weight: f64,
    extension: &str,
    overlap_samples_per_triangle: Option<f64>,
    threads: Option<i64>,
) -> PyResult<DevelopableSeamsResult> {
    let start = Instant::now();
    let buffers = mesh(positions, corner_vertices, polygon_starts, triangles, seam_edges)?;
    let extra = optional_indices("extra_terminals", extra_terminals)?;
    let cut = cut_options(
        energy, epsilon, defect_tolerance, length_weight, alignment_weight, extension,
        overlap_samples_per_triangle, threads,
    )?;
    let options = dv::UnwrapDevOptions { cut, margin: dv::DEFAULT_PACK_MARGIN };
    let out = py
        .detach(move || {
            let input = UnwrapInput::from_buffers(buffers)?;
            dv::developable_unwrap(&input, &extra, &options)
        })
        .map_err(unwrap_error)?;
    let r = &out.report;
    let report = PyDict::new(py);
    report.set_item("energy", r.energy)?;
    report.set_item("n_seam_edges", r.n_seam_edges)?;
    report.set_item("n_terminals", r.n_terminals)?;
    report.set_item("n_added_terminals", r.n_added_terminals)?;
    report.set_item("max_vertex_energy", r.max_vertex_energy)?;
    report.set_item("seams_on_boundary", r.seams_on_boundary)?;
    report.set_item("seconds", r.seconds)?;
    report.set_item("seconds_total", start.elapsed().as_secs_f64())?;
    Ok(DevelopableSeamsResult {
        seam_edges: index_array(py, "seam_edges", out.seam_edges.as_flattened().iter().copied())?.unbind(),
        terminals: index_array(py, "terminals", out.terminals.iter().copied())?.unbind(),
        vertex_energy: real_array(py, &out.vertex_energy)?.unbind(),
        vertex_defect: real_array(py, &out.vertex_defect)?.unbind(),
        star_embeddedness: real_array(py, &out.star_embeddedness)?.unbind(),
        report: report.unbind(),
    })
}
