//! The PyO3 glue of the Python module `ddg_tools`.
//!
//! No logic lives here. Every array argument is taken through the buffer protocol, described to
//! [`crate::buffers`] as a [`BufferView`] and copied into Rust vectors while the GIL is held;
//! the work then runs with the GIL released; results go back as `array.array` objects and a
//! plain `dict`. The module is built for the stable ABI (`abi3`, Python 3.11 or newer), so one
//! binary serves every Blender, whichever Python it embeds, and numpy is never required.

use std::time::Instant;

use pyo3::buffer::PyUntypedBuffer;
use pyo3::exceptions::{PyException, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList};

use crate::buffers::{BufferView, decode_indices, decode_reals};
use crate::cones::{ConesInput, DEFAULT_LAMBDA, reoptimize};
use crate::contract::{
    MeshBuffers, ReportValue, UnwrapInput, UnwrapOptions, polygon_offsets_from_starts,
};
use crate::error::{InputError as RustInputError, UnwrapError};
use crate::pipeline::{IS_PLACEHOLDER, PIPELINE_NAME};

mod cut;
mod developability;
mod sphere;
mod intrinsic_apps;
mod intrinsic_vector;

pyo3::create_exception!(
    ddg_tools,
    DdgError,
    PyException,
    "Base class of the errors raised by ddg_tools."
);
pyo3::create_exception!(
    ddg_tools,
    InputError,
    DdgError,
    "An argument is malformed: wrong buffer element type, wrong length, an index out of range, \
     contradictory steering, or an option outside its domain. The message starts with the name \
     of the parameter."
);
pyo3::create_exception!(
    ddg_tools,
    MeshError,
    DdgError,
    "The arrays are well formed but the mesh cannot be unwrapped as it is (non-manifold, \
     inconsistently oriented, ...). The message names the offending vertices and faces by their \
     indices in the caller's mesh."
);
pyo3::create_exception!(
    ddg_tools,
    SolveError,
    DdgError,
    "The mesh was accepted but no valid UV map could be produced."
);

fn input_error(e: RustInputError) -> PyErr {
    InputError::new_err(e.to_string())
}

fn unwrap_error(e: UnwrapError) -> PyErr {
    match e {
        UnwrapError::Input(e) => input_error(e),
        UnwrapError::Mesh(message) => MeshError::new_err(message),
        UnwrapError::Solve(message) => SolveError::new_err(message),
    }
}

/// Takes the buffer of `obj` and hands a view of it to `decode`.
///
/// A missing buffer protocol is a `TypeError`; everything `decode` rejects is an `InputError`.
fn with_view<T>(
    name: &'static str,
    obj: &Bound<'_, PyAny>,
    decode: impl FnOnce(&BufferView<'_>) -> Result<T, RustInputError>,
) -> PyResult<T> {
    let buffer = PyUntypedBuffer::get(obj).map_err(|_| {
        let type_name = obj
            .get_type()
            .name()
            .map(|n| n.to_string())
            .unwrap_or_else(|_| "?".to_owned());
        PyTypeError::new_err(format!(
            "{name}: expected an object that supports the buffer protocol (array.array, bytes, \
             bytearray, memoryview, numpy.ndarray), got {type_name}"
        ))
    })?;
    if !buffer.is_c_contiguous() {
        return Err(input_error(RustInputError::NotContiguous { name }));
    }
    let format = buffer.format().to_string_lossy();
    let len = buffer.len_bytes();
    let bytes: &[u8] = if len == 0 {
        &[]
    } else {
        // SAFETY: the exporter keeps `len` bytes at `buf_ptr()` alive and in place until the
        // buffer is released, which happens when `buffer` is dropped at the end of this
        // function; C-contiguity (checked above) means those bytes are the items in C order.
        // The GIL is held and no Python code runs while the slice is alive, so nothing mutates
        // the memory through the interpreter in the meantime.
        unsafe { std::slice::from_raw_parts(buffer.buf_ptr().cast::<u8>(), len) }
    };
    let view = BufferView {
        name,
        format: &format,
        item_size: buffer.item_size(),
        shape: buffer.shape(),
        bytes,
    };
    decode(&view).map_err(input_error)
}

fn reals(
    name: &'static str,
    obj: &Bound<'_, PyAny>,
    width: usize,
    row: &'static str,
) -> PyResult<Vec<f64>> {
    with_view(name, obj, |view| decode_reals(view, width, row))
}

fn indices(
    name: &'static str,
    obj: &Bound<'_, PyAny>,
    width: usize,
    row: &'static str,
) -> PyResult<Vec<u32>> {
    with_view(name, obj, |view| decode_indices(view, width, row))
}

/// `array.array(typecode, bytes)`.
fn new_array<'py>(py: Python<'py>, typecode: &str, bytes: &[u8]) -> PyResult<Bound<'py, PyAny>> {
    py.import("array")?
        .getattr("array")?
        .call1((typecode, PyBytes::new(py, bytes)))
}

/// `array.array('i')` from indices; an index beyond the int32 range is an error.
fn index_array<'py>(
    py: Python<'py>,
    name: &'static str,
    values: impl ExactSizeIterator<Item = u32>,
) -> PyResult<Bound<'py, PyAny>> {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for (position, value) in values.enumerate() {
        let value = i32::try_from(value).map_err(|_| {
            SolveError::new_err(format!(
                "internal error: {name}[{position}] = {value} does not fit the int32 output array"
            ))
        })?;
        bytes.extend_from_slice(&value.to_ne_bytes());
    }
    new_array(py, "i", &bytes)
}

/// `array.array('f')` from UVs, rounded to float32; a value that overflows float32 is an error.
fn uv_array<'py>(py: Python<'py>, uvs: &[[f64; 2]]) -> PyResult<Bound<'py, PyAny>> {
    let mut bytes = Vec::with_capacity(uvs.len() * 8);
    for (corner, uv) in uvs.iter().enumerate() {
        for &x in uv {
            let rounded = x as f32;
            if !rounded.is_finite() {
                return Err(SolveError::new_err(format!(
                    "the UV coordinate {x} at polygon corner {corner} does not fit float32"
                )));
            }
            bytes.extend_from_slice(&rounded.to_ne_bytes());
        }
    }
    new_array(py, "f", &bytes)
}

/// The result of `unwrap`. All arrays are `array.array` objects.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct UnwrapResult {
    /// `array('f')`, 2 per polygon corner, in corner ("loop") order: `uvs[2*c]`, `uvs[2*c+1]`.
    /// Not packed: charts overlap until `pack_charts` lays them out.
    uvs: Py<PyAny>,
    /// `array('i')`, one per polygon: the chart the polygon belongs to.
    polygon_charts: Py<PyAny>,
    /// `array('i')`: the vertices that carry a cone, ascending.
    cone_vertices: Py<PyAny>,
    /// `array('d')`: the cone angles in radians, parallel to `cone_vertices`.
    cone_angles: Py<PyAny>,
    /// `array('i')`, 2 per edge: the edges the surface was cut along, as vertex pairs.
    cut_edges: Py<PyAny>,
    /// `dict` of plain values (`bool`, `int`, `float`, `str`, and `list[str]` for `warnings`).
    report: Py<PyDict>,
}

#[pymethods]
impl UnwrapResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let len = |obj: &Py<PyAny>| obj.bind(py).len();
        Ok(format!(
            "UnwrapResult(corners={}, polygons={}, cones={}, cut_edges={}, pipeline={:?})",
            len(&self.uvs)? / 2,
            len(&self.polygon_charts)?,
            len(&self.cone_vertices)?,
            len(&self.cut_edges)? / 2,
            PIPELINE_NAME
        ))
    }
}

/// Compute a UV map for a polygon mesh. Returns an UnwrapResult.
///
/// Every array argument is any object with the buffer protocol (array.array, bytes, bytearray,
/// memoryview, numpy.ndarray), C-contiguous, flat or of shape (n, width). Real numbers are
/// float32 or float64; indices are 16-, 32- or 64-bit integers; a byte buffer is read as raw
/// native-endian float32 (real numbers) or int32 (indices).
///
/// positions         3 reals per vertex.
/// corner_vertices   1 index per polygon corner ("loop"): its vertex.
/// polygon_starts    1 index per polygon: its first corner (Blender's loop_start). Polygons
///                   own consecutive corners, the first starts at 0, each has at least 3.
/// triangles         3 corner indices per triangle: Blender's mesh.loop_triangles[...].loops.
/// seam_edges        2 vertex indices per edge to cut along before anything else.
/// pinned_vertices   vertices where a cone costs nothing.
/// prescribed_vertices, prescribed_angles
///                   vertices with a fixed cone angle, and those angles in radians (reals).
/// forbidden_vertices  vertices where no cone may be placed.
/// lambda_           trade-off between distortion and cones, >= 0; larger means fewer cones.
/// positive_cones_only  allow only cones with a positive angle.
/// optimize_boundary  let the solver choose the scale factor along the boundary of every piece
///                   (the mesh's boundary and the seams) instead of keeping boundary lengths:
///                   lower distortion, but the two sides of a seam may get different lengths
///                   (report["max_seam_length_mismatch"]), and the flattening can fall short of
///                   the optimized values on some meshes (report["max_boundary_length_error"]).
///                   Default False.
/// split_overlaps    split charts that overlap themselves: a cut that separates the two layers
///                   is added to the seams and the charts are flattened again, a cut is kept
///                   only if it lowers the overlap without mirroring a triangle
///                   (report["n_splits"], ["overlapped_fraction_before_split"]). Default False.
/// split_max         at most this many cuts are kept (default 32).
/// split_threshold   a chart is split while more than this fraction of its footprint is covered
///                   more than once (default 0.001).
/// cone_feedback     compare the cone solver's linearized scale factor with the one the
///                   flattening realizes, lower the distortion weight where they differ by more
///                   than feedback_threshold (thin spikes) and solve again, until nothing changes
///                   (report["cone_feedback_rounds"], ["n_cone_feedback_reweighted"],
///                   ["cone_feedback_vertices"], ["cone_feedback_max_mismatch"]). Default False.
/// feedback_threshold  the mismatch threshold in log units (default 20.0).
/// feedback_max_rounds the most solves, the first included (default 4).
///
/// Raises InputError for malformed arguments, MeshError for a mesh that cannot be unwrapped
/// (for example "edge (12, 57) is shared by 3 faces: 4, 9, 31"), SolveError when no valid map
/// could be produced, and TypeError for an argument without the buffer protocol.
#[pyfunction]
#[pyo3(signature = (
    positions,
    corner_vertices,
    polygon_starts,
    triangles,
    *,
    seam_edges = None,
    pinned_vertices = None,
    prescribed_vertices = None,
    prescribed_angles = None,
    forbidden_vertices = None,
    lambda_ = 0.1,
    positive_cones_only = true,
    optimize_boundary = false,
    split_overlaps = false,
    split_max = 32,
    split_threshold = 0.001,
    cone_feedback = false,
    feedback_threshold = 20.0,
    feedback_max_rounds = 4,
))]
#[allow(clippy::too_many_arguments)]
fn unwrap<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    corner_vertices: &Bound<'py, PyAny>,
    polygon_starts: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    seam_edges: Option<&Bound<'py, PyAny>>,
    pinned_vertices: Option<&Bound<'py, PyAny>>,
    prescribed_vertices: Option<&Bound<'py, PyAny>>,
    prescribed_angles: Option<&Bound<'py, PyAny>>,
    forbidden_vertices: Option<&Bound<'py, PyAny>>,
    lambda_: f64,
    positive_cones_only: bool,
    optimize_boundary: bool,
    split_overlaps: bool,
    split_max: usize,
    split_threshold: f64,
    cone_feedback: bool,
    feedback_threshold: f64,
    feedback_max_rounds: usize,
) -> PyResult<UnwrapResult> {
    let start = Instant::now();
    let optional_indices = |name, obj: Option<&Bound<'py, PyAny>>, width, row| match obj {
        Some(obj) => indices(name, obj, width, row),
        None => Ok(Vec::new()),
    };
    let buffers = MeshBuffers {
        positions: reals("positions", positions, 3, "vertex")?,
        corner_vertices: indices("corner_vertices", corner_vertices, 1, "corner")?,
        polygon_starts: indices("polygon_starts", polygon_starts, 1, "polygon")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
        seam_edges: optional_indices("seam_edges", seam_edges, 2, "edge")?,
        pinned_vertices: optional_indices("pinned_vertices", pinned_vertices, 1, "vertex")?,
        prescribed_vertices: optional_indices(
            "prescribed_vertices",
            prescribed_vertices,
            1,
            "vertex",
        )?,
        prescribed_angles: match prescribed_angles {
            Some(obj) => reals("prescribed_angles", obj, 1, "cone")?,
            None => Vec::new(),
        },
        forbidden_vertices: optional_indices(
            "forbidden_vertices",
            forbidden_vertices,
            1,
            "vertex",
        )?,
    };
    let options = UnwrapOptions {
        lambda: lambda_,
        positive_cones_only,
        optimize_boundary,
        split_overlaps,
        split_max,
        split_threshold,
        cone_feedback,
        feedback_threshold,
        feedback_max_rounds,
    };
    let seconds_decode = start.elapsed().as_secs_f64();

    // The arguments are copied; from here on the interpreter is free to run other threads.
    let (output, seconds_validate) = py
        .detach(move || {
            let validate = Instant::now();
            let input = UnwrapInput::from_buffers(buffers)?;
            let seconds_validate = validate.elapsed().as_secs_f64();
            crate::pipeline::unwrap(&input, &options).map(|output| (output, seconds_validate))
        })
        .map_err(unwrap_error)?;

    let encode = Instant::now();
    let uvs = uv_array(py, &output.uvs)?;
    let polygon_charts = index_array(py, "polygon_charts", output.polygon_charts.iter().copied())?;
    let cone_vertices = index_array(py, "cone_vertices", output.cones.iter().map(|c| c.vertex))?;
    let cone_angles = {
        let bytes: Vec<u8> = output
            .cones
            .iter()
            .flat_map(|c| c.angle.to_ne_bytes())
            .collect();
        new_array(py, "d", &bytes)?
    };
    let cut_edges = index_array(
        py,
        "cut_edges",
        output.cut_edges.as_flattened().iter().copied(),
    )?;

    let report = PyDict::new(py);
    for (key, value) in output.report.entries() {
        match value {
            ReportValue::Bool(v) => report.set_item(key, *v)?,
            ReportValue::Int(v) => report.set_item(key, *v)?,
            ReportValue::Float(v) => report.set_item(key, *v)?,
            ReportValue::Text(v) => report.set_item(key, v)?,
            ReportValue::TextList(v) => report.set_item(key, PyList::new(py, v)?)?,
        }
    }
    report.set_item("seconds_decode", seconds_decode)?;
    report.set_item("seconds_validate", seconds_validate)?;
    report.set_item("seconds_encode", encode.elapsed().as_secs_f64())?;
    report.set_item("seconds_total", start.elapsed().as_secs_f64())?;
    report.set_item("warnings", PyList::new(py, output.report.warnings())?)?;

    Ok(UnwrapResult {
        uvs: uvs.unbind(),
        polygon_charts: polygon_charts.unbind(),
        cone_vertices: cone_vertices.unbind(),
        cone_angles: cone_angles.unbind(),
        cut_edges: cut_edges.unbind(),
        report: report.unbind(),
    })
}

/// Lay the charts of a UV map out in the unit square. Returns (uvs, info): a new array('f') with
/// 2 values per polygon corner, and a dict with n_charts, scale, margin and coverage.
///
/// uvs             2 reals per polygon corner (UnwrapResult.uvs).
/// polygon_starts  1 index per polygon: its first corner, as passed to unwrap.
/// polygon_charts  1 index per polygon: its chart (UnwrapResult.polygon_charts).
/// margin          free space around every chart's bounding box, as a fraction of the unit
///                 square, 0 <= margin < 0.5.
///
/// Charts are translated and scaled by one common factor, never rotated. Their bounding boxes,
/// padded by the margin, are placed on shelves, tallest first, from the lower left corner; the
/// scale is the largest for which the shelves fit. With many charts the margin is capped at
/// 1 / (4 * ceil(sqrt(n_charts))); the dict says which margin was used.
///
/// Raises InputError for malformed arguments.
#[pyfunction]
#[pyo3(signature = (uvs, polygon_starts, polygon_charts, *, margin = 0.0))]
fn pack_charts<'py>(
    py: Python<'py>,
    uvs: &Bound<'py, PyAny>,
    polygon_starts: &Bound<'py, PyAny>,
    polygon_charts: &Bound<'py, PyAny>,
    margin: f64,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyDict>)> {
    let flat = reals("uvs", uvs, 2, "corner")?;
    let starts = indices("polygon_starts", polygon_starts, 1, "polygon")?;
    let charts = indices("polygon_charts", polygon_charts, 1, "polygon")?;
    let (packed, report) = py
        .detach(move || {
            let mut uvs: Vec<[f64; 2]> = flat.chunks_exact(2).map(|c| [c[0], c[1]]).collect();
            let offsets = polygon_offsets_from_starts(starts, uvs.len())?;
            let report = crate::pack::pack_charts(&mut uvs, &offsets, &charts, margin)?;
            Ok((uvs, report))
        })
        .map_err(input_error)?;
    let info = PyDict::new(py);
    info.set_item("n_charts", report.n_charts)?;
    info.set_item("scale", report.scale)?;
    info.set_item("margin", report.margin)?;
    info.set_item("coverage", report.coverage)?;
    Ok((uv_array(py, &packed)?, info))
}

/// The result of `reoptimize_cones`. All arrays are `array.array` objects.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct ReoptimizeResult {
    /// `array('i')`: the cone vertices, as given.
    cone_vertices: Py<PyAny>,
    /// `array('d')`: the optimal cone angles in radians, parallel to `cone_vertices`.
    cone_angles: Py<PyAny>,
    /// `array('d')`: the log scale factor `u`, one per vertex.
    log_scale: Py<PyAny>,
    /// `array('d')`: the boundary curvature, one per vertex.
    boundary_curvature: Py<PyAny>,
    /// `dict` of plain values (`bool`, `int`, `float`, `str`, `None`).
    report: Py<PyDict>,
}

#[pymethods]
impl ReoptimizeResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let status = match self.report.bind(py).get_item("status")? {
            Some(status) => status.to_string(),
            None => String::new(),
        };
        Ok(format!(
            "ReoptimizeResult(cones={}, vertices={}, status={:?})",
            self.cone_vertices.bind(py).len()?,
            self.log_scale.bind(py).len()?,
            status
        ))
    }
}

/// `array.array('d')` from reals.
fn real_array<'py>(py: Python<'py>, values: &[f64]) -> PyResult<Bound<'py, PyAny>> {
    let bytes: Vec<u8> = values.iter().flat_map(|x| x.to_ne_bytes()).collect();
    new_array(py, "d", &bytes)
}

/// Find the least-distortion cone angles for given cone positions. Returns a ReoptimizeResult.
///
/// Array arguments are read as for `unwrap`.
///
/// positions         3 reals per vertex.
/// triangles         3 vertex indices per triangle: a manifold, consistently oriented triangle
///                   mesh in which every vertex is used. Several connected components are allowed.
/// cone_vertices     the vertices where cones are placed; their angles are optimized. Every other
///                   interior vertex gets no cone. Interior vertices only, each listed once.
/// boundary          "fixed" (log scale u = 0 on the boundary, the default), "free", "convex" or
///                   "penalized" (a cost on the boundary curvature, weight lambda_).
/// distortion_weights  one weight >= 0 per vertex (default 1).
/// lower_bounds, upper_bounds
///                   one bound per vertex on the cone angles, <= 0 and >= 0 (default none).
/// prescribed_vertices, prescribed_angles
///                   vertices with a fixed cone angle, and those angles in radians.
/// lambda_           cost of the boundary curvature with boundary="penalized" (finite, >= 0,
///                   default 0.1); ignored otherwise.
/// method            "auto" (default), "direct" (the linear optimality system wherever it
///                   applies) or "solve" (always the generic solver).
///
/// The result has cone_vertices (as given), cone_angles (radians, parallel, may be 0),
/// log_scale and boundary_curvature (one per vertex) and a report dict: status ("converged",
/// "not_converged", "linear_solver_failed"), method, distortion, objective, total_cone_angle,
/// n_components, certificate_holds, state_residual, duality_gap, seconds, message, seconds_total.
///
/// Raises InputError for malformed arguments, a cone vertex that is out of range, repeated, on
/// the boundary or prescribed, and a closed component whose cones cannot carry 2*pi*chi;
/// MeshError for a mesh that is not a valid manifold triangle mesh; TypeError for an argument
/// without the buffer protocol. A result that did not converge is not an error: see
/// report["status"].
#[pyfunction]
#[pyo3(signature = (
    positions,
    triangles,
    cone_vertices,
    *,
    boundary = "fixed",
    distortion_weights = None,
    lower_bounds = None,
    upper_bounds = None,
    prescribed_vertices = None,
    prescribed_angles = None,
    method = "auto",
    lambda_ = DEFAULT_LAMBDA,
))]
#[allow(clippy::too_many_arguments)]
fn reoptimize_cones<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    cone_vertices: &Bound<'py, PyAny>,
    boundary: &str,
    distortion_weights: Option<&Bound<'py, PyAny>>,
    lower_bounds: Option<&Bound<'py, PyAny>>,
    upper_bounds: Option<&Bound<'py, PyAny>>,
    prescribed_vertices: Option<&Bound<'py, PyAny>>,
    prescribed_angles: Option<&Bound<'py, PyAny>>,
    method: &str,
    lambda_: f64,
) -> PyResult<ReoptimizeResult> {
    let start = Instant::now();
    let optional_reals = |name, obj: Option<&Bound<'py, PyAny>>| match obj {
        Some(obj) => reals(name, obj, 1, "vertex").map(Some),
        None => Ok(None),
    };
    let input = ConesInput {
        positions: reals("positions", positions, 3, "vertex")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
        cone_vertices: indices("cone_vertices", cone_vertices, 1, "cone")?,
        boundary: boundary.to_owned(),
        method: method.to_owned(),
        lambda: lambda_,
        distortion_weights: optional_reals("distortion_weights", distortion_weights)?,
        lower_bounds: optional_reals("lower_bounds", lower_bounds)?,
        upper_bounds: optional_reals("upper_bounds", upper_bounds)?,
        prescribed_vertices: match prescribed_vertices {
            Some(obj) => indices("prescribed_vertices", obj, 1, "cone")?,
            None => Vec::new(),
        },
        prescribed_angles: match prescribed_angles {
            Some(obj) => reals("prescribed_angles", obj, 1, "cone")?,
            None => Vec::new(),
        },
    };
    let output = py
        .detach(move || reoptimize(&input))
        .map_err(unwrap_error)?;

    let report = PyDict::new(py);
    let r = &output.report;
    report.set_item("status", r.status)?;
    report.set_item("method", r.method)?;
    report.set_item("distortion", r.distortion)?;
    report.set_item("objective", r.objective)?;
    report.set_item("total_cone_angle", r.total_cone_angle)?;
    report.set_item("n_components", r.n_components)?;
    report.set_item("certificate_holds", r.certificate_holds)?;
    report.set_item("state_residual", r.state_residual)?;
    report.set_item("duality_gap", r.duality_gap)?;
    report.set_item("seconds", r.seconds)?;
    report.set_item("message", r.message.as_deref())?;
    report.set_item("seconds_total", start.elapsed().as_secs_f64())?;

    Ok(ReoptimizeResult {
        cone_vertices: index_array(py, "cone_vertices", output.cone_vertices.iter().copied())?
            .unbind(),
        cone_angles: real_array(py, &output.cone_angles)?.unbind(),
        log_scale: real_array(py, &output.log_scale)?.unbind(),
        boundary_curvature: real_array(py, &output.boundary_curvature)?.unbind(),
        report: report.unbind(),
    })
}

/// Mesh buffers in, per-corner UVs out: the Python side of the ddg-tools Rust workspace.
///
/// unwrap(...)       compute a UV map (see its docstring for the array formats)
/// pack_charts(...)  lay the charts of a UV map out in the unit square
/// reoptimize_cones(...)  the least-distortion angles of cones at given vertices
/// map_to_sphere(...)  a conformal map of a closed genus-0 mesh to the unit sphere
/// developability(...), make_developable(...), developable_seams(...), developable_unwrap(...)
///                   the developability energy, flow and cut of Stein, Grinspun, Crane (2018)
///
/// PIPELINE          name of the unwrap pipeline this build contains
/// IS_PLACEHOLDER    False: the linked pipeline is the real unwrap (`ddg-uv`), not a stand-in
///                   (the flag exists so that a build with a stand-in would label itself)
///
/// CutSession(...)   the cut flow on a polygon mesh (see its docstring)
///
/// Errors: DdgError, with subclasses InputError, MeshError, SolveError and NotAvailableError.
#[pymodule(name = "ddg_tools")]
fn ddg_tools_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add_function(wrap_pyfunction!(unwrap, m)?)?;
    m.add_function(wrap_pyfunction!(pack_charts, m)?)?;
    m.add_function(wrap_pyfunction!(reoptimize_cones, m)?)?;
    m.add_class::<ReoptimizeResult>()?;
    m.add_class::<UnwrapResult>()?;
    cut::register(m)?;
    sphere::register(m)?;
    developability::register(m)?;
    intrinsic_apps::register(m)?;
    intrinsic_vector::register(m)?;
    m.add("DdgError", py.get_type::<DdgError>())?;
    m.add("InputError", py.get_type::<InputError>())?;
    m.add("MeshError", py.get_type::<MeshError>())?;
    m.add("SolveError", py.get_type::<SolveError>())?;
    m.add("PIPELINE", PIPELINE_NAME)?;
    m.add("IS_PLACEHOLDER", IS_PLACEHOLDER)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
