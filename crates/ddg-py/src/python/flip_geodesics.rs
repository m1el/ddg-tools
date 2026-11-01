//! The PyO3 glue of `geodesic_path`, `geodesic_network` and `geodesic_bezier`
//! ([`crate::flip_geodesics`]). As in the parent module no logic lives here.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use super::{index_array, indices, real_array, reals, unwrap_error};
use crate::flip_geodesics::{
    GeodesicArguments, GeodesicOutput, GeodesicTask, InputCertificateOutput, run,
};
use crate::integer_coords::Places;
use crate::intrinsic_apps::MeshInput;

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(geodesic_path, m)?)?;
    m.add_function(wrap_pyfunction!(geodesic_network, m)?)?;
    m.add_function(wrap_pyfunction!(geodesic_bezier, m)?)?;
    m.add_class::<GeodesicResult>()?;
    Ok(())
}

/// The result of `geodesic_path`, `geodesic_network` and `geodesic_bezier`.
#[pyclass(frozen, module = "ddg_tools")]
pub struct GeodesicResult {
    /// `list` of `array('d')`: per path, 3 reals per point on the input surface (reshape to
    /// `(n, 3)` with numpy). A closed path repeats its first point at the end.
    #[pyo3(get)]
    polylines: Py<PyAny>,
    /// `list` of `dict`: per path, the places of its points on the input mesh: `kind`
    /// (`array('i')`: 0 vertex, 1 point of an edge, 2 point of a face), `id` (vertex, edge or
    /// face), `params` (`array('d')`, 3 per point: `(0, 0, 0)`; `(t, 0, 0)` with `t` running from
    /// the lower-numbered end of the input edge; barycentric coordinates).
    #[pyo3(get)]
    surface_points: Py<PyAny>,
    /// `array('d')`: the length of every path.
    #[pyo3(get)]
    lengths: Py<PyAny>,
    /// `list` of `bool`: whether every path is closed.
    #[pyo3(get)]
    closed: Py<PyAny>,
    /// `dict`: representation, status, iterations, flips, ..., `certificate`,
    /// `input_certificate` and for a Bézier curve `bezier`.
    #[pyo3(get)]
    report: Py<PyDict>,
    /// Bézier only: `array('d')`, 3 reals per control point of the result (the control polygon
    /// of the `2**rounds` pieces); otherwise `None`.
    #[pyo3(get)]
    control_polygon: Option<Py<PyAny>>,
    /// `dict` of the edges of the final triangulation as polylines on the input mesh, as
    /// `IntrinsicTriangulation.edge_polylines` (edges, points, offsets, along_input_edge,
    /// n_crossings, max_relative_mismatch), if `triangulation_edges=True`; otherwise `None`.
    #[pyo3(get)]
    triangulation_edges: Option<Py<PyDict>>,
}

#[pymethods]
impl GeodesicResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let report = self.report.bind(py);
        Ok(format!(
            "GeodesicResult(paths={}, length={:.6}, status={})",
            self.lengths.bind(py).len()?,
            report
                .get_item("length")?
                .map_or(f64::NAN, |x| x.extract().unwrap_or(f64::NAN)),
            report
                .get_item("status")?
                .map_or("?".into(), |x| x.to_string()),
        ))
    }
}

fn places<'py>(py: Python<'py>, p: &Places) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let kinds: Vec<u32> = p.kind.iter().map(|&k| k as u32).collect();
    d.set_item("kind", index_array(py, "kind", kinds.into_iter())?)?;
    d.set_item("id", index_array(py, "id", p.id.iter().copied())?)?;
    d.set_item("params", real_array(py, &p.params)?)?;
    Ok(d)
}

fn input_certificate<'py>(
    py: Python<'py>,
    c: &InputCertificateOutput,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("passed", c.passed)?;
    d.set_item("n_crossings", c.n_crossings)?;
    d.set_item("max_crossing_deviation", c.max_crossing_deviation)?;
    d.set_item("n_joints", c.n_joints)?;
    d.set_item("n_exempt", c.n_exempt)?;
    d.set_item("min_joint_angle", c.min_joint_angle)?;
    d.set_item("max_wedge_mismatch", c.max_wedge_mismatch)?;
    d.set_item("max_length_error", c.max_length_error)?;
    d.set_item("n_faults", c.n_faults)?;
    d.set_item("faults", c.faults.clone())?;
    Ok(d)
}

fn result(py: Python<'_>, out: GeodesicOutput) -> PyResult<GeodesicResult> {
    let polylines = PyList::empty(py);
    for p in &out.polylines {
        polylines.append(real_array(py, p)?)?;
    }
    let surface_points = PyList::empty(py);
    for p in &out.surface_points {
        surface_points.append(places(py, p)?)?;
    }
    let report = PyDict::new(py);
    let s = &out.shortening;
    report.set_item("representation", out.representation)?;
    report.set_item("status", s.status)?;
    report.set_item("length", out.lengths.iter().sum::<f64>())?;
    report.set_item("initial_length", s.initial_length)?;
    report.set_item("iterations", s.iterations)?;
    report.set_item("flips", s.flips)?;
    report.set_item("refused_flips", s.refused_flips)?;
    report.set_item("stalled", s.stalled)?;
    report.set_item("n_paths", out.lengths.len())?;
    report.set_item("n_crossings", out.n_crossings)?;
    report.set_item("max_relative_mismatch", out.max_relative_mismatch)?;
    report.set_item("n_triangulation_vertices", out.n_triangulation_vertices)?;
    let c = &out.certificate;
    let cert = PyDict::new(py);
    cert.set_item("geodesic", c.geodesic)?;
    cert.set_item("tolerance", c.tolerance)?;
    cert.set_item("n_joints", c.n_joints)?;
    cert.set_item("n_pinned", c.n_pinned)?;
    cert.set_item("n_straight", c.n_straight)?;
    cert.set_item("n_blocked", c.n_blocked)?;
    cert.set_item("min_angle", c.min_angle)?;
    cert.set_item("failures", c.failures.clone())?;
    report.set_item("certificate", cert)?;
    report.set_item(
        "input_certificate",
        out.input_certificate
            .as_ref()
            .map(|c| input_certificate(py, c))
            .transpose()?,
    )?;
    let mut control_polygon = None;
    if let Some(b) = &out.bezier {
        let d = PyDict::new(py);
        d.set_item("rounds", b.rounds)?;
        d.set_item(
            "control_points",
            index_array(py, "control_points", b.control_points.iter().copied())?,
        )?;
        d.set_item("n_splits", b.n_splits)?;
        d.set_item("n_snapped", b.n_snapped)?;
        d.set_item("n_delaunay_retries", b.n_delaunay_retries)?;
        d.set_item("n_nudged", b.n_nudged)?;
        d.set_item("n_shortenings", b.n_shortenings)?;
        d.set_item("n_not_geodesic", b.n_not_geodesic)?;
        report.set_item("bezier", d)?;
        control_polygon = Some(real_array(py, &b.control_polygon)?.unbind());
    } else {
        report.set_item("bezier", py.None())?;
    }
    let triangulation_edges = match &out.triangulation_edges {
        None => None,
        Some(o) => {
            let d = PyDict::new(py);
            d.set_item("edges", index_array(py, "edges", o.edges.iter().copied())?)?;
            d.set_item("points", real_array(py, &o.points)?)?;
            d.set_item(
                "offsets",
                index_array(py, "offsets", o.offsets.iter().copied())?,
            )?;
            d.set_item("along_input_edge", o.along_input_edge.clone())?;
            d.set_item("n_crossings", o.n_crossings.clone())?;
            d.set_item("max_relative_mismatch", o.max_relative_mismatch)?;
            Some(d.unbind())
        }
    };
    Ok(GeodesicResult {
        polylines: polylines.into_any().unbind(),
        surface_points: surface_points.into_any().unbind(),
        lengths: real_array(py, &out.lengths)?.unbind(),
        closed: PyList::new(py, &out.closed)?.into_any().unbind(),
        report: report.unbind(),
        control_polygon,
        triangulation_edges,
    })
}

fn mesh_input(positions: &Bound<'_, PyAny>, triangles: &Bound<'_, PyAny>) -> PyResult<MeshInput> {
    Ok(MeshInput {
        positions: reals("positions", positions, 3, "vertex")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
    })
}

#[allow(clippy::too_many_arguments)]
fn arguments(
    representation: Option<&str>,
    mollify: Option<f64>,
    tolerance: f64,
    max_iterations: Option<usize>,
    stop_length_ratio: Option<f64>,
    certify_input: bool,
    triangulation_edges: bool,
) -> GeodesicArguments {
    GeodesicArguments {
        representation: representation.map(str::to_owned),
        mollify,
        tolerance,
        max_iterations,
        stop_length_ratio,
        certify_input,
        triangulation_edges,
    }
}

fn execute(
    py: Python<'_>,
    input: MeshInput,
    task: GeodesicTask,
    args: GeodesicArguments,
) -> PyResult<GeodesicResult> {
    // The extraction runs on rayon's global pool: release the GIL for all of it.
    let out = py
        .detach(move || run(&input, &task, &args))
        .map_err(unwrap_error)?;
    result(py, out)
}

/// A geodesic path or loop through vertices of a mesh, by edge flips (Sharp, Crane 2020, *You
/// Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges*). Returns a
/// GeodesicResult with one path.
///
/// positions       3 reals per vertex.
/// triangles       3 vertex indices per triangle: a manifold, consistently oriented mesh.
/// vertices        the waypoints (at least two; no two consecutive equal): the initial path is the
///                 shortest edge path between consecutive ones.
/// closed          join the last waypoint back to the first (a loop).
/// pin_waypoints   keep the path through every waypoint (default True); False: they only seed the
///                 path, which may leave them (the ends of an open path stay).
/// representation  None (default: integer coordinates, Gillespie-Sharp-Crane 2021), "signposts" or "integer"
///                 (Gillespie-Sharp-Crane 2021).
/// mollify         relative intrinsic mollification (integer only; e.g. 1e-4).
/// tolerance       a joint is straight when its smaller side exceeds pi - tolerance (1e-5).
/// max_iterations  limit on shortened joints (default: 100 (edges + segments)).
/// stop_length_ratio  stop once the length is below this fraction of the initial length.
/// certify_input   also check the path on the input mesh (default True; False is faster).
/// triangulation_edges  also return the edges of the final triangulation as polylines on the
///                 input mesh (default False).
///
/// Raises InputError for malformed arguments (message names the parameter), MeshError for a mesh
/// that is not a manifold triangle mesh, SolveError when the path crosses itself or cannot be
/// traced on the input mesh (the message says what to change).
#[pyfunction]
#[pyo3(signature = (positions, triangles, vertices, *, closed = false, pin_waypoints = true, representation = None, mollify = None, tolerance = 1e-5, max_iterations = None, stop_length_ratio = None, certify_input = true, triangulation_edges = false))]
#[allow(clippy::too_many_arguments)]
fn geodesic_path(
    py: Python<'_>,
    positions: &Bound<'_, PyAny>,
    triangles: &Bound<'_, PyAny>,
    vertices: &Bound<'_, PyAny>,
    closed: bool,
    pin_waypoints: bool,
    representation: Option<&str>,
    mollify: Option<f64>,
    tolerance: f64,
    max_iterations: Option<usize>,
    stop_length_ratio: Option<f64>,
    certify_input: bool,
    triangulation_edges: bool,
) -> PyResult<GeodesicResult> {
    let input = mesh_input(positions, triangles)?;
    let task = GeodesicTask::Path {
        vertices: indices("vertices", vertices, 1, "vertex")?,
        closed,
        pin_waypoints,
    };
    let args = arguments(
        representation,
        mollify,
        tolerance,
        max_iterations,
        stop_length_ratio,
        certify_input,
        triangulation_edges,
    );
    execute(py, input, task, args)
}

/// Geodesic paths and loops along a network of mesh edges (seams, cut graphs), shortened
/// together without crossing, junctions staying put (Sharp, Crane 2020, Sec. 6.1). Returns a
/// GeodesicResult with one path per branch and loop of the network.
///
/// positions, triangles  as for geodesic_path.
/// edges           the edges of the network as vertex pairs, shape (m, 2) or flat; each pair must
///                 be an edge of the mesh. Vertices where other than two edges meet are pinned.
/// pinned          vertices to pin in addition (default none).
///
/// The other options are those of geodesic_path.
#[pyfunction]
#[pyo3(signature = (positions, triangles, edges, *, pinned = None, representation = None, mollify = None, tolerance = 1e-5, max_iterations = None, stop_length_ratio = None, certify_input = true, triangulation_edges = false))]
#[allow(clippy::too_many_arguments)]
fn geodesic_network(
    py: Python<'_>,
    positions: &Bound<'_, PyAny>,
    triangles: &Bound<'_, PyAny>,
    edges: &Bound<'_, PyAny>,
    pinned: Option<&Bound<'_, PyAny>>,
    representation: Option<&str>,
    mollify: Option<f64>,
    tolerance: f64,
    max_iterations: Option<usize>,
    stop_length_ratio: Option<f64>,
    certify_input: bool,
    triangulation_edges: bool,
) -> PyResult<GeodesicResult> {
    let input = mesh_input(positions, triangles)?;
    let task = GeodesicTask::Network {
        edges: indices("edges", edges, 2, "edge")?,
        pinned: pinned
            .map(|p| indices("pinned", p, 1, "vertex"))
            .transpose()?
            .unwrap_or_default(),
    };
    let args = arguments(
        representation,
        mollify,
        tolerance,
        max_iterations,
        stop_length_ratio,
        certify_input,
        triangulation_edges,
    );
    execute(py, input, task, args)
}

/// A geodesic Bézier curve with control vertices on the surface, by midpoint subdivision (Sharp,
/// Crane 2020, Sec. 6.2; Morera et al. 2008). Returns a GeodesicResult with one path and the
/// control polygon of the result.
///
/// positions, triangles  as for geodesic_path.
/// control_vertices  the control points (at least two; no two consecutive equal; none repeated).
/// rounds          subdivision rounds, 0 to 24 (default 4): the curve is made of 2**rounds pieces;
///                 1 to 8 suit interactive use.
///
/// The other options are those of geodesic_path. The report has a `bezier` dict with the
/// control points (vertex ids of the final triangulation) and the counts of the subdivision.
///
/// Raises SolveError when the shortest paths between the control points touch or cross ("move a
/// control point") or when a control point is repeated.
#[pyfunction]
#[pyo3(signature = (positions, triangles, control_vertices, *, rounds = 4, representation = None, mollify = None, tolerance = 1e-5, max_iterations = None, stop_length_ratio = None, certify_input = true, triangulation_edges = false))]
#[allow(clippy::too_many_arguments)]
fn geodesic_bezier(
    py: Python<'_>,
    positions: &Bound<'_, PyAny>,
    triangles: &Bound<'_, PyAny>,
    control_vertices: &Bound<'_, PyAny>,
    rounds: usize,
    representation: Option<&str>,
    mollify: Option<f64>,
    tolerance: f64,
    max_iterations: Option<usize>,
    stop_length_ratio: Option<f64>,
    certify_input: bool,
    triangulation_edges: bool,
) -> PyResult<GeodesicResult> {
    let input = mesh_input(positions, triangles)?;
    let task = GeodesicTask::Bezier {
        control_vertices: indices("control_vertices", control_vertices, 1, "vertex")?,
        rounds,
    };
    let args = arguments(
        representation,
        mollify,
        tolerance,
        max_iterations,
        stop_length_ratio,
        certify_input,
        triangulation_edges,
    );
    execute(py, input, task, args)
}
