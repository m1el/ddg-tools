//! The PyO3 glue of the cut flow: `ddg_tools.CutSession` and its results.
//!
//! As in the parent module no logic lives here: buffers are decoded, [`crate::cutting`] does the
//! work with the GIL released, results go back as `array.array` objects, plain dictionaries and
//! small result classes.

use std::sync::{Mutex, MutexGuard};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList, PyTuple};

use super::{
    DdgError, InputError, MeshError, SolveError, index_array, indices, new_array, reals, uv_array,
    with_view,
};
use crate::buffers::decode_reals;
use crate::contract::{MeshBuffers, UnwrapInput};
use crate::cutting::{
    BindError, Controls, CutSession as RustCutSession, DescriptionArrays, Dict, EnergyKind,
    FluxChoice, InsertedArrays, StageSet, StartSpec, StepRule, TriangulationSpec, UnwrapArrays,
    UnwrapSettings, Value,
};

pyo3::create_exception!(
    ddg_tools,
    StaleStateError,
    DdgError,
    "A stored state belongs to another triangulation: the mesh or its seams changed since it \
     was made. Start again from corner values, labels or a new start."
);

pyo3::create_exception!(
    ddg_tools,
    NotAvailableError,
    DdgError,
    "The request is valid but not implemented yet: for example a cut session with more than two \
     regions. The message says what is missing."
);

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add_class::<CutSession>()?;
    m.add_class::<CutDescription>()?;
    m.add_class::<InsertedMesh>()?;
    m.add_class::<CutUnwrapResult>()?;
    m.add_function(wrap_pyfunction!(visibility_penalty, m)?)?;
    m.add_function(wrap_pyfunction!(length_weight_for_ripples, m)?)?;
    m.add("NotAvailableError", py.get_type::<NotAvailableError>())?;
    m.add("StaleStateError", py.get_type::<StaleStateError>())?;
    Ok(())
}

fn bind_error(e: BindError) -> PyErr {
    match e {
        BindError::Input(m) => InputError::new_err(m),
        BindError::Mesh(m) => MeshError::new_err(m),
        BindError::NotAvailable(m) => NotAvailableError::new_err(m),
        BindError::Stale(m) => StaleStateError::new_err(m),
        BindError::Solve(m) => SolveError::new_err(m),
    }
}

fn value_to_py<'py>(py: Python<'py>, value: &Value) -> PyResult<Bound<'py, PyAny>> {
    Ok(match value {
        Value::Bool(v) => v.into_pyobject(py)?.to_owned().into_any(),
        Value::Int(v) => v.into_pyobject(py)?.into_any(),
        Value::Float(v) => v.into_pyobject(py)?.into_any(),
        Value::Text(v) => v.into_pyobject(py)?.into_any(),
        Value::TextList(v) => PyList::new(py, v)?.into_any(),
        Value::List(v) => {
            let items = v
                .iter()
                .map(|x| value_to_py(py, x))
                .collect::<PyResult<Vec<_>>>()?;
            PyList::new(py, items)?.into_any()
        }
        Value::Dict(d) => dict_to_py(py, d)?.into_any(),
    })
}

fn dict_to_py<'py>(py: Python<'py>, dict: &Dict) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new(py);
    for (key, value) in &dict.0 {
        out.set_item(key, value_to_py(py, value)?)?;
    }
    Ok(out)
}

/// `array.array('d')`.
fn real_array<'py>(py: Python<'py>, values: &[f64]) -> PyResult<Bound<'py, PyAny>> {
    let bytes: Vec<u8> = values.iter().flat_map(|x| x.to_ne_bytes()).collect();
    new_array(py, "d", &bytes)
}

/// `array.array('i')` from signed values; a value beyond int32 is an error.
fn signed_array<'py>(
    py: Python<'py>,
    name: &'static str,
    values: impl IntoIterator<Item = i64>,
) -> PyResult<Bound<'py, PyAny>> {
    let mut bytes = Vec::new();
    for (position, value) in values.into_iter().enumerate() {
        let value = i32::try_from(value).map_err(|_| {
            SolveError::new_err(format!(
                "internal error: {name}[{position}] = {value} does not fit the int32 output array"
            ))
        })?;
        bytes.extend_from_slice(&value.to_ne_bytes());
    }
    new_array(py, "i", &bytes)
}

fn u32_array<'py>(
    py: Python<'py>,
    name: &'static str,
    values: &[u32],
) -> PyResult<Bound<'py, PyAny>> {
    index_array(py, name, values.iter().copied())
}

/// The keywords of `step` and `run`.
const CONTROL_KEYS: &[&str] = &[
    "energy",
    "scale_weights",
    "feature_r_eps",
    "feature_smoothing",
    "length_weight",
    "ripple_wavelength",
    "penalty",
    "penalty_rule",
    "distortion_weights",
    "area_weights",
    "symmetry",
    "curvature_uniform",
    "curvature_per_triangle",
    "stages",
    "flux_recovery",
    "flux_smoothing",
    "sliver_aspect_ratio",
    "smoothing_passes",
    "length_step",
    "extension",
    "speed_scaling",
    "redistance",
    "step_rule",
    "time_scale",
    "max_motion",
    "tau_max",
    "tau",
    "min_crossing_offset",
];
/// The keywords of `run` only.
const STOP_KEYS: &[&str] = &[
    "max_steps",
    "motion_tolerance",
    "energy_tolerance",
    "energy_window",
];

fn choice<T: Copy>(key: &str, value: &Bound<'_, PyAny>, table: &[(&str, T)]) -> PyResult<T> {
    let text: String = value
        .extract()
        .map_err(|_| InputError::new_err(format!("{key}: expected a string")))?;
    table
        .iter()
        .find(|(name, _)| *name == text)
        .map(|(_, v)| *v)
        .ok_or_else(|| {
            let names: Vec<String> = table.iter().map(|(n, _)| format!("{n:?}")).collect();
            InputError::new_err(format!(
                "{key}: {text:?} is not one of {}",
                names.join(", ")
            ))
        })
}

fn number(key: &str, value: &Bound<'_, PyAny>) -> PyResult<f64> {
    value
        .extract::<f64>()
        .map_err(|_| InputError::new_err(format!("{key}: expected a number")))
}

fn count(key: &str, value: &Bound<'_, PyAny>) -> PyResult<usize> {
    value
        .extract::<usize>()
        .map_err(|_| InputError::new_err(format!("{key}: expected a non-negative integer")))
}

/// A number per region: a list or tuple of numbers, or any real buffer.
fn region_weights(name: &'static str, value: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    if value.is_instance_of::<PyList>() || value.is_instance_of::<PyTuple>() {
        return value
            .extract::<Vec<f64>>()
            .map_err(|_| InputError::new_err(format!("{name}: expected numbers")));
    }
    reals(name, value, 1, "region")
}

/// Real numbers of a buffer of any shape `(.., width)` or flat (the rows are checked later).
fn reals_any_shape(name: &'static str, value: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    with_view(name, value, |view| {
        let width = view.shape.last().copied().unwrap_or(1).max(1);
        decode_reals(view, width, "corner")
    })
}

/// The keyword options of `step` / `run`. An unknown keyword is a `TypeError`, as for a function
/// with a fixed signature; a value of `None` means "not given".
fn parse_controls(
    function: &str,
    kwargs: Option<&Bound<'_, PyDict>>,
    with_stop: bool,
) -> PyResult<Controls> {
    use ddg_variational_cuts as vc;
    let mut c = Controls::default();
    let Some(kwargs) = kwargs else {
        return Ok(c);
    };
    for (key, value) in kwargs.iter() {
        let key: String = key.extract()?;
        let known = CONTROL_KEYS.contains(&key.as_str())
            || (with_stop && STOP_KEYS.contains(&key.as_str()));
        if !known {
            return Err(PyTypeError::new_err(format!(
                "{function}() got an unexpected keyword argument '{key}'"
            )));
        }
        if value.is_none() {
            continue;
        }
        let k = key.as_str();
        match k {
            "energy" => {
                c.energy = choice(
                    k,
                    &value,
                    &[
                        ("dirichlet", EnergyKind::Dirichlet),
                        ("scaled_dirichlet", EnergyKind::ScaledDirichlet),
                        ("hencky", EnergyKind::Hencky),
                    ],
                )?
            }
            "scale_weights" => c.scale_weights = Some(reals("scale_weights", &value, 1, "vertex")?),
            "feature_r_eps" => c.feature_r_eps = number(k, &value)?,
            "feature_smoothing" => c.feature_smoothing = number(k, &value)?,
            "length_weight" => c.length_weight = Some(number(k, &value)?),
            "ripple_wavelength" => c.ripple_wavelength = Some(number(k, &value)?),
            "penalty" => c.penalty = Some(reals("penalty", &value, 1, "vertex")?),
            "penalty_rule" => {
                c.penalty_rule = Some(choice(
                    k,
                    &value,
                    &[
                        ("gradient", vc::PenaltyRule::Gradient),
                        ("paper_diagonal", vc::PenaltyRule::PaperDiagonal),
                    ],
                )?)
            }
            "distortion_weights" => {
                c.distortion_weights = Some(region_weights("distortion_weights", &value)?)
            }
            "area_weights" => c.area_weights = Some(region_weights("area_weights", &value)?),
            "symmetry" => c.symmetry = Some(indices("symmetry", &value, 1, "vertex")?),
            "curvature_uniform" => c.curvature_uniform = Some(number(k, &value)?),
            "curvature_per_triangle" => {
                c.curvature_per_triangle =
                    Some(reals("curvature_per_triangle", &value, 1, "triangle")?)
            }
            "stages" => {
                c.stages = Some(choice(
                    k,
                    &value,
                    &[
                        ("default", StageSet::Default),
                        ("paper", StageSet::Paper),
                        ("recommended", StageSet::Recommended),
                    ],
                )?)
            }
            "flux_recovery" => {
                c.flux_recovery = Some(choice(
                    k,
                    &value,
                    &[
                        ("nodal", FluxChoice::Nodal),
                        ("smoothed", FluxChoice::Smoothed),
                    ],
                )?)
            }
            "flux_smoothing" => c.flux_smoothing = Some(number(k, &value)?),
            "sliver_aspect_ratio" => c.sliver_aspect_ratio = Some(number(k, &value)?),
            "smoothing_passes" => c.smoothing_passes = Some(count(k, &value)?),
            "length_step" => {
                c.length_step = Some(choice(
                    k,
                    &value,
                    &[
                        ("on_cut", vc::LengthStep::OnCut),
                        ("on_mesh", vc::LengthStep::OnMesh),
                    ],
                )?)
            }
            "extension" => {
                c.extension = Some(choice(
                    k,
                    &value,
                    &[
                        ("local", vc::Extension::Local),
                        ("harmonic", vc::Extension::Harmonic),
                    ],
                )?)
            }
            "speed_scaling" => {
                c.speed_scaling = Some(choice(
                    k,
                    &value,
                    &[
                        ("mean_gradient", vc::SpeedScaling::MeanGradient),
                        ("unit", vc::SpeedScaling::Unit),
                    ],
                )?)
            }
            "redistance" => {
                c.redistance = Some(choice(
                    k,
                    &value,
                    &[
                        ("keep_adjacent", vc::Redistance::KeepAdjacent),
                        ("off", vc::Redistance::Off),
                    ],
                )?)
            }
            "step_rule" => {
                c.step_rule = choice(
                    k,
                    &value,
                    &[
                        ("adaptive", StepRule::Adaptive),
                        ("capped", StepRule::Capped),
                        ("fixed", StepRule::Fixed),
                    ],
                )?
            }
            "time_scale" => c.time_scale = Some(number(k, &value)?),
            "max_motion" => c.max_motion = Some(number(k, &value)?),
            "tau_max" => c.tau_max = Some(number(k, &value)?),
            "tau" => c.tau = Some(number(k, &value)?),
            "min_crossing_offset" => c.min_crossing_offset = Some(number(k, &value)?),
            "max_steps" => c.max_steps = Some(count(k, &value)?),
            "motion_tolerance" => c.motion_tolerance = Some(number(k, &value)?),
            "energy_tolerance" => c.energy_tolerance = Some(number(k, &value)?),
            "energy_window" => c.energy_window = Some(count(k, &value)?),
            _ => unreachable!("checked against the key lists"),
        }
    }
    Ok(c)
}

/// The start from the keywords of the constructor / `restart`: exactly one of a circle
/// (`circle_vertex` with `circle_radius`), `polygon_labels` or `corner_values`.
fn parse_start(
    circle_vertex: Option<u32>,
    circle_radius: Option<f64>,
    polygon_labels: Option<&Bound<'_, PyAny>>,
    corner_values: Option<&Bound<'_, PyAny>>,
    cluster_iterations: Option<usize>,
) -> PyResult<StartSpec> {
    let n = usize::from(circle_vertex.is_some() || circle_radius.is_some())
        + usize::from(polygon_labels.is_some())
        + usize::from(corner_values.is_some())
        + usize::from(cluster_iterations.is_some());
    if n != 1 {
        return Err(InputError::new_err(
            "give exactly one start: circle_vertex with circle_radius, polygon_labels, \
             corner_values, or cluster_iterations (normal clusters)",
        ));
    }
    if let Some(max_iterations) = cluster_iterations {
        return Ok(StartSpec::NormalClusters { max_iterations });
    }
    if let Some(labels) = polygon_labels {
        return Ok(StartSpec::Labels(indices(
            "polygon_labels",
            labels,
            1,
            "polygon",
        )?));
    }
    if let Some(values) = corner_values {
        return Ok(StartSpec::Corners(reals_any_shape(
            "corner_values",
            values,
        )?));
    }
    match (circle_vertex, circle_radius) {
        (Some(vertex), Some(radius)) => Ok(StartSpec::Circle { vertex, radius }),
        _ => Err(InputError::new_err(
            "circle_vertex and circle_radius belong together",
        )),
    }
}

/// The cut on the caller's polygons: points on polygon edges and, per polygon, the paths of the
/// cut through it. Plain arrays, all flat (see the crate README for the layout).
#[pyclass(module = "ddg_tools", frozen)]
struct CutDescription {
    cut: ddg_uv::cutting::CutOnPolygons,
    /// `array('i')`, 2 per edge point: the vertex pair `[min, max]` of its polygon edge.
    #[pyo3(get)]
    edge_vertices: Py<PyAny>,
    /// `array('d')`: the parameter from the first vertex of the edge.
    #[pyo3(get)]
    edge_t: Py<PyAny>,
    /// `array('d')`, 3 per edge point.
    #[pyo3(get)]
    edge_positions: Py<PyAny>,
    /// `array('i')`: for a point on a seam the polygon whose side it is on, else -1.
    #[pyo3(get)]
    edge_side: Py<PyAny>,
    /// `array('i')`: the polygon of every path.
    #[pyo3(get)]
    path_polygon: Py<PyAny>,
    /// `array('i')`: the start of every path: an edge point, or `-(j + 1)` for junction `j`.
    #[pyo3(get)]
    path_start: Py<PyAny>,
    /// `array('i')`: the end of every path, coded like the start.
    #[pyo3(get)]
    path_end: Py<PyAny>,
    /// `array('i')`, 2 per path: the regions on its left and on its right.
    #[pyo3(get)]
    path_regions: Py<PyAny>,
    /// `array('i')`, `n_paths + 1` offsets into the interior points.
    #[pyo3(get)]
    path_interior_starts: Py<PyAny>,
    /// `array('i')`, 2 per interior point: the polygon corners of the diagonal it lies on.
    #[pyo3(get)]
    interior_corners: Py<PyAny>,
    /// `array('d')`: the parameter from the first corner.
    #[pyo3(get)]
    interior_t: Py<PyAny>,
    /// `array('d')`, 3 per interior point.
    #[pyo3(get)]
    interior_positions: Py<PyAny>,
    /// `array('d')`, 3 per junction of three or more regions (empty for two).
    #[pyo3(get)]
    junction_positions: Py<PyAny>,
    /// `array('i')`: the polygon of every junction.
    #[pyo3(get)]
    junction_polygon: Py<PyAny>,
    /// `array('i')`: the triangle of the supplied triangulation that contains each junction.
    #[pyo3(get)]
    junction_triangle: Py<PyAny>,
    /// `array('i')`, 3 per junction: the polygon corners of its triangle.
    #[pyo3(get)]
    junction_corners: Py<PyAny>,
    /// `array('d')`, 3 per junction: barycentric weights of those corners (the junction is
    /// `sum(w[k] * x[corner k])`).
    #[pyo3(get)]
    junction_weights: Py<PyAny>,
    /// Number of crossings that were moved away from a vertex.
    #[pyo3(get)]
    n_clamped_crossings: usize,
    /// What the simplification of the cut did (intrinsic triangulations; empty otherwise):
    /// `n_merged_nodes` (bends removed), `max_displacement` (a length, an upper bound on how far
    /// the described cut is from the flow's), `max_displacement_relative` (as a fraction of the
    /// longest edge of the triangle), `n_moved_points`.
    #[pyo3(get)]
    report: Py<PyDict>,
}

impl CutDescription {
    fn new(
        py: Python<'_>,
        cut: ddg_uv::cutting::CutOnPolygons,
        report: Option<&Dict>,
    ) -> PyResult<CutDescription> {
        let a = DescriptionArrays::from_cut(&cut);
        let report = match report {
            Some(d) => dict_to_py(py, d)?,
            None => PyDict::new(py),
        };
        Ok(CutDescription {
            report: report.unbind(),
            edge_vertices: u32_array(py, "edge_vertices", &a.edge_vertices)?.unbind(),
            edge_t: real_array(py, &a.edge_t)?.unbind(),
            edge_positions: real_array(py, &a.edge_positions)?.unbind(),
            edge_side: signed_array(py, "edge_side", a.edge_side.iter().copied())?.unbind(),
            path_polygon: u32_array(py, "path_polygon", &a.path_polygon)?.unbind(),
            path_start: signed_array(py, "path_start", a.path_start.iter().copied())?.unbind(),
            path_end: signed_array(py, "path_end", a.path_end.iter().copied())?.unbind(),
            path_regions: u32_array(py, "path_regions", &a.path_regions)?.unbind(),
            path_interior_starts: u32_array(py, "path_interior_starts", &a.path_interior_starts)?
                .unbind(),
            interior_corners: u32_array(py, "interior_corners", &a.interior_corners)?.unbind(),
            interior_t: real_array(py, &a.interior_t)?.unbind(),
            interior_positions: real_array(py, &a.interior_positions)?.unbind(),
            junction_positions: real_array(py, &a.junction_positions)?.unbind(),
            junction_polygon: u32_array(py, "junction_polygon", &a.junction_polygon)?.unbind(),
            junction_triangle: u32_array(py, "junction_triangle", &a.junction_triangle)?.unbind(),
            junction_corners: u32_array(py, "junction_corners", &a.junction_corners)?.unbind(),
            junction_weights: real_array(py, &a.junction_weights)?.unbind(),
            n_clamped_crossings: a.n_clamped_crossings,
            cut,
        })
    }
}

#[pymethods]
impl CutDescription {
    /// The cut as lines for display, joined at the points where exactly two paths meet:
    /// `(positions, starts, closed, regions)` as `CutSession.preview` returns them.
    fn lines<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        lines_tuple(py, &self.cut)
    }
    /// Number of points where the cut crosses a polygon edge.
    #[getter]
    fn n_edge_points(&self) -> usize {
        self.cut.edge_points.len()
    }
    /// Number of paths (pieces of the cut inside one polygon).
    #[getter]
    fn n_paths(&self) -> usize {
        self.cut.paths.len()
    }
    fn __repr__(&self) -> String {
        format!(
            "CutDescription(edge_points={}, paths={}, junctions={})",
            self.cut.edge_points.len(),
            self.cut.paths.len(),
            self.cut.junctions.len()
        )
    }
}

/// The polygon mesh with a cut inserted: the polygons the cut crosses are split along it.
#[pyclass(module = "ddg_tools", frozen)]
struct InsertedMesh {
    /// `array('d')`, 3 per vertex.
    #[pyo3(get)]
    positions: Py<PyAny>,
    /// `array('i')`: the vertex of every polygon corner.
    #[pyo3(get)]
    corner_vertices: Py<PyAny>,
    /// `array('i')`: the first corner of every polygon.
    #[pyo3(get)]
    polygon_starts: Py<PyAny>,
    /// `array('i')`, 3 corner indices per triangle of the new polygons.
    #[pyo3(get)]
    triangles: Py<PyAny>,
    /// `array('i')`, 2 vertices `[min, max]` per edge of the cut, ascending.
    #[pyo3(get)]
    new_seams: Py<PyAny>,
    /// `array('i')`: the old polygon each new polygon is a part of.
    #[pyo3(get)]
    polygon_parent: Py<PyAny>,
    /// `array('i')`, 3 per new vertex: kind (0 old vertex, 1 edge point, 2 interior point of a
    /// path, 3 junction), then the old vertex / the edge point / the path / the junction, and for kind 2 the position in
    /// the path's interior points (else 0).
    #[pyo3(get)]
    vertex_source: Py<PyAny>,
    /// `array('i')`: the vertex of every edge point of the description.
    #[pyo3(get)]
    edge_point_vertex: Py<PyAny>,
    /// `array('i')`: the vertex of every junction of the description.
    #[pyo3(get)]
    junction_vertex: Py<PyAny>,
    /// `array('i')`: the region of every new polygon where the description says so (a part
    /// next to a piece of the cut), else -1.
    #[pyo3(get)]
    polygon_region: Py<PyAny>,
    /// `array('i')`: the new polygons that were split or got extra vertices and have more than
    /// three corners (the ones whose triangulation an editor should keep), ascending.
    #[pyo3(get)]
    changed_polygons: Py<PyAny>,
    /// `array('i')`, `len + 1` offsets into `changed_vertices`.
    #[pyo3(get)]
    changed_vertex_starts: Py<PyAny>,
    /// `array('i')`: the vertices of every changed polygon, in order.
    #[pyo3(get)]
    changed_vertices: Py<PyAny>,
    /// `array('i')`, `len + 1` offsets, in triangles, into `changed_triangles`.
    #[pyo3(get)]
    changed_triangle_starts: Py<PyAny>,
    /// `array('i')`: 3 vertices per triangle of the changed polygons (the library's
    /// triangulation).
    #[pyo3(get)]
    changed_triangles: Py<PyAny>,
    /// Number of vertices of the old mesh (vertices below it keep their index).
    #[pyo3(get)]
    n_old_vertices: usize,
}

impl InsertedMesh {
    fn new(py: Python<'_>, a: &InsertedArrays) -> PyResult<InsertedMesh> {
        Ok(InsertedMesh {
            positions: real_array(py, &a.positions)?.unbind(),
            corner_vertices: u32_array(py, "corner_vertices", &a.corner_vertices)?.unbind(),
            polygon_starts: u32_array(py, "polygon_starts", &a.polygon_starts)?.unbind(),
            triangles: u32_array(py, "triangles", &a.triangles)?.unbind(),
            new_seams: u32_array(py, "new_seams", &a.new_seams)?.unbind(),
            polygon_parent: u32_array(py, "polygon_parent", &a.polygon_parent)?.unbind(),
            vertex_source: u32_array(py, "vertex_source", &a.vertex_source)?.unbind(),
            edge_point_vertex: u32_array(py, "edge_point_vertex", &a.edge_point_vertex)?.unbind(),
            junction_vertex: u32_array(py, "junction_vertex", &a.junction_vertex)?.unbind(),
            polygon_region: signed_array(py, "polygon_region", a.polygon_region.iter().copied())?
                .unbind(),
            changed_polygons: u32_array(py, "changed_polygons", &a.changed_polygons)?.unbind(),
            changed_vertex_starts: u32_array(
                py,
                "changed_vertex_starts",
                &a.changed_vertex_starts,
            )?
            .unbind(),
            changed_vertices: u32_array(py, "changed_vertices", &a.changed_vertices)?.unbind(),
            changed_triangle_starts: u32_array(
                py,
                "changed_triangle_starts",
                &a.changed_triangle_starts,
            )?
            .unbind(),
            changed_triangles: u32_array(py, "changed_triangles", &a.changed_triangles)?.unbind(),
            n_old_vertices: a.n_old_vertices,
        })
    }
}

#[pymethods]
impl InsertedMesh {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "InsertedMesh(vertices={}, polygons={}, new_seams={})",
            self.positions.bind(py).len()? / 3,
            self.polygon_starts.bind(py).len()?,
            self.new_seams.bind(py).len()? / 2
        ))
    }
}

/// The result of `CutSession.unwrap`: the inserted mesh and the UV map of the new mesh.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct CutUnwrapResult {
    /// The cut on the polygons that was inserted.
    cut: Py<CutDescription>,
    /// The polygon mesh with the cut inserted.
    inserted: Py<InsertedMesh>,
    /// `array('i')`, 2 per seam of the new mesh: the old seams and the cut.
    seams: Py<PyAny>,
    /// `array('f')`, 2 per corner of the new mesh. Not packed.
    uvs: Py<PyAny>,
    /// `array('i')`: the chart of every polygon of the new mesh.
    polygon_charts: Py<PyAny>,
    /// `array('i')`: vertices of the new mesh that carry a cone.
    cone_vertices: Py<PyAny>,
    /// `array('d')`: their angles.
    cone_angles: Py<PyAny>,
    /// `array('i')`, 2 per edge: the edges the unwrap cut along, seams included.
    cut_edges: Py<PyAny>,
    /// `dict` of the unwrap report.
    report: Py<PyDict>,
}

fn lines_tuple<'py>(
    py: Python<'py>,
    cut: &ddg_uv::cutting::CutOnPolygons,
) -> PyResult<Bound<'py, PyTuple>> {
    let p = crate::cutting::PreviewArrays::from_cut(cut);
    PyTuple::new(
        py,
        [
            real_array(py, &p.positions)?,
            u32_array(py, "starts", &p.starts)?,
            PyList::new(py, &p.closed)?.into_any(),
            u32_array(py, "regions", &p.regions)?,
        ],
    )
}

/// The bytes of a state: `bytes`, `bytearray`, `memoryview` of bytes.
fn state_bytes(name: &'static str, obj: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    obj.extract::<Vec<u8>>()
        .map_err(|_| InputError::new_err(format!("{name}: expected bytes (from snapshot())")))
}

fn poisoned() -> PyErr {
    SolveError::new_err("internal error: a previous call panicked inside this session")
}

/// The cut flow of "Variational Surface Cutting" (Sharp, Crane 2018) on a polygon mesh.
///
/// A session owns the state of the flow: one function per region and per polygon corner (for
/// two regions `phi` and `-phi`, region 0 where `phi <= 0`). Build it from the mesh buffers of
/// `unwrap` and a start; advance it with `step` / `run`; read the cut with `preview`,
/// `describe_cut`, `insert` and `unwrap`. Weights, penalties and step sizes refer to the surface
/// scaled to unit area (`length_scale`); the state and the positions are in the unit of the mesh.
///
/// Mesh arguments (as for `unwrap`): positions, corner_vertices, polygon_starts, triangles.
/// seam_edges, pinned_vertices, forbidden_vertices, prescribed_vertices, prescribed_angles:
/// steering as for `unwrap`; seams are cut open before the flow runs (the flow cannot cross
/// them) and stay seams. n_regions: 2 is the only value the flow handles so far; more raises
/// NotAvailableError.
///
/// Exactly one start: circle_vertex with circle_radius (the geodesic disk, region 0);
/// polygon_labels (one label per polygon below n_regions); corner_values (row-major
/// n_regions x n_corners, as `state` returns them).
///
/// A general pointwise energy given as a Python callable is not offered.
#[pyclass(module = "ddg_tools", frozen)]
struct CutSession {
    inner: Mutex<RustCutSession>,
}

impl CutSession {
    fn lock(&self) -> PyResult<MutexGuard<'_, RustCutSession>> {
        self.inner.lock().map_err(|_| poisoned())
    }
}

#[pymethods]
impl CutSession {
    #[new]
    #[pyo3(signature = (
        positions,
        corner_vertices,
        polygon_starts,
        triangles,
        *,
        n_regions = None,
        triangulation = "intrinsic",
        refine_bound = None,
        optimal_delaunay = false,
        representation = None,
        min_offset = None,
        merge_tolerance = None,
        state = None,
        seam_edges = None,
        pinned_vertices = None,
        forbidden_vertices = None,
        prescribed_vertices = None,
        prescribed_angles = None,
        circle_vertex = None,
        circle_radius = None,
        polygon_labels = None,
        corner_values = None,
        cluster_iterations = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new<'py>(
        py: Python<'py>,
        positions: &Bound<'py, PyAny>,
        corner_vertices: &Bound<'py, PyAny>,
        polygon_starts: &Bound<'py, PyAny>,
        triangles: &Bound<'py, PyAny>,
        n_regions: Option<usize>,
        triangulation: &str,
        refine_bound: Option<f64>,
        optimal_delaunay: bool,
        representation: Option<&str>,
        min_offset: Option<f64>,
        merge_tolerance: Option<f64>,
        state: Option<&Bound<'py, PyAny>>,
        seam_edges: Option<&Bound<'py, PyAny>>,
        pinned_vertices: Option<&Bound<'py, PyAny>>,
        forbidden_vertices: Option<&Bound<'py, PyAny>>,
        prescribed_vertices: Option<&Bound<'py, PyAny>>,
        prescribed_angles: Option<&Bound<'py, PyAny>>,
        circle_vertex: Option<u32>,
        circle_radius: Option<f64>,
        polygon_labels: Option<&Bound<'py, PyAny>>,
        corner_values: Option<&Bound<'py, PyAny>>,
        cluster_iterations: Option<usize>,
    ) -> PyResult<CutSession> {
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
        let spec = match (triangulation, refine_bound) {
            ("caller", None) => TriangulationSpec::Caller,
            ("intrinsic", None) => TriangulationSpec::Flips,
            ("refined", bound) => TriangulationSpec::Refined {
                bound_degrees: bound.unwrap_or(25.0),
            },
            ("caller" | "intrinsic", Some(_)) => {
                return Err(InputError::new_err(
                    "refine_bound only applies to triangulation=\"refined\"",
                ));
            }
            (other, _) => {
                return Err(InputError::new_err(format!(
                    "triangulation: {other:?} is not one of \"intrinsic\", \"caller\", \
                     \"refined\""
                )));
            }
        };
        let representation = match representation {
            None => None,
            Some("integer") => Some(ddg_uv::cutting::Representation::Integer),
            Some("signposts") => Some(ddg_uv::cutting::Representation::Signposts),
            Some(other) => {
                return Err(InputError::new_err(format!(
                    "representation: {other:?} is not one of \"integer\", \"signposts\""
                )));
            }
        };
        if state.is_some() && representation.is_some() {
            return Err(InputError::new_err(
                "state carries its own triangulation, its representation included: do not pass \
                 representation with it (read CutSession.representation afterwards)",
            ));
        }
        if state.is_some() && optimal_delaunay {
            return Err(InputError::new_err(
                "state carries its own triangulation, optimal_delaunay included: do not pass \
                 optimal_delaunay with it (read CutSession.optimal_delaunay afterwards)",
            ));
        }
        let session = if let Some(state) = state {
            if circle_vertex.is_some()
                || circle_radius.is_some()
                || polygon_labels.is_some()
                || corner_values.is_some()
                || cluster_iterations.is_some()
            {
                return Err(InputError::new_err(
                    "state replaces the start: give no other start with it",
                ));
            }
            let bytes = state_bytes("state", state)?;
            py.detach(move || {
                let input = UnwrapInput::from_buffers(buffers)
                    .map_err(|e| BindError::Input(e.to_string()))?;
                let session =
                    RustCutSession::from_state(&input, &bytes, min_offset, merge_tolerance)?;
                if n_regions.is_some_and(|n| n != session.n_regions()) {
                    return Err(BindError::Input(format!(
                        "n_regions={} but the state has {} regions",
                        n_regions.unwrap_or(0),
                        session.n_regions()
                    )));
                }
                Ok(session)
            })
            .map_err(bind_error)?
        } else {
            let n_regions = n_regions.unwrap_or(2);
            let start = parse_start(
                circle_vertex,
                circle_radius,
                polygon_labels,
                corner_values,
                cluster_iterations,
            )?;
            py.detach(move || {
                let input = UnwrapInput::from_buffers(buffers)
                    .map_err(|e| BindError::Input(e.to_string()))?;
                RustCutSession::new_on(
                    &input,
                    n_regions,
                    &start,
                    &spec,
                    optimal_delaunay,
                    representation,
                    min_offset,
                    merge_tolerance,
                )
            })
            .map_err(bind_error)?
        };
        Ok(CutSession {
            inner: Mutex::new(session),
        })
    }

    /// The triangulation option the session runs on: "intrinsic" (Delaunay flips), "refined" or
    /// "caller".
    #[getter]
    fn triangulation(&self) -> PyResult<&'static str> {
        Ok(self.lock()?.triangulation_name())
    }
    /// The representation of the intrinsic triangulation: "integer" (the default; Gillespie,
    /// Sharp, Crane 2021) or "signposts" (Sharp, Soliman, Crane 2019); None on "caller". After
    /// `state=` this is the state's.
    #[getter]
    fn representation(&self) -> PyResult<Option<&'static str>> {
        Ok(self.lock()?.representation_name())
    }
    /// Whether the intrinsic optimal Delaunay triangulation is on (False on "caller"). After
    /// `state=` this is the flag the state was made with.
    #[getter]
    fn optimal_delaunay(&self) -> PyResult<bool> {
        Ok(self.lock()?.optimal_delaunay())
    }
    /// Number of patches of the current state: connected pieces of the surface in one region
    /// (corners joined along polygons and across edges that are not seams).
    #[getter]
    fn n_patches(&self) -> PyResult<usize> {
        Ok(self.lock()?.n_patches())
    }
    /// Whether the per-corner `state` is the whole state (the triangulation added no vertex).
    /// If not, keep `snapshot()` as well.
    #[getter]
    fn state_is_per_corner(&self) -> PyResult<bool> {
        Ok(self.lock()?.state_is_per_corner())
    }
    /// Hex fingerprint of the triangulation the flow runs on (what a stored state must match).
    #[getter]
    fn fingerprint(&self) -> PyResult<String> {
        Ok(self.lock()?.fingerprint())
    }
    /// What the construction of the intrinsic triangulation did, as a dict (n_added,
    /// n_initial_flips, n_splits, min_angle_degrees, max_angle_degrees, is_delaunay, seconds,
    /// ...); None on the caller's triangulation.
    #[getter]
    fn intrinsic_report<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let report = self.lock()?.intrinsic_report();
        report.map(|d| dict_to_py(py, &d)).transpose()
    }
    /// The whole state as `bytes` (versioned, with the fingerprint of the triangulation, the
    /// number of regions and the triangulation option). Needed when `state_is_per_corner` is
    /// false; always valid.
    fn snapshot<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let bytes = self.lock()?.snapshot();
        Ok(PyBytes::new(py, &bytes))
    }
    /// Sets the state from `bytes` made by `snapshot()`. Raises StaleStateError if the mesh or
    /// its seams are not the ones the state was made for.
    fn restore(&self, py: Python<'_>, state: &Bound<'_, PyAny>) -> PyResult<()> {
        let bytes = state_bytes("state", state)?;
        py.detach(|| self.lock()?.restore(&bytes).map_err(bind_error))
    }

    /// Number of regions.
    #[getter]
    fn n_regions(&self) -> PyResult<usize> {
        Ok(self.lock()?.n_regions())
    }
    /// Number of vertices of the polygon mesh.
    #[getter]
    fn n_vertices(&self) -> PyResult<usize> {
        Ok(self.lock()?.n_vertices())
    }
    /// Number of polygon corners ("loops"): the length of every row of `state`.
    #[getter]
    fn n_corners(&self) -> PyResult<usize> {
        Ok(self.lock()?.n_corners())
    }
    /// Number of polygons.
    #[getter]
    fn n_polygons(&self) -> PyResult<usize> {
        Ok(self.lock()?.n_polygons())
    }
    /// Number of triangles of the supplied triangulation.
    #[getter]
    fn n_triangles(&self) -> PyResult<usize> {
        Ok(self.lock()?.n_triangles())
    }
    /// The factor from a length of the mesh to a length on the unit-area surface that the
    /// weights, the penalty, the step sizes and the ripple wavelength refer to.
    #[getter]
    fn length_scale(&self) -> PyResult<f64> {
        Ok(self.lock()?.length_scale())
    }
    /// Total area of the mesh, in the unit of the mesh squared.
    #[getter]
    fn total_area(&self) -> PyResult<f64> {
        Ok(self.lock()?.total_area())
    }
    /// Flow time since the last start / `state` assignment.
    #[getter]
    fn time(&self) -> PyResult<f64> {
        Ok(self.lock()?.time())
    }
    /// Steps taken since the last start / `state` assignment.
    #[getter]
    fn n_steps(&self) -> PyResult<usize> {
        Ok(self.lock()?.n_steps())
    }
    /// Warnings of the validation of the mesh (ignored or extended seams), `list[str]`.
    #[getter]
    fn warnings(&self) -> PyResult<Vec<String>> {
        Ok(self.lock()?.warnings().to_vec())
    }
    /// `array('i')`, 2 per edge: the seams that were cut open before the flow ran.
    #[getter]
    fn seams_cut<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let seams: Vec<u32> = self.lock()?.seams_cut().iter().flatten().copied().collect();
        u32_array(py, "seams_cut", &seams)
    }

    /// The state: a `memoryview` of float64 with shape (n_regions, n_corners), a copy
    /// (`numpy.asarray(session.state)` gives the 2-D array). In the unit of the mesh; for two
    /// regions row 1 is minus row 0. Assigning accepts any float buffer of that shape (or flat,
    /// row-major) and resets time and step count.
    #[getter(state)]
    fn get_state<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let (flat, n_regions, n_corners) = {
            let g = self.lock()?;
            (g.state(), g.n_regions(), g.n_corners())
        };
        let array = real_array(py, &flat)?;
        let view = py
            .import("builtins")?
            .getattr("memoryview")?
            .call1((array,))?;
        let bytes = view.call_method1("cast", ("B",))?;
        bytes.call_method1("cast", ("d", (n_regions, n_corners)))
    }

    #[setter(state)]
    fn set_state(&self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let n_corners = self.lock()?.n_corners();
        let flat = reals("state", value, n_corners, "corner")?;
        py.detach(|| self.lock()?.set_state(&flat).map_err(bind_error))
    }

    /// One row of the state as `array('d')`.
    fn state_row<'py>(&self, py: Python<'py>, region: usize) -> PyResult<Bound<'py, PyAny>> {
        let (flat, n_regions, n_corners) = {
            let g = self.lock()?;
            (g.state(), g.n_regions(), g.n_corners())
        };
        if region >= n_regions {
            return Err(InputError::new_err(format!(
                "region {region} is out of range ({n_regions} regions)"
            )));
        }
        real_array(py, &flat[region * n_corners..(region + 1) * n_corners])
    }

    /// The region of every polygon corner, `array('i')`.
    fn corner_regions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let regions = self.lock()?.corner_regions();
        u32_array(py, "corner_regions", &regions)
    }

    /// Replaces the state by a new start; the keywords are those of the constructor.
    #[pyo3(signature = (*, circle_vertex = None, circle_radius = None, polygon_labels = None, corner_values = None, cluster_iterations = None))]
    fn restart(
        &self,
        py: Python<'_>,
        circle_vertex: Option<u32>,
        circle_radius: Option<f64>,
        polygon_labels: Option<&Bound<'_, PyAny>>,
        corner_values: Option<&Bound<'_, PyAny>>,
        cluster_iterations: Option<usize>,
    ) -> PyResult<()> {
        let start = parse_start(
            circle_vertex,
            circle_radius,
            polygon_labels,
            corner_values,
            cluster_iterations,
        )?;
        py.detach(|| self.lock()?.restart(&start).map_err(bind_error))
    }

    /// The roughening rate of the current state: the scale of the length weight. Takes the
    /// energy keywords of `step`.
    #[pyo3(signature = (**kwargs))]
    fn roughening_rate(&self, py: Python<'_>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<f64> {
        let controls = parse_controls("roughening_rate", kwargs, false)?;
        py.detach(|| self.lock()?.roughening_rate(&controls).map_err(bind_error))
    }

    /// One step of the flow. Returns a dict (energy terms of the state the step started from,
    /// `tau`, speeds, counters, timings, `length_weight`).
    ///
    /// Keywords (defaults are the library's; see the crate README for the table):
    /// energy ("dirichlet", "scaled_dirichlet", "hencky"; default "hencky"), scale_weights
    /// (per vertex, with "scaled_dirichlet"; else the curvature-based weight with feature_r_eps
    /// and feature_smoothing), length_weight or ripple_wavelength (exclusive; a wavelength sets
    /// the weight from the roughening rate of the state; neither: ripple_wavelength 0.25),
    /// penalty (eta per vertex) and penalty_rule ("gradient", "paper_diagonal"),
    /// distortion_weights and area_weights (per region), symmetry (class per vertex),
    /// curvature_uniform or curvature_per_triangle, stages ("default" or "paper") and the
    /// individual stages flux_recovery, flux_smoothing, sliver_aspect_ratio, smoothing_passes,
    /// length_step, extension, speed_scaling, redistance, step_rule ("adaptive", "capped",
    /// "fixed") with time_scale, max_motion, tau_max, tau, and min_crossing_offset.
    #[pyo3(signature = (**kwargs))]
    fn step<'py>(
        &self,
        py: Python<'py>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let controls = parse_controls("step", kwargs, false)?;
        let report = py.detach(|| self.lock()?.step(&controls).map_err(bind_error))?;
        dict_to_py(py, &report)
    }

    /// Steps until the stop rule says so. Takes the keywords of `step` and max_steps (default
    /// 1000), motion_tolerance (1e-5), energy_tolerance (1e-4), energy_window (25); a tolerance
    /// of 0 never stops for that reason. A ripple_wavelength sets the length weight once, at
    /// the start of the run. Returns a dict: steps_taken, stop_reason ("max_steps",
    /// "converged", "energy_stalled", "cut_vanished"), energy_before, energy_after,
    /// cut_length_before/after (in the unit of the mesh), per-region distortion and area,
    /// length_weight, and `steps`, the list of the dicts that `step` returns.
    #[pyo3(signature = (**kwargs))]
    fn run<'py>(
        &self,
        py: Python<'py>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let controls = parse_controls("run", kwargs, true)?;
        let report = py.detach(|| self.lock()?.run(&controls).map_err(bind_error))?;
        dict_to_py(py, &report)
    }

    /// The cut as polylines for display: `(positions, starts, closed)` with `positions` an
    /// `array('d')` of 3 per point, `starts` an `array('i')` of `n_polylines + 1` offsets into
    /// the points, and `closed` a `list[bool]`. An open polyline ends on the boundary or on a
    /// seam.
    #[pyo3(signature = (*, min_crossing_offset = None))]
    fn preview<'py>(
        &self,
        py: Python<'py>,
        min_crossing_offset: Option<f64>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let p = py.detach(|| {
            self.lock()?
                .preview(min_crossing_offset)
                .map_err(bind_error)
        })?;
        PyTuple::new(
            py,
            [
                real_array(py, &p.positions)?,
                u32_array(py, "starts", &p.starts)?,
                PyList::new(py, &p.closed)?.into_any(),
                u32_array(py, "regions", &p.regions)?,
            ],
        )
    }

    /// The current cut as lines for display, from the description (`CutDescription.lines()`).
    #[pyo3(signature = (*, min_crossing_offset = None))]
    fn lines<'py>(
        &self,
        py: Python<'py>,
        min_crossing_offset: Option<f64>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let cut = py.detach(|| {
            self.lock()?
                .describe_cut(min_crossing_offset)
                .map_err(bind_error)
        })?;
        lines_tuple(py, &cut)
    }

    /// The current cut expressed on the polygons of the mesh: a CutDescription.
    #[pyo3(signature = (*, min_crossing_offset = None))]
    fn describe_cut(
        &self,
        py: Python<'_>,
        min_crossing_offset: Option<f64>,
    ) -> PyResult<CutDescription> {
        let (cut, report) = py.detach(|| {
            self.lock()?
                .describe_cut_report(min_crossing_offset)
                .map_err(bind_error)
        })?;
        CutDescription::new(py, cut, Some(&report))
    }

    /// Inserts a cut into the polygon mesh (the Rust reference of what an editor must produce):
    /// `cut` is a CutDescription, or the current cut if omitted. Returns an InsertedMesh.
    #[pyo3(signature = (cut = None, *, min_crossing_offset = None))]
    fn insert(
        &self,
        py: Python<'_>,
        cut: Option<&Bound<'_, CutDescription>>,
        min_crossing_offset: Option<f64>,
    ) -> PyResult<InsertedMesh> {
        let cut = cut.map(|c| c.get().cut.clone());
        let inserted = py.detach(|| {
            self.lock()?
                .insert(cut.as_ref(), min_crossing_offset)
                .map_err(bind_error)
        })?;
        let old = self.lock()?.polygon_offsets();
        InsertedMesh::new(
            py,
            &InsertedArrays::from_inserted(&inserted).with_changes(&old, &inserted),
        )
    }

    /// Describes the current cut, inserts it, and unwraps the new mesh with the old seams and
    /// the cut as seams. cones=False (default): no cone is placed besides the cut; True: the
    /// cone solver may place cones (lambda_, positive_cones_only, optimize_boundary as for
    /// `unwrap`). The UVs are not packed (`pack_charts`).
    #[pyo3(signature = (*, min_crossing_offset = None, lambda_ = 0.1, positive_cones_only = true, optimize_boundary = false, cones = false))]
    fn unwrap(
        &self,
        py: Python<'_>,
        min_crossing_offset: Option<f64>,
        lambda_: f64,
        positive_cones_only: bool,
        optimize_boundary: bool,
        cones: bool,
    ) -> PyResult<CutUnwrapResult> {
        let settings = UnwrapSettings {
            lambda: lambda_,
            positive_cones_only,
            optimize_boundary,
            cones,
        };
        let (cut, result) = py.detach(|| {
            self.lock()?
                .unwrap(min_crossing_offset, &settings)
                .map_err(bind_error)
        })?;
        let a = UnwrapArrays::from_unwrap(&result);
        let report = dict_to_py(py, &a.report)?;
        report.set_item("warnings", PyList::new(py, &a.warnings)?)?;
        Ok(CutUnwrapResult {
            cut: Py::new(py, CutDescription::new(py, cut, None)?)?,
            inserted: Py::new(py, InsertedMesh::new(py, &a.inserted)?)?,
            seams: u32_array(py, "seams", &a.seams)?.unbind(),
            uvs: uv_array(py, &a.uvs)?.unbind(),
            polygon_charts: u32_array(py, "polygon_charts", &a.polygon_charts)?.unbind(),
            cone_vertices: u32_array(py, "cone_vertices", &a.cone_vertices)?.unbind(),
            cone_angles: real_array(py, &a.cone_angles)?.unbind(),
            cut_edges: u32_array(py, "cut_edges", &a.cut_edges)?.unbind(),
            report: report.unbind(),
        })
    }

    fn __repr__(&self) -> PyResult<String> {
        let g = self.lock()?;
        Ok(format!(
            "CutSession(regions={}, vertices={}, corners={}, steps={})",
            g.n_regions(),
            g.n_vertices(),
            g.n_corners(),
            g.n_steps()
        ))
    }
}

/// The visibility penalty of the paper (Sec. 7.1), `eta_i = nu_i * max(n_i . e_i, 0)`: `n_i` is
/// the area-weighted unit normal of the triangulation at vertex i, `e_i` the unit vector from the
/// vertex to `camera` (3 numbers), `nu_i` the visibility of the vertex (1: the camera sees it;
/// 0: something is in the way; the function casts no rays; `visibility=None` means 1
/// everywhere). Returns `array('d')`, one number per vertex in [0, nu_i], without unit: scale
/// it to the length weight before passing it as `penalty`.
#[pyfunction]
#[pyo3(signature = (positions, corner_vertices, polygon_starts, triangles, camera, *, visibility = None))]
fn visibility_penalty<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    corner_vertices: &Bound<'py, PyAny>,
    polygon_starts: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    camera: [f64; 3],
    visibility: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyAny>> {
    let buffers = MeshBuffers {
        positions: reals("positions", positions, 3, "vertex")?,
        corner_vertices: indices("corner_vertices", corner_vertices, 1, "corner")?,
        polygon_starts: indices("polygon_starts", polygon_starts, 1, "polygon")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
        ..MeshBuffers::default()
    };
    let visibility = match visibility {
        Some(v) => Some(reals("visibility", v, 1, "vertex")?),
        None => None,
    };
    let eta = py
        .detach(move || {
            let input =
                UnwrapInput::from_buffers(buffers).map_err(|e| BindError::Input(e.to_string()))?;
            crate::cutting::visibility_penalty(&input, camera, visibility.as_deref())
        })
        .map_err(bind_error)?;
    real_array(py, &eta)
}

/// The length weight `alpha_L = omega * wavelength / (2 pi)` at which ripples of the cut longer
/// than `wavelength` (on the unit-area surface) grow and shorter ones decay, for the roughening
/// rate `omega` (`CutSession.roughening_rate`).
#[pyfunction]
fn length_weight_for_ripples(roughening_rate: f64, wavelength: f64) -> f64 {
    crate::cutting::length_weight_for_ripples(roughening_rate, wavelength)
}
