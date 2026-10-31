//! The PyO3 glue of `ddg_tools.smoothest_direction_field` and `ddg_tools.adaptive_solve`
//! (`crate::intrinsic_vector`). No logic: buffers decoded, work with the GIL released, results as
//! `array.array` objects and a plain `dict`.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use super::{indices, real_array, reals, unwrap_error};
use crate::contract::ReportValue;
use crate::intrinsic_vector::{
    AdaptiveOptions, FieldOptions, TriangleArrays, adaptive, direction_field,
};

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(smoothest_direction_field, m)?)?;
    m.add_function(wrap_pyfunction!(adaptive_solve, m)?)?;
    m.add_class::<DirectionFieldResult>()?;
    m.add_class::<AdaptiveResult>()?;
    Ok(())
}

fn report_dict<'py>(
    py: Python<'py>,
    items: &[(&'static str, ReportValue)],
) -> PyResult<Bound<'py, PyDict>> {
    let report = PyDict::new(py);
    for (key, value) in items {
        match value {
            ReportValue::Bool(v) => report.set_item(key, *v)?,
            ReportValue::Int(v) => report.set_item(key, *v)?,
            ReportValue::Float(v) => report.set_item(key, *v)?,
            ReportValue::Text(v) => report.set_item(key, v)?,
            ReportValue::TextList(v) => report.set_item(key, PyList::new(py, v)?)?,
        }
    }
    Ok(report)
}

/// The result of `smoothest_direction_field`.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct DirectionFieldResult {
    /// `array('d')`, 3 per vertex: one of the field's n directions, a unit vector tangent to an
    /// input face at the vertex.
    vectors: Py<PyAny>,
    /// `array('d')`, 3 per singularity: its place on the surface.
    singularity_points: Py<PyAny>,
    /// `array('i')`, per singularity: its index in units of 1/n.
    singularity_indices: Py<PyAny>,
    /// `dict` of plain values.
    report: Py<PyDict>,
}

#[pymethods]
impl DirectionFieldResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "DirectionFieldResult(vertices={}, singularities={})",
            self.vectors.bind(py).len()? / 3,
            self.singularity_indices.bind(py).len()?
        ))
    }
}

/// The smoothest n-direction field of a triangle mesh (Knöppel et al. 2013, on an intrinsic
/// triangulation: Sharp, Soliman, Crane 2019, Sec. 5.5). Returns a DirectionFieldResult.
///
/// positions       3 reals per vertex.
/// triangles       3 vertex indices per triangle (manifold, consistently oriented).
/// n               the symmetry: 1 vectors, 2 lines, 4 crosses.
/// triangulation   "delaunay" (default: the intrinsic Delaunay triangulation), "refined"
///                 (intrinsic Delaunay refinement to min_angle_degrees) or "input".
/// laplacian       "fem" (finite elements, default) or "fd" (cotan weights and transport).
///
/// The report has status ("converged" or "iteration_limit"), iterations, eigenvalue, energy
/// (eigenvalue times area), residual, n_triangulation_vertices, index_sum (n times the Euler
/// characteristic on a closed surface), euler_characteristic, seconds.
///
/// Raises InputError for malformed arguments, MeshError for a mesh that is not a manifold
/// triangle mesh, SolveError if the eigenproblem fails.
#[pyfunction]
#[pyo3(signature = (positions, triangles, *, n = 1, triangulation = "delaunay", laplacian = "fem", min_angle_degrees = 25.0))]
fn smoothest_direction_field<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    n: u32,
    triangulation: &str,
    laplacian: &str,
    min_angle_degrees: f64,
) -> PyResult<DirectionFieldResult> {
    let mesh = TriangleArrays {
        positions: reals("positions", positions, 3, "vertex")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
    };
    let options = FieldOptions {
        n,
        triangulation: triangulation.to_owned(),
        laplacian: laplacian.to_owned(),
        min_angle_degrees,
    };
    let out = py
        .detach(move || direction_field(&mesh, &options))
        .map_err(unwrap_error)?;
    let signed: Vec<u8> = out
        .singularity_indices
        .iter()
        .flat_map(|x| x.to_ne_bytes())
        .collect();
    Ok(DirectionFieldResult {
        vectors: real_array(py, &out.vectors)?.unbind(),
        singularity_points: real_array(py, &out.singularity_points)?.unbind(),
        singularity_indices: super::new_array(py, "i", &signed)?.unbind(),
        report: report_dict(py, &out.report)?.unbind(),
    })
}

/// The result of `adaptive_solve`.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct AdaptiveResult {
    /// `array('d')`, n_components per input vertex (interleaved).
    values: Py<PyAny>,
    /// 1 (green, heat) or 2 (disk: coordinates in the unit disk).
    n_components: usize,
    /// `dict` of plain values.
    report: Py<PyDict>,
}

#[pymethods]
impl AdaptiveResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "AdaptiveResult(vertices={}, components={})",
            self.values.bind(py).len()? / self.n_components.max(1),
            self.n_components
        ))
    }
}

/// Solve a model problem with adaptive refinement of an intrinsic triangulation (Sharp,
/// Soliman, Crane 2019, Sec. 5.4). Returns an AdaptiveResult with the solution at the vertices.
///
/// problem      "green" (harmonic Green's function of `source`, mean zero), "heat" (one backward
///              Euler step of length `time` of the heat equation from a unit point source at
///              `source`; default time: the squared mean edge length) or "disk" (a harmonic map
///              of a disk-shaped mesh onto the unit disk, boundary by arc length).
/// max_vertices the vertex budget of the refinement.
/// bulk         the fraction of the error estimate the marked triangles carry (default 0.5).
/// representation "signposts" (default) or "integer" (Gillespie, Sharp, Crane 2021).
///
/// The report has status ("vertex_limit", "tolerance", "iteration_limit", "stalled"),
/// iterations, n_triangulation_vertices, estimate, oscillation, seconds.
///
/// Raises InputError for malformed arguments, MeshError for a mesh that is not a manifold
/// triangle mesh (or not a disk for "disk"), SolveError if refinement or a solve fails.
#[pyfunction]
#[pyo3(signature = (positions, triangles, problem, *, source = 0, time = None, max_vertices = 20000, bulk = 0.5, min_angle_degrees = 25.0, representation = "signposts"))]
#[allow(clippy::too_many_arguments)]
fn adaptive_solve<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    problem: &str,
    source: u32,
    time: Option<f64>,
    max_vertices: usize,
    bulk: f64,
    min_angle_degrees: f64,
    representation: &str,
) -> PyResult<AdaptiveResult> {
    let mesh = TriangleArrays {
        positions: reals("positions", positions, 3, "vertex")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
    };
    let options = AdaptiveOptions {
        problem: problem.to_owned(),
        source,
        time,
        max_vertices,
        bulk,
        min_angle_degrees,
        representation: representation.to_owned(),
    };
    let out = py
        .detach(move || adaptive(&mesh, &options))
        .map_err(unwrap_error)?;
    Ok(AdaptiveResult {
        values: real_array(py, &out.values)?.unbind(),
        n_components: out.n_components,
        report: report_dict(py, &out.report)?.unbind(),
    })
}
