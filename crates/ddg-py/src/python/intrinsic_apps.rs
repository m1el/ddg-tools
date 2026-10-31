//! The PyO3 glue of `geodesic_distance` and `steiner_tree` (applications of intrinsic
//! triangulations, [`crate::intrinsic_apps`]). As in the parent module no logic lives here.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use super::{index_array, indices, real_array, reals, unwrap_error};
use crate::intrinsic_apps::{
    HeatArguments, MeshInput, SteinerArguments, geodesic_distance as rust_geodesic_distance,
    steiner_tree as rust_steiner_tree,
};

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(geodesic_distance, m)?)?;
    m.add_function(wrap_pyfunction!(steiner_tree, m)?)?;
    m.add_class::<GeodesicDistanceResult>()?;
    m.add_class::<SteinerTreeResult>()?;
    Ok(())
}

fn mesh_input(positions: &Bound<'_, PyAny>, triangles: &Bound<'_, PyAny>) -> PyResult<MeshInput> {
    Ok(MeshInput {
        positions: reals("positions", positions, 3, "vertex")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
    })
}

/// The result of `geodesic_distance`.
#[pyclass(frozen, module = "ddg_tools")]
pub struct GeodesicDistanceResult {
    /// `array('d')`, the distance at every vertex.
    #[pyo3(get)]
    distances: Py<PyAny>,
    /// `dict`: n_triangulation_vertices, time, seconds_factor, seconds_solve.
    #[pyo3(get)]
    report: Py<PyDict>,
}

/// Geodesic distance from source vertices by the heat method on an intrinsic triangulation
/// (Sharp, Soliman, Crane 2019, Sec. 5.3; Crane, Weischedel, Wardetzky 2013). Returns a
/// GeodesicDistanceResult.
///
/// positions      3 reals per vertex.
/// triangles      3 vertex indices per triangle: a manifold, consistently oriented mesh.
/// sources        vertex indices (at least one).
/// triangulation  "refined" (intrinsic Delaunay refinement, the default), "delaunay" (intrinsic
///                Delaunay triangulation) or "input" (the mesh as it is).
/// time_factor    m in the time step t = m h^2, h the mean edge length (default 1).
/// boundary       "neumann" (default), "dirichlet" or "average": conditions of the heat step.
/// representation "signposts" (default) or "integer": the intrinsic triangulation's representation
///                (Gillespie, Sharp, Crane 2021); the distances agree to rounding.
///
/// Raises InputError for malformed arguments, MeshError for a mesh that is not a valid
/// manifold triangle mesh, SolveError when the triangulation or a solve fails.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (positions, triangles, sources, *, triangulation = "refined", time_factor = 1.0, boundary = "neumann", representation = "signposts"))]
fn geodesic_distance<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    sources: &Bound<'py, PyAny>,
    triangulation: &str,
    time_factor: f64,
    boundary: &str,
    representation: &str,
) -> PyResult<GeodesicDistanceResult> {
    let input = mesh_input(positions, triangles)?;
    let args = HeatArguments {
        sources: indices("sources", sources, 1, "source")?,
        triangulation: triangulation.to_owned(),
        time_factor,
        boundary: boundary.to_owned(),
        representation: representation.to_owned(),
    };
    let out = py
        .detach(move || rust_geodesic_distance(&input, &args))
        .map_err(unwrap_error)?;
    let report = PyDict::new(py);
    report.set_item("n_triangulation_vertices", out.n_triangulation_vertices)?;
    report.set_item("time", out.time)?;
    report.set_item("seconds_factor", out.seconds_factor)?;
    report.set_item("seconds_solve", out.seconds_solve)?;
    Ok(GeodesicDistanceResult {
        distances: real_array(py, &out.distances)?.unbind(),
        report: report.unbind(),
    })
}

/// The result of `steiner_tree`.
#[pyclass(frozen, module = "ddg_tools")]
pub struct SteinerTreeResult {
    /// `array('d')`, 3 reals per point of the polylines.
    #[pyo3(get)]
    points: Py<PyAny>,
    /// `array('i')`, where the polyline of every tree edge starts in `points`, plus the end.
    #[pyo3(get)]
    offsets: Py<PyAny>,
    /// `dict`: length, mst_length, n_steiner_points, is_disk, max_interior_defect.
    #[pyo3(get)]
    report: Py<PyDict>,
}

/// A short cut through all vertices: an approximate Steiner tree over the intrinsic Delaunay
/// triangulation with Fermat points (Sharp, Soliman, Crane 2019, Sec. 5.1). Returns a
/// SteinerTreeResult whose polylines lie on the surface.
///
/// positions       3 reals per vertex.
/// triangles       3 vertex indices per triangle.
/// intrinsic       over the intrinsic Delaunay triangulation (default True) or the input edges.
/// steiner_points  insert Fermat points (default True); False gives the minimal spanning tree.
/// representation "signposts" (default) or "integer" (Gillespie, Sharp, Crane 2021).
///
/// Raises InputError for malformed arguments, MeshError for a mesh that is not a valid
/// manifold triangle mesh, SolveError when the triangulation or a trace fails.
#[pyfunction]
#[pyo3(signature = (positions, triangles, *, intrinsic = true, steiner_points = true, representation = "signposts"))]
fn steiner_tree<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    intrinsic: bool,
    steiner_points: bool,
    representation: &str,
) -> PyResult<SteinerTreeResult> {
    let input = mesh_input(positions, triangles)?;
    let args = SteinerArguments {
        intrinsic,
        steiner_points,
        representation: representation.to_owned(),
    };
    let out = py
        .detach(move || rust_steiner_tree(&input, &args))
        .map_err(unwrap_error)?;
    let report = PyDict::new(py);
    report.set_item("length", out.length)?;
    report.set_item("mst_length", out.mst_length)?;
    report.set_item("n_steiner_points", out.n_steiner_points)?;
    report.set_item("is_disk", out.is_disk)?;
    report.set_item("max_interior_defect", out.max_interior_defect)?;
    Ok(SteinerTreeResult {
        points: real_array(py, &out.points)?.unbind(),
        offsets: index_array(py, "offsets", out.offsets.iter().copied())?.unbind(),
        report: report.unbind(),
    })
}
