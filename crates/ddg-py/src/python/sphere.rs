//! The PyO3 glue of `ddg_tools.map_to_sphere`. As in the parent module no logic lives here:
//! buffers are decoded, [`crate::sphere`] does the work with the GIL released, the result goes
//! back as an `array.array` and a plain dictionary.

use std::time::Instant;

use pyo3::prelude::*;
use pyo3::types::PyDict;

use super::{InputError, indices, real_array, reals, unwrap_error};
use crate::sphere::{SphereInput, map_to_sphere as rust_map_to_sphere};

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(map_to_sphere, m)?)?;
    m.add_class::<SphereResult>()?;
    Ok(())
}

/// The result of `map_to_sphere`.
#[pyclass(module = "ddg_tools", frozen, get_all)]
struct SphereResult {
    /// `array('d')`: 3 per vertex, a point of the unit sphere.
    positions: Py<PyAny>,
    /// `dict` of plain values (`bool`, `int`, `float`, `str`, `None`).
    report: Py<PyDict>,
}

#[pymethods]
impl SphereResult {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let get = |key: &str| -> PyResult<String> {
            Ok(match self.report.bind(py).get_item(key)? {
                Some(value) => value.to_string(),
                None => String::new(),
            })
        };
        Ok(format!(
            "SphereResult(vertices={}, pole={}, is_injective={})",
            self.positions.bind(py).len()? / 3,
            get("pole")?,
            get("is_injective")?
        ))
    }
}

/// Map a closed genus-0 triangle mesh conformally to the unit sphere. Returns a SphereResult.
///
/// Boundary First Flattening (Sawhney, Crane 2017), Sec. 7.5: one vertex star is removed, the
/// rest is flattened to a round disk, projected stereographically to the sphere and then
/// Moebius-centered (Baden, Crane, Kazhdan 2018) so that the surface's area is balanced around
/// the origin. The map is determined up to a rotation of the sphere.
///
/// positions   3 reals per vertex.
/// triangles   3 vertex indices per triangle: a connected, closed, consistently oriented
///             manifold triangle mesh of genus 0 in which every vertex is used.
/// pole        the vertex to remove (it lands on the sphere like every other vertex); None
///             (default) picks the vertex with the most regular star.
/// center      Moebius-center the result (default True). Without it the pole is at (0, 0, -1).
///
/// The result has positions (array('d'), 3 per vertex, on the unit sphere) and a report dict:
/// pole, is_injective, flipped_triangles, degenerate_triangles, covering_degree,
/// quasi_conformal_mean, quasi_conformal_max, area_distortion, log_scale_min, log_scale_max,
/// max_radius_error, centering_status ("converged", "max_iterations", "stalled",
/// "singular_jacobian", or None without centering), centering_iterations, center_of_mass
/// (float or None), seconds, seconds_total.
///
/// Raises InputError for malformed arguments and a pole that is not a vertex; MeshError for a
/// mesh that is not a valid manifold triangle mesh, has degenerate triangles, has a boundary,
/// several components or a genus other than 0; SolveError if the computation fails; TypeError
/// for an argument without the buffer protocol.
#[pyfunction]
#[pyo3(signature = (positions, triangles, *, pole = None, center = true))]
fn map_to_sphere<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    pole: Option<i64>,
    center: bool,
) -> PyResult<SphereResult> {
    let start = Instant::now();
    let pole = match pole {
        None => None,
        Some(v) => Some(u32::try_from(v).map_err(|_| {
            InputError::new_err(format!("pole: {v} is not a vertex index (must be >= 0)"))
        })?),
    };
    let input = SphereInput {
        positions: reals("positions", positions, 3, "vertex")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
        pole,
        center,
    };
    let output = py
        .detach(move || rust_map_to_sphere(&input))
        .map_err(unwrap_error)?;

    let report = PyDict::new(py);
    let r = &output.report;
    report.set_item("pole", r.pole)?;
    report.set_item("is_injective", r.is_injective)?;
    report.set_item("flipped_triangles", r.flipped_triangles)?;
    report.set_item("degenerate_triangles", r.degenerate_triangles)?;
    report.set_item("covering_degree", r.covering_degree)?;
    report.set_item("quasi_conformal_mean", r.quasi_conformal_mean)?;
    report.set_item("quasi_conformal_max", r.quasi_conformal_max)?;
    report.set_item("area_distortion", r.area_distortion)?;
    report.set_item("log_scale_min", r.log_scale_min)?;
    report.set_item("log_scale_max", r.log_scale_max)?;
    report.set_item("max_radius_error", r.max_radius_error)?;
    report.set_item("centering_status", r.centering_status)?;
    report.set_item("centering_iterations", r.centering_iterations)?;
    report.set_item("center_of_mass", r.center_of_mass)?;
    report.set_item("seconds", r.seconds)?;
    report.set_item("seconds_total", start.elapsed().as_secs_f64())?;
    Ok(SphereResult {
        positions: real_array(py, &output.positions)?.unbind(),
        report: report.unbind(),
    })
}
