//! The PyO3 glue of `ddg_tools.IntrinsicTriangulation` ([`crate::integer_coords`]). As in the
//! parent module no logic lives here.

use std::sync::{Mutex, MutexGuard};

use pyo3::prelude::*;
use pyo3::types::PyDict;

use super::{SolveError, index_array, indices, real_array, reals, unwrap_error};
use crate::integer_coords::{
    CertifyOutput, Counts, FlipOutput, Places, PolylinesOutput, RefineArguments, RefineOutput,
    Session, SubdivisionOutput,
};
use crate::intrinsic_apps::MeshInput;

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<IntrinsicTriangulation>()?;
    Ok(())
}

fn places<'py>(py: Python<'py>, p: &Places) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let kinds: Vec<u32> = p.kind.iter().map(|&k| k as u32).collect();
    d.set_item("kind", index_array(py, "kind", kinds.into_iter())?)?;
    d.set_item("id", index_array(py, "id", p.id.iter().copied())?)?;
    d.set_item("params", real_array(py, &p.params)?)?;
    Ok(d)
}

fn counts<'py>(py: Python<'py>, c: &Counts) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("n_vertices", c.n_vertices)?;
    d.set_item("n_edges", c.n_edges)?;
    d.set_item("n_faces", c.n_faces)?;
    d.set_item("n_input_vertices", c.n_input_vertices)?;
    d.set_item("n_input_faces", c.n_input_faces)?;
    d.set_item("min_angle", c.min_angle)?;
    d.set_item("max_angle", c.max_angle)?;
    d.set_item("n_crossings", c.n_crossings)?;
    Ok(d)
}

/// An intrinsic triangulation of a mesh: integer coordinates (Gillespie, Sharp, Crane 2021, the
/// default) or signposts (Sharp, Soliman, Crane 2019).
///
/// IntrinsicTriangulation(positions, triangles, representation="integer", mollify=None)
///
/// positions       3 reals per vertex.
/// triangles       3 vertex indices per triangle: a manifold, consistently oriented mesh.
/// representation  "integer" (default) or "signposts".
/// mollify         relative margin of intrinsic mollification (integer only; the paper uses
///                 1e-4); None: no mollification.
///
/// The input mesh is never changed. Vertex ids 0..n_input_vertices are those of the input.
///
/// Raises InputError for malformed arguments, MeshError for a mesh that is not a valid
/// manifold triangle mesh, SolveError when an operation fails.
#[pyclass(module = "ddg_tools", frozen)]
pub struct IntrinsicTriangulation {
    inner: Mutex<Session>,
}

impl IntrinsicTriangulation {
    fn lock(&self) -> PyResult<MutexGuard<'_, Session>> {
        self.inner
            .lock()
            .map_err(|_| SolveError::new_err("the triangulation is poisoned by an earlier panic"))
    }
}

#[pymethods]
impl IntrinsicTriangulation {
    #[new]
    #[pyo3(signature = (positions, triangles, representation = "integer", mollify = None))]
    fn new(
        py: Python<'_>,
        positions: &Bound<'_, PyAny>,
        triangles: &Bound<'_, PyAny>,
        representation: &str,
        mollify: Option<f64>,
    ) -> PyResult<Self> {
        let input = MeshInput {
            positions: reals("positions", positions, 3, "vertex")?,
            triangles: indices("triangles", triangles, 3, "triangle")?,
        };
        let representation = representation.to_owned();
        let session = py
            .detach(move || Session::new(&input, &representation, mollify))
            .map_err(unwrap_error)?;
        Ok(IntrinsicTriangulation {
            inner: Mutex::new(session),
        })
    }

    /// "integer" or "signposts".
    #[getter]
    fn representation(&self) -> PyResult<&'static str> {
        Ok(self.lock()?.representation())
    }

    /// dict: n_vertices, n_edges, n_faces, n_input_vertices, n_input_faces, min_angle,
    /// max_angle (radians), n_crossings (integer coordinates; 0 for signposts).
    fn counts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let c = self.lock()?.counts();
        counts(py, &c)
    }

    /// Flips edges until every edge is Delaunay (paper Sec. 4.1). Returns a dict: status
    /// ("delaunay", "flip_limit", "not_flippable"), n_flips, n_not_delaunay.
    fn flip_to_delaunay<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let mut guard = self.lock()?;
        let s: &mut Session = &mut guard;
        let out: FlipOutput = py.detach(|| s.flip_to_delaunay()).map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item("status", out.status)?;
        d.set_item("n_flips", out.n_flips)?;
        d.set_item("n_not_delaunay", out.n_not_delaunay)?;
        Ok(d)
    }

    /// Delaunay refinement to an angle bound. rules: "chew" (default; the integer paper's
    /// rules) or "ruppert" (the signpost paper's). Returns a dict: status ("done",
    /// "insertion_limit"), n_inserted, n_face_insertions, n_edge_splits, n_flips, n_removals,
    /// n_removals_refused, min_angle_valid (radians), n_below_bound, certified.
    #[pyo3(signature = (min_angle_degrees = 25.0, rules = "chew", max_insertions = None, max_circumradius = None))]
    fn refine<'py>(
        &self,
        py: Python<'py>,
        min_angle_degrees: f64,
        rules: &str,
        max_insertions: Option<usize>,
        max_circumradius: Option<f64>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let args = RefineArguments {
            min_angle_degrees,
            rules: rules.to_owned(),
            max_insertions,
            max_circumradius,
        };
        let mut guard = self.lock()?;
        let s: &mut Session = &mut guard;
        let out: RefineOutput = py.detach(|| s.refine(&args)).map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item("status", out.status)?;
        d.set_item("n_inserted", out.n_inserted)?;
        d.set_item("n_face_insertions", out.n_face_insertions)?;
        d.set_item("n_edge_splits", out.n_edge_splits)?;
        d.set_item("n_flips", out.n_flips)?;
        d.set_item("n_removals", out.n_removals)?;
        d.set_item("n_removals_refused", out.n_removals_refused)?;
        d.set_item("min_angle_valid", out.min_angle_valid)?;
        d.set_item("n_below_bound", out.n_below_bound)?;
        d.set_item("certified", out.certified)?;
        Ok(d)
    }

    /// (triangles, edge_lengths): array('i') 3 vertex ids per face of the triangulation, and
    /// array('d') one length per edge.
    fn triangles<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
        let (t, l) = self.lock()?.triangles();
        Ok((
            index_array(py, "triangles", t.into_iter())?,
            real_array(py, &l)?,
        ))
    }

    /// dict of the places on the input mesh of the vertices of the triangulation: kind
    /// (array('i'): 0 vertex, 1 edge, 2 face), id, params (array('d'), 3 per vertex) and
    /// positions (array('d'), 3 per vertex, in space).
    fn vertex_places<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let (p, xyz) = self.lock()?.vertex_places();
        let d = places(py, &p)?;
        d.set_item("positions", real_array(py, &xyz)?)?;
        Ok(d)
    }

    /// The common subdivision of the input mesh and the triangulation, as a dict:
    /// positions (array('d'), 3 per vertex), on_input and on_triangulation (dicts as in
    /// vertex_places: kind, id, params), polygon_starts (array('i'), one more than polygons),
    /// corners (array('i')), input_face and intrinsic_face (array('i'), per polygon),
    /// n_intrinsic_vertices (the vertices of the triangulation come first).
    fn common_subdivision<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let s = self.lock()?;
        let o: SubdivisionOutput = py.detach(|| s.subdivision()).map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item("positions", real_array(py, &o.positions)?)?;
        d.set_item("on_input", places(py, &o.on_input)?)?;
        d.set_item("on_triangulation", places(py, &o.on_triangulation)?)?;
        d.set_item(
            "polygon_starts",
            index_array(py, "polygon_starts", o.polygon_starts.into_iter())?,
        )?;
        d.set_item(
            "corners",
            index_array(py, "corners", o.corners.into_iter())?,
        )?;
        d.set_item(
            "input_face",
            index_array(py, "input_face", o.input_face.into_iter())?,
        )?;
        d.set_item(
            "intrinsic_face",
            index_array(py, "intrinsic_face", o.intrinsic_face.into_iter())?,
        )?;
        d.set_item("n_intrinsic_vertices", o.n_intrinsic_vertices)?;
        Ok(d)
    }

    /// The intrinsic edges as polylines on the input mesh (edges: ids, default all), as a dict:
    /// edges (array('i')), points (array('d'), 3 per point, on the input surface), offsets
    /// (array('i'), where each polyline starts in points, one more entry at the end),
    /// along_input_edge (list: the input edge a polyline lies along or -1), n_crossings (list),
    /// max_relative_mismatch.
    #[pyo3(signature = (edges = None))]
    fn edge_polylines<'py>(
        &self,
        py: Python<'py>,
        edges: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let ids = edges.map(|e| indices("edges", e, 1, "edge")).transpose()?;
        let s = self.lock()?;
        let o: PolylinesOutput = py
            .detach(|| s.edge_polylines(ids.as_deref()))
            .map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item("edges", index_array(py, "edges", o.edges.into_iter())?)?;
        d.set_item("points", real_array(py, &o.points)?)?;
        d.set_item(
            "offsets",
            index_array(py, "offsets", o.offsets.into_iter())?,
        )?;
        d.set_item("along_input_edge", o.along_input_edge)?;
        d.set_item("n_crossings", o.n_crossings)?;
        d.set_item("max_relative_mismatch", o.max_relative_mismatch)?;
        Ok(d)
    }

    /// The certificate of the integer coordinates (SPEC 2.4): a dict n_faces, n_along,
    /// n_crossings, n_input_edges, n_inserted_on_edges, max_crossings_per_input_edge. Raises
    /// InputError for signposts and SolveError naming the item that fails.
    fn certify<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let s = self.lock()?;
        let c: CertifyOutput = py.detach(|| s.certify()).map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item("n_faces", c.n_faces)?;
        d.set_item("n_along", c.n_along)?;
        d.set_item("n_crossings", c.n_crossings)?;
        d.set_item("n_input_edges", c.n_input_edges)?;
        d.set_item("n_inserted_on_edges", c.n_inserted_on_edges)?;
        d.set_item(
            "max_crossings_per_input_edge",
            c.max_crossings_per_input_edge,
        )?;
        Ok(d)
    }

    /// Transfers a function on the vertices of the triangulation (one real per vertex) to the
    /// input vertices by L2 projection (Eq. 13-14); integer coordinates only.
    fn to_input<'py>(
        &self,
        py: Python<'py>,
        values: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let values = reals("values", values, 1, "vertex")?;
        let mut guard = self.lock()?;
        let s: &mut Session = &mut guard;
        let out = py.detach(|| s.to_input(&values)).map_err(unwrap_error)?;
        real_array(py, &out)
    }

    /// Transfers a function on the input vertices to the vertices of the triangulation;
    /// integer coordinates only.
    fn to_intrinsic<'py>(
        &self,
        py: Python<'py>,
        values: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let values = reals("values", values, 1, "vertex")?;
        let mut guard = self.lock()?;
        let s: &mut Session = &mut guard;
        let out = py
            .detach(|| s.to_intrinsic(&values))
            .map_err(unwrap_error)?;
        real_array(py, &out)
    }
}
