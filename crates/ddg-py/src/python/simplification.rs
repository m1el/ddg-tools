//! The PyO3 glue of `simplify`, `build_hierarchy` and `qem_simplify` ([`crate::simplification`]).
//! As in the parent module no logic lives here: buffers in, `array.array` objects and dicts out.

use std::sync::{Mutex, MutexGuard};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict, PyString, PyTuple};

use ddg_intrinsic_simplification::multigrid::SolveReport;
use ddg_intrinsic_simplification::refine::RefinementReport;

use super::{SolveError, index_array, indices, new_array, real_array, reals, unwrap_error};
use crate::intrinsic_apps::MeshInput;
use crate::simplification::{
    AnisotropyArguments, AnisotropyDirection, Coo, HierarchyArguments, HierarchyHandle,
    MultigridArguments, MultigridHandle, QemArguments, QemOutput,
    Simplification as RustSimplification, SimplifyArguments, qem as rust_qem,
};

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(simplify, m)?)?;
    m.add_function(wrap_pyfunction!(build_hierarchy, m)?)?;
    m.add_function(wrap_pyfunction!(qem_simplify, m)?)?;
    m.add_class::<Simplification>()?;
    m.add_class::<Hierarchy>()?;
    m.add_class::<MultigridSolver>()?;
    m.add_class::<QemResult>()?;
    Ok(())
}

fn mesh_input(positions: &Bound<'_, PyAny>, triangles: &Bound<'_, PyAny>) -> PyResult<MeshInput> {
    Ok(MeshInput {
        positions: reals("positions", positions, 3, "vertex")?,
        triangles: indices("triangles", triangles, 3, "triangle")?,
    })
}

fn opt_indices(name: &'static str, v: Option<&Bound<'_, PyAny>>) -> PyResult<Vec<u32>> {
    Ok(match v {
        Some(v) if !v.is_none() => indices(name, v, 1, "vertex")?,
        _ => Vec::new(),
    })
}

fn opt_reals(name: &'static str, v: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec<f64>>> {
    Ok(match v {
        Some(v) if !v.is_none() => Some(reals(name, v, 1, "value")?),
        _ => None,
    })
}

/// `array('i')` of signed values (-1 for none).
fn signed_array<'py>(
    py: Python<'py>,
    name: &'static str,
    values: &[i64],
) -> PyResult<Bound<'py, PyAny>> {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for (p, &v) in values.iter().enumerate() {
        let v = i32::try_from(v).map_err(|_| {
            SolveError::new_err(format!(
                "internal error: {name}[{p}] = {v} does not fit int32"
            ))
        })?;
        bytes.extend_from_slice(&v.to_ne_bytes());
    }
    new_array(py, "i", &bytes)
}

fn coo_tuple<'py>(py: Python<'py>, c: &Coo) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(
        py,
        [
            index_array(py, "rows", c.rows.iter().copied())?,
            index_array(py, "cols", c.cols.iter().copied())?,
            real_array(py, &c.values)?,
        ],
    )
}

/// `scipy.sparse.coo_matrix((values, (rows, cols)), shape)`; ImportError if scipy is missing.
fn to_scipy<'py>(
    py: Python<'py>,
    values: Bound<'py, PyAny>,
    rows: Bound<'py, PyAny>,
    cols: Bound<'py, PyAny>,
    shape: (usize, usize),
) -> PyResult<Bound<'py, PyAny>> {
    let sparse = py.import("scipy.sparse")?;
    let np = py.import("numpy")?;
    let rows = np.getattr("asarray")?.call1((rows,))?;
    let cols = np.getattr("asarray")?.call1((cols,))?;
    let values = np.getattr("asarray")?.call1((values,))?;
    let coo = sparse
        .getattr("coo_matrix")?
        .call1(((values, (rows, cols)), shape))?;
    coo.call_method0("tocsr")
}

fn solve_report<'py>(py: Python<'py>, r: &SolveReport) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("converged", r.converged)?;
    d.set_item("cycles", r.cycles)?;
    d.set_item("residuals", r.residuals.clone())?;
    d.set_item("rhs_norm", r.rhs_norm)?;
    d.set_item("relative_residual", r.relative_residual)?;
    d.set_item("convergence_factor", r.convergence_factor)?;
    d.set_item("level_sizes", r.level_sizes.clone())?;
    d.set_item("solve_seconds", r.solve_seconds)?;
    Ok(d)
}

fn refinement_report<'py>(py: Python<'py>, r: &RefinementReport) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item(
        "representation",
        format!("{:?}", r.representation).to_lowercase(),
    )?;
    d.set_item("n_input_vertices", r.n_input_vertices)?;
    d.set_item("n_vertices", r.n_vertices)?;
    d.set_item("n_added", r.n_added)?;
    d.set_item("min_angle", r.min_angle)?;
    d.set_item("max_angle", r.max_angle)?;
    d.set_item("closed", r.closed)?;
    d.set_item("guarantee_applies", r.guarantee_applies)?;
    d.set_item("n_located", r.n_located)?;
    d.set_item("n_clamped", r.n_clamped)?;
    d.set_item("seconds_refine", r.seconds_refine)?;
    d.set_item("seconds_locate", r.seconds_locate)?;
    match &r.repair {
        Some(p) => {
            let e = PyDict::new(py);
            e.set_item("defects_before", p.defects_before)?;
            e.set_item("defects_after", p.defects_after)?;
            e.set_item("n_flips", p.n_flips)?;
            d.set_item("repair", e)?;
        }
        None => d.set_item("repair", py.None())?,
    }
    Ok(d)
}

fn string_map<'py>(
    py: Python<'py>,
    m: &std::collections::BTreeMap<&'static str, usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    for (k, v) in m {
        d.set_item(*k, *v)?;
    }
    Ok(d)
}

/// A coarse mesh of a surface: the result of `simplify`.
///
/// Fine vertices are the input's vertices; coarse vertices keep their input ids (`vertices`),
/// except vertices that refinement inserted, whose ids start at the input vertex count
/// (`n_input_vertices`). Tracked points are the input vertices (ids 0..n_input_vertices-1)
/// followed by the extra points of `points=`.
///
/// Attributes (arrays are `array.array`; 3-wide arrays are flat):
///
/// n_input_vertices, n_vertices, n_faces   counts
/// vertices            fine id of every coarse vertex (int32)
/// faces               3 fine ids per coarse face
/// faces_compact       3 coarse vertex ids (0..n_vertices-1) per face
/// lengths             3 reals per face: the side from corner k to corner k+1 (a Δ-complex face may
///                     repeat a vertex)
/// point_face          coarse face of every tracked point (a coarse vertex: its lowest incident
///                     face), -1 for none
/// point_bary          3 barycentric coordinates per tracked point in that face's corners
/// point_vertex        compacted coarse vertex a tracked point is, -1 inside a face
/// face_colours, n_colours   greedy colouring of the coarse faces (faces that share a vertex
///                     differ); point_colours: colour of every tracked point's coarse face
/// report              dict (see below)
///
/// report: stop ("targetreached", "maxcost", "nothingremovable", "queueempty"), n_vertices_before,
/// n_vertices_after, n_interior, n_boundary, n_ear, n_initial_delaunay_flips, n_delaunay_flips,
/// n_fallback_flips, n_degree_flips, n_pops, n_stale_pops, n_evaluations, n_clamped,
/// flatten_iterations_max, flatten_iterations_mean, skips, commit_failures, infinite_evaluations
/// (dicts reason -> count), components (list of (before, after)), keys (dict), timings (dict of
/// seconds), seconds (the whole call), anisotropy (dict or None), pre_refinement and
/// post_refinement (dicts or None).
#[pyclass(module = "ddg_tools", frozen)]
pub struct Simplification {
    inner: Mutex<RustSimplification>,
}

impl Simplification {
    fn lock(&self) -> PyResult<MutexGuard<'_, RustSimplification>> {
        self.inner
            .lock()
            .map_err(|_| SolveError::new_err("the simplification is poisoned by an earlier panic"))
    }
}

#[pymethods]
impl Simplification {
    /// Vertices of the input mesh.
    #[getter]
    fn n_input_vertices(&self) -> PyResult<usize> {
        Ok(self.lock()?.n_input)
    }

    /// Vertices of the coarse mesh.
    #[getter]
    fn n_vertices(&self) -> PyResult<usize> {
        Ok(self.lock()?.output.n_vertices())
    }

    /// Faces of the coarse mesh.
    #[getter]
    fn n_faces(&self) -> PyResult<usize> {
        Ok(self.lock()?.output.n_faces())
    }

    /// array('i'): the fine id of every coarse vertex.
    #[getter]
    fn vertices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let s = self.lock()?;
        index_array(py, "vertices", s.output.vertices.iter().copied())
    }

    /// array('i'): 3 fine vertex ids per coarse face.
    #[getter]
    fn faces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let v = self.lock()?.faces_fine();
        index_array(py, "faces", v.into_iter())
    }

    /// array('i'): 3 compacted vertex ids per coarse face.
    #[getter]
    fn faces_compact<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let v = self.lock()?.faces_compact();
        index_array(py, "faces_compact", v.into_iter())
    }

    /// array('d'): 3 side lengths per face (corner k to corner k+1).
    #[getter]
    fn lengths<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let v = self.lock()?.lengths();
        real_array(py, &v)
    }

    /// array('i'): the coarse face of every tracked point, -1 for none.
    #[getter]
    fn point_face<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let l = self.lock()?.point_locations();
        signed_array(py, "point_face", &l.face)
    }

    /// array('d'): 3 barycentric coordinates per tracked point.
    #[getter]
    fn point_bary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let l = self.lock()?.point_locations();
        real_array(py, &l.bary)
    }

    /// array('i'): the compacted coarse vertex a tracked point is, -1 inside a face.
    #[getter]
    fn point_vertex<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let l = self.lock()?.point_locations();
        signed_array(py, "point_vertex", &l.vertex)
    }

    /// array('i'): the colour of every coarse face.
    #[getter]
    fn face_colours<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let (c, _) = self.lock()?.face_colours();
        index_array(py, "face_colours", c.into_iter())
    }

    /// The number of colours used.
    #[getter]
    fn n_colours(&self) -> PyResult<u32> {
        Ok(self.lock()?.face_colours().1)
    }

    /// array('i'): the colour of the coarse face of every tracked point, -1 for none.
    #[getter]
    fn point_colours<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let c = self.lock()?.point_colours();
        signed_array(py, "point_colours", &c)
    }

    /// The report: dict (see the class docstring).
    #[getter]
    fn report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let s = self.lock()?;
        let r = &s.report;
        let d = PyDict::new(py);
        d.set_item("stop", format!("{:?}", r.stop).to_lowercase())?;
        d.set_item("n_vertices_before", r.n_vertices_before)?;
        d.set_item("n_vertices_after", r.n_vertices_after)?;
        d.set_item("n_initial_delaunay_flips", r.n_initial_delaunay_flips)?;
        d.set_item("n_interior", r.n_interior)?;
        d.set_item("n_boundary", r.n_boundary)?;
        d.set_item("n_ear", r.n_ear)?;
        d.set_item("n_delaunay_flips", r.n_delaunay_flips)?;
        d.set_item("n_delaunay_refused", r.n_delaunay_refused)?;
        d.set_item("n_fallback_flips", r.n_fallback_flips)?;
        d.set_item("n_degree_flips", r.n_degree_flips)?;
        d.set_item("n_pops", r.n_pops)?;
        d.set_item("n_stale_pops", r.n_stale_pops)?;
        d.set_item("n_evaluations", r.n_evaluations)?;
        d.set_item("n_clamped", r.n_clamped)?;
        d.set_item("flatten_iterations_max", r.flatten_iterations_max)?;
        d.set_item("flatten_iterations_mean", r.flatten_iterations_mean)?;
        d.set_item("skips", string_map(py, &r.skips)?)?;
        d.set_item("commit_failures", string_map(py, &r.commit_failures)?)?;
        d.set_item(
            "infinite_evaluations",
            string_map(py, &r.infinite_evaluations)?,
        )?;
        d.set_item(
            "components",
            r.components
                .iter()
                .map(|c| (c.before, c.after))
                .collect::<Vec<_>>(),
        )?;
        let k = PyDict::new(py);
        k.set_item("n_clean", r.keys.n_clean)?;
        k.set_item("n_clean_mismatch", r.keys.n_clean_mismatch)?;
        k.set_item("n_differ", r.keys.n_differ)?;
        k.set_item("max_relative_difference", r.keys.max_relative_difference)?;
        d.set_item("keys", k)?;
        let t = PyDict::new(py);
        t.set_item("initial_delaunay", r.timings.initial_delaunay)?;
        t.set_item("channels", r.timings.channels)?;
        t.set_item("initial_costs", r.timings.initial_costs)?;
        t.set_item("removals", r.timings.removals)?;
        t.set_item("requeue", r.timings.requeue)?;
        t.set_item("total", r.timings.total)?;
        d.set_item("timings", t)?;
        d.set_item("seconds", s.seconds)?;
        match &s.anisotropy {
            Some(a) => {
                let e = PyDict::new(py);
                e.set_item("tau_requested", a.tau_requested)?;
                e.set_item("tau_used", a.tau_used)?;
                e.set_item("clamped", a.clamped)?;
                e.set_item("tau_max", a.tau_max)?;
                e.set_item("min_factor", a.min_factor)?;
                e.set_item("max_factor", a.max_factor)?;
                e.set_item("mean_factor", a.mean_factor)?;
                d.set_item("anisotropy", e)?;
            }
            None => d.set_item("anisotropy", py.None())?,
        }
        for (key, r) in [("pre_refinement", &s.pre), ("post_refinement", &s.post)] {
            match r {
                Some(r) => d.set_item(key, refinement_report(py, r)?)?,
                None => d.set_item(key, py.None())?,
            }
        }
        Ok(d)
    }

    /// P as COO triplets (rows, cols, values): fine vertices by coarse vertices (SPEC 7.2); rows
    /// sum to one.
    #[pyo3(name = "prolongation")]
    fn prolongation<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let s = self.lock()?;
        let c = s.prolongation().map_err(unwrap_error)?;
        coo_tuple(py, &c)
    }

    /// P as a scipy.sparse CSR matrix (ImportError without scipy and numpy).
    fn prolongation_scipy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let c = self.lock()?.prolongation().map_err(unwrap_error)?;
        to_scipy(
            py,
            real_array(py, &c.values)?,
            index_array(py, "rows", c.rows.iter().copied())?,
            index_array(py, "cols", c.cols.iter().copied())?,
            (c.nrows, c.ncols),
        )
    }

    /// The same as `prolongation` for every tracked point (the `points=` included), rows in
    /// tracked-point order.
    fn point_prolongation<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let c = self.lock()?.point_prolongation().map_err(unwrap_error)?;
        coo_tuple(py, &c)
    }

    /// P^vec for n-direction fields (n=1 vector fields, 2 line fields, 4 cross fields) as a
    /// dict: rows, cols (array('i')), re, im (array('d')), shape and report (how the rows were
    /// found). Fine and coarse values are in the normalized frames of the input mesh.
    /// Not available with refine_before, refine_after or anisotropy.
    #[pyo3(signature = (n = 1))]
    fn vector_prolongation<'py>(&self, py: Python<'py>, n: u32) -> PyResult<Bound<'py, PyDict>> {
        let s = self.lock()?;
        let c = py
            .detach(|| s.vector_prolongation(n))
            .map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item("rows", index_array(py, "rows", c.rows.iter().copied())?)?;
        d.set_item("cols", index_array(py, "cols", c.cols.iter().copied())?)?;
        d.set_item("re", real_array(py, &c.re)?)?;
        d.set_item("im", real_array(py, &c.im)?)?;
        d.set_item("shape", (c.nrows, c.ncols))?;
        let r = &c.report;
        let e = PyDict::new(py);
        e.set_item("n_same_face", r.n_same_face)?;
        e.set_item("n_adjacent_face", r.n_adjacent_face)?;
        e.set_item("n_kept_from_images", r.n_kept_from_images)?;
        e.set_item("n_averaged", r.n_averaged)?;
        e.set_item("n_kept_averaged", r.n_kept_averaged)?;
        e.set_item("n_default", r.n_default)?;
        e.set_item("n_cancelled", r.n_cancelled)?;
        d.set_item("report", e)?;
        Ok(d)
    }

    /// P^vec as a complex scipy.sparse CSR matrix (ImportError without scipy and numpy).
    #[pyo3(signature = (n = 1))]
    fn vector_prolongation_scipy<'py>(
        &self,
        py: Python<'py>,
        n: u32,
    ) -> PyResult<Bound<'py, PyAny>> {
        let s = self.lock()?;
        let c = py
            .detach(|| s.vector_prolongation(n))
            .map_err(unwrap_error)?;
        let np = py.import("numpy")?;
        let re = np.getattr("asarray")?.call1((real_array(py, &c.re)?,))?;
        let im = np.getattr("asarray")?.call1((real_array(py, &c.im)?,))?;
        let z = re.add(im.mul(PyComplex::from_doubles(py, 0.0, 1.0))?)?;
        to_scipy(
            py,
            z,
            index_array(py, "rows", c.rows.iter().copied())?,
            index_array(py, "cols", c.cols.iter().copied())?,
            (c.nrows, c.ncols),
        )
    }

    /// Poisson's equation on the coarse mesh, prolongated (paper Fig. 21): solves L x = M f
    /// with f sampled from `f` (one real per input vertex) at the coarse vertices. `fixed`:
    /// input vertices that are coarse vertices (keep them with simplify(fixed=...)) with
    /// Dirichlet `values`; without them the mean-zero problem on a closed surface. Returns a
    /// dict: coarse (one real per coarse vertex), fine (one per input vertex, P x), report.
    #[pyo3(signature = (f, fixed = None, values = None))]
    fn poisson<'py>(
        &self,
        py: Python<'py>,
        f: &Bound<'py, PyAny>,
        fixed: Option<&Bound<'py, PyAny>>,
        values: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let f = reals("f", f, 1, "vertex")?;
        let fixed = opt_indices("fixed", fixed)?;
        let values = opt_reals("values", values)?.unwrap_or_else(|| vec![0.0; fixed.len()]);
        let mut guard = self.lock()?;
        let s: &mut RustSimplification = &mut guard;
        let out = py
            .detach(|| s.poisson(&f, &fixed, &values))
            .map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item("coarse", real_array(py, &out.coarse)?)?;
        d.set_item("fine", real_array(py, &out.fine)?)?;
        let r = PyDict::new(py);
        r.set_item("n_free", out.report.n_free)?;
        r.set_item("n_fixed", out.report.n_fixed)?;
        r.set_item("rhs_shift", out.report.rhs_shift)?;
        r.set_item("residual_max", out.report.residual_max)?;
        d.set_item("report", r)?;
        Ok(d)
    }

    /// Geodesic distance from tracked point `source` (an input vertex id, or an extra point
    /// n_input_vertices + k) computed on the coarse mesh and prolongated to the input
    /// vertices. method: "exact" (default; exact polyhedral distance on the Δ-complex) or "heat"
    /// (heat method; time_factor scales the time step). Returns a dict: coarse, fine.
    #[pyo3(signature = (source, method = "exact", time_factor = 1.0))]
    fn geodesic<'py>(
        &self,
        py: Python<'py>,
        source: usize,
        method: &str,
        time_factor: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let mut guard = self.lock()?;
        let s: &mut RustSimplification = &mut guard;
        let (coarse, fine) = py
            .detach(|| s.geodesic(source, method, time_factor))
            .map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item("coarse", real_array(py, &coarse)?)?;
        d.set_item("fine", real_array(py, &fine)?)?;
        Ok(d)
    }

    /// Geodesic Voronoi labels from tracked points `sources` (paper Fig. 25), exact distance on
    /// the coarse mesh, one run per source in parallel (threads, at most 16). Returns a dict:
    /// fine_labels (index into sources, per input vertex), coarse_labels, fine_distance.
    #[pyo3(signature = (sources, threads = 4))]
    fn voronoi<'py>(
        &self,
        py: Python<'py>,
        sources: &Bound<'py, PyAny>,
        threads: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        let sources: Vec<usize> = indices("sources", sources, 1, "source")?
            .into_iter()
            .map(|s| s as usize)
            .collect();
        let mut guard = self.lock()?;
        let s: &mut RustSimplification = &mut guard;
        let v = py
            .detach(|| s.voronoi(&sources, threads))
            .map_err(unwrap_error)?;
        let d = PyDict::new(py);
        d.set_item(
            "fine_labels",
            index_array(py, "fine_labels", v.fine_labels.into_iter())?,
        )?;
        d.set_item(
            "coarse_labels",
            index_array(py, "coarse_labels", v.coarse_labels.into_iter())?,
        )?;
        d.set_item("fine_distance", real_array(py, &v.fine_distance)?)?;
        Ok(d)
    }

    fn __repr__(&self) -> PyResult<String> {
        let s = self.lock()?;
        Ok(format!(
            "<ddg_tools.Simplification {} -> {} vertices, {} faces>",
            s.n_input,
            s.output.n_vertices(),
            s.output.n_faces()
        ))
    }
}

/// Simplifies a triangle mesh by intrinsic vertex removal with the intrinsic error metric
/// (Liu, Gillespie, Chislett, Sharp, Jacobson, Crane, ACM TOG 42(4), 2023, Alg. 1). The input is
/// never changed; the coarse mesh is an intrinsic triangulation (a Δ-complex) of the same
/// surface: its vertices are a subset of the input's, its faces are given by edge lengths, and
/// every input vertex and extra point is mapped onto it. Returns a Simplification.
///
/// positions, triangles  3 reals per vertex, 3 vertex indices per triangle: a manifold,
///                consistently oriented mesh (boundaries allowed).
/// target         stop at this many vertices, or ...
/// ratio          ... at round(ratio * input vertices) (exactly one of them).
/// weights        w in [0, 1] (default 0): the metric's mass is curvature (1 - w) plus area (w).
/// memory         memory-based metric (default True) or the memoryless one (False).
/// fixed          input vertices that stay (indices).
/// masses         one real per vertex: an extra mass channel (the sign splits it in two
///                channels), weight mass_weight (default 1).
/// points         extra tracked points: a tuple (faces, bary), faces one input face index per
///                point and bary 3 barycentric coordinates each. They are tracked points
///                n_input_vertices + k (the Blender view uses the face centroids).
/// refine_before  Delaunay-refine the input first (paper Fig. 12); inserted vertices have ids
///                >= n_input_vertices and are removable.
/// refine_after   Delaunay-refine the coarse mesh at the end (Fig. 13) to a 25 degree minimum
///                angle; repair_simplicial (default True) first flips a Δ-complex into a
///                simplicial mesh where the refinement needs one.
/// anisotropy     strength tau in [0, 1] of anisotropic coarsening (Sec. 8.6.2), None: off;
///                anisotropy_direction "max" (default), "min" (principal curvature direction
///                lines) or 3 reals per vertex; anisotropy_clamp lowers tau to the largest valid
///                value instead of failing.
/// representation "integer" (default) or "signposts": the refinements' representation.
/// projective_sign "negative" (default) or "positive": the sign of the exponent in how tracked
///                points move through a flattening (SPEC S25): "negative" (e^-u) is the
///                circumcircle-preserving map of Springborn et al. 2008 that the paper cites,
///                "positive" (e^+u) is Eq. 10 as printed (the reference's). Same coarse mesh;
///                "negative" distorts the map much less.
///
/// Raises InputError for malformed arguments, MeshError for a mesh that is not a valid manifold
/// triangle mesh, SolveError when a refinement or the simplification fails.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (positions, triangles, *, target = None, ratio = None, weights = 0.0, memory = true, fixed = None, masses = None, mass_weight = 1.0, points = None, refine_before = false, refine_after = false, repair_simplicial = true, anisotropy = None, anisotropy_direction = None, anisotropy_clamp = false, representation = "integer", projective_sign = "negative"))]
fn simplify<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    target: Option<usize>,
    ratio: Option<f64>,
    weights: f64,
    memory: bool,
    fixed: Option<&Bound<'py, PyAny>>,
    masses: Option<&Bound<'py, PyAny>>,
    mass_weight: f64,
    points: Option<&Bound<'py, PyAny>>,
    refine_before: bool,
    refine_after: bool,
    repair_simplicial: bool,
    anisotropy: Option<f64>,
    anisotropy_direction: Option<&Bound<'py, PyAny>>,
    anisotropy_clamp: bool,
    representation: &str,
    projective_sign: &str,
) -> PyResult<Simplification> {
    let input = mesh_input(positions, triangles)?;
    let (point_faces, point_bary) = match points {
        Some(p) if !p.is_none() => {
            let t = p.cast::<PyTuple>().map_err(|_| {
                PyTypeError::new_err("points: expected a tuple (faces, bary) of buffers")
            })?;
            if t.len() != 2 {
                return Err(PyTypeError::new_err(
                    "points: expected a tuple (faces, bary)",
                ));
            }
            (
                indices("points (faces)", &t.get_item(0)?, 1, "point")?,
                reals("points (bary)", &t.get_item(1)?, 3, "point")?,
            )
        }
        _ => (Vec::new(), Vec::new()),
    };
    let direction = match anisotropy_direction {
        None => AnisotropyDirection::Max,
        Some(d) if d.is_none() => AnisotropyDirection::Max,
        Some(d) => {
            if let Ok(s) = d.cast::<PyString>() {
                match s.to_str()? {
                    "max" => AnisotropyDirection::Max,
                    "min" => AnisotropyDirection::Min,
                    other => {
                        return Err(super::InputError::new_err(format!(
                            "anisotropy_direction = {other:?}: must be \"max\", \"min\" or 3 reals per vertex"
                        )));
                    }
                }
            } else {
                AnisotropyDirection::Given(reals("anisotropy_direction", d, 3, "vertex")?)
            }
        }
    };
    let args = SimplifyArguments {
        target,
        ratio,
        weights,
        memory,
        fixed: opt_indices("fixed", fixed)?,
        masses: opt_reals("masses", masses)?,
        mass_weight,
        point_faces,
        point_bary,
        refine_before,
        refine_after,
        repair_simplicial,
        anisotropy: anisotropy.map(|tau| AnisotropyArguments {
            tau,
            direction,
            clamp: anisotropy_clamp,
        }),
        representation: representation.to_owned(),
        projective_sign: projective_sign.to_owned(),
    };
    let out = py
        .detach(move || RustSimplification::new(&input, &args))
        .map_err(unwrap_error)?;
    Ok(Simplification {
        inner: Mutex::new(out),
    })
}

/// An intrinsic hierarchy of a mesh (paper Sec. 8.7): levels of decreasing vertex count from
/// one simplification run (or one per level), with the cotan Laplacian and lumped masses of
/// every level and the prolongations between them.
///
/// n_levels, level_sizes            counts (finest first; level 0 is the input)
/// level_vertices(k)                input vertex ids of level k
/// faces(k)                         3 vertex indices (into level k) per face
/// laplacian(k), prolongation(k)    COO triplets (rows, cols, values); prolongation(k) is
///                                  level k by level k+1
/// mass(k)                          lumped masses of level k
/// report                           dict: targets, counts, fresh_levels, prolongation_nnz,
///                                  simplify_seconds, operators_seconds
/// multigrid(...)                   a MultigridSolver on the finest level
#[pyclass(module = "ddg_tools", frozen)]
pub struct Hierarchy {
    inner: HierarchyHandle,
}

#[pymethods]
impl Hierarchy {
    /// Number of levels.
    #[getter]
    fn n_levels(&self) -> usize {
        self.inner.n_levels()
    }

    /// Vertices per level (list).
    #[getter]
    fn level_sizes(&self) -> Vec<usize> {
        self.inner.level_sizes()
    }

    /// Input vertex ids of level k.
    fn level_vertices<'py>(&self, py: Python<'py>, level: usize) -> PyResult<Bound<'py, PyAny>> {
        let v = self.inner.level_vertices(level).map_err(unwrap_error)?;
        index_array(py, "level_vertices", v.into_iter())
    }

    /// 3 vertex indices (into the level) per face of level k.
    fn faces<'py>(&self, py: Python<'py>, level: usize) -> PyResult<Bound<'py, PyAny>> {
        let v = self.inner.faces(level).map_err(unwrap_error)?;
        index_array(py, "faces", v.into_iter())
    }

    /// The cotan Laplacian of level k as COO triplets (rows, cols, values).
    fn laplacian<'py>(&self, py: Python<'py>, level: usize) -> PyResult<Bound<'py, PyTuple>> {
        coo_tuple(py, &self.inner.laplacian(level).map_err(unwrap_error)?)
    }

    /// The lumped masses of level k.
    fn mass<'py>(&self, py: Python<'py>, level: usize) -> PyResult<Bound<'py, PyAny>> {
        real_array(py, &self.inner.mass(level).map_err(unwrap_error)?)
    }

    /// P_k (level k by level k+1) as COO triplets (rows, cols, values).
    fn prolongation<'py>(&self, py: Python<'py>, level: usize) -> PyResult<Bound<'py, PyTuple>> {
        coo_tuple(py, &self.inner.prolongation(level).map_err(unwrap_error)?)
    }

    /// The report: dict.
    #[getter]
    fn report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let r = self.inner.report();
        let d = PyDict::new(py);
        d.set_item("targets", r.targets.clone())?;
        d.set_item("counts", r.counts.clone())?;
        d.set_item("fresh_levels", r.fresh_levels)?;
        d.set_item("prolongation_nnz", r.prolongation_nnz.clone())?;
        d.set_item("simplify_seconds", r.simplify_seconds)?;
        d.set_item("operators_seconds", r.operators_seconds)?;
        Ok(d)
    }

    /// The Galerkin multigrid solver (Liu et al. 2021) for a system on level 0.
    ///
    /// matrix       the system as a tuple (rows, cols, values) (square, symmetric positive
    ///              definite once `known` is eliminated); default L + shift M of level 0
    /// shift        the shift of the default matrix (default 0; a pure Laplacian is only
    ///              semidefinite: give `known`, a shift, or rely on the coarsest-level shift)
    /// known        level-0 vertices with prescribed values
    /// pre_sweeps, post_sweeps, tolerance, max_cycles   smoother and stopping rule
    #[pyo3(signature = (matrix = None, *, shift = 0.0, known = None, pre_sweeps = 2, post_sweeps = 2, tolerance = 1e-10, max_cycles = 50))]
    #[allow(clippy::too_many_arguments)]
    fn multigrid<'py>(
        &self,
        py: Python<'py>,
        matrix: Option<&Bound<'py, PyAny>>,
        shift: f64,
        known: Option<&Bound<'py, PyAny>>,
        pre_sweeps: usize,
        post_sweeps: usize,
        tolerance: f64,
        max_cycles: usize,
    ) -> PyResult<MultigridSolver> {
        let matrix = match matrix {
            Some(m) if !m.is_none() => {
                let t = m.cast::<PyTuple>().map_err(|_| {
                    PyTypeError::new_err("matrix: expected a tuple (rows, cols, values) of buffers")
                })?;
                if t.len() != 3 {
                    return Err(PyTypeError::new_err(
                        "matrix: expected a tuple (rows, cols, values)",
                    ));
                }
                let rows = indices("matrix (rows)", &t.get_item(0)?, 1, "entry")?;
                let cols = indices("matrix (cols)", &t.get_item(1)?, 1, "entry")?;
                let values = reals("matrix (values)", &t.get_item(2)?, 1, "entry")?;
                if rows.len() != cols.len() || rows.len() != values.len() {
                    return Err(super::InputError::new_err(format!(
                        "matrix: rows, cols and values must have the same length, got {}, {} and {}",
                        rows.len(),
                        cols.len(),
                        values.len()
                    )));
                }
                let n = self.inner.level_sizes()[0];
                Some(Coo {
                    nrows: n,
                    ncols: n,
                    rows,
                    cols,
                    values,
                })
            }
            _ => None,
        };
        let args = MultigridArguments {
            matrix,
            shift,
            known: opt_indices("known", known)?,
            pre_sweeps,
            post_sweeps,
            tolerance,
            max_cycles,
        };
        let solver = py
            .detach(|| self.inner.multigrid(&args))
            .map_err(unwrap_error)?;
        Ok(MultigridSolver {
            inner: Mutex::new(solver),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "<ddg_tools.Hierarchy levels {:?}>",
            self.inner.level_sizes()
        )
    }
}

/// A hierarchy of the mesh, for multigrid solves (paper Sec. 8.7).
///
/// ratio          vertex count of a level over the one above, in (0, 1) (default 0.25)
/// n_min          levels are added while the current one has more vertices (default 500)
/// fresh_levels   one simplification run per level with fresh masses (the reference's variant)
///                instead of one run with snapshots (default False)
/// weights, memory, fixed   as simplify
///
/// Raises InputError, MeshError and SolveError as simplify.
#[pyfunction]
#[pyo3(signature = (positions, triangles, *, ratio = 0.25, n_min = 500, fresh_levels = false, weights = 0.0, memory = true, fixed = None))]
#[allow(clippy::too_many_arguments)]
fn build_hierarchy<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    ratio: f64,
    n_min: usize,
    fresh_levels: bool,
    weights: f64,
    memory: bool,
    fixed: Option<&Bound<'py, PyAny>>,
) -> PyResult<Hierarchy> {
    let input = mesh_input(positions, triangles)?;
    let args = HierarchyArguments {
        ratio,
        n_min,
        fresh_levels,
        weights,
        memory,
        fixed: opt_indices("fixed", fixed)?,
    };
    let inner = py
        .detach(move || HierarchyHandle::new(&input, &args))
        .map_err(unwrap_error)?;
    Ok(Hierarchy { inner })
}

/// The multigrid solver of a Hierarchy.
///
/// solve(b, known_values=None, x0=None)   V-cycles until the relative residual is below the
///                tolerance or max_cycles: returns (x, report)
/// pcg(b, known_values=None, x0=None, max_iterations=100)   conjugate gradients preconditioned
///                by a V-cycle: returns (x, report)
/// info           dict: level_sizes, level_nnz, dropped_columns, coarsest_shift, n_known,
///                setup_seconds
///
/// b has one real per level-0 vertex; known_values one per `known` vertex (default zeros); x is
/// the full solution. The report is a dict: converged, cycles, residuals, rhs_norm,
/// relative_residual, convergence_factor, level_sizes, solve_seconds.
#[pyclass(module = "ddg_tools", frozen)]
pub struct MultigridSolver {
    inner: Mutex<MultigridHandle>,
}

impl MultigridSolver {
    fn lock(&self) -> PyResult<MutexGuard<'_, MultigridHandle>> {
        self.inner
            .lock()
            .map_err(|_| SolveError::new_err("the solver is poisoned by an earlier panic"))
    }
}

#[pymethods]
impl MultigridSolver {
    /// V-cycles. Returns (x, report).
    #[pyo3(signature = (b, known_values = None, x0 = None))]
    fn solve<'py>(
        &self,
        py: Python<'py>,
        b: &Bound<'py, PyAny>,
        known_values: Option<&Bound<'py, PyAny>>,
        x0: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyDict>)> {
        let b = reals("b", b, 1, "vertex")?;
        let kv = opt_reals("known_values", known_values)?;
        let x0 = opt_reals("x0", x0)?;
        let s = self.lock()?;
        let (x, r) = py
            .detach(|| s.solve(&b, kv.as_deref(), x0.as_deref()))
            .map_err(unwrap_error)?;
        Ok((real_array(py, &x)?, solve_report(py, &r)?))
    }

    /// Preconditioned conjugate gradients. Returns (x, report).
    #[pyo3(signature = (b, known_values = None, x0 = None, max_iterations = 100))]
    fn pcg<'py>(
        &self,
        py: Python<'py>,
        b: &Bound<'py, PyAny>,
        known_values: Option<&Bound<'py, PyAny>>,
        x0: Option<&Bound<'py, PyAny>>,
        max_iterations: usize,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyDict>)> {
        let b = reals("b", b, 1, "vertex")?;
        let kv = opt_reals("known_values", known_values)?;
        let x0 = opt_reals("x0", x0)?;
        let s = self.lock()?;
        let (x, r) = py
            .detach(|| s.pcg(&b, kv.as_deref(), x0.as_deref(), max_iterations))
            .map_err(unwrap_error)?;
        Ok((real_array(py, &x)?, solve_report(py, &r)?))
    }

    /// What the setup built: dict.
    #[getter]
    fn info<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let s = self.lock()?;
        let i = s.info();
        let d = PyDict::new(py);
        d.set_item("level_sizes", i.level_sizes.clone())?;
        d.set_item("level_nnz", i.level_nnz.clone())?;
        d.set_item("dropped_columns", i.dropped_columns.clone())?;
        d.set_item("coarsest_shift", i.coarsest_shift)?;
        d.set_item("n_known", i.n_known)?;
        d.set_item("setup_seconds", i.setup_seconds)?;
        Ok(d)
    }
}

/// The result of `qem_simplify`.
///
/// positions, triangles   the coarse mesh (3 reals per vertex, 3 indices per triangle)
/// vertex_map             per input vertex, the coarse vertex it collapsed into
/// projection_face, projection_bary, projection_distance   per input vertex, its closest point
///                        on the coarse mesh (face, 3 barycentric coordinates, distance)
/// report                 dict
#[pyclass(module = "ddg_tools", frozen)]
pub struct QemResult {
    out: QemOutput,
}

#[pymethods]
impl QemResult {
    /// array('d'), 3 per coarse vertex.
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        real_array(py, &self.out.positions)
    }
    /// array('i'), 3 per face.
    #[getter]
    fn triangles<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        index_array(py, "triangles", self.out.triangles.iter().copied())
    }
    /// array('i'), the coarse vertex of every input vertex.
    #[getter]
    fn vertex_map<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        index_array(py, "vertex_map", self.out.vertex_map.iter().copied())
    }
    /// array('i').
    #[getter]
    fn projection_face<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        index_array(
            py,
            "projection_face",
            self.out.projection_face.iter().copied(),
        )
    }
    /// array('d'), 3 per input vertex.
    #[getter]
    fn projection_bary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        real_array(py, &self.out.projection_bary)
    }
    /// array('d').
    #[getter]
    fn projection_distance<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        real_array(py, &self.out.projection_distance)
    }
    /// dict: n_input_vertices, n_output_vertices, n_output_faces, reached_target, n_collapses,
    /// n_rejected_link, n_rejected_flip, n_rejected_other, total_cost, fine_area, coarse_area,
    /// area_error_percent, max_projection_distance, mean_projection_distance.
    #[getter]
    fn report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let r = &self.out.report;
        let d = PyDict::new(py);
        d.set_item("n_input_vertices", r.n_input_vertices)?;
        d.set_item("n_input_faces", r.n_input_faces)?;
        d.set_item("n_output_vertices", r.n_output_vertices)?;
        d.set_item("n_output_faces", r.n_output_faces)?;
        d.set_item("reached_target", r.reached_target)?;
        d.set_item("n_collapses", r.n_collapses)?;
        d.set_item("n_rejected_link", r.n_rejected_link)?;
        d.set_item("n_rejected_flip", r.n_rejected_flip)?;
        d.set_item("n_rejected_other", r.n_rejected_other)?;
        d.set_item("total_cost", r.total_cost)?;
        d.set_item("fine_area", r.fine_area)?;
        d.set_item("coarse_area", r.coarse_area)?;
        d.set_item("area_error_percent", r.area_error_percent)?;
        d.set_item("max_projection_distance", r.max_projection_distance)?;
        d.set_item("mean_projection_distance", r.mean_projection_distance)?;
        Ok(d)
    }
}

/// Quadric error metric simplification (Garland, Heckbert 1997): the extrinsic baseline of the
/// paper's comparisons (Sec. 8.4). Collapses edges to `target` vertices; every input vertex is
/// also projected to its closest point on the result (the map used for the comparisons).
///
/// area_weighted  weight each face plane by its area (default False: weight 1)
/// boundary_weight  weight of the boundary constraint planes (default 1000)
///
/// Raises InputError, MeshError as simplify.
#[pyfunction]
#[pyo3(signature = (positions, triangles, *, target, area_weighted = false, boundary_weight = 1000.0))]
fn qem_simplify<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    triangles: &Bound<'py, PyAny>,
    target: usize,
    area_weighted: bool,
    boundary_weight: f64,
) -> PyResult<QemResult> {
    let input = mesh_input(positions, triangles)?;
    let args = QemArguments {
        target,
        area_weighted,
        boundary_weight,
    };
    let out = py
        .detach(move || rust_qem(&input, &args))
        .map_err(unwrap_error)?;
    Ok(QemResult { out })
}
