//! The PyO3 glue of the Points as Tori bindings: `PointsAsTori`, `PatCsgField`, `PatMeshTruth`,
//! `Pat2d` and the `pat_*` functions (also reachable as the module `ddg_tools.points_as_tori`
//! with the prefix dropped).
//!
//! No logic lives here: buffers (or plain Python sequences) are decoded, keyword arguments become
//! a [`Kw`] list, [`crate::points_as_tori`] does the work with the GIL released, results go back
//! as `array.array` objects (flat; `numpy.asarray(a).reshape(-1, 3)` for rows) and plain
//! dictionaries. Every error is a `PointsAsToriError`, which is both a `ddg_tools.DdgError` and a
//! `ValueError`, carrying the library's message.

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt, PyList, PyModule, PyString, PyTuple, PyType};

use super::{DdgError, InputError, index_array, indices, real_array, reals};
use crate::points_as_tori as core;
use crate::points_as_tori::{Kw, O, PErr, Pat, Shared, V};

static ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();

/// A progress counter and cancel flag for `PointsAsTori(...)` builds and evaluations: pass it as
/// `progress=`, poll `stage`, `done`, `total` or `fraction` from another thread, and `cancel()` to
/// make the running call raise `PointsAsToriError("cancelled")` at its next check.
#[pyclass(name = "PatProgress", module = "ddg_tools", frozen)]
struct PyProgress(ddg_points_as_tori::Progress);

#[pymethods]
impl PyProgress {
    #[new]
    fn new() -> Self {
        PyProgress(ddg_points_as_tori::Progress::new())
    }

    /// Ask the running call to stop (sticky until `reset()`).
    fn cancel(&self) {
        self.0.cancel();
    }

    /// Clear the cancel flag and the counters.
    fn reset(&self) {
        self.0.reset();
    }

    /// Whether `cancel()` was called.
    #[getter]
    fn cancelled(&self) -> bool {
        self.0.is_cancelled()
    }

    /// The current stage: "idle", "neighborhoods", "network", "fit", "tree" or "evaluate".
    #[getter]
    fn stage(&self) -> &'static str {
        self.0.stage()
    }

    /// Items finished in the current stage.
    #[getter]
    fn done(&self) -> u64 {
        self.0.done()
    }

    /// Items in the current stage.
    #[getter]
    fn total(&self) -> u64 {
        self.0.total()
    }

    /// `done / total` of the current stage.
    #[getter]
    fn fraction(&self) -> f64 {
        self.0.fraction()
    }

    fn __repr__(&self) -> String {
        format!("PatProgress(stage={:?}, done={}, total={}, cancelled={})", self.0.stage(), self.0.done(), self.0.total(), self.0.is_cancelled())
    }
}

fn err(e: PErr) -> PyErr {
    Python::attach(|py| match ERROR.get(py) {
        Some(t) => PyErr::from_type(t.bind(py).clone(), e.0),
        None => DdgError::new_err(e.0),
    })
}

fn bad(name: &str, message: impl std::fmt::Display) -> PyErr {
    err(PErr(format!("{name}: {message}")))
}

/// A buffer-decoding `InputError` becomes a `PointsAsToriError`.
fn as_input(e: PyErr, py: Python<'_>) -> PyErr {
    if e.is_instance_of::<InputError>(py) { err(PErr(e.value(py).to_string())) } else { e }
}

fn is_type_error(e: &PyErr, py: Python<'_>) -> bool {
    e.is_instance_of::<pyo3::exceptions::PyTypeError>(py)
}

/// Reals from a buffer, or from a (nested) Python sequence of numbers.
fn nums(name: &'static str, obj: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    // flat buffers, or (n, 3) / (n, 2) ones
    let mut first = None;
    for width in [1, 3, 2] {
        match reals(name, obj, width, "number") {
            Ok(v) => return Ok(v),
            Err(e) => {
                first.get_or_insert(e);
            }
        }
    }
    match Err::<Vec<f64>, PyErr>(first.unwrap_or_else(|| bad(name, "unreadable"))) {
        Ok(v) => Ok(v),
        Err(e) if is_type_error(&e, obj.py()) => {
            let flat: Vec<f64> = if let Ok(v) = obj.extract::<Vec<f64>>() {
                v
            } else if let Ok(v) = obj.extract::<Vec<Vec<f64>>>() {
                v.into_iter().flatten().collect()
            } else {
                return Err(bad(name, "expected an object with the buffer protocol or a sequence of numbers"));
            };
            Ok(flat)
        }
        Err(e) => Err(as_input(e, obj.py())),
    }
}

/// Indices from a buffer, or from a (nested) sequence of integers.
fn idx(name: &'static str, obj: &Bound<'_, PyAny>) -> PyResult<Vec<u32>> {
    let mut first = None;
    for width in [1, 3] {
        match indices(name, obj, width, "index") {
            Ok(v) => return Ok(v),
            Err(e) => {
                first.get_or_insert(e);
            }
        }
    }
    match Err::<Vec<u32>, PyErr>(first.unwrap_or_else(|| bad(name, "unreadable"))) {
        Ok(v) => Ok(v),
        Err(e) if is_type_error(&e, obj.py()) => {
            let flat: Vec<i64> = if let Ok(v) = obj.extract::<Vec<i64>>() {
                v
            } else if let Ok(v) = obj.extract::<Vec<Vec<i64>>>() {
                v.into_iter().flatten().collect()
            } else {
                return Err(bad(name, "expected an object with the buffer protocol or a sequence of integers"));
            };
            flat.into_iter()
                .map(|i| u32::try_from(i).map_err(|_| bad(name, format!("index {i} out of range"))))
                .collect()
        }
        Err(e) => Err(as_input(e, obj.py())),
    }
}

fn leak(name: &str) -> &'static str {
    static NAMES: std::sync::Mutex<Vec<&'static str>> = std::sync::Mutex::new(Vec::new());
    let mut names = NAMES.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(n) = names.iter().find(|n| **n == name) {
        return n;
    }
    let n: &'static str = Box::leak(name.to_owned().into_boxed_str());
    names.push(n);
    n
}

fn to_opt(name: &str, obj: &Bound<'_, PyAny>) -> PyResult<O> {
    if obj.is_none() {
        return Ok(O::None);
    }
    if obj.is_instance_of::<PyBool>() {
        return Ok(O::Bool(obj.extract()?));
    }
    if obj.is_instance_of::<PyInt>() {
        return obj.extract::<i64>().map(O::Int).map_err(|_| bad(name, "integer out of range"));
    }
    if obj.is_instance_of::<PyFloat>() {
        return Ok(O::Float(obj.extract()?));
    }
    if obj.is_instance_of::<PyString>() {
        return Ok(O::Str(obj.extract()?));
    }
    if let Ok(p) = obj.cast::<PyProgress>() {
        return Ok(O::Progress(p.get().0.clone()));
    }
    if let Ok(v) = nums(leak(name), obj) {
        return Ok(O::Array(v));
    }
    if let Ok(i) = obj.extract::<i64>()
        && obj.hasattr("__index__")?
    {
        return Ok(O::Int(i));
    }
    if let Ok(x) = obj.extract::<f64>() {
        return Ok(O::Float(x));
    }
    Err(bad(name, format!("unsupported value of type {}", obj.get_type().name().map(|n| n.to_string()).unwrap_or_default())))
}

/// Keyword arguments plus explicit extras (`None` values of extras are left out); the
/// `return_report` flag is split off.
fn kwargs(kw: Option<&Bound<'_, PyDict>>, extra: &[(&str, Option<&Bound<'_, PyAny>>)]) -> PyResult<(Kw, bool)> {
    let mut items = Vec::new();
    let mut report = false;
    if let Some(kw) = kw {
        for (k, v) in kw.iter() {
            let k: String = k.extract()?;
            if k == "return_report" {
                report = v.extract::<bool>().map_err(|_| bad("return_report", "expected a bool"))?;
                continue;
            }
            let o = to_opt(&k, &v)?;
            items.push((k, o));
        }
    }
    for (k, v) in extra {
        if let Some(v) = v {
            if items.iter().any(|(n, _)| n == k) {
                return Err(bad(k, "given twice"));
            }
            items.push(((*k).to_owned(), to_opt(k, v)?));
        }
    }
    Ok((Kw::new(items), report))
}

fn to_py<'py>(py: Python<'py>, v: &V) -> PyResult<Bound<'py, PyAny>> {
    Ok(match v {
        V::None => py.None().into_bound(py),
        V::B(b) => PyBool::new(py, *b).to_owned().into_any(),
        V::I(i) => i.into_pyobject(py)?.into_any(),
        V::F(x) => x.into_pyobject(py)?.into_any(),
        V::S(s) => s.into_pyobject(py)?.into_any(),
        V::L(items) => {
            let l = PyList::empty(py);
            for i in items {
                l.append(to_py(py, i)?)?;
            }
            l.into_any()
        }
        V::D(items) => {
            let d = PyDict::new(py);
            for (k, i) in items {
                d.set_item(k, to_py(py, i)?)?;
            }
            d.into_any()
        }
    })
}

fn dict<'py>(py: Python<'py>, v: &V) -> PyResult<Bound<'py, PyDict>> {
    to_py(py, v)?.cast_into::<PyDict>().map_err(PyErr::from)
}

fn arr<'py>(py: Python<'py>, v: &[f64]) -> PyResult<Bound<'py, PyAny>> {
    real_array(py, v)
}

fn tuple<'py>(py: Python<'py>, items: Vec<Bound<'py, PyAny>>) -> PyResult<Bound<'py, PyAny>> {
    Ok(PyTuple::new(py, items)?.into_any())
}

/// `(a, report)` or just `a`.
fn with_report<'py>(py: Python<'py>, mut items: Vec<Bound<'py, PyAny>>, report: &V, want: bool) -> PyResult<Bound<'py, PyAny>> {
    if want {
        items.push(to_py(py, report)?);
    }
    if items.len() == 1 { Ok(items.remove(0)) } else { tuple(py, items) }
}

// -------------------------------------------------------------------------------------------
// Registration
// -------------------------------------------------------------------------------------------

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let bases = PyTuple::new(py, [py.get_type::<DdgError>(), py.get_type::<PyValueError>()])?;
    let d = PyDict::new(py);
    d.set_item("__module__", "ddg_tools")?;
    d.set_item("__doc__", "A Points as Tori call failed: an invalid argument or option, a cloud that does not validate, a degenerate query, a missing GPU, ... Both a ddg_tools.DdgError and a ValueError; the message is the library's.")?;
    let cls = py.import("builtins")?.getattr("type")?.call1(("PointsAsToriError", bases, d))?.cast_into::<PyType>().map_err(PyErr::from)?;
    m.add("PointsAsToriError", &cls)?;
    let _ = ERROR.set(py, cls.unbind());
    m.add_class::<PyPat>()?;
    m.add_class::<PyProgress>()?;
    m.add_class::<PyCsg>()?;
    m.add_class::<PyTruth>()?;
    m.add_class::<PyPat2>()?;
    m.add_function(wrap_pyfunction!(pat_load_tori, m)?)?;
    m.add_function(wrap_pyfunction!(pat_from_tori_dict, m)?)?;
    m.add_function(wrap_pyfunction!(pat_boolean, m)?)?;
    m.add_function(wrap_pyfunction!(pat_offset, m)?)?;
    m.add_function(wrap_pyfunction!(pat_reconstruct, m)?)?;
    m.add_function(wrap_pyfunction!(pat_sample_level_set, m)?)?;
    m.add_function(wrap_pyfunction!(pat_morphology, m)?)?;
    m.add_function(wrap_pyfunction!(pat_orient_normals, m)?)?;
    m.add_function(wrap_pyfunction!(pat_point_areas, m)?)?;
    m.add_function(wrap_pyfunction!(pat_gwn, m)?)?;
    m.add_function(wrap_pyfunction!(pat_fwn, m)?)?;
    m.add_function(wrap_pyfunction!(pat_fwn_mesh, m)?)?;
    m.add_function(wrap_pyfunction!(pat_shc, m)?)?;
    m.add_function(wrap_pyfunction!(pat_sspd, m)?)?;
    m.add_function(wrap_pyfunction!(pat_lse, m)?)?;
    m.add_function(wrap_pyfunction!(pat_defaults, m)?)?;

    // The submodule `ddg_tools.points_as_tori`: the same objects with the short names.
    let sub = PyModule::new(py, "points_as_tori")?;
    sub.setattr("__doc__", "Points as Tori (Feng, Gkioulekas, Crane 2026): signed distance to oriented point clouds. The same objects as the pat_* names of ddg_tools, without the prefix.")?;
    let pairs: [(&str, &str); 23] = [
        ("PointsAsTori", "PointsAsTori"),
        ("CsgField", "PatCsgField"),
        ("Progress", "PatProgress"),
        ("MeshTruth", "PatMeshTruth"),
        ("Pat2d", "Pat2d"),
        ("PointsAsToriError", "PointsAsToriError"),
        ("load_tori", "pat_load_tori"),
        ("from_tori_dict", "pat_from_tori_dict"),
        ("boolean", "pat_boolean"),
        ("offset", "pat_offset"),
        ("reconstruct", "pat_reconstruct"),
        ("sample_level_set", "pat_sample_level_set"),
        ("morphology", "pat_morphology"),
        ("orient_normals", "pat_orient_normals"),
        ("point_areas", "pat_point_areas"),
        ("gwn", "pat_gwn"),
        ("fwn", "pat_fwn"),
        ("fwn_mesh", "pat_fwn_mesh"),
        ("shc", "pat_shc"),
        ("sspd", "pat_sspd"),
        ("lse", "pat_lse"),
        ("defaults", "pat_defaults"),
        ("gpu_info", "gpu_info"),
    ];
    for (short, long) in pairs {
        sub.setattr(short, m.getattr(long)?)?;
    }
    m.add_submodule(&sub)?;
    py.import("sys")?.getattr("modules")?.set_item("ddg_tools.points_as_tori", &sub)?;
    Ok(())
}

// -------------------------------------------------------------------------------------------
// Fields as arguments
// -------------------------------------------------------------------------------------------

fn shared(obj: &Bound<'_, PyAny>, name: &str) -> PyResult<Shared> {
    if let Ok(p) = obj.extract::<PyRef<'_, PyPat>>() {
        let f: Shared = p.pat.field.clone();
        return Ok(f);
    }
    if let Ok(c) = obj.extract::<PyRef<'_, PyCsg>>() {
        return Ok(c.field.clone());
    }
    Err(bad(name, "expected a PointsAsTori or a CsgField"))
}

// -------------------------------------------------------------------------------------------
// PointsAsTori
// -------------------------------------------------------------------------------------------

/// Signed distance to an oriented point cloud by blending fitted tori (Feng, Gkioulekas, Crane,
/// ACM TOG 45(4), 2026). Building runs the network and fits one torus per point; the object is
/// immutable and evaluates in parallel (the GIL is released). See `ddg_tools.pyi` for the options.
#[pyclass(name = "PointsAsTori", module = "ddg_tools", frozen)]
struct PyPat {
    pat: Arc<Pat>,
}

fn queries(q: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    nums("queries", q)
}

#[pymethods]
impl PyPat {
    #[new]
    #[pyo3(signature = (points, normals, **options))]
    fn new(py: Python<'_>, points: &Bound<'_, PyAny>, normals: &Bound<'_, PyAny>, options: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let (p, n) = (nums("points", points)?, nums("normals", normals)?);
        let (mut kw, _) = kwargs(options, &[])?;
        let pat = py
            .detach(move || {
                let spec = core::parse_build(&mut kw, false)?;
                kw.finish("PointsAsTori")?;
                Pat::build(&p, &n, &spec)
            })
            .map_err(err)?;
        Ok(PyPat { pat: Arc::new(pat) })
    }

    /// Number of points.
    #[getter]
    fn n(&self) -> usize {
        self.pat.cloud.len()
    }

    /// The sign rule the tori were fitted with ("reference" or "geometric"), or None when not
    /// known (tori imported from a dict or from a file written without the marker).
    #[getter]
    fn sign_rule(&self) -> Option<&'static str> {
        self.pat.field.sign_rule().map(|r| if r == ddg_points_as_tori::SignRule::Reference { "reference" } else { "geometric" })
    }

    /// The (validated) points, flat `array('d')`, 3 per point.
    #[getter]
    fn points<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        arr(py, self.pat.cloud.positions().as_flattened())
    }

    /// The (validated, unit) normals, flat `array('d')`, 3 per point.
    #[getter]
    fn normals<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        arr(py, self.pat.cloud.normals().as_flattened())
    }

    /// `φ(x)` for flat queries (3 numbers each); `array('d')`.
    #[pyo3(signature = (queries))]
    fn signed_distance<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, p) = (self::queries(queries)?, self.pat.clone());
        let (d, _, _) = py.detach(move || p.evaluate(&q, true, false, false)).map_err(err)?;
        arr(py, &d)
    }

    /// Alias of `signed_distance`.
    #[pyo3(signature = (queries))]
    fn distance<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        self.signed_distance(py, queries)
    }

    /// `∇φ(x)`, flat `array('d')`, 3 per query.
    #[pyo3(signature = (queries))]
    fn gradient<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, p) = (self::queries(queries)?, self.pat.clone());
        let (_, g, _) = py.detach(move || p.evaluate(&q, false, true, false)).map_err(err)?;
        arr(py, &g)
    }

    /// `(φ, ∇φ)` in one pass.
    #[pyo3(signature = (queries))]
    fn distance_and_gradient<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, p) = (self::queries(queries)?, self.pat.clone());
        let (d, g, _) = py.detach(move || p.evaluate(&q, true, true, false)).map_err(err)?;
        tuple(py, vec![arr(py, &d)?, arr(py, &g)?])
    }

    /// `Δφ(x)` (frozen λ).
    #[pyo3(signature = (queries))]
    fn laplacian<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, p) = (self::queries(queries)?, self.pat.clone());
        let (_, _, l) = py.detach(move || p.evaluate(&q, false, false, true)).map_err(err)?;
        arr(py, &l)
    }

    /// Write our binary tori file (tied to the cloud by a hash).
    fn save_tori(&self, py: Python<'_>, path: &str) -> PyResult<()> {
        let (p, path) = (self.pat.clone(), path.to_owned());
        py.detach(move || p.save_tori(&path)).map_err(err)
    }

    /// The tori in the reference program's layout: a dict `centers` (3n), `axes` (2n, the
    /// angles theta, phi), `major` (n), `minor` (n, signed) of flat `array('d')`.
    fn export_tori<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let t = core::export_tori(self.pat.tori());
        let d = PyDict::new(py);
        d.set_item("centers", arr(py, &t.centers)?)?;
        d.set_item("axes", arr(py, &t.axes)?)?;
        d.set_item("major", arr(py, &t.major)?)?;
        d.set_item("minor", arr(py, &t.minor)?)?;
        Ok(d)
    }

    /// All the tori: the keys of `export_tori` plus `axis_vectors` (3n unit axes), `sigma` (n)
    /// and `coefficients` (6n, the network's quadric coefficients).
    #[getter]
    fn tori<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let t = core::export_tori(self.pat.tori());
        let d = self.export_tori(py)?;
        d.set_item("axis_vectors", arr(py, &t.axis_vectors)?)?;
        d.set_item("sigma", arr(py, &t.sigma)?)?;
        d.set_item("coefficients", arr(py, &t.coefficients)?)?;
        Ok(d)
    }

    /// What building did (timings, devices, fit counters); None for a field loaded from tori.
    #[getter]
    fn build_report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        match &self.pat.build_report {
            Some(r) => to_py(py, &core::build_report_v(r)),
            None => Ok(py.None().into_bound(py)),
        }
    }

    /// What the last evaluation did on this object (None before the first).
    #[getter]
    fn eval_report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        match self.pat.last_eval.lock().map_err(|_| bad("eval_report", "poisoned"))?.as_ref() {
            Some(r) => to_py(py, &core::eval_report_v(r)),
            None => Ok(py.None().into_bound(py)),
        }
    }

    /// What validating the cloud found.
    #[getter]
    fn cloud_report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        to_py(py, &core::cloud_report_v(&self.pat.cloud_report))
    }

    /// The evaluation options in force, as a dict.
    #[getter]
    fn options<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        to_py(py, &self.pat.options())
    }

    /// The same tori with other evaluation options (neighborhood, lambda_rule, gradient, ...).
    #[pyo3(signature = (**options))]
    fn with_options(&self, py: Python<'_>, options: Option<&Bound<'_, PyDict>>) -> PyResult<PyPat> {
        let (mut kw, _) = kwargs(options, &[])?;
        let p = self.pat.clone();
        let q = py
            .detach(move || {
                let r = p.with_options(&mut kw)?;
                kw.finish("with_options")?;
                Ok(r)
            })
            .map_err(err)?;
        Ok(PyPat { pat: Arc::new(q) })
    }

    /// Marching cubes on `φ = level`: `(vertices, faces)` (flat `array('d')` and `array('i')`),
    /// plus normals with `normals=True` and the report with `return_report=True`.
    #[pyo3(signature = (level=None, resolution=None, bbox=None, **options))]
    fn reconstruct<'py>(
        &self,
        py: Python<'py>,
        level: Option<&Bound<'py, PyAny>>,
        resolution: Option<&Bound<'py, PyAny>>,
        bbox: Option<&Bound<'py, PyAny>>,
        options: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let f: Shared = self.pat.field.clone();
        reconstruct_impl(py, f, Some(self.pat.clone()), level, resolution, bbox, options)
    }

    /// Points and unit normals on `φ = level` (area-uniform, projected): `(points, normals)`.
    #[pyo3(signature = (level=None, count=None, **options))]
    fn sample_level_set<'py>(&self, py: Python<'py>, level: Option<&Bound<'py, PyAny>>, count: Option<&Bound<'py, PyAny>>, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
        let f: Shared = self.pat.field.clone();
        sample_impl(py, f, level, count, options)
    }

    /// The offset surface `φ - c` as a CSG field.
    fn offset(&self, c: f64) -> PyResult<PyCsg> {
        let f: Shared = self.pat.field.clone();
        Ok(PyCsg { field: Arc::new(core::offset(f, c).map_err(err)?) })
    }

    /// `min(self, other)` (or a smooth minimum with radius `smooth`).
    #[pyo3(signature = (other, smooth=None))]
    fn union(&self, other: &Bound<'_, PyAny>, smooth: Option<f64>) -> PyResult<PyCsg> {
        boolean_impl(self.pat.field.clone(), other, "union", smooth)
    }

    /// `max(self, other)`.
    #[pyo3(signature = (other, smooth=None))]
    fn intersect(&self, other: &Bound<'_, PyAny>, smooth: Option<f64>) -> PyResult<PyCsg> {
        boolean_impl(self.pat.field.clone(), other, "intersect", smooth)
    }

    /// `max(self, -other)`.
    #[pyo3(signature = (other, smooth=None))]
    fn subtract(&self, other: &Bound<'_, PyAny>, smooth: Option<f64>) -> PyResult<PyCsg> {
        boolean_impl(self.pat.field.clone(), other, "subtract", smooth)
    }

    fn __repr__(&self) -> String {
        format!("PointsAsTori(n={})", self.pat.cloud.len())
    }
}

fn boolean_impl(a: Shared, other: &Bound<'_, PyAny>, op: &str, smooth: Option<f64>) -> PyResult<PyCsg> {
    let b = shared(other, "other")?;
    Ok(PyCsg { field: Arc::new(core::boolean(a, b, op, smooth).map_err(err)?) })
}

/// Marching cubes; with `sink` (the `PointsAsTori` object), `eval_report` afterwards describes the
/// evaluations of the run.
fn reconstruct_impl<'py>(
    py: Python<'py>,
    f: Shared,
    sink: Option<Arc<Pat>>,
    level: Option<&Bound<'py, PyAny>>,
    resolution: Option<&Bound<'py, PyAny>>,
    bbox: Option<&Bound<'py, PyAny>>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyAny>> {
    let (mut kw, want) = kwargs(options, &[("level", level), ("resolution", resolution), ("bbox", bbox)])?;
    let want_normals = options.is_some_and(|o| o.get_item("normals").ok().flatten().is_some_and(|v| v.extract::<bool>().unwrap_or(false)));
    let m = py.detach(move || core::reconstruct(&*f, &mut kw)).map_err(err)?;
    if let (Some(p), Some(e)) = (sink, m.eval.clone())
        && let Ok(mut last) = p.last_eval.lock()
    {
        *last = Some(e);
    }
    let mut items = vec![arr(py, &m.positions)?, index_array(py, "faces", m.faces.iter().copied())?];
    if want_normals && let Some(n) = &m.normals {
        items.push(arr(py, n)?);
    }
    with_report(py, items, &m.report, want)
}

fn sample_impl<'py>(
    py: Python<'py>,
    f: Shared,
    level: Option<&Bound<'py, PyAny>>,
    count: Option<&Bound<'py, PyAny>>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyAny>> {
    let (mut kw, want) = kwargs(options, &[("count", count)])?;
    let level = match level {
        Some(l) if !l.is_none() => l.extract::<f64>().map_err(|_| bad("level", "expected a number"))?,
        _ => 0.0,
    };
    let (p, n, rep) = py.detach(move || core::level_set(&*f, level, &mut kw)).map_err(err)?;
    with_report(py, vec![arr(py, &p)?, arr(py, &n)?], &rep, want)
}

// -------------------------------------------------------------------------------------------
// Tori I/O
// -------------------------------------------------------------------------------------------

/// A field from a file written by `PointsAsTori.save_tori` and the cloud it belongs to (a
/// different cloud is an error). Keyword options as in `PointsAsTori` (evaluation options,
/// `device`, `threads`; no network is run).
#[pyfunction]
#[pyo3(signature = (path, points, normals, **options))]
fn pat_load_tori(py: Python<'_>, path: &str, points: &Bound<'_, PyAny>, normals: &Bound<'_, PyAny>, options: Option<&Bound<'_, PyDict>>) -> PyResult<PyPat> {
    let (p, n, path) = (nums("points", points)?, nums("normals", normals)?, path.to_owned());
    let (mut kw, _) = kwargs(options, &[])?;
    let pat = py
        .detach(move || {
            let spec = core::parse_build(&mut kw, false)?;
            kw.finish("load_tori")?;
            let (tori, rule) = core::load_tori_for(&path, &p, &n, spec.policy)?;
            Pat::from_tori(&p, &n, tori, rule, &spec)
        })
        .map_err(err)?;
    Ok(PyPat { pat: Arc::new(pat) })
}

/// A field from tori in the reference program's layout (`export_tori()`): a dict with
/// `centers`, `axes` (angles), `major`, `minor`; and the cloud they belong to.
#[pyfunction]
#[pyo3(signature = (points, normals, tori, **options))]
fn pat_from_tori_dict(py: Python<'_>, points: &Bound<'_, PyAny>, normals: &Bound<'_, PyAny>, tori: &Bound<'_, PyDict>, options: Option<&Bound<'_, PyDict>>) -> PyResult<PyPat> {
    let (p, n) = (nums("points", points)?, nums("normals", normals)?);
    let get = |k: &'static str| -> PyResult<Vec<f64>> {
        match tori.get_item(k)? {
            Some(v) => nums(k, &v),
            None => Err(bad("tori", format!("missing key {k:?} (centers, axes, major, minor)"))),
        }
    };
    let (c, a, r, s) = (get("centers")?, get("axes")?, get("major")?, get("minor")?);
    let (mut kw, _) = kwargs(options, &[])?;
    let pat = py
        .detach(move || {
            let spec = core::parse_build(&mut kw, false)?;
            kw.finish("from_tori_dict")?;
            let t = core::import_tori(&c, &a, &r, &s)?;
            Pat::from_tori(&p, &n, t, None, &spec)
        })
        .map_err(err)?;
    Ok(PyPat { pat: Arc::new(pat) })
}

// -------------------------------------------------------------------------------------------
// CSG
// -------------------------------------------------------------------------------------------

/// A Boolean or offset of PointsAsTori fields: an implicit field with the same evaluation and
/// reconstruction methods.
#[pyclass(name = "PatCsgField", module = "ddg_tools", frozen)]
struct PyCsg {
    field: Shared,
}

#[pymethods]
impl PyCsg {
    /// `φ(x)`.
    fn signed_distance<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, f) = (self::queries(queries)?, self.field.clone());
        let (d, _) = py.detach(move || core::eval_shared(&*f, &q, false)).map_err(err)?;
        arr(py, &d)
    }

    /// `(φ, ∇φ)`.
    fn distance_and_gradient<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, f) = (self::queries(queries)?, self.field.clone());
        let (d, g) = py.detach(move || core::eval_shared(&*f, &q, true)).map_err(err)?;
        tuple(py, vec![arr(py, &d)?, arr(py, &g)?])
    }

    /// `∇φ`.
    fn gradient<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, f) = (self::queries(queries)?, self.field.clone());
        let (_, g) = py.detach(move || core::eval_shared(&*f, &q, true)).map_err(err)?;
        arr(py, &g)
    }

    /// A box containing the surface, `(xmin, ymin, zmin, xmax, ymax, zmax)`.
    #[getter]
    fn bbox<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        arr(py, &core::bbox(&*self.field))
    }

    /// Marching cubes, as `PointsAsTori.reconstruct`.
    #[pyo3(signature = (level=None, resolution=None, bbox=None, **options))]
    fn reconstruct<'py>(
        &self,
        py: Python<'py>,
        level: Option<&Bound<'py, PyAny>>,
        resolution: Option<&Bound<'py, PyAny>>,
        bbox: Option<&Bound<'py, PyAny>>,
        options: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        reconstruct_impl(py, self.field.clone(), None, level, resolution, bbox, options)
    }

    /// Points and normals on `φ = level`, as `PointsAsTori.sample_level_set` (`count` defaults to 4096).
    #[pyo3(signature = (level=None, count=None, **options))]
    fn sample_level_set<'py>(&self, py: Python<'py>, level: Option<&Bound<'py, PyAny>>, count: Option<&Bound<'py, PyAny>>, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
        sample_impl(py, self.field.clone(), level, count, options)
    }

    /// The offset surface `φ - c`.
    fn offset(&self, c: f64) -> PyResult<PyCsg> {
        Ok(PyCsg { field: Arc::new(core::offset(self.field.clone(), c).map_err(err)?) })
    }

    /// `min(self, other)`.
    #[pyo3(signature = (other, smooth=None))]
    fn union(&self, other: &Bound<'_, PyAny>, smooth: Option<f64>) -> PyResult<PyCsg> {
        boolean_impl(self.field.clone(), other, "union", smooth)
    }

    /// `max(self, other)`.
    #[pyo3(signature = (other, smooth=None))]
    fn intersect(&self, other: &Bound<'_, PyAny>, smooth: Option<f64>) -> PyResult<PyCsg> {
        boolean_impl(self.field.clone(), other, "intersect", smooth)
    }

    /// `max(self, -other)`.
    #[pyo3(signature = (other, smooth=None))]
    fn subtract(&self, other: &Bound<'_, PyAny>, smooth: Option<f64>) -> PyResult<PyCsg> {
        boolean_impl(self.field.clone(), other, "subtract", smooth)
    }
}

/// A Boolean of two fields: `op` is `"union"`, `"intersect"` or `"subtract"`; `smooth` a radius
/// for the smooth variants.
#[pyfunction]
#[pyo3(signature = (a, b, op, smooth=None))]
fn pat_boolean(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>, op: &str, smooth: Option<f64>) -> PyResult<PyCsg> {
    boolean_impl(shared(a, "a")?, b, op, smooth)
}

/// The offset field `φ - c` (outward for `c > 0`).
#[pyfunction]
fn pat_offset(a: &Bound<'_, PyAny>, c: f64) -> PyResult<PyCsg> {
    Ok(PyCsg { field: Arc::new(core::offset(shared(a, "a")?, c).map_err(err)?) })
}

// -------------------------------------------------------------------------------------------
// Functions on fields
// -------------------------------------------------------------------------------------------

/// Marching cubes of a PointsAsTori or CsgField; see `PointsAsTori.reconstruct`.
#[pyfunction]
#[pyo3(signature = (field, level=None, resolution=None, bbox=None, **options))]
fn pat_reconstruct<'py>(
    py: Python<'py>,
    field: &Bound<'py, PyAny>,
    level: Option<&Bound<'py, PyAny>>,
    resolution: Option<&Bound<'py, PyAny>>,
    bbox: Option<&Bound<'py, PyAny>>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyAny>> {
    reconstruct_impl(py, shared(field, "field")?, None, level, resolution, bbox, options)
}

/// Points and normals on a level set; see `PointsAsTori.sample_level_set`.
#[pyfunction]
#[pyo3(signature = (field, level=None, count=None, **options))]
fn pat_sample_level_set<'py>(py: Python<'py>, field: &Bound<'py, PyAny>, level: Option<&Bound<'py, PyAny>>, count: Option<&Bound<'py, PyAny>>, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    sample_impl(py, shared(field, "field")?, level, count, options)
}

/// Morphology (Fig. 24): `op` is `"erode"`, `"dilate"`, `"open"` or `"close"` with radius
/// `radius`; the level set is sampled and tori are refitted. Returns `(PointsAsTori, report)`.
/// Keyword options: those of `sample_level_set` and of `PointsAsTori` (with `paper=True` by default).
#[pyfunction]
#[pyo3(signature = (field, op, radius, **options))]
fn pat_morphology<'py>(py: Python<'py>, field: &Bound<'py, PyAny>, op: &str, radius: f64, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    let (f, op) = (shared(field, "field")?, op.to_owned());
    let (mut kw, _) = kwargs(options, &[])?;
    let (p, rep) = py.detach(move || core::morphology(&*f, &op, radius, &mut kw)).map_err(err)?;
    tuple(py, vec![Bound::new(py, PyPat { pat: Arc::new(p) })?.into_any(), to_py(py, &rep)?])
}

/// Orient normals by the iteration of Fig. 27. Returns `(normals, report)`.
#[pyfunction]
#[pyo3(signature = (points, normals, **options))]
fn pat_orient_normals<'py>(py: Python<'py>, points: &Bound<'py, PyAny>, normals: &Bound<'py, PyAny>, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    let (p, n) = (nums("points", points)?, nums("normals", normals)?);
    let (mut kw, _) = kwargs(options, &[])?;
    let (out, rep) = py.detach(move || core::orient(&p, &n, &mut kw)).map_err(err)?;
    tuple(py, vec![arr(py, &out)?, to_py(py, &rep)?])
}

// -------------------------------------------------------------------------------------------
// Baselines
// -------------------------------------------------------------------------------------------

/// Barill point areas (`k`, `filter_opposite`); `areas`, or `(areas, report)` with `return_report=True`.
#[pyfunction]
#[pyo3(signature = (points, normals, **options))]
fn pat_point_areas<'py>(py: Python<'py>, points: &Bound<'py, PyAny>, normals: &Bound<'py, PyAny>, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    let (p, n) = (nums("points", points)?, nums("normals", normals)?);
    let (mut kw, want) = kwargs(options, &[])?;
    let (a, rep) = py.detach(move || core::areas(&p, &n, &mut kw)).map_err(err)?;
    with_report(py, vec![arr(py, &a)?], &rep, want)
}

/// Generalized winding number (direct sum; `regularization=eps` for the regularized one).
#[pyfunction]
#[pyo3(signature = (points, normals, queries, **options))]
fn pat_gwn<'py>(py: Python<'py>, points: &Bound<'py, PyAny>, normals: &Bound<'py, PyAny>, queries: &Bound<'py, PyAny>, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    let (p, n, q) = (nums("points", points)?, nums("normals", normals)?, self::queries(queries)?);
    let (mut kw, want) = kwargs(options, &[])?;
    let (w, rep) = py.detach(move || core::gwn(&p, &n, &q, &mut kw)).map_err(err)?;
    with_report(py, vec![arr(py, &w)?], &rep, want)
}

/// Barill's fast winding number of a point cloud (`beta`, `order`, `leaf_size`).
#[pyfunction]
#[pyo3(signature = (points, normals, queries, **options))]
fn pat_fwn<'py>(py: Python<'py>, points: &Bound<'py, PyAny>, normals: &Bound<'py, PyAny>, queries: &Bound<'py, PyAny>, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    let (p, n, q) = (nums("points", points)?, nums("normals", normals)?, self::queries(queries)?);
    let (mut kw, want) = kwargs(options, &[])?;
    let (w, rep) = py.detach(move || core::fwn(&p, &n, &q, &mut kw)).map_err(err)?;
    with_report(py, vec![arr(py, &w)?], &rep, want)
}

/// Barill's fast winding number of a triangle mesh.
#[pyfunction]
#[pyo3(signature = (vertices, faces, queries, **options))]
fn pat_fwn_mesh<'py>(py: Python<'py>, vertices: &Bound<'py, PyAny>, faces: &Bound<'py, PyAny>, queries: &Bound<'py, PyAny>, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    let (v, f, q) = (nums("vertices", vertices)?, idx("faces", faces)?, self::queries(queries)?);
    let (mut kw, want) = kwargs(options, &[])?;
    let (w, rep) = py.detach(move || core::fwn_mesh(&v, &f, &q, &mut kw)).map_err(err)?;
    with_report(py, vec![arr(py, &w)?], &rep, want)
}

fn kernel_impl<'py>(
    py: Python<'py>,
    method: &'static str,
    points: &Bound<'py, PyAny>,
    normals: &Bound<'py, PyAny>,
    queries: &Bound<'py, PyAny>,
    gradient: bool,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyAny>> {
    let (p, n, q) = (nums("points", points)?, nums("normals", normals)?, self::queries(queries)?);
    let (mut kw, want) = kwargs(options, &[])?;
    let (v, g, rep) = py.detach(move || core::kernel(method, &p, &n, &q, gradient, &mut kw)).map_err(err)?;
    let mut items = vec![arr(py, &v)?];
    if gradient {
        items.push(arr(py, &g)?);
    }
    with_report(py, items, &rep, want)
}

/// Signed Hopf-Cole distance (Eq. 28). With `gradient=True` returns `(values, gradients)`.
#[pyfunction]
#[pyo3(signature = (points, normals, queries, gradient=false, **options))]
fn pat_shc<'py>(py: Python<'py>, points: &Bound<'py, PyAny>, normals: &Bound<'py, PyAny>, queries: &Bound<'py, PyAny>, gradient: bool, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    kernel_impl(py, "shc", points, normals, queries, gradient, options)
}

/// Smoothed signed planar distance (Eq. 29).
#[pyfunction]
#[pyo3(signature = (points, normals, queries, gradient=false, **options))]
fn pat_sspd<'py>(py: Python<'py>, points: &Bound<'py, PyAny>, normals: &Bound<'py, PyAny>, queries: &Bound<'py, PyAny>, gradient: bool, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    kernel_impl(py, "sspd", points, normals, queries, gradient, options)
}

/// LogSumExp unsigned distance (Madan and Levin 2022).
#[pyfunction]
#[pyo3(signature = (points, normals, queries, gradient=false, **options))]
fn pat_lse<'py>(py: Python<'py>, points: &Bound<'py, PyAny>, normals: &Bound<'py, PyAny>, queries: &Bound<'py, PyAny>, gradient: bool, options: Option<&Bound<'py, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
    kernel_impl(py, "lse", points, normals, queries, gradient, options)
}

/// The defaults of the `PointsAsTori` options, read from the library.
#[pyfunction]
fn pat_defaults(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    dict(py, &core::defaults())
}

// -------------------------------------------------------------------------------------------
// Exact mesh signed distance
// -------------------------------------------------------------------------------------------

/// The exact signed distance to a closed triangle mesh (closest point on a BVH, sign by the
/// winding number), the ground truth of the evaluation.
#[pyclass(name = "PatMeshTruth", module = "ddg_tools", frozen)]
struct PyTruth {
    t: Arc<core::Truth>,
}

#[pymethods]
impl PyTruth {
    #[new]
    #[pyo3(signature = (vertices, faces, **options))]
    fn new(py: Python<'_>, vertices: &Bound<'_, PyAny>, faces: &Bound<'_, PyAny>, options: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let (v, f) = (nums("vertices", vertices)?, idx("faces", faces)?);
        let (mut kw, _) = kwargs(options, &[])?;
        let t = py.detach(move || core::truth(&v, &f, &mut kw)).map_err(err)?;
        Ok(PyTruth { t: Arc::new(t) })
    }

    /// Topology of the mesh (closed, boundary edges, signed volume, ...).
    #[getter]
    fn mesh_report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        to_py(py, &self.t.report)
    }

    /// Signed distance (negative inside).
    fn signed_distance<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, t) = (self::queries(queries)?, self.t.clone());
        let (d, _, _) = py.detach(move || core::truth_eval(&t.truth, &q, false)).map_err(err)?;
        arr(py, &d)
    }

    /// `(signed distance, winding number)`.
    fn evaluate<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, t) = (self::queries(queries)?, self.t.clone());
        let (d, w, _) = py.detach(move || core::truth_eval(&t.truth, &q, true)).map_err(err)?;
        tuple(py, vec![arr(py, &d)?, arr(py, &w)?])
    }
}

// -------------------------------------------------------------------------------------------
// The classical 2D variant
// -------------------------------------------------------------------------------------------

/// The classical 2D variant (Fig. 10 analogue): a circle per point from a weighted quadratic
/// fit of its neighbors, blended like the tori. Points and normals are 2 numbers per row.
#[pyclass(name = "Pat2d", module = "ddg_tools", frozen)]
struct PyPat2 {
    p: Arc<core::Pat2>,
}

#[pymethods]
impl PyPat2 {
    #[new]
    #[pyo3(signature = (points, normals, **options))]
    fn new(py: Python<'_>, points: &Bound<'_, PyAny>, normals: &Bound<'_, PyAny>, options: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let (p, n) = (nums("points", points)?, nums("normals", normals)?);
        let (mut kw, _) = kwargs(options, &[])?;
        let f = py.detach(move || core::build_2d(&p, &n, &mut kw)).map_err(err)?;
        Ok(PyPat2 { p: Arc::new(f) })
    }

    /// What building did.
    #[getter]
    fn build_report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        to_py(py, &self.p.report)
    }

    /// `φ(x)` for flat queries (2 numbers each).
    fn signed_distance<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, f) = (self::queries(queries)?, self.p.clone());
        let (d, _) = py.detach(move || core::eval_2d(&f.field, &q, false)).map_err(err)?;
        arr(py, &d)
    }

    /// `(φ, ∇φ)`.
    fn distance_and_gradient<'py>(&self, py: Python<'py>, queries: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let (q, f) = (self::queries(queries)?, self.p.clone());
        let (d, g) = py.detach(move || core::eval_2d(&f.field, &q, true)).map_err(err)?;
        tuple(py, vec![arr(py, &d)?, arr(py, &g)?])
    }

    /// Contour polylines of `φ = level` on an `nx` by `ny` grid over `bbox` =
    /// `(xmin, ymin, xmax, ymax)`: a list of `(points, closed)`.
    #[pyo3(signature = (bbox, nx, ny, level=0.0))]
    fn contours<'py>(&self, py: Python<'py>, bbox: &Bound<'py, PyAny>, nx: usize, ny: usize, level: f64) -> PyResult<Bound<'py, PyAny>> {
        let b = nums("bbox", bbox)?;
        let b: [f64; 4] = b.try_into().map_err(|_| bad("bbox", "expected 4 numbers (xmin, ymin, xmax, ymax)"))?;
        let f = self.p.clone();
        let lines = py.detach(move || core::contours_2d(&f.field, b, nx, ny, level)).map_err(err)?;
        let out = PyList::empty(py);
        for (pts, closed) in lines {
            out.append(tuple(py, vec![arr(py, &pts)?, PyBool::new(py, closed).to_owned().into_any()])?)?;
        }
        Ok(out.into_any())
    }
}
