//! The PyO3 glue of `ddg_tools.RepulsiveSurfaceSession`, `repulsive_surface_flow`,
//! `surface_energy`, `repulsive_surface_defaults` and `repulsive_surface_scene`.
//!
//! No logic lives here: buffers (or plain Python sequences) are decoded, option values become a
//! [`Val`] tree, [`crate::repulsive_surfaces`] does the work with the GIL released (taking it
//! back only to call the progress function and to check for Ctrl-C), results go back as
//! `array.array` objects and plain dictionaries. Errors in the arguments or options are
//! `RepulsiveSurfaceInputError`, which is both a `ddg_tools.InputError` and a `ValueError`; a
//! failed computation is `SolveError`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt, PyList, PyString, PyTuple, PyType};

use super::{InputError, SolveError, index_array, indices, new_array, real_array, reals};
use crate::repulsive_surfaces::{Engine, RunEnd, SError, Val, make_surface};

static INPUT_ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let bases = PyTuple::new(
        py,
        [py.get_type::<InputError>(), py.get_type::<PyValueError>()],
    )?;
    let dict = PyDict::new(py);
    dict.set_item("__module__", "ddg_tools")?;
    dict.set_item(
        "__doc__",
        "An argument or an option of the repulsive-surfaces functions is invalid. Both a \
         ddg_tools.InputError and a ValueError; the message is the library's.",
    )?;
    let cls = py
        .import("builtins")?
        .getattr("type")?
        .call1(("RepulsiveSurfaceInputError", bases, dict))?
        .cast_into::<PyType>()
        .map_err(PyErr::from)?;
    m.add("RepulsiveSurfaceInputError", &cls)?;
    let _ = INPUT_ERROR.set(py, cls.unbind());
    m.add_class::<RepulsiveSurfaceSession>()?;
    m.add_class::<RepulsiveSurfaceHandle>()?;
    m.add_function(wrap_pyfunction!(repulsive_surface_flow, m)?)?;
    m.add_function(wrap_pyfunction!(surface_energy, m)?)?;
    m.add_function(wrap_pyfunction!(repulsive_surface_defaults, m)?)?;
    m.add_function(wrap_pyfunction!(gpu_info, m)?)?;
    m.add_function(wrap_pyfunction!(repulsive_surface_scene, m)?)?;
    Ok(())
}

fn input_err(message: String) -> PyErr {
    Python::attach(|py| match INPUT_ERROR.get(py) {
        Some(t) => PyErr::from_type(t.bind(py).clone(), message),
        None => InputError::new_err(message),
    })
}

fn s_err(e: SError) -> PyErr {
    match e {
        SError::Input(m) => input_err(m),
        SError::Solve(m) => SolveError::new_err(m),
    }
}

fn bad(name: &str, message: impl std::fmt::Display) -> PyErr {
    input_err(format!("{name}: {message}"))
}

// -------------------------------------------------------------------------------------------
// Decoding
// -------------------------------------------------------------------------------------------

/// A buffer-decoding `InputError` becomes a `RepulsiveSurfaceInputError`.
fn as_input(e: PyErr, py: Python<'_>) -> PyErr {
    if e.is_instance_of::<InputError>(py) {
        input_err(e.value(py).to_string())
    } else {
        e
    }
}

fn is_type_error(e: &PyErr, py: Python<'_>) -> bool {
    e.is_instance_of::<pyo3::exceptions::PyTypeError>(py)
}

/// Reals from a buffer, or from a (nested) Python sequence of numbers.
fn reals_any(
    name: &'static str,
    obj: &Bound<'_, PyAny>,
    width: usize,
    row: &'static str,
) -> PyResult<Vec<f64>> {
    match reals(name, obj, width, row) {
        Ok(v) => Ok(v),
        Err(e) if is_type_error(&e, obj.py()) => {
            let flat: Vec<f64> = if let Ok(v) = obj.extract::<Vec<f64>>() {
                v
            } else if let Ok(v) = obj.extract::<Vec<Vec<f64>>>() {
                v.into_iter().flatten().collect()
            } else {
                return Err(bad(
                    name,
                    "expected an object with the buffer protocol or a sequence of numbers",
                ));
            };
            if !flat.len().is_multiple_of(width) {
                return Err(bad(
                    name,
                    format!(
                        "{} numbers are not a whole number of {row}s of {width}",
                        flat.len()
                    ),
                ));
            }
            if let Some(i) = flat.iter().position(|x| !x.is_finite()) {
                return Err(bad(name, format!("entry {i} is not finite")));
            }
            Ok(flat)
        }
        Err(e) => Err(as_input(e, obj.py())),
    }
}

/// Indices from a buffer, or from a (nested) Python sequence of non-negative integers.
fn indices_any(
    name: &'static str,
    obj: &Bound<'_, PyAny>,
    width: usize,
    row: &'static str,
) -> PyResult<Vec<u32>> {
    match indices(name, obj, width, row) {
        Ok(v) => Ok(v),
        Err(e) if is_type_error(&e, obj.py()) => {
            let flat: Vec<i64> = if let Ok(v) = obj.extract::<Vec<i64>>() {
                v
            } else if let Ok(v) = obj.extract::<Vec<Vec<i64>>>() {
                v.into_iter().flatten().collect()
            } else {
                return Err(bad(
                    name,
                    "expected an object with the buffer protocol or a sequence of integers",
                ));
            };
            if !flat.len().is_multiple_of(width) {
                return Err(bad(
                    name,
                    format!(
                        "{} indices are not a whole number of {row}s of {width}",
                        flat.len()
                    ),
                ));
            }
            flat.iter()
                .enumerate()
                .map(|(i, &x)| {
                    u32::try_from(x)
                        .map_err(|_| bad(name, format!("entry {i} = {x} is not a valid index")))
                })
                .collect()
        }
        Err(e) => Err(as_input(e, obj.py())),
    }
}

type Mesh = (Vec<[f64; 3]>, Vec<[u32; 3]>);

fn mesh(positions: &Bound<'_, PyAny>, faces: &Bound<'_, PyAny>) -> PyResult<Mesh> {
    let p = reals_any("positions", positions, 3, "vertex")?;
    let f = indices_any("faces", faces, 3, "face")?;
    Ok((
        p.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect(),
        f.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect(),
    ))
}

/// An opaque penalty or obstacle that a scene file produced. It goes back into the `penalties`
/// or `ccd_obstacles` option unchanged.
#[pyclass(module = "ddg_tools", frozen)]
struct RepulsiveSurfaceHandle {
    val: Val,
}

#[pymethods]
impl RepulsiveSurfaceHandle {
    fn __repr__(&self) -> String {
        format!("RepulsiveSurfaceHandle({:?})", self.val)
    }
}

/// An option value or a nested value of one.
fn to_val(name: &str, obj: &Bound<'_, PyAny>) -> PyResult<Val> {
    let py = obj.py();
    if obj.is_none() {
        return Ok(Val::None);
    }
    if obj.is_instance_of::<PyBool>() {
        return Ok(Val::Bool(obj.extract()?));
    }
    if obj.is_instance_of::<PyInt>() {
        return obj
            .extract::<i64>()
            .map(Val::Int)
            .map_err(|_| bad(name, "integer out of range"));
    }
    if obj.is_instance_of::<PyFloat>() {
        return Ok(Val::Float(obj.extract()?));
    }
    if obj.is_instance_of::<PyString>() {
        return Ok(Val::Str(obj.extract()?));
    }
    if let Ok(d) = obj.cast::<PyDict>() {
        let mut items = Vec::new();
        for (k, v) in d.iter() {
            let k: String = k
                .extract()
                .map_err(|_| bad(name, "dict keys must be strings"))?;
            let v = to_val(&format!("{name}.{k}"), &v)?;
            items.push((k, v));
        }
        return Ok(Val::Dict(items));
    }
    if obj.is_instance_of::<PyList>() || obj.is_instance_of::<PyTuple>() {
        let mut items = Vec::new();
        for (i, v) in obj.try_iter()?.enumerate() {
            items.push(to_val(&format!("{name}[{i}]"), &v?)?);
        }
        return Ok(Val::List(items));
    }
    if let Ok(h) = obj.extract::<PyRef<'_, RepulsiveSurfaceHandle>>() {
        return Ok(h.val.clone());
    }
    // A buffer: floats, or integers (indices travel as reals); 2-d arrays are read flat.
    for width in [1, 3, 2] {
        if let Ok(v) = reals(leak(name), obj, width, "number") {
            return Ok(Val::Array(v));
        }
    }
    for width in [1, 3, 2] {
        if let Ok(v) = indices(leak(name), obj, width, "number") {
            return Ok(Val::Array(v.into_iter().map(f64::from).collect()));
        }
    }
    // numpy scalars and other numbers.
    if let Ok(i) = obj.extract::<i64>()
        && obj.hasattr("__index__")?
    {
        return Ok(Val::Int(i));
    }
    if let Ok(x) = obj.extract::<f64>() {
        return Ok(Val::Float(x));
    }
    let _ = py;
    Err(bad(
        name,
        format!(
            "unsupported value of type {}",
            obj.get_type().name().map(|n| n.to_string()).unwrap_or_default()
        ),
    ))
}

/// The buffer helpers want a `&'static str` for the name in their messages.
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

fn options(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<(String, Val)>> {
    let Some(kwargs) = kwargs else {
        return Ok(Vec::new());
    };
    let mut out = Vec::new();
    for (k, v) in kwargs.iter() {
        let k: String = k.extract()?;
        let v = to_val(&k, &v)?;
        out.push((k, v));
    }
    Ok(out)
}

// -------------------------------------------------------------------------------------------
// Encoding
// -------------------------------------------------------------------------------------------

fn to_py<'py>(py: Python<'py>, v: &Val) -> PyResult<Bound<'py, PyAny>> {
    Ok(match v {
        Val::None => py.None().into_bound(py),
        Val::Bool(b) => PyBool::new(py, *b).to_owned().into_any(),
        Val::Int(i) => i.into_pyobject(py)?.into_any(),
        Val::Float(x) => x.into_pyobject(py)?.into_any(),
        Val::Str(s) => s.into_pyobject(py)?.into_any(),
        Val::List(items) => {
            let l = PyList::empty(py);
            for i in items {
                l.append(to_py(py, i)?)?;
            }
            l.into_any()
        }
        Val::Array(a) => real_array(py, a)?,
        Val::Dict(items) => {
            let d = PyDict::new(py);
            for (k, i) in items {
                d.set_item(k, to_py(py, i)?)?;
            }
            d.into_any()
        }
        Val::Penalty(_) | Val::Obstacle(_) => Bound::new(py, RepulsiveSurfaceHandle { val: v.clone() })?.into_any(),
    })
}

fn dict<'py>(py: Python<'py>, v: &Val) -> PyResult<Bound<'py, PyDict>> {
    to_py(py, v)?.cast_into::<PyDict>().map_err(PyErr::from)
}

fn flat_points<'py>(py: Python<'py>, p: &[[f64; 3]]) -> PyResult<Bound<'py, PyAny>> {
    real_array(py, p.as_flattened())
}

fn flat_faces<'py>(py: Python<'py>, f: &[[u32; 3]]) -> PyResult<Bound<'py, PyAny>> {
    index_array(py, "faces", f.as_flattened().iter().copied())
}

fn index_vec<'py>(py: Python<'py>, name: &'static str, v: &[u32]) -> PyResult<Bound<'py, PyAny>> {
    index_array(py, name, v.iter().copied())
}

fn _unused(py: Python<'_>) -> PyResult<()> {
    let _ = new_array(py, "i", &[])?;
    Ok(())
}

/// Runs `f` with a progress closure for the library. The closure takes the GIL back, checks
/// for Ctrl-C and calls the Python function with the step's report; `False` stops. An exception
/// stops the run and is returned in `raised`.
fn with_progress<T>(
    py: Python<'_>,
    progress: Option<Py<PyAny>>,
    f: impl FnOnce(Option<&mut dyn FnMut(&Val) -> bool>) -> T + Send,
) -> (T, Option<PyErr>)
where
    T: Send,
{
    let mut raised: Option<PyErr> = None;
    let out = {
        let raised = &mut raised;
        py.detach(move || {
            let mut cb = |report: &Val| -> bool {
                Python::attach(|py| {
                    if let Err(e) = py.check_signals() {
                        *raised = Some(e);
                        return false;
                    }
                    let Some(p) = &progress else { return true };
                    let arg = match to_py(py, report) {
                        Ok(a) => a,
                        Err(e) => {
                            *raised = Some(e);
                            return false;
                        }
                    };
                    match p.call1(py, (arg,)) {
                        Ok(value) => !matches!(value.extract::<bool>(py), Ok(false)),
                        Err(e) => {
                            *raised = Some(e);
                            false
                        }
                    }
                })
            };
            f(Some(&mut cb))
        })
    };
    (out, raised)
}

// -------------------------------------------------------------------------------------------
// The session
// -------------------------------------------------------------------------------------------

/// A repulsive-surface flow that is advanced one step at a time (the engine of a live Blender
/// operator).
///
/// RepulsiveSurfaceSession(positions, faces, **options)
///
/// positions   3 reals per vertex (any buffer, or a sequence of numbers).
/// faces       3 vertex indices per triangle: a manifold triangle mesh (orientation is optional
///             unless a volume constraint or penalty needs it).
///
/// The options are listed in `repulsive_surface_defaults()` and in the .pyi stub; an unknown
/// name is a ValueError. `step()` computes one step with the GIL released and returns a dict.
/// `snapshot()` and the getters read the last completed step and never wait for a running one,
/// so another thread (a timer) may poll while a worker steps. `cancel()` is safe from any
/// thread; a step in progress stops between line-search trials.
///
/// Remeshing (the default) changes the vertex and face counts in every step: always take
/// positions, faces, pinned and face_origin from the same snapshot().
#[pyclass(module = "ddg_tools", frozen)]
struct RepulsiveSurfaceSession {
    engine: Engine,
}

#[pymethods]
impl RepulsiveSurfaceSession {
    #[new]
    #[pyo3(signature = (positions, faces, **options))]
    fn new(
        py: Python<'_>,
        positions: &Bound<'_, PyAny>,
        faces: &Bound<'_, PyAny>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let (p, f) = mesh(positions, faces)?;
        let o = self::options(options)?;
        let engine = py.detach(move || Engine::new(p, f, &o)).map_err(s_err)?;
        Ok(RepulsiveSurfaceSession { engine })
    }

    /// One step. Releases the GIL. Returns a dict: iteration, status ("accepted", "stuck", or
    /// "cancelled" if cancel() was called: nothing is computed then), session_status ("running",
    /// "converged", "max_iterations", "time_limit", "stuck", "cancelled", "error"), n_vertices,
    /// n_faces, energy_before, energy_after, tangent_point_after, penalties_after (list),
    /// gradient_norm, sobolev_norm_sq, tau, tau0, tau_ccd, halvings, newton_steps, projection,
    /// max_violation, dropped_rows, linear (dict), remesh (None or dict: splits, collapses,
    /// flips, ...), seconds (dict of phases).
    fn step<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let report = py.detach(|| self.engine.step()).map_err(s_err)?;
        dict(py, &report)
    }

    /// Steps until the session is final (converged, stuck, max_iterations steps in all),
    /// max_steps steps of this call were taken, cancel() was called, time_limit seconds passed,
    /// or progress returned False. progress(report) is called from the calling thread after every
    /// step with the step's dict (the GIL is released between steps). Returns a dict: end
    /// ("max_steps", "final", "callback", "cancelled", "time_limit"), status, steps_taken,
    /// iteration, energy_before, energy_after, seconds, steps (the list of the step dicts).
    #[pyo3(signature = (max_steps = None, progress = None, time_limit = None))]
    fn run<'py>(
        &self,
        py: Python<'py>,
        max_steps: Option<i64>,
        progress: Option<Bound<'py, PyAny>>,
        time_limit: Option<f64>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let max_steps = max_steps
            .map(|n| usize::try_from(n).map_err(|_| bad("max_steps", format!("{n} must be >= 0"))))
            .transpose()?;
        let progress = progress.map(Bound::unbind);
        let (out, raised) = with_progress(py, progress, |cb| {
            self.engine.run(max_steps, time_limit, cb)
        });
        if let Some(e) = raised {
            return Err(e);
        }
        let summary = out.map_err(s_err)?;
        let snap = self.engine.snapshot().map_err(s_err)?;
        let d = PyDict::new(py);
        d.set_item(
            "end",
            match summary.end {
                RunEnd::MaxSteps
                | RunEnd::Final
                | RunEnd::Callback
                | RunEnd::Cancelled
                | RunEnd::TimeLimit => summary.end.name(),
            },
        )?;
        d.set_item("status", snap.status)?;
        d.set_item("steps_taken", summary.steps.len())?;
        d.set_item("iteration", snap.iteration)?;
        d.set_item("energy_before", summary.energy_before)?;
        d.set_item("energy_after", snap.energy)?;
        d.set_item("seconds", summary.seconds)?;
        d.set_item("steps", to_py(py, &Val::List(summary.steps))?)?;
        Ok(d)
    }

    /// The last completed state as a dict: positions (array('d'), 3 per vertex), faces
    /// (array('i'), 3 per face), pinned (array('i'): the pinned vertices, current ids),
    /// face_origin (array('i'), one per face: the face of the initial mesh it descends from),
    /// iteration, energy, status ("running", "converged", "max_iterations", "time_limit",
    /// "stuck", "cancelled", "error"), error (None or the message), mean_edge_length0, report
    /// (the dict of the last step, or None before the first). Does not wait for a running step.
    fn snapshot<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let s = self.engine.snapshot().map_err(s_err)?;
        let d = PyDict::new(py);
        d.set_item("positions", flat_points(py, &s.positions)?)?;
        d.set_item("faces", flat_faces(py, &s.faces)?)?;
        d.set_item("pinned", index_vec(py, "pinned", &s.pinned)?)?;
        d.set_item("face_origin", index_vec(py, "face_origin", &s.face_origin)?)?;
        d.set_item("iteration", s.iteration)?;
        d.set_item("energy", s.energy)?;
        d.set_item("status", s.status)?;
        d.set_item("error", s.error)?;
        d.set_item("mean_edge_length0", s.mean_edge_length0)?;
        d.set_item(
            "report",
            match &s.last_report {
                Some(v) => to_py(py, v)?,
                None => py.None().into_bound(py),
            },
        )?;
        Ok(d)
    }

    /// Asks the session to stop: a running step ends between line-search trials (the mesh stays
    /// that of the last completed step), later steps do nothing and report status "cancelled".
    /// Safe from any thread; never waits. reset() makes the session usable again.
    fn cancel(&self) -> PyResult<()> {
        self.engine.cancel().map_err(s_err)
    }

    /// Back to the initial mesh with the same options (clears cancel). Waits for a running step.
    fn reset(&self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| self.engine.reset()).map_err(s_err)
    }

    /// array('d'), 3 per vertex: the current positions (last completed step).
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        flat_points(py, &self.engine.snapshot().map_err(s_err)?.positions)
    }

    /// array('i'), 3 per face.
    #[getter]
    fn faces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        flat_faces(py, &self.engine.snapshot().map_err(s_err)?.faces)
    }

    /// The energy (tangent-point energy plus penalties) of the last completed step.
    #[getter]
    fn energy(&self) -> PyResult<f64> {
        Ok(self.engine.snapshot().map_err(s_err)?.energy)
    }

    /// Steps taken.
    #[getter]
    fn iteration(&self) -> PyResult<usize> {
        Ok(self.engine.snapshot().map_err(s_err)?.iteration)
    }

    /// "running", "converged", "max_iterations", "time_limit", "stuck", "cancelled" or "error".
    #[getter]
    fn status(&self) -> PyResult<&'static str> {
        Ok(self.engine.snapshot().map_err(s_err)?.status)
    }

    /// The message of the error when status is "error", else None.
    #[getter]
    fn error(&self) -> PyResult<Option<String>> {
        Ok(self.engine.snapshot().map_err(s_err)?.error)
    }

    /// Number of vertices now.
    #[getter]
    fn n_vertices(&self) -> PyResult<usize> {
        Ok(self.engine.snapshot().map_err(s_err)?.positions.len())
    }

    /// Number of faces now.
    #[getter]
    fn n_faces(&self) -> PyResult<usize> {
        Ok(self.engine.snapshot().map_err(s_err)?.faces.len())
    }

    /// The dict of the last step, or None before the first.
    #[getter]
    fn report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        match self.engine.snapshot().map_err(s_err)?.last_report {
            Some(v) => to_py(py, &v),
            None => Ok(py.None().into_bound(py)),
        }
    }

    fn __repr__(&self) -> PyResult<String> {
        let s = self.engine.snapshot().map_err(s_err)?;
        Ok(format!(
            "RepulsiveSurfaceSession(vertices={}, faces={}, iteration={}, energy={:.6e}, status={:?})",
            s.positions.len(),
            s.faces.len(),
            s.iteration,
            s.energy,
            s.status
        ))
    }
}

// -------------------------------------------------------------------------------------------
// Functions
// -------------------------------------------------------------------------------------------

/// The repulsive-surface flow in one call. Returns (positions, faces, report): positions an
/// array('d') of 3 per vertex and faces an array('i') of 3 per face (both of the final mesh,
/// which differs from the input with remeshing), and a dict: status, iterations,
/// energy_initial, energy_final, seconds, pinned and face_origin (array('i'), as in
/// RepulsiveSurfaceSession.snapshot), steps (the list of the step dicts).
///
/// steps     the most steps (default: the option max_iterations, 1000).
/// progress  progress(step_report) after every step, called from the calling thread; False stops
///           the flow (status "cancelled"). An exception (and Ctrl-C) stops it and propagates.
/// options   as for RepulsiveSurfaceSession.
///
/// A RepulsiveSurfaceSession stepped the same number of times gives bitwise the same mesh.
#[pyfunction]
#[pyo3(signature = (positions, faces, steps = None, progress = None, **options))]
fn repulsive_surface_flow<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    faces: &Bound<'py, PyAny>,
    steps: Option<i64>,
    progress: Option<Bound<'py, PyAny>>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let (p, f) = mesh(positions, faces)?;
    let mut o = self::options(options)?;
    if let Some(n) = steps {
        let n = usize::try_from(n).map_err(|_| bad("steps", format!("{n} must be >= 0")))?;
        o.retain(|(k, _)| k != "max_iterations");
        o.push(("max_iterations".into(), Val::Int(n as i64)));
    }
    let progress = progress.map(Bound::unbind);
    let (out, raised) = with_progress(py, progress, |cb| {
        crate::repulsive_surfaces::flow_once(p, f, &o, cb)
    });
    if let Some(e) = raised {
        return Err(e);
    }
    let out = out.map_err(s_err)?;
    let report = dict(py, &out.report)?;
    report.set_item("pinned", index_vec(py, "pinned", &out.pinned)?)?;
    report.set_item("face_origin", index_vec(py, "face_origin", &out.face_origin)?)?;
    PyTuple::new(
        py,
        [
            flat_points(py, &out.positions)?,
            flat_faces(py, &out.faces)?,
            report.into_any(),
        ],
    )
}

/// The tangent-point energy of a mesh plus the penalties of the options, and the differential
/// of the total: (energy, differential, report). energy is the total, differential an
/// array('d') of 3 per vertex; report is a dict: energy (tangent-point part), total, penalties
/// (a list of {name, value} in the order given), subcritical, pairs_exact, nodes_lumped,
/// seconds. The options read are alpha, beta, energy_eval ("barnes_hut" with theta and
/// leaf_size, or "exact"), penalties and threads; the others are validated and ignored.
#[pyfunction]
#[pyo3(signature = (positions, faces, **options))]
fn surface_energy<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    faces: &Bound<'py, PyAny>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<(f64, Bound<'py, PyAny>, Bound<'py, PyDict>)> {
    let (p, f) = mesh(positions, faces)?;
    let o = self::options(options)?;
    let (differential, report) = py
        .detach(move || crate::repulsive_surfaces::surface_energy(&p, f, &o))
        .map_err(s_err)?;
    let total = match report.get("total") {
        Some(Val::Float(x)) => *x,
        _ => f64::NAN,
    };
    Ok((total, flat_points(py, &differential)?, dict(py, &report)?))
}

/// The defaults of all options of RepulsiveSurfaceSession and repulsive_surface_flow, as a dict.
/// They are the defaults of the Rust library, read from it; the .pyi stub and the tests check
/// that they agree. Options that are absent by default (penalties, pins, ...) are None. Passing
/// the whole dict back as keywords gives the default flow.
#[pyfunction]
fn repulsive_surface_defaults(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    dict(py, &crate::repulsive_surfaces::defaults())
}

/// What this wheel and machine offer for the `device` option of the Repulsive Surfaces flow, as
/// a dict: `cuda_compiled` (the wheel was built with CUDA support), `available` (a CUDA device
/// can be opened), `devices` (a list of dicts `ordinal`, `name`, `compute_capability`,
/// `memory_bytes`) and `reason` (why no device is usable, or None). Never raises.
#[pyfunction]
fn gpu_info(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    let v = py.detach(crate::repulsive_surfaces::gpu_info);
    dict(py, &v)
}

/// Loads a scene file of the reference program (repel_mesh, constrain, minimize, obstacle,
/// implicit, method, ...) and returns (positions, faces, options): the mesh and the keyword
/// options that run the scene, as a dict to pass as `**options`. The scene's penalties and
/// obstacles are opaque RepulsiveSurfaceHandle objects inside the dict. Paths in the scene are
/// relative to its directory. ValueError for a file that cannot be read or parsed and for
/// normal constraints, which the library does not support.
#[pyfunction]
fn repulsive_surface_scene<'py>(
    py: Python<'py>,
    path: &str,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>, Bound<'py, PyDict>)> {
    let path = path.to_owned();
    let (p, f, o) = py
        .detach(move || crate::repulsive_surfaces::scene(&path))
        .map_err(s_err)?;
    let _ = make_surface;
    Ok((flat_points(py, &p)?, flat_faces(py, &f)?, dict(py, &o)?))
}
