//! The PyO3 glue of `ddg_tools.RepulsiveSession`, `repulsive_flow` and `curve_energy`.
//!
//! As in the parent module no logic lives here: buffers (or plain Python sequences) are decoded,
//! [`crate::repulsive`] does the work with the GIL released (taking it back only to call the
//! progress function and to check for Ctrl-C), results go back as `array.array` objects and
//! plain dictionaries. Errors in the arguments or options are `RepulsiveInputError`, which is
//! both a `ddg_tools.InputError` and a `ValueError`; a failed computation is `SolveError`.

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyBool, PyDict, PyList, PyTuple, PyType};

use ddg_repulsive_curves::{Axis, BhOptions};

use super::{InputError, SolveError, index_array, indices, new_array, real_array, reals};
use crate::repulsive::{
    Engine, FieldSpec, LengthMode, OPTION_NAMES, RepError, RepulsiveOptions, RunEnd, SurfaceSpec,
    Val,
};

static INPUT_ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    // class RepulsiveInputError(InputError, ValueError)
    let bases = PyTuple::new(
        py,
        [py.get_type::<InputError>(), py.get_type::<PyValueError>()],
    )?;
    let dict = PyDict::new(py);
    dict.set_item("__module__", "ddg_tools")?;
    dict.set_item(
        "__doc__",
        "An argument or an option of the repulsive-curves functions is invalid. Both a \
         ddg_tools.InputError and a ValueError; the message is the library's.",
    )?;
    let cls = py
        .import("builtins")?
        .getattr("type")?
        .call1(("RepulsiveInputError", bases, dict))?
        .cast_into::<PyType>()
        .map_err(PyErr::from)?;
    m.add("RepulsiveInputError", &cls)?;
    let _ = INPUT_ERROR.set(py, cls.unbind());
    m.add_class::<RepulsiveSession>()?;
    m.add_function(wrap_pyfunction!(repulsive_flow, m)?)?;
    m.add_function(wrap_pyfunction!(curve_energy, m)?)?;
    m.add_function(wrap_pyfunction!(repulsive_defaults, m)?)?;
    m.add("REPULSIVE_BACKENDS", PyList::new(py, REPULSIVE_BACKENDS)?)?;
    Ok(())
}

/// The values of the `backend` option.
const REPULSIVE_BACKENDS: [&str; 5] = ["auto", "dense", "barnes_hut", "multigrid", "accelerated"];

fn input_err(message: String) -> PyErr {
    Python::attach(|py| match INPUT_ERROR.get(py) {
        Some(t) => PyErr::from_type(t.bind(py).clone(), message),
        None => InputError::new_err(message),
    })
}

fn rep_err(e: RepError) -> PyErr {
    match e {
        RepError::Input(m) => input_err(m),
        RepError::Solve(m) => SolveError::new_err(m),
    }
}

// -------------------------------------------------------------------------------------------
// Decoding
// -------------------------------------------------------------------------------------------

fn bad(name: &str, message: impl std::fmt::Display) -> PyErr {
    input_err(format!("{name}: {message}"))
}

/// A buffer-decoding `InputError` becomes a `RepulsiveInputError` (also a `ValueError`).
fn as_input(e: PyErr, py: Python<'_>) -> PyErr {
    if e.is_instance_of::<InputError>(py) {
        input_err(e.value(py).to_string())
    } else {
        e
    }
}

fn is_type_error(e: &PyErr, py: Python<'_>) -> bool {
    e.is_instance_of::<PyTypeError>(py)
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

fn points(v: Vec<f64>) -> Vec<[f64; 3]> {
    v.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

fn pairs(v: Vec<u32>) -> Vec<[u32; 2]> {
    v.chunks_exact(2).map(|c| [c[0], c[1]]).collect()
}

fn triples(v: Vec<u32>) -> Vec<[u32; 3]> {
    v.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

fn float(name: &str, v: &Bound<'_, PyAny>) -> PyResult<f64> {
    if v.is_instance_of::<PyBool>() {
        return Err(bad(name, "expected a number, got a bool"));
    }
    v.extract::<f64>()
        .map_err(|_| bad(name, format!("expected a number, got {}", type_name(v))))
}

fn count(name: &str, v: &Bound<'_, PyAny>) -> PyResult<usize> {
    let x = v
        .extract::<i64>()
        .map_err(|_| bad(name, format!("expected an integer, got {}", type_name(v))))?;
    usize::try_from(x).map_err(|_| bad(name, format!("{x} must be >= 0")))
}

fn text(name: &str, v: &Bound<'_, PyAny>) -> PyResult<String> {
    v.extract::<String>()
        .map_err(|_| bad(name, format!("expected a string, got {}", type_name(v))))
}

fn flag(name: &str, v: &Bound<'_, PyAny>) -> PyResult<bool> {
    v.extract::<bool>()
        .map_err(|_| bad(name, format!("expected a bool, got {}", type_name(v))))
}

fn vec3(name: &str, v: &Bound<'_, PyAny>) -> PyResult<[f64; 3]> {
    let x: Vec<f64> = v.extract().map_err(|_| bad(name, "expected 3 numbers"))?;
    <[f64; 3]>::try_from(x).map_err(|x| bad(name, format!("expected 3 numbers, got {}", x.len())))
}

fn type_name(v: &Bound<'_, PyAny>) -> String {
    v.get_type()
        .name()
        .map(|n| n.to_string())
        .unwrap_or_else(|_| "?".into())
}

fn index_list(name: &'static str, v: &Bound<'_, PyAny>) -> PyResult<Vec<u32>> {
    indices_any(name, v, 1, "index")
}

/// The entries of a dict that names a surface or a field: a strict reader.
struct Spec<'py> {
    what: String,
    dict: Bound<'py, PyDict>,
    used: Vec<String>,
}

impl<'py> Spec<'py> {
    fn new(what: &str, obj: &Bound<'py, PyAny>) -> PyResult<Spec<'py>> {
        let dict = obj
            .cast::<PyDict>()
            .map_err(|_| {
                bad(
                    what,
                    format!("expected a dict with a 'type', got {}", type_name(obj)),
                )
            })?
            .clone();
        Ok(Spec {
            what: what.to_owned(),
            dict,
            used: Vec::new(),
        })
    }

    fn get(&mut self, key: &str) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.used.push(key.to_owned());
        Ok(self.dict.get_item(key)?.filter(|v| !v.is_none()))
    }

    fn require(&mut self, key: &str) -> PyResult<Bound<'py, PyAny>> {
        self.get(key)?
            .ok_or_else(|| bad(&self.what, format!("missing '{key}'")))
    }

    fn finish(&self) -> PyResult<()> {
        for (k, _) in self.dict.iter() {
            let k = k.to_string();
            if k != "type" && !self.used.contains(&k) {
                return Err(bad(&self.what, format!("unknown key '{k}'")));
            }
        }
        Ok(())
    }

    fn kind(&self) -> PyResult<String> {
        match self.dict.get_item("type")? {
            Some(t) => text(&self.what, &t),
            None => Err(bad(&self.what, "missing 'type'")),
        }
    }
}

const SURFACES: &str =
    "sphere, torus, plane, double_torus, rounded_box, metaballs, union, smooth_union";

fn surface_spec(what: &str, obj: &Bound<'_, PyAny>) -> PyResult<SurfaceSpec> {
    let mut s = Spec::new(what, obj)?;
    let kind = s.kind()?;
    let v3 = |s: &mut Spec<'_>, key: &str, default: Option<[f64; 3]>| -> PyResult<[f64; 3]> {
        match s.get(key)? {
            Some(v) => vec3(&format!("{what}.{key}"), &v),
            None => default.ok_or_else(|| bad(what, format!("missing '{key}'"))),
        }
    };
    let num = |s: &mut Spec<'_>, key: &str, default: Option<f64>| -> PyResult<f64> {
        match s.get(key)? {
            Some(v) => float(&format!("{what}.{key}"), &v),
            None => default.ok_or_else(|| bad(what, format!("missing '{key}'"))),
        }
    };
    let spec = match kind.as_str() {
        "sphere" => SurfaceSpec::Sphere {
            center: v3(&mut s, "center", Some([0.0; 3]))?,
            radius: num(&mut s, "radius", Some(1.0))?,
        },
        "torus" => SurfaceSpec::Torus {
            center: v3(&mut s, "center", Some([0.0; 3]))?,
            major: num(&mut s, "major", None)?,
            minor: num(&mut s, "minor", None)?,
        },
        "plane" => SurfaceSpec::Plane {
            point: v3(&mut s, "point", Some([0.0; 3]))?,
            normal: v3(&mut s, "normal", None)?,
        },
        "double_torus" => SurfaceSpec::DoubleTorus {
            radius2: num(&mut s, "radius2", None)?,
        },
        "rounded_box" => SurfaceSpec::RoundedBox {
            center: v3(&mut s, "center", Some([0.0; 3]))?,
            half: v3(&mut s, "half", None)?,
            radius: num(&mut s, "radius", Some(0.0))?,
        },
        "metaballs" => {
            let balls = s.require("balls")?;
            let list: Vec<Bound<'_, PyAny>> = balls
                .extract()
                .map_err(|_| bad(what, "'balls' must be a list of (center, amplitude, width)"))?;
            let mut out = Vec::new();
            for b in list {
                let t: (Bound<'_, PyAny>, f64, f64) = b
                    .extract()
                    .map_err(|_| bad(what, "each ball is (center, amplitude, width)"))?;
                out.push((vec3(&format!("{what}.balls"), &t.0)?, t.1, t.2));
            }
            SurfaceSpec::Metaballs {
                balls: out,
                threshold: num(&mut s, "threshold", None)?,
            }
        }
        "union" => {
            let parts = s.require("parts")?;
            let list: Vec<Bound<'_, PyAny>> = parts
                .extract()
                .map_err(|_| bad(what, "'parts' must be a list of surface dicts"))?;
            SurfaceSpec::Union(
                list.iter()
                    .map(|p| surface_spec(&format!("{what}.parts"), p))
                    .collect::<PyResult<_>>()?,
            )
        }
        "smooth_union" => {
            let a = s.require("a")?;
            let b = s.require("b")?;
            SurfaceSpec::SmoothUnion {
                a: Box::new(surface_spec(&format!("{what}.a"), &a)?),
                b: Box::new(surface_spec(&format!("{what}.b"), &b)?),
                k: num(&mut s, "k", None)?,
            }
        }
        other => {
            return Err(bad(
                what,
                format!("unknown surface type '{other}'; use one of {SURFACES}"),
            ));
        }
    };
    s.finish()?;
    Ok(spec)
}

fn field_spec(obj: &Bound<'_, PyAny>) -> PyResult<FieldSpec> {
    let mut s = Spec::new("field", obj)?;
    let spec = match s.kind()?.as_str() {
        "constant" => {
            let vector = vec3("field.vector", &s.require("vector")?)?;
            let normalize = match s.get("normalize")? {
                Some(v) => flag("field.normalize", &v)?,
                None => false,
            };
            FieldSpec::Constant { vector, normalize }
        }
        "circular" => FieldSpec::Circular,
        "sphere_tangent" => FieldSpec::SphereTangent,
        other => {
            return Err(bad(
                "field",
                format!("unknown field type '{other}'; use constant, circular or sphere_tangent"),
            ));
        }
    };
    s.finish()?;
    Ok(spec)
}

/// Fills the options from the keyword arguments; an unknown name is an error.
fn parse_options(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<RepulsiveOptions> {
    let mut o = RepulsiveOptions::default();
    let Some(kwargs) = kwargs else { return Ok(o) };
    for (key, v) in kwargs.iter() {
        let key: String = key.extract()?;
        let k = key.as_str();
        let none = v.is_none();
        match k {
            "alpha" => o.alpha = float(k, &v)?,
            "beta" => o.beta = float(k, &v)?,
            "metric" => o.metric = text(k, &v)?,
            "quadrature" => o.quadrature = text(k, &v)?,
            "length" => {
                o.length = match text(k, &v)?.as_str() {
                    "edges" => LengthMode::Edges,
                    "total" => LengthMode::Total,
                    "free" => LengthMode::Free,
                    other => {
                        return Err(bad(
                            k,
                            format!("'{other}' is not valid; use 'edges', 'total' or 'free'"),
                        ));
                    }
                }
            }
            "grow_factor" => o.grow_factor = float(k, &v)?,
            "grow_rate" => o.grow_rate = float(k, &v)?,
            "max_subdivisions" => {
                o.max_subdivisions =
                    u32::try_from(count(k, &v)?).map_err(|_| bad(k, "too large"))?;
            }
            "barycenter" => o.barycenter = if none { None } else { Some(flag(k, &v)?) },
            "pinned" => {
                o.pinned = if none {
                    Vec::new()
                } else {
                    index_list("pinned", &v)?
                }
            }
            "tangent_vertices" => {
                o.tangent_vertices = if none {
                    Vec::new()
                } else {
                    index_list("tangent_vertices", &v)?
                }
            }
            "tangent_edges" => {
                o.tangent_edges = if none {
                    Vec::new()
                } else {
                    index_list("tangent_edges", &v)?
                }
            }
            "surface" => {
                o.surface = if none {
                    None
                } else {
                    let vertices = match v.cast::<PyDict>() {
                        Ok(d) => match d.get_item("vertices")? {
                            Some(vs) if !vs.is_none() => Some(index_list("surface.vertices", &vs)?),
                            _ => None,
                        },
                        Err(_) => None,
                    };
                    // `vertices` is not a surface parameter.
                    let stripped = match v.cast::<PyDict>() {
                        Ok(d) => {
                            let c = d.copy()?;
                            c.del_item("vertices").ok();
                            c.into_any()
                        }
                        Err(_) => v.clone(),
                    };
                    Some((surface_spec("surface", &stripped)?, vertices))
                }
            }
            "planar_axis" => {
                o.planar_axis = if none {
                    None
                } else {
                    Some(match text(k, &v)?.as_str() {
                        "x" => Axis::X,
                        "y" => Axis::Y,
                        "z" => Axis::Z,
                        other => {
                            return Err(bad(
                                k,
                                format!("'{other}' is not valid; use 'x', 'y' or 'z'"),
                            ));
                        }
                    })
                }
            }
            "obstacles" => {
                o.obstacles = Vec::new();
                if !none {
                    let list: Vec<(Bound<'_, PyAny>, Bound<'_, PyAny>)> = v
                        .extract()
                        .map_err(|_| bad(k, "expected a list of (vertices, triangles) pairs"))?;
                    for (p, t) in list {
                        o.obstacles.push((
                            points(reals_any("obstacles.vertices", &p, 3, "vertex")?),
                            triples(indices_any("obstacles.triangles", &t, 3, "triangle")?),
                        ));
                    }
                }
            }
            "obstacle_weight" => o.obstacle_weight = float(k, &v)?,
            "obstacle_exponent" => {
                o.obstacle_exponent = if none { None } else { Some(float(k, &v)?) }
            }
            "obstacle_theta" => o.obstacle_theta = if none { None } else { Some(float(k, &v)?) },
            "total_length_weight" => o.total_length_weight = float(k, &v)?,
            "length_difference_weight" => o.length_difference_weight = float(k, &v)?,
            "field" => o.field = if none { None } else { Some(field_spec(&v)?) },
            "field_weight" => o.field_weight = float(k, &v)?,
            "line_search" => o.line_search = text(k, &v)?,
            "line_search_start" => o.line_search_start = text(k, &v)?,
            "line_search_start_factor" => o.line_search_start_factor = float(k, &v)?,
            "armijo" => o.armijo = float(k, &v)?,
            "shrink" => o.shrink = float(k, &v)?,
            "max_halvings" => {
                o.max_halvings = u32::try_from(count(k, &v)?).map_err(|_| bad(k, "too large"))?;
            }
            "gap" => o.gap = float(k, &v)?,
            "backend" => o.backend = text(k, &v)?,
            "auto_dense_max_vertices" => o.auto_dense_max_vertices = count(k, &v)?,
            "dense_fallback_max_vertices" => o.dense_fallback_max_vertices = count(k, &v)?,
            "mg_coarsest_vertices" => o.mg_coarsest_vertices = count(k, &v)?,
            "mg_max_levels" => o.mg_max_levels = count(k, &v)?,
            "mg_smoothing_steps" => o.mg_smoothing_steps = count(k, &v)?,
            "mg_initial_smoothing_steps" => o.mg_initial_smoothing_steps = count(k, &v)?,
            "mg_tolerance" => o.mg_tolerance = float(k, &v)?,
            "mg_max_cycles" => o.mg_max_cycles = count(k, &v)?,
            "mg_flexible_cg" => o.mg_flexible_cg = flag(k, &v)?,
            "mg_near_field_factor" => o.mg_near_field_factor = float(k, &v)?,
            "mg_smoother_factor" => o.mg_smoother_factor = float(k, &v)?,
            "mg_max_fixed_fraction" => o.mg_max_fixed_fraction = float(k, &v)?,
            "mg_energy_check" => o.mg_energy_check = flag(k, &v)?,
            "mg_failure_residual" => o.mg_failure_residual = float(k, &v)?,
            "bh_theta_x" => o.bh_theta_x = float(k, &v)?,
            "bh_theta_t" => o.bh_theta_t = float(k, &v)?,
            "bh_leaf_size" => o.bh_leaf_size = count(k, &v)?,
            "threads" => o.threads = if none { None } else { Some(count(k, &v)?) },
            "dense_threads" => o.dense_threads = if none { None } else { Some(count(k, &v)?) },
            "projection_tolerance" => o.projection_tolerance = float(k, &v)?,
            "projection_absolute" => o.projection_absolute = flag(k, &v)?,
            "projection_iterations" => {
                o.projection_iterations =
                    u32::try_from(count(k, &v)?).map_err(|_| bad(k, "too large"))?;
            }
            "stop" => o.stop = text(k, &v)?,
            "stop_threshold" => o.stop_threshold = if none { None } else { Some(float(k, &v)?) },
            "stall_window" => o.stall_window = count(k, &v)?,
            "gradient_norm" => o.gradient_norm = text(k, &v)?,
            "max_steps" => o.max_steps = count(k, &v)?,
            "time_limit" => o.time_limit = if none { None } else { Some(float(k, &v)?) },
            other => {
                return Err(input_err(format!(
                    "unknown option '{other}'; the options are {}",
                    OPTION_NAMES.join(", ")
                )));
            }
        }
    }
    Ok(o)
}

type Curve = (Vec<[f64; 3]>, Vec<[u32; 2]>);

fn curve(positions: &Bound<'_, PyAny>, edges: &Bound<'_, PyAny>) -> PyResult<Curve> {
    Ok((
        points(reals_any("positions", positions, 3, "vertex")?),
        pairs(indices_any("edges", edges, 2, "edge")?),
    ))
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
        Val::Dict(items) => {
            let d = PyDict::new(py);
            for (k, i) in items {
                d.set_item(k, to_py(py, i)?)?;
            }
            d.into_any()
        }
    })
}

fn dict<'py>(py: Python<'py>, v: &Val) -> PyResult<Bound<'py, PyDict>> {
    to_py(py, v)?.cast_into::<PyDict>().map_err(PyErr::from)
}

fn flat_points<'py>(py: Python<'py>, p: &[[f64; 3]]) -> PyResult<Bound<'py, PyAny>> {
    real_array(py, p.as_flattened())
}

fn flat_edges<'py>(py: Python<'py>, e: &[[u32; 2]]) -> PyResult<Bound<'py, PyAny>> {
    index_array(py, "edges", e.as_flattened().iter().copied())
}

fn signed_array<'py>(py: Python<'py>, values: &[i32]) -> PyResult<Bound<'py, PyAny>> {
    let bytes: Vec<u8> = values.iter().flat_map(|x| x.to_ne_bytes()).collect();
    new_array(py, "i", &bytes)
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

/// A repulsive-curves flow that is advanced one step at a time (the engine of a live Blender
/// operator).
///
/// RepulsiveSession(positions, edges, **options)
///
/// positions   3 reals per vertex (any buffer, or a sequence of numbers).
/// edges       2 vertex indices per edge. A closed curve is a cycle; junctions are allowed.
///
/// The options are listed in `repulsive_defaults()` and in the .pyi stub; an unknown name is a
/// ValueError. `step()` computes one step with the GIL released and returns a dict. `snapshot()`,
/// `positions`, `edges`, `energy`, `status` ... read the last completed step and never wait for
/// a running one, so another thread (a timer) may poll while a worker steps. `cancel()` is safe
/// from any thread; a step in progress stops between line-search trials.
///
/// Subdivision (option max_subdivisions): when it happens the vertex count changes (see
/// `n_vertices` or the report key `subdivided`). Old vertices keep their indices, the midpoint of
/// edge I is vertex |V| + I and edge I becomes edges I and |E| + I. `vertex_source` and
/// `edge_source` give, for every vertex and edge, the input vertex or edge it is or descends
/// from (-1 for a created vertex).
#[pyclass(module = "ddg_tools", frozen)]
struct RepulsiveSession {
    engine: Engine,
}

#[pymethods]
impl RepulsiveSession {
    #[new]
    #[pyo3(signature = (positions, edges, **options))]
    fn new(
        py: Python<'_>,
        positions: &Bound<'_, PyAny>,
        edges: &Bound<'_, PyAny>,
        options: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let (p, e) = curve(positions, edges)?;
        let o = parse_options(options)?;
        let engine = py.detach(move || Engine::new(p, e, &o)).map_err(rep_err)?;
        Ok(RepulsiveSession { engine })
    }

    /// One step. Releases the GIL. Returns a dict: iteration, status ("accepted", "stuck",
    /// "converged", "retarget_stalled", or "cancelled" if cancel() was called: nothing is
    /// computed then), session_status, n_vertices, n_edges, energy_before, energy_after,
    /// tangent_point_after, penalties_after, gradient_norm, tau, tau_max, halvings,
    /// max_violation, projection_iterations, linear (dict), retarget (None or dict), subdivided,
    /// energy_backend, solver_backend, seconds (dict of phases).
    fn step<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let report = py.detach(|| self.engine.step()).map_err(rep_err)?;
        dict(py, &report)
    }

    /// Steps until the session is final (converged, stuck, max_steps), max_steps steps of this
    /// call were taken, cancel() was called, time_limit seconds passed, or progress returned
    /// False. progress(report) is called from the calling thread after every step with the
    /// step's dict (the GIL is released between steps). Returns a dict: end ("max_steps",
    /// "final", "callback", "cancelled", "time_limit"), status, steps_taken, iteration,
    /// energy_before, energy_after, seconds, steps (the list of the step dicts).
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
        let summary = out.map_err(rep_err)?;
        let snap = self.engine.snapshot().map_err(rep_err)?;
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

    /// The last completed state: (positions, edges, iteration, energy, status). positions is an
    /// array('d') of 3 per vertex, edges an array('i') of 2 per edge, status one of "running",
    /// "converged", "max_steps", "time_limit", "stuck", "cancelled", "error". Does not wait for
    /// a running step.
    fn snapshot<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let s = self.engine.snapshot().map_err(rep_err)?;
        PyTuple::new(
            py,
            [
                flat_points(py, &s.positions)?,
                flat_edges(py, &s.edges)?,
                s.iteration.into_pyobject(py)?.into_any(),
                s.energy.into_pyobject(py)?.into_any(),
                s.status.into_pyobject(py)?.into_any(),
            ],
        )
    }

    /// Asks the session to stop: a running step ends between line-search trials (the curve stays
    /// that of the last completed step), later steps do nothing and report status "cancelled".
    /// Safe from any thread; never waits. reset() makes the session usable again.
    fn cancel(&self) -> PyResult<()> {
        self.engine.cancel().map_err(rep_err)
    }

    /// Back to the initial curve with the same options (clears cancel). Waits for a running step.
    fn reset(&self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| self.engine.reset()).map_err(rep_err)
    }

    /// Changes the length schedule: the lengths move to `factor` times the current targets, by
    /// `rate` of the current total per step. factor = 1 removes the schedule. ValueError if
    /// there is no length constraint to move (length="free").
    #[pyo3(signature = (factor, rate = 0.01))]
    fn set_schedule(&self, py: Python<'_>, factor: f64, rate: f64) -> PyResult<()> {
        py.detach(|| self.engine.set_schedule(factor, rate))
            .map_err(rep_err)
    }

    /// array('d'), 3 per vertex: the current curve (last completed step).
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        flat_points(py, &self.engine.snapshot().map_err(rep_err)?.positions)
    }

    /// array('i'), 2 per edge.
    #[getter]
    fn edges<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        flat_edges(py, &self.engine.snapshot().map_err(rep_err)?.edges)
    }

    /// The energy (tangent-point energy plus penalties) of the last completed step.
    #[getter]
    fn energy(&self) -> PyResult<f64> {
        Ok(self.engine.snapshot().map_err(rep_err)?.energy)
    }

    /// Steps taken.
    #[getter]
    fn iteration(&self) -> PyResult<usize> {
        Ok(self.engine.snapshot().map_err(rep_err)?.iteration)
    }

    /// "running", "converged", "max_steps", "time_limit", "stuck", "cancelled" or "error".
    #[getter]
    fn status(&self) -> PyResult<&'static str> {
        Ok(self.engine.snapshot().map_err(rep_err)?.status)
    }

    /// The message of the error when status is "error", else None.
    #[getter]
    fn error(&self) -> PyResult<Option<String>> {
        Ok(self.engine.snapshot().map_err(rep_err)?.error)
    }

    /// Number of vertices now (grows with subdivision).
    #[getter]
    fn n_vertices(&self) -> PyResult<usize> {
        Ok(self.engine.snapshot().map_err(rep_err)?.positions.len())
    }

    /// Number of edges now.
    #[getter]
    fn n_edges(&self) -> PyResult<usize> {
        Ok(self.engine.snapshot().map_err(rep_err)?.edges.len())
    }

    /// Number of subdivisions so far.
    #[getter]
    fn subdivisions(&self) -> PyResult<u32> {
        Ok(self.engine.snapshot().map_err(rep_err)?.subdivisions)
    }

    /// Mean edge length of the start; the unit of the relative tolerances.
    #[getter]
    fn mean_edge_length0(&self) -> PyResult<f64> {
        Ok(self.engine.snapshot().map_err(rep_err)?.mean_edge_length0)
    }

    /// Total length of the curve now.
    #[getter]
    fn total_length(&self) -> PyResult<f64> {
        Ok(self.engine.snapshot().map_err(rep_err)?.total_length)
    }

    /// array('i'), one per vertex: the input vertex it is, -1 for a vertex created by subdivision.
    #[getter]
    fn vertex_source<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        signed_array(
            py,
            &self
                .engine
                .snapshot()
                .map_err(rep_err)?
                .sources
                .vertex_source,
        )
    }

    /// array('i'), one per edge: the input edge it descends from.
    #[getter]
    fn edge_source<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let s = self.engine.snapshot().map_err(rep_err)?;
        index_array(py, "edge_source", s.sources.edge_source.iter().copied())
    }

    /// The dict of the last step, or None before the first.
    #[getter]
    fn report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        match self.engine.snapshot().map_err(rep_err)?.last_report {
            Some(v) => to_py(py, &v),
            None => Ok(py.None().into_bound(py)),
        }
    }

    fn __repr__(&self) -> PyResult<String> {
        let s = self.engine.snapshot().map_err(rep_err)?;
        Ok(format!(
            "RepulsiveSession(vertices={}, edges={}, iteration={}, energy={:.6e}, status={:?})",
            s.positions.len(),
            s.edges.len(),
            s.iteration,
            s.energy,
            s.status
        ))
    }
}

// -------------------------------------------------------------------------------------------
// Functions
// -------------------------------------------------------------------------------------------

/// The repulsive flow in one call. Returns (positions, edges, report): positions an array('d') of
/// 3 per vertex and edges an array('i') of 2 per edge (both of the final curve, which has more
/// vertices after a subdivision), and a dict: status, iterations, energy_initial, energy_final,
/// seconds, n_subdivisions, vertex_source, edge_source (lists, as in RepulsiveSession), steps (the
/// list of the step dicts of RepulsiveSession.step).
///
/// steps     the most steps (default: the option max_steps, 1000).
/// progress  progress(step_report) after every step, called from the calling thread; False stops
///           the flow (status "cancelled"). An exception (and Ctrl-C) stops it and propagates.
/// options   as for RepulsiveSession.
///
/// A RepulsiveSession stepped the same number of times gives bitwise the same curve.
#[pyfunction]
#[pyo3(signature = (positions, edges, steps = None, progress = None, **options))]
fn repulsive_flow<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    edges: &Bound<'py, PyAny>,
    steps: Option<i64>,
    progress: Option<Bound<'py, PyAny>>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let (p, e) = curve(positions, edges)?;
    let mut o = parse_options(options)?;
    if let Some(n) = steps {
        o.max_steps = usize::try_from(n).map_err(|_| bad("steps", format!("{n} must be >= 0")))?;
    }
    let progress = progress.map(Bound::unbind);
    let (out, raised) = with_progress(py, progress, |cb| crate::repulsive::flow_once(p, e, &o, cb));
    if let Some(e) = raised {
        return Err(e);
    }
    let out = out.map_err(rep_err)?;
    PyTuple::new(
        py,
        [
            flat_points(py, &out.positions)?,
            flat_edges(py, &out.edges)?,
            to_py(py, &out.report)?,
        ],
    )
}

/// The tangent-point energy of a curve and its differential: (energy, differential), the
/// differential an array('d') of 3 per vertex (dE/dx_v).
///
/// backend   "exact" (all pairs, default) or "barnes_hut" (with bh_theta_x, bh_theta_t,
///           bh_leaf_size).
#[pyfunction]
#[pyo3(signature = (
    positions,
    edges,
    alpha = 3.0,
    beta = 6.0,
    *,
    backend = "exact",
    bh_theta_x = 0.125,
    bh_theta_t = 0.25,
    bh_leaf_size = 8,
    threads = None,
))]
#[allow(clippy::too_many_arguments)]
fn curve_energy<'py>(
    py: Python<'py>,
    positions: &Bound<'py, PyAny>,
    edges: &Bound<'py, PyAny>,
    alpha: f64,
    beta: f64,
    backend: &str,
    bh_theta_x: f64,
    bh_theta_t: f64,
    bh_leaf_size: i64,
    threads: Option<i64>,
) -> PyResult<(f64, Bound<'py, PyAny>)> {
    let (p, e) = curve(positions, edges)?;
    let bh = match backend {
        "exact" => None,
        "barnes_hut" => Some(BhOptions {
            theta_x: bh_theta_x,
            theta_t: bh_theta_t,
            leaf_size: usize::try_from(bh_leaf_size)
                .map_err(|_| bad("bh_leaf_size", format!("{bh_leaf_size} must be >= 0")))?,
        }),
        other => {
            return Err(bad(
                "backend",
                format!("'{other}' is not valid; use 'exact' or 'barnes_hut'"),
            ));
        }
    };
    let threads = threads
        .map(|n| usize::try_from(n).map_err(|_| bad("threads", format!("{n} must be >= 1"))))
        .transpose()?;
    let (energy, differential) = py
        .detach(move || crate::repulsive::curve_energy(&p, e, alpha, beta, bh, threads))
        .map_err(rep_err)?;
    Ok((energy, flat_points(py, &differential)?))
}

/// The defaults of all options of RepulsiveSession and repulsive_flow, as a dict. They are the
/// defaults of the Rust library (SPEC.md section 10), read from it; the .pyi stub and the tests
/// check that they agree. List-valued and absent options (pinned, obstacles, surface, ...) are
/// None; stop_threshold None means the rule's own default (1e-4 gradient, 1e-3 relative).
#[pyfunction]
fn repulsive_defaults(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    let o = RepulsiveOptions::default();
    let d = PyDict::new(py);
    d.set_item("alpha", o.alpha)?;
    d.set_item("beta", o.beta)?;
    d.set_item("metric", &o.metric)?;
    d.set_item("quadrature", &o.quadrature)?;
    d.set_item("length", "edges")?;
    d.set_item("grow_factor", o.grow_factor)?;
    d.set_item("grow_rate", o.grow_rate)?;
    d.set_item("max_subdivisions", o.max_subdivisions)?;
    d.set_item("barycenter", py.None())?;
    d.set_item("pinned", py.None())?;
    d.set_item("tangent_vertices", py.None())?;
    d.set_item("tangent_edges", py.None())?;
    d.set_item("surface", py.None())?;
    d.set_item("planar_axis", py.None())?;
    d.set_item("obstacles", py.None())?;
    d.set_item("obstacle_weight", o.obstacle_weight)?;
    d.set_item("obstacle_exponent", o.obstacle_exponent)?;
    d.set_item("obstacle_theta", o.obstacle_theta)?;
    d.set_item("total_length_weight", o.total_length_weight)?;
    d.set_item("length_difference_weight", o.length_difference_weight)?;
    d.set_item("field", py.None())?;
    d.set_item("field_weight", o.field_weight)?;
    d.set_item("line_search", &o.line_search)?;
    d.set_item("line_search_start", &o.line_search_start)?;
    d.set_item("line_search_start_factor", o.line_search_start_factor)?;
    d.set_item("armijo", o.armijo)?;
    d.set_item("shrink", o.shrink)?;
    d.set_item("max_halvings", o.max_halvings)?;
    d.set_item("gap", o.gap)?;
    d.set_item("backend", &o.backend)?;
    d.set_item("auto_dense_max_vertices", o.auto_dense_max_vertices)?;
    d.set_item("dense_fallback_max_vertices", o.dense_fallback_max_vertices)?;
    d.set_item("mg_coarsest_vertices", o.mg_coarsest_vertices)?;
    d.set_item("mg_max_levels", o.mg_max_levels)?;
    d.set_item("mg_smoothing_steps", o.mg_smoothing_steps)?;
    d.set_item("mg_initial_smoothing_steps", o.mg_initial_smoothing_steps)?;
    d.set_item("mg_tolerance", o.mg_tolerance)?;
    d.set_item("mg_max_cycles", o.mg_max_cycles)?;
    d.set_item("mg_flexible_cg", o.mg_flexible_cg)?;
    d.set_item("mg_near_field_factor", o.mg_near_field_factor)?;
    d.set_item("mg_smoother_factor", o.mg_smoother_factor)?;
    d.set_item("mg_max_fixed_fraction", o.mg_max_fixed_fraction)?;
    d.set_item("mg_energy_check", o.mg_energy_check)?;
    d.set_item("mg_failure_residual", o.mg_failure_residual)?;
    d.set_item("bh_theta_x", o.bh_theta_x)?;
    d.set_item("bh_theta_t", o.bh_theta_t)?;
    d.set_item("bh_leaf_size", o.bh_leaf_size)?;
    d.set_item("threads", o.threads)?;
    d.set_item("dense_threads", o.dense_threads)?;
    d.set_item("min_edge_length_ratio", o.min_edge_length_ratio)?;
    d.set_item("projection_tolerance", o.projection_tolerance)?;
    d.set_item("projection_absolute", o.projection_absolute)?;
    d.set_item("projection_iterations", o.projection_iterations)?;
    d.set_item("stop", &o.stop)?;
    d.set_item("stop_threshold", o.stop_threshold)?;
    d.set_item("stall_window", o.stall_window)?;
    d.set_item("gradient_norm", &o.gradient_norm)?;
    d.set_item("max_steps", o.max_steps)?;
    d.set_item("time_limit", o.time_limit)?;
    Ok(d)
}
