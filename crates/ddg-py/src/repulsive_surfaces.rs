//! Repulsive Surfaces: the binding of `ddg-repulsive-surfaces`.
//!
//! Yu, Brakensiek, Schumacher, Crane, *Repulsive Surfaces*, ACM Transactions on Graphics 40(6),
//! 2021: the tangent-point energy of a triangle mesh and its fractional Sobolev descent with
//! constraints, penalties, obstacles, continuous collision detection and remeshing. This
//! module holds everything that is not PyO3 (the glue is `python/repulsive_surfaces.rs`):
//!
//! * [`Val`]: one value tree for option values going in and reports going out (the glue turns
//!   Python objects into it and back), so that the option parser is plain Rust and tested by
//!   `cargo test`. Buffers become [`Val::Array`]; penalties and CCD obstacles that cannot be
//!   written as options (from a scene file) travel as opaque [`Val::Penalty`] / [`Val::Obstacle`];
//! * [`flow_options`]: options as a dict-like value -> the library's [`FlowOptions`] (the
//!   targets that mean "the initial value" are resolved against the start mesh here);
//!   [`options_val`] is its inverse, and [`defaults`] is `options_val` of the library's defaults;
//! * [`Engine`]: a session that can be driven from several threads. `step` takes the session's
//!   lock; `snapshot`, `cancel` never wait for a running step (the last completed state is stored
//!   separately and the cancel flag is an atomic);
//! * [`flow_once`], [`surface_energy`], [`scene`].
//!
//! # Remeshing
//!
//! With remeshing on (the default) the vertex and face counts change in every step: the
//! snapshot carries the current faces, the current pinned vertex ids and, per face, the face of
//! the *initial* mesh it descends from (`face_origin`, for a material index).
//!
//! # Errors
//!
//! Everything the caller can get wrong is [`SError::Input`] carrying the library's message
//! (the glue raises `ValueError`); a failed computation is [`SError::Solve`].

use std::collections::BTreeSet;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};
use std::time::{Duration, Instant};

use ddg_repulsive_curves::obstacle::Obstacle;
use ddg_repulsive_surfaces as rs;
use rs::implicit::{Cylinder, FramedTorus, ImplicitSurface, Plane, Sphere};
use rs::penalties::simple::{area_and_gradient, volume_and_gradient};
use rs::remesh::{EnabledOps, RemeshOptions, RemeshOrder};
use rs::gpu::{CudaOptions, Determinism, Device};
use rs::solve::iterative::HierarchicalOptions;
use rs::{
    ArmijoRule, Backend, BarycenterMode, BhOptions, CcdOptions, Constraint, Control, Energy,
    BcqnVariant, EnergyEval, Exponents, FlowOptions, GmresSide, L2Kind, LbfgsVariant, LineSearch, Method,
    Metric, NewtonSteps, ObstacleKind, PinSchedule, Penalty, Quadrature, RsError, Schedule,
    ScheduleTiming, Session, StartRule, StaticSet, Status, StepReport, StepStatus, StopRule,
    Surface,
};

/// Leaf size of the BVHs of mesh obstacles.
const OBSTACLE_LEAF_SIZE: usize = 8;

// ---------------------------------------------------------------------------------------------
// Errors and values
// ---------------------------------------------------------------------------------------------

/// What can go wrong.
#[derive(Debug, Clone, PartialEq)]
pub enum SError {
    /// The caller's input or options are wrong; the message is the library's or names the option.
    Input(String),
    /// The input was accepted but a computation failed.
    Solve(String),
}

impl std::fmt::Display for SError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SError::Input(m) | SError::Solve(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for SError {}

impl From<RsError> for SError {
    fn from(e: RsError) -> Self {
        match e {
            RsError::NonFiniteEnergy { .. }
            | RsError::SingularSaddle { .. }
            | RsError::SolverDiverged { .. }
            | RsError::Internal { .. }
            | RsError::Cancelled => SError::Solve(e.to_string()),
            e => SError::Input(e.to_string()),
        }
    }
}

impl From<ddg_repulsive_curves::RcError> for SError {
    fn from(e: ddg_repulsive_curves::RcError) -> Self {
        SError::Input(e.to_string())
    }
}

fn bad(name: &str, message: impl std::fmt::Display) -> SError {
    SError::Input(format!("{name}: {message}"))
}

/// A value of an option or a report. The glue turns it into `None`, `bool`, `int`, `float`,
/// `str`, `list`, `dict`, `array.array('d')` or an opaque handle.
#[derive(Clone)]
pub enum Val {
    /// `None`.
    None,
    /// `bool`.
    Bool(bool),
    /// `int`.
    Int(i64),
    /// `float`.
    Float(f64),
    /// `str`.
    Str(String),
    /// `list`.
    List(Vec<Val>),
    /// A flat numeric buffer (`array`, numpy); indices are stored as reals.
    Array(Vec<f64>),
    /// `dict` (insertion order).
    Dict(Vec<(String, Val)>),
    /// A penalty built elsewhere (a scene file); accepted in the `penalties` list.
    Penalty(Penalty),
    /// A triangle obstacle built elsewhere; accepted in the `ccd_obstacles` list.
    Obstacle(Arc<Obstacle>),
}

impl std::fmt::Debug for Val {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Val::None => f.write_str("None"),
            Val::Bool(b) => write!(f, "{b}"),
            Val::Int(i) => write!(f, "{i}"),
            Val::Float(x) => write!(f, "{x:?}"),
            Val::Str(s) => write!(f, "{s:?}"),
            Val::List(l) => f.debug_list().entries(l).finish(),
            Val::Array(a) => write!(f, "array[{}]", a.len()),
            Val::Dict(d) => f.debug_map().entries(d.iter().map(|(k, v)| (k, v))).finish(),
            Val::Penalty(p) => write!(f, "{p:?}"),
            Val::Obstacle(_) => f.write_str("obstacle"),
        }
    }
}

impl Val {
    /// The entry `key` of a dict.
    pub fn get(&self, key: &str) -> Option<&Val> {
        match self {
            Val::Dict(d) => d.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    fn s(x: &str) -> Val {
        Val::Str(x.into())
    }

    fn uint(x: usize) -> Val {
        Val::Int(x as i64)
    }

    fn opt_f(x: Option<f64>) -> Val {
        x.map_or(Val::None, Val::Float)
    }

    fn floats(x: &[f64]) -> Val {
        Val::List(x.iter().map(|&v| Val::Float(v)).collect())
    }

    fn dict(items: Vec<(&str, Val)>) -> Val {
        Val::Dict(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
    }

    fn type_name(&self) -> &'static str {
        match self {
            Val::None => "None",
            Val::Bool(_) => "bool",
            Val::Int(_) => "int",
            Val::Float(_) => "float",
            Val::Str(_) => "str",
            Val::List(_) => "list",
            Val::Array(_) => "array",
            Val::Dict(_) => "dict",
            Val::Penalty(_) => "penalty",
            Val::Obstacle(_) => "obstacle",
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Reading values
// ---------------------------------------------------------------------------------------------

fn as_f64(name: &str, v: &Val) -> Result<f64> {
    match v {
        Val::Float(x) => Ok(*x),
        Val::Int(i) => Ok(*i as f64),
        other => Err(bad(name, format!("expected a number, got {}", other.type_name()))),
    }
}

fn as_usize(name: &str, v: &Val) -> Result<usize> {
    match v {
        Val::Int(i) => usize::try_from(*i).map_err(|_| bad(name, format!("{i} must be >= 0"))),
        Val::Float(x) if x.fract() == 0.0 && *x >= 0.0 && *x < 1e15 => Ok(*x as usize),
        other => Err(bad(
            name,
            format!("expected a non-negative integer, got {}", other.type_name()),
        )),
    }
}

fn as_bool(name: &str, v: &Val) -> Result<bool> {
    match v {
        Val::Bool(b) => Ok(*b),
        other => Err(bad(name, format!("expected a bool, got {}", other.type_name()))),
    }
}

fn as_str<'a>(name: &str, v: &'a Val) -> Result<&'a str> {
    match v {
        Val::Str(s) => Ok(s),
        other => Err(bad(name, format!("expected a string, got {}", other.type_name()))),
    }
}

fn flatten_into(name: &str, v: &Val, out: &mut Vec<f64>) -> Result<()> {
    match v {
        Val::Float(x) => out.push(*x),
        Val::Int(i) => out.push(*i as f64),
        Val::Array(a) => out.extend_from_slice(a),
        Val::List(l) => {
            for i in l {
                flatten_into(name, i, out)?;
            }
        }
        other => {
            return Err(bad(
                name,
                format!("expected numbers, got {}", other.type_name()),
            ));
        }
    }
    Ok(())
}

fn reals(name: &str, v: &Val) -> Result<Vec<f64>> {
    let mut out = Vec::new();
    flatten_into(name, v, &mut out)?;
    if let Some(i) = out.iter().position(|x| !x.is_finite()) {
        return Err(bad(name, format!("entry {i} is not finite")));
    }
    Ok(out)
}

fn points(name: &str, v: &Val) -> Result<Vec<[f64; 3]>> {
    let flat = reals(name, v)?;
    if !flat.len().is_multiple_of(3) {
        return Err(bad(
            name,
            format!("{} numbers are not a whole number of points of 3", flat.len()),
        ));
    }
    Ok(flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect())
}

fn vec3(name: &str, v: &Val) -> Result<[f64; 3]> {
    let flat = reals(name, v)?;
    <[f64; 3]>::try_from(flat)
        .map_err(|x| bad(name, format!("expected 3 numbers, got {}", x.len())))
}

type Result<T> = std::result::Result<T, SError>;

fn indices(name: &str, v: &Val) -> Result<Vec<u32>> {
    let flat = reals(name, v)?;
    flat.iter()
        .enumerate()
        .map(|(i, &x)| {
            if x >= 0.0 && x.fract() == 0.0 && x < f64::from(u32::MAX) {
                Ok(x as u32)
            } else {
                Err(bad(name, format!("entry {i} = {x} is not a valid index")))
            }
        })
        .collect()
}

fn triangles(name: &str, v: &Val) -> Result<Vec<[u32; 3]>> {
    let flat = indices(name, v)?;
    if !flat.len().is_multiple_of(3) {
        return Err(bad(
            name,
            format!("{} indices are not a whole number of triangles", flat.len()),
        ));
    }
    Ok(flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect())
}

/// A dict whose keys must all be read; `finish` reports the first unknown key.
struct Dict<'a> {
    what: String,
    items: &'a [(String, Val)],
    used: BTreeSet<&'a str>,
}

impl<'a> Dict<'a> {
    fn new(what: &str, v: &'a Val) -> Result<Dict<'a>> {
        match v {
            Val::Dict(items) => Ok(Dict {
                what: what.to_owned(),
                items,
                used: BTreeSet::new(),
            }),
            other => Err(bad(
                what,
                format!("expected a dict, got {}", other.type_name()),
            )),
        }
    }

    /// The value of `key`, `None` if absent or `None`.
    fn get(&mut self, key: &str) -> Option<&'a Val> {
        let (k, v) = self.items.iter().find(|(k, _)| k == key)?;
        self.used.insert(k.as_str());
        (!matches!(v, Val::None)).then_some(v)
    }

    fn name(&self, key: &str) -> String {
        format!("{}.{key}", self.what)
    }

    fn require(&mut self, key: &str) -> Result<&'a Val> {
        let n = self.what.clone();
        self.get(key).ok_or_else(|| bad(&n, format!("missing '{key}'")))
    }

    fn f64(&mut self, key: &str, default: f64) -> Result<f64> {
        match self.get(key) {
            Some(v) => as_f64(&self.name(key), v),
            None => Ok(default),
        }
    }

    fn opt_f64(&mut self, key: &str) -> Result<Option<f64>> {
        self.get(key).map(|v| as_f64(&self.name(key), v)).transpose()
    }

    fn req_f64(&mut self, key: &str) -> Result<f64> {
        let v = self.require(key)?;
        as_f64(&self.name(key), v)
    }

    fn vec3(&mut self, key: &str, default: Option<[f64; 3]>) -> Result<[f64; 3]> {
        match self.get(key) {
            Some(v) => vec3(&self.name(key), v),
            None => default.ok_or_else(|| bad(&self.what, format!("missing '{key}'"))),
        }
    }

    fn string(&mut self, key: &str, default: &'static str) -> Result<&'a str> {
        match self.get(key) {
            Some(v) => as_str(&self.name(key), v),
            None => Ok(default),
        }
    }

    fn bool(&mut self, key: &str, default: bool) -> Result<bool> {
        match self.get(key) {
            Some(v) => as_bool(&self.name(key), v),
            None => Ok(default),
        }
    }

    fn finish(&self) -> Result<()> {
        for (k, _) in self.items {
            if !self.used.contains(k.as_str()) {
                return Err(bad(&self.what, format!("unknown key '{k}'")));
            }
        }
        Ok(())
    }
}

fn one_of(name: &str, value: &str, allowed: &[&str]) -> Result<()> {
    if allowed.contains(&value) {
        Ok(())
    } else {
        Err(bad(
            name,
            format!("'{value}' is not valid; use one of {}", allowed.join(", ")),
        ))
    }
}

// ---------------------------------------------------------------------------------------------
// Implicit surfaces and penalties
// ---------------------------------------------------------------------------------------------

const IMPLICIT_TYPES: &[&str] = &["sphere", "cylinder", "torus", "plane"];

fn implicit_surface(what: &str, v: &Val) -> Result<Arc<dyn ImplicitSurface>> {
    let mut d = Dict::new(what, v)?;
    let kind = d.string("type", "")?;
    if kind.is_empty() {
        return Err(bad(what, "missing 'type'"));
    }
    one_of(&format!("{what}.type"), kind, IMPLICIT_TYPES)?;
    let s: Arc<dyn ImplicitSurface> = match kind {
        "sphere" => Arc::new(Sphere::new(d.vec3("center", Some([0.0; 3]))?, d.f64("radius", 1.0)?)?),
        "cylinder" => Arc::new(Cylinder::new(
            d.vec3("point", Some([0.0; 3]))?,
            d.vec3("axis", Some([0.0, 1.0, 0.0]))?,
            d.req_f64("radius")?,
        )?),
        "torus" => Arc::new(FramedTorus::new(
            d.vec3("center", Some([0.0; 3]))?,
            d.vec3("axis", Some([0.0, 1.0, 0.0]))?,
            d.req_f64("major")?,
            d.req_f64("minor")?,
        )?),
        _ => Arc::new(Plane::new(
            d.vec3("point", Some([0.0; 3]))?,
            d.vec3("normal", None)?,
        )?),
    };
    d.finish()?;
    Ok(s)
}

const PENALTY_TYPES: &[&str] = &[
    "area_deviation",
    "volume_deviation",
    "total_area",
    "total_volume",
    "squared_error",
    "boundary_length",
    "boundary_curvature",
    "willmore",
    "implicit_obstacle",
    "implicit_attractor",
    "mesh_obstacle",
    "point_cloud_obstacle",
];

fn quadrature(name: &str, s: &str) -> Result<Quadrature> {
    one_of(name, s, &["dual_area", "per_vertex"])?;
    Ok(if s == "dual_area" {
        Quadrature::DualArea
    } else {
        Quadrature::PerVertex
    })
}

/// One penalty dict; `Some(obstacle)` for a mesh obstacle that also guards the line search.
fn penalty(
    i: usize,
    v: &Val,
    surface: &Surface,
    positions: &[[f64; 3]],
) -> Result<(Penalty, Option<Arc<Obstacle>>)> {
    if let Val::Penalty(p) = v {
        return Ok((p.clone(), None));
    }
    let what = format!("penalties[{i}]");
    let mut d = Dict::new(&what, v)?;
    let kind = d.string("type", "")?;
    if kind.is_empty() {
        return Err(bad(&what, "missing 'type'"));
    }
    one_of(&format!("{what}.type"), kind, PENALTY_TYPES)?;
    let weight = d.f64("weight", 1.0)?;
    let mut ccd = None;
    let p = match kind {
        "area_deviation" => Penalty::AreaDeviation {
            weight,
            target: match d.opt_f64("target")? {
                Some(t) => t,
                None => area_and_gradient(surface, positions).0,
            },
        },
        "volume_deviation" => Penalty::VolumeDeviation {
            weight,
            target: match d.opt_f64("target")? {
                Some(t) => t,
                None => volume_and_gradient(surface, positions).0,
            },
        },
        "total_area" => Penalty::TotalArea {
            weight,
            power: d.f64("power", 1.0)?,
        },
        "total_volume" => Penalty::TotalVolume {
            weight,
            power: d.f64("power", 1.0)?,
        },
        "squared_error" => Penalty::SquaredError {
            weight,
            targets: match d.get("targets") {
                Some(t) => points(&d.name("targets"), t)?,
                None => positions.to_vec(),
            },
        },
        "boundary_length" => Penalty::BoundaryLength {
            weight,
            target: d.f64("target", 0.0)?,
        },
        "boundary_curvature" => Penalty::BoundaryCurvature {
            weight,
            dual_length: d.bool("dual_length", true)?,
        },
        "willmore" => Penalty::Willmore { weight },
        "implicit_obstacle" | "implicit_attractor" => {
            let surf = implicit_surface(&d.name("surface"), d.require("surface")?)?;
            let power = d.f64("power", 2.0)?;
            let q = quadrature(&d.name("quadrature"), d.string("quadrature", "dual_area")?)?;
            if kind == "implicit_obstacle" {
                Penalty::ImplicitObstacle {
                    weight,
                    surface: surf,
                    power,
                    quadrature: q,
                }
            } else {
                Penalty::ImplicitAttractor {
                    weight,
                    surface: surf,
                    power,
                    quadrature: q,
                    vertices: d
                        .get("vertices")
                        .map(|v| indices(&d.name("vertices"), v))
                        .transpose()?,
                }
            }
        }
        "mesh_obstacle" => {
            let vertices = points(&d.name("vertices"), d.require("vertices")?)?;
            let tris = triangles(&d.name("triangles"), d.require("triangles")?)?;
            let obstacle = Arc::new(Obstacle::new(vertices, tris)?);
            let kind = match d.string("kind", "tangent_point")? {
                "tangent_point" => ObstacleKind::TangentPoint,
                "coulomb" => ObstacleKind::Coulomb {
                    power: d.f64("power", 2.0)?,
                },
                other => {
                    return Err(bad(
                        &d.name("kind"),
                        format!("'{other}' is not valid; use tangent_point or coulomb"),
                    ));
                }
            };
            if d.bool("ccd", true)? {
                ccd = Some(obstacle.clone());
            }
            Penalty::MeshObstacle {
                weight,
                obstacle: Arc::new(StaticSet::from_obstacle(&obstacle, OBSTACLE_LEAF_SIZE)?),
                kind,
            }
        }
        _ => Penalty::PointCloudObstacle {
            weight,
            points: Arc::new(StaticSet::from_points(
                &points(&d.name("points"), d.require("points")?)?,
                OBSTACLE_LEAF_SIZE,
            )?),
        },
    };
    d.finish()?;
    Ok((p, ccd))
}

// ---------------------------------------------------------------------------------------------
// Constraints
// ---------------------------------------------------------------------------------------------

#[derive(Default)]
struct Constraints {
    list: Vec<Constraint>,
    schedules: Vec<Schedule>,
    pin_schedules: Vec<PinSchedule>,
}

fn pins_constraint(
    what: &str,
    vertices: Vec<u32>,
    targets: Option<&Val>,
    movement: Option<&Val>,
    steps: Option<usize>,
    out: &mut Constraints,
) -> Result<()> {
    let n = vertices.len();
    let targets = targets
        .map(|t| {
            let p = points(&format!("{what}.targets"), t)?;
            if p.len() != n {
                return Err(bad(
                    &format!("{what}.targets"),
                    format!("{} targets for {n} pins", p.len()),
                ));
            }
            Ok(p)
        })
        .transpose()?;
    if let Some(m) = movement {
        let name = format!("{what}.move");
        let flat = reals(&name, m)?;
        let displacement = if flat.len() == 3 {
            vec![[flat[0], flat[1], flat[2]]; n]
        } else {
            let p = points(&name, m)?;
            if p.len() != n {
                return Err(bad(&name, format!("{} displacements for {n} pins", p.len())));
            }
            p
        };
        let steps = steps.ok_or_else(|| bad(what, "'move' needs 'steps'"))?;
        out.pin_schedules.push(PinSchedule {
            constraint: out.list.len(),
            displacement,
            steps,
        });
    } else if steps.is_some() {
        return Err(bad(what, "'steps' without 'move'"));
    }
    out.list.push(Constraint::Pins { vertices, targets });
    Ok(())
}

fn constraint(i: usize, v: &Val, out: &mut Constraints) -> Result<()> {
    let what = format!("constraints[{i}]");
    let mut d = Dict::new(&what, v)?;
    let kind = d.string("type", "")?;
    one_of(
        &format!("{what}.type"),
        kind,
        &["barycenter", "pins", "area", "volume"],
    )?;
    match kind {
        "barycenter" => {
            let component = d.get("component").map(|v| as_usize(&d.name("component"), v)).transpose()?;
            let target = d.get("target").map(|v| vec3(&d.name("target"), v)).transpose()?;
            out.list.push(Constraint::Barycenter {
                component: component.map(|c| c as u32),
                target,
            });
        }
        "pins" => {
            let vertices = indices(&d.name("vertices"), d.require("vertices")?)?;
            let steps = d.get("steps").map(|v| as_usize(&d.name("steps"), v)).transpose()?;
            let (t, m) = (d.get("targets"), d.get("move"));
            pins_constraint(&what, vertices, t, m, steps, out)?;
        }
        _ => {
            let component = d.get("component").map(|v| as_usize(&d.name("component"), v)).transpose()?;
            let target = d.opt_f64("target")?;
            let (factor, offset) = (d.opt_f64("factor")?, d.opt_f64("offset")?);
            let steps = d.get("steps").map(|v| as_usize(&d.name("steps"), v)).transpose()?;
            let component = component.map(|c| c as u32);
            if factor.is_some() || offset.is_some() {
                out.schedules.push(Schedule {
                    constraint: out.list.len(),
                    factor: factor.unwrap_or(1.0),
                    offset: offset.unwrap_or(0.0),
                    steps: steps.ok_or_else(|| bad(&what, "a growth schedule needs 'steps'"))?,
                });
            } else if steps.is_some() {
                return Err(bad(&what, "'steps' without 'factor' or 'offset'"));
            }
            out.list.push(if kind == "area" {
                Constraint::Area { component, target }
            } else {
                Constraint::Volume { component, target }
            });
        }
    }
    d.finish()
}

// ---------------------------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------------------------

/// The names of all options.
pub const OPTION_NAMES: &[&str] = &[
    "alpha",
    "beta",
    "energy",
    "coulomb_alpha",
    "tp_weight",
    "metric",
    "low_order",
    "l2_kind",
    "backend",
    "dense_max_faces",
    "energy_eval",
    "theta",
    "leaf_size",
    "chi",
    "cluster_leaf_size",
    "gmres_restart",
    "gmres_tolerance",
    "gmres_max_iterations",
    "gmres_side",
    "device",
    "cuda_ordinal",
    "gpu_determinism",
    "method",
    "aqp_kappa",
    "lbfgs_memory",
    "lbfgs_variant",
    "constraints",
    "barycenter",
    "barycenter_mode",
    "barycenter_shift",
    "pins",
    "pin_targets",
    "pin_move",
    "pin_move_steps",
    "schedule_timing",
    "penalties",
    "ccd_obstacles",
    "armijo",
    "armijo_c",
    "start",
    "start_lengths",
    "start_factor",
    "start_tau",
    "tau_min_rel",
    "max_halvings",
    "ccd",
    "ccd_gap",
    "ccd_max_iterations",
    "ccd_close_factor",
    "newton_steps",
    "newton_max",
    "remesh",
    "remesh_split",
    "remesh_collapse",
    "remesh_rounds",
    "remesh_rho",
    "remesh_flip_cap",
    "remesh_order",
    "remesh_ops",
    "stop",
    "stop_relative",
    "stop_window",
    "stop_eps",
    "max_iterations",
    "time_limit",
    "stuck_limit",
    "threads",
    "reference_compat",
];

/// Names of the option values that are enumerations, for the tests and the stub.
pub const REMESH_OPS: &[&str] = &["split", "collapse", "flip", "smooth"];

const DEFAULT_AQP_KAPPA: f64 = 100.0;
const DEFAULT_LBFGS_MEMORY: usize = 10;

/// Reads the options of a flow. Unknown names are errors naming the option; the start mesh
/// resolves the targets that default to "the initial value" of penalties.
///
/// # Errors
/// [`SError::Input`] with the option's name or the library's message.
pub fn flow_options(
    options: &[(String, Val)],
    surface: &Surface,
    positions: &[[f64; 3]],
) -> Result<FlowOptions> {
    for (k, _) in options {
        if !OPTION_NAMES.contains(&k.as_str()) {
            return Err(SError::Input(format!(
                "unknown option '{k}'; repulsive_surface_defaults() lists the options"
            )));
        }
    }
    let mut d = Dict {
        what: "options".to_owned(),
        items: options,
        used: BTreeSet::new(),
    };
    let mut o = FlowOptions::default();
    let n = |k: &str| k.to_owned();

    o.exponents = Exponents::new(d.f64("alpha", 6.0)?, d.f64("beta", 12.0)?)?;
    let coulomb = d.opt_f64("coulomb_alpha")?;
    let energy = d.string("energy", "tangent_point")?;
    one_of("energy", energy, &["tangent_point", "coulomb", "willmore"])?;
    o.energy = match energy {
        "tangent_point" => Energy::TangentPoint,
        "willmore" => Energy::Willmore,
        _ => Energy::Coulomb {
            alpha: coulomb.ok_or_else(|| bad("coulomb_alpha", "required with energy='coulomb'"))?,
        },
    };
    o.tp_weight = d.f64("tp_weight", 1.0)?;

    let low_order = d.bool("low_order", true)?;
    let l2 = d.string("l2_kind", "mass")?;
    one_of("l2_kind", l2, &["mass", "identity"])?;
    let metric = d.string("metric", "sobolev")?;
    one_of(&n("metric"), metric, &["sobolev", "l2", "h1", "h2"])?;
    o.metric = match metric {
        "sobolev" => Metric::Sobolev { low_order },
        "l2" => Metric::L2(if l2 == "mass" { L2Kind::Mass } else { L2Kind::Identity }),
        "h1" => Metric::H1,
        _ => Metric::H2,
    };

    // Back end and energy evaluation.
    let mut h = HierarchicalOptions::default();
    h.cluster.chi = d.f64("chi", h.cluster.chi)?;
    h.cluster.leaf_size = d
        .get("cluster_leaf_size")
        .map(|v| as_usize("cluster_leaf_size", v))
        .transpose()?
        .unwrap_or(h.cluster.leaf_size);
    if let Some(v) = d.get("gmres_restart") {
        h.gmres.restart = as_usize("gmres_restart", v)?;
    }
    h.gmres.tolerance = d.f64("gmres_tolerance", h.gmres.tolerance)?;
    if let Some(v) = d.get("gmres_max_iterations") {
        h.gmres.max_iterations = as_usize("gmres_max_iterations", v)?;
    }
    let side = d.string("gmres_side", "right")?;
    one_of("gmres_side", side, &["right", "left"])?;
    h.gmres.side = if side == "right" { GmresSide::Right } else { GmresSide::Left };
    let device = d.string("device", "cpu")?;
    one_of("device", device, &["cpu", "cuda", "auto"])?;
    let ordinal = d
        .get("cuda_ordinal")
        .map(|v| as_usize("cuda_ordinal", v))
        .transpose()?
        .unwrap_or(0);
    let det = d.string("gpu_determinism", "run_to_run")?;
    one_of("gpu_determinism", det, &["run_to_run", "match_cpu"])?;
    let cuda = CudaOptions {
        ordinal,
        determinism: if det == "match_cpu" { Determinism::MatchCpu } else { Determinism::RunToRun },
    };
    h.device = match device {
        "cuda" => Device::Cuda(cuda),
        "auto" => Device::Auto(cuda),
        _ => Device::Cpu,
    };
    let dense_max = d
        .get("dense_max_faces")
        .map(|v| as_usize("dense_max_faces", v))
        .transpose()?
        .unwrap_or(rs::solve::iterative::DEFAULT_DENSE_MAX_FACES);
    let backend = d.string("backend", "auto")?;
    one_of("backend", backend, &["auto", "dense", "hierarchical"])?;
    o.backend = match backend {
        "auto" => Backend::Auto {
            dense_max_faces: dense_max,
            hierarchical: h,
        },
        "dense" => Backend::Dense,
        _ => Backend::Hierarchical(h),
    };
    let mut bh = BhOptions::default();
    bh.theta = d.f64("theta", bh.theta)?;
    if let Some(v) = d.get("leaf_size") {
        bh.leaf_size = as_usize("leaf_size", v)?;
    }
    let eval = d.string("energy_eval", "barnes_hut")?;
    one_of("energy_eval", eval, &["barnes_hut", "exact"])?;
    o.energy_eval = if eval == "exact" { EnergyEval::Exact } else { EnergyEval::BarnesHut(bh) };

    // Method.
    let kappa = d.f64("aqp_kappa", DEFAULT_AQP_KAPPA)?;
    let memory = d
        .get("lbfgs_memory")
        .map(|v| as_usize("lbfgs_memory", v))
        .transpose()?
        .unwrap_or(DEFAULT_LBFGS_MEMORY);
    let variant = d.string("lbfgs_variant", "standard")?;
    one_of("lbfgs_variant", variant, &["standard", "reference"])?;
    let variant = if variant == "standard" { LbfgsVariant::Standard } else { LbfgsVariant::Reference };
    let method = d.string("method", "projected_gradient")?;
    one_of(
        "method",
        method,
        &["projected_gradient", "aqp", "h1_lbfgs", "bcqn"],
    )?;
    o.method = match method {
        "projected_gradient" => Method::ProjectedGradient,
        "aqp" => Method::Aqp { kappa },
        "h1_lbfgs" => Method::H1Lbfgs { memory, variant },
        _ => Method::Bcqn {
            memory,
            variant: if variant == LbfgsVariant::Standard {
                BcqnVariant::Standard
            } else {
                BcqnVariant::Reference
            },
        },
    };

    // Constraints.
    let mut cs = Constraints::default();
    let pins = d.get("pins").map(|v| indices("pins", v)).transpose()?;
    let explicit = d.get("constraints");
    let bary_opt = d.get("barycenter");
    match explicit {
        Some(Val::List(items)) => {
            for (i, c) in items.iter().enumerate() {
                constraint(i, c, &mut cs)?;
            }
        }
        Some(other) => {
            return Err(bad(
                "constraints",
                format!("expected a list of dicts, got {}", other.type_name()),
            ));
        }
        None => {
            let bary = match bary_opt {
                Some(v) => as_bool("barycenter", v)?,
                None => pins.is_none(),
            };
            if bary {
                cs.list.push(Constraint::Barycenter {
                    component: None,
                    target: None,
                });
            }
        }
    }
    let pin_targets = d.get("pin_targets");
    let pin_move = d.get("pin_move");
    let pin_steps = d.get("pin_move_steps").map(|v| as_usize("pin_move_steps", v)).transpose()?;
    match pins {
        Some(p) => pins_constraint("pins", p, pin_targets, pin_move, pin_steps, &mut cs)?,
        None if pin_targets.is_some() || pin_move.is_some() || pin_steps.is_some() => {
            return Err(bad("pin_targets", "pin_targets, pin_move and pin_move_steps need 'pins'"));
        }
        None => {}
    }
    o.constraints = cs.list;
    o.schedules = cs.schedules;
    o.pin_schedules = cs.pin_schedules;
    o.barycenter_mode = match d.string("barycenter_mode", "auto")? {
        "auto" => None,
        "global" => Some(BarycenterMode::Global),
        "per_component" => Some(BarycenterMode::PerComponent),
        other => {
            return Err(bad(
                "barycenter_mode",
                format!("'{other}' is not valid; use auto, global or per_component"),
            ));
        }
    };
    o.barycenter_shift = d.bool("barycenter_shift", false)?;
    o.schedule_timing = match d.string("schedule_timing", "before_projection")? {
        "before_projection" => ScheduleTiming::BeforeProjection,
        "after_projection" => ScheduleTiming::AfterProjection,
        other => {
            return Err(bad(
                "schedule_timing",
                format!("'{other}' is not valid; use before_projection or after_projection"),
            ));
        }
    };

    // Penalties and obstacles.
    match d.get("penalties") {
        Some(Val::List(items)) => {
            for (i, p) in items.iter().enumerate() {
                let (p, ccd) = penalty(i, p, surface, positions)?;
                o.penalties.push(p);
                o.ccd_obstacles.extend(ccd);
            }
        }
        Some(other) => {
            return Err(bad(
                "penalties",
                format!("expected a list of dicts, got {}", other.type_name()),
            ));
        }
        None => {}
    }
    match d.get("ccd_obstacles") {
        Some(Val::List(items)) => {
            for v in items {
                match v {
                    Val::Obstacle(ob) => o.ccd_obstacles.push(ob.clone()),
                    other => {
                        return Err(bad(
                            "ccd_obstacles",
                            format!("expected obstacle handles, got {}", other.type_name()),
                        ));
                    }
                }
            }
        }
        Some(other) => {
            return Err(bad(
                "ccd_obstacles",
                format!("expected a list, got {}", other.type_name()),
            ));
        }
        None => {}
    }

    // Line search.
    let armijo_c = d.f64("armijo_c", 1e-2)?;
    let armijo = d.string("armijo", "standard")?;
    one_of("armijo", armijo, &["standard", "reference"])?;
    let lengths = d.f64("start_lengths", 1.0)?;
    let factor = d.f64("start_factor", 2.0)?;
    let tau = d.opt_f64("start_tau")?;
    let start = d.string("start", "grow_previous")?;
    one_of(
        "start",
        start,
        &["grow_previous", "max_displacement", "reference", "absolute"],
    )?;
    let mut ls = LineSearch::default();
    ls.armijo = if armijo == "standard" {
        ArmijoRule::Standard { c: armijo_c }
    } else {
        ArmijoRule::Reference { c: armijo_c }
    };
    ls.start = match start {
        "grow_previous" => StartRule::GrowPrevious {
            factor,
            first: Box::new(StartRule::MaxDisplacement { lengths }),
        },
        "max_displacement" => StartRule::MaxDisplacement { lengths },
        "reference" => StartRule::Reference,
        _ => StartRule::Absolute {
            tau: tau.ok_or_else(|| bad("start_tau", "required with start='absolute'"))?,
        },
    };
    ls.tau_min_rel = d.f64("tau_min_rel", ls.tau_min_rel)?;
    ls.max_halvings = d
        .get("max_halvings")
        .map(|v| as_usize("max_halvings", v).map(|x| x as u32))
        .transpose()?;
    let mut ccd = CcdOptions::default();
    ccd.gap = d.f64("ccd_gap", ccd.gap)?;
    if let Some(v) = d.get("ccd_max_iterations") {
        ccd.max_iterations = as_usize("ccd_max_iterations", v)? as u32;
    }
    ccd.close_factor = d.f64("ccd_close_factor", ccd.close_factor)?;
    ls.ccd = d.bool("ccd", true)?.then_some(ccd);
    o.line_search = ls;

    // Projection.
    let newton_max = d
        .get("newton_max")
        .map(|v| as_usize("newton_max", v))
        .transpose()?
        .unwrap_or(5) as u32;
    o.newton_steps = match d.get("newton_steps") {
        Some(v) => NewtonSteps::Fixed(as_usize("newton_steps", v)? as u32),
        None => NewtonSteps::UntilTolerance { max: newton_max },
    };

    // Remeshing.
    let mut r = RemeshOptions::default();
    r.split = d.f64("remesh_split", r.split)?;
    r.collapse = d.f64("remesh_collapse", r.collapse)?;
    if let Some(v) = d.get("remesh_rounds") {
        r.rounds = as_usize("remesh_rounds", v)?;
    }
    r.rho = d.f64("remesh_rho", r.rho)?;
    if let Some(v) = d.get("remesh_flip_cap") {
        r.flip_cap_per_vertex = as_usize("remesh_flip_cap", v)?;
    }
    r.order = match d.string("remesh_order", "flip_then_smooth")? {
        "flip_then_smooth" => RemeshOrder::FlipThenSmooth,
        "smooth_then_flip" => RemeshOrder::SmoothThenFlip,
        other => {
            return Err(bad(
                "remesh_order",
                format!("'{other}' is not valid; use flip_then_smooth or smooth_then_flip"),
            ));
        }
    };
    if let Some(Val::List(ops)) = d.get("remesh_ops") {
        let mut e = EnabledOps {
            split: false,
            collapse: false,
            flip: false,
            smooth: false,
        };
        for op in ops {
            match as_str("remesh_ops", op)? {
                "split" => e.split = true,
                "collapse" => e.collapse = true,
                "flip" => e.flip = true,
                "smooth" => e.smooth = true,
                other => {
                    return Err(bad(
                        "remesh_ops",
                        format!("'{other}' is not valid; use split, collapse, flip, smooth"),
                    ));
                }
            }
        }
        r.enabled_ops = e;
    } else if let Some(other) = d.get("remesh_ops") {
        return Err(bad(
            "remesh_ops",
            format!("expected a list of names, got {}", other.type_name()),
        ));
    }
    o.remesh = d.bool("remesh", true)?.then_some(r);

    // Stop rules and limits.
    let rel = d.f64("stop_relative", 1e-7)?;
    let window = d
        .get("stop_window")
        .map(|v| as_usize("stop_window", v))
        .transpose()?
        .unwrap_or(10);
    let eps = d.opt_f64("stop_eps")?;
    let stop = d.string("stop", "energy_stall")?;
    one_of("stop", stop, &["energy_stall", "sobolev", "off"])?;
    o.stop = match stop {
        "energy_stall" => StopRule::EnergyStall { relative: rel, window },
        "sobolev" => StopRule::SobolevNorm {
            eps: eps.ok_or_else(|| bad("stop_eps", "required with stop='sobolev'"))?,
        },
        _ => StopRule::Off,
    };
    if let Some(v) = d.get("max_iterations") {
        o.max_iterations = as_usize("max_iterations", v)?;
    }
    if let Some(t) = d.opt_f64("time_limit")? {
        o.time_limit = Some(
            Duration::try_from_secs_f64(t)
                .map_err(|_| bad("time_limit", format!("{t} must be finite and >= 0 seconds")))?,
        );
    }
    if let Some(v) = d.get("stuck_limit") {
        o.stuck_limit = as_usize("stuck_limit", v)?;
    }
    if let Some(v) = d.get("threads") {
        o.threads = Some(as_usize("threads", v)?);
    }
    o.reference_compat = d.bool("reference_compat", false)?;
    d.finish()?;
    o.validate()?;
    Ok(o)
}

// ---------------------------------------------------------------------------------------------
// Options back to values; the defaults
// ---------------------------------------------------------------------------------------------

fn vec3_val(p: [f64; 3]) -> Val {
    Val::floats(&p)
}

fn flat_val(p: &[[f64; 3]]) -> Val {
    Val::Array(p.iter().flatten().copied().collect())
}

/// The options as the dict [`flow_options`] reads: `flow_options(options_val(o))` equals `o`
/// (the penalties and obstacles of a scene travel as opaque handles).
pub fn options_val(o: &FlowOptions) -> Val {
    let mut items: Vec<(&str, Val)> = Vec::new();
    items.push(("alpha", Val::Float(o.exponents.alpha())));
    items.push(("beta", Val::Float(o.exponents.beta())));
    let (energy, coulomb) = match o.energy {
        Energy::TangentPoint => ("tangent_point", None),
        Energy::Willmore => ("willmore", None),
        Energy::Coulomb { alpha } => ("coulomb", Some(alpha)),
    };
    items.push(("energy", Val::s(energy)));
    items.push(("coulomb_alpha", Val::opt_f(coulomb)));
    items.push(("tp_weight", Val::Float(o.tp_weight)));
    let (metric, low_order, l2) = match o.metric {
        Metric::Sobolev { low_order } => ("sobolev", low_order, "mass"),
        Metric::L2(L2Kind::Mass) => ("l2", true, "mass"),
        Metric::L2(_) => ("l2", true, "identity"),
        Metric::H1 => ("h1", true, "mass"),
        _ => ("h2", true, "mass"),
    };
    items.push(("metric", Val::s(metric)));
    items.push(("low_order", Val::Bool(low_order)));
    items.push(("l2_kind", Val::s(l2)));
    let (backend, dense_max, h) = match o.backend {
        Backend::Dense => ("dense", rs::solve::iterative::DEFAULT_DENSE_MAX_FACES, HierarchicalOptions::default()),
        Backend::Hierarchical(h) => ("hierarchical", rs::solve::iterative::DEFAULT_DENSE_MAX_FACES, h),
        Backend::Auto {
            dense_max_faces,
            hierarchical,
        } => ("auto", dense_max_faces, hierarchical),
    };
    items.push(("backend", Val::s(backend)));
    items.push(("dense_max_faces", Val::uint(dense_max)));
    let (eval, bh) = match &o.energy_eval {
        EnergyEval::Exact => ("exact", BhOptions::default()),
        EnergyEval::BarnesHut(b) => ("barnes_hut", *b),
    };
    items.push(("energy_eval", Val::s(eval)));
    items.push(("theta", Val::Float(bh.theta)));
    items.push(("leaf_size", Val::uint(bh.leaf_size)));
    items.push(("chi", Val::Float(h.cluster.chi)));
    items.push(("cluster_leaf_size", Val::uint(h.cluster.leaf_size)));
    items.push(("gmres_restart", Val::uint(h.gmres.restart)));
    items.push(("gmres_tolerance", Val::Float(h.gmres.tolerance)));
    items.push(("gmres_max_iterations", Val::uint(h.gmres.max_iterations)));
    items.push((
        "gmres_side",
        Val::s(if h.gmres.side == GmresSide::Right { "right" } else { "left" }),
    ));
    let (device, cuda) = match h.device {
        Device::Cpu => ("cpu", CudaOptions::default()),
        Device::Cuda(c) => ("cuda", c),
        Device::Auto(c) => ("auto", c),
    };
    items.push(("device", Val::s(device)));
    items.push(("cuda_ordinal", Val::uint(cuda.ordinal)));
    items.push((
        "gpu_determinism",
        Val::s(if cuda.determinism == Determinism::MatchCpu { "match_cpu" } else { "run_to_run" }),
    ));
    let (method, kappa, memory, variant) = match o.method {
        Method::ProjectedGradient => ("projected_gradient", DEFAULT_AQP_KAPPA, DEFAULT_LBFGS_MEMORY, LbfgsVariant::Standard),
        Method::Aqp { kappa } => ("aqp", kappa, DEFAULT_LBFGS_MEMORY, LbfgsVariant::Standard),
        Method::H1Lbfgs { memory, variant } => ("h1_lbfgs", DEFAULT_AQP_KAPPA, memory, variant),
        Method::Bcqn { memory, variant } => (
            "bcqn",
            DEFAULT_AQP_KAPPA,
            memory,
            if variant == BcqnVariant::Standard { LbfgsVariant::Standard } else { LbfgsVariant::Reference },
        ),
    };
    items.push(("method", Val::s(method)));
    items.push(("aqp_kappa", Val::Float(kappa)));
    items.push(("lbfgs_memory", Val::uint(memory)));
    items.push((
        "lbfgs_variant",
        Val::s(if variant == LbfgsVariant::Standard { "standard" } else { "reference" }),
    ));
    // Constraints: None is the default (one barycenter).
    let default_constraints = FlowOptions::default().constraints;
    if o.constraints == default_constraints && o.schedules.is_empty() && o.pin_schedules.is_empty() {
        items.push(("constraints", Val::None));
    } else {
        let list = o
            .constraints
            .iter()
            .enumerate()
            .map(|(i, c)| constraint_val(o, i, c))
            .collect();
        items.push(("constraints", Val::List(list)));
    }
    items.push(("barycenter", Val::None));
    items.push((
        "barycenter_mode",
        Val::s(match o.barycenter_mode {
            None => "auto",
            Some(BarycenterMode::Global) => "global",
            Some(_) => "per_component",
        }),
    ));
    items.push(("barycenter_shift", Val::Bool(o.barycenter_shift)));
    items.push(("pins", Val::None));
    items.push(("pin_targets", Val::None));
    items.push(("pin_move", Val::None));
    items.push(("pin_move_steps", Val::None));
    items.push((
        "schedule_timing",
        Val::s(match o.schedule_timing {
            ScheduleTiming::BeforeProjection => "before_projection",
            _ => "after_projection",
        }),
    ));
    items.push((
        "penalties",
        if o.penalties.is_empty() {
            Val::None
        } else {
            Val::List(o.penalties.iter().cloned().map(Val::Penalty).collect())
        },
    ));
    items.push((
        "ccd_obstacles",
        if o.ccd_obstacles.is_empty() {
            Val::None
        } else {
            Val::List(o.ccd_obstacles.iter().cloned().map(Val::Obstacle).collect())
        },
    ));
    let ls = &o.line_search;
    let (armijo, c) = match ls.armijo {
        ArmijoRule::Standard { c } => ("standard", c),
        ArmijoRule::Reference { c } => ("reference", c),
    };
    items.push(("armijo", Val::s(armijo)));
    items.push(("armijo_c", Val::Float(c)));
    let (start, lengths, factor, tau) = match &ls.start {
        StartRule::GrowPrevious { factor, first } => (
            "grow_previous",
            match **first {
                StartRule::MaxDisplacement { lengths } => lengths,
                _ => 1.0,
            },
            *factor,
            None,
        ),
        StartRule::MaxDisplacement { lengths } => ("max_displacement", *lengths, 2.0, None),
        StartRule::Reference => ("reference", 1.0, 2.0, None),
        StartRule::Absolute { tau } => ("absolute", 1.0, 2.0, Some(*tau)),
    };
    items.push(("start", Val::s(start)));
    items.push(("start_lengths", Val::Float(lengths)));
    items.push(("start_factor", Val::Float(factor)));
    items.push(("start_tau", Val::opt_f(tau)));
    items.push(("tau_min_rel", Val::Float(ls.tau_min_rel)));
    items.push((
        "max_halvings",
        ls.max_halvings.map_or(Val::None, |h| Val::Int(i64::from(h))),
    ));
    let ccd = ls.ccd.unwrap_or_default();
    items.push(("ccd", Val::Bool(ls.ccd.is_some())));
    items.push(("ccd_gap", Val::Float(ccd.gap)));
    items.push(("ccd_max_iterations", Val::Int(i64::from(ccd.max_iterations))));
    items.push(("ccd_close_factor", Val::Float(ccd.close_factor)));
    let (fixed, max) = match o.newton_steps {
        NewtonSteps::Fixed(n) => (Val::Int(i64::from(n)), 5),
        NewtonSteps::UntilTolerance { max } => (Val::None, max),
    };
    items.push(("newton_steps", fixed));
    items.push(("newton_max", Val::Int(i64::from(max))));
    let r = o.remesh.unwrap_or_default();
    items.push(("remesh", Val::Bool(o.remesh.is_some())));
    items.push(("remesh_split", Val::Float(r.split)));
    items.push(("remesh_collapse", Val::Float(r.collapse)));
    items.push(("remesh_rounds", Val::uint(r.rounds)));
    items.push(("remesh_rho", Val::Float(r.rho)));
    items.push(("remesh_flip_cap", Val::uint(r.flip_cap_per_vertex)));
    items.push((
        "remesh_order",
        Val::s(match r.order {
            RemeshOrder::FlipThenSmooth => "flip_then_smooth",
            _ => "smooth_then_flip",
        }),
    ));
    let e = r.enabled_ops;
    let ops: Vec<Val> = [
        (e.split, "split"),
        (e.collapse, "collapse"),
        (e.flip, "flip"),
        (e.smooth, "smooth"),
    ]
    .into_iter()
    .filter(|(on, _)| *on)
    .map(|(_, n)| Val::s(n))
    .collect();
    items.push(("remesh_ops", Val::List(ops)));
    let (stop, rel, window, eps) = match o.stop {
        StopRule::EnergyStall { relative, window } => ("energy_stall", relative, window, None),
        StopRule::SobolevNorm { eps } => ("sobolev", 1e-7, 10, Some(eps)),
        _ => ("off", 1e-7, 10, None),
    };
    items.push(("stop", Val::s(stop)));
    items.push(("stop_relative", Val::Float(rel)));
    items.push(("stop_window", Val::uint(window)));
    items.push(("stop_eps", Val::opt_f(eps)));
    items.push(("max_iterations", Val::uint(o.max_iterations)));
    items.push(("time_limit", Val::opt_f(o.time_limit.map(|t| t.as_secs_f64()))));
    items.push(("stuck_limit", Val::uint(o.stuck_limit)));
    items.push(("threads", o.threads.map_or(Val::None, Val::uint)));
    items.push(("reference_compat", Val::Bool(o.reference_compat)));
    Val::dict(items)
}

fn constraint_val(o: &FlowOptions, i: usize, c: &Constraint) -> Val {
    let comp = |c: &Option<u32>| c.map_or(Val::None, |c| Val::Int(i64::from(c)));
    match c {
        Constraint::Barycenter { component, target } => Val::dict(vec![
            ("type", Val::s("barycenter")),
            ("component", comp(component)),
            ("target", target.map_or(Val::None, vec3_val)),
        ]),
        Constraint::Pins { vertices, targets } => {
            let sch = o.pin_schedules.iter().find(|s| s.constraint == i);
            Val::dict(vec![
                ("type", Val::s("pins")),
                (
                    "vertices",
                    Val::Array(vertices.iter().map(|&v| f64::from(v)).collect()),
                ),
                ("targets", targets.as_ref().map_or(Val::None, |t| flat_val(t))),
                ("move", sch.map_or(Val::None, |s| flat_val(&s.displacement))),
                ("steps", sch.map_or(Val::None, |s| Val::uint(s.steps))),
            ])
        }
        Constraint::Area { component, target } | Constraint::Volume { component, target } => {
            let sch = o.schedules.iter().find(|s| s.constraint == i);
            Val::dict(vec![
                (
                    "type",
                    Val::s(if matches!(c, Constraint::Area { .. }) { "area" } else { "volume" }),
                ),
                ("component", comp(component)),
                ("target", Val::opt_f(*target)),
                ("factor", sch.map_or(Val::None, |s| Val::Float(s.factor))),
                ("offset", sch.map_or(Val::None, |s| Val::Float(s.offset))),
                ("steps", sch.map_or(Val::None, |s| Val::uint(s.steps))),
            ])
        }
    }
}

/// Every option with its default (the library's, read from [`FlowOptions::default`]).
pub fn defaults() -> Val {
    options_val(&FlowOptions::default())
}

// ---------------------------------------------------------------------------------------------
// Reports
// ---------------------------------------------------------------------------------------------

/// The name of a session status.
pub fn status_name(s: &Status) -> &'static str {
    match s {
        Status::Running => "running",
        Status::Converged => "converged",
        Status::MaxIterations => "max_iterations",
        Status::TimeLimit => "time_limit",
        Status::Stuck => "stuck",
        Status::Cancelled => "cancelled",
        Status::Error(_) => "error",
    }
}

/// The report of one step as a dict (`status`: the outcome of the step, `session_status`: that
/// of the session after it).
pub fn step_val(r: &StepReport, session_status: &Status, device: Option<&str>) -> Val {
    let solves = r.linear.solves.len();
    let linear = Val::dict(vec![
        ("solver", Val::s(r.linear.solver)),
        ("dimension", Val::uint(r.linear.dimension)),
        ("solves", Val::uint(solves)),
        ("iterations", Val::uint(r.linear.total_iterations())),
        ("relative_residual", Val::Float(r.linear.max_relative_residual())),
    ]);
    let remesh = r.remesh.as_ref().map_or(Val::None, |m| {
        Val::dict(vec![
            ("splits", Val::Int(m.splits as i64)),
            ("collapses", Val::Int(m.collapses as i64)),
            ("flips", Val::Int(m.flips as i64)),
            ("smoothing_moves", Val::Int(m.smoothing_moves as i64)),
            ("n_vertices", Val::uint(m.n_vertices)),
            ("n_faces", Val::uint(m.n_faces)),
            ("seconds", Val::Float(m.seconds)),
        ])
    });
    let t = &r.seconds;
    let seconds = Val::dict(vec![
        ("bvh", Val::Float(t.bvh)),
        ("energy", Val::Float(t.energy)),
        ("factor", Val::Float(t.factor)),
        ("solve", Val::Float(t.solve)),
        ("line_search", Val::Float(t.line_search)),
        ("ccd", Val::Float(t.ccd)),
        ("projection", Val::Float(t.projection)),
        ("remesh", Val::Float(t.remesh)),
        ("total", Val::Float(t.total)),
    ]);
    Val::dict(vec![
        ("iteration", Val::uint(r.iteration)),
        (
            "status",
            Val::s(match r.status {
                StepStatus::Accepted => "accepted",
                _ => "stuck",
            }),
        ),
        ("session_status", Val::s(status_name(session_status))),
        ("device", Val::s(if device.is_some() { "cuda" } else { "cpu" })),
        ("device_name", device.map_or(Val::None, Val::s)),
        ("n_vertices", Val::uint(r.vertices)),
        ("n_faces", Val::uint(r.faces)),
        ("energy_before", Val::Float(r.energy_before)),
        ("energy_after", Val::Float(r.energy_after)),
        ("tangent_point_after", Val::Float(r.energy_term_after)),
        ("penalties_after", Val::floats(&r.penalties_after)),
        ("gradient_norm", Val::Float(r.gradient_norm)),
        ("sobolev_norm_sq", Val::Float(r.sobolev_norm_sq)),
        ("tau", Val::Float(r.tau)),
        ("tau0", Val::Float(r.tau0)),
        ("tau_ccd", Val::opt_f(r.tau_ccd)),
        ("halvings", Val::Int(i64::from(r.halvings))),
        ("newton_steps", Val::Int(i64::from(r.newton_steps))),
        ("projection", Val::s(&format!("{:?}", r.projection_status).to_lowercase())),
        ("max_violation", Val::Float(r.max_violation)),
        (
            "dropped_rows",
            Val::List(r.dropped_rows.iter().map(|&i| Val::uint(i)).collect()),
        ),
        ("linear", linear),
        ("remesh", remesh),
        ("seconds", seconds),
    ])
}

fn cancelled_val(iteration: usize, energy: f64, n_vertices: usize, n_faces: usize) -> Val {
    Val::dict(vec![
        ("iteration", Val::uint(iteration)),
        ("status", Val::s("cancelled")),
        ("session_status", Val::s("cancelled")),
        ("n_vertices", Val::uint(n_vertices)),
        ("n_faces", Val::uint(n_faces)),
        ("energy_before", Val::Float(energy)),
        ("energy_after", Val::Float(energy)),
    ])
}

// ---------------------------------------------------------------------------------------------
// The engine
// ---------------------------------------------------------------------------------------------

/// The last completed state of a session, readable without waiting for a step.
#[derive(Debug, Clone)]
pub struct Snapshot {
    /// Positions.
    pub positions: Vec<[f64; 3]>,
    /// Faces (they change with remeshing).
    pub faces: Vec<[u32; 3]>,
    /// For each face, the face of the initial mesh it descends from.
    pub face_origin: Vec<u32>,
    /// The pinned vertices (current ids).
    pub pinned: Vec<u32>,
    /// Steps taken.
    pub iteration: usize,
    /// `F`.
    pub energy: f64,
    /// Status name.
    pub status: &'static str,
    /// The message of the error when the status is `"error"`.
    pub error: Option<String>,
    /// `l0`, the mean edge length of the start.
    pub mean_edge_length0: f64,
    /// The report of the last step.
    pub last_report: Option<Val>,
}

/// Why a [`Engine::run`] ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunEnd {
    /// `max_steps` steps were taken.
    MaxSteps,
    /// The session reached a final status.
    Final,
    /// The callback asked to stop.
    Callback,
    /// `cancel()` was called.
    Cancelled,
    /// The time limit of the call.
    TimeLimit,
}

impl RunEnd {
    /// The name.
    pub fn name(self) -> &'static str {
        match self {
            RunEnd::MaxSteps => "max_steps",
            RunEnd::Final => "final",
            RunEnd::Callback => "callback",
            RunEnd::Cancelled => "cancelled",
            RunEnd::TimeLimit => "time_limit",
        }
    }
}

/// The result of [`Engine::run`].
#[derive(Debug, Clone)]
pub struct RunSummary {
    /// Why the run ended.
    pub end: RunEnd,
    /// The reports of the steps of this run.
    pub steps: Vec<Val>,
    /// `F` at the start of the run.
    pub energy_before: f64,
    /// Wall time.
    pub seconds: f64,
}

/// A session that several threads can use (see the module docs).
pub struct Engine {
    session: Mutex<Session>,
    cancel: Mutex<Arc<AtomicBool>>,
    snapshot: Mutex<Snapshot>,
    start: (Surface, Vec<[f64; 3]>),
    options: FlowOptions,
}

fn poisoned() -> SError {
    SError::Solve("internal error: a previous call panicked inside this session".into())
}

fn make_snapshot(session: &Session) -> Snapshot {
    Snapshot {
        positions: session.positions().to_vec(),
        faces: session.surface().faces().to_vec(),
        face_origin: session.face_origin().to_vec(),
        pinned: session.pinned().to_vec(),
        iteration: session.iteration(),
        energy: session.energy(),
        status: status_name(&session.status()),
        error: None,
        mean_edge_length0: session.mean_edge_length0(),
        last_report: None,
    }
}

/// Builds the surface of the input mesh.
///
/// # Errors
/// [`SError::Input`] with the library's message.
pub fn make_surface(positions: &[[f64; 3]], faces: Vec<[u32; 3]>) -> Result<Surface> {
    let s = Surface::new(positions.len(), faces)?;
    s.validate_positions(positions)?;
    Ok(s)
}

impl Engine {
    /// Starts a session on the mesh.
    ///
    /// # Errors
    /// [`SError::Input`]: the mesh, the positions or an option is invalid (the library's message).
    pub fn new(
        positions: Vec<[f64; 3]>,
        faces: Vec<[u32; 3]>,
        options: &[(String, Val)],
    ) -> Result<Engine> {
        let surface = make_surface(&positions, faces)?;
        let flow = flow_options(options, &surface, &positions)?;
        let session = Session::new(surface.clone(), positions.clone(), flow.clone())?;
        let snapshot = make_snapshot(&session);
        let cancel = session.cancel_flag();
        Ok(Engine {
            session: Mutex::new(session),
            cancel: Mutex::new(cancel),
            snapshot: Mutex::new(snapshot),
            start: (surface, positions),
            options: flow,
        })
    }

    fn session(&self) -> Result<MutexGuard<'_, Session>> {
        self.session.lock().map_err(|_| poisoned())
    }

    fn snap(&self) -> Result<MutexGuard<'_, Snapshot>> {
        self.snapshot.lock().map_err(|_| poisoned())
    }

    fn cancel_flag(&self) -> Result<Arc<AtomicBool>> {
        Ok(self.cancel.lock().map_err(|_| poisoned())?.clone())
    }

    /// The last completed state. Never waits for a running step.
    ///
    /// # Errors
    /// [`SError::Solve`] if a call panicked inside the session.
    pub fn snapshot(&self) -> Result<Snapshot> {
        Ok(self.snap()?.clone())
    }

    /// Asks the running step (between line-search trials) or the next one to stop. The snapshot
    /// reads `"cancelled"` at once. Never waits.
    ///
    /// # Errors
    /// [`SError::Solve`] if a call panicked inside the session.
    pub fn cancel(&self) -> Result<()> {
        self.cancel_flag()?.store(true, Ordering::Relaxed);
        self.snap()?.status = "cancelled";
        Ok(())
    }

    /// Whether `cancel` was called since the start or the last `reset`.
    ///
    /// # Errors
    /// [`SError::Solve`] if a call panicked inside the session.
    pub fn is_cancelled(&self) -> Result<bool> {
        Ok(self.cancel_flag()?.load(Ordering::Relaxed))
    }

    /// Starts again from the initial mesh with the same options (clears the cancel flag).
    ///
    /// # Errors
    /// The errors of [`Engine::new`].
    pub fn reset(&self) -> Result<()> {
        let session = Session::new(
            self.start.0.clone(),
            self.start.1.clone(),
            self.options.clone(),
        )?;
        let snapshot = make_snapshot(&session);
        let mut s = self.session()?;
        *self.cancel.lock().map_err(|_| poisoned())? = session.cancel_flag();
        *s = session;
        *self.snap()? = snapshot;
        Ok(())
    }

    /// One step. A cancelled session does not step and returns a report with status
    /// `"cancelled"`.
    ///
    /// # Errors
    /// [`SError::Solve`] / [`SError::Input`] for a failed step; the session's status is then
    /// `"error"` and the mesh that of the last completed step.
    pub fn step(&self) -> Result<Val> {
        let mut s = self.session()?;
        let finish_cancelled = |s: &Session| {
            self.snap()?.status = "cancelled";
            Ok(cancelled_val(
                s.iteration(),
                s.energy(),
                s.surface().n_vertices(),
                s.surface().n_faces(),
            ))
        };
        if self.cancel_flag()?.load(Ordering::Relaxed) {
            return finish_cancelled(&s);
        }
        match s.step() {
            Ok(report) => {
                let status = s.status();
                let val = step_val(&report, &status, s.gpu_device_name().as_deref());
                let mut snap = make_snapshot(&s);
                snap.last_report = Some(val.clone());
                *self.snap()? = snap;
                Ok(val)
            }
            Err(RsError::Cancelled) => finish_cancelled(&s),
            Err(e) => {
                let status = s.status();
                let mut snap = self.snap()?;
                snap.status = status_name(&status);
                snap.error = Some(e.to_string());
                Err(e.into())
            }
        }
    }

    /// Steps until the session reaches a final status (converged, stuck, `max_iterations`
    /// steps in all), `max_steps` steps of this call are taken, `cancel` is called,
    /// `time_limit` seconds have passed, or `progress` returns `false`.
    ///
    /// Between steps the lock is free, so `snapshot` and `cancel` work during a run.
    ///
    /// # Errors
    /// The errors of [`Engine::step`].
    pub fn run(
        &self,
        max_steps: Option<usize>,
        time_limit: Option<f64>,
        mut progress: Option<&mut dyn FnMut(&Val) -> bool>,
    ) -> Result<RunSummary> {
        let t0 = Instant::now();
        let energy_before = self.snap()?.energy;
        let mut steps = Vec::new();
        let end = loop {
            if max_steps.is_some_and(|n| steps.len() >= n) {
                break RunEnd::MaxSteps;
            }
            if time_limit.is_some_and(|l| t0.elapsed().as_secs_f64() >= l) {
                break RunEnd::TimeLimit;
            }
            if self.is_cancelled()? {
                break RunEnd::Cancelled;
            }
            {
                let mut snap = self.snap()?;
                if snap.status != "running" {
                    break RunEnd::Final;
                }
                if snap.iteration >= self.options.max_iterations {
                    snap.status = "max_iterations";
                    break RunEnd::Final;
                }
            }
            let report = self.step()?;
            if report.get("status") == Some(&Val::s("cancelled")) {
                break RunEnd::Cancelled;
            }
            let go = progress.as_mut().is_none_or(|p| p(&report));
            steps.push(report);
            if !go {
                break RunEnd::Callback;
            }
        };
        Ok(RunSummary {
            end,
            steps,
            energy_before,
            seconds: t0.elapsed().as_secs_f64(),
        })
    }
}

impl PartialEq for Val {
    fn eq(&self, other: &Val) -> bool {
        match (self, other) {
            (Val::None, Val::None) => true,
            (Val::Bool(a), Val::Bool(b)) => a == b,
            (Val::Int(a), Val::Int(b)) => a == b,
            (Val::Float(a), Val::Float(b)) => a == b,
            (Val::Str(a), Val::Str(b)) => a == b,
            (Val::List(a), Val::List(b)) => a == b,
            (Val::Array(a), Val::Array(b)) => a == b,
            (Val::Dict(a), Val::Dict(b)) => a == b,
            _ => false,
        }
    }
}

// ---------------------------------------------------------------------------------------------
// One-shot flow, energy, scenes
// ---------------------------------------------------------------------------------------------

/// The result of [`flow_once`].
#[derive(Debug, Clone)]
pub struct FlowOutput {
    /// Final positions.
    pub positions: Vec<[f64; 3]>,
    /// Final faces.
    pub faces: Vec<[u32; 3]>,
    /// For each final face, the face of the input it descends from.
    pub face_origin: Vec<u32>,
    /// The pinned vertices (final ids).
    pub pinned: Vec<u32>,
    /// The report (see the glue for the keys).
    pub report: Val,
}

/// `Session::new` and `run` in one call, with reports as [`Val`]s. `progress` gets the report of
/// every step and returns `false` to stop (the status is then `"cancelled"`).
///
/// # Errors
/// As [`Engine::new`].
pub fn flow_once(
    positions: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    options: &[(String, Val)],
    mut progress: Option<&mut dyn FnMut(&Val) -> bool>,
) -> Result<FlowOutput> {
    let surface = make_surface(&positions, faces)?;
    let flow = flow_options(options, &surface, &positions)?;
    let mut session = Session::new(surface, positions, flow)?;
    let device = session.gpu_device_name();
    let mut cb = |r: &StepReport| {
        let go = progress
            .as_mut()
            .is_none_or(|p| p(&step_val(r, &Status::Running, device.as_deref())));
        if go { Control::Continue } else { Control::Stop }
    };
    let r = session.run(Some(&mut cb))?;
    if let Status::Error(e) = &r.status {
        return Err(e.clone().into());
    }
    let report = Val::dict(vec![
        ("status", Val::s(status_name(&r.status))),
        ("iterations", Val::uint(r.iterations)),
        (
            "energy_initial",
            r.steps
                .first()
                .map_or(Val::None, |s| Val::Float(s.energy_before)),
        ),
        ("energy_final", Val::Float(r.final_energy)),
        ("seconds", Val::Float(r.seconds)),
        (
            "steps",
            Val::List(r.steps.iter().map(|s| step_val(s, &r.status, device.as_deref())).collect()),
        ),
    ]);
    Ok(FlowOutput {
        positions: session.positions().to_vec(),
        faces: session.surface().faces().to_vec(),
        face_origin: session.face_origin().to_vec(),
        pinned: session.pinned().to_vec(),
        report,
    })
}

/// What this build and machine offer for the `device` option: `cuda_compiled` (the wheel was
/// built with the `cuda` feature), `available` (a device can be opened), `devices` (a list of
/// `{ordinal, name, compute_capability, memory_bytes}`) and `reason` (why none is usable, or
/// `None`). Opens the driver only; compiles nothing.
pub fn gpu_info() -> Val {
    let compiled = rs::gpu::cuda_compiled();
    let (devices, reason) = match rs::gpu::devices() {
        Ok(d) if d.is_empty() => (Vec::new(), Some("the CUDA driver reports no device".to_owned())),
        Ok(d) => (d, None),
        Err(e) => (Vec::new(), Some(e.to_string())),
    };
    Val::dict(vec![
        ("cuda_compiled", Val::Bool(compiled)),
        ("available", Val::Bool(compiled && !devices.is_empty())),
        (
            "devices",
            Val::List(
                devices
                    .iter()
                    .map(|d| {
                        Val::dict(vec![
                            ("ordinal", Val::uint(d.ordinal)),
                            ("name", Val::s(&d.name)),
                            (
                                "compute_capability",
                                Val::s(&format!("{}.{}", d.compute_capability.0, d.compute_capability.1)),
                            ),
                            ("memory_bytes", Val::Int(d.total_memory as i64)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("reason", reason.map_or(Val::None, |r| Val::s(&r))),
    ])
}

/// The tangent-point energy with the penalties of the options, and the differential of the
/// total: `(differential, report)`; the report has `energy`, `total`, `penalties` (a list of
/// `{name, value}` in the order given), `subcritical`, `pairs_exact`, `nodes_lumped`, `seconds`.
/// Reads the options `alpha`, `beta`, `energy_eval`, `theta`, `leaf_size`, `penalties`,
/// `threads`; every other option is validated and ignored.
///
/// # Errors
/// [`SError::Input`] for a bad mesh or option.
pub fn surface_energy(
    positions: &[[f64; 3]],
    faces: Vec<[u32; 3]>,
    options: &[(String, Val)],
) -> Result<(Vec<[f64; 3]>, Val)> {
    let surface = make_surface(positions, faces)?;
    let o = flow_options(options, &surface, positions)?;
    let r = rs::energy(
        &surface,
        positions,
        o.exponents,
        &o.energy_eval,
        &o.penalties,
        o.threads,
    )?;
    let report = Val::dict(vec![
        ("energy", Val::Float(r.energy)),
        ("total", Val::Float(r.total())),
        (
            "penalties",
            Val::List(
                o.penalties
                    .iter()
                    .zip(&r.penalties)
                    .map(|(p, &v)| {
                        Val::dict(vec![("name", Val::s(p.name())), ("value", Val::Float(v))])
                    })
                    .collect(),
            ),
        ),
        ("subcritical", Val::Bool(r.subcritical)),
        ("pairs_exact", Val::Int(r.pairs_exact as i64)),
        ("nodes_lumped", Val::Int(r.nodes_lumped as i64)),
        ("seconds", Val::Float(r.seconds)),
    ]);
    Ok((r.differential, report))
}

/// A scene file loaded: `(positions, faces, options)`; the options are the dict
/// [`flow_options`] reads (the scene's penalties and obstacles as opaque handles).
///
/// # Errors
/// [`SError::Input`] for a file that cannot be read or parsed, or for a scene feature the
/// library does not support (normal constraints).
pub fn scene(path: &str) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>, Val)> {
    let spec = rs::SceneSpec::load(path)?;
    let o = rs::scene_options(&spec)?;
    Ok((spec.positions.clone(), spec.surface.faces().to_vec(), options_val(&o)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rs::generate::primitives;

    fn dict(items: Vec<(&str, Val)>) -> Vec<(String, Val)> {
        items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect()
    }

    fn mesh() -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        let m = primitives::torus(12, 8, 2.0, 0.7).unwrap();
        (m.positions, m.faces)
    }

    fn items(v: &Val) -> Vec<(String, Val)> {
        match v {
            Val::Dict(d) => d.clone(),
            _ => panic!("not a dict"),
        }
    }

    #[test]
    fn defaults_cover_all_options_and_equal_the_library() {
        let d = defaults();
        let keys: Vec<&str> = items(&d).iter().map(|(k, _)| k.as_str()).collect::<Vec<_>>().into_iter().map(|k| OPTION_NAMES.iter().find(|n| **n == k).copied().unwrap()).collect();
        for n in OPTION_NAMES {
            assert!(keys.contains(n), "defaults lack {n}");
        }
        assert_eq!(keys.len(), OPTION_NAMES.len());
        let (p, f) = mesh();
        let s = make_surface(&p, f).unwrap();
        let o = flow_options(&items(&d), &s, &p).unwrap();
        assert_eq!(format!("{o:?}"), format!("{:?}", FlowOptions::default()));
    }

    #[test]
    fn unknown_option_is_named() {
        let (p, f) = mesh();
        let e = Engine::new(p, f, &dict(vec![("alfa", Val::Float(1.0))])).err().unwrap();
        assert!(e.to_string().contains("alfa"), "{e}");
    }

    #[test]
    fn bad_values_carry_the_library_message() {
        let (p, f) = mesh();
        let s = make_surface(&p, f).unwrap();
        let e = flow_options(&dict(vec![("alpha", Val::Float(0.5))]), &s, &p).unwrap_err();
        assert!(e.to_string().to_lowercase().contains("alpha"), "{e}");
        let e = flow_options(&dict(vec![("remesh_split", Val::Float(-1.0))]), &s, &p).unwrap_err();
        assert!(e.to_string().contains("split"), "{e}");
        let e = flow_options(&dict(vec![("metric", Val::s("l3"))]), &s, &p).unwrap_err();
        assert!(e.to_string().contains("metric"), "{e}");
    }

    #[test]
    fn options_round_trip() {
        let (p, f) = mesh();
        let s = make_surface(&p, f).unwrap();
        let opts = dict(vec![
            ("metric", Val::s("h1")),
            ("remesh", Val::Bool(false)),
            ("ccd", Val::Bool(false)),
            ("pins", Val::List(vec![Val::Int(0), Val::Int(3)])),
            ("pin_move", Val::floats(&[0.0, 0.0, 0.5])),
            ("pin_move_steps", Val::Int(10)),
            (
                "constraints",
                Val::List(vec![Val::dict(vec![
                    ("type", Val::s("volume")),
                    ("factor", Val::Float(1.5)),
                    ("steps", Val::Int(20)),
                ])]),
            ),
            (
                "penalties",
                Val::List(vec![
                    Val::dict(vec![
                        ("type", Val::s("implicit_obstacle")),
                        (
                            "surface",
                            Val::dict(vec![("type", Val::s("sphere")), ("radius", Val::Float(0.5))]),
                        ),
                    ]),
                    Val::dict(vec![("type", Val::s("area_deviation"))]),
                ]),
            ),
        ]);
        let o = flow_options(&opts, &s, &p).unwrap();
        assert_eq!(o.pin_schedules.len(), 1);
        assert_eq!(o.schedules.len(), 1);
        assert_eq!(o.penalties.len(), 2);
        let o2 = flow_options(&items(&options_val(&o)), &s, &p).unwrap();
        assert_eq!(format!("{o:?}"), format!("{o2:?}"));
    }

    #[test]
    fn engine_steps_equal_the_one_call_flow() {
        let (p, f) = mesh();
        let opts = dict(vec![("max_iterations", Val::Int(4)), ("threads", Val::Int(2))]);
        let e = Engine::new(p.clone(), f.clone(), &opts).unwrap();
        for _ in 0..4 {
            e.step().unwrap();
        }
        let out = flow_once(p, f, &opts, None).unwrap();
        let snap = e.snapshot().unwrap();
        assert_eq!(snap.positions, out.positions);
        assert_eq!(snap.faces, out.faces);
        assert_eq!(snap.face_origin, out.face_origin);
    }

    #[test]
    fn cancel_makes_steps_report_cancelled() {
        let (p, f) = mesh();
        let e = Engine::new(p, f, &[]).unwrap();
        e.cancel().unwrap();
        let r = e.step().unwrap();
        assert_eq!(r.get("status"), Some(&Val::s("cancelled")));
        e.reset().unwrap();
        assert_eq!(e.snapshot().unwrap().status, "running");
    }

    #[test]
    fn energy_report_lists_penalties() {
        let (p, f) = mesh();
        let (d, r) = surface_energy(
            &p,
            f,
            &dict(vec![(
                "penalties",
                Val::List(vec![Val::dict(vec![("type", Val::s("willmore"))])]),
            )]),
        )
        .unwrap();
        assert_eq!(d.len(), p.len());
        assert!(matches!(r.get("energy"), Some(Val::Float(x)) if *x > 0.0));
    }
}
