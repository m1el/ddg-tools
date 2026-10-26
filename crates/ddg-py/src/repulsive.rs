//! Repulsive Curves: the binding of `ddg-repulsive-curves`.
//!
//! Yu, Schumacher, Crane, *Repulsive Curves*, ACM Transactions on Graphics 40(2), 2021: the
//! tangent-point energy of a curve network, its fractional Sobolev descent, constraints and
//! penalties. This module holds everything that is not PyO3 (the Python glue is
//! `python/repulsive.rs`):
//!
//! * [`RepulsiveOptions`]: the options of `RepulsiveSession` and `repulsive_flow` as plain
//!   Rust, their defaults (equal to those of the library, checked by the contract tests below),
//!   and [`RepulsiveOptions::flow_options`], which validates them and builds the library's
//!   [`FlowOptions`] (constraints, penalties, obstacles, back end);
//! * [`Engine`]: a session that can be driven from several threads. `step` takes the session's
//!   lock; `snapshot`, `cancel` and the property getters never wait for a running step (the last
//!   completed state is stored separately, and the cancel flag is an atomic). This is what a
//!   modal Blender operator needs: a worker thread steps, the timer reads the snapshot;
//! * [`flow_once`]: the one-shot flow;
//! * [`Val`]: a small value tree for the reports, turned into dicts by the glue.
//!
//! # Subdivision
//!
//! The library subdivides every edge when the mean edge length exceeds twice its start value
//! (at most `max_subdivisions` times). Old vertices keep their indices; the midpoint of edge
//! `I` becomes vertex `|V| + I`, and edge `I` becomes edges `I` and `|E| + I`
//! (`ddg_repulsive_curves::SubdivisionMap`). [`SourceMap`] follows that through any number of
//! subdivisions: `vertex_source[v]` is the vertex of the input that `v` is (`-1` for a vertex
//! created by subdivision) and `edge_source[e]` the input edge that `e` descends from.
//!
//! # Errors
//!
//! Everything the caller can get wrong is a [`RepError::Input`] carrying the library's message
//! (the glue raises `ValueError`); a failed solve is [`RepError::Solve`].

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};
use std::time::Instant;

use ddg_repulsive_curves as rc;
use ddg_repulsive_curves::energy::exact::make_pool;
use ddg_repulsive_curves::multigrid::{MgOptions, MultigridSolvers};
use ddg_repulsive_curves::surfaces::{
    CircularField, ConstantField, DoubleTorus, Metaballs, Plane, RoundedBox, SmoothUnion, Sphere,
    SphereTangentField, Torus, Union, VectorField,
};
use rc::{
    Axis, BhOptions, Constraint, Control, CurveNetwork, EnergyEval, Exponents, FlowOptions,
    GradNorm, GramQuadrature, ImplicitSurface, LengthDifference, LengthSchedule, LineSearch,
    Metric, Obstacle, PenaltyTerm, ProjectionOptions, RcError, Session, Solver, StartRule, Status,
    StepReport, StepStatus, StopRule, SurfacePotential, TangentSite, Tolerance, TotalLengthPenalty,
};

// ---------------------------------------------------------------------------------------------
// Errors and values
// ---------------------------------------------------------------------------------------------

/// What can go wrong.
#[derive(Debug, Clone, PartialEq)]
pub enum RepError {
    /// The caller's input or options are wrong; the message is the library's, or names the
    /// option. Raised as `ValueError` (and `ddg_tools.InputError`).
    Input(String),
    /// The input was accepted but a computation failed (singular system, non-finite energy).
    Solve(String),
}

impl std::fmt::Display for RepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RepError::Input(m) | RepError::Solve(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for RepError {}

impl From<RcError> for RepError {
    fn from(e: RcError) -> Self {
        match e {
            RcError::SingularMetric { .. }
            | RcError::SingularSaddle { .. }
            | RcError::MultigridDiverged { .. }
            | RcError::NonFiniteEnergy { .. }
            | RcError::Cancelled => RepError::Solve(e.to_string()),
            e => RepError::Input(e.to_string()),
        }
    }
}

fn input(name: &str, message: impl std::fmt::Display) -> RepError {
    RepError::Input(format!("{name}: {message}"))
}

/// A value of a report: the glue turns it into `None`, `bool`, `int`, `float`, `str`, `list`
/// or `dict`.
#[derive(Debug, Clone, PartialEq)]
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
    /// `dict` (insertion order).
    Dict(Vec<(&'static str, Val)>),
}

impl Val {
    /// The entry `key` of a dict.
    pub fn get(&self, key: &str) -> Option<&Val> {
        match self {
            Val::Dict(d) => d.iter().find(|(k, _)| *k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    fn opt_f(x: Option<f64>) -> Val {
        x.map_or(Val::None, Val::Float)
    }

    fn uint(x: usize) -> Val {
        Val::Int(x as i64)
    }
}

// ---------------------------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------------------------

/// The names of all options, in the order of the `.pyi` stub.
pub const OPTION_NAMES: &[&str] = &[
    "alpha",
    "beta",
    "metric",
    "quadrature",
    "length",
    "grow_factor",
    "grow_rate",
    "max_subdivisions",
    "barycenter",
    "pinned",
    "tangent_vertices",
    "tangent_edges",
    "surface",
    "planar_axis",
    "obstacles",
    "obstacle_weight",
    "obstacle_exponent",
    "obstacle_theta",
    "total_length_weight",
    "length_difference_weight",
    "field",
    "field_weight",
    "line_search",
    "line_search_start",
    "armijo",
    "shrink",
    "max_halvings",
    "gap",
    "backend",
    "bh_theta_x",
    "bh_theta_t",
    "bh_leaf_size",
    "threads",
    "projection_tolerance",
    "projection_absolute",
    "projection_iterations",
    "stop",
    "stop_threshold",
    "stall_window",
    "gradient_norm",
    "max_steps",
    "time_limit",
];

/// An obstacle: vertices and triangles.
pub type ObstacleSoup = (Vec<[f64; 3]>, Vec<[u32; 3]>);

/// Which lengths the flow keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthMode {
    /// Every edge keeps its length (the default).
    Edges,
    /// The total length is kept; edges may trade length.
    Total,
    /// Nothing is kept (only the barycenter is fixed).
    Free,
}

/// An implicit surface by name (`f(x) = 0`), for the surface constraint.
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceSpec {
    /// `|x - center| - radius`.
    Sphere {
        /// Centre.
        center: [f64; 3],
        /// Radius.
        radius: f64,
    },
    /// Torus about the `y` axis through `center`.
    Torus {
        /// Centre.
        center: [f64; 3],
        /// Major radius.
        major: f64,
        /// Minor radius.
        minor: f64,
    },
    /// The plane through `point` with normal `normal`.
    Plane {
        /// A point of the plane.
        point: [f64; 3],
        /// The normal.
        normal: [f64; 3],
    },
    /// The double torus of the reference, with `radius2` the squared radius.
    DoubleTorus {
        /// The squared radius.
        radius2: f64,
    },
    /// A box with rounded edges.
    RoundedBox {
        /// Centre.
        center: [f64; 3],
        /// Half extents.
        half: [f64; 3],
        /// Rounding radius.
        radius: f64,
    },
    /// Metaballs `T - sum a exp(-|x-c|^2/b^2)`.
    Metaballs {
        /// `(centre, amplitude, width)`.
        balls: Vec<([f64; 3], f64, f64)>,
        /// The threshold `T`.
        threshold: f64,
    },
    /// The union (pointwise minimum) of surfaces.
    Union(Vec<SurfaceSpec>),
    /// The smooth union of two surfaces with blend width `k`.
    SmoothUnion {
        /// First part.
        a: Box<SurfaceSpec>,
        /// Second part.
        b: Box<SurfaceSpec>,
        /// Blend width.
        k: f64,
    },
}

impl SurfaceSpec {
    /// The surface.
    ///
    /// # Errors
    /// [`RepError::Input`] with the library's message.
    pub fn build(&self) -> Result<Arc<dyn ImplicitSurface>, RepError> {
        Ok(match self {
            SurfaceSpec::Sphere { center, radius } => Arc::new(Sphere::new(*center, *radius)?),
            SurfaceSpec::Torus {
                center,
                major,
                minor,
            } => Arc::new(Torus::new(*center, *major, *minor)?),
            SurfaceSpec::Plane { point, normal } => Arc::new(Plane::new(*point, *normal)?),
            SurfaceSpec::DoubleTorus { radius2 } => Arc::new(DoubleTorus::new(*radius2)?),
            SurfaceSpec::RoundedBox {
                center,
                half,
                radius,
            } => Arc::new(RoundedBox::new(*center, *half, *radius)?),
            SurfaceSpec::Metaballs { balls, threshold } => {
                Arc::new(Metaballs::new(balls.clone(), *threshold)?)
            }
            SurfaceSpec::Union(parts) => {
                let parts = parts
                    .iter()
                    .map(SurfaceSpec::build)
                    .collect::<Result<Vec<_>, _>>()?;
                Arc::new(Union::new(parts)?)
            }
            SurfaceSpec::SmoothUnion { a, b, k } => {
                Arc::new(SmoothUnion::new(a.build()?, b.build()?, *k)?)
            }
        })
    }
}

/// A vector field by name, for the field-alignment penalty.
#[derive(Debug, Clone, PartialEq)]
pub enum FieldSpec {
    /// The constant field `vector` (divided by its length with `normalize`).
    Constant {
        /// The vector.
        vector: [f64; 3],
        /// Use `vector / |vector|`.
        normalize: bool,
    },
    /// The unit field circulating around the `y` axis.
    Circular,
    /// The reference's `interesting` field on the unit sphere.
    SphereTangent,
}

impl FieldSpec {
    /// The field.
    ///
    /// # Errors
    /// [`RepError::Input`] with the library's message.
    pub fn build(&self) -> Result<Arc<dyn VectorField>, RepError> {
        Ok(match self {
            FieldSpec::Constant { vector, normalize } => {
                if *normalize {
                    Arc::new(ConstantField::unit(*vector)?)
                } else {
                    Arc::new(ConstantField::new(*vector)?)
                }
            }
            FieldSpec::Circular => Arc::new(CircularField),
            FieldSpec::SphereTangent => Arc::new(SphereTangentField),
        })
    }
}

/// The options of `RepulsiveSession` and `repulsive_flow`.
#[derive(Debug, Clone, PartialEq)]
pub struct RepulsiveOptions {
    /// `alpha` of the tangent-point energy (library default 3).
    pub alpha: f64,
    /// `beta` (library default 6).
    pub beta: f64,
    /// `"sobolev"`, `"l2"`, `"h1"`, `"h2"` (the last three are for comparisons).
    pub metric: String,
    /// `"midpoint"` or `"trapezoid"` (the Sobolev metric's quadrature).
    pub quadrature: String,
    /// Which lengths are kept.
    pub length: LengthMode,
    /// The final length relative to the start (`1`: no schedule).
    pub grow_factor: f64,
    /// Progress of the schedule per step, as a fraction of the start length.
    pub grow_rate: f64,
    /// Maximal number of subdivisions of the whole network.
    pub max_subdivisions: u32,
    /// Fix the barycenter; `None`: yes unless a vertex is pinned.
    pub barycenter: Option<bool>,
    /// Vertices that do not move.
    pub pinned: Vec<u32>,
    /// Vertices (degree 1 or 2) whose tangent is kept.
    pub tangent_vertices: Vec<u32>,
    /// Edges whose tangent is kept.
    pub tangent_edges: Vec<u32>,
    /// Keep the vertices on a surface: the surface and the vertices (`None`: all).
    pub surface: Option<(SurfaceSpec, Option<Vec<u32>>)>,
    /// Planar mode: the held axis.
    pub planar_axis: Option<Axis>,
    /// Obstacles: `(vertices, triangles)`. The curve never crosses them (collision-guarded
    /// line search).
    pub obstacles: Vec<ObstacleSoup>,
    /// Weight of the repulsive potential of the obstacles (0: none, only the collision check).
    pub obstacle_weight: f64,
    /// Exponent `p` of that potential (`None`: `beta - alpha`).
    pub obstacle_exponent: Option<f64>,
    /// Far-field opening of that potential (`None`: the library's).
    pub obstacle_theta: Option<f64>,
    /// Weight of the total-length penalty (0: none).
    pub total_length_weight: f64,
    /// Weight of the length-difference penalty (0: none).
    pub length_difference_weight: f64,
    /// A vector field the tangents are pushed to align with.
    pub field: Option<FieldSpec>,
    /// Weight of the field penalty.
    pub field_weight: f64,
    /// `"collision"` or `"armijo"`.
    pub line_search: String,
    /// `"two_thirds"` or `"grow_previous"` (the start of the collision-guarded search).
    pub line_search_start: String,
    /// Armijo constant.
    pub armijo: f64,
    /// Shrink factor.
    pub shrink: f64,
    /// Maximal number of shrinks.
    pub max_halvings: u32,
    /// The collision-guarded search's gap.
    pub gap: f64,
    /// `"auto"`, `"dense"`, `"barnes_hut"`, `"multigrid"`, `"accelerated"`.
    pub backend: String,
    /// Barnes-Hut `theta_x`.
    pub bh_theta_x: f64,
    /// Barnes-Hut `theta_t`.
    pub bh_theta_t: f64,
    /// Barnes-Hut leaf size.
    pub bh_leaf_size: usize,
    /// Threads (`None`: all cores).
    pub threads: Option<usize>,
    /// Tolerance of the constraint projection.
    pub projection_tolerance: f64,
    /// The tolerance is an absolute length, not a fraction of the mean edge length.
    pub projection_absolute: bool,
    /// Maximal number of chord iterations of the projection.
    pub projection_iterations: u32,
    /// `"gradient"`, `"relative"` or `"stall"`.
    pub stop: String,
    /// The threshold of the rule (`None`: 1e-4 for gradient, 1e-3 for relative; required for stall).
    pub stop_threshold: Option<f64>,
    /// Window of the stall rule.
    pub stall_window: usize,
    /// `"euclidean"` or `"mass"`: the norm of the gradient rule.
    pub gradient_norm: String,
    /// Maximal number of steps of the session.
    pub max_steps: usize,
    /// Wall-time limit of a run, in seconds.
    pub time_limit: Option<f64>,
}

impl Default for RepulsiveOptions {
    /// The library's defaults (`SPEC.md` §10).
    fn default() -> Self {
        let lib = FlowOptions::default();
        let (armijo, shrink, max_halvings, gap) = match LineSearch::default() {
            LineSearch::CollisionGuarded {
                armijo,
                shrink,
                gap,
                max_halvings,
                ..
            } => (armijo, shrink, max_halvings, gap),
            LineSearch::NormalizedArmijo {
                armijo,
                shrink,
                max_halvings,
            } => (armijo, shrink, max_halvings, 0.1),
        };
        let (projection_tolerance, projection_absolute) = match lib.projection.tolerance {
            Tolerance::MeanEdgeLength(t) => (t, false),
            Tolerance::Absolute(t) => (t, true),
        };
        let (_stop_eps, gradient_norm) = match lib.stop {
            StopRule::GradientNorm { norm, eps } => (
                eps,
                match norm {
                    GradNorm::Euclidean => "euclidean",
                    GradNorm::MassL2 => "mass",
                },
            ),
            _ => (1e-4, "euclidean"),
        };
        let bh = BhOptions::default();
        RepulsiveOptions {
            alpha: lib.exponents.alpha(),
            beta: lib.exponents.beta(),
            metric: "sobolev".into(),
            quadrature: "midpoint".into(),
            length: LengthMode::Edges,
            grow_factor: LengthSchedule::default().factor,
            grow_rate: LengthSchedule::default().rate,
            max_subdivisions: lib.max_subdivisions,
            barycenter: None,
            pinned: Vec::new(),
            tangent_vertices: Vec::new(),
            tangent_edges: Vec::new(),
            surface: None,
            planar_axis: None,
            obstacles: Vec::new(),
            obstacle_weight: 0.0,
            obstacle_exponent: None,
            obstacle_theta: None,
            total_length_weight: 0.0,
            length_difference_weight: 0.0,
            field: None,
            field_weight: 1.0,
            line_search: "collision".into(),
            line_search_start: "two_thirds".into(),
            armijo,
            shrink,
            max_halvings,
            gap,
            backend: "auto".into(),
            bh_theta_x: bh.theta_x,
            bh_theta_t: bh.theta_t,
            bh_leaf_size: bh.leaf_size,
            threads: lib.threads,
            projection_tolerance,
            projection_absolute,
            projection_iterations: lib.projection.max_iterations,
            stop: "gradient".into(),
            stop_threshold: None,
            stall_window: 10,
            gradient_norm: gradient_norm.into(),
            max_steps: lib.max_iterations,
            time_limit: lib.time_limit.map(|d| d.as_secs_f64()),
        }
    }
}

fn one_of(name: &str, value: &str, allowed: &[&str]) -> Result<(), RepError> {
    if allowed.contains(&value) {
        Ok(())
    } else {
        Err(input(
            name,
            format!(
                "'{value}' is not valid; use {}",
                allowed
                    .iter()
                    .map(|a| format!("'{a}'"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        ))
    }
}

impl RepulsiveOptions {
    /// Validates the options and builds the library's [`FlowOptions`].
    ///
    /// # Errors
    /// [`RepError::Input`]: an option outside its domain (the message names it, or is the
    /// library's), an obstacle that is not a valid triangle soup.
    pub fn flow_options(&self) -> Result<FlowOptions, RepError> {
        let exponents = Exponents::new(self.alpha, self.beta)?;
        one_of("metric", &self.metric, &["sobolev", "l2", "h1", "h2"])?;
        one_of("quadrature", &self.quadrature, &["midpoint", "trapezoid"])?;
        let metric = match self.metric.as_str() {
            "l2" => Metric::L2,
            "h1" => Metric::H1,
            "h2" => Metric::H2,
            _ => Metric::Sobolev {
                quadrature: if self.quadrature == "trapezoid" {
                    GramQuadrature::Trapezoid
                } else {
                    GramQuadrature::Midpoint
                },
            },
        };

        // Constraints.
        let mut constraints = Vec::new();
        if self.barycenter.unwrap_or(self.pinned.is_empty()) {
            constraints.push(Constraint::Barycenter { target: None });
        }
        match self.length {
            LengthMode::Edges => constraints.push(Constraint::EdgeLengths {
                edges: None,
                targets: None,
            }),
            LengthMode::Total => constraints.push(Constraint::TotalLength { target: None }),
            LengthMode::Free => {}
        }
        for &vertex in &self.pinned {
            constraints.push(Constraint::Point {
                vertex,
                target: None,
            });
        }
        for &v in &self.tangent_vertices {
            constraints.push(Constraint::Tangent {
                site: TangentSite::Vertex(v),
                direction: None,
            });
        }
        for &e in &self.tangent_edges {
            constraints.push(Constraint::Tangent {
                site: TangentSite::Edge(e),
                direction: None,
            });
        }
        if let Some((spec, vertices)) = &self.surface {
            constraints.push(Constraint::Surface {
                vertices: vertices.clone(),
                surface: spec.build()?,
            });
        }

        // Penalties and obstacles.
        let mut penalties: Vec<Arc<dyn PenaltyTerm>> = Vec::new();
        if self.total_length_weight != 0.0 {
            penalties.push(Arc::new(TotalLengthPenalty::new(self.total_length_weight)?));
        }
        if self.length_difference_weight != 0.0 {
            penalties.push(Arc::new(LengthDifference::new(
                self.length_difference_weight,
            )?));
        }
        let mut obstacles: Vec<Arc<dyn rc::CollisionObstacle>> = Vec::new();
        for (vertices, triangles) in &self.obstacles {
            let obstacle = Arc::new(Obstacle::new(vertices.clone(), triangles.clone())?);
            if self.obstacle_weight != 0.0 {
                let p = self.obstacle_exponent.unwrap_or(self.beta - self.alpha);
                let theta = self
                    .obstacle_theta
                    .unwrap_or(rc::penalties::surface::DEFAULT_THETA);
                penalties.push(Arc::new(SurfacePotential::new(
                    obstacle,
                    self.obstacle_weight,
                    p,
                    theta,
                )?));
            } else {
                obstacles.push(obstacle);
            }
        }
        if let Some(field) = &self.field {
            penalties.push(Arc::new(rc::penalties::FieldAlignment::new(
                field.build()?,
                self.field_weight,
            )?));
        }

        // Line search.
        one_of("line_search", &self.line_search, &["collision", "armijo"])?;
        one_of(
            "line_search_start",
            &self.line_search_start,
            &["two_thirds", "grow_previous"],
        )?;
        let line_search = if self.line_search == "armijo" {
            LineSearch::NormalizedArmijo {
                armijo: self.armijo,
                shrink: self.shrink,
                max_halvings: self.max_halvings,
            }
        } else {
            LineSearch::CollisionGuarded {
                armijo: self.armijo,
                shrink: self.shrink,
                start: if self.line_search_start == "grow_previous" {
                    StartRule::grow_previous()
                } else {
                    StartRule::TwoThirdsTauMax
                },
                gap: self.gap,
                max_halvings: self.max_halvings,
            }
        };

        // Stop rule.
        one_of("stop", &self.stop, &["gradient", "relative", "stall"])?;
        one_of("gradient_norm", &self.gradient_norm, &["euclidean", "mass"])?;
        let threshold = |default: Option<f64>| {
            self.stop_threshold
                .or(default)
                .ok_or_else(|| input("stop_threshold", "required with stop='stall'"))
        };
        let stop = match self.stop.as_str() {
            "relative" => StopRule::RelativeGradient {
                ratio: threshold(Some(1e-3))?,
            },
            "stall" => StopRule::EnergyStall {
                relative: threshold(None)?,
                window: self.stall_window,
            },
            _ => StopRule::GradientNorm {
                norm: if self.gradient_norm == "mass" {
                    GradNorm::MassL2
                } else {
                    GradNorm::Euclidean
                },
                eps: threshold(Some(1e-4))?,
            },
        };

        // Back end.
        one_of(
            "backend",
            &self.backend,
            &["auto", "dense", "barnes_hut", "multigrid", "accelerated"],
        )?;
        let bh = EnergyEval::BarnesHut(BhOptions {
            theta_x: self.bh_theta_x,
            theta_t: self.bh_theta_t,
            leaf_size: self.bh_leaf_size,
        });
        let multigrid = || -> Result<Solver, RepError> {
            let pool = Arc::new(make_pool(self.threads)?);
            let dense = Arc::new(make_pool(Some(
                rc::DEFAULT_DENSE_THREADS.min(pool.current_num_threads()),
            ))?);
            Ok(Solver::Custom(Arc::new(MultigridSolvers::new(
                pool,
                dense,
                MgOptions::default(),
            )?)))
        };
        let backend = match self.backend.as_str() {
            "dense" => rc::Backend::Explicit {
                energy: EnergyEval::Exact,
                solver: Solver::Dense,
            },
            "barnes_hut" => rc::Backend::Explicit {
                energy: bh,
                solver: Solver::Dense,
            },
            "multigrid" => rc::Backend::Explicit {
                energy: EnergyEval::Exact,
                solver: multigrid()?,
            },
            "accelerated" => rc::Backend::Explicit {
                energy: bh,
                solver: multigrid()?,
            },
            _ => rc::Backend::default(),
        };

        let tolerance = if self.projection_absolute {
            Tolerance::Absolute(self.projection_tolerance)
        } else {
            Tolerance::MeanEdgeLength(self.projection_tolerance)
        };
        let time_limit = match self.time_limit {
            None => None,
            Some(s) if s.is_finite() && s >= 0.0 => Some(std::time::Duration::from_secs_f64(s)),
            Some(s) => {
                return Err(input(
                    "time_limit",
                    format!("{s} must be a finite number of seconds >= 0 (or None)"),
                ));
            }
        };
        let schedule = if self.grow_factor == 1.0 {
            None
        } else {
            Some(LengthSchedule {
                factor: self.grow_factor,
                rate: self.grow_rate,
            })
        };
        let o = FlowOptions {
            exponents,
            metric,
            backend,
            constraints,
            penalties,
            obstacles,
            planar: self.planar_axis,
            line_search,
            projection: ProjectionOptions {
                tolerance,
                max_iterations: self.projection_iterations,
            },
            stop,
            max_iterations: self.max_steps,
            time_limit,
            schedule,
            max_subdivisions: self.max_subdivisions,
            threads: self.threads,
            ..FlowOptions::default()
        };
        o.validate()?;
        Ok(o)
    }
}

// ---------------------------------------------------------------------------------------------
// Subdivision bookkeeping
// ---------------------------------------------------------------------------------------------

/// Where the vertices and edges of a subdivided network come from (see the module docs).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceMap {
    /// For each vertex, the input vertex it is, `-1` for a vertex created by subdivision.
    pub vertex_source: Vec<i32>,
    /// For each edge, the input edge it descends from.
    pub edge_source: Vec<u32>,
}

impl SourceMap {
    /// The identity map of a network with `n_vertices` vertices and `n_edges` edges.
    pub fn identity(n_vertices: usize, n_edges: usize) -> SourceMap {
        SourceMap {
            vertex_source: (0..n_vertices as i32).collect(),
            edge_source: (0..n_edges as u32).collect(),
        }
    }

    /// One subdivision: `|E|` midpoints are appended, and edge `I` becomes `I` and `|E| + I`.
    pub fn subdivide(&mut self) {
        let e = self.edge_source.len();
        self.vertex_source.extend(std::iter::repeat_n(-1, e));
        self.edge_source.extend_from_within(..);
    }
}

// ---------------------------------------------------------------------------------------------
// Status and reports
// ---------------------------------------------------------------------------------------------

/// The name of a session status.
pub fn status_name(s: &Status) -> &'static str {
    match s {
        Status::Running => "running",
        Status::Converged => "converged",
        Status::MaxIterations => "max_steps",
        Status::TimeLimit => "time_limit",
        Status::Stuck => "stuck",
        Status::Cancelled => "cancelled",
        Status::Error(_) => "error",
    }
}

fn step_status_name(s: StepStatus) -> &'static str {
    match s {
        StepStatus::Accepted => "accepted",
        StepStatus::Stuck => "stuck",
        StepStatus::Converged => "converged",
        StepStatus::RetargetStalled => "retarget_stalled",
    }
}

/// The report of one step as a dict. `status` is the outcome of the step (`"accepted"`,
/// `"stuck"`, `"converged"`, `"retarget_stalled"`); `session_status` that of the session after
/// it.
pub fn step_val(r: &StepReport, session_status: &Status) -> Val {
    let linear = Val::Dict(vec![
        ("solver", Val::Str(r.linear.solver.into())),
        ("dimension", Val::uint(r.linear.dimension)),
        ("iterations", Val::uint(r.linear.iterations)),
        ("relative_residual", Val::Float(r.linear.relative_residual)),
        ("fallback_used", Val::Bool(r.linear.fallback_used)),
        ("levels", Val::uint(r.linear.levels)),
        ("pivot_ratio", Val::opt_f(r.linear.pivot_ratio)),
        ("seconds", Val::Float(r.linear.seconds)),
    ]);
    let retarget = match &r.retarget {
        None => Val::None,
        Some(t) => Val::Dict(vec![
            (
                "status",
                Val::Str(
                    match t.status {
                        rc::RetargetStatus::Moved => "moved",
                        rc::RetargetStatus::Stalled => "stalled",
                    }
                    .into(),
                ),
            ),
            ("progress_from", Val::Float(t.progress_from)),
            ("progress_to", Val::Float(t.progress_to)),
            ("attempts", Val::Int(i64::from(t.attempts))),
            (
                "projection_iterations",
                Val::Int(i64::from(t.projection_iterations)),
            ),
            ("max_violation", Val::Float(t.max_violation)),
            ("collision", Val::Bool(t.collision)),
        ]),
    };
    let s = &r.seconds;
    let seconds = Val::Dict(vec![
        ("retarget", Val::Float(s.retarget)),
        ("energy", Val::Float(s.energy)),
        ("metric", Val::Float(s.metric)),
        ("solve", Val::Float(s.solve)),
        ("line_search", Val::Float(s.line_search)),
        ("projection", Val::Float(s.projection)),
        ("ccd", Val::Float(s.ccd)),
        ("subdivision", Val::Float(s.subdivision)),
        ("total", Val::Float(s.total)),
    ]);
    Val::Dict(vec![
        ("iteration", Val::uint(r.iteration)),
        ("status", Val::Str(step_status_name(r.status).into())),
        (
            "session_status",
            Val::Str(status_name(session_status).into()),
        ),
        ("n_vertices", Val::uint(r.n_vertices)),
        ("n_edges", Val::uint(r.n_edges)),
        ("energy_before", Val::Float(r.energy_before)),
        ("energy_after", Val::Float(r.energy_after)),
        ("tangent_point_after", Val::Float(r.tangent_point_after)),
        (
            "penalties_after",
            Val::List(r.penalties_after.iter().map(|&x| Val::Float(x)).collect()),
        ),
        ("gradient_norm", Val::Float(r.gradient_norm)),
        ("gradient_mass_norm", Val::Float(r.gradient_mass_norm)),
        ("sobolev_norm_sq", Val::Float(r.sobolev_norm_sq)),
        ("tau", Val::Float(r.tau)),
        ("tau_start", Val::Float(r.tau_start)),
        ("tau_max", Val::opt_f(r.tau_max)),
        ("halvings", Val::Int(i64::from(r.halvings))),
        ("trials", Val::uint(r.trials.len())),
        (
            "projection_iterations",
            Val::Int(i64::from(r.projection_iterations)),
        ),
        ("max_violation", Val::Float(r.max_violation)),
        (
            "dropped_rows",
            Val::List(r.dropped_rows.iter().map(|&i| Val::uint(i)).collect()),
        ),
        ("linear", linear),
        ("ccd_checks", Val::Int(i64::from(r.ccd_checks))),
        ("ccd_contact", Val::Bool(r.ccd.contact.is_some())),
        ("retarget", retarget),
        ("subdivided", Val::Bool(r.subdivided)),
        ("energy_backend", Val::Str(r.energy_backend.into())),
        ("solver_backend", Val::Str(r.solver_backend.into())),
        ("seconds", seconds),
    ])
}

/// The report `step` returns for a cancelled session (nothing was computed).
fn cancelled_val(iteration: usize, energy: f64, n_vertices: usize, n_edges: usize) -> Val {
    Val::Dict(vec![
        ("iteration", Val::uint(iteration)),
        ("status", Val::Str("cancelled".into())),
        ("session_status", Val::Str("cancelled".into())),
        ("n_vertices", Val::uint(n_vertices)),
        ("n_edges", Val::uint(n_edges)),
        ("energy_before", Val::Float(energy)),
        ("energy_after", Val::Float(energy)),
        ("subdivided", Val::Bool(false)),
    ])
}

// ---------------------------------------------------------------------------------------------
// The engine
// ---------------------------------------------------------------------------------------------

/// The last completed state of a session, readable without waiting for a step.
#[derive(Debug, Clone)]
pub struct Snapshot {
    /// Positions of the network.
    pub positions: Vec<[f64; 3]>,
    /// Edges of the network.
    pub edges: Vec<[u32; 2]>,
    /// Steps taken.
    pub iteration: usize,
    /// `F` (tangent-point energy plus penalties).
    pub energy: f64,
    /// Status name (see [`status_name`]).
    pub status: &'static str,
    /// The message of the error when the status is `"error"`.
    pub error: Option<String>,
    /// Where vertices and edges come from.
    pub sources: SourceMap,
    /// Number of subdivisions so far.
    pub subdivisions: u32,
    /// `ℓ̄₀`, the mean edge length of the start.
    pub mean_edge_length0: f64,
    /// The report of the last step.
    pub last_report: Option<Val>,
    /// `Σ ℓ_I` now.
    pub total_length: f64,
}

/// Why a [`Engine::run`] ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunEnd {
    /// `max_steps` steps were taken.
    MaxSteps,
    /// The session reached a final status (converged, stuck, step cap, ...).
    Final,
    /// The callback asked to stop.
    Callback,
    /// `cancel()` was called.
    Cancelled,
    /// The time limit of the options passed.
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
    network: CurveNetwork,
    positions: Vec<[f64; 3]>,
    options: FlowOptions,
}

fn poisoned() -> RepError {
    RepError::Solve("internal error: a previous call panicked inside this session".into())
}

fn make_snapshot(session: &Session, sources: SourceMap) -> Snapshot {
    let net = session.network();
    let positions = session.positions().to_vec();
    Snapshot {
        total_length: net.edge_geometry(&positions).length.iter().sum(),
        positions,
        edges: net.edges().to_vec(),
        iteration: session.iteration(),
        energy: session.energy(),
        status: status_name(&session.status()),
        error: None,
        sources,
        subdivisions: session.subdivisions(),
        mean_edge_length0: session.mean_edge_length0(),
        last_report: None,
    }
}

impl Engine {
    /// Starts a session on the curve network `edges` at `positions`.
    ///
    /// # Errors
    /// [`RepError::Input`]: the network, the positions or an option is invalid (the library's
    /// message).
    pub fn new(
        positions: Vec<[f64; 3]>,
        edges: Vec<[u32; 2]>,
        options: &RepulsiveOptions,
    ) -> Result<Engine, RepError> {
        let flow = options.flow_options()?;
        let network = CurveNetwork::new(positions.len(), edges)?;
        let session = Session::new(network.clone(), positions.clone(), flow.clone())?;
        let sources = SourceMap::identity(network.n_vertices(), network.n_edges());
        let snapshot = make_snapshot(&session, sources);
        let cancel = session.cancel_flag();
        Ok(Engine {
            session: Mutex::new(session),
            cancel: Mutex::new(cancel),
            snapshot: Mutex::new(snapshot),
            network,
            positions,
            options: flow,
        })
    }

    fn session(&self) -> Result<MutexGuard<'_, Session>, RepError> {
        self.session.lock().map_err(|_| poisoned())
    }

    fn snap(&self) -> Result<MutexGuard<'_, Snapshot>, RepError> {
        self.snapshot.lock().map_err(|_| poisoned())
    }

    fn cancel_flag(&self) -> Result<Arc<AtomicBool>, RepError> {
        Ok(self.cancel.lock().map_err(|_| poisoned())?.clone())
    }

    /// The last completed state. Never waits for a running step.
    ///
    /// # Errors
    /// [`RepError::Solve`] if a call panicked inside the session.
    pub fn snapshot(&self) -> Result<Snapshot, RepError> {
        Ok(self.snap()?.clone())
    }

    /// Asks the running step (between line-search trials) or the next one to stop. The snapshot
    /// reads `"cancelled"` at once. Never waits.
    ///
    /// # Errors
    /// [`RepError::Solve`] if a call panicked inside the session.
    pub fn cancel(&self) -> Result<(), RepError> {
        self.cancel_flag()?.store(true, Ordering::Relaxed);
        self.snap()?.status = "cancelled";
        Ok(())
    }

    /// Whether `cancel` was called since the start or the last `reset`.
    ///
    /// # Errors
    /// [`RepError::Solve`] if a call panicked inside the session.
    pub fn is_cancelled(&self) -> Result<bool, RepError> {
        Ok(self.cancel_flag()?.load(Ordering::Relaxed))
    }

    /// Starts again from the initial curve with the same options (the cancel flag is cleared).
    ///
    /// # Errors
    /// The errors of [`Engine::new`].
    pub fn reset(&self) -> Result<(), RepError> {
        let session = Session::new(
            self.network.clone(),
            self.positions.clone(),
            self.options.clone(),
        )?;
        let snapshot = make_snapshot(
            &session,
            SourceMap::identity(self.network.n_vertices(), self.network.n_edges()),
        );
        let mut s = self.session()?;
        *self.cancel.lock().map_err(|_| poisoned())? = session.cancel_flag();
        *s = session;
        *self.snap()? = snapshot;
        Ok(())
    }

    /// Changes the length schedule: the final length relative to the *current* targets and the
    /// progress per step. `factor = 1` removes it.
    ///
    /// # Errors
    /// [`RepError::Input`]: no length constraint, or a bad factor or rate.
    pub fn set_schedule(&self, factor: f64, rate: f64) -> Result<(), RepError> {
        let mut s = self.session()?;
        let schedule = (factor != 1.0).then_some(LengthSchedule { factor, rate });
        s.set_schedule(schedule)?;
        Ok(())
    }

    /// One step. A cancelled session does not step and returns a report with status
    /// `"cancelled"`.
    ///
    /// # Errors
    /// [`RepError::Solve`] / [`RepError::Input`] for a failed step; the session's status is
    /// then `"error"` and the curve that of the last completed step.
    pub fn step(&self) -> Result<Val, RepError> {
        let mut s = self.session()?;
        let finish_cancelled = |s: &Session| {
            let mut snap = self.snap()?;
            snap.status = "cancelled";
            Ok(cancelled_val(
                s.iteration(),
                s.energy(),
                s.network().n_vertices(),
                s.network().n_edges(),
            ))
        };
        if self.cancel_flag()?.load(Ordering::Relaxed) {
            return finish_cancelled(&s);
        }
        match s.step() {
            Ok(report) => {
                let status = s.status();
                let val = step_val(&report, &status);
                let mut snap = self.snap()?;
                if report.subdivided {
                    snap.sources.subdivide();
                    snap.edges = s.network().edges().to_vec();
                    debug_assert_eq!(snap.sources.vertex_source.len(), s.positions().len());
                }
                snap.positions.clear();
                snap.positions.extend_from_slice(s.positions());
                snap.total_length = s.network().edge_geometry(s.positions()).length.iter().sum();
                snap.iteration = s.iteration();
                snap.energy = s.energy();
                snap.status = status_name(&status);
                snap.subdivisions = s.subdivisions();
                snap.last_report = Some(val.clone());
                Ok(val)
            }
            Err(RcError::Cancelled) => finish_cancelled(&s),
            Err(e) => {
                let status = s.status();
                let mut snap = self.snap()?;
                snap.status = status_name(&status);
                snap.error = Some(e.to_string());
                Err(e.into())
            }
        }
    }

    /// Steps until the session reaches a final status, `max_steps` steps of this run are taken,
    /// `cancel` is called, `time_limit` seconds have passed, or `progress` returns `false`. The
    /// callback gets the report of every step. A session that is final already does not step.
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
    ) -> Result<RunSummary, RepError> {
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
            if self.snap()?.status != "running" {
                break RunEnd::Final;
            }
            let report = self.step()?;
            let cancelled = report.get("status") == Some(&Val::Str("cancelled".into()));
            if cancelled {
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

// ---------------------------------------------------------------------------------------------
// One-shot flow and energy
// ---------------------------------------------------------------------------------------------

/// The result of [`flow_once`].
#[derive(Debug, Clone)]
pub struct FlowOutput {
    /// Final positions.
    pub positions: Vec<[f64; 3]>,
    /// Final edges.
    pub edges: Vec<[u32; 2]>,
    /// The report (see the glue for the keys).
    pub report: Val,
}

/// `ddg_repulsive_curves::flow` with reports as [`Val`]s. `progress` gets the report of every
/// step and returns `false` to stop (the status is then `"cancelled"`, as in the library).
///
/// # Errors
/// As [`Engine::new`].
pub fn flow_once(
    positions: Vec<[f64; 3]>,
    edges: Vec<[u32; 2]>,
    options: &RepulsiveOptions,
    mut progress: Option<&mut dyn FnMut(&Val) -> bool>,
) -> Result<FlowOutput, RepError> {
    let flow = options.flow_options()?;
    let network = CurveNetwork::new(positions.len(), edges)?;
    let mut sources = SourceMap::identity(network.n_vertices(), network.n_edges());
    let mut subdivisions = 0u32;
    let mut cb = |r: &StepReport| {
        if r.subdivided {
            sources.subdivide();
            subdivisions += 1;
        }
        let go = progress
            .as_mut()
            .is_none_or(|p| p(&step_val(r, &Status::Running)));
        if go { Control::Continue } else { Control::Stop }
    };
    let result = rc::flow(network, positions, flow, Some(&mut cb))?;
    let r = &result.report;
    let status = match &r.status {
        Status::Error(e) => return Err(e.clone().into()),
        s => s,
    };
    let report = Val::Dict(vec![
        ("status", Val::Str(status_name(status).into())),
        ("iterations", Val::uint(r.iterations)),
        (
            "energy_initial",
            r.steps
                .first()
                .map_or(Val::None, |s| Val::Float(s.energy_before)),
        ),
        ("energy_final", Val::Float(r.final_energy)),
        ("seconds", Val::Float(r.seconds)),
        ("n_subdivisions", Val::Int(i64::from(subdivisions))),
        (
            "vertex_source",
            Val::List(
                sources
                    .vertex_source
                    .iter()
                    .map(|&v| Val::Int(i64::from(v)))
                    .collect(),
            ),
        ),
        (
            "edge_source",
            Val::List(
                sources
                    .edge_source
                    .iter()
                    .map(|&e| Val::Int(i64::from(e)))
                    .collect(),
            ),
        ),
        (
            "steps",
            Val::List(r.steps.iter().map(|s| step_val(s, status)).collect()),
        ),
    ]);
    Ok(FlowOutput {
        edges: result.network.edges().to_vec(),
        positions: result.positions,
        report,
    })
}

/// The tangent-point energy and its differential (one array of 3 per vertex).
///
/// # Errors
/// [`RepError::Input`] for a bad network, positions or option.
pub fn curve_energy(
    positions: &[[f64; 3]],
    edges: Vec<[u32; 2]>,
    alpha: f64,
    beta: f64,
    barnes_hut: Option<BhOptions>,
    threads: Option<usize>,
) -> Result<(f64, Vec<[f64; 3]>), RepError> {
    let network = CurveNetwork::new(positions.len(), edges)?;
    let eval = barnes_hut.map_or(EnergyEval::Exact, EnergyEval::BarnesHut);
    let report = rc::energy(
        &network,
        positions,
        Exponents::new(alpha, beta)?,
        &eval,
        threads,
    )?;
    Ok((report.energy, report.differential))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A trefoil with `n` vertices.
    pub(crate) fn trefoil(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 2]>) {
        let positions = (0..n)
            .map(|i| {
                let t = std::f64::consts::TAU * i as f64 / n as f64;
                [
                    t.sin() + 2.0 * (2.0 * t).sin(),
                    t.cos() - 2.0 * (2.0 * t).cos(),
                    -(3.0 * t).sin(),
                ]
            })
            .collect();
        let edges = (0..n as u32).map(|i| [i, (i + 1) % n as u32]).collect();
        (positions, edges)
    }

    fn debug(o: &FlowOptions) -> String {
        format!("{o:?}")
    }

    #[test]
    fn defaults_equal_the_library() {
        let o = RepulsiveOptions::default().flow_options().unwrap();
        let lib = FlowOptions::default();
        assert_eq!(debug(&o), debug(&lib));
        let d = RepulsiveOptions::default();
        assert_eq!((d.alpha, d.beta), (3.0, 6.0));
        assert_eq!(d.max_steps, 1000);
        assert_eq!(d.max_subdivisions, 0);
        assert_eq!(d.stop_threshold, None);
        assert_eq!(d.grow_rate, 0.01);
        assert_eq!(d.grow_factor, 1.0);
    }

    #[test]
    fn bad_options_are_input_errors() {
        let bad: Vec<(&str, RepulsiveOptions)> = vec![
            (
                "alpha",
                RepulsiveOptions {
                    alpha: 0.5,
                    ..Default::default()
                },
            ),
            (
                "beta",
                RepulsiveOptions {
                    beta: f64::NAN,
                    ..Default::default()
                },
            ),
            (
                "line_search",
                RepulsiveOptions {
                    line_search: "magic".into(),
                    ..Default::default()
                },
            ),
            (
                "armijo",
                RepulsiveOptions {
                    armijo: 2.0,
                    ..Default::default()
                },
            ),
            (
                "backend",
                RepulsiveOptions {
                    backend: "gpu".into(),
                    ..Default::default()
                },
            ),
            (
                "threads",
                RepulsiveOptions {
                    threads: Some(0),
                    ..Default::default()
                },
            ),
            (
                "stall",
                RepulsiveOptions {
                    stop: "stall".into(),
                    stop_threshold: None,
                    ..Default::default()
                },
            ),
            (
                "stop_threshold",
                RepulsiveOptions {
                    stop_threshold: Some(-1.0),
                    ..Default::default()
                },
            ),
            (
                "time_limit",
                RepulsiveOptions {
                    time_limit: Some(-1.0),
                    ..Default::default()
                },
            ),
            (
                "weight",
                RepulsiveOptions {
                    total_length_weight: f64::NAN,
                    ..Default::default()
                },
            ),
            (
                "obstacle",
                RepulsiveOptions {
                    obstacles: vec![(vec![[0.0; 3]; 3], vec![[0, 1, 5]])],
                    ..Default::default()
                },
            ),
            (
                "surface",
                RepulsiveOptions {
                    surface: Some((
                        SurfaceSpec::Sphere {
                            center: [0.0; 3],
                            radius: -1.0,
                        },
                        None,
                    )),
                    ..Default::default()
                },
            ),
        ];
        for (what, o) in bad {
            match o.flow_options() {
                Err(RepError::Input(m)) => assert!(!m.is_empty(), "{what}"),
                other => panic!(
                    "{what}: expected an input error, got {:?}",
                    other.map(|_| ())
                ),
            }
        }
        // A schedule without a length constraint is the library's error.
        let o = RepulsiveOptions {
            length: LengthMode::Free,
            grow_factor: 2.0,
            ..Default::default()
        };
        let (p, e) = trefoil(40);
        assert!(matches!(Engine::new(p, e, &o), Err(RepError::Input(_))));
    }

    #[test]
    fn invalid_network_is_an_input_error() {
        let (p, _) = trefoil(10);
        let o = RepulsiveOptions::default();
        for edges in [vec![[0, 99]], vec![[3, 3]], vec![[0, 1], [1, 0]]] {
            assert!(
                matches!(Engine::new(p.clone(), edges, &o), Err(RepError::Input(_))),
                "bad edges accepted"
            );
        }
    }

    #[test]
    fn session_matches_flow_bitwise() {
        let (p, e) = trefoil(48);
        let o = RepulsiveOptions {
            max_steps: 6,
            threads: Some(2),
            ..Default::default()
        };
        let engine = Engine::new(p.clone(), e.clone(), &o).unwrap();
        for _ in 0..6 {
            engine.step().unwrap();
        }
        let flow = flow_once(p, e, &o, None).unwrap();
        let snap = engine.snapshot().unwrap();
        assert_eq!(snap.positions, flow.positions);
        assert_eq!(snap.status, "max_steps");
    }

    #[test]
    fn cancel_and_reset() {
        let (p, e) = trefoil(40);
        let engine = Engine::new(p.clone(), e, &RepulsiveOptions::default()).unwrap();
        engine.step().unwrap();
        engine.cancel().unwrap();
        let r = engine.step().unwrap();
        assert_eq!(r.get("status"), Some(&Val::Str("cancelled".into())));
        assert_eq!(engine.snapshot().unwrap().status, "cancelled");
        engine.reset().unwrap();
        let s = engine.snapshot().unwrap();
        assert_eq!(
            (s.status, s.iteration, s.positions.clone()),
            ("running", 0, p)
        );
        engine.step().unwrap();
    }

    #[test]
    fn source_map_follows_subdivision() {
        let mut m = SourceMap::identity(3, 3);
        m.subdivide();
        assert_eq!(m.vertex_source, vec![0, 1, 2, -1, -1, -1]);
        assert_eq!(m.edge_source, vec![0, 1, 2, 0, 1, 2]);
        m.subdivide();
        assert_eq!(m.vertex_source.len(), 12);
        assert_eq!(m.edge_source.len(), 12);
        assert_eq!(&m.edge_source[6..], &[0, 1, 2, 0, 1, 2]);
    }
}
