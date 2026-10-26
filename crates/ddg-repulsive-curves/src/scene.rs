//! Scene description for the applications of `SPEC.md` §14.4 (work package WP13).
//!
//! A [`SceneSpec`] says what to flow (which constraints, penalties, schedule, back end); a
//! [`Scene`] is a spec bound to a curve network and its positions. Two ways to get one:
//!
//! * [`parse_scene`] / [`read_scene`]: an importer for the scene files of the authors' program
//!   (`scenes/FORMATS.md` of `repulsive-curves`, MIT, Copyright (c) 2019 Christopher Yu; the
//!   licence notice is in `LICENSES/`).
//! * the generators of [`builtin`]: small procedural versions of each application, used by the
//!   tests and as fast demonstrations (`examples/scenes.rs`).
//!
//! [`run_scene`] runs a scene to convergence or a step budget, optionally writing OBJ frames, and
//! returns a [`RunReport`] (steps, time, energy, minimal distance between non-adjacent edges,
//! constraint residuals) with a JSON writer; [`audit`] is the geometric audit it uses.
//!
//! # Directives of the reference format
//!
//! Supported: `curve`, `repel_curve [alpha beta [weight]]`, `repel_surface file [weight]`,
//! `repel_plane cx cy cz nx ny nz [weight]` (a large square obstacle, not an infinite plane),
//! `optimize_length [w]`, `optimize_length_diffs [w]`, `optimize_field circular|interesting [w]`,
//! `fix_barycenter`, `fix_length [s]`, `fix_totallength [s]`, `fix_edgelengths [s]`,
//! `fix_vertex i`, `fix_tangent i`, `fix_special_vertices`, `fix_endpoint_vertices`,
//! `fix_special_tangents`, `constraint_surface sphere|torus|yplane|doubletorus [params]`,
//! `union_type disjoint|smooth`, `constrain_vertex i`, `constrain_all`, `constrain_endpoints`,
//! `iteration_limit n`, `subdivide_limit n`, comment lines (`#`).
//!
//! Not supported (recorded in [`Scene::notes`], never silently): `optimize_pin_angles`,
//! `optimize_area` (no such penalty in the paper's Sec. 8.1 list that this crate implements),
//! `show_surface` (display only). Differences from the reference, also noted in the scene:
//!
//! * The weight of `repel_curve` scales the energy; since the flow's energy has weight one, the
//!   other penalty weights are divided by it (same minimizers).
//! * `fix_edgelengths s` / `fix_totallength s` / `fix_length s` set a [`LengthSchedule`] with
//!   `factor = s` and `rate = 0.01` (`SPEC.md` §5.4): one common factor for all edges.
//! * A curve constrained to the plane `y = 0` cannot also have a barycenter constraint (its `y`
//!   row would be dependent on the surface rows, `SPEC.md` §5.2): the importer pins vertex 0
//!   instead and drops that vertex's surface row.
//! * `fix_special_tangents` adds tangents at endpoints (degree 1) only; at junctions
//!   (degree 3 or more) a vertex tangent is undefined (`InvalidTangentSite`).
//! * No constraint directive at all means `default_constraints()` (barycenter and all edge
//!   lengths), as in the reference.
//! * The reference stops at an absolute Sobolev dot product; the stop rule here is the crate's
//!   default (`SPEC.md` §10), usually not reached before `iteration_limit`.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use rayon::prelude::*;

use crate::ccd::segment_segment_distance;
use crate::constraints::{Constraint, ImplicitSurface, TangentSite, default_constraints};
use crate::energy::{BhOptions, EnergyEval};
use crate::error::RcError;
use crate::io::{parse_curve_obj, write_curve_obj};
use crate::kernel::Exponents;
use crate::multigrid::MgOptions;
use crate::network::CurveNetwork;
use crate::obstacle::Obstacle;
use crate::penalties::{FieldAlignment, LengthDifference, SurfacePotential, TotalLengthPenalty};
use crate::step::schedule::LengthSchedule;
use crate::step::stop::Status;
use crate::step::{Backend, FlowOptions, PenaltyTerm, Session, Solver};
use crate::surfaces::{
    CircularField, DoubleTorus, Plane, SmoothUnion, Sphere, SphereTangentField, Torus, Union,
    VectorField,
};

type P = [f64; 3];

// =============================================================================================
// Spec
// =============================================================================================

/// How the edge lengths (or the total length) are constrained.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LengthMode {
    /// Not constrained.
    Free,
    /// Every edge keeps its length; `scale` is the final factor of a growing schedule.
    EdgeLengths {
        /// Final factor relative to the start (1: constant).
        scale: f64,
    },
    /// The total length is kept; `scale` as above.
    TotalLength {
        /// Final factor relative to the start (1: constant).
        scale: f64,
    },
}

/// Which vertices are constrained to the surface.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum SurfaceVertices {
    /// None.
    #[default]
    None,
    /// All vertices (`constrain_all`).
    All,
    /// The endpoints, vertices of degree 1 (`constrain_endpoints`).
    Endpoints,
    /// The listed vertices (`constrain_vertex`).
    List(Vec<u32>),
}

/// A repulsive obstacle: a triangle mesh and its weight.
#[derive(Debug, Clone)]
pub struct ObstacleTerm {
    /// The mesh.
    pub mesh: Arc<Obstacle>,
    /// The weight of the surface potential.
    pub weight: f64,
    /// Where it came from (a file name), for reports.
    pub source: String,
}

/// What to do with a curve network (the reference's `SceneData`).
#[derive(Clone)]
pub struct SceneSpec {
    /// `(α, β)`, default (3, 6).
    pub exponents: Exponents,
    /// Weight of the tangent-point energy (the others are divided by it).
    pub repel_weight: f64,
    /// Repulsive meshes.
    pub obstacles: Vec<ObstacleTerm>,
    /// The length constraint.
    pub length: LengthMode,
    /// Barycenter constraint.
    pub barycenter: bool,
    /// Pinned vertices (to their start positions).
    pub pins: Vec<u32>,
    /// Pin all vertices of degree other than 2.
    pub pin_special: bool,
    /// Pin all vertices of degree 1.
    pub pin_endpoints: bool,
    /// Tangent constraints at these vertices (degree 1 or 2).
    pub tangent_pins: Vec<u32>,
    /// Tangent constraints at the endpoints of degree 1 (of the special vertices).
    pub special_tangents: bool,
    /// The constraint surface and its name.
    pub surface: Option<(Arc<dyn ImplicitSurface>, String)>,
    /// Which vertices are constrained to it.
    pub surface_vertices: SurfaceVertices,
    /// Weight of the total-length penalty (0: none).
    pub length_penalty: f64,
    /// Weight of the length-difference penalty (0: none).
    pub length_diff_penalty: f64,
    /// Field alignment penalty: field name and weight.
    pub field: Option<(String, f64)>,
    /// Stop after this many steps (`None`: the crate's default cap).
    pub iteration_limit: Option<usize>,
    /// Number of subdivisions allowed.
    pub subdivide_limit: u32,
    /// Growth of the target length per step as a fraction of the start length.
    pub schedule_rate: f64,
    /// Constraint list given explicitly (overrides all of the above constraint fields).
    pub explicit_constraints: Option<Vec<Constraint>>,
    /// Extra penalty terms (built-in scenes).
    pub extra_penalties: Vec<Arc<dyn PenaltyTerm>>,
}

impl std::fmt::Debug for SceneSpec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SceneSpec")
            .field("exponents", &self.exponents)
            .field("length", &self.length)
            .field("barycenter", &self.barycenter)
            .field("pins", &self.pins.len())
            .field("surface", &self.surface.as_ref().map(|s| &s.1))
            .field("obstacles", &self.obstacles.len())
            .finish_non_exhaustive()
    }
}

impl Default for SceneSpec {
    fn default() -> Self {
        SceneSpec {
            exponents: Exponents::default(),
            repel_weight: 1.0,
            obstacles: Vec::new(),
            length: LengthMode::Free,
            barycenter: false,
            pins: Vec::new(),
            pin_special: false,
            pin_endpoints: false,
            tangent_pins: Vec::new(),
            special_tangents: false,
            surface: None,
            surface_vertices: SurfaceVertices::None,
            length_penalty: 0.0,
            length_diff_penalty: 0.0,
            field: None,
            iteration_limit: None,
            subdivide_limit: 0,
            schedule_rate: 0.01,
            explicit_constraints: None,
            extra_penalties: Vec::new(),
        }
    }
}

/// A scene: a curve, what to do with it, and bookkeeping for the audit.
#[derive(Clone)]
pub struct Scene {
    /// A name for reports and file names.
    pub name: String,
    /// The curve network.
    pub network: CurveNetwork,
    /// Its start positions.
    pub positions: Vec<P>,
    /// The spec.
    pub spec: SceneSpec,
    /// What the importer did or did not do (see the module documentation).
    pub notes: Vec<String>,
    /// Obstacle meshes with their file names (for renderers).
    pub obstacle_sources: Vec<String>,
    /// The constraints that [`Scene::options`] builds, resolved against the network
    /// (pins carry their start positions).
    pub resolved: ResolvedScene,
}

/// The constraints of a scene in terms the audit understands.
#[derive(Clone, Default)]
pub struct ResolvedScene {
    /// Pinned vertices and their positions.
    pub pins: Vec<(u32, P)>,
    /// Vertices with a tangent constraint and the direction they must keep.
    pub tangents: Vec<(u32, P)>,
    /// Vertices constrained to the surface.
    pub surface_vertices: Vec<u32>,
    /// The constraint surface.
    pub surface: Option<Arc<dyn ImplicitSurface>>,
    /// The final length factor of the schedule (1: none).
    pub length_factor: f64,
    /// Total length of the start curve.
    pub length0: f64,
}

/// The back end of a run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackendChoice {
    /// The session's [`Backend::Auto`]: exact energy and dense solver up to the measured
    /// crossover, Barnes–Hut and multigrid above (`SPEC.md` §17h); re-chosen after a subdivision.
    Auto,
    /// Exact energy, dense solver.
    Exact,
    /// Barnes–Hut energy, dense solver.
    BarnesHutDense,
    /// Barnes–Hut energy, multigrid solver.
    BarnesHutMultigrid,
    /// Exact energy, multigrid solver.
    ExactMultigrid,
}

impl BackendChoice {
    /// Parses `auto`, `exact`, `bh` (Barnes–Hut + dense), `mg` or `bh-mg` (Barnes–Hut +
    /// multigrid), `emg` or `exact-mg` (exact + multigrid).
    pub fn parse(s: &str) -> Option<BackendChoice> {
        Some(match s {
            "auto" => BackendChoice::Auto,
            "exact" => BackendChoice::Exact,
            "bh" => BackendChoice::BarnesHutDense,
            "mg" | "bh-mg" => BackendChoice::BarnesHutMultigrid,
            "emg" | "exact-mg" => BackendChoice::ExactMultigrid,
            _ => return None,
        })
    }
}

fn unit(v: P) -> P {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dist(a: P, b: P) -> f64 {
    let d = sub(a, b);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// The tangent at a vertex of degree 1 or 2 as the constraint defines it (`SPEC.md` §5.1).
fn vertex_tangent(network: &CurveNetwork, positions: &[P], v: u32) -> Option<P> {
    let inc = network.incident_edges(v);
    let other = |e: u32| {
        let [a, b] = network.edges()[e as usize];
        if a == v { b } else { a }
    };
    match inc.len() {
        1 => Some(unit(sub(
            positions[other(inc[0]) as usize],
            positions[v as usize],
        ))),
        2 => {
            let (a, b) = (other(inc[0]) as usize, other(inc[1]) as usize);
            let x = positions[v as usize];
            let s1 = unit(sub(x, positions[a]));
            let s2 = unit(sub(positions[b], x));
            Some(unit([s1[0] + s2[0], s1[1] + s2[1], s1[2] + s2[2]]))
        }
        _ => None,
    }
}

impl Scene {
    /// Binds a spec to a network. Resolves the pins, the tangents and the surface vertices.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for an index out of range or a tangent at a vertex of degree
    /// three or more; the errors of the penalty constructors.
    pub fn new(
        name: impl Into<String>,
        network: CurveNetwork,
        positions: Vec<P>,
        spec: SceneSpec,
    ) -> Result<Scene, RcError> {
        let mut scene = Scene {
            name: name.into(),
            obstacle_sources: spec.obstacles.iter().map(|o| o.source.clone()).collect(),
            network,
            positions,
            spec,
            notes: Vec::new(),
            resolved: ResolvedScene::default(),
        };
        scene.resolve()?;
        Ok(scene)
    }

    fn resolve(&mut self) -> Result<(), RcError> {
        let nv = self.network.n_vertices();
        let bad = |what: &str, i: u32| RcError::InvalidOption {
            name: "scene",
            value: format!("{what} {i}"),
            requirement: "a vertex index below the number of vertices",
        };
        let sp = &self.spec;
        let mut pins: Vec<u32> = sp.pins.clone();
        for v in 0..nv as u32 {
            let d = self.network.degree(v);
            if (sp.pin_special && d != 2) || (sp.pin_endpoints && d == 1) {
                pins.push(v);
            }
        }
        pins.sort_unstable();
        pins.dedup();
        let mut tangents: Vec<u32> = sp.tangent_pins.clone();
        if sp.special_tangents {
            for v in 0..nv as u32 {
                let d = self.network.degree(v);
                if d == 1 && (sp.pin_special || sp.pin_endpoints || pins.contains(&v)) {
                    tangents.push(v);
                } else if d >= 3 && (sp.pin_special || sp.pin_endpoints) {
                    // junctions have no vertex tangent
                }
            }
        }
        tangents.sort_unstable();
        tangents.dedup();
        for &v in pins.iter().chain(tangents.iter()) {
            if v as usize >= nv {
                return Err(bad("vertex", v));
            }
        }
        let mut tan = Vec::new();
        for &v in &tangents {
            match vertex_tangent(&self.network, &self.positions, v) {
                Some(t) => tan.push((v, t)),
                None => {
                    return Err(RcError::InvalidTangentSite {
                        vertex: v,
                        degree: self.network.degree(v),
                    });
                }
            }
        }
        let mut svs: Vec<u32> = match &sp.surface_vertices {
            SurfaceVertices::None => Vec::new(),
            SurfaceVertices::All => (0..nv as u32).collect(),
            SurfaceVertices::Endpoints => (0..nv as u32)
                .filter(|&v| self.network.degree(v) == 1)
                .collect(),
            SurfaceVertices::List(l) => l.clone(),
        };
        for &v in &svs {
            if v as usize >= nv {
                return Err(bad("surface vertex", v));
            }
        }
        // A pinned point on the surface carries a redundant surface row: drop the row.
        svs.retain(|v| !pins.contains(v));
        let scale = match sp.length {
            LengthMode::Free => 1.0,
            LengthMode::EdgeLengths { scale } | LengthMode::TotalLength { scale } => scale,
        };
        self.resolved = ResolvedScene {
            pins: pins
                .iter()
                .map(|&v| (v, self.positions[v as usize]))
                .collect(),
            tangents: tan,
            surface_vertices: svs,
            surface: sp.surface.as_ref().map(|s| s.0.clone()),
            length_factor: scale,
            length0: self
                .network
                .edge_geometry(&self.positions)
                .length
                .iter()
                .sum(),
        };
        Ok(())
    }

    /// The flow options of the scene for `backend` and `threads`. Returns the options and the
    /// name of the back end used.
    ///
    /// # Errors
    /// Errors of the penalty and solver constructors.
    pub fn options(
        &self,
        backend: BackendChoice,
        threads: Option<usize>,
    ) -> Result<FlowOptions, RcError> {
        let sp = &self.spec;
        let mut o = FlowOptions {
            exponents: sp.exponents,
            threads,
            max_subdivisions: sp.subdivide_limit,
            ..FlowOptions::default()
        };
        if let Some(n) = sp.iteration_limit {
            o.max_iterations = n;
        }
        // constraints
        if let Some(c) = &sp.explicit_constraints {
            o.constraints = c.clone();
        } else {
            o.constraints = self.constraints();
        }
        // schedule
        let scale = self.resolved.length_factor;
        if (scale - 1.0).abs() > 1e-12 {
            o.schedule = Some(LengthSchedule {
                factor: scale,
                rate: sp.schedule_rate,
            });
        }
        // penalties
        let w = sp.repel_weight;
        for ob in &sp.obstacles {
            o.penalties.push(Arc::new(SurfacePotential::with_exponents(
                ob.mesh.clone(),
                ob.weight / w,
                &sp.exponents,
            )?));
        }
        if sp.length_penalty != 0.0 {
            o.penalties
                .push(Arc::new(TotalLengthPenalty::new(sp.length_penalty / w)?));
        }
        if sp.length_diff_penalty != 0.0 {
            o.penalties
                .push(Arc::new(LengthDifference::new(sp.length_diff_penalty / w)?));
        }
        if let Some((name, fw)) = &sp.field {
            let field: Arc<dyn VectorField> = match name.as_str() {
                "circular" => Arc::new(CircularField),
                "interesting" => Arc::new(SphereTangentField),
                other => {
                    return Err(RcError::InvalidOption {
                        name: "optimize_field",
                        value: other.to_string(),
                        requirement: "circular or interesting",
                    });
                }
            };
            o.penalties
                .push(Arc::new(FieldAlignment::new(field, fw / w)?));
        }
        o.penalties.extend(sp.extra_penalties.iter().cloned());
        // back end
        o.backend = match backend {
            BackendChoice::Auto => Backend::default(),
            BackendChoice::Exact => Backend::Explicit {
                energy: EnergyEval::Exact,
                solver: Solver::Dense,
            },
            BackendChoice::BarnesHutDense => Backend::Explicit {
                energy: EnergyEval::BarnesHut(BhOptions::default()),
                solver: Solver::Dense,
            },
            BackendChoice::BarnesHutMultigrid => Backend::Explicit {
                energy: EnergyEval::BarnesHut(BhOptions::default()),
                solver: Solver::Multigrid(MgOptions::default()),
            },
            BackendChoice::ExactMultigrid => Backend::Explicit {
                energy: EnergyEval::Exact,
                solver: Solver::Multigrid(MgOptions::default()),
            },
        };
        Ok(o)
    }

    /// The constraint list of the spec (the reference's order: barycenter, length, pins,
    /// tangents, surface).
    pub fn constraints(&self) -> Vec<Constraint> {
        let sp = &self.spec;
        let r = &self.resolved;
        let none = !sp.barycenter
            && matches!(sp.length, LengthMode::Free)
            && r.pins.is_empty()
            && r.tangents.is_empty()
            && r.surface_vertices.is_empty();
        if none {
            return default_constraints();
        }
        let mut c = Vec::new();
        let mut pins = r.pins.clone();
        let mut barycenter = sp.barycenter;
        if barycenter && sp.surface.as_ref().is_some_and(|s| s.1 == "yplane") {
            if matches!(sp.surface_vertices, SurfaceVertices::All) && pins.is_empty() {
                // the barycenter's y row is dependent on the surface rows (SPEC §5.2)
                pins.push((0, self.positions[0]));
            }
            if !pins.is_empty() {
                barycenter = false;
            }
        }
        if barycenter {
            c.push(Constraint::Barycenter { target: None });
        }
        match sp.length {
            LengthMode::Free => {}
            LengthMode::EdgeLengths { .. } => c.push(Constraint::EdgeLengths {
                edges: None,
                targets: None,
            }),
            LengthMode::TotalLength { .. } => c.push(Constraint::TotalLength { target: None }),
        }
        for (v, p) in &pins {
            c.push(Constraint::Point {
                vertex: *v,
                target: Some(*p),
            });
        }
        for (v, t) in &r.tangents {
            c.push(Constraint::Tangent {
                site: TangentSite::Vertex(*v),
                direction: Some(*t),
            });
        }
        let mut svs = r.surface_vertices.clone();
        svs.retain(|v| !pins.iter().any(|(p, _)| p == v));
        if let (Some((f, _)), false) = (&sp.surface, svs.is_empty()) {
            c.push(Constraint::Surface {
                vertices: Some(svs),
                surface: f.clone(),
            });
        }
        c
    }

    /// Starts a session.
    ///
    /// # Errors
    /// Errors of [`Scene::options`] and [`Session::new`].
    pub fn session(
        &self,
        backend: BackendChoice,
        threads: Option<usize>,
    ) -> Result<Session, RcError> {
        Session::new(
            self.network.clone(),
            self.positions.clone(),
            self.options(backend, threads)?,
        )
    }
}

// =============================================================================================
// Importer
// =============================================================================================

fn perr(line: usize, message: impl Into<String>) -> RcError {
    RcError::Parse {
        line,
        message: message.into(),
    }
}

fn num(parts: &[&str], i: usize, line: usize, what: &str) -> Result<f64, RcError> {
    parts
        .get(i)
        .ok_or_else(|| perr(line, format!("{what}: missing argument {i}")))?
        .parse::<f64>()
        .map_err(|_| perr(line, format!("{what}: `{}` is not a number", parts[i])))
}

fn idx(parts: &[&str], i: usize, line: usize, what: &str) -> Result<u32, RcError> {
    parts
        .get(i)
        .ok_or_else(|| perr(line, format!("{what}: missing argument")))?
        .parse::<u32>()
        .map_err(|_| {
            perr(
                line,
                format!("{what}: `{}` is not a vertex index", parts[i]),
            )
        })
}

/// A large square obstacle standing for the plane `repel_plane` (the reference uses the
/// infinite plane). Side `2 r` around the foot of the origin.
fn plane_obstacle(center: P, normal: P, r: f64) -> Result<Obstacle, RcError> {
    let n = unit(normal);
    let a = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = unit([
        a[1] * n[2] - a[2] * n[1],
        a[2] * n[0] - a[0] * n[2],
        a[0] * n[1] - a[1] * n[0],
    ]);
    let v = [
        n[1] * u[2] - n[2] * u[1],
        n[2] * u[0] - n[0] * u[2],
        n[0] * u[1] - n[1] * u[0],
    ];
    let k = 8usize;
    let mut verts = Vec::new();
    for i in 0..=k {
        for j in 0..=k {
            let (s, t) = (
                r * (2.0 * i as f64 / k as f64 - 1.0),
                r * (2.0 * j as f64 / k as f64 - 1.0),
            );
            verts.push([
                center[0] + s * u[0] + t * v[0],
                center[1] + s * u[1] + t * v[1],
                center[2] + s * u[2] + t * v[2],
            ]);
        }
    }
    let mut tris = Vec::new();
    for i in 0..k {
        for j in 0..k {
            let a = (i * (k + 1) + j) as u32;
            let b = a + 1;
            let c = a + (k + 1) as u32;
            let d = c + 1;
            tris.push([a, b, d]);
            tris.push([a, d, c]);
        }
    }
    Obstacle::new(verts, tris)
}

/// Parses the text of a scene file of the reference format. `base` is the directory the paths
/// are relative to; `name` names the scene.
///
/// # Errors
/// [`RcError::Parse`] with the line of a malformed directive or an unknown surface type;
/// [`RcError::Io`] for a missing curve or obstacle file; the errors of the curve reader.
pub fn parse_scene(text: &str, base: &Path, name: &str) -> Result<Scene, RcError> {
    let mut spec = SceneSpec::default();
    let mut curve: Option<PathBuf> = None;
    let mut notes = Vec::new();
    let mut smooth_union = false;
    let mut any_constraint = false;
    let mut planes: Vec<(P, P, f64)> = Vec::new();
    let mut surface: Option<Arc<dyn ImplicitSurface>> = None;
    let mut surface_names: Vec<String> = Vec::new();
    for (ln, raw) in text.lines().enumerate() {
        let line = ln + 1;
        let t = raw.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = t.split_whitespace().collect();
        let key = parts[0];
        let n = parts.len();
        match key {
            "curve" => {
                if n != 2 {
                    return Err(perr(line, "curve takes one path"));
                }
                curve = Some(base.join(parts[1]));
            }
            "repel_curve" => {
                if n == 3 || n == 4 {
                    let (a, b) = (num(&parts, 1, line, key)?, num(&parts, 2, line, key)?);
                    spec.exponents = Exponents::new(a, b)?;
                    if n == 4 {
                        spec.repel_weight = num(&parts, 3, line, key)?;
                    }
                } else if n != 1 {
                    return Err(perr(line, "repel_curve takes [alpha beta [weight]]"));
                }
            }
            "repel_surface" => {
                if !(2..=3).contains(&n) {
                    return Err(perr(line, "repel_surface takes path [weight]"));
                }
                let w = if n == 3 {
                    num(&parts, 2, line, key)?
                } else {
                    1.0
                };
                let path = base.join(parts[1]);
                spec.obstacles.push(ObstacleTerm {
                    mesh: Arc::new(Obstacle::read_obj(&path)?),
                    weight: w,
                    source: parts[1].to_string(),
                });
            }
            "repel_plane" => {
                if n != 7 && n != 8 {
                    return Err(perr(line, "repel_plane takes cx cy cz nx ny nz [weight]"));
                }
                let c = [
                    num(&parts, 1, line, key)?,
                    num(&parts, 2, line, key)?,
                    num(&parts, 3, line, key)?,
                ];
                let nn = [
                    num(&parts, 4, line, key)?,
                    num(&parts, 5, line, key)?,
                    num(&parts, 6, line, key)?,
                ];
                let w = if n == 8 {
                    num(&parts, 7, line, key)?
                } else {
                    1.0
                };
                planes.push((c, nn, w));
            }
            "show_surface" => {
                notes.push(format!("line {line}: show_surface ignored (display only)"))
            }
            "optimize_length" | "optimize_length_diffs" => {
                let w = if n == 2 {
                    num(&parts, 1, line, key)?
                } else {
                    1.0
                };
                if key == "optimize_length" {
                    spec.length_penalty += w;
                } else {
                    spec.length_diff_penalty += w;
                }
            }
            "optimize_pin_angles" | "optimize_area" => notes.push(format!(
                "line {line}: {key} is not supported and was ignored"
            )),
            "optimize_field" => {
                if n != 2 && n != 3 {
                    return Err(perr(line, "optimize_field takes fieldType [weight]"));
                }
                let w = if n == 3 {
                    num(&parts, 2, line, key)?
                } else {
                    1.0
                };
                spec.field = Some((parts[1].to_string(), w));
            }
            "fix_barycenter" => {
                spec.barycenter = true;
                any_constraint = true;
            }
            "fix_length" | "fix_totallength" | "fix_edgelengths" => {
                any_constraint = true;
                if !matches!(spec.length, LengthMode::Free) {
                    notes.push(format!(
                        "line {line}: {key} ignored, a length constraint is already set"
                    ));
                    continue;
                }
                let scale = if n == 2 {
                    num(&parts, 1, line, key)?
                } else {
                    1.0
                };
                if !(scale.is_finite() && scale > 0.0) {
                    return Err(perr(line, "the scale must be positive"));
                }
                spec.length = if key == "fix_edgelengths" {
                    LengthMode::EdgeLengths { scale }
                } else {
                    LengthMode::TotalLength { scale }
                };
            }
            "fix_vertex" => {
                spec.pins.push(idx(&parts, 1, line, key)?);
                any_constraint = true;
            }
            "fix_special_vertices" => {
                spec.pin_special = true;
                any_constraint = true;
            }
            "fix_endpoint_vertices" => {
                spec.pin_endpoints = true;
                any_constraint = true;
            }
            "fix_special_tangents" => {
                spec.special_tangents = true;
                any_constraint = true;
                notes.push(
                    "fix_special_tangents: tangents at endpoints only (a junction has no vertex tangent)"
                        .to_string(),
                );
            }
            "fix_tangent" => {
                spec.tangent_pins.push(idx(&parts, 1, line, key)?);
                any_constraint = true;
            }
            "iteration_limit" => {
                let k = idx(&parts, 1, line, key)? as usize;
                spec.iteration_limit = if k == 0 { None } else { Some(k) };
            }
            "subdivide_limit" => spec.subdivide_limit = idx(&parts, 1, line, key)?,
            "union_type" => match parts.get(1).copied() {
                Some("smooth") => smooth_union = true,
                Some("disjoint") => smooth_union = false,
                _ => return Err(perr(line, "union_type is `disjoint` or `smooth`")),
            },
            "constraint_surface" => {
                let kind = *parts
                    .get(1)
                    .ok_or_else(|| perr(line, "constraint_surface takes a type"))?;
                let f: Arc<dyn ImplicitSurface> = match (kind, n) {
                    ("sphere", 2) => Arc::new(Sphere::new([0.0; 3], 1.0)?),
                    ("sphere", 5) | ("sphere", 6) => {
                        let c = [
                            num(&parts, 2, line, key)?,
                            num(&parts, 3, line, key)?,
                            num(&parts, 4, line, key)?,
                        ];
                        let r = if n == 6 {
                            num(&parts, 5, line, key)?
                        } else {
                            1.0
                        };
                        Arc::new(Sphere::new(c, r)?)
                    }
                    ("torus", 2) => Arc::new(Torus::new([0.0; 3], 1.0, 0.25)?),
                    ("torus", 5) => {
                        let c = [
                            num(&parts, 2, line, key)?,
                            num(&parts, 3, line, key)?,
                            num(&parts, 4, line, key)?,
                        ];
                        Arc::new(Torus::new(c, 1.0, 0.25)?)
                    }
                    ("torus", 7) => {
                        let c = [
                            num(&parts, 2, line, key)?,
                            num(&parts, 3, line, key)?,
                            num(&parts, 4, line, key)?,
                        ];
                        Arc::new(Torus::new(
                            c,
                            num(&parts, 5, line, key)?,
                            num(&parts, 6, line, key)?,
                        )?)
                    }
                    ("yplane", 2) => Arc::new(Plane::y_zero()),
                    ("doubletorus", 2) => Arc::new(DoubleTorus::new(0.04)?),
                    ("sphere" | "torus" | "yplane" | "doubletorus", _) => {
                        return Err(perr(line, format!("wrong number of arguments to {kind}")));
                    }
                    _ => return Err(perr(line, format!("unknown surface type `{kind}`"))),
                };
                surface = Some(match surface.take() {
                    None => f,
                    Some(prev) if smooth_union => Arc::new(SmoothUnion::new(f, prev, 1.0)?),
                    Some(prev) => Arc::new(Union::new(vec![f, prev])?),
                });
                surface_names.push(kind.to_string());
            }
            "constrain_vertex" => {
                let v = idx(&parts, 1, line, key)?;
                match &mut spec.surface_vertices {
                    SurfaceVertices::List(l) => l.push(v),
                    s @ SurfaceVertices::None => *s = SurfaceVertices::List(vec![v]),
                    _ => {}
                }
                any_constraint = true;
            }
            "constrain_all" => {
                spec.surface_vertices = SurfaceVertices::All;
                any_constraint = true;
            }
            "constrain_endpoints" => {
                if !matches!(spec.surface_vertices, SurfaceVertices::All) {
                    spec.surface_vertices = SurfaceVertices::Endpoints;
                }
                any_constraint = true;
            }
            other => notes.push(format!("line {line}: unknown directive `{other}` ignored")),
        }
    }
    let curve = curve.ok_or_else(|| perr(0, "the scene has no `curve` line"))?;
    let (network, positions) = read_scene_curve(&curve)?;
    if let Some(f) = surface {
        let label = if surface_names.len() == 1 {
            surface_names[0].clone()
        } else {
            format!(
                "{}({})",
                if smooth_union {
                    "smooth-union"
                } else {
                    "union"
                },
                surface_names.join("+")
            )
        };
        // the plane special case is keyed on the plain name
        spec.surface = Some((f, label));
    }
    if !planes.is_empty() {
        let (lo, hi) = positions
            .iter()
            .fold(([f64::MAX; 3], [f64::MIN; 3]), |(lo, hi), p| {
                (
                    [lo[0].min(p[0]), lo[1].min(p[1]), lo[2].min(p[2])],
                    [hi[0].max(p[0]), hi[1].max(p[1]), hi[2].max(p[2])],
                )
            });
        let r = 20.0 * dist(lo, hi).max(1e-9);
        for (c, nn, w) in planes {
            spec.obstacles.push(ObstacleTerm {
                mesh: Arc::new(plane_obstacle(c, nn, r)?),
                weight: w,
                source: "repel_plane".to_string(),
            });
        }
        notes.push("repel_plane is a square of side 40 times the curve's diameter".to_string());
    }
    if spec.repel_weight != 1.0 {
        notes.push(format!(
            "repel_curve weight {}: the other penalty weights are divided by it",
            spec.repel_weight
        ));
    }
    if !any_constraint {
        notes.push("no constraint directive: barycenter and edge lengths (default)".into());
    }
    if spec.surface.is_none() && !matches!(spec.surface_vertices, SurfaceVertices::None) {
        return Err(perr(0, "constrain_* without a constraint_surface"));
    }
    // The authors' curves lie on their surfaces only to the precision of the file; the
    // constraints must hold at the start (`SPEC.md` §5.1), so move the vertices onto it.
    let mut positions = positions;
    if let Some((f, _)) = &spec.surface {
        let ids: Vec<u32> = match &spec.surface_vertices {
            SurfaceVertices::None => Vec::new(),
            SurfaceVertices::All => (0..network.n_vertices() as u32).collect(),
            SurfaceVertices::Endpoints => (0..network.n_vertices() as u32)
                .filter(|&v| network.degree(v) == 1)
                .collect(),
            SurfaceVertices::List(l) => l
                .iter()
                .copied()
                .filter(|&v| (v as usize) < positions.len())
                .collect(),
        };
        let mut worst = 0.0f64;
        for v in ids {
            let x0 = positions[v as usize];
            let mut x = x0;
            for _ in 0..60 {
                let (fx, g) = (f.value(x), f.gradient(x));
                let g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
                if g2 < 1e-300 || fx.abs() < 1e-14 * g2.sqrt() {
                    break;
                }
                x = [
                    x[0] - fx * g[0] / g2,
                    x[1] - fx * g[1] / g2,
                    x[2] - fx * g[2] / g2,
                ];
            }
            worst = worst.max(dist(x, x0));
            positions[v as usize] = x;
        }
        if worst > 0.0 {
            notes.push(format!(
                "vertices moved onto the constraint surface by at most {worst:.3e}"
            ));
        }
    }
    let mut scene = Scene::new(name, network, positions, spec)?;
    scene.notes = notes;
    Ok(scene)
}

/// Reads a scene file; the name is the parent directory's name plus the file stem when the stem
/// is not `scene`.
///
/// # Errors
/// [`RcError::Io`] and the errors of [`parse_scene`].
pub fn read_scene(path: impl AsRef<Path>) -> Result<Scene, RcError> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path).map_err(|e| RcError::Io {
        path: path.display().to_string(),
        message: e.to_string(),
    })?;
    let base = path.parent().unwrap_or_else(|| Path::new("."));
    let dir = base
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_default();
    let stem = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_default();
    let name = if stem == "scene" {
        dir
    } else {
        format!("{dir}-{stem}")
    };
    parse_scene(&text, base, &name)
}

// =============================================================================================
// Audit
// =============================================================================================

/// Geometric and constraint audit of a curve.
#[derive(Debug, Clone, PartialEq)]
pub struct Audit {
    /// Smallest distance between two edges that share no vertex (`INFINITY` if none).
    pub min_edge_distance: f64,
    /// Mean edge length.
    pub mean_edge_length: f64,
    /// Total length.
    pub total_length: f64,
    /// Largest distance of a pinned vertex from its start position.
    pub max_pin_error: f64,
    /// Largest `sin` of the angle between a constrained tangent and its direction.
    pub max_tangent_error: f64,
    /// Largest `|f| / |∇f|` (a distance estimate) over the surface-constrained vertices.
    pub max_surface_error: f64,
    /// Largest absolute value of a row of the library's constraint residual `Φ`
    /// (`None` when no session was given).
    pub max_constraint_residual: Option<f64>,
    /// Smallest distance of a vertex or edge midpoint to an obstacle (`None` without obstacles).
    pub min_obstacle_distance: Option<f64>,
}

/// Smallest distance between two non-adjacent edges, in parallel; with the pair.
pub fn min_edge_distance(network: &CurveNetwork, positions: &[P]) -> (f64, (usize, usize)) {
    let e = network.edges();
    let m = e.len();
    (0..m)
        .into_par_iter()
        .map(|i| {
            let [a, b] = e[i];
            let (pa, pb) = (positions[a as usize], positions[b as usize]);
            let mut best = (f64::INFINITY, (i, i));
            for j in i + 1..m {
                let [c, d] = e[j];
                if a == c || a == d || b == c || b == d {
                    continue;
                }
                let dd =
                    segment_segment_distance(pa, pb, positions[c as usize], positions[d as usize]);
                if dd < best.0 {
                    best = (dd, (i, j));
                }
            }
            best
        })
        .reduce(
            || (f64::INFINITY, (0, 0)),
            |a, b| if b.0 < a.0 { b } else { a },
        )
}

/// Audits `positions` against the scene's constraints; `session` adds the library's own
/// residual.
pub fn audit(
    scene: &Scene,
    network: &CurveNetwork,
    positions: &[P],
    obstacles: &[Arc<Obstacle>],
    session: Option<&Session>,
) -> Audit {
    let g = network.edge_geometry(positions);
    let total: f64 = g.length.iter().sum();
    let r = &scene.resolved;
    let max_pin_error = r
        .pins
        .iter()
        .map(|(v, p)| dist(positions[*v as usize], *p))
        .fold(0.0, f64::max);
    let max_tangent_error = r
        .tangents
        .iter()
        .filter_map(|(v, x)| {
            vertex_tangent(network, positions, *v).map(|s| {
                let c = [
                    s[1] * x[2] - s[2] * x[1],
                    s[2] * x[0] - s[0] * x[2],
                    s[0] * x[1] - s[1] * x[0],
                ];
                // sine of the angle, or 2 for an antiparallel tangent
                let dot = s[0] * x[0] + s[1] * x[1] + s[2] * x[2];
                if dot < 0.0 {
                    2.0
                } else {
                    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
                }
            })
        })
        .fold(0.0, f64::max);
    let max_surface_error = match &r.surface {
        Some(f) => r
            .surface_vertices
            .iter()
            .map(|&v| {
                let x = positions[v as usize];
                let gr = f.gradient(x);
                let gn = (gr[0] * gr[0] + gr[1] * gr[1] + gr[2] * gr[2])
                    .sqrt()
                    .max(1e-300);
                f.value(x).abs() / gn
            })
            .fold(0.0, f64::max),
        None => 0.0,
    };
    let max_constraint_residual = session.and_then(|s| {
        s.constraints()
            .values(s.network(), s.positions())
            .ok()
            .map(|v| v.iter().fold(0.0, |m: f64, x| m.max(x.abs())))
    });
    let min_obstacle_distance = if obstacles.is_empty() {
        None
    } else {
        let e = network.edges();
        let mut pts: Vec<P> = positions.to_vec();
        for [a, b] in e {
            let (p, q) = (positions[*a as usize], positions[*b as usize]);
            pts.push([
                0.5 * (p[0] + q[0]),
                0.5 * (p[1] + q[1]),
                0.5 * (p[2] + q[2]),
            ]);
        }
        Some(
            obstacles
                .iter()
                .map(|o| {
                    pts.par_iter()
                        .map(|&x| o.distance_to(x).0)
                        .reduce(|| f64::INFINITY, f64::min)
                })
                .fold(f64::INFINITY, f64::min),
        )
    };
    Audit {
        min_edge_distance: min_edge_distance(network, positions).0,
        mean_edge_length: total / network.n_edges().max(1) as f64,
        total_length: total,
        max_pin_error,
        max_tangent_error,
        max_surface_error,
        max_constraint_residual,
        min_obstacle_distance,
    }
}

// =============================================================================================
// Running a scene
// =============================================================================================

/// Options of [`run_scene`].
#[derive(Debug, Clone)]
pub struct RunOptions {
    /// Step budget (`None`: the scene's iteration limit, else the crate default).
    pub steps: Option<usize>,
    /// Back end.
    pub backend: BackendChoice,
    /// Threads of the main pool (`None`: all cores).
    pub threads: Option<usize>,
    /// Write `frame%05d.obj` every this many steps (0: no frames; the first and last are
    /// always written when `out_dir` is set).
    pub frame_every: usize,
    /// Output directory (frames, `report.json`).
    pub out_dir: Option<PathBuf>,
    /// Audit every this many steps (0: only first and last). Audits cost `O(|E|²)`.
    pub audit_every: usize,
}

impl Default for RunOptions {
    fn default() -> Self {
        RunOptions {
            steps: None,
            backend: BackendChoice::Auto,
            threads: None,
            frame_every: 0,
            out_dir: None,
            audit_every: 0,
        }
    }
}

/// One line of the step log.
#[derive(Debug, Clone, PartialEq)]
pub struct StepLog {
    /// Step number (1-based).
    pub step: usize,
    /// Vertices after the step.
    pub vertices: usize,
    /// Energy after the step (tangent-point plus penalties).
    pub energy: f64,
    /// Step size.
    pub tau: f64,
    /// Largest constraint violation after projection.
    pub max_violation: f64,
    /// Whether the step subdivided the curve.
    pub subdivided: bool,
    /// Seconds of the step.
    pub seconds: f64,
}

/// What [`run_scene`] reports.
#[derive(Debug, Clone)]
pub struct RunReport {
    /// Scene name.
    pub scene: String,
    /// Back end names used (energy, solver).
    pub backend: String,
    /// Threads.
    pub threads: Option<usize>,
    /// Steps taken.
    pub steps: usize,
    /// Final status.
    pub status: String,
    /// Wall seconds of the flow (frames and audits excluded).
    pub seconds: f64,
    /// Vertices and edges at the start and the end.
    pub vertices: (usize, usize),
    /// Edges at the start and the end.
    pub edges: (usize, usize),
    /// Energy before the first and after the last step.
    pub energy: (f64, f64),
    /// Audit of the start curve.
    pub audit_start: Audit,
    /// Audit of the end curve.
    pub audit_end: Audit,
    /// Smallest `min_edge_distance` over the audited steps.
    pub min_distance_over_run: f64,
    /// Total length of the start curve and of the end curve.
    pub length: (f64, f64),
    /// The log (one entry per step).
    pub log: Vec<StepLog>,
    /// The notes of the scene.
    pub notes: Vec<String>,
    /// Frames written.
    pub frames: usize,
}

fn jf(x: f64) -> String {
    if x.is_finite() {
        format!("{x:e}")
    } else if x.is_nan() {
        "null".into()
    } else if x > 0.0 {
        "1e999".into()
    } else {
        "-1e999".into()
    }
}

fn jo(x: Option<f64>) -> String {
    x.map_or("null".to_string(), jf)
}

fn js(s: &str) -> String {
    let mut o = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            c if (c as u32) < 0x20 => {
                let _ = write!(o, "\\u{:04x}", c as u32);
            }
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

impl Audit {
    /// JSON object.
    pub fn to_json(&self) -> String {
        format!(
            "{{\"min_edge_distance\":{},\"mean_edge_length\":{},\"total_length\":{},\"max_pin_error\":{},\"max_tangent_error\":{},\"max_surface_error\":{},\"max_constraint_residual\":{},\"min_obstacle_distance\":{}}}",
            jf(self.min_edge_distance),
            jf(self.mean_edge_length),
            jf(self.total_length),
            jf(self.max_pin_error),
            jf(self.max_tangent_error),
            jf(self.max_surface_error),
            jo(self.max_constraint_residual),
            jo(self.min_obstacle_distance)
        )
    }
}

impl RunReport {
    /// The report as JSON (hand-written; the crate has no serializer dependency).
    pub fn to_json(&self) -> String {
        let mut s = String::from("{\n");
        let _ = writeln!(s, "  \"scene\": {},", js(&self.scene));
        let _ = writeln!(s, "  \"backend\": {},", js(&self.backend));
        let _ = writeln!(
            s,
            "  \"threads\": {},",
            self.threads.map_or("null".to_string(), |t| t.to_string())
        );
        let _ = writeln!(s, "  \"steps\": {},", self.steps);
        let _ = writeln!(s, "  \"status\": {},", js(&self.status));
        let _ = writeln!(s, "  \"seconds\": {},", jf(self.seconds));
        let _ = writeln!(
            s,
            "  \"vertices\": [{}, {}],\n  \"edges\": [{}, {}],",
            self.vertices.0, self.vertices.1, self.edges.0, self.edges.1
        );
        let _ = writeln!(
            s,
            "  \"energy\": [{}, {}],",
            jf(self.energy.0),
            jf(self.energy.1)
        );
        let _ = writeln!(
            s,
            "  \"length\": [{}, {}],",
            jf(self.length.0),
            jf(self.length.1)
        );
        let _ = writeln!(
            s,
            "  \"min_distance_over_run\": {},",
            jf(self.min_distance_over_run)
        );
        let _ = writeln!(s, "  \"audit_start\": {},", self.audit_start.to_json());
        let _ = writeln!(s, "  \"audit_end\": {},", self.audit_end.to_json());
        let _ = writeln!(s, "  \"frames\": {},", self.frames);
        let notes: Vec<String> = self.notes.iter().map(|n| js(n)).collect();
        let _ = writeln!(s, "  \"notes\": [{}],", notes.join(", "));
        s.push_str("  \"log\": [\n");
        for (i, l) in self.log.iter().enumerate() {
            let _ = writeln!(
                s,
                "    [{}, {}, {}, {}, {}, {}, {}]{}",
                l.step,
                l.vertices,
                jf(l.energy),
                jf(l.tau),
                jf(l.max_violation),
                l.subdivided as u8,
                jf(l.seconds),
                if i + 1 < self.log.len() { "," } else { "" }
            );
        }
        s.push_str("  ],\n  \"log_columns\": [\"step\", \"vertices\", \"energy\", \"tau\", \"max_violation\", \"subdivided\", \"seconds\"]\n}\n");
        s
    }
}

/// The result of [`run_scene`].
pub struct RunOutcome {
    /// The report.
    pub report: RunReport,
    /// Final network (subdivision changes it).
    pub network: CurveNetwork,
    /// Final positions.
    pub positions: Vec<P>,
}

fn io_err(path: &Path, e: std::io::Error) -> RcError {
    RcError::Io {
        path: path.display().to_string(),
        message: e.to_string(),
    }
}

/// Runs a scene to convergence or the step budget.
///
/// # Errors
/// Errors of the session; [`RcError::Io`] for output files.
pub fn run_scene(scene: &Scene, opts: &RunOptions) -> Result<RunOutcome, RcError> {
    let mut options = scene.options(opts.backend, opts.threads)?;
    if let Some(n) = opts.steps {
        options.max_iterations = n;
    }
    let max_steps = options.max_iterations;
    let mut session = Session::new(scene.network.clone(), scene.positions.clone(), options)?;
    let obstacles: Vec<Arc<Obstacle>> = scene
        .spec
        .obstacles
        .iter()
        .map(|o| o.mesh.clone())
        .collect();
    if let Some(d) = &opts.out_dir {
        std::fs::create_dir_all(d).map_err(|e| io_err(d, e))?;
    }
    let mut frames = 0usize;
    let write_frame = |s: &Session, frames: &mut usize| -> Result<(), RcError> {
        if let Some(d) = &opts.out_dir {
            write_curve_obj(
                d.join(format!("frame{:05}.obj", *frames)),
                s.network(),
                s.positions(),
            )?;
            *frames += 1;
        }
        Ok(())
    };
    let audit_start = audit(
        scene,
        session.network(),
        session.positions(),
        &obstacles,
        Some(&session),
    );
    let mut min_over = audit_start.min_edge_distance;
    let v0 = (scene.network.n_vertices(), scene.network.n_edges());
    let e_start = session.energy();
    write_frame(&session, &mut frames)?;
    let mut log = Vec::new();
    let mut seconds = 0.0;
    let mut status = Status::Running;
    let mut last_frame_step = 0usize;
    let mut error: Option<String> = None;
    for k in 1..=max_steps {
        let t0 = Instant::now();
        let r = match session.step() {
            Ok(r) => r,
            Err(e) => {
                // keep what was computed: the report says why the run stopped
                error = Some(e.to_string());
                break;
            }
        };
        seconds += t0.elapsed().as_secs_f64();
        log.push(StepLog {
            step: k,
            vertices: r.n_vertices,
            energy: r.energy_after,
            tau: r.tau,
            max_violation: r.max_violation,
            subdivided: r.subdivided,
            seconds: r.seconds.total,
        });
        if opts.audit_every > 0 && k % opts.audit_every == 0 {
            min_over = min_over.min(min_edge_distance(session.network(), session.positions()).0);
        }
        if opts.frame_every > 0 && k % opts.frame_every == 0 {
            write_frame(&session, &mut frames)?;
            last_frame_step = k;
        }
        status = session.status();
        if status.is_terminal() {
            break;
        }
    }
    let steps = log.len();
    if opts.out_dir.is_some() && last_frame_step != steps {
        write_frame(&session, &mut frames)?;
    }
    let audit_end = audit(
        scene,
        session.network(),
        session.positions(),
        &obstacles,
        Some(&session),
    );
    min_over = min_over.min(audit_end.min_edge_distance);
    let (eb, sb) = log
        .last()
        .map(|_| {
            let r = session
                .options()
                .backend
                .resolve(session.network().n_vertices());
            (format!("{:?}", r.0), format!("{:?}", r.1))
        })
        .unwrap_or_default();
    let report = RunReport {
        scene: scene.name.clone(),
        backend: format!("{eb} + {sb}"),
        threads: opts.threads,
        steps,
        status: match &error {
            Some(e) => format!("Error: {e}"),
            None => format!("{status:?}"),
        },
        seconds,
        vertices: (v0.0, session.network().n_vertices()),
        edges: (v0.1, session.network().n_edges()),
        energy: (e_start, session.energy()),
        length: (audit_start.total_length, audit_end.total_length),
        audit_start,
        audit_end,
        min_distance_over_run: min_over,
        log,
        notes: scene.notes.clone(),
        frames,
    };
    if let Some(d) = &opts.out_dir {
        let p = d.join("report.json");
        std::fs::write(&p, report.to_json()).map_err(|e| io_err(&p, e))?;
    }
    Ok(RunOutcome {
        report,
        network: session.network().clone(),
        positions: session.positions().to_vec(),
    })
}

pub mod builtin;

/// Reads the curve file of a scene: [`crate::io::read_curve_obj`], except that a file without `l` lines and
/// with exactly one polygon face (`f i j k ...`, the authors' `DoubleTorus/InitialLoop.obj` and
/// `FreedmanUnknot/FreedmanKnotPlot.obj`) is read as the closed loop through the face's vertices
/// in that order.
///
/// # Errors
/// [`RcError::Io`]; as [`parse_curve_obj`]; [`RcError::Parse`] for a malformed face.
pub fn read_scene_curve(
    path: impl AsRef<std::path::Path>,
) -> Result<(CurveNetwork, Vec<P>), RcError> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path).map_err(|e| RcError::Io {
        path: path.display().to_string(),
        message: e.to_string(),
    })?;
    parse_scene_curve(&text).map_err(|e| match e {
        RcError::Parse { line, message } => RcError::Parse {
            line,
            message: format!("{}: {message}", path.display()),
        },
        e => e,
    })
}

/// [`read_scene_curve`] on text.
///
/// # Errors
/// As [`read_scene_curve`], without the I/O error.
pub fn parse_scene_curve(text: &str) -> Result<(CurveNetwork, Vec<P>), RcError> {
    let is_tag = |l: &str, t: &str| l.split_whitespace().next() == Some(t);
    let has_lines = text.lines().any(|l| is_tag(l, "l"));
    let faces: Vec<(usize, &str)> = text
        .lines()
        .enumerate()
        .filter(|(_, l)| is_tag(l, "f"))
        .collect();
    if has_lines || faces.len() != 1 {
        return parse_curve_obj(text);
    }
    let (ln, face) = faces[0];
    let first = face
        .split_whitespace()
        .nth(1)
        .ok_or_else(|| RcError::Parse {
            line: ln + 1,
            message: "polygon face without vertices".into(),
        })?;
    let mut loop_line = String::from("l");
    for tok in face.split_whitespace().skip(1) {
        loop_line.push(' ');
        loop_line.push_str(tok);
    }
    loop_line.push(' ');
    loop_line.push_str(first);
    let mut out = String::with_capacity(text.len() + loop_line.len() + 1);
    for (k, l) in text.lines().enumerate() {
        // Keep the line numbers: the face line becomes the loop.
        out.push_str(if k == ln { &loop_line } else { l });
        out.push('\n');
    }
    parse_curve_obj(&out)
}

/// Parses a curve OBJ given as text (convenience for generated scenes).
///
/// # Errors
/// As [`parse_curve_obj`].
pub fn curve_from_text(text: &str) -> Result<(CurveNetwork, Vec<P>), RcError> {
    parse_curve_obj(text)
}

impl std::fmt::Debug for Scene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Scene")
            .field("name", &self.name)
            .field("vertices", &self.network.n_vertices())
            .field("edges", &self.network.n_edges())
            .field("spec", &self.spec)
            .field("notes", &self.notes)
            .finish_non_exhaustive()
    }
}
