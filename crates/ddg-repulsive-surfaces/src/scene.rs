//! Reader for the reference's scene-file format (`SPEC.md` §14.4, WP12).
//!
//! Two levels:
//!
//! * [`SceneFile`] is the syntax of a scene file, one field per keyword, with the line of every
//!   statement. [`SceneFile::parse`] never touches the file system; [`SceneFile::to_text`]
//!   writes a canonical scene that parses back to an equal value.
//! * [`SceneSpec`] is the scene resolved against its files: the mesh is read and validated,
//!   vertex sets are computed (boundary vertices, UV-selected attractor vertices), every
//!   potential is turned into a [`Penalty`] with its targets, every obstacle into a
//!   [`Penalty::MeshObstacle`] or [`Penalty::PointCloudObstacle`], every implicit surface into
//!   a penalty with [`Quadrature::PerVertex`] (the reference's convention, D21). The flow
//!   session (WP5) turns it into `FlowOptions`.
//!
//! # Grammar
//!
//! `FORMAT.txt` of the reference plus its undocumented keywords. Lines are split at white
//! space (the reference splits at single spaces); empty lines and lines whose first
//! non-blank character is `#` are skipped. Unknown keywords, unknown arguments, missing
//! arguments and surplus arguments are errors ([`RsError::Parse`], 1-based line; line 0 for
//! a problem of the whole file). Where the reference silently ignores or leaves values
//! uninitialized, this reader is stricter:
//!
//! | statement | notes |
//! |---|---|
//! | `repel_mesh <obj> <α> <β>` | required, once; the reference leaves `α β` uninitialized when they are missing (`Boundary/scene-tall.txt`), here an error; `(α, β)` must satisfy [`Exponents::new`] |
//! | `constrain barycenter` | `Constraint::Barycenter { component: None, target: None }` |
//! | `constrain area\|volume [mult iters [add]]` | growth to `mult·X₀ + add` over `iters` steps; `mult` without `iters` is an error (the reference ignores it) |
//! | `constrain vertices\|vertex i j … [move dx dy dz iters]` | `Pins` to the start positions; `move` gives a [`ScheduleSpec::Move`] (the reference's moving pins) |
//! | `constrain boundary_vertices` | `Pins` of all boundary vertices, computed at load time |
//! | `constrain boundary_normals`, `constrain normal\|normals i j …` | [`NormalPins`]; the library has no normal constraint (the reference falls through and never implemented the vertex form) |
//! | `minimize squared_error\|area\|volume\|area_deviation\|volume_deviation\|boundary_length\|boundary_curvature\|willmore [w [target]]` | `w` default 1, `target` default 0 (only `boundary_length` uses it). `area` is `TotalArea { power: 2 }` and `volume` `TotalVolume { power: 2 }` (the reference squares them); `boundary_curvature` is `dual_length: false` |
//! | `obstacle <obj> [w [centered]]`, `point_cloud_obstacle <obj> [w [centered]]` | `w` default 1; `centered` subtracts the area-weighted barycenter; mesh obstacles are [`ObstacleKind::TangentPoint`] |
//! | `implicit sphere\|cylinder\|torus\|plane repel\|attract <power> <weight> <params>` | exact parameter counts; the torus has the reference's fixed `y` axis |
//! | `iteration_limit n`, `time_limit ms` | `0` means no limit (the reference's default) |
//! | `allow_barycenter_shift`, `disable_near_field` | flags |
//! | `method hs\|aqp\|bqn\|h1\|h1-lbfgs\|h2\|willmore\|l2` | [`SceneMethod`] |
//! | `log <csv>` | path relative to the scene file |
//! | `autotarget_volume [ratio]` | `ratio` default 1 (the reference reads it unconditionally) |
//!
//! Paths (mesh, obstacles, log) are relative to the directory of the scene file.
//!
//! The reference collects all `constrain vertices` lines of a scene into one pin constraint;
//! here every line is one `Pins` constraint with its own motion.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use ddg_repulsive_curves::obstacle::Obstacle;

use crate::constraints::Constraint;
use crate::error::RsError;
use crate::geometry::v3::{cross, dot, sub};
use crate::implicit::{Cylinder, FramedTorus, ImplicitSurface, Plane, Sphere};
use crate::kernel::Exponents;
use crate::obstacle::{ObstacleKind, StaticSet};
use crate::penalties::simple::{area_and_gradient, volume_and_gradient};
use crate::penalties::{Penalty, Quadrature};
use crate::surface::Surface;

/// Leaf size of the obstacle BVHs created by the loader.
pub const OBSTACLE_LEAF_SIZE: usize = 8;

fn parse_err(path: &str, line: usize, message: impl Into<String>) -> RsError {
    RsError::Parse {
        path: Some(path.to_string()),
        line,
        message: message.into(),
    }
}

/// The reference's `method` names.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SceneMethod {
    /// `hs`: the fractional Sobolev gradient (default).
    #[default]
    Hs,
    /// `aqp`.
    Aqp,
    /// `bqn`: BCQN.
    Bqn,
    /// `h1`.
    H1,
    /// `h1-lbfgs`.
    H1Lbfgs,
    /// `h2`.
    H2,
    /// `willmore`: the tangent-point energy is replaced by the Willmore energy (H² flow).
    Willmore,
    /// `l2`.
    L2,
}

impl SceneMethod {
    /// The scene-file spelling.
    pub fn name(self) -> &'static str {
        match self {
            SceneMethod::Hs => "hs",
            SceneMethod::Aqp => "aqp",
            SceneMethod::Bqn => "bqn",
            SceneMethod::H1 => "h1",
            SceneMethod::H1Lbfgs => "h1-lbfgs",
            SceneMethod::H2 => "h2",
            SceneMethod::Willmore => "willmore",
            SceneMethod::L2 => "l2",
        }
    }

    fn from_name(s: &str) -> Option<SceneMethod> {
        Some(match s {
            "hs" => SceneMethod::Hs,
            "aqp" => SceneMethod::Aqp,
            "bqn" => SceneMethod::Bqn,
            "h1" => SceneMethod::H1,
            "h1-lbfgs" => SceneMethod::H1Lbfgs,
            "h2" => SceneMethod::H2,
            "willmore" => SceneMethod::Willmore,
            "l2" => SceneMethod::L2,
            _ => return None,
        })
    }
}

/// A growth schedule `constrain area 3 200 1.5`: the target moves linearly from `X₀` to
/// `multiplier·X₀ + addition` over `iterations` steps (`SPEC.md` D11).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Growth {
    /// The factor on the initial value.
    pub multiplier: f64,
    /// Number of steps, at least 1.
    pub iterations: u64,
    /// Added to the final target (0 if absent).
    pub addition: f64,
}

impl Growth {
    /// `multiplier·x0 + addition`.
    pub fn final_target(&self, initial: f64) -> f64 {
        self.multiplier * initial + self.addition
    }
}

/// `constrain vertices … move dx dy dz iters`: the pins move by `displacement` over `iterations`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PinMotion {
    /// The total displacement of every pin of the statement.
    pub displacement: [f64; 3],
    /// Number of steps, at least 1.
    pub iterations: u64,
}

/// One `constrain …` statement.
#[derive(Debug, Clone, PartialEq)]
pub enum ConstrainStmt {
    /// `constrain barycenter`.
    Barycenter,
    /// `constrain area [growth]`.
    Area(Option<Growth>),
    /// `constrain volume [growth]`.
    Volume(Option<Growth>),
    /// `constrain vertices i j … [move …]`.
    Vertices {
        /// The pinned vertices, as written.
        vertices: Vec<u32>,
        /// The moving-pin clause.
        motion: Option<PinMotion>,
    },
    /// `constrain boundary_vertices`.
    BoundaryVertices,
    /// `constrain boundary_normals`.
    BoundaryNormals,
    /// `constrain normal i j …`.
    Normals(Vec<u32>),
}

/// The potentials of `minimize …`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MinimizeKind {
    /// `squared_error`.
    SquaredError,
    /// `area` (squared total area).
    Area,
    /// `volume` (squared total volume).
    Volume,
    /// `area_deviation`.
    AreaDeviation,
    /// `volume_deviation`.
    VolumeDeviation,
    /// `boundary_length`.
    BoundaryLength,
    /// `boundary_curvature`.
    BoundaryCurvature,
    /// `willmore`.
    Willmore,
}

impl MinimizeKind {
    fn name(self) -> &'static str {
        match self {
            MinimizeKind::SquaredError => "squared_error",
            MinimizeKind::Area => "area",
            MinimizeKind::Volume => "volume",
            MinimizeKind::AreaDeviation => "area_deviation",
            MinimizeKind::VolumeDeviation => "volume_deviation",
            MinimizeKind::BoundaryLength => "boundary_length",
            MinimizeKind::BoundaryCurvature => "boundary_curvature",
            MinimizeKind::Willmore => "willmore",
        }
    }

    fn from_name(s: &str) -> Option<MinimizeKind> {
        Some(match s {
            "squared_error" => MinimizeKind::SquaredError,
            "area" => MinimizeKind::Area,
            "volume" => MinimizeKind::Volume,
            "area_deviation" => MinimizeKind::AreaDeviation,
            "volume_deviation" => MinimizeKind::VolumeDeviation,
            "boundary_length" => MinimizeKind::BoundaryLength,
            "boundary_curvature" => MinimizeKind::BoundaryCurvature,
            "willmore" => MinimizeKind::Willmore,
            _ => return None,
        })
    }
}

/// `minimize kind [weight [target]]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Minimize {
    /// Which potential.
    pub kind: MinimizeKind,
    /// The weight, default 1.
    pub weight: f64,
    /// The target value, default 0 (used by `boundary_length`).
    pub target: f64,
}

/// `obstacle` / `point_cloud_obstacle`.
#[derive(Debug, Clone, PartialEq)]
pub struct ObstacleStmt {
    /// The OBJ path as written (relative to the scene file).
    pub path: String,
    /// The weight, default 1.
    pub weight: f64,
    /// Recentre on the area-weighted barycenter.
    pub centered: bool,
    /// `point_cloud_obstacle`: only the vertices are used.
    pub point_cloud: bool,
}

/// The parameters of an `implicit` statement.
#[derive(Debug, Clone, PartialEq)]
pub enum ImplicitDesc {
    /// `sphere <radius> <cx cy cz>`.
    Sphere {
        /// Radius.
        radius: f64,
        /// Centre.
        center: [f64; 3],
    },
    /// `cylinder <radius> <cx cy cz> <ax ay az>` (infinite).
    Cylinder {
        /// Radius.
        radius: f64,
        /// A point on the axis.
        center: [f64; 3],
        /// The axis direction.
        axis: [f64; 3],
    },
    /// `torus <major> <minor> <cx cy cz>` (axis `y`).
    Torus {
        /// Major radius.
        major: f64,
        /// Minor radius.
        minor: f64,
        /// Centre.
        center: [f64; 3],
    },
    /// `plane <px py pz> <nx ny nz>`.
    Plane {
        /// A point of the plane.
        point: [f64; 3],
        /// The normal.
        normal: [f64; 3],
    },
}

/// `implicit …`.
#[derive(Debug, Clone, PartialEq)]
pub struct ImplicitStmt {
    /// The surface.
    pub surface: ImplicitDesc,
    /// `repel` (obstacle, `d^{-p}`) or `attract` (`|d|^p`).
    pub repel: bool,
    /// The power.
    pub power: f64,
    /// The weight.
    pub weight: f64,
}

/// A statement with the line it came from.
#[derive(Debug, Clone, PartialEq)]
pub struct At<T> {
    /// 1-based line.
    pub line: usize,
    /// The statement.
    pub value: T,
}

/// The syntax of a scene file (the reference's `SceneData`, before any file is read).
#[derive(Debug, Clone, PartialEq)]
pub struct SceneFile {
    /// The path of the scene file (for messages and relative paths); empty for text.
    pub path: PathBuf,
    /// `repel_mesh`: the OBJ as written.
    pub mesh: String,
    /// `α`.
    pub alpha: f64,
    /// `β`.
    pub beta: f64,
    /// `constrain …`, in file order.
    pub constraints: Vec<At<ConstrainStmt>>,
    /// `minimize …`, in file order.
    pub potentials: Vec<At<Minimize>>,
    /// `obstacle` and `point_cloud_obstacle`, in file order.
    pub obstacles: Vec<At<ObstacleStmt>>,
    /// `implicit …`, in file order.
    pub implicits: Vec<At<ImplicitStmt>>,
    /// `iteration_limit` (`None`: 0 or absent).
    pub iteration_limit: Option<u64>,
    /// `time_limit` in milliseconds (`None`: 0 or absent).
    pub time_limit_ms: Option<u64>,
    /// `allow_barycenter_shift`.
    pub allow_barycenter_shift: bool,
    /// `disable_near_field`.
    pub disable_near_field: bool,
    /// `method`.
    pub method: SceneMethod,
    /// `log`, as written.
    pub log: Option<String>,
    /// `autotarget_volume`, the ratio.
    pub autotarget_volume: Option<f64>,
}

struct Tokens<'a> {
    path: &'a str,
    line: usize,
    keyword: &'a str,
    rest: Vec<&'a str>,
    next: usize,
}

impl<'a> Tokens<'a> {
    fn err(&self, message: impl Into<String>) -> RsError {
        parse_err(self.path, self.line, message)
    }

    fn remaining(&self) -> usize {
        self.rest.len() - self.next
    }

    fn peek(&self) -> Option<&'a str> {
        self.rest.get(self.next).copied()
    }

    fn word(&mut self, what: &str) -> Result<&'a str, RsError> {
        let w = self
            .rest
            .get(self.next)
            .copied()
            .ok_or_else(|| self.err(format!("`{}` needs {what}", self.keyword)))?;
        self.next += 1;
        Ok(w)
    }

    fn num(&mut self, what: &str) -> Result<f64, RsError> {
        let w = self.word(what)?;
        match w.parse::<f64>() {
            Ok(x) if x.is_finite() => Ok(x),
            _ => Err(self.err(format!("{what}: `{w}` is not a finite number"))),
        }
    }

    fn uint(&mut self, what: &str) -> Result<u64, RsError> {
        let w = self.word(what)?;
        w.parse::<u64>()
            .map_err(|_| self.err(format!("{what}: `{w}` is not a non-negative integer")))
    }

    fn vec3(&mut self, what: &str) -> Result<[f64; 3], RsError> {
        Ok([self.num(what)?, self.num(what)?, self.num(what)?])
    }

    fn finish(&self) -> Result<(), RsError> {
        if self.next < self.rest.len() {
            Err(self.err(format!(
                "`{}`: unexpected `{}`",
                self.keyword, self.rest[self.next]
            )))
        } else {
            Ok(())
        }
    }
}

fn growth_of(t: &mut Tokens) -> Result<Option<Growth>, RsError> {
    if t.remaining() == 0 {
        return Ok(None);
    }
    let multiplier = t.num("growth multiplier")?;
    if multiplier <= 0.0 {
        return Err(t.err("the growth multiplier must be positive"));
    }
    if t.remaining() == 0 {
        return Err(t.err(
            "a growth multiplier needs the number of growth iterations (`constrain area 3 200`)",
        ));
    }
    let iterations = t.uint("growth iterations")?;
    if iterations == 0 {
        return Err(t.err("the number of growth iterations must be at least 1"));
    }
    let addition = if t.remaining() > 0 {
        t.num("growth addition")?
    } else {
        0.0
    };
    Ok(Some(Growth {
        multiplier,
        iterations,
        addition,
    }))
}

fn vertex_list(t: &mut Tokens) -> Result<(Vec<u32>, Option<PinMotion>), RsError> {
    let mut vertices = Vec::new();
    let mut motion = None;
    while let Some(w) = t.peek() {
        if w == "move" {
            t.next += 1;
            let displacement = t.vec3("the displacement after `move`")?;
            let iterations = t.uint("the number of iterations of `move`")?;
            if iterations == 0 {
                return Err(t.err("`move` needs at least 1 iteration"));
            }
            motion = Some(PinMotion {
                displacement,
                iterations,
            });
            break;
        }
        let v = t.word("a vertex")?;
        vertices.push(v.parse::<u32>().map_err(|_| {
            t.err(format!("`{v}` is not a vertex index (non-negative integer)"))
        })?);
    }
    Ok((vertices, motion))
}

impl SceneFile {
    /// Parses scene text. `path` is only used in error messages and for relative paths.
    ///
    /// # Errors
    /// [`RsError::Parse`] with the line; see the module documentation for what is rejected.
    pub fn parse(text: &str, path: &Path) -> Result<SceneFile, RsError> {
        let pstr = path.display().to_string();
        let pstr = pstr.as_str();
        let mut s = SceneFile {
            path: path.to_path_buf(),
            mesh: String::new(),
            alpha: f64::NAN,
            beta: f64::NAN,
            constraints: vec![],
            potentials: vec![],
            obstacles: vec![],
            implicits: vec![],
            iteration_limit: None,
            time_limit_ms: None,
            allow_barycenter_shift: false,
            disable_near_field: false,
            method: SceneMethod::Hs,
            log: None,
            autotarget_volume: None,
        };
        let mut have_mesh = false;
        let mut seen: BTreeMap<&'static str, usize> = BTreeMap::new();
        for (i, raw) in text.lines().enumerate() {
            let line = i + 1;
            let trimmed = raw.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let mut words = trimmed.split_whitespace();
            let keyword = words.next().unwrap_or("");
            let mut t = Tokens {
                path: pstr,
                line,
                keyword,
                rest: words.collect(),
                next: 0,
            };
            let mut once = |name: &'static str| -> Result<(), RsError> {
                if let Some(prev) = seen.insert(name, line) {
                    return Err(parse_err(
                        pstr,
                        line,
                        format!("`{name}` is given twice (first on line {prev})"),
                    ));
                }
                Ok(())
            };
            match keyword {
                "repel_mesh" => {
                    once("repel_mesh")?;
                    s.mesh = t.word("an OBJ path")?.to_string();
                    if t.remaining() == 0 {
                        return Err(t.err(
                            "`repel_mesh` needs the exponents `alpha beta` (the reference leaves them uninitialized)",
                        ));
                    }
                    s.alpha = t.num("alpha")?;
                    s.beta = t.num("beta")?;
                    Exponents::new(s.alpha, s.beta).map_err(|e| t.err(e.to_string()))?;
                    have_mesh = true;
                }
                "constrain" => {
                    let kind = t.word("a constraint type")?;
                    let c = match kind {
                        "barycenter" => ConstrainStmt::Barycenter,
                        "area" => ConstrainStmt::Area(growth_of(&mut t)?),
                        "volume" => ConstrainStmt::Volume(growth_of(&mut t)?),
                        "vertices" | "vertex" => {
                            let (vertices, motion) = vertex_list(&mut t)?;
                            if vertices.is_empty() {
                                return Err(t.err("`constrain vertices` needs at least one vertex"));
                            }
                            ConstrainStmt::Vertices { vertices, motion }
                        }
                        "boundary_vertices" => ConstrainStmt::BoundaryVertices,
                        "boundary_normals" => ConstrainStmt::BoundaryNormals,
                        "normal" | "normals" => {
                            let (vertices, motion) = vertex_list(&mut t)?;
                            if motion.is_some() {
                                return Err(t.err("`move` applies to vertex pins only"));
                            }
                            ConstrainStmt::Normals(vertices)
                        }
                        other => {
                            return Err(t.err(format!("unknown constraint type `{other}`")));
                        }
                    };
                    t.finish()?;
                    s.constraints.push(At { line, value: c });
                }
                "minimize" => {
                    let k = t.word("a potential")?;
                    let kind = MinimizeKind::from_name(k)
                        .ok_or_else(|| t.err(format!("unknown potential `{k}`")))?;
                    let weight = if t.remaining() > 0 { t.num("weight")? } else { 1.0 };
                    let target = if t.remaining() > 0 {
                        t.num("target value")?
                    } else {
                        0.0
                    };
                    t.finish()?;
                    s.potentials.push(At {
                        line,
                        value: Minimize {
                            kind,
                            weight,
                            target,
                        },
                    });
                }
                "obstacle" | "point_cloud_obstacle" => {
                    let path = t.word("an OBJ path")?.to_string();
                    let weight = if t.remaining() > 0 { t.num("weight")? } else { 1.0 };
                    let mut centered = false;
                    if t.remaining() > 0 {
                        let w = t.word("`centered`")?;
                        if w != "centered" {
                            return Err(t.err(format!("expected `centered`, found `{w}`")));
                        }
                        centered = true;
                    }
                    t.finish()?;
                    s.obstacles.push(At {
                        line,
                        value: ObstacleStmt {
                            path,
                            weight,
                            centered,
                            point_cloud: keyword == "point_cloud_obstacle",
                        },
                    });
                }
                "implicit" => {
                    let kind = t.word("a surface name")?;
                    let action = t.word("`repel` or `attract`")?;
                    let repel = match action {
                        "repel" => true,
                        "attract" => false,
                        o => {
                            return Err(t.err(format!(
                                "the implicit action must be `repel` or `attract`, found `{o}`"
                            )));
                        }
                    };
                    let power = t.num("power")?;
                    let weight = t.num("weight")?;
                    let surface = match kind {
                        "sphere" => ImplicitDesc::Sphere {
                            radius: t.num("radius")?,
                            center: t.vec3("centre")?,
                        },
                        "cylinder" => ImplicitDesc::Cylinder {
                            radius: t.num("radius")?,
                            center: t.vec3("centre")?,
                            axis: t.vec3("axis")?,
                        },
                        "torus" => ImplicitDesc::Torus {
                            major: t.num("major radius")?,
                            minor: t.num("minor radius")?,
                            center: t.vec3("centre")?,
                        },
                        "plane" => ImplicitDesc::Plane {
                            point: t.vec3("point")?,
                            normal: t.vec3("normal")?,
                        },
                        o => return Err(t.err(format!("unrecognized implicit surface `{o}`"))),
                    };
                    t.finish()?;
                    s.implicits.push(At {
                        line,
                        value: ImplicitStmt {
                            surface,
                            repel,
                            power,
                            weight,
                        },
                    });
                }
                "iteration_limit" => {
                    once("iteration_limit")?;
                    let n = t.uint("the iteration limit")?;
                    t.finish()?;
                    s.iteration_limit = (n > 0).then_some(n);
                }
                "time_limit" => {
                    once("time_limit")?;
                    let n = t.uint("the time limit in milliseconds")?;
                    t.finish()?;
                    s.time_limit_ms = (n > 0).then_some(n);
                }
                "allow_barycenter_shift" => {
                    t.finish()?;
                    s.allow_barycenter_shift = true;
                }
                "disable_near_field" => {
                    t.finish()?;
                    s.disable_near_field = true;
                }
                "method" => {
                    once("method")?;
                    let m = t.word("a method name")?;
                    s.method = SceneMethod::from_name(m)
                        .ok_or_else(|| t.err(format!("unknown method `{m}`")))?;
                    t.finish()?;
                }
                "log" => {
                    once("log")?;
                    s.log = Some(t.word("a file name")?.to_string());
                    t.finish()?;
                }
                "autotarget_volume" => {
                    once("autotarget_volume")?;
                    let r = if t.remaining() > 0 { t.num("ratio")? } else { 1.0 };
                    t.finish()?;
                    s.autotarget_volume = Some(r);
                }
                other => {
                    return Err(parse_err(pstr, line, format!("unknown keyword `{other}`")));
                }
            }
        }
        if !have_mesh {
            return Err(parse_err(pstr, 0, "the scene has no `repel_mesh` statement"));
        }
        Ok(s)
    }

    /// Reads and parses a scene file.
    ///
    /// # Errors
    /// [`RsError::Io`]; then as [`SceneFile::parse`].
    pub fn read(path: impl AsRef<Path>) -> Result<SceneFile, RsError> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| RsError::Io {
            path: path.display().to_string(),
            message: e.to_string(),
        })?;
        SceneFile::parse(&text, path)
    }

    /// The directory that relative paths refer to.
    pub fn directory(&self) -> PathBuf {
        match self.path.parent() {
            Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
            _ => PathBuf::from("."),
        }
    }

    /// The path of `name` relative to the scene file.
    pub fn resolve_path(&self, name: &str) -> PathBuf {
        self.directory().join(name)
    }

    /// Every file the scene refers to (mesh, obstacles), resolved.
    pub fn mesh_files(&self) -> Vec<PathBuf> {
        let mut v = vec![self.resolve_path(&self.mesh)];
        v.extend(self.obstacles.iter().map(|o| self.resolve_path(&o.value.path)));
        v
    }

    /// The scene as text in the canonical form; `parse(to_text())` is equal to `self` up to
    /// the line numbers and the path.
    pub fn to_text(&self) -> String {
        let mut o = String::new();
        let n = |x: f64| format!("{x:?}");
        let _ = writeln!(o, "repel_mesh {} {} {}", self.mesh, n(self.alpha), n(self.beta));
        for c in &self.constraints {
            match &c.value {
                ConstrainStmt::Barycenter => o.push_str("constrain barycenter\n"),
                ConstrainStmt::Area(g) | ConstrainStmt::Volume(g) => {
                    let name = if matches!(c.value, ConstrainStmt::Area(_)) {
                        "area"
                    } else {
                        "volume"
                    };
                    match g {
                        None => {
                            let _ = writeln!(o, "constrain {name}");
                        }
                        Some(g) => {
                            let _ = writeln!(
                                o,
                                "constrain {name} {} {} {}",
                                n(g.multiplier),
                                g.iterations,
                                n(g.addition)
                            );
                        }
                    }
                }
                ConstrainStmt::Vertices { vertices, motion } => {
                    let _ = write!(o, "constrain vertices");
                    for v in vertices {
                        let _ = write!(o, " {v}");
                    }
                    if let Some(m) = motion {
                        let _ = write!(
                            o,
                            " move {} {} {} {}",
                            n(m.displacement[0]),
                            n(m.displacement[1]),
                            n(m.displacement[2]),
                            m.iterations
                        );
                    }
                    o.push('\n');
                }
                ConstrainStmt::BoundaryVertices => o.push_str("constrain boundary_vertices\n"),
                ConstrainStmt::BoundaryNormals => o.push_str("constrain boundary_normals\n"),
                ConstrainStmt::Normals(vs) => {
                    let _ = write!(o, "constrain normals");
                    for v in vs {
                        let _ = write!(o, " {v}");
                    }
                    o.push('\n');
                }
            }
        }
        for p in &self.potentials {
            let m = &p.value;
            let _ = writeln!(o, "minimize {} {} {}", m.kind.name(), n(m.weight), n(m.target));
        }
        for ob in &self.obstacles {
            let b = &ob.value;
            let _ = writeln!(
                o,
                "{} {} {}{}",
                if b.point_cloud { "point_cloud_obstacle" } else { "obstacle" },
                b.path,
                n(b.weight),
                if b.centered { " centered" } else { "" }
            );
        }
        for im in &self.implicits {
            let i = &im.value;
            let v3 = |v: &[f64; 3]| format!("{} {} {}", n(v[0]), n(v[1]), n(v[2]));
            let (name, params) = match &i.surface {
                ImplicitDesc::Sphere { radius, center } => {
                    ("sphere", format!("{} {}", n(*radius), v3(center)))
                }
                ImplicitDesc::Cylinder { radius, center, axis } => (
                    "cylinder",
                    format!("{} {} {}", n(*radius), v3(center), v3(axis)),
                ),
                ImplicitDesc::Torus { major, minor, center } => (
                    "torus",
                    format!("{} {} {}", n(*major), n(*minor), v3(center)),
                ),
                ImplicitDesc::Plane { point, normal } => {
                    ("plane", format!("{} {}", v3(point), v3(normal)))
                }
            };
            let _ = writeln!(
                o,
                "implicit {name} {} {} {} {params}",
                if i.repel { "repel" } else { "attract" },
                n(i.power),
                n(i.weight)
            );
        }
        if let Some(k) = self.iteration_limit {
            let _ = writeln!(o, "iteration_limit {k}");
        }
        if let Some(k) = self.time_limit_ms {
            let _ = writeln!(o, "time_limit {k}");
        }
        if self.allow_barycenter_shift {
            o.push_str("allow_barycenter_shift\n");
        }
        if self.disable_near_field {
            o.push_str("disable_near_field\n");
        }
        if self.method != SceneMethod::Hs {
            let _ = writeln!(o, "method {}", self.method.name());
        }
        if let Some(l) = &self.log {
            let _ = writeln!(o, "log {l}");
        }
        if let Some(r) = self.autotarget_volume {
            let _ = writeln!(o, "autotarget_volume {}", n(r));
        }
        o
    }
}

/// What a constraint does over time (`SPEC.md` D11).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScheduleSpec {
    /// Constant.
    None,
    /// `constrain area|volume m n a`: the target grows (D11).
    Growth(Growth),
    /// `constrain vertices … move`: every pin moves by the displacement.
    Move(PinMotion),
}

/// A constraint of the library plus its schedule and origin.
#[derive(Debug, Clone, PartialEq)]
pub struct SceneConstraint {
    /// The line of the statement.
    pub line: usize,
    /// The constraint (targets `None` = initial values).
    pub constraint: Constraint,
    /// Its schedule.
    pub schedule: ScheduleSpec,
}

/// A request to pin the normals of vertices (`constrain boundary_normals`, `constrain normals`).
/// The library has no such constraint: the flow session must report
/// [`RsError::Unsupported`] or implement it.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalPins {
    /// The line of the statement.
    pub line: usize,
    /// The vertices (the boundary vertices for `boundary_normals`).
    pub vertices: Vec<u32>,
}

/// A scene resolved against its files; the description WP5 turns into `FlowOptions`.
#[derive(Clone)]
pub struct SceneSpec {
    /// The syntax this was made from.
    pub file: SceneFile,
    /// The mesh path.
    pub mesh_path: PathBuf,
    /// The mesh topology.
    pub surface: Surface,
    /// The mesh positions.
    pub positions: Vec<[f64; 3]>,
    /// `(α, β)`.
    pub exponents: Exponents,
    /// Hard constraints in file order.
    pub constraints: Vec<SceneConstraint>,
    /// Normal pins (not supported by the library).
    pub normal_pins: Vec<NormalPins>,
    /// Penalties in the reference's order: `minimize …`, obstacles, implicit surfaces.
    pub penalties: Vec<Penalty>,
    /// The line of each entry of `penalties`.
    pub penalty_lines: Vec<usize>,
    /// Sum of the signed volumes of the mesh obstacles (not the point clouds).
    pub obstacle_volume: f64,
    /// `autotarget_volume r`: the final target of the growing volume constraint is
    /// `r · obstacle_volume`.
    pub autotarget_volume: Option<f64>,
    /// `method`.
    pub method: SceneMethod,
    /// `allow_barycenter_shift`.
    pub allow_barycenter_shift: bool,
    /// `disable_near_field` (`low_order: false`, D19).
    pub disable_near_field: bool,
    /// Iteration limit, if any.
    pub iteration_limit: Option<u64>,
    /// Wall-clock limit in milliseconds, if any.
    pub time_limit_ms: Option<u64>,
    /// The CSV log, resolved against the scene directory.
    pub log: Option<PathBuf>,
}

impl std::fmt::Debug for SceneSpec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SceneSpec")
            .field("mesh", &self.mesh_path)
            .field("vertices", &self.surface.n_vertices())
            .field("faces", &self.surface.n_faces())
            .field("exponents", &self.exponents)
            .field("constraints", &self.constraints)
            .field("penalties", &self.penalties)
            .field("method", &self.method)
            .finish()
    }
}

fn signed_volume(vertices: &[[f64; 3]], triangles: impl Iterator<Item = [u32; 3]>) -> f64 {
    let mut v = 0.0;
    for t in triangles {
        let [a, b, c] = t.map(|i| vertices[i as usize]);
        v += dot(a, cross(b, c)) / 6.0;
    }
    v
}

/// The area-weighted barycenter of a triangle soup (the reference's `meshBarycenter`).
fn soup_barycenter(vertices: &[[f64; 3]], triangles: &[[u32; 3]]) -> [f64; 3] {
    let mut c = [0.0; 3];
    let mut w = 0.0;
    for t in triangles {
        let [a, b, d] = t.map(|i| vertices[i as usize]);
        let area = 0.5 * crate::geometry::v3::norm(cross(sub(b, a), sub(d, a)));
        for k in 0..3 {
            c[k] += area * (a[k] + b[k] + d[k]) / 3.0;
        }
        w += area;
    }
    if w > 0.0 { c.map(|x| x / w) } else { [0.0; 3] }
}

fn read_obj_mesh(path: &Path) -> Result<ddg_core::io::ObjMesh, RsError> {
    ddg_core::io::read_obj(path).map_err(|e| match RsError::from(e) {
        RsError::Io { path: _, message } => RsError::Io {
            path: path.display().to_string(),
            message,
        },
        other => other,
    })
}

impl SceneSpec {
    /// Reads a scene file and everything it refers to.
    ///
    /// # Errors
    /// As [`SceneFile::read`] and [`SceneSpec::resolve`].
    pub fn load(path: impl AsRef<Path>) -> Result<SceneSpec, RsError> {
        SceneSpec::resolve(SceneFile::read(path)?)
    }

    /// Loads the mesh and obstacles of a parsed scene and builds the penalties.
    ///
    /// # Errors
    /// [`RsError::Io`] for a missing file (the message names it); [`RsError::Parse`] (with the
    /// scene line) for a vertex index outside the mesh, `boundary_vertices` on a closed mesh,
    /// an invalid implicit surface or a penalty that does not fit the mesh (for instance a
    /// volume on a non-orientable mesh); the mesh errors of [`Surface::new`].
    pub fn resolve(file: SceneFile) -> Result<SceneSpec, RsError> {
        let pstr = file.path.display().to_string();
        let pstr = pstr.as_str();
        let mesh_path = file.resolve_path(&file.mesh);
        let obj = read_obj_mesh(&mesh_path)?;
        let triangles = obj.triangles().ok_or_else(|| {
            let f = (0..obj.n_faces())
                .find(|&f| obj.face_corners(f).len() != 3)
                .unwrap_or(0);
            RsError::PolygonFace {
                face: f as u32,
                corners: obj.face_corners(f).len(),
            }
        })?;
        let n = obj.positions.len();
        let surface = Surface::new(n, triangles)?;
        let positions = obj.positions.clone();
        surface.validate_positions(&positions)?;
        let exponents = Exponents::new(file.alpha, file.beta)?;

        // Attractor vertices: the reference attracts the vertices with a corner UV u > 1e-10
        // when the mesh has any such corner (D22).
        let attracted: Option<Vec<u32>> = {
            let mut flag = vec![false; n];
            for (c, &t) in obj.corner_texcoords.iter().enumerate() {
                if t != ddg_core::topology::NONE && obj.texcoords[t as usize][0] > 1e-10 {
                    flag[obj.corner_verts[c] as usize] = true;
                }
            }
            flag.iter().any(|&f| f).then(|| (0..n as u32).filter(|&v| flag[v as usize]).collect())
        };

        let boundary_vertices = || -> Vec<u32> {
            let mut v: Vec<u32> = surface.boundary_loops().iter().flatten().copied().collect();
            v.sort_unstable();
            v.dedup();
            v
        };

        let mut constraints = Vec::new();
        let mut normal_pins = Vec::new();
        for c in &file.constraints {
            let line = c.line;
            let check = |vs: &[u32]| -> Result<(), RsError> {
                match vs.iter().find(|&&v| v as usize >= n) {
                    Some(v) => Err(parse_err(
                        pstr,
                        line,
                        format!("vertex {v} does not exist (the mesh has {n} vertices)"),
                    )),
                    None => Ok(()),
                }
            };
            let need_boundary = || -> Result<Vec<u32>, RsError> {
                let b = boundary_vertices();
                if b.is_empty() {
                    Err(parse_err(pstr, line, "the mesh has no boundary vertices"))
                } else {
                    Ok(b)
                }
            };
            let sc = |constraint, schedule| SceneConstraint {
                line,
                constraint,
                schedule,
            };
            match &c.value {
                ConstrainStmt::Barycenter => constraints.push(sc(
                    Constraint::Barycenter {
                        component: None,
                        target: None,
                    },
                    ScheduleSpec::None,
                )),
                ConstrainStmt::Area(g) => constraints.push(sc(
                    Constraint::Area {
                        component: None,
                        target: None,
                    },
                    g.map_or(ScheduleSpec::None, ScheduleSpec::Growth),
                )),
                ConstrainStmt::Volume(g) => constraints.push(sc(
                    Constraint::Volume {
                        component: None,
                        target: None,
                    },
                    g.map_or(ScheduleSpec::None, ScheduleSpec::Growth),
                )),
                ConstrainStmt::Vertices { vertices, motion } => {
                    check(vertices)?;
                    constraints.push(sc(
                        Constraint::Pins {
                            vertices: vertices.clone(),
                            targets: None,
                        },
                        motion.map_or(ScheduleSpec::None, ScheduleSpec::Move),
                    ));
                }
                ConstrainStmt::BoundaryVertices => constraints.push(sc(
                    Constraint::Pins {
                        vertices: need_boundary()?,
                        targets: None,
                    },
                    ScheduleSpec::None,
                )),
                ConstrainStmt::BoundaryNormals => normal_pins.push(NormalPins {
                    line,
                    vertices: need_boundary()?,
                }),
                ConstrainStmt::Normals(vs) => {
                    check(vs)?;
                    normal_pins.push(NormalPins {
                        line,
                        vertices: vs.clone(),
                    });
                }
            }
        }

        let (area0, _) = area_and_gradient(&surface, &positions);
        let mut penalties = Vec::new();
        let mut penalty_lines = Vec::new();
        for p in &file.potentials {
            let m = p.value;
            let w = m.weight;
            let pen = match m.kind {
                MinimizeKind::SquaredError => Penalty::SquaredError {
                    weight: w,
                    targets: positions.clone(),
                },
                MinimizeKind::Area => Penalty::TotalArea {
                    weight: w,
                    power: 2.0,
                },
                MinimizeKind::Volume => Penalty::TotalVolume {
                    weight: w,
                    power: 2.0,
                },
                MinimizeKind::AreaDeviation => Penalty::AreaDeviation {
                    weight: w,
                    target: area0,
                },
                MinimizeKind::VolumeDeviation => {
                    let (v0, _) = volume_and_gradient(&surface, &positions);
                    Penalty::VolumeDeviation {
                        weight: w,
                        target: v0,
                    }
                }
                MinimizeKind::BoundaryLength => Penalty::BoundaryLength {
                    weight: w,
                    target: m.target,
                },
                MinimizeKind::BoundaryCurvature => Penalty::BoundaryCurvature {
                    weight: w,
                    dual_length: false,
                },
                MinimizeKind::Willmore => Penalty::Willmore { weight: w },
            };
            pen.validate(&surface)
                .map_err(|e| parse_err(pstr, p.line, e.to_string()))?;
            penalties.push(pen);
            penalty_lines.push(p.line);
        }

        let mut obstacle_volume = 0.0;
        for o in &file.obstacles {
            let path = file.resolve_path(&o.value.path);
            let mut mesh = read_obj_mesh(&path)?;
            if o.value.centered {
                // Triangles of the (fan-triangulated) polygons; the reference recentres the
                // vertices by the area-weighted barycenter.
                let tris = fan_triangles(&mesh);
                let c = soup_barycenter(&mesh.positions, &tris);
                for p in &mut mesh.positions {
                    *p = sub(*p, c);
                }
            }
            let pen = if o.value.point_cloud {
                let set = StaticSet::from_points(&mesh.positions, OBSTACLE_LEAF_SIZE)
                    .map_err(|e| parse_err(pstr, o.line, e.to_string()))?;
                Penalty::PointCloudObstacle {
                    weight: o.value.weight,
                    points: Arc::new(set),
                }
            } else {
                let ob = Obstacle::from_obj_mesh(&mesh)
                    .map_err(|e| parse_err(pstr, o.line, format!("{}: {e}", path.display())))?;
                obstacle_volume += signed_volume(
                    ob.vertices(),
                    ob.triangle_indices().iter().copied(),
                );
                let set = StaticSet::from_obstacle(&ob, OBSTACLE_LEAF_SIZE)
                    .map_err(|e| parse_err(pstr, o.line, e.to_string()))?;
                Penalty::MeshObstacle {
                    weight: o.value.weight,
                    obstacle: Arc::new(set),
                    kind: ObstacleKind::TangentPoint,
                }
            };
            penalties.push(pen);
            penalty_lines.push(o.line);
        }

        for im in &file.implicits {
            let i = &im.value;
            let bad = |e: &dyn std::fmt::Display| parse_err(pstr, im.line, e.to_string());
            let imp: Arc<dyn ImplicitSurface> = match &i.surface {
                ImplicitDesc::Sphere { radius, center } => {
                    Arc::new(Sphere::new(*center, *radius).map_err(|e| bad(&e))?)
                }
                ImplicitDesc::Cylinder {
                    radius,
                    center,
                    axis,
                } => Arc::new(Cylinder::new(*center, *axis, *radius).map_err(|e| bad(&e))?),
                ImplicitDesc::Torus {
                    major,
                    minor,
                    center,
                } => Arc::new(
                    FramedTorus::new(*center, [0.0, 1.0, 0.0], *major, *minor)
                        .map_err(|e| bad(&e))?,
                ),
                ImplicitDesc::Plane { point, normal } => {
                    Arc::new(Plane::new(*point, *normal).map_err(|e| bad(&e))?)
                }
            };
            let pen = if i.repel {
                Penalty::ImplicitObstacle {
                    weight: i.weight,
                    surface: imp,
                    power: i.power,
                    quadrature: Quadrature::PerVertex,
                }
            } else {
                Penalty::ImplicitAttractor {
                    weight: i.weight,
                    surface: imp,
                    power: i.power,
                    quadrature: Quadrature::PerVertex,
                    vertices: attracted.clone(),
                }
            };
            pen.validate(&surface)
                .map_err(|e| bad(&e))?;
            penalties.push(pen);
            penalty_lines.push(im.line);
        }

        let dir = file.directory();
        Ok(SceneSpec {
            mesh_path,
            exponents,
            constraints,
            normal_pins,
            penalties,
            penalty_lines,
            obstacle_volume,
            autotarget_volume: file.autotarget_volume,
            method: file.method,
            allow_barycenter_shift: file.allow_barycenter_shift,
            disable_near_field: file.disable_near_field,
            iteration_limit: file.iteration_limit,
            time_limit_ms: file.time_limit_ms,
            log: file.log.as_ref().map(|l| dir.join(l)),
            surface,
            positions,
            file,
        })
    }
}

fn fan_triangles(mesh: &ddg_core::io::ObjMesh) -> Vec<[u32; 3]> {
    let mut out = Vec::new();
    for f in 0..mesh.n_faces() {
        let c = &mesh.corner_verts[mesh.face_corners(f)];
        for i in 1..c.len().saturating_sub(1) {
            out.push([c[0], c[i], c[i + 1]]);
        }
    }
    out
}
