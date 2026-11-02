//! Running a reference scene file (`SPEC.md` §14.4): [`SceneSpec`] → [`FlowOptions`] → a run,
//! with a per-step CSV log in the spirit of the reference's `performance.csv`.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Duration;

use ddg_repulsive_curves::obstacle::Obstacle;

use crate::ccd::CcdOptions;
use crate::constraints::Constraint;
use crate::energy::Energy;
use crate::error::RsError;
use crate::metric::{L2Kind, Metric};
use crate::scene::{SceneMethod, SceneSpec, ScheduleSpec};
use crate::surface::Surface;

use super::options::{BcqnVariant, FlowOptions, LbfgsVariant, Method, PinSchedule, Schedule};
use super::session::{Control, FlowReport, Session, StepReport};

/// The [`FlowOptions`] of a scene: exponents, constraints with their growth schedules and pin
/// motions, penalties, the method (`hs` = the Sobolev metric, `disable_near_field` = no `B⁰`;
/// `l2`/`h1`/`h2` the comparison metrics; `willmore` = the Willmore energy with the `H²` metric;
/// `aqp`/`bqn`/`h1_lbfgs` the comparison optimizers of WP9), the barycenter shift, the limits.
/// `autotarget_volume r` makes the growing volume constraint end at `r ×` the obstacles' volume.
/// The mesh obstacles become the CCD's obstacles; the CCD stays on for every exponent (an input that
/// intersects itself, such as the Klein bottle, is handled by the CCD's intersection band).
/// Everything else keeps its default.
///
/// # Errors
/// [`RsError::Unsupported`] for normal constraints (`constrain normals`, `boundary_normals`);
/// [`RsError::Io`] / [`RsError::Curves`] for an obstacle that cannot be read again.
pub fn scene_options(spec: &SceneSpec) -> Result<FlowOptions, RsError> {
    if !spec.normal_pins.is_empty() {
        return Err(RsError::Unsupported {
            what: "normal constraints (`constrain normals`, `constrain boundary_normals`)",
        });
    }
    let mut o = FlowOptions {
        exponents: spec.exponents,
        constraints: Vec::new(),
        penalties: spec.penalties.clone(),
        barycenter_shift: spec.allow_barycenter_shift,
        ..FlowOptions::default()
    };
    for (i, c) in spec.constraints.iter().enumerate() {
        o.constraints.push(c.constraint.clone());
        match &c.schedule {
            ScheduleSpec::None => {}
            ScheduleSpec::Growth(g) => {
                let (factor, offset) = match (&c.constraint, spec.autotarget_volume) {
                    (Constraint::Volume { .. }, Some(r)) => (0.0, r * spec.obstacle_volume),
                    _ => (g.multiplier, g.addition),
                };
                o.schedules.push(Schedule {
                    constraint: i,
                    factor,
                    offset,
                    steps: g.iterations as usize,
                });
            }
            ScheduleSpec::Move(m) => {
                let n = match &c.constraint {
                    Constraint::Pins { vertices, .. } => vertices.len(),
                    _ => 0,
                };
                o.pin_schedules.push(PinSchedule {
                    constraint: i,
                    displacement: vec![m.displacement; n],
                    steps: m.iterations as usize,
                });
            }
        }
    }
    let low_order = !spec.disable_near_field;
    match spec.method {
        SceneMethod::Hs => o.metric = Metric::Sobolev { low_order },
        SceneMethod::L2 => o.metric = Metric::L2(L2Kind::Mass),
        SceneMethod::H1 => o.metric = Metric::H1,
        SceneMethod::H2 => o.metric = Metric::H2,
        SceneMethod::Willmore => {
            o.energy = Energy::Willmore;
            o.metric = Metric::H2;
        }
        SceneMethod::Aqp => o.method = Method::Aqp { kappa: 100.0 },
        SceneMethod::Bqn => o.method = Method::Bcqn { memory: 20, variant: BcqnVariant::Standard },
        SceneMethod::H1Lbfgs => {
            o.method = Method::H1Lbfgs {
                memory: 20,
                variant: LbfgsVariant::Standard,
            }
        }
    }
    for ob in &spec.file.obstacles {
        if !ob.value.point_cloud {
            o.ccd_obstacles.push(Arc::new(ccd_obstacle(spec, &ob.value.path, ob.value.centered)?));
        }
    }
    if let Some(n) = spec.iteration_limit {
        o.max_iterations = n as usize;
    }
    o.time_limit = spec.time_limit_ms.map(Duration::from_millis);
    Ok(o)
}

/// The triangle obstacle of a scene's `obstacle` statement, read again for the CCD (the
/// penalty keeps only centroids and areas), recentred like the scene reader does.
fn ccd_obstacle(spec: &SceneSpec, path: &str, centered: bool) -> Result<Obstacle, RsError> {
    let path = spec.file.resolve_path(path);
    let mut mesh = ddg_core::io::read_obj(&path).map_err(|e| match RsError::from(e) {
        RsError::Io { message, .. } => RsError::Io {
            path: path.display().to_string(),
            message,
        },
        e => e,
    })?;
    if centered {
        let (mut c, mut w) = ([0.0; 3], 0.0);
        for f in 0..mesh.n_faces() {
            let cv = &mesh.corner_verts[mesh.face_corners(f)];
            for i in 1..cv.len().saturating_sub(1) {
                let [a, b, d] = [cv[0], cv[i], cv[i + 1]].map(|v| mesh.positions[v as usize]);
                let (area, _, x) = crate::geometry::face_quantities([a, b, d]);
                for k in 0..3 {
                    c[k] += area * x[k];
                }
                w += area;
            }
        }
        if w > 0.0 {
            for p in &mut mesh.positions {
                for k in 0..3 {
                    p[k] -= c[k] / w;
                }
            }
        }
    }
    Ok(Obstacle::from_obj_mesh(&mesh)?)
}

/// Options of [`run_scene`].
#[derive(Debug, Clone, Default)]
pub struct RunSceneOptions {
    /// Overrides the scene's `iteration_limit`.
    pub max_steps: Option<usize>,
    /// Thread count (`None`: all cores).
    pub threads: Option<usize>,
    /// Turns remeshing off.
    pub no_remesh: bool,
    /// Turns the CCD of the line search on with these options.
    pub ccd: Option<CcdOptions>,
    /// Turns the CCD off (it is on by default; wins over `ccd`).
    pub no_ccd: bool,
    /// Turns the energy-stall stop rule off (run the full step count).
    pub no_stop_rule: bool,
    /// The per-step CSV log (default: the scene's `log`, if any).
    pub csv: Option<PathBuf>,
    /// Writes the final mesh here.
    pub output_obj: Option<PathBuf>,
}

/// The result of [`run_scene`].
#[derive(Debug, Clone)]
pub struct RunReport {
    /// The scene file.
    pub scene: PathBuf,
    /// `|V|`, `|F|` at the start.
    pub initial_size: (usize, usize),
    /// The options used.
    pub options: FlowOptions,
    /// The run.
    pub flow: FlowReport,
    /// The final surface.
    pub surface: Surface,
    /// The final positions.
    pub positions: Vec<[f64; 3]>,
    /// The CSV log written, if any.
    pub csv: Option<PathBuf>,
}

/// The header of the per-step CSV log.
pub const CSV_HEADER: &str = "iteration,status,vertices,faces,energy_before,energy_after,energy_term,\
gradient_norm,sobolev_norm_sq,tau0,tau,tau_ccd,halvings,solver,solves,gmres_iterations,max_residual,\
newton,max_violation,splits,collapses,flips,ms_total,ms_bvh,ms_energy,ms_factor,ms_solve,\
ms_line_search,ms_ccd,ms_projection,ms_remesh,fallbacks,projection_blocked";

/// One CSV line of a step (see [`CSV_HEADER`]).
pub fn csv_line(r: &StepReport) -> String {
    let mut s = String::new();
    let ms = |x: f64| x * 1e3;
    let (sp, co, fl) = r
        .remesh
        .as_ref()
        .map_or((0, 0, 0), |m| (m.splits, m.collapses, m.flips));
    let _ = write!(
        s,
        "{},{:?},{},{},{:.17e},{:.17e},{:.17e},{:.6e},{:.6e},{:.6e},{:.6e},{},{},{},{},{},{:.3e},{},{:.3e},{},{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{},{}",
        r.iteration,
        r.status,
        r.vertices,
        r.faces,
        r.energy_before,
        r.energy_after,
        r.energy_term_after,
        r.gradient_norm,
        r.sobolev_norm_sq,
        r.tau0,
        r.tau,
        r.tau_ccd.map_or(String::new(), |t| format!("{t:.6e}")),
        r.halvings,
        r.linear.solver,
        r.linear.solves.len(),
        r.linear.total_iterations(),
        r.linear.max_relative_residual(),
        r.newton_steps,
        r.max_violation,
        sp,
        co,
        fl,
        ms(r.seconds.total),
        ms(r.seconds.bvh),
        ms(r.seconds.energy),
        ms(r.seconds.factor),
        ms(r.seconds.solve),
        ms(r.seconds.line_search),
        ms(r.seconds.ccd),
        ms(r.seconds.projection),
        ms(r.seconds.remesh),
        r.linear.fallbacks.len(),
        r.projection_blocked,
    );
    s
}

/// Loads a scene file ([`SceneSpec::load`]), builds its options ([`scene_options`]) with the
/// overrides of `run`, runs it, writes the per-step CSV log (flushed after every step) and the
/// final OBJ if asked.
///
/// # Errors
/// The errors of [`SceneSpec::load`], [`scene_options`], [`Session::new`]; [`RsError::Io`] for
/// the log and the OBJ. Errors of steps end the run with `Status::Error` in the report.
pub fn run_scene(path: impl AsRef<Path>, run: &RunSceneOptions) -> Result<RunReport, RsError> {
    use std::io::Write;
    let path = path.as_ref();
    let spec = SceneSpec::load(path)?;
    let mut options = scene_options(&spec)?;
    if let Some(n) = run.max_steps {
        options.max_iterations = n;
    }
    options.threads = run.threads;
    if run.no_remesh {
        options.remesh = None;
    }
    if run.ccd.is_some() {
        options.line_search.ccd = run.ccd;
    }
    if run.no_ccd {
        options.line_search.ccd = None;
    }
    if run.no_stop_rule {
        options.stop = super::options::StopRule::Off;
    }
    let initial_size = (spec.surface.n_vertices(), spec.surface.n_faces());
    let csv = run.csv.clone().or_else(|| spec.log.clone());
    let io_err = |p: &Path, e: std::io::Error| RsError::Io {
        path: p.display().to_string(),
        message: e.to_string(),
    };
    let mut writer = match &csv {
        Some(p) => {
            let f = std::fs::File::create(p).map_err(|e| io_err(p, e))?;
            let mut w = std::io::BufWriter::new(f);
            writeln!(w, "{CSV_HEADER}").map_err(|e| io_err(p, e))?;
            Some(w)
        }
        None => None,
    };
    let mut session = Session::new(spec.surface, spec.positions, options.clone())?;
    let mut write_error: Option<RsError> = None;
    let mut progress = |r: &StepReport| -> Control {
        if let (Some(w), Some(p)) = (writer.as_mut(), csv.as_ref())
            && let Err(e) = writeln!(w, "{}", csv_line(r)).and_then(|()| w.flush())
        {
            write_error = Some(io_err(p, e));
            return Control::Stop;
        }
        Control::Continue
    };
    let flow = session.run(Some(&mut progress))?;
    if let Some(e) = write_error {
        return Err(e);
    }
    if let Some(p) = &run.output_obj {
        crate::io::write_obj(p, session.surface(), session.positions())?;
    }
    Ok(RunReport {
        scene: path.to_path_buf(),
        initial_size,
        options,
        flow,
        surface: session.surface().clone(),
        positions: session.positions().to_vec(),
        csv,
    })
}
