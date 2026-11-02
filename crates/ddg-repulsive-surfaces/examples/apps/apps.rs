//! The applications of Sec. 9 of "Repulsive Surfaces" (Yu, Brakensiek, Schumacher, Crane 2021),
//! `SPEC.md` §14.4. Shared by `examples/apps/main.rs` (full runs, OBJ sequences) and
//! `tests/applications.rs` (smoke sizes). Every application returns [`Outcome`]s.
#![allow(dead_code, clippy::too_many_arguments)]

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use ddg_repulsive_curves::obstacle::Obstacle;
use ddg_repulsive_surfaces::generate::{Mesh, primitives, shapes};
use ddg_repulsive_surfaces::implicit::Plane;
use ddg_repulsive_surfaces::scene::{SceneFile, SceneSpec};
use ddg_repulsive_surfaces::step::{
    FlowOptions, Schedule, Session, StopRule, scene_options,
};
use ddg_repulsive_surfaces::{
    Constraint, Exponents, ImplicitSurface, Metric, ObstacleKind, Penalty, Quadrature,
    RemeshOptions, RsError, StaticSet, Surface, ccd, read_obj, write_obj,
};

pub type R<T> = Result<T, RsError>;
type P = [f64; 3];

/// Run configuration.
#[derive(Debug, Clone)]
pub struct Cfg {
    /// Reduced sizes and step counts (the tests).
    pub smoke: bool,
    pub threads: usize,
    /// Output directory (OBJ frames, CSV); `None`: nothing is written.
    pub out: Option<PathBuf>,
    /// Frame spacing in steps (0: only the first and last).
    pub frame_every: usize,
    /// Overrides the step count of every flow.
    pub steps: Option<usize>,
    /// The reference scenes (`data/repulsive-surfaces/scenes`).
    pub scenes: PathBuf,
}

impl Cfg {
    pub fn smoke(threads: usize) -> Cfg {
        Cfg { smoke: true, threads, out: None, frame_every: 0, steps: Some(3), scenes: scenes_dir() }
    }
    fn n(&self, full: usize) -> usize {
        self.steps.unwrap_or(full)
    }
    fn faces(&self, full: usize, smoke: usize) -> usize {
        if self.smoke { smoke } else { full }
    }
}

pub fn scenes_dir() -> PathBuf {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    for base in [root.clone(), root.join("../../..")] {
        let d = base.join("data/repulsive-surfaces/scenes");
        if d.is_dir() {
            return d;
        }
    }
    // worktrees: the main checkout may hold the data
    PathBuf::from("/var/data/discrete-geometry/ddg-tools/data/repulsive-surfaces/scenes")
}

/// What an application reports.
#[derive(Debug, Clone, Default)]
pub struct Outcome {
    pub name: String,
    pub faces0: usize,
    pub faces1: usize,
    pub steps: usize,
    pub status: String,
    pub seconds: f64,
    pub e0: f64,
    pub e1: f64,
    pub area0: f64,
    pub area1: f64,
    pub volume0: f64,
    pub volume1: f64,
    pub max_violation: f64,
    pub intersects: Option<bool>,
    pub extra: String,
}

impl Outcome {
    pub fn line(&self) -> String {
        format!(
            "RESULT name={} faces={}->{} steps={} status={} seconds={:.1} energy={:.6e}->{:.6e} area={:.5}->{:.5} volume={:.5}->{:.5} max_violation={:.2e} intersects={} {}",
            self.name, self.faces0, self.faces1, self.steps, self.status, self.seconds, self.e0, self.e1,
            self.area0, self.area1, self.volume0, self.volume1, self.max_violation,
            self.intersects.map_or("n/a".to_string(), |b| b.to_string()), self.extra
        )
    }
}

// ------------------------------------------------------------------------------------------
// Geometry helpers
// ------------------------------------------------------------------------------------------

fn sub(a: P, b: P) -> P { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
fn cross(a: P, b: P) -> P { [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]] }
fn dot(a: P, b: P) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
fn norm(a: P) -> f64 { dot(a, a).sqrt() }

pub fn area_volume(faces: &[[u32; 3]], x: &[P]) -> (f64, f64) {
    let (mut a, mut v) = (0.0, 0.0);
    for f in faces {
        let [p, q, r] = f.map(|i| x[i as usize]);
        a += 0.5 * norm(cross(sub(q, p), sub(r, p)));
        v += dot(p, cross(q, r)) / 6.0;
    }
    (a, v)
}

pub fn vertex_normals(faces: &[[u32; 3]], x: &[P]) -> Vec<P> {
    let mut n = vec![[0.0; 3]; x.len()];
    for f in faces {
        let [p, q, r] = f.map(|i| x[i as usize]);
        let c = cross(sub(q, p), sub(r, p));
        for &i in f {
            for k in 0..3 {
                n[i as usize][k] += c[k];
            }
        }
    }
    n.iter().map(|v| { let l = norm(*v).max(1e-300); [v[0] / l, v[1] / l, v[2] / l] }).collect()
}

/// Eigenvalues (ascending) of a symmetric 3x3 matrix (Jacobi).
pub fn eig3(mut a: [[f64; 3]; 3]) -> [f64; 3] {
    for _ in 0..60 {
        let (mut p, mut q, mut m) = (0, 1, 0.0);
        for i in 0..3 {
            for j in i + 1..3 {
                if a[i][j].abs() > m { m = a[i][j].abs(); p = i; q = j; }
            }
        }
        if m < 1e-14 * (a[0][0].abs() + a[1][1].abs() + a[2][2].abs()).max(1e-300) { break; }
        let th = 0.5 * (2.0 * a[p][q]).atan2(a[q][q] - a[p][p]);
        let (c, s) = (th.cos(), th.sin());
        for k in 0..3 {
            let (akp, akq) = (a[k][p], a[k][q]);
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
        }
        for k in 0..3 {
            let (apk, aqk) = (a[p][k], a[q][k]);
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
        }
    }
    let mut e = [a[0][0], a[1][1], a[2][2]];
    e.sort_by(|x, y| x.partial_cmp(y).unwrap());
    e
}

/// Symmetry indicators of a closed surface: eigenvalues of the area-weighted second-moment
/// tensor about the area centroid (equal eigenvalues: tetrahedral/cubic/icosahedral-type or
/// spherical symmetry; two equal: an axis of order >= 3; none: at most dihedral/mirror), the
/// multiplicity at 3 % relative tolerance, and the ratios.
pub fn inertia_indicators(faces: &[[u32; 3]], x: &[P]) -> ([f64; 3], usize, String) {
    let (mut c, mut aw) = ([0.0; 3], 0.0);
    let mut cen = Vec::new();
    for f in faces {
        let [p, q, r] = f.map(|i| x[i as usize]);
        let a = 0.5 * norm(cross(sub(q, p), sub(r, p)));
        let m = [0, 1, 2].map(|k| (p[k] + q[k] + r[k]) / 3.0);
        for k in 0..3 { c[k] += a * m[k]; }
        aw += a;
        cen.push((a, m));
    }
    for k in 0..3 { c[k] /= aw; }
    let mut t = [[0.0; 3]; 3];
    for (a, m) in &cen {
        let d = sub(*m, c);
        for i in 0..3 { for j in 0..3 { t[i][j] += a * d[i] * d[j] / aw; } }
    }
    let e = eig3(t);
    let tol = 0.03;
    let close = |a: f64, b: f64| (a - b).abs() <= tol * a.abs().max(b.abs());
    let mult = if close(e[0], e[1]) && close(e[1], e[2]) { 3 } else if close(e[0], e[1]) || close(e[1], e[2]) { 2 } else { 1 };
    let desc = format!("moment_eigs=[{:.4},{:.4},{:.4}] ratios={:.3},{:.3} multiplicity={}", e[0], e[1], e[2], e[1] / e[0], e[2] / e[1], mult);
    (e, mult, desc)
}

// ------------------------------------------------------------------------------------------
// Running a flow
// ------------------------------------------------------------------------------------------

pub fn mesh_pair(m: &Mesh) -> R<(Surface, Vec<P>)> {
    Ok((m.surface()?, m.positions.clone()))
}

fn plane_attractor(z: f64, weight: f64, power: f64) -> R<Penalty> {
    let s: Arc<dyn ImplicitSurface> = Arc::new(Plane::new([0.0, 0.0, z], [0.0, 0.0, 1.0])?);
    Ok(Penalty::ImplicitAttractor { weight, surface: s, power, quadrature: Quadrature::DualArea, vertices: None })
}

pub fn mesh_obstacle(m: &Mesh, weight: f64) -> R<(Penalty, Arc<Obstacle>)> {
    let ob = Obstacle::new(m.positions.clone(), m.faces.clone())?;
    let set = StaticSet::from_obstacle(&ob, 8)?;
    Ok((Penalty::MeshObstacle { weight, obstacle: Arc::new(set), kind: ObstacleKind::TangentPoint }, Arc::new(ob)))
}

/// Frames collected by a run: every `frame_every` steps, plus first and last.
pub struct Frames {
    pub positions: Vec<Vec<P>>,
}

/// Runs a flow, writes OBJ frames and a CSV, checks the final surface.
pub fn run(
    cfg: &Cfg,
    name: &str,
    surface: Surface,
    positions: Vec<P>,
    mut opts: FlowOptions,
    keep_frames: bool,
) -> R<(Outcome, Session, Frames)> {
    opts.threads = Some(cfg.threads);
    opts.max_iterations = cfg.n(opts.max_iterations);
    let (a0, v0) = area_volume(surface.faces(), &positions);
    let faces0 = surface.n_faces();
    let io = |p: &Path, e: std::io::Error| RsError::Io { path: p.display().to_string(), message: e.to_string() };
    let dir = cfg.out.as_ref().map(|o| o.join(name));
    let mut csv = None;
    if let Some(d) = &dir {
        std::fs::create_dir_all(d).map_err(|e| io(d, e))?;
        write_obj(d.join("frame_00000.obj"), &surface, &positions)?;
        use std::io::Write;
        let p = d.join("steps.csv");
        let mut f = std::io::BufWriter::new(std::fs::File::create(&p).map_err(|e| io(&p, e))?);
        writeln!(f, "{}", ddg_repulsive_surfaces::step::CSV_HEADER).map_err(|e| io(&p, e))?;
        csv = Some(f);
    }
    let mut frames = Frames { positions: Vec::new() };
    if keep_frames {
        frames.positions.push(positions.clone());
    }
    let t0 = Instant::now();
    let mut session = Session::new(surface, positions, opts)?;
    let total = session.options().max_iterations;
    let every = cfg.frame_every;
    let (mut e0, mut e1, mut viol) = (f64::NAN, f64::NAN, 0.0);
    let mut steps = 0usize;
    let mut status = "MaxIterations".to_string();
    let mut last_written = 0usize;
    while session.iteration() < total {
        match session.step() {
            Ok(r) => {
                if e0.is_nan() { e0 = r.energy_before; }
                e1 = r.energy_after;
                viol = r.max_violation;
                steps += 1;
                if let Some(f) = csv.as_mut() {
                    use std::io::Write;
                    let _ = writeln!(f, "{}", ddg_repulsive_surfaces::step::csv_line(&r));
                    let _ = f.flush();
                }
                let it = session.iteration();
                if every > 0 && it % every == 0 {
                    if let Some(d) = &dir {
                        write_obj(d.join(format!("frame_{it:05}.obj")), session.surface(), session.positions())?;
                        last_written = it;
                    }
                    if keep_frames { frames.positions.push(session.positions().to_vec()); }
                }
                if session.status().is_terminal() {
                    status = format!("{:?}", session.status());
                    break;
                }
            }
            Err(e) => {
                status = format!("Error({e})");
                break;
            }
        }
    }
    if let Some(d) = &dir {
        let it = session.iteration();
        if it != last_written || it == 0 {
            write_obj(d.join(format!("frame_{it:05}.obj")), session.surface(), session.positions())?;
        }
        write_obj(d.join("final.obj"), session.surface(), session.positions())?;
    }
    if keep_frames && every == 0 { frames.positions.push(session.positions().to_vec()); }
    let (a1, v1) = area_volume(session.surface().faces(), session.positions());
    let inter = ccd::intersects(session.surface(), session.positions()).ok();
    let o = Outcome {
        name: name.to_string(),
        faces0,
        faces1: session.surface().n_faces(),
        steps,
        status,
        seconds: t0.elapsed().as_secs_f64(),
        e0,
        e1,
        area0: a0,
        area1: a1,
        volume0: v0,
        volume1: v1,
        max_violation: viol,
        intersects: inter,
        extra: String::new(),
    };
    Ok((o, session, frames))
}

fn base_opts() -> FlowOptions {
    FlowOptions { max_iterations: 1000, stop: StopRule::Off, ..FlowOptions::default() }
}

fn area_constraints() -> Vec<Constraint> {
    vec![Constraint::Barycenter { component: None, target: None }, Constraint::Area { component: None, target: None }]
}

// ------------------------------------------------------------------------------------------
// 1. Fig. 15, 17: genus g
// ------------------------------------------------------------------------------------------

/// Fig. 15 (`plane = None`), Fig. 17 (attracting plane of the given weight).
pub fn genus(cfg: &Cfg, g: u32, plane: Option<f64>) -> R<Outcome> {
    let faces = cfg.faces(3000 + 600 * g as usize, 2500);
    let mesh = shapes::genus_surface(g, faces)?;
    let (s, x) = mesh_pair(&mesh)?;
    let mut o = base_opts();
    o.max_iterations = 800;
    o.constraints = area_constraints();
    if let Some(w) = plane {
        o.penalties.push(plane_attractor(0.0, w, 2.0)?);
    }
    let tag = match plane { None => format!("genus{g:02}"), Some(w) => format!("genus{g:02}_plane{w}") };
    let (mut out, sess, _) = run(cfg, &tag, s, x, o, false)?;
    let (_, _, d) = inertia_indicators(sess.surface().faces(), sess.positions());
    out.extra = format!("genus={} {}", 1 - sess.surface().euler_characteristic() / 2, d);
    Ok(out)
}

// ------------------------------------------------------------------------------------------
// 3. Isotopies
// ------------------------------------------------------------------------------------------

/// Unlinked handcuffs: two rings joined by a bar (same tube radius as the linked ones).
fn unlinked_handcuffs(radius: f64, faces: usize) -> R<Mesh> {
    let n = 160;
    let ring = |cx: f64| -> Vec<P> {
        (0..n).map(|i| { let t = std::f64::consts::TAU * i as f64 / n as f64; [cx + 1.5 * t.cos(), 1.5 * t.sin(), 0.0] }).collect()
    };
    let bar: Vec<P> = vec![[-0.5, 0.0, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
    shapes::graph_tube(&[(ring(-2.0), true), (ring(2.0), true), (bar, false)], radius, faces)
}

/// Fig. 14 / Fig. 1: linked handcuffs and unlinked handcuffs flow (plane attractor, fixed area
/// and barycenter) to the same minimizer. Returns two outcomes; the comparison is in `extra`.
pub fn handcuffs(cfg: &Cfg) -> R<Vec<Outcome>> {
    let faces = cfg.faces(9000, 5000);
    let linked = shapes::handcuffs(1, faces)?;
    let unlinked = unlinked_handcuffs(linked.radius, faces)?;
    // the energy is not scale invariant at (6, 12): same total area for both
    let unlinked = unlinked.clone().transformed((linked.mesh.area() / unlinked.area()).sqrt(), [0.0; 3]);
    let mut res = Vec::new();
    for (tag, mesh) in [("handcuffs_linked", &linked.mesh), ("handcuffs_unlinked", &unlinked)] {
        let (s, x) = mesh_pair(mesh)?;
        let mut o = base_opts();
        o.max_iterations = 1500;
        o.constraints = area_constraints();
        o.penalties.push(plane_attractor(0.0, 5.0, 2.0)?);
        let (mut out, sess, _) = run(cfg, tag, s, x, o, false)?;
        out.extra = format!("genus={}", 1 - sess.surface().euler_characteristic() / 2);
        res.push(out);
    }
    let (e0, e1) = (res[0].e1, res[1].e1);
    let note = format!("energy_linked_final={e0:.5e} energy_unlinked_final={e1:.5e} rel_diff={:.3e}", (e0 - e1).abs() / e1.abs().max(1e-300));
    for r in &mut res { r.extra.push(' '); r.extra.push_str(&note); }
    Ok(res)
}

/// Loads a reference scene with the mesh replaced (the scene's relative paths still resolve
/// against the scene's directory).
fn load_scene(cfg: &Cfg, scene: &str, mesh: Option<&str>) -> R<SceneSpec> {
    let path = cfg.scenes.join(scene);
    let text = std::fs::read_to_string(&path).map_err(|e| RsError::Io { path: path.display().to_string(), message: e.to_string() })?;
    let mut f = SceneFile::parse(&text, &path)?;
    if let Some(m) = mesh { f.mesh = m.to_string(); }
    SceneSpec::resolve(f)
}

/// A flow from a scene: the scene's constraints, penalties, exponents; the surface and positions
/// are the scene's, or given.
fn scene_run(
    cfg: &Cfg,
    name: &str,
    spec: SceneSpec,
    start: Option<(Surface, Vec<P>)>,
    steps: usize,
    tweak: impl FnOnce(&mut FlowOptions),
    keep: bool,
) -> R<(Outcome, Session, Frames)> {
    let mut o = scene_options(&spec)?;
    o.max_iterations = steps;
    o.stop = StopRule::Off;
    tweak(&mut o);
    let (s, x) = start.unwrap_or((spec.surface.clone(), spec.positions.clone()));
    run(cfg, name, s, x, o, keep)
}

/// Fig. 1/18: a handle taken off a pole. Stage 1: the scene (implicit cylinder, plane attractor,
/// fixed area, free barycenter); stage 2: the same mesh with the cylinder weakened as in
/// `endscene.txt`.
pub fn handcuffs_pole(cfg: &Cfg) -> R<Vec<Outcome>> {
    let n1 = if cfg.smoke { 3 } else { 1500 };
    let n2 = if cfg.smoke { 3 } else { 600 };
    let spec = load_scene(cfg, "HandcuffsPole/scene.txt", None)?;
    let (mut o1, s1, _) = scene_run(cfg, "pole_stage1", spec, None, n1, |_| {}, false)?;
    let end = load_scene(cfg, "HandcuffsPole/endscene.txt", None)?;
    let start = (s1.surface().clone(), s1.positions().to_vec());
    let (mut o2, _s2, _) = scene_run(cfg, "pole_stage2", end, Some(start), n2, |_| {}, false)?;
    o1.extra = "cylinder weight 0.1".into();
    o2.extra = "cylinder weight 0.000125 (weakened), from stage 1".into();
    Ok(vec![o1, o2])
}

/// The SlideOffBar variant: a handle on a torus bar (implicit torus), then the torus weakened and
/// a plane attractor added (as `scene_end.txt`; its `end.obj` is not in the reference, so stage 2
/// continues from stage 1's result).
pub fn slide_off_bar(cfg: &Cfg) -> R<Vec<Outcome>> {
    let n1 = if cfg.smoke { 3 } else { 1500 };
    let n2 = if cfg.smoke { 3 } else { 600 };
    let spec = load_scene(cfg, "SlideOffBar/scene_beginning.txt", None)?;
    let (mut o1, s1, _) = scene_run(cfg, "slideoff_stage1", spec, None, n1, |_| {}, false)?;
    let end = load_scene(cfg, "SlideOffBar/scene_end.txt", Some("start.obj"))?;
    let start = (s1.surface().clone(), s1.positions().to_vec());
    let (mut o2, _, _) = scene_run(cfg, "slideoff_stage2", end, Some(start), n2, |_| {}, false)?;
    o1.extra = "torus bar weight 1".into();
    o2.extra = "torus weight 0.001 + plane attractor, from stage 1".into();
    Ok(vec![o1, o2])
}

// ------------------------------------------------------------------------------------------
// 4. Fig. 20: punctured torus eversion
// ------------------------------------------------------------------------------------------

fn transform(x: &[P], refl: usize, quarter: usize) -> Vec<P> {
    x.iter().map(|p| {
        let mut q = *p;
        q[refl] = -q[refl];
        for _ in 0..quarter { q = [-q[1], q[0], q[2]]; }
        q
    }).collect()
}

fn chamfer(a: &[P], b: &[P]) -> f64 {
    let one = |a: &[P], b: &[P]| a.iter().map(|p| b.iter().map(|q| norm(sub(*p, *q))).fold(f64::MAX, f64::min)).sum::<f64>() / a.len() as f64;
    0.5 * (one(a, b) + one(b, a))
}

/// Fig. 20: minimize the tangent-point energy of a punctured torus with the total area fixed and
/// a zero-volume penalty (plus boundary length and curvature penalties); reflect and rotate the
/// midsurface so that the flow can be run in reverse; writes the eversion as an OBJ sequence
/// (forward frames, then the transformed frames reversed with flipped orientation).
pub fn eversion(cfg: &Cfg) -> R<Outcome> {
    let (nm, nn) = if cfg.smoke { (32, 16) } else { (64, 32) };
    let mesh = primitives::punctured_torus(nm, nn, 2.0, 0.8, 0.55)?;
    let (s, x) = mesh_pair(&mesh)?;
    let blen: f64 = s.edges().iter().enumerate().filter(|(e, _)| s.is_boundary_edge(*e as u32)).map(|(_, e)| norm(sub(x[e[0] as usize], x[e[1] as usize]))).sum();
    let mut o = base_opts();
    o.max_iterations = 800;
    o.remesh = None;
    o.constraints = area_constraints();
    o.penalties = vec![
        Penalty::TotalVolume { weight: 1.0e4, power: 2.0 },
        Penalty::BoundaryLength { weight: 30.0, target: blen },
        Penalty::BoundaryCurvature { weight: 1.0, dual_length: true },
    ];
    let (mut out, sess, frames) = run(cfg, "eversion", s, x.clone(), o, true)?;
    let faces = sess.surface().faces().to_vec();
    let mid = sess.positions().to_vec();
    let (a, v) = area_volume(&faces, &mid);
    let (_, v0) = area_volume(&faces, &x);
    // The transformation T (a reflection and a quarter-turn rotation) with T(M) = M as oriented
    // surfaces: small chamfer distance, and the winding-based normals of the (unflipped) image
    // agree with the midsurface's. Candidates whose normals disagree would need a flipped
    // winding: the end of the path would then not be inside-out (reported as `flipped=true`).
    let n_mid = vertex_normals(&faces, &mid);
    let diam = mid.iter().map(|p| norm(*p)).fold(0.0, f64::max);
    let mut best: Option<(f64, usize, usize, f64)> = None;
    let mut best_any: Option<(f64, usize, usize, f64)> = None;
    for refl in 0..3 { for q in 0..4 {
        let tm = transform(&mid, refl, q);
        let c = chamfer(&mid, &tm);
        let n_t = vertex_normals(&faces, &tm);
        let agree: f64 = (0..mid.len()).map(|i| {
            let j = (0..mid.len()).min_by(|&p, &q| norm(sub(tm[p], mid[i])).partial_cmp(&norm(sub(tm[q], mid[i]))).unwrap()).unwrap();
            dot(n_mid[i], n_t[j])
        }).sum::<f64>() / mid.len() as f64;
        let cand = (c, refl, q, agree);
        if best_any.is_none_or(|b| c < b.0) { best_any = Some(cand); }
        if agree > 0.0 && best.is_none_or(|b| c < b.0) { best = Some(cand); }
    } }
    let (cham, brefl, bq, agree) = best.or(best_any).unwrap();
    let best = (cham, brefl, bq);
    let flip = agree < 0.0;
    let last = transform(&x, best.1, best.2);
    let fl: Vec<[u32; 3]> = if flip { faces.iter().map(|f| [f[0], f[2], f[1]]).collect() } else { faces.clone() };
    let (_, v_end) = area_volume(&fl, &last);
    if let Some(d) = cfg.out.as_ref().map(|o| o.join("eversion")) {
        let n1 = frames.positions.len();
        for (k, p) in frames.positions.iter().enumerate() {
            write_obj(d.join(format!("evert_{k:05}.obj")), sess.surface(), p)?;
        }
        for (k, p) in frames.positions.iter().rev().enumerate() {
            let t = transform(p, best.1, best.2);
            let surf = Surface::new(t.len(), fl.clone())?;
            write_obj(d.join(format!("evert_{:05}.obj", n1 + k)), &surf, &t)?;
        }
    }
    out.extra = format!(
        "V0={v0:.4} V_mid={v:.2e} (|V_mid|/|V0|={:.2e}) area_mid={a:.4} boundary_length0={blen:.4} symmetry_T=reflect{}+rot{}x90 chamfer={:.3e} (radius {diam:.3}) normals_agree={agree:.3} flipped={flip} V_end={v_end:.4}",
        (v / v0).abs(), best.1, best.2, best.0
    );
    Ok(out)
}

// ------------------------------------------------------------------------------------------
// 5. Fig. 21-24
// ------------------------------------------------------------------------------------------

type EdgeSet = std::collections::HashSet<(u32, u32)>;

/// 1-to-4 subdivision keeping the old vertices' indices; `sharp` edges (pairs of old vertices)
/// propagate; returns the mesh and the vertices on sharp edges.
fn subdivide(m: &Mesh, levels: u32, sharp: &[(u32, u32)]) -> (Mesh, Vec<u32>) {
    let mut m = m.clone();
    let mut sharp: EdgeSet = sharp.iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
    for _ in 0..levels {
        let mut mid: std::collections::HashMap<(u32, u32), u32> = Default::default();
        let mut pos = m.positions.clone();
        let mut new_sharp: EdgeSet = EdgeSet::new();
        let mut faces = Vec::new();
        for f in &m.faces {
            let mut ids = [0u32; 3];
            for k in 0..3 {
                let (a, b) = (f[k], f[(k + 1) % 3]);
                let key = (a.min(b), a.max(b));
                ids[k] = *mid.entry(key).or_insert_with(|| {
                    let p = [0, 1, 2].map(|c| 0.5 * (pos[a as usize][c] + pos[b as usize][c]));
                    pos.push(p);
                    let i = (pos.len() - 1) as u32;
                    if sharp.contains(&key) {
                        new_sharp.insert((a.min(i), a.max(i)));
                        new_sharp.insert((b.min(i), b.max(i)));
                    }
                    i
                });
            }
            let [a, b, c] = *f;
            let [ab, bc, ca] = ids;
            faces.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
        }
        m = Mesh { positions: pos, faces };
        sharp = new_sharp;
    }
    let mut v: Vec<u32> = sharp.iter().flat_map(|&(a, b)| [a, b]).collect();
    v.sort();
    v.dedup();
    (m, v)
}

/// The prism edges of `h_shape_coarse`: outline edges at both caps and the vertical edges.
fn h_sharp_edges() -> Vec<(u32, u32)> {
    let mut e = Vec::new();
    for i in 0..12u32 {
        let j = (i + 1) % 12;
        e.push((i, j));
        e.push((i + 12, j + 12));
        e.push((i, i + 12));
    }
    e
}

/// Variants of Fig. 22.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HVariant {
    SparsePinsVolume,
    DensePinsVolume,
    NoPinsVolume,
    SparsePinsArea,
    SparsePinsSquaredError,
    SparsePinsVolumeBox,
    DensePinsVolumeBox,
}

impl HVariant {
    pub const ALL: [HVariant; 7] = [
        HVariant::SparsePinsVolume, HVariant::DensePinsVolume, HVariant::NoPinsVolume, HVariant::SparsePinsArea,
        HVariant::SparsePinsSquaredError, HVariant::SparsePinsVolumeBox, HVariant::DensePinsVolumeBox,
    ];
    pub fn name(self) -> &'static str {
        match self {
            HVariant::SparsePinsVolume => "h_sparse_volume",
            HVariant::DensePinsVolume => "h_dense_volume",
            HVariant::NoPinsVolume => "h_nopins_volume",
            HVariant::SparsePinsArea => "h_sparse_area",
            HVariant::SparsePinsSquaredError => "h_sparse_sqerr_volume",
            HVariant::SparsePinsVolumeBox => "h_sparse_volume_box",
            HVariant::DensePinsVolumeBox => "h_dense_volume_box",
        }
    }
}

/// Fig. 22: the H of text with pins; volume growth (x2) or area growth (x1.5); remeshing off.
pub fn h_shape(cfg: &Cfg, v: HVariant) -> R<Outcome> {
    let levels = if cfg.smoke { 1 } else { 3 };
    let (mesh, dense) = subdivide(&primitives::h_shape_coarse(), levels, &h_sharp_edges());
    let (s, x) = mesh_pair(&mesh)?;
    let sparse: Vec<u32> = (0..24).collect();
    let (a0, v0) = area_volume(s.faces(), &x);
    let steps = 300;
    let mut o = base_opts();
    o.max_iterations = steps;
    o.remesh = None;
    if std::env::var("RS_AUTO_BACKEND").is_err() { o.backend = ddg_repulsive_surfaces::Backend::Dense; }
    o.constraints = vec![];
    let pins = |idx: &[u32]| Constraint::Pins { vertices: idx.to_vec(), targets: None };
    use HVariant::*;
    match v {
        SparsePinsVolume | SparsePinsVolumeBox | SparsePinsArea => o.constraints.push(pins(&sparse)),
        DensePinsVolume | DensePinsVolumeBox => o.constraints.push(pins(&dense)),
        SparsePinsSquaredError => {
            o.penalties.push(Penalty::SquaredError { weight: 100.0, targets: x.clone() });
            o.constraints.push(Constraint::Barycenter { component: None, target: None });
        }
        NoPinsVolume => o.constraints.push(Constraint::Barycenter { component: None, target: None }),
    }
    let ci = o.constraints.len();
    match v {
        SparsePinsArea => {
            o.constraints.push(Constraint::Area { component: None, target: None });
            o.schedules.push(Schedule { constraint: ci, factor: 1.5, offset: 0.0, steps });
        }
        _ => {
            o.constraints.push(Constraint::Volume { component: None, target: None });
            o.schedules.push(Schedule { constraint: ci, factor: 2.0, offset: 0.0, steps });
        }
    }
    if matches!(v, SparsePinsVolumeBox | DensePinsVolumeBox) {
        let b = primitives::box_mesh([1.9, 2.4, 0.9], 6)?.transformed(1.0, [1.5, 2.0, 0.5]);
        let (pen, ob) = mesh_obstacle(&b, 0.1)?;
        o.penalties.push(pen);
        o.ccd_obstacles.push(ob);
    }
    let (mut out, sess, _) = run(cfg, v.name(), s, x.clone(), o, false)?;
    let zmax = sess.positions().iter().map(|p| p[2]).fold(f64::MIN, f64::max);
    let zmin = sess.positions().iter().map(|p| p[2]).fold(f64::MAX, f64::min);
    let moved = sparse.iter().map(|&i| norm(sub(sess.positions()[i as usize], x[i as usize]))).fold(0.0, f64::max);
    let npins = if matches!(v, NoPinsVolume) { 0 } else if matches!(v, DensePinsVolume | DensePinsVolumeBox) { dense.len() } else { 24 };
    out.extra = format!("pins={npins} area_factor={:.3} volume_factor={:.3} z_range=[{zmin:.3},{zmax:.3}] max_pin_drift={moved:.2e}", out.area1 / a0, out.volume1 / v0);
    Ok(out)
}

fn read_mesh(path: &Path) -> R<Mesh> {
    let (s, x) = read_obj(path)?;
    Ok(Mesh { positions: x, faces: s.faces().to_vec() })
}

/// Fig. 23: shrink-wrap a point cloud (`scene_bunny` / `scene_volume`: the sphere with the volume
/// constraint falling, a one-sided point-cloud obstacle); the sphere is the 20k-face one.
pub fn shrink_wrap_points(cfg: &Cfg, which: &str) -> R<Outcome> {
    let (scene, tag) = if which == "spot" { ("Reconstruction/scene_volume.txt", "shrinkwrap_spot_points") } else { ("Reconstruction/scene_bunny.txt", "shrinkwrap_bunny_points") };
    let spec = load_scene(cfg, scene, Some("sphere20k.obj"))?;
    let steps = if cfg.smoke { 3 } else { 220 };
    let (mut out, _, _) = scene_run(cfg, tag, spec, None, steps, |_| {}, false)?;
    out.extra = "start=sphere20k".to_string();
    Ok(out)
}

/// Fig. 23: shrink-wrap a polygon soup with holes (the bunny with 35 % of its faces deleted) as
/// a mesh obstacle.
pub fn shrink_wrap_soup(cfg: &Cfg) -> R<Outcome> {
    let spec = load_scene(cfg, "Reconstruction/scene_bunny.txt", Some("sphere20k.obj"))?;
    let bunny = read_mesh(&cfg.scenes.join("Reconstruction/bunny-scaled.obj"))?;
    let soup = primitives::soup_with_holes(&bunny, 0.35, 7)?;
    let verts: Vec<P> = soup.iter().flatten().copied().collect();
    let tris: Vec<[u32; 3]> = (0..soup.len() as u32).map(|i| [3 * i, 3 * i + 1, 3 * i + 2]).collect();
    let ob = Obstacle::new(verts, tris)?;
    let set = StaticSet::from_obstacle(&ob, 8)?;
    let steps = if cfg.smoke { 3 } else { 220 };
    let (mut out, _, _) = scene_run(cfg, "shrinkwrap_bunny_soup", spec, None, steps, |o| {
        o.penalties.retain(|p| !matches!(p, Penalty::PointCloudObstacle { .. }));
        o.penalties.push(Penalty::MeshObstacle { weight: 1.0, obstacle: Arc::new(set), kind: ObstacleKind::TangentPoint });
        o.ccd_obstacles.push(Arc::new(ob));
    }, false)?;
    out.extra = format!("soup: {} of {} faces kept", soup.len(), bunny.faces.len());
    Ok(out)
}

/// Even-odd ray test.
fn inside(faces: &[[u32; 3]], x: &[P], p: P) -> bool {
    let d = [0.5773502691896258, 0.5773502691896258 + 1e-3, 0.5773502691896258 - 1.3e-3];
    let mut hits = 0;
    for f in faces {
        let [a, b, c] = f.map(|i| x[i as usize]);
        let (e1, e2) = (sub(b, a), sub(c, a));
        let h = cross(d, e2);
        let det = dot(e1, h);
        if det.abs() < 1e-14 { continue; }
        let s = sub(p, a);
        let u = dot(s, h) / det;
        let q = cross(s, e1);
        let v = dot(d, q) / det;
        let t = dot(e2, q) / det;
        if u >= 0.0 && v >= 0.0 && u + v <= 1.0 && t > 0.0 { hits += 1; }
    }
    hits % 2 == 1
}

/// Fig. 24: nested envelopes. Level 0 is spot; level k+1 starts from level k pushed outwards by a
/// fraction of an edge, grows to the volume of level k times 1.15, has the previous level as a
/// mesh obstacle (and for the CCD) and is remeshed coarser (`collapse` 1.0, `split` 2.5 of the
/// level's initial mean edge length).
pub fn nested_envelopes(cfg: &Cfg, levels: usize) -> R<Vec<Outcome>> {
    let input = if cfg.smoke { "Walnut/inside.obj" } else { "Reconstruction/spot_triangulated.obj" };
    let (mut surf, mut pos) = read_obj(cfg.scenes.join(input))?;
    let mut res = Vec::new();
    for k in 1..=levels {
        let (_, vprev) = area_volume(surf.faces(), &pos);
        let nrm = vertex_normals(surf.faces(), &pos);
        let mean_e = surf.edges().iter().map(|e| norm(sub(pos[e[0] as usize], pos[e[1] as usize]))).sum::<f64>() / surf.n_edges() as f64;
        let prev_mesh = Mesh { positions: pos.clone(), faces: surf.faces().to_vec() };
        let start: Vec<P> = pos.iter().zip(&nrm).map(|(p, n)| [0, 1, 2].map(|a| p[a] + 0.3 * mean_e * n[a])).collect();
        let (pen, ob) = mesh_obstacle(&prev_mesh, 1.0)?;
        let mut o = base_opts();
        let steps = if cfg.smoke { 3 } else { 150 };
        o.max_iterations = steps;
        o.constraints = vec![Constraint::Barycenter { component: None, target: None }, Constraint::Volume { component: None, target: None }];
        o.schedules.push(Schedule { constraint: 1, factor: 0.0, offset: 1.15 * vprev, steps: (steps * 2 / 3).max(1) });
        o.penalties.push(pen);
        o.ccd_obstacles.push(ob);
        o.remesh = Some(RemeshOptions { collapse: 1.0, split: 2.5, ..RemeshOptions::default() });
        if k >= 2 && std::env::var("RS_AUTO_BACKEND").is_err() { o.backend = ddg_repulsive_surfaces::Backend::Dense; }
        let (mut out, sess, _) = run(cfg, &format!("envelope_level{k}"), surf.clone(), start, o, false)?;
        let (_, v1) = area_volume(sess.surface().faces(), sess.positions());
        let inside_frac = pos.iter().filter(|p| inside(sess.surface().faces(), sess.positions(), **p)).count() as f64 / pos.len() as f64;
        out.extra = format!("volume_ratio_to_previous={:.4} (target 1.15) previous_level_vertices_inside={inside_frac:.4}", v1 / vprev);
        res.push(out);
        surf = sess.surface().clone();
        pos = sess.positions().to_vec();
    }
    Ok(res)
}

// ------------------------------------------------------------------------------------------
// 6. Fig. 25, 26
// ------------------------------------------------------------------------------------------

/// Fig. 25 top: the walnut scene (area x8 inside a shell obstacle, volume deviation penalty).
pub fn walnut(cfg: &Cfg) -> R<Outcome> {
    let spec = load_scene(cfg, "Walnut/scene.txt", None)?;
    let steps = if cfg.smoke { 3 } else { 220 };
    let (mut out, sess, _) = scene_run(cfg, "walnut", spec, None, steps, |_| {}, false)?;
    out.extra = format!("area_ratio={:.3} vertices={}", out.area1 / out.area0, sess.surface().n_vertices());
    Ok(out)
}

/// Fig. 25 bottom: spheres in a slab. Every sphere's area grows while the slab (a box obstacle,
/// also for the CCD) confines them.
pub fn cobblestones(cfg: &Cfg) -> R<Outcome> {
    let (n, level) = if cfg.smoke { (4, 1) } else { (14, 2) };
    let half = if cfg.smoke { [1.6, 1.6, 0.9] } else { [3.2, 3.2, 0.9] };
    let c = primitives::cobblestones(n, 0.5, 0.15, half, level, 11)?;
    let (s, x) = mesh_pair(&c.spheres)?;
    let nc = s.n_components();
    let steps = if cfg.smoke { 3 } else { 400 };
    let growth = 2.5;
    let mut o = base_opts();
    o.max_iterations = steps;
    o.constraints = vec![Constraint::Barycenter { component: None, target: None }];
    for k in 0..nc {
        o.constraints.push(Constraint::Area { component: Some(k as u32), target: None });
        o.schedules.push(Schedule { constraint: 1 + k, factor: growth, offset: 0.0, steps: (steps * 3 / 4).max(1) });
    }
    let slab_w: f64 = std::env::var("RS_SLAB_W").ok().and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let (pen, ob) = mesh_obstacle(&c.slab, slab_w)?;
    if std::env::var("RS_AUTO_BACKEND").is_err() { o.backend = ddg_repulsive_surfaces::Backend::Dense; }
    o.penalties.push(pen);
    o.ccd_obstacles.push(ob);
    let no_remesh = std::env::var("RS_NO_REMESH").is_ok();
    if no_remesh {
        o.remesh = None;
    }
    let tag = if no_remesh { "cobblestones_noremesh".to_string() } else if slab_w == 1.0 { "cobblestones".to_string() } else { format!("cobblestones_w{slab_w}") };
    let (mut out, sess, _) = run(cfg, &tag, s, x, o, false)?;
    let areas: Vec<f64> = sess.surface().components().iter().map(|cm| cm.faces().iter().map(|&f| {
        let t = sess.surface().faces()[f as usize].map(|i| sess.positions()[i as usize]);
        0.5 * norm(cross(sub(t[1], t[0]), sub(t[2], t[0])))
    }).sum()).collect();
    let zmax = sess.positions().iter().map(|p| p[2].abs()).fold(0.0, f64::max);
    let a_each = out.area0 / nc as f64;
    let nout = sess.positions().iter().filter(|p| p[2].abs() > half[2]).count();
    out.extra = format!(
        "components={nc} component_area_ratio_min={:.4} max={:.4} (schedule x{growth}) max|z|={zmax:.3} (slab half {}) vertices_outside_slab={nout}",
        areas.iter().cloned().fold(f64::MAX, f64::min) / a_each,
        areas.iter().cloned().fold(0.0, f64::max) / a_each,
        half[2]
    );
    Ok(out)
}

fn sub_opts(alpha: f64, beta: f64) -> R<FlowOptions> {
    let mut o = base_opts();
    o.exponents = Exponents::new(alpha, beta)?;
    o.metric = Metric::Sobolev { low_order: false };
    o.line_search.ccd = None;
    if std::env::var("RS_AUTO_BACKEND").is_err() { o.backend = ddg_repulsive_surfaces::Backend::Dense; }
    if std::env::var("RS_REMESH").is_err() { o.remesh = None; }
    Ok(o)
}

/// Fig. 26 top: a torus whose tube is slightly too fat (a spindle torus: the tube self-intersects
/// near the axis), flowed with exponents (3.5, 7) and (3, 6), `low_order: false`.
pub fn subcritical_small(cfg: &Cfg) -> R<Vec<Outcome>> {
    let (nm, nn) = if cfg.smoke { (32, 16) } else { (64, 32) };
    let r: f64 = std::env::var("RS_SPINDLE_R").ok().and_then(|s| s.parse().ok()).unwrap_or(1.02);
    let rr = 1.0f64;
    let mut positions = Vec::new();
    for i in 0..nm {
        for j in 0..nn {
            let (u, v) = (std::f64::consts::TAU * i as f64 / nm as f64, std::f64::consts::TAU * j as f64 / nn as f64);
            positions.push([(rr + r * v.cos()) * u.cos(), (rr + r * v.cos()) * u.sin(), r * v.sin()]);
        }
    }
    let id = |i: usize, j: usize| ((i % nm) * nn + (j % nn)) as u32;
    let mut faces = Vec::new();
    for i in 0..nm { for j in 0..nn {
        faces.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
        faces.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
    } }
    let surf = Surface::new(positions.len(), faces)?;
    let pairs0 = ccd::intersecting_pairs(&surf, &positions)?.len();
    let mut res = Vec::new();
    for (tag, a, b) in [("subcritical_spindle_p3.5", 3.5, 7.0), ("subcritical_spindle_p3.0", 3.0, 6.0)] {
        let mut o = sub_opts(a, b)?;
        o.max_iterations = 600;
        o.constraints = area_constraints();
        let tag = if r == 1.02 { tag.to_string() } else { format!("{tag}_r{r}") };
        let (mut out, sess, _) = run(cfg, &tag, surf.clone(), positions.clone(), o, false)?;
        let pairs1 = ccd::intersecting_pairs(sess.surface(), sess.positions())?.len();
        out.extra = format!("exponents=({a},{b}) intersecting_pairs {pairs0}->{pairs1}");
        res.push(out);
    }
    Ok(res)
}

/// Fig. 26 bottom: the Klein bottle, which cannot be embedded: it keeps intersecting.
pub fn klein(cfg: &Cfg) -> R<Outcome> {
    let spec = load_scene(cfg, "KleinBottle/scene.txt", None)?;
    let (s, x) = (spec.surface.clone(), spec.positions.clone());
    let pairs0 = ccd::intersecting_pairs(&s, &x)?.len();
    let steps = if cfg.smoke { 3 } else { 700 };
    let (mut out, sess, _) = scene_run(cfg, "klein_subcritical", spec, None, steps, |o| {
        o.metric = Metric::Sobolev { low_order: false };
        if std::env::var("RS_AUTO_BACKEND").is_err() { o.backend = ddg_repulsive_surfaces::Backend::Dense; }
        o.line_search.ccd = None;
    }, false)?;
    let pairs1 = ccd::intersecting_pairs(sess.surface(), sess.positions())?.len();
    out.extra = format!("exponents=(3.5,7) intersecting_pairs {pairs0}->{pairs1}");
    Ok(out)
}

/// All applications by name (for the CLI): returns the outcomes.
pub fn by_name(cfg: &Cfg, name: &str) -> R<Vec<Outcome>> {
    if let Some(g) = name.strip_prefix("genus") {
        let (g, plane) = match g.split_once("_plane") {
            Some((g, w)) => (g, Some(w.parse::<f64>().unwrap_or(5.0))),
            None => (g, None),
        };
        let g: u32 = g.parse().map_err(|_| RsError::Unsupported { what: "genus<N>[_plane<W>]" })?;
        return Ok(vec![genus(cfg, g, plane)?]);
    }
    for v in HVariant::ALL {
        if v.name() == name { return Ok(vec![h_shape(cfg, v)?]); }
    }
    Ok(match name {
        "handcuffs" => handcuffs(cfg)?,
        "pole" => handcuffs_pole(cfg)?,
        "slideoff" => slide_off_bar(cfg)?,
        "eversion" => vec![eversion(cfg)?],
        "shrinkwrap_bunny" => vec![shrink_wrap_points(cfg, "bunny")?],
        "shrinkwrap_spot" => vec![shrink_wrap_points(cfg, "spot")?],
        "shrinkwrap_soup" => vec![shrink_wrap_soup(cfg)?],
        "envelopes" => nested_envelopes(cfg, if cfg.smoke { 2 } else { 4 })?,
        "walnut" => vec![walnut(cfg)?],
        "cobblestones" => vec![cobblestones(cfg)?],
        "subcritical" => subcritical_small(cfg)?,
        "klein" => vec![klein(cfg)?],
        "control_torus" => vec![control_torus(cfg)?],
        _ => return Err(RsError::Unsupported { what: "unknown application name" }),
    })
}

pub fn all_names() -> Vec<String> {
    let mut v: Vec<String> = (0..=11).map(|g| format!("genus{g:02}")).collect();
    v.extend((2..=6).map(|g| format!("genus{g:02}_plane5")));
    v.extend((2..=6).map(|g| format!("genus{g:02}_plane100")));
    v.extend(HVariant::ALL.iter().map(|h| h.name().to_string()));
    for n in ["handcuffs", "pole", "slideoff", "eversion", "shrinkwrap_bunny", "shrinkwrap_spot", "shrinkwrap_soup", "envelopes", "walnut", "cobblestones", "subcritical", "klein", "control_torus"] {
        v.push(n.to_string());
    }
    v
}

/// Control for Fig. 15, g = 1: a torus of revolution (R/r = 2.5) with the area of `genus01`,
/// flowed with the same options; shows whether the genus-1 run's plateau is a minimizer.
pub fn control_torus(cfg: &Cfg) -> R<Outcome> {
    let m = primitives::torus(64, 32, 1.0555, 0.4222)?;
    let (s, x) = mesh_pair(&m)?;
    let mut o = base_opts();
    o.max_iterations = 800;
    o.constraints = area_constraints();
    let (mut out, sess, _) = run(cfg, "control_torus", s, x, o, false)?;
    let (_, _, d) = inertia_indicators(sess.surface().faces(), sess.positions());
    out.extra = d;
    Ok(out)
}
