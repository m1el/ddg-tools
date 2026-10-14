//! STUDY: does more resolution for the cut (remeshing, paper Sec. 6.1) improve the convergence
//! and the quality of the cut flow?
//!
//! The flow ("Variational Surface Cutting", Sharp, Crane 2018, crate `ddg-variational-cuts`)
//! runs on an intrinsic triangulation of the input mesh (crate `ddg-intrinsic`), so vertices
//! can be added without changing the surface or the caller's mesh; the cut is carried back to
//! the input mesh. Treatments of the triangulation the flow runs on:
//!
//! ```text
//! input      the caller's triangles (no intrinsic triangulation)
//! flip       Delaunay flips (ddg_uv::cutting::IntrinsicOptions::flips_only, the default), run
//!            through ddg_uv::cutting::Session so that the bend merge of the way back applies
//! u2, u4     uniform refinement: every edge longer than 1/2 (1/4) of the mean input edge is
//!            split at its midpoint, then flips to Delaunay, repeated (no vertex moves)
//! ka2, ka4   curvature-adaptive refinement in the spirit of paper Sec. 6.1: the longest edge of
//!            every face whose integrated curvature |k| A (k: the curvature density of the input
//!            face under the face's centroid, unit-area surface) exceeds a bound is split, then
//!            flips, repeated; the bound is chosen so that as many vertices are added as by u2
//!            (u4)
//! odt2, odt4 the intrinsic optimal Delaunay triangulation with target length 1/2 (1/4) of the
//!            mean input edge (IntrinsicOptions::optimal_delaunay), through the Session
//! suffix S   the flow uses, per face of the triangulation, the curvature density of the input
//!            faces it lies on ("spread"; CurvatureSource::PerFace) instead of the Gaussian
//!            curvature of the triangulation itself, which on an intrinsic triangulation sits
//!            only at the input vertices (an inserted vertex is flat)
//! ```
//!
//! **The yardstick.** Every final cut is a set of straight pieces inside the input faces (the
//! way back of `ddg-intrinsic`; for the Session treatments the description on the caller's
//! triangles after the bend merge, with `merge_tolerance` 0, 0.02 and the default 0.1). It is
//! evaluated on one common reference: the regular 1:n² subdivision of every input face (`n` odd,
//! so that no lattice line follows the dyadic points the refinements insert; the same
//! polyhedron), twice: with the curvature density of the input face each reference face lies in
//! (`ref_*`, "spread": the problem the caller's triangles pose, at a fine resolution), and with
//! the reference's own Gaussian curvature (`poly_*`: the cones of the polyhedron at the input
//! vertices). The level-set function on the reference has the side of the cut (propagated over
//! the lattice by the parity of crossings) times the distance to the pieces, so its zero set
//! crosses every lattice edge where the cut does. The exact length of the cut (sum of the
//! pieces) is printed next to the reference's length term as a check.
//!
//! ```text
//! cargo run --release -p ddg-uv --example remesh_study -- real MESH      one input mesh, all cases
//! cargo run --release -p ddg-uv --example remesh_study -- speed MESH...  timings, run alone
//! cargo run --release -p ddg-uv --example remesh_study -- laws           exact laws (flat circle, sphere latitude)
//! ```
//!
//! `MESH` is a name (`spot`, `bob`, `blub`, `camel`, `spot_control`, `bob_control`,
//! `blub_control`) or an OBJ path. Environment: `REAL_ROOT` (the `data` directory, default the
//! main checkout's), `STUDY_THREADS` (default 9), `REF_FACES` (the largest number of reference
//! faces, default 1 000 000; `n` is the largest odd number up to 33 that fits), `REF_N` (forces
//! `n`), `TREATS` (comma-separated subset), `ALPHAS` (factors of the
//! automatic length weight, default `1,0.5,0.25`), `ENERGIES` (`hencky,dirichlet`), `STARTS`
//! (`0` the plane of the longest axis, `1` of the second longest; default `0`), `CUTS_DIR`
//! (write every final cut as an OBJ of line pieces there).
//!
//! Lines starting with `ROW` are one run each, as `key=value` pairs, for
//! `scripts/remesh_tables.py`.

// Index loops follow the formulas; the result rows are wide tuples.
#![allow(clippy::needless_range_loop, clippy::type_complexity, clippy::too_many_arguments)]

use std::collections::HashMap;
use std::f64::consts::PI;
use std::panic::AssertUnwindSafe;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use ddg_core::generators::{SplitMix64, icosphere};
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, ExportOptions, IntrinsicError, IntrinsicTriangulation,
    OptimalDelaunayOptions, Segment, SimplicialExport, SurfacePoint, TargetLength,
};
use ddg_uv::PolygonMesh;
use ddg_uv::cutting::{
    CutOnPolygons, DEFAULT_CROSSING_OFFSET, DescribeReport, FlowKind, IntrinsicOptions, PathEnd,
    Session, Start, Triangulation,
};
use ddg_variational_cuts::{
    CurvatureSource, Cut, CutProblem, Energy, EvalOptions, Evaluator, Flow, FlowOptions,
    RunReport, StepSize, StopRule, length_weight_for_ripples,
};

// ---------------------------------------------------------------------------------------------
// small helpers

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn env_list(name: &str, default: &str) -> Vec<String> {
    std::env::var(name)
        .unwrap_or_else(|_| default.to_string())
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

/// Runs `f` on every item on `threads` threads; a panic becomes an `Err`.
fn par_map<T: Sync, R: Send>(
    items: &[T],
    threads: usize,
    f: impl Fn(&T) -> R + Sync,
) -> Vec<Result<R, String>> {
    let threads = threads.max(1).min(items.len().max(1));
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<Option<Result<R, String>>>> =
        Mutex::new((0..items.len()).map(|_| None).collect());
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| {
                loop {
                    let k = next.fetch_add(1, Ordering::Relaxed);
                    if k >= items.len() {
                        break;
                    }
                    let r =
                        std::panic::catch_unwind(AssertUnwindSafe(|| f(&items[k]))).map_err(|p| {
                            p.downcast_ref::<String>()
                                .cloned()
                                .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                                .unwrap_or_else(|| "panic".into())
                        });
                    out.lock().unwrap()[k] = Some(r);
                }
            });
        }
    });
    out.into_inner()
        .unwrap()
        .into_iter()
        .map(|x| x.unwrap())
        .collect()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Distance from `p` to the segment `a b` in space.
fn segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(b, a);
    let len2 = dot(d, d);
    let t = if len2 > 0.0 {
        (dot(sub(p, a), d) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    norm(sub(p, [a[0] + t * d[0], a[1] + t * d[1], a[2] + t * d[2]]))
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

// ---------------------------------------------------------------------------------------------
// the input mesh

struct Input {
    name: String,
    positions: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    mesh: TriMesh,
    geometry: Geometry,
    mean_edge: f64,
    /// The length the refinement factors multiply (`TARGET`: `mean` edge, default; `median`
    /// edge; `area`: the side of the equilateral triangle of the mean face area, which does not
    /// depend on how the input's area is cut into slivers).
    base_edge: f64,
    /// `1 / √(total area)`: lengths times this are lengths on the unit-area surface.
    scale: f64,
    /// The Gaussian curvature density of every input face on the unit-area surface (SPEC §3.1
    /// of the cut flow): `CurvatureSource::Gaussian` of the input mesh.
    density: Vec<f64>,
    vertex_faces: Vec<Vec<u32>>,
    /// The planes of the starts: (axis, offset), longest axis first.
    planes: Vec<(usize, f64)>,
}

fn load_obj(path: &str) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>), String> {
    let obj = ddg_core::io::read_obj(path).map_err(|e| format!("read_obj {path}: {e}"))?;
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let r = obj.face_corners(f);
        let c = &obj.corner_verts[r.clone()];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    Ok((obj.positions, faces))
}

/// The planes through the centroid normal to the coordinate axes, by decreasing extent.
fn planes_of(positions: &[[f64; 3]]) -> Vec<(usize, f64)> {
    let n = positions.len() as f64;
    let mut axes: Vec<(usize, f64, f64)> = (0..3)
        .map(|a| {
            let lo = positions.iter().map(|p| p[a]).fold(f64::INFINITY, f64::min);
            let hi = positions
                .iter()
                .map(|p| p[a])
                .fold(f64::NEG_INFINITY, f64::max);
            (a, hi - lo, positions.iter().map(|p| p[a]).sum::<f64>() / n)
        })
        .collect();
    axes.sort_by(|x, y| y.1.total_cmp(&x.1));
    axes.into_iter().map(|(a, _, c)| (a, c)).collect()
}

fn mesh_path(name: &str) -> String {
    let root = std::env::var("REAL_ROOT")
        .unwrap_or_else(|_| "/var/data/discrete-geometry/ddg-tools/data".into());
    let models = format!("{root}/models");
    let cones = format!("{root}/optimal-cones/optimal-cones-data/data");
    match name {
        "spot" => format!("{models}/spot/spot/spot_triangulated.obj"),
        "bob" => format!("{models}/bob/bob/bob_tri.obj"),
        "blub" => format!("{models}/blub/blub/blub_triangulated.obj"),
        "spot_control" => format!("{models}/spot/spot/spot_control_mesh.obj"),
        "bob_control" => format!("{models}/bob/bob/bob_controlmesh.obj"),
        "blub_control" => format!("{models}/blub/blub/blub_control_mesh.obj"),
        n if n.starts_with("camel") => format!("{cones}/{n}.obj"),
        other => other.to_string(),
    }
}

impl Input {
    fn new(name: &str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Result<Input, String> {
        let mesh = TriMesh::new(positions.len(), faces.clone()).map_err(|e| e.to_string())?;
        let geometry = Geometry::from_positions(&mesh, &positions).map_err(|e| e.to_string())?;
        let lengths = geometry.edge_lengths();
        let mean_edge = lengths.iter().sum::<f64>() / lengths.len() as f64;
        let base_edge = match std::env::var("TARGET").as_deref() {
            Ok("median") => median(lengths.to_vec()),
            Ok("area") => {
                (4.0 * geometry.total_area() / (3f64.sqrt() * faces.len() as f64)).sqrt()
            }
            _ => mean_edge,
        };
        let scale = 1.0 / geometry.total_area().sqrt();
        let density = Evaluator::new(&mesh, &geometry)
            .map_err(|e| e.to_string())?
            .gaussian_curvature()
            .to_vec();
        let mut vertex_faces = vec![Vec::new(); positions.len()];
        for (f, face) in faces.iter().enumerate() {
            for &v in face {
                vertex_faces[v as usize].push(f as u32);
            }
        }
        let planes = planes_of(&positions);
        Ok(Input {
            name: name.to_string(),
            positions,
            faces,
            mesh,
            geometry,
            mean_edge,
            base_edge,
            scale,
            density,
            vertex_faces,
            planes,
        })
    }

    fn load(name: &str) -> Result<Input, String> {
        let (positions, faces) = load_obj(&mesh_path(name))?;
        let mut input = Input::new(name, positions, faces)?;
        if name.starts_with("camel") {
            // All camels start from the plane of camel.obj (as in flow_on_bad_meshes.rs).
            if let Ok((p, _)) = load_obj(&mesh_path("camel")) {
                input.planes = planes_of(&p);
            }
        }
        Ok(input)
    }

    fn position(&self, p: SurfacePoint) -> [f64; 3] {
        match p {
            SurfacePoint::Vertex(v) => self.positions[v as usize],
            SurfacePoint::Edge { edge, t } => {
                let [a, b] = self.mesh.edge_vertices(edge as usize);
                let (pa, pb) = (self.positions[a], self.positions[b]);
                [0, 1, 2].map(|k| (1.0 - t) * pa[k] + t * pb[k])
            }
            SurfacePoint::Face { face, bary } => {
                let f = self.faces[face as usize];
                [0, 1, 2].map(|k| (0..3).map(|c| bary[c] * self.positions[f[c] as usize][k]).sum())
            }
        }
    }

    /// The curvature density at a point of the input mesh: that of its face, the mean of the
    /// faces of its edge or of its vertex.
    fn density_at(&self, p: SurfacePoint) -> f64 {
        match p {
            SurfacePoint::Face { face, .. } => self.density[face as usize],
            SurfacePoint::Edge { edge, .. } => {
                let h = self.mesh.edge_halfedge(edge as usize);
                let mut sum = self.density[self.mesh.face(h)];
                let mut n = 1.0;
                if let Some(t) = self.mesh.twin(h) {
                    sum += self.density[self.mesh.face(t)];
                    n += 1.0;
                }
                sum / n
            }
            SurfacePoint::Vertex(v) => {
                let fs = &self.vertex_faces[v as usize];
                fs.iter().map(|&f| self.density[f as usize]).sum::<f64>() / fs.len() as f64
            }
        }
    }

    /// `|k_f| A_f` of every input face on the unit-area surface.
    fn integrated_abs(&self) -> Vec<f64> {
        let s2 = self.scale * self.scale;
        self.geometry
            .face_areas()
            .iter()
            .zip(&self.density)
            .map(|(a, k)| a * s2 * k.abs())
            .collect()
    }

    /// The start function: the signed distance to the plane `which`.
    fn start(&self, which: usize, p: [f64; 3]) -> f64 {
        let (axis, offset) = self.planes[which.min(self.planes.len() - 1)];
        p[axis] - offset
    }
}

// ---------------------------------------------------------------------------------------------
// TREATMENTS

#[derive(Clone, Copy, PartialEq, Debug)]
enum Refine {
    /// The caller's triangles.
    Caller,
    /// Delaunay flips only.
    Flips,
    /// Split edges longer than this multiple of the mean input edge, flip, repeat.
    Uniform(f64),
    /// Curvature-adaptive, with as many added vertices as `Uniform(factor)`.
    Adaptive(f64),
    /// The optimal Delaunay triangulation with this target (multiple of the mean input edge).
    Odt(f64),
}

#[derive(Clone, Copy, Debug)]
struct Treat {
    tag: &'static str,
    refine: Refine,
    /// Curvature density transferred from the input faces.
    spread: bool,
    /// Run through `ddg_uv::cutting::Session` (bend merge measured).
    session: bool,
}

fn all_treats() -> Vec<Treat> {
    let t = |tag, refine, spread, session| Treat {
        tag,
        refine,
        spread,
        session,
    };
    vec![
        t("input", Refine::Caller, false, false),
        t("flip", Refine::Flips, false, true),
        t("flipS", Refine::Flips, true, false),
        t("u2", Refine::Uniform(0.5), false, false),
        t("u2S", Refine::Uniform(0.5), true, false),
        t("u4", Refine::Uniform(0.25), false, false),
        t("u4S", Refine::Uniform(0.25), true, false),
        t("ka2", Refine::Adaptive(0.5), false, false),
        t("ka2S", Refine::Adaptive(0.5), true, false),
        t("ka4", Refine::Adaptive(0.25), false, false),
        t("ka4S", Refine::Adaptive(0.25), true, false),
        t("odt2", Refine::Odt(0.5), false, true),
        t("odt4", Refine::Odt(0.25), false, true),
    ]
}

fn selected_treats() -> Vec<Treat> {
    let all = all_treats();
    match std::env::var("TREATS") {
        Ok(list) => {
            let wanted: Vec<&str> = list.split(',').map(str::trim).collect();
            all.into_iter().filter(|t| wanted.contains(&t.tag)).collect()
        }
        Err(_) => all,
    }
}

/// What a refinement did.
#[derive(Clone, Debug, Default)]
struct BuildInfo {
    added: usize,
    rounds: usize,
    refused: usize,
    seconds: f64,
    /// The bound of the adaptive refinement (integrated curvature per face, unit-area surface).
    bound: f64,
    /// Attempts to match the vertex budget (adaptive).
    attempts: usize,
    /// The uniform refinement stopped at `MAX_VERTICES`.
    capped: bool,
}

/// A triangulation the flow runs on directly (not through the Session).
struct Built {
    tri: Option<IntrinsicTriangulation>,
    export: Option<SimplicialExport>,
    mesh: TriMesh,
    geometry: Geometry,
    /// Positions of the vertices on the input mesh (the polyhedron).
    positions: Vec<[f64; 3]>,
    /// Per face: the curvature density of the input faces it lies on (spread).
    spread: Vec<f64>,
    /// `Σ |k| A` of the spread density against that of the input, minus one.
    spread_mismatch: f64,
    info: BuildInfo,
}

fn flips(tri: &mut IntrinsicTriangulation) -> Result<usize, IntrinsicError> {
    Ok(tri.flip_to_delaunay(&DelaunayOptions::default())?.n_flips)
}

/// Splits the given edges at their midpoints (refusals counted), then flips to Delaunay.
fn split_and_flip(
    tri: &mut IntrinsicTriangulation,
    edges: &[u32],
    info: &mut BuildInfo,
) -> Result<usize, IntrinsicError> {
    let mut n = 0;
    for &e in edges {
        match tri.split_edge(e as usize, 0.5) {
            Ok(_) => n += 1,
            Err(IntrinsicError::InsertionRefused { .. }) => info.refused += 1,
            Err(other) => return Err(other),
        }
    }
    info.added += n;
    flips(tri)?;
    Ok(n)
}

/// Uniform refinement: every edge longer than `h` is split, longest first, then flips; repeated.
fn refine_uniform(
    tri: &mut IntrinsicTriangulation,
    h: f64,
    info: &mut BuildInfo,
) -> Result<(), IntrinsicError> {
    let cap = env_usize("MAX_VERTICES", usize::MAX);
    for _ in 0..80 {
        if tri.n_vertices() >= cap {
            info.capped = true;
            break;
        }
        let mut long: Vec<(f64, u32)> = (0..tri.n_edges())
            .filter(|&e| tri.edge_length(e) > h)
            .map(|e| (tri.edge_length(e), e as u32))
            .collect();
        if long.is_empty() {
            break;
        }
        long.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        info.rounds += 1;
        let room = cap.saturating_sub(tri.n_vertices());
        let edges: Vec<u32> = long.into_iter().map(|x| x.1).take(room).collect();
        if split_and_flip(tri, &edges, info)? == 0 {
            break;
        }
    }
    Ok(())
}

/// The curvature density under the centroid of face `f` of the triangulation.
fn centroid_density(
    tri: &IntrinsicTriangulation,
    input: &Input,
    f: usize,
) -> Result<f64, IntrinsicError> {
    let p = tri.point_on_input(SurfacePoint::Face {
        face: f as u32,
        bary: [1.0 / 3.0; 3],
    })?;
    Ok(input.density_at(p))
}

/// Curvature-adaptive refinement (paper Sec. 6.1): the longest edge of every face whose
/// integrated curvature `|k| A` on the unit-area surface exceeds `bound` is split; flips; repeat.
fn refine_adaptive(
    tri: &mut IntrinsicTriangulation,
    input: &Input,
    bound: f64,
    info: &mut BuildInfo,
) -> Result<(), IntrinsicError> {
    let s2 = input.scale * input.scale;
    let kmax = input.density.iter().fold(0.0f64, |m, k| m.max(k.abs()));
    for _ in 0..80 {
        let mut marked: Vec<(f64, u32)> = Vec::new();
        for f in 0..tri.n_faces() {
            let area = tri.face_area(f) * s2;
            if kmax * area <= bound {
                continue;
            }
            if centroid_density(tri, input, f)?.abs() * area <= bound {
                continue;
            }
            let (mut best, mut best_len) = (0usize, -1.0);
            for k in 0..3 {
                let e = tri.edge(3 * f + k);
                if tri.edge_length(e) > best_len {
                    best_len = tri.edge_length(e);
                    best = e;
                }
            }
            marked.push((best_len, best as u32));
        }
        if marked.is_empty() {
            break;
        }
        marked.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        marked.dedup_by_key(|x| x.1);
        let mut edges: Vec<u32> = marked.into_iter().map(|x| x.1).collect();
        let mut seen = std::collections::HashSet::new();
        edges.retain(|e| seen.insert(*e));
        info.rounds += 1;
        if split_and_flip(tri, &edges, info)? == 0 {
            break;
        }
    }
    Ok(())
}

/// Builds the triangulation of a direct treatment. `budget` is the number of vertices the
/// uniform refinement of the same factor added (for `Adaptive`).
fn build(input: &Input, refine: Refine, budget: Option<usize>) -> Result<Built, String> {
    let start = Instant::now();
    if refine == Refine::Caller {
        let n = input.faces.len();
        return Ok(Built {
            tri: None,
            export: None,
            mesh: input.mesh.clone(),
            geometry: input.geometry.clone(),
            positions: input.positions.clone(),
            spread: input.density.clone(),
            spread_mismatch: 0.0,
            info: BuildInfo {
                seconds: start.elapsed().as_secs_f64(),
                ..BuildInfo::default()
            },
        })
        .inspect(|b| debug_assert_eq!(b.spread.len(), n));
    }
    let err = |e: IntrinsicError| e.to_string();
    let mut tri = IntrinsicTriangulation::new(&input.mesh, &input.geometry).map_err(err)?;
    flips(&mut tri).map_err(err)?;
    let mut info = BuildInfo::default();
    match refine {
        Refine::Caller | Refine::Flips => {}
        Refine::Uniform(factor) => {
            refine_uniform(&mut tri, factor * input.base_edge, &mut info).map_err(err)?;
        }
        Refine::Odt(factor) => {
            let report = tri
                .optimal_delaunay(&OptimalDelaunayOptions {
                    target_length: TargetLength::Absolute(factor * input.base_edge),
                    ..OptimalDelaunayOptions::default()
                })
                .map_err(err)?;
            info.added = report.n_inserted;
            info.rounds = report.n_iterations;
        }
        Refine::Adaptive(_) => {
            let want = budget.ok_or("the adaptive refinement needs a vertex budget")? as f64;
            // First guess: faces split until their share of Σ|k|A is the bound; then a
            // multiplicative secant on the actual count.
            let total: f64 = input.integrated_abs().iter().sum();
            let mut bound = total / (input.faces.len() as f64 + 2.0 * want);
            let base = tri.clone();
            let mut best: Option<(f64, IntrinsicTriangulation, BuildInfo)> = None;
            for attempt in 0..8 {
                let mut t = base.clone();
                let mut i = BuildInfo {
                    bound,
                    attempts: attempt + 1,
                    ..BuildInfo::default()
                };
                refine_adaptive(&mut t, input, bound, &mut i).map_err(err)?;
                let ratio = i.added.max(1) as f64 / want;
                let better = best
                    .as_ref()
                    .is_none_or(|(r, _, _)| (ratio.ln()).abs() < r.ln().abs());
                if better {
                    best = Some((ratio, t, i));
                }
                if (ratio - 1.0).abs() < 0.1 {
                    break;
                }
                // added ≈ C / bound: scale the bound by the ratio (damped).
                bound *= ratio.powf(0.9);
            }
            let (_, t, i) = best.expect("at least one attempt");
            tri = t;
            info = i;
        }
    }
    let export = tri
        .export_simplicial(&ExportOptions::default())
        .map_err(err)?;
    info.added = export.mesh.n_vertices() - input.positions.len();
    info.seconds = start.elapsed().as_secs_f64();
    let positions: Vec<[f64; 3]> = export
        .locations
        .iter()
        .map(|&l| input.position(l))
        .collect();
    // Spread curvature: per face, the mean of the input density at four points.
    let samples = [
        [1.0 / 3.0; 3],
        [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
    ];
    let mut spread = Vec::with_capacity(export.mesh.n_faces());
    for f in 0..export.mesh.n_faces() {
        let mut sum = 0.0;
        for bary in samples {
            let p = tri
                .point_on_input(SurfacePoint::Face {
                    face: f as u32,
                    bary,
                })
                .map_err(err)?;
            sum += input.density_at(p);
        }
        spread.push(sum / samples.len() as f64);
    }
    let s2 = input.scale * input.scale;
    let ours: f64 = spread
        .iter()
        .zip(export.geometry.face_areas())
        .map(|(k, a)| k.abs() * a * s2)
        .sum();
    let theirs: f64 = input.integrated_abs().iter().sum();
    Ok(Built {
        mesh: export.mesh.clone(),
        geometry: export.geometry.clone(),
        positions,
        spread,
        spread_mismatch: ours / theirs - 1.0,
        tri: Some(tri),
        export: Some(export),
        info,
    })
}

// ---------------------------------------------------------------------------------------------
// the cut as pieces on the input faces

/// An end of a piece: barycentric coordinates in the piece's face, and the input edge with the
/// parameter from its lower vertex if the point lies on one.
#[derive(Clone, Copy, Debug)]
struct End {
    bary: [f64; 3],
    edge: Option<(u32, f64)>,
}

#[derive(Clone, Copy, Debug)]
struct Piece {
    face: u32,
    ends: [End; 2],
}

fn end_in_face(input: &Input, face: u32, p: SurfacePoint) -> Result<End, String> {
    let f = input.faces[face as usize];
    let slot = |v: usize| f.iter().position(|&w| w as usize == v);
    match p {
        SurfacePoint::Vertex(v) => {
            let k = slot(v as usize).ok_or("a vertex of a piece is not in its face")?;
            let mut bary = [0.0; 3];
            bary[k] = 1.0;
            Ok(End { bary, edge: None })
        }
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = input.mesh.edge_vertices(edge as usize);
            let (ka, kb) = (
                slot(a).ok_or("an edge point of a piece is not on its face")?,
                slot(b).ok_or("an edge point of a piece is not on its face")?,
            );
            let mut bary = [0.0; 3];
            bary[ka] = 1.0 - t;
            bary[kb] = t;
            Ok(End {
                bary,
                edge: Some((edge, t)),
            })
        }
        SurfacePoint::Face { face: g, bary } => {
            if g != face {
                return Err(format!("a face point of face {g} ends a piece in face {face}"));
            }
            Ok(End { bary, edge: None })
        }
    }
}

/// The cut of `phi` on the caller's triangles (the flow ran on the input mesh).
fn pieces_on_caller(input: &Input, phi: &[f64]) -> Result<Vec<Piece>, String> {
    let cut = Cut::from_level_set(&input.mesh, phi, 1e-3).map_err(|e| e.to_string())?;
    cut.segments
        .iter()
        .map(|s| {
            let e = |k: usize| {
                let p = cut.points[s.points[k] as usize];
                end_in_face(
                    input,
                    s.face,
                    SurfacePoint::Edge {
                        edge: p.edge,
                        t: p.t,
                    },
                )
            };
            Ok(Piece {
                face: s.face,
                ends: [e(0)?, e(1)?],
            })
        })
        .collect()
}

/// The cut of `phi` on an export, carried back to the input mesh (no simplification).
fn pieces_way_back(input: &Input, built: &Built, phi: &[f64]) -> Result<(Vec<Piece>, f64), String> {
    let (tri, export) = (built.tri.as_ref().unwrap(), built.export.as_ref().unwrap());
    let cut = Cut::from_level_set(&export.mesh, phi, 1e-3).map_err(|e| e.to_string())?;
    let nodes: Vec<SurfacePoint> = cut
        .points
        .iter()
        .map(|p| export.edge_point(p.edge as usize, p.t))
        .collect();
    let segments: Vec<Segment> = cut
        .segments
        .iter()
        .map(|s| Segment {
            face: s.face,
            nodes: s.points,
        })
        .collect();
    let back = tri
        .segments_on_input(&nodes, &segments, &EdgeTraceOptions::default())
        .map_err(|e| e.to_string())?;
    let mut pieces = Vec::new();
    for path in &back.paths {
        for (i, &face) in path.faces.iter().enumerate() {
            pieces.push(Piece {
                face,
                ends: [
                    end_in_face(input, face, path.points[i])?,
                    end_in_face(input, face, path.points[i + 1])?,
                ],
            });
        }
    }
    Ok((pieces, back.max_relative_mismatch))
}

/// The cut described on the caller's triangles by the Session (after its bend merge).
fn pieces_described(input: &Input, cut: &CutOnPolygons) -> Result<Vec<Piece>, String> {
    let end = |face: u32, e: PathEnd| -> Result<End, String> {
        match e {
            PathEnd::EdgePoint(i) => {
                let p = &cut.edge_points[i as usize];
                let edge = input
                    .mesh
                    .find_edge(p.edge[0] as usize, p.edge[1] as usize)
                    .ok_or("an edge point is on no edge")?;
                let [a, _] = input.mesh.edge_vertices(edge);
                let t = if a == p.edge[0] as usize { p.t } else { 1.0 - p.t };
                end_in_face(
                    input,
                    face,
                    SurfacePoint::Edge {
                        edge: edge as u32,
                        t,
                    },
                )
            }
            PathEnd::Junction(i) => {
                let j = &cut.junctions[i as usize];
                // corners are polygon corners; for `PolygonMesh::from_triangles` corner c of
                // triangle p is 3 p + k.
                let f = input.faces[face as usize];
                let mut bary = [0.0; 3];
                for k in 0..3 {
                    let v = input.faces[(j.corners[k] / 3) as usize][(j.corners[k] % 3) as usize];
                    let slot = f
                        .iter()
                        .position(|&w| w == v)
                        .ok_or("a junction corner is not in its face")?;
                    bary[slot] = j.barycentric[k];
                }
                Ok(End { bary, edge: None })
            }
        }
    };
    let mut pieces = Vec::new();
    for path in &cut.paths {
        if !path.interior.is_empty() {
            return Err("a path crosses a diagonal; the input has polygons".into());
        }
        pieces.push(Piece {
            face: path.polygon,
            ends: [end(path.polygon, path.start)?, end(path.polygon, path.end)?],
        });
    }
    Ok(pieces)
}

fn piece_points(input: &Input, p: &Piece) -> [[f64; 3]; 2] {
    p.ends.map(|e| {
        input.position(SurfacePoint::Face {
            face: p.face,
            bary: e.bary,
        })
    })
}

/// The exact length of the cut on the unit-area surface.
fn pieces_length(input: &Input, pieces: &[Piece]) -> f64 {
    pieces
        .iter()
        .map(|p| {
            let [a, b] = piece_points(input, p);
            norm(sub(b, a))
        })
        .sum::<f64>()
        * input.scale
}

fn write_cut_obj(input: &Input, pieces: &[Piece], name: &str) {
    let Ok(dir) = std::env::var("CUTS_DIR") else {
        return;
    };
    let _ = std::fs::create_dir_all(&dir);
    let mut s = String::new();
    for p in pieces {
        for q in piece_points(input, p) {
            s.push_str(&format!("v {} {} {}\n", q[0], q[1], q[2]));
        }
    }
    for i in 0..pieces.len() {
        s.push_str(&format!("l {} {}\n", 2 * i + 1, 2 * i + 2));
    }
    let _ = std::fs::write(format!("{dir}/{name}.obj"), s);
}

// ---------------------------------------------------------------------------------------------
// the yardstick: the midpoint subdivision of the input mesh

struct Reference {
    n: usize,
    mesh: TriMesh,
    positions: Vec<[f64; 3]>,
    /// Per reference vertex: an input face it lies in, and its barycentric coordinates there.
    loc: Vec<(u32, [f64; 3])>,
    /// Per input face: the reference vertex of lattice point `(i, j)` at `row_offset(n, j) + i`.
    lattice: Vec<Vec<u32>>,
    /// Per input vertex: its reference vertex.
    of_input_vertex: Vec<u32>,
    /// Per reference face: the input curvature density.
    curvature: Vec<f64>,
    /// Per reference edge (index of `mesh.edges()`), position in `edge_index`.
    edge_index: HashMap<(u32, u32), u32>,
    adjacency: Vec<Vec<(u32, u32)>>,
    evaluators: Vec<Mutex<Evaluator>>,
    next_evaluator: AtomicUsize,
    seconds: f64,
}

/// Offset of row `j` in the triangular lattice of side `n`: `Σ_{j' < j} (n + 1 − j')`.
fn row_offset(n: usize, j: usize) -> usize {
    j * (n + 1) - j * (j.max(1) - 1) / 2
}

impl Reference {
    /// The number of segments per input edge: odd, so that no lattice line of the reference
    /// is a line of dyadic points (midpoints, quarter points) that the refinements insert and
    /// that edges of their triangulations follow; the largest odd `n ≤ 33` whose `F n²` faces
    /// fit `REF_FACES` (`REF_N` forces it).
    fn segments(input: &Input) -> usize {
        if let Ok(k) = std::env::var("REF_N")
            && let Ok(k) = k.parse::<usize>()
        {
            return k.max(1);
        }
        let max_faces = env_usize("REF_FACES", 1_000_000);
        let mut n = 1;
        while n + 2 <= 33 && input.faces.len() * (n + 2) * (n + 2) <= max_faces {
            n += 2;
        }
        n
    }

    fn new(input: &Input, n: usize) -> Result<Reference, String> {
        let start = Instant::now();
        let idx = |i: usize, j: usize| row_offset(n, j) + i;
        let mut positions: Vec<[f64; 3]> = Vec::new();
        let mut loc: Vec<(u32, [f64; 3])> = Vec::new();
        let mut of_input_vertex = vec![u32::MAX; input.positions.len()];
        let mut edge_points: HashMap<(u32, u32, u32), u32> = HashMap::new();
        let mut lattice = Vec::with_capacity(input.faces.len());
        let mut faces: Vec<[u32; 3]> = Vec::new();
        let mut curvature = Vec::new();
        for (fi, f) in input.faces.iter().enumerate() {
            let mut ids = vec![u32::MAX; (n + 1) * (n + 2) / 2];
            for j in 0..=n {
                for i in 0..=(n - j) {
                    let bary = [
                        (n - i - j) as f64 / n as f64,
                        i as f64 / n as f64,
                        j as f64 / n as f64,
                    ];
                    // Which vertex or side of the face the point is on.
                    let corner = if i == 0 && j == 0 {
                        Some(f[0])
                    } else if i == n {
                        Some(f[1])
                    } else if j == n {
                        Some(f[2])
                    } else {
                        None
                    };
                    let side = if corner.is_some() {
                        None
                    } else if j == 0 {
                        Some((f[0], f[1], i))
                    } else if i == 0 {
                        Some((f[0], f[2], j))
                    } else if i + j == n {
                        Some((f[1], f[2], j))
                    } else {
                        None
                    };
                    let mut fresh = || {
                        let id = positions.len() as u32;
                        positions.push(input.position(SurfacePoint::Face {
                            face: fi as u32,
                            bary,
                        }));
                        loc.push((fi as u32, bary));
                        id
                    };
                    let id = if let Some(v) = corner {
                        if of_input_vertex[v as usize] == u32::MAX {
                            of_input_vertex[v as usize] = fresh();
                        }
                        of_input_vertex[v as usize]
                    } else if let Some((a, b, k)) = side {
                        let key = if a < b {
                            (a, b, k as u32)
                        } else {
                            (b, a, (n - k) as u32)
                        };
                        match edge_points.get(&key) {
                            Some(&id) => id,
                            None => {
                                let id = fresh();
                                edge_points.insert(key, id);
                                id
                            }
                        }
                    } else {
                        fresh()
                    };
                    ids[idx(i, j)] = id;
                }
            }
            for j in 0..n {
                for i in 0..(n - j) {
                    faces.push([ids[idx(i, j)], ids[idx(i + 1, j)], ids[idx(i, j + 1)]]);
                    curvature.push(input.density[fi]);
                    if i + j + 1 < n {
                        faces.push([ids[idx(i + 1, j)], ids[idx(i + 1, j + 1)], ids[idx(i, j + 1)]]);
                        curvature.push(input.density[fi]);
                    }
                }
            }
            lattice.push(ids);
        }
        let mesh = TriMesh::new(positions.len(), faces).map_err(|e| e.to_string())?;
        let geometry = Geometry::from_positions(&mesh, &positions).map_err(|e| e.to_string())?;
        let mut edge_index = HashMap::with_capacity(mesh.n_edges());
        let mut adjacency = vec![Vec::new(); positions.len()];
        for (e, &[a, b]) in mesh.edges().iter().enumerate() {
            edge_index.insert((a.min(b), a.max(b)), e as u32);
            adjacency[a as usize].push((b, e as u32));
            adjacency[b as usize].push((a, e as u32));
        }
        let mut evaluators = Vec::new();
        for _ in 0..env_usize("REF_POOL", 2).max(1) {
            evaluators.push(Mutex::new(
                Evaluator::new(&mesh, &geometry).map_err(|e| e.to_string())?,
            ));
        }
        Ok(Reference {
            n,
            mesh,
            positions,
            loc,
            lattice,
            of_input_vertex,
            curvature,
            edge_index,
            adjacency,
            evaluators,
            next_evaluator: AtomicUsize::new(0),
            seconds: start.elapsed().as_secs_f64(),
        })
    }

    fn edge(&self, a: u32, b: u32) -> usize {
        self.edge_index[&(a.min(b), a.max(b))] as usize
    }

    /// The level-set function of the cut on the reference: the side (parity of crossings,
    /// propagated from input vertex 0 with the sign `sign0`) times the distance to the pieces.
    /// Returns it with the number of lattice edges whose parity disagrees with the sides and
    /// the number of input vertices whose side disagrees with `input_sides`.
    fn level_set(
        &self,
        input: &Input,
        pieces: &[Piece],
        input_sides: &[bool],
    ) -> (Vec<f64>, usize, usize) {
        let n = self.n;
        let idx = |i: usize, j: usize| row_offset(n, j) + i;
        let mut parity = vec![false; self.mesh.n_edges()];
        let mut by_face: HashMap<u32, Vec<usize>> = HashMap::new();
        for (k, p) in pieces.iter().enumerate() {
            by_face.entry(p.face).or_default().push(k);
        }
        // Lattice edges inside a face: proper intersections with the pieces of the face, in
        // the affine coordinates (b1, b2) (intersection is affine invariant).
        let orient = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| {
            (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
        };
        let crosses = |a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]| {
            let (o1, o2) = (orient(a, b, c), orient(a, b, d));
            let (o3, o4) = (orient(c, d, a), orient(c, d, b));
            // Half-open: a point on a line counts as on its positive side, the same for every
            // lattice edge and piece that share it, so the parity stays consistent when a piece
            // ends on a lattice edge or passes through a lattice point.
            ((o1 >= 0.0) != (o2 >= 0.0)) && ((o3 >= 0.0) != (o4 >= 0.0))
        };
        for (&face, list) in &by_face {
            let ids = &self.lattice[face as usize];
            let segs: Vec<([f64; 2], [f64; 2])> = list
                .iter()
                .map(|&k| {
                    let e = pieces[k].ends;
                    ([e[0].bary[1], e[0].bary[2]], [e[1].bary[1], e[1].bary[2]])
                })
                .collect();
            let pt = |i: usize, j: usize| [i as f64 / n as f64, j as f64 / n as f64];
            let mut visit = |(i0, j0): (usize, usize), (i1, j1): (usize, usize)| {
                let on_side = (j0 == 0 && j1 == 0)
                    || (i0 == 0 && i1 == 0)
                    || (i0 + j0 == n && i1 + j1 == n);
                if on_side {
                    return;
                }
                let (a, b) = (pt(i0, j0), pt(i1, j1));
                let count = segs.iter().filter(|(c, d)| crosses(a, b, *c, *d)).count();
                if count % 2 == 1 {
                    let e = self.edge(ids[idx(i0, j0)], ids[idx(i1, j1)]);
                    parity[e] ^= true;
                }
            };
            for j in 0..n {
                for i in 0..(n - j) {
                    visit((i, j), (i + 1, j));
                    visit((i, j), (i, j + 1));
                    visit((i + 1, j), (i, j + 1));
                }
            }
        }
        // Lattice edges on an input edge: the crossings of that edge seen from its first face.
        let mut on_edge: HashMap<u32, Vec<f64>> = HashMap::new();
        for p in pieces {
            for e in p.ends {
                if let Some((edge, t)) = e.edge {
                    let first = input.mesh.face(input.mesh.edge_halfedge(edge as usize)) as u32;
                    if first == p.face {
                        on_edge.entry(edge).or_default().push(t);
                    }
                }
            }
        }
        if std::env::var("DEBUG").is_ok() {
            let mut other: HashMap<u32, Vec<f64>> = HashMap::new();
            let mut kinds = [0usize; 3];
            for p in pieces {
                for e in p.ends {
                    kinds[usize::from(e.edge.is_some())] += 1;
                    if e.bary.contains(&1.0) {
                        kinds[2] += 1;
                    }
                    if let Some((edge, t)) = e.edge {
                        let first = input.mesh.face(input.mesh.edge_halfedge(edge as usize)) as u32;
                        if first != p.face {
                            other.entry(edge).or_default().push(t);
                        }
                    }
                }
            }
            let mut bad = 0;
            for (edge, ts) in &other {
                let mut a = on_edge.get(edge).cloned().unwrap_or_default();
                let mut b = ts.clone();
                a.sort_by(f64::total_cmp);
                b.sort_by(f64::total_cmp);
                if a.len() != b.len() || a.iter().zip(&b).any(|(x, y)| (x - y).abs() > 1e-9) {
                    bad += 1;
                    if bad < 5 {
                        eprintln!("edge {edge}: first {a:?} other {b:?}");
                    }
                }
            }
            eprintln!(
                "DEBUG pieces {} ends inside/on edge/at vertex {:?}, edges seen from one side only or differing {bad}",
                pieces.len(),
                kinds
            );
        }
        for (&edge, ts) in &on_edge {
            let [a, b] = input.mesh.edge_vertices(edge as usize);
            // Reference vertices along the edge from a to b.
            let h = input.mesh.edge_halfedge(edge as usize);
            let f = input.mesh.face(h);
            let face = input.faces[f];
            let ids = &self.lattice[f];
            let ka = face.iter().position(|&w| w as usize == a).unwrap();
            let kb = face.iter().position(|&w| w as usize == b).unwrap();
            let lattice_point = |s: usize| {
                // The point at s/n from a towards b, as lattice coordinates (i, j).
                let mut bary = [0usize; 3];
                bary[ka] = n - s;
                bary[kb] = s;
                (bary[1], bary[2])
            };
            for s in 0..n {
                let (lo, hi) = (s as f64 / n as f64, (s + 1) as f64 / n as f64);
                let count = ts.iter().filter(|&&t| t > lo && t <= hi).count();
                if count % 2 == 1 {
                    let (i0, j0) = lattice_point(s);
                    let (i1, j1) = lattice_point(s + 1);
                    let e = self.edge(ids[idx(i0, j0)], ids[idx(i1, j1)]);
                    parity[e] ^= true;
                }
            }
        }
        // Sides by propagation.
        let nv = self.positions.len();
        let mut side = vec![0i8; nv];
        let seed = self.of_input_vertex[0] as usize;
        side[seed] = if input_sides[0] { 1 } else { -1 };
        let mut queue = std::collections::VecDeque::from([seed]);
        while let Some(v) = queue.pop_front() {
            for &(w, e) in &self.adjacency[v] {
                if side[w as usize] == 0 {
                    side[w as usize] = if parity[e as usize] { -side[v] } else { side[v] };
                    queue.push_back(w as usize);
                }
            }
        }
        let conflicts = self
            .mesh
            .edges()
            .iter()
            .enumerate()
            .filter(|&(e, &[a, b])| {
                (side[a as usize] != side[b as usize]) != parity[e] || side[a as usize] == 0
            })
            .count();
        let wrong_inputs = (0..input.positions.len())
            .filter(|&v| (side[self.of_input_vertex[v] as usize] > 0) != input_sides[v])
            .count();
        // Distances to the pieces of the faces around the vertex's face.
        let points: Vec<[[f64; 3]; 2]> = pieces.iter().map(|p| piece_points(input, p)).collect();
        let far = input.mean_edge;
        let floor = 1e-9 * input.mean_edge;
        let phi: Vec<f64> = (0..nv)
            .map(|v| {
                let (f, _) = self.loc[v];
                let x = self.positions[v];
                let mut d = far;
                for &c in &input.faces[f as usize] {
                    for &g in &input.vertex_faces[c as usize] {
                        if let Some(list) = by_face.get(&g) {
                            for &k in list {
                                d = d.min(segment_distance(x, points[k][0], points[k][1]));
                            }
                        }
                    }
                }
                f64::from(side[v]) * d.max(floor)
            })
            .collect();
        (phi, conflicts, wrong_inputs)
    }

    /// The energy of the cut on the reference.
    /// `polyhedral`: the Gaussian curvature of the reference itself (the cones of the
    /// polyhedron at the input vertices, spread over their faces of the reference); otherwise
    /// the curvature density of the input face under every reference face.
    fn energy(
        &self,
        phi: &[f64],
        energy: &Energy,
        alpha: f64,
        polyhedral: bool,
    ) -> Result<ddg_variational_cuts::EnergyTerms, String> {
        let problem = CutProblem {
            curvature: if polyhedral {
                CurvatureSource::Gaussian
            } else {
                CurvatureSource::PerFace(self.curvature.clone())
            },
            ..CutProblem::new(energy.clone(), alpha)
        };
        let k = self.next_evaluator.fetch_add(1, Ordering::Relaxed) % self.evaluators.len();
        let mut ev = self.evaluators[k].lock().unwrap();
        ev.evaluate(phi, &problem, &EvalOptions::default())
            .map(|e| e.energy)
            .map_err(|e| e.to_string())
    }
}

// ---------------------------------------------------------------------------------------------
// runs

#[derive(Clone, Debug)]
struct Case {
    energy: Energy,
    energy_name: &'static str,
    alpha_factor: f64,
    alpha: f64,
    start: usize,
}

impl Case {
    fn tag(&self) -> String {
        format!("{}-a{}-s{}", self.energy_name, self.alpha_factor, self.start)
    }
}

/// What one run gave.
struct Outcome {
    pieces: Vec<Piece>,
    /// Pieces after the bend merge with the tolerances 0.02 and 0.1 (Session runs).
    merged: Vec<(f64, Vec<Piece>, DescribeReport)>,
    sides: Vec<bool>,
    /// The final level-set function at the input vertices (the first vertices of every
    /// triangulation), for the energy of the same cut restricted to the input mesh.
    phi_input: Vec<f64>,
    report: RunReport,
    seconds: f64,
    build_seconds: f64,
    n_vertices: usize,
    added: usize,
    /// Mean edge of the triangulation the flow ran on.
    own_mean_edge: f64,
    mismatch: f64,
    note: String,
}

fn summarize_run(r: &RunReport, input: &Input) -> String {
    let e: Vec<f64> = r
        .steps
        .iter()
        .map(|s| s.energy.total)
        .chain(std::iter::once(r.last_energy))
        .collect();
    let rises = e.windows(2).filter(|w| w[1] > w[0]).count();
    let n = r.steps.len();
    let unit = input.scale * input.mean_edge;
    let tail: Vec<f64> = r.steps[n.saturating_sub(10)..]
        .iter()
        .map(|s| s.cut_motion / unit)
        .collect();
    let pred: Vec<f64> = r.steps[n.saturating_sub(10)..]
        .iter()
        .map(|s| s.predicted_motion / unit)
        .collect();
    format!(
        "steps={} stop={:?} returned={} rises={} motion_end={:.3e} predicted_end={:.3e} e_first={:.6e} e_last={:.6e} e_lowest={:.6e} e_final={:.6e} len_final={:.6e}",
        n,
        r.stop,
        r.returned_step,
        rises,
        median(tail),
        median(pred),
        e.first().copied().unwrap_or(f64::NAN),
        r.last_energy,
        r.lowest_energy,
        r.final_energy.total,
        r.final_energy.length,
    )
}

fn stop_rule() -> StopRule {
    StopRule::default()
}

fn flow_options() -> FlowOptions<'static> {
    FlowOptions::new(StepSize::default())
}

/// A run on a direct treatment.
fn run_direct(input: &Input, built: &Built, treat: &Treat, case: &Case) -> Result<Outcome, String> {
    let phi0: Vec<f64> = built
        .positions
        .iter()
        .map(|&p| input.start(case.start, p))
        .collect();
    let mut flow = Flow::new(&built.mesh, &built.geometry, &phi0).map_err(|e| e.to_string())?;
    flow.redistance().map_err(|e| e.to_string())?;
    let mut problem = CutProblem::new(case.energy.clone(), case.alpha);
    if treat.spread {
        problem.curvature = CurvatureSource::PerFace(built.spread.clone());
    }
    let start = Instant::now();
    let report = flow
        .run(&problem, &flow_options(), &stop_rule())
        .map_err(|e| e.to_string())?;
    let seconds = start.elapsed().as_secs_f64();
    let phi = flow.level_set();
    let phi_input = phi[..input.positions.len()].to_vec();
    let sides: Vec<bool> = phi_input.iter().map(|&x| x > 0.0).collect();
    let (pieces, mismatch) = if built.tri.is_some() {
        pieces_way_back(input, built, &phi)?
    } else {
        (pieces_on_caller(input, &phi)?, 0.0)
    };
    let lengths = built.geometry.edge_lengths();
    Ok(Outcome {
        pieces,
        merged: Vec::new(),
        sides,
        phi_input,
        report,
        seconds,
        build_seconds: built.info.seconds,
        n_vertices: built.mesh.n_vertices(),
        added: built.info.added,
        own_mean_edge: lengths.iter().sum::<f64>() / lengths.len() as f64,
        mismatch,
        note: String::new(),
    })
}

fn session_options(refine: Refine, merge: f64, base_edge: f64) -> IntrinsicOptions {
    let mut o = match refine {
        Refine::Odt(factor) => IntrinsicOptions {
            optimal_delaunay: Some(OptimalDelaunayOptions {
                target_length: TargetLength::Absolute(factor * base_edge),
                ..OptimalDelaunayOptions::default()
            }),
            ..IntrinsicOptions::flips_only()
        },
        _ => IntrinsicOptions::flips_only(),
    };
    o.merge_tolerance = merge;
    o
}

/// The Sessions of one treatment: the one the flow runs in (merge tolerance 0, so its
/// description is the cut as it came back) and one per other merge tolerance, which describe
/// the same state. Built once (the optimal Delaunay triangulation takes seconds) and restarted
/// for every case.
struct SessionSet {
    poly: PolygonMesh,
    main: Session,
    others: Vec<(f64, Session)>,
    build_seconds: f64,
}

impl SessionSet {
    fn new(input: &Input, treat: &Treat) -> Result<SessionSet, String> {
        let poly = PolygonMesh::from_triangles(input.positions.clone(), &input.faces);
        let zero = Start::Corners(vec![vec![0.0; poly.corner_vertices.len()]; 2]);
        let make = |merge: f64| {
            Session::with_triangulation(
                &poly.input(),
                2,
                &zero,
                FlowKind::TwoRegion,
                &Triangulation::Intrinsic(session_options(treat.refine, merge, input.base_edge)),
            )
            .map_err(|e| e.to_string())
        };
        let start = Instant::now();
        let main = make(0.0)?;
        let build_seconds = start.elapsed().as_secs_f64();
        let mut others = Vec::new();
        for merge in [0.02, 0.1] {
            others.push((merge, make(merge)?));
        }
        Ok(SessionSet {
            poly,
            main,
            others,
            build_seconds,
        })
    }

    /// A run through the Session, and the same state described with the other tolerances.
    fn run(&mut self, input: &Input, case: &Case) -> Result<Outcome, String> {
        let corner_values: Vec<f64> = self
            .poly
            .corner_vertices
            .iter()
            .map(|&v| input.start(case.start, input.positions[v as usize]))
            .collect();
        let neg: Vec<f64> = corner_values.iter().map(|x| -x).collect();
        let session = &mut self.main;
        session
            .restart(&Start::Corners(vec![corner_values, neg]))
            .map_err(|e| e.to_string())?;
        session.flow_mut().redistance().map_err(|e| e.to_string())?;
        let problem = session.new_problem(case.energy.clone(), case.alpha);
        let start = Instant::now();
        let report = session
            .run(&problem, &flow_options(), &stop_rule())
            .map_err(|e| e.to_string())?;
        let seconds = start.elapsed().as_secs_f64();
        let state = session.state();
        let mut first_corner = vec![usize::MAX; input.positions.len()];
        for (c, &v) in self.poly.corner_vertices.iter().enumerate() {
            if first_corner[v as usize] == usize::MAX {
                first_corner[v as usize] = c;
            }
        }
        let phi_input: Vec<f64> = first_corner.iter().map(|&c| state[0][c]).collect();
        let sides: Vec<bool> = phi_input.iter().map(|&x| x > 0.0).collect();
        let (cut0, _) = session
            .describe_cut_with_report(DEFAULT_CROSSING_OFFSET)
            .map_err(|e| e.to_string())?;
        let pieces = pieces_described(input, &cut0)?;
        let snapshot = session.snapshot();
        let ir = session.intrinsic_report().cloned();
        let n_vertices = session.engine_mesh().n_vertices();
        let mut merged = Vec::new();
        for (merge, other) in &mut self.others {
            other.restore(&snapshot).map_err(|e| e.to_string())?;
            let (cut, rep) = other
                .describe_cut_with_report(DEFAULT_CROSSING_OFFSET)
                .map_err(|e| e.to_string())?;
            merged.push((*merge, pieces_described(input, &cut)?, rep));
        }
        Ok(Outcome {
            pieces,
            merged,
            sides,
            phi_input,
            report,
            seconds,
            build_seconds: self.build_seconds,
            n_vertices,
            added: ir.as_ref().map_or(0, |r| r.n_added),
            own_mean_edge: f64::NAN,
            mismatch: 0.0,
            note: ir.map_or(String::new(), |r| format!("odt={:?}", r.optimal_delaunay)),
        })
    }
}

// ---------------------------------------------------------------------------------------------
// the study on one input mesh

fn energy_of(name: &str) -> (Energy, &'static str) {
    match name {
        "dirichlet" => (Energy::Dirichlet, "dirichlet"),
        _ => (Energy::Hencky, "hencky"),
    }
}

/// The automatic length weight of the bindings (ripple wavelength 0.25 from the roughening
/// rate of the start state), computed on the flips triangulation for every treatment alike.
fn alpha_auto(input: &Input, flips: &Built, energy: &Energy, start: usize) -> Result<f64, String> {
    let phi0: Vec<f64> = flips
        .positions
        .iter()
        .map(|&p| input.start(start, p))
        .collect();
    let mut flow = Flow::new(&flips.mesh, &flips.geometry, &phi0).map_err(|e| e.to_string())?;
    flow.redistance().map_err(|e| e.to_string())?;
    let problem = CutProblem::new(energy.clone(), 0.0);
    let omega = flow
        .roughening_rate(&problem, &flow_options().eval)
        .map_err(|e| e.to_string())?;
    Ok(length_weight_for_ripples(omega, 0.25))
}

fn cost_rank(t: &Treat) -> u32 {
    match t.refine {
        Refine::Uniform(f) | Refine::Adaptive(f) | Refine::Odt(f) if f < 0.3 => 0,
        Refine::Uniform(_) | Refine::Adaptive(_) | Refine::Odt(_) => 1,
        _ => 2,
    }
}

fn real(name: &str) {
    let threads = env_usize("STUDY_THREADS", 9);
    let input = match Input::load(name) {
        Ok(i) => i,
        Err(e) => {
            println!("{name}: load failed: {e}");
            return;
        }
    };
    let total_k: f64 = input
        .density
        .iter()
        .zip(input.geometry.face_areas())
        .map(|(k, a)| k * a * input.scale * input.scale)
        .sum();
    println!(
        "TARGET {}: base length {:.4e} (mean edge {:.4e}, median {:.4e}), MAX_VERTICES {}",
        std::env::var("TARGET").unwrap_or_else(|_| "mean".into()),
        input.base_edge,
        input.mean_edge,
        median(input.geometry.edge_lengths().to_vec()),
        env_usize("MAX_VERTICES", usize::MAX)
    );
    println!(
        "MESH {}: V {}, F {}, mean edge {:.4e}, total curvature {:.4} (2πχ = {:.4}), Σ|k|A {:.4}",
        input.name,
        input.positions.len(),
        input.faces.len(),
        input.mean_edge,
        total_k,
        2.0 * PI * input.mesh.euler_characteristic() as f64,
        input.integrated_abs().iter().sum::<f64>()
    );
    let n_ref = Reference::segments(&input);
    let reference = match Reference::new(&input, n_ref) {
        Ok(r) => r,
        Err(e) => {
            println!("reference failed: {e}");
            return;
        }
    };
    println!(
        "REFERENCE: 1:{} subdivision ({} segments per input edge), V {}, F {}, built in {:.2} s",
        reference.n * reference.n,
        reference.n,
        reference.positions.len(),
        reference.mesh.n_faces(),
        reference.seconds
    );
    let treats = selected_treats();
    // Direct triangulations: uniform first (their counts are the budgets of the adaptive ones).
    let mut wanted: Vec<Refine> = vec![Refine::Flips];
    for t in &treats {
        if !t.session && !wanted.contains(&t.refine) {
            wanted.push(t.refine);
        }
        if let Refine::Adaptive(f) = t.refine
            && !wanted.contains(&Refine::Uniform(f))
        {
            wanted.push(Refine::Uniform(f));
        }
    }
    let first: Vec<Refine> = wanted
        .iter()
        .copied()
        .filter(|r| !matches!(r, Refine::Adaptive(_)))
        .collect();
    let mut builds: Vec<(Refine, Built)> = Vec::new();
    for (r, b) in first
        .iter()
        .zip(par_map(&first, threads, |&r| build(&input, r, None)))
    {
        match b.and_then(|x| x) {
            Ok(b) => builds.push((*r, b)),
            Err(e) => println!("BUILD {r:?} FAILED: {e}"),
        }
    }
    let second: Vec<(Refine, Option<usize>)> = wanted
        .iter()
        .filter_map(|&r| match r {
            Refine::Adaptive(f) => Some((
                r,
                builds
                    .iter()
                    .find(|(q, _)| *q == Refine::Uniform(f))
                    .map(|(_, b)| b.info.added),
            )),
            _ => None,
        })
        .collect();
    for ((r, _), b) in second.iter().zip(par_map(&second, threads, |&(r, budget)| {
        build(&input, r, budget)
    })) {
        match b.and_then(|x| x) {
            Ok(b) => builds.push((*r, b)),
            Err(e) => println!("BUILD {r:?} FAILED: {e}"),
        }
    }
    for (r, b) in &builds {
        let (lo, hi) = b
            .export
            .as_ref()
            .map_or((f64::NAN, f64::NAN), |e| (e.min_angle, e.max_angle));
        println!(
            "BUILD {:?}: V {} (+{}){}, rounds {}, refused {}, {:.3} s (threads shared), bound {:.3e} after {} attempt(s), export angles {:.1}..{:.1} deg, spread Σ|k|A mismatch {:+.2e}",
            r,
            b.mesh.n_vertices(),
            b.info.added,
            if b.info.capped {
                " CAPPED at MAX_VERTICES"
            } else {
                ""
            },
            b.info.rounds,
            b.info.refused,
            b.info.seconds,
            b.info.bound,
            b.info.attempts,
            lo.to_degrees(),
            hi.to_degrees(),
            b.spread_mismatch
        );
    }
    let Some(flips) = builds
        .iter()
        .find(|(r, _)| *r == Refine::Flips)
        .map(|x| &x.1)
    else {
        println!("no flips triangulation");
        return;
    };
    // Cases.
    let mut cases = Vec::new();
    let starts: Vec<usize> = env_list("STARTS", "0")
        .iter()
        .filter_map(|s| s.parse().ok())
        .collect();
    for e in env_list("ENERGIES", "hencky,dirichlet") {
        let (energy, energy_name) = energy_of(&e);
        for &start in &starts {
            let auto = match alpha_auto(&input, flips, &energy, start) {
                Ok(a) => a,
                Err(err) => {
                    println!("alpha_auto {energy_name} start {start} failed: {err}");
                    continue;
                }
            };
            println!("ALPHA_AUTO {energy_name} start {start}: {auto:.5e}");
            for f in env_list("ALPHAS", "1,0.5,0.25") {
                let factor: f64 = f.parse().unwrap_or(1.0);
                cases.push(Case {
                    energy: energy.clone(),
                    energy_name,
                    alpha_factor: factor,
                    alpha: factor * auto,
                    start,
                });
            }
        }
    }
    // Jobs: one per case and direct treatment; one per Session treatment and energy (its
    // Sessions are built once and restarted for each of its cases).
    let mut jobs: Vec<(Vec<usize>, Treat)> = Vec::new();
    for t in &treats {
        if t.session {
            let mut by_energy: Vec<(&'static str, Vec<usize>)> = Vec::new();
            for (c, case) in cases.iter().enumerate() {
                match by_energy.iter_mut().find(|(e, _)| *e == case.energy_name) {
                    Some((_, list)) => list.push(c),
                    None => by_energy.push((case.energy_name, vec![c])),
                }
            }
            for (_, list) in by_energy {
                jobs.push((list, *t));
            }
        } else {
            for c in 0..cases.len() {
                jobs.push((vec![c], *t));
            }
        }
    }
    jobs.sort_by_key(|(c, t)| (cost_rank(t), !t.session, c[0]));
    let results = par_map(&jobs, threads, |(list, treat)| {
        let mut lines = Vec::new();
        let mut set = if treat.session {
            match SessionSet::new(&input, treat) {
                Ok(s) => Some(s),
                Err(e) => {
                    for &c in list {
                        lines.push(error_row(&input, &cases[c], treat, &e));
                    }
                    return lines;
                }
            }
        } else {
            None
        };
        for &c in list {
            let case = &cases[c];
            let outcome = match &mut set {
                Some(set) => set.run(&input, case),
                None => match builds.iter().find(|(r, _)| *r == treat.refine) {
                    Some((_, b)) => run_direct(&input, b, treat, case),
                    None => Err("not built".to_string()),
                },
            };
            match outcome {
                Ok(o) => lines.extend(emit_rows(&input, &reference, case, treat, &o)),
                Err(e) => lines.push(error_row(&input, case, treat, &e)),
            }
        }
        lines
    });
    let failed = results.iter().filter(|r| r.is_err()).count();
    for (k, r) in results.iter().enumerate() {
        if let Err(e) = r {
            for &c in &jobs[k].0 {
                error_row(&input, &cases[c], &jobs[k].1, &format!("PANIC {e}"));
            }
        }
    }
    println!("DONE {} jobs, {failed} panicked", jobs.len());
}

fn error_row(input: &Input, case: &Case, treat: &Treat, e: &str) -> String {
    let line = format!(
        "ROW mesh={} case={} treat={} variant=cut error={}",
        input.name,
        case.tag(),
        treat.tag,
        e.replace(' ', "_")
    );
    println!("{line}");
    line
}

/// The rows of one run: the cut as it came back, and the merged descriptions.
fn emit_rows(
    input: &Input,
    reference: &Reference,
    case: &Case,
    treat: &Treat,
    o: &Outcome,
) -> Vec<String> {
    let mut lines = Vec::new();
    let run = summarize_run(&o.report, input);
    // The energy of the final level-set function restricted to the input vertices, on the
    // input mesh with its own Gaussian curvature (the measure of flow_on_bad_meshes.rs).
    let on_input = Evaluator::new(&input.mesh, &input.geometry)
        .map_err(|e| e.to_string())
        .and_then(|mut ev| {
            ev.evaluate(
                &o.phi_input,
                &CutProblem::new(case.energy.clone(), case.alpha),
                &EvalOptions::default(),
            )
            .map(|e| e.energy.total)
            .map_err(|e| e.to_string())
        })
        .unwrap_or(f64::NAN);
    let mut variants: Vec<(String, &Vec<Piece>, Option<&DescribeReport>)> =
        vec![("cut".to_string(), &o.pieces, None)];
    for (m, p, r) in &o.merged {
        variants.push((format!("m{m}"), p, Some(r)));
    }
    for (variant, pieces, rep) in variants {
        write_cut_obj(
            input,
            pieces,
            &format!("{}_{}_{}_{}", input.name, case.tag(), treat.tag, variant),
        );
        let start = Instant::now();
        let (phi, conflicts, wrong) = reference.level_set(input, pieces, &o.sides);
        let exact = pieces_length(input, pieces);
        let terms = reference.energy(&phi, &case.energy, case.alpha, false);
        let poly = reference.energy(&phi, &case.energy, case.alpha, true);
        let secs = start.elapsed().as_secs_f64();
        let (total, length, distortion) = match &terms {
            Ok(t) => (t.total, t.length, t.distortion.iter().sum::<f64>()),
            Err(_) => (f64::NAN, f64::NAN, f64::NAN),
        };
        let (poly_total, poly_distortion) = match &poly {
            Ok(t) => (t.total, t.distortion.iter().sum::<f64>()),
            Err(_) => (f64::NAN, f64::NAN),
        };
        let merge = rep.map_or(String::new(), |r| {
            format!(
                " merged={} disp_rel={:.3e} moved_points={}",
                r.n_merged_nodes, r.max_displacement_relative, r.n_moved_points
            )
        });
        let line = format!(
            "ROW mesh={} case={} energy={} alpha_factor={} alpha={:.5e} start={} treat={} variant={} V={} added={} build_s={:.3} run_s={:.3} ms_step={:.2} {} ref_total={:.6e} ref_length={:.6e} ref_distortion={:.6e} exact_length={:.6e} ref_total_exact_len={:.6e} poly_total={:.6e} poly_distortion={:.6e} on_input={:.6e} own_edge={:.4e} pieces={} conflicts={} wrong_sides={} mismatch={:.1e} ref_s={:.2}{}{}{}",
            input.name,
            case.tag(),
            case.energy_name,
            case.alpha_factor,
            case.alpha,
            case.start,
            treat.tag,
            variant,
            o.n_vertices,
            o.added,
            o.build_seconds,
            o.seconds,
            1e3 * o.seconds / o.report.steps.len().max(1) as f64,
            run,
            total,
            length,
            distortion,
            exact,
            distortion + case.alpha * exact,
            poly_total,
            poly_distortion,
            on_input,
            o.own_mean_edge,
            pieces.len(),
            conflicts,
            wrong,
            o.mismatch,
            secs,
            merge,
            if o.note.is_empty() {
                String::new()
            } else {
                format!(" note={}", o.note.replace(' ', ""))
            },
            match (&terms, &poly) {
                (Err(e), _) | (_, Err(e)) => format!(" ref_error={}", e.replace(' ', "_")),
                _ => String::new(),
            }
        );
        println!("{line}");
        lines.push(line);
    }
    lines
}

// ---------------------------------------------------------------------------------------------
// timings, run alone

fn speed(names: &[String]) {
    println!(
        "SPEED (one thread, run alone): build = triangulation and export (Session: the whole session); step = median of steps 5..25, Hencky, automatic alpha, start 0."
    );
    println!(
        "{:<14} {:<6} {:>8} {:>8} {:>9} {:>9}",
        "mesh", "treat", "V", "added", "build s", "step ms"
    );
    for name in names {
        let Ok(input) = Input::load(name) else {
            println!("{name}: load failed");
            continue;
        };
        let Ok(flips) = build(&input, Refine::Flips, None) else {
            continue;
        };
        let alpha = alpha_auto(&input, &flips, &Energy::Hencky, 0).unwrap_or(0.5);
        let mut budgets: HashMap<String, usize> = HashMap::new();
        for treat in selected_treats() {
            if treat.spread && treat.refine != Refine::Flips {
                // Same triangulation and step cost as without S.
                continue;
            }
            let mut times = Vec::new();
            let (v, added, build_s);
            if treat.session {
                let poly = PolygonMesh::from_triangles(input.positions.clone(), &input.faces);
                let values: Vec<f64> = poly
                    .corner_vertices
                    .iter()
                    .map(|&v| input.start(0, input.positions[v as usize]))
                    .collect();
                let neg: Vec<f64> = values.iter().map(|x| -x).collect();
                let start = Instant::now();
                let Ok(mut session) = Session::with_triangulation(
                    &poly.input(),
                    2,
                    &Start::Corners(vec![values, neg]),
                    FlowKind::TwoRegion,
                    &Triangulation::Intrinsic(session_options(treat.refine, 0.1, input.base_edge)),
                ) else {
                    continue;
                };
                build_s = start.elapsed().as_secs_f64();
                let _ = session.flow_mut().redistance();
                let problem = session.new_problem(Energy::Hencky, alpha);
                for _ in 0..25 {
                    let s = Instant::now();
                    if session.step(&problem, &flow_options()).is_err() {
                        break;
                    }
                    times.push(s.elapsed().as_secs_f64());
                }
                v = session.engine_mesh().n_vertices();
                added = session.intrinsic_report().map_or(0, |r| r.n_added);
            } else {
                let budget = match treat.refine {
                    Refine::Adaptive(f) => budgets.get(&format!("{f}")).copied(),
                    _ => None,
                };
                let Ok(b) = build(&input, treat.refine, budget) else {
                    continue;
                };
                if let Refine::Uniform(f) = treat.refine {
                    budgets.insert(format!("{f}"), b.info.added);
                }
                build_s = b.info.seconds;
                let phi0: Vec<f64> = b.positions.iter().map(|&p| input.start(0, p)).collect();
                let Ok(mut flow) = Flow::new(&b.mesh, &b.geometry, &phi0) else {
                    continue;
                };
                let _ = flow.redistance();
                let mut problem = CutProblem::new(Energy::Hencky, alpha);
                if treat.spread {
                    problem.curvature = CurvatureSource::PerFace(b.spread.clone());
                }
                for _ in 0..25 {
                    let s = Instant::now();
                    if flow.step(&problem, &flow_options()).is_err() {
                        break;
                    }
                    times.push(s.elapsed().as_secs_f64());
                }
                v = b.mesh.n_vertices();
                added = b.info.added;
            }
            let step = median(times.iter().skip(5).copied().collect());
            println!(
                "{:<14} {:<6} {:>8} {:>8} {:>9.3} {:>9.2}",
                input.name,
                treat.tag,
                v,
                added,
                build_s,
                1e3 * step
            );
        }
    }
}

// ---------------------------------------------------------------------------------------------
// exact laws

fn polar(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2])
}

/// Positions of the cut points of `phi` on a built triangulation.
fn cut_positions(b: &Built, phi: &[f64]) -> Vec<[f64; 3]> {
    match Cut::from_level_set(&b.mesh, phi, 1e-3) {
        Ok(cut) => cut
            .point_positions(&b.mesh, &b.positions)
            .unwrap_or_default(),
        Err(_) => Vec::new(),
    }
}

fn mean_std(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    (
        mean,
        (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n).sqrt(),
    )
}

fn jittered_grid(n: usize, jitter: f64, seed: u64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let h = 1.0 / n as f64;
    let mut rng = SplitMix64::new(seed);
    let id = |i: usize, j: usize| (i * (n + 1) + j) as u32;
    let mut positions = Vec::new();
    for i in 0..=n {
        for j in 0..=n {
            let mut p = [i as f64 * h, j as f64 * h, 0.0];
            let inside = |x: f64| x > 0.5 * h && x < 1.0 - 0.5 * h;
            if inside(p[0]) && inside(p[1]) {
                p[0] += jitter * h * rng.uniform(-1.0, 1.0);
                p[1] += jitter * h * rng.uniform(-1.0, 1.0);
            }
            positions.push(p);
        }
    }
    let mut faces = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    (positions, faces)
}

fn latitude_exact(psi0: f64, alpha: f64, t: f64) -> f64 {
    let r = 1.0 / (4.0 * PI).sqrt();
    let rhs = |psi: f64| psi.cos() / (r * r * psi.sin()) * (4.0 / (r * psi.sin()) - alpha);
    let n = 20000;
    let h = t / n as f64;
    let mut psi = psi0;
    for _ in 0..n {
        let k1 = rhs(psi);
        let k2 = rhs(psi + 0.5 * h * k1);
        let k3 = rhs(psi + 0.5 * h * k2);
        let k4 = rhs(psi + h * k3);
        psi += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    psi
}

/// The treatments of the laws: all but the adaptive refinement (the flat meshes have no
/// curvature, the sphere a uniform one), run directly (flips and the optimal Delaunay
/// triangulation too, without the Session).
fn law_treats() -> Vec<Treat> {
    let mut out: Vec<Treat> = all_treats()
        .into_iter()
        .filter(|t| !matches!(t.refine, Refine::Adaptive(_) | Refine::Flips))
        .map(|t| Treat {
            session: false,
            ..t
        })
        .collect();
    out.push(Treat {
        tag: "flipS",
        refine: Refine::Flips,
        spread: true,
        session: false,
    });
    out.insert(
        1,
        Treat {
            tag: "flip",
            refine: Refine::Flips,
            spread: false,
            session: false,
        },
    );
    out
}

fn laws() {
    let threads = env_usize("STUDY_THREADS", 32);
    // Flat circle law: only the length term acts; the curvature is zero, so S changes nothing.
    let (positions, faces) = jittered_grid(32, 0.3, 7);
    let grid = Input::new("grid32", positions, faces).expect("grid");
    let center = [0.511, 0.493];
    let (rho0, alpha, t_end, steps) = (0.3f64, 1.0f64, 0.03f64, 96usize);
    let exact = (rho0 * rho0 - 2.0 * alpha * t_end).sqrt();
    let flat: Vec<Treat> = law_treats().into_iter().filter(|t| !t.spread).collect();
    let builds: Vec<Result<Built, String>> =
        par_map(&flat, threads, |t| build(&grid, t.refine, None))
            .into_iter()
            .map(|r| r.and_then(|x| x))
            .collect();
    println!(
        "CIRCLE LAW (jittered grid 32, rho0 = 0.3, alpha_L = 1, t = 0.03, {steps} fixed steps, Dirichlet): exact radius {exact:.5}"
    );
    let ks: Vec<usize> = (0..flat.len()).collect();
    let res = par_map(&ks, threads, |&k| -> Result<String, String> {
        let b = builds[k].as_ref().map_err(|e| e.clone())?;
        let radius =
            |p: [f64; 3]| ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2)).sqrt();
        let phi0: Vec<f64> = b.positions.iter().map(|&p| radius(p) - rho0).collect();
        let mut flow = Flow::new(&b.mesh, &b.geometry, &phi0).map_err(|e| e.to_string())?;
        let problem = CutProblem::new(Energy::Dirichlet, alpha);
        let options = FlowOptions::new(StepSize::Fixed(t_end / steps as f64));
        for _ in 0..steps {
            flow.step(&problem, &options).map_err(|e| e.to_string())?;
        }
        let rs: Vec<f64> = cut_positions(b, &flow.level_set())
            .into_iter()
            .map(radius)
            .collect();
        let (mean, std) = mean_std(&rs);
        Ok(format!(
            "LAW circle treat={} V={} radius_error={:+.3e} roundness={:.2e}",
            flat[k].tag,
            b.mesh.n_vertices(),
            mean / exact - 1.0,
            std / mean
        ))
    });
    for r in res {
        println!(
            "{}",
            r.and_then(|x| x)
                .unwrap_or_else(|e| format!("LAW circle ERROR {e}"))
        );
    }
    // Sphere: the latitude law (ODE) and the balanced latitude, on the polyhedral icosphere.
    for level in [3u32, 4] {
        let ico = icosphere(level).expect("icosphere");
        let sphere = Input::new(&format!("ico{level}"), ico.positions, ico.faces).expect("sphere");
        let treats = law_treats();
        let builds: Vec<Result<Built, String>> =
            par_map(&treats, threads, |t| build(&sphere, t.refine, None))
                .into_iter()
                .map(|r| r.and_then(|x| x))
                .collect();
        let r = 1.0 / (4.0 * PI).sqrt();
        let star = (4.0 / (20.0 * r)).asin();
        println!(
            "\nSPHERE icosphere({level}): ODE (alpha_L = 8, psi0 = 0.8, t = 4e-3, 50 fixed steps; exact {:.5}) and balance (alpha_L = 20, psi* = {star:.5}, adaptive steps to rest from psi0 = 0.5 and 1.1), Dirichlet",
            latitude_exact(0.8, 8.0, 4e-3)
        );
        let jobs: Vec<(usize, u8)> = (0..treats.len())
            .flat_map(|k| [(k, 0u8), (k, 1), (k, 2)])
            .collect();
        let res = par_map(&jobs, threads, |&(k, which)| -> Result<String, String> {
            let b = builds[k].as_ref().map_err(|e| e.clone())?;
            let treat = &treats[k];
            let psi0 = [0.8, 0.5, 1.1][which as usize];
            let alpha = if which == 0 { 8.0 } else { 20.0 };
            let phi0: Vec<f64> = b.positions.iter().map(|&p| polar(p) - psi0).collect();
            let mut flow = Flow::new(&b.mesh, &b.geometry, &phi0).map_err(|e| e.to_string())?;
            let mut problem = CutProblem::new(Energy::Dirichlet, alpha);
            if treat.spread {
                problem.curvature = CurvatureSource::PerFace(b.spread.clone());
            }
            let (steps, stop) = if which == 0 {
                let options = FlowOptions::new(StepSize::Fixed(4e-3 / 50.0));
                for _ in 0..50 {
                    flow.step(&problem, &options).map_err(|e| e.to_string())?;
                }
                (50, "fixed".to_string())
            } else {
                let stop = StopRule {
                    max_steps: 3000,
                    motion_tolerance: 1e-7,
                    energy_tolerance: 0.0,
                    energy_window: 1,
                };
                let rep = flow
                    .run(&problem, &flow_options(), &stop)
                    .map_err(|e| e.to_string())?;
                (rep.steps.len(), format!("{:?}", rep.stop))
            };
            let ps: Vec<f64> = cut_positions(b, &flow.level_set())
                .into_iter()
                .map(polar)
                .collect();
            if ps.is_empty() {
                return Err(format!("{}: the cut vanished", treat.tag));
            }
            let (mean, std) = mean_std(&ps);
            let target = if which == 0 {
                latitude_exact(0.8, 8.0, 4e-3)
            } else {
                star
            };
            Ok(format!(
                "LAW {} level={level} treat={} V={} psi0={psi0} error={:+.3e} spread={:.2e} steps={steps} stop={stop}",
                if which == 0 { "ode" } else { "balance" },
                treat.tag,
                b.mesh.n_vertices(),
                mean - target,
                std
            ))
        });
        for r in res {
            println!(
                "{}",
                r.and_then(|x| x)
                    .unwrap_or_else(|e| format!("LAW sphere ERROR {e}"))
            );
        }
    }
}

// ---------------------------------------------------------------------------------------------
// a reproducer of a failed way back

/// `repro MESH TREAT ENERGY ALPHA_FACTOR START`: the run of one case of `real`, and if the way
/// back fails, every segment that fails alone, with its nodes and their places on the input.
fn repro(args: &[String]) {
    let [name, tag, energy_name, factor, start] = args else {
        eprintln!("usage: remesh_study repro MESH TREAT ENERGY ALPHA_FACTOR START");
        return;
    };
    let input = Input::load(name).expect("load");
    let treat = *all_treats()
        .iter()
        .find(|t| t.tag == tag)
        .expect("unknown treatment");
    let budget = match treat.refine {
        Refine::Adaptive(f) => Some(
            build(&input, Refine::Uniform(f), None)
                .expect("uniform build")
                .info
                .added,
        ),
        _ => None,
    };
    let flips = build(&input, Refine::Flips, None).expect("flips");
    let built = build(&input, treat.refine, budget).expect("build");
    println!(
        "REPRO mesh={name} treat={tag} V={} added={} bound={:.6e}",
        built.mesh.n_vertices(),
        built.info.added,
        built.info.bound
    );
    let (energy, energy_name) = energy_of(energy_name);
    let start: usize = start.parse().expect("start");
    let factor: f64 = factor.parse().expect("alpha factor");
    let auto = alpha_auto(&input, &flips, &energy, start).expect("alpha");
    let case = Case {
        energy,
        energy_name,
        alpha_factor: factor,
        alpha: factor * auto,
        start,
    };
    println!("REPRO case={} alpha={:.6e}", case.tag(), case.alpha);
    let phi0: Vec<f64> = built
        .positions
        .iter()
        .map(|&p| input.start(case.start, p))
        .collect();
    let mut flow = Flow::new(&built.mesh, &built.geometry, &phi0).expect("flow");
    flow.redistance().expect("redistance");
    let mut problem = CutProblem::new(case.energy.clone(), case.alpha);
    if treat.spread {
        problem.curvature = CurvatureSource::PerFace(built.spread.clone());
    }
    let report = flow
        .run(&problem, &flow_options(), &stop_rule())
        .expect("run");
    println!("REPRO {}", summarize_run(&report, &input));
    let phi = flow.level_set();
    let Some(tri) = built.tri.as_ref() else {
        println!("REPRO no way back for the caller's triangles");
        return;
    };
    let export = built.export.as_ref().unwrap();
    let cut = Cut::from_level_set(&export.mesh, &phi, 1e-3).expect("cut");
    let nodes: Vec<SurfacePoint> = cut
        .points
        .iter()
        .map(|p| export.edge_point(p.edge as usize, p.t))
        .collect();
    let segments: Vec<Segment> = cut
        .segments
        .iter()
        .map(|s| Segment {
            face: s.face,
            nodes: s.points,
        })
        .collect();
    match tri.segments_on_input(&nodes, &segments, &EdgeTraceOptions::default()) {
        Ok(b) => {
            println!(
                "REPRO way back OK: {} segments, max relative mismatch {:.2e}",
                segments.len(),
                b.max_relative_mismatch
            );
            return;
        }
        Err(e) => println!("REPRO way back FAILED: {e}"),
    }
    let mut n_failed = 0;
    for (k, seg) in segments.iter().enumerate() {
        let one = [*seg];
        let Err(e) = tri.segments_on_input(&nodes, &one, &EdgeTraceOptions::default()) else {
            continue;
        };
        n_failed += 1;
        println!("REPRO segment {k} in face {} of the triangulation: {e}", seg.face);
        for &n in &seg.nodes {
            let p = cut.points[n as usize];
            let [a, b] = export.mesh.edge_vertices(p.edge as usize);
            println!(
                "  node {n}: export edge {} ({a}, {b}) t={:.17e} -> {:?} on the triangulation -> {:?} on the input",
                p.edge,
                p.t,
                nodes[n as usize],
                tri.point_on_input(nodes[n as usize])
            );
        }
        let corners: Vec<String> = (0..3)
            .map(|c| {
                let e = tri.edge(3 * seg.face as usize + c);
                format!("edge {e} length {:.6e}", tri.edge_length(e))
            })
            .collect();
        println!("  face {} sides: {}", seg.face, corners.join(", "));
        for tol in [1e-4, 1e-2, 1e-1] {
            let r = tri.segments_on_input(
                &nodes,
                &one,
                &EdgeTraceOptions {
                    tolerance: tol,
                    ..EdgeTraceOptions::default()
                },
            );
            println!(
                "  with tolerance {tol:.0e}: {}",
                match r {
                    Ok(b) => format!("OK, relative mismatch {:.2e}", b.max_relative_mismatch),
                    Err(e) => format!("{e}"),
                }
            );
        }
    }
    println!("REPRO {n_failed} of {} segments fail alone", segments.len());
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("real") => {
            for name in &args[1..] {
                real(name);
            }
        }
        Some("speed") => speed(&args[1..]),
        Some("laws") => laws(),
        Some("repro") => repro(&args[1..]),
        _ => eprintln!("usage: remesh_study real MESH... | speed MESH... | laws"),
    }
}
