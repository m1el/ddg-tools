//! The published UV maps of the optimal-cones dataset against our flattening and our cones.
//!
//! ```text
//! cargo run --release -p ddg-uv --example published_uvs [name ...]
//! ```
//!
//! Reads the 23 meshes of `data/optimal-cones/optimal-cones-data/data` (or `DDG_DATASET`); the
//! arguments restrict the meshes to those whose name equals one of them. Prints three tables.
//!
//! * A: the published data as it is: cones, total cone angle, Gauss-Bonnet of the cone file,
//!   the dataset's area distortion of the `vt` coordinates, flipped triangles, charts.
//! * B: our Boundary First Flattening with exactly the file's cones (positions and angles), for
//!   `u = 0` on the boundary and for optimized boundary values, compared triangle by triangle
//!   with the published UVs.
//! * C: our algorithm against theirs at equal or smaller total cone angle: C1 the penalized
//!   optimum, C2 the recipe the authors evidently used (solve, merge clusters, re-optimize the
//!   angles at the merged vertices), each with `lambda` found by a coarse scan plus bisection
//!   (cap: the published total cone angle `sum |angle|` times 1.001; the candidate with the
//!   smallest solver distortion among those within the cap is flattened), C3 the same recipes
//!   at the number of cones of the file (no cap on the angle there: read the angle column).
//!   For meshes with a boundary Table C runs for both boundary treatments (`u = 0` and
//!   optimized boundary values); the first block is the treatment under which the file's
//!   angles are best reproduced.
//!
//! The measure, used everywhere: `D = sum_t A_t u_t^2`, `u_t = 1/2 log(B_t/A_t)`, `A_t` the 3D
//! area (the weights `A_t / sum A` sum to 1; the ratio `B_t/A_t` uses the areas as they are),
//! `B_t` the area in the UV plane. `D*` is `D` after the best global scale of the UV plane
//! (the area-weighted mean of `u` removed). Flipped (orientation opposite to the majority) and
//! degenerate UV triangles are left out of the sums and counted.
//!
//! Environment: `DDG_DATASET` (directory), `DDG_BUDGET_S` (seconds of searching per mesh and
//! variant, default 1500), `DDG_ITERS` (bisection steps, default 20), `DDG_SLACK` (relative slack on the cap, default 1e-3), `DDG_GRID` (values of
//! `lambda` in the coarse scan, default 12; no scan on meshes of more than 100 000 faces),
//! `DDG_RECIPES` (comma-separated recipes of Table C, default `c1,c2`; `c1fb` is C1 with the
//! feedback loop of the cone solver at every probe: solve, flatten, lower the distortion weight
//! where the solver's `u` is far from the flattening's, repeat; its note lists the rounds and
//! the reweighted vertices), `DDG_FB_THRESHOLD`, `DDG_FB_ROUNDS` (the loop's threshold and
//! round limit, defaults those of `FeedbackOptions`).

use std::f64::consts::PI;
use std::path::PathBuf;
use std::time::Instant;

use ddg_bff::{
    ConeFlattenInput, ConeFlattenOptions, CutSource, GaussBonnetPolicy, ScaleFactorOptions,
    flatten_with_cones,
};
use ddg_core::io::{Cone, read_cones, read_obj};
use ddg_core::operators::angle_defects;
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    BoundaryTreatment, ConeProblem, FeedbackOptions, Parallelism, Solution, SolveOptions, Status,
    evaluate_cones, solve_with_boundary, solve_with_feedback,
};
use ddg_uv::realized_face_log_scale;

/// Relative slack on the published total cone angle (`DDG_SLACK`, default 1e-3) (the file's angles have six digits and the
/// same positions re-optimized sum to it only to about 4e-5, `bump`).
fn total_slack() -> f64 {
    env_f64("DDG_SLACK", 1e-3)
}
/// Two `D*` are called equal if they differ by less than this relative amount.
const EQUAL: f64 = 1e-3;

fn dataset_dir() -> PathBuf {
    std::env::var_os("DDG_DATASET").map_or_else(
        || {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../data/optimal-cones/optimal-cones-data/data")
        },
        PathBuf::from,
    )
}

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

struct Data {
    name: String,
    mesh: TriMesh,
    geometry: Geometry,
    cones: Vec<Cone>,
    /// Published UV per corner (`3 f + k`).
    uvs: Vec<[f64; 2]>,
    boundary: bool,
    /// Decimals of the published `vt` coordinates.
    decimals: usize,
}

fn load(name: &str) -> Result<Data, String> {
    let dir = dataset_dir();
    let obj = read_obj(dir.join(format!("{name}.obj"))).map_err(|e| e.to_string())?;
    let faces = obj.triangles().ok_or("not a triangle mesh")?;
    let mesh = TriMesh::new(obj.positions.len(), faces).map_err(|e| e.to_string())?;
    let geometry = Geometry::from_positions(&mesh, &obj.positions).map_err(|e| e.to_string())?;
    let cones = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices())
        .map_err(|e| e.to_string())?;
    let uvs = obj.corner_uvs().ok_or("a corner has no vt")?;
    let boundary = mesh.n_boundary_edges() > 0;
    let text =
        std::fs::read_to_string(dir.join(format!("{name}.obj"))).map_err(|e| e.to_string())?;
    let decimals = text
        .lines()
        .filter(|l| l.starts_with("vt "))
        .flat_map(|l| l.split_whitespace().skip(1))
        .map(|w| {
            w.split('.')
                .nth(1)
                .map_or(0, |f| f.trim_end_matches(char::is_alphabetic).len())
        })
        .max()
        .unwrap_or(0);
    Ok(Data {
        decimals,
        name: name.to_string(),
        mesh,
        geometry,
        cones,
        uvs,
        boundary,
    })
}

// ---------------------------------------------------------------- the measure

struct Meas {
    /// `u_t`, NaN where flipped or degenerate.
    u: Vec<f64>,
    d: f64,
    dstar: f64,
    flipped: usize,
    degenerate: usize,
    /// The majority orientation of the UV triangles is opposite to the mesh's.
    mirrored: bool,
    /// The area-weighted mean of `u` over the triangles that count.
    mean: f64,
}

fn measure(d: &Data, uvs: &[[f64; 2]]) -> Meas {
    let areas = d.geometry.face_areas();
    let total: f64 = areas.iter().sum();
    let nf = areas.len();
    let signed: Vec<f64> = (0..nf)
        .map(|f| {
            let (p, q, r) = (uvs[3 * f], uvs[3 * f + 1], uvs[3 * f + 2]);
            0.5 * ((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
        })
        .collect();
    let positive = signed.iter().filter(|b| **b > 0.0).count();
    let negative = signed.iter().filter(|b| **b < 0.0).count();
    let s = if negative > positive { -1.0 } else { 1.0 };
    let (mut u, mut flipped, mut degenerate) = (vec![f64::NAN; nf], 0, 0);
    let (mut sw, mut swu, mut swuu) = (0.0, 0.0, 0.0);
    for f in 0..nf {
        let b = s * signed[f];
        if !b.is_finite()
            || b == 0.0
            || areas[f].partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
        {
            degenerate += 1;
        } else if b < 0.0 {
            flipped += 1;
        } else {
            let v = 0.5 * (b / areas[f]).ln();
            if !v.is_finite() {
                degenerate += 1;
                continue;
            }
            u[f] = v;
            let w = areas[f] / total;
            sw += w;
            swu += w * v;
            swuu += w * v * v;
        }
    }
    let dstar = if sw > 0.0 {
        swuu - swu * swu / sw
    } else {
        f64::NAN
    };
    Meas {
        u,
        d: swuu,
        dstar: dstar.max(0.0),
        flipped,
        degenerate,
        mirrored: s < 0.0,
        mean: if sw > 0.0 { swu / sw } else { f64::NAN },
    }
}

/// Triangle by triangle: `u_a - u_b` after removing each one's area-weighted mean over the
/// triangles that count in both. Returns (largest |diff|, weighted rms, weighted correlation,
/// number of triangles compared).
fn compare(d: &Data, a: &[f64], b: &[f64]) -> (f64, f64, f64, usize) {
    let areas = d.geometry.face_areas();
    let total: f64 = areas.iter().sum();
    let idx: Vec<usize> = (0..a.len())
        .filter(|&f| a[f].is_finite() && b[f].is_finite())
        .collect();
    let sw: f64 = idx.iter().map(|&f| areas[f] / total).sum();
    let ma: f64 = idx.iter().map(|&f| areas[f] / total * a[f]).sum::<f64>() / sw;
    let mb: f64 = idx.iter().map(|&f| areas[f] / total * b[f]).sum::<f64>() / sw;
    let (mut mx, mut ss, mut saa, mut sbb, mut sab) = (0.0f64, 0.0, 0.0, 0.0, 0.0);
    for &f in &idx {
        let w = areas[f] / total;
        let (x, y) = (a[f] - ma, b[f] - mb);
        mx = mx.max((x - y).abs());
        ss += w * (x - y) * (x - y);
        saa += w * x * x;
        sbb += w * y * y;
        sab += w * x * y;
    }
    (mx, (ss / sw).sqrt(), sab / (saa * sbb).sqrt(), idx.len())
}

/// Where the squared difference of two `u` (means removed) sits: its share on triangles with a
/// boundary vertex and on triangles with a cone vertex, and the area share of those triangles.
/// Returns `[(error share, area share); 2]`.
fn error_shares(d: &Data, a: &[f64], b: &[f64]) -> [(f64, f64); 2] {
    let areas = d.geometry.face_areas();
    let total: f64 = areas.iter().sum();
    let idx: Vec<usize> = (0..a.len())
        .filter(|&f| a[f].is_finite() && b[f].is_finite())
        .collect();
    let sw: f64 = idx.iter().map(|&f| areas[f] / total).sum();
    let ma: f64 = idx.iter().map(|&f| areas[f] / total * a[f]).sum::<f64>() / sw;
    let mb: f64 = idx.iter().map(|&f| areas[f] / total * b[f]).sum::<f64>() / sw;
    let mut is_cone = vec![false; d.mesh.n_vertices()];
    for c in &d.cones {
        is_cone[c.vertex as usize] = true;
    }
    let mut out = [(0.0, 0.0); 2];
    let mut all = 0.0;
    for &f in &idx {
        let w = areas[f] / total;
        let e = w * ((a[f] - ma) - (b[f] - mb)).powi(2);
        all += e;
        let vs = d.mesh.faces()[f];
        let sets = [
            vs.iter().any(|&v| d.mesh.is_boundary_vertex(v as usize)),
            vs.iter().any(|&v| is_cone[v as usize]),
        ];
        for (o, hit) in out.iter_mut().zip(sets) {
            if hit {
                o.0 += e;
                o.1 += w;
            }
        }
    }
    for o in &mut out {
        o.0 /= all;
        o.1 /= sw;
    }
    out
}

/// Connected pieces of the surface in the UV plane: faces are joined across an edge if both
/// endpoints have the same UV in both faces.
fn count_charts(d: &Data, uvs: &[[f64; 2]]) -> usize {
    let m = &d.mesh;
    let mut parent: Vec<usize> = (0..m.n_faces()).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for h in 0..m.n_halfedges() {
        if let Some(t) = m.twin(h)
            && h < t
            && uvs[h] == uvs[m.next(t)]
            && uvs[m.next(h)] == uvs[t]
        {
            let (a, b) = (find(&mut parent, m.face(h)), find(&mut parent, m.face(t)));
            parent[a] = b;
        }
    }
    (0..m.n_faces())
        .filter(|&f| find(&mut parent, f) == f)
        .count()
}

/// The cone angle `2 pi - (sum of the UV corner angles)` at every cone of the file, and its
/// largest difference from the file's angle.
fn implied_angle_difference(d: &Data, uvs: &[[f64; 2]]) -> f64 {
    let mut sum = vec![0.0f64; d.mesh.n_vertices()];
    for (f, face) in d.mesh.faces().iter().enumerate() {
        for k in 0..3 {
            let (p, q, r) = (
                uvs[3 * f + k],
                uvs[3 * f + (k + 1) % 3],
                uvs[3 * f + (k + 2) % 3],
            );
            let (a, b) = ([q[0] - p[0], q[1] - p[1]], [r[0] - p[0], r[1] - p[1]]);
            let cross = a[0] * b[1] - a[1] * b[0];
            let dot = a[0] * b[0] + a[1] * b[1];
            sum[face[k] as usize] += cross.abs().atan2(dot);
        }
    }
    d.cones
        .iter()
        .map(|c| ((2.0 * PI - sum[c.vertex as usize]) - c.angle).abs())
        .fold(0.0, f64::max)
}

// ---------------------------------------------------------------- solving and flattening

fn options(merge: bool) -> SolveOptions<'static> {
    SolveOptions {
        merge_cone_clusters: merge,
        max_merge_cluster_size: usize::MAX,
        ..SolveOptions::default()
    }
}

type Weights = Option<Vec<f64>>;

fn penalized(
    d: &Data,
    lambda: f64,
    bt: BoundaryTreatment,
    merge: bool,
    dw: &Weights,
) -> Result<Solution, String> {
    let problem = ConeProblem {
        distortion_weights: dw.clone(),
        ..ConeProblem::new(&d.mesh, &d.geometry, lambda)
    };
    solve_with_boundary(&problem, bt, &options(merge)).map_err(|e| e.to_string())
}

/// Angles of least distortion at the given vertices (free cone sites there, none elsewhere).
fn refit(
    d: &Data,
    vertices: &[u32],
    bt: BoundaryTreatment,
    dw: &Weights,
) -> Result<Solution, String> {
    let mut costs = vec![f64::INFINITY; d.mesh.n_vertices()];
    for &v in vertices {
        costs[v as usize] = 0.0;
    }
    let problem = ConeProblem {
        cone_costs: Some(costs),
        distortion_weights: dw.clone(),
        ..ConeProblem::new(&d.mesh, &d.geometry, 1.0)
    };
    solve_with_boundary(&problem, bt, &options(false)).map_err(|e| e.to_string())
}

/// The solution with the given cones fixed (angles from the file), nothing else free except the
/// boundary values for `Free`.
fn prescribed(
    d: &Data,
    cones: &[Cone],
    bt: BoundaryTreatment,
    dw: &Weights,
) -> Result<Solution, String> {
    let problem = ConeProblem {
        cone_costs: Some(vec![f64::INFINITY; d.mesh.n_vertices()]),
        prescribed_cones: cones.to_vec(),
        distortion_weights: dw.clone(),
        ..ConeProblem::new(&d.mesh, &d.geometry, 0.0)
    };
    solve_with_boundary(&problem, bt, &options(false)).map_err(|e| e.to_string())
}

struct Flat {
    uvs: Vec<[f64; 2]>,
    charts: usize,
    secs: f64,
    note: String,
}

fn flatten(d: &Data, cones: &[Cone], u: Option<&[f64]>) -> Result<Flat, String> {
    let start = Instant::now();
    let run = |policy| {
        flatten_with_cones(
            &ConeFlattenInput {
                mesh: &d.mesh,
                geometry: &d.geometry,
                cones,
                scale_factors: u,
                cut: CutSource::Plan {
                    may_cut: None,
                    must_cut: &[],
                },
                anchors: &[],
            },
            &ConeFlattenOptions {
                scale_factors: ScaleFactorOptions {
                    gauss_bonnet: policy,
                    ..ScaleFactorOptions::default()
                },
                overlap_samples_per_triangle: 0.0,
                ..ConeFlattenOptions::default()
            },
        )
    };
    let mut note = String::new();
    let flat = match run(GaussBonnetPolicy::Reject) {
        Err(ddg_bff::BffError::GaussBonnetViolated { .. }) if u.is_none() => {
            note = "GB spread over vertices".into();
            run(GaussBonnetPolicy::SpreadOverVertices)
        }
        other => other,
    }
    .map_err(|e| e.to_string())?;
    Ok(Flat {
        charts: flat.n_charts(),
        uvs: flat.corner_uvs,
        secs: start.elapsed().as_secs_f64(),
        note,
    })
}

/// The options of the feedback loop (`DDG_FB_THRESHOLD`, `DDG_FB_ROUNDS`).
fn feedback_options() -> FeedbackOptions {
    let d = FeedbackOptions::default();
    FeedbackOptions {
        threshold: env_f64("DDG_FB_THRESHOLD", d.threshold),
        max_rounds: env_f64("DDG_FB_ROUNDS", d.max_rounds as f64) as usize,
        max_reweighted_fraction: env_f64("DDG_FB_FRACTION", d.max_reweighted_fraction),
        reduced_weight: env_f64("DDG_FB_REDUCED", d.reduced_weight),
        cone_free_probe: env_f64("DDG_FB_PROBE", 1.0) != 0.0,
    }
}

/// The recipes of Table C to run (`DDG_RECIPES`).
fn recipe(name: &str) -> bool {
    std::env::var("DDG_RECIPES")
        .unwrap_or_else(|_| "c1,c2".to_owned())
        .split(',')
        .any(|r| r.trim() == name)
}

fn total_angle(cones: &[Cone]) -> f64 {
    cones.iter().map(|c| c.angle.abs()).sum()
}

fn short(s: String) -> String {
    s.chars().take(90).collect()
}

// ---------------------------------------------------------------- Table C: searches

struct Probe {
    lambda: f64,
    cones: Vec<Cone>,
    u: Vec<f64>,
    total: f64,
    /// The solver's `1/2 sum m u^2` of these cones.
    dist: f64,
    converged: bool,
    /// What else to print with the probe (the feedback loop's rounds).
    note: String,
}

struct Search {
    best: Option<Probe>,
    probes: Vec<Probe>,
    notes: Vec<String>,
}

/// `eval(lambda)` gives the cones for that `lambda`. The total cone angle is assumed to fall
/// with `lambda`, which holds for the penalized optimum but not for merging and re-optimizing
/// the angles; so a coarse scan of `grid` values of `lambda` precedes the bisection, and the
/// feasible probe (total at most `cap`) with the smallest solver distortion is kept.
fn search(
    cap: f64,
    iters: usize,
    grid: usize,
    budget: f64,
    start: Instant,
    mut eval: impl FnMut(f64) -> Result<Probe, String>,
) -> Search {
    let mut out = Search {
        best: None,
        probes: Vec::new(),
        notes: Vec::new(),
    };
    let (lo, hi) = (1e-9f64, 1e6f64);
    let feasible = |p: Probe, out: &mut Search| -> bool {
        let ok = p.total <= cap && !p.cones.is_empty();
        let better = out.best.as_ref().is_none_or(|b| p.dist < b.dist);
        let free = p.total <= cap;
        if !p.converged {
            out.notes
                .push(format!("lambda {:.3e}: not converged", p.lambda));
        }
        if ok && better {
            out.best = Some(Probe {
                lambda: p.lambda,
                cones: p.cones.clone(),
                u: p.u.clone(),
                total: p.total,
                dist: p.dist,
                converged: p.converged,
                note: p.note.clone(),
            });
        }
        out.probes.push(p);
        free
    };
    let mut run = |l: f64, out: &mut Search| -> Option<bool> {
        let t0 = Instant::now();
        let r = eval(l);
        if let Ok(p) = &r {
            eprintln!(
                "  probe lambda {l:.3e}: {} cones, total {:.3}, {:.1}s",
                p.cones.len(),
                p.total,
                t0.elapsed().as_secs_f64()
            );
        }
        match r {
            Ok(p) => Some(feasible(p, out)),
            Err(e) => {
                out.notes.push(format!("lambda {l:.3e}: {}", short(e)));
                None
            }
        }
    };
    let Some(hi_ok) = run(hi, &mut out) else {
        return out;
    };
    if !hi_ok {
        out.notes.push(format!(
            "even lambda = {hi:e} has a total angle above the cap"
        ));
        return out;
    }
    let lo_ok = run(lo, &mut out);
    if lo_ok == Some(true) {
        out.notes.push(format!(
            "lambda = {lo:e} already within the cap: fewest restrictions"
        ));
        return out;
    }
    for g in 0..grid {
        if start.elapsed().as_secs_f64() > budget {
            out.notes.push("time budget used up".into());
            break;
        }
        let l = 10f64.powf(-5.0 + 5.5 * g as f64 / (grid.max(2) - 1) as f64);
        run(l, &mut out);
    }
    let (mut a, mut b) = (lo, hi);
    for _ in 0..iters {
        if start.elapsed().as_secs_f64() > budget {
            out.notes.push("time budget used up".into());
            break;
        }
        let mid = (a * b).sqrt();
        match run(mid, &mut out) {
            Some(true) => b = mid,
            Some(false) => a = mid,
            None => break,
        }
    }
    out
}

/// A probe whose number of cones is `target`: bisect `lambda` between probes that bracket it.
fn equal_count(
    search: &mut Search,
    target: usize,
    iters: usize,
    budget: f64,
    start: Instant,
    mut eval: impl FnMut(f64) -> Result<Probe, String>,
) -> Option<usize> {
    let exact = |s: &Search| s.probes.iter().position(|p| p.cones.len() == target);
    if let Some(i) = exact(search) {
        return Some(i);
    }
    let mut ps: Vec<(f64, usize)> = search
        .probes
        .iter()
        .map(|p| (p.lambda, p.cones.len()))
        .collect();
    ps.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    // Adjacent pair (low lambda has more cones than target, higher lambda fewer).
    let pair = ps
        .windows(2)
        .find(|w| w[0].1 > target && w[1].1 < target)
        .map(|w| (w[0].0, w[1].0));
    if let Some((mut a, mut b)) = pair {
        for _ in 0..iters {
            if start.elapsed().as_secs_f64() > budget {
                break;
            }
            let mid = (a * b).sqrt();
            let Ok(p) = eval(mid) else { break };
            let n = p.cones.len();
            search.probes.push(p);
            if n == target {
                return Some(search.probes.len() - 1);
            }
            if n > target {
                a = mid;
            } else {
                b = mid;
            }
        }
    }
    // The closest count attained, fewer preferred on ties.
    (0..search.probes.len()).min_by_key(|&i| {
        let n = search.probes[i].cones.len();
        (n.abs_diff(target), n > target)
    })
}

// ---------------------------------------------------------------- one mesh

fn verdict(ours: f64, theirs: f64) -> String {
    let r = ours / theirs;
    if !r.is_finite() {
        "?".into()
    } else if (r - 1.0).abs() <= EQUAL {
        format!("equal (within {:.1}%; ratio {r:.4})", 100.0 * EQUAL)
    } else if r < 1.0 {
        format!("LOWER by {:.2}% (ratio {r:.4})", 100.0 * (1.0 - r))
    } else {
        format!("HIGHER by {:.2}% (ratio {r:.4})", 100.0 * (r - 1.0))
    }
}

struct Tables {
    a: Vec<String>,
    b: Vec<String>,
    c: Vec<String>,
    /// The boundary treatment whose re-optimized angles are closest to the file's.
    family: BoundaryTreatment,
}

const A_HEAD: &str = "mesh                  verts  faces | cones  sum|angle| | Gauss-Bonnet of file |           D          D* | flipped degen | charts";
const B_HEAD: &str = "mesh                  boundary |   D*(ours)  D*(pub)  ratio | flip degen charts | max|du|  rms du     corr | D* rounded, rms floor of the rounding | where the squared difference sits | status, time";
const C_HEAD: &str = "mesh / recipe                  | cones  sum|angle|       lambda | solver D |          D          D* | flip degen | D* vs published | time";

fn row_c(d: &Data, label: &str, p: &Probe, pub_dstar: f64, note: &str) -> String {
    let (cones, solver) = (&p.cones[..], p.dist);
    let lam = format!("{:.4e}", p.lambda);
    let head = format!(
        "{:<30} | {:>5} {:>9.4} {:>12} | {:>8.5} |",
        label,
        cones.len(),
        total_angle(cones),
        lam,
        solver
    );
    let t = Instant::now();
    match flatten(d, cones, Some(&p.u)) {
        Err(e) => format!("{head} flattening FAILED: {} {note}", short(e)),
        Ok(f) => {
            let m = measure(d, &f.uvs);
            format!(
                "{head} {:>11.4e} {:>11.4e} | {:>4} {:>5} | {:<15} | {:.1}s {note}{}",
                m.d,
                m.dstar,
                m.flipped,
                m.degenerate,
                verdict(m.dstar, pub_dstar),
                f.secs.max(t.elapsed().as_secs_f64()),
                if f.charts > 1 {
                    " [several charts]"
                } else {
                    ""
                }
            ) + &p.note
        }
    }
}

/// With `treatment` the recipes of Table C run for that boundary treatment only (no Tables A, B);
/// without it, A and B are produced (unless `tips_variant`) and Table C runs for the boundary
/// family of the file.
fn run_mesh(d: &Data, tips_variant: bool, treatment: Option<BoundaryTreatment>) -> Tables {
    let mut t = Tables {
        a: vec![],
        b: vec![],
        c: vec![],
        family: BoundaryTreatment::Fixed,
    };
    let start = Instant::now();
    let budget = env_f64("DDG_BUDGET_S", 1500.0);
    let iters = env_f64("DDG_ITERS", 20.0) as usize;
    // The scan of lambda is skipped on the very large meshes (it costs a solve per value).
    let grid = if d.mesh.n_faces() > 100_000 {
        0
    } else {
        env_f64("DDG_GRID", 12.0) as usize
    };
    let nv = d.mesh.n_vertices();
    let pub_total = total_angle(&d.cones);
    let tag = &d.name;

    // ---- Table A
    let pm = measure(d, &d.uvs);
    let chi = d.mesh.euler_characteristic();
    let gb = if d.boundary {
        "n/a (boundary)".to_string()
    } else {
        let sum: f64 = d.cones.iter().map(|c| c.angle).sum();
        format!("sum - 2 pi chi = {:+.2e}", sum - 2.0 * PI * chi as f64)
    };
    let charts = count_charts(d, &d.uvs);
    let tables_ab = !tips_variant && treatment.is_none();
    if tables_ab {
        t.a.push(format!(
        "{:<20} {:>7} {:>6} | {:>5} {:>10.4} | {:<20} | {:>11.4e} {:>11.4e} | {:>7} {:>5} | {}{}",
        tag,
        nv,
        d.mesh.n_faces(),
        d.cones.len(),
        pub_total,
        gb,
        pm.d,
        pm.dstar,
        pm.flipped,
        pm.degenerate,
        if charts == 1 { "one chart".to_string() } else { format!("{charts} charts") },
        if pm.mirrored { " (UV orientation mirrored)" } else { "" }
    ));
    }

    // ---- Table B (not repeated for the variant)
    let treatments: Vec<(&str, BoundaryTreatment)> = if d.boundary {
        vec![
            ("u=0", BoundaryTreatment::Fixed),
            ("optimized", BoundaryTreatment::Free),
        ]
    } else {
        vec![("closed", BoundaryTreatment::Fixed)]
    };
    let mut best_rms = (f64::INFINITY, "");
    for (name, bt) in treatments.iter().filter(|_| tables_ab) {
        let s = Instant::now();
        let flat = if *name == "optimized" {
            prescribed(d, &d.cones, *bt, &None).and_then(|sol| {
                let note = format!("{:?}", sol.report.status);
                flatten(d, &d.cones, Some(&sol.log_scale_optimal)).map(|f| (f, note))
            })
        } else {
            flatten(d, &d.cones, None).map(|f| {
                let n = f.note.clone();
                (f, n)
            })
        };
        match flat {
            Err(e) => t.b.push(format!(
                "{:<20} {:<9} | FAILED after {:.1}s: {}",
                tag,
                name,
                s.elapsed().as_secs_f64(),
                short(e)
            )),
            Ok((f, note)) => {
                let m = measure(d, &f.uvs);
                let (mx, rms, corr, _) = compare(d, &m.u, &pm.u);
                // Our UVs scaled like the published ones and rounded to their digits.
                let scale = (pm.mean - m.mean).exp();
                let q = 10f64.powi(d.decimals as i32);
                let rounded: Vec<[f64; 2]> = f
                    .uvs
                    .iter()
                    .map(|p| {
                        [
                            (scale * p[0] * q).round() / q,
                            (scale * p[1] * q).round() / q,
                        ]
                    })
                    .collect();
                let mr = measure(d, &rounded);
                if rms < best_rms.0 {
                    best_rms = (rms, name);
                }
                let sh = error_shares(d, &m.u, &pm.u);
                let floor = compare(d, &mr.u, &m.u).1;
                t.b.push(format!(
                    "{:<20} {:<9} | {:>10.4e} {:>8.4e} {:>6.4} | {:>4} {:>5} {:>6} | {:>8.2e} {:>8.2e} {:>8.5} | {:>10.4e} rms {:>8.2e} | bdry {:>3.0}% of error ({:.1}% of area), cones {:>3.0}% ({:.2}%) | {} {:.1}s",
                    tag,
                    name,
                    m.dstar,
                    pm.dstar,
                    m.dstar / pm.dstar,
                    m.flipped,
                    m.degenerate,
                    f.charts,
                    mx,
                    rms,
                    corr,
                    mr.dstar,
                    floor,
                    100.0 * sh[0].0,
                    100.0 * sh[0].1,
                    100.0 * sh[1].0,
                    100.0 * sh[1].1,
                    note,
                    f.secs
                ));
            }
        }
    }
    if tables_ab {
        t.b.push(format!(
            "{:<20} {}| cone angles implied by the published UVs (2 pi - corner angle sum) against the file: largest difference {:.2e} rad",
            tag,
            if d.boundary {
                format!("best match to the published u: {} (rms {:.2e}) ", best_rms.1, best_rms.0)
            } else {
                String::new()
            },
            implied_angle_difference(d, &d.uvs)
        ));
    }

    // ---- Table C
    let omega = angle_defects(&d.mesh, &d.geometry).unwrap_or_default();
    let dw: Weights = if tips_variant {
        let mut w = vec![1.0; nv];
        for v in 0..nv {
            if !d.mesh.is_boundary_vertex(v) && omega[v] > 1.5 * PI {
                w[v] = 0.0;
            }
        }
        Some(w)
    } else {
        None
    };
    let mut label_base = if tips_variant {
        format!("{tag}[wE=0 at tips]")
    } else {
        tag.to_string()
    };
    if d.boundary {
        label_base.push_str(match treatment {
            Some(BoundaryTreatment::Free) => " [optimized]",
            Some(_) => " [u=0]",
            None => "",
        });
    }
    // Family: the boundary treatment whose re-optimized angles at the file's positions are
    // closer to the file's.
    let verts: Vec<u32> = d.cones.iter().map(|c| c.vertex).collect();
    let angle_diff = |bt| -> f64 {
        refit(d, &verts, bt, &None).map_or(f64::NAN, |s| {
            d.cones
                .iter()
                .map(|c| {
                    (s.cones
                        .iter()
                        .find(|x| x.vertex == c.vertex)
                        .map_or(0.0, |x| x.angle)
                        - c.angle)
                        .abs()
                })
                .fold(0.0, f64::max)
        })
    };
    let (family, fam_note) = if let Some(bt) = treatment {
        (
            bt,
            format!("boundary treatment {bt:?} (forced, the other family of this mesh)"),
        )
    } else if d.boundary {
        let (f0, f1) = (
            angle_diff(BoundaryTreatment::Fixed),
            angle_diff(BoundaryTreatment::Free),
        );
        let bt = if f1 < f0 || f0.is_nan() {
            BoundaryTreatment::Free
        } else {
            BoundaryTreatment::Fixed
        };
        (
            bt,
            format!(
                "boundary family {bt:?} (largest angle difference of the re-optimized file cones: u=0 {f0:.2e}, optimized {f1:.2e})"
            ),
        )
    } else {
        let f0 = angle_diff(BoundaryTreatment::Fixed);
        (
            BoundaryTreatment::Fixed,
            format!(
                "closed mesh (largest angle difference of the re-optimized file cones: {f0:.2e})"
            ),
        )
    };
    t.family = family;
    t.c.push(format!("{label_base}: {fam_note}"));
    // Published cones: solver distortion and published D*.
    let pub_sol = if family == BoundaryTreatment::Fixed {
        let problem = ConeProblem {
            distortion_weights: dw.clone(),
            ..ConeProblem::new(&d.mesh, &d.geometry, 0.0)
        };
        evaluate_cones(&problem, &d.cones, Parallelism::Sequential)
            .map(|e| (e.distortion, e.log_scale))
            .map_err(|e| e.to_string())
    } else {
        prescribed(d, &d.cones, family, &dw).map(|s| (s.distortion, s.log_scale_optimal))
    };
    match &pub_sol {
        Ok((dist, _)) => t.c.push(format!(
            "{:<30} | {:>5} {:>9.4} {:>12} | {:>8.5} | {:>11.4e} {:>11.4e} | {:>4} {:>5} | (published UVs)",
            format!("{label_base} published"),
            d.cones.len(),
            pub_total,
            "-",
            dist,
            pm.d,
            pm.dstar,
            pm.flipped,
            pm.degenerate
        )),
        Err(e) => t.c.push(format!("{label_base} published cones: solver FAILED {}", short(e.clone()))),
    }
    let cap = pub_total * (1.0 + total_slack());
    let sol_probe = |lambda: f64, s: Solution| Probe {
        lambda,
        total: total_angle(&s.cones),
        cones: s.cones,
        u: s.log_scale,
        dist: s.distortion,
        converged: s.report.status == Status::Converged,
        note: String::new(),
    };
    // C1
    let c1_eval = |l: f64| -> Result<Probe, String> {
        penalized(d, l, family, false, &dw).map(|s| sol_probe(l, s))
    };
    // C1 with the feedback loop: every probe is a loop of solves and flattenings.
    let c1fb_eval = |l: f64| -> Result<Probe, String> {
        let problem = ConeProblem {
            distortion_weights: dw.clone(),
            ..ConeProblem::new(&d.mesh, &d.geometry, l)
        };
        let result = solve_with_feedback(
            &problem,
            family,
            &options(false),
            &feedback_options(),
            |cones: &[Cone], u: &[f64]| {
                let f = flatten(d, cones, Some(u))?;
                Ok::<_, String>((realized_face_log_scale(&d.geometry, &f.uvs), ()))
            },
        )
        .map_err(|e| e.to_string())?;
        let report = &result.feedback;
        let rounds: Vec<String> = report
            .rounds
            .iter()
            .map(|r| {
                format!(
                    "{}max|d| {:.2e} (eligible {:.2e}), {} over, solve {:.2}s, flatten {:.2}s",
                    if r.probe { "probe: " } else { "" },
                    r.max_abs_mismatch,
                    r.max_abs_eligible,
                    r.n_over_threshold,
                    r.seconds_solve,
                    r.seconds_realize
                )
            })
            .collect();
        let note = format!(
            "[feedback: {} rounds, {} reweighted{}, {:?}, {:.2}s; {}]",
            report.rounds.len(),
            report.reweighted.len(),
            if report.reweighted.is_empty() {
                String::new()
            } else {
                format!(
                    " ({})",
                    report
                        .reweighted
                        .iter()
                        .take(32)
                        .map(u32::to_string)
                        .collect::<Vec<_>>()
                        .join(" ")
                )
            },
            report.stop,
            report.seconds,
            rounds.join(" | ")
        );
        let mut probe = sol_probe(l, result.solution);
        probe.note = note;
        Ok(probe)
    };
    let c2_eval = |l: f64| -> Result<Probe, String> {
        let s = penalized(d, l, family, true, &dw)?;
        let vs: Vec<u32> = s.cones.iter().map(|c| c.vertex).collect();
        if vs.is_empty() {
            return Ok(Probe {
                lambda: l,
                cones: vec![],
                u: vec![],
                total: 0.0,
                dist: s.distortion,
                converged: true,
                note: String::new(),
            });
        }
        refit(d, &vs, family, &dw).map(|r| sol_probe(l, r))
    };
    let pub_dstar = pm.dstar;
    let report = |name: &str,
                  s: &mut Search,
                  eval: &dyn Fn(f64) -> Result<Probe, String>,
                  t: &mut Tables| {
        let st = Instant::now();
        match &s.best {
            None => t.c.push(format!(
                "{:<30} | no feasible configuration found: {}",
                format!("{label_base} {name}"),
                s.notes.join("; ")
            )),
            Some(b) => {
                let note = if s.notes.is_empty() {
                    String::new()
                } else {
                    format!("[{}]", s.notes.join("; "))
                };
                let r = row_c(d, &format!("{label_base} {name}"), b, pub_dstar, &note);
                t.c.push(r);
            }
        }
        // equal number of cones
        let target = d.cones.len();
        let sb = Instant::now();
        let idx = equal_count(s, target, 10, budget, start, eval);
        match idx {
            None => t.c.push(format!(
                "{:<30} | equal count: no probe",
                format!("{label_base} {name} n={target}")
            )),
            Some(i) => {
                let p = &s.probes[i];
                if p.cones.is_empty() {
                    t.c.push(format!(
                        "{:<30} | equal count: closest probe has no cone",
                        format!("{label_base} {name} n={target}")
                    ));
                } else {
                    let mut note = if p.cones.len() == target {
                        String::new()
                    } else {
                        format!("[closest count attained: {} of {target}]", p.cones.len())
                    };
                    if p.total > cap {
                        note.push_str(&format!(
                            "[total angle {:.2}x the published: not an equal-angle comparison]",
                            p.total / pub_total
                        ));
                    }
                    let r = row_c(
                        d,
                        &format!("{label_base} {name} n={}", p.cones.len()),
                        p,
                        pub_dstar,
                        &note,
                    );
                    t.c.push(r);
                }
            }
        }
        let _ = (st, sb);
    };
    {
        if recipe("c1") {
            let mut s1 = search(cap, iters, grid, budget, start, c1_eval);
            report("C1 penalized", &mut s1, &c1_eval, &mut t);
        }
        // With the weights set by hand (the tips variant) the loop is not run.
        if recipe("c1fb") && !tips_variant {
            let mut s3 = search(cap, iters, grid, budget, start, c1fb_eval);
            report("C1 + feedback", &mut s3, &c1fb_eval, &mut t);
        }
        if recipe("c2") {
            let mut s2 = search(cap, iters, grid, budget, start, c2_eval);
            report("C2 merge+reopt", &mut s2, &c2_eval, &mut t);
        }
    }
    t.c.push(format!(
        "{label_base}: total {:.1}s",
        start.elapsed().as_secs_f64()
    ));
    t
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut names: Vec<String> = std::fs::read_dir(dataset_dir())
        .map(|rd| {
            rd.filter_map(|e| e.ok())
                .filter_map(|e| e.file_name().into_string().ok())
                .filter_map(|n| n.strip_suffix(".obj").map(str::to_owned))
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    if !args.is_empty() {
        names.retain(|n| args.contains(n));
    }
    if names.is_empty() {
        eprintln!("no dataset meshes found in {}", dataset_dir().display());
        return;
    }
    let mut all = Tables {
        a: vec![],
        b: vec![],
        c: vec![],
        family: BoundaryTreatment::Fixed,
    };
    for name in &names {
        eprintln!("== {name}");
        let d = match load(name) {
            Ok(d) => d,
            Err(e) => {
                all.a.push(format!("{name}: cannot load: {e}"));
                continue;
            }
        };
        let mut t = run_mesh(&d, false, None);
        let other = |bt| {
            if bt == BoundaryTreatment::Free {
                BoundaryTreatment::Fixed
            } else {
                BoundaryTreatment::Free
            }
        };
        if d.boundary {
            let t2 = run_mesh(&d, false, Some(other(t.family)));
            t.c.extend(t2.c);
        }
        if name == "camelshotnoise" {
            let t3 = run_mesh(&d, true, None);
            let fam = t3.family;
            t.c.extend(t3.c);
            let t4 = run_mesh(&d, true, Some(other(fam)));
            t.c.extend(t4.c);
        }
        all.a.extend(t.a);
        all.b.extend(t.b);
        all.c.extend(t.c);
    }
    println!("TABLE A: the published data as it is");
    println!("{A_HEAD}");
    for l in &all.a {
        println!("{l}");
    }
    println!(
        "\nTABLE B: our Boundary First Flattening with the file's cones against the published UVs"
    );
    println!("{B_HEAD}");
    for l in &all.b {
        println!("{l}");
    }
    println!(
        "\nTABLE C: our cones against the published cones at equal or smaller total cone angle"
    );
    println!(
        "(cap = published total * (1 + {}); 'equal' means D* within {:.1}%)",
        total_slack(),
        100.0 * EQUAL
    );
    println!("{C_HEAD}");
    for l in &all.c {
        println!("{l}");
    }
}
