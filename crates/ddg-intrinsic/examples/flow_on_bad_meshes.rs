//! FOLLOW-UP STUDY: the same measurements as `crates/ddg-variational-cuts/examples/bad_mesh_study.rs`
//! (commit 4b1985e of this repository; the mesh builders, the stage sets, the exact laws and the
//! tables are COPIED from it and only the parts marked "TREATMENT" are new), but on the input mesh
//! and on three intrinsic triangulations of it, selected by the environment variable `TREAT`:
//! `input` (default), `flip` (intrinsic Delaunay by flips; the export splits an edge only where
//! the result is not a simplicial complex), `r25`, `r30` (refined to 25 / 30 degrees, exported with
//! `ExportOptions::refine`). The flow runs unchanged on the export. For the sphere, the
//! intrinsic triangulation has the geometry of the POLYHEDRON; the exact law refers to the
//! smooth sphere; a vertex that is not an input vertex is placed on the polyhedron and its polar
//! angle is that of that point. Each treatment is compared with the input mesh, which has the
//! same polyhedron error.
//!
//! Original header of `bad_mesh_study.rs`:
//!
//! What obtuse and badly shaped triangles do to the cut flow: every stage and the whole flow
//! measured against the exact laws of the other studies, on meshes with the same geometry and a
//! worse triangulation (and on real meshes).
//!
//! ```text
//! cargo run --release -p ddg-variational-cuts --example bad_mesh_study -- quality
//! cargo run --release -p ddg-variational-cuts --example bad_mesh_study -- eval        # stage 1: evaluation
//! cargo run --release -p ddg-variational-cuts --example bad_mesh_study -- redist      # stage 2: redistancing
//! cargo run --release -p ddg-variational-cuts --example bad_mesh_study -- laws        # stage 3: the three laws
//! cargo run --release -p ddg-variational-cuts --example bad_mesh_study -- energy      # stage 4: energy history
//! cargo run --release -p ddg-variational-cuts --example bad_mesh_study -- attr        # one stage option at a time
//! cargo run --release -p ddg-variational-cuts --example bad_mesh_study -- real [ATTR] # real meshes
//! ```
//!
//! Environment: `FLAT_N` (cells per side of the flat meshes, default 32), `SPHERE_LEVEL`
//! (icosphere level, default 4), `STUDY_THREADS` (default 12), `LAW_STEPS` (steps of the circle
//! law, default 96), `REAL_ROOT` (the `data` directory with the models, default the main
//! checkout), `EXTRA_OBJ` (extra OBJ files for `real`, separated by `:`).
//!
//! Flat meshes are the unit square (sheared ones: the sheared lattice clipped to the square, so
//! the boundary is ragged), spheres are the unit sphere; all vertices of the sphere meshes lie
//! exactly on it. Three sets of stages are compared: `D` (the defaults), `P` (the paper as
//! printed, `FlowOptions::paper`: nodal normal derivative with sliver smoothing, length step on
//! the mesh, harmonic extension, `φ_t = s`) and `PS` (the structure of the paper with the
//! smoothed flux recovery, as in `redistance_study` and the `VARIANTS[0]` of `flow_study`).

// Index loops follow the formulas they implement; the task table and the loader return what
// they need.
#![allow(clippy::needless_range_loop, clippy::type_complexity)]
// Parts copied from the study this one follows are kept although nothing here calls them.
#![allow(dead_code)]

use std::collections::HashMap;
use std::f64::consts::PI;
use std::panic::AssertUnwindSafe;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use ddg_core::distance::{DistanceMode, upwind_distance_with};
use ddg_core::generators::{SplitMix64, icosphere};
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    EdgeTraceOptions, ExportOptions, IntrinsicTriangulation, RefineOptions, Segment,
    SimplicialExport, SurfacePoint,
};
use ddg_variational_cuts::{
    CurvatureSource, Cut, CutProblem, Energy, EvalOptions, Evaluation, Extension, Flow,
    FlowOptions, FluxRecovery, LengthStep, RegionWeights, SpeedScaling, StepReport, StepSize,
    StopRule, length_weight_for_ripples, region_of,
};

// ---------------------------------------------------------------------------------------------
// small helpers

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

/// Runs `f` on every item on a few threads; a panic becomes an `Err`.
fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<Result<R, String>> {
    let threads = env_usize("STUDY_THREADS", 12)
        .max(1)
        .min(items.len().max(1));
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

fn sci(x: f64) -> String {
    if x.is_nan() {
        "nan".into()
    } else {
        format!("{x:+.1e}")
    }
}

fn short(msg: &str) -> String {
    let one: String = msg.chars().take(110).collect();
    format!("ERR[{}]", one.replace('\n', " "))
}

fn flatten<T>(r: Result<Result<T, String>, String>) -> Result<T, String> {
    match r {
        Ok(Ok(x)) => Ok(x),
        Ok(Err(e)) | Err(e) => Err(e),
    }
}

// ---------------------------------------------------------------------------------------------
// mesh quality

#[derive(Clone, Copy, PartialEq, Debug)]
enum Kind {
    Flat,
    Sphere,
    Real,
}

#[derive(Clone, Debug)]
struct Quality {
    nv: usize,
    nf: usize,
    largest_median: f64,
    largest_p95: f64,
    largest_max: f64,
    smallest: f64,
    neg_cotan: f64,
    inverted: usize,
    degenerate: usize,
    h_median: f64,
    edge_ratio: f64,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn quality(positions: &[[f64; 3]], faces: &[[u32; 3]], kind: Kind) -> Quality {
    let mut largest = Vec::with_capacity(faces.len());
    let mut smallest = f64::INFINITY;
    let mut inverted = 0;
    let mut degenerate = 0;
    let mut edges: HashMap<(u32, u32), (f64, u8)> = HashMap::new();
    let mut lengths = Vec::new();
    for f in faces {
        let p = f.map(|v| positions[v as usize]);
        let mut a = [0.0; 3];
        for k in 0..3 {
            let u = sub(p[(k + 1) % 3], p[k]);
            let v = sub(p[(k + 2) % 3], p[k]);
            a[k] = norm(cross(u, v)).atan2(dot(u, v));
            if !a[k].is_finite() {
                a[k] = 0.0;
            }
        }
        let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
        if norm(n) == 0.0 {
            degenerate += 1;
        }
        let bad = match kind {
            Kind::Flat => n[2] <= 0.0,
            Kind::Sphere => {
                dot(
                    n,
                    [
                        p[0][0] + p[1][0] + p[2][0],
                        p[0][1] + p[1][1] + p[2][1],
                        p[0][2] + p[1][2] + p[2][2],
                    ],
                ) <= 0.0
            }
            Kind::Real => false,
        };
        if bad {
            inverted += 1;
        }
        largest.push(a.iter().cloned().fold(0.0, f64::max));
        smallest = smallest.min(a.iter().cloned().fold(f64::INFINITY, f64::min));
        for k in 0..3 {
            let (i, j) = (f[(k + 1) % 3], f[(k + 2) % 3]);
            let key = (i.min(j), i.max(j));
            let e = edges.entry(key).or_insert((0.0, 0));
            e.0 += a[k];
            e.1 += 1;
        }
    }
    let mut n_interior = 0;
    let mut neg = 0;
    for (&(i, j), &(sum, count)) in &edges {
        lengths.push(norm(sub(positions[i as usize], positions[j as usize])));
        if count == 2 {
            n_interior += 1;
            if sum > PI + 1e-12 {
                neg += 1;
            }
        }
    }
    largest.sort_by(f64::total_cmp);
    lengths.sort_by(f64::total_cmp);
    let at = |q: f64| -> f64 {
        if largest.is_empty() {
            f64::NAN
        } else {
            largest[((q * largest.len() as f64) as usize).min(largest.len() - 1)].to_degrees()
        }
    };
    Quality {
        nv: positions.len(),
        nf: faces.len(),
        largest_median: at(0.5),
        largest_p95: at(0.95),
        largest_max: at(1.0),
        smallest: smallest.to_degrees(),
        neg_cotan: neg as f64 / n_interior.max(1) as f64,
        inverted,
        degenerate,
        h_median: lengths.get(lengths.len() / 2).copied().unwrap_or(f64::NAN),
        edge_ratio: lengths.last().copied().unwrap_or(f64::NAN)
            / lengths.first().copied().unwrap_or(f64::NAN),
    }
}

fn quality_header() -> String {
    format!(
        "{:<34} {:>7} {:>7} | {:>6} {:>6} {:>6} | {:>6} | {:>8} | {:>5} {:>5}",
        "mesh", "verts", "faces", "med", "p95", "max", "min", "neg cot", "flip", "degen"
    )
}

fn quality_row(name: &str, q: &Quality) -> String {
    format!(
        "{:<34} {:>7} {:>7} | {:>6.1} {:>6.1} {:>6.1} | {:>6.2} | {:>7.2}% | {:>5} {:>5}",
        name,
        q.nv,
        q.nf,
        q.largest_median,
        q.largest_p95,
        q.largest_max,
        q.smallest,
        100.0 * q.neg_cotan,
        q.inverted,
        q.degenerate
    )
}

// ---------------------------------------------------------------------------------------------
// meshes

/// A mesh before it is validated: positions and faces, and for flat meshes the centre of the
/// test circle.
#[derive(Clone)]
struct Spec {
    name: String,
    kind: Kind,
    positions: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    center: [f64; 2],
}

/// TREATMENT: what the intrinsic triangulation cost and gave.
#[derive(Clone, Default, Debug)]
struct TreatInfo {
    tag: &'static str,
    /// Vertices added to the input mesh (refinement and splits).
    added: usize,
    n_input: usize,
    /// Seconds: `IntrinsicTriangulation::new`, flips, refinement, export (its own flips,
    /// splits and refinements included).
    t_new: f64,
    t_flip: f64,
    t_refine: f64,
    t_export: f64,
    n_splits: usize,
    n_export_flips: usize,
    export_min: f64,
    export_max: f64,
}

struct Shape {
    spec: Spec,
    mesh: TriMesh,
    geometry: Geometry,
    q: Quality,
    info: TreatInfo,
}

impl Spec {
    fn quality(&self) -> Quality {
        quality(&self.positions, &self.faces, self.kind)
    }
    fn build(&self) -> Result<Shape, String> {
        let mesh = TriMesh::new(self.positions.len(), self.faces.clone())
            .map_err(|e| format!("TriMesh::new: {e}"))?;
        let geometry = Geometry::from_positions(&mesh, &self.positions)
            .map_err(|e| format!("Geometry::from_positions: {e}"))?;
        let q = quality_geom(&mesh, &geometry);
        Ok(Shape {
            spec: self.clone(),
            mesh,
            geometry,
            q,
            info: TreatInfo {
                tag: "in",
                n_input: self.positions.len(),
                ..TreatInfo::default()
            },
        })
    }
}

impl Shape {
    fn positions(&self) -> &[[f64; 3]] {
        &self.spec.positions
    }
}

/// Rescales a flat spec to total area 1 (about the origin).
fn unit_area(mut spec: Spec) -> Spec {
    let area: f64 = spec
        .faces
        .iter()
        .map(|f| {
            let p = f.map(|v| spec.positions[v as usize]);
            0.5 * cross(sub(p[1], p[0]), sub(p[2], p[0]))[2]
        })
        .sum();
    let s = 1.0 / area.sqrt();
    for p in &mut spec.positions {
        for x in p.iter_mut() {
            *x *= s;
        }
    }
    spec.center = [spec.center[0] * s, spec.center[1] * s];
    spec
}

/// Quad lattice `(nx+1) × (ny+1)`, vertex `(i, j)` at `pos(i, j)`, each cell split along the
/// diagonal `(i, j)–(i+1, j+1)`.
fn lattice(
    nx: usize,
    ny: usize,
    pos: impl Fn(usize, usize) -> [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let id = |i: usize, j: usize| (i * (ny + 1) + j) as u32;
    let mut positions = Vec::new();
    for i in 0..=nx {
        for j in 0..=ny {
            positions.push(pos(i, j));
        }
    }
    let mut faces = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    (positions, faces)
}

fn jittered_grid(n: usize, jitter: f64, seed: u64, name: &str) -> Spec {
    let h = 1.0 / n as f64;
    let mut rng = SplitMix64::new(seed);
    let (mut positions, faces) = lattice(n, n, |i, j| [i as f64 * h, j as f64 * h, 0.0]);
    for p in &mut positions {
        let inside = |x: f64| x > 0.5 * h && x < 1.0 - 0.5 * h;
        if inside(p[0]) && inside(p[1]) {
            p[0] += jitter * h * rng.uniform(-1.0, 1.0);
            p[1] += jitter * h * rng.uniform(-1.0, 1.0);
        }
    }
    Spec {
        name: name.into(),
        kind: Kind::Flat,
        positions,
        faces,
        center: [0.511, 0.493],
    }
}

fn anisotropic(nx: usize, ratio: usize, name: &str) -> Spec {
    let ny = nx * ratio;
    let (positions, faces) = lattice(nx, ny, |i, j| {
        [i as f64 / nx as f64, j as f64 / ny as f64, 0.0]
    });
    Spec {
        name: name.into(),
        kind: Kind::Flat,
        positions,
        faces,
        center: [0.511, 0.493],
    }
}

/// Removes faces until no vertex is a pinch point (its faces form one fan).
fn prune_pinches(mut faces: Vec<[u32; 3]>) -> Vec<[u32; 3]> {
    loop {
        let nf = faces.len();
        let mut parent: Vec<usize> = (0..3 * nf).collect();
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        let mut directed: HashMap<(u32, u32), usize> = HashMap::new();
        for (f, face) in faces.iter().enumerate() {
            for k in 0..3 {
                directed.insert((face[k], face[(k + 1) % 3]), f);
            }
        }
        let corner = |faces: &[[u32; 3]], f: usize, v: u32| -> usize {
            3 * f + faces[f].iter().position(|&x| x == v).unwrap()
        };
        for f in 0..nf {
            for k in 0..3 {
                let (a, b) = (faces[f][k], faces[f][(k + 1) % 3]);
                if let Some(&g) = directed.get(&(b, a)) {
                    for v in [a, b] {
                        let (x, y) = (corner(&faces, f, v), corner(&faces, g, v));
                        let (rx, ry) = (find(&mut parent, x), find(&mut parent, y));
                        parent[rx] = ry;
                    }
                }
            }
        }
        // Components of corners per vertex.
        let mut comps: HashMap<u32, HashMap<usize, usize>> = HashMap::new();
        for f in 0..nf {
            for k in 0..3 {
                let root = find(&mut parent, 3 * f + k);
                *comps
                    .entry(faces[f][k])
                    .or_default()
                    .entry(root)
                    .or_insert(0) += 1;
            }
        }
        let mut drop = vec![false; nf];
        let mut any = false;
        for (&v, c) in &comps {
            if c.len() > 1 {
                any = true;
                let best = *c
                    .iter()
                    .max_by_key(|(r, n)| (**n, std::cmp::Reverse(**r)))
                    .unwrap()
                    .0;
                for f in 0..nf {
                    if let Some(k) = faces[f].iter().position(|&x| x == v)
                        && find(&mut parent, 3 * f + k) != best
                    {
                        drop[f] = true;
                    }
                }
            }
        }
        if !any {
            return faces;
        }
        faces = faces
            .into_iter()
            .zip(drop)
            .filter(|(_, d)| !d)
            .map(|(f, _)| f)
            .collect();
    }
}

fn compact(positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut map = vec![u32::MAX; positions.len()];
    let mut out = Vec::new();
    let faces: Vec<[u32; 3]> = faces
        .iter()
        .map(|f| {
            f.map(|v| {
                if map[v as usize] == u32::MAX {
                    map[v as usize] = out.len() as u32;
                    out.push(positions[v as usize]);
                }
                map[v as usize]
            })
        })
        .collect();
    (out, faces)
}

/// The lattice `x = (i + σ j) h, y = j h` (the unit-cell lattice sheared by `σ`) clipped to the
/// unit square: every triangle has the largest angle `180° − atan(1/σ)`, the boundary is ragged.
fn sheared(n: usize, sigma: f64, name: &str) -> Spec {
    let h = 1.0 / n as f64;
    let mut index: HashMap<(i64, i64), u32> = HashMap::new();
    let mut positions = Vec::new();
    let eps = 1e-9;
    for j in 0..=n as i64 {
        let lo = (-sigma * j as f64).floor() as i64 - 1;
        let hi = (n as f64 - sigma * j as f64).ceil() as i64 + 1;
        for i in lo..=hi {
            let x = (i as f64 + sigma * j as f64) * h;
            if x >= -eps && x <= 1.0 + eps {
                index.insert((i, j), positions.len() as u32);
                positions.push([x, j as f64 * h, 0.0]);
            }
        }
    }
    let mut faces = Vec::new();
    for (&(i, j), &a) in &index {
        let (b, c, d) = (
            index.get(&(i + 1, j)),
            index.get(&(i + 1, j + 1)),
            index.get(&(i, j + 1)),
        );
        if let (Some(&b), Some(&c), Some(&d)) = (b, c, d) {
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    faces.sort();
    let faces = prune_pinches(faces);
    let (positions, faces) = compact(positions, faces);
    unit_area(Spec {
        name: name.into(),
        kind: Kind::Flat,
        positions,
        faces,
        center: [0.511, 0.493],
    })
}

/// Rows alternate between thin strips of height `0.05 d` (caps: one angle of 168°) and tall ones
/// of height `3 d` (needles: apex 19°), the rows shifted by half a cell.
fn needles_and_caps(n: usize, name: &str) -> Spec {
    let d = 1.0 / n as f64;
    let (ha, hb) = (0.05 * d, 3.0 * d);
    let pairs = ((1.0 / (ha + hb)).round() as usize).max(1);
    let scale = 1.0 / (pairs as f64 * (ha + hb));
    let rows = 2 * pairs + 1;
    let mut y = vec![0.0f64];
    for j in 0..2 * pairs {
        let last = y[j];
        y.push(last + scale * if j % 2 == 0 { ha } else { hb });
    }
    *y.last_mut().unwrap() = 1.0;
    let mut positions = Vec::new();
    let mut row_x: Vec<Vec<f64>> = Vec::new();
    let mut start = Vec::new();
    for j in 0..rows {
        let xs: Vec<f64> = if j % 2 == 0 {
            (0..=n).map(|i| i as f64 * d).collect()
        } else {
            std::iter::once(0.0)
                .chain((0..n).map(|i| (i as f64 + 0.5) * d))
                .chain(std::iter::once(1.0))
                .collect()
        };
        start.push(positions.len() as u32);
        for &x in &xs {
            positions.push([x, y[j], 0.0]);
        }
        row_x.push(xs);
    }
    let mut faces = Vec::new();
    for j in 0..rows - 1 {
        let (xa, xb) = (&row_x[j], &row_x[j + 1]);
        let (mut i, mut m) = (0usize, 0usize);
        while i < xa.len() - 1 || m < xb.len() - 1 {
            let bottom = if i == xa.len() - 1 {
                false
            } else if m == xb.len() - 1 {
                true
            } else {
                xa[i + 1] <= xb[m + 1]
            };
            let (a, b) = (start[j] + i as u32, start[j + 1] + m as u32);
            if bottom {
                faces.push([a, a + 1, b]);
                i += 1;
            } else {
                faces.push([a, b + 1, b]);
                m += 1;
            }
        }
    }
    Spec {
        name: name.into(),
        kind: Kind::Flat,
        positions,
        faces,
        center: [0.511, 0.493],
    }
}

fn ico(level: u32, amplitude: f64, name: &str) -> Spec {
    let g = icosphere(level).unwrap();
    let mut positions = g.positions.clone();
    if amplitude > 0.0 {
        let mut lengths = Vec::new();
        for f in &g.faces {
            for k in 0..3 {
                lengths.push(norm(sub(
                    g.positions[f[k] as usize],
                    g.positions[f[(k + 1) % 3] as usize],
                )));
            }
        }
        let h = lengths.iter().sum::<f64>() / lengths.len() as f64;
        let mut rng = SplitMix64::new(11);
        for p in &mut positions {
            let n = norm(*p);
            let nrm = [p[0] / n, p[1] / n, p[2] / n];
            let a = if nrm[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            let t1 = cross(nrm, a);
            let l1 = norm(t1);
            let t1 = [t1[0] / l1, t1[1] / l1, t1[2] / l1];
            let t2 = cross(nrm, t1);
            let (u, v) = (
                amplitude * h * rng.uniform(-1.0, 1.0),
                amplitude * h * rng.uniform(-1.0, 1.0),
            );
            let q = [
                nrm[0] + u * t1[0] + v * t2[0],
                nrm[1] + u * t1[1] + v * t2[1],
                nrm[2] + u * t1[2] + v * t2[2],
            ];
            let l = norm(q);
            *p = [q[0] / l, q[1] / l, q[2] / l];
        }
    }
    Spec {
        name: name.into(),
        kind: Kind::Sphere,
        positions,
        faces: g.faces,
        center: [0.0, 0.0],
    }
}

/// `rings − 1` rings of `n_lon` vertices between two pole vertices.
fn latlon(rings: usize, n_lon: usize, name: &str) -> Spec {
    let mut positions = vec![[0.0, 0.0, 1.0]];
    for k in 1..rings {
        let th = PI * k as f64 / rings as f64;
        for j in 0..n_lon {
            let ph = 2.0 * PI * j as f64 / n_lon as f64;
            positions.push([th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()]);
        }
    }
    positions.push([0.0, 0.0, -1.0]);
    let south = (positions.len() - 1) as u32;
    let ring = |k: usize, j: usize| (1 + (k - 1) * n_lon + j % n_lon) as u32;
    let mut faces = Vec::new();
    for j in 0..n_lon {
        faces.push([0, ring(1, j), ring(1, j + 1)]);
        faces.push([south, ring(rings - 1, j + 1), ring(rings - 1, j)]);
    }
    for k in 1..rings - 1 {
        for j in 0..n_lon {
            let (a, b, c, d) = (
                ring(k, j),
                ring(k, j + 1),
                ring(k + 1, j + 1),
                ring(k + 1, j),
            );
            faces.push([a, d, c]);
            faces.push([a, c, b]);
        }
    }
    Spec {
        name: name.into(),
        kind: Kind::Sphere,
        positions,
        faces,
        center: [0.0, 0.0],
    }
}

/// Flips edges whose opposite angles sum to more than 180 degrees until none is left (the
/// Delaunay triangulation of the same vertices: in the plane with the boundary edges fixed, on the
/// sphere the convex hull). A flip is made only if both new triangles keep their orientation.
/// Returns the faces and the number of flips.
fn delaunay(spec: &Spec) -> (Spec, usize) {
    let pos = &spec.positions;
    let mut faces = spec.faces.clone();
    let mut directed: HashMap<(u32, u32), usize> = HashMap::new();
    for (f, face) in faces.iter().enumerate() {
        for k in 0..3 {
            directed.insert((face[k], face[(k + 1) % 3]), f);
        }
    }
    let angle_at = |v: u32, a: u32, b: u32| -> f64 {
        let (p, q, r) = (pos[v as usize], pos[a as usize], pos[b as usize]);
        let (u, w) = (sub(q, p), sub(r, p));
        norm(cross(u, w)).atan2(dot(u, w))
    };
    let oriented = |a: u32, b: u32, c: u32| -> bool {
        let (p, q, r) = (pos[a as usize], pos[b as usize], pos[c as usize]);
        let n = cross(sub(q, p), sub(r, p));
        match spec.kind {
            Kind::Flat => n[2] > 1e-14,
            _ => {
                dot(
                    n,
                    [p[0] + q[0] + r[0], p[1] + q[1] + r[1], p[2] + q[2] + r[2]],
                ) > 1e-14
            }
        }
    };
    let mut stack: Vec<(u32, u32)> = directed.keys().copied().collect();
    stack.sort();
    let mut flips = 0;
    while let Some((a, b)) = stack.pop() {
        let (Some(&f), Some(&g)) = (directed.get(&(a, b)), directed.get(&(b, a))) else {
            continue;
        };
        let c = *faces[f].iter().find(|&&v| v != a && v != b).unwrap();
        let d = *faces[g].iter().find(|&&v| v != a && v != b).unwrap();
        if c == d || directed.contains_key(&(c, d)) || directed.contains_key(&(d, c)) {
            continue;
        }
        if angle_at(c, a, b) + angle_at(d, a, b) <= PI + 1e-9 {
            continue;
        }
        if !(oriented(a, d, c) && oriented(d, b, c)) {
            continue;
        }
        for e in [(a, b), (b, c), (c, a), (b, a), (a, d), (d, b)] {
            directed.remove(&e);
        }
        faces[f] = [a, d, c];
        faces[g] = [d, b, c];
        for (face, idx) in [([a, d, c], f), ([d, b, c], g)] {
            for k in 0..3 {
                directed.insert((face[k], face[(k + 1) % 3]), idx);
            }
        }
        flips += 1;
        stack.extend([(a, d), (d, b), (b, c), (c, a)]);
    }
    let mut out = spec.clone();
    out.faces = faces;
    out.name = format!("{} [Del]", spec.name);
    (out, flips)
}

fn maybe_delaunay(specs: Vec<Spec>) -> Vec<Spec> {
    if std::env::var("DELAUNAY").is_ok() {
        specs.iter().map(|s| delaunay(s).0).collect()
    } else {
        specs
    }
}

fn flat_family() -> Vec<Spec> {
    let n = env_usize("FLAT_N", 32);
    maybe_delaunay(vec![
        jittered_grid(n, 0.0, 5, &format!("flat: regular grid {n}^2")),
        jittered_grid(n, 0.3, 5, "flat: jitter 0.30 (baseline)"),
        jittered_grid(n, 0.33, 5, "flat: jitter 0.33"),
        sheared(n, 0.5, "flat: sheared 117"),
        sheared(n, 1.0, "flat: sheared 135"),
        sheared(n, 2.0, "flat: sheared 153"),
        sheared(n, 1.0 / (15.0f64).to_radians().tan(), "flat: sheared 165"),
        anisotropic(n, 4, "flat: cells 1:4"),
        anisotropic(n, 16, "flat: cells 1:16"),
        needles_and_caps(n, "flat: needles and caps"),
    ])
}

fn sphere_family() -> Vec<Spec> {
    let l = env_usize("SPHERE_LEVEL", 4) as u32;
    maybe_delaunay(vec![
        ico(l, 0.0, &format!("sphere: icosphere({l})")),
        ico(l, 0.2, "sphere: ico tangential 0.20 h"),
        ico(l, 0.3, "sphere: ico tangential 0.30 h"),
        ico(l, 0.35, "sphere: ico tangential 0.35 h"),
        latlon(41, 82, "sphere: latlon 41x82"),
        latlon(41, 246, "sphere: latlon 41x246"),
        latlon(81, 162, "sphere: latlon 81x162 (finer)"),
    ])
}

// ---------------------------------------------------------------------------------------------
// TREATMENT: the intrinsic triangulation and its export as a `Shape`

#[derive(Clone, Copy, PartialEq, Debug)]
enum Treat {
    Input,
    Flip,
    R25,
    R30,
}

impl Treat {
    fn from_env() -> Treat {
        match std::env::var("TREAT").as_deref() {
            Ok("flip") => Treat::Flip,
            Ok("r25") => Treat::R25,
            Ok("r30") => Treat::R30,
            _ => Treat::Input,
        }
    }
    fn tag(self) -> &'static str {
        match self {
            Treat::Input => "in",
            Treat::Flip => "flip",
            Treat::R25 => "r25",
            Treat::R30 => "r30",
        }
    }
    fn refine(self) -> Option<RefineOptions> {
        let degrees = match self {
            Treat::R25 => 25.0f64,
            Treat::R30 => 30.0,
            _ => return None,
        };
        Some(RefineOptions {
            min_angle: degrees.to_radians(),
            ..RefineOptions::default()
        })
    }
}

/// The position in space of a point of the input mesh (copied from `ddg-intrinsic/tests/common`).
fn position_on_input(positions: &[[f64; 3]], mesh: &TriMesh, point: SurfacePoint) -> [f64; 3] {
    match point {
        SurfacePoint::Vertex(v) => positions[v as usize],
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = mesh.edge_vertices(edge as usize);
            [0, 1, 2].map(|k| (1.0 - t) * positions[a][k] + t * positions[b][k])
        }
        SurfacePoint::Face { face, bary } => {
            let f = mesh.faces()[face as usize];
            [0, 1, 2].map(|k| (0..3).map(|c| bary[c] * positions[f[c] as usize][k]).sum())
        }
    }
}

/// Quality numbers from the angles of a geometry (the intrinsic ones for an export).
fn quality_geom(mesh: &TriMesh, geometry: &Geometry) -> Quality {
    let angles = geometry.corner_angles();
    let mut largest = Vec::new();
    let mut smallest = f64::INFINITY;
    let mut edges: HashMap<(u32, u32), (f64, u8)> = HashMap::new();
    for (f, face) in mesh.faces().iter().enumerate() {
        let a = angles[f];
        largest.push(a.iter().cloned().fold(0.0, f64::max));
        smallest = smallest.min(a.iter().cloned().fold(f64::INFINITY, f64::min));
        for k in 0..3 {
            let (i, j) = (face[(k + 1) % 3], face[(k + 2) % 3]);
            let e = edges.entry((i.min(j), i.max(j))).or_insert((0.0, 0));
            e.0 += a[k];
            e.1 += 1;
        }
    }
    let (mut n_interior, mut neg) = (0, 0);
    for &(sum, count) in edges.values() {
        if count == 2 {
            n_interior += 1;
            if sum > PI + 1e-12 {
                neg += 1;
            }
        }
    }
    largest.sort_by(f64::total_cmp);
    let at =
        |q: f64| largest[((q * largest.len() as f64) as usize).min(largest.len() - 1)].to_degrees();
    let mut lengths = geometry.edge_lengths().to_vec();
    lengths.sort_by(f64::total_cmp);
    Quality {
        nv: mesh.n_vertices(),
        nf: mesh.n_faces(),
        largest_median: at(0.5),
        largest_p95: at(0.95),
        largest_max: at(1.0),
        smallest: smallest.to_degrees(),
        neg_cotan: neg as f64 / n_interior.max(1) as f64,
        inverted: 0,
        degenerate: 0,
        h_median: lengths[lengths.len() / 2],
        edge_ratio: lengths[lengths.len() - 1] / lengths[0],
    }
}

/// Everything a treatment produces, kept for the way back.
struct Treated {
    shape: Shape,
    tri: Option<IntrinsicTriangulation>,
    export: Option<SimplicialExport>,
}

fn treat_full(spec: &Spec, treat: Treat) -> Result<Treated, String> {
    let input = spec.build()?;
    if treat == Treat::Input {
        return Ok(Treated {
            shape: input,
            tri: None,
            export: None,
        });
    }
    let mut info = TreatInfo {
        tag: treat.tag(),
        n_input: spec.positions.len(),
        ..TreatInfo::default()
    };
    let start = Instant::now();
    let mut tri = IntrinsicTriangulation::new(&input.mesh, &input.geometry)
        .map_err(|e| format!("IntrinsicTriangulation::new: {e}"))?;
    info.t_new = start.elapsed().as_secs_f64();
    let export = match treat.refine() {
        None => {
            let start = Instant::now();
            tri.flip_to_delaunay(&Default::default())
                .map_err(|e| format!("flip_to_delaunay: {e}"))?;
            info.t_flip = start.elapsed().as_secs_f64();
            let start = Instant::now();
            let e = tri
                .export_simplicial(&ExportOptions::default())
                .map_err(|e| format!("export_simplicial: {e}"))?;
            info.t_export = start.elapsed().as_secs_f64();
            e
        }
        Some(options) => {
            let start = Instant::now();
            tri.refine(&options).map_err(|e| format!("refine: {e}"))?;
            info.t_refine = start.elapsed().as_secs_f64();
            let start = Instant::now();
            let e = tri
                .export_simplicial(&ExportOptions {
                    refine: Some(options),
                    ..ExportOptions::default()
                })
                .map_err(|e| format!("export_simplicial: {e}"))?;
            info.t_export = start.elapsed().as_secs_f64();
            e
        }
    };
    info.added = export.mesh.n_vertices() - spec.positions.len();
    info.n_splits = export.n_splits;
    info.n_export_flips = export.n_flips;
    info.export_min = export.min_angle.to_degrees();
    info.export_max = export.max_angle.to_degrees();
    let positions: Vec<[f64; 3]> = export
        .locations
        .iter()
        .map(|&l| position_on_input(&spec.positions, &input.mesh, l))
        .collect();
    let mut new_spec = spec.clone();
    new_spec.name = format!("{} [{}]", spec.name, treat.tag());
    new_spec.positions = positions;
    new_spec.faces = export.mesh.faces().to_vec();
    let q = quality_geom(&export.mesh, &export.geometry);
    let shape = Shape {
        spec: new_spec,
        mesh: export.mesh.clone(),
        geometry: export.geometry.clone(),
        q,
        info,
    };
    Ok(Treated {
        shape,
        tri: Some(tri),
        export: Some(export),
    })
}

/// A spec as a `Shape` under the treatment of the environment.
fn make(spec: &Spec) -> Result<Shape, String> {
    treat_full(spec, Treat::from_env()).map(|t| t.shape)
}

fn build_all(specs: &[Spec]) -> Vec<(Spec, Result<Shape, String>)> {
    // Independent meshes are treated in parallel.
    let res = par_map(specs, make);
    specs
        .iter()
        .zip(res)
        .map(|(s, r)| {
            let r = r.and_then(|x| x);
            let mut s = s.clone();
            if let Ok(sh) = &r {
                s.name = sh.spec.name.clone();
            }
            (s, r)
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// sets of stages and the stepping driver

#[derive(Clone, Copy, PartialEq, Debug)]
enum Redist {
    /// The library's (`Redistance::KeepAdjacent`, single-pass fast marching).
    Keep,
    /// None.
    Off,
    /// Cut-adjacent values kept, the rest by the fixed point of the same scheme (done here).
    FixedPoint,
    /// Cut-adjacent values kept, every other vertex gets the exact distance to the circle /
    /// latitude through the mean position of the cut (done here; ideal far field).
    Oracle,
}

#[derive(Clone, Copy)]
enum Oracle {
    None,
    Circle([f64; 2]),
    Latitude,
}

#[derive(Clone)]
struct Stages {
    name: String,
    flux: FluxRecovery,
    passes: usize,
    length: LengthStep,
    ext: Extension,
    speed: SpeedScaling,
    redist: Redist,
}

impl Stages {
    fn d() -> Stages {
        Stages {
            name: "D".into(),
            flux: FluxRecovery::default(),
            passes: 0,
            length: LengthStep::OnCut,
            ext: Extension::Local,
            speed: SpeedScaling::MeanGradient,
            redist: Redist::Keep,
        }
    }
    fn p() -> Stages {
        Stages {
            name: "P".into(),
            flux: FluxRecovery::Nodal,
            passes: 1,
            length: LengthStep::OnMesh,
            ext: Extension::Harmonic,
            speed: SpeedScaling::Unit,
            redist: Redist::Keep,
        }
    }
    fn ps() -> Stages {
        Stages {
            name: "PS".into(),
            flux: FluxRecovery::default(),
            passes: 0,
            ..Stages::p()
        }
    }
    fn with(mut self, name: &str, f: impl FnOnce(&mut Stages)) -> Stages {
        f(&mut self);
        self.name = name.into();
        self
    }
    fn options(&self, step: StepSize) -> FlowOptions<'static> {
        let mut o = FlowOptions::new(step);
        o.eval.flux_recovery = self.flux;
        o.eval.smoothing_passes = self.passes;
        o.length_step = self.length;
        o.extension = self.ext;
        o.speed_scaling = self.speed;
        o.redistance = if self.redist == Redist::Keep {
            ddg_variational_cuts::Redistance::KeepAdjacent
        } else {
            ddg_variational_cuts::Redistance::Off
        };
        o
    }
    fn eval_options(&self) -> EvalOptions<'static> {
        self.options(StepSize::default()).eval
    }
}

fn paper_check() {
    // `P` must be `FlowOptions::paper`.
    let a = Stages::p().options(StepSize::default());
    let b = FlowOptions::paper(StepSize::default());
    assert_eq!(a.eval.flux_recovery, b.eval.flux_recovery);
    assert_eq!(a.eval.smoothing_passes, b.eval.smoothing_passes);
    assert_eq!(a.length_step, b.length_step);
    assert_eq!(a.extension, b.extension);
    assert_eq!(a.speed_scaling, b.speed_scaling);
    assert_eq!(a.redistance, b.redistance);
    let d = Stages::d().options(StepSize::default());
    let e = FlowOptions::new(StepSize::default());
    assert_eq!(d.eval.flux_recovery, e.eval.flux_recovery);
    assert_eq!(d.length_step, e.length_step);
    assert_eq!(d.extension, e.extension);
    assert_eq!(d.speed_scaling, e.speed_scaling);
}

fn polar(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1]).sqrt().atan2(p[2])
}

fn cut_points(sh: &Shape, phi: &[f64]) -> Vec<[f64; 3]> {
    match Cut::from_level_set(&sh.mesh, phi, 1e-3) {
        Ok(cut) => cut
            .point_positions(&sh.mesh, sh.positions())
            .unwrap_or_default(),
        Err(_) => Vec::new(),
    }
}

/// Cut-adjacent vertices keep their values, the others get the fixed point of the upwind scheme.
fn redistance_fixed_point(sh: &Shape, phi: &[f64]) -> Result<(Vec<f64>, usize, bool), String> {
    let cut = Cut::from_level_set(&sh.mesh, phi, 1e-3).map_err(|e| e.to_string())?;
    if cut.points.is_empty() {
        return Ok((phi.to_vec(), 0, true));
    }
    let adjacent = cut.cut_adjacent_vertices(&sh.mesh);
    let region: Vec<u32> = phi.iter().map(|x| region_of(*x)).collect();
    let seeds: Vec<(u32, f64)> = (0..phi.len())
        .filter(|&v| adjacent[v])
        .map(|v| (v as u32, phi[v].abs()))
        .collect();
    let r = upwind_distance_with(
        &sh.mesh,
        &sh.geometry,
        &seeds,
        Some(&region),
        &DistanceMode::fixed_point(phi.len()),
    )
    .map_err(|e| e.to_string())?;
    let out = (0..phi.len())
        .map(|v| {
            if adjacent[v] || !r.values[v].is_finite() {
                phi[v]
            } else if region[v] == 0 {
                -r.values[v]
            } else {
                r.values[v]
            }
        })
        .collect();
    let (reopened, converged) = r
        .fixed_point
        .map_or((0, true), |f| (f.n_reopened, f.converged));
    Ok((out, reopened, converged))
}

fn oracle_phi(sh: &Shape, phi: &[f64], o: Oracle) -> Option<Vec<f64>> {
    let pts = cut_points(sh, phi);
    if pts.is_empty() {
        return None;
    }
    match o {
        Oracle::None => None,
        Oracle::Circle(c) => {
            let r = |p: [f64; 3]| ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt();
            let mean = pts.iter().map(|p| r(*p)).sum::<f64>() / pts.len() as f64;
            Some(sh.positions().iter().map(|p| r(*p) - mean).collect())
        }
        Oracle::Latitude => {
            let mean = pts.iter().map(|p| polar(*p)).sum::<f64>() / pts.len() as f64;
            Some(sh.positions().iter().map(|p| polar(*p) - mean).collect())
        }
    }
}

fn do_step(
    flow: &mut Flow,
    sh: &Shape,
    problem: &CutProblem,
    st: &Stages,
    step: StepSize,
    oracle: Oracle,
) -> Result<StepReport, String> {
    let r = flow
        .step(problem, &st.options(step))
        .map_err(|e| e.to_string())?;
    match st.redist {
        Redist::Keep | Redist::Off => {}
        Redist::FixedPoint => {
            let (phi, _, _) = redistance_fixed_point(sh, &flow.level_set())?;
            flow.set_level_set(&phi).map_err(|e| e.to_string())?;
        }
        Redist::Oracle => {
            // The cut-adjacent vertices keep their values, every other vertex gets the exact
            // distance to the circle (or latitude) through the mean position of the cut points.
            let now = flow.level_set();
            if let Some(target) = oracle_phi(sh, &now, oracle) {
                let adjacent = Cut::from_level_set(&sh.mesh, &now, 1e-3)
                    .map_err(|e| e.to_string())?
                    .cut_adjacent_vertices(&sh.mesh);
                let phi: Vec<f64> = (0..now.len())
                    .map(|v| if adjacent[v] { now[v] } else { target[v] })
                    .collect();
                flow.set_level_set(&phi).map_err(|e| e.to_string())?;
            }
        }
    }
    Ok(r)
}

struct RunOut {
    steps: Vec<StepReport>,
    stop: String,
    final_energy: Option<f64>,
    seconds: f64,
    returned_step: Option<usize>,
}

/// `Flow::run` of the crate (the new stop rule: stalled / rose, set back to the lowest energy
/// state). `final_energy` is the energy of the last state before a set-back.
fn run_loop(
    flow: &mut Flow,
    _sh: &Shape,
    problem: &CutProblem,
    st: &Stages,
    step: StepSize,
    stop: &StopRule,
    _oracle: Oracle,
) -> RunOut {
    assert!(matches!(st.redist, Redist::Keep | Redist::Off));
    let start = Instant::now();
    match flow.run(problem, &st.options(step), stop) {
        Ok(r) => RunOut {
            steps: r.steps,
            stop: format!("{:?}", r.stop),
            final_energy: Some(r.last_energy),
            seconds: start.elapsed().as_secs_f64(),
            returned_step: Some(r.returned_step),
        },
        Err(e) => RunOut {
            steps: Vec::new(),
            stop: short(&e.to_string()),
            final_energy: None,
            seconds: start.elapsed().as_secs_f64(),
            returned_step: None,
        },
    }
}

// ---------------------------------------------------------------------------------------------
// stage 1: evaluation only

#[derive(Default, Debug, Clone)]
struct EvalErr {
    u_max: f64,
    d_max: [f64; 2],
    d_mean: [f64; 2],
    s_max: f64,
    s_mean: f64,
    j: f64,
    n_cut: usize,
    n_sliver: usize,
    smallest_offset: f64,
    clamped: usize,
}

fn rel_stats(values: impl Iterator<Item = f64>, exact: f64) -> (f64, f64) {
    let (mut max, mut sum, mut n) = (0.0f64, 0.0, 0usize);
    for v in values {
        let e = (v / exact - 1.0).abs();
        max = max.max(if e.is_nan() { f64::INFINITY } else { e });
        sum += e;
        n += 1;
    }
    (max, sum / n.max(1) as f64)
}

fn per_side(e: &Evaluation, r: usize) -> Vec<f64> {
    (0..e.cut.points.len())
        .map(|k| e.normal_derivative[e.triangles.point_vertices_by_region(k)[r]])
        .collect()
}

/// Sphere cap of polar angle 1: `u`, `∂u/∂n` on both sides, `s`, `J_0`.
fn cap_eval(sh: &Shape, opts: &EvalOptions<'_>) -> Result<EvalErr, String> {
    let psi0 = 1.0f64;
    let r = 1.0 / (4.0 * PI).sqrt();
    let mut ev =
        ddg_variational_cuts::Evaluator::new(&sh.mesh, &sh.geometry).map_err(|e| e.to_string())?;
    let phi: Vec<f64> = sh.positions().iter().map(|p| polar(*p) - psi0).collect();
    let problem = CutProblem::new(Energy::Dirichlet, 0.0);
    let e = ev
        .evaluate(&phi, &problem, opts)
        .map_err(|e| e.to_string())?;
    let exact_d = [(0.5 * psi0).tan() / r, 1.0 / ((0.5 * psi0).tan() * r)];
    let exact_s = exact_d[0] * exact_d[0] - exact_d[1] * exact_d[1];
    let mut u_max = 0.0f64;
    for (v, &p) in sh.positions().iter().enumerate() {
        let psi = polar(p);
        let exact = if e.cut.vertex_region[v] == 0 {
            2.0 * ((0.5 * psi0).cos() / (0.5 * psi).cos()).ln()
        } else {
            2.0 * ((0.5 * psi0).sin() / (0.5 * psi).sin()).ln()
        };
        u_max = u_max.max((e.log_scale[v] - exact).abs());
    }
    let exact_j = -8.0 * PI * (0.5 * psi0).cos().ln() - 4.0 * PI * (0.5 * psi0).sin().powi(2);
    let s0 = rel_stats(per_side(&e, 0).into_iter(), exact_d[0]);
    let s1 = rel_stats(per_side(&e, 1).into_iter(), exact_d[1]);
    let s = rel_stats(e.speed.iter().copied(), exact_s);
    Ok(EvalErr {
        u_max,
        d_max: [s0.0, s1.0],
        d_mean: [s0.1, s1.1],
        s_max: s.0,
        s_mean: s.1,
        j: (e.energy.distortion[0] / exact_j - 1.0).abs(),
        n_cut: e.cut.points.len(),
        n_sliver: e.report.n_sliver_points,
        smallest_offset: e.report.smallest_crossing_offset,
        clamped: e.report.n_clamped_crossings,
    })
}

/// Flat disk of radius 0.3 with the constant source `K = 5` in region 0 only.
fn disk_eval(sh: &Shape, opts: &EvalOptions<'_>) -> Result<EvalErr, String> {
    let (k, rho) = (5.0f64, 0.3f64);
    let c = sh.spec.center;
    let radius = |p: [f64; 3]| ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt();
    let mut ev =
        ddg_variational_cuts::Evaluator::new(&sh.mesh, &sh.geometry).map_err(|e| e.to_string())?;
    let phi: Vec<f64> = sh.positions().iter().map(|p| radius(*p) - rho).collect();
    let problem = CutProblem {
        curvature: CurvatureSource::Uniform(k),
        regions: vec![
            RegionWeights {
                distortion: 1.0,
                area: 0.0,
            },
            RegionWeights {
                distortion: 0.0,
                area: 0.0,
            },
        ],
        ..CutProblem::new(Energy::Dirichlet, 0.0)
    };
    let e = ev
        .evaluate(&phi, &problem, opts)
        .map_err(|e| e.to_string())?;
    let mut u_max = 0.0f64;
    for (v, &p) in sh.positions().iter().enumerate() {
        if e.cut.vertex_region[v] == 0 {
            let exact = 0.25 * k * (radius(p).powi(2) - rho * rho);
            u_max = u_max.max((e.log_scale[v] - exact).abs());
        }
    }
    let exact_d = 0.5 * k * rho;
    let s0 = rel_stats(per_side(&e, 0).into_iter(), exact_d);
    let s = rel_stats(e.speed.iter().copied(), exact_d * exact_d);
    let exact_j = PI * k * k * rho.powi(4) / 8.0;
    Ok(EvalErr {
        u_max,
        d_max: [s0.0, f64::NAN],
        d_mean: [s0.1, f64::NAN],
        s_max: s.0,
        s_mean: s.1,
        j: (e.energy.distortion[0] / exact_j - 1.0).abs(),
        n_cut: e.cut.points.len(),
        n_sliver: e.report.n_sliver_points,
        smallest_offset: e.report.smallest_crossing_offset,
        clamped: e.report.n_clamped_crossings,
    })
}

fn eval_study() {
    let sets = [
        ("D smoothed c=.5", EvalOptions::default()),
        (
            "nodal",
            EvalOptions {
                flux_recovery: FluxRecovery::Nodal,
                ..EvalOptions::default()
            },
        ),
        ("P nodal+sliver", EvalOptions::paper()),
    ];
    println!(
        "STAGE 1: evaluation only against exact solutions. Relative errors (max / mean) of du/dn and of s, max |u-u*|, relative error of J_0."
    );
    println!(
        "Sphere cap psi0 = 1 (du/dn on region 0 / region 1 sides); flat disk rho = 0.3, K = 5 in region 0 (du/dn of region 0)."
    );
    println!(
        "columns: u_max | du/dn max,mean (region 0) | du/dn max,mean (region 1) | s max,mean | J0 | cut pts, sliver pts, smallest offset, clamped\n"
    );
    for (title, specs, cap) in [
        ("SPHERE CAP", sphere_family(), true),
        ("FLAT DISK", flat_family(), false),
    ] {
        println!("{title}");
        let built = build_all(&specs);
        let jobs: Vec<(usize, usize)> = (0..built.len())
            .flat_map(|m| (0..sets.len()).map(move |s| (m, s)))
            .collect();
        let res = par_map(&jobs, |&(m, s)| match &built[m].1 {
            Ok(sh) => {
                if cap {
                    cap_eval(sh, &sets[s].1)
                } else {
                    disk_eval(sh, &sets[s].1)
                }
            }
            Err(e) => Err(e.clone()),
        });
        for (m, (spec, _)) in built.iter().enumerate() {
            for s in 0..sets.len() {
                let r = flatten(res[m * sets.len() + s].clone());
                match r {
                    Ok(x) => println!(
                        "{:<32} {:<16} u {:.1e} | d0 {:.1e}, {:.1e} | d1 {:.1e}, {:.1e} | s {:.1e}, {:.1e} | J0 {:.1e} | {} {} {:.1e} {}",
                        spec.name,
                        sets[s].0,
                        x.u_max,
                        x.d_max[0],
                        x.d_mean[0],
                        x.d_max[1],
                        x.d_mean[1],
                        x.s_max,
                        x.s_mean,
                        x.j,
                        x.n_cut,
                        x.n_sliver,
                        x.smallest_offset,
                        x.clamped
                    ),
                    Err(e) => println!("{:<32} {:<16} {}", spec.name, sets[s].0, short(&e)),
                }
            }
        }
        println!();
    }
}

// ---------------------------------------------------------------------------------------------
// stage 2: redistancing only

#[derive(Clone)]
struct RedistOut {
    abs_mean: f64,
    abs_max: f64,
    rel_mean: f64,
    rel_max: f64,
    fp_abs_mean: f64,
    fp_abs_max: f64,
    fp_rel_max: f64,
    fp_reopened: usize,
    fp_converged: bool,
    unreached: usize,
    crossed: (f64, f64, f64),
    near: (f64, f64, f64),
}

/// `|grad phi|` of the linear interpolant in face `fi` from the intrinsic geometry:
/// `|grad|^2 A = 1/2 sum_k cot_k (phi_i - phi_j)^2` (as `init::redistanced` of the crate).
fn face_gradient_geometry(sh: &Shape, fi: usize, v: [f64; 3]) -> f64 {
    let cot = sh.geometry.corner_cotans()[fi];
    let s: f64 = (0..3)
        .map(|c| cot[c] * (v[(c + 1) % 3] - v[(c + 2) % 3]).powi(2))
        .sum();
    (s / (2.0 * sh.geometry.face_areas()[fi])).max(0.0).sqrt()
}

#[allow(dead_code)]
fn face_gradient(p: [[f64; 3]; 3], f: [f64; 3]) -> f64 {
    let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
    let l = norm(n);
    if l == 0.0 {
        return f64::NAN;
    }
    let (ea, eb, ec) = (sub(p[2], p[1]), sub(p[0], p[2]), sub(p[1], p[0]));
    let v = [
        f[0] * ea[0] + f[1] * eb[0] + f[2] * ec[0],
        f[0] * ea[1] + f[1] * eb[1] + f[2] * ec[1],
        f[0] * ea[2] + f[1] * eb[2] + f[2] * ec[2],
    ];
    let nn = [n[0] / l, n[1] / l, n[2] / l];
    norm(cross(nn, v)) / l
}

fn redist_test(sh: &Shape, exact: &dyn Fn([f64; 3]) -> f64) -> Result<RedistOut, String> {
    let pos = sh.positions();
    let d: Vec<f64> = pos.iter().map(|p| exact(*p)).collect();
    let mut flow = Flow::new(&sh.mesh, &sh.geometry, &d).map_err(|e| e.to_string())?;
    let unreached = flow.redistance().map_err(|e| e.to_string())?;
    let phi = flow.level_set();
    let h = (sh.geometry.total_area() / pos.len() as f64).sqrt();
    let adjacent = Cut::from_level_set(&sh.mesh, &d, 1e-3)
        .map_err(|e| e.to_string())?
        .cut_adjacent_vertices(&sh.mesh);
    let stats = |phi: &[f64]| -> (f64, f64, f64, f64) {
        let (mut sa, mut ma, mut sr, mut mr, mut na, mut nr) =
            (0.0, 0.0f64, 0.0, 0.0f64, 0usize, 0usize);
        for v in 0..phi.len() {
            if adjacent[v] || d[v].abs() > 0.25 {
                continue;
            }
            let e = (phi[v] - d[v]).abs();
            sa += e / h;
            ma = ma.max(e / h);
            na += 1;
            if d[v].abs() > 3.0 * h {
                sr += e / d[v].abs();
                mr = mr.max(e / d[v].abs());
                nr += 1;
            }
        }
        (sa / na.max(1) as f64, ma, sr / nr.max(1) as f64, mr)
    };
    let (abs_mean, abs_max, rel_mean, rel_max) = stats(&phi);
    let (fp, reopened, converged) = redistance_fixed_point(sh, &d)?;
    let (fa_mean, fa_max, _, fr_max) = stats(&fp);
    // |grad phi| of the linear interpolant after the redistancing.
    let mut crossed = (f64::INFINITY, 0.0, 0.0f64);
    let mut near = (f64::INFINITY, 0.0, 0.0f64);
    let (mut nc, mut nn) = (0usize, 0usize);
    for (fi, f) in sh.mesh.faces().iter().enumerate() {
        let _ = pos;
        let vals = f.map(|v| phi[v as usize]);
        let g = face_gradient_geometry(sh, fi, vals);
        if !g.is_finite() {
            continue;
        }
        let regions = vals.map(region_of);
        let mixed = regions[0] != regions[1] || regions[1] != regions[2];
        let closest = vals.iter().map(|x| x.abs()).fold(f64::INFINITY, f64::min);
        if mixed {
            crossed = (crossed.0.min(g), crossed.1 + g, crossed.2.max(g));
            nc += 1;
        }
        if closest <= 2.0 * h {
            near = (near.0.min(g), near.1 + g, near.2.max(g));
            nn += 1;
        }
    }
    Ok(RedistOut {
        abs_mean,
        abs_max,
        rel_mean,
        rel_max,
        fp_abs_mean: fa_mean,
        fp_abs_max: fa_max,
        fp_rel_max: fr_max,
        fp_reopened: reopened,
        fp_converged: converged,
        unreached,
        crossed: (crossed.0, crossed.1 / nc.max(1) as f64, crossed.2),
        near: (near.0, near.1 / nn.max(1) as f64, near.2),
    })
}

fn redist_study() {
    println!(
        "STAGE 2: redistancing only. Start: the exact signed distance; the cut-adjacent vertices keep it, the rest is recomputed by the single-pass fast marching of the flow (SP) and, for comparison, by the fixed point of the same scheme (FP)."
    );
    println!(
        "abs = |phi - d| / h with h = sqrt(area / vertices), over non-adjacent vertices with |d| < 0.25; rel = |phi - d| / |d| for |d| > 3 h; grad = |grad phi| of the linear interpolant in crossed faces (min/mean/max) and in faces within 2 h of the cut (near).\n"
    );
    for (title, specs, flat) in [
        ("FLAT", flat_family(), true),
        ("SPHERE", sphere_family(), false),
    ] {
        let built = build_all(&specs);
        let names: [&str; 2] = if flat {
            ["line", "circle"]
        } else {
            ["equator", "latitude 1.0"]
        };
        let exact_for = |s: &Shape, which: usize| -> Box<dyn Fn([f64; 3]) -> f64> {
            let c = s.spec.center;
            match (flat, which) {
                (true, 0) => Box::new(move |p: [f64; 3]| p[0] - c[0]),
                (true, _) => Box::new(move |p: [f64; 3]| {
                    ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt() - 0.3
                }),
                (false, 0) => Box::new(|p: [f64; 3]| polar(p) - 0.5 * PI),
                (false, _) => Box::new(|p: [f64; 3]| polar(p) - 1.0),
            }
        };
        println!("{title}");
        println!(
            "{:<32} {:<12} | SP abs mean,max  rel mean,max | FP abs mean,max rel max (reopened, conv) | unreached | grad crossed min/mean/max | near min/mean/max",
            "mesh", "shape"
        );
        let jobs: Vec<(usize, usize)> = (0..built.len())
            .flat_map(|m| (0..2).map(move |s| (m, s)))
            .collect();
        let res = par_map(&jobs, |&(m, s)| match &built[m].1 {
            Ok(sh) => {
                let f = exact_for(sh, s);
                redist_test(sh, &*f)
            }
            Err(e) => Err(e.clone()),
        });
        for (k, &(m, s)) in jobs.iter().enumerate() {
            let name = &built[m].0.name;
            match flatten(res[k].clone()) {
                Ok(r) => println!(
                    "{:<32} {:<12} | {:.3} {:.3}  {:.1e} {:.1e} | {:.3} {:.3} {:.1e} ({}, {}) | {} | {:.3}/{:.3}/{:.3} | {:.3}/{:.3}/{:.3}",
                    name,
                    names[s],
                    r.abs_mean,
                    r.abs_max,
                    r.rel_mean,
                    r.rel_max,
                    r.fp_abs_mean,
                    r.fp_abs_max,
                    r.fp_rel_max,
                    r.fp_reopened,
                    r.fp_converged,
                    r.unreached,
                    r.crossed.0,
                    r.crossed.1,
                    r.crossed.2,
                    r.near.0,
                    r.near.1,
                    r.near.2
                ),
                Err(e) => println!("{:<32} {:<12} {}", name, names[s], short(&e)),
            }
        }
        println!();
    }
}

// ---------------------------------------------------------------------------------------------
// stage 3: the laws

#[derive(Clone, Debug)]
struct CircleOut {
    err: f64,
    round: f64,
    g: (f64, f64, f64),
    unreached: usize,
    clamped: usize,
}

fn circle_law(sh: &Shape, st: &Stages, steps: usize) -> Result<CircleOut, String> {
    let (rho0, alpha, t_end) = (0.3f64, 1.0f64, 0.03f64);
    let c = sh.spec.center;
    let exact = (rho0 * rho0 - 2.0 * alpha * t_end).sqrt();
    let radius = |p: [f64; 3]| ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt();
    let phi0: Vec<f64> = sh.positions().iter().map(|p| radius(*p) - rho0).collect();
    let mut flow = Flow::new(&sh.mesh, &sh.geometry, &phi0).map_err(|e| e.to_string())?;
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    let tau = t_end / steps as f64;
    let (mut unreached, mut clamped) = (0, 0);
    let mut last = None;
    for _ in 0..steps {
        let r = do_step(
            &mut flow,
            sh,
            &problem,
            st,
            StepSize::Fixed(tau),
            Oracle::Circle(c),
        )?;
        unreached = unreached.max(r.n_unreached);
        clamped = clamped.max(r.n_clamped_crossings);
        last = Some(r);
    }
    let pts = cut_points(sh, &flow.level_set());
    if pts.is_empty() {
        return Err("the cut vanished".into());
    }
    let rs: Vec<f64> = pts.iter().map(|p| radius(*p)).collect();
    let mean = rs.iter().sum::<f64>() / rs.len() as f64;
    let std = (rs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / rs.len() as f64).sqrt();
    let g = last
        .and_then(|r| r.gradient_norm)
        .map_or((f64::NAN, f64::NAN, f64::NAN), |g| (g.min, g.mean, g.max));
    Ok(CircleOut {
        err: mean / exact - 1.0,
        round: std / mean,
        g,
        unreached,
        clamped,
    })
}

/// Flat disk with the uniform source `K` in region 0 only (Dirichlet, `w = (1, 0)`):
/// `ρ' = −(Kρ/2)² − α_L/ρ`. Returns (relative error of the mean radius, roundness, last
/// |grad phi| min/mean/max).
fn disk_law(sh: &Shape, st: &Stages, steps: usize) -> Result<CircleOut, String> {
    let (k, rho0, alpha, t_end) = (5.0f64, 0.3f64, 0.3f64, 0.05f64);
    let rhs = |r: f64| -(0.5 * k * r).powi(2) - alpha / r;
    let n = 20000;
    let h = t_end / n as f64;
    let mut rho = rho0;
    for _ in 0..n {
        let k1 = rhs(rho);
        let k2 = rhs(rho + 0.5 * h * k1);
        let k3 = rhs(rho + 0.5 * h * k2);
        let k4 = rhs(rho + h * k3);
        rho += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    let c = sh.spec.center;
    let radius = |p: [f64; 3]| ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt();
    let phi0: Vec<f64> = sh.positions().iter().map(|p| radius(*p) - rho0).collect();
    let mut flow = Flow::new(&sh.mesh, &sh.geometry, &phi0).map_err(|e| e.to_string())?;
    let problem = CutProblem {
        curvature: CurvatureSource::Uniform(k),
        regions: vec![
            RegionWeights {
                distortion: 1.0,
                area: 0.0,
            },
            RegionWeights {
                distortion: 0.0,
                area: 0.0,
            },
        ],
        ..CutProblem::new(Energy::Dirichlet, alpha)
    };
    let (mut unreached, mut clamped) = (0, 0);
    let mut last = None;
    for _ in 0..steps {
        let r = do_step(
            &mut flow,
            sh,
            &problem,
            st,
            StepSize::Fixed(t_end / steps as f64),
            Oracle::Circle(c),
        )?;
        unreached = unreached.max(r.n_unreached);
        clamped = clamped.max(r.n_clamped_crossings);
        last = Some(r);
    }
    let pts = cut_points(sh, &flow.level_set());
    if pts.is_empty() {
        return Err("the cut vanished".into());
    }
    let rs: Vec<f64> = pts.iter().map(|p| radius(*p)).collect();
    let (mean, std) = mean_std(&rs);
    let g = last
        .and_then(|r| r.gradient_norm)
        .map_or((f64::NAN, f64::NAN, f64::NAN), |g| (g.min, g.mean, g.max));
    Ok(CircleOut {
        err: mean / rho - 1.0,
        round: std / mean,
        g,
        unreached,
        clamped,
    })
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

fn mean_std(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    (
        mean,
        (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n).sqrt(),
    )
}

/// Latitude circle against the ODE: (error of the mean latitude, spread).
fn lat_ode(
    sh: &Shape,
    st: &Stages,
    alpha: f64,
    psi0: f64,
    t_end: f64,
    steps: usize,
) -> Result<(f64, f64), String> {
    let target = latitude_exact(psi0, alpha, t_end);
    let phi0: Vec<f64> = sh.positions().iter().map(|p| polar(*p) - psi0).collect();
    let mut flow = Flow::new(&sh.mesh, &sh.geometry, &phi0).map_err(|e| e.to_string())?;
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    for _ in 0..steps {
        do_step(
            &mut flow,
            sh,
            &problem,
            st,
            StepSize::Fixed(t_end / steps as f64),
            Oracle::Latitude,
        )?;
    }
    let pts = cut_points(sh, &flow.level_set());
    if pts.is_empty() {
        return Err("the cut vanished".into());
    }
    let (mean, std) = mean_std(&pts.iter().map(|p| polar(*p)).collect::<Vec<_>>());
    Ok((mean - target, std))
}

/// Latitude circle with `alpha = 20` run to rest with the adaptive step: (final latitude minus
/// the balanced latitude, spread, steps, stop).
fn balance(sh: &Shape, st: &Stages, psi0: f64) -> Result<(f64, f64, usize, String), String> {
    let r = 1.0 / (4.0 * PI).sqrt();
    let alpha = 20.0;
    let star = (4.0 / (alpha * r)).asin();
    let phi0: Vec<f64> = sh.positions().iter().map(|p| polar(*p) - psi0).collect();
    let mut flow = Flow::new(&sh.mesh, &sh.geometry, &phi0).map_err(|e| e.to_string())?;
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    let stop = StopRule {
        max_steps: 3000,
        motion_tolerance: 1e-7,
        energy_tolerance: 0.0,
        energy_window: 1,
    };
    let out = run_loop(
        &mut flow,
        sh,
        &problem,
        st,
        StepSize::default(),
        &stop,
        Oracle::Latitude,
    );
    let pts = cut_points(sh, &flow.level_set());
    if pts.is_empty() {
        return Err(format!(
            "the cut vanished after {} steps ({})",
            out.steps.len(),
            out.stop
        ));
    }
    let (mean, std) = mean_std(&pts.iter().map(|p| polar(*p)).collect::<Vec<_>>());
    Ok((mean - star, std, out.steps.len(), out.stop))
}

fn mode_amplitude(pts: &[[f64; 3]], m: usize) -> (f64, f64) {
    let mut samples: Vec<(f64, f64)> = pts.iter().map(|p| (p[1].atan2(p[0]), polar(*p))).collect();
    samples.sort_by(|a, b| a.0.total_cmp(&b.0));
    let n = samples.len();
    if n < 3 {
        return (f64::NAN, f64::NAN);
    }
    let weight = |k: usize| {
        let before = samples[(k + n - 1) % n].0;
        let after = samples[(k + 1) % n].0;
        0.5 * (after - before).rem_euclid(2.0 * PI) / (2.0 * PI)
    };
    let mean: f64 = (0..n).map(|k| weight(k) * samples[k].1).sum();
    let amp: f64 = (0..n)
        .map(|k| 2.0 * weight(k) * (samples[k].1 - mean) * (m as f64 * samples[k].0).cos())
        .sum();
    let rest: f64 = (0..n)
        .map(|k| {
            let r = samples[k].1 - mean - amp * (m as f64 * samples[k].0).cos();
            weight(k) * r * r
        })
        .sum();
    (amp, rest.sqrt())
}

/// Growth rate of the amplitude of the ripple `ε cos(mθ)` of the equator: (measured, exact, rest).
fn ripple(sh: &Shape, st: &Stages, alpha: f64, m: usize) -> Result<(f64, f64, f64), String> {
    let r = 1.0 / (4.0 * PI).sqrt();
    let exact = (m as f64 - 1.0) * (4.0 / r.powi(3) - alpha * (m as f64 + 1.0) / (r * r));
    let horizon = (0.7 / exact.abs().max(30.0)).min(4e-3);
    let steps = 160;
    let eps = 0.02;
    let phi0: Vec<f64> = sh
        .positions()
        .iter()
        .map(|p| polar(*p) - 0.5 * PI - eps * (m as f64 * p[1].atan2(p[0])).cos())
        .collect();
    let mut flow = Flow::new(&sh.mesh, &sh.geometry, &phi0).map_err(|e| e.to_string())?;
    flow.redistance().map_err(|e| e.to_string())?;
    let (a0, _) = mode_amplitude(&cut_points(sh, &flow.level_set()), m);
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    for _ in 0..steps {
        do_step(
            &mut flow,
            sh,
            &problem,
            st,
            StepSize::Fixed(horizon / steps as f64),
            Oracle::None,
        )?;
    }
    let (a1, rest) = mode_amplitude(&cut_points(sh, &flow.level_set()), m);
    Ok(((a1 / a0).ln() / horizon, exact, rest))
}

fn ripple_cases() -> Vec<(f64, usize)> {
    vec![(0.0, 4), (2.0, 2), (2.0, 4), (6.0, 4), (10.0, 4)]
}

fn laws_study() {
    let steps = env_usize("LAW_STEPS", 96);
    let sets = [Stages::d(), Stages::p(), Stages::ps()];
    println!(
        "STAGE 3: the laws, for the stage sets D (defaults), P (paper as printed), PS (paper structure, smoothed flux)."
    );
    println!("A cell shows the error; ERR[...] is a failed run.\n");
    // Circle law.
    println!(
        "CIRCLE LAW (flat, rho0 = 0.3, alpha_L = 1, t = 0.03, {steps} steps): relative error of the mean radius / roundness (std/mean) / |grad phi| min,mean,max in the crossed faces at the last step / unreached, clamped (max over steps)"
    );
    let flat = build_all(&flat_family());
    let jobs: Vec<(usize, usize)> = (0..flat.len())
        .flat_map(|m| (0..sets.len()).map(move |s| (m, s)))
        .collect();
    let res = par_map(&jobs, |&(m, s)| match &flat[m].1 {
        Ok(sh) => circle_law(sh, &sets[s], steps),
        Err(e) => Err(e.clone()),
    });
    for (k, &(m, s)) in jobs.iter().enumerate() {
        match flatten(res[k].clone()) {
            Ok(r) => println!(
                "{:<32} {:<3} {} / {:.1e} / {:.2}-{:.2}-{:.2} / {} {}",
                flat[m].0.name,
                sets[s].name,
                sci(r.err),
                r.round,
                r.g.0,
                r.g.1,
                r.g.2,
                r.unreached,
                r.clamped
            ),
            Err(e) => println!("{:<32} {:<3} {}", flat[m].0.name, sets[s].name, short(&e)),
        }
    }
    println!();
    println!(
        "DISK LAW (flat, uniform source K = 5 in region 0 only, Dirichlet, alpha_L = 0.3, rho0 = 0.3, t = 0.05, {steps} steps, rho' = -(K rho/2)^2 - alpha/rho): same columns"
    );
    let res = par_map(&jobs, |&(m, s)| match &flat[m].1 {
        Ok(sh) => disk_law(sh, &sets[s], steps),
        Err(e) => Err(e.clone()),
    });
    for (k, &(m, s)) in jobs.iter().enumerate() {
        match flatten(res[k].clone()) {
            Ok(r) => println!(
                "{:<32} {:<3} {} / {:.1e} / {:.2}-{:.2}-{:.2} / {} {}",
                flat[m].0.name,
                sets[s].name,
                sci(r.err),
                r.round,
                r.g.0,
                r.g.1,
                r.g.2,
                r.unreached,
                r.clamped
            ),
            Err(e) => println!("{:<32} {:<3} {}", flat[m].0.name, sets[s].name, short(&e)),
        }
    }
    println!();
    let sph = build_all(&sphere_family());
    // Latitude ODE.
    let lat_cases = [(20.0, 0.5, 4e-3, 200usize), (8.0, 0.8, 4e-3, 200)];
    println!(
        "LATITUDE CIRCLE vs ODE (sphere): error of the mean latitude (spread); cases (alpha_L, psi0, t, steps)"
    );
    let jobs: Vec<(usize, usize, usize)> = (0..sph.len())
        .flat_map(|m| {
            (0..sets.len()).flat_map(move |s| (0..lat_cases.len()).map(move |c| (m, s, c)))
        })
        .collect();
    let res = par_map(&jobs, |&(m, s, c)| match &sph[m].1 {
        Ok(sh) => lat_ode(
            sh,
            &sets[s],
            lat_cases[c].0,
            lat_cases[c].1,
            lat_cases[c].2,
            lat_cases[c].3,
        ),
        Err(e) => Err(e.clone()),
    });
    for (k, &(m, s, c)) in jobs.iter().enumerate() {
        match flatten(res[k].clone()) {
            Ok((e, sp)) => println!(
                "{:<32} {:<3} case {:?}: {} ({:.0e})",
                sph[m].0.name,
                sets[s].name,
                lat_cases[c],
                sci(e),
                sp
            ),
            Err(e) => println!(
                "{:<32} {:<3} case {:?}: {}",
                sph[m].0.name,
                sets[s].name,
                lat_cases[c],
                short(&e)
            ),
        }
    }
    println!();
    // Balanced latitude.
    let r = 1.0 / (4.0 * PI).sqrt();
    println!(
        "BALANCED LATITUDE (alpha_L = 20, sin psi* = 4/(alpha R), psi* = {:.5}); run to rest with the adaptive step (motion <= 1e-7, <= 3000 steps): final psi - psi* (spread), steps, stop",
        (4.0 / (20.0 * r)).asin()
    );
    let starts = [0.5, 1.1];
    let jobs: Vec<(usize, usize, usize)> = (0..sph.len())
        .flat_map(|m| (0..sets.len()).flat_map(move |s| (0..2).map(move |c| (m, s, c))))
        .collect();
    let res = par_map(&jobs, |&(m, s, c)| match &sph[m].1 {
        Ok(sh) => balance(sh, &sets[s], starts[c]),
        Err(e) => Err(e.clone()),
    });
    for (k, &(m, s, c)) in jobs.iter().enumerate() {
        match flatten(res[k].clone()) {
            Ok((e, sp, n, stop)) => println!(
                "{:<32} {:<3} from {}: {} ({:.0e}) {} steps {}",
                sph[m].0.name,
                sets[s].name,
                starts[c],
                sci(e),
                sp,
                n,
                stop
            ),
            Err(e) => println!(
                "{:<32} {:<3} from {}: {}",
                sph[m].0.name,
                sets[s].name,
                starts[c],
                short(&e)
            ),
        }
    }
    println!();
    // Ripples.
    let cases = ripple_cases();
    println!(
        "RIPPLES OF THE EQUATOR: measured growth rate / exact - 1 (and the measured rate); cases (alpha_L, m)"
    );
    let jobs: Vec<(usize, usize, usize)> = (0..sph.len())
        .flat_map(|m| {
            (0..sets.len()).flat_map({
                let n = cases.len();
                move |s| (0..n).map(move |c| (m, s, c))
            })
        })
        .collect();
    let res = par_map(&jobs, |&(m, s, c)| match &sph[m].1 {
        Ok(sh) => ripple(sh, &sets[s], cases[c].0, cases[c].1),
        Err(e) => Err(e.clone()),
    });
    for (k, &(m, s, c)) in jobs.iter().enumerate() {
        match flatten(res[k].clone()) {
            Ok((meas, exact, rest)) => println!(
                "{:<32} {:<3} (alpha {}, m {}): {:+.3} (measured {:+.1}, exact {:+.1}, rest {:.0e})",
                sph[m].0.name,
                sets[s].name,
                cases[c].0,
                cases[c].1,
                meas / exact - 1.0,
                meas,
                exact,
                rest
            ),
            Err(e) => println!(
                "{:<32} {:<3} (alpha {}, m {}): {}",
                sph[m].0.name,
                sets[s].name,
                cases[c].0,
                cases[c].1,
                short(&e)
            ),
        }
    }
}

// ---------------------------------------------------------------------------------------------
// stage 4: the energy history

struct EnergySummary {
    steps: usize,
    stop: String,
    e0: f64,
    e_end: f64,
    rises: usize,
    max_rise: f64,
    max_rise_frac: f64,
    unreached_max: usize,
    clamped_mean: f64,
    clamped_max: usize,
    floor: f64,
    secs_per_step: f64,
    motion_last: f64,
}

fn summarize(out: &RunOut, flow_scale: f64, mean_edge: f64) -> EnergySummary {
    let mut e: Vec<f64> = out.steps.iter().map(|r| r.energy.total).collect();
    if let Some(f) = out.final_energy {
        e.push(f);
    }
    let (e0, e_end) = (
        e.first().copied().unwrap_or(f64::NAN),
        e.last().copied().unwrap_or(f64::NAN),
    );
    let drop = e0 - e_end;
    let rises: Vec<f64> = e
        .windows(2)
        .map(|w| w[1] - w[0])
        .filter(|d| *d > 0.0)
        .collect();
    let max_rise = rises.iter().cloned().fold(0.0, f64::max);
    let n = out.steps.len();
    let tail = (n / 4).clamp(1, 50).min(n);
    let mut motions: Vec<f64> = out.steps[n - tail..].iter().map(|r| r.cut_motion).collect();
    motions.sort_by(f64::total_cmp);
    let floor = if motions.is_empty() {
        f64::NAN
    } else {
        motions[motions.len() / 2]
    } / (flow_scale * mean_edge);
    EnergySummary {
        steps: n,
        stop: out.stop.clone(),
        e0,
        e_end,
        rises: rises.len(),
        max_rise,
        max_rise_frac: if drop > 0.0 {
            max_rise / drop
        } else {
            f64::NAN
        },
        unreached_max: out.steps.iter().map(|r| r.n_unreached).max().unwrap_or(0),
        clamped_mean: out
            .steps
            .iter()
            .map(|r| r.n_clamped_crossings as f64)
            .sum::<f64>()
            / n.max(1) as f64,
        clamped_max: out
            .steps
            .iter()
            .map(|r| r.n_clamped_crossings)
            .max()
            .unwrap_or(0),
        floor,
        secs_per_step: out.seconds / n.max(1) as f64,
        motion_last: out.steps.last().map_or(f64::NAN, |r| r.cut_motion),
    }
}

fn mean_edge(sh: &Shape) -> f64 {
    let l = sh.geometry.edge_lengths();
    l.iter().sum::<f64>() / l.len() as f64
}

fn energy_line(name: &str, st: &str, s: &EnergySummary) -> String {
    format!(
        "{:<32} {:<3} steps {:>4} stop {:<14} E {:.5e} -> {:.5e} | rises {:>3} (max {:.1e} = {:.1e} of drop) | unreached max {} | clamped mean {:.1} max {} | floor {:.2e} edge | last motion {:.1e} | {:.1} ms/step",
        name,
        st,
        s.steps,
        s.stop,
        s.e0,
        s.e_end,
        s.rises,
        s.max_rise,
        s.max_rise_frac,
        s.unreached_max,
        s.clamped_mean,
        s.clamped_max,
        s.floor,
        s.motion_last,
        1e3 * s.secs_per_step
    )
}

/// Flat: a bump of height 0.15 in the middle, Hencky, a straight cut at the left, small `α_L`.
fn flat_energy_run(spec: &Spec, st: &Stages) -> Result<String, String> {
    let c = spec.center;
    let mut spec = spec.clone();
    for p in &mut spec.positions {
        let r2 = (p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2);
        p[2] = 0.15 * (-r2 / 0.02).exp();
    }
    let sh = make(&spec)?;
    let phi0: Vec<f64> = sh.positions().iter().map(|p| p[0] - (c[0] - 0.3)).collect();
    let mut flow = Flow::new(&sh.mesh, &sh.geometry, &phi0).map_err(|e| e.to_string())?;
    flow.redistance().map_err(|e| e.to_string())?;
    let problem = CutProblem::new(Energy::Hencky, 0.004);
    // The cut starts far from the bump and the energy hardly changes for many steps there, so the
    // energy criterion would stop the run at once: only the motion criterion is used.
    let stop = StopRule {
        max_steps: 1000,
        motion_tolerance: 1e-5,
        energy_tolerance: 0.0,
        energy_window: 1,
    };
    let out = run_loop(
        &mut flow,
        &sh,
        &problem,
        st,
        StepSize::default(),
        &stop,
        Oracle::None,
    );
    let scale = flow.scale().length_scale;
    Ok(energy_line(
        &spec.name,
        &st.name,
        &summarize(&out, scale, mean_edge(&sh)),
    ))
}

/// Sphere: the equator with the perturbation of Fig. 1 (0.01 in the modes 2..8), Dirichlet.
fn sphere_energy_run(sh: &Shape, st: &Stages, alpha: f64) -> Result<String, String> {
    let mut rng = SplitMix64::new(5);
    let phases: Vec<f64> = (0..9).map(|_| rng.uniform(0.0, 2.0 * PI)).collect();
    let phi0: Vec<f64> = sh
        .positions()
        .iter()
        .map(|p| {
            let theta = p[1].atan2(p[0]);
            let bump: f64 = (2..=8)
                .map(|m| 0.01 * (m as f64 * theta + phases[m]).cos())
                .sum();
            polar(*p) - 0.5 * PI - bump
        })
        .collect();
    let mut flow = Flow::new(&sh.mesh, &sh.geometry, &phi0).map_err(|e| e.to_string())?;
    flow.redistance().map_err(|e| e.to_string())?;
    let problem = CutProblem::new(Energy::Dirichlet, alpha);
    let out = run_loop(
        &mut flow,
        sh,
        &problem,
        st,
        StepSize::default(),
        &StopRule::default(),
        Oracle::None,
    );
    let scale = flow.scale().length_scale;
    let line = energy_line(
        &sh.spec.name,
        &st.name.to_string(),
        &summarize(&out, scale, mean_edge(sh)),
    );
    Ok(format!("alpha {alpha}: {line}"))
}

fn energy_study() {
    let sets = [Stages::d(), Stages::p()];
    println!(
        "STAGE 4: energy history of whole runs, adaptive step. Flat: motion <= 1e-5 only, <= 1000 steps. Sphere: the default StopRule (<= 1000 steps, motion <= 1e-5 or energy fall < 1e-4 over 25 steps)."
    );
    println!(
        "Flat: bump 0.15 exp(-r^2/0.02), Hencky, alpha_L = 0.004, straight cut 0.3 left of the bump. Sphere: equator + Fig. 1 perturbation, Dirichlet.\n"
    );
    let flat = flat_family();
    let jobs: Vec<(usize, usize)> = (0..flat.len())
        .flat_map(|m| (0..sets.len()).map(move |s| (m, s)))
        .collect();
    let res = par_map(&jobs, |&(m, s)| flat_energy_run(&flat[m], &sets[s]));
    for (k, &(m, s)) in jobs.iter().enumerate() {
        match flatten(res[k].clone()) {
            Ok(l) => println!("{l}"),
            Err(e) => println!("{:<32} {:<3} {}", flat[m].name, sets[s].name, short(&e)),
        }
    }
    println!();
    let sph = build_all(&sphere_family());
    let alphas = [6.0, 3.0];
    let jobs: Vec<(usize, usize, usize)> = (0..sph.len())
        .flat_map(|m| (0..sets.len()).flat_map(move |s| (0..2).map(move |a| (m, s, a))))
        .collect();
    let res = par_map(&jobs, |&(m, s, a)| match &sph[m].1 {
        Ok(sh) => sphere_energy_run(sh, &sets[s], alphas[a]),
        Err(e) => Err(e.clone()),
    });
    for (k, &(m, s, a)) in jobs.iter().enumerate() {
        match flatten(res[k].clone()) {
            Ok(l) => println!("{l}"),
            Err(e) => println!(
                "{:<32} {:<3} alpha {}: {}",
                sph[m].0.name,
                sets[s].name,
                alphas[a],
                short(&e)
            ),
        }
    }
}

// ---------------------------------------------------------------------------------------------
// attribution: one stage option at a time

fn variants() -> Vec<Stages> {
    vec![
        Stages::d(),
        Stages::d().with("flux:Nodal", |s| s.flux = FluxRecovery::Nodal),
        Stages::d().with("flux:Nodal+sl", |s| {
            s.flux = FluxRecovery::Nodal;
            s.passes = 1;
        }),
        Stages::d().with("flux:c=0.25", |s| {
            s.flux = FluxRecovery::Smoothed { length: 0.25 }
        }),
        Stages::d().with("flux:c=1", |s| {
            s.flux = FluxRecovery::Smoothed { length: 1.0 }
        }),
        Stages::d().with("len:OnMesh", |s| s.length = LengthStep::OnMesh),
        Stages::d().with("ext:Harmonic", |s| s.ext = Extension::Harmonic),
        Stages::d().with("speed:Unit", |s| s.speed = SpeedScaling::Unit),
        Stages::d().with("redist:Off", |s| s.redist = Redist::Off),
        Stages::d().with("redist:FixedPt", |s| s.redist = Redist::FixedPoint),
        Stages::d().with("redist:Oracle", |s| s.redist = Redist::Oracle),
        Stages::ps(),
        Stages::p(),
    ]
}

fn attr_study() {
    let vars = variants();
    let steps = env_usize("LAW_STEPS", 96);
    println!(
        "ATTRIBUTION: the default stages with one option changed at a time (PS = all structural options of the paper, P = the paper as printed)."
    );
    println!(
        "redist:Off = no redistancing; FixedPt = cut-adjacent values kept, the rest by the fixed point of the upwind scheme; Oracle = cut-adjacent values kept, every other vertex gets the exact distance to the circle / latitude through the mean position of the cut (ideal far field).\n"
    );
    println!("CIRCLE LAW (relative error of the radius; roundness in brackets)");
    let flat = build_all(&flat_family());
    let jobs: Vec<(usize, usize)> = (0..flat.len())
        .flat_map(|m| (0..vars.len()).map(move |v| (m, v)))
        .collect();
    let res = par_map(&jobs, |&(m, v)| match &flat[m].1 {
        Ok(sh) => circle_law(sh, &vars[v], steps),
        Err(e) => Err(e.clone()),
    });
    print!("{:<30}", "mesh");
    for v in &vars {
        print!(" {:>16}", v.name);
    }
    println!();
    for m in 0..flat.len() {
        print!("{:<30}", flat[m].0.name);
        for v in 0..vars.len() {
            match flatten(res[m * vars.len() + v].clone()) {
                Ok(r) => print!(" {:>8}({:.0e})", sci(r.err), r.round),
                Err(e) => print!(" {:>16}", &short(&e)[..9.min(short(&e).len())]),
            }
        }
        println!();
    }
    println!();
    println!(
        "DISK LAW (flat, uniform source K = 5 in region 0, alpha_L = 0.3; relative error of the radius; roundness in brackets)"
    );
    let res = par_map(&jobs, |&(m, v)| match &flat[m].1 {
        Ok(sh) => disk_law(sh, &vars[v], steps),
        Err(e) => Err(e.clone()),
    });
    print!("{:<30}", "mesh");
    for v in &vars {
        print!(" {:>16}", v.name);
    }
    println!();
    for m in 0..flat.len() {
        print!("{:<30}", flat[m].0.name);
        for v in 0..vars.len() {
            match flatten(res[m * vars.len() + v].clone()) {
                Ok(r) => print!(" {:>8}({:.0e})", sci(r.err), r.round),
                Err(e) => print!(" {:>16}", &short(&e)[..9.min(short(&e).len())]),
            }
        }
        println!();
    }
    println!();
    let sph = build_all(&sphere_family());
    let tasks: Vec<(
        &str,
        Box<dyn Fn(&Shape, &Stages) -> Result<String, String> + Sync>,
    )> = vec![
        (
            "LATITUDE vs ODE alpha 20, psi0 0.5, t 4e-3, 200 steps: error of psi (spread)",
            Box::new(|sh, st| {
                lat_ode(sh, st, 20.0, 0.5, 4e-3, 200).map(|(e, s)| format!("{}({:.0e})", sci(e), s))
            }),
        ),
        (
            "BALANCED LATITUDE alpha 20 from 0.5: psi - psi* (steps)",
            Box::new(|sh, st| balance(sh, st, 0.5).map(|(e, _, n, _)| format!("{}({n})", sci(e)))),
        ),
        (
            "RIPPLE alpha 0, m 4 (distortion term only): measured/exact - 1",
            Box::new(|sh, st| {
                ripple(sh, st, 0.0, 4).map(|(m, e, _)| format!("{:+.3}", m / e - 1.0))
            }),
        ),
        (
            "RIPPLE alpha 10, m 4 (length term dominates): measured/exact - 1",
            Box::new(|sh, st| {
                ripple(sh, st, 10.0, 4).map(|(m, e, _)| format!("{:+.3}", m / e - 1.0))
            }),
        ),
    ];
    for (title, f) in &tasks {
        println!("{title}");
        let jobs: Vec<(usize, usize)> = (0..sph.len())
            .flat_map(|m| (0..vars.len()).map(move |v| (m, v)))
            .collect();
        let res = par_map(&jobs, |&(m, v)| match &sph[m].1 {
            Ok(sh) => f(sh, &vars[v]),
            Err(e) => Err(e.clone()),
        });
        print!("{:<30}", "mesh");
        for v in &vars {
            print!(" {:>16}", v.name);
        }
        println!();
        for m in 0..sph.len() {
            print!("{:<30}", sph[m].0.name);
            for v in 0..vars.len() {
                match flatten(res[m * vars.len() + v].clone()) {
                    Ok(s) => print!(" {:>16}", s),
                    Err(e) => print!(" {:>16}", &short(&e)[..9.min(short(&e).len())]),
                }
            }
            println!();
        }
        println!();
    }
}

// ---------------------------------------------------------------------------------------------
// real meshes

fn load_obj(path: &str) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>, usize), String> {
    let obj = ddg_core::io::read_obj(path).map_err(|e| format!("read_obj: {e}"))?;
    let mut faces = Vec::new();
    let mut polygons = 0;
    for f in 0..obj.n_faces() {
        let r = obj.face_corners(f);
        let c = &obj.corner_verts[r.clone()];
        if c.len() != 3 {
            polygons += 1;
        }
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    Ok((obj.positions, faces, polygons))
}

struct RealMesh {
    name: String,
    path: String,
    /// The plane: axis and offset.
    plane: Option<(usize, f64)>,
}

fn plane_of(positions: &[[f64; 3]]) -> (usize, f64) {
    let n = positions.len() as f64;
    let mut best = (0, 0.0, 0.0);
    for a in 0..3 {
        let lo = positions.iter().map(|p| p[a]).fold(f64::INFINITY, f64::min);
        let hi = positions
            .iter()
            .map(|p| p[a])
            .fold(f64::NEG_INFINITY, f64::max);
        if hi - lo > best.1 {
            best = (a, hi - lo, positions.iter().map(|p| p[a]).sum::<f64>() / n);
        }
    }
    (best.0, best.2)
}

fn real_run(
    path: &str,
    plane: Option<(usize, f64)>,
    vars: &[Stages],
    alpha_rules: &[(&str, Option<f64>)],
) -> Vec<String> {
    let mut lines = Vec::new();
    let (positions, faces, polygons) = match load_obj(path) {
        Ok(x) => x,
        Err(e) => return vec![format!("  load failed: {e}")],
    };
    let spec = Spec {
        name: path.into(),
        kind: Kind::Real,
        positions: positions.clone(),
        faces: faces.clone(),
        center: [0.0; 2],
    };
    let q = spec.quality();
    lines.push(format!(
        "  {}{}",
        quality_row("quality", &q),
        if polygons > 0 {
            format!("  ({polygons} polygons fan-triangulated)")
        } else {
            String::new()
        }
    ));
    lines.push(format!(
        "  median edge {:.4e}, longest/shortest edge {:.1e}",
        q.h_median, q.edge_ratio
    ));
    let sh = match spec.build() {
        Ok(s) => s,
        Err(e) => {
            lines.push(format!("  REJECTED by the crate: {e}"));
            return lines;
        }
    };
    let (axis, offset) = plane.unwrap_or_else(|| plane_of(&positions));
    let phi0: Vec<f64> = positions.iter().map(|p| p[axis] - offset).collect();
    let nbad = phi0.iter().filter(|x| **x == 0.0).count();
    lines.push(format!("  mesh accepted; components {}, boundary edges {}, start plane: axis {axis} offset {offset:.4} ({nbad} vertices exactly on it)", sh.mesh.n_components(), sh.mesh.n_boundary_edges()));
    for st in vars {
        for (rule_name, alpha) in alpha_rules {
            let mut flow = match Flow::new(&sh.mesh, &sh.geometry, &phi0) {
                Ok(f) => f,
                Err(e) => {
                    lines.push(format!("  Flow::new failed: {e}"));
                    return lines;
                }
            };
            if let Err(e) = flow.redistance() {
                lines.push(format!("  redistance failed: {e}"));
                continue;
            }
            let mut problem = CutProblem::new(Energy::Hencky, 0.0);
            let alpha_value = match alpha {
                Some(a) => *a,
                None => match flow.roughening_rate(&problem, &st.eval_options()) {
                    Ok(w) => length_weight_for_ripples(w, 0.25),
                    Err(e) => {
                        lines.push(format!(
                            "  {} {}: roughening_rate failed: {e}",
                            st.name, rule_name
                        ));
                        continue;
                    }
                },
            };
            problem.length_weight = alpha_value;
            let out = run_loop(
                &mut flow,
                &sh,
                &problem,
                st,
                StepSize::default(),
                &StopRule::default(),
                Oracle::None,
            );
            let scale = flow.scale().length_scale;
            let s = summarize(&out, scale, mean_edge(&sh));
            let first = out.steps.first();
            lines.push(format!(
                "  {:<14} alpha {:<6} ({:.3}): {}  [cut pts at start {}]",
                st.name,
                rule_name,
                alpha_value,
                energy_line("", "", &s).trim_start(),
                first.map_or(0, |r| r.n_cut_points)
            ));
        }
    }
    lines
}

fn real_study(attr: bool) {
    let root = std::env::var("REAL_ROOT")
        .unwrap_or_else(|_| "/var/data/discrete-geometry/ddg-tools/data".into());
    let cones = format!("{root}/optimal-cones/optimal-cones-data/data");
    let models = format!("{root}/models");
    let camel_plane = load_obj(&format!("{cones}/camel.obj"))
        .ok()
        .map(|(p, _, _)| plane_of(&p));
    let mut meshes = vec![
        RealMesh {
            name: "spot (triangulated)".into(),
            path: format!("{models}/spot/spot/spot_triangulated.obj"),
            plane: None,
        },
        RealMesh {
            name: "bob (tri)".into(),
            path: format!("{models}/bob/bob/bob_tri.obj"),
            plane: None,
        },
        RealMesh {
            name: "blub (triangulated)".into(),
            path: format!("{models}/blub/blub/blub_triangulated.obj"),
            plane: None,
        },
    ];
    for n in [
        "camel",
        "camelbad",
        "camelbadnonuniform",
        "camel_noise1",
        "camel_noise2",
        "camelshotnoise",
    ] {
        meshes.push(RealMesh {
            name: n.into(),
            path: format!("{cones}/{n}.obj"),
            plane: camel_plane,
        });
    }
    if let Ok(extra) = std::env::var("EXTRA_OBJ") {
        for p in extra.split(':').filter(|s| !s.is_empty()) {
            meshes.push(RealMesh {
                name: p.rsplit('/').next().unwrap().into(),
                path: p.into(),
                plane: None,
            });
        }
    }
    let vars: Vec<Stages> = if attr {
        // The fixed-point redistancing is left out on real meshes: with obtuse triangles it
        // re-opens vertices millions of times per step (see `redist`, "reopened").
        variants()
            .into_iter()
            .filter(|v| !matches!(v.redist, Redist::Oracle | Redist::FixedPoint))
            .collect()
    } else {
        vec![Stages::d()]
    };
    let rules: Vec<(&str, Option<f64>)> = if attr {
        vec![("0.5", Some(0.5))]
    } else {
        vec![("auto", None), ("0.5", Some(0.5))]
    };
    println!(
        "REAL MESHES: Hencky, planar start (plane through the centroid across the longest axis; camel variants use the plane of camel.obj), default StopRule, adaptive step."
    );
    println!(
        "alpha 'auto' = length_weight_for_ripples(omega, 0.25) from the roughening rate of the start; floor = median cut motion of the last quarter of the steps (<= 50) in median edges of the unit-area surface.\n"
    );
    println!("{}", quality_header());
    if attr {
        // One job per (mesh, variant).
        let jobs: Vec<(usize, usize)> = (0..meshes.len())
            .flat_map(|m| (0..vars.len()).map(move |v| (m, v)))
            .collect();
        let res = par_map(&jobs, |&(m, v)| {
            real_run(
                &meshes[m].path,
                meshes[m].plane,
                std::slice::from_ref(&vars[v]),
                &rules,
            )
        });
        for (m, mesh) in meshes.iter().enumerate() {
            println!("{} [{}]", mesh.name, mesh.path);
            for v in 0..vars.len() {
                match &res[m * vars.len() + v] {
                    Ok(lines) => {
                        for l in lines {
                            let keep = l.trim_start().starts_with(&vars[v].name)
                                || (v == 0
                                    && !l.contains("quality")
                                    && !l.contains("median edge")
                                    && !l.contains("mesh accepted"));
                            if keep {
                                println!("{l}");
                            }
                        }
                    }
                    Err(e) => println!("  {} PANIC: {e}", vars[v].name),
                }
            }
            println!();
        }
        return;
    }
    let res = par_map(&meshes, |m| real_run(&m.path, m.plane, &vars, &rules));
    for (m, r) in meshes.iter().zip(res) {
        println!("{} [{}]", m.name, m.path);
        match r {
            Ok(lines) => {
                for l in lines {
                    println!("{l}");
                }
            }
            Err(e) => println!("  PANIC: {e}"),
        }
        println!();
    }
}

fn quality_study() {
    println!(
        "MESH QUALITY (TREAT={}): largest angle per triangle: median / 95% / max (degrees), min angle, fraction of interior edges with a negative cotan weight. Angles of the intrinsic geometry. +V: vertices added; splits / flips of the export; seconds new/flip/refine/export.\n",
        Treat::from_env().tag()
    );
    println!(
        "{}  {:>6} {:>6} {:>6}  s new/flip/refine/export",
        quality_header(),
        "+V",
        "splits",
        "flips",
    );
    let specs: Vec<Spec> = flat_family().into_iter().chain(sphere_family()).collect();
    for (s, r) in build_all(&specs) {
        match r {
            Ok(sh) => println!(
                "{}  {:>6} {:>6} {:>6}  {:.3}/{:.3}/{:.3}/{:.3}",
                quality_row(&s.name, &sh.q),
                sh.info.added,
                sh.info.n_splits,
                sh.info.n_export_flips,
                sh.info.t_new,
                sh.info.t_flip,
                sh.info.t_refine,
                sh.info.t_export
            ),
            Err(e) => println!("{:<34} {}", s.name, short(&e)),
        }
    }
}

// ---------------------------------------------------------------------------------------------
// TREATMENT: real meshes on the intrinsic triangulation, the way back, speed

fn real_list() -> Vec<RealMesh> {
    let root = std::env::var("REAL_ROOT")
        .unwrap_or_else(|_| "/var/data/discrete-geometry/ddg-tools/data".into());
    let cones = format!("{root}/optimal-cones/optimal-cones-data/data");
    let models = format!("{root}/models");
    let camel_plane = load_obj(&format!("{cones}/camel.obj"))
        .ok()
        .map(|(p, _, _)| plane_of(&p));
    let mut meshes = vec![
        RealMesh {
            name: "spot".into(),
            path: format!("{models}/spot/spot/spot_triangulated.obj"),
            plane: None,
        },
        RealMesh {
            name: "bob".into(),
            path: format!("{models}/bob/bob/bob_tri.obj"),
            plane: None,
        },
        RealMesh {
            name: "blub".into(),
            path: format!("{models}/blub/blub/blub_triangulated.obj"),
            plane: None,
        },
    ];
    for n in [
        "camel",
        "camelbad",
        "camelbadnonuniform",
        "camel_noise1",
        "camel_noise2",
        "camelshotnoise",
    ] {
        meshes.push(RealMesh {
            name: n.into(),
            path: format!("{cones}/{n}.obj"),
            plane: camel_plane,
        });
    }
    meshes
}

/// Connected pieces of a cut from its segments (union-find over the node indices).
fn n_pieces(n_nodes: usize, segments: &[Segment]) -> usize {
    let mut parent: Vec<usize> = (0..n_nodes).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut used = vec![false; n_nodes];
    for s in segments {
        let (a, b) = (s.nodes[0] as usize, s.nodes[1] as usize);
        used[a] = true;
        used[b] = true;
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent[ra] = rb;
    }
    (0..n_nodes)
        .filter(|&i| used[i] && find(&mut parent, i) == i)
        .count()
}

/// The cut of `phi` on the export, carried back to the input mesh: (nodes, segments, pieces of
/// the input path, connected pieces, largest relative mismatch, largest absolute mismatch,
/// seconds).
fn way_back(
    tri: &IntrinsicTriangulation,
    export: &SimplicialExport,
    phi: &[f64],
) -> Result<(usize, usize, usize, usize, f64, f64, f64), String> {
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
    let start = Instant::now();
    let paths = tri
        .segments_on_input(&nodes, &segments, &EdgeTraceOptions::default())
        .map_err(|e| e.to_string())?;
    let seconds = start.elapsed().as_secs_f64();
    let pieces: usize = paths.paths.iter().map(|p| p.faces.len()).sum();
    let abs = paths.paths.iter().map(|p| p.mismatch).fold(0.0, f64::max);
    Ok((
        nodes.len(),
        segments.len(),
        pieces,
        n_pieces(nodes.len(), &segments),
        paths.max_relative_mismatch,
        abs,
        seconds,
    ))
}

fn real_job(m: &RealMesh, treat: Treat) -> Vec<String> {
    let mut lines = Vec::new();
    let (positions, faces, _) = match load_obj(&m.path) {
        Ok(x) => x,
        Err(e) => return vec![format!("  load failed: {e}")],
    };
    let spec = Spec {
        name: m.name.clone(),
        kind: Kind::Real,
        positions: positions.clone(),
        faces,
        center: [0.0; 2],
    };
    let input = match spec.build() {
        Ok(s) => s,
        Err(e) => return vec![format!("  input mesh REJECTED by the crate: {e}")],
    };
    let tr = match treat_full(&spec, treat) {
        Ok(t) => t,
        Err(e) => return vec![format!("  treatment FAILED: {e}")],
    };
    let info = tr.shape.info.clone();
    lines.push(format!(
        "  {} V {} (+{}), largest angle med/p95/max {:.1}/{:.1}/{:.1}, min {:.2}, neg-cot {:.2}%, splits {}, export flips {}; seconds new {:.3} flip {:.3} refine {:.3} export {:.3} (sum {:.3})",
        treat.tag(), tr.shape.q.nv, info.added, tr.shape.q.largest_median, tr.shape.q.largest_p95, tr.shape.q.largest_max,
        tr.shape.q.smallest, 100.0 * tr.shape.q.neg_cotan, info.n_splits, info.n_export_flips,
        info.t_new, info.t_flip, info.t_refine, info.t_export, info.t_new + info.t_flip + info.t_refine + info.t_export
    ));
    let sh = &tr.shape;
    let (axis, offset) = m.plane.unwrap_or_else(|| plane_of(&positions));
    let phi0: Vec<f64> = sh.positions().iter().map(|p| p[axis] - offset).collect();
    let mut flow = match Flow::new(&sh.mesh, &sh.geometry, &phi0) {
        Ok(f) => f,
        Err(e) => return vec![format!("  Flow::new failed: {e}")],
    };
    if let Err(e) = flow.redistance() {
        lines.push(format!("  redistance failed: {e}"));
        return lines;
    }
    let problem = CutProblem::new(Energy::Hencky, 0.5);
    let st = Stages::d();
    let out = run_loop(
        &mut flow,
        sh,
        &problem,
        &st,
        StepSize::default(),
        &StopRule::default(),
        Oracle::None,
    );
    let scale = flow.scale().length_scale;
    let s = summarize(&out, scale, mean_edge(sh));
    lines.push(format!(
        "  {}  returned_step {:?}, total {:.2} s",
        energy_line("", "", &s).trim_start(),
        out.returned_step,
        out.seconds
    ));
    // The final cut, evaluated by the same code on the input mesh: the values at the input
    // vertices (the first vertices of the export) as a level-set function of the input mesh.
    let phi = flow.level_set();
    let own = flow
        .evaluate(&problem, &st.eval_options())
        .map(|e| e.energy.total)
        .unwrap_or(f64::NAN);
    let on_input = (|| -> Result<f64, String> {
        let mut ev = ddg_variational_cuts::Evaluator::new(&input.mesh, &input.geometry)
            .map_err(|e| e.to_string())?;
        let e = ev
            .evaluate(&phi[..positions.len()], &problem, &st.eval_options())
            .map_err(|e| e.to_string())?;
        Ok(e.energy.total)
    })();
    lines.push(format!(
        "  final energy: own (on the mesh the flow ran on) {:.5e}; the same cut restricted to the input vertices, evaluated on the input mesh {}",
        own,
        on_input.map_or_else(|e| short(&e), |x| format!("{x:.5e}"))
    ));
    if let (Some(tri), Some(export)) = (&tr.tri, &tr.export) {
        match way_back(tri, export, &phi) {
            Ok((nodes, segs, pieces, comps, rel, abs, secs)) => lines.push(format!(
                "  way back: {nodes} nodes, {segs} segments, {pieces} pieces on input faces, {comps} connected piece(s), largest mismatch {rel:.1e} of a segment ({abs:.1e} absolute), {secs:.4} s"
            )),
            Err(e) => lines.push(format!("  way back FAILED: {e}")),
        }
    }
    lines
}

fn real_intr_study() {
    let treat = Treat::from_env();
    println!(
        "REAL MESHES, TREAT={}: Hencky, alpha_L = 0.5, planar start (plane of the longest axis through the centroid; camels: the plane of camel.obj), default StopRule (new rule), adaptive step, default stages.\n",
        treat.tag()
    );
    let meshes = real_list();
    let res = par_map(&meshes, |m| real_job(m, treat));
    for (m, r) in meshes.iter().zip(res) {
        println!("{}", m.name);
        match r {
            Ok(lines) => lines.iter().for_each(|l| println!("{l}")),
            Err(e) => println!("  PANIC {e}"),
        }
        println!();
    }
}

/// Timings, run alone (one process, one thread): treatment stages, flow per step, way back.
fn speed_study() {
    println!(
        "SPEED (single thread, run alone): seconds. new/flip/refine/export: the triangulation; step: median of steps 5..25 of the flow on the mesh it runs on; back: the cut of the last step carried to the input mesh.\n"
    );
    println!(
        "{:<28} {:<5} {:>7} {:>7} {:>7} {:>7} {:>8} {:>9} {:>9} {:>8}",
        "mesh", "treat", "V", "new", "flip", "refine", "export", "step", "10 steps", "back"
    );
    let mut specs: Vec<(String, Spec)> = Vec::new();
    for m in real_list() {
        if let Ok((positions, faces, _)) = load_obj(&m.path) {
            specs.push((
                m.name.clone(),
                Spec {
                    name: m.name,
                    kind: Kind::Real,
                    positions,
                    faces,
                    center: [0.0; 2],
                },
            ));
        }
    }
    for s in flat_family()
        .into_iter()
        .filter(|s| s.name.contains("165") || s.name.contains("1:16"))
    {
        specs.push((s.name.clone(), s));
    }
    for (name, spec) in &specs {
        let (axis, offset) = if spec.kind == Kind::Real {
            plane_of(&spec.positions)
        } else {
            (0, 0.5)
        };
        for treat in [Treat::Input, Treat::Flip, Treat::R25, Treat::R30] {
            let tr = match treat_full(spec, treat) {
                Ok(t) => t,
                Err(e) => {
                    println!("{:<28} {:<5} {}", name, treat.tag(), short(&e));
                    continue;
                }
            };
            let sh = &tr.shape;
            let phi0: Vec<f64> = sh.positions().iter().map(|p| p[axis] - offset).collect();
            let mut flow = Flow::new(&sh.mesh, &sh.geometry, &phi0).unwrap();
            let _ = flow.redistance();
            let problem = CutProblem::new(Energy::Hencky, 0.5);
            let st = Stages::d();
            let mut times = Vec::new();
            for _ in 0..25 {
                let start = Instant::now();
                match do_step(
                    &mut flow,
                    sh,
                    &problem,
                    &st,
                    StepSize::default(),
                    Oracle::None,
                ) {
                    Ok(_) => times.push(start.elapsed().as_secs_f64()),
                    Err(_) => break,
                }
            }
            let mut mid: Vec<f64> = times.iter().skip(5).cloned().collect();
            mid.sort_by(f64::total_cmp);
            let step = mid.get(mid.len() / 2).copied().unwrap_or(f64::NAN);
            let ten: f64 = times.iter().skip(5).take(10).sum();
            let back = match (&tr.tri, &tr.export) {
                (Some(tri), Some(export)) => {
                    way_back(tri, export, &flow.level_set()).map_or(f64::NAN, |r| r.6)
                }
                _ => 0.0,
            };
            let i = &sh.info;
            println!(
                "{:<28} {:<5} {:>7} {:>7.3} {:>7.3} {:>7.3} {:>8.3} {:>9.4} {:>9.3} {:>8.4}",
                name,
                treat.tag(),
                sh.mesh.n_vertices(),
                i.t_new,
                i.t_flip,
                i.t_refine,
                i.t_export,
                step,
                ten,
                back
            );
        }
    }
}

fn main() {
    paper_check();
    let what = std::env::args().nth(1).unwrap_or_default();
    match what.as_str() {
        "quality" => quality_study(),
        "eval" => eval_study(),
        "redist" => redist_study(),
        "laws" => laws_study(),
        "energy" => energy_study(),
        "attr" => attr_study(),
        "real" => real_intr_study(),
        "speed" => speed_study(),
        _ => {
            quality_study();
            eval_study();
            redist_study();
            laws_study();
            energy_study();
            attr_study();
            real_intr_study();
        }
    }
}
