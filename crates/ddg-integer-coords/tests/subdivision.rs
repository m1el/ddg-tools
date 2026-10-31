//! The common subdivision, the interpolation matrices and the point queries (WP5; SPEC §4.8,
//! §4.9, §8 layer 2).
//!
//! The states are those of `tests/states.rs`: signpost triangulations of the standard shapes, spot
//! and camelbad in 7 scenarios and 2 seeds, converted to integer coordinates. The integer
//! subdivision is compared with the signpost crate's overlay of the triangulation it was converted
//! from. Run with `--nocapture` to see the tables.

#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord, clippy::manual_contains)]

mod common;

use std::sync::OnceLock;

use common::{SCENARIOS, SEEDS, Scenario, Shape, all_shapes, convert, state};
use ddg_core::generators::SplitMix64;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::{IntegerTriangulation, Parts};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, IntrinsicTriangulation, Overlay, SurfacePoint,
};

struct Job {
    name: String,
    shape: usize,
    signposts: IntrinsicTriangulation,
    overlay: Result<Overlay, String>,
    tri: IntegerTriangulation,
}

fn shapes() -> &'static [Shape] {
    static SHAPES: OnceLock<Vec<Shape>> = OnceLock::new();
    SHAPES.get_or_init(all_shapes)
}

fn jobs() -> &'static [Job] {
    static JOBS: OnceLock<Vec<Job>> = OnceLock::new();
    JOBS.get_or_init(|| {
        let shapes = shapes();
        let mut specs = Vec::new();
        for s in 0..shapes.len() {
            for scenario in SCENARIOS {
                for seed in SEEDS {
                    specs.push((s, scenario, seed));
                }
            }
        }
        par_map(&specs, |&(s, scenario, seed): &(usize, Scenario, u64)| {
            let signposts = state(&shapes[s], scenario, seed);
            let (tri, _) = convert(&signposts).ok()?;
            let overlay = signposts
                .trace_edges(&EdgeTraceOptions::default())
                .map_err(|e| e.to_string())
                .and_then(|t| signposts.overlay(&t).map_err(|e| e.to_string()));
            Some(Job {
                name: format!("{} / {scenario:?} / {seed}", shapes[s].name),
                shape: s,
                signposts,
                overlay,
                tri,
            })
        })
        .into_iter()
        .flatten()
        .collect()
    })
}

fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(48);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let out: Vec<std::sync::Mutex<Option<R>>> =
        items.iter().map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= items.len() {
                        break;
                    }
                    *out[i].lock().unwrap() = Some(f(&items[i]));
                }
            });
        }
    });
    out.into_iter()
        .map(|m| m.into_inner().unwrap().unwrap())
        .collect()
}

#[derive(Debug, Clone, Default)]
struct Worst {
    value: f64,
    place: String,
}

impl Worst {
    fn note(&mut self, value: f64, place: impl FnOnce() -> String) {
        if value > self.value || value.is_nan() {
            self.value = value;
            self.place = place();
        }
    }
    fn merge(&mut self, other: &Worst) {
        if other.value > self.value || other.value.is_nan() {
            *self = other.clone();
        }
    }
}

fn position(shape: &Shape, point: SurfacePoint) -> [f64; 3] {
    let p = &shape.positions;
    match point {
        SurfacePoint::Vertex(v) => p[v as usize],
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = shape.mesh.edge_vertices(edge as usize);
            [0, 1, 2].map(|k| (1.0 - t) * p[a][k] + t * p[b][k])
        }
        SurfacePoint::Face { face, bary } => {
            let f = shape.mesh.faces()[face as usize];
            [0, 1, 2].map(|k| (0..3).map(|c| bary[c] * p[f[c] as usize][k]).sum())
        }
    }
}

fn minus(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
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
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    dot(minus(a, b), minus(a, b)).sqrt()
}

fn in_input_face(shape: &Shape, point: SurfacePoint, f: usize) -> bool {
    match point {
        SurfacePoint::Vertex(v) => shape.mesh.faces()[f].contains(&v),
        SurfacePoint::Edge { edge, .. } => {
            (0..3).any(|k| shape.mesh.edge(3 * f + k) == edge as usize)
        }
        SurfacePoint::Face { face, .. } => face as usize == f,
    }
}

fn in_face(tri: &IntrinsicTriangulation, point: SurfacePoint, f: usize) -> bool {
    match point {
        SurfacePoint::Vertex(v) => tri.faces()[f].contains(&v),
        SurfacePoint::Edge { edge, .. } => (0..3).any(|k| tri.edge(3 * f + k) == edge as usize),
        SurfacePoint::Face { face, .. } => face as usize == f,
    }
}

/// The polygon rotated to its lexicographically smallest rotation.
fn canonical(p: &[u32]) -> Vec<u32> {
    (0..p.len())
        .map(|r| {
            (0..p.len())
                .map(|i| p[(r + i) % p.len()])
                .collect::<Vec<_>>()
        })
        .min()
        .unwrap()
}

fn kind(p: SurfacePoint) -> (u8, u32) {
    match p {
        SurfacePoint::Vertex(v) => (0, v),
        SurfacePoint::Edge { edge, .. } => (1, edge),
        SurfacePoint::Face { face, .. } => (2, face),
    }
}

#[derive(Default)]
struct Stats {
    worst_raw: Worst,
    worst_conditioned: Worst,
    worst_u3d: Worst,
    worst_area_input: Worst,
    worst_area_tri: Worst,
    worst_area_total: Worst,
    crossings: usize,
    faces: usize,
    refused: usize,
    problems: Vec<String>,
}

fn check_job(job: &Job) -> Stats {
    let mut st = Stats::default();
    let shape = &shapes()[job.shape];
    let mesh = &shape.mesh;
    let t = &job.tri;
    let c = t.complex();
    let subd = match t.common_subdivision() {
        Ok(s) => s,
        Err(e) => {
            st.problems
                .push(format!("{}: common_subdivision: {e}", job.name));
            return st;
        }
    };
    let overlay = &subd.overlay;
    let nv = c.n_vertices();
    // The signposts' overlay, if it exists, is the reference.
    match &job.overlay {
        Err(_) => st.refused += 1,
        Ok(reference) => {
            if reference.vertices.len() != overlay.vertices.len()
                || reference.n_faces() != overlay.n_faces()
            {
                st.problems.push(format!(
                    "{}: sizes: {} vertices and {} faces, signposts {} and {}",
                    job.name,
                    overlay.vertices.len(),
                    overlay.n_faces(),
                    reference.vertices.len(),
                    reference.n_faces()
                ));
                return st;
            }
            // Sines of the crossing angles, by overlay vertex.
            let mut offsets = vec![0usize; c.n_edges() + 1];
            for e in 0..c.n_edges() {
                offsets[e + 1] = offsets[e] + t.normal_coordinate(e).max(0) as usize;
            }
            let mut sin = vec![f64::NAN; offsets[c.n_edges()]];
            for ie in 0..mesh.n_edges() {
                if let Ok(g) = t.extract_edge_geometry(mesh.edge_halfedge(ie)) {
                    for x in g.crossings.iter().flatten() {
                        let h = x.crossing.halfedge as usize;
                        let e = c.edge(h);
                        let n = t.normal_coordinate(e);
                        let i = if c.is_forward(h) {
                            x.crossing.index
                        } else {
                            n - 1 - x.crossing.index
                        };
                        sin[offsets[e] + i as usize] = x.sin_angle;
                    }
                }
            }
            for (i, (mine, theirs)) in overlay.vertices.iter().zip(&reference.vertices).enumerate()
            {
                if i < nv {
                    if mine != theirs {
                        st.problems
                            .push(format!("{}: vertex {i}: {mine:?} vs {theirs:?}", job.name));
                    }
                    continue;
                }
                st.crossings += 1;
                if kind(mine.on_input) != kind(theirs.on_input)
                    || kind(mine.on_triangulation) != kind(theirs.on_triangulation)
                {
                    st.problems.push(format!(
                        "{}: crossing vertex {i}: {mine:?} vs {theirs:?}",
                        job.name
                    ));
                    continue;
                }
                let (SurfacePoint::Edge { t: u1, .. }, SurfacePoint::Edge { t: u2, .. }) =
                    (mine.on_input, theirs.on_input)
                else {
                    continue;
                };
                let (SurfacePoint::Edge { t: v1, .. }, SurfacePoint::Edge { t: v2, .. }) =
                    (mine.on_triangulation, theirs.on_triangulation)
                else {
                    continue;
                };
                let d = (u1 - u2).abs().max((v1 - v2).abs());
                let s = sin[i - nv];
                st.worst_raw.note(d, || format!("{} vertex {i}", job.name));
                st.worst_conditioned
                    .note(d * s, || format!("{} vertex {i}", job.name));
                st.worst_u3d.note(
                    dist(
                        position(shape, mine.on_input),
                        position(shape, theirs.on_input),
                    ) * s,
                    || format!("{} vertex {i}", job.name),
                );
                if !(d * s <= 1e-9) {
                    st.problems.push(format!(
                        "{}: crossing vertex {i}: |Δ| {d:.2e}, sin θ {s:.2e}",
                        job.name
                    ));
                }
            }
            // Faces as labelled polygons, per face of the triangulation.
            let polygons = |o: &Overlay| {
                let mut per_face: Vec<Vec<(u32, Vec<u32>)>> = vec![Vec::new(); c.n_faces()];
                for f in 0..o.n_faces() {
                    per_face[o.face_triangulation[f] as usize]
                        .push((o.face_input[f], canonical(o.face(f))));
                }
                for v in &mut per_face {
                    v.sort();
                }
                per_face
            };
            let (a, b) = (polygons(overlay), polygons(reference));
            for (f, (x, y)) in a.iter().zip(&b).enumerate() {
                if x != y {
                    st.problems.push(format!(
                        "{}: face {f} of the triangulation: {x:?} vs signposts {y:?}",
                        job.name
                    ));
                    break;
                }
            }
        }
    }
    // The checks of the signpost crate's overlay (IT-SPEC N12) on ours.
    let positions: Vec<[f64; 3]> = overlay
        .vertices
        .iter()
        .map(|v| position(shape, v.on_input))
        .collect();
    let mut area_in_input = vec![0.0; mesh.n_faces()];
    let mut area_in_tri = vec![0.0; c.n_faces()];
    let mut pieces_of = vec![0usize; c.n_faces()];
    for f in 0..overlay.n_faces() {
        let corners = overlay.face(f);
        let (g, tf) = (
            overlay.face_input[f] as usize,
            overlay.face_triangulation[f] as usize,
        );
        if !(3..=6).contains(&corners.len())
            || g >= mesh.n_faces()
            || tf >= c.n_faces()
            || (f > 0 && overlay.face_triangulation[f - 1] > overlay.face_triangulation[f])
        {
            st.problems.push(format!(
                "{}: piece {f}: {corners:?}, faces {g} / {tf}",
                job.name
            ));
            continue;
        }
        pieces_of[tf] += 1;
        for &k in corners {
            let v = overlay.vertices[k as usize];
            if !in_input_face(shape, v.on_input, g)
                || !in_face(&job.signposts, v.on_triangulation, tf)
            {
                st.problems.push(format!(
                    "{}: piece {f}: corner {k} {v:?} is not in faces {g} / {tf}",
                    job.name
                ));
            }
        }
        let normal = {
            let [a, b, cc] = mesh.faces()[g].map(|v| shape.positions[v as usize]);
            cross(minus(b, a), minus(cc, a))
        };
        let scale = dot(normal, normal).sqrt();
        let n = corners.len();
        let mut area = 0.0;
        for k in 0..n {
            let (a, b, cc) = (
                positions[corners[k] as usize],
                positions[corners[(k + 1) % n] as usize],
                positions[corners[(k + 2) % n] as usize],
            );
            let turn = dot(cross(minus(b, a), minus(cc, b)), normal) / scale;
            if turn <= 0.0 {
                st.problems.push(format!(
                    "{}: piece {f} is not convex at corner {k}: {turn}",
                    job.name
                ));
            }
            if k + 2 < n {
                let o = positions[corners[0] as usize];
                area += 0.5 * dot(cross(minus(b, o), minus(cc, o)), normal) / scale;
            }
        }
        if area <= 0.0 {
            st.problems
                .push(format!("{}: piece {f} has area {area}", job.name));
        }
        area_in_input[g] += area;
        area_in_tri[tf] += area;
    }
    if pieces_of.iter().any(|&n| n == 0) {
        st.problems.push(format!(
            "{}: a face of the triangulation has no piece",
            job.name
        ));
    }
    let total = shape.geometry.total_area();
    for (g, &a) in area_in_input.iter().enumerate() {
        let want = shape.geometry.face_areas()[g];
        st.worst_area_input.note((a - want).abs() / want, || {
            format!("{} input face {g}", job.name)
        });
        if (a - want).abs() >= 1e-9 * total {
            st.problems.push(format!(
                "{}: pieces in input face {g}: {a}, face {want}",
                job.name
            ));
        }
    }
    for (f, &a) in area_in_tri.iter().enumerate() {
        let want = job.signposts.face_area(f);
        st.worst_area_tri
            .note((a - want).abs() / want, || format!("{} face {f}", job.name));
        if (a - want).abs() >= 1e-9 * total {
            st.problems.push(format!(
                "{}: pieces in face {f}: {a}, face {want}",
                job.name
            ));
        }
    }
    let sum: f64 = area_in_input.iter().sum();
    st.worst_area_total
        .note((sum - total).abs() / total, || job.name.clone());
    let (triangles, _) = overlay.triangles();
    match TriMesh::new(overlay.vertices.len(), triangles) {
        Err(e) => st
            .problems
            .push(format!("{}: not a valid mesh: {e}", job.name)),
        Ok(m) => {
            if m.euler_characteristic() != mesh.euler_characteristic()
                || m.boundary_loops().len() != mesh.boundary_loops().len()
                || m.n_components() != mesh.n_components()
            {
                st.problems
                    .push(format!("{}: topology differs from the input", job.name));
            }
            if let Err(e) = Geometry::from_positions(&m, &positions) {
                st.problems
                    .push(format!("{}: degenerate triangle: {e}", job.name));
            }
        }
    }
    if subd.fingerprint != t.fingerprint() || subd.report.n_vertices != overlay.vertices.len() {
        st.problems
            .push(format!("{}: fingerprint or report", job.name));
    }
    st.faces += overlay.n_faces();
    st
}

#[test]
fn subdivision_matches_the_signpost_overlay() {
    let jobs = jobs();
    assert!(jobs.len() >= 200, "{} states", jobs.len());
    let results = par_map(jobs, check_job);
    let mut all = Stats::default();
    for r in &results {
        all.worst_raw.merge(&r.worst_raw);
        all.worst_conditioned.merge(&r.worst_conditioned);
        all.worst_u3d.merge(&r.worst_u3d);
        all.worst_area_input.merge(&r.worst_area_input);
        all.worst_area_tri.merge(&r.worst_area_tri);
        all.worst_area_total.merge(&r.worst_area_total);
        all.crossings += r.crossings;
        all.faces += r.faces;
        all.refused += r.refused;
        all.problems.extend(r.problems.iter().cloned());
    }
    eprintln!(
        "subdivision vs signposts: {} states ({} without a signpost overlay), {} pieces, {} crossings\n  \
         raw max(|Δu|, |Δv|) {:.2e} at {}\n  conditioned (· sin θ) {:.2e} at {}; 3-D distance of on_input · sin θ {:.2e}\n  \
         worst relative area error: input faces {:.2e} ({}), faces of T¹ {:.2e} ({}), total {:.2e}",
        jobs.len(),
        all.refused,
        all.faces,
        all.crossings,
        all.worst_raw.value,
        all.worst_raw.place,
        all.worst_conditioned.value,
        all.worst_conditioned.place,
        all.worst_u3d.value,
        all.worst_area_input.value,
        all.worst_area_input.place,
        all.worst_area_tri.value,
        all.worst_area_tri.place,
        all.worst_area_total.value
    );
    for p in all.problems.iter().take(30) {
        eprintln!("  PROBLEM {p}");
    }
    assert!(all.problems.is_empty(), "{} problems", all.problems.len());
}

/// Row sums and entry counts of `P₀` and `P₁`, `P₀ x⁰` against the vertices' positions, and on
/// flat shapes `P₁ x¹ = P₀ x⁰` for `x¹ = P₀ x⁰` restricted to the vertices of `T¹`.
#[test]
fn interpolation_matrices() {
    let jobs = jobs();
    let results = par_map(jobs, |job| {
        let shape = &shapes()[job.shape];
        let t = &job.tri;
        let mut problems = Vec::new();
        let (mut worst_sum, mut worst_pos, mut worst_lin) = (0.0f64, 0.0f64, 0.0f64);
        let sub = t.common_subdivision().unwrap();
        let (p0, p1) = t.interpolation_matrices(&sub).unwrap();
        let ns = sub.overlay.vertices.len();
        if (p0.nrows(), p0.ncols(), p1.nrows(), p1.ncols())
            != (ns, t.n_input_vertices(), ns, t.complex().n_vertices())
        {
            problems.push(format!("{}: shapes", job.name));
            return (problems, 0.0, 0.0, 0.0, 0, 0);
        }
        let ones0 = vec![1.0; p0.ncols()];
        let ones1 = vec![1.0; p1.ncols()];
        let (mut s0, mut s1) = (vec![0.0; ns], vec![0.0; ns]);
        p0.mul_vec(&ones0, &mut s0).unwrap();
        p1.mul_vec(&ones1, &mut s1).unwrap();
        for i in 0..ns {
            worst_sum = worst_sum.max((s0[i] - 1.0).abs()).max((s1[i] - 1.0).abs());
        }
        let mut count0 = vec![0usize; ns];
        let mut count1 = vec![0usize; ns];
        for (r, _, _) in p0.triplets() {
            count0[r as usize] += 1;
        }
        for (r, _, _) in p1.triplets() {
            count1[r as usize] += 1;
        }
        for i in 0..ns {
            let want0 = match sub.overlay.vertices[i].on_input {
                SurfacePoint::Vertex(_) => 1,
                SurfacePoint::Edge { .. } => 2,
                SurfacePoint::Face { .. } => 3,
            };
            let ok1 = if i < t.complex().n_vertices() {
                count1[i] == 1
            } else {
                (1..=2).contains(&count1[i])
            }; // a self-edge gives one entry
            if count0[i] != want0 || !ok1 {
                problems.push(format!(
                    "{}: row {i}: {} / {} entries, expected {want0}",
                    job.name, count0[i], count1[i]
                ));
                break;
            }
        }
        // Positions.
        let n0 = t.n_input_vertices();
        let size = shape.mean_edge_length();
        let mut at0: Vec<Vec<f64>> = vec![vec![0.0; ns]; 3];
        for k in 0..3 {
            let x: Vec<f64> = (0..n0).map(|v| shape.positions[v][k]).collect();
            p0.mul_vec(&x, &mut at0[k]).unwrap();
        }
        for i in 0..ns {
            let want = position(shape, sub.overlay.vertices[i].on_input);
            for k in 0..3 {
                worst_pos = worst_pos.max((at0[k][i] - want[k]).abs() / size);
            }
        }
        if shape.flat {
            for k in 0..3 {
                let x1: Vec<f64> = at0[k][..t.complex().n_vertices()].to_vec();
                let mut y = vec![0.0; ns];
                p1.mul_vec(&x1, &mut y).unwrap();
                for i in 0..ns {
                    worst_lin = worst_lin.max((y[i] - at0[k][i]).abs() / size);
                }
            }
        }
        (
            problems,
            worst_sum,
            worst_pos,
            worst_lin,
            ns,
            usize::from(shape.flat),
        )
    });
    let (mut sum, mut pos, mut lin) = (0.0f64, 0.0f64, 0.0f64);
    let (mut rows, mut flat_states) = (0, 0);
    let mut problems = Vec::new();
    for (p, a, b, c, n, f) in results {
        problems.extend(p);
        sum = sum.max(a);
        pos = pos.max(b);
        lin = lin.max(c);
        rows += n;
        flat_states += f;
    }
    eprintln!(
        "P0, P1: {rows} rows in {} states; worst row sum defect {sum:.2e}; P0 x0 vs positions {pos:.2e} (edge lengths); \
         flat shapes ({flat_states} states): P1 (P0 x0 on V1) vs P0 x0 {lin:.2e}",
        jobs.len()
    );
    for p in problems.iter().take(20) {
        eprintln!("  PROBLEM {p}");
    }
    assert!(problems.is_empty());
    assert!(
        sum < 1e-14 && pos < 1e-12 && lin < 1e-8,
        "sum {sum}, pos {pos}, lin {lin}"
    );
}

/// A stale subdivision is refused by the matrices (SPEC I19), and the fingerprint changes with
/// every flip.
#[test]
fn a_stale_subdivision_is_detected() {
    let job = &jobs()[jobs().len() / 2];
    let mut t = job.tri.clone();
    let sub = t.common_subdivision().unwrap();
    assert_eq!(sub.fingerprint, t.fingerprint());
    let mut flipped = false;
    for e in 0..t.complex().n_edges() {
        if t.flip_edge(e).is_ok() {
            flipped = true;
            break;
        }
    }
    assert!(flipped);
    assert_ne!(sub.fingerprint, t.fingerprint());
    assert!(t.interpolation_matrices(&sub).is_err());
    let fresh = t.common_subdivision().unwrap();
    assert!(t.interpolation_matrices(&fresh).is_ok());
}

/// Point queries `T⁰ → T¹` against the signposts' `point_on_triangulation`.
#[test]
fn point_queries_match_the_signposts() {
    let jobs = jobs();
    let results = par_map(jobs, |job| {
        let shape = &shapes()[job.shape];
        let t = &job.tri;
        let mut rng = SplitMix64::new(77);
        let mut locator = t.point_locator();
        let size = shape.mean_edge_length();
        let (mut worst_back, mut worst_theirs, mut worst_bary) = (0.0f64, 0.0f64, 0.0f64);
        let (mut n, mut skipped) = (0usize, 0usize);
        let mut problems = Vec::new();
        for q in 0..40 {
            let point = if q % 4 == 3 {
                SurfacePoint::Edge {
                    edge: (rng.next_u64() % shape.mesh.n_edges() as u64) as u32,
                    t: rng.uniform(0.01, 0.99),
                }
            } else {
                let f = (rng.next_u64() % shape.mesh.n_faces() as u64) as u32;
                let (mut a, mut b) = (rng.next_f64(), rng.next_f64());
                if a + b > 1.0 {
                    (a, b) = (1.0 - a, 1.0 - b);
                }
                SurfacePoint::Face {
                    face: f,
                    bary: [a, b, 1.0 - a - b],
                }
            };
            let mine = match locator.point_on_triangulation(point) {
                Ok(p) => p,
                Err(e) => {
                    problems.push(format!("{}: {point}: {e}", job.name));
                    continue;
                }
            };
            n += 1;
            match job.signposts.point_on_input(mine) {
                Ok(back) => {
                    let d = dist(position(shape, back), position(shape, point)) / size;
                    worst_back = worst_back.max(d);
                    if d > 1e-8 {
                        problems.push(format!(
                            "{}: {point} -> {mine}, back at distance {d:.2e}",
                            job.name
                        ));
                    }
                }
                Err(e) => problems.push(format!("{}: {point} -> {mine}: back: {e}", job.name)),
            }
            match job.signposts.point_on_triangulation(point) {
                Ok(theirs) => {
                    if let (
                        SurfacePoint::Face { face: f1, bary: b1 },
                        SurfacePoint::Face { face: f2, bary: b2 },
                    ) = (mine, theirs)
                        && f1 == f2
                    {
                        let d = (0..3).map(|k| (b1[k] - b2[k]).abs()).fold(0.0, f64::max);
                        worst_bary = worst_bary.max(d);
                    }
                    if let Ok(back) = job.signposts.point_on_input(theirs) {
                        worst_theirs = worst_theirs
                            .max(dist(position(shape, back), position(shape, point)) / size);
                    }
                }
                Err(_) => skipped += 1,
            }
        }
        (problems, worst_back, worst_theirs, worst_bary, n, skipped)
    });
    let (mut back, mut theirs, mut bary) = (0.0f64, 0.0f64, 0.0f64);
    let (mut n, mut skipped) = (0, 0);
    let mut problems = Vec::new();
    for (p, a, b, c, m, s) in results {
        problems.extend(p);
        back = back.max(a);
        theirs = theirs.max(b);
        bary = bary.max(c);
        n += m;
        skipped += s;
    }
    eprintln!(
        "point queries: {n} points, {skipped} the signposts could not locate; round trip through the signposts' point_on_input: worst {back:.2e} edge lengths (signposts' own result: {theirs:.2e}); \
         same face, barycentric difference {bary:.2e}"
    );
    for p in problems.iter().take(20) {
        eprintln!("  PROBLEM {p}");
    }
    assert!(problems.is_empty());
}

#[test]
fn point_queries_refuse_bad_points() {
    let job = &jobs()[0];
    let t = &job.tri;
    let mesh = t.input().mesh();
    let nan = f64::NAN;
    let bad = [
        SurfacePoint::Vertex(mesh.n_vertices() as u32),
        SurfacePoint::Edge {
            edge: mesh.n_edges() as u32,
            t: 0.5,
        },
        SurfacePoint::Edge { edge: 0, t: nan },
        SurfacePoint::Edge { edge: 0, t: 1.5 },
        SurfacePoint::Face {
            face: mesh.n_faces() as u32,
            bary: [0.2, 0.3, 0.5],
        },
        SurfacePoint::Face {
            face: 0,
            bary: [nan, 0.3, 0.5],
        },
        SurfacePoint::Face {
            face: 0,
            bary: [0.5, 0.5, 0.5],
        },
    ];
    for p in bad {
        assert!(t.point_on_triangulation(p).is_err(), "{p}");
    }
    assert_eq!(
        t.point_on_triangulation(SurfacePoint::Vertex(1)).unwrap(),
        SurfacePoint::Vertex(1)
    );
    assert_eq!(
        t.point_on_triangulation(SurfacePoint::Edge { edge: 0, t: 0.0 })
            .unwrap(),
        SurfacePoint::Vertex(mesh.edge_vertices(0)[0] as u32)
    );
}

/// Corrupted integers never panic: the subdivision and the queries return a result.
#[test]
fn bad_data_never_panics() {
    let jobs = jobs();
    let picks: Vec<&Job> = jobs.iter().step_by(5).collect();
    let results = par_map(&picks, |job| {
        let mut rng = SplitMix64::new(3);
        let (mut ok, mut err) = (0, 0);
        for round in 0..30 {
            let Parts {
                input,
                complex,
                lengths,
                mut normal,
                mut roundabouts,
                locations,
            } = job.tri.clone().into_parts();
            match round % 3 {
                0 => {
                    for _ in 0..1 + round / 3 {
                        let e = (rng.next_u64() % normal.len() as u64) as usize;
                        normal[e] += [-2, -1, 1, 2][(rng.next_u64() % 4) as usize];
                        normal[e] = normal[e].max(-1);
                    }
                }
                1 => {
                    for _ in 0..1 + round / 3 {
                        let h = (rng.next_u64() % roundabouts.len() as u64) as usize;
                        if roundabouts[h] != u32::MAX {
                            roundabouts[h] = roundabouts[h].wrapping_add(1);
                        }
                    }
                }
                _ => {
                    let e = (rng.next_u64() % normal.len() as u64) as usize;
                    normal[e] = (rng.next_u64() % 50) as i64;
                }
            }
            let Ok(t) = IntegerTriangulation::from_parts(Parts {
                input,
                complex,
                lengths,
                normal,
                roundabouts,
                locations,
            }) else {
                err += 1;
                continue;
            };
            match t.common_subdivision() {
                Ok(s) => {
                    ok += 1;
                    let _ = t.interpolation_matrices(&s);
                }
                Err(_) => err += 1,
            }
            let _ = t.point_on_triangulation(SurfacePoint::Face {
                face: 0,
                bary: [0.3, 0.3, 0.4],
            });
        }
        (ok, err)
    });
    let (ok, err) = results.iter().fold((0, 0), |a, r| (a.0 + r.0, a.1 + r.1));
    eprintln!("corrupted states: {ok} gave a subdivision, {err} an error");
    assert!(err > 0);
}

/// Timings: single thread, release, run alone (`--ignored --nocapture`).
#[test]
#[ignore]
fn timings() {
    use std::time::Instant;
    for shape in shapes()
        .iter()
        .filter(|s| s.name.contains("spot") || s.name.contains("camelbad"))
    {
        let mut sp = shape.triangulation();
        sp.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let (tri, _) = convert(&sp).unwrap();
        let best = |f: &mut dyn FnMut()| {
            (0..5)
                .map(|_| {
                    let t0 = Instant::now();
                    f();
                    t0.elapsed().as_secs_f64() * 1e3
                })
                .fold(f64::MAX, f64::min)
        };
        let t_sp = best(&mut || {
            let traces = sp.trace_edges(&EdgeTraceOptions::default()).unwrap();
            std::hint::black_box(sp.overlay(&traces).unwrap());
        });
        let t_ic = best(&mut || {
            std::hint::black_box(tri.common_subdivision().unwrap());
        });
        let t_batch = best(&mut || {
            std::hint::black_box(tri.extract_all_edges().unwrap());
        });
        let mesh = tri.input().mesh();
        let t_seq = best(&mut || {
            for e in 0..mesh.n_edges() {
                std::hint::black_box(tri.extract_edge(mesh.edge_halfedge(e)).unwrap());
            }
        });
        let sub = tri.common_subdivision().unwrap();
        let t_p = best(&mut || {
            std::hint::black_box(tri.interpolation_matrices(&sub).unwrap());
        });
        let mut locator_time = 0.0;
        {
            let mut rng = SplitMix64::new(1);
            let t0 = Instant::now();
            let mut loc = tri.point_locator();
            for _ in 0..1000 {
                let f = (rng.next_u64() % mesh.n_faces() as u64) as u32;
                let _ = loc.point_on_triangulation(SurfacePoint::Face {
                    face: f,
                    bary: [0.3, 0.3, 0.4],
                });
            }
            locator_time += t0.elapsed().as_secs_f64() * 1e3;
        }
        eprintln!(
            "{}: |F0| {}, |V^S| {}, pieces {}: signposts trace+overlay {t_sp:.2} ms; integer subdivision {t_ic:.2} ms \
             (extract_all_edges {t_batch:.2} ms, one by one {t_seq:.2} ms); P0+P1 {t_p:.2} ms; 1000 point queries {locator_time:.1} ms",
            shape.name,
            mesh.n_faces(),
            sub.overlay.vertices.len(),
            sub.overlay.n_faces()
        );
    }
}
