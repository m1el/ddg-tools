//! Replay of the Lean oracle's fixtures for this paper (`lean/oracle/fixtures/intrinsic`).
//!
//! The oracle is an independent implementation of single operations, written from the paper's
//! pseudocode and `SPEC.md` without this crate's code; this test is written from the oracle's
//! `README.md` (section on the third paper) without the oracle's code. Every case is run here
//! and everything its output gives is compared: the input's angle sums, the triangulation
//! before and after a flip or Delaunay flipping (faces as a set, each up to rotation, with
//! lengths and signposts), the unfolded diamond of a flip, refusals, the report and measures
//! of Delaunay flipping, traces from vertices over the input mesh with every crossing, the
//! trace of every halfedge of a flipped triangulation, and circumcentres with their traces
//! from the barycentre.
//!
//! Run with `cargo test --release -p ddg-intrinsic --test oracle -- --nocapture` for the
//! table.

use std::collections::HashMap;
use std::f64::consts::TAU;
use std::path::PathBuf;

use ddg_core::generators::SplitMix64;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::formulas::{
    circumcenter_barycentric, corner_distances, triangle_layout, unfold_diamond,
};
use ddg_intrinsic::{
    DelaunayCriterion, DelaunayOptions, DelaunayStatus, EdgeTraceOptions, FlipLimit, FlipRefusal,
    IntrinsicError, IntrinsicTriangulation, SurfacePoint, Trace, TraceOptions, TraceStop,
};

/// The tolerance of the comparison: that of the oracle's own check of its fixtures. Lengths
/// are compared relative to their size, angles, barycentric coordinates and edge parameters
/// absolutely.
const TOLERANCE: f64 = 1e-9;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/intrinsic")
}

// ----- the text format (oracle README, "Text format") -----

type Line = Vec<String>;

fn read_lines(path: &std::path::Path) -> Vec<Line> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .map(|line| {
            line.split_whitespace()
                .map(str::to_string)
                .collect::<Line>()
        })
        .filter(|tokens| !tokens.is_empty() && !tokens[0].starts_with('#'))
        .collect()
}

fn float(token: &str) -> f64 {
    assert_eq!(token.len(), 16, "a float is 16 hexadecimal digits: {token}");
    f64::from_bits(u64::from_str_radix(token, 16).unwrap_or_else(|e| panic!("{token}: {e}")))
}

fn int(token: &str) -> usize {
    token.parse().unwrap_or_else(|e| panic!("{token}: {e}"))
}

/// A case: the input mesh and the lines of its task.
struct Case {
    name: String,
    task: String,
    mesh: TriMesh,
    geometry: Geometry,
    lines: Vec<Line>,
}

fn read_case(name: &str) -> Case {
    let lines = read_lines(&fixtures().join(format!("{name}.in")));
    let mut task = String::new();
    let mut n_vertices = 0;
    let mut faces = Vec::new();
    let mut edges = Vec::new();
    let mut rest = Vec::new();
    for line in lines {
        match line[0].as_str() {
            "case" => assert_eq!(line[1], name),
            "task" => task = line[1].clone(),
            "mesh" => n_vertices = int(&line[1]),
            "f" => faces.push([1, 2, 3].map(|k| int(&line[k]) as u32)),
            "e" => edges.push((int(&line[1]), int(&line[2]), float(&line[3]))),
            "end" => break,
            _ => rest.push(line),
        }
    }
    let mesh = TriMesh::new(n_vertices, faces).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(edges.len(), mesh.n_edges(), "{name}: every edge once");
    let mut lengths = vec![f64::NAN; mesh.n_edges()];
    for (a, b, length) in edges {
        let e = mesh
            .find_edge(a, b)
            .unwrap_or_else(|| panic!("{name}: no edge {a} {b}"));
        lengths[e] = length;
    }
    let geometry =
        Geometry::from_edge_lengths(&mesh, lengths).unwrap_or_else(|e| panic!("{name}: {e}"));
    Case {
        name: name.to_string(),
        task,
        mesh,
        geometry,
        lines: rest,
    }
}

/// One face of a printed state: vertices, signposts and lengths of its three halfedges.
#[derive(Debug, Clone, Copy)]
struct Corner3 {
    vertices: [u32; 3],
    signposts: [f64; 3],
    lengths: [f64; 3],
}

/// The output of a case.
struct Output {
    lines: Vec<Line>,
}

impl Output {
    fn read(name: &str) -> Output {
        Output {
            lines: read_lines(&fixtures().join(format!("{name}.out"))),
        }
    }

    fn first(&self, tag: &str) -> Option<&Line> {
        self.lines.iter().find(|line| line[0] == tag)
    }

    fn all<'a>(&'a self, tag: &'a str) -> impl Iterator<Item = &'a Line> + 'a {
        self.lines.iter().filter(move |line| line[0] == tag)
    }

    /// The faces of the state with this label, if it was printed.
    fn state(&self, label: &str) -> Option<Vec<Corner3>> {
        let at = self
            .lines
            .iter()
            .position(|line| line[0] == "state" && line[1] == label)?;
        let n_faces = int(&self.lines[at][2]);
        let faces = self.lines[at + 1..at + 1 + n_faces]
            .iter()
            .map(|line| {
                assert_eq!(line[0], "cf");
                Corner3 {
                    vertices: [2, 3, 4].map(|k| int(&line[k]) as u32),
                    signposts: [5, 6, 7].map(|k| float(&line[k])),
                    lengths: [8, 9, 10].map(|k| float(&line[k])),
                }
            })
            .collect();
        Some(faces)
    }
}

// ----- what a replay found -----

/// The largest differences of one case, by kind, and everything that disagrees.
#[derive(Default)]
struct Tally {
    n_compared: usize,
    angle_sum: f64,
    length: f64,
    signpost: f64,
    point: f64,
    parameter: f64,
    other: f64,
    notes: Vec<String>,
    findings: Vec<String>,
}

impl Tally {
    fn check(&mut self, what: &str, difference: f64, kind: fn(&mut Tally) -> &mut f64) {
        self.n_compared += 1;
        let slot = kind(self);
        *slot = slot.max(difference);
        if difference > TOLERANCE || difference.is_nan() {
            self.findings
                .push(format!("{what}: difference {difference:e}"));
        }
    }

    fn equal<T: PartialEq + std::fmt::Debug>(&mut self, what: &str, mine: T, oracle: T) {
        self.n_compared += 1;
        if mine != oracle {
            self.findings
                .push(format!("{what}: here {mine:?}, oracle {oracle:?}"));
        }
    }

    fn relative(&mut self, what: &str, mine: f64, oracle: f64, kind: fn(&mut Tally) -> &mut f64) {
        let scale = oracle.abs().max(mine.abs()).max(f64::MIN_POSITIVE);
        self.check(what, (mine - oracle).abs() / scale, kind);
    }
}

/// The difference of two normalized angles at a vertex: on the circle at an interior vertex.
fn angle_difference(boundary: bool, a: f64, b: f64) -> f64 {
    let d = (a - b).abs();
    if boundary { d } else { d.min((TAU - d).abs()) }
}

fn build(case: &Case) -> IntrinsicTriangulation {
    IntrinsicTriangulation::new(&case.mesh, &case.geometry)
        .unwrap_or_else(|e| panic!("{}: {e}", case.name))
}

/// The header: counts, and the angle sum and kind of every vertex.
fn compare_header(case: &Case, tri: &IntrinsicTriangulation, out: &Output, tally: &mut Tally) {
    let mesh_line = out.first("mesh").unwrap();
    let counts: Vec<i64> = mesh_line[1..].iter().map(|t| t.parse().unwrap()).collect();
    tally.equal(
        "mesh counts (vertices, faces, edges, boundary edges, Euler characteristic)",
        vec![
            case.mesh.n_vertices() as i64,
            case.mesh.n_faces() as i64,
            case.mesh.n_edges() as i64,
            case.mesh.n_boundary_edges() as i64,
            case.mesh.euler_characteristic(),
        ],
        counts,
    );
    let mut n = 0;
    for line in out.all("vt") {
        let v = int(&line[1]);
        tally.check(
            &format!("angle sum of vertex {v}"),
            (tri.vertex_angle_sum(v) - float(&line[2])).abs(),
            |t| &mut t.angle_sum,
        );
        tally.equal(
            &format!("vertex {v} on the boundary"),
            tri.is_boundary_vertex(v),
            line[3] == "1",
        );
        n += 1;
    }
    tally.equal("number of vt lines", n, case.mesh.n_vertices());
}

fn faces_of(tri: &IntrinsicTriangulation) -> Vec<Corner3> {
    (0..tri.n_faces())
        .map(|f| Corner3 {
            vertices: tri.faces()[f],
            signposts: [0, 1, 2].map(|k| tri.signpost(3 * f + k)),
            lengths: tri.face_lengths(f),
        })
        .collect()
}

/// The smallest rotation of a vertex triple.
fn rotation_key(v: [u32; 3]) -> [u32; 3] {
    (0..3)
        .map(|r| [v[r], v[(r + 1) % 3], v[(r + 2) % 3]])
        .min()
        .unwrap()
}

/// A printed state against the triangulation: the faces as a set, each up to rotation.
fn compare_state(tri: &IntrinsicTriangulation, out: &Output, label: &str, tally: &mut Tally) {
    let Some(theirs) = out.state(label) else {
        tally
            .findings
            .push(format!("no state {label} in the output"));
        return;
    };
    let mine = faces_of(tri);
    tally.equal(
        &format!("state {label}: number of faces"),
        mine.len(),
        theirs.len(),
    );
    let mut by_key: HashMap<[u32; 3], Vec<usize>> = HashMap::new();
    for (f, face) in mine.iter().enumerate() {
        by_key
            .entry(rotation_key(face.vertices))
            .or_default()
            .push(f);
    }
    let mut used = vec![false; mine.len()];
    for (n, face) in theirs.iter().enumerate() {
        // The face of this triangulation with these vertices that fits best, and its rotation.
        let mut best: Option<(f64, f64, usize)> = None;
        for &f in by_key
            .get(&rotation_key(face.vertices))
            .into_iter()
            .flatten()
        {
            if used[f] {
                continue;
            }
            for r in 0..3 {
                let rotated = [0, 1, 2].map(|k| (k + r) % 3);
                if rotated.map(|k| mine[f].vertices[k]) != face.vertices {
                    continue;
                }
                let mut length = 0.0f64;
                let mut signpost = 0.0f64;
                for (k, &m) in rotated.iter().enumerate() {
                    length =
                        length.max((mine[f].lengths[m] - face.lengths[k]).abs() / face.lengths[k]);
                    signpost = signpost.max(angle_difference(
                        tri.is_boundary_vertex(face.vertices[k] as usize),
                        mine[f].signposts[m],
                        face.signposts[k],
                    ));
                }
                if best.is_none_or(|(l, s, _)| length.max(signpost) < l.max(s)) {
                    best = Some((length, signpost, f));
                }
            }
        }
        match best {
            Some((length, signpost, f)) => {
                used[f] = true;
                let what = format!("state {label}: face {n} {:?}", face.vertices);
                tally.check(&format!("{what}, lengths"), length, |t| &mut t.length);
                tally.check(&format!("{what}, signposts"), signpost, |t| &mut t.signpost);
            }
            None => tally.findings.push(format!(
                "state {label}: the oracle's face {n} {:?} is not a face here",
                face.vertices
            )),
        }
    }
    let left = used.iter().filter(|u| !**u).count();
    if left > 0 && mine.len() == theirs.len() {
        tally.findings.push(format!(
            "state {label}: {left} faces here are not in the oracle's state"
        ));
    }
}

/// The one edge of the triangulation between two vertices.
fn edge_between(tri: &IntrinsicTriangulation, a: usize, b: usize) -> usize {
    let found: Vec<usize> = (0..tri.n_edges())
        .filter(|&e| {
            let [p, q] = tri.edge_vertices(e);
            (p, q) == (a, b) || (p, q) == (b, a)
        })
        .collect();
    assert_eq!(found.len(), 1, "the pair {a} {b} must identify one edge");
    found[0]
}

fn refusal_name(reason: &FlipRefusal) -> &'static str {
    match reason {
        FlipRefusal::BoundaryEdge => "boundary",
        FlipRefusal::DegreeOneEndpoint { .. } => "degree_one",
        FlipRefusal::NotConvex { .. } => "not_convex",
        FlipRefusal::TooThin { .. } => "too_thin",
    }
}

// ----- flip -----

fn replay_flip(case: &Case, out: &Output, tally: &mut Tally) {
    let mut tri = build(case);
    compare_header(case, &tri, out, tally);
    let flip = case.lines.iter().find(|l| l[0] == "flip").unwrap();
    let again = case.lines.iter().any(|l| l[0] == "again" && l[1] == "1");
    let e = edge_between(&tri, int(&flip[1]), int(&flip[2]));
    let status = out.first("status").unwrap();
    if out.state("before").is_some() {
        compare_state(&tri, out, "before", tally);
    } else {
        tally
            .notes
            .push("no state before in the output of a refused flip".to_string());
    }
    let before = tri.clone();
    match (tri.flip_edge(e), status[1].as_str()) {
        (Err(IntrinsicError::FlipRefused { reason, .. }), "refused") => {
            tally.equal(
                "the reason of the refusal",
                refusal_name(&reason),
                status[2].as_str(),
            );
            tally.equal("a refused flip changes nothing", tri == before, true);
        }
        (Ok(done), "ok") => {
            let info = out.first("flipinfo").unwrap();
            let [i, j, k, l] = [1, 2, 3, 4].map(|n| int(&info[n]) as u32);
            // The halfedge from i to j, and the apex k on its left, l on its right.
            let h = {
                let h = before.edge_halfedge(e);
                if before.from_vertex(h) as u32 == i {
                    h
                } else {
                    before.twin(h).unwrap()
                }
            };
            let t = before.twin(h).unwrap();
            tally.equal(
                "the diamond i j k l of the flip",
                [
                    before.from_vertex(h) as u32,
                    before.to_vertex(h) as u32,
                    before.from_vertex(before.prev(h)) as u32,
                    before.from_vertex(before.prev(t)) as u32,
                ],
                [i, j, k, l],
            );
            let mut ends = done.new_vertices;
            ends.sort_unstable();
            let mut theirs = [k, l];
            theirs.sort_unstable();
            tally.equal("the end points of the flipped edge", ends, theirs);
            let length = |h: usize| before.edge_length(before.edge(h));
            // a = l_ij, b = l_ik, c = l_jk, d = l_il, e = l_jl.
            let (a, b, c, d, e5) = (
                length(h),
                length(before.prev(h)),
                length(before.next(h)),
                length(before.next(t)),
                length(before.prev(t)),
            );
            for (name, mine, n) in [
                ("a", a, 5),
                ("b", b, 6),
                ("c", c, 7),
                ("d", d, 8),
                ("e", e5, 9),
            ] {
                tally.relative(
                    &format!("flipinfo length {name}"),
                    mine,
                    float(&info[n]),
                    |t| &mut t.length,
                );
            }
            let diamond = unfold_diamond(a, c, b, d, e5).unwrap();
            for (name, mine, n) in [
                ("x_k", diamond.x_k, 10),
                ("h_k", diamond.h_k, 11),
                ("x_l", diamond.x_l, 12),
                ("h_l", diamond.h_l, 13),
            ] {
                tally.check(
                    &format!("flipinfo {name}"),
                    (mine - float(&info[n])).abs() / a,
                    |t| &mut t.other,
                );
            }
            tally.relative(
                "the length of the flipped edge, from the diamond",
                diamond.flipped_length(),
                float(&info[14]),
                |t| &mut t.length,
            );
            tally.relative(
                "the length of the flipped edge, as stored",
                done.new_length,
                float(&info[14]),
                |t| &mut t.length,
            );
            compare_state(&tri, out, "after", tally);
            if again {
                match (tri.flip_edge(e), out.first("again")) {
                    (Ok(_), None) => compare_state(&tri, out, "again", tally),
                    (Err(IntrinsicError::FlipRefused { reason, .. }), Some(line)) => {
                        tally.equal("the second flip's refusal", refusal_name(&reason), &line[2]);
                    }
                    (mine, theirs) => tally
                        .findings
                        .push(format!("the second flip: here {mine:?}, oracle {theirs:?}")),
                }
            }
        }
        (mine, theirs) => tally
            .findings
            .push(format!("flip: here {mine:?}, oracle status {theirs}")),
    }
}

// ----- Delaunay flipping -----

fn delaunay_options(case: &Case) -> DelaunayOptions {
    let value = |tag: &str| {
        case.lines
            .iter()
            .find(|l| l[0] == tag)
            .map(|l| l[1].clone())
    };
    let tolerance = value("tolerance").map_or(1e-10, |t| float(&t));
    let limit = value("limit").map_or(64, |t| int(&t));
    let criterion = match value("criterion").as_deref() {
        // The reference's test has the fixed tolerance 1e-6 on the weight.
        Some("cotan") => DelaunayCriterion::REFERENCE,
        _ => DelaunayCriterion::AngleSum { tolerance },
    };
    DelaunayOptions {
        criterion,
        max_flips: FlipLimit::PerEdge(limit),
    }
}

/// The queue of SPEC N6 with a given first order of the edges, through the public interface:
/// the first edge is taken off and flipped if it fails the test, and the four other edges of
/// its two triangles are put back unless they are queued. Returns the number of flips.
fn plain_queue(
    mut tri: IntrinsicTriangulation,
    order: Vec<usize>,
    options: &DelaunayOptions,
) -> usize {
    let fails = |tri: &IntrinsicTriangulation, e: usize| match options.criterion {
        DelaunayCriterion::AngleSum { tolerance } => {
            tri.edge_angle_excess(e).is_some_and(|x| x > tolerance)
        }
        DelaunayCriterion::CotanWeight { tolerance } => {
            !tri.is_boundary_edge(e) && tri.edge_cotan_weight(e) < -tolerance
        }
    };
    let mut queued = vec![true; tri.n_edges()];
    let mut queue: std::collections::VecDeque<usize> = order.into();
    let mut n_flips = 0;
    while let Some(e) = queue.pop_front() {
        queued[e] = false;
        if !fails(&tri, e) || tri.flip_edge(e).is_err() {
            continue;
        }
        n_flips += 1;
        let h = tri.edge_halfedge(e);
        let t = tri.twin(h).unwrap();
        for g in [tri.next(h), tri.prev(h), tri.next(t), tri.prev(t)] {
            let neighbour = tri.edge(g);
            if !std::mem::replace(&mut queued[neighbour], true) {
                queue.push_back(neighbour);
            }
        }
    }
    n_flips
}

fn replay_delaunay(case: &Case, out: &Output, tally: &mut Tally) {
    let mut tri = build(case);
    compare_header(case, &tri, out, tally);
    compare_state(&tri, out, "input", tally);
    let report = tri.flip_to_delaunay(&delaunay_options(case)).unwrap();
    let theirs = out.first("report").unwrap();
    // The number of flips is informational on the oracle's side, and rightly: it depends on
    // the order in which the queue first holds the edges (here by edge id, there by halfedge),
    // while the triangulation it ends with does not. The queue of SPEC N6 replayed through
    // the public interface, in the order of the edge ids and in 200 random orders:
    let options = delaunay_options(case);
    let in_order: Vec<usize> = (0..tri.n_edges()).collect();
    tally.equal(
        "flips of the plain queue in the order of the edge ids",
        plain_queue(build(case), in_order.clone(), &options),
        report.n_flips,
    );
    let mut rng = SplitMix64::new(0x0c1e);
    let (mut fewest, mut most) = (usize::MAX, 0);
    for _ in 0..200 {
        let mut order = in_order.clone();
        for i in (1..order.len()).rev() {
            order.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let n = plain_queue(build(case), order, &options);
        fewest = fewest.min(n);
        most = most.max(n);
    }
    let oracle_flips = int(&theirs[1]);
    tally.notes.push(format!(
        "flips: {} here, {oracle_flips} in the oracle; {fewest} to {most} over 200 random orders of the queue",
        report.n_flips
    ));
    tally.equal(
        "the flips of both lie in the range of the random orders",
        (fewest..=most).contains(&report.n_flips) && (fewest..=most).contains(&oracle_flips),
        true,
    );
    tally.equal("refused flips", report.n_refused, int(&theirs[2]));
    tally.equal("ties", report.n_ties, int(&theirs[3]));
    tally.equal(
        "status",
        match report.status {
            DelaunayStatus::Delaunay => "done",
            DelaunayStatus::FlipLimit => "flip_limit",
            DelaunayStatus::Blocked => "blocked",
        },
        theirs[4].as_str(),
    );
    let measures = out.first("measures").unwrap();
    tally.check(
        "largest angle excess",
        (report.max_angle_excess - float(&measures[1])).abs(),
        |t| &mut t.other,
    );
    // SPEC N6: the cotan weight of an edge is w = ½ (cot α + cot β). The oracle's MINWEIGHT is
    // the sum of the two cotangents, 2 w (the formula tagged `it-cotan-weight` on both
    // sides). Both are accepted here, so that the fixtures can be corrected to print w.
    let theirs_weight = float(&measures[2]);
    let as_weight = (report.min_cotan_weight - theirs_weight).abs() / theirs_weight.abs();
    let as_sum = (2.0 * report.min_cotan_weight - theirs_weight).abs() / theirs_weight.abs();
    tally.notes.push(format!(
        "smallest cotan weight {:.6} here; the oracle prints {:.6}, {}",
        report.min_cotan_weight,
        theirs_weight,
        if as_weight <= as_sum {
            "the weight"
        } else {
            "twice the weight (cot + cot)"
        }
    ));
    tally.check("smallest cotan weight", as_weight.min(as_sum), |t| {
        &mut t.other
    });
    tally.relative("total area", tri.total_area(), float(&measures[3]), |t| {
        &mut t.other
    });
    compare_state(&tri, out, "final", tally);
    tally.notes.push(format!(
        "{} flips, {} self-edges and {} edges in groups afterwards",
        report.n_flips,
        tri.irregularities().self_edges.len(),
        tri.irregularities().multi_edges.len()
    ));
}

// ----- traces -----

/// A trace as the oracle prints it.
struct OracleTrace {
    stop: String,
    face: usize,
    bary: [f64; 3],
    left: f64,
    vertex: i64,
    clamp: bool,
    crossings: Vec<(usize, usize, f64)>,
}

/// Reads `TAG STOP FACE B0 B1 B2 LEFT VERTEX CLAMP NCROSS` and its `TAGx U V T` lines.
fn oracle_trace(out: &Output, tag: &str) -> OracleTrace {
    let line = out
        .first(tag)
        .unwrap_or_else(|| panic!("no line {tag} in the output"));
    let cross = format!("{tag}x");
    let crossings: Vec<(usize, usize, f64)> = out
        .all(&cross)
        .map(|l| (int(&l[1]), int(&l[2]), float(&l[3])))
        .collect();
    assert_eq!(crossings.len(), int(&line[9]), "{tag}: number of crossings");
    OracleTrace {
        stop: line[1].clone(),
        face: int(&line[2]),
        bary: [3, 4, 5].map(|k| float(&line[k])),
        left: float(&line[6]),
        vertex: line[7].parse().unwrap(),
        clamp: line[8] == "1",
        crossings,
    }
}

/// A point of the mesh in the barycentric coordinates of face `f`, if it lies in its closure.
fn in_face(mesh: &TriMesh, point: SurfacePoint, f: usize) -> Option<[f64; 3]> {
    let corners = mesh.faces()[f];
    let slot = |v: usize| corners.iter().position(|&c| c as usize == v);
    let mut bary = [0.0; 3];
    match point {
        SurfacePoint::Vertex(v) => bary[slot(v as usize)?] = 1.0,
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = mesh.edge_vertices(edge as usize);
            bary[slot(a)?] = 1.0 - t;
            bary[slot(b)?] = t;
        }
        SurfacePoint::Face { face, bary: b } => {
            if face as usize != f {
                return None;
            }
            bary = b;
        }
    }
    Some(bary)
}

/// What the oracle calls the way a trace stopped.
fn stop_name(trace: &Trace) -> &'static str {
    match (trace.stop, trace.end) {
        (TraceStop::Boundary, _) => "boundary",
        (TraceStop::Vertex, _) | (TraceStop::Completed, SurfacePoint::Vertex(_)) => "vertex",
        (TraceStop::Completed, SurfacePoint::Edge { .. }) => "edge",
        (TraceStop::Completed, SurfacePoint::Face { .. }) => "end",
    }
}

/// The crossings of a trace as the oracle lists them: the input edge by its vertices and the
/// parameter from the lower one. The oracle lists the boundary edge at which a path leaves
/// the surface as its last crossing; here that edge is the end point, so it is added.
fn crossings_of(mesh: &TriMesh, trace: &Trace) -> Vec<(usize, usize, f64)> {
    let mut list: Vec<(usize, usize, f64)> = trace
        .crossings
        .iter()
        .map(|c| {
            let [a, b] = mesh.edge_vertices(c.edge as usize);
            (a, b, c.t)
        })
        .collect();
    assert_eq!(list.len(), trace.n_crossings);
    if let (TraceStop::Boundary, SurfacePoint::Edge { edge, t }) = (trace.stop, trace.end) {
        let [a, b] = mesh.edge_vertices(edge as usize);
        if list.last().is_none_or(|&(p, q, _)| (p, q) != (a, b)) {
            list.push((a, b, t));
        }
    }
    list
}

fn compare_crossings(
    what: &str,
    mine: &[(usize, usize, f64)],
    theirs: &[(usize, usize, f64)],
    tally: &mut Tally,
) {
    let edges = |list: &[(usize, usize, f64)]| -> Vec<(usize, usize)> {
        list.iter().map(|&(a, b, _)| (a, b)).collect()
    };
    if edges(mine) != edges(theirs) {
        tally.findings.push(format!(
            "{what}: crossed edges here {:?}, oracle {:?}",
            edges(mine),
            edges(theirs)
        ));
        return;
    }
    for (m, o) in mine.iter().zip(theirs) {
        tally.check(
            &format!("{what}: parameter on edge {} {}", o.0, o.1),
            (m.2 - o.2).abs(),
            |t| &mut t.parameter,
        );
    }
}

/// A trace over the input mesh for `requested` against the oracle's.
fn compare_trace(
    case: &Case,
    what: &str,
    trace: &Trace,
    requested: f64,
    theirs: &OracleTrace,
    tally: &mut Tally,
) {
    tally.equal(
        &format!("{what}: stop"),
        stop_name(trace),
        theirs.stop.as_str(),
    );
    match in_face(&case.mesh, trace.end, theirs.face) {
        Some(bary) => {
            let difference = (0..3)
                .map(|k| (bary[k] - theirs.bary[k]).abs())
                .fold(0.0, f64::max);
            tally.check(
                &format!("{what}: end point in face {}", theirs.face),
                difference,
                |t| &mut t.point,
            );
        }
        None => tally.findings.push(format!(
            "{what}: the end {} is not in the oracle's face {}",
            trace.end, theirs.face
        )),
    }
    tally.check(
        &format!("{what}: length left"),
        ((requested - trace.length) - theirs.left).abs() / requested.max(1.0),
        |t| &mut t.other,
    );
    if theirs.vertex >= 0 {
        tally.equal(
            &format!("{what}: the vertex reached"),
            trace.end,
            SurfacePoint::Vertex(theirs.vertex as u32),
        );
    }
    tally.equal(
        &format!("{what}: clamped"),
        trace.clamp != 0.0,
        theirs.clamp,
    );
    compare_crossings(
        what,
        &crossings_of(&case.mesh, trace),
        &theirs.crossings,
        tally,
    );
}

fn recording() -> TraceOptions {
    TraceOptions {
        record_crossings: true,
        ..TraceOptions::default()
    }
}

fn replay_trace(case: &Case, out: &Output, tally: &mut Tally) {
    let tri = build(case);
    compare_header(case, &tri, out, tally);
    let queries: Vec<&Line> = case.lines.iter().filter(|l| l[0] == "trace").collect();
    tally.equal("number of queries", out.all("query").count(), queries.len());
    for (n, query) in queries.iter().enumerate() {
        let (v, phi, distance) = (int(&query[1]), float(&query[2]), float(&query[3]));
        let trace = tri
            .input()
            .trace(SurfacePoint::Vertex(v as u32), phi, distance, &recording())
            .unwrap_or_else(|e| panic!("{}: query {n}: {e}", case.name));
        let theirs = oracle_trace(out, &format!("tr{n}"));
        compare_trace(
            case,
            &format!("query {n} from vertex {v}"),
            &trace,
            distance,
            &theirs,
            tally,
        );
    }
}

// ----- the trace of every halfedge -----

/// The distance of a point of the input mesh from vertex `w`, if both lie in one face.
fn distance_from_vertex(case: &Case, point: SurfacePoint, w: usize) -> Option<f64> {
    let lengths = case.geometry.edge_lengths();
    match point {
        SurfacePoint::Vertex(v) => (v as usize == w).then_some(0.0),
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = case.mesh.edge_vertices(edge as usize);
            let length = lengths[edge as usize];
            if a == w {
                Some(t * length)
            } else if b == w {
                Some((1.0 - t) * length)
            } else {
                None
            }
        }
        SurfacePoint::Face { face, bary } => {
            let f = face as usize;
            let slot = case.mesh.faces()[f].iter().position(|&c| c as usize == w)?;
            let l = [0, 1, 2].map(|k| lengths[case.mesh.edge(3 * f + k)]);
            Some(corner_distances(l, bary)?[slot])
        }
    }
}

/// The crossings of a path to vertex `w` without those of its arrival: where a path ends at a
/// vertex, rounding decides whether it still crosses an edge at that vertex, and the two
/// implementations need not agree on that.
fn before_arrival(list: &[(usize, usize, f64)], w: usize) -> Vec<(usize, usize, f64)> {
    let mut list = list.to_vec();
    while let Some(&(a, b, t)) = list.last() {
        let at_w = (a == w && t < 1e-7) || (b == w && t > 1.0 - 1e-7);
        if !at_w {
            break;
        }
        list.pop();
    }
    list
}

fn replay_edgetrace(case: &Case, out: &Output, tally: &mut Tally) {
    let mut tri = build(case);
    compare_header(case, &tri, out, tally);
    for line in case.lines.iter().filter(|l| l[0] == "preflip") {
        let e = edge_between(&tri, int(&line[1]), int(&line[2]));
        tri.flip_edge(e)
            .unwrap_or_else(|e| panic!("{}: preflip: {e}", case.name));
    }
    if case
        .lines
        .iter()
        .any(|l| l[0] == "predelaunay" && l[1] == "1")
    {
        tri.flip_to_delaunay(&delaunay_options(case)).unwrap();
    }
    compare_state(&tri, out, "state", tally);

    // The oracle's rows, each with its crossings; a row is identified by (U, W, PHI).
    let rows: Vec<&Line> = out.all("te").collect();
    tally.equal("number of traced halfedges", rows.len(), tri.n_halfedges());
    let mut crossings: HashMap<(usize, usize, u64), Vec<(usize, usize, f64)>> = HashMap::new();
    for l in out.all("tx") {
        crossings
            .entry((int(&l[1]), int(&l[2]), float(&l[3]).to_bits()))
            .or_default()
            .push((int(&l[4]), int(&l[5]), float(&l[6])));
    }
    let mut used = vec![false; tri.n_halfedges()];
    let checked = tri
        .trace_edges(&EdgeTraceOptions::default())
        .unwrap_or_else(|e| panic!("{}: {e}", case.name));
    let (mut same_stop, mut arrivals_differ) = (0usize, 0usize);
    for row in rows {
        let (u, w, phi, length) = (int(&row[1]), int(&row[2]), float(&row[3]), float(&row[4]));
        let theirs = crossings
            .get(&(u, w, phi.to_bits()))
            .cloned()
            .unwrap_or_default();
        assert_eq!(
            theirs.len(),
            int(&row[11]),
            "te {u} {w}: number of crossings"
        );
        // The halfedge here: from u to w with that signpost.
        let boundary = tri.is_boundary_vertex(u);
        let h = (0..tri.n_halfedges())
            .filter(|&h| !used[h] && tri.from_vertex(h) == u && tri.to_vertex(h) == w)
            .min_by(|&a, &b| {
                angle_difference(boundary, tri.signpost(a), phi).total_cmp(&angle_difference(
                    boundary,
                    tri.signpost(b),
                    phi,
                ))
            });
        let Some(h) = h else {
            tally
                .findings
                .push(format!("te {u} {w}: no such halfedge here"));
            continue;
        };
        used[h] = true;
        let what = format!("halfedge {u} -> {w}");
        tally.check(
            &format!("{what}: signpost"),
            angle_difference(boundary, tri.signpost(h), phi),
            |t| &mut t.signpost,
        );
        let e = tri.edge(h);
        tally.relative(
            &format!("{what}: length"),
            tri.edge_length(e),
            length,
            |t| &mut t.length,
        );
        // The same trace as the oracle's: over the input mesh, along the signpost, for the
        // length, with no target.
        let trace = tri
            .input()
            .trace(
                SurfacePoint::Vertex(u as u32),
                tri.direction_on_input(u, tri.signpost(h)),
                tri.edge_length(e),
                &recording(),
            )
            .unwrap_or_else(|e| panic!("{}: {what}: {e}", case.name));
        // It arrives: the end is as far from w as the oracle's, which is nowhere.
        let miss = float(&row[10]);
        match distance_from_vertex(case, trace.end, w) {
            Some(distance) => tally.check(
                &format!("{what}: distance of the end from {w}"),
                (distance - miss).abs() / length,
                |t| &mut t.point,
            ),
            None => tally.findings.push(format!(
                "{what}: the path ends at {}, not next to {w}",
                trace.end
            )),
        }
        tally.check(&format!("{what}: the oracle's miss"), miss / length, |t| {
            &mut t.point
        });
        let mine = crossings_of(&case.mesh, &trace);
        same_stop += usize::from(stop_name(&trace) == row[5]);
        arrivals_differ += usize::from(mine.len() != theirs.len());
        compare_crossings(
            &what,
            &before_arrival(&mine, w),
            &before_arrival(&theirs, w),
            tally,
        );
        // And the checked trace of the edge (SPEC N8), for the halfedge that gives the edge
        // its direction: the same crossings, and it arrives.
        if tri.edge_halfedge(e) == h {
            let listed: Vec<(usize, usize, f64)> = checked
                .crossings(e)
                .iter()
                .map(|c| {
                    let [a, b] = case.mesh.edge_vertices(c.edge as usize);
                    (a, b, c.t)
                })
                .collect();
            compare_crossings(
                &format!("{what}, checked"),
                &before_arrival(&listed, w),
                &before_arrival(&theirs, w),
                tally,
            );
            tally.check(
                &format!("{what}, checked: mismatch"),
                checked.mismatch(e),
                |t| &mut t.point,
            );
        }
    }
    tally.notes.push(format!(
        "{} halfedges; the same stop in {same_stop}; {arrivals_differ} differ in a crossing at the far vertex",
        tri.n_halfedges()
    ));
}

// ----- circumcentres -----

fn replay_circumcenter(case: &Case, out: &Output, tally: &mut Tally) {
    let tri = build(case);
    compare_header(case, &tri, out, tally);
    let mut faces: Vec<usize> = Vec::new();
    for line in case.lines.iter().filter(|l| l[0] == "circumcenter") {
        if line[1] == "all" {
            faces.extend(0..case.mesh.n_faces());
        } else {
            faces.push(int(&line[1]));
        }
    }
    let rows: HashMap<usize, &Line> = out.all("cc").map(|l| (int(&l[1]), l)).collect();
    tally.equal("number of circumcentres", rows.len(), faces.len());
    let (mut n_inside, mut n_boundary) = (0usize, 0usize);
    for f in faces {
        let Some(row) = rows.get(&f) else {
            tally.findings.push(format!("no line cc {f}"));
            continue;
        };
        let lengths = case.geometry.edge_lengths();
        let l = [0, 1, 2].map(|k| lengths[case.mesh.edge(3 * f + k)]);
        let weights = circumcenter_barycentric(l).unwrap();
        let theirs = [2, 3, 4].map(|k| float(&row[k]));
        let scale = theirs.iter().fold(1.0f64, |s, w| s.max(w.abs()));
        for k in 0..3 {
            tally.check(
                &format!("face {f}: weight {k} of the circumcentre"),
                (weights[k] - theirs[k]).abs() / scale,
                |t| &mut t.point,
            );
        }
        let radii = corner_distances(l, weights).unwrap();
        for k in 0..3 {
            tally.relative(
                &format!("face {f}: distance of the circumcentre from corner {k}"),
                radii[k],
                float(&row[5 + k]),
                |t| &mut t.length,
            );
        }
        let inside = weights.iter().all(|&w| w > 0.0);
        tally.equal(&format!("face {f}: inside"), inside, row[8] == "1");
        n_inside += usize::from(inside);
        // The trace from the barycentre straight to the circumcentre, in the frame of the
        // face: corner 0 at the origin, corner 1 on the x-axis.
        let layout = triangle_layout(l).unwrap();
        let point = |b: [f64; 3]| {
            [0, 1].map(|i| b[0] * layout[0][i] + b[1] * layout[1][i] + b[2] * layout[2][i])
        };
        let (from, to) = (point([1.0 / 3.0; 3]), point(weights));
        let d = [to[0] - from[0], to[1] - from[1]];
        let distance = d[0].hypot(d[1]);
        let theirs = oracle_trace(out, &format!("cr{f}"));
        if distance == 0.0 {
            tally.equal(&format!("face {f}: stop"), "end", theirs.stop.as_str());
            continue;
        }
        let direction = d[1].atan2(d[0]).rem_euclid(TAU);
        let trace = tri
            .input()
            .trace(
                SurfacePoint::Face {
                    face: f as u32,
                    bary: [1.0 / 3.0; 3],
                },
                if direction >= TAU { 0.0 } else { direction },
                distance,
                &recording(),
            )
            .unwrap_or_else(|e| panic!("{}: face {f}: {e}", case.name));
        n_boundary += usize::from(trace.stop == TraceStop::Boundary);
        compare_trace(
            case,
            &format!("face {f}: trace to the circumcentre"),
            &trace,
            distance,
            &theirs,
            tally,
        );
    }
    tally.notes.push(format!(
        "{n_inside} circumcentres inside their face, {n_boundary} traces reach the boundary"
    ));
}

// ----- the replay -----

#[test]
fn the_oracle_s_fixtures_replay() {
    let manifest = fixtures().join("MANIFEST.tsv");
    let text = std::fs::read_to_string(&manifest)
        .unwrap_or_else(|e| panic!("{}: {e}", manifest.display()));
    let cases: Vec<(String, String)> = text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut columns = line.split('\t');
            (
                columns.next().unwrap().to_string(),
                columns.next().unwrap().to_string(),
            )
        })
        .collect();
    assert_eq!(cases.len(), 18, "the oracle's README lists 18 cases");
    // One job per case, in parallel.
    let tallies: Vec<(String, String, Tally)> = std::thread::scope(|scope| {
        let handles: Vec<_> = cases
            .iter()
            .map(|(name, expected)| {
                scope.spawn(move || {
                    let case = read_case(name);
                    let out = Output::read(name);
                    let mut tally = Tally::default();
                    let status = out.first("status").unwrap();
                    tally.equal(
                        "the status of the manifest",
                        status[1].as_str(),
                        expected.as_str(),
                    );
                    tally.equal(
                        "the task",
                        out.first("task").unwrap()[1].as_str(),
                        case.task.as_str(),
                    );
                    match case.task.as_str() {
                        "flip" => replay_flip(&case, &out, &mut tally),
                        "delaunay" => replay_delaunay(&case, &out, &mut tally),
                        "trace" => replay_trace(&case, &out, &mut tally),
                        "edgetrace" => replay_edgetrace(&case, &out, &mut tally),
                        "circumcenter" => replay_circumcenter(&case, &out, &mut tally),
                        other => panic!("{name}: unknown task {other}"),
                    }
                    (name.clone(), case.task.clone(), tally)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    println!(
        "largest differences per case (tolerance {TOLERANCE:e}); lengths relative, the rest absolute"
    );
    println!(
        "{:<28} {:<12} {:>8} | {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} | {:>8}",
        "case",
        "task",
        "compared",
        "angle sum",
        "length",
        "signpost",
        "point",
        "parameter",
        "other",
        "findings"
    );
    let mut n_findings = 0;
    for (name, task, tally) in &tallies {
        println!(
            "{:<28} {:<12} {:>8} | {:>9.1e} {:>9.1e} {:>9.1e} {:>9.1e} {:>9.1e} {:>9.1e} | {:>8}",
            name,
            task,
            tally.n_compared,
            tally.angle_sum,
            tally.length,
            tally.signpost,
            tally.point,
            tally.parameter,
            tally.other,
            tally.findings.len()
        );
        for note in &tally.notes {
            println!("    {note}");
        }
        for finding in tally.findings.iter().take(12) {
            println!("    FINDING: {finding}");
        }
        if tally.findings.len() > 12 {
            println!("    ... and {} more", tally.findings.len() - 12);
        }
        n_findings += tally.findings.len();
    }
    assert_eq!(n_findings, 0, "the replay disagrees with the oracle");
}

/// The oracle's reading of a direction that equals a signpost exactly: the path runs along
/// the edge, and ends on it (`edge`) if the distance is shorter than the edge, at the far
/// vertex otherwise (`vertex`, with the length that is left). There is no fixture for the
/// short and the long case, so they are checked here on a fixture's mesh, for every halfedge.
/// At a boundary vertex the direction `π` is the boundary edge that arrives (SPEC N2).
#[test]
fn a_direction_equal_to_a_signpost_runs_along_the_edge() {
    for name in ["it-trace-grid", "it-trace-ico2"] {
        let case = read_case(name);
        let tri = build(&case);
        let input = tri.input();
        let mesh = &case.mesh;
        let lengths = case.geometry.edge_lengths();
        let along = |v: usize, phi: f64, edge: usize, far: usize| {
            let length = lengths[edge];
            let from_lower = mesh.edge_vertices(edge)[0] == v;
            for (factor, expected) in [(0.25, "edge"), (1.0, "vertex"), (1.75, "vertex")] {
                let distance = factor * length;
                let trace = input
                    .trace(SurfacePoint::Vertex(v as u32), phi, distance, &recording())
                    .unwrap_or_else(|e| panic!("{name}: vertex {v}, direction {phi}: {e}"));
                let what = format!("{name}: from {v} along edge {edge} for {factor} of it");
                assert_eq!(stop_name(&trace), expected, "{what}");
                assert_eq!(trace.along_edge, Some(edge as u32), "{what}");
                assert_eq!(trace.n_crossings, 0, "{what}");
                assert_eq!(trace.clamp, 0.0, "{what}");
                if factor < 1.0 {
                    // A point of the edge, at that distance from the start.
                    let t = if from_lower { factor } else { 1.0 - factor };
                    match trace.end {
                        SurfacePoint::Edge { edge: e, t: at } => {
                            assert_eq!(e as usize, edge, "{what}");
                            assert!((at - t).abs() < 1e-15, "{what}: {at}");
                        }
                        other => panic!("{what}: {other}"),
                    }
                    assert_eq!(trace.stop, TraceStop::Completed, "{what}");
                    assert_eq!(trace.length, distance, "{what}");
                } else {
                    // The far vertex; what is left of the distance is not travelled.
                    assert_eq!(trace.end, SurfacePoint::Vertex(far as u32), "{what}");
                    assert_eq!(trace.length, length, "{what}");
                    assert_eq!(
                        trace.stop,
                        if factor > 1.0 {
                            TraceStop::Vertex
                        } else {
                            TraceStop::Completed
                        },
                        "{what}"
                    );
                }
            }
        };
        for h in 0..mesh.n_halfedges() {
            along(
                mesh.from_vertex(h),
                input.signpost(h),
                mesh.edge(h),
                mesh.to_vertex(h),
            );
        }
        // The boundary edge that arrives at a boundary vertex: the direction π.
        let mut n_boundary = 0;
        for h in 0..mesh.n_halfedges() {
            if mesh.twin(h).is_none() {
                let v = mesh.to_vertex(h);
                along(v, std::f64::consts::PI, mesh.edge(h), mesh.from_vertex(h));
                n_boundary += 1;
            }
        }
        assert_eq!(n_boundary, mesh.n_boundary_edges());
    }
}
