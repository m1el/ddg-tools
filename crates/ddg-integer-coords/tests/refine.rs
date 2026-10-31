//! Delaunay refinement on integer coordinates (WP6, SPEC §5): the shared driver against the
//! signpost representation, the paper's angle bounds, termination.
//!
//! Run with `--nocapture` for the tables; `timing` is ignored (single thread, pinned).

mod common;

use std::collections::HashMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use common::{Scenario, Shape, all_shapes, bad_meshes, convert, load_data, models, state};
use ddg_core::generators;
use ddg_integer_coords::{
    IntegerError, IntegerRefineOptions, IntegerRefineReport, IntegerTriangulation,
};
use ddg_intrinsic::{
    BelowReason, Complex, EdgeTrace, EdgeTraceOptions, Flip, FlipRefusal, InputSurface,
    InsertedVertexRemoval, Insertion, InsertionLimit, IntrinsicError, IntrinsicTriangulation,
    RefineOptions, RefineReport, RefineRules, RefineStatus, SurfacePoint, Trace, TraceOptions,
    Triangulation,
};

fn par_map<S: Sync, R: Send>(specs: &[S], f: impl Fn(&S) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(16);
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<Option<R>>> = Mutex::new((0..specs.len()).map(|_| None).collect());
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= specs.len() {
                        break;
                    }
                    let r = f(&specs[i]);
                    out.lock().unwrap()[i] = Some(r);
                }
            });
        }
    });
    out.into_inner()
        .unwrap()
        .into_iter()
        .map(|r| r.unwrap())
        .collect()
}

fn options(degrees: f64) -> RefineOptions {
    RefineOptions {
        min_angle: degrees.to_radians(),
        ..Default::default()
    }
}

fn integer_options(degrees: f64, rules: RefineRules) -> IntegerRefineOptions {
    IntegerRefineOptions {
        refine: options(degrees),
        rules,
        ..Default::default()
    }
}

fn rules_name(rules: RefineRules) -> &'static str {
    match rules {
        RefineRules::Ruppert => "Ruppert",
        RefineRules::Chew => "Chew",
    }
}

// ----- shapes -----

/// A closed double cone over a regular `n`-gon of radius 1: the top apex has the angle sum
/// `top` degrees, the bottom apex `bottom` degrees, the ring vertices whatever follows.
fn bicone(n: usize, top: f64, bottom: f64) -> Shape {
    let chord = 2.0 * (std::f64::consts::PI / n as f64).sin();
    let height = |degrees: f64| {
        let corner = (degrees / n as f64).to_radians();
        let slant = chord / (2.0 * (0.5 * corner).sin());
        (slant * slant - 1.0).sqrt()
    };
    let mut positions: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            [a.cos(), a.sin(), 0.0]
        })
        .collect();
    positions.push([0.0, 0.0, height(top)]);
    positions.push([0.0, 0.0, -height(bottom)]);
    let (t, b) = (n as u32, n as u32 + 1);
    let mut faces = Vec::new();
    for i in 0..n as u32 {
        let j = (i + 1) % n as u32;
        faces.push([i, j, t]);
        faces.push([j, i, b]);
    }
    Shape::new(
        &format!("bicone {n} {top}/{bottom}"),
        positions,
        faces,
        false,
    )
}

/// An open cone (a disk with a cone point): `n` triangles around an apex of angle sum
/// `degrees`, rings of `rings` vertices, boundary at the outer ring.
fn open_cone(n: usize, degrees: f64) -> Shape {
    let corner = (degrees / n as f64).to_radians();
    // Lay the cone out flat: the apex at the origin, ring vertices at radius r along a fan of
    // total angle `degrees`, then roll it up into space.
    let chord = 2.0 * (std::f64::consts::PI / n as f64).sin();
    let slant = chord / (2.0 * (0.5 * corner).sin());
    let height = (slant * slant - 1.0).max(0.0).sqrt();
    let mut positions = vec![[0.0, 0.0, height]];
    for ring in 1..=2u32 {
        let s = ring as f64 / 2.0;
        for i in 0..n {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            positions.push([s * a.cos(), s * a.sin(), height * (1.0 - s)]);
        }
    }
    let v = |ring: u32, i: usize| 1 + (ring - 1) * n as u32 + (i % n) as u32;
    let mut faces = Vec::new();
    for i in 0..n {
        faces.push([0, v(1, i), v(1, i + 1)]);
        faces.push([v(1, i), v(2, i), v(2, i + 1)]);
        faces.push([v(1, i), v(2, i + 1), v(1, i + 1)]);
    }
    Shape::new(&format!("open cone {n} {degrees}"), positions, faces, false)
}

/// Cones with small angle sums (narrow vertices, IC-SPEC I23) and one above 60°.
fn cone_shapes() -> Vec<Shape> {
    vec![
        bicone(6, 20.0, 200.0),
        bicone(6, 45.0, 100.0),
        bicone(8, 59.0, 59.0),
        bicone(8, 70.0, 300.0),
        open_cone(6, 30.0),
        open_cone(8, 50.0),
    ]
}

/// More meshes with boundary.
fn boundary_shapes() -> Vec<Shape> {
    let mut out = vec![Shape::generated(
        "bumpy grid 20x20",
        generators::bumpy_grid_disk(20, 20, 0.3).unwrap(),
        false,
    )];
    out.extend(load_data("optimal-cones/optimal-cones-data/data/bump.obj"));
    out.extend(load_data("optimal-cones/optimal-cones-data/data/hand.obj"));
    out
}

/// Every shape of the differential test.
fn shapes() -> Vec<Shape> {
    let mut out = all_shapes();
    out.extend(cone_shapes());
    out.extend(boundary_shapes());
    out.extend(
        models()
            .into_iter()
            .filter(|s| s.name != "spot_triangulated"),
    );
    out.extend(
        bad_meshes()
            .into_iter()
            .filter(|s| s.name != "camelbad")
            .take(2),
    );
    out
}

// ----- a log of the operations the driver makes -----

#[derive(Debug, Clone, PartialEq)]
enum Op {
    Flip(u32, bool),
    Insert {
        face: u32,
        bary: [u64; 3],
        result: Result<u32, String>,
    },
    Split {
        edge: u32,
        t: u64,
        result: Result<u32, String>,
    },
    Remove(u32, Result<usize, String>),
}

/// What an operation's error was, for the log.
fn kind(e: &IntrinsicError) -> String {
    match e {
        IntrinsicError::InsertionRefused { reason, .. } => format!("refused: {reason}"),
        IntrinsicError::RemovalRefused { reason, .. } => format!("refused: {reason}"),
        IntrinsicError::FlipRefused { .. } => "flip refused".to_string(),
        other => format!("{other}"),
    }
}

/// A representation whose operations are logged; integer face splits note their distance to
/// the nearest arc (I6), signpost insertions whether the place was snapped onto an input edge.
struct Log<'a, T> {
    t: &'a mut T,
    ops: Vec<Op>,
    notes: HashMap<usize, String>,
    note: fn(&T, usize, [f64; 3]) -> Option<String>,
    integer: bool,
}

impl<T: Triangulation> Triangulation for Log<'_, T> {
    fn complex(&self) -> &Complex {
        self.t.complex()
    }
    fn input(&self) -> &InputSurface {
        self.t.input()
    }
    fn edge_lengths(&self) -> &[f64] {
        self.t.edge_lengths()
    }
    fn vertex_angle_sum(&self, v: usize) -> f64 {
        self.t.vertex_angle_sum(v)
    }
    fn vertex_location(&self, v: usize) -> SurfacePoint {
        self.t.vertex_location(v)
    }
    fn can_flip(&self, e: usize) -> Result<Result<(), FlipRefusal>, IntrinsicError> {
        self.t.can_flip(e)
    }
    fn flip_edge(&mut self, e: usize) -> Result<Flip, IntrinsicError> {
        let r = self.t.flip_edge(e);
        self.ops.push(Op::Flip(e as u32, r.is_ok()));
        r
    }
    fn insert_vertex(&mut self, f: usize, bary: [f64; 3]) -> Result<Insertion, IntrinsicError> {
        if let Some(note) = (self.note)(self.t, f, bary) {
            self.notes.insert(self.ops.len(), note);
        }
        let mut r = self.t.insert_vertex(f, bary);
        // As the integer crate's refinement: a place App. A cannot recover is a refusal.
        if self.integer
            && let Err(IntrinsicError::Representation { message, .. }) = &r
            && message.contains("could not be")
        {
            r = Err(IntrinsicError::InsertionRefused {
                at: SurfacePoint::Face {
                    face: f as u32,
                    bary,
                },
                reason: ddg_intrinsic::InsertRefusal::NotLocated {
                    detail: message.clone(),
                },
            });
        }
        if let Ok(ins) = &r
            && ins.snapped.is_some()
        {
            self.notes
                .insert(self.ops.len(), format!("snapped {:?}", ins.snapped));
        }
        self.ops.push(Op::Insert {
            face: f as u32,
            bary: bary.map(f64::to_bits),
            result: r.as_ref().map(|i| i.vertex).map_err(kind),
        });
        r
    }
    fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntrinsicError> {
        let r = self.t.split_edge(e, t);
        self.ops.push(Op::Split {
            edge: e as u32,
            t: t.to_bits(),
            result: r.as_ref().map(|i| i.vertex).map_err(kind),
        });
        r
    }
    fn set_constrained(&mut self, e: usize, constrained: bool) -> Result<(), IntrinsicError> {
        self.t.set_constrained(e, constrained)
    }
    fn trace_edge_on_input(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, IntrinsicError> {
        self.t.trace_edge_on_input(e, options)
    }
    fn trace_geodesic(
        &self,
        start: SurfacePoint,
        direction: f64,
        length: f64,
        options: &TraceOptions,
    ) -> Result<Trace, IntrinsicError> {
        self.t.trace_geodesic(start, direction, length, options)
    }
    fn remove_inserted_vertex(
        &mut self,
        v: usize,
    ) -> Result<InsertedVertexRemoval, IntrinsicError> {
        let r = self.t.remove_inserted_vertex(v);
        self.ops.push(Op::Remove(
            v as u32,
            r.as_ref().map(|x| x.flips.len()).map_err(kind),
        ));
        r
    }
}

fn no_note<T>(_: &T, _: usize, _: [f64; 3]) -> Option<String> {
    None
}

fn arc_note(t: &IntegerTriangulation, f: usize, bary: [f64; 3]) -> Option<String> {
    let d = t.locate_in_face(f, bary).ok()?.choice.curve_distance?;
    (d < 1e-6).then(|| format!("{d:.1e} shortest sides from an arc"))
}

/// The first operation at which the two representations' refinements part, with what each did
/// and the notes there, or `None` if they made the same operations.
fn first_divergence(
    start: &IntrinsicTriangulation,
    integer_start: &IntegerTriangulation,
    options: &RefineOptions,
    rules: RefineRules,
) -> Option<String> {
    let mut sig = start.clone();
    let mut int = integer_start.clone();
    let mut a = Log {
        t: &mut sig,
        ops: Vec::new(),
        notes: HashMap::new(),
        note: no_note,
        integer: false,
    };
    let mut b = Log {
        t: &mut int,
        ops: Vec::new(),
        notes: HashMap::new(),
        note: arc_note,
        integer: true,
    };
    let ra = ddg_intrinsic::refine::refine(&mut a, options, rules, &|_, _| false);
    let rb = ddg_intrinsic::refine::refine(&mut b, options, rules, &|_, _| false);
    let n = a.ops.len().min(b.ops.len());
    let i = (0..n).find(|&i| a.ops[i] != b.ops[i]);
    let i = match i {
        Some(i) => i,
        None if a.ops.len() != b.ops.len() => n,
        None => {
            return match (ra, rb) {
                (Ok(_), Ok(_)) => None,
                (x, y) => Some(format!(
                    "same operations, results {:?} / {:?}",
                    x.err(),
                    y.err()
                )),
            };
        }
    };
    // The nearest notes at or before the divergence.
    let near = |notes: &HashMap<usize, String>| {
        notes
            .iter()
            .filter(|(k, _)| **k <= i)
            .max_by_key(|(k, _)| **k)
            .map(|(k, v)| format!("op {k}: {v}"))
    };
    Some(format!(
        "op {i} of {}/{}: signposts {:?} [{:?}], integers {:?} [{:?}]",
        a.ops.len(),
        b.ops.len(),
        a.ops.get(i),
        near(&a.notes),
        b.ops.get(i),
        near(&b.notes),
    ))
}

// ----- the differential test -----

/// The first difference between two integer states, places compared within `place_tol`
/// barycentric, or `None`.
fn differs(a: &IntegerTriangulation, b: &IntegerTriangulation) -> Option<String> {
    if a.complex() != b.complex() {
        return Some("connectivity".into());
    }
    if let Some(e) = (0..a.edge_lengths().len())
        .find(|&e| a.edge_lengths()[e].to_bits() != b.edge_lengths()[e].to_bits())
    {
        return Some(format!("length of edge {e}"));
    }
    let na = a.normal_coordinates();
    let nb = b.normal_coordinates();
    let n_normal = (0..na.len()).filter(|&e| na[e] != nb[e]).count();
    if n_normal > 0 {
        return Some(format!("n of {n_normal} edges"));
    }
    if let Some(h) = (0..a.roundabouts().len()).find(|&h| a.roundabouts()[h] != b.roundabouts()[h])
    {
        return Some(format!("roundabout of halfedge {h}"));
    }
    None
}

fn position(shape: &Shape, p: SurfacePoint) -> [f64; 3] {
    let x = |v: usize| shape.positions[v];
    match p {
        SurfacePoint::Vertex(v) => x(v as usize),
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = shape.mesh.edge_vertices(edge as usize);
            let (pa, pb) = (x(a), x(b));
            [0, 1, 2].map(|c| (1.0 - t) * pa[c] + t * pb[c])
        }
        SurfacePoint::Face { face, bary } => {
            let f = shape.mesh.faces()[face as usize];
            let p = f.map(|v| x(v as usize));
            [0, 1, 2].map(|c| bary[0] * p[0][c] + bary[1] * p[1][c] + bary[2] * p[2][c])
        }
    }
}

/// The counts the two reports must share.
fn counts(r: &RefineReport) -> [usize; 9] {
    [
        r.n_inserted,
        r.n_face_insertions,
        r.n_edge_splits,
        r.n_boundary_splits,
        r.n_flips,
        r.n_removals,
        r.n_removals_refused,
        r.n_trace_failures,
        r.n_narrow_cone_faces,
    ]
}

struct Outcome {
    label: String,
    sig: RefineReport,
    int: Result<IntegerRefineReport, IntegerError>,
    /// The two triangulations: connectivity and lengths (bitwise) equal.
    same_triangulation: bool,
    /// Largest distance of the places of the same vertex, over the mean edge length.
    place_error: f64,
    /// The integer state against `from_intrinsic` of the refined signposts: `None` if equal.
    against_conversion: Option<String>,
    divergence: Option<String>,
}

fn run_one(shape: &Shape, scenario: Scenario, rules: RefineRules, degrees: f64) -> Outcome {
    let label = format!(
        "{} / {scenario:?} / {} / {degrees}°",
        shape.name,
        rules_name(rules)
    );
    let start = state(shape, scenario, 11);
    let integer_start = match scenario {
        Scenario::Initial => IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap(),
        _ => {
            convert(&start)
                .unwrap_or_else(|e| panic!("{label}: conversion {e}"))
                .0
        }
    };
    let opts = options(degrees);
    let mut sig = start.clone();
    let sig_report = match rules {
        RefineRules::Ruppert => sig.refine(&opts),
        RefineRules::Chew => sig.refine_with_rules(&opts, rules),
    }
    .unwrap_or_else(|e| panic!("{label}: signposts {e}"));
    let mut int = integer_start.clone();
    let int_report = int.refine(&integer_options(degrees, rules));
    let same_triangulation = int.complex() == sig.complex()
        && int
            .edge_lengths()
            .iter()
            .zip(sig.edge_lengths())
            .all(|(a, b)| a.to_bits() == b.to_bits());
    let mut place_error = 0.0f64;
    if same_triangulation {
        for v in 0..int.complex().n_vertices() {
            let d = common::dist(
                position(shape, int.vertex_location(v)),
                position(shape, sig.vertex_location(v)),
            );
            place_error = place_error.max(d / shape.mean_edge_length());
        }
    }
    let against_conversion = if same_triangulation {
        match convert(&sig) {
            Ok((c, _)) => differs(&int, &c),
            Err(e) => Some(format!("conversion failed: {e}")),
        }
    } else {
        Some("triangulations differ".into())
    };
    let same_counts = int_report
        .as_ref()
        .is_ok_and(|r| counts(&r.refine) == counts(&sig_report));
    let divergence = if !same_triangulation || !same_counts {
        first_divergence(&start, &integer_start, &opts, rules)
    } else {
        None
    };
    Outcome {
        label,
        sig: sig_report,
        int: int_report,
        same_triangulation,
        place_error,
        against_conversion,
        divergence,
    }
}

/// IC-SPEC I11, I23: the same driver on both representations inserts the same points, makes
/// the same flips and removals, and gives the same triangulation (connectivity, lengths bitwise)
/// wherever the floating-point decisions of the two agree; the integer result certifies and
/// equals `from_intrinsic` of the refined signposts exactly unless a place was decided
/// differently (a signpost place snapped onto an input edge, a point within rounding of an
/// arc).
#[test]
fn refinement_matches_the_signposts() {
    let shapes = shapes();
    let mut jobs = Vec::new();
    for (i, shape) in shapes.iter().enumerate() {
        let big = shape.mesh.n_faces() > 5000;
        let scenarios: &[Scenario] = if big {
            &[Scenario::Initial]
        } else {
            &[
                Scenario::Initial,
                Scenario::RandomFlips,
                Scenario::Insertions,
                Scenario::EdgeSplits,
            ]
        };
        for &scenario in scenarios {
            for rules in [RefineRules::Ruppert, RefineRules::Chew] {
                jobs.push((i, scenario, rules, 25.0));
                if !big && scenario == Scenario::Initial {
                    jobs.push((i, scenario, rules, 30.0));
                }
            }
        }
    }
    let outcomes = par_map(&jobs, |&(i, scenario, rules, degrees)| {
        run_one(&shapes[i], scenario, rules, degrees)
    });
    let mut unexplained = Vec::new();
    let (mut n_same, mut n_exact, mut n_inserted, mut n_removed, mut n_refused) = (0, 0, 0, 0, 0);
    let mut worst_place = 0.0f64;
    let mut n_trace_failures = 0;
    let (mut n_near_arc, mut min_curve, mut n_clamped, mut n_moved, mut n_recover) =
        (0, f64::INFINITY, 0, 0, 0);
    let mut n_snapped = 0;
    println!(
        "{:<52} {:>6} {:>6} {:>5} {:>5} {:>6} {:>6} {:>8}  result",
        "run", "ins", "bnd", "rem", "stuck", "trace", "clamp", "near-arc"
    );
    for o in &outcomes {
        let int = match &o.int {
            Ok(r) => r,
            Err(e) => {
                unexplained.push(format!("{}: integer refinement failed: {e}", o.label));
                continue;
            }
        };
        assert!(int.certificate.is_some(), "{}: not certified", o.label);
        n_inserted += int.refine.n_inserted;
        n_removed += int.refine.n_removals;
        n_refused += int.refine.n_removals_refused;
        n_trace_failures += int.refine.n_trace_failures;
        n_near_arc += int.n_near_arc;
        min_curve = min_curve.min(int.min_curve_distance);
        n_clamped += int.clamped_crossings;
        n_moved += int.n_places_moved;
        n_recover += int.n_recover_failed;
        n_snapped += o.sig.n_snapped;
        let result = if o.same_triangulation && counts(&int.refine) == counts(&o.sig) {
            n_same += 1;
            worst_place = worst_place.max(o.place_error);
            match &o.against_conversion {
                None => {
                    n_exact += 1;
                    "exact".to_string()
                }
                Some(d) => {
                    // Explained by a place decided differently.
                    let why = if o.sig.n_snapped > 0 {
                        format!("{} signpost places snapped", o.sig.n_snapped)
                    } else if int.n_near_arc > 0 {
                        format!("{} splits near an arc", int.n_near_arc)
                    } else {
                        unexplained.push(format!("{}: same triangulation, {d}", o.label));
                        "UNEXPLAINED".to_string()
                    };
                    format!("same triangulation; vs conversion: {d} ({why})")
                }
            }
        } else {
            let d = o.divergence.clone().unwrap_or_default();
            // Floating point decided a place differently: a split point within rounding of a
            // crossing (the integers refuse it, I20; the signposts snap it onto the input
            // edge), a point within rounding of an arc (I6), a snapped signpost place.
            let explained =
                d.contains("within rounding") || d.contains("from an arc") || d.contains("snapped");
            if !explained {
                unexplained.push(format!("{}: {d}", o.label));
            }
            format!("DIVERGE {d}")
        };
        println!(
            "{:<52} {:>6} {:>6} {:>5} {:>5} {:>6} {:>6} {:>8}  {result}",
            o.label,
            int.refine.n_inserted,
            int.refine.n_boundary_splits,
            int.refine.n_removals,
            int.refine.n_removals_refused,
            int.refine.n_trace_failures,
            int.clamped_crossings,
            int.n_near_arc,
        );
    }
    println!(
        "{} runs: {n_same} with the same operations and triangulation ({n_exact} equal to the \
         conversion exactly), places within {worst_place:.1e} mean edge lengths; {n_inserted} \
         insertions, {n_removed} removals, {n_refused} refused, {n_trace_failures} trace \
         failures; integer face splits: {n_near_arc} within 1e-9 of an arc (nearest \
         {min_curve:.1e}), {n_clamped} clamped crossings, {n_moved} places moved inside, \
         {n_recover} refused by App. A; signpost places snapped: {n_snapped}",
        outcomes.len()
    );
    assert!(unexplained.is_empty(), "{unexplained:#?}");
}

// ----- the paper's statements -----

/// Paper Sec. 4.2, Thm. C.1: refinement to 25° reaches the bound everywhere but in the faces
/// the narrow-cone rule leaves (and at corners that cannot be improved); on closed surfaces
/// whose angle sums are all at least 60°, 30° is reached.
#[test]
fn the_angle_bounds_are_reached() {
    let mut shapes = all_shapes();
    shapes.extend(cone_shapes());
    shapes.extend(boundary_shapes());
    let mut jobs = Vec::new();
    for (i, shape) in shapes.iter().enumerate() {
        jobs.push((i, 25.0));
        let closed = (0..shape.mesh.n_vertices()).all(|v| !shape.mesh.is_boundary_vertex(v));
        let sums = ddg_core::operators::vertex_angle_sums(&shape.mesh, &shape.geometry).unwrap();
        let wide = sums.iter().all(|&s| s >= 60f64.to_radians());
        if closed && wide {
            jobs.push((i, 30.0));
        }
    }
    let reports = par_map(&jobs, |&(i, degrees)| {
        let shape = &shapes[i];
        let mut t = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let r = t
            .refine(&integer_options(degrees, RefineRules::Chew))
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        (i, degrees, r)
    });
    println!(
        "{:<28} {:>4} {:>7} {:>7} {:>8} {:>8} {:>6} {:>6}",
        "shape", "bound", "min", "valid", "narrow", "inserted", "below", "status"
    );
    for (i, degrees, r) in &reports {
        let shape = &shapes[*i];
        let rr = &r.refine;
        println!(
            "{:<28} {:>4} {:>7.3} {:>7.3} {:>8} {:>8} {:>6} {:?}",
            shape.name,
            degrees,
            rr.min_angle.to_degrees(),
            rr.min_angle_valid.to_degrees(),
            rr.n_narrow_cone_faces,
            rr.n_inserted,
            rr.below_bound.len(),
            rr.status,
        );
        assert_eq!(rr.status, RefineStatus::Done, "{}", shape.name);
        let bound = degrees.to_radians();
        for b in &rr.below_bound {
            assert!(
                matches!(
                    b.reason,
                    BelowReason::AngleSum { .. } | BelowReason::NeedleVertex { .. }
                ),
                "{}: {b:?}",
                shape.name
            );
        }
        let theorem = *degrees == 30.0;
        if theorem || rr.n_narrow_cone_faces == 0 {
            // Every face is valid; corners below the bound only at vertices whose angle sum
            // is below it.
            assert!(
                rr.min_angle_reachable >= bound * (1.0 - 1e-9),
                "{}: {}",
                shape.name,
                rr.min_angle_reachable.to_degrees()
            );
        }
        if theorem {
            assert!(rr.min_angle_valid >= bound * (1.0 - 1e-9), "{}", shape.name);
            assert!(rr.min_angle >= bound * (1.0 - 1e-9), "{}", shape.name);
        }
    }
}

// ----- termination -----

/// The insertion limit is honoured and reported; the state certifies; refused removals are
/// counted, not fatal.
#[test]
fn insertion_limits_are_honoured_and_refused_removals_counted() {
    let mut shapes = all_shapes();
    shapes.extend(cone_shapes());
    let jobs: Vec<usize> = (0..shapes.len()).collect();
    let results = par_map(&jobs, |&i| {
        let shape = &shapes[i];
        let full = {
            let mut t = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
            t.refine(&integer_options(30.0, RefineRules::Chew)).unwrap()
        };
        let mut lines = Vec::new();
        for limit in [0usize, 1, 7, 40] {
            let mut t = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
            let mut o = integer_options(30.0, RefineRules::Chew);
            o.refine.max_insertions = InsertionLimit::Vertices(limit);
            let r = t.refine(&o).unwrap();
            assert!(r.refine.n_inserted <= limit, "{}", shape.name);
            assert!(r.certificate.is_some());
            if full.refine.n_inserted > limit {
                assert_eq!(
                    r.refine.status,
                    RefineStatus::InsertionLimit,
                    "{}",
                    shape.name
                );
                assert_eq!(r.refine.n_inserted, limit, "{}", shape.name);
            } else if full.refine.n_inserted < limit {
                // (At exactly the limit the driver stops before it looks: `InsertionLimit`.)
                assert_eq!(r.refine.status, RefineStatus::Done, "{}", shape.name);
            }
            lines.push((limit, r.refine.n_inserted, r.refine.status));
        }
        // From a state full of vertices made by edge splits (straight stars, I27), whose
        // removals are often refused.
        let start = state(shape, Scenario::EdgeSplits, 7);
        let mut t = convert(&start).unwrap().0;
        let r = t.refine(&integer_options(25.0, RefineRules::Chew)).unwrap();
        (i, full.refine.n_inserted, lines, r)
    });
    let (mut removed, mut refused) = (0, 0);
    for (i, full, lines, r) in &results {
        println!(
            "{:<28} full {full:>6}  limits {lines:?}  from splits: removed {} refused {}",
            shapes[*i].name, r.refine.n_removals, r.refine.n_removals_refused
        );
        removed += r.refine.n_removals;
        refused += r.refine.n_removals_refused;
    }
    println!("removals from split states: {removed} removed, {refused} refused");
}

/// Bad options and constrained edges are refused before anything changes.
#[test]
fn bad_options_change_nothing() {
    let shape = Shape::generated("grid 6x4", generators::grid_disk(6, 4).unwrap(), true);
    let mut t = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    let before = t.clone();
    for degrees in [-1.0, 36.0, f64::NAN] {
        let err = t
            .refine(&integer_options(degrees, RefineRules::Chew))
            .unwrap_err();
        assert!(
            matches!(
                err,
                IntegerError::Intrinsic(IntrinsicError::InvalidOptions { .. })
            ),
            "{err:?}"
        );
        assert_eq!(t, before);
    }
    Triangulation::set_constrained(&mut t, 3, true).unwrap();
    let constrained = t.clone();
    let err = t.refine(&IntegerRefineOptions::default()).unwrap_err();
    assert!(
        matches!(
            err,
            IntegerError::Intrinsic(IntrinsicError::ConstrainedEdges { .. })
        ),
        "{err:?}"
    );
    assert_eq!(t, constrained);
}

// ----- timings -----

/// Refinement to 25° on spot, camelbad, bunny and meshes with boundary, signposts against
/// integers, both rules; single thread, best of five. Run alone, pinned:
/// `taskset -c N cargo test --release --test refine -- --ignored --nocapture timing`.
#[test]
#[ignore]
fn timing() {
    let mut shapes: Vec<Shape> = Vec::new();
    for path in [
        "models/spot/spot/spot_triangulated.obj",
        "optimal-cones/optimal-cones-data/data/camelbad.obj",
        "optimal-cones/optimal-cones-data/data/bunny.obj",
        "optimal-cones/optimal-cones-data/data/bump.obj",
        "optimal-cones/optimal-cones-data/data/hand.obj",
        "optimal-cones/optimal-cones-data/data/face.obj",
    ] {
        shapes.extend(load_data(path));
    }
    let best = |f: &mut dyn FnMut() -> f64| (0..5).map(|_| f()).fold(f64::INFINITY, f64::min);
    println!(
        "{:<20} {:>7} {:>8} {:>9} {:>6} | {:>10} {:>10} {:>10} {:>10}",
        "mesh", "rules", "|V0|", "inserted", "rem", "signposts", "integer", "+certify", "ratio"
    );
    for shape in &shapes {
        for rules in [RefineRules::Ruppert, RefineRules::Chew] {
            let opts = options(25.0);
            let mut inserted = 0;
            let mut removed = 0;
            let sig = best(&mut || {
                let mut t = shape.triangulation();
                let start = std::time::Instant::now();
                let r = t.refine_with_rules(&opts, rules).unwrap();
                let s = start.elapsed().as_secs_f64();
                inserted = r.n_inserted;
                s
            });
            let int = best(&mut || {
                let mut t = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
                let start = std::time::Instant::now();
                let r = t
                    .refine(&IntegerRefineOptions {
                        certify: false,
                        ..integer_options(25.0, rules)
                    })
                    .unwrap();
                let s = start.elapsed().as_secs_f64();
                removed = r.refine.n_removals;
                assert_eq!(r.refine.n_inserted, inserted, "{}", shape.name);
                s
            });
            let certified = best(&mut || {
                let mut t = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
                let start = std::time::Instant::now();
                t.refine(&integer_options(25.0, rules)).unwrap();
                start.elapsed().as_secs_f64()
            });
            println!(
                "{:<20} {:>7} {:>8} {:>9} {:>6} | {:>8.1}ms {:>8.1}ms {:>8.1}ms {:>10.2}",
                shape.name,
                rules_name(rules),
                shape.mesh.n_vertices(),
                inserted,
                removed,
                1e3 * sig,
                1e3 * int,
                1e3 * certified,
                int / sig,
            );
        }
    }
}
