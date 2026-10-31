//! The operations of WP3 (SPEC §4.2, §4.3, §4.5, §4.6): flips, Delaunay flipping, edge splits
//! and vertex removal, differentially against the signpost representation: the same operation
//! is done in an `IntrinsicTriangulation` and in the `IntegerTriangulation` converted from its
//! state, and afterwards the integer state must **equal**
//! `IntegerTriangulation::from_intrinsic` of the signposts exactly (connectivity, lengths
//! bitwise, every `n`, every roundabout, every place) and certify.
//!
//! Run with `--nocapture` to see the tables.

mod common;

use std::collections::BTreeMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use common::{SCENARIOS, Scenario, Shape, all_shapes, convert, standard_shapes, state};
use ddg_core::generators::SplitMix64;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::formulas::{self, FlipNormalError};
use ddg_integer_coords::{IntegerError, IntegerTriangulation, Parts};
use ddg_intrinsic::{
    Complex, DelaunayOptions, EdgeTrace, EdgeTraceOptions, Flip, FlipRefusal, InputSurface,
    Insertion, IntrinsicError, IntrinsicTriangulation, RemovalRefusal, SurfacePoint, Triangulation,
};

// ----- helpers -----

/// Runs `f` over `specs` on all cores, results in order.
fn par_map<S: Sync, R: Send>(specs: &[S], f: impl Fn(&S) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(48);
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

/// The first difference between the integer state and the conversion of the signposts, or
/// `None` if they are equal (lengths bitwise). `Err` if the conversion fails.
fn difference(
    it: &IntegerTriangulation,
    s: &IntrinsicTriangulation,
) -> Result<Option<String>, String> {
    let (c, _) = convert(s).map_err(|e| format!("conversion failed: {e}"))?;
    Ok(differs(it, &c, true))
}

/// The first difference between two integer states (`with_locations`: also the places).
fn differs(
    a: &IntegerTriangulation,
    b: &IntegerTriangulation,
    with_locations: bool,
) -> Option<String> {
    if a.complex() != b.complex() {
        return Some("connectivity differs".into());
    }
    for (e, (x, y)) in a.edge_lengths().iter().zip(b.edge_lengths()).enumerate() {
        if x.to_bits() != y.to_bits() {
            return Some(format!("length of edge {e}: {x:e} vs {y:e}"));
        }
    }
    for (e, (x, y)) in a
        .normal_coordinates()
        .iter()
        .zip(b.normal_coordinates())
        .enumerate()
    {
        if x != y {
            return Some(format!("n of edge {e}: {x} vs {y}"));
        }
    }
    for (h, (x, y)) in a.roundabouts().iter().zip(b.roundabouts()).enumerate() {
        if x != y {
            return Some(format!("roundabout of halfedge {h}: {x} vs {y}"));
        }
    }
    if with_locations {
        for (v, (x, y)) in a
            .vertex_locations()
            .iter()
            .zip(b.vertex_locations())
            .enumerate()
        {
            if x != y {
                return Some(format!("place of vertex {v}: {x:?} vs {y:?}"));
            }
        }
    }
    None
}

/// Number of edges from a vertex to itself and of faces that use a vertex twice (I21).
fn delta_features(c: &Complex) -> (usize, usize) {
    let loops = (0..c.n_edges())
        .filter(|&e| {
            let [a, b] = c.edge_vertices(e);
            a == b
        })
        .count();
    let faces = c
        .faces()
        .iter()
        .filter(|f| f[0] == f[1] || f[1] == f[2] || f[2] == f[0])
        .count();
    (loops, faces)
}

/// A triangulation whose flips are logged (for the identical-sequence check of SPEC §4.3).
struct Recorder<'a, T> {
    t: &'a mut T,
    flips: Vec<u32>,
}

impl<T: Triangulation> Triangulation for Recorder<'_, T> {
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
        if r.is_ok() {
            self.flips.push(e as u32);
        }
        r
    }
    fn insert_vertex(&mut self, f: usize, bary: [f64; 3]) -> Result<Insertion, IntrinsicError> {
        self.t.insert_vertex(f, bary)
    }
    fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntrinsicError> {
        self.t.split_edge(e, t)
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
}

fn shape_scale(shape: &Shape) -> usize {
    shape.mesh.n_edges()
}

/// The jobs of the differential tests: every shape, scenario and seed.
fn specs(shapes: &[Shape]) -> Vec<(usize, Scenario, u64)> {
    let mut out = Vec::new();
    for i in 0..shapes.len() {
        for scenario in SCENARIOS {
            for seed in [5u64, 77] {
                out.push((i, scenario, seed));
            }
        }
    }
    out
}

// ----- flips and Delaunay flipping -----

#[derive(Default, Debug)]
struct FlipStats {
    states: usize,
    batches: usize,
    attempts: usize,
    flips: usize,
    refusals_agree: usize,
    delaunay_flips: usize,
    loops_max: usize,
    repeated_faces_max: usize,
    mismatches: Vec<String>,
    conversion_failures: Vec<String>,
}

#[test]
fn flips_and_delaunay_match_the_signposts_exactly() {
    let shapes = all_shapes();
    let specs = specs(&shapes);
    let results = par_map(&specs, |&(i, scenario, seed)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?} / {seed}", shape.name);
        let mut st = FlipStats::default();
        let mut s = state(shape, scenario, seed);
        let mut it = match convert(&s) {
            Ok((it, _)) => it,
            Err(e) => {
                st.conversion_failures.push(format!("{tag}: initial: {e}"));
                return st;
            }
        };
        st.states = 1;
        let mut rng = SplitMix64::new(seed ^ 0x9e37_79b9);
        let big = shape_scale(shape) > 5000;
        let (batches, per_batch) = if big {
            (2, 1500)
        } else {
            (4, (3 * s.n_edges()).min(3000))
        };
        for b in 0..batches {
            for _ in 0..per_batch {
                let e = (rng.next_u64() % s.n_edges() as u64) as usize;
                st.attempts += 1;
                let rs = s.flip_edge(e);
                let ri = it.flip_edge(e);
                match (&rs, &ri) {
                    (Ok(fs), Ok(fi)) => {
                        st.flips += 1;
                        if fs != fi {
                            st.mismatches
                                .push(format!("{tag}: flip {e}: {fs:?} vs {fi:?}"));
                        }
                    }
                    (Err(es), Err(IntegerError::Intrinsic(ei))) if es == ei => {
                        st.refusals_agree += 1
                    }
                    _ => st
                        .mismatches
                        .push(format!("{tag}: flip {e}: {rs:?} vs {ri:?}")),
                }
            }
            st.batches += 1;
            let (lp, rf) = delta_features(it.complex());
            st.loops_max = st.loops_max.max(lp);
            st.repeated_faces_max = st.repeated_faces_max.max(rf);
            match difference(&it, &s) {
                Ok(None) => {}
                Ok(Some(d)) => st.mismatches.push(format!("{tag}: batch {b}: {d}")),
                Err(e) => st
                    .conversion_failures
                    .push(format!("{tag}: batch {b}: {e}")),
            }
            if let Err(e) = it.certify() {
                st.mismatches
                    .push(format!("{tag}: batch {b}: certify: {e}"));
            }
        }
        // Delaunay flipping: identical sequences.
        let opts = DelaunayOptions::default();
        let mut rs = Recorder {
            t: &mut s,
            flips: Vec::new(),
        };
        let report_s = rs.flip_to_delaunay(&opts).unwrap();
        let flips_s = rs.flips;
        let mut ri = Recorder {
            t: &mut it,
            flips: Vec::new(),
        };
        let report_i = ri.flip_to_delaunay(&opts).unwrap();
        let flips_i = ri.flips;
        st.delaunay_flips += flips_i.len();
        if flips_s != flips_i || report_s != report_i {
            st.mismatches.push(format!(
                "{tag}: Delaunay: {} vs {} flips, reports {report_s:?} vs {report_i:?}",
                flips_s.len(),
                flips_i.len()
            ));
        }
        match difference(&it, &s) {
            Ok(None) => {}
            Ok(Some(d)) => st.mismatches.push(format!("{tag}: after Delaunay: {d}")),
            Err(e) => st
                .conversion_failures
                .push(format!("{tag}: after Delaunay: {e}")),
        }
        if let Err(e) = it.certify() {
            st.mismatches
                .push(format!("{tag}: after Delaunay: certify: {e}"));
        }
        st
    });
    let mut total = FlipStats::default();
    for r in results {
        total.states += r.states;
        total.batches += r.batches;
        total.attempts += r.attempts;
        total.flips += r.flips;
        total.refusals_agree += r.refusals_agree;
        total.delaunay_flips += r.delaunay_flips;
        total.loops_max = total.loops_max.max(r.loops_max);
        total.repeated_faces_max = total.repeated_faces_max.max(r.repeated_faces_max);
        total.mismatches.extend(r.mismatches);
        total.conversion_failures.extend(r.conversion_failures);
    }
    eprintln!(
        "flip differential: {} states, {} batches, {} attempts, {} flips, {} refusals (all \
         identical), {} Delaunay flips; up to {} self-edges and {} faces with a repeated vertex \
         in one state; {} mismatches, {} conversion failures",
        total.states,
        total.batches,
        total.attempts,
        total.flips,
        total.refusals_agree,
        total.delaunay_flips,
        total.loops_max,
        total.repeated_faces_max,
        total.mismatches.len(),
        total.conversion_failures.len()
    );
    for m in total.conversion_failures.iter().chain(&total.mismatches) {
        eprintln!("  {m}");
    }
    assert!(total.mismatches.is_empty(), "{:#?}", total.mismatches);
    assert!(
        total.conversion_failures.is_empty(),
        "{:#?}",
        total.conversion_failures
    );
}

/// Per edge, `n` and the sorted list of (tail, roundabout) of its halfedges: what a flip and
/// its inverse must restore (the halfedge ids of the two faces change).
fn integers_by_edge(t: &IntegerTriangulation) -> Vec<(i64, Vec<(usize, u32)>)> {
    let c = t.complex();
    (0..c.n_edges())
        .map(|e| {
            let h = c.edge_halfedge(e);
            let mut r: Vec<(usize, u32)> = std::iter::once(h)
                .chain(c.twin(h))
                .map(|g| (c.from_vertex(g), t.roundabouts()[g]))
                .collect();
            r.sort_unstable();
            (t.normal_coordinate(e), r)
        })
        .collect()
}

#[test]
fn flipping_back_restores_the_integers_exactly() {
    let shapes = standard_shapes();
    let mut n_pairs = 0usize;
    for shape in &shapes {
        for scenario in [Scenario::RandomFlips, Scenario::Insertions, Scenario::Mixed] {
            let s = state(shape, scenario, 3);
            let (it, _) = convert(&s).unwrap();
            let before = integers_by_edge(&it);
            for e in 0..it.complex().n_edges() {
                let mut t = it.clone();
                if t.flip_edge(e).is_err() {
                    continue;
                }
                let flipped_n = t.normal_coordinate(e);
                // The involution on the integers alone (ic-flip-involution).
                t.flip_edge(e).unwrap_or_else(|err| {
                    panic!("{}: edge {e} cannot be flipped back: {err}", shape.name)
                });
                assert_eq!(
                    integers_by_edge(&t),
                    before,
                    "{} / {scenario:?}: edge {e} (n after the flip {flipped_n})",
                    shape.name
                );
                let l0 = it.edge_lengths()[e];
                assert!((t.edge_lengths()[e] - l0).abs() <= 1e-12 * l0);
                n_pairs += 1;
            }
        }
    }
    eprintln!("flip and flip back: {n_pairs} edges, integers restored exactly");
    assert!(n_pairs > 1000);
}

// ----- edge splits -----

#[derive(Default, Debug)]
struct SplitStats {
    along_known: usize,
    along_traced: usize,
    along_traced_max_distance: f64,
    crossed_given: usize,
    crossed_own: usize,
    crossed_own_max_distance: f64,
    crossed_own_refused: Vec<String>,
    crossed_skipped: usize,
    refusals_agree: usize,
    mismatches: Vec<String>,
}

/// The distance of two places on the same input edge (relative to its length), or of two
/// places in the same input face (largest barycentric difference); infinite otherwise.
fn location_distance(a: SurfacePoint, b: SurfacePoint) -> f64 {
    match (a, b) {
        (SurfacePoint::Edge { edge: x, t: s }, SurfacePoint::Edge { edge: y, t }) if x == y => {
            (s - t).abs()
        }
        (SurfacePoint::Face { face: x, bary: p }, SurfacePoint::Face { face: y, bary: q })
            if x == y =>
        {
            (0..3).map(|k| (p[k] - q[k]).abs()).fold(0.0, f64::max)
        }
        _ => f64::INFINITY,
    }
}

/// The insertions agree, except for `snapped`, which only the signposts set.
fn same_insertion(a: &Insertion, b: &Insertion) -> bool {
    a.vertex == b.vertex
        && a.location == b.location
        && a.new_edges == b.new_edges
        && a.new_faces == b.new_faces
        && a.halves == b.halves
}

#[test]
fn edge_splits_match_the_signposts_exactly() {
    let shapes = all_shapes();
    let specs = specs(&shapes);
    let results = par_map(&specs, |&(i, scenario, seed)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?} / {seed}", shape.name);
        let mut st = SplitStats::default();
        let mut s = state(shape, scenario, seed);
        let Ok((mut it, _)) = convert(&s) else {
            st.mismatches
                .push(format!("{tag}: initial conversion failed"));
            return st;
        };
        let mut rng = SplitMix64::new(seed ^ 0x5eed);
        let attempts = if shape_scale(shape) > 5000 { 30 } else { 40 };
        let opts = EdgeTraceOptions::default();
        for a in 0..attempts {
            let e = (rng.next_u64() % s.n_edges() as u64) as usize;
            let t = rng.uniform(0.15, 0.85);
            let n = it.normal_coordinate(e);
            let before_state = it.clone();
            if n == -1 && s.edge_on_input(e).is_none() {
                // Along an input edge that the signposts do not know: they place the vertex by
                // tracing. Compare everything but the place, then put the state back.
                let (mut s2, mut it2) = (s.clone(), it.clone());
                match (s2.split_edge(e, t), it2.split_edge(e, t)) {
                    (Ok(xs), Ok(xi)) => {
                        let (c, _) = convert(&s2).unwrap();
                        if let Some(d) = differs(&it2, &c, false) {
                            st.mismatches.push(format!("{tag}: traced split {e}: {d}"));
                        }
                        let dist = location_distance(xs.location, xi.location);
                        st.along_traced_max_distance = st.along_traced_max_distance.max(dist);
                        st.along_traced += 1;
                    }
                    (Err(_), Err(_)) => st.refusals_agree += 1,
                    (rs, ri) => st
                        .mismatches
                        .push(format!("{tag}: traced split {e}: {rs:?} vs {ri:?}")),
                }
                continue;
            }
            if n >= 0 {
                // With the geometric decisions taken from the signposts.
                let trace = s.trace_edge(e, &opts).unwrap();
                let at = t * s.edge_lengths()[e];
                let tol = 1e-9 * s.edge_lengths()[e];
                if trace
                    .crossings
                    .iter()
                    .any(|x| (x.distance - at).abs() < tol)
                {
                    st.crossed_skipped += 1;
                    continue;
                }
                let before = trace.crossings.iter().filter(|x| x.distance < at).count() as i64;
                // The integer representation's own split (I20 with the crossings of
                // ExtractIntrinsicEdge): everything but the place exact, on copies.
                let (mut s2, mut it2) = (s.clone(), it.clone());
                match (s2.split_edge(e, t), it2.split_edge(e, t)) {
                    (Ok(xs), Ok(xi)) => {
                        st.crossed_own += 1;
                        let (c, _) = convert(&s2).unwrap();
                        if let Some(d) = differs(&it2, &c, false) {
                            st.mismatches.push(format!("{tag}: own split {e}: {d}"));
                        }
                        if let Err(err) = it2.certify() {
                            st.mismatches
                                .push(format!("{tag}: own split {e}: certify: {err}"));
                        }
                        let dist = location_distance(xs.location, xi.location);
                        st.crossed_own_max_distance = st.crossed_own_max_distance.max(dist);
                    }
                    (Err(_), Err(IntegerError::Intrinsic(_))) => {}
                    (Ok(_), Err(err)) => st.crossed_own_refused.push(format!("{tag}: {e}: {err}")),
                    (rs, ri) => st
                        .mismatches
                        .push(format!("{tag}: own split {e}: {rs:?} vs {ri:?}")),
                }
                // With the geometric decisions taken from the signposts: exact.
                match s.split_edge(e, t) {
                    Ok(xs) => match it.split_edge_given(e, t, before, xs.location) {
                        Ok(xi) => {
                            st.crossed_given += 1;
                            if !same_insertion(&xs, &xi) {
                                st.mismatches
                                    .push(format!("{tag}: split {e}: {xs:?} vs {xi:?}"));
                            }
                        }
                        Err(err) => st.mismatches.push(format!("{tag}: split {e}: {err}")),
                    },
                    Err(_) => {
                        st.refusals_agree += 1;
                        continue;
                    }
                }
            } else {
                match (s.split_edge(e, t), it.split_edge(e, t)) {
                    (Ok(xs), Ok(xi)) => {
                        st.along_known += 1;
                        if !same_insertion(&xs, &xi) {
                            st.mismatches
                                .push(format!("{tag}: split {e}: {xs:?} vs {xi:?}"));
                        }
                    }
                    (Err(es), Err(IntegerError::Intrinsic(ei))) if es == ei => {
                        st.refusals_agree += 1;
                        if it != before_state {
                            st.mismatches
                                .push(format!("{tag}: refused split changed the state"));
                        }
                        continue;
                    }
                    (rs, ri) => {
                        st.mismatches
                            .push(format!("{tag}: split {e}: {rs:?} vs {ri:?}"));
                        return st;
                    }
                }
            }
            match difference(&it, &s) {
                Ok(None) => {}
                Ok(Some(d)) => {
                    st.mismatches
                        .push(format!("{tag}: after split {a} (edge {e}, n {n}): {d}"));
                    return st;
                }
                Err(err) => {
                    st.mismatches.push(format!("{tag}: after split {a}: {err}"));
                    return st;
                }
            }
            if let Err(err) = it.certify() {
                st.mismatches
                    .push(format!("{tag}: after split {a}: certify: {err}"));
                return st;
            }
        }
        st
    });
    let mut total = SplitStats::default();
    for r in results {
        total.along_known += r.along_known;
        total.along_traced += r.along_traced;
        total.along_traced_max_distance = total
            .along_traced_max_distance
            .max(r.along_traced_max_distance);
        total.crossed_given += r.crossed_given;
        total.crossed_own += r.crossed_own;
        total.crossed_own_max_distance = total
            .crossed_own_max_distance
            .max(r.crossed_own_max_distance);
        total.crossed_own_refused.extend(r.crossed_own_refused);
        total.crossed_skipped += r.crossed_skipped;
        total.refusals_agree += r.refusals_agree;
        total.mismatches.extend(r.mismatches);
    }
    eprintln!(
        "split differential: {} splits along input edges known to the signposts (exact), {} \
         along input edges the signposts trace (all but the place exact; places within {:e}), \
         {} splits of crossed edges with the signposts' crossings and place (exact), {} with \
         this crate's own crossings (all but the place exact; places within {:e}; {} refused \
         where the signposts split), {} skipped (a crossing within 1e-9 of the point), {} \
         refusals identical; {} mismatches",
        total.along_known,
        total.along_traced,
        total.along_traced_max_distance,
        total.crossed_given,
        total.crossed_own,
        total.crossed_own_max_distance,
        total.crossed_own_refused.len(),
        total.crossed_skipped,
        total.refusals_agree,
        total.mismatches.len()
    );
    for m in &total.mismatches {
        eprintln!("  {m}");
    }
    for m in &total.crossed_own_refused {
        eprintln!("  refused: {m}");
    }
    assert!(total.mismatches.is_empty());
    assert!(total.along_known > 100 && total.crossed_given > 100 && total.crossed_own > 100);
    assert!(total.along_traced_max_distance < 1e-9);
    assert!(total.crossed_own_max_distance < 1e-9);
    assert!(total.crossed_own_refused.is_empty());
}

#[test]
fn split_side_is_eq6_along_input_edges() {
    for a in -1..=6 {
        for b in -1..=6 {
            for (s, n) in [(0, [-1, a, b]), (1, [b, -1, a]), (2, [a, b, -1])] {
                let x = formulas::split_side(n, s, 0).unwrap();
                // n[s + 2] = n_ca (= n_ki), n[s + 1] = n_bc (= n_jk).
                assert_eq!(
                    x.to_opposite,
                    formulas::split_edge_normal(n[(s + 2) % 3], n[(s + 1) % 3])
                );
                assert_eq!(x.emanating_before, 0);
            }
        }
    }
}

// ----- removal -----

#[derive(Default, Debug)]
struct RemovalStats {
    removed: usize,
    flips: usize,
    max_flips: usize,
    by_degree: BTreeMap<usize, usize>,
    stuck: usize,
    stuck_by_degree: BTreeMap<usize, usize>,
    /// Stuck at a star with two spokes on one straight line (angle π between them).
    stuck_straight: usize,
    refused: BTreeMap<String, usize>,
    split_vertices: usize,
    split_removed: usize,
    split_stuck: usize,
    /// Of `split_vertices`, those on edges crossing input edges.
    split_crossed: usize,
    mismatches: Vec<String>,
}

impl RemovalStats {
    fn add(&mut self, o: RemovalStats) {
        self.removed += o.removed;
        self.flips += o.flips;
        self.max_flips = self.max_flips.max(o.max_flips);
        for (k, v) in o.by_degree {
            *self.by_degree.entry(k).or_default() += v;
        }
        self.stuck += o.stuck;
        self.stuck_straight += o.stuck_straight;
        for (k, v) in o.stuck_by_degree {
            *self.stuck_by_degree.entry(k).or_default() += v;
        }
        for (k, v) in o.refused {
            *self.refused.entry(k).or_default() += v;
        }
        self.split_vertices += o.split_vertices;
        self.split_removed += o.split_removed;
        self.split_stuck += o.split_stuck;
        self.split_crossed += o.split_crossed;
        self.mismatches.extend(o.mismatches);
    }
}

/// Whether two spokes of `v` lie on one straight line: the corner angles between them sum to
/// `π` (within `1e-9`), as at a vertex made by an edge split (I27).
fn straight_star(t: &IntegerTriangulation, v: usize) -> bool {
    let angles: Vec<f64> = t
        .complex()
        .outgoing_halfedges(v)
        .map(|h| Triangulation::corner_angle(t, h))
        .collect();
    let n = angles.len();
    (0..n).any(|start| {
        let mut sum = 0.0;
        (0..n - 1).any(|k| {
            sum += angles[(start + k) % n];
            (sum - std::f64::consts::PI).abs() < 1e-9
        })
    })
}

/// Removes inserted vertices from the integer state one by one (random order), doing the same
/// flips and the same removal in the signposts, and compares after each.
fn removals(
    tag: &str,
    s: &mut IntrinsicTriangulation,
    it: &mut IntegerTriangulation,
    rng: &mut SplitMix64,
    st: &mut RemovalStats,
    max: usize,
) {
    let n_input = it.n_input_vertices();
    let mut stuck: Vec<usize> = Vec::new();
    for _ in 0..max {
        let candidates: Vec<usize> = (n_input..it.complex().n_vertices())
            .filter(|&v| !it.complex().is_boundary_vertex(v) && !stuck.contains(&v))
            .collect();
        if candidates.is_empty() {
            break;
        }
        let v = candidates[(rng.next_u64() % candidates.len() as u64) as usize];
        let before = it.clone();
        match it.remove_vertex(v) {
            Ok(r) => {
                st.removed += 1;
                st.flips += r.flips.len();
                st.max_flips = st.max_flips.max(r.flips.len());
                *st.by_degree.entry(r.degree).or_default() += 1;
                for &e in &r.flips {
                    if let Err(err) = s.flip_edge(e as usize) {
                        st.mismatches
                            .push(format!("{tag}: signposts refuse flip {e}: {err}"));
                        return;
                    }
                }
                match s.remove_vertex(v) {
                    Ok(x) if x == r.removal => {}
                    other => {
                        st.mismatches
                            .push(format!("{tag}: signposts' removal of {v}: {other:?}"));
                        return;
                    }
                }
                match difference(it, s) {
                    Ok(None) => {}
                    Ok(Some(d)) => {
                        st.mismatches
                            .push(format!("{tag}: after removing {v}: {d}"));
                        return;
                    }
                    Err(e) => {
                        st.mismatches
                            .push(format!("{tag}: after removing {v}: {e}"));
                        return;
                    }
                }
                if let Err(e) = it.certify() {
                    st.mismatches
                        .push(format!("{tag}: after removing {v}: certify: {e}"));
                    return;
                }
                stuck.clear();
            }
            Err(err) => {
                if *it != before {
                    st.mismatches
                        .push(format!("{tag}: failed removal of {v} changed the state"));
                    return;
                }
                match err {
                    IntegerError::RemovalStuck { degree, .. } => {
                        st.stuck += 1;
                        *st.stuck_by_degree.entry(degree).or_default() += 1;
                        if straight_star(it, v) {
                            st.stuck_straight += 1;
                        }
                    }
                    other => {
                        let key = match &other {
                            IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                                reason,
                                ..
                            }) => format!("{reason:?}")
                                .split([' ', '{', '('])
                                .next()
                                .unwrap()
                                .to_string(),
                            x => format!("{x}"),
                        };
                        *st.refused.entry(key).or_default() += 1;
                    }
                }
                stuck.push(v);
            }
        }
    }
}

#[test]
fn removals_match_the_signposts_exactly() {
    let shapes = all_shapes();
    let mut specs = Vec::new();
    for i in 0..shapes.len() {
        for scenario in [
            Scenario::Insertions,
            Scenario::EdgeSplits,
            Scenario::SplitInputEdgesThenFlips,
            Scenario::Mixed,
        ] {
            for seed in [5u64, 77] {
                specs.push((i, scenario, seed));
            }
        }
    }
    let results = par_map(&specs, |&(i, scenario, seed)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?} / {seed}", shape.name);
        let mut st = RemovalStats::default();
        let mut s = state(shape, scenario, seed);
        let Ok((mut it, _)) = convert(&s) else {
            st.mismatches.push(format!("{tag}: conversion failed"));
            return st;
        };
        let mut rng = SplitMix64::new(seed ^ 0xdead);
        let max = if shape_scale(shape) > 5000 { 25 } else { 60 };
        removals(&tag, &mut s, &mut it, &mut rng, &mut st, max);
        st
    });
    let mut total = RemovalStats::default();
    let mut per_scenario: BTreeMap<String, [usize; 3]> = BTreeMap::new();
    for (r, &(_, scenario, _)) in results.into_iter().zip(&specs) {
        let x = per_scenario.entry(format!("{scenario:?}")).or_default();
        x[0] += r.removed;
        x[1] += r.stuck;
        x[2] += r.stuck_straight;
        total.add(r);
    }
    eprintln!("removal per scenario [removed, stuck, stuck at a straight star]: {per_scenario:?}");
    eprintln!(
        "removal differential: {} removed ({} flips, {:.2} per removal, at most {}; by degree \
         before: {:?}), {} stuck (state restored; degree when stuck: {:?}; {} of them at a star with two \
         spokes on a straight line), refused: {:?}; {} mismatches",
        total.removed,
        total.flips,
        total.flips as f64 / total.removed.max(1) as f64,
        total.max_flips,
        total.by_degree,
        total.stuck,
        total.stuck_by_degree,
        total.stuck_straight,
        total.refused,
        total.mismatches.len()
    );
    for m in &total.mismatches {
        eprintln!("  {m}");
    }
    assert!(total.mismatches.is_empty());
    assert!(total.removed > 500);
}

/// I27: vertices made by splitting edges along input edges (degree four, two angles of π),
/// removed right after the split and after Delaunay flipping.
#[test]
fn removal_of_split_vertices_i27() {
    let shapes = all_shapes();
    let mut specs = Vec::new();
    for i in 0..shapes.len() {
        for scenario in [Scenario::Initial, Scenario::Delaunay, Scenario::RandomFlips] {
            for seed in [5u64, 77] {
                specs.push((i, scenario, seed));
            }
        }
    }
    let results = par_map(&specs, |&(i, scenario, seed)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?} / {seed}", shape.name);
        let mut st = RemovalStats::default();
        let s0 = state(shape, scenario, seed);
        let Ok((it0, _)) = convert(&s0) else {
            st.mismatches.push(format!("{tag}: conversion failed"));
            return st;
        };
        let mut rng = SplitMix64::new(seed ^ 0x127);
        let (mut s, mut it) = (s0.clone(), it0.clone());
        let attempts = if shape_scale(shape) > 5000 { 20 } else { 30 };
        // Right after the split.
        for _ in 0..attempts {
            let e = (rng.next_u64() % it.complex().n_edges() as u64) as usize;
            let n = it.normal_coordinate(e);
            if it.complex().is_boundary_edge(e) || (n == -1 && s.edge_on_input(e).is_none()) {
                continue;
            }
            let t = rng.uniform(0.2, 0.8);
            let xs = if n == -1 {
                let (Ok(xs), Ok(_)) = (s.split_edge(e, t), it.split_edge(e, t)) else {
                    continue;
                };
                xs
            } else {
                // A crossed edge, with the signposts' decisions so that the states stay equal.
                let trace = s.trace_edge(e, &EdgeTraceOptions::default()).unwrap();
                let at = t * s.edge_lengths()[e];
                let before = trace.crossings.iter().filter(|x| x.distance < at).count() as i64;
                let Ok(xs) = s.split_edge(e, t) else { continue };
                it.split_edge_given(e, t, before, xs.location).unwrap();
                st.split_crossed += 1;
                xs
            };
            st.split_vertices += 1;
            let v = xs.vertex as usize;
            let before = it.clone();
            match it.remove_vertex(v) {
                Ok(r) => {
                    st.split_removed += 1;
                    for &f in &r.flips {
                        s.flip_edge(f as usize).unwrap();
                    }
                    s.remove_vertex(v).unwrap();
                    if let Ok(Some(d)) | Err(d) = difference(&it, &s) {
                        st.mismatches.push(format!("{tag}: split vertex {v}: {d}"));
                        return st;
                    }
                    it.certify().unwrap();
                }
                Err(IntegerError::RemovalStuck { .. }) => {
                    st.split_stuck += 1;
                    assert!(it == before, "{tag}: stuck removal changed the state");
                }
                Err(e) => st.mismatches.push(format!("{tag}: split vertex {v}: {e}")),
            }
        }
        // After Delaunay flipping: split some more, flip, remove all inserted vertices.
        for _ in 0..attempts {
            let e = (rng.next_u64() % it.complex().n_edges() as u64) as usize;
            if it.normal_coordinate(e) == -1 && s.edge_on_input(e).is_some() {
                let t = rng.uniform(0.2, 0.8);
                let _ = (s.split_edge(e, t), it.split_edge(e, t));
            }
        }
        let opts = DelaunayOptions::default();
        s.flip_to_delaunay(&opts).unwrap();
        Triangulation::flip_to_delaunay(&mut it, &opts).unwrap();
        removals(&tag, &mut s, &mut it, &mut rng, &mut st, 200);
        st
    });
    let mut total = RemovalStats::default();
    for r in results {
        total.add(r);
    }
    eprintln!(
        "I27: {} vertices split from interior edges ({} of them crossing input edges) and \
         removed at once: {} removed, {} stuck \
         (state restored); after Delaunay flipping: {} removed ({} flips), {} stuck (degree: \
         {:?}; {} at a straight star), refused {:?}; {} mismatches",
        total.split_vertices,
        total.split_crossed,
        total.split_removed,
        total.split_stuck,
        total.removed,
        total.flips,
        total.stuck,
        total.stuck_by_degree,
        total.stuck_straight,
        total.refused,
        total.mismatches.len()
    );
    for m in &total.mismatches {
        eprintln!("  {m}");
    }
    assert!(total.mismatches.is_empty());
    assert!(total.split_vertices > 100);
}

#[test]
fn refused_removals_change_nothing() {
    let shapes = standard_shapes();
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for shape in &shapes {
        let s = state(shape, Scenario::Mixed, 9);
        let (mut it, _) = convert(&s).unwrap();
        // An input vertex.
        let before = it.clone();
        match it.remove_vertex(0) {
            Err(IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                reason: RemovalRefusal::InputVertex,
                ..
            })) => *counts.entry("input vertex").or_default() += 1,
            other => panic!("{}: {other:?}", shape.name),
        }
        assert!(it == before);
        // Out of range.
        assert!(matches!(
            it.remove_vertex(it.complex().n_vertices()),
            Err(IntegerError::IndexOutOfRange { .. })
        ));
        for v in it.n_input_vertices()..it.complex().n_vertices() {
            if it.complex().is_boundary_vertex(v) {
                match it.remove_vertex(v) {
                    Err(IntegerError::Unsupported { .. }) => {
                        *counts.entry("boundary (unsupported)").or_default() += 1
                    }
                    other => panic!("{}: {other:?}", shape.name),
                }
                assert!(it == before);
            } else {
                // A constrained edge at the vertex.
                let e = it.complex().edge(it.complex().outgoing_halfedge(v));
                let mut t = it.clone();
                Triangulation::set_constrained(&mut t, e, true).unwrap();
                let before_c = t.clone();
                match t.remove_vertex(v) {
                    Err(IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                        reason: RemovalRefusal::ConstrainedEdge { .. },
                        ..
                    })) => *counts.entry("constrained edge").or_default() += 1,
                    other => panic!("{}: {other:?}", shape.name),
                }
                assert!(t == before_c);
            }
        }
    }
    eprintln!("refused removals (state unchanged): {counts:?}");
    assert!(counts.len() == 3);
}

// ----- the shared interface -----

#[test]
fn the_interface_delegates_and_refuses_face_insertion() {
    let shape = &standard_shapes()[0];
    let mut it = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    let before = it.clone();
    assert!(matches!(
        Triangulation::insert_vertex(&mut it, 0, [0.2, 0.3, 0.5]),
        Err(IntrinsicError::Unsupported {
            representation: "integer coordinates",
            ..
        })
    ));
    // WP2: ExtractIntrinsicEdge; edge 0 of the input itself runs along input edge 0.
    assert_eq!(
        it.trace_edge_on_input(0, &EdgeTraceOptions::default())
            .unwrap()
            .along_input_edge,
        Some(0)
    );
    assert!(matches!(
        Triangulation::flip_edge(&mut it, usize::MAX),
        Err(IntrinsicError::IndexOutOfRange { .. })
    ));
    assert!(it == before);
    // Lengths, angle sums and places as the signposts.
    let s = shape.triangulation();
    for v in 0..s.n_vertices() {
        assert_eq!(
            Triangulation::vertex_angle_sum(&it, v),
            Triangulation::vertex_angle_sum(&s, v)
        );
        assert_eq!(
            Triangulation::vertex_location(&it, v),
            Triangulation::vertex_location(&s, v)
        );
    }
}

// ----- overflow -----

/// A flat unit square split along its diagonal 0–2 (edge `diag`), with hand-set integers.
fn square(normal_of: impl Fn(usize, bool) -> i64) -> (IntegerTriangulation, usize) {
    let mesh = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
    let pos = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let geometry = Geometry::from_positions(&mesh, &pos).unwrap();
    let it = IntegerTriangulation::new(&mesh, &geometry).unwrap();
    let mut parts: Parts = it.into_parts();
    let diag = (0..parts.complex.n_edges())
        .find(|&e| !parts.complex.is_boundary_edge(e))
        .unwrap();
    for e in 0..parts.complex.n_edges() {
        parts.normal[e] = normal_of(e, e == diag);
    }
    (IntegerTriangulation::from_parts(parts).unwrap(), diag)
}

#[test]
fn overflow_is_reported_never_wrapped() {
    for m in [i64::MAX, i64::MAX - 1, i64::MAX / 2 + 1, i64::MAX / 4 + 7] {
        // Every side of both triangles crossed m times, the diagonal not at all: m arcs cut
        // each corner at the ends of the new diagonal, which would cross 2m of them.
        let (mut it, diag) = square(|_, d| if d { 0 } else { m });
        let before = it.clone();
        let expect_overflow = 2 * (m as i128) > i64::MAX as i128;
        match it.flip_edge(diag) {
            Err(IntegerError::Overflow { edge }) => {
                assert!(expect_overflow);
                assert_eq!(edge as usize, diag);
            }
            Ok(_) => assert!(!expect_overflow),
            other => panic!("{m}: {other:?}"),
        }
        if expect_overflow {
            assert!(it == before);
            assert!(matches!(
                Triangulation::flip_edge(&mut it, diag),
                Err(IntrinsicError::Representation {
                    representation: "integer coordinates",
                    ..
                })
            ));
        }
        assert_eq!(
            formulas::flip_normal(0, m, m, m, m),
            if expect_overflow {
                Err(FlipNormalError::Overflow)
            } else {
                Ok(2 * m)
            }
        );
    }
    // Eq. 4 without reduction near the top of the range.
    assert_eq!(formulas::roundabout_update(i64::MAX, 1, 0, 0), None);
    assert_eq!(
        formulas::roundabout_update(i64::MAX, i64::MAX, 1, 7),
        Some(((i64::MAX as i128 * 2 + 1) % 7) as i64)
    );
    // Eq. 6 and the crossed split at the top of the range.
    assert_eq!(formulas::split_edge_normal(i64::MAX, -1), i64::MAX);
    // Corner 1 is cut by i64::MAX arcs; with the point beyond all of them every arc crosses
    // the new edge, before all of them none does.
    let side = |before| formulas::split_side([i64::MAX, i64::MAX, 0], 0, before);
    assert_eq!(side(i64::MAX).map(|x| x.to_opposite), Some(i64::MAX));
    assert_eq!(side(0).map(|x| x.to_opposite), Some(0));
    // Both sides at corner 2 at the top: c_2 = i64::MAX arcs, all cross the new edge.
    let x = formulas::split_side([0, i64::MAX, i64::MAX], 0, 0).unwrap();
    assert_eq!(x.to_opposite, i64::MAX);
}

#[test]
fn eq3_below_minus_one_is_invalid() {
    // Quad 2,0,0,0,0 (SPEC §4.2): two arcs from k to l; Eq. 3 gives −2.
    assert_eq!(formulas::flip_normal(2, 0, 0, 0, 0), Ok(-2));
    let (mut it, diag) = square(|_, d| if d { 2 } else { 0 });
    let before = it.clone();
    match it.flip_edge(diag) {
        Err(IntegerError::InvalidCoordinates(x)) => {
            assert_eq!(x.place, ddg_integer_coords::Place::Edge(diag as u32))
        }
        other => panic!("{other:?}"),
    }
    assert!(it == before);
}

// ----- timing (SPEC §10.1) -----

/// Delaunay flipping and random flips, integer coordinates against signposts, on spot and
/// camelbad. Run alone, single thread, release:
/// `cargo test --release -p ddg-integer-coords --test ops -- --ignored --nocapture timing`.
#[test]
#[ignore]
fn timing() {
    use std::time::Instant;
    let load = std::fs::read_to_string("/proc/loadavg").unwrap_or_default();
    eprintln!("load average before: {}", load.trim());
    let best = |runs: usize, mut f: Box<dyn FnMut() -> f64 + '_>| {
        (0..runs).map(|_| f()).fold(f64::INFINITY, f64::min)
    };
    for path in [
        "models/spot/spot/spot_triangulated.obj",
        "optimal-cones/optimal-cones-data/data/camelbad.obj",
    ] {
        let Some(shape) = common::load_data(path) else {
            eprintln!("{path}: not present");
            continue;
        };
        let s0 = shape.triangulation();
        let i0 = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let opts = DelaunayOptions::default();
        let mut n_delaunay = 0;
        let t_s = best(
            5,
            Box::new(|| {
                let mut s = s0.clone();
                let t = Instant::now();
                let r = s.flip_to_delaunay(&opts).unwrap();
                n_delaunay = r.n_flips;
                t.elapsed().as_secs_f64()
            }),
        );
        let t_i = best(
            5,
            Box::new(|| {
                let mut it = i0.clone();
                let t = Instant::now();
                Triangulation::flip_to_delaunay(&mut it, &opts).unwrap();
                t.elapsed().as_secs_f64()
            }),
        );
        let mut it = i0.clone();
        Triangulation::flip_to_delaunay(&mut it, &opts).unwrap();
        let t_c = best(
            5,
            Box::new(|| {
                let t = Instant::now();
                it.certify().unwrap();
                t.elapsed().as_secs_f64()
            }),
        );
        // Random flips: the same 10 |E| attempts on both.
        let ne = s0.n_edges();
        let mut rng = SplitMix64::new(1);
        let edges: Vec<usize> = (0..10 * ne)
            .map(|_| (rng.next_u64() % ne as u64) as usize)
            .collect();
        let mut done_s = 0;
        let r_s = best(
            3,
            Box::new(|| {
                let mut s = s0.clone();
                let t = Instant::now();
                done_s = edges.iter().filter(|&&e| s.flip_edge(e).is_ok()).count();
                t.elapsed().as_secs_f64()
            }),
        );
        let mut done_i = 0;
        let mut sum_n = 0u128;
        let r_i = best(
            3,
            Box::new(|| {
                let mut it = i0.clone();
                let t = Instant::now();
                done_i = edges.iter().filter(|&&e| it.flip_edge(e).is_ok()).count();
                let dt = t.elapsed().as_secs_f64();
                sum_n = it.total_crossings();
                dt
            }),
        );
        assert_eq!(done_s, done_i);
        eprintln!(
            "{}: |E| {ne}; Delaunay flipping ({n_delaunay} flips): signposts {:.2} ms, integer \
             {:.2} ms ({:.2}x), certify after it {:.2} ms; random flips ({done_i} of {} \
             attempts, Σn⁺ after {sum_n}): signposts {:.2} M/s, integer {:.2} M/s ({:.0} / {:.0} \
             ns per flip)",
            shape.name,
            1e3 * t_s,
            1e3 * t_i,
            t_i / t_s,
            1e3 * t_c,
            edges.len(),
            done_s as f64 / r_s / 1e6,
            done_i as f64 / r_i / 1e6,
            1e9 * r_s / done_s as f64,
            1e9 * r_i / done_i as f64,
        );
        // Removal: insertions by the signposts, converted, then every inserted vertex removed.
        let s1 = state(&shape, Scenario::Insertions, 5);
        let (i1, _) = convert(&s1).unwrap();
        let mut n_removed = 0;
        let mut n_flips = 0;
        let t_r = best(
            3,
            Box::new(|| {
                let mut it = i1.clone();
                let n_input = it.n_input_vertices();
                let t = Instant::now();
                let (mut removed, mut flips) = (0, 0);
                let mut v = it.complex().n_vertices();
                while v > n_input {
                    v -= 1;
                    if let Ok(r) = it.remove_vertex(v) {
                        removed += 1;
                        flips += r.flips.len();
                    }
                }
                n_removed = removed;
                n_flips = flips;
                t.elapsed().as_secs_f64()
            }),
        );
        eprintln!(
            "{}: removal of {n_removed} inserted vertices ({n_flips} flips): {:.1} µs per removal",
            shape.name,
            1e6 * t_r / n_removed.max(1) as f64
        );
    }
    let load = std::fs::read_to_string("/proc/loadavg").unwrap_or_default();
    eprintln!("load average after: {}", load.trim());
}
