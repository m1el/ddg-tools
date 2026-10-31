//! The geometry of the crossings (WP2): `ExtractIntrinsicEdge` against the signposts' traces,
//! the two directions of a crossing against each other (FlipOut contract R2), the positions of
//! the input edges against the signposts' overlay, clamps (SPEC I5), and bad data.
//!
//! The states are those of `tests/states.rs`: signpost triangulations of the standard shapes,
//! spot and camelbad in 7 scenarios and 2 seeds, converted to integer coordinates. Run with
//! `--nocapture` to see the tables.

mod common;

use std::sync::OnceLock;

use common::{SCENARIOS, SEEDS, Scenario, all_shapes, convert, state};
use ddg_core::generators::SplitMix64;
use ddg_integer_coords::{
    Crossing, CurveSegment, IntegerError, IntegerTriangulation, Parts, PathPiece, StripReport,
};
use ddg_intrinsic::{
    EdgeTrace, EdgeTraceOptions, EdgeTraces, IntrinsicTriangulation, SurfacePoint,
};

/// Positions must agree to this, as parameters (fractions of the edge's length).
const TOLERANCE: f64 = 1e-9;

struct Job {
    name: String,
    signposts: IntrinsicTriangulation,
    traces: EdgeTraces,
    tri: IntegerTriangulation,
}

/// Every (shape, scenario, seed) that converts, computed once per test binary, in parallel.
fn jobs() -> &'static [Job] {
    static JOBS: OnceLock<Vec<Job>> = OnceLock::new();
    JOBS.get_or_init(|| {
        let shapes = all_shapes();
        let mut specs = Vec::new();
        for s in 0..shapes.len() {
            for scenario in SCENARIOS {
                for seed in SEEDS {
                    specs.push((s, scenario, seed));
                }
            }
        }
        par_map(&specs, |&(s, scenario, seed): &(usize, Scenario, u64)| {
            let shape = &shapes[s];
            let signposts = state(shape, scenario, seed);
            let (tri, _) = convert(&signposts).ok()?;
            let traces = signposts
                .trace_edges(&EdgeTraceOptions::default())
                .expect("the conversion traced every edge");
            Some(Job {
                name: format!("{} / {scenario:?} / {seed}", shape.name),
                signposts,
                traces,
                tri,
            })
        })
        .into_iter()
        .flatten()
        .collect()
    })
}

/// `f` over `items` on all cores, results in order.
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

/// The largest value seen and where.
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

/// Number of self-edges and of faces that use a vertex twice.
fn delta_complex_features(t: &IntegerTriangulation) -> (usize, usize) {
    let c = t.complex();
    let self_edges = (0..c.n_edges())
        .filter(|&e| {
            let [a, b] = c.edge_vertices(e);
            a == b
        })
        .count();
    let repeated = c
        .faces()
        .iter()
        .filter(|f| f[0] == f[1] || f[1] == f[2] || f[2] == f[0])
        .count();
    (self_edges, repeated)
}

/// The signposts' crossings of edge `e`, without a last crossing that is the arrival at a far
/// vertex lying on that input edge (as the conversion reads it); and whether one was dropped.
fn signpost_crossings(job: &Job, e: usize) -> (&[ddg_intrinsic::Crossing], bool) {
    let c = job.traces.crossings(e);
    let b = job
        .tri
        .complex()
        .to_vertex(job.tri.complex().edge_halfedge(e));
    if let (SurfacePoint::Edge { edge, .. }, Some(last)) = (job.tri.vertex_location(b), c.last())
        && last.edge == edge
    {
        return (&c[..c.len() - 1], true);
    }
    (c, false)
}

#[derive(Default)]
struct EdgeStats {
    n_edges: usize,
    n_along: usize,
    n_crossings: usize,
    dropped: usize,
    worst_t: Worst,
    worst_distance: Worst,
    worst_mismatch: Worst,
    /// `|Δt| sin θ`: the distance between the two crossing points normal to the intrinsic edge,
    /// relative to the input edge's length.
    worst_conditioned: Worst,
    report: StripReport,
    /// Per crossing: (edge, index, |Δt|), judged once the crossing angles are known.
    raw: Vec<(usize, usize, f64)>,
    /// Crossings with `|Δt| > 1e-9` whose conditioned difference is within the tolerance.
    ill_conditioned: Vec<String>,
    disagreements: Vec<String>,
}

impl EdgeStats {
    /// Judges the position differences with the crossing angles `sins[e][i]`.
    fn judge(&mut self, name: &str, sins: &[Vec<f64>]) {
        for &(e, i, dt) in &self.raw {
            let Some(&sin) = sins.get(e).and_then(|s| s.get(i)) else {
                continue;
            };
            let conditioned = dt * sin;
            self.worst_conditioned
                .note(conditioned, || format!("{name}: edge {e}, crossing {i}"));
            let line = format!(
                "{name}: edge {e}, crossing {i}: |Δt| = {dt:.2e}, sin θ = {sin:.2e}, |Δt| sin θ = \
                 {conditioned:.2e}"
            );
            if conditioned > TOLERANCE {
                self.disagreements.push(line);
            } else if dt > TOLERANCE {
                self.ill_conditioned.push(line);
            }
        }
    }
}

/// Every edge's `extract_intrinsic_edge` against the signposts' `trace_edge` of the same edge.
fn compare_edges(job: &Job) -> (EdgeStats, Vec<EdgeTrace>) {
    let options = EdgeTraceOptions::default();
    let t = &job.tri;
    let mut s = EdgeStats::default();
    let mut ours_all = Vec::with_capacity(t.complex().n_edges());
    for e in 0..t.complex().n_edges() {
        s.n_edges += 1;
        let place = || format!("{}: edge {e}", job.name);
        let (ours, report) = match t.extract_intrinsic_edge_with_report(e, &options) {
            Ok(x) => x,
            Err(err) => {
                s.disagreements.push(format!("{}: error {err}", place()));
                continue;
            }
        };
        s.report.merge(&report);
        s.worst_mismatch
            .note(ours.mismatch / t.edge_lengths()[e], place);
        let (theirs, dropped) = signpost_crossings(job, e);
        s.dropped += usize::from(dropped);
        if ours.vertices != job.signposts.edge_vertices(e).map(|v| v as u32) {
            s.disagreements
                .push(format!("{}: end points differ", place()));
        }
        if ours.along_input_edge != job.traces.along_input_edge(e) {
            s.disagreements.push(format!(
                "{}: along {:?}, signposts {:?}",
                place(),
                ours.along_input_edge,
                job.traces.along_input_edge(e)
            ));
        }
        if ours.along_input_edge.is_some() {
            s.n_along += 1;
        } else if ours.first_face != job.traces.first_face(e) {
            s.disagreements.push(format!(
                "{}: first face {:?}, signposts {:?}",
                place(),
                ours.first_face,
                job.traces.first_face(e)
            ));
        }
        let key = |c: &ddg_intrinsic::Crossing| (c.halfedge, c.edge);
        let same = ours.crossings.len() == theirs.len()
            && ours
                .crossings
                .iter()
                .zip(theirs)
                .all(|(a, b)| key(a) == key(b));
        if !same {
            s.disagreements.push(format!(
                "{}: crossed input halfedges {:?}, signposts {:?}",
                place(),
                ours.crossings.iter().map(key).collect::<Vec<_>>(),
                theirs.iter().map(key).collect::<Vec<_>>()
            ));
        } else {
            let length = t.edge_lengths()[e];
            for (i, (a, b)) in ours.crossings.iter().zip(theirs).enumerate() {
                s.n_crossings += 1;
                let dt = (a.t - b.t).abs();
                s.worst_t.note(dt, || format!("{}, crossing {i}", place()));
                s.raw.push((e, i, dt));
                s.worst_distance
                    .note((a.distance - b.distance).abs() / length, || {
                        format!("{}, crossing {i}", place())
                    });
            }
        }
        ours_all.push(ours);
    }
    (s, ours_all)
}

#[test]
fn intrinsic_edges_match_the_signpost_traces_and_both_directions_agree() {
    let jobs = jobs();
    assert!(jobs.len() >= 200, "{} states converted", jobs.len());
    let results = par_map(jobs, |job| {
        let (mut edges, traces) = compare_edges(job);
        let (directions, sins) = compare_directions(job, &traces);
        edges.judge(&job.name, &sins);
        (edges, directions)
    });
    let mut total = EdgeStats::default();
    let mut dir = DirectionStats::default();
    let mut n_self_edge_states = 0;
    let (mut n_self_edges, mut n_repeated) = (0, 0);
    eprintln!(
        "{:<52} {:>6} {:>7} {:>8} {:>9} {:>9} {:>9} {:>4} {:>4}",
        "state", "|E|", "along", "cross", "Δt", "Δu,Δv", "mismatch", "self", "clmp"
    );
    for (job, (s, d)) in jobs.iter().zip(&results) {
        let (se, rep) = delta_complex_features(&job.tri);
        n_self_edges += se;
        n_repeated += rep;
        n_self_edge_states += usize::from(se > 0);
        eprintln!(
            "{:<52} {:>6} {:>7} {:>8} {:>9.1e} {:>9.1e} {:>9.1e} {:>4} {:>4}",
            job.name,
            s.n_edges,
            s.n_along,
            s.n_crossings,
            s.worst_t.value,
            d.worst.value,
            s.worst_mismatch.value,
            se,
            s.report.clamped_crossings + d.report.clamped_crossings,
        );
        total.n_edges += s.n_edges;
        total.n_along += s.n_along;
        total.n_crossings += s.n_crossings;
        total.dropped += s.dropped;
        total.worst_t.merge(&s.worst_t);
        total.worst_distance.merge(&s.worst_distance);
        total.worst_mismatch.merge(&s.worst_mismatch);
        total.worst_conditioned.merge(&s.worst_conditioned);
        total.report.merge(&s.report);
        total.disagreements.extend(s.disagreements.iter().cloned());
        total
            .ill_conditioned
            .extend(s.ill_conditioned.iter().cloned());
        dir.merge(d);
    }
    eprintln!(
        "{} states ({} with self-edges: {} self-edges, {} faces using a vertex twice)",
        jobs.len(),
        n_self_edge_states,
        n_self_edges,
        n_repeated
    );
    eprintln!(
        "ExtractIntrinsicEdge vs signpost trace_edge: {} edges ({} along input edges), {} \
         crossings, {} signpost arrival crossings dropped; worst |Δt| {:.2e} at {}; worst \
         |Δt| sin θ {:.2e} at {}; worst |Δdistance|/ℓ {:.2e} at {}; worst relative length \
         mismatch {:.2e} at {}; clamped {} (largest {:.1e}); {} walk steps",
        total.n_edges,
        total.n_along,
        total.n_crossings,
        total.dropped,
        total.worst_t.value,
        total.worst_t.place,
        total.worst_conditioned.value,
        total.worst_conditioned.place,
        total.worst_distance.value,
        total.worst_distance.place,
        total.worst_mismatch.value,
        total.worst_mismatch.place,
        total.report.clamped_crossings,
        total.report.max_clamp,
        total.report.walk_steps,
    );
    eprintln!(
        "Alg. 2 (strip over T1) vs ExtractIntrinsicEdge (strip over T0): {} crossings of {} \
         input edges, worst |Δu| {:.2e} at {}, worst |Δv| {:.2e} at {}; worst max(|Δu|, \
         |Δv|) sin θ {:.2e} at {}; worst relative length mismatch of a piece {:.2e}; clamped {}",
        dir.n_crossings,
        dir.n_curves,
        dir.worst_u.value,
        dir.worst_u.place,
        dir.worst_v.value,
        dir.worst_v.place,
        dir.worst_conditioned.value,
        dir.worst_conditioned.place,
        dir.report.max_relative_mismatch,
        dir.report.clamped_crossings,
    );
    for d in total.ill_conditioned.iter().chain(&dir.ill_conditioned) {
        eprintln!("  above 1e-9 only by conditioning: {d}");
    }
    for d in total
        .disagreements
        .iter()
        .chain(&dir.disagreements)
        .take(40)
    {
        eprintln!("  DISAGREES {d}");
    }
    assert!(
        total.disagreements.is_empty() && dir.disagreements.is_empty(),
        "{} disagreements with the signposts, {} between the directions",
        total.disagreements.len(),
        dir.disagreements.len()
    );
    assert_eq!(total.report.clamped_crossings, 0);
    assert_eq!(dir.report.clamped_crossings, 0);
    assert!(n_self_edge_states > 0, "no Δ-complex state was tested");
}

#[derive(Default)]
struct DirectionStats {
    n_curves: usize,
    n_crossings: usize,
    worst_u: Worst,
    worst_v: Worst,
    worst: Worst,
    worst_conditioned: Worst,
    report: StripReport,
    ill_conditioned: Vec<String>,
    disagreements: Vec<String>,
}

impl DirectionStats {
    fn merge(&mut self, o: &DirectionStats) {
        self.n_curves += o.n_curves;
        self.n_crossings += o.n_crossings;
        self.worst_u.merge(&o.worst_u);
        self.worst_v.merge(&o.worst_v);
        self.worst.merge(&o.worst);
        self.worst_conditioned.merge(&o.worst_conditioned);
        self.report.merge(&o.report);
        self.ill_conditioned
            .extend(o.ill_conditioned.iter().cloned());
        self.disagreements.extend(o.disagreements.iter().cloned());
    }
}

/// R2: every crossing of an input edge with an intrinsic edge, from Alg. 2's strip over `T¹`
/// (all input edges, `extract_edge_geometry`) and from `extract_intrinsic_edge`'s strip over
/// `T⁰`: the same input edge, the same point (`u` along the input edge, `v` along the intrinsic
/// edge). Every crossing of every intrinsic edge is met exactly once. Also returns the sine of
/// the crossing angle of every crossing, by intrinsic edge.
fn compare_directions(job: &Job, traces: &[EdgeTrace]) -> (DirectionStats, Vec<Vec<f64>>) {
    let t = &job.tri;
    let c = t.complex();
    let mesh = t.input().mesh();
    let mut s = DirectionStats::default();
    if traces.len() != c.n_edges() {
        s.disagreements
            .push(format!("{}: an intrinsic edge failed", job.name));
        return (s, Vec::new());
    }
    let mut met: Vec<Vec<bool>> = traces
        .iter()
        .map(|x| vec![false; x.crossings.len()])
        .collect();
    let mut sins: Vec<Vec<f64>> = traces
        .iter()
        .map(|x| vec![0.0; x.crossings.len()])
        .collect();
    for ie in 0..mesh.n_edges() {
        let g = match t.extract_edge_geometry(mesh.edge_halfedge(ie)) {
            Ok(g) => g,
            Err(err) => {
                s.disagreements
                    .push(format!("{}: input edge {ie}: {err}", job.name));
                continue;
            }
        };
        s.report.merge(&g.report);
        s.n_curves += 1;
        for (piece, crossings) in g.path.pieces.iter().zip(&g.crossings) {
            if let PathPiece::Crossings(seg) = piece {
                assert_eq!(seg.crossings.len(), crossings.len());
            }
            for x in crossings {
                s.n_crossings += 1;
                let h = x.crossing.halfedge as usize;
                let e = c.edge(h);
                let n = traces[e].crossings.len() as i64;
                let i = if c.is_forward(h) {
                    x.crossing.index
                } else {
                    n - 1 - x.crossing.index
                } as usize;
                let place = || format!("{}: input edge {ie}, edge {e} crossing {i}", job.name);
                let Some(y) = traces[e].crossings.get(i) else {
                    s.disagreements
                        .push(format!("{}: no such crossing", place()));
                    continue;
                };
                if met[e][i] {
                    s.disagreements.push(format!("{}: met twice", place()));
                }
                met[e][i] = true;
                sins[e][i] = x.sin_angle;
                if y.edge != x.input_edge {
                    s.disagreements.push(format!(
                        "{}: input edge {} from the intrinsic edge",
                        place(),
                        y.edge
                    ));
                    continue;
                }
                let du = (x.u - y.t).abs();
                let SurfacePoint::Edge { t: v_edge, .. } = x.on_triangulation(t) else {
                    unreachable!()
                };
                let dv = (v_edge - y.distance / t.edge_lengths()[e]).abs();
                s.worst_u.note(du, place);
                s.worst_v.note(dv, place);
                s.worst.note(du.max(dv), place);
                let conditioned = du.max(dv) * x.sin_angle;
                s.worst_conditioned.note(conditioned, place);
                let line = format!(
                    "{}: |Δu| = {du:.2e}, |Δv| = {dv:.2e}, sin θ = {:.2e}",
                    place(),
                    x.sin_angle
                );
                if conditioned > TOLERANCE {
                    s.disagreements.push(line);
                } else if du.max(dv) > TOLERANCE {
                    s.ill_conditioned.push(line);
                }
            }
        }
    }
    let unmet = met.iter().flatten().filter(|&&m| !m).count();
    if unmet > 0 {
        s.disagreements.push(format!(
            "{}: {unmet} crossings of intrinsic edges not met",
            job.name
        ));
    }
    (s, sins)
}

/// The positions of the input edges' crossings (Alg. 2) against the signposts' overlay, whose
/// crossing vertices are the points where the traced intrinsic edges cross input edges.
#[test]
fn input_edge_positions_match_the_signpost_overlay() {
    let jobs = jobs();
    let results = par_map(jobs, |job| {
        let mut worst = Worst::default();
        let mut worst_conditioned = Worst::default();
        let mut bad = Vec::new();
        let mut ill = Vec::new();
        let mut matched = 0usize;
        let overlay = match job.signposts.overlay(&job.traces) {
            Ok(o) => o,
            Err(err) => {
                let bad = vec![format!("{}: overlay: {err}", job.name)];
                return (worst, worst_conditioned, bad, ill, 0, true);
            }
        };
        let t = &job.tri;
        let c = t.complex();
        let mut offsets = vec![0usize; c.n_edges() + 1];
        for e in 0..c.n_edges() {
            offsets[e + 1] = offsets[e] + job.traces.crossings(e).len();
        }
        let mesh = t.input().mesh();
        for ie in 0..mesh.n_edges() {
            let g = match t.extract_edge_geometry(mesh.edge_halfedge(ie)) {
                Ok(g) => g,
                Err(err) => {
                    bad.push(format!("{}: input edge {ie}: {err}", job.name));
                    continue;
                }
            };
            for x in g.crossings.iter().flatten() {
                let h = x.crossing.halfedge as usize;
                let e = c.edge(h);
                let n = job.traces.crossings(e).len() as i64;
                let i = if c.is_forward(h) {
                    x.crossing.index
                } else {
                    n - 1 - x.crossing.index
                } as usize;
                let vertex = &overlay.vertices[c.n_vertices() + offsets[e] + i];
                let (
                    SurfacePoint::Edge { edge: oe, t: ou },
                    SurfacePoint::Edge { edge: te, t: tv },
                ) = (vertex.on_input, vertex.on_triangulation)
                else {
                    bad.push(format!(
                        "{}: edge {e} crossing {i}: not an edge point",
                        job.name
                    ));
                    continue;
                };
                let SurfacePoint::Edge { t: v_edge, .. } = x.on_triangulation(t) else {
                    unreachable!()
                };
                if oe != x.input_edge || te as usize != e {
                    bad.push(format!(
                        "{}: edge {e} crossing {i}: overlay edges ({oe}, {te}), ours ({}, {e})",
                        job.name, x.input_edge
                    ));
                    continue;
                }
                matched += 1;
                let d = (ou - x.u).abs().max((tv - v_edge).abs());
                worst.note(d, || format!("{}: edge {e} crossing {i}", job.name));
                worst_conditioned.note(d * x.sin_angle, || {
                    format!("{}: edge {e} crossing {i}", job.name)
                });
                let line = format!(
                    "{}: edge {e} crossing {i}: (u, v) = ({}, {v_edge}), overlay ({ou}, {tv}), \
                     sin θ = {:.2e}",
                    job.name, x.u, x.sin_angle
                );
                if d * x.sin_angle > TOLERANCE {
                    bad.push(line);
                } else if d > TOLERANCE {
                    ill.push(line);
                }
            }
        }
        (worst, worst_conditioned, bad, ill, matched, false)
    });
    let mut worst = Worst::default();
    let mut worst_conditioned = Worst::default();
    let (mut matched, mut no_overlay) = (0, 0);
    let (mut bad, mut ill) = (Vec::new(), Vec::new());
    for (w, wc, b, i, m, refused) in &results {
        worst.merge(w);
        worst_conditioned.merge(wc);
        matched += m;
        no_overlay += usize::from(*refused);
        bad.extend(b.iter().cloned());
        ill.extend(i.iter().cloned());
    }
    eprintln!(
        "Alg. 2 positions vs the signposts' overlay: {matched} crossings matched in {} states \
         ({no_overlay} overlays refused); worst max(|Δu|, |Δv|) {:.2e} at {}; worst \
         max(|Δu|, |Δv|) sin θ {:.2e} at {}",
        jobs.len() - no_overlay,
        worst.value,
        worst.place,
        worst_conditioned.value,
        worst_conditioned.place
    );
    for b in &ill {
        eprintln!("  above 1e-9 only by conditioning: {b}");
    }
    for b in bad.iter().take(30) {
        eprintln!("  DISAGREES {b}");
    }
    assert!(bad.is_empty(), "{} disagreements", bad.len());
}

/// `extract_geometric_crossing` returns the crossing of `extract_geometric_curve`, and the
/// curve's ends lie at the parameters of their places.
#[test]
fn geometric_crossings_and_curves_are_consistent() {
    let jobs = jobs();
    let mut checked = 0;
    for job in jobs.iter().step_by(7) {
        let t = &job.tri;
        let c = t.complex();
        let mut rng = SplitMix64::new(5);
        for _ in 0..20 {
            let h = (rng.next_u64() % c.n_halfedges() as u64) as usize;
            let n = t.halfedge_normal(h);
            if n <= 0 {
                continue;
            }
            let z = Crossing::new(h, (rng.next_u64() % n as u64) as i64);
            let curve = t.extract_geometric_curve(z).unwrap();
            let (x, _) = t.extract_geometric_crossing(z).unwrap();
            assert!(curve.crossings.contains(&x));
            assert_eq!(x.crossing, z);
            assert_eq!(x.input_edge, curve.traced.input_edge);
            assert!((0.0..=1.0).contains(&x.u) && (0.0..=1.0).contains(&x.v));
            // u runs monotonically from u_start to u_end.
            let (lo, hi) = (
                curve.u_start.min(curve.u_end),
                curve.u_start.max(curve.u_end),
            );
            assert!(curve.crossings.iter().all(|y| y.u >= lo && y.u <= hi));
            checked += 1;
        }
    }
    eprintln!("{checked} random crossings");
    assert!(checked > 100);
}

/// Mutated integers: every extraction returns a result or an error, never panics.
#[test]
fn bad_data_never_panics() {
    let jobs = jobs();
    let options = EdgeTraceOptions::default();
    let mut rng = SplitMix64::new(77);
    let (mut n_ok, mut n_err) = (0usize, 0usize);
    // Small states only: every mutation re-extracts everything.
    let small: Vec<&Job> = jobs
        .iter()
        .filter(|j| j.tri.complex().n_edges() < 400)
        .step_by(3)
        .collect();
    assert!(!small.is_empty());
    for job in &small {
        for _ in 0..6 {
            let mut p: Parts = job.tri.clone().into_parts();
            let ne = p.normal.len();
            let e = (rng.next_u64() % ne as u64) as usize;
            match rng.next_u64() % 4 {
                0 => p.normal[e] += 1,
                1 => p.normal[e] = (p.normal[e] - 1).max(-1),
                2 => p.normal[e] += 2,
                _ => {
                    let h = (rng.next_u64() % p.roundabouts.len() as u64) as usize;
                    if p.roundabouts[h] != ddg_integer_coords::NO_ROUNDABOUT {
                        p.roundabouts[h] += 1;
                    }
                }
            }
            let m = IntegerTriangulation::from_parts(p).unwrap();
            for e in 0..m.complex().n_edges() {
                match m.extract_intrinsic_edge(e, &options) {
                    Ok(_) => n_ok += 1,
                    Err(_) => n_err += 1,
                }
            }
            let mesh = m.input().mesh();
            for ie in 0..mesh.n_edges() {
                match m.extract_edge_geometry(mesh.edge_halfedge(ie)) {
                    Ok(_) => n_ok += 1,
                    Err(_) => n_err += 1,
                }
            }
        }
    }
    eprintln!("mutated states: {n_ok} extractions returned, {n_err} errors");
    assert!(n_err > 0);
    // Out of range, bad options, huge counts.
    let t = &small
        .iter()
        .find(|j| j.tri.total_crossings() > 0)
        .expect("a small state with crossings")
        .tri;
    assert!(matches!(
        t.extract_intrinsic_edge(t.complex().n_edges(), &options),
        Err(IntegerError::IndexOutOfRange { .. })
    ));
    let nan = EdgeTraceOptions {
        tolerance: f64::NAN,
        max_crossings: None,
    };
    assert!(matches!(
        t.extract_intrinsic_edge(0, &nan),
        Err(IntegerError::Intrinsic(_))
    ));
    let mut p = t.clone().into_parts();
    let e = (0..p.normal.len())
        .find(|&e| !p.complex.is_boundary_edge(e))
        .unwrap();
    p.normal[e] = 1 << 40;
    let m = IntegerTriangulation::from_parts(p).unwrap();
    assert!(m.extract_intrinsic_edge(e, &options).is_err());
    // A hand-built segment that is not a chain of faces.
    let h = (0..t.complex().n_halfedges())
        .find(|&h| t.halfedge_normal(h) > 0)
        .unwrap();
    let curve = t.extract_geometric_curve(Crossing::new(h, 0)).unwrap();
    let mut broken: CurveSegment = curve.segment.clone();
    broken.start = (broken.start + 3) % t.complex().n_halfedges() as u32;
    assert!(
        t.segment_geometry(&broken, curve.traced.input_edge as usize, 0.0, 1.0)
            .is_err()
    );
    let mut broken = curve.segment.clone();
    broken.end.corner = (broken.end.corner + 3) % t.complex().n_halfedges() as u32;
    assert!(
        t.segment_geometry(&broken, curve.traced.input_edge as usize, 0.0, 1.0)
            .is_err()
    );
    let mut empty = curve.segment.clone();
    empty.crossings.clear();
    assert!(t.segment_geometry(&empty, 0, 0.0, 1.0).is_err());
}

/// The shared interface's `trace_edge_on_input` is `extract_intrinsic_edge`.
#[test]
fn the_interface_traces_with_extract_intrinsic_edge() {
    use ddg_intrinsic::{IntrinsicError, Triangulation};
    let options = EdgeTraceOptions::default();
    let mut n = 0;
    for job in jobs().iter().step_by(5) {
        let t = &job.tri;
        for e in (0..t.complex().n_edges()).step_by(7) {
            let via_trait = Triangulation::trace_edge_on_input(t, e, &options).unwrap();
            assert_eq!(via_trait, t.extract_intrinsic_edge(e, &options).unwrap());
            n += 1;
        }
        let ne = t.complex().n_edges();
        assert!(matches!(
            Triangulation::trace_edge_on_input(t, ne, &options),
            Err(IntrinsicError::IndexOutOfRange { .. })
        ));
    }
    assert!(n > 1000, "{n}");
}
