//! Face splits, App. A and moves (WP4; SPEC §4.4, §4.7, §4.9, I6, I15, I16), differentially
//! against the signpost representation: the same insertion is done in an
//! `IntrinsicTriangulation` and in the `IntegerTriangulation` converted from its state, and the
//! integer state must equal `IntegerTriangulation::from_intrinsic` of the signposts exactly
//! (connectivity, lengths bitwise, every `n`, every roundabout); the places `q⁰` are compared
//! by distance. Where they differ, the distance of the point to the nearest arc (I6) says why.
//!
//! Run with `--nocapture` to see the tables.

mod common;

use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use common::{SCENARIOS, Scenario, Shape, all_shapes, convert, inner_bary, standard_shapes, state};
use ddg_core::generators::SplitMix64;
use ddg_integer_coords::insert::FaceArrangement;
use ddg_integer_coords::{
    FaceSplitOptions, IntegerError, IntegerTriangulation, MoveInsertion, Region,
    recover_barycentric,
};
use ddg_intrinsic::{DelaunayOptions, IntrinsicError, MoveRefusal, SurfacePoint, Triangulation};

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

/// The first difference between two integer states (places not compared), or `None`.
fn differs(a: &IntegerTriangulation, b: &IntegerTriangulation) -> Option<String> {
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
    None
}

/// The position in space of a point of the input mesh of `shape`.
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

/// The distance in space of two points of the input mesh of `shape`, over its mean edge length.
fn place_distance(shape: &Shape, a: SurfacePoint, b: SurfacePoint) -> f64 {
    common::dist(position(shape, a), position(shape, b)) / shape.mean_edge_length()
}

/// The largest difference of the barycentric coordinates of two points of the same input face
/// (or parameters on the same input edge); infinite for points of different elements.
fn bary_difference(a: SurfacePoint, b: SurfacePoint) -> f64 {
    match (a, b) {
        (SurfacePoint::Face { face: f, bary: p }, SurfacePoint::Face { face: g, bary: q })
            if f == g =>
        {
            (0..3).map(|c| (p[c] - q[c]).abs()).fold(0.0, f64::max)
        }
        (SurfacePoint::Edge { edge: e, t: s }, SurfacePoint::Edge { edge: g, t }) if e == g => {
            (s - t).abs()
        }
        (SurfacePoint::Vertex(v), SurfacePoint::Vertex(w)) if v == w => 0.0,
        _ => f64::INFINITY,
    }
}

/// The jobs: every shape, scenario and seed.
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

/// `Σ n⁺` of the sides of face `f`.
fn face_crossings(t: &IntegerTriangulation, f: usize) -> i64 {
    t.face_normal(f).iter().map(|&n| n.max(0)).sum()
}

/// What kind of point to insert.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PointKind {
    Inner,
    NearEdge,
    NearCorner,
    Crowded,
    NearArc,
}

/// The parameters of the crossings of side `k` of face `f` from the side's tail, from the
/// integer representation's own `ExtractIntrinsicEdge`.
fn side_params(t: &IntegerTriangulation, f: usize, k: usize) -> Vec<f64> {
    let h = 3 * f + k;
    let e = t.complex().edge(h);
    if t.normal_coordinate(e) < 0 {
        return Vec::new();
    }
    let tr = t.extract_intrinsic_edge(e, &Default::default()).unwrap();
    let l = t.edge_lengths()[e];
    let p: Vec<f64> = tr.crossings.iter().map(|x| x.distance / l).collect();
    if t.complex().edge_halfedge(e) == h {
        p
    } else {
        p.iter().rev().map(|x| 1.0 - x).collect()
    }
}

/// A point at the relative distance `delta` (towards the corner or away from it) from the
/// middle of an arc around a corner of a random face that has one.
fn near_arc_point(t: &IntegerTriangulation, rng: &mut SplitMix64) -> Option<(usize, [f64; 3])> {
    let nf = t.complex().n_faces();
    for _ in 0..40 {
        let f = (rng.next_u64() % nf as u64) as usize;
        let Some(a) = FaceArrangement::of(t.face_normal(f)) else {
            continue;
        };
        let ks: Vec<usize> = (0..3).filter(|&k| a.coords.c[k] > 0).collect();
        if ks.is_empty() {
            continue;
        }
        let k = ks[(rng.next_u64() % ks.len() as u64) as usize];
        let i = (rng.next_u64() % a.coords.c[k] as u64) as usize;
        let (pk, pk2) = (side_params(t, f, k), side_params(t, f, (k + 2) % 3));
        let (u, w) = (pk[i], pk2[pk2.len() - 1 - i]);
        let mut p = [0.0; 3];
        p[k] += 0.5 * (1.0 - u) + 0.5 * w;
        p[(k + 1) % 3] += 0.5 * u;
        p[(k + 2) % 3] += 0.5 * (1.0 - w);
        let delta = 10f64.powf(-rng.uniform(3.0, 16.0))
            * if rng.next_u64().is_multiple_of(2) {
                1.0
            } else {
                -1.0
            };
        let mut q = [0.0; 3];
        for c in 0..3 {
            let corner = if c == k { 1.0 } else { 0.0 };
            q[c] = p[c] + delta * (corner - p[c]);
        }
        if q.iter().all(|&x| x > 0.0) {
            return Some((f, q));
        }
    }
    None
}

/// A face and a point in it: random inner points, points near a side or a corner (down to the
/// thinness `plan_face_split` refuses), and points in the face with the most crossings among a
/// few random ones.
fn pick_point(t: &IntegerTriangulation, rng: &mut SplitMix64) -> (usize, [f64; 3], PointKind) {
    let nf = t.complex().n_faces();
    let kind = match rng.next_u64() % 9 {
        0..=2 => PointKind::Inner,
        3 | 4 => PointKind::NearEdge,
        5 => PointKind::NearCorner,
        6 => PointKind::NearArc,
        _ => PointKind::Crowded,
    };
    if kind == PointKind::NearArc
        && let Some((f, q)) = near_arc_point(t, rng)
    {
        return (f, q, kind);
    }
    let mut f = (rng.next_u64() % nf as u64) as usize;
    let bary = match kind {
        PointKind::Inner | PointKind::NearArc => inner_bary(rng),
        PointKind::Crowded => {
            for _ in 0..30 {
                let g = (rng.next_u64() % nf as u64) as usize;
                if face_crossings(t, g) > face_crossings(t, f) {
                    f = g;
                }
            }
            inner_bary(rng)
        }
        PointKind::NearEdge => {
            let small = 10f64.powf(-rng.uniform(2.0, 7.0));
            let a = rng.uniform(0.1, 0.9) * (1.0 - small);
            let k = (rng.next_u64() % 3) as usize;
            let mut b = [0.0; 3];
            b[k] = small;
            b[(k + 1) % 3] = a;
            b[(k + 2) % 3] = 1.0 - small - a;
            b
        }
        PointKind::NearCorner => {
            let small = 10f64.powf(-rng.uniform(1.5, 5.0));
            let a = rng.uniform(0.1, 0.9) * small;
            let k = (rng.next_u64() % 3) as usize;
            let mut b = [0.0; 3];
            b[k] = 1.0 - small;
            b[(k + 1) % 3] = a;
            b[(k + 2) % 3] = small - a;
            b
        }
    };
    (f, bary, kind)
}

/// The number of corners App. A must use for `region` of `a` (SPEC I16).
fn expected_corners(a: &FaceArrangement, region: Region) -> usize {
    let c = a.coords.c;
    match region {
        Region::Corner { index: 0, .. } => 3,
        Region::Corner { .. } => 4,
        Region::Middle { .. } if a.fan == 0 => (0..3).map(|k| if c[k] > 0 { 2 } else { 1 }).sum(),
        Region::Middle { index } => {
            let (j1, j2) = ((a.kappa + 1) % 3, (a.kappa + 2) % 3);
            let mut n = 1;
            if index == 0 {
                n += if c[j1] > 0 { 2 } else { 1 };
            }
            if index > 0 {
                n += 1;
            }
            if index < a.fan {
                n += 1;
            }
            if index == a.fan {
                n += if c[j2] > 0 { 2 } else { 1 };
            }
            n
        }
    }
}

// ----- differential insertion -----

#[derive(Default, Debug)]
struct InsertStats {
    attempts: usize,
    inserted: usize,
    by_kind: [usize; 5],
    refusals_agree: usize,
    /// Signposts refused (their trace over the input mesh failed), integers inserted.
    signposts_refused: Vec<String>,
    /// Integers refused, signposts inserted.
    integer_refused: Vec<String>,
    /// Integer states that differ from the conversion, with the distance to the nearest arc.
    near_arc: Vec<(f64, String)>,
    /// The signposts snapped the place onto an input edge (the states then differ by design).
    snapped: usize,
    mismatches: Vec<String>,
    worst_place: f64,
    worst_place_where: String,
    min_curve_distance_equal: f64,
    corners_hist: [usize; 7],
    corners_wrong: Vec<String>,
    crowded_max: i64,
    in_repeated_faces: usize,
    moved: usize,
    worst_residual: f64,
    min_coordinate: f64,
    runner_up_max: f64,
    /// Insertions by the decade of their distance to the nearest arc (`10^-k`, `k < 17`):
    /// `[equal, differ]`.
    arc_decades: [[usize; 2]; 18],
}

fn decade(d: f64) -> usize {
    if d.is_nan() || d <= 0.0 || !d.is_finite() {
        return 17;
    }
    ((-d.log10()).floor().max(0.0) as usize).min(17)
}

#[test]
fn insertions_match_the_signposts_exactly() {
    let shapes = all_shapes();
    let specs = specs(&shapes);
    let results = par_map(&specs, |&(i, scenario, seed)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?} / {seed}", shape.name);
        let mut st = InsertStats {
            min_curve_distance_equal: f64::INFINITY,
            min_coordinate: f64::INFINITY,
            runner_up_max: f64::NEG_INFINITY,
            ..InsertStats::default()
        };
        let mut s = state(shape, scenario, seed);
        let Ok((mut it, _)) = convert(&s) else {
            st.mismatches
                .push(format!("{tag}: initial conversion failed"));
            return st;
        };
        let mut rng = SplitMix64::new(seed ^ 0x1a5e47);
        let large = shape.mesh.n_edges() > 5000;
        let attempts = if large { 40 } else { 120 };
        for a in 0..attempts {
            st.attempts += 1;
            let (f, bary, kind) = pick_point(&it, &mut rng);
            let before = it.clone();
            let crossings = face_crossings(&it, f);
            let [x, y, z] = it.complex().faces()[f];
            let repeated = x == y || y == z || z == x;
            let rs = s.insert_vertex(f, bary);
            let ri = it.split_face(f, bary, &FaceSplitOptions::default());
            match (rs, ri) {
                (Ok(xs), Ok(xi)) => {
                    st.inserted += 1;
                    st.by_kind[kind as usize] += 1;
                    if kind == PointKind::Crowded {
                        st.crowded_max = st.crowded_max.max(crossings);
                    }
                    st.in_repeated_faces += usize::from(repeated);
                    let loc = &xi.location;
                    let rho = xi.recover.n_corners;
                    st.corners_hist[rho.min(6)] += 1;
                    let want = expected_corners(&loc.arrangement, loc.choice.region);
                    if rho != want {
                        st.corners_wrong.push(format!(
                            "{tag}: region {} used {rho} corners, expected {want}",
                            loc.choice.region
                        ));
                    }
                    if xi.recover.moved > 0.0 {
                        st.moved += 1;
                    }
                    st.worst_residual = st.worst_residual.max(xi.recover.residual);
                    st.min_coordinate = st.min_coordinate.min(xi.recover.min_coordinate);
                    if let Some((_, g)) = loc.choice.runner_up {
                        st.runner_up_max = st.runner_up_max.max(g);
                    }
                    let distance = loc.choice.curve_distance.unwrap_or(f64::INFINITY);
                    if xs.snapped.is_some() {
                        st.snapped += 1;
                        it = convert(&s).unwrap().0;
                        continue;
                    }
                    let d = place_distance(shape, xs.location, xi.insertion.location);
                    if d > st.worst_place {
                        st.worst_place = d;
                        st.worst_place_where = format!(
                            "{tag}: face {f} {bary:?} ({kind:?}), region {}, {crossings} \
                             crossings, arc at {distance:e}",
                            loc.choice.region
                        );
                    }
                    match convert(&s) {
                        Ok((c, _)) => match differs(&it, &c) {
                            None => {
                                st.min_curve_distance_equal =
                                    st.min_curve_distance_equal.min(distance);
                                st.arc_decades[decade(distance)][0] += 1;
                            }
                            Some(diff) => {
                                st.arc_decades[decade(distance)][1] += 1;
                                st.near_arc.push((
                                    distance,
                                    format!(
                                        "{tag}: insertion {a} in face {f} at {bary:?} \
                                         ({kind:?}): {diff}; region {} margin {:e}, runner-up \
                                         {:?}",
                                        loc.choice.region, loc.choice.margin, loc.choice.runner_up
                                    ),
                                ));
                                it = c;
                            }
                        },
                        Err(err) => {
                            st.mismatches
                                .push(format!("{tag}: conversion after insertion {a}: {err}"));
                            return st;
                        }
                    }
                    if (!large || a % 10 == 9)
                        && let Err(err) = it.certify()
                    {
                        st.mismatches
                            .push(format!("{tag}: certify after insertion {a}: {err}"));
                        return st;
                    }
                }
                (Err(es), Err(IntegerError::Intrinsic(ei))) if es == ei => {
                    st.refusals_agree += 1;
                    if it != before {
                        st.mismatches
                            .push(format!("{tag}: a refused insertion changed the state"));
                    }
                }
                (Err(es), Ok(xi)) => {
                    st.signposts_refused.push(format!(
                        "{tag}: face {f} {bary:?} ({kind:?}): signposts: {es}; integer region {}",
                        xi.location.choice.region
                    ));
                    // Put the integer state back in step with the signposts.
                    it = before;
                }
                (Ok(_), Err(ei)) => {
                    st.integer_refused
                        .push(format!("{tag}: face {f} {bary:?} ({kind:?}): {ei}"));
                    if it != before {
                        st.mismatches
                            .push(format!("{tag}: a refused insertion changed the state"));
                    }
                    it = convert(&s).unwrap().0;
                }
                (Err(es), Err(ei)) => {
                    st.mismatches
                        .push(format!("{tag}: face {f}: refusals differ: {es} vs {ei}"));
                }
            }
        }
        if let Err(err) = it.certify() {
            st.mismatches.push(format!("{tag}: final certify: {err}"));
        }
        st
    });
    let mut t = InsertStats {
        min_curve_distance_equal: f64::INFINITY,
        min_coordinate: f64::INFINITY,
        runner_up_max: f64::NEG_INFINITY,
        ..InsertStats::default()
    };
    for r in results {
        t.attempts += r.attempts;
        t.inserted += r.inserted;
        for k in 0..5 {
            t.by_kind[k] += r.by_kind[k];
        }
        t.refusals_agree += r.refusals_agree;
        t.signposts_refused.extend(r.signposts_refused);
        t.integer_refused.extend(r.integer_refused);
        t.near_arc.extend(r.near_arc);
        t.snapped += r.snapped;
        t.mismatches.extend(r.mismatches);
        if r.worst_place > t.worst_place {
            t.worst_place = r.worst_place;
            t.worst_place_where = r.worst_place_where;
        }
        t.min_curve_distance_equal = t.min_curve_distance_equal.min(r.min_curve_distance_equal);
        for k in 0..7 {
            t.corners_hist[k] += r.corners_hist[k];
        }
        t.corners_wrong.extend(r.corners_wrong);
        t.crowded_max = t.crowded_max.max(r.crowded_max);
        t.in_repeated_faces += r.in_repeated_faces;
        t.moved += r.moved;
        t.worst_residual = t.worst_residual.max(r.worst_residual);
        t.min_coordinate = t.min_coordinate.min(r.min_coordinate);
        t.runner_up_max = t.runner_up_max.max(r.runner_up_max);
        for k in 0..18 {
            t.arc_decades[k][0] += r.arc_decades[k][0];
            t.arc_decades[k][1] += r.arc_decades[k][1];
        }
    }
    t.near_arc.sort_by(|a, b| a.0.total_cmp(&b.0));
    eprintln!(
        "insertion differential: {} attempts, {} inserted in both (inner {}, near a side {}, \
         near a corner {}, crowded {} with up to {} crossings in the face, near an arc {}; {} in \
         faces with a repeated vertex), {} refusals identical; {} states differ from the conversion, {} \
         signpost places snapped onto an input edge; signposts refused {} that the integers \
         inserted, integers refused {}; smallest arc distance with equal states {:e}; worst \
         place difference {:e} mean edge lengths ({}); App. A corners used 3..6: {:?}, {} wrong; \
         {} places moved inside (I16), smallest coordinate before {:e}, worst residual {:e}; \
         largest runner-up margin {:e}; {} mismatches",
        t.attempts,
        t.inserted,
        t.by_kind[0],
        t.by_kind[1],
        t.by_kind[2],
        t.by_kind[3],
        t.crowded_max,
        t.by_kind[4],
        t.in_repeated_faces,
        t.refusals_agree,
        t.near_arc.len(),
        t.snapped,
        t.signposts_refused.len(),
        t.integer_refused.len(),
        t.min_curve_distance_equal,
        t.worst_place,
        t.worst_place_where,
        &t.corners_hist[3..],
        t.corners_wrong.len(),
        t.moved,
        t.min_coordinate,
        t.worst_residual,
        t.runner_up_max,
        t.mismatches.len()
    );
    eprintln!("  by distance to the nearest arc (decade 10^-k: equal / differ):");
    for (k, [eq, df]) in t.arc_decades.iter().enumerate() {
        if eq + df > 0 {
            eprintln!("    1e-{k:<2} {eq:>6} {df:>4}");
        }
    }
    for (d, m) in &t.near_arc {
        eprintln!("  differs (arc at {d:e}): {m}");
    }
    for m in t.signposts_refused.iter().take(10) {
        eprintln!("  signposts refused: {m}");
    }
    for m in &t.integer_refused {
        eprintln!("  integer refused: {m}");
    }
    for m in t.corners_wrong.iter().take(10) {
        eprintln!("  {m}");
    }
    for m in &t.mismatches {
        eprintln!("  {m}");
    }
    assert!(t.mismatches.is_empty());
    assert!(t.corners_wrong.is_empty());
    assert!(t.integer_refused.is_empty());
    assert!(t.inserted > 3000);
    // A difference from the signposts is allowed only within rounding of an arc.
    for (d, m) in &t.near_arc {
        assert!(*d < 1e-6, "differs far from every arc: {m}");
    }
    assert!(t.worst_place < 1e-8, "{}", t.worst_place_where);
    // The runner-up's margin is never positive (zero only for a point on an arc as computed).
    assert!(t.runner_up_max <= 0.0);
}

// ----- many insertions, Delaunay flipping, removal -----

#[test]
fn many_insertions_then_delaunay_then_removal() {
    let shapes = all_shapes();
    let specs: Vec<(usize, Scenario, u64)> = (0..shapes.len())
        .flat_map(|i| {
            [Scenario::Initial, Scenario::RandomFlips, Scenario::Mixed]
                .into_iter()
                .map(move |sc| (i, sc, 9u64))
        })
        .collect();
    #[derive(Default)]
    struct Stats {
        inserted: usize,
        refused: usize,
        delaunay_flips: usize,
        removed: usize,
        stuck: usize,
        boundary: usize,
        failures: Vec<String>,
    }
    let results = par_map(&specs, |&(i, scenario, seed)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?}", shape.name);
        let mut st = Stats::default();
        let s = state(shape, scenario, seed);
        let Ok((mut it, _)) = convert(&s) else {
            st.failures.push(format!("{tag}: conversion failed"));
            return st;
        };
        let mut rng = SplitMix64::new(seed ^ 0xbadc0de);
        let first = it.complex().n_vertices();
        let batches = 3;
        let per_batch = (shape.mesh.n_faces() / 4).clamp(10, 200);
        for b in 0..batches {
            for _ in 0..per_batch {
                let (f, bary, _) = pick_point(&it, &mut rng);
                match it.insert_vertex(f, bary) {
                    Ok(_) => st.inserted += 1,
                    Err(IntegerError::Intrinsic(_)) => st.refused += 1,
                    Err(err) => {
                        st.failures.push(format!("{tag}: insertion: {err}"));
                        return st;
                    }
                }
            }
            if let Err(err) = it.certify() {
                st.failures.push(format!("{tag}: batch {b}: {err}"));
                return st;
            }
            match it.flip_to_delaunay_certified(&DelaunayOptions::default()) {
                Ok((r, _)) => st.delaunay_flips += r.n_flips,
                Err(err) => {
                    st.failures
                        .push(format!("{tag}: batch {b}: Delaunay: {err}"));
                    return st;
                }
            }
        }
        // Remove the inserted vertices, last first (no vertex is renumbered then).
        let mut v = it.complex().n_vertices();
        while v > first {
            v -= 1;
            match it.remove_vertex(v) {
                Ok(_) => st.removed += 1,
                Err(IntegerError::RemovalStuck { .. }) => st.stuck += 1,
                Err(IntegerError::Unsupported { .. }) => st.boundary += 1,
                Err(err) => {
                    st.failures.push(format!("{tag}: removal of {v}: {err}"));
                    return st;
                }
            }
        }
        if let Err(err) = it.certify() {
            st.failures.push(format!("{tag}: after removals: {err}"));
        }
        st
    });
    let mut t = Stats::default();
    for r in results {
        t.inserted += r.inserted;
        t.refused += r.refused;
        t.delaunay_flips += r.delaunay_flips;
        t.removed += r.removed;
        t.stuck += r.stuck;
        t.boundary += r.boundary;
        t.failures.extend(r.failures);
    }
    eprintln!(
        "many insertions: {} inserted ({} refused by plan_face_split), {} Delaunay flips, then {} \
         removed, {} stuck (RemovalStuck), {} on the boundary; every batch certified; {} failures",
        t.inserted,
        t.refused,
        t.delaunay_flips,
        t.removed,
        t.stuck,
        t.boundary,
        t.failures.len()
    );
    for m in &t.failures {
        eprintln!("  {m}");
    }
    assert!(t.failures.is_empty());
    assert!(t.inserted > 5000 && t.removed > 4000);
}

// ----- App. A -----

/// Eq. 15–16 reproduce any affine map exactly (to rounding), for every number of corners.
#[test]
fn recover_barycentric_is_exact_for_affine_maps() {
    let mut rng = SplitMix64::new(15);
    let mut worst = 0.0f64;
    for trial in 0..2000 {
        let rho = 3 + trial % 4;
        // A random affine change: columns are the input coordinates of the corners of T¹'s face.
        let mut a = [[0.0; 3]; 3];
        for col in 0..3 {
            let b = inner_bary(&mut rng);
            for (row, x) in a.iter_mut().enumerate() {
                x[col] = b[row];
            }
        }
        let map = |v: [f64; 3]| [0, 1, 2].map(|r| (0..3).map(|c| a[r][c] * v[c]).sum::<f64>());
        let corners: Vec<[f64; 3]> = (0..rho).map(|_| inner_bary(&mut rng)).collect();
        let on_input: Vec<[f64; 3]> = corners.iter().map(|&v| map(v)).collect();
        let v = inner_bary(&mut rng);
        let (u, residual) = recover_barycentric(&corners, &on_input, v).unwrap();
        let want = map(v);
        let err = (0..3).map(|c| (u[c] - want[c]).abs()).fold(0.0, f64::max);
        worst = worst.max(err);
        assert!(residual < 1e-12, "residual {residual}");
    }
    eprintln!("App. A on random affine maps (3 to 6 corners): worst error {worst:e}");
    assert!(worst < 1e-10);
    // Fewer than three corners, or corners on a line: no answer.
    assert!(
        recover_barycentric(
            &[[1.0, 0.0, 0.0]; 2],
            &[[1.0, 0.0, 0.0]; 2],
            [1.0, 0.0, 0.0]
        )
        .is_none()
    );
    let line = [[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 1.0, 0.0]];
    assert!(recover_barycentric(&line, &line, [0.3, 0.3, 0.4]).is_none());
}

/// On a flat mesh every place is a point of the plane: the inserted vertex must be where the
/// barycentric combination of the places of its face's corners is (App. A reproduces the point
/// known in both, SPEC §9.2 `ic-recover-barycentric`), whatever the region, with every corner
/// of the region used (I16). On the input itself (`T¹ = T⁰`) the place is the point asked for.
#[test]
fn recover_reproduces_known_points() {
    let shapes: Vec<Shape> = standard_shapes().into_iter().filter(|s| s.flat).collect();
    assert!(shapes.len() >= 4);
    let mut worst_flat = 0.0f64;
    let mut worst_identity = 0.0f64;
    let mut n = 0;
    let mut hist = [0usize; 7];
    for shape in &shapes {
        // T¹ = T⁰: the place is the point itself.
        let input = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let mut rng = SplitMix64::new(3);
        for _ in 0..50 {
            let f = (rng.next_u64() % input.complex().n_faces() as u64) as usize;
            let bary = inner_bary(&mut rng);
            let x = input.clone().insert_vertex(f, bary).unwrap();
            let SurfacePoint::Face { face, bary: q } = x.location else {
                panic!("{x:?}");
            };
            assert_eq!(face as usize, f);
            worst_identity =
                worst_identity.max((0..3).map(|c| (q[c] - bary[c]).abs()).fold(0.0, f64::max));
        }
        for scenario in SCENARIOS {
            let s = state(shape, scenario, 1);
            let (mut it, _) = convert(&s).unwrap();
            let mut rng = SplitMix64::new(4);
            for _ in 0..60 {
                let (f, bary, _) = pick_point(&it, &mut rng);
                let corners = it.complex().faces()[f];
                let corner_pos = corners.map(|v| position(shape, it.vertex_location(v as usize)));
                let want: [f64; 3] = [0, 1, 2].map(|c| {
                    bary[0] * corner_pos[0][c]
                        + bary[1] * corner_pos[1][c]
                        + bary[2] * corner_pos[2][c]
                });
                let Ok(r) = it.split_face(f, bary, &FaceSplitOptions::default()) else {
                    continue;
                };
                n += 1;
                hist[r.recover.n_corners] += 1;
                let got = position(shape, r.insertion.location);
                worst_flat = worst_flat.max(common::dist(got, want) / shape.mean_edge_length());
            }
            it.certify().unwrap();
        }
    }
    eprintln!(
        "App. A on flat meshes: {n} insertions, worst distance from the known point {worst_flat:e} \
         mean edge lengths, regions with 3..6 corners {:?}; on T¹ = T⁰ worst coordinate error \
         {worst_identity:e}",
        &hist[3..]
    );
    assert!(n > 500);
    assert!(worst_flat < 1e-12);
    assert!(worst_identity < 1e-15);
    assert!(hist[3] > 0 && hist[4] > 0 && hist[5] + hist[6] > 0);
}

/// The point query `T¹ → T⁰` (SPEC §4.9) against the signposts' `point_on_input`, for points
/// of faces, of edges and vertices; and the place of a vertex inserted at a point is the
/// query's answer.
#[test]
fn point_queries_agree_with_the_signposts() {
    let shapes = all_shapes();
    let specs = specs(&shapes);
    let results = par_map(&specs, |&(i, scenario, seed)| {
        let shape = &shapes[i];
        let s = state(shape, scenario, seed);
        let (it, _) = convert(&s).unwrap();
        let mut rng = SplitMix64::new(seed ^ 0x9e77);
        let mut worst = (0.0f64, String::new());
        let mut n = 0usize;
        let mut failures = Vec::new();
        for q in 0..40 {
            let point = match q % 4 {
                0 | 1 => {
                    let f = (rng.next_u64() % it.complex().n_faces() as u64) as usize;
                    SurfacePoint::Face {
                        face: f as u32,
                        bary: inner_bary(&mut rng),
                    }
                }
                2 => SurfacePoint::Edge {
                    edge: (rng.next_u64() % it.complex().n_edges() as u64) as u32,
                    t: rng.uniform(0.05, 0.95),
                },
                _ => {
                    SurfacePoint::Vertex((rng.next_u64() % it.complex().n_vertices() as u64) as u32)
                }
            };
            let (Ok(a), b) = (s.point_on_input(point), it.point_on_input(point)) else {
                continue;
            };
            match b {
                Ok(b) => {
                    n += 1;
                    let d = place_distance(shape, a, b);
                    if d > worst.0 {
                        worst = (
                            d,
                            format!("{} / {scenario:?}: {point}: {a} vs {b}", shape.name),
                        );
                    }
                    if let SurfacePoint::Face { face, bary } = point {
                        let mut it2 = it.clone();
                        if let Ok(x) = it2.insert_vertex(face as usize, bary)
                            && bary_difference(x.location, b) > 1e-14
                        {
                            failures
                                .push(format!("{point}: inserted at {} vs query {b}", x.location));
                        }
                    }
                }
                Err(err) => failures.push(format!("{} / {scenario:?}: {point}: {err}", shape.name)),
            }
        }
        (n, worst, failures)
    });
    let mut n = 0;
    let mut worst = (0.0f64, String::new());
    let mut failures = Vec::new();
    for (k, w, f) in results {
        n += k;
        if w.0 > worst.0 {
            worst = w;
        }
        failures.extend(f);
    }
    eprintln!(
        "point queries: {n} agree with the signposts' point_on_input within {:e} mean edge \
         lengths ({}); {} failures",
        worst.0,
        worst.1,
        failures.len()
    );
    for m in &failures {
        eprintln!("  {m}");
    }
    assert!(failures.is_empty());
    assert!(n > 5000);
    assert!(worst.0 < 1e-8);
}

// ----- I6 -----

/// SPEC I6: the distance to the nearest arc is reported; a limit above it refuses the split and
/// changes nothing, a limit below it does not; points on an arc (the middle of a chord) still give
/// a valid state.
#[test]
fn too_close_to_a_curve_is_refused_and_points_on_arcs_are_valid() {
    let shapes = standard_shapes();
    let mut refused = 0;
    let mut accepted = 0;
    let mut on_arc = 0;
    let mut on_arc_differs = 0;
    let mut on_arc_snapped = 0;
    let mut on_arc_worst_margin = 0.0f64;
    for shape in &shapes {
        let s = state(shape, Scenario::RandomFlips, 2);
        let (it, _) = convert(&s).unwrap();
        let mut rng = SplitMix64::new(6);
        for _ in 0..40 {
            let (f, bary, _) = pick_point(&it, &mut rng);
            let Ok(loc) = it.locate_in_face(f, bary) else {
                continue;
            };
            let Some(d) = loc.choice.curve_distance else {
                continue;
            };
            let mut a = it.clone();
            match a.split_face(f, bary, &FaceSplitOptions { too_close: d * 1.5 }) {
                Err(IntegerError::TooCloseToCurve { distance, .. }) => {
                    assert_eq!(distance, d);
                    assert!(a == it, "a refused split changed the state");
                    refused += 1;
                }
                Err(IntegerError::Intrinsic(_)) => continue,
                other => panic!("{}: expected TooCloseToCurve, got {other:?}", shape.name),
            }
            let mut b = it.clone();
            if b.split_face(f, bary, &FaceSplitOptions { too_close: d * 0.5 })
                .is_ok()
            {
                accepted += 1;
            }
        }
        // Points on an arc: the midpoint of the two crossings of an arc around a corner.
        let mut s2 = s.clone();
        let mut it2 = it.clone();
        for f in 0..it.complex().n_faces() {
            if on_arc >= 400 {
                break;
            }
            let Some(a) = FaceArrangement::of(it.face_normal(f)) else {
                continue;
            };
            let Some(k) = (0..3).find(|&k| a.coords.c[k] > 0) else {
                continue;
            };
            // Arc 0 around corner k: crossing 0 of side k and the last crossing of side k + 2,
            // read from the signposts' traces of the two sides.
            let side = |slot: usize| {
                let h = 3 * f + slot;
                let e = it.complex().edge(h);
                let tr = s.trace_edge(e, &Default::default()).unwrap();
                let l = s.edge_lengths()[e];
                let forward = it.complex().edge_halfedge(e) == h;
                let params: Vec<f64> = tr.crossings.iter().map(|x| x.distance / l).collect();
                if forward {
                    params
                } else {
                    params.iter().rev().map(|p| 1.0 - p).collect()
                }
            };
            let (pk, pk2) = (side(k), side((k + 2) % 3));
            let (Some(&u), Some(&w)) = (pk.first(), pk2.last()) else {
                continue;
            };
            let mut p = [0.0; 3];
            p[k] += 0.5 * (1.0 - u) + 0.5 * w;
            p[(k + 1) % 3] += 0.5 * u;
            p[(k + 2) % 3] += 0.5 * (1.0 - w);
            let (rs, ri) = (
                s2.clone().insert_vertex(f, p),
                it2.split_face(f, p, &FaceSplitOptions::default()),
            );
            let Ok(ri) = ri else { continue };
            on_arc += 1;
            on_arc_worst_margin = on_arc_worst_margin.max(ri.location.choice.margin.abs());
            it2.certify()
                .unwrap_or_else(|e| panic!("{}: on an arc: {e}", shape.name));
            if let Ok(xs) = rs {
                if xs.snapped.is_some() {
                    on_arc_snapped += 1;
                } else {
                    let mut s3 = s2.clone();
                    s3.insert_vertex(f, p).unwrap();
                    if let Ok((c, _)) = convert(&s3)
                        && differs(&it2, &c).is_some()
                    {
                        on_arc_differs += 1;
                    }
                }
            }
            // Back to the state without the vertex.
            it2 = it.clone();
            s2 = s.clone();
        }
    }
    eprintln!(
        "I6: {refused} splits refused with a limit above the reported distance (state unchanged), \
         {accepted} accepted below it; {on_arc} insertions at the middle of an arc certify \
         (|margin| up to {on_arc_worst_margin:e} shortest sides); the signposts snap \
         {on_arc_snapped} of them onto the input edge, {on_arc_differs} of the others differ from \
         the signposts' conversion"
    );
    assert!(refused > 200 && accepted == refused);
    assert!(on_arc > 100);
}

// ----- moves -----

#[derive(Default)]
struct MoveStats {
    compared: usize,
    worst_place: f64,
    worst_lengths: f64,
    normal_multiset_differs: Vec<String>,
    beyond_star: usize,
    refused_stuck: usize,
    refused_other: std::collections::BTreeMap<String, usize>,
    signposts_refused: usize,
    ties: usize,
    restore_failures: Vec<String>,
    failures: Vec<String>,
}

/// A short name for the kind of an error.
fn kind_of(err: &IntegerError) -> String {
    match err {
        IntegerError::Intrinsic(IntrinsicError::MoveRefused { reason, .. }) => {
            format!("MoveRefused {reason:?}")
                .split([' ', '{', '('])
                .take(2)
                .collect::<Vec<_>>()
                .join(" ")
        }
        IntegerError::Intrinsic(IntrinsicError::MoveNotLocated { .. }) => "MoveNotLocated".into(),
        IntegerError::Intrinsic(IntrinsicError::InsertionRefused { reason, .. }) => {
            format!("InsertionRefused {reason:?}")
                .split([' ', '{', '('])
                .take(2)
                .collect::<Vec<_>>()
                .join(" ")
        }
        IntegerError::Intrinsic(other) => format!("{other:?}")
            .split([' ', '{', '('])
            .next()
            .unwrap_or("?")
            .to_string(),
        other => format!("{other:?}")
            .split([' ', '{', '('])
            .next()
            .unwrap_or("?")
            .to_string(),
    }
}

/// Moves of inserted vertices: inside their star against the signposts' local
/// `move_vertex` (same new place; after Delaunay flipping both, the same lengths and the same
/// normal coordinates as multisets), beyond their star (the signposts refuse; the integers
/// move, certify); every refused move restores the state exactly.
#[test]
fn moves_match_the_signposts_and_refusals_restore() {
    let shapes = all_shapes();
    let specs: Vec<(usize, Scenario, u64)> = (0..shapes.len())
        .flat_map(|i| {
            [Scenario::Insertions, Scenario::Mixed, Scenario::EdgeSplits]
                .into_iter()
                .flat_map(move |sc| [3u64, 8].map(move |seed| (i, sc, seed)))
        })
        .collect();
    let results = par_map(&specs, |&(i, scenario, seed)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?} / {seed}", shape.name);
        let mut st = MoveStats::default();
        let s = state(shape, scenario, seed);
        let (it, _) = convert(&s).unwrap();
        let mut rng = SplitMix64::new(seed ^ 0x3047e);
        let inserted: Vec<usize> = (it.n_input_vertices()..it.complex().n_vertices())
            .filter(|&v| !it.complex().is_boundary_vertex(v))
            .collect();
        for &v in inserted
            .iter()
            .take(if shape.mesh.n_edges() > 5000 { 8 } else { 20 })
        {
            let Some(min_shape) = s.star_min_shape(v) else {
                continue;
            };
            let _ = min_shape;
            // A displacement in the signposts' frame of v, a fraction of the shortest edge at v.
            let shortest = s
                .complex()
                .outgoing_halfedges(v)
                .map(|h| s.edge_lengths()[s.complex().edge(h)])
                .fold(f64::INFINITY, f64::min);
            let psi = rng.uniform(0.0, std::f64::consts::TAU);
            let r = rng.uniform(0.05, 0.3) * shortest;
            // The same direction in the integer frame (angle frames, SPEC I9).
            let out = it.complex().outgoing_halfedge(v);
            let phi = psi - s.signpost(out);
            let d_s = [r * psi.cos(), r * psi.sin()];
            let d_i = [r * phi.cos(), r * phi.sin()];
            let (mut s2, mut it2) = (s.clone(), it.clone());
            let before = it2.clone();
            match (
                s2.move_vertex(v, d_s),
                it2.move_vertex(v, d_i, &FaceSplitOptions::default()),
            ) {
                (Ok(ms), Ok(mi)) => {
                    st.compared += 1;
                    st.worst_place = st.worst_place.max(place_distance(shape, ms.to, mi.to));
                    if let Err(err) = it2.certify() {
                        st.failures
                            .push(format!("{tag}: move of {v}: certify: {err}"));
                        continue;
                    }
                    // The same Delaunay triangulation of the same points.
                    let ds = s2.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
                    let di = it2.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
                    let mut ls: Vec<f64> = s2.edge_lengths().to_vec();
                    let mut li: Vec<f64> = it2.edge_lengths().to_vec();
                    ls.sort_by(f64::total_cmp);
                    li.sort_by(f64::total_cmp);
                    if ls.len() != li.len() {
                        st.failures
                            .push(format!("{tag}: move of {v}: edge counts differ"));
                        continue;
                    }
                    let dl = ls
                        .iter()
                        .zip(&li)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0, f64::max)
                        / shape.mean_edge_length();
                    st.worst_lengths = st.worst_lengths.max(dl);
                    if let Ok((c, _)) = convert(&s2) {
                        let mut ns = c.normal_coordinates().to_vec();
                        let mut ni = it2.normal_coordinates().to_vec();
                        ns.sort();
                        ni.sort();
                        if ns != ni && (ds.n_ties > 0 || di.n_ties > 0) {
                            // Not unique: a tie of the Delaunay criterion (cocircular corners)
                            // was broken differently.
                            st.ties += 1;
                        } else if ns != ni {
                            let how = match &mi.insertion {
                                Some(MoveInsertion::Face(r)) => format!(
                                    "face split, arc at {:?}",
                                    r.location.choice.curve_distance
                                ),
                                Some(MoveInsertion::Edge(_)) => "edge split".to_string(),
                                None => "no move".to_string(),
                            };
                            st.normal_multiset_differs.push(format!(
                                "{tag}: move of {v} ({how}; signposts' place {} snapped {:?}, \
                                 integer {}): Σn {} vs {}",
                                ms.to,
                                ms.snapped,
                                mi.to,
                                ns.iter().map(|&x| x.max(0)).sum::<i64>(),
                                ni.iter().map(|&x| x.max(0)).sum::<i64>()
                            ));
                        }
                    }
                }
                (_, Err(err)) => {
                    if it2 != before {
                        st.restore_failures
                            .push(format!("{tag}: move of {v}: {err}"));
                    }
                    match err {
                        IntegerError::RemovalStuck { .. } => st.refused_stuck += 1,
                        other => *st.refused_other.entry(kind_of(&other)).or_default() += 1,
                    }
                }
                (Err(_), Ok(_)) => {
                    // The star is thin and the move leaves it (NotInStar, TooThin): the
                    // signposts' local update does not apply; the integer move must be valid.
                    st.signposts_refused += 1;
                    if let Err(err) = it2.certify() {
                        st.failures
                            .push(format!("{tag}: move of {v}: certify: {err}"));
                    }
                }
            }
            // Beyond the star: three times the longest edge at v.
            let longest = it
                .complex()
                .outgoing_halfedges(v)
                .map(|h| it.edge_lengths()[it.complex().edge(h)])
                .fold(0.0, f64::max);
            let far = 3.0 * longest;
            let mut it3 = it.clone();
            match it3.move_vertex(
                v,
                [far * phi.cos(), far * phi.sin()],
                &FaceSplitOptions::default(),
            ) {
                Ok(_) => {
                    st.beyond_star += 1;
                    if let Err(err) = it3.certify() {
                        st.failures
                            .push(format!("{tag}: far move of {v}: certify: {err}"));
                    }
                }
                Err(err) => {
                    if it3 != it {
                        st.restore_failures
                            .push(format!("{tag}: far move of {v}: {err}"));
                    }
                    match err {
                        IntegerError::RemovalStuck { .. } => st.refused_stuck += 1,
                        other => *st.refused_other.entry(kind_of(&other)).or_default() += 1,
                    }
                }
            }
        }
        st
    });
    let mut t = MoveStats::default();
    for r in results {
        t.compared += r.compared;
        t.worst_place = t.worst_place.max(r.worst_place);
        t.worst_lengths = t.worst_lengths.max(r.worst_lengths);
        t.normal_multiset_differs.extend(r.normal_multiset_differs);
        t.beyond_star += r.beyond_star;
        t.refused_stuck += r.refused_stuck;
        for (k, n) in r.refused_other {
            *t.refused_other.entry(k).or_default() += n;
        }
        t.signposts_refused += r.signposts_refused;
        t.ties += r.ties;
        t.restore_failures.extend(r.restore_failures);
        t.failures.extend(r.failures);
    }
    eprintln!(
        "moves: {} inside the star compared with the signposts' move_vertex: places within {:e} \
         mean edge lengths, after Delaunay flipping both lengths within {:e}, normal coordinates \
         as multisets differ in {} ({} more after Delaunay ties); {} moves in thin stars the signposts refused; {} moves beyond \
         the star (signposts cannot); refused: {} RemovalStuck, others {:?}, all restored ({} \
         not); {} failures",
        t.compared,
        t.worst_place,
        t.worst_lengths,
        t.normal_multiset_differs.len(),
        t.ties,
        t.signposts_refused,
        t.beyond_star,
        t.refused_stuck,
        t.refused_other,
        t.restore_failures.len(),
        t.failures.len()
    );
    for m in t.normal_multiset_differs.iter().take(10) {
        eprintln!("  n multiset differs: {m}");
    }
    for m in t.restore_failures.iter().chain(&t.failures) {
        eprintln!("  {m}");
    }
    assert!(t.restore_failures.is_empty() && t.failures.is_empty());
    assert!(t.normal_multiset_differs.is_empty());
    assert!(t.compared > 300);
    assert!(t.worst_place < 1e-8);
    assert!(t.worst_lengths < 1e-8);
}

/// Refused moves: an input vertex, a vertex on the boundary, a displacement that is not finite;
/// a removal that gets stuck (the straight star of an edge-split vertex, I27) restores the
/// whole state.
#[test]
fn refused_moves_change_nothing() {
    let shape = &standard_shapes()[0];
    let mut it = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    let before = it.clone();
    assert!(matches!(
        it.move_vertex(0, [0.01, 0.0], &FaceSplitOptions::default()),
        Err(IntegerError::Intrinsic(IntrinsicError::MoveRefused {
            reason: MoveRefusal::InputVertex,
            ..
        }))
    ));
    it.insert_vertex(0, [0.3, 0.3, 0.4]).unwrap();
    let v = it.complex().n_vertices() - 1;
    let after_insert = it.clone();
    assert!(matches!(
        it.move_vertex(v, [f64::NAN, 0.0], &FaceSplitOptions::default()),
        Err(IntegerError::Intrinsic(IntrinsicError::MoveRefused {
            reason: MoveRefusal::NotFinite,
            ..
        }))
    ));
    assert!(it == after_insert);
    let zero = it
        .move_vertex(v, [0.0, 0.0], &FaceSplitOptions::default())
        .unwrap();
    assert!(zero.insertion.is_none() && it == after_insert);
    let _ = before;
    // A boundary vertex.
    let grid = standard_shapes()
        .into_iter()
        .find(|s| s.name.starts_with("grid"))
        .unwrap();
    let mut g = IntegerTriangulation::new(&grid.mesh, &grid.geometry).unwrap();
    let e = (0..g.complex().n_edges())
        .find(|&e| g.complex().twin(g.complex().edge_halfedge(e)).is_none())
        .unwrap();
    g.split_edge(e, 0.5).unwrap();
    let b = g.complex().n_vertices() - 1;
    let gb = g.clone();
    assert!(matches!(
        g.move_vertex(b, [0.01, 0.0], &FaceSplitOptions::default()),
        Err(IntegerError::Unsupported { .. })
    ));
    assert!(g == gb);
    // Stuck removals: vertices made by splitting interior edges (straight stars, I27), moved a
    // little; whenever the removal is stuck the state must be exactly as before.
    let mut stuck = 0;
    let mut moved = 0;
    for shape in standard_shapes() {
        let s = state(&shape, Scenario::Delaunay, 4);
        let (mut it, _) = convert(&s).unwrap();
        let mut rng = SplitMix64::new(12);
        for _ in 0..30 {
            let e = (rng.next_u64() % it.complex().n_edges() as u64) as usize;
            if it.complex().twin(it.complex().edge_halfedge(e)).is_none() {
                continue;
            }
            let mut t = it.clone();
            if t.split_edge(e, rng.uniform(0.3, 0.7)).is_err() {
                continue;
            }
            let v = t.complex().n_vertices() - 1;
            let before = t.clone();
            let h = t.complex().outgoing_halfedge(v);
            let l = t.edge_lengths()[t.complex().edge(h)];
            let phi = rng.uniform(0.0, std::f64::consts::TAU);
            let r = 0.1 * l;
            match t.move_vertex(
                v,
                [r * phi.cos(), r * phi.sin()],
                &FaceSplitOptions::default(),
            ) {
                Ok(m) => {
                    moved += 1;
                    assert_eq!(m.vertex as usize, v);
                    t.certify().unwrap();
                }
                Err(IntegerError::RemovalStuck { .. }) => {
                    stuck += 1;
                    assert!(
                        t == before,
                        "{}: a stuck move changed the state",
                        shape.name
                    );
                }
                Err(err) => {
                    assert!(
                        t == before,
                        "{}: a refused move changed the state: {err}",
                        shape.name
                    );
                }
            }
            it = t;
        }
    }
    eprintln!(
        "moves of edge-split vertices: {moved} done, {stuck} refused as RemovalStuck (state restored)"
    );
    assert!(stuck > 0 && moved > 0);
}

// ----- the shared interface -----

#[test]
fn the_interface_inserts() {
    let shape = &standard_shapes()[1];
    let s = state(shape, Scenario::RandomFlips, 1);
    let (mut it, _) = convert(&s).unwrap();
    let mut s2 = s.clone();
    let a = Triangulation::insert_vertex(&mut it, 3, [0.2, 0.3, 0.5]).unwrap();
    let b = s2.insert_vertex(3, [0.2, 0.3, 0.5]).unwrap();
    assert_eq!(
        (a.vertex, &a.new_edges, &a.new_faces),
        (b.vertex, &b.new_edges, &b.new_faces)
    );
    it.certify().unwrap();
    // A refusal of the shared kind comes through as itself.
    let before = it.clone();
    assert!(matches!(
        Triangulation::insert_vertex(&mut it, 0, [1.0 - 2e-9, 1e-9, 1e-9]),
        Err(IntrinsicError::InsertionRefused { .. })
    ));
    assert!(it == before);
}

// ----- timing -----

/// Insertion vs the signposts' insertion (single thread, release, run alone):
/// `cargo test --release -p ddg-integer-coords --test insert -- --ignored --nocapture timing`.
#[test]
#[ignore]
fn timing() {
    let shapes: Vec<Shape> = all_shapes()
        .into_iter()
        .filter(|s| s.name.starts_with("spot") || s.name == "camelbad")
        .collect();
    for shape in &shapes {
        let mut s = shape.triangulation();
        s.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let (it, _) = convert(&s).unwrap();
        let mut rng = SplitMix64::new(1);
        let n = 2000;
        let points: Vec<(usize, [f64; 3])> = (0..n)
            .map(|_| {
                (
                    (rng.next_u64() % s.n_faces() as u64) as usize,
                    inner_bary(&mut rng),
                )
            })
            .collect();
        let mut best = [f64::INFINITY; 3];
        let mut crossings = 0i64;
        for round in 0..5 {
            let mut s2 = s.clone();
            let t0 = std::time::Instant::now();
            for &(f, b) in &points {
                let _ = s2.insert_vertex(f, b);
            }
            best[0] = best[0].min(t0.elapsed().as_secs_f64() / n as f64);
            let mut i2 = it.clone();
            let t0 = std::time::Instant::now();
            for &(f, b) in &points {
                let _ = i2.insert_vertex(f, b);
            }
            best[1] = best[1].min(t0.elapsed().as_secs_f64() / n as f64);
            let t0 = std::time::Instant::now();
            for &(f, b) in &points {
                let _ = it.locate_in_face(f, b);
            }
            best[2] = best[2].min(t0.elapsed().as_secs_f64() / n as f64);
            if round == 0 {
                crossings = points.iter().map(|&(f, _)| face_crossings(&it, f)).sum();
            }
        }
        eprintln!(
            "{}: |E| {}, {n} insertions after Delaunay flipping (mean {:.1} crossings per face): \
             signposts {:.2} µs, integer {:.2} µs, point query {:.2} µs",
            shape.name,
            s.n_edges(),
            crossings as f64 / n as f64,
            best[0] * 1e6,
            best[1] * 1e6,
            best[2] * 1e6
        );
    }
}
