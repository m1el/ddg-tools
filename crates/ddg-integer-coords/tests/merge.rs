//! The merge of the two halves of a split edge (SPEC §4.6, I27): the integers of the merged edge
//! (`ic-merge-edge-normal`), split then merge as the identity, and removal of every split vertex
//! after Delaunay flipping, differentially against the signposts, back to the intrinsic Delaunay
//! triangulation.
//!
//! Run with `--nocapture` to see the counts.

mod common;

use std::collections::BTreeMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use common::{SCENARIOS, Scenario, all_shapes, convert, standard_shapes, state};
use ddg_core::generators::SplitMix64;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::{IntegerError, IntegerTriangulation, Removal};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, IntrinsicError, IntrinsicTriangulation, RemovalStep,
    SurfacePoint, Triangulation,
};

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

/// The first difference between two integer states, the lengths of the edges in `loose`
/// within `ulps` units in the last place and all others bitwise.
fn differs(
    a: &IntegerTriangulation,
    b: &IntegerTriangulation,
    loose: &[usize],
    ulps: u64,
) -> Option<String> {
    if a.complex() != b.complex() {
        return Some("connectivity differs".into());
    }
    for (e, (x, y)) in a.edge_lengths().iter().zip(b.edge_lengths()).enumerate() {
        let d = x.to_bits().abs_diff(y.to_bits());
        if d > if loose.contains(&e) { ulps } else { 0 } {
            return Some(format!("length of edge {e}: {x:e} vs {y:e}"));
        }
    }
    if a.normal_coordinates() != b.normal_coordinates() {
        return Some("normal coordinates differ".into());
    }
    if a.roundabouts() != b.roundabouts() {
        return Some("roundabouts differ".into());
    }
    if a.vertex_locations() != b.vertex_locations() {
        return Some("places differ".into());
    }
    None
}

/// The first difference between the integer state and the conversion of the signposts.
fn difference(it: &IntegerTriangulation, s: &IntrinsicTriangulation) -> Option<String> {
    match convert(s) {
        Ok((c, _)) => differs(it, &c, &[], 0),
        Err(e) => Some(format!("conversion failed: {e}")),
    }
}

/// Splits edge `e` at `t` in both representations, with the signposts' crossings and place
/// for a crossed edge so that the states stay equal (as `tests/ops.rs`). `None` if either
/// refuses.
fn split_both(
    s: &mut IntrinsicTriangulation,
    it: &mut IntegerTriangulation,
    e: usize,
    t: f64,
) -> Option<usize> {
    let n = it.normal_coordinate(e);
    if n == -1 && s.edge_on_input(e).is_none() {
        return None;
    }
    if n == -1 {
        let (sb, ib) = (s.clone(), it.clone());
        match (s.split_edge(e, t), it.split_edge(e, t)) {
            (Ok(x), Ok(_)) => Some(x.vertex as usize),
            _ => {
                *s = sb;
                *it = ib;
                None
            }
        }
    } else {
        let trace = s.trace_edge(e, &EdgeTraceOptions::default()).ok()?;
        let at = t * s.edge_lengths()[e];
        let before = trace.crossings.iter().filter(|x| x.distance < at).count() as i64;
        let near = trace
            .crossings
            .iter()
            .any(|x| (x.distance - at).abs() < 1e-9 * s.edge_lengths()[e]);
        if near {
            return None;
        }
        let x = s.split_edge(e, t).ok()?;
        it.split_edge_given(e, t, before, x.location).unwrap();
        Some(x.vertex as usize)
    }
}

/// Does in the signposts what the integer removal `r` of vertex `v` did.
fn replay_removal(s: &mut IntrinsicTriangulation, v: usize, r: &Removal) -> Result<(), String> {
    for &e in &r.flips {
        s.flip_edge(e as usize)
            .map_err(|err| format!("signposts refuse flip {e}: {err}"))?;
    }
    match &r.removal {
        RemovalStep::Star(x) => match s.remove_vertex(v) {
            Ok(y) if &y == x => Ok(()),
            other => Err(format!("signposts' removal of {v}: {other:?}")),
        },
        RemovalStep::Merge(m) => match s.merge_edge(v, m.halves[0] as usize) {
            Ok(y) if &y == m => Ok(()),
            other => Err(format!("signposts' merge at {v}: {other:?}")),
        },
    }
}

/// A split followed at once by the merge gives back the integer state: connectivity, every
/// `n`, every roundabout and every place exactly, every length bitwise except the merged one
/// (`t ℓ + (1 − t) ℓ`, within one ulp of `ℓ`): along input edges (Eq. 6), across them
/// (`ic-split-crossed-edge`), on the boundary, in every scenario.
#[test]
fn merge_right_after_a_split_is_the_identity() {
    let shapes = all_shapes();
    let specs: Vec<(usize, Scenario)> = (0..shapes.len())
        .flat_map(|i| SCENARIOS.into_iter().map(move |sc| (i, sc)))
        .collect();
    #[derive(Default)]
    struct Stats {
        // along an input edge (interior), along the boundary, across input edges
        kinds: [usize; 3],
        bitwise: usize,
        worst_ulps: u64,
        refused_split: usize,
        failures: Vec<String>,
    }
    let results = par_map(&specs, |&(i, scenario)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?}", shape.name);
        let mut st = Stats::default();
        let s = state(shape, scenario, 31);
        let Ok((mut it, _)) = convert(&s) else {
            st.failures.push(format!("{tag}: conversion failed"));
            return st;
        };
        let mut rng = SplitMix64::new(0x1227);
        let attempts = if it.complex().n_edges() > 5000 {
            60
        } else {
            150
        };
        for _ in 0..attempts {
            let e = (rng.next_u64() % it.complex().n_edges() as u64) as usize;
            let t = rng.uniform(0.05, 0.95);
            let before = it.clone();
            let n = it.normal_coordinate(e);
            let boundary = it.complex().is_boundary_edge(e);
            let Ok(split) = it.split_edge(e, t) else {
                st.refused_split += 1;
                continue;
            };
            let v = split.vertex as usize;
            let half = split.halves.unwrap()[(rng.next_u64() % 2) as usize] as usize;
            match it.merge_edge(v, half) {
                Ok(m) => {
                    assert_eq!(m.edge as usize, e);
                    if let Some(d) = differs(&it, &before, &[e], 1) {
                        st.failures
                            .push(format!("{tag}: split and merge of {e}: {d}"));
                        return st;
                    }
                    let d = it.edge_lengths()[e]
                        .to_bits()
                        .abs_diff(before.edge_lengths()[e].to_bits());
                    st.bitwise += usize::from(d == 0);
                    st.worst_ulps = st.worst_ulps.max(d);
                    st.kinds[if n >= 0 {
                        2
                    } else if boundary {
                        1
                    } else {
                        0
                    }] += 1;
                }
                Err(err) => {
                    st.failures
                        .push(format!("{tag}: merge after the split of {e}: {err}"));
                    return st;
                }
            }
            // Go on from the split state (deeper Δ-complex configurations).
            it = before;
            let _ = it.split_edge(e, t);
        }
        if let Err(err) = it.certify() {
            st.failures.push(format!("{tag}: certify: {err}"));
        }
        st
    });
    let mut t = Stats::default();
    for r in results {
        for k in 0..3 {
            t.kinds[k] += r.kinds[k];
        }
        t.bitwise += r.bitwise;
        t.worst_ulps = t.worst_ulps.max(r.worst_ulps);
        t.refused_split += r.refused_split;
        t.failures.extend(r.failures);
    }
    let total: usize = t.kinds.iter().sum();
    eprintln!(
        "split then merge: {total} ({} along interior input edges, {} along the boundary, {} \
         across input edges; {} splits refused), integers and connectivity exact, merged length \
         bitwise in {}, otherwise within {} ulp; {} failures",
        t.kinds[0],
        t.kinds[1],
        t.kinds[2],
        t.refused_split,
        t.bitwise,
        t.worst_ulps,
        t.failures.len()
    );
    for f in &t.failures {
        eprintln!("  {f}");
    }
    assert!(t.failures.is_empty());
    assert!(t.kinds.iter().all(|&k| k > 50));
    assert!(t.worst_ulps <= 1);
}

/// Many splits merged in the reverse order give back the state (lengths within one ulp).
#[test]
fn many_splits_merged_in_reverse_order_restore_the_state() {
    for shape in standard_shapes() {
        let s = state(&shape, Scenario::Mixed, 3);
        let (mut it, _) = convert(&s).unwrap();
        let original = it.clone();
        let mut rng = SplitMix64::new(77);
        let mut done = Vec::new();
        for _ in 0..it.complex().n_edges() {
            let e = (rng.next_u64() % it.complex().n_edges() as u64) as usize;
            if let Ok(x) = it.split_edge(e, rng.uniform(0.1, 0.9)) {
                done.push((x.vertex as usize, e));
            }
        }
        it.certify().unwrap();
        let touched: Vec<usize> = done.iter().map(|&(_, e)| e).collect();
        for &(v, e) in done.iter().rev() {
            it.merge_edge(v, e).unwrap();
        }
        if let Some(d) = differs(&it, &original, &touched, 4) {
            panic!("{}: {d}", shape.name);
        }
        it.certify().unwrap();
    }
}

/// A canonical form of a state without vertex insertions, independent of edge, face and
/// halfedge ids: per edge its end points (sorted), `n` and length; per halfedge leaving an
/// input vertex its end points and roundabout.
#[allow(clippy::type_complexity)]
fn canonical(
    t: &IntegerTriangulation,
) -> (Vec<(usize, usize, i64, f64)>, Vec<(usize, usize, u32)>) {
    let c = t.complex();
    let mut edges: Vec<(usize, usize, i64, f64)> = (0..c.n_edges())
        .map(|e| {
            let [a, b] = c.edge_vertices(e);
            (
                a.min(b),
                a.max(b),
                t.normal_coordinate(e),
                t.edge_lengths()[e],
            )
        })
        .collect();
    edges.sort_by(|x, y| {
        (x.0, x.1, x.2)
            .cmp(&(y.0, y.1, y.2))
            .then(x.3.total_cmp(&y.3))
    });
    let mut halfedges: Vec<(usize, usize, u32)> = (0..c.n_halfedges())
        .map(|h| (c.from_vertex(h), c.to_vertex(h), t.roundabouts()[h]))
        .collect();
    halfedges.sort_unstable();
    (edges, halfedges)
}

/// Splits many edges of the intrinsic Delaunay triangulation (along input edges, across them,
/// on the boundary), flips to Delaunay, removes every split vertex (last first) and compares:
/// after each removal the integers equal `from_intrinsic` of the signposts that did the same
/// flips and the same removal or merge, and certify; at the end, flipped to Delaunay again, the
/// state is the intrinsic Delaunay triangulation it started from (unique where no edge is
/// within `1e-9` of a Delaunay tie): connectivity, `n`, roundabouts, lengths within `1e-12`.
#[test]
fn split_delaunay_remove_every_split_vertex_gives_back_the_idt() {
    let shapes = all_shapes();
    let specs: Vec<(usize, u64)> = (0..shapes.len())
        .flat_map(|i| [5u64, 77].into_iter().map(move |seed| (i, seed)))
        .collect();
    #[derive(Default)]
    struct Stats {
        split: usize,
        removed: [usize; 2], // by flips (Alg. 8), by a merge
        stuck: BTreeMap<String, usize>,
        idt_equal: usize,
        idt_skipped_ties: usize,
        idt_skipped_stuck: usize,
        worst_length: f64,
        failures: Vec<String>,
    }
    let results = par_map(&specs, |&(i, seed)| {
        let shape = &shapes[i];
        let tag = format!("{} / {seed}", shape.name);
        let mut st = Stats::default();
        let opts = DelaunayOptions::default();
        let mut s = shape.triangulation();
        s.flip_to_delaunay(&opts).unwrap();
        let Ok((mut it, _)) = convert(&s) else {
            st.failures.push(format!("{tag}: conversion failed"));
            return st;
        };
        let idt = it.clone();
        let ties = (0..idt.complex().n_edges())
            .filter_map(|e| idt.edge_angle_excess(e))
            .any(|x| x > -1e-9);
        let small = it.complex().n_edges() < 3000;
        let mut rng = SplitMix64::new(seed ^ 0x5911);
        let n_input = it.n_input_vertices();
        let target = (it.complex().n_edges() / 5).clamp(10, 200);
        let mut tries = 0;
        while st.split < target && tries < 20 * target {
            tries += 1;
            let e = (rng.next_u64() % it.complex().n_edges() as u64) as usize;
            if split_both(&mut s, &mut it, e, rng.uniform(0.15, 0.85)).is_some() {
                st.split += 1;
            }
        }
        s.flip_to_delaunay(&opts).unwrap();
        Triangulation::flip_to_delaunay(&mut it, &opts).unwrap();
        if let Some(d) = difference(&it, &s) {
            st.failures.push(format!("{tag}: after the splits: {d}"));
            return st;
        }
        let mut v = it.complex().n_vertices();
        while v > n_input {
            v -= 1;
            let before = it.clone();
            match it.remove_vertex(v) {
                Ok(r) => {
                    st.removed[usize::from(matches!(r.removal, RemovalStep::Merge(_)))] += 1;
                    if let Err(m) = replay_removal(&mut s, v, &r) {
                        st.failures.push(format!("{tag}: {m}"));
                        return st;
                    }
                    if small {
                        if let Some(d) = difference(&it, &s) {
                            st.failures.push(format!("{tag}: after removing {v}: {d}"));
                            return st;
                        }
                        if let Err(err) = it.certify() {
                            st.failures
                                .push(format!("{tag}: after removing {v}: {err}"));
                            return st;
                        }
                    }
                }
                Err(err) => {
                    if it != before {
                        st.failures
                            .push(format!("{tag}: failed removal changed the state"));
                        return st;
                    }
                    let key = match err {
                        IntegerError::RemovalStuck { degree, .. } => {
                            let self_edge = it
                                .complex()
                                .outgoing_halfedges(v)
                                .any(|h| it.complex().to_vertex(h) == v);
                            format!(
                                "stuck at degree {degree}{}",
                                if self_edge { " (self-edge)" } else { "" }
                            )
                        }
                        IntegerError::Intrinsic(IntrinsicError::RemovalRefused {
                            reason, ..
                        }) => format!("{reason:?}")
                            .split([' ', '{'])
                            .next()
                            .unwrap()
                            .to_string(),
                        other => {
                            st.failures.push(format!("{tag}: removal of {v}: {other}"));
                            return st;
                        }
                    };
                    *st.stuck.entry(key).or_default() += 1;
                }
            }
        }
        if let Err(err) = it.certify() {
            st.failures
                .push(format!("{tag}: after the removals: {err}"));
            return st;
        }
        if let Some(d) = difference(&it, &s) {
            st.failures.push(format!("{tag}: after the removals: {d}"));
            return st;
        }
        if it.complex().n_vertices() > n_input {
            st.idt_skipped_stuck += 1;
            return st;
        }
        Triangulation::flip_to_delaunay(&mut it, &opts).unwrap();
        if ties {
            st.idt_skipped_ties += 1;
            return st;
        }
        let (ea, ha) = canonical(&it);
        let (eb, hb) = canonical(&idt);
        let same_ints = ea.len() == eb.len()
            && ea
                .iter()
                .zip(&eb)
                .all(|(x, y)| (x.0, x.1, x.2) == (y.0, y.1, y.2))
            && ha == hb;
        if !same_ints {
            st.failures.push(format!(
                "{tag}: not the intrinsic Delaunay triangulation again"
            ));
            return st;
        }
        for (x, y) in ea.iter().zip(&eb) {
            st.worst_length = st.worst_length.max((x.3 - y.3).abs() / y.3);
        }
        st.idt_equal += 1;
        st
    });
    let mut t = Stats::default();
    for r in results {
        t.split += r.split;
        t.removed[0] += r.removed[0];
        t.removed[1] += r.removed[1];
        for (k, v) in r.stuck {
            *t.stuck.entry(k).or_default() += v;
        }
        t.idt_equal += r.idt_equal;
        t.idt_skipped_ties += r.idt_skipped_ties;
        t.idt_skipped_stuck += r.idt_skipped_stuck;
        t.worst_length = t.worst_length.max(r.worst_length);
        t.failures.extend(r.failures);
    }
    eprintln!(
        "split, Delaunay, remove: {} vertices split, {} removed by flips (Alg. 8), {} by a merge, \
         not removed {:?}; every removal equal to the signposts' and certified; back to the iDT \
         in {} runs (lengths within {:e} relative), {} not compared (Delaunay ties), {} not \
         compared (a vertex left); {} failures",
        t.split,
        t.removed[0],
        t.removed[1],
        t.stuck,
        t.idt_equal,
        t.worst_length,
        t.idt_skipped_ties,
        t.idt_skipped_stuck,
        t.failures.len()
    );
    for f in &t.failures {
        eprintln!("  {f}");
    }
    assert!(t.failures.is_empty());
    assert!(t.idt_equal > 5 && t.removed[1] > 100);
    assert!(t.worst_length < 1e-12);
}

/// A vertex on an input edge whose two straight edges cross that input edge there: the merged
/// edge crosses it once at the vertex, `n = n₁ + n₂ + 1` (the third case of
/// `ic-merge-edge-normal`). Flat square `a c b d` with the diagonal `a b` an input edge; split
/// it where the other diagonal `c d` crosses it, merge the halves towards `c` and `d`: the
/// result is the square with the diagonal `c d` (n = 1), as the flip of `a b` gives.
#[test]
fn a_merge_across_an_input_edge_counts_the_crossing_at_the_vertex() {
    let mesh = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
    let pos = [
        [0.0, 0.0, 0.0],
        [1.0, -0.4, 0.0],
        [2.0, 0.0, 0.0],
        [0.8, 0.7, 0.0],
    ];
    let geometry = Geometry::from_positions(&mesh, &pos).unwrap();
    let mut it = IntegerTriangulation::new(&mesh, &geometry).unwrap();
    let mut s = IntrinsicTriangulation::new(&mesh, &geometry).unwrap();
    let diag = (0..it.complex().n_edges())
        .find(|&e| !it.complex().is_boundary_edge(e))
        .unwrap();
    // Where the segment 1–3 crosses 0–2 (y = 0): 1 + (0.4 / 1.1)(0.8 − 1) along x.
    let x = 1.0 + (0.4 / 1.1) * (0.8 - 1.0);
    let [a, _] = it.complex().edge_vertices(diag);
    let t = if a == 0 { x / 2.0 } else { 1.0 - x / 2.0 };
    let split = it.split_edge(diag, t).unwrap();
    s.split_edge(diag, t).unwrap();
    let v = split.vertex as usize;
    let to_opposite = split.new_edges[1] as usize;
    let mut flipped = IntegerTriangulation::new(&mesh, &geometry).unwrap();
    flipped.flip_edge(diag).unwrap();
    assert_eq!(flipped.normal_coordinate(diag), 1);
    // The edges to the opposite corners are straight by construction (within rounding).
    let plan =
        ddg_intrinsic::plan_edge_merge(it.complex(), it.edge_lengths(), v, to_opposite).unwrap();
    assert!(
        plan.excess.iter().all(|x| x.abs() < 1e-12),
        "{:?}",
        plan.excess
    );
    let m = it.merge_edge(v, to_opposite).unwrap();
    s.merge_edge(v, to_opposite).unwrap();
    assert_eq!(it.normal_coordinate(m.edge as usize), 1);
    it.certify().unwrap();
    assert!(matches!(it.vertex_locations()[0], SurfacePoint::Vertex(0)));
    assert_eq!(difference(&it, &s), None);
    let ends = it.complex().edge_vertices(m.edge as usize);
    let mut ends = ends.to_vec();
    ends.sort_unstable();
    assert_eq!(ends, vec![1, 3]);
    assert!((it.edge_lengths()[m.edge as usize] - flipped.edge_lengths()[diag]).abs() < 1e-14);
    // Merging the halves along the input edge instead gives the input back.
    let mut back = IntegerTriangulation::new(&mesh, &geometry).unwrap();
    let split = back.split_edge(diag, t).unwrap();
    back.merge_edge(split.vertex as usize, diag).unwrap();
    assert_eq!(back.normal_coordinate(diag), -1);
    back.certify().unwrap();
}
