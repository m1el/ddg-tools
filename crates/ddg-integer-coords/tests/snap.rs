//! The edge split at a crossing (SPEC I28, `ic-split-at-crossing`): the integer split of a point
//! within rounding of a crossing of the edge with an input edge puts the vertex on that input
//! edge, as the signposts do, and the integers equal those of `from_intrinsic` of the signposts
//! after the same split (connectivity, lengths bitwise, every `n`, every roundabout, the places
//! up to the position along the input edge). `certify` passes after every split.
//!
//! Run with `--nocapture` to see the counts.

mod common;

use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use common::{SCENARIOS, Scenario, all_shapes, convert, state};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{EdgeTraceOptions, IntrinsicTriangulation, SurfacePoint};

fn par_map<S: Sync, R: Send>(specs: &[S], f: impl Fn(&S) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(20);
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

/// The first difference between two integer states; places agree in kind and input edge, and
/// the position along it within `1e-9`.
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
            return Some(format!("normal coordinate of edge {e}: {x} vs {y}"));
        }
    }
    if a.roundabouts() != b.roundabouts() {
        return Some("roundabouts differ".into());
    }
    for (v, (p, q)) in a
        .vertex_locations()
        .iter()
        .zip(b.vertex_locations())
        .enumerate()
    {
        let same = match (p, q) {
            (SurfacePoint::Edge { edge: x, t: s }, SurfacePoint::Edge { edge: y, t: u }) => {
                x == y && (s - u).abs() < 1e-9
            }
            _ => p == q,
        };
        if !same {
            return Some(format!("place of vertex {v}: {p:?} vs {q:?}"));
        }
    }
    None
}

#[derive(Default)]
struct Stats {
    splits: usize,
    snapped: usize,
    first: usize,
    last: usize,
    next_to_vertex: usize,
    boundary: usize,
    spokes_along: usize,
    chained: usize,
    refused: usize,
    disagree: usize,
    skipped: usize,
    failures: Vec<String>,
}

/// Splits edge `e` of both at parameter `t`; compares. `Ok(None)` if either refuses.
fn split_and_compare(
    s: &mut IntrinsicTriangulation,
    it: &mut IntegerTriangulation,
    e: usize,
    t: f64,
    st: &mut Stats,
) -> Result<bool, String> {
    let (sb, ib) = (s.clone(), it.clone());
    let n = it.normal_coordinate(e);
    let (rs, ri) = (s.split_edge(e, t), it.split_edge(e, t));
    let (xs, xi) = match (rs, ri) {
        (Ok(a), Ok(b)) => (a, b),
        (Err(_), Err(_)) => {
            st.refused += 1;
            *s = sb;
            *it = ib;
            return Ok(false);
        }
        (a, b) => {
            *s = sb;
            *it = ib;
            return Err(format!(
                "edge {e} at {t}: refusals differ (signposts {:?}, integers {:?})",
                a.map(|x| x.snapped),
                b.map(|x| x.snapped)
            ));
        }
    };
    st.splits += 1;
    if xs.snapped.is_some() != xi.snapped.is_some() {
        // Within rounding of the snapping tolerance the two traces decide differently.
        if xs.snapped.or(xi.snapped).is_some_and(|c| c < 1e-10) {
            st.disagree += 1;
            return Ok(false);
        }
        return Err(format!(
            "edge {e} at {t}: snapped {:?} vs {:?}",
            xs.snapped, xi.snapped
        ));
    }
    if xi.snapped.is_some() {
        st.snapped += 1;
        if let Err(err) = it.certify() {
            return Err(format!("edge {e} (n = {n}) at {t}: certify: {err}"));
        }
        match convert(s) {
            Ok((c, _)) => {
                if let Some(d) = differs(it, &c) {
                    return Err(format!("edge {e} (n = {n}) at {t}: {d}"));
                }
            }
            Err(err) => return Err(format!("edge {e}: conversion of the signposts: {err}")),
        }
        let spokes = &xi.new_edges[1..];
        if spokes.iter().any(|&x| it.normal_coordinate(x as usize) == -1) {
            st.spokes_along += 1;
        }
        Ok(true)
    } else {
        Ok(false)
    }
}

/// For every state, every edge across input edges, every crossing: split there (and a hair to
/// either side) in both representations.
#[test]
fn split_at_a_crossing_equals_the_signposts() {
    let shapes = all_shapes();
    let specs: Vec<(usize, Scenario)> = (0..shapes.len())
        .flat_map(|i| SCENARIOS.into_iter().map(move |sc| (i, sc)))
        .collect();
    let results = par_map(&specs, |&(i, scenario)| {
        let shape = &shapes[i];
        let tag = format!("{} / {scenario:?}", shape.name);
        let mut st = Stats::default();
        let s0 = state(shape, scenario, 31);
        let Ok((it0, _)) = convert(&s0) else {
            st.skipped += 1;
            return st;
        };
        let n_edges = it0.complex().n_edges();
        let stride = if n_edges > 3000 { 7 } else { 1 };
        for e in (0..n_edges).step_by(stride) {
            let n = it0.normal_coordinate(e);
            if n <= 0 {
                continue;
            }
            let Ok(trace) = s0.trace_edge(e, &EdgeTraceOptions::default()) else {
                continue;
            };
            let len = s0.edge_lengths()[e];
            let nx = trace.crossings.len();
            for (j, x) in trace.crossings.iter().enumerate() {
                for shift in [0.0, 1e-14, -3e-14] {
                    let t = (x.distance / len) + shift;
                    if !(t > 0.0 && t < 1.0) {
                        continue;
                    }
                    let (mut s, mut it) = (s0.clone(), it0.clone());
                    match split_and_compare(&mut s, &mut it, e, t, &mut st) {
                        Ok(true) => {
                            st.first += usize::from(j == 0);
                            st.last += usize::from(j + 1 == nx);
                            let dv = it0.complex().is_boundary_edge(e)
                                || [0, 1].iter().any(|&k| {
                                    let h = it0.complex().edge_halfedge(e);
                                    let h = if k == 0 {
                                        h
                                    } else {
                                        it0.complex().twin(h).unwrap_or(h)
                                    };
                                    let f = h / 3;
                                    (0..3).any(|c| {
                                        it0.complex().is_boundary_edge(it0.complex().edge(3 * f + c))
                                    })
                                });
                            st.boundary += usize::from(dv);
                            let near = x.distance < 1e-3 * len || len - x.distance < 1e-3 * len;
                            st.next_to_vertex += usize::from(near);
                            // A second split on the result, at a crossing of another edge.
                            if let Some(e2) = (0..it.complex().n_edges())
                                .skip(e % 7)
                                .step_by(11)
                                .find(|&e2| it.normal_coordinate(e2) > 0)
                                && let Ok(tr) =
                                    s.trace_edge(e2, &EdgeTraceOptions::default())
                                && let Some(x2) = tr.crossings.get(tr.crossings.len() / 2)
                            {
                                let l2 = s.edge_lengths()[e2];
                                if let Ok(true) = split_and_compare(
                                    &mut s,
                                    &mut it,
                                    e2,
                                    x2.distance / l2,
                                    &mut st,
                                ) {
                                    st.chained += 1;
                                }
                            }
                        }
                        Ok(false) => {}
                        Err(err) => {
                            st.failures.push(format!("{tag}: {err}"));
                        }
                    }
                }
            }
        }
        st
    });
    let mut t = Stats::default();
    for r in results {
        t.splits += r.splits;
        t.snapped += r.snapped;
        t.first += r.first;
        t.last += r.last;
        t.next_to_vertex += r.next_to_vertex;
        t.boundary += r.boundary;
        t.spokes_along += r.spokes_along;
        t.chained += r.chained;
        t.refused += r.refused;
        t.disagree += r.disagree;
        t.skipped += r.skipped;
        t.failures.extend(r.failures);
    }
    eprintln!(
        "splits at crossings: {} splits, {} snapped (first crossing {}, last {}, within 1e-3 of an \
         end {}, on a boundary face {}, a spoke along the input edge {}, chained second splits \
         {}); {} refused by both, {} within rounding of the tolerance decided differently, {} states \
         skipped; {} failures",
        t.splits,
        t.snapped,
        t.first,
        t.last,
        t.next_to_vertex,
        t.boundary,
        t.spokes_along,
        t.chained,
        t.refused,
        t.disagree,
        t.skipped,
        t.failures.len()
    );
    for f in t.failures.iter().take(20) {
        eprintln!("  {f}");
    }
    assert!(t.failures.is_empty());
    assert!(t.snapped > 100);
}
