//! `Complex::merge_edge` and the merge of a straight star (SPEC N53; IC-SPEC §4.6, I27): the
//! inverse of an edge split, its renumbering, and removal of split vertices by the signposts.

mod common;

use common::{inner_bary, standard_shapes};
use ddg_core::generators::SplitMix64;
use ddg_core::topology::NONE;
use ddg_intrinsic::{
    Complex, DelaunayOptions, IntrinsicError, MERGE_STRAIGHT_TOLERANCE, RemovalRefusal,
    RemovalStep, plan_edge_merge, straight_pair,
};

/// The difference in ulps between two positive lengths.
fn ulps(a: f64, b: f64) -> u64 {
    a.to_bits().abs_diff(b.to_bits())
}

/// Splits followed at once by the merge give back the complex exactly, ids included,
/// whichever half is given and on interior and boundary edges, also after flips, face splits
/// and earlier splits (Δ-complex configurations).
#[test]
fn merge_right_after_a_split_restores_the_complex_exactly() {
    let mut n = [0usize; 2];
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut c = Complex::from_tri_mesh(&shape.mesh);
        let mut rng = SplitMix64::new(1300 + i as u64);
        for step in 0..300 {
            match rng.next_u64() % 4 {
                0 => {
                    let e = (rng.next_u64() % c.n_edges() as u64) as usize;
                    let _ = c.flip_edge(e);
                }
                1 => {
                    let f = (rng.next_u64() % c.n_faces() as u64) as usize;
                    c.split_face(f).unwrap();
                }
                _ => {
                    let e = (rng.next_u64() % c.n_edges() as u64) as usize;
                    let before = c.clone();
                    let Ok(split) = c.split_edge(e) else { continue };
                    let half = split.halves[(rng.next_u64() % 2) as usize] as usize;
                    let merge = c.merge_edge(split.vertex as usize, half).unwrap();
                    assert_eq!(c, before, "{}: step {step}, edge {e}", shape.name);
                    assert_eq!(merge.edge as usize, e);
                    assert!(merge.halfedge_moves.len() <= 3);
                    assert!(merge.edge_moves.is_empty() && merge.face_moves.is_empty());
                    assert_eq!(merge.vertex_move, None);
                    n[usize::from(split.faces[1] == NONE)] += 1;
                    // Keep the split for the next steps.
                    c.split_edge(e).unwrap();
                }
            }
        }
        c.check().unwrap();
    }
    eprintln!(
        "split then merge: {} interior, {} boundary edges, all exact",
        n[0], n[1]
    );
    assert!(n[0] > 1000 && n[1] > 50);
}

/// Many splits undone in the reverse order give back the complex exactly.
#[test]
fn many_splits_merged_in_reverse_order_restore_the_complex() {
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut c = Complex::from_tri_mesh(&shape.mesh);
        let mut rng = SplitMix64::new(1400 + i as u64);
        for _ in 0..c.n_edges() {
            let e = (rng.next_u64() % c.n_edges() as u64) as usize;
            let _ = c.flip_edge(e);
        }
        let original = c.clone();
        let mut splits = Vec::new();
        for _ in 0..2 * c.n_edges() {
            let e = (rng.next_u64() % c.n_edges() as u64) as usize;
            if let Ok(s) = c.split_edge(e) {
                splits.push(s);
            }
        }
        for s in splits.iter().rev() {
            c.merge_edge(s.vertex as usize, s.halves[1] as usize)
                .unwrap();
        }
        assert_eq!(c, original, "{}", shape.name);
    }
}

/// Per-halfedge and per-edge data carried through merges of arbitrary vertices of degree four
/// (any pair of opposite edges, whatever the geometry; the complex does not look at it) and of
/// boundary vertices of degree three, in random order: every halfedge keeps its tail, every
/// edge its end points, the merged edge the direction of the half whose id it keeps.
#[test]
fn merge_renumbers_and_brings_the_data_along() {
    let mut merged = [0usize; 2];
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut c = Complex::from_tri_mesh(&shape.mesh);
        let mut rng = SplitMix64::new(1500 + i as u64);
        let chi = c.euler_characteristic();
        for _ in 0..60 {
            match rng.next_u64() % 3 {
                0 => {
                    let e = (rng.next_u64() % c.n_edges() as u64) as usize;
                    let _ = c.flip_edge(e);
                }
                1 => {
                    let f = (rng.next_u64() % c.n_faces() as u64) as usize;
                    c.split_face(f).unwrap();
                }
                _ => {
                    let e = (rng.next_u64() % c.n_edges() as u64) as usize;
                    let _ = c.split_edge(e);
                }
            }
        }
        let n_input = shape.mesh.n_vertices();
        let mut tails: Vec<u32> = (0..c.n_halfedges())
            .map(|h| c.from_vertex(h) as u32)
            .collect();
        let mut ends: Vec<[u32; 2]> = (0..c.n_edges())
            .map(|e| c.edge_vertices(e).map(|v| v as u32))
            .collect();
        let mut ids: Vec<u32> = (0..c.n_vertices() as u32).collect();
        loop {
            let candidates: Vec<usize> = (n_input..c.n_vertices())
                .filter(|&v| {
                    let d = c.outgoing_halfedges(v).count();
                    let ok = if c.is_boundary_vertex(v) {
                        d == 2
                    } else {
                        d == 4
                    };
                    ok && c.outgoing_halfedges(v).all(|h| c.to_vertex(h) != v)
                })
                .collect();
            if candidates.is_empty() {
                break;
            }
            let v = candidates[(rng.next_u64() % candidates.len() as u64) as usize];
            let out: Vec<usize> = c.outgoing_halfedges(v).collect();
            let e = if c.is_boundary_vertex(v) {
                c.edge(out[0])
            } else {
                c.edge(out[(rng.next_u64() % 4) as usize])
            };
            let star = c.clone();
            let merge = match c.merge_edge(v, e) {
                Ok(m) => m,
                Err(IntrinsicError::RemovalRefused { reason, .. }) => {
                    panic!("{}: merge of {v} at {e}: {reason}", shape.name)
                }
                Err(err) => panic!("{err}"),
            };
            // The expected ends of the merged edge: the kept half with v replaced by the far
            // end of the other half.
            let [keep, dead] = merge.halves.map(|x| x as usize);
            let far = |g: usize| {
                let [a, b] = star.edge_vertices(g);
                (if a == v { b } else { a }) as u32
            };
            let mut kept = ends[keep];
            for x in kept.iter_mut() {
                if *x == v as u32 {
                    *x = far(dead);
                }
            }
            ends[keep] = kept;
            merge.apply_halfedge_data(&mut tails);
            merge.apply_edge_data(&mut ends);
            merge.apply_vertex_data(&mut ids);
            if let Some((old, new)) = merge.vertex_move {
                for t in tails.iter_mut().chain(ends.iter_mut().flatten()) {
                    if *t == old {
                        *t = new;
                    }
                }
            }
            c.check()
                .unwrap_or_else(|err| panic!("{}: merge of {v}: {err}", shape.name));
            assert_eq!(tails.len(), c.n_halfedges());
            for (h, &t) in tails.iter().enumerate() {
                assert_eq!(t, c.from_vertex(h) as u32, "{}: halfedge {h}", shape.name);
            }
            for (g, &x) in ends.iter().enumerate() {
                assert_eq!(
                    x,
                    c.edge_vertices(g).map(|y| y as u32),
                    "{}: edge {g}",
                    shape.name
                );
            }
            assert_eq!(ends[merge.edge as usize], merge.ends);
            assert_eq!(ids.len(), c.n_vertices());
            assert_eq!(c.euler_characteristic(), chi);
            assert_eq!(
                merge.counts,
                [
                    c.n_halfedges() as u32,
                    c.n_edges() as u32,
                    c.n_faces() as u32,
                    c.n_vertices() as u32
                ]
            );
            for f in merge.faces {
                if f != NONE {
                    let h = (0..3)
                        .map(|k| 3 * f as usize + k)
                        .find(|&h| c.edge(h) == merge.edge as usize);
                    assert!(
                        h.is_some(),
                        "{}: face {f} lacks the merged edge",
                        shape.name
                    );
                }
            }
            merged[usize::from(merge.faces[1] == NONE)] += 1;
        }
    }
    eprintln!(
        "merged {} interior and {} boundary vertices",
        merged[0], merged[1]
    );
    assert!(merged[0] > 100 && merged[1] > 10);
}

#[test]
fn merges_that_do_not_apply_are_refused_and_change_nothing() {
    let shape = &standard_shapes()[6]; // grid 6x4, with boundary
    let mut tri = shape.triangulation();
    let refused = |r: Result<_, IntrinsicError>| match r {
        Err(IntrinsicError::RemovalRefused { reason, .. }) => reason,
        other => panic!("{other:?}"),
    };
    let interior = (0..tri.n_edges())
        .find(|&e| !tri.is_boundary_edge(e))
        .unwrap();
    let split = tri.split_edge(interior, 0.3).unwrap();
    let p = split.vertex as usize;
    let [a, b] = split.halves.unwrap().map(|x| x as usize);
    let before = tri.clone();
    // An input vertex; an edge that is not at the vertex; the edges to the opposite corners,
    // which are not straight.
    assert_eq!(refused(tri.merge_edge(0, a)), RemovalRefusal::InputVertex);
    let far = (0..tri.n_edges())
        .find(|&e| !tri.edge_vertices(e).contains(&p))
        .unwrap();
    assert_eq!(
        refused(tri.merge_edge(p, far)),
        RemovalRefusal::NotAHalf { edge: far as u32 }
    );
    let [c, _] = split.new_edges[1..].try_into().unwrap();
    assert!(matches!(
        refused(tri.merge_edge(p, c as usize)),
        RemovalRefusal::NotStraight { .. }
    ));
    // A constrained edge to an opposite corner; one constrained half.
    tri.set_constrained(c as usize, true).unwrap();
    let constrained = tri.clone();
    assert_eq!(
        refused(tri.merge_edge(p, a)),
        RemovalRefusal::ConstrainedEdge { edge: c }
    );
    assert!(tri == constrained);
    tri.set_constrained(c as usize, false).unwrap();
    tri.set_constrained(b, true).unwrap();
    assert_eq!(
        refused(tri.merge_edge(p, a)),
        RemovalRefusal::ConstrainedEdge { edge: b as u32 }
    );
    tri.set_constrained(b, false).unwrap();
    assert!(tri == before);
    // Degree five.
    let f = tri.complex().face(tri.complex().outgoing_halfedge(p));
    let q = tri.insert_vertex(f, [0.2, 0.3, 0.5]).unwrap();
    let mut after = tri.clone();
    assert!(matches!(
        refused(after.merge_edge(p, a)),
        RemovalRefusal::MergeDegree {
            degree: 5,
            boundary: false
        }
    ));
    assert!(after == tri);
    assert!(matches!(
        tri.merge_edge(p, tri.n_edges()),
        Err(IntrinsicError::IndexOutOfRange { .. })
    ));
    let _ = q;
    // Both halves constrained: merged into one constrained edge.
    let mut tri = before.clone();
    tri.set_constrained(a, true).unwrap();
    tri.set_constrained(b, true).unwrap();
    let merge = tri.merge_edge(p, b).unwrap();
    assert!(tri.is_constrained(merge.edge as usize));
}

/// The signposts: a split undone by the merge gives back the triangulation bitwise, except
/// the merged length (`t ℓ + (1 − t) ℓ`), within an ulp or two; boundary edges included.
#[test]
fn signpost_merge_after_a_split_is_the_identity() {
    let (mut exact, mut total, mut worst, mut worst_excess) = (0usize, 0usize, 0u64, 0f64);
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(1600 + i as u64);
        for _ in 0..tri.n_edges() {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let _ = tri.flip_edge(e);
        }
        for _ in 0..200 {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let t = 0.05 + 0.9 * (rng.next_u64() % 10_000) as f64 / 10_000.0;
            let before = tri.clone();
            let Ok(split) = tri.split_edge(e, t) else {
                continue;
            };
            let p = split.vertex as usize;
            let plan = plan_edge_merge(tri.complex(), tri.edge_lengths(), p, e).unwrap();
            worst_excess = worst_excess
                .max(plan.excess[0].abs())
                .max(plan.excess[1].abs());
            assert!(straight_pair(tri.complex(), tri.edge_lengths(), p).is_some());
            tri.merge_edge(p, e).unwrap();
            assert_eq!(tri.complex(), before.complex(), "{}", shape.name);
            assert_eq!(tri.signposts(), before.signposts(), "{}", shape.name);
            assert_eq!(tri.vertex_locations(), before.vertex_locations());
            for g in 0..tri.n_edges() {
                assert_eq!(tri.edge_on_input(g), before.edge_on_input(g));
                if g != e {
                    assert_eq!(
                        tri.edge_length(g).to_bits(),
                        before.edge_length(g).to_bits()
                    );
                }
            }
            let d = ulps(tri.edge_length(e), before.edge_length(e));
            worst = worst.max(d);
            exact += usize::from(d == 0);
            total += 1;
            tri.check().unwrap();
            // Go on from the split state.
            tri = before;
            let _ = tri.split_edge(e, t);
        }
    }
    eprintln!(
        "signpost split then merge: {total} merges, {exact} with the length bitwise, the others \
         within {worst} ulp; largest |angle sum - pi| at a split vertex {worst_excess:e}"
    );
    assert!(total > 1000 && worst <= 2);
    assert!(worst_excess < 0.1 * MERGE_STRAIGHT_TOLERANCE);
}

/// Removal of edge-split vertices by the signposts, after Delaunay flipping: Alg. 8 where it
/// completes, the merge where it was stuck; every state consistent, nothing refused as stuck at
/// a straight star.
#[test]
fn signpost_removal_of_split_vertices_merges_straight_stars() {
    let mut counts = [0usize; 4]; // flips only, merge, stuck, boundary merged
    let mut refusals = std::collections::BTreeMap::new();
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(1700 + i as u64);
        let n_input = tri.n_input_vertices();
        for _ in 0..40 {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let t = 0.2 + 0.6 * (rng.next_u64() % 1000) as f64 / 1000.0;
            let _ = tri.split_edge(e, t);
            if rng.next_u64().is_multiple_of(4) {
                let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
                let _ = tri.insert_vertex(f, inner_bary(&mut rng));
            }
        }
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let mut v = tri.n_vertices();
        while v > n_input {
            v -= 1;
            let before = tri.clone();
            match tri.remove_inserted_vertex(v) {
                Ok(r) => {
                    match r.removal {
                        RemovalStep::Star(_) => counts[0] += 1,
                        RemovalStep::Merge(ref m) if m.faces[1] == NONE => counts[3] += 1,
                        RemovalStep::Merge(_) => counts[1] += 1,
                    }
                    let report = tri.check().unwrap();
                    assert!(report.max_signpost_error < 1e-9, "{}", shape.name);
                    assert!(report.relative_area_error.abs() < 1e-12, "{}", shape.name);
                }
                Err(IntrinsicError::RemovalRefused { reason, .. }) => {
                    assert!(tri == before);
                    let key = format!("{reason:?}");
                    *refusals
                        .entry(key.split([' ', '{']).next().unwrap().to_string())
                        .or_insert(0usize) += 1;
                    assert!(
                        !matches!(reason, RemovalRefusal::Stuck { .. })
                            || straight_pair(tri.complex(), tri.edge_lengths(), v).is_none(),
                        "{}: stuck at a straight star at {v}",
                        shape.name
                    );
                    counts[2] += 1;
                }
                Err(err) => panic!("{}: {err}", shape.name),
            }
        }
    }
    eprintln!(
        "signpost removal after Delaunay flipping: {} by flips (Alg. 8), {} by a merge, {} \
         boundary vertices merged, {} refused {refusals:?}",
        counts[0], counts[1], counts[3], counts[2]
    );
    assert!(counts[1] > 50 && counts[3] > 10);
}
