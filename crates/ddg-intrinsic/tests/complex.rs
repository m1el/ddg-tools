//! The shared connectivity (`Complex`, SPEC N50) and the shared interface (`Triangulation`):
//! the index bookkeeping of every topological operation, removal of a degree-3 vertex,
//! constrained edges, and the guarantees of the trait.

mod common;

use common::{inner_bary, standard_shapes};
use ddg_core::generators::SplitMix64;
use ddg_core::topology::NONE;
use ddg_intrinsic::{
    Complex, DelaunayOptions, DelaunayStatus, EdgeTraceOptions, FlipRefusal, IntrinsicError,
    IntrinsicTriangulation, OptimalDelaunayOptions, RefineOptions, RemovalRefusal, Triangulation,
};

/// Per-halfedge and per-edge data a representation would carry: the tail of every halfedge
/// and the end points of every edge, which the operations' bookkeeping must keep right.
struct Tracked {
    tails: Vec<u32>,
    ends: Vec<[u32; 2]>,
    vertices: Vec<u32>,
}

impl Tracked {
    fn new(c: &Complex) -> Tracked {
        Tracked {
            tails: (0..c.n_halfedges())
                .map(|h| c.from_vertex(h) as u32)
                .collect(),
            ends: (0..c.n_edges())
                .map(|e| c.edge_vertices(e).map(|v| v as u32))
                .collect(),
            vertices: (0..c.n_vertices() as u32).collect(),
        }
    }

    /// Every halfedge leaves the vertex recorded for it, every edge has the recorded end
    /// points in the recorded order.
    fn assert_matches(&self, c: &Complex, what: &str) {
        c.check().unwrap_or_else(|e| panic!("{what}: {e}"));
        assert_eq!(self.tails.len(), c.n_halfedges(), "{what}");
        for h in 0..c.n_halfedges() {
            assert_eq!(
                self.tails[h],
                c.from_vertex(h) as u32,
                "{what}: halfedge {h}"
            );
        }
        for e in 0..c.n_edges() {
            let ends = c.edge_vertices(e).map(|v| v as u32);
            assert_eq!(self.ends[e], ends, "{what}: edge {e}");
        }
    }
}

#[test]
fn bookkeeping_of_flips_and_splits_keeps_tails_and_edge_directions() {
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut c = Complex::from_tri_mesh(&shape.mesh);
        let mut rng = SplitMix64::new(100 + i as u64);
        let chi = c.euler_characteristic();
        let mut tracked = Tracked::new(&c);
        for step in 0..200 {
            let op = rng.next_u64() % 3;
            apply_one(&mut c, &mut tracked, &mut rng, op);
            tracked.assert_matches(&c, &format!("{}: step {step}", shape.name));
            assert_eq!(c.euler_characteristic(), chi, "{}", shape.name);
        }
    }
}

fn apply_one(c: &mut Complex, tracked: &mut Tracked, rng: &mut SplitMix64, op: u64) {
    match op {
        0 => {
            let e = (rng.next_u64() % c.n_edges() as u64) as usize;
            let directions: Vec<[usize; 2]> =
                (0..c.n_edges()).map(|g| c.edge_vertices(g)).collect();
            if let Ok(flip) = c.flip_edge(e) {
                flip.move_halfedge_data(&mut tracked.tails);
                let [l, k] = flip.new_vertices;
                let [lk, kl] = flip.new_halfedges;
                tracked.tails[lk as usize] = l;
                tracked.tails[kl as usize] = k;
                tracked.ends[e] = [l, k];
                // A flip keeps the direction of every edge it does not flip.
                for (g, &d) in directions.iter().enumerate() {
                    if g != e {
                        assert_eq!(c.edge_vertices(g), d, "edge {g} after flipping {e}");
                    }
                }
            }
        }
        1 => {
            let f = (rng.next_u64() % c.n_faces() as u64) as usize;
            let corners = c.faces()[f];
            let split = c.split_face(f).unwrap();
            split.move_halfedge_data(&mut tracked.tails, NONE);
            let p = split.vertex;
            for (&corner, &out) in corners.iter().zip(&split.outgoing) {
                tracked.ends.push([corner, p]);
                let out = out as usize;
                assert_eq!(c.to_vertex(out) as u32, corner);
                tracked.tails[out] = p;
                tracked.tails[c.twin(out).unwrap()] = corner;
            }
            tracked.vertices.push(p);
        }
        _ => {
            let e = (rng.next_u64() % c.n_edges() as u64) as usize;
            let [a, b] = tracked.ends[e];
            let Ok(split) = c.split_edge(e) else {
                return;
            };
            split.move_halfedge_data(&mut tracked.tails, NONE);
            let p = split.vertex;
            tracked.ends.resize(c.n_edges(), [NONE; 2]);
            tracked.ends[e] = [a, p];
            tracked.ends[split.halves[1] as usize] = [p, b];
            for &g in &split.to_opposite {
                if g != NONE {
                    let h = c.edge_halfedge(g as usize);
                    assert_eq!(c.to_vertex(h) as u32, p, "edge {g} runs towards p");
                    tracked.ends[g as usize] = [c.from_vertex(h) as u32, p];
                }
            }
            let outgoing = split.outgoing();
            assert_eq!(outgoing.len(), if split.faces[1] == NONE { 2 } else { 4 });
            for h in outgoing {
                assert_eq!(c.from_vertex(h as usize) as u32, p);
                tracked.tails[h as usize] = p;
                if let Some(t) = c.twin(h as usize) {
                    tracked.tails[t] = c.from_vertex(t) as u32;
                }
            }
            // The boundary halfedge p → b comes first at a vertex on the boundary.
            assert_eq!(
                c.outgoing_halfedge(p as usize),
                3 * split.new_faces[0] as usize
            );
            tracked.vertices.push(p);
            assert_eq!(c.edge_vertices(e), [a as usize, p as usize]);
            assert_eq!(
                c.edge_vertices(split.halves[1] as usize),
                [p as usize, b as usize]
            );
        }
    }
}

#[test]
fn removal_renumbers_and_brings_the_data_along() {
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut c = Complex::from_tri_mesh(&shape.mesh);
        let mut rng = SplitMix64::new(200 + i as u64);
        let chi = c.euler_characteristic();
        let mut tracked = Tracked::new(&c);
        // Insert some vertices, mix with flips and edge splits, then remove every vertex of
        // degree three that was inserted, in random order.
        for _ in 0..40 {
            let op = rng.next_u64() % 3;
            apply_one(&mut c, &mut tracked, &mut rng, op);
        }
        let n_input = shape.mesh.n_vertices();
        let mut n_removed = 0;
        loop {
            let candidates: Vec<usize> = (n_input..c.n_vertices())
                .filter(|&v| !c.is_boundary_vertex(v) && c.vertex_degree(v) == 3)
                .filter(|&v| c.outgoing_halfedges(v).all(|h| c.to_vertex(h) != v))
                .collect();
            if candidates.is_empty() {
                break;
            }
            let v = candidates[(rng.next_u64() % candidates.len() as u64) as usize];
            let star: Vec<u32> = c
                .outgoing_halfedges(v)
                .map(|h| c.to_vertex(h) as u32)
                .collect();
            let removal = c.remove_vertex(v).unwrap();
            removal.apply_halfedge_data(&mut tracked.tails);
            removal.apply_edge_data(&mut tracked.ends);
            removal.apply_vertex_data(&mut tracked.vertices);
            // Vertex ids stored as data are renamed too.
            if let Some((old, new)) = removal.vertex_move {
                for t in tracked.tails.iter_mut() {
                    if *t == old {
                        *t = new;
                    }
                }
                for ends in tracked.ends.iter_mut() {
                    for x in ends.iter_mut() {
                        if *x == old {
                            *x = new;
                        }
                    }
                }
            }
            tracked.assert_matches(&c, &format!("{}: removal of {v}", shape.name));
            let rename = |x: u32| match removal.vertex_move {
                Some((old, new)) if x == old => new,
                _ => x,
            };
            let mut corners = removal.corners.to_vec();
            corners.sort_unstable();
            let mut expected: Vec<u32> = star.iter().map(|&x| rename(x)).collect();
            expected.sort_unstable();
            assert_eq!(corners, expected, "{}", shape.name);
            assert_eq!(c.faces()[removal.face as usize], removal.corners);
            assert_eq!(c.euler_characteristic(), chi);
            n_removed += 1;
        }
        // Every inserted vertex that is left has degree other than three or is on the
        // boundary.
        assert!(n_removed > 0, "{}: nothing removed", shape.name);
    }
}

#[test]
fn insertion_then_removal_restores_the_triangulation_exactly() {
    for shape in standard_shapes() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(5);
        for _ in 0..tri.n_edges() {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let _ = tri.flip_edge(e);
        }
        let before = tri.clone();
        let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
        let Ok(insertion) = tri.insert_vertex(f, inner_bary(&mut rng)) else {
            continue;
        };
        assert_eq!(insertion.halves, None);
        let removal = tri.remove_vertex(insertion.vertex as usize).unwrap();
        assert_eq!(removal.face as usize, f);
        // The two halfedges the insertion moved out of face f move back.
        assert_eq!(removal.halfedge_moves.len(), 2, "{}", shape.name);
        assert_eq!(tri, before, "{}", shape.name);
    }
}

#[test]
fn removal_of_an_earlier_vertex_keeps_the_triangulation_consistent() {
    for shape in standard_shapes() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(9);
        let mut inserted = Vec::new();
        for _ in 0..8 {
            let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
            if let Ok(insertion) = tri.insert_vertex(f, inner_bary(&mut rng)) {
                inserted.push(insertion);
            }
        }
        let Some(first) = inserted.first() else {
            continue;
        };
        let v = first.vertex as usize;
        let locations: Vec<_> = (0..tri.n_vertices())
            .filter(|&u| u != v)
            .map(|u| tri.vertex_location(u))
            .collect();
        match tri.remove_vertex(v) {
            Ok(removal) => {
                let mut after: Vec<_> = tri.vertex_locations().to_vec();
                let mut expected = locations;
                let key = |p: &ddg_intrinsic::SurfacePoint| format!("{p}");
                after.sort_by_key(key);
                expected.sort_by_key(key);
                assert_eq!(after, expected, "{}", shape.name);
                assert_eq!(removal.counts[3] as usize, tri.n_vertices());
                let report = tri.check().unwrap();
                assert!(report.max_signpost_error < 1e-9, "{}", shape.name);
                assert!(report.relative_area_error.abs() < 1e-12, "{}", shape.name);
                let traces = tri.trace_edges(&EdgeTraceOptions::default()).unwrap();
                tri.overlay(&traces).unwrap();
            }
            Err(IntrinsicError::RemovalRefused { reason, .. }) => {
                // The first inserted vertex has degree three unless a later insertion landed
                // in one of its faces.
                assert!(
                    matches!(reason, RemovalRefusal::Degree { .. }),
                    "{}: {reason}",
                    shape.name
                );
            }
            Err(e) => panic!("{}: {e}", shape.name),
        }
    }
}

#[test]
fn removal_is_refused_where_it_does_not_apply() {
    let shape = &standard_shapes()[6]; // grid 6x4, with boundary
    let mut tri = shape.triangulation();
    let before = tri.clone();
    let refused = |r: Result<_, IntrinsicError>| match r {
        Err(IntrinsicError::RemovalRefused { reason, .. }) => reason,
        other => panic!("{other:?}"),
    };
    assert_eq!(refused(tri.remove_vertex(0)), RemovalRefusal::InputVertex);
    assert!(matches!(
        tri.remove_vertex(tri.n_vertices()),
        Err(IntrinsicError::IndexOutOfRange { .. })
    ));
    // A vertex on a boundary edge.
    let boundary_edge = (0..tri.n_edges())
        .find(|&e| tri.is_boundary_edge(e))
        .unwrap();
    let on_boundary = tri.split_edge(boundary_edge, 0.5).unwrap().vertex as usize;
    assert_eq!(
        refused(tri.remove_vertex(on_boundary)),
        RemovalRefusal::BoundaryVertex
    );
    // A vertex of degree four (an interior edge split), and one with a constrained edge.
    let interior = (0..tri.n_edges())
        .find(|&e| !tri.is_boundary_edge(e))
        .unwrap();
    let split = tri.split_edge(interior, 0.5).unwrap();
    let halves = split.halves.unwrap();
    assert_eq!(
        refused(tri.remove_vertex(split.vertex as usize)),
        RemovalRefusal::Degree { degree: 4 }
    );
    let f = 0;
    let p = tri.insert_vertex(f, [0.3, 0.3, 0.4]).unwrap();
    tri.set_constrained(p.new_edges[0] as usize, true).unwrap();
    assert!(matches!(
        refused(tri.remove_vertex(p.vertex as usize)),
        RemovalRefusal::ConstrainedEdge { .. }
    ));
    let _ = (halves, before);
}

#[test]
fn constrained_edges_are_left_alone() {
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut rng = SplitMix64::new(300 + i as u64);
        let mut tri = shape.triangulation();
        // Constrain a third of the interior edges.
        let mut constrained = Vec::new();
        for e in 0..tri.n_edges() {
            if !tri.is_boundary_edge(e) && rng.next_u64().is_multiple_of(3) {
                tri.set_constrained(e, true).unwrap();
                constrained.push((e, tri.edge_vertices(e), tri.edge_length(e)));
            }
        }
        for &(e, _, _) in &constrained {
            assert_eq!(tri.can_flip(e).unwrap(), Err(FlipRefusal::Constrained));
            assert!(matches!(
                tri.flip_edge(e),
                Err(IntrinsicError::FlipRefused {
                    reason: FlipRefusal::Constrained,
                    ..
                })
            ));
        }
        let report = tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_eq!(report.n_constrained, constrained.len(), "{}", shape.name);
        assert_ne!(report.status, DelaunayStatus::FlipLimit);
        for &(e, ends, length) in &constrained {
            assert!(tri.is_constrained(e));
            assert_eq!(tri.edge_vertices(e), ends, "{}: edge {e}", shape.name);
            assert_eq!(tri.edge_length(e), length);
        }
        assert_eq!(
            tri.is_delaunay(&DelaunayOptions::default().criterion),
            report.status == DelaunayStatus::Delaunay
        );
        // A split of a constrained edge gives two constrained halves, directed like it.
        let splittable = constrained
            .iter()
            .find(|&&(e, _, _)| {
                let h = tri.complex().edge_halfedge(e);
                tri.complex().twin(h).is_some_and(|t| t / 3 != h / 3)
            })
            .copied();
        if let Some((e, [a, b], length)) = splittable {
            let split = Triangulation::split_edge(&mut tri, e, 0.25).unwrap();
            let [h0, h1] = split.halves.unwrap().map(|g| g as usize);
            assert_eq!(h0, e);
            assert!(tri.is_constrained(h0) && tri.is_constrained(h1));
            let p = split.vertex as usize;
            assert_eq!(tri.edge_vertices(h0), [a, p]);
            assert_eq!(tri.edge_vertices(h1), [p, b]);
            assert_eq!(tri.edge_length(h0), 0.25 * length);
            assert_eq!(tri.edge_length(h1), (1.0 - 0.25) * length);
            for &g in &split.new_edges[1..] {
                assert!(!tri.is_constrained(g as usize));
            }
        }
        // Refinement and the optimal Delaunay iteration refuse, and change nothing.
        let before = tri.clone();
        assert!(matches!(
            tri.refine(&RefineOptions::default()),
            Err(IntrinsicError::ConstrainedEdges { .. })
        ));
        assert!(matches!(
            tri.optimal_delaunay(&OptimalDelaunayOptions::default()),
            Err(IntrinsicError::ConstrainedEdges { .. })
        ));
        assert_eq!(tri, before);
        // Releasing the constraints lets Delaunay flipping finish the job, as from scratch.
        for e in 0..tri.n_edges() {
            tri.set_constrained(e, false).unwrap();
        }
        let report = tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_eq!(report.n_constrained, 0);
    }
}

/// Code written against the trait only.
fn delaunay_through_the_trait<T: Triangulation>(t: &mut T) -> (usize, Vec<f64>) {
    let report = t.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let c = t.complex();
    let angles = (0..c.n_halfedges()).map(|h| t.corner_angle(h)).collect();
    (report.n_flips, angles)
}

#[test]
fn the_trait_runs_the_same_code_as_the_signpost_representation() {
    for shape in standard_shapes() {
        let mut a = shape.triangulation();
        let mut b = a.clone();
        let report = a.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let (n_flips, angles) = delaunay_through_the_trait(&mut b);
        assert_eq!(a, b, "{}", shape.name);
        assert_eq!(report.n_flips, n_flips);
        for (h, angle) in angles.into_iter().enumerate() {
            assert_eq!(angle.to_bits(), a.corner_angle(h).to_bits());
        }
        for e in 0..a.n_edges().min(20) {
            let t =
                Triangulation::trace_edge_on_input(&a, e, &EdgeTraceOptions::default()).unwrap();
            assert_eq!(t, a.trace_edge(e, &EdgeTraceOptions::default()).unwrap());
        }
        for v in 0..a.n_vertices() {
            assert_eq!(Triangulation::vertex_location(&a, v), a.vertex_location(v));
        }
        let _: &dyn Fn(&IntrinsicTriangulation) -> &Complex = &|t| Triangulation::complex(t);
    }
}
