//! `Complex::remove_boundary_vertex` (SPEC N54): removal of a boundary vertex on two faces, its
//! renumbering, the round trip with a boundary edge split, Δ-complex cases and the refusals.

mod common;

use common::standard_shapes;
use ddg_core::TriMesh;
use ddg_core::generators::SplitMix64;
use ddg_core::topology::NONE;
use ddg_intrinsic::{Complex, IntrinsicError, RemovalRefusal};

fn refused<T: std::fmt::Debug>(r: Result<T, IntrinsicError>) -> RemovalRefusal {
    match r {
        Err(IntrinsicError::RemovalRefused { reason, .. }) => reason,
        other => panic!("expected a refusal, got {other:?}"),
    }
}

/// The boundary vertices on exactly two faces with no self-edge at them.
fn candidates(c: &Complex) -> Vec<usize> {
    (0..c.n_vertices())
        .filter(|&v| {
            c.is_boundary_vertex(v)
                && c.outgoing_halfedges(v).count() == 2
                && c.outgoing_halfedges(v).all(|h| c.to_vertex(h) != v)
        })
        .collect()
}

/// A boundary split followed at once by the removal of the new vertex gives back the complex
/// exactly, ids included, also after flips, face splits and earlier splits (Δ-complexes).
#[test]
fn removal_right_after_a_boundary_split_restores_the_complex_exactly() {
    let mut n = 0usize;
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut c = Complex::from_tri_mesh(&shape.mesh);
        let boundary: Vec<usize> = (0..c.n_edges())
            .filter(|&e| c.is_boundary_edge(e))
            .collect();
        if boundary.is_empty() {
            continue;
        }
        let mut rng = SplitMix64::new(5400 + i as u64);
        for step in 0..300 {
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
                    let boundary: Vec<usize> = (0..c.n_edges())
                        .filter(|&e| c.is_boundary_edge(e))
                        .collect();
                    let e = boundary[(rng.next_u64() % boundary.len() as u64) as usize];
                    let before = c.clone();
                    let Ok(split) = c.split_edge(e) else { continue };
                    let removal = c.remove_boundary_vertex(split.vertex as usize).unwrap();
                    assert_eq!(c, before, "{}: step {step}, edge {e}", shape.name);
                    assert_eq!(removal.edge as usize, e);
                    assert_eq!(removal.faces[1], NONE);
                    assert!(removal.edge_moves.is_empty() && removal.face_moves.is_empty());
                    assert_eq!(removal.vertex_move, None);
                    n += 1;
                    // Keep a split for the next steps.
                    c.split_edge(e).unwrap();
                }
            }
        }
        c.check().unwrap();
    }
    eprintln!("boundary split then removal: {n}, all exact");
    assert!(n > 300);
}

/// Per-halfedge, per-edge and per-vertex data carried through removals of arbitrary boundary
/// vertices on two faces (the complex does not look at the geometry), in random order, with
/// flips in between: every halfedge keeps its tail, every edge its end points (the merged one
/// runs from `b` to `c`), the merged face is `abc` with its outer halfedges, χ is unchanged.
#[test]
fn removal_renumbers_and_brings_the_data_along() {
    let mut removed = 0usize;
    let mut loops = 0usize;
    for (i, shape) in standard_shapes().iter().enumerate() {
        let mut c = Complex::from_tri_mesh(&shape.mesh);
        if !(0..c.n_vertices()).any(|v| c.is_boundary_vertex(v)) {
            continue;
        }
        let mut rng = SplitMix64::new(5500 + i as u64);
        let chi = c.euler_characteristic();
        let mut tails: Vec<u32> = (0..c.n_halfedges())
            .map(|h| c.from_vertex(h) as u32)
            .collect();
        let mut ends: Vec<[u32; 2]> = (0..c.n_edges())
            .map(|e| c.edge_vertices(e).map(|v| v as u32))
            .collect();
        let mut ids: Vec<u32> = (0..c.n_vertices() as u32).collect();
        for _ in 0..4 * shape.mesh.n_vertices() {
            // Flips bring more boundary vertices down to two faces.
            for _ in 0..3 {
                let e = (rng.next_u64() % c.n_edges() as u64) as usize;
                if let Ok(flip) = c.flip_edge(e) {
                    flip.move_halfedge_data(&mut tails);
                    for h in flip.new_halfedges {
                        tails[h as usize] = c.from_vertex(h as usize) as u32;
                    }
                    ends[e] = c.edge_vertices(e).map(|v| v as u32);
                }
            }
            let cand = candidates(&c);
            if cand.is_empty() {
                continue;
            }
            let v = cand[(rng.next_u64() % cand.len() as u64) as usize];
            let out: Vec<usize> = c.outgoing_halfedges(v).collect();
            // b → v arrives, v → c leaves, v → a is the spoke.
            let (vc, va) = (out[0], out[1]);
            let (cc, a) = (c.to_vertex(vc) as u32, c.to_vertex(va) as u32);
            let b = c.from_vertex(c.prev(va)) as u32;
            let star = c.clone();
            let removal = c
                .remove_boundary_vertex(v)
                .unwrap_or_else(|e| panic!("{}: removal of {v}: {e}", shape.name));
            let [keep, _] = removal.halves.map(|x| x as usize);
            ends[keep] = [b, cc];
            removal.apply_halfedge_data(&mut tails);
            removal.apply_edge_data(&mut ends);
            removal.apply_vertex_data(&mut ids);
            let rename = |x: u32| match removal.vertex_move {
                Some((old, new)) if x == old => new,
                _ => x,
            };
            if removal.vertex_move.is_some() {
                for t in tails.iter_mut().chain(ends.iter_mut().flatten()) {
                    *t = rename(*t);
                }
            }
            c.check()
                .unwrap_or_else(|err| panic!("{}: removal of {v}: {err}", shape.name));
            assert_eq!(c.euler_characteristic(), chi, "{}", shape.name);
            assert_eq!(
                (c.n_vertices(), c.n_edges(), c.n_faces()),
                (
                    star.n_vertices() - 1,
                    star.n_edges() - 2,
                    star.n_faces() - 1
                )
            );
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
            assert_eq!(removal.ends, [rename(b), rename(cc)]);
            assert!(c.is_boundary_edge(removal.edge as usize));
            let f = removal.faces[0] as usize;
            let mut corners = c.faces()[f];
            let want = [rename(a), rename(b), rename(cc)];
            // The face is abc up to rotation.
            let r = (0..3).find(|&k| corners[k] == want[0] && corners[(k + 1) % 3] == want[1]);
            let r = r.unwrap_or_else(|| panic!("{}: face {corners:?}, want {want:?}", shape.name));
            corners.rotate_left(r);
            assert_eq!(corners, want);
            assert_eq!(ids.len(), c.n_vertices());
            if b == cc {
                loops += 1;
            }
            removed += 1;
        }
    }
    eprintln!("removed {removed} boundary vertices ({loops} leaving a boundary self-edge)");
    assert!(removed > 200);
    assert!(loops > 0);
}

/// Two triangles `(i, c, a)`, `(i, a, b)`: removing `i` leaves the triangle `abc`; the ids are
/// filled by the last ones.
#[test]
fn removal_of_a_quad_corner() {
    // Vertices: i = 0, a = 1, b = 2, c = 3; faces (i, c, a), (i, a, b).
    let mesh = TriMesh::new(4, vec![[0, 3, 1], [0, 1, 2]]).unwrap();
    let mut c = Complex::from_tri_mesh(&mesh);
    let removal = c.remove_boundary_vertex(0).unwrap();
    c.check().unwrap();
    assert_eq!(removal.vertex, 0);
    assert_eq!(removal.vertex_move, Some((3, 0)));
    assert_eq!(c.n_faces(), 1);
    assert_eq!(c.n_edges(), 3);
    assert_eq!(c.n_vertices(), 3);
    // c (= 3) was renamed to 0.
    let mut f = c.faces()[0];
    let r = (0..3).find(|&k| f[k] == 1).unwrap();
    f.rotate_left(r);
    assert_eq!(f, [1, 2, 0]);
    assert_eq!(removal.ends, [2, 0]);
}

#[test]
fn removals_that_do_not_apply_are_refused_and_change_nothing() {
    // A fan of three triangles around the boundary vertex 0, and an interior vertex.
    let mesh = TriMesh::new(5, vec![[0, 1, 2], [0, 2, 3], [0, 3, 4]]).unwrap();
    let mut c = Complex::from_tri_mesh(&mesh);
    let before = c.clone();
    assert_eq!(
        refused(c.remove_boundary_vertex(0)),
        RemovalRefusal::BoundaryFaces { faces: 3 }
    );
    assert_eq!(
        refused(c.remove_boundary_vertex(1)),
        RemovalRefusal::BoundaryFaces { faces: 1 }
    );
    assert!(matches!(
        c.remove_boundary_vertex(9),
        Err(IntrinsicError::IndexOutOfRange { .. })
    ));
    assert_eq!(c, before);
    // Vertex 2 is on two faces; its spoke is the edge 0-2.
    let spoke = (0..c.n_edges())
        .find(|&e| {
            let mut v = c.edge_vertices(e);
            v.sort_unstable();
            v == [0, 2]
        })
        .unwrap();
    c.set_constrained(spoke, true).unwrap();
    let before = c.clone();
    assert_eq!(
        refused(c.remove_boundary_vertex(2)),
        RemovalRefusal::ConstrainedEdge { edge: spoke as u32 }
    );
    assert_eq!(c, before);
    c.set_constrained(spoke, false).unwrap();
    assert_eq!(c.boundary_removal_refusal(2), None);

    // Interior vertex: the centre of an icosahedron-like closed shape.
    let shape = &standard_shapes()[0];
    let mut closed = Complex::from_tri_mesh(&shape.mesh);
    let before = closed.clone();
    assert_eq!(
        refused(closed.remove_boundary_vertex(0)),
        RemovalRefusal::InteriorVertex
    );
    assert_eq!(closed, before);
}

/// A boundary self-edge: the only vertex of a boundary loop is refused (removing it would
/// leave a boundary loop without vertices).
#[test]
fn the_vertex_of_a_boundary_loop_is_refused() {
    // Random flips and removals on an annulus shrink its boundary loops until some loop has
    // one vertex; that vertex is refused.
    let shape = standard_shapes()
        .into_iter()
        .find(|s| s.name.starts_with("annulus"))
        .unwrap();
    let mut c = Complex::from_tri_mesh(&shape.mesh);
    let mut rng = SplitMix64::new(5600);
    let mut seen = 0usize;
    for _ in 0..20_000 {
        let e = (rng.next_u64() % c.n_edges() as u64) as usize;
        let _ = c.flip_edge(e);
        let cand = candidates(&c);
        if !cand.is_empty() {
            let v = cand[(rng.next_u64() % cand.len() as u64) as usize];
            c.remove_boundary_vertex(v).unwrap();
        }
        for v in 0..c.n_vertices() {
            let h = c.outgoing_halfedge(v);
            if c.is_boundary_vertex(v) && c.to_vertex(h) == v {
                let before = c.clone();
                let reason = refused(c.remove_boundary_vertex(v));
                assert!(
                    matches!(
                        reason,
                        RemovalRefusal::SelfEdge { .. } | RemovalRefusal::BoundaryFaces { .. }
                    ),
                    "{reason:?}"
                );
                assert_eq!(c, before);
                seen += 1;
            }
        }
        c.check().unwrap();
    }
    eprintln!("boundary self-edges refused: {seen} times");
    assert!(seen > 0);
}
