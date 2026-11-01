//! `Surface` validation, orientation classes, components, fans, boundary loops, geometry and
//! OBJ I/O (`SPEC.md` §2.1, §2.3, §9.1).

#![allow(clippy::needless_range_loop)]

use ddg_core::generators::{annulus, bumpy_grid_disk, icosphere, perturb, torus};
use ddg_core::operators::{cotan_laplacian, lumped_mass as core_mass};
use ddg_core::{Geometry, TriMesh};
use ddg_repulsive_surfaces::geometry::{
    FaceGeometry, cotan_stiffness, face_area_gradient, face_normal_variation, face_quantities,
    lumped_mass, vertex_areas,
};
use ddg_repulsive_surfaces::{
    NONE, Orientation, RsError, Surface, obj_to_string, parse_obj, read_obj, write_obj,
};

const KLEIN: &str = "/var/data/discrete-geometry/RepulsiveSurfaces/repulsive-surfaces/scenes/KleinBottle/klein_bottle_blocky.obj";

/// A Möbius band with `n` quads (2 triangles each) around the core circle of radius 2,
/// half-width 0.5, one half twist.
fn moebius(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut pos = Vec::new();
    for i in 0..n {
        let th = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let (c, s) = (th.cos(), th.sin());
        let (ch, sh) = ((th / 2.0).cos(), (th / 2.0).sin());
        for w in [0.5, -0.5] {
            pos.push([(2.0 + w * ch) * c, (2.0 + w * ch) * s, w * sh]);
        }
    }
    let t = |i: usize| (2 * (i % n)) as u32;
    let b = |i: usize| (2 * (i % n) + 1) as u32;
    let mut faces = Vec::new();
    for i in 0..n {
        let (ti, bi) = (t(i), b(i));
        // Across the seam the top and bottom rows swap.
        let (tj, bj) = if i + 1 == n {
            (b(0), t(0))
        } else {
            (t(i + 1), b(i + 1))
        };
        faces.push([ti, bi, bj]);
        faces.push([ti, bj, tj]);
    }
    (pos, faces)
}

fn err(n: usize, faces: Vec<[u32; 3]>) -> RsError {
    Surface::new(n, faces).expect_err("must be rejected")
}

#[test]
fn bad_topology_is_rejected_with_located_errors() {
    assert_eq!(err(3, vec![]), RsError::NoFaces);
    assert!(matches!(
        err(u32::MAX as usize, vec![[0, 1, 2]]),
        RsError::TooLarge { .. }
    ));
    let e = err(3, vec![[0, 1, 5]]);
    assert_eq!(
        e,
        RsError::VertexOutOfRange {
            face: 0,
            corner: 2,
            vertex: 5,
            n_vertices: 3
        }
    );
    assert!(
        e.to_string()
            .contains("face 0 corner 2 references vertex 5")
    );
    let e = err(4, vec![[0, 1, 2], [1, 3, 1]]);
    assert_eq!(
        e,
        RsError::RepeatedVertex {
            face: 1,
            vertices: [1, 3, 1]
        }
    );
    // Duplicates in either orientation.
    let e = err(4, vec![[0, 1, 2], [1, 2, 3], [2, 1, 0]]);
    assert_eq!(
        e,
        RsError::DuplicateFace {
            first: 0,
            second: 2,
            vertices: [0, 1, 2]
        }
    );
    assert!(e.to_string().contains("faces 0 and 2 are duplicates"));
    let e = err(5, vec![[0, 1, 2], [1, 0, 3], [0, 1, 4]]);
    assert_eq!(
        e,
        RsError::NonManifoldEdge {
            edge: [0, 1],
            faces: vec![0, 1, 2]
        }
    );
    assert_eq!(e.to_string(), "edge (0, 1) is shared by 3 faces: 0, 1, 2");
    // Bow tie: two triangles touching at vertex 2.
    let e = err(5, vec![[0, 1, 2], [2, 3, 4]]);
    assert_eq!(e, RsError::NonManifoldVertex { vertex: 2, fans: 2 });
    assert!(e.to_string().contains("vertex 2 is non-manifold"));
    // Two closed fans sharing a vertex (two tetrahedra glued at a vertex).
    let mut faces = vec![[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]];
    faces.extend([[0, 4, 5], [0, 5, 6], [0, 6, 4], [4, 6, 5]]);
    assert_eq!(
        err(7, faces),
        RsError::NonManifoldVertex { vertex: 0, fans: 2 }
    );
    let e = err(6, vec![[0, 1, 2], [0, 2, 4]]);
    assert_eq!(
        e,
        RsError::UnreferencedVertex {
            vertex: 3,
            count: 2
        }
    );
}

#[test]
fn bad_positions_are_rejected() {
    let s = Surface::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
    let mut p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    assert!(s.validate_positions(&p).is_ok());
    assert_eq!(
        s.validate_positions(&p[..3]),
        Err(RsError::PositionCountMismatch {
            expected: 4,
            found: 3
        })
    );
    p[2][1] = f64::NAN;
    assert_eq!(
        s.validate_positions(&p),
        Err(RsError::NonFinitePosition { vertex: 2 })
    );
    // Vertex 3 on the line through 0 and 2: face 1 collapses.
    p[2] = [1.0, 1.0, 0.0];
    p[3] = [2.0, 2.0, 0.0];
    let e = s.validate_positions(&p).unwrap_err();
    assert!(
        matches!(e, RsError::DegenerateFace { face: 1, .. }),
        "{e:?}"
    );
    assert!(e.to_string().contains("face 1 is degenerate"));
    // Coincident vertices.
    p[3] = p[0];
    assert!(matches!(
        s.validate_positions(&p),
        Err(RsError::DegenerateFace { face: 1, .. })
    ));
}

#[test]
fn icosphere_is_consistent_and_matches_trimesh() {
    let m = icosphere(2).unwrap();
    let s = Surface::new(m.positions.len(), m.faces.clone()).unwrap();
    let t = TriMesh::new(m.positions.len(), m.faces.clone()).unwrap();
    assert_eq!(s.n_components(), 1);
    assert_eq!(s.orientation(0), &Orientation::Consistent);
    assert_eq!(s.components()[0].inconsistent_edge(), None);
    assert_eq!(s.euler_characteristic(), 2);
    assert_eq!(s.components()[0].euler_characteristic(), 2);
    assert!(s.boundary_loops().is_empty());
    assert_eq!(s.edges(), t.edges());
    for v in 0..s.n_vertices() as u32 {
        assert!(!s.is_boundary_vertex(v));
        assert_eq!(s.fan(v).len(), t.vertex_face_count(v as usize));
        check_fan(&s, v);
    }
    for (e, &[a, b]) in s.edges().iter().enumerate() {
        assert_eq!(s.find_edge(b, a), Some(e as u32));
        assert!(!s.is_boundary_edge(e as u32));
        for f in s.edge_faces()[e] {
            assert!(s.faces()[f as usize].contains(&a) && s.faces()[f as usize].contains(&b));
        }
    }
    for (f, fe) in s.face_edges().iter().enumerate() {
        let t = s.faces()[f];
        for k in 0..3 {
            let [a, b] = s.edges()[fe[k] as usize];
            let (x, y) = (t[k], t[(k + 1) % 3]);
            assert_eq!([a, b], [x.min(y), x.max(y)]);
        }
    }
}

/// The fan of `v` is a chain of faces containing `v` in which consecutive faces share an edge
/// at `v`; a boundary fan starts and ends with boundary edges, an interior fan closes.
fn check_fan(s: &Surface, v: u32) {
    let fan = s.fan(v);
    let edges_at = |f: u32| -> Vec<u32> {
        s.face_edges()[f as usize]
            .iter()
            .copied()
            .filter(|&e| s.edges()[e as usize].contains(&v))
            .collect()
    };
    let shared = |f: u32, g: u32| edges_at(f).iter().any(|e| edges_at(g).contains(e));
    for w in fan.windows(2) {
        assert!(shared(w[0], w[1]), "fan of {v}: {fan:?}");
    }
    let has_boundary = |f: u32| edges_at(f).iter().any(|&e| s.is_boundary_edge(e));
    if s.is_boundary_vertex(v) {
        assert!(has_boundary(fan[0]) && has_boundary(fan[fan.len() - 1]));
    } else {
        assert!(shared(fan[0], fan[fan.len() - 1]));
    }
}

#[test]
fn flipped_faces_give_orientable_with_flip_set() {
    let m = icosphere(1).unwrap();
    let mut faces = m.faces.clone();
    for f in [3usize, 7, 40] {
        faces[f].swap(1, 2);
    }
    let s = Surface::new(m.positions.len(), faces).unwrap();
    assert_eq!(
        s.orientation(0),
        &Orientation::Orientable {
            flips: vec![3, 7, 40]
        }
    );
    let [a, b] = s.components()[0].inconsistent_edge().unwrap();
    // The witness edge belongs to a flipped face.
    assert!([3usize, 7, 40].iter().any(|&f| {
        let t = s.faces()[f];
        t.contains(&a) && t.contains(&b)
    }));
    // Lowest face flipped: the lowest face keeps its orientation, all others flip.
    let mut faces = m.faces.clone();
    faces[0].swap(0, 1);
    let s = Surface::new(m.positions.len(), faces).unwrap();
    let Orientation::Orientable { flips } = s.orientation(0) else {
        panic!("expected Orientable");
    };
    assert_eq!(flips.len(), m.faces.len() - 1);
    assert!(!flips.contains(&0));
}

#[test]
fn moebius_band_is_non_orientable_with_one_boundary_loop() {
    let (pos, faces) = moebius(12);
    assert!(TriMesh::new(pos.len(), faces.clone()).is_err());
    let s = Surface::new(pos.len(), faces).unwrap();
    s.validate_positions(&pos).unwrap();
    assert_eq!(s.n_components(), 1);
    assert_eq!(s.orientation(0), &Orientation::NonOrientable);
    assert!(s.components()[0].inconsistent_edge().is_some());
    assert_eq!(s.euler_characteristic(), 0);
    assert_eq!(s.boundary_loops().len(), 1);
    assert_eq!(s.boundary_loops()[0].len(), 24);
    assert_eq!(s.components()[0].boundary_loops(), &[0]);
    assert_eq!(s.boundary_loops()[0][0], 0);
    for v in 0..s.n_vertices() as u32 {
        assert!(s.is_boundary_vertex(v));
        check_fan(&s, v);
    }
    check_loops(&s);
}

/// Every boundary loop is a cycle of boundary edges, each boundary edge on exactly one loop.
fn check_loops(s: &Surface) {
    let mut seen = vec![false; s.n_edges()];
    for lp in s.boundary_loops() {
        assert_eq!(lp[0], *lp.iter().min().unwrap());
        for i in 0..lp.len() {
            let e = s
                .find_edge(lp[i], lp[(i + 1) % lp.len()])
                .expect("loop edge");
            assert!(s.is_boundary_edge(e));
            assert!(!seen[e as usize]);
            seen[e as usize] = true;
        }
    }
    let n_boundary = (0..s.n_edges() as u32)
        .filter(|&e| s.is_boundary_edge(e))
        .count();
    assert_eq!(seen.iter().filter(|&&x| x).count(), n_boundary);
}

#[test]
fn klein_bottle_from_the_reference_is_accepted() {
    let (s, pos) = read_obj(KLEIN).unwrap();
    assert_eq!((s.n_vertices(), s.n_faces()), (900, 1800));
    s.validate_positions(&pos).unwrap();
    assert_eq!(s.n_components(), 1);
    assert_eq!(s.orientation(0), &Orientation::NonOrientable);
    assert_eq!(s.euler_characteristic(), 0);
    assert!(s.boundary_loops().is_empty());
    for v in 0..s.n_vertices() as u32 {
        assert!(!s.is_boundary_vertex(v));
        check_fan(&s, v);
    }
    // ddg_core refuses it (it needs an orientation).
    assert!(TriMesh::new(s.n_vertices(), s.faces().to_vec()).is_err());
}

#[test]
fn components_and_boundaries() {
    // Annulus (two loops) + a translated icosphere + a bumpy disk (one loop).
    let a = annulus(12, 3, 0.5, 1.0).unwrap();
    let b = icosphere(1).unwrap();
    let c = bumpy_grid_disk(4, 3, 0.2).unwrap();
    let mut pos = a.positions.clone();
    let mut faces = a.faces.clone();
    for m in [&b, &c] {
        let off = pos.len() as u32;
        let shift = 3.0 * (off as f64).cbrt();
        pos.extend(m.positions.iter().map(|p| [p[0] + shift, p[1], p[2]]));
        faces.extend(m.faces.iter().map(|t| [t[0] + off, t[1] + off, t[2] + off]));
    }
    let s = Surface::new(pos.len(), faces.clone()).unwrap();
    let t = TriMesh::new(pos.len(), faces.clone()).unwrap();
    assert_eq!(s.n_components(), 3);
    assert_eq!(s.face_components(), t.face_components());
    let fa = a.faces.len();
    let fb = b.faces.len();
    assert_eq!(s.components()[0].faces()[0], 0);
    assert_eq!(s.components()[1].faces()[0] as usize, fa);
    assert_eq!(s.components()[2].faces()[0] as usize, fa + fb);
    let chis: Vec<i64> = s
        .components()
        .iter()
        .map(|c| c.euler_characteristic())
        .collect();
    assert_eq!(chis, vec![0, 2, 1]);
    assert_eq!(s.components()[0].boundary_loops().len(), 2);
    assert_eq!(s.components()[1].boundary_loops().len(), 0);
    assert_eq!(s.components()[2].boundary_loops().len(), 1);
    for c in s.components() {
        assert_eq!(c.orientation(), &Orientation::Consistent);
        for &v in c.vertices() {
            assert_eq!(
                s.vertex_components()[v as usize] as usize,
                s.face_components()[c.faces()[0] as usize] as usize
            );
        }
    }
    for v in 0..s.n_vertices() as u32 {
        assert_eq!(s.is_boundary_vertex(v), t.is_boundary_vertex(v as usize));
        check_fan(&s, v);
    }
    check_loops(&s);
    // On consistently oriented components the loops run as ddg_core's (surface on the left).
    assert_eq!(s.boundary_loops().len(), t.boundary_loops().len());
    for lp in s.boundary_loops() {
        let other = t
            .boundary_loops()
            .iter()
            .find(|l| l.contains(&lp[0]))
            .expect("same loop");
        assert_eq!(other.len(), lp.len());
        let k = other.iter().position(|&x| x == lp[0]).unwrap();
        for i in 0..lp.len() {
            assert_eq!(lp[i], other[(k + i) % other.len()]);
        }
    }
    // Boundary edges have NONE as their second face.
    let nb = s.edge_faces().iter().filter(|f| f[1] == NONE).count();
    assert_eq!(nb, t.n_boundary_edges());
}

#[test]
fn geometry_matches_ddg_core_and_is_orientation_free() {
    for (name, m) in [
        ("icosphere", icosphere(2).unwrap()),
        ("torus", torus(10, 7, 2.0, 0.7).unwrap()),
        ("disk", bumpy_grid_disk(5, 4, 0.3).unwrap()),
    ] {
        let mut pos = m.positions.clone();
        perturb(&mut pos, [0.02; 3], 7);
        let s = Surface::new(pos.len(), m.faces.clone()).unwrap();
        let t = TriMesh::new(pos.len(), m.faces.clone()).unwrap();
        let g = Geometry::from_positions(&t, &pos).unwrap();
        let w_core = cotan_laplacian(&t, &g).unwrap();
        let w = cotan_stiffness(&s, &pos).unwrap();
        let mass = lumped_mass(&s, &pos).unwrap();
        let m_core = core_mass(&t, &g).unwrap();
        let n = pos.len();
        let scale = (0..n).map(|i| w_core.get(i, i).abs()).fold(0.0, f64::max);
        for i in 0..n {
            assert!(
                (mass[i] - m_core[i]).abs() <= 1e-14 * m_core[i],
                "{name}: mass {i}"
            );
            let mut row = 0.0;
            for j in 0..n {
                let d = (w.get(i, j) - w_core.get(i, j)).abs();
                assert!(d <= 1e-12 * scale, "{name}: W[{i},{j}]");
                assert_eq!(w.get(i, j), w.get(j, i));
                row += w.get(i, j);
            }
            assert!(row.abs() <= 1e-12 * scale);
        }
        // PSD (diagonal dominance is not guaranteed; check x^T W x >= 0 on a few vectors).
        let mut rng = ddg_core::generators::SplitMix64::new(3);
        for _ in 0..5 {
            let x: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let mut q = 0.0;
            for (r, c, val) in w.triplets() {
                q += x[r as usize] * val * x[c as usize];
            }
            assert!(q >= -1e-12 * scale, "{name}: x^T W x = {q}");
        }
        // Flip every other face: W, M, areas and barycenters unchanged, normals negated.
        let flipped: Vec<[u32; 3]> = m
            .faces
            .iter()
            .enumerate()
            .map(|(f, t)| if f % 2 == 0 { [t[0], t[2], t[1]] } else { *t })
            .collect();
        let sf = Surface::new(pos.len(), flipped).unwrap();
        assert!(!matches!(sf.orientation(0), Orientation::Consistent));
        let wf = cotan_stiffness(&sf, &pos).unwrap();
        for i in 0..n {
            for j in 0..n {
                assert!((wf.get(i, j) - w.get(i, j)).abs() <= 1e-13 * scale);
            }
        }
        let g0 = FaceGeometry::new(&s, &pos).unwrap();
        let g1 = FaceGeometry::new(&sf, &pos).unwrap();
        for f in 0..s.n_faces() {
            assert!((g0.areas[f] - g1.areas[f]).abs() <= 1e-15 * g0.areas[f]);
            let sign = if f % 2 == 0 { -1.0 } else { 1.0 };
            for k in 0..3 {
                assert!((g0.normals[f][k] - sign * g1.normals[f][k]).abs() <= 1e-15);
                assert!((g0.barycenters[f][k] - g1.barycenters[f][k]).abs() <= 1e-15);
            }
        }
        let va = vertex_areas(&s, &g0.areas);
        let total: f64 = va.iter().sum();
        assert!((total - g0.total_area()).abs() <= 1e-13 * total);
    }
}

#[test]
fn face_derivatives_match_finite_differences() {
    let p0 = [[0.1, -0.3, 0.2], [1.2, 0.1, -0.4], [0.3, 0.9, 0.5]];
    let ga = face_area_gradient(p0);
    let jn = face_normal_variation(p0);
    let h = 1e-6;
    for c in 0..3 {
        for s in 0..3 {
            let mut pp = p0;
            let mut pm = p0;
            pp[c][s] += h;
            pm[c][s] -= h;
            let (ap, np, _) = face_quantities(pp);
            let (am, nm, _) = face_quantities(pm);
            let da = (ap - am) / (2.0 * h);
            assert!((da - ga[c][s]).abs() < 1e-9, "area c{c} s{s}");
            for r in 0..3 {
                let dn = (np[r] - nm[r]) / (2.0 * h);
                assert!((dn - jn[c][r][s]).abs() < 1e-9, "normal c{c} r{r} s{s}");
            }
        }
    }
    // Barycenter and area of a known triangle.
    let (a, n, x) = face_quantities([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]]);
    assert_eq!(a, 3.0);
    assert_eq!(n, [0.0, 0.0, 1.0]);
    assert_eq!(x, [2.0 / 3.0, 1.0, 0.0]);
}

#[test]
fn obj_round_trip_and_polygons() {
    let (pos, faces) = moebius(8);
    let s = Surface::new(pos.len(), faces).unwrap();
    let text = obj_to_string(&s, &pos).unwrap();
    let (s2, pos2) = parse_obj(&text).unwrap();
    assert_eq!(s2, s);
    assert_eq!(pos2, pos);
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("rs_moebius.obj");
    write_obj(&path, &s, &pos).unwrap();
    let (s3, pos3) = read_obj(&path).unwrap();
    assert_eq!((s3, pos3), (s, pos.clone()));
    let quad = "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 0 0 1\nf 1 2 5\nf 1 2 3 4\n";
    let e = parse_obj(quad).unwrap_err();
    assert_eq!(
        e,
        RsError::PolygonFace {
            face: 1,
            corners: 4
        }
    );
    assert!(e.to_string().contains("face 1 has 4 corners"));
    assert!(matches!(
        parse_obj("v 0 0 0\nv 1 0\nf 1 2 3\n"),
        Err(RsError::Parse { line: 2, .. })
    ));
    assert!(matches!(
        read_obj("/nonexistent/x.obj"),
        Err(RsError::Io { .. })
    ));
    let s4 = Surface::new(3, vec![[0, 1, 2]]).unwrap();
    assert_eq!(
        obj_to_string(&s4, &pos[..2]),
        Err(RsError::PositionCountMismatch {
            expected: 3,
            found: 2
        })
    );
}
