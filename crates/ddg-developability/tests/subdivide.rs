//! Tests of 4-1 midpoint subdivision (SPEC §5.2, §10.6).

use ddg_core::generators::{annulus, bumpy_grid_disk, grid_disk, icosphere, perturb, torus, GeneratedMesh};
use ddg_core::topology::TriMesh;
use ddg_developability::subdivide::subdivide;
use ddg_developability::DevError;

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn normal(p: &[[f64; 3]], f: [u32; 3]) -> [f64; 3] {
    cross(sub(p[f[1] as usize], p[f[0] as usize]), sub(p[f[2] as usize], p[f[0] as usize]))
}

fn meshes() -> Vec<(&'static str, GeneratedMesh)> {
    let mut sphere = icosphere(1).unwrap();
    perturb(&mut sphere.positions, [0.05, 0.05, 0.05], 3);
    vec![
        ("icosahedron", icosphere(0).unwrap()),
        ("sphere", sphere),
        ("torus", torus(12, 7, 2.0, 0.7).unwrap()),
        ("disk", grid_disk(5, 4).unwrap()),
        ("bumpy disk", bumpy_grid_disk(6, 6, 0.3).unwrap()),
        ("annulus", annulus(9, 3, 0.5, 1.0).unwrap()),
    ]
}

#[test]
fn counts_topology_and_geometry() {
    for (name, g) in meshes() {
        let m: TriMesh = g.topology().unwrap();
        let s = subdivide(&m, &g.positions).unwrap();
        let (nv, ne, nf) = (m.n_vertices(), m.n_edges(), m.n_faces());
        assert_eq!(s.mesh.n_faces(), 4 * nf, "{name}");
        assert_eq!(s.mesh.n_vertices(), nv + ne, "{name}");
        assert_eq!(s.positions.len(), nv + ne);
        assert_eq!(s.mesh.euler_characteristic(), m.euler_characteristic(), "{name}");
        assert_eq!(s.mesh.boundary_loops().len(), m.boundary_loops().len(), "{name}");
        assert_eq!(s.mesh.n_boundary_edges(), 2 * m.n_boundary_edges(), "{name}");
        assert_eq!(s.mesh.n_components(), m.n_components(), "{name}");
        assert_eq!(&s.positions[..nv], &g.positions[..]);
        // midpoints
        for (e, &[a, b]) in m.edges().iter().enumerate() {
            let q = s.positions[s.edge_midpoint[e] as usize];
            let (pa, pb) = (g.positions[a as usize], g.positions[b as usize]);
            assert_eq!(q, [0, 1, 2].map(|k| 0.5 * (pa[k] + pb[k])));
        }
        // children are coplanar with their parent, same orientation, a quarter of the area
        for (t, f) in s.mesh.faces().iter().enumerate() {
            let parent = m.faces()[s.face_parent[t] as usize];
            let np = normal(&g.positions, parent);
            let nc = normal(&s.positions, *f);
            let cr = cross(np, nc);
            let scale = dot(np, np).sqrt() * dot(nc, nc).sqrt();
            assert!(dot(cr, cr).sqrt() <= 1e-12 * scale, "{name}: child {t} not parallel");
            assert!((dot(np, nc) / (dot(np, np)) - 0.25).abs() < 1e-12, "{name}: child {t} area/orientation");
        }
        // valences: old unchanged, interior midpoints 6, boundary midpoints 3
        for v in 0..nv {
            assert_eq!(s.mesh.vertex_face_count(v), m.vertex_face_count(v));
        }
        for e in 0..m.n_edges() {
            let v = s.edge_midpoint[e] as usize;
            let expect = if m.is_boundary_edge(e) { 3 } else { 6 };
            assert_eq!(s.mesh.vertex_face_count(v), expect, "{name}");
            assert_eq!(s.mesh.is_boundary_vertex(v), m.is_boundary_edge(e));
        }
    }
}

#[test]
fn repeated_subdivision_stays_valid() {
    let g = bumpy_grid_disk(4, 4, 0.2).unwrap();
    let mut m = g.topology().unwrap();
    let mut p = g.positions.clone();
    for _ in 0..3 {
        let s = subdivide(&m, &p).unwrap();
        m = s.mesh;
        p = s.positions;
    }
    assert_eq!(m.n_faces(), 64 * g.faces.len());
    assert_eq!(m.boundary_loops().len(), 1);
}

#[test]
fn pinned_ids_follow() {
    let g = grid_disk(3, 3).unwrap();
    let m = g.topology().unwrap();
    let s = subdivide(&m, &g.positions).unwrap();
    let e = m.find_edge(0, m.vertex_neighbors(0).next().unwrap() as usize).unwrap();
    let [a, b] = m.edge_vertices(e);
    let pinned = s.pinned_after(&m, &[a as u32, b as u32]).unwrap();
    assert_eq!(pinned, vec![a.min(b) as u32, a.max(b) as u32, s.edge_midpoint[e]]);
    let one = s.pinned_after(&m, &[a as u32]).unwrap();
    assert_eq!(one, vec![a as u32]);
    assert!(matches!(
        s.pinned_after(&m, &[m.n_vertices() as u32]),
        Err(DevError::VertexOutOfRange { .. })
    ));
}

#[test]
fn bad_input_is_an_error() {
    let g = icosphere(0).unwrap();
    let m = g.topology().unwrap();
    assert!(matches!(subdivide(&m, &g.positions[1..]), Err(DevError::LengthMismatch { .. })));
    let mut p = g.positions.clone();
    p[3][1] = f64::NAN;
    assert!(matches!(subdivide(&m, &p), Err(DevError::NonFinitePosition { vertex: 3 })));
}

#[test]
fn all_six_energies_are_preserved() {
    use ddg_developability::{Energy, EnergyOptions, evaluate_energy};
    let mut sphere = icosphere(2).unwrap();
    perturb(&mut sphere.positions, [0.03, 0.03, 0.03], 9);
    let mut disk = bumpy_grid_disk(6, 5, 0.15).unwrap();
    perturb(&mut disk.positions, [0.0, 0.0, 0.02], 2);
    let mut tor = torus(10, 6, 2.0, 0.7).unwrap();
    perturb(&mut tor.positions, [0.02, 0.02, 0.02], 5);
    for (name, g) in [("sphere", sphere), ("disk", disk), ("torus", tor)] {
        let m = g.topology().unwrap();
        let s = subdivide(&m, &g.positions).unwrap();
        for energy in [
            Energy::Combinatorial,
            Energy::CombinatorialMax,
            Energy::Covariance,
            Energy::CovarianceMax,
            Energy::Intrinsic,
            Energy::IntrinsicMax,
        ] {
            let o = EnergyOptions { energy, ..Default::default() };
            let (e0, _) = evaluate_energy(&m, &g.positions, &o).unwrap();
            let (e1, _) = evaluate_energy(&s.mesh, &s.positions, &o).unwrap();
            assert!(e0 > 1e-3, "{name} {energy:?}: test energy too small ({e0})");
            assert!((e0 - e1).abs() <= 1e-12 * e0.max(1.0), "{name} {energy:?}: {e0} vs {e1}");
        }
    }
}
