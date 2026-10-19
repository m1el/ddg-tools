//! Tests of small-angle remeshing (SPEC §5.3, §10.6). The mesh is validated after every operation
//! (`validate_each_operation`) and again here.

use ddg_core::generators::{GeneratedMesh, SplitMix64, grid_disk, icosphere, perturb};
use ddg_core::topology::TriMesh;
use ddg_developability::DevError;
use ddg_developability::remesh::{RemeshOptions, remesh};

fn opts() -> RemeshOptions {
    RemeshOptions { validate_each_operation: true, ..Default::default() }
}

fn angles_min(m: &TriMesh, p: &[[f64; 3]]) -> f64 {
    let mut best = f64::INFINITY;
    for f in m.faces() {
        let q = f.map(|v| p[v as usize]);
        for k in 0..3 {
            let (a, b, c) = (q[k], q[(k + 1) % 3], q[(k + 2) % 3]);
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            let n = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
            best = best.min(n.atan2(u[0] * v[0] + u[1] * v[1] + u[2] * v[2]));
        }
    }
    best
}

fn n_bad(m: &TriMesh, p: &[[f64; 3]], small: f64) -> usize {
    m.faces()
        .iter()
        .filter(|f| {
            let t = f.map(|v| p[v as usize]);
            angles_min(&TriMesh::new(3, vec![[0, 1, 2]]).unwrap(), &t) < small
        })
        .count()
}

#[test]
fn constructed_cap_is_flipped() {
    let pos = vec![[0., 0., 0.], [10., 0., 0.], [5., 0.05, 0.], [5., -3., 0.]];
    let m = TriMesh::new(4, vec![[0, 1, 2], [1, 0, 3]]).unwrap();
    let r = remesh(&m, &pos, &opts(), &[]).unwrap();
    assert_eq!((r.report.flips, r.report.collapses), (1, 0));
    assert_eq!(r.report.unresolved_bad_triangles, 0);
    assert!(r.report.min_angle_after > r.report.min_angle_before);
    assert_eq!(r.mesh.n_faces(), 2);
    assert_eq!(r.positions, pos);
    assert!(r.changed());
}

/// Grid vertex `(i, j)` of `grid_disk(n, n)`.
fn gv(n: usize, i: usize, j: usize) -> usize {
    i * (n + 1) + j
}

#[test]
fn constructed_needle_is_collapsed() {
    let n = 5;
    let mut g = grid_disk(n, n).unwrap();
    // interior vertex (2,2) moved next to its neighbour (3,2)
    let (u, w) = (gv(n, 2, 2), gv(n, 3, 2));
    let h = 1.0 / n as f64;
    g.positions[u] = [g.positions[w][0] - 1e-3 * h, g.positions[w][1] + 2e-4 * h, 0.0];
    let m = g.topology().unwrap();
    let before = n_bad(&m, &g.positions, RemeshOptions::default().small_angle);
    assert!(before >= 2);
    let r = remesh(&m, &g.positions, &opts(), &[]).unwrap();
    assert!(r.report.collapses >= 1);
    assert_eq!(r.report.unresolved_bad_triangles, 0);
    assert_eq!(r.mesh.n_vertices(), m.n_vertices() - r.report.collapses);
    assert_eq!(r.mesh.n_faces(), m.n_faces() - 2 * r.report.collapses);
    assert_eq!(r.mesh.boundary_loops().len(), 1);
    assert_eq!(r.mesh.euler_characteristic(), m.euler_characteristic());
    assert_eq!(r.vertex_to_old.len(), r.mesh.n_vertices());
    assert!(r.report.min_angle_after > r.report.min_angle_before);
    // untouched vertices keep their positions
    for (new, &old) in r.vertex_to_old.iter().enumerate() {
        if old as usize != u.min(w) {
            assert_eq!(r.positions[new], g.positions[old as usize]);
        }
    }
}

#[test]
fn boundary_needle_respects_the_option() {
    let n = 5;
    let mut g = grid_disk(n, n).unwrap();
    let (u, w) = (gv(n, 0, 2), gv(n, 0, 3)); // two boundary vertices
    let h = 1.0 / n as f64;
    g.positions[u] = [g.positions[w][0], g.positions[w][1] - 1e-3 * h, 0.0];
    let m = g.topology().unwrap();
    let off = remesh(&m, &g.positions, &RemeshOptions { collapse_boundary_edges: false, ..opts() }, &[]).unwrap();
    assert!(off.report.refused_collapses >= 1);
    assert_eq!(off.report.collapses, 0);
    assert!(off.report.unresolved_bad_triangles >= 1);
    assert_eq!(off.mesh, m);
    let on = remesh(&m, &g.positions, &opts(), &[]).unwrap();
    assert!(on.report.collapses >= 1);
    assert_eq!(on.report.unresolved_bad_triangles, 0);
    assert_eq!(on.mesh.boundary_loops().len(), 1);
    assert_eq!(on.mesh.n_boundary_edges(), m.n_boundary_edges() - on.report.collapses.min(1));
}

#[test]
fn good_mesh_is_unchanged() {
    let g = icosphere(2).unwrap();
    let m = g.topology().unwrap();
    let r = remesh(&m, &g.positions, &opts(), &[3]).unwrap();
    assert!(!r.changed());
    assert_eq!(r.mesh, m);
    assert_eq!(r.positions, g.positions);
    assert_eq!(r.vertex_to_old, (0..m.n_vertices() as u32).collect::<Vec<_>>());
    assert_eq!(r.pinned, vec![3]);
    assert_eq!(r.report.min_angle_before, r.report.min_angle_after);
}

/// Moves random vertices next to a neighbour (needles) or onto the segment between two
/// neighbours of a neighbour (caps).
fn inject(g: &mut GeneratedMesh, m: &TriMesh, count: usize, seed: u64) {
    let mut rng = SplitMix64::new(seed);
    for _ in 0..count {
        let v = (rng.next_u64() % m.n_vertices() as u64) as usize;
        let nb: Vec<usize> = m.vertex_neighbors(v).collect();
        let w = nb[(rng.next_u64() % nb.len() as u64) as usize];
        let t = rng.uniform(0.0, 0.01);
        if rng.next_f64() < 0.5 {
            for k in 0..3 {
                g.positions[v][k] = (1.0 - t) * g.positions[w][k] + t * g.positions[v][k];
            }
        } else {
            let w2 = nb[(rng.next_u64() % nb.len() as u64) as usize];
            for k in 0..3 {
                g.positions[v][k] = (0.5 - t) * g.positions[w][k] + (0.5 + t) * g.positions[w2][k];
            }
        }
    }
}

#[test]
fn random_meshes_with_injected_bad_triangles_stay_valid() {
    let mut total = (0, 0, 0, 0);
    for seed in 0..40u64 {
        let mut g = match seed % 3 {
            0 => icosphere(2).unwrap(),
            1 => {
                let mut g = grid_disk(8, 7).unwrap();
                perturb(&mut g.positions, [0.02, 0.02, 0.05], seed);
                g
            }
            _ => {
                let mut g = icosphere(1).unwrap();
                perturb(&mut g.positions, [0.1, 0.1, 0.1], seed);
                g
            }
        };
        let m = g.topology().unwrap();
        inject(&mut g, &m, 8, 100 + seed);
        let small = RemeshOptions::default().small_angle;
        let before = n_bad(&m, &g.positions, small);
        let r = remesh(&m, &g.positions, &opts(), &[]).unwrap_or_else(|e| panic!("seed {seed}: {e}"));
        // topology invariants
        assert_eq!(r.mesh.euler_characteristic(), m.euler_characteristic(), "seed {seed}");
        assert_eq!(r.mesh.boundary_loops().len(), m.boundary_loops().len(), "seed {seed}");
        assert_eq!(r.mesh.n_components(), m.n_components());
        assert_eq!(r.mesh.n_vertices(), m.n_vertices() - r.report.collapses);
        // each collapse removes two faces, or one if the edge is on the boundary
        let c = r.report.collapses;
        if m.n_boundary_edges() == 0 {
            assert_eq!(r.mesh.n_faces(), m.n_faces() - 2 * c);
        } else {
            assert!(r.mesh.n_faces() + c <= m.n_faces() && m.n_faces() - 2 * c <= r.mesh.n_faces());
        }
        assert!(r.positions.iter().all(|p| p.iter().all(|c| c.is_finite())));
        let after = n_bad(&r.mesh, &r.positions, small);
        assert_eq!(after, r.report.unresolved_bad_triangles);
        assert!(after <= before, "seed {seed}: {before} -> {after}");
        assert!((angles_min(&r.mesh, &r.positions) - r.report.min_angle_after).abs() < 1e-12);
        total.0 += before;
        total.1 += after;
        total.2 += r.report.flips;
        total.3 += r.report.collapses;
    }
    eprintln!("bad before {} after {}, flips {}, collapses {}", total.0, total.1, total.2, total.3);
    assert!(total.0 > 0 && total.1 * 4 < total.0 && total.2 > 0 && total.3 > 0);
}

#[test]
fn pinned_vertices_never_move_and_follow_merges() {
    let n = 5;
    let mut g = grid_disk(n, n).unwrap();
    let (u, w) = (gv(n, 2, 2), gv(n, 3, 2));
    let h = 1.0 / n as f64;
    g.positions[u] = [g.positions[w][0] - 1e-3 * h, g.positions[w][1] + 2e-4 * h, 0.0];
    let m = g.topology().unwrap();
    // pin the higher id of the short edge: it must survive with its position
    let hi = u.max(w);
    let r = remesh(&m, &g.positions, &opts(), &[hi as u32, 0]).unwrap();
    assert!(r.report.collapses >= 1);
    for &p in &r.pinned {
        let old = r.vertex_to_old[p as usize] as usize;
        assert_eq!(r.positions[p as usize], g.positions[old]);
    }
    assert!(r.pinned.iter().any(|&p| r.vertex_to_old[p as usize] as usize == hi));
    assert_eq!(r.pinned.len(), 2);
    // both ends pinned: refused
    let r = remesh(&m, &g.positions, &opts(), &[u as u32, w as u32]).unwrap();
    assert!(r.report.refused_collapses >= 1);
    assert_eq!(r.positions[u], g.positions[u]);
}

#[test]
fn bad_input_and_options() {
    let g = icosphere(1).unwrap();
    let m = g.topology().unwrap();
    let bad = |o: RemeshOptions| matches!(remesh(&m, &g.positions, &o, &[]), Err(DevError::InvalidOption { .. }));
    assert!(bad(RemeshOptions { small_angle: 0.0, ..Default::default() }));
    assert!(bad(RemeshOptions { small_angle: std::f64::consts::FRAC_PI_3, ..Default::default() }));
    assert!(bad(RemeshOptions { small_angle: f64::NAN, ..Default::default() }));
    assert!(bad(RemeshOptions { max_passes: 0, ..Default::default() }));
    assert!(bad(RemeshOptions { check_interval: 0, ..Default::default() }));
    assert!(matches!(remesh(&m, &g.positions[1..], &opts(), &[]), Err(DevError::LengthMismatch { .. })));
    assert!(matches!(
        remesh(&m, &g.positions, &opts(), &[m.n_vertices() as u32]),
        Err(DevError::VertexOutOfRange { .. })
    ));
    let mut p = g.positions.clone();
    p[2][0] = f64::INFINITY;
    assert!(matches!(remesh(&m, &p, &opts(), &[]), Err(DevError::NonFinitePosition { vertex: 2 })));
}
