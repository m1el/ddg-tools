//! Flattening, vertex removal and the undo log (SPEC §3.4–3.6; WP2).

mod common;

use std::f64::consts::PI;

use common::{Shape, crane_models, dist, planar_grid, reference_meshes, small_shapes};
use ddg_core::generators::SplitMix64;
use ddg_intrinsic_simplification::{
    CoarseMesh, IemError, RemovalKind, RemovalRefusal, RemovalReport,
};

#[derive(Debug, Default)]
struct Tally {
    removed: [usize; 3],
    skipped: std::collections::BTreeMap<String, usize>,
    fallback_flips: usize,
    degree_flips: usize,
    delaunay_flips: usize,
    max_iterations: usize,
    max_conservation_error: f64,
    clamped: usize,
    /// Interior edges that fail the Delaunay test at the end.
    not_delaunay: usize,
}

/// Checks one removal: curvature conservation (`iem-flatten-conserves-curvature`), the flat
/// star.
fn check_removal(shape: &Shape, r: &RemovalReport, t: &mut Tally) {
    let p = &r.prepared;
    let before: f64 = p.neighbours.iter().map(|n| n.curvature_before).sum();
    let after: f64 = p.neighbours.iter().map(|n| n.curvature_after).sum();
    let scale = 1.0
        + p.curvature.abs()
        + p.neighbours
            .iter()
            .map(|n| n.curvature_before.abs())
            .sum::<f64>();
    let err = (after - before - p.curvature).abs() / scale;
    assert!(
        err < 1e-9,
        "{}: vertex {}: curvature not conserved: {before} + {} != {after}",
        shape.name,
        p.vertex,
        p.curvature
    );
    t.max_conservation_error = t.max_conservation_error.max(err);
    t.removed[match p.kind {
        RemovalKind::Interior => 0,
        RemovalKind::Boundary => 1,
        RemovalKind::Ear => 2,
    }] += 1;
    t.fallback_flips += p.fallback_flips.len();
    t.degree_flips += p.degree_flips.len();
    t.delaunay_flips += r.delaunay.as_ref().map_or(0, |d| d.flips.len());
    t.max_iterations = t.max_iterations.max(p.flatten.iterations);
    t.clamped += p.n_clamped + r.delaunay.as_ref().map_or(0, |d| d.n_clamped);
}

/// Removes every removable vertex in passes by ascending fine id until a pass removes nothing;
/// certifies every `every` removals.
fn remove_all(shape: &Shape, mesh: &mut CoarseMesh, every: usize) -> Tally {
    let mut t = Tally::default();
    let mut n = 0usize;
    loop {
        let mut progress = false;
        for v in 0..mesh.n_fine_vertices() as u32 {
            if !mesh.is_vertex(v) {
                continue;
            }
            match mesh.remove_vertex(v) {
                Ok(r) => {
                    check_removal(shape, &r, &mut t);
                    progress = true;
                    n += 1;
                    if n.is_multiple_of(every) {
                        mesh.certify(false).unwrap_or_else(|e| {
                            panic!("{}: after removing {v} ({n} removals): {e}", shape.name)
                        });
                    }
                }
                Err(IemError::Removal(f)) => {
                    if let RemovalRefusal::Internal { what } = &f.reason {
                        panic!("{}: {what}", shape.name);
                    }
                    let key = format!("{:?}", f.reason)
                        .split([' ', '{', '('])
                        .next()
                        .unwrap_or("")
                        .to_string();
                    *t.skipped.entry(key).or_default() += 1;
                }
                Err(e) => panic!("{}: removal of {v}: {e}", shape.name),
            }
        }
        if !progress {
            break;
        }
    }
    let cert = mesh
        .certify(false)
        .unwrap_or_else(|e| panic!("{}: at the end: {e}", shape.name));
    t.not_delaunay = cert.n_not_delaunay;
    t
}

#[test]
fn removing_every_removable_vertex_of_the_small_shapes() {
    for shape in small_shapes() {
        let points = shape.random_points(3);
        let mut mesh = shape.coarse(&points);
        mesh.flip_to_delaunay().unwrap();
        let t = remove_all(&shape, &mut mesh, 1);
        eprintln!(
            "{}: {} -> {} vertices, removed {:?} (interior, boundary, ear), skips {:?}, \
             fallback {}, degree flips {}, Delaunay flips {}, max flatten steps {}, \
             conservation {:e}, clamped {}",
            shape.name,
            shape.positions.len(),
            mesh.n_vertices(),
            t.removed,
            t.skipped,
            t.fallback_flips,
            t.degree_flips,
            t.delaunay_flips,
            t.max_iterations,
            t.max_conservation_error,
            t.clamped
        );
    }
}

#[test]
fn removing_every_removable_vertex_of_the_models() {
    for shape in crane_models().iter().chain(&reference_meshes(3_000)) {
        let points = shape.centroids();
        let mut mesh = shape.coarse(&points);
        mesh.flip_to_delaunay().unwrap();
        let every = if shape.positions.len() > 1000 { 7 } else { 1 };
        let t = remove_all(shape, &mut mesh, every);
        eprintln!(
            "{}: {} -> {} vertices ({} faces), removed {:?} (interior, boundary, ear), skips \
             {:?}, fallback {}, degree flips {}, Delaunay flips {}, max flatten steps {}, \
             conservation {:e}, clamped {}, not Delaunay at the end {}",
            shape.name,
            shape.positions.len(),
            mesh.n_vertices(),
            mesh.n_faces(),
            t.removed,
            t.skipped,
            t.fallback_flips,
            t.degree_flips,
            t.delaunay_flips,
            t.max_iterations,
            t.max_conservation_error,
            t.clamped,
            t.not_delaunay
        );
        if shape.name == "CrumpledDevelopable" {
            assert!(t.removed[2] > 0, "no ear was removed");
        }
    }
}

/// The larger reference meshes (capsule, spot_subdiv, bs_rest, dragon_fat; up to 25 000
/// vertices), every removable vertex by id, certified every 101 removals. Slow in a debug
/// build: `cargo test --release -- --ignored`.
#[test]
#[ignore]
fn removing_every_removable_vertex_of_the_large_reference_meshes() {
    for shape in reference_meshes(30_000)
        .iter()
        .filter(|s| s.positions.len() > 3_000)
    {
        let mut mesh = shape.coarse(&[]);
        mesh.flip_to_delaunay().unwrap();
        let t = remove_all(shape, &mut mesh, 101);
        eprintln!(
            "{}: {} -> {} vertices, removed {:?}, skips {:?}, fallback {}, degree flips {}, \
             Delaunay flips {}, max flatten steps {}, conservation {:e}, clamped {}, not \
             Delaunay at the end {}",
            shape.name,
            shape.positions.len(),
            mesh.n_vertices(),
            t.removed,
            t.skipped,
            t.fallback_flips,
            t.degree_flips,
            t.delaunay_flips,
            t.max_iterations,
            t.max_conservation_error,
            t.clamped,
            t.not_delaunay
        );
    }
}

/// A planar grid with its four corners fixed: every removal flattens nothing (`u = 0`), so the
/// coarse metric stays the planar one and the tracked map reproduces the positions of every
/// fine vertex and extra point (SPEC §12.2; the reference reaches 1.9e-8).
#[test]
fn removals_on_a_planar_grid_keep_the_points_in_place() {
    for (n, jitter, seed) in [(8usize, 0.0, 1u64), (16, 0.3, 2), (24, 0.6, 3)] {
        let shape = planar_grid(n, jitter, seed);
        let points = shape.random_points(seed);
        let mut mesh = shape.coarse(&points);
        for corner in [0, n, n * (n + 1), (n + 1) * (n + 1) - 1] {
            mesh.set_fixed(corner as u32, true).unwrap();
        }
        mesh.flip_to_delaunay().unwrap();
        let t = remove_all(&shape, &mut mesh, 1);
        let mut worst = 0.0f64;
        for v in 0..shape.positions.len() {
            let got = shape.reconstruct(mesh.point_location(v).unwrap());
            worst = worst.max(dist(got, shape.positions[v]));
        }
        for (k, p) in points.iter().enumerate() {
            let got = shape.reconstruct(mesh.point_location(shape.positions.len() + k).unwrap());
            worst = worst.max(dist(got, shape.input_position(p)));
        }
        eprintln!(
            "grid {n} jitter {jitter}: {} -> {} vertices, removed {:?}, skips {:?}, worst point \
             error {worst:e}, clamped {}",
            shape.positions.len(),
            mesh.n_vertices(),
            t.removed,
            t.skipped,
            t.clamped
        );
        assert!(worst < 1e-9, "{worst}");
        assert_eq!(mesh.n_vertices(), 4, "only the fixed corners should remain");
    }
}

/// Reverting a removal (after the first phase, with or without tracking) restores the mesh
/// bitwise; a tentative removal computes what the committed removal computes (S5).
#[test]
fn reverting_restores_the_mesh_bitwise() {
    let mut checked = 0;
    let mut models = small_shapes();
    models.extend(crane_models());
    for (i, shape) in models.iter().enumerate() {
        let points = shape.random_points(11);
        let mut mesh = shape.coarse(&points);
        mesh.flip_to_delaunay().unwrap();
        let mut rng = SplitMix64::new(500 + i as u64);
        for step in 0..300 {
            let v = (rng.next_u64() % mesh.n_fine_vertices() as u64) as u32;
            if !mesh.is_vertex(v) {
                continue;
            }
            let before = mesh.clone();
            let tentative = mesh.tentative_removal(v);
            assert!(mesh.bitwise_eq(&before), "{}: tentative {v}", shape.name);
            match mesh.begin_removal(v) {
                Ok(prepared) => {
                    let mut t = tentative.expect("tentative failed where the removal works");
                    t.n_clamped = prepared.n_clamped;
                    assert_eq!(t, prepared, "{}: vertex {v}", shape.name);
                    mesh.revert_removal().unwrap();
                    assert!(mesh.bitwise_eq(&before), "{}: revert of {v}", shape.name);
                    checked += 1;
                    // Commit every third one, so that the reverts happen on coarsened meshes.
                    if step % 3 == 0 {
                        mesh.remove_vertex(v).unwrap();
                    }
                }
                Err(f) => {
                    assert_eq!(tentative.unwrap_err(), f);
                    assert!(mesh.bitwise_eq(&before), "{}: failed {v}", shape.name);
                }
            }
        }
    }
    eprintln!("{checked} reverts bitwise exact");
    assert!(checked > 500);
}

/// Removal refusals: fixed vertices, the vertices left on a coarse torus, and the mesh is
/// unchanged by each.
#[test]
fn refusals_change_nothing() {
    let shape = planar_grid(4, 0.2, 9);
    let mut mesh = shape.coarse(&[]);
    mesh.set_fixed(6, true).unwrap();
    let before = mesh.clone();
    match mesh.remove_vertex(6) {
        Err(IemError::Removal(f)) => assert_eq!(f.reason, RemovalRefusal::Fixed),
        other => panic!("{other:?}"),
    }
    assert!(mesh.bitwise_eq(&before));
    mesh.set_fixed(6, false).unwrap();
    mesh.remove_vertex(6).unwrap();
    match mesh.remove_vertex(6) {
        Err(IemError::Removal(f)) => assert_eq!(f.reason, RemovalRefusal::NotAVertex),
        other => panic!("{other:?}"),
    }
    assert!(matches!(
        mesh.revert_removal(),
        Err(IemError::NoRemovalPending)
    ));
    mesh.begin_removal(7).unwrap();
    assert!(matches!(mesh.flip_edge(0), Err(IemError::RemovalPending)));
    mesh.revert_removal().unwrap();
}

/// Coarsened small shapes are Δ-complexes; random flips and removals on them keep the
/// certificate.
#[test]
fn delta_complexes_from_coarsening() {
    let mut delta = 0;
    for (i, shape) in small_shapes().iter().enumerate() {
        let points = shape.random_points(21);
        let mut mesh = shape.coarse(&points);
        remove_all(shape, &mut mesh, 1);
        let mut rng = SplitMix64::new(700 + i as u64);
        for _ in 0..200 {
            let e = (rng.next_u64() % mesh.n_edges() as u64) as usize;
            match mesh.flip_edge(e) {
                Ok(_) | Err(IemError::FlipRefused { .. }) => {}
                Err(err) => panic!("{}: {err}", shape.name),
            }
            mesh.certify(false).unwrap();
        }
        let c = mesh.complex();
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
        eprintln!(
            "{}: {} vertices left, {self_edges} self-edges, {repeated} faces with a repeated \
             vertex",
            shape.name,
            mesh.n_vertices()
        );
        delta += self_edges + repeated;
    }
    assert!(delta > 0);
}

#[test]
fn flattening_reaches_the_target_angle_sum() {
    for shape in crane_models().iter().take(2) {
        let mut mesh = shape.coarse(&[]);
        mesh.flip_to_delaunay().unwrap();
        for v in (0..mesh.n_fine_vertices() as u32).step_by(37) {
            let Ok(p) = mesh.begin_removal(v) else {
                continue;
            };
            let c = mesh.complex_vertex(v).unwrap();
            let target = if mesh.complex().is_boundary_vertex(c) {
                PI
            } else {
                2.0 * PI
            };
            // The solve stops within its tolerance; the flips to degree three that follow keep
            // the angle sum up to the rounding of their new corner angles.
            assert!(
                p.flatten.residual <= 1e-12 * target,
                "{}: {v}: {:?}",
                shape.name,
                p.flatten
            );
            let theta = mesh.angle_sum(c);
            assert!(
                (theta - target).abs() <= 1e-9,
                "{}: {v}: {theta}",
                shape.name
            );
            assert!(p.flatten.iterations <= 20, "{:?}", p.flatten);
            mesh.revert_removal().unwrap();
        }
    }
}
