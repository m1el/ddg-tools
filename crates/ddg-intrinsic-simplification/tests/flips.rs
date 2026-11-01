//! The data model, flips with point tracking and Delaunay flipping (SPEC §2, §3.1–3.3; WP1).

mod common;

use common::{Shape, crane_models, dist, planar_grid, reference_meshes, small_shapes};
use ddg_core::generators::SplitMix64;
use ddg_intrinsic_simplification::{
    CoarseMesh, CoarseMeshOptions, IemError, InputPoint, PointLocation,
};

/// Random flips with the certificate after every flip.
fn random_flips(shape: &Shape, mesh: &mut CoarseMesh, n: usize, seed: u64, every: usize) -> usize {
    let mut rng = SplitMix64::new(seed);
    let mut done = 0;
    for step in 0..n {
        let e = (rng.next_u64() % mesh.n_edges() as u64) as usize;
        match mesh.flip_edge(e) {
            Ok(_) => done += 1,
            Err(IemError::FlipRefused { .. }) => {}
            Err(err) => panic!("{}: step {step}: {err}", shape.name),
        }
        if step % every == 0 {
            mesh.certify(false).unwrap_or_else(|err| {
                panic!("{}: after flip {step} of edge {e}: {err}", shape.name)
            });
        }
    }
    mesh.certify(false)
        .unwrap_or_else(|err| panic!("{}: at the end: {err}", shape.name));
    done
}

#[test]
fn construction_is_certified_and_keeps_the_input_ids() {
    for shape in small_shapes().iter().chain(&crane_models()) {
        let points = shape.centroids();
        let mesh = shape.coarse(&points);
        let cert = mesh.certify(false).unwrap();
        assert!(cert.max_angle_error < 1e-12, "{}: {cert:?}", shape.name);
        assert_eq!(mesh.faces(), shape.faces, "{}", shape.name);
        assert_eq!(mesh.n_points(), shape.positions.len() + shape.faces.len());
        for v in 0..shape.positions.len() as u32 {
            assert_eq!(
                mesh.point_location(v as usize).unwrap(),
                PointLocation::Vertex(v)
            );
            assert_eq!(mesh.complex_vertex(v), Some(v as usize));
        }
    }
}

#[test]
fn bad_input_is_an_error_naming_the_element() {
    let shape = planar_grid(3, 0.0, 1);
    let opts = CoarseMeshOptions::default();
    // A degenerate face.
    let mut positions = shape.positions.clone();
    positions[5] = positions[4];
    let err = CoarseMesh::from_positions(&positions, &shape.faces, &[], opts).unwrap_err();
    assert!(matches!(err, IemError::Geometry(_)), "{err}");
    // A non-manifold edge.
    let mut faces = shape.faces.clone();
    faces.push([0, 1, 9]);
    let err = CoarseMesh::from_positions(&shape.positions, &faces, &[], opts).unwrap_err();
    assert!(matches!(err, IemError::Mesh(_)), "{err}");
    // Bad points.
    for (p, what) in [
        (
            InputPoint {
                face: 99,
                bary: [1.0, 0.0, 0.0],
            },
            "out of range",
        ),
        (
            InputPoint {
                face: 0,
                bary: [-0.1, 0.6, 0.5],
            },
            "non-negative",
        ),
        (
            InputPoint {
                face: 0,
                bary: [0.5, 0.6, 0.0],
            },
            "sum",
        ),
        (
            InputPoint {
                face: 0,
                bary: [f64::NAN, 0.5, 0.5],
            },
            "finite",
        ),
    ] {
        let err =
            CoarseMesh::from_positions(&shape.positions, &shape.faces, &[p], opts).unwrap_err();
        match err {
            IemError::InvalidPoint { index: 0, what: w } => assert!(w.contains(what), "{w}"),
            other => panic!("{other}"),
        }
    }
}

/// Random flips on the small shapes keep the certificate after every flip.
#[test]
fn random_flips_on_small_shapes_keep_the_certificate() {
    let mut delta = 0;
    for (i, shape) in small_shapes().iter().enumerate() {
        let points = shape.random_points(10 + i as u64);
        let mut mesh = shape.coarse(&points);
        let flips = random_flips(shape, &mut mesh, 3000, 100 + i as u64, 1);
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
            "{}: {flips} flips, {self_edges} self-edges, {repeated} faces with a repeated vertex",
            shape.name
        );
        delta += self_edges + repeated;
    }
    // Geometric flips rarely make Δ-complexes out of simplicial meshes; the coarsened shapes
    // of `tests/removal.rs` are Δ-complexes and are flipped there.
    eprintln!("Δ-configurations reached by flips alone: {delta}");
}

/// Random flips on the models (Crane's and the reference's), with tracked face centroids.
#[test]
fn random_flips_on_models_keep_the_certificate() {
    for (i, shape) in crane_models()
        .iter()
        .chain(&reference_meshes(12_000))
        .enumerate()
    {
        let points = shape.centroids();
        let mut mesh = shape.coarse(&points);
        let flips = random_flips(shape, &mut mesh, 400, 200 + i as u64, 1);
        eprintln!("{}: {flips} flips certified", shape.name);
    }
}

/// On a planar grid every tracked point keeps its position through any flips: rebuilt from
/// its coarse face and the fine positions of the face's corners.
#[test]
fn flips_on_a_planar_grid_keep_the_points_in_place() {
    for (seed, jitter) in [(1u64, 0.0), (2, 0.3), (3, 0.6)] {
        let shape = planar_grid(12, jitter, seed);
        let points = shape.random_points(seed);
        let mut mesh = shape.coarse(&points);
        random_flips(&shape, &mut mesh, 5000, seed, 50);
        let report = mesh.flip_to_delaunay().unwrap();
        mesh.certify(true).unwrap();
        let mut worst = 0.0f64;
        for (k, p) in points.iter().enumerate() {
            let id = shape.positions.len() + k;
            let got = shape.reconstruct(mesh.point_location(id).unwrap());
            worst = worst.max(dist(got, shape.input_position(p)));
        }
        eprintln!(
            "grid jitter {jitter}: {} Delaunay flips, worst point error {worst:e}, {} clamped",
            report.flips.len(),
            mesh.n_clamped_points()
        );
        assert!(worst < 1e-12, "{worst}");
    }
}

/// Delaunay flipping from all edges reaches a certified Delaunay mesh on the models; tracked
/// points stay valid.
#[test]
fn delaunay_flipping_on_models() {
    for shape in crane_models().iter().chain(&reference_meshes(30_000)) {
        let points = shape.centroids();
        let mut mesh = shape.coarse(&points);
        let before = mesh.certify(false).unwrap();
        let report = mesh.flip_to_delaunay().unwrap();
        let cert = mesh
            .certify(true)
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        eprintln!(
            "{}: {} non-Delaunay edges before, {} flips, {:?}, max excess after {:e}, {} clamped",
            shape.name,
            before.n_not_delaunay,
            report.flips.len(),
            report.status,
            cert.max_delaunay_excess,
            report.n_clamped
        );
    }
}

/// Flips undone by flipping back: the flip of the new edge restores the end points and, up to
/// rounding, the length.
#[test]
fn a_flip_twice_restores_the_edge() {
    let shape = planar_grid(5, 0.3, 4);
    let mut mesh = shape.coarse(&[]);
    for e in 0..mesh.n_edges() {
        let l = mesh.edge_lengths()[e];
        let Ok(f) = mesh.flip_edge(e) else { continue };
        let g = mesh.flip_edge(e).unwrap();
        assert_eq!(g.new, [f.old[1], f.old[0]]);
        assert!((mesh.edge_lengths()[e] - l).abs() <= 1e-14 * l);
        mesh.certify(false).unwrap();
    }
}
