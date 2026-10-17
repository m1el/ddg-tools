//! The heat method on intrinsic triangulations (paper Sec. 5.3).

mod common;

use common::{Shape, dist, jittered_grid, position, sheared_torus, vertex_positions};
use ddg_core::generators::{grid_disk, icosphere};
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation, RefineOptions, SurfacePoint};
use ddg_intrinsic_apps::exact::{ExactOptions, exact_distance};
use ddg_intrinsic_apps::heat::{
    HeatBoundary, HeatError, HeatOptions, HeatSolver, HeatSource, geodesic_distance,
};

fn idt(shape: &Shape) -> IntrinsicTriangulation {
    let mut tri = shape.triangulation();
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    tri
}

fn idr(shape: &Shape) -> IntrinsicTriangulation {
    let mut tri = idt(shape);
    tri.refine(&RefineOptions::default()).unwrap();
    tri
}

/// Mean of `|d − e| / e` over the vertices with `e > 0`.
fn mean_relative_error(d: &[f64], e: &[f64]) -> f64 {
    let pairs: Vec<f64> = d
        .iter()
        .zip(e)
        .filter(|(_, e)| **e > 0.0)
        .map(|(d, e)| (d - e).abs() / e)
        .collect();
    pairs.iter().sum::<f64>() / pairs.len() as f64
}

#[test]
fn on_a_flat_square_the_distance_is_nearly_euclidean() {
    let shape = jittered_grid(24, 0.3, 11);
    for (mode, tri) in [("idt", idt(&shape)), ("idr", idr(&shape))] {
        let positions = vertex_positions(&shape, &tri);
        let source = 12 * 25 + 12; // the centre vertex
        let d = geodesic_distance(&tri, &[HeatSource::Vertex(source)], &HeatOptions::default())
            .unwrap();
        assert_eq!(d.report.n_flat_faces, 0);
        assert!(d.report.inconsistency < 1e-10, "{:?}", d.report);
        assert!(d.distances[source].abs() < 1e-12);
        let euclid: Vec<f64> = positions
            .iter()
            .map(|&p| dist(p, positions[source]))
            .collect();
        let err = mean_relative_error(&d.distances, &euclid);
        assert!(err < 0.03, "{mode}: mean relative error {err}");
        // Far from the source the error is small in absolute terms too.
        let worst = d
            .distances
            .iter()
            .zip(&euclid)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(worst < 0.03, "{mode}: largest error {worst}");
    }
}

#[test]
fn symmetric_input_gives_symmetric_distance() {
    // A grid symmetric about the vertical line through its centre column, whose triangles are
    // mirrored too: the distance from a vertex on the axis is symmetric to rounding.
    let n = 10;
    let g = grid_disk(n, n).unwrap();
    let mut faces = Vec::new();
    let id = |i: usize, j: usize| (i * (n + 1) + j) as u32;
    for i in 0..n {
        for j in 0..n {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            if (i < n / 2) == (j % 2 == 0) {
                faces.push([a, b, c]);
                faces.push([a, c, d]);
            } else {
                faces.push([a, b, d]);
                faces.push([b, c, d]);
            }
        }
    }
    let shape = Shape::new("mirrored grid", g.positions, faces);
    let tri = shape.triangulation();
    let source = id(n / 2, 3) as usize;
    let d = geodesic_distance(&tri, &[HeatSource::Vertex(source)], &HeatOptions::default())
        .unwrap()
        .distances;
    for i in 0..=n {
        for j in 0..=n {
            let (a, b) = (d[id(i, j) as usize], d[id(n - i, j) as usize]);
            assert!((a - b).abs() < 1e-10, "({i}, {j}): {a} != {b}");
        }
    }
    // Swapping source and target: not exact for the heat method, close on a good mesh.
    let sphere = Shape::generated("icosphere 3", icosphere(3).unwrap());
    let solver = HeatSolver::new(&sphere.triangulation(), &HeatOptions::default()).unwrap();
    let (a, b) = (0, 300);
    let ab = solver.distance(&[HeatSource::Vertex(a)]).unwrap().distances[b];
    let ba = solver.distance(&[HeatSource::Vertex(b)]).unwrap().distances[a];
    assert!((ab - ba).abs() < 0.01 * ab, "{ab} {ba}");
}

#[test]
fn the_intrinsic_delaunay_triangulation_and_refinement_reduce_the_error() {
    // Paper Fig. 18 on a closed mesh of long obtuse triangles: input > iDT > iDR.
    let shape = sheared_torus(60, 16, 1.0, 0.45, 5);
    let source = 0;
    let exact = exact_distance(&shape.mesh, &shape.geometry, &[source], &ExactOptions::default())
        .unwrap()
        .distances;
    let n = shape.mesh.n_vertices();
    let error = |tri: &IntrinsicTriangulation| {
        let d = geodesic_distance(tri, &[HeatSource::Vertex(source)], &HeatOptions::default())
            .unwrap()
            .distances;
        mean_relative_error(&d[..n], &exact)
    };
    // Refinement with a bound on the circumradius (half the mean edge of the iDT): an angle
    // bound alone does not make the triangles smaller (ddg-intrinsic SPEC N16).
    let mut fine = idt(&shape);
    let mean_edge = fine.edge_lengths().iter().sum::<f64>() / fine.n_edges() as f64;
    fine.refine(&RefineOptions {
        max_circumradius: Some(0.5 * mean_edge),
        ..RefineOptions::default()
    })
    .unwrap();
    let (e_input, e_idt, e_idr, e_fine) = (
        error(&shape.triangulation()),
        error(&idt(&shape)),
        error(&idr(&shape)),
        error(&fine),
    );
    println!("mean relative error: input {e_input}, iDT {e_idt}, iDR {e_idr}, iDR/r {e_fine}");
    assert!(e_input > 5.0 * e_idt, "{e_input} {e_idt}");
    assert!(e_idr < 1.1 * e_idt, "{e_idt} {e_idr}");
    assert!(e_fine < 0.8 * e_idt, "{e_idt} {e_fine}");
    assert!(e_fine < 0.02, "{e_fine}");
}

#[test]
fn factorizations_are_reused_and_sources_can_be_points() {
    let shape = jittered_grid(12, 0.3, 3);
    let tri = idr(&shape);
    let solver = HeatSolver::new(&tri, &HeatOptions::default()).unwrap();
    assert!(solver.setup().time > 0.0);
    let first = solver.distance(&[HeatSource::Vertex(20)]).unwrap();
    let _other = solver.distance(&[HeatSource::Vertex(70)]).unwrap();
    let again = solver.distance(&[HeatSource::Vertex(20)]).unwrap();
    assert_eq!(first.distances, again.distances);
    let fresh = geodesic_distance(&tri, &[HeatSource::Vertex(20)], &HeatOptions::default())
        .unwrap();
    assert_eq!(first.distances, fresh.distances);
    // An input vertex as a point of the input mesh.
    let as_point = solver
        .distance(&[HeatSource::Point(SurfacePoint::Vertex(20))])
        .unwrap();
    assert_eq!(first.distances, as_point.distances);
    // A point inside an input face. The shift makes the value interpolated at the source the
    // barycentric mean of the distances of the corners of the face of the triangulation it
    // lies in: positive and below the longest edge. At the corners of the input face the
    // value is close to the straight distance.
    let face = 37;
    let bary = [0.2, 0.5, 0.3];
    let point = SurfacePoint::Face { face, bary };
    let d = solver.distance(&[HeatSource::Point(point)]).unwrap();
    let at = solver.value_at(&d.distances, point).unwrap();
    let longest = tri.edge_lengths().iter().fold(0.0f64, |m, &x| m.max(x));
    assert!(at > 0.0 && at < longest, "value at the source point: {at}");
    let p = position(&shape, point);
    for &v in &shape.mesh.faces()[face as usize] {
        let e = dist(p, shape.positions[v as usize]);
        assert!((d.distances[v as usize] - e).abs() < 0.5 * longest, "{}", d.distances[v as usize]);
    }
    // Values at inserted vertices: their place on the input mesh, queried, gives the value.
    let positions = vertex_positions(&shape, &tri);
    for v in shape.mesh.n_vertices()..tri.n_vertices() {
        let value = solver
            .value_at(&first.distances, tri.vertex_location(v))
            .unwrap();
        assert!((value - first.distances[v]).abs() < 1e-9, "vertex {v}");
        let e = dist(positions[v], positions[20]);
        assert!((first.distances[v] - e).abs() < 0.05, "vertex {v}");
    }
    // Several sources: the distance to the nearest (smoothed where the two meet).
    let two = solver
        .distance(&[HeatSource::Vertex(20), HeatSource::Vertex(150)])
        .unwrap();
    let far = solver.distance(&[HeatSource::Vertex(150)]).unwrap();
    for v in 0..tri.n_vertices() {
        let nearest = first.distances[v].min(far.distances[v]);
        assert!((two.distances[v] - nearest).abs() < 0.1 * nearest + 0.02, "vertex {v}: {} {nearest}", two.distances[v]);
    }
}

#[test]
fn boundary_conditions() {
    let shape = jittered_grid(16, 0.25, 5);
    let tri = idt(&shape);
    let source = 8 * 17 + 8;
    let positions = vertex_positions(&shape, &tri);
    let euclid: Vec<f64> = positions
        .iter()
        .map(|&p| dist(p, positions[source]))
        .collect();
    for boundary in [HeatBoundary::Neumann, HeatBoundary::Dirichlet, HeatBoundary::Average] {
        let options = HeatOptions {
            boundary,
            ..HeatOptions::default()
        };
        let d = geodesic_distance(&tri, &[HeatSource::Vertex(source)], &options).unwrap();
        let err = mean_relative_error(&d.distances, &euclid);
        assert!(err < 0.05, "{boundary:?}: {err}");
    }
    let dirichlet = HeatOptions {
        boundary: HeatBoundary::Dirichlet,
        ..HeatOptions::default()
    };
    assert_eq!(
        geodesic_distance(&tri, &[HeatSource::Vertex(0)], &dirichlet).map(|_| ()),
        Err(HeatError::SourceOnDirichletBoundary { vertex: 0 })
    );
}

#[test]
fn bad_input_is_an_error() {
    let shape = Shape::generated("icosphere 1", icosphere(1).unwrap());
    let tri = shape.triangulation();
    let solver = HeatSolver::new(&tri, &HeatOptions::default()).unwrap();
    assert_eq!(solver.distance(&[]).map(|_| ()), Err(HeatError::NoSource));
    assert!(matches!(
        solver.distance(&[HeatSource::Vertex(99)]),
        Err(HeatError::SourceOutOfRange { vertex: 99, .. })
    ));
    assert!(matches!(
        solver.value_at(&[0.0; 3], SurfacePoint::Vertex(0)),
        Err(HeatError::ValueCount { got: 3, .. })
    ));
    for options in [
        HeatOptions {
            time_factor: 0.0,
            ..HeatOptions::default()
        },
        HeatOptions {
            time: Some(f64::NAN),
            ..HeatOptions::default()
        },
    ] {
        assert!(matches!(
            HeatSolver::new(&tri, &options),
            Err(HeatError::InvalidOptions(_))
        ));
    }
}
