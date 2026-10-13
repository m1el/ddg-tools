//! SPEC §5 item 4: after Delaunay flipping every interior edge satisfies the Delaunay condition
//! and has a non-negative cotan weight, and the result does not depend on the order in which
//! edges are flipped, except at cocircular ties.

mod common;

use common::{
    Shape, edge_list, edge_list_difference, long_torus, needle_tetrahedron, sheared_parallelogram,
    standard_shapes, thin_cylinder, tube,
};
use ddg_core::generators::{self, SplitMix64};
use ddg_intrinsic::{
    DelaunayCriterion, DelaunayOptions, DelaunayStatus, FlipLimit, IntrinsicError,
    IntrinsicTriangulation,
};

const TOLERANCE: f64 = 1e-10;

fn delaunay(shape: &Shape) -> (IntrinsicTriangulation, ddg_intrinsic::DelaunayReport) {
    let mut tri = shape.triangulation();
    let report = tri
        .flip_to_delaunay(&DelaunayOptions::default())
        .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
    (tri, report)
}

#[test]
fn every_interior_edge_is_delaunay_afterwards() {
    for shape in standard_shapes() {
        let (tri, report) = delaunay(&shape);
        let context = format!("{}: {report:?}", shape.name);
        assert_eq!(report.status, DelaunayStatus::Delaunay, "{context}");
        assert_eq!(report.n_not_delaunay, 0, "{context}");
        assert_eq!(report.n_refused, 0, "{context}");
        assert!(report.n_tested >= tri.n_edges(), "{context}");
        assert!(report.max_angle_excess <= TOLERANCE, "{context}");
        assert!(tri.is_delaunay(&DelaunayCriterion::default()), "{context}");
        let mut largest = f64::NEG_INFINITY;
        let mut smallest = f64::INFINITY;
        let mut ties = 0;
        for e in 0..tri.n_edges() {
            let weight = tri.edge_cotan_weight(e);
            assert!(weight.is_finite(), "{context}");
            let Some(excess) = tri.edge_angle_excess(e) else {
                assert!(tri.is_boundary_edge(e));
                continue;
            };
            assert!(
                excess <= TOLERANCE,
                "{context}: edge {e} has excess {excess}"
            );
            // The weight has the sign of minus the excess; where the excess is within the
            // tolerance of zero the weight is zero up to the tolerance over the sines.
            assert!(
                weight >= -1e-8,
                "{context}: edge {e} has cotan weight {weight} (excess {excess})"
            );
            if excess < -TOLERANCE {
                assert!(weight > 0.0, "{context}: edge {e}: {weight} {excess}");
            }
            largest = largest.max(excess);
            smallest = smallest.min(weight);
            ties += usize::from(excess.abs() <= TOLERANCE);
        }
        assert_eq!(report.max_angle_excess, largest, "{context}");
        assert_eq!(report.min_cotan_weight, smallest, "{context}");
        assert_eq!(report.n_ties, ties, "{context}");
        // No positive off-diagonal entry in the Laplacian between vertices that are not both on
        // the boundary: the maximum principle. A boundary edge cannot be flipped and keeps the
        // weight ½ cot of its one opposite angle, which is negative if that angle is obtuse.
        let laplacian = tri.cotan_laplacian();
        for (i, j, value) in laplacian.triplets() {
            let both_on_boundary =
                tri.is_boundary_vertex(i as usize) && tri.is_boundary_vertex(j as usize);
            assert!(
                i == j || both_on_boundary || value <= 1e-8,
                "{context}: L[{i}, {j}] = {value}"
            );
        }
        // A second run has nothing to do.
        let mut again = tri.clone();
        let second = again.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_eq!(second.n_flips, 0, "{context}");
        assert_eq!(second.n_tested, tri.n_edges(), "{context}");
        assert_eq!(again, tri, "{context}");
    }
}

/// Flips, in a random order, edges that fail the criterion until none is left.
fn delaunay_in_random_order(shape: &Shape, seed: u64) -> (IntrinsicTriangulation, usize) {
    let mut tri = shape.triangulation();
    let mut rng = SplitMix64::new(seed);
    let mut flips = 0;
    loop {
        let failing: Vec<usize> = (0..tri.n_edges())
            .filter(|&e| tri.edge_angle_excess(e).is_some_and(|x| x > TOLERANCE))
            .collect();
        if failing.is_empty() {
            return (tri, flips);
        }
        // Flip a random half of them (those that still fail when their turn comes).
        for _ in 0..failing.len().div_ceil(2) {
            let e = failing[(rng.next_u64() % failing.len() as u64) as usize];
            if tri.edge_angle_excess(e).is_some_and(|x| x > TOLERANCE) {
                tri.flip_edge(e)
                    .unwrap_or_else(|err| panic!("{}: {err}", shape.name));
                flips += 1;
            }
        }
        assert!(
            flips < 100 * tri.n_edges(),
            "{}: no end in sight",
            shape.name
        );
    }
}

#[test]
fn the_result_does_not_depend_on_the_order_of_flips() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let (queue, report) = delaunay(shape);
        let reference = edge_list(&queue);
        for seed in 0..3 {
            let (random, flips) = delaunay_in_random_order(shape, 1000 * s as u64 + seed);
            let other = edge_list(&random);
            let unmatched = edge_list_difference(&reference, &other, 1e-10);
            let context = format!(
                "{}: {} flips by the queue, {flips} in random order {seed}, {} ties, {unmatched} \
                 edges differ",
                shape.name, report.n_flips, report.n_ties
            );
            // Without ties the triangulation is unique. A tie is a quadrilateral with both
            // diagonals Delaunay; two orders can differ in the diagonal of each tie, which is
            // one edge of each list.
            assert!(unmatched <= 2 * report.n_ties, "{context}");
            if report.n_ties == 0 {
                assert_eq!(unmatched, 0, "{context}");
            }
            // The triangulations may differ; the Dirichlet energy of any function on the
            // vertices does not (both are Delaunay, ties do not change it).
            let x: Vec<f64> = (0..queue.n_vertices())
                .map(|v| ((v * 37 % 101) as f64).sin())
                .collect();
            let energy = |tri: &IntrinsicTriangulation| {
                let l = tri.cotan_laplacian();
                let mut y = vec![0.0; x.len()];
                l.mul_vec(&x, &mut y).unwrap();
                x.iter().zip(&y).map(|(a, b)| a * b).sum::<f64>()
            };
            let (a, b) = (energy(&queue), energy(&random));
            assert!(
                (a - b).abs() <= 1e-9 * a.abs().max(1.0),
                "{context}: energies {a} and {b}"
            );
        }
    }
}

/// Where the ties are: on the generators that are regular by construction.
#[test]
fn cocircular_ties_are_reported_and_left_alone() {
    // A grid of squares: every interior diagonal is a tie, nothing is flipped.
    let grid = Shape::generated("grid", generators::grid_disk(5, 4).unwrap(), true);
    let (tri, report) = delaunay(&grid);
    assert_eq!(report.n_flips, 0);
    assert_eq!(report.n_ties, 20, "{report:?}");
    assert!(report.max_angle_excess.abs() < 1e-14);
    assert_eq!(tri, grid.triangulation());
    // The sheared lattice has the same Delaunay triangulation as the square lattice it is a
    // shear of only for the interior; here all that is checked is that ties are counted and
    // that the result is Delaunay.
    let sheared = sheared_parallelogram(8, 1.0);
    let (_, report) = delaunay(&sheared);
    assert_eq!(report.status, DelaunayStatus::Delaunay);
    assert!(report.n_flips > 0);
    // The icosphere is Delaunay as it is; its ties come from the symmetric subdivision.
    let sphere = Shape::generated("icosphere", generators::icosphere(2).unwrap(), false);
    let (_, report) = delaunay(&sphere);
    assert_eq!(report.n_flips, 0, "{report:?}");
}

/// On a convex planar domain the intrinsic Delaunay triangulation is the Delaunay
/// triangulation of the points: no vertex lies inside the circumcircle of a triangle.
#[test]
fn on_a_convex_flat_domain_circumcircles_are_empty() {
    for shape in [
        sheared_parallelogram(8, 1.0),
        sheared_parallelogram(7, 2.5),
        common::jittered_grid(9, 0.3, 12),
    ] {
        let (tri, report) = delaunay(&shape);
        assert!(report.n_flips > 0, "{}", shape.name);
        let p = &shape.positions;
        common::assert_not_folded(&shape);
        for (f, face) in tri.faces().iter().enumerate() {
            let [a, b, c] = face.map(|v| [p[v as usize][0], p[v as usize][1]]);
            // Circumcentre and squared radius.
            let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
            assert!(
                d > 0.0,
                "{}: face {f} is not counter-clockwise in the plane",
                shape.name
            );
            let (a2, b2, c2) = (
                a[0] * a[0] + a[1] * a[1],
                b[0] * b[0] + b[1] * b[1],
                c[0] * c[0] + c[1] * c[1],
            );
            let centre = [
                (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d,
                (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d,
            ];
            let radius = (a[0] - centre[0]).hypot(a[1] - centre[1]);
            for (v, q) in p.iter().enumerate() {
                let distance = (q[0] - centre[0]).hypot(q[1] - centre[1]);
                assert!(
                    distance >= radius * (1.0 - 1e-9),
                    "{}: vertex {v} is inside the circumcircle of face {f} = {face:?}",
                    shape.name
                );
            }
        }
    }
}

#[test]
fn the_flip_limit_is_reported() {
    let shape = sheared_parallelogram(8, 3.0);
    let (_, full) = delaunay(&shape);
    assert!(full.n_flips > 20);
    let mut tri = shape.triangulation();
    let limited = DelaunayOptions {
        max_flips: FlipLimit::Flips(7),
        ..DelaunayOptions::default()
    };
    let report = tri.flip_to_delaunay(&limited).unwrap();
    assert_eq!(report.status, DelaunayStatus::FlipLimit);
    assert_eq!(report.n_flips, 7);
    assert!(report.n_not_delaunay > 0 && report.max_angle_excess > TOLERANCE);
    assert!(!tri.is_delaunay(&DelaunayCriterion::default()));
    tri.check().unwrap();
    // A limit of zero flips changes nothing.
    let mut untouched = shape.triangulation();
    let none = DelaunayOptions {
        max_flips: FlipLimit::Flips(0),
        ..DelaunayOptions::default()
    };
    let report = untouched.flip_to_delaunay(&none).unwrap();
    assert_eq!(
        (report.status, report.n_flips),
        (DelaunayStatus::FlipLimit, 0)
    );
    assert_eq!(untouched, shape.triangulation());
    // Carrying on reaches the same triangulation as an unlimited run.
    let rest = tri
        .flip_to_delaunay(&DelaunayOptions {
            max_flips: FlipLimit::Unlimited,
            ..DelaunayOptions::default()
        })
        .unwrap();
    assert_eq!(rest.status, DelaunayStatus::Delaunay);
    let (unlimited, _) = delaunay(&shape);
    let difference = edge_list_difference(&edge_list(&tri), &edge_list(&unlimited), 1e-10);
    assert!(difference <= 2 * rest.n_ties, "{difference} {rest:?}");
    // A limit per edge.
    let mut tri = shape.triangulation();
    let per_edge = DelaunayOptions {
        max_flips: FlipLimit::PerEdge(1),
        ..DelaunayOptions::default()
    };
    assert_eq!(
        tri.flip_to_delaunay(&per_edge).unwrap().status,
        DelaunayStatus::Delaunay
    );
}

#[test]
fn the_reference_criterion_bounds_the_cotan_weights() {
    for shape in standard_shapes() {
        let mut tri = shape.triangulation();
        let options = DelaunayOptions {
            criterion: DelaunayCriterion::REFERENCE,
            ..DelaunayOptions::default()
        };
        let report = tri.flip_to_delaunay(&options).unwrap();
        assert_eq!(
            report.status,
            DelaunayStatus::Delaunay,
            "{}: {report:?}",
            shape.name
        );
        assert!(
            report.min_cotan_weight >= -1e-6,
            "{}: {report:?}",
            shape.name
        );
        assert!(tri.is_delaunay(&DelaunayCriterion::REFERENCE));
        // It is the same triangulation as with the angle criterion, up to edges whose weight
        // is within 1e-6 of zero.
        let (ours, _) = delaunay(&shape);
        let near_ties = (0..ours.n_edges())
            .filter(|&e| !ours.is_boundary_edge(e) && ours.edge_cotan_weight(e).abs() <= 1e-6)
            .count();
        let difference = edge_list_difference(&edge_list(&tri), &edge_list(&ours), 1e-10);
        assert!(
            difference <= 2 * near_ties,
            "{}: {difference} > 2 * {near_ties}",
            shape.name
        );
    }
}

#[test]
fn an_invalid_tolerance_is_an_error() {
    let shape = sheared_parallelogram(3, 1.0);
    for tolerance in [-1e-12, f64::NAN, f64::INFINITY] {
        for criterion in [
            DelaunayCriterion::AngleSum { tolerance },
            DelaunayCriterion::CotanWeight { tolerance },
        ] {
            let mut tri = shape.triangulation();
            let error = tri
                .flip_to_delaunay(&DelaunayOptions {
                    criterion,
                    ..DelaunayOptions::default()
                })
                .unwrap_err();
            assert!(
                matches!(error, IntrinsicError::InvalidOptions { .. }),
                "{error}"
            );
            assert!(error.to_string().contains("tolerance"), "{error}");
            assert_eq!(tri, shape.triangulation());
        }
    }
}

/// Meshes on which Delaunay flipping leaves the class of simplicial complexes.
#[test]
fn delaunay_flipping_produces_self_edges_and_double_edges() {
    // The coarse torus with a thin tube of the reference runs: 6 self-edges, 16 vertex pairs
    // joined by two edges (data/intrinsic-reference/SUMMARY.md).
    let torus = long_torus(8, 3, 4.0, 0.2, 0.37);
    let (tri, report) = delaunay(&torus);
    let irregular = tri.irregularities();
    println!("{}: {report:?}\n  {irregular:?}", torus.name);
    assert_eq!(irregular.self_edges.len(), 6);
    let mut pairs: Vec<[usize; 2]> = irregular
        .multi_edges
        .iter()
        .map(|&e| {
            let [a, b] = tri.edge_vertices(e as usize);
            [a.min(b), a.max(b)]
        })
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    assert_eq!(pairs.len(), 16, "vertex pairs joined by more than one edge");
    assert!(irregular.multi_edges.len() >= 32);
    assert!(!irregular.faces_with_repeated_vertex.is_empty());
    for &e in &irregular.self_edges {
        let [a, b] = tri.edge_vertices(e as usize);
        assert_eq!(a, b);
        // A loop around the tube: about as long as the tube's circumference.
        let length = tri.edge_length(e as usize);
        assert!(
            length > 0.9 && length < 1.3,
            "self-edge {e} has length {length}"
        );
    }
    assert!(matches!(
        tri.to_tri_mesh().unwrap_err(),
        IntrinsicError::NotSimplicial { .. }
    ));
    tri.check().unwrap();

    // The tilted needle: a vertex of degree one.
    let needle = needle_tetrahedron(8.0, 4.0);
    let (tri, report) = delaunay(&needle);
    let irregular = tri.irregularities();
    println!("{}: {report:?}\n  {irregular:?}", needle.name);
    assert_eq!(irregular.degree_one_vertices.len(), 1);
    assert_eq!(irregular.self_edges.len(), 1);
    assert_eq!(irregular.multi_edges.len(), 2);

    // A thin tube wound with few vertices: double edges, no self-edge. The tube of the
    // reference runs stays simplicial; a thinner one does not.
    let reference_tube = thin_cylinder(4, 12, 0.1, 10.0, 0.37);
    let (tri, report) = delaunay(&reference_tube);
    println!(
        "{}: {report:?}\n  {:?}",
        reference_tube.name,
        tri.irregularities()
    );
    assert!(tri.irregularities().is_empty());
    assert_eq!(report.n_flips, 48, "as the reference (stats.txt)");
    // With the vertices on a helix a thin tube gets loops around it, on a surface with
    // boundary.
    let helix = tube(3, 6, 0.02, 10.0, 0.37, 1.0);
    let (tri, report) = delaunay(&helix);
    let irregular = tri.irregularities();
    println!("{}: {report:?}\n  {irregular:?}", helix.name);
    assert!(irregular.self_edges.len() >= 10 && irregular.multi_edges.len() >= 10);
    for &e in &irregular.self_edges {
        // A loop around the tube has the length of its circumference, here the perimeter of
        // the triangle the three vertices of a ring would span.
        let length = tri.edge_length(e as usize);
        let around = 0.02 * 2.0 * std::f64::consts::PI;
        assert!(
            length > 0.8 * around && length < 1.05 * around,
            "self-edge {e}: {length}"
        );
    }
    tri.check().unwrap();
}
