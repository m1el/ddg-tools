//! The cotan Laplacian, the lumped vertex areas and the quality numbers of an intrinsic
//! triangulation, with the conventions of `ddg_core::operators` (the repository's `CLAUDE.md`).

mod common;

use common::{Shape, long_torus, needle_tetrahedron, sheared_parallelogram, standard_shapes, tube};
use ddg_core::SparseMatrix;
use ddg_core::generators::{self, SplitMix64};
use ddg_core::operators::{cotan_laplacian, lumped_mass, lumped_mass_normalized};
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation};

fn energy(l: &SparseMatrix, x: &[f64]) -> f64 {
    let mut y = vec![0.0; x.len()];
    l.mul_vec(x, &mut y).unwrap();
    x.iter().zip(&y).map(|(a, b)| a * b).sum()
}

fn random_vector(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = SplitMix64::new(seed);
    (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect()
}

/// The Laplacian assembled edge by edge from the cotan weights, as a dense matrix: an
/// independent path to the same entries.
fn dense_from_edges(tri: &IntrinsicTriangulation) -> Vec<Vec<f64>> {
    let n = tri.n_vertices();
    let mut dense = vec![vec![0.0; n]; n];
    for e in 0..tri.n_edges() {
        let [a, b] = tri.edge_vertices(e);
        let w = tri.edge_cotan_weight(e);
        // (u_a − u_b)² w: nothing for a self-edge.
        dense[a][a] += w;
        dense[b][b] += w;
        dense[a][b] -= w;
        dense[b][a] -= w;
    }
    dense
}

fn assert_laplacian(shape: &Shape, tri: &IntrinsicTriangulation, what: &str) {
    let l = tri.cotan_laplacian();
    let n = tri.n_vertices();
    let context = format!("{} ({what})", shape.name);
    assert_eq!((l.nrows(), l.ncols()), (n, n), "{context}");
    assert_eq!(
        l.symmetry_defect().unwrap(),
        0.0,
        "{context}: not exactly symmetric"
    );
    let scale = l.values().iter().fold(0.0f64, |a, &b| a.max(b.abs()));
    assert!(scale.is_finite() && scale > 0.0, "{context}");
    // Zero row sums: constants are in the kernel.
    let ones = vec![1.0; n];
    let mut rows = vec![0.0; n];
    l.mul_vec(&ones, &mut rows).unwrap();
    for (v, r) in rows.iter().enumerate() {
        assert!(r.abs() <= 1e-13 * scale, "{context}: row {v} sums to {r}");
    }
    // Positive semidefinite.
    for seed in 0..5 {
        let x = random_vector(n, seed);
        let e = energy(&l, &x);
        assert!(e >= -1e-12 * scale, "{context}: energy {e}");
    }
    // Entry by entry against the edge weights; the pattern is adjacency plus diagonal.
    let dense = dense_from_edges(tri);
    let mut stored = vec![vec![false; n]; n];
    for (i, j, value) in l.triplets() {
        let (i, j) = (i as usize, j as usize);
        stored[i][j] = true;
        assert!(
            (value - dense[i][j]).abs() <= 1e-12 * scale,
            "{context}: L[{i}, {j}] = {value}, from the edges {}",
            dense[i][j]
        );
    }
    for e in 0..tri.n_edges() {
        let [a, b] = tri.edge_vertices(e);
        assert!(
            stored[a][b] && stored[b][a] && stored[a][a],
            "{context}: edge {e}"
        );
    }
    for i in 0..n {
        for j in 0..n {
            assert!(
                stored[i][j] || dense[i][j] == 0.0,
                "{context}: ({i}, {j}) is missing"
            );
        }
    }
}

#[test]
fn a_fresh_triangulation_has_the_operators_of_the_input() {
    for shape in standard_shapes() {
        let tri = shape.triangulation();
        // Bit for bit: the same formulas on the same lengths, summed in the same order.
        assert_eq!(
            tri.cotan_laplacian(),
            cotan_laplacian(&shape.mesh, &shape.geometry).unwrap(),
            "{}",
            shape.name
        );
        assert_eq!(
            tri.lumped_mass(),
            lumped_mass(&shape.mesh, &shape.geometry).unwrap()
        );
        let normalized = lumped_mass_normalized(&shape.mesh, &shape.geometry).unwrap();
        for (ours, theirs) in tri.lumped_mass_normalized().iter().zip(&normalized) {
            assert!((ours - theirs).abs() < 1e-15);
        }
        assert!((tri.total_area() - shape.geometry.total_area()).abs() < 1e-14);
        assert_laplacian(&shape, &tri, "fresh");
    }
}

#[test]
fn the_laplacian_of_a_flipped_triangulation_is_a_cotan_laplacian() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let mut tri = shape.triangulation();
        let mut rng = SplitMix64::new(77 + s as u64);
        for _ in 0..2 * tri.n_edges() {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let _ = tri.flip_edge(e);
        }
        assert_laplacian(shape, &tri, "random flips");
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_laplacian(shape, &tri, "Delaunay");
        // Where the result is a simplicial complex, ddg_core computes the same operators on it.
        if let Ok((mesh, geometry)) = tri.to_tri_mesh() {
            let theirs = cotan_laplacian(&mesh, &geometry).unwrap();
            let ours = tri.cotan_laplacian();
            assert_eq!(ours.col_ptr(), theirs.col_ptr(), "{}", shape.name);
            assert_eq!(ours.row_idx(), theirs.row_idx(), "{}", shape.name);
            let scale = ours.values().iter().fold(0.0f64, |a, &b| a.max(b.abs()));
            for (a, b) in ours.values().iter().zip(theirs.values()) {
                assert!((a - b).abs() <= 1e-13 * scale, "{}: {a} {b}", shape.name);
            }
            assert_eq!(tri.lumped_mass(), lumped_mass(&mesh, &geometry).unwrap());
        }
    }
}

/// Double edges add up in one entry; a self-edge does not appear at all.
#[test]
fn double_edges_and_self_edges_in_the_laplacian() {
    for shape in [
        long_torus(8, 3, 4.0, 0.2, 0.37),
        needle_tetrahedron(8.0, 4.0),
        tube(3, 6, 0.02, 10.0, 0.37, 1.0),
    ] {
        let mut tri = shape.triangulation();
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let irregular = tri.irregularities();
        assert!(!irregular.self_edges.is_empty() && !irregular.multi_edges.is_empty());
        assert_laplacian(&shape, &tri, "Delaunay");
        let l = tri.cotan_laplacian();
        // A pair of vertices joined by two edges: the entry is minus the sum of both weights.
        let e = irregular.multi_edges[0] as usize;
        let [a, b] = tri.edge_vertices(e);
        let between: Vec<f64> = (0..tri.n_edges())
            .filter(|&o| {
                let [p, q] = tri.edge_vertices(o);
                (p, q) == (a, b) || (p, q) == (b, a)
            })
            .map(|o| tri.edge_cotan_weight(o))
            .collect();
        assert!(between.len() >= 2, "{}", shape.name);
        assert!((l.get(a, b) + between.iter().sum::<f64>()).abs() < 1e-12);
        // A vertex with a self-edge: its diagonal entry is the sum over its other edges.
        let loop_edge = irregular.self_edges[0] as usize;
        let [v, _] = tri.edge_vertices(loop_edge);
        let others: f64 = (0..tri.n_edges())
            .filter(|&o| {
                let [p, q] = tri.edge_vertices(o);
                (p == v) != (q == v)
            })
            .map(|o| tri.edge_cotan_weight(o))
            .sum();
        assert!(
            (l.get(v, v) - others).abs() < 1e-12 * others.abs().max(1.0),
            "{}",
            shape.name
        );
        // The lumped areas: a face counts once for each of its corners at a vertex.
        let mass = tri.lumped_mass();
        assert!((mass.iter().sum::<f64>() - tri.total_area()).abs() < 1e-13 * tri.total_area());
        let f = irregular.faces_with_repeated_vertex[0] as usize;
        let face = tri.faces()[f];
        let repeated = if face[0] == face[1] || face[0] == face[2] {
            face[0]
        } else {
            face[1]
        };
        let corners = face.iter().filter(|&&x| x == repeated).count();
        assert!(corners >= 2);
        let expected: f64 = (0..tri.n_faces())
            .map(|g| {
                let n = tri.faces()[g].iter().filter(|&&x| x == repeated).count();
                n as f64 * tri.face_area(g) / 3.0
            })
            .sum();
        assert!((mass[repeated as usize] - expected).abs() < 1e-14 * expected.max(1.0));
        let normalized = tri.lumped_mass_normalized();
        assert!((normalized.iter().sum::<f64>() - 1.0).abs() < 1e-14);
    }
}

/// Rippa's theorem, edge by edge: a Delaunay flip does not increase the Dirichlet energy of
/// the piecewise linear interpolant of any values at the vertices.
#[test]
fn delaunay_flipping_does_not_increase_the_dirichlet_energy() {
    for (s, shape) in standard_shapes().iter().enumerate() {
        let before = shape.triangulation();
        let mut after = before.clone();
        let report = after.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let (l0, l1) = (before.cotan_laplacian(), after.cotan_laplacian());
        let mut strictly = 0;
        for seed in 0..8 {
            let x = random_vector(before.n_vertices(), 500 + 10 * s as u64 + seed);
            let (e0, e1) = (energy(&l0, &x), energy(&l1, &x));
            assert!(e1 <= e0 * (1.0 + 1e-12), "{}: {e0} -> {e1}", shape.name);
            strictly += usize::from(e1 < e0 * (1.0 - 1e-9));
        }
        if report.n_flips > 0 {
            assert!(
                strictly > 0,
                "{}: {} flips left every energy as it was",
                shape.name,
                report.n_flips
            );
        }
        // Linear functions on a flat shape have the same energy on every triangulation: the
        // integral of |∇f|² over the shape.
        if shape.flat {
            let x: Vec<f64> = shape
                .positions
                .iter()
                .map(|p| 2.0 * p[0] - 3.0 * p[1])
                .collect();
            let exact = 13.0 * shape.geometry.total_area();
            assert!(
                (energy(&l0, &x) - exact).abs() < 1e-11 * exact,
                "{}",
                shape.name
            );
            assert!(
                (energy(&l1, &x) - exact).abs() < 1e-11 * exact,
                "{}",
                shape.name
            );
        }
    }
}

#[test]
fn quality_numbers_of_known_meshes() {
    let deg = f64::to_degrees;
    // Right isosceles triangles.
    let grid = Shape::generated("grid", generators::grid_disk(4, 3).unwrap(), true);
    let q = grid.triangulation().quality(1e-12);
    // 43 edges, 14 of them on the boundary.
    assert_eq!((q.n_faces, q.n_interior_edges), (24, 29));
    assert!((deg(q.min_angle) - 45.0).abs() < 1e-12 && (deg(q.max_angle) - 90.0).abs() < 1e-12);
    for a in q.largest_angle_quantiles {
        assert!((deg(a) - 90.0).abs() < 1e-12);
    }
    assert_eq!(q.largest_angle_above, [0.0, 0.0, 0.0]);
    assert_eq!(
        (q.n_negative_weights, q.negative_weight_fraction()),
        (0, 0.0)
    );
    assert!(q.min_cotan_weight.abs() < 1e-15, "the diagonals are ties");
    // With a threshold of zero the ties whose weight is negative by rounding are counted.
    assert!(grid.triangulation().quality(0.0).n_negative_weights <= 12);

    // Sheared by one cell per row: angles 45°, 18.43° (atan 1/3)... every face has an angle of
    // 135° or more.
    let sheared = sheared_parallelogram(6, 1.0);
    let mut tri = sheared.triangulation();
    let q = tri.quality(1e-12);
    assert!(
        (deg(q.largest_angle_quantiles[0]) - 135.0).abs() < 1e-10,
        "{q:?}"
    );
    assert_eq!(q.largest_angle_above, [1.0, 1.0, 0.0]);
    assert!(
        q.negative_weight_fraction() > 0.3 && q.min_cotan_weight < -0.4,
        "{q:?}"
    );
    // Its Delaunay triangulation is the one of the square lattice, right isosceles triangles,
    // except along the two slanted sides.
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    let q = tri.quality(1e-12);
    assert_eq!(q.n_negative_weights, 0, "{q:?}");
    assert!(
        (deg(q.largest_angle_quantiles[0]) - 90.0).abs() < 1e-10,
        "{q:?}"
    );
    assert!(q.min_cotan_weight > -1e-12);
    assert_eq!(
        q.n_interior_edges,
        tri.n_edges() - sheared.mesh.n_boundary_edges()
    );
}
