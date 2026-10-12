//! Property tests of the mesh, geometry, operator, sparse and I/O layers on generated meshes.

// Index loops mirror the subscripts of the formulas being checked.
#![allow(clippy::needless_range_loop)]

use std::f64::consts::PI;

use ddg_core::generators::{self, GeneratedMesh, SplitMix64};
use ddg_core::geometry::{Geometry, corner_angle, corner_cotan, triangle_area};
use ddg_core::io::{obj_to_string, parse_obj};
use ddg_core::operators::{
    angle_defects, boundary_geodesic_curvatures, cotan_laplacian, gauss_bonnet_residual,
    lumped_mass, lumped_mass_normalized,
};
use ddg_core::sparse::{
    Method, Parallelism, SparseMatrix, Symbolic, ThreadPool, relative_residual,
};
use ddg_core::topology::{MeshError, TriMesh, compact_vertices, orient_faces};
use proptest::prelude::*;

/// What a generated mesh is expected to be topologically.
#[derive(Debug, Clone, Copy)]
struct Expect {
    euler: i64,
    genus: usize,
    boundary_loops: usize,
    /// Whether all vertices lie in the plane z = 0.
    flat: bool,
}

/// A mesh of one of the four topology classes, with positions perturbed by less than a
/// quarter of the shortest edge so that every triangle stays non-degenerate.
fn arbitrary_mesh() -> impl Strategy<Value = (GeneratedMesh, Expect)> {
    let sphere = (0u32..3, any::<u64>(), 0.0..0.2f64).prop_map(|(level, seed, amp)| {
        let mut m = generators::icosphere(level).unwrap();
        perturb_relative(&mut m, [amp, amp, amp], seed);
        (
            m,
            Expect {
                euler: 2,
                genus: 0,
                boundary_loops: 0,
                flat: false,
            },
        )
    });
    let torus =
        (3usize..14, 3usize..10, any::<u64>(), 0.0..0.2f64).prop_map(|(a, b, seed, amp)| {
            let mut m = generators::torus(a, b, 2.0, 0.7).unwrap();
            perturb_relative(&mut m, [amp, amp, amp], seed);
            (
                m,
                Expect {
                    euler: 0,
                    genus: 1,
                    boundary_loops: 0,
                    flat: false,
                },
            )
        });
    let disk = (
        1usize..9,
        1usize..9,
        any::<u64>(),
        0.0..0.2f64,
        prop::bool::ANY,
        0.0..0.3f64,
    )
        .prop_map(|(nx, ny, seed, amp, flat, bump)| {
            let mut m = if flat {
                generators::grid_disk(nx, ny).unwrap()
            } else {
                generators::bumpy_grid_disk(nx, ny, bump).unwrap()
            };
            perturb_relative(&mut m, [amp, amp, if flat { 0.0 } else { amp }], seed);
            (
                m,
                Expect {
                    euler: 1,
                    genus: 0,
                    boundary_loops: 1,
                    flat,
                },
            )
        });
    let annulus = (
        3usize..14,
        1usize..5,
        any::<u64>(),
        0.0..0.2f64,
        prop::bool::ANY,
    )
        .prop_map(|(around, across, seed, amp, flat)| {
            let mut m = if flat {
                generators::annulus(around, across, 0.6, 1.7).unwrap()
            } else {
                generators::open_cylinder(around, across, 1.0, 1.5).unwrap()
            };
            perturb_relative(&mut m, [amp, amp, if flat { 0.0 } else { amp }], seed);
            (
                m,
                Expect {
                    euler: 0,
                    genus: 0,
                    boundary_loops: 2,
                    flat,
                },
            )
        });
    prop_oneof![sphere, torus, disk, annulus]
}

/// Perturbs positions by `fraction` of the shortest edge length (per axis).
fn perturb_relative(m: &mut GeneratedMesh, fraction: [f64; 3], seed: u64) {
    let mut shortest = f64::INFINITY;
    for f in &m.faces {
        for k in 0..3 {
            let (p, q) = (
                m.positions[f[k] as usize],
                m.positions[f[(k + 1) % 3] as usize],
            );
            shortest = shortest.min(ddg_core::geometry::edge_length(p, q));
        }
    }
    generators::perturb(&mut m.positions, fraction.map(|a| a * shortest), seed);
}

fn random_vector(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = SplitMix64::new(seed);
    (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    #[test]
    fn topology_invariants((shape, expect) in arbitrary_mesh()) {
        let mesh = shape.topology().unwrap();
        prop_assert_eq!(mesh.euler_characteristic(), expect.euler);
        prop_assert_eq!(mesh.genus(), expect.genus);
        prop_assert_eq!(mesh.boundary_loops().len(), expect.boundary_loops);
        prop_assert_eq!(mesh.n_components(), 1);
        prop_assert_eq!(mesh.n_halfedges(), 3 * mesh.n_faces());

        // Twins are an involution that reverses direction and preserves the edge.
        let mut boundary_halfedges = 0;
        for h in 0..mesh.n_halfedges() {
            prop_assert_eq!(mesh.from_vertex(mesh.next(h)), mesh.to_vertex(h));
            prop_assert_eq!(mesh.next(mesh.prev(h)), h);
            let e = mesh.edge(h);
            let [a, b] = mesh.edge_vertices(e);
            prop_assert!(a < b);
            let (from, to) = (mesh.from_vertex(h), mesh.to_vertex(h));
            prop_assert_eq!((from.min(to), from.max(to)), (a, b));
            prop_assert_eq!(mesh.find_edge(from, to), Some(e));
            prop_assert_eq!(mesh.find_halfedge(from, to), Some(h));
            match mesh.twin(h) {
                Some(t) => {
                    prop_assert_eq!(mesh.twin(t), Some(h));
                    prop_assert_eq!((mesh.from_vertex(t), mesh.to_vertex(t)), (to, from));
                    prop_assert_eq!(mesh.edge(t), e);
                    prop_assert!(!mesh.is_boundary_edge(e));
                }
                None => {
                    boundary_halfedges += 1;
                    prop_assert!(mesh.is_boundary_edge(e));
                    prop_assert!(mesh.is_boundary_vertex(from) && mesh.is_boundary_vertex(to));
                }
            }
            prop_assert!(mesh.edge_halfedge(e) == h || mesh.twin(mesh.edge_halfedge(e)) == Some(h));
        }
        prop_assert_eq!(boundary_halfedges, mesh.n_boundary_edges());
        prop_assert_eq!(2 * mesh.n_edges(), mesh.n_halfedges() + boundary_halfedges);
        // Edges are sorted and unique.
        prop_assert!(mesh.edges().windows(2).all(|w| w[0] < w[1]));

        // One-rings: every outgoing halfedge starts at the vertex; degrees add up.
        let mut total_faces = 0;
        for v in 0..mesh.n_vertices() {
            let out: Vec<usize> = mesh.outgoing_halfedges(v).collect();
            prop_assert_eq!(out[0], mesh.outgoing_halfedge(v));
            prop_assert!(out.iter().all(|&h| mesh.from_vertex(h) == v));
            total_faces += out.len();
            prop_assert_eq!(out.len(), mesh.vertex_face_count(v));
            let neighbors: Vec<usize> = mesh.vertex_neighbors(v).collect();
            let boundary = mesh.is_boundary_vertex(v);
            prop_assert_eq!(neighbors.len(), out.len() + usize::from(boundary));
            let mut sorted = neighbors.clone();
            sorted.sort_unstable();
            sorted.dedup();
            prop_assert_eq!(sorted.len(), neighbors.len());
            prop_assert!(neighbors.iter().all(|&w| mesh.find_edge(v, w).is_some()));
            // Consecutive neighbours span a face with v, counter-clockwise.
            for (i, &h) in out.iter().enumerate() {
                prop_assert_eq!(mesh.to_vertex(h), neighbors[i]);
                prop_assert_eq!(mesh.opposite_vertex(h), neighbors[(i + 1) % neighbors.len()]);
            }
        }
        prop_assert_eq!(total_faces, 3 * mesh.n_faces());

        // Boundary loops are closed chains of boundary halfedges and cover the boundary once.
        let mut on_loops = 0;
        for l in mesh.boundary_loops() {
            for (i, &v) in l.iter().enumerate() {
                let h = mesh.outgoing_halfedge(v as usize);
                prop_assert_eq!(mesh.twin(h), None);
                prop_assert_eq!(mesh.to_vertex(h), l[(i + 1) % l.len()] as usize);
            }
            on_loops += l.len();
        }
        prop_assert_eq!(on_loops, mesh.n_boundary_edges());
        prop_assert_eq!(mesh.boundary_vertices().len(), mesh.n_boundary_edges());
        prop_assert_eq!(mesh.boundary_vertices().len() + mesh.interior_vertices().len(), mesh.n_vertices());
    }

    #[test]
    fn gauss_bonnet_holds((shape, expect) in arbitrary_mesh()) {
        let mesh = shape.topology().unwrap();
        let g = Geometry::from_positions(&mesh, &shape.positions).unwrap();
        let residual = gauss_bonnet_residual(&mesh, &g).unwrap();
        let tolerance = 1e-12 * (1.0 + mesh.n_vertices() as f64);
        prop_assert!(residual.abs() < tolerance, "residual {residual}");
        let omega: f64 = angle_defects(&mesh, &g).unwrap().iter().sum();
        let kappa: f64 = boundary_geodesic_curvatures(&mesh, &g).unwrap().iter().sum();
        prop_assert!((omega + kappa - 2.0 * PI * expect.euler as f64).abs() < tolerance);
        if expect.flat {
            // A flat mesh has no interior curvature at all.
            for (v, &w) in angle_defects(&mesh, &g).unwrap().iter().enumerate() {
                prop_assert!(w.abs() < 1e-12, "vertex {v}: defect {w}");
            }
        }
        // Triangle angles sum to pi, and the three scalar formulas agree with each other.
        for f in 0..mesh.n_faces() {
            let a = g.corner_angles()[f];
            prop_assert!((a[0] + a[1] + a[2] - PI).abs() < 1e-13);
            for k in 0..3 {
                prop_assert!(a[k] > 0.0 && a[k] < PI);
                let cot = g.corner_cotans()[f][k];
                prop_assert!((cot * a[k].sin() - a[k].cos()).abs() < 1e-12 * (1.0 + cot.abs()));
            }
            prop_assert!(g.face_areas()[f] > 0.0);
        }
    }

    /// Gauss–Bonnet is a statement about the metric alone: it also holds for edge lengths that
    /// come from no embedding, here a conformal rescaling of a sphere.
    #[test]
    fn gauss_bonnet_holds_intrinsically(level in 0u32..3, seed in any::<u64>(), strength in 0.0..0.3f64) {
        let shape = generators::icosphere(level).unwrap();
        let mesh = shape.topology().unwrap();
        let embedded = Geometry::from_positions(&mesh, &shape.positions).unwrap();
        let u: Vec<f64> = random_vector(mesh.n_vertices(), seed).iter().map(|x| strength * x).collect();
        let lengths: Vec<f64> = mesh
            .edges()
            .iter()
            .zip(embedded.edge_lengths())
            .map(|(&[i, j], &l)| l * ((u[i as usize] + u[j as usize]) / 2.0).exp())
            .collect();
        let g = Geometry::from_edge_lengths(&mesh, lengths).unwrap();
        let residual = gauss_bonnet_residual(&mesh, &g).unwrap();
        prop_assert!(residual.abs() < 1e-12 * (1.0 + mesh.n_vertices() as f64), "residual {residual}");
    }

    #[test]
    fn laplacian_is_symmetric_psd_with_zero_row_sums((shape, expect) in arbitrary_mesh(), seed in any::<u64>()) {
        let mesh = shape.topology().unwrap();
        let g = Geometry::from_positions(&mesh, &shape.positions).unwrap();
        let l = cotan_laplacian(&mesh, &g).unwrap();
        let n = mesh.n_vertices();
        prop_assert_eq!((l.nrows(), l.ncols()), (n, n));
        prop_assert_eq!(l.nnz(), n + 2 * mesh.n_edges());
        prop_assert_eq!(l.symmetry_defect().unwrap(), 0.0);

        // Entries: -1/2 (cot a + cot b) on edges, checked against the corner cotangents.
        for e in 0..mesh.n_edges() {
            let [i, j] = mesh.edge_vertices(e);
            let h = mesh.edge_halfedge(e);
            let mut w = 0.0;
            for h in std::iter::once(h).chain(mesh.twin(h)) {
                let c = mesh.prev(h); // the corner opposite halfedge h
                w += 0.5 * g.corner_cotans()[c / 3][c % 3];
            }
            prop_assert!((l.get(i, j) + w).abs() <= 1e-14 * (1.0 + w.abs()));
        }

        // Zero row sums and scale for tolerances.
        let ones = vec![1.0; n];
        let mut row_sums = vec![0.0; n];
        l.mul_vec(&ones, &mut row_sums).unwrap();
        let scale = l.values().iter().fold(0.0f64, |m, v| m.max(v.abs()));
        for (i, s) in row_sums.iter().enumerate() {
            prop_assert!(s.abs() <= 1e-13 * scale, "row {i} sums to {s} (scale {scale})");
        }

        // Positive semidefinite: x' L x >= 0, and it is the Dirichlet energy, which the
        // constant vector does not have.
        let x = random_vector(n, seed);
        let mut lx = vec![0.0; n];
        l.mul_vec(&x, &mut lx).unwrap();
        let energy: f64 = x.iter().zip(&lx).map(|(a, b)| a * b).sum();
        prop_assert!(energy >= -1e-12 * scale * n as f64, "x'Lx = {energy}");

        // Linear precision: on a flat mesh, L annihilates linear functions at interior vertices.
        if expect.flat {
            let mut rng = SplitMix64::new(seed ^ 0xABCD);
            let (a, b, c) = (rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0));
            let f: Vec<f64> = shape.positions.iter().map(|p| a * p[0] + b * p[1] + c).collect();
            let mut lf = vec![0.0; n];
            l.mul_vec(&f, &mut lf).unwrap();
            for v in mesh.interior_vertices() {
                let v = v as usize;
                prop_assert!(lf[v].abs() <= 1e-11 * scale, "(L f)[{v}] = {} (scale {scale})", lf[v]);
            }
        }
    }

    #[test]
    fn masses_sum_to_the_area((shape, _) in arbitrary_mesh()) {
        let mesh = shape.topology().unwrap();
        let g = Geometry::from_positions(&mesh, &shape.positions).unwrap();
        let m = lumped_mass(&mesh, &g).unwrap();
        let area = g.total_area();
        prop_assert!(m.iter().all(|&x| x > 0.0));
        prop_assert!((m.iter().sum::<f64>() - area).abs() <= 1e-13 * area);
        let mn = lumped_mass_normalized(&mesh, &g).unwrap();
        prop_assert!((mn.iter().sum::<f64>() - 1.0).abs() <= 1e-13);
        // Area from Heron agrees with the cross product.
        let mut cross_area = 0.0;
        for f in &shape.faces {
            let [p, q, r] = f.map(|v| shape.positions[v as usize]);
            let u = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
            let w = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];
            let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
            cross_area += 0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        }
        prop_assert!((cross_area - area).abs() <= 1e-12 * area);
    }

    /// The scalar formulas against an independent evaluation from coordinates.
    #[test]
    fn scalar_formulas_match_coordinates(bx in 0.1..5.0f64, cx in -5.0..5.0f64, cy in 0.01..5.0f64) {
        // Triangle A = (0,0), B = (bx,0), C = (cx,cy).
        let a = ((cx - bx).powi(2) + cy * cy).sqrt(); // |BC|, opposite A
        let b = (cx * cx + cy * cy).sqrt(); // |CA|, opposite B
        let c = bx; // |AB|, opposite C
        let area = triangle_area(a, b, c).unwrap();
        prop_assert!((area - 0.5 * bx * cy).abs() <= 1e-12 * (1.0 + area));
        let angle_a = corner_angle(a, b, c).unwrap();
        prop_assert!((angle_a - cy.atan2(cx)).abs() < 1e-11);
        let angle_b = corner_angle(b, c, a).unwrap();
        prop_assert!((angle_b - cy.atan2(bx - cx)).abs() < 1e-11);
        let angle_c = corner_angle(c, a, b).unwrap();
        prop_assert!((angle_a + angle_b + angle_c - PI).abs() < 1e-11);
        // cot A = (AB . AC) / |AB x AC|.
        let cot_a = corner_cotan(a, b, c).unwrap();
        prop_assert!((cot_a - cx / cy).abs() <= 1e-10 * (1.0 + (cx / cy).abs()));
        // Symmetric in the two adjacent sides.
        prop_assert_eq!(corner_cotan(a, c, b).unwrap(), cot_a);
        prop_assert_eq!(corner_angle(a, c, b).unwrap(), angle_a);
    }

    #[test]
    fn orient_faces_undoes_random_flips((shape, _) in arbitrary_mesh(), seed in any::<u64>()) {
        let mut rng = SplitMix64::new(seed);
        let n = shape.positions.len();
        let flips: Vec<bool> = shape.faces.iter().map(|_| rng.next_f64() < 0.3).collect();
        let scrambled: Vec<[u32; 3]> = shape
            .faces
            .iter()
            .zip(&flips)
            .map(|(&[a, b, c], &flip)| if flip { [a, c, b] } else { [a, b, c] })
            .collect();
        let n_flips = flips.iter().filter(|&&f| f).count();
        if n_flips != 0 && n_flips != shape.faces.len() {
            let rejected = matches!(TriMesh::new(n, scrambled.clone()), Err(MeshError::InconsistentOrientation { .. }));
            prop_assert!(rejected);
        }
        let oriented = orient_faces(n, &scrambled).unwrap();
        prop_assert_eq!(oriented.n_components, 1);
        prop_assert_eq!(oriented.n_flipped, oriented.flipped.iter().filter(|&&f| f).count());
        // The result is valid and is the original orientation or its mirror image.
        TriMesh::new(n, oriented.faces.clone()).unwrap();
        let same = oriented.flipped == flips;
        let mirrored = oriented.flipped.iter().zip(&flips).all(|(a, b)| a != b);
        prop_assert!(same || mirrored);
        // The fewest faces are flipped.
        prop_assert!(2 * oriented.n_flipped <= shape.faces.len());
        for (f, (&[a, b, c], &flip)) in scrambled.iter().zip(&oriented.flipped).enumerate() {
            prop_assert_eq!(oriented.faces[f], if flip { [a, c, b] } else { [a, b, c] });
        }
    }

    #[test]
    fn components_split_and_compact((a, _) in arbitrary_mesh(), (b, _) in arbitrary_mesh(), extra in 0usize..4) {
        // Two meshes side by side, with `extra` unreferenced vertices in between.
        let offset = (a.positions.len() + extra) as u32;
        let mut faces = a.faces.clone();
        faces.extend(b.faces.iter().map(|f| f.map(|v| v + offset)));
        let n = a.positions.len() + extra + b.positions.len();
        if extra > 0 {
            let err = TriMesh::new(n, faces.clone()).unwrap_err();
            prop_assert_eq!(err, MeshError::UnreferencedVertex { vertex: a.positions.len() as u32, count: extra });
        }
        let compact = compact_vertices(n, &faces).unwrap();
        prop_assert_eq!(compact.new_to_old.len(), n - extra);
        for (new, &old) in compact.new_to_old.iter().enumerate() {
            prop_assert_eq!(compact.old_to_new[old as usize], new as u32);
        }
        let mesh = TriMesh::new(compact.new_to_old.len(), compact.faces).unwrap();
        prop_assert_eq!(mesh.n_components(), 2);
        let (ta, tb) = (a.topology().unwrap(), b.topology().unwrap());
        prop_assert_eq!(mesh.euler_characteristic(), ta.euler_characteristic() + tb.euler_characteristic());
        prop_assert_eq!(mesh.genus(), ta.genus() + tb.genus());
        let parts = mesh.split_components();
        prop_assert_eq!(parts.len(), 2);
        prop_assert_eq!(&parts[0].mesh, &ta);
        prop_assert_eq!(&parts[1].mesh, &tb);
        for part in &parts {
            for (f, &pf) in part.face_to_parent.iter().enumerate() {
                for k in 0..3 {
                    let v = part.mesh.faces()[f][k];
                    prop_assert_eq!(part.vertex_to_parent[v as usize], mesh.faces()[pf as usize][k]);
                }
            }
        }
    }

    #[test]
    fn sparse_solves_on_mesh_matrices((shape, _) in arbitrary_mesh(), seed in any::<u64>()) {
        let mesh = shape.topology().unwrap();
        let g = Geometry::from_positions(&mesh, &shape.positions).unwrap();
        let l = cotan_laplacian(&mesh, &g).unwrap();
        let m = lumped_mass_normalized(&mesh, &g).unwrap();
        let n = mesh.n_vertices();
        let pool = ThreadPool::new(2).unwrap();

        // Positive definite: L + M.
        let mut a = l.clone();
        a.add_to_diagonal(&m).unwrap();
        let x_true = random_vector(n, seed);
        let mut b = vec![0.0; n];
        a.mul_vec(&x_true, &mut b).unwrap();
        let symbolic = Symbolic::analyze(&a, Method::Llt).unwrap();
        for par in [Parallelism::Sequential, Parallelism::Pool(&pool)] {
            let f = symbolic.factor(&a, par).unwrap();
            let mut x = b.clone();
            f.solve_in_place(&mut x, par).unwrap();
            prop_assert!(relative_residual(&a, &x, &b).unwrap() < 1e-10);
        }

        // Quasi-definite saddle point system [[-M, L], [L, D]] with D > 0.
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i as u32, i as u32, -m[i]));
            triplets.push(((n + i) as u32, (n + i) as u32, 0.01 + (i % 3) as f64));
        }
        for (i, j, v) in l.triplets() {
            triplets.push((n as u32 + i, j, v));
            triplets.push((i, n as u32 + j, v));
        }
        let k = SparseMatrix::from_triplets(2 * n, 2 * n, &triplets).unwrap();
        prop_assert_eq!(k.symmetry_defect().unwrap(), 0.0);
        let x_true = random_vector(2 * n, seed ^ 1);
        let mut b = vec![0.0; 2 * n];
        k.mul_vec(&x_true, &mut b).unwrap();
        for method in [Method::Ldlt, Method::Lu] {
            let f = Symbolic::analyze(&k, method).unwrap().factor(&k, Parallelism::Sequential).unwrap();
            let mut x = b.clone();
            f.solve_in_place(&mut x, Parallelism::Sequential).unwrap();
            let residual = relative_residual(&k, &x, &b).unwrap();
            prop_assert!(residual < 1e-8, "{method:?}: residual {residual}");
        }
    }

    #[test]
    fn triplets_agree_with_a_dense_matrix(
        nrows in 1usize..8,
        ncols in 1usize..8,
        entries in prop::collection::vec((0u32..8, 0u32..8, -4i32..5), 0..40),
        seed in any::<u64>(),
    ) {
        let triplets: Vec<(u32, u32, f64)> = entries
            .iter()
            .map(|&(i, j, v)| (i % nrows as u32, j % ncols as u32, f64::from(v)))
            .collect();
        let a = SparseMatrix::from_triplets(nrows, ncols, &triplets).unwrap();
        let mut dense = vec![vec![0.0; ncols]; nrows];
        for &(i, j, v) in &triplets {
            dense[i as usize][j as usize] += v;
        }
        for i in 0..nrows {
            for j in 0..ncols {
                prop_assert_eq!(a.get(i, j), dense[i][j]);
            }
        }
        // Structure invariants and a round trip through the raw arrays.
        let rebuilt = SparseMatrix::from_csc(nrows, ncols, a.col_ptr().to_vec(), a.row_idx().to_vec(), a.values().to_vec());
        prop_assert_eq!(rebuilt.unwrap(), a.clone());
        let x = random_vector(ncols, seed);
        let mut y = vec![0.0; nrows];
        a.mul_vec(&x, &mut y).unwrap();
        for i in 0..nrows {
            let expect: f64 = (0..ncols).map(|j| dense[i][j] * x[j]).sum();
            prop_assert!((y[i] - expect).abs() < 1e-12);
        }
        let t = a.transpose();
        for i in 0..nrows {
            for j in 0..ncols {
                prop_assert_eq!(t.get(j, i), dense[i][j]);
            }
        }
    }

    #[test]
    fn obj_round_trip(
        coords in prop::collection::vec(prop::array::uniform3(-1e6..1e6f64), 3..30),
        sizes in prop::collection::vec(3usize..7, 1..12),
        with_uvs in prop::bool::ANY,
        seed in any::<u64>(),
    ) {
        let mut rng = SplitMix64::new(seed);
        let mut corner_verts = Vec::new();
        let mut face_starts = Vec::new();
        for &size in &sizes {
            face_starts.push(corner_verts.len() as u32);
            for _ in 0..size {
                corner_verts.push((rng.next_u64() % coords.len() as u64) as u32);
            }
        }
        let uvs: Vec<[f64; 2]> = corner_verts.iter().map(|_| [rng.uniform(-3.0, 3.0), rng.next_f64()]).collect();
        let text = obj_to_string(&coords, &corner_verts, &face_starts, with_uvs.then_some(&uvs[..])).unwrap();
        let read = parse_obj(&text).unwrap();
        prop_assert_eq!(&read.positions, &coords);
        prop_assert_eq!(&read.corner_verts, &corner_verts);
        prop_assert_eq!(&read.face_starts, &face_starts);
        if with_uvs {
            prop_assert_eq!(read.corner_uvs().unwrap(), uvs);
        } else {
            prop_assert!(read.texcoords.is_empty());
            prop_assert_eq!(read.corner_uvs(), None);
        }
    }
}

#[cfg(feature = "blender-tessellation")]
mod triangulation {
    use super::*;
    use ddg_core::triangulate::triangulate_polygons;

    /// A star-shaped polygon (convex when `r_min` is 1) in a random plane, with coordinates
    /// that are exact in single precision.
    fn star_polygon(n: usize, r_min: f64, rng: &mut SplitMix64) -> Vec<[f64; 3]> {
        // Angles are kept apart so that no two corners nearly coincide.
        let mut angles: Vec<f64> = (0..n)
            .map(|i| 2.0 * PI * (i as f64 + rng.uniform(0.1, 0.9)) / n as f64)
            .collect();
        angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // An orthonormal frame from two random vectors.
        let u = loop {
            let v = [
                rng.uniform(-1.0, 1.0),
                rng.uniform(-1.0, 1.0),
                rng.uniform(-1.0, 1.0),
            ];
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if len > 0.2 {
                break v.map(|x| x / len);
            }
        };
        let w = loop {
            let v = [
                rng.uniform(-1.0, 1.0),
                rng.uniform(-1.0, 1.0),
                rng.uniform(-1.0, 1.0),
            ];
            let d = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
            let v = [v[0] - d * u[0], v[1] - d * u[1], v[2] - d * u[2]];
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if len > 0.2 {
                break v.map(|x| x / len);
            }
        };
        angles
            .iter()
            .map(|&a| {
                let r = rng.uniform(r_min, 1.0);
                let (x, y) = (r * a.cos(), r * a.sin());
                [0, 1, 2].map(|k| f64::from((x * u[k] + y * w[k]) as f32))
            })
            .collect()
    }

    fn vector_area(points: &[[f64; 3]]) -> [f64; 3] {
        let mut n = [0.0; 3];
        for i in 0..points.len() {
            let (p, q) = (points[i], points[(i + 1) % points.len()]);
            n[0] += 0.5 * (p[1] * q[2] - p[2] * q[1]);
            n[1] += 0.5 * (p[2] * q[0] - p[0] * q[2]);
            n[2] += 0.5 * (p[0] * q[1] - p[1] * q[0]);
        }
        n
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(256))]

        #[test]
        fn triangulation_preserves_area_and_corners(
            sizes in prop::collection::vec(3usize..14, 1..8),
            convex in prop::bool::ANY,
            seed in any::<u64>(),
        ) {
            let mut rng = SplitMix64::new(seed);
            let mut positions = Vec::new();
            let mut corner_verts: Vec<u32> = Vec::new();
            let mut face_starts = Vec::new();
            for &n in &sizes {
                face_starts.push(corner_verts.len() as u32);
                let start = positions.len() as u32;
                positions.extend(star_polygon(n, if convex { 1.0 } else { 0.35 }, &mut rng));
                // Corners in a rotated order, so that corner index != vertex index.
                let shift = (rng.next_u64() % n as u64) as u32;
                corner_verts.extend((0..n as u32).map(|i| start + (i + shift) % n as u32));
            }
            let t = triangulate_polygons(&positions, &corner_verts, &face_starts).unwrap();
            prop_assert_eq!(t.n_triangles(), corner_verts.len() - 2 * sizes.len());
            prop_assert_eq!(t.triangle_corners.len(), t.n_triangles());
            prop_assert_eq!(t.triangle_polygons.len(), t.n_triangles());
            prop_assert_eq!(t.is_diagonal.len(), t.n_triangles());

            let mut tri = 0;
            for (p, &n) in sizes.iter().enumerate() {
                let start = face_starts[p];
                let polygon: Vec<[f64; 3]> =
                    (0..n).map(|i| positions[corner_verts[start as usize + i] as usize]).collect();
                let normal = vector_area(&polygon);
                let polygon_area = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
                let mut area_sum = 0.0;
                // Count how often each polygon edge and each diagonal is used, per direction.
                let mut edge_uses = vec![0u32; n];
                let mut diagonal_uses = std::collections::HashMap::new();
                for _ in 0..n - 2 {
                    prop_assert_eq!(t.triangle_polygons[tri], p as u32);
                    let corners = t.triangle_corners[tri];
                    for k in 0..3 {
                        // Corner mapping: inside the polygon, and consistent with the vertices.
                        prop_assert!(corners[k] >= start && corners[k] < start + n as u32);
                        prop_assert_eq!(t.triangles[tri][k], corner_verts[corners[k] as usize]);
                        let (a, b) = (corners[k] - start, corners[(k + 1) % 3] - start);
                        if (a + 1) % n as u32 == b {
                            prop_assert!(!t.is_diagonal[tri][k]);
                            edge_uses[a as usize] += 1;
                        } else {
                            // Orientation is preserved, so a polygon edge never appears reversed.
                            prop_assert!((b + 1) % n as u32 != a);
                            prop_assert!(t.is_diagonal[tri][k]);
                            *diagonal_uses.entry((a, b)).or_insert(0u32) += 1;
                        }
                    }
                    let pts = t.triangles[tri].map(|v| positions[v as usize]);
                    let a = vector_area(&pts);
                    // Each triangle is oriented like its polygon and has positive area.
                    let along = (a[0] * normal[0] + a[1] * normal[1] + a[2] * normal[2]) / polygon_area;
                    prop_assert!(along > 0.0, "triangle {tri} is flipped or degenerate: {along}");
                    area_sum += along;
                    tri += 1;
                }
                prop_assert!((area_sum - polygon_area).abs() <= 1e-9 * polygon_area, "{area_sum} vs {polygon_area}");
                prop_assert!(edge_uses.iter().all(|&c| c == 1));
                prop_assert_eq!(diagonal_uses.len(), 2 * (n - 3));
                for (&(a, b), &count) in &diagonal_uses {
                    prop_assert_eq!(count, 1);
                    prop_assert_eq!(diagonal_uses.get(&(b, a)), Some(&1));
                }
            }
            prop_assert_eq!(tri, t.n_triangles());

            // Per-corner values travel back to the polygon corners.
            let values: Vec<[u32; 3]> = t.triangle_corners.clone();
            let back = t.to_polygon_corners(&values, corner_verts.len(), u32::MAX).unwrap();
            prop_assert_eq!(back, (0..corner_verts.len() as u32).collect::<Vec<_>>());
        }

        /// A grid of (non-planar) quads becomes a valid disk, whose interior edges are exactly
        /// the quad edges plus one diagonal per quad.
        #[test]
        fn quad_grid_becomes_a_valid_mesh(nx in 1usize..7, ny in 1usize..7, seed in any::<u64>(), bend in 0.0..0.4f64) {
            let mut rng = SplitMix64::new(seed);
            let mut positions = Vec::new();
            for i in 0..=nx {
                for j in 0..=ny {
                    positions.push([i as f64, j as f64, bend * rng.uniform(-1.0, 1.0)]);
                }
            }
            let index = |i: usize, j: usize| (i * (ny + 1) + j) as u32;
            let mut corner_verts = Vec::new();
            let mut face_starts = Vec::new();
            for i in 0..nx {
                for j in 0..ny {
                    face_starts.push(corner_verts.len() as u32);
                    corner_verts.extend([index(i, j), index(i + 1, j), index(i + 1, j + 1), index(i, j + 1)]);
                }
            }
            let t = triangulate_polygons(&positions, &corner_verts, &face_starts).unwrap();
            let mesh = TriMesh::new(positions.len(), t.triangles.clone()).unwrap();
            prop_assert_eq!(mesh.euler_characteristic(), 1);
            prop_assert_eq!(mesh.boundary_loops().len(), 1);
            let mut diagonal_edges = std::collections::HashSet::new();
            for (tri, flags) in t.is_diagonal.iter().enumerate() {
                for k in 0..3 {
                    let h = 3 * tri + k;
                    if flags[k] {
                        diagonal_edges.insert(mesh.edge(h));
                        // Both sides of a diagonal are flagged and belong to the same quad.
                        let twin = mesh.twin(h).unwrap();
                        prop_assert!(t.is_diagonal[twin / 3][twin % 3]);
                        prop_assert_eq!(t.triangle_polygons[twin / 3], t.triangle_polygons[tri]);
                    }
                }
            }
            prop_assert_eq!(diagonal_edges.len(), nx * ny);
            Geometry::from_positions(&mesh, &positions).unwrap();
        }
    }
}
