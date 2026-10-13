//! Tests of `cotan_laplacian_weighted`: constants, symmetry, row sums, positive semidefiniteness
//! on obtuse meshes, and agreement with the face-by-face Dirichlet energy.

#![allow(clippy::needless_range_loop)]

use ddg_core::generators::{self, GeneratedMesh, SplitMix64};
use ddg_core::geometry::GeometryError;
use ddg_core::operators::{cotan_laplacian, cotan_laplacian_weighted};
use ddg_core::{Geometry, SparseMatrix, TriMesh};

fn build(m: &GeneratedMesh) -> (TriMesh, Geometry) {
    let mesh = m.topology().unwrap();
    let geometry = Geometry::from_positions(&mesh, &m.positions).unwrap();
    (mesh, geometry)
}

fn sheared(n: usize, s: f64) -> GeneratedMesh {
    let mut m = generators::grid_disk(n, n).unwrap();
    for p in &mut m.positions {
        p[0] += s * p[1];
    }
    m
}

fn meshes() -> Vec<(&'static str, GeneratedMesh)> {
    let mut bumpy = generators::bumpy_grid_disk(9, 7, 0.4).unwrap();
    generators::perturb(&mut bumpy.positions, [0.01, 0.01, 0.01], 5);
    vec![
        ("icosphere", generators::icosphere(2).unwrap()),
        ("sheared 3", sheared(8, 3.0)),
        ("sheared 6", sheared(6, 6.0)),
        ("bumpy", bumpy),
    ]
}

fn quad_form(l: &SparseMatrix, x: &[f64]) -> f64 {
    let n = x.len();
    let mut s = 0.0;
    for j in 0..n {
        for k in l.col_ptr()[j]..l.col_ptr()[j + 1] {
            s += x[l.row_idx()[k] as usize] * l.values()[k] * x[j];
        }
    }
    s
}

fn random_vec(rng: &mut SplitMix64, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    (0..n).map(|_| rng.uniform(lo, hi)).collect()
}

#[test]
fn constant_weight_gives_scaled_laplacian() {
    for (name, m) in meshes() {
        let (mesh, geometry) = build(&m);
        let l = cotan_laplacian(&mesh, &geometry).unwrap();
        for c in [0.0, 1.0, 0.37, 12.5] {
            let lw =
                cotan_laplacian_weighted(&mesh, &geometry, &vec![c; mesh.n_vertices()]).unwrap();
            // Same pattern, entries equal to c L to rounding.
            assert_eq!(lw.col_ptr(), l.col_ptr(), "{name}: pattern");
            assert_eq!(lw.row_idx(), l.row_idx(), "{name}: pattern");
            let scale = l.values().iter().fold(0.0f64, |a, x| a.max(x.abs()));
            for (a, b) in lw.values().iter().zip(l.values()) {
                assert!(
                    (a - c * b).abs() <= 1e-14 * c.max(1.0) * scale,
                    "{name} c = {c}"
                );
            }
        }
        // w = 1 reproduces the unweighted matrix.
        let one =
            cotan_laplacian_weighted(&mesh, &geometry, &vec![1.0; mesh.n_vertices()]).unwrap();
        for (a, b) in one.values().iter().zip(l.values()) {
            assert!((a - b).abs() <= 1e-15 * (1.0 + b.abs()));
        }
    }
}

#[test]
fn symmetric_with_zero_row_sums_and_explicit_zeros() {
    for (name, m) in meshes() {
        let (mesh, geometry) = build(&m);
        let mut rng = SplitMix64::new(1);
        let w = random_vec(&mut rng, mesh.n_vertices(), 0.0, 3.0);
        let lw = cotan_laplacian_weighted(&mesh, &geometry, &w).unwrap();
        assert_eq!(lw.symmetry_defect().unwrap(), 0.0, "{name}");
        let l = cotan_laplacian(&mesh, &geometry).unwrap();
        assert_eq!(lw.nnz(), l.nnz(), "{name}");
        let mut row_sum = vec![0.0; mesh.n_vertices()];
        let mut row_abs = vec![0.0; mesh.n_vertices()];
        for j in 0..mesh.n_vertices() {
            for k in lw.col_ptr()[j]..lw.col_ptr()[j + 1] {
                row_sum[lw.row_idx()[k] as usize] += lw.values()[k];
                row_abs[lw.row_idx()[k] as usize] += lw.values()[k].abs();
            }
        }
        for i in 0..mesh.n_vertices() {
            assert!(
                row_sum[i].abs() <= 1e-14 * row_abs[i].max(1.0),
                "{name}: row {i} sums to {}",
                row_sum[i]
            );
        }
        // Zero weights everywhere keep the pattern, with explicit zeros.
        let zero =
            cotan_laplacian_weighted(&mesh, &geometry, &vec![0.0; mesh.n_vertices()]).unwrap();
        assert_eq!(zero.nnz(), l.nnz());
        assert!(zero.values().iter().all(|&x| x == 0.0 || x == -0.0));
    }
}

#[test]
fn positive_semidefinite_for_nonnegative_weights() {
    for (name, m) in meshes() {
        let (mesh, geometry) = build(&m);
        let n = mesh.n_vertices();
        let mut rng = SplitMix64::new(2);
        let mut worst = f64::INFINITY;
        for trial in 0..200 {
            // Mix of dense weights and weights that vanish on most of the mesh.
            let mut w = random_vec(&mut rng, n, 0.0, 5.0);
            if trial % 3 == 0 {
                for x in &mut w {
                    if rng.next_f64() < 0.7 {
                        *x = 0.0;
                    }
                }
            }
            let x = random_vec(&mut rng, n, -1.0, 1.0);
            let lw = cotan_laplacian_weighted(&mesh, &geometry, &w).unwrap();
            let q = quad_form(&lw, &x);
            // Rounding scale: sum of |entries| times |x|².
            let scale: f64 = lw.values().iter().map(|v| v.abs()).sum::<f64>() * 1.0;
            worst = worst.min(q / scale.max(1e-300));
            assert!(
                q >= -1e-13 * scale,
                "{name}: x^T L_w x = {q} (scale {scale})"
            );
        }
        println!("{name}: smallest x^T L_w x / sum|L_w| over 200 trials = {worst:.3e}");
    }
}

#[test]
fn quadratic_form_is_the_weighted_face_by_face_dirichlet_energy() {
    for (name, m) in meshes() {
        let (mesh, geometry) = build(&m);
        let n = mesh.n_vertices();
        let mut rng = SplitMix64::new(3);
        let w = random_vec(&mut rng, n, 0.0, 4.0);
        let x = random_vec(&mut rng, n, -2.0, 2.0);
        let lw = cotan_laplacian_weighted(&mesh, &geometry, &w).unwrap();
        // Dirichlet energy of the linear interpolant on each face, from 3D positions:
        // |grad f|^2 area, with grad f = sum_i f_i (n x e_i) / (2 A) on a triangle.
        let mut energy = 0.0;
        let mut scale = 0.0;
        for (f, face) in mesh.faces().iter().enumerate() {
            let p = face.map(|v| m.positions[v as usize]);
            let area = geometry.face_areas()[f];
            let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            let cross = |a: [f64; 3], b: [f64; 3]| {
                [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ]
            };
            let nrm = cross(sub(p[1], p[0]), sub(p[2], p[0])); // 2 A * unit normal
            let mut grad = [0.0; 3];
            for c in 0..3 {
                // Edge opposite corner c, counter-clockwise.
                let e = sub(p[(c + 2) % 3], p[(c + 1) % 3]);
                let r = cross(nrm, e);
                let fc = x[face[c] as usize];
                for d in 0..3 {
                    grad[d] += fc * r[d] / (2.0 * area * 2.0 * area);
                }
            }
            let g2 = grad.iter().map(|g| g * g).sum::<f64>();
            let wbar = face.iter().map(|&v| w[v as usize]).sum::<f64>() / 3.0;
            energy += wbar * g2 * area;
            scale += wbar * g2.abs() * area;
        }
        let q = quad_form(&lw, &x);
        assert!(
            (q - energy).abs() <= 1e-12 * scale.max(1e-300),
            "{name}: {q} vs {energy}"
        );
    }
}

#[test]
fn bad_weights_are_errors() {
    let m = generators::grid_disk(2, 2).unwrap();
    let (mesh, geometry) = build(&m);
    let n = mesh.n_vertices();
    let e = cotan_laplacian_weighted(&mesh, &geometry, &vec![1.0; n - 1]).unwrap_err();
    assert_eq!(
        e,
        GeometryError::WeightCountMismatch {
            expected: n,
            got: n - 1
        }
    );
    assert!(e.to_string().contains(&n.to_string()));
    for bad in [-1e-9, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut w = vec![1.0; n];
        w[4] = bad;
        let e = cotan_laplacian_weighted(&mesh, &geometry, &w).unwrap_err();
        assert!(
            matches!(e, GeometryError::InvalidWeight { vertex: 4, .. }),
            "{e:?}"
        );
        assert!(e.to_string().contains("vertex 4"));
    }
    // A geometry of another mesh.
    let other = generators::grid_disk(3, 3).unwrap();
    let (_, other_geometry) = build(&other);
    assert!(matches!(
        cotan_laplacian_weighted(&mesh, &other_geometry, &vec![1.0; n]),
        Err(GeometryError::MeshMismatch { .. })
    ));
}
