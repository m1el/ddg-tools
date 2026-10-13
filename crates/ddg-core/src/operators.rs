//! Discrete operators on a triangle mesh: cotan Laplacian, lumped mass, angle defect and boundary
//! geodesic curvature (Sec. 5.1 of Soliman, Slepčev, Crane, *Optimal Cone Singularities for
//! Conformal Flattening*, 2018), with the sign and normalization conventions of the repository's
//! `CLAUDE.md`:
//!
//! * the Laplacian `L` (and its per-vertex weighted variant) is symmetric **positive semidefinite** with zero row sums;
//! * the lumped mass is offered both raw ([`lumped_mass`], `Σ m_i` = total area) and area
//!   normalized ([`lumped_mass_normalized`], `Σ m_i = 1`) — normalization is the caller's
//!   explicit choice;
//! * the angle defect `Ω` lives on interior vertices and the geodesic curvature `k` on boundary
//!   vertices; each is returned as a full per-vertex vector that is zero where it is not defined.
//!
//! Every function takes the mesh and a [`Geometry`] built for it and fails with
//! [`GeometryError::MeshMismatch`] if the two do not fit together.

use std::f64::consts::PI;

use crate::geometry::{Geometry, GeometryError};
use crate::sparse::SparseMatrix;
use crate::topology::TriMesh;

/// The cotan Laplacian as an `n_vertices × n_vertices` sparse matrix:
/// `L_ij = −½ Σ cot θ` over the one or two corners opposite edge `ij`, and
/// `L_ii = −Σ_{j≠i} L_ij`.
///
/// The result is exactly symmetric, stores both triangles, and its pattern is the vertex
/// adjacency plus the diagonal — it depends only on the mesh, never on the values (an edge
/// whose two cotangents cancel is stored as an explicit zero).
///
/// `L` is the stiffness matrix of piecewise linear finite elements, `xᵀ L x` being the
/// Dirichlet energy of the interpolant of `x`; it is therefore positive semidefinite for every
/// valid geometry (individual off-diagonal entries can still be positive on non-Delaunay
/// meshes), and its kernel consists of the vectors that are constant on each connected
/// component. It does not change when all lengths are scaled by a common factor.
///
/// [formula: cotan-laplacian]
pub fn cotan_laplacian(mesh: &TriMesh, geometry: &Geometry) -> Result<SparseMatrix, GeometryError> {
    geometry.check_mesh(mesh)?;
    Ok(assemble_cotan(mesh, geometry, |_| 1.0))
}

/// The cotan Laplacian with a per-vertex weight `w`, as the sum of the stiffness matrices of the
/// single faces, each scaled by the mean of `w` over its three vertices:
///
/// ```text
/// L_w = Σ_f w̄_f · L_f,    w̄_f = (w_i + w_j + w_k) / 3,
/// ```
///
/// where `L_f` is the cotan stiffness matrix of face `f` alone (the contribution of `f` to
/// [`cotan_laplacian`]). This is the exact finite element matrix of `∫ w ∇f · ∇g` for
/// piecewise linear `w`, `f`, `g`: the gradients are constant on a face and the integral of
/// the linear function `w` over a triangle is its area times the mean of its vertex values.
/// Since each `L_f` is positive semidefinite (it is the Dirichlet form of the interpolant on
/// one triangle), `L_w` is positive semidefinite for every `w ≥ 0`, obtuse triangles
/// included; the rule "multiply the cotangent of edge `ij` by `(w_i + w_j)/2`" has no such
/// guarantee. With `w ≡ c` the result is `c L`.
///
/// Storage, symmetry (exact), zero row sums, sparsity pattern and explicit zeros are those of
/// [`cotan_laplacian`].
///
/// # Errors
/// [`GeometryError::MeshMismatch`] for a geometry of another mesh,
/// [`GeometryError::WeightCountMismatch`] unless `w` has one entry per vertex, and
/// [`GeometryError::InvalidWeight`] for the first weight that is negative, NaN or infinite.
///
/// [formula: cut-weighted-laplacian]
pub fn cotan_laplacian_weighted(
    mesh: &TriMesh,
    geometry: &Geometry,
    w: &[f64],
) -> Result<SparseMatrix, GeometryError> {
    geometry.check_mesh(mesh)?;
    if w.len() != mesh.n_vertices() {
        return Err(GeometryError::WeightCountMismatch {
            expected: mesh.n_vertices(),
            got: w.len(),
        });
    }
    if let Some(vertex) = w.iter().position(|x| !(x.is_finite() && *x >= 0.0)) {
        return Err(GeometryError::InvalidWeight {
            vertex: vertex as u32,
            weight: w[vertex],
        });
    }
    let faces = mesh.faces();
    Ok(assemble_cotan(mesh, geometry, |f| {
        let [a, b, c] = faces[f];
        (w[a as usize] + w[b as usize] + w[c as usize]) / 3.0
    }))
}

/// Assembles `Σ_f scale(f) · L_f` (see [`cotan_laplacian`]); the mesh and geometry must match.
fn assemble_cotan(
    mesh: &TriMesh,
    geometry: &Geometry,
    scale: impl Fn(usize) -> f64,
) -> SparseMatrix {
    let n = mesh.n_vertices();
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(6 * mesh.n_faces() + n);
    for v in 0..n as u32 {
        triplets.push((v, v, 0.0));
    }
    for (f, (face, cotans)) in mesh
        .faces()
        .iter()
        .zip(geometry.corner_cotans())
        .enumerate()
    {
        let s = scale(f);
        for k in 0..3 {
            // The edge opposite corner k joins corners k + 1 and k + 2.
            let (i, j) = (face[(k + 1) % 3], face[(k + 2) % 3]);
            let w = 0.5 * s * cotans[k];
            triplets.push((i, j, -w));
            triplets.push((j, i, -w));
        }
    }
    let mut l = SparseMatrix::from_triplets(n, n, &triplets)
        .expect("vertex indices of a valid mesh are in range");
    // Diagonal: minus the sum of the off-diagonal entries of the column (= row, by symmetry).
    let col_ptr = l.col_ptr().to_vec();
    for j in 0..n {
        let range = col_ptr[j]..col_ptr[j + 1];
        let mut sum = 0.0;
        let mut diagonal_slot = range.start;
        for k in range {
            if l.row_idx()[k] as usize == j {
                diagonal_slot = k;
            } else {
                sum += l.values()[k];
            }
        }
        l.values_mut()[diagonal_slot] = -sum;
    }
    l
}

/// The lumped (barycentric) vertex masses `m_i = ⅓ Σ A_f` over the faces `f` at vertex `i`.
/// They sum to the total area and scale with the square of the length unit.
///
/// [formula: lumped-mass]
pub fn lumped_mass(mesh: &TriMesh, geometry: &Geometry) -> Result<Vec<f64>, GeometryError> {
    geometry.check_mesh(mesh)?;
    let mut mass = vec![0.0; mesh.n_vertices()];
    for (face, &area) in mesh.faces().iter().zip(geometry.face_areas()) {
        for &v in face {
            mass[v as usize] += area / 3.0;
        }
    }
    Ok(mass)
}

/// The lumped masses divided by the total area, `m_i / Σ_f A_f`, so that they sum to 1 and do
/// not depend on the length unit. This is the normalization under which the cone-count
/// parameter `λ` of the optimal cones paper is dimensionless.
///
/// [formula: lumped-mass]
pub fn lumped_mass_normalized(
    mesh: &TriMesh,
    geometry: &Geometry,
) -> Result<Vec<f64>, GeometryError> {
    let mut mass = lumped_mass(mesh, geometry)?;
    let total = geometry.total_area();
    for m in &mut mass {
        *m /= total;
    }
    Ok(mass)
}

/// For every vertex, the sum of the corner angles of its incident faces.
pub fn vertex_angle_sums(mesh: &TriMesh, geometry: &Geometry) -> Result<Vec<f64>, GeometryError> {
    geometry.check_mesh(mesh)?;
    let mut sums = vec![0.0; mesh.n_vertices()];
    for (face, angles) in mesh.faces().iter().zip(geometry.corner_angles()) {
        for k in 0..3 {
            sums[face[k] as usize] += angles[k];
        }
    }
    Ok(sums)
}

/// The angle defect (integrated Gaussian curvature) `Ω_i = 2π − Σ θ` at every interior vertex
/// `i`, where the sum runs over the corner angles at `i`. Entries of boundary vertices are `0.0`.
///
/// [formula: angle-defect]
pub fn angle_defects(mesh: &TriMesh, geometry: &Geometry) -> Result<Vec<f64>, GeometryError> {
    let mut out = vertex_angle_sums(mesh, geometry)?;
    for (v, x) in out.iter_mut().enumerate() {
        *x = if mesh.is_boundary_vertex(v) {
            0.0
        } else {
            2.0 * PI - *x
        };
    }
    Ok(out)
}

/// The integrated geodesic curvature `k_i = π − Σ θ` at every boundary vertex `i`, where the
/// sum runs over the corner angles at `i`. Entries of interior vertices are `0.0`.
///
/// [formula: boundary-geodesic-curvature]
pub fn boundary_geodesic_curvatures(
    mesh: &TriMesh,
    geometry: &Geometry,
) -> Result<Vec<f64>, GeometryError> {
    let mut out = vertex_angle_sums(mesh, geometry)?;
    for (v, x) in out.iter_mut().enumerate() {
        *x = if mesh.is_boundary_vertex(v) {
            PI - *x
        } else {
            0.0
        };
    }
    Ok(out)
}

/// The discrete Gauss–Bonnet residual `Σ_I Ω_i + Σ_B k_i − 2π χ`, which is zero in exact
/// arithmetic for every valid geometry. Its size measures accumulated rounding error.
///
/// [formula: gauss-bonnet]
pub fn gauss_bonnet_residual(mesh: &TriMesh, geometry: &Geometry) -> Result<f64, GeometryError> {
    let omega: f64 = angle_defects(mesh, geometry)?.iter().sum();
    let kappa: f64 = boundary_geodesic_curvatures(mesh, geometry)?.iter().sum();
    Ok(omega + kappa - 2.0 * PI * mesh.euler_characteristic() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generators;

    fn square() -> (TriMesh, Geometry) {
        let mesh = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let g = Geometry::from_positions(&mesh, &p).unwrap();
        (mesh, g)
    }

    #[test]
    fn laplacian_of_a_square() {
        let (mesh, g) = square();
        let l = cotan_laplacian(&mesh, &g).unwrap();
        // Right triangles: the diagonal 0-2 is opposite two right angles (weight 0), each side
        // is opposite one 45 degree angle (weight 1/2).
        let expect = [
            [1.0, -0.5, 0.0, -0.5],
            [-0.5, 1.0, -0.5, 0.0],
            [0.0, -0.5, 1.0, -0.5],
            [-0.5, 0.0, -0.5, 1.0],
        ];
        for (i, row) in expect.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                assert!(
                    (l.get(i, j) - value).abs() < 1e-15,
                    "L[{i}][{j}] = {}",
                    l.get(i, j)
                );
            }
        }
        // Pattern: adjacency plus diagonal. 0-2 is an edge (stored, ~0); 1-3 is not.
        assert!(l.find(0, 2).is_some() && l.find(2, 0).is_some());
        assert!(l.find(1, 3).is_none());
        assert_eq!(l.nnz(), 4 + 2 * 5);
        assert_eq!(l.symmetry_defect().unwrap(), 0.0);
    }

    #[test]
    fn mass_and_curvature_of_a_square() {
        let (mesh, g) = square();
        let m = lumped_mass(&mesh, &g).unwrap();
        let expect = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 3.0, 1.0 / 6.0];
        for i in 0..4 {
            assert!((m[i] - expect[i]).abs() < 1e-15);
        }
        let mn = lumped_mass_normalized(&mesh, &g).unwrap();
        assert!((mn.iter().sum::<f64>() - 1.0).abs() < 1e-15);
        let sums = vertex_angle_sums(&mesh, &g).unwrap();
        for s in sums {
            assert!((s - PI / 2.0).abs() < 1e-15);
        }
        // All four vertices are on the boundary: Ω = 0 everywhere, k = π/2 at each corner.
        assert_eq!(angle_defects(&mesh, &g).unwrap(), vec![0.0; 4]);
        for k in boundary_geodesic_curvatures(&mesh, &g).unwrap() {
            assert!((k - PI / 2.0).abs() < 1e-15);
        }
        assert!(gauss_bonnet_residual(&mesh, &g).unwrap().abs() < 1e-14);
    }

    #[test]
    fn icosahedron_curvature() {
        let ico = generators::icosphere(0).unwrap();
        let mesh = ico.topology().unwrap();
        let g = Geometry::from_positions(&mesh, &ico.positions).unwrap();
        // Five equilateral triangles meet at every vertex.
        for omega in angle_defects(&mesh, &g).unwrap() {
            assert!((omega - (2.0 * PI - 5.0 * PI / 3.0)).abs() < 1e-14);
        }
        assert_eq!(
            boundary_geodesic_curvatures(&mesh, &g).unwrap(),
            vec![0.0; 12]
        );
        assert!(gauss_bonnet_residual(&mesh, &g).unwrap().abs() < 1e-13);
        // Each edge has two opposite 60 degree angles: weight cot(60°) = 1/√3, degree 5.
        let l = cotan_laplacian(&mesh, &g).unwrap();
        let w = 1.0 / 3f64.sqrt();
        for i in 0..12 {
            assert!((l.get(i, i) - 5.0 * w).abs() < 1e-14);
            for j in mesh.vertex_neighbors(i) {
                assert!((l.get(i, j) + w).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn scale_invariance() {
        let base = generators::bumpy_grid_disk(5, 4, 0.2).unwrap();
        let mesh = base.topology().unwrap();
        let g1 = Geometry::from_positions(&mesh, &base.positions).unwrap();
        let s = 8.0; // a power of two: scaling is exact in floating point
        let scaled: Vec<[f64; 3]> = base.positions.iter().map(|p| p.map(|x| s * x)).collect();
        let g2 = Geometry::from_positions(&mesh, &scaled).unwrap();
        assert_eq!(
            cotan_laplacian(&mesh, &g1).unwrap(),
            cotan_laplacian(&mesh, &g2).unwrap()
        );
        assert_eq!(
            angle_defects(&mesh, &g1).unwrap(),
            angle_defects(&mesh, &g2).unwrap()
        );
        assert_eq!(
            lumped_mass_normalized(&mesh, &g1).unwrap(),
            lumped_mass_normalized(&mesh, &g2).unwrap()
        );
        let (m1, m2) = (
            lumped_mass(&mesh, &g1).unwrap(),
            lumped_mass(&mesh, &g2).unwrap(),
        );
        for (a, b) in m1.iter().zip(&m2) {
            assert_eq!(a * s * s, *b);
        }
    }

    #[test]
    fn mismatched_geometry_is_an_error() {
        let (mesh, _) = square();
        let other = generators::grid_disk(2, 2).unwrap();
        let other_mesh = other.topology().unwrap();
        let g = Geometry::from_positions(&other_mesh, &other.positions).unwrap();
        let e = cotan_laplacian(&mesh, &g).unwrap_err();
        assert!(matches!(e, GeometryError::MeshMismatch { .. }));
        assert!(lumped_mass(&mesh, &g).is_err());
        assert!(angle_defects(&mesh, &g).is_err());
        assert!(boundary_geodesic_curvatures(&mesh, &g).is_err());
        assert!(gauss_bonnet_residual(&mesh, &g).is_err());
    }
}
