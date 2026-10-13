//! Curvature per face (SPEC §3.1) and the curvature-based weight of the scaled Dirichlet energy
//! (SPEC §3.6).

use std::f64::consts::PI;

use ddg_core::operators::vertex_angle_sums;
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::{Geometry, GeometryError, SparseMatrix, TriMesh};

use crate::error::CutError;

/// The integrated Gaussian curvature `Ω_f` of every face (paper Sec. 5.2, Eq. 24 with the sign
/// of SPEC D1): with the corner angles of an interior vertex `i` rescaled to sum to `2π`,
/// `θ̃_i^f = θ_i^f · 2π/Θ_i`, and those of a boundary vertex left alone,
///
/// ```text
/// Ω_f = θ̃_i^f + θ̃_j^f + θ̃_k^f − π
/// ```
///
/// It is positive on a sphere and scale invariant. Because the corner angles of a Euclidean
/// triangle sum to `π`, the same number is `Σ_corners θ (2π − Θ)/Θ` over the interior corners of
/// the face; that form is evaluated, since it has no cancellation where the surface is flat.
///
/// # Errors
/// [`GeometryError::MeshMismatch`].
///
/// [formula: cut-face-curvature]
pub fn face_curvatures(mesh: &TriMesh, geometry: &Geometry) -> Result<Vec<f64>, GeometryError> {
    let sums = vertex_angle_sums(mesh, geometry)?;
    Ok(face_curvatures_from(mesh, geometry.corner_angles(), &sums))
}

/// [`face_curvatures`] from corner angles and vertex angle sums that are already at hand.
pub(crate) fn face_curvatures_from(
    mesh: &TriMesh,
    corner_angles: &[[f64; 3]],
    angle_sums: &[f64],
) -> Vec<f64> {
    mesh.faces()
        .iter()
        .zip(corner_angles)
        .map(|(face, angles)| {
            let mut omega = 0.0;
            for k in 0..3 {
                let v = face[k] as usize;
                if !mesh.is_boundary_vertex(v) {
                    omega += angles[k] * (2.0 * PI - angle_sums[v]) / angle_sums[v];
                }
            }
            omega
        })
        .collect()
}

/// The residual `Σ_f Ω_f − Σ_{i ∉ ∂M} (2π − Θ_i)` of the identity that makes the face curvatures
/// a redistribution of the vertex angle defects. It is zero in exact arithmetic; on a closed
/// mesh the second sum is `2πχ`.
///
/// # Errors
/// [`GeometryError::MeshMismatch`].
///
/// [formula: cut-face-curvature-sum]
pub fn face_curvature_sum_residual(
    mesh: &TriMesh,
    geometry: &Geometry,
) -> Result<f64, GeometryError> {
    let sums = vertex_angle_sums(mesh, geometry)?;
    let faces: f64 = face_curvatures_from(mesh, geometry.corner_angles(), &sums)
        .iter()
        .sum();
    let vertices: f64 = sums
        .iter()
        .enumerate()
        .filter(|&(v, _)| !mesh.is_boundary_vertex(v))
        .map(|(_, theta)| 2.0 * PI - theta)
        .sum();
    Ok(faces - vertices)
}

/// The curvature-based weight `a = r³` of the scaled Dirichlet energy (paper Sec. 3.5.1, 5.6;
/// SPEC §3.6, D5), on the unit-area surface:
///
/// ```text
/// K_i  = (2π − Θ_i)/A_i  (interior),  0 on ∂M
/// r⁰_i = 1 / ( √|K_i| + 1/r_ε )
/// (W_V + c L) r = W_V r⁰
/// a_i  = r_i³
/// ```
///
/// `vertex_area` are the lumped vertex areas `A_i` of the unit-area surface, `face_cot` the
/// corner cotangents of every face.
///
/// [formula: cut-feature-size]
pub(crate) fn feature_size_weight(
    mesh: &TriMesh,
    face_cot: &[[f64; 3]],
    vertex_area: &[f64],
    angle_sums: &[f64],
    r_eps: f64,
    smoothing: f64,
    par: Parallelism<'_>,
) -> Result<Vec<f64>, CutError> {
    crate::problem::validate_feature_size(r_eps, smoothing)?;
    let n = mesh.n_vertices();
    let mut r: Vec<f64> = (0..n)
        .map(|v| {
            let k = if mesh.is_boundary_vertex(v) {
                0.0
            } else {
                (2.0 * PI - angle_sums[v]) / vertex_area[v]
            };
            1.0 / (k.abs().sqrt() + 1.0 / r_eps)
        })
        .collect();
    if smoothing > 0.0 {
        // (W_V + c L), lower triangle; duplicates are summed by `from_triplets`.
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(n + 9 * mesh.n_faces());
        for (v, &area) in vertex_area.iter().enumerate() {
            triplets.push((v as u32, v as u32, area));
        }
        for (face, cot) in mesh.faces().iter().zip(face_cot) {
            for k in 0..3 {
                let (i, j) = (face[(k + 1) % 3], face[(k + 2) % 3]);
                let w = 0.5 * smoothing * cot[k];
                triplets.push((i, i, w));
                triplets.push((j, j, w));
                triplets.push((i.max(j), i.min(j), -w));
            }
        }
        let matrix = SparseMatrix::from_triplets(n, n, &triplets)?;
        let factorization = Symbolic::analyze(&matrix, Method::Llt)?.factor(&matrix, par)?;
        for (x, &area) in r.iter_mut().zip(vertex_area) {
            *x *= area;
        }
        factorization.solve_in_place(&mut r, par)?;
    }
    if let Some(v) = r.iter().position(|x| !(x.is_finite() && *x > 0.0)) {
        return Err(CutError::NonPositiveFeatureSize {
            vertex: v as u32,
            value: r[v],
        });
    }
    Ok(r.into_iter().map(|x| x * x * x).collect())
}
