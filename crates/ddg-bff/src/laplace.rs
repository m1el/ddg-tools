//! The cotan Laplacian in the forms this crate needs: applied to a vector, and assembled for a
//! subset of the vertices, on a mesh whose faces carry given corner cotangents.
//!
//! A cut mesh has the faces and corners of the mesh it was cut from, so its Laplacian is
//! assembled from the *original* corner cotangents and the vertex indices of the cut mesh; no
//! geometry has to be rebuilt after cutting.

use ddg_core::sparse::{Factorization, Method, Parallelism, SparseMatrix, Symbolic};
use ddg_core::topology::NONE;

use crate::error::BffError;

/// Computes `out = L x` for the cotan Laplacian `L` of the faces with the given corner
/// cotangents: `(L x)_i = Σ_j ½ (cot α_ij + cot β_ij) (x_i − x_j)`.
///
/// `cotans[f][k]` belongs to the corner of face `f` at `faces[f][k]` and weights the edge
/// opposite that corner. `L` is symmetric positive semidefinite.
///
/// [formula: cotan-laplacian]
pub(crate) fn laplacian_apply(faces: &[[u32; 3]], cotans: &[[f64; 3]], x: &[f64], out: &mut [f64]) {
    out.fill(0.0);
    for (face, cot) in faces.iter().zip(cotans) {
        for k in 0..3 {
            let (i, j) = (face[(k + 1) % 3] as usize, face[(k + 2) % 3] as usize);
            let flow = 0.5 * cot[k] * (x[i] - x[j]);
            out[i] += flow;
            out[j] -= flow;
        }
    }
}

/// Assembles the principal submatrix `L_FF + shift · I` of the cotan Laplacian for the vertex
/// subset `F` given by `free_index`: `free_index[v]` is the row of vertex `v`, or [`NONE`] if
/// `v` is not in `F`. `n_free` is the number of rows.
///
/// [formula: cotan-laplacian]
pub(crate) fn laplacian_block(
    faces: &[[u32; 3]],
    cotans: &[[f64; 3]],
    free_index: &[u32],
    n_free: usize,
    shift: f64,
) -> SparseMatrix {
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(12 * faces.len() + n_free);
    for i in 0..n_free as u32 {
        triplets.push((i, i, shift));
    }
    for (face, cot) in faces.iter().zip(cotans) {
        for k in 0..3 {
            let (i, j) = (
                free_index[face[(k + 1) % 3] as usize],
                free_index[face[(k + 2) % 3] as usize],
            );
            let w = 0.5 * cot[k];
            if i != NONE {
                triplets.push((i, i, w));
            }
            if j != NONE {
                triplets.push((j, j, w));
            }
            if i != NONE && j != NONE {
                triplets.push((i, j, -w));
                triplets.push((j, i, -w));
            }
        }
    }
    SparseMatrix::from_triplets(n_free, n_free, &triplets)
        .expect("row indices are below n_free by construction")
}

/// Cholesky factorization of a symmetric positive definite matrix, on the calling thread.
pub(crate) fn factor_spd(
    system: &'static str,
    matrix: &SparseMatrix,
) -> Result<Factorization, BffError> {
    let solve = |error| BffError::Solve { system, error };
    Symbolic::analyze(matrix, Method::Llt)
        .map_err(solve)?
        .factor(matrix, Parallelism::Sequential)
        .map_err(solve)
}

/// Solves with `nrhs` right-hand sides stored column after column, in place.
pub(crate) fn solve_spd(
    system: &'static str,
    factorization: &Factorization,
    rhs: &mut [f64],
    nrhs: usize,
) -> Result<(), BffError> {
    factorization
        .solve_many_in_place(rhs, nrhs, Parallelism::Sequential)
        .map_err(|error| BffError::Solve { system, error })
}

/// `‖r‖₂ / ‖b‖₂`, or `‖r‖₂` if `b = 0`.
pub(crate) fn relative_norm(r: impl Iterator<Item = f64>, b: impl Iterator<Item = f64>) -> f64 {
    let rr: f64 = r.map(|x| x * x).sum();
    let bb: f64 = b.map(|x| x * x).sum();
    if bb == 0.0 {
        rr.sqrt()
    } else {
        (rr / bb).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::generators;
    use ddg_core::geometry::Geometry;
    use ddg_core::operators::cotan_laplacian;

    #[test]
    fn apply_and_block_agree_with_ddg_core() {
        let g = generators::bumpy_grid_disk(5, 4, 0.2).unwrap();
        let mesh = g.topology().unwrap();
        let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
        let l = cotan_laplacian(&mesh, &geometry).unwrap();
        let n = mesh.n_vertices();
        let x: Vec<f64> = (0..n).map(|i| (0.7 * i as f64).sin()).collect();
        let mut expect = vec![0.0; n];
        l.mul_vec(&x, &mut expect).unwrap();
        let mut got = vec![0.0; n];
        laplacian_apply(mesh.faces(), geometry.corner_cotans(), &x, &mut got);
        for i in 0..n {
            assert!((got[i] - expect[i]).abs() < 1e-12, "row {i}");
        }
        // Interior block.
        let interior = mesh.interior_vertices();
        let mut index = vec![NONE; n];
        for (row, &v) in interior.iter().enumerate() {
            index[v as usize] = row as u32;
        }
        let block = laplacian_block(
            mesh.faces(),
            geometry.corner_cotans(),
            &index,
            interior.len(),
            0.25,
        );
        let expect = l.submatrix(&interior, &interior).unwrap();
        for (a, &i) in interior.iter().enumerate() {
            for (b, &j) in interior.iter().enumerate() {
                let shift = if a == b { 0.25 } else { 0.0 };
                let want = l.get(i as usize, j as usize) + shift;
                assert!((block.get(a, b) - want).abs() < 1e-12);
            }
        }
        assert_eq!(block.nnz(), expect.nnz());
    }
}
