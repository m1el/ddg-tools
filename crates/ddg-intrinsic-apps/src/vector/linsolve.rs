//! Sparse solves shared by the modules of this crate.

use ddg_core::sparse::{Factorization, Method, Parallelism, Symbolic};
use ddg_core::{SparseError, SparseMatrix};

/// Which factorization a solve used.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Factorized {
    /// Cholesky: the matrix was positive definite.
    Cholesky,
    /// LU with partial pivoting, after Cholesky reported a matrix that is not positive
    /// definite (a connection Laplacian with negative cotan weights can be indefinite).
    Lu,
}

/// Cholesky of a symmetric matrix, or LU if Cholesky finds it not positive definite.
pub(crate) fn factor(a: &SparseMatrix) -> Result<(Factorization, Factorized), SparseError> {
    let symbolic = Symbolic::analyze(a, Method::Llt)?;
    match symbolic.factor(a, Parallelism::Sequential) {
        Ok(f) => Ok((f, Factorized::Cholesky)),
        Err(SparseError::NotPositiveDefinite { .. }) => {
            let symbolic = Symbolic::analyze(a, Method::Lu)?;
            Ok((symbolic.factor(a, Parallelism::Sequential)?, Factorized::Lu))
        }
        Err(other) => Err(other),
    }
}

/// Cholesky only: an error if the matrix is not positive definite.
pub(crate) fn factor_spd(a: &SparseMatrix) -> Result<Factorization, SparseError> {
    Symbolic::analyze(a, Method::Llt)?.factor(a, Parallelism::Sequential)
}

/// Solves in place.
pub(crate) fn solve(f: &Factorization, rhs: &mut [f64]) -> Result<(), SparseError> {
    f.solve_in_place(rhs, Parallelism::Sequential)
}

/// `y = A x`.
pub(crate) fn mul(a: &SparseMatrix, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; a.nrows()];
    a.mul_vec(x, &mut y)
        .expect("dimensions are consistent by construction");
    y
}

pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

pub(crate) fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
