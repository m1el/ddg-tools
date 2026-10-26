//! Sparse saddle-point solver for the classic metrics (`SPEC.md` §4.6, §6.4, §8.5).
//!
//! `[M ⊗ I_d  Cᵀ; C  0]` assembled as a sparse matrix and factored once by sparse LU with
//! partial pivoting (`ddg_core::sparse`, `Method::Lu`), sequentially (deterministic). Used for
//! the L², H¹, H² projected gradients of the Sec. 7 evaluation.

use std::time::Instant;

use ddg_core::sparse::{Factorization, Method, Parallelism, SparseMatrix, Symbolic};

use super::{LinearReport, SaddleSolution, SaddleSolver, check_rhs, sparse_error};
use crate::error::RcError;
use crate::metric::DofLayout;

/// The factored sparse saddle matrix of one step.
#[derive(Debug)]
pub struct SparseSaddle {
    k: usize,
    nd: usize,
    matrix: SparseMatrix,
    factor: Factorization,
    factor_seconds: f64,
}

impl SparseSaddle {
    /// Builds and factors `[M ⊗ I_d  Cᵀ; C 0]` for a sparse symmetric metric `m`
    /// (`|V| × |V|`) and the Jacobian `c` (`k × d|V|`).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for mismatched dimensions; [`RcError::SingularSaddle`] if the
    /// LU finds no pivot.
    pub fn new(m: &SparseMatrix, c: &SparseMatrix, layout: DofLayout) -> Result<Self, RcError> {
        let t0 = Instant::now();
        let d = layout.dim();
        let n = m.nrows();
        let nd = n * d;
        if m.ncols() != n || c.ncols() != nd {
            return Err(RcError::InvalidOption {
                name: "sparse saddle",
                value: format!(
                    "metric {}x{}, Jacobian {}x{}",
                    n,
                    m.ncols(),
                    c.nrows(),
                    c.ncols()
                ),
                requirement: "a square metric and a Jacobian with d|V| columns",
            });
        }
        let k = c.nrows();
        let mut t = Vec::with_capacity(m.nnz() * d + 2 * c.nnz());
        for (r, col, v) in m.triplets() {
            for s in 0..d as u32 {
                t.push((r * d as u32 + s, col * d as u32 + s, v));
            }
        }
        for (r, col, v) in c.triplets() {
            t.push((nd as u32 + r, col, v));
            t.push((col, nd as u32 + r, v));
        }
        let matrix = SparseMatrix::from_triplets(nd + k, nd + k, &t)
            .map_err(|e| sparse_error(e, "sparse saddle"))?;
        let sym =
            Symbolic::analyze(&matrix, Method::Lu).map_err(|e| sparse_error(e, "sparse saddle"))?;
        let factor = sym
            .factor(&matrix, Parallelism::Sequential)
            .map_err(|e| sparse_error(e, "sparse saddle"))?;
        Ok(SparseSaddle {
            k,
            nd,
            matrix,
            factor,
            factor_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// Seconds spent assembling and factoring.
    pub fn factor_seconds(&self) -> f64 {
        self.factor_seconds
    }
}

impl SaddleSolver for SparseSaddle {
    fn n_dofs(&self) -> usize {
        self.nd
    }

    fn n_constraints(&self) -> usize {
        self.k
    }

    /// [formula: rc-saddle-gradient]
    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError> {
        check_rhs(self, a, b)?;
        let t0 = Instant::now();
        let rhs: Vec<f64> = a.iter().chain(b).copied().collect();
        let mut sol = rhs.clone();
        self.factor
            .solve_in_place(&mut sol, Parallelism::Sequential)
            .map_err(|e| sparse_error(e, "sparse saddle"))?;
        let mut r = vec![0.0; sol.len()];
        self.matrix
            .mul_vec(&sol, &mut r)
            .map_err(|e| sparse_error(e, "sparse saddle"))?;
        let r2: f64 = r.iter().zip(&rhs).map(|(p, q)| (p - q) * (p - q)).sum();
        let rhs2: f64 = rhs.iter().map(|v| v * v).sum();
        let relative_residual = if rhs2 > 0.0 {
            (r2 / rhs2).sqrt()
        } else {
            r2.sqrt()
        };
        let lambda = sol.split_off(self.nd);
        Ok(SaddleSolution {
            x: sol,
            lambda,
            report: LinearReport {
                solver: "sparse LU",
                dimension: self.nd + self.k,
                iterations: 0,
                relative_residual,
                fallback_used: false,
                levels: 0,
                pivot_ratio: None,
                seconds: t0.elapsed().as_secs_f64(),
            },
        })
    }
}
