//! Dense saddle-point solver (`SPEC.md` §6.1–6.2, §8.5).
//!
//! One LU factorization with partial pivoting (`faer`, on an explicit thread pool) of
//!
//! ```text
//! K = [ Ā/s  Cᵀ ]      Ā = A ⊗ I_d (interleaved, D29),  s = max_v A_vv,
//!     [ C    0  ]
//! ```
//!
//! reused for Eq. 22 and every chord iteration of the projection (Sec. 7.4). The scaling by
//! `s` only balances the two blocks for pivoting (`λ` is scaled back). The factorization is
//! `O((d|V| + k)³)`.

use std::sync::Arc;
use std::time::Instant;

use ddg_core::sparse::SparseMatrix;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::lu::partial_pivoting::{factor as lu_factor, solve as lu_solve};
use faer::perm::PermRef;
use faer::{Mat, Par};
use rayon::prelude::*;

use super::{LinearReport, SaddleSolution, SaddleSolver, check_rhs};
use crate::error::RcError;
use crate::metric::{DenseMatrix, DofLayout};

/// Pivots with `|u_ii| ≤ PIVOT_TOLERANCE · max |u_jj|` mean the saddle matrix is singular.
pub const PIVOT_TOLERANCE: f64 = 1e-14;

/// The factored dense saddle matrix of one step.
#[derive(Debug)]
pub struct DenseSaddle {
    pool: Arc<rayon::ThreadPool>,
    layout: DofLayout,
    n_vertices: usize,
    k: usize,
    scale: f64,
    lu: Mat<f64>,
    perm: Vec<usize>,
    perm_inv: Vec<usize>,
    /// Copies for the residual check.
    a: DenseMatrix,
    c: SparseMatrix,
    ct: SparseMatrix,
    pivot_ratio: f64,
    factor_seconds: f64,
}

impl DenseSaddle {
    /// Builds and factors the saddle matrix of the metric `a` (`|V| × |V|`, symmetric) and the
    /// constraint Jacobian `c` (`k × d|V|`) on `pool`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] if the dimensions do not match; [`RcError::SingularMetric`]
    /// (`components: 0`) if `a` has no positive diagonal entry; [`RcError::SingularSaddle`] with
    /// the index of the smallest pivot if a pivot is below [`PIVOT_TOLERANCE`] relative to the
    /// largest (or not finite).
    pub fn new(
        pool: Arc<rayon::ThreadPool>,
        a: &DenseMatrix,
        c: &SparseMatrix,
        layout: DofLayout,
    ) -> Result<DenseSaddle, RcError> {
        let t0 = Instant::now();
        let n = a.n();
        let d = layout.dim();
        let nd = n * d;
        if c.ncols() != nd {
            return Err(RcError::InvalidOption {
                name: "constraint Jacobian",
                value: format!("{} columns", c.ncols()),
                requirement: "d|V| columns, matching the metric",
            });
        }
        let k = c.nrows();
        let dim = nd + k;
        let scale = (0..n).map(|v| a.get(v, v)).fold(0.0f64, f64::max);
        if !(scale > 0.0 && scale.is_finite()) {
            return Err(RcError::SingularMetric { components: 0 });
        }
        let inv = 1.0 / scale;
        let ct = c.transpose();
        let mut lu = Mat::<f64>::zeros(dim, dim);
        let n_threads = pool.current_num_threads();
        let par = if n_threads > 1 {
            Par::rayon(n_threads)
        } else {
            Par::Seq
        };
        let (ccp, cri, cval) = (c.col_ptr(), c.row_idx(), c.values());
        let (tcp, tri, tval) = (ct.col_ptr(), ct.row_idx(), ct.values());
        pool.install(|| {
            lu.par_col_iter_mut().enumerate().for_each(|(j, mut col)| {
                if j < nd {
                    let (v, s) = (j / d, j % d);
                    let row = a.row(v); // A is symmetric: column v = row v
                    for w in 0..n {
                        col[w * d + s] = row[w] * inv;
                    }
                    for p in ccp[j]..ccp[j + 1] {
                        col[nd + cri[p] as usize] = cval[p];
                    }
                } else {
                    let r = j - nd;
                    for p in tcp[r]..tcp[r + 1] {
                        col[tri[p] as usize] = tval[p];
                    }
                }
            });
        });
        let mut perm = vec![0usize; dim];
        let mut perm_inv = vec![0usize; dim];
        pool.install(|| {
            let mut mem = MemBuffer::new(lu_factor::lu_in_place_scratch::<usize, f64>(
                dim,
                dim,
                par,
                Default::default(),
            ));
            lu_factor::lu_in_place(
                lu.as_mut(),
                &mut perm,
                &mut perm_inv,
                par,
                MemStack::new(&mut mem),
                Default::default(),
            );
        });
        let mut max_p: f64 = 0.0;
        let mut min_p = f64::INFINITY;
        let mut argmin = 0;
        for i in 0..dim {
            let p = lu[(i, i)].abs();
            if !p.is_finite() {
                return Err(RcError::SingularSaddle { pivot: i });
            }
            max_p = max_p.max(p);
            if p < min_p {
                min_p = p;
                argmin = i;
            }
        }
        let pivot_ratio = if dim == 0 { 1.0 } else { min_p / max_p };
        if pivot_ratio.is_nan() || pivot_ratio <= PIVOT_TOLERANCE {
            return Err(RcError::SingularSaddle { pivot: argmin });
        }
        Ok(DenseSaddle {
            pool,
            layout,
            n_vertices: n,
            k,
            scale,
            lu,
            perm,
            perm_inv,
            a: a.clone(),
            c: c.clone(),
            ct,
            pivot_ratio,
            factor_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// Seconds spent assembling and factoring.
    pub fn factor_seconds(&self) -> f64 {
        self.factor_seconds
    }

    /// `min |u_ii| / max |u_ii|` of the LU factor.
    pub fn pivot_ratio(&self) -> f64 {
        self.pivot_ratio
    }

    /// The metric `A` the solver was built with.
    pub fn metric(&self) -> &DenseMatrix {
        &self.a
    }

    /// The constraint Jacobian `C` the solver was built with.
    pub fn jacobian(&self) -> &SparseMatrix {
        &self.c
    }

    /// The layout of the unknowns.
    pub fn layout(&self) -> DofLayout {
        self.layout
    }

    /// `Ā x` (`A ⊗ I_d` times an interleaved vector), parallel over vertices, fixed order.
    pub fn apply_metric(&self, x: &[f64]) -> Vec<f64> {
        let d = self.layout.dim();
        let n = self.n_vertices;
        let mut y = vec![0.0; n * d];
        self.pool.install(|| {
            y.par_chunks_mut(d).enumerate().for_each(|(v, yv)| {
                let row = self.a.row(v);
                for (w, &aw) in row.iter().enumerate() {
                    for s in 0..d {
                        yv[s] += aw * x[w * d + s];
                    }
                }
            });
        });
        y
    }
}

impl SaddleSolver for DenseSaddle {
    fn n_dofs(&self) -> usize {
        self.n_vertices * self.layout.dim()
    }

    fn n_constraints(&self) -> usize {
        self.k
    }

    /// [formula: rc-saddle-gradient]
    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError> {
        check_rhs(self, a, b)?;
        let t0 = Instant::now();
        let nd = self.n_dofs();
        let dim = nd + self.k;
        let inv = 1.0 / self.scale;
        let mut rhs =
            Mat::<f64>::from_fn(dim, 1, |i, _| if i < nd { a[i] * inv } else { b[i - nd] });
        // The triangular solves of one right-hand side run sequentially: they are memory bound
        // (no gain from threads) and faer's parallel triangular solve sums in an order that
        // depends on the thread count, which would break D25.
        let par = Par::Seq;
        let perm = PermRef::new_checked(&self.perm, &self.perm_inv, dim);
        let mut mem = MemBuffer::new(lu_solve::solve_in_place_scratch::<usize, f64>(dim, 1, par));
        lu_solve::solve_in_place(
            self.lu.as_ref(),
            self.lu.as_ref(),
            perm,
            rhs.as_mut(),
            par,
            MemStack::new(&mut mem),
        );
        let x: Vec<f64> = (0..nd).map(|i| rhs[(i, 0)]).collect();
        let lambda: Vec<f64> = (0..self.k).map(|r| rhs[(nd + r, 0)] * self.scale).collect();
        if let Some(i) = x.iter().chain(&lambda).position(|v| !v.is_finite()) {
            return Err(RcError::SingularSaddle { pivot: i });
        }
        // Residual of the unscaled system.
        let ax = self.apply_metric(&x);
        let mut ctl = vec![0.0; nd];
        self.ct
            .mul_vec(&lambda, &mut ctl)
            .map_err(|e| super::sparse_error(e, "C^T"))?;
        let mut cx = vec![0.0; self.k];
        self.c
            .mul_vec(&x, &mut cx)
            .map_err(|e| super::sparse_error(e, "C"))?;
        let mut r2 = 0.0;
        for i in 0..nd {
            let r = ax[i] + ctl[i] - a[i];
            r2 += r * r;
        }
        for i in 0..self.k {
            let r = cx[i] - b[i];
            r2 += r * r;
        }
        let rhs2: f64 = a.iter().chain(b).map(|v| v * v).sum();
        let relative_residual = if rhs2 > 0.0 {
            (r2 / rhs2).sqrt()
        } else {
            r2.sqrt()
        };
        Ok(SaddleSolution {
            x,
            lambda,
            report: LinearReport {
                solver: "dense LU",
                dimension: dim,
                iterations: 0,
                relative_residual,
                fallback_used: false,
                levels: 0,
                pivot_ratio: Some(self.pivot_ratio),
                seconds: t0.elapsed().as_secs_f64(),
            },
        })
    }
}
