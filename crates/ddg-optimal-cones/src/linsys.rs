//! The linear systems of the Newton iterations, behind one seam.
//!
//! Every Newton step of both stages, for whatever choice of pieces, is a linear system in
//! `(u, φ)` of the same shape. With `W = diag(wE_i m_i)` and `A` the Laplacian block of the
//! unknowns:
//!
//! ```text
//! W u + A φ = 0                                    (adjoint equation, every vertex)
//! (A u)_i                   = µ_i − Ω_i            row i is Flat:    µ_i is a given number
//! (A u)_i − k_i (φ_i − o_i) = b_i − Ω_i            row i is Sloped:  µ_i = b_i + k_i (φ_i − o_i)
//! φ_i                       = p_i                  row i is Pinned:  µ_i is free, φ_i given
//! ```
//!
//! A *flat* row is a vertex whose cone angle is fixed by its piece (no cone, an active bound,
//! a forbidden or prescribed vertex). A *sloped* row is an active vertex of a regularized
//! problem (paper Eq. 22: the entry `(1/γ) D` of Eq. 25). A *pinned* row is an active vertex of
//! the exact problem (paper Eq. 26, where the cone row reads `φ_i = ±λ`).
//!
//! [`NewtonLinearSolver`] is the seam: the Newton code only asks for the solution of this
//! system. [`SaddleSolver`] is the reference implementation, a sparse factorization of the
//! symmetric `2n × 2n` matrix `[[W, A], [A, −K]]` (pinned rows and columns replaced by the
//! identity). A solver that factorizes `A` once and works on the active set alone can be added
//! as a second implementation.
//!
//! The unknown of a sloped row is the excess `ψ_i = φ_i − o_i` over the threshold `o_i = ±t_i`,
//! not `φ_i` itself. The stiffness `k_i` reaches `10¹⁰`; with `φ_i` as unknown the right-hand
//! side would contain `k_i o_i` and every equation would lose that many digits, while
//! `k_i ψ_i = µ_i − b_i` is of order one.

use std::time::{Duration, Instant};

use ddg_core::SparseMatrix;
use ddg_core::sparse::{
    EliminationOrder, Factorization, Method, Parallelism, SparseError, Symbolic,
    fill_reducing_ordering,
};

use crate::problem::Discrete;

/// How one vertex enters a Newton system. See the [module documentation](self).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Row {
    /// The cone angle is the given number.
    Flat {
        /// The cone angle `µ_i`.
        mu: f64,
    },
    /// The cone angle is the affine function `base + k (φ_i − offset)` of the adjoint.
    Sloped {
        /// The slope `k > 0`.
        k: f64,
        /// The adjoint value at which the cone angle equals `base` (the threshold `±t_i`).
        offset: f64,
        /// The cone angle at `φ_i = offset`.
        base: f64,
    },
    /// The adjoint is the given number and the cone angle is free.
    Pinned {
        /// The adjoint value `φ_i`.
        phi: f64,
    },
}

/// Which factorization the reference solver uses for the saddle-point matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LinearSolver {
    /// `LDLᵀ` in the paired elimination order if all distortion weights are positive, `LU`
    /// otherwise; and `LU` for any system on which `LDLᵀ` breaks down or is not accurate
    /// enough after iterative refinement. Once `LU` was needed it is used for the rest of the
    /// solve.
    #[default]
    Auto,
    /// `LDLᵀ` without pivoting, eliminating `u_i` immediately before `φ_i`, the pairs in a
    /// fill-reducing order of the Laplacian. Fast; needs all distortion weights positive.
    SaddleLdlt,
    /// `LU` with partial pivoting. Several times slower, no assumptions.
    SaddleLu,
}

/// What the linear solver did, accumulated over a solve.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LinearSolverReport {
    /// Numeric `LDLᵀ` factorizations of the saddle-point matrix (including failed attempts).
    pub ldlt_factorizations: usize,
    /// Numeric `LU` factorizations of the saddle-point matrix.
    pub lu_factorizations: usize,
    /// Triangular solve pairs with a saddle-point factor (initial solves plus refinement steps).
    pub solves: usize,
    /// Steps of iterative refinement among those solves.
    pub refinement_steps: usize,
    /// Number of systems on which `LDLᵀ` failed or was too inaccurate and `LU` took over.
    pub ldlt_fallbacks: usize,
    /// The largest componentwise backward error `max_i |K x − b|_i / (|K| |x| + |b|)_i` of any
    /// accepted solution: how well each equation holds relative to the size of its terms.
    pub worst_backward_error: f64,
    /// Time spent in symbolic analysis.
    pub symbolic_time: Duration,
    /// Time spent in numeric factorization.
    pub factorization_time: Duration,
    /// Time spent in solves, including refinement.
    pub solve_time: Duration,
}

/// A Newton system could not be solved.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct LinearFailure {
    /// What went wrong, for the report.
    pub message: String,
}

/// The seam between the Newton iterations and the linear algebra.
pub(crate) trait NewtonLinearSolver {
    /// Solves the system described in the [module documentation](self) for the given rows and
    /// returns `(u, φ)`. Pinned entries of `φ` are returned exactly as given.
    fn solve(
        &mut self,
        rows: &[Row],
        par: Parallelism<'_>,
    ) -> Result<(Vec<f64>, Vec<f64>), LinearFailure>;

    /// What has been done so far.
    fn report(&self) -> LinearSolverReport;
}

/// Componentwise backward error up to which a solution is accepted: every equation must hold
/// up to this fraction of the magnitude of its own terms.
const ACCEPT_BACKWARD_ERROR: f64 = 1e-11;
/// Iterative refinement stops once the componentwise backward error is at most this.
const TARGET_BACKWARD_ERROR: f64 = 1e-14;
/// Limit on the refinement steps of one solve.
const MAX_REFINEMENT_STEPS: usize = 4;

/// Reference implementation: sparse factorization of the symmetric saddle-point matrix
/// `[[W, A], [A, −K]]` in the unknowns `(u, φ)`.
///
/// [formula: saddle-point-newton-system]
pub(crate) struct SaddleSolver<'d> {
    d: &'d Discrete,
    strategy: LinearSolver,
    /// The `2n × 2n` matrix; its pattern is fixed (blocks of zeros are stored explicitly).
    k: SparseMatrix,
    /// The paired elimination order for `LDLᵀ`.
    paired_order: Vec<u32>,
    ldlt_symbolic: Option<Symbolic>,
    ldlt: Option<Factorization>,
    lu_symbolic: Option<Symbolic>,
    lu: Option<Factorization>,
    /// Set once `LDLᵀ` has failed in `Auto` mode.
    use_lu: bool,
    report: LinearSolverReport,
}

impl<'d> SaddleSolver<'d> {
    /// Builds the matrix pattern and the elimination order. No factorization happens here.
    pub fn new(d: &'d Discrete, strategy: LinearSolver) -> Result<SaddleSolver<'d>, SparseError> {
        let n = d.n;
        let a = &d.a;
        let start = Instant::now();
        // Column j < n (unknown u_j): W_jj, then A_ij in rows n + i.
        // Column n + j (unknown φ_j): A_ij in rows i, then the diagonal entry.
        let nnz = 2 * a.nnz() + 2 * n;
        let mut col_ptr = Vec::with_capacity(2 * n + 1);
        let mut row_idx = Vec::with_capacity(nnz);
        col_ptr.push(0);
        for j in 0..n {
            row_idx.push(j as u32);
            for p in a.col_ptr()[j]..a.col_ptr()[j + 1] {
                row_idx.push((n + a.row_idx()[p] as usize) as u32);
            }
            col_ptr.push(row_idx.len());
        }
        for j in 0..n {
            for p in a.col_ptr()[j]..a.col_ptr()[j + 1] {
                row_idx.push(a.row_idx()[p]);
            }
            row_idx.push((n + j) as u32);
            col_ptr.push(row_idx.len());
        }
        let k = SparseMatrix::from_csc(2 * n, 2 * n, col_ptr, row_idx, vec![0.0; nnz])?;
        let order = fill_reducing_ordering(a)?;
        let mut paired_order = Vec::with_capacity(2 * n);
        for &i in &order {
            paired_order.push(i);
            paired_order.push(n as u32 + i);
        }
        let report = LinearSolverReport {
            symbolic_time: start.elapsed(),
            ..LinearSolverReport::default()
        };
        Ok(SaddleSolver {
            d,
            strategy,
            k,
            paired_order,
            ldlt_symbolic: None,
            ldlt: None,
            lu_symbolic: None,
            lu: None,
            // With a zero distortion weight the paired order meets a zero pivot (or one that
            // is non-zero only through fill, which is inaccurate): use LU from the start.
            use_lu: strategy == LinearSolver::SaddleLu
                || (strategy == LinearSolver::Auto && d.w.contains(&0.0)),
            report,
        })
    }

    /// Writes the matrix values and the right-hand side for the given rows.
    fn assemble(&mut self, rows: &[Row]) -> Vec<f64> {
        let d = self.d;
        let n = d.n;
        let a = &d.a;
        let a_ptr = a.col_ptr();
        let a_rows = a.row_idx();
        let a_vals = a.values();
        let mut rhs = vec![0.0; 2 * n];
        let pinned = |i: usize| matches!(rows[i], Row::Pinned { .. });
        let values = self.k.values_mut();
        let mut q = 0;
        for j in 0..n {
            values[q] = d.w[j];
            q += 1;
            // Block (φ rows, u_j column): the state rows, which a pinned vertex does not have.
            for p in a_ptr[j]..a_ptr[j + 1] {
                values[q] = if pinned(a_rows[p] as usize) {
                    0.0
                } else {
                    a_vals[p]
                };
                q += 1;
            }
        }
        for j in 0..n {
            let range = a_ptr[j]..a_ptr[j + 1];
            match rows[j] {
                Row::Pinned { phi } => {
                    // The known φ_j moves to the right-hand side of the adjoint equation.
                    for p in range {
                        values[q] = 0.0;
                        q += 1;
                        rhs[a_rows[p] as usize] -= a_vals[p] * phi;
                    }
                    values[q] = 1.0;
                    rhs[n + j] = phi;
                }
                Row::Flat { mu } => {
                    for p in range {
                        values[q] = a_vals[p];
                        q += 1;
                    }
                    values[q] = 0.0;
                    rhs[n + j] = mu - d.omega[j];
                }
                Row::Sloped { k, offset, base } => {
                    // Unknown ψ_j = φ_j − offset: the offset moves to the right-hand side of
                    // the adjoint equation, like a pinned value.
                    for p in range {
                        values[q] = a_vals[p];
                        q += 1;
                        rhs[a_rows[p] as usize] -= a_vals[p] * offset;
                    }
                    values[q] = -k;
                    rhs[n + j] = base - d.omega[j];
                }
            }
            q += 1;
        }
        rhs
    }

    fn factor(&mut self, lu: bool, par: Parallelism<'_>) -> Result<(), SparseError> {
        let (method, symbolic, numeric) = if lu {
            (Method::Lu, &mut self.lu_symbolic, &mut self.lu)
        } else {
            (Method::Ldlt, &mut self.ldlt_symbolic, &mut self.ldlt)
        };
        if symbolic.is_none() {
            let start = Instant::now();
            let order = if lu {
                EliminationOrder::FillReducing
            } else {
                EliminationOrder::Custom(&self.paired_order)
            };
            *symbolic = Some(Symbolic::analyze_with_order(&self.k, method, order)?);
            self.report.symbolic_time += start.elapsed();
        }
        let start = Instant::now();
        if lu {
            self.report.lu_factorizations += 1;
        } else {
            self.report.ldlt_factorizations += 1;
        }
        let result = match numeric {
            Some(f) => f.refactor(&self.k, par),
            None => symbolic
                .as_ref()
                .expect("analysed above")
                .factor(&self.k, par)
                .map(|f| *numeric = Some(f)),
        };
        self.report.factorization_time += start.elapsed();
        result
    }

    /// Factorizes with the chosen method and solves with refinement. `Err` carries a reason.
    fn try_solve(
        &mut self,
        lu: bool,
        rhs: &[f64],
        par: Parallelism<'_>,
    ) -> Result<Vec<f64>, String> {
        let name = if lu { "LU" } else { "LDLt" };
        self.factor(lu, par)
            .map_err(|e| format!("{name} factorization failed: {e}"))?;
        let factorization = if lu { &self.lu } else { &self.ldlt };
        let factorization = factorization.as_ref().expect("factorized above");
        // Solve, then refine until every equation holds to rounding level relative to its own
        // terms. A step that does not improve the solution ends the refinement.
        let start = Instant::now();
        let mut x = rhs.to_vec();
        let outcome = factorization
            .solve_in_place(&mut x, par)
            .map_err(|e| format!("{name} solve failed: {e}"));
        let mut solves = 1;
        let mut residual = vec![0.0; rhs.len()];
        let mut error = f64::NAN;
        if outcome.is_ok() {
            error = backward_error(&self.k, &x, rhs, &mut residual);
            let mut candidate = vec![0.0; rhs.len()];
            let mut candidate_residual = vec![0.0; rhs.len()];
            while error > TARGET_BACKWARD_ERROR && solves <= MAX_REFINEMENT_STEPS {
                let mut correction = residual.clone();
                if factorization.solve_in_place(&mut correction, par).is_err() {
                    break;
                }
                solves += 1;
                for i in 0..x.len() {
                    candidate[i] = x[i] + correction[i];
                }
                let refined = backward_error(&self.k, &candidate, rhs, &mut candidate_residual);
                // A NaN counts as no improvement.
                if refined.is_nan() || refined >= error {
                    break;
                }
                std::mem::swap(&mut x, &mut candidate);
                std::mem::swap(&mut residual, &mut candidate_residual);
                error = refined;
            }
        }
        self.report.solve_time += start.elapsed();
        self.report.solves += solves;
        self.report.refinement_steps += solves - 1;
        outcome?;
        // A NaN is rejected.
        if error.is_nan() || error > ACCEPT_BACKWARD_ERROR {
            return Err(format!(
                "{name} solve is inaccurate: componentwise backward error {error:.3e} after {} \
                 refinement steps",
                solves - 1
            ));
        }
        self.report.worst_backward_error = self.report.worst_backward_error.max(error);
        Ok(x)
    }
}

/// The componentwise backward error `max_i |b − K x|_i / (|K| |x| + |b|)_i` (Oettli–Prager):
/// the smallest relative perturbation of the entries of `K` and `b` for which `x` is an exact
/// solution. Unlike a residual norm it does not depend on how the equations are scaled.
/// `residual` receives `b − K x`. Returns NaN if the residual is not finite.
fn backward_error(k: &SparseMatrix, x: &[f64], b: &[f64], residual: &mut [f64]) -> f64 {
    residual.copy_from_slice(b);
    let mut magnitude: Vec<f64> = b.iter().map(|v| v.abs()).collect();
    let (col_ptr, row_idx, values) = (k.col_ptr(), k.row_idx(), k.values());
    for (j, &xj) in x.iter().enumerate() {
        for p in col_ptr[j]..col_ptr[j + 1] {
            let term = values[p] * xj;
            residual[row_idx[p] as usize] -= term;
            magnitude[row_idx[p] as usize] += term.abs();
        }
    }
    let mut worst = 0.0f64;
    for (r, m) in residual.iter().zip(&magnitude) {
        if !r.is_finite() {
            return f64::NAN;
        }
        if *m > 0.0 {
            worst = worst.max(r.abs() / m);
        }
    }
    worst
}

impl NewtonLinearSolver for SaddleSolver<'_> {
    fn solve(
        &mut self,
        rows: &[Row],
        par: Parallelism<'_>,
    ) -> Result<(Vec<f64>, Vec<f64>), LinearFailure> {
        let n = self.d.n;
        debug_assert_eq!(rows.len(), n);
        let rhs = self.assemble(rows);
        let mut x = if self.use_lu {
            self.try_solve(true, &rhs, par)
                .map_err(|message| LinearFailure { message })?
        } else {
            match self.try_solve(false, &rhs, par) {
                Ok(x) => x,
                Err(ldlt_message) if self.strategy == LinearSolver::Auto => {
                    self.report.ldlt_fallbacks += 1;
                    self.use_lu = true;
                    self.try_solve(true, &rhs, par)
                        .map_err(|lu_message| LinearFailure {
                            message: format!("{ldlt_message}; {lu_message}"),
                        })?
                }
                Err(message) => return Err(LinearFailure { message }),
            }
        };
        let mut phi = x.split_off(n);
        for (value, row) in phi.iter_mut().zip(rows) {
            match *row {
                Row::Pinned { phi } => *value = phi,
                Row::Sloped { offset, .. } => *value += offset,
                Row::Flat { .. } => {}
            }
        }
        Ok((x, phi))
    }

    fn report(&self) -> LinearSolverReport {
        self.report
    }
}
