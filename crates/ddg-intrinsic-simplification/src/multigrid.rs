//! Galerkin multigrid solver for symmetric positive definite systems on a hierarchy of
//! prolongation matrices (SPEC §9.2).
//!
//! Hsueh-Ti Derek Liu, Jiayi Eris Zhang, Mirela Ben-Chen, Alec Jacobson, *Surface Multigrid via
//! Intrinsic Prolongation*, ACM TOG 40(4), 2021 (`papers/liu2021_surface_multigrid.pdf`,
//! Alg. 1 and 2, App. D), and the `min_quad_with_fixed_mg` of the authors' reference code of the
//! intrinsic-simplification paper (MIT, Copyright (c) 2023 Hsueh-Ti Derek Liu;
//! `src/mg_precompute.cpp`, `mg_VCycle.cpp`, `min_quad_with_fixed_mg.cpp`; SPEC §9.2, map §10).
//! This is a port of the algorithm to the workspace's sparse types, not of the code.
//!
//! # The method
//!
//! Given the fine matrix `A` (`n × n`, symmetric positive definite once the fixed rows are
//! removed) and prolongations `P_k ∈ ℝ^{n_k × n_{k+1}}` (`n_0 = n`), the Galerkin hierarchy is
//! `A_0 = A_uu`, `A_{k+1} = P_kᵀ A_k P_k` ([`Multigrid::new`]). Every level is symmetric positive
//! definite when `A_0` is and every `P_k` has full column rank (`iem-galerkin-spd`, checked on
//! the matrices by the tests; on a rank-deficient `P` the setup fails at the coarsest Cholesky
//! factorization or at a non-positive diagonal, and says so).
//!
//! * **Dirichlet elimination** (as `min_quad_with_fixed_mg_precompute`): the fixed vertices
//!   (`known`) are removed from the system (`A_uu`, `b_u − A_uk x_k`) and their rows from `P_0`;
//!   then at every level the columns of `P_k` without any entry `> 1e-15` (coarse vertices
//!   supported only on removed rows) are dropped, which removes the same coarse vertices from
//!   the rows of `P_{k+1}`, and so on up the hierarchy. (The reference stops the cascade at the
//!   first level where nothing is dropped; this one always runs through all levels.)
//! * **Smoothing**: Gauss–Seidel, `ν₁ = ν₂ = 2` by default. The pre-smoother is `ν₁` *forward*
//!   sweeps and the post-smoother `ν₂` *backward* sweeps: the post-smoother is the transpose of
//!   the pre-smoother, so the V-cycle operator is symmetric positive definite and can be used as
//!   a preconditioner of conjugate gradients ([`Multigrid::pcg`]). (The reference uses forward
//!   sweeps twice, which is not symmetric.)
//! * **Coarsest level**: sparse Cholesky (the workspace's) of `A_L`. If the matrix is not
//!   positive definite (a pure Neumann coarse problem, say), a diagonal shift
//!   `coarsest_shift · max_i (A_L)_ii` is added and the factorization retried once; the shift
//!   used is reported ([`MultigridInfo::coarsest_shift`], `0` when none was needed). The
//!   reference adds `1e-12` unconditionally.
//! * **Stopping**: V-cycles until `‖b_u − A_0 x_u‖₂ ≤ tolerance · ‖b_u − A_uk x_k‖₂` (relative to
//!   the reduced right-hand side; to the initial residual if that is zero) or `max_cycles`.
//!   The reference's criterion is absolute with at most 20 cycles.
//!
//! # Reports
//!
//! [`MultigridInfo`] (level sizes, dropped columns, shift, setup time) and [`SolveReport`]
//! (residual after each cycle, convergence factor, solve time).

use std::time::Instant;

use ddg_core::sparse::{Factorization, Method, Parallelism, Symbolic};
use ddg_core::{SparseError, SparseMatrix, TriMesh};

/// Errors of the multigrid setup and solves.
#[derive(Debug, Clone, PartialEq)]
pub enum MultigridError {
    /// The fine matrix is not square.
    NotSquare {
        /// Rows.
        rows: usize,
        /// Columns.
        cols: usize,
    },
    /// The fine matrix is not symmetric.
    NotSymmetric {
        /// Largest `|a_ij − a_ji|`.
        defect: f64,
    },
    /// A size does not match.
    SizeMismatch {
        /// What is wrong.
        what: String,
        /// Expected size.
        expected: usize,
        /// Size found.
        got: usize,
    },
    /// A fixed vertex is out of range or listed twice, or every vertex is fixed.
    BadKnown {
        /// What is wrong.
        reason: String,
    },
    /// A level has no unknowns (zero columns in its prolongation, possibly after the dropping
    /// of the columns that only the fixed vertices support).
    EmptyLevel {
        /// The prolongation `P_level` (from level `level + 1` to `level`).
        level: usize,
        /// Why.
        reason: String,
    },
    /// A diagonal entry is not positive (the matrix is not positive definite) or not finite.
    BadDiagonal {
        /// The level.
        level: usize,
        /// The row.
        row: usize,
        /// The value.
        value: f64,
    },
    /// An input contains a NaN or infinity.
    NonFinite {
        /// Which input.
        what: String,
    },
    /// An option is out of range.
    BadOption {
        /// What is wrong.
        reason: String,
    },
    /// The coarsest matrix could not be factorized, even with a shift.
    Coarsest {
        /// The shift of the last attempt.
        shift: f64,
        /// The failure.
        reason: String,
    },
    /// A failure of the sparse backend.
    Sparse(SparseError),
}

impl std::fmt::Display for MultigridError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MultigridError::NotSquare { rows, cols } => {
                write!(f, "the system matrix is {rows} x {cols}, not square")
            }
            MultigridError::NotSymmetric { defect } => {
                write!(f, "the system matrix is not symmetric (max |a_ij - a_ji| = {defect:e})")
            }
            MultigridError::SizeMismatch {
                what,
                expected,
                got,
            } => write!(f, "{what}: expected {expected}, got {got}"),
            MultigridError::BadKnown { reason } => write!(f, "bad fixed vertices: {reason}"),
            MultigridError::EmptyLevel { level, reason } => {
                write!(f, "the level below prolongation P_{level} has no unknowns: {reason}")
            }
            MultigridError::BadDiagonal { level, row, value } => write!(
                f,
                "level {level}: diagonal entry {row} is {value:e}; the matrix is not positive definite"
            ),
            MultigridError::NonFinite { what } => write!(f, "{what} contains a non-finite value"),
            MultigridError::BadOption { reason } => write!(f, "bad multigrid option: {reason}"),
            MultigridError::Coarsest { shift, reason } => write!(
                f,
                "the coarsest matrix cannot be factorized (diagonal shift {shift:e}): {reason}"
            ),
            MultigridError::Sparse(e) => write!(f, "sparse backend: {e}"),
        }
    }
}

impl std::error::Error for MultigridError {}

impl From<SparseError> for MultigridError {
    fn from(e: SparseError) -> Self {
        MultigridError::Sparse(e)
    }
}

/// Options of [`Multigrid::new`] and of the solves.
#[derive(Debug, Clone, PartialEq)]
pub struct MultigridOptions {
    /// Forward Gauss–Seidel sweeps before the coarse-grid correction (`ν₁`, default 2).
    pub pre_sweeps: usize,
    /// Backward Gauss–Seidel sweeps after it (`ν₂`, default 2).
    pub post_sweeps: usize,
    /// Relative residual tolerance of [`Multigrid::solve`] (default `1e-10`).
    pub tolerance: f64,
    /// Cap on the V-cycles of [`Multigrid::solve`] (default 50).
    pub max_cycles: usize,
    /// Relative diagonal shift tried if the coarsest matrix is not positive definite
    /// (default `1e-12`; `0` disables the fallback: the setup then fails).
    pub coarsest_shift: f64,
}

impl Default for MultigridOptions {
    fn default() -> Self {
        MultigridOptions {
            pre_sweeps: 2,
            post_sweeps: 2,
            tolerance: 1e-10,
            max_cycles: 50,
            coarsest_shift: 1e-12,
        }
    }
}

/// What the setup built.
#[derive(Debug, Clone, PartialEq)]
pub struct MultigridInfo {
    /// Unknowns per level, finest first (after the elimination of the fixed vertices).
    pub level_sizes: Vec<usize>,
    /// Stored entries of the level matrices.
    pub level_nnz: Vec<usize>,
    /// Columns dropped from each prolongation because only fixed vertices supported them.
    pub dropped_columns: Vec<usize>,
    /// The relative diagonal shift added to the coarsest matrix (`0.0` if none).
    pub coarsest_shift: f64,
    /// Number of fixed vertices.
    pub n_known: usize,
    /// Seconds spent in the setup (Galerkin products, factorization).
    pub setup_seconds: f64,
}

/// What a solve did.
#[derive(Debug, Clone, PartialEq)]
pub struct SolveReport {
    /// Whether the relative residual reached the tolerance.
    pub converged: bool,
    /// V-cycles (or CG iterations) performed.
    pub cycles: usize,
    /// `‖b_u − A_0 x_u‖₂` before the first cycle and after each cycle (`cycles + 1` values).
    pub residuals: Vec<f64>,
    /// `‖b_u − A_uk x_k‖₂`, the scale of the relative residual (the initial residual if zero).
    pub rhs_norm: f64,
    /// Last residual over `rhs_norm`.
    pub relative_residual: f64,
    /// Geometric mean of the residual reduction per cycle,
    /// `(r_N / r_0)^{1/N}` (`0` if no cycle was needed).
    pub convergence_factor: f64,
    /// Unknowns per level.
    pub level_sizes: Vec<usize>,
    /// Seconds spent in the solve.
    pub solve_seconds: f64,
}

/// The Galerkin hierarchy of a system with optional fixed values, ready to solve.
#[derive(Debug)]
pub struct Multigrid {
    n: usize,
    options: MultigridOptions,
    unknown: Vec<u32>,
    known: Vec<u32>,
    a_uk: SparseMatrix,
    a: Vec<SparseMatrix>,
    p: Vec<SparseMatrix>,
    pt: Vec<SparseMatrix>,
    diag: Vec<Vec<f64>>,
    coarse: Factorization,
    info: MultigridInfo,
}

/// `C = A B` for compressed-column matrices.
fn mat_mul(a: &SparseMatrix, b: &SparseMatrix) -> Result<SparseMatrix, SparseError> {
    let (m, n) = (a.nrows(), b.ncols());
    let (acp, ari, av) = (a.col_ptr(), a.row_idx(), a.values());
    let (bcp, bri, bv) = (b.col_ptr(), b.row_idx(), b.values());
    let mut col_ptr = Vec::with_capacity(n + 1);
    col_ptr.push(0usize);
    let mut row_idx: Vec<u32> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    let mut acc = vec![0.0f64; m];
    let mut mark = vec![usize::MAX; m];
    let mut touched: Vec<u32> = Vec::new();
    for j in 0..n {
        touched.clear();
        for kb in bcp[j]..bcp[j + 1] {
            let (r, s) = (bri[kb] as usize, bv[kb]);
            for ka in acp[r]..acp[r + 1] {
                let i = ari[ka] as usize;
                if mark[i] != j {
                    mark[i] = j;
                    acc[i] = 0.0;
                    touched.push(i as u32);
                }
                acc[i] += av[ka] * s;
            }
        }
        touched.sort_unstable();
        for &i in &touched {
            row_idx.push(i);
            values.push(acc[i as usize]);
        }
        col_ptr.push(row_idx.len());
    }
    SparseMatrix::from_csc(m, n, col_ptr, row_idx, values)
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|t| t * t).sum::<f64>().sqrt()
}

/// The prolongation of one level of regular midpoint subdivision (4-to-1 split of every
/// triangle) to its parent: `(n + n_edges) × n`, identity on the old vertices, and the new
/// vertex `n + e` on edge `e = [a, b]` of `mesh` is the average of `a` and `b` (linear
/// interpolation, exact). The numbering is that of `ddg_intrinsic_apps::fem::regular_subdivision`.
///
/// # Errors
/// Only if the sizes exceed the sparse index range.
pub fn midpoint_subdivision_prolongation(mesh: &TriMesh) -> Result<SparseMatrix, SparseError> {
    let n = mesh.n_vertices();
    let mut t: Vec<(u32, u32, f64)> = (0..n as u32).map(|i| (i, i, 1.0)).collect();
    for (e, &[a, b]) in mesh.edges().iter().enumerate() {
        t.push(((n + e) as u32, a, 0.5));
        t.push(((n + e) as u32, b, 0.5));
    }
    SparseMatrix::from_triplets(n + mesh.n_edges(), n, &t)
}

impl Multigrid {
    /// Builds the hierarchy for `A x = b` with the vertices `known` fixed.
    ///
    /// `prolongations[k]` is `P_k` (`n_k × n_{k+1}`, with `n_0 = a.nrows()`); an empty slice is a
    /// single level (a direct solve). See the module documentation for the elimination, the
    /// smoothers and the coarsest solve.
    ///
    /// [formula: iem-galerkin-spd]
    ///
    /// # Errors
    /// [`MultigridError`]: `A` not square or not symmetric, a prolongation whose rows do not
    /// match the level above, a level with zero columns, `known` out of range or repeated or
    /// covering every vertex, non-finite entries, a non-positive diagonal at some level, a
    /// coarsest matrix that cannot be factorized.
    pub fn new(
        a: &SparseMatrix,
        prolongations: &[SparseMatrix],
        known: &[u32],
        options: MultigridOptions,
    ) -> Result<Multigrid, MultigridError> {
        let t0 = Instant::now();
        let n = a.nrows();
        if a.ncols() != n {
            return Err(MultigridError::NotSquare {
                rows: n,
                cols: a.ncols(),
            });
        }
        if !(options.tolerance >= 0.0 && options.tolerance.is_finite()) {
            return Err(MultigridError::BadOption {
                reason: format!("tolerance must be finite and non-negative, got {}", options.tolerance),
            });
        }
        if !(options.coarsest_shift >= 0.0 && options.coarsest_shift.is_finite()) {
            return Err(MultigridError::BadOption {
                reason: format!(
                    "coarsest_shift must be finite and non-negative, got {}",
                    options.coarsest_shift
                ),
            });
        }
        if a.values().iter().any(|v| !v.is_finite()) {
            return Err(MultigridError::NonFinite {
                what: "the system matrix".into(),
            });
        }
        let max_abs = a.values().iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let defect = a.symmetry_defect()?;
        if defect > 1e-8 * max_abs.max(f64::MIN_POSITIVE) {
            return Err(MultigridError::NotSymmetric { defect });
        }
        // Prolongation shapes.
        let mut rows = n;
        for (k, p) in prolongations.iter().enumerate() {
            if p.nrows() != rows {
                return Err(MultigridError::SizeMismatch {
                    what: format!("rows of P_{k} (must equal the unknowns of level {k})"),
                    expected: rows,
                    got: p.nrows(),
                });
            }
            if p.ncols() == 0 {
                return Err(MultigridError::EmptyLevel {
                    level: k,
                    reason: "P has zero columns".into(),
                });
            }
            if p.values().iter().any(|v| !v.is_finite()) {
                return Err(MultigridError::NonFinite {
                    what: format!("P_{k}"),
                });
            }
            rows = p.ncols();
        }
        // Fixed vertices.
        let mut is_known = vec![false; n];
        for &k in known {
            let k = k as usize;
            if k >= n {
                return Err(MultigridError::BadKnown {
                    reason: format!("vertex {k} is out of range (the system has {n} rows)"),
                });
            }
            if is_known[k] {
                return Err(MultigridError::BadKnown {
                    reason: format!("vertex {k} is listed twice"),
                });
            }
            is_known[k] = true;
        }
        let unknown: Vec<u32> = (0..n as u32).filter(|&i| !is_known[i as usize]).collect();
        if unknown.is_empty() {
            return Err(MultigridError::BadKnown {
                reason: "every vertex is fixed; there is nothing to solve".into(),
            });
        }
        // `known` may be unordered; the solves take its values in the order given.
        let a0 = a.submatrix(&unknown, &unknown)?;
        let a_uk = a.submatrix(&unknown, known)?;

        // Prolongations with the fixed rows removed and the unsupported columns dropped.
        let mut p: Vec<SparseMatrix> = Vec::with_capacity(prolongations.len());
        let mut dropped = Vec::with_capacity(prolongations.len());
        let mut rows_keep: Option<Vec<u32>> = Some(unknown.clone());
        for (k, pk) in prolongations.iter().enumerate() {
            let all_cols: Vec<u32> = (0..pk.ncols() as u32).collect();
            let rk = rows_keep.take().expect("rows set at every iteration");
            let pr = if rk.len() == pk.nrows() {
                pk.clone()
            } else {
                pk.submatrix(&rk, &all_cols)?
            };
            // Columns with an entry > 1e-15.
            let (cp, vals) = (pr.col_ptr(), pr.values());
            let keep: Vec<u32> = (0..pr.ncols())
                .filter(|&j| vals[cp[j]..cp[j + 1]].iter().any(|v| v.abs() > 1e-15))
                .map(|j| j as u32)
                .collect();
            if keep.is_empty() {
                return Err(MultigridError::EmptyLevel {
                    level: k,
                    reason: "no coarse vertex is supported by a free fine vertex".into(),
                });
            }
            dropped.push(pr.ncols() - keep.len());
            let pr = if keep.len() == pr.ncols() {
                pr
            } else {
                let all_rows: Vec<u32> = (0..pr.nrows() as u32).collect();
                pr.submatrix(&all_rows, &keep)?
            };
            p.push(pr);
            rows_keep = Some(keep);
        }
        // The kept coarse vertices of level k+1 select the rows of P_{k+1} (done above through
        // `rows_keep` at the next iteration); the matrices of different levels are consistent
        // because `keep` is the row set for the next prolongation.

        // Galerkin hierarchy.
        let pt: Vec<SparseMatrix> = p.iter().map(|m| m.transpose()).collect();
        let mut mats = vec![a0];
        for k in 0..p.len() {
            let ap = mat_mul(&mats[k], &p[k])?;
            mats.push(mat_mul(&pt[k], &ap)?);
        }
        let mut diag = Vec::with_capacity(mats.len());
        for (lv, m) in mats.iter().enumerate() {
            let d = m.diagonal();
            if let Some((row, &value)) = d.iter().enumerate().find(|(_, v)| !(**v > 0.0 && v.is_finite())) {
                return Err(MultigridError::BadDiagonal {
                    level: lv,
                    row,
                    value,
                });
            }
            diag.push(d);
        }
        // Coarsest factorization.
        let last = mats.last().expect("at least one level");
        let mut shift = 0.0;
        let coarse = {
            let first = Symbolic::analyze(last, Method::Llt)?.factor(last, Parallelism::Sequential);
            match first {
                Ok(f) => f,
                Err(SparseError::NotPositiveDefinite { .. }) if options.coarsest_shift > 0.0 => {
                    let dmax = diag.last().expect("level").iter().cloned().fold(0.0, f64::max);
                    shift = options.coarsest_shift;
                    let mut shifted = last.clone();
                    shifted.add_to_diagonal(&vec![shift * dmax; last.nrows()])?;
                    Symbolic::analyze(&shifted, Method::Llt)?
                        .factor(&shifted, Parallelism::Sequential)
                        .map_err(|e| MultigridError::Coarsest {
                            shift,
                            reason: e.to_string(),
                        })?
                }
                Err(e) => {
                    return Err(MultigridError::Coarsest {
                        shift: 0.0,
                        reason: e.to_string(),
                    });
                }
            }
        };
        let info = MultigridInfo {
            level_sizes: mats.iter().map(|m| m.nrows()).collect(),
            level_nnz: mats.iter().map(|m| m.nnz()).collect(),
            dropped_columns: dropped,
            coarsest_shift: shift,
            n_known: known.len(),
            setup_seconds: t0.elapsed().as_secs_f64(),
        };
        Ok(Multigrid {
            n,
            options,
            unknown,
            known: known.to_vec(),
            a_uk,
            a: mats,
            p,
            pt,
            diag,
            coarse,
            info,
        })
    }

    /// What the setup built.
    pub fn info(&self) -> &MultigridInfo {
        &self.info
    }

    /// The Galerkin matrices `A_0, …, A_L` (finest first; `A_0` has the fixed rows removed).
    pub fn level_matrices(&self) -> &[SparseMatrix] {
        &self.a
    }

    /// The prolongations after the elimination of the fixed vertices and the dropping of the
    /// unsupported columns.
    pub fn prolongations(&self) -> &[SparseMatrix] {
        &self.p
    }

    /// Number of unknowns on the finest level.
    pub fn n_unknown(&self) -> usize {
        self.unknown.len()
    }

    /// The free vertices of the fine system, ascending.
    pub fn unknown(&self) -> &[u32] {
        &self.unknown
    }

    /// One Gauss–Seidel sweep on level `lv`: forward or backward.
    fn sweep(&self, lv: usize, b: &[f64], x: &mut [f64], forward: bool) {
        let m = &self.a[lv];
        let (cp, ri, v) = (m.col_ptr(), m.row_idx(), m.values());
        let d = &self.diag[lv];
        let n = b.len();
        let mut body = |i: usize| {
            let mut s = 0.0;
            for k in cp[i]..cp[i + 1] {
                let r = ri[k] as usize;
                if r != i {
                    s += v[k] * x[r];
                }
            }
            x[i] = (b[i] - s) / d[i];
        };
        if forward {
            for i in 0..n {
                body(i);
            }
        } else {
            for i in (0..n).rev() {
                body(i);
            }
        }
    }

    fn residual(&self, lv: usize, b: &[f64], x: &[f64]) -> Vec<f64> {
        let mut r = vec![0.0; b.len()];
        // Dimensions are consistent by construction.
        self.a[lv]
            .mul_vec(x, &mut r)
            .expect("level matrices match the vectors");
        for (ri, bi) in r.iter_mut().zip(b) {
            *ri = bi - *ri;
        }
        r
    }

    fn vcycle_level(&self, lv: usize, b: &[f64], x: &mut [f64]) -> Result<(), MultigridError> {
        if lv + 1 == self.a.len() {
            let mut r = self.residual(lv, b, x);
            self.coarse.solve_in_place(&mut r, Parallelism::Sequential)?;
            for (xi, di) in x.iter_mut().zip(&r) {
                *xi += di;
            }
            return Ok(());
        }
        for _ in 0..self.options.pre_sweeps {
            self.sweep(lv, b, x, true);
        }
        let r = self.residual(lv, b, x);
        let mut rc = vec![0.0; self.a[lv + 1].nrows()];
        self.pt[lv].mul_vec(&r, &mut rc)?;
        let mut xc = vec![0.0; rc.len()];
        self.vcycle_level(lv + 1, &rc, &mut xc)?;
        let mut corr = vec![0.0; x.len()];
        self.p[lv].mul_vec(&xc, &mut corr)?;
        for (xi, ci) in x.iter_mut().zip(&corr) {
            *xi += ci;
        }
        for _ in 0..self.options.post_sweeps {
            self.sweep(lv, b, x, false);
        }
        Ok(())
    }

    /// One V-cycle on the reduced system `A_0 x = b` (both of length
    /// [`n_unknown`](Self::n_unknown)), improving `x` in place.
    ///
    /// # Errors
    /// [`MultigridError::SizeMismatch`] for wrong lengths; backend failures.
    pub fn vcycle(&self, b: &[f64], x: &mut [f64]) -> Result<(), MultigridError> {
        self.check_len("right-hand side", b.len())?;
        self.check_len("solution", x.len())?;
        self.vcycle_level(0, b, x)
    }

    /// The preconditioner of conjugate gradients: `z = V(r)`, one V-cycle from the zero guess.
    /// The operator `r ↦ z` is symmetric positive definite.
    ///
    /// # Errors
    /// As [`vcycle`](Self::vcycle).
    pub fn apply_preconditioner(&self, r: &[f64], z: &mut [f64]) -> Result<(), MultigridError> {
        self.check_len("residual", r.len())?;
        self.check_len("output", z.len())?;
        z.fill(0.0);
        self.vcycle_level(0, r, z)
    }

    fn check_len(&self, what: &str, got: usize) -> Result<(), MultigridError> {
        if got != self.unknown.len() {
            return Err(MultigridError::SizeMismatch {
                what: format!("{what} of the reduced system"),
                expected: self.unknown.len(),
                got,
            });
        }
        Ok(())
    }

    /// The reduced right-hand side `b_u − A_uk x_k`, and the initial guess of the unknowns.
    fn reduce(
        &self,
        b: &[f64],
        known_values: &[f64],
        x0: Option<&[f64]>,
    ) -> Result<(Vec<f64>, Vec<f64>), MultigridError> {
        if b.len() != self.n {
            return Err(MultigridError::SizeMismatch {
                what: "right-hand side".into(),
                expected: self.n,
                got: b.len(),
            });
        }
        if known_values.len() != self.known.len() {
            return Err(MultigridError::SizeMismatch {
                what: "fixed values (one per fixed vertex, in the order given to `new`)".into(),
                expected: self.known.len(),
                got: known_values.len(),
            });
        }
        if let Some(x0) = x0
            && x0.len() != self.n
        {
            return Err(MultigridError::SizeMismatch {
                what: "initial guess".into(),
                expected: self.n,
                got: x0.len(),
            });
        }
        if b.iter().chain(known_values).any(|v| !v.is_finite())
            || x0.is_some_and(|x| x.iter().any(|v| !v.is_finite()))
        {
            return Err(MultigridError::NonFinite {
                what: "the right-hand side, fixed values or initial guess".into(),
            });
        }
        let mut rhs: Vec<f64> = self.unknown.iter().map(|&i| b[i as usize]).collect();
        if !self.known.is_empty() {
            let mut t = vec![0.0; rhs.len()];
            self.a_uk.mul_vec(known_values, &mut t)?;
            for (r, ti) in rhs.iter_mut().zip(&t) {
                *r -= ti;
            }
        }
        let z0 = match x0 {
            Some(x) => self.unknown.iter().map(|&i| x[i as usize]).collect(),
            None => vec![0.0; rhs.len()],
        };
        Ok((rhs, z0))
    }

    fn assemble(&self, z: &[f64], known_values: &[f64]) -> Vec<f64> {
        let mut x = vec![0.0; self.n];
        for (&i, &v) in self.unknown.iter().zip(z) {
            x[i as usize] = v;
        }
        for (&i, &v) in self.known.iter().zip(known_values) {
            x[i as usize] = v;
        }
        x
    }

    fn make_report(
        &self,
        residuals: Vec<f64>,
        rhs_norm: f64,
        converged: bool,
        t0: Instant,
    ) -> SolveReport {
        let cycles = residuals.len() - 1;
        let (r0, rn) = (residuals[0], *residuals.last().expect("non-empty"));
        let factor = if cycles > 0 && r0 > 0.0 {
            (rn / r0).powf(1.0 / cycles as f64)
        } else {
            0.0
        };
        SolveReport {
            converged,
            cycles,
            relative_residual: if rhs_norm > 0.0 { rn / rhs_norm } else { rn },
            residuals,
            rhs_norm,
            convergence_factor: factor,
            level_sizes: self.info.level_sizes.clone(),
            solve_seconds: t0.elapsed().as_secs_f64(),
        }
    }

    /// Solves `A x = b` with the fixed values `known_values[i]` at `known[i]` (the order given
    /// to [`new`](Self::new)) by V-cycles. `x0` is an optional initial guess (full length; its
    /// fixed entries are ignored). Returns the full solution and the report; `converged` is
    /// `false` (not an error) if `max_cycles` ran out first.
    ///
    /// # Errors
    /// [`MultigridError`] for wrong lengths or non-finite inputs.
    pub fn solve(
        &self,
        b: &[f64],
        known_values: &[f64],
        x0: Option<&[f64]>,
    ) -> Result<(Vec<f64>, SolveReport), MultigridError> {
        let t0 = Instant::now();
        let (rhs, mut z) = self.reduce(b, known_values, x0)?;
        let mut res = vec![norm(&self.residual(0, &rhs, &z))];
        let rhs_norm = {
            let r = norm(&rhs);
            if r > 0.0 { r } else { res[0] }
        };
        let tol = self.options.tolerance * rhs_norm;
        let mut converged = res[0] <= tol;
        while !converged && res.len() <= self.options.max_cycles {
            self.vcycle_level(0, &rhs, &mut z)?;
            let r = norm(&self.residual(0, &rhs, &z));
            if !r.is_finite() {
                return Err(MultigridError::NonFinite {
                    what: "the residual (the iteration diverged)".into(),
                });
            }
            res.push(r);
            converged = r <= tol;
        }
        let report = self.make_report(res, rhs_norm, converged, t0);
        Ok((self.assemble(&z, known_values), report))
    }

    /// Conjugate gradients on the reduced system, preconditioned by one V-cycle per iteration
    /// (stopping at the relative residual `options.tolerance` or after `max_iterations`).
    /// Same arguments and results as [`solve`](Self::solve); `cycles` counts the iterations.
    ///
    /// # Errors
    /// As [`solve`](Self::solve).
    pub fn pcg(
        &self,
        b: &[f64],
        known_values: &[f64],
        x0: Option<&[f64]>,
        max_iterations: usize,
    ) -> Result<(Vec<f64>, SolveReport), MultigridError> {
        let t0 = Instant::now();
        let (rhs, mut z) = self.reduce(b, known_values, x0)?;
        let n = rhs.len();
        let mut r = self.residual(0, &rhs, &z);
        let mut res = vec![norm(&r)];
        let rhs_norm = {
            let t = norm(&rhs);
            if t > 0.0 { t } else { res[0] }
        };
        let tol = self.options.tolerance * rhs_norm;
        let mut converged = res[0] <= tol;
        let mut s = vec![0.0; n];
        let mut p = vec![0.0; n];
        let mut ap = vec![0.0; n];
        let mut rho = 0.0;
        if !converged {
            self.apply_preconditioner(&r, &mut s)?;
            p.copy_from_slice(&s);
            rho = r.iter().zip(&s).map(|(a, b)| a * b).sum();
        }
        while !converged && res.len() <= max_iterations {
            self.a[0].mul_vec(&p, &mut ap)?;
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if !(pap > 0.0 && pap.is_finite()) {
                return Err(MultigridError::NonFinite {
                    what: "the CG curvature p'Ap (the matrix or the preconditioner is not positive definite)".into(),
                });
            }
            let alpha = rho / pap;
            for i in 0..n {
                z[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            let rn = norm(&r);
            res.push(rn);
            if rn <= tol {
                converged = true;
                break;
            }
            self.apply_preconditioner(&r, &mut s)?;
            let rho_new: f64 = r.iter().zip(&s).map(|(a, b)| a * b).sum();
            let beta = rho_new / rho;
            rho = rho_new;
            for i in 0..n {
                p[i] = s[i] + beta * p[i];
            }
        }
        let report = self.make_report(res, rhs_norm, converged, t0);
        Ok((self.assemble(&z, known_values), report))
    }
}
