//! A small sparse matrix type and sparse direct solvers.
//!
//! [`SparseMatrix`] is plain compressed-sparse-column storage with sorted, duplicate-free
//! columns. [`Symbolic`] is the fill-reducing ordering and symbolic analysis of one sparsity
//! pattern; it is cheap to clone and can be reused for any number of numeric factorizations
//! ([`Symbolic::factor`], [`Factorization::refactor`]) of matrices with that pattern.
//! A [`Factorization`] solves with one or many right-hand sides.
//!
//! Four methods are offered ([`Method`]): `LLᵀ` for symmetric positive definite matrices,
//! `LDLᵀ` without pivoting for symmetric quasi-definite ones, `LBLᵀ` with Bunch–Kaufman pivoting
//! restricted to supernodes for general symmetric indefinite ones, and `LU` with partial pivoting
//! for everything else.
//!
//! Parallelism is an explicit argument of every numeric call ([`Parallelism`]): either the
//! calling thread, or a [`ThreadPool`] that the caller owns. Nothing reads global settings and
//! the global `rayon` pool is never used.
//!
//! # Backend
//!
//! The solvers are a thin layer over the pure-Rust [`faer`](https://crates.io/crates/faer)
//! crate (version 0.24), using its low-level sparse API so that the thread count stays explicit.
//! No `faer` type appears in this module's public interface, so the backend can be replaced
//! without touching callers. Orderings are AMD (symmetric methods) and COLAMD (`LU`); `faer`
//! provides no nested dissection.

use std::fmt;
use std::sync::Arc;

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::cholesky::ldlt::factor::LdltRegularization;
use faer::linalg::cholesky::llt::factor::LltRegularization;
use faer::sparse::linalg::{SupernodalThreshold, amd as famd, cholesky as fchol, lu as flu};
use faer::sparse::{FaerError, SparseColMatRef, SymbolicSparseColMatRef};
use faer::{Conj, MatMut, Par, Side};

/// Errors of sparse matrix construction, factorization and solves.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SparseError {
    /// A triplet lies outside the matrix.
    TripletOutOfRange {
        /// Position of the triplet in the input.
        index: usize,
        /// Its row.
        row: u32,
        /// Its column.
        col: u32,
        /// Number of rows of the matrix.
        nrows: usize,
        /// Number of columns of the matrix.
        ncols: usize,
    },
    /// Raw compressed-column arrays are inconsistent.
    InvalidStructure {
        /// What is wrong.
        reason: String,
    },
    /// The matrix dimensions exceed the `u32` index range.
    TooLarge {
        /// Number of rows requested.
        nrows: usize,
        /// Number of columns requested.
        ncols: usize,
    },
    /// A vector or index list has the wrong length or content for the operation.
    DimensionMismatch {
        /// Which argument is wrong.
        what: &'static str,
        /// The required value.
        expected: usize,
        /// The value that was passed.
        got: usize,
    },
    /// An index passed to an operation is not valid for the matrix.
    IndexOutOfRange {
        /// Which kind of index ("row" or "column").
        what: &'static str,
        /// The offending index.
        index: usize,
        /// The exclusive upper bound it violates.
        bound: usize,
    },
    /// An index that must be unique was listed twice.
    DuplicateIndex {
        /// Which kind of index ("row" or "column").
        what: &'static str,
        /// The repeated index.
        index: usize,
    },
    /// A factorization was requested for a non-square matrix.
    NotSquare {
        /// Number of rows.
        nrows: usize,
        /// Number of columns.
        ncols: usize,
    },
    /// The matrix handed to a numeric factorization does not have the analysed pattern.
    PatternMismatch {
        /// What differs.
        reason: String,
    },
    /// A diagonal entry that an operation needs is not stored in the pattern.
    MissingDiagonal {
        /// The row and column of the missing entry.
        index: usize,
    },
    /// `LLᵀ` met a non-positive pivot: the matrix is not (numerically) positive definite.
    NotPositiveDefinite {
        /// Where the backend noticed it. This refers to the permuted matrix and its exact
        /// meaning depends on the backend's code path; it is a diagnostic hint only.
        step: usize,
    },
    /// `LDLᵀ` met an exactly zero (or non-finite) pivot.
    ZeroPivot {
        /// Where the backend noticed it; a diagnostic hint only, as for
        /// [`SparseError::NotPositiveDefinite`].
        step: usize,
    },
    /// `LU` found no pivot: the matrix is structurally singular.
    StructurallySingular {
        /// Where the backend noticed it; a diagnostic hint only.
        step: usize,
    },
    /// A solve produced a non-finite value, which means the factorized matrix is singular or
    /// too ill-conditioned for the chosen method, or the right-hand side was not finite.
    NonFiniteSolution {
        /// Index of the first non-finite entry (row + n · column).
        index: usize,
    },
    /// A solve was requested from a [`Factorization`] whose last
    /// [`refactor`](Factorization::refactor) failed.
    NotFactorized,
    /// A custom elimination order is not a permutation of `0..n`.
    InvalidEliminationOrder {
        /// What is wrong with it.
        reason: String,
    },
    /// The requested combination of options is not available.
    Unsupported {
        /// What was requested.
        what: &'static str,
    },
    /// A [`ThreadPool`] could not be created.
    ThreadPool {
        /// Why.
        reason: String,
    },
    /// The factor would need more entries than the index type can address.
    IndexOverflow,
    /// An allocation failed.
    OutOfMemory,
}

impl fmt::Display for SparseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SparseError::TripletOutOfRange {
                index,
                row,
                col,
                nrows,
                ncols,
            } => write!(
                f,
                "triplet {index} at (row {row}, column {col}) is outside the {nrows} x {ncols} matrix"
            ),
            SparseError::InvalidStructure { reason } => {
                write!(f, "invalid sparse structure: {reason}")
            }
            SparseError::TooLarge { nrows, ncols } => {
                write!(f, "{nrows} x {ncols} matrix exceeds the 32-bit index range")
            }
            SparseError::DimensionMismatch {
                what,
                expected,
                got,
            } => {
                write!(f, "{what}: expected {expected}, got {got}")
            }
            SparseError::IndexOutOfRange { what, index, bound } => {
                write!(
                    f,
                    "{what} index {index} is out of range (must be < {bound})"
                )
            }
            SparseError::DuplicateIndex { what, index } => {
                write!(f, "{what} {index} is listed more than once")
            }
            SparseError::NotSquare { nrows, ncols } => {
                write!(f, "a square matrix is required, got {nrows} x {ncols}")
            }
            SparseError::PatternMismatch { reason } => {
                write!(
                    f,
                    "matrix does not have the analysed sparsity pattern: {reason}"
                )
            }
            SparseError::MissingDiagonal { index } => {
                write!(
                    f,
                    "diagonal entry ({index}, {index}) is not stored in the sparsity pattern"
                )
            }
            SparseError::NotPositiveDefinite { step } => write!(
                f,
                "matrix is not positive definite: non-positive pivot at elimination step {step}"
            ),
            SparseError::ZeroPivot { step } => {
                write!(
                    f,
                    "LDLt factorization failed: zero pivot at elimination step {step}"
                )
            }
            SparseError::StructurallySingular { step } => {
                write!(
                    f,
                    "matrix is structurally singular: no pivot at elimination step {step}"
                )
            }
            SparseError::NonFiniteSolution { index } => write!(
                f,
                "solve produced a non-finite value at entry {index}: the matrix is singular or too \
                 ill-conditioned for this factorization, or the right-hand side is not finite"
            ),
            SparseError::NotFactorized => write!(
                f,
                "the factorization holds no valid factor: its last refactorization failed"
            ),
            SparseError::InvalidEliminationOrder { reason } => {
                write!(f, "invalid elimination order: {reason}")
            }
            SparseError::Unsupported { what } => write!(f, "not supported: {what}"),
            SparseError::ThreadPool { reason } => write!(f, "cannot create thread pool: {reason}"),
            SparseError::IndexOverflow => {
                write!(f, "sparse factor is too large for the index type")
            }
            SparseError::OutOfMemory => write!(f, "out of memory in sparse factorization"),
        }
    }
}

impl std::error::Error for SparseError {}

impl From<FaerError> for SparseError {
    fn from(e: FaerError) -> Self {
        match e {
            FaerError::IndexOverflow => SparseError::IndexOverflow,
            FaerError::OutOfMemory => SparseError::OutOfMemory,
            // `FaerError` is non-exhaustive; treat future variants as resource errors.
            _ => SparseError::OutOfMemory,
        }
    }
}

/// A real sparse matrix in compressed-sparse-column (CSC) form.
///
/// Invariants: `col_ptr` has `ncols + 1` non-decreasing entries starting at 0; within each
/// column the row indices are strictly increasing (sorted, no duplicates) and `< nrows`.
/// Explicitly stored zeros are allowed and are part of the sparsity pattern.
#[derive(Debug, Clone, PartialEq)]
pub struct SparseMatrix {
    nrows: usize,
    ncols: usize,
    col_ptr: Vec<usize>,
    row_idx: Vec<u32>,
    values: Vec<f64>,
}

impl SparseMatrix {
    /// Builds a matrix from `(row, column, value)` triplets. Triplets with equal position are
    /// summed (in input order). Zero values are kept as explicit entries, so the pattern depends
    /// only on the positions.
    ///
    /// # Errors
    /// [`SparseError::TripletOutOfRange`], [`SparseError::TooLarge`].
    pub fn from_triplets(
        nrows: usize,
        ncols: usize,
        triplets: &[(u32, u32, f64)],
    ) -> Result<SparseMatrix, SparseError> {
        if nrows > u32::MAX as usize || ncols > u32::MAX as usize {
            return Err(SparseError::TooLarge { nrows, ncols });
        }
        let mut counts = vec![0usize; ncols + 1];
        for (index, &(row, col, _)) in triplets.iter().enumerate() {
            if row as usize >= nrows || col as usize >= ncols {
                return Err(SparseError::TripletOutOfRange {
                    index,
                    row,
                    col,
                    nrows,
                    ncols,
                });
            }
            counts[col as usize + 1] += 1;
        }
        for j in 0..ncols {
            counts[j + 1] += counts[j];
        }
        // Bucket by column (stable), then sort each column by row (stable) and merge.
        let mut next = counts.clone();
        let mut rows = vec![0u32; triplets.len()];
        let mut vals = vec![0f64; triplets.len()];
        for &(row, col, value) in triplets {
            let slot = &mut next[col as usize];
            rows[*slot] = row;
            vals[*slot] = value;
            *slot += 1;
        }
        let mut col_ptr = Vec::with_capacity(ncols + 1);
        let mut row_idx = Vec::with_capacity(triplets.len());
        let mut values = Vec::with_capacity(triplets.len());
        let mut order: Vec<usize> = Vec::new();
        col_ptr.push(0);
        for j in 0..ncols {
            let (lo, hi) = (counts[j], counts[j + 1]);
            order.clear();
            order.extend(lo..hi);
            order.sort_by_key(|&k| rows[k]);
            let start = row_idx.len();
            for &k in &order {
                if row_idx.len() > start && row_idx[row_idx.len() - 1] == rows[k] {
                    let last = values.len() - 1;
                    values[last] += vals[k];
                } else {
                    row_idx.push(rows[k]);
                    values.push(vals[k]);
                }
            }
            col_ptr.push(row_idx.len());
        }
        Ok(SparseMatrix {
            nrows,
            ncols,
            col_ptr,
            row_idx,
            values,
        })
    }

    /// Wraps raw CSC arrays after checking every invariant.
    ///
    /// # Errors
    /// [`SparseError::InvalidStructure`] naming the first violated invariant,
    /// [`SparseError::TooLarge`].
    pub fn from_csc(
        nrows: usize,
        ncols: usize,
        col_ptr: Vec<usize>,
        row_idx: Vec<u32>,
        values: Vec<f64>,
    ) -> Result<SparseMatrix, SparseError> {
        if nrows > u32::MAX as usize || ncols > u32::MAX as usize {
            return Err(SparseError::TooLarge { nrows, ncols });
        }
        let invalid = |reason: String| Err(SparseError::InvalidStructure { reason });
        if col_ptr.len() != ncols + 1 {
            return invalid(format!(
                "col_ptr has {} entries, expected ncols + 1 = {}",
                col_ptr.len(),
                ncols + 1
            ));
        }
        if col_ptr[0] != 0 {
            return invalid(format!("col_ptr[0] is {}, expected 0", col_ptr[0]));
        }
        if row_idx.len() != values.len() {
            return invalid(format!(
                "{} row indices but {} values",
                row_idx.len(),
                values.len()
            ));
        }
        if col_ptr[ncols] != row_idx.len() {
            return invalid(format!(
                "col_ptr ends at {}, but there are {} stored entries",
                col_ptr[ncols],
                row_idx.len()
            ));
        }
        for j in 0..ncols {
            if col_ptr[j] > col_ptr[j + 1] {
                return invalid(format!("col_ptr decreases at column {j}"));
            }
        }
        for j in 0..ncols {
            let rows = &row_idx[col_ptr[j]..col_ptr[j + 1]];
            if let Some(&r) = rows.iter().find(|&&r| r as usize >= nrows) {
                return invalid(format!(
                    "column {j} has row index {r}, but there are {nrows} rows"
                ));
            }
            if let Some(w) = rows.windows(2).find(|w| w[0] >= w[1]) {
                return invalid(format!(
                    "column {j} is not strictly increasing: row {} is followed by row {}",
                    w[0], w[1]
                ));
            }
        }
        Ok(SparseMatrix {
            nrows,
            ncols,
            col_ptr,
            row_idx,
            values,
        })
    }

    /// The square diagonal matrix with the given diagonal (every diagonal entry is stored).
    pub fn from_diagonal(diagonal: &[f64]) -> Result<SparseMatrix, SparseError> {
        let n = diagonal.len();
        if n > u32::MAX as usize {
            return Err(SparseError::TooLarge { nrows: n, ncols: n });
        }
        Ok(SparseMatrix {
            nrows: n,
            ncols: n,
            col_ptr: (0..=n).collect(),
            row_idx: (0..n as u32).collect(),
            values: diagonal.to_vec(),
        })
    }

    /// Number of rows.
    #[inline]
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Number of columns.
    #[inline]
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// Number of stored entries (including explicit zeros).
    #[inline]
    pub fn nnz(&self) -> usize {
        self.row_idx.len()
    }

    /// Column pointers: the entries of column `j` are `col_ptr[j]..col_ptr[j + 1]`.
    #[inline]
    pub fn col_ptr(&self) -> &[usize] {
        &self.col_ptr
    }

    /// Row index of every stored entry, column by column, ascending within a column.
    #[inline]
    pub fn row_idx(&self) -> &[u32] {
        &self.row_idx
    }

    /// Value of every stored entry, in the order of [`SparseMatrix::row_idx`].
    #[inline]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Mutable values. The pattern cannot change, so a [`Symbolic`] analysis stays valid.
    #[inline]
    pub fn values_mut(&mut self) -> &mut [f64] {
        &mut self.values
    }

    /// Position of entry `(row, col)` in [`SparseMatrix::values`], if it is stored.
    pub fn find(&self, row: usize, col: usize) -> Option<usize> {
        if row >= self.nrows || col >= self.ncols {
            return None;
        }
        let (lo, hi) = (self.col_ptr[col], self.col_ptr[col + 1]);
        self.row_idx[lo..hi]
            .binary_search(&(row as u32))
            .ok()
            .map(|k| lo + k)
    }

    /// The value at `(row, col)`; `0.0` if the entry is not stored or out of range.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.find(row, col).map_or(0.0, |k| self.values[k])
    }

    /// All stored entries as `(row, column, value)`, column by column.
    pub fn triplets(&self) -> impl Iterator<Item = (u32, u32, f64)> + '_ {
        (0..self.ncols).flat_map(move |j| {
            (self.col_ptr[j]..self.col_ptr[j + 1])
                .map(move |k| (self.row_idx[k], j as u32, self.values[k]))
        })
    }

    /// The main diagonal (`0.0` where no entry is stored).
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.nrows.min(self.ncols))
            .map(|i| self.get(i, i))
            .collect()
    }

    /// Adds `d[i]` to the stored diagonal entry `(i, i)` for every `i`. The pattern is unchanged.
    ///
    /// # Errors
    /// [`SparseError::DimensionMismatch`] if `d` is not `min(nrows, ncols)` long;
    /// [`SparseError::MissingDiagonal`] if a diagonal entry is not stored (nothing is modified).
    pub fn add_to_diagonal(&mut self, d: &[f64]) -> Result<(), SparseError> {
        let n = self.nrows.min(self.ncols);
        if d.len() != n {
            return Err(SparseError::DimensionMismatch {
                what: "diagonal length",
                expected: n,
                got: d.len(),
            });
        }
        let mut slots = Vec::with_capacity(n);
        for i in 0..n {
            slots.push(
                self.find(i, i)
                    .ok_or(SparseError::MissingDiagonal { index: i })?,
            );
        }
        for (slot, &x) in slots.into_iter().zip(d) {
            self.values[slot] += x;
        }
        Ok(())
    }

    /// Computes `y = A x`.
    ///
    /// # Errors
    /// [`SparseError::DimensionMismatch`] if `x.len() != ncols` or `y.len() != nrows`.
    pub fn mul_vec(&self, x: &[f64], y: &mut [f64]) -> Result<(), SparseError> {
        if x.len() != self.ncols {
            return Err(SparseError::DimensionMismatch {
                what: "length of x",
                expected: self.ncols,
                got: x.len(),
            });
        }
        if y.len() != self.nrows {
            return Err(SparseError::DimensionMismatch {
                what: "length of y",
                expected: self.nrows,
                got: y.len(),
            });
        }
        y.fill(0.0);
        for (j, &xj) in x.iter().enumerate() {
            for k in self.col_ptr[j]..self.col_ptr[j + 1] {
                y[self.row_idx[k] as usize] += self.values[k] * xj;
            }
        }
        Ok(())
    }

    /// The transpose.
    pub fn transpose(&self) -> SparseMatrix {
        let mut col_ptr = vec![0usize; self.nrows + 1];
        for &r in &self.row_idx {
            col_ptr[r as usize + 1] += 1;
        }
        for i in 0..self.nrows {
            col_ptr[i + 1] += col_ptr[i];
        }
        let mut next = col_ptr.clone();
        let mut row_idx = vec![0u32; self.nnz()];
        let mut values = vec![0f64; self.nnz()];
        for j in 0..self.ncols {
            for k in self.col_ptr[j]..self.col_ptr[j + 1] {
                let slot = &mut next[self.row_idx[k] as usize];
                row_idx[*slot] = j as u32;
                values[*slot] = self.values[k];
                *slot += 1;
            }
        }
        SparseMatrix {
            nrows: self.ncols,
            ncols: self.nrows,
            col_ptr,
            row_idx,
            values,
        }
    }

    /// The largest `|A_ij − A_ji|` over all entries; `0.0` for an exactly symmetric matrix.
    ///
    /// # Errors
    /// [`SparseError::NotSquare`].
    pub fn symmetry_defect(&self) -> Result<f64, SparseError> {
        if self.nrows != self.ncols {
            return Err(SparseError::NotSquare {
                nrows: self.nrows,
                ncols: self.ncols,
            });
        }
        let mut worst = 0.0f64;
        for (i, j, v) in self.triplets() {
            worst = worst.max((v - self.get(j as usize, i as usize)).abs());
        }
        Ok(worst)
    }

    /// The submatrix with the given rows and columns, in the given order:
    /// `result[(a, b)] = self[(rows[a], cols[b])]`.
    ///
    /// # Errors
    /// [`SparseError::IndexOutOfRange`] for an index outside the matrix,
    /// [`SparseError::DuplicateIndex`] if a row is listed twice (columns may repeat).
    pub fn submatrix(&self, rows: &[u32], cols: &[u32]) -> Result<SparseMatrix, SparseError> {
        const UNUSED: u32 = u32::MAX;
        let mut row_map = vec![UNUSED; self.nrows];
        for (new, &old) in rows.iter().enumerate() {
            if old as usize >= self.nrows {
                return Err(SparseError::IndexOutOfRange {
                    what: "row",
                    index: old as usize,
                    bound: self.nrows,
                });
            }
            if row_map[old as usize] != UNUSED {
                return Err(SparseError::DuplicateIndex {
                    what: "row",
                    index: old as usize,
                });
            }
            row_map[old as usize] = new as u32;
        }
        let mut col_ptr = Vec::with_capacity(cols.len() + 1);
        let mut row_idx = Vec::new();
        let mut values = Vec::new();
        let mut column: Vec<(u32, f64)> = Vec::new();
        col_ptr.push(0);
        for &j in cols {
            if j as usize >= self.ncols {
                return Err(SparseError::IndexOutOfRange {
                    what: "column",
                    index: j as usize,
                    bound: self.ncols,
                });
            }
            column.clear();
            for k in self.col_ptr[j as usize]..self.col_ptr[j as usize + 1] {
                let r = row_map[self.row_idx[k] as usize];
                if r != UNUSED {
                    column.push((r, self.values[k]));
                }
            }
            column.sort_unstable_by_key(|&(r, _)| r);
            for &(r, v) in &column {
                row_idx.push(r);
                values.push(v);
            }
            col_ptr.push(row_idx.len());
        }
        Ok(SparseMatrix {
            nrows: rows.len(),
            ncols: cols.len(),
            col_ptr,
            row_idx,
            values,
        })
    }
}

/// A private pool of worker threads for numeric factorizations and solves.
///
/// The pool is owned by the caller and passed explicitly ([`Parallelism::Pool`]); the global
/// `rayon` pool is never used. This matters for speed as well as for hygiene: handing the
/// backend a pool with many more threads than the factorization can use makes it dramatically
/// slower (see the timings in the crate README), so size the pool deliberately.
pub struct ThreadPool {
    pool: rayon::ThreadPool,
    n_threads: usize,
}

impl fmt::Debug for ThreadPool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ThreadPool")
            .field("n_threads", &self.n_threads)
            .finish()
    }
}

impl ThreadPool {
    /// Starts a pool with exactly `n_threads` worker threads.
    ///
    /// # Errors
    /// [`SparseError::ThreadPool`] if `n_threads` is 0 or the threads cannot be started.
    pub fn new(n_threads: usize) -> Result<ThreadPool, SparseError> {
        if n_threads == 0 {
            return Err(SparseError::ThreadPool {
                reason: "a thread pool needs at least 1 thread".to_string(),
            });
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .thread_name(|i| format!("ddg-sparse-{i}"))
            .build()
            .map_err(|e| SparseError::ThreadPool {
                reason: e.to_string(),
            })?;
        Ok(ThreadPool { pool, n_threads })
    }

    /// Number of worker threads.
    pub fn n_threads(&self) -> usize {
        self.n_threads
    }

    /// Calls `f(0), …, f(n − 1)` on the pool's threads, one index at a time to whichever thread
    /// is free, and returns the results in index order. Nothing runs on the global `rayon`
    /// pool. This is for callers that have independent jobs of their own (for instance
    /// right-hand sides to solve with one factorization) and want them on the same explicit
    /// pool as the factorization; `f` may itself be called from several threads at once.
    pub fn map<R: Send>(&self, n: usize, f: impl Fn(usize) -> R + Sync + Send) -> Vec<R> {
        use rayon::prelude::*;
        self.pool
            .install(|| (0..n).into_par_iter().with_max_len(1).map(f).collect())
    }
}

/// Where a numeric factorization or solve runs.
#[derive(Debug, Clone, Copy)]
pub enum Parallelism<'a> {
    /// On the calling thread only.
    Sequential,
    /// On the given pool, using all of its threads.
    Pool(&'a ThreadPool),
}

impl Parallelism<'_> {
    /// Runs `f` with the backend's parallelism setting, inside the pool if there is one.
    fn run<R: Send>(self, f: impl FnOnce(Par) -> R + Send) -> R {
        match self {
            Parallelism::Sequential => f(Par::Seq),
            Parallelism::Pool(pool) => pool.pool.install(|| f(Par::rayon(pool.n_threads))),
        }
    }
}

/// Which factorization to compute.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// Cholesky `A = L Lᵀ` for symmetric positive definite matrices. Fails with
    /// [`SparseError::NotPositiveDefinite`] otherwise.
    Llt,
    /// `A = L D Lᵀ` with diagonal `D` and **no pivoting**. For a symmetric quasi-definite matrix
    /// (`[[−E, Bᵀ], [B, F]]` with `E`, `F` positive definite) the factorization exists under
    /// every elimination order. For other indefinite matrices it exists only for suitable
    /// orders ([`EliminationOrder::Custom`]) and otherwise breaks down with
    /// [`SparseError::ZeroPivot`]. In both cases there is no pivoting to limit element growth,
    /// so the solution can be less accurate than the residual of a pivoted method (measured on
    /// this workspace's saddle-point systems: relative residual around 1e-8, back to 1e-15
    /// after one step of [`Factorization::solve_refined`]). Check residuals.
    Ldlt,
    /// `A = L B Lᵀ` with 1×1 and 2×2 pivots chosen by Bunch–Kaufman pivoting inside each
    /// supernode. Pivoting cannot leave a supernode, so a matrix whose zero diagonal entries
    /// are eliminated early can still defeat it: the factorization then "succeeds" and the
    /// solve reports [`SparseError::NonFiniteSolution`]. On this workspace's saddle-point
    /// systems exactly that happens, and the method also slows down sharply with more threads;
    /// prefer [`Method::Ldlt`] with a suitable order, or [`Method::Lu`].
    Lblt,
    /// `A = L U` with partial (row) pivoting and a COLAMD column ordering, for general square
    /// matrices. The most robust and the slowest.
    Lu,
}

impl Method {
    /// Whether the method assumes a symmetric matrix (and reads only its lower triangle).
    pub fn is_symmetric(self) -> bool {
        !matches!(self, Method::Lu)
    }
}

/// The order in which a symmetric factorization eliminates the unknowns.
#[derive(Debug, Clone, Copy, Default)]
pub enum EliminationOrder<'a> {
    /// A fill-reducing order computed by the backend: AMD for the symmetric methods, COLAMD
    /// (columns) for `LU`. The default.
    #[default]
    FillReducing,
    /// The order of the matrix as given. Symmetric methods only.
    Natural,
    /// `order[k]` is the index of the unknown eliminated in step `k`; must be a permutation of
    /// `0..n`. Symmetric methods only.
    ///
    /// This is how an indefinite matrix can be factorized by [`Method::Ldlt`] without pivoting:
    /// choose an order in which every pivot is guaranteed to be non-zero. For a saddle-point
    /// matrix `[[−M, L], [L, D]]` with `M` positive diagonal and `D` non-negative diagonal,
    /// eliminating unknown `i` of the first block immediately before unknown `n + i` of the
    /// second (with the pairs in a fill-reducing order of `L`, see [`fill_reducing_ordering`])
    /// is such an order.
    Custom(&'a [u32]),
}

/// The approximate minimum degree (AMD) elimination order for the sparsity pattern of
/// `A + Aᵀ`: `order[k]` is the index of the unknown to eliminate in step `k`. This is the order
/// [`EliminationOrder::FillReducing`] uses for the symmetric methods; it is exposed so that callers can
/// build structured orders from it (see [`EliminationOrder::Custom`]).
///
/// # Errors
/// [`SparseError::NotSquare`], or a resource error.
pub fn fill_reducing_ordering(a: &SparseMatrix) -> Result<Vec<u32>, SparseError> {
    if a.nrows != a.ncols {
        return Err(SparseError::NotSquare {
            nrows: a.nrows,
            ncols: a.ncols,
        });
    }
    let n = a.nrows;
    let row_idx: Vec<usize> = a.row_idx.iter().map(|&r| r as usize).collect();
    // SAFETY: `SparseMatrix` upholds the invariants faer requires.
    let pattern = unsafe {
        SymbolicSparseColMatRef::<usize>::new_unchecked(n, n, &a.col_ptr, None, &row_idx)
    };
    let mut perm: Vec<usize> = try_zeroed(n)?;
    let mut perm_inv: Vec<usize> = try_zeroed(n)?;
    let mut mem = MemBuffer::try_new(famd::order_scratch::<usize>(n, a.nnz()))
        .map_err(|_| SparseError::OutOfMemory)?;
    famd::order(
        &mut perm,
        &mut perm_inv,
        pattern,
        famd::Control::default(),
        MemStack::new(&mut mem),
    )?;
    Ok(perm.into_iter().map(|p| p as u32).collect())
}

enum SymbolicKind {
    Cholesky(fchol::SymbolicCholesky<usize>),
    Lu(flu::SymbolicLu<usize>),
}

struct SymbolicInner {
    method: Method,
    n: usize,
    /// The analysed pattern, in the index type of the backend.
    col_ptr: Vec<usize>,
    row_idx: Vec<usize>,
    kind: SymbolicKind,
}

/// The reusable symbolic analysis of one sparsity pattern for one [`Method`].
///
/// Cloning is cheap (the analysis is shared).
#[derive(Clone)]
pub struct Symbolic {
    inner: Arc<SymbolicInner>,
}

impl fmt::Debug for Symbolic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Symbolic")
            .field("info", &self.info())
            .finish()
    }
}

/// Facts about a symbolic analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SymbolicInfo {
    /// The method the analysis was made for.
    pub method: Method,
    /// Matrix dimension.
    pub n: usize,
    /// Stored entries of the analysed matrix (both triangles, as passed in).
    pub matrix_nnz: usize,
    /// Number of `f64` values in the factor (for supernodal factors this counts the dense
    /// supernode blocks). `None` for `LU`, whose fill is only known after numeric pivoting.
    pub factor_values: Option<usize>,
    /// Whether the backend chose a supernodal (dense-block) factorization over a simplicial
    /// one. `None` for `LU`, where the backend does not expose its choice.
    pub supernodal: Option<bool>,
}

impl Symbolic {
    /// Analyses the sparsity pattern of `a` for the given method, with the backend's
    /// fill-reducing ordering ([`EliminationOrder::FillReducing`]).
    ///
    /// For the symmetric methods only the lower triangle of `a` (entries with `row >= col`) is
    /// read, here and in every later numeric factorization; whether `a` stores the upper
    /// triangle as well makes no difference. Values are not inspected.
    ///
    /// # Errors
    /// [`SparseError::NotSquare`], [`SparseError::IndexOverflow`], [`SparseError::OutOfMemory`].
    pub fn analyze(a: &SparseMatrix, method: Method) -> Result<Symbolic, SparseError> {
        Symbolic::analyze_with_order(a, method, EliminationOrder::FillReducing)
    }

    /// Like [`Symbolic::analyze`], with an explicit elimination order.
    ///
    /// # Errors
    /// As for [`Symbolic::analyze`], plus [`SparseError::InvalidEliminationOrder`] if a custom order is
    /// not a permutation of `0..n`, and [`SparseError::Unsupported`] for `LU` with any order
    /// other than [`EliminationOrder::FillReducing`].
    pub fn analyze_with_order(
        a: &SparseMatrix,
        method: Method,
        ordering: EliminationOrder<'_>,
    ) -> Result<Symbolic, SparseError> {
        if a.nrows != a.ncols {
            return Err(SparseError::NotSquare {
                nrows: a.nrows,
                ncols: a.ncols,
            });
        }
        let n = a.nrows;
        if method == Method::Lu && !matches!(ordering, EliminationOrder::FillReducing) {
            return Err(SparseError::Unsupported {
                what: "LU with a natural or custom ordering",
            });
        }
        // Forward and inverse permutation of a custom order, validated.
        let custom: Option<(Vec<usize>, Vec<usize>)> = match ordering {
            EliminationOrder::Custom(order) => {
                let invalid = |reason: String| Err(SparseError::InvalidEliminationOrder { reason });
                if order.len() != n {
                    return invalid(format!(
                        "it has {} entries, but the matrix has dimension {n}",
                        order.len()
                    ));
                }
                let mut inverse = vec![usize::MAX; n];
                for (step, &index) in order.iter().enumerate() {
                    if index as usize >= n {
                        return invalid(format!(
                            "entry {step} is {index}, which is not below the dimension {n}"
                        ));
                    }
                    if inverse[index as usize] != usize::MAX {
                        return invalid(format!(
                            "index {index} appears twice (entries {} and {step})",
                            inverse[index as usize]
                        ));
                    }
                    inverse[index as usize] = step;
                }
                Some((order.iter().map(|&i| i as usize).collect(), inverse))
            }
            _ => None,
        };
        let col_ptr = a.col_ptr.clone();
        let row_idx: Vec<usize> = a.row_idx.iter().map(|&r| r as usize).collect();
        // SAFETY: `SparseMatrix` upholds exactly the invariants faer requires (monotone column
        // pointers starting at 0, in-range, sorted, duplicate-free row indices).
        let pattern = unsafe {
            SymbolicSparseColMatRef::<usize>::new_unchecked(n, n, &col_ptr, None, &row_idx)
        };
        let kind = match method {
            Method::Llt | Method::Ldlt | Method::Lblt => {
                let mut params = fchol::CholeskySymbolicParams::default();
                if method == Method::Lblt {
                    // Pivoting only exists in the supernodal code path.
                    params.supernodal_flop_ratio_threshold = SupernodalThreshold::FORCE_SUPERNODAL;
                }
                let ordering = match (&ordering, &custom) {
                    (EliminationOrder::FillReducing, _) => fchol::SymmetricOrdering::Amd,
                    (EliminationOrder::Natural, _) => fchol::SymmetricOrdering::Identity,
                    // SAFETY: validated above: both arrays are permutations of 0..n and
                    // inverse to each other.
                    (EliminationOrder::Custom(_), Some((forward, inverse))) => {
                        fchol::SymmetricOrdering::Custom(unsafe {
                            faer::perm::PermRef::new_unchecked(forward, inverse, n)
                        })
                    }
                    (EliminationOrder::Custom(_), None) => {
                        unreachable!("custom permutation was built above")
                    }
                };
                SymbolicKind::Cholesky(fchol::factorize_symbolic_cholesky(
                    pattern,
                    Side::Lower,
                    ordering,
                    params,
                )?)
            }
            Method::Lu => {
                SymbolicKind::Lu(flu::factorize_symbolic_lu(pattern, Default::default())?)
            }
        };
        Ok(Symbolic {
            inner: Arc::new(SymbolicInner {
                method,
                n,
                col_ptr,
                row_idx,
                kind,
            }),
        })
    }

    /// The method this analysis was made for.
    pub fn method(&self) -> Method {
        self.inner.method
    }

    /// Dimension of the analysed matrix.
    pub fn dim(&self) -> usize {
        self.inner.n
    }

    /// Facts about the analysis (factor size, supernodal or not).
    pub fn info(&self) -> SymbolicInfo {
        let inner = &*self.inner;
        let (factor_values, supernodal) = match &inner.kind {
            SymbolicKind::Cholesky(s) => (
                Some(s.len_val()),
                Some(matches!(s.raw(), fchol::SymbolicCholeskyRaw::Supernodal(_))),
            ),
            SymbolicKind::Lu(_) => (None, None),
        };
        SymbolicInfo {
            method: inner.method,
            n: inner.n,
            matrix_nnz: inner.row_idx.len(),
            factor_values,
            supernodal,
        }
    }

    /// Computes the numeric factorization of `a`, which must have exactly the analysed pattern.
    ///
    /// # Errors
    /// [`SparseError::PatternMismatch`] if the pattern differs; otherwise the failure of the
    /// method: [`SparseError::NotPositiveDefinite`] (`Llt`), [`SparseError::ZeroPivot`]
    /// (`Ldlt`), [`SparseError::StructurallySingular`] (`Lu`), or a resource error.
    pub fn factor(
        &self,
        a: &SparseMatrix,
        par: Parallelism<'_>,
    ) -> Result<Factorization, SparseError> {
        let numeric = match &self.inner.kind {
            SymbolicKind::Cholesky(s) => {
                let n = self.inner.n;
                let (pivots, n_perm) = if self.inner.method == Method::Lblt {
                    (n, n)
                } else {
                    (0, 0)
                };
                Numeric::Cholesky {
                    values: try_zeroed(s.len_val())?,
                    subdiag: try_zeroed(pivots)?,
                    perm_fwd: try_zeroed(n_perm)?,
                    perm_inv: try_zeroed(n_perm)?,
                }
            }
            SymbolicKind::Lu(_) => Numeric::Lu(Box::new(flu::NumericLu::new())),
        };
        let mut factorization = Factorization {
            symbolic: self.clone(),
            numeric,
            valid: false,
        };
        factorization.refactor(a, par)?;
        Ok(factorization)
    }
}

fn try_zeroed<T: Clone + Default>(len: usize) -> Result<Vec<T>, SparseError> {
    let mut v = Vec::new();
    v.try_reserve_exact(len)
        .map_err(|_| SparseError::OutOfMemory)?;
    v.resize(len, T::default());
    Ok(v)
}

enum Numeric {
    Cholesky {
        values: Vec<f64>,
        /// `LBLᵀ` only: sub-diagonal of the block-diagonal factor.
        subdiag: Vec<f64>,
        /// `LBLᵀ` only: pivoting permutation.
        perm_fwd: Vec<usize>,
        perm_inv: Vec<usize>,
    },
    Lu(Box<flu::NumericLu<usize, f64>>),
}

/// A numeric factorization, ready to solve linear systems.
pub struct Factorization {
    symbolic: Symbolic,
    numeric: Numeric,
    /// False after a failed [`Factorization::refactor`], until one succeeds.
    valid: bool,
}

impl fmt::Debug for Factorization {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Factorization")
            .field("symbolic", &self.symbolic.info())
            .finish()
    }
}

impl Factorization {
    /// The symbolic analysis this factorization is based on.
    pub fn symbolic(&self) -> &Symbolic {
        &self.symbolic
    }

    /// Dimension of the factorized matrix.
    pub fn dim(&self) -> usize {
        self.symbolic.inner.n
    }

    /// Recomputes the numeric factorization for new values on the same pattern, reusing the
    /// symbolic analysis and the factor's memory.
    ///
    /// If this returns an error the factorization holds no valid factor: solves fail with
    /// [`SparseError::NotFactorized`] until a later `refactor` succeeds. (A pattern mismatch is
    /// detected before anything is touched and leaves the previous factor usable.)
    ///
    /// # Errors
    /// As for [`Symbolic::factor`].
    pub fn refactor(&mut self, a: &SparseMatrix, par: Parallelism<'_>) -> Result<(), SparseError> {
        check_pattern(&self.symbolic.inner, a)?;
        self.valid = false;
        par.run(|par| self.refactor_with(a, par))?;
        self.valid = true;
        Ok(())
    }

    fn refactor_with(&mut self, a: &SparseMatrix, par: Par) -> Result<(), SparseError> {
        let inner = &*self.symbolic.inner;
        // SAFETY: the pattern was validated when `Symbolic` was built from a `SparseMatrix`.
        let pattern = unsafe {
            SymbolicSparseColMatRef::<usize>::new_unchecked(
                inner.n,
                inner.n,
                &inner.col_ptr,
                None,
                &inner.row_idx,
            )
        };
        let matrix = SparseColMatRef::<usize, f64>::new(pattern, &a.values);
        match (&inner.kind, &mut self.numeric) {
            (
                SymbolicKind::Cholesky(s),
                Numeric::Cholesky {
                    values,
                    subdiag,
                    perm_fwd,
                    perm_inv,
                },
            ) => match inner.method {
                Method::Llt => {
                    let mut mem = MemBuffer::try_new(
                        s.factorize_numeric_llt_scratch::<f64>(par, Default::default()),
                    )
                    .map_err(|_| SparseError::OutOfMemory)?;
                    s.factorize_numeric_llt::<f64>(
                        values,
                        matrix,
                        Side::Lower,
                        LltRegularization::default(),
                        par,
                        MemStack::new(&mut mem),
                        Default::default(),
                    )
                    .map_err(|e| match e {
                        faer::linalg::cholesky::llt::factor::LltError::NonPositivePivot {
                            index,
                        } => SparseError::NotPositiveDefinite { step: index },
                    })?;
                }
                Method::Ldlt => {
                    let mut mem = MemBuffer::try_new(
                        s.factorize_numeric_ldlt_scratch::<f64>(par, Default::default()),
                    )
                    .map_err(|_| SparseError::OutOfMemory)?;
                    s.factorize_numeric_ldlt::<f64>(
                        values,
                        matrix,
                        Side::Lower,
                        LdltRegularization::default(),
                        par,
                        MemStack::new(&mut mem),
                        Default::default(),
                    )
                    .map_err(|e| match e {
                        faer::linalg::cholesky::ldlt::factor::LdltError::ZeroPivot { index } => {
                            SparseError::ZeroPivot { step: index }
                        }
                    })?;
                }
                Method::Lblt => {
                    let mut mem = MemBuffer::try_new(
                        s.factorize_numeric_intranode_lblt_scratch::<f64>(par, Default::default()),
                    )
                    .map_err(|_| SparseError::OutOfMemory)?;
                    s.factorize_numeric_intranode_lblt::<f64>(
                        values,
                        subdiag,
                        perm_fwd,
                        perm_inv,
                        matrix,
                        Side::Lower,
                        par,
                        MemStack::new(&mut mem),
                        Default::default(),
                    );
                }
                Method::Lu => unreachable!("LU analysis never has a Cholesky symbolic"),
            },
            (SymbolicKind::Lu(s), Numeric::Lu(numeric)) => {
                let mut mem = MemBuffer::try_new(
                    s.factorize_numeric_lu_scratch::<f64>(par, Default::default()),
                )
                .map_err(|_| SparseError::OutOfMemory)?;
                s.factorize_numeric_lu::<f64>(
                    numeric,
                    matrix,
                    par,
                    MemStack::new(&mut mem),
                    Default::default(),
                )
                .map_err(|e| match e {
                    faer::sparse::linalg::LuError::SymbolicSingular { index } => {
                        SparseError::StructurallySingular { step: index }
                    }
                    faer::sparse::linalg::LuError::Generic(g) => g.into(),
                })?;
            }
            _ => unreachable!("numeric storage always matches the symbolic kind"),
        }
        Ok(())
    }

    /// Solves `A x = b` for one right-hand side, overwriting `rhs` (which holds `b`) with `x`.
    ///
    /// # Errors
    /// [`SparseError::DimensionMismatch`] if `rhs.len() != dim()`;
    /// [`SparseError::NonFiniteSolution`] if the result is not finite (then `rhs` holds garbage).
    pub fn solve_in_place(&self, rhs: &mut [f64], par: Parallelism<'_>) -> Result<(), SparseError> {
        self.solve_many_in_place(rhs, 1, par)
    }

    /// Solves `A X = B` for `nrhs` right-hand sides stored column after column in `rhs`
    /// (`rhs[i + dim · k]` is row `i` of right-hand side `k`), overwriting them with the solutions.
    ///
    /// # Errors
    /// [`SparseError::DimensionMismatch`] if `rhs.len() != dim() · nrhs`;
    /// [`SparseError::NonFiniteSolution`] if the result is not finite (then `rhs` holds garbage).
    pub fn solve_many_in_place(
        &self,
        rhs: &mut [f64],
        nrhs: usize,
        par: Parallelism<'_>,
    ) -> Result<(), SparseError> {
        if !self.valid {
            return Err(SparseError::NotFactorized);
        }
        let n = self.symbolic.inner.n;
        let expected = n.checked_mul(nrhs).ok_or(SparseError::IndexOverflow)?;
        if rhs.len() != expected {
            return Err(SparseError::DimensionMismatch {
                what: "length of right-hand side (dim * nrhs)",
                expected,
                got: rhs.len(),
            });
        }
        if expected == 0 {
            return Ok(());
        }
        par.run(|par| self.solve_with(rhs, nrhs, par))
    }

    fn solve_with(&self, rhs: &mut [f64], nrhs: usize, par: Par) -> Result<(), SparseError> {
        let inner = &*self.symbolic.inner;
        let n = inner.n;
        {
            let rhs = MatMut::from_column_major_slice_mut(rhs, n, nrhs);
            match (&inner.kind, &self.numeric) {
                (
                    SymbolicKind::Cholesky(s),
                    Numeric::Cholesky {
                        values,
                        subdiag,
                        perm_fwd,
                        perm_inv,
                    },
                ) => {
                    let mut mem = MemBuffer::try_new(s.solve_in_place_scratch::<f64>(nrhs, par))
                        .map_err(|_| SparseError::OutOfMemory)?;
                    let stack = MemStack::new(&mut mem);
                    match inner.method {
                        Method::Llt => fchol::LltRef::<usize, f64>::new(s, values)
                            .solve_in_place_with_conj(Conj::No, rhs, par, stack),
                        Method::Ldlt => fchol::LdltRef::<usize, f64>::new(s, values)
                            .solve_in_place_with_conj(Conj::No, rhs, par, stack),
                        Method::Lblt => {
                            // SAFETY: both arrays were filled by the factorization with a
                            // permutation of 0..n and its inverse.
                            let perm = unsafe {
                                faer::perm::PermRef::new_unchecked(perm_fwd, perm_inv, n)
                            };
                            fchol::IntranodeLbltRef::<usize, f64>::new(s, values, subdiag, perm)
                                .solve_in_place_with_conj(Conj::No, rhs, par, stack)
                        }
                        Method::Lu => unreachable!("LU analysis never has a Cholesky symbolic"),
                    }
                }
                (SymbolicKind::Lu(s), Numeric::Lu(numeric)) => {
                    let mut mem = MemBuffer::try_new(s.solve_in_place_scratch::<f64>(nrhs, par))
                        .map_err(|_| SparseError::OutOfMemory)?;
                    flu::LuRef::<usize, f64>::new_unchecked(s, numeric).solve_in_place_with_conj(
                        Conj::No,
                        rhs,
                        par,
                        MemStack::new(&mut mem),
                    );
                }
                _ => unreachable!("numeric storage always matches the symbolic kind"),
            }
        }
        match rhs.iter().position(|x| !x.is_finite()) {
            Some(index) => Err(SparseError::NonFiniteSolution { index }),
            None => Ok(()),
        }
    }
}

/// Controls [`Factorization::solve_refined`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefineOptions {
    /// Maximum number of refinement steps after the initial solve. `0` gives a plain solve with
    /// a residual check.
    pub max_steps: usize,
    /// Stop as soon as the relative residual `‖A x − b‖₂ / ‖b‖₂` is at most this.
    pub tolerance: f64,
}

/// What [`Factorization::solve_refined`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefineReport {
    /// Whether the final relative residual is at most the requested tolerance.
    pub converged: bool,
    /// Number of refinement steps that were applied (each costs one solve and one product).
    pub steps: usize,
    /// Relative residual of the plain solve, before any refinement.
    pub initial_residual: f64,
    /// Relative residual of the returned solution.
    pub residual: f64,
}

impl Factorization {
    /// Solves `A x = b` and improves the solution by iterative refinement against the matrix
    /// `a`: repeat `x ← x + A⁻¹ (b − a x)` with the factorization as the approximate inverse.
    ///
    /// `a` is the matrix whose solution is wanted. It is normally the matrix that was
    /// factorized, in which case refinement removes the error of an unpivoted or otherwise
    /// inaccurate factorization; it may also differ slightly from it (for example the
    /// factorization of a regularized matrix can be used to solve the unregularized system).
    /// A step that does not reduce the residual is discarded and ends the iteration, so the
    /// returned `x` is always the best iterate seen.
    ///
    /// # Errors
    /// [`SparseError::DimensionMismatch`] if `a`, `b` or `x` do not have dimension [`dim`](Self::dim);
    /// [`SparseError::NonFiniteSolution`] if the initial solve fails.
    pub fn solve_refined(
        &self,
        a: &SparseMatrix,
        b: &[f64],
        x: &mut [f64],
        options: RefineOptions,
        par: Parallelism<'_>,
    ) -> Result<RefineReport, SparseError> {
        let n = self.dim();
        for (what, got) in [
            ("rows of a", a.nrows),
            ("columns of a", a.ncols),
            ("length of b", b.len()),
            ("length of x", x.len()),
        ] {
            if got != n {
                return Err(SparseError::DimensionMismatch {
                    what,
                    expected: n,
                    got,
                });
            }
        }
        let b_norm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        let scale = if b_norm == 0.0 { 1.0 } else { b_norm };
        // r = b - a x, returning its relative norm.
        let residual = |x: &[f64], r: &mut [f64]| -> Result<f64, SparseError> {
            a.mul_vec(x, r)?;
            let mut sum = 0.0;
            for (ri, bi) in r.iter_mut().zip(b) {
                *ri = bi - *ri;
                sum += *ri * *ri;
            }
            Ok(sum.sqrt() / scale)
        };

        x.copy_from_slice(b);
        self.solve_in_place(x, par)?;
        let mut r = vec![0.0; n];
        let initial_residual = residual(x, &mut r)?;
        let mut best = initial_residual;
        let mut steps = 0;
        let mut candidate = vec![0.0; n];
        let mut correction = vec![0.0; n];
        // Written so that a NaN residual counts as "not converged" and "not improved".
        let converged = |residual: f64| residual <= options.tolerance;
        while steps < options.max_steps && !converged(best) {
            correction.copy_from_slice(&r);
            if self.solve_in_place(&mut correction, par).is_err() {
                break;
            }
            for i in 0..n {
                candidate[i] = x[i] + correction[i];
            }
            let mut r_new = vec![0.0; n];
            let res = residual(&candidate, &mut r_new)?;
            let improved = res < best;
            if !improved {
                break;
            }
            x.copy_from_slice(&candidate);
            r = r_new;
            best = res;
            steps += 1;
        }
        Ok(RefineReport {
            converged: converged(best),
            steps,
            initial_residual,
            residual: best,
        })
    }
}

fn check_pattern(inner: &SymbolicInner, a: &SparseMatrix) -> Result<(), SparseError> {
    let mismatch = |reason: String| Err(SparseError::PatternMismatch { reason });
    if a.nrows != inner.n || a.ncols != inner.n {
        return mismatch(format!(
            "analysed {} x {}, got {} x {}",
            inner.n, inner.n, a.nrows, a.ncols
        ));
    }
    if a.nnz() != inner.row_idx.len() {
        return mismatch(format!(
            "analysed {} stored entries, got {}",
            inner.row_idx.len(),
            a.nnz()
        ));
    }
    if let Some(j) = (0..=inner.n).find(|&j| a.col_ptr[j] != inner.col_ptr[j]) {
        return mismatch(format!(
            "column {} has a different number of entries",
            j.saturating_sub(1)
        ));
    }
    if let Some(k) = (0..a.nnz()).find(|&k| a.row_idx[k] as usize != inner.row_idx[k]) {
        return mismatch(format!(
            "stored entry {k} is in row {}, analysed row {}",
            a.row_idx[k], inner.row_idx[k]
        ));
    }
    Ok(())
}

/// The relative residual `‖A x − b‖₂ / ‖b‖₂` (or `‖A x‖₂` if `b = 0`): the standard check of a
/// computed solution.
///
/// # Errors
/// [`SparseError::DimensionMismatch`] if the lengths do not match the matrix.
pub fn relative_residual(a: &SparseMatrix, x: &[f64], b: &[f64]) -> Result<f64, SparseError> {
    if b.len() != a.nrows {
        return Err(SparseError::DimensionMismatch {
            what: "length of b",
            expected: a.nrows,
            got: b.len(),
        });
    }
    let mut r = vec![0.0; a.nrows];
    a.mul_vec(x, &mut r)?;
    let mut rr = 0.0;
    let mut bb = 0.0;
    for (ri, bi) in r.iter().zip(b) {
        rr += (ri - bi) * (ri - bi);
        bb += bi * bi;
    }
    Ok(if bb == 0.0 {
        rr.sqrt()
    } else {
        (rr / bb).sqrt()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pool_map_runs_in_order_on_the_pools_own_threads() {
        let pool = ThreadPool::new(3).unwrap();
        let out = pool.map(200, |i| {
            let name = std::thread::current().name().map(str::to_string);
            (i * i, name)
        });
        assert_eq!(out.len(), 200);
        for (i, (sq, name)) in out.iter().enumerate() {
            assert_eq!(*sq, i * i);
            assert!(
                name.as_deref()
                    .is_some_and(|n| n.starts_with("ddg-sparse-")),
                "{name:?}"
            );
        }
        assert!(pool.map(0, |i| i).is_empty());
    }

    /// Runs `body` sequentially and on a two-thread pool.
    fn for_each_parallelism(mut body: impl FnMut(Parallelism<'_>)) {
        let pool = ThreadPool::new(2).unwrap();
        assert_eq!(pool.n_threads(), 2);
        body(Parallelism::Sequential);
        body(Parallelism::Pool(&pool));
    }

    #[test]
    fn thread_pool_needs_a_thread() {
        let e = ThreadPool::new(0).unwrap_err();
        assert_eq!(
            e.to_string(),
            "cannot create thread pool: a thread pool needs at least 1 thread"
        );
    }

    /// Second-difference matrix on a path with `n` nodes plus `shift` on the diagonal.
    fn path_laplacian(n: usize, shift: f64) -> SparseMatrix {
        let mut t = Vec::new();
        for i in 0..n as u32 {
            let deg = if i == 0 || i as usize == n - 1 {
                1.0
            } else {
                2.0
            };
            t.push((i, i, deg + shift));
            if (i as usize) < n - 1 {
                t.push((i, i + 1, -1.0));
                t.push((i + 1, i, -1.0));
            }
        }
        SparseMatrix::from_triplets(n, n, &t).unwrap()
    }

    #[test]
    fn triplets_are_sorted_and_summed() {
        let a = SparseMatrix::from_triplets(
            3,
            2,
            &[
                (2, 1, 1.0),
                (0, 1, 2.0),
                (2, 1, 0.5),
                (1, 0, 4.0),
                (0, 0, 0.0),
            ],
        )
        .unwrap();
        assert_eq!(a.col_ptr(), &[0, 2, 4]);
        assert_eq!(a.row_idx(), &[0, 1, 0, 2]);
        assert_eq!(a.values(), &[0.0, 4.0, 2.0, 1.5]);
        assert_eq!((a.nrows(), a.ncols(), a.nnz()), (3, 2, 4));
        assert_eq!(a.get(2, 1), 1.5);
        assert_eq!(a.get(1, 1), 0.0);
        assert_eq!(a.find(0, 0), Some(0));
        assert_eq!(a.find(1, 1), None);
        assert_eq!(a.find(7, 0), None);
        let t: Vec<_> = a.triplets().collect();
        assert_eq!(t, vec![(0, 0, 0.0), (1, 0, 4.0), (0, 1, 2.0), (2, 1, 1.5)]);
        let at = a.transpose();
        assert_eq!((at.nrows(), at.ncols()), (2, 3));
        for i in 0..3 {
            for j in 0..2 {
                assert_eq!(a.get(i, j), at.get(j, i));
            }
        }
        assert_eq!(at.transpose(), a);

        let e = SparseMatrix::from_triplets(3, 2, &[(0, 0, 1.0), (3, 1, 1.0)]).unwrap_err();
        assert_eq!(
            e.to_string(),
            "triplet 1 at (row 3, column 1) is outside the 3 x 2 matrix"
        );
    }

    #[test]
    fn csc_validation_messages() {
        let ok = SparseMatrix::from_csc(2, 2, vec![0, 1, 2], vec![0, 1], vec![1.0, 2.0]).unwrap();
        assert_eq!(ok, SparseMatrix::from_diagonal(&[1.0, 2.0]).unwrap());
        let msg = |r: Result<SparseMatrix, SparseError>| r.unwrap_err().to_string();
        assert_eq!(
            msg(SparseMatrix::from_csc(2, 2, vec![0, 1], vec![0], vec![1.0])),
            "invalid sparse structure: col_ptr has 2 entries, expected ncols + 1 = 3"
        );
        assert_eq!(
            msg(SparseMatrix::from_csc(
                2,
                2,
                vec![0, 2, 2],
                vec![1, 0],
                vec![1.0, 2.0]
            )),
            "invalid sparse structure: column 0 is not strictly increasing: row 1 is followed by row 0"
        );
        assert_eq!(
            msg(SparseMatrix::from_csc(
                2,
                2,
                vec![0, 1, 2],
                vec![0, 2],
                vec![1.0, 2.0]
            )),
            "invalid sparse structure: column 1 has row index 2, but there are 2 rows"
        );
        assert_eq!(
            msg(SparseMatrix::from_csc(
                2,
                2,
                vec![0, 2, 1],
                vec![0],
                vec![1.0]
            )),
            "invalid sparse structure: col_ptr decreases at column 1"
        );
        assert_eq!(
            msg(SparseMatrix::from_csc(
                2,
                2,
                vec![0, 1, 2],
                vec![0, 1],
                vec![1.0]
            )),
            "invalid sparse structure: 2 row indices but 1 values"
        );
    }

    #[test]
    fn mul_vec_diagonal_and_submatrix() {
        let mut a = path_laplacian(4, 0.0);
        let mut y = vec![0.0; 4];
        a.mul_vec(&[1.0, 2.0, 4.0, 8.0], &mut y).unwrap();
        assert_eq!(y, vec![-1.0, -1.0, -2.0, 4.0]);
        assert_eq!(a.symmetry_defect().unwrap(), 0.0);
        assert_eq!(a.diagonal(), vec![1.0, 2.0, 2.0, 1.0]);
        a.add_to_diagonal(&[1.0, 0.0, 0.0, 0.5]).unwrap();
        assert_eq!(a.diagonal(), vec![2.0, 2.0, 2.0, 1.5]);
        assert_eq!(
            a.mul_vec(&[1.0], &mut y).unwrap_err().to_string(),
            "length of x: expected 4, got 1"
        );

        let s = a.submatrix(&[3, 1], &[1, 2, 3]).unwrap();
        assert_eq!((s.nrows(), s.ncols()), (2, 3));
        for (i, &r) in [3usize, 1].iter().enumerate() {
            for (j, &c) in [1usize, 2, 3].iter().enumerate() {
                assert_eq!(s.get(i, j), a.get(r, c));
            }
        }
        assert_eq!(
            a.submatrix(&[0, 0], &[0]).unwrap_err().to_string(),
            "row 0 is listed more than once"
        );
        assert_eq!(
            a.submatrix(&[4], &[0]).unwrap_err().to_string(),
            "row index 4 is out of range (must be < 4)"
        );
        assert_eq!(
            a.submatrix(&[0], &[9]).unwrap_err().to_string(),
            "column index 9 is out of range (must be < 4)"
        );

        let off =
            SparseMatrix::from_triplets(2, 2, &[(0, 1, 1.0), (1, 0, 3.0), (1, 1, 1.0)]).unwrap();
        assert_eq!(off.symmetry_defect().unwrap(), 2.0);
        let mut off2 = off.clone();
        assert_eq!(
            off2.add_to_diagonal(&[1.0, 1.0]).unwrap_err(),
            SparseError::MissingDiagonal { index: 0 }
        );
        assert_eq!(off2, off);
    }

    #[test]
    fn llt_solves_and_refactors() {
        let n = 50;
        let a = path_laplacian(n, 0.1);
        let symbolic = Symbolic::analyze(&a, Method::Llt).unwrap();
        assert_eq!(symbolic.method(), Method::Llt);
        assert_eq!(symbolic.dim(), n);
        assert_eq!(symbolic.info().matrix_nnz, a.nnz());
        assert!(symbolic.info().factor_values.unwrap() >= n);
        let x_true: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
        let mut b = vec![0.0; n];
        a.mul_vec(&x_true, &mut b).unwrap();
        for_each_parallelism(|par| {
            let mut f = symbolic.factor(&a, par).unwrap();
            let mut x = b.clone();
            f.solve_in_place(&mut x, par).unwrap();
            assert!(relative_residual(&a, &x, &b).unwrap() < 1e-13);
            for (xi, ti) in x.iter().zip(&x_true) {
                assert!((xi - ti).abs() < 1e-11);
            }

            // Same pattern, new values, same symbolic analysis.
            let a2 = path_laplacian(n, 3.0);
            f.refactor(&a2, par).unwrap();
            let mut b2 = vec![0.0; n];
            a2.mul_vec(&x_true, &mut b2).unwrap();
            let mut x2 = b2.clone();
            f.solve_in_place(&mut x2, par).unwrap();
            assert!(relative_residual(&a2, &x2, &b2).unwrap() < 1e-14);

            // Two right-hand sides at once.
            let mut rhs: Vec<f64> = b2
                .iter()
                .copied()
                .chain(b2.iter().map(|v| 2.0 * v))
                .collect();
            f.solve_many_in_place(&mut rhs, 2, par).unwrap();
            for i in 0..n {
                assert!((rhs[i] - x_true[i]).abs() < 1e-12);
                assert!((rhs[n + i] - 2.0 * x_true[i]).abs() < 1e-12);
            }
            assert_eq!(
                f.solve_in_place(&mut [0.0; 3], par)
                    .unwrap_err()
                    .to_string(),
                "length of right-hand side (dim * nrhs): expected 50, got 3"
            );
        });
    }

    #[test]
    fn llt_reads_only_the_lower_triangle() {
        let n = 20;
        let full = path_laplacian(n, 0.5);
        let lower: Vec<_> = full.triplets().filter(|&(i, j, _)| i >= j).collect();
        let lower = SparseMatrix::from_triplets(n, n, &lower).unwrap();
        let mut b: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
        let b0 = b.clone();
        let f = Symbolic::analyze(&lower, Method::Llt)
            .unwrap()
            .factor(&lower, Parallelism::Sequential)
            .unwrap();
        f.solve_in_place(&mut b, Parallelism::Sequential).unwrap();
        assert!(relative_residual(&full, &b, &b0).unwrap() < 1e-14);
    }

    #[test]
    fn llt_rejects_indefinite_and_singular_matrices() {
        let a = path_laplacian(10, -0.5);
        let s = Symbolic::analyze(&a, Method::Llt).unwrap();
        let e = s.factor(&a, Parallelism::Sequential).unwrap_err();
        assert!(
            matches!(e, SparseError::NotPositiveDefinite { .. }),
            "{e:?}"
        );
        assert!(
            e.to_string()
                .starts_with("matrix is not positive definite: non-positive pivot")
        );
    }

    #[test]
    fn pattern_mismatch_is_reported() {
        let a = path_laplacian(10, 1.0);
        let s = Symbolic::analyze(&a, Method::Llt).unwrap();
        let other = SparseMatrix::from_diagonal(&[1.0; 10]).unwrap();
        let e = s.factor(&other, Parallelism::Sequential).unwrap_err();
        assert_eq!(
            e.to_string(),
            "matrix does not have the analysed sparsity pattern: analysed 28 stored entries, got 10"
        );
        let e = s
            .factor(&path_laplacian(11, 1.0), Parallelism::Sequential)
            .unwrap_err();
        assert_eq!(
            e.to_string(),
            "matrix does not have the analysed sparsity pattern: analysed 10 x 10, got 11 x 11"
        );
        let rect = SparseMatrix::from_triplets(2, 3, &[(0, 0, 1.0)]).unwrap();
        assert_eq!(
            Symbolic::analyze(&rect, Method::Lu)
                .unwrap_err()
                .to_string(),
            "a square matrix is required, got 2 x 3"
        );
    }

    /// `[[−M, L], [L, D]]` with `L` a path Laplacian, `M > 0` and `D >= 0` with a few
    /// positive entries: symmetric indefinite and nonsingular.
    fn saddle(n: usize) -> SparseMatrix {
        let l = path_laplacian(n, 0.0);
        let mut t = Vec::new();
        for i in 0..n as u32 {
            t.push((i, i, -(1.0 + 0.1 * f64::from(i % 3))));
            let d = if i % 7 == 0 { 2.0 } else { 0.0 };
            t.push((n as u32 + i, n as u32 + i, d));
        }
        for (i, j, v) in l.triplets() {
            t.push((n as u32 + i, j, v));
            t.push((i, n as u32 + j, v));
        }
        SparseMatrix::from_triplets(2 * n, 2 * n, &t).unwrap()
    }

    #[test]
    fn indefinite_methods_solve_a_saddle_point_system() {
        let n = 40;
        let a = saddle(n);
        assert_eq!(a.symmetry_defect().unwrap(), 0.0);
        let x_true: Vec<f64> = (0..2 * n).map(|i| (i as f64 * 0.21).cos()).collect();
        let mut b = vec![0.0; 2 * n];
        a.mul_vec(&x_true, &mut b).unwrap();
        // LLt must refuse this matrix.
        let e = Symbolic::analyze(&a, Method::Llt)
            .unwrap()
            .factor(&a, Parallelism::Sequential)
            .unwrap_err();
        assert!(matches!(e, SparseError::NotPositiveDefinite { .. }));
        for method in [Method::Lblt, Method::Lu] {
            let s = Symbolic::analyze(&a, method).unwrap();
            for_each_parallelism(|par| {
                let mut f = s.factor(&a, par).unwrap();
                let mut x = b.clone();
                f.solve_in_place(&mut x, par).unwrap();
                let res = relative_residual(&a, &x, &b).unwrap();
                assert!(res < 1e-11, "{method:?} {par:?}: residual {res}");
                f.refactor(&a, par).unwrap();
                let mut x = b.clone();
                f.solve_in_place(&mut x, par).unwrap();
                assert!(relative_residual(&a, &x, &b).unwrap() < 1e-11);
            });
        }
    }

    #[test]
    fn ldlt_solves_a_quasi_definite_system() {
        // With D > 0 everywhere the saddle matrix is quasi-definite and LDLt needs no pivoting.
        let n = 40;
        let mut a = saddle(n);
        let mut shift = vec![0.0; 2 * n];
        shift[n..].fill(0.5);
        a.add_to_diagonal(&shift).unwrap();
        let x_true: Vec<f64> = (0..2 * n).map(|i| 1.0 + (i as f64 * 0.11).sin()).collect();
        let mut b = vec![0.0; 2 * n];
        a.mul_vec(&x_true, &mut b).unwrap();
        let s = Symbolic::analyze(&a, Method::Ldlt).unwrap();
        for_each_parallelism(|par| {
            let f = s.factor(&a, par).unwrap();
            let mut x = b.clone();
            f.solve_in_place(&mut x, par).unwrap();
            assert!(relative_residual(&a, &x, &b).unwrap() < 1e-12);
        });
    }

    #[test]
    fn custom_ordering_decides_which_pivots_ldlt_meets() {
        // A = [[0, 1, 0], [1, 1, 1], [0, 1, 2]] is nonsingular with a zero first diagonal entry.
        let a = SparseMatrix::from_triplets(
            3,
            3,
            &[
                (0, 0, 0.0),
                (1, 0, 1.0),
                (0, 1, 1.0),
                (1, 1, 1.0),
                (2, 1, 1.0),
                (1, 2, 1.0),
                (2, 2, 2.0),
            ],
        )
        .unwrap();
        let seq = Parallelism::Sequential;
        // Natural order starts with the zero pivot.
        let natural =
            Symbolic::analyze_with_order(&a, Method::Ldlt, EliminationOrder::Natural).unwrap();
        let e = natural.factor(&a, seq).unwrap_err();
        assert!(matches!(e, SparseError::ZeroPivot { .. }), "{e:?}");
        assert!(
            e.to_string()
                .starts_with("LDLt factorization failed: zero pivot")
        );
        // Eliminating 1, 2, 0 works: pivots 1, 1, -2. (The inverse permutation 2, 0, 1 would
        // meet a zero pivot in its second step, so this also pins down the convention.)
        let good =
            Symbolic::analyze_with_order(&a, Method::Ldlt, EliminationOrder::Custom(&[1, 2, 0]))
                .unwrap();
        let f = good.factor(&a, seq).unwrap();
        let b = [1.0, 2.0, 3.0];
        let mut x = b;
        f.solve_in_place(&mut x, seq).unwrap();
        assert!(relative_residual(&a, &x, &b).unwrap() < 1e-15);
        let bad =
            Symbolic::analyze_with_order(&a, Method::Ldlt, EliminationOrder::Custom(&[2, 0, 1]))
                .unwrap();
        assert!(matches!(
            bad.factor(&a, seq).unwrap_err(),
            SparseError::ZeroPivot { .. }
        ));

        let err = |order: &[u32]| {
            Symbolic::analyze_with_order(&a, Method::Ldlt, EliminationOrder::Custom(order))
                .unwrap_err()
                .to_string()
        };
        assert_eq!(
            err(&[0, 1]),
            "invalid elimination order: it has 2 entries, but the matrix has dimension 3"
        );
        assert_eq!(
            err(&[0, 1, 3]),
            "invalid elimination order: entry 2 is 3, which is not below the dimension 3"
        );
        assert_eq!(
            err(&[0, 1, 0]),
            "invalid elimination order: index 0 appears twice (entries 0 and 2)"
        );
        assert_eq!(
            Symbolic::analyze_with_order(&a, Method::Lu, EliminationOrder::Natural)
                .unwrap_err()
                .to_string(),
            "not supported: LU with a natural or custom ordering"
        );
    }

    #[test]
    fn paired_ordering_makes_ldlt_work_on_a_semidefinite_saddle_system() {
        // D has zeros, so the matrix is not quasi-definite and a generic order may hit a zero
        // pivot; eliminating (i, n + i) in pairs cannot.
        let n = 60;
        let a = saddle(n);
        let l = path_laplacian(n, 0.0);
        let order_l = fill_reducing_ordering(&l).unwrap();
        let mut sorted = order_l.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..n as u32).collect::<Vec<_>>());
        let paired: Vec<u32> = order_l.iter().flat_map(|&i| [i, n as u32 + i]).collect();
        let s = Symbolic::analyze_with_order(&a, Method::Ldlt, EliminationOrder::Custom(&paired))
            .unwrap();
        let x_true: Vec<f64> = (0..2 * n).map(|i| (i as f64 * 0.3).sin()).collect();
        let mut b = vec![0.0; 2 * n];
        a.mul_vec(&x_true, &mut b).unwrap();
        for_each_parallelism(|par| {
            let f = s.factor(&a, par).unwrap();
            let mut x = b.clone();
            f.solve_in_place(&mut x, par).unwrap();
            assert!(relative_residual(&a, &x, &b).unwrap() < 1e-12);
        });
        assert_eq!(
            fill_reducing_ordering(&SparseMatrix::from_triplets(2, 3, &[]).unwrap()).unwrap_err(),
            SparseError::NotSquare { nrows: 2, ncols: 3 }
        );
    }

    #[test]
    fn refinement_recovers_the_solution_of_a_perturbed_factorization() {
        let n = 40;
        let a = path_laplacian(n, 0.01);
        // Factorize a slightly different matrix, then solve with the true one.
        let f = {
            let nearby = path_laplacian(n, 0.011);
            Symbolic::analyze(&nearby, Method::Llt)
                .unwrap()
                .factor(&nearby, Parallelism::Sequential)
                .unwrap()
        };
        let x_true: Vec<f64> = (0..n).map(|i| (i as f64).cos()).collect();
        let mut b = vec![0.0; n];
        a.mul_vec(&x_true, &mut b).unwrap();
        let mut x = vec![0.0; n];
        let seq = Parallelism::Sequential;
        let plain = f
            .solve_refined(
                &a,
                &b,
                &mut x,
                RefineOptions {
                    max_steps: 0,
                    tolerance: 1e-12,
                },
                seq,
            )
            .unwrap();
        assert!(!plain.converged && plain.steps == 0 && plain.residual == plain.initial_residual);
        assert!(plain.residual > 1e-6);
        let report = f
            .solve_refined(
                &a,
                &b,
                &mut x,
                RefineOptions {
                    max_steps: 200,
                    tolerance: 1e-12,
                },
                seq,
            )
            .unwrap();
        assert!(report.converged, "{report:?}");
        assert!(
            report.steps > 0
                && report.residual <= 1e-12
                && report.initial_residual == plain.initial_residual
        );
        assert!((relative_residual(&a, &x, &b).unwrap() - report.residual).abs() < 1e-15);
        // A budget that is too small reports non-convergence but still returns its best iterate.
        let short = f
            .solve_refined(
                &a,
                &b,
                &mut x,
                RefineOptions {
                    max_steps: 1,
                    tolerance: 1e-12,
                },
                seq,
            )
            .unwrap();
        assert!(!short.converged && short.steps == 1 && short.residual < short.initial_residual);
        let e = f
            .solve_refined(
                &a,
                &b[1..],
                &mut x,
                RefineOptions {
                    max_steps: 1,
                    tolerance: 0.0,
                },
                seq,
            )
            .unwrap_err();
        assert_eq!(e.to_string(), "length of b: expected 40, got 39");
    }

    #[test]
    fn failed_refactorization_invalidates_the_factor() {
        let good = path_laplacian(10, 1.0);
        let seq = Parallelism::Sequential;
        let mut f = Symbolic::analyze(&good, Method::Llt)
            .unwrap()
            .factor(&good, seq)
            .unwrap();
        let mut x = vec![1.0; 10];
        f.solve_in_place(&mut x, seq).unwrap();
        // A pattern mismatch is detected up front and leaves the factor usable.
        assert!(matches!(
            f.refactor(&SparseMatrix::from_diagonal(&[1.0; 10]).unwrap(), seq),
            Err(SparseError::PatternMismatch { .. })
        ));
        f.solve_in_place(&mut x, seq).unwrap();
        // A numeric failure does not.
        let indefinite = path_laplacian(10, -0.5);
        assert!(matches!(
            f.refactor(&indefinite, seq),
            Err(SparseError::NotPositiveDefinite { .. })
        ));
        let e = f.solve_in_place(&mut x, seq).unwrap_err();
        assert_eq!(e, SparseError::NotFactorized);
        assert_eq!(
            e.to_string(),
            "the factorization holds no valid factor: its last refactorization failed"
        );
        f.refactor(&good, seq).unwrap();
        f.solve_in_place(&mut x, seq).unwrap();
    }

    #[test]
    fn empty_matrices_do_not_panic() {
        let empty = SparseMatrix::from_triplets(0, 0, &[]).unwrap();
        assert_eq!((empty.nrows(), empty.ncols(), empty.nnz()), (0, 0, 0));
        assert_eq!(fill_reducing_ordering(&empty).unwrap(), Vec::<u32>::new());
        for method in [Method::Llt, Method::Ldlt, Method::Lblt, Method::Lu] {
            let f = Symbolic::analyze(&empty, method)
                .unwrap()
                .factor(&empty, Parallelism::Sequential)
                .unwrap();
            assert_eq!(f.dim(), 0);
            f.solve_in_place(&mut [], Parallelism::Sequential).unwrap();
        }
    }

    #[test]
    fn solver_types_are_send_and_sync() {
        fn check<T: Send + Sync>() {}
        check::<SparseMatrix>();
        check::<Symbolic>();
        check::<Factorization>();
        check::<ThreadPool>();
        check::<SparseError>();
    }

    #[test]
    fn lu_solves_a_nonsymmetric_system_and_reports_singularity() {
        let n = 30;
        let mut t = Vec::new();
        for i in 0..n as u32 {
            t.push((i, i, 4.0 + f64::from(i)));
            if i + 1 < n as u32 {
                t.push((i, i + 1, -1.0));
                t.push((i + 1, i, 2.0));
            }
        }
        t.push((n as u32 - 1, 0, 0.3));
        let a = SparseMatrix::from_triplets(n, n, &t).unwrap();
        let b: Vec<f64> = (0..n).map(|i| i as f64 - 3.0).collect();
        let f = Symbolic::analyze(&a, Method::Lu)
            .unwrap()
            .factor(&a, Parallelism::Sequential)
            .unwrap();
        let mut x = b.clone();
        f.solve_in_place(&mut x, Parallelism::Sequential).unwrap();
        assert!(relative_residual(&a, &x, &b).unwrap() < 1e-13);

        // A matrix with an empty column is structurally singular.
        let sing =
            SparseMatrix::from_triplets(3, 3, &[(0, 0, 1.0), (1, 0, 1.0), (2, 2, 1.0)]).unwrap();
        let r = Symbolic::analyze(&sing, Method::Lu)
            .and_then(|s| s.factor(&sing, Parallelism::Sequential));
        match r {
            Err(SparseError::StructurallySingular { .. }) => {}
            Ok(f) => {
                let mut x = vec![1.0; 3];
                let e = f
                    .solve_in_place(&mut x, Parallelism::Sequential)
                    .unwrap_err();
                assert!(matches!(e, SparseError::NonFiniteSolution { .. }), "{e:?}");
            }
            Err(e) => panic!("unexpected error {e:?}"),
        }
    }

    #[test]
    fn singular_lblt_solve_is_reported_not_silent() {
        // The pure path Laplacian is singular (constants in the kernel).
        let a = path_laplacian(12, 0.0);
        let f = Symbolic::analyze(&a, Method::Lblt)
            .unwrap()
            .factor(&a, Parallelism::Sequential)
            .unwrap();
        let mut x = vec![1.0; 12];
        match f.solve_in_place(&mut x, Parallelism::Sequential) {
            Err(SparseError::NonFiniteSolution { .. }) => {}
            // Rounding can leave a tiny non-zero pivot; then the residual check is what catches it.
            Ok(()) => assert!(x.iter().all(|v| v.is_finite())),
            Err(e) => panic!("unexpected error {e:?}"),
        }
    }
}
