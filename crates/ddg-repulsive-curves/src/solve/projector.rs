//! Null-space projector (`SPEC.md` §6.3, B.3.3, D4).
//!
//! ```text
//! C† = Cᵀ (C Cᵀ)⁻¹,     P = I − C† C
//! ```
//!
//! for `C` of full row rank: `C P = 0`, `P² = P`, `Pᵀ = P`, `range P = ker C`. The paper prints
//! `P := C C†` (the `k × k` identity) and `C† := (C Cᵀ)⁻¹ Cᵀ` (incompatible shapes); D4.
//!
//! `P` is never formed: it is applied through a sparse Cholesky factorization
//! (`ddg_core::sparse`, `Llt`) of the row-normalized `S̃ = D^{−½} C Cᵀ D^{−½}`,
//! `D = diag(C Cᵀ)`, so that `(C Cᵀ)⁻¹ = D^{−½} S̃⁻¹ D^{−½}`. The normalization makes the
//! rank test independent of the rows' units. Factorization and solves run sequentially
//! (`C Cᵀ` is small and sparse; the result is deterministic).
//!
//! Uses by the multigrid (WP6): Eq. 30 (`y` solves `Pᵀ Ā P y = Pᵀ a` ⇔ `x = P y` solves the
//! saddle system with `b = 0`) and Eq. 33 (`z = C† b`, `y` from the `b = 0` system with right
//! side `Ā z`, then `x = z − y` solves it with general `b`). Both are proved in
//! `lean/proofs` (`rc-projected-system`, `rc-projection-transform`) and tested in
//! `tests/constraints.rs`.

use ddg_core::sparse::{Factorization, Method, Parallelism, SparseError, SparseMatrix, Symbolic};

use super::sparse_error;
use crate::error::RcError;

/// Rows of `C Cᵀ` whose normalized matrix has an eigenvalue below this are rank deficient.
pub const ROW_RANK_TOLERANCE: f64 = 1e-10;

/// Number of inverse-iteration steps for the smallest eigenvalue estimate.
const INVERSE_ITERATIONS: usize = 8;

/// The projector `P = I − Cᵀ (C Cᵀ)⁻¹ C` onto `ker C`, applied matrix-free.
#[derive(Debug)]
pub struct Projector {
    c: SparseMatrix,
    ct: SparseMatrix,
    /// `D^{−½}`.
    scale: Vec<f64>,
    factor: Option<Factorization>,
    min_scaled_eigenvalue: f64,
}

/// Deterministic pseudo-random start vector in `[−1, 1]` (splitmix64).
fn start_vector(k: usize) -> Vec<f64> {
    (0..k as u64)
        .map(|i| {
            let mut z = i
                .wrapping_add(0x9E37_79B9_7F4A_7C15)
                .wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            (z >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        })
        .collect()
}

fn normalize(x: &mut [f64]) -> f64 {
    let n = x.iter().map(|v| v * v).sum::<f64>().sqrt();
    if n > 0.0 {
        for v in x.iter_mut() {
            *v /= n;
        }
    }
    n
}

/// Inverse iteration with `f`: returns the Rayleigh quotient `yᵀ S y / yᵀ y` of the last
/// iterate and the iterate (normalized).
fn inverse_iteration(s: &SparseMatrix, f: &Factorization) -> Result<(f64, Vec<f64>), RcError> {
    let k = s.nrows();
    let mut x = start_vector(k);
    normalize(&mut x);
    for _ in 0..INVERSE_ITERATIONS {
        match f.solve_in_place(&mut x, Parallelism::Sequential) {
            Ok(()) => {}
            Err(_) => {
                // A non-finite solve: the matrix is singular to working precision.
                return Ok((0.0, x));
            }
        }
        if normalize(&mut x) == 0.0 || x.iter().any(|v| !v.is_finite()) {
            return Ok((0.0, x));
        }
    }
    let mut sx = vec![0.0; k];
    s.mul_vec(&x, &mut sx)
        .map_err(|e| sparse_error(e, "C C^T"))?;
    let rq = x.iter().zip(&sx).map(|(a, b)| a * b).sum::<f64>();
    Ok((rq, x))
}

impl Projector {
    /// Factors `C Cᵀ` (rows normalized) and checks the row rank.
    ///
    /// `row_labels[r]` names row `r` in errors (e.g. the row index in
    /// [`crate::constraints::ConstraintSet::rows`]); `None` uses `r` itself.
    ///
    /// # Errors
    /// [`RcError::RedundantConstraints`] if a row is zero or the normalized `C Cᵀ` has an
    /// eigenvalue below [`ROW_RANK_TOLERANCE`]; the rows named are those of a (numerical)
    /// dependency, i.e. the support of a null vector of `C Cᵀ`. [`RcError::InvalidOption`] for
    /// labels of the wrong length or resource errors.
    pub fn new(c: &SparseMatrix, row_labels: Option<&[usize]>) -> Result<Projector, RcError> {
        let k = c.nrows();
        if let Some(l) = row_labels
            && l.len() != k
        {
            return Err(RcError::InvalidOption {
                name: "row_labels",
                value: l.len().to_string(),
                requirement: "one label per row of C",
            });
        }
        let label = |r: usize| row_labels.map_or(r, |l| l[r]);
        let ct = c.transpose();
        if k == 0 {
            return Ok(Projector {
                c: c.clone(),
                ct,
                scale: Vec::new(),
                factor: None,
                min_scaled_eigenvalue: 1.0,
            });
        }
        // C Cᵀ as the sum over columns of C of outer products.
        let (cp, ri, vals) = (c.col_ptr(), c.row_idx(), c.values());
        let mut diag = vec![0.0; k];
        let mut t = Vec::new();
        for j in 0..c.ncols() {
            for p in cp[j]..cp[j + 1] {
                diag[ri[p] as usize] += vals[p] * vals[p];
                for q in cp[j]..cp[j + 1] {
                    t.push((ri[p], ri[q], vals[p] * vals[q]));
                }
            }
        }
        if let Some(r) = diag.iter().position(|&d| d.is_nan() || d <= 0.0) {
            return Err(RcError::RedundantConstraints {
                rows: vec![label(r)],
            });
        }
        let scale: Vec<f64> = diag.iter().map(|d| 1.0 / d.sqrt()).collect();
        for e in &mut t {
            e.2 *= scale[e.0 as usize] * scale[e.1 as usize];
        }
        // Make sure every diagonal entry is in the pattern (it is: d > 0).
        let s = SparseMatrix::from_triplets(k, k, &t).map_err(|e| sparse_error(e, "C C^T"))?;
        let sym = Symbolic::analyze(&s, Method::Llt).map_err(|e| sparse_error(e, "C C^T"))?;
        let redundant = |x: &[f64]| {
            let m = x.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let mut rows: Vec<usize> = (0..k)
                .filter(|&r| x[r].abs() >= 1e-6 * m)
                .map(label)
                .collect();
            if rows.is_empty() {
                rows = (0..k).map(label).collect();
            }
            rows.sort_unstable();
            RcError::RedundantConstraints { rows }
        };
        match sym.factor(&s, Parallelism::Sequential) {
            Ok(f) => {
                let (lmin, x) = inverse_iteration(&s, &f)?;
                if lmin.is_nan() || lmin <= ROW_RANK_TOLERANCE {
                    return Err(redundant(&x));
                }
                Ok(Projector {
                    c: c.clone(),
                    ct,
                    scale,
                    factor: Some(f),
                    min_scaled_eigenvalue: lmin,
                })
            }
            Err(SparseError::NotPositiveDefinite { .. }) => {
                // Singular: find a null vector with a slightly shifted matrix.
                let mut shifted = s.clone();
                shifted
                    .add_to_diagonal(&vec![1e-3 * ROW_RANK_TOLERANCE; k])
                    .map_err(|e| sparse_error(e, "C C^T"))?;
                let f = sym
                    .factor(&shifted, Parallelism::Sequential)
                    .map_err(|e| sparse_error(e, "C C^T"))?;
                let (_, x) = inverse_iteration(&shifted, &f)?;
                Err(redundant(&x))
            }
            Err(e) => Err(sparse_error(e, "C C^T")),
        }
    }

    /// Number of rows `k`.
    pub fn n_rows(&self) -> usize {
        self.c.nrows()
    }

    /// Number of columns `d|V|`.
    pub fn n_dofs(&self) -> usize {
        self.c.ncols()
    }

    /// Estimate (inverse iteration) of the smallest eigenvalue of the row-normalized `C Cᵀ`.
    pub fn min_scaled_eigenvalue(&self) -> f64 {
        self.min_scaled_eigenvalue
    }

    /// `(C Cᵀ)⁻¹ b`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a wrong length; solve failures.
    pub fn solve_gram(&self, b: &[f64]) -> Result<Vec<f64>, RcError> {
        if b.len() != self.n_rows() {
            return Err(RcError::InvalidOption {
                name: "projector right-hand side",
                value: b.len().to_string(),
                requirement: "one entry per constraint row",
            });
        }
        let Some(f) = &self.factor else {
            return Ok(Vec::new());
        };
        let mut y: Vec<f64> = b.iter().zip(&self.scale).map(|(x, s)| x * s).collect();
        f.solve_in_place(&mut y, Parallelism::Sequential)
            .map_err(|e| sparse_error(e, "C C^T"))?;
        for (v, s) in y.iter_mut().zip(&self.scale) {
            *v *= s;
        }
        Ok(y)
    }

    fn mul_c(&self, x: &[f64]) -> Result<Vec<f64>, RcError> {
        let mut y = vec![0.0; self.n_rows()];
        self.c
            .mul_vec(x, &mut y)
            .map_err(|e| sparse_error(e, "C"))?;
        Ok(y)
    }

    fn mul_ct(&self, y: &[f64]) -> Result<Vec<f64>, RcError> {
        let mut x = vec![0.0; self.n_dofs()];
        self.ct
            .mul_vec(y, &mut x)
            .map_err(|e| sparse_error(e, "C^T"))?;
        Ok(x)
    }

    /// `C† b = Cᵀ (C Cᵀ)⁻¹ b`: the smallest (Euclidean) `z` with `C z = b`.
    ///
    /// [formula: rc-projector]
    ///
    /// # Errors
    /// As [`Projector::solve_gram`].
    pub fn pseudo_inverse(&self, b: &[f64]) -> Result<Vec<f64>, RcError> {
        if self.n_rows() == 0 {
            return Ok(vec![0.0; self.n_dofs()]);
        }
        let y = self.solve_gram(b)?;
        self.mul_ct(&y)
    }

    /// `P x = x − C† C x`, the orthogonal projection of `x` onto `ker C`.
    ///
    /// [formula: rc-projector]
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a wrong length; solve failures.
    pub fn project(&self, x: &[f64]) -> Result<Vec<f64>, RcError> {
        if x.len() != self.n_dofs() {
            return Err(RcError::InvalidOption {
                name: "projector argument",
                value: x.len().to_string(),
                requirement: "d|V| entries",
            });
        }
        if self.n_rows() == 0 {
            return Ok(x.to_vec());
        }
        let z = self.pseudo_inverse(&self.mul_c(x)?)?;
        Ok(x.iter().zip(&z).map(|(a, b)| a - b).collect())
    }

    /// `Pᵀ Ā P y`, the operator of Eq. 30 (`P` is symmetric), for a metric product
    /// `apply_metric(v) = Ā v`. With `Ā` positive definite on `ker C`, `y` solves
    /// `Pᵀ Ā P y = Pᵀ a` iff `x = P y` solves the saddle system with `b = 0` (`x` is unique,
    /// `y` is not).
    ///
    /// [formula: rc-projected-system]
    ///
    /// # Errors
    /// As [`Projector::project`].
    pub fn projected_operator(
        &self,
        apply_metric: impl Fn(&[f64]) -> Vec<f64>,
        y: &[f64],
    ) -> Result<Vec<f64>, RcError> {
        let py = self.project(y)?;
        self.project(&apply_metric(&py))
    }

    /// Eq. 33, for a general right side: the `x` of `[Ā Cᵀ; C 0] [x; λ] = [a; b]` from a solver
    /// of the `b = 0` system (`solve_kernel(r)` returns the `x` of `[Ā Cᵀ; C 0] [x; µ] = [r; 0]`,
    /// e.g. Eq. 30 by multigrid): `z = C† b`, `y = solve_kernel(a − Ā z)`, `x = z + y`. With
    /// `a = 0` this is the paper's `x = z − y'`, `y'` the solution for the right side `Ā z`.
    ///
    /// [formula: rc-projection-transform]
    ///
    /// # Errors
    /// As [`Projector::pseudo_inverse`], and the errors of `solve_kernel`.
    pub fn solve_with_kernel_solver(
        &self,
        apply_metric: impl Fn(&[f64]) -> Vec<f64>,
        mut solve_kernel: impl FnMut(&[f64]) -> Result<Vec<f64>, RcError>,
        a: &[f64],
        b: &[f64],
    ) -> Result<Vec<f64>, RcError> {
        if a.len() != self.n_dofs() {
            return Err(RcError::InvalidOption {
                name: "saddle right-hand side",
                value: a.len().to_string(),
                requirement: "d|V| entries",
            });
        }
        let z = self.pseudo_inverse(b)?;
        let az = apply_metric(&z);
        let r: Vec<f64> = a.iter().zip(&az).map(|(p, q)| p - q).collect();
        let y = solve_kernel(&r)?;
        Ok(z.iter().zip(&y).map(|(p, q)| p + q).collect())
    }
}
