//! Translation and rank checks of the constraint set (`SPEC.md` §5.2, D5, D9).
//!
//! 1. **Translation.** `C [t_x t_y t_z]` (two columns in planar mode) must have full column
//!    rank, otherwise `Ā` (whose kernel is the translations, D5) is singular on `ker C`:
//!    [`RcError::TranslationNotFixed`]. For the classic metrics (§4.6), whose kernel is the
//!    per-component translations, the same test runs per component
//!    ([`check_component_translations`], [`RcError::SingularMetric`]).
//! 2. **Row rank.** `C Cᵀ` (rows normalized to unit length) is factored by sparse Cholesky and
//!    its smallest eigenvalue estimated by inverse iteration ([`Projector`]). The known
//!    redundancies of D9 are dropped beforehand by [`ConstraintSet::new`] and reported in
//!    [`RankReport::dropped_rows`]; any other rank deficiency is
//!    [`RcError::RedundantConstraints`] naming the rows of the dependency.
//!
//! The rank decisions use the row-normalized Jacobian, so they do not depend on the units of
//! the rows (lengths, `f`, tangents) or on the scale of the curve.

use ddg_core::sparse::SparseMatrix;

use super::ConstraintSet;
use crate::error::RcError;
use crate::metric::DofLayout;
use crate::network::CurveNetwork;
use crate::solve::projector::Projector;

/// Threshold on the smallest eigenvalue of the Gram matrix of the row-normalized translation
/// images (§5.2 step 1).
pub const TRANSLATION_RANK_TOLERANCE: f64 = 1e-10;

/// Result of [`ConstraintSet::analyze`].
#[derive(Debug, Clone, PartialEq)]
pub struct RankReport {
    /// Rows dropped by D9 (indices into [`ConstraintSet::rows`]).
    pub dropped_rows: Vec<usize>,
    /// Number of active rows.
    pub n_active_rows: usize,
    /// Smallest eigenvalue of the row-normalized `C Cᵀ` (estimate; `1` with no rows).
    pub min_scaled_eigenvalue: f64,
    /// Smallest eigenvalue of the Gram matrix of the row-normalized translation images.
    pub translation_eigenvalue: f64,
}

impl ConstraintSet {
    /// The checks of §5.2 at `positions`: translations fixed, active rows independent.
    ///
    /// # Errors
    /// [`RcError::TranslationNotFixed`], [`RcError::RedundantConstraints`],
    /// [`RcError::PositionCountMismatch`].
    pub fn analyze(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<RankReport, RcError> {
        let c = self.jacobian(network, positions)?;
        let translation_eigenvalue = check_translations(&c, self.layout())?;
        let p = Projector::new(&c, Some(self.active_rows()))?;
        Ok(RankReport {
            dropped_rows: self.dropped_rows().to_vec(),
            n_active_rows: self.n_active_rows(),
            min_scaled_eigenvalue: p.min_scaled_eigenvalue(),
            translation_eigenvalue,
        })
    }
}

/// Eigenvalues (increasing) of a symmetric `m × m` row-major matrix, sequential.
pub(crate) fn symmetric_eigenvalues(m: usize, data: &[f64]) -> Result<Vec<f64>, RcError> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::evd::{ComputeEigenvectors, self_adjoint_evd, self_adjoint_evd_scratch};
    if m == 0 {
        return Ok(Vec::new());
    }
    let a = faer::Mat::<f64>::from_fn(m, m, |i, j| data[i * m + j]);
    let mut s = faer::diag::Diag::<f64>::zeros(m);
    let par = faer::Par::Seq;
    let mut mem = MemBuffer::new(self_adjoint_evd_scratch::<f64>(
        m,
        ComputeEigenvectors::No,
        par,
        Default::default(),
    ));
    self_adjoint_evd(
        a.as_ref(),
        s.as_mut(),
        None,
        par,
        MemStack::new(&mut mem),
        Default::default(),
    )
    .map_err(|_| RcError::InvalidOption {
        name: "eigenvalue solver",
        value: format!("{m} x {m} matrix"),
        requirement: "the symmetric eigenvalue iteration did not converge",
    })?;
    let mut ev: Vec<f64> = (0..m).map(|i| s[i]).collect();
    ev.sort_by(f64::total_cmp);
    Ok(ev)
}

/// Smallest eigenvalue of `G = M̃ᵀ M̃`, `M̃[r][(g, s)] = Σ_{v ∈ g} C[r, (v, s)] / ‖C_r‖`, where
/// `group[v]` is the group of vertex `v` (`n_groups` groups).
fn translation_gram_min(
    c: &SparseMatrix,
    layout: DofLayout,
    group: &[usize],
    n_groups: usize,
) -> Result<f64, RcError> {
    let d = layout.dim();
    let m = n_groups * d;
    let ct = c.transpose(); // column r = row r of C
    let (cp, ri, vals) = (ct.col_ptr(), ct.row_idx(), ct.values());
    let mut g = vec![0.0; m * m];
    let mut acc: Vec<(usize, f64)> = Vec::new();
    for r in 0..c.nrows() {
        let range = cp[r]..cp[r + 1];
        let norm2: f64 = vals[range.clone()].iter().map(|x| x * x).sum();
        if norm2 == 0.0 {
            continue;
        }
        let inv = 1.0 / norm2.sqrt();
        acc.clear();
        for p in range {
            let j = ri[p] as usize;
            let (v, s) = (j / d, j % d);
            acc.push((group[v] * d + s, vals[p] * inv));
        }
        acc.sort_by_key(|&(k, _)| k);
        let mut merged: Vec<(usize, f64)> = Vec::with_capacity(acc.len());
        for &(k, x) in &acc {
            match merged.last_mut() {
                Some(last) if last.0 == k => last.1 += x,
                _ => merged.push((k, x)),
            }
        }
        for &(a, x) in &merged {
            for &(b, y) in &merged {
                g[a * m + b] += x * y;
            }
        }
    }
    let ev = symmetric_eigenvalues(m, &g)?;
    Ok(ev.first().copied().unwrap_or(f64::INFINITY))
}

/// §5.2 step 1: `C [t_x t_y t_z]` has full column rank (`d` columns). Returns the smallest
/// eigenvalue of the Gram matrix of the row-normalized images.
///
/// # Errors
/// [`RcError::TranslationNotFixed`] if it is `≤` [`TRANSLATION_RANK_TOLERANCE`].
pub fn check_translations(c: &SparseMatrix, layout: DofLayout) -> Result<f64, RcError> {
    let n = c.ncols() / layout.dim();
    let lmin = translation_gram_min(c, layout, &vec![0; n], 1)?;
    if lmin.is_nan() || lmin <= TRANSLATION_RANK_TOLERANCE {
        return Err(RcError::TranslationNotFixed);
    }
    Ok(lmin)
}

/// The translation test per connected component, for metrics whose kernel is the
/// per-component translations (L-type metrics of §4.6): `C` must be injective on them.
///
/// # Errors
/// [`RcError::SingularMetric`] (with the number of components) otherwise.
pub fn check_component_translations(
    c: &SparseMatrix,
    network: &CurveNetwork,
    layout: DofLayout,
) -> Result<f64, RcError> {
    let n = network.n_vertices();
    let group: Vec<usize> = (0..n as u32).map(|v| network.component_of(v)).collect();
    let nc = network.components().len();
    let lmin = translation_gram_min(c, layout, &group, nc)?;
    if lmin.is_nan() || lmin <= TRANSLATION_RANK_TOLERANCE {
        return Err(RcError::SingularMetric { components: nc });
    }
    Ok(lmin)
}
