//! Sparse L², H¹, H² metrics (`SPEC.md` §4.6, Supp. 2.1, D23), used by the evaluation of
//! Sec. 7 only.
//!
//! ```text
//! L²:  M = diag(m_v),  m_v = ½ Σ_{I ∋ v} ℓ_I                    (lumped mass; reference: Σℓ/deg)
//! H¹:  L_ij = −1/ℓ_ij for every edge (i, j),  L_ii = Σ_{I ∋ i} 1/ℓ_I
//! H²:  L M⁻¹ L
//! ```
//!
//! `ker L = ker L M⁻¹ L` are the per-component constants, so with several components the
//! constraints must fix every component's translation
//! ([`crate::constraints::rank::check_component_translations`]).
//!
//! [formula: rc-classic-metrics]

use ddg_core::sparse::SparseMatrix;

use super::Metric;
use crate::error::RcError;
use crate::network::CurveNetwork;
use crate::solve::sparse_error;

/// `M = diag(m_v)` with `m_v = ½ Σ_{I ∋ v} ℓ_I`.
///
/// [formula: rc-classic-metrics]
///
/// # Errors
/// Only if the matrix is too large for `u32` indices.
pub fn mass_matrix(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
) -> Result<SparseMatrix, RcError> {
    let m = network.dual_lengths(positions);
    SparseMatrix::from_diagonal(&m).map_err(|e| sparse_error(e, "mass matrix"))
}

/// The weak curve Laplacian `L` with weights `1/ℓ_ij` (positive semidefinite).
///
/// [formula: rc-classic-metrics]
///
/// # Errors
/// Only if the matrix is too large for `u32` indices.
pub fn curve_laplacian(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
) -> Result<SparseMatrix, RcError> {
    let n = network.n_vertices();
    let geom = network.edge_geometry(positions);
    let mut t = Vec::with_capacity(4 * network.n_edges());
    for (i, &[a, b]) in network.edges().iter().enumerate() {
        let w = 1.0 / geom.length[i];
        t.push((a, a, w));
        t.push((b, b, w));
        t.push((a, b, -w));
        t.push((b, a, -w));
    }
    SparseMatrix::from_triplets(n, n, &t).map_err(|e| sparse_error(e, "curve Laplacian"))
}

/// `L M⁻¹ L`: `(L M⁻¹ L)_uv = Σ_w L_uw L_wv / m_w`, summed over `w` in increasing order.
///
/// [formula: rc-classic-metrics]
///
/// # Errors
/// Only if the matrix is too large for `u32` indices.
pub fn bilaplacian(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
) -> Result<SparseMatrix, RcError> {
    let n = network.n_vertices();
    let l = curve_laplacian(network, positions)?;
    let mass = network.dual_lengths(positions);
    let (cp, ri, vals) = (l.col_ptr(), l.row_idx(), l.values());
    let mut t = Vec::new();
    // L is symmetric: column w lists the entries L_uw.
    for w in 0..n {
        let range = cp[w]..cp[w + 1];
        for p in range.clone() {
            for q in range.clone() {
                t.push((ri[p], ri[q], vals[p] * vals[q] / mass[w]));
            }
        }
    }
    SparseMatrix::from_triplets(n, n, &t).map_err(|e| sparse_error(e, "bilaplacian"))
}

/// The metric matrix of `metric` (`L2`, `H1`, `H2`).
///
/// # Errors
/// [`RcError::InvalidOption`] for [`Metric::Sobolev`] (dense, see
/// [`super::sobolev_dense::sobolev_gram`]).
pub fn classic_metric(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    metric: &Metric,
) -> Result<SparseMatrix, RcError> {
    match metric {
        Metric::L2 => mass_matrix(network, positions),
        Metric::H1 => curve_laplacian(network, positions),
        Metric::H2 => bilaplacian(network, positions),
        Metric::Sobolev { .. } => Err(RcError::InvalidOption {
            name: "metric",
            value: "Sobolev".into(),
            requirement: "classic_metric takes L2, H1 or H2; the Sobolev metric is dense",
        }),
    }
}
