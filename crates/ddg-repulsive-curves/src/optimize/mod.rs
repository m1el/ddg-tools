//! Soft-constraint optimizers of the Sec. 7 evaluation (`SPEC.md` §6.4, D21, D22).
//!
//! The objective is `F = Ê + Hencky (+ other penalties)` without hard constraints. Translations
//! are removed from the direction only (Supp. 2.2: the gradient is taken "orthogonal (in the L²
//! sense) to the space of constant vector fields"): the preconditioned gradient `z` solves Eq. 22
//! with `C` replaced by the rows `M_c = [m₁I₃ … m_|V|I₃]` ([`mean_free_rows`]), with the metric of
//! the current curve at every iteration ("dynamic", Supp. 3.3).
//!
//! * [`lbfgs`]: L-BFGS with memory 30 and the dynamic metric as the initial inverse Hessian.
//! * [`ncg`]: Polak–Ribière (PR+) nonlinear conjugate gradients with the same preconditioner.
//!
//! The session (`step`) runs both with the line searches of §7. This module also holds the
//! Hencky strain penalty that the evaluation needs ([`Hencky`], §5.3, D21); work package WP7
//! owns the other penalties.

use std::sync::Arc;

use ddg_core::sparse::SparseMatrix;

use crate::error::RcError;
use crate::kernel::v3::{norm, sub};
use crate::metric::DofLayout;
use crate::network::CurveNetwork;
use crate::solve::sparse_error;
use crate::step::PenaltyTerm;
use crate::subdivide::SubdivisionMap;

pub mod lbfgs;
pub mod ncg;

pub use lbfgs::Lbfgs;
pub use ncg::Ncg;

/// The descent method (`SPEC.md` §9.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Method {
    /// Projected gradient in the chosen metric (§6.1), with hard constraints.
    #[default]
    ProjectedGradient,
    /// L-BFGS on a soft problem (no hard constraints), memory `memory` (Supp. 3.3: 30).
    Lbfgs {
        /// Number of stored pairs.
        memory: usize,
    },
    /// Polak–Ribière nonlinear CG on a soft problem.
    Ncg,
}

impl Method {
    /// `true` for the soft-problem methods.
    pub fn is_soft(&self) -> bool {
        !matches!(self, Method::ProjectedGradient)
    }

    /// Checks the parameters.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for an L-BFGS memory of 0.
    pub fn validate(&self) -> Result<(), RcError> {
        if let Method::Lbfgs { memory: 0 } = self {
            return Err(RcError::InvalidOption {
                name: "method.memory",
                value: "0".into(),
                requirement: "at least one stored pair",
            });
        }
        Ok(())
    }
}

/// The rows `M_c = [m₁ I_d … m_|V| I_d]` (`d` rows, interleaved columns, held axis dropped) that
/// remove the translations from the direction of the soft problem: `M_c z = 0` means
/// `Σ_v m_v z_v = 0` (the direction is L²-orthogonal to the constant fields).
///
/// [formula: rc-mean-free-direction]
///
/// # Errors
/// Resource errors of the sparse matrix.
pub fn mean_free_rows(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    layout: DofLayout,
) -> Result<SparseMatrix, RcError> {
    let m = network.dual_lengths(positions);
    let d = layout.dim();
    let mut t = Vec::with_capacity(d * m.len());
    for (v, &mv) in m.iter().enumerate() {
        for s in 0..d {
            t.push((s as u32, layout.index(v, s) as u32, mv));
        }
    }
    SparseMatrix::from_triplets(d, layout.n_dofs(m.len()), &t)
        .map_err(|e| sparse_error(e, "mean-free rows"))
}

/// The Hencky strain penalty `w Σ_I ℓ⁰_I ln²(ℓ_I/ℓ⁰_I)` (§5.3; Supp. 2.2 "Hencky strain terms …
/// for each edge, weighted by edge lengths of the initial configuration"), the soft replacement
/// of `EdgeLengths` in the evaluation (D21). `∂/∂ℓ_I = 2w (ℓ⁰_I/ℓ_I) ln(ℓ_I/ℓ⁰_I)`, and
/// `∂ℓ_I/∂γ_{i₂} = T_I = −∂ℓ_I/∂γ_{i₁}`.
///
/// [formula: rc-penalty-hencky]
#[derive(Debug, Clone, PartialEq)]
pub struct Hencky {
    weight: f64,
    rest: Vec<f64>,
}

impl Hencky {
    /// The penalty with rest lengths `rest` (one per edge).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a negative or non-finite weight or non-positive rest
    /// lengths.
    pub fn new(weight: f64, rest: Vec<f64>) -> Result<Hencky, RcError> {
        if !(weight.is_finite() && weight >= 0.0) {
            return Err(RcError::InvalidOption {
                name: "Hencky.weight",
                value: weight.to_string(),
                requirement: "a finite non-negative weight",
            });
        }
        if let Some(x) = rest.iter().find(|x| !(x.is_finite() && **x > 0.0)) {
            return Err(RcError::InvalidOption {
                name: "Hencky.rest_lengths",
                value: x.to_string(),
                requirement: "finite positive lengths",
            });
        }
        Ok(Hencky { weight, rest })
    }

    /// The penalty with the current edge lengths as rest lengths.
    ///
    /// # Errors
    /// As [`Hencky::new`]; [`RcError::PositionCountMismatch`].
    pub fn at_rest(
        weight: f64,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
    ) -> Result<Hencky, RcError> {
        if positions.len() != network.n_vertices() {
            return Err(RcError::PositionCountMismatch {
                expected: network.n_vertices(),
                found: positions.len(),
            });
        }
        let rest = network
            .edges()
            .iter()
            .map(|&[a, b]| norm(sub(positions[b as usize], positions[a as usize])))
            .collect();
        Hencky::new(weight, rest)
    }

    /// The weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }

    /// The rest lengths.
    pub fn rest_lengths(&self) -> &[f64] {
        &self.rest
    }

    /// The largest `|ln(ℓ_I/ℓ⁰_I)|` (the calibration measure of §14.2).
    pub fn max_log_strain(&self, network: &CurveNetwork, positions: &[[f64; 3]]) -> f64 {
        network
            .edges()
            .iter()
            .zip(&self.rest)
            .map(|(&[a, b], r)| {
                (norm(sub(positions[b as usize], positions[a as usize])) / r)
                    .ln()
                    .abs()
            })
            .fold(0.0, f64::max)
    }
}

impl PenaltyTerm for Hencky {
    fn name(&self) -> &str {
        "Hencky"
    }

    fn check(&self, network: &CurveNetwork) -> Result<(), RcError> {
        if self.rest.len() != network.n_edges() {
            return Err(RcError::InvalidOption {
                name: "Hencky.rest_lengths",
                value: format!("{} lengths", self.rest.len()),
                requirement: "one rest length per edge",
            });
        }
        Ok(())
    }

    /// [formula: rc-penalty-hencky]
    fn evaluate(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        mut gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<f64, RcError> {
        self.check(network)?;
        let mut value = 0.0;
        for (&[a, b], &r) in network.edges().iter().zip(&self.rest) {
            let d = sub(positions[b as usize], positions[a as usize]);
            let l = norm(d);
            let s = (l / r).ln();
            value += self.weight * r * s * s;
            if let Some(g) = gradient.as_deref_mut() {
                let dl = 2.0 * self.weight * (r / l) * s;
                for k in 0..3 {
                    let t = dl * d[k] / l;
                    g[b as usize][k] += t;
                    g[a as usize][k] -= t;
                }
            }
        }
        Ok(value)
    }

    fn subdivided(&self, map: &SubdivisionMap) -> Result<Arc<dyn PenaltyTerm>, RcError> {
        Ok(Arc::new(Hencky {
            weight: self.weight,
            rest: map.map_edge_values(&self.rest, |x| 0.5 * x),
        }))
    }
}

/// `Σ a_i b_i`, sequential (fixed order).
pub(crate) fn dotv(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
