//! Total length and length difference penalties (`SPEC.md` §5.3, D20).
//!
//! Both weights may be negative (the reference's scene format allows it).

use std::sync::Arc;

use super::{check_args, check_weight};
use crate::error::RcError;
use crate::kernel::v3::{norm, sub};
use crate::network::CurveNetwork;
use crate::step::PenaltyTerm;
use crate::subdivide::SubdivisionMap;

/// `w Σ_I ℓ_I`; gradient `i₂ += w T_I`, `i₁ −= w T_I`.
///
/// The paper's penalty is the total length itself (the reference uses `½ w L²`, D20). Splitting
/// edges at their midpoints keeps the value.
///
/// [formula: rc-penalty-total-length]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TotalLengthPenalty {
    weight: f64,
}

impl TotalLengthPenalty {
    /// The penalty with weight `w` (any finite sign).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-finite weight.
    pub fn new(weight: f64) -> Result<Self, RcError> {
        check_weight("TotalLength.weight", weight, true)?;
        Ok(TotalLengthPenalty { weight })
    }

    /// The weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
}

impl PenaltyTerm for TotalLengthPenalty {
    fn name(&self) -> &str {
        "TotalLength"
    }

    /// [formula: rc-penalty-total-length]
    fn evaluate(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        mut gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<f64, RcError> {
        check_args(network, positions, &gradient)?;
        let mut total = 0.0;
        for &[a, b] in network.edges() {
            let d = sub(positions[b as usize], positions[a as usize]);
            let l = norm(d);
            total += l;
            if let Some(g) = gradient.as_deref_mut() {
                for k in 0..3 {
                    let t = self.weight * d[k] / l;
                    g[b as usize][k] += t;
                    g[a as usize][k] -= t;
                }
            }
        }
        Ok(self.weight * total)
    }

    fn subdivided(&self, _map: &SubdivisionMap) -> Result<Arc<dyn PenaltyTerm>, RcError> {
        Ok(Arc::new(*self))
    }
}

/// `w Σ_{v : deg v = 2} (ℓ_{I_v} − ℓ_{J_v})²` with `I_v, J_v` the two edges at `v`; gradient
/// `2w (ℓ_I − ℓ_J)(∇ℓ_I − ∇ℓ_J)` (the paper's formula; the reference has an extra ½, D20). It
/// encourages equal lengths of adjacent edges without a hard constraint. Vertices of other
/// degrees contribute nothing.
///
/// [formula: rc-penalty-length-difference]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LengthDifference {
    weight: f64,
}

impl LengthDifference {
    /// The penalty with weight `w` (any finite sign).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-finite weight.
    pub fn new(weight: f64) -> Result<Self, RcError> {
        check_weight("LengthDifference.weight", weight, true)?;
        Ok(LengthDifference { weight })
    }

    /// The weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
}

impl PenaltyTerm for LengthDifference {
    fn name(&self) -> &str {
        "LengthDifference"
    }

    /// [formula: rc-penalty-length-difference]
    fn evaluate(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        mut gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<f64, RcError> {
        check_args(network, positions, &gradient)?;
        let edges = network.edges();
        let geo = network.edge_geometry(positions);
        let mut total = 0.0;
        for v in 0..network.n_vertices() as u32 {
            let inc = network.incident_edges(v);
            if inc.len() != 2 {
                continue;
            }
            let (i, j) = (inc[0] as usize, inc[1] as usize);
            let diff = geo.length[i] - geo.length[j];
            total += diff * diff;
            if let Some(g) = gradient.as_deref_mut() {
                let f = 2.0 * self.weight * diff;
                // ∇ℓ_I: i2 += T_I, i1 -= T_I; the J term enters with the opposite sign
                for (e, s) in [(i, f), (j, -f)] {
                    let [a, b] = edges[e];
                    for k in 0..3 {
                        let t = s * geo.tangent[e][k];
                        g[b as usize][k] += t;
                        g[a as usize][k] -= t;
                    }
                }
            }
        }
        Ok(self.weight * total)
    }

    fn subdivided(&self, _map: &SubdivisionMap) -> Result<Arc<dyn PenaltyTerm>, RcError> {
        Ok(Arc::new(*self))
    }
}
