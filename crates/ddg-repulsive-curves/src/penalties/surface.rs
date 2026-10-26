//! Surface potential penalty (`SPEC.md` §5.3, D19).

use std::sync::Arc;

use rayon::prelude::*;

use super::{check_args, check_weight};
use crate::ccd::CollisionObstacle;
use crate::error::RcError;
use crate::kernel::Exponents;
use crate::kernel::v3::{norm, sub};
use crate::network::CurveNetwork;
use crate::obstacle::Obstacle;
use crate::step::PenaltyTerm;
use crate::subdivide::SubdivisionMap;

/// The default far-field opening `θ_M` of the lumping (`r_N <= θ_M |γ_v − c_N|`).
pub const DEFAULT_THETA: f64 = 0.1;

/// `Ê_M = w Σ_v m_v Σ_t A_t δ(γ_v, t)^{−p}` for a rigid triangle-soup [`Obstacle`]: a
/// Coulomb-like potential of the same order as the tangent-point energy (`p = β − α` by
/// default), with `m_v = ½ Σ_{I ∋ v} ℓ_I` the dual length of the curve vertex, `A_t` the
/// triangle area, and `δ` the point–triangle distance, so that the barrier sits on the faces
/// and not on the mesh vertices (D19). The near field is evaluated exactly, the far field by
/// centroid lumping of BVH nodes ([`Obstacle::potential_at`], opening `theta`).
///
/// The gradient is `w (m_v ∇φ(γ_v) + ½ (φ_{a} + φ_{b}) ± T_I)` per vertex and edge, with
/// `φ(x) = Σ_t A_t δ^{−p}` and `∇δ^{−p} = −p δ^{−p−2}(x − π_t(x))`. The reference's potential
/// is unweighted by the curve masses and its value omits the node mass that its force
/// includes (D19, D20); here value and gradient are one function.
///
/// If a vertex lies on the obstacle the value is `+∞` (which every line search rejects) and a
/// gradient request is an error. The obstacle is also given to the collision detection
/// ([`PenaltyTerm::obstacle`]), so that guarded steps never cross it.
///
/// Deterministic: per-vertex terms are independent and summed in vertex order, whatever the
/// size of the rayon pool in which `evaluate` runs.
///
/// [formula: rc-penalty-surface]
#[derive(Debug, Clone)]
pub struct SurfacePotential {
    obstacle: Arc<Obstacle>,
    weight: f64,
    p: f64,
    theta: f64,
}

impl SurfacePotential {
    /// The potential with weight `w >= 0`, exponent `p > 0` and far-field opening `theta >= 0`
    /// (`0`: exact, no lumping).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a negative or non-finite weight or `theta`, or a
    /// non-positive or non-finite `p`.
    pub fn new(obstacle: Arc<Obstacle>, weight: f64, p: f64, theta: f64) -> Result<Self, RcError> {
        check_weight("SurfacePotential.weight", weight, false)?;
        if !(p.is_finite() && p > 0.0) {
            return Err(RcError::InvalidOption {
                name: "SurfacePotential.p",
                value: p.to_string(),
                requirement: "a finite positive exponent",
            });
        }
        if !(theta.is_finite() && theta >= 0.0) {
            return Err(RcError::InvalidOption {
                name: "SurfacePotential.theta",
                value: theta.to_string(),
                requirement: "a finite non-negative opening",
            });
        }
        Ok(SurfacePotential {
            obstacle,
            weight,
            p,
            theta,
        })
    }

    /// The potential with `p = β − α` ("a Coulomb potential of the same order") and
    /// [`DEFAULT_THETA`].
    ///
    /// # Errors
    /// As [`SurfacePotential::new`].
    pub fn with_exponents(
        obstacle: Arc<Obstacle>,
        weight: f64,
        exponents: &Exponents,
    ) -> Result<Self, RcError> {
        SurfacePotential::new(
            obstacle,
            weight,
            exponents.beta() - exponents.alpha(),
            DEFAULT_THETA,
        )
    }

    /// The obstacle.
    pub fn obstacle_mesh(&self) -> &Arc<Obstacle> {
        &self.obstacle
    }

    /// The weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }

    /// The exponent `p`.
    pub fn exponent(&self) -> f64 {
        self.p
    }

    /// The far-field opening.
    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// The same potential with a different opening (`0`: exact).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a negative or non-finite `theta`.
    pub fn with_theta(&self, theta: f64) -> Result<Self, RcError> {
        SurfacePotential::new(self.obstacle.clone(), self.weight, self.p, theta)
    }

    /// Value and optionally gradient with all triangles exact (the reference for tests and the
    /// brute-force timing).
    ///
    /// # Errors
    /// As [`PenaltyTerm::evaluate`].
    pub fn evaluate_brute(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<f64, RcError> {
        self.evaluate_with(network, positions, gradient, true)
    }

    fn evaluate_with(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        gradient: Option<&mut [[f64; 3]]>,
        brute: bool,
    ) -> Result<f64, RcError> {
        check_args(network, positions, &gradient)?;
        let want = gradient.is_some();
        let obs = &*self.obstacle;
        let (p, theta) = (self.p, self.theta);
        let q: Vec<_> = positions
            .par_iter()
            .map(|&x| {
                if brute {
                    obs.potential_at_brute(x, p, want)
                } else {
                    obs.potential_at(x, p, theta, want)
                }
            })
            .collect();
        if want && let Some(v) = q.iter().position(|r| r.value.is_infinite()) {
            return Err(RcError::InvalidOption {
                name: "SurfacePotential.positions",
                value: format!("vertex {v} lies on an obstacle triangle"),
                requirement: "curve vertices off the obstacle (the potential is infinite there)",
            });
        }
        let m = network.dual_lengths(positions);
        let mut total = 0.0;
        for (v, r) in q.iter().enumerate() {
            total += m[v] * r.value;
        }
        if let Some(g) = gradient {
            for (v, r) in q.iter().enumerate() {
                for k in 0..3 {
                    g[v][k] += self.weight * m[v] * r.gradient[k];
                }
            }
            // ∂m_v/∂γ: each edge contributes ℓ_I (φ_a + φ_b)/2
            for &[a, b] in network.edges() {
                let d = sub(positions[b as usize], positions[a as usize]);
                let l = norm(d);
                let s = 0.5 * self.weight * (q[a as usize].value + q[b as usize].value) / l;
                for k in 0..3 {
                    g[b as usize][k] += s * d[k];
                    g[a as usize][k] -= s * d[k];
                }
            }
        }
        Ok(self.weight * total)
    }
}

impl PenaltyTerm for SurfacePotential {
    fn name(&self) -> &str {
        "SurfacePotential"
    }

    /// [formula: rc-penalty-surface]
    fn evaluate(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<f64, RcError> {
        self.evaluate_with(network, positions, gradient, false)
    }

    fn obstacle(&self) -> Option<Arc<dyn CollisionObstacle>> {
        Some(self.obstacle.clone())
    }

    fn subdivided(&self, _map: &SubdivisionMap) -> Result<Arc<dyn PenaltyTerm>, RcError> {
        Ok(Arc::new(self.clone()))
    }
}
