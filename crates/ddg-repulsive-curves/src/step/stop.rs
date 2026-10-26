//! Stop rules and status (`SPEC.md` §7.5, D14).
//!
//! * [`StopRule::GradientNorm`]: `‖g̃‖ < eps`, by default the Euclidean 2-norm over all
//!   coordinates with `eps = 1e-4` (paper Sec. 5.4 and Supp. 1.3, D14), or the mass-weighted
//!   `√(Σ_v m_v |g̃_v|²)` ([`GradNorm::MassL2`]). Both are absolute and depend on the scale of the
//!   curve (§4.5), as in the paper.
//! * [`StopRule::RelativeGradient`]: `‖g̃_k‖ / ‖g̃_0‖ < ratio` (scale free, for Blender).
//! * [`StopRule::EnergyStall`]: the objective decreased by less than `relative · |F|` over the
//!   last `window` accepted steps (scale free).
//!
//! `max_iterations` and `time_limit` are options of the session; a stuck line search ends a run
//! with [`Status::Stuck`] unless a length schedule is still moving the targets.

use std::collections::VecDeque;

use crate::error::RcError;

/// Which norm of the projected gradient `g̃` the stop rule measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GradNorm {
    /// `‖g̃‖₂` over all coordinates (D14, the default).
    #[default]
    Euclidean,
    /// `√(Σ_v m_v |g̃_v|²)` with the dual lengths `m_v`.
    MassL2,
}

/// When a run has converged.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StopRule {
    /// `‖g̃‖ < eps` in the chosen norm.
    GradientNorm {
        /// The norm.
        norm: GradNorm,
        /// The threshold.
        eps: f64,
    },
    /// `‖g̃_k‖₂ / ‖g̃_0‖₂ < ratio` (`g̃_0` of the first step of the session).
    RelativeGradient {
        /// The ratio.
        ratio: f64,
    },
    /// `F_{k−window} − F_k ≤ relative · |F_k|` over the last `window` accepted steps.
    EnergyStall {
        /// The relative decrease.
        relative: f64,
        /// Number of accepted steps.
        window: usize,
    },
}

impl Default for StopRule {
    /// `GradientNorm { Euclidean, 1e-4 }` (`SPEC.md` §10, D14).
    fn default() -> Self {
        StopRule::GradientNorm {
            norm: GradNorm::Euclidean,
            eps: 1e-4,
        }
    }
}

impl StopRule {
    /// Checks the parameters.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for non-positive or non-finite thresholds or a zero window.
    pub fn validate(&self) -> Result<(), RcError> {
        let (name, v) = match *self {
            StopRule::GradientNorm { eps, .. } => ("stop.eps", eps),
            StopRule::RelativeGradient { ratio } => ("stop.ratio", ratio),
            StopRule::EnergyStall { relative, window } => {
                if window == 0 {
                    return Err(RcError::InvalidOption {
                        name: "stop.window",
                        value: "0".into(),
                        requirement: "at least one step",
                    });
                }
                ("stop.relative", relative)
            }
        };
        if !(v.is_finite() && v > 0.0) {
            return Err(RcError::InvalidOption {
                name,
                value: v.to_string(),
                requirement: "a finite positive threshold",
            });
        }
        Ok(())
    }
}

/// The state of a session or a run.
#[derive(Debug, Clone, PartialEq)]
pub enum Status {
    /// Not finished.
    Running,
    /// A stop rule holds.
    Converged,
    /// `max_iterations` steps were taken.
    MaxIterations,
    /// The run's time limit was reached.
    TimeLimit,
    /// The line search found no acceptable step (and no schedule is moving the targets).
    Stuck,
    /// Cancelled through the cancel flag or the progress callback.
    Cancelled,
    /// A step failed with an error; the positions are those of the last completed step.
    Error(RcError),
}

impl Status {
    /// `true` for every status but [`Status::Running`].
    pub fn is_terminal(&self) -> bool {
        !matches!(self, Status::Running)
    }
}

/// The outcome of one step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepStatus {
    /// A step was accepted (the curve moved).
    Accepted,
    /// The line search found no acceptable step; the curve did not move (except for a
    /// successful re-target).
    Stuck,
    /// A gradient stop rule holds; no line search was done.
    Converged,
    /// The re-target projection failed even with the smallest increment, and the descent did
    /// not move the curve either.
    RetargetStalled,
}

/// `√(Σ_v m_v |g_v|²)`.
pub fn mass_norm(direction: &[[f64; 3]], masses: &[f64]) -> f64 {
    direction
        .iter()
        .zip(masses)
        .map(|(g, m)| m * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]))
        .sum::<f64>()
        .sqrt()
}

/// What the stop rules remember between steps.
#[derive(Debug, Clone, Default)]
pub struct StopState {
    first_norm: Option<f64>,
    energies: VecDeque<f64>,
}

impl StopState {
    /// Forgets the reference norm and the energy history (after a subdivision the norms are not
    /// comparable).
    pub fn reset(&mut self) {
        self.first_norm = None;
        self.energies.clear();
    }

    /// Whether a gradient rule holds for the direction with Euclidean norm `euclidean` and
    /// mass norm `mass` (computed by the caller). `EnergyStall` never holds here.
    pub fn gradient_converged(&mut self, rule: &StopRule, euclidean: f64, mass: f64) -> bool {
        let first = *self.first_norm.get_or_insert(euclidean);
        match *rule {
            StopRule::GradientNorm { norm, eps } => match norm {
                GradNorm::Euclidean => euclidean < eps,
                GradNorm::MassL2 => mass < eps,
            },
            StopRule::RelativeGradient { ratio } => {
                euclidean == 0.0 || (first > 0.0 && euclidean / first < ratio)
            }
            StopRule::EnergyStall { .. } => euclidean == 0.0,
        }
    }

    /// Records the objective at the base of the first step after a (re)start, so that the
    /// first window of `EnergyStall` spans `window` steps from there.
    pub fn seed_energy(&mut self, energy: f64) {
        if self.energies.is_empty() {
            self.energies.push_back(energy);
        }
    }

    /// Records the objective after an accepted step; whether `EnergyStall` holds.
    pub fn energy_stalled(&mut self, rule: &StopRule, energy: f64) -> bool {
        let StopRule::EnergyStall { relative, window } = *rule else {
            return false;
        };
        self.energies.push_back(energy);
        while self.energies.len() > window + 1 {
            self.energies.pop_front();
        }
        if self.energies.len() < window + 1 {
            return false;
        }
        let old = self.energies[0];
        old - energy <= relative * energy.abs()
    }
}
