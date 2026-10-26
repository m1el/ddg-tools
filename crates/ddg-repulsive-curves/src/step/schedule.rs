//! Target-length schedules (`SPEC.md` §5.4, D13).
//!
//! A [`LengthSchedule`] `{ factor, rate }` moves the target of every `TotalLength` constraint and
//! all targets of every `EdgeLengths` constraint by one common factor
//!
//! ```text
//! c(p) = 1 + sign(factor − 1) · min(p, |factor − 1|),     L*_k = c(p_k) L⁰,   ℓ*_I = c(p_k) ℓ⁰_I
//! ```
//!
//! where `p_k` is the progress after step `k` (`k · rate` when every re-target succeeds; a failed
//! re-target retries with half the increment, up to 5 times, D13). With one common factor the
//! ratios of the edge targets are kept (the reference adds equal absolute amounts instead).
//! The bases `L⁰`, `ℓ⁰_I` are the targets when the schedule was set (subdivision halves the edge
//! bases with the targets).
//!
//! [formula: rc-length-schedule]

use crate::constraints::{ConstraintSet, ResolvedConstraint};
use crate::error::RcError;
use crate::subdivide::SubdivisionMap;

/// A length schedule (`SPEC.md` §5.4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LengthSchedule {
    /// The final length relative to the start (`> 1` grows, `< 1` shrinks, `1` does nothing).
    pub factor: f64,
    /// Progress per step, as a fraction of the start length (the reference: 0.01).
    pub rate: f64,
}

impl Default for LengthSchedule {
    /// `factor 1` (no change), `rate 0.01` (`SPEC.md` §10).
    fn default() -> Self {
        LengthSchedule {
            factor: 1.0,
            rate: 0.01,
        }
    }
}

impl LengthSchedule {
    /// Checks the parameters.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-positive or non-finite factor or rate.
    pub fn validate(&self) -> Result<(), RcError> {
        for (name, v) in [
            ("schedule.factor", self.factor),
            ("schedule.rate", self.rate),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(RcError::InvalidOption {
                    name,
                    value: v.to_string(),
                    requirement: "a finite positive number",
                });
            }
        }
        Ok(())
    }

    /// The common factor `c(p) = 1 + sign(factor − 1) · min(p, |factor − 1|)`.
    ///
    /// [formula: rc-length-schedule]
    pub fn scale_at(&self, progress: f64) -> f64 {
        let span = (self.factor - 1.0).abs();
        let p = progress.clamp(0.0, span);
        if self.factor >= 1.0 { 1.0 + p } else { 1.0 - p }
    }

    /// The total target after `k` successful steps: `L⁰ + sign(f − 1) · min(k · rate · L⁰,
    /// |f − 1| L⁰)`.
    ///
    /// [formula: rc-length-schedule]
    pub fn total_target(&self, l0: f64, k: usize) -> f64 {
        l0 * self.scale_at(k as f64 * self.rate)
    }

    /// `|factor − 1|`, the progress at which the schedule ends.
    pub fn span(&self) -> f64 {
        (self.factor - 1.0).abs()
    }
}

/// The base targets of one constraint.
#[derive(Debug, Clone, PartialEq)]
enum Base {
    Total { constraint: usize, l0: f64 },
    Edges { constraint: usize, l0: Vec<f64> },
}

/// A schedule with its progress and base targets.
#[derive(Debug, Clone, PartialEq)]
pub struct ScheduleState {
    schedule: LengthSchedule,
    progress: f64,
    bases: Vec<Base>,
}

impl ScheduleState {
    /// Starts `schedule` from the current targets of `constraints`.
    ///
    /// # Errors
    /// [`RcError::ScheduleWithoutLengthConstraint`] if there is neither a `TotalLength` nor an
    /// `EdgeLengths` constraint; the errors of [`LengthSchedule::validate`].
    pub fn new(
        schedule: LengthSchedule,
        constraints: &ConstraintSet,
    ) -> Result<ScheduleState, RcError> {
        schedule.validate()?;
        let mut bases = Vec::new();
        for (ci, c) in constraints.constraints().iter().enumerate() {
            match c {
                ResolvedConstraint::TotalLength { target } => bases.push(Base::Total {
                    constraint: ci,
                    l0: *target,
                }),
                ResolvedConstraint::EdgeLengths { targets, .. } => bases.push(Base::Edges {
                    constraint: ci,
                    l0: targets.clone(),
                }),
                _ => {}
            }
        }
        if bases.is_empty() {
            return Err(RcError::ScheduleWithoutLengthConstraint);
        }
        Ok(ScheduleState {
            schedule,
            progress: 0.0,
            bases,
        })
    }

    /// The schedule.
    pub fn schedule(&self) -> LengthSchedule {
        self.schedule
    }

    /// The progress reached so far.
    pub fn progress(&self) -> f64 {
        self.progress
    }

    /// The current common factor.
    pub fn scale(&self) -> f64 {
        self.schedule.scale_at(self.progress)
    }

    /// `true` when the targets have reached `factor` times their bases.
    pub fn finished(&self) -> bool {
        self.progress >= self.schedule.span()
    }

    /// The progress after one step whose increment is halved `halvings` times.
    pub fn next_progress(&self, halvings: u32) -> f64 {
        (self.progress + self.schedule.rate * 0.5f64.powi(halvings as i32))
            .min(self.schedule.span())
    }

    /// Writes the targets of progress `progress` into `constraints`.
    ///
    /// # Errors
    /// The errors of the target setters of [`ConstraintSet`].
    pub fn apply(&self, constraints: &mut ConstraintSet, progress: f64) -> Result<(), RcError> {
        let c = self.schedule.scale_at(progress);
        for b in &self.bases {
            match b {
                Base::Total { constraint, l0 } => {
                    constraints.set_total_length_target(*constraint, c * l0)?
                }
                Base::Edges { constraint, l0 } => {
                    let t: Vec<f64> = l0.iter().map(|x| c * x).collect();
                    constraints.set_edge_length_targets(*constraint, &t)?
                }
            }
        }
        Ok(())
    }

    /// Records a successful re-target to `progress`.
    pub fn commit(&mut self, progress: f64) {
        self.progress = progress;
    }

    /// The state after a subdivision: constraint indices through `index` (old → new indices, as
    /// returned by [`crate::subdivide::subdivide_constraints`]), edge bases halved on both
    /// children.
    pub fn subdivided(&self, map: &SubdivisionMap, index: &[Vec<usize>]) -> ScheduleState {
        let bases = self
            .bases
            .iter()
            .map(|b| match b {
                Base::Total { constraint, l0 } => Base::Total {
                    constraint: index[*constraint][0],
                    l0: *l0,
                },
                Base::Edges { constraint, l0 } => Base::Edges {
                    constraint: index[*constraint][0],
                    l0: map.map_edge_values(l0, |x| 0.5 * x),
                },
            })
            .collect();
        ScheduleState {
            schedule: self.schedule,
            progress: self.progress,
            bases,
        }
    }
}
