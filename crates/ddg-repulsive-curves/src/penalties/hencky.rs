//! Hencky strain penalty (`SPEC.md` §5.3, D21).
//!
//! The penalty `w Σ_I ℓ⁰_I ln²(ℓ_I/ℓ⁰_I)` was implemented by work package WP3 with the soft
//! methods and stays where it is, as [`crate::optimize::Hencky`] (the session, the soft
//! problem and its tests use it there); it is re-exported here so that all penalties are
//! reachable from [`crate::penalties`]. This module adds the weight rule of D21 (§14.2).

use std::ops::RangeInclusive;

use crate::error::RcError;
pub use crate::optimize::Hencky;

/// The calibration rule of `SPEC.md` §14.2 (D21): the smallest weight `w_H = 10^k`, `k` in
/// `exponents`, for which `median_strain(w_H) <= target`, where `median_strain` runs the
/// reference method on the calibration knots with the Hencky weight `w_H` and returns the median
/// over the knots of the final largest `|ln(ℓ/ℓ⁰)|` ([`Hencky::max_log_strain`]). `None` if no
/// weight in the range reaches the target.
///
/// # Errors
/// [`RcError::InvalidOption`] for a non-positive or non-finite `target`; the errors of
/// `median_strain`.
pub fn calibrate_hencky_weight(
    exponents: RangeInclusive<i32>,
    target: f64,
    mut median_strain: impl FnMut(f64) -> Result<f64, RcError>,
) -> Result<Option<f64>, RcError> {
    if !(target.is_finite() && target > 0.0) {
        return Err(RcError::InvalidOption {
            name: "calibrate_hencky_weight.target",
            value: target.to_string(),
            requirement: "a finite positive strain",
        });
    }
    for k in exponents {
        let w = 10f64.powi(k);
        if median_strain(w)? <= target {
            return Ok(Some(w));
        }
    }
    Ok(None)
}
