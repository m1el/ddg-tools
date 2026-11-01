//! Penalties of `SPEC.md` §5.2 (area/volume deviation, total area/volume, squared error,
//! boundary length and curvature, Willmore, implicit obstacles and attractors, mesh and
//! point-cloud obstacles).
//!
//! Stub created by WP1; filled by WP7 (SPEC §15). [`Penalty`] is empty until then; the hook
//! `evaluate` is called by [`crate::energy::energy`], and WP7 replaces both (and may change
//! everything else in this file).

use crate::error::RsError;
use crate::surface::Surface;

/// A soft penalty `w · P(f)` (`SPEC.md` §5.2). No variants until WP7.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Penalty {}

/// Hook for [`crate::energy::energy`]: returns the weighted value `w_i P_i(f)` of every
/// penalty (in order) and adds the weighted gradients to `differential`. Positions are
/// validated.
///
/// # Errors
/// None until WP7.
pub(crate) fn evaluate(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    penalties: &[Penalty],
    differential: &mut [[f64; 3]],
) -> Result<Vec<f64>, RsError> {
    let _ = (pool, surface, positions, differential);
    penalties.iter().map(|p| match *p {}).collect()
}
