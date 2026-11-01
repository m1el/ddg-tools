//! Barnes–Hut energy and the exact gradient of the approximated energy for frozen lists
//! (`SPEC.md` §3.4, D4, D5).
//!
//! Stub created by WP1; filled by WP3 (SPEC §15). The hook below is called by
//! [`crate::energy::energy`]; WP3 replaces its body (and may change everything else in this
//! file).

use crate::energy::{BhOptions, EnergyEvaluation};
use crate::error::RsError;
use crate::geometry::FaceGeometry;
use crate::kernel::Exponents;
use crate::surface::Surface;

/// Hook for [`crate::energy::energy`]: `Ẽ`, `∂Ẽ/∂f_v` and the pair/node counts.
/// `geom` holds the validated face quantities of `positions`; `options` are validated.
///
/// # Errors
/// [`RsError::Unsupported`] until WP3 fills it in.
pub(crate) fn energy_and_differential(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    geom: &FaceGeometry,
    exponents: &Exponents,
    options: &BhOptions,
) -> Result<EnergyEvaluation, RsError> {
    let _ = (pool, surface, positions, geom, exponents, options);
    Err(RsError::Unsupported {
        what: "the Barnes–Hut energy (SPEC §15 WP3)",
    })
}
