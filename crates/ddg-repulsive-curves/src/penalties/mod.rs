//! Penalties of `SPEC.md` §5.3 (Sec. 8.1 of the paper), implementing
//! [`PenaltyTerm`](crate::step::PenaltyTerm):
//!
//! | penalty | module |
//! |---|---|
//! | [`TotalLengthPenalty`], [`LengthDifference`] | [`length`] |
//! | [`FieldAlignment`] | [`field`] |
//! | [`SurfacePotential`] | [`surface`] |
//! | Hencky (kept in [`crate::optimize::Hencky`]) | [`hencky`] |
//!
//! **D20.** The authors' versions (`src/extra_potentials.cpp`, MIT, Copyright (c) 2019
//! Christopher Yu) are not copied where they are inconsistent: their total length is `½ w L²`
//! with a gradient that indexes `GetVertex(j)` (the neighbour's *local* index) as a global
//! vertex; their length difference carries a factor ½ that the paper's `w Σ (ℓ_I − ℓ_J)²` does
//! not; their field energy `w Σ |X × T|²` omits the factor `ℓ_I` that its gradient contains.
//! Here value and gradient are the derivative of the formula printed in the paper and in
//! `SPEC.md` §5.3, and every gradient is tested against central differences.

pub mod field;
pub mod hencky;
pub mod length;
pub mod surface;

pub use field::FieldAlignment;
pub use hencky::{Hencky, calibrate_hencky_weight};
pub use length::{LengthDifference, TotalLengthPenalty};
pub use surface::SurfacePotential;

use crate::error::RcError;
use crate::network::CurveNetwork;

/// The checks every penalty's `evaluate` makes on its arguments.
pub(crate) fn check_args(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    gradient: &Option<&mut [[f64; 3]]>,
) -> Result<(), RcError> {
    if positions.len() != network.n_vertices() {
        return Err(RcError::PositionCountMismatch {
            expected: network.n_vertices(),
            found: positions.len(),
        });
    }
    if let Some(g) = gradient
        && g.len() != network.n_vertices()
    {
        return Err(RcError::PositionCountMismatch {
            expected: network.n_vertices(),
            found: g.len(),
        });
    }
    Ok(())
}

/// A finite weight (any sign when `allow_negative`).
pub(crate) fn check_weight(
    name: &'static str,
    w: f64,
    allow_negative: bool,
) -> Result<(), RcError> {
    if w.is_finite() && (allow_negative || w >= 0.0) {
        Ok(())
    } else {
        Err(RcError::InvalidOption {
            name,
            value: w.to_string(),
            requirement: if allow_negative {
                "a finite weight"
            } else {
                "a finite non-negative weight"
            },
        })
    }
}
