//! Intrinsic mollification: making a degenerate metric usable.
//!
//! The function and its record live in `ddg_core::mollify` (moved there so that crates that do
//! not depend on this one can use it); this module re-exports them and keeps
//! [`mollify_edge_lengths`] with this crate's error type, so callers are unchanged.
//!
//! **Nothing in this crate mollifies on its own.** The flattenings take a
//! `ddg_core::Geometry` and use it as it is; a mesh with a degenerate face is refused when
//! that geometry is built (`GeometryError::DegenerateFaces` names the faces). A caller who
//! wants the reference's behaviour calls [`mollify_edge_lengths`] and gets back, besides the
//! lengths, how much was added and how many faces made it necessary ([`Mollified`]).

use ddg_core::mollify::MollifyError;
use ddg_core::topology::TriMesh;

pub use ddg_core::mollify::{Mollified, REFERENCE_MOLLIFY_DELTA};

use crate::error::BffError;

/// [`ddg_core::mollify::mollify_edge_lengths`] with the errors of this crate: lengthens all
/// edges by the smallest common amount `ε ≥ 0` such that every triangle satisfies
/// `ℓ_a + ℓ_b ≥ ℓ_c + δ` for each of its sides `c`, where `δ = relative_delta · (mean edge
/// length)`.
///
/// # Errors
/// [`BffError::LengthMismatch`] or [`BffError::NonFiniteInput`] for bad lengths or a negative
/// or non-finite `relative_delta`.
pub fn mollify_edge_lengths(
    mesh: &TriMesh,
    edge_lengths: &[f64],
    relative_delta: f64,
) -> Result<Mollified, BffError> {
    ddg_core::mollify::mollify_edge_lengths(mesh, edge_lengths, relative_delta).map_err(|e| match e
    {
        MollifyError::LengthCountMismatch { expected, got } => BffError::LengthMismatch {
            what: "edge lengths",
            expected,
            got,
        },
        MollifyError::NonFiniteLength { edge, value } => BffError::NonFiniteInput {
            what: "edge lengths",
            index: edge,
            value,
        },
        MollifyError::InvalidDelta { value } => BffError::NonFiniteInput {
            what: "relative_delta (must be finite and non-negative)",
            index: 0,
            value,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn errors_are_those_of_this_crate() {
        let mesh = TriMesh::new(3, vec![[0, 1, 2]]).unwrap();
        assert!(matches!(
            mollify_edge_lengths(&mesh, &[1.0, 1.0], 1e-6),
            Err(BffError::LengthMismatch {
                what: "edge lengths",
                expected: 3,
                got: 2
            })
        ));
        assert!(matches!(
            mollify_edge_lengths(&mesh, &[1.0, f64::INFINITY, 1.0], 1e-6),
            Err(BffError::NonFiniteInput { index: 1, .. })
        ));
        assert!(matches!(
            mollify_edge_lengths(&mesh, &[1.0, 1.0, 1.0], f64::NAN),
            Err(BffError::NonFiniteInput { index: 0, .. })
        ));
        let m = mollify_edge_lengths(&mesh, &[1.0, 1.0, 2.0], 1e-3).unwrap();
        assert!(m.epsilon > 0.0);
    }
}
