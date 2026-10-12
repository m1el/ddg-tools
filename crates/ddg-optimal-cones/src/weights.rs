//! The curvature-adapted weights of the paper's Sec. 8.1 ("nonuniform importance").
//!
//! This is a convenience helper. Nothing else in the crate uses it, and the solver never
//! applies such weights on its own.

use ddg_core::operators::{angle_defects, lumped_mass_normalized};
use ddg_core::{Geometry, TriMesh};

use crate::error::ConeError;

/// Which of the paper's two contradictory statements about the local feature size `r` to use.
///
/// Sec. 8.1 prints `r_i := 1/(|Ω_i|/M_ii + ε)` and then says "this value will be small in flat
/// regions and large in highly curved regions. We then set `(w_E)_i := r_i` and
/// `(w_R)_i := 1/r_i`, emphasizing the importance of small features, and decreasing the cost of
/// placing cones in those same regions." The printed formula is *large* in flat regions, so
/// with it the weights do the opposite of what the text describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureSizeConvention {
    /// `r_i = |Ω_i|/m_i + ε`: small in flat regions and large in curved ones, as the paper's
    /// text describes. `wE = r` then emphasizes distortion at small features and `wR = 1/r`
    /// makes cones cheaper there, which is the stated purpose.
    AsDescribed,
    /// `r_i = 1/(|Ω_i|/m_i + ε)`: the formula exactly as printed. `wE = r` then emphasizes
    /// flat regions and `wR = 1/r` makes cones cheaper there.
    AsPrinted,
}

/// Per-vertex weights for [`ConeProblem`](crate::ConeProblem).
#[derive(Debug, Clone, PartialEq)]
pub struct FeatureSizeWeights {
    /// The feature-size value `r_i` per vertex.
    pub feature_size: Vec<f64>,
    /// Distortion weights `wE_i = r_i`.
    pub distortion_weights: Vec<f64>,
    /// Cone costs `wR_i = 1/r_i`.
    pub cone_costs: Vec<f64>,
}

/// The weights `wE_i = r_i`, `wR_i = 1/r_i` of the paper's Sec. 8.1, where `r` is built from the
/// curvature density `|Ω_i|/m_i` according to `convention` (see [`FeatureSizeConvention`] for
/// why there are two).
///
/// `Ω_i` is the angle defect and `m_i` the lumped mass **divided by the total area**, so
/// `|Ω_i|/m_i` and `epsilon` are dimensionless and the weights do not depend on the length
/// unit. `Ω_i` is zero at boundary vertices, so `r_i` is `ε` resp. `1/ε` there (the weights of
/// boundary vertices have no effect on the problem). `epsilon` must be positive and finite; it
/// keeps `r` and `1/r` finite at flat vertices.
///
/// The weights are returned as they are, without any rescaling. Their overall magnitude
/// multiplies into the problem: scaling all `wE` by `a` and all `wR` by `b` is the same as
/// replacing `λ` by `λ·b/a`.
///
/// # Errors
/// [`ConeError::Geometry`] if the geometry does not belong to the mesh;
/// [`ConeError::InvalidArgument`] if `epsilon` is not positive and finite.
///
/// [formula: feature-size-weights]
pub fn feature_size_weights(
    mesh: &TriMesh,
    geometry: &Geometry,
    epsilon: f64,
    convention: FeatureSizeConvention,
) -> Result<FeatureSizeWeights, ConeError> {
    if !(epsilon.is_finite() && epsilon > 0.0) {
        return Err(ConeError::InvalidArgument {
            reason: format!("epsilon must be positive and finite, got {epsilon}"),
        });
    }
    let omega = angle_defects(mesh, geometry)?;
    let mass = lumped_mass_normalized(mesh, geometry)?;
    let feature_size: Vec<f64> = omega
        .iter()
        .zip(&mass)
        .map(|(omega, mass)| {
            let density = omega.abs() / mass + epsilon;
            match convention {
                FeatureSizeConvention::AsDescribed => density,
                FeatureSizeConvention::AsPrinted => 1.0 / density,
            }
        })
        .collect();
    Ok(FeatureSizeWeights {
        distortion_weights: feature_size.clone(),
        cone_costs: feature_size.iter().map(|r| 1.0 / r).collect(),
        feature_size,
    })
}
