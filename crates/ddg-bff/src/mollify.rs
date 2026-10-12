//! Intrinsic mollification: making a degenerate metric usable.
//!
//! `ddg_core::Geometry` rejects faces whose edge lengths do not satisfy the strict triangle
//! inequality. The reference implementation of Boundary First Flattening instead lengthens
//! every edge by a common small amount until all triangle inequalities hold with a margin
//! (`BFFData::computeEdgeLengths`; the technique is from N. Sharp and K. Crane, *A Laplacian
//! for Nonmanifold Triangle Meshes*, SGP 2020, Sec. 4.5). This module offers the same as an
//! explicit step: compute lengths, mollify them, then build the geometry from the result.

use ddg_core::topology::TriMesh;

use crate::error::{BffError, check_values};

/// The relative margin the reference implementation uses: `1e-6` of the mean edge length.
pub const REFERENCE_MOLLIFY_DELTA: f64 = 1e-6;

/// Mollified edge lengths.
#[derive(Debug, Clone, PartialEq)]
pub struct Mollified {
    /// The new length of every edge: the old one plus `epsilon`.
    pub edge_lengths: Vec<f64>,
    /// The amount that was added to every edge. Zero if all triangle inequalities already
    /// held with the requested margin; the lengths are then unchanged.
    pub epsilon: f64,
    /// The absolute margin `δ = relative_delta · mean edge length`.
    pub delta: f64,
}

/// Lengthens all edges by the smallest common amount `ε ≥ 0` such that every triangle
/// satisfies `ℓ_a + ℓ_b ≥ ℓ_c + δ` for each of its sides `c`, where
/// `δ = relative_delta · (mean edge length)`:
/// `ε = max(0, max over corners (ℓ_c − ℓ_a − ℓ_b) + δ)`.
///
/// # Errors
/// [`BffError::LengthMismatch`] or [`BffError::NonFiniteInput`] for bad lengths or a negative
/// or non-finite `relative_delta`.
///
/// [formula: intrinsic-mollification]
pub fn mollify_edge_lengths(
    mesh: &TriMesh,
    edge_lengths: &[f64],
    relative_delta: f64,
) -> Result<Mollified, BffError> {
    check_values("edge lengths", edge_lengths, mesh.n_edges())?;
    if !relative_delta.is_finite() || relative_delta < 0.0 {
        return Err(BffError::NonFiniteInput {
            what: "relative_delta (must be finite and non-negative)",
            index: 0,
            value: relative_delta,
        });
    }
    let mean = edge_lengths.iter().sum::<f64>() / edge_lengths.len() as f64;
    let delta = relative_delta * mean;
    let mut epsilon = 0.0f64;
    for f in 0..mesh.n_faces() {
        let l = [0, 1, 2].map(|k| edge_lengths[mesh.edge(3 * f + k)]);
        for k in 0..3 {
            epsilon = epsilon.max(l[k] - l[(k + 1) % 3] - l[(k + 2) % 3] + delta);
        }
    }
    Ok(Mollified {
        edge_lengths: edge_lengths.iter().map(|l| l + epsilon).collect(),
        epsilon,
        delta,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::geometry::Geometry;

    #[test]
    fn valid_triangles_are_left_alone_and_degenerate_ones_repaired() {
        let mesh = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
        // Edges in id order: (0,1), (0,2), (0,3), (1,2), (2,3).
        let good = [1.0, 1.2, 1.0, 1.0, 1.1];
        let m = mollify_edge_lengths(&mesh, &good, REFERENCE_MOLLIFY_DELTA).unwrap();
        assert_eq!(m.epsilon, 0.0);
        assert_eq!(m.edge_lengths, good);
        // Face (0, 1, 2) with lengths 1, 1, 2 is degenerate.
        let bad = [1.0, 2.0, 1.5, 1.0, 1.5];
        assert!(Geometry::from_edge_lengths(&mesh, bad.to_vec()).is_err());
        let m = mollify_edge_lengths(&mesh, &bad, 1e-3).unwrap();
        assert!((m.epsilon - m.delta).abs() < 1e-15 && m.epsilon > 0.0);
        let geometry = Geometry::from_edge_lengths(&mesh, m.edge_lengths.clone()).unwrap();
        assert!(geometry.face_areas().iter().all(|&a| a > 0.0));
        // Every triangle inequality now holds with margin delta.
        for f in 0..2 {
            let l = [0, 1, 2].map(|k| m.edge_lengths[mesh.edge(3 * f + k)]);
            for k in 0..3 {
                assert!(l[(k + 1) % 3] + l[(k + 2) % 3] - l[k] >= m.delta * (1.0 - 1e-9));
            }
        }
    }
}
