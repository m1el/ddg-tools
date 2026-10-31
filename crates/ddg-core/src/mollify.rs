//! Intrinsic mollification: making a degenerate metric usable.
//!
//! [`Geometry`](crate::Geometry) rejects faces whose edge lengths do not satisfy the strict
//! triangle inequality. The reference implementation of Boundary First Flattening (Sawhney and
//! Crane 2017) instead lengthens every edge by a common small amount until all triangle
//! inequalities hold with a margin (`BFFData::computeEdgeLengths`; the technique is from
//! N. Sharp and K. Crane, *A Laplacian for Nonmanifold Triangle Meshes*, SGP 2020, Sec. 4.5).
//! This module offers the same as an explicit step: compute lengths, mollify them, then build
//! the geometry from the result. Integer coordinates (Gillespie, Sharp and Crane 2021, Sec. 3.9)
//! use it as their robustness preprocessing.
//!
//! **Nothing mollifies on its own.** A caller who wants it calls [`mollify_edge_lengths`] and
//! gets back, besides the lengths, how much was added and how many faces made it necessary
//! ([`Mollified`]) — that is the record to put into the caller's report. With
//! [`Mollified::epsilon`]` == 0` the lengths are returned bit for bit.
//!
//! The reference mollifies always. Its margin (`1e-6` of the mean edge length) is also
//! exceeded by valid but very thin triangles, so it changes the input of some meshes that
//! `Geometry` accepts: of the 23 meshes of the optimal-cones dataset, `camelbad` and
//! `camelbadnonuniform`, by `2e-8` of the mesh radius (see the README of `ddg-bff`).
//!
//! It lived in `ddg-bff` first; `ddg_bff::mollify_edge_lengths` is a wrapper around this one
//! with that crate's error type.

use std::fmt;

use crate::topology::TriMesh;

/// The relative margin the reference implementation uses: `1e-6` of the mean edge length.
pub const REFERENCE_MOLLIFY_DELTA: f64 = 1e-6;

/// Mollified edge lengths, and the record of what was done.
#[derive(Debug, Clone, PartialEq)]
pub struct Mollified {
    /// The new length of every edge: the old one plus `epsilon`.
    pub edge_lengths: Vec<f64>,
    /// The amount that was added to every edge. Zero if all triangle inequalities already
    /// held with the requested margin; the lengths are then unchanged.
    pub epsilon: f64,
    /// The absolute margin `δ = relative_delta · mean edge length`.
    pub delta: f64,
    /// Number of faces that missed the margin before mollification, i.e. the faces because of
    /// which `epsilon > 0`. Zero exactly when `epsilon == 0`.
    pub n_faces_below_margin: usize,
    /// The face that missed the margin by most (the one that determines `epsilon`), if any.
    pub worst_face: Option<u32>,
}

/// Why [`mollify_edge_lengths`] refused its input.
#[derive(Debug, Clone, PartialEq)]
pub enum MollifyError {
    /// There is not one length per edge.
    LengthCountMismatch {
        /// The number of edges.
        expected: usize,
        /// The number of lengths passed.
        got: usize,
    },
    /// A length is not finite.
    NonFiniteLength {
        /// The edge.
        edge: usize,
        /// Its length.
        value: f64,
    },
    /// The relative margin is negative or not finite.
    InvalidDelta {
        /// The margin passed.
        value: f64,
    },
}

impl fmt::Display for MollifyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MollifyError::LengthCountMismatch { expected, got } => {
                write!(f, "edge lengths has {got} entries, expected {expected}")
            }
            MollifyError::NonFiniteLength { edge, value } => {
                write!(f, "edge lengths[{edge}] = {value} is not finite")
            }
            MollifyError::InvalidDelta { value } => {
                write!(f, "relative_delta = {value} is not finite and non-negative")
            }
        }
    }
}

impl std::error::Error for MollifyError {}

/// Lengthens all edges by the smallest common amount `ε ≥ 0` such that every triangle
/// satisfies `ℓ_a + ℓ_b ≥ ℓ_c + δ` for each of its sides `c`, where
/// `δ = relative_delta · (mean edge length)`:
/// `ε = max(0, max over corners (ℓ_c − ℓ_a − ℓ_b) + δ)`.
///
/// # Errors
/// [`MollifyError`] for a number of lengths other than the number of edges, a length that is
/// not finite, or a negative or non-finite `relative_delta`.
///
/// [formula: intrinsic-mollification]
pub fn mollify_edge_lengths(
    mesh: &TriMesh,
    edge_lengths: &[f64],
    relative_delta: f64,
) -> Result<Mollified, MollifyError> {
    if edge_lengths.len() != mesh.n_edges() {
        return Err(MollifyError::LengthCountMismatch {
            expected: mesh.n_edges(),
            got: edge_lengths.len(),
        });
    }
    if let Some(edge) = edge_lengths.iter().position(|x| !x.is_finite()) {
        return Err(MollifyError::NonFiniteLength {
            edge,
            value: edge_lengths[edge],
        });
    }
    if !relative_delta.is_finite() || relative_delta < 0.0 {
        return Err(MollifyError::InvalidDelta {
            value: relative_delta,
        });
    }
    let mean = edge_lengths.iter().sum::<f64>() / edge_lengths.len() as f64;
    let delta = relative_delta * mean;
    let mut epsilon = 0.0f64;
    let mut n_faces_below_margin = 0usize;
    let mut worst_face = None;
    for f in 0..mesh.n_faces() {
        let l = [0, 1, 2].map(|k| edge_lengths[mesh.edge(3 * f + k)]);
        let mut needed = 0.0f64;
        for k in 0..3 {
            needed = needed.max(l[k] - l[(k + 1) % 3] - l[(k + 2) % 3] + delta);
        }
        if needed > 0.0 {
            n_faces_below_margin += 1;
            if needed > epsilon {
                epsilon = needed;
                worst_face = Some(f as u32);
            }
        }
    }
    Ok(Mollified {
        edge_lengths: edge_lengths.iter().map(|l| l + epsilon).collect(),
        epsilon,
        delta,
        n_faces_below_margin,
        worst_face,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Geometry;

    #[test]
    fn valid_triangles_are_left_alone_and_degenerate_ones_repaired() {
        let mesh = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
        // Edges in id order: (0,1), (0,2), (0,3), (1,2), (2,3).
        let good = [1.0, 1.2, 1.0, 1.0, 1.1];
        let m = mollify_edge_lengths(&mesh, &good, REFERENCE_MOLLIFY_DELTA).unwrap();
        assert_eq!(m.epsilon, 0.0);
        assert_eq!(m.edge_lengths, good);
        assert_eq!((m.n_faces_below_margin, m.worst_face), (0, None));
        // Face (0, 1, 2) with lengths 1, 1, 2 is degenerate.
        let bad = [1.0, 2.0, 1.5, 1.0, 1.5];
        assert!(Geometry::from_edge_lengths(&mesh, bad.to_vec()).is_err());
        let m = mollify_edge_lengths(&mesh, &bad, 1e-3).unwrap();
        assert!((m.epsilon - m.delta).abs() < 1e-15 && m.epsilon > 0.0);
        assert_eq!((m.n_faces_below_margin, m.worst_face), (1, Some(0)));
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

    #[test]
    fn bad_input_is_refused() {
        let mesh = TriMesh::new(3, vec![[0, 1, 2]]).unwrap();
        assert_eq!(
            mollify_edge_lengths(&mesh, &[1.0, 1.0], 1e-6),
            Err(MollifyError::LengthCountMismatch {
                expected: 3,
                got: 2
            })
        );
        assert!(matches!(
            mollify_edge_lengths(&mesh, &[1.0, f64::NAN, 1.0], 1e-6),
            Err(MollifyError::NonFiniteLength { edge: 1, .. })
        ));
        assert_eq!(
            mollify_edge_lengths(&mesh, &[1.0, 1.0, 1.0], -1.0),
            Err(MollifyError::InvalidDelta { value: -1.0 })
        );
    }
}
