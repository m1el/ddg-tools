//! Discrete operators and quality measures on an intrinsic triangulation.
//!
//! The triangulation is not a `ddg_core::TriMesh`, so `ddg_core::operators` cannot be called on
//! it; these are the same operators with the same conventions (the repository's `CLAUDE.md`),
//! assembled from the same per-triangle formulas of `ddg_core::geometry`. Double edges add
//! their weights into one matrix entry; a self-edge contributes nothing to the Laplacian.

use std::f64::consts::{PI, TAU};

use ddg_core::SparseMatrix;
use ddg_core::geometry::corner_cotan;

use crate::triangulation::IntrinsicTriangulation;

/// Angle and cotan-weight statistics of a triangulation. See
/// [`IntrinsicTriangulation::quality`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quality {
    /// Number of faces.
    pub n_faces: usize,
    /// Number of interior edges.
    pub n_interior_edges: usize,
    /// The smallest corner angle, in radians.
    pub min_angle: f64,
    /// The largest corner angle, in radians.
    pub max_angle: f64,
    /// The median, the 90 % and the 99 % quantile and the maximum of the largest angle of a
    /// face, in radians.
    pub largest_angle_quantiles: [f64; 4],
    /// The fractions of faces whose largest angle exceeds 90°, 120° and 150° by more than
    /// `1e-9` radians (so that right triangles are not counted as obtuse by rounding).
    pub largest_angle_above: [f64; 3],
    /// Number of interior edges whose cotan weight is below minus the threshold passed to
    /// [`IntrinsicTriangulation::quality`].
    pub n_negative_weights: usize,
    /// The smallest cotan weight of an interior edge (`+∞` if there is none).
    pub min_cotan_weight: f64,
}

impl Quality {
    /// The fraction of interior edges with a negative cotan weight (zero if there is no
    /// interior edge).
    pub fn negative_weight_fraction(&self) -> f64 {
        if self.n_interior_edges == 0 {
            0.0
        } else {
            self.n_negative_weights as f64 / self.n_interior_edges as f64
        }
    }
}

impl IntrinsicTriangulation {
    /// The cotangent of the angle opposite halfedge `h` in its face.
    pub(crate) fn opposite_cotan(&self, h: usize) -> f64 {
        let l = self.face_lengths(h / 3);
        let k = h % 3;
        corner_cotan(l[k], l[(k + 1) % 3], l[(k + 2) % 3])
            .expect("every face keeps the triangle inequalities")
    }

    /// The cotan Laplacian of the triangulation as an `n_vertices × n_vertices` matrix,
    /// symmetric **positive semidefinite** with zero row sums:
    /// `L_ij = −½ Σ cot θ` over the corners opposite the edges between `i` and `j`, and
    /// `L_ii = −Σ_{j≠i} L_ij`.
    ///
    /// It is the stiffness matrix of piecewise linear finite elements on the triangulation.
    /// Several edges between `i` and `j` add up in the entry `(i, j)`; an edge from a vertex
    /// to itself contributes nothing (its two ends have the same value). As in
    /// `ddg_core::operators::cotan_laplacian` the matrix is exactly symmetric and its pattern
    /// is the vertex adjacency plus the diagonal, with explicit zeros. After Delaunay flipping
    /// the off-diagonal entries are not positive (up to the tolerance of the criterion), which
    /// is the maximum principle the paper's Sec. 4 refers to.
    ///
    /// [formula: cotan-laplacian]
    pub fn cotan_laplacian(&self) -> SparseMatrix {
        let n = self.n_vertices();
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(6 * self.n_faces() + n);
        for v in 0..n as u32 {
            triplets.push((v, v, 0.0));
        }
        for h in 0..self.n_halfedges() {
            // The corner opposite halfedge h weighs the edge of h.
            let (i, j) = (self.from_vertex(h) as u32, self.to_vertex(h) as u32);
            if i == j {
                continue;
            }
            let w = 0.5 * self.opposite_cotan(h);
            triplets.push((i, j, -w));
            triplets.push((j, i, -w));
        }
        let mut l = SparseMatrix::from_triplets(n, n, &triplets)
            .expect("vertex indices of a valid triangulation are in range");
        let col_ptr = l.col_ptr().to_vec();
        for j in 0..n {
            let range = col_ptr[j]..col_ptr[j + 1];
            let mut sum = 0.0;
            let mut diagonal_slot = range.start;
            for k in range {
                if l.row_idx()[k] as usize == j {
                    diagonal_slot = k;
                } else {
                    sum += l.values()[k];
                }
            }
            l.values_mut()[diagonal_slot] = -sum;
        }
        l
    }

    /// The lumped vertex areas `m_i = ⅓ Σ A_f` over the corners of faces `f` at vertex `i` (a
    /// face that uses a vertex twice counts twice for it). They sum to the total area, which
    /// is that of the input mesh.
    ///
    /// [formula: lumped-mass]
    pub fn lumped_mass(&self) -> Vec<f64> {
        let mut mass = vec![0.0; self.n_vertices()];
        for (f, face) in self.faces.iter().enumerate() {
            let third = self.face_area(f) / 3.0;
            for &v in face {
                mass[v as usize] += third;
            }
        }
        mass
    }

    /// The lumped vertex areas divided by the total area, so that they sum to one (the area
    /// normalization of the repository's `CLAUDE.md`).
    ///
    /// [formula: lumped-mass]
    pub fn lumped_mass_normalized(&self) -> Vec<f64> {
        let mut mass = self.lumped_mass();
        let total = self.total_area();
        for m in &mut mass {
            *m /= total;
        }
        mass
    }

    /// For every vertex, the sum of the corner angles at it, computed from the edge lengths.
    /// Up to rounding it is the stored [`vertex_angle_sum`](Self::vertex_angle_sum): flips do
    /// not change the cone angle of a vertex.
    pub fn corner_angle_sums(&self) -> Vec<f64> {
        let mut sums = vec![0.0; self.n_vertices()];
        for h in 0..self.n_halfedges() {
            sums[self.from_vertex(h)] += self.corner_angle(h);
        }
        sums
    }

    /// The angle defect `Ω_i = 2π − Θ_i` of every interior vertex and `0.0` for boundary
    /// vertices, from the stored angle sums.
    ///
    /// [formula: angle-defect]
    pub fn angle_defects(&self) -> Vec<f64> {
        (0..self.n_vertices())
            .map(|v| {
                if self.is_boundary_vertex(v) {
                    0.0
                } else {
                    TAU - self.angle_sums[v]
                }
            })
            .collect()
    }

    /// The geodesic curvature `k_i = π − Θ_i` of every boundary vertex and `0.0` for interior
    /// vertices, from the stored angle sums.
    ///
    /// [formula: boundary-geodesic-curvature]
    pub fn boundary_geodesic_curvatures(&self) -> Vec<f64> {
        (0..self.n_vertices())
            .map(|v| {
                if self.is_boundary_vertex(v) {
                    PI - self.angle_sums[v]
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Angle and cotan-weight statistics. An interior edge is counted as having a negative
    /// weight when its weight is below `−weight_threshold`; pass `0.0` to count every negative
    /// weight, or a small positive number to leave out cocircular ties whose weight is
    /// negative by rounding.
    pub fn quality(&self, weight_threshold: f64) -> Quality {
        let mut min_angle = f64::INFINITY;
        let mut max_angle = 0.0f64;
        let mut largest: Vec<f64> = Vec::with_capacity(self.n_faces());
        for f in 0..self.n_faces() {
            let angles = [0, 1, 2].map(|k| self.corner_angle(3 * f + k));
            let big = angles[0].max(angles[1]).max(angles[2]);
            min_angle = min_angle.min(angles[0].min(angles[1]).min(angles[2]));
            max_angle = max_angle.max(big);
            largest.push(big);
        }
        largest.sort_unstable_by(f64::total_cmp);
        let quantile = |p: f64| {
            let n = largest.len();
            largest[(((n - 1) as f64) * p).round() as usize]
        };
        let above = |degrees: f64| {
            let limit = degrees.to_radians() + 1e-9;
            let below = largest.partition_point(|&a| a <= limit);
            (largest.len() - below) as f64 / largest.len() as f64
        };
        let mut n_interior_edges = 0;
        let mut n_negative_weights = 0;
        let mut min_cotan_weight = f64::INFINITY;
        for e in 0..self.n_edges() {
            if self.is_boundary_edge(e) {
                continue;
            }
            n_interior_edges += 1;
            let w = self.edge_cotan_weight(e);
            min_cotan_weight = min_cotan_weight.min(w);
            if w < -weight_threshold {
                n_negative_weights += 1;
            }
        }
        Quality {
            n_faces: self.n_faces(),
            n_interior_edges,
            min_angle,
            max_angle,
            largest_angle_quantiles: [quantile(0.5), quantile(0.9), quantile(0.99), quantile(1.0)],
            largest_angle_above: [above(90.0), above(120.0), above(150.0)],
            n_negative_weights,
            min_cotan_weight,
        }
    }
}
