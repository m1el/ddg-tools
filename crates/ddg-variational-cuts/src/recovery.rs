//! Recovery of the normal derivative on the cut from the integrated normal derivatives
//! (SPEC §3.5, D12).
//!
//! The integrated normal derivative `v_p = (L u)_p + g_p` of a cut vertex is accurate when it is
//! summed over a stretch of the cut that is as long as the triangles are wide, but not per cut
//! vertex: the hat function it is tested with reaches one mesh width into the region, so its
//! error is of the order of the mesh width, whatever the length `ℓ_p` of the cut that belongs to
//! the vertex. Where consecutive crossings crowd around one mesh vertex, `ℓ_p` is as small as
//! their distance and `v_p / ℓ_p` is meaningless. [`FluxRecovery::Smoothed`] therefore smooths
//! the fluxes along the cut, implicitly and with the mass `ℓ`, before dividing.

use ddg_core::SparseMatrix;
use ddg_core::sparse::{Factorization, Method, Parallelism, Symbolic};

use crate::cut_mesh::{BaseGeometry, CutTriangles, EdgeKind};
use crate::error::CutError;
use crate::problem::FluxRecovery;

/// The coupling `κ_σ = (c · h_f)² / |σ|` of the two ends of a cut segment `σ` in the smoothed
/// flux recovery (SPEC §3.5): `c` is the smoothing length in units of the local edge length,
/// `h_f²` the mean squared edge length of the face that contains the segment, `|σ|` the length
/// of the segment. It is the stiffness of a 1D linear element of length `|σ|`, times the local
/// smoothing time `(c h_f)²`.
///
/// [formula: cut-flux-recovery]
#[inline]
pub fn flux_smoothing_stiffness(length: f64, mean_squared_edge: f64, segment_length: f64) -> f64 {
    length * length * mean_squared_edge / segment_length
}

/// Maps integrated normal derivatives at the boundary vertices of the cut mesh to normal
/// derivatives.
pub(crate) struct FluxSolver {
    n_original: usize,
    /// The factorization of `diag(ℓ) + S` over the copies of the cut points; `None` for the
    /// nodal recovery or when there is no cut.
    factorization: Option<Factorization>,
}

impl FluxSolver {
    /// Sets up the recovery for the cut mesh `tri` with dual lengths `dual_length`: for
    /// [`FluxRecovery::Smoothed`] the matrix
    ///
    /// ```text
    /// diag(ℓ_p) + Σ_{cut segments σ = pq of each region} κ_σ [[1, −1], [−1, 1]]
    /// ```
    ///
    /// over the cut vertices of the cut mesh (each cut point once per side; the two sides are
    /// not coupled) is assembled and factorized. It is symmetric positive definite, an
    /// M-matrix, and has the row sums `ℓ_p`.
    ///
    /// [formula: cut-flux-recovery]
    pub(crate) fn new(
        tri: &CutTriangles,
        base: &BaseGeometry,
        dual_length: &[f64],
        recovery: FluxRecovery,
        par: Parallelism<'_>,
    ) -> Result<FluxSolver, CutError> {
        let nv = tri.n_original;
        let n = 2 * tri.n_points;
        let length = match recovery {
            FluxRecovery::Smoothed { length } if length > 0.0 && n > 0 => length,
            _ => {
                return Ok(FluxSolver {
                    n_original: nv,
                    factorization: None,
                });
            }
        };
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(4 * n);
        for (i, &dual) in dual_length[nv..].iter().enumerate() {
            triplets.push((i as u32, i as u32, dual));
        }
        for (t, kinds) in tri.edge_kind.iter().enumerate() {
            for k in 0..3 {
                if kinds[k] != EdgeKind::Cut {
                    continue;
                }
                let f = tri.parent[t] as usize;
                let mean_squared_edge = (base.len2[f][0] + base.len2[f][1] + base.len2[f][2]) / 3.0;
                let kappa =
                    flux_smoothing_stiffness(length, mean_squared_edge, tri.len2[t][k].sqrt());
                let i = tri.faces[t][k] - nv as u32;
                let j = tri.faces[t][(k + 1) % 3] - nv as u32;
                triplets.push((i, i, kappa));
                triplets.push((j, j, kappa));
                triplets.push((i.max(j), i.min(j), -kappa));
            }
        }
        let matrix = SparseMatrix::from_triplets(n, n, &triplets)?;
        let factorization = Symbolic::analyze(&matrix, Method::Llt)?.factor(&matrix, par)?;
        Ok(FluxSolver {
            n_original: nv,
            factorization: Some(factorization),
        })
    }

    /// The normal derivative at every vertex of the cut mesh from the integrated normal
    /// derivative `flux` (`v` or `q`): at the copies of the cut points by the recovery this
    /// solver was set up for, at the other boundary vertices `flux / dual_length`, and zero at
    /// interior vertices.
    ///
    /// [formula: cut-flux-recovery]
    pub(crate) fn recover(
        &self,
        flux: &[f64],
        dual_length: &[f64],
        par: Parallelism<'_>,
    ) -> Result<Vec<f64>, CutError> {
        let mut out: Vec<f64> = flux
            .iter()
            .zip(dual_length)
            .map(|(&v, &dual)| if dual > 0.0 { v / dual } else { 0.0 })
            .collect();
        if let Some(factorization) = &self.factorization {
            let cut = &mut out[self.n_original..];
            cut.copy_from_slice(&flux[self.n_original..]);
            factorization.solve_in_place(cut, par)?;
        }
        Ok(out)
    }
}
