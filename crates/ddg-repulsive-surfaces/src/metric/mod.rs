//! The fractional Sobolev metric `A = B + B⁰` and the comparison metrics (`SPEC.md` §4).
//!
//! Vertex functions `u ∈ R^{|V|}` enter the metric through two per-face operators (§4.1):
//!
//! ```text
//! ū(S)     = (u_i + u_j + u_k)/3                                   (U u = a ⊙ ū)
//! D_f u(S) = (1/(2a_S)) N_S × (u_i e_jk + u_j e_ki + u_k e_ij)       (V u = a ⊙ D_f u)
//! ```
//!
//! with the **un-normalized** edge vectors (D2: the paper's "unit edge vectors" is a misprint).
//! [`FaceOperators`] holds the per-face coefficients; [`dense`] assembles `L_τ`, `B`, `B⁰`, `A`
//! and the comparison metrics densely (small meshes, tests, the oracle and `Backend::Dense`);
//! `hierarchical` (WP6) multiplies by them through the block cluster tree.
//!
//! Everything here is orientation free: flipping a face flips `N_S` and the cyclic order of its
//! edges, so `D_f` is unchanged.

pub mod dense;
pub mod hierarchical;

pub use ddg_repulsive_curves::metric::DenseMatrix;

use crate::error::RsError;
use crate::geometry::{FaceGeometry, corners, opposite_edges};
use crate::surface::Surface;

/// Which metric the descent uses (`SPEC.md` §4.4, §4.5, §9.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Metric {
    /// The fractional Sobolev metric `A = w_TP (B + c_lo B⁰) + …` of §4.4; `low_order: false`
    /// drops `B⁰` (`c_lo = 0`, Sec. 9.2.4, the reference's `disable_near_field`).
    Sobolev {
        /// Keep the low-order term `B⁰`.
        low_order: bool,
    },
    /// The L² metric (§4.5).
    L2(L2Kind),
    /// `W ⊗ I₃` (cotan stiffness).
    H1,
    /// `W M⁻¹ W ⊗ I₃`.
    H2,
}

impl Default for Metric {
    /// `Sobolev { low_order: true }`.
    fn default() -> Self {
        Metric::Sobolev { low_order: true }
    }
}

/// The two L² metrics of §4.5 (D25).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum L2Kind {
    /// The lumped mass `M = diag(a_v)` (the paper).
    #[default]
    Mass,
    /// The identity (the reference).
    Identity,
}

/// The weights of the Sobolev metric `A = w_TP (B + c_lo B⁰) + w_W W M⁻¹ W (+ O_diag)`
/// (`SPEC.md` §4.4, D6, D7).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SobolevWeights {
    /// `w_TP`, the tangent-point weight (default 1).
    pub tp_weight: f64,
    /// `c_lo ∈ {0, 1}`: keep `B⁰` (default `true`).
    pub low_order: bool,
    /// `w_W`, the weight of the Willmore term `W M⁻¹ W` (0 = none, the default).
    pub willmore_weight: f64,
}

impl Default for SobolevWeights {
    fn default() -> Self {
        SobolevWeights {
            tp_weight: 1.0,
            low_order: true,
            willmore_weight: 0.0,
        }
    }
}

impl SobolevWeights {
    /// `c_lo` as a number.
    pub fn c_lo(&self) -> f64 {
        if self.low_order { 1.0 } else { 0.0 }
    }

    /// Checks the weights.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a non-finite or non-positive `tp_weight`, or a non-finite
    /// or negative `willmore_weight`.
    pub fn validate(&self) -> Result<(), RsError> {
        if !(self.tp_weight.is_finite() && self.tp_weight > 0.0) {
            return Err(RsError::InvalidOption {
                name: "tp_weight",
                value: self.tp_weight.to_string(),
                requirement: "a finite positive weight",
            });
        }
        if !(self.willmore_weight.is_finite() && self.willmore_weight >= 0.0) {
            return Err(RsError::InvalidOption {
                name: "willmore_weight",
                value: self.willmore_weight.to_string(),
                requirement: "a finite non-negative weight",
            });
        }
        Ok(())
    }
}

/// The coefficients of `D_f` on the triangle `p = [f_i, f_j, f_k]`: `g_c = (1/(2a)) N × e_opp(c)`,
/// so that `D_f u = Σ_c u_c g_c` (`SPEC.md` §4.1). `g_c` is the gradient of the hat function of
/// corner `c`: it is orthogonal to `N`, `g_c·(f_c − f_b) = 1` and `g_c·(f_b − f_d) = 0` for the
/// other two corners `b, d`; `Σ_c g_c = 0` (`D_f 1 = 0`). Un-normalized edges (D2).
///
/// [formula: rs-df-gradient]
#[inline]
pub fn df_coefficients(p: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    use crate::geometry::v3::{cross, norm, scale, sub};
    let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
    // N/(2a) = n/|n|², since 2a = |n|.
    let len = norm(n);
    let w = scale(n, 1.0 / (len * len));
    let e = opposite_edges(p);
    [cross(w, e[0]), cross(w, e[1]), cross(w, e[2])]
}

/// `D_f u(S)` on one triangle from its three corner values (`SPEC.md` §4.1).
///
/// [formula: rs-df-gradient]
#[inline]
pub fn df_face(p: [[f64; 3]; 3], u: [f64; 3]) -> [f64; 3] {
    let g = df_coefficients(p);
    let mut out = [0.0; 3];
    for c in 0..3 {
        for r in 0..3 {
            out[r] += u[c] * g[c][r];
        }
    }
    out
}

/// Per-face data of the operators `U` and `V` (`SPEC.md` §4.1): areas, barycenters, normals and
/// the `D_f` coefficients ([`df_coefficients`]) of every face.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceOperators {
    /// `a_S`.
    pub areas: Vec<f64>,
    /// `X_S`.
    pub barycenters: Vec<[f64; 3]>,
    /// `N_S`.
    pub normals: Vec<[f64; 3]>,
    /// `g_c(S)` for the corners `c = 0, 1, 2` of face `S` (in the face's vertex order).
    pub df: Vec<[[f64; 3]; 3]>,
    /// The faces (copied, so the operators can be applied without the surface).
    pub faces: Vec<[u32; 3]>,
    /// `|V|`.
    pub n_vertices: usize,
}

impl FaceOperators {
    /// Computes the per-face data after [`Surface::validate_positions`].
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`].
    pub fn new(surface: &Surface, positions: &[[f64; 3]]) -> Result<FaceOperators, RsError> {
        let geom = FaceGeometry::new(surface, positions)?;
        Ok(Self::from_geometry(surface, positions, geom))
    }

    /// As [`FaceOperators::new`] from already computed (validated) face geometry.
    pub fn from_geometry(
        surface: &Surface,
        positions: &[[f64; 3]],
        geom: FaceGeometry,
    ) -> FaceOperators {
        let df = surface
            .faces()
            .iter()
            .map(|t| df_coefficients(corners(positions, t)))
            .collect();
        FaceOperators {
            areas: geom.areas,
            barycenters: geom.barycenters,
            normals: geom.normals,
            df,
            faces: surface.faces().to_vec(),
            n_vertices: surface.n_vertices(),
        }
    }

    /// `|F|`.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// The face averages `ū_S` of a vertex function.
    pub fn average(&self, u: &[f64]) -> Vec<f64> {
        self.faces
            .iter()
            .map(|t| (u[t[0] as usize] + u[t[1] as usize] + u[t[2] as usize]) / 3.0)
            .collect()
    }

    /// `D_f u(S)` for every face (`SPEC.md` §4.1).
    ///
    /// [formula: rs-df-gradient]
    pub fn gradient(&self, u: &[f64]) -> Vec<[f64; 3]> {
        self.faces
            .iter()
            .zip(&self.df)
            .map(|(t, g)| {
                let mut out = [0.0; 3];
                for c in 0..3 {
                    let uc = u[t[c] as usize];
                    for r in 0..3 {
                        out[r] += uc * g[c][r];
                    }
                }
                out
            })
            .collect()
    }

    /// The sparse row of `X = U` (`channel = Average`: one channel, coefficients `a_S/3`) or
    /// `X = V` (`Gradient`: three channels, `a_S g_c(S)`) of face `s`: for each corner `c`, the
    /// vertex and the (up to 3) channel coefficients.
    #[inline]
    pub fn row(&self, s: usize, channel: Channel) -> [(usize, [f64; 3]); 3] {
        let t = self.faces[s];
        let a = self.areas[s];
        let mut out = [(0usize, [0.0; 3]); 3];
        for c in 0..3 {
            out[c].0 = t[c] as usize;
            match channel {
                Channel::Average => out[c].1 = [a / 3.0, 0.0, 0.0],
                Channel::Gradient => {
                    let g = self.df[s][c];
                    out[c].1 = [a * g[0], a * g[1], a * g[2]];
                }
            }
        }
        out
    }
}

/// Which per-face operator a kernel acts through (`SPEC.md` §4.3): `U` (face averages, `L_τ`,
/// `B⁰`) or `V` (face gradients `D_f`, `B`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Channel {
    /// `X = U = diag(a) · averaging`; one channel.
    Average,
    /// `X = V = diag(a) D_f`; three channels (the components of `D_f u`).
    Gradient,
}

impl Channel {
    /// Number of channels (1 or 3).
    pub fn width(self) -> usize {
        match self {
            Channel::Average => 1,
            Channel::Gradient => 3,
        }
    }
}

/// The hook for the mesh-obstacle term `O_diag` of §4.4 (D23; filled by WP6/WP7): per face,
/// `d_S = Σ_𝒪 w_𝒪 Σ_{T∈𝒪} h(S, T) a_T` for the high-order kernel (`hi`, acting through `V`)
/// and the low-order kernel (`lo`, through `U`). The metric adds
/// `2 Vᵀ diag(hi/a) V + 2 Uᵀ diag(lo/a) U` (the diagonal part of the kernel-matrix form of
/// §4.3 for mesh ∪ obstacle with the obstacle fixed), unscaled by `w_TP`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FaceDiagonal {
    /// `d_S` of the high-order kernel (empty = none).
    pub hi: Vec<f64>,
    /// `d_S` of the low-order kernel (empty = none).
    pub lo: Vec<f64>,
}
