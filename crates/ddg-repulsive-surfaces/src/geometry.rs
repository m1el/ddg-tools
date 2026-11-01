//! Per-face geometry, its derivatives, and the orientation-free cotan stiffness and lumped mass
//! (`SPEC.md` §2.1, §3.2).
//!
//! For a face `S = (i, j, k)` with positions `f_i, f_j, f_k`:
//!
//! ```text
//! n_S = (f_j − f_i) × (f_k − f_i),   a_S = ½|n_S|,   N_S = n_S/|n_S|,   X_S = (f_i + f_j + f_k)/3
//! e_jk = f_k − f_j,   e_ki = f_i − f_k,   e_ij = f_j − f_i          (not normalized, D2)
//! ∂a_S/∂f_i = ½ N_S × e_jk,     ∂N_S/∂f_i = (1/(2a_S)) (e_jk × N_S) N_Sᵀ
//! ```
//!
//! (cyclically for `j`, `k`). The vertex area is `a_v = ⅓ Σ_{S∋v} a_S`. Flipping a face's
//! orientation flips `N_S` and leaves `a_S`, `X_S`, `a_v`, `W` and `M` unchanged.

use ddg_core::SparseMatrix;

use crate::error::RsError;
use crate::surface::Surface;

/// Small 3-vector helpers shared by the modules of the crate.
pub mod v3 {
    pub use ddg_repulsive_curves::kernel::v3::{add, axpy, cross, dot, norm, scale, sub};
}

use v3::{add, cross, dot, norm, scale, sub};

/// `(a_S, N_S, X_S)` of the triangle `p = [f_i, f_j, f_k]`: area, unit normal, barycenter.
/// A degenerate triangle gives area 0 and a non-finite normal (callers validate first, see
/// [`Surface::validate_positions`]).
///
/// [formula: rs-face-quantities]
#[inline]
pub fn face_quantities(p: [[f64; 3]; 3]) -> (f64, [f64; 3], [f64; 3]) {
    let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
    let len = norm(n);
    let x = scale(add(add(p[0], p[1]), p[2]), 1.0 / 3.0);
    (0.5 * len, scale(n, 1.0 / len), x)
}

/// The edge vectors opposite the three corners: `[e_jk, e_ki, e_ij]` (not normalized, D2).
#[inline]
pub fn opposite_edges(p: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [sub(p[2], p[1]), sub(p[0], p[2]), sub(p[1], p[0])]
}

/// `∂a_S/∂f_c = ½ N_S × e_opp(c)` for the three corners `c` of the triangle `p`.
///
/// [formula: rs-face-area-gradient]
#[inline]
pub fn face_area_gradient(p: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let (_, nrm, _) = face_quantities(p);
    let e = opposite_edges(p);
    [
        scale(cross(nrm, e[0]), 0.5),
        scale(cross(nrm, e[1]), 0.5),
        scale(cross(nrm, e[2]), 0.5),
    ]
}

/// The Jacobians `J_c = ∂N_S/∂f_c = (1/(2a_S)) (e_opp(c) × N_S) N_Sᵀ` of the unit normal with
/// respect to the three corners, as row-major matrices: `J_c[r][s] = ∂(N_S)_r / ∂(f_c)_s`.
///
/// [formula: rs-face-normal-variation]
pub fn face_normal_variation(p: [[f64; 3]; 3]) -> [[[f64; 3]; 3]; 3] {
    let (a, nrm, _) = face_quantities(p);
    let e = opposite_edges(p);
    let mut out = [[[0.0; 3]; 3]; 3];
    for c in 0..3 {
        let w = scale(cross(e[c], nrm), 1.0 / (2.0 * a));
        for r in 0..3 {
            for s in 0..3 {
                out[c][r][s] = w[r] * nrm[s];
            }
        }
    }
    out
}

/// Pull-back of the per-face partials of a function `F(a_S, X_S, N_S)` to the corners of the
/// triangle `p`: `∂F/∂f_c = (∂F/∂a) ∂a/∂f_c + ⅓ ∂F/∂X + J_cᵀ ∂F/∂N`, with
/// `J_cᵀ g = N_S ((e_opp(c) × N_S)·g)/(2a_S)` (the transposes of
/// [`face_normal_variation`] and the gradient of [`face_area_gradient`]).
///
/// [formula: rs-face-area-gradient]
/// [formula: rs-face-normal-variation]
#[inline]
pub fn face_pull_back(
    p: [[f64; 3]; 3],
    d_area: f64,
    d_bary: [f64; 3],
    d_normal: [f64; 3],
) -> [[f64; 3]; 3] {
    let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
    let len = norm(n);
    let nrm = scale(n, 1.0 / len);
    let e = opposite_edges(p);
    let xb = scale(d_bary, 1.0 / 3.0);
    // 2a = |n|.
    let inv_2a = 1.0 / len;
    let mut out = [[0.0; 3]; 3];
    for c in 0..3 {
        let ga = scale(cross(nrm, e[c]), 0.5 * d_area);
        let gn = scale(nrm, dot(cross(e[c], nrm), d_normal) * inv_2a);
        out[c] = add(add(ga, xb), gn);
    }
    out
}

/// Per-face areas, unit normals and barycenters of a surface (`SPEC.md` §2.1).
#[derive(Debug, Clone, PartialEq)]
pub struct FaceGeometry {
    /// `a_S`.
    pub areas: Vec<f64>,
    /// `N_S` (sign from the caller's vertex order; irrelevant to the energy).
    pub normals: Vec<[f64; 3]>,
    /// `X_S`.
    pub barycenters: Vec<[f64; 3]>,
}

impl FaceGeometry {
    /// Computes the face quantities after [`Surface::validate_positions`].
    ///
    /// [formula: rs-face-quantities]
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`].
    pub fn new(surface: &Surface, positions: &[[f64; 3]]) -> Result<FaceGeometry, RsError> {
        surface.validate_positions(positions)?;
        Ok(Self::new_unchecked(surface, positions))
    }

    /// As [`FaceGeometry::new`] without validation (positions must have been validated; a
    /// degenerate face gives a non-finite normal).
    pub fn new_unchecked(surface: &Surface, positions: &[[f64; 3]]) -> FaceGeometry {
        let nf = surface.n_faces();
        let mut g = FaceGeometry {
            areas: Vec::with_capacity(nf),
            normals: Vec::with_capacity(nf),
            barycenters: Vec::with_capacity(nf),
        };
        for t in surface.faces() {
            let (a, n, x) = face_quantities(corners(positions, t));
            g.areas.push(a);
            g.normals.push(n);
            g.barycenters.push(x);
        }
        g
    }

    /// `Σ_S a_S`.
    pub fn total_area(&self) -> f64 {
        self.areas.iter().sum()
    }
}

/// The three corner positions of face `t`.
#[inline]
pub fn corners(positions: &[[f64; 3]], t: &[u32; 3]) -> [[f64; 3]; 3] {
    [
        positions[t[0] as usize],
        positions[t[1] as usize],
        positions[t[2] as usize],
    ]
}

/// Pull-back of per-face partials `(∂F/∂a_S, ∂F/∂X_S, ∂F/∂N_S)` to the vertices: the
/// differential `∂F/∂f_v` (see [`face_pull_back`]). Faces are added in face order, so the
/// result is deterministic.
///
/// [formula: rs-face-area-gradient]
/// [formula: rs-face-normal-variation]
pub fn pull_back(
    surface: &Surface,
    positions: &[[f64; 3]],
    d_area: &[f64],
    d_bary: &[[f64; 3]],
    d_normal: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let mut out = vec![[0.0; 3]; surface.n_vertices()];
    for (f, t) in surface.faces().iter().enumerate() {
        let g = face_pull_back(corners(positions, t), d_area[f], d_bary[f], d_normal[f]);
        for c in 0..3 {
            let v = t[c] as usize;
            out[v] = add(out[v], g[c]);
        }
    }
    out
}

/// Barycentric dual vertex areas `a_v = ⅓ Σ_{S∋v} a_S` from the face areas (the diagonal of
/// the lumped mass `M`). They sum to the total area.
///
/// [formula: rs-face-quantities]
pub fn vertex_areas(surface: &Surface, face_areas: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; surface.n_vertices()];
    for (f, t) in surface.faces().iter().enumerate() {
        let a = face_areas[f] / 3.0;
        for &v in t {
            out[v as usize] += a;
        }
    }
    out
}

/// The lumped mass `M = diag(a_v)` (`SPEC.md` §2.1; un-normalized: `Σ a_v` is the total area).
///
/// # Errors
/// The errors of [`Surface::validate_positions`].
pub fn lumped_mass(surface: &Surface, positions: &[[f64; 3]]) -> Result<Vec<f64>, RsError> {
    let g = FaceGeometry::new(surface, positions)?;
    Ok(vertex_areas(surface, &g.areas))
}

/// The cotan stiffness matrix `W` (`|V| × |V|`, symmetric positive semidefinite, `CLAUDE.md`
/// convention): `W_ij = −½ Σ cot(angle opposite ij)` over the one or two faces of edge `ij`,
/// `W_ii = −Σ_{j≠i} W_ij`. Computed from the face list, so it does not depend on the faces'
/// orientation; equal to `ddg_core::operators::cotan_laplacian` on oriented meshes.
///
/// # Errors
/// The errors of [`Surface::validate_positions`].
pub fn cotan_stiffness(surface: &Surface, positions: &[[f64; 3]]) -> Result<SparseMatrix, RsError> {
    surface.validate_positions(positions)?;
    let n = surface.n_vertices();
    let mut trip: Vec<(u32, u32, f64)> = Vec::with_capacity(12 * surface.n_faces());
    for t in surface.faces() {
        let p = corners(positions, t);
        for c in 0..3 {
            // Angle at corner c, opposite the edge (c+1, c+2).
            let a = t[(c + 1) % 3];
            let b = t[(c + 2) % 3];
            let u = sub(p[(c + 1) % 3], p[c]);
            let w = sub(p[(c + 2) % 3], p[c]);
            let cot = dot(u, w) / norm(cross(u, w));
            let h = 0.5 * cot;
            trip.push((a, b, -h));
            trip.push((b, a, -h));
            trip.push((a, a, h));
            trip.push((b, b, h));
        }
    }
    SparseMatrix::from_triplets(n, n, &trip).map_err(|_| RsError::TooLarge {
        n_vertices: n,
        n_faces: surface.n_faces(),
    })
}

/// Mean edge length `ℓ̄` (the length scale of finite differences, line-search starts and
/// remeshing).
pub fn mean_edge_length(surface: &Surface, positions: &[[f64; 3]]) -> f64 {
    let sum: f64 = surface
        .edges()
        .iter()
        .map(|&[a, b]| norm(sub(positions[a as usize], positions[b as usize])))
        .sum();
    sum / surface.n_edges() as f64
}
