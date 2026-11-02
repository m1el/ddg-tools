//! The Willmore term `fᵀ W M⁻¹ W f` (`SPEC.md` §5.2, Sec. 7.3.5) with its exact gradient, and the
//! sparse matrix `W M⁻¹ W` for the metric.
//!
//! With `r_i = (W f)_i ∈ R³` and `m_i = a_v` (a third of the incident face areas), the value is
//! `P = Σ_{i interior} |r_i|² / m_i`: rows of boundary vertices are dropped (as the reference).
//! Since `|(W f)_i| ≈ m_i |Δf| = 2 m_i |H|`, on a smooth closed surface `P ≈ ∫ |Δf|² dA =
//! 4 ∫ H² dA` (`H` the mean curvature); the unit sphere has `P ≈ 16π`. (The reference's value is
//! `¼` of this.) `P` is invariant under scaling and rigid motions.
//!
//! The gradient, with `y_i = r_i / m_i` (zero on boundary rows):
//!
//! ```text
//! dP = 2 Σ_i y_i·dr_i − Σ_i |y_i|² dm_i,     dr_i = Σ_j (dW_ij f_j + W_ij df_j)
//! ∂P/∂f_j  = 2 (W y)_j + (terms of dW) − ...
//! ```
//!
//! The variation of `W` goes through the cotangents: `Σ_ij dW_ij y_i·f_j = Σ_{faces, corners c}
//! ½ dcot_c (y_a − y_b)·(f_a − f_b)` (`a, b` the other two corners), and
//! `∂m_i/∂f = ⅓ Σ_{S∋i} ∂a_S/∂f`.
//!
//! [formula: rs-willmore]

use ddg_core::SparseMatrix;

use crate::error::RsError;
use crate::geometry::v3::{add, axpy, cross, dot, norm, scale, sub};
use crate::geometry::{corners, cotan_stiffness, face_area_gradient, vertex_areas, FaceGeometry};
use crate::surface::Surface;

/// `W f` row by row (vectors) from the cotan weights, and the vertex areas.
fn apply_w(
    surface: &Surface,
    positions: &[[f64; 3]],
    x: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let mut out = vec![[0.0; 3]; surface.n_vertices()];
    for t in surface.faces() {
        let p = corners(positions, t);
        for c in 0..3 {
            let (a, b) = (t[(c + 1) % 3] as usize, t[(c + 2) % 3] as usize);
            let u = sub(p[(c + 1) % 3], p[c]);
            let w = sub(p[(c + 2) % 3], p[c]);
            let h = 0.5 * dot(u, w) / norm(cross(u, w));
            let d = sub(x[a], x[b]);
            out[a] = axpy(out[a], h, d);
            out[b] = axpy(out[b], -h, d);
        }
    }
    out
}

/// The Willmore value `P` and its exact gradient.
///
/// # Errors
/// The errors of [`Surface::validate_positions`].
pub fn willmore_and_gradient(
    surface: &Surface,
    positions: &[[f64; 3]],
) -> Result<(f64, Vec<[f64; 3]>), RsError> {
    let geom = FaceGeometry::new(surface, positions)?;
    let m = vertex_areas(surface, &geom.areas);
    let n = surface.n_vertices();
    let r = apply_w(surface, positions, positions);
    let mut y = vec![[0.0; 3]; n];
    let mut value = 0.0;
    for i in 0..n {
        if !surface.is_boundary_vertex(i as u32) {
            y[i] = scale(r[i], 1.0 / m[i]);
            value += dot(r[i], y[i]);
        }
    }
    // 2 W y (W symmetric).
    let wy = apply_w(surface, positions, &y);
    let mut g: Vec<[f64; 3]> = wy.iter().map(|v| scale(*v, 2.0)).collect();
    for t in surface.faces() {
        let p = corners(positions, t);
        // Mass variation: −Σ_{i ∈ S} |y_i|² / 3 · ∂a_S.
        let ya: f64 = t.iter().map(|&v| dot(y[v as usize], y[v as usize])).sum();
        let ga = face_area_gradient(p);
        for c in 0..3 {
            g[t[c] as usize] = axpy(g[t[c] as usize], -ya / 3.0, ga[c]);
        }
        // Cotangent variation.
        for c in 0..3 {
            let (ia, ib) = ((c + 1) % 3, (c + 2) % 3);
            let (a, b) = (t[ia] as usize, t[ib] as usize);
            let coef = dot(sub(y[a], y[b]), sub(positions[a], positions[b]));
            let u = sub(p[ia], p[c]);
            let w = sub(p[ib], p[c]);
            let cr = cross(u, w);
            let s = norm(cr);
            let nh = scale(cr, 1.0 / s);
            let cot = dot(u, w) / s;
            // ∂cot/∂u, ∂cot/∂w.
            let du = sub(scale(w, 1.0 / s), scale(cross(w, nh), cot / s));
            let dw = sub(scale(u, 1.0 / s), scale(cross(nh, u), cot / s));
            g[a] = axpy(g[a], coef, du);
            g[b] = axpy(g[b], coef, dw);
            let vc = t[c] as usize;
            g[vc] = axpy(g[vc], -coef, add(du, dw));
        }
    }
    Ok((value, g))
}

/// The sparse symmetric PSD matrix `W M⁻¹ W` (`|V| × |V|`, to be used as `⊗ I₃`): with
/// `interior_only` the rows of boundary vertices of the first `W` are dropped, i.e.
/// `Σ_{i interior} w_i w_iᵀ / m_i` (the Hessian-like term of the penalty, `SPEC.md` §4.4);
/// otherwise the full `W M⁻¹ W` (the `H²` comparison metric, `SPEC.md` §4.5).
///
/// # Errors
/// The errors of [`Surface::validate_positions`]; [`RsError::TooLarge`] from the sparse
/// assembly.
pub fn willmore_matrix(
    surface: &Surface,
    positions: &[[f64; 3]],
    interior_only: bool,
) -> Result<SparseMatrix, RsError> {
    let geom = FaceGeometry::new(surface, positions)?;
    let m = vertex_areas(surface, &geom.areas);
    let w = cotan_stiffness(surface, positions)?;
    let n = surface.n_vertices();
    // W is symmetric: its column i is the row i.
    let mut trip: Vec<(u32, u32, f64)> = Vec::new();
    for i in 0..n {
        if interior_only && surface.is_boundary_vertex(i as u32) {
            continue;
        }
        let col: Vec<(u32, f64)> = (w.col_ptr()[i]..w.col_ptr()[i + 1])
            .map(|k| (w.row_idx()[k], w.values()[k]))
            .collect();
        for &(a, wa) in &col {
            for &(b, wb) in &col {
                trip.push((a, b, wa * wb / m[i]));
            }
        }
    }
    SparseMatrix::from_triplets(n, n, &trip).map_err(|_| RsError::TooLarge {
        n_vertices: n,
        n_faces: surface.n_faces(),
    })
}
