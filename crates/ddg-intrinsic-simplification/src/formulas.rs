//! The closed-form formulas of intrinsic simplification, as plain functions of lengths, angles
//! and barycentric coordinates. Each carries a `[formula: iem-…]` tag (SPEC §13). Areas,
//! angles and cotangents of single triangles come from [`ddg_core::geometry`]; flips and
//! the Delaunay test from [`ddg_intrinsic::formulas`].

use std::f64::consts::{PI, TAU};

use ddg_core::geometry::{corner_angle, triangle_area};

/// An angular coordinate brought to its range (SPEC S4): `[0, 2π)` at an interior vertex,
/// unchanged at a boundary vertex (where the coordinates lie in `[0, π]` by construction and
/// do not wrap). One function for every place that stores a coordinate.
///
/// [formula: iem-wrap]
#[inline]
pub fn wrap_angle(phi: f64, boundary: bool) -> f64 {
    if boundary {
        return phi;
    }
    let r = phi.rem_euclid(TAU);
    // rem_euclid can round up to TAU for a tiny negative input.
    if r >= TAU { 0.0 } else { r }
}

/// The corner angle opposite the side `opposite` in a triangle with sides `opposite`,
/// `side1`, `side2`, *extended* to lengths that violate the triangle inequality (Springborn,
/// Schröder, Pinkall 2008, Sec. 3.1): `π` when `opposite ≥ side1 + side2`, `0` when another
/// side is at least the sum of the two others, the law-of-cosines angle otherwise. With it
/// the angle sum at a vertex is a continuous function of the conformal factor
/// (`iem-flatten-monotone`). Returns `None` for a non-positive or non-finite length.
///
/// [formula: iem-flatten-monotone]
pub fn extended_corner_angle(opposite: f64, side1: f64, side2: f64) -> Option<f64> {
    let ok = |x: f64| x > 0.0 && x.is_finite();
    if !(ok(opposite) && ok(side1) && ok(side2)) {
        return None;
    }
    if opposite >= side1 + side2 {
        return Some(PI);
    }
    if side1 >= opposite + side2 || side2 >= opposite + side1 {
        return Some(0.0);
    }
    // Strictly inside the triangle inequality: an angle in (0, π), or the limit value where
    // rounding leaves no area.
    Some(corner_angle(opposite, side1, side2).unwrap_or_else(|| {
        let longest = opposite.max(side1).max(side2);
        if opposite == longest { PI } else { 0.0 }
    }))
}

/// The angles of a face at its three slots (`angles[k]` at slot `k`, opposite the side
/// `l[(k + 1) % 3]`, where `l[k]` is the side from slot `k` to slot `k + 1`), extended as in
/// [`extended_corner_angle`], and whether the face is a valid triangle (strict triangle
/// inequality and positive area).
///
/// [formula: iem-flatten-monotone]
pub fn extended_face_angles(l: [f64; 3]) -> Option<([f64; 3], bool)> {
    let valid = triangle_area(l[0], l[1], l[2]).is_some();
    let mut a = [0.0; 3];
    for (k, slot) in a.iter_mut().enumerate() {
        *slot = extended_corner_angle(l[(k + 1) % 3], l[(k + 2) % 3], l[k])?;
    }
    Some((a, valid))
}

/// The contribution of one face to the angle sum `Θ_i` of a vertex `i` and to its derivative
/// in the conformal factor `u` of `i` (paper Eq. 3, 4; Springborn et al. 2008 Eq. 9, 10
/// restricted to one vertex), for a face whose slots `at_i[k]` are `i`. With `l` the scaled
/// side lengths (`l[k]` from slot `k` to `k + 1`): the sum of the extended angles at the slots
/// of `i`, and `Σ_{a, s ∈ slots of i} ∂θ_a/∂u_s` with `∂θ_a/∂u_a = −½ (cot θ_b + cot θ_c)`
/// and `∂θ_a/∂u_b = ½ cot θ_c` (`c` the third slot). For a face with one slot at `i` this is
/// Eq. 4's `−½ (cot θ^k_ij + cot θ^j_ki)`. A face that violates the triangle inequality
/// contributes `0` to the derivative (its extended angles are locally constant). Returns
/// `None` for a non-positive or non-finite length.
///
/// [formula: iem-flatten-newton]
pub fn face_angle_sum_and_derivative(l: [f64; 3], at_i: [bool; 3]) -> Option<(f64, f64, bool)> {
    // Only the angles at the slots of i are needed (bitwise those of `extended_face_angles`).
    let mut sum = 0.0;
    for k in 0..3 {
        if at_i[k] {
            sum += extended_corner_angle(l[(k + 1) % 3], l[(k + 2) % 3], l[k])?;
        }
    }
    let Some(area) = triangle_area(l[0], l[1], l[2]) else {
        return Some((sum, 0.0, false));
    };
    // cot at slot k (opposite side l[(k + 1) % 3]), with the one area of the face: bitwise
    // `ddg_core::geometry::corner_cotan`.
    let mut cot = [0.0; 3];
    for (k, c) in cot.iter_mut().enumerate() {
        let (a, s1, s2) = (l[(k + 1) % 3], l[(k + 2) % 3], l[k]);
        let (b, cc) = if s1 >= s2 { (s1, s2) } else { (s2, s1) };
        let x = (cc * cc + (b - a) * (b + a)) / (4.0 * area);
        *c = if x.is_finite() { x } else { 0.0 };
    }
    let mut d = 0.0;
    for a in 0..3 {
        if !at_i[a] {
            continue;
        }
        let (b, c) = ((a + 1) % 3, (a + 2) % 3);
        d -= 0.5 * (cot[b] + cot[c]);
        if at_i[b] {
            d += 0.5 * cot[c];
        }
        if at_i[c] {
            d += 0.5 * cot[b];
        }
    }
    Some((sum, d, true))
}

/// Which sign of the exponent the projective interpolation of Eq. 10 uses (SPEC S25).
///
/// `Negative` (default): `b̃_k ∝ e^{−u} b_k` at the slots of the flattened vertex, the
/// circumcircle-preserving projective interpolation of Springborn, Schröder, Pinkall 2008
/// (Sec. 3.4, App. D), which the paper cites for this step. `Positive`: `e^{+u}`, Eq. 10 as the
/// paper prints it and as the reference implements it — the inverse of that map. Both are continuous
/// bijections of every face onto itself ([`projective_interpolation_signed`]), so the tracking
/// is correct with either; they differ in where interior points land (`Negative`: much lower
/// anisotropic distortion of the map; `Positive`: 3–19 % lower area distortion; SPEC §10.5).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ProjectiveSign {
    /// `e^{+u}` (the paper's Eq. 10 as printed and the reference's).
    Positive,
    /// `e^{−u}` (Springborn et al. 2008; default).
    #[default]
    Negative,
}

/// Eq. 10 of the paper, as printed: the barycentric coordinates of a point of a face after the
/// vertex `i` with conformal factor `u` is flattened: `b̃_k ∝ e^{u} b_k` at the slots of `i`,
/// `b_k` elsewhere, normalized to sum one. It maps the face bijectively onto itself, fixes its
/// corners and agrees along a shared side with the neighbouring face. It is the *inverse* of
/// the circumcircle-preserving projective interpolation of Springborn et al. 2008 (Sec. 3.4,
/// App. D), which uses `e^{−u}`; both are valid tracking maps (SPEC S25, option
/// [`ProjectiveSign`]).
///
/// [formula: iem-projective-interpolation]
#[inline]
pub fn projective_interpolation(b: [f64; 3], at_i: [bool; 3], u: f64) -> [f64; 3] {
    projective_interpolation_signed(b, at_i, u, ProjectiveSign::Positive)
}

/// [`projective_interpolation`] with the sign of the exponent chosen by `sign` (SPEC S25).
///
/// [formula: iem-projective-interpolation]
#[inline]
pub fn projective_interpolation_signed(
    b: [f64; 3],
    at_i: [bool; 3],
    u: f64,
    sign: ProjectiveSign,
) -> [f64; 3] {
    let s = match sign {
        ProjectiveSign::Positive => u.exp(),
        ProjectiveSign::Negative => (-u).exp(),
    };
    let w = [0, 1, 2].map(|k| if at_i[k] { s * b[k] } else { b[k] });
    let sum = w[0] + w[1] + w[2];
    w.map(|x| x / sum)
}

/// The barycentric coordinates of the removed vertex `i` of a flat star of degree three in the
/// merged face `x₀ x₁ x₂` (paper Sec. 7.1, "standard formulas"): the coordinate of corner `m`
/// is the area of the sub-triangle `i x_{m+1} x_{m+2}` over the sum of the three, each area
/// from its lengths (`areas[m]` is that of the face of the star opposite `x_m`). Exact for a
/// flat star, where the three faces tile the merged one. `None` if an area is not positive.
///
/// [formula: iem-degree-three-barycentric]
pub fn degree_three_barycentric(areas: [f64; 3]) -> Option<[f64; 3]> {
    if !areas.iter().all(|&a| a > 0.0 && a.is_finite()) {
        return None;
    }
    let sum = areas[0] + areas[1] + areas[2];
    Some(areas.map(|a| a / sum))
}

/// Where a flat boundary vertex `i` of degree two lies on the merged boundary edge `bc` (paper
/// Sec. 4.2): the parameter `t = ℓ_bi / (ℓ_bi + ℓ_ic)` from `b`, so that `i = (1 − t) b + t c`,
/// and the merged length `ℓ_bc = ℓ_bi + ℓ_ic`. Exact when the corner angles at `i` sum to `π`.
///
/// [formula: iem-boundary-merge]
#[inline]
pub fn boundary_merge(l_bi: f64, l_ic: f64) -> (f64, f64) {
    let l = l_bi + l_ic;
    (l_bi / l, l)
}
