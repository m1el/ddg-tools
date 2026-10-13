//! The closed-form formulas of the signpost data structure, as plain functions of lengths and
//! angles.
//!
//! Each carries a `[formula: it-…]` tag (see the repository's `CLAUDE.md`, "Formula
//! traceability"). Areas, angles and cotangents of single triangles come from
//! [`ddg_core::geometry`] (Kahan's forms, accurate for needle and cap triangles).
//!
//! Functions that take the lengths of a triangle return `None` when they do not form a triangle
//! of positive area in floating point arithmetic; nothing is clamped.

use std::f64::consts::{PI, TAU};

use ddg_core::geometry::{corner_angle, corner_cotan, triangle_area};

/// The factor that turns an angle measured in the triangles around a vertex into a normalized
/// angle: `2π / Θ` at an interior vertex and `π / Θ` at a boundary vertex, `Θ` being the sum of
/// the corner angles at the vertex (paper Eq. 1 and the text below it; SPEC N2 for the
/// boundary). With it the directions around an interior vertex fill `[0, 2π)` and those around
/// a boundary vertex `[0, π]`, whatever the cone angle.
///
/// [formula: it-angle-normalization]
#[inline]
pub fn angle_scale(theta: f64, boundary: bool) -> f64 {
    (if boundary { PI } else { TAU }) / theta
}

/// The signpost of the next halfedge counter-clockwise around a vertex: `φ_ik = φ_ij +
/// s · θ_i^{jk}` with `s` the [`angle_scale`] of the vertex and `θ_i^{jk}` the corner angle
/// between the two halfedges (paper Eq. 3, Alg. 1). At an interior vertex the result is reduced
/// to `[0, 2π)`; at a boundary vertex directions do not wrap.
///
/// [formula: it-signpost-update]
#[inline]
pub fn signpost_update(phi_ij: f64, corner: f64, theta: f64, boundary: bool) -> f64 {
    let phi = phi_ij + angle_scale(theta, boundary) * corner;
    if !boundary && phi >= TAU {
        phi - TAU
    } else {
        phi
    }
}

/// The two triangles `ijk` and `jil` on an edge `ij`, unfolded into the plane with `i` at the
/// origin and `j` at `(ℓ_ij, 0)`: `k = (x_k, h_k)` above the edge and `l = (x_l, −h_l)` below.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Diamond {
    /// Length of the shared edge `ij`.
    pub l_ij: f64,
    /// Abscissa of `k`, measured from `i` along `ij` (negative if the angle at `i` is obtuse).
    pub x_k: f64,
    /// Distance from `j` back to the foot of `k`: `ℓ_ij − x_k`, computed without that
    /// subtraction.
    pub x_k_from_j: f64,
    /// Height of `k` over `ij`, positive.
    pub h_k: f64,
    /// Abscissa of `l`, measured from `i`.
    pub x_l: f64,
    /// `ℓ_ij − x_l`, computed without that subtraction.
    pub x_l_from_j: f64,
    /// Height of `l` under `ij`, positive.
    pub h_l: f64,
}

/// Unfolds the triangles `ijk` and `jil` (see [`Diamond`]). With `a = ℓ_ij`:
/// `x_k = (a² + ℓ_ki² − ℓ_jk²) / 2a`, evaluated as `(a² + (ℓ_ki − ℓ_jk)(ℓ_ki + ℓ_jk)) / 2a`,
/// and `h_k = 2 A_ijk / a`; the same for `l` with `ℓ_il`, `ℓ_lj`.
///
/// [formula: it-unfold-diamond]
pub fn unfold_diamond(l_ij: f64, l_jk: f64, l_ki: f64, l_il: f64, l_lj: f64) -> Option<Diamond> {
    let area_k = triangle_area(l_ij, l_jk, l_ki)?;
    let area_l = triangle_area(l_ij, l_il, l_lj)?;
    let aa = l_ij * l_ij;
    let dk = (l_ki - l_jk) * (l_ki + l_jk);
    let dl = (l_il - l_lj) * (l_il + l_lj);
    let two_a = 2.0 * l_ij;
    Some(Diamond {
        l_ij,
        x_k: (aa + dk) / two_a,
        x_k_from_j: (aa - dk) / two_a,
        h_k: 2.0 * area_k / l_ij,
        x_l: (aa + dl) / two_a,
        x_l_from_j: (aa - dl) / two_a,
        h_l: 2.0 * area_l / l_ij,
    })
}

impl Diamond {
    /// Length of the other diagonal `kl`, the edge a flip creates:
    /// `ℓ_kl² = (x_k − x_l)² + (h_k + h_l)²`. The heights add without cancellation, and the
    /// difference of the abscissas is evaluated as
    /// `((ℓ_ki − ℓ_jk)(ℓ_ki + ℓ_jk) − (ℓ_il − ℓ_lj)(ℓ_il + ℓ_lj)) / 2 ℓ_ij`.
    /// It equals the paper's `√(ℓ_ik² + ℓ_il² − 2 ℓ_ik ℓ_il cos θ_i^{lk})` (Fig. 9).
    ///
    /// [formula: it-flip-length]
    pub fn flipped_length(&self) -> f64 {
        let dx = self.x_k - self.x_l;
        let dy = self.h_k + self.h_l;
        dx.hypot(dy)
    }

    /// Twice the signed areas of the triangles `ilk` and `jkl` that a flip creates:
    /// `x_k h_l + x_l h_k` and `(ℓ_ij − x_k) h_l + (ℓ_ij − x_l) h_k`. The quadrilateral
    /// `i l j k` is strictly convex, i.e. the diagonal `kl` crosses `ij` in its interior, exactly
    /// when both are positive.
    ///
    /// [formula: it-flip-convexity]
    pub fn flipped_areas(&self) -> [f64; 2] {
        [
            self.x_k * self.h_l + self.x_l * self.h_k,
            self.x_k_from_j * self.h_l + self.x_l_from_j * self.h_k,
        ]
    }
}

/// The smallest shape number a flip may give a triangle it creates (SPEC N5).
///
/// The *shape number* of a new triangle `T` is `2 A_T / (ℓ_new · ℓ_max(T))`: twice its area
/// over the product of the new edge and its longest side. Rounding the new length to a double
/// (a relative error `ε ≈ 1.1e-16`) moves the angles of `T`, as they follow from its three
/// lengths, by up to `ε` divided by the shape number: `∂γ/∂c = c / 2A` for the angle opposite
/// a side `c`, and `|∂α/∂c| ≤ a / 2A` for the others. With the limit `1e-6` a flip changes the
/// surface that the edge lengths describe by at most about `1e-10` radians per angle. A
/// triangle with a shape number below `√ε ≈ 1e-8` is not described by its edge lengths at
/// all: the slack of its triangle inequality is smaller than the rounding error of a length.
///
/// The shape number is small for a triangle that the flip makes nearly degenerate at an end
/// point of the old edge (a quadrilateral that is barely convex) and is not small for a needle
/// whose short side is the new edge. During Delaunay flipping it stayed above `4e-4` on every
/// mesh this crate was tested on.
pub const FLIP_SHAPE_LIMIT: f64 = 1e-6;

/// How far an edge `ij` with the triangles `ijk` and `jil` is from being Delaunay, as an
/// angle: `θ_k^{ij} + θ_l^{ji} − π`, the sum of the two angles opposite the edge minus `π`
/// (paper Eq. 4). The edge is Delaunay exactly when this is not positive.
///
/// [formula: it-delaunay-condition]
pub fn delaunay_angle_excess(l_ij: f64, l_jk: f64, l_ki: f64, l_il: f64, l_lj: f64) -> Option<f64> {
    let at_k = corner_angle(l_ij, l_jk, l_ki)?;
    let at_l = corner_angle(l_ij, l_il, l_lj)?;
    Some((at_k - PI) + at_l)
}

/// The sum of the cotangents of the two angles opposite an edge, `cot θ_k^{ij} + cot θ_l^{ji}
/// = sin(θ_k + θ_l) / (sin θ_k sin θ_l)`: twice the cotan weight of the edge. Since the sines
/// are positive it has the sign of `π − (θ_k + θ_l)`, so an edge is Delaunay exactly when its
/// cotan weight is not negative (paper Sec. 4).
///
/// [formula: it-cotan-weight]
pub fn opposite_cotan_sum(l_ij: f64, l_jk: f64, l_ki: f64, l_il: f64, l_lj: f64) -> Option<f64> {
    Some(corner_cotan(l_ij, l_jk, l_ki)? + corner_cotan(l_ij, l_il, l_lj)?)
}

/// The position of the third vertex `c` of a triangle `abc` (counter-clockwise) whose vertices
/// `a` and `b` are already placed in the plane: with `d = b − a` and `d⊥ = (−d_y, d_x)`,
///
/// ```text
/// c = a + α d + β d⊥,   α = (ℓ_ab² + ℓ_ca² − ℓ_bc²) / (2 ℓ_ab²),   β = 2 A / ℓ_ab²,
/// ```
///
/// with `A` the area of the triangle. `α` is evaluated as
/// `(ℓ_ab² + (ℓ_ca − ℓ_bc)(ℓ_ca + ℓ_bc)) / (2 ℓ_ab²)`. This is the step that unfolds the next
/// triangle along a traced path (paper Sec. 3.2.2); the given lengths are used, not the
/// distance of the placed points, so errors of the layout are not fed back into the shape.
///
/// [formula: it-unfold-apex]
#[inline]
pub fn unfold_apex(a: [f64; 2], b: [f64; 2], l_ab: f64, l_bc: f64, l_ca: f64) -> Option<[f64; 2]> {
    let area = triangle_area(l_ab, l_bc, l_ca)?;
    let ll = l_ab * l_ab;
    let alpha = (ll + (l_ca - l_bc) * (l_ca + l_bc)) / (2.0 * ll);
    let beta = 2.0 * area / ll;
    let d = [b[0] - a[0], b[1] - a[1]];
    Some([
        a[0] + alpha * d[0] - beta * d[1],
        a[1] + alpha * d[1] + beta * d[0],
    ])
}

/// Where the ray along the positive x-axis meets the segment from `tail` to `head`, given that
/// the end points lie on opposite sides of the axis (`tail_y ≤ 0 ≤ head_y`, not both zero, or
/// the reverse): the parameter `t = tail_y / (tail_y − head_y) ∈ [0, 1]` along the segment and
/// the abscissa `x = tail_x + t (head_x − tail_x)`. In floating point `t` is in `[0, 1]`
/// without clamping, because the numerator and the denominator have the same sign and the
/// numerator is not larger in magnitude.
///
/// [formula: it-ray-edge-intersection]
#[inline]
pub fn axis_crossing(tail: [f64; 2], head: [f64; 2]) -> (f64, f64) {
    let t = tail[1] / (tail[1] - head[1]);
    (t, tail[0] + t * (head[0] - tail[0]))
}

/// Barycentric coordinates of the point `p` in the triangle `q₀ q₁ q₂`, as ratios of signed
/// areas: `b_k = [p, q_{k+1}, q_{k+2}] / ([p, q₁, q₂] + [p, q₂, q₀] + [p, q₀, q₁])`, where
/// `[p, u, v]` is the cross product `(u − p) × (v − p)`. They sum to one by construction; a
/// coordinate is negative when `p` lies outside the triangle.
///
/// [formula: it-barycentric-coordinates]
pub fn barycentric(q: &[[f64; 2]; 3], p: [f64; 2]) -> [f64; 3] {
    let d = [
        [q[0][0] - p[0], q[0][1] - p[1]],
        [q[1][0] - p[0], q[1][1] - p[1]],
        [q[2][0] - p[0], q[2][1] - p[1]],
    ];
    let w = [
        d[1][0] * d[2][1] - d[1][1] * d[2][0],
        d[2][0] * d[0][1] - d[2][1] * d[0][0],
        d[0][0] * d[1][1] - d[0][1] * d[1][0],
    ];
    let sum = w[0] + w[1] + w[2];
    [w[0] / sum, w[1] / sum, w[2] / sum]
}

/// A triangle with halfedge lengths `l` placed in the plane in the frame of a face: corner 0
/// at the origin, corner 1 at `(l[0], 0)`, corner 2 above the x-axis ([`unfold_apex`]).
/// `l[k]` is the length of the side from corner `k` to corner `k + 1`.
pub fn triangle_layout(l: [f64; 3]) -> Option<[[f64; 2]; 3]> {
    let v1 = [l[0], 0.0];
    let v2 = unfold_apex([0.0, 0.0], v1, l[0], l[1], l[2])?;
    Some([[0.0, 0.0], v1, v2])
}

/// The distances from the point with barycentric coordinates `bary` of a triangle to its three
/// corners (paper Alg. 6), evaluated in the plane: with the corners `v_k` of
/// [`triangle_layout`] and `p = Σ bary_k v_k`, `d_k = |p − v_k|`. They are the lengths of the
/// edges that a vertex inserted at the point gets (Alg. 7), and with a point on a side the
/// lengths an edge split creates (Alg. 8). The paper's closed form
/// `d² = −(ℓ₁₂² u₁u₂ + ℓ₂₃² u₂u₃ + ℓ₃₁² u₃u₁)` (`u` the difference of the barycentric
/// coordinates) is the same number; it is not used because its three terms cancel in a thin
/// triangle, where the coordinates in the plane add in quadrature.
///
/// [formula: it-point-distances]
pub fn corner_distances(l: [f64; 3], bary: [f64; 3]) -> Option<[f64; 3]> {
    let v = triangle_layout(l)?;
    let p = [
        bary[1] * v[1][0] + bary[2] * v[2][0],
        bary[1] * v[1][1] + bary[2] * v[2][1],
    ];
    Some([0, 1, 2].map(|k| (p[0] - v[k][0]).hypot(p[1] - v[k][1])))
}

/// The circumcentre of a triangle in barycentric coordinates: with `a_k` the side opposite
/// corner `k`, the weight of corner `k` is `a_k² (a_{k+1}² + a_{k+2}² − a_k²)`, divided by the
/// sum of the three (which is `16 A²`). A coordinate is negative when the angle at that corner
/// is obtuse: the circumcentre then lies beyond the opposite side. `l[k]` is the length of the
/// side from corner `k` to corner `k + 1`, so `a_k = l[(k + 1) % 3]`.
///
/// [formula: it-circumcenter]
pub fn circumcenter_barycentric(l: [f64; 3]) -> Option<[f64; 3]> {
    triangle_area(l[0], l[1], l[2])?;
    let weight = |k: usize| {
        let (a, b, c) = (l[(k + 1) % 3], l[(k + 2) % 3], l[k]);
        let (hi, lo) = if b >= c { (b, c) } else { (c, b) };
        a * a * (lo * lo + (hi - a) * (hi + a))
    };
    let w = [weight(0), weight(1), weight(2)];
    let sum = w[0] + w[1] + w[2];
    (sum > 0.0 && sum.is_finite()).then(|| [w[0] / sum, w[1] / sum, w[2] / sum])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diamond_of_a_unit_square() {
        // i = (0,0), j = (1,1) is the diagonal of the unit square, k = (0,1), l = (1,0).
        let s = 2f64.sqrt();
        let d = unfold_diamond(s, 1.0, 1.0, 1.0, 1.0).unwrap();
        assert!((d.flipped_length() - s).abs() < 1e-15);
        let [a, b] = d.flipped_areas();
        assert!((a - 1.0).abs() < 1e-15 && (b - 1.0).abs() < 1e-15);
        assert!(delaunay_angle_excess(s, 1.0, 1.0, 1.0, 1.0).unwrap().abs() < 1e-15);
        assert!(opposite_cotan_sum(s, 1.0, 1.0, 1.0, 1.0).unwrap().abs() < 1e-15);
    }

    #[test]
    fn diamond_against_coordinates() {
        // i = (0,0), j = (4,0), k = (1,2), l = (3.5,-1).
        let dist = |p: [f64; 2], q: [f64; 2]| (p[0] - q[0]).hypot(p[1] - q[1]);
        let (i, j, k, l) = ([0.0, 0.0], [4.0, 0.0], [1.0, 2.0], [3.5, -1.0]);
        let d = unfold_diamond(dist(i, j), dist(j, k), dist(k, i), dist(i, l), dist(l, j)).unwrap();
        assert!((d.x_k - 1.0).abs() < 1e-14 && (d.h_k - 2.0).abs() < 1e-14);
        assert!((d.x_l - 3.5).abs() < 1e-14 && (d.h_l - 1.0).abs() < 1e-14);
        assert!((d.x_k + d.x_k_from_j - 4.0).abs() < 1e-14);
        assert!((d.flipped_length() - dist(k, l)).abs() < 1e-14);
    }

    #[test]
    fn nonconvex_diamond_has_a_negative_area() {
        // k = (-1, 1) and l = (-2, -1): the diagonal kl passes left of i = (0, 0).
        let dist = |p: [f64; 2], q: [f64; 2]| (p[0] - q[0]).hypot(p[1] - q[1]);
        let (i, j, k, l) = ([0.0, 0.0], [4.0, 0.0], [-1.0, 1.0], [-2.0, -1.0]);
        let d = unfold_diamond(dist(i, j), dist(j, k), dist(k, i), dist(i, l), dist(l, j)).unwrap();
        let [at_i, at_j] = d.flipped_areas();
        assert!(at_i < 0.0 && at_j > 0.0);
    }

    #[test]
    fn apex_and_crossing() {
        let c = unfold_apex([1.0, 1.0], [4.0, 1.0], 3.0, 5.0, 4.0).unwrap();
        assert!((c[0] - 1.0).abs() < 1e-15 && (c[1] - 5.0).abs() < 1e-15);
        let (t, x) = axis_crossing([1.0, -1.0], [3.0, 3.0]);
        assert!((t - 0.25).abs() < 1e-16 && (x - 1.5).abs() < 1e-16);
        let b = barycentric(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], [0.25, 0.5]);
        assert!((b[0] - 0.25).abs() < 1e-16 && (b[1] - 0.25).abs() < 1e-16);
        assert!((b[2] - 0.5).abs() < 1e-16);
    }

    #[test]
    fn distances_and_circumcentre_of_a_right_triangle() {
        // Corners (0,0), (4,0), (0,3): sides 4, 5, 3; the circumcentre is the midpoint of the
        // hypotenuse, from corner 1 to corner 2.
        let l = [4.0, 5.0, 3.0];
        let c = circumcenter_barycentric(l).unwrap();
        assert!(c[0].abs() < 1e-16 && (c[1] - 0.5).abs() < 1e-16 && (c[2] - 0.5).abs() < 1e-16);
        let d = corner_distances(l, c).unwrap();
        for x in d {
            assert!((x - 2.5).abs() < 1e-15);
        }
        // An obtuse triangle: the coordinate of the obtuse corner is negative.
        let obtuse = circumcenter_barycentric([1.0, 1.0, 1.9]).unwrap();
        assert!(obtuse[1] < 0.0 && obtuse[0] > 0.0 && obtuse[2] > 0.0);
        let r = corner_distances([1.0, 1.0, 1.9], obtuse).unwrap();
        assert!((r[0] - r[1]).abs() < 1e-14 && (r[1] - r[2]).abs() < 1e-14);
        // The point at a third of side 0 is 1/3 and 2/3 of it from its ends.
        let on_side = corner_distances(l, [2.0 / 3.0, 1.0 / 3.0, 0.0]).unwrap();
        assert!((on_side[0] - 4.0 / 3.0).abs() < 1e-15 && (on_side[1] - 8.0 / 3.0).abs() < 1e-15);
        assert!((on_side[2] - (16.0f64 / 9.0 + 9.0).sqrt()).abs() < 1e-15);
    }

    #[test]
    fn signposts_wrap_at_interior_vertices_only() {
        let theta = 5.0; // cone angle
        let phi = signpost_update(6.0, 1.0, theta, false);
        assert!((phi - (6.0 + TAU / 5.0 - TAU)).abs() < 1e-15);
        assert_eq!(signpost_update(1.0, 1.0, 2.0, true), 1.0 + PI / 2.0);
        assert_eq!(angle_scale(TAU, false), 1.0);
        assert_eq!(angle_scale(PI, true), 1.0);
    }
}
