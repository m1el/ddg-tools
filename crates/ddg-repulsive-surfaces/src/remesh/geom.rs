//! Geometric predicates of the remesher: circumcenters, the Delaunay test, the foldover
//! (normal-turn) test (`SPEC.md` §10).

use crate::geometry::v3::{add, cross, dot, norm, scale, sub};

/// The foldover threshold `π − 2` rad (`SPEC.md` §10; the reference's `checkFoldover(…, 2)`).
pub const FOLDOVER_LIMIT: f64 = std::f64::consts::PI - 2.0;

/// The angle between two vectors in `[0, π]`, `atan2(|u × v|, u · v)` (accurate near 0 and π).
#[inline]
pub fn angle_between(u: [f64; 3], v: [f64; 3]) -> f64 {
    norm(cross(u, v)).atan2(dot(u, v))
}

/// The (unnormalized) normal `(p₁ − p₀) × (p₂ − p₀)` of a triangle; its length is twice the area.
#[inline]
pub fn tri_normal(p: [[f64; 3]; 3]) -> [f64; 3] {
    cross(sub(p[1], p[0]), sub(p[2], p[0]))
}

/// Does a face whose (unnormalized) normal was `before` and is `after` fold over: the new
/// triangle is degenerate (zero or non-finite area) or its normal turned by more than `π − 2`?
/// Both normals must come from the same vertex order (the test is then orientation-free).
/// Evaluated as `before · after < cos(π − 2) |before| |after|` (no `atan2`).
#[inline]
pub fn folds(before: [f64; 3], after: [f64; 3]) -> bool {
    let a = norm(after);
    !(a > 0.0 && a.is_finite()) || dot(before, after) < COS_FOLDOVER_LIMIT * norm(before) * a
}

/// `cos(π − 2) = −cos 2`.
pub const COS_FOLDOVER_LIMIT: f64 = 0.4161468365471424;

/// The corner angle at `v` between `w − v` and `x − v`.
#[inline]
pub fn corner_angle(v: [f64; 3], w: [f64; 3], x: [f64; 3]) -> f64 {
    angle_between(sub(w, v), sub(x, v))
}

/// The Delaunay test of an interior edge `p q` with opposite vertices `c` and `d`: the edge
/// must be flipped iff the two opposite angles sum to more than `π`, equivalently iff its
/// cotan weight `½(cot ∠c + cot ∠d)` is negative (`cot α + cot β = sin(α + β)/(sin α sin β)`
/// and both sines are positive). Evaluated in `f64` (the reference uses `float`, D16) in the
/// cotan form without divisions: with `cot ∠c = d_c/s_c` (`d_c` the dot product and `s_c > 0`
/// the cross-product length of the two edge vectors at `c`), the test is
/// `d_c s_d + d_d s_c < 0`. A degenerate corner (`s = 0`) never asks for a flip.
///
/// [formula: rs-delaunay-test]
#[inline]
pub fn violates_delaunay(p: [f64; 3], q: [f64; 3], c: [f64; 3], d: [f64; 3]) -> bool {
    let (uc, vc) = (sub(p, c), sub(q, c));
    let (ud, vd) = (sub(p, d), sub(q, d));
    let (dc, sc) = (dot(uc, vc), norm(cross(uc, vc)));
    let (dd, sd) = (dot(ud, vd), norm(cross(ud, vd)));
    sc > 0.0 && sd > 0.0 && dc * sd + dd * sc < 0.0
}

/// The circumcenter of a triangle from its barycentric weights
/// `a²(b² + c² − a²) : b²(c² + a² − b²) : c²(a² + b² − c²)` (`a, b, c` the lengths of the sides
/// opposite `p₀, p₁, p₂`), normalized to sum 1; the sum is `16 · area²` and positive for a
/// non-degenerate triangle (`ref: remeshing.cpp:73–87`).
///
/// [formula: rs-circumcenter]
pub fn circumcenter(p: [[f64; 3]; 3]) -> [f64; 3] {
    let a2 = norm2(sub(p[2], p[1]));
    let b2 = norm2(sub(p[2], p[0]));
    let c2 = norm2(sub(p[1], p[0]));
    let w = [
        a2 * (b2 + c2 - a2),
        b2 * (c2 + a2 - b2),
        c2 * (a2 + b2 - c2),
    ];
    let s = w[0] + w[1] + w[2];
    add(
        add(scale(p[0], w[0] / s), scale(p[1], w[1] / s)),
        scale(p[2], w[2] / s),
    )
}

/// The barycenter `(p₀ + p₁ + p₂)/3`.
#[inline]
pub fn barycenter(p: [[f64; 3]; 3]) -> [f64; 3] {
    scale(add(add(p[0], p[1]), p[2]), 1.0 / 3.0)
}

#[inline]
fn norm2(v: [f64; 3]) -> f64 {
    dot(v, v)
}

/// The tangential smoothing displacement of a vertex at `x` (`SPEC.md` §10 step 2)
///
/// ```text
/// u = ρ Σ_S a_S (c_S − x) / Σ_S a_S
/// ```
///
/// over the faces `S` of its star with areas `a_S` and centers `c_S` (circumcenter, or
/// barycenter for faces with a boundary edge). Since `Σ_S a_S = 3 a_v`, this is the
/// reference's `½ Σ a_S (c_S − x) / (3 a_v)` for `ρ = ½` (`ref: remeshing.cpp:569–578`). Not yet
/// projected onto the tangent plane. Returns `None` if the total area is not positive.
///
/// [formula: rs-smoothing-weights]
pub fn smoothing_displacement(
    x: [f64; 3],
    rho: f64,
    faces: impl Iterator<Item = (f64, [f64; 3])>,
) -> Option<[f64; 3]> {
    let mut sum = [0.0; 3];
    let mut area = 0.0;
    for (a, c) in faces {
        sum = add(sum, scale(sub(c, x), a));
        area += a;
    }
    if area > 0.0 && area.is_finite() {
        Some(scale(sum, rho / area))
    } else {
        None
    }
}
