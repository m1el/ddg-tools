//! Per-face quantities, their derivatives and the adjoint (back-propagation) form (SPEC §1–§2).
//!
//! A face `σ = (i, j, k)` (corners in the stored CCW order) has the unnormalized normal
//! `n = (f_j − f_i) × (f_k − f_i)`, area `A = |n|/2`, unit normal `N = n/|n|` and corner angles
//! `θ = atan2(|a × b|, a · b)`. The edge opposite corner `p` is `e_p` (CCW):
//! `e_i = f_k − f_j`, `e_j = f_i − f_k`, `e_k = f_j − f_i`.
//!
//! These functions are building blocks of [`crate::energy`]; they are public so that they can be
//! tested (finite differences) and reused, and take plain arrays.

/// 3-vector helpers (small, inlined).
pub(crate) mod v3 {
    #[inline(always)]
    pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }
    #[inline(always)]
    pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }
    #[inline(always)]
    pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
        [a[0] * s, a[1] * s, a[2] * s]
    }
    /// `a + s b`
    #[inline(always)]
    pub fn axpy(a: [f64; 3], s: f64, b: [f64; 3]) -> [f64; 3] {
        [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]]
    }
    #[inline(always)]
    pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
    #[inline(always)]
    pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }
    #[inline(always)]
    pub fn norm(a: [f64; 3]) -> f64 {
        dot(a, a).sqrt()
    }
    #[inline(always)]
    pub fn norm2(a: [f64; 3]) -> f64 {
        dot(a, a)
    }
    #[inline(always)]
    pub fn is_finite(a: [f64; 3]) -> bool {
        a[0].is_finite() && a[1].is_finite() && a[2].is_finite()
    }
}

use v3::*;

/// Normal, area and unit normal of one face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceFrame {
    /// `n = (f_j − f_i) × (f_k − f_i)` (twice the vector area).
    pub n: [f64; 3],
    /// `|n| = 2 A`.
    pub n_len: f64,
    /// Unit normal `N = n / |n|`.
    pub normal: [f64; 3],
}

impl FaceFrame {
    /// Area `A = |n| / 2`.
    pub fn area(&self) -> f64 {
        0.5 * self.n_len
    }
}

/// Normal `n`, `|n| = 2A` and unit normal `N` of the face with corners `p = [f_i, f_j, f_k]`
/// (CCW). `None` if `n` is zero or not finite (degenerate face: the unit normal is undefined).
///
/// [formula: dev-face-frame]
pub fn face_frame(p: [[f64; 3]; 3]) -> Option<FaceFrame> {
    let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
    let n_len = norm(n);
    if !n_len.is_finite() || n_len <= 0.0 {
        return None;
    }
    let normal = scale(n, 1.0 / n_len);
    if !is_finite(normal) {
        return None;
    }
    Some(FaceFrame { n, n_len, normal })
}

/// Angle between `a` and `b`: `atan2(|a × b|, a · b)` (accurate for small and near-π angles).
///
/// [formula: dev-corner-angle]
pub fn corner_angle(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(cross(a, b)).atan2(dot(a, b))
}

/// The three corner angles of the face `p = [f_i, f_j, f_k]`, angle `m` at corner `p[m]`.
pub fn corner_angles(p: [[f64; 3]; 3]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for m in 0..3 {
        let a = sub(p[(m + 1) % 3], p[m]);
        let b = sub(p[(m + 2) % 3], p[m]);
        out[m] = corner_angle(a, b);
    }
    out
}

/// Edge vector opposite corner `m` of `p`, CCW: `e_m = p[m+2] − p[m+1]` (indices mod 3).
#[inline]
pub fn opposite_edge(p: [[f64; 3]; 3], m: usize) -> [f64; 3] {
    sub(p[(m + 2) % 3], p[(m + 1) % 3])
}

/// Gradient of the area `A` with respect to corner `m` of `p`: `½ N × e_m` (Eq. 7).
///
/// [formula: dev-area-gradient]
pub fn area_gradient(p: [[f64; 3]; 3], normal: [f64; 3], m: usize) -> [f64; 3] {
    scale(cross(normal, opposite_edge(p, m)), 0.5)
}

/// Jacobian `∂N/∂f_m = (e_m × N) Nᵀ / (2A)` of the unit normal with respect to corner `m` (Eq. 8
/// with the corrected factor `1/(2A)`, D2). Row `r` is `∂N_r/∂f_m`.
///
/// [formula: dev-normal-gradient]
pub fn normal_jacobian(p: [[f64; 3]; 3], frame: &FaceFrame, m: usize) -> [[f64; 3]; 3] {
    let t = scale(cross(opposite_edge(p, m), frame.normal), 1.0 / frame.n_len);
    let nn = frame.normal;
    [scale(nn, t[0]), scale(nn, t[1]), scale(nn, t[2])]
}

/// Gradients of the angle `θ_m` at corner `m` of `p` with respect to `[f_m, f_{m+1}, f_{m+2}]`
/// (Eq. 9 with the corrected `|e|²`, D3):
/// `∇_{f_j} θ_i = N × (f_i − f_j)/|f_i − f_j|²`, `∇_{f_k} θ_i = N × (f_k − f_i)/|f_k − f_i|²`,
/// `∇_{f_i} θ_i = −(both)`, with `i = m`, `j = m+1`, `k = m+2` (mod 3).
///
/// [formula: dev-angle-gradient]
pub fn angle_gradients(p: [[f64; 3]; 3], normal: [f64; 3], m: usize) -> [[f64; 3]; 3] {
    let (fi, fj, fk) = (p[m], p[(m + 1) % 3], p[(m + 2) % 3]);
    let dij = sub(fi, fj);
    let dki = sub(fk, fi);
    let gj = scale(cross(normal, dij), 1.0 / norm2(dij));
    let gk = scale(cross(normal, dki), 1.0 / norm2(dki));
    [scale(add(gj, gk), -1.0), gj, gk]
}

/// Vector area `S_v = ½ Σ n_σ` of a star, summed in the given (star) order.
/// The unit vertex normal is `S_v / |S_v|` (undefined if `S_v = 0`); for an interior vertex
/// `S_v` does not depend on the position of `v`.
///
/// [formula: dev-vertex-normal]
pub fn vector_area<'a>(face_normals: impl IntoIterator<Item = &'a [f64; 3]>) -> [f64; 3] {
    let mut s = [0.0; 3];
    for n in face_normals {
        s = add(s, *n);
    }
    scale(s, 0.5)
}

/// Adjoint form of SPEC §2 for one face: given the sensitivity `y_σ = Σ_c ∂E/∂N_σ` (summed over
/// the three corners, i.e. three stars), `z_sum = Σ_{u∈σ} ∂E/∂S_u` (vertex-normal sensitivities
/// of the three vertices) and the angle sensitivities `s[m] = ∂E/∂θ_m`, return the gradient
/// contributions `q[m] = ∂E/∂f_m` of the three corners:
///
/// ```text
/// w = (y − (y·N) N) / |n|  +  ½ z_sum              (sensitivity w.r.t. n)
/// q_m = w × e_m                                     (normal term)
/// then for corners m = 0, 1, 2 the angle terms of Eq. 9 times s[m]
/// ```
///
/// [formula: dev-backprop-face]
pub fn backprop_face(
    p: [[f64; 3]; 3],
    frame: &FaceFrame,
    y: [f64; 3],
    z_sum: [f64; 3],
    s: [f64; 3],
) -> [[f64; 3]; 3] {
    let nn = frame.normal;
    let tang = axpy(y, -dot(y, nn), nn);
    let w = axpy(scale(tang, 1.0 / frame.n_len), 0.5, z_sum);
    let mut q = [
        cross(w, opposite_edge(p, 0)),
        cross(w, opposite_edge(p, 1)),
        cross(w, opposite_edge(p, 2)),
    ];
    for (m, &sm) in s.iter().enumerate() {
        if sm == 0.0 {
            continue;
        }
        let (i, j, k) = (m, (m + 1) % 3, (m + 2) % 3);
        let dij = sub(p[i], p[j]);
        let dki = sub(p[k], p[i]);
        let gj = scale(cross(nn, dij), sm / norm2(dij));
        let gk = scale(cross(nn, dki), sm / norm2(dki));
        q[j] = add(q[j], gj);
        q[k] = add(q[k], gk);
        q[i] = sub(q[i], add(gj, gk));
    }
    q
}
