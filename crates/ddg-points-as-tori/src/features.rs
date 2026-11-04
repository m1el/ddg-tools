//! Network input features (SPEC.md §4.1; NOTES-REFERENCE.md §2.2), ported from the reference's
//! `InputFeatures.__call__` and `orthonormal_basis_with_rotation` (`network.py`, points-as-tori,
//! MIT, Copyright (c) 2026 Nicole Feng).

/// `max(‖n‖, NAN_EPSILON)` guard of the reference.
pub const NORMAL_EPS: f64 = 1e-10;

/// How the features are computed (decision D3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FeatureOrder {
    /// Ours (default): `d = (p_j − p_i)/σ_i` difference-first in f64, frame and dot products in
    /// f64, then rounded to f32. Translation invariant to round-off.
    #[default]
    DifferenceFirst,
    /// The reference's arithmetic: f32 throughout, `p_j/σ − p_i/σ` on absolute coordinates.
    /// Used for parity with the authors' program.
    Reference,
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross32(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot32(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The local frame `(ŝ, t̂)` of a unit normal `c`: `u = e_y` if `|c_x| > 0.9` else `e_x`,
/// `ŝ = (c × u)/‖c × u‖`, `t̂ = c × ŝ`.
///
/// `c` is used as given (callers pass a unit vector). `‖c × u‖ ≥ √(1 − 0.9²) > 0.43` for unit `c`,
/// so no guard is needed.
///
/// [formula: pat-local-frame]
pub fn local_frame(c: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let u = if c[0].abs() > 0.9 { [0.0, 1.0, 0.0] } else { [1.0, 0.0, 0.0] };
    let s = cross(c, u);
    let ls = dot(s, s).sqrt();
    let s = [s[0] / ls, s[1] / ls, s[2] / ls];
    (s, cross(c, s))
}

/// [`local_frame`] in f32 with the reference's `max(‖s‖, 1e-10)` guard (the Rodrigues rotation
/// by 0 is the identity, so it is omitted).
///
/// [formula: pat-local-frame]
pub fn local_frame_f32(c: [f32; 3]) -> ([f32; 3], [f32; 3]) {
    let u = if c[0].abs() > 0.9 { [0.0, 1.0, 0.0] } else { [1.0, 0.0, 0.0] };
    let s = cross32(c, u);
    let ls = dot32(s, s).sqrt().max(NORMAL_EPS as f32);
    let s = [s[0] / ls, s[1] / ls, s[2] / ls];
    (s, cross32(c, s))
}

/// Median of `v` (sorted in place): the middle value, or the mean of the two middle values for
/// an even count. `σ(c·P) = c·σ(P)` for `c > 0`, translation invariant.
///
/// [formula: pat-median-scale]
pub fn median(v: &mut [f64]) -> f64 {
    v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let m = v.len();
    if m % 2 == 1 { v[m / 2] } else { 0.5 * (v[m / 2 - 1] + v[m / 2]) }
}

/// [`median`] in f32 (`jnp.median`).
///
/// [formula: pat-median-scale]
pub fn median_f32(v: &mut [f32]) -> f32 {
    v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let m = v.len();
    if m % 2 == 1 { v[m / 2] } else { 0.5 * (v[m / 2 - 1] + v[m / 2]) }
}

/// Scratch buffers for [`neighborhood_features`].
#[derive(Debug, Default, Clone)]
pub struct FeatureScratch {
    d64: Vec<f64>,
    d32: Vec<f32>,
    p32: Vec<[f32; 3]>,
}

/// The features of one neighborhood `row = (j_0 = i, j_1, …)`, written to `out` (`6·row.len()`
/// values, token-major: `(d·ŝ, d·t̂, d·c, n̂·ŝ, n̂·t̂, n̂·c)` per token), and its scale `σ_i`.
///
/// `σ_i = median_m ‖p_{j_m} − p_i‖` over all tokens including `m = 0`; `d = (p_j − p_i)/σ_i`;
/// `n̂_j = n_j / max(‖n_j‖, 1e-10)`; the frame is that of `c = n̂_i`. In
/// [`FeatureOrder::Reference`] the returned σ is the f32 value promoted.
///
/// σ = 0 (all tokens at one position) gives non-finite features; the network then reports the
/// row as non-finite.
///
/// [formula: pat-features]
/// [formula: pat-feature-invariance]
pub fn neighborhood_features(
    positions: &[[f64; 3]],
    normals: &[[f64; 3]],
    row: &[u32],
    order: FeatureOrder,
    scratch: &mut FeatureScratch,
    out: &mut [f32],
) -> f64 {
    let k = row.len();
    debug_assert_eq!(out.len(), 6 * k);
    match order {
        FeatureOrder::DifferenceFirst => {
            let pi = positions[row[0] as usize];
            let c = unit(normals[row[0] as usize]);
            let (s, t) = local_frame(c);
            scratch.d64.clear();
            scratch.d64.extend(row.iter().map(|&j| {
                let p = positions[j as usize];
                let d = [p[0] - pi[0], p[1] - pi[1], p[2] - pi[2]];
                dot(d, d).sqrt()
            }));
            let sigma = median(&mut scratch.d64);
            for (m, &j) in row.iter().enumerate() {
                let p = positions[j as usize];
                let d = [(p[0] - pi[0]) / sigma, (p[1] - pi[1]) / sigma, (p[2] - pi[2]) / sigma];
                let nj = unit(normals[j as usize]);
                let o = &mut out[6 * m..6 * m + 6];
                o[0] = dot(d, s) as f32;
                o[1] = dot(d, t) as f32;
                o[2] = dot(d, c) as f32;
                o[3] = dot(nj, s) as f32;
                o[4] = dot(nj, t) as f32;
                o[5] = dot(nj, c) as f32;
            }
            sigma
        }
        FeatureOrder::Reference => {
            let f = |v: [f64; 3]| [v[0] as f32, v[1] as f32, v[2] as f32];
            scratch.p32.clear();
            scratch.p32.extend(row.iter().map(|&j| f(positions[j as usize])));
            let p0 = scratch.p32[0];
            scratch.d32.clear();
            scratch.d32.extend(scratch.p32.iter().map(|p| {
                let d = [p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]];
                dot32(d, d).sqrt()
            }));
            let sigma = median_f32(&mut scratch.d32);
            let c = unit32(f(normals[row[0] as usize]));
            let (s, t) = local_frame_f32(c);
            let q0 = [p0[0] / sigma, p0[1] / sigma, p0[2] / sigma];
            for (m, &j) in row.iter().enumerate() {
                let p = scratch.p32[m];
                let d = [p[0] / sigma - q0[0], p[1] / sigma - q0[1], p[2] / sigma - q0[2]];
                let nj = unit32(f(normals[j as usize]));
                let o = &mut out[6 * m..6 * m + 6];
                o[0] = dot32(d, s);
                o[1] = dot32(d, t);
                o[2] = dot32(d, c);
                o[3] = dot32(nj, s);
                o[4] = dot32(nj, t);
                o[5] = dot32(nj, c);
            }
            sigma as f64
        }
    }
}

#[inline]
fn unit(n: [f64; 3]) -> [f64; 3] {
    let l = dot(n, n).sqrt().max(NORMAL_EPS);
    [n[0] / l, n[1] / l, n[2] / l]
}

#[inline]
fn unit32(n: [f32; 3]) -> [f32; 3] {
    let l = dot32(n, n).sqrt().max(NORMAL_EPS as f32);
    [n[0] / l, n[1] / l, n[2] / l]
}
