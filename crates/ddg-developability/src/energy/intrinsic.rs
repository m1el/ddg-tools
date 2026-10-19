//! Intrinsic energies in the log map of the Gauss sphere at the vertex normal (App. B.5, B.5.1),
//! SPEC §3.6–§3.7.

use super::{Scratch, Star, StarOut};
use crate::eig::sym_eig2;
use crate::geometry::v3::*;

/// Faces whose normal is within this angle of the antipode of the vertex normal make the log map
/// singular; the intrinsic energies are then undefined (D12).
pub const LOG_MAP_SINGULAR_ANGLE: f64 = 1e-6;

/// Reason reported when the vector area of a star vanishes.
pub const REASON_ZERO_VECTOR_AREA: &str =
    "the vector area of the star is zero (or not finite), so the vertex normal is undefined";
/// Reason reported when a face normal is nearly opposite the vertex normal.
pub const REASON_ANTIPODAL_NORMAL: &str = "a face normal is within 1e-6 rad of the antipode of the \
     vertex normal, where the log map of the Gauss sphere is singular";

/// `G(φ) = φ / sin φ` (series `1 + φ²/6 + 7φ⁴/360` below `1e-3`).
#[inline]
pub fn log_map_g(phi: f64) -> f64 {
    if phi < 1e-3 {
        let p2 = phi * phi;
        1.0 + p2 / 6.0 + 7.0 * p2 * p2 / 360.0
    } else {
        phi / phi.sin()
    }
}

/// `H(φ) = G′(φ) / sin φ = (sin φ − φ cos φ) / sin³ φ` (series `1/3 + 2φ²/15 + 2φ⁴/63` below
/// `1e-2`, where the closed form loses digits to cancellation).
#[inline]
pub fn log_map_h(phi: f64) -> f64 {
    if phi < 1e-2 {
        let p2 = phi * phi;
        1.0 / 3.0 + 2.0 * p2 / 15.0 + 2.0 * p2 * p2 / 63.0
    } else {
        let (sn, cs) = phi.sin_cos();
        (sn - phi * cs) / (sn * sn * sn)
    }
}

/// The log map of the unit sphere at `n0`, applied to the unit vector `normal`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LogMap {
    /// `Ñ = log_{n0}(N) = G(φ) (N − c n0)`, tangent at `n0`, `|Ñ| = φ`.
    pub tilde: [f64; 3],
    /// `c = ⟨n0, N⟩`.
    pub c: f64,
    /// `φ = atan2(|n0 × N|, c)`, the geodesic distance.
    pub phi: f64,
    /// `G(φ)`.
    pub g: f64,
    /// `H(φ)`.
    pub h: f64,
}

/// `log_{n0}(N) = (φ / sin φ)(N − ⟨n0,N⟩ n0)`, `φ = atan2(|n0 × N|, ⟨n0,N⟩)`.
///
/// [formula: dev-log-map]
pub fn log_map(n0: [f64; 3], normal: [f64; 3]) -> LogMap {
    let c = dot(n0, normal);
    let phi = norm(cross(n0, normal)).atan2(c);
    let g = log_map_g(phi);
    let h = log_map_h(phi);
    LogMap {
        tilde: scale(axpy(normal, -c, n0), g),
        c,
        phi,
        g,
        h,
    }
}

/// Adjoint of the log map: for a sensitivity `z = ∂E/∂Ñ`, returns `(∂E/∂N, ∂E/∂n0)`:
/// `K = G (z·n0) + H z·(N − c n0)`, `∂/∂N = G z − K n0`, `∂/∂n0 = −G c z − K N`
/// (the smooth, corrected form of Eq. 11, D6; no `1/sin φ` at `φ = 0`).
///
/// [formula: dev-log-map-gradient]
pub fn log_map_adjoint(
    n0: [f64; 3],
    normal: [f64; 3],
    lm: &LogMap,
    z: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let t = axpy(normal, -lm.c, n0);
    let kk = lm.g * dot(z, n0) + lm.h * dot(z, t);
    let d_n = axpy(scale(z, lm.g), -kk, n0);
    let d_n0 = axpy(scale(z, -lm.g * lm.c), -kk, normal);
    (d_n, d_n0)
}

/// Deterministic orthonormal basis `(t1, t2)` of the plane orthogonal to the unit vector `n`,
/// with `t1 × t2 = n`.
pub fn tangent_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let ax = [n[0].abs(), n[1].abs(), n[2].abs()];
    let e = if ax[0] <= ax[1] && ax[0] <= ax[2] {
        [1.0, 0.0, 0.0]
    } else if ax[1] <= ax[2] {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let t1 = cross(n, e);
    let t1 = scale(t1, 1.0 / norm(t1));
    let t2 = cross(n, t1);
    (t1, t2)
}

/// Vertex normal and log-mapped normals of a star; fails where the energy is undefined (D12).
fn prepare(star: &Star, scratch: &mut Scratch) -> Result<([f64; 3], f64), &'static str> {
    let sv = star.vector_area;
    let len = norm(sv);
    if !len.is_finite() || len <= 0.0 {
        return Err(REASON_ZERO_VECTOR_AREA);
    }
    let n0 = scale(sv, 1.0 / len);
    scratch.logs.clear();
    for n in star.normals {
        let lm = log_map(n0, *n);
        if lm.phi.is_nan() || std::f64::consts::PI - lm.phi < LOG_MAP_SINGULAR_ANGLE {
            return Err(REASON_ANTIPODAL_NORMAL);
        }
        scratch.logs.push(lm);
    }
    Ok((n0, len))
}

/// Back-propagate `z_σ = ∂e/∂Ñ_σ` (in `scratch.v1`) through the log map into `y_c` and the
/// vertex-normal sensitivity `z_v = (y_n0 − (y_n0·n0) n0)/|S_v|`.
fn finish(
    star: &Star,
    n0: [f64; 3],
    s_len: f64,
    y: &mut [[f64; 3]],
    scratch: &Scratch,
) -> [f64; 3] {
    let mut y_n0 = [0.0; 3];
    for ((n, lm), (z, yc)) in star
        .normals
        .iter()
        .zip(&scratch.logs)
        .zip(scratch.v1.iter().zip(y.iter_mut()))
    {
        let (d_n, d_n0) = log_map_adjoint(n0, *n, lm, *z);
        *yc = d_n;
        y_n0 = add(y_n0, d_n0);
    }
    scale(axpy(y_n0, -dot(y_n0, n0), n0), 1.0 / s_len)
}

/// Intrinsic width `e_v = λ_min(Ã_v)`, `Ã_v = Σ θ Ñ_σ Ñ_σᵀ` restricted to the tangent plane at the
/// vertex normal `n0 = S_v/|S_v|` (App. B.5), by [`sym_eig2`] in the basis [`tangent_basis`]. With
/// `x ⊥ n0` the unit eigenvector, `dλ = xᵀ dÃ x` (the basis terms vanish), so
/// `∂e/∂Ñ_σ = 2θ (x·Ñ_σ) x`, `s_c = (x·Ñ_σ)²`, then the log-map adjoint.
///
/// Margin: `(λ_2 − λ_1) / Σθ` (`Σθ = tr A_v` of §3.4: an absolute gap on the scale of `E_λ`;
/// `tr Ã_v = Σ θ φ²` would make the margin of a nearly flat star scale-free rounding noise).
///
/// [formula: dev-tangent-eigen-derivative]
pub fn intrinsic(
    star: &Star,
    y: &mut [[f64; 3]],
    s: &mut [f64],
    scratch: &mut Scratch,
) -> Result<StarOut, &'static str> {
    let (n0, s_len) = prepare(star, scratch)?;
    let (t1, t2) = tangent_basis(n0);
    let mut m = [[0.0; 2]; 2];
    for (lm, &th) in scratch.logs.iter().zip(star.angles) {
        let a = [dot(lm.tilde, t1), dot(lm.tilde, t2)];
        m[0][0] += th * a[0] * a[0];
        m[0][1] += th * a[0] * a[1];
        m[1][1] += th * a[1] * a[1];
    }
    let eig = sym_eig2(m).map_err(|_| "the intrinsic covariance matrix has a non-finite entry")?;
    let yv = eig.vectors[0];
    let x = axpy(scale(t1, yv[0]), yv[1], t2);
    scratch.v1.clear();
    for ((lm, &th), sc) in scratch.logs.iter().zip(star.angles).zip(s.iter_mut()) {
        let xn = dot(x, lm.tilde);
        scratch.v1.push(scale(x, 2.0 * th * xn));
        *sc = xn * xn;
    }
    let z = finish(star, n0, s_len, y, scratch);
    let theta_sum: f64 = star.angles.iter().sum();
    let margin = (eig.values[1] - eig.values[0]) / theta_sum;
    Ok(StarOut {
        value: eig.values[0],
        margin,
        z,
    })
}

/// The width direction of a hull-edge candidate: `u = n0 × ŵ`, `ŵ = w/|w|`. Its derivative in
/// the plane is `du = −(u·dw) ŵ/|w|`.
///
/// [formula: dev-width-altitude-gradient]
pub fn width_direction(n0: [f64; 3], w: [f64; 3]) -> Option<([f64; 3], [f64; 3], f64)> {
    let len = norm(w);
    if !len.is_finite() || len <= 0.0 {
        return None;
    }
    let wh = scale(w, 1.0 / len);
    Some((cross(n0, wh), wh, len))
}

/// An intrinsic-max candidate: the pair `(a, b)`, `u`, `ŵ`, `|w|`.
type WidthCandidate = (usize, usize, [f64; 3], [f64; 3], f64);

/// Branching-free intrinsic width (App. B.5.1): `e_v = (min width of conv{Ñ_σ})²`,
/// `ψ(u) = (max_σ⟨Ñ_σ,u⟩ − min_σ⟨Ñ_σ,u⟩)²` minimized over the hull-edge normals
/// `u_ab = n0 × (Ñ_b − Ñ_a)/|Ñ_b − Ñ_a|`, `a < b` (rotating calipers; the paper's `±` enumeration
/// is a superset with the same minimum, D7). Chosen: lowest `(ψ, a, b)`; then `s = argmax`,
/// `t = argmin` of `⟨Ñ,u⟩` (lowest index on ties). No candidate gives `e_v = 0`.
///
/// Gradient (`d = Ñ_s − Ñ_t`, `Δ = d·u`): `z_s += 2Δu`, `z_t −= 2Δu`,
/// `z_a += 2Δ (d·ŵ)/|w| u`, `z_b −= 2Δ (d·ŵ)/|w| u`, `s_c = 0`, then the log-map adjoint.
/// Margin: second-best minus best `ψ`.
///
/// [formula: dev-intrinsic-width]
/// [formula: dev-min-width-hull-edge]
/// [formula: dev-min-width-gradient]
pub fn intrinsic_max(
    star: &Star,
    y: &mut [[f64; 3]],
    s: &mut [f64],
    scratch: &mut Scratch,
) -> Result<StarOut, &'static str> {
    let (n0, s_len) = prepare(star, scratch)?;
    s.fill(0.0);
    let k = star.normals.len();
    let logs = &scratch.logs;
    let width = |u: [f64; 3]| {
        let (mut hi, mut lo) = (f64::NEG_INFINITY, f64::INFINITY);
        for lm in logs {
            let p = dot(lm.tilde, u);
            hi = hi.max(p);
            lo = lo.min(p);
        }
        (hi - lo) * (hi - lo)
    };
    let mut best: Option<WidthCandidate> = None;
    let mut best_psi = f64::INFINITY;
    let mut second = f64::INFINITY;
    for a in 0..k {
        for b in a + 1..k {
            let Some((u, wh, wl)) = width_direction(n0, sub(logs[b].tilde, logs[a].tilde)) else {
                continue;
            };
            let psi = width(u);
            if psi < best_psi {
                second = best_psi;
                best_psi = psi;
                best = Some((a, b, u, wh, wl));
            } else if psi < second {
                second = psi;
            }
        }
    }
    scratch.v1.clear();
    scratch.v1.resize(k, [0.0; 3]);
    let Some((a, b, u, wh, wl)) = best else {
        y.fill([0.0; 3]);
        return Ok(StarOut {
            value: 0.0,
            margin: f64::INFINITY,
            z: [0.0; 3],
        });
    };
    let (mut si, mut ti) = (0, 0);
    let (mut hi, mut lo) = (f64::NEG_INFINITY, f64::INFINITY);
    for (i, lm) in logs.iter().enumerate() {
        let p = dot(lm.tilde, u);
        if p > hi {
            hi = p;
            si = i;
        }
        if p < lo {
            lo = p;
            ti = i;
        }
    }
    let d = sub(logs[si].tilde, logs[ti].tilde);
    let delta = dot(d, u);
    let zs = scale(u, 2.0 * delta);
    let za = scale(u, 2.0 * delta * dot(d, wh) / wl);
    let z = &mut scratch.v1;
    z[si] = add(z[si], zs);
    z[ti] = sub(z[ti], zs);
    z[a] = add(z[a], za);
    z[b] = sub(z[b], za);
    let zv = finish(star, n0, s_len, y, scratch);
    Ok(StarOut {
        value: best_psi,
        margin: second - best_psi,
        z: zv,
    })
}
