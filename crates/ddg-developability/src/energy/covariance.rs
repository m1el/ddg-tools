//! Covariance energy `E_λ` (Eq. 3–4, 10) and max covariance energy `λ_max` (Eq. 6, App. B.4),
//! SPEC §3.4–§3.5.

use super::{Star, StarOut};
use crate::eig::sym_eig3;
use crate::geometry::v3::*;

/// Covariance energy `e_v = λ_min(A_v)`, `A_v = Σ_σ θ_v^σ N_σ N_σᵀ` (Eq. 3 = Eq. 4 by the
/// Rayleigh quotient), by the cyclic Jacobi solver [`sym_eig3`] (D29); the result is not clamped
/// at 0 (D14). With `x` the unit eigenvector of `λ_min`, the (super)gradient is
/// `y_c = 2 θ (x·N_σ) x`, `s_c = (x·N_σ)²` (Eq. 10 with the corrected Eq. 8–9): `∇_A λ = x xᵀ` at
/// a simple eigenvalue (Magnus 1985) and `λ_min(B) ≤ λ_min(A) + xᵀ(B − A)x` always.
///
/// Angle weights make `A_v` invariant under retessellation of a planar region: splitting a corner
/// into two corners with the same normal leaves `A_v` unchanged.
///
/// Margin: `(λ_2 − λ_1) / tr(A_v)`.
///
/// [formula: dev-covariance-matrix]
/// [formula: dev-covariance-gradient]
/// [formula: dev-simple-eigenvalue-derivative]
/// [formula: dev-eigenvalue-supergradient]
/// [formula: dev-angle-weight-retessellation]
pub fn covariance(star: &Star, y: &mut [[f64; 3]], s: &mut [f64]) -> Result<StarOut, &'static str> {
    let mut a = [[0.0; 3]; 3];
    for (n, &t) in star.normals.iter().zip(star.angles) {
        for i in 0..3 {
            let tn = t * n[i];
            for j in i..3 {
                a[i][j] += tn * n[j];
            }
        }
    }
    let eig = sym_eig3(a).map_err(|_| "the covariance matrix has a non-finite entry")?;
    let lam = eig.values[0];
    let x = eig.vectors[0];
    for ((n, &t), (yc, sc)) in star.normals.iter().zip(star.angles).zip(y.iter_mut().zip(s.iter_mut()))
    {
        let xn = dot(x, *n);
        *yc = scale(x, 2.0 * t * xn);
        *sc = xn * xn;
    }
    let tr = a[0][0] + a[1][1] + a[2][2];
    let margin = if tr > 0.0 { (eig.values[1] - eig.values[0]) / tr } else { 0.0 };
    Ok(StarOut { value: lam, margin, z: [0.0; 3] })
}

/// Spherical circumcenter of the sites `a, b, c`: `n = (b − a) × (c − a)`, `v = n/|n|`, and
/// `g = ⟨v,a⟩ = ⟨v,b⟩ = ⟨v,c⟩ = det(a,b,c)/|n|`. `None` if `n` is exactly zero or not finite.
///
/// [formula: dev-spherical-circumcenter]
pub fn spherical_circumcenter(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Option<([f64; 3], f64)> {
    let n = cross(sub(b, a), sub(c, a));
    let len = norm(n);
    if !(len > 0.0) || !len.is_finite() {
        return None;
    }
    let v = scale(n, 1.0 / len);
    Some((v, len))
}

/// Derivatives of `e = g²`, `g = ⟨v, a⟩`, `v` the circumcenter of `a, b, c` (SPEC §3.5):
/// `∂e/∂q = 2g [ [q = a] v + ((e_q × v)·a / |n|) v ]` with `e_a = c − b`, `e_b = a − c`,
/// `e_c = b − a`. Returns `[∂e/∂a, ∂e/∂b, ∂e/∂c]`.
///
/// [formula: dev-max-covariance-gradient]
pub fn circumcenter_value_gradient(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    v: [f64; 3],
    n_len: f64,
) -> [[f64; 3]; 3] {
    let g = dot(v, a);
    let edges = [sub(c, b), sub(a, c), sub(b, a)];
    let mut out = [[0.0; 3]; 3];
    for (q, e) in edges.into_iter().enumerate() {
        let mut coef = dot(cross(e, v), a) / n_len;
        if q == 0 {
            coef += 1.0;
        }
        out[q] = scale(v, 2.0 * g * coef);
    }
    out
}

/// Max covariance energy `e_v = min_{|u|=1} max_σ ⟨u, N_σ⟩²` (Eq. 6), the largest empty circle
/// among the sites `±N_σ`, by the candidates of SPEC §3.5 (D11):
///
/// 1. for face triples `σ1 < σ2 < σ3` and signs `(s2, s3)` in the order `++, +−, −+, −−`: the
///    circumcenter of `(N_σ1, s2 N_σ2, s3 N_σ3)`;
/// 2. for pairs `σ1 < σ2`: the circumcenter of `(N_σ1, −N_σ1, N_σ2)` (the pole `∥ N_σ1 × N_σ2`).
///
/// Candidates with an exactly zero (or non-finite) `n` are skipped. The minimum of `φ` is attained
/// at a spherical Voronoi vertex of the sites, which is a candidate (proof sketch in SPEC §3.5;
/// every unit vector is an upper bound, so extra candidates are harmless). The chosen candidate
/// minimizes `(φ(v), φ(v) − ⟨v,a⟩², enumeration index)`. No candidate (all normals parallel or
/// antiparallel) gives `e_v = 0`.
///
/// Subgradient: `y(σ_q) += s_q ∂(g²)/∂q` for the three sites `q = s_q N_σq` of the chosen
/// candidate ([`circumcenter_value_gradient`]; the site sign `s_q` is missing in the paper, D5),
/// `s_c = 0`. Margin: second-smallest minus smallest `φ` over the candidates.
///
/// [formula: dev-max-covariance-voronoi]
pub fn covariance_max(star: &Star, y: &mut [[f64; 3]], s: &mut [f64]) -> StarOut {
    let nrm = star.normals;
    let k = nrm.len();
    s.fill(0.0);
    y.fill([0.0; 3]);
    let phi = |v: [f64; 3]| {
        let mut m = 0.0f64;
        for n in nrm {
            let t = dot(v, *n);
            m = m.max(t * t);
        }
        m
    };
    // Candidate: (faces, signs) of the sites (a, b, c).
    let mut best: Option<([usize; 3], [f64; 3], [f64; 3], f64)> = None; // faces, signs, v, |n|
    let mut best_key = (f64::INFINITY, f64::INFINITY);
    let mut second = f64::INFINITY;
    let mut consider = |faces: [usize; 3], signs: [f64; 3]| {
        let a = scale(nrm[faces[0]], signs[0]);
        let b = scale(nrm[faces[1]], signs[1]);
        let c = scale(nrm[faces[2]], signs[2]);
        let Some((v, len)) = spherical_circumcenter(a, b, c) else {
            return;
        };
        let p = phi(v);
        let g = dot(v, a);
        let slack = p - g * g;
        if p < best_key.0 || (p == best_key.0 && slack < best_key.1) {
            if p < best_key.0 {
                second = best_key.0;
            } else {
                second = second.min(p);
            }
            best_key = (p, slack);
            best = Some((faces, signs, v, len));
        } else {
            second = second.min(p);
        }
    };
    for i in 0..k {
        for j in i + 1..k {
            for l in j + 1..k {
                for (s2, s3) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
                    consider([i, j, l], [1.0, s2, s3]);
                }
            }
        }
    }
    for i in 0..k {
        for j in i + 1..k {
            consider([i, i, j], [1.0, -1.0, 1.0]);
        }
    }
    let Some((faces, signs, v, len)) = best else {
        return StarOut { value: 0.0, margin: f64::INFINITY, z: [0.0; 3] };
    };
    let sites = [0, 1, 2].map(|q| scale(nrm[faces[q]], signs[q]));
    let grad = circumcenter_value_gradient(sites[0], sites[1], sites[2], v, len);
    for q in 0..3 {
        y[faces[q]] = axpy(y[faces[q]], signs[q], grad[q]);
    }
    StarOut { value: best_key.0, margin: second - best_key.0, z: [0.0; 3] }
}
