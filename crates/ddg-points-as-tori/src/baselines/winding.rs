//! Winding numbers of oriented point clouds by direct summation (SPEC.md §9).
//!
//! * The generalized winding number of a cloud (Barill et al. 2018, Eq. 5; Jacobson et al. 2013):
//!   `w(x) = Σ A_i ⟨n_i, p_i − x⟩/(4π‖p_i − x‖³)`, about 1 inside a closed outward-oriented
//!   surface and 0 outside. Terms with `p_i = x` are left out and counted.
//! * The regularized winding number (Chen et al. 2024): each dipole kernel multiplied by
//!   `S(t) = erf t − (2t/√π) e^{−t²}`, `t = ρ/ε`; smooth everywhere (the term vanishes at
//!   `ρ = 0`). The paper's Fig. 8 uses `ε = 0.02`.
//! * The winding-number-signed distance: `|d(x)|` with the sign of `½ − w(x)` (negative inside),
//!   the "GWN sign + distance" method of SPEC §16.2.
//!
//! `S(t)` is the regularized lower incomplete gamma function `P(3/2, t²)`, evaluated by its power
//! series for `t² < 5/2` and by the Legendre continued fraction of `Q = 1 − P` otherwise, so that
//! `S(t)/ρ³` has no cancellation at small `t`.

use super::{CHUNK, check_len, check_queries, dot, norm, sub};
use crate::{PatError, PointCloud};
use rayon::prelude::*;
use std::f64::consts::PI;

/// What a batch evaluation did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct WindingReport {
    /// Queries.
    pub queries: usize,
    /// Terms with `p_i = x` left out (direct winding number only).
    pub skipped: u64,
    /// Wall time.
    pub seconds: f64,
}

fn check_areas(cloud: &PointCloud, areas: &[f64]) -> Result<(), PatError> {
    if areas.len() != cloud.len() {
        return Err(PatError::SizeMismatch { what: format!("{} areas for {} points", areas.len(), cloud.len()) });
    }
    if let Some(i) = areas.iter().position(|a| !a.is_finite()) {
        return Err(PatError::NonFinite { index: i, what: "point area" });
    }
    Ok(())
}

/// The point-cloud winding number at `x` and the number of left-out terms (`p_i = x`).
///
/// [formula: pat-winding-dipole]
pub fn winding_number(positions: &[[f64; 3]], normals: &[[f64; 3]], areas: &[f64], x: [f64; 3]) -> (f64, u32) {
    let mut w = 0.0;
    let mut skipped = 0;
    for i in 0..positions.len() {
        let d = sub(positions[i], x);
        let r2 = dot(d, d);
        if r2 == 0.0 {
            skipped += 1;
            continue;
        }
        let r = r2.sqrt();
        w += areas[i] * dot(normals[i], d) / (4.0 * PI * r2 * r);
    }
    (w, skipped)
}

/// `w` at every query on an explicit pool (`threads`; `None`: all cores).
pub fn winding_numbers(
    cloud: &PointCloud,
    areas: &[f64],
    xs: &[[f64; 3]],
    out: &mut [f64],
    threads: Option<usize>,
) -> Result<WindingReport, PatError> {
    let t0 = std::time::Instant::now();
    check_areas(cloud, areas)?;
    check_queries(xs)?;
    check_len("winding", out.len(), xs.len())?;
    let pool = crate::pool::make_pool(threads)?;
    let (p, n) = (cloud.positions(), cloud.normals());
    let skipped: u64 = pool.install(|| {
        out.par_chunks_mut(CHUNK)
            .zip(xs.par_chunks(CHUNK))
            .map(|(o, q)| {
                let mut s = 0u64;
                for (o, &x) in o.iter_mut().zip(q) {
                    let (w, k) = winding_number(p, n, areas, x);
                    *o = w;
                    s += k as u64;
                }
                s
            })
            .collect::<Vec<u64>>()
            .iter()
            .sum()
    });
    Ok(WindingReport { queries: xs.len(), skipped, seconds: t0.elapsed().as_secs_f64() })
}

const SQRT_PI: f64 = 1.772_453_850_905_516;

/// `Σ_{k≥0} x^k / ((a+1)(a+2)…(a+k))` for `a = 3/2` (the series of `P(3/2, x)` without its
/// prefactor `x^{3/2} e^{−x}/Γ(5/2)`).
fn p32_series(x: f64) -> f64 {
    let a = 1.5;
    let mut term = 1.0;
    let mut sum = 1.0;
    for k in 1..400 {
        term *= x / (a + k as f64);
        sum += term;
        if term < sum * 1e-17 {
            break;
        }
    }
    sum
}

/// `Q(3/2, x) = Γ(3/2, x)/Γ(3/2)` by the Legendre continued fraction (modified Lentz), `x ≥ 5/2`.
fn q32_cf(x: f64) -> f64 {
    let a = 1.5;
    let tiny = 1e-300;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / tiny;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..400 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < tiny {
            d = tiny;
        }
        c = b + an / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-16 {
            break;
        }
    }
    // Γ(3/2) = √π/2
    (-x + a * x.ln()).exp() * h / (0.5 * SQRT_PI)
}

/// The regularized dipole factor `K_ε(ρ) = S(ρ/ε)/ρ³`, `S(t) = erf t − (2t/√π)e^{−t²}`
/// (`= 4/(3√π ε³)` at `ρ = 0`). `ε > 0`.
///
/// [formula: pat-regularized-wn]
pub fn regularized_kernel(rho: f64, eps: f64) -> f64 {
    let t = rho / eps;
    let x = t * t;
    if x < 2.5 {
        // S(t) = t³ e^{−t²} Σ / Γ(5/2), Γ(5/2) = 3√π/4
        (-x).exp() * p32_series(x) / (0.75 * SQRT_PI) / (eps * eps * eps)
    } else {
        (1.0 - q32_cf(x)) / (rho * rho * rho)
    }
}

/// `S(t) = erf t − (2t/√π) e^{−t²}` for `t ≥ 0`.
///
/// [formula: pat-regularized-wn]
pub fn regularization_s(t: f64) -> f64 {
    let x = t * t;
    if x < 2.5 { t * x * (-x).exp() * p32_series(x) / (0.75 * SQRT_PI) } else { 1.0 - q32_cf(x) }
}

/// The regularized winding number at `x`: `Σ A_i ⟨n_i, p_i − x⟩ K_ε(ρ_i)/(4π)`.
///
/// [formula: pat-regularized-wn]
pub fn regularized_winding_number(positions: &[[f64; 3]], normals: &[[f64; 3]], areas: &[f64], x: [f64; 3], eps: f64) -> f64 {
    let mut w = 0.0;
    for i in 0..positions.len() {
        let d = sub(positions[i], x);
        w += areas[i] * dot(normals[i], d) * regularized_kernel(norm(d), eps) / (4.0 * PI);
    }
    w
}

/// The regularized winding number at every query (`ε > 0`).
pub fn regularized_winding_numbers(
    cloud: &PointCloud,
    areas: &[f64],
    eps: f64,
    xs: &[[f64; 3]],
    out: &mut [f64],
    threads: Option<usize>,
) -> Result<WindingReport, PatError> {
    let t0 = std::time::Instant::now();
    if !(eps.is_finite() && eps > 0.0) {
        return Err(PatError::InvalidOption { what: format!("regularized winding number: ε = {eps} (must be > 0)") });
    }
    check_areas(cloud, areas)?;
    check_queries(xs)?;
    check_len("winding", out.len(), xs.len())?;
    let pool = crate::pool::make_pool(threads)?;
    let (p, n) = (cloud.positions(), cloud.normals());
    pool.install(|| {
        out.par_chunks_mut(CHUNK).zip(xs.par_chunks(CHUNK)).for_each(|(o, q)| {
            for (o, &x) in o.iter_mut().zip(q) {
                *o = regularized_winding_number(p, n, areas, x, eps);
            }
        })
    });
    Ok(WindingReport { queries: xs.len(), skipped: 0, seconds: t0.elapsed().as_secs_f64() })
}

/// `|d|` signed by the winding number: negative where `w > threshold` (inside), `|d|` otherwise
/// (`threshold = ½` for the smooth theory).
pub fn winding_signed_distance(unsigned: f64, w: f64, threshold: f64) -> f64 {
    if w > threshold { -unsigned.abs() } else { unsigned.abs() }
}
