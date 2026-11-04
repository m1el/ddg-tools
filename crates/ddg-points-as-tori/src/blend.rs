//! The blend of the tori's SDFs and its derivatives (SPEC.md §7.4; paper Eq. 25–26, App. D; the
//! reference's `evaluate_distance_gradient_laplacian_single`, points-as-tori, MIT, Copyright (c)
//! 2026 Nicole Feng). f64.
//!
//! For a query `x` and a neighbor set `𝒩`: `ρ_i = ‖x − p_i‖`, `σ_x = ½ max_𝒩 ρ_i`, `λ` by
//! [`LambdaRule`], `w_i = A_i exp(−λ(ρ_i − s))` (`A_i = 1` without area weights) and
//! `φ = Σ w_i g_i / Σ w_i − isovalue`. Sums run over `𝒩` in its given order.
//!
//! The exponent shift `s` cancels in `φ` and its derivatives. It is `σ_x` for
//! [`LambdaRule::PerQuery`] (the reference) and `min_𝒩 ρ` for the constant rules: with
//! `λ = 10³/D` and a 32-NN fallback far from the cloud, `λσ_x` exceeds 709 and `exp` overflows
//! (SPEC resolution from WP3).

use crate::sdf::{Singular, torus_sdf, torus_sdf_gradient, torus_sdf_laplacian};
use crate::Torus;

/// How `λ` is chosen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LambdaRule {
    /// `λ = c/σ_x` per query (the released reference; `c = 64`).
    PerQuery {
        /// The constant `C`.
        c: f64,
    },
    /// A fixed `λ`.
    Constant {
        /// The value.
        lambda: f64,
    },
    /// Sec. 4.2 as printed: `λ = factor/D` with `D` the mean distance of every point to its
    /// `k_mean` nearest points (self included, as the reference's helper), constant over queries.
    /// Without the reference's dead `λ·R/√3` factor. Its evaluation radius `R_eval = 2C/λ` is set
    /// through [`crate::Neighborhood::Radius`] (see [`crate::TorusField::paper_options`]).
    PaperHeuristic {
        /// The constant `C` (for `R_eval = 2C/λ`).
        c: f64,
        /// The numerator (1e3).
        factor: f64,
        /// Neighbors in the mean distance `D` (64).
        k_mean: usize,
    },
}

impl Default for LambdaRule {
    fn default() -> Self {
        LambdaRule::PerQuery { c: 64.0 }
    }
}

/// Which gradient is returned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GradientMode {
    /// App. D with `λ` and `σ` held fixed (the reference).
    #[default]
    Frozen,
    /// The exact gradient of `φ` where the neighbor set and the farthest neighbor are locally
    /// constant: adds the `∂λ/∂x` term for [`LambdaRule::PerQuery`].
    Exact,
}

/// Counters of one blend evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BlendFlags {
    /// Torus terms evaluated on an axis.
    pub singular_axis: u32,
    /// Torus terms evaluated on a core circle.
    pub singular_core: u32,
    /// Data points coinciding with the query (`r̂ := 0`).
    pub at_data_point: u32,
}

/// Value, gradient and Laplacian of one blend (gradient and Laplacian zero when not requested).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BlendValue {
    /// `φ(x) − isovalue`.
    pub distance: f64,
    /// `∇φ`.
    pub gradient: [f64; 3],
    /// `Δφ`.
    pub laplacian: f64,
    /// `λ` used.
    pub lambda: f64,
    /// `σ_x = ½ max ρ`.
    pub sigma: f64,
    /// The exponent shift `s` (`σ_x` for PerQuery, `min ρ` for the constant rules).
    pub shift: f64,
    /// Counters.
    pub flags: BlendFlags,
}

/// Inputs shared by all queries.
#[derive(Debug, Clone, Copy)]
pub struct BlendInputs<'a> {
    /// Data positions.
    pub positions: &'a [[f64; 3]],
    /// One torus per data point.
    pub tori: &'a [Torus],
    /// Optional per-point weights `A_i`.
    pub area_weights: Option<&'a [f64]>,
    /// The λ rule.
    pub lambda: LambdaRule,
    /// For [`LambdaRule::PaperHeuristic`]: the precomputed constant λ.
    pub lambda_global: f64,
    /// Subtracted from `φ`.
    pub isovalue: f64,
}

/// What to compute.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendWant {
    /// `φ` only.
    Distance,
    /// `φ` and `∇φ` (with the given mode).
    Gradient(GradientMode),
    /// `φ`, the frozen `∇φ` and `Δφ` (App. D divergence, frozen `λ`).
    Laplacian,
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `σ_x = ½ max_𝒩 ‖x − p_i‖`, `min_𝒩 ‖x − p_i‖`, and the position in `nbrs` of the member
/// with the largest `(ρ², index)` (the last member of the `(ρ², index)`-sorted set; its `r̂`
/// enters the exact gradient).
pub fn shift(x: [f64; 3], positions: &[[f64; 3]], nbrs: &[u32]) -> (f64, f64, usize) {
    let mut m2 = f64::NEG_INFINITY;
    let mut min2 = f64::INFINITY;
    let mut arg = 0;
    for (t, &i) in nbrs.iter().enumerate() {
        let r = sub(x, positions[i as usize]);
        let d2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        if d2 > m2 || (d2 == m2 && i > nbrs[arg]) {
            m2 = d2;
            arg = t;
        }
        min2 = min2.min(d2);
    }
    (0.5 * m2.sqrt(), min2.sqrt(), arg)
}

/// `w_i = A_i exp(−λ(ρ_i − σ))`. With `λ = C/σ` and `0 ≤ ρ ≤ 2σ` the exponent lies in `[−C, C]`.
///
/// [formula: pat-blend-weight-range]
#[inline]
pub fn blend_weight(lambda: f64, rho: f64, shift: f64, area: f64) -> f64 {
    area * (-lambda * (rho - shift)).exp()
}

/// The neighbor set degenerates: [`LambdaRule::PerQuery`] with `σ_x = 0` (every member at `x`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Degenerate;

/// One blend over `nbrs` (non-empty).
///
/// [formula: pat-blend-shift]
/// [formula: pat-blend-convex]
/// [formula: pat-blend-gradient]
/// [formula: pat-blend-gradient-exact]
/// [formula: pat-blend-laplacian]
pub fn blend(inp: &BlendInputs<'_>, x: [f64; 3], nbrs: &[u32], want: BlendWant) -> Result<BlendValue, Degenerate> {
    let (sigma, rho_min, arg) = shift(x, inp.positions, nbrs);
    let (lambda, shift) = match inp.lambda {
        LambdaRule::PerQuery { c } => {
            if sigma == 0.0 {
                return Err(Degenerate);
            }
            (c / sigma, sigma)
        }
        LambdaRule::Constant { lambda } => (lambda, rho_min),
        LambdaRule::PaperHeuristic { .. } => (inp.lambda_global, rho_min),
    };
    let mut flags = BlendFlags::default();
    let mut wsum = 0.0;
    let mut gsum = 0.0;
    match want {
        BlendWant::Distance => {
            for &i in nbrs {
                let i = i as usize;
                let rho = norm(sub(x, inp.positions[i]));
                let a = inp.area_weights.map_or(1.0, |a| a[i]);
                let w = blend_weight(lambda, rho, shift, a);
                gsum += w * torus_sdf(&inp.tori[i], x);
                wsum += w;
            }
            Ok(BlendValue { distance: gsum / wsum - inp.isovalue, lambda, sigma, shift, flags, ..Default::default() })
        }
        BlendWant::Gradient(_) | BlendWant::Laplacian => {
            let lap = matches!(want, BlendWant::Laplacian);
            let mut asum = [0.0; 3]; // Σ w (∇g − λ r̂ g)
            let mut nsum = [0.0; 3]; // Σ w r̂
            let mut numer1 = 0.0; // Σ w (Δg − 2λ r̂·∇g + g(λ² − 2λ/ρ))
            let mut numer2 = 0.0; // Σ w (2/ρ − λ)
            let mut esum = 0.0; // Σ w ρ g (exact-gradient term)
            let mut ersum = 0.0; // Σ w ρ
            for &i in nbrs {
                let i = i as usize;
                let rv = sub(x, inp.positions[i]);
                let rho = norm(rv);
                let rhat = if rho > 0.0 {
                    [rv[0] / rho, rv[1] / rho, rv[2] / rho]
                } else {
                    flags.at_data_point += 1;
                    [0.0; 3]
                };
                let a = inp.area_weights.map_or(1.0, |a| a[i]);
                let w = blend_weight(lambda, rho, shift, a);
                let (g, gr, lg, sing) = if lap {
                    torus_sdf_laplacian(&inp.tori[i], x)
                } else {
                    let (g, gr, s) = torus_sdf_gradient(&inp.tori[i], x);
                    (g, gr, 0.0, s)
                };
                match sing {
                    Singular::Axis => flags.singular_axis += 1,
                    Singular::Core => flags.singular_core += 1,
                    Singular::None => {}
                }
                for a in 0..3 {
                    asum[a] += w * (gr[a] - lambda * rhat[a] * g);
                    nsum[a] += w * rhat[a];
                }
                gsum += w * g;
                wsum += w;
                esum += w * rho * g;
                ersum += w * rho;
                if lap {
                    let inv = if rho > 0.0 { 2.0 / rho } else { 0.0 };
                    numer1 += w * (lg - 2.0 * lambda * dot(rhat, gr) + g * lambda * (lambda - inv));
                    numer2 += w * (inv - lambda);
                }
            }
            let a = [asum[0] / wsum, asum[1] / wsum, asum[2] / wsum];
            let b = [lambda * nsum[0] / wsum, lambda * nsum[1] / wsum, lambda * nsum[2] / wsum];
            let phi = gsum / wsum;
            let mut grad = [a[0] + phi * b[0], a[1] + phi * b[1], a[2] + phi * b[2]];
            if let (BlendWant::Gradient(GradientMode::Exact), LambdaRule::PerQuery { c }) = (want, inp.lambda)
                && sigma > 0.0
            {
                // ∂φ/∂λ = Σ w(σ − ρ)(g − φ)/W = −(Σ w ρ g − φ Σ w ρ)/W; ∂λ/∂x = −(C/σ²)·½ r̂_last
                let dphi_dlambda = -(esum - phi * ersum) / wsum;
                let far = sub(x, inp.positions[nbrs[arg] as usize]);
                let lf = norm(far);
                if lf > 0.0 {
                    let k = -(c / (sigma * sigma)) * 0.5 / lf * dphi_dlambda;
                    for t in 0..3 {
                        grad[t] += k * far[t];
                    }
                }
            }
            let laplacian = if lap {
                let first = numer1 / wsum + dot(a, b);
                let second = dot(grad, b) + phi * (lambda * numer2 / wsum + dot(b, b));
                first + second
            } else {
                0.0
            };
            Ok(BlendValue { distance: phi - inp.isovalue, gradient: grad, laplacian, lambda, sigma, shift, flags })
        }
    }
}
