//! The generalized tangent-point kernel and its gradients (`SPEC.md` §3.1–3.2, Eq. 3).
//!
//! ```text
//! k(p, q, T) = |T × d|^α / |d|^β,    d = p − q,    |T × d| = |P_T d|,    P_T = I − T Tᵀ
//! ```
//!
//! The exponents are validated once, in [`Exponents::new`] (D24).

use crate::error::RcError;

/// Small 3-vector helpers shared by the modules of the crate.
pub mod v3 {
    /// `a + b`.
    #[inline(always)]
    pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }
    /// `a - b`.
    #[inline(always)]
    pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }
    /// `s a`.
    #[inline(always)]
    pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
        [a[0] * s, a[1] * s, a[2] * s]
    }
    /// `a + s b`.
    #[inline(always)]
    pub fn axpy(a: [f64; 3], s: f64, b: [f64; 3]) -> [f64; 3] {
        [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]]
    }
    /// `a · b`.
    #[inline(always)]
    pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
    /// `a × b`.
    #[inline(always)]
    pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }
    /// `|a|`.
    #[inline(always)]
    pub fn norm(a: [f64; 3]) -> f64 {
        dot(a, a).sqrt()
    }
}

use v3::{axpy, dot, scale, sub};

/// The exponent pair `(α, β)` of the energy `ℰ^α_β` (Eq. 3).
///
/// Invariant (checked by [`Exponents::new`]): `α > 1` and `α + 2 ≤ β < 2α + 1`, both finite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Exponents {
    alpha: f64,
    beta: f64,
}

impl Exponents {
    /// Validates `(α, β)`: the energy is finite on embedded `C¹` curves of the energy space
    /// `W^{s,α}` exactly when `α > 1` and `β ∈ [α + 2, 2α + 1)` (Blatt 2013), `s = (β − 1)/α`.
    /// Nothing is clamped (D24).
    ///
    /// [formula: rc-exponent-range]
    ///
    /// # Errors
    /// [`RcError::InvalidExponents`] naming the violated requirement.
    pub fn new(alpha: f64, beta: f64) -> Result<Exponents, RcError> {
        let bad = |requirement: &'static str| {
            Err(RcError::InvalidExponents {
                alpha,
                beta,
                requirement,
            })
        };
        if !alpha.is_finite() || !beta.is_finite() {
            return bad("both exponents must be finite");
        }
        if alpha <= 1.0 {
            return bad("alpha must be greater than 1");
        }
        if beta < alpha + 2.0 {
            return bad("beta must be at least alpha + 2");
        }
        if beta >= 2.0 * alpha + 1.0 {
            return bad("beta must be smaller than 2 alpha + 1");
        }
        Ok(Exponents { alpha, beta })
    }

    /// `α`.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// `β`.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// `s = (β − 1)/α ∈ (1, 2)`, the Sobolev–Slobodeckij order of the energy space.
    pub fn s(&self) -> f64 {
        (self.beta - 1.0) / self.alpha
    }

    /// `σ = s − 1 ∈ (0, 1)`, the fractional order of the metric of SPEC §4.
    pub fn sigma(&self) -> f64 {
        self.s() - 1.0
    }

    /// `2 + α − β`: the energy scales as `c^{2+α−β}` under `γ → cγ` (SPEC §3.3).
    pub fn energy_scaling(&self) -> f64 {
        2.0 + self.alpha - self.beta
    }
}

impl Default for Exponents {
    /// `(3, 6)`, the paper's and the reference's default.
    fn default() -> Self {
        Exponents {
            alpha: 3.0,
            beta: 6.0,
        }
    }
}

/// `x^e` with `powi` for small integer exponents (exact repeated multiplication, much faster
/// than `powf`). `x` must be `>= 0` when `e` is not an integer.
#[inline(always)]
fn pw(x: f64, e: f64) -> f64 {
    if e == e.trunc() && e.abs() <= 32.0 {
        x.powi(e as i32)
    } else {
        x.powf(e)
    }
}

/// Value and gradients of the kernel at `(p, q, T)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KernelEval {
    /// `k(p, q, T)`.
    pub value: f64,
    /// `∇_p k`; `∇_q k = −∇_p k`.
    pub grad_p: [f64; 3],
    /// `∇_T k` projected onto `T^⊥` (the part unit-`T` variations see).
    pub grad_t: [f64; 3],
}

/// `k(p, q, T) = |T × (p − q)|^α / |p − q|^β` for a unit tangent `T` (Eq. 3). Infinite if
/// `p = q`; the caller (network validation) excludes that.
///
/// [formula: rc-kernel]
#[inline]
pub fn kernel(p: [f64; 3], q: [f64; 3], t: [f64; 3], e: Exponents) -> f64 {
    let d = sub(p, q);
    let td = dot(t, d);
    let n = axpy(d, -td, t);
    let nu = dot(n, n).sqrt();
    let rho = dot(d, d).sqrt();
    pw(nu, e.alpha) * pw(rho, -e.beta)
}

/// The kernel together with its gradients (`SPEC.md` §3.2):
///
/// ```text
/// ∇_p k = α ν^{α−1} n̂ / ρ^β − β k d / ρ²,   ∇_q k = −∇_p k,   ∇_T k = −α ν^{α−1} (T·d) n̂ / ρ^β
/// ```
///
/// with `n = P_T d`, `ν = |n|`, `ρ = |d|`. The term `ν^{α−1} n̂` is `0` at `ν = 0`, the
/// continuous extension for `α > 1` (D24); no cutoff is applied.
///
/// [formula: rc-kernel-gradient]
#[inline]
pub fn kernel_with_gradient(p: [f64; 3], q: [f64; 3], t: [f64; 3], e: Exponents) -> KernelEval {
    let d = sub(p, q);
    let rho2 = dot(d, d);
    let td = dot(t, d);
    let n = axpy(d, -td, t);
    let nu = dot(n, n).sqrt();
    let rho = rho2.sqrt();
    let rho_mb = pw(rho, -e.beta);
    let value = pw(nu, e.alpha) * rho_mb;
    // a = α ν^{α−1} n̂ / ρ^β, computed without dividing by ν where possible.
    let a = if nu == 0.0 {
        [0.0; 3]
    } else if e.alpha >= 2.0 {
        scale(n, e.alpha * pw(nu, e.alpha - 2.0) * rho_mb)
    } else {
        scale(n, e.alpha * pw(nu, e.alpha - 1.0) / nu * rho_mb)
    };
    let grad_p = axpy(a, -e.beta * value / rho2, d);
    let grad_t = scale(a, -td);
    KernelEval {
        value,
        grad_p,
        grad_t,
    }
}
