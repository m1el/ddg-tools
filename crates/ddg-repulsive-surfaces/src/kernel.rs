//! The exponents and the tangent-point kernel of surfaces with its gradient (`SPEC.md` §2.2,
//! §3.1, §3.2; paper Eq. 9, 16).
//!
//! ```text
//! K(X, N; Y) = |N·(X − Y)|^α |X − Y|^{−β}
//! ```
//!
//! With `v = X − Y`, `ν = N·v`, `r = |v|`:
//!
//! ```text
//! ∂K/∂v = α|ν|^{α−2}ν r^{−β} N − β|ν|^α r^{−β−2} v,   ∂K/∂N = α|ν|^{α−2}ν r^{−β} v,   ∂K/∂Y = −∂K/∂v
//! ```
//!
//! `|ν|^{α−2}ν = sign(ν)|ν|^{α−1}` is continuous at `ν = 0` for `α > 1` (and 0 there). Integer
//! exponents (`α`, `β` integers up to 64) use exact repeated multiplication (binary powering),
//! the others `powf`; both are deterministic.

use crate::error::RsError;
use crate::geometry::v3::{axpy, dot, scale};

/// The exponent pair `(α, β)` of the tangent-point energy (`SPEC.md` §2.2, D3).
///
/// Invariant (checked by [`Exponents::new`]): `α > 1` finite and `s = (β − 2)/α ∈ (1, 2)`,
/// i.e. `α + 2 < β < 2α + 2`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Exponents {
    alpha: f64,
    beta: f64,
    /// `α` as an integer when it is one (≤ 64).
    alpha_int: Option<u32>,
    /// `β` as an integer when it is one (≤ 128).
    beta_int: Option<u32>,
}

fn as_small_int(x: f64, max: f64) -> Option<u32> {
    (x == x.trunc() && x > 0.0 && x <= max).then_some(x as u32)
}

impl Exponents {
    /// Validates `(α, β)`: finite, `α > 1`, `α + 2 < β < 2α + 2` (`s ∈ (1, 2)`). The
    /// subcritical range `β − α ≤ 4` is allowed (paper Fig. 26) and reported by
    /// [`Exponents::subcritical`]. Nothing is clamped.
    ///
    /// [formula: rs-exponents]
    /// [formula: rs-exponent-range]
    ///
    /// # Errors
    /// [`RsError::InvalidExponents`] naming the violated requirement.
    pub fn new(alpha: f64, beta: f64) -> Result<Exponents, RsError> {
        let bad = |requirement: &'static str| {
            Err(RsError::InvalidExponents {
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
        if beta <= alpha + 2.0 {
            return bad("beta must be greater than alpha + 2 (s = (beta - 2)/alpha > 1)");
        }
        if beta >= 2.0 * alpha + 2.0 {
            return bad("beta must be smaller than 2 alpha + 2 (s = (beta - 2)/alpha < 2)");
        }
        Ok(Exponents {
            alpha,
            beta,
            alpha_int: as_small_int(alpha, 64.0),
            beta_int: as_small_int(beta, 128.0),
        })
    }

    /// The paper's family `α = p`, `β = 2p` (`p > 2`; `s = 2 − 2/p`).
    ///
    /// [formula: rs-exponents]
    ///
    /// # Errors
    /// As [`Exponents::new`].
    pub fn from_p(p: f64) -> Result<Exponents, RsError> {
        Exponents::new(p, 2.0 * p)
    }

    /// `α`.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// `β`.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// `s = (β − 2)/α ∈ (1, 2)`: the differential has order `2s` (`= 2 − 2/p` for `β = 2p`).
    ///
    /// [formula: rs-exponents]
    pub fn s(&self) -> f64 {
        (self.beta - 2.0) / self.alpha
    }

    /// `σ = s − 1 ∈ (0, 1)`.
    ///
    /// [formula: rs-exponents]
    pub fn sigma(&self) -> f64 {
        self.s() - 1.0
    }

    /// The kernel exponent of the high- and low-order terms `B`, `B⁰`: `2s = 2σ + 2`
    /// (`ρ^{−2s}`).
    ///
    /// [formula: rs-exponents]
    pub fn high_order_kernel_exponent(&self) -> f64 {
        2.0 * self.sigma() + 2.0
    }

    /// The kernel exponent of the preconditioner's `L_{2−s}`: `6 − 2s = 2(2 − s) + 2`.
    ///
    /// [formula: rs-exponents]
    pub fn preconditioner_kernel_exponent(&self) -> f64 {
        2.0 * (2.0 - self.s()) + 2.0
    }

    /// `β − α ≤ 4`: below the barrier threshold (paper Fig. 26; `p ≤ 4` for `β = 2p`). The
    /// energy is then finite on intersecting surfaces.
    ///
    /// [formula: rs-exponent-range]
    pub fn subcritical(&self) -> bool {
        self.beta - self.alpha <= 4.0
    }

    /// `4 + α − β`: `ℰ̂(cf) = c^{4+α−β} ℰ̂(f)` for `c > 0` (`SPEC.md` §3.3).
    ///
    /// [formula: rs-energy-scaling]
    pub fn energy_scaling(&self) -> f64 {
        4.0 + self.alpha - self.beta
    }

    /// `(|ν|^α, sign(ν)|ν|^{α−1})`.
    #[inline(always)]
    fn nu_powers(&self, nu: f64) -> (f64, f64) {
        let a = nu.abs();
        let p1 = match self.alpha_int {
            Some(k) => ipow(a, k - 1),
            None => a.powf(self.alpha - 1.0),
        };
        (p1 * a, p1.copysign(nu))
    }

    /// `r^{−β}` from `r² `.
    #[inline(always)]
    fn r_power(&self, r2: f64) -> f64 {
        match self.beta_int {
            Some(k) if k % 2 == 0 => 1.0 / ipow(r2, k / 2),
            Some(k) => 1.0 / (ipow(r2, k / 2) * r2.sqrt()),
            None => r2.powf(-0.5 * self.beta),
        }
    }
}

impl Default for Exponents {
    /// `(6, 12)`, the paper's `p = 6`.
    fn default() -> Self {
        Exponents {
            alpha: 6.0,
            beta: 12.0,
            alpha_int: Some(6),
            beta_int: Some(12),
        }
    }
}

/// `x^n` by binary powering (exact operation order, independent of the platform's `powi`).
#[inline(always)]
pub fn ipow(mut x: f64, mut n: u32) -> f64 {
    let mut acc = 1.0;
    while n > 0 {
        if n & 1 == 1 {
            acc *= x;
        }
        n >>= 1;
        if n > 0 {
            x *= x;
        }
    }
    acc
}

/// `K(X, N; Y) = |N·(X − Y)|^α |X − Y|^{−β}` (Eq. 16; `N` a unit normal, its sign is
/// irrelevant). `X = Y` gives a non-finite value.
///
/// [formula: rs-kernel]
#[inline]
pub fn kernel(e: &Exponents, x: [f64; 3], n: [f64; 3], y: [f64; 3]) -> f64 {
    kernel_of_v(e, crate::geometry::v3::sub(x, y), n)
}

/// [`kernel`] as a function of `v = X − Y`.
///
/// [formula: rs-kernel]
#[inline(always)]
pub fn kernel_of_v(e: &Exponents, v: [f64; 3], n: [f64; 3]) -> f64 {
    let (pa, _) = e.nu_powers(dot(n, v));
    pa * e.r_power(dot(v, v))
}

/// The kernel and its partial derivatives (`SPEC.md` §3.2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KernelEval {
    /// `K`.
    pub value: f64,
    /// `∂K/∂v = ∂K/∂X = −∂K/∂Y`.
    pub d_v: [f64; 3],
    /// `∂K/∂N` (the unconstrained partial; a unit-normal variation is applied by the caller).
    pub d_n: [f64; 3],
}

/// `K(X, N; Y)` with `∂K/∂X = d_v`, `∂K/∂Y = −d_v`, `∂K/∂N = d_n`.
///
/// [formula: rs-kernel-gradient]
#[inline]
pub fn kernel_with_gradient(e: &Exponents, x: [f64; 3], n: [f64; 3], y: [f64; 3]) -> KernelEval {
    kernel_with_gradient_of_v(e, crate::geometry::v3::sub(x, y), n)
}

/// [`kernel_with_gradient`] as a function of `v = X − Y`.
///
/// [formula: rs-kernel-gradient]
#[inline(always)]
pub fn kernel_with_gradient_of_v(e: &Exponents, v: [f64; 3], n: [f64; 3]) -> KernelEval {
    let r2 = dot(v, v);
    let rb = e.r_power(r2);
    let (pa, sp1) = e.nu_powers(dot(n, v));
    let value = pa * rb;
    let c_n = e.alpha * sp1 * rb;
    let d_v = axpy(scale(n, c_n), -e.beta * value / r2, v);
    KernelEval {
        value,
        d_v,
        d_n: scale(v, c_n),
    }
}

/// Kernel values and `v`-partials of both directions of a face pair with a shared `r^{−β}`:
/// for `v = X_S − X_T`, returns `(K(S→T), K(T→S), ∂K(S→T)/∂v, ∂K(T→S)/∂v', ∂K(S→T)/∂N_S)`
/// with `v' = −v` (crate-internal; used by the exact energy).
#[inline(always)]
pub(crate) fn pair_terms(
    e: &Exponents,
    v: [f64; 3],
    n_s: [f64; 3],
    n_t: [f64; 3],
) -> (f64, f64, [f64; 3], [f64; 3], [f64; 3]) {
    let r2 = dot(v, v);
    let rb = e.r_power(r2);
    let (pa_s, sp_s) = e.nu_powers(dot(n_s, v));
    // ν_T = N_T·(−v).
    let (pa_t, sp_t) = e.nu_powers(-dot(n_t, v));
    let k_st = pa_s * rb;
    let k_ts = pa_t * rb;
    let c_s = e.alpha * sp_s * rb;
    let c_t = e.alpha * sp_t * rb;
    let inv_r2 = 1.0 / r2;
    let dv_st = axpy(scale(n_s, c_s), -e.beta * k_st * inv_r2, v);
    let dv_ts = axpy(scale(n_t, c_t), e.beta * k_ts * inv_r2, v);
    (k_st, k_ts, dv_st, dv_ts, scale(v, c_s))
}
