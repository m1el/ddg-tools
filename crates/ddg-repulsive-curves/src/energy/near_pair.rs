//! Near-pair quadrature of the pair integral (WP16; `SPEC.md` §3.4, §17o). **Not in the paper.**
//!
//! The paper's discrete energy (Eq. 17–18, [`super::exact`]) evaluates the pair integral
//! `∫_I ∫_J k(x, y, T_I) dx dy` with the four endpoint pairs only. When two disjoint edges come
//! close in their interiors (a crossing whose closest points are far from all four endpoints),
//! that rule stays bounded while the pair integral of the polygon diverges like
//! `d^{2+α−β}`; the flow can then push two edges together until the contact floor of the
//! line search stops it (`SPEC.md` §17l.6, §17o) [formula: rc-near-divergence]. This module evaluates the same pair integral
//! accurately for close pairs and blends it into the paper's term:
//!
//! ```text
//! e_IJ = e^paper_IJ + w(r_IJ) (e^near_IJ − e^paper_IJ),    r_IJ = d_IJ / ℓ_IJ,
//! ℓ_IJ = √((ℓ_I² + ℓ_J²)/2),   w = 1 for r ≤ c₁, 0 for r ≥ c₂, C² quintic in between
//! ```
//!
//! with `d_IJ` the distance of the two segments. Pairs with `r ≥ c₂` are untouched; when no
//! pair of a curve has `r < c₂` the energy and its differential are bitwise those of the
//! paper's formula.
//!
//! **Inner integral (along `I`).** For a point `y` on `J`, `x = γ_{i₁} + u T_I` and `T_I`
//! fixed, `|T_I × (x − y)| = h(y)` does not depend on `u` [formula: rc-near-h-constant]: it is
//! the distance of `y` from the line of `I`. With `u₀ = T_I·(y − γ_{i₁})` and the substitution
//! `u − u₀ = h tan θ` [formula: rc-near-inner]
//!
//! ```text
//! F(y) = ∫_0^{ℓ_I} h^α ((u − u₀)² + h²)^{−β/2} du = h^{1+α−β} (S(θ₁) − S(θ₀)),
//! θ₁ = atan2(ℓ_I − u₀, h),  θ₀ = atan2(−u₀, h),   S(θ) = ∫_0^θ cos^{β−2} φ dφ.
//! ```
//!
//! `S` is evaluated to rounding by two binomial series, each with ratio at most ½
//! (`|θ| ≤ π/4`: `S = Σ_k C(a,k)(−1)^k sin^{2k+1}θ/(2k+1)`, `a = (β−3)/2`; `|θ| > π/4`:
//! `S = ±(S(π/2) − R(cos θ))`, `R(c) = c^{β−1} Σ_k C(−½,k)(−1)^k c^{2k}/(β−1+2k)`), 60 terms in
//! Horner form [formula: rc-near-cos-series]. So `∂S/∂θ = cos^{β−2}θ` holds to rounding and the
//! differential of `F` is analytic.
//!
//! **Outer integral (along `J`).** `F` peaks at the point of `J` closest to `I`, with width
//! about `d_IJ`, and has steps of the same width where `y` passes the ends of `I` (nearly
//! parallel pairs). `J = [0, 1]` (parameter `t`) is split at the closest-point parameter `t*` and
//! at the projections of the two endpoints of `I` onto `J` (sorted); every panel is halved and each
//! half is graded toward its end `c` by `t = c ± δ sinh w`, `w ∈ [0, asinh(L/δ)]`, `δ = d_IJ/ℓ_J`
//! (Johnston & Elliott 2005), with `M` sub-panels of 8-point Gauss–Legendre in `w`
//! [formula: rc-near-outer]. The node count is fixed (`8 · 8 M` per ordered pair), so the
//! quadrature is a smooth function of the positions wherever the breakpoints are; the energy
//! grows without bound as `d → 0` at bounded cost.
//!
//! **Differential.** Exact for the evaluated quadrature: the nodes and weights are
//! differentiated through `t*` (implicit differentiation of the closest-point conditions),
//! the endpoint projections, `δ` and `ℓ_J`, and `d_IJ` by the envelope theorem
//! (`∇d = n̂ ⊗ ((1−s*), s*, −(1−t*), −t*)`) [formula: rc-near-differential]. `t*`, the
//! projections and the order of the sorted breakpoints are piecewise smooth (clamping,
//! re-ordering, the parallel case); at their switches the quadrature has a kink of the size of
//! its error (about 1e-7 relative), nothing more.

use rayon::prelude::*;

use crate::error::RcError;
use crate::kernel::v3::{add, axpy, dot, norm, scale, sub};
use crate::kernel::{Exponents, kernel_with_gradient};
use crate::network::CurveNetwork;

type P = [f64; 3];
/// A gradient with respect to the four endpoints `(p₀, p₁, q₀, q₁)` of an edge pair.
pub type Grad4 = [[f64; 3]; 4];

const ZERO4: Grad4 = [[0.0; 3]; 4];

/// The 8-point Gauss–Legendre rule on `[−1, 1]` (nodes, weights).
const GL8: [(f64, f64); 8] = [
    (-0.9602898564975362, 0.10122853629037669),
    (-0.7966664774136267, 0.22238103445337434),
    (-0.525532409916329, 0.31370664587788705),
    (-0.18343464249564978, 0.36268378337836177),
    (0.18343464249564978, 0.36268378337836177),
    (0.525532409916329, 0.31370664587788705),
    (0.7966664774136267, 0.22238103445337434),
    (0.9602898564975362, 0.10122853629037669),
];

/// Number of terms of the two series of `S` (ratio at most ½: remainder below `2⁻⁶⁰`).
const SERIES_TERMS: usize = 60;

/// Options of the near-pair quadrature (`SPEC.md` §3.4). Not part of the paper: `None` in
/// [`FlowOptions::near_pair`](crate::step::FlowOptions::near_pair) gives the paper's energy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NearPairOptions {
    /// `c₁`: pairs with `d_IJ ≤ c₁ ℓ_IJ` use the accurate pair integral only.
    pub blend_start: f64,
    /// `c₂`: pairs with `d_IJ ≥ c₂ ℓ_IJ` use the paper's term only.
    pub blend_end: f64,
    /// `M`: Gauss–Legendre sub-panels (8 nodes each) per graded half-panel of the outer integral.
    pub subpanels: usize,
    /// Also correct pairs of edges joined by one edge (`I`, `K`, `J` consecutive at a vertex
    /// each). Default `true` with [`SecondNeighbourRule::fold_gated`] (`SPEC.md` §17p, §17r,
    /// §18 Q10, Q11): the paper's rule cannot see an accordion fold (`I` and `J` overlapping
    /// across `K`), while at a corner, where the closest points are at the ends of `K` and Eq. 18
    /// overestimates (a barrier the paper's flows and the authors' scenes rely on), the one-sided
    /// rule keeps Eq. 18; hairpins are corrected only near contact. `false`: the pair is the
    /// paper's term (WP16, §17o).
    pub second_neighbours: bool,
    /// How second neighbours are corrected when [`NearPairOptions::second_neighbours`] is set
    /// (WP16b, `SPEC.md` §17p). Default [`SecondNeighbourRule::fold_gated`] (never below Eq. 18;
    /// §17r, §18 Q11).
    /// Ignored when second neighbours are excluded.
    pub second_neighbour_rule: SecondNeighbourRule,
    /// The band `(c₁, c₂)` in `r = d/ℓ_IJ` for second neighbours; `None` (the default): the
    /// band of every other pair. Ignored when second neighbours are excluded.
    pub second_neighbour_band: Option<(f64, f64)>,
}

/// The correction of a pair `I, K, J` joined by one edge `K` (WP16b, `SPEC.md` §17p). Every rule
/// keeps the blend weight `w(r)` of §3.4 (with [`NearPairOptions::second_neighbour_band`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SecondNeighbourRule {
    /// `w (E^n − E^p)`, as for every other pair (§3.4). Lowers the energy at corners, where Eq. 18
    /// overestimates.
    TwoSided,
    /// Never below Eq. 18: `w E^p ψ(E^n/E^p − 1)` with `E^{n,p}` the sums over the two orders and
    /// `ψ(z) = 0` for `z ≤ 0`, `6z³ − 8z⁴ + 3z⁵` on `[0, 1]`, `z` beyond (C², nondecreasing,
    /// `0 ≤ ψ(z) ≤ max(z, 0)`). So `E^p ≤ e ≤ max(E^p, E^n)`, and the pair is the paper's wherever
    /// the accurate integral is not larger. [formula: rc-near-onesided]
    OneSided,
    /// `w g (E^n − E^p)` with the gate `g = G(|s* − s_K|/σ) G(|t* − t_K|/σ)`, `s_K`, `t_K` the
    /// parameters of the ends of `I`, `J` at `K`, `G(x) = W(min(x, 1))` the quintic of the blend:
    /// zero when a closest point is at `K` (a corner), one when both are at least `σ` away from it
    /// (a fold). [formula: rc-near-gate]
    Gated {
        /// `σ ∈ (0, 1]`, in segment parameters.
        width: f64,
    },
    /// The gate of [`SecondNeighbourRule::Gated`] on the one-sided term:
    /// `w g E^p ψ(E^n/E^p − 1)`.
    OneSidedGated {
        /// `σ ∈ (0, 1]`, in segment parameters.
        width: f64,
    },
    /// The one-sided term with the "either" gate `g = 1 − (1 − G_s)(1 − G_t)`
    /// (`G_s = G(|s* − s_K|/σ)`, `G_t` likewise): shut only when **both** closest points are at
    /// the ends of `K` (a corner), open as soon as one of them is at least `σ` away. The product
    /// gate of [`SecondNeighbourRule::OneSidedGated`] lets a flow slide the contact onto an end
    /// of `K` (one closest point at `K`, the other in the interior of its edge) and shut the
    /// barrier (`SPEC.md` §17p). [formula: rc-near-gate]
    OneSidedGatedEither {
        /// `σ ∈ (0, 1]`, in segment parameters.
        width: f64,
    },
    /// The one-sided term with the weight `W = 1 − (1 − w_c(r))(1 − w(r) G(cos φ))` in place of
    /// `w(r)` (WP16c, `SPEC.md` §17r): `w` the band of the rule (§3.4, or
    /// [`NearPairOptions::second_neighbour_band`]), `w_c` the same quintic on the narrower
    /// `contact_band`, `φ` the angle between `I` oriented toward `K` and `J` oriented away from
    /// `K` (0 for a Z-fold or a straight line, π for a hairpin), `G` the quintic of the blend
    /// rising from 0 at `cos φ = cos_shut` to 1 at `cos φ = cos_open`. Z-folds are corrected as by
    /// [`SecondNeighbourRule::OneSided`]; hairpins only within the contact band. `0 ≤ W ≤ 1`, so
    /// `E^p ≤ e ≤ max(E^p, E^n)` still holds. The default since §18 Q11, with the parameters of
    /// [`SecondNeighbourRule::fold_gated`]. [formula: rc-near-foldgate]
    OneSidedFoldGated {
        /// The band `(c₁, c₂)` in `r` within which every second-neighbour pair is corrected
        /// (`0 < c₁ < c₂ ≤` the rule's `c₂`).
        contact_band: (f64, f64),
        /// `cos φ` at and below which the fold gate is shut.
        cos_shut: f64,
        /// `cos φ` at and above which the fold gate is open (`cos_shut < cos_open ≤ 1`).
        cos_open: f64,
    },
}

impl SecondNeighbourRule {
    /// The default rule (`SPEC.md` §17r, §18 Q11): contact band `(0.001, 0.005)`, the fold gate shut at
    /// `φ ≥ 90°` and open at `φ ≤ 60°`.
    pub const fn fold_gated() -> SecondNeighbourRule {
        SecondNeighbourRule::OneSidedFoldGated {
            contact_band: (0.001, 0.005),
            cos_shut: 0.0,
            cos_open: 0.5,
        }
    }
}

impl Default for NearPairOptions {
    /// `c₁ = 0.25, c₂ = 0.75, M = 4`, second neighbours corrected by the fold-gated one-sided
    /// rule [`SecondNeighbourRule::fold_gated`] (`SPEC.md` §17o, §17p, §17r, §18 Q10, Q11).
    fn default() -> Self {
        NearPairOptions {
            blend_start: 0.25,
            blend_end: 0.75,
            subpanels: 4,
            second_neighbours: true,
            second_neighbour_rule: SecondNeighbourRule::fold_gated(),
            second_neighbour_band: None,
        }
    }
}

impl NearPairOptions {
    /// Checks `0 < c₁ < c₂ ≤ 4` (finite) and `1 ≤ M ≤ 64`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] naming the field.
    pub fn validate(&self) -> Result<(), RcError> {
        let (c1, c2) = (self.blend_start, self.blend_end);
        if !(c1.is_finite() && c1 > 0.0 && c1 < c2) {
            return Err(RcError::InvalidOption {
                name: "near_pair.blend_start",
                value: c1.to_string(),
                requirement: "finite, > 0 and < blend_end",
            });
        }
        if !(c2.is_finite() && c2 <= 4.0) {
            return Err(RcError::InvalidOption {
                name: "near_pair.blend_end",
                value: c2.to_string(),
                requirement: "finite, > blend_start and ≤ 4",
            });
        }
        if let Some((c1, c2)) = self.second_neighbour_band
            && !(c1.is_finite() && c1 > 0.0 && c1 < c2 && c2.is_finite() && c2 <= 4.0)
        {
            return Err(RcError::InvalidOption {
                name: "near_pair.second_neighbour_band",
                value: format!("({c1}, {c2})"),
                requirement: "finite, 0 < c1 < c2 ≤ 4",
            });
        }
        if let SecondNeighbourRule::Gated { width }
        | SecondNeighbourRule::OneSidedGated { width }
        | SecondNeighbourRule::OneSidedGatedEither { width } = self.second_neighbour_rule
            && !(width > 0.0 && width <= 1.0)
        {
            return Err(RcError::InvalidOption {
                name: "near_pair.second_neighbour_rule.width",
                value: width.to_string(),
                requirement: "in (0, 1]",
            });
        }
        if let SecondNeighbourRule::OneSidedFoldGated {
            contact_band: (c1, c2),
            cos_shut,
            cos_open,
        } = self.second_neighbour_rule
        {
            let wide = self
                .second_neighbour_band
                .map_or(self.blend_end, |(_, c2)| c2);
            if !(c1.is_finite() && c1 > 0.0 && c1 < c2 && c2 <= wide) {
                return Err(RcError::InvalidOption {
                    name: "near_pair.second_neighbour_rule.contact_band",
                    value: format!("({c1}, {c2})"),
                    requirement: "finite, 0 < c1 < c2 ≤ the second neighbours' c2",
                });
            }
            if !(cos_shut.is_finite() && cos_shut >= -1.0 && cos_shut < cos_open && cos_open <= 1.0)
            {
                return Err(RcError::InvalidOption {
                    name: "near_pair.second_neighbour_rule.cos_shut",
                    value: format!("({cos_shut}, {cos_open})"),
                    requirement: "−1 ≤ cos_shut < cos_open ≤ 1",
                });
            }
        }
        if !(1..=64).contains(&self.subpanels) {
            return Err(RcError::InvalidOption {
                name: "near_pair.subpanels",
                value: self.subpanels.to_string(),
                requirement: "between 1 and 64",
            });
        }
        Ok(())
    }
}

// =============================================================================================
// The blend
// =============================================================================================

/// The blend weight `w(r)` and `dw/dr`: with `τ = (c₂ − r)/(c₂ − c₁)` clamped to `[0, 1]`,
/// `w = τ³(10 − 15τ + 6τ²)` (C², `w' = w'' = 0` at both ends), `dw/dr = −30 τ²(1 − τ)²/(c₂ − c₁)`.
///
/// [formula: rc-near-blend]
pub fn blend_weight(r: f64, c1: f64, c2: f64) -> (f64, f64) {
    let tau = (c2 - r) / (c2 - c1);
    if tau >= 1.0 {
        (1.0, 0.0)
    } else if tau <= 0.0 {
        (0.0, 0.0)
    } else {
        let w = tau * tau * tau * (10.0 - 15.0 * tau + 6.0 * tau * tau);
        let dw = -30.0 * tau * tau * (1.0 - tau) * (1.0 - tau) / (c2 - c1);
        (w, dw)
    }
}

/// The smooth positive part of [`SecondNeighbourRule::OneSided`] and its derivative:
/// `ψ(z) = 0` (`z ≤ 0`), `6z³ − 8z⁴ + 3z⁵` (`0 < z < 1`), `z` (`z ≥ 1`). `ψ` and `ψ'` are
/// continuous, `ψ'' ` too (`ψ''(0) = ψ''(1) = 0`); `ψ' = z²(18 − 32z + 15z²) > 0` inside.
///
/// [formula: rc-near-onesided]
pub fn smooth_positive_part(z: f64) -> (f64, f64) {
    if z <= 0.0 {
        (0.0, 0.0)
    } else if z >= 1.0 {
        (z, 1.0)
    } else {
        let z2 = z * z;
        (
            z2 * z * (6.0 - 8.0 * z + 3.0 * z2),
            z2 * (18.0 - 32.0 * z + 15.0 * z2),
        )
    }
}

/// The gate factor `G(x) = W(min(x, 1))` of [`SecondNeighbourRule::Gated`] (`W` the quintic of
/// [`blend_weight`]) and `dG/dx`, for `x ≥ 0`.
///
/// [formula: rc-near-gate]
pub fn gate_factor(x: f64) -> (f64, f64) {
    if x <= 0.0 {
        (0.0, 0.0)
    } else if x >= 1.0 {
        (1.0, 0.0)
    } else {
        (
            x * x * x * (10.0 - 15.0 * x + 6.0 * x * x),
            30.0 * x * x * (1.0 - x) * (1.0 - x),
        )
    }
}

// =============================================================================================
// Inner integral
// =============================================================================================

/// The exponent-dependent tables of the inner integral.
#[derive(Debug, Clone)]
pub struct InnerTables {
    alpha: f64,
    /// `m = β − 2`.
    m: f64,
    /// `γ = 1 + α − β`.
    gamma: f64,
    /// `C(a, k)(−1)^k/(2k + 1)`, `a = (m − 1)/2`.
    series_a: [f64; SERIES_TERMS],
    /// `C(−½, k)(−1)^k/(m + 1 + 2k)`.
    series_r: [f64; SERIES_TERMS],
    /// `S(π/2)`.
    s_half: f64,
}

impl InnerTables {
    /// The tables for `exponents`.
    pub fn new(exponents: Exponents) -> InnerTables {
        let (alpha, beta) = (exponents.alpha(), exponents.beta());
        let m = beta - 2.0;
        let a = 0.5 * (m - 1.0);
        let mut series_a = [0.0; SERIES_TERMS];
        let mut series_r = [0.0; SERIES_TERMS];
        let (mut ca, mut cr) = (1.0f64, 1.0f64);
        for k in 0..SERIES_TERMS {
            let kf = k as f64;
            series_a[k] = ca / (2.0 * kf + 1.0);
            series_r[k] = cr / (m + 1.0 + 2.0 * kf);
            ca *= -(a - kf) / (kf + 1.0);
            cr *= (kf + 0.5) / (kf + 1.0);
        }
        let mut t = InnerTables {
            alpha,
            m,
            gamma: 1.0 + alpha - beta,
            series_a,
            series_r,
            s_half: 0.0,
        };
        let c = std::f64::consts::FRAC_1_SQRT_2;
        t.s_half = t.series_low(c) + t.series_high(c);
        t
    }

    /// `S(θ)` for `|θ| ≤ π/4` from `s = sin θ`: `Σ_k C(a,k)(−1)^k s^{2k+1}/(2k+1)`.
    #[inline]
    fn series_low(&self, s: f64) -> f64 {
        let s2 = s * s;
        let mut acc = 0.0;
        for k in (0..SERIES_TERMS).rev() {
            acc = acc * s2 + self.series_a[k];
        }
        acc * s
    }

    /// `R(c) = ∫_0^ψ sin^m φ dφ` with `c = sin ψ`, `ψ ≤ π/4`:
    /// `c^{m+1} Σ_k C(−½,k)(−1)^k c^{2k}/(m + 1 + 2k)`.
    #[inline]
    fn series_high(&self, c: f64) -> f64 {
        let c2 = c * c;
        let mut acc = 0.0;
        for k in (0..SERIES_TERMS).rev() {
            acc = acc * c2 + self.series_r[k];
        }
        acc * c.powf(self.m + 1.0)
    }

    /// `S(θ) = ∫_0^θ cos^m φ dφ` for `θ = atan2(sn, cs)`, `cs ≥ 0`, `sn² + cs² = 1`.
    ///
    /// [formula: rc-near-cos-series]
    #[inline]
    pub fn cos_power_integral(&self, sn: f64, cs: f64) -> f64 {
        if sn.abs() <= cs {
            self.series_low(sn)
        } else if sn > 0.0 {
            self.s_half - self.series_high(cs)
        } else {
            self.series_high(cs) - self.s_half
        }
    }

    /// `S(θ₁) − S(θ₀)` for `θ_i = atan2(X_i, h)`, `X₁ > X₀`; when both angles are beyond `π/4`
    /// on the same side, the difference of the two `R` terms (no cancellation of `S(π/2)`).
    #[inline]
    fn angle_integral(&self, x1: f64, r1: f64, x0: f64, r0: f64, h: f64) -> f64 {
        if x0 >= h {
            self.series_high(h / r0) - self.series_high(h / r1)
        } else if x1 <= -h {
            self.series_high(h / r1) - self.series_high(h / r0)
        } else {
            self.cos_power_integral(x1 / r1, h / r1) - self.cos_power_integral(x0 / r0, h / r0)
        }
    }

    /// `F(y) = ∫_I k(x, y, T_I) dx` for the edge `(a, a1)`, and its gradients with respect to
    /// `y`, `a`, `a1` (`SPEC.md` §3.4). `h = |P_T (y − a)|` is the value of `|T × (x − y)|` at
    /// every `x` of the edge.
    ///
    /// [formula: rc-near-h-constant]
    /// [formula: rc-near-inner]
    #[inline]
    pub fn inner(&self, y: P, a: P, a1: P) -> InnerEval {
        let e = sub(a1, a);
        let l = norm(e);
        let t = scale(e, 1.0 / l);
        let r = sub(y, a);
        let u0 = dot(t, r);
        let n = axpy(r, -u0, t);
        let h = norm(n);
        if h == 0.0 {
            // On the line of I, outside I (inside would be an intersection): F = 0, and its
            // gradient vanishes there because α > 1.
            return InnerEval {
                value: 0.0,
                grad_y: [0.0; 3],
                grad_a: [0.0; 3],
                grad_a1: [0.0; 3],
            };
        }
        let (x1, x0) = (l - u0, -u0);
        let (r1, r0) = (x1.hypot(h), x0.hypot(h));
        let c = self.angle_integral(x1, r1, x0, r0, h);
        let hg = h.powf(self.gamma);
        let value = hg * c;
        // ∂C/∂X_i = ±cos^m θ_i · h/R_i², ∂C/∂h = −cos^m θ₁ X₁/R₁² + cos^m θ₀ X₀/R₀².
        let (cm1, cm0) = ((h / r1).powf(self.m), (h / r0).powf(self.m));
        let (q1, q0) = (cm1 / (r1 * r1), cm0 / (r0 * r0));
        let dc_dx1 = q1 * h;
        let dc_dx0 = -q0 * h;
        let dc_dh = -q1 * x1 + q0 * x0;
        let f_u = hg * (-dc_dx1 - dc_dx0);
        let f_l = hg * dc_dx1;
        let f_h = self.gamma * value / h + hg * dc_dh;
        let nh = scale(n, 1.0 / h);
        let grad_y = axpy(scale(t, f_u), f_h, nh);
        // δu₀ = n·(δa1 − δa)/ℓ + T·(δy − δa),  δh = n̂·(δy − δa) − u₀ n̂·(δa1 − δa)/ℓ,
        // δℓ = T·(δa1 − δa).
        let grad_a1 = add(
            add(scale(n, f_u / l), scale(nh, -f_h * u0 / l)),
            scale(t, f_l),
        );
        let grad_a = scale(add(grad_y, grad_a1), -1.0);
        InnerEval {
            value,
            grad_y,
            grad_a,
            grad_a1,
        }
    }

    /// `F(y)` only.
    #[inline]
    pub fn inner_value(&self, y: P, a: P, a1: P) -> f64 {
        let e = sub(a1, a);
        let l = norm(e);
        let t = scale(e, 1.0 / l);
        let r = sub(y, a);
        let u0 = dot(t, r);
        let h = norm(axpy(r, -u0, t));
        if h == 0.0 {
            return 0.0;
        }
        let (x1, x0) = (l - u0, -u0);
        let c = self.angle_integral(x1, x1.hypot(h), x0, x0.hypot(h), h);
        h.powf(self.gamma) * c
    }

    /// `α`.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }
}

/// Value and gradients of the inner integral [`InnerTables::inner`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InnerEval {
    /// `F(y)`.
    pub value: f64,
    /// `∂F/∂y`.
    pub grad_y: [f64; 3],
    /// `∂F/∂γ_{i₁}`.
    pub grad_a: [f64; 3],
    /// `∂F/∂γ_{i₂}`.
    pub grad_a1: [f64; 3],
}

// =============================================================================================
// Closest points of two segments
// =============================================================================================

/// Closest points `p₀ + s (p₁ − p₀)`, `q₀ + t (q₁ − q₀)` of two segments, their distance and the
/// derivatives of `d`, `s`, `t` with respect to `(p₀, p₁, q₀, q₁)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Closest {
    /// `d = |P(s) − Q(t)|`.
    pub d: f64,
    /// Parameter on the first segment.
    pub s: f64,
    /// Parameter on the second segment.
    pub t: f64,
    /// `∇d` (envelope theorem).
    pub grad_d: Grad4,
    /// `∇s` (0 where `s` is clamped).
    pub grad_s: Grad4,
    /// `∇t` (0 where `t` is clamped).
    pub grad_t: Grad4,
}

/// The closest points of the segments `p₀p₁`, `q₀q₁` (Ericson 2004, §5.1.9; the parallel case
/// takes `s = 0`), with derivatives: `∇d = n̂ ⊗ ((1−s), s, −(1−t), −t)` by the envelope theorem;
/// `∇s, ∇t` by implicit differentiation of `A·D = 0` (if `s` is free) and `B·D = 0` (if `t` is
/// free), `D = p₀ − q₀ + sA − tB`, `A = p₁ − p₀`, `B = q₁ − q₀`. The caller guarantees `d > 0`.
///
/// [formula: rc-near-closest]
pub fn closest_points(p0: P, p1: P, q0: P, q1: P) -> Closest {
    let aa_v = sub(p1, p0);
    let bb_v = sub(q1, q0);
    let r = sub(p0, q0);
    let aa = dot(aa_v, aa_v);
    let ee = dot(bb_v, bb_v);
    let f = dot(bb_v, r);
    let c = dot(aa_v, r);
    let bb = dot(aa_v, bb_v);
    let denom = aa * ee - bb * bb;
    // `s_solved`: `s` satisfies `A·D = 0` (it is not clamped and not the parallel choice).
    let (mut s, mut s_solved) = if denom > 1e-14 * aa * ee {
        let s = (bb * f - c * ee) / denom;
        (s.clamp(0.0, 1.0), s > 0.0 && s < 1.0)
    } else {
        (0.0, false)
    };
    let mut t = (bb * s + f) / ee;
    let mut t_free = t > 0.0 && t < 1.0;
    if t <= 0.0 || t >= 1.0 {
        t = t.clamp(0.0, 1.0);
        t_free = false;
        let s1 = (bb * t - c) / aa;
        s = s1.clamp(0.0, 1.0);
        s_solved = s1 > 0.0 && s1 < 1.0;
    }
    let dv = sub(axpy(r, s, aa_v), scale(bb_v, t));
    let d = norm(dv);
    let nh = scale(dv, 1.0 / d);
    let grad_d = [
        scale(nh, 1.0 - s),
        scale(nh, s),
        scale(nh, -(1.0 - t)),
        scale(nh, -t),
    ];
    let s_free = s_solved;
    // r1 = −∂(A·D)/∂(p₀, p₁, q₀, q₁) at fixed (s, t); r2 = −∂(B·D)/∂(…).
    let r1: Grad4 = [
        scale(sub(scale(aa_v, 1.0 - s), dv), -1.0),
        scale(axpy(dv, s, aa_v), -1.0),
        scale(aa_v, 1.0 - t),
        scale(aa_v, t),
    ];
    let r2: Grad4 = [
        scale(bb_v, -(1.0 - s)),
        scale(bb_v, -s),
        add(scale(bb_v, 1.0 - t), dv),
        sub(scale(bb_v, t), dv),
    ];
    let mut grad_s = ZERO4;
    let mut grad_t = ZERO4;
    if s_free && t_free {
        // [aa −bb; bb −ee] [δs; δt] = [r1; r2].
        for k in 0..4 {
            grad_s[k] = scale(sub(scale(r1[k], ee), scale(r2[k], bb)), 1.0 / denom);
            grad_t[k] = scale(sub(scale(r1[k], bb), scale(r2[k], aa)), 1.0 / denom);
        }
    } else if t_free {
        for k in 0..4 {
            grad_t[k] = scale(r2[k], -1.0 / ee);
        }
    } else if s_free {
        for k in 0..4 {
            grad_s[k] = scale(r1[k], 1.0 / aa);
        }
    }
    Closest {
        d,
        s,
        t,
        grad_d,
        grad_s,
        grad_t,
    }
}

/// The parameter of the projection of `x` onto the segment `q₀q₁`, clamped to `[0, 1]`, and its
/// gradient with respect to `(x, q₀, q₁)`.
fn projection_parameter(x: P, q0: P, q1: P) -> (f64, [[f64; 3]; 3]) {
    let b = sub(q1, q0);
    let ee = dot(b, b);
    let w = sub(x, q0);
    let q = dot(w, b) / ee;
    if q <= 0.0 {
        return (0.0, [[0.0; 3]; 3]);
    }
    if q >= 1.0 {
        return (1.0, [[0.0; 3]; 3]);
    }
    let gx = scale(b, 1.0 / ee);
    let gq1 = scale(axpy(w, -2.0 * q, b), 1.0 / ee);
    let gq0 = scale(add(gx, gq1), -1.0);
    (q, [gx, gq0, gq1])
}

// =============================================================================================
// Outer integral: one ordered pair
// =============================================================================================

/// A breakpoint of the outer integral with its gradient.
#[derive(Clone, Copy)]
struct Breakpoint {
    z: f64,
    grad: Grad4,
}

fn grad4_axpy(acc: &mut Grad4, s: f64, g: &Grad4) {
    for k in 0..4 {
        acc[k] = axpy(acc[k], s, g[k]);
    }
}

/// `e^near_IJ = ∫_I ∫_J k(x, y, T_I) dx dy` for the edge `I = p₀p₁` and `J = q₀q₁` by the graded
/// quadrature of the module docs, with its gradient with respect to `(p₀, p₁, q₀, q₁)` if
/// `with_gradient`. `t_star` is the closest-point parameter on `J` with its gradient and
/// `d` the segment distance with its gradient (both in the order `(p₀, p₁, q₀, q₁)`).
///
/// [formula: rc-near-outer]
/// [formula: rc-near-quadrature-error]
#[allow(clippy::too_many_arguments)]
pub fn ordered_near_energy(
    tables: &InnerTables,
    p: [P; 2],
    q: [P; 2],
    t_star: (f64, &Grad4),
    d: (f64, &Grad4),
    subpanels: usize,
    with_gradient: bool,
) -> (f64, Grad4) {
    let bv = sub(q[1], q[0]);
    let lj = norm(bv);
    let tj = scale(bv, 1.0 / lj);
    let delta = d.0 / lj;
    // ∇δ = ∇d/ℓ_J − d/ℓ_J² ∇ℓ_J, ∇ℓ_J = (0, 0, −T_J, T_J).
    let mut grad_delta = ZERO4;
    if with_gradient {
        grad4_axpy(&mut grad_delta, 1.0 / lj, d.1);
        let s = d.0 / (lj * lj);
        grad_delta[2] = axpy(grad_delta[2], s, tj);
        grad_delta[3] = axpy(grad_delta[3], -s, tj);
    }
    // Breakpoints: t*, projections of p₀ and p₁ (gradients re-ordered to (p₀, p₁, q₀, q₁)).
    let mut bps = [Breakpoint {
        z: t_star.0,
        grad: *t_star.1,
    }; 3];
    for (k, &x) in p.iter().enumerate() {
        let (z, g) = projection_parameter(x, q[0], q[1]);
        let mut grad = ZERO4;
        grad[k] = g[0];
        grad[2] = g[1];
        grad[3] = g[2];
        bps[k + 1] = Breakpoint { z, grad };
    }
    // Stable sort of three (ties keep the order t*, p₀, p₁).
    let mut order = [0usize, 1, 2];
    for i in 1..3 {
        let mut j = i;
        while j > 0 && bps[order[j - 1]].z > bps[order[j]].z {
            order.swap(j - 1, j);
            j -= 1;
        }
    }
    let zs = [0.0, bps[order[0]].z, bps[order[1]].z, bps[order[2]].z, 1.0];
    // Adjoints of the five panel ends (the first and last are constants) and of δ.
    let mut adj_z = [0.0f64; 5];
    let mut adj_delta = 0.0f64;
    let mut sum = 0.0f64;
    let mut direct = ZERO4;
    let m = subpanels as f64;
    // A panel end moves with the curve unless it is 0, 1 or a clamped breakpoint (zero gradient).
    let moving = |slot: usize| {
        (1..4).contains(&slot)
            && bps[order[slot - 1]]
                .grad
                .iter()
                .flatten()
                .any(|&g| g != 0.0)
    };
    for panel in 0..4 {
        let (lo, hi) = (zs[panel], zs[panel + 1]);
        let half = 0.5 * (hi - lo);
        // An empty panel adds nothing to the value, but its derivative with respect to a moving
        // end is not zero (∂/∂hi ∫_lo^hi f = f(hi)): keep it unless both ends are fixed.
        if half <= 0.0 && (!with_gradient || !(moving(panel) || moving(panel + 1))) {
            continue;
        }
        let q_len = half.hypot(delta);
        let w_max = (half / delta).asinh();
        let dw_dl = 1.0 / q_len;
        let dw_dd = -half / (delta * q_len);
        for (end, sigma) in [(panel, 1.0f64), (panel + 1, -1.0f64)] {
            let c = zs[end];
            let (mut adj_c, mut adj_l, mut adj_d) = (0.0, 0.0, 0.0);
            for j in 0..subpanels {
                for &(xi, wk) in &GL8 {
                    let eta = (j as f64 + 0.5 * (1.0 + xi)) / m;
                    let om = wk / (2.0 * m);
                    let we = w_max * eta;
                    let (sh, ch) = (we.sinh(), we.cosh());
                    let t = c + sigma * delta * sh;
                    let wt = om * delta * w_max * ch;
                    let y = axpy(q[0], t, bv);
                    if !with_gradient {
                        sum += wt * tables.inner_value(y, p[0], p[1]);
                        continue;
                    }
                    let fe = tables.inner(y, p[0], p[1]);
                    sum += wt * fe.value;
                    // Direct dependence on the four endpoints at fixed (t, wt).
                    direct[0] = axpy(direct[0], wt, fe.grad_a);
                    direct[1] = axpy(direct[1], wt, fe.grad_a1);
                    direct[2] = axpy(direct[2], wt * (1.0 - t), fe.grad_y);
                    direct[3] = axpy(direct[3], wt * t, fe.grad_y);
                    // Through the node t and the weight wt.
                    let g_t = wt * dot(fe.grad_y, bv);
                    let g_w = fe.value;
                    let dt_dd = sigma * (sh + delta * ch * eta * dw_dd);
                    let dt_dl = sigma * delta * ch * eta * dw_dl;
                    let dwt_dd = om * (w_max * ch + delta * dw_dd * (ch + we * sh));
                    let dwt_dl = om * delta * dw_dl * (ch + we * sh);
                    adj_c += g_t;
                    adj_d += g_t * dt_dd + g_w * dwt_dd;
                    adj_l += g_t * dt_dl + g_w * dwt_dl;
                }
            }
            if with_gradient {
                adj_z[end] += adj_c;
                adj_z[panel + 1] += 0.5 * adj_l;
                adj_z[panel] -= 0.5 * adj_l;
                adj_delta += adj_d;
            }
        }
    }
    let value = lj * sum;
    if !with_gradient {
        return (value, ZERO4);
    }
    // ∇E = ℓ_J (direct + Σ adj_z ∇z + adj_δ ∇δ) + Q ∇ℓ_J.
    let mut g = direct;
    for (k, &o) in order.iter().enumerate() {
        grad4_axpy(&mut g, adj_z[k + 1], &bps[o].grad);
    }
    grad4_axpy(&mut g, adj_delta, &grad_delta);
    for gk in g.iter_mut() {
        *gk = scale(*gk, lj);
    }
    g[2] = axpy(g[2], -sum, tj);
    g[3] = axpy(g[3], sum, tj);
    (value, g)
}

/// The paper's ordered pair term `e_IJ = k̂_IJ ℓ_I ℓ_J` (Eq. 17–18) with its gradient with
/// respect to `(p₀, p₁, q₀, q₁)` (`SPEC.md` §3.2).
///
/// [formula: rc-energy-differential]
pub fn paper_pair(p: [P; 2], q: [P; 2], e: Exponents) -> (f64, Grad4) {
    let ev = sub(p[1], p[0]);
    let li = norm(ev);
    let ti = scale(ev, 1.0 / li);
    let fv = sub(q[1], q[0]);
    let lj = norm(fv);
    let tj = scale(fv, 1.0 / lj);
    let mut khat = 0.0;
    let mut gt = [0.0; 3];
    let mut gp = [[0.0; 3]; 2];
    let mut gq = [[0.0; 3]; 2];
    for a in 0..2 {
        for b in 0..2 {
            let ke = kernel_with_gradient(p[a], q[b], ti, e);
            khat += ke.value;
            gt = add(gt, ke.grad_t);
            gp[a] = add(gp[a], ke.grad_p);
            gq[b] = sub(gq[b], ke.grad_p);
        }
    }
    khat *= 0.25;
    let ll = 0.25 * li * lj;
    let first = axpy(scale(ti, khat), 0.25, gt);
    let signs = [-1.0, 1.0];
    let mut g = ZERO4;
    for a in 0..2 {
        g[a] = axpy(scale(first, signs[a] * lj), ll, gp[a]);
        g[2 + a] = axpy(scale(tj, signs[a] * li * khat), ll, gq[a]);
    }
    (khat * li * lj, g)
}

/// The paper's ordered pair term only.
fn paper_pair_value(p: [P; 2], q: [P; 2], e: Exponents) -> f64 {
    let ev = sub(p[1], p[0]);
    let li = norm(ev);
    let ti = scale(ev, 1.0 / li);
    let lj = norm(sub(q[1], q[0]));
    let mut khat = 0.0;
    for a in p {
        for b in q {
            khat += crate::kernel::kernel(a, b, ti, e);
        }
    }
    0.25 * khat * li * lj
}

// =============================================================================================
// One unordered pair: the correction
// =============================================================================================

/// The correction of one unordered edge pair `{I, J}` and its gradient with respect to
/// `(p₀, p₁, q₀, q₁)` (`I = p₀p₁`, `J = q₀q₁`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairCorrection {
    /// `w(r) Σ (e^near − e^paper)` over the included orders.
    pub value: f64,
    /// Its gradient.
    pub grad: Grad4,
    /// `w(r_IJ)`.
    pub weight: f64,
    /// `r_IJ = d_IJ/ℓ_IJ`.
    pub ratio: f64,
}

/// The near-pair correction `w(r_IJ) Σ_o (e^near_o − e^paper_o)` of the unordered pair
/// `{I, J}` over the included orders `o ∈ {(I, J), (J, I)}` (`include_ij`, `include_ji`), with
/// its gradient if `with_gradient` (`SPEC.md` §3.4). Zero (and no work) for `r ≥ c₂`.
///
/// [formula: rc-near-differential]
#[allow(clippy::too_many_arguments)]
pub fn pair_correction(
    tables: &InnerTables,
    exponents: Exponents,
    options: &NearPairOptions,
    p: [P; 2],
    q: [P; 2],
    include_ij: bool,
    include_ji: bool,
    with_gradient: bool,
) -> PairCorrection {
    pair_correction_of_kind(
        tables,
        exponents,
        options,
        p,
        q,
        None,
        include_ij,
        include_ji,
        with_gradient,
    )
}

/// [`pair_correction`] for a pair of the given kind: `joined = Some((a, b))` for second
/// neighbours, whose joining edge `K` meets `I` at its end `a ∈ {0, 1}` (`p_a`) and `J` at its
/// end `b` (`q_b`); the pair then follows [`NearPairOptions::second_neighbour_rule`] and
/// [`NearPairOptions::second_neighbour_band`] (`SPEC.md` §17p). `None`: every other pair (§3.4).
///
/// [formula: rc-near-differential]
#[allow(clippy::too_many_arguments)]
pub fn pair_correction_of_kind(
    tables: &InnerTables,
    exponents: Exponents,
    options: &NearPairOptions,
    p: [P; 2],
    q: [P; 2],
    joined: Option<(usize, usize)>,
    include_ij: bool,
    include_ji: bool,
    with_gradient: bool,
) -> PairCorrection {
    let cl = closest_points(p[0], p[1], q[0], q[1]);
    let ev = sub(p[1], p[0]);
    let fv = sub(q[1], q[0]);
    let (li, lj) = (norm(ev), norm(fv));
    let lm = (0.5 * (li * li + lj * lj)).sqrt();
    let ratio = cl.d / lm;
    let rule = match joined {
        Some(_) => options.second_neighbour_rule,
        None => SecondNeighbourRule::TwoSided,
    };
    let (c1, c2) = match joined {
        Some(_) => options
            .second_neighbour_band
            .unwrap_or((options.blend_start, options.blend_end)),
        None => (options.blend_start, options.blend_end),
    };
    let (w, dw) = blend_weight(ratio, c1, c2);
    if w == 0.0 || !(include_ij || include_ji) {
        return PairCorrection {
            value: 0.0,
            grad: ZERO4,
            weight: w,
            ratio,
        };
    }
    let mut diff = 0.0;
    let mut gdiff = ZERO4;
    // The sums E^n, E^p over the orders and their gradients (the one-sided rules only).
    let (mut en_sum, mut ep_sum) = (0.0f64, 0.0f64);
    let (mut gen_sum, mut gep_sum) = (ZERO4, ZERO4);
    if include_ij {
        let (en, gn) = ordered_near_energy(
            tables,
            p,
            q,
            (cl.t, &cl.grad_t),
            (cl.d, &cl.grad_d),
            options.subpanels,
            with_gradient,
        );
        if with_gradient {
            let (ep, gp) = paper_pair(p, q, exponents);
            diff += en - ep;
            grad4_axpy(&mut gdiff, 1.0, &gn);
            grad4_axpy(&mut gdiff, -1.0, &gp);
            en_sum += en;
            ep_sum += ep;
            grad4_axpy(&mut gen_sum, 1.0, &gn);
            grad4_axpy(&mut gep_sum, 1.0, &gp);
        } else {
            let ep = paper_pair_value(p, q, exponents);
            diff += en - ep;
            en_sum += en;
            ep_sum += ep;
        }
    }
    if include_ji {
        // (J, I): the outer integral runs over I; re-order the gradients to (q₀, q₁, p₀, p₁).
        let sw = |g: &Grad4| [g[2], g[3], g[0], g[1]];
        let (gs, gd) = (sw(&cl.grad_s), sw(&cl.grad_d));
        let (en, gn) = ordered_near_energy(
            tables,
            q,
            p,
            (cl.s, &gs),
            (cl.d, &gd),
            options.subpanels,
            with_gradient,
        );
        if with_gradient {
            let (ep, gp) = paper_pair(q, p, exponents);
            diff += en - ep;
            grad4_axpy(&mut gdiff, 1.0, &sw(&gn));
            grad4_axpy(&mut gdiff, -1.0, &sw(&gp));
            en_sum += en;
            ep_sum += ep;
            grad4_axpy(&mut gen_sum, 1.0, &sw(&gn));
            grad4_axpy(&mut gep_sum, 1.0, &sw(&gp));
        } else {
            let ep = paper_pair_value(q, p, exponents);
            diff += en - ep;
            en_sum += en;
            ep_sum += ep;
        }
    }
    if rule != SecondNeighbourRule::TwoSided {
        return second_neighbour_term(
            rule,
            joined.unwrap_or((1, 0)),
            &cl,
            (w, dw),
            ratio,
            (lm, li, lj, ev, fv),
            (diff, &gdiff),
            (en_sum, &gen_sum),
            (ep_sum, &gep_sum),
            with_gradient,
        );
    }
    let mut grad = ZERO4;
    if with_gradient {
        grad4_axpy(&mut grad, w, &gdiff);
        if dw != 0.0 {
            // ∇r = ∇d/ℓ_m − r/ℓ_m ∇ℓ_m,  ∇ℓ_m = (ℓ_I ∇ℓ_I + ℓ_J ∇ℓ_J)/(2 ℓ_m).
            let mut gr = ZERO4;
            grad4_axpy(&mut gr, 1.0 / lm, &cl.grad_d);
            let ci = -ratio / lm * li / (2.0 * lm);
            let cj = -ratio / lm * lj / (2.0 * lm);
            let (ti, tj) = (scale(ev, 1.0 / li), scale(fv, 1.0 / lj));
            gr[0] = axpy(gr[0], -ci, ti);
            gr[1] = axpy(gr[1], ci, ti);
            gr[2] = axpy(gr[2], -cj, tj);
            gr[3] = axpy(gr[3], cj, tj);
            grad4_axpy(&mut grad, diff * dw, &gr);
        }
    }
    PairCorrection {
        value: w * diff,
        grad,
        weight: w,
        ratio,
    }
}

/// `∇r` of `r = d/ℓ_IJ` (`SPEC.md` §3.4).
fn ratio_gradient(cl: &Closest, ratio: f64, lm: f64, li: f64, lj: f64, ev: P, fv: P) -> Grad4 {
    let mut gr = ZERO4;
    grad4_axpy(&mut gr, 1.0 / lm, &cl.grad_d);
    let ci = -ratio / lm * li / (2.0 * lm);
    let cj = -ratio / lm * lj / (2.0 * lm);
    let (ti, tj) = (scale(ev, 1.0 / li), scale(fv, 1.0 / lj));
    gr[0] = axpy(gr[0], -ci, ti);
    gr[1] = axpy(gr[1], ci, ti);
    gr[2] = axpy(gr[2], -cj, tj);
    gr[3] = axpy(gr[3], cj, tj);
    gr
}

/// The second-neighbour rules other than [`SecondNeighbourRule::TwoSided`] (`SPEC.md` §17p).
#[allow(clippy::too_many_arguments)]
fn second_neighbour_term(
    rule: SecondNeighbourRule,
    (a, b): (usize, usize),
    cl: &Closest,
    (w, dw): (f64, f64),
    ratio: f64,
    (lm, li, lj, ev, fv): (f64, f64, f64, P, P),
    (diff, gdiff): (f64, &Grad4),
    (en, gen_): (f64, &Grad4),
    (ep, gep): (f64, &Grad4),
    with_gradient: bool,
) -> PairCorrection {
    // The gate g = G(|s* − s_K|/σ) G(|t* − t_K|/σ) and its gradient.
    let (g, grad_g) = match rule {
        SecondNeighbourRule::Gated { width }
        | SecondNeighbourRule::OneSidedGated { width }
        | SecondNeighbourRule::OneSidedGatedEither { width } => {
            let either = matches!(rule, SecondNeighbourRule::OneSidedGatedEither { .. });
            let (sk, tk) = (a as f64, b as f64);
            let (xs, xt) = ((cl.s - sk).abs() / width, (cl.t - tk).abs() / width);
            let (gs, dgs) = gate_factor(xs);
            let (gt, dgt) = gate_factor(xt);
            let mut gg = ZERO4;
            if with_gradient {
                // d|s − s_K|/ds = +1 for s_K = 0, −1 for s_K = 1.
                let ss = if a == 0 { 1.0 } else { -1.0 } / width;
                let st = if b == 0 { 1.0 } else { -1.0 } / width;
                // Product gate: ∇(G_s G_t); either gate: ∇(1 − (1 − G_s)(1 − G_t)).
                let (fs, ft) = if either {
                    (1.0 - gt, 1.0 - gs)
                } else {
                    (gt, gs)
                };
                grad4_axpy(&mut gg, dgs * ss * fs, &cl.grad_s);
                grad4_axpy(&mut gg, ft * dgt * st, &cl.grad_t);
            }
            if either {
                (1.0 - (1.0 - gs) * (1.0 - gt), gg)
            } else {
                (gs * gt, gg)
            }
        }
        _ => (1.0, ZERO4),
    };
    let one_sided = matches!(
        rule,
        SecondNeighbourRule::OneSided
            | SecondNeighbourRule::OneSidedGated { .. }
            | SecondNeighbourRule::OneSidedGatedEither { .. }
            | SecondNeighbourRule::OneSidedFoldGated { .. }
    );
    // The un-weighted term c (w and g factored out) and its gradient.
    let (c, grad_c) = if one_sided {
        if ep.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            (0.0, ZERO4)
        } else {
            let rho = en / ep;
            let (psi, dpsi) = smooth_positive_part(rho - 1.0);
            let mut gc = ZERO4;
            if with_gradient && psi != 0.0 {
                // ∇(E^p ψ(E^n/E^p − 1)) = ψ ∇E^p + ψ' (∇E^n − ρ ∇E^p).
                grad4_axpy(&mut gc, psi - dpsi * rho, gep);
                grad4_axpy(&mut gc, dpsi, gen_);
            }
            (ep * psi, gc)
        }
    } else {
        (diff, *gdiff)
    };
    if let SecondNeighbourRule::OneSidedFoldGated {
        contact_band: (cc1, cc2),
        cos_shut,
        cos_open,
    } = rule
    {
        return fold_gated_term(
            (a, b),
            cl,
            (w, dw),
            (cc1, cc2, cos_shut, cos_open),
            ratio,
            (lm, li, lj, ev, fv),
            (c, &grad_c),
            with_gradient,
        );
    }
    let mut grad = ZERO4;
    if with_gradient && c != 0.0 {
        grad4_axpy(&mut grad, w * g, &grad_c);
        grad4_axpy(&mut grad, w * c, &grad_g);
        if dw != 0.0 {
            let gr = ratio_gradient(cl, ratio, lm, li, lj, ev, fv);
            grad4_axpy(&mut grad, g * c * dw, &gr);
        }
    } else if with_gradient && !one_sided {
        // c = 0 with a nonzero gradient is possible for the gated two-sided rule.
        grad4_axpy(&mut grad, w * g, &grad_c);
    }
    PairCorrection {
        value: w * g * c,
        grad,
        weight: w,
        ratio,
    }
}

/// [`SecondNeighbourRule::OneSidedFoldGated`]: `W c` with `W = 1 − (1 − w_c)(1 − w G)` and
/// `c = E^p ψ(E^n/E^p − 1)` (`SPEC.md` §17r);
/// `∇W = (1 − w G) w_c' ∇r + (1 − w_c)(G w' ∇r + w G' ∇cos φ)`,
/// `∇_u cos φ = (v̂ − cos φ û)/|u|` for `u = ±(p₁ − p₀)` (toward `K`), `v = ±(q₁ − q₀)` (away).
///
/// [formula: rc-near-foldgate]
#[allow(clippy::too_many_arguments)]
fn fold_gated_term(
    (a, b): (usize, usize),
    cl: &Closest,
    (w, dw): (f64, f64),
    (cc1, cc2, cos_shut, cos_open): (f64, f64, f64, f64),
    ratio: f64,
    (lm, li, lj, ev, fv): (f64, f64, f64, P, P),
    (c, grad_c): (f64, &Grad4),
    with_gradient: bool,
) -> PairCorrection {
    let (wc, dwc) = blend_weight(ratio, cc1, cc2);
    // I oriented toward K (K meets I at p_a), J away from K (K meets J at q_b).
    let sa = if a == 1 { 1.0 } else { -1.0 };
    let sb = if b == 0 { 1.0 } else { -1.0 };
    let (uh, vh) = (scale(ev, sa / li), scale(fv, sb / lj));
    let cosphi = dot(uh, vh);
    let (gf, dgf) = gate_factor((cosphi - cos_shut) / (cos_open - cos_shut));
    let dgf = dgf / (cos_open - cos_shut);
    let wg = w * gf;
    let weight = 1.0 - (1.0 - wc) * (1.0 - wg);
    let mut grad = ZERO4;
    if with_gradient && c != 0.0 {
        grad4_axpy(&mut grad, weight, grad_c);
        let gr = ratio_gradient(cl, ratio, lm, li, lj, ev, fv);
        let coef_r = (1.0 - wg) * dwc + (1.0 - wc) * gf * dw;
        if coef_r != 0.0 {
            grad4_axpy(&mut grad, c * coef_r, &gr);
        }
        let coef_cos = (1.0 - wc) * w * dgf;
        if coef_cos != 0.0 {
            // ∇cos φ with respect to ev and fv, then to the endpoints.
            let gu = scale(sub(vh, scale(uh, cosphi)), sa / li);
            let gv = scale(sub(uh, scale(vh, cosphi)), sb / lj);
            let s = c * coef_cos;
            grad[0] = axpy(grad[0], -s, gu);
            grad[1] = axpy(grad[1], s, gu);
            grad[2] = axpy(grad[2], -s, gv);
            grad[3] = axpy(grad[3], s, gv);
        }
    }
    PairCorrection {
        value: weight * c,
        grad,
        weight,
        ratio,
    }
}

/// The ends `(a, b)` of edges `i`, `j` joined by one edge `K` (`p_a` of `i` and `q_b` of `j` are
/// the two ends of `K`), the first such edge in the incidence lists; `None` if none.
pub fn joining_ends(network: &CurveNetwork, i: usize, j: usize) -> Option<(usize, usize)> {
    let e = network.edges();
    for (a, &va) in e[i].iter().enumerate() {
        for &k in network.incident_edges(va) {
            let [x, y] = e[k as usize];
            let other = if x == va { y } else { x };
            if let Some(b) = e[j].iter().position(|&v| v == other) {
                return Some((a, b));
            }
        }
    }
    None
}

// =============================================================================================
// Broad phase and the whole correction
// =============================================================================================

/// `true` if some edge joins an endpoint of edge `i` to an endpoint of edge `j`.
pub fn joined_by_edge(network: &CurveNetwork, i: usize, j: usize) -> bool {
    let e = network.edges();
    let [c, d] = e[j];
    e[i].iter().any(|&a| {
        network.incident_edges(a).iter().any(|&k| {
            let [x, y] = e[k as usize];
            let other = if x == a { y } else { x };
            other == c || other == d
        })
    })
}

/// All unordered pairs `(I, J)`, `I < J`, of disjoint edges with `d_IJ < c₂ ℓ_IJ`, sorted
/// (without pairs joined by one edge unless [`NearPairOptions::second_neighbours`]).
/// Broad phase: a uniform grid of cell size `(1 + 2c₂) max ℓ`, each edge in the cells its box
/// grown by `c₂ max ℓ` touches (sorted cell entries, so no hashing order enters); `O(|E|)` for
/// curves of bounded density.
pub fn near_pairs(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    options: &NearPairOptions,
) -> Vec<(u32, u32)> {
    let edges = network.edges();
    if edges.len() < 2 {
        return Vec::new();
    }
    let ends: Vec<[P; 2]> = edges
        .iter()
        .map(|&[a, b]| [positions[a as usize], positions[b as usize]])
        .collect();
    let lens: Vec<f64> = ends.iter().map(|e| norm(sub(e[1], e[0]))).collect();
    let lmax = lens.iter().copied().fold(0.0f64, f64::max);
    // The second neighbours' band may reach further (`SPEC.md` §17p).
    let c2_second = match (options.second_neighbours, options.second_neighbour_band) {
        (true, Some((_, c2))) => c2,
        _ => options.blend_end,
    };
    let tol = options.blend_end.max(c2_second) * lmax;
    let h = lmax + 2.0 * tol;
    if !(h > 0.0 && h.is_finite()) {
        return Vec::new();
    }
    let cell = |x: f64| (x / h).floor() as i64;
    let mut entries: Vec<([i64; 3], u32)> = Vec::with_capacity(edges.len() * 8);
    for (e, [p, q]) in ends.iter().enumerate() {
        let lo = [0, 1, 2].map(|c| cell(p[c].min(q[c]) - tol));
        let hi = [0, 1, 2].map(|c| cell(p[c].max(q[c]) + tol));
        for x in lo[0]..=hi[0] {
            for y in lo[1]..=hi[1] {
                for z in lo[2]..=hi[2] {
                    entries.push(([x, y, z], e as u32));
                }
            }
        }
    }
    entries.sort_unstable();
    let mut cand: Vec<(u32, u32)> = Vec::new();
    let mut k = 0;
    while k < entries.len() {
        let mut k1 = k + 1;
        while k1 < entries.len() && entries[k1].0 == entries[k].0 {
            k1 += 1;
        }
        for a in k..k1 {
            for b in a + 1..k1 {
                let (i, j) = (entries[a].1, entries[b].1);
                if !network.edges_adjacent(i as usize, j as usize)
                    && (options.second_neighbours
                        || !joined_by_edge(network, i as usize, j as usize))
                {
                    cand.push((i.min(j), i.max(j)));
                }
            }
        }
        k = k1;
    }
    cand.sort_unstable();
    cand.dedup();
    cand.retain(|&(i, j)| {
        let (pi, pj) = (&ends[i as usize], &ends[j as usize]);
        let lm = (0.5 * (lens[i as usize].powi(2) + lens[j as usize].powi(2))).sqrt();
        let c2 =
            if c2_second != options.blend_end && joined_by_edge(network, i as usize, j as usize) {
                c2_second
            } else {
                options.blend_end
            };
        closest_points(pi[0], pi[1], pj[0], pj[1]).d < c2 * lm
    });
    cand
}

/// The result of [`correction`].
#[derive(Debug, Clone, PartialEq)]
pub struct NearCorrection {
    /// The correction of the energy.
    pub energy: f64,
    /// Its differential (one entry per vertex), if asked for.
    pub differential: Option<Vec<[f64; 3]>>,
    /// Totals.
    pub stats: NearPairStats,
}

/// Totals of a near-pair pass.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct NearPairStats {
    /// Unordered pairs with `d_IJ < c₂ ℓ_IJ` (the ones evaluated).
    pub pairs: u64,
    /// Of those, pairs with `w = 1` (`d_IJ ≤ c₁ ℓ_IJ`).
    pub full: u64,
    /// The smallest `r_IJ` (∞ when there is none).
    pub min_ratio: f64,
}

/// The near-pair correction of the whole curve, `Σ_{pairs} w (Σ_o e^near_o − e^paper_o)`, and
/// optionally its differential (one entry per vertex, added in pair order). `include(i, j)`
/// says whether the ordered pair `(i, j)` is evaluated exactly by the caller (always for the
/// exact energy; only the exact entries of the interaction lists for Barnes–Hut). Pairs are
/// evaluated in parallel and summed in their sorted order: the result does not depend on the
/// thread count (D25). Returns `None` when no pair is near (the caller then leaves its energy
/// untouched, bit for bit).
///
/// [formula: rc-near-differential]
#[allow(clippy::too_many_arguments)]
pub fn correction(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
    options: &NearPairOptions,
    include: &(dyn Fn(usize, usize) -> bool + Sync),
    with_gradient: bool,
) -> Option<NearCorrection> {
    let pairs = near_pairs(network, positions, options);
    if pairs.is_empty() {
        return None;
    }
    let tables = InnerTables::new(exponents);
    let edges = network.edges();
    let results: Vec<PairCorrection> = pool.install(|| {
        pairs
            .par_iter()
            .map(|&(i, j)| {
                let [a, b] = edges[i as usize];
                let [c, d] = edges[j as usize];
                let p = [positions[a as usize], positions[b as usize]];
                let q = [positions[c as usize], positions[d as usize]];
                let joined = if options.second_neighbours {
                    joining_ends(network, i as usize, j as usize)
                } else {
                    None
                };
                pair_correction_of_kind(
                    &tables,
                    exponents,
                    options,
                    p,
                    q,
                    joined,
                    include(i as usize, j as usize),
                    include(j as usize, i as usize),
                    with_gradient,
                )
            })
            .collect()
    });
    let mut total = 0.0;
    let mut stats = NearPairStats {
        pairs: pairs.len() as u64,
        full: 0,
        min_ratio: f64::INFINITY,
    };
    let mut diff = if with_gradient {
        Some(vec![[0.0; 3]; network.n_vertices()])
    } else {
        None
    };
    for (&(i, j), r) in pairs.iter().zip(&results) {
        total += r.value;
        if r.weight == 1.0 {
            stats.full += 1;
        }
        stats.min_ratio = stats.min_ratio.min(r.ratio);
        if let Some(g) = diff.as_mut() {
            let [a, b] = edges[i as usize];
            let [c, d] = edges[j as usize];
            for (v, gv) in [a, b, c, d].into_iter().zip(r.grad) {
                g[v as usize] = add(g[v as usize], gv);
            }
        }
    }
    Some(NearCorrection {
        energy: total,
        differential: diff,
        stats,
    })
}
