//! Torus fitting from the predicted coefficients (SPEC.md §7.1; paper Sec. 4.1, Eq. 18–22,
//! App. C; NOTES-REFERENCE.md §6).
//!
//! Ported from `get_torus_parameters_python` and `fit_torus_at_point` of the reference's
//! `torus_distance.cpp` (points-as-tori, MIT, Copyright (c) 2026 Nicole Feng), in f64, with the
//! JAX training rule's `|κ| ≥ 1e-10` floor added (decision D9) and every clamp counted.

use crate::features::local_frame;
use crate::{PatError, PointCloud};
use rayon::prelude::*;

/// A torus. `sign(T) = sign(minor)`: the SDF is `sign⁺(minor)·(‖d‖ − |minor|)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Torus {
    /// Center.
    pub center: [f64; 3],
    /// Unit axis.
    pub axis: [f64; 3],
    /// Major radius `R` (may be smaller than `|minor|`: a spindle torus).
    pub major: f64,
    /// Signed minor radius `s·r`.
    pub minor: f64,
}

impl Torus {
    /// The axis as the reference's spherical angles `(θ, φ) = (atan2(u_y, u_x), acos(u_z))`, for
    /// export only (`u_z` is clamped to `[−1, 1]` first, as the JAX version does).
    pub fn axis_angles(&self) -> [f64; 2] {
        let u = self.axis;
        [u[1].atan2(u[0]), u[2].clamp(-1.0, 1.0).acos()]
    }

    /// The unit axis from the reference's angles: `(cos θ sin φ, sin θ sin φ, cos φ)`.
    pub fn axis_from_angles(theta: f64, phi: f64) -> [f64; 3] {
        let sp = phi.sin();
        [theta.cos() * sp, theta.sin() * sp, phi.cos()]
    }
}

/// The tori of a cloud with the coefficients and scales they came from.
#[derive(Debug, Clone, PartialEq)]
pub struct Tori {
    /// One torus per point.
    pub tori: Vec<Torus>,
    /// The coefficients `(a00, a01, a10, a11, a02, a20)` per point.
    pub coefficients: Vec<[f32; 6]>,
    /// The scales `σ_i` per point.
    pub sigma: Vec<f64>,
}

/// Which clamps fired for one point (one flag per clamp, as the Lean oracle prints them).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ClampFlags {
    /// The discriminant argument was below `eig_eps` (also the definition of "umbilic").
    pub disc: bool,
    /// `ev₊x = v1 + disc` or `ev₋y = v2` was raised to `eig_eps` in magnitude.
    pub eigvec: bool,
    /// A Rayleigh denominator `ev·I·ev` was raised to `eig_eps` in magnitude.
    pub denominator: bool,
    /// `|κmin|` or `|κmax|` was below `curvature_floor`.
    pub curvature: bool,
    /// `κmin κmax < 0` (a saddle; see [`SignRule`]).
    pub saddle: bool,
    /// A principal curvature is exactly zero (see [`SignRule`]).
    pub zero: bool,
    /// Under [`EigenRule::Robust`]: the complemented eigenvector differs from the reference's
    /// `(v1 ∓ disc, v2)` by more than 1e-6 in direction (`|sin|`), i.e. the rule changed the fit.
    pub eigen_switched: bool,
}

/// The intermediate quantities of one fit (orientation, diagnostics).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FitDiagnostics {
    /// Signed curvature of smaller magnitude.
    pub kappa_min: f64,
    /// Signed curvature of larger magnitude.
    pub kappa_max: f64,
    /// Unit principal direction of `kappa_min` (in space, tangent to `Q*` at 0).
    pub v_min: [f64; 3],
    /// Unit principal direction of `kappa_max`.
    pub v_max: [f64; 3],
    /// Unit normal of `Q*` at 0 (Eq. 21).
    pub n_star: [f64; 3],
    /// Clamps that fired.
    pub clamps: ClampFlags,
}

/// How the torus's sign and center follow from the curvatures.
///
/// **Reference** (default; the paper's Eq. 22 and center formula, the JAX training rule and the
/// C++ inference code): `s` by Eq. 22, `c = o − (s/|κmin|) n*`, `R = 1/|κmin| − sgn(κ₊κ₋) r`, with
/// `|κ|` floored but `s` and `sgn` taken from the unfloored values. Two defects, measured in
/// WP1 (NOTES-WP1.md):
/// - at a **saddle** (`κmin κmax < 0`) the center is put on the wrong side: the torus still passes
///   through `o`, but it is the mirror image of the osculating one, and `∂g/∂n* = −1` at `o`
///   (inside and outside swapped locally). The golden tori of the authors' program show it (33 of
///   12,917 tori of holey-bunny, 4 of femur, 1 of snail);
/// - at an **exactly zero** curvature `s = −1` and `sgn = 0`: a plane gives a horn torus whose
///   core circle passes through `o`; a cylinder bending away from `n` is put on the wrong side.
///
/// **Geometric** (ours): `|κ|` below the floor becomes `±floor` with its sign (`+floor` at 0);
/// `s = +1` iff `κmax < 0` (the tube center `o + n*/κmax` is on the inner side); the center is the
/// center of the parallel circle, `c = o + n*/κmin` (signed). Identical to `Reference` except at
/// saddles and exact zeros; the torus then osculates at `o` in every case (outward normal, both
/// principal curvatures).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SignRule {
    /// The reference's rule, as the network was trained (SPEC.md D9).
    #[default]
    Reference,
    /// The osculating torus in every case.
    Geometric,
}

/// How the generalized eigenvector pair of `(II, I)` is formed.
///
/// **Reference** (default, parity): `ev± = (v1 ± disc, v2)`. One of the two cancels when
/// `|v2| ≪ |v1|` (`disc² − v1² = O(v2²)`): `ev₋` for `v1 > 0`, `ev₊` for `v1 < 0`. At `v2 = 0` it is
/// the zero vector, the 1e-10 clamp then returns the *other* eigenvector and both curvatures come
/// out equal (NOTES-WP2.md item 1). In f64 this does not happen on real data; in f32 (training) the
/// radii lose digits below `|v2/v1| ≈ 3e-4` and are wrong by 0.4–0.75 below `3e-7` (NOTES-WP7.md).
///
/// **Robust** (ours): keep the non-cancelling vector of the pair (`ev₊` for `v1 ≥ 0`, else `ev₋`) and
/// take the other as its `I`-orthogonal complement `J·I·e` (`J` the rotation by 90°): the
/// eigenvectors of the symmetric pencil `(II, I)` are `I`-orthogonal, so this is exact whenever the
/// curvatures differ, needs no clamp, and both forms agree (up to sign) at the switch `v1 = 0`, so
/// tori and their derivatives are continuous there.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EigenRule {
    /// The reference's pair (as trained).
    #[default]
    Reference,
    /// The non-cancelling vector and its `I`-orthogonal complement.
    Robust,
}

/// The robust eigenvector pair `(ev₊, ev₋)` (unnormalized) of [`EigenRule::Robust`].
///
/// [formula: pat-principal-pair]
pub fn robust_pair(v1: f64, v2: f64, disc: f64, i: [f64; 3]) -> ([f64; 2], [f64; 2]) {
    let comp = |e: [f64; 2]| [-(i[1] * e[0] + i[2] * e[1]), i[0] * e[0] + i[1] * e[1]];
    if v1 >= 0.0 {
        let e = [v1 + disc, v2];
        (e, comp(e))
    } else {
        let e = [v1 - disc, v2];
        (comp(e), e)
    }
}

/// Options of [`fit_tori`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FitOptions {
    /// Floor of the discriminant, eigenvector components and Rayleigh denominators (the
    /// reference's `NAN_EPSILON`, 1e-10).
    pub eig_eps: f64,
    /// Floor of `|κmin|`, `|κmax|` before inversion (JAX rule, 1e-10).
    pub curvature_floor: f64,
    /// Sign and center rule (see [`SignRule`]).
    pub sign_rule: SignRule,
    /// Eigenvector rule (see [`EigenRule`]).
    pub eigen_rule: EigenRule,
    /// Threads (`None`: all cores). An explicit pool; results do not depend on it.
    pub threads: Option<usize>,
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions { eig_eps: 1e-10, curvature_floor: 1e-10, sign_rule: SignRule::Reference, eigen_rule: EigenRule::Reference, threads: None }
    }
}

/// What fitting did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FitReport {
    /// Points whose discriminant was clamped.
    pub clamped_disc: usize,
    /// Points with an eigenvector-component clamp.
    pub clamped_eigvec: usize,
    /// Points with a Rayleigh-denominator clamp.
    pub clamped_denominator: usize,
    /// Points with a curvature floor.
    pub clamped_curvature: usize,
    /// Umbilic points: the discriminant clamp is active (equal to `clamped_disc`; kept as its own
    /// name for the oracle's flag).
    pub umbilic: usize,
    /// Saddle points (`κmin κmax < 0`; see [`SignRule`]).
    pub saddle: usize,
    /// Points with an exactly zero principal curvature (see [`SignRule`]).
    pub zero_curvature: usize,
    /// Points where [`EigenRule::Robust`] changed the eigenvector (see [`ClampFlags::eigen_switched`]).
    pub eigen_switched: usize,
    /// Points with a non-finite torus (an error is returned when non-empty).
    pub nonfinite: Vec<u32>,
    /// Neighborhoods with fewer than the model's `k` tokens (`n < k`).
    pub tokens_lt_k: usize,
    /// Wall time.
    pub seconds: f64,
}

/// C++ `sign`: +1 for `x ≥ 0` (NaN gives −1).
#[inline]
fn sign_pos(x: f64) -> f64 {
    if x >= 0.0 { 1.0 } else { -1.0 }
}

/// `jnp.sign` / C++ `sgn`: −1, 0 or +1.
#[inline]
fn sgn(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Fundamental forms of the Monge patch `Q*(s,t) = p + s ŝ + t t̂ + Q(s,t) n` at 0:
/// `A = √(1 + a01² + a10²)`, `I = [[1 + a10², a10 a01], [a10 a01, 1 + a01²]]`,
/// `II = [[2 a20, a11], [a11, 2 a02]] / A`. Returns `(A, I, II)` with symmetric matrices as
/// `[m11, m12, m22]`.
///
/// [formula: pat-fundamental-forms]
pub fn fundamental_forms(a: [f64; 6]) -> (f64, [f64; 3], [f64; 3]) {
    let [_a00, a01, a10, a11, a02, a20] = a;
    let big_a = (1.0 + a01 * a01 + a10 * a10).sqrt();
    let i = [1.0 + a10 * a10, a10 * a01, 1.0 + a01 * a01];
    let ii = [2.0 * a20 / big_a, a11 / big_a, 2.0 * a02 / big_a];
    (big_a, i, ii)
}

/// Mean and Gaussian curvature `H = tr(I⁻¹ II)/2`, `K = det II / det I` (Eq. 19).
///
/// [formula: pat-mean-gauss]
pub fn mean_gauss(a: [f64; 6]) -> (f64, f64) {
    let (_, i, ii) = fundamental_forms(a);
    let det_i = i[0] * i[2] - i[1] * i[1];
    let h = 0.5 * (i[2] * ii[0] - 2.0 * i[1] * ii[1] + i[0] * ii[2]) / det_i;
    let k = (ii[0] * ii[2] - ii[1] * ii[1]) / det_i;
    (h, k)
}

/// The unnormalized generalized eigenvector pair of `(II, I)` before clamping, the raw
/// discriminant (before the floor) and `(v1, v2)`: `ev± = (v1 ± disc, v2)` with
/// `v1 = a20 a01² − a02 a10² − a02 + a20`, `v2 = a11(a10² + 1) − 2 a01 a10 a20`,
/// `disc² = (1 + a01² + a10²)(a11² − 4 a02 a20) + term²`,
/// `term = a02(a10² + 1) − a01 a10 a11 + a20(a01² + 1)`.
///
/// [formula: pat-principal-pair]
pub fn principal_pair_terms(a: [f64; 6]) -> (f64, f64, f64) {
    let [_a00, a01, a10, a11, a02, a20] = a;
    let v1 = a20 * a01 * a01 - a02 * a10 * a10 - a02 + a20;
    let v2 = a11 * (a10 * a10 + 1.0) - 2.0 * a01 * a10 * a20;
    let term = a02 * (a10 * a10 + 1.0) - a01 * a10 * a11 + a20 * (a01 * a01 + 1.0);
    let disc2 = (1.0 + a01 * a01 + a10 * a10) * (a11 * a11 - 4.0 * a02 * a20) + term * term;
    (v1, v2, disc2)
}

/// Rayleigh quotient `κ = (evᵀ II ev)/(evᵀ I ev)` with the denominator's magnitude floored at
/// `eps` (C++ sign rule). Returns `(κ, clamped)`.
///
/// [formula: pat-curvature-rayleigh]
pub fn rayleigh(ev: [f64; 2], i: [f64; 3], ii: [f64; 3], eps: f64) -> (f64, bool) {
    let num = ev[0] * (ii[0] * ev[0] + ii[1] * ev[1]) + ev[1] * (ii[1] * ev[0] + ii[2] * ev[1]);
    let den = ev[0] * (i[0] * ev[0] + i[1] * ev[1]) + ev[1] * (i[1] * ev[0] + i[2] * ev[1]);
    let clamped = den.abs() < eps;
    let den = sign_pos(den) * den.abs().max(eps);
    (num / den, clamped)
}

/// The space direction `v = w_x(ŝ + a10 n) + w_y(t̂ + a01 n)` of a parameter direction `w`
/// (tangent to `Q*` at 0).
///
/// [formula: pat-principal-directions]
pub fn principal_direction(w: [f64; 2], s: [f64; 3], t: [f64; 3], n: [f64; 3], a10: f64, a01: f64) -> [f64; 3] {
    let ds = [s[0] + a10 * n[0], s[1] + a10 * n[1], s[2] + a10 * n[2]];
    let dt = [t[0] + a01 * n[0], t[1] + a01 * n[1], t[2] + a01 * n[2]];
    [w[0] * ds[0] + w[1] * dt[0], w[0] * ds[1] + w[1] * dt[1], w[0] * ds[2] + w[1] * dt[2]]
}

/// The normal of `Q*` at 0: `n* = (n − a10 ŝ − a01 t̂)/A` (Eq. 21).
///
/// [formula: pat-surface-normal]
pub fn surface_normal(s: [f64; 3], t: [f64; 3], n: [f64; 3], a10: f64, a01: f64, big_a: f64) -> [f64; 3] {
    [
        (n[0] - a10 * s[0] - a01 * t[0]) / big_a,
        (n[1] - a10 * s[1] - a01 * t[1]) / big_a,
        (n[2] - a10 * s[2] - a01 * t[2]) / big_a,
    ]
}

/// Torus radii from the curvatures: `r = 1/|κmax|`, `R = 1/|κmin| − sgn(κmax κmin) r`, with the
/// magnitudes floored at `floor`. Returns `(R, r, clamped)`.
///
/// [formula: pat-torus-radii]
pub fn torus_radii(kappa_min: f64, kappa_max: f64, floor: f64) -> (f64, f64, bool) {
    let clamped = kappa_max.abs() < floor || kappa_min.abs() < floor;
    let akmax = kappa_max.abs().max(floor);
    let akmin = kappa_min.abs().max(floor);
    let r = 1.0 / akmax;
    let big_r = 1.0 / akmin - sgn(kappa_max * kappa_min) * r;
    (big_r, r, clamped)
}

/// Torus sign (Eq. 22): `+1` iff `κmax < 0 ∧ (κmin < 0 ∨ κmin > 0)`, else `−1`.
///
/// [formula: pat-torus-sign]
#[allow(clippy::nonminimal_bool)] // written as Eq. 22
pub fn torus_sign(kappa_min: f64, kappa_max: f64) -> f64 {
    if (kappa_max < 0.0 && kappa_min < 0.0) || (kappa_max < 0.0 && kappa_min > 0.0) { 1.0 } else { -1.0 }
}

/// Fit one torus: `p`, raw normal `n` (unit by validation), coefficients `a`.
/// Returns the torus and its diagnostics.
///
/// Center `c = p + a00 n − (s/|κmin|) n*`, axis `u = (n* × vmin)/max(‖n* × vmin‖, eps)`,
/// minor radius `s·r`.
///
/// [formula: pat-torus-frame]
/// [formula: pat-torus-contact]
pub fn fit_torus(p: [f64; 3], n: [f64; 3], a: [f64; 6], opts: &FitOptions) -> (Torus, FitDiagnostics) {
    let eps = opts.eig_eps;
    let ln = norm(n);
    let c = [n[0] / ln, n[1] / ln, n[2] / ln];
    let (s_hat, t_hat) = local_frame(c);
    let [a00, a01, a10, _a11, _a02, _a20] = a;
    let (big_a, i, ii) = fundamental_forms(a);
    let (v1, v2, disc2) = principal_pair_terms(a);
    let mut clamps = ClampFlags { disc: disc2 < eps, ..Default::default() };
    let disc = disc2.max(eps).sqrt();

    let mut ev_pos_0 = v1 + disc;
    if ev_pos_0.abs() < eps {
        clamps.eigvec = true;
    }
    ev_pos_0 = sign_pos(ev_pos_0) * ev_pos_0.abs().max(eps);
    let mut ev_neg_1 = v2;
    if ev_neg_1.abs() < eps {
        clamps.eigvec = true;
    }
    ev_neg_1 = sign_pos(ev_neg_1) * ev_neg_1.abs().max(eps);
    let (mut ev_pos, mut ev_neg) = match opts.eigen_rule {
        EigenRule::Reference => ([ev_pos_0, v2], [v1 - disc, ev_neg_1]),
        EigenRule::Robust => {
            clamps.eigvec = false;
            let (p, q) = robust_pair(v1, v2, disc, i);
            // the reference's vector that the rule replaced, for the counter
            let (old, new) = if v1 >= 0.0 { ([v1 - disc, ev_neg_1], q) } else { ([ev_pos_0, v2], p) };
            let (lo, ln) = ((old[0] * old[0] + old[1] * old[1]).sqrt(), (new[0] * new[0] + new[1] * new[1]).sqrt());
            clamps.eigen_switched = ((old[0] * new[1] - old[1] * new[0]) / (lo * ln)).abs() > 1e-6;
            (p, q)
        }
    };
    let lp = (ev_pos[0] * ev_pos[0] + ev_pos[1] * ev_pos[1]).sqrt();
    let lq = (ev_neg[0] * ev_neg[0] + ev_neg[1] * ev_neg[1]).sqrt();
    ev_pos = [ev_pos[0] / lp, ev_pos[1] / lp];
    ev_neg = [ev_neg[0] / lq, ev_neg[1] / lq];

    let v_pos = principal_direction(ev_pos, s_hat, t_hat, n, a10, a01);
    let v_neg = principal_direction(ev_neg, s_hat, t_hat, n, a10, a01);
    let n_star = surface_normal(s_hat, t_hat, n, a10, a01, big_a);

    let (k_pos, cp) = rayleigh(ev_pos, i, ii, eps);
    let (k_neg, cn) = rayleigh(ev_neg, i, ii, eps);
    clamps.denominator = cp || cn;

    let (kappa_min, kappa_max, v_min, v_max) =
        if k_pos.abs() < k_neg.abs() { (k_pos, k_neg, v_pos, v_neg) } else { (k_neg, k_pos, v_neg, v_pos) };

    clamps.zero = kappa_min == 0.0 || kappa_max == 0.0;
    clamps.saddle = kappa_min * kappa_max < 0.0;
    let floor = opts.curvature_floor;
    clamps.curvature = kappa_min.abs() < floor || kappa_max.abs() < floor;
    let (big_r, r, s, center) = match opts.sign_rule {
        SignRule::Reference => {
            let (big_r, r, _) = torus_radii(kappa_min, kappa_max, floor);
            let s = torus_sign(kappa_min, kappa_max);
            let akmin = kappa_min.abs().max(floor);
            let c = [
                p[0] + a00 * n[0] - s / akmin * n_star[0],
                p[1] + a00 * n[1] - s / akmin * n_star[1],
                p[2] + a00 * n[2] - s / akmin * n_star[2],
            ];
            (big_r, r, s, c)
        }
        SignRule::Geometric => {
            let f = |k: f64| if k.abs() >= floor { k } else if k < 0.0 { -floor } else { floor };
            let (km, kx) = (f(kappa_min), f(kappa_max));
            let (big_r, r, _) = torus_radii(km, kx, floor);
            let s = if kx < 0.0 { 1.0 } else { -1.0 };
            // c = o + n*/κmin, written as the reference's expression with the side −sgn(κmin) so
            // that the two rules agree bitwise away from saddles and zeros
            let side = if km < 0.0 { 1.0 } else { -1.0 };
            let akmin = km.abs();
            let c = [
                p[0] + a00 * n[0] - side / akmin * n_star[0],
                p[1] + a00 * n[1] - side / akmin * n_star[1],
                p[2] + a00 * n[2] - side / akmin * n_star[2],
            ];
            (big_r, r, s, c)
        }
    };
    let u = cross(n_star, v_min);
    let lu = norm(u).max(eps);
    let axis = [u[0] / lu, u[1] / lu, u[2] / lu];
    let unitize = |v: [f64; 3]| {
        let l = norm(v);
        [v[0] / l, v[1] / l, v[2] / l]
    };
    let diag = FitDiagnostics { kappa_min, kappa_max, v_min: unitize(v_min), v_max: unitize(v_max), n_star, clamps };
    (Torus { center, axis, major: big_r, minor: s * r }, diag)
}

fn check_inputs(cloud: &PointCloud, coefficients: &[[f32; 6]], opts: &FitOptions) -> Result<(), PatError> {
    if coefficients.len() != cloud.len() {
        return Err(PatError::SizeMismatch {
            what: format!("{} coefficient rows for {} points", coefficients.len(), cloud.len()),
        });
    }
    if !(opts.eig_eps.is_finite() && opts.eig_eps > 0.0) {
        return Err(PatError::InvalidOption { what: format!("FitOptions::eig_eps = {}", opts.eig_eps) });
    }
    if !(opts.curvature_floor.is_finite() && opts.curvature_floor > 0.0) {
        return Err(PatError::InvalidOption { what: format!("FitOptions::curvature_floor = {}", opts.curvature_floor) });
    }
    for (i, a) in coefficients.iter().enumerate() {
        if !a.iter().all(|x| x.is_finite()) {
            return Err(PatError::NonFinite { index: i, what: "coefficient row" });
        }
    }
    Ok(())
}

const FIT_CHUNK: usize = 1024;

/// Fit one torus per point (SPEC.md §7.1). `sigma` and `tokens` are carried into [`Tori`] and the
/// report (`tokens < k` is counted); pass the values from prediction.
pub fn fit_tori_with(
    cloud: &PointCloud,
    coefficients: &[[f32; 6]],
    sigma: &[f64],
    tokens_lt_k: usize,
    opts: &FitOptions,
) -> Result<(Tori, FitReport), PatError> {
    let t0 = std::time::Instant::now();
    check_inputs(cloud, coefficients, opts)?;
    if sigma.len() != cloud.len() {
        return Err(PatError::SizeMismatch { what: format!("{} sigma values for {} points", sigma.len(), cloud.len()) });
    }
    let pool = crate::pool::make_pool(opts.threads)?;
    let (p, nr) = (cloud.positions(), cloud.normals());
    let res: Vec<(Torus, ClampFlags)> = pool.install(|| {
        coefficients
            .par_chunks(FIT_CHUNK)
            .enumerate()
            .flat_map_iter(|(c, rows)| {
                rows.iter().enumerate().map(move |(r, a)| {
                    let i = c * FIT_CHUNK + r;
                    let a = a.map(|x| x as f64);
                    let (t, d) = fit_torus(p[i], nr[i], a, opts);
                    (t, d.clamps)
                })
            })
            .collect()
    });
    let mut report = FitReport { tokens_lt_k, ..Default::default() };
    let mut tori = Vec::with_capacity(res.len());
    for (i, (t, cf)) in res.into_iter().enumerate() {
        report.clamped_disc += cf.disc as usize;
        report.clamped_eigvec += cf.eigvec as usize;
        report.clamped_denominator += cf.denominator as usize;
        report.umbilic += cf.disc as usize;
        report.saddle += cf.saddle as usize;
        report.clamped_curvature += cf.curvature as usize;
        report.zero_curvature += cf.zero as usize;
        report.eigen_switched += cf.eigen_switched as usize;
        let finite = t.center.iter().chain(t.axis.iter()).all(|x| x.is_finite()) && t.major.is_finite() && t.minor.is_finite();
        if !finite {
            report.nonfinite.push(i as u32);
        }
        tori.push(t);
    }
    report.seconds = t0.elapsed().as_secs_f64();
    if !report.nonfinite.is_empty() {
        return Err(PatError::NonFiniteTorus { indices: report.nonfinite });
    }
    Ok((Tori { tori, coefficients: coefficients.to_vec(), sigma: sigma.to_vec() }, report))
}

/// Fit one torus per point from given coefficients (σ unknown: stored as NaN).
pub fn fit_tori(cloud: &PointCloud, coefficients: &[[f32; 6]], opts: &FitOptions) -> Result<(Tori, FitReport), PatError> {
    let sigma = vec![f64::NAN; cloud.len()];
    fit_tori_with(cloud, coefficients, &sigma, 0, opts)
}

/// The fit's intermediate quantities at point `i` (curvatures, directions, normal, clamps).
pub fn fit_diagnostics(cloud: &PointCloud, coefficients: &[[f32; 6]], i: usize, opts: &FitOptions) -> Result<FitDiagnostics, PatError> {
    if i >= cloud.len() || i >= coefficients.len() {
        return Err(PatError::SizeMismatch {
            what: format!("point {i} requested; cloud has {}, coefficients {}", cloud.len(), coefficients.len()),
        });
    }
    let a = coefficients[i].map(|x| x as f64);
    Ok(fit_torus(cloud.positions()[i], cloud.normals()[i], a, opts).1)
}

const MAGIC: &[u8; 8] = b"PATTORI\0";
const VERSION: u32 = 1;

impl Tori {
    /// Our binary format (SPEC.md §2.2): magic `PATTORI\0`, `u32` version 1, `u32` sign-rule marker
    /// (0: not recorded, 1: [`SignRule::Reference`], 2: [`SignRule::Geometric`]; files written
    /// before the marker existed have 0 there), `u64` n,
    /// `u64` hash of the source cloud ([`PointCloud::hash`]), then little-endian arrays: centers
    /// `f64[3n]`, axes `f64[3n]`, major `f64[n]`, minor `f64[n]`, σ `f64[n]`, coefficients `f32[6n]`.
    pub fn to_bytes(&self, cloud_hash: u64) -> Vec<u8> {
        self.to_bytes_with(cloud_hash, None)
    }

    /// [`Tori::to_bytes`] recording the sign rule the tori were fitted with (the tori of the two
    /// rules differ at saddles, so a cache keyed by the cloud alone would mix them).
    pub fn to_bytes_with(&self, cloud_hash: u64, sign_rule: Option<SignRule>) -> Vec<u8> {
        let n = self.tori.len();
        let mut b = Vec::with_capacity(32 + n * (8 * 9 + 24));
        b.extend_from_slice(MAGIC);
        b.extend_from_slice(&VERSION.to_le_bytes());
        let marker: u32 = match sign_rule {
            None => 0,
            Some(SignRule::Reference) => 1,
            Some(SignRule::Geometric) => 2,
        };
        b.extend_from_slice(&marker.to_le_bytes());
        b.extend_from_slice(&(n as u64).to_le_bytes());
        b.extend_from_slice(&cloud_hash.to_le_bytes());
        let mut f = |x: f64| b.extend_from_slice(&x.to_le_bytes());
        self.tori.iter().for_each(|t| t.center.iter().for_each(|&x| f(x)));
        self.tori.iter().for_each(|t| t.axis.iter().for_each(|&x| f(x)));
        self.tori.iter().for_each(|t| f(t.major));
        self.tori.iter().for_each(|t| f(t.minor));
        self.sigma.iter().for_each(|&x| f(x));
        for c in &self.coefficients {
            for x in c {
                b.extend_from_slice(&x.to_le_bytes());
            }
        }
        b
    }

    /// Read [`Tori::to_bytes`]. With `expected_hash`, a different source-cloud hash is an error.
    /// Returns the tori and the stored hash.
    pub fn from_bytes(bytes: &[u8], expected_hash: Option<u64>) -> Result<(Tori, u64), PatError> {
        Self::from_bytes_with(bytes, expected_hash).map(|(t, h, _)| (t, h))
    }

    /// [`Tori::from_bytes`] also returning the recorded sign rule (`None`: not recorded, as in
    /// files written by [`Tori::to_bytes`] or before the marker existed; an unknown marker is an
    /// error).
    pub fn from_bytes_with(bytes: &[u8], expected_hash: Option<u64>) -> Result<(Tori, u64, Option<SignRule>), PatError> {
        let ferr = |reason: String| PatError::Format { reason };
        if bytes.len() < 32 || &bytes[..8] != MAGIC {
            return Err(ferr("not a tori file (bad magic)".into()));
        }
        let u32at = |i: usize| u32::from_le_bytes([bytes[i], bytes[i + 1], bytes[i + 2], bytes[i + 3]]);
        let u64at = |i: usize| u64::from_le_bytes(bytes[i..i + 8].try_into().unwrap_or([0; 8]));
        if u32at(8) != VERSION {
            return Err(ferr(format!("tori file version {} (expected {VERSION})", u32at(8))));
        }
        let sign_rule = match u32at(12) {
            0 => None,
            1 => Some(SignRule::Reference),
            2 => Some(SignRule::Geometric),
            m => return Err(ferr(format!("tori file has sign-rule marker {m} (expected 0, 1 or 2)"))),
        };
        let n = u64at(16) as usize;
        let hash = u64at(24);
        if let Some(h) = expected_hash
            && h != hash
        {
            return Err(ferr(format!("tori file belongs to cloud {hash:016x}, not {h:016x}")));
        }
        let need = n.checked_mul(8 * 9 + 24).and_then(|x| x.checked_add(32));
        if need != Some(bytes.len()) {
            return Err(ferr(format!("tori file has {} bytes, {n} tori need {:?}", bytes.len(), need)));
        }
        let f = |k: usize| f64::from_le_bytes(bytes[32 + 8 * k..40 + 8 * k].try_into().unwrap_or([0; 8]));
        let tori = (0..n)
            .map(|i| Torus {
                center: [f(3 * i), f(3 * i + 1), f(3 * i + 2)],
                axis: [f(3 * n + 3 * i), f(3 * n + 3 * i + 1), f(3 * n + 3 * i + 2)],
                major: f(6 * n + i),
                minor: f(7 * n + i),
            })
            .collect();
        let sigma = (0..n).map(|i| f(8 * n + i)).collect();
        let o = 32 + 8 * 9 * n;
        let g = |k: usize| f32::from_le_bytes([bytes[o + 4 * k], bytes[o + 4 * k + 1], bytes[o + 4 * k + 2], bytes[o + 4 * k + 3]]);
        let coefficients = (0..n).map(|i| std::array::from_fn(|c| g(6 * i + c))).collect();
        Ok((Tori { tori, coefficients, sigma }, hash, sign_rule))
    }
}
