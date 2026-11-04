//! Torus fitting as framework ops (SPEC.md §10.4): the reference's `torus_precompute_masked`
//! (`src/pointsastori/network.py`, points-as-tori, MIT, Copyright (c) 2026 Nicole Feng) with the JAX
//! clamp rule, differentiable in the coefficients, on any burn device and in its float dtype (f32
//! for training, f64 for the finite-difference checks).
//!
//! The math is that of [`crate::fit::fit_torus`] (formula ids there), with the same
//! [`SignRule`] and [`EigenRule`] options. Two differences from the JAX code, neither changing a
//! value: the axis is kept as a unit vector (JAX converts it to angles and back, whose `acos` has an
//! infinite derivative at `|u_z| = 1`), and `sgn⁺` (C++) replaces `jnp.sign` in the clamps, as in
//! [`crate::fit`] (they differ only at exact zeros).

use crate::fit::{EigenRule, SignRule};
use burn::tensor::{Bool, DType, Device, Tensor, TensorData};

/// A device and float dtype.
#[derive(Debug, Clone)]
pub struct Ctx {
    /// The device (plain or autodiff).
    pub dev: Device,
    /// The float dtype of everything created here.
    pub dtype: DType,
}

impl Ctx {
    /// A tensor from f64 host values (rounded to the dtype).
    pub fn t2(&self, v: Vec<f64>, shape: [usize; 2]) -> Tensor<2> {
        Tensor::<2>::from_data(TensorData::new(v, shape), (&self.dev, self.dtype))
    }

    /// A tensor from f32 host values (exact in either dtype).
    pub fn t2_f32(&self, v: Vec<f32>, shape: [usize; 2]) -> Tensor<2> {
        Tensor::<2>::from_data(TensorData::new(v, shape), (&self.dev, self.dtype))
    }

    /// A 3-D tensor from f64 values.
    pub fn t3(&self, v: Vec<f64>, shape: [usize; 3]) -> Tensor<3> {
        Tensor::<3>::from_data(TensorData::new(v, shape), (&self.dev, self.dtype))
    }
}

/// Read a tensor back as f64 values.
pub fn to_f64<const D: usize>(t: Tensor<D>) -> Vec<f64> {
    let d = t.into_data();
    d.convert::<f64>().try_to_vec::<f64>().unwrap_or_default()
}

/// The constant per-point geometry of the fit: positions, the normal used by the fitting formulas
/// (raw or unit, D23), and the local frame `ŝ, t̂` of the unit normal; all `[n, 3]`.
#[derive(Debug, Clone)]
pub struct FitGeometry {
    /// Points `p`.
    pub p: Tensor<2>,
    /// Normals `n` as the fitting uses them.
    pub n: Tensor<2>,
    /// `ŝ`.
    pub s: Tensor<2>,
    /// `t̂`.
    pub t: Tensor<2>,
}

impl FitGeometry {
    /// From host points and normals; the frame is [`crate::features::local_frame`] of `n/‖n‖`.
    pub fn new(ctx: &Ctx, p: &[[f64; 3]], n: &[[f64; 3]]) -> FitGeometry {
        let m = p.len();
        let mut s = Vec::with_capacity(3 * m);
        let mut t = Vec::with_capacity(3 * m);
        for v in n {
            let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(crate::features::NORMAL_EPS);
            let (a, b) = crate::features::local_frame([v[0] / l, v[1] / l, v[2] / l]);
            s.extend(a);
            t.extend(b);
        }
        let flat = |x: &[[f64; 3]]| x.iter().flatten().copied().collect::<Vec<f64>>();
        FitGeometry { p: ctx.t2(flat(p), [m, 3]), n: ctx.t2(flat(n), [m, 3]), s: ctx.t2(s, [m, 3]), t: ctx.t2(t, [m, 3]) }
    }
}

/// The fitted tori: `[n, 9]` rows `(center 3, axis 3, R, r = |minor|, sign)`. The sign is a
/// constant (piecewise constant in the coefficients).
#[derive(Debug, Clone)]
pub struct TorusTensor {
    /// `[n, 9]`.
    pub packed: Tensor<2>,
}

/// Options of [`fit_tensor`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensorFitOptions {
    /// Sign and center rule. Training defaults to [`SignRule::Reference`] (the released weights were
    /// trained with it), independent of the inference default [`SignRule::default`] (`Geometric`).
    pub sign_rule: SignRule,
    /// Eigenvector rule.
    pub eigen_rule: EigenRule,
    /// Clamp of the discriminant, eigenvector components and denominators (1e-10).
    pub eig_eps: f64,
    /// Floor of `|κ|` (1e-10).
    pub curvature_floor: f64,
}

impl Default for TensorFitOptions {
    /// Explicitly `SignRule::Reference` and `EigenRule::Reference` (as the authors trained), not
    /// `SignRule::default()`.
    fn default() -> Self {
        TensorFitOptions { sign_rule: SignRule::Reference, eigen_rule: EigenRule::Reference, eig_eps: 1e-10, curvature_floor: 1e-10 }
    }
}

type T = Tensor<2>;

/// `+1` for `x ≥ 0`, else `−1` (a constant).
fn sgn_pos(x: &T) -> T {
    let x = x.clone().detach();
    let m = x.clone().greater_equal_elem(0.0);
    x.ones_like().neg().mask_fill(m, 1.0)
}

/// `sgn⁺(x)·max(|x|, eps)`.
fn clamp_sign(x: T, eps: f64) -> T {
    sgn_pos(&x) * x.abs().clamp_min(eps)
}

/// `jnp.sign` (0 at 0), a constant.
fn sgn(x: &T) -> T {
    x.clone().detach().sign()
}

fn where3(mask: &Tensor<2, Bool>, yes: T, no: T) -> T {
    let [n, c] = no.dims();
    no.mask_where(mask.clone().expand([n, c]), yes)
}

/// Generalized Rayleigh quotient `(eᵀ II e)/(eᵀ I e)` with the denominator clamped.
fn rayleigh(ex: &T, ey: &T, i: &(T, T, T), ii: &(T, T, T), eps: f64) -> T {
    let num = ex.clone() * (ii.0.clone() * ex.clone() + ii.1.clone() * ey.clone()) + ey.clone() * (ii.1.clone() * ex.clone() + ii.2.clone() * ey.clone());
    let den = ex.clone() * (i.0.clone() * ex.clone() + i.1.clone() * ey.clone()) + ey.clone() * (i.1.clone() * ex.clone() + i.2.clone() * ey.clone());
    num / clamp_sign(den, eps)
}

/// Fit the tori of coefficients `a` (`[n, 6]`, `(a00, a01, a10, a11, a02, a20)`).
///
/// [formula: pat-fundamental-forms]
/// [formula: pat-principal-pair]
/// [formula: pat-curvature-rayleigh]
/// [formula: pat-surface-normal]
/// [formula: pat-torus-radii]
/// [formula: pat-torus-sign]
/// [formula: pat-torus-frame]
pub fn fit_tensor(a: T, g: &FitGeometry, o: &TensorFitOptions) -> TorusTensor {
    let eps = o.eig_eps;
    let floor = o.curvature_floor;
    let col = |i: usize| a.clone().narrow(1, i, 1);
    let (a00, a01, a10, a11, a02, a20) = (col(0), col(1), col(2), col(3), col(4), col(5));
    let a01s = a01.clone().square();
    let a10s = a10.clone().square();
    let det_i = a01s.clone() + a10s.clone() + 1.0; // A²
    let big_a = det_i.clone().sqrt();
    let v1 = a20.clone() * a01s.clone() - a02.clone() * a10s.clone() - a02.clone() + a20.clone();
    let v2 = a11.clone() * (a10s.clone() + 1.0) - a01.clone() * a10.clone() * a20.clone() * 2.0;
    let term = a02.clone() * (a10s.clone() + 1.0) - a01.clone() * a10.clone() * a11.clone() + a20.clone() * (a01s.clone() + 1.0);
    let disc2 = det_i.clone() * (a11.clone().square() - a02.clone() * a20.clone() * 4.0) + term.square();
    let disc = disc2.clamp_min(eps).sqrt();

    let i_f = (a10s.clone() + 1.0, a10.clone() * a01.clone(), a01s.clone() + 1.0);
    let ii_f = (a20.clone() * 2.0 / big_a.clone(), a11.clone() / big_a.clone(), a02.clone() * 2.0 / big_a.clone());

    let (pex, pey, nex, ney) = match o.eigen_rule {
        EigenRule::Reference => {
            (clamp_sign(v1.clone() + disc.clone(), eps), v2.clone(), v1.clone() - disc.clone(), clamp_sign(v2.clone(), eps))
        }
        EigenRule::Robust => {
            let comp = |ex: &T, ey: &T| {
                (
                    -(i_f.1.clone() * ex.clone() + i_f.2.clone() * ey.clone()),
                    i_f.0.clone() * ex.clone() + i_f.1.clone() * ey.clone(),
                )
            };
            // v1 ≥ 0: e₊ = (v1 + disc, v2), e₋ = J I e₊; else e₋ = (v1 − disc, v2), e₊ = J I e₋
            let m = v1.clone().detach().greater_equal_elem(0.0);
            let (p1x, p1y) = (v1.clone() + disc.clone(), v2.clone());
            let (n1x, n1y) = comp(&p1x, &p1y);
            let (n2x, n2y) = (v1.clone() - disc.clone(), v2.clone());
            let (p2x, p2y) = comp(&n2x, &n2y);
            (where3(&m, p1x, p2x), where3(&m, p1y, p2y), where3(&m, n1x, n2x), where3(&m, n1y, n2y))
        }
    };
    let lp = (pex.clone().square() + pey.clone().square()).sqrt();
    let ln = (nex.clone().square() + ney.clone().square()).sqrt();
    let (pex, pey) = (pex / lp.clone(), pey / lp);
    let (nex, ney) = (nex / ln.clone(), ney / ln);

    let ds = g.s.clone() + a10.clone() * g.n.clone();
    let dt = g.t.clone() + a01.clone() * g.n.clone();
    let v_pos = pex.clone() * ds.clone() + pey.clone() * dt.clone();
    let v_neg = nex.clone() * ds + ney.clone() * dt;
    let n_star = (g.n.clone() - a10 * g.s.clone() - a01 * g.t.clone()) / big_a;

    let k_pos = rayleigh(&pex, &pey, &i_f, &ii_f, eps);
    let k_neg = rayleigh(&nex, &ney, &i_f, &ii_f, eps);
    let use_pos = k_pos.clone().detach().abs().lower(k_neg.clone().detach().abs());
    let kmin = where3(&use_pos, k_pos.clone(), k_neg.clone());
    let kmax = where3(&use_pos, k_neg, k_pos);
    let vmin = where3(&use_pos, v_pos, v_neg);

    let ones = kmin.clone().detach().ones_like();
    let (big_r, r, s, side, akmin) = match o.sign_rule {
        SignRule::Reference => {
            let akmax = kmax.clone().abs().clamp_min(floor);
            let akmin = kmin.clone().abs().clamp_min(floor);
            let r = akmax.recip();
            let big_r = akmin.clone().recip() - sgn(&(kmax.clone() * kmin.clone())) * r.clone();
            let kd = (kmax.clone().detach(), kmin.clone().detach());
            let neg_max = kd.0.lower_elem(0.0);
            let pos = neg_max.clone().bool_and(kd.1.clone().lower_elem(0.0)).bool_or(neg_max.bool_and(kd.1.greater_elem(0.0)));
            let s = ones.clone().neg().mask_fill(pos, 1.0);
            (big_r, r, s.clone(), s, akmin)
        }
        SignRule::Geometric => {
            let km = clamp_sign(kmin.clone(), floor);
            let kx = clamp_sign(kmax.clone(), floor);
            let r = kx.clone().abs().recip();
            let akmin = km.clone().abs();
            let big_r = akmin.clone().recip() - sgn(&(km.clone() * kx.clone())) * r.clone();
            let s = ones.clone().neg().mask_fill(kx.detach().lower_elem(0.0), 1.0);
            let side = ones.neg().mask_fill(km.detach().lower_elem(0.0), 1.0);
            (big_r, r, s, side, akmin)
        }
    };
    let center = g.p.clone() + a00 * g.n.clone() - (side / akmin) * n_star.clone();
    let u = n_star.cross(vmin, 1);
    let lu = u.clone().square().sum_dim(1).sqrt().clamp_min(eps);
    let axis = u / lu;
    TorusTensor { packed: Tensor::cat(vec![center, axis, big_r, r, s], 1) }
}
