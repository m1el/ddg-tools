//! The training loss (Eq. 27 as coded; SPEC.md §10.1, NOTES-REFERENCE.md §9.1), ported from
//! `SDF_loss` and `torus_distance_gradient_masked` of the reference (points-as-tori, MIT,
//! Copyright (c) 2026 Nicole Feng):
//!
//! ```text
//! L = mean_q |φ(q) − φ_true(q)| + mean_q |1 − ‖∇φ(q)‖|
//! φ(q) = Σ_i w_i g_i(q) / max(W, 1e-10),  w_i = exp(−λ(ρ_i − σ)),  W = Σ_i w_i
//! σ = ½ max of the 64 nearest valid ρ_i,  λ = 64/σ,  ρ_i = ‖q − p_i‖
//! ∇φ = A + φB,  A = Σ w_i(∇g_i − λ r̂_i g_i)/W',  B = λ Σ w_i r̂_i/W'    (frozen λ, σ; App. D)
//! ```
//!
//! over all valid tori of the cloud. The weights depend on `q` and the points only, not on the
//! network, so they are computed here on the host (f64) as constants; with `ŵ_i = w_i/W'` and
//! `c_i = λ(Σ_j ŵ_j r̂_j − r̂_i)`
//!
//! ```text
//! φ = Σ_i ŵ_i g_i,   ∇φ = Σ_i ŵ_i (∇g_i + c_i g_i),
//! ```
//!
//! which is the reference's expression regrouped. Only `g_i` and `∇g_i` depend on the coefficients;
//! they are framework ops ([`super::tfit`] then the closed-form torus SDF and gradient), so the
//! eikonal term needs first-order autodiff only.
//!
//! **Truncation.** A torus with `λ(ρ_i − ρ_min) > T` has `ŵ_i ≤ e^{−T}`; such pairs are dropped
//! (`T = 40` by default: `e^{−40} = 4e−18`, far below f32 round-off). The reference sums all 2048 tori
//! densely; the f32 weights of the dropped pairs are below its own rounding of the sums. `T = ∞`
//! keeps every pair. The dropped weight mass is reported.

use super::tfit::{Ctx, FitGeometry, TensorFitOptions, TorusTensor, fit_tensor, to_f64};
use crate::{PatError, Torus};
use burn::tensor::{Int, Tensor, TensorData};
use rayon::prelude::*;

/// Options of the loss.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LossOptions {
    /// The fit.
    pub fit: TensorFitOptions,
    /// `λσ` (64).
    pub lambda_c: f64,
    /// Nearest valid points that define `σ` (64).
    pub sigma_k: usize,
    /// Truncation `T` of the pair list (40; `f64::INFINITY` keeps all pairs).
    pub cutoff: f64,
    /// Padded pairs per autodiff chunk (memory bound). Queries are sorted by their pair count and
    /// cut into chunks of at most this many padded pairs ([`chunks`]), so padding stays small; the
    /// result depends on it only through the round-off of the final sum.
    pub pair_chunk: usize,
}

impl Default for LossOptions {
    fn default() -> Self {
        LossOptions { fit: TensorFitOptions::default(), lambda_c: 64.0, sigma_k: 64, cutoff: 40.0, pair_chunk: 1 << 16 }
    }
}

/// The host-side constants of the blend for a set of queries (CSR over queries).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BlendPlan {
    /// Query positions (valid queries only).
    pub q: Vec<[f64; 3]>,
    /// Ground truth per query.
    pub truth: Vec<f64>,
    /// Index of each kept query in the caller's list.
    pub source: Vec<u32>,
    /// `σ` per query.
    pub sigma: Vec<f64>,
    /// CSR offsets (`q.len() + 1`).
    pub offset: Vec<u32>,
    /// Torus index per pair (into the valid-torus list given to [`build_plan`]).
    pub idx: Vec<u32>,
    /// `ŵ` per pair.
    pub w: Vec<f64>,
    /// `c` per pair.
    pub c: Vec<[f64; 3]>,
    /// Queries dropped: `σ = 0` or non-finite truth.
    pub dropped: usize,
    /// Largest dropped weight mass `Σ_dropped ŵ` over the queries.
    pub max_dropped_mass: f64,
    /// Largest number of pairs of a query.
    pub max_pairs: usize,
    /// Queries whose `W` was raised to 1e-10 (the reference's guard).
    pub guarded: usize,
}

impl BlendPlan {
    /// Number of queries.
    pub fn len(&self) -> usize {
        self.q.len()
    }

    /// True without queries.
    pub fn is_empty(&self) -> bool {
        self.q.is_empty()
    }

    /// Pairs in all.
    pub fn pairs(&self) -> usize {
        self.idx.len()
    }
}

#[inline]
fn d3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// The truncated plan of one query through the kd-tree: `σ` from the `k` nearest, then only the
/// points within `ρ_min + T/λ` (every other weight is below `e^{−T}` times the largest, so `W` and
/// `Σ ŵ r̂` over this set differ from the full sums by less than `n·e^{−T}` relatively).
fn plan_query_tree(q: [f64; 3], p: &[[f64; 3]], tree: &crate::spatial::KdTree, o: &LossOptions, buf: &mut Vec<crate::spatial::Neighbor>) -> Option<QueryPlan> {
    let n = p.len();
    tree.knn(q, o.sigma_k.min(n), buf);
    let sigma = 0.5 * buf.last()?.d2.sqrt();
    if !(sigma > 0.0 && sigma.is_finite()) {
        return None;
    }
    let lam = o.lambda_c / sigma;
    let rmin = buf[0].d2.sqrt();
    tree.radius(q, rmin + o.cutoff / lam, buf);
    buf.sort_unstable_by_key(|nb| nb.index);
    let rho: Vec<f64> = buf.iter().map(|nb| nb.d2.sqrt()).collect();
    let w: Vec<f64> = rho.iter().map(|&r| (-lam * (r - sigma)).exp()).collect();
    let big_w: f64 = w.iter().sum();
    let guarded = big_w < 1e-10;
    let ws = big_w.max(1e-10);
    let rhat: Vec<[f64; 3]> = buf
        .iter()
        .zip(&rho)
        .map(|(nb, &r)| {
            let d = d3(q, p[nb.index as usize]);
            let l = r.max(1e-10);
            [d[0] / l, d[1] / l, d[2] / l]
        })
        .collect();
    let mut bbar = [0.0; 3];
    for (wi, rh) in w.iter().zip(&rhat) {
        for k in 0..3 {
            bbar[k] += wi / ws * rh[k];
        }
    }
    let wmax = w.iter().copied().fold(0.0, f64::max);
    let dropped_bound = (n - buf.len()) as f64 * wmax * (-o.cutoff).exp() / ws;
    Some(QueryPlan {
        sigma,
        idx: buf.iter().map(|nb| nb.index).collect(),
        c: rhat.iter().map(|rh| std::array::from_fn(|k| lam * (bbar[k] - rh[k]))).collect(),
        w: w.iter().map(|wi| wi / ws).collect(),
        dropped_mass: dropped_bound,
        guarded,
    })
}

struct QueryPlan {
    sigma: f64,
    idx: Vec<u32>,
    w: Vec<f64>,
    c: Vec<[f64; 3]>,
    dropped_mass: f64,
    guarded: bool,
}

fn plan_query(q: [f64; 3], p: &[[f64; 3]], o: &LossOptions, frozen_sigma: Option<f64>, scratch: &mut Vec<f64>) -> Option<QueryPlan> {
    let n = p.len();
    scratch.clear();
    scratch.extend(p.iter().map(|&pi| {
        let d = d3(q, pi);
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }));
    let rho = scratch.clone();
    let sigma = match frozen_sigma {
        Some(s) => s,
        None => {
            let k = o.sigma_k.min(n);
            let (_, kth, _) = scratch.select_nth_unstable_by(k - 1, |a, b| a.total_cmp(b));
            0.5 * *kth
        }
    };
    if !(sigma > 0.0 && sigma.is_finite()) {
        return None;
    }
    let lam = o.lambda_c / sigma;
    let rmin = rho.iter().copied().fold(f64::INFINITY, f64::min);
    let w: Vec<f64> = rho.iter().map(|&r| (-lam * (r - sigma)).exp()).collect();
    let big_w: f64 = w.iter().sum();
    let guarded = big_w < 1e-10;
    let ws = big_w.max(1e-10);
    let mut bbar = [0.0; 3];
    let rhat: Vec<[f64; 3]> = p
        .iter()
        .zip(&rho)
        .map(|(&pi, &r)| {
            let d = d3(q, pi);
            let l = r.max(1e-10);
            [d[0] / l, d[1] / l, d[2] / l]
        })
        .collect();
    for (wi, rh) in w.iter().zip(&rhat) {
        let a = wi / ws;
        for k in 0..3 {
            bbar[k] += a * rh[k];
        }
    }
    let mut out = QueryPlan { sigma, idx: vec![], w: vec![], c: vec![], dropped_mass: 0.0, guarded };
    for i in 0..n {
        let a = w[i] / ws;
        if lam * (rho[i] - rmin) <= o.cutoff {
            out.idx.push(i as u32);
            out.w.push(a);
            out.c.push(std::array::from_fn(|k| lam * (bbar[k] - rhat[i][k])));
        } else {
            out.dropped_mass += a;
        }
    }
    Some(out)
}

/// Build the blend constants for queries `q` (with ground truth `truth`) over the valid tori at
/// points `p`. `frozen_sigma` (tests) fixes `σ` per query. Runs on the current rayon pool.
///
/// [formula: pat-blend-gradient]
/// [formula: pat-blend-weight-range]
pub fn build_plan(q: &[[f64; 3]], truth: &[f64], p: &[[f64; 3]], o: &LossOptions, frozen_sigma: Option<&[f64]>) -> Result<BlendPlan, PatError> {
    if q.len() != truth.len() || frozen_sigma.is_some_and(|s| s.len() != q.len()) {
        return Err(PatError::SizeMismatch { what: format!("{} queries, {} truth values", q.len(), truth.len()) });
    }
    if p.is_empty() {
        return Err(PatError::TooFewPoints { n: 0, need: 1 });
    }
    if !(o.lambda_c > 0.0 && o.sigma_k > 0 && o.cutoff >= 0.0 && o.pair_chunk > 0) {
        return Err(PatError::InvalidOption { what: format!("LossOptions {o:?}") });
    }
    let tree = if o.cutoff.is_finite() && frozen_sigma.is_none() {
        Some(crate::spatial::KdTree::new(p, crate::spatial::KdTreeOptions::default())?)
    } else {
        None
    };
    let per: Vec<Option<QueryPlan>> = q
        .par_iter()
        .enumerate()
        .map_init(|| (Vec::new(), Vec::new()), |(s, b), (j, &x)| {
            if !truth[j].is_finite() || !x.iter().all(|v| v.is_finite()) {
                return None;
            }
            match &tree {
                Some(t) => plan_query_tree(x, p, t, o, b),
                None => plan_query(x, p, o, frozen_sigma.map(|f| f[j]), s),
            }
        })
        .collect();
    let mut plan = BlendPlan { offset: vec![0], ..Default::default() };
    for (j, qp) in per.iter().enumerate() {
        let Some(qp) = qp else {
            plan.dropped += 1;
            continue;
        };
        plan.q.push(q[j]);
        plan.truth.push(truth[j]);
        plan.source.push(j as u32);
        plan.sigma.push(qp.sigma);
        plan.max_pairs = plan.max_pairs.max(qp.idx.len());
        plan.max_dropped_mass = plan.max_dropped_mass.max(qp.dropped_mass);
        plan.guarded += qp.guarded as usize;
        plan.offset.push(plan.offset[plan.offset.len() - 1] + qp.idx.len() as u32);
    }
    // the pair arrays, filled in parallel (same values and order as appending query by query)
    let total = plan.offset[plan.offset.len() - 1] as usize;
    plan.idx = vec![0; total];
    plan.w = vec![0.0; total];
    plan.c = vec![[0.0; 3]; total];
    let kept: Vec<&QueryPlan> = per.iter().flatten().collect();
    let mut parts = Vec::with_capacity(kept.len());
    let (mut ri, mut rw, mut rc) = (&mut plan.idx[..], &mut plan.w[..], &mut plan.c[..]);
    for qp in &kept {
        let l = qp.idx.len();
        let (a, b) = std::mem::take(&mut ri).split_at_mut(l);
        let (x, y) = std::mem::take(&mut rw).split_at_mut(l);
        let (u, v) = std::mem::take(&mut rc).split_at_mut(l);
        parts.push((a, x, u));
        (ri, rw, rc) = (b, y, v);
    }
    parts.into_par_iter().zip(kept.par_iter()).for_each(|((a, x, u), qp)| {
        a.copy_from_slice(&qp.idx);
        x.copy_from_slice(&qp.w);
        u.copy_from_slice(&qp.c);
    });
    Ok(plan)
}

/// φ and ∇φ of the plan's queries from tori on the host (f64, [`crate::sdf`]): the check of the
/// framework path.
pub fn host_blend(tori: &[Torus], plan: &BlendPlan) -> Vec<(f64, [f64; 3])> {
    (0..plan.len())
        .map(|j| {
            let (lo, hi) = (plan.offset[j] as usize, plan.offset[j + 1] as usize);
            let mut phi = 0.0;
            let mut gr = [0.0; 3];
            for k in lo..hi {
                let (g, dg, _) = crate::sdf::torus_sdf_gradient(&tori[plan.idx[k] as usize], plan.q[j]);
                phi += plan.w[k] * g;
                for a in 0..3 {
                    gr[a] += plan.w[k] * (dg[a] + plan.c[k][a] * g);
                }
            }
            (phi, gr)
        })
        .collect()
}

/// One chunk of queries laid out densely: `Q × K` pairs, padded pairs have `ŵ = c = 0` and the
/// torus index of the query's first pair.
struct Dense {
    qn: usize,
    kk: usize,
    q: Vec<f64>,
    idx: Vec<i64>,
    w: Vec<f64>,
    c: Vec<f64>,
    truth: Vec<f64>,
}

fn dense(plan: &BlendPlan, qs: &[u32]) -> Dense {
    let kk = qs.iter().map(|&j| (plan.offset[j as usize + 1] - plan.offset[j as usize]) as usize).max().unwrap_or(1).max(1);
    let qn = qs.len();
    let mut d = Dense {
        qn,
        kk,
        q: Vec::with_capacity(qn * kk * 3),
        idx: Vec::with_capacity(qn * kk),
        w: Vec::with_capacity(qn * kk),
        c: Vec::with_capacity(qn * kk * 3),
        truth: qs.iter().map(|&j| plan.truth[j as usize]).collect(),
    };
    for j in qs.iter().map(|&j| j as usize) {
        let (a, b) = (plan.offset[j] as usize, plan.offset[j + 1] as usize);
        let first = if b > a { plan.idx[a] } else { 0 };
        for k in 0..kk {
            d.q.extend(plan.q[j]);
            if a + k < b {
                d.idx.push(plan.idx[a + k] as i64);
                d.w.push(plan.w[a + k]);
                d.c.extend(plan.c[a + k]);
            } else {
                d.idx.push(first as i64);
                d.w.push(0.0);
                d.c.extend([0.0; 3]);
            }
        }
    }
    d
}

/// φ `[Q, 1]` and ∇φ `[Q, 3]` of a dense chunk from packed tori (framework ops).
///
/// Per pair: `x = q − c`, `v = ‖x × u‖`, `d = (v − R, x·u)`, `g = s(‖d‖ − r)`,
/// `∇g = s(w₁ (x − (x·u)u)/max(v, 1e-10) + w₂ u)`, `w = d/max(‖d‖, 1e-10)` (the reference's guards;
/// the norms are taken of `max(·², 1e-30)` so an exact zero has a finite derivative).
///
/// [formula: pat-torus-sdf]
/// [formula: pat-torus-sdf-gradient]
fn blend_dense(ctx: &Ctx, tori: &TorusTensor, d: &Dense) -> (Tensor<2>, Tensor<2>) {
    let m = d.qn * d.kk;
    let idx = Tensor::<1, Int>::from_data(TensorData::new(d.idx.clone(), [m]), &ctx.dev);
    let t = tori.packed.clone().select(0, idx);
    let c = t.clone().narrow(1, 0, 3);
    let u = t.clone().narrow(1, 3, 3);
    let big_r = t.clone().narrow(1, 6, 1);
    let r = t.clone().narrow(1, 7, 1);
    let s = t.narrow(1, 8, 1);
    let q = ctx.t2(d.q.clone(), [m, 3]);
    let w = ctx.t2(d.w.clone(), [m, 1]);
    let cc = ctx.t2(d.c.clone(), [m, 3]);
    let x = q - c;
    let cr = x.clone().cross(u.clone(), 1);
    let v = cr.square().sum_dim(1).clamp_min(1e-30).sqrt();
    let d1 = v.clone() - big_r;
    let d2 = (x.clone() * u.clone()).sum_dim(1);
    let dn = (d1.clone().square() + d2.clone().square()).clamp_min(1e-30).sqrt();
    let g = s.clone() * (dn.clone() - r);
    let dsafe = dn.clamp_min(1e-10);
    let g1 = (x - d2.clone() * u.clone()) / v.clamp_min(1e-10);
    let grad_g = s * ((d1 / dsafe.clone()) * g1 + (d2 / dsafe) * u);
    let phi = (w.clone() * g.clone()).reshape([d.qn, d.kk]).sum_dim(1);
    let gphi = (w * (grad_g + cc * g)).reshape([d.qn, d.kk, 3]).sum_dim(1).reshape([d.qn, 3]);
    (phi, gphi)
}

/// The queries sorted by `(pairs, index)` and cut into chunks of at most `budget` padded pairs
/// (a chunk holds at least one query).
pub fn chunks(plan: &BlendPlan, budget: usize) -> Vec<Vec<u32>> {
    let cnt = |j: u32| (plan.offset[j as usize + 1] - plan.offset[j as usize]).max(1) as usize;
    let mut order: Vec<u32> = (0..plan.len() as u32).collect();
    order.sort_by_key(|&j| (cnt(j), j));
    let mut out: Vec<Vec<u32>> = Vec::new();
    let mut cur: Vec<u32> = Vec::new();
    for j in order {
        if !cur.is_empty() && (cur.len() + 1) * cnt(j) > budget {
            out.push(std::mem::take(&mut cur));
        }
        cur.push(j);
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// Values of one evaluation of the loss.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LossValue {
    /// `L` (the two means).
    pub loss: f64,
    /// Mean `|φ − φ_true|`.
    pub distance: f64,
    /// Mean `|1 − ‖∇φ‖|`.
    pub eikonal: f64,
    /// Queries in the means.
    pub queries: usize,
    /// Non-finite φ or ∇φ (the reference drops them from the means; a training step with any is
    /// skipped, SPEC.md §10.2).
    pub nonfinite: usize,
}

/// The fitting inputs of a cloud's valid tori: points, normals (as the fit uses them).
pub struct CloudFit<'a> {
    /// Points.
    pub p: &'a [[f64; 3]],
    /// Normals.
    pub n: &'a [[f64; 3]],
}

/// The loss and, if `ctx.dev` is an autodiff device, `∂L/∂a` (`[n, 6]` row major, f64) for
/// coefficients `a` (`n × 6`, f64 values in the dtype of `ctx`). Chunks of at most
/// [`LossOptions::pair_chunk`] padded pairs ([`chunks`]) run on the current rayon pool, each on its
/// own graph; chunk sums are combined in chunk order.
///
/// [formula: pat-loss]
pub fn loss_and_grad(ctx: &Ctx, cloud: &CloudFit<'_>, a: &[f64], plan: &BlendPlan, o: &LossOptions, want_grad: bool) -> Result<(LossValue, Option<Vec<f64>>), PatError> {
    let n = cloud.p.len();
    if a.len() != 6 * n || cloud.n.len() != n {
        return Err(PatError::SizeMismatch { what: format!("{} coefficients for {n} points", a.len()) });
    }
    let nq = plan.len();
    let mut val = LossValue { queries: nq, ..Default::default() };
    if nq == 0 {
        return Ok((val, want_grad.then(|| vec![0.0; 6 * n])));
    }
    let inv = 1.0 / nq as f64;
    let chunks = chunks(plan, o.pair_chunk);
    let outs: Vec<(f64, f64, usize, Option<Vec<f64>>)> = chunks
        .par_iter()
        .map(|qs| {
            let geo = FitGeometry::new(ctx, cloud.p, cloud.n);
            let mut leaf = ctx.t2(a.to_vec(), [n, 6]);
            if want_grad {
                leaf = leaf.autodiff().require_grad();
            }
            let tori = fit_tensor(leaf.clone(), &geo, &o.fit);
            let d = dense(plan, qs);
            let (phi, gphi) = blend_dense(ctx, &tori, &d);
            let truth = ctx.t2(d.truth.clone(), [d.qn, 1]);
            let err = (phi - truth).abs();
            let gn = gphi.square().sum_dim(1).clamp_min(1e-30).sqrt();
            let eik = (gn.neg() + 1.0).abs();
            let ev = to_f64(err.clone());
            let ek = to_f64(eik.clone());
            let bad = ev.iter().zip(&ek).filter(|(x, y)| !x.is_finite() || !y.is_finite()).count();
            let sd: f64 = ev.iter().sum();
            let se: f64 = ek.iter().sum();
            let grad = if want_grad {
                let l = (err.sum() + eik.sum()).mul_scalar(inv);
                let grads = l.backward();
                Some(leaf.grad(&grads).map(to_f64).unwrap_or_else(|| vec![0.0; 6 * n]))
            } else {
                None
            };
            (sd, se, bad, grad)
        })
        .collect();
    let mut grad = want_grad.then(|| vec![0.0; 6 * n]);
    let (mut sd, mut se) = (0.0, 0.0);
    for (a_, b_, bad, g) in outs {
        sd += a_;
        se += b_;
        val.nonfinite += bad;
        if let (Some(acc), Some(g)) = (grad.as_mut(), g) {
            for (x, y) in acc.iter_mut().zip(g) {
                *x += y;
            }
        }
    }
    val.distance = sd * inv;
    val.eikonal = se * inv;
    val.loss = val.distance + val.eikonal;
    Ok((val, grad))
}

/// The tori of coefficients `a` computed by the framework fit, read back (tests, diagnostics):
/// `(center, axis, R, signed minor)`.
pub fn tensor_tori(ctx: &Ctx, cloud: &CloudFit<'_>, a: &[f64], o: &TensorFitOptions) -> Vec<Torus> {
    let n = cloud.p.len();
    let geo = FitGeometry::new(ctx, cloud.p, cloud.n);
    let t = fit_tensor(ctx.t2(a.to_vec(), [n, 6]), &geo, o);
    let v = to_f64(t.packed);
    (0..n)
        .map(|i| {
            let r = &v[9 * i..9 * i + 9];
            Torus { center: [r[0], r[1], r[2]], axis: [r[3], r[4], r[5]], major: r[6], minor: r[8] * r[7] }
        })
        .collect()
}
