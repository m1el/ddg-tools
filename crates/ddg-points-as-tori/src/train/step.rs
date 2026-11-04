//! One training step (SPEC.md §10.4), in passes:
//!
//! 1. the blend constants of the step's queries over the cloud's valid tori (host, f64,
//!    [`super::loss::build_plan`]); they fix which tori the loss can see;
//! 2. coefficients of those tori without the tape (the network on their neighborhoods);
//! 3. `∂L/∂a` from the loss with the tape on the coefficients only ([`super::loss::loss_and_grad`]);
//! 4. the network again, with the tape, in chunks of [`StepOptions::network_chunk`] neighborhoods,
//!    back-propagating `⟨∂L/∂y, y(θ)⟩` (`y` the raw outputs, `∂L/∂y = ∂L/∂a ⊙ ∂a/∂y`); rows whose
//!    `∂L/∂a` is exactly zero are skipped (their contribution is exactly zero).
//!
//! This is exact (the loss depends on θ only through `a`), needs the activations of one chunk at a
//! time, and on the CPU runs the chunks on a fixed number of groups ([`StepOptions::groups`]) whose
//! gradients are summed in group order, so the result does not depend on the thread count.
//!
//! Padded clouds (D19): by default the network runs only on valid rows, short neighborhoods are
//! padded to `k` tokens with an attention pad mask. `reference_padding` reproduces the reference:
//! the network runs on all `P` rows, `−1` neighbors read row `P − 1`, no mask, and non-finite
//! gradient entries are zeroed.

use super::data::TrainShape;
use super::loss::{BlendPlan, CloudFit, LossOptions, LossValue, build_plan, loss_and_grad};
use super::net::{NetParams, PT, RowFeatures, forward};
use super::queries::{QueryOptions, RowQueries, sample_rows};
use super::tfit::{Ctx, to_f64};
use super::truth::MeshTruth;
use crate::features::{FeatureOrder, FeatureScratch, neighborhood_features, neighborhood_features_f64};
use crate::network::{scale_outputs, scale_outputs_f64};
use crate::PatError;
use burn::tensor::DType;
use rayon::prelude::*;

/// Options of a step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepOptions {
    /// The loss.
    pub loss: LossOptions,
    /// The queries.
    pub queries: QueryOptions,
    /// Feature arithmetic (D3; `Reference` for parity with the trained regime).
    pub feature_order: FeatureOrder,
    /// Reproduce the reference's padding (D19).
    pub reference_padding: bool,
    /// Neighborhoods per taped network chunk (pass 4; on the GPU also the fixed batch shape).
    pub network_chunk: usize,
    /// Neighborhoods per launch of the untaped forward pass on the GPU (pass 2); 0: `network_chunk`.
    /// Ignored on the CPU.
    pub forward_chunk: usize,
    /// Fixed number of gradient groups (pass 4); results depend on it, not on the threads.
    pub groups: usize,
}

impl Default for StepOptions {
    fn default() -> Self {
        StepOptions {
            loss: LossOptions::default(),
            queries: QueryOptions::default(),
            feature_order: FeatureOrder::DifferenceFirst,
            reference_padding: false,
            network_chunk: 4,
            forward_chunk: 0,
            groups: 64,
        }
    }
}

/// A shape prepared for an epoch: features of the network rows, the fit inputs, the queries and
/// their ground truth.
#[derive(Debug, Clone)]
pub struct PreparedShape {
    /// Rows the network runs on.
    pub net_rows: Vec<u32>,
    /// Their features.
    pub feats: RowFeatures,
    /// Point and normal per network row (the fit's inputs).
    pub fit_p: Vec<[f64; 3]>,
    /// Normals as the fit uses them.
    pub fit_n: Vec<[f64; 3]>,
    /// Network-row positions (`net_rows` order) of the valid tori.
    pub valid: Vec<u32>,
    /// Queries of every valid row (rows `0..n`).
    pub queries: RowQueries,
    /// Points `P`.
    pub rows: usize,
    /// Valid points `n`.
    pub n: usize,
}

/// Features of `rows` of a shape (f64 values; the reference order's are its f32 values).
pub fn row_features(shape: &TrainShape, rows: &[u32], order: FeatureOrder, reference_padding: bool) -> RowFeatures {
    let k = shape.k;
    let b = rows.len();
    let mut f = RowFeatures { b, k, feats: vec![0.0; b * k * 6], sigma: vec![0.0; b], pad: None };
    let mut pad = vec![false; b * k];
    let mut any = false;
    let per: Vec<(f64, Vec<f64>, usize)> = rows
        .par_iter()
        .map_init(FeatureScratch::default, |scr, &i| {
            let row = shape.neighbor_row(i as usize);
            let idx: Vec<u32> = if reference_padding {
                row.iter().map(|&j| if j < 0 { (shape.rows - 1) as u32 } else { j as u32 }).collect()
            } else {
                row.iter().take_while(|&&j| j >= 0).map(|&j| j as u32).collect()
            };
            let t = idx.len();
            let mut out = vec![0.0; k * 6];
            let s = match order {
                FeatureOrder::DifferenceFirst => neighborhood_features_f64(&shape.positions, &shape.normals, &idx, scr, &mut out[..t * 6]),
                FeatureOrder::Reference => {
                    let mut o32 = vec![0f32; t * 6];
                    let s = neighborhood_features(&shape.positions, &shape.normals, &idx, order, scr, &mut o32);
                    for (o, v) in out.iter_mut().zip(o32) {
                        *o = v as f64;
                    }
                    s
                }
            };
            (s, out, t)
        })
        .collect();
    for (r, (s, out, t)) in per.into_iter().enumerate() {
        f.sigma[r] = s;
        f.feats[r * k * 6..(r + 1) * k * 6].copy_from_slice(&out);
        if t < k {
            any = true;
            pad[r * k + t..(r + 1) * k].iter_mut().for_each(|x| *x = true);
        }
    }
    if any {
        f.pad = Some(pad);
    }
    f
}

/// Prepare a shape for an epoch: features and the queries and ground truth of every valid row,
/// keyed by `(seed, epoch, shape_id)`.
pub fn prepare_shape(shape: &TrainShape, o: &StepOptions, key: (u64, u64, u64)) -> Result<PreparedShape, PatError> {
    let n = shape.n;
    if n == 0 {
        return Err(PatError::TooFewPoints { n: 0, need: 1 });
    }
    let net_rows: Vec<u32> = if o.reference_padding { (0..shape.rows as u32).collect() } else { (0..n as u32).collect() };
    let feats = row_features(shape, &net_rows, o.feature_order, o.reference_padding);
    let fit_p: Vec<[f64; 3]> = net_rows.iter().map(|&i| shape.positions[i as usize]).collect();
    let fit_n: Vec<[f64; 3]> = net_rows.iter().map(|&i| shape.normals[i as usize]).collect();
    let valid: Vec<u32> = (0..net_rows.len() as u32).filter(|&r| (net_rows[r as usize] as usize) < n).collect();
    let mesh = MeshTruth::new(&shape.vertices, &shape.faces);
    let rows_q: Vec<u32> = (0..n as u32).collect();
    let padding = o.reference_padding;
    let columns = |i: usize| -> Vec<usize> {
        let row = shape.neighbor_row(i);
        if padding {
            row.iter().map(|&j| if j < 0 { shape.rows - 1 } else { j as usize }).collect()
        } else {
            row.iter().take_while(|&&j| j >= 0).map(|&j| j as usize).collect()
        }
    };
    let queries = sample_rows(&shape.positions, &rows_q, &columns, &mesh, &shape.vertices, &shape.faces, &o.queries, key)?;
    Ok(PreparedShape { net_rows, feats, fit_p, fit_n, valid, queries, rows: shape.rows, n })
}

/// The output of a step.
#[derive(Debug, Clone, Default)]
pub struct StepOut {
    /// The loss before the update.
    pub value: LossValue,
    /// Parameter gradients (canonical order), `None` if not requested or the step is empty.
    pub grads: Option<Vec<PT>>,
    /// Gradient entries that were non-finite (zeroed under `reference_padding`; otherwise the step
    /// is to be skipped).
    pub nonfinite_grads: usize,
    /// Network rows run without / with the tape.
    pub rows_forward: usize,
    /// Network rows run with the tape.
    pub rows_backward: usize,
    /// Pairs in the blend.
    pub pairs: usize,
    /// The plan's diagnostics.
    pub plan_dropped: usize,
    /// Largest truncated weight mass.
    pub max_dropped_mass: f64,
    /// Seconds: plan, forward, loss, backward.
    pub seconds: [f64; 4],
}

/// Coefficients (f64 values of the network dtype) and `∂a/∂y` of network rows `sel`
/// (positions in `prep.net_rows`), without the tape.
fn coefficients(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, sel: &[u32], order: FeatureOrder, chunk: usize) -> Result<(Vec<f64>, Vec<[f64; 6]>), PatError> {
    let (net, _) = params.build(ctx, false)?;
    // gather the selected rows' features into one buffer
    let sub = select_rows(&prep.feats, sel);
    let chunks: Vec<(usize, usize)> = (0..sel.len()).step_by(chunk).map(|lo| (lo, (lo + chunk).min(sel.len()))).collect();
    let raw: Vec<Vec<f64>> = chunks.par_iter().map(|&(lo, hi)| to_f64(forward(&net, ctx, &sub, lo, hi))).collect();
    Ok(scale_raw(&raw.concat(), &sub.sigma, ctx, order))
}

/// Coefficients `a` and `∂a/∂y` from raw outputs `y` (6 per row) and the rows' `σ`.
pub(crate) fn scale_raw(raw: &[f64], sigma: &[f64], ctx: &Ctx, order: FeatureOrder) -> (Vec<f64>, Vec<[f64; 6]>) {
    let mut a = Vec::with_capacity(raw.len());
    let mut da = Vec::with_capacity(raw.len() / 6);
    for (r, y) in raw.chunks_exact(6).enumerate() {
        let s = sigma[r];
        let y6 = [y[0], y[1], y[2], y[3], y[4], y[5]];
        let v = if ctx.dtype == DType::F64 {
            scale_outputs_f64(y6, s)
        } else {
            scale_outputs(y6.map(|v| v as f32), s, order).map(|v| v as f64)
        };
        a.extend(v);
        da.push([s, 1.0, 1.0, 1.0 / s, 1.0 / s, 1.0 / s]);
    }
    (a, da)
}

/// The features of network rows `sel` (positions in the shape's network rows), in that order.
pub(crate) fn select_rows(f: &RowFeatures, sel: &[u32]) -> RowFeatures {
    let k = f.k;
    let mut out = RowFeatures { b: sel.len(), k, feats: Vec::with_capacity(sel.len() * k * 6), sigma: Vec::with_capacity(sel.len()), pad: None };
    let mut pad = Vec::with_capacity(sel.len() * k);
    for &r in sel {
        let r = r as usize;
        out.feats.extend_from_slice(&f.feats[r * k * 6..(r + 1) * k * 6]);
        out.sigma.push(f.sigma[r]);
        if let Some(p) = &f.pad {
            pad.extend_from_slice(&p[r * k..(r + 1) * k]);
        }
    }
    if f.pad.is_some() {
        out.pad = Some(pad);
    }
    out
}

/// Parameter gradients of `Σ_rows ⟨g_y, y(θ)⟩` for network rows `sel` with output weights `gy`.
fn backprop(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, sel: &[u32], gy: &[f64], o: &StepOptions) -> Result<Vec<PT>, PatError> {
    let sub = select_rows(&prep.feats, sel);
    let nrow = sel.len();
    let groups = o.groups.max(1).min(nrow.max(1));
    let bounds: Vec<(usize, usize)> = (0..groups).map(|g| (g * nrow / groups, (g + 1) * nrow / groups)).collect();
    let per: Vec<Result<Option<Vec<PT>>, PatError>> = bounds
        .par_iter()
        .map(|&(glo, ghi)| {
            if glo == ghi {
                return Ok(None);
            }
            let (net, leaves) = params.build(ctx, true)?;
            let mut acc: Option<Vec<PT>> = None;
            let mut lo = glo;
            while lo < ghi {
                let hi = (lo + o.network_chunk).min(ghi);
                let y = forward(&net, ctx, &sub, lo, hi);
                let w = ctx.t2(gy[6 * lo..6 * hi].to_vec(), [hi - lo, 6]);
                let grads = (y * w).sum().backward();
                let g: Vec<PT> = leaves
                    .iter()
                    .map(|l| match l {
                        PT::D1(t) => PT::D1(t.grad(&grads).unwrap_or_else(|| t.clone().detach().without_autodiff().zeros_like())),
                        PT::D2(t) => PT::D2(t.grad(&grads).unwrap_or_else(|| t.clone().detach().without_autodiff().zeros_like())),
                    })
                    .collect();
                acc = Some(match acc {
                    None => g,
                    Some(a) => a.iter().zip(&g).map(|(x, y)| x.add(y)).collect(),
                });
                lo = hi;
            }
            Ok(acc)
        })
        .collect();
    let mut total: Option<Vec<PT>> = None;
    for g in per {
        if let Some(g) = g? {
            total = Some(match total {
                None => g,
                Some(t) => t.iter().zip(&g).map(|(x, y)| x.add(y)).collect(),
            });
        }
    }
    total.ok_or_else(|| PatError::Backend { reason: "no rows to back-propagate".into() })
}

/// Queries of the batch of rows `[lo, hi)` (rows ≥ n have none).
fn batch_queries(prep: &PreparedShape, lo: usize, hi: usize) -> (Vec<[f64; 3]>, Vec<f64>) {
    let m = prep.queries.points.len() / prep.n.max(1);
    let (lo, hi) = (lo.min(prep.n), hi.min(prep.n));
    (prep.queries.points[lo * m..hi * m].to_vec(), prep.queries.truth[lo * m..hi * m].to_vec())
}

/// The host part of a step (parameter independent, so it can be computed ahead on CPU threads
/// while the device trains): the blend plan of the slice's queries and the network rows it needs.
#[derive(Debug, Clone)]
pub struct SlicePlan {
    /// The blend constants.
    pub plan: BlendPlan,
    /// Network rows (positions in `prep.net_rows`) whose tori the plan references.
    pub sel_net: Vec<u32>,
    /// The pair layout of the GPU loss (feature `cuda`, built when asked for).
    #[cfg(feature = "cuda")]
    pub gpu: Option<super::gpu::PairLayout>,
    /// Seconds spent building it.
    pub seconds: f64,
}

/// The host part of the step on query rows `[lo, hi)`; `None` when the slice has no queries.
/// `gpu_layout` also builds the GPU loss's pair layout (ignored without the feature `cuda`).
/// Runs on the current rayon pool.
pub fn plan_slice(prep: &PreparedShape, lo: usize, hi: usize, o: &StepOptions, gpu_layout: bool) -> Result<Option<SlicePlan>, PatError> {
    let t0 = std::time::Instant::now();
    let (q, truth) = batch_queries(prep, lo, hi);
    if q.is_empty() {
        return Ok(None);
    }
    let vp: Vec<[f64; 3]> = prep.valid.iter().map(|&r| prep.fit_p[r as usize]).collect();
    let plan: BlendPlan = build_plan(&q, &truth, &vp, &o.loss, None)?;
    // tori the loss sees (valid-list positions) → network rows
    let mut used = vec![false; vp.len()];
    plan.idx.iter().for_each(|&i| used[i as usize] = true);
    let sel_valid: Vec<u32> = if o.reference_padding { (0..vp.len() as u32).collect() } else { (0..vp.len() as u32).filter(|&i| used[i as usize]).collect() };
    let sel_net: Vec<u32> = if o.reference_padding { (0..prep.net_rows.len() as u32).collect() } else { sel_valid.iter().map(|&i| prep.valid[i as usize]).collect() };
    #[cfg(feature = "cuda")]
    let gpu = if gpu_layout { Some(super::gpu::pair_layout(&plan, vp.len(), &o.loss)) } else { None };
    #[cfg(not(feature = "cuda"))]
    let _ = gpu_layout;
    Ok(Some(SlicePlan {
        plan,
        sel_net,
        #[cfg(feature = "cuda")]
        gpu,
        seconds: t0.elapsed().as_secs_f64(),
    }))
}

/// True when `ctx` is a CUDA device (feature `cuda`): the step then runs the GPU passes of
/// [`super::gpu`] (fixed batch shapes, one dispatching thread, gathers only).
pub fn on_gpu(ctx: &Ctx) -> bool {
    #[cfg(feature = "cuda")]
    {
        super::gpu::is_cuda(&ctx.dev)
    }
    #[cfg(not(feature = "cuda"))]
    {
        let _ = ctx;
        false
    }
}

/// The loss of query rows `[lo, hi)` with the current parameters, and the parameter gradients if
/// `want_grad`. `ctx` is the plain device and dtype (autodiff is enabled on the leaves).
/// [`plan_slice`] followed by [`step_planned`].
pub fn step(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, lo: usize, hi: usize, o: &StepOptions, want_grad: bool) -> Result<StepOut, PatError> {
    match plan_slice(prep, lo, hi, o, on_gpu(ctx))? {
        Some(sp) => step_planned(params, ctx, prep, &sp, o, want_grad),
        None => Ok(StepOut::default()),
    }
}

/// [`step`] with its host part computed ahead. On a CUDA device ([`on_gpu`]) the network passes
/// and the loss run through [`super::gpu`]; call it from one thread (the dispatching thread).
pub fn step_planned(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, sp: &SlicePlan, o: &StepOptions, want_grad: bool) -> Result<StepOut, PatError> {
    let mut out = StepOut::default();
    let plan = &sp.plan;
    let sel_net = &sp.sel_net;
    let gpu = on_gpu(ctx);
    out.pairs = plan.pairs();
    out.plan_dropped = plan.dropped;
    out.max_dropped_mass = plan.max_dropped_mass;
    out.seconds[0] = sp.seconds;
    let vp: Vec<[f64; 3]> = prep.valid.iter().map(|&r| prep.fit_p[r as usize]).collect();

    let t1 = std::time::Instant::now();
    let (a_net, da_net) = if gpu {
        gpu_coefficients(params, ctx, prep, sel_net, o)?
    } else {
        coefficients(params, ctx, prep, sel_net, o.feature_order, o.network_chunk.max(1))?
    };
    out.rows_forward = sel_net.len();
    out.seconds[1] = t1.elapsed().as_secs_f64();

    // coefficients of all valid tori in valid-list order (unused tori get any finite value; the
    // plan never reads them)
    let t2 = std::time::Instant::now();
    let pos_in_sel: std::collections::HashMap<u32, usize> = sel_net.iter().enumerate().map(|(s, &r)| (r, s)).collect();
    let mut a_valid = vec![0.0; 6 * vp.len()];
    let sphere = [0.0, 0.0, 0.0, 0.0, -0.5, -0.5];
    for (j, &r) in prep.valid.iter().enumerate() {
        match pos_in_sel.get(&r) {
            Some(&s) => a_valid[6 * j..6 * j + 6].copy_from_slice(&a_net[6 * s..6 * s + 6]),
            None => a_valid[6 * j..6 * j + 6].copy_from_slice(&sphere),
        }
    }
    let vn: Vec<[f64; 3]> = prep.valid.iter().map(|&r| prep.fit_n[r as usize]).collect();
    let cloud = CloudFit { p: &vp, n: &vn };
    let (value, ga) = if gpu { gpu_loss(ctx, &cloud, &a_valid, sp, &o.loss, want_grad)? } else { loss_and_grad(ctx, &cloud, &a_valid, plan, &o.loss, want_grad)? };
    out.value = value;
    out.seconds[2] = t2.elapsed().as_secs_f64();
    let Some(ga) = ga else {
        return Ok(out);
    };

    let t3 = std::time::Instant::now();
    // ∂L/∂y per selected network row
    let mut gy_net = vec![0.0; 6 * sel_net.len()];
    for (j, &r) in prep.valid.iter().enumerate() {
        if let Some(&s) = pos_in_sel.get(&r) {
            for c in 0..6 {
                gy_net[6 * s + c] = ga[6 * j + c] * da_net[s][c];
            }
        }
    }
    let keep: Vec<usize> = if o.reference_padding {
        (0..sel_net.len()).collect()
    } else {
        (0..sel_net.len()).filter(|&s| gy_net[6 * s..6 * s + 6].iter().any(|&g| g != 0.0)).collect()
    };
    let rows: Vec<u32> = keep.iter().map(|&s| sel_net[s]).collect();
    let bad_gy = gy_net.iter().filter(|g| !g.is_finite()).count();
    let gy: Vec<f64> = keep.iter().flat_map(|&s| gy_net[6 * s..6 * s + 6].to_vec()).map(|g| if g.is_finite() { g } else { 0.0 }).collect();
    out.rows_backward = rows.len();
    if rows.is_empty() {
        out.grads = Some(params.items.iter().map(|p| p.value.scale(0.0)).collect());
        out.seconds[3] = t3.elapsed().as_secs_f64();
        return Ok(out);
    }
    let mut grads = if gpu { gpu_backprop(params, ctx, prep, &rows, &gy, o)? } else { backprop(params, ctx, prep, &rows, &gy, o)? };
    // non-finite entries
    let mut bad = bad_gy;
    for g in grads.iter_mut() {
        let v = g.to_f64();
        let nb = v.iter().filter(|x| !x.is_finite()).count();
        if nb > 0 {
            bad += nb;
            if o.reference_padding {
                *g = g.like(v.into_iter().map(|x| if x.is_finite() { x } else { 0.0 }).collect(), ctx);
            }
        }
    }
    out.nonfinite_grads = bad;
    out.grads = Some(grads);
    out.seconds[3] = t3.elapsed().as_secs_f64();
    Ok(out)
}

/// The loss of a shape as the reference's validation computes it: the network on all rows, the
/// loss per batch of `batch` rows that has valid rows; returns the per-batch values.
pub fn evaluate_shape(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, batch: usize, o: &StepOptions) -> Result<Vec<LossValue>, PatError> {
    let mut out = Vec::new();
    let mut lo = 0;
    while lo < prep.rows {
        let hi = lo + batch;
        if lo < prep.n {
            out.push(step(params, ctx, prep, lo, hi, o, false)?.value);
        }
        lo = hi;
    }
    Ok(out)
}

#[cfg(feature = "cuda")]
fn gpu_coefficients(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, sel: &[u32], o: &StepOptions) -> Result<(Vec<f64>, Vec<[f64; 6]>), PatError> {
    super::gpu::coefficients(params, ctx, prep, sel, o)
}

#[cfg(feature = "cuda")]
fn gpu_loss(ctx: &Ctx, cloud: &CloudFit<'_>, a: &[f64], sp: &SlicePlan, o: &LossOptions, want_grad: bool) -> Result<(LossValue, Option<Vec<f64>>), PatError> {
    match &sp.gpu {
        Some(l) => super::gpu::loss_and_grad(ctx, cloud, a, &sp.plan, l, o, want_grad),
        None => super::gpu::loss_and_grad(ctx, cloud, a, &sp.plan, &super::gpu::pair_layout(&sp.plan, cloud.p.len(), o), o, want_grad),
    }
}

#[cfg(feature = "cuda")]
fn gpu_backprop(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, rows: &[u32], gy: &[f64], o: &StepOptions) -> Result<Vec<PT>, PatError> {
    super::gpu::backprop(params, ctx, prep, rows, gy, o)
}

#[cfg(not(feature = "cuda"))]
fn gpu_coefficients(_: &NetParams, _: &Ctx, _: &PreparedShape, _: &[u32], _: &StepOptions) -> Result<(Vec<f64>, Vec<[f64; 6]>), PatError> {
    Err(PatError::Gpu { reason: "built without the feature cuda".into() })
}

#[cfg(not(feature = "cuda"))]
fn gpu_loss(_: &Ctx, _: &CloudFit<'_>, _: &[f64], _: &SlicePlan, _: &LossOptions, _: bool) -> Result<(LossValue, Option<Vec<f64>>), PatError> {
    Err(PatError::Gpu { reason: "built without the feature cuda".into() })
}

#[cfg(not(feature = "cuda"))]
fn gpu_backprop(_: &NetParams, _: &Ctx, _: &PreparedShape, _: &[u32], _: &[f64], _: &StepOptions) -> Result<Vec<PT>, PatError> {
    Err(PatError::Gpu { reason: "built without the feature cuda".into() })
}
