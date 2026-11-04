//! Training on a CUDA device (feature `cuda` + `train`; SPEC.md §10.4, §18 WP7 item 9; NOTES-WP10.md).
//!
//! The step of [`super::step`] with three passes replaced, all through burn-cuda (CubeCL kernels,
//! fusion on, TF32 matmuls: the user's decision of 2026-11-04 for training):
//!
//! - **Coefficients** (pass 2): the network without the tape on the plan's rows, in batches of the
//!   fixed shape [`super::StepOptions::network_chunk`] × `k` tokens (the last batch is padded with
//!   copies of its first row), one batch after the other from the calling thread.
//! - **Loss and `∂L/∂a`** (pass 3): the tori are fitted once; the per-pair gather of the tori from the
//!   detached fit feeds one leaf per dense chunk of pairs, so the backward of the blend stops at the
//!   pairs. burn-cubecl's backward of `select` is `select_assign` with one thread per output column
//!   (9 threads walking all pairs serially; deterministic but ~1 s per step at 3.7M pairs), and the
//!   backward of `gather` is an atomic scatter-add (not repeatable). Instead, the pair gradients are
//!   summed per torus by a **gather** over the transposed pair lists ([`PairLayout`], built on the
//!   host ahead of time) and a fixed-shape `sum_dim`, then pushed through the fit by a second
//!   backward of `⟨T(a), ∂L/∂T⟩`. Exact (the chain rule regrouped) and free of atomics.
//! - **Backprop** (pass 4): the network with the tape on the rows with a nonzero `∂L/∂y`, in batches
//!   of the same fixed shape, accumulating the parameter gradients on the device in batch order.
//!
//! No full `sum()` is used where its value matters: burn-cubecl's `float_sum` is a one-shot atomic
//! kernel (`SumStrategy::OneShot` without `autotune`). Loss values are read back per query and
//! summed on the host in f64, as on the CPU.

use super::loss::{BlendPlan, CloudFit, LossOptions, LossValue};
use super::net::{NetParams, PT, forward};
use super::step::{PreparedShape, StepOptions, scale_raw, select_rows};
use super::tfit::{Ctx, FitGeometry, fit_tensor, to_f64};
use crate::PatError;
use burn::tensor::{DType, Device, Int, Tensor, TensorData};

/// Largest CUDA ordinal [`is_cuda`] recognizes.
const MAX_ORDINAL: usize = 16;

/// True when `dev` (with or without autodiff) is a CUDA device of burn's dispatch backend.
pub fn is_cuda(dev: &Device) -> bool {
    let d = dev.clone().without_autodiff();
    (0..MAX_ORDINAL).any(|i| d == Device::cuda(i))
}

/// The training context on CUDA device `ordinal` (f32). Errors when the driver, the ordinal or
/// NVRTC is missing.
pub fn cuda_ctx(ordinal: usize) -> Result<Ctx, PatError> {
    crate::gpu::probe(ordinal)?;
    let dev = std::panic::catch_unwind(|| {
        let dev = Device::cuda(ordinal);
        // open the device (and compile nothing yet): one tiny upload and read back
        let t = Tensor::<1>::from_data(TensorData::new(vec![1f32], [1]), &dev);
        let _ = t.into_data();
        dev
    })
    .map_err(|_| PatError::Gpu { reason: format!("burn-cuda could not open CUDA device {ordinal}") })?;
    Ok(Ctx { dev, dtype: DType::F32 })
}

fn gpu_err(what: &str) -> impl FnOnce(Box<dyn std::any::Any + Send>) -> PatError + '_ {
    move |e| {
        let msg = e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_default();
        PatError::Gpu { reason: format!("burn-cuda panicked in {what}: {msg}") }
    }
}

/// One dense chunk of the GPU loss: `qn` queries × `kk` pairs.
#[derive(Debug, Clone, Default)]
pub struct DenseChunk {
    /// Queries.
    pub qn: usize,
    /// Pairs per query (padded).
    pub kk: usize,
    /// First padded pair of the chunk in the global padded order.
    pub base: usize,
    /// Query positions `[qn, 3]`.
    pub q: Vec<f32>,
    /// Ground truth `[qn]`.
    pub truth: Vec<f32>,
    /// Torus per pair `[qn·kk]` (padding: the query's first torus, with `ŵ = c = 0`).
    pub idx: Vec<i32>,
    /// `ŵ` per pair.
    pub w: Vec<f32>,
    /// `c` per pair `[qn·kk, 3]`.
    pub c: Vec<f32>,
}

/// One dense chunk of the per-torus reduction: `tn` tori × `kk` pair positions.
#[derive(Debug, Clone, Default)]
pub struct TorusChunk {
    /// Tori.
    pub tn: usize,
    /// Positions per torus (padded with the zero row `m_pad`).
    pub kk: usize,
    /// Global padded pair positions `[tn·kk]`.
    pub pos: Vec<i32>,
}

/// The host layout of the GPU loss for a [`BlendPlan`]: the queries in dense chunks (sorted by their
/// pair count, as [`super::loss::chunks`]) and, transposed, the pair positions of every torus in
/// dense chunks (sorted by count), with the permutation back to torus order.
#[derive(Debug, Clone, Default)]
pub struct PairLayout {
    /// Query chunks.
    pub chunks: Vec<DenseChunk>,
    /// Padded pairs in all.
    pub m_pad: usize,
    /// Tori.
    pub n: usize,
    /// Torus chunks.
    pub tchunks: Vec<TorusChunk>,
    /// `inv[t]` = row of torus `t` in the concatenated torus chunks.
    pub inv: Vec<i32>,
    /// Indices of the plan's queries in chunk order (the loss's summation order).
    pub order: Vec<u32>,
}

/// Build the [`PairLayout`] of `plan` over `n` tori (host, f32 values of the plan's constants; runs
/// on the current rayon pool). Query chunks are those of [`super::loss::chunks`] with
/// [`LossOptions::pair_chunk`] padded pairs at most; torus chunks likewise.
pub fn pair_layout(plan: &BlendPlan, n: usize, o: &LossOptions) -> PairLayout {
    use rayon::prelude::*;
    let budget = o.pair_chunk.max(1);
    let groups = super::loss::chunks(plan, budget);
    let cnt_q = |j: u32| (plan.offset[j as usize + 1] - plan.offset[j as usize]) as usize;
    let kks: Vec<usize> = groups.iter().map(|qs| qs.iter().map(|&j| cnt_q(j)).max().unwrap_or(1).max(1)).collect();
    let mut bases = Vec::with_capacity(groups.len());
    let mut base = 0usize;
    for (qs, &kk) in groups.iter().zip(&kks) {
        bases.push(base);
        base += qs.len() * kk;
    }
    let m_pad = base;
    let chunks: Vec<DenseChunk> = groups
        .par_iter()
        .zip(kks.par_iter().zip(bases.par_iter()))
        .map(|(qs, (&kk, &base))| {
            let qn = qs.len();
            let mut d = DenseChunk {
                qn,
                kk,
                base,
                q: Vec::with_capacity(3 * qn),
                truth: Vec::with_capacity(qn),
                idx: Vec::with_capacity(qn * kk),
                w: Vec::with_capacity(qn * kk),
                c: Vec::with_capacity(3 * qn * kk),
            };
            for &j in qs {
                let j = j as usize;
                d.q.extend(plan.q[j].map(|v| v as f32));
                d.truth.push(plan.truth[j] as f32);
                let (a, b) = (plan.offset[j] as usize, plan.offset[j + 1] as usize);
                let first = if b > a { plan.idx[a] as i32 } else { 0 };
                for k in 0..kk {
                    if a + k < b {
                        d.idx.push(plan.idx[a + k] as i32);
                        d.w.push(plan.w[a + k] as f32);
                        d.c.extend(plan.c[a + k].map(|v| v as f32));
                    } else {
                        d.idx.push(first);
                        d.w.push(0.0);
                        d.c.extend([0f32; 3]);
                    }
                }
            }
            d
        })
        .collect();
    let order: Vec<u32> = groups.iter().flatten().copied().collect();
    // transposed (counting sort): the padded positions of every torus's real pairs, ascending
    let mut count = vec![0usize; n + 1];
    for &t in &plan.idx {
        if (t as usize) < n {
            count[t as usize + 1] += 1;
        }
    }
    for t in 0..n {
        count[t + 1] += count[t];
    }
    let mut flat = vec![0i32; count[n]];
    let mut fill = count.clone();
    for (d, qs) in chunks.iter().zip(&groups) {
        for (r, &j) in qs.iter().enumerate() {
            let a = plan.offset[j as usize] as usize;
            for k in 0..cnt_q(j) {
                let t = plan.idx[a + k] as usize;
                if t < n {
                    flat[fill[t]] = (d.base + r * d.kk + k) as i32;
                    fill[t] += 1;
                }
            }
        }
    }
    let cnt = |t: usize| (count[t + 1] - count[t]).max(1);
    let mut torder: Vec<usize> = (0..n).collect();
    torder.sort_by_key(|&t| (cnt(t), t));
    let mut inv = vec![0i32; n];
    let mut tgroups: Vec<Vec<usize>> = Vec::new();
    let mut cur: Vec<usize> = Vec::new();
    for (row, &t) in torder.iter().enumerate() {
        if !cur.is_empty() && (cur.len() + 1) * cnt(t) > budget {
            tgroups.push(std::mem::take(&mut cur));
        }
        cur.push(t);
        inv[t] = row as i32;
    }
    if !cur.is_empty() {
        tgroups.push(cur);
    }
    let tchunks: Vec<TorusChunk> = tgroups
        .par_iter()
        .map(|ts| {
            let kk = ts.iter().map(|&t| cnt(t)).max().unwrap_or(1);
            let mut pos = Vec::with_capacity(ts.len() * kk);
            for &t in ts {
                let l = &flat[count[t]..count[t + 1]];
                pos.extend_from_slice(l);
                pos.extend(std::iter::repeat_n(m_pad as i32, kk - l.len()));
            }
            TorusChunk { tn: ts.len(), kk, pos }
        })
        .collect();
    PairLayout { chunks, m_pad, n, tchunks, inv, order }
}

fn int1(ctx: &Ctx, v: &[i32]) -> Tensor<1, Int> {
    Tensor::<1, Int>::from_data(TensorData::new(v.to_vec(), [v.len()]), &ctx.dev)
}

fn f2(ctx: &Ctx, v: &[f32], shape: [usize; 2]) -> Tensor<2> {
    Tensor::<2>::from_data(TensorData::new(v.to_vec(), shape), (&ctx.dev, ctx.dtype))
}

/// φ `[Q, 1]` and ∇φ `[Q, 3]` of a dense chunk from its gathered tori `t` (`[Q·K, 9]`): the
/// framework ops of [`super::loss`] (same guards), with the query positions broadcast per pair.
///
/// [formula: pat-torus-sdf]
/// [formula: pat-torus-sdf-gradient]
fn blend_chunk(ctx: &Ctx, t: Tensor<2>, d: &DenseChunk) -> (Tensor<2>, Tensor<2>) {
    let (qn, kk) = (d.qn, d.kk);
    let m = qn * kk;
    let c = t.clone().narrow(1, 0, 3);
    let u = t.clone().narrow(1, 3, 3);
    let big_r = t.clone().narrow(1, 6, 1);
    let r = t.clone().narrow(1, 7, 1);
    let s = t.narrow(1, 8, 1);
    let q = f2(ctx, &d.q, [qn, 3]).reshape([qn, 1, 3]).expand([qn, kk, 3]).reshape([m, 3]);
    let w = f2(ctx, &d.w, [m, 1]);
    let cc = f2(ctx, &d.c, [m, 3]);
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
    let phi = (w.clone() * g.clone()).reshape([qn, kk]).sum_dim(1);
    let gphi = (w * (grad_g + cc * g)).reshape([qn, kk, 3]).sum_dim(1).reshape([qn, 3]);
    (phi, gphi)
}

/// The loss and `∂L/∂a` on the GPU (see the module docs); `a` as in
/// [`super::loss::loss_and_grad`], whose values and gradient this equals up to f32 round-off.
///
/// [formula: pat-loss]
pub fn loss_and_grad(
    ctx: &Ctx,
    cloud: &CloudFit<'_>,
    a: &[f64],
    plan: &BlendPlan,
    lay: &PairLayout,
    o: &LossOptions,
    want_grad: bool,
) -> Result<(LossValue, Option<Vec<f64>>), PatError> {
    let n = cloud.p.len();
    if a.len() != 6 * n || cloud.n.len() != n || lay.n != n || lay.order.len() != plan.len() {
        return Err(PatError::SizeMismatch {
            what: format!("{} coefficients for {n} points, layout of {} tori and {} queries for a plan of {}", a.len(), lay.n, lay.order.len(), plan.len()),
        });
    }
    let nq = plan.len();
    let mut val = LossValue { queries: nq, ..Default::default() };
    if nq == 0 {
        return Ok((val, want_grad.then(|| vec![0.0; 6 * n])));
    }
    let inv = 1.0 / nq as f64;
    let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let geo = FitGeometry::new(ctx, cloud.p, cloud.n);
        let mut leaf = ctx.t2(a.to_vec(), [n, 6]);
        if want_grad {
            leaf = leaf.autodiff().require_grad();
        }
        let tori = fit_tensor(leaf.clone(), &geo, &o.fit).packed;
        let td = tori.clone().detach().without_autodiff();
        let mut leaves = Vec::with_capacity(lay.chunks.len());
        let mut terms = Vec::with_capacity(lay.chunks.len());
        let mut obj: Option<Tensor<2>> = None;
        for d in &lay.chunks {
            let mut tp = td.clone().select(0, int1(ctx, &d.idx));
            if want_grad {
                tp = tp.autodiff().require_grad();
            }
            let (phi, gphi) = blend_chunk(ctx, tp.clone(), d);
            let truth = f2(ctx, &d.truth, [d.qn, 1]);
            let err = (phi - truth).abs();
            let gn = gphi.square().sum_dim(1).clamp_min(1e-30).sqrt();
            let eik = (gn.neg() + 1.0).abs();
            if want_grad {
                // the value of this sum is never read; only its gradient (ones) is used
                let s = (err.clone() + eik.clone()).sum_dim(0);
                obj = Some(match obj {
                    None => s,
                    Some(x) => x + s,
                });
            }
            terms.push(Tensor::cat(vec![err.detach().without_autodiff(), eik.detach().without_autodiff()], 1));
            leaves.push(tp);
        }
        let ev = to_f64(Tensor::cat(terms, 0));
        let grad = match obj {
            Some(obj) if want_grad => {
                let g1 = obj.mul_scalar(inv).backward();
                let mut gp: Vec<Tensor<2>> = leaves
                    .iter()
                    .zip(&lay.chunks)
                    .map(|(l, d)| l.grad(&g1).unwrap_or_else(|| Tensor::<2>::zeros([d.qn * d.kk, 9], (&ctx.dev, ctx.dtype))))
                    .collect();
                gp.push(Tensor::<2>::zeros([1, 9], (&ctx.dev, ctx.dtype)));
                let gp = Tensor::cat(gp, 0);
                let per: Vec<Tensor<2>> =
                    lay.tchunks.iter().map(|tc| gp.clone().select(0, int1(ctx, &tc.pos)).reshape([tc.tn, tc.kk, 9]).sum_dim(1).reshape([tc.tn, 9])).collect();
                let gt = Tensor::cat(per, 0).select(0, int1(ctx, &lay.inv));
                let g2 = (tori * gt).sum_dim(1).sum_dim(0).backward();
                Some(leaf.grad(&g2).map(to_f64).unwrap_or_else(|| vec![0.0; 6 * n]))
            }
            _ => None,
        };
        (ev, grad)
    }));
    let (ev, grad) = run.map_err(gpu_err("the loss"))?;
    let (mut sd, mut se) = (0.0, 0.0);
    for t in ev.chunks_exact(2) {
        if !t[0].is_finite() || !t[1].is_finite() {
            val.nonfinite += 1;
        }
        sd += t[0];
        se += t[1];
    }
    val.distance = sd * inv;
    val.eikonal = se * inv;
    val.loss = val.distance + val.eikonal;
    Ok((val, grad))
}

/// `sel` padded to a multiple of `batch` with copies of its first row, and the batch count.
fn padded(sel: &[u32], batch: usize) -> (Vec<u32>, usize) {
    let nb = sel.len().div_ceil(batch);
    let mut v = sel.to_vec();
    v.resize(nb * batch, sel.first().copied().unwrap_or(0));
    (v, nb)
}

/// Pass 2 on the GPU: coefficients and `∂a/∂y` of network rows `sel`, batches of the fixed shape
/// `o.network_chunk` from the calling thread.
pub fn coefficients(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, sel: &[u32], o: &StepOptions) -> Result<(Vec<f64>, Vec<[f64; 6]>), PatError> {
    if sel.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let batch = o.network_chunk.max(1);
    let (rows, nb) = padded(sel, batch);
    let sub = select_rows(&prep.feats, &rows);
    let raw = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| -> Result<Vec<f64>, PatError> {
        let (net, _) = params.build(ctx, false)?;
        let mut raw = Vec::with_capacity(6 * rows.len());
        for b in 0..nb {
            raw.extend(to_f64(forward(&net, ctx, &sub, b * batch, (b + 1) * batch)));
        }
        Ok(raw)
    }))
    .map_err(gpu_err("the forward pass"))??;
    Ok(scale_raw(&raw[..6 * sel.len()], &sub.sigma[..sel.len()], ctx, o.feature_order))
}

/// Pass 4 on the GPU: parameter gradients of `Σ_rows ⟨g_y, y(θ)⟩` for network rows `sel` with
/// output weights `gy`, batches of the fixed shape `o.network_chunk` (padded rows have weight 0),
/// summed on the device in batch order.
pub fn backprop(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, sel: &[u32], gy: &[f64], o: &StepOptions) -> Result<Vec<PT>, PatError> {
    if sel.is_empty() || gy.len() != 6 * sel.len() {
        return Err(PatError::SizeMismatch { what: format!("{} output weights for {} rows", gy.len(), sel.len()) });
    }
    let batch = o.network_chunk.max(1);
    let (rows, nb) = padded(sel, batch);
    let sub = select_rows(&prep.feats, &rows);
    let mut w = gy.to_vec();
    w.resize(6 * rows.len(), 0.0);
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| -> Result<Vec<PT>, PatError> {
        let (net, leaves) = params.build(ctx, true)?;
        let mut acc: Option<Vec<PT>> = None;
        for b in 0..nb {
            let (lo, hi) = (b * batch, (b + 1) * batch);
            let y = forward(&net, ctx, &sub, lo, hi);
            let wt = ctx.t2(w[6 * lo..6 * hi].to_vec(), [batch, 6]);
            let grads = (y * wt).sum_dim(1).sum_dim(0).backward();
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
        }
        acc.ok_or_else(|| PatError::Backend { reason: "no rows to back-propagate".into() })
    }))
    .map_err(gpu_err("the backward pass"))?
}
