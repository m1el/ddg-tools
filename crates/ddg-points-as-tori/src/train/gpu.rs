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
use super::net::{NetParams, PT};
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
    /// Queries (padded to [`ladder`]; the padding has `ŵ = c = 0` and is masked out of the loss).
    pub qn: usize,
    /// Real queries (the first `qr` rows).
    pub qr: usize,
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

/// The smallest size of a fixed ladder (1…8, then multiples of 8 growing by ~1.25×) that holds `x`:
/// chunk shapes are rounded up to it so that the loss's shapes repeat from step to step (burn-fusion
/// re-explores every new shape: 0.29 s per step in phase 0 against 0.065 s with repeated shapes).
pub fn ladder(x: usize) -> usize {
    if x <= 8 {
        return x.max(1);
    }
    let mut v = 8.0f64;
    loop {
        v *= 1.25;
        let r = (v / 8.0).ceil() as usize * 8;
        if r >= x {
            return r;
        }
    }
}

/// Build the [`PairLayout`] of `plan` over `n` tori (host, f32 values of the plan's constants; runs
/// on the current rayon pool). Query chunks are those of [`super::loss::chunks`] with
/// [`LossOptions::pair_chunk`] padded pairs at most; torus chunks likewise.
pub fn pair_layout(plan: &BlendPlan, n: usize, o: &LossOptions) -> PairLayout {
    use rayon::prelude::*;
    let budget = o.pair_chunk.max(1);
    let groups = super::loss::chunks(plan, budget);
    let cnt_q = |j: u32| (plan.offset[j as usize + 1] - plan.offset[j as usize]) as usize;
    let kks: Vec<usize> = groups.iter().map(|qs| ladder(qs.iter().map(|&j| cnt_q(j)).max().unwrap_or(1).max(1))).collect();
    let mut bases = Vec::with_capacity(groups.len());
    let mut base = 0usize;
    for (qs, &kk) in groups.iter().zip(&kks) {
        bases.push(base);
        base += ladder(qs.len()) * kk;
    }
    let m_pad = base;
    let chunks: Vec<DenseChunk> = groups
        .par_iter()
        .zip(kks.par_iter().zip(bases.par_iter()))
        .map(|(qs, (&kk, &base))| {
            let qr = qs.len();
            let qn = ladder(qr);
            let mut d = DenseChunk {
                qn,
                qr,
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
            // padded queries: torus 0, no weight, truth 0 (masked out of the loss)
            d.q.resize(3 * qn, 0.0);
            d.truth.resize(qn, 0.0);
            d.idx.resize(qn * kk, 0);
            d.w.resize(qn * kk, 0.0);
            d.c.resize(3 * qn * kk, 0.0);
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
    for &t in &torder {
        if !cur.is_empty() && (cur.len() + 1) * cnt(t) > budget {
            tgroups.push(std::mem::take(&mut cur));
        }
        cur.push(t);
    }
    if !cur.is_empty() {
        tgroups.push(cur);
    }
    // rows of the concatenated (padded) torus chunks
    let mut row0 = 0usize;
    for ts in &tgroups {
        for (r, &t) in ts.iter().enumerate() {
            inv[t] = (row0 + r) as i32;
        }
        row0 += ladder(ts.len());
    }
    let tchunks: Vec<TorusChunk> = tgroups
        .par_iter()
        .map(|ts| {
            let kk = ladder(ts.iter().map(|&t| cnt(t)).max().unwrap_or(1));
            let tn = ladder(ts.len());
            let mut pos = Vec::with_capacity(tn * kk);
            for &t in ts {
                let l = &flat[count[t]..count[t + 1]];
                pos.extend_from_slice(l);
                pos.extend(std::iter::repeat_n(m_pad as i32, kk - l.len()));
            }
            pos.resize(tn * kk, m_pad as i32);
            TorusChunk { tn, kk, pos }
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
                // the value of this sum is never read; only its gradient is used; padded queries
                // (rows qr..qn) are masked out
                let valid = f2(ctx, &(0..d.qn).map(|r| if r < d.qr { 1.0 } else { 0.0 }).collect::<Vec<f32>>(), [d.qn, 1]);
                let s = ((err.clone() + eik.clone()) * valid).sum_dim(0);
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
    let mut row = 0usize;
    for d in &lay.chunks {
        for t in ev[2 * row..2 * (row + d.qr)].chunks_exact(2) {
            if !t[0].is_finite() || !t[1].is_finite() {
                val.nonfinite += 1;
            }
            sd += t[0];
            se += t[1];
        }
        row += d.qn;
    }
    val.distance = sd * inv;
    val.eikonal = se * inv;
    val.loss = val.distance + val.eikonal;
    Ok((val, grad))
}

/// Row groups of the training forward ([`forward_split`]); `PAT_ROW_GROUPS` overrides it
/// (measurements).
pub const ROW_GROUPS: usize = 64;

fn row_groups() -> usize {
    std::env::var("PAT_ROW_GROUPS").ok().and_then(|v| v.parse().ok()).filter(|&g| g > 0).unwrap_or(ROW_GROUPS)
}

/// A per-column parameter `[c]` broadcast over `[g, rows, c]` through `[g, 1, c]`: its gradient
/// is then reduced in two short stages (over the rows of each group by the broadcast, over the
/// groups by `expand`'s backward), each with many parallel units. Reducing the `[16384, c]`
/// gradient over its leading axis in one go, as burn does for a bias or a LayerNorm scale, runs
/// CubeCL's unit strategy on `c/4` threads that each walk all rows (0.24–0.9 ms per parameter,
/// measured; NOTES-WP10.md).
fn bcast(p: Tensor<1>, g: usize) -> Tensor<3> {
    let c = p.dims()[0];
    p.reshape([1, 1, c]).expand([g, 1, c])
}

fn groups_for(n: usize, groups: usize) -> usize {
    if groups > 1 && n.is_multiple_of(groups) { groups } else { 1 }
}

/// `x W + b` for `x` `[n, din]` as a matmul broadcast over `g` row groups
/// (`[g, n/g, din] · [1, din, dout]`), so that the weight gradient is `g` partial products summed
/// by a short reduction (a split-K by hand; burn's folded `[din, n]·[n, dout]` weight gradient runs
/// on a few SMs at these shapes), and the bias through [`bcast`].
fn lin_split(x: Tensor<2>, w: &Tensor<2>, bias: &Tensor<1>, groups: usize) -> Tensor<2> {
    let [n, din] = x.dims();
    let dout = w.dims()[1];
    let g = groups_for(n, groups);
    let y = x.reshape([g, n / g, din]).matmul(w.clone().unsqueeze::<3>());
    (y + bcast(bias.clone(), g)).reshape([n, dout])
}

/// burn's LayerNorm (the backend's default `layer_norm`: mean, centered, two-pass variance,
/// `√(var + ε)`, scale, shift; the same operations) on `x` `[n, e]` viewed as `[g, n/g, e]`, with
/// the scale and shift through [`bcast`].
fn ln_split(x: Tensor<2>, gamma: &Tensor<1>, beta: &Tensor<1>, groups: usize) -> Tensor<2> {
    let [n, e] = x.dims();
    let g = groups_for(n, groups);
    let x = x.reshape([g, n / g, e]);
    let mean = x.clone().mean_dim(2);
    let c = x - mean;
    let var = (c.clone() * c.clone()).mean_dim(2);
    let y = c / (var + crate::network::LAYER_NORM_EPS).sqrt() * bcast(gamma.clone(), g) + bcast(beta.clone(), g);
    y.reshape([n, e])
}

/// The parameters as plain tensors in the canonical order of [`NetParams`] (input, input norm,
/// blocks: norm1, q, k, v, out, norm2, ffn1, ffn2; final norm, head; kernel before bias, scale
/// before shift).
pub struct FlatParams<'a> {
    d2: Vec<&'a Tensor<2>>,
    d1: Vec<&'a Tensor<1>>,
}

impl<'a> FlatParams<'a> {
    /// From parameter tensors in canonical order (`8 + 16·layers` of them).
    pub fn new(items: &'a [PT]) -> Result<FlatParams<'a>, PatError> {
        let (mut d2, mut d1) = (Vec::new(), Vec::new());
        for p in items {
            match p {
                PT::D2(t) => d2.push(t),
                PT::D1(t) => d1.push(t),
            }
        }
        // per block 6 kernels and 10 vectors; input, head: 1 kernel 1 bias; norms: 2 vectors
        if d2.len() < 2 || (d2.len() - 2) % 6 != 0 || d1.len() != 2 + 4 + 10 * ((d2.len() - 2) / 6) {
            return Err(PatError::Weights { reason: format!("{} kernels and {} vectors are not a network", d2.len(), d1.len()) });
        }
        Ok(FlatParams { d2, d1 })
    }

    fn layers(&self) -> usize {
        (self.d2.len() - 2) / 6
    }
}

/// The network of [`crate::network::burn::PatNet::forward`] (same parameters, same operations up to
/// round-off: pre-LN blocks, attention scores divided by `√d_k` after the product, padded keys
/// filled with burn's `−1e4`, softmax, GELU-tanh) with the linears and norms of [`lin_split`] and
/// [`ln_split`] on `[groups, rows, ·]` views, whose parameter gradients reduce fast on the GPU.
/// Raw outputs `[B, 6]` of features `[B, T, 6]`. The pad mask, if given, fills padded keys; an
/// all-false mask gives the same values as none.
///
/// [formula: pat-network]
pub fn forward_flat(p: &FlatParams<'_>, heads: usize, x: Tensor<3>, mask_pad: Option<Tensor<2, burn::tensor::Bool>>, groups: usize) -> Tensor<2> {
    use burn::tensor::activation::{gelu_approximate, softmax};
    let [b, t, c] = x.dims();
    let n = b * t;
    let g = groups;
    let (k2, k1) = (&p.d2, &p.d1);
    let h0 = lin_split(x.reshape([n, c]), k2[0], k1[0], g);
    let e = h0.dims()[1];
    let d = e / heads;
    let mut h = ln_split(h0, k1[1], k1[2], g);
    let mask = mask_pad.map(|m| m.reshape([b, 1, 1, t]));
    for l in 0..p.layers() {
        let (w, v) = (&k2[1 + 6 * l..7 + 6 * l], &k1[3 + 10 * l..13 + 10 * l]);
        // v: norm1 scale, shift, q b, k b, v b, out b, norm2 scale, shift, ffn1 b, ffn2 b
        let a = ln_split(h.clone(), v[0], v[1], g);
        let split = |y: Tensor<2>| y.reshape([b, t, heads, d]).swap_dims(1, 2);
        let q = split(lin_split(a.clone(), w[0], v[2], g));
        let kk = split(lin_split(a.clone(), w[1], v[3], g));
        let vv = split(lin_split(a, w[2], v[4], g));
        let mut sc = q.matmul(kk.transpose()).div_scalar((d as f32).sqrt());
        if let Some(m) = &mask {
            sc = sc.mask_fill(m.clone(), -1.0e4);
        }
        let ctx = softmax(sc, 3).matmul(vv).swap_dims(1, 2).reshape([n, e]);
        h = h + lin_split(ctx, w[3], v[5], g);
        let f = ln_split(h.clone(), v[6], v[7], g);
        h = h + lin_split(gelu_approximate(lin_split(f, w[4], v[8], g)), w[5], v[9], g);
    }
    let nl = p.layers();
    let h = ln_split(h, k1[3 + 10 * nl], k1[4 + 10 * nl], g);
    let center = h.reshape([b, t, e]).narrow(1, 0, 1).reshape([b, e]);
    lin_split(center, k2[1 + 6 * nl], k1[5 + 10 * nl], 1)
}

/// [`forward_flat`] on the parameters of a built [`crate::network::burn::PatNet`] (tests).
pub fn forward_split(
    net: &crate::network::burn::PatNet,
    heads: usize,
    x: Tensor<3>,
    mask_pad: Option<Tensor<2, burn::tensor::Bool>>,
    groups: usize,
) -> Tensor<2> {
    let lin = |l: &burn::nn::Linear| -> [PT; 2] {
        [
            PT::D2(l.weight.val()),
            PT::D1(l.bias.as_ref().map(|b| b.val()).unwrap_or_else(|| l.weight.val().narrow(0, 0, 1).reshape([l.weight.dims()[1]]).zeros_like())),
        ]
    };
    let norm = |n: &burn::nn::LayerNorm| -> [PT; 2] {
        [PT::D1(n.gamma.val()), PT::D1(n.beta.as_ref().map(|b| b.val()).unwrap_or_else(|| n.gamma.val().zeros_like()))]
    };
    let mut items: Vec<PT> = Vec::new();
    items.extend(lin(&net.input));
    items.extend(norm(&net.input_norm));
    for blk in &net.blocks {
        items.extend(norm(&blk.norm1));
        for l in [&blk.attention.query, &blk.attention.key, &blk.attention.value, &blk.attention.output] {
            items.extend(lin(l));
        }
        items.extend(norm(&blk.norm2));
        items.extend(lin(&blk.ffn1));
        items.extend(lin(&blk.ffn2));
    }
    items.extend(norm(&net.final_norm));
    items.extend(lin(&net.head));
    match FlatParams::new(&items) {
        Ok(fp) => forward_flat(&fp, heads, x, mask_pad, groups),
        Err(_) => net.forward(x, None),
    }
}

/// Network inputs of rows `[lo, hi)` (features and, if any row is short, the pad mask) on `ctx`.
fn row_inputs(ctx: &Ctx, f: &super::net::RowFeatures, lo: usize, hi: usize) -> (Tensor<3>, Option<Tensor<2, burn::tensor::Bool>>) {
    let (k, b) = (f.k, hi - lo);
    let x = ctx.t3(f.feats[lo * k * 6..hi * k * 6].to_vec(), [b, k, 6]);
    let mask = f.pad.as_ref().map(|p| Tensor::<2, burn::tensor::Bool>::from_data(TensorData::new(p[lo * k..hi * k].to_vec(), [b, k]), &ctx.dev));
    (x, mask)
}

// ------------------------------------------------------------------------------------------------
// Ground truth: the winding-number sums on the GPU
// ------------------------------------------------------------------------------------------------

/// The solid-angle sum of `super::truth::MeshTruth::solid_angle_sum` (f64), one thread per query
/// over the faces in their order, with the CPU's operation order (`sub`, `dot` left to right,
/// `cross`, Van Oosterom–Strackee) and no FMA contraction: only `atan2` may differ from the
/// platform libm, by ulps.
const GWN_SRC: &str = r#"
#define TILE 128
extern "C" __global__ void gwn_sum(const double* __restrict__ q, int nq, const double* __restrict__ tri, int nf, double* __restrict__ out) {
    __shared__ double s[TILE * 9];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    double qx = 0.0, qy = 0.0, qz = 0.0;
    if (i < nq) { qx = q[3 * i]; qy = q[3 * i + 1]; qz = q[3 * i + 2]; }
    double acc = 0.0;
    for (int base = 0; base < nf; base += TILE) {
        int m = nf - base < TILE ? nf - base : TILE;
        for (int t = threadIdx.x; t < m * 9; t += blockDim.x) s[t] = tri[base * 9 + t];
        __syncthreads();
        if (i < nq) {
            for (int k = 0; k < m; k++) {
                const double* T = s + 9 * k;
                double ax = T[0] - qx, ay = T[1] - qy, az = T[2] - qz;
                double bx = T[3] - qx, by = T[4] - qy, bz = T[5] - qz;
                double cx = T[6] - qx, cy = T[7] - qy, cz = T[8] - qz;
                double la = sqrt(ax * ax + ay * ay + az * az);
                double lb = sqrt(bx * bx + by * by + bz * bz);
                double lc = sqrt(cx * cx + cy * cy + cz * cz);
                double c0 = by * cz - bz * cy, c1 = bz * cx - bx * cz, c2 = bx * cy - by * cx;
                double num = ax * c0 + ay * c1 + az * c2;
                double ab = ax * bx + ay * by + az * bz;
                double bc = bx * cx + by * cy + bz * cz;
                double ca = cx * ax + cy * ay + cz * az;
                double den = la * lb * lc + ab * lc + bc * la + ca * lb;
                acc += 2.0 * atan2(num, den);
            }
        }
        __syncthreads();
    }
    if (i < nq) out[i] = acc;
}
"#;

struct GwnGpu {
    _ctx: std::sync::Arc<cudarc::driver::CudaContext>,
    stream: std::sync::Arc<cudarc::driver::CudaStream>,
    f: cudarc::driver::CudaFunction,
}

static GWN: std::sync::Mutex<Option<GwnGpu>> = std::sync::Mutex::new(None);

fn gwn_open(ordinal: usize) -> Result<GwnGpu, PatError> {
    use cudarc::nvrtc::{CompileOptions, compile_ptx_with_opts};
    let info = crate::gpu::probe(ordinal)?;
    let e = |w: &'static str| move |err: cudarc::driver::DriverError| PatError::Gpu { reason: format!("{w}: {err:?}") };
    let ctx = cudarc::driver::CudaContext::new(ordinal).map_err(e("context"))?;
    let (major, minor) = info.compute_capability;
    let opts = CompileOptions {
        fmad: Some(false),
        prec_div: Some(true),
        prec_sqrt: Some(true),
        ftz: Some(false),
        options: vec![format!("--gpu-architecture=sm_{major}{minor}")],
        name: Some("ddg_points_as_tori_gwn.cu".into()),
        ..Default::default()
    };
    let ptx = std::panic::catch_unwind(|| compile_ptx_with_opts(GWN_SRC, opts))
        .map_err(|_| PatError::Gpu { reason: "NVRTC (libnvrtc.so.12) could not be loaded".into() })?
        .map_err(|err| PatError::Gpu { reason: format!("NVRTC: {err:?}") })?;
    let module = ctx.load_module(ptx).map_err(e("module"))?;
    let f = module.load_function("gwn_sum").map_err(e("function"))?;
    // its own stream: no implicit synchronization with burn's streams
    let stream = ctx.new_stream().map_err(e("stream"))?;
    Ok(GwnGpu { _ctx: ctx, stream, f })
}

/// Solid-angle sums at `q` over the triangles `tri` (9 coordinates each, in the mesh's face order)
/// on CUDA device 0 (the context is opened on first use and shared by the calling threads).
pub fn solid_angle_sums(q: &[[f64; 3]], tri: &[[f64; 9]]) -> Result<Vec<f64>, PatError> {
    use cudarc::driver::{LaunchConfig, PushKernelArg};
    if q.is_empty() {
        return Ok(Vec::new());
    }
    let mut guard = GWN.lock().map_err(|_| PatError::Backend { reason: "poisoned GPU ground-truth lock".into() })?;
    if guard.is_none() {
        *guard = Some(gwn_open(0)?);
    }
    let g = guard.as_ref().ok_or_else(|| PatError::Backend { reason: "GPU ground truth".into() })?;
    let e = |w: &'static str| move |err: cudarc::driver::DriverError| PatError::Gpu { reason: format!("ground truth {w}: {err:?}") };
    let qf: Vec<f64> = q.iter().flatten().copied().collect();
    let tf: Vec<f64> = if tri.is_empty() { vec![0.0; 9] } else { tri.iter().flatten().copied().collect() };
    let dq = g.stream.clone_htod(&qf).map_err(e("upload"))?;
    let dt = g.stream.clone_htod(&tf).map_err(e("upload"))?;
    let mut dout = g.stream.alloc_zeros::<f64>(q.len()).map_err(e("alloc"))?;
    let (nq, nf) = (q.len() as i32, tri.len() as i32);
    let block = 128u32;
    let cfg = LaunchConfig { grid_dim: ((q.len() as u32).div_ceil(block), 1, 1), block_dim: (block, 1, 1), shared_mem_bytes: 0 };
    // Safety: the buffers have the sizes the kernel reads and writes (3·nq, 9·nf, nq).
    unsafe { g.stream.launch_builder(&g.f).arg(&dq).arg(&nq).arg(&dt).arg(&nf).arg(&mut dout).launch(cfg) }.map_err(e("launch"))?;
    g.stream.clone_dtoh(&dout).map_err(e("download"))
}

// ------------------------------------------------------------------------------------------------
// CUDA graphs of the network batches
// ------------------------------------------------------------------------------------------------

type Shared<T> = std::rc::Rc<std::cell::RefCell<T>>;
type BoxFn<T> = Box<dyn FnMut() -> T>;

/// A captured network batch and its stable input buffers.
struct Slot<T> {
    graph: burn::tensor::Graph<T, BoxFn<T>>,
    x: Shared<Tensor<3>>,
    w: Option<Shared<Tensor<2>>>,
    mask: Option<Shared<Tensor<2, burn::tensor::Bool>>>,
}

/// The graphs of one device and network (per thread: GPU steps run on one dispatching thread).
struct GraphState {
    key: String,
    params: Shared<Vec<PT>>,
    fwd: std::collections::HashMap<bool, Slot<Tensor<2>>>,
    bwd: std::collections::HashMap<bool, Slot<Vec<PT>>>,
}

thread_local! {
    static GRAPHS: std::cell::RefCell<Option<GraphState>> = const { std::cell::RefCell::new(None) };
}

/// Overwrite `dst` in place (its buffer is the one a graph reads) with `src` of the same shape.
fn write_in_place<const D: usize>(dst: &mut Tensor<D>, src: Tensor<D>) {
    let dims = dst.dims();
    let ranges: [std::ops::Range<usize>; D] = std::array::from_fn(|i| 0..dims[i]);
    dst.inplace(|t| t.slice_assign(ranges, src));
}

fn write_mask(dst: &mut Tensor<2, burn::tensor::Bool>, src: Tensor<2, burn::tensor::Bool>) {
    let [a, b] = dst.dims();
    dst.inplace(|t| t.slice_assign([0..a, 0..b], src));
}

/// A uniquely owned device copy of `p` (a stable buffer).
fn fresh(ctx: &Ctx, p: &PT) -> PT {
    let v = p.to_f32();
    match p {
        PT::D1(t) => PT::D1(Tensor::<1>::from_data(TensorData::new(v, t.dims()), (&ctx.dev, ctx.dtype))),
        PT::D2(t) => PT::D2(Tensor::<2>::from_data(TensorData::new(v, t.dims()), (&ctx.dev, ctx.dtype))),
    }
}

fn graph_key(params: &NetParams, ctx: &Ctx, o: &StepOptions, k: usize) -> String {
    format!("{:?}|{:?}|{}|{}|{}|{}", ctx.dev, params.cfg, o.network_chunk, o.forward_chunk, row_groups(), k)
}

/// Run `f` with the thread's graph state for `params` (created on first use, recreated when the
/// device, network or batch shapes change), with the stable parameter buffers set to `params`.
fn with_graphs<R>(params: &NetParams, ctx: &Ctx, o: &StepOptions, k: usize, f: impl FnOnce(&mut GraphState) -> R) -> R {
    let key = graph_key(params, ctx, o, k);
    GRAPHS.with(|cell| {
        let mut st = cell.borrow_mut();
        if st.as_ref().is_none_or(|s| s.key != key) {
            *st = Some(GraphState {
                key,
                params: std::rc::Rc::new(std::cell::RefCell::new(params.items.iter().map(|p| fresh(ctx, &p.value)).collect())),
                fwd: Default::default(),
                bwd: Default::default(),
            });
        } else if let Some(s) = st.as_mut() {
            let mut ps = s.params.borrow_mut();
            for (dst, src) in ps.iter_mut().zip(&params.items) {
                match (dst, &src.value) {
                    (PT::D1(d), PT::D1(v)) => write_in_place(d, v.clone()),
                    (PT::D2(d), PT::D2(v)) => write_in_place(d, v.clone()),
                    _ => {}
                }
            }
        }
        f(st.as_mut().expect("graph state was just set"))
    })
}

fn batch_inputs(ctx: &Ctx, f: &super::net::RowFeatures, lo: usize, hi: usize) -> (Tensor<3>, Tensor<2, burn::tensor::Bool>) {
    let (x, mask) = row_inputs(ctx, f, lo, hi);
    let mask = mask.unwrap_or_else(|| Tensor::<2, burn::tensor::Bool>::from_data(TensorData::new(vec![false; (hi - lo) * f.k], [hi - lo, f.k]), &ctx.dev));
    (x, mask)
}

/// Raw outputs of `nb` forward batches through a captured graph.
fn coefficients_graph(params: &NetParams, ctx: &Ctx, sub: &super::net::RowFeatures, batch: usize, nb: usize, o: &StepOptions) -> Result<Vec<f64>, PatError> {
    let heads = params.cfg.heads;
    let has_mask = sub.pad.is_some();
    with_graphs(params, ctx, o, sub.k, |st| {
        let mut raw = Vec::with_capacity(6 * batch * nb);
        for b in 0..nb {
            let (x, mask) = batch_inputs(ctx, sub, b * batch, (b + 1) * batch);
            if !st.fwd.contains_key(&has_mask) {
                let (xs, ms) =
                    (std::rc::Rc::new(std::cell::RefCell::new(x.clone())), has_mask.then(|| std::rc::Rc::new(std::cell::RefCell::new(mask.clone()))));
                let (p, xc, mc) = (st.params.clone(), xs.clone(), ms.clone());
                let body: BoxFn<Tensor<2>> = Box::new(move || {
                    let ps = p.borrow();
                    let fp = FlatParams::new(&ps).expect("checked parameter layout");
                    forward_flat(&fp, heads, xc.borrow().clone(), mc.as_ref().map(|m| m.borrow().clone()), row_groups())
                });
                let graph = burn::tensor::capture(&ctx.dev, body);
                st.fwd.insert(has_mask, Slot { graph, x: xs, w: None, mask: ms });
            } else if let Some(slot) = st.fwd.get(&has_mask) {
                write_in_place(&mut slot.x.borrow_mut(), x);
                if let Some(m) = &slot.mask {
                    write_mask(&mut m.borrow_mut(), mask);
                }
            }
            let slot = st.fwd.get_mut(&has_mask).ok_or_else(|| PatError::Backend { reason: "graph slot".into() })?;
            // Safety: the graph's buffers (parameters, inputs, intermediates, output) are owned by
            // the slot and the thread's state, written only on this thread, before the replay.
            let y = unsafe { slot.graph.replay() };
            raw.extend(to_f64(y.clone()));
        }
        FlatParams::new(&st.params.borrow()).map(|_| raw)
    })
}

/// Pass 4 through a captured graph: the batch gradients, summed in batch order.
fn backprop_graph(
    params: &NetParams,
    ctx: &Ctx,
    sub: &super::net::RowFeatures,
    w: &[f64],
    batch: usize,
    nb: usize,
    o: &StepOptions,
) -> Result<Vec<PT>, PatError> {
    let heads = params.cfg.heads;
    let has_mask = sub.pad.is_some();
    with_graphs(params, ctx, o, sub.k, |st| {
        let mut acc: Option<Vec<PT>> = None;
        for b in 0..nb {
            let (lo, hi) = (b * batch, (b + 1) * batch);
            let (x, mask) = batch_inputs(ctx, sub, lo, hi);
            let wt = ctx.t2(w[6 * lo..6 * hi].to_vec(), [batch, 6]);
            if !st.bwd.contains_key(&has_mask) {
                let xs = std::rc::Rc::new(std::cell::RefCell::new(x.clone()));
                let ws = std::rc::Rc::new(std::cell::RefCell::new(wt.clone()));
                let ms = has_mask.then(|| std::rc::Rc::new(std::cell::RefCell::new(mask.clone())));
                let (p, xc, wc, mc) = (st.params.clone(), xs.clone(), ws.clone(), ms.clone());
                let body: BoxFn<Vec<PT>> = Box::new(move || {
                    let leaves: Vec<PT> = p.borrow().iter().map(|q| q.map(|t| t.autodiff().require_grad(), |t| t.autodiff().require_grad())).collect();
                    let fp = FlatParams::new(&leaves).expect("checked parameter layout");
                    let y = forward_flat(&fp, heads, xc.borrow().clone(), mc.as_ref().map(|m| m.borrow().clone()), row_groups());
                    let g = (y * wc.borrow().clone()).sum_dim(1).sum_dim(0).backward();
                    leaves
                        .iter()
                        .map(|l| match l {
                            PT::D1(t) => PT::D1(t.grad(&g).unwrap_or_else(|| t.clone().detach().without_autodiff().zeros_like())),
                            PT::D2(t) => PT::D2(t.grad(&g).unwrap_or_else(|| t.clone().detach().without_autodiff().zeros_like())),
                        })
                        .collect()
                });
                let graph = burn::tensor::capture(&ctx.dev, body);
                st.bwd.insert(has_mask, Slot { graph, x: xs, w: Some(ws), mask: ms });
            } else if let Some(slot) = st.bwd.get(&has_mask) {
                write_in_place(&mut slot.x.borrow_mut(), x);
                if let Some(ws) = &slot.w {
                    write_in_place(&mut ws.borrow_mut(), wt);
                }
                if let Some(m) = &slot.mask {
                    write_mask(&mut m.borrow_mut(), mask);
                }
            }
            let slot = st.bwd.get_mut(&has_mask).ok_or_else(|| PatError::Backend { reason: "graph slot".into() })?;
            // Safety: as in `coefficients_graph`; the outputs are consumed (copied into `acc`) before
            // the next replay, which drains the queued copies first.
            let g = unsafe { slot.graph.replay() };
            acc = Some(match acc {
                None => g.iter().map(|t| t.scale(1.0)).collect(),
                Some(a) => a.iter().zip(g.iter()).map(|(x, y)| x.add(y)).collect(),
            });
        }
        acc.ok_or_else(|| PatError::Backend { reason: "no rows to back-propagate".into() })
    })
}

/// `sel` padded to a multiple of `batch` with copies of its first row, and the batch count.
fn padded(sel: &[u32], batch: usize) -> (Vec<u32>, usize) {
    let nb = sel.len().div_ceil(batch);
    let mut v = sel.to_vec();
    v.resize(nb * batch, sel.first().copied().unwrap_or(0));
    (v, nb)
}

/// Pass 2 on the GPU: coefficients and `∂a/∂y` of network rows `sel`, batches of the fixed shape
/// `o.forward_chunk` (or `o.network_chunk`) from the calling thread.
pub fn coefficients(params: &NetParams, ctx: &Ctx, prep: &PreparedShape, sel: &[u32], o: &StepOptions) -> Result<(Vec<f64>, Vec<[f64; 6]>), PatError> {
    if sel.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let batch = if o.forward_chunk > 0 { o.forward_chunk } else { o.network_chunk.max(1) };
    let (rows, nb) = padded(sel, batch);
    let sub = select_rows(&prep.feats, &rows);
    let raw = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| -> Result<Vec<f64>, PatError> {
        if o.gpu_graphs {
            return coefficients_graph(params, ctx, &sub, batch, nb, o);
        }
        let vals: Vec<PT> = params.items.iter().map(|p| p.value.clone()).collect();
        let fp = FlatParams::new(&vals)?;
        let mut raw = Vec::with_capacity(6 * rows.len());
        for b in 0..nb {
            let (x, mask) = row_inputs(ctx, &sub, b * batch, (b + 1) * batch);
            raw.extend(to_f64(forward_flat(&fp, params.cfg.heads, x, mask, row_groups())));
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
        if o.gpu_graphs {
            return backprop_graph(params, ctx, &sub, &w, batch, nb, o);
        }
        let leaves: Vec<PT> = params.items.iter().map(|p| p.value.map(|t| t.autodiff().require_grad(), |t| t.autodiff().require_grad())).collect();
        let fp = FlatParams::new(&leaves)?;
        let mut acc: Option<Vec<PT>> = None;
        for b in 0..nb {
            let (lo, hi) = (b * batch, (b + 1) * batch);
            let (x, mask) = row_inputs(ctx, &sub, lo, hi);
            let y = forward_flat(&fp, params.cfg.heads, x, mask, row_groups());
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

#[cfg(test)]
#[allow(clippy::type_complexity, clippy::single_range_in_vec_init)]
mod tests {
    use super::*;
    use crate::fit::{EigenRule, SignRule};
    use crate::train::loss::build_plan;

    fn lcg(s: &mut u64) -> f64 {
        *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (*s >> 11) as f64 / (1u64 << 53) as f64
    }

    /// An ellipsoid-like cloud with outward normals, random coefficients near the tangent sphere,
    /// and queries around it.
    /// Points, normals, coefficients, queries, ground truth.
    type Synthetic = (Vec<[f64; 3]>, Vec<[f64; 3]>, Vec<f64>, Vec<[f64; 3]>, Vec<f64>);

    fn synthetic(n: usize, nq: usize, seed: u64) -> Synthetic {
        let mut s = seed;
        let mut p = Vec::new();
        let mut nr = Vec::new();
        for _ in 0..n {
            let (u, v) = (2.0 * std::f64::consts::PI * lcg(&mut s), (2.0 * lcg(&mut s) - 1.0).acos());
            let d = [v.sin() * u.cos(), v.sin() * u.sin(), v.cos()];
            p.push([0.6 * d[0], 0.45 * d[1], 0.3 * d[2]]);
            let g = [d[0] / 0.6, d[1] / 0.45, d[2] / 0.3];
            let l = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            nr.push([g[0] / l, g[1] / l, g[2] / l]);
        }
        let a: Vec<f64> = (0..6 * n).map(|i| if i % 6 >= 4 { -0.5 + 0.3 * (lcg(&mut s) - 0.5) } else { 0.1 * (lcg(&mut s) - 0.5) }).collect();
        let q: Vec<[f64; 3]> = (0..nq).map(|_| [1.6 * lcg(&mut s) - 0.8, 1.2 * lcg(&mut s) - 0.6, 0.8 * lcg(&mut s) - 0.4]).collect();
        let truth: Vec<f64> = q.iter().map(|x| (x[0] * x[0] / 0.36 + x[1] * x[1] / 0.2025 + x[2] * x[2] / 0.09).sqrt() - 1.0).collect();
        (p, nr, a, q, truth)
    }

    #[test]
    fn pair_layout_transposes_the_plan() {
        let (p, _, _, q, truth) = synthetic(500, 900, 3);
        for budget in [1usize, 700, 1 << 20] {
            let o = LossOptions { pair_chunk: budget, ..Default::default() };
            let plan = build_plan(&q, &truth, &p, &o, None).unwrap();
            let l = pair_layout(&plan, p.len(), &o);
            assert_eq!(l.order.len(), plan.len());
            // every real pair appears once, at its torus, and the padded ones carry zero weight
            let mut seen = vec![Vec::new(); p.len()];
            for d in &l.chunks {
                assert!(d.qr * d.kk <= 2 * budget.max(d.kk) && d.qn == ladder(d.qr));
                for (k, (&t, &w)) in d.idx.iter().zip(&d.w).enumerate() {
                    if w != 0.0 {
                        seen[t as usize].push((d.base + k) as i32);
                    }
                }
            }
            let mut back = vec![Vec::new(); p.len()];
            let mut r0 = 0usize;
            let mut rows: Vec<Vec<i32>> = Vec::new();
            for tc in &l.tchunks {
                for r in 0..tc.tn {
                    rows.push(tc.pos[r * tc.kk..(r + 1) * tc.kk].iter().copied().filter(|&x| x != l.m_pad as i32).collect());
                }
                r0 += tc.tn;
            }
            assert!(r0 >= p.len());
            for t in 0..p.len() {
                back[t] = rows[l.inv[t] as usize].clone();
            }
            assert_eq!(seen, back, "budget {budget}");
            assert_eq!(seen.iter().map(Vec::len).sum::<usize>(), plan.pairs() - plan.w.iter().filter(|&&w| (w as f32) == 0.0).count());
        }
    }

    #[test]
    fn gpu_loss_matches_cpu_and_repeats() {
        let Ok(gctx) = cuda_ctx(0) else {
            assert!(std::env::var_os("DDG_REQUIRE_GPU").is_none(), "DDG_REQUIRE_GPU is set but no GPU");
            eprintln!("SKIP: no usable CUDA device");
            return;
        };
        let cctx = Ctx { dev: crate::network::burn::cpu_device(), dtype: DType::F32 };
        let (p, nr, a, q, truth) = synthetic(800, 3000, 7);
        // the f32 reference eigenvector form cancels near v2 = 0 (NOTES-WP7.md), so CPU and GPU f32
        // differ there by more than round-off of a well-conditioned fit; the robust form is the
        // training default
        for (sign, eigen, tol) in [
            (SignRule::Reference, EigenRule::Robust, 1e-4),
            (SignRule::Geometric, EigenRule::Robust, 1e-4),
            (SignRule::Reference, EigenRule::Reference, 0.1),
            (SignRule::Geometric, EigenRule::Reference, 0.1),
        ] {
            let mut o = LossOptions::default();
            o.fit.sign_rule = sign;
            o.fit.eigen_rule = eigen;
            let plan = build_plan(&q, &truth, &p, &o, None).unwrap();
            let cloud = CloudFit { p: &p, n: &nr };
            let (vc, gc) = super::super::loss::loss_and_grad(&cctx, &cloud, &a, &plan, &o, true).unwrap();
            let mut og = o;
            og.pair_chunk = 1 << 14;
            let lay = pair_layout(&plan, p.len(), &og);
            let (vg, gg) = loss_and_grad(&gctx, &cloud, &a, &plan, &lay, &og, true).unwrap();
            let (vg2, gg2) = loss_and_grad(&gctx, &cloud, &a, &plan, &lay, &og, true).unwrap();
            let (gc, gg, gg2) = (gc.unwrap(), gg.unwrap(), gg2.unwrap());
            assert!(vg.loss.to_bits() == vg2.loss.to_bits() && gg.iter().zip(&gg2).all(|(x, y)| x.to_bits() == y.to_bits()), "GPU loss not repeatable");
            let gmax = gc.iter().fold(0f64, |m, x| m.max(x.abs()));
            let dmax = gc.iter().zip(&gg).fold(0f64, |m, (x, y)| m.max((x - y).abs()));
            let rl = (vc.loss - vg.loss).abs() / vc.loss;
            eprintln!(
                "{sign:?}/{eigen:?}: loss CPU {} GPU {} (rel {rl:.2e}), ∂L/∂a max|Δ|/max|g| {:.2e}, {} pairs",
                vc.loss,
                vg.loss,
                dmax / gmax,
                plan.pairs()
            );
            assert!(rl < 1e-5 && dmax / gmax < tol, "GPU loss differs from the CPU");
            assert_eq!(vc.nonfinite, vg.nonfinite);
        }
    }

    /// Forward + backward time of the network's parts at the training batch shape (measurement).
    #[test]
    #[ignore]
    fn profile_network_parts() {
        use burn::nn::LayerNormConfig;
        use burn::tensor::activation::{gelu_approximate, softmax};
        let Ok(ctx) = cuda_ctx(0) else { return };
        let b: usize = std::env::var("PAT_B").ok().and_then(|v| v.parse().ok()).unwrap_or(256);
        let (t, e, h, m) = (64usize, 128usize, 8usize, 512usize);
        let d = e / h;
        let mut s = 1u64;
        let mut rnd = |n: usize| -> Vec<f32> { (0..n).map(|_| (lcg(&mut s) - 0.5) as f32).collect() };
        let x0 = Tensor::<3>::from_data(TensorData::new(rnd(b * t * e), [b, t, e]), &ctx.dev);
        let w_e = Tensor::<2>::from_data(TensorData::new(rnd(e * e), [e, e]), &ctx.dev);
        let w_m = Tensor::<2>::from_data(TensorData::new(rnd(e * m), [e, m]), &ctx.dev);
        let time = |name: &str, f: &dyn Fn(Tensor<3>) -> Tensor<3>| {
            for rep in 0..3 {
                let t0 = std::time::Instant::now();
                for _ in 0..8 {
                    let x = x0.clone().autodiff().require_grad();
                    let y = f(x.clone());
                    let g = y.sum_dim(2).sum_dim(1).sum_dim(0).backward();
                    let _ = x.grad(&g).map(|g| g.slice([0..1, 0..1, 0..1]).into_data());
                }
                if rep == 2 {
                    eprintln!("{name}: {:.2} ms per batch of {b} (fwd + bwd)", t0.elapsed().as_secs_f64() * 1e3 / 8.0);
                }
            }
        };
        // reductions over the leading rows (bias and LayerNorm parameter gradients)
        for cols in [128usize, 512] {
            let g = Tensor::<2>::from_data(TensorData::new(rnd(b * t * cols), [b * t, cols]), &ctx.dev);
            let g3 = g.clone().reshape([b, t, cols]);
            let variants: Vec<(&str, Box<dyn Fn() -> Tensor<2>>)> = vec![
                ("sum_dim(0)", Box::new(|| g.clone().sum_dim(0))),
                ("[B,T,C] sum_dim(0).sum_dim(1)", Box::new(|| g3.clone().sum_dim(0).sum_dim(1).reshape([1, cols]))),
                ("swap_dims then sum_dim(1)", Box::new(|| g.clone().swap_dims(0, 1).sum_dim(1).reshape([1, cols]))),
                ("contiguous transpose then sum_dim(1)", Box::new(|| (g.clone().swap_dims(0, 1) * 1.0).sum_dim(1).reshape([1, cols]))),
                (
                    "[128, N/128, C] sum_dim(1) then sum_dim(0)",
                    Box::new(|| g.clone().reshape([128, b * t / 128, cols]).sum_dim(1).sum_dim(0).reshape([1, cols])),
                ),
            ];
            for (name, f) in &variants {
                let _ = f().into_data();
                let t0 = std::time::Instant::now();
                for _ in 0..20 {
                    let _ = f().into_data();
                }
                eprintln!("reduce [{}, {cols}] {name}: {:.3} ms", b * t, t0.elapsed().as_secs_f64() * 1e3 / 20.0);
            }
        }
        // the whole network, parameters tracked (as pass 4) and not
        let params = NetParams::from_model(&crate::Model::reference().unwrap(), &ctx);
        for track in [false, true] {
            let (net, leaves) = params.build(&ctx, track).unwrap();
            let feats = Tensor::<3>::from_data(TensorData::new(rnd(b * t * 6), [b, t, 6]), &ctx.dev);
            for rep in 0..3 {
                let t0 = std::time::Instant::now();
                for _ in 0..8 {
                    let x = feats.clone().autodiff().require_grad();
                    let y = if std::env::var_os("PAT_SPLIT").is_some() {
                        forward_split(&net, 8, x.clone(), None, ROW_GROUPS)
                    } else {
                        net.forward(x.clone(), None)
                    };
                    let g = y.sum_dim(1).sum_dim(0).backward();
                    let gx = x.grad(&g).map(|g| g.slice([0..1, 0..1, 0..1]).into_data());
                    if track && let Some(PT::D2(t)) = leaves.first() {
                        let _ = t.grad(&g).map(|g| g.slice([0..1, 0..1]).into_data());
                    }
                    let _ = gx;
                }
                if rep == 2 {
                    eprintln!("network (params tracked: {track}): {:.2} ms per batch of {b} (fwd + bwd)", t0.elapsed().as_secs_f64() * 1e3 / 8.0);
                }
            }
            if track {
                continue;
            }
            let t0 = std::time::Instant::now();
            for _ in 0..8 {
                let y = net.forward(feats.clone(), None);
                let _ = y.slice([0..1, 0..1]).into_data();
            }
            eprintln!("network forward only: {:.2} ms per batch of {b}", t0.elapsed().as_secs_f64() * 1e3 / 8.0);
        }
        // bias adds with tracked biases, 8 in a row, several formulations
        for (name, mode) in [("x + b[1,E]", 0), ("grouped [32, N/32, E] + b", 1), ("x + b, b untracked", 2), ("(x + b) then * 1.0001", 3)] {
            let bs: Vec<Tensor<1>> = (0..8)
                .map(|_| {
                    let b0 = Tensor::<1>::from_data(TensorData::new(rnd(e), [e]), &ctx.dev);
                    if mode == 2 { b0 } else { b0.autodiff().require_grad() }
                })
                .collect();
            for rep in 0..3 {
                let t0 = std::time::Instant::now();
                for _ in 0..8 {
                    let x = x0.clone().reshape([b * t, e]).autodiff().require_grad();
                    let mut y = x.clone();
                    for bb in &bs {
                        y = match mode {
                            1 => (y.reshape([32, b * t / 32, e]) + bb.clone().unsqueeze::<3>()).reshape([b * t, e]),
                            3 => (y + bb.clone().unsqueeze::<2>()).mul_scalar(1.0001),
                            _ => y + bb.clone().unsqueeze::<2>(),
                        };
                    }
                    let g = y.sum_dim(1).sum_dim(0).backward();
                    let _ = x.grad(&g).map(|g| g.slice([0..1, 0..1]).into_data());
                }
                if rep == 2 {
                    eprintln!("8 bias adds, {name}: {:.2} ms per batch of {b} (fwd + bwd)", t0.elapsed().as_secs_f64() * 1e3 / 8.0);
                }
            }
        }
        // the network with one category of parameters tracked
        for cat in ["linear weights", "linear biases", "norm params", "attention q/k/v/o", "ffn", "head+input"] {
            let (mut net, _) = params.build(&ctx, false).unwrap();
            let tr2 = |p: &mut burn::module::Param<Tensor<2>>| *p = burn::module::Param::from_tensor(p.val().autodiff().require_grad());
            let tr1 = |p: &mut burn::module::Param<Tensor<1>>| *p = burn::module::Param::from_tensor(p.val().autodiff().require_grad());
            let lin = |l: &mut burn::nn::Linear, w: bool, bias: bool| {
                if w {
                    tr2(&mut l.weight);
                }
                if bias && let Some(b) = l.bias.as_mut() {
                    tr1(b);
                }
            };
            let (w, bi) = (cat == "linear weights", cat == "linear biases");
            for blk in &mut net.blocks {
                let att = cat == "attention q/k/v/o";
                let ffn = cat == "ffn";
                for l in [&mut blk.attention.query, &mut blk.attention.key, &mut blk.attention.value, &mut blk.attention.output] {
                    lin(l, w || att, bi || att);
                }
                for l in [&mut blk.ffn1, &mut blk.ffn2] {
                    lin(l, w || ffn, bi || ffn);
                }
                if cat == "norm params" {
                    for n in [&mut blk.norm1, &mut blk.norm2] {
                        tr1(&mut n.gamma);
                        if let Some(b) = n.beta.as_mut() {
                            tr1(b);
                        }
                    }
                }
            }
            if cat == "head+input" {
                lin(&mut net.head, true, true);
                lin(&mut net.input, true, true);
            }
            let feats = Tensor::<3>::from_data(TensorData::new(rnd(b * t * 6), [b, t, 6]), &ctx.dev);
            for rep in 0..3 {
                let t0 = std::time::Instant::now();
                for _ in 0..8 {
                    let x = feats.clone().autodiff().require_grad();
                    let y = net.forward(x.clone(), None);
                    let g = y.sum_dim(1).sum_dim(0).backward();
                    let _ = x.grad(&g).map(|g| g.slice([0..1, 0..1, 0..1]).into_data());
                }
                if rep == 2 {
                    eprintln!("network, tracked {cat}: {:.2} ms per batch of {b} (fwd + bwd)", t0.elapsed().as_secs_f64() * 1e3 / 8.0);
                }
            }
        }
        let mut lnt = LayerNormConfig::new(e).init(&ctx.dev);
        lnt.gamma = burn::module::Param::from_tensor(lnt.gamma.val().autodiff().require_grad());
        lnt.beta = lnt.beta.map(|b| burn::module::Param::from_tensor(b.val().autodiff().require_grad()));
        time("LayerNorm x17, gamma/beta tracked", &|x| {
            let mut x = x;
            for _ in 0..17 {
                x = lnt.forward(x);
            }
            x
        });
        let ln = LayerNormConfig::new(e).init(&ctx.dev);
        time("LayerNorm x17", &|x| {
            let mut x = x;
            for _ in 0..17 {
                x = ln.forward(x);
            }
            x
        });
        time("linear E->E x8", &|x| {
            let mut x = x;
            for _ in 0..8 {
                x = x.matmul(w_e.clone().unsqueeze::<3>());
            }
            x
        });
        time("ffn E->M gelu M->E x8", &|x| {
            let mut x = x;
            for _ in 0..8 {
                x = gelu_approximate(x.matmul(w_m.clone().unsqueeze::<3>())).matmul(w_m.clone().transpose().unsqueeze::<3>());
            }
            x
        });
        time("gelu on [B,T,M] x8", &|x| {
            let mut y = x.matmul(w_m.clone().unsqueeze::<3>());
            for _ in 0..8 {
                y = gelu_approximate(y);
            }
            y.matmul(w_m.clone().transpose().unsqueeze::<3>())
        });
        time("attention core x8", &|x| {
            let mut x = x;
            for _ in 0..8 {
                let q = x.clone().reshape([b, t, h, d]).swap_dims(1, 2);
                let sc = q.clone().matmul(q.clone().transpose()).div_scalar((d as f32).sqrt());
                let w = softmax(sc, 3);
                x = w.matmul(q).swap_dims(1, 2).reshape([b, t, e]);
            }
            x
        });
    }

    /// [`forward_split`] equals [`crate::network::burn::PatNet::forward`] (f64 on the CPU), outputs
    /// and parameter gradients, with and without a pad mask.
    #[test]
    fn forward_split_equals_patnet() {
        let ctx = Ctx { dev: crate::network::burn::cpu_device(), dtype: DType::F64 };
        let cfg = crate::ModelConfig { embed: 16, layers: 2, heads: 2, mlp: 32, k: 8, n_freq: 0 };
        let mut params = NetParams::init(&cfg, 5, &ctx);
        // a non-trivial head (the init's head kernel is zero) and LayerNorm parameters
        let mut s = 11u64;
        params = params
            .with_values(params.items.iter().map(|p| p.value.like(p.value.to_f64().iter().map(|v| v + 0.3 * (lcg(&mut s) - 0.5)).collect(), &ctx)).collect());
        let (b, t) = (8usize, 8usize);
        let feats: Vec<f64> = (0..b * t * 6).map(|_| lcg(&mut s) - 0.5).collect();
        let pad: Vec<bool> = (0..b * t).map(|i| i % t >= t - (i / t) % 3).collect();
        for mask in [false, true] {
            let mk = || mask.then(|| Tensor::<2, burn::tensor::Bool>::from_data(TensorData::new(pad.clone(), [b, t]), &ctx.dev));
            let run = |split: bool| -> (Vec<f64>, Vec<f64>) {
                let (net, leaves) = params.build(&ctx, true).unwrap();
                let x = ctx.t3(feats.clone(), [b, t, 6]);
                let y = if split { forward_split(&net, cfg.heads, x, mk(), 4) } else { net.forward(x, mk()) };
                let w = ctx.t2((0..b * 6).map(|i| (i % 7) as f64 - 3.0).collect(), [b, 6]);
                let yv = to_f64(y.clone());
                let g = (y * w).sum_dim(1).sum_dim(0).backward();
                let gv = leaves
                    .iter()
                    .flat_map(|l| match l {
                        PT::D1(t) => to_f64(t.grad(&g).unwrap()),
                        PT::D2(t) => to_f64(t.grad(&g).unwrap()),
                    })
                    .collect();
                (yv, gv)
            };
            let (y0, g0) = run(false);
            let (y1, g1) = run(true);
            let dy = y0.iter().zip(&y1).fold(0f64, |m, (a, b)| m.max((a - b).abs()));
            let gmax = g0.iter().fold(0f64, |m, a| m.max(a.abs()));
            let dg = g0.iter().zip(&g1).fold(0f64, |m, (a, b)| m.max((a - b).abs()));
            eprintln!("mask {mask}: outputs max |Δ| {dy:.2e}, gradients max |Δ|/max|g| {:.2e}", dg / gmax);
            assert!(dy < 1e-12 && dg / gmax < 1e-12, "forward_split differs from PatNet");
        }
    }

    /// CUDA-graph capture of one network batch (forward with tape, backward, gradients) on stable
    /// buffers refreshed in place: replays equal eager runs bitwise (measurement and check).
    #[test]
    #[ignore]
    fn graph_replay_matches_eager() {
        use std::cell::RefCell;
        let Ok(ctx) = cuda_ctx(0) else { return };
        let b: usize = std::env::var("PAT_B").ok().and_then(|v| v.parse().ok()).unwrap_or(256);
        let t = 64usize;
        let model = crate::Model::reference().unwrap();
        let host = NetParams::from_model(&model, &Ctx { dev: crate::network::burn::cpu_device(), dtype: DType::F32 });
        let upload = |p: &PT, scale: f64| -> PT {
            let v: Vec<f32> = p.to_f32().iter().map(|x| (*x as f64 * scale) as f32).collect();
            match p {
                PT::D1(x) => PT::D1(Tensor::<1>::from_data(TensorData::new(v, x.dims()), &ctx.dev)),
                PT::D2(x) => PT::D2(Tensor::<2>::from_data(TensorData::new(v, x.dims()), &ctx.dev)),
            }
        };
        let mut s = 3u64;
        let mut rnd = |n: usize| -> Vec<f32> { (0..n).map(|_| (lcg(&mut s) - 0.5) as f32).collect() };
        let params = RefCell::new(host.items.iter().map(|p| upload(&p.value, 1.0)).collect::<Vec<PT>>());
        let x = RefCell::new(Tensor::<3>::from_data(TensorData::new(rnd(b * t * 6), [b, t, 6]), &ctx.dev));
        let w = RefCell::new(Tensor::<2>::from_data(TensorData::new(rnd(b * 6), [b, 6]), &ctx.dev));
        let mask = RefCell::new(Tensor::<2, burn::tensor::Bool>::from_data(TensorData::new(vec![false; b * t], [b, t]), &ctx.dev));
        let body = || -> Vec<PT> {
            let leaves: Vec<PT> = params.borrow().iter().map(|p| p.map(|t| t.autodiff().require_grad(), |t| t.autodiff().require_grad())).collect();
            let fp = FlatParams::new(&leaves).unwrap();
            let y = forward_flat(&fp, 8, x.borrow().clone(), Some(mask.borrow().clone()), ROW_GROUPS);
            let g = (y * w.borrow().clone()).sum_dim(1).sum_dim(0).backward();
            leaves
                .iter()
                .map(|l| match l {
                    PT::D1(t) => PT::D1(t.grad(&g).unwrap()),
                    PT::D2(t) => PT::D2(t.grad(&g).unwrap()),
                })
                .collect()
        };
        let hash =
            |g: &[PT]| -> u64 { g.iter().flat_map(|t| t.to_f32()).fold(0xcbf29ce484222325u64, |h, x| (h ^ x.to_bits() as u64).wrapping_mul(0x100000001b3)) };
        let t0 = std::time::Instant::now();
        let mut graph = burn::tensor::capture(&ctx.dev, body);
        eprintln!("capture: {:.2} s, hardware graph: {}", t0.elapsed().as_secs_f64(), graph.is_hardware());
        let refresh3 = |dst: &RefCell<Tensor<3>>, v: Vec<f32>| {
            let d = dst.borrow().dims();
            let src = Tensor::<3>::from_data(TensorData::new(v, d), &ctx.dev);
            dst.borrow_mut().inplace(|t| t.slice_assign([0..d[0], 0..d[1], 0..d[2]], src));
        };
        let refresh2 = |dst: &RefCell<Tensor<2>>, v: Vec<f32>| {
            let d = dst.borrow().dims();
            let src = Tensor::<2>::from_data(TensorData::new(v, d), &ctx.dev);
            dst.borrow_mut().inplace(|t| t.slice_assign([0..d[0], 0..d[1]], src));
        };
        for it in 0..4 {
            refresh3(&x, rnd(b * t * 6));
            refresh2(&w, rnd(b * 6));
            // parameters changed in place too
            {
                let scale = 1.0 + 0.01 * it as f64;
                let mut ps = params.borrow_mut();
                for (dst, src) in ps.iter_mut().zip(&host.items) {
                    let up = upload(&src.value, scale);
                    match (dst, up) {
                        (PT::D1(d), PT::D1(u)) => {
                            let n = d.dims()[0];
                            d.inplace(|t| t.slice_assign([0..n], u));
                        }
                        (PT::D2(d), PT::D2(u)) => {
                            let [r, c] = d.dims();
                            d.inplace(|t| t.slice_assign([0..r, 0..c], u));
                        }
                        _ => unreachable!(),
                    }
                }
            }
            let tr = std::time::Instant::now();
            let hg = hash(unsafe { graph.replay() });
            let t_replay = tr.elapsed().as_secs_f64();
            let te = std::time::Instant::now();
            let he = hash(&body());
            let t_eager = te.elapsed().as_secs_f64();
            eprintln!("iteration {it}: replay {hg:016x} ({:.1} ms), eager {he:016x} ({:.1} ms)", t_replay * 1e3, t_eager * 1e3);
            assert_eq!(hg, he, "replay differs from eager");
        }
    }
}
