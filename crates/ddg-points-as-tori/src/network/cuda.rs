//! GPU inference through burn-cuda (feature `cuda`, SPEC.md §5.2, §5.4, §18.2 D4).
//!
//! The network runs on a CUDA device of burn's dispatch backend (CubeCL kernels compiled by NVRTC
//! at run time). The crate builds burn without `autotune` and without `fusion`: kernel selection
//! is a fixed heuristic of the shape, never a timing, so a fixed shape always runs the same
//! kernels.
//!
//! **TF32.** CubeCL's default f32 matmul on `sm_80+` is a tensor-core (`wmma …tf32`) kernel: the
//! operands are rounded to TF32 (10 explicit mantissa bits) inside the kernel, which moves the raw
//! outputs by ~1e-4 relative against the CPU (NOTES-GPU.md). That is the default,
//! [`GpuMatmul::Tf32`] (user decision 2026-11-04). The option [`GpuMatmul::Split3`] recovers f32
//! accuracy with framework ops only: every operand is split exactly into a TF32 head and an f32
//! tail, `a = a_h + a_l`, and `a·b ≈ a_h·b_h + (a_h·b_l + a_l·b_h)` (the `a_l·b_l`
//! term is ~2⁻²² relative). Weights are split on the host by truncation, activations on the device
//! by Veltkamp's split `c = 8193 a`, `a_h = c − (c − a)` (separate kernels, no FMA contraction).
//!
//! Batches have the fixed shape [`GPU_BATCH`] × tokens (the last batch is padded with copies of
//! its first neighborhood and the padding is dropped). The outputs were measured to be bitwise
//! independent of the batch shape anyway (NOTES-GPU.md); the constant makes it hold by
//! construction.
//!
//! [formula: pat-network]

use super::{GpuMatmul, LinearW, Model, NetWeights, NormW};
use crate::PatError;
use burn::nn::{LayerNorm, LayerNormConfig};
use burn::module::Param;
use burn::tensor::activation::{gelu_approximate, softmax};
use burn::tensor::{Device, Tensor, TensorData};

/// Neighborhoods per GPU batch (a constant: the shape of every launch, NOTES-GPU.md).
pub const GPU_BATCH: usize = 256;

/// A linear layer on the device: `W` as is (TF32 mode), its exact split `W = W_h + W_l`
/// (Split3) and the bias.
struct Lin {
    full: Tensor<3>,
    hi: Tensor<3>,
    lo: Tensor<3>,
    b: Tensor<3>,
}

struct Block {
    norm1: LayerNorm,
    /// q, k and v side by side (`E → 3E`): one matmul for the three projections.
    qkv: Lin,
    o: Lin,
    norm2: LayerNorm,
    ffn1: Lin,
    ffn2: Lin,
}

/// The network's weights on a CUDA device, split for [`GpuMatmul::Split3`].
struct GpuNet {
    input: Lin,
    input_norm: LayerNorm,
    blocks: Vec<Block>,
    final_norm: LayerNorm,
    head: Lin,
    heads: usize,
    mode: GpuMatmul,
}

/// `w = hi + lo` with `hi` the TF32 truncation of `w` (exact: `lo` has at most 13 bits).
fn split_host(w: &[f32]) -> (Vec<f32>, Vec<f32>) {
    let hi: Vec<f32> = w.iter().map(|x| f32::from_bits(x.to_bits() & 0xFFFF_E000)).collect();
    let lo = w.iter().zip(&hi).map(|(x, h)| x - h).collect();
    (hi, lo)
}

/// The layers `ws` (same `din`) side by side along the outputs, uploaded.
fn lin(ws: &[&LinearW], dev: &Device) -> Lin {
    let din = ws[0].din;
    let dout: usize = ws.iter().map(|w| w.dout).sum();
    let mut w = Vec::with_capacity(din * dout);
    for r in 0..din {
        for l in ws {
            w.extend_from_slice(&l.w[r * l.dout..(r + 1) * l.dout]);
        }
    }
    let bias: Vec<f32> = ws.iter().flat_map(|l| l.b.iter().copied()).collect();
    let (hi, lo) = split_host(&w);
    let t = |v: Vec<f32>| Tensor::<3>::from_data(TensorData::new(v, [1, din, dout]), dev);
    Lin { full: t(w), hi: t(hi), lo: t(lo), b: Tensor::<3>::from_data(TensorData::new(bias, [1, 1, dout]), dev) }
}

fn norm(w: &NormW, dev: &Device) -> LayerNorm {
    let e = w.scale.len();
    let mut ln = LayerNormConfig::new(e).with_epsilon(super::LAYER_NORM_EPS).init(dev);
    ln.gamma = Param::from_tensor(Tensor::<1>::from_data(TensorData::new(w.scale.clone(), [e]), dev));
    ln.beta = Some(Param::from_tensor(Tensor::<1>::from_data(TensorData::new(w.bias.clone(), [e]), dev)));
    ln
}

/// Veltkamp's split of an f32 tensor: `a = a_h + a_l` exactly, `a_h` with 11 significant bits.
///
/// `c = 8193 a` is formed as `8192 a + a`: the product by a power of two is exact, so the one
/// rounding is the same whether or not the compiler contracts it into an FMA, and the remaining
/// steps are additions only. The split therefore survives burn's kernel fusion (which compiles the
/// whole chain into one kernel with FMA contraction on).
fn split_dev<const D: usize>(a: Tensor<D>) -> (Tensor<D>, Tensor<D>) {
    let c = a.clone() * 8192.0 + a.clone();
    let hi = c.clone() - (c - a.clone());
    let lo = a - hi.clone();
    (hi, lo)
}

impl GpuNet {
    fn new(w: &NetWeights, heads: usize, mode: GpuMatmul, dev: &Device) -> GpuNet {
        let blocks = w
            .blocks
            .iter()
            .map(|b| Block {
                norm1: norm(&b.norm1, dev),
                qkv: lin(&[&b.q, &b.k, &b.v], dev),
                o: lin(&[&b.o], dev),
                norm2: norm(&b.norm2, dev),
                ffn1: lin(&[&b.ffn1], dev),
                ffn2: lin(&[&b.ffn2], dev),
            })
            .collect();
        GpuNet {
            input: lin(&[&w.input], dev),
            input_norm: norm(&w.input_norm, dev),
            blocks,
            final_norm: norm(&w.final_norm, dev),
            head: lin(&[&w.head], dev),
            heads,
            mode,
        }
    }

    /// `a·b` with `b = b_h + b_l` given split: Split3 `a_h b_h + (a_h b_l + a_l b_h)`.
    fn mm<const D: usize>(&self, a: Tensor<D>, b_hi: Tensor<D>, b_lo: Tensor<D>) -> Tensor<D> {
        let (ah, al) = split_dev(a);
        let small = ah.clone().matmul(b_lo) + al.matmul(b_hi.clone());
        ah.matmul(b_hi) + small
    }

    /// `x W + b`.
    fn linear(&self, x: Tensor<3>, l: &Lin) -> Tensor<3> {
        match self.mode {
            GpuMatmul::Tf32 => x.matmul(l.full.clone()) + l.b.clone(),
            GpuMatmul::Split3 => self.mm(x, l.hi.clone(), l.lo.clone()) + l.b.clone(),
        }
    }

    /// `a·b` for activations.
    fn mm_act(&self, a: Tensor<4>, b: Tensor<4>) -> Tensor<4> {
        match self.mode {
            GpuMatmul::Tf32 => a.matmul(b),
            GpuMatmul::Split3 => {
                let (bh, bl) = split_dev(b);
                self.mm(a, bh, bl)
            }
        }
    }

    /// The attention of burn's `MultiHeadAttention` (no mask), with our matmuls: scores
    /// `q kᵀ / √d`, softmax over keys, context, output projection.
    fn attention(&self, x: Tensor<3>, blk: &Block) -> Tensor<3> {
        let [b, t, e] = x.dims();
        let (h, d) = (self.heads, e / self.heads);
        let qkv = self.linear(x, &blk.qkv);
        let heads = |i: usize| qkv.clone().narrow(2, i * e, e).reshape([b, t, h, d]).swap_dims(1, 2);
        let (q, k, v) = (heads(0), heads(1), heads(2));
        let scores = self.mm_act(q, k.transpose()).div_scalar((d as f32).sqrt());
        let w = softmax(scores, 3);
        let ctx = self.mm_act(w, v).swap_dims(1, 2).reshape([b, t, e]);
        self.linear(ctx, &blk.o)
    }

    /// Features `[B, T, 6]` → raw outputs `[B, 6]`.
    fn forward(&self, features: Tensor<3>) -> Tensor<2> {
        let [b, _, _] = features.dims();
        let mut x = self.input_norm.forward(self.linear(features, &self.input));
        for blk in &self.blocks {
            let a = self.attention(blk.norm1.forward(x.clone()), blk);
            x = x + a;
            let h = blk.norm2.forward(x.clone());
            let f = self.linear(gelu_approximate(self.linear(h, &blk.ffn1)), &blk.ffn2);
            x = x + f;
        }
        let center = self.final_norm.forward(x.narrow(1, 0, 1)); // per token: narrow first
        self.linear(center, &self.head).reshape([b, 6])
    }
}

/// The network on one CUDA device.
pub(crate) struct CudaRunner {
    net: GpuNet,
    dev: Device,
    batch: usize,
    ordinal: usize,
}

impl CudaRunner {
    /// Build the network on device `ordinal` with the default batch and matmul mode.
    pub(crate) fn new(model: &Model, ordinal: usize, mode: GpuMatmul) -> Result<Self, PatError> {
        Self::with(model, ordinal, GPU_BATCH, mode)
    }

    /// As [`CudaRunner::new`] with another batch shape or matmul mode (measurements only).
    pub(crate) fn with(model: &Model, ordinal: usize, batch: usize, mode: GpuMatmul) -> Result<Self, PatError> {
        crate::gpu::probe(ordinal)?;
        let built = std::panic::catch_unwind(|| {
            let dev = Device::cuda(ordinal);
            let net = GpuNet::new(&model.weights, model.config.heads, mode, &dev);
            (net, dev)
        });
        let (net, dev) = built.map_err(|_| PatError::Gpu { reason: format!("burn-cuda could not open CUDA device {ordinal}") })?;
        Ok(CudaRunner { net, dev, batch: batch.max(1), ordinal })
    }

    /// The device ordinal.
    pub(crate) fn ordinal(&self) -> usize {
        self.ordinal
    }

    /// Queue the forward pass of `b ≤ batch` neighborhoods of `t` tokens as one padded batch; the
    /// returned tensor is read back later ([`CudaRunner::read`]).
    fn launch(&self, features: &[f32], b: usize, t: usize) -> Result<Tensor<2>, PatError> {
        let row = t * 6;
        if b == 0 || b > self.batch || features.len() != b * row {
            return Err(PatError::SizeMismatch { what: format!("{} feature values for {b} neighborhoods of {t} tokens", features.len()) });
        }
        let mut data = Vec::with_capacity(self.batch * row);
        data.extend_from_slice(features);
        for _ in b..self.batch {
            data.extend_from_slice(&features[..row]);
        }
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let x = Tensor::<3>::from_data(TensorData::new(data, [self.batch, t, 6]), &self.dev);
            self.net.forward(x)
        }))
        .map_err(|_| PatError::Gpu { reason: "burn-cuda panicked during the forward pass".into() })
    }

    /// Read back one queued batch (`b` real rows).
    fn read(y: Tensor<2>, b: usize) -> Result<Vec<f32>, PatError> {
        let out = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| y.into_data().try_into_vec::<f32>().map_err(|e| format!("{e:?}"))));
        let mut y = match out {
            Ok(Ok(y)) => y,
            Ok(Err(reason)) => return Err(PatError::Backend { reason }),
            Err(_) => return Err(PatError::Gpu { reason: "burn-cuda panicked reading the outputs".into() }),
        };
        y.truncate(b * 6);
        Ok(y)
    }

    /// Raw outputs of `b ≤ batch` neighborhoods of `t` tokens, run as one padded batch.
    pub(crate) fn run_batch(&self, features: &[f32], b: usize, t: usize) -> Result<Vec<f32>, PatError> {
        Self::read(self.launch(features, b, t)?, b)
    }

    /// Raw outputs of any number of neighborhoods, in batches of the runner's batch shape, one
    /// batch at a time with a read-back after each (NOTES-GPU.md: queueing several batches before
    /// reading back, or dispatching from several host threads, gave no gain, and queueing changed
    /// burn's fusion segments and with them the bits of the Split3 outputs).
    pub(crate) fn run(&self, features: &[f32], t: usize) -> Result<Vec<f32>, PatError> {
        let row = t * 6;
        let mut out = Vec::with_capacity(features.len() / row.max(1) * 6);
        for chunk in features.chunks(self.batch * row) {
            out.extend(self.run_batch(chunk, chunk.len() / row, t)?);
        }
        Ok(out)
    }
}

/// Raw outputs of precomputed features (`b·t·6` values) on CUDA device `ordinal`, in batches of
/// `batch` neighborhoods (`None`: [`GPU_BATCH`]) with matmul `mode`. For parity tests and
/// measurements; the pipeline goes through [`super::predict_coefficients`] with
/// [`super::Device::Cuda`].
pub fn forward_features_cuda(
    model: &Model,
    features: &[f32],
    tokens: usize,
    ordinal: usize,
    batch: Option<usize>,
    mode: GpuMatmul,
) -> Result<Vec<f32>, PatError> {
    if tokens == 0 || !features.len().is_multiple_of(6 * tokens) {
        return Err(PatError::SizeMismatch { what: format!("{} feature values for {tokens} tokens", features.len()) });
    }
    let r = CudaRunner::with(model, ordinal, batch.unwrap_or(GPU_BATCH), mode)?;
    r.run(features, tokens)
}

/// Neighborhoods whose features are computed on the CPU while the GPU runs the previous group.
const FEATURE_GROUP: usize = 8 * GPU_BATCH;

/// [`super::predict_coefficients`] on CUDA device `ordinal` (inputs already validated by the
/// caller): features on the CPU pool in groups, overlapped with the forward pass on the GPU, then
/// the σ scaling on the host. The report's `device` is the GPU; `features_seconds` is CPU wall time
/// of the feature groups, `forward_seconds` the wall time of the GPU batches.
pub(crate) fn predict_coefficients_cuda(
    model: &Model,
    cloud: &crate::PointCloud,
    nbhds: &crate::spatial::Neighborhoods,
    opts: &super::PredictOptions,
    ordinal: usize,
) -> Result<super::Prediction, PatError> {
    use crate::features::{FeatureScratch, neighborhood_features};
    use rayon::prelude::*;
    let t0 = std::time::Instant::now();
    let n = cloud.len();
    let t = nbhds.k;
    let pool = crate::pool::make_pool(opts.threads)?;
    let runner = CudaRunner::new(model, ordinal, opts.gpu_matmul)?;
    let (pos, nrm) = (cloud.positions(), cloud.normals());
    let order = opts.feature_order;
    let groups: Vec<(usize, usize)> = (0..n).step_by(FEATURE_GROUP).map(|lo| (lo, (lo + FEATURE_GROUP).min(n))).collect();
    let features = |lo: usize, hi: usize| -> (Vec<f32>, Vec<f64>, f64) {
        let tf = std::time::Instant::now();
        let mut feat = vec![0f32; (hi - lo) * t * 6];
        let mut sig = vec![0f64; hi - lo];
        pool.install(|| {
            feat.par_chunks_mut(t * 6).zip(sig.par_iter_mut()).enumerate().for_each_init(FeatureScratch::default, |scr, (r, (f, s))| {
                *s = neighborhood_features(pos, nrm, nbhds.row(lo + r), order, scr, f);
            })
        });
        (feat, sig, tf.elapsed().as_secs_f64())
    };
    let mut coeffs = Vec::with_capacity(n);
    let mut sigma = Vec::with_capacity(n);
    let mut report = super::PredictReport {
        device: super::Device::Cuda { ordinal: runner.ordinal() },
        tokens: t,
        chunks: n.div_ceil(GPU_BATCH),
        ..Default::default()
    };
    std::thread::scope(|s| -> Result<(), PatError> {
        let (tx, rx) = std::sync::mpsc::sync_channel::<(Vec<f32>, Vec<f64>, f64)>(1);
        let groups = &groups;
        let features = &features;
        s.spawn(move || {
            for &(lo, hi) in groups {
                if tx.send(features(lo, hi)).is_err() {
                    return;
                }
            }
        });
        for _ in 0..groups.len() {
            let (feat, sig, tf) = rx.recv().map_err(|_| PatError::Backend { reason: "the feature thread stopped".into() })?;
            report.features_seconds += tf;
            let tn = std::time::Instant::now();
            let raw = runner.run(&feat, t)?;
            report.forward_seconds += tn.elapsed().as_secs_f64();
            for (r, &sg) in sig.iter().enumerate() {
                let y = [raw[6 * r], raw[6 * r + 1], raw[6 * r + 2], raw[6 * r + 3], raw[6 * r + 4], raw[6 * r + 5]];
                coeffs.push(super::scale_outputs(y, sg, order));
            }
            sigma.extend(sig);
        }
        Ok(())
    })?;
    let bad: Vec<u32> = coeffs.iter().enumerate().filter(|(_, a)| !a.iter().all(|x| x.is_finite())).map(|(i, _)| i as u32).collect();
    if !bad.is_empty() {
        return Err(PatError::NonFiniteOutput { indices: bad });
    }
    report.seconds = t0.elapsed().as_secs_f64();
    Ok((coeffs, sigma, report))
}

/// A network kept on a CUDA device (weights uploaded once) for repeated calls; see
/// [`forward_features_cuda`]. Measurements only.
pub struct CudaNetwork(CudaRunner);

impl CudaNetwork {
    /// Upload `model` to device `ordinal`.
    pub fn new(model: &Model, ordinal: usize, mode: GpuMatmul) -> Result<Self, PatError> {
        Ok(CudaNetwork(CudaRunner::new(model, ordinal, mode)?))
    }

    /// Raw outputs of precomputed features (`b·t·6` values).
    pub fn forward(&self, features: &[f32], tokens: usize) -> Result<Vec<f32>, PatError> {
        if tokens == 0 || !features.len().is_multiple_of(6 * tokens) {
            return Err(PatError::SizeMismatch { what: format!("{} feature values for {tokens} tokens", features.len()) });
        }
        self.0.run(features, tokens)
    }
}
