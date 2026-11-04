//! The network's parameters for training (SPEC.md §6, §10.4): a named list in the reference's Flax
//! order on a device, built into WP1's burn module [`PatNet`] for each pass, initialized as the
//! reference (lecun-normal kernels, zero biases, unit LayerNorm scales, a zero head kernel with
//! bias `(0, 0, 0, 0, −½, −½)`: the tangent sphere), and saved in the crate's safetensors layout
//! (the Flax names and shapes of `weights/pat-reference.safetensors`, so
//! [`crate::Model::from_safetensors`] reads a checkpoint). The reference's `nnx.Rngs(42)` stream is
//! not reproduced (D22).

use super::queries::uniform;
use super::tfit::Ctx;
use crate::network::burn::{PatBlock, PatNet};
use crate::network::{BlockW, LinearW, ModelConfig, NetWeights, NormW};
use crate::{Model, PatError};
use burn::module::Param;
use burn::nn::attention::MultiHeadAttentionConfig;
use burn::nn::{LayerNorm, LayerNormConfig, Linear};
use burn::tensor::{DType, Tensor, TensorData};
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;
use std::collections::HashMap;

/// A parameter tensor of rank 1 or 2.
#[derive(Debug, Clone)]
pub enum PT {
    /// Biases, LayerNorm scales and biases.
    D1(Tensor<1>),
    /// Kernels `[in, out]`.
    D2(Tensor<2>),
}

impl PT {
    /// Values as f64 (row major).
    pub fn to_f64(&self) -> Vec<f64> {
        let d = match self {
            PT::D1(t) => t.clone().into_data(),
            PT::D2(t) => t.clone().into_data(),
        };
        d.convert::<f64>().try_to_vec::<f64>().unwrap_or_default()
    }

    /// Values as f32 (row major).
    pub fn to_f32(&self) -> Vec<f32> {
        let d = match self {
            PT::D1(t) => t.clone().into_data(),
            PT::D2(t) => t.clone().into_data(),
        };
        d.convert::<f32>().try_to_vec::<f32>().unwrap_or_default()
    }

    /// Apply `f` to the tensor whatever its rank.
    pub fn map(&self, f1: impl FnOnce(Tensor<1>) -> Tensor<1>, f2: impl FnOnce(Tensor<2>) -> Tensor<2>) -> PT {
        match self {
            PT::D1(t) => PT::D1(f1(t.clone())),
            PT::D2(t) => PT::D2(f2(t.clone())),
        }
    }

    /// Number of values.
    pub fn len(&self) -> usize {
        match self {
            PT::D1(t) => t.dims()[0],
            PT::D2(t) => t.dims()[0] * t.dims()[1],
        }
    }

    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `self + o` (same rank).
    pub fn add(&self, o: &PT) -> PT {
        match (self, o) {
            (PT::D1(a), PT::D1(b)) => PT::D1(a.clone() + b.clone()),
            (PT::D2(a), PT::D2(b)) => PT::D2(a.clone() + b.clone()),
            _ => self.clone(),
        }
    }

    /// `c·self`.
    pub fn scale(&self, c: f64) -> PT {
        self.map(|t| t.mul_scalar(c), |t| t.mul_scalar(c))
    }

    /// The same shape from f64 values on `ctx`.
    pub fn like(&self, v: Vec<f64>, ctx: &Ctx) -> PT {
        match self {
            PT::D1(t) => PT::D1(Tensor::<1>::from_data(TensorData::new(v, t.dims()), (&ctx.dev, ctx.dtype))),
            PT::D2(t) => PT::D2(Tensor::<2>::from_data(TensorData::new(v, t.dims()), (&ctx.dev, ctx.dtype))),
        }
    }
}

/// One named parameter.
#[derive(Debug, Clone)]
pub struct ParamItem {
    /// Flax name (`blocks.3.attention.query.kernel`, …).
    pub name: String,
    /// Flax shape (as stored in the safetensors file).
    pub flax_shape: Vec<usize>,
    /// The values on the device, as `[in, out]` or `[len]`.
    pub value: PT,
}

/// All parameters, in the canonical order (input, input norm, blocks, final norm, head; per block
/// norm1, q, k, v, out, norm2, ffn1, ffn2; kernel before bias, scale before bias).
#[derive(Debug, Clone)]
pub struct NetParams {
    /// Architecture.
    pub cfg: ModelConfig,
    /// The 8 + 16·layers tensors.
    pub items: Vec<ParamItem>,
}

/// One entry of [`layout`]: name, our shape, Flax shape, host values.
type Entry = (String, Vec<usize>, Vec<usize>, Vec<f32>);

/// Flax names and shapes, in the canonical order, with the host values of `w`.
fn layout(cfg: &ModelConfig, w: &NetWeights) -> Vec<Entry> {
    let (e, h) = (cfg.embed, cfg.heads);
    let d = e / h;
    let mut v: Vec<Entry> = Vec::new();
    fn lin(v: &mut Vec<Entry>, p: &str, l: &LinearW) {
        v.push((format!("{p}.kernel"), vec![l.din, l.dout], vec![l.din, l.dout], l.w.clone()));
        v.push((format!("{p}.bias"), vec![l.dout], vec![l.dout], l.b.clone()));
    }
    fn nrm(v: &mut Vec<Entry>, p: &str, n: &NormW) {
        v.push((format!("{p}.scale"), vec![n.scale.len()], vec![n.scale.len()], n.scale.clone()));
        v.push((format!("{p}.bias"), vec![n.bias.len()], vec![n.bias.len()], n.bias.clone()));
    }
    lin(&mut v, "input_encoder.mlp.layers.0", &w.input);
    nrm(&mut v, "input_encoder.mlp.layers.1", &w.input_norm);
    for (l, b) in w.blocks.iter().enumerate() {
        let p = format!("blocks.{l}");
        nrm(&mut v, &format!("{p}.norm1"), &b.norm1);
        for (nm, x) in [("query", &b.q), ("key", &b.k), ("value", &b.v)] {
            v.push((format!("{p}.attention.{nm}.kernel"), vec![e, e], vec![e, h, d], x.w.clone()));
            v.push((format!("{p}.attention.{nm}.bias"), vec![e], vec![h, d], x.b.clone()));
        }
        v.push((format!("{p}.attention.out.kernel"), vec![e, e], vec![h, d, e], b.o.w.clone()));
        v.push((format!("{p}.attention.out.bias"), vec![e], vec![e], b.o.b.clone()));
        nrm(&mut v, &format!("{p}.norm2"), &b.norm2);
        lin(&mut v, &format!("{p}.ffn.layers.0"), &b.ffn1);
        lin(&mut v, &format!("{p}.ffn.layers.3"), &b.ffn2);
    }
    nrm(&mut v, "final_norm", &w.final_norm);
    lin(&mut v, "mlp.layers.0", &w.head);
    v
}

/// Reads parameters in canonical order while building the module.
struct Cursor<'a> {
    items: std::slice::Iter<'a, ParamItem>,
    track: bool,
    leaves: Vec<PT>,
}

impl Cursor<'_> {
    fn next(&mut self) -> Result<PT, PatError> {
        let v = self.items.next().map(|p| p.value.clone()).ok_or_else(|| PatError::Weights { reason: "too few parameters".into() })?;
        let v = if self.track { v.map(|t| t.autodiff().require_grad(), |t| t.autodiff().require_grad()) } else { v };
        self.leaves.push(v.clone());
        Ok(v)
    }

    fn d1(&mut self) -> Result<Param<Tensor<1>>, PatError> {
        match self.next()? {
            PT::D1(t) => Ok(Param::from_tensor(t)),
            PT::D2(_) => Err(PatError::Weights { reason: "expected a 1-D parameter".into() }),
        }
    }

    fn d2(&mut self) -> Result<Param<Tensor<2>>, PatError> {
        match self.next()? {
            PT::D2(t) => Ok(Param::from_tensor(t)),
            PT::D1(_) => Err(PatError::Weights { reason: "expected a 2-D parameter".into() }),
        }
    }

    fn linear(&mut self) -> Result<Linear, PatError> {
        let weight = self.d2()?;
        Ok(Linear { weight, bias: Some(self.d1()?) })
    }

    fn norm(&mut self, ctx: &Ctx, e: usize) -> Result<LayerNorm, PatError> {
        let mut ln = LayerNormConfig::new(e).with_epsilon(crate::network::LAYER_NORM_EPS).init(&ctx.dev);
        ln.gamma = self.d1()?;
        ln.beta = Some(self.d1()?);
        Ok(ln)
    }
}

impl NetParams {
    /// Upload framework-neutral weights to `ctx` (in its dtype).
    pub fn from_weights(cfg: &ModelConfig, w: &NetWeights, ctx: &Ctx) -> NetParams {
        let items = layout(cfg, w)
            .into_iter()
            .map(|(name, shape, flax_shape, vals)| {
                let value = if shape.len() == 2 {
                    PT::D2(Tensor::<2>::from_data(TensorData::new(vals, [shape[0], shape[1]]), (&ctx.dev, ctx.dtype)))
                } else {
                    PT::D1(Tensor::<1>::from_data(TensorData::new(vals, [shape[0]]), (&ctx.dev, ctx.dtype)))
                };
                ParamItem { name, flax_shape, value }
            })
            .collect();
        NetParams { cfg: *cfg, items }
    }

    /// From a model.
    pub fn from_model(m: &Model, ctx: &Ctx) -> NetParams {
        NetParams::from_weights(m.config(), m.weights(), ctx)
    }

    /// The same parameters cast to `dtype`.
    pub fn cast(&self, dtype: DType) -> NetParams {
        let items = self.items.iter().map(|it| ParamItem { value: it.value.map(|t| t.cast(dtype), |t| t.cast(dtype)), ..it.clone() }).collect();
        NetParams { cfg: self.cfg, items }
    }

    /// Replace the values (same order and shapes).
    pub fn with_values(&self, values: Vec<PT>) -> NetParams {
        let items = self.items.iter().zip(values).map(|(it, value)| ParamItem { value, ..it.clone() }).collect();
        NetParams { cfg: self.cfg, items }
    }

    /// Back to host weights (f32).
    pub fn to_weights(&self) -> NetWeights {
        let mut vals: Vec<Vec<f32>> = self.items.iter().map(|p| p.value.to_f32()).collect();
        vals.reverse();
        let mut pop = || vals.pop().unwrap_or_default();
        let (e, m) = (self.cfg.embed, self.cfg.mlp);
        fn l2(pop: &mut dyn FnMut() -> Vec<f32>, din: usize, dout: usize) -> LinearW {
            let w = pop();
            LinearW { din, dout, w, b: pop() }
        }
        fn n2(pop: &mut dyn FnMut() -> Vec<f32>) -> NormW {
            let scale = pop();
            NormW { scale, bias: pop() }
        }
        let input = l2(&mut pop, 6, e);
        let input_norm = n2(&mut pop);
        let blocks = (0..self.cfg.layers)
            .map(|_| BlockW {
                norm1: n2(&mut pop),
                q: l2(&mut pop, e, e),
                k: l2(&mut pop, e, e),
                v: l2(&mut pop, e, e),
                o: l2(&mut pop, e, e),
                norm2: n2(&mut pop),
                ffn1: l2(&mut pop, e, m),
                ffn2: l2(&mut pop, m, e),
            })
            .collect();
        let final_norm = n2(&mut pop);
        let head = l2(&mut pop, e, 6);
        NetWeights { input, input_norm, blocks, final_norm, head }
    }

    /// As a [`Model`].
    pub fn to_model(&self) -> Model {
        Model::from_weights(self.cfg, self.to_weights())
    }

    /// Number of scalar parameters.
    pub fn count(&self) -> usize {
        self.items.iter().map(|p| p.value.len()).sum()
    }

    /// Build the burn module on `ctx`. With `track`, every parameter is a fresh autodiff leaf that
    /// requires its gradient; the leaves are returned in canonical order.
    pub fn build(&self, ctx: &Ctx, track: bool) -> Result<(PatNet, Vec<PT>), PatError> {
        let (e, h) = (self.cfg.embed, self.cfg.heads);
        let mut c = Cursor { items: self.items.iter(), track, leaves: Vec::with_capacity(self.items.len()) };
        let input = c.linear()?;
        let input_norm = c.norm(ctx, e)?;
        let mut blocks = Vec::with_capacity(self.cfg.layers);
        for _ in 0..self.cfg.layers {
            let norm1 = c.norm(ctx, e)?;
            let mut attention = MultiHeadAttentionConfig::new(e, h).with_dropout(0.0).init(&ctx.dev);
            attention.query = c.linear()?;
            attention.key = c.linear()?;
            attention.value = c.linear()?;
            attention.output = c.linear()?;
            let norm2 = c.norm(ctx, e)?;
            let ffn1 = c.linear()?;
            let ffn2 = c.linear()?;
            blocks.push(PatBlock { norm1, attention, norm2, ffn1, ffn2 });
        }
        let final_norm = c.norm(ctx, e)?;
        let head = c.linear()?;
        Ok((PatNet { input, input_norm, blocks, final_norm, head }, c.leaves))
    }

    /// Fresh parameters as the reference initializes them, seeded (D22): kernels lecun-normal
    /// (truncated normal at ±2 standard deviations, variance `1/fan_in`, Flax's
    /// `variance_scaling(1, fan_in, truncated_normal)`), biases 0, LayerNorm scale 1 and bias 0, head
    /// kernel 0 and bias `(0, 0, 0, 0, −½, −½)`, so that the initial output is the tangent sphere of
    /// radius `σ_i` at every point.
    pub fn init(cfg: &ModelConfig, seed: u64, ctx: &Ctx) -> NetParams {
        let (e, m) = (cfg.embed, cfg.mlp);
        let mut rng = ChaCha8Rng::seed_from_u64(seed ^ 0x494e_4954); // "INIT"
        // truncated standard normal on [−2, 2], rescaled to unit variance (Flax's constant)
        let tn = |rng: &mut ChaCha8Rng| loop {
            let (u1, u2) = (uniform(rng), uniform(rng));
            let z = (-2.0 * (1.0 - u1).ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            if z.abs() <= 2.0 {
                break z / 0.879_625_661_034_239_8;
            }
        };
        let kernel = |din: usize, dout: usize, rng: &mut ChaCha8Rng| {
            let s = (1.0 / din as f64).sqrt();
            LinearW { din, dout, w: (0..din * dout).map(|_| (s * tn(rng)) as f32).collect(), b: vec![0.0; dout] }
        };
        let norm = || NormW { scale: vec![1.0; e], bias: vec![0.0; e] };
        let input = kernel(6, e, &mut rng);
        let blocks = (0..cfg.layers)
            .map(|_| BlockW {
                norm1: norm(),
                q: kernel(e, e, &mut rng),
                k: kernel(e, e, &mut rng),
                v: kernel(e, e, &mut rng),
                o: kernel(e, e, &mut rng),
                norm2: norm(),
                ffn1: kernel(e, m, &mut rng),
                ffn2: kernel(m, e, &mut rng),
            })
            .collect();
        let head = LinearW { din: e, dout: 6, w: vec![0.0; e * 6], b: vec![0.0, 0.0, 0.0, 0.0, -0.5, -0.5] };
        let w = NetWeights { input, input_norm: norm(), blocks, final_norm: norm(), head };
        NetParams::from_weights(cfg, &w, ctx)
    }

    /// Serialize as safetensors (f32, Flax names and shapes) with `metadata` added to the config.
    pub fn to_safetensors(&self, metadata: &HashMap<String, String>) -> Result<Vec<u8>, PatError> {
        let bytes: Vec<(String, Vec<usize>, Vec<u8>)> = self
            .items
            .iter()
            .map(|p| (p.name.clone(), p.flax_shape.clone(), p.value.to_f32().iter().flat_map(|x| x.to_le_bytes()).collect()))
            .collect();
        let views: Vec<(String, safetensors::tensor::TensorView<'_>)> = bytes
            .iter()
            .map(|(n, s, b)| {
                safetensors::tensor::TensorView::new(safetensors::Dtype::F32, s.clone(), b)
                    .map(|v| (n.clone(), v))
                    .map_err(|e| PatError::Weights { reason: format!("{n}: {e}") })
            })
            .collect::<Result<_, _>>()?;
        let mut meta = metadata.clone();
        meta.insert(
            "config".into(),
            format!(
                "{{\"embed_dim\": {}, \"n_layers\": {}, \"n_heads\": {}, \"mlp_dim\": {}, \"k_neighbors\": {}, \"n_freq\": 0}}",
                self.cfg.embed, self.cfg.layers, self.cfg.heads, self.cfg.mlp, self.cfg.k
            ),
        );
        meta.insert("format".into(), "flax-nnx paths, flax shapes, float32 little-endian".into());
        safetensors::serialize(views, Some(meta)).map_err(|e| PatError::Weights { reason: format!("serialize: {e}") })
    }
}

/// Network inputs of a batch of rows: features `[B, k, 6]` (f64 values), scales `σ`, and the pad
/// mask `[B, k]` (true = padding) when any row is short.
#[derive(Debug, Clone, Default)]
pub struct RowFeatures {
    /// Rows.
    pub b: usize,
    /// Tokens per row.
    pub k: usize,
    /// `b·k·6` features.
    pub feats: Vec<f64>,
    /// `σ_i` per row.
    pub sigma: Vec<f64>,
    /// `b·k` pad flags, `None` when no row is short.
    pub pad: Option<Vec<bool>>,
}

/// Run `net` (built on `ctx`) on features: raw outputs `[B, 6]`.
pub fn forward(net: &PatNet, ctx: &Ctx, f: &RowFeatures, lo: usize, hi: usize) -> Tensor<2> {
    let (k, b) = (f.k, hi - lo);
    let x = ctx.t3(f.feats[lo * k * 6..hi * k * 6].to_vec(), [b, k, 6]);
    let mask = f.pad.as_ref().map(|p| Tensor::<2, burn::tensor::Bool>::from_data(TensorData::new(p[lo * k..hi * k].to_vec(), [b, k]), &ctx.dev));
    net.forward(x, mask)
}
