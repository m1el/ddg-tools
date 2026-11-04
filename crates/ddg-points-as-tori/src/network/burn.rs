//! The network as a burn 0.22 module (SPEC.md §5.2, §6.1).
//!
//! Architecture of the reference's `FundamentalFormPredictor` (points-as-tori, MIT, Copyright (c)
//! 2026 Nicole Feng; NOTES-REFERENCE.md §2.3): `LayerNorm(Linear 6→E)`, `L` pre-LN blocks
//! `x += MHA(LN₁ x); x += Linear(GELU_tanh(Linear(LN₂ x)))`, final LayerNorm, the center token,
//! `Linear E→6`.
//!
//! Differences in round-off from Flax, all within the parity tolerance:
//! - burn's LayerNorm uses the two-pass variance; Flax the "fast" `E[x²] − E[x]²`.
//! - burn divides the attention scores by `√d_k = 4` after the product; Flax scales the query
//!   first. Division by 4 is exact, so the two agree up to the product's own rounding.
//!
//! [formula: pat-network]

use super::{LinearW, ModelConfig, NetWeights, NormW};
use burn::module::{Module, Param};
use burn::nn::attention::{MhaInput, MultiHeadAttention, MultiHeadAttentionConfig};
use burn::nn::{LayerNorm, LayerNormConfig, Linear};
use burn::tensor::{Bool, Device, Tensor, TensorData, activation};

/// One pre-LN transformer block.
#[derive(Module, Debug)]
pub struct PatBlock {
    /// LayerNorm before attention.
    pub norm1: LayerNorm,
    /// Self-attention (biases on q, k, v, out).
    pub attention: MultiHeadAttention,
    /// LayerNorm before the FFN.
    pub norm2: LayerNorm,
    /// FFN expansion `E → mlp`.
    pub ffn1: Linear,
    /// FFN projection `mlp → E`.
    pub ffn2: Linear,
}

/// The whole network.
#[derive(Module, Debug)]
pub struct PatNet {
    /// Input projection `6 → E`.
    pub input: Linear,
    /// LayerNorm after the input projection.
    pub input_norm: LayerNorm,
    /// The blocks.
    pub blocks: Vec<PatBlock>,
    /// Final LayerNorm.
    pub final_norm: LayerNorm,
    /// Output head `E → 6` on the center token.
    pub head: Linear,
}

fn linear(w: &LinearW, dev: &Device) -> Linear {
    let weight = Tensor::<2>::from_data(TensorData::new(w.w.clone(), [w.din, w.dout]), dev);
    let bias = Tensor::<1>::from_data(TensorData::new(w.b.clone(), [w.dout]), dev);
    Linear { weight: Param::from_tensor(weight), bias: Some(Param::from_tensor(bias)) }
}

fn norm(w: &NormW, eps: f64, dev: &Device) -> LayerNorm {
    let mut ln = LayerNormConfig::new(w.scale.len()).with_epsilon(eps).init(dev);
    ln.gamma = Param::from_tensor(Tensor::<1>::from_data(TensorData::new(w.scale.clone(), [w.scale.len()]), dev));
    ln.beta = Some(Param::from_tensor(Tensor::<1>::from_data(TensorData::new(w.bias.clone(), [w.bias.len()]), dev)));
    ln
}

impl PatNet {
    /// Build the module from framework-neutral weights on `dev`.
    pub fn from_weights(cfg: &ModelConfig, w: &NetWeights, dev: &Device) -> PatNet {
        let eps = super::LAYER_NORM_EPS;
        let blocks = w
            .blocks
            .iter()
            .map(|b| {
                let mut attention = MultiHeadAttentionConfig::new(cfg.embed, cfg.heads).with_dropout(0.0).init(dev);
                attention.query = linear(&b.q, dev);
                attention.key = linear(&b.k, dev);
                attention.value = linear(&b.v, dev);
                attention.output = linear(&b.o, dev);
                PatBlock {
                    norm1: norm(&b.norm1, eps, dev),
                    attention,
                    norm2: norm(&b.norm2, eps, dev),
                    ffn1: linear(&b.ffn1, dev),
                    ffn2: linear(&b.ffn2, dev),
                }
            })
            .collect();
        PatNet {
            input: linear(&w.input, dev),
            input_norm: norm(&w.input_norm, eps, dev),
            blocks,
            final_norm: norm(&w.final_norm, eps, dev),
            head: linear(&w.head, dev),
        }
    }

    /// Forward: features `[B, T, 6]` → raw outputs `[B, 6]` (before the σ scaling).
    /// `mask_pad` (`[B, T]`, true = padding) masks padded tokens in attention (training on padded
    /// clouds, D19); `None` for the pretrained weights.
    pub fn forward(&self, features: Tensor<3>, mask_pad: Option<Tensor<2, Bool>>) -> Tensor<2> {
        let [b, _t, _] = features.dims();
        let mut x = self.input_norm.forward(self.input.forward(features));
        for blk in &self.blocks {
            let h = blk.norm1.forward(x.clone());
            let mut inp = MhaInput::self_attn(h);
            if let Some(m) = &mask_pad {
                inp = inp.mask_pad(m.clone());
            }
            x = x + blk.attention.forward(inp).context;
            let h = blk.norm2.forward(x.clone());
            x = x + blk.ffn2.forward(activation::gelu_approximate(blk.ffn1.forward(h)));
        }
        let x = self.final_norm.forward(x);
        let e = x.dims()[2];
        let center: Tensor<2> = x.narrow(1, 0, 1).reshape([b, e]);
        self.head.forward(center)
    }
}

/// The CPU device of the Flex backend (built without its rayon feature: single-threaded; our pool
/// runs one chunk per worker, SPEC.md §5.4).
pub fn cpu_device() -> Device {
    Device::flex()
}

/// Run `net` on `features` (`b·t·6` values) and return `b·6` raw outputs.
pub fn forward_chunk(net: &PatNet, dev: &Device, features: &[f32], b: usize, t: usize) -> Result<Vec<f32>, String> {
    let x = Tensor::<3>::from_data(TensorData::new(features.to_vec(), [b, t, 6]), dev);
    let y = net.forward(x, None);
    y.into_data().try_into_vec::<f32>().map_err(|e| format!("{e:?}"))
}
