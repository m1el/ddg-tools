//! The coefficient network: weights, the inference API and the CPU back ends (SPEC.md §5, §6).
//!
//! The weights are kept framework neutral ([`NetWeights`], f32, Flax layouts flattened to
//! `[in, out]` row-major matrices). The CPU forward pass is the burn 0.22 module of
//! [`burn`](self::burn) on the Flex backend, built without its rayon feature: every call is
//! single-threaded, and our explicit pool runs fixed chunks of [`CHUNK`] neighborhoods, one per
//! task. Results are bitwise independent of the thread count (SPEC.md §5.4) and, as measured, of
//! the chunk size too (every op of the forward pass is row-wise per neighborhood).

pub mod burn;
#[cfg(feature = "cuda")]
pub mod cuda;

use crate::features::{FeatureOrder, FeatureScratch, neighborhood_features};
use crate::spatial::Neighborhoods;
use crate::{PatError, PointCloud};
use rayon::prelude::*;
use std::collections::HashMap;

/// Neighborhoods per network chunk (a constant, SPEC.md §5.4).
///
/// Measured in WP1 (NOTES-WP1.md): burn-flex runs 4,050/s at chunk 1, 3,300/s at 2, 2,700–3,000/s
/// at 4–8, 2,100/s at 16, 960/s at 32 and 560/s at 128 on 32 threads; the activations of a small
/// chunk stay in cache. The outputs were bitwise identical for every chunk size from 1 to 512.
pub const CHUNK: usize = 1;

/// LayerNorm ε of the reference (Flax `nnx.LayerNorm(epsilon = 1e-6)`).
pub const LAYER_NORM_EPS: f64 = 1e-6;

/// Network hyper-parameters (from the weights' metadata, checked against the tensors).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ModelConfig {
    /// Embedding width `E` (128).
    pub embed: usize,
    /// Transformer blocks (8).
    pub layers: usize,
    /// Attention heads (8).
    pub heads: usize,
    /// FFN width (512).
    pub mlp: usize,
    /// Neighborhood size the network was trained with (64).
    pub k: usize,
    /// Positional-encoding frequencies (0; other values are not supported, NOTES §2.2).
    pub n_freq: usize,
}

impl Default for ModelConfig {
    fn default() -> Self {
        ModelConfig { embed: 128, layers: 8, heads: 8, mlp: 512, k: 64, n_freq: 0 }
    }
}

/// A linear layer `y = x W + b`, `W` row-major `[din, dout]`.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearW {
    /// Input width.
    pub din: usize,
    /// Output width.
    pub dout: usize,
    /// `din·dout` weights, row-major `[din, dout]`.
    pub w: Vec<f32>,
    /// `dout` biases.
    pub b: Vec<f32>,
}

/// LayerNorm scale and bias.
#[derive(Debug, Clone, PartialEq)]
pub struct NormW {
    /// Scale (γ).
    pub scale: Vec<f32>,
    /// Bias (β).
    pub bias: Vec<f32>,
}

/// One block's weights. Attention projections are `[E, E]` with column `h·d + j` for head `h`,
/// head dimension `j` (Flax `[E, H, d]` flattened); the output projection's rows likewise.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct BlockW {
    pub norm1: NormW,
    pub q: LinearW,
    pub k: LinearW,
    pub v: LinearW,
    pub o: LinearW,
    pub norm2: NormW,
    pub ffn1: LinearW,
    pub ffn2: LinearW,
}

/// All weights.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct NetWeights {
    pub input: LinearW,
    pub input_norm: NormW,
    pub blocks: Vec<BlockW>,
    pub final_norm: NormW,
    pub head: LinearW,
}

/// A loaded model.
#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    config: ModelConfig,
    weights: NetWeights,
    metadata: HashMap<String, String>,
}

/// The bundled reference weights (`weights/pat-reference.safetensors`, MIT, see
/// `weights/README.md`).
#[cfg(feature = "bundled-weights")]
pub static REFERENCE_WEIGHTS: &[u8] = include_bytes!("../../weights/pat-reference.safetensors");

fn json_int(s: &str, key: &str) -> Option<usize> {
    let pat = format!("\"{key}\"");
    let at = s.find(&pat)? + pat.len();
    let rest = s[at..].trim_start().strip_prefix(':')?.trim_start();
    let end = rest.find(|c: char| !c.is_ascii_digit()).unwrap_or(rest.len());
    rest[..end].parse().ok()
}

impl Model {
    /// The authors' pretrained weights, bundled with the crate.
    #[cfg(feature = "bundled-weights")]
    pub fn reference() -> Result<Model, PatError> {
        Model::from_safetensors(REFERENCE_WEIGHTS)
    }

    /// Load from safetensors bytes (names and shapes as written by `tools/convert_weights.py`).
    pub fn from_safetensors(bytes: &[u8]) -> Result<Model, PatError> {
        let werr = |reason: String| PatError::Weights { reason };
        let (_, meta) = safetensors::SafeTensors::read_metadata(bytes).map_err(|e| werr(format!("header: {e}")))?;
        let metadata = meta.metadata().clone().unwrap_or_default();
        let st = safetensors::SafeTensors::deserialize(bytes).map_err(|e| werr(format!("{e}")))?;
        let mut config = ModelConfig::default();
        if let Some(c) = metadata.get("config") {
            let get = |k: &str| json_int(c, k).ok_or_else(|| werr(format!("metadata config lacks {k}: {c}")));
            config = ModelConfig {
                embed: get("embed_dim")?,
                layers: get("n_layers")?,
                heads: get("n_heads")?,
                mlp: get("mlp_dim")?,
                k: get("k_neighbors")?,
                n_freq: json_int(c, "n_freq").unwrap_or(0),
            };
        }
        if config.n_freq != 0 {
            return Err(werr(format!("n_freq = {} is not supported (the reference never computes it)", config.n_freq)));
        }
        if config.heads == 0 || config.embed % config.heads != 0 || config.k == 0 || config.layers == 0 {
            return Err(werr(format!("inconsistent config {config:?}")));
        }
        let get = |name: &str, shape: &[usize]| -> Result<Vec<f32>, PatError> {
            let t = st.tensor(name).map_err(|_| werr(format!("missing tensor {name}")))?;
            if t.dtype() != safetensors::Dtype::F32 {
                return Err(werr(format!("tensor {name} has dtype {:?}, expected F32", t.dtype())));
            }
            if t.shape() != shape {
                return Err(werr(format!("tensor {name} has shape {:?}, expected {shape:?}", t.shape())));
            }
            let v: Vec<f32> = t.data().chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
            if let Some(p) = v.iter().position(|x| !x.is_finite()) {
                return Err(werr(format!("tensor {name} has a non-finite value at flat index {p}")));
            }
            Ok(v)
        };
        let (e, h, m) = (config.embed, config.heads, config.mlp);
        let d = e / h;
        let lin = |p: &str, din: usize, dout: usize| -> Result<LinearW, PatError> {
            Ok(LinearW { din, dout, w: get(&format!("{p}.kernel"), &[din, dout])?, b: get(&format!("{p}.bias"), &[dout])? })
        };
        let nrm = |p: &str| -> Result<NormW, PatError> {
            Ok(NormW { scale: get(&format!("{p}.scale"), &[e])?, bias: get(&format!("{p}.bias"), &[e])? })
        };
        let qkv = |p: &str| -> Result<LinearW, PatError> {
            Ok(LinearW { din: e, dout: e, w: get(&format!("{p}.kernel"), &[e, h, d])?, b: get(&format!("{p}.bias"), &[h, d])? })
        };
        let mut blocks = Vec::with_capacity(config.layers);
        for l in 0..config.layers {
            let p = format!("blocks.{l}");
            blocks.push(BlockW {
                norm1: nrm(&format!("{p}.norm1"))?,
                q: qkv(&format!("{p}.attention.query"))?,
                k: qkv(&format!("{p}.attention.key"))?,
                v: qkv(&format!("{p}.attention.value"))?,
                o: LinearW {
                    din: e,
                    dout: e,
                    w: get(&format!("{p}.attention.out.kernel"), &[h, d, e])?,
                    b: get(&format!("{p}.attention.out.bias"), &[e])?,
                },
                norm2: nrm(&format!("{p}.norm2"))?,
                ffn1: lin(&format!("{p}.ffn.layers.0"), e, m)?,
                ffn2: lin(&format!("{p}.ffn.layers.3"), m, e)?,
            });
        }
        let weights = NetWeights {
            input: lin("input_encoder.mlp.layers.0", 6, e)?,
            input_norm: nrm("input_encoder.mlp.layers.1")?,
            blocks,
            final_norm: nrm("final_norm")?,
            head: lin("mlp.layers.0", e, 6)?,
        };
        let expected = 8 + 16 * config.layers;
        if st.len() != expected {
            let mut names: Vec<&str> = st.names();
            names.sort_unstable();
            return Err(werr(format!("{} tensors, expected {expected}: {}", st.len(), names.join(", "))));
        }
        Ok(Model { config, weights, metadata })
    }

    /// Build a model from weights (e.g. after training).
    pub fn from_weights(config: ModelConfig, weights: NetWeights) -> Model {
        Model { config, weights, metadata: HashMap::new() }
    }

    /// The configuration.
    pub fn config(&self) -> &ModelConfig {
        &self.config
    }

    /// The weights.
    pub fn weights(&self) -> &NetWeights {
        &self.weights
    }

    /// The safetensors metadata (source, licence, hashes).
    pub fn metadata(&self) -> &HashMap<String, String> {
        &self.metadata
    }
}

/// Where the network runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Device {
    /// The CPU.
    #[default]
    Cpu,
    /// A CUDA device (feature `cuda`, WP4).
    Cuda {
        /// Device ordinal.
        ordinal: usize,
    },
    /// CUDA when available, else the CPU.
    Auto,
}

/// Options of [`predict_coefficients`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PredictOptions {
    /// Device.
    pub device: Device,
    /// Threads (`None`: all cores).
    pub threads: Option<usize>,
    /// Feature arithmetic (D3).
    pub feature_order: FeatureOrder,
}

/// What prediction did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PredictReport {
    /// The device used.
    pub device: Device,
    /// Tokens per neighborhood.
    pub tokens: usize,
    /// Number of chunks.
    pub chunks: usize,
    /// Seconds in features (summed over workers).
    pub features_seconds: f64,
    /// Seconds in the forward pass (summed over workers).
    pub forward_seconds: f64,
    /// Wall time.
    pub seconds: f64,
}

/// Scale raw outputs `y` to coefficients `a = (σ y0, y1, y2, y3/σ, y4/σ, y5/σ)`. In
/// [`FeatureOrder::Reference`] the arithmetic is f32 (the reference's); otherwise f64, then
/// rounded.
///
/// [formula: pat-output-scaling]
pub fn scale_outputs(y: [f32; 6], sigma: f64, order: FeatureOrder) -> [f32; 6] {
    match order {
        FeatureOrder::Reference => {
            let s = sigma as f32;
            [y[0] * s, y[1], y[2], y[3] / s, y[4] / s, y[5] / s]
        }
        FeatureOrder::DifferenceFirst => scale_outputs_f64(y.map(|v| v as f64), sigma).map(|v| v as f32),
    }
}

/// [`scale_outputs`] in f64 without rounding (the default order's arithmetic before the cast).
///
/// [formula: pat-output-scaling]
pub fn scale_outputs_f64(y: [f64; 6], sigma: f64) -> [f64; 6] {
    [y[0] * sigma, y[1], y[2], y[3] / sigma, y[4] / sigma, y[5] / sigma]
}

/// Raw network outputs for precomputed features (`b·t·6` values, `b` neighborhoods of `t`
/// tokens), chunked as in [`predict_coefficients`]. For parity tests and benchmarks.
pub fn forward_features(model: &Model, features: &[f32], tokens: usize, opts: &PredictOptions) -> Result<Vec<f32>, PatError> {
    if tokens == 0 || features.len() % (6 * tokens) != 0 {
        return Err(PatError::SizeMismatch { what: format!("{} feature values for {tokens} tokens", features.len()) });
    }
    check_device(opts.device)?;
    let b = features.len() / (6 * tokens);
    let pool = crate::pool::make_pool(opts.threads)?;
    let runner = Runner::new(model);
    let out: Result<Vec<Vec<f32>>, PatError> = pool.install(|| {
        features
            .par_chunks(CHUNK * tokens * 6)
            .map(|f| runner.run(f, f.len() / (6 * tokens), tokens))
            .collect()
    });
    let out: Vec<f32> = out?.concat();
    debug_assert_eq!(out.len(), b * 6);
    Ok(out)
}

fn check_device(d: Device) -> Result<(), PatError> {
    match d {
        Device::Cpu | Device::Auto => Ok(()),
        Device::Cuda { .. } => Err(PatError::Gpu { reason: "the cuda back end is not built yet (WP4)".into() }),
    }
}

/// The framework module and device, shared by the workers (the module is `Sync`).
pub(crate) struct Runner {
    net: self::burn::PatNet,
    dev: ::burn::tensor::Device,
}

impl Runner {
    pub(crate) fn new(model: &Model) -> Self {
        let dev = self::burn::cpu_device();
        Runner { net: self::burn::PatNet::from_weights(&model.config, &model.weights, &dev), dev }
    }

    /// Raw outputs of `b` neighborhoods of `t` tokens.
    pub(crate) fn run(&self, features: &[f32], b: usize, t: usize) -> Result<Vec<f32>, PatError> {
        self::burn::forward_chunk(&self.net, &self.dev, features, b, t).map_err(|reason| PatError::Backend { reason })
    }
}

/// Coefficients and scales `σ_i` for every neighborhood (SPEC.md §6.3).
///
/// Rows of `nbhds` must start with their own index (as [`crate::spatial::fitting_neighborhoods`]
/// builds them). Non-finite outputs are an error naming the rows.
pub fn predict_coefficients(
    model: &Model,
    cloud: &PointCloud,
    nbhds: &Neighborhoods,
    opts: &PredictOptions,
) -> Result<(Vec<[f32; 6]>, Vec<f64>, PredictReport), PatError> {
    let t0 = std::time::Instant::now();
    check_device(opts.device)?;
    let n = cloud.len();
    if nbhds.len() != n {
        return Err(PatError::SizeMismatch { what: format!("{} neighborhoods for {n} points", nbhds.len()) });
    }
    let t = nbhds.k;
    if t > model.config.k {
        return Err(PatError::SizeMismatch {
            what: format!("neighborhoods of {t} tokens, the model was trained with {}", model.config.k),
        });
    }
    for i in 0..n {
        let row = nbhds.row(i);
        if row[0] as usize != i {
            return Err(PatError::SizeMismatch { what: format!("neighborhood {i} starts with {}, not with itself", row[0]) });
        }
        if let Some(&j) = row.iter().find(|&&j| j as usize >= n) {
            return Err(PatError::SizeMismatch { what: format!("neighborhood {i} has index {j} ≥ {n}") });
        }
    }
    let pool = crate::pool::make_pool(opts.threads)?;
    let runner = Runner::new(model);
    let (pos, nrm) = (cloud.positions(), cloud.normals());
    let n_chunks = n.div_ceil(CHUNK);
    type ChunkOut = (Vec<[f32; 6]>, Vec<f64>, f64, f64);
    let res: Result<Vec<ChunkOut>, PatError> = pool.install(|| {
        (0..n_chunks)
            .into_par_iter()
            .map_init(
                || (FeatureScratch::default(), vec![0f32; CHUNK * t * 6]),
                |(scr, feat), c| {
                    let lo = c * CHUNK;
                    let hi = (lo + CHUNK).min(n);
                    let b = hi - lo;
                    let tf = std::time::Instant::now();
                    let mut sig = Vec::with_capacity(b);
                    for (r, i) in (lo..hi).enumerate() {
                        let s = neighborhood_features(pos, nrm, nbhds.row(i), opts.feature_order, scr, &mut feat[r * t * 6..(r + 1) * t * 6]);
                        sig.push(s);
                    }
                    let tf = tf.elapsed().as_secs_f64();
                    let tn = std::time::Instant::now();
                    let raw = runner.run(&feat[..b * t * 6], b, t)?;
                    let tn = tn.elapsed().as_secs_f64();
                    let co = (0..b)
                        .map(|r| {
                            let y = [raw[6 * r], raw[6 * r + 1], raw[6 * r + 2], raw[6 * r + 3], raw[6 * r + 4], raw[6 * r + 5]];
                            scale_outputs(y, sig[r], opts.feature_order)
                        })
                        .collect();
                    Ok((co, sig, tf, tn))
                },
            )
            .collect()
    });
    let mut coeffs = Vec::with_capacity(n);
    let mut sigma = Vec::with_capacity(n);
    let mut report = PredictReport { device: Device::Cpu, tokens: t, chunks: n_chunks, ..Default::default() };
    for (co, sg, tf, tn) in res? {
        coeffs.extend(co);
        sigma.extend(sg);
        report.features_seconds += tf;
        report.forward_seconds += tn;
    }
    let bad: Vec<u32> = coeffs.iter().enumerate().filter(|(_, a)| !a.iter().all(|x| x.is_finite())).map(|(i, _)| i as u32).collect();
    if !bad.is_empty() {
        return Err(PatError::NonFiniteOutput { indices: bad });
    }
    report.seconds = t0.elapsed().as_secs_f64();
    Ok((coeffs, sigma, report))
}
