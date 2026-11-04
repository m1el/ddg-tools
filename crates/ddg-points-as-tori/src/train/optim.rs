//! The optimizer and the learning-rate schedules (SPEC.md §10.1, D21; NOTES-REFERENCE.md §9.4).
//!
//! The reference uses `optax.adamw(warmup_cosine_decay_schedule(0, lr, warmup, 1e8, 0.01·lr),
//! weight_decay = 0)`: Adam (β = 0.9/0.999, ε = 1e-8) with a linear warmup from 0 and a cosine over
//! `1e8` steps, i.e. practically constant after the warmup. The moment update is burn's
//! [`Adam`](burn::optim::Adam) applied per parameter tensor (its bias-corrected step equals optax's
//! `m̂/(√v̂ + ε)`); AdamW's decoupled decay `θ ← θ − lr·wd·θ` is added on top when `weight_decay > 0`
//! (optax's `add_decayed_weights` before the learning rate).
//!
//! [formula: pat-adam]

use super::net::{NetParams, PT};
use super::tfit::Ctx;
use crate::PatError;
use burn::optim::{Adam, AdamConfig, AdamState, AdaptiveMomentumState, Optimizer};
use burn::tensor::{Tensor, TensorData};
use std::collections::HashMap;

/// Adam/AdamW options.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdamOptions {
    /// β₁ (0.9).
    pub beta1: f32,
    /// β₂ (0.999).
    pub beta2: f32,
    /// ε (1e-8, optax's; burn's default is 1e-5).
    pub epsilon: f32,
    /// Decoupled weight decay (0: Adam, as the reference).
    pub weight_decay: f64,
}

impl Default for AdamOptions {
    fn default() -> Self {
        AdamOptions { beta1: 0.9, beta2: 0.999, epsilon: 1e-8, weight_decay: 0.0 }
    }
}

/// The schedule after the warmup (D21).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Schedule {
    /// The code's: cosine over `1e8` steps, i.e. constant (default).
    #[default]
    Code,
    /// A real cosine from the peak to `0.01·peak` over the phase's steps (the paper's intent).
    Cosine,
}

/// optax's `warmup_cosine_decay_schedule(0, peak, warmup, decay_steps, 0.01·peak)` at `count`.
pub fn warmup_cosine(count: u64, peak: f64, warmup: u64, decay_steps: u64) -> f64 {
    let end = 0.01 * peak;
    if warmup > 0 && count < warmup {
        return peak * count as f64 / warmup as f64;
    }
    let total = decay_steps.saturating_sub(warmup).max(1) as f64;
    let t = ((count - warmup) as f64).min(total);
    let alpha = end / peak;
    let cos = 0.5 * (1.0 + (std::f64::consts::PI * t / total).cos());
    peak * ((1.0 - alpha) * cos + alpha)
}

/// The learning rate of step `count` (0-based) of a phase of `phase_steps` steps.
pub fn learning_rate(s: Schedule, count: u64, peak: f64, warmup: u64, phase_steps: u64) -> f64 {
    match s {
        Schedule::Code => warmup_cosine(count, peak, warmup, 100_000_000),
        Schedule::Cosine => warmup_cosine(count, peak, warmup, warmup + phase_steps),
    }
}

/// A named f32 tensor: name, shape, values.
pub type NamedTensor = (String, Vec<usize>, Vec<f32>);

#[derive(Clone)]
enum St {
    D1(AdamState<1>),
    D2(AdamState<2>),
}

/// Adam state for every parameter (keyed by the parameter's Flax name).
#[derive(Clone)]
pub struct AdamTrainer {
    adam: Adam,
    opts: AdamOptions,
    state: Vec<Option<St>>,
    /// Updates done.
    pub count: u64,
}

impl AdamTrainer {
    /// Fresh state for `params`.
    pub fn new(params: &NetParams, opts: AdamOptions) -> AdamTrainer {
        let adam = AdamConfig::new().with_beta_1(opts.beta1).with_beta_2(opts.beta2).with_epsilon(opts.epsilon).build();
        AdamTrainer { adam, opts, state: vec![None; params.items.len()], count: 0 }
    }

    /// Apply one update with `grads` (canonical order) at learning rate `lr`.
    pub fn update(&mut self, params: &NetParams, grads: &[PT], lr: f64) -> Result<NetParams, PatError> {
        if grads.len() != params.items.len() {
            return Err(PatError::SizeMismatch { what: format!("{} gradients for {} parameters", grads.len(), params.items.len()) });
        }
        let mut vals = Vec::with_capacity(grads.len());
        for (i, (p, g)) in params.items.iter().zip(grads).enumerate() {
            let wd = self.opts.weight_decay;
            let v = match (&p.value, g, self.state[i].take()) {
                (PT::D1(t), PT::D1(gr), st) => {
                    let st = match st {
                        Some(St::D1(s)) => Some(s),
                        _ => None,
                    };
                    let (mut nt, ns) = self.adam.step(lr, t.clone(), gr.clone(), st);
                    if wd > 0.0 {
                        nt = nt - t.clone().mul_scalar(lr * wd);
                    }
                    self.state[i] = ns.map(St::D1);
                    PT::D1(nt)
                }
                (PT::D2(t), PT::D2(gr), st) => {
                    let st = match st {
                        Some(St::D2(s)) => Some(s),
                        _ => None,
                    };
                    let (mut nt, ns) = self.adam.step(lr, t.clone(), gr.clone(), st);
                    if wd > 0.0 {
                        nt = nt - t.clone().mul_scalar(lr * wd);
                    }
                    self.state[i] = ns.map(St::D2);
                    PT::D2(nt)
                }
                _ => return Err(PatError::SizeMismatch { what: format!("gradient rank of {}", p.name) }),
            };
            vals.push(v);
        }
        self.count += 1;
        Ok(params.with_values(vals))
    }

    /// The moments as named f32 tensors (`m.<name>`, `v.<name>`) and the per-tensor step count.
    pub fn to_named(&self, params: &NetParams) -> (Vec<NamedTensor>, u64) {
        let mut out = Vec::new();
        let mut time = 0u64;
        for (p, s) in params.items.iter().zip(&self.state) {
            let (m1, m2, t, shape): (Vec<f32>, Vec<f32>, usize, Vec<usize>) = match s {
                Some(St::D1(s)) => (f32v(s.momentum.moment_1.clone()), f32v(s.momentum.moment_2.clone()), s.momentum.time, s.momentum.moment_1.dims().to_vec()),
                Some(St::D2(s)) => (f32v(s.momentum.moment_1.clone()), f32v(s.momentum.moment_2.clone()), s.momentum.time, s.momentum.moment_1.dims().to_vec()),
                None => continue,
            };
            time = t as u64;
            out.push((format!("m.{}", p.name), shape.clone(), m1));
            out.push((format!("v.{}", p.name), shape, m2));
        }
        (out, time)
    }

    /// Restore the moments written by [`AdamTrainer::to_named`] (`time` the moment count).
    pub fn from_named(params: &NetParams, opts: AdamOptions, named: &HashMap<String, Vec<f32>>, time: u64, count: u64, ctx: &Ctx) -> Result<AdamTrainer, PatError> {
        let mut t = AdamTrainer::new(params, opts);
        t.count = count;
        if time == 0 {
            return Ok(t);
        }
        for (i, p) in params.items.iter().enumerate() {
            let (Some(m1), Some(m2)) = (named.get(&format!("m.{}", p.name)), named.get(&format!("v.{}", p.name))) else {
                return Err(PatError::Format { reason: format!("optimizer state lacks the moments of {}", p.name) });
            };
            if m1.len() != p.value.len() || m2.len() != p.value.len() {
                return Err(PatError::Format { reason: format!("optimizer moments of {} have the wrong size", p.name) });
            }
            t.state[i] = Some(match &p.value {
                PT::D1(x) => {
                    let mk = |v: &Vec<f32>| Tensor::<1>::from_data(TensorData::new(v.clone(), x.dims()), (&ctx.dev, ctx.dtype));
                    St::D1(AdamState::new(AdaptiveMomentumState::new(time as usize, mk(m1), mk(m2))))
                }
                PT::D2(x) => {
                    let mk = |v: &Vec<f32>| Tensor::<2>::from_data(TensorData::new(v.clone(), x.dims()), (&ctx.dev, ctx.dtype));
                    St::D2(AdamState::new(AdaptiveMomentumState::new(time as usize, mk(m1), mk(m2))))
                }
            });
        }
        Ok(t)
    }
}

fn f32v<const D: usize>(t: Tensor<D>) -> Vec<f32> {
    t.into_data().convert::<f32>().try_to_vec::<f32>().unwrap_or_default()
}
