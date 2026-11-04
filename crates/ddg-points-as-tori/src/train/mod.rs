//! Training (feature `train`, SPEC.md §10; NOTES-REFERENCE.md §9; NOTES-WP7.md).
//!
//! A port of the reference's training code (`training/train.py`, `mesh_sampling_jax.py` and the loss
//! parts of `src/pointsastori/network.py` of points-as-tori, MIT License, Copyright (c) 2026 Nicole
//! Feng; notice in `LICENSES/points-as-tori-MIT.txt`) to burn 0.22 with autodiff. Everything runs on
//! a burn [`Device`](burn::tensor::Device) given by the caller: the CPU (`Device::flex()`) here, a
//! CUDA device in WP10. Autodiff is enabled per leaf tensor, so the same code serves both.
//!
//! - [`npz`], [`data`]: the training files, validated (D19, D23).
//! - [`queries`], [`truth`]: per-row query sampling (D20, D22) and the mesh ground truth.
//! - [`tfit`], [`loss`]: the torus fit and the loss of Eq. 27 as framework ops, with the blend
//!   weights as host constants.
//! - [`net`]: the parameters, the initialization, the burn module per pass, safetensors.
//! - [`step`]: one training step in two passes (SPEC.md §10.4) and the evaluation of the loss.
//! - [`optim`]: Adam/AdamW (burn's Adam per tensor) and the schedules (D21).
//! - [`run`]: the curriculum, the loop, checkpoints, logs and resumption.

pub mod data;
pub mod loss;
pub mod net;
pub mod npz;
pub mod optim;
pub mod queries;
pub mod run;
pub mod step;
pub mod tfit;
pub mod truth;

#[cfg(feature = "cuda")]
pub mod gpu;

pub use data::{DataOptions, DataReport, TrainShape, load_file};
pub use loss::{BlendPlan, LossOptions, LossValue};
pub use net::NetParams;
pub use optim::{AdamOptions, Schedule};
pub use queries::{QueryBoxes, QueryOptions};
pub use run::{Curriculum, Phase, StepLog, TrainOptions, TrainReport, train};
pub use step::{StepOptions, StepOut};
pub use tfit::{Ctx, TensorFitOptions};
