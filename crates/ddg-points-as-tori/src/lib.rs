//! Points as Tori: fast pointwise signed distance for oriented point clouds.
//!
//! An implementation of
//!
//! > Nicole Feng, Ioannis Gkioulekas, Keenan Crane. *Points as Tori: Fast Pointwise Signed
//! > Distance for Point Clouds.* ACM Transactions on Graphics 45(4), Article 53, 2026.
//! > <https://doi.org/10.1145/3811385>
//!
//! For every point of an oriented point cloud a small transformer predicts the six coefficients
//! of a local quadratic height function; a torus is fitted to its shift, principal curvatures and
//! directions, and the signed distance at any query point is a distance-weighted blend of the
//! tori's closed-form signed distance functions.
//!
//! The authors' reference implementation (`points-as-tori`, MIT, "Copyright (c) 2026 Nicole
//! Feng") may be read and ported; its licence notice is kept in `LICENSES/points-as-tori-MIT.txt`
//! and covers the converted pretrained weights as well.
//!
//! The contract of this crate is `SPEC.md` next to its `Cargo.toml` (formula ids `pat-…`, the
//! decisions, the public API, the verification and measurement plans and the work packages);
//! `NOTES-REFERENCE.md` records what the reference does where the paper is silent or differs.
//!
//! # Pipeline
//!
//! 1. [`PointCloud::new`] validates the cloud (finite, unit normals).
//! 2. [`spatial::fitting_neighborhoods`]: the 64 nearest points of every point, self first, by our
//!    deterministic kd-tree.
//! 3. [`network::predict_coefficients`]: features in the local frame scaled by the median distance,
//!    the pretrained transformer (f32), the coefficients of the height function.
//! 4. [`fit::fit_tori_with`]: a torus per point (f64), with every clamp counted.
//! 5. [`TorusField`]: the blended SDF, its gradient and Laplacian (f64).
//!
//! [`PointsAsTori::build`] runs 2–5. Every parallel stage uses an explicit pool; results are
//! bitwise independent of the thread count.
//!
//! ```no_run
//! use ddg_points_as_tori::{BuildOptions, Model, NormalPolicy, PointCloud, PointsAsTori};
//! # fn main() -> Result<(), ddg_points_as_tori::PatError> {
//! # let (positions, normals): (Vec<[f64; 3]>, Vec<[f64; 3]>) = (vec![], vec![]);
//! let (cloud, _) = PointCloud::new(positions, normals, NormalPolicy::default())?;
//! let model = Model::reference()?;
//! let (field, report) = PointsAsTori::build(&cloud, &model, &BuildOptions::default())?;
//! let (phi, grad) = field.distance_gradient([0.1, 0.2, 0.3]);
//! # let _ = (phi, grad, report); Ok(()) }
//! ```

#![warn(missing_docs)]

pub mod blend;
pub mod cloud;
pub mod error;
pub mod features;
pub mod field;
pub mod fit;
pub mod network;
pub mod npy;
pub mod pool;
pub mod sdf;
pub mod spatial;

// Modules of the later work packages (SPEC.md §17); stubs until they land.
#[allow(missing_docs)]
pub mod analytic;
#[allow(missing_docs)]
pub mod apps;
#[allow(missing_docs)]
pub mod baselines;
#[allow(missing_docs)]
pub mod csg;
#[allow(missing_docs)]
pub mod generate;
#[cfg(feature = "cuda")]
#[allow(missing_docs)]
pub mod gpu;
#[allow(missing_docs)]
pub mod mc;
#[allow(missing_docs)]
pub mod pat2d;
#[allow(missing_docs)]
pub mod render;
#[cfg(feature = "train")]
#[allow(missing_docs)]
pub mod train;
#[allow(missing_docs)]
pub mod truth;

pub use blend::{GradientMode, LambdaRule};
pub use cloud::{CloudReport, NormalPolicy, PointCloud};
pub use error::PatError;
pub use features::FeatureOrder;
pub use field::{BuildOptions, BuildReport, EvalOptions, EvalOut, EvalPrecision, EvalReport, Neighborhood, PointsAsTori, TorusField};
pub use fit::{FitDiagnostics, FitOptions, FitReport, Torus, Tori, fit_diagnostics, fit_tori};
pub use network::{Device, Model, ModelConfig, PredictOptions, PredictReport, predict_coefficients};
pub use spatial::Neighborhoods;
