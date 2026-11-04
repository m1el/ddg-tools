//! Blob generator, simulated scanner, farthest-point sampling, training files (SPEC.md §10.3,
//! §11). Owned by WP8.
//!
//! * [`rng`]: seeded generator; [`perlin`]: improved noise (our own);
//! * [`blob`]: Perlin-noise implicit blobs, marching cubes, watertight check ([`simplify`] brings
//!   them to 2048 faces);
//! * [`bvh`], [`scanner`]: ray casting and the simulated scanner;
//! * [`sampling`]: area-uniform mesh sampling, farthest points, centering;
//! * [`npz`], [`dataset`]: training/validation files in the authors' layout;
//! * [`d2`]: the 2D generators for the trained 2D network (WP9b).
#![allow(clippy::type_complexity, clippy::neg_cmp_op_on_partial_ord)]

pub mod blob;
pub mod bvh;
pub mod d2;
pub mod dataset;
pub mod npz;
pub mod perlin;
pub mod rng;
pub mod sampling;
pub mod scanner;
pub mod simplify;
