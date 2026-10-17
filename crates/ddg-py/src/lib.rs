//! Python module `ddg_tools`: the bridge between Blender and the Rust crates of this workspace.
//!
//! A mesh goes in as flat arrays taken straight from Blender (`foreach_get`), a UV per polygon
//! corner comes out as an array Blender takes back (`foreach_set`). No numpy is needed on either
//! side: arguments are read through the buffer protocol, results are `array.array` objects.
//!
//! The unwrap pipeline behind this binding is the `ddg-uv` crate: seams, optimal cone
//! singularities, cutting, conformal flattening. [`pipeline::run_pipeline`] is the single place
//! that calls it.
//!
//! # Layout
//!
//! * [`buffers`] — decoding buffer-protocol memory into typed vectors.
//! * [`cones`] — `reoptimize_cones`: the optimal angles of cones at given vertices.
//! * [`contract`] — [`contract::UnwrapInput`], [`contract::UnwrapOptions`],
//!   [`contract::UnwrapOutput`] and [`contract::Report`]: the data that crosses the boundary,
//!   validated.
//! * [`error`] — [`error::InputError`] and [`error::UnwrapError`].
//! * [`pipeline`] — the seam to the pipeline and the driver around it.
//! * [`measure`] — distortion measured on the resulting UV map.
//! * [`pack`] — packing charts into the unit square.
//! * `python` (cargo feature `python`) — the PyO3 glue. It contains no logic of its own:
//!   it turns Python buffers into [`buffers::BufferView`]s, calls the functions above and
//!   turns the results into Python objects.
//!
//! Everything except `python` is plain Rust and is what `cargo test -p ddg-py` exercises; no
//! Python interpreter is needed to build or test it. The Python-facing documentation (signature,
//! buffer formats, report keys) is in the crate README.
//!
//! # References
//!
//! The pipeline this binding is meant for implements
//! Yousuf Soliman, Dejan Slepčev, Keenan Crane, *Optimal Cone Singularities for Conformal
//! Flattening*, ACM Transactions on Graphics 37(4), 2018,
//! <https://doi.org/10.1145/3197517.3201367>, and
//! Rohan Sawhney, Keenan Crane, *Boundary First Flattening*, ACM Transactions on Graphics 37(1),
//! 2017, <https://doi.org/10.1145/3132705>. Neither algorithm is in this crate.

#![warn(missing_docs)]

pub mod buffers;
pub mod cones;
pub mod contract;
pub mod cutting;
pub mod error;
pub mod intrinsic_apps;
pub mod measure;
pub mod pack;
pub mod pipeline;

#[cfg(feature = "python")]
mod python;
