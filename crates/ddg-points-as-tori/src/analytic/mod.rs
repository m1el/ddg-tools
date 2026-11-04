//! Analytic test shapes (App. E, SPEC.md §11). Owned by WP8.
//!
//! * [`sdf2d`]: the 16 Quilez 2D signed distance functions with exact gradients
//!   `[formula: pat-quilez-2d]`;
//! * [`solid`]: revolution and extrusion with gradients `[formula: pat-revolution]`,
//!   `[formula: pat-extrusion]`, `[formula: pat-extrusion-gradient]`, area-uniform surface
//!   sampling with exact normals, and the [`ImplicitField`](crate::csg::ImplicitField) impl;
//! * [`boundary`]: traced 2D boundaries (also used by the 2D generators);
//! * [`set`]: the 70-instance set 14 × 5, seeds 0..69, and its JSON record;
//! * [`saddle`]: a small saddle-rich set and curvature statistics for the sign-rule question.
//!
//! The 2D formulas follow Inigo Quilez (<https://iquilezles.org/articles/distfunctions2d/>) and
//! the reference's `shape_2d.*` (MIT, Copyright (c) 2026 Nicole Feng, notice in
//! `LICENSES/points-as-tori-MIT.txt`).

pub mod boundary;
pub mod saddle;
pub mod sdf2d;
pub mod set;
pub mod solid;

pub use boundary::Boundary2d;
pub use saddle::{SaddleShape, SaddleStats, principal_curvatures, saddle_set, saddle_stats};
pub use sdf2d::{ALL_KINDS, Kind, PAPER_KINDS, Shape2d};
pub use set::{Instance, analytic_json, analytic_set, write_analytic_json};
pub use solid::{Solid, SolidKind};
