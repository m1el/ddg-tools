//! Linear solves for the constrained direction and the projection (`SPEC.md` §6).
//!
//! Stub created by WP1; filled by WP2 (`dense`, `schur`, `gmres`) and WP6 (`iterative`,
//! `precond`, SPEC §15).

pub mod dense;
pub mod gmres;
pub mod iterative;
pub mod precond;
pub mod schur;
