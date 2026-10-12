//! Experiments around the optimal-cones solver (not published).
//!
//! * [`prepare`]: the reduced problem `µ -> u` through one fixed Cholesky factor of the Laplacian.
//! * [`lasso`]: a factor-once solver (working set + dense active set) for the default problem.
//! * [`dense`]: the dense active-set method.
//! * [`meshes`]: loading and generating test meshes.

pub mod dense;
pub mod lasso;
pub mod meshes;
pub mod prepare;
