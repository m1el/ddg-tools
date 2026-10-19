//! Optimizers for the developability flow (SPEC §4): L-BFGS with the Lewis–Overton weak Wolfe
//! line search, gradient descent with a line search, and a fixed step.
//!
//! The optimizer is generic over an [`Objective`] on a flat `&[f64]` (groups of `group`
//! coordinates are "vertices"; `group = 3` for meshes), so it can be tested on classic smooth and
//! nonsmooth functions. [`minimize`] is the generic entry point. Lengths are measured in units of
//! the problem's `length_scale` (the mean edge length `ℓ̄` for meshes, D17); positions are never
//! rescaled.
//!
//! References: Lewis & Overton, *Nonsmooth optimization via quasi-Newton methods*, Math.
//! Program. 141 (2013), Algorithm 2.6 and Sec. 5.1; Nocedal & Wright, *Numerical Optimization*
//! (2006), Alg. 3.1, 7.4, 7.5.

mod core;
mod options;

pub use self::core::{
    ClosureObjective, EvalError, Minimized, Objective, OptimizeError, Problem, Undefined, minimize,
};
pub use self::options::{
    Control, FlatProgress, InvalidOption, LineSearch, Method, OptimizeOptions, OptimizeReport,
    OptimizeStatus,
};
