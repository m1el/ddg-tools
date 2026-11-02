//! Comparison optimizers: AQP, H¹ L-BFGS, BCQN (`SPEC.md` §14.2, D26–D28; Sec. 8.2.1 of the
//! paper).
//!
//! All three use the cotan Laplacian `W` of the current geometry as the proxy (the session sets
//! `metric = H1` for them) with the positional constraints (barycenters, pins) in saddle form,
//! `W′ = [[W ⊗ I₃, C_Qᵀ], [C_Q, 0]]`, and need a soft problem: nonlinear constraints (area,
//! volume) become penalties ([`soft_problem`]). They share the session's energy and gradient,
//! line search, CCD, projection, remeshing and stop rules; only the direction differs
//! (and, for AQP, the point at which the gradient is taken).
//!
//! * **AQP** (D28; Kovalsky, Galun, Lipman, *Accelerated Quadratic Proxy for Geometric
//!   Optimization*, SIGGRAPH 2016): Nesterov momentum `y = (1+θ) x_k − θ x_{k−1}`,
//!   `θ = (1 − 1/√κ)/(1 + 1/√κ)`, then an H¹ step from `y`. The momentum is reset when the
//!   energy increases and on remeshing. [formula: rs-aqp]
//! * **H¹ L-BFGS** (D26): the two-loop recursion with `H₀ = γ W′⁻¹` (the curves crate's
//!   `Lbfgs`) — equal to the BFGS inverse update [formula: rs-lbfgs-two-loop]. Memory is reset
//!   on remeshing and on a pair without positive curvature. [`LbfgsVariant::Reference`](crate::step::LbfgsVariant): the
//!   reference's recursion with all inner products in `W` (`ρ = 1/yᵀWs`, …), which treats dual
//!   vectors as primal.
//! * **BCQN** (D27; Zhu, Bouaziz, Kavan *et al.*, *Blended Cured Quasi-Newton for Distortion
//!   Optimization*, SIGGRAPH 2018, without barriers or inversion-aware line search, as the
//!   paper): the stored gradient difference is blended `ŷ = (1−β) y + β W s`,
//!   `β = clamp(‖W‖ · yᵀWs / B, 0, 1)`, with `‖W‖` the spectral norm by power iteration and
//!   `B = (total area)^{4/3}`. [`BcqnVariant::Reference`](crate::step::BcqnVariant): Frobenius norm of the constrained
//!   Laplacian and the `W`-inner-product recursion. [formula: rs-bcqn]
//!
//! Determinism: everything here is serial and sparse, so trajectories do not depend on the
//! thread count.

mod soft;
mod state;

pub use soft::soft_problem;
pub use state::{MethodDiagnostics, aqp_theta};
#[doc(hidden)]
pub use state::{apply_k3, spectral_norm, two_loop_reference};
pub(crate) use state::MethodState;
