//! The optimality certificate: residuals of the optimality system and a duality gap that
//! bounds the suboptimality, computed from `(u, φ, µ)` and the problem data alone.
//!
//! Nothing here uses the solver's state: the Laplacian is applied face by face from the corner
//! cotangents (not through the assembled matrix or its interior block), the cone equation is
//! evaluated from the problem's weights and bounds (not from the solver's active sets), and the
//! dual point is rebuilt from `φ`.

use ddg_core::operators::{angle_defects, lumped_mass_normalized};
use ddg_core::{Geometry, TriMesh};

use crate::error::ConeError;
use crate::problem::{ConeProblem, validate};

/// Residuals of the optimality conditions and the duality gap of a triple `(u, φ, µ)`.
///
/// With `w_i = wE_i m_i`, `t_i = λ wR_i` and `g_i(µ) = t_i |µ| + ι_[lo_i, hi_i](µ)` the
/// optimality system of the canonical problem is
///
/// ```text
/// (L u)_i = µ_i − Ω_i        (L φ)_i = −w_i u_i        φ_i ∈ ∂g_i(µ_i)        (i interior)
/// ```
///
/// with `u = φ = µ = 0` on the boundary. A triple that satisfies it is a global minimizer
/// (machine-checked: `optimality-system-sufficient` in `lean/proofs`).
///
/// The **duality gap** makes this quantitative. Let `φ̂` be `φ` clipped to the domain of the
/// conjugates (`φ̂_i ≤ t_i` where `hi_i = +∞`, `φ̂_i ≥ −t_i` where `lo_i = −∞`), and
/// `ũ_i = −(L φ̂)_i / w_i`. Then `(ũ, φ̂)` is a feasible point of the dual problem, its value
///
/// ```text
/// D = −½ Σ (L φ̂)_i² / w_i + Σ φ̂_i Ω_i − Σ g_i*(φ̂_i)
/// ```
///
/// is a lower bound for the objective of every feasible point, and so
/// `objective(u, µ) − optimum ≤ primal_objective − dual_objective = duality_gap` provided
/// `(u, µ)` is feasible, i.e. up to `state_residual` and `bound_violation`
/// (`duality-gap-certificate` in `lean/proofs`, theorems `duality_gap_certificate` and
/// `adjoint_gap_certificate`). At a vertex with `w_i = 0` the dual point needs `(L φ̂)_i = 0`,
/// which cannot be enforced by construction; there the bound holds up to `adjoint_residual`.
///
/// All statements about bounds are for exact arithmetic; the values are computed in `f64`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Certificate {
    /// `max |(L u)_i − µ_i + Ω_i|` over interior vertices (and `|u_i|` over boundary vertices).
    pub state_residual: f64,
    /// `max |(L φ)_i + wE_i m_i u_i|` over interior vertices (and `|φ_i|` over boundary
    /// vertices).
    pub adjoint_residual: f64,
    /// `max |µ_i − clamp(P_{t_i}(φ_i + µ_i), lo_i, hi_i)|` over interior vertices: zero exactly
    /// when the inclusion `φ_i ∈ ∂g_i(µ_i)` holds everywhere.
    pub cone_residual: f64,
    /// Largest distance of a cone angle from its admissible set: `[lo_i, hi_i]`, the prescribed
    /// angle, `0` at a forbidden or boundary vertex.
    pub bound_violation: f64,
    /// `max |φ̂_i − φ_i|`: how far `φ` had to be clipped to become dual feasible.
    pub dual_infeasibility: f64,
    /// Closed meshes: `|Σ µ_i − Σ Ω_i|`, the violation of Gauss–Bonnet. Zero on meshes with
    /// boundary.
    pub gauss_bonnet_residual: f64,
    /// The objective `½ Σ wE m u² + λ Σ wR |µ|` of `(u, µ)`.
    pub primal_objective: f64,
    /// The dual value `D` of `(ũ, φ̂)`: a lower bound of the optimal objective.
    pub dual_objective: f64,
    /// `primal_objective − dual_objective`: the bound on the suboptimality of `(u, µ)`.
    pub duality_gap: f64,
    /// The same gap evaluated without cancellation, as the sum of the nonnegative terms
    /// `½ Σ w_i (u_i − ũ_i)² + Σ (g_i(µ_i) + g_i*(φ̂_i) − φ̂_i µ_i)`. It differs from
    /// `duality_gap` by `Σ φ̂_i · (state residual)_i` and by rounding.
    pub complementarity_gap: f64,
}

impl Certificate {
    /// Whether every residual is at most `tolerance` and both forms of the gap are at most
    /// `tolerance · max(1, |primal_objective|)`.
    pub fn holds(&self, tolerance: f64) -> bool {
        let gap_tolerance = tolerance * self.primal_objective.abs().max(1.0);
        // Written so that a NaN fails.
        self.state_residual <= tolerance
            && self.adjoint_residual <= tolerance
            && self.cone_residual <= tolerance
            && self.bound_violation <= tolerance
            && self.dual_infeasibility <= tolerance
            && self.gauss_bonnet_residual <= tolerance
            && self.duality_gap.abs() <= gap_tolerance
            && self.complementarity_gap.abs() <= gap_tolerance
    }
}

/// Applies the cotan Laplacian to a per-vertex vector, face by face:
/// `(L x)_i = ½ Σ_faces Σ_{corners k of the face} cot θ_k (x_i − x_j)` for the edge `ij`
/// opposite corner `k`.
///
/// [formula: cotan-laplacian]
pub(crate) fn apply_laplacian(mesh: &TriMesh, geometry: &Geometry, x: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; x.len()];
    for (face, cotans) in mesh.faces().iter().zip(geometry.corner_cotans()) {
        for k in 0..3 {
            let (i, j) = (face[(k + 1) % 3] as usize, face[(k + 2) % 3] as usize);
            let flow = 0.5 * cotans[k] * (x[i] - x[j]);
            out[i] += flow;
            out[j] -= flow;
        }
    }
    out
}

fn check_len(what: &'static str, x: &[f64], n: usize) -> Result<(), ConeError> {
    if x.len() != n {
        return Err(ConeError::WrongLength {
            what,
            expected: n,
            got: x.len(),
        });
    }
    Ok(())
}

/// Computes the [`Certificate`] of a triple `(u, φ, µ)` for a problem. Each argument has one
/// entry per mesh vertex.
///
/// This is the check [`solve`](crate::solve) runs on its own result; it can equally be used on
/// a triple from any other source.
///
/// # Errors
/// The errors of the problem's validation, and [`ConeError::WrongLength`].
///
/// [formula: optimality-system-sufficient]
/// [formula: duality-gap-certificate]
/// [formula: cone-fixed-point-bounded]
pub fn certify(
    problem: &ConeProblem<'_>,
    log_scale: &[f64],
    adjoint: &[f64],
    measure: &[f64],
) -> Result<Certificate, ConeError> {
    let data = validate(problem)?;
    let mesh = problem.mesh;
    let geometry = problem.geometry;
    let nv = mesh.n_vertices();
    check_len("log_scale", log_scale, nv)?;
    check_len("adjoint", adjoint, nv)?;
    check_len("measure", measure, nv)?;
    let (u, phi, mu) = (log_scale, adjoint, measure);

    let omega = angle_defects(mesh, geometry)?;
    let mass = lumped_mass_normalized(mesh, geometry)?;
    let closed = data.unknown.iter().all(|&b| b);

    // φ clipped to the domain of the conjugates: a dual-feasible point.
    let phi_hat: Vec<f64> = (0..nv)
        .map(|v| {
            if data.unknown[v] {
                let (a, b) = data.sites[v].dual_domain();
                phi[v].clamp(a, b)
            } else {
                phi[v]
            }
        })
        .collect();
    let l_u = apply_laplacian(mesh, geometry, u);
    let l_phi = apply_laplacian(mesh, geometry, phi);
    let l_phi_hat = apply_laplacian(mesh, geometry, &phi_hat);

    let mut c = Certificate {
        state_residual: 0.0,
        adjoint_residual: 0.0,
        cone_residual: 0.0,
        bound_violation: 0.0,
        dual_infeasibility: 0.0,
        gauss_bonnet_residual: 0.0,
        primal_objective: 0.0,
        dual_objective: 0.0,
        duality_gap: 0.0,
        complementarity_gap: 0.0,
    };
    let mut distortion = 0.0;
    let mut cone_cost = 0.0;
    let mut dual_quadratic = 0.0;
    let mut dual_linear = 0.0;
    let mut dual_conjugate = 0.0;
    let mut mismatch = 0.0;
    let mut slack = 0.0;
    let mut sum_mu = 0.0;
    let mut sum_omega = 0.0;
    for v in 0..nv {
        let site = &data.sites[v];
        if !data.unknown[v] {
            c.state_residual = c.state_residual.max(u[v].abs());
            c.adjoint_residual = c.adjoint_residual.max(phi[v].abs());
            c.bound_violation = c.bound_violation.max(mu[v].abs());
            continue;
        }
        let w = data.w_e[v] * mass[v];
        c.state_residual = c.state_residual.max((l_u[v] - mu[v] + omega[v]).abs());
        c.adjoint_residual = c.adjoint_residual.max((l_phi[v] + w * u[v]).abs());
        c.cone_residual = c
            .cone_residual
            .max((mu[v] - site.prox(phi[v] + mu[v])).abs());
        c.bound_violation = c.bound_violation.max(site.bound_violation(mu[v]));
        c.dual_infeasibility = c.dual_infeasibility.max((phi_hat[v] - phi[v]).abs());
        sum_mu += mu[v];
        sum_omega += omega[v];

        distortion += w * u[v] * u[v];
        cone_cost += site.cost(mu[v]);
        if w > 0.0 {
            dual_quadratic += l_phi_hat[v] * l_phi_hat[v] / w;
            let u_tilde = -l_phi_hat[v] / w;
            mismatch += w * (u[v] - u_tilde) * (u[v] - u_tilde);
        }
        dual_linear += phi_hat[v] * omega[v];
        let conjugate = site.conjugate(phi_hat[v]);
        dual_conjugate += conjugate;
        slack += site.cost(mu[v]) + conjugate - phi_hat[v] * mu[v];
    }
    if closed {
        c.gauss_bonnet_residual = (sum_mu - sum_omega).abs();
    }
    c.primal_objective = 0.5 * distortion + cone_cost;
    c.dual_objective = -0.5 * dual_quadratic + dual_linear - dual_conjugate;
    c.duality_gap = c.primal_objective - c.dual_objective;
    c.complementarity_gap = 0.5 * mismatch + slack;
    Ok(c)
}
