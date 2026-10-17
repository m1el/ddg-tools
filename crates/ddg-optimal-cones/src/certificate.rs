//! The optimality certificate: residuals of the optimality system and a duality gap that
//! bounds the suboptimality, computed from `(u, φ, µ)` and the problem data alone.
//!
//! Nothing here uses the solver's state: the Laplacian is applied face by face from the corner
//! cotangents (not through the assembled matrix or its interior block), the cone equation is
//! evaluated from the problem's weights and bounds (not from the solver's active sets), and the
//! dual point is rebuilt from `φ`.

use ddg_core::{Geometry, TriMesh};

use crate::error::ConeError;
use crate::problem::{BoundaryTreatment, ConeProblem, Instance};

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
/// (machine-checked: `optimality-system-sufficient` in `lean/proofs`). The first six fields
/// measure how far the triple is from satisfying it.
///
/// # The duality gap, and exactly what is computed
///
/// The gap makes "how far" quantitative. From the triple the following are formed, with
/// sums over the unknown (interior, or all on a closed mesh) vertices:
///
/// ```text
/// r_i  = (L u)_i − µ_i + Ω_i                         the residual of the state equation
/// φ̂_i  = φ_i clipped to the domain of g_i*:          φ̂_i ≤ t_i if hi_i = +∞,  φ̂_i ≥ −t_i if lo_i = −∞
/// ũ_i  = −(L φ̂)_i / w_i                              where w_i > 0
/// J    = ½ Σ w_i u_i² + Σ t_i |µ_i|                   primal_objective
/// D    = −½ Σ_{w_i>0} (L φ̂)_i² / w_i + Σ φ̂_i Ω_i − Σ g_i*(φ̂_i)        dual_objective
/// ```
///
/// (`g_i*(p) = hi_i·max(0, p − t_i) + lo_i·min(0, p + t_i)`, a term with an infinite bound
/// absent; `a·p` at a vertex whose angle is fixed to `a`, where also `t_i |µ_i|` is zero.) The
/// dual point `(ũ, φ̂)` is built so that the adjoint equation `L φ̂ = −w ũ` holds **by
/// construction** and `φ̂` is dual feasible **exactly**; nothing about it relies on the
/// residuals of the triple. For any numbers `u, φ, µ` whatsoever (zero on the boundary) the
/// following is an algebraic identity (`duality-gap-inexact`):
///
/// ```text
/// J − D  =  ½ Σ_{w_i>0} w_i (u_i − ũ_i)²  +  Σ (t_i|µ_i| + g_i*(φ̂_i) − φ̂_i µ_i)
///           +  Σ_{w_i=0} (L φ̂)_i u_i  −  Σ φ̂_i r_i
/// ```
///
/// `duality_gap` is the left side. `complementarity_gap` is the first line of the right side;
/// each of its summands is nonnegative when `µ` respects its bounds (Fenchel–Young).
/// `gap_uncertainty` is `Σ_{w_i=0} |(L φ̂)_i u_i| + Σ |φ̂_i r_i|`, a bound for the magnitude of
/// the second line.
///
/// # What this proves
///
/// The statements are machine-checked in `lean/proofs` (`DdgProofs/Certificate.lean`,
/// `DdgProofs/Duality.lean`) for exact real arithmetic; `φ̂` is dual feasible by construction.
///
/// * The identity above is `objective_sub_dualObjectiveOfAdjoint`.
/// * **Lower bound** (`dual_lower_bound_of_clip`): if `(L φ̂)_i = 0` wherever `w_i = 0` (no
///   condition if all `w_i > 0`), then `D ≤ objective(v, ν)` for every feasible pair `(v, ν)`,
///   whatever the residuals of the triple are. So `D` is a lower bound of the optimal value.
/// * **Backward statement** (`backward_certificate`): if in addition `µ` respects its bounds
///   (`bound_violation = 0`), then `(u, µ)` is exactly feasible for the problem whose angle
///   defects are `Ω − r` instead of `Ω`, and its objective exceeds the optimum *of that
///   problem* by at most `complementarity_gap ≥ 0`. `|r_i| ≤ state_residual`.
/// * **Exactly feasible triples**: if `r = 0`, then
///   `objective(u, µ) − optimum ≤ duality_gap = complementarity_gap`
///   (`adjoint_gap_certificate`).
///
/// What it does **not** prove:
///
/// * Nothing bounds the suboptimality of the returned `(u, µ)` in the *given* problem,
///   because that pair satisfies the state equation only up to `state_residual` and a pair
///   that is not feasible can have an objective *below* the optimum. `duality_gap` can
///   therefore be negative; it is `complementarity_gap` (which cannot) shifted by at most
///   `gap_uncertainty`. Read the gap as "`complementarity_gap`, give or take
///   `gap_uncertainty`". (An enclosure of the optimal value of the given problem,
///   `D ≤ optimum ≤ J + Σ t_i |r_i|`, holds if `µ + r` is admissible at every vertex,
///   `certificate_enclosure`; that fails as soon as a forbidden or prescribed vertex has a
///   nonzero residual, so the code does not evaluate it.)
/// * Where `w_i = 0` the dual problem requires `(L φ̂)_i = 0` exactly, which cannot be
///   enforced by construction. The vertex contributes nothing to `D`, and its violation
///   `(L φ̂)_i` is in `adjoint_residual` and, times `u_i`, in `gap_uncertainty`. With some
///   `w_i = 0`, `D` is a lower bound of the optimum only up to `Σ_{w_i=0} (L φ̂)_i v_i`,
///   which depends on the competitor `v` (`dualObjectiveOfAdjoint_add_le`).
/// * All of this is about exact arithmetic; the values are computed in `f64`.
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
    /// boundary. On a mesh with several components: the largest such violation over the
    /// closed components (each must satisfy Gauss–Bonnet on its own), zero if none is closed.
    pub gauss_bonnet_residual: f64,
    /// The objective `J = ½ Σ wE m u² + λ Σ wR |µ|` of `(u, µ)`.
    pub primal_objective: f64,
    /// The dual value `D` of `(ũ, φ̂)`: a lower bound of the optimal objective (exact
    /// arithmetic, all `wE_i m_i > 0`).
    pub dual_objective: f64,
    /// `primal_objective − dual_objective`. A bound on the suboptimality of `(u, µ)` if that
    /// pair is exactly feasible; otherwise see `gap_uncertainty`.
    pub duality_gap: f64,
    /// `½ Σ_{w_i>0} w_i (u_i − ũ_i)² + Σ (t_i|µ_i| + g_i*(φ̂_i) − φ̂_i µ_i)`: the gap as a sum
    /// of terms that are nonnegative when `µ` respects its bounds. It differs from
    /// `duality_gap` by at most `gap_uncertainty` (and rounding).
    pub complementarity_gap: f64,
    /// `Σ_{w_i=0} |(L φ̂)_i u_i| + Σ |φ̂_i r_i|` with `r` the residual of the state equation:
    /// how much of `duality_gap` is due to the triple not satisfying the state equation
    /// exactly (and to dual infeasibility at vertices of zero distortion weight) rather than
    /// to suboptimality.
    pub gap_uncertainty: f64,
    /// The repair cost `ε = λ Σ wR_i |r_i|` of absorbing the state residual `r = L u − µ + Ω`
    /// into the cone angles, over the penalized unknown vertices (a vertex with a fixed angle
    /// and a free cone site cost nothing). See [`Certificate::optimum_enclosure`].
    pub repair_cost: f64,
    /// How far `µ + r` is from being admissible: the largest, over the unknown vertices, distance
    /// of `µ_i + r_i` from `[lo_i, hi_i]` (from the fixed angle at a forbidden or prescribed
    /// vertex, from `0` at a forbidden one), and `max(|u_i|, |µ_i|)` over the vertices that are
    /// not unknowns (boundary vertices with `u = 0`), which must be exactly `0`. It is `0` when
    /// `µ + r` is admissible everywhere. With a forbidden or prescribed vertex and a nonzero
    /// residual there it is positive, and nothing is proved about the optimum (see
    /// [`Certificate::optimum_enclosure`]).
    pub repair_bound_violation: f64,
    /// Whether `wE_i m_i > 0` at every unknown vertex, which the lower bound `D ≤ optimum`
    /// needs.
    pub weights_positive: bool,
}

impl Certificate {
    /// An enclosure `(lower, upper)` of the optimal objective value of the **given** problem, or
    /// `None` when none is proved.
    ///
    /// * `lower = dual_objective = D`. Proved (`dual_lower_bound_of_clip`) when
    ///   [`weights_positive`](Certificate::weights_positive): `D` is at most the objective of every
    ///   feasible pair, whatever the residuals of the triple are.
    /// * `upper = primal_objective + repair_cost`. Proved (`certificate_enclosure`) when
    ///   `repair_bound_violation == 0`: then `µ + r` is admissible, `(u, µ + r)` is exactly
    ///   feasible, and its objective is at most `J(u, µ) + λ Σ wR_i |r_i|`.
    ///
    /// Both statements are about exact real arithmetic applied to the `f64` values of the
    /// triple and of `Ω`, `L`, `m`; rounding in evaluating `D`, `J`, `r` is not covered.
    ///
    /// When is `repair_bound_violation` zero? Always at a penalized vertex without bounds. At a
    /// vertex with a bound it needs `µ_i + r_i` inside the bounds, and at a forbidden or
    /// prescribed vertex it needs `r_i = 0` exactly (the angle must stay the fixed one). So with a
    /// forbidden or prescribed vertex and a nonzero residual there, or an angle at its bound
    /// with a residual pushing it out, the result is `None`: `(u, µ + r)` is then not feasible,
    /// no feasible point with this `u` exists, and nothing bounds the optimum from above. (The
    /// lower bound is still valid; see [`Certificate::dual_objective`].) Not finite values give
    /// `None` too.
    ///
    /// When it is `Some((lower, upper))`, the optimal value lies in it, the pair `(u, µ + r)`
    /// is feasible, and its objective exceeds the optimum by at most
    /// `duality_gap + repair_cost` (`certificate_of_computed`).
    ///
    /// [formula: duality-gap-certificate]
    pub fn optimum_enclosure(&self) -> Option<(f64, f64)> {
        // Written so that a NaN fails.
        let finite = self.dual_objective.is_finite()
            && self.primal_objective.is_finite()
            && self.repair_cost.is_finite()
            && self.duality_gap.is_finite();
        if finite && self.weights_positive && self.repair_bound_violation == 0.0 {
            Some((
                self.dual_objective,
                self.primal_objective + self.repair_cost,
            ))
        } else {
            None
        }
    }

    /// The proved bound `duality_gap + repair_cost` on how much the objective of the feasible
    /// pair `(u, µ + r)` exceeds the optimum, under the conditions of
    /// [`Certificate::optimum_enclosure`] (`None` otherwise).
    pub fn suboptimality_bound(&self) -> Option<f64> {
        self.optimum_enclosure()
            .map(|_| self.duality_gap + self.repair_cost)
    }

    /// Whether every residual is at most `tolerance` and both forms of the gap and the gap
    /// uncertainty are at most `tolerance · max(1, |primal_objective|)`.
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
            && self.gap_uncertainty <= gap_tolerance
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
/// The mesh may have several connected components; the certificate is then that of the whole
/// problem (the sums run over all components, the residuals are maxima over all vertices, and
/// the Gauss–Bonnet residual is the largest over the closed components).
///
/// This is the check [`solve`](crate::solve) runs on its own result; it can equally be used on
/// a triple from any other source.
///
/// # Errors
/// The errors of the problem's validation, and [`ConeError::WrongLength`].
///
/// [formula: optimality-system-sufficient]
/// [formula: duality-gap-certificate]
/// [formula: duality-gap-inexact]
/// [formula: cone-fixed-point-bounded]
pub fn certify(
    problem: &ConeProblem<'_>,
    log_scale: &[f64],
    adjoint: &[f64],
    measure: &[f64],
) -> Result<Certificate, ConeError> {
    certify_with_boundary(
        problem,
        BoundaryTreatment::Fixed,
        log_scale,
        adjoint,
        measure,
    )
}

/// [`certify`] for a boundary treatment: the certificate of a triple for the problem that
/// [`solve_with_boundary`](crate::solve_with_boundary) solves. With a treatment other than
/// `Fixed` the boundary vertices are unknowns: `measure` holds the curvature `ν` of the
/// flattened boundary at them, `k_i` (the boundary geodesic curvature) plays the role of the
/// angle defect there, and all residuals and sums run over all vertices; the Gauss–Bonnet
/// residual is `|Σ_I µ + Σ_B ν − 2πχ|`.
///
/// # Errors
/// As [`certify`].
///
/// [formula: optimality-system-sufficient]
/// [formula: duality-gap-certificate]
/// [formula: duality-gap-inexact]
/// [formula: cone-fixed-point-bounded]
/// [formula: boundary-optimality-system]
pub fn certify_with_boundary(
    problem: &ConeProblem<'_>,
    boundary: BoundaryTreatment,
    log_scale: &[f64],
    adjoint: &[f64],
    measure: &[f64],
) -> Result<Certificate, ConeError> {
    let instance = Instance::new(problem, boundary)?;
    let nv = problem.mesh.n_vertices();
    check_len("log_scale", log_scale, nv)?;
    check_len("adjoint", adjoint, nv)?;
    check_len("measure", measure, nv)?;
    Ok(certify_instance(&instance, log_scale, adjoint, measure))
}

/// [`certify_with_boundary`] on validated data: the whole problem, or one component of it
/// (then the certificate is that of the component's part of the problem: its share of the
/// objective, with the masses of the whole mesh). The vectors have one entry per vertex of
/// `instance.mesh`.
///
/// [formula: duality-gap-inexact]
/// [formula: component-separation]
pub(crate) fn certify_instance(
    instance: &Instance<'_>,
    log_scale: &[f64],
    adjoint: &[f64],
    measure: &[f64],
) -> Certificate {
    let data = &instance.data;
    let mesh = instance.mesh;
    let geometry = instance.geometry;
    let nv = mesh.n_vertices();
    let (u, phi, mu) = (log_scale, adjoint, measure);
    let omega = &instance.omega;
    let mass = &instance.mass;
    // Gauss–Bonnet holds per closed component: the sums of µ and Ω, and whether every vertex is
    // an unknown, per component.
    let component: Option<Vec<u32>> = (mesh.n_components() > 1).then(|| mesh.vertex_components());
    let n_components = mesh.n_components();
    let mut sums = vec![(0.0f64, 0.0f64, true); n_components];

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
        gap_uncertainty: 0.0,
        repair_cost: 0.0,
        repair_bound_violation: 0.0,
        weights_positive: true,
    };
    let mut distortion = 0.0;
    let mut cone_cost = 0.0;
    let mut dual_quadratic = 0.0;
    let mut dual_linear = 0.0;
    let mut dual_conjugate = 0.0;
    let mut mismatch = 0.0;
    let mut slack = 0.0;
    for v in 0..nv {
        let site = &data.sites[v];
        let sum = &mut sums[component.as_ref().map_or(0, |c| c[v] as usize)];
        if !data.unknown[v] {
            sum.2 = false;
            c.state_residual = c.state_residual.max(u[v].abs());
            c.adjoint_residual = c.adjoint_residual.max(phi[v].abs());
            c.bound_violation = c.bound_violation.max(mu[v].abs());
            c.repair_bound_violation = c.repair_bound_violation.max(u[v].abs()).max(mu[v].abs());
            continue;
        }
        let w = data.w_e[v] * mass[v];
        let state = l_u[v] - mu[v] + omega[v];
        c.state_residual = c.state_residual.max(state.abs());
        c.gap_uncertainty += (phi_hat[v] * state).abs();
        c.repair_cost += site.cost(state);
        c.repair_bound_violation = c
            .repair_bound_violation
            .max(site.bound_violation(mu[v] + state));
        c.weights_positive &= w > 0.0;
        c.adjoint_residual = c.adjoint_residual.max((l_phi[v] + w * u[v]).abs());
        c.cone_residual = c
            .cone_residual
            .max((mu[v] - site.prox(phi[v] + mu[v])).abs());
        c.bound_violation = c.bound_violation.max(site.bound_violation(mu[v]));
        c.dual_infeasibility = c.dual_infeasibility.max((phi_hat[v] - phi[v]).abs());
        sum.0 += mu[v];
        sum.1 += omega[v];

        distortion += w * u[v] * u[v];
        cone_cost += site.cost(mu[v]);
        if w > 0.0 {
            dual_quadratic += l_phi_hat[v] * l_phi_hat[v] / w;
            let u_tilde = -l_phi_hat[v] / w;
            mismatch += w * (u[v] - u_tilde) * (u[v] - u_tilde);
        } else {
            c.gap_uncertainty += (l_phi_hat[v] * u[v]).abs();
        }
        dual_linear += phi_hat[v] * omega[v];
        let conjugate = site.conjugate(phi_hat[v]);
        dual_conjugate += conjugate;
        slack += site.cost(mu[v]) + conjugate - phi_hat[v] * mu[v];
    }
    for &(sum_mu, sum_omega, closed) in &sums {
        if closed {
            c.gauss_bonnet_residual = c.gauss_bonnet_residual.max((sum_mu - sum_omega).abs());
        }
    }
    c.primal_objective = 0.5 * distortion + cone_cost;
    c.dual_objective = -0.5 * dual_quadratic + dual_linear - dual_conjugate;
    c.duality_gap = c.primal_objective - c.dual_objective;
    c.complementarity_gap = 0.5 * mismatch + slack;
    c
}
