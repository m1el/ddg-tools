//! The solver: options, result, report, and the two-stage driver.

use std::time::Instant;

use ddg_core::TriMesh;
use ddg_core::io::Cone;
use ddg_core::sparse::Parallelism;

use crate::certificate::{Certificate, apply_laplacian, certify_instance};
use crate::error::ConeError;
use crate::evaluate::{LaplaceSolver, cones_to_measure, log_scale_of_measure};
use crate::extract::{ExtractOptions, ExtractionReport, extract_cones};
use crate::factor_once::{
    self, FactorOnceOptions, FactorOnceReport, FallbackReason, MethodUsed, SolveMethod, Unsupported,
};
use crate::linsys::{LinearSolver, LinearSolverReport, NewtonLinearSolver, SaddleSolver};
use crate::newton::{
    Iterate, NewtonOutcome, StageStop, adjoint_constant_interval, solve_exact, solve_regularized,
};
use crate::problem::{BoundaryTreatment, ConeProblem, Discrete, Instance, split_pieces};
use ddg_core::operators::boundary_geodesic_curvatures;

/// The quantity whose square the regularized problems of Stage I penalize.
///
/// Stage I only initializes Stage II, so this choice affects the number of Newton steps and
/// nothing else: both variants lead to the same solution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RegularizationMetric {
    /// `(γ/2) Σ µ_i²`: the regularized cone angle is `µ_i = clamp(P_t(φ_i)/γ)`, the update the
    /// paper prints (Eqs. 21–22), which has no mass factor. With area-normalized masses
    /// (`m_i ≈ 1/n`) this is about `n` times stiffer than `Mass` at the same `γ`, so the last
    /// regularized solution has its cones on single vertices and the exact iteration that
    /// follows typically needs one step. The default, because it needs fewer Newton steps
    /// (measured on the dataset, see the crate README).
    #[default]
    Unweighted,
    /// `(γ/2) Σ µ_i² / m_i`: the squared `L²` norm of the curvature *density* `µ_i / m_i`,
    /// which is what paper Eq. 12 discretizes to. The regularized cone angle is
    /// `µ_i = clamp((m_i/γ)·P_t(φ_i))`.
    Mass,
}

/// Options of [`solve`]. `SolveOptions::default()` is the configuration the crate is tested
/// with.
#[derive(Debug, Clone)]
pub struct SolveOptions<'a> {
    /// Stage I: the regularization parameters `γ`, strictly decreasing and positive. One
    /// regularized problem is solved for each, starting from the solution of the previous one.
    /// Empty skips Stage I, which leaves the plain semismooth Newton method started from zero.
    /// Default: `1, 10⁻¹, …, 10⁻¹⁰` (paper Fig. 12 with `N = 10`;
    /// [`SolveOptions::paper_regularization`]). The shorter `10⁻², 10⁻⁶, 10⁻¹⁰`
    /// ([`SolveOptions::short_regularization`]) reaches the same optimum with 1.2 to 3.9 times
    /// fewer factorizations on the unbounded dataset problems, but needs up to 2.5 times more
    /// on problems with positive-only capped angles (what the UV pipeline asks for), which is
    /// why it is not the default (README, "The regularization schedule").
    pub regularization: Vec<f64>,
    /// What the regularization measures. Default: [`RegularizationMetric::Unweighted`].
    pub regularization_metric: RegularizationMetric,
    /// Limit on the Newton steps of one regularized problem. A regularized problem that hits
    /// the limit is not an error (its result only initializes what follows); the stage report
    /// records it. Default 50.
    pub max_regularized_iterations: usize,
    /// Limit on the Newton steps of one run of the exact iteration (Stage II). Default 50.
    pub max_exact_iterations: usize,
    /// Safeguard: how many extra regularized problems the solver may solve when the paper's
    /// two stages do not finish. If the exact iteration cycles, hits its limit or meets a
    /// singular system, a *proximal step* follows (a regularized problem whose penalty is
    /// centred at the current cone angles instead of zero) and the exact iteration is retried.
    /// Proximal steps converge to a minimizer whether or not it is unique; once a step no
    /// longer moves the cone angles the iterate is accepted as it is. If a regularized problem
    /// itself hits its iteration limit, it is continued before anything else. `0` disables
    /// the safeguard. Default 12.
    pub max_safeguard_rounds: usize,
    /// The exact iteration also stops when the cone residual
    /// `max_i |µ_i − clamp(P_t(φ_i + µ_i))|` of a solved step is at most this (paper
    /// Sec. 5.3.3, "norm of the residual below a small tolerance"). This only matters when two
    /// pieces tie at rounding level. Cone angles of at most this magnitude are set to exactly
    /// zero at the end ([`SolveReport::zeroed_angles`]). Default `10⁻¹²`.
    pub newton_tolerance: f64,
    /// The result is reported as [`Status::Converged`] iff its [`Certificate`] holds with this
    /// tolerance ([`Certificate::holds`]). Default `10⁻⁹`.
    pub certificate_tolerance: f64,
    /// A vertex is a cone if `|µ_i|` exceeds this. Default `10⁻¹²` (paper Sec. 6.1).
    pub cone_threshold: f64,
    /// Whether to replace each cluster of cones by a single cone (paper Sec. 6.1). A cluster
    /// is a set of cones **of the same sign** connected through mesh edges between its
    /// members; a positive and a negative cone are never merged, because summing them cancels
    /// cone angle that the optimum needs. Default `false`: the returned cones are the optimal
    /// measure. What merging cost in distortion is reported
    /// ([`ExtractionReport::distortion_change`]; `distortion_optimal` against `distortion`).
    pub merge_cone_clusters: bool,
    /// Clusters with more cones than this are never merged. The paper merges "tiny clusters";
    /// for very small `λ`, and on convex parts of closed meshes, the cones cover whole
    /// regions, and merging such a region into one vertex would be meaningless. Default 8.
    pub max_merge_cluster_size: usize,
    /// Clusters whose diameter exceeds this are never merged. The diameter is the largest
    /// distance between two cones of the cluster along mesh edges inside the cluster, in the
    /// length unit of the geometry ([`MergedCluster::diameter`](crate::MergedCluster)); a few
    /// edge lengths is what "tiny" means in the paper. Default `f64::INFINITY`: only the
    /// number of cones limits merging.
    pub max_merge_cluster_diameter: f64,
    /// How the Newton systems are solved. Default [`LinearSolver::Auto`].
    pub linear_solver: LinearSolver,
    /// Where factorizations and solves run. Default: on the calling thread. Threads gain
    /// little for these factorizations (see the `ddg-core` README).
    pub parallelism: Parallelism<'a>,
    /// Which method solves the problem. Default [`SolveMethod::Auto`]: the factor-once method
    /// where it applies and is fast, the reference method otherwise. The report says which one
    /// produced the result ([`SolveReport::method`]).
    pub method: SolveMethod,
    /// Options of the factor-once method (thresholds for giving up, tolerance).
    pub factor_once: FactorOnceOptions,
}

impl Default for SolveOptions<'_> {
    fn default() -> Self {
        SolveOptions {
            regularization: SolveOptions::paper_regularization(),
            regularization_metric: RegularizationMetric::default(),
            max_regularized_iterations: 50,
            max_exact_iterations: 50,
            max_safeguard_rounds: 12,
            newton_tolerance: 1e-12,
            certificate_tolerance: 1e-9,
            cone_threshold: 1e-12,
            merge_cone_clusters: false,
            max_merge_cluster_size: 8,
            max_merge_cluster_diameter: f64::INFINITY,
            linear_solver: LinearSolver::default(),
            parallelism: Parallelism::Sequential,
            method: SolveMethod::default(),
            factor_once: FactorOnceOptions::default(),
        }
    }
}

impl SolveOptions<'_> {
    /// The paper's schedule of Stage I (Fig. 12 with `N = 10`): `γ = 1, 10⁻¹, …, 10⁻¹⁰`. The
    /// default.
    pub fn paper_regularization() -> Vec<f64> {
        (0..=10).map(|k| 10f64.powi(-k)).collect()
    }

    /// The short schedule `γ = 10⁻², 10⁻⁶, 10⁻¹⁰`: the same optimum with fewer factorizations
    /// on problems without bounds, more on problems with positive-only capped angles (see
    /// [`SolveOptions::regularization`]).
    pub fn short_regularization() -> Vec<f64> {
        vec![1e-2, 1e-6, 1e-10]
    }

    pub(crate) fn validate(&self) -> Result<(), ConeError> {
        let bad = |reason: String| Err(ConeError::InvalidOptions { reason });
        for (i, &gamma) in self.regularization.iter().enumerate() {
            if !(gamma.is_finite() && gamma > 0.0) {
                return bad(format!(
                    "regularization[{i}] is {gamma}; every entry must be positive and finite"
                ));
            }
            // Both entries are finite here, so this is the negation of "strictly decreasing".
            if i > 0 && gamma >= self.regularization[i - 1] {
                return bad(format!(
                    "regularization must be strictly decreasing, but entry {i} ({gamma}) is not \
                     below entry {} ({})",
                    i - 1,
                    self.regularization[i - 1]
                ));
            }
        }
        // `+∞` is allowed; NaN and negative values are not.
        if self.max_merge_cluster_diameter.is_nan() || self.max_merge_cluster_diameter < 0.0 {
            return bad(format!(
                "max_merge_cluster_diameter is {}; it must be non-negative",
                self.max_merge_cluster_diameter
            ));
        }
        if self.max_exact_iterations == 0 {
            return bad("max_exact_iterations must be at least 1".to_string());
        }
        if !self.regularization.is_empty() && self.max_regularized_iterations == 0 {
            return bad("max_regularized_iterations must be at least 1".to_string());
        }
        for (name, value) in [("factor_once.tolerance", self.factor_once.tolerance)] {
            if !(value.is_finite() && value >= 0.0) {
                return bad(format!(
                    "{name} is {value}; it must be finite and non-negative"
                ));
            }
        }
        if self.factor_once.max_rounds == 0 {
            return bad("factor_once.max_rounds must be at least 1".to_string());
        }
        for (name, value) in [
            ("newton_tolerance", self.newton_tolerance),
            ("certificate_tolerance", self.certificate_tolerance),
            ("cone_threshold", self.cone_threshold),
        ] {
            if !(value.is_finite() && value >= 0.0) {
                return bad(format!(
                    "{name} is {value}; it must be finite and non-negative"
                ));
            }
        }
        Ok(())
    }
}

/// Which method produced a result and why.
#[derive(Debug, Clone, PartialEq)]
pub struct MethodReport {
    /// What [`SolveOptions::method`] asked for.
    pub requested: SolveMethod,
    /// The method that produced the returned result. With [`SolveMethod::FactorOnce`] this is
    /// always the factor-once method, also when it gave up (then `fallback` says why and the
    /// status is not [`Status::Converged`]).
    pub used: MethodUsed,
    /// Why the factor-once method did not produce a certified result, if it was tried or
    /// excluded and did not. `None` when it succeeded, and when the reference method was
    /// requested. With `requested == Auto` and `used == Reference` it is the reason for the
    /// fallback.
    pub fallback: Option<FallbackReason>,
    /// What the factor-once method did, if it ran (also when it gave up).
    pub factor_once: Option<FactorOnceReport>,
}

/// Whether the result is a certified optimum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// The optimality certificate holds within
    /// [`certificate_tolerance`](SolveOptions::certificate_tolerance): the returned measure is
    /// optimal up to the certified duality gap.
    Converged,
    /// The iterations ended at their limits without a point whose certificate holds. The
    /// returned measure is the best available iterate and **is not** known to be optimal; its
    /// [`Certificate`] says how far off it is.
    NotConverged,
    /// A Newton system could not be solved accurately by any available factorization, and the
    /// certificate of the last iterate does not hold. See [`SolveReport::message`].
    LinearSolverFailed,
}

/// Which problem a stage solved.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StageKind {
    /// Stage I: the regularized problem with parameter `gamma`.
    Regularized {
        /// The regularization parameter `γ`.
        gamma: f64,
    },
    /// Stage II: the exact problem.
    Exact,
    /// Safeguard: a proximal step, i.e. the regularized problem with parameter `gamma` whose
    /// penalty is centred at the current cone angles instead of zero.
    Proximal {
        /// The regularization parameter `γ`.
        gamma: f64,
    },
}

/// What one stage of the solve did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StageReport {
    /// Which problem.
    pub kind: StageKind,
    /// Whether the stage ran in a safeguard round, i.e. after the paper's two stages (Stage I,
    /// then one run of the exact iteration) had failed to finish. A safeguard stage of kind
    /// [`StageKind::Regularized`] is a regularized problem centred at zero that the safeguard
    /// solved (or continued); it is not part of Stage I.
    pub safeguard: bool,
    /// Newton steps. Each is one factorization of the Newton matrix and one linear solve
    /// (plus refinement).
    pub newton_iterations: usize,
    /// Steps that the line search shortened (regularized and proximal stages only).
    pub damped_steps: usize,
    /// Closed meshes: how often the constant of `φ` was shifted to restore Gauss–Bonnet.
    pub constant_shifts: usize,
    /// Why the stage stopped.
    pub stop: StageStop,
    /// Wall time of the stage in seconds.
    pub seconds: f64,
}

/// Wall times of a solve, in seconds.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Timings {
    /// Validation and assembly of the operators.
    pub setup: f64,
    /// Stage I.
    pub regularized: f64,
    /// Stage II including safeguard rounds.
    pub exact: f64,
    /// The factor-once method (also when it gave up and the reference method followed).
    pub factor_once: f64,
    /// The certificate.
    pub certificate: f64,
    /// Cone extraction, including the recomputation of `u` when the cones changed.
    pub extraction: f64,
    /// Everything.
    pub total: f64,
}

/// What a solve did and how good its result is.
#[derive(Debug, Clone, PartialEq)]
pub struct SolveReport {
    /// Whether the result is a certified optimum.
    pub status: Status,
    /// The stages in the order they ran.
    pub stages: Vec<StageReport>,
    /// Total Newton steps of Stage I: the regularized problems of
    /// [`SolveOptions::regularization`], and nothing else (zero if that list is empty).
    pub regularized_iterations: usize,
    /// Total Newton steps of the exact iteration (all runs, including those of safeguard
    /// rounds).
    pub exact_iterations: usize,
    /// Number of safeguard rounds that were needed (0 when the paper's two stages sufficed).
    pub safeguard_rounds: usize,
    /// Total Newton steps of the regularized and proximal problems that the safeguard solved
    /// (the stages with [`StageReport::safeguard`] set that are not exact). Every Newton step
    /// of a solve is counted in exactly one of `regularized_iterations`, `exact_iterations`
    /// and `safeguard_iterations`.
    pub safeguard_iterations: usize,
    /// Factorization and solve counts and times of the linear solver.
    pub linear_solver: LinearSolverReport,
    /// Wall times.
    pub timings: Timings,
    /// The largest amount by which a cone angle was moved onto its bound at the very end.
    /// Zero whenever the last Newton step left the active sets unchanged; otherwise at most
    /// of the order of [`SolveOptions::newton_tolerance`]. The certificate is computed after
    /// this correction.
    pub bound_correction: f64,
    /// Number of cone angles of magnitude at most [`SolveOptions::newton_tolerance`] that were
    /// set to exactly zero at the very end. An active vertex carries such an angle when it
    /// absorbs the rounding error of the data (for instance of `Σ Ω` on a closed mesh). The
    /// certificate is computed after this.
    pub zeroed_angles: usize,
    /// The largest magnitude among them (zero if there were none).
    pub largest_zeroed_angle: f64,
    /// The optimality certificate of `(u, φ, µ)` as returned.
    pub certificate: Certificate,
    /// Which method produced the result, and why.
    pub method: MethodReport,
    /// Explanation when something did not go to plan (a linear solver failure, or why the
    /// result is not converged). On a mesh with several components: the messages of the
    /// components, each prefixed with `component <c>: `.
    pub message: Option<String>,
    /// One report per connected component of the mesh, in the order of
    /// `TriMesh::face_components` (`TriMesh::vertex_components` says which vertex belongs to
    /// which), when the mesh has more than one; **empty for a connected mesh**, whose report is
    /// this one. With several components the other fields of this report describe the whole
    /// solve: `certificate` is the certificate of the whole problem and alone decides `status`;
    /// `stages` are those of the components one after another; the counts and the times of the
    /// phases are sums over the components (`timings.total` is wall time);
    /// `bound_correction` and `largest_zeroed_angle` are maxima; `method.used` is
    /// [`MethodUsed::FactorOnce`] only if every component with unknowns used it, `method.fallback`
    /// is the first component's reason, `method.factor_once` is `None` (see the components).
    pub components: Vec<ComponentReport>,
}

/// What happened on one connected component of a disconnected mesh
/// ([`SolveReport::components`]).
///
/// The component is solved as its part of the whole problem: masses divided by the total area
/// of the whole mesh, the given `λ`. Its `objective`, `distortion_optimal`, `cone_cost` and its
/// certificate are about that part (they sum, over the components, to those of the whole).
#[derive(Debug, Clone, PartialEq)]
pub struct ComponentReport {
    /// The number of the component (`TriMesh::face_components`).
    pub component: u32,
    /// Its vertices.
    pub n_vertices: usize,
    /// Its faces.
    pub n_faces: usize,
    /// Its unknowns (interior vertices; all vertices if it is closed or the boundary is
    /// optimized). `0` for a component with nothing to solve, such as a single triangle with
    /// `u = 0` on its boundary: it gets `u = µ = 0` and no stages.
    pub n_unknowns: usize,
    /// Whether every vertex of the component is an unknown (a closed component, or a boundary
    /// treatment other than `Fixed`): then `u` has its own additive constant on it, and its
    /// angles (and `ν`) sum to its own `2πχ`.
    pub closed: bool,
    /// Its share of the total area, `A_c / A` (the sum of its normalized masses).
    pub area_fraction: f64,
    /// The `λ` for which the component solved alone, with its masses divided by its own area,
    /// has the same `u` and `µ`: `λ · A / A_c` (its adjoint is then `A / A_c` times the one
    /// returned here).
    pub equivalent_lambda: f64,
    /// Its part of the objective.
    pub objective: f64,
    /// Its part of the distortion of the optimum.
    pub distortion_optimal: f64,
    /// Its part of the cone cost of the optimum.
    pub cone_cost: f64,
    /// The interval of constants that may be added to the adjoint on this component alone
    /// (`(0, 0)` unless it is closed).
    pub adjoint_shift: (f64, f64),
    /// The solve of the component: status and certificate of its part of the problem, stages,
    /// iterations, method, times (`timings.total` is its wall time; setup and extraction are
    /// those of the whole). Its `components` is empty.
    pub report: SolveReport,
}

/// The result of [`solve`]. It holds two things that must not be confused:
///
/// * **The optimum** of the canonical problem: `measure`, `log_scale_optimal`, `adjoint`,
///   `objective`, `distortion_optimal`, `cone_cost`. The [`Certificate`] in `report` is
///   computed for exactly this triple, as returned: after the rounding-level clean-up of the
///   cone angles ([`SolveReport::bound_correction`], [`SolveReport::zeroed_angles`]) and after
///   the adjoint constant of a closed mesh was chosen. `certify(problem, &log_scale_optimal,
///   &adjoint, &measure)` reproduces it bit for bit. The scale factor is not recomputed after
///   the clean-up, so what the clean-up changed is visible in the certificate's
///   `state_residual`.
/// * **The extracted cones**: `cones`, `log_scale`, `distortion`, `total_cone_angle`. They
///   are the optimum unless [`ExtractionReport::cones_differ_from_measure`] is set (a
///   threshold dropped small angles, or clusters were merged). Then they are deliberately
///   *not* optimal, nothing certifies them, and `extraction` lists every difference.
///
/// Every per-vertex vector has one entry per mesh vertex and is zero at boundary vertices.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// The optimal measure `µ`: the cone angle at every vertex, exactly zero where there is no
    /// cone.
    pub measure: Vec<f64>,
    /// The log scale factor `u` of the optimal measure: `(L u)_i = µ_i − Ω_i`. Zero at the
    /// boundary vertices for [`BoundaryTreatment::Fixed`]; with an optimized boundary it is
    /// returned there too.
    pub log_scale_optimal: Vec<f64>,
    /// The geodesic curvature `ν_i` of the flattened boundary at every boundary vertex, zero at
    /// interior vertices: `ν_i = (L u)_i + k_i` with `k_i = π −` (angle sum) the boundary
    /// curvature of the mesh. For the treatments other than `Fixed` it is the optimized
    /// variable, a copy of `measure` at the boundary vertices; for `Fixed` it is computed from
    /// `u` (which vanishes on the boundary), so it is what the flattening with `u = 0` has
    /// there. All zero on a closed mesh.
    pub boundary_curvature: Vec<f64>,
    /// The adjoint state `φ`: `(L φ)_i = −wE_i m_i u_i`, with `φ_i = ±λ wR_i` at cones strictly
    /// inside their bounds and `|φ_i| ≤ λ wR_i` at unbounded vertices without a cone.
    pub adjoint: Vec<f64>,
    /// Closed meshes: the interval `(lo, hi)` of constants that may be added to `adjoint`
    /// without violating the optimality system; `lo ≤ 0 ≤ hi`. It is `(0, 0)` whenever some
    /// cone lies strictly inside its bounds (there `φ = ±λ wR` pins the constant) or a free
    /// site does, and always on meshes with boundary. If the interval is not a point, the
    /// returned `adjoint` is the admissible one whose mass-weighted mean `Σ m φ` is closest to
    /// zero.
    pub adjoint_shift: (f64, f64),
    /// The objective `½ Σ wE m u² + λ Σ wR |µ|` of the optimal measure.
    pub objective: f64,
    /// The distortion `½ Σ wE m u²` of the optimal measure.
    pub distortion_optimal: f64,
    /// The cone cost `λ Σ wR |µ|` of the optimal measure (prescribed cones cost nothing).
    pub cone_cost: f64,
    /// The cones: vertex and angle, by ascending vertex.
    ///
    /// Range of the angles: every angle lies within the bounds given for its vertex
    /// ([`ConeProblem::lower_bounds`], [`ConeProblem::upper_bounds`]), exactly, also after
    /// merging. Without bounds nothing limits an angle. In particular the linearized problem
    /// does not know that a cone angle must be below `2π`: a free cone site on a closed convex
    /// mesh takes up to `4π`, and a merged cone carries the sum of its cluster.
    /// [`ExtractionReport::n_beyond_two_pi`] counts the returned angles outside `(−2π, 2π)`.
    pub cones: Vec<Cone>,
    /// The log scale factor `u` of the returned `cones`. Equal to `log_scale_optimal` unless
    /// extraction changed the cones.
    pub log_scale: Vec<f64>,
    /// The distortion `½ Σ wE m u²` of the returned `cones`.
    pub distortion: f64,
    /// The total cone angle `Σ |angle|` of the returned `cones`.
    pub total_cone_angle: f64,
    /// What cone extraction did.
    pub extraction: ExtractionReport,
    /// What the solver did, and the certificate.
    pub report: SolveReport,
}

/// Largest stiffness `k_i` of a proximal step. Stiffer problems converge in fewer proximal
/// steps but their Newton systems are worse conditioned, and the regularized cone angle
/// `k (φ − t)` is evaluated from `φ` with an error of about `k · 10⁻¹⁶`.
const MAX_PROXIMAL_STIFFNESS: f64 = 1e8;
/// The stiffness cap is never lowered below this.
const MIN_PROXIMAL_STIFFNESS: f64 = 1e2;
/// A proximal step that moves no cone angle by more than this (relative to the largest angle,
/// at least 1) is a fixed point: its cone angles solve the exact problem.
const FIXED_POINT_TOLERANCE: f64 = 1e-10;

/// The stiffness `k_i` of a regularized problem, i.e. the slope of its cone angle
/// `µ_i = clamp(k_i · P_{t_i}(φ_i))` in the adjoint: `m_i/γ` for the mass metric, `1/γ` for the
/// unweighted one, at most `cap`.
///
/// [formula: regularized-cone-update]
fn stiffness(d: &Discrete, metric: RegularizationMetric, gamma: f64, cap: f64) -> Vec<f64> {
    d.mass
        .iter()
        .map(|&m| {
            let weight = match metric {
                RegularizationMetric::Mass => m,
                RegularizationMetric::Unweighted => 1.0,
            };
            (weight / gamma).min(cap)
        })
        .collect()
}

fn stage(kind: StageKind, safeguard: bool, outcome: &NewtonOutcome, start: Instant) -> StageReport {
    StageReport {
        kind,
        safeguard,
        newton_iterations: outcome.iterations,
        damped_steps: outcome.damped_steps,
        constant_shifts: outcome.constant_shifts,
        stop: outcome.stop,
        seconds: start.elapsed().as_secs_f64(),
    }
}

fn max_difference(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Runs one regularized problem, records its stage, and returns how it ended
/// ([`StageStop::LinearSolverFailure`] if a Newton system could not be solved accurately).
#[allow(clippy::too_many_arguments)]
fn run_regularized(
    d: &Discrete,
    solver: &mut dyn NewtonLinearSolver,
    (kind, safeguard): (StageKind, bool),
    (k, center): (&[f64], &[f64]),
    x: &mut Iterate,
    max_iterations: usize,
    par: Parallelism<'_>,
    (stages, message): (&mut Vec<StageReport>, &mut Option<String>),
) -> StageStop {
    let start = Instant::now();
    match solve_regularized(d, solver, k, center, x, max_iterations, par) {
        Ok(outcome) => {
            match outcome.stop {
                StageStop::Stalled => {
                    *message = Some(
                        "a regularized problem could not be solved: its Newton systems are too \
                         ill-conditioned"
                            .to_string(),
                    )
                }
                StageStop::IterationLimit => {
                    *message = Some(
                        "a regularized problem reached its iteration limit before it was solved"
                            .to_string(),
                    )
                }
                _ => {}
            }
            stages.push(stage(kind, safeguard, &outcome, start));
            outcome.stop
        }
        Err(failure) => {
            *message = Some(failure.message);
            stages.push(failed_stage(kind, safeguard, start));
            StageStop::LinearSolverFailure
        }
    }
}

fn failed_stage(kind: StageKind, safeguard: bool, start: Instant) -> StageReport {
    StageReport {
        kind,
        safeguard,
        newton_iterations: 0,
        damped_steps: 0,
        constant_shifts: 0,
        stop: StageStop::LinearSolverFailure,
        seconds: start.elapsed().as_secs_f64(),
    }
}

/// What the two-stage iteration produced.
struct Driven {
    x: Iterate,
    stages: Vec<StageReport>,
    message: Option<String>,
    linear_failure: bool,
    safeguard_rounds: usize,
    linear_report: LinearSolverReport,
    regularized_seconds: f64,
    exact_seconds: f64,
}

impl Driven {
    /// The result of the factor-once method, in the shape of a driver result.
    fn fast(x: Iterate, message: Option<String>) -> Driven {
        Driven {
            x,
            stages: Vec::new(),
            message,
            linear_failure: false,
            safeguard_rounds: 0,
            linear_report: LinearSolverReport::default(),
            regularized_seconds: 0.0,
            exact_seconds: 0.0,
        }
    }
}

/// The optimum as returned, after the rounding-level clean-up, with its certificate.
struct Optimum {
    /// The cleaned iterate.
    x: Iterate,
    measure: Vec<f64>,
    log_scale_optimal: Vec<f64>,
    adjoint: Vec<f64>,
    adjoint_shift: (f64, f64),
    distortion_optimal: f64,
    cone_cost: f64,
    certificate: Certificate,
    bound_correction: f64,
    zeroed_angles: usize,
    largest_zeroed_angle: f64,
    seconds: f64,
}

/// Cleans up the iterate (bounds, tiny angles, the adjoint constant of a closed mesh) and
/// computes the certificate of what is returned.
fn optimum(
    instance: &Instance<'_>,
    d: &Discrete,
    options: &SolveOptions<'_>,
    mut x: Iterate,
) -> Optimum {
    // Rounding-level clean-up of the cone angles; the certificate is computed afterwards.
    // The bounds are constraints: a cone angle that a step ending by tolerance left beyond
    // its bound by rounding is moved onto it. And a cone angle below the Newton tolerance is
    // no cone: an active vertex can carry such an angle when it absorbs the rounding error of
    // the data, for instance of Σ Ω on a closed mesh.
    let mut bound_correction = 0.0f64;
    let mut zeroed_angles = 0;
    let mut largest_zeroed_angle = 0.0f64;
    for (mu, site) in x.mu.iter_mut().zip(&d.sites) {
        let clamped = mu.clamp(site.lo, site.hi);
        bound_correction = bound_correction.max((clamped - *mu).abs());
        *mu = clamped;
        if !site.fixed && *mu != 0.0 && mu.abs() <= options.newton_tolerance {
            zeroed_angles += 1;
            largest_zeroed_angle = largest_zeroed_angle.max(mu.abs());
            *mu = 0.0;
        }
    }

    // On a closed mesh the adjoint is determined up to the constants that keep it in the
    // subdifferentials. Return the admissible adjoint whose mass-weighted mean is closest to
    // zero, and the interval of the others.
    let mut adjoint_shift = (0.0, 0.0);
    if d.closed {
        let (lo, hi) = adjoint_constant_interval(&d.sites, &x.phi, &x.mu);
        if lo <= hi {
            let mean = x.phi.iter().zip(&d.mass).map(|(p, m)| p * m).sum::<f64>()
                / d.mass.iter().sum::<f64>();
            let c = (-mean).clamp(lo, hi);
            for phi in &mut x.phi {
                *phi += c;
            }
            adjoint_shift = (lo - c, hi - c);
        }
    }

    // The optimum, per vertex.
    let measure = d.to_vertices(&x.mu);
    let log_scale_optimal = d.to_vertices(&x.u);
    let adjoint = d.to_vertices(&x.phi);
    let distortion_optimal = d.distortion(&x.u);
    let cone_cost = d.cone_cost(&x.mu);

    let start = Instant::now();
    let certificate = certify_instance(instance, &log_scale_optimal, &adjoint, &measure);
    Optimum {
        x,
        measure,
        log_scale_optimal,
        adjoint,
        adjoint_shift,
        distortion_optimal,
        cone_cost,
        certificate,
        bound_correction,
        zeroed_angles,
        largest_zeroed_angle,
        seconds: start.elapsed().as_secs_f64(),
    }
}

/// Stages I and II with the safeguards, on an assembled problem.
///
/// [formula: proximal-fixed-point]
fn drive(d: &Discrete, options: &SolveOptions<'_>) -> Result<Driven, ConeError> {
    let n = d.n;
    let par = options.parallelism;
    let mut regularized_seconds = 0.0;
    let mut exact_seconds = 0.0;
    let mut x = Iterate::zeros(n);
    let mut stages: Vec<StageReport> = Vec::new();
    let mut message: Option<String> = None;
    let mut linear_failure = false;
    let mut safeguard_rounds = 0;
    let mut linear_report = LinearSolverReport::default();

    if n > 0 {
        let mut solver = SaddleSolver::new(d, options.linear_solver)?;
        // The regularized problem that ran last: its parameter, stiffness and centre, and how
        // its Newton iteration ended.
        let mut gamma = 1.0;
        let mut k: Vec<f64> = Vec::new();
        let mut center = vec![0.0; n];
        let mut last_stop = StageStop::ActiveSetsUnchanged;

        // Stage I: the regularized problems.
        let stage_start = Instant::now();
        for &next in &options.regularization {
            gamma = next;
            k = stiffness(d, options.regularization_metric, gamma, f64::INFINITY);
            last_stop = run_regularized(
                d,
                &mut solver,
                (StageKind::Regularized { gamma }, false),
                (&k, &center),
                &mut x,
                options.max_regularized_iterations,
                par,
                (&mut stages, &mut message),
            );
            if last_stop == StageStop::LinearSolverFailure {
                // Too stiff for an accurate solve: the safeguard below lowers the stiffness.
                break;
            }
        }
        regularized_seconds = stage_start.elapsed().as_secs_f64();

        // Stage II: the exact problem, with further regularized problems as a safeguard.
        let stage_start = Instant::now();
        let mut stiffness_cap = MAX_PROXIMAL_STIFFNESS;
        // Whether the last proximal step did not move the cone angles: then they solve the
        // exact problem.
        let mut fixed_point = false;
        loop {
            // The exact iteration is only worth starting from a solved regularized problem.
            if last_stop.is_success() {
                let start = Instant::now();
                match solve_exact(
                    d,
                    &mut solver,
                    &mut x,
                    options.max_exact_iterations,
                    options.newton_tolerance,
                    par,
                ) {
                    Ok(outcome) => {
                        stages.push(stage(
                            StageKind::Exact,
                            safeguard_rounds > 0,
                            &outcome,
                            start,
                        ));
                        if outcome.stop.is_success() {
                            message = None;
                            break;
                        }
                        let what = match outcome.stop {
                            StageStop::Cycle => "returned to active sets it had already tried",
                            _ => "reached its iteration limit",
                        };
                        message = Some(format!("the exact Newton iteration {what}"));
                    }
                    Err(failure) => {
                        stages.push(failed_stage(StageKind::Exact, safeguard_rounds > 0, start));
                        message = Some(format!(
                            "a Newton system of the exact problem could not be solved (the \
                             minimizer may not be unique): {}",
                            failure.message
                        ));
                    }
                }
                if fixed_point {
                    // The proximal iteration has converged: the iterate is a minimizer even
                    // though the exact iteration cannot confirm it. The certificate decides.
                    break;
                }
            }
            let ill_conditioned = matches!(
                last_stop,
                StageStop::Stalled | StageStop::LinearSolverFailure
            );
            if safeguard_rounds == options.max_safeguard_rounds
                || (ill_conditioned && stiffness_cap <= MIN_PROXIMAL_STIFFNESS)
            {
                linear_failure = last_stop == StageStop::LinearSolverFailure;
                break;
            }
            safeguard_rounds += 1;
            match last_stop {
                // Out of iterations: continue the same regularized problem.
                StageStop::IterationLimit if !k.is_empty() => {}
                // Too ill-conditioned: a proximal step with a hundred times less stiffness.
                StageStop::Stalled | StageStop::LinearSolverFailure => {
                    stiffness_cap = (stiffness_cap * 1e-2).max(MIN_PROXIMAL_STIFFNESS);
                    k = stiffness(d, options.regularization_metric, gamma, stiffness_cap);
                    center.copy_from_slice(&x.mu);
                }
                // Solved, but the exact iteration failed: a proximal step, ten times stiffer.
                _ => {
                    gamma *= 0.1;
                    k = stiffness(d, options.regularization_metric, gamma, stiffness_cap);
                    center.copy_from_slice(&x.mu);
                }
            }
            let kind = if center.iter().all(|&c| c == 0.0) {
                StageKind::Regularized { gamma }
            } else {
                StageKind::Proximal { gamma }
            };
            last_stop = run_regularized(
                d,
                &mut solver,
                (kind, true),
                (&k, &center),
                &mut x,
                options.max_regularized_iterations.max(1),
                par,
                (&mut stages, &mut message),
            );
            let moved = max_difference(&x.mu, &center);
            let size = x.mu.iter().fold(1.0f64, |m, v| m.max(v.abs()));
            fixed_point = last_stop.is_success() && moved <= FIXED_POINT_TOLERANCE * size;
        }
        exact_seconds = stage_start.elapsed().as_secs_f64();
        linear_report = solver.report();
    }
    Ok(Driven {
        x,
        stages,
        message,
        linear_failure,
        safeguard_rounds,
        linear_report,
        regularized_seconds,
        exact_seconds,
    })
}

/// Whether the factor-once method supports a problem: `Ok(None)` if it does, `Ok(Some(why))` if
/// [`SolveMethod::Auto`] would use the reference method for it (and
/// [`SolveMethod::FactorOnce`] would refuse it).
///
/// On a mesh with several components each component is a problem of its own: `Ok(None)` if the
/// method supports every component that has unknowns, otherwise the reason of the first
/// component (in the order of `TriMesh::face_components`) that it does not support; vertex
/// numbers in the reason are those of the whole mesh. A component without unknowns (a single
/// triangle with `u = 0` on its boundary) has nothing to solve and does not count;
/// [`Unsupported::NoUnknowns`] only if no component has unknowns.
///
/// # Errors
/// As [`solve_with_boundary`] for an invalid problem.
pub fn factor_once_support(
    problem: &ConeProblem<'_>,
    boundary: BoundaryTreatment,
) -> Result<Option<Unsupported>, ConeError> {
    let whole = Instance::new(problem, boundary)?;
    if problem.mesh.n_components() == 1 {
        let d = Discrete::from_instance(&whole)?;
        return Ok(factor_once::check_support(&d, problem.lambda).err());
    }
    let pieces = split_pieces(problem.mesh, problem.geometry)?;
    let mut any_unknowns = false;
    for (c, piece) in pieces.iter().enumerate() {
        let d = Discrete::from_instance(&whole.part(piece, c as u32))?;
        if d.n == 0 {
            continue;
        }
        any_unknowns = true;
        if let Err(why) = factor_once::check_support(&d, problem.lambda) {
            return Ok(Some(why));
        }
    }
    Ok((!any_unknowns).then_some(Unsupported::NoUnknowns))
}

/// Computes optimal cone singularities: solves the canonical problem and extracts the cones.
///
/// The method is the paper's (Fig. 12): Stage I solves a sequence of regularized problems by a
/// semismooth Newton method, Stage II solves the exact optimality system by the same method
/// started from the last regularized solution. If Stage II does not terminate, proximal steps
/// are inserted ([`SolveOptions::max_safeguard_rounds`]). The result always carries its
/// optimality [`Certificate`], and [`SolveReport::status`] is [`Status::Converged`] only if the
/// certificate holds.
///
/// The safeguard rests on one fact: a minimizer `(u, µ)` of
/// `objective + Σ (µ_i − µ̄_i)² / (2 k_i)` with `µ = µ̄` is a minimizer of the objective itself
/// (the penalty and its derivative vanish there), and the iteration `µ̄ ← µ` decreases the
/// objective. It is used only to reach a point; whether that point is optimal is decided by the
/// certificate.
///
/// **Disconnected meshes.** The mesh may have several connected components. The result is the
/// exact solution of the canonical problem on the whole mesh: the masses are divided by the total
/// area of the whole mesh, `λ` is the given one, and every closed component has its own additive
/// constant of `u` and its own Gauss–Bonnet constraint. The problem separates: each component is
/// solved on its own with these masses (concurrently on the pool of
/// [`SolveOptions::parallelism`], each component on one thread), and the results are assembled.
/// [`SolveReport::components`] reports every component (status, iterations, method, and the
/// certificate of its part of the problem); [`SolveReport::certificate`] is the certificate of the
/// whole problem, computed from the assembled triple, and it alone decides
/// [`SolveReport::status`]. A component solved alone with the masses normalized by its own area
/// needs `λ · A / A_c` for the same result ([`ComponentReport::equivalent_lambda`]).
///
/// # Errors
/// [`ConeError`] if the problem or the options are invalid, or the problem is infeasible. A
/// solve that does not converge is **not** an error; check `report.status`.
///
/// [formula: proximal-fixed-point]
/// [formula: component-separation]
pub fn solve(problem: &ConeProblem<'_>, options: &SolveOptions<'_>) -> Result<Solution, ConeError> {
    solve_with_boundary(problem, BoundaryTreatment::Fixed, options)
}

/// What the methods produced for one problem on a connected mesh (the whole problem, or one
/// component).
struct Ran {
    driven: Driven,
    opt: Optimum,
    method: MethodReport,
    factor_once_seconds: f64,
}

/// Solves one problem on a connected mesh: the factor-once method first unless the reference
/// method is forced, the reference method when that gives up.
fn run(
    instance: &Instance<'_>,
    d: &Discrete,
    lambda: f64,
    options: &SolveOptions<'_>,
) -> Result<Ran, ConeError> {
    let par = options.parallelism;
    let mut method = MethodReport {
        requested: options.method,
        used: MethodUsed::Reference,
        fallback: None,
        factor_once: None,
    };
    let mut factor_once_seconds = 0.0;
    let mut accepted: Option<(Driven, Optimum)> = None;
    if d.n == 0 && instance.component.is_some() {
        // A component without unknowns (a single triangle with u = 0 on its boundary): there is
        // nothing to solve, with either method.
        if options.method != SolveMethod::Reference {
            method.fallback = Some(FallbackReason::Unsupported(Unsupported::NoUnknowns));
        }
    } else if options.method != SolveMethod::Reference {
        let forced = options.method == SolveMethod::FactorOnce;
        let start = Instant::now();
        let attempt = factor_once::attempt(d, lambda, &options.factor_once, par, forced);
        factor_once_seconds = start.elapsed().as_secs_f64();
        method.factor_once = Some(attempt.report.clone());
        let mut message = None;
        let mut reason = attempt.gave_up;
        if let Some(x) = attempt.x {
            let opt = optimum(instance, d, options, x);
            if reason.is_none() && !opt.certificate.holds(options.certificate_tolerance) {
                reason = Some(FallbackReason::CertificateFailed {
                    certificate: Box::new(opt.certificate),
                });
            }
            if forced {
                if let Some(r) = &reason {
                    message = Some(format!("the factor-once method did not finish: {r}"));
                }
                accepted = Some((Driven::fast(opt.x.clone(), message), opt));
                method.used = MethodUsed::FactorOnce;
            } else if reason.is_none() {
                accepted = Some((Driven::fast(opt.x.clone(), None), opt));
                method.used = MethodUsed::FactorOnce;
            }
        } else if forced {
            let why = reason.map_or_else(|| "no result".to_string(), |r| r.to_string());
            return Err(ConeError::InvalidOptions {
                reason: format!("SolveMethod::FactorOnce cannot solve this problem: {why}"),
            });
        }
        method.fallback = reason;
    }
    let (driven, opt) = match accepted {
        Some(pair) => pair,
        None => {
            let driven = drive(d, options)?;
            let opt = optimum(instance, d, options, driven.x.clone());
            (driven, opt)
        }
    };
    Ok(Ran {
        driven,
        opt,
        method,
        factor_once_seconds,
    })
}

/// The status of a result: the certificate decides.
fn status_of(
    certificate: &Certificate,
    linear_failure: bool,
    options: &SolveOptions<'_>,
) -> Status {
    if certificate.holds(options.certificate_tolerance) {
        Status::Converged
    } else if linear_failure {
        Status::LinearSolverFailed
    } else {
        Status::NotConverged
    }
}

/// `(regularized, safeguard, exact)` Newton steps of a list of stages.
fn iteration_counts(stages: &[StageReport]) -> (usize, usize, usize) {
    let regularized = stages
        .iter()
        .filter(|s| !s.safeguard && matches!(s.kind, StageKind::Regularized { .. }))
        .map(|s| s.newton_iterations)
        .sum();
    let safeguard = stages
        .iter()
        .filter(|s| s.safeguard && !matches!(s.kind, StageKind::Exact))
        .map(|s| s.newton_iterations)
        .sum();
    let exact = stages
        .iter()
        .filter(|s| matches!(s.kind, StageKind::Exact))
        .map(|s| s.newton_iterations)
        .sum();
    (regularized, safeguard, exact)
}

fn extract_options(options: &SolveOptions<'_>) -> ExtractOptions {
    ExtractOptions {
        threshold: options.cone_threshold,
        merge: options.merge_cone_clusters,
        max_cluster_size: options.max_merge_cluster_size,
        max_cluster_diameter: options.max_merge_cluster_diameter,
    }
}

/// The measure with the boundary entries zeroed: only interior vertices carry cone angles (with
/// an optimized boundary the boundary entries of `measure` are the curvatures `ν`).
fn interior_measure(mesh: &TriMesh, boundary: BoundaryTreatment, measure: &[f64]) -> Vec<f64> {
    let mut cone_measure = measure.to_vec();
    if boundary != BoundaryTreatment::Fixed {
        for (v, m) in cone_measure.iter_mut().enumerate() {
            if mesh.is_boundary_vertex(v) {
                *m = 0.0;
            }
        }
    }
    cone_measure
}

/// The scale factor (per vertex of `d`'s mesh) and distortion of given cones, with the boundary
/// values optimized again for them if the boundary is unknown; and a message if that did not
/// reach the tolerance.
fn scale_of_cones(
    d: &Discrete,
    boundary: BoundaryTreatment,
    cones: &[Cone],
    options: &SolveOptions<'_>,
) -> Result<(Vec<f64>, f64, Option<String>), ConeError> {
    let par = options.parallelism;
    if boundary == BoundaryTreatment::Fixed || d.n < d.n_vertices {
        let mu = cones_to_measure(d, cones)?;
        let laplace = LaplaceSolver::new(d, par)?;
        let (u, _) = log_scale_of_measure(d, &laplace, &mu, par)?;
        Ok((d.to_vertices(&u), d.distortion(&u), None))
    } else {
        // Unknown boundary: the interior angles are those of the cones; the boundary values are
        // optimized again for them, with the same treatment.
        let fixed = d.with_fixed_interior(cones)?;
        let again = drive(&fixed, options)?;
        let message = again.message.map(|why| {
            format!(
                "the boundary values for the extracted cones were not obtained to tolerance: {why}"
            )
        });
        Ok((d.to_vertices(&again.x.u), d.distortion(&again.x.u), message))
    }
}

/// `ν` at the boundary vertices: `(L u)_i + k_i` for `Fixed`, the optimized variable otherwise.
pub(crate) fn boundary_curvature_of(
    problem: &ConeProblem<'_>,
    boundary: BoundaryTreatment,
    log_scale_optimal: &[f64],
    measure: &[f64],
) -> Result<Vec<f64>, ConeError> {
    let nv = problem.mesh.n_vertices();
    let mut boundary_curvature = vec![0.0; nv];
    if problem.mesh.n_boundary_edges() > 0 {
        if boundary == BoundaryTreatment::Fixed {
            let k = boundary_geodesic_curvatures(problem.mesh, problem.geometry)?;
            let l_u = apply_laplacian(problem.mesh, problem.geometry, log_scale_optimal);
            for v in 0..nv {
                if problem.mesh.is_boundary_vertex(v) {
                    boundary_curvature[v] = l_u[v] + k[v];
                }
            }
        } else {
            for v in 0..nv {
                if problem.mesh.is_boundary_vertex(v) {
                    boundary_curvature[v] = measure[v];
                }
            }
        }
    }
    Ok(boundary_curvature)
}

/// [`solve`] with a choice of what happens at the boundary: [`BoundaryTreatment::Fixed`] is
/// `solve` itself; the other treatments make the boundary values and curvatures unknowns (paper
/// Sec. 8.3; see [`BoundaryTreatment`]).
///
/// With a treatment other than `Fixed` on a mesh with boundary:
///
/// * `u` is returned at the boundary vertices too, and `measure[v]` at a boundary vertex `v`
///   is the curvature `ν_v` of the flattened boundary there (the cone angles `µ` stay at the
///   interior vertices; `cones` lists interior vertices only);
/// * the objective, `distortion_optimal` and `cone_cost` include the boundary vertices (with
///   their masses and distortion weights; the cost of `ν` for `Penalized`);
/// * `u` is determined up to a constant, which is optimized (`Σ wE m u = 0`), as on a closed
///   mesh; `adjoint` is determined up to a constant unless some boundary site pins it
///   (`adjoint_shift`);
/// * the certificate is that of the problem with the boundary as unknowns, and its
///   `gauss_bonnet_residual` is `|Σ_I µ + Σ_B ν − 2πχ|`;
/// * if extraction changes the cones (a threshold, merging), the scale factor of the returned
///   cones, `log_scale`, is the one whose boundary values are re-optimized for the new
///   interior cones (the same treatment of the boundary, interior angles fixed).
///
/// [`Solution::boundary_curvature`] reports `ν` for every treatment. Disconnected meshes are
/// handled as in [`solve`]: per component (with an unknown boundary every component with
/// boundary has its own constant of `u` and its own Gauss–Bonnet sum).
///
/// # Errors
/// As [`solve`]. `ConeError::Infeasible` can occur with `Convex` (and prescribed or
/// forbidden vertices).
///
/// [formula: proximal-fixed-point]
/// [formula: boundary-optimality-system]
pub fn solve_with_boundary(
    problem: &ConeProblem<'_>,
    boundary: BoundaryTreatment,
    options: &SolveOptions<'_>,
) -> Result<Solution, ConeError> {
    let total_start = Instant::now();
    options.validate()?;
    if problem.mesh.n_components() > 1 {
        return solve_components(problem, boundary, options, total_start);
    }
    let instance = Instance::new(problem, boundary)?;
    let d = Discrete::from_instance(&instance)?;
    let mut timings = Timings {
        setup: total_start.elapsed().as_secs_f64(),
        ..Timings::default()
    };

    let Ran {
        driven,
        opt,
        method,
        factor_once_seconds,
    } = run(&instance, &d, problem.lambda, options)?;
    timings.factor_once = factor_once_seconds;
    let Driven {
        x: _,
        stages,
        mut message,
        linear_failure,
        safeguard_rounds,
        linear_report,
        regularized_seconds,
        exact_seconds,
    } = driven;
    timings.regularized = regularized_seconds;
    timings.exact = exact_seconds;
    let Optimum {
        x: _,
        measure,
        log_scale_optimal,
        adjoint,
        adjoint_shift,
        distortion_optimal,
        cone_cost,
        certificate,
        bound_correction,
        zeroed_angles,
        largest_zeroed_angle,
        seconds: certificate_seconds,
    } = opt;
    timings.certificate = certificate_seconds;
    let status = status_of(&certificate, linear_failure, options);
    if status == Status::Converged && safeguard_rounds == 0 {
        message = None;
    }

    // Cone extraction.
    let start = Instant::now();
    let cone_measure = interior_measure(problem.mesh, boundary, &measure);
    let (cones, mut extraction) = extract_cones(
        problem.mesh,
        problem.geometry,
        &instance.data,
        &problem.prescribed_cones,
        &cone_measure,
        &extract_options(options),
    );
    let total_cone_angle: f64 = cones.iter().map(|c| c.angle.abs()).sum();
    let (log_scale, distortion) = if !extraction.cones_differ_from_measure {
        (log_scale_optimal.clone(), distortion_optimal)
    } else {
        let (u, distortion, why) = scale_of_cones(&d, boundary, &cones, options)?;
        if why.is_some() {
            message = why;
        }
        (u, distortion)
    };
    extraction.distortion_change = distortion - distortion_optimal;
    extraction.total_angle_change = if extraction.cones_differ_from_measure {
        total_cone_angle - measure.iter().map(|m| m.abs()).sum::<f64>()
    } else {
        0.0
    };
    timings.extraction = start.elapsed().as_secs_f64();
    timings.total = total_start.elapsed().as_secs_f64();

    let (regularized_iterations, safeguard_iterations, exact_iterations) =
        iteration_counts(&stages);
    let boundary_curvature =
        boundary_curvature_of(problem, boundary, &log_scale_optimal, &measure)?;
    Ok(Solution {
        measure,
        log_scale_optimal,
        boundary_curvature,
        adjoint,
        adjoint_shift,
        objective: distortion_optimal + cone_cost,
        distortion_optimal,
        cone_cost,
        cones,
        log_scale,
        distortion,
        total_cone_angle,
        extraction,
        report: SolveReport {
            status,
            stages,
            regularized_iterations,
            exact_iterations,
            safeguard_rounds,
            safeguard_iterations,
            linear_solver: linear_report,
            timings,
            bound_correction,
            zeroed_angles,
            largest_zeroed_angle,
            certificate,
            method,
            message,
            components: Vec::new(),
        },
    })
}

/// [`solve_with_boundary`] on a mesh with several components: each component is solved as its
/// part of the whole problem (masses of the whole mesh, the given `λ`), then the results are
/// assembled and the whole problem is certified.
///
/// [formula: component-separation]
fn solve_components(
    problem: &ConeProblem<'_>,
    boundary: BoundaryTreatment,
    options: &SolveOptions<'_>,
    total_start: Instant,
) -> Result<Solution, ConeError> {
    let mesh = problem.mesh;
    let whole = Instance::new(problem, boundary)?;
    let pieces = split_pieces(mesh, problem.geometry)?;
    let parts: Vec<Instance<'_>> = pieces
        .iter()
        .enumerate()
        .map(|(c, piece)| whole.part(piece, c as u32))
        .collect();
    let discretes = parts
        .iter()
        .map(Discrete::from_instance)
        .collect::<Result<Vec<_>, _>>()?;
    let n_parts = pieces.len();
    let mut timings = Timings {
        setup: total_start.elapsed().as_secs_f64(),
        ..Timings::default()
    };

    // The components are independent problems: with a pool they run concurrently, each on one
    // thread (so the result of a component does not depend on the number of threads).
    let part_options = SolveOptions {
        parallelism: Parallelism::Sequential,
        ..options.clone()
    };
    let solve_one = |c: usize| solve_part(&parts[c], &discretes[c], problem.lambda, &part_options);
    let results: Vec<Result<PartSolve, ConeError>> = match options.parallelism {
        Parallelism::Pool(pool) => pool.map(n_parts, solve_one),
        Parallelism::Sequential => (0..n_parts).map(solve_one).collect(),
    };
    let mut solved = Vec::with_capacity(n_parts);
    for result in results {
        solved.push(result?);
    }

    // Assemble.
    let nv = mesh.n_vertices();
    let mut measure = vec![0.0; nv];
    let mut log_scale_optimal = vec![0.0; nv];
    let mut adjoint = vec![0.0; nv];
    let mut stages = Vec::new();
    let mut linear_report = LinearSolverReport::default();
    let mut distortion_optimal = 0.0;
    let mut cone_cost = 0.0;
    // A constant added to the whole adjoint must be admissible on every component; a component
    // that is not closed admits none.
    let mut adjoint_shift = (f64::NEG_INFINITY, f64::INFINITY);
    let mut bound_correction = 0.0f64;
    let mut zeroed_angles = 0;
    let mut largest_zeroed_angle = 0.0f64;
    let mut safeguard_rounds = 0;
    let mut linear_failure = false;
    let mut messages = Vec::new();
    let mut all_factor_once = true;
    let mut fallback = None;
    let mut components = Vec::with_capacity(n_parts);
    for (c, part) in solved.into_iter().enumerate() {
        let PartSolve {
            measure: part_measure,
            log_scale_optimal: part_u,
            adjoint: part_phi,
            linear_failure: part_failure,
            report: component,
        } = part;
        for (local, &v) in pieces[c].sub.vertex_to_parent.iter().enumerate() {
            let v = v as usize;
            measure[v] = part_measure[local];
            log_scale_optimal[v] = part_u[local];
            adjoint[v] = part_phi[local];
        }
        let r = &component.report;
        distortion_optimal += component.distortion_optimal;
        cone_cost += component.cone_cost;
        adjoint_shift = (
            adjoint_shift.0.max(component.adjoint_shift.0),
            adjoint_shift.1.min(component.adjoint_shift.1),
        );
        bound_correction = bound_correction.max(r.bound_correction);
        zeroed_angles += r.zeroed_angles;
        largest_zeroed_angle = largest_zeroed_angle.max(r.largest_zeroed_angle);
        safeguard_rounds += r.safeguard_rounds;
        linear_failure |= part_failure;
        linear_report.accumulate(&r.linear_solver);
        timings.regularized += r.timings.regularized;
        timings.exact += r.timings.exact;
        timings.factor_once += r.timings.factor_once;
        timings.certificate += r.timings.certificate;
        if component.n_unknowns > 0 {
            all_factor_once &= r.method.used == MethodUsed::FactorOnce;
            if fallback.is_none() {
                fallback.clone_from(&r.method.fallback);
            }
        }
        if let Some(m) = &r.message {
            messages.push(format!("component {c}: {m}"));
        }
        stages.extend(r.stages.iter().copied());
        components.push(component);
    }

    // The certificate of the whole problem decides the status.
    let start = Instant::now();
    let certificate = certify_instance(&whole, &log_scale_optimal, &adjoint, &measure);
    timings.certificate += start.elapsed().as_secs_f64();
    let status = status_of(&certificate, linear_failure, options);
    let mut message = (!messages.is_empty()).then(|| messages.join("; "));
    let method = MethodReport {
        requested: options.method,
        used: if all_factor_once && components.iter().any(|c| c.n_unknowns > 0) {
            MethodUsed::FactorOnce
        } else {
            MethodUsed::Reference
        },
        fallback,
        factor_once: None,
    };

    // Cone extraction, on the whole mesh (a cluster never spans two components); the scale factor
    // of changed cones per component.
    let start = Instant::now();
    let cone_measure = interior_measure(mesh, boundary, &measure);
    let (cones, mut extraction) = extract_cones(
        mesh,
        problem.geometry,
        &whole.data,
        &problem.prescribed_cones,
        &cone_measure,
        &extract_options(options),
    );
    let total_cone_angle: f64 = cones.iter().map(|c| c.angle.abs()).sum();
    let (log_scale, distortion) = if !extraction.cones_differ_from_measure {
        (log_scale_optimal.clone(), distortion_optimal)
    } else {
        let component_of = mesh.vertex_components();
        let mut local = vec![0u32; nv];
        for piece in &pieces {
            for (l, &v) in piece.sub.vertex_to_parent.iter().enumerate() {
                local[v as usize] = l as u32;
            }
        }
        let mut per_part: Vec<Vec<Cone>> = vec![Vec::new(); n_parts];
        let mut extracted = vec![0.0; nv];
        for cone in &cones {
            let v = cone.vertex as usize;
            extracted[v] = cone.angle;
            per_part[component_of[v] as usize].push(Cone {
                vertex: local[v],
                angle: cone.angle,
            });
        }
        let mut u_all = log_scale_optimal.clone();
        let mut distortion = 0.0;
        for c in 0..n_parts {
            let vertices = &pieces[c].sub.vertex_to_parent;
            let unchanged = vertices
                .iter()
                .all(|&v| extracted[v as usize] == cone_measure[v as usize]);
            if unchanged {
                distortion += components[c].distortion_optimal;
                continue;
            }
            let (u, part_distortion, why) =
                scale_of_cones(&discretes[c], boundary, &per_part[c], &part_options)?;
            for (l, &v) in vertices.iter().enumerate() {
                u_all[v as usize] = u[l];
            }
            distortion += part_distortion;
            if let Some(why) = why {
                let why = format!("component {c}: {why}");
                message = Some(match message {
                    Some(m) => format!("{m}; {why}"),
                    None => why,
                });
            }
        }
        (u_all, distortion)
    };
    extraction.distortion_change = distortion - distortion_optimal;
    extraction.total_angle_change = if extraction.cones_differ_from_measure {
        total_cone_angle - measure.iter().map(|m| m.abs()).sum::<f64>()
    } else {
        0.0
    };
    timings.extraction = start.elapsed().as_secs_f64();
    timings.total = total_start.elapsed().as_secs_f64();

    let (regularized_iterations, safeguard_iterations, exact_iterations) =
        iteration_counts(&stages);
    let boundary_curvature =
        boundary_curvature_of(problem, boundary, &log_scale_optimal, &measure)?;
    Ok(Solution {
        measure,
        log_scale_optimal,
        boundary_curvature,
        adjoint,
        adjoint_shift,
        objective: distortion_optimal + cone_cost,
        distortion_optimal,
        cone_cost,
        cones,
        log_scale,
        distortion,
        total_cone_angle,
        extraction,
        report: SolveReport {
            status,
            stages,
            regularized_iterations,
            exact_iterations,
            safeguard_rounds,
            safeguard_iterations,
            linear_solver: linear_report,
            timings,
            bound_correction,
            zeroed_angles,
            largest_zeroed_angle,
            certificate,
            method,
            message,
            components,
        },
    })
}

/// The `λ` with which a component, solved alone with its masses divided by its own area, has the
/// minimizers of its part of the whole problem: `λ · A / A_c = λ / area_fraction` (masses
/// multiplied by `c = A / A_c` need `λ · c`).
///
/// [formula: component-lambda-rescaling]
/// [formula: lambda-area-normalization]
pub(crate) fn equivalent_lambda(lambda: f64, area_fraction: f64) -> f64 {
    lambda / area_fraction
}

/// The solve of one component of a disconnected mesh: its optimum per vertex of the component,
/// and its report.
pub(crate) struct PartSolve {
    pub measure: Vec<f64>,
    pub log_scale_optimal: Vec<f64>,
    pub adjoint: Vec<f64>,
    /// Whether a linear solver failure ended the reference method.
    pub linear_failure: bool,
    pub report: ComponentReport,
}

/// Solves one component (an [`Instance::part`]) of a disconnected mesh, as its part of the whole
/// problem.
///
/// [formula: component-separation]
pub(crate) fn solve_part(
    instance: &Instance<'_>,
    d: &Discrete,
    lambda: f64,
    options: &SolveOptions<'_>,
) -> Result<PartSolve, ConeError> {
    let start = Instant::now();
    let Ran {
        driven,
        opt,
        method,
        factor_once_seconds,
    } = run(instance, d, lambda, options)?;
    let status = status_of(&opt.certificate, driven.linear_failure, options);
    let mut message = driven.message;
    if status == Status::Converged && driven.safeguard_rounds == 0 {
        message = None;
    }
    let (regularized_iterations, safeguard_iterations, exact_iterations) =
        iteration_counts(&driven.stages);
    let area_fraction: f64 = instance.mass.iter().sum();
    let shift = if d.closed {
        opt.adjoint_shift
    } else {
        (0.0, 0.0)
    };
    let report = ComponentReport {
        component: instance.component.unwrap_or(0),
        n_vertices: instance.mesh.n_vertices(),
        n_faces: instance.mesh.n_faces(),
        closed: d.closed,
        n_unknowns: d.n,
        area_fraction,
        equivalent_lambda: equivalent_lambda(lambda, area_fraction),
        objective: opt.distortion_optimal + opt.cone_cost,
        distortion_optimal: opt.distortion_optimal,
        cone_cost: opt.cone_cost,
        adjoint_shift: shift,
        report: SolveReport {
            status,
            stages: driven.stages,
            regularized_iterations,
            exact_iterations,
            safeguard_rounds: driven.safeguard_rounds,
            safeguard_iterations,
            linear_solver: driven.linear_report,
            timings: Timings {
                regularized: driven.regularized_seconds,
                exact: driven.exact_seconds,
                factor_once: factor_once_seconds,
                certificate: opt.seconds,
                total: start.elapsed().as_secs_f64(),
                ..Timings::default()
            },
            bound_correction: opt.bound_correction,
            zeroed_angles: opt.zeroed_angles,
            largest_zeroed_angle: opt.largest_zeroed_angle,
            certificate: opt.certificate,
            method,
            message,
            components: Vec::new(),
        },
    };
    Ok(PartSolve {
        measure: opt.measure,
        log_scale_optimal: opt.log_scale_optimal,
        adjoint: opt.adjoint,
        linear_failure: driven.linear_failure,
        report,
    })
}
