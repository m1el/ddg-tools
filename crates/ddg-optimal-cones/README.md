# ddg-optimal-cones

Where to put cone singularities, and with which angles, so that a triangle mesh can be flattened
conformally with the least area distortion for the total cone angle spent. An implementation of

> Yousuf Soliman, Dejan Slepčev, Keenan Crane.
> *Optimal Cone Singularities for Conformal Flattening.*
> ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201367>

("the paper"). The method, the optimality system and the two-stage semismooth Newton algorithm
are theirs. What this crate adds is listed under [Deviations from the paper](#deviations-from-the-paper):
a certificate on every solve, safeguards, closed meshes, angle bounds, prescribed cones, and
[a second method](#the-factor-once-method) that factorizes the Laplacian once and is 3 to 50 times
faster (the default whenever it applies; the paper's method is the fallback and can be forced).

Contents: [example](#example) · [the problem](#the-problem) · [the boundary](#the-boundary) · [interface](#interface) ·
[steering](#steering) · [what comes back](#what-comes-back) · [the certificate](#the-certificate) ·
[algorithm](#algorithm) · [the factor-once method](#the-factor-once-method) ·
[deviations from the paper](#deviations-from-the-paper) ·
[the published dataset](#the-published-dataset) · [choosing λ](#choosing-λ) ·
[closed meshes](#closed-meshes) · [disconnected meshes](#disconnected-meshes) ·
[re-optimizing given cones](#re-optimizing-given-cones) · [needle-like spikes](#needle-like-spikes) ·
[the feedback loop](#the-feedback-loop) ·
[performance baseline](#performance-baseline) · [performance of the factor-once method](#performance-of-the-factor-once-method) · [formula tags](#formula-tags) · [tests](#tests) ·
[not implemented](#not-implemented)

## Example

```rust
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{ConeProblem, SolveOptions, Status, solve};

let mesh = TriMesh::new(positions.len(), triangles)?;        // manifold, oriented; components allowed
let geometry = Geometry::from_positions(&mesh, &positions)?; // or Geometry::from_edge_lengths

let problem = ConeProblem::new(&mesh, &geometry, 0.5);       // λ = 0.5
let solution = solve(&problem, &SolveOptions::default())?;

assert_eq!(solution.report.status, Status::Converged);       // not converging is not an Err
for cone in &solution.cones {
    println!("vertex {} angle {}", cone.vertex, cone.angle); // 0-based vertex, radians
}
```

The crate documentation has this example and a steering example as tested doctests.

## The problem

Conventions are those of the repository's `CLAUDE.md`. `L` is the positive semidefinite cotan
Laplacian, `Ω_i` the angle defect of interior vertex `i`, `m_i` the lumped vertex mass **divided
by the total area** (`Σ m_i = 1`), and `u` the true log conformal scale factor (flattened edge
length `≈ e^{(u_i+u_j)/2} ℓ_ij`). `solve` finds the cone angles `µ`, and the `u` they determine,
that

```
minimize     ½ Σ_i wE_i m_i u_i²  +  λ Σ_i wR_i |µ_i|
subject to   (L u)_i = µ_i − Ω_i  (i interior),   u_i = 0  (i boundary),   lo_i ≤ µ_i ≤ hi_i
```

A cone angle is the curvature of the flattened metric, `µ_i = 2π − (angle sum around i after
flattening)`: positive for a cone tip. This is the convention of the dataset's `*-cones.txt`.
On a closed mesh `u` is defined up to a constant, which is optimized too, and `Σ µ = 2πχ`.

`λ` is dimensionless here. If masses are left as they are (not divided by the area), the same
problem has `λ_raw = λ · area` ([choosing λ](#choosing-λ)).

## The boundary

By default `u = 0` at the boundary vertices (paper Fig. 10 b). [`solve_with_boundary`] also
optimizes the boundary (paper Sec. 8.3, Fig. 10 c–e); the treatment is an argument
(`BoundaryTreatment`), not a field of `ConeProblem`, so that existing struct literals of
`ConeProblem` keep compiling (`solve` is `solve_with_boundary(.., Fixed, ..)`, bit for bit).

The model. Let `B` be the boundary vertices, `I` the interior ones, `k_i = π − Σ(corner angles at i)`
the geodesic curvature of the boundary of the mesh at `i ∈ B`, and `ν_i` the geodesic curvature
of the flattened boundary there. The boundary vertices become unknowns:

```
minimize     ½ Σ_{i∈V} wE_i m_i u_i²  +  λ Σ_{i∈I} wR_i |µ_i|  +  Σ_{i∈B} cost_i(ν_i)
subject to   (L u)_i = µ_i − Ω_i   (i ∈ I),      (L u)_i = ν_i − k_i   (i ∈ B)
```

so every boundary vertex has a *site* for `ν_i`, as every interior vertex has one for `µ_i`:

| `BoundaryTreatment` | paper | `cost_i(ν_i)` |
|---|---|---|
| `Fixed` (default) | Fig. 10 b | not a variable: `u_i = 0` |
| `Free` | Fig. 10 c | `0`, `ν_i` unconstrained (optimized boundary values) |
| `Penalized` | Fig. 10 d | `λ wR_i \|ν_i\|`, `wR_i` from `cone_costs` (`+∞`: `ν_i = 0`); a polygonal boundary |
| `Convex` | Fig. 10 e | `0` with `ν_i ≥ 0` |

The distortion runs over all vertices (their masses and `wE`). With no vertex pinned `u` is
defined up to a constant that is optimized too (`Σ wE m u = 0`), as on a closed mesh, and
`Σ_I µ + Σ_B ν = 2πχ` holds automatically. In the other treatments the arrays `lower_bounds`,
`upper_bounds` have no effect at boundary vertices, `prescribed_cones` there is still an error.

**Optimality conditions** (derived from the statement above; the same as the closed-mesh system with
`Ω̂ = (Ω on I, k on B)` and a site per vertex): `(L u)_i = s_i − Ω̂_i`, `(L φ)_i = −wE_i m_i u_i`
for all vertices, `φ_i ∈ ∂g_i(s_i)` with `s = (µ, ν)`. At a boundary vertex: `Free` gives `φ_i = 0`
(a natural boundary condition: the adjoint vanishes on the boundary), `Penalized` gives
`φ_i ∈ λ wR_i ∂|ν_i|`, `Convex` gives `φ_i ≤ 0` with `φ_i ν_i = 0`. This is
`optimality-system-sufficient` for these data, hence sufficient for a global minimum; the solver is the
closed-mesh solver of the other sections (nothing new is iterated), so the Newton systems, the
safeguards and the certificate are the same. The modelling statement `ν_i = (L u)_i + k_i` (the
first-order change of the boundary geodesic curvature) is not proved here; it is checked
against two routes that do not share the boundary code: the generic convex solver, and the
double of the mesh (below).

**What comes back.** `Solution::log_scale_optimal` holds `u` at the boundary vertices,
`Solution::measure` holds `ν` at the boundary vertices (and `µ` at the interior ones), and
`Solution::boundary_curvature` holds `ν` at the boundary vertices in every treatment (for `Fixed`
it is `(L u)_i + k_i` with `u_B = 0`, which is not a variable there). `objective`,
`distortion_optimal` and `cone_cost` include the boundary vertices, and the penalty of `ν` for
`Penalized`. `cones` lists interior vertices only. `certify_with_boundary` certifies a triple
for the treatment (`gauss_bonnet_residual` is `|Σ_I µ + Σ_B ν − 2πχ|`). If extraction changes the
cones (a threshold, merging), `log_scale` and `distortion` are those of the new interior cones
with the boundary optimized again for them (same treatment).

**Verified against.**

* The generic convex solver Clarabel, with the problem written out again from the statement
  above (`tests/convex_reference.rs`): 60 random meshes per row (disk, annulus, cylinder) with
  all steering features, with and without zero distortion weights, for `Free`, `Penalized`,
  `Convex`, and 200 random cases of a property test. Every case agrees (objective to `10⁻¹¹`,
  `µ` and `ν` to `10⁻⁵` and better, the limit is Clarabel's), and the verdicts on infeasibility
  (`Convex` with the other constraints) agree. One degenerate case does not: a `Free` boundary
  vertex with zero distortion weight has a free site and no weight, its `u` is determined by
  nothing, and the minimum need not be attained (Clarabel returns an objective of `3·10⁻¹⁴` with
  `u` of order 200). The solver reports `NotConverged` (objective `9·10⁻⁶` after 12 safeguard
  rounds, `5·10⁻⁷` after 300); the sweeps keep boundary vertices weighted for `Free`, and a
  test pins down that this case is never reported as optimal.
* The double of the mesh (`tests/boundary.rs`, `tests/dataset.rs`): the mesh glued to a mirror
  copy along the boundary is a closed mesh, on which the closed-mesh path (the seam vertices free
  sites, `λ/2`) solves the same problem. Random meshes: agreement to `10⁻¹³` for the three
  treatments (objective, `µ`, `u`, `ν` = half the seam angle). Dataset, the 16 meshes with boundary
  at `λ = 0.5` that have a double (`face`, `hand`, `mannequin` have a triangle with three boundary
  vertices, whose double is not a mesh): relative objective to `10⁻¹¹` (`camelbad`; `10⁻¹³` or
  better on the others), `u` to `10⁻¹¹`, `ν` to `2·10⁻¹⁰`, `µ` to `10⁻¹²` except
  `camelshotnoise` (`3·10⁻⁹`).
* The optimized boundary is never worse than the others (the fixed-boundary optimum is feasible
  for it at no cost), tested on random meshes (strictly better, by more than `10⁻⁶`, in 93 of 120 cases) and on the dataset.

**Distortion with `u = 0` on the boundary against optimized boundary values**
(`cargo run --release -p ddg-optimal-cones --example boundary -- table`; `λ = 0.5`, the optimum of
each problem; `D = ½ Σ wE m u²`, objective `= D +` cone cost; the optimized boundary uses fewer
cones, so `D` alone can be larger):

| mesh | cones `u = 0` | `D`, `u = 0` | objective, `u = 0` | cones free | `D`, free | objective, free |
|---|---:|---:|---:|---:|---:|---:|
| body | 4 | 1.249 | 6.027 | 4 | 1.695 | 5.635 |
| bump | 0 | 0.1785 | 0.1785 | 0 | 0.06192 | 0.06192 |
| camel | 2 | 1.361 | 3.1 | 0 | 0.7725 | 0.7725 |
| camel_noise1 | 5 | 1.399 | 3.251 | 0 | 0.76 | 0.76 |
| camel_noise2 | 8 | 1.607 | 3.837 | 0 | 0.6146 | 0.6146 |
| camelbad | 2 | 1.36 | 3.098 | 0 | 0.768 | 0.768 |
| camelbadnonuniform | 1 | 1.348 | 3.082 | 0 | 0.7731 | 0.7731 |
| camelshotnoise | 452 | 3.493 | 122.6 | 63 | 14.69 | 102.2 |
| cow (5 files) | 2 | 1.89 | 3.525 | 0 | 1.321 | 1.321 |
| face | 0 | 0.01957 | 0.01957 | 0 | 0.00991 | 0.00991 |
| hand | 6 | 1.576 | 5.085 | 4 | 2.811 | 3.63 |
| mannequin | 4 | 0.8287 | 2.442 | 0 | 0.606 | 0.606 |
| rhino | 2 | 1.323 | 2.721 | 0 | 1.063 | 1.063 |
| tree | 4 | 2.017 | 2.29 | 0 | 0.5162 | 0.5162 |

The optimized boundary lowers the objective on every mesh with boundary (it can only), by a factor
of 1.07 (`body`) to 4.4 (`tree`). It lowers `D` too except on `airplane`, `body`, `hand` and
`camelshotnoise`, where it trades distortion for fewer cones. At `λ = 0.5` it places no cone at all
on most meshes: the curvature is absorbed by the boundary values. The example also prints
`Penalized` and `Convex` at the same `λ`. `Penalized` puts hundreds of cones on `bump` (833),
`mannequin` (529) and `camelshotnoise` (870) at `λ = 0.5`; I did not investigate why (the solver
agrees with Clarabel on small meshes, so it is the problem, not an error of the solver that I
know of).

## Interface

```rust
// problem.rs
pub struct ConeProblem<'a> {
    pub mesh: &'a TriMesh,
    pub geometry: &'a Geometry,
    pub lambda: f64,                             // λ ≥ 0
    pub distortion_weights: Option<Vec<f64>>,    // wE ≥ 0 per vertex; None = 1
    pub cone_costs: Option<Vec<f64>>,            // wR ∈ [0, +∞] per vertex; None = 1
    pub lower_bounds: Option<Vec<f64>>,          // lo ≤ 0 per vertex; None = −∞
    pub upper_bounds: Option<Vec<f64>>,          // hi ≥ 0 per vertex; None = +∞
    pub prescribed_cones: Vec<Cone>,             // fixed angle at no cost
}
impl ConeProblem<'_> { pub fn new(mesh, geometry, lambda) -> ConeProblem }

// solve.rs
pub fn solve(problem: &ConeProblem<'_>, options: &SolveOptions<'_>) -> Result<Solution, ConeError>;
pub fn solve_with_boundary(problem, boundary: BoundaryTreatment, options) -> Result<Solution, ConeError>;
pub enum BoundaryTreatment { Fixed /* default: u = 0 */, Free, Penalized, Convex }   // see "The boundary"
pub struct SolveOptions<'a> {
    pub regularization: Vec<f64>,                // γ of Stage I (reference method); default 1, 1e-1, …, 1e-10 (paper_regularization());
                                                 // short_regularization(): 1e-2, 1e-6, 1e-10
    pub regularization_metric: RegularizationMetric,  // Unweighted (default) | Mass
    pub max_regularized_iterations: usize,       // 50
    pub max_exact_iterations: usize,             // 50
    pub max_safeguard_rounds: usize,             // 12
    pub newton_tolerance: f64,                   // 1e-12
    pub certificate_tolerance: f64,              // 1e-9
    pub cone_threshold: f64,                     // 1e-12
    pub merge_cone_clusters: bool,               // false
    pub max_merge_cluster_size: usize,           // 8
    pub max_merge_cluster_diameter: f64,         // +∞ (a length)
    pub linear_solver: LinearSolver,             // Auto | SaddleLdlt | SaddleLu
    pub parallelism: Parallelism<'a>,            // Sequential | Pool(&ThreadPool)
    pub method: SolveMethod,                     // Auto (default) | Reference | FactorOnce
    pub factor_once: FactorOnceOptions,          // tolerance 1e-12, max_working_set None, memory_limit_bytes 4 GiB, max_rounds 200
}
pub fn factor_once_support(problem, boundary) -> Result<Option<Unsupported>, ConeError>;  // None = supported
pub struct Solution {
    // the optimum of the problem (this is what the certificate is about)
    pub measure: Vec<f64>,            // µ per vertex, exactly 0 where there is no cone
    pub log_scale_optimal: Vec<f64>,  // u per vertex
    pub boundary_curvature: Vec<f64>, // ν at boundary vertices, 0 elsewhere (every treatment)
    pub adjoint: Vec<f64>,            // φ per vertex
    pub adjoint_shift: (f64, f64),    // closed meshes: constants that may be added to φ
    pub objective: f64, pub distortion_optimal: f64, pub cone_cost: f64,
    // the cones extracted from it
    pub cones: Vec<Cone>,             // by ascending vertex
    pub log_scale: Vec<f64>,          // u of `cones`
    pub distortion: f64, pub total_cone_angle: f64,
    pub extraction: ExtractionReport,
    pub report: SolveReport,
}
pub struct SolveReport {
    pub status: Status,               // Converged | NotConverged | LinearSolverFailed
    pub stages: Vec<StageReport>,     // reference method only (empty for the factor-once method):
                                      // kind, safeguard, newton_iterations, damped_steps,
                                      // constant_shifts, stop: StageStop, seconds
    pub regularized_iterations: usize, pub exact_iterations: usize,
    pub safeguard_rounds: usize, pub safeguard_iterations: usize,
    pub linear_solver: LinearSolverReport,   // factorizations, solves, refinement, times
    pub method: MethodReport,                // requested, used: MethodUsed, fallback: Option<FallbackReason>,
                                             // factor_once: Option<FactorOnceReport>
    pub timings: Timings,                    // setup, regularized, exact, factor_once, certificate, extraction, total
    pub bound_correction: f64, pub zeroed_angles: usize, pub largest_zeroed_angle: f64,
    pub certificate: Certificate,
    pub message: Option<String>,
    pub components: Vec<ComponentReport>,    // one per connected component (empty if connected)
}
pub enum StageKind { Regularized { gamma }, Exact, Proximal { gamma } }
pub enum StageStop { ActiveSetsUnchanged, Tolerance, Stalled, IterationLimit, Cycle, LinearSolverFailure }

// certificate.rs
pub fn certify(problem, log_scale: &[f64], adjoint: &[f64], measure: &[f64]) -> Result<Certificate, ConeError>;
pub fn certify_with_boundary(problem, boundary, log_scale, adjoint, measure) -> Result<Certificate, ConeError>;
pub struct Certificate {
    pub state_residual: f64, pub adjoint_residual: f64, pub cone_residual: f64,
    pub bound_violation: f64, pub dual_infeasibility: f64, pub gauss_bonnet_residual: f64,
    pub primal_objective: f64, pub dual_objective: f64,
    pub duality_gap: f64, pub complementarity_gap: f64, pub gap_uncertainty: f64,
    pub repair_cost: f64, pub repair_bound_violation: f64, pub weights_positive: bool,
}
impl Certificate {
    pub fn holds(&self, tolerance: f64) -> bool;
    pub fn optimum_enclosure(&self) -> Option<(f64, f64)>;   // proved [D, J + repair_cost], or None
    pub fn suboptimality_bound(&self) -> Option<f64>;        // duality_gap + repair_cost, or None
}

// evaluate.rs: scale factor and adjoint of a given cone configuration
pub fn evaluate_cones(problem, cones: &[Cone], par: Parallelism<'_>) -> Result<ConeEvaluation, ConeError>;
pub struct ConeEvaluation { log_scale, adjoint, distortion, total_cone_angle, gauss_bonnet_defect,
    component_defects }

// reoptimize.rs: the optimal angles of given cones (see "Re-optimizing given cones")
pub fn reoptimize_cones(problem, cone_vertices: &[u32], options: &ReoptimizeOptions<'_>) -> Result<Reoptimization, ConeError>;
pub fn reoptimization_problem(problem, cone_vertices) -> Result<ConeProblem, ConeError>;

// extract.rs: what extraction did
pub struct ExtractionReport { threshold, n_dropped, dropped_angle, merge_requested, n_clusters,
    merged: Vec<MergedCluster>,       // members, cone, diameter
    unmerged: Vec<UnmergedCluster>,   // vertices, angle, diameter, reason: TooLarge | TooWide | BoundsViolated
    cones_differ_from_measure, distortion_change, total_angle_change, n_beyond_two_pi }

// feedback.rs: lower wE where the solver's u is far from a flattening's (see "The feedback loop")
pub struct FeedbackOptions { threshold: f64 /* 4 */, reduced_weight: f64 /* 0 */,
    max_rounds: usize /* 4 */, max_reweighted_fraction: f64 /* 0.01 */ }
pub fn solve_with_feedback<R, E: Display>(problem, boundary, options: &SolveOptions, feedback: &FeedbackOptions,
    realize: impl FnMut(&Solution) -> Result<(Vec<f64> /* r_t per face */, R), E>)
    -> Result<FeedbackSolution<R>, FeedbackError<E>>;    // solution, realized: R, feedback: FeedbackReport
pub struct Feedback;                                     // the rule as a state machine (new, weights, observe,
                                                         // round_failed, finish), for callers with their own solve
pub struct FeedbackReport { rounds: Vec<FeedbackRound>, reweighted, pending, weights, stop: FeedbackStop, seconds }
pub struct FeedbackRound { max_abs_mismatch, rms_mismatch, max_abs_eligible, offsets, n_unmeasured,
    n_over_threshold, reweighted, seconds_solve, seconds_realize, seconds_mismatch }
pub enum FeedbackStop { Stable, RoundLimit, TooManyVertices { wanted, allowed }, NoWeightLeft { component },
    RoundFailed { round, message }, NotContinued }
pub fn scale_mismatch(mesh, geometry, log_scale: &[f64], realized_face_log_scale: &[f64])
    -> Result<ScaleMismatch, ConeError>;                 // δ per vertex, offsets, max, rms

// weights.rs: the curvature-adapted weights of the paper's Sec. 8.1
pub fn feature_size_weights(mesh, geometry, epsilon: f64, convention: FeatureSizeConvention)
    -> Result<FeatureSizeWeights, ConeError>;      // AsDescribed | AsPrinted

pub fn soft_threshold(t: f64, x: f64) -> f64;      // P_t(x), paper Eq. 18
pub enum ConeError { Geometry, NotConnected /* no longer returned */, InvalidLambda, WrongLength, InvalidDistortionWeight,
    InvalidConeCost, InvalidBounds, InvalidPrescribedCone, NoDistortionWeight, Infeasible,
    InvalidCone, InvalidArgument, InvalidOptions, Sparse }
```

Nothing is normalized, clamped or repaired silently. Invalid input is an `Err` whose message
says what is wrong and where; a solve that runs but does not reach a certified optimum is an
`Ok` with `report.status != Converged`.

## Steering

| wish | how |
|---|---|
| pin a cone (vertex and angle) | `prescribed_cones: vec![Cone { vertex, angle }]`; costs nothing, overrides cost and bounds of its vertex |
| a cone must be possible at a vertex, angle free | `cone_costs[v] = 0.0` (a free cone site, paper Sec. 8.1) |
| forbid cones in a region | `cone_costs[v] = f64::INFINITY` there |
| make cones cheaper or dearer in a region | `cone_costs[v]` smaller or larger than 1 |
| make distortion matter more in a region, or not at all | `distortion_weights[v]` larger than 1, or `0.0` |
| positive cones only | `lower_bounds = Some(vec![0.0; nv])` |
| angles in `[−π/2, π/2]` | both bounds |
| re-optimize the angles of a given set of cones | `reoptimize_cones(problem, &vertices, ..)`: `cone_costs = 0` at those vertices and `+∞` elsewhere, solved directly ([re-optimizing given cones](#re-optimizing-given-cones)) |
| keep every angle realizable | bounds inside `(−2π, 2π)`; `extraction.n_beyond_two_pi` counts what is outside |
| ignore needle-like spikes | `distortion_weights[v] = 0.0` at their tips ([needle-like spikes](#needle-like-spikes)) |
| find such vertices from a flattening | `solve_with_feedback(.., realize)` with a flattener in `realize` ([the feedback loop](#the-feedback-loop)); `ddg-uv` does it with `UnwrapOptions::cone_feedback` |
| a closed convex mesh | several free sites, cones forbidden elsewhere, bounds ([closed meshes](#closed-meshes)) |

Every array has one entry per mesh vertex; entries of boundary vertices are validated and have
no effect. On a closed mesh the constraints must allow `Σ µ = 2πχ`; otherwise
`ConeError::Infeasible` says which sum is required and which sums the bounds, forbidden vertices
and prescribed cones allow. With a zero distortion weight somewhere the minimizer need not be
unique; a certified minimizer is still returned.

## What comes back

A `Solution` holds two things that must not be confused.

* **The optimum**: `measure`, `log_scale_optimal`, `adjoint`, `objective`. The certificate in
  the report is computed for exactly this triple as returned, after the rounding-level clean-up
  of the cone angles (an angle that a step ending by tolerance left beyond its bound is moved
  onto it; an angle of magnitude at most `newton_tolerance` is set to zero; both are reported in
  `bound_correction`, `zeroed_angles`, `largest_zeroed_angle`) and after the adjoint constant of a
  closed mesh was chosen. `u` is not recomputed after the clean-up, so what the clean-up changed
  shows in the certificate's `state_residual`. `certify(problem, &log_scale_optimal, &adjoint,
  &measure)` reproduces the report's certificate bit for bit (tested).
* **The extracted cones**: `cones`, `log_scale`, `distortion`, `total_cone_angle`. A vertex is
  a cone if `|µ_i| > cone_threshold`. With `merge_cone_clusters` every cluster of cones is
  replaced by one cone with the summed angle at the cluster's Fréchet mean vertex. A cluster
  is a set of cones **of the same sign** connected through mesh edges between its members; it
  is left alone if it has more than `max_merge_cluster_size` cones or is wider than
  `max_merge_cluster_diameter`. If extraction changed anything
  (`extraction.cones_differ_from_measure`), the cones are deliberately not optimal, `log_scale`
  and `distortion` are recomputed for them, nothing certifies them, and `extraction` lists every
  difference. What merging cost is always visible: `distortion_optimal` is the distortion
  before, `distortion` after, `extraction.distortion_change` the difference.

**The range of the returned angles.** Every angle lies within the bounds given for its vertex,
exactly, also after merging (a merged cone is only placed at a vertex whose bounds admit it).
Without bounds nothing limits an angle, and the linearized problem does not know that a cone
angle must be below `2π`: a free cone site on a closed convex mesh takes up to `4π`, the base
of a needle-like spike gets angles below `−2π`, and a merged cone carries the sum of its
cluster. `extraction.n_beyond_two_pi` counts the returned angles that are `≥ 2π` (no such cone
exists) or `≤ −2π` (an angle sum of `4π` or more, which no layout realizes without overlap).
Bounds inside `(−2π, 2π)` exclude both.

## The certificate

The first six fields are residuals of the optimality system
`(L u)_i = µ_i − Ω_i`, `(L φ)_i = −wE_i m_i u_i`, `φ_i ∈ ∂(λ wR_i|·| + ι_[lo_i,hi_i])(µ_i)`,
which is sufficient for a global minimum (machine-checked in `lean/proofs`:
`optimality-system-sufficient`). The gap fields are computed as follows, with
`w_i = wE_i m_i`, `t_i = λ wR_i` and sums over the unknown vertices:

```
r_i = (L u)_i − µ_i + Ω_i
φ̂_i = φ_i clipped to the domain of g_i*          (φ̂_i ≤ t_i if hi_i = +∞, φ̂_i ≥ −t_i if lo_i = −∞)
ũ_i = −(L φ̂)_i / w_i                             where w_i > 0
J   = ½ Σ w_i u_i² + Σ t_i|µ_i|                                      primal_objective
D   = −½ Σ_{w_i>0} (L φ̂)_i²/w_i + Σ φ̂_i Ω_i − Σ g_i*(φ̂_i)           dual_objective
J − D = ½ Σ_{w_i>0} w_i (u_i − ũ_i)² + Σ (t_i|µ_i| + g_i*(φ̂_i) − φ̂_i µ_i)     complementarity_gap
        + Σ_{w_i=0} (L φ̂)_i u_i − Σ φ̂_i r_i                                    |·| ≤ gap_uncertainty
```

The last two lines are an algebraic identity for any numbers `u, φ, µ` that vanish on the
boundary (`duality-gap-inexact`; checked on random non-solutions in `tests/api.rs`). The dual
point `(ũ, φ̂)` satisfies the adjoint equation by construction and is dual feasible exactly.

What is proved, in `lean/proofs` (`DdgProofs/Certificate.lean`), for exact real arithmetic:

* The identity above (`objective_sub_dualObjectiveOfAdjoint`).
* **Lower bound** (`dual_lower_bound_of_clip`): if all `wE_i m_i > 0`, `D` is at most the
  objective of every feasible pair, whatever the residuals of the triple are.
* **Backward statement** (`backward_certificate`): if moreover `µ` respects its bounds, then
  `(u, µ)` is exactly feasible for the problem whose angle defects are `Ω − r` instead of `Ω`
  (`|r_i| ≤ state_residual`) and is suboptimal *there* by at most `complementarity_gap ≥ 0`.
* If `r = 0`: `objective − optimum ≤ duality_gap = complementarity_gap`.
* **Enclosure from a computed triple** (`certificate_enclosure`, `certificate_of_computed`): if all
  `wE_i m_i > 0` and `µ_i + r_i` is admissible at every vertex (inside `[lo_i, hi_i]`; equal to the
  fixed angle at a forbidden or prescribed vertex), then `(u, µ + r)` is exactly feasible and
  `D ≤ optimum ≤ J(u, µ + r) ≤ J(u, µ) + λ Σ wR_i |r_i|`; its objective exceeds the optimum by at
  most `duality_gap + repair_cost`. The certificate reports the numbers this needs:
  `repair_cost = λ Σ wR_i |r_i|` (penalized unknown vertices; a fixed vertex and a free site
  cost nothing), `repair_bound_violation` (the largest distance of `µ_i + r_i` from the
  admissible set, and `max(|u|, |µ|)` at vertices that are not unknowns, which must be exactly
  `0`) and `weights_positive`. `optimum_enclosure()` returns `Some((D, J + repair_cost))` exactly
  when `weights_positive` and `repair_bound_violation == 0` (and the values are finite), else
  `None`; `suboptimality_bound()` is then `duality_gap + repair_cost`. In particular: at a
  penalized vertex without bounds every residual is admissible; at a forbidden or prescribed
  vertex the residual must be exactly `0`, which a computed solution has only by luck (it is
  rounding noise: on problems with such features even a converged solve usually gets `None`,
  215 of 366 in `tests/enclosure.rs`), and then nothing is proved from above; likewise for an
  angle at its bound with a residual pushing it out. Without positive weights the lower bound is
  not proved. Tested: on 339 deliberately perturbed triples (`u`, `µ`, `φ` moved by up to 0.1)
  the enclosure contained the optimal value of the generic convex solver every time. All of this
  is about exact arithmetic on the `f64` values of the triple; rounding in `D`, `J`, `r` is not
  covered.

What is not:

* Nothing bounds the suboptimality of the returned pair in the *given* problem: it satisfies
  the state equation only up to `state_residual`, and an infeasible pair can have an objective
  below the optimum, so `duality_gap` can be negative. It is `complementarity_gap` (a sum of
  nonnegative terms) shifted by at most
  `gap_uncertainty = Σ_{w_i=0} |(L φ̂)_i u_i| + Σ |φ̂_i r_i|`. Read the gap as
  "`complementarity_gap`, give or take `gap_uncertainty`". (The proved enclosure
  `D ≤ optimum ≤ J + λ Σ wR_i |r_i|` needs `µ + r` to be admissible at every vertex; see the
  next item and `Certificate::optimum_enclosure`.)
* **Where `wE_i = 0`** the dual problem requires `(L φ̂)_i = 0` exactly, which is not enforced by
  construction. Such a vertex contributes nothing to `D`; its violation is in `adjoint_residual`
  and, times `u_i`, in `gap_uncertainty`. `D` is then a lower bound only up to
  `Σ_{w_i=0} (L φ̂)_i v_i`, which depends on the competitor `v`.
* Everything is evaluated in `f64`; no statement is made about rounding.

`Status::Converged` means `certificate.holds(certificate_tolerance)`: the six residuals are at
most the tolerance, and `|duality_gap|`, `|complementarity_gap|` and `gap_uncertainty` are at
most `tolerance · max(1, |J|)`. Nothing else decides the status.

## Disconnected meshes

Every function that takes a mesh (`solve`, `solve_with_boundary`, `certify`,
`certify_with_boundary`, `evaluate_cones`, `factor_once_support`, `reoptimize_cones`) accepts a
mesh with several connected components and returns the **exact solution of the canonical
problem on the whole mesh**: the masses divided by the total area `A` of the whole mesh, the
given `λ`, one additive constant of `u` per closed component.

**The problem separates.** Vertices of different components share no edge, so `L` (and its
block on the unknowns) is block diagonal; the objective `½ Σ wE m u² + λ Σ wR |µ|` is a sum over
the components; every constraint involves one component. So `(u, µ)` is a minimizer iff its
restriction to every component minimizes that component's part of the objective under that
component's constraints, and the optimality systems separate the same way
(`isMinimizer_sum_iff`, `kkt_sum_iff` in `lean/proofs/DdgProofs/Components.lean`,
`component-separation`). In particular every closed component has its own additive constant of
`u` (the competitors' constants are independent), and its own Gauss–Bonnet constraint
`Σ_C µ = 2πχ_C` (its block of `L` annihilates its constants; `sum_cone_eq_sum_defect_inl`). A
closed component whose bounds, forbidden vertices or prescribed cones cannot sum to its `2πχ` is
`ConeError::Infeasible { component: Some(c), .. }`, with the component's `χ` and sums; a closed
component without distortion weight is `NoDistortionWeight { component: Some(c) }`.

**Which normalization a component gets.** The part of the objective that belongs to component
`c` has the masses of the whole problem, `a_i / A`:

```
½ Σ_C wE (a_i/A) u² + λ Σ_C wR |µ|  =  (A_c/A) · [ ½ Σ_C wE (a_i/A_c) u² + (λ A/A_c) Σ_C wR |µ| ]
```

A positive factor does not change the minimizers, so the component solved **alone**, with its
masses normalized by its own area `A_c`, has the same `u` and `µ` exactly when it gets
`λ_c = λ · A / A_c`, and its adjoint is then `A / A_c` times the adjoint of the whole problem
(`isMinimizer_normalizeMass_alone_iff`, `kkt_normalizeMass_alone_iff`; tag
`component-lambda-rescaling`; the special case of one component is `lambda-area-normalization`).
The solver uses the first form: each component keeps the masses of the whole mesh and the given
`λ` (`Instance::part`), so no `λ` is rescaled and the adjoint needs no conversion;
`ComponentReport::equivalent_lambda` reports `λ · A / A_c` for whoever wants to solve a
component alone. The consequence is that a component's cones depend on the rest of the mesh: a
small component counts little in the total distortion and gets fewer cones than alone (in
`tests/components.rs` a disk and a copy scaled by 0.2: 5 cones on each when solved alone, 5 and 0
in the union at `λ = 0.004`).

**Interface.** The `Solution` is that of the whole mesh. `SolveReport::components` has one
`ComponentReport` per component (in the order of `TriMesh::face_components`;
`TriMesh::vertex_components` maps vertices to them) when there is more than one, and is empty
for a connected mesh:

```rust
pub struct ComponentReport {
    pub component: u32, pub n_vertices: usize, pub n_faces: usize, pub n_unknowns: usize,
    pub closed: bool,                 // every vertex an unknown: own constant of u, own 2πχ
    pub area_fraction: f64,           // A_c / A
    pub equivalent_lambda: f64,       // λ · A / A_c: the same result when solved alone
    pub objective: f64, pub distortion_optimal: f64, pub cone_cost: f64,   // its part
    pub adjoint_shift: (f64, f64),    // constants admissible on this component alone
    pub report: SolveReport,          // its status, certificate (of its part), stages, method, times
}
```

The certificate in the top-level report is computed from the assembled triple for the whole
problem (`certify` reproduces it bit for bit), and it alone decides `status`; its
`gauss_bonnet_residual` is the largest over the closed components. The other top-level fields
aggregate: stages one component after another, counts and phase times summed (`timings.total`
is wall time), `bound_correction` and `largest_zeroed_angle` maxima, `method.used` is
`FactorOnce` only if every component with unknowns used it, `method.fallback` the first
component's reason, `adjoint_shift` the constants admissible on every component at once
(`(0, 0)` unless all are closed). Extraction runs on the whole mesh (a cluster never spans two
components); if it changes cones, `log_scale` is recomputed per component.
`ConeEvaluation::component_defects` gives the Gauss–Bonnet defect of every component.

**Solving.** Each component is a connected problem with everything of the other sections: the
factor-once method factorizes the Laplacian block of that component (a closed component pins
one of its own vertices), with its own working set and limits (`factor_once_support` answers per
component, naming vertices of the whole mesh); the reference method runs per component. With
`SolveOptions::parallelism = Pool(..)` the components run concurrently on the pool, each on one
thread, so the result does not depend on the number of threads (tested bit for bit); a mesh with
one large component and small ones then loses the inner threading of the large one (worth 1.4×
to 2× only with hundreds of cones, see [performance](#performance-of-the-factor-once-method)).
A connected mesh takes the path it took before, bit for bit (`examples/fingerprint.rs`: 44 of 44
hashes unchanged).

**Inputs at the edge.** Isolated vertices cannot reach the solver: `TriMesh::new` rejects a
vertex that no face references (`MeshError::UnreferencedVertex`; `compact_vertices` removes
them). A component that is a single triangle has three boundary vertices and, with `u = 0` on the
boundary, no unknowns: it gets `u = µ = 0`, contributes nothing, and its report has
`n_unknowns = 0` and no stages (with any `SolveMethod`, also `FactorOnce`); with an optimized
boundary it is a problem with three unknowns like any other.

**Verified against** (`tests/components.rs`): the union against separate solves with
`λ · A / A_c` (two closed meshes, a closed one and one with boundary, three components, every
boundary treatment; `µ`, `u` to `10⁻⁷`, the adjoint up to the factor `A_c / A`, the objective
as the area-weighted sum, per-component reports); 40 random disconnected problems with every
steering feature against the generic convex solver Clarabel (40 of 40 compared); steering per
component; random relabeling of vertices, faces and corners; both methods and a thread pool;
`evaluate_cones`, merging, `factor_once_support`.

**`ddg-uv` makes a different choice.** `crates/ddg-uv/src/prepare.rs` cuts the mesh along its
seams and splits it into pieces (its connected components after cutting), and
`crates/ddg-uv/src/cones.rs` solves every piece as a `ConeProblem` of its own with the user's
`λ`, i.e. with the masses normalized by the **piece's** area. By the identity above, that is the
canonical problem of the whole (cut) mesh with `λ · A_c / A` on piece `c` instead of `λ`: every
piece gets the cones it would get alone, independent of the other pieces, where the canonical
problem of the whole would give a small piece fewer cones (and a large one slightly more). This
is a defensible choice for UV charts (a trade-off per chart, unaffected by unrelated parts of the
model), and the seams make the pieces a different problem from the uncut mesh anyway. It was
not changed.

## Re-optimizing given cones

```rust
pub fn reoptimize_cones(problem: &ConeProblem<'_>, cone_vertices: &[u32], options: &ReoptimizeOptions<'_>)
    -> Result<Reoptimization, ConeError>;
pub fn reoptimization_problem<'a>(problem: &ConeProblem<'a>, cone_vertices: &[u32]) -> Result<ConeProblem<'a>, ConeError>;
pub struct ReoptimizeOptions<'a> { pub boundary: BoundaryTreatment, pub method: ReoptimizeMethod /* Auto | DirectKkt | Solve */,
                                   pub solve: SolveOptions<'a> }
pub struct Reoptimization { pub cones: Vec<Cone> /* the given vertices, given order, optimal angles */,
    pub measure, pub log_scale, pub adjoint, pub boundary_curvature: Vec<f64>,
    pub distortion: f64, pub objective: f64, pub total_cone_angle: f64, pub report: ReoptimizeReport }
pub struct ReoptimizeReport { pub status: Status, pub method: ReoptimizeMethodUsed /* DirectKkt | Solve */,
    pub certificate: Certificate, pub linear_solver: LinearSolverReport,
    pub components: Vec<ReoptimizeComponentReport>, pub seconds: f64, pub message: Option<String> }
pub const REDUCED_SITES_LIMIT: usize = 32;
```

Given the cone **positions**, it finds the cone **angles** of least distortion: the canonical
problem with cone cost `0` at the given vertices and `+∞` at every other interior vertex
(`reoptimization_problem` builds it; the distortion weights, bounds and prescribed cones of
`problem` are kept, the cone costs of `problem` matter only at boundary vertices with
`Penalized`). `λ` then plays no role (except for `Penalized`). Works with every boundary
treatment and on disconnected meshes. Errors (`InvalidCone` with the position in the list): a
vertex out of range, listed twice, on the boundary, or carrying a prescribed cone. An empty list
is allowed (no cones: on a mesh with boundary `µ = 0`; on a torus too; on a sphere it is
`Infeasible`). On a closed component the given cones must be able to carry `2πχ`: one unbounded
cone suffices (a single cone on a sphere gets `4π`), bounds that cannot is `Infeasible` naming
the component.

**Method.** Without the `ℓ¹` term the problem is a quadratic program. If no bound is active,
its optimality system is linear: `W u + A φ = 0`, `(A u)_i = a_i − Ω_i` at the fixed vertices,
`φ_i = 0` at the given cones (a free site's subdifferential is `{0}`). That is the Newton system
of the reference method with every piece known beforehand (`saddle-point-newton-system`), and a
solution is a global minimizer (`optimality-system-sufficient`). Two ways to solve it, both
exact:

* **direct** (`DirectKkt`): one sparse `LDLᵀ` factorization of the `2n × 2n` saddle-point matrix
  and one solve with refinement; `µ_i = (A u)_i + Ω_i` at the cones. The cost does not depend on
  the number of cones. On a closed component the pinned rows `φ_i = 0` fix the constant of `φ`,
  and the adjoint equation summed gives `Σ w u = 0`, the optimal constant of `u`.
* **reduced**: the generic solver's factor-once method on `reoptimization_problem`: one Cholesky
  factorization of the `n × n` Laplacian block, one solve per free site, and the dense system of
  the free sites, solved exactly by its active-set method (bounds included).

`Auto` takes the reduced system up to `REDUCED_SITES_LIMIT = 32` free sites per component (the
given cones, plus the boundary vertices with an optimized boundary) and the direct solve beyond.
The direct solution is accepted only if every angle is inside its bounds and its certificate
holds; otherwise (an active bound), and with `Penalized` (an `ℓ¹` term remains) or a closed
component whose angles are all fixed (a singular system), the component goes to the generic
solver. Every component's report says which method ran and why
(`ReoptimizeComponentReport::direct_not_used`), and the result carries the certificate of the
whole problem of `reoptimization_problem`.

Single-thread times on the dataset (`cargo run --release -p ddg-optimal-cones --example
reoptimize -- --repeat 3`, fastest of 3, alone on the machine; "solve" is `solve` on
`reoptimization_problem`, i.e. the factor-once method): the reduced system is 1.1 to 2.6 times
faster than the direct solve up to 29 cones (`body` 0.025 s against 0.066 s, forced
`DirectKkt`), the two are even at 36 (`rhino` 0.075 / 0.081 s) and 45 (`brain`, 153 k vertices,
1.06 / 1.00 s), and the direct solve wins by 6.6 times at 439 cones (`cow01`, 0.052 / 0.342 s)
and 34 times at 13 836 (`cow00001`, 0.051 / 1.72 s).

**The published cone files** (`examples/reoptimize.rs`; asserted in
`tests/reoptimize.rs::published_cones_reoptimized`): at the published positions, the
re-optimized angles never give more distortion than the published angles. `D = ½ Σ m u²` with
`λ = ½` irrelevant; "published angles" is `evaluate_cones` with the file's angles (on a closed
mesh shifted by a common constant to sum to exactly `4π`, a change of at most `3·10⁻⁴`);
"boundary optimized" is `BoundaryTreatment::Free` (the file's angles with the boundary
optimized for them, against `reoptimize_cones` with `Free`).

| file (cones) | D, published angles | D, re-optimized | ratio | largest angle change | D, published, boundary optimized | D, re-optimized, boundary optimized | ratio | `reoptimize_cones` (s) | method | `solve`, costs 0/∞ (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| airplane (13) | 2.5099e-1 | 7.6469e-2 | 0.3047 | 6.3e-1 | 6.2628e-2 | 6.2385e-2 | 0.9961 | 0.004 | Solve | 0.004 |
| axes, closed (4) | 3.9905e-1 | 3.9905e-1 | 1.0000 | 1.3e-4 | – | – | – | 0.009 | Solve | 0.010 |
| body (10) | 2.0317e-1 | 2.0092e-1 | 0.9889 | 8.6e-2 | 1.9649e-1 | 1.8291e-1 | 0.9309 | 0.035 | Solve | 0.030 |
| brain, closed (45) | 1.7153e-2 | 1.7153e-2 | 1.0000 | 2.6e-3 | – | – | – | 1.134 | DirectKkt | 1.070 |
| bump (5) | 9.9707e-3 | 9.9707e-3 | 1.0000 | 2.4e-4 | 9.9198e-3 | 9.7577e-3 | 0.9837 | 0.006 | Solve | 0.006 |
| bunny, closed (4) | 1.1511e-1 | 1.1485e-1 | 0.9977 | 9.0e-2 | – | – | – | 0.040 | Solve | 0.038 |
| camel (6) | 5.0121e-2 | 2.0273e-2 | 0.4045 | 3.6e-1 | 1.4915e-2 | 1.4915e-2 | 1.0000 | 0.027 | Solve | 0.027 |
| camel_noise1 (7) | 6.6811e-2 | 2.7369e-2 | 0.4097 | 3.9e-1 | 2.1570e-2 | 2.1570e-2 | 1.0000 | 0.027 | Solve | 0.026 |
| camel_noise2 (7) | 2.6944e-1 | 1.1042e-1 | 0.4098 | 6.4e-1 | 9.6985e-2 | 9.6983e-2 | 1.0000 | 0.025 | Solve | 0.025 |
| camelbad (6) | 4.9319e-2 | 1.9947e-2 | 0.4044 | 3.5e-1 | 1.4574e-2 | 1.4573e-2 | 1.0000 | 0.028 | Solve | 0.028 |
| camelbadnonuniform (6) | 4.7527e-2 | 2.0221e-2 | 0.4255 | 3.5e-1 | 1.4966e-2 | 1.4964e-2 | 0.9998 | 0.020 | Solve | 0.020 |
| camelshotnoise (6) | 8.4785e4 | 8.4712e4 | 0.9991 | 3.9e1 | 8.4706e4 | 8.4696e4 | 0.9999 | 0.030 | Solve | 0.031 |
| chair, closed (14) | 6.5701e-3 | 6.5701e-3 | 1.0000 | 7.8e-4 | – | – | – | 0.017 | Solve | 0.016 |
| cow00001 (13836) | 1.6034e-6 | 5.2731e-10 | 0.0003 | 2.5e-3 | 5.2362e-7 | 4.6731e-10 | 0.0009 | 0.052 | DirectKkt | 1.791 |
| cow01 (439) | 1.5395e-2 | 5.3820e-3 | 0.3496 | 3.5e0 | 1.5098e-2 | 5.3226e-3 | 0.3525 | 0.052 | DirectKkt | 0.342 |
| cow1 (2) | 5.0953e-1 | 4.7997e-1 | 0.9420 | 4.1e-1 | 4.4975e-1 | 4.4975e-1 | 1.0000 | 0.038 | Solve | 0.037 |
| cow5e-1 (5) | 7.4871e-1 | 2.6767e-1 | 0.3575 | 1.0e0 | 1.8475e-1 | 1.8475e-1 | 1.0000 | 0.040 | Solve | 0.039 |
| cow75e-2 (3) | 3.8681e-1 | 3.6560e-1 | 0.9452 | 4.1e-1 | 3.5835e-1 | 3.5835e-1 | 1.0000 | 0.039 | Solve | 0.038 |
| face (2) | 2.1981e-3 | 2.1981e-3 | 1.0000 | 4.6e-5 | 2.0880e-3 | 2.0824e-3 | 0.9974 | 0.030 | Solve | 0.029 |
| hand (29) | 7.4499e-3 | 6.5589e-3 | 0.8804 | 1.2e-1 | 6.4206e-3 | 6.4205e-3 | 1.0000 | 0.008 | Solve | 0.008 |
| mannequin (8) | 9.1555e-3 | 9.1555e-3 | 1.0000 | 3.6e-4 | 8.8263e-3 | 7.3326e-3 | 0.8308 | 0.016 | Solve | 0.015 |
| rhino (36) | 1.9547e-2 | 1.8860e-2 | 0.9648 | 2.2e-1 | 1.8585e-2 | 1.8582e-2 | 0.9999 | 0.076 | DirectKkt | 0.085 |
| tree (8) | 2.1784e-2 | 2.1784e-2 | 1.0000 | 2.7e-3 | 2.1765e-2 | 2.1658e-2 | 0.9950 | 0.046 | Solve | 0.045 |

The ratios sort the files as [the published dataset](#the-published-dataset) did: where the
published angles are re-optimized angles with `u = 0` on the boundary (`axes`, `brain`, `bump`,
`chair`, `face`, `mannequin`, `tree`) the first ratio is 1.0000 (`body` 0.989 and `bunny` 0.998,
whose angles were matched less closely there); where they belong to optimized
boundary values (the camels, the sparse cows, `hand`, `rhino`, `airplane`) the second ratio is
0.9961 to 1.0000 and the first 0.30 to 0.97; `cow01` and `cow00001` are penalized optima, whose
angles are shrunk by `λ`, so re-optimizing lowers their distortion by 3 and 3 000 times;
`camelshotnoise` is dominated by its spikes ([needle-like spikes](#needle-like-spikes)).

**Verified against** (`tests/reoptimize.rs`): 164 random problems (every topology, weights,
bounds, prescribed cones, the four boundary treatments) against `solve_with_boundary` on
`reoptimization_problem` with both methods, and with `Auto` and `Solve` (`µ`, `u` to `10⁻⁷`,
objective to `10⁻⁹`); 80 of them were solved directly, 84 had an active bound; 40 against the
generic convex solver (38 compared, the rest infeasible for both); moving one cone by one vertex
on an ellipsoid (the distortion lies between that of the cones with both vertices and that
without the moved one, and changes by less than 25 %); the edge cases above; a disconnected mesh
(the angles of each component equal those of the component alone, whatever its normalization);
the dataset.

## Algorithm

The paper's (Fig. 12). **Stage I** solves the regularized problems for `γ = 1, 10⁻¹, …, 10⁻¹⁰` (a shorter schedule
is selectable; see [the schedule](#the-regularization-schedule)),
each started from the previous solution; **Stage II** solves the exact optimality system
started from the last regularized solution. Both are semismooth Newton iterations in which a
step chooses a *piece* for every vertex (no cone, positive cone, negative cone, at the lower
bound, at the upper bound, fixed) and solves the linear optimality system of these pieces.
Every such system is

```
W u + A φ = 0                         W = diag(wE m), A = L on the unknown vertices
(A u)_i                   = µ_i − Ω_i          vertex with a fixed angle on its piece
(A u)_i − k_i (φ_i − o_i) = b_i − Ω_i          active vertex of a regularized problem
φ_i = ±λ wR_i                                   active vertex of the exact problem
```

and is solved by a sparse `LDLᵀ` factorization of the symmetric `2n × 2n` matrix (unknown `u_i`
eliminated immediately before `φ_i`), with iterative refinement on the componentwise backward
error, and by `LU` where a distortion weight is zero or `LDLᵀ` is not accurate enough.

## The factor-once method

The reference method refactorizes a `2n × 2n` saddle-point matrix at every Newton step (12 to 57
times per solve). The Laplacian never changes, and
the optimum has few cones, so a second method factorizes it once. It is the default
([`SolveMethod::Auto`](#choosing-the-method)); the reference method is the fallback.

**The reduced problem.** Let `A` be the Laplacian block of the unknowns, `W = diag(wE m)` (all
`> 0`), `S` the free cone sites, `F` the sites with a fixed angle (prescribed, or forbidden:
angle 0) with angles `a_F`, `Z` the columns `z_j = A⁻¹ e_j` of the free sites and
`w0 = A⁻¹ (Ω − a_F)`. Then `u = Z µ_S − w0`, and the canonical problem is

```
minimize   ½ µ_Sᵀ G µ_S + qᵀ µ_S + Σ_j t_j |µ_j|     subject to  lo_j ≤ µ_j ≤ hi_j   [and Σ µ_S = c]
G = Zᵀ W Z,    q = −Zᵀ W w0,    t_j = λ wR_j                                [formula: reduced-lasso-problem, reduced-gram-matrix]
```

(`G` is positive definite when all weights are positive; the bracket is for a singular `A`, below.)
With `y = A⁻¹ W u`, the gradient of the distortion is `∂/∂µ_j = y_j`, and `µ` is optimal iff for
every free site `j` (and the multiplier `ν` of the sum constraint, `0` without it)
`−(y_j + ν) ∈ ∂g_j(µ_j)`; the adjoint is `φ = −(y + ν)`
[formula: reduced-optimality-conditions]. So one solve for `u` and one for `y` check the
optimality conditions of **all** sites.

**The loop.** Keep a *working set* of sites whose columns `z_j` have been computed (one solve
each; independent, so they run on the thread pool in chunks of 4 to 16 right-hand sides) and the
Gram matrix `G_WW` (a product of the stored columns, row-chunked over the pool). A round:

1. solve the problem restricted to the working set with a dense active-set method (an ℓ¹
   method of Osborne, Presnell and Turlach, extended here to box constraints: every variable is
   zero, free with a sign, or clamped at a bound; a Cholesky factor of `G_BB` of the free ones
   is extended and shrunk in `O(a²)`);
2. compute `u` (from the stored columns, or one solve) and `y` (one solve) and re-anchor the
   dense gradient to the exact one (one step of iterative refinement; without it the clustered
   cones of a noisy mesh leave residuals of `10⁻⁷` in the Gram matrix's badly conditioned
   directions);
3. every site that violates its optimality condition by more than `tolerance · λ` joins the
   working set: the worst ones, but at most one per neighbourhood of 3 edges, up to half the
   number of nonzero angles (at least 16) per round, because neighbouring sites have nearly
   dependent columns and adding both is wasted.

It ends when no site outside the working set violates its condition. There is no convergence
proof and no bound on the number of rounds; the certificate is what makes the result safe, and
the method gives up (below) rather than iterate forever.

**Closed meshes, and all the boundary treatments.** `A = L` is singular exactly when every
vertex is an unknown (a closed mesh, or a boundary treatment other than `Fixed`: then the
problem *is* the closed-mesh problem with one site per vertex, so this path serves them all).
One free site, the *pin* `p`, loses its equation and its `u` is set to zero before shifting:
`A_p` (without row and column `p`) is positive definite, `z_j = P A_p⁻¹ e_j` (zero for `j = p`),
`w0 = P A_p⁻¹ (Ω − a_F)_p̂`, where `P x = x − (Σ w_i x_i / Σ w_i) 1` chooses the additive constant of
`u` that minimizes the distortion. The equation of the pin is the sum constraint
`Σ_S µ = Σ (Ω − a_F)` [formula: reduced-sum-constraint-pin], which the dense method keeps with a
multiplier `ν`; the pin has a zero row and column in `G` and starts free with the whole sum.
The pin is any free site whose bounds allow the whole sum; if there is none (every angle capped
below the sum, as the UV pipeline does on a sphere: `0 ≤ µ ≤ 0.95·2π`, sum `4π`) the method
starts from a feasible point with the sites of largest curvature clamped at their bounds, as many as
needed (three there), and the pin carrying the rest.

**What is supported.** Everything in the canonical problem with *positive* distortion weights:
per-vertex cone costs (free sites `wR = 0` and forbidden ones `wR = ∞` included), bounds `lo`
and `hi` (finite or not, positive-only and negative-only angles), prescribed cones (they move
to the right-hand side), meshes with boundary (`u = 0` there), closed meshes, and the boundary
treatments of `solve_with_boundary`. Each was checked against the reference method
(`tests/factor_once.rs`: 969 random problems of every topology, steering and boundary treatment, and 5 closed spheres whose angles are all capped) and the
rest of the suite was run with the method forced.

**What falls back** (`MethodReport::fallback = Unsupported(..)`, nothing is lost but the
speed-up):

| situation | why |
|---|---|
| a distortion weight `wE_i m_i = 0` at an unknown (boundary vertices excluded with `Fixed`) | `G` is singular and the minimizer may not be unique; I have no treatment that is shown to be correct |
| `λ = 0` | every site is free, there is nothing sparse; the tolerance is relative to `λ` |
| a singular `A` and the bounds of the free sites do not add up to the required sum `Σ (Ω − a_F)` (up to the feasibility tolerance of the problem) | no feasible starting point |
| an unknown set that is empty | nothing to solve |

**Giving up.** In [`SolveMethod::Auto`] the method stops and the reference method solves the
problem when one of these happens, and the report says which
([`FallbackReason`]; the failed attempt's counters stay in `report.method.factor_once`):

| reason | limit | default |
|---|---|---|
| `WorkingSetLimit` | the working set would exceed `max_working_set` | `1300 (n/15000)^0.2 · (1 + (threads−1)/7)`, threads ≤ 8 (justified under [performance](#performance-of-the-factor-once-method)) |
| `MemoryLimit` | the stored columns (`8 n` bytes each) and the dense matrices would exceed `memory_limit_bytes` | 4 GiB (`brain` needs 0.7 GB at `λ = ½`) |
| `IterationLimit` | `max_rounds` rounds, or the dense method hit its iteration cap three rounds in a row | 200 |
| `CertificateFailed` | the result does not satisfy `certificate.holds(certificate_tolerance)` (the failed certificate is in the reason) | `10⁻⁹` |
| `FactorizationFailed` | the factorization or a solve failed | – |

With `SolveMethod::FactorOnce` there is no fallback: a problem outside the support is
`Err(ConeError::InvalidOptions)`; the other reasons give a result with `status = NotConverged`
whose `report.message` and `report.method.fallback` say why. `SolveMethod::Reference` never tries
the fast method. In every case the result carries the certificate computed by `certify`, exactly
as for the reference method; the factor-once method has no say in whether it is `Converged`.

**Choosing the method.** `SolveOptions::method` (`Auto`, `Reference`, `FactorOnce`),
`SolveOptions::factor_once` (the limits above, and `tolerance`), and `factor_once_support(problem,
boundary)` to ask whether a problem is supported. `report.method` says what happened:
`used` (the method that produced the result), `fallback` (why the factor-once method did not,
or why it ended without a certified result), `factor_once` (rounds, solves, working set, nonzero
angles at the end and at the peak, memory estimate, counters of the dense method, times by phase).
`report.stages` is empty when the factor-once method produced the result (there are no Newton
stages), and `report.timings.factor_once` is its time.

**Threads.** `SolveOptions::parallelism` with a pool is used for the column solves (chunks of
right-hand sides), the Gram product (row blocks) and the evaluation of `u` from the stored
columns, through `ThreadPool::map` of `ddg-core`; the factorization and single solves run on the
calling thread (threads do not help them, see the `ddg-core` README), faer is called with
`Par::Seq`, and nothing touches the global `rayon` pool.

**Tolerance.** `FactorOnceOptions::tolerance = 10⁻¹²` (relative to `λ`) is not negotiable:
in the prototype a looser `10⁻⁸` left cone residuals of `5·10⁻⁹` on `camelshotnoise` (452 cones,
many nearly zero) and failed the certificate at `10⁻⁹`, `10⁻¹⁰` missed one near-zero cone.
`10⁻¹²` costs nothing. (`tests/factor_once.rs` shows a loose `10⁻⁴` failing the certificate on a
noisy sphere and the reference method taking over.)

**Not the same bits.** The two methods return the same optimum up to what their certificates
allow, not the same bits. `tests/factor_once.rs` compares them with *derived* tolerances: the
certificate's `suboptimality_bound g` bounds `J − J*`; the objective is `½ ‖u − u*‖²_W`-strongly
convex in `u`, so `‖u_a − u_b‖_W ≤ √(2 g_a) + √(2 g_b)`, and `µ = L u + Ω` gives
`|Δµ_i| ≤ Σ_j |L_ij| |Δu_j|`. On those 969 problems the largest observed ratio to the allowed
difference was `1.8·10⁻⁷`, the largest `|Δµ|` `2.2·10⁻¹¹`. On the dataset at `λ = ½` (below)
`|Δµ| ≤ 4·10⁻⁸` (`camelshotnoise`; `4·10⁻⁷` on `brain` with 8 threads), `|Δu| ≤ 9·10⁻⁹`
(`8·10⁻⁸` on `brain` with 8 threads), relative difference of the objective `≤ 3·10⁻¹²`.

## The regularization schedule

The default stays the paper's `γ = 1, 10⁻¹, …, 10⁻¹⁰` (`SolveOptions::paper_regularization()`).
The shorter `γ = 10⁻², 10⁻⁶, 10⁻¹⁰` (`SolveOptions::short_regularization()`) was evaluated as a
replacement and **rejected**. Stage I only initializes Stage II, so the optimum does not depend on
the schedule, and the evidence is that it does not: with the reference method on every dataset file
at `λ = ½` (23 files, `brain` included) and on 17 meshes at `λ = 0.1, 0.02, 0.005`
(`cargo run --release -p ddg-optimal-cones --example speedup -- --schedules`), both schedules
converged, certified, to **bit-identical** `µ` and `u` (largest difference `0`: the exact stage
ends on the same active sets, hence on the same linear system). Factorizations fell in all but one
of these cases, `1.2×` on `brain` (50 → 41, 32.8 → 27.0 s) to `3.9×`. The whole test suite also
passed with the short schedule as the default (reference method forced, debug and release,
including the release dataset tests; and in automatic mode).

What decided against it is the benchmark of `ddg-uv` (`cargo bench -p ddg-uv`, before and after),
which asks for positive-only cone angles capped below `2π` (`lo = 0`, `hi = 0.95·2π`), and no test
or dataset run above has such bounds. On that problem the short schedule needs more
factorizations, not fewer (`--schedules --cap 5.969 --lambda 0.1`, reference method):

| mesh | paper schedule | short schedule `10⁻², 10⁻⁶, 10⁻¹⁰` |
|---|---:|---:|
| bunny | 68 factorizations, 3.1 s | 171, 7.6 s (0.40×) |
| chair | 43, 0.50 s | 67, 0.80 s (0.63×) |
| camel | 49, 1.45 s | 48, 1.44 s |
| hand | 58, 0.46 s | 32, 0.26 s (1.8×) |

and the end-to-end `unwrap/brain` of the benchmark went from 35.1 s to 98.7 s, `unwrap/bunny`
from 3.1 s to 7.5 s. Two other schedules tried on the same four problems (`1, 10⁻², 10⁻⁴, …,
10⁻¹⁰` and `1, 10⁻¹, 10⁻³, 10⁻⁵, 10⁻⁷, 10⁻¹⁰`) were no better than the default overall (bunny 73 and
96 factorizations, chair 59 and 70, camel 61 and 51). Without bounds the exact stage ends on the
same active sets whatever Stage I did, so the schedule is only speed; with bounds the last
regularized solution is a worse start for it. Nothing here explains why, and I did not look into
the safeguard rounds that the extra factorizations are. The schedule stays the paper's, and the
short one is available for unbounded problems, where it is up to `3.9×` fewer factorizations
(`camelshotnoise` at `λ = 0.02` is the exception: 62 → 81).

## Deviations from the paper

Each is derived from the problem statement above; the proofs are in `lean/proofs`, an
independent implementation is `lean/oracle`.

1. **Sign of `u`.** The paper prints `L u = Ω − µ` (Eq. 17) with a positive semidefinite `L`,
   which makes its `u` minus the log scale factor. Here `L u = µ − Ω` and `L φ = −wE m u`; `µ`,
   `φ` and the cones are the same (`discrete-yamabe-linearization`).
2. **Normalization of `λ`.** Masses are divided by the total area. The paper says `λ` is
   normalized by the surface area without giving the formula; the one published file that is
   a plain optimum (`cow01`) corresponds to masses that are *not* normalized.
3. **Weights.** `wE` multiplies the mass in the adjoint equation, `wR` scales the threshold per
   vertex (`λ wR_i`); App. B of the paper drops the `λ`.
4. **Angle bounds** (Sec. 8.2 says only "a simple projection"): the cone equation is
   `µ = clamp(P_t(φ + µ), lo, hi)`, the clamp after the soft threshold, and the Newton method
   gets two more pieces (`cone-fixed-point-bounded`).
5. **Forbidden vertices and prescribed cones** are vertices with a fixed angle (`0`, or the
   prescribed one) and no cost. Prescribing an angle is not in the paper (its "fixed cone
   points" have free angles: `wR = 0` here).
6. **Closed meshes** are not treated in the paper's algorithm. Here the constant of `u` is a
   variable (`Σ wE m u = 0` at the optimum), `Σ µ = 2πχ` is checked for feasibility, the
   constant of `φ` (the multiplier of `Σ µ = 2πχ`) is fixed by the active vertices, and when
   no vertex is active it is chosen so that the cone angles the pieces ask for sum to `2πχ`
   (`gauss-bonnet-adjoint-shift`); the interval of admissible constants is returned as
   `adjoint_shift`.
7. **Regularized stage.** The paper's update `µ = P_λ(φ)/γ` (Eqs. 21–22) has no mass factor,
   while the `L²` norm of Eq. 12 gives `µ_i = (m_i/γ) P_λ(φ_i)`. Both are available
   (`RegularizationMetric`); the default is the paper's printed form, because it needs fewer
   Newton steps: on the 23 dataset meshes at `λ = 0.5` it took fewer factorizations on 20
   (for example `body` 33 against 52, `camel` 36 against 51, the cow 28 against 41), the same
   number on 2 and more on 1 (`chair`, 39 against 37), and never a safeguard round, where the
   mass form needed two on `camelshotnoise`. With normalized masses the printed form is about
   `n` times stiffer at the same `γ`. The choice only affects the starting point of Stage II;
   both give the same optimum. The unknown of an active row is `φ_i − (±t_i)` rather than
   `φ_i`, which keeps the digits when the stiffness is `10¹⁰`.
8. **Start of Stage II.** The last regularized `(u, φ, µ)` is used as it is; the paper's
   `µ ← P_λ(φ)` (Fig. 12) would drop the factor `1/γ`.
9. **Line search in Stage I.** A regularized problem has a convex, continuously differentiable
   dual objective; a Newton step that would increase it is cut at its minimizer along the step
   (exact, the derivative is piecewise linear). The step is the paper's whenever it leaves the
   pieces unchanged or still descends.
10. **Safeguards in Stage II.** Plain semismooth Newton cycles on ordinary inputs. Here a
    vertex moves to a *neighbouring* piece per step (never from a positive to a negative cone
    at once), a repeated set of pieces or the iteration limit ends the run, and then a proximal
    step (a regularized problem centred at the current cone angles) is inserted before the
    exact iteration is retried, at most `max_safeguard_rounds` times. None of this decides the
    status; the certificate does.
11. **Ties.** The sloped pieces are closed intervals: a vertex exactly on a kink counts as
    active, so a free cone site is always active. The paper and `lean/proofs` (`Site.pieceReg`)
    use strict inequalities at the thresholds, the Lean oracle gives both kinks to the flat
    pieces. The classifiers differ only where the adjacent affine formulas agree; the stopping
    rule "pieces unchanged ⟹ exact solution" holds for each (`newton_stageI_exact`,
    `newton_stageII_exact` for any piecewise affine classifier). For this crate's classifier
    the required property is unit-tested, not machine-checked.
12. **Linear solver.** Sparse `LDLᵀ`/`LU` of the saddle-point matrix instead of sparse QR.
13. **Rounding-level clean-up** of the cone angles before the certificate (see
    [what comes back](#what-comes-back)).
14. **Cone extraction.** Threshold `10⁻¹²` as in the paper. Merging (Sec. 6.1) is off by
    default, and differs from the paper's one sentence in four ways. (a) Only cones **of the
    same sign** form a cluster. The paper's reason for merging is that one cone of the smooth
    optimum shows up as a tiny ring of cones; a positive cone next to a negative one is a
    feature of the surface (the tip and the base of a spike), and summing them cancels cone
    angle that the optimum needs. On `camelshotnoise` at `λ = 0.5` merging clusters of any sign
    (of at most 8 cones) raises the distortion from 3.49 to 1.1·10⁴; merging by sign leaves it
    at 3.49. (b) Clusters of more than `max_merge_cluster_size` cones (default 8) or wider than
    `max_merge_cluster_diameter` (default: no limit) are left alone: on convex parts of closed
    meshes the cones cover whole regions. (c) A prescribed cone is never moved, and a merged
    cone must respect the angle bounds of its vertex. (d) The Fréchet mean, which the paper
    does not define, is the member minimizing `Σ_k |µ_k| dist(v, c_k)²` with distances along
    edges inside the cluster (on `axes` one of four merged cones lands on a different vertex of
    its cluster than the authors').
15. **Feature-size weights** (Sec. 8.1): the printed formula `r = 1/(|Ω|/M + ε)` contradicts the
    sentence after it; `feature_size_weights` offers both readings and never applies either on
    its own.
16. **Disconnected meshes** are solved per component with the masses of the whole mesh
    ([disconnected meshes](#disconnected-meshes)); the paper does not discuss them.

## The published dataset

The authors published 23 meshes with the cones used to flatten them
(`data/optimal-cones/optimal-cones-data`). `examples/published.rs` investigates them;
`tests/dataset.rs` asserts what was established. Everything below is measured with those two;
the logs are `runs_tmp/cones-published-*.log` and `runs_tmp/cones-sweep-*.log`.

**The published angles are not the angles of the penalized optimum.** In an optimum of the
problem above the adjoint is `φ = ±λ wR` at every cone. In the published files it is zero at
every cone (on closed meshes: one constant, whatever the sign of the cone). That is the
optimality condition of a *free cone site*: the published angles are the angles that minimize
the distortion for the published positions when cones cost nothing there (in the terms of
this crate: `cone_costs = 0` at the cones, `+∞` elsewhere). The paper does not mention such a
step. Largest difference between the re-optimized and the published angle, in radians:

| file (cones) | `u = 0` on the boundary | `u` = one free constant | boundary values optimized |
|---|---|---|---|
| bump (5) | **2.4e-4** | 1.4e-1 | 1.4e-1 |
| face (2) | **4.6e-5** | 1.2e-2 | 4.9e-2 |
| mannequin (8) | **3.6e-4** | 5.0e-1 | 3.6e-1 |
| tree (8) | **2.7e-3** | 1.1e-1 | 1.1e-1 |
| body (10) | **8.6e-2** (6 of 10 within 6e-3) | 8.7e-1 | 8.8e-1 |
| camel (6) | 3.6e-1 | 1.8e-1 | **3.4e-3** |
| camel_noise1 (7) | 3.9e-1 | 1.9e-1 | **5.8e-3** |
| camel_noise2 (7) | 6.4e-1 | 8.6e-2 | **9.4e-3** |
| camelbad (6) | 3.5e-1 | 1.8e-1 | **4.6e-3** |
| camelbadnonuniform (6) | 3.6e-1 | 1.6e-1 | **9.8e-3** |
| cow1 (2) | 4.1e-1 | 9.5e-2 | **4.9e-3** |
| cow75e-2 (3) | 4.1e-1 | 3.3e-3 | **3.3e-3** |
| cow5e-1 (5) | 1.0e0 | 3.2e-1 | **1.0e-2** |
| hand (29) | 1.2e-1 | 6.1e-2 | **3.9e-3** |
| rhino (36) | 2.2e-1 | 1.7e-1 | **8.2e-3** |
| airplane (13) | 6.3e-1 | 1.2e-1 | **1.1e-1** (12 of 13 within 4.4e-3) |
| axes (4), closed | **1.3e-4** | | |
| chair (14), closed | **7.8e-4** | | |
| brain (45), closed | **2.6e-3** | | |
| bunny (4), closed | 9.0e-2 | | |

So the files fall into two families. One treats the boundary as this crate does (`u = 0`,
paper Fig. 10 b); for it the agreement is at the level that the six-digit positions of the OBJ
files allow. The other evidently used **optimized boundary values** (paper Sec. 8.3, Fig. 10 c),
which `BoundaryTreatment::Free` implements (see [the boundary](#the-boundary)); the agreement
there is `10⁻²`, not `10⁻⁴`, so some detail
of the authors' boundary treatment is not captured (uniform instead of lumped masses, and
other weights of the boundary mass, were tried and are worse or no better).

Not explained by re-optimized angles: `bunny`; `camelshotnoise` (its 6 cones are the plain
camel's, while in this crate's problem its spikes dominate everything, see
[needle-like spikes](#needle-like-spikes)); and the two
dense cow files, which are plain optima: `cow01` (439 cones) is the optimum of this crate's
problem for `λ = 0.01 / area` to `8·10⁻⁴` rad, 439 of 439 vertices, and `cow00001` (13 836
cones) is essentially `µ = Ω` (our optimum for `λ = 10⁻⁹ / area` has 14 725 cones, 13 832 of
them the file's, within `2.5·10⁻³` rad), not the optimum for `10⁻⁵ / area` (5 829 cones).

**Where the published positions come from** is only partly reproduced.

* `bump` is reproduced from start to finish: this crate's optimum for `λ = 0.01 / area` (any
  `λ·area` in `[0.0095, 0.012]`) has 32 cones in clusters; merged they are the 5 vertices of
  the file; the re-optimized angles at these vertices are the file's (total 5.1349 against
  5.1347, largest difference `2.4·10⁻⁴`). Re-optimizing before merging does not give the file
  (`6·10⁻²` and more), so the order is: penalized problem, merge, re-optimize.
* `body`: the file's 10 vertices are exactly the cones of this crate's optimum for `λ·area` in
  `[0.08, 0.09]`.
* The three sparse cow files have exactly the cones of this crate's optimum (`u = 0` on the
  boundary) for `λ·area` in `[1.39, 3.52]` (`cow1`), `[1.12, 1.39]` (`cow75e-2`) and
  `[0.97, 1.10]` (`cow5e-1`), although their angles belong to the other boundary family. With
  optimized boundary values the penalized optimum has no cone at or next to the file's vertex 9
  for any `λ·area` between 0.11 and 1.1 (and no cone at all above 0.96).
  The file names (1, 0.75, 0.5) are not these `λ` under any single factor.
* `axes` (closed, positive cones): 14 cones in 4 clusters for every `λ ≥ 0.045`; merged, 3 are
  at the file's vertices and the 4th is one edge away.
* Approximately: `face` (3 merged cones for `λ·area ≈ 0.011`, the file has 2 of them, one at
  the same vertex and one an edge away), `tree` (the 4 positive cones for `λ·area` in
  `[0.022, 0.048]`; the file's 4 negative ones never appear together, we get 5 or 9 cones),
  `brain` (55 merged cones in the positive regime, 28 of the file's 45 at the same vertex).
* Not reproduced at any `λ`, with either boundary model that could be solved: `mannequin`,
  `bunny`, `chair`, `airplane`, the camels, `hand`, `rhino`.

A certified optimum cannot be beaten in its own problem: for every mesh with boundary the
file's cones, evaluated with `evaluate_cones`, have an objective at `λ = 0.5` that is not
below ours (asserted in `tests/dataset.rs`).

**Published cones with optimized boundary values, re-run with `BoundaryTreatment::Free`**
(`tests/dataset.rs::published_cones_with_optimized_boundary_values`; largest difference between the
re-optimized and the published angle, radians). The numbers are the same as the doubled mesh gave,
to `10⁻¹¹` (the two routes were compared cone by cone; `hand` has no double): it does **not**
get closer than the `10⁻²` the doubled mesh gave. Best agreement `3.3·10⁻³`, worst `1.0·10⁻²`.

| file | `u = 0` | optimized (this mode) | doubled mesh |
|---|---|---|---|
| camel (6) | 3.6e-1 | 3.43e-3 | 3.43e-3 |
| camel_noise1 (7) | 3.9e-1 | 5.75e-3 | 5.75e-3 |
| camel_noise2 (7) | 6.4e-1 | 9.40e-3 | 9.40e-3 |
| camelbad (6) | 3.5e-1 | 4.56e-3 | 4.56e-3 |
| camelbadnonuniform (6) | 3.6e-1 | 9.80e-3 | 9.80e-3 |
| cow1 (2) | 4.1e-1 | 4.93e-3 | 4.93e-3 |
| cow75e-2 (3) | 4.1e-1 | 3.27e-3 | 3.27e-3 |
| cow5e-1 (5) | 1.0e0 | 1.02e-2 | 1.02e-2 |
| hand (29) | 1.2e-1 | 3.85e-3 | (no double) |
| rhino (36) | 2.2e-1 | 8.24e-3 | 8.24e-3 |

The residual `3·10⁻³ … 10⁻²` is thus not due to the way the boundary was optimized here; it is
still an unknown detail of the authors' computation (it was not removed by the mass variants
mentioned above).

What this means for a user: to get cones like the authors', solve, merge the clusters, and
re-optimize the angles at the merged vertices (the last row of the [steering](#steering)
table; `tests/dataset.rs::least_distortion_angles` is that recipe). The re-optimized cones are
larger than the penalized ones (they are not shrunk by `λ`) and give less distortion for the
same positions.

### Published UVs: measured comparison

`cargo run --release -p ddg-uv --example published_uvs` (one process per mesh is faster:
run it with a mesh name as argument, in parallel). Measured on 2026-10-13; logs
`runs_tmp/published-uvs/`. The measure is the dataset README's: `D = Σ_t A_t u_t²` with
`u_t = ½ log(B_t/A_t)`, `A_t` the triangle's share of the 3D area, `B_t` its UV area; `D*` is `D`
with the best global scale of the UV plane (the area-weighted mean of `u` removed). Flipped and
degenerate UV triangles are counted and left out.

What it shows:

* **The flattening is reproduced on 11 of the 23 meshes** (our Boundary First Flattening with
  the file's cones against the published `vt`, triangle by triangle: rms difference of `u` at
  most about `10⁻³`, correlation at least 0.9999), approximately on 5 more (rms up to `6·10⁻³`),
  and not on 7 (`airplane`, `body`, `camel_noise1`, `camelbadnonuniform`, `camelshotnoise`,
  `rhino`: rms `1.2·10⁻²` to `5·10⁻²`; `cow00001`). Rounding of the published coordinates is not
  the cause (its floor is `6·10⁻⁵`); what is, was not found.
* **The published UVs realize the file's cone angles** to `5·10⁻⁶` rad on all 23 meshes.
* **Every published flattening with a boundary has `u = 0` on the boundary**, also for the files
  whose cone angles belong to the optimized-boundary family (the camels, the sparse cows, `hand`,
  `rhino`, `airplane`): with `u = 0` our `D*` is within 0.969 to 1.006 of the published one, with
  optimized boundary values the match is clearly worse. So in those files the angles and the
  flattening do not use the same boundary treatment.
* **The dataset's claim** (least distortion among all conformal flattenings with the same or
  smaller total cone angle), tested with a cap of the published total plus 0.1 % (the file's
  angles have six digits): our penalized optimum (C1) has a lower `D*` on 20 meshes, by 1.6 % to
  83 %, is within 1 % on 2 (`cow01`, `face`), and is higher on `camelshotnoise` as it is (the
  spikes, see below; lower by 20 % with zero distortion weight at the 22 tips). It gets there with
  many small cones in clusters (254 on `camel` against 6 in the file). The authors' evident recipe
  at our positions (C2: solve, merge, re-optimize the angles) is lower on 10, within 1 % on 5 and
  higher on 8.
* **At the same number of cones as the file** our merged and re-optimized cones are mostly
  higher (`bump` +3.3 %, `body` +0.9 %, `camel` +4.4 %, `rhino` +15 %, `camelbad` +67 %), lower or
  equal on a few (`cow5e-1` −13 %, `axes` −0.2 %, `cow75e-2` +0.1 %): the published positions are
  better than ours for few cones. These rows are in the logs; the search for an equal count is
  approximate.

Not settled: the C2 search can miss better configurations (its total angle is not monotone in
`λ`); a run with a cap without the 0.1 % was stopped before it finished and is not used (on closed
meshes such a cap is infeasible where the file's total is below `4π` by rounding), so the two
verdicts within 1 % are checked with the slack only; several C2 rows fail to flatten because a
re-optimized cone angle reaches `2π`.

Table A and B: the published data, and our flattening with the file's cones against it.

| mesh | cones | Σ\|angle\| | Gauss–Bonnet of the file | D | D* | flipped + degenerate | charts |
|---|---:|---:|---|---:|---:|---:|---|
| airplane | 13 | 32.1502 | n/a (has boundary) | 1.2502e-1 | 1.2502e-1 | 2 + 0 | one chart |
| axes | 4 | 12.5663 | -9.06e-5 | 8.0896e-1 | 8.0892e-1 | 0 + 0 | one chart |
| body | 10 | 27.9356 | n/a (has boundary) | 3.9685e-1 | 3.9323e-1 | 1 + 56 | one chart |
| brain | 45 | 32.0843 | -2.71e-3 | 3.2931e-2 | 3.2843e-2 | 35 + 2 | one chart |
| bump | 5 | 5.1347 | n/a (has boundary) | 1.9776e-2 | 1.9693e-2 | 0 + 0 | one chart |
| bunny | 4 | 12.5656 | -7.61e-4 | 2.3244e-1 | 2.3229e-1 | 0 + 0 | one chart |
| camel | 6 | 20.3953 | n/a (has boundary) | 3.7399e-2 | 3.7399e-2 | 0 + 0 | one chart |
| camel_noise1 | 7 | 21.0340 | n/a (has boundary) | 5.9617e-2 | 5.0959e-2 | 19 + 0 | one chart |
| camel_noise2 | 7 | 13.6965 | n/a (has boundary) | 3.0651e0 | 1.7762e-1 | 113 + 0 | one chart |
| camelbad | 6 | 20.2976 | n/a (has boundary) | 4.4093e-2 | 3.6806e-2 | 1 + 0 | one chart |
| camelbadnonuniform | 6 | 19.9219 | n/a (has boundary) | 3.5926e-2 | 3.5926e-2 | 4 + 0 | one chart |
| camelshotnoise | 6 | 20.0649 | n/a (has boundary) | 6.8172e-2 | 6.8138e-2 | 70 + 0 | one chart |
| chair | 14 | 20.4193 | +7.83e-4 | 1.2758e-2 | 1.2752e-2 | 0 + 0 | one chart |
| cow00001 | 13836 | 123.9842 | n/a (has boundary) | 1.0883e-5 | 1.9331e-6 | 0 + 0 | one chart |
| cow01 | 439 | 39.3188 | n/a (has boundary) | 3.0575e-2 | 3.0076e-2 | 0 + 0 | one chart |
| cow1 | 2 | 8.0929 | n/a (has boundary) | 9.0993e-1 | 9.0993e-1 | 0 + 0 | one chart |
| cow5e-1 | 5 | 14.3826 | n/a (has boundary) | 4.2601e-1 | 4.2601e-1 | 0 + 0 | one chart |
| cow75e-2 | 3 | 10.4061 | n/a (has boundary) | 7.2186e-1 | 7.2038e-1 | 0 + 0 | one chart |
| face | 2 | 1.4673 | n/a (has boundary) | 3.8059e-3 | 3.8059e-3 | 0 + 0 | one chart |
| hand | 29 | 52.5442 | n/a (has boundary) | 1.1477e-2 | 1.0504e-2 | 0 + 0 | one chart |
| mannequin | 8 | 9.6786 | n/a (has boundary) | 2.4814e0 | 1.7559e-2 | 2 + 0 | one chart |
| rhino | 36 | 65.4459 | n/a (has boundary) | 3.7754e-2 | 3.7539e-2 | 0 + 0 | one chart |
| tree | 8 | 19.0045 | n/a (has boundary) | 4.3301e-2 | 4.3161e-2 | 0 + 0 | one chart |

| mesh | boundary that matches | D* ours / D* published | rms Δu | max \|Δu\| | corr | flipped ours / published | rms Δu, other boundary | file angle vs angle implied by the UVs | reproduces |
|---|---|---:|---:|---:|---:|---|---:|---:|---|
| airplane | u=0 | 1.2432e-1 / 1.2502e-1 = 0.9944 | 3.38e-2 | 1.07e0 | 0.99540 | 4 / 2 | 4.02e-2 | 4.63e-6 | no |
| axes | closed | 8.0895e-1 / 8.0892e-1 = 1.0000 | 2.52e-4 | 3.46e-3 | 1.00000 | 0 / 0 | - | 4.94e-6 | yes (rms ≈1e-3 or less) |
| body | u=0 | 3.9399e-1 / 3.9323e-1 = 1.0019 | 5.12e-2 | 1.06e0 | 0.99666 | 0 / 1 | 5.75e-2 | 4.72e-6 | no |
| brain | closed | 3.2837e-2 / 3.2843e-2 = 0.9998 | 1.46e-3 | 3.82e-1 | 0.99997 | 38 / 35 | - | 4.78e-6 | approximately (rms ≤1e-2) |
| bump | u=0 | 1.9693e-2 / 1.9693e-2 = 1.0000 | 3.58e-5 | 2.01e-4 | 1.00000 | 0 / 0 | 5.30e-3 | 4.47e-6 | yes (rms ≈1e-3 or less) |
| bunny | closed | 2.3230e-1 / 2.3229e-1 = 1.0001 | 4.88e-4 | 1.01e-2 | 1.00000 | 0 / 0 | - | 3.95e-6 | yes (rms ≈1e-3 or less) |
| camel | u=0 | 3.7311e-2 / 3.7399e-2 = 0.9977 | 5.57e-3 | 6.77e-1 | 0.99959 | 0 / 0 | 8.89e-2 | 4.90e-6 | approximately (rms ≤1e-2) |
| camel_noise1 | u=0 | 5.0891e-2 / 5.0959e-2 = 0.9987 | 3.47e-2 | 8.81e-1 | 0.98820 | 19 / 19 | 8.97e-2 | 4.26e-6 | no |
| camel_noise2 | u=0 | 1.7723e-1 / 1.7762e-1 = 0.9978 | 4.34e-3 | 1.67e-1 | 0.99995 | 113 / 113 | 6.31e-2 | 1.70e-6 | approximately (rms ≤1e-2) |
| camelbad | u=0 | 3.7014e-2 / 3.6806e-2 = 1.0056 | 6.07e-3 | 6.29e-1 | 0.99950 | 1 / 1 | 8.82e-2 | 4.87e-6 | approximately (rms ≤1e-2) |
| camelbadnonuniform | u=0 | 3.5782e-2 / 3.5926e-2 = 0.9960 | 1.18e-2 | 1.74e0 | 0.99806 | 1 / 4 | 8.73e-2 | 4.37e-6 | no |
| camelshotnoise | u=0 | 6.6023e-2 / 6.8138e-2 = 0.9690 | 3.55e-2 | 7.52e-1 | 0.99045 | 72 / 70 | 2.31e0 | 4.49e-6 | no |
| chair | closed | 1.2718e-2 / 1.2752e-2 = 0.9974 | 9.86e-4 | 2.30e-2 | 0.99996 | 0 / 0 | - | 4.18e-6 | yes (rms ≈1e-3 or less) |
| cow00001 | optimized | 1.0436e-6 / 1.9331e-6 = 0.5398 | 6.79e-4 | 3.30e-3 | 0.88571 | 0 / 0 | 7.32e-4 | 5.00e-7 | no |
| cow01 | u=0 | 3.0084e-2 / 3.0076e-2 = 1.0002 | 2.42e-4 | 1.20e-2 | 1.00000 | 0 / 0 | 9.37e-3 | 4.99e-7 | yes (rms ≈1e-3 or less) |
| cow1 | u=0 | 9.0999e-1 / 9.0993e-1 = 1.0001 | 1.05e-3 | 2.73e-2 | 1.00000 | 0 / 0 | 1.05e-1 | 3.84e-6 | yes (rms ≈1e-3 or less) |
| cow5e-1 | u=0 | 4.2604e-1 / 4.2601e-1 = 1.0001 | 3.44e-4 | 1.53e-2 | 1.00000 | 0 / 0 | 2.38e-1 | 4.52e-6 | yes (rms ≈1e-3 or less) |
| cow75e-2 | u=0 | 7.2040e-1 / 7.2038e-1 = 1.0000 | 1.04e-3 | 2.55e-2 | 1.00000 | 0 / 0 | 5.50e-2 | 2.60e-6 | yes (rms ≈1e-3 or less) |
| face | u=0 | 3.8056e-3 / 3.8059e-3 = 0.9999 | 1.10e-4 | 6.95e-3 | 1.00000 | 0 / 0 | 1.47e-2 | 1.50e-6 | yes (rms ≈1e-3 or less) |
| hand | u=0 | 1.0534e-2 / 1.0504e-2 = 1.0028 | 4.87e-3 | 4.62e-2 | 0.99887 | 0 / 0 | 1.34e-2 | 4.87e-6 | approximately (rms ≤1e-2) |
| mannequin | u=0 | 1.7538e-2 / 1.7559e-2 = 0.9988 | 7.42e-4 | 1.60e-2 | 0.99998 | 2 / 2 | 1.83e-2 | 4.50e-6 | yes (rms ≈1e-3 or less) |
| rhino | u=0 | 3.7122e-2 / 3.7539e-2 = 0.9889 | 1.50e-2 | 1.59e-1 | 0.99700 | 0 / 0 | 3.00e-2 | 4.96e-6 | no |
| tree | u=0 | 4.3162e-2 / 4.3161e-2 = 1.0000 | 2.61e-4 | 5.04e-2 | 1.00000 | 0 / 0 | 4.32e-3 | 4.75e-6 | yes (rms ≈1e-3 or less) |

Table C: distortion at equal or smaller total cone angle (first block of each mesh: the boundary
treatment of the file's family).

| mesh (boundary of the family) | published: cones, Σ\|µ\|, solver ½Σmu², D* | C1: cones, Σ\|µ\|, λ, solver, D* | C1 vs published D* | C2: cones, Σ\|µ\|, λ, solver, D* | C2 vs published D* |
|---|---|---|---:|---|---:|
| airplane | 13, 32.1502, 0.06263, 1.2502e-1 | 116, 32.1823, 6.5017e-3, 0.04484, 8.8078e-2 | -29.5% | 10, 31.9137, 4.5443e-2, 0.08268, 1.6560e-1 | +32.5% |
| axes | 4, 12.5663, 0.39905, 8.0892e-1 | 15, 12.5786, 3.4709e-2, 0.39823, 7.9630e-1 | -1.6% | 4, 12.5664, 1.0000e6, 0.39901, 8.0698e-1 | -0.2% |
| body | 10, 27.9356, 0.20317, 3.9323e-1 | 409, 27.9635, 9.6354e-3, 0.16586, 3.2797e-1 | -16.6% | 10, 27.9207, 4.6953e-2, 0.19864, 3.9239e-1 | -0.2% |
| brain | 45, 32.0843, 0.01715, 3.2843e-2 | 287, 32.1160, 1.8161e-4, 0.01199, 2.2696e-2 | -30.9% | 65, 26.5667, 1.1751e-3, 0.01625, 3.1045e-2 | -5.5% |
| bump | 5, 5.1347, 0.00997, 1.9693e-2 | 40, 5.1397, 1.7724e-3, 0.00959, 1.8260e-2 | -7.3% | 5, 5.1349, 2.1288e-3, 0.00997, 1.9693e-2 | +0.0% (=) |
| bunny | 4, 12.5656, 0.11512, 2.3229e-1 | 254, 12.5780, 1.3555e-2, 0.11045, 2.2036e-1 | -5.1% | 15, 12.5664, 1.0000e6, 0.11102, 2.2319e-1 | -3.9% |
| camel | 6, 20.3953, 0.01491, 3.7399e-2 | 254, 20.4156, 1.0728e-3, 0.00858, 1.7116e-2 | -54.2% | 9, 17.3144, 3.3003e-3, 0.01681, 3.3497e-2 | -10.4% |
| camel_noise1 | 7, 21.0340, 0.02157, 5.0959e-2 | 115, 21.0548, 1.2310e-3, 0.01475, 2.8397e-2 | -44.3% | 22, 18.0495, 3.0080e-3, 0.02287, 4.4215e-2 | -13.2% |
| camel_noise2 | 7, 13.6965, 0.09698, 1.7762e-1 | 35, 13.7101, 3.3666e-3, 0.09299, 1.6565e-1 | -6.7% | 15, 13.5681, 2.4144e-2, 0.09563, 1.7269e-1 | -2.8% |
| camelbad | 6, 20.2976, 0.01457, 3.6806e-2 | 210, 20.3179, 1.0788e-3, 0.00855, 1.7411e-2 | -52.7% | 12, 17.9253, 3.1908e-3, 0.01658, 3.3467e-2 | -9.1% |
| camelbadnonuniform | 6, 19.9219, 0.01497, 3.5926e-2 | 140, 19.9418, 1.2510e-3, 0.00929, 1.7600e-2 | -51.0% | 9, 18.2938, 3.1623e-3, 0.01598, 3.1694e-2 | -11.8% |
| camelshotnoise | 6, 20.0649, 84784.65243, 6.8138e-2 | 22, 20.0834, 1.0432e3, 61463.19595, 7.3762e-1 | +982.5% | 3, 18.8637, 1.3917e3, 71275.77592, no flattening | failed (cone angle ≥ 2π) |
| camelshotnoise[wE=0 at tips] | 6, 20.0649, 0.13550, 6.8138e-2 | 151, 20.0848, 1.6352e-3, 0.01840, 5.4351e-2 | -20.2% | 10, 19.8262, 6.0566e-3, 0.02426, 6.0572e-2 | -11.1% |
| chair | 14, 20.4193, 0.00657, 1.2752e-2 | 243, 20.4394, 4.5409e-4, 0.00305, 5.9280e-3 | -53.5% | 10, 20.3272, 5.5232e-4, 0.00456, 8.3847e-3 | -34.2% |
| cow00001 (u=0) | 13836, 123.9842, 0.00000, 1.9331e-6 | 14236, 122.1582, 1.0000e-9, 0.00000, 8.1615e-10 | -100.0% | 337, 124.0967, 1.0698e-9, 0.00183, 5.6212e-3 | +290685.3% |
| cow01 (u=0) | 439, 39.3188, 0.01540, 3.0076e-2 | 441, 39.3579, 2.2349e-3, 0.01531, 2.9913e-2 | -0.5% | 8, 29.3834, 1.1262e-2, 0.08028, 1.6008e-1 | +432.2% |
| cow1 | 2, 8.0929, 0.44975, 9.0993e-1 | 5, 8.1010, 5.6065e-2, 0.35996, 7.2054e-1 | -20.8% | 1, 4.3525, 2.1500e-1, 0.85271, 1.6999e0 | +86.8% |
| cow5e-1 | 5, 14.3826, 0.18475, 4.2601e-1 | 18, 14.3967, 1.4802e-2, 0.16941, 3.3992e-1 | -20.2% | 5, 14.3936, 3.1623e-2, 0.18360, 3.6965e-1 | -13.2% |
| cow75e-2 | 3, 10.4061, 0.35835, 7.2038e-1 | 7, 10.4165, 3.5420e-2, 0.25405, 5.0910e-1 | -29.3% | 3, 10.1345, 8.7723e-2, 0.36031, 7.2115e-1 | +0.1% |
| face | 2, 1.4673, 0.00220, 3.8059e-3 | 16, 1.4687, 6.4352e-4, 0.00214, 3.8153e-3 | +0.2% | 3, 1.4667, 3.1623e-3, 0.00218, 3.7996e-3 | -0.2% |
| hand | 29, 52.5442, 0.00642, 1.0504e-2 | 536, 52.5967, 2.2657e-4, 0.00102, 1.7605e-3 | -83.2% | 22, 52.5002, 1.8455e-4, 0.01606, 7.6183e-2 | +625.2% |
| mannequin | 8, 9.6786, 0.00916, 1.7559e-2 | 673, 9.6882, 1.4974e-3, 0.00596, 1.1463e-2 | -34.7% | 6, 9.6540, 2.8838e-3, 0.01706, 3.1095e-2 | +77.1% |
| rhino | 36, 65.4459, 0.01858, 3.7539e-2 | 735, 65.5109, 2.4770e-4, 0.00512, 1.0011e-2 | -73.3% | 49, 65.1052, 6.3476e-4, 0.00977, 1.8595e-2 | -50.5% |
| tree | 8, 19.0045, 0.02178, 4.3161e-2 | 482, 19.0233, 1.9094e-3, 0.01146, 2.2589e-2 | -47.7% | 5, 16.5022, 1.0000e-4, 0.02786, 5.3293e-2 | +23.5% |

| mesh | other boundary: C1 vs published | C2 vs published |
|---|---:|---:|
| airplane | -15.6% | +54.2% |
| body | -24.0% | -7.0% |
| bump | -22.1% | +2.5% |
| camel | -50.2% | -0.6% |
| camel_noise1 | -39.4% | +9.1% |
| camel_noise2 | +1.8% | +12.2% |
| camelbad | -48.5% | -18.6% |
| camelbadnonuniform | -46.8% | +55.6% |
| camelshotnoise | +57943.7% | +7334.9% |
| camelshotnoise[wE=0 at tips] | -27.1% | -7.8% |
| cow00001 (u=0) | -100.0% | +289097.3% |
| cow01 (u=0) | -17.1% | +172.2% |
| cow1 | -20.1% | failed (cone angle ≥ 2π) |
| cow5e-1 | -9.8% | +20.6% |
| cow75e-2 | -27.3% | -14.2% |
| face | -9.6% | -6.4% |
| hand | -82.0% | +424.9% |
| mannequin | -49.7% | +56.4% |
| rhino | -72.1% | -46.5% |
| tree | -53.3% | -13.5% |

| mesh | C1 at n = published | C2 at n = published |
|---|---|---|
| airplane | n=13, Σ=14.7044, D*=9.6996e-1: +675.8% | n=13, Σ=50.0242, D*=2.5320e-1: +102.5% [Σ>cap] |
| axes | n=14, Σ=12.5664, D*=7.9753e-1: -1.4% (n differs) | n=4, Σ=12.5664, D*=8.0698e-1: -0.2% |
| body | n=10, Σ=12.0846, D*=1.2339e0: +213.8% | n=10, Σ=27.8317, D*=3.9688e-1: +0.9% |
| brain | n=123, Σ=12.5664, D*=4.6882e-2: +42.8% (n differs) | n=55, Σ=12.5664, D*=4.6886e-2: +42.8% (n differs) |
| bump | n=5, Σ=2.6444, D*=2.5621e-2: +30.1% | n=5, Σ=4.9740, D*=2.0338e-2: +3.3% |
| bunny | n=229, Σ=15.3832, D*=1.5551e-1: -33.1% (n differs) [Σ>cap] | n=15, Σ=12.5664, D*=2.2319e-1: -3.9% (n differs) |
| camel | n=6, Σ=4.1317, D*=6.1295e-1: +1538.9% | n=6, Σ=15.8866, D*=3.9047e-2: +4.4% |
| camel_noise1 | n=7, Σ=3.9127, D*=2.0003e-1: +292.5% | n=7, Σ=12.0367, D*=7.5908e-2: +49.0% |
| camel_noise2 | n=7, Σ=2.4138, D*=6.8552e-1: +285.9% | n=7, Σ=12.5917, D*=3.9050e-1: +119.9% |
| camelbad | n=6, Σ=4.0910, D*=7.0258e-1: +1808.9% | n=6, Σ=11.1575, D*=6.1524e-2: +67.2% |
| camelbadnonuniform | n=6, Σ=4.1362, D*=6.3507e-1: +1667.7% | n=6, Σ=11.0813, D*=5.8888e-2: +63.9% |
| camelshotnoise | n=6, Σ=2.7520, D*=1.4323e0: +2002.1% | n=6, Σ=37.6815, D*=1.4052e0: +1962.3% [Σ>cap] |
| camelshotnoise[wE=0 at tips] | n=6, Σ=3.7052, D*=7.4853e-1: +998.5% | n=6, Σ=15.3416, D*=9.5257e-2: +39.8% |
| chair | n=73, Σ=12.5664, D*=2.3976e-2: +88.0% (n differs) | n=13, Σ=26.9170, D*=4.9870e-3: -60.9% (n differs) [Σ>cap] |
| cow00001 (u=0) | n=13837, Σ=120.5789, D*=2.7875e-9: -99.9% (n differs) | n=341, Σ=124.3630, D*=5.5931e-3: +289232.9% (n differs) [Σ>cap] |
| cow01 (u=0) | n=439, Σ=39.3030, D*=3.0153e-2: +0.2% | n=341, Σ=125.1340, D*=5.6236e-3: -81.3% (n differs) [Σ>cap] |
| cow1 | n=2, Σ=4.3338, D*=1.2697e0: +39.5% | n=2, Σ=8.1016, D*=8.9933e-1: -1.2% [Σ>cap] |
| cow5e-1 | n=5, Σ=8.0821, D*=7.2266e-1: +69.6% | n=5, Σ=14.3936, D*=3.6965e-1: -13.2% |
| cow75e-2 | n=3, Σ=4.8028, D*=1.1974e0: +66.2% | n=3, Σ=10.1345, D*=7.2115e-1: +0.1% |
| face | n=2, Σ=0.4125, D*=1.1597e-2: +204.7% | n=2, Σ=1.3444, D*=4.2932e-3: +12.8% |
| hand | n=29, Σ=21.7938, D*=5.6625e-1: +5290.5% | n=29, Σ=55.2651, D*=8.2852e-2: +688.7% [Σ>cap] |
| mannequin | n=8, Σ=3.9050, D*=2.5909e-1: +1375.5% | n=8, Σ=13.3909, D*=2.7356e-2: +55.8% [Σ>cap] |
| rhino | n=36, Σ=15.4330, D*=1.5627e-1: +316.3% | n=36, Σ=56.3305, D*=4.3284e-2: +15.3% |
| tree | n=8, Σ=6.8537, D*=1.5862e-1: +267.5% | n=7, Σ=21.2411, D*=4.2963e-2: -0.5% (n differs) [Σ>cap] |

### Published UVs with the feedback loop

C1 of Table C with [the feedback loop](#the-feedback-loop) at every probe of the search
(`DDG_RECIPES=c1,c1fb cargo run --release -p ddg-uv --example published_uvs`), against C1 without
it and, on `camelshotnoise`, against zero weight at the 22 tips set by hand. Not measured yet
(phase 2 of the work).

## Choosing λ

`λ` here is for masses divided by the area. For masses as they are in a file,
`λ_raw = λ · area`. No single rule maps the published examples to a `λ`: the values that
reproduce the published positions, where any do, span three orders of magnitude.

| file | published cones (total angle) | `λ` here that gives the published positions | `λ · area` | at `λ = 0.5`: cones (merged) |
|---|---|---|---|---|
| bump | 5 (5.13) | 0.0019 – 0.0024 | 0.0095 – 0.012 | 0 |
| face | 2 (1.47) | about 0.003 (approximately) | about 0.011 | 0 |
| tree | 8 (19.0) | 0.0036 – 0.0079 (the 4 positive ones only) | 0.022 – 0.048 | 4 (4) |
| body | 10 (27.9) | 0.050 – 0.056 | 0.08 – 0.09 | 4 (4) |
| cow5e-1 | 5 (14.4) | 0.217 – 0.247 | 0.97 – 1.10 | 2 (2) |
| cow75e-2 | 3 (10.4) | 0.251 – 0.311 | 1.12 – 1.39 | 2 (2) |
| cow1 | 2 (8.09) | 0.311 – 0.790 | 1.39 – 3.52 | 2 (2) |
| cow01 | 439 (39.3) | 0.002244 | 0.0100 | 2 (2) |
| cow00001 | 13 836 (124) | about 2·10⁻¹⁰ or less (`µ ≈ Ω`) | about 10⁻⁹ or less | 2 (2) |
| axes (closed) | 4 (4π) | any ≥ 0.045 | ≥ 0.14 | 14 (4) |
| brain (closed) | 45 (32.1) | partly, any ≥ 0.002 | | 123 (55) |
| bunny (closed) | 4 (4π) | none | | 254 (15) |
| chair (closed) | 14 (20.4) | none | | 73 (7) |
| mannequin | 8 (9.68) | none | | 4 (1) |
| airplane | 13 (32.2) | none | | 6 (4) |
| camel, camel_noise1/2, camelbad, camelbadnonuniform | 6 – 7 (13.7 – 21.0) | none | | 1 – 8 (1 – 6) |
| camelshotnoise | 6 (20.1) | none | | 452 (51) |
| hand | 29 (52.5) | none | | 6 (3) |
| rhino | 36 (65.4) | none | | 2 (1) |

Remarks for choosing a default:

* `λ = 0.5`, the paper's stated default, gives few cones here: 0 to 8 on the meshes with
  boundary (452 on `camelshotnoise`, at and around its 22 spikes). On nearly developable disks (`bump`, `face`) it gives none, and cones appear only
  below `λ ≈ 0.1` and `0.03`.
* On a closed genus-0 mesh the optimum does not depend on `λ` above a mesh-dependent value
  (`axes` 0.045, `bunny` 0.016, `chair` 0.003, `brain` 0.002): all cones are then positive, the
  cone cost is `4πλ` whatever the configuration, and the problem is the search for the best
  positive measure. That optimum is not sparse: it comes in clusters (`bunny` 254 cones in 12
  clusters, `brain` 123, `chair` 73) and needs merging, with a cluster limit above the default
  8 (`bunny`: 228 cones are left with the limit, 15 without). See
  [closed meshes](#closed-meshes).
* The distortion-minimizing angles for a set of positions do not depend on `λ` at all.

## Closed meshes

On a closed mesh `Σ µ = 2πχ`. On a genus-0 mesh whose optimum has only positive cones (above
a mesh-dependent `λ`, or with `µ ≥ 0` imposed) the cone cost is the constant `4πλ`: **`λ` steers
nothing**, and the optimum is the best positive measure, which is not sparse.
`examples/steering.rs` measures what the options of this crate give on an icosphere, an
ellipsoid and a capped cylinder (642, 642 and 314 vertices; `closed`), and on the dataset's
closed meshes (`closed-dataset`). Logs: `runs_tmp/cones-steering-*.log`.

| setting | icosphere | ellipsoid | capped cylinder |
|---|---|---|---|
| plain, or `µ ≥ 0`, or `µ ≤ π` | a cone at every vertex, `µ = Ω`, any `λ`; merged: one cone of `4π` | the same | 48 cones on the two rims; merged: two cones of `2π` |
| one free site (with or without `µ ≥ 0`), `λ ≥ 0.3` | the free site alone, with `4π` | the same | the same (`λ = 1`) |
| the same, `λ = 0.1` | 8 cones; merged: `3.50π` (the free site) and `0.50π` | 24 cones; merged `2.69π`, `1.31π` | 26 cones; merged `2.57π`, `1.43π` |
| the same, `λ = 0.03` | 143 cones; merged `2.10π`, `1.90π` | 166; merged `2.02π`, `1.98π` | 27; merged 3, largest `2.08π` |
| the same, `λ = 0.01` | 316 cones; merged `1.27π`, `2.73π` | 334; merged `1.52π`, `2.48π` | 39; merged `0.72π`, `1.32π`, `1.97π` |
| one free site and `µ ≤ π` everywhere | 388 cones, any `λ ≥ 0.01` | 496 | 45 |
| 8 free sites spread over the surface, cones forbidden elsewhere, `0 ≤ µ ≤ π` | 8 cones, `0.24π … 0.59π`, distortion 0.0084 | 8 cones, `0.27π … π`, 0.0091 | 7 cones, `0.01π … 0.93π`, 0.0155 |

What this says:

* A single free cone site (paper Sec. 8.1) acts as a puncture: it holds more than `2π`, which
  is not a cone, unless `λ` is small, and then the rest of the curvature sits in one cluster
  opposite to it whose sum again exceeds `2π` or nearly does. Two cones on a sphere must have
  `2π` each. Only the cylinder at `λ ≤ 0.01` gave a configuration that exists (three cones).
  Bounding the free site as well just spreads the remainder over hundreds of vertices.
* A sparse configuration that can be flattened needs at least three cones below `2π`. On a
  round shape nothing in the geometry singles out positions, so the caller has to: several
  free sites, cones forbidden elsewhere, and bounds `0 ≤ µ ≤ hi` with `hi < 2π` and
  `hi · (number of sites) ≥ 4π`. The solver then returns the least-distortion angles for these
  sites. `λ` is not involved.
* On closed meshes that are not round the positive optimum comes in a few clusters, and
  merging them is enough. With `µ ≥ 0`, merging without a size limit:

| mesh | cones of the optimum | clusters | after merging | largest angle | distortion before → after | angles re-optimized at the merged positions, `0 ≤ µ ≤ 1.9π` |
|---|---:|---:|---:|---:|---|---|
| axes | 14 | 4 | 4 | `1.02π` | 0.39866 → 0.39901 | 4 cones, `0.99π … 1.02π`, 0.39901 |
| bunny | 254 | 12 | 15 | `1.09π` | 0.11061 → 0.11119 | 15 cones, up to `1.09π`, 0.11102 |
| chair | 73 | 7 | 7 | `0.87π` | 0.012171 → 0.012435 | 7 cones, `0.18π … 0.87π`, 0.012431 |
| brain | 123 | 31 | 55 | `0.31π` | 0.024134 → 0.024136 | 55 cones, up to `0.31π`, 0.024135 |

  The number of cones is then set by the geometry, not by a parameter; fewer cones means
  choosing a subset of the merged positions and re-optimizing the angles there.

## Needle-like spikes

`camelshotnoise` is the camel with 22 needles: slant height 0.52, base radius `8·10⁻⁵`
(`examples/published.rs spikes`). The paper (Sec. 7.2, Fig. 20) says such outliers are ignored
because they carry almost no area. In the linearized problem they are not negligible. For a
vertex lifted by `h` over a base of radius `ρ` the cotan weights at the tip are of order `ρ/h`,
so the cone-free scale factor there is of order `−h/ρ`, on a mass of order `ρh`: a distortion of
order `h³/ρ`, growing as the needle gets thinner. (The smooth problem agrees: a cone of tip
angle `θ` flattens with `u = (2π/θ − 1) log r`.) Measured on a flat unit grid with one vertex
lifted:

| `h/ρ` | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|---:|
| `E · ρ / h³` (`ρ = 0.05`) | 0.40 | 0.75 | 1.13 | 1.40 | 1.49 | 1.42 |
| `u_tip · ρ / h` | −0.71 | −1.10 | −1.41 | −1.61 | −1.72 | −1.78 |

and at fixed `h = 0.4`, halving `ρ` from 0.1 to 0.0125 multiplies `E` by 2.7, 2.3 and 2.2. On
`camelshotnoise` (`h/ρ ≈ 6200`) the cone-free scale factor at the tips is `−1.3·10⁴` to
`−1.7·10⁴`, and the 22 tips carry `8.4822·10⁴` of the `8.4826·10⁴` cone-free distortion. What
the solver does with that, at `λ = 0.5` (0.1 and 0.02 behave alike):

| | cones | negative | angles | distortion |
|---|---:|---:|---|---:|
| the optimum | 452 | 430 | −6.84 … 6.28 | 3.49 |
| merged by sign, clusters ≤ 8 (this crate's default limit) | 417 | 395 | −6.84 … 6.28 | 3.49 |
| merged by sign, all clusters | 51 | 29 | −34.8 … 6.28 | 3.85 |
| merged regardless of sign, clusters ≤ 8 (not what this crate does) | 420 | 399 | −4.28 … 6.28 | 1.1·10⁴ |
| merged regardless of sign, all clusters | 43 | 24 | −28.5 … 6.28 | 1.3·10⁵ |
| `µ ≥ 0` | 22 | 0 | 6.19 … 6.27 | 493 |
| `wE = 0` at the 22 tips | 2 | 0 | 0.31, 2.79 | 1.37 (weighted) |

The certified optimum is to put a cone of almost `2π` on every tip and negative cones around
its base; cancelling the two by merging restores the distortion of the bare spikes, which is
why this crate never merges cones of opposite sign. With positive cones only the base cannot
be compensated and the distortion stays at 493 whatever `λ` is. To have the solver ignore the
spikes, as the paper describes, exclude their tips from the distortion measure
(`distortion_weights = 0` there): the result is then that of the plain camel (2 cones at
`λ = 0.5`, 9 at `0.1`, 10 at `0.02`, all positive).

The published file has 6 cones at the positions of the plain camel's, with nearly its angles.
How the authors arrived there on this mesh was not found: re-optimizing the angles at the 6
positions gives nonsense (the spikes dominate), with the tips and their rings as additional
free sites the angles are off by 0.1 to 1.0, and with zero weight at the tips by 0.35.

**The hypothesis that the authors' weights made the tips weightless, tested with optimized boundary
values.** The first of those experiments used a constant boundary value, under which even the plain
`camel` is off by 0.36. Redone with both boundary models (`example boundary -- spikes`, log
`runs_tmp/boundary-spikes*.log`): the angles at the 6 published positions re-optimized (free
sites there, no cones elsewhere, so independent of `λ` and of the overall scale of the weights),
largest difference from the file:

| distortion weights | `camel` `u = 0` | `camel` optimized | `camelshotnoise` `u = 0` | `camelshotnoise` optimized |
|---|---:|---:|---:|---:|
| uniform | 0.36 | **3.4e-3** | 39 | 44 |
| zero at the 22 tips | | | 1.04 | **0.55** |
| zero at the tips and their 1-rings | | | 1.04 | **0.55** |
| feature size as printed, `ε` = 0.01, 1, 100, 10⁴ × median `\|Ω\|/M`, and at the maximum | 0.37 – 0.45 | 5.7e-3 – 0.30 | 1.3 – 31 | 0.87 – 34 |
| feature size as described, same `ε` | 0.26 – 0.98 | 1.4e-3 – 1.4 | 49 – 1.1e4 (not converged, so unreliable) | 57 – 6.8e3 (not converged) |

(`\|Ω\|/M` on `camelshotnoise`: median 56, 99 % 2.8·10⁴, max 2.5·10⁵.) On the control `camel`,
feature-size weights with `ε` of the order of the maximum make the weights nearly uniform and
give about the uniform value (the best, 1.4·10⁻³ at `ε` = 5·10⁵, is as good as uniform, 3.4·10⁻³);
smaller `ε` is worse. On `camelshotnoise` nothing matches: the best is zero weight at the tips, 0.55 with
optimized boundary values (and rings make no difference), against 3.4·10⁻³ for the plain camel.
Penalized solves with zero weight at the tips (and with the rings, the same), `λ` from 2 to 0.001:
with `u = 0`: no cone for `λ ≥ 1`; for `λ ≤ 0.5`, 1 to 3 of the 6 published vertices
and 1 to 245 other cones, one cone next to a tip (at `λ = 0.002`); with optimized boundary values: no
cone for `λ ≥ 0.2`, for `λ ≤ 0.1` 0 to 3 of the 6 and 2 to 231 other cones, none at a tip. The
as-described feature-size weights were not solved penalized, since their refit does not converge. So this
hypothesis is **not** confirmed: the weights do not reproduce the published angles on this mesh.

## The feedback loop

The canonical problem measures distortion with the solution `u` of the **linearized** Yamabe
equation. Where that linearization is poor the solver optimizes a distortion that no flattening
has: at the tips of needles ([needle-like spikes](#needle-like-spikes)) the linearized `u` is of
order `−h/ρ` (`−1.5·10⁴` on `camelshotnoise`) while a flattening places the tip between its
neighbours (about `−4`), and the solver spends its cones there (C1 at `+982 %` `D*`, [published
UVs](#published-uvs-measured-comparison)). Zero distortion weight at the 22 tips, set by hand,
gives `−20 %`. The feedback loop finds such vertices from a flattening instead of by hand
(ROADMAP, "Scheduled" 1 of 2026-10-14).

**The rule** (`src/feedback.rs`).

1. Solve with the current weights `wE`.
2. Flatten the solution's `cones` with its `log_scale` (the caller does this; this crate does not
   depend on a flattener) and measure the realized log area factor of every face,
   `r_t = ½ ln(|B_t|/A_t)` (`A_t` on the surface, `B_t` in the plane; absolute value, so a flipped
   face still says how much it was scaled; `ddg_uv::realized_face_log_scale`).
3. Compare at every vertex (`scale_mismatch`):
   `δ_i = (u_i − R_i) − c`, `R_i = Σ_{t∋i} A_t r_t / Σ_{t∋i} A_t`, `c` the mass-weighted
   median of `u − R` over the measured vertices of the connected component.
4. Every interior vertex with `wE_i > 0`, not reweighted before, with `|δ_i| > threshold` gets
   `wE_i ← reduced_weight · wE_i`.
5. Repeat until a round picks nothing (`Stable`). Stop without applying a round if `max_rounds`
   solves were done (`RoundLimit`), if more than `max_reweighted_fraction` of the vertices would
   be reweighted in total (`TooManyVertices`), or if a component would be left without a positive
   weight (`NoWeightLeft`). A solve or flattening that fails after a reweighting undoes it and
   keeps the previous solution (`RoundFailed`).

**Design decisions, and why.**

| question | decision | reason |
|---|---|---|
| where it lives | the rule (`Feedback`, `scale_mismatch`) and a driver with a flattening closure (`solve_with_feedback`) here; the pipeline's driver in `ddg-uv` (`UnwrapOptions::cone_feedback`, `ddg_uv::cone_feedback`) | this crate does not depend on `ddg-bff`; `ddg-uv` depends on both and has its own solve (policies for closed pieces, merge guard, slits), so it drives the same state machine with its own solve and flattening instead of a second copy of the rule |
| what is compared | the solver's vertex value `u_i` with the area-weighted mean `R_i` of the realized face values around it | comparing the solver's own face values `(u_a + u_b + u_c)/3` with `r_t` would spread a tip's `−1.5·10⁴` over its neighbours (the ring of a needle would be flagged too); the price is a bias where `u` has a kink: at a cone `u_i` lies `(µ − Ω_i)/Σ_j w_ij` below its neighbours, about `1.7` for `µ = 0.95 · 2π` on a regular mesh |
| scale invariance | `δ` is a difference of logs (dimensionless), and the offset `c` removes the free constant of `u` (closed meshes) and the scale of the plane; `δ` is invariant under both (`lean/proofs`, `cone-feedback-mismatch`) | a closed mesh has a free constant, and a flattening may be scaled |
| offset | weighted median, not mean | the tips' `\|δ\|` of `10⁴` on `0.26 %` of the area would move a mean by about 13 and flag every other vertex |
| "far" | `\|δ_i\| > threshold`, default `4` (a factor `e⁸ ≈ 3000` in area), provisional until the measurements below | a needle tip has `\|δ\|` of order `h/ρ` (`10⁴` on `camelshotnoise`, about 50 on the synthetic spikes of the tests); a cone near `2π` about 1.7 |
| per vertex or per region | per vertex | zero weight at the tips and at the tips plus their 1-rings give the same angles ([needle-like spikes](#needle-like-spikes)); with the vertex value the ring is not flagged |
| how the weight drops | to `0` by default (`reduced_weight`, any factor in `[0, 1)`) | what the hand-set experiment did; a factor `f > 0` keeps the minimizer unique, but the tip's `m u²` of `10⁴` times `f` must be negligible, so `f` would have to be tiny anyway |
| which vertices | interior vertices with a positive weight, each at most once | a zero weight at a `Free` boundary vertex can make the problem ill posed ([the boundary](#the-boundary)); weights only decrease, so the loop ends |
| round cap | `max_rounds = 4` solves, the first included; `1` measures only | after the tips are excluded the next round should be stable; a needle whose tip got a cone in round 0 shows up in round 1 when that cone moves away |
| guard | at most `1 %` of the vertices (`max_reweighted_fraction`); a round that would exceed it is not applied | a mismatch on many vertices is not a spike; lowering weights everywhere would change the problem, not repair it. With `1 %` a mesh of fewer than 100 vertices is never reweighted (the rule of `ddg-uv`'s spike check) |
| stop | a round that picks no vertex | then the solution agrees with its flattening at the threshold |
| default | off in this crate (it needs a flattener); in the pipeline, Python and Blender decided by the measurements below | the user's rule: the more stable and reliable option is the default |

**What is proved and what is not.** `cone-feedback-mismatch` (`DdgProofs/ConeFeedback.lean`):
`faceMean` shifts with its values, the weighted median is shift-equivariant (the quantile of
`cut-weighted-quantile`, the same scan in Rust), `δ` does not change when a constant is added
to `u` and another to every `r_t`, and `δ = 0` wherever the flattening realizes `u` up to such
constants. The rule itself is a design choice (`CITED.md`, `cone-feedback-rule`): there is
nothing to prove about a threshold. Each solve of the loop is certified as any other; the result
is the certified optimum of the canonical problem **with the returned weights**
(`FeedbackReport::weights`).

**Cost.** On a mesh where nothing is reweighted: one comparison (`O(F + V log V)`, a sort for
the median) after the solve and the flattening that were made anyway. Per extra round: one solve
and one flattening. A zero weight disables the factor-once method
([not implemented](#not-implemented)), so the extra solves run the reference method. The
measured numbers are below.

**Measured.** (Not measured yet: phase 2 of the work.)

## Performance baseline

`cargo run --release -p ddg-optimal-cones --example baseline -- --merge --repeat 3`: every
dataset mesh at `λ = 0.5`, the **reference method** (`baseline` forces it) with default options, one thread, fastest of three runs. AMD Ryzen
Threadripper 3970X; wall-clock seconds inside `solve`.

| mesh | triangles | total | numeric factorization | triangular solves | rest | Newton steps I / II | factorizations | solves | cones | merged (≤ 8) | merged (all) | distortion | after merging all |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| airplane | 4 656 | 0.12 | 0.09 | 0.01 | 0.01 | 31 / 1 | 32 | 64 | 6 | 4 | 4 | 1.499 | 1.499 |
| axes | 7 466 | 0.40 | 0.34 | 0.04 | 0.02 | 31 / 1 | 32 | 64 | 14 | 4 | 4 | 0.3987 | 0.3990 |
| body | 23 892 | 1.28 | 1.10 | 0.14 | 0.04 | 32 / 1 | 33 | 66 | 4 | 4 | 4 | 1.249 | 1.249 |
| brain | 305 782 | 33.5 | 25.9 | 6.20 | 1.35 | 49 / 1 | 50 | 100 | 123 | 63 | 55 | 0.02413 | 0.02414 |
| bump | 6 528 | 0.06 | 0.04 | 0.01 | 0.01 | 11 / 1 | 12 | 24 | 0 | 0 | 0 | 0.1785 | 0.1785 |
| bunny | 28 576 | 2.00 | 1.43 | 0.49 | 0.08 | 43 / 1 | 44 | 88 | 254 | 228 | 15 | 0.1106 | 0.1112 |
| camel | 22 704 | 1.08 | 0.77 | 0.27 | 0.05 | 35 / 1 | 36 | 73 | 2 | 1 | 1 | 1.361 | 1.361 |
| camel_noise1 | 22 704 | 1.06 | 0.76 | 0.26 | 0.05 | 34 / 1 | 35 | 70 | 5 | 1 | 1 | 1.399 | 1.399 |
| camel_noise2 | 22 704 | 1.03 | 0.73 | 0.26 | 0.05 | 33 / 1 | 34 | 68 | 8 | 6 | 6 | 1.607 | 1.607 |
| camelbad | 22 725 | 0.96 | 0.60 | 0.32 | 0.03 | 34 / 1 | 35 | 105 | 2 | 2 | 2 | 1.360 | 1.360 |
| camelbadnonuniform | 18 125 | 0.62 | 0.41 | 0.18 | 0.03 | 35 / 1 | 36 | 108 | 1 | 1 | 1 | 1.348 | 1.348 |
| camelshotnoise | 22 968 | 1.82 | 1.22 | 0.52 | 0.08 | 56 / 1 | 57 | 148 | 452 | 417 | 51 | 3.493 | 3.847 |
| chair | 11 448 | 0.49 | 0.36 | 0.09 | 0.03 | 38 / 1 | 39 | 78 | 73 | 69 | 7 | 0.01217 | 0.01243 |
| cow (5 files, one mesh) | 29 999 | 0.91 | 0.62 | 0.25 | 0.04 | 27 / 1 | 28 | 56 | 2 | 2 | 2 | 1.890 | 1.890 |
| face | 27 150 | 0.42 | 0.29 | 0.11 | 0.02 | 11 / 1 | 12 | 24 | 0 | 0 | 0 | 0.01957 | 0.01957 |
| hand | 6 660 | 0.24 | 0.20 | 0.03 | 0.01 | 32 / 1 | 33 | 66 | 6 | 3 | 3 | 1.576 | 1.576 |
| mannequin | 12 195 | 0.43 | 0.32 | 0.09 | 0.03 | 30 / 1 | 31 | 62 | 4 | 1 | 1 | 0.8287 | 0.8288 |
| rhino | 40 304 | 1.74 | 1.18 | 0.47 | 0.10 | 32 / 1 | 33 | 66 | 2 | 1 | 1 | 1.323 | 1.323 |
| tree | 33 592 | 1.34 | 0.91 | 0.39 | 0.04 | 25 / 1 | 26 | 62 | 4 | 4 | 4 | 2.017 | 2.017 |

("rest" includes symbolic analysis, assembly, line searches, certificate, extraction.
"merged (≤ 8)": clusters of at most 8 cones, the default limit; "merged (all)": no limit. The
last two columns are the distortion of the optimum and of the cones after merging all
clusters; merging costs 0.1 % on `axes`, 0.5 % on `bunny`, 2 % on `chair`, 10 % on
`camelshotnoise`, and less than 0.01 % elsewhere.) All 23 solves are certified optima reached
without a safeguard round; state residuals are `10⁻¹⁵` to `2·10⁻¹¹`. `brain` took 33.3 and
33.5 s in two runs of this table.

Where the time goes: 63 to 85 % is numeric factorization. Every Newton step
refactorizes the `2n × 2n` saddle-point matrix (about 0.5 s on `brain`), and there are 12 to 57
steps, all but one in Stage I: 3 to 8 per regularized problem for the first six or seven values
of `γ`, usually one for the last ones. The paper
reports 20 – 25 s for 100 – 150k triangles with 50 – 100 linear solves on a 2018 laptop; the
dataset has no mesh of that size, and on this faster machine 306k triangles take 33 s, so this
reference implementation is in the same range, not far below it. The Laplacian never changes
and can be factorized once; that is [the factor-once method](#the-factor-once-method), which is
the default when it applies.

## Performance of the factor-once method

`cargo run --release -p ddg-optimal-cones --example speedup -- --repeat 3 [--threads 8]`: every
dataset mesh at `λ = ½`, fastest of three runs, wall-clock seconds inside `solve`
(`report.timings.total`); AMD Ryzen Threadripper 3970X, shared with other jobs (load about 2 during
the runs). *Before* is `SolveMethod::Reference` with the schedule of the table above (the
default); *after* is `SolveMethod::Auto` with default options, which used the factor-once method
on every row. The reference method runs on one thread (threads do not help it: `rhino` 1.77 s on
one thread, 2.73 s on 16, prototype measurement). "Working set" is the number of columns the method computed.

| mesh | vertices | cones | working set | before: reference (s) | factorizations | after, 1 thread (s) | speedup | after, 8 threads (s) | speedup | max abs diff µ | max abs diff u | certificate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| airplane | 2349 | 6 | 34 | 0.117 | 32 | 0.005 | 23.4x | 0.005 | 22.8x | 7.2e-12 | 2.6e-12 | ok |
| axes | 3735 | 14 | 32 | 0.404 | 32 | 0.012 | 33.7x | 0.011 | 33.4x | 3.3e-11 | 9.2e-12 | ok |
| body | 11971 | 4 | 16 | 1.341 | 33 | 0.026 | 51.6x | 0.026 | 48.7x | 3.4e-12 | 6.6e-12 | ok |
| brain | 152893 | 123 | 546 | 32.882 | 50 | 7.162 | 4.6x | 3.563 | 9.1x | 1.8e-9 | 4.8e-10 | ok |
| bump | 3329 | 0 | 0 | 0.061 | 12 | 0.006 | 10.2x | 0.006 | 10.2x | 0.0e0 | 1.2e-14 | ok |
| bunny | 14290 | 254 | 365 | 1.852 | 44 | 0.333 | 5.6x | 0.212 | 8.9x | 1.9e-10 | 4.7e-11 | ok |
| camel | 11381 | 2 | 17 | 1.017 | 36 | 0.031 | 32.8x | 0.035 | 29.1x | 1.1e-9 | 1.5e-10 | ok |
| camel_noise1 | 11381 | 5 | 22 | 0.990 | 35 | 0.042 | 23.6x | 0.045 | 21.8x | 2.6e-9 | 3.7e-10 | ok |
| camel_noise2 | 11381 | 8 | 35 | 0.980 | 34 | 0.044 | 22.3x | 0.043 | 22.0x | 1.0e-10 | 3.8e-11 | ok |
| camelbad | 11392 | 2 | 17 | 0.859 | 35 | 0.032 | 26.8x | 0.034 | 25.5x | 1.8e-10 | 6.1e-11 | ok |
| camelbadnonuniform | 9092 | 1 | 16 | 0.596 | 36 | 0.023 | 25.9x | 0.026 | 22.9x | 1.2e-11 | 1.9e-11 | ok |
| camelshotnoise | 11513 | 452 | 651 | 1.690 | 57 | 0.482 | 3.5x | 0.270 | 6.5x | 4.2e-8 | 8.7e-9 | ok |
| chair | 5726 | 73 | 156 | 0.426 | 39 | 0.055 | 7.7x | 0.040 | 10.8x | 3.6e-10 | 5.1e-11 | ok |
| cow1 | 15049 | 2 | 16 | 0.863 | 28 | 0.043 | 20.1x | 0.045 | 19.1x | 4.1e-13 | 7.3e-13 | ok |
| face | 13746 | 0 | 0 | 0.409 | 12 | 0.030 | 13.6x | 0.030 | 14.0x | 0.0e0 | 8.7e-15 | ok |
| hand | 3372 | 6 | 36 | 0.250 | 33 | 0.009 | 27.8x | 0.009 | 27.6x | 6.6e-9 | 1.2e-9 | ok |
| mannequin | 6130 | 4 | 20 | 0.391 | 31 | 0.021 | 18.6x | 0.020 | 19.7x | 1.2e-11 | 4.0e-12 | ok |
| rhino | 20205 | 2 | 17 | 1.654 | 33 | 0.069 | 24.0x | 0.065 | 26.0x | 6.8e-11 | 1.7e-11 | ok |
| tree | 16917 | 4 | 16 | 1.329 | 26 | 0.050 | 26.6x | 0.051 | 26.3x | 4.9e-13 | 2.6e-12 | ok |

"max abs diff" compares the automatic result with the reference result (one thread). All
results are certified (`certificate.holds(1e-9)`). The speedup is `3.5×` (`camelshotnoise`, 452 cones
and a working set of 651) to `52×` (`body`); `brain` takes 7.2 s on one thread and 3.6 s on 8
(32.9 s before). Threads help only where there are many columns to solve (`brain`, `bunny`,
`camelshotnoise`, `chair`: `1.4×` to `2×`); elsewhere the time is the factorization and the
setup. The difference on `brain` at 8 threads is `|Δµ| = 4·10⁻⁷`, `|Δu| = 8·10⁻⁸` (the single-thread run
`2·10⁻⁹` and `5·10⁻¹⁰`: the batches of columns differ, and the certificates allow differences of
this size). Where the memory goes: `brain` stores 546 columns of 153 k doubles, 0.67 GB.

**Where the method stops paying: `λ` sweep on `cow1`** (15 049 vertices), reference method, the
factor-once method with its limits lifted, and the automatic mode, fastest of three, seconds:

| λ | cones | working set | reference | factor-once only, 1 thread | automatic, 1 thread | 8 threads: factor-once only | 8 threads: automatic |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 2 | 16 | 0.88 | 0.043 | 0.049 | 0.043 | 0.045 |
| 0.1 | 7 | 51 | 1.12 | 0.062 | 0.062 | 0.050 | 0.055 |
| 0.03 | 8 | 51 | 1.16 | 0.065 | 0.065 | 0.057 | 0.058 |
| 0.01 | 109 | 196 | 1.34 | 0.173 | 0.165 | 0.103 | 0.102 |
| 0.003 | 402 | 591 | 1.36 | 0.534 | 0.555 | 0.237 | 0.249 |
| 0.001 | 691 | 956 | 1.37 | 1.01 | 1.05 | 0.455 | 0.449 |
| 0.0003 | 1 132 | 1 539 | 1.39 | 2.08 | **2.80** (fallback) | 0.96 | 0.95 |
| 0.0001 | 1 682 | 2 213 | 1.47 | 4.18 | **3.05** (fallback) | 2.10 | 2.09 |

On one thread the factor-once method stops winning at a working set of about 1 200 (4 % of the
vertices are cones; the dense active-set part is `O(a²)` per iteration and sequential, the Gram
product `O(n a²)`), on 8 threads at about 1 900 (10 % cones) at `cow1`. On one thread `rhino`
(20 205 vertices), `bunny` (14 290) and `camelshotnoise` (11 513) were still ahead of the
reference method at the largest working sets tried (1 148, 1 194, 1 134; extrapolated crossovers
about 1 300 to 1 400), and `brain` (152 893) took 25 s at a working set of 1 613 and 59.5 s at
3 008 against 32.9 s for the reference method (crossover about 1 900); `brain` at 8 threads, in an
earlier run of the same method, took 17.5 s at 3 008, still ahead. The limit
`1300 · (n/15000)^0.2 · (1 + (threads − 1)/7)` (at most twice) is a fit of these points, chosen so
that the working set stops growing about where the two methods take equal time. It is **not** a
guarantee that automatic mode is never slower than the better method: a problem whose working
set grows just past the limit has paid for the attempt and then runs the reference method.
Measured worst cases: automatic mode is `2.0×` to `2.1×` the reference time on one thread right
beyond the limit (`cow1` at `λ = 3·10⁻⁴` and `10⁻⁴`: 2.80 s and 3.05 s against 1.39 s and 1.47 s;
`brain` at `λ = 3·10⁻⁵`: 65 s against 33 s), and `1.4×` at 8 threads beyond the larger limit
(`cow1`, `λ = 10⁻⁴`: 2.09 s against 1.53 s); one thread, working set 1 194 on `bunny`, just below
the limit, 1.46 s against 1.80 s for the reference method. The `λ` range that the UV pipeline and the
paper use (`0.1` to `1`) is far from the limit on every dataset mesh.

Other meshes, one thread (reference / automatic, seconds), `λ = 0.01 … 0.0003`: `rhino`
2.18 / 0.18, 2.33 / 0.86, 2.40 / 0.95, 2.55 / 1.82; `bunny` 1.88 / 0.38, 1.84 / 0.51, 1.89 / 0.82,
1.80 / 1.46; `camelshotnoise` 1.79 / 0.81, 1.95 / 0.73, 1.81 / 0.80, 1.96 / 1.03. Automatic mode
used the factor-once method in all twelve.

**The UV pipeline** (`cargo bench -p ddg-uv`, Criterion, one thread, `DDG_BENCH_LARGE=1`; *before* is
the tree at commit `584ddef`, i.e. this crate without the factor-once method, *after* the tree with
it; median of 10 samples, `--save-baseline before` / `--baseline before`): the cone solves of the
pipeline ask for positive-only angles capped at `0.95·2π`, with free sites on round shapes.

| benchmark | before | after | change |
|---|---:|---:|---:|
| `unwrap/bumpy_disk_64` (8 192 triangles, boundary) | 87.9 ms | 27.7 ms | −69 % |
| `unwrap/uv_sphere_64x32` (closed, round, free-site fallback) | 374.6 ms | 468.6 ms | **+25 %** |
| `unwrap/grid_cube_16` | 125.2 ms | 12.7 ms | −90 % |
| `unwrap/torus_96x24` (no cone solve) | 15.1 ms | 14.3 ms | −8 % |
| `unwrap/hand` | 492 ms | 33.2 ms | −93 % |
| `unwrap/camel` | 1.542 s | 157 ms | −90 % |
| `unwrap/bunny` | 3.111 s | 536 ms | −83 % |
| `unwrap/brain` (305 782 triangles) | 35.07 s | 9.74 s | −72 % |
| `stages/camel/2_cones` | 1.413 s | 100 ms | −93 % |
| `stages/brain/2_cones` | 34.16 s | 8.85 s | −74 % |

(the other stages are unchanged within noise). The slower one is `uv_sphere_64x32`: its first cone
solve (the penalized problem, positive angles only) has a cone at every one of its 1 986 vertices,
the working set reaches the limit (868 for this size), the method gives up and the reference method
runs; the attempt costs about 90 ms. The second solve there (8 free sites, cones forbidden
elsewhere) is fast. Closed meshes in this pipeline only work because the method starts from sites
clamped at their bounds when no single site can carry the angle sum `4π` (otherwise `bunny` and
`brain` would have fallen back as unsupported). A caller that knows its first solve will fill the
surface can ask for `SolveMethod::Reference` there. The default schedule matters here too: with
the short schedule `10⁻², 10⁻⁶, 10⁻¹⁰` as the default, this benchmark's `unwrap/bunny` went from
3.1 s to 7.5 s and `unwrap/brain` from 35 s to 99 s (the reference method, with capped positive
angles; see [the schedule](#the-regularization-schedule)), which is why it is not the default.

## Formula tags

Functions that implement a formula carry `[formula: <id>]`; `scripts/check-formulas.sh` matches
the ids with `lean/proofs`. Ids taken from there: `soft-threshold-fixed-point`,
`soft-threshold-active-sets`, `abs-subdifferential-bounded`, `cone-fixed-point-bounded`,
`regularized-cone-update`, `optimality-system-sufficient`, `duality-gap-certificate`,
`newton-stopping-rule`, `closed-mesh-cone-sum`, `closed-mesh-optimal-constant`,
`lambda-area-normalization`, `cotan-laplacian`, `discrete-yamabe-linearization`.

Ids introduced by this crate, with the statement the function implements. `lean/proofs` has
theorems for them by now (`DdgProofs/Solver.lean`, `DdgProofs/Certificate.lean`; one is cited
in `CITED.md`), except for the four `reduced-…` ids of the factor-once method, which are new and
have neither a theorem nor an entry in `CITED.md` yet (`lean/` was not touched when they were
added, so `scripts/check-formulas.sh` reports them as missing). The last column says
what the theorems cover, which is not always the whole statement.

| id | function | statement | in `lean/proofs` |
|---|---|---|---|
| `duality-gap-inexact` | `certify` | For any `u, φ, µ` vanishing on the boundary, with `r = L u − µ + Ω`, `φ̂` in the domain of every `g_i*`, `ũ_i = −(L φ̂)_i/w_i` where `w_i > 0`, `J = ½ Σ w_i u_i² + Σ t_i\|µ_i\|` and `D = −½ Σ_{w_i>0} (L φ̂)_i²/w_i + Σ φ̂_i Ω_i − Σ g_i*(φ̂_i)`: `J − D = ½ Σ_{w_i>0} w_i (u_i − ũ_i)² + Σ (t_i\|µ_i\| + g_i*(φ̂_i) − φ̂_i µ_i) + Σ_{w_i=0} (L φ̂)_i u_i − Σ φ̂_i r_i`. (At a vertex with fixed angle `a`: `t_i\|µ_i\|` is `0` and `g_i*(p) = a p`.) | the identity: `gap_identity_general`, `objective_sub_dualObjectiveOfAdjoint` |
| `boundary-optimality-system` | `Discrete::with_boundary`, `solve_with_boundary`, `certify_with_boundary` | With `Ω̂_i = Ω_i` (`i ∈ I`), `Ω̂_i = k_i = π − Σ corner angles` (`i ∈ B`), unknowns `u ∈ R^V`, `s = (µ, ν)` and one site `g_i` per vertex (`t_i\|·\|` with bounds): `minimize ½ Σ_V wE_i m_i u_i² + Σ_V g_i(s_i)` subject to `(L u)_i = s_i − Ω̂_i` for all `i`. Its optimality system is `(L u)_i = s_i − Ω̂_i`, `(L φ)_i = −wE_i m_i u_i`, `φ_i ∈ ∂g_i(s_i)`; it is sufficient for a global minimum, and `Σ s = Σ Ω̂ = 2πχ`. Sites on `B`: `Free` `g = 0` (so `φ_i = 0`), `Penalized` `g = λ wR_i\|·\|`, `Convex` `g = ι_[0,∞)`. | the sufficiency is `optimality-system-sufficient` for these data (a `ConeProblem` with this `Ω̂` and sites); that `ν_i = (L u)_i + k_i` is the geodesic curvature of the flattened boundary (paper Sec. 8.3) is a modelling statement, checked numerically, not proved |
| `saddle-point-newton-system` | `SaddleSolver` | The linear system of a choice of pieces is `W u + A φ = 0` together with, per vertex, `(A u)_i = µ_i − Ω_i` (angle fixed), `(A u)_i − k_i (φ_i − o_i) = b_i − Ω_i` (angle `b_i + k_i(φ_i − o_i)`), or `φ_i = p_i` (angle free); in the unknowns `(u, φ)` its matrix is the symmetric `[[W, A], [A, −K]]` with the rows and columns of the pinned `φ_i` replaced by the identity. | `isSolution_iff`, `isAssembledSolution_iff`; a solution with every vertex on its piece is a minimizer (`isMinimizer_of_isSolution`, `isProxMinimizer_of_isSolution`) |
| `gauss-bonnet-adjoint-shift` | `gauss_bonnet_shift` | For `k_i > 0`, `h(c) = Σ_i clamp(k_i P_{t_i}(x_i + c), lo_i, hi_i)` is continuous, nondecreasing and piecewise linear, with breakpoints among `κ − x_i` for the kinks `κ ∈ {t_i, t_i + hi_i/k_i, −t_i, −t_i + lo_i/k_i}`. The function returns the `c` of smallest magnitude with `h(c) = target`, or, if the target is not attained, the `c` of smallest magnitude at which `h` is nearest to it. | monotonicity (`sum_proxReg_add_monotone`) and that the pinned system is the full one once `h(c) = Σ Ω` (`isSolution_of_pinned_of_sum_eq`); not the search on the breakpoints |
| `adjoint-constant-interval` | `adjoint_constant_interval` | With `∂g_i(µ_i) = [a_i, b_i]`: `{c : φ_i + c ∈ ∂g_i(µ_i) for all i} = [max_i (a_i − φ_i), min_i (b_i − φ_i)]`. | `kkt_add_const_iff`, `Site.stationary_iff_mem_subdiff` |
| `regularized-dual-objective` | `DualSlope` | For `w > 0`, stiffness `k_i > 0`, centre `µ̄`, with `Φ(φ) = ½ Σ w_i u_i² − Σ Ω_i φ_i + Σ e_i(φ_i)`, `W u = −A φ`, `e_i'(p) = clamp(k_i P_{t_i}(p + µ̄_i/k_i), lo_i, hi_i)` (the intended reading: `Φ` is the convex dual objective of `minimize J(u, µ) + Σ (µ_i − µ̄_i)²/(2 k_i)`, `∇Φ = −(A u − µ + Ω)`; this reading is not claimed as proved). Along `(du, dφ)` with `W du = −A dφ`: `d/dα Φ = Σ w_i du_i (u_i + α du_i) − Σ Ω_i dφ_i + Σ e_i'(φ_i + α dφ_i) dφ_i`, nondecreasing and piecewise linear in `α`. | the slope along a direction and its monotonicity (`dualSlope_eq`, `dualSlope_monotone`); `Φ` itself is not formalized |
| `cone-cluster-frechet-mean` | `frechet_order` | For a set of cone vertices connected through edges between its members, with angles `µ_k`: the member `v` that minimizes `Σ_k \|µ_k\| dist(v, c_k)²`, `dist` the length of the shortest path along edges with both ends in the set; ties go to the smaller index. | existence of a minimizer, and that merging at any member keeps `Σ µ` and does not increase `Σ \|µ\|` (`exists_frechetMean`, `sum_mergeCones`, `sum_abs_mergeCones_le`) |
| `reduced-lasso-problem` | `Reduced::new`, `factor_once::attempt` | Let `A` be positive definite (the Laplacian block of the unknowns, or `L` without the row and column of the pin), `w_i > 0`, `S` the free sites with `z_j = A⁻¹ e_j`, `F` the sites with a fixed angle `a_F`, `w0 = A⁻¹ (Ω − a_F)`. Then `u = Z µ_S − w0` solves `A u = µ − Ω` (with `µ_F = a_F`), and the canonical problem is equivalent to `minimize ½ µ_Sᵀ G µ_S + qᵀ µ_S + Σ_{j∈S} t_j \|µ_j\|` subject to `lo_j ≤ µ_j ≤ hi_j` (and, for a pin, `Σ_S µ = Σ (Ω − a_F)`), with `G = Zᵀ W Z`, `q = −Zᵀ W w0`, up to the constant `½ w0ᵀ W w0`. | not in `lean/proofs` (see below) |
| `reduced-gram-matrix` | `Reduced::gram_cross` | `G_jk = z_jᵀ W z_k` and `q_j = −z_jᵀ W w0` with `z_j = P A_p⁻¹ e_j` (`A_p⁻¹` padded with a zero at the pin `p`, `P x = x − (Σ w_i x_i / Σ w_i) 1`; without a pin `P = I`); `z_p = 0`. `G` is positive semidefinite, positive definite on every set of linearly independent columns. | not in `lean/proofs` |
| `reduced-optimality-conditions` | `Reduced::gradient`, `violation_at_zero`, `Dense::exact_violation` | With `u = Z µ_S − w0` (so `Σ w u = 0` with a pin) and `y = A⁻¹ (W u)` (zero at the pin): `∂/∂µ_j (½ uᵀ W u) = z_jᵀ W u = y_j`. `µ_S` minimizes the reduced problem iff there is `ν ∈ R` (`ν = 0` without the sum constraint) with `−(y_j + ν) ∈ ∂g_j(µ_j)` for every `j ∈ S`, `g_j(m) = t_j\|m\| + ι_[lo_j, hi_j](m)`; then `φ = −(y + ν)` satisfies `L φ = −W u` and `φ_j ∈ ∂g_j(µ_j)`. At `µ_j = 0` the condition is `\|y_j + ν\| ≤ t_j` (restricted to the allowed signs); at `µ_j = hi_j > 0` it is `y_j + ν + t_j ≤ 0`; at `µ_j = lo_j < 0` it is `y_j + ν − t_j ≥ 0`. | not in `lean/proofs` (it is `optimality-system-sufficient` for the reduced variables) |
| `reduced-sum-constraint-pin` | `Reduced::new`, `Reduced::project` | If `A` is positive semidefinite with kernel exactly the constants, `p` an index, `A_p` positive definite, and `Σ b = 0`, then `x = P pad(A_p⁻¹ b_p̂)` solves `A x = b` (the equation of `p` follows from the others because `1ᵀ A = 0` and `Σ b = 0`) and, among `x + c 1`, `x` is the one with `Σ w x = 0`, which minimizes `½ Σ w (x + c)²`. For `b = µ − Ω` the condition is `Σ µ = Σ Ω`, i.e. `Σ_S µ_S = Σ (Ω − a_F)`. | not in `lean/proofs` |
| `proximal-fixed-point` | `solve` | If `(u, µ)` minimizes `J(v, ν) + Σ (ν_i − µ̄_i)²/(2 k_i)` over the feasible pairs and `µ = µ̄`, then `(u, µ)` minimizes `J`. | `IsProxMinimizer.isMinimizer_of_fixedPoint`; also that a proximal step does not increase the objective (`IsProxMinimizer.objective_add_le`) |
| `componentwise-backward-error` | `backward_error` | `max_i \|b − K x\|_i / (\|K\| \|x\| + \|b\|)_i` is the smallest `ε` such that `(K + δK) x = b + δb` for some `\|δK\| ≤ ε\|K\|`, `\|δb\| ≤ ε\|b\|` entrywise (Oettli and Prager 1964). | `oettli_prager` |
| `component-separation` | `Instance::part`, `Discrete::from_instance`, `certify_instance`, `solve`, `solve_part`, `evaluate_cones`, `reoptimize_cones` | On a mesh with several components `L` is block diagonal and the objective a sum: `(u, µ)` is a minimizer (a solution of the optimality system) of the whole problem iff every restriction is one of its component's part, with the masses of the whole problem; every closed component has its own constant of `u` and its own `Σ µ = 2πχ`. | `isMinimizer_sum_iff`, `kkt_sum_iff`, `objective_sum`, `feasible_sum_iff`, `isMinimizer_sum_iff_alone` (`DdgProofs/Components.lean`) |
| `component-lambda-rescaling` | `equivalent_lambda` | A component with its masses normalized by its own area `A_c` instead of the total `A` has the same minimizers with `λ · A / A_c`, and the adjoint of its optimality system is `A / A_c` times that of the whole problem. | `isMinimizer_normalizeMass_alone_iff`, `kkt_normalizeMass_alone_iff`, `isMinimizer_scaleMass_withLam_iff`, `kkt_scaleMass_withLam_iff`, `Site.stationary_mul_iff` |
| `feature-size-weights` | `feature_size_weights` | `wE_i = r_i`, `wR_i = 1/r_i` with `r_i = \|Ω_i\|/m_i + ε` ("as described") or its reciprocal ("as printed"), `m` area-normalized (paper Sec. 8.1). | cited in `CITED.md` (a choice of input data, not a derived formula) |
| `cone-feedback-mismatch` | `scale_mismatch` | `R_i = Σ_{t∋i} A_t r_t / Σ_{t∋i} A_t` over the measured faces, `δ_i = (u_i − R_i) − c` with `c` the mass-weighted median of `u − R` over the measured vertices of the component; `δ` is invariant under `u + a`, `r + b`, and zero where `u − R` is constant. | `faceMean_shift`, `mismatch_shift`, `mismatch_eq_zero_of_const` (`DdgProofs/ConeFeedback.lean`) |
| `cut-weighted-quantile` (taken from the cut flow) | `weighted_quantile` (feedback) | the first value, in ascending order, at which the weights summed so far reach `q` times their sum; here `q = ½`. Shifting the values shifts it. | the cut flow's theorems, and `isWeightedQuantile_shift`, `weightedQuantile_shift` |
| `cone-feedback-rule` | `Feedback::observe` | the reweighting rule of [the feedback loop](#the-feedback-loop) | cited in `CITED.md` as a design choice: a threshold rule, nothing to prove |

## Tests

```
cargo test --release -p ddg-optimal-cones            # everything; the dataset tests need --release
cargo test --release -p ddg-optimal-cones --test dataset -- --nocapture
cargo run  --release -p ddg-optimal-cones --example published -- refit     # or facts, adjoint, sweep, iterate
cargo run  --release -p ddg-optimal-cones --example baseline -- --merge --repeat 3   # the reference method
cargo run  --release -p ddg-optimal-cones --example speedup -- --repeat 3               # factor-once vs reference table
cargo run  --release -p ddg-optimal-cones --example boundary -- table     # or spikes
cargo run  --release -p ddg-optimal-cones --example fingerprint            # bit hashes of default-mode solves
```

| file | what it checks |
|---|---|
| unit tests in `src/` | the scalar formulas (pieces, proximal map, subdifferential, conjugate), the adjoint shift |
| `tests/convex_reference.rs` | agreement with the generic convex solver Clarabel on random problems with every feature |
| `tests/oracle.rs` | replay of the Lean oracle's golden fixtures and live differential fuzzing against its executable (skipped with a notice if it is not built) |
| `tests/properties.rs` | invariance under scaling and relabeling, monotonicity in `λ`, the guarantee of Sec. 4.5, exact constraints, the certificate of the returned triple |
| `tests/api.rs` | errors, options, tiny meshes, extraction and merging, the certificate identity, non-converged status |
| `tests/dataset.rs` | the published dataset: every mesh solved to a certified optimum, the findings above, the boundary modes against the double and the published refits |
| `tests/boundary.rs` | the boundary treatments against the double of the mesh, orderings (free ≤ the others), reporting, certificate, extraction |
| `tests/factor_once.rs` | the factor-once method against the reference method on random problems of every topology, steering and boundary treatment, with tolerances derived from the certificates; every fallback reason forced; what falls back |
| `tests/enclosure.rs` | `optimum_enclosure` against the generic convex solver on perturbed triples; when it is withheld |
| `tests/components.rs` | disconnected meshes: the union against separate solves with `λ · A / A_c`, against Clarabel, relabeling, steering per component, methods, threads, every entry point |
| `tests/reoptimize.rs` | `reoptimize_cones` against `solve` on `reoptimization_problem` and against Clarabel, moving a cone, edge cases, disconnected meshes, the published cone files |
| unit tests in `src/feedback.rs` | the weighted median, the mismatch (no mismatch for a consistent realization, invariance under constants, unmeasured faces, errors), the rule (a far vertex is reweighted once, the limits end the loop without applying the round, a failed round is undone, factors and initial weights, boundary vertices, a component keeps a weight), invalid options |
| `crates/ddg-uv/tests/feedback.rs` | the loop with Boundary First Flattening: it finds the tips of thin spikes on a grid and on a sphere, ends where the hand-set exclusion starts, changes nothing on smooth meshes, invalid options are errors |

The dataset is not part of the repository; without it the dataset tests print a notice and pass.

**Which method the tests use.** Everything runs in automatic mode by default (so the factor-once
method on all small meshes). The files that compare with an independent solver (`convex_reference`,
`oracle`, `properties`, `enclosure`, `dataset`) solve through `tests/common::solve_t`, which reads
`DDG_CONES_METHOD`: `reference` forces the reference method, `fast` forces the factor-once
method wherever `factor_once_support` allows it (with the limits for giving up lifted, so that
the test checks correctness and not thresholds; the problems it does not support go to the
reference method), anything else leaves the options alone:

```
DDG_CONES_METHOD=fast      cargo test --release -p ddg-optimal-cones
DDG_CONES_METHOD=reference cargo test --release -p ddg-optimal-cones
```

The tests in `api.rs` that read the reference method's stages, counters and options
(`report_counts_what_was_done`, `starved_solver_reports_not_converged`,
`every_configuration_finds_the_same_optimum`) force `SolveMethod::Reference` themselves.

## Not implemented

* The boundary treatment is not a field of `ConeProblem` (that would break struct literals in
  other crates); `evaluate_cones` has no boundary variant (use `solve_with_boundary` with
  `cone_costs` of `0` and `+∞`); `Penalized` and `Convex` have no per-vertex bounds of their own;
  a zero distortion weight at a `Free` boundary vertex can make the problem ill posed (the solver
  then says `NotConverged`); the `Infeasible` error text speaks of a "closed mesh" also for a
  boundary treatment other than `Fixed`.
* Multiresolution (Sec. 6.3).
* The factor-once method does not handle zero distortion weights (the Gram matrix is singular),
  `λ = 0`, and it has no proof of convergence (the certificate and the fallback replace it); the
  re-optimization of the boundary for extracted cones (`with_fixed_interior`) uses the reference
  method; the dense part is sequential and `O(a²)` per iteration, which is why it loses to the
  reference method above about 5 % of the vertices being cones; there is no GPU path (assessed,
  not worth it).
* Non-manifold meshes, intrinsic Delaunay refinement.
* Disconnected meshes: the components are solved concurrently with one thread each; a large
  component next to small ones does not get the pool's threads inside its own solve.
* `reoptimize_cones` with a `Penalized` boundary has an `ℓ¹` term left and always runs the generic
  solver; with an active bound it solves the unconstrained system first and then the generic
  solver (the work of the direct solve is lost). The limit of 32 free sites between the reduced
  and the direct method is a constant measured on the dataset, not adapted to the mesh.
