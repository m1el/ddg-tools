# ddg-optimal-cones

Where to put cone singularities, and with which angles, so that a triangle mesh can be flattened
conformally with the least area distortion for the total cone angle spent. An implementation of

> Yousuf Soliman, Dejan Slepčev, Keenan Crane.
> *Optimal Cone Singularities for Conformal Flattening.*
> ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201367>

("the paper"). The method, the optimality system and the two-stage semismooth Newton algorithm
are theirs. What this crate adds is listed under [Deviations from the paper](#deviations-from-the-paper):
a certificate on every solve, safeguards, closed meshes, angle bounds, prescribed cones.

Contents: [example](#example) · [the problem](#the-problem) · [the boundary](#the-boundary) · [interface](#interface) ·
[steering](#steering) · [what comes back](#what-comes-back) · [the certificate](#the-certificate) ·
[algorithm](#algorithm) · [deviations from the paper](#deviations-from-the-paper) ·
[the published dataset](#the-published-dataset) · [choosing λ](#choosing-λ) ·
[closed meshes](#closed-meshes) · [needle-like spikes](#needle-like-spikes) ·
[performance baseline](#performance-baseline) · [formula tags](#formula-tags) · [tests](#tests) ·
[not implemented](#not-implemented)

## Example

```rust
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{ConeProblem, SolveOptions, Status, solve};

let mesh = TriMesh::new(positions.len(), triangles)?;        // connected, manifold, oriented
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
    pub regularization: Vec<f64>,                // γ of Stage I; default 1, 1e-1, …, 1e-10
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
}
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
    pub stages: Vec<StageReport>,     // kind, safeguard, newton_iterations, damped_steps,
                                      // constant_shifts, stop: StageStop, seconds
    pub regularized_iterations: usize, pub exact_iterations: usize,
    pub safeguard_rounds: usize, pub safeguard_iterations: usize,
    pub linear_solver: LinearSolverReport,   // factorizations, solves, refinement, times
    pub timings: Timings,                    // setup, regularized, exact, certificate, extraction, total
    pub bound_correction: f64, pub zeroed_angles: usize, pub largest_zeroed_angle: f64,
    pub certificate: Certificate,
    pub message: Option<String>,
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
pub struct ConeEvaluation { log_scale, adjoint, distortion, total_cone_angle, gauss_bonnet_defect }

// extract.rs: what extraction did
pub struct ExtractionReport { threshold, n_dropped, dropped_angle, merge_requested, n_clusters,
    merged: Vec<MergedCluster>,       // members, cone, diameter
    unmerged: Vec<UnmergedCluster>,   // vertices, angle, diameter, reason: TooLarge | TooWide | BoundsViolated
    cones_differ_from_measure, distortion_change, total_angle_change, n_beyond_two_pi }

// weights.rs: the curvature-adapted weights of the paper's Sec. 8.1
pub fn feature_size_weights(mesh, geometry, epsilon: f64, convention: FeatureSizeConvention)
    -> Result<FeatureSizeWeights, ConeError>;      // AsDescribed | AsPrinted

pub fn soft_threshold(t: f64, x: f64) -> f64;      // P_t(x), paper Eq. 18
pub enum ConeError { Geometry, NotConnected, InvalidLambda, WrongLength, InvalidDistortionWeight,
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
| re-optimize the angles of a given set of cones | `cone_costs = 0` at those vertices and `+∞` elsewhere ([the published dataset](#the-published-dataset)) |
| keep every angle realizable | bounds inside `(−2π, 2π)`; `extraction.n_beyond_two_pi` counts what is outside |
| ignore needle-like spikes | `distortion_weights[v] = 0.0` at their tips ([needle-like spikes](#needle-like-spikes)) |
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

## Algorithm

The paper's (Fig. 12). **Stage I** solves the regularized problems for `γ = 1, 10⁻¹, …, 10⁻¹⁰`,
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
16. **Connected meshes only**; split components first.

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

## Performance baseline

`cargo run --release -p ddg-optimal-cones --example baseline -- --merge --repeat 3`: every
dataset mesh at `λ = 0.5`, default options, one thread, fastest of three runs. AMD Ryzen
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
and could be factorized once; that path is not implemented.

## Formula tags

Functions that implement a formula carry `[formula: <id>]`; `scripts/check-formulas.sh` matches
the ids with `lean/proofs`. Ids taken from there: `soft-threshold-fixed-point`,
`soft-threshold-active-sets`, `abs-subdifferential-bounded`, `cone-fixed-point-bounded`,
`regularized-cone-update`, `optimality-system-sufficient`, `duality-gap-certificate`,
`newton-stopping-rule`, `closed-mesh-cone-sum`, `closed-mesh-optimal-constant`,
`lambda-area-normalization`, `cotan-laplacian`, `discrete-yamabe-linearization`.

Ids introduced by this crate, with the statement the function implements. `lean/proofs` has
theorems for them by now (`DdgProofs/Solver.lean`, `DdgProofs/Certificate.lean`; one is cited
in `CITED.md`), and `scripts/check-formulas.sh` reports no missing id. The last column says
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
| `proximal-fixed-point` | `solve` | If `(u, µ)` minimizes `J(v, ν) + Σ (ν_i − µ̄_i)²/(2 k_i)` over the feasible pairs and `µ = µ̄`, then `(u, µ)` minimizes `J`. | `IsProxMinimizer.isMinimizer_of_fixedPoint`; also that a proximal step does not increase the objective (`IsProxMinimizer.objective_add_le`) |
| `componentwise-backward-error` | `backward_error` | `max_i \|b − K x\|_i / (\|K\| \|x\| + \|b\|)_i` is the smallest `ε` such that `(K + δK) x = b + δb` for some `\|δK\| ≤ ε\|K\|`, `\|δb\| ≤ ε\|b\|` entrywise (Oettli and Prager 1964). | `oettli_prager` |
| `feature-size-weights` | `feature_size_weights` | `wE_i = r_i`, `wR_i = 1/r_i` with `r_i = \|Ω_i\|/m_i + ε` ("as described") or its reciprocal ("as printed"), `m` area-normalized (paper Sec. 8.1). | cited in `CITED.md` (a choice of input data, not a derived formula) |

## Tests

```
cargo test --release -p ddg-optimal-cones            # everything; the dataset tests need --release
cargo test --release -p ddg-optimal-cones --test dataset -- --nocapture
cargo run  --release -p ddg-optimal-cones --example published -- refit     # or facts, adjoint, sweep, iterate
cargo run  --release -p ddg-optimal-cones --example baseline -- --merge --repeat 3
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
| `tests/enclosure.rs` | `optimum_enclosure` against the generic convex solver on perturbed triples; when it is withheld |

The dataset is not part of the repository; without it the dataset tests print a notice and pass.

## Not implemented

* The boundary treatment is not a field of `ConeProblem` (that would break struct literals in
  other crates); `evaluate_cones` has no boundary variant (use `solve_with_boundary` with
  `cone_costs` of `0` and `+∞`); `Penalized` and `Convex` have no per-vertex bounds of their own;
  a zero distortion weight at a `Free` boundary vertex can make the problem ill posed (the solver
  then says `NotConverged`); the `Infeasible` error text speaks of a "closed mesh" also for a
  boundary treatment other than `Fixed`.
* A function that re-optimizes the angles of given cones; it is a call of `solve` with
  `cone_costs` of `0` and `+∞`.
* Multiresolution (Sec. 6.3), and a solver that factorizes the Laplacian once.
* Disconnected meshes, non-manifold meshes, intrinsic Delaunay refinement.
