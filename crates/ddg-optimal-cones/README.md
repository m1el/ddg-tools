# ddg-optimal-cones

Where to put cone singularities, and with which angles, so that a triangle mesh can be flattened
conformally with the least area distortion for the total cone angle spent. An implementation of

> Yousuf Soliman, Dejan Slepčev, Keenan Crane.
> *Optimal Cone Singularities for Conformal Flattening.*
> ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201367>

("the paper"). The method, the optimality system and the two-stage semismooth Newton algorithm
are theirs. What this crate adds is listed under [Deviations from the paper](#deviations-from-the-paper):
a certificate on every solve, safeguards, closed meshes, angle bounds, prescribed cones.

Contents: [example](#example) · [the problem](#the-problem) · [interface](#interface) ·
[steering](#steering) · [what comes back](#what-comes-back) · [the certificate](#the-certificate) ·
[algorithm](#algorithm) · [deviations from the paper](#deviations-from-the-paper) ·
[the published dataset](#the-published-dataset) · [choosing λ](#choosing-λ) ·
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
`ddg-core` is used with `default-features = false`, so nothing GPL-derived is linked.

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
    pub linear_solver: LinearSolver,             // Auto | SaddleLdlt | SaddleLu
    pub parallelism: Parallelism<'a>,            // Sequential | Pool(&ThreadPool)
}
pub struct Solution {
    // the optimum of the problem (this is what the certificate is about)
    pub measure: Vec<f64>,            // µ per vertex, exactly 0 where there is no cone
    pub log_scale_optimal: Vec<f64>,  // u per vertex
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
pub struct Certificate {
    pub state_residual: f64, pub adjoint_residual: f64, pub cone_residual: f64,
    pub bound_violation: f64, pub dual_infeasibility: f64, pub gauss_bonnet_residual: f64,
    pub primal_objective: f64, pub dual_objective: f64,
    pub duality_gap: f64, pub complementarity_gap: f64, pub gap_uncertainty: f64,
}
impl Certificate { pub fn holds(&self, tolerance: f64) -> bool }

// evaluate.rs: scale factor and adjoint of a given cone configuration
pub fn evaluate_cones(problem, cones: &[Cone], par: Parallelism<'_>) -> Result<ConeEvaluation, ConeError>;
pub struct ConeEvaluation { log_scale, adjoint, distortion, total_cone_angle, gauss_bonnet_defect }

// extract.rs: what extraction did
pub struct ExtractionReport { threshold, n_dropped, dropped_angle, merge_requested, n_clusters,
    merged: Vec<MergedCluster>, unmerged: Vec<UnmergedCluster>, cones_differ_from_measure,
    distortion_change, total_angle_change }

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
  a cone if `|µ_i| > cone_threshold`. With `merge_cone_clusters` every edge-connected cluster of
  cones (of at most `max_merge_cluster_size` cones) is replaced by one cone with the summed
  angle at the cluster's Fréchet mean vertex. If extraction changed anything
  (`extraction.cones_differ_from_measure`), the cones are deliberately not optimal, `log_scale`
  and `distortion` are recomputed for them, nothing certifies them, and `extraction` lists every
  difference including the cost in distortion.

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

* **Proved** (exact arithmetic, all `wE_i m_i > 0`): `D` is a lower bound of the optimal value
  whatever the residuals are (`duality-gap-certificate`, theorem `adjoint_gap_certificate`). If
  `(u, µ)` is exactly feasible, `objective − optimum ≤ duality_gap = complementarity_gap`.
* **Not proved**: the returned pair satisfies the state equation only up to `state_residual`, and
  an infeasible pair can have an objective below the optimum, so `duality_gap` can be negative.
  It is `complementarity_gap` (a sum of nonnegative terms) shifted by at most
  `gap_uncertainty = Σ_{w_i=0} |(L φ̂)_i u_i| + Σ |φ̂_i r_i|`. Read the gap as
  "`complementarity_gap`, give or take `gap_uncertainty`".
* **Where `wE_i = 0`** the dual problem requires `(L φ̂)_i = 0` exactly, which is not enforced by
  construction. Such a vertex contributes nothing to `D`; its violation is in `adjoint_residual`
  and, times `u_i`, in `gap_uncertainty`. `D` is then a lower bound only up to that term.
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
    default, leaves clusters of more than `max_merge_cluster_size` cones alone, never moves a
    prescribed cone, and keeps the merged cone within the angle bounds of its vertex. The
    Fréchet mean is the member minimizing `Σ_k |µ_k| dist(v, c_k)²` with distances along edges
    inside the cluster; the paper does not define it (on `axes` one of four merged cones lands
    on a different vertex of its cluster than the authors').
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
which this crate does not implement; the agreement there is `10⁻²`, not `10⁻⁴`, so some detail
of the authors' boundary treatment is not captured (uniform instead of lumped masses, and
other weights of the boundary mass, were tried and are worse or no better).

Not explained by re-optimized angles: `bunny`; `camelshotnoise` (with our operators its
spikes dominate the energy, the distortion of the file's cones is about `10⁵`); and the two
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
| camelshotnoise | 6 (20.1) | none | | 452 (43) |
| hand | 29 (52.5) | none | | 6 (3) |
| rhino | 36 (65.4) | none | | 2 (1) |

Remarks for choosing a default:

* `λ = 0.5`, the paper's stated default, gives few cones here: 0 to 8 on the meshes with
  boundary (452 on `camelshotnoise`, one at every spike). On nearly developable disks (`bump`, `face`) it gives none, and cones appear only
  below `λ ≈ 0.1` and `0.03`.
* On a closed genus-0 mesh the optimum does not depend on `λ` above a mesh-dependent value
  (`axes` 0.045, `bunny` 0.016, `chair` 0.003, `brain` 0.002): all cones are then positive, the
  cone cost is `4πλ` whatever the configuration, and the problem is the search for the best
  positive measure. That optimum is not sparse: it comes in clusters (`bunny` 254 cones in 12
  clusters, `brain` 123, `chair` 73) and needs merging, with a cluster limit above the default
  8 (`bunny`: 228 cones are left with the limit, 15 without).
* The distortion-minimizing angles for a set of positions do not depend on `λ` at all.

## Performance baseline

`cargo run --release -p ddg-optimal-cones --example baseline -- --merge --repeat 3`: every
dataset mesh at `λ = 0.5`, default options, one thread, fastest of three runs. AMD Ryzen
Threadripper 3970X; wall-clock seconds inside `solve`.

| mesh | triangles | total | numeric factorization | triangular solves | rest | Newton steps I / II | factorizations | solves | cones | merged |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| airplane | 4 656 | 0.12 | 0.09 | 0.02 | 0.01 | 31 / 1 | 32 | 64 | 6 | 4 |
| axes | 7 466 | 0.39 | 0.33 | 0.04 | 0.02 | 31 / 1 | 32 | 64 | 14 | 4 |
| body | 23 892 | 1.26 | 1.08 | 0.14 | 0.05 | 32 / 1 | 33 | 66 | 4 | 4 |
| brain | 305 782 | 33.3 | 25.8 | 6.3 | 1.2 | 49 / 1 | 50 | 100 | 123 | 55 |
| bump | 6 528 | 0.06 | 0.05 | 0.01 | 0.01 | 11 / 1 | 12 | 24 | 0 | 0 |
| bunny | 28 576 | 1.94 | 1.41 | 0.46 | 0.08 | 43 / 1 | 44 | 88 | 254 | 15 |
| camel | 22 704 | 1.02 | 0.76 | 0.21 | 0.05 | 35 / 1 | 36 | 73 | 2 | 1 |
| camel_noise1 | 22 704 | 1.00 | 0.74 | 0.22 | 0.05 | 34 / 1 | 35 | 70 | 5 | 1 |
| camel_noise2 | 22 704 | 1.03 | 0.75 | 0.22 | 0.05 | 33 / 1 | 34 | 68 | 8 | 6 |
| camelbad | 22 725 | 0.86 | 0.58 | 0.25 | 0.03 | 34 / 1 | 35 | 105 | 2 | 2 |
| camelbadnonuniform | 18 125 | 0.64 | 0.43 | 0.18 | 0.03 | 35 / 1 | 36 | 108 | 1 | 1 |
| camelshotnoise | 22 968 | 1.79 | 1.24 | 0.48 | 0.07 | 56 / 1 | 57 | 148 | 452 | 43 |
| chair | 11 448 | 0.46 | 0.36 | 0.08 | 0.03 | 38 / 1 | 39 | 78 | 73 | 7 |
| cow (5 files, one mesh) | 29 999 | 0.87 – 0.91 | 0.61 – 0.63 | 0.21 – 0.24 | 0.04 | 27 / 1 | 28 | 56 | 2 | 2 |
| face | 27 150 | 0.41 | 0.29 | 0.09 | 0.02 | 11 / 1 | 12 | 24 | 0 | 0 |
| hand | 6 660 | 0.23 | 0.19 | 0.03 | 0.01 | 32 / 1 | 33 | 66 | 6 | 3 |
| mannequin | 12 195 | 0.39 | 0.30 | 0.07 | 0.03 | 30 / 1 | 31 | 62 | 4 | 1 |
| rhino | 40 304 | 1.69 | 1.16 | 0.43 | 0.10 | 32 / 1 | 33 | 66 | 2 | 1 |
| tree | 33 592 | 1.39 | 0.96 | 0.39 | 0.05 | 25 / 1 | 26 | 62 | 4 | 4 |

("rest" includes symbolic analysis, assembly, line searches, certificate, extraction. "merged":
all clusters, no size limit.) All 23 solves are certified optima reached without a safeguard
round; state residuals are `10⁻¹⁵` to `2·10⁻¹¹`.

Where the time goes: two thirds to 86 % is numeric factorization. Every Newton step
refactorizes the `2n × 2n` saddle-point matrix (about 0.5 s on `brain`), and there are 12 to 57
steps, all but one in Stage I: 3 to 8 per regularized problem for the first six or seven values
of `γ`, usually one for the last ones. The paper
reports 20 – 25 s for 100 – 150k triangles with 50 – 100 linear solves on a 2018 laptop; the
dataset has no mesh of that size, and on this faster machine 306k triangles take 33 s, so this
reference implementation is in the same range, not far below it. The Laplacian never changes
and could be factorized once; that path is not implemented.

## Formula tags

Functions that implement a formula carry `[formula: <id>]`. Ids that `lean/proofs` proves:
`soft-threshold-fixed-point`, `soft-threshold-active-sets`, `abs-subdifferential-bounded`,
`cone-fixed-point-bounded`, `regularized-cone-update`, `optimality-system-sufficient`,
`duality-gap-certificate`, `newton-stopping-rule`, `closed-mesh-cone-sum`,
`closed-mesh-optimal-constant`, `lambda-area-normalization`, `cotan-laplacian`,
`discrete-yamabe-linearization`. Ids introduced here, with their statements (not yet in
`lean/proofs` or its `CITED.md`, so `scripts/check-formulas.sh` lists them as missing):

| id | function | statement |
|---|---|---|
| `duality-gap-inexact` | `certify` | For any `u, φ, µ` vanishing on the boundary, with `r = L u − µ + Ω`, `φ̂` in the domain of every `g_i*`, `ũ_i = −(L φ̂)_i/w_i` where `w_i > 0`, `J = ½ Σ w_i u_i² + Σ t_i\|µ_i\|` and `D = −½ Σ_{w_i>0} (L φ̂)_i²/w_i + Σ φ̂_i Ω_i − Σ g_i*(φ̂_i)`: `J − D = ½ Σ_{w_i>0} w_i (u_i − ũ_i)² + Σ (t_i\|µ_i\| + g_i*(φ̂_i) − φ̂_i µ_i) + Σ_{w_i=0} (L φ̂)_i u_i − Σ φ̂_i r_i`. (At a vertex with fixed angle `a`: `t_i\|µ_i\|` is `0` and `g_i*(p) = a p`.) |
| `gauss-bonnet-adjoint-shift` | `gauss_bonnet_shift` | For `k_i > 0`, `h(c) = Σ_i clamp(k_i P_{t_i}(x_i + c), lo_i, hi_i)` is continuous, nondecreasing and piecewise linear, with breakpoints among `κ − x_i` for the kinks `κ ∈ {t_i, t_i + hi_i/k_i, −t_i, −t_i + lo_i/k_i}`. The function returns the `c` of smallest magnitude with `h(c) = target`, or, if the target is not attained, the `c` of smallest magnitude at which `h` is nearest to it. |
| `adjoint-constant-interval` | `adjoint_constant_interval` | With `∂g_i(µ_i) = [a_i, b_i]`: `{c : φ_i + c ∈ ∂g_i(µ_i) for all i} = [max_i (a_i − φ_i), min_i (b_i − φ_i)]`. |
| `regularized-dual-objective` | `DualSlope` | For `w > 0`, stiffness `k_i > 0`, centre `µ̄`: `Φ(φ) = ½ Σ w_i u_i² − Σ Ω_i φ_i + Σ e_i(φ_i)` with `W u = −A φ` and `e_i'(p) = clamp(k_i P_{t_i}(p + µ̄_i/k_i), lo_i, hi_i)` is convex and continuously differentiable, `∇Φ = −(A u − µ + Ω)` with `µ_i = e_i'(φ_i)`, and its minimizers are the adjoint states of the problem `minimize J(u, µ) + Σ (µ_i − µ̄_i)²/(2 k_i)`. Along `(du, dφ)` with `W du = −A dφ`: `d/dα Φ = Σ w_i du_i (u_i + α du_i) − Σ Ω_i dφ_i + Σ e_i'(φ_i + α dφ_i) dφ_i`, nondecreasing and piecewise linear in `α`. |
| `proximal-fixed-point` | `solve` | If `(u, µ)` minimizes `J(v, ν) + Σ (ν_i − µ̄_i)²/(2 k_i)` over the feasible pairs and `µ = µ̄`, then `(u, µ)` minimizes `J`. |
| `saddle-point-newton-system` | `SaddleSolver` | The linear system of a choice of pieces is `W u + A φ = 0` together with, per vertex, `(A u)_i = µ_i − Ω_i` (angle fixed), `(A u)_i − k_i (φ_i − o_i) = b_i − Ω_i` (angle `b_i + k_i(φ_i − o_i)`), or `φ_i = p_i` (angle free); in the unknowns `(u, φ)` its matrix is the symmetric `[[W, A], [A, −K]]` with the rows and columns of the pinned `φ_i` replaced by the identity. |
| `componentwise-backward-error` | `backward_error` | `max_i \|b − K x\|_i / (\|K\| \|x\| + \|b\|)_i` is the smallest `ε` such that `(K + δK) x = b + δb` for some `\|δK\| ≤ ε\|K\|`, `\|δb\| ≤ ε\|b\|` entrywise (Oettli and Prager 1964). |
| `cone-cluster-frechet-mean` | `frechet_order` | For a set of cone vertices connected through edges between its members, with angles `µ_k`: the member `v` that minimizes `Σ_k \|µ_k\| dist(v, c_k)²`, `dist` the length of the shortest path along edges with both ends in the set; ties go to the smaller index. |
| `feature-size-weights` | `feature_size_weights` | `wE_i = r_i`, `wR_i = 1/r_i` with `r_i = \|Ω_i\|/m_i + ε` ("as described") or its reciprocal ("as printed"), `m` area-normalized (paper Sec. 8.1). |

## Tests

```
cargo test --release -p ddg-optimal-cones            # everything; the dataset tests need --release
cargo test --release -p ddg-optimal-cones --test dataset -- --nocapture
cargo run  --release -p ddg-optimal-cones --example published -- refit     # or facts, adjoint, sweep, iterate
cargo run  --release -p ddg-optimal-cones --example baseline -- --merge --repeat 3
```

| file | what it checks |
|---|---|
| unit tests in `src/` | the scalar formulas (pieces, proximal map, subdifferential, conjugate), the adjoint shift |
| `tests/convex_reference.rs` | agreement with the generic convex solver Clarabel on random problems with every feature |
| `tests/oracle.rs` | replay of the Lean oracle's golden fixtures and live differential fuzzing against its executable (skipped with a notice if it is not built) |
| `tests/properties.rs` | invariance under scaling and relabeling, monotonicity in `λ`, the guarantee of Sec. 4.5, exact constraints, the certificate of the returned triple |
| `tests/api.rs` | errors, options, tiny meshes, extraction and merging, the certificate identity, non-converged status |
| `tests/dataset.rs` | the published dataset: every mesh solved to a certified optimum, and the findings above |

The dataset is not part of the repository; without it the dataset tests print a notice and pass.

## Not implemented

* Optimized boundary values, boundary curvature penalties, convex boundaries (paper Sec. 8.3).
  Half of the published files used the first. In the terms of this crate it is the same
  problem with the boundary vertices as unknowns that carry their geodesic curvature and are
  free cone sites; `tests/dataset.rs` evaluates it on the double of the mesh.
* A function that re-optimizes the angles of given cones; it is a call of `solve` with
  `cone_costs` of `0` and `+∞`.
* Multiresolution (Sec. 6.3), and a solver that factorizes the Laplacian once.
* Disconnected meshes, non-manifold meshes, intrinsic Delaunay refinement.
