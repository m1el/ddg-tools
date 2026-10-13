# ddg-oracle

A small, dense, dependency-free Lean 4 implementation of

> Yousuf Soliman, Dejan Slepčev, Keenan Crane.
> *Optimal Cone Singularities for Conformal Flattening.* ACM Trans. Graph. 37(4), 2018.

It exists to be compared against: the Rust crates of this repository are fuzzed against it, and
`lean/proofs` proves theorems about the formulas it runs. It is written for meshes of up to a
few hundred vertices and for being obviously right, not for speed.

It solves the **canonical discrete problem** of the repository's `CLAUDE.md` and nothing else:

```
minimize     ½ Σ_i wE_i m_i u_i²  +  λ Σ_i wR_i |µ_i|
subject to   (L u)_i = µ_i − Ω_i  (i interior),   u_i = 0  (i boundary),   lo_i ≤ µ_i ≤ hi_i
```

A second part of the package, [the cut flow](#cut-flow-second-paper) (`ddg-oracle cuts`), implements one
step of *Variational Surface Cutting* (Sharp, Crane 2018) from `crates/ddg-variational-cuts/SPEC.md`;
everything below up to that section is about the first paper.

Contents: [build and run](#build-and-run) · [text format](#text-format) ·
[conventions](#conventions) · [algorithm](#algorithm) · [certificate](#certificate) ·
[how to compare against it](#how-to-compare-against-it) · [code map](#code-map) ·
[tests, fixtures, throughput](#tests-fixtures-throughput) · [limitations](#limitations)

## Build and run

```
cd lean/oracle
lake build                      # Lean 4.34.1 (lean-toolchain), no dependencies
.lake/build/bin/ddg-oracle selftest
```

| command | |
|---|---|
| `ddg-oracle operators [IN [OUT]]` | curvatures, masses, cotan Laplacian of every case in `IN` |
| `ddg-oracle solve [IN [OUT]]` | solve the canonical problem for every case in `IN` |
| `ddg-oracle batch [IN [OUT]]` | the same, each case naming its task in a `task` line |
| `ddg-oracle selftest [N [SIZE]]` | internal checks, with a sweep over `N` random cases (default 300) of up to about `SIZE` vertices (default 41) |
| `ddg-oracle fixtures write DIR` / `check DIR` | regenerate / verify the golden fixtures |
| `ddg-oracle bench` | throughput |
| `ddg-oracle cuts [IN [OUT]]` | one step of the cut flow, [see below](#cut-flow-second-paper) |
| `ddg-oracle cuts-selftest` | checks of the cut-flow oracle |
| `ddg-oracle cuts-fixtures write DIR` / `check DIR` | regenerate / verify `fixtures/cuts` |

`IN`/`OUT` default to `-` (standard input/output). The three processing commands are all batch
commands: they read any number of cases from one stream, write one output block per case in the
same order, and **flush after every block**, so a fuzzer can keep one process open, write a case,
and read lines until `end`. A malformed case produces an error block and never disturbs the
cases after it. The exit code is 0 unless the command line is wrong or I/O fails; the outcome of
each case is in its `status` line.

## Text format

Version 1. Line oriented, ASCII. A line is a tag followed by tokens separated by spaces or tabs;
lines end with LF (a preceding CR is ignored). Blank lines and lines whose first token starts
with `#` are ignored everywhere.

* **Integers** are decimal.
* **Floats are never written in decimal.** A float is the IEEE-754 binary64 bit pattern as
  exactly 16 hexadecimal digits, most significant first: `3ff0000000000000` is `1.0`,
  `7ff0000000000000` is `+∞`, `fff0000000000000` is `−∞`, `8000000000000000` is `−0.0`.
  Output uses lower case; input accepts both cases. In Rust:
  `format!("{:016x}", x.to_bits())` and `f64::from_bits(u64::from_str_radix(tok, 16)?)`.
* **Vertex and face indices are 0-based**, in the order of the `v` and `f` lines.

### Input

```
case NAME                     NAME: any token without whitespace, not starting with '#'
[task operators|solve]        optional; required by `batch`; must agree with the command otherwise
mesh NV NF
v X Y Z                       × NV   vertex positions
f A B C                       × NF   triangles, consistently oriented
[lambda L                     the parameter section: required by `solve`, optional (and then
 w WE WR LO HI                × NV   only syntax-checked) for `operators`
 cones K
 c VERTEX ANGLE]              × K    prescribed cones
end
```

Parameter section, one `w` line per vertex in vertex order:

| token | meaning | allowed values |
|---|---|---|
| `L` | `λ` | finite, `≥ 0` |
| `WE` | distortion weight `wE_i` | finite, `≥ 0` |
| `WR` | cone weight `wR_i` | `[0, +∞]`; `0` is a free cone site, `+∞` (`7ff0000000000000`) forbids a cone |
| `LO`, `HI` | bounds on the cone angle | `−∞ ≤ LO ≤ 0 ≤ HI ≤ +∞`; unbounded is `fff0000000000000` / `7ff0000000000000` |
| `c VERTEX ANGLE` | prescribes `µ_VERTEX = ANGLE` at no cost | interior vertex, finite angle, each vertex at most once |

A `w` line must be present for every vertex. Its values are validated at every vertex but have
an effect only at interior vertices. A prescribed cone overrides `WR`, `LO`, `HI` of its vertex.

The mesh must be a connected, manifold, consistently oriented triangle mesh with at most 2000
vertices, with or without boundary, with finite positions and no degenerate triangle (the three
edge lengths, as computed in floating point, must satisfy the strict triangle inequality).

### Output of `operators`

```
case NAME
task operators
status ok
mesh NV NF NE NB CHI          NE edges, NB boundary vertices, CHI = NV − NE + NF
area A                        total area
o I B CURV MASS               × NV   B = 1 for a boundary vertex, else 0
t F AREA ANG_A ANG_B ANG_C COT_A COT_B COT_C    × NF
laplacian NNZ
l I J VALUE                   × NNZ  entries L_IJ with I ≤ J, sorted by (I, J)
end
```

* `CURV` is the angle defect `Ω_I = 2π − Σ angles` at an interior vertex and the geodesic
  curvature `k_I = π − Σ angles` at a boundary vertex. `Σ CURV = 2π·CHI`.
* `MASS` is the lumped mass `m_I = (⅓ Σ area of the triangles at I) / A`; `Σ MASS = 1`.
* `t` lines give, for face `F = (A, B, C)`, its area and the interior angle and its cotangent at
  each of its three corners, in the order of the face's vertices.
* `l` lines list the upper triangle of the symmetric positive semidefinite cotan Laplacian of
  **all** vertices (boundary rows and columns included): every diagonal entry and one entry per
  edge. `L_IJ = −½ Σ cot(angle opposite edge IJ)`, `L_II = −Σ_{J≠I} L_IJ`.

### Output of `solve`

```
case NAME
task solve
status converged
mesh NV NF NE NB CHI
interior N                    number of interior vertices
closed 0|1                    1 if the mesh has no boundary
iterations S1 S2 STOP         Newton steps of Stage I (all γ together) and of Stage II; see below
objective TOTAL DISTORTION CONE_COST
residuals PRIMAL ADJOINT INCLUSION BOUNDS GAP
gauss_bonnet SUM_MU SUM_OMEGA
phi_shift LO HI
s I U PHI MU                  × NV
end
```

* `s` lines: the log scale factor `u_I`, the adjoint `φ_I` and the cone angle `µ_I` of every
  vertex; all three are exactly `0` at boundary vertices. Where the cone equation fixes `µ_I`
  (no cone, a bound attained, a prescribed or forbidden cone) the value is exactly `0`, `LO`,
  `HI` or the prescribed angle, not a rounded version of it; where it fixes `φ_I` (a cone
  strictly inside its bounds) `φ_I` is exactly `±λ·wR_I`.
* `objective`: `TOTAL = DISTORTION + CONE_COST`, `DISTORTION = ½ Σ wE m u²`,
  `CONE_COST = λ Σ wR |µ|` over the vertices that are neither forbidden nor prescribed.
* `residuals`: the [certificate](#certificate), all computed from the `s` lines.
* `gauss_bonnet`: `Σ µ_I` and `Σ Ω_I` over interior vertices; equal on a closed mesh.
* `phi_shift`: on a closed mesh, `φ + c` is an equally valid adjoint for every constant
  `c ∈ [LO, HI]` (see [closed meshes](#closed-meshes)); `0 0` on a mesh with boundary.
* `STOP` says how Stage II ended: `residual` (all residuals at rounding level), `active_sets`
  (active sets unchanged), `regularized` (see [safeguards](#convergence-and-safeguards)),
  `cap`, `singular`.

### Status

| `status` line | meaning | rest of the block |
|---|---|---|
| `status ok` | `operators` succeeded | as above |
| `status converged` | the certificate holds | as above |
| `status not_converged MESSAGE` | no certified solution was found | as above **without `s` lines** |
| `status infeasible MESSAGE` | closed mesh whose bounds contradict Gauss–Bonnet | nothing |
| `status error MESSAGE` | malformed case, invalid mesh or invalid parameters | nothing |

`MESSAGE` is free text to the end of the line and says what is wrong and where, e.g.
`directed edge (0, 1) occurs in faces 0 and 2: the mesh is non-manifold or not consistently oriented`.
A block always starts with `case` and ends with `end`. A case that could not even be framed is
reported under the name `?`.

## Conventions

Those of `CLAUDE.md`. In particular:

* `L` is positive semidefinite, `u` is the **true** log scale factor, `L u = µ − Ω`.
* A positive `µ_i` is a cone tip: `µ_i = 2π − (angle sum of the flattened cone metric at i)`.
* Masses are normalized by the total area, so `λ` is dimensionless and the result is invariant
  under uniform scaling.
* `φ` is the multiplier of the constraint `L u − µ + Ω = 0` with the sign for which
  `φ_i ∈ ∂g_i(µ_i)`, `g_i(µ) = λ wR_i |µ| + indicator of [lo_i, hi_i]`. Then **`L φ = −wE m u`**.
  Without bounds: `|φ_i| ≤ λ wR_i` everywhere and `φ_i = λ wR_i · sign µ_i` at cones.
  This `φ` is the paper's adjoint; the paper's `u` is the negative of ours.

## Algorithm

The paper's two-stage method (Fig. 12): Stage I solves the regularized problems for
`γ = 1, 10⁻¹, …, 10⁻¹⁰`, each started from the previous solution; Stage II solves the exact
problem started from the last regularized solution. Every problem is solved by semismooth
Newton steps with dense LU factorization with partial pivoting. What follows is everything the
paper leaves open or gets inconsistent, and what the oracle does about it.

### Signs and the optimality system

With `W = diag(wE_i m_i)`, `t_i = λ wR_i` and the Lagrangian
`½ uᵀWu + Σ g_i(µ_i) + φᵀ(L u − µ + Ω)`, stationarity in `u`, `µ` and `φ` gives

```
L u = µ − Ω          L φ = −W u          φ_i ∈ ∂g_i(µ_i)
```

Nothing is copied from the paper's Eq. 20 or App. B (which disagree with each other in the sign
of `µ`): with the paper's `u_paper = −u` the first two equations are the paper's
`L u = Ω − µ`, `L φ = M u` (Eq. 17), so Eq. 20 (`L u − Ω + µ`) is right and App. B
(`L u − Ω − µ`) is a misprint. The weights enter exactly as `CLAUDE.md` says: `wE` multiplies the
mass in the adjoint equation, `wR` scales the threshold per vertex.

### Bounds, forbidden and prescribed cones

In one dimension the proximal map of `f + indicator of [lo, hi]` is the proximal map of `f`
followed by the projection to `[lo, hi]`. Hence

```
φ_i ∈ ∂g_i(µ_i)   ⟺   µ_i = clamp(P_{t_i}(φ_i + µ_i), lo_i, hi_i),     P_t(x) = max(0, x−t) + min(0, x+t)
```

which is the paper's Eq. 19 with the projection added. A forbidden cone (`wR = +∞`) is
`t = 0, lo = hi = 0`; a prescribed cone `c` is `t = 0, lo = hi = c`: both have a constant
proximal map and no cost. The graph of `∂g_i` has up to five pieces, ordered

```
lower (µ = lo, φ ≤ −t)  –  neg (lo < µ < 0, φ = −t)  –  zero (µ = 0, |φ| ≤ t)  –  pos (0 < µ < hi, φ = t)  –  upper (µ = hi, φ ≥ t)
```

The **active sets** are `pos ∪ neg` (the paper's `A⁺ ∪ A⁻`, Eq. 24); on them the Newton system
has the row `φ_i = ±t_i`, elsewhere the row `µ_i = 0`, `lo_i` or `hi_i`. At a free site
(`t = 0`) the pieces `pos` and `neg` are one line and `zero` is a point, so a free site is
always active (the paper's strict inequalities would make it inactive at `φ + µ = 0`).

### The regularized problems (Stage I)

The paper's `µ = P_λ(φ)/γ` (Eq. 21) equates an integrated quantity with a pointwise one. The
oracle regularizes with the squared `L²` norm of the *density* `µ_i / m_i`:

```
minimize  ½ uᵀWu + Σ g_i(µ_i) + (γ/2) Σ µ_i² / m_i        ⟹        φ_i − (γ/m_i) µ_i ∈ ∂g_i(µ_i)
```

Solved for `µ_i` this is `µ_i = clamp((m_i/γ) P_{t_i}(φ_i))`: the paper's formula with the mass
factor. The oracle does not eliminate `µ` like the paper's `F_γ` (Eq. 22) but keeps it as an
unknown and writes the same condition as the proximal fixed point with step `c_i`:

```
µ_i = clamp(c_i · P_{t_i}(φ_i + µ_i), lo_i, hi_i),        c_i = m_i / (m_i + γ) ∈ (0, 1]
```

(on the piece `pos` both say `µ_i = (m_i/γ)(φ_i − t_i)`). For `γ = 0` this is the exact
condition, so **Stage II is the same code with `γ = 0`**, and `µ` is not recomputed from `φ`
between the stages (the paper's `µ ← P_λ(φ)` in Fig. 12 would lose the factor `1/γ`).
The eliminated form is also numerically poor for small `γ`: it recovers `µ` from `φ − t`, a
difference of order `γ`.

### Convergence and safeguards

A Newton step chooses a piece for every vertex, solves the linear optimality system of these
pieces exactly, and adds the step. It stops when all residuals are at rounding level
(`|r| ≤ 10⁻¹³ ·` the magnitude of the terms of `r`) or when the pieces did not change
(the paper's two criteria, Sec. 5.3.3). In Stage II up to two more steps are then taken with the
same LU factorization (iterative refinement), and the values that the cone equation fixes
(`µ_i ∈ {0, lo_i, hi_i}` off the active sets, `φ_i = ±t_i` on them) are assigned exactly.

Plain semismooth Newton takes the piece each vertex currently lies on. That iteration was
observed to cycle on ordinary random inputs in both forms of the regularized problem: in the
eliminated form with bounds, `clamp((m/γ) P(φ))` jumps from `0` to `hi` over an interval of
width `γ·hi/m`, and vertices flip between "no cone" and "cone at the bound" forever; in the form
with `µ` kept, a vertex with a large angle of the wrong sign flips between `pos` and `neg`.
Three safeguards, in this order:

1. **Adjacent pieces only** (always on). A vertex moves from the piece of the previous step to
   a neighbouring piece of the graph above, never further: not from `zero` to `upper`, not from
   `pos` to `neg`. The next piece equals the previous one exactly when the point lies on it, so
   "pieces unchanged" still means "solved". With this rule the eliminated form (Eq. 22) and the
   form with `µ` kept choose the same pieces and solve the same linear systems (in exact
   arithmetic), and whenever no vertex would skip a piece the iteration is the paper's.
2. **Single change** (after 60 steps for a regularized problem, 100 for the exact one). The
   problem is restarted with only one vertex changing its piece per step, the one whose piece is
   most wrong; at most 400 such steps.
3. **Singular exact system.** A Stage II matrix with an active vertex can be singular only if
   some `wE_i = 0`; it typically means that the minimizer is not unique. The oracle then solves
   the regularized problems for `γ = 10⁻¹¹, 10⁻¹², 10⁻¹³`, which are strictly convex in `µ`,
   and returns the last one that converged (`STOP = regularized`); as `γ → 0` their minimizers
   tend to the minimizer of the exact problem with the smallest `Σ µ_i²/m_i`.

None of this decides the status: the result is `converged` if and only if the
[certificate](#certificate) holds, whatever path led to it. If it does not hold the status is
`not_converged` and no solution is printed.

### Closed meshes

Not treated in the paper. `L` annihilates constants, so:

* summing the state equation: `Σ µ = Σ Ω = 2πχ` (Gauss–Bonnet) is necessary for feasibility;
  if `Σ lo > 2πχ` or `Σ hi < 2πχ` (tolerance `10⁻⁹`) the status is `infeasible`;
* summing the adjoint equation: `Σ wE m u = 0`, which fixes the additive constant of `u`
  (it is the optimality condition for that constant); at least one `wE_i` must be positive;
* the adjoint equation determines `φ` only up to a constant `c`, the multiplier of the
  Gauss–Bonnet constraint. It is fixed by the inclusion: every active vertex has `φ_i = ±t_i`.

In the Newton system this means: **with at least one active vertex the `3n × 3n` matrix is
regular** (for `wE > 0`) and nothing special is needed; a free site is always active, which is
why it "removes the degeneracy". With **no active vertex** the matrix is singular, and the start
`u = φ = µ = 0` is such a point. Then `c` is determined from the constraint it belongs to:

1. If the cone angles the current pieces prescribe violate Gauss–Bonnet, `φ` is shifted by the
   constant `c` of smallest magnitude such that `Σ_i clamp(c_i P_{t_i}(φ_i + µ_i + c)) = 2πχ`.
   The left side is nondecreasing and piecewise linear in `c`; the equation is solved exactly
   on its breakpoints. This is the proximal step of `g + indicator{Σ µ = 2πχ}`, and for a
   regularized problem the exact minimization of the dual objective along the null direction
   of its Hessian. Afterwards some vertex is active, or the bounds alone satisfy Gauss–Bonnet.
2. In the latter case (for example the cube with `0 ≤ µ ≤ π/2`: eight cones of exactly `π/2`)
   the system is singular but consistent. It is bordered with one column (a slack on the state
   equation, which absorbs a violation of Gauss–Bonnet at rounding level) and one row
   (`Σ m_i · δφ_i = 0`), and then `φ` is shifted into the interval of constants for which the
   inclusion holds at every vertex.

If at the solution no vertex is active, that interval does not shrink to a point and `φ` is
genuinely not unique. The oracle reports the admissible `φ` whose mean `Σ m_i φ_i` is closest
to `0`, and the remaining freedom as `phi_shift LO HI`.

## Certificate

After every solve the optimality conditions of the canonical problem are checked directly on
the returned `(u, φ, µ)`, independently of the Newton code: `L` is applied face by face from
corner cotangents (not through the dense matrix), and the inclusion is evaluated from the
input parameters (not from the solver's active sets).

| residual | definition |
|---|---|
| `PRIMAL` | `max_i |(L u)_i − µ_i + Ω_i|` over interior vertices |
| `ADJOINT` | `max_i |(L φ)_i + wE_i m_i u_i|` over interior vertices |
| `INCLUSION` | `max_i dist(φ_i, ∂g_i(µ̃_i))`, where `µ̃_i` is `µ_i` snapped to `0`, `lo_i`, `hi_i` if within `10⁻⁹` |
| `BOUNDS` | `max_i dist(µ_i, [lo_i, hi_i])`; `|µ_i − c|` for a prescribed, `|µ_i|` for a forbidden cone |
| `GAP` | primal objective − dual objective, `dual = −½ Σ wE m u² + Σ Ω_i φ_i − Σ g_i*(φ_i)` |

`g_i*(φ) = max(0, hi_i (φ − t_i), lo_i (φ + t_i))` (terms with an infinite bound omitted; their
violation is in `INCLUSION`), `c·φ` for a prescribed cone, `0` for a forbidden one. By
Fenchel–Young the gap of an exactly feasible pair is `Σ (g_i(µ_i) + g_i*(φ_i) − µ_i φ_i) ≥ 0`
and vanishes exactly at the optimum.

The status is `converged` iff the first four residuals are `≤ 10⁻⁹` and
`|GAP| ≤ 10⁻⁹ · max(1, |objective|)`. On the fixtures all five are below `2·10⁻¹⁴`, except on
the one returned by the regularized fallback (`3·10⁻¹²` in `INCLUSION`, `2·10⁻¹¹` in `GAP`).

## How to compare against it

* **`u`, `µ`, objective**: unique whenever all `wE_i > 0` at interior vertices; compare with an
  absolute tolerance. As a guide to the oracle's own error: on the tiny problems of `selftest`
  it differs from the brute-force enumeration by less than `10⁻¹³`, and rescaling or relabeling
  a mesh changes its results by about `10⁻¹⁴`; `fixtures check` uses `10⁻⁹`. Compare `µ`
  numerically, not as a set of cone vertices: at a degenerate vertex (`µ_i = 0` and
  `|φ_i| = t_i`), and at a free site whose optimal angle is `0`, one correct solver may report
  `0` and another `10⁻¹⁷`.
* **`φ`**: unique on meshes with boundary. On closed meshes compare modulo a constant in
  `[phi_shift LO, phi_shift HI]`; the interval is `[0, 0]` up to rounding when some vertex is
  active.
* **Some `wE_i = 0`**: only the objective and the products `wE_i u_i` are unique. The oracle
  still returns a certified minimizer (`STOP = regularized` when it noticed the
  non-uniqueness).
* **Iteration counts** are informational.
* **Operators**: angles, cotangents and areas are functions of the *rounded edge lengths*
  (`CLAUDE.md`: "corner angles come from edge lengths"). For well-shaped triangles these are
  well-conditioned. For nearly degenerate triangles they are ill-conditioned functions of the
  lengths, however accurately the formulas are evaluated: expect two correct implementations
  to differ by roughly `ε · (longest edge / shortest altitude)` in a cotangent and by roughly
  `ε / sin(angle)` in an angle close to `π` (an estimate, not a measured bound).

## Code map

Scalar-generic definitions (`Ddg.Scalar`: arithmetic, order with decidability, `max`, `min`,
`sqrt`, `atan`, `π`, numerals; instance at `Float`). `lean/proofs` imports
`DdgOracle.Formulas` and proves theorems about these definitions at `ℝ`.

| module | definitions |
|---|---|
| `DdgOracle.Scalar` | `Scalar`, `Scalar.abs`, the `Float` instance (`floatPi`, `floatOfNat`) |
| `DdgOracle.Formulas.Sum` | `sumTo`, `sumList` |
| `DdgOracle.Formulas.Triangle` | `Point`, `dist`, `cornerAngle`, `cornerCot`, `triArea` |
| `DdgOracle.Formulas.Mesh` | `Face` (`angleAt`, `cotOpposite`, `area`, …), `angleSum`, `angleDefect`, `geodesicCurvature`, `laplacianOffDiag`, `laplacianDiag`, `laplacian`, `areaAt`, `totalArea`, `massOf`, `lumpedMass` |
| `DdgOracle.Formulas.Threshold` | `softThreshold` (Eq. 18), `clamp`, `boundedThreshold`, `regularizedMeasure`, `ActivePos`/`ActiveNeg` (Eq. 24), `ConeState` (`coneState`, `value`, `indicator`, `isActive`, `canon`, `toward`), `nextState`, `InSubdiff`, `subdiffViolation`, `subdiffViolationTol`, `boundViolation` |
| `DdgOracle.Formulas.Optimality` | `Problem`, `applyL`, `stateResidual`, `adjointResidual`, `coneResidual`, `F` and `DF` (Eqs. 20, 26), `regMeasure`, `Freg` and `DFreg` (Eqs. 22, 25, with the mass factor), `coneScale`, `Fimp`, `DFimp`, `Fpiece`, `DFpiece`, `distortion`, `coneCost`, `objective` |

`Float` only:

| module | |
|---|---|
| `DdgOracle.Hex` | hex encoding of floats, tokenizer |
| `DdgOracle.Topology` | mesh validation, adjacency, boundary, Euler characteristic |
| `DdgOracle.Operators` | gathers the generic mesh formulas into arrays and a dense `L` |
| `DdgOracle.Dense` | LU with partial pivoting; Gauss–Jordan with complete pivoting (brute force only) |
| `DdgOracle.Solver` | problem setup, Newton stages, safeguards, closed-mesh treatment |
| `DdgOracle.Certificate` | certificate, status, the complete `solve` |
| `DdgOracle.BruteForce` | enumeration of all piece assignments for tiny problems |
| `DdgOracle.Format` | parser, printers, streaming driver |
| `DdgOracle.Generators`, `Cases`, `Fixtures`, `Selftest`, `Bench` | test meshes, parameter builders, fixtures, checks |

What runs is the generic code: every operator value is produced by a formula of
`Formulas.Mesh`; the Newton matrix is `Problem.DFpiece` (or `Problem.DF` when, in Stage II, every
vertex stays on its piece) and the right-hand side is `Problem.Fpiece` (or `Problem.F`,
`Problem.Fimp`), entry by entry. `Freg`/`DFreg` (the paper's eliminated form) are definitions
for the proofs only; `selftest` checks that they describe the same solutions.

## Tests, fixtures, throughput

`ddg-oracle selftest` (under 10 s) checks, without any external reference:

* hex round trips; triangle formulas against closed forms; the proximal identities
  (`µ = P(φ+µ) ⟺ φ ∈ ∂g(µ)`, pieces, eliminated vs. implicit regularization) on exact grids
  covering every sign and bound case;
* on a fan, grid, disk, annulus, tetrahedron, octahedron, cube, sphere and torus: Gauss–Bonnet,
  symmetry / zero row sums / positive semidefiniteness of `L`, agreement of the dense `L` with
  the face-wise application, `Σ m = 1`; known curvatures; rejection of eight kinds of invalid
  mesh with their messages;
* the solver against a closed form (one interior vertex), against **brute-force enumeration**
  of all assignments of vertices to pieces on 38 tiny problems (4–8 unknown vertices; disks,
  annulus, tetrahedron, octahedron, bipyramid, cube; all features), and against **random
  feasible perturbations** (none may lower the objective) on medium meshes;
* the paper's qualitative claims: a large `λ` on a disk gives no cones; a closed convex surface
  without free site has a cone at every vertex with `µ = Ω`, `u = 0`; one free site removes
  most of them; the cube with `0 ≤ µ ≤ π/2` gets eight cones of `π/2`;
* infeasible bounds, fully prescribed cones, forbidden-everywhere torus;
* invariance under uniform scaling and under relabeling of the vertices;
* the text format: round trip, batch mode, error blocks and resynchronization;
* a sweep over random meshes and parameters with every feature: all must be `converged` or
  `infeasible`.

`fixtures/` holds 57 golden cases (`MANIFEST.tsv` lists them): for each `NAME` the input
`NAME.in` and the outputs `NAME.operators.out`, `NAME.solve.out`. They cover closed sphere-like
(convex and strongly non-convex), torus, disk, annulus; uniform, random and curvature-adapted
weights; free sites; forbidden regions and the visibility set-up of Sec. 8.1; positive-only,
`[−π/2, π/2]`, `[0, π/2]` and tight bounds; prescribed cones; `λ` from `0` to `0.5`; a problem
with a non-unique minimizer (`STOP = regularized`) and one on which the plain active-set update
cycles; and the outcomes `infeasible` and `error`. `ddg-oracle fixtures check fixtures` recomputes them and
compares with tolerance `10⁻⁹` (the C library's `atan` may differ in the last bit between
platforms, so the files are not expected to be bit-identical everywhere).

Throughput of the whole text pipeline (parse, operators, solve, certificate, print) on one
thread of an AMD Threadripper 3970X, measured with `ddg-oracle bench`:

| mesh | 50 vertices | 200 vertices |
|---|---|---|
| closed, sphere-like (all vertices unknown) | 33 cases/s | 1.3 cases/s |
| torus | 78 cases/s | 3.8 cases/s |
| disk (25 / 144 interior vertices) | 230 cases/s | 5.1 cases/s |
| annulus (30 / 160 interior vertices) | 190 cases/s | 6.1 cases/s |
| `operators` only | about 800 cases/s | about 160 cases/s |

The process is single-threaded; a fuzzer can run several.

## Limitations

* Dense: `O(n³)` per Newton step, at most 2000 vertices; meant for a few hundred.
* If some `wE_i = 0` the minimizer need not be unique; the oracle returns one certified
  minimizer and cannot say whether another solver's different one is wrong. On a closed mesh at
  least one `wE_i` must be positive.
* Connected meshes only.
* Positions only (no intrinsic edge-length input), although all formulas are intrinsic.
* Boundary conditions are `u = 0`; the paper's boundary optimization (Sec. 8.3) is not
  implemented.
* Cone merging (Sec. 6.1) and multiresolution (Sec. 6.3) are not part of the canonical problem
  and are not implemented.
* No global convergence proof for the Newton iteration; this is why every result is certified
  and uncertified results are not printed.

# Cut flow (second paper)

> Nicholas Sharp, Keenan Crane. *Variational Surface Cutting.* ACM Trans. Graph. 37(4), 2018.

`ddg-oracle cuts` computes **one step of the flow** of `crates/ddg-variational-cuts/SPEC.md` (two
regions) on a small mesh with dense linear algebra, and prints every intermediate quantity. It
implements SPEC §2, §3 (including the flux recovery of §3.5, D12, and D13), the seven stages of §4 and the
energy of §6, written from SPEC.md alone; it follows SPEC.md as of commit f1d56b5
(clean room: nothing was read from the Rust crate or from the authors' code). Every sentence of
SPEC that was ambiguous is listed under [choices](#choices-where-the-spec-is-ambiguous).
`ddg-oracle cuts-selftest` checks it (below); the golden cases are `fixtures/cuts/*.in|out`.

The closed-form formulas of SPEC §10 are scalar-generic in `DdgOracle.Formulas.Cuts` (face
curvature, edge crossing, barycentric length, sub-triangle area, weighted-Laplacian entries,
length/penalty step entries, gradient densities, feature size, speed, aspect ratio, smoothing,
explicit step, step cap, class mean, local Hopf-Lax update). The Hopf-Lax update needs one square
root, `min`, `max` and `<`, all in `Ddg.Scalar`; no new operation was added to the class. What
runs is that generic code instantiated at `Float`.

## Cut flow: input

One case per `case ... end` block, any number per stream, same lexical conventions as above
(floats are 16 hex digits, indices are decimal and 0-based, `#` comments). The mesh is **intrinsic**:
faces and edge lengths, no positions. The mesh may be disconnected, with or without boundary,
at most 1500 vertices. All lines after the `e` lines are keyed, may come in any order, and each
occurs once (`eta` and `classes` are optional, everything else is required).

```
case NAME
[task cuts]                       optional; must be 'cuts'
mesh NV NF NE
f A B C                           x NF   oriented triangles
e A B LEN                         x NE   every edge once, A < B, any order; LEN > 0 in the caller's unit
phi P_0 ... P_{NV-1}              level set function in the caller's unit (floats)
energy dirichlet
energy hencky                     pointwise b(u) = u^2
energy scaled_a A_0 ... A_{NV-1}  scaled Dirichlet, weight a > 0 per vertex
energy scaled_curv R_EPS C        scaled Dirichlet, a = r^3 of SPEC 3.6 with r_eps = R_EPS, c = C
curvature gaussian                k_f = Omega_f / A_f
curvature uniform C               k_f = C on every face
curvature per_face K_0 ... K_{NF-1}
alpha ALPHA_L                     >= 0
eta E_0 ... E_{NV-1}              optional, each >= 0 (absent = 0)
penalty gradient | paper          SPEC 3.10: (W_V + tau L_w) rule | PenaltyRule::PaperDiagonal
w W0 W1                           >= 0
beta B0 B1
step fixed TAU                    TAU > 0
step capped TAU_MAX MAX_MOTION    tau = min(TAU_MAX, MAX_MOTION h / max|s|)
flux nodal                        d = v / l (the paper)
flux smoothed C                   SPEC 3.5 flux recovery with c = C >= 0 (SPEC default C = 0.5; C = 0 equals nodal)
delta DELTA                       min_crossing_offset, 0 <= DELTA < 1/2 (0: nothing is moved; a crossing on a vertex is an error)
rho0 RHO0                         sliver threshold > 0 (SPEC default 5)
passes N                          sliver-smoothing passes, decimal integer (SPEC default 0 = off)
classes C_0 ... C_{NV-1}          optional, decimal labels; vertices with equal labels form a class
end
```

The mesh is scaled to unit total area (SPEC §2): every length is multiplied by
`s = 1/sqrt(total area)`. **All parameters (`alpha`, `eta`, `a`'s units (none), `R_EPS`, `C`, `delta` is
dimensionless, curvature `uniform`/`per_face`, `step`) refer to the scaled surface, `phi` is in the
caller's unit** and is multiplied by `s` on the way in. `a` is a dimensionless weight. Every output
quantity is in the scaled unit, except where a column says "caller". Region 0 is `phi <= 0`.

Invalid input gives `status error MESSAGE` (bad mesh or edge table, non-finite or out-of-range
parameter, missing or repeated line, wrong number of values, the paper's penalty rule with some but
not all `alpha + eta_i = 0`, a singular system).

## Cut flow: output

`status ok` (a full step), `status no_cut` (no edge is crossed: only the evaluation at the input
is printed, from `scale` to `energy`), or `status error MESSAGE` (then nothing else but `end`).
Lines in this order; `K` is a cut-vertex index, node number `NV + K`; **nodes** are the original
vertices `0..NV-1` followed by the cut vertices.

```
case NAME
task cuts
status ok|no_cut|error MESSAGE
scale S TOTAL_AREA                  S = 1/sqrt(TOTAL_AREA); TOTAL_AREA in the caller's unit
mesh NV NF NE NB CHI                NB boundary vertices, CHI = NV - NE + NF
sum_omega X                         sum_f Omega_f
counts NC NSEG NTRI NMOVED NSLIVERV NSLIVERT NTIES
f F OMEGA KGAUSS KUSED AREA         x NF   Omega_f, Omega_f/A_f, density used (after an override), A_f
v I REGION THETA AREA R A W PHI SIGMA   x NV
x K A B TRAW T MOVED LEN            x NC   crossed edges, sorted by (A, B)
t IDX PARENT REGION N0 N1 N2 L01 L12 L20 AREA   x NTRI
g K1 K2 FACE LEN                    x NSEG
r R NTRI NNODES NCOMP NUNSOLVED_COMP NUNSOLVED_NODES J AREA AREA_SOLVED     x 2
u R NODE SOLVED U P                 for every node of S_R, regions 0 then 1, nodes increasing
b R NODE V ELL D Q E G              for every node of B_S_R
s K SRAW SSMOOTH SLIVER             x NC
tau TAU HMIN SMAX                   (status ok)
energy TOTAL J0 J1 A0 A1 LENGTH     energy at the input phi
n I PHITILDE PHILEN PHISYM REGION2 ADJ T PHINEW PHINEW_CALLER    x NV
redist SWEEPS NADJ UNREACHED
energy_new TOTAL J0 J1 A0 A1 LENGTH energy at PHINEW
cuts_new NC_NEW                     cut vertices of PHINEW
end
```

* `counts`: `NMOVED` crossings whose `t` was moved into `[delta, 1-delta]` (strict: `t < delta` or
  `t > 1-delta`); `NSLIVERV` cut vertices that have an incident cut-mesh triangle (either region) with
  aspect ratio `> rho0`; `NSLIVERT` such triangles; `NTIES` quads in which the squared lengths of the two diagonals agree to a
  relative `1e-9` (the choice `p-c` / `b-q` is then decided by rounding; the fixtures have none).
* `v`: `REGION` of `phi` (0/1), `THETA` = `Theta_i`, `AREA` = `A_i` (lumped, on `M`), `R` the feature
  size `r_i` and `A` the weight `a_i` (both `0` unless the energy is scaled; `R` only for
  `scaled_curv`), `W` = `alpha + eta_i`, `PHI` the scaled input, `SIGMA` the harmonic extension
  (`0` on `no_cut`).
* `x`: the cut vertex lies at `(1-T) x_A + T x_B`; `TRAW = phi_A/(phi_A-phi_B)`, `T` after the clamp,
  `MOVED` 0/1, `LEN` the length of the edge.
* `t`: the triangles of the cut mesh, faces in index order. An uncrossed face gives one triangle
  `(a, b, c)` in its own order. A crossed face with corners `v0 v1 v2` and alone vertex `a = v_k`,
  `b = v_{k+1}`, `c = v_{k+2}` gives `(a, p, q)` (`p` on `ab`, `q` on `ac`), then, if `|pc| <= |bq|`,
  `(p, b, c), (p, c, q)`, otherwise `(p, b, q), (b, c, q)`; `p-c` is chosen iff `|pc|^2 <= |bq|^2`, both squared
  lengths computed directly from the barycentric formula. The cotangents of a sub-triangle come from the
  barycentric inner product (SPEC §3.3), those of an uncut face are the ones of `M`. `Lij` is the length between corner `i` and
  corner `j` (`L01 = |N0 N1|`, `L12`, `L20 = |N2 N0|`), `AREA` the area (`|det| A_f`).
* `g`: cut segments, one per crossed face (the `(a, p, q)` triangle), `K1 = p`, `K2 = q`, in face order.
* `r`: per region: number of triangles, nodes, connected components of `S_R` (vertex-connected),
  components without a boundary vertex, nodes in them, `J`, area of `S_R`, area of the components
  that were solved.
* `u`: `SOLVED` 0 for nodes of an unsolved component (then `U = P = 0`), the scale factor and the
  adjoint (`P = 0` for the Dirichlet energy).
* `b`: `V` the integrated outward normal derivative, `ELL` the lumped boundary length, `D` the recovered
  normal derivative of `u` (SPEC §3.5: the solution of the flux-recovery system at cut vertices, `V/ELL` at
  the other nodes of `B_S` and for `flux nodal`), `Q` the integrated derivative of the adjoint, `E` the
  recovered normal derivative of the adjoint (`Q` and `E` are `0` for Dirichlet), `G` the gradient density
  before `w_R` and `beta_R`, computed from `D`, `E` (SPEC §3.6). `B_S` = cut vertices and the original boundary vertices in `S_R`.
* `s`: `SRAW = w_0 G_0 + beta_0 - w_1 G_1 - beta_1` at cut vertex `K`, `SSMOOTH` after smoothing,
  `SLIVER` 0/1.
* `tau`: the step used, `HMIN` the shortest crossed edge, `SMAX = max_c |SSMOOTH_c|`.
* `energy`: SPEC §6: `TOTAL = w_0 J0 + beta_0 A0 + w_1 J1 + beta_1 A1 + LENGTH`, where `A_R` is the
  area of the *solved* part of `S_R` (see choices) and `LENGTH = sum_segments (w_p + w_q)/2 |pq|`
  with `w` interpolated at the cut vertices with the clamped `T`.
* `n`: `PHITILDE = phi + tau sigma`; `PHILEN` after SPEC §3.10; `PHISYM` after the symmetry
  projection (equal to `PHILEN` without `classes`); `REGION2`, `ADJ` (0/1) and `T` are the regions, the
  cut-adjacent vertices and the Hopf-Lax values of the redistancing (`T = |PHISYM|` on `ADJ`,
  `+inf` = `7ff0000000000000` where no front arrives); `PHINEW` is the new level set function after
  redistancing and the second projection (scaled unit), `PHINEW_CALLER = PHINEW / S`.
* `redist`: number of accepted non-seed vertices, `NADJ` cut-adjacent vertices, `UNREACHED` vertices that kept
  `PHISYM`.
* `energy_new`/`cuts_new`: the same evaluation at `PHINEW`.

### How to compare

Everything is a deterministic function of the input; the oracle and a correct implementation
agree to rounding **provided the discrete decisions agree** (regions, crossings, the diagonal of
every quad, which vertices are slivers, the sign classes in the redistancing). Fixtures avoid
near-ties (`NTIES = 0`, no `phi` within rounding of 0 except where intended). With `flux nodal`, quantities
that divide by the lumped boundary length (`G`, `SRAW`, `SIGMA` and what follows) are ill-conditioned when a
crossing lies very close to a vertex (`ico1-phi-zero-vertex-nodal`: `max|s| = 10^6`); the default recovery
is not. Compare with a relative tolerance (`1e-9` is what `cuts-fixtures check` uses).

## Cut flow: choices where the SPEC is ambiguous

Each item quotes the sentence and states the choice made (these are also in the delivery report).

1. SPEC §3.4 "A component of `S` with no boundary vertex has no solution ... It is reported and left
   out of everything (D9c)". *Components* are the vertex-connected components of the submesh.
   An unsolved component gets `u = p = 0`, no `v, ell, q, G` (it has no `B_S` node), `sigma = 0`, and
   is excluded from `J` (it would be `0` anyway), **and from `Area(S_r)` in the energy**. The `r`
   line gives both `AREA` and `AREA_SOLVED`; `energy` uses `AREA_SOLVED`. (Reading "everything" as
   not including the area would add `beta_r * (AREA - AREA_SOLVED)`.) A component that touches `dM`
   (e.g. a disk component) has a boundary vertex and is solved.
2. SPEC §3.6 "A region with `w_r = 0` needs no solve." The oracle solves every region anyway and
   prints `u, v, G`; the energy uses `w_r`, so nothing changes.
3. SPEC §3.7 (smoothing optional, `passes` default 0; `passes` and `rho0` are explicit input lines). "a cut
   vertex is a sliver vertex if an incident cut-mesh triangle has aspect ratio `> rho0`": any triangle of the cut mesh (either region, after the quad is split) that has the
   cut vertex as a corner. Aspect ratio from the *area `|det| A_f`* and the barycentric lengths.
   "every sliver vertex simultaneously": every pass reads the speeds of the previous pass.
   Neighbours are the other ends of the cut segments; a vertex with fewer than two neighbours is left
   alone.
4. SPEC §4 "Step size ... `max_c |s_c|`": the *smoothed* speed (the one that is extended), `h` is the
   shortest crossed original edge in the scaled unit. If `max|s| = 0` then `tau = tau_max`.
5. SPEC §3.10 "`w ≡ 0` skips the step": all `alpha + eta_i` exactly `0`, for both rules (then
   `phi' = phi~`). The paper rule with `alpha + eta_i = 0` at some but not all vertices is an error.
6. SPEC §3.12 "Vertices in `Adj` keep their values": the values of `phi'` after the symmetry projection.
   "Signs stay": the sign class of a vertex is the region of that same `phi'` (`phi' <= 0` is class 0).
   A vertex of class `c` ignores neighbours of the other class and neighbours whose `T` is still `+inf`.
   Redistancing is the single-pass fast marching of the current SPEC §3.12 (commit 736f580): seeds are
   accepted; repeatedly the unaccepted vertex with the smallest tentative value (ties: lower index) is
   accepted, by a naive O(n^2) selection. `redistanceFixedPoint` (sweeps to the fixed point, the earlier
   definition) is kept only as a self-test cross-check on a mesh without obtuse angles.
7. SPEC §3.3 "Sub-triangles get ... linearly interpolated `a` of their parent": `a` at a cut vertex is
   `(1-t) a_x + t a_y` (clamped `t`) and `a_t` of a sub-triangle is the mean of the three node values.
8. SPEC §3.3 "the quad is split by the shorter of its diagonals ... in floating point `p-c` iff
   `|pc|^2 <= |bq|^2`": both squared lengths are `baryLength2` of the barycentric differences (no root).
   `NTIES` counts relative gaps below `1e-9`.
9. SPEC §6 "`w` interpolated linearly to the cut vertices": with the clamped `t`; each segment once.
10. SPEC §3.2: `NMOVED` counts `t < delta` or `t > 1 - delta` strictly; a `t` exactly equal to
    `delta` is not counted.
11. SPEC §3.1 "A curvature override replaces `k_f`": it replaces the density that drives the Yamabe
    equation (`g`); `OMEGA` and `KGAUSS` in the output stay the Gaussian ones. `uniform`/`per_face` are
    in the scaled unit (1/length^2). The curvature-based weight (SPEC §3.6) always uses the vertex
    curvature `K_i = (2 pi - Theta_i)/A_i` of `M`, never the override.
12. SPEC §3.6 curvature-based weight: `(W_V + c L) r = W_V r0` is solved on the whole of `M` with
    `L` the cotan Laplacian of `M` (scale invariant) and `A_i` the scaled lumped area.
13. SPEC §3.6 pointwise: the Hencky case has `b(0) = 0`, `b'(u) = 2u`; the adjoint right-hand side at a
    boundary node is `-A_i^S b'(0) = -0`.
14. SPEC §3.11: a class mean is the sum in increasing vertex index divided by the class size.
15. SPEC §2 "`E(φ)`" in `energy_new` is the same function at `PHINEW`, not part of the SPEC's one
    step; it is there for the energy-decrease check.

16. SPEC §3.5 flux recovery: one linear system per region over *all* cut vertices of the cut mesh (it
    decouples over the pieces by itself). `kappa = c^2 h_f^2 / |pq|` with `|pq|` the segment length from the
    `(a, p, q)` sub-triangle and `h_f^2` from the scaled edge lengths of the parent face. The same
    system, with the same `ELL`, is solved for `V` and for `Q` (scaled and Hencky energies). For
    `flux nodal` no system is solved and `D = V/ELL` exactly; `smoothed 0` solves a diagonal system and
    agrees with it to rounding (measured: identical).
17. SPEC §3.6 "`G` is needed at cut vertices only": `G` is nevertheless printed at every node of `B_S`, with
    `D = V/ELL`, `E = Q/ELL` at the original boundary vertices.
18. SPEC §3.2 / D13f "with `delta = 0` ... a crossing on a vertex is an error": the oracle raises the error
    when any sub-triangle has non-positive area (this is also the only check of D13a), message
    "sub-triangle of face F has no area". With `delta > 0` crossings are moved and no such error occurs.
19. SPEC §3.1 / D13e: `Omega_f` is evaluated in the cancellation-free form (sum over interior corners of
    `theta (2 pi - Theta)/Theta`); it differs from the first form by rounding only.
20. SPEC §3.4 / D13b: the "identity rows" are realized by excluding the unsolved component from the
    systems; `u = p = 0` and `sigma = 0` there. The report is the `r` line (`NUNSOLVED_COMP`, `NUNSOLVED_NODES`,
    `AREA - AREA_SOLVED`).

## Cut flow: tests and fixtures

`ddg-oracle cuts-selftest` (about 1 s, 37 checks) runs, without external reference:

* icosphere generator (levels 0-3: outward, `chi = 2`, vertex counts);
* `sum_f Omega_f` equals the sum of `2 pi - Theta_i` over interior vertices to `4e-14` on a sphere,
  torus, bumpy disk, annulus, lumpy sphere (`2 pi chi` on the closed ones);
* cut mesh: areas sum to 1 (`2e-15`), `chi(S_0) + chi(S_1) - chi(cut) = chi(M)` on a sphere, torus
  and disk;
* flat disk: `u = 0` (`2e-15`), `v = 0`, `s = 0` (`2e-28`);
* scaled Dirichlet with `a = 1` equals Dirichlet (`J`, `s` to `2e-15`; `p = -2u`, `q = -2v`);
* flat disk with a constant source `K = 5` and a circle of radius `0.3`: Dirichlet `J` is within `3e-4`
  of `pi K^2 rho^4/8`, mean `G` within `1e-3` of `(K rho/2)^2`; Hencky `J` within `1.2e-2` of
  `pi K^2 rho^6/48`, `G` within `6e-4` of `K^2 rho^4/16`; `max |u - exact| = 1e-3`;
* sphere cap (SPEC §6.2) on the icosphere of level 3 (642 vertices), polar angle 1: `J = 0.391` against
  `0.394` (0.55%), `du/dn = 1.9368` against `1.9370` (1e-4), `max|u - exact| = 1e-3`, and the exact
  identity `sum_B v_i = sum_f k_f A_f` to rounding;
* the stage equations on the output (Yamabe, Hencky adjoint, harmonic extension) to `3e-14`;
* finite differences (SPEC §6.3): shifting `phi` by `+-e` changes `E` by `-2 e sum_c s_c l_c` within 0.5%
  for Dirichlet, Hencky with `w = (1, 2)`, `beta = (0.5, 0.1)`, and curvature-scaled with `w_0 = 1.5`;
* redistancing (fast marching): every value equals the update from the vertices accepted before it,
  cut-adjacent values are kept; on a flat grid with a straight cut it reproduces the exact distance
  and equals the fixed point of the scheme (`5e-17`);
* flux recovery on three cuts (icosphere 2 tilted plane, icosphere 1 with `phi = z`, a geodesic cap): `sum l d = sum v` on
  every piece of the cut (`1e-14` relative), `d` within the range of `v/l` of the piece, `smoothed 0` equals
  nodal; maximal `|s|` (Hencky) nodal -> default recovery: `3.8e5 -> 0.9`, `8.7e4 -> 2.9e-3`, `7.6 -> 6.2`;
  `delta = 0` runs when no crossing is on a vertex and is an error otherwise;
* the closed form of the local Hopf-Lax update against a ternary search over lambda (`1e-12`);
* a change of unit (lengths and `phi` times 3) leaves the scaled quantities unchanged;
* energy decrease along a small step; the text format round trip; error blocks.

`fixtures/cuts/` has 26 cases (`MANIFEST.tsv`): `NAME.in` and `NAME.out` (the output of `ddg-oracle cuts
NAME.in`): icosphere levels 1 and 2 with a tilted plane, geodesic caps with each energy, a torus, a grid
disk and an annulus whose cuts meet the boundary, an annulus ring, `eta` with both penalty rules, `w_r`,
`beta_r`, symmetry classes (both projections active), `phi = 0` exactly at vertices, `delta = 0.2` clamps,
most cases use the default recovery (`flux smoothed 0.5`, no sliver smoothing); the paper's procedure (`flux nodal`,
sliver smoothing 1-3 passes) is pinned by `ico2-tilted-dirichlet`, `ico1-delta-clamp`, `ico2-sliver-smoothing`,
`ico1-phi-zero-vertex-nodal`; also `flux smoothed 1.0`, recovery plus sliver passes, `delta = 0`; the fixtures
`ico2-tilted-dirichlet(-recovery)` and `ico1-phi-zero-vertex(-nodal)` show the blow-up (`max |SSMOOTH|` after
one pass: `63 -> 10` and `1.15e6 -> 1.14`), curvature overrides (uniform and per face), a flat disk
and two disjoint spheres where one component is not touched. `ddg-oracle cuts-fixtures check
fixtures/cuts` recomputes and compares (tolerance `1e-9`, as for the first paper).

## Cut flow: code map

| module | |
|---|---|
| `DdgOracle.Formulas.Cuts` | scalar-generic formulas of SPEC §10 (see above) |
| `DdgOracle.Cuts.Case`, `Mat` | input data model; dense matrices with per-triangle cotan assembly and principal-subsystem solves |
| `DdgOracle.Cuts.Base` | validation, scaling, `Theta`, `Omega_f`, vertex areas, `L`, feature size |
| `DdgOracle.Cuts.CutMesh` | regions, crossings, cut mesh |
| `DdgOracle.Cuts.Region` | `u, v, ell, p, q, G, J` of one region |
| `DdgOracle.Cuts.Step` | speed, smoothing, extension, energy, length step, symmetry, redistancing, `runStep` |
| `DdgOracle.Cuts.Format` | parser, printers, stream driver |
| `DdgOracle.Cuts.Gen`, `Selftest`, `Fixtures` | icosphere and test cases, checks, golden fixtures |

`buildTopology` got an optional argument `connected := true`; the cut flow passes `false`.

## Cut flow: limitations

* Two regions only (SPEC §8 is not implemented); `b` is Hencky only (no general `b, b'`).
* Dense: at most 1500 vertices, practical limit a few hundred (LU of the full `n x n` systems).
* One step; `cuts` does not iterate and has no stopping rule.
* Not implemented: the visibility helper `eta_i = nu_i n_i . E_i`, initial-`phi` helpers.
