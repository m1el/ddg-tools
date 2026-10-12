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
| `ddg-oracle selftest [N]` | internal checks, with a sweep over `N` random cases (default 300) |
| `ddg-oracle fixtures write DIR` / `check DIR` | regenerate / verify the golden fixtures |
| `ddg-oracle bench` | throughput |

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
  vertex; all three are exactly `0` at boundary vertices. `µ_I` is exactly `0` where there is
  no cone and exactly `LO`/`HI`/the prescribed angle where a bound is attained.
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

Plain semismooth Newton takes the piece each vertex currently lies on. That iteration cycles on
ordinary inputs: with bounds, `clamp((m/γ) P(φ))` jumps from `0` to `hi` over an interval of
width `γ·hi/m`, and vertices flip between "no cone" and "cone at the bound" forever; without
bounds, a vertex with a large angle of the wrong sign flips between `pos` and `neg`. Three
safeguards, in this order:

1. **Adjacent pieces only** (always on). A vertex moves from the piece of the previous step to
   a neighbouring piece of the graph above, never further: not from `zero` to `upper`, not from
   `pos` to `neg`. The next piece equals the previous one exactly when the point lies on it, so
   "pieces unchanged" still means "solved". Whenever no vertex wants to skip a piece this is the
   paper's iteration, and with this rule the eliminated form (Eq. 22) and the form with `µ` kept
   generate the same pieces and the same linear systems.
2. **Single change** (after 60 steps for a regularized problem, 100 for the exact one). The
   problem is restarted with only one vertex changing its piece per step, the one whose piece is
   most wrong; at most 400 such steps.
3. **Singular exact system.** If the Stage II matrix is singular, the minimizer is not unique
   (possible only if some `wE_i = 0`). The oracle then solves the regularized problems for
   `γ = 10⁻¹¹, 10⁻¹², 10⁻¹³` and returns the last one (`STOP = regularized`); these converge
   to the minimizer with the smallest `Σ µ_i²/m_i`.

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
`|GAP| ≤ 10⁻⁹ · max(1, |objective|)`. On the fixtures the residuals are below `2·10⁻¹⁴`.

## How to compare against it

* **`u`, `µ`, objective**: unique whenever all `wE_i > 0` at interior vertices; compare with an
  absolute tolerance (the fixtures reproduce to `10⁻⁹` across platforms; the oracle's own error
  is typically `10⁻¹³`). Compare `µ` numerically, not as a set of cone vertices: at a
  degenerate vertex (`µ_i = 0` and `|φ_i| = t_i`) one correct solver may report `0` and
  another `10⁻¹⁷`.
* **`φ`**: unique on meshes with boundary. On closed meshes compare modulo a constant in
  `[phi_shift LO, phi_shift HI]`; the interval is `[0, 0]` up to rounding when some vertex is
  active.
* **Some `wE_i = 0`**: only the objective and the products `wE_i u_i` are unique. The oracle
  still returns a certified minimizer (`STOP = regularized` when it noticed the
  non-uniqueness).
* **Iteration counts** are informational.
* **Operators**: angles, cotangents and areas are functions of the *rounded edge lengths*
  (`CLAUDE.md`: "corner angles come from edge lengths"). For well-shaped triangles the results
  are accurate to a few ulp. For nearly degenerate triangles they are ill-conditioned functions
  of the lengths, so two correct implementations may differ by about
  `ε · (longest edge / shortest altitude)` in a cotangent and by `ε / sin(angle)` in an angle
  close to `π`.

## Code map

Scalar-generic definitions (`Ddg.Scalar`: arithmetic, order with decidability, `max`, `min`,
`sqrt`, `atan`, `π`, numerals; instance at `Float`). `lean/proofs` imports
`DdgOracle.Formulas` and proves theorems about these definitions at `ℝ`.

| module | definitions |
|---|---|
| `DdgOracle.Scalar` | `Scalar`, `Scalar.abs`, the `Float` instance |
| `DdgOracle.Formulas.Sum` | `sumTo`, `sumList` |
| `DdgOracle.Formulas.Triangle` | `Point`, `dist`, `cornerAngle`, `cornerCot`, `triArea` |
| `DdgOracle.Formulas.Mesh` | `Face` (`angleAt`, `cotOpposite`, `area`, …), `angleSum`, `angleDefect`, `geodesicCurvature`, `laplacianOffDiag`, `laplacianDiag`, `laplacian`, `areaAt`, `totalArea`, `massOf`, `lumpedMass` |
| `DdgOracle.Formulas.Threshold` | `softThreshold` (Eq. 18), `clamp`, `boundedThreshold`, `regularizedMeasure`, `ActivePos`/`ActiveNeg` (Eq. 24), `ConeState` (`coneState`, `value`, `indicator`, `canon`, `toward`), `nextState`, `InSubdiff`, `subdiffViolation`, `subdiffViolationTol`, `boundViolation` |
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

`ddg-oracle selftest` (about 15 s) checks, without any external reference:

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
thread, measured with `ddg-oracle bench`: about 30–220 cases/s at 50 vertices and 1.2–6 cases/s
at 200 vertices; closed meshes are the slow end because all their vertices are unknowns.
`operators` alone: about 800 cases/s at 50 vertices, 150 cases/s at 200. The process is
single-threaded; a fuzzer can run several.

## Limitations

* Dense: `O(n³)` per Newton step, at most 2000 vertices; meant for a few hundred.
* Connected meshes only.
* Positions only (no intrinsic edge-length input), although all formulas are intrinsic.
* Boundary conditions are `u = 0`; the paper's boundary optimization (Sec. 8.3) is not
  implemented.
* Cone merging (Sec. 6.1) and multiresolution (Sec. 6.3) are not part of the canonical problem
  and are not implemented.
* No global convergence proof for the Newton iteration; this is why every result is certified
  and uncertified results are not printed.
