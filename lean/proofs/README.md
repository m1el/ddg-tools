# lean/proofs — machine-checked derivations of the formulas

Lean 4 / Mathlib proofs of the discrete formulas used by `ddg-tools`:

* the canonical cone problem and its solver, `crates/ddg-optimal-cones` — *Optimal Cone
  Singularities for Conformal Flattening* (Soliman, Slepčev, Crane, ACM TOG 2018; "the paper");
* the flattener, `crates/ddg-bff` — *Boundary First Flattening* (Sawhney, Crane, ACM TOG 2017;
  "BFF");
* the cut flow, `crates/ddg-variational-cuts` — *Variational Surface Cutting* (Sharp, Crane,
  ACM TOG 2018; "the cuts paper"); its contract is `crates/ddg-variational-cuts/SPEC.md` ("SPEC");
* intrinsic triangulations, `crates/ddg-intrinsic` — *Navigating Intrinsic Triangulations*
  (Sharp, Soliman, Crane, ACM TOG 2019; "the signpost paper"); its contract is
  `crates/ddg-intrinsic/SPEC.md` ("the intrinsic SPEC", decisions `N1`–`N17`);
* the developability flow, `crates/ddg-developability` — *Developability of Triangle Meshes*
  (Stein, Grinspun, Crane, ACM TOG 2018); its contract is `crates/ddg-developability/SPEC.md`;
* the operators of `crates/ddg-core` all of them use.

The ground truth is the **canonical discrete problem** and the math conventions of the
repository's `CLAUDE.md`; where a paper is silent or inconsistent the statements here are derived
from that problem (see [Discrepancies](#discrepancies-with-the-papers-as-printed)).

**Scope in one paragraph.** Every theorem is a statement in exact real arithmetic about a formula;
rounding is not covered anywhere, and the Rust code is tied to the statements by reading, by the
`[formula: id]` tags and by differential tests, not mechanically. Results that are used but not
proved here (the smooth theory, convergence of the Newton method, three topological facts, …) are
listed with exact references in [`CITED.md`](CITED.md), which also states precisely what a
run-time certificate of the cone solver does and does not establish. For the cut flow there is no
discrete problem with a certifiable minimizer: the theorems cover the formulas and linear systems
of one flow step and three closed-form solutions of the smooth problem; `CITED.md` says what that
leaves open. For intrinsic triangulations the theorems cover the closed forms of single
operations (unfolding, flip, Delaunay test, signposts, one tracing step, barycentric coordinates,
circumcentre); the data structure and the algorithms built on them are not modelled, and the
global results (existence, uniqueness and termination for the Delaunay triangulation, the
refinement guarantees) are cited. For developability the theorems cover the per-face derivatives (with the
corrections of SPEC D2, D3), the gradients of the active pieces of the six energies, the line
search's acceptance conditions, subdivision and two closed forms of the regular cone; the hinge
theory of the paper, the Voronoi and rotating-calipers arguments and the derivative of a simple
eigenvalue are cited.

## Pins and how to check

| | |
|---|---|
| Lean | `leanprover/lean4:v4.34.1` (`lean-toolchain`) |
| Mathlib | tag `v4.34.1`, revision `d13f23b723b8a846827a245b89c10fc7d3f11612` (`lake-manifest.json`) |
| Oracle | path dependency `../oracle` (package `ddgOracle`, no Mathlib); only `DdgOracle.Formulas` is imported |

```sh
cd lean/proofs
lake exe cache get                      # prebuilt Mathlib; do not compile Mathlib from source
lake build > ../../runs_tmp/lean-proofs-build.log 2>&1
lake env lean DdgProofs/Audit.lean      # prints the axioms each main theorem depends on
../../scripts/check-formulas.sh         # formula ids: Rust <-> Lean theorem / CITED.md; no sorry/axiom
../../scripts/check-formulas.sh --build # the same, plus build and axiom audit of the main theorems
../../scripts/check-formulas.sh --audit-all   # the same, plus the axioms of every declaration (a few seconds more)
```

There is no `sorry` and no `axiom` declaration in this package. `DdgProofs/Audit.lean` runs
`#print axioms` on the main theorems, by name (the count is printed by `check-formulas.sh --build`);
every one reports exactly `[propext, Classical.choice, Quot.sound]`.
`lake env lean DdgProofs/AuditAll.lean` (a few seconds) collects the axioms of *every*
declaration of the package and of the oracle's formula modules and reports any other axiom.

The lakefile sets `autoImplicit = false`; `lake env lean File.lean` on a single file does **not**
apply that option (pass `-DautoImplicit=false`), `lake build` does.

## Layout

| File | Content |
|---|---|
| `DdgProofs/SoftThreshold.lean` | soft threshold, clamp, subdifferential of `t·\|x\|` on an interval (item A.1) |
| `DdgProofs/Problem.lean` | the canonical discrete problem: `Site`, `ConeProblem`, objective, feasibility, optimality system |
| `DdgProofs/Optimality.lean` | sufficiency, quadratic growth, uniqueness, necessity (invertible `L`), guarantee of paper Sec. 4.5 (A.2, A.3) |
| `DdgProofs/Duality.lean` | discrete dual, Fenchel–Young, weak/strong duality, duality-gap certificate (A.4) |
| `DdgProofs/Certificate.lean` | the certificate for inexact data: gap identity for arbitrary points, clipped dual lower bound, repaired upper bound, enclosure, backward statements, counterexamples (A.4′) |
| `DdgProofs/Newton.lean` | residual maps, semismooth differentials, stopping rule, regularized problem (A.5) |
| `DdgProofs/Solver.lean` | the linear systems `crates/ddg-optimal-cones` solves, closed-mesh treatment, adjoint constant, Gauss–Bonnet shift, line-search slope, proximal stage, backward error, cone merging (C) |
| `DdgProofs/Closed.lean` | closed meshes, optimal additive constant, area normalization of `λ` (A.6, A.7) |
| `DdgProofs/Components.lean` | the problem of a mesh with several connected components is the sum of the problems of its components (block diagonal `L`, common `λ`, masses normalized by the total area): minimizers, optimality system and Gauss–Bonnet per component; a component solved alone with its own area normalization and `λ_c = (A/A_c)·λ` (A.8) |
| `DdgProofs/Triangle.lean` | corner angle, cotangent, Heron area, angle sum, scaling, link to Mathlib's Euclidean geometry (B.8) |
| `DdgProofs/Cotan.lean` | cotan formula, cotan Laplacian, lumped mass (B.9, B.10) |
| `DdgProofs/GaussBonnet.lean` | angle defect, boundary curvature, discrete Gauss–Bonnet (B.11) |
| `DdgProofs/Yamabe.lean` | linearized discrete Yamabe equation (B.12) |
| `DdgProofs/Evaluated.lean` | the same formulas in the arrangements `crates/ddg-core` evaluates (Kahan's area, `atan2` corner angle, curvature fields, Laplacian entries) |
| `DdgProofs/MeshProblem.lean` | the canonical problem of a mesh (`TriMesh.coneProblem`) and the translation between vectors on the unknowns and per-vertex vectors that vanish on the boundary |
| `DdgProofs/FactorOnce.lean` | the factor-once method of `crates/ddg-optimal-cones/src/factor_once/`: the reduced (LASSO-type) problem in the cone vector, its Gram matrix, its optimality conditions and their link to `KKT`, the pinned solve of a singular `A` (C) |
| `DdgProofs/BoundaryProblem.lean` | the mesh problem with boundary unknowns (`TriMesh.boundaryProblem`, `Ω̂`): sufficiency of the optimality system, `Σ s = 2πχ`, the free / penalized / convex boundary sites, first-order curvature of the scaled boundary |
| `DdgProofs/Bff/Laplace.lean` | BFF: Poisson solves, Dirichlet-to-Neumann, Neumann problem, exterior-angle sum, first-order Cherrier formula (D) |
| `DdgProofs/Bff/Cut.lean` | BFF: cut meshes, angle sums at cones, Euler characteristic after cutting (D) |
| `DdgProofs/Bff/Curve.lean` | BFF: boundary tangents and curve, closing lengths, exterior-angle normalizations, Hilbert-transform data (D) |
| `DdgProofs/Bff/Quality.lean` | BFF: distortion measures, signed area and angle, segment crossing, similarity fit, mollification (D) |
| `DdgProofs/Bff/Sphere.lean` | BFF map to the sphere: stereographic projection (on the sphere, inverse, conformal and orientation preserving), inversions of the sphere and the derivative behind the centering Jacobian (D) |
| `DdgProofs/Cuts/Crossing.lean` | cut flow: crossing parameter and its clamp, barycentric length and inner product, areas and cotangents of sub-triangles (SPEC §3.2, §3.3) |
| `DdgProofs/Cuts/Curvature.lean` | cut flow: rescaled angles, face curvature `Ω_f`, its sum and Gauss–Bonnet, own and oracle definitions (§3.1) |
| `DdgProofs/Cuts/Weighted.lean` | cut flow: the weighted cotan Laplacian `L_w = Σ_f w̄_f L_f`, the paper's edge rule and a counterexample to its semidefiniteness (§3.6, §3.10, D11) |
| `DdgProofs/Cuts/Dirichlet.lean` | cut flow: Dirichlet solves on an index set split into `I` and `B`: Yamabe system, boundary flux, Dirichlet energy, scaled and pointwise adjoints, gradient densities, harmonic extension (§3.4–§3.6, §3.8) |
| `DdgProofs/Cuts/Step.lean` | cut flow: sign of the flow, explicit step, length step, the paper's penalty rule, symmetry projection, step cap, feature size, sliver smoothing (§1, §3.6–§3.11, §4) |
| `DdgProofs/Cuts/Recovery.lean` | cut flow: recovery of the normal derivative from the integrated fluxes (§3.5, D12) |
| `DdgProofs/Cuts/Redistance.lean` | cut flow: the Hopf–Lax triangle update in closed form is the minimum over the edge (§3.12) |
| `DdgProofs/Cuts/OnCut.lean` | cut flow, default stages: the two systems on the cut (length force, rate), the extension from the cut, `\|∇φ\|` in a crossed face and its mean, the per-step factor of a ripple, the adaptive step (§3.8–§3.10, §4, D14–D16) |
| `DdgProofs/Cuts/Ripple.lean` | cut flow: the roughening rate `ωk − wk²`, the length weight for a cutoff wavelength, the ingredients of the linearization on closed forms, the sphere rates (§1, §4, D18) |
| `DdgProofs/Cuts/Quantile.lean` | cut flow: the weighted quantile, the roughening rate of a state as a quantile, the robust step rule and the local limit per cut point (§1, §4, D22–D25) |
| `DdgProofs/Cuts/VertexLimit.lean` | cut flow: the local limit at the vertices, on the change of the whole step, for two and for `N` regions (§4, §8.4, D26) |
| `DdgProofs/Cuts/Stop.lean` | cut flow: the stop rule of a run: the three cases of the energy criterion, the state of lowest energy, the motion criterion with the predicted motion (§4, D27, D29, D31; `cut-stop-energy`, `cut-run-lowest`, `cut-stop-motion`) |
| `DdgProofs/Cuts/Regions.lean` | cut flow with `N` regions: region of a vertex, crossings between regions, triple points, the projection of the redistancing, rates of the functions, the systems on the boundary of a region (§8, D20, D21) |
| `DdgProofs/Cuts/Analytic.lean` | cut flow: closed-form solutions (flat disk, spherical cap, shrinking circle) and the shape derivatives of the energies on them (§1, §6) |
| `DdgProofs/Intrinsic/Diamond.lean` | intrinsic triangulations: the layout of a triangle and of the two triangles on an edge from their lengths, uniqueness up to a rigid motion, the length of a flipped edge, strict convexity (intrinsic SPEC §4, N5) |
| `DdgProofs/Intrinsic/Delaunay.lean` | intrinsic triangulations: the Delaunay condition and the cotan weight, the link to the cotan Laplacian, what a flip does to the angles and to the condition (§4, §5, N6, N17) |
| `DdgProofs/Intrinsic/Signpost.lean` | intrinsic triangulations: angle normalization, the signpost update around interior and boundary vertices and under a flip (§3, N2) |
| `DdgProofs/Intrinsic/Trace.lean` | intrinsic triangulations: one step of a tracing query: crossing of an edge with the path, choice of the exit edge, forward progress (§4, N7) |
| `DdgProofs/Intrinsic/Barycentric.lean` | intrinsic triangulations: barycentric coordinates, distances of a point to the corners, edge split, circumcentre (§4, N10, N14, N17) |
| `DdgProofs/Intrinsic/Move.lean` | intrinsic triangulations: moving a vertex (lengths, areas of its star), the circumcentre and barycentre times the area, the energy of the optimal Delaunay triangulation of a star, its gradient, the step of Chen and Holst's Eq. 4.13 (§4, N19, N21) |
| `DdgProofs/Developability/*.lean` | developability: vectors and derivatives along curves (`Vec`), per-face quantities and Eqs. 7–9 (`Basic`), invariances, the six energies (`Combinatorial`, `Covariance`, `MaxCovariance`, `Intrinsic`, `Width`), Jacobi rotation (`Eig`), line search (`Optimizer`), subdivision, remeshing classification, regular cone (developability SPEC §1–§5, §10.3) |
| `DdgProofs/Oracle/*.lean` | the real-number instance of the oracle's `Ddg.Scalar` and the theorems restated about the oracle's definitions; `Oracle/Unbounded.lean`: the same with absent bounds |
| `DdgProofs/Audit.lean`, `DdgProofs/AuditAll.lean` | axiom audits (not imported by the library) |

## Conventions

- Namespace `DdgProofs` holds the math-level definitions of this package. Namespace `Ddg` belongs
  to `lean/oracle`; theorems in `namespace Ddg.Problem` / `DdgProofs.Oracle` are about the oracle's
  scalar-generic definitions instantiated at `ℝ` (the same definitions the `Float` executable runs).
- Signs are those of `CLAUDE.md`: `L` is positive semidefinite, `u` is the true log scale factor,
  the state equation is `L u = µ − Ω`. The **optimality system** for `(u, φ, µ)` is

  ```
  L u  = µ − Ω                       (state)
  Lᵀφ  = −(wE ∘ m ∘ u)               (adjoint:  (Lᵀφ)_i = −wE_i m_i u_i)
  φ_i ∈ ∂(λ wR_i |·| + ι_[lo_i,hi_i])(µ_i)      (inclusion)   ⟺   µ_i = clamp_[lo_i,hi_i](P_{λ wR_i}(φ_i + µ_i))
  ```

  It is paper Eq. 17 / System B.2 with `u_paper = −u`; `φ` and `µ` are the same in both conventions.
- `ConeProblem ι`: `ι` is the set of unknowns — the interior vertices on a mesh with boundary
  (`L = L_II`, `u = 0` on the boundary already eliminated, see `laplacian_mulVec_of_zero_off`), or
  all vertices on a closed mesh. On a closed mesh the competitors in `IsMinimizer` are *all*
  solutions of `L v = ν − Ω`, so the additive constant of `u` is optimized.
  `TriMesh.coneProblem M p wE λ site` (`DdgProofs/MeshProblem.lean`) is this problem for a mesh
  `M` with unknowns `p`, with the masses divided by the total area. For per-vertex vectors that
  vanish off the unknowns (as the Rust arrays do; an exact test), its state residual is
  `(L u)_v − µ_v + Ω_v` with the full Laplacian (`coneProblem_stateResidual`), its adjoint product
  is `(L φ)_v` (`coneProblem_adjoint`), feasibility is the constraint set of `CLAUDE.md`
  (`coneProblem_feasible_iff`), and the objective is the per-vertex sum (`coneProblem_objective`);
  the matrix is symmetric, positive definite in the Dirichlet case and annihilates constants on
  a closed mesh (`coneProblem_transpose`, `coneProblem_posDef`, `coneProblem_mulVec_one`).
- `Site` is the per-vertex data: `penalized wR lo hi` (cost `λ·wR·|µ|`, `µ ∈ [lo, hi]`,
  `lo, hi : EReal` with `⊥/⊤` = unbounded, `wR = 0` = free site) or `fixed a` (`µ = a` at no cost:
  `fixed 0` is a forbidden vertex, i.e. `wR = +∞`; `fixed a` a prescribed cone).
- The subdifferential is written out: `IsSubgradientWithinAt f K g x` means `x ∈ K` and
  `f x + g·(y − x) ≤ f y` for all `y ∈ K`, i.e. `g ∈ ∂(f + ι_K)(x)`.
- In `cornerAngle a b c`, `cornerCot a b c` the corner is the one opposite the **first** argument.

## Formula ids and their exact statements

This is the list to take `[formula: id]` tags from. Each id names one formula; a Rust function
carrying the tag is expected to compute or check exactly this. Notation: `t = λ·wR ≥ 0`,
`w_i = wE_i m_i`, `P_t(x) = max(0, x − t) + min(0, x + t)`,
`clamp(x) = max(lo, min(hi, x))` with `lo ≤ 0 ≤ hi`, `J(u, µ) = ½ Σ w_i u_i² + Σ t_i|µ_i|`,
`θ_a` the corner angle opposite the side `a` of a triangle with sides `a, b, c`.
All statements are for exact real arithmetic. Ids listed in [`CITED.md`](CITED.md) are not here.

### Geometry and operators

| id | exact statement | theorems |
|---|---|---|
| `edge-length-euclidean` | `ℓ(p, q) = √((p₀−q₀)² + (p₁−q₁)² + (p₂−q₂)²)` is the Euclidean distance in `ℝ³`. | `dist_eq_sqrt_sum_sq` |
| `triangle-area` | For `a < b + c`, `b < c + a`, `c < a + b`: `A = ¼√((a+b+c)(−a+b+c)(a−b+c)(a+b−c))`, so `16A² = (a+b+c)(−a+b+c)(a−b+c)(a+b−c)`, and `A = ½·b·c·sin θ_a`. With the sides sorted `x ≥ y ≥ z`: the three inequalities hold iff `z − (x − y) > 0`, and `A = ¼√((x+(y+z))(z−(x−y))(z+(x−y))(x+(y−z)))`. | `sixteen_mul_triArea_sq_eq`, `sin_cornerAngle`, `half_mul_sin_angle_eq_triArea`, `kahanArea_of_sorted_rearrangement` |
| `corner-angle` | `θ_a = arccos((b² + c² − a²)/(2bc)) ∈ (0, π)`; equivalently `θ_a = 2·atan2(√((a−b+c)(a+b−c)), √((a+b+c)(−a+b+c)))` (`atan2(y, x)` = argument of `x + iy`), also in Kahan's arrangement with the two adjacent sides sorted. For points of a Euclidean space it is the angle at the vertex opposite `a`. | `cos_cornerAngle`, `cornerAngle_pos`, `cornerAngle_lt_pi`, `cornerAngle_eq_two_arctan`, `two_mul_atan2_sqrt_eq_cornerAngle`, `kahanCornerAngle_eq_cornerAngle`, `angle_eq_cornerAngle` |
| `corner-cot` | `cot θ_a = (b² + c² − a²)/(4A)`; the numerator may be evaluated as `c² + (b − a)(b + a)`. | `cornerCot_eq_cot`, `kahanCornerCot_eq_cornerCot` |
| `triangle-angle-sum` | `θ_a + θ_b + θ_c = π`. | `cornerAngle_add_cornerAngle_add_cornerAngle` |
| `cotan-formula` | For an affine `f` on a triangle with vertex values `f₁, f₂, f₃` and gradient `g`: `A·‖g‖² = ½(cot θ₁ (f₂−f₃)² + cot θ₂ (f₃−f₁)² + cot θ₃ (f₁−f₂)²)`. | `cotan_formula` |
| `cotan-laplacian` | `L_ij = −½ Σ_f cot(corner of f opposite edge ij)` for `i ≠ j` (faces not containing the edge contribute 0), `L_ii = −Σ_{j≠i} L_ij`. Then `L = Lᵀ`, `L·1 = 0`, and `uᵀLu = Σ_f ½(cot_a (u_b−u_c)² + cot_b (u_c−u_a)² + cot_c (u_a−u_b)²)`. | `TriMesh.laplacian_apply_of_ne`, `laplacian_apply_self`, `laplacian_isSymm`, `laplacian_mulVec_one`, `dotProduct_laplacian_mulVec_eq_sum_faceEnergy` |
| `cotan-laplacian-psd` | `uᵀLu ≥ 0` for every valid mesh (no Delaunay condition). | `cotEnergy_nonneg`, `TriMesh.laplacian_posSemidef` |
| `cotan-laplacian-kernel` | Connected valid mesh: `L u = 0 ⟺ u` constant. | `TriMesh.laplacian_mulVec_eq_zero_iff` |
| `cotan-laplacian-interior-posdef` | Connected valid mesh with at least one non-interior vertex: `L_II` is positive definite. If `u = 0` off `I` then `(L u)_i = (L_II u_I)_i` for `i ∈ I`. | `TriMesh.laplacian_interior_posDef`, `laplacian_mulVec_of_zero_off` |
| `lumped-mass` | `m_i = ⅓ Σ_{f ∋ i} A_f`; `Σ_i m_i = Σ_f A_f`; hence `Σ_i m_i/(Σ_f A_f) = 1`. | `TriMesh.sum_mass`, `sum_mass_div_totalArea` |
| `angle-defect` | `Ω_i = 2π − Σ(corner angles at i)` at interior vertices (stored as 0 at boundary vertices). | `TriMesh.gauss_bonnet`, `gauss_bonnet_fields` |
| `boundary-geodesic-curvature` | `k_i = π − Σ(corner angles at i)` at boundary vertices (stored as 0 at interior vertices). | `TriMesh.gauss_bonnet`, `gauss_bonnet_fields` |
| `euler-characteristic` | `χ = \|V\| − \|E\| + \|F\|`. | `TriMesh.gauss_bonnet` |
| `gauss-bonnet` | Valid mesh, every edge in at most two faces, every boundary vertex on exactly two boundary edges: `Σ_{i interior} Ω_i + Σ_{i boundary} k_i = 2πχ`. Closed: `Σ_i Ω_i = 2πχ`. | `TriMesh.gauss_bonnet`, `gauss_bonnet_closed`, `gauss_bonnet_fields` |
| `scale-invariance` | Under `ℓ ↦ s·ℓ` (`s ≠ 0`): `θ`, `cot θ`, `L`, `Ω` unchanged; `A ↦ s²A`, `m ↦ s²m`. | `cornerAngle_smul`, `cornerCot_smul`, `triArea_smul`, `TriMesh.laplacian_scale`, `mass_scale`, `angleDefect_scale` |
| `discrete-yamabe-linearization` | For `ℓ̃_ij(t) = e^{t(v_i+v_j)/2} ℓ_ij`: `d/dt θ_i = ½cot θ_k (v_j − v_i) + ½cot θ_j (v_k − v_i)` in a triangle `ijk`, and `d/dt Ω_i = (L v)_i`, at `t = 0`. | `hasDerivAt_cornerAngle_conformal`, `TriMesh.hasDerivAt_angleSum_conformalScale`, `hasDerivAt_angleDefect_conformalScale` |

### The optimization

| id | exact statement | theorems |
|---|---|---|
| `soft-threshold-fixed-point` | For `t ≥ 0`: `µ = P_t(φ + µ)` ⟺ `\|φ\| ≤ t ∧ φ·µ = t·\|µ\|` ⟺ (`\|φ\| ≤ t`, `φ = t` if `µ > 0`, `φ = −t` if `µ < 0`) ⟺ `t\|µ\| + φ(y − µ) ≤ t\|y\|` for all `y`. For `t = 0`: ⟺ `φ = 0`. | `eq_softThreshold_add_iff`, `eq_softThreshold_add_iff_cases`, `eq_softThreshold_add_iff_isSubgradient`, `eq_softThreshold_zero_add_iff` |
| `abs-subdifferential-bounded` | `φ ∈ ∂(t\|·\| + ι_[lo,hi])(µ)` ⟺ `lo ≤ µ ≤ hi` ∧ (`µ < hi` ⟹ `φ ≤ (t if µ ≥ 0 else −t)`) ∧ (`lo < µ` ⟹ `(t if µ > 0 else −t) ≤ φ`). | `isSubgradientWithinAt_abs_iff` |
| `cone-fixed-point-bounded` | The inclusion above ⟺ `µ = clamp(P_t(φ + µ))` (cone residual `µ − clamp(P_t(φ + µ)) = 0`). At a fixed vertex (forbidden: `a = 0`; prescribed: `a`): ⟺ `µ = a`. | `isSubgradientWithinAt_abs_iff_eq_clampE`, `Site.stationary_iff_eq_prox` |
| `regularized-cone-update` | For `γ' > 0`: `φ − γ'µ ∈ ∂(t\|·\| + ι_[lo,hi])(µ)` ⟺ `µ = clamp(P_t(φ)/γ')`, the unique minimizer of `y ↦ t\|y\| + (γ'/2)y² − φy` on `[lo, hi]`. With `γ' = γ/m`: `µ = clamp((m/γ)·P_t(φ))`. Equivalent implicit form: `µ = clamp(c·P_t(φ + µ))`, `c = m/(m + γ)`. | `isSubgradientWithinAt_abs_sub_iff_eq_clampE`, `minimizes_reg_iff_eq_clampE`, `Oracle.eq_regularizedMeasure_add_iff` |
| `soft-threshold-active-sets` | `x ↦ clamp(k·P_t(x))` (`k > 0`; `k = 1` in Stage II, `k = m/γ` in Stage I) equals on each piece an affine function: clamped-high `hi` (slope 0), positive-active `k(x − t)` (slope `k`), inactive `0`, negative-active `k(x + t)` (slope `k`), clamped-low `lo` (slope 0). For the unbounded soft threshold: `P_t(x) = A⁺(x)(x − t) + A⁻(x)(x + t)`, slope `A⁺ + A⁻`, `A⁺ = [x > t]`, `A⁻ = [x < −t]`. | `softThreshold_isPiecewiseAffine`, `Site.proxReg_isPiecewiseAffine`, `Oracle.regularizedMeasure_eq_value` |
| `optimality-system-sufficient` | If `L u = µ − Ω`, `Lᵀφ = −w∘u`, and at every vertex `φ_i ∈ ∂(t_i\|·\| + ι_[lo_i,hi_i])(µ_i)` (resp. `µ_i = a_i` at a fixed vertex), then for every `(v, ν)` with `L v = ν − Ω` and admissible `ν`: `J(u, µ) + ½ Σ w_i (v_i − u_i)² ≤ J(v, ν)`. Any `L`; `w ≥ 0`. | `ConeProblem.KKT.objective_add_le`, `KKT.isMinimizer`, `KKT.minimizes_reduced_of_posDef`, `KKT.minimizes_closed` |
| `optimality-system-necessary` | If `det L ≠ 0` and `(u, µ)` is a minimizer of a well-formed problem, then `φ = −(Lᵀ)⁻¹(w∘u)` is the unique vector completing it to a solution of the optimality system. | `ConeProblem.IsMinimizer.existsUnique_kkt`, `isMinimizer_iff_exists_kkt_of_posDef` |
| `minimizer-unique` | If all `w_i > 0`, two minimizers have equal `u` and equal `µ` (any `L`). | `ConeProblem.IsMinimizer.unique` |
| `least-distortion-guarantee` | For a minimizer `(u, µ)` and feasible `(v, ν)`: `E(u) ≤ E(v) + λ(R(ν) − R(µ))`, `E = ½Σ w u²`, `R = Σ wR\|·\|`; so `E(u) ≤ E(v)` whenever `R(ν) ≤ R(µ)`. | `ConeProblem.IsMinimizer.distortion_le_add`, `distortion_le` |
| `duality-gap-certificate` | Dual value `D(ũ, φ) = −½Σ w_i ũ_i² + Σ φ_i Ω_i − Σ c_i(φ_i)`, where `Lᵀφ = −w∘ũ`, `c_i(p) = hi_i·max(0, p − t_i) + lo_i·min(0, p + t_i)` (a term with an infinite bound is dropped and then `p ≤ t_i`, resp. `−t_i ≤ p`, is required), `c_i(p) = a_i p` at a fixed vertex. For `L u = µ − Ω`: `J(u, µ) − D(ũ, φ) = ½Σ w_i(u_i − ũ_i)² + Σ_i (t_i\|µ_i\| + c_i(φ_i) − φ_i µ_i)`, every summand `≥ 0` at feasible points; hence `D ≤ J(v, ν)` for all feasible `(v, ν)`, and `J(u, µ) ≤ J(v, ν) + gap` whenever `J(u, µ) − D ≤ gap`. If all `w_i > 0`: `ũ_i = −(Lᵀφ)_i/w_i`, `D = −½Σ (Lᵀφ)_i²/w_i + Σ φ_i Ω_i − Σ c_i(φ_i)`. At a solution of the optimality system `D = J`. | `ConeProblem.gap_identity`, `weak_duality`, `duality_gap_certificate`, `adjoint_gap_certificate`, `KKT.dualObjective_eq` |
| `duality-gap-inexact` | For **arbitrary** `u, ũ, φ, µ` with `r = L u − µ + Ω`, `s = Lᵀφ + w∘ũ`: `J(u, µ) − D(ũ, φ) = ½Σ w_i(u_i − ũ_i)² + Σ_i (t_i\|µ_i\| + c_i(φ_i) − φ_i µ_i) − Σ φ_i r_i + Σ s_i u_i`. With `ũ_i = −(Lᵀφ)_i/w_i` (`0` where `w_i = 0`): `J(u, µ) − D = G − Σ φ_i r_i + Σ_{w_i=0} (Lᵀφ)_i u_i`, `G` the first two sums (the "complementarity gap"). | `ConeProblem.gap_identity_general`, `objective_sub_dualObjectiveOfAdjoint` |
| `duality-gap-certificate` (inexact data) | `φ̂ = dualClip φ` (`φ_i` clamped to `[−t_i` if `lo_i = −∞`, `t_i` if `hi_i = +∞]`), `D = −½Σ_{w_i>0} (Lᵀφ̂)_i²/w_i + Σ φ̂_i Ω_i − Σ c_i(φ̂_i)`. If `(Lᵀφ̂)_i = 0` wherever `w_i = 0`: `D ≤ J(v, ν)` for all feasible `(v, ν)`. If moreover `µ_i + r_i` is admissible for all `i`: `(u, µ + r)` is feasible, `J(u, µ + r) ≤ J(u, µ) + Σ t_i\|r_i\|` (sum over penalized vertices), and `J(u, µ + r) ≤ J(v, ν) + (J(u, µ) − D + Σ t_i\|r_i\|)` for all feasible `(v, ν)`. Always (for admissible `µ`): `(u, µ)` is feasible for the problem with `Ω − r` and within `G` of its optimum. | `ConeProblem.dual_lower_bound_of_clip`, `certificate_enclosure`, `certificate_of_computed`, `IsMinimizer.weighted_sq_sub_le_width`, `repaired_certificate`, `backward_certificate`, `KKT.objective_add_le_of_residual` |
| `newton-stopping-rule` | Stage II: `F(u, φ, µ) = (L u − µ + Ω, Lᵀφ + w∘u, µ − clamp(P_t(φ + µ)))`, `DF·(v, q, ν) = (L v − ν, Lᵀq + w∘v, ν − D(q + ν))` with `D_i` the slope of the piece of `φ_i + µ_i`. If `DF·(v, q, ν) = −F(u, φ, µ)` and every `φ_i + q_i + µ_i + ν_i` lies on the same piece as `φ_i + µ_i`, then `F(u + v, φ + q, µ + ν) = 0` and `(u + v, µ + ν)` is a global minimizer. Stage I: the same for `F_γ(u, φ) = (L u + Ω − clamp((m/γ)P_t(φ)), Lᵀφ + w∘u)`, `DF_γ·(v, q) = (L v − D_γ q, Lᵀq + w∘v)`, pieces of `φ_i`. | `paper_stageII_stopping_rule`, `paper_stageI_stopping_rule`, `ConeProblem.isMinimizer_of_newton_step`, `isRegMinimizer_of_newton_step`, `Ddg.Problem.Fimp_eq_zero_of_newton_step` |
| `regularized-optimality-sufficient` | If `L u = µ − Ω`, `Lᵀφ = −w∘u` and `φ_i − γ_i µ_i ∈ ∂(t_i\|·\| + ι)(µ_i)` (`γ_i ≥ 0`), then `(u, µ)` minimizes `J(u, µ) + Σ (γ_i/2) µ_i²` over all feasible pairs. (`γ_i = γ/m_i` is the mass-consistent regularization.) | `ConeProblem.isRegMinimizer_of_stationary`, `isRegMinimizer_of_residualReg_eq_zero` |
| `regularized-suboptimality` | If `(u_γ, µ_γ)` minimizes the regularized objective and `(u, µ)` is feasible: `J(u_γ, µ_γ) ≤ J(u, µ) + Σ (γ_i/2) µ_i²`. | `ConeProblem.IsRegMinimizer.objective_le` |
| `closed-mesh-cone-sum` | If `1ᵀL = 0` and `L u = µ − Ω` then `Σ µ_i = Σ Ω_i`; on a closed valid manifold mesh this is `2πχ`. | `ConeProblem.sum_cone_eq_sum_defect`, `TriMesh.sum_cone_eq_two_pi_eulerChar` |
| `closed-mesh-optimal-constant` | `E(u + c) = E(u) + c·Σ w_i u_i + (c²/2)·Σ w_i`; if `Σ w_i > 0` the unique minimizing constant is `c₀ = −Σ w_i u_i / Σ w_i`. If `L·1 = 0`, every minimizer (and every solution of the optimality system) has `Σ w_i u_i = 0`. | `ConeProblem.distortion_add_const`, `distortion_add_const_eq`, `KKT.weighted_mean_eq_zero`, `IsMinimizer.weighted_mean_eq_zero` |
| `lambda-area-normalization` | For `c > 0`: the minimizers with masses `c·m` and parameter `λ` are the minimizers with masses `m` and parameter `λ/c`. With masses divided by the total area `A` the problem is unchanged by `(m, A) ↦ (c·m, c·A)`; dividing the masses by `A` is `λ ↦ λ·A`. | `ConeProblem.isMinimizer_scaleMass_iff`, `normalizeMass_scaleMass`, `isMinimizer_normalizeMass_iff` |
| `component-separation` | Two problems `P` on `ι`, `Q` on `κ` with a common `λ`; `P.sum Q` on `ι ⊕ κ` has `L = [[L_P, 0], [0, L_Q]]` and the data of `P`, `Q` on the two blocks (the problem of a mesh with the two as its connected components, masses of both normalized by the same total area). Then `J(u ⊕ v, µ ⊕ ν) = J_P(u, µ) + J_Q(v, ν)`; `(u ⊕ v, µ ⊕ ν)` is feasible ⟺ `(u, µ)` and `(v, ν)` are; it is a global minimizer ⟺ `(u, µ)` minimizes `P` and `(v, ν)` minimizes `Q` (on each closed component the additive constant of `u` is a separate free variable); `(u ⊕ v, φ ⊕ ψ, µ ⊕ ν)` solves the optimality system ⟺ each part does. If `1ᵀL_P = 0`, every feasible pair has `Σ_ι µ = Σ_ι Ω`. Combined with the rescaling below: minimizer of the whole ⟺ each part minimizes its component with masses `a·m`, `λ_a = a·λ` (any `a > 0` per component). | `ConeProblem.objective_sum`, `feasible_sum_iff`, `isMinimizer_sum_iff`, `kkt_sum_iff`, `isMinimizer_sum_iff_alone` (Gauss–Bonnet: `Feasible.sum_cone_eq_sum_defect_inl`, `…_inr`, tagged `closed-mesh-cone-sum`) |
| `component-lambda-rescaling` | For `c > 0`: `c·φ ∈ ∂(c·λ·cost + ι)(µ)` ⟺ `φ ∈ ∂(λ·cost + ι)(µ)`; masses `c·m` with `c·λ` have the minimizers of masses `m` with `λ`, and `(u, φ, µ)` solves the optimality system of the latter ⟺ `(u, c·φ, µ)` solves that of the former. With raw masses `a_i`, total area `A > 0` and component area `A_c > 0`: masses `a_i/A_c` with `λ_c = (A/A_c)·λ` have the minimizers of masses `a_i/A` with `λ`, and `φ_global = (A_c/A)·φ_alone`. | `Site.stationary_mul_iff`, `ConeProblem.isMinimizer_scaleMass_withLam_iff`, `kkt_scaleMass_withLam_iff`, `kkt_scaleMass_withLam_iff'`, `isMinimizer_normalizeMass_alone_iff`, `kkt_normalizeMass_alone_iff` |

### The cone solver (`crates/ddg-optimal-cones`)

`A` is the Laplacian block of the unknowns, `w_i = wE_i m_i`, `g_i = t_i|·| + ι_[lo_i,hi_i]`.

| id | exact statement | theorems |
|---|---|---|
| `saddle-point-newton-system` | A Newton system has one row per vertex: *flat* `(A u)_i = µ_i − Ω_i`, *sloped* `(A u)_i − k_i(φ_i − o_i) = b_i − Ω_i`, *pinned* `φ_i = p_i`; plus `w∘u + A φ = 0`. `(u, φ)` solves it ⟺ `A φ = −w∘u` and `µ := A u + Ω` satisfies `µ_i` given / `µ_i = b_i + k_i(φ_i − o_i)` / `φ_i = p_i`. The assembled matrix form (unknown `ψ_i = φ_i − o_i` on sloped rows, pinned columns removed, known values on the right) has the same solutions. If `A` is symmetric and `φ_i ∈ ∂g_i(µ_i)` for all `i`, the solution is a solution of the optimality system, hence a global minimizer; with `φ_i − (µ_i − µ̄_i)/k_i ∈ ∂g_i(µ_i)` it minimizes `J + Σ(µ_i − µ̄_i)²/(2k_i)`. | `Solver.isSolution_iff`, `Solver.isAssembledSolution_iff`, `ConeProblem.kkt_of_isSolution`, `isMinimizer_of_isSolution`, `isProxMinimizer_of_isSolution` |
| `soft-threshold-active-sets` (rows) | "The vertex is on the piece of its row": positive cone `g = t`, `0 ≤ µ ≤ hi`; negative cone `g = −t`, `lo ≤ µ ≤ 0`; no cone `µ = 0`, `g ≤ t` unless `hi = 0`, `−t ≤ g` unless `lo = 0`; upper bound `µ = hi`, `t ≤ g` (`−t ≤ g` if `hi = 0`); lower bound symmetric; fixed `µ = a`. Each implies `g ∈ ∂g_i(µ)`. | `Site.stationary_of_pos`, `…_neg`, `…_zero`, `…_upper`, `…_lower`, `Site.stationary_fixed_iff` |
| `proximal-fixed-point` | If `(u, µ)` minimizes `J(v, ν) + Σ(ν_i − µ_i)²/(2k_i)` (`k_i > 0`) over the feasible pairs, it minimizes `J`. If `(u, µ)` minimizes the problem with centre `µ̄` and `(ū, µ̄)` is feasible: `J(u, µ) + Σ(µ_i − µ̄_i)²/(2k_i) ≤ J(ū, µ̄)`. | `ConeProblem.IsProxMinimizer.isMinimizer_of_fixedPoint`, `IsProxMinimizer.objective_add_le` |
| `regularized-cone-update` (proximal form) | For `k > 0`: `φ − (µ − µ̄)/k ∈ ∂g(µ)` ⟺ `µ = clamp(k·P_t(φ + µ̄/k))`. | `ConeProblem.stationary_prox_iff` |
| `adjoint-constant-interval` | If `Lᵀ·1 = 0` and `(u, φ, µ)` satisfies the state and adjoint equations, then `(u, φ + c, µ)` solves the optimality system ⟺ `a_i ≤ φ_i + c ≤ b_i` for all `i`, `[a_i, b_i] = ∂g_i(µ_i)` (`a_i = −∞` at the lower bound and at a fixed vertex, else `t_i` if `µ_i > 0`, else `−t_i`; `b_i = +∞` at the upper bound and at a fixed vertex, else `t_i` if `µ_i ≥ 0`, else `−t_i`). | `ConeProblem.kkt_add_const_iff`, `kkt_add_const_iff_interval`, `Site.stationary_iff_mem_subdiff` |
| `gauss-bonnet-adjoint-shift` | `c ↦ Σ_i clamp(k_i·P_{t_i}(x_i + c))` is nondecreasing. If `1ᵀA = 0`, `Σ b = 0` and `(A u)_i = b_i` for `i ≠ i₀`, then `(A u)_{i₀} = b_{i₀}`: once the flat cone angles sum to `Σ Ω`, the system with one adjoint value pinned and that vertex's state equation dropped is the full system. | `ConeProblem.sum_proxReg_add_monotone`, `Solver.mulVec_eq_of_sum_eq`, `Solver.isSolution_of_pinned_of_sum_eq` |
| `regularized-dual-objective` | For symmetric `A` and a direction with `w∘du + A dφ = 0`: `g(α) := Σ w_i du_i(u_i + α du_i) − Σ Ω_i dφ_i + Σ m_i(x_i + α dφ_i) dφ_i = −Σ dφ_i·((A(u + α du))_i − m_i(x_i + α dφ_i) + Ω_i)`; `g` is nondecreasing if `w ≥ 0` and every `m_i` is. | `ConeProblem.dualSlope_eq`, `dualSlope_monotone` |
| `componentwise-backward-error` | (Oettli–Prager) `\|b − Kx\|_i ≤ ε(\|K\|\|x\| + \|b\|)_i` for all `i` ⟺ `(K + δK)x = b + δb` for some `\|δK_ij\| ≤ ε\|K_ij\|`, `\|δb_i\| ≤ ε\|b_i\|`. | `oettli_prager`, `residual_eq_zero_of_magnitude_eq_zero` |
| `cone-cluster-frechet-mean` | Merging a cluster `C` into one cone at any `v ∈ C` (`µ'_v = Σ_C µ`, `µ' = 0` on `C∖{v}`): `Σ µ' = Σ µ`, `Σ\|µ'\| ≤ Σ\|µ\|`, and `µ'` is admissible ⟺ `Σ_C µ` is admissible at `v` and `0` is admissible on `C∖{v}`. A vertex of least Fréchet variance exists. | `sum_mergeCones`, `sum_abs_mergeCones_le`, `mergeCones_mem_feasible_iff`, `exists_frechetMean` |
| `closed-mesh-cone-sum` (evaluation) | `Σ_i (µ_i − Ω_i − (Σµ − ΣΩ)·m_i/Σm) = 0`. | `sum_sub_defect_eq_zero` |
| `closed-mesh-optimal-constant` (solver) | If `1ᵀA = 0`, every solution of a Newton system has `Σ w_i u_i = 0`. | `Solver.IsSolution.sum_weighted_eq_zero` |
| `boundary-optimality-system` | Let `Ω̂_i = Ω_i` (`i` interior), `Ω̂_i = k_i = π − Σ(corner angles)` (`i` boundary), `s = (µ on I, ν on B)`, one site `g_i` per vertex (`λ·wR_i|·|` with bounds), `m_i` the lumped masses divided by the total area, `wE ≥ 0`, `M` valid. **(1)** If at every vertex `(L u)_i = s_i − Ω̂_i`, `(L φ)_i = −wE_i m_i u_i` and `φ_i ∈ ∂g_i(s_i)` (`Site.Stationary`), then `(u, s)` is feasible and minimizes `½ Σ_V wE_i m_i u_i² + λ Σ_V cost_i(s_i)` over all `(v, ν)` with `L v = ν − Ω̂`, `ν_i` admissible, including all shifts `v + c` (closed-mesh case: `L·1 = 0`). **(2)** On a valid manifold mesh every `(u, s)` with `L u = s − Ω̂` has `Σ_I µ + Σ_B ν = Σ_V s = 2πχ`; also `Σ Ω̂ = 2πχ`. **(3)** Boundary sites: free (`g = 0`): `φ_i ∈ ∂g(ν_i) ⟺ φ_i = 0` (`λ ≥ 0`); penalized (`g = λ·wR·\|·\|`, `λ, wR ≥ 0`): `⟺ ν = P_{λ wR}(φ + ν)` `⟺ \|φ\| ≤ λ wR ∧ (ν > 0 → φ = λ wR) ∧ (ν < 0 → φ = −λ wR)`; convex (`g` the indicator of `[0, ∞)`, `Site.convex = penalized 0 0 ⊤`, no extension of `Site` needed): `⟺ ν ≥ 0 ∧ φ ≤ 0 ∧ φ ν = 0`. **(4)** Modelling step: for `t ↦ Ω̂_i` of the metric `ℓ̃_ij = e^{t(u_i+u_j)/2} ℓ_ij` (geodesic curvature at a boundary vertex, angle defect at an interior one): value `Ω̂_i` at `0`, derivative `(L u)_i` at `0`; hence for feasible `(u, s)`, `s_i = Ω̂_i + (L u)_i` is the first-order value (only first order is claimed). | `TriMesh.boundaryProblem`, `boundaryProblem_feasible_iff`, `boundaryProblem_objective`, `boundaryProblem_isMinimizer`, `sum_hatCurvature`, `sum_eq_of_feasible`, `Site.stationary_boundary_free_iff`, `Site.stationary_penalized_iff`, `Site.stationary_convex_iff`, `hasDerivAt_hatCurvature_conformalScale`, `boundaryProblem_firstOrder_boundary`, `boundaryProblem_firstOrder_interior` |
| `reduced-lasso-problem` | `ConeProblem P` with `A = P.L`, `w_i = wE_i m_i`, `S` the free (penalized) sites, `F` the fixed sites (forbidden `a = 0`, prescribed `a`), `a_F` the fixed angles (zero on `S`), `R` a linear solver (`R = A⁻¹` or the pinned solver below), `z_j = R e_j`, `w0 = R(Ω − a_F)`, `G_jk = z_jᵀWz_k`, `q_j = −z_jᵀWw0`. Then for every `µ`: `u := Σ_{j∈S} µ_j z_j − w0 = R(µ_S − (Ω − a_F))`, `½ Σ w u² = ½ µ_SᵀGµ_S + qᵀµ_S + ½ w0ᵀWw0`, hence `J(u, µ) = [½ µ_SᵀGµ_S + qᵀµ_S + λ Σ_S wR_j\|µ_j\|] + ½ w0ᵀWw0` (fixed sites cost 0). **Dirichlet** (`A` positive definite, `R = A⁻¹`): `(u, µ)` is a minimizer of the canonical problem ⟺ `µ` is feasible (`µ_F = a_F`, `lo_j ≤ µ_j ≤ hi_j` on `S`), `u = Z µ_S − w0`, and `µ` minimizes the bracket over all feasible cone vectors. **Closed** (`A` PSD, kernel the constants, `w > 0`, `R` = pinned solver): the same with the extra constraint `Σ_S µ_S = Σ_i (Ω − a_F)_i` on both `µ` and the competitors. `t_j = λ wR_j` is `λ·Site.cost`. | `ConeProblem.reducedU_eq`, `distortion_reducedU`, `objective_reducedU`, `isMinimizer_iff_reduced_dirichlet`, `isMinimizer_iff_reduced_pinned` |
| `reduced-gram-matrix` | Same notation, `G = (z_jᵀ W z_k)_{jk}`, `q_j = −z_jᵀ W w0`, `z_j = R e_j`, where in the closed case `R b = P·pad(A_p⁻¹ b_p̂)` (so `z_p = 0`). If all `w_i ≥ 0` then `G` is positive semidefinite; if all `w_i > 0` and the columns `z_j`, `j ∈ T`, are linearly independent then `xᵀGx > 0` for every nonzero `x` supported on `T`. The columns `z_j`, `j ∈ T`, are linearly independent for every `T` in the Dirichlet case (`A` invertible) and for every `T` not containing the pin in the closed case (so `G` is positive definite on all of `S` resp. `S ∖ {p}`). For a free site `j` and any `t`, `½ Σ w (u + t z_j)² = ½ Σ w u² + t·z_jᵀWu + (t²/2) G_jj`, i.e. `∂/∂µ_j (½ uᵀWu) = z_jᵀ W u`. | `ConeProblem.gramMatrix_posSemidef`, `gram_quadratic_pos`, `linearIndepOn_column_dirichlet`, `linearIndepOn_column_pinned`, `distortion_reducedU_update` |
| `reduced-optimality-conditions` | **Gradient.** `z_jᵀ W u = y_j`, `y = A⁻¹(W u)` (Dirichlet, `A` symmetric); closed case: `y = pad(A_p⁻¹ (W u)_p̂)` (so `y_p = 0`) for every `u` with `Σ w u = 0`. **Reduced conditions** `ReducedOptimal µ y ν`: `µ_F = a_F` and `−(y_j + ν) ∈ ∂g_j(µ_j)` for `j ∈ S`, `g_j(m) = t_j\|m\| + ι_[lo_j,hi_j](m)` (`Site.Stationary`): at `µ_j = 0` it says `\|y_j + ν\| ≤ t_j` restricted to the signs the bounds allow, at `µ_j = hi_j > 0`: `y_j + ν + t_j ≤ 0`, at `µ_j = lo_j < 0`: `y_j + ν − t_j ≥ 0`. **Link to `KKT`** (the direction the code relies on): Dirichlet (`A` positive definite): `µ` with `ReducedOptimal µ y 0`, `u = Z µ_S − w0`, `y = A⁻¹(Wu)` ⟹ `KKT u (−y) µ`, hence `(u, µ)` is a global minimizer (`optimality-system-sufficient`; needs `WellFormed`). Closed (`A` PSD with kernel the constants, `w > 0`, any pin `p`): `µ` with `Σ_S µ_S = Σ(Ω − a_F)`, `ReducedOptimal µ y ν`, `u = Z µ_S − w0`, `y = pad(A_p⁻¹ (Wu)_p̂)` ⟹ `KKT u (−(y + ν)) µ`; so `φ = −(y + ν)` satisfies `Lφ = −Wu` and `φ_j ∈ ∂g_j(µ_j)`. **Converse** (every minimizer satisfies the conditions): proved for the Dirichlet case only (`u = Z µ_S − w0`, `y = A⁻¹ W u`, `ν = 0`); for the closed case it is the cited `optimality-system-necessary-closed`. | `ConeProblem.column_dot_dirichlet`, `column_dot_pinned`, `column_dot_weighted_pinned`, `kkt_of_reducedOptimal_dirichlet`, `isMinimizer_of_reducedOptimal_dirichlet`, `kkt_of_reducedOptimal_pinned`, `isMinimizer_of_reducedOptimal_pinned`, `reducedOptimal_of_isMinimizer_dirichlet` |
| `reduced-sum-constraint-pin` | `A` PSD whose kernel is exactly the constants (`1ᵀA = 0`), `p` any index, `A_p` = `A` without row and column `p` (then invertible: proved), `Σ b = 0`: `x = pad(A_p⁻¹ b_p̂)` has `(Ax)_i = b_i` for `i ≠ p`, and `A x = b` (the pin's equation follows from `1ᵀA = 0`, `Σ b = 0`); with `Σ w ≠ 0`, `x̄ = P x`, `P x = x − (Σ w x/Σ w)·1`, also solves `A x̄ = b` and `Σ w x̄ = 0`; every solution of `A v = b` is `x̄ + c·1`, `x̄` is the only one with `Σ w v = 0`, and `E(x̄ + c) = E(x̄) + ½(Σ w)c²` for `E = ½ Σ w (·)²`, so `x̄` is the unique minimizer of the distortion among the solutions. For `b = µ − Ω` (fixed sites at `a_F`): `Σ b = 0 ⟺ Σ_S µ_S = Σ_i (Ω − a_F)_i`, and then `u = Z µ_S − w0` solves `L u = µ − Ω` with `Σ w u = 0`. For the mesh problem with every vertex an unknown on a connected valid mesh the hypotheses on `A` hold. | `ConeProblem.mulVec_pinSolve`, `mulVec_pinnedSolver`, `weighted_pinnedSolver`, `pinnedSolver_unique`, `distortion_pinned`, `reducedU_closed_spec`, `KernelIsConstants.pinBlock_isUnit`, `TriMesh.coneProblem_kernelIsConstants` |

### Boundary First Flattening (`crates/ddg-bff`)

`L` is the cotan Laplacian of a mesh, `I` its interior and `B` its boundary vertices; boundary
loops are indexed with the surface on the left, edge `j` from boundary vertex `j` to `j+1`.

| id | exact statement | theorems |
|---|---|---|
| `bff-harmonic-extension` | Connected valid mesh, `F` the free vertices, at least one vertex prescribed: for every source `φ` and prescribed values `g` there is exactly one `a` with `a = g` off `F` and `(L a)_v = φ_v` on `F` (`L_FF a_F = φ_F − L_FX a_X`). | `TriMesh.existsUnique_poisson` |
| `bff-dirichlet-to-neumann` | If `(L a)_v = φ_v` on `F`, then with `h_v := φ_v − (L a)_v` off `F` (`0` on `F`): `L a = φ − h` everywhere and `Σ h = Σ φ`. | `TriMesh.neumannData_spec`, `sum_neumannData` |
| `bff-neumann-to-dirichlet` | Connected valid mesh: `L a = b` is solvable ⟺ `Σ b = 0`; then for each vertex `v₀` exactly one solution has `a(v₀) = 0`, any two solutions differ by a constant, and for weights with `Σ m ≠ 0` exactly one solution has `Σ m a = 0`. | `TriMesh.exists_laplacian_mulVec_eq_iff`, `existsUnique_pinned`, `sub_const_of_laplacian_mulVec_eq`, `existsUnique_neumann` |
| `bff-cherrier-boundary-curvature` | Under `ℓ̃_ij = e^{t(v_i+v_j)/2}ℓ_ij`: `d/dt k_i = (L v)_i` at `t = 0` for `k_i = π − Σθ` (any vertex). The code's `k̃ = k − h`, `h = −(L u)_B`, is `k̃ = k + (L u)_B`. | `TriMesh.hasDerivAt_geodesicCurvature_conformalScale` |
| `bff-exterior-angle-sum` | Valid manifold mesh, any `u`: `Σ_I (Ω + L u) + Σ_B (k + L u) = 2πχ`. If `(L u)_I = −Ω_I`: `Σ_B (k + L u) = 2πχ` (`= 2π` for a disk). If `(L u)_I = µ_I − Ω_I`: `= 2πχ − Σ_I µ`. | `TriMesh.gauss_bonnet_conformal`, `sum_exteriorAngle`, `sum_exteriorAngle_cones` |
| `bff-cone-angle-sum` | `M'` a cut of the valid mesh `M` (same faces, vertex copies, `ρ` copy ↦ original), `u` on `M`: `Σ_{p : ρ p = v} (π − (k'_p + (L'(u∘ρ))_p)) = 2π − (Ω_v + (L u)_v)`; `= 2π − µ_v` if `(L u)_v = µ_v − Ω_v`. At a vertex with one copy `Ω' + L'(u∘ρ) = Ω + L u`. | `TriMesh.sum_fiber_interior_angle`, `sum_fiber_interior_angle_cone`, `curvature_copy_eq`, `sum_fiber_angleSum`, `sum_fiber_laplacian` |
| `cut-euler-characteristic` | `χ(M') = χ(M) + Σ_v (s_v − 1) − Σ_e (c_e − 1)`, `s_v` copies of vertex `v`, `c_e` copies of edge `e`; `c_e ≤` number of faces at `e` (so `c_e − 1 ∈ {0, 1}` on a manifold mesh, and the last sum is the number of cut edges). | `TriMesh.cut_euler_characteristic`, `card_vertices_cut`, `card_edges_cut`, `card_edge_copies_le` |
| `bff-hilbert-transform` | Boundary loop with successor `next`, `𝒜(a, b) = ½Σ_j (a_j b_{j+1} − a_{j+1} b_j)` (signed area of the polygon `(a_j, b_j)`, counter-clockwise positive), `h_j = ½(a_{j−1} − a_{j+1})`: `𝒜(a, b) = Σ_j h_j b_j`, `Σ_j h_j = 0`. For symmetric positive semidefinite `L` and `rhs = h` at the boundary vertices, `0` elsewhere: `L b = rhs` ⟺ `b` minimizes `½(aᵀLa + bᵀLb) − 𝒜(a, b)`. For the cotan Laplacian of a connected valid mesh such a `b` exists. | `loopSignedArea_eq_sum_hilbertRhs_mul`, `sum_hilbertRhs`, `hilbertTransform_minimizes`, `hilbertTransform_of_minimizes`, `hilbertTransform_minimizes_conformalEnergy_iff`, `TriMesh.hilbertTransform_cotan` |
| `bff-dual-boundary-length` | `d_j = ½(ℓ_{j−1} + ℓ_j)`: `Σ_j d_j = Σ_i ℓ_i`; `d_j > 0` if all `ℓ > 0`. | `sum_dualBoundaryLength`, `dualBoundaryLength_pos` |
| `bff-close-exterior-angles` | `k'_j = k_j + (d_j/Σd)(2π − Σk)`: if `Σd ≠ 0` then `Σk' = 2π`; nothing changes if `Σk = 2π`. | `sum_closedExteriorAngle`, `closedExteriorAngle_of_defect_eq_zero` |
| `bff-uniformization-angles` | `k_j = 2π d_j/Σd`: `Σk = 2π` if `Σd ≠ 0`, `k_j > 0` if all `d_j > 0`; in particular for `d` the dual lengths of the target lengths of positive `ℓ`. | `sum_uniformizationAngle`, `uniformizationAngle_pos`, `flattenToDisk_angles` |
| `bff-boundary-tangents` | `φ_0 = 0`, `φ_i = φ_{i−1} + k_i`, `T_i = e^{iφ_i}`: `\|T_i\| = 1`; `T_{i+1} = e^{ik_{i+1}} T_i`; `T_0 = e^{ik_0} T_{n−1}` ⟺ `Σ_{j<n} k_j ∈ 2πℤ` (in particular if the sum is `2π`). | `norm_boundaryTangent`, `boundaryTangent_succ`, `boundaryTangent_closes_iff`, `boundaryTangent_closes` |
| `bff-boundary-curve` | `γ_p = Σ_{i<p} ℓ_i T_i`: `γ_{p+1} − γ_p = ℓ_p T_p`, of length `\|ℓ_p\|` for unit `T_p`; the closing edge `γ_0 − γ_{n−1}` equals `ℓ_{n−1}T_{n−1}` ⟺ `Σ_{i<n} ℓ_i T_i = 0`. Gap spreading `p'_j = γ_j − (arc_j/total)γ_n`, `total ≠ 0`: `p'_n = 0 = p'_0` and `p'_{j+1} − p'_j = ℓ_jT_j − (ℓ_j/total)γ_n`. | `boundaryCurve_succ_sub`, `norm_boundaryCurve_succ_sub`, `boundaryCurve_closingEdge_iff`, `spreadCurve_end`, `spreadCurve_succ_sub` |
| `bff-closing-lengths` | `T : κ×δ`, weights `N_d > 0`, `S = T N⁻¹Tᵀ`: for any `λ`, `ℓ̃ = ℓ* − N⁻¹Tᵀλ` satisfies `Tℓ̃ = 0` ⟺ `Sλ = Tℓ*`; then `Σ N_d(ℓ̃_d − ℓ*_d)² ≤ Σ N_d(ℓ'_d − ℓ*_d)²` for every `ℓ'` with `Tℓ' = 0`, with equality only for `ℓ' = ℓ̃`. A solution `λ` always exists and all solutions give the same `ℓ̃`; if `S` is invertible, `ℓ̃ = ℓ* − N⁻¹Tᵀ(TN⁻¹Tᵀ)⁻¹Tℓ*` (BFF Eq. 20). One loop, `N_d = 1/ℓ_d`: the polygon with the lengths `ℓ̃` closes. | `closedLengths_closes_iff`, `closedLengths_optimal`, `closedLengths_unique`, `exists_closing_multiplier`, `closedLengths_eq_of_solutions`, `closedLengthsPaper_isMinimizer`, `boundaryCurve_loopClosedLength`, `loopClosedLength_optimal` |
| `signed-area-2d` | `sda(a,b,c) = (b−a)×(c−a)`: antisymmetric, cyclic, translation invariant; `4·sda² = 2xy + 2yz + 2zx − x² − y² − z²` in the squared side lengths; `\|sda\| = 2·(Heron area of the side lengths)`. | `sda_swap_left`, `sda_rotate`, `sda_translate`, `sda_sq_mul_four`, `abs_sda_eq_two_mul_triArea` |
| `quasi-conformal-error` | For a valid triangle `(ℓ01, ℓ12, ℓ20)` with its Heron area: `p0 = 0`, `p1 = (ℓ01, 0)`, `p2 = (x, y)` is an isometric layout with `y > 0`; `J = [J_s J_t]` maps `p1−p0 ↦ q1−q0`, `p2−p0 ↦ q2−q0`. With `S± = ½(a + c ± √((a−c)² + 4b²))`: `S₊ + S₋ = a + c`, `S₊S₋ = ac − b² = det²`, `0 ≤ S₋ ≤ S₊`, and `S₋ ≤ ‖J e‖² ≤ S₊` for unit `e`, both attained. `σ_min` as evaluated is `√S₋`; `σ_min > 0` ⟺ `det ≠ 0`; `Q ≥ 1`, `Q = 1` ⟺ `a = c ∧ b = 0` ⟺ `J` is a similarity. Dividing `J` by its largest entry and scaling back (as the code does) changes nothing. `sda(q) = det·2A`. | `layout_isometric`, `jacobian_maps_triangle`, `sqStretchMax_add_sqStretchMin`, `sqStretchMax_mul_sqStretchMin`, `sqStretch_bounds`, `sqStretch_attained`, `sigmaMin_eq_sqrt`, `sigmaMin_pos_iff`, `one_le_qcError`, `qcError_eq_one_iff`, `similarity_iff`, `singularValuesEval_eq`, `qcErrorEval_eq`, `sda_eq_jacDet_mul_area` |
| `face-log-scale-factor` | `u = ½ln(B/A)`: `= ln s` for `B = s²A`; `= ½ln\|det J\| = ½(ln σ_max + ln σ_min)` for the triangle map. | `logScaleFactor_sq_mul`, `logScaleFactor_abs_det_mul`, `logScaleFactor_eq_half_log_sigma`, `logScaleFactor_triangle` |
| `signed-angle-2d` | `θ = atan2(u×v, u·v) = arg(ū·v) ∈ (−π, π]`; `\|θ\|` is the angle between `u` and `v` (both nonzero); `0 < θ < π` ⟺ `u×v > 0`, `θ < 0` ⟺ `u×v < 0`, `θ = 0` resp. `π` for parallel resp. antiparallel vectors; the three signed corner angles of a counter-clockwise triangle are positive and sum to `π`. | `signedCornerAngle_eq_arg`, `signedCornerAngle_mem_Ioc`, `abs_signedCornerAngle_eq_angle`, `signedCornerAngle_sign`, `signedCornerAngle_sum_eq_pi` |
| `segment-proper-crossing` | `d1·d2 < 0 ∧ d3·d4 < 0` ⟺ the segments are not parallel and have a common point interior to both; that point is `a0 + t(a1−a0)`, `t = d3/(d3−d4)`. With a margin `≥ 0` the code's test still implies it (the tolerance itself is not modelled). | `properCrossing_iff`, `properCrossing_sound`, `properCrossing_complete`, `strictlyOpposite_imp` |
| `similarity-fit` | `a = Σ(q_i − q̄)conj(p_i − p̄)/Σ\|p_i − p̄\|²`, `b = q̄ − a p̄` is the unique minimizer of `Σ\|a p_i + b − q_i\|²` (points not all equal); the real sums the code forms are the real and imaginary parts of these. | `simFit_residual_le`, `simFit_residual_eq_iff`, `fit_similarity_optimal`, `fit_similarity_unique`, `toC_fitCoeff`, `toC_fitTranslation`, `fitNorm_eq_zero_iff` |
| `intrinsic-mollification` | `ε* = max(0, max over faces and sides (ℓ_c − ℓ_a − ℓ_b + δ))`, computed by the code's loop, is the least `ε ≥ 0` with `(ℓ_a + ε) + (ℓ_b + ε) ≥ (ℓ_c + ε) + δ` for all sides of all faces; if `δ > 0` every mollified face satisfies the strict triangle inequalities (for `δ = 0` it need not). | `mollifyEpsilon_le_iff`, `mollifyEpsilonLoop_eq`, `mollifyEpsilon_feasible`, `mollifyEpsilon_le_of_feasible`, `isTriangle_mollified` |

### Variational Surface Cutting (`crates/ddg-variational-cuts`)

Notation of SPEC: `φ` the level-set function (region 0 where `φ ≤ 0`), `Θ_v` the angle sum at a
vertex, `L` a square matrix on an index set `V` split by a predicate into the unknowns `I` and
the rest `B` (for a mesh: the cotan Laplacian of a region, `B` = cut vertices and vertices of
`∂M`), `u_I` the restriction to `I`, `L_II`, `L_BI` the blocks. "Oracle" names are the
definitions of `lean/oracle/DdgOracle/Formulas/Cuts.lean` at `ℝ`; the theorems are about them
unless a definition of this package is named. All theorems are in `namespace DdgProofs.Cuts`
(mesh statements in `DdgProofs.TriMesh`).

| id | exact statement | theorems |
|---|---|---|
| `cut-face-curvature` | `θ̃ = θ·2π/Θ_v` at `v ∉ B`, `θ̃ = θ` at `v ∈ B`; `Ω_f = θ̃_a + θ̃_b + θ̃_c − π`. If `Θ_v ≠ 0` and `v ∉ B`: `Σ_f θ̃_v^f = 2π`. A face with its three vertices in `B` has `Ω_f = 0`. On a valid mesh `Ω_f = Σ_{corners i ∉ B} θ_i (2π − Θ_i)/Θ_i` (the cancellation-free form, D13e), and the oracle's `Face.omegaCancellationFree` equals `Face.omega`. `k_f·A_f = Ω_f`. | `TriMesh.sum_normAngleAt_of_not_mem`, `faceCurvature_of_mem`, `faceCurvature_eq_sum_defect`, `sum_cornerNormalized`, `face_omega_eq`, `face_omegaCancellationFree_eq`, `faceCurvatureDensity_mul_area` |
| `cut-face-curvature-sum` | Valid mesh, `Θ_v ≠ 0` at every `v ∉ B` (true if `v` lies in a face): `Σ_f Ω_f = Σ_{v ∉ B} (2π − Θ_v)`. With `B` the boundary vertices of a valid manifold mesh: `Σ_f Ω_f + Σ_B k_v = 2πχ`; closed: `Σ_f Ω_f = 2πχ`. | `TriMesh.sum_faceCurvature`, `sum_faceCurvature_eq_sum_angleDefect`, `sum_faceCurvature_add_geodesicCurvature`, `sum_faceCurvature_closed`, `sum_face_omega`, `sum_face_omega_closed` |
| `cut-edge-crossing` | `t = φ_a/(φ_a − φ_b)` (`φ_a ≠ φ_b`): `(1−t)φ_a + tφ_b = 0`, and `t` is the only such number. `φ_a ≤ 0 < φ_b` ⟹ `0 ≤ t < 1`; `φ_b ≤ 0 < φ_a` ⟹ `0 < t ≤ 1`; `t = 0 ⟺ φ_a = 0`, `t = 1 ⟺ φ_b = 0`. Clamp `c(t) = max(δ, min(1−δ, t))`, `δ ≤ ½`: `δ ≤ c(t) ≤ 1−δ`; `c(t) = t ⟺ δ ≤ t ≤ 1−δ`; `crossingMoved δ t ⟺ c(t) ≠ t`; `\|c(t) − t\| ≤ δ` for `t ∈ [0,1]`, `δ ≥ 0`. | `interp_crossingParam`, `interp_eq_zero_iff`, `crossingParam_mem_Ico`, `crossingParam_mem_Ioc`, `crossingParam_mem_Icc`, `clampCrossing_mem`, `clampCrossing_eq_self_iff`, `crossingMoved_eq_true_iff`, `abs_clampCrossing_sub_le` |
| `cut-barycentric-length` | Points `P_a, P_b, P_c` of a real inner product space, `d_a + d_b + d_c = 0`: `‖Σ d_i P_i‖² = −(d_a d_b ℓ_ab² + d_b d_c ℓ_bc² + d_c d_a ℓ_ca²)` (`= baryLength2`; any finite family: `‖Σ d_i P_i‖² = −½ Σ_{i,j} d_i d_j ‖P_i − P_j‖²`); for `x = Σ x_i P_i`, `y = Σ y_i P_i` with coefficient sums 1, `baryLength(ℓ; x − y) = dist(x, y)`. If `ℓ_ab, ℓ_bc, ℓ_ca` satisfy the strict triangle inequalities, `baryLength2 ≥ 0` (the floor in `baryLength` is inactive). Polarization: for sum-zero `d`, `e`, `baryInner(ℓ; d, e) = ⟨Σ d_i P_i, Σ e_i P_i⟩`; and `cotFromInner(baryInner(ℓ; y−x, z−x), \|det\|·A_f)` is `cornerCot` of the corner `x` of the sub-triangle with the `baryLength` sides (D13a). | `baryLength2_eq_norm_sq`, `norm_sq_sum_smul_of_sum_eq_zero`, `baryLength_eq_dist`, `baryLength2_nonneg`, `baryInner_eq_inner`, `cotFromInner_eq_cornerCot` |
| `cut-subtriangle-area` | Parent side lengths with strict triangle inequalities, area `A_f` (Heron); barycentric rows `x, y, z` with sums 1: the Heron area of the triangle with sides `baryLength(x−y)`, `baryLength(y−z)`, `baryLength(z−x)` is `subTriangleArea(det3(x; y; z), A_f) = \|det\|·A_f`. For points: `triArea(dists of X, Y, Z) = \|det\|·triArea(dists of P_a, P_b, P_c)`. Cut face with `p = (1−t)a + tb`, `q = (1−s)a + sc`: the determinants of `(a,p,q)`, `(p,b,c)`, `(p,c,q)` are `ts`, `1−t`, `t(1−s)`, of `(p,b,q)`, `(b,c,q)`: `(1−t)s`, `1−s`; both triples sum to 1. | `subTriangleArea_eq_triArea`, `heronPoly_baryLength2`, `triArea_bary_points`, `cutFace_det_sum_pc`, `cutFace_det_sum_bq` |
| `cut-weighted-laplacian` | `(L_w)_ij = −½ Σ_f w̄_f cot(corner of f opposite ij)` (`i ≠ j`), zero row sums (`TriMesh.weightedLaplacian`; the oracle's `weightedOffDiag`, `weightedDiag` are its entries). `L_w = L_wᵀ`, `L_w·1 = 0`, `xᵀL_w x = Σ_f w̄_f·(cotan energy of x on f)` on a valid mesh, hence positive semidefinite for `w̄ ≥ 0`, in particular for `w̄_f = triangleMean` of a vertex weight `≥ 0`. `w̄ ≡ c` (or a constant vertex weight `c`): `L_w = c·L`. The paper's edge rule also gives `c·L` for a constant weight; on the triangle `4, 13, 15` with vertex weights `(1,1,0)` its energy at `(1,−1,−1)` is negative. | `TriMesh.weightedLaplacian_isSymm`, `weightedLaplacian_mulVec_one`, `dotProduct_weightedLaplacian_mulVec`, `weightedLaplacian_posSemidef`, `weightedLaplacian_faceMean_posSemidef`, `weightedLaplacian_const`, `weightedLaplacian_faceMean_const`, `edgeRule_const`, `edgeRuleEnergy_neg_example`, `weightedEntry_eq`, `weighted_quadratic_nonneg`, `weighted_row_sum`, `weightedEntry_symm` |
| `cut-yamabe-system` | `IsDirichletSolution L I f u`: `u = 0` off `I`, `(L u)_i = f_i` on `I` ⟺ `u_B = 0`, `L_II u_I = f_I`. If `det L_II ≠ 0` (e.g. `L_II` positive definite) there is exactly one such `u`, for every `f`; the Yamabe system is `f = −g`. For the cotan Laplacian of a connected valid mesh with a vertex outside `I` the hypothesis holds. | `existsUnique_isDirichletSolution`, `isDirichletSolution_iff_block`, `yamabe_existsUnique`, `yamabe_iff_block`, `isUnit_det_of_posDef`, `TriMesh.yamabe_existsUnique` |
| `cut-normal-derivative` | For a solution of the Yamabe system, `v_i := (L u)_i + g_i` at `i ∉ I`, `0` on `I` (`boundaryFlux L I (−g) u`): at `i ∉ I`, `v_i = normalFlux((L_BI u_I)_i, g_i)`; `L u = −g + v`, and `v` is the only vector with this property. If the columns of `L` sum to zero (for a symmetric `L`: the rows): `Σ v = Σ g` (flux balance); on a mesh with `g = Σ_{t∋i} c_t/3`: `Σ_B v = Σ_t c_t`. `normalDerivative(v, ℓ)·ℓ = v`; `dualLength [ℓ₁, ℓ₂] = (ℓ₁+ℓ₂)/2`. | `normalFlux_apply`, `normalFlux_spec`, `sum_normalFlux`, `sum_normalFlux_of_isSymm`, `TriMesh.sum_lumped`, `TriMesh.sum_normalFlux_lumped`, `normalDerivative_mul`, `dualLength_pair` |
| `cut-flux-recovery` | Cut vertices `C`, `ℓ_p > 0`, symmetric couplings `k(p,q) ≥ 0` (`κ_σ` for a segment `σ = pq`, else 0). `fluxMatrix` has the entries `fluxDiagEntry(ℓ_p, Σ_{q≠p} k(p,q))` and `fluxOffEntry(k(p,q))`; row `p` of `fluxMatrix·d` is `ℓ_p d_p + Σ_q k(p,q)(d_p − d_q)`. It is symmetric positive definite, so `fluxMatrix·d = v` has exactly one solution. For a set `S ⊆ C` that no coupling leaves (a piece of the cut; `S = C` included): `Σ_S ℓ_p d_p = Σ_S v_p`, and `m ≤ v_p/ℓ_p ≤ M` on `S` ⟹ `m ≤ d_p ≤ M` on `S`; `v_p = aℓ_p` on `S` ⟹ `d = a` on `S`. With `k ≡ 0`: `d_p = v_p/ℓ_p`. `v ↦ d` is linear. `fluxCoupling(c, h², s) = c²h²/s ≥ 0`, `= 0` for `c = 0`. | `fluxMatrix_mulVec_apply`, `fluxMatrix_isSymm`, `fluxMatrix_posDef`, `fluxRecovery_existsUnique`, `fluxRecovery_sum_piece`, `fluxRecovery_sum`, `fluxRecovery_ge_piece`, `fluxRecovery_le_piece`, `fluxRecovery_bounds`, `fluxRecovery_const_piece`, `fluxRecovery_nodal`, `fluxRecovery_smul`, `fluxRecovery_add`, `fluxCoupling_nonneg`, `fluxCoupling_zero` |
| `cut-dirichlet-energy` | For a solution of the Yamabe system: `uᵀLu = −Σ_i u_i g_i = −u_I·g_I`; `≥ 0` if `L` is positive semidefinite; on a valid mesh `= Σ_f (cotan energy of u on f)`. | `dirichletEnergy_eq`, `dirichletEnergy_eq_interior`, `dirichletEnergy_nonneg`, `TriMesh.dirichletEnergy_eq_sum_faceEnergy` |
| `cut-gradient-dirichlet` | `gradDirichletD(d) = d² ≥ 0`; with `d = v/ℓ` it is `gradDirichlet(v, ℓ) = (v/ℓ)²`. Closed-form instances: flat disk `d/dρ ∫₀^ρ u′² 2πr dr = gradDirichletD(u′(ρ))·2πρ`; spherical cap `d/dψ₀ ∫₀^{ψ₀} u′² 2π sin ψ dψ = gradDirichletD(u′(ψ₀))·2π sin ψ₀`. | `gradDirichletD_eq`, `gradDirichletD_nonneg`, `gradDirichletD_nodal`, `gradDirichlet_eq`, `hasDerivAt_disk_dirichletEnergy`, `hasDerivAt_disk_dirichletEnergy_density`, `hasDerivAt_capEnergy_integral`, `hasDerivAt_capEnergy_density` |
| `cut-adjoint-scaled` | `r_i = −2(L_a u)_i` on `I`, `r_i = −2((L_a u)_i − a_i v_i)` off `I` (`scaledAdjointRhs`, from the oracle's `adjointRhsScaledInterior/Boundary`); `q = boundaryFlux L I r p`, at `i ∉ I`: `q_i = adjointFlux((L_BI p_I)_i, r_i)`. If `L_a = cL`, `a ≡ c` and `u = 0` off `I`: `p = −2c·u` solves `p_B = 0`, `L_II p_I = r_I`, it is the only solution when `det L_II ≠ 0`, and `q_i = −2c·v_i` (`c = 1`: `p = −2u`, `q = −2v`). Exactness for symmetric `L`, `L_a`: if `δu_B = 0`, `L_II δu_I = −δg_I`, `p_B = 0`, `(Lp)_i = −2(L_a u)_i` on `I`, then `(u + tδu)ᵀL_a(u + tδu) = uᵀL_a u + t·Σ_I p_i δg_i + t²·δuᵀL_aδu`. | `adjointFlux_apply`, `scaledAdjoint_const`, `scaledAdjoint_const_unique`, `scaledEnergy_perturb`, `IsDirichletSolution.reciprocity` |
| `cut-gradient-scaled` | `gradScaledD(a, d, e) = −a d² − d e`; `gradScaledD(c, d, −2cd) = c·gradDirichletD(d)`, `gradScaledD(1, d, −2d) = gradDirichletD(d)`; with `d = v/ℓ`, `e = q/ℓ` it is `gradScaled(a, v, q, ℓ)`, and `gradScaled(1, v, −2v, ℓ) = gradDirichlet(v, ℓ)`. Flat disk, constant weight `c`, `p = −2cu`: `d/dρ [c ∫₀^ρ u′² 2πr dr] = gradScaledD(c, u′(ρ), p′(ρ))·2πρ`. | `gradScaledD_eq`, `gradScaledD_const`, `gradScaledD_one`, `gradScaledD_nodal`, `gradScaled_const`, `gradScaled_one`, `hasDerivAt_disk_scaledEnergy`, `hasDerivAt_disk_scaledEnergy_density` |
| `cut-adjoint-pointwise` | Symmetric `L`; `b` differentiable at every `u_i` with derivative `b′(u_i)`; `δu_B = 0`, `L_II δu_I = −δg_I`; `p_B = 0`, `(L p)_i = adjointRhsPointwise(A_i, b′(u_i)) = −A_i b′(u_i)` on `I`. Then `d/dt Σ_i A_i b(u_i + tδu_i)` at `t = 0` is `Σ_I p_i δg_i`. Hencky (`b = u²`): the same with `r_i = −2A_i u_i`. | `hasDerivAt_pointwiseEnergy`, `hasDerivAt_henckyEnergy` |
| `cut-gradient-pointwise` | `gradPointwiseD(b₀, d, e) = b₀ − d e`; with `d = v/ℓ`, `e = q/ℓ` it is `gradPointwise(b₀, v, q, ℓ)`. Flat disk, `b = u²`, adjoint `p = K r⁴/32 − Kρ²r²/8 + 3Kρ⁴/32`: `d/dρ ∫₀^ρ u² 2πr dr = gradPointwiseD(b(0), u′(ρ), p′(ρ))·2πρ = πK²ρ⁵/8`. | `gradPointwiseD_eq`, `gradPointwiseD_nodal`, `gradPointwise_eq`, `hasDerivAt_disk_henckyEnergy`, `hasDerivAt_disk_henckyEnergy_density`, `disk_hencky_derivative_value` |
| `cut-feature-size` | `featureSizeRaw(r_ε, K) = 1/(√\|K\| + 1/r_ε)`; for `r_ε > 0`: `0 < r⁰ ≤ r_ε`, `r⁰ = r_ε` at `K = 0`, `r⁰ < 1/√\|K\|` for `K ≠ 0`, `r⁰` decreasing in `\|K\|`; `featureWeight r = r³ > 0` for `r > 0`. Smoothing `(W_V + cL) r = W_V r⁰` with positive areas, `c ≥ 0`, zero row sums and **nonpositive off-diagonal entries** of `L`: every lower (upper) bound of `r⁰` is one of `r`; hence `0 < r ≤ r_ε`, `r³ > 0`. | `featureSizeRaw_pos`, `featureSizeRaw_le`, `featureSizeRaw_zero`, `featureSizeRaw_lt_inv_sqrt`, `featureSizeRaw_antitone`, `featureWeight_pos`, `smoothing_ge`, `smoothing_le`, `smoothedFeatureSize_bounds`, `smoothedFeatureSize_bounds_positivePart` (the Laplacian `L⁺` of the edge weights `max(w_ij, 0)` of SPEC D17 satisfies the sign condition on every mesh) |
| `cut-speed` | `cutSpeed(G₀, G₁) = G₀ − G₁`, `regionDensity(w, G, β) = wG + β`. With `dE(σ) := Σ_c g_c σ_c ℓ_c`: if `φ_t = −σ` and `φ_t = g` then `dE = −Σ g_c² ℓ_c` (`≤ 0` for `ℓ ≥ 0`, `= 0` only for `g = 0` when `ℓ > 0`); with `φ_t = −g` instead, `dE = +Σ g_c² ℓ_c`; `dE(σ) ≥ −Σ g²ℓ` for every `σ` with `Σσ²ℓ ≤ Σg²ℓ`. For `g = G₀ − G₁`, `σ = −cutSpeed(G₀, G₁)`: `dE = −Σ s_c² ℓ_c ≤ 0`. Closed forms: disk `d/dρ [wJ + βπρ²] = regionDensity(w, u′(ρ)², β)·2πρ`; two caps of a sphere `d/dψ₀ [J(ψ₀) + J(π − ψ₀)] = cutSpeed(tan²(ψ₀/2), tan²((π−ψ₀)/2))·2π sin ψ₀`. | `firstVariation_levelSet`, `firstVariation_levelSet_nonpos`, `firstVariation_levelSet_eq_zero_iff`, `firstVariation_ascent`, `firstVariation_ge`, `firstVariation_cutSpeed`, `hasDerivAt_disk_regionEnergy`, `hasDerivAt_twoCaps` |
| `cut-sliver-smoothing` | `smoothedSpeed(ℓ₁, ℓ₂, s₁, s₂) = (ℓ₂s₁ + ℓ₁s₂)/(ℓ₁+ℓ₂) = ℓ₂/(ℓ₁+ℓ₂)·s₁ + ℓ₁/(ℓ₁+ℓ₂)·s₂`; for `ℓ₁, ℓ₂ ≥ 0`, `ℓ₁+ℓ₂ > 0` the weights are `≥ 0` with sum 1 and `min(s₁,s₂) ≤ smoothedSpeed ≤ max(s₁,s₂)`; `smoothedSpeed(ℓ₁, ℓ₂, α − βℓ₁, α + βℓ₂) = α` (exact for a speed linear in arclength). `aspectRatio(longest, area) = longest²√3/(4·area)`: `= 1` for an equilateral triangle, `≥ 1` for every triangle with `longest` its longest side. | `smoothedSpeed_eq_convex`, `smoothedSpeed_weights`, `smoothedSpeed_mem`, `smoothedSpeed_linear`, `aspectRatio_equilateral`, `one_le_aspectRatio` |
| `cut-harmonic-extension` | If `det L_OO ≠ 0`: for all data `s` (and every source `φ`) there is exactly one `σ` with `σ = s` off `O` and `(Lσ)_i = 0` (`= φ_i`) on `O`; `(Lσ)_O = 0 ⟺ L_OO σ_O = −L_OC σ_C`. For symmetric `L`, such a `σ`, and `τ = σ` off `O`: `τᵀLτ = σᵀLσ + (τ−σ)ᵀL(τ−σ)`; so `σᵀLσ ≤ τᵀLτ` if `L` is positive semidefinite, with `<` for `τ ≠ σ` if `L_OO` is positive definite. On a connected valid mesh with a cut vertex the hypotheses hold and `σ` minimizes `Σ_f (cotan energy)`. | `existsUnique_dirichlet_data`, `existsUnique_harmonicExtension`, `harmonic_iff_block`, `harmonic_energy_identity`, `harmonicExtension_energy_le`, `harmonicExtension_energy_lt`, `TriMesh.existsUnique_harmonicExtension`, `TriMesh.harmonicExtension_faceEnergy_le` |
| `cut-explicit-step` | `explicitStep(φ, τ, σ) = φ + τσ`, so `φ̃ − φ = τσ`. With `dE` as in `cut-speed`, the displacement `−(φ̃_c − φ_c)` of the cut for the speed `s`, `τ ≥ 0`, `ℓ ≥ 0`: `dE = −τ Σ s_c² ℓ_c ≤ 0`. | `explicitStep_sub`, `firstVariation_explicitStep` |
| `cut-length-step` | `A_i > 0`, `L` positive semidefinite, `τ ≥ 0`: `diag(A) + τL` is positive definite (entries `lengthStepEntry`), `(diag(A) + τL)φ′ = diag(A)φ̃` has exactly one solution, and for it `φ′ᵀLφ′ ≤ φ̃ᵀLφ̃`, `ΣA_iφ′_i² ≤ ΣA_iφ̃_i²`, `ΣA_iφ′_i = ΣA_iφ̃_i` (zero column sums), and `φ′` minimizes `½ΣA_i(ψ_i − φ̃_i)² + (τ/2)ψᵀLψ` over all `ψ`. | `lengthStep_posDef`, `lengthStep_existsUnique`, `lengthStep_entry`, `lengthStep_energy_le`, `lengthStep_mass_le`, `lengthStep_sum`, `lengthStep_minimizes` |
| `cut-penalty-step-paper` | `d_i > 0`: `(diag(A/d) + τL)φ′ = diag(A/d)φ̃ ⟺ (diag(A) + τ·diag(d)·L)φ′ = diag(A)φ̃` (entries `penaltyStepEntry`); exactly one solution for `A > 0`, `L` positive semidefinite, `τ ≥ 0`; for constant `d` it is the length step with `L_w = dL`. | `penaltyStep_iff`, `penaltyStep_entry`, `penaltyStep_existsUnique`, `penaltyStep_const_iff` |
| `cut-symmetry-projection` | `classMean cls φ v` = mean of `φ` over `{w : cls w = cls v}` (`P`; the oracle's `classMean` of the listed class is this). `Pψ = ψ ⟺ ψ` is constant on classes; `P(Pφ) = Pφ`; `Σ(Pφ)_vχ_v = Σφ_v(Pχ)_v`; for class-constant `ψ`: `Σ(φ − Pφ)_vψ_v = 0` and `Σ(φ−ψ)² = Σ(φ−Pφ)² + Σ(Pφ−ψ)²`, so `Pφ` is the nearest class-constant vector. | `classMean_eq_self_iff`, `classMean_idem`, `classMean_selfAdjoint`, `sum_sub_classMean_mul`, `classMean_pythagoras`, `classMean_nearest`, `oracle_classMean_eq` |
| `cut-redistance` | `hopfLaxEdge(T, ℓ) = T + ℓ`. For points `x_i, x_j, x_k` of a real inner product space (any, also `x_j = x_k`): `hopfLaxTriangle(\|x_i−x_j\|, \|x_i−x_k\|, \|x_j−x_k\|, T_j, T_k)` is the least value of `λT_j + (1−λ)T_k + ‖x_i − (λx_j + (1−λ)x_k)‖` over `λ ∈ [0,1]` (attained, and a lower bound). For side lengths with the strict triangle inequalities: the same with the distance `baryLength(ℓ_ij, ℓ_jk, ℓ_ik; 1, −λ, −(1−λ))`. Consequences: `≤ hopfLaxEdge(T_j, ℓ_ij)`, `≤ hopfLaxEdge(T_k, ℓ_ik)`, `≥ min(T_j, T_k)`. | `hopfLaxTriangle_isLeast`, `hopfLaxTriangle_le`, `hopfLaxTriangle_mem`, `hopfLaxTriangle_isLeast_intrinsic`, `hopfLaxTriangle_le_edge`, `min_le_hopfLaxTriangle`, `hl_stationary_le`, `clamp_min_of_convex`, `hopfLaxObjective_convex` |
| `cut-step-cap` | `stepCap(τ_max, c, h, m) = min(τ_max, c·h/m)`: `τ ≤ τ_max`; for `m > 0`: `τ·m ≤ c·h`, `τ\|s\| ≤ c·h` whenever `\|s\| ≤ m` (and `τ_max, c, h ≥ 0`), and every `τ′ ≤ τ_max` with `τ′m ≤ ch` has `τ′ ≤ τ`. `adaptiveStep(τ_max, time_scale, ω, max_motion, h, m)` and `robustStep(τ_max, time_scale, ω, max_motion, m̂) = adaptiveStep(…, h = 1, m̂)` are `min(τ_max, time_scale/ω, max_motion·h/m)` **without the terms whose denominator is not positive**: `τ_max` for `ω ≤ 0`, `m ≤ 0`; `stepCap` (the capped rule) for `ω ≤ 0 < m`; `min(τ_max, time_scale/ω)` for `m ≤ 0 < ω`. `robustStep`: `τ ≤ τ_max`; `τω ≤ time_scale` if `ω > 0`; `τm̂ ≤ max_motion` if `m̂ > 0`; every `τ′` with these three properties has `τ′ ≤ τ` (it is the largest such step). For `τ_max, time_scale, max_motion ≥ 0`: `τ ≥ 0` and `τm̂ ≤ max_motion` whatever the sign of `m̂`. Robust statistics (`h_c > 0`): `speedInEdges(ρ_c, h_c) = \|ρ_c\|/h_c ≤ m̂` ⟹ `τ\|ρ_c\| ≤ max_motion·h_c` (for `m̂ = weightedQuantile 0.9 (m_c) (ℓ_c)` that is at least nine tenths of the length of the cut, `cut-weighted-quantile`); `τ·m ≤ max_motion` for every `m ≤ m̂`, `m̂ > 0`; `robustStep(max(m₁, m₂)) = min(robustStep(m₁), robustStep(m₂))` for `max_motion ≥ 0` (the larger of the two `m̂` is the smaller of the two steps). Extremes: `robustStep(ρ_max/h_min) = adaptiveStep(h_min, ρ_max)` for `h_min > 0`; with `\|ρ_c\| ≤ ρ_max`, `0 ≤ h_min ≤ h_c`: `τ\|ρ_c\| ≤ max_motion·h_c` at every cut point. | `stepCap_le_max`, `stepCap_mul_le`, `stepCap_mul_abs_le`, `le_stepCap`, `stepCap_nonneg`, `adaptiveStep_le_max`, `adaptiveStep_mul_omega_le`, `adaptiveStep_mul_rho_le`, `adaptiveStep_skip_both`, `adaptiveStep_skip_omega`, `adaptiveStep_skip_rho`, `adaptiveStep_full`, `le_adaptiveStep`, `adaptiveStep_nonneg`, `robustStep_eq_adaptiveStep`, `robustStep_spec`, `robustStep_capped`, `robustStep_mul_mhat_le`, `robustStep_motion_le'`, `robustStep_mul_le_edge`, `robustStep_mul_le_of_le`, `robustStep_max`, `robustStep_extremes`, `extremesStep_motion_le` |
| `cut-length-step-on-cut` | Cut points `C`, `ℓ_p > 0`, symmetric couplings `≥ 0`. `fluxMatrix ℓ k` has off-diagonal entries `−k(p,q) ≤ 0`, row sums `ℓ_p`, is symmetric positive definite, and `fluxMatrix·d ≥ 0 ⟹ d ≥ 0`. Force system `fluxMatrix ℓ κ·λ⁰ = Θ` and rate system `fluxMatrix ℓ (τ·m)·r = ℓ(s + λ⁰)` (`m_σ = w̄_σ/\|σ\|`, `τ ≥ 0`; the oracle's `rateCoupling`, `rateRhs`): each has exactly one solution; `Σℓλ⁰ = ΣΘ`; `Σ_S ℓ r = Σ_S ℓ(s+λ⁰)` and `lo ≤ s+λ⁰ ≤ hi` on `S` ⟹ `lo ≤ r ≤ hi` on `S`, for every piece `S` of the cut; `s + λ⁰ = 0 ⟹ r = 0` for every `τ`; without coupling `r = s + λ⁰`. `rippleFactor = 1 + τ(ωk − wk²)/(1 + τwk²) = (1 + τωk)/(1 + τwk²)`: in `(0,1)` for `τ, w > 0`, `ω ≥ 0`, `k > ω/w`; `> 1` for `0 < k < ω/w`; `< 0` for `(τ, ω, w, k) = (1, −10, 1, 1)`. Turning: for a valid manifold mesh (the cut mesh of a region) `Σ_{v∈∂}(π − β_v) = 2πχ − Σ_{interior}Ω_v` (`β_v` the angle sum), `= 2π` for a flat disk; with `Θ_p = lengthForce(w, π, β_p, segs)` and constant `w`: `Σℓλ⁰ = w(2πχ − ΣΩ)`. `angleFromCot(cornerCot a b c) = cornerAngle a b c`; `footFraction` of the cotangents at `p`, `q` is the parameter of the foot of the height (the vector to the third vertex is orthogonal to `pq`); `normalDerivativeOfLinear` of an affine `w` is `⟨∇w, f − z⟩/h`. | `fluxMatrix_apply_of_ne_nonpos`, `fluxMatrix_row_sum`, `fluxMatrix_posDef`, `fluxMatrix_monotone`, `onCut_force_existsUnique`, `onCut_rate_existsUnique`, `onCut_force_sum`, `onCut_rate_sum_piece`, `onCut_rate_bounds_piece`, `onCut_rest`, `onCut_rate_nodal`, `onCut_rate_oracle`, `ripple_step`, `rippleFactor_eq`, `rippleFactor_mem_Ioo`, `one_lt_rippleFactor`, `rippleFactor_neg_example`, `TriMesh.sum_turningAngle`, `TriMesh.sum_turningAngle_flat_disk`, `TriMesh.onCut_lengthForce_sum`, `TriMesh.onCut_force_sum_flat_disk`, `angleFromCot_cornerCot`, `footFraction_orthogonal`, `normalDerivativeOfLinear_eq` |
| `cut-extension-local` | `wmean s c x = Σ_s c_j x_j/Σ_s c_j` with `c ≥ 0`, `Σc > 0` lies between every lower and upper bound of the `x_j` with `c_j ≠ 0`. `σ_i = wmean(ℓ_pβ_ip; r_p)` (`localExtension`; the oracle's `weightedMean` of the two lists) lies between the smallest and largest rate at the cut points of `i`; the outward rule `Σ(σ_j/ℓ_ij)/Σ(1/ℓ_ij)` (the oracle's `ringMean`) between the smallest and largest `σ_j`; by induction over the rings every vertex **that a ring reaches** has `lo ≤ σ ≤ hi`. If `ℓ, β ≥ 0` and `Σ_i β_ip = 1` for every `p`: `Σ_p ℓ_p Σ_i β_ip σ_i = Σ_p ℓ_p r_p`. | `wmean_mem`, `localExtension_mem`, `outwardExtension_mem`, `extension_mem_of_rings`, `localExtension_conservation`, `oracle_weightedMean_eq`, `oracle_weightedMean_mem`, `oracle_ringMean_mem`; for any two lists `weightedMean ws vs = Σ_zip w·v / Σ_zip w` and lies between the bounds of the zipped values when the zipped weights are `≥ 0` with positive sum (`weightedMean_eq`, `weightedMean_mem`) |
| `cut-gradient-norm` | Triangle `a, p = a + e₁, q = a + e₂` (nondegenerate) in a real inner product space, `g = αe₁ + βe₂` with `⟨g, e₁⟩ = ⟨g, e₂⟩ = −φ_a` (the gradient of the affine function with value `φ_a` at `a` that vanishes at `p` and `q`): `‖g‖ = gradientNormLone(φ_a, \|pq\|, area) = \|φ_a\|·\|pq\|/(2·area)`, area by Heron. `ḡ = meanGradient(ΣL, Σ L/g)`: for `L ≥ 0`, `ΣL > 0`, `0 < lo ≤ g_σ ≤ hi`: `lo ≤ ḡ ≤ hi`; `Σ L_σ·(ḡ/g_σ) = ΣL`; `ḡ = 1` if all `g_σ = 1`. `gradientNormLone(0, \|pq\|, area) = 0`: a face whose lone vertex has `φ_a = 0` has slope 0 (SPEC §3.9, D23a leaves it out of the mean; the bounds above are for the faces that count, `g_σ > 0`). **Which faces count** (SPEC §3.9, D23, D30): every statement on `ḡ` is over an arbitrary finite set `s` of counted faces, so leaving out the faces with a moved crossing changes the set, not the theorems. For the set the code uses, `s.filter (¬moved ∧ g > 0)`, the slopes are positive by the choice of the set, and with `L ≥ 0`, `ΣL > 0` on it `lo ≤ ḡ ≤ hi` for bounds of the counted slopes. `pieceSlope s L g` (own; `meanGradient` if that is `> 0`, else `1`) is the scale the flow uses: it is `> 0`; `= 1` for `s = ∅`, in particular when every face has a moved crossing or the slope `0`; `= meanGradient` when `L ≥ 0`, `ΣL > 0`, `g > 0` on `s`. Why such a face must go: for `σ₀ ∈ s` with `L_{σ₀} > 0` (and `L ≥ 0`, `g > 0` on `s`), `ḡ ≤ (ΣL/L_{σ₀})·g_{σ₀}`. | `norm_grad_eq_gradientNormLone`, `norm_grad_eq_gradNorm`, `meanGradient_mem`, `oracle_meanGradient_mem`, `meanGradient_normalizes`, `meanGradient_one`, `gradientNormLone_zero`, `meanGradient_le_of_face`, `pieceSlope_pos`, `pieceSlope_empty`, `pieceSlope_of_no_counted`, `pieceSlope_eq_meanGradient`, `pieceSlope_counted_mem` |
| `cut-roughening-rate` | `rippleRate ω w k = ωk − wk²`; `roughRate(Q₀, Q₁) = 2(Q₀ + Q₁)`; for the Dirichlet energy `roughQ = w_r d²`, so `ω ≥ 0`. For `w, k > 0`: `rate > 0 ⟺ k < ω/w`; `rate ≤ ω²/(4w)`, attained at `k = ω/(2w)`. With `α_L = lengthWeightForRipples(ω, ℓ) = ωℓ/(2π)`: rate `= 0` at `k = 2π/ℓ`, `> 0` for `0 < k < 2π/ℓ`, `< 0` for `k > 2π/ℓ` (`ω, ℓ > 0`). Ingredients of the cited linearization: `d/dδ [Q₀(1−kδ)² − Q₁(1+kδ)²]` at `0` is `−roughRate·k`; `e^{−ky}cos(kx)` is harmonic and `−∂_y` at `y = 0` is `k` times its value; `sin ψ·(tan(ψ/2)^m)′ = m·tan(ψ/2)^m`, derivative `m` at `ψ = π/2`; `(tan(ψ/2))′ = 1` at `π/2`. Sphere formula `(m−1)(4/R³ − α(m+1)/R²)` (cited): `= 0` for `m = 1`; `= rippleRate(4/R², α, m/R) − (4/R³ − α/R²)`; `< 0` for all integers `m ≥ 2` ⟺ `α > 4/(3R)`; `4.72 < 4/(3R) < 4.73` for `R = 1/√(4π)`. **Rate of the state** (SPEC §1, D22): `ω = weightedQuantile 0.9 (ω_c) (ℓ_c)` (oracle), `ℓ_c ≥ 0`, `Σℓ > 0`: `ω` is one of the `ω_c`; `Σ{ℓ_c : ω_c ≤ ω} ≥ 0.9·Σℓ`; `Σ{ℓ_c : ω_c < ω} < 0.9·Σℓ`; `ω ≥ 0` if all `ω_c ≥ 0` (any weights). `rippleRate` is nondecreasing in `ω` for `k ≥ 0`. At a cut point with `ω_c ≤ ω`: `rate(ω_c, w, k) < 0` for `w, k > 0`, `k > ω/w`; with `w = lengthWeightForRipples(ω, ℓ)`, `ω, ℓ > 0`: for `k > 2π/ℓ`. At a cut point with `ω_c ≥ ω`: `rate > 0` for `0 < k < 2π/ℓ`. For any `ω_c`: `rate < 0` for `k > (ω_c/ω)·2π/ℓ`, `k > 0`. `lengthWeightForRipples(ω, ℓ) ≤ 0` for `ω ≤ 0`, `ℓ ≥ 0`; `rate(ω, w, k) < 0` for `ω < 0`, `w ≥ 0`, `k > 0`. | `rippleRate_pos_iff`, `rippleRate_le_max`, `rippleRate_cutoff`, `rippleRate_pos_of_longer`, `rippleRate_neg_of_shorter`, `roughQ_dirichlet`, `roughRate_dirichlet_nonneg`, `hasDerivAt_ripple_density`, `halfPlane_mode_harmonic`, `halfPlane_mode_dtn`, `sphere_mode_ode`, `sphere_mode_dtn`, `capU_second_derivative_equator`, `sphereRippleRate_eq_planar`, `sphere_equator_stable_iff`, `sphere_unit_area_threshold`, `roughening_quantile_spec`, `roughening_quantile_above`, `roughening_quantile_nonneg`, `rippleRate_mono_omega`, `rippleRate_neg_of_le_of_lt`, `rippleRate_neg_of_le`, `rippleRate_pos_of_ge`, `rippleRate_neg_of_lt_longer`, `lengthWeightForRipples_nonpos`, `rippleRate_neg_of_neg` |
| `cut-weighted-quantile` | Pairs `(x_c, ℓ_c)` = `xs.zip ls`; this package's `totalWeight = Σℓ`, `weightLE v = Σ{ℓ_c : x_c ≤ v}`, `weightLT v = Σ{ℓ_c : x_c < v}`, and `IsWeightedQuantile pairs q v`: `v` is one of the `x_c`, `q·Σℓ ≤ weightLE v`, and every `x_c < v` has `weightLE x_c < q·Σℓ` (at most one `v` has this property). For `ℓ_c ≥ 0`, `Σℓ > 0`, `0 < q ≤ 1` the oracle's `weightedQuantile q xs ls` has it. Hence: it is one of the values; `lo ≤ · ≤ hi` for all bounds `lo ≤ x_c ≤ hi`; `weightLE ≥ q·Σℓ`; `weightLT < q·Σℓ`; `weightedQuantile q ≤ weightedQuantile q′` for `0 < q ≤ q′ ≤ 1`; the value `c` if all `x_c = c`; it equals `quantileScan` of **every** ascending arrangement of the pairs (the order of equal values does not matter) and this package's `Cuts.weightedQuantile` (Mathlib's insertion sort, then the scan). For `ℓ_c > 0` and a nonempty list: `weightedQuantile 1` is one of the values and `≥` every `x_c`. `weightedQuantile q [] ls = 0`. For any weights and any `q` the result is `0` or one of the values, so `≥ 0` if all `x_c ≥ 0`. Oracle's sort: `sortAscending l` is a permutation of `l`, ascending; `sortAscending [(1,10),(1,20)] = [(1,20),(1,10)]`. | `oracle_weightedQuantile_spec`, `IsWeightedQuantile.unique`, `oracle_weightedQuantile_bounds`, `oracle_le_weightLE_weightedQuantile`, `oracle_weightLT_weightedQuantile_lt`, `oracle_weightedQuantile_mono`, `oracle_weightedQuantile_const`, `isWeightedQuantile_scan`, `oracle_weightedQuantile_eq_scan`, `oracle_weightedQuantile_eq`, `oracle_weightedQuantile_one`, `oracle_weightedQuantile_nil`, `oracle_weightedQuantile_mem_or`, `oracle_weightedQuantile_nonneg`, `sortAscending_perm`, `sortAscending_pairwise`, `sortAscending_reverses_ties` |
| `cut-step-local` | **At the cut points** (the length step on the cut). `f = localLimitFactor(max_motion, h_c, τ, v_c)` (oracle): `min(1, max_motion·h_c/(τv_c))` if `τv_c > 0`, else `1`; for `max_motion·h_c ≥ 0` it is this package's `localFactor(allowed, moved)` (`allowed/moved` if `allowed < moved`, else `1`: the form of the Rust code) at `allowed = max_motion·h_c`, `moved = τv_c`. `max_motion·h_c > 0` ⟹ `0 < f ≤ 1`. `τv_c ≤ max_motion·h_c` ⟹ `f = 1`. `max_motion·h_c ≥ 0`, `τ ≥ 0`, `\|r_c\| ≤ v_c` ⟹ `τ\|r_c f\| ≤ max_motion·h_c`; in particular for `v_c = localSpeed(r_c, ω_c, ω) = max(\|r_c\|, max(0, ω_c − ω)/2)`, which is `≥ \|r_c\|`, `≥ (ω_c − ω)/2`, and `= \|r_c\|` for `ω_c ≤ ω`. `τ > 0`, `τv_c > max_motion·h_c ≥ 0` ⟹ `τ·f = max_motion·h_c/v_c`. For `max_motion·h_c > 0`: `r f > 0 ⟺ r > 0`, `r f < 0 ⟺ r < 0`, `r f = 0 ⟺ r = 0`, `\|r f\| ≤ \|r\|`. The two cases agree (`τ ≥ 0`, `h_c > 0`): `τm̂ ≤ max_motion` and `v_c/h_c ≤ m̂` ⟹ `τv_c ≤ max_motion·h_c` and `f = 1`; so `(1 if v_c/h_c ≤ m̂, else f) = f`; the robust step satisfies `τm̂ ≤ max_motion` (`τ_max, time_scale, max_motion ≥ 0`) and is `> 0` for positive parameters; `0/h_c ≤ m̂` for a quantile `m̂` of nonnegative speeds. Descent, with `dE(σ) = Σ_c g_c σ_c ℓ_c` (`cut-speed`), `ℓ ≥ 0`: a gradient density scaled point by point, `σ_c = −g_c f_c`, `f_c ≥ 0`, has `dE = −Σ f_c g_c² ℓ_c ≤ 0`; in the flow that is the limited rate when there is no coupling along the cut (then `r = s + λ⁰`): `dE(−r·f) ≤ 0`. For the rate `r` of the system with coupling (`cut-length-step-on-cut`), `g = s + λ⁰`: `dE(−r) = −rᵀ(diag ℓ + τK)r ≤ 0`; but `dE(−r·f) > 0` in the example `ℓ = (9, 1)`, `h = (1, 31/110)`, `g = (1, −20)`, `m = 36`, `Capped{1, 1/4}` (then `m̂ = 1`, `τ = 1/4`, `r = (−10/11, −31/11)`, `f = (1, 1/10)`, `dE = −530/11` unlimited and `+28/11` limited). `\|cutSpeed(Q₀, Q₁)\| ≤ roughRate(Q₀, Q₁)/2` for `Q₀, Q₁ ≥ 0`. **At the vertices** (the length step on the mesh, D26), all about the oracle's definitions: `vertexLimitFactor = localLimitFactor`; `limitedValue(φ, φ′, f) = φ + f(φ′ − φ)`; `vertexSpeedCapped(Δ, τ) = \|Δ\|/τ = \|Δ/τ\|` and `vertexSpeedAdaptive(Δ, τ, ω_i, ω) = max(\|Δ\|/τ, max(0, ω_i − ω)/2) = localSpeed(Δ/τ, ω_i, ω)` for `τ > 0`; the adaptive speed is the capped one for `ω_i ≤ ω` (`τ ≥ 0`), and is the same for `max(ω, M)` as for `M` in place of `ω_i`. With `φ″ = limitedValue(φ, φ′, vertexLimitFactor(max_motion, h_i, τ, v_i))`: `max_motion·h_i ≥ 0`, `τ > 0`, `v_i ≥ \|φ′ − φ\|/τ` (both speeds) ⟹ `\|φ″ − φ\| ≤ max_motion·h_i`. `max_motion·h_i > 0` ⟹ `φ″ − φ > 0 ⟺ φ′ − φ > 0`, `< 0 ⟺ < 0`, `φ″ = φ ⟺ φ′ = φ`, `\|φ″ − φ\| ≤ \|φ′ − φ\|`, `min(φ, φ′) ≤ φ″ ≤ max(φ, φ′)`. `τv_i ≤ max_motion·h_i` ⟹ `φ″ = φ′`; capped (`τ > 0`): `\|φ′ − φ\| ≤ max_motion·h_i` ⟹ `φ″ = φ′`; adaptive: if also `τ(ω_i − ω)/2 ≤ max_motion·h_i`. `v_i = 0` ⟹ factor `1`. Capped rule (`max_motion·h_i ≥ 0`, `τ > 0`): `φ″ = φ + max(−max_motion·h_i, min(max_motion·h_i, φ′ − φ))`, independent of `τ`; a limited vertex changes by exactly `max_motion·h_i`. **`N` regions**: `Δ = vertexChange(φ_b′, φ_b, φ_a′, φ_a) = ½[(φ_b′ − φ_b) − (φ_a′ − φ_a)]`; with every function taking the factor `f`, `vertexChange(φ_b″, φ_b, φ_a″, φ_a) = Δ·f`, so under the hypotheses above (`v ≥ \|Δ\|/τ`, both speeds) `\|vertexChange(φ_b″, φ_b, φ_a″, φ_a)\| ≤ max_motion·h_i`, with its sign kept and `0` only if `Δ = 0`; `φ_k″ − φ_l″ = (φ_k − φ_l) + f·[(φ_k′ − φ_k) − (φ_l′ − φ_l)]`; `\|φ_k″ − φ_k\| ≤ \|φ_k′ − φ_k\|` for `0 ≤ f ≤ 1`. On `[φ, −φ]` with the whole step `[φ′, −φ′]` and `a ≠ b` in `Fin 2`: `Δ = −(φ′ − φ)` (`a = 0`) or `φ′ − φ` (`a = 1`), both speeds equal those of the two-region rule for `φ′ − φ`, and for every `f` the limited state is `[limitedValue(φ, φ′, f), −limitedValue(φ, φ′, f)]`. | at the cut points: `localLimitFactor_def`, `localLimitFactor_eq`, `localLimitFactor_mem_Ioc`, `localLimitFactor_eq_one`, `localLimit_motion_le`, `localLimit_motion_le_adaptive`, `abs_le_localSpeed`, `localSpeed_of_le`, `localLimit_step_eq`, `localLimit_sign`, `localLimit_cases_agree`, `localLimit_twoCase_eq`, `robustStep_mul_mhat_le`, `robustStep_pos`, `localLimitFactor_eq_one_of_speed_le`, `localLimit_first_case_of_zero`, `firstVariation_localLimit`, `firstVariation_limited`, `onCut_limited_descent_of_no_coupling`, `onCut_rate_descent`, `limited_rate_not_descent_example`, `abs_cutSpeed_le_half_roughRate`; on `localFactor`: `localFactor_eq_min`, `localFactor_eq_one`, `localFactor_mem_Ioc`, `limited_motion_le`, `limited_sign`; at the vertices: `vertexLimitFactor_eq`, `limitedValue_eq`, `vertexSpeedCapped_eq_abs_rate`, `vertexSpeedAdaptive_eq_localSpeed`, `vertexSpeedAdaptive_of_le`, `vertexSpeedAdaptive_max_omega`, `vertexLimit_change_le`, `vertexLimit_change_le_capped`, `vertexLimit_change_le_adaptive`, `vertexLimit_sign`, `vertexLimit_between`, `vertexLimit_unchanged`, `vertexLimit_unchanged_capped`, `vertexLimit_unchanged_adaptive`, `vertexLimitFactor_zero_speed`, `vertexLimit_capped_clamp`, `vertexLimit_capped_limited`; `N` regions: `vertexChange_limited`, `vertexLimit_regions_change_le`, `vertexLimit_regions_change_le_capped`, `vertexLimit_regions_change_le_adaptive`, `vertexLimit_regions_sign`, `limitedValue_sub_limitedValue`, `abs_limitedValue_sub_le`, `vertexChange_two_lower`, `vertexChange_two_upper`, `abs_vertexChange_two`, `vertexLimit_two_regions` |
| `cut-edge-crossing-regions` | `IsRegion x k`: `x_k ≤ x_l` for all `l` and `x_k < x_l` for `l < k` (unique; exists). For `i` the region of `a`, `j ≠ i` the region of `b`: `ψ_a = φ_i(a) − φ_j(a) ≤ 0 ≤ ψ_b`, `ψ_b > 0` if `i < j`, `ψ_a < 0` if `j < i`, so `ψ_a ≠ ψ_b`; `t = crossingParamRegions = crossingParam(ψ_a, ψ_b)` satisfies `(1−t)φ_i(a) + tφ_i(b) = (1−t)φ_j(a) + tφ_j(b)`, `0 ≤ t ≤ 1`, `t < 1` if `i < j`, `t > 0` if `j < i`. On `[φ, −φ]`: region 0 ⟺ `φ ≤ 0`, region 1 ⟺ `φ > 0`, and `t = crossingParam(φ_a, φ_b)` for either orientation. | `IsRegion.unique`, `exists_isRegion`, `region_differences`, `crossingParamRegions_interp`, `crossingParamRegions_mem`, `isRegion_two_zero_iff`, `isRegion_two_one_iff`, `crossingParamRegions_two` |
| `cut-triple-point` | `m_k = triplePointCoord(s_k, s_{k−1}) = ((1 − s_k) + s_{k−1})/3`: `m₀ + m₁ + m₂ = 1`; `m_k > 0` if `s_k < 1`, `s_{k−1} > 0`; `(p₀+p₁+p₂)/3 = Σ m_k c_k` for `p_k = (1−s_k)c_k + s_k c_{k+1}`. Barycentric determinants: `(c_k, p_k, m)`: `s_k m_{k+2}`; `(c_k, m, p_{k−1})`: `m_{k+1}(1 − s_{k−1})`; the six are positive for `0 < s < 1` and sum to 1. Convexity of the quad `(c_k, p_k, m, p_{k−1})`: `(c_k, p_k, p_{k−1})` has `s_k(1 − s_{k−1}) > 0` and `(p_k, m, p_{k−1})` has `((1−s₀)(1−s₁)(1−s₂) + s₀s₁s₂)/3 > 0`, a third of `det(p₀, p₁, p₂)`. | `triplePoint_coord_sum`, `triplePoint_coord_pos`, `triplePoint_eq_centroid`, `tripleFace_dets`, `tripleFace_dets_pos`, `tripleFace_convex` |
| `cut-region-projection` | `IsTwoSmallest x i₀ i₁`: `i₀ ≠ i₁`, `x_{i₀} ≤ x_{i₁} ≤ x_k` for `k ≠ i₀` (exists for at least two functions). With `m = projectionShift(x_{i₀}, x_{i₁})`, `y_k = projectedValue(x_k, m)`, `d = (x_{i₁} − x_{i₀})/2 ≥ 0`: `y_{i₀} = −d`, `y_{i₁} = d`, `y_k ≥ d` for `k ≠ i₀`; if `x_{i₀} < x_{i₁}`: `y_k < 0 ⟺ k = i₀`; `y_k − y_l = x_k − x_l`; the region is unchanged; `crossingParamRegions` of an edge is unchanged by the shifts at its two ends; `IsTwoSmallest y i₀ i₁` and its shift is `0`. On `[φ, −φ]`: `m = 0`. | `exists_isTwoSmallest`, `projected_two_smallest`, `projected_rest`, `projected_neg_iff`, `projected_sub`, `isRegion_projected`, `crossingParamRegions_projected`, `projected_isTwoSmallest`, `projectionShift_projected`, `projectionShift_two` |
| `cut-region-extension` | `tripleRate(\|mp\|, \|mp′\|, s_p, s_{p′}) = smoothedSpeed` of the same arguments: between `min` and `max` of the two rates, exact for a rate linear in arclength, odd. `regionRate own e`: `e` for the own function, `−e` for the others; `φ_own − φ_j` changes at `2e`, the others keep their differences. Two regions: `regionRate 0 σ = regionRate 1 (−σ) = [σ, −σ]`; if `σ` is the harmonic extension of `s` then `−σ` is that of `−s`; `explicitStep(−φ, τ, −σ) = −explicitStep(φ, τ, σ)`; if `φ′` solves the length step for `φ̃` then `−φ′` solves it for `−φ̃`; `classMean(−φ) = −classMean φ`. | `tripleRate_eq_smoothedSpeed`, `tripleRate_mem`, `tripleRate_linear`, `tripleRate_neg`, `regionRate_sub`, `regionRate_two`, `harmonicExtension_neg`, `explicitStep_neg`, `lengthStep_neg`, `classMean_neg` |
| `cut-length-step-regions` | The two systems of `cut-length-step-on-cut` over the boundary copies of all regions: exactly one `λ⁰` and one `ρ`; for the boundary nodes `X` of one region (no coupling leaves `X`): `Σ_X ℓλ⁰ = Σ_X Θ`, `Σ_X ℓρ = Σ_X ℓ(s+λ⁰)`, `ρ` between the bounds of `s + λ⁰` on `X`. Two chains with the same `ℓ` and couplings and opposite `s`, `λ⁰` have opposite `ρ`; opposite right-hand sides give opposite solutions. If `β_i + β_j = 2π`: `turningAngle(π, β_j) = −turningAngle(π, β_i)`. If `β₀ + β₁ + β₂ = 2π`: `Σ turningAngle(2π/3, β_k) = 0`; `turningAngle(2π/3, 2π/3) = 0`; `> 0` for `β < 2π/3`. Explicit step: `explicitStep(φ_k, τ, e) − explicitStep(φ_j, τ, e − δ) = (φ_k − φ_j) + τδ`; with opposite rates on the two sides `δ = 2e`, and `explicitStep(−φ, τ, e − 2e) = −explicitStep(φ, τ, e)`. | `regionChain_existsUnique`, `regionChain_sums`, `regionChain_bounds`, `onCut_rates_opposite`, `fluxRecovery_neg`, `curve_loads_opposite`, `junction_loads_sum`, `junction_load_balanced`, `junction_load_pos`, `regionStep_difference`, `wmean_sub_neg`, `regionStep_two` |
| `cut-disk-radial-solution` | `u(r) = K(r² − ρ²)/4`: `u″ + u′/r = K` (`r ≠ 0`), `∂²/∂x² + ∂²/∂y²` of `K(x²+y²−ρ²)/4` is `K`, `u(ρ) = 0`, `u′(ρ)·2πρ = Kπρ²`; `∫₀^ρ u′² 2πr dr = πK²ρ⁴/8` with `ρ`-derivative `u′(ρ)²·2πρ`; `∫₀^ρ u² 2πr dr = πK²ρ⁶/48`; `p(r) = Kr⁴/32 − Kρ²r²/8 + 3Kρ⁴/32`: `p″ + p′/r = 2u`, `p(ρ) = 0`, `p′(ρ) = −Kρ³/8`; `d/dρ ∫₀^ρ u² 2πr dr = (u(ρ)² − u′(ρ)p′(ρ))·2πρ = πK²ρ⁵/8`. | `diskU_ode`, `diskU_laplacian`, `diskU_boundary`, `disk_flux_balance`, `disk_dirichletEnergy`, `hasDerivAt_disk_dirichletEnergy`, `disk_henckyEnergy`, `diskP_ode`, `diskP_boundary`, `deriv_diskP_boundary`, `hasDerivAt_disk_henckyEnergy`, `disk_hencky_derivative_value` |
| `cut-sphere-cap-solution` | `u(ψ) = 2 ln(cos(ψ₀/2)/cos(ψ/2))`: `u(ψ₀) = 0`; `u′(ψ) = tan(ψ/2)` where `cos(ψ₀/2), cos(ψ/2) ≠ 0`; `(sin ψ·u′)′ = sin ψ` on `(−π, π)`; `tan(ψ₀/2)·2π sin ψ₀ = 2π(1 − cos ψ₀)`; for `\|ψ₀\| < π`: `∫₀^{ψ₀} u′² 2π sin ψ dψ = −8π ln cos(ψ₀/2) − 4π sin²(ψ₀/2)` with `ψ₀`-derivative `u′(ψ₀)²·2π sin ψ₀`; two caps: see `cut-speed`; `cutSpeed = 0` at `ψ₀ = π/2`. | `capU_boundary`, `hasDerivAt_capU`, `capU_ode`, `cap_flux_balance`, `capEnergy_eq`, `hasDerivAt_capEnergy`, `hasDerivAt_capEnergy_integral`, `hasDerivAt_twoCaps`, `twoCaps_equator`, `cap_integrand_scale` |
| `cut-curve-shortening-circle` | `ρ(t) = √(ρ₀² − 2αt)`, `ρ₀² − 2αt > 0`: `ρ′(t) = −α/ρ(t)`; `ρ(0) = ρ₀` (`ρ₀ ≥ 0`); `πρ(t)² = πρ₀² − 2παt`; for `φ(t, r) = r − ρ(t)`: `φ_t = α(φ_rr + φ_r/r)` at `r = ρ(t)`; `d/dt [α·2πρ(t)] = −(α/ρ)²·2πρ`. | `hasDerivAt_shrinkRadius`, `shrinkRadius_zero`, `shrinkRadius_area`, `circle_levelSet_flow`, `hasDerivAt_circle_lengthEnergy` |

Cited for the cut flow (`CITED.md`): `cut-yamabe-equation`, `cut-shape-derivative`,
`cut-cea-method`, `cut-level-set-evolution`, `cut-ripple-linearization`,
`cut-network-first-variation`, `cut-weighted-fem-stiffness`, `cut-fast-marching-convergence`.

**The stop rule of a run** (SPEC §4 "Stopping", D27, D29, D31): the ids `cut-stop-energy` (the
Rust `energy_stop`; theorem `energy_cases`), `cut-run-lowest` (`run_outcome`; theorem
`runOutcome_rose`) and `cut-stop-motion` (`motion_criterion`; theorem `motionCriterion_spec`),
in `DdgProofs/Cuts/Stop.lean`; the definitions are this package's own, written after the Rust
functions.

`cut-stop-motion` (D31). `motionCriterion tol moved predicted` is
`tol > 0 ∧ moved ≤ tol ∧ predicted ≤ tol`:

- it implies `moved ≤ tol` (the criterion before D31); for `tol ≤ 0` it never holds; it is kept
  when `moved` and `predicted` get smaller (`motionCriterion_spec`, `motionCriterion_mono`);
- length step on the cut: `predictedMotionOnCut τ A c r = τ · max_{i ∈ A} |Σ_p c_ip r_ip / Σ_p
  c_ip|` (`A` the cut-adjacent vertices, `c_ip = ℓ_p β_ip`, `r_ip` the rate of the cut point `p`
  as the vertex `i` sees it; the maximum over an empty set is `0`). For `τ ≥ 0`, and at every
  `i ∈ A` weights `c_ip ≥ 0` with `Σ_p c_ip > 0`: it is `≤ τ·R` for every `R ≥ 0` with
  `|r_ip| ≤ R` where `c_ip ≠ 0`, in particular `≤ τ · max_p |r_p|`
  (`predictedMotionOnCut_le`, `predictedMotionOnCut_le_max`); it is `≥ 0`; it is `0` when all
  rates are `0`, whatever the weights, for instance at a rest state `s + λ⁰ = 0` of the system
  on the cut (`predictedMotionOnCut_zero`, `predictedMotionOnCut_rest`). With all rates `0` the
  criterion is `tol > 0 ∧ moved ≤ tol`, and holds for `moved = 0` and every `tol > 0`
  (`motionCriterion_zero_rates`, `motionCriterion_rest`);
- length step on the mesh: `predictedMotionOnMesh A φ φ′ = max_{i ∈ A} |φ′(i) − φ(i)|`; it is
  `0` if `φ′ = φ` on `A`; for factors `f_i ∈ [0, 1]`,
  `max_{i ∈ A} |limitedValue(φ_i, φ′_i, f_i) − φ_i| ≤` it: the local limit at the vertices does
  not exceed the prediction (`predictedMotionOnMesh_zero`, `limited_le_predictedMotionOnMesh`).

`cut-stop-energy`, `cut-run-lowest`. With
`t = ε\|E_k\|`, `Eb = E_{k−W}`:

- for `ε ≥ 0` exactly one of `Eb − E_k > t`, `\|Eb − E_k\| ≤ t`, `E_k − Eb > t` holds
  (`energy_cases`); they are exhaustive for every `ε`, and for `ε < 0` the first and the third
  can hold together (`energy_cases_exhaustive`, `energy_cases_overlap_of_neg`);
- `energyVerdict ε Eb E_k` (the chain of comparisons of `energy_stop`) is `goesOn` iff
  `Eb − E_k > t`, `stalled` iff `\|Eb − E_k\| ≤ t`, and for `ε ≥ 0` `rose` iff `E_k − Eb > t`; it
  is not `goesOn` iff `Eb − E_k ≤ t`, the stop condition before D27 (`energyVerdict_goesOn_iff`,
  `energyVerdict_stalled_iff`, `energyVerdict_rose_iff`, `energyVerdict_stops_iff`,
  `energyVerdict_zero_stops_iff`);
- `lowestState E n` (replace the state held only by one of strictly lower energy) returns
  `(E_j, j)` with `j ≤ n`, `E_j ≤ E_i` for all `i ≤ n` and `E_j < E_i` for all `i < j`, and `j`
  is the only such state (`lowestState_spec`, `lowestState_unique`);
- `runOutcome` (`run_outcome`): without a vanished cut, after `n + 1` steps, the lowest energy
  and its state are `lowestState E (n + 1)`, and the state returned is that one for
  `EnergyRose` and `n + 1` otherwise; with no step it is the start state; with a vanished cut
  the last state does not count and the run stays in it (`runOutcome_eq_lowestState`,
  `runOutcome_rose`, `runOutcome_not_rose`, `runOutcome_no_steps`, `runOutcome_vanished`);
- after `EnergyRose` (`ε ≥ 0`, `k ≤ n`) the lowest energy of the run is `< E_k`
  (`lowestState_lt_of_rose`).

### Navigating Intrinsic Triangulations (`crates/ddg-intrinsic`)

Notation of the oracle (`lean/oracle/DdgOracle/Formulas/Intrinsic.lean`, at `ℝ`; the theorems are
about its definitions unless a definition of this package is named). The *diamond* of an edge
`ij` of length `a`: the triangles `ijk` (`b = ℓ_ik`, `c = ℓ_jk`) and `jil` (`d = ℓ_il`,
`e = ℓ_jl`); `x_k = diamondX a b c`, `h_k = diamondH a b c`, `x_l = diamondX a d e`,
`h_l = diamondH a d e`; the *standard layout* `i = (0,0)`, `j = (a,0)`, `k = (x_k, h_k)` (this
package's `apexK a b c`), `l = (x_l, −h_l)` (`apexL a d e`); `f = flippedLengthFromLengths a b c
d e`. `T(x,y,z)` is `IsTriangle x y z` (the strict triangle inequalities). Plane points are pairs
of reals; `|pq| = dist2 p q` and `[p,q,r] = sda p q r` (twice the signed area) are those of
`DdgProofs/Bff/Quality.lean`. Angles: `θ_k = cornerAngle a b c`, `θ_l = cornerAngle a d e`
(opposite the edge); at `i`: `α_i = cornerAngle c a b`, `β_i = cornerAngle e a d`; at `j`:
`α_j = cornerAngle b c a`, `β_j = cornerAngle d e a`. `Convex` is this package's
`DiamondConvex a b c d e`: the oracle's `flipConvex` of the standard layout. After a flip the
edge `kl` has the diamond `(f, c, e, b, d)` (apexes `j` and `i`). All theorems are in
`namespace DdgProofs.Intrinsic`.

| id | exact statement | theorems |
|---|---|---|
| `it-unfold-diamond` | `T(a,b,c)`, `T(a,d,e)`: `h_k, h_l > 0`, and in the standard layout `\|ij\| = a`, `\|ik\| = b`, `\|jk\| = c`, `\|il\| = d`, `\|jl\| = e`; `[i,j,k] = a·h_k`, `[i,j,l] = −a·h_l` (the apexes on opposite sides). Uniqueness: for any plane points with `[i,j,k] > 0 > [i,j,l]` there are `cs, sn` with `cs² + sn² = 1` such that `i, j, k, l` are the images of `(0,0)`, `(\|ij\|, 0)`, `apexK(\|ij\|, \|ik\|, \|jk\|)`, `apexL(\|ij\|, \|il\|, \|jl\|)` under `rigid cs sn i : p ↦ i + (cs·p₁ − sn·p₂, sn·p₁ + cs·p₂)`; these maps keep `dist2` and `sda`. In every layout `diamondX(\|ij\|, \|ik\|, \|jk\|) = ⟨j−i, k−i⟩/\|ij\|` (`\|ij\| > 0`) and `diamondH(\|ij\|, \|ik\|, \|jk\|) = \|[i,j,k]\|/\|ij\|`. `diamondX a b c = layoutX a c b`, `diamondH a b c = layoutY a (triArea a c b)` (the reference layout of BFF). A plane triangle with `[p,q,r] ≠ 0` has `T(\|qr\|, \|rp\|, \|pq\|)`. | `diamond_isometric`, `diamond_sda`, `diamond_unique`, `dist2_rigid`, `sda_rigid`, `diamondX_realized`, `diamondH_realized`, `diamondX_eq_layoutX`, `diamondH_eq_layoutY`, `isTriangle_of_sda_ne` |
| `it-unfold-apex` | `C = unfoldApex A (B − A) ℓ_ab ℓ_bc ℓ_ca`. If `\|AB\| = ℓ_ab` and `T(ℓ_ab, ℓ_bc, ℓ_ca)`: `\|AC\| = ℓ_ca`, `\|BC\| = ℓ_bc`, `[A,B,C] = 2·triArea ℓ_ab ℓ_bc ℓ_ca > 0`. With `T` alone, whatever `\|AB\|`: `\|AC\|² = (\|AB\|/ℓ_ab)²·ℓ_ca²`, `\|BC\|² = (\|AB\|/ℓ_ab)²·ℓ_bc²` (the given triangle, scaled). Uniqueness: `[A,B,C′] > 0 ⟹ unfoldApex A (B − A) \|AB\| \|BC′\| \|C′A\| = C′`. `unfoldApex (0,0) (a,0) a c b = (x_k, h_k)` for `a ≠ 0`. The exit rule: `apexAboveAxis y = true ⟺ y > 0`; with this package's `MeetsAxis p q` (the closed segment has a point with `y = 0`; `⟺` its ends are not strictly on one side) and `y_up > 0 > y_down`: `MeetsAxis down apex ⟺ y_apex ≥ 0`; `MeetsAxis apex up ⟺ y_apex ≤ 0`; one of the two holds; both hold `⟺ y_apex = 0`; `apexAboveAxis y_apex = true ⟹ MeetsAxis down apex ∧ ¬MeetsAxis apex up`. | `unfoldApex_isometric`, `unfoldApex_similar`, `unfoldApex_unique`, `apex_unique`, `unfoldApex_eq_rigid`, `unfoldApex_layout`, `apexAboveAxis_iff`, `meetsAxis_iff`, `exit_cases` |
| `it-flip-length` | `flippedLength x_k h_k x_l h_l = \|(x_k, h_k)(x_l, −h_l)\|`. For plane points with `[i,j,k] > 0 > [i,j,l]`: `flippedLengthFromLengths \|ij\| \|ik\| \|jk\| \|il\| \|jl\| = \|kl\|`. `T(a,b,c)`, `T(a,d,e)`: `f² = b² + d² − 2bd·cos(α_i + β_i)` (the paper's Fig. 9 and Alg. 5) `= c² + e² − 2ce·cos(α_j + β_j) = b² + d² − [(a²+b²−c²)(a²+d²−e²) − 16·A_ijk·A_jil]/(2a²)`; `flippedLengthFromLengths a c b e d = f` (`a ≠ 0`; mirrored) and `flippedLengthFromLengths a d e b c = f` (the triangles exchanged). If also `Convex`: `flippedLengthFromLengths f c e b d = a` (flipping back); `T(f,b,d)` and `T(f,e,c)` (the new triangles `ilk`, `jkl`); `triArea f b d + triArea f e c = triArea a b c + triArea a d e`. | `flippedLength_eq_dist2`, `flippedLengthFromLengths_realized`, `flippedLengthFromLengths_sq`, `flippedLengthFromLengths_sq'`, `flippedLengthFromLengths_sq_lengths`, `flippedLengthFromLengths_mirror`, `flippedLengthFromLengths_comm`, `flip_back`, `flip_isTriangle`, `flip_area` |
| `it-flip-convexity` | `flipConvex a x_k h_k x_l h_l = true ⟺ flipArea1 > 0 ∧ flipArea2 > 0`, with `flipArea1 = x_k h_l + x_l h_k = [i,l,k]` and `flipArea2 = (a − x_k) h_l + (a − x_l) h_k = [j,k,l]` in the standard layout, and in every layout with `[i,j,k] > 0 > [i,j,l]` (standard values of its five distances). For `a, h_k, h_l > 0`: `flipConvex = true ⟺ ∃ s, t ∈ (0,1): (s·a, 0) = k + t·(l − k)` (the diagonals cross in their interiors) `⟺ [i,l,j], [l,j,k], [j,k,i], [k,i,l] > 0` (the quadrilateral `i, l, j, k` is strictly convex). `T(a,b,c)`, `T(a,d,e)`: `flipArea1 = b·d·sin(α_i + β_i)`, `flipArea2 = c·e·sin(α_j + β_j)`, `Convex ⟺ α_i + β_i < π ∧ α_j + β_j < π`, and `Convex ⟹ DiamondConvex f c e b d`. | `flipConvex_iff`, `flipArea1_eq_sda`, `flipArea2_eq_sda`, `flipArea_realized`, `flipConvex_iff_diagonals_cross`, `flipConvex_iff_turns`, `flipArea1_eq_sin`, `flipArea2_eq_sin`, `diamondConvex_iff_angles`, `flip_convex` |
| `it-delaunay-condition` | `delaunayExcess θ_k θ_l = θ_k + θ_l − π = (θ_k − π) + θ_l`; `delaunayHolds tol θ_k θ_l = true ⟺ delaunayExcess θ_k θ_l ≤ tol`; `oppositeAngles a b c d e = (θ_k, θ_l)` for `T(a,b,c)`, `T(a,d,e)`. Then: `delaunayHolds 0 (oppositeAngles a b c d e) = true ⟺ cotanWeightFromLengths a b c d e ≥ 0 ⟺ θ_k + θ_l ≤ π`, and `θ_k + θ_l > π ⟺ cotanWeightFromLengths < 0`. `θ_k + θ_l ≥ π ⟹ Convex` (a violating edge, and a tie, can be flipped). Two corners of one triangle: `delaunayExcess (cornerAngle a b c) (cornerAngle b c a) < 0`. `(α_i + β_i) + (α_j + β_j) = 2π − (θ_k + θ_l)`. `Convex ⟹ delaunayExcess (cornerAngle f c e) (cornerAngle f b d) = −delaunayExcess θ_k θ_l`. `θ_k + θ_l > π ⟹ Convex`, `delaunayExcess (oppositeAngles f c e b d) < 0` and `cotanWeightFromLengths f c e b d > 0`. `θ_k + θ_l = π ⟹ Convex` and the excess after the flip is `0`. | `delaunayExcess_eq`, `delaunayExcess_eq'`, `delaunayHolds_iff`, `oppositeAngles_eq`, `delaunay_iff_cotanWeight_nonneg`, `angle_sum_le_pi_iff_cotanWeight_nonneg`, `pi_lt_angle_sum_iff_cotanWeight_neg`, `diamondConvex_of_pi_le`, `diamondConvex_of_not_delaunay`, `delaunay_of_same_triangle`, `quad_angle_sum`, `flip_delaunayExcess`, `flip_of_not_delaunay`, `flip_of_tie` |
| `it-cotan-weight` | `cotanWeight x y = x + y`; `cotanWeightFromLengths a b c d e = cornerCot a b c + cornerCot a d e`, which for `T(a,b,c)`, `T(a,d,e)` is `cot θ_k + cot θ_l`: **twice** the cotan weight `w = ½(cot θ_k + cot θ_l)` of the edge (intrinsic SPEC N6, N17). For `θ_k, θ_l ∈ (0, π)`: `cotanWeight (cot θ_k) (cot θ_l) = sin(θ_k + θ_l)/(sin θ_k·sin θ_l)` with a positive denominator; it is `≥ 0 ⟺ θ_k + θ_l ≤ π`, `> 0 ⟺ θ_k + θ_l < π`, `= 0 ⟺ θ_k + θ_l = π`, and the same three hold for `w = cotanWeight/2`. Link to `cotan-laplacian`: in a `TriMesh` in which only the faces `f₁ ≠ f₂` have the edge `{i,j}`, `L_ij = offDiag i j = −cotanWeight (cotOpp f₁ i j) (cotOpp f₂ i j)/2 = −w`, so `cotanWeight = −2·L_ij` and `L_ij ≤ 0 ⟺ cotanWeight ≥ 0`. `cotanWeightFromLengths a b c d e < 0 ⟹ Convex`. | `cotanWeightFromLengths_eq`, `cotanWeightFromLengths_eq_cot`, `cotanWeight_eq_sin_div`, `cotanWeight_nonneg_iff`, `cotanWeight_pos_iff`, `cotanWeight_eq_zero_iff`, `half_cotanWeight_sign`, `offDiag_eq_neg_half_cotanWeight`, `cotanWeight_eq_neg_two_mul_offDiag`, `offDiag_nonpos_iff`, `diamondConvex_of_cotanWeight_neg` |
| `it-angle-normalization` | `angleScale false Θ = 2π/Θ`, `angleScale true Θ = π/Θ`. For corner angles `θ_n` with `Σθ ≠ 0`: `Σ_n angleScale false (Σθ)·θ_n = 2π`, `Σ_n angleScale true (Σθ)·θ_n = π`. `angleScale false (2π) = 1 = angleScale true π`; `angleScale b Θ > 0` for `Θ > 0`; `Θ/(2π)·(angleScale false Θ·θ) = θ` for `Θ ≠ 0` (the paper's Alg. 2, line 7). A flip keeps the angle sums at its four vertices (`T(a,b,c)`, `T(a,d,e)`, `Convex`): `cornerAngle f b d = α_i + β_i`, `cornerAngle f c e = α_j + β_j`, `cornerAngle e f c + cornerAngle d f b = θ_k`, `cornerAngle c e f + cornerAngle b d f = θ_l`. | `angleScale_interior`, `angleScale_boundary`, `sum_normalized_interior`, `sum_normalized_boundary`, `angleScale_flat`, `angleScale_pos`, `angleScale_undo`, `flip_angle_sums`, `flip_angle_i`, `flip_angle_j` |
| `it-signpost-update` | `signpostUpdate false φ θ Θ = reduceFullTurn(φ + (2π/Θ)·θ)`, `reduceFullTurn x = x − 2π` if `2π ≤ x`, else `x`; `signpostUpdate true φ θ Θ = φ + (π/Θ)·θ`. `0 ≤ x < 4π ⟹ 0 ≤ reduceFullTurn x < 2π`, and it is `x` or `x − 2π`. For `0 ≤ φ < 2π`, `Θ > 0`, `0 ≤ θ ≤ Θ`: the update is in `[0, 2π)`, equals `φ + (2π/Θ)θ` or that minus `2π`, and has its cosine and sine. For `φ < 2π`, `Θ > 0`, `θ₁, θ₂ ≥ 0`, `θ₁ + θ₂ ≤ Θ`: two updates in a row are the update by `θ₁ + θ₂`. The update by `Θ ≠ 0` is the identity for `φ ≥ 0`. `List.foldl` of the update over corners `θs` (`Θ > 0`, each `≥ 0`, sum `≤ Θ`, start in `[0, 2π)`) is the update by `θs.sum`; with `θs.sum = Θ` it returns the start; if `φ ≥ 0` and `φ + (2π/Θ)·θs.sum < 2π` it is `φ + (2π/Θ)·θs.sum` (the cumulative sums of the paper's Eq. 2). Boundary vertex: the fold is `φ + (π/Θ)·θs.sum`; from `φ = 0` with `θs.sum = Θ ≠ 0` it is `π`. Under a flip (`T(a,b,c)`, `T(a,d,e)`, `Convex`, `φ < 2π`, `Θ > 0`, the old corner `≤ Θ`): the updates by `cornerAngle e f c` and then `cornerAngle d f b` give the update by `θ_k`; by `cornerAngle c e f` and then `cornerAngle b d f`, the update by `θ_l`. | `signpostUpdate_interior`, `signpostUpdate_boundary`, `reduceFullTurn_mem`, `signpostUpdate_mem`, `signpostUpdate_direction`, `signpostUpdate_add`, `signpostUpdate_full_turn`, `signpost_fold`, `signpost_around_interior`, `signpost_fold_no_wrap`, `signpost_fold_boundary`, `signpost_boundary_end`, `flip_signpost_at_k`, `flip_signpost_at_l` |
| `it-ray-edge-intersection` | `t = axisCrossingT y_tail y_head = y_tail/(y_tail − y_head)` (the `crossingParam` of `cut-edge-crossing`), `axisCrossingX x_tail x_head t = x_tail + t·(x_head − x_tail)`. For `y_tail ≠ y_head`: `tail + t·(head − tail) = (axisCrossingX x_tail x_head t, 0)`, and `y_tail + s·(y_head − y_tail) = 0 ⟺ s = t`. For `y_tail ≤ 0 ≤ y_head`, `y_tail < y_head` (and with the roles exchanged): `0 ≤ t ≤ 1`, `t = 0 ⟺ y_tail = 0`, `t = 1 ⟺ y_head = 0`; for `y_tail < 0 < y_head`: `0 < t < 1`. Forward progress, for `y_up > 0 > y_down` and `[up, down, apex] > 0`: if `y_apex > 0`, `x_exit − x_entry = t·[up,down,apex]/(y_up − y_down) > 0` (exit on `down–apex`, `t` its parameter); if `y_apex < 0`, `x_exit − x_entry = (1 − t)·[up,down,apex]/(y_up − y_down) > 0` (exit on `apex–up`). | `axisCrossingT_eq_crossingParam`, `axisCrossing_on_axis`, `axisCrossing_unique`, `axisCrossingT_mem_Icc`, `axisCrossingT_mem_Icc'`, `axisCrossingT_mem_Ioo`, `exit_ahead_above`, `exit_ahead_below` |
| `it-barycentric-coordinates` | `signedArea2` is `sda`; `barycentric a b c p = ([p,b,c], [a,p,c], [a,b,p])/[a,b,c]`; `[p,b,c] + [a,p,c] + [a,b,p] = [a,b,c]`. For `[a,b,c] ≠ 0`: the three coordinates sum to `1`; `b₀·a + b₁·b + b₂·c = p`; and if `u + v + w = 1`, `u·a + v·b + w·c = p` then `barycentric a b c p = (u, v, w)`. The form of the Rust code, `([p,b,c], [p,c,a], [p,a,b])` divided by their sum, is equal. `T(ℓ_ab, ℓ_bc, ℓ_ca)`: `pointFromDistances ℓ_ab ℓ_bc ℓ_ca \|v₀p\| \|v₁p\| \|v₂p\| = p` for every `p`, in the layout `v₀ = (0,0)`, `v₁ = (ℓ_ab, 0)`, `v₂ = apexK ℓ_ab ℓ_ca ℓ_bc`. | `signedArea2_eq_sda`, `barycentric_eq`, `sda_split`, `barycentric_sum`, `barycentric_reproduces`, `barycentric_unique`, `barycentric_eq_normalized`, `pointFromDistances_eq_self` |
| `it-point-distances` | `cornerDistances ℓ₀₁ ℓ₁₂ ℓ₂₀ b₀ b₁ b₂ = (\|v₀P\|, \|v₁P\|, \|v₂P\|)` with `v₀ = (0,0)`, `v₁ = (ℓ₀₁, 0)`, `v₂ = apexK ℓ₀₁ ℓ₂₀ ℓ₁₂` and `P = b₁·v₁ + b₂·v₂` (this package's `layoutPoint`; `b₀` is not used). For `T(ℓ₀₁, ℓ₁₂, ℓ₂₀)` and `b₀ + b₁ + b₂ = 1` it is `(baryLength ℓ₀₁ ℓ₁₂ ℓ₂₀ (b₀−1) b₁ b₂, baryLength … b₀ (b₁−1) b₂, baryLength … b₀ b₁ (b₂−1))`: the paper's Alg. 6, by `baryLength_eq_dist` of `cut-barycentric-length` (used, not reproved, through `dist2 p q = dist (toC p) (toC q)`); for any plane points and coordinates `x`, `y` with sum `1`: `\|Σ y_k v_k  Σ x_k v_k\| = baryLength \|v₀v₁\| \|v₁v₂\| \|v₂v₀\| (x − y)`. Edge split: `edgeSplitLengths ℓ b_i b_j = (ℓ·b_j, ℓ·b_i)`; the two sum to `ℓ` if `b_i + b_j = 1`; they are the first two of `cornerDistances ℓ₀₁ ℓ₁₂ ℓ₂₀ b₀ b₁ 0` (`b₀, b₁ ≥ 0`, `b₀ + b₁ = 1`); `baryLength ℓ_ij · · (b_i − 1) b_j 0 = ℓ_ij·b_j` (Alg. 6 from the corner `i`); for `0 ≤ t ≤ 1` and `P = (1−t)·tail + t·head`: `edgeSplitLengths \|tail head\| (1−t) t = (\|tail P\|, \|head P\|) = (t·ℓ, (1−t)·ℓ)` (intrinsic SPEC N17). The paper's Alg. 8 as printed: `baryLength 1 · · (3/4 − 1) (1/4) 0 = 1/4`, while its line 4 gives `1·(3/4)`. | `cornerDistances_eq`, `cornerDistances_eq_baryLength`, `dist2_bary_eq_baryLength`, `dist2_eq_dist_toC`, `edgeSplitLengths_eq`, `edgeSplitLengths_sum`, `edgeSplitLengths_eq_cornerDistances`, `edgeSplitLengths_param`, `dist2_lerp_tail`, `dist2_lerp_head`, `alg6_distance_on_edge`, `alg6_distance_on_edge'`, `alg8_as_printed_differs` |
| `it-circumcenter` | `circumWeight a_k a_{k+1} a_{k+2} = a_k²·(a_{k+1}² + a_{k+2}² − a_k²)`, `a_k` the side opposite corner `k`; `circumcenterBary a₀ a₁ a₂ = (w₀, w₁, w₂)/S`, `S = w₀ + w₁ + w₂` (this package's `circumSum`). For `T(a₀,a₁,a₂)`: `S = 16·triArea² > 0`; the coordinates sum to `1`; `cornerDistances a₂ a₀ a₁ (circumcenterBary a₀ a₁ a₂) = (R, R, R)` with `R = a₀a₁a₂/(4·triArea) = a₀/(2 sin θ₀)` (equidistant from the three corners); the coordinate of corner `k` is `< 0 ⟺` the angle there is `> π/2`, and `= 0 ⟺` it is `π/2` (angles `cornerAngle a₀ a₁ a₂`, `cornerAngle a₁ a₂ a₀`, `cornerAngle a₂ a₀ a₁`); if the coordinate of corner `0` is negative the other two are positive; `θ₀ < π/6 ⟹ a₀ < R`. | `circumWeight_sum`, `circumSum_pos`, `circumcenterBary_sum`, `circumcenter_equidistant`, `circumradius_eq_div_sin`, `circumcenterBary_neg_iff`, `circumcenterBary_eq_zero_iff`, `circumcenterBary_pos_of_neg`, `circumradius_gt_of_angle_lt` |
| `it-move-length` | `movedLength ℓ r c s = √((ℓ − r c)² + (r s)²)`; `movedLength ℓ r (cos α) (sin α) = √(r² + ℓ² − 2rℓ cos α)` (the paper's Alg. 10); `movedLength ℓ r (cos(ψ − φ)) (sin(ψ − φ)) = \|(r cos φ, r sin φ)(ℓ cos ψ, ℓ sin ψ)\|`; `movedLength ℓ 0 c s = ℓ` for `ℓ ≥ 0`. (The scalar interface has no `cos`/`sin`; the oracle passes `cos α`, `sin α`.) | `movedLength_eq`, `movedLength_eq_law_of_cosines`, `movedLength_eq_dist2`, `movedLength_zero` |
| `it-move-area` | `movedDoubleArea A₂ v d = A₂ − (v₁d₂ − v₂d₁)`; for all plane points `[v,j,k] = movedDoubleArea [0,j,k] v (k − j)`; for `p : ℕ → ℝ²` with `p n = p 0`: `Σ_{k<n} [v, p_k, p_{k+1}] = Σ_{k<n} [w, p_k, p_{k+1}]` for all `v, w` (the star keeps its area). | `movedDoubleArea_eq`, `movedDoubleArea_eq_sda`, `sda_change_centre`, `sum_sda_closed_link` |
| `it-circumcenter-vector` | `weightedCircumcenter ℓ_ab ℓ_ac ℓ_bc A = (A ℓ_ab/2, ℓ_ab(ℓ_ac² + ℓ_bc² − ℓ_ab²)/8)`. With `A ≠ 0`, `ℓ_ab ≠ 0`, `c_x = (ℓ_ab² + ℓ_ac² − ℓ_bc²)/(2ℓ_ab)`, `h = 2A/ℓ_ab`, `c_x² + h² = ℓ_ac²`, and `o = weightedCircumcenter … / A` (componentwise): `\|0 o\|² = \|(ℓ_ab,0) o\|² = \|(c_x,h) o\|²`. Heron form: `T(ℓ_ab, ℓ_ac, ℓ_bc)`, `A = triArea ℓ_ab ℓ_ac ℓ_bc`, `c = (diamondX ℓ_ab ℓ_ac ℓ_bc, diamondH ℓ_ab ℓ_ac ℓ_bc)`: the same. | `weightedCircumcenter_eq`, `weightedCircumcenter_equidistant`, `weightedCircumcenter_equidistant_heron` |
| `it-odt-update` | `weightedBarycenter ℓ_ab ℓ_ac ℓ_bc A = A·⅓((0,0) + (ℓ_ab,0) + (diamondX ℓ_ab ℓ_ac ℓ_bc, 2A/ℓ_ab))`; `odtDisplacement w total α = (α/total)·w`. A triangle `(x,b,c)`: `circumVec x b c = ¼[\|x\|² rot(c−b) − \|b\|² rot(c−x) + \|c\|² rot(b−x)]`, `rot u = (−u₂, u₁)`; for `[x,b,c] ≠ 0` the point `circumVec/(½[x,b,c])` is equidistant from `x, b, c`; `circumVec (0,0) (ℓ_ab,0) (c_x,h) = weightedCircumcenter ℓ_ab ℓ_ac ℓ_bc (ℓ_ab h/2)` (`c_x` as above, `c_x² + h² = ℓ_ac²`). A closed link `p n = p 0`, `\|T_k\| = [x,p_k,p_{k+1}]/2`, `S(x) = Σ_k (\|x\|² + \|p_k\|² + \|p_{k+1}\|²)\|T_k\|/3` (`odtEnergy`), `G = ⅔ Σ_k (\|T_k\| x − circumVec x p_k p_{k+1})` (`odtGradient`), `\|ω\| = Σ_k \|T_k\|` (`starArea`): `HasFDerivAt S (v ↦ ⟨G, v⟩) x` and `t ↦ S(x + t e)` has the derivative `⟨G, e⟩` at `0` (Chen–Holst Eq. 4.5, `d = 2`); `\|ω\|` does not depend on `x`; `S(x + v) = S(x) + ⟨G,v⟩ + \|ω\|·\|v\|²/3`; with `w = Σ_k (circumVec − \|T_k\| x)`: `odtDisplacement w \|ω\| α = −(3α/(2\|ω\|))·G` and, for `\|ω\| ≠ 0`, `S(x + odtDisplacement w \|ω\| α) = S(x) − 3\|G\|²/(4\|ω\|)·α(2 − α)`; for `\|ω\| > 0`, `0 ≤ α ≤ 2` the step does not increase `S`, and `α = 1` gives the minimum of `S` over all `y`. | `weightedBarycenter_eq`, `odtDisplacement_eq`, `circumVec_fst`, `circumVec_snd`, `circumVec_equidistant`, `circumVec_layout`, `hasFDerivAt_triTerm`, `triGrad_split`, `sum_triGrad`, `hasFDerivAt_odtEnergy`, `hasDerivAt_odtEnergy_line`, `starArea_indep`, `odtEnergy_add`, `odtDisplacement_eq_gradient`, `odtEnergy_step`, `odtEnergy_step_le` |

Cited for intrinsic triangulations (`CITED.md`): `it-delaunay-existence-uniqueness`,
`it-delaunay-flip-termination`, `it-rippa-theorem`, `it-chew-second-algorithm`,
`it-ruppert-encroachment`, `it-straightest-geodesic`, `it-odt-energy`,
`it-odt-target-length-splits`.

### Developability of Triangle Meshes (`crates/ddg-developability`)

Stein, Grinspun, Crane, ACM TOG 2018 ("the developability paper"); contract
`crates/ddg-developability/SPEC.md` (decisions `D1`–`D29`). Clean room: written from the paper's
text and the SPEC only. Notation of the oracle (`lean/oracle/DdgOracle/Formulas/Developability.lean`,
namespace `Ddg.Dev`, at `ℝ`); vectors are `Ddg.Point ℝ` (`V3`), a face is `(f_i, f_j, f_k)`,
`n = (f_j − f_i) × (f_k − f_i)`, `N = n/|n|`, `A = |n|/2`, `θ_i = atan2(|a × b|, a·b)`,
`e_i = f_k − f_j`, `e_j = f_i − f_k`, `e_k = f_j − f_i`. **Derivatives** are along arbitrary
differentiable curves (`HasDerivAt3`: each coordinate has a derivative; `DdgProofs/Developability/Vec.lean`),
so "`dF = Σ_p g_p·δ_p`" contains the gradient with respect to each argument and the chain rule.
All theorems are in `namespace DdgProofs.Dev`.

| id | exact statement | theorems |
|---|---|---|
| `dev-face-frame` | `n ≠ 0`: `n = f_i×f_j + f_j×f_k + f_k×f_i`, `\|N\| = 1`, `N ⊥` the three edges, `A = triArea` (Heron) of the three edge lengths. | `face_frame` |
| `dev-corner-angle` | `a, b ≠ 0`: `θ ∈ [0, π]`, `cos θ = a·b/(\|a\|\|b\|)`, `sin θ = \|a×b\|/(\|a\|\|b\|)`, `θ = arccos(a·b/(\|a\|\|b\|))`; for `n ≠ 0` it is the law-of-cosines `cornerAngle` of the lengths and the three angles sum to `π`. | `cornerAngleAt_spec`, `cornerAngleAt_eq_arccos`, `cornerAngleAt_eq_cornerAngle`, `cornerAngleAt_add` |
| `dev-area-gradient` | `n ≠ 0`: `dA = Σ_p (½ N × e_p)·δ_p` (Eq. 7). | `hasDerivAt_faceArea` |
| `dev-normal-gradient` | `n ≠ 0`: `dN = Σ_p (N·δ_p)(e_p × N)/(2A)` (Eq. 8 with `1/(2A)`, D2). | `hasDerivAt3_faceUnitNormal` |
| `dev-angle-gradient` | `n ≠ 0`: `dθ_i = ∇_iθ_i·δ_i + ∇_jθ_i·δ_j + ∇_kθ_i·δ_k`, `∇_jθ_i = N×(f_i−f_j)/\|f_i−f_j\|²`, `∇_kθ_i = N×(f_k−f_i)/\|f_k−f_i\|²`, `∇_iθ_i = −(∇_j + ∇_k)` (Eq. 9 with the squares, D3). | `hasDerivAt_cornerAngleAt` |
| `dev-backprop-face` | `n ≠ 0`: `F = y·N + ζ·n + s_iθ_i + s_jθ_j + s_kθ_k` has `dF = Σ_p g_p·δ_p`, `(g_i, g_j, g_k) = backpropFace`, i.e. `g_p = w × e_p + angle terms`, `w = (y − (y·N)N)/\|n\| + ζ`. | `hasDerivAt_backpropFace` |
| `dev-vertex-normal` | Interior star with cyclic ring `w_m`: `S_v = ½ Σ_m w_m × w_{m+1}` (independent of `f_v`); `N_v = S_v/\|S_v\|` is a unit vector. | `ringVectorArea_eq`, `vertexNormal_unit` |
| `dev-gradient-invariance` | `E : (Fin n → Fin 3 → ℝ) → ℝ`, `HasFDerivAt E E' f`, `E'δ = Σ g·δ`, `E` constant along the translations, the rotations about the three axes and the scalings of `f`: `Σ g = 0`, `Σ f × g = 0`, `Σ f·g = 0`. | `gradient_invariance` |
| `dev-star-bipartitions` | `interiorBipartitions k` = the pairs `a < b < k`, `k(k−1)/2` of them; `F₁ = arcPart a b`, `F₂ = arcComplement k a b` are nonempty, disjoint, cover `0 … k−1`; distinct pairs give distinct unordered partitions; every cyclic interval of length `1 ≤ l < k` is `F₁` or `F₂` of a pair; `boundarySplits k` has `k − 1` entries. ("Edge-connected" = cyclic interval is taken as the definition.) | `length_interiorBipartitions`, `mem_interiorBipartitions`, `arcPart_arcComplement`, `arc_bipartition_injective`, `cycInterval_mem_bipartitions`, `length_boundarySplits` |
| `dev-combinatorial-width` | `n ≥ 1`, unit normals: `(1/n)Σ\|N−N̄\|² = (1/n²)Σ_{σ1<σ2}\|N_σ1−N_σ2\|² = 1 − \|N̄\|²`; over ordered pairs the sum is twice (D1). | `partWidth_eq_pairs_eq_mean`, `sum_ordered_pairs` |
| `dev-combinatorial-gradient` | `n ≥ 1`, any curves `N_σ(t)`: `d partWidth = Σ_σ (2/n)(N_σ − N̄)·N_σ'`. | `hasDerivAt_partWidth` |
| `dev-combinatorial-max` | `\|N_i − N_j\|² ≤ partMaxWidth` for all `i, j ≤ n`, with equality for some pair; `d\|N₁−N₂\|² = y·N₁' − y·N₂'`, `y = 2(N₁ − N₂)`. | `partMaxWidth_spec`, `hasDerivAt_pairTerm` |
| `dev-covariance-matrix` | `uᵀ(covMatrix)u = Σθ⟨u,N⟩²`; for an orthonormal eigenbasis with `λ₀ ≤ λ₁, λ₂`, `λ₀` is the least value of Eq. 3 over unit `u` (Rayleigh). | `quad_covMatrix`, `covObjective_isLeast` |
| `dev-eigenvalue-supergradient` | `\|x\| = 1`, `xᵀAx = λ_A`, `λ_B ≤ uᵀBu` for all unit `u`: `λ_B ≤ λ_A + xᵀ(B−A)x`. | `eigenvalue_supergradient`, `eigenvalue_supergradient_of_eigenbasis` |
| `dev-covariance-gradient` | Curves `θ_m(t)`, `N_m(t)`, fixed `x`: `d(xᵀA x) = Σ ((x·N)² θ' + 2θ(x·N) x·N')` (Eq. 10 up to the cited `dλ = xᵀ dA x`). | `hasDerivAt_quad_covMatrix` |
| `dev-angle-weight-retessellation` | Splitting corner `j` into two corners with the same normal and angles summing to `θ_j` keeps `covMatrix`. | `covMatrix_split_corner` |
| `dev-spherical-circumcenter` | `n = (b−a)×(c−a) ≠ 0`, `v = n/\|n\|`: `v·a = v·b = v·c = det(a,b,c)/\|n\|`, `\|v\| = 1`, `\|n\| = 2A_abc`; a unit `u` with `u·a = u·b = u·c` is `±v`. | `circumcenter_equidistant`, `eq_circumcenter_of_equidistant` |
| `dev-max-covariance-gradient` | Curves `a, b, c` (no unit length needed), `n ≠ 0`: `d((v·a)²) = Σ_q ∂_q·q'`, `∂_q = 2g([q=a]v + ((e_q×v)·a/\|n\|) v)`; for sites `s_q N_q` the normal gets `s_q ∂_q` (D5). | `hasDerivAt_circumG_sq`, `hasDerivAt_circumG_sq_signed` |
| `dev-log-map` | Unit `n0, N`, `s = \|n0×N\| > 0`: `cos φ = n0·N`, `sin φ = s`, `φ > 0`, `Ñ ⊥ n0`, `\|Ñ\| = φ`, `exp_{n0}(Ñ) = N`; `s = 0 < n0·N ⟹ N = n0, Ñ = 0`. | `logMap_spec`, `logMap_of_parallel` |
| `dev-log-map-gradient` | At `t`: `n0, N` unit, `n0·n0' = N·N' = 0`, `s > 0`: `d(z·Ñ) = y_N·N' + y_{n0}·n0'` with `logMapAdjoint` (exact `G`, `H`; corrected Eq. 11, D6). | `hasDerivAt_logMap` |
| `dev-tangent-eigen-derivative` | `t1, t2` orthonormal to first order, `Ã(t)` symmetric, `Ãx = λx`, `x = y₁t₁ + y₂t₂`: `d quad(Ã, y₁t₁ + y₂t₂) = quad(Ã', x)` (basis terms vanish); a 2x2 eigenvector of `TᵀÃT` lifts to a 3D eigenvector when all `Ñ ⊥ t₁×t₂`. | `hasDerivAt_quad_tangent`, `apply_eq_smul_of_tangent_eigen` |
| `dev-intrinsic-width` | `pairWidthSq n P u = (max⟨P,u⟩ − min⟨P,u⟩)²`. | `pairWidthSq_eq_sq_widthAlong` |
| `dev-width-altitude-gradient` | `n0` fixed unit, `w = x_b − x_a ≠ 0`, `w, w' ⊥ n0` at `t`: `du = −(u·w')ŵ/\|w\|`, `u = n0 × ŵ`. | `hasDerivAt3_altitudeDir` |
| `dev-min-width-gradient` | `W = \|d\|² − (d·w)²/\|w\|²` equals `(d·u)²` for `d, w ⊥` unit `n0`; along arbitrary curves of the four points (`d, w ⊥ n0` at `t`) `dW = Σ z·x'` with `minWidthSensitivities`. | `gramWidth_eq_sq`, `hasDerivAt_gramWidth` |
| `dev-jacobi-rotation` | `a_pq ≠ 0`: `t² + 2θt − 1 = 0`, `\|t\| ≤ 1`, `jacobiStep = JᵀMJ` (`J` the rotation by `c, s = tc` in the plane `(p,q)`), `a_pq ← 0`, trace and Frobenius norm kept; 2x2: `(c, −s)`, `(s, c)` orthonormal eigenvectors for `a_pp − t a_pq`, `a_qq + t a_pq`. | `jacobiStep_eq_rotated`, `jacobiT_quadratic`, `abs_jacobiT_le_one`, `jacobiStep_trace`, `jacobiStep_frobenius`, `jacobi2_eigen` |
| `dev-weak-wolfe-curvature` | `t > 0`, `gᵀd < 0`, `c2 < 1`, `C(t)`: `sᵀy > t(c2 − 1)gᵀd > 0`. | `weakWolfe_curvature` |
| `dev-lewis-overton-line-search` | The outcome of `lineSearch`: `accepted t` satisfies `S(t)` and `C(t)`; `capped t` satisfies `S(t)`; `weak t` satisfies `S(t)`, `t > 0`. | `lineSearch_spec`, `lineSearch_accepted` |
| `dev-midpoint-subdivision` | The four children have `n/4` (same `N`); old corners keep their angles; old stars keep their vertex normal (`S → S/4`); a star whose normals (and face normals) are `⊥ w` has Eq. 3, Eq. 6, the intrinsic objectives `= 0` at `w`, `n0 ⊥ w`, `Ñ ⊥ w`; a part with one normal has width `0`. | `child_normals`, `child_unit_normals`, `child_angles`, `vectorArea_children`, `coplanar_energies_vanish`, `partWidth_const` |
| `dev-small-angle-classification` | `n ≠ 0`, `θ_s ≤ π/3`: at most two angles `< θ_s`; if `θ_i, θ_j < θ_s` then `θ_k > π − 2θ_s` and the edge opposite `k` is strictly the longest. | `smallAngle_classification` |
| `dev-regular-cone-values` | Regular cone, `k ≥ 3`: `Σθ⟨u,N_m⟩² = θk(h²/(2ρ²)(u_x²+u_y²) + cos²(π/k)/ρ² u_z²)`, least value over unit `u` `θk min(sin²α/2, cos²α)`; for `h > 0` the vertex normal is `e_z` and the intrinsic objective is `θkα²/2`, `α = atan2(h, cos π/k)`, in every tangent direction. (`E_P`, `λ_max`, intrinsic max: tests only.) | `cone_covObjective`, `cone_covariance_energy`, `cone_intrinsic_energy` |

Cited for the developability paper (`CITED.md`): `dev-hinge-coplanar`, `dev-hinge-antiparallel`,
`dev-hinge-ruled`, `dev-valence-three-flat`, `dev-simple-eigenvalue-derivative`,
`dev-max-covariance-voronoi`, `dev-min-width-hull-edge`, `dev-lbfgs-two-loop`,
`dev-lewis-overton-termination`, `dev-armijo-backtracking`, `dev-flip-manifold`,
`dev-link-condition`, `dev-cut-edge-weight`. Not covered: the Taylor series of the log-map
coefficients (`logGSeries`, `logHSeries`) and the overflow guard of the Jacobi rotation, which the
oracle file marks as numerical devices.

## Theorems

Sections A and B below index the theorems about the canonical problem and the operators. The
theorems about the cone solver (`DdgProofs/Solver.lean`), about Boundary First Flattening
(`DdgProofs/Bff/*.lean`) and about the cut flow (`DdgProofs/Cuts/*.lean`) are indexed by the
three id tables just above. The first two groups are about this package's own definitions, which
transcribe the Rust expressions where the statement says "as the code evaluates it"; the cut-flow
theorems are about the oracle's definitions wherever one exists (see below).

"Def." says which definitions the statement is about: **own** = math-level definitions of this
package, **oracle** = the definitions of `lean/oracle` at `ℝ`. Hypotheses abbreviated: *valid mesh*
= every face has three distinct vertices and edge lengths satisfying the strict triangle
inequalities (`TriMesh.Valid`); *well-formed* = `wE, m, λ, wR ≥ 0`, `lo ≤ 0 ≤ hi`.

### A.1 Scalar convex facts (paper Eqs. 18–19, 21, Props. B.2–B.3)

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `soft-threshold-fixed-point` | `eq_softThreshold_add_iff_isSubgradient` | For `t ≥ 0`: `µ = P_t(φ + µ)` iff `φ ∈ t·∂\|µ\|`. | Eq. 19, App. B.2 | own |
| `soft-threshold-fixed-point` | `eq_softThreshold_add_iff_cases` | … iff `\|φ\| ≤ t`, `φ = t` where `µ > 0`, `φ = −t` where `µ < 0`. | Eq. 19, Prop. 4.1/C.1 (discrete) | own |
| `soft-threshold-fixed-point` | `eq_softThreshold_add_iff` | … iff `\|φ\| ≤ t` and `φ·µ = t·\|µ\|`. | Eq. 19 | own |
| `soft-threshold-fixed-point` | `eq_softThreshold_add_iff_complementarity` | … iff `\|φ\| ≤ t`, `µ(\|φ\| − t) = 0` and `φ·µ ≥ 0`. | Eq. 19 | own |
| — | `eq_softThreshold_zero_add_iff` | Free site (`t = 0`): the inclusion is `φ = 0`. | Sec. 8.1 | own |
| `abs-subdifferential-bounded` | `isSubgradientWithinAt_abs_iff` | `g ∈ ∂(t\|·\| + ι_[lo,hi])(µ)` iff `µ ∈ [lo,hi]`, `g ≤` right slope of `t\|·\|` unless `µ = hi`, `g ≥` left slope unless `µ = lo`. | Sec. 8.2 (derived) | own |
| `cone-fixed-point-bounded` | `isSubgradientWithinAt_abs_iff_eq_clampE` | With bounds: `φ ∈ ∂(t\|·\| + ι_[lo,hi])(µ)` iff `µ = clamp(P_t(φ + µ))`. | Sec. 8.2 (derived) | own |
| `cone-fixed-point-bounded` | `Site.stationary_iff_eq_prox` | The same for every vertex kind (weighted, bounded, free, forbidden, prescribed). | Eq. 19, Sec. 8 | own |
| `cone-fixed-point-bounded` | `Oracle.inSubdiff_iff_eq_boundedThreshold` | `InSubdiff t lo hi µ φ` iff `µ = boundedThreshold t lo hi (φ + µ)`. | Eq. 19 | oracle |
| `regularized-cone-update` | `isSubgradientWithinAt_abs_sub_iff_eq_clampE` | `ψ − γµ ∈ ∂(t\|·\| + ι_[lo,hi])(µ)` iff `µ = clamp(P_t(ψ)/γ)` (`γ > 0`). | Eq. 21, System B.4 | own |
| `regularized-cone-update` | `eq_softThreshold_div_iff_isSubgradient` | Unbounded case: `µ = P_t(φ)/γ` iff `φ − γµ ∈ t∂\|µ\|`. | Eq. 21 | own |
| `regularized-cone-update` | `minimizes_reg_iff_eq_clampE` | `clamp(P_t(ψ)/γ)` is the unique minimizer of `y ↦ t\|y\| + (γ/2)y² − ψy` on `[lo,hi]`. | Prop. B.3 | own |
| `regularized-cone-update` | `Site.stationary_sub_iff_eq_proxReg` | The regularized inclusion in fixed-point form for every vertex kind. | Eq. 21 | own |
| `regularized-cone-update` | `Oracle.inSubdiff_sub_div_iff_eq_regularizedMeasure` | `InSubdiff t lo hi µ (φ − µ/k)` iff `µ = regularizedMeasure k t lo hi φ` (`k = m/γ`). | Eq. 21 with mass factor | oracle |
| `regularized-cone-update` | `Oracle.eq_regularizedMeasure_add_iff` | `µ = clamp(c·P_t(φ + µ))`, `c = m/(m+γ)`, iff `φ − (γ/m)µ ∈ ∂g(µ)`. | not in the paper (oracle's form) | oracle |

### A.2 The optimality system (paper Eq. 17, System B.2)

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `optimality-system-sufficient` | `ConeProblem.KKT.isMinimizer` | If `(u, φ, µ)` solves the optimality system (and `wE, m ≥ 0`) then `(u, µ)` is a global minimizer. Holds for every matrix `L`. | Sec. 4.4, App. A.3 (discrete) | own |
| `optimality-system-sufficient` | `ConeProblem.KKT.objective_add_le` | Quadratic growth: every feasible `(v, ν)` has `J(v,ν) ≥ J(u,µ) + ½ Σ wE m (v − u)²`. | — | own |
| `optimality-system-sufficient` | `ConeProblem.KKT.minimizes_reduced_of_posDef` | Dirichlet case (`L` positive definite): `µ` minimizes `ν ↦ J(L⁻¹(ν − Ω), ν)`. | Sec. 5.2.1 | own |
| `optimality-system-sufficient` | `ConeProblem.KKT.minimizes_closed` | Closed case (`L·1 = 0`): `(u, µ)` beats `(v + c, ν)` for every solution `v` and every constant `c`. | (closed meshes: not in the paper) | own |
| `optimality-system-sufficient` | `Ddg.Problem.objective_le_of_F_eq_zero` | A zero of the oracle's residual map `F` is feasible and minimizes the oracle's `objective`. | Eq. 20 | oracle |
| `minimizer-unique` | `ConeProblem.IsMinimizer.unique` | If all `wE_i m_i > 0`, two minimizers have the same `u` and the same `µ` (any `L`). | — | own |
| — | `ConeProblem.IsMinimizer.weighted_sq_sub_eq_zero` | In general two minimizers agree in `u` wherever `wE_i m_i > 0`. | — | own |
| `minimizer-unique` | `Ddg.Problem.unique_of_minimizers` | The same for the oracle's problem. | — | oracle |
| `optimality-system-necessary` | `ConeProblem.IsMinimizer.existsUnique_kkt` | If `det L ≠ 0`, every minimizer has exactly one adjoint state `φ` completing it to a solution of the optimality system. | Sec. 4.4 (discrete) | own |
| `optimality-system-necessary` | `ConeProblem.isMinimizer_iff_exists_kkt_of_posDef` | Dirichlet case: `(u, µ)` is a minimizer iff some `φ` solves the optimality system with it. | — | own |

Necessity on closed meshes is **not** proved (see `CITED.md`, `optimality-system-necessary-closed`).

### A.3 Guarantee of paper Sec. 4.5

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `least-distortion-guarantee` | `ConeProblem.IsMinimizer.distortion_le_add` | `E(u) ≤ E(v) + λ(R(ν) − R(µ))` for every feasible `(v, ν)`. | Sec. 4.5 | own |
| `least-distortion-guarantee` | `ConeProblem.IsMinimizer.distortion_le` | A minimizer has the least distortion among all feasible `ν` with `Σ wR\|ν\| ≤ Σ wR\|µ\|`. | Sec. 4.5 | own |

### A.4 Dual problem and certificate (paper Problem 10 / A.1, discrete)

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `duality-gap-certificate` | `ConeProblem.gap_identity` | `J(u,µ) − D(ũ,φ) = ½ Σ wE m (u − ũ)² + Σ_i (λ cost_i(µ_i) + cost_i*(φ_i) − φ_i µ_i)`. | — | own |
| `duality-gap-certificate` | `ConeProblem.weak_duality` | Every dual-feasible value is ≤ every primal-feasible value. | Sec. 3.4, 4.3 (discrete) | own |
| `duality-gap-certificate` | `ConeProblem.duality_gap_certificate` | Given a primal-feasible and a dual-feasible point, if `J − D ≤ gap` then `J` is within `gap` of the objective of every feasible point. | — | own |
| `duality-gap-certificate` | `ConeProblem.adjoint_gap_certificate` | The same with the dual point built from `φ` alone (`ũ = −Lᵀφ/(wE m)`), all `wE m > 0`. | — | own |
| `duality-gap-certificate` | `ConeProblem.KKT.dualObjective_eq` | At a solution of the optimality system the gap is zero (strong duality). | Sec. 3.4 | own |
| — | `Site.fenchel_young`, `Site.exists_conj_eq`, `Site.Stationary.fenchel_young_eq` | Fenchel–Young for the per-vertex cost, attainment, and equality at a stationary pair. | Prop. A.3 | own |
| `duality-gap-certificate` | `Ddg.Problem.objective_le_add_gap` | The certificate for the oracle's problem (finite bounds: every `φ` is dual feasible). | — | oracle (dual value: own formula `dualValue`) |

### A.4′ The certificate for inexact data (`DdgProofs/Certificate.lean`)

No assumption on `(u, φ, µ)` unless stated. `r = L u − µ + Ω` (`stateResidual`),
`φ̂ = dualClip φ`, `D(φ) = dualObjectiveOfAdjoint φ`, `G = complementarityGap`,
`cost_i(x) = wR_i|x|` (`0` at a fixed vertex). "Admissible" = in `[lo_i, hi_i]`, resp. equal to the
fixed angle.

| id | theorem | statement | Def. |
|---|---|---|---|
| `duality-gap-inexact` | `ConeProblem.gap_identity_general` | `J(u,µ) − D(ũ,φ) = ½Σw(u−ũ)² + Σ(λcost_i(µ_i) + conj_i(φ_i) − φ_iµ_i) − Σφ_i r_i + Σ s_i u_i` for arbitrary `u, ũ, φ, µ`. | own |
| `duality-gap-inexact` | `ConeProblem.objective_sub_dualObjectiveOfAdjoint` | `J(u,µ) − D(φ) = G(u,φ,µ) − Σφ_i r_i + Σ_{w_i=0}(Lᵀφ)_i u_i` (the relation between the code's `duality_gap` and `complementarity_gap`). | own |
| `duality-gap-inexact` | `ConeProblem.abs_gap_sub_complementarityGap_le` | `\|J(u,µ) − D(φ) − G(u,φ,µ)\| ≤ Σ\|φ_i r_i\| + Σ_{w_i=0}\|(Lᵀφ)_i u_i\|` (the code's `gap_uncertainty`; it bounds the difference of the two forms of the gap, not the suboptimality). | own |
| `duality-gap-certificate` | `ConeProblem.dualObjectiveOfAdjoint_add_le` | `D(φ) + Σ_{w_i=0}(Lᵀφ)_i v_i ≤ J(v,ν)` for dual-feasible `φ` and feasible `(v,ν)`. | own |
| `duality-gap-certificate` | `ConeProblem.dual_lower_bound`, `dual_lower_bound_of_clip` | `D(φ̂) ≤ J(v,ν)` for every feasible `(v,ν)`, for **any** `φ`, if `(Lᵀφ̂)_i = 0` wherever `w_i = 0`. | own |
| — | `ConeProblem.feasible_iff_eq_inducedCone` | `(u,ν)` is feasible ⟺ `ν = L u + Ω` and `L u + Ω` is admissible: `µ + r` is the only repair that keeps `u`. | own |
| `duality-gap-certificate` | `ConeProblem.objective_add_stateResidual_le`, `repairCost_le_of_abs_le` | `J(u, µ+r) ≤ J(u,µ) + λΣcost_i(r_i) ≤ J(u,µ) + ε·λΣ_i wR_i` if `\|r_i\| ≤ ε`. | own |
| `duality-gap-certificate` | `ConeProblem.certificate_enclosure` | If `φ` is dual feasible, `(Lᵀφ)_i = 0` wherever `w_i = 0`, and `µ + r` is admissible: `(u, µ+r)` is feasible; `D(φ) ≤ J(v,ν)` for all feasible `(v,ν)`; `J(u,µ+r) ≤ J(u,µ) + λΣcost_i(r_i) =: U`; `J(u,µ+r) ≤ J(v,ν) + (U − D(φ))`; `0 ≤ U − D(φ)`. | own |
| `duality-gap-certificate` | `ConeProblem.certificate_of_computed` | The same from a computed triple, all `w_i > 0`: with `J(u,µ) − D(φ̂) ≤ gap` and `λΣcost_i(r_i) ≤ ε`, `(u, µ+r)` is within `gap + ε` of the optimum. | own |
| `duality-gap-certificate` | `ConeProblem.certificate_of_computed_unbounded` | If no vertex has a bound or a fixed angle (the default problem), `µ + r` is always admissible: the previous statement holds for **every** triple. | own |
| `duality-gap-certificate` | `ConeProblem.IsMinimizer.objective_mem_enclosure`, `IsMinimizer.weighted_sq_sub_le_width` | Every minimizer `(u*,µ*)` has `D(φ) ≤ J(u*,µ*) ≤ U` and `¼Σw(u−u*)² ≤ U − D(φ)`. | own |
| `duality-gap-certificate` | `ConeProblem.repaired_certificate` | For admissible `µ`: with `µ̄` the admissible vector nearest to `µ+r` and `ρ = µ+r−µ̄`, `(u,µ̄)` is feasible for the problem with defects `Ω − ρ` and within `J(u,µ) + λΣcost_i(r_i) − D(φ) + Σφ_iρ_i` of its optimum. | own |
| `duality-gap-certificate` | `ConeProblem.backward_certificate` | For admissible `µ`: `(u,µ)` is feasible for the problem with defects `Ω − r` and within `G(u,φ,µ) ≥ 0` of its optimum. | own |
| `duality-gap-certificate` | `ConeProblem.KKT.objective_add_le_of_residual` | If `(u*,φ*,µ*)` solves the optimality system: `J(u*,µ*) + ½Σw(u−u*)² ≤ J(u,µ) + Σφ*_i r_i` for admissible `µ`. | own |
| `minimizer-unique` | `ConeProblem.IsMinimizer.weighted_sq_sub_le` | `¼Σw(u−v)² ≤ J(u,µ) − J(v,ν)` for a minimizer `(v,ν)` and feasible `(u,µ)`. | own |

Five `example`s at the end of the file: the hypotheses of `certificate_of_computed` hold with
nonzero residuals in a case where `J(u,µ) − D < 0`; `J(u,µ)` can lie below the optimum; the
enclosure is false at a forbidden vertex with `r ≠ 0` and at a vertex on its bound; the lower
bound is false when `w_i = 0` and `(Lᵀφ)_i ≠ 0`.

### A.5 Semismooth Newton method (paper Sec. 5.3)

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `soft-threshold-active-sets` | `softThreshold_isPiecewiseAffine` | `P_t` is affine on each configuration of the active sets, with slope `A⁺ + A⁻`. | Eq. 24, Sec. 5.3.2 | own |
| `soft-threshold-active-sets` | `Site.proxReg_isPiecewiseAffine` | The same for `clamp(P_t(·)/γ)` with five pieces (two bound-active ones). | Sec. 5.3.2, 8.2 | own |
| `soft-threshold-active-sets` | `Oracle.regularizedMeasure_eq_value`, `Oracle.regularizedMeasure_isPiecewiseAffine` | The oracle's `coneState` pieces describe `clamp(k·P_t(·))` exactly. | — | oracle |
| `newton-stopping-rule` | `paper_stageII_stopping_rule` | Stage II, paper's form: if a Newton step leaves `A⁺(φ+µ)`, `A⁻(φ+µ)` unchanged, the new iterate is an exact zero of `F`. | Eqs. 20, 26, Sec. 5.3.3 | own |
| `newton-stopping-rule` | `paper_stageI_stopping_rule` | Stage I, paper's form (`A^±(φ)`, per-vertex `γ_i`). | Eqs. 22, 25, Sec. 5.3.3 | own |
| `newton-stopping-rule` | `ConeProblem.residual_eq_zero_of_newton_step` | Stage II for the canonical problem with all vertex kinds. | Sec. 5.3.3 | own |
| `newton-stopping-rule` | `ConeProblem.isMinimizer_of_newton_step` | … and then the new iterate is a global minimizer. | Sec. 5.3.3 | own |
| `newton-stopping-rule` | `ConeProblem.residualReg_eq_zero_of_newton_step`, `ConeProblem.isRegMinimizer_of_newton_step` | Stage I for the canonical problem; the new iterate minimizes the regularized objective. | Sec. 5.3.3 | own |
| `newton-stopping-rule` | `Ddg.Problem.F_eq_zero_of_newton_step`, `Ddg.Problem.objective_le_of_newton_step` | Stage II on the oracle's `stateResidual`, `adjointResidual`, `coneResidual`, `state`: a step that keeps the **canonical** state (`ConeState.canon`) of every vertex lands on a zero of `F`. (Raw states flip at a prescribed vertex without any change of the equation; an `example` shows it.) | Sec. 5.3.3 | oracle |
| `newton-stopping-rule` | `Ddg.Problem.Freg_eq_zero_of_newton_step` | Stage I on the oracle's `Freg`, `regState`, canonical states, for `γ > 0`. | Sec. 5.3.3 | oracle |
| `soft-threshold-active-sets` | `Oracle.regularizedMeasure_affine_of_canon_eq` | `clamp(k·P_t(·))` is affine between two points whose pieces have the same canonical name. | — | oracle |
| `newton-stopping-rule` | `Ddg.Problem.zero_off_of_Fpiece_eq_zero` | The oracle's piecewise linear system with absent bounds (see below). | — | oracle |
| `newton-stopping-rule` | `Ddg.Problem.Fimp_eq_zero_of_Fpiece_eq_zero`, `Ddg.Problem.Fimp_eq_zero_of_newton_step` | The oracle's own iteration: a solution of `DFpiece γ σ · y = −Fpiece γ σ` whose new point lies on the pieces `σ` is a zero of `Fimp γ`. | — | oracle |
| — | `Ddg.Problem.Fpiece_add`, `Ddg.Problem.DFpiece_mulVec` | `DFpiece` is the Jacobian of the affine map `Fpiece`. | — | oracle |
| — | `Ddg.Problem.Fimp_zero_eq_F` | `Fimp 0 = F`. | — | oracle |
| `regularized-optimality-sufficient` | `ConeProblem.isRegMinimizer_of_stationary`, `ConeProblem.isRegMinimizer_of_residualReg_eq_zero` | A zero of the Stage I residual minimizes `J(u,µ) + Σ (γ_i/2) µ_i²`. | Eq. 12–13, System B.3 | own |
| `regularized-optimality-sufficient` | `Ddg.Problem.regObjective_le_of_Fimp_eq_zero` | A zero of the oracle's `Fimp γ` minimizes `J + Σ (γ/(2 m_i)) µ_i²`. | Eq. 12 with mass factor | oracle |
| `regularized-suboptimality` | `ConeProblem.IsRegMinimizer.objective_le` | `J(u_γ, µ_γ) ≤ J(u, µ) + Σ (γ_i/2) µ_i²` for every feasible `(u, µ)`. | Sec. 4.6 | own |

Convergence of the iteration is not proved (`CITED.md`, `semismooth-newton-convergence`).

### A.6, A.7 Closed meshes; scaling

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `closed-mesh-cone-sum` | `ConeProblem.sum_cone_eq_sum_defect` | If `1ᵀL = 0` and `L u = µ − Ω` then `Σµ = ΣΩ`. | Sec. 8.1 (Gauss–Bonnet remark) | own |
| `closed-mesh-cone-sum` | `TriMesh.sum_cone_eq_two_pi_eulerChar` | On a closed valid manifold mesh `Σµ = 2πχ`. | — | own |
| `closed-mesh-optimal-constant` | `ConeProblem.distortion_add_const`, `distortion_add_const_eq` | `E(u + c)` is a parabola in `c`; its unique minimizer is `c₀ = −Σ wE m u / Σ wE m`. | — | own |
| `closed-mesh-optimal-constant` | `ConeProblem.KKT.weighted_mean_eq_zero`, `IsMinimizer.weighted_mean_eq_zero` | At the optimum on a closed mesh `Σ wE m u = 0`. | — | own |
| `lambda-area-normalization` | `ConeProblem.isMinimizer_scaleMass_iff` | Multiplying all masses by `c > 0` equals replacing `λ` by `λ/c`. | Sec. 6.2 | own |
| `lambda-area-normalization` | `ConeProblem.normalizeMass_scaleMass` | `normalizeMass P A` divides the masses by `A`, to be taken as the **total area** (the sum of the lumped masses over *all* vertices, boundary included — not the sum over the unknowns). If masses and `A` are multiplied by the same `c ≠ 0`, the normalized problem is unchanged. | Sec. 6.2 | own |
| `lambda-area-normalization` | `ConeProblem.isMinimizer_normalizeMass_iff` | Dividing the masses by `A > 0` is the same as replacing `λ` by `λ·A`. | Sec. 6.2 | own |
| `scale-invariance` | `cornerAngle_smul`, `cornerCot_smul`, `triArea_smul` | Angles and cotangents are scale invariant; areas scale by `s²`. | Sec. 6.2 | own |
| `scale-invariance` | `TriMesh.laplacian_scale`, `TriMesh.mass_scale`, `TriMesh.angleDefect_scale` | `L` and `Ω` are scale invariant; masses scale by `s²`. | Sec. 6.2 | own |

### A.8 Connected components (`DdgProofs/Components.lean`)

`P.sum Q`: block diagonal `L = fromBlocks P.L 0 0 Q.L`, `Ω, m, wE, site` by `Sum.elim`, `λ = P.lam`;
every theorem that involves `λ` assumes `P.lam = Q.lam`. Several components: iterate.

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `component-separation` | `ConeProblem.objective_sum` | `J(u ⊕ v, µ ⊕ ν) = J_P(u, µ) + J_Q(v, ν)`. | — | own |
| `component-separation` | `ConeProblem.feasible_sum_iff` | Feasible for the sum ⟺ each part feasible (no hypothesis on `λ`). | — | own |
| `component-separation` | `ConeProblem.isMinimizer_sum_iff` | Global minimizer of the sum ⟺ each part a global minimizer of its component (no other hypothesis; the additive constant of each closed component is free). | — | own |
| `component-separation` | `ConeProblem.kkt_sum_iff` | Optimality system of the sum ⟺ that of each part. | — | own |
| `component-separation` | `ConeProblem.isMinimizer_sum_iff_alone` | For `a, b > 0`: minimizer of the sum ⟺ each part minimizes its component with masses and `λ` multiplied by `a` (resp. `b`). | — | own |
| `closed-mesh-cone-sum` | `ConeProblem.Feasible.sum_cone_eq_sum_defect_inl`, `…_inr` | A closed component (`1ᵀL_P = 0`) of a feasible pair of the sum has `Σ µ = Σ Ω` on that component. | — | own |
| `component-lambda-rescaling` | `Site.stationary_mul_iff` | `c > 0`: `Stationary (c·λ) s (c·φ) µ ⟺ Stationary λ s φ µ`. | — | own |
| `component-lambda-rescaling` | `ConeProblem.isMinimizer_scaleMass_withLam_iff` | `c > 0`: masses `c·m` with `c·λ` have the minimizers of `P`. | — | own |
| `component-lambda-rescaling` | `ConeProblem.kkt_scaleMass_withLam_iff`, `kkt_scaleMass_withLam_iff'` | `c > 0`: `KKT_c(u, c·φ, µ) ⟺ KKT(u, φ, µ)`; `KKT_c(u, φ, µ) ⟺ KKT(u, c⁻¹·φ, µ)`. | — | own |
| `component-lambda-rescaling` | `ConeProblem.isMinimizer_normalizeMass_alone_iff`, `kkt_normalizeMass_alone_iff` | `A, A_c > 0`: masses `a_i/A_c` with `(A/A_c)·λ` ⟺ masses `a_i/A` with `λ`; adjoints `φ_global = (A_c/A)·φ_alone`. | — | own |

### B.8 Triangles from edge lengths

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `corner-angle` | `cos_cornerAngle` | Law of cosines: `cos(angle opposite a) = (b² + c² − a²)/(2bc)`. | Sec. 5.1.2 | own |
| `corner-angle` | `cornerAngle_eq_two_arctan` | The half-angle form `2·atan √(((a−b+c)(a+b−c))/((a+b+c)(b+c−a)))` is the same angle. | — | own |
| `corner-angle` | `Oracle.cornerAngle_eq` | The oracle's `cornerAngle` is that angle, for every nondegenerate triangle. | — | oracle |
| `corner-angle` | `angle_eq_cornerAngle` | For points of a Euclidean space it is Mathlib's angle `∠ p₁ p₂ p₃`. | — | own |
| `corner-angle` | `two_mul_atan2_sqrt_eq_cornerAngle`, `kahanCornerAngle_eq_cornerAngle` | `2·atan2(√num, √den)`, also in Kahan's arrangement as `corner_angle` in `ddg-core` evaluates it, is that angle. | — | own (mirrors the Rust expression) |
| `triangle-area` | `sin_cornerAngle` | `sin(angle opposite a) = 2·Area/(bc)` with `Area` from Heron's formula. | — | own |
| `triangle-area` | `half_mul_sin_angle_eq_triArea` | Heron: for points of a Euclidean space `½·\|p₁p₂\|·\|p₃p₂\|·sin ∠ = triArea`. | Eq. 16 | own |
| `triangle-area` | `kahanArea_of_sorted_rearrangement`, `sixteen_mul_triArea_sq_eq` | With sides sorted `x ≥ y ≥ z`, the test `z − (x − y) > 0` is the triangle inequality and Kahan's expression is the Heron area, as `triangle_area` in `ddg-core` evaluates it. | — | own (mirrors the Rust expression) |
| `corner-cot` | `cornerCot_eq_cot` | `(b² + c² − a²)/(4·Area)` is the cotangent of the angle opposite `a`. | Sec. 5.1.3 | own (`Oracle.cornerCot_eq`, `Oracle.triArea_eq` identify the oracle's) |
| `corner-cot` | `kahanCornerCot_eq_cornerCot` | The numerator `c² + (b − a)(b + a)` with the sides sorted, as `corner_cotan` evaluates it, gives the same value. | — | own (mirrors the Rust expression) |
| `triangle-angle-sum` | `cornerAngle_add_cornerAngle_add_cornerAngle` | The three angles computed from the lengths sum to `π`. | — | own |
| `edge-length-euclidean` | `dist_eq_sqrt_sum_sq`, `Oracle.dist_eq` | Distance in `ℝ³` is `√(Δx² + Δy² + Δz²)`. | Sec. 5.1.1 | own / oracle |

### B.9, B.10 Cotan Laplacian and mass

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `cotan-formula` | `cotan_formula` | On a triangle, `Area·‖∇f‖²` of an affine `f` equals `½ Σ cot θ_k (f_i − f_j)²` over its three edges. | Sec. 5.1.3 | own |
| `cotan-laplacian` | `TriMesh.laplacian_apply_of_ne`, `TriMesh.laplacian_apply_self` | Entries: `L_ij = −½ Σ_f cot(corner opposite ij)` for `i ≠ j`, `L_ii = −Σ_{j≠i} L_ij`. | Sec. 5.1.3 | own |
| `cotan-laplacian` | `TriMesh.laplacian_isSymm`, `Oracle.laplacian_symm` | `L` is symmetric. | Sec. 5.1.3 | own / oracle |
| `cotan-laplacian` | `TriMesh.laplacian_mulVec_one`, `Oracle.laplacian_row_sum` | Rows of `L` sum to zero. | Sec. 5.1.3 | own / oracle |
| `cotan-laplacian` | `TriMesh.dotProduct_laplacian_mulVec_eq_sum_faceEnergy` | `uᵀLu = Σ_faces ½ Σ cot θ_k (u_i − u_j)²`. | Sec. 5.1.3 | own |
| `cotan-laplacian-psd` | `cotEnergy_nonneg`, `TriMesh.laplacian_posSemidef`, `Oracle.laplacian_quadratic_nonneg` | `L` is positive semidefinite on every valid mesh, Delaunay or not. | Sec. 5.1.3 | own / oracle |
| `cotan-laplacian-kernel` | `TriMesh.laplacian_mulVec_eq_zero_iff` | On a connected valid mesh `L u = 0` iff `u` is constant. | — | own |
| `cotan-laplacian-interior-posdef` | `TriMesh.laplacian_interior_posDef` | On a connected valid mesh with a boundary vertex the interior block `L_II` is positive definite. | Sec. 5.1.3 | own |
| — | `TriMesh.laplacian_mulVec_of_zero_off` | If `u = 0` on the boundary, `(L u)_i = (L_II u_I)_i` at interior `i`. | Sec. 5.1.3 | own |
| `lumped-mass` | `TriMesh.sum_mass`, `TriMesh.sum_mass_div_totalArea`, `Oracle.sum_lumpedMass` | Lumped masses sum to the total area; normalized masses sum to 1. | Eq. 16 | own / oracle |

### B.11 Curvature and Gauss–Bonnet

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `gauss-bonnet`, `angle-defect`, `boundary-geodesic-curvature`, `euler-characteristic` | `TriMesh.gauss_bonnet` | `Σ_interior Ω + Σ_boundary k = 2π(\|V\| − \|E\| + \|F\|)` for a valid mesh in which every edge lies in at most two faces and every boundary vertex lies on exactly two boundary edges. | Eqs. 14–15 | own |
| `gauss-bonnet` | `TriMesh.gauss_bonnet_closed` | Closed mesh: `Σ Ω = 2πχ`. | Sec. 8.1 | own |
| `gauss-bonnet`, `angle-defect`, `boundary-geodesic-curvature` | `TriMesh.gauss_bonnet_fields` | The same for the per-vertex fields `ddg-core` stores (`Ω` set to 0 on boundary vertices, `k` set to 0 on interior vertices). | Eqs. 14–15 | own (mirrors the Rust fields) |
| `gauss-bonnet`, `angle-defect`, `boundary-geodesic-curvature` | `Oracle.gauss_bonnet`, `Oracle.angleDefect_eq`, `Oracle.geodesicCurvature_eq` | The same for the oracle's `angleDefect`, `geodesicCurvature` (boundary taken from the `TriMesh`). | Eqs. 14–15 | oracle |

`angle-defect`, `boundary-geodesic-curvature` and `euler-characteristic` are definitions; the
theorem carrying their ids is the one that ties them together.

### B.12 Linearized discrete Yamabe equation

| id | theorem | statement | paper | Def. |
|---|---|---|---|---|
| `discrete-yamabe-linearization` | `hasDerivAt_cornerAngle_conformal` | Under `ℓ̃_ij = e^{t(v_i+v_j)/2} ℓ_ij` the angle at `i` of a triangle `ijk` has derivative `½cot θ_k (v_j − v_i) + ½cot θ_j (v_k − v_i)` at `t = 0`. | (Springborn et al. 2008; BFF) | own |
| `discrete-yamabe-linearization` | `TriMesh.hasDerivAt_angleSum_conformalScale` | The derivative of the angle sum at a vertex is `−(L v)_i`. | Sec. 5.2.1 | own |
| `discrete-yamabe-linearization` | `TriMesh.hasDerivAt_angleDefect_conformalScale`, `…_add` | Along the ray `t ↦ t·v` the angle defect of the scaled metric has derivative `(L v)_i` at `t = 0`; the same at any base point `u` with the Laplacian of the scaled metric. (A derivative along each ray; no expansion `Ω + L u + o(‖u‖)` uniform in the direction is proved.) | Eq. 17 | own |

This confirms the sign convention of `CLAUDE.md` by proof: with a positive semidefinite `L` and
`u` the true log scale factor, the curvature after scaling by `t·v` is `Ω + t·L v` to first order
in `t`, so prescribing `µ` gives the linearized equation `L u = µ − Ω`, not the paper's
`L u = Ω − µ`.

## Which definitions are shared with the code that runs

The oracle (`lean/oracle`, `DdgOracle/Formulas/*`) defines the formulas once for a generic scalar;
`DdgProofs/Oracle/Scalar.lean` provides the instance `Ddg.Scalar ℝ` (every field the genuine real
operation). The instance is `scoped`: it is active inside `namespace DdgProofs.Oracle` and after
`open DdgProofs.Oracle`, and does not affect how arithmetic on `ℝ` is elaborated elsewhere.
Theorems **about the oracle's definitions** (at `ℝ`, i.e. exact arithmetic):

- scalar: `softThreshold`, `clamp`, `boundedThreshold`, `regularizedMeasure`, `InSubdiff`,
  `coneState`, `ConeState.value/indicator/canon`, `cornerAngle`, `triArea`, `cornerCot`, `dist`;
- mesh: `laplacian`, `angleSum`, `angleDefect`, `geodesicCurvature`, `areaAt`, `totalArea`,
  `lumpedMass`;
- problem: `Problem.applyL`, `stateResidual`, `adjointResidual`, `coneResidual`, `F`, `state`,
  `regMeasure`, `regStateResidual`, `Freg`, `regState`, `coneScale`, `regConeResidual`, `impState`,
  `Fimp`, `pieceConeResidual`, `Fpiece`, `DFpiece`, `objective`;
- cut flow (`DdgOracle/Formulas/Cuts.lean`, every definition of that file): `normalizedAngle`,
  `cornerNormalized`, `faceCurvature`, `faceCurvatureDensity`, `Face.omega`,
  `Face.omegaCancellationFree`, `crossingParam`, `clampCrossing`, `crossingMoved`, `baryLength2`,
  `baryLength`, `baryInner`, `cotFromInner`, `det3`, `subTriangleArea`, `triangleMean`,
  `weightedOffDiag`, `weightedDiag`, `lengthStepEntry`, `penaltyStepEntry`, `normalDerivative`,
  `dualLength`, `normalFlux`, `adjointFlux`, `adjointRhsScaledInterior`,
  `adjointRhsScaledBoundary`, `adjointRhsPointwise`, `gradDirichlet`, `gradScaled`,
  `gradPointwise`, `gradDirichletD`, `gradScaledD`, `gradPointwiseD`, `faceMeanSquareEdge`,
  `fluxCoupling`, `fluxDiagEntry`, `fluxOffEntry`, `fluxRhs`, `vertexCurvatureDensity`,
  `featureSizeRaw`, `featureWeight`, `cutSpeed`, `regionDensity`, `aspectRatio`, `smoothedSpeed`,
  `explicitStep`, `stepCap`, `classMean`, `hopfLaxEdge`, `hopfLaxTriangle`; and of the revision
  "default stages, N regions": `angleFromCot`, `turningAngle`, `normalDerivativeOfLinear`,
  `footFraction`, `heightOver`, `lengthForce`, `segmentMeanWeight`, `rateCoupling`, `rateRhs`,
  `gradientNormLone`, `meanGradient`, `weightedMean`, `ringMean`, `roughQ`, `roughRate`,
  `lengthWeightForRipples`, `adaptiveStep`, `crossingParamRegions`, `triplePointCoord`,
  `tripleRate`, `projectionShift`, `projectedValue`; and of the revision "robust step
  statistics": `insertAscending`, `sortAscending`, `sumLeft`, `firstReaching`,
  `weightedQuantile`, `speedInEdges`, `robustStep`, `localSpeed`, `localLimitFactor`; and of
  the revision "D26": `vertexChange`, `vertexSpeedCapped`, `vertexSpeedAdaptive`,
  `vertexLimitFactor`, `limitedValue`;
- intrinsic triangulations (`DdgOracle/Formulas/Intrinsic.lean`, every definition of that file):
  `angleScale`, `reduceFullTurn`, `signpostUpdate`, `diamondX`, `diamondH`, `flippedLength`,
  `flippedLengthFromLengths`, `flipArea1`, `flipArea2`, `flipConvex`, `delaunayExcess`,
  `delaunayHolds`, `cotanWeight`, `cotanWeightFromLengths`, `oppositeAngles`, `unfoldApex`,
  `axisCrossingT`, `axisCrossingX`, `apexAboveAxis`, `signedArea2`, `barycentric`,
  `pointFromDistances`, `planeDistance`, `cornerDistances`, `edgeSplitLengths`, `circumWeight`,
  `circumcenterBary`, `movedLength`, `movedDoubleArea`, `weightedCircumcenter`,
  `weightedBarycenter`, `odtDisplacement`; they reuse `triArea`, `cornerAngle`, `cornerCot` of the first paper
  (`Oracle.triArea_eq`, `Oracle.cornerAngle_eq`, `Oracle.cornerCot_eq`).

Own definitions of this package for intrinsic triangulations are shorthands and the vehicle of
the proofs, never the subject of a tagged claim on their own: `Intrinsic.apexK`, `apexL` (the
two apexes of the standard layout), `DiamondConvex` (the oracle's `flipConvex` of the standard
layout of five lengths), `rigid` (a rotation and a translation of the plane), `MeetsAxis` (a
closed segment has a point on the x-axis), `layoutPoint` (`b₁v₁ + b₂v₂`), `circumSum` (the sum
of the three weights), `rot`, `nsq`, `toDual`, `circumVec`, `odtEnergy`, `odtGradient`,
`starArea`, `triGrad`, `forceF` (the vehicle of Eq. 4.5: a quarter turn, the squared norm, a
vector as a linear form, `|T| c_T` as a polynomial, the energy and the gradient of one star, one
term of the gradient and the telescoping force of Chen and Holst's proof); plane points, `dist2` and `sda` are those of `DdgProofs/Bff/Quality.lean`,
`IsTriangle`, `triArea`, `cornerAngle`, `cornerCot` those of `DdgProofs/Triangle.lean`.

Own definitions of this package for the default stages and `N` regions (what SPEC states in
words or as a sum over an index set, for which the oracle has only the scalar pieces):
`Cuts.wmean`, `localExtension`, `gradNorm`, `Cuts.meanGradient` (identified with the oracle's
`gradientNormLone`, `meanGradient`, `weightedMean`, `ringMean` by `gradientNormLone_eq`,
`oracle_meanGradient_eq`, `oracle_weightedMean_eq`, `oracle_ringMean_mem`), `rippleRate`,
`rippleFactor`, `sphereRippleRate` (no oracle counterpart), `IsRegion`, `IsTwoSmallest`,
`regionRate` (the argmin with its tie rule, the two smallest values, the substitution
`φ̇_j := −φ̇_k`); for the step statistics `totalWeight`, `weightLE`, `weightLT`,
`IsWeightedQuantile` (the quantile of SPEC §1 stated without an order; the oracle's
`weightedQuantile` satisfies it, `oracle_weightedQuantile_spec`), `quantileScan`,
`Cuts.weightedQuantile` (equal to the oracle's, `oracle_weightedQuantile_eq`) and `localFactor`
(the Rust form of the factor; equal to the oracle's `localLimitFactor` for
`max_motion · h_c ≥ 0`, `localLimitFactor_eq`); for the stop rule `EnergyVerdict`,
`energyVerdict`, `lowestState`, `runOutcome`, `motionCriterion`, `maxAbs`,
`predictedMotionOnCut`, `predictedMotionOnMesh` (no oracle counterpart: the oracle replays
single steps); for the mean slope `pieceSlope` (the harmonic mean over the counted faces, or
`1`: the rule of the Rust `gradient_norms`).

For the cut flow the matrices are not oracle definitions (the executable assembles them from the
entries above and solves densely), so the linear-algebra theorems are about matrices of this
package built from those entries: `TriMesh.weightedLaplacian` (entries identified with
`weightedOffDiag`/`weightedDiag` by `Cuts.weightedEntry_eq`), `diag(A) + τL`
(`Cuts.lengthStep_entry`, `penaltyStep_entry`), `Cuts.fluxMatrix` (defined from `fluxDiagEntry`,
`fluxOffEntry`), and an arbitrary matrix `L` with a predicate `I` for the Dirichlet solves
(`Cuts.IsDirichletSolution`, `boundaryFlux`, with `normalFlux`/`adjointFlux` identified
entrywise). The oracle's dense implementation of a flow step (`DdgOracle/Cuts/*`: cut extraction,
cut mesh, assembly, the marching loop) is outside `DdgOracle/Formulas` and is not covered.

Not covered on the oracle's side: the flattened matrices `DF` and `DFreg` (the Stage I/II
theorems are stated block by block; only `DFpiece`, of which `DFimp` is an instance by definition,
is identified with a Jacobian); `nextState`/`toward` (the safeguard is covered only through the
hypothesis "the new point lies on the pieces `σ`"); the bordered linear system used on closed
meshes; everything outside `DdgOracle/Formulas` (solver, certificate code, topology).
**Absent bounds.** At `ℝ` the oracle's bounds `lo_i`, `hi_i` are finite numbers, while the
executable's default is `lo_i = −∞`, `hi_i = +∞`. `DdgProofs/Oracle/Unbounded.lean` extends the
oracle-side theorems to that case: two predicates `Lo`, `Hi` mark the vertices whose lower resp.
upper bound is absent, the cone map there is `clampOff ∘ P_t` (`clampOff` leaves out the absent
side of `min hi (max lo ·)`), and

| theorem | statement |
|---|---|
| `Ddg.Problem.objective_le_of_zero_off` | zero state, adjoint and cone residuals (`µ_i − clampOff(P_{t_i}(φ_i + µ_i))`) ⟹ `(u, µ)` is feasible and minimizes the oracle's `objective` over all pairs feasible with the same bounds absent |
| `Ddg.Problem.unique_of_minimizers_off` | uniqueness when all `wE_i·mass_i > 0` |
| `Ddg.Problem.objective_le_add_gap_off` | the duality-gap certificate; at an absent bound the dual point must satisfy `φ_i ≤ t_i` resp. `−t_i ≤ φ_i` and the conjugate loses that term |
| `Ddg.Problem.regObjective_le_of_zero_off` | a zero of the regularized residuals minimizes the regularized objective |
| `Ddg.Problem.zero_off_of_Fpiece_eq_zero` | a solution of the oracle's piecewise linear system `Fpiece γ σ = 0` whose vertices lie on their pieces is a zero of the residuals |
| `Oracle.clamp_box_eq`, `Oracle.clampOff_eq_value` | the reduction: a finite bound that is out of the way acts like an absent one; the pieces `stateOff` are correct |

The proofs reduce to the finite theorems: the points that occur in a statement lie in a finite
box (`Ddg.Problem.boxed`), on which the oracle's own `coneResidual`, `Feasible`, `objective`,
`dualValue` for the boxed problem agree with the unbounded ones. What is **not** machine-checked
is the remark that, in IEEE arithmetic, the oracle's `clamp lo hi y = min hi (max lo y)` with
`lo = −∞` and/or `hi = +∞` evaluates to `clampOff` at a number `y` (nothing here is about
`Float`). The stopping rule for absent bounds takes the hypothesis "the affine function of `σ_i`
has the value of the point's piece" rather than an equality of canonical names.

Everything else (all of `ConeProblem`, `Site`, `TriMesh`, the dual objective, Gauss–Bonnet's
combinatorics, the Yamabe linearization) is about this package's own definitions.

The Rust code is not connected to Lean mechanically. For `crates/ddg-core`,
`DdgProofs/Evaluated.lean` transcribes the expressions the Rust functions evaluate
(`kahanArea`, `kahanCornerAngle`, `kahanCornerCot`, `angleDefectField`, `geodesicCurvatureField`)
and proves them equal to the math-level quantities; `DdgProofs/Bff/Quality.lean` and
`DdgProofs/Bff/Curve.lean` do the same for the per-triangle and per-loop expressions of
`crates/ddg-bff`, and `DdgProofs/Solver.lean` transcribes the rows of the cone solver's linear
systems. That a transcription matches the Rust source is checked by reading, and by the
differential tests against the oracle.

### Where the hypotheses of the theorems come from in the Rust code

| hypothesis (Lean) | guaranteed by |
|---|---|
| `TriMesh.Valid.ab/bc/ca` (three distinct vertices per face) | `ddg_core::TriMesh::new` (`MeshError::RepeatedVertexInFace`) |
| `TriMesh.Valid.tri` (strict triangle inequalities of the stored lengths) | `ddg_core::Geometry::from_edge_lengths` / `from_positions` (`GeometryError::DegenerateFaces`), **not** `TriMesh::new`; the test `z − (x − y) > 0` on sorted sides is exact for the stored `f64` values (`kahanArea_of_sorted_rearrangement`) |
| `TriMesh.Manifold.edge_le_two` | `TriMesh::new` (`MeshError::NonManifoldEdge`) |
| `TriMesh.Manifold.boundary_two` | `TriMesh::new` (one fan per vertex, `MeshError::NonManifoldVertex`) |
| `TriMesh.Connected` (kernel of `L`, positive definiteness of `L_II`, existence/uniqueness of the Poisson solves) | **not** `TriMesh::new`, which only counts components; `crates/ddg-optimal-cones` solves each connected component as its own problem (with the masses of the whole mesh and the common `λ`; `component-separation`, `ConeProblem.isMinimizer_sum_iff`, `kkt_sum_iff`), so the theorems that need connectedness are applied per component; `crates/ddg-bff` works per connected component / chart. (Lean's `Connected` also requires every vertex to lie in a face.) |
| `ConeProblem.WellFormed` (`wE, m, λ, wR ≥ 0`, `lo ≤ 0 ≤ hi`) | `ddg_optimal_cones::problem::validate` |
| `P.Lᵀ = P.L`, `1ᵀL = 0` on a closed mesh | `TriMesh.laplacian_isSymm`, `laplacian_vecMul_one` (proved) |
| `boundaryProblem_isMinimizer` hypotheses: `(L u)_i = s_i − Ω̂_i`, `(L φ)_i = −wE_i m_i u_i`, `φ_i ∈ ∂g_i(s_i)` at all vertices (`wE ≥ 0`, valid mesh) | `ddg_optimal_cones::certify_with_boundary` (residuals of the state and adjoint equations, and the per-site inclusion with `Free` / `Penalized` / `Convex` sites); `Discrete::with_boundary` builds `Ω̂ = (Ω on I, k on B)` and `m` normalized by the total area (`curvatures`, `lumped_mass_normalized`). The solver and the certificate check these only up to tolerance; the exact statement is for exact solutions (`certificate_enclosure` and related give the inexact version for `ConeProblem`). |
| `sum_eq_of_feasible` (`Valid`, `Manifold`, the constraint rows) | `TriMesh::new` (manifoldness), `gauss_bonnet_residual` of `certify_with_boundary` |
| `Site.convex = penalized 0 0 ⊤`, `Site.free`, `Site.boundaryPenalized wR` | `BoundaryTreatment::{Convex, Free, Penalized}` (`wR_i = +∞` is not a penalized site: it is `Site.forbidden`, i.e. `ν_i = 0`) |
| modelling step `ν_i = k_i + (L u)_i` | proved only to first order in the scaling parameter (`boundaryProblem_firstOrder_boundary`); the Rust crate additionally checks it against the generic convex solver and the double of the mesh, which are tests, not proofs |
| `ReducedOptimal`, `RedFeasible`, the sum constraint (factor-once KKT link) | the conditions that `factor_once::attempt` establishes by its working-set loop and `Dense::exact_violation`; the tolerance `tolerance·λ` is **not** exact, so the link is a statement about exact conditions, and the run-time `certify` (`optimality-system-sufficient`, `duality-gap-certificate`) remains what bounds a computed answer |
| `u = Z µ_S − w0`, `y = A⁻¹(Wu)` | `Reduced::u_from_columns` / `u_of_mu` and `Reduced::gradient` (floating point; the stored `u` is checked by `certify`, not by this theorem) |
| `∀ i, 0 < w_i` (`gram_quadratic_pos`, `sum_weight_ne_zero`) | `Reduced::check` (`Unsupported::ZeroDistortionWeight`) |
| `P.L.PosDef` (factor-once Dirichlet) | `TriMesh.coneProblem_posDef` (connected, a non-unknown vertex); the Cholesky factorization of `Reduced::new` (`FallbackReason::FactorizationFailed`) fails otherwise |
| `KernelIsConstants P.L` (factor-once closed) | `TriMesh.coneProblem_kernelIsConstants` (connected valid mesh, every vertex an unknown); connectedness: one component at a time (row `TriMesh.Connected` above); `A_p` positive definite then follows (`pinBlock_isUnit`), which the Cholesky factorization of `A_p` in `Reduced::new` also detects |
| the pin `p` | `Reduced::check`: a free site whose bounds allow the sum, or the plan with clamped sites; the Lean theorems hold for **any** index `p` (the pin need not be a free site), the code's choice only matters for feasibility of the starting point |
| `Σ_S µ_S = Σ(Ω − a_F)` | `Reduced::check` (`c`), kept by the dense method with the multiplier `ν` (`Dense::solve`) |
| `IsCutOf M' M ρ` (same faces, copies of vertices, lengths pulled back) | the construction in `ddg_bff::cut::cut_mesh` (by reading; `CutMesh::vertex_to_original` is `ρ`) |
| disk topology of a chart | checked at run time by counting (`ComponentTopology::is_disk`), which rests on the cited `genus-from-euler` |
| cut flow: `δ ≤ ½` (`clampCrossing_mem` …) | `0 ≤ min_crossing_offset < ½` is checked (`CutError::InvalidOptions`) |
| cut flow: `Θ_v ≠ 0` at every `v ∉ B` (`sum_faceCurvature`) | `ddg_core::TriMesh::new` rejects a vertex that lies in no face (`MeshError::UnreferencedVertex`), and the angle sum at a vertex of a face of a valid mesh is positive (`TriMesh.angleSum_pos`) |
| cut flow: strict triangle inequalities of the parent face (`subTriangleArea_eq_triArea`, `cotFromInner_eq_cornerCot`, `baryLength2_nonneg`) | `ddg_core::Geometry` (`GeometryError::DegenerateFaces`), as for `TriMesh.Valid.tri` |
| cut flow: `det L_II ≠ 0` / `L_II` positive definite (`yamabe_existsUnique`, `existsUnique_harmonicExtension`) | proved for a connected valid mesh with a vertex outside `I` (`TriMesh.yamabe_existsUnique`); a component of a region without boundary vertex is detected and left out (SPEC §3.4, D9c), and the Cholesky factorization fails otherwise. The cut mesh of a region satisfying `TriMesh.Valid` (positive sub-triangle areas) is by construction (`\|det\| > 0` after the clamp with `δ > 0`), not proved here |
| cut flow: `L` symmetric with zero column sums (`sum_normalFlux`, `reciprocity`) | `TriMesh.laplacian_isSymm`, `laplacian_vecMul_one`; for the weighted matrix `weightedLaplacian_isSymm`, `weightedLaplacian_vecMul_one` (proved) |
| cut flow: `A_i > 0`, `L_w` positive semidefinite, `τ ≥ 0` (`lengthStep_*`) | lumped areas of a valid mesh; `w = α_L + η ≥ 0` is validated input, then `weightedLaplacian_faceMean_posSemidef` |
| cut flow: `ℓ_p > 0`, `κ_σ ≥ 0`, symmetric couplings (`fluxMatrix_posDef` …) | dual lengths of cut segments of positive length; `fluxCoupling_nonneg`; one `κ_σ` per segment added to both of its ends (by reading `recovery.rs`) |
| cut flow: `a = r³ > 0` for the scaled energy | the smoothing uses the Laplacian `L⁺` of the edge weights `max(w_ij, 0)` (SPEC §3.6, D17), whose off-diagonal entries are `≤ 0`: `smoothedFeatureSize_bounds_positivePart`; a caller-supplied `a` is checked for `a_i > 0` |
| cut flow: `ω ≥ 0` (`rippleFactor_mem_Ioo`) | proved for the Dirichlet energy (`roughRate_dirichlet_nonneg`); **not guaranteed** for the scaled and pointwise energies, where `Q_r` is not a square |
| cut flow: the interior angles of the regions at a cut point / triple point sum to `2π` (`curve_loads_opposite`, `junction_loads_sum`) | the cut-mesh triangles around the point lie in one planar face; taken as a hypothesis, not derived from the cut mesh |
| cut flow: `Σ_i β_ip = 1`, `β ≥ 0` (`localExtension_conservation`) | the two weights `1 − t`, `t` of a cut point with `t ∈ [δ, 1 − δ]` (`clampCrossing_mem`) |
| cut flow: `IsRegion`, `IsTwoSmallest` | the scans of `regions_of` (`cut.rs`: a strict `<`, so the lowest index among equal values) and of `redistance_regions` (`regions.rs`: smallest and second smallest, equal values allowed), by reading |
| cut flow: `ℓ_c ≥ 0`, `Σℓ > 0`, `0 < q ≤ 1` (`oracle_weightedQuantile_*`, `roughening_quantile_*`) | `ℓ_c` is a sum of halves or wholes of segment lengths of the cut (`≥ 0`; positive sum for a cut of positive length), by reading `flow.rs`; `q = ROUGHENING_QUANTILE = 0.9`; a cut without points gives 0 (`oracle_weightedQuantile_nil`) |
| cut flow: `τ_max, time_scale, max_motion ≥ 0` (`robustStep_*`, `localLimit_*`) | `FlowOptions::validate`: positive and finite (`CutError::InvalidOptions`) |
| cut flow: `h_c > 0`, `max_motion·h_c > 0`, `\|r_c\| ≤ v_c` (`localLimit_*`) | `h_c` is the length of a crossed edge of a valid geometry (`local_edges`); `v_c` is `\|r_c\|` or `max(\|r_c\|, ½·max(excess, 0))` in `StepSize::local_factor` (`abs_le_localSpeed` for the oracle's `localSpeed`); by reading |
| cut flow: `τ > 0`, `max_motion·h_i > 0` (`vertexLimit_*`) | `Fixed(τ)` is validated positive, and nothing is limited with it; the rule's `τ` is positive for positive parameters (`robustStep_pos`); `h_i` is the length of an edge of a valid geometry (`vertex_limits`), by reading |
| cut flow: `ε ≥ 0` (`energy_cases`, `energyVerdict_rose_iff`) | `energy_stop` evaluates the criterion only for `energy_tolerance > 0` and more than `energy_window` steps |
| intrinsic: `IsTriangle` of the three lengths of a face (`diamond_isometric`, `unfoldApex_isometric`, `circumcenter_equidistant`, …) | every function of `formulas.rs` that takes lengths returns `None` unless `ddg_core::geometry::triangle_area` finds a positive area; a flip keeps the property in exact arithmetic (`flip_isTriangle`), and the shape limit of intrinsic SPEC N5 guards it in floating point |
| intrinsic: `DiamondConvex` (`flip_back`, `flip_area`, `flip_angle_sums`, …) | `Flip` is refused unless both `Diamond::flipped_areas` are positive, and above the shape limit (N5) |
| intrinsic: `0 ≤ φ < 2π`, `Θ > 0`, `0 ≤ θ ≤ Θ` (`signpostUpdate_*`) | stored signposts of an interior vertex are reduced (N2); `θ` is a corner angle at the vertex and `Θ` its angle sum (of the input vertex; `2π` or `π` at an inserted one); by reading |
| intrinsic: `y_up > 0 > y_down`, the apex to the left of `up → down` (`exit_cases`, `exit_ahead_*`) | `Tracer::walk` and `exit_after_entry` of `trace.rs` (the entered halfedge has its tail above and its head below the path, the apex is unfolded to its left); by reading |
| intrinsic: `b₀ + b₁ + b₂ = 1` (`cornerDistances_eq_baryLength`), `[a,b,c] ≠ 0` (`barycentric_*`) | barycentric coordinates of a point of a face; a face of positive area laid out by `triangle_layout`; by reading |
| cut flow: `φ_t = −σ` and `dE = Σ g σ ℓ` (`firstVariation_levelSet`) | the smooth theory (`cut-level-set-evolution`, `cut-cea-method`, cited); checked on closed forms (`hasDerivAt_twoCaps`, `circle_levelSet_flow`) and by the finite-difference test of SPEC §6.3 |

## Discrepancies with the papers as printed

Items 1–10 and 14 concern the cone paper, 11–13 BFF, 15–26 the cuts paper, its SPEC and the
oracle's definitions for it. (BFF Proposition 1 is stated for a disk;
`TriMesh.sum_exteriorAngle` gives `2πχ` for any valid manifold mesh.)

1. **Sign of the state equation.** Eq. 17 / B.2 print `L u = Ω − µ` with a positive semidefinite
   `L` (Sec. 5.1.3). By `TriMesh.hasDerivAt_angleDefect_conformalScale` the curvature of
   `e^{(u_i+u_j)/2} ℓ_ij` is `Ω + L u` to first order, so the paper's `u` is minus the log scale
   factor. Here: `L u = µ − Ω`, `Lᵀφ = −wE m u`; `φ`, `µ` and the cones are unchanged.
2. **First row of `F`.** Eq. 20 has `L u − Ω + µ` (consistent with Eq. 17); the last display of
   App. B.2 and the boundary-control residual of Sec. 5.3.4 have `L u − Ω − µ`.
3. **Mass factor in the regularized stage.** Eqs. 21–22 and System B.4 print `µ = P_λ(φ)/γ`.
   That is the optimality condition for the regularizer `(γ/2) Σ µ_i²` on the integrated vertex
   values, not for `(γ/2)‖µ‖²_{L²}` of Eq. 12, which gives `µ_i = (m_i/γ) P_λ(φ_i)`. The theorems
   take a per-vertex `γ_i > 0` and cover both (`γ_i = γ/m_i` is the mass-consistent choice).
4. **Initialization of Stage II** (Fig. 12): `µ ← P_λ(φ)` drops the factor `1/γ` (and the mass) of
   Eq. 21. It is only a starting point.
5. **Weights in the threshold.** App. B.2 writes `P_{w_R}`, `∂1_{w_R}` and `Proj_{w_R}` where the
   threshold is `λ·w_R` (as in App. A.2–A.3); System B.4 and the display before it write `x` for `φ`.
6. **Pre-dual with weights.** Problem A.1 prints the linear term as `∫ (Kφ) w_R dA`; from
   `G(φ) = 1_{λ w_R}(φ) − ∫ Kφ dA` (same appendix) there is no factor `w_R`. The discrete dual
   here has the term `Σ φ_i Ω_i`.
7. **Eq. 12** writes the constraint as `Δu_γ = Ω − µ_γ` (`Ω` for `K`).
8. **Bounds** (Sec. 8.2) are described only as "a simple projection". Derived here: the cone
   condition becomes `µ = clamp_[lo,hi](P_t(φ + µ))` (clamp *after* the soft threshold; valid
   because `0 ∈ [lo, hi]`), the regularized one `µ = clamp_[lo,hi]((m/γ) P_t(φ))`, and the Newton
   method gets two more active sets (clamped at `lo`, clamped at `hi`) with slope `0`.
9. **Closed meshes** are not treated in the paper's algorithm. Derived here: the additive constant
   of `u` is a free variable, optimality forces `Σ wE m u = 0` and `Σµ = ΣΩ = 2πχ`; the additive
   constant of `φ` is the multiplier of the constraint `Σµ = ΣΩ`.
10. Thesis only: Eq. (5.6) has the sign of item 2 and writes the cone row as `µ − Proj(·)` (it
    should be `µ − P(·)`, since `Proj(v) = v − P(v)`); Eq. (5.8), third row, prints
    `[0, 1_{|Rᵀp+µ|>λw}, 1 − 1_{|Rᵀp|>λw}]` for `[0, −D(p+µ), I − D(p+µ)]`.

11. **BFF, sign of the source term.** BFF Eq. 11 prints `A u = Ω − Ω̃` and the proofs of its
    Propositions 1 and 2 use `(A u)_i + h_i = Ω_i` and `Δu = Ω − Θ`, with the positive
    semidefinite cotan matrix `A` of its Sec. 4.1. Its Alg. 1 (lines 6, 8) and the reference
    implementation pass the source term `−Ω`, i.e. `A u = Ω̃ − Ω` in the interior. The latter is
    what `CLAUDE.md` fixes (`L u = µ − Ω`) and what is proved (item 1). The boundary relation
    `k̃ = k − h = k + (A u)_B` and both propositions are unaffected
    (`TriMesh.sum_exteriorAngle`, `TriMesh.sum_fiber_interior_angle_cone`).
12. **BFF, weights of the closing problem.** Eq. 19 and the text after it weight the length change
    of a boundary *edge* by `1/ℓ_ij` (`N_ii = 1/ℓ_ij`); Alg. 5 line 7 comments "`ℓ_i` are dual
    lengths". The reference implementation and `crates/ddg-bff` use the edge lengths
    (`N_d = Σ 1/ℓ` over the one or two sides of a length degree of freedom).
13. **BFF, sign in the Hilbert transform; index of the tangent angles.** Sec. 4.3.3 defines
    `h := Uᵀa`, `h_j = ½(a_k − a_i)` for consecutive boundary vertices `i, j, k` (Eq. 17), and
    solves `A b = −Uᵀa`, i.e. `(A b)_j = ½(a_i − a_k)`. That right-hand side, `½(a_{j−1} − a_{j+1})`
    with the surface on the left, is what `hilbert_transform` in `crates/ddg-bff` returns, and it
    is the gradient of the counter-clockwise area of the boundary polygon with respect to `b_j`
    (proved in `Bff/Curve.lean`; the paper's `aᵀUb` of Eq. 16 is minus that area). Alg. 6
    (lines 2–3) instead writes `h_i ← ½(a_{i+1} − a_{i−1})` and solves `A b = h`, which is the
    same system only for the opposite direction of the walk. Separately, Sec. 4.5.1 prints
    `φ_p := Σ_{i=1}^{p−1} k̃_i` for the direction of the edge leaving vertex `p`, which omits the
    turn at `p`; Alg. 5 (lines 2–4) accumulates `φ_{i,i+1} = φ_{i−1,i} + k̃_i`, which is what the
    code does (`boundaryTangent_succ`).
14. **Cone paper, Sec. 8.1** prints `r_i := 1/(|Ω_i|/M_ii + ε)` and describes it as "small in flat
    regions and large in highly curved regions"; the printed formula is the opposite. See
    `feature-size-weights` in `CITED.md`.

15. **Cuts paper: the list of differences is SPEC §9 (D1–D31).** The ones a theorem here bears
    on: D1 (sign of the face curvature: `Σ_f Ω_f = +2πχ` on a closed mesh,
    `sum_faceCurvature_closed`); D2 (sign of the speed: `firstVariation_levelSet` against
    `firstVariation_ascent`, and `hasDerivAt_twoCaps` on a closed-form example); D4 (`W_V` on the
    left of the length step: the matrix that is positive definite and whose step is the
    minimizing movement, `lengthStep_posDef`, `lengthStep_minimizes`); D8 (the paper's penalty
    rule equals the default one exactly when the penalty is constant, `penaltyStep_const_iff`);
    D11 (the edge rule is indefinite on an obtuse triangle, `edgeRuleEnergy_neg_example`); D12
    (the recovery keeps the total flux and satisfies a maximum principle, `fluxRecovery_*`); D13a,
    D13e (`cotFromInner_eq_cornerCot`, `face_omegaCancellationFree_eq`); D14 (the systems on the
    cut: `onCut_*`, and the total turning by Gauss–Bonnet); D15 (the extension is a convex
    combination that keeps the mean rate, `localExtension_*`); D16 (`norm_grad_eq_gradientNormLone`,
    `meanGradient_normalizes`); D17 (item 18); D18 (`rippleRate_*`, `sphere_equator_stable_iff`;
    the linearization is cited); D20, D21 (`DdgProofs/Cuts/Regions.lean`); D22 (the quantile,
    the robust step and the local limit: `DdgProofs/Cuts/Quantile.lean`); D23a
    (`gradientNormLone_zero`), D23c (`lengthWeightForRipples_nonpos`, `adaptiveStep_skip_omega`,
    `rippleRate_neg_of_neg`); D24a, D24e (`localLimit_cases_agree`,
    `localLimit_first_case_of_zero`); D25 (items 27–31 below, taken into SPEC); D26
    (`DdgProofs/Cuts/VertexLimit.lean`); D27, D31 (`DdgProofs/Cuts/Stop.lean`); D30
    (`pieceSlope_*`, `meanGradient_le_of_face`).
16. **SPEC §3.1, a hypothesis that is not stated.** `Σ_f Ω_f = Σ_{i ∉ B_M} (2π − Θ_i)` needs
    every vertex outside `B_M` to lie in a face (`Θ_i ≠ 0`). A vertex in no face adds `2π` on the
    right and nothing on the left; with such vertices the closed-mesh value `2πχ` (they count in
    `χ`) is also not `Σ_f Ω_f`. `ddg_core::TriMesh::new` rejects such meshes, so the identity
    holds for every mesh the crate accepts; SPEC does not say that it relies on this.
17. **SPEC §3.5, "the unique vector".** `v` is determined by `L u = −g + v` alone; that it
    vanishes on `I` is the Yamabe system, not an extra condition (`normalFlux_spec`).
18. **SPEC §3.6, positivity of the curvature-based weight** (now SPEC D17). `r⁰ > 0`, but with
    the cotan Laplacian the smoothed `r` (`(W_V + cL) r = W_V r⁰`) is only proved to stay in the
    range of `r⁰` when the off-diagonal entries of `L` are nonpositive (`smoothing_ge`); with
    negative cotangent weights the discrete maximum principle can fail and `a = r³` could be
    negative. SPEC therefore smooths with the Laplacian `L⁺` of the edge weights `max(w_ij, 0)`,
    for which the bounds hold on every mesh (`smoothedFeatureSize_bounds_positivePart`). The flux
    recovery and the systems on the cut have the maximum principle by construction (their
    off-diagonal entries are `−κ ≤ 0`).
19. **SPEC §3.6, what the adjoint is exact for.** The discrete adjoints are the exact first
    variation of the discrete objectives with respect to the *source* `g` at fixed mesh and fixed
    `I` (`scaledEnergy_perturb`, `hasDerivAt_pointwiseEnergy`). The gradient density with respect
    to the *cut* is not an exact derivative of the discrete energy; it is the discretization of
    the cited Table 1, exact on the closed-form examples.
20. **SPEC §3.12, the closed form.** `hopfLaxTriangle` is the minimum over `λ ∈ [0, 1]` also for a
    degenerate edge (`x_j = x_k`) and without any acuteness assumption; "a convex function of `λ`
    in the plane of the triangle" holds in every inner product space
    (`hopfLaxObjective_convex`). The floors `max 0 (·)` of the oracle's definition are inactive
    for side lengths of a Euclidean triangle (`hlQ_nonneg`, `baryLength2_nonneg`).

21. **SPEC §3.10, the factor of a ripple.** `1 + τ(ωk − wk²)/(1 + τwk²)` "lies in `(0, 1)` for
    every `k > ω/w` and every `τ`" needs `ω ≥ 0`: for `ω < 0` it is still `< 1` but can be
    negative (`rippleFactor_neg_example`: `−4.5` for `τ = 1`, `ω = −10`, `w = k = 1`). `ω ≥ 0` is
    proved for the Dirichlet energy only (`roughRate_dirichlet_nonneg`).
22. **SPEC §3.8, "`σ` lies between the smallest and the largest `r`".** True for every vertex
    that a ring reaches (`extension_mem_of_rings`); a vertex that no ring reaches gets `0`, which
    need not lie in that range (if all rates have one sign).
23. **SPEC §8.2, range of the crossing parameter.** `t ∈ [0, 1]`; more precisely `[0, 1)` if
    `i < j` and `(0, 1]` if `j < i` (`crossingParamRegions_mem`): the end value is possible only
    on the side whose region wins the tie.
24. **SPEC §8.5, junction loads.** The loads `2π/3 − β^X` of the three regions sum to zero at
    *every* junction whose three interior angles sum to `2π`, balanced or not
    (`junction_loads_sum`); what characterizes a balanced junction is that each of them vanishes
    (`junction_load_balanced`). SPEC §10 words it as "the loads of a balanced junction sum to
    zero".
25. **SPEC §3.10, the turning of a closed cut.** "A closed cut on a flat surface loses area at
    the exact rate `2πw`" is `Σ ℓ λ⁰ = 2πw`; it holds because the turning angles `π − β_p` are
    the boundary curvatures of the cut mesh of the region, by discrete Gauss–Bonnet
    (`TriMesh.sum_turningAngle`), with the hypotheses that the region is a disk (`χ = 1`) and
    flat; no theorem about plane polygons is needed. On a curved region the total is
    `w(2πχ − Σ_interior Ω)`.
26. **Oracle, `weightedMean`** (resolved in the oracle, revision `3884a02`). It zipped the two
    lists but divided by the sum of *all* weights, so for lists of different lengths the result
    was not a weighted mean. It now divides by the sum of the weights of the zipped pairs
    (`weightedMean_eq`) and is a weighted mean of the zipped values for any two lists
    (`weightedMean_mem`).
27. **Oracle, `sortAscending` and the order of equal values** (SPEC D25e: left as it is). The docstrings say "after all
    entries with an equal key (so that equal values stay in their original order)" and "stable
    ascending sort", and SPEC §1 says "equal values in the order of the cut points". The
    function inserts the *first* element last (`sortAscending (e :: es) = insertAscending e
    (sortAscending es)`), behind the equal keys of all later elements, so equal values come out
    in the **reverse** of their original order: `sortAscending [(1,10),(1,20)] =
    [(1,20),(1,10)]` (`sortAscending_reverses_ties`). In exact arithmetic the quantile is the
    same for every ascending arrangement (`isWeightedQuantile_scan`), so no theorem is affected;
    in floating point the running sums of a group of equal values are accumulated in another
    order than in the Rust code (a stable sort), which can only matter when the threshold
    `q·Σℓ` falls within rounding of the sum at the end of such a group.
28. **SPEC §1, "`q = 1` returns the largest value"** (SPEC D25b). This needs positive weights
    (`oracle_weightedQuantile_one`): a largest value of weight zero is not returned, because the
    sum reaches `Σℓ` before it (the oracle's docstring says so). The same happens in floating
    point when the weight of the largest values is absorbed by the running sum: the oracle's
    `Float` instance gives `weightedQuantile 1 [1, 2] [1e20, 1] = 1` (evaluated, not a theorem).
    So "accumulated in the ascending order, so that `q = 1` returns the largest value in
    floating point too" holds when no trailing weight is lost in the sum, which is the case for
    lengths of comparable size.
29. **SPEC §4, the factor of the local limit** (SPEC D25c). "The factor lies in `(0, 1]`" needs
    `max_motion · h_c > 0` (`localLimitFactor_mem_Ioc`). `min(1, max_motion · h_c/(τ v_c))` is
    not defined at `τ v_c = 0`; the oracle and the Rust code return `1` there, and since D24 the
    first case of the rule (`v_c/h_c ≤ m̂`) covers it, because `m̂ ≥ 0`
    (`localLimit_first_case_of_zero`). "The two cases agree" (§10) needs `τ ≥ 0` and `h_c > 0`
    in addition to the two inequalities it names (`localLimit_cases_agree`).
30. **SPEC §4, "ripples longer than `ℓ` grow and shorter ones decay on nine tenths of the
    cut"** (SPEC D25d). Exactly: the shorter ones decay where `ω_c ≤ ω`, at least nine tenths of the length
    (`rippleRate_neg_of_le`, `roughening_quantile_spec`); the longer ones grow where `ω_c ≥ ω`,
    which is more than one tenth of the length (`rippleRate_pos_of_ge`,
    `roughening_quantile_above`) and nine tenths only where the local rates are equal; where
    `ω_c < ω`, ripples with wavelengths between `ℓ` and `ℓ·ω/ω_c` decay as well
    (`rippleRate_neg_of_lt_longer`). §1 states it exactly ("on nine tenths of the cut, ripples
    shorter than `2πw/ω` decay").
31. **SPEC §4, the limited rate as a descent direction** (SPEC D25a; §4 now words it as below).
    A gradient density scaled point by point by factors `≥ 0` is a descent direction
    (`firstVariation_localLimit`). In the flow this is the limit at the cut points when the
    rate on the cut is the gradient density itself, which with the length step on the cut means
    no coupling along the cut: `w ≡ 0`, `r = s` (`onCut_limited_descent_of_no_coupling`). With
    coupling the unlimited rate is a descent direction for `g = s + λ⁰` (`onCut_rate_descent`:
    `dE = −rᵀ(diag ℓ + τK)r`), the limited one need not be: in
    `limited_rate_not_descent_example` (two cut points, all quantities as the capped rule with
    the robust statistics produces them) `dE = −530/11` without the limit and `+28/11` with it.
    The coupling makes one point move against its own gradient density, and the limit removes
    most of the neighbour's share that paid for it. "Rest states are unchanged" holds
    (`localLimit_sign`: `r f = 0 ⟺ r = 0`). Until D26 this theorem was also read as covering
    the length step on the mesh (where the rate on the cut is `s`); since D26 the rate is not
    limited there, the limit acts at the vertices on the whole step, and no descent statement
    is made or proved for that form.
32. **Oracle, tags.** `robustStep` and `speedInEdges` carry `[formula: cut-step-local]`; SPEC §4
    puts the step rule and `m_c` under `cut-step-cap`. The theorems about `robustStep` are
    tagged `cut-step-cap` here.
33. **SPEC §4, the limit at the vertices (D26): hypotheses, and the formula in floating point.**
    "No value changes by more than `max_motion · h_i`" needs `τ > 0` (the speed is `|Δ_i|/τ`)
    and `max_motion · h_i ≥ 0`; "the direction is kept" needs `max_motion · h_i > 0`
    (`vertexLimit_change_le`, `vertexLimit_sign`). Both hold for every accepted option
    (`robustStep_pos`). For the capped rule `τ` cancels: the rule is the clamp of `Δ_i` to
    `±max_motion · h_i` (`vertexLimit_capped_clamp`). In exact arithmetic a vertex with the
    factor `1` keeps `φ′` (`limitedValue_one`); evaluated literally in floating point,
    `φ + 1·(φ′ − φ)` differs from `φ′` in the last bits for most values (doubles:
    `0.9 + (0.001 − 0.9) = 0.0010000000000000009`). The oracle evaluates the formula at every
    vertex, the Rust code only where the factor is below `1`.
34. **SPEC §4, §8.4, `N` regions: what the limit bounds.** Half the difference of the two
    functions that are smallest at the vertex changes by at most `max_motion · h_i`
    (`vertexLimit_regions_change_le`), as §10 says. The change of each single function, and of
    its difference to a third function, is multiplied by the factor
    (`limitedValue_sub_limitedValue`, `abs_limitedValue_sub_le`) but is not bounded by
    `max_motion · h_i`; "no value changes by more than `max_motion · h_i`" (§4) is the
    two-region statement.
35. **SPEC §4, the stop rule (D27).** The three cases are exclusive only for `ε ≥ 0`
    (`energy_cases`, `energy_cases_overlap_of_neg`). SPEC does not say what `ε = 0` means: read
    literally the table stops every run whose energy has not strictly fallen over the window
    (`energyVerdict_zero_stops_iff`), while the Rust code does not evaluate the criterion for
    `ε ≤ 0`. "A caller that takes fewer than `W` steps per run is never set back" also holds
    for exactly `W` steps (the criterion needs `E_{k−W}` with `k` the number of steps minus
    one). (Taken into SPEC as D29; the rule now has the ids `cut-stop-energy` and
    `cut-run-lowest`.)
36. **Signpost paper, Alg. 8 (edge split): the two lengths are exchanged.** With `b_i` the
    barycentric weight of the corner `i`, the point on the edge is `b_i·i + b_j·j` and its
    distance from `i` is `ℓ_ij·b_j`. That is also what the paper's own Alg. 6 gives for
    `Distance(ℓ_ijk, i, (b_i, b_j, 0))` (`alg6_distance_on_edge`), and lines 6–7 of Alg. 8 use
    `(b_i, b_j, 0)` in this sense; lines 4–5 print `ℓ_ip ← ℓ_ij b_i`, `ℓ_jp ← ℓ_ij b_j`
    (`alg8_as_printed_differs`: `1/4` against `3/4`). The oracle (`edgeSplitLengths`), the Rust
    code (`la = t·ℓ` for the parameter `t` from the tail) and the intrinsic SPEC (§3, N17) have
    the right one (`edgeSplitLengths_param`).
37. **Signpost paper, Alg. 1 (signpost update): no reduction.** Sec. 3.1 has `φ : H → [0, 2π)`,
    but `φ_ik ← φ_ij + 2πθ/Θ` leaves that range whenever the two halfedges lie on different
    sides of the reference direction (in Alg. 5 and Alg. 3; not in Alg. 4, which never passes
    it: `signpost_fold_no_wrap`). The crate and the oracle subtract `2π` once; that is enough
    for `0 ≤ φ < 2π` and a corner `0 ≤ θ ≤ Θ` (`signpostUpdate_mem`), and makes the updates
    around a vertex consistent (`signpostUpdate_add`, `signpost_around_interior`). (Recorded
    as N18 of the intrinsic SPEC.)
38. **`it-cotan-weight`: the function is `2w`, and the oracle's docstring has the sign of the
    Laplacian wrong.** The oracle's `cotanWeight` and the Rust `opposite_cotan_sum` return
    `cot θ_k + cot θ_l`, twice the weight `w = ½(cot θ_k + cot θ_l)` (now stated in SPEC N17 and
    the crate's README). The oracle's docstring calls it "twice the entry of the Laplacian"; in
    the convention of `CLAUDE.md` (`L` positive semidefinite) the entry is `L_ij = −w`, so the
    function is `−2·L_ij` (`cotanWeight_eq_neg_two_mul_offDiag`). (The oracle's docstring was
    corrected after this finding.)
39. **Oracle and Rust, `cornerDistances`: `b₀` is not used.** The function returns the
    distances of `b₁v₁ + b₂v₂` (the Lean build of the oracle prints "unused variable `b0`").
    That is the point with the coordinates `(b₀, b₁, b₂)` only when they sum to `1`; the link
    to the paper's Alg. 6 needs that hypothesis (`cornerDistances_eq_baryLength`). The README of
    the crate states `d_k = |Σ_j b_j v_j − v_k|` without it.
40. **`it-ray-edge-intersection`: what the formulas intersect.** `t = y_tail/(y_tail − y_head)`
    and `x = x_tail + t(x_head − x_tail)` are the intersection of the *line* of the edge with
    the *line* of the path. That the point is on the edge needs the ends on opposite sides
    (`axisCrossingT_mem_Icc`); that it lies ahead on the ray is not in the formula: it follows
    from the orientation of the unfolded triangle (`exit_ahead_above`, `exit_ahead_below`). The
    table of the crate's README gives the two formulas without hypotheses.
41. **`it-unfold-apex`: an edge that does not have the given length.** The rustdoc says the
    given lengths are used "so errors of the layout are not fed back into the shape". Exactly:
    the result is the triangle of the three given lengths scaled by `|b − a|/ℓ_ab`
    (`unfoldApex_similar`): its angles are those of the lengths, its size follows the placed
    edge; the two new sides have the given lengths only if `|b − a| = ℓ_ab`
    (`unfoldApex_isometric`).
42. **Flip: the paper's formula and its range.** The law of cosines of the paper's Fig. 9 and the
    layout formula of the crate are the same number for every pair of nondegenerate triangles
    (`flippedLengthFromLengths_sq`). It is the distance of the two apexes in the plane; it is
    the length of the new edge only when the diamond is strictly convex, which is when the
    crate flips. The paper's remark (Sec. 4.1) that an edge which needs a flip can be flipped
    holds, ties included (`diamondConvex_of_pi_le`); in floating point the shape limit of the
    intrinsic SPEC N5 can still refuse it (N6 counts that case).
43. **Order of evaluation.** The Rust code evaluates the Delaunay excess as `(θ_k − π) + θ_l`
    and normalizes barycentric coordinates by the sum of the three areas; the oracle evaluates
    `θ_k + θ_l − π` and divides by the area of the triangle. They are equal in exact arithmetic
    (`delaunayExcess_eq'`, `barycentric_eq_normalized`) and can differ in the last bits.
44. **`CITED.md`, intrinsic triangulations.** The six entries were checked against the sources
    (2026-10-14); each has a "Source check" line. Read: Bobenko–Springborn (arXiv v3), Indermitte
    et al., Chew 1993, Shewchuk 1997 and 2002, Polthier–Schmies (author's preprint). Not read
    (paywalled): Rivin, Rippa, Fortune and Ruppert's journal paper; their data are from Crossref
    and what they are cited for is taken from the works that quote them (B–S for Rippa and
    Rivin, the Signpost paper for Fortune, Shewchuk for Ruppert). Corrections: the termination
    of flipping is Indermitte et al.'s (B–S give a second proof); `30°` is Chew's bound and
    `26.5°` is Shewchuk's; the `90°` condition on input angles for Ruppert's rule is dropped
    (only Shewchuk's `60°` was found); the pages of Fortune were added.
45. **SPEC §4, "Stopping", after D31.** The first sentence still lists "when the cut motion of
    a step … is at most `motion_tolerance`" among the reasons to stop; with D31 that is
    necessary and not sufficient (`motionCriterion_spec`), as the formula two lines below says.
    "The predicted motion is at most `τ · max |r_p|`" (§10) needs `τ ≥ 0` and, at every
    cut-adjacent vertex, weights `ℓ_p β_ip ≥ 0` that are not all zero
    (`predictedMotionOnCut_le_max`); on a network the bound is over the rates of the copies. No
    order between the cut motion and the predicted motion is claimed or proved: the cut motion
    is measured after the symmetry projection and the redistancing and also over the vertices
    that are cut-adjacent after the step. "A piece without any face that counts has `ḡ = 1`"
    (§3.9) is a rule of the code, not a property of the formula; the code also takes `1` when
    the harmonic mean is not positive (`pieceSlope`).

## Caveats on faithfulness

- All statements are about exact real arithmetic; nothing is proved about rounding.
- **The certificate.** `duality_gap_certificate` and `adjoint_gap_certificate` assume an exactly
  feasible primal point. For a computed triple the statements of section A.4′ apply instead; they
  need (a) the clipped `φ̂` (a comparison and a `min`/`max`, exact in floating point), (b) where
  `wE_i m_i = 0`, `(Lᵀφ̂)_i = 0` exactly, (c) for the forward enclosure, `µ + r` admissible, which
  at a forbidden or prescribed vertex means `r_i = 0` exactly and at a vertex on a bound that the
  residual points inward. When (c) fails only the backward statements hold
  (`backward_certificate`, `repaired_certificate`). The quantities `r`, `D`, `J` themselves are
  evaluated in floating point by the Rust code, and that evaluation is outside every theorem.
- The theorems of `DdgProofs/Solver.lean` say what an *exact* solution of the solver's linear
  system is. The solver accepts a computed solution by its componentwise backward error
  (`oettli_prager` says what that number means); the effect of that error on the iterate is not
  analysed. The solver's final answer is judged by the certificate, not by these theorems.
- The Rust classifier `Site::piece` (pieces closed at their kinks, a free site labelled `Pos`
  throughout) is not transcribed; what is proved is, for each kind of row, the test under which a
  vertex lies on the piece of its row (`Site.stationary_of_pos` …), and the package's own
  classifier `Site.pieceReg`.
- A forbidden vertex (`wR = +∞`) is modelled as the constraint `µ_i = 0` (`Site.fixed 0`), not
  as an infinite weight.
- The Dirichlet energy of an affine function on a triangle is taken to be `Area·‖∇f‖²`; the
  integral itself is not formalized.
- Gauss–Bonnet's "boundary vertex" is an end point of an edge that lies in exactly one face;
  vertices that belong to no face count as interior (with `Ω = 2π`) and contribute to `χ`.
- The Newton theorems say that a step which keeps the active sets yields an exact solution; they
  do not say that such a step occurs.
- BFF: the theorems cover the linear algebra and the bookkeeping (existence and uniqueness of the
  solves, angle sums, closing of the curve, Euler characteristic). They do not say that the
  result is conformal, injective or close to the smooth map: quality is measured by the
  distortion formulas (whose definitions are proved to be what they are called) and reported.
  The rank-deficient branch of the closing problem (dropped directions), the optional diagonal
  shift of the reference implementation, and all tolerances are not covered.
- `IsCutOf` describes a cut by its result (same faces, vertex copies); that `cut_mesh` produces
  such a mesh, and that the copies are the sectors around each vertex, is by reading the code.
- `atan2 y x` is defined as the argument of the complex number `x + iy`; the Kahan arrangements of
  `DdgProofs/Evaluated.lean` are proved equal to the plain formulas as real numbers, which says
  nothing about their (better) rounding behaviour.
- **Division by zero.** In Lean `x/0 = 0` for real numbers, while IEEE arithmetic gives `±∞` or
  NaN. A statement about a quotient that carries no hypothesis on its denominator (for example
  `smoothedSpeed_eq_convex`, `gradScaledD_nodal`, `cotFromInner_eq_cornerCot` for `det = 0`,
  `featureSizeRaw_zero`) holds at a zero denominator only by this convention and says nothing
  about the floating-point value there; wherever a zero denominator would change the meaning the
  hypothesis is stated (`φ_a ≠ φ_b`, `Θ_v ≠ 0`, `ℓ_p > 0`, `m > 0`, `r ≠ 0`, …).
- **Cut flow.** The theorems are per formula and per linear system. They do not say that a flow
  step decreases the discrete energy `E(φ)` of SPEC §6, that the cut converges, or that the
  discrete gradient density is consistent on a mesh; those are tested, not proved. The
  first-variation statements (`firstVariation_*`) take the formula `dE = Σ g σ ℓ` and the relation
  `φ_t = −σ` as hypotheses. The closed-form theorems (`DdgProofs/Cuts/Analytic.lean`) are
  one-variable calculus on explicit functions: they verify the smooth formulas on radially
  symmetric examples (the Laplace–Beltrami operator of the sphere enters only as the ordinary
  differential equation `(sin ψ·u′)′ = K sin ψ`; for the disk the planar Laplacian is also
  computed in Cartesian coordinates), not the smooth theory.
- **Cut flow, combinatorics.** Which faces are cut, the choice of the quad diagonal, the regions
  and their components, the cut mesh as a `TriMesh`, and the marching order of the redistancing
  are not modelled. The sub-triangle statements are for one face and arbitrary barycentric rows;
  `cutFace_det_sum_*` covers the two splittings of SPEC §3.3.
- **Cut flow, the flux recovery.** Proved: well-posedness, conservation and the maximum principle
  for any symmetric nonnegative couplings. Not proved: any accuracy statement (SPEC D12 reports
  measurements), and that `κ_σ = c²h_f²/|pq|` is a good choice.
- **Cut flow, default stages and `N` regions.** The theorems about the systems on the cut are
  for arbitrary symmetric nonnegative couplings on an index set; that the Rust code assembles
  these couplings from the segments of the cut (one chain per region, through triple points) is
  by reading. The length force `Θ` is taken as given, except for its total on a closed cut with
  a constant weight. The angle sums `2π` around a cut point or a triple point are hypotheses.
  The roughening rate and the junction rule rest on cited smooth results
  (`cut-ripple-linearization`, `cut-network-first-variation`); SPEC itself reports the junction
  rule as empirical (D20). `IsRegion`, `IsTwoSmallest` and `regionRate` are this package's
  formulations of rules that SPEC states in words. Nothing is proved about the marching per
  function of the redistancing, the numbering of nodes, or the curves between junctions.
- **Intrinsic triangulations.** The theorems are about scalar closed forms and about points of
  the plane. Not modelled: the halfedge mesh, its signposts as stored data and their update by
  the operations (the theorems say what one update and one flip do to the numbers); that the
  unfolding of a path of triangles composes (each step is covered, the path is not); the
  tracing query as a whole, edge traces and their checks, point queries, the overlay, the
  refinement loop, the simplicial export. Moving a vertex and the optimal Delaunay update are
  covered for one planar star with fixed connectivity and signed areas (that the star of an
  intrinsic vertex unfolds into one is SPEC N19); the iteration, the shortened step, the flips
  between iterations and the interpolation-error meaning of the energy are not (cited). The Delaunay *triangulation* never appears: what is
  proved is local (one edge and its two triangles), and the global facts are cited. Strict
  convexity is the oracle's `flipConvex`; the shape limit of N5, the tolerance of the Delaunay
  test and of the traces (N6, N8) and every other tolerance are outside. A corner angle is
  `cornerAngle` of three lengths (`arccos` of the law of cosines); the oracle and the Rust code
  evaluate it by a half-angle formula, which is the same number for a nondegenerate triangle
  (`Oracle.cornerAngle_eq`), so theorems stated with `cornerAngle` need `IsTriangle` to be
  about the oracle's `oppositeAngles`.
- **Cut flow, the limit at the vertices and the stop rule (D26, D27).** The limit is covered as
  a function of `φ(i)`, `φ′(i)`, `h_i`, `ω_i`, `τ`; how `h_i` and `ω_i` are collected
  (`vertex_limits`), which vertices are cut-adjacent and counted, which function is `b`, and the
  order of the stages are not modelled. That the paper's stages keep the form `[φ, −φ]` is
  proved stage by stage (`explicitStep_neg`, `lengthStep_neg`, `vertexLimit_two_regions`,
  `classMean_neg`), not for the composed step. Nothing is proved about the energy after a
  limited step. The stop rule is two functions of this package written after `energy_stop` and
  `run_outcome`; that the run loop feeds them the energies of the states in order, and holds
  the earliest lowest state, is by reading. NaN energies are outside the statements. The
  motion criterion (D31) is covered as a function of three numbers, and the predicted motion as
  a function of the weights and the rates at the cut-adjacent vertices; which vertices are
  cut-adjacent, which rate (before the local limit and before `ḡ`) the code passes, and the
  cut motion itself are by reading. That a step with all rates zero moves nothing is not
  proved (the symmetry projection and the redistancing are not modelled); the theorems take
  the cut motion as given.
- **Cut flow, step statistics (D22–D24).** The quantile, the step rule and the factor of the
  local limit are covered as functions of their arguments. Not modelled: which crossed edges
  `h_c` ranges over (`local_edges`) and how `ℓ_c` is assembled; which speeds enter the two `m̂`
  of the length step on the mesh and on a network (the theorem is `robustStep_max` for any two
  numbers); the copies of triple points; the count of limited cut points; `Fixed` (nothing is
  limited, nothing to prove); "a piece without any face that counts has `ḡ = 1`" (the theorems
  on `ḡ` need `ΣL > 0`). Everything about the quantile in floating point (the order of the
  running sums, equal values, `v_c/h_c` "equal to `m̂` to the last bit", the rounding of `τ`)
  is outside the theorems; the Rust test `v_c/h_c ≤ m̂` before the comparison of `τ v_c` with
  `max_motion · h_c` changes nothing in exact arithmetic (`localLimit_twoCase_eq`). The descent
  statements are first variations in the model `dE = Σ g σ ℓ`; nothing is proved about the
  energy after a step.
