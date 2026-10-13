# lean/proofs — machine-checked derivations of the formulas

Lean 4 / Mathlib proofs of the discrete formulas used by `ddg-tools`:

* the canonical cone problem and its solver, `crates/ddg-optimal-cones` — *Optimal Cone
  Singularities for Conformal Flattening* (Soliman, Slepčev, Crane, ACM TOG 2018; "the paper");
* the flattener, `crates/ddg-bff` — *Boundary First Flattening* (Sawhney, Crane, ACM TOG 2017;
  "BFF");
* the operators of `crates/ddg-core` both use.

The ground truth is the **canonical discrete problem** and the math conventions of the
repository's `CLAUDE.md`; where a paper is silent or inconsistent the statements here are derived
from that problem (see [Discrepancies](#discrepancies-with-the-papers-as-printed)).

**Scope in one paragraph.** Every theorem is a statement in exact real arithmetic about a formula;
rounding is not covered anywhere, and the Rust code is tied to the statements by reading, by the
`[formula: id]` tags and by differential tests, not mechanically. Results that are used but not
proved here (the smooth theory, convergence of the Newton method, three topological facts, …) are
listed with exact references in [`CITED.md`](CITED.md), which also states precisely what a
run-time certificate of the cone solver does and does not establish.

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
../../scripts/check-formulas.sh --audit-all   # the same, plus the axioms of every declaration (minutes)
```

There is no `sorry` and no `axiom` declaration in this package. `DdgProofs/Audit.lean` runs
`#print axioms` on the main theorems, by name (the count is printed by `check-formulas.sh --build`);
every one reports exactly `[propext, Classical.choice, Quot.sound]`.
`lake env lean DdgProofs/AuditAll.lean` (several minutes) collects the axioms of *every*
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

## Theorems

Sections A and B below index the theorems about the canonical problem and the operators. The
theorems about the cone solver (`DdgProofs/Solver.lean`) and about Boundary First Flattening
(`DdgProofs/Bff/*.lean`) are indexed by the two id tables just above; they are all about this
package's own definitions, which transcribe the Rust expressions where the statement says "as
the code evaluates it".

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
  `Fimp`, `pieceConeResidual`, `Fpiece`, `DFpiece`, `objective`.

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
| `TriMesh.Connected` (kernel of `L`, positive definiteness of `L_II`, existence/uniqueness of the Poisson solves) | **not** `TriMesh::new`, which only counts components; `crates/ddg-optimal-cones` requires one component (`ConeError::NotConnected`), `crates/ddg-bff` works per connected component / chart. (Lean's `Connected` also requires every vertex to lie in a face.) |
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
| `KernelIsConstants P.L` (factor-once closed) | `TriMesh.coneProblem_kernelIsConstants` (connected valid mesh, every vertex an unknown); connectedness is `ConeError::NotConnected`; `A_p` positive definite then follows (`pinBlock_isUnit`), which the Cholesky factorization of `A_p` in `Reduced::new` also detects |
| the pin `p` | `Reduced::check`: a free site whose bounds allow the sum, or the plan with clamped sites; the Lean theorems hold for **any** index `p` (the pin need not be a free site), the code's choice only matters for feasibility of the starting point |
| `Σ_S µ_S = Σ(Ω − a_F)` | `Reduced::check` (`c`), kept by the dense method with the multiplier `ν` (`Dense::solve`) |
| `IsCutOf M' M ρ` (same faces, copies of vertices, lengths pulled back) | the construction in `ddg_bff::cut::cut_mesh` (by reading; `CutMesh::vertex_to_original` is `ρ`) |
| disk topology of a chart | checked at run time by counting (`ComponentTopology::is_disk`), which rests on the cited `genus-from-euler` |

## Discrepancies with the papers as printed

Items 1–10 and 14 concern the cone paper, 11–13 BFF. (BFF Proposition 1 is stated for a disk;
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
