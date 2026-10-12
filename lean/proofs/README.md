# lean/proofs — machine-checked derivations of the formulas

Lean 4 / Mathlib proofs of the discrete formulas used by the `ddg-tools` implementation of
*Optimal Cone Singularities for Conformal Flattening* (Soliman, Slepčev, Crane, ACM TOG 2018;
"the paper"). The ground truth is the **canonical discrete problem** and the math conventions of
the repository's `CLAUDE.md`; where the paper is silent or inconsistent the statements here are
derived from that problem (see [Discrepancies](#discrepancies-with-the-paper-as-printed)).

Results that are used but not proved here (the paper's smooth theory, convergence of the Newton
method, …) are listed with exact references in [`CITED.md`](CITED.md).

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
../../scripts/check-formulas.sh --build # the same, plus build and axiom audit
```

There is no `sorry` and no `axiom` declaration in this package. `DdgProofs/Audit.lean` runs
`#print axioms` on 89 theorems (all the main ones); every one reports exactly
`[propext, Classical.choice, Quot.sound]`.

## Layout

| File | Content |
|---|---|
| `DdgProofs/SoftThreshold.lean` | soft threshold, clamp, subdifferential of `t·\|x\|` on an interval (item A.1) |
| `DdgProofs/Problem.lean` | the canonical discrete problem: `Site`, `ConeProblem`, objective, feasibility, optimality system |
| `DdgProofs/Optimality.lean` | sufficiency, quadratic growth, uniqueness, necessity (invertible `L`), guarantee of paper Sec. 4.5 (A.2, A.3) |
| `DdgProofs/Duality.lean` | discrete dual, Fenchel–Young, weak/strong duality, duality-gap certificate (A.4) |
| `DdgProofs/Newton.lean` | residual maps, semismooth differentials, stopping rule, regularized problem (A.5) |
| `DdgProofs/Closed.lean` | closed meshes, optimal additive constant, area normalization of `λ` (A.6, A.7) |
| `DdgProofs/Triangle.lean` | corner angle, cotangent, Heron area, angle sum, scaling, link to Mathlib's Euclidean geometry (B.8) |
| `DdgProofs/Cotan.lean` | cotan formula, cotan Laplacian, lumped mass (B.9, B.10) |
| `DdgProofs/GaussBonnet.lean` | angle defect, boundary curvature, discrete Gauss–Bonnet (B.11) |
| `DdgProofs/Yamabe.lean` | linearized discrete Yamabe equation (B.12) |
| `DdgProofs/Evaluated.lean` | the same formulas in the arrangements `crates/ddg-core` evaluates (Kahan's area, `atan2` corner angle, curvature fields, Laplacian entries) |
| `DdgProofs/Oracle/*.lean` | the real-number instance of the oracle's `Ddg.Scalar` and the theorems restated about the oracle's definitions |
| `DdgProofs/Audit.lean` | axiom audit (not imported by the library) |

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
| `newton-stopping-rule` | Stage II: `F(u, φ, µ) = (L u − µ + Ω, Lᵀφ + w∘u, µ − clamp(P_t(φ + µ)))`, `DF·(v, q, ν) = (L v − ν, Lᵀq + w∘v, ν − D(q + ν))` with `D_i` the slope of the piece of `φ_i + µ_i`. If `DF·(v, q, ν) = −F(u, φ, µ)` and every `φ_i + q_i + µ_i + ν_i` lies on the same piece as `φ_i + µ_i`, then `F(u + v, φ + q, µ + ν) = 0` and `(u + v, µ + ν)` is a global minimizer. Stage I: the same for `F_γ(u, φ) = (L u + Ω − clamp((m/γ)P_t(φ)), Lᵀφ + w∘u)`, `DF_γ·(v, q) = (L v − D_γ q, Lᵀq + w∘v)`, pieces of `φ_i`. | `paper_stageII_stopping_rule`, `paper_stageI_stopping_rule`, `ConeProblem.isMinimizer_of_newton_step`, `isRegMinimizer_of_newton_step`, `Ddg.Problem.Fimp_eq_zero_of_newton_step` |
| `regularized-optimality-sufficient` | If `L u = µ − Ω`, `Lᵀφ = −w∘u` and `φ_i − γ_i µ_i ∈ ∂(t_i\|·\| + ι)(µ_i)` (`γ_i ≥ 0`), then `(u, µ)` minimizes `J(u, µ) + Σ (γ_i/2) µ_i²` over all feasible pairs. (`γ_i = γ/m_i` is the mass-consistent regularization.) | `ConeProblem.isRegMinimizer_of_stationary`, `isRegMinimizer_of_residualReg_eq_zero` |
| `regularized-suboptimality` | If `(u_γ, µ_γ)` minimizes the regularized objective and `(u, µ)` is feasible: `J(u_γ, µ_γ) ≤ J(u, µ) + Σ (γ_i/2) µ_i²`. | `ConeProblem.IsRegMinimizer.objective_le` |
| `closed-mesh-cone-sum` | If `1ᵀL = 0` and `L u = µ − Ω` then `Σ µ_i = Σ Ω_i`; on a closed valid manifold mesh this is `2πχ`. | `ConeProblem.sum_cone_eq_sum_defect`, `TriMesh.sum_cone_eq_two_pi_eulerChar` |
| `closed-mesh-optimal-constant` | `E(u + c) = E(u) + c·Σ w_i u_i + (c²/2)·Σ w_i`; if `Σ w_i > 0` the unique minimizing constant is `c₀ = −Σ w_i u_i / Σ w_i`. If `L·1 = 0`, every minimizer (and every solution of the optimality system) has `Σ w_i u_i = 0`. | `ConeProblem.distortion_add_const`, `distortion_add_const_eq`, `KKT.weighted_mean_eq_zero`, `IsMinimizer.weighted_mean_eq_zero` |
| `lambda-area-normalization` | For `c > 0`: the minimizers with masses `c·m` and parameter `λ` are the minimizers with masses `m` and parameter `λ/c`. With masses divided by their sum the problem is unchanged by `m ↦ c·m`. | `ConeProblem.isMinimizer_scaleMass_iff`, `normalizeMass_scaleMass` |

## Theorems

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
| `newton-stopping-rule` | `Ddg.Problem.F_eq_zero_of_newton_step`, `Ddg.Problem.objective_le_of_newton_step` | Stage II on the oracle's `stateResidual`, `adjointResidual`, `coneResidual`, `state`. | Sec. 5.3.3 | oracle |
| `newton-stopping-rule` | `Ddg.Problem.Freg_eq_zero_of_newton_step` | Stage I on the oracle's `Freg`, `regState`. | Sec. 5.3.3 | oracle |
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
| `lambda-area-normalization` | `ConeProblem.normalizeMass_scaleMass` | With masses normalized to sum 1 the problem does not change under scaling. | Sec. 6.2 | own |
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
| `discrete-yamabe-linearization` | `TriMesh.hasDerivAt_angleDefect_conformalScale`, `…_add` | The angle defect of the scaled metric is `Ω + L u + o(u)`; the same at any `u` with the Laplacian of the scaled metric. | Eq. 17 | own |

This confirms the sign convention of `CLAUDE.md` by proof: with a positive semidefinite `L` and
`u` the true log scale factor, the curvature after scaling is `Ω + L u`, so prescribing `µ` gives
`L u = µ − Ω`, not the paper's `L u = Ω − µ`.

## Which definitions are shared with the code that runs

The oracle (`lean/oracle`, `DdgOracle/Formulas/*`) defines the formulas once for a generic scalar;
`DdgProofs/Oracle/Scalar.lean` provides `instance : Ddg.Scalar ℝ`. Theorems **about the oracle's
definitions** (at `ℝ`, i.e. exact arithmetic):

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
At `ℝ` the oracle's bounds are finite numbers, so its unbounded default (`±∞` in `Float`) is not
an instance of these statements; the unbounded problem is covered by the **own** theorems
(`Site.penalized wR ⊥ ⊤`).

Everything else (all of `ConeProblem`, `Site`, `TriMesh`, the dual objective, Gauss–Bonnet's
combinatorics, the Yamabe linearization) is about this package's own definitions.

The Rust code is not connected to Lean mechanically. For `crates/ddg-core`,
`DdgProofs/Evaluated.lean` transcribes the expressions the Rust functions evaluate
(`kahanArea`, `kahanCornerAngle`, `kahanCornerCot`, `angleDefectField`, `geodesicCurvatureField`)
and proves them equal to the math-level quantities; that the transcription matches the Rust source
is checked by reading, and by the differential tests against the oracle.

## Discrepancies with the paper as printed

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

## Caveats on faithfulness

- All statements are about exact real arithmetic. The certificate theorems assume exact primal
  feasibility (`L u = µ − Ω`) and an exact adjoint equation; nothing is proved about rounding. A
  dual point that satisfies the adjoint equation by construction is `ũ := −Lᵀφ/(wE m)`
  (`adjoint_gap_certificate`).
- At a vertex with an infinite bound the dual point must satisfy `φ_i ≤ λ wR_i` (no upper bound) /
  `−λ wR_i ≤ φ_i` (no lower bound) exactly; a computed `φ` generally needs to be scaled or clipped
  first.
- A forbidden vertex (`wR = +∞`) is modelled as the constraint `µ_i = 0` (`Site.fixed 0`), not
  as an infinite weight.
- The Dirichlet energy of an affine function on a triangle is taken to be `Area·‖∇f‖²`; the
  integral itself is not formalized.
- Gauss–Bonnet's "boundary vertex" is an end point of an edge that lies in exactly one face;
  vertices that belong to no face count as interior (with `Ω = 2π`) and contribute to `χ`.
- The Newton theorems say that a step which keeps the active sets yields an exact solution; they
  do not say that such a step occurs.
- `atan2 y x` is defined as the argument of the complex number `x + iy`; the Kahan arrangements of
  `DdgProofs/Evaluated.lean` are proved equal to the plain formulas as real numbers, which says
  nothing about their (better) rounding behaviour.
