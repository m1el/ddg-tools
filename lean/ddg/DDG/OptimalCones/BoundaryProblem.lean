import DDG.OptimalCones.MeshProblem
import DDG.Conformal.Yamabe
import DDG.Conformal.Cherrier
import DDG.Laplace.Kernel
import DDG.Laplace.Cotan

/-!
# The mesh problem with boundary unknowns

`crates/ddg-optimal-cones` (`Discrete::with_boundary`, `solve_with_boundary`,
`certify_with_boundary`) optimizes the boundary values of `u` as well. Let `B` be the boundary
vertices, `I` the interior ones, `Ω_i` the angle defect (`i ∈ I`) and `k_i = π − Σ(corner angles)`
the geodesic curvature of the boundary (`i ∈ B`). Put `Ω̂_i = Ω_i` on `I` and `Ω̂_i = k_i` on `B`
(`TriMesh.hatCurvature`). The unknowns are `u ∈ ℝ^V` and `s = (µ on I, ν on B)`, with one site
`g_i` per vertex (`λ·wR_i·|·|` with bounds, as in `Site`):

```
minimize     ½ Σ_V wE_i m_i u_i²  +  λ Σ_V cost_i(s_i)
subject to   (L u)_i = s_i − Ω̂_i     for all vertices i,   s_i ∈ feasible_i.
```

This is `ConeProblem V` for `L` the full cotan Laplacian (`TriMesh.boundaryProblem`), so the
results of `DDG.OptimalCones.Optimality` and `DDG.OptimalCones.Closed` apply: the full Laplacian annihilates
constants, so `u` is determined only up to an additive constant, which is optimized too.

* `boundaryProblem_isMinimizer`: the optimality system is sufficient for a global minimum;
* `sum_eq_of_feasible`: `Σ s = 2πχ` for every feasible point;
* the three boundary sites: `Site.stationary_free_iff` (`φ = 0`), `Site.stationary_penalized_iff`
  (soft threshold), `Site.stationary_convex_iff` (`ν ≥ 0`, `φ ≤ 0`, `φ ν = 0`); the convex site is
  `Site.penalized 0 0 ⊤` of the existing `Site` type, no extension was needed;
* `hasDerivAt_hatCurvature_conformalScale`, `boundaryProblem_firstOrder_boundary`: the modelling
  step, `ν_i = k_i + (L u)_i` is the first-order boundary curvature of the scaled mesh.
-/

noncomputable section

namespace DDG

open Matrix Real

namespace Site

/-- The convex boundary site: `ν ≥ 0` at no cost (`g` is the indicator of `[0, ∞)`). It is the
penalized site with weight `0` and lower bound `0`, so it is well-formed and `Site` needed no
extension. -/
def convex : Site := penalized 0 0 ⊤

/-- The penalized boundary site: cost `λ·wR·|ν|`, no bounds. -/
def boundaryPenalized (wR : ℝ) : Site := penalized wR ⊥ ⊤

theorem convex_wellFormed : convex.WellFormed :=
  ⟨le_rfl, EReal.coe_zero.symm.le, le_top⟩

theorem boundaryPenalized_wellFormed {wR : ℝ} (h : 0 ≤ wR) : (boundaryPenalized wR).WellFormed :=
  ⟨h, bot_le, le_top⟩

theorem convex_feasible : convex.feasible = {x | 0 ≤ x} := by
  ext x
  simp [convex, feasible, boundSet, EReal.coe_nonneg]

/-- **Free boundary site** (`g = 0`): the optimality inclusion is `φ_i = 0`, the natural
boundary condition.

[formula: boundary-optimality-system] -/
theorem stationary_boundary_free_iff {lam : ℝ} (hlam : 0 ≤ lam) (φ ν : ℝ) :
    free.Stationary lam φ ν ↔ φ = 0 :=
  stationary_free_iff hlam φ ν

/-- **Penalized boundary site** (`g = λ·wR·|·|`): the inclusion `φ ∈ ∂g(ν)` is the soft-threshold
fixed point `ν = P_{λ wR}(φ + ν)`, equivalently `|φ| ≤ λ wR`, `φ = λ wR` if `ν > 0`,
`φ = −λ wR` if `ν < 0`.

[formula: boundary-optimality-system] -/
theorem stationary_penalized_iff {lam : ℝ} (hlam : 0 ≤ lam) {wR : ℝ} (hwR : 0 ≤ wR) (φ ν : ℝ) :
    (boundaryPenalized wR).Stationary lam φ ν ↔ ν = softThreshold (lam * wR) (φ + ν) ∧
      (|φ| ≤ lam * wR ∧ (0 < ν → φ = lam * wR) ∧ (ν < 0 → φ = -(lam * wR))) := by
  have h1 := stationary_iff_eq_prox hlam (boundaryPenalized_wellFormed hwR) φ ν
  have hp : (boundaryPenalized wR).prox lam (φ + ν) = softThreshold (lam * wR) (φ + ν) := by
    simp [boundaryPenalized, prox, clampE_bot_top]
  rw [h1, hp]
  exact ⟨fun h => ⟨h, (eq_softThreshold_add_iff_cases (mul_nonneg hlam hwR) φ ν).mp h⟩,
    fun h => h.1⟩

/-- **Convex boundary site** (`g` the indicator of `[0, ∞)`): the inclusion `φ ∈ ∂g(ν)` is
`ν ≥ 0`, `φ ≤ 0` and `φ·ν = 0` (complementarity), for every `λ`.

[formula: boundary-optimality-system] -/
theorem stationary_convex_iff (lam φ ν : ℝ) :
    convex.Stationary lam φ ν ↔ 0 ≤ ν ∧ φ ≤ 0 ∧ φ * ν = 0 := by
  have hc : ∀ x, lam * convex.cost x = 0 := fun x => by simp [convex, cost]
  simp only [Stationary, IsSubgradientWithinAt, hc, convex_feasible,
    zero_add]
  constructor
  · rintro ⟨hν', h⟩
    have hν : 0 ≤ ν := by simpa using hν'
    have h0 := h 0 (by simp)
    have h2 := h (ν + 1) (by simp; linarith)
    refine ⟨hν, ?_, ?_⟩
    · linarith
    · nlinarith
  · rintro ⟨hν, hφ, hc⟩
    refine ⟨by simpa using hν, fun y hy => ?_⟩
    have hy' : 0 ≤ y := by simpa using hy
    nlinarith

end Site

namespace TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype V] [Fintype F] (M : TriMesh V F)

/-- The curvature `Ω̂`: the angle defect `Ω_i` at an interior vertex, the boundary geodesic
curvature `k_i = π − Σ(corner angles)` at a boundary vertex. -/
def hatCurvature (i : V) : ℝ :=
  if i ∈ M.boundaryVertices then M.geodesicCurvature i else M.angleDefect i

/-- **The mesh problem with boundary unknowns.** Every vertex is an unknown (`u ∈ ℝ^V`,
`s ∈ ℝ^V`; `s = µ` on interior vertices and `ν` on boundary vertices), `L` is the full cotan
Laplacian, `Ω̂ = hatCurvature`, `m` the lumped masses divided by the total area, and `site i`
is the site of `s_i` (use `Site.free`, `Site.boundaryPenalized wR`, `Site.convex` at boundary
vertices).

[formula: boundary-optimality-system] -/
def boundaryProblem (wE : V → ℝ) (lam : ℝ) (site : V → Site) : ConeProblem V where
  L := M.laplacian
  Ω := M.hatCurvature
  m := fun i => M.mass i / M.totalArea
  wE := wE
  lam := lam
  site := site

variable (wE : V → ℝ) (lam : ℝ) (site : V → Site)

omit [Fintype V] in
theorem mass_div_totalArea_nonneg (hM : M.Valid) (i : V) : 0 ≤ M.mass i / M.totalArea := by
  have hA : 0 ≤ M.totalArea := Finset.sum_nonneg fun f _ => (triArea_pos (hM.tri f)).le
  have hm : 0 ≤ M.mass i := by
    unfold mass vertexArea
    refine div_nonneg (Finset.sum_nonneg fun f _ => ?_) (by norm_num)
    split_ifs
    · exact (triArea_pos (hM.tri f)).le
    · exact le_rfl
  exact div_nonneg hm hA

theorem boundaryProblem_wellFormed (hM : M.Valid) (hwE : ∀ v, 0 ≤ wE v) (hlam : 0 ≤ lam)
    (hsite : ∀ v, (site v).WellFormed) : (M.boundaryProblem wE lam site).WellFormed where
  wE_nonneg := hwE
  m_nonneg := M.mass_div_totalArea_nonneg hM
  lam_nonneg := hlam
  site := hsite

/-- **Feasibility is the constraint of the statement**: `(L u)_i = s_i − Ω̂_i` at every vertex
and `s_i` admissible at its site.

[formula: boundary-optimality-system] -/
theorem boundaryProblem_feasible_iff (u s : V → ℝ) :
    (M.boundaryProblem wE lam site).Feasible u s ↔
      (∀ i, (M.laplacian *ᵥ u) i = s i - M.hatCurvature i) ∧ ∀ i, s i ∈ (site i).feasible := by
  unfold ConeProblem.Feasible
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨fun i => congrFun h1 i, h2⟩
  · rintro ⟨h1, h2⟩
    exact ⟨funext h1, h2⟩

/-- **The objective**: `½ Σ_V wE_i (m_i/A) u_i² + λ Σ_V cost_i(s_i)`.

[formula: boundary-optimality-system] -/
theorem boundaryProblem_objective (u s : V → ℝ) :
    (M.boundaryProblem wE lam site).objective u s =
      1 / 2 * ∑ i, wE i * (M.mass i / M.totalArea) * u i ^ 2 +
        lam * ∑ i, (site i).cost (s i) := by
  simp [ConeProblem.objective, ConeProblem.distortion, ConeProblem.coneCost, boundaryProblem]

/-- **The optimality system of the problem with boundary unknowns is sufficient for a global
minimum.** Let `(L u)_i = s_i − Ω̂_i`, `(L φ)_i = −wE_i (m_i/A) u_i` and `φ_i ∈ ∂g_i(s_i)` (the
site's `Stationary`) hold at every vertex. Then `(u, s)` is feasible and no feasible `(v, ν)` has
a smaller objective; in particular not `(v + c, ν)` for any solution `v` of `L v = ν − Ω̂` and any
constant `c`: the additive constant of `u` is optimized too.

Hypotheses: `M` valid (positive areas, so the masses are nonnegative), `wE ≥ 0`. Nothing is
assumed about `λ` or the sites (the inclusion already contains feasibility).

[formula: boundary-optimality-system] -/
theorem boundaryProblem_isMinimizer (hM : M.Valid) (hwE : ∀ v, 0 ≤ wE v) {u φ s : V → ℝ}
    (hstate : ∀ i, (M.laplacian *ᵥ u) i = s i - M.hatCurvature i)
    (hadj : ∀ i, (M.laplacian *ᵥ φ) i = -(wE i * (M.mass i / M.totalArea) * u i))
    (hincl : ∀ i, (site i).Stationary lam (φ i) (s i)) :
    (M.boundaryProblem wE lam site).IsMinimizer u s ∧
      ∀ (v ν : V → ℝ), (∀ i, (M.laplacian *ᵥ v) i = ν i - M.hatCurvature i) →
        (∀ i, ν i ∈ (site i).feasible) → ∀ c : ℝ,
        (M.boundaryProblem wE lam site).objective u s ≤
          (M.boundaryProblem wE lam site).objective (fun i => v i + c) ν := by
  have hK : (M.boundaryProblem wE lam site).KKT u φ s := by
    refine ⟨funext hstate, ?_, hincl⟩
    have hT : (M.boundaryProblem wE lam site).Lᵀ = M.laplacian := M.laplacian_isSymm.eq
    rw [hT]
    exact funext hadj
  have hm := M.mass_div_totalArea_nonneg hM
  refine ⟨hK.isMinimizer hwE hm, fun v ν hv hν c => ?_⟩
  exact hK.minimizes_closed M.laplacian_mulVec_one hwE hm (funext hv) hν c

/-- The curvatures `Ω̂` sum to `2πχ` (discrete Gauss–Bonnet, `Σ_I Ω + Σ_B k = 2πχ`). -/
theorem sum_hatCurvature (hM : M.Valid) (hman : M.Manifold) :
    ∑ i, M.hatCurvature i = 2 * π * (M.eulerChar : ℝ) := by
  rw [← Finset.sum_compl_add_sum M.boundaryVertices, ← M.gauss_bonnet hM hman]
  congr 1
  · exact Finset.sum_congr rfl fun i hi => by simp [hatCurvature, Finset.mem_compl.mp hi]
  · exact Finset.sum_congr rfl fun i hi => by simp [hatCurvature, hi]

/-- **The total curvature is `2πχ`.** For every `(u, s)` with `(L u)_i = s_i − Ω̂_i`:
`Σ_I µ + Σ_B ν = Σ_V s = 2πχ`.

[formula: boundary-optimality-system] -/
theorem sum_eq_of_feasible (hM : M.Valid) (hman : M.Manifold) {u s : V → ℝ}
    (hstate : ∀ i, (M.laplacian *ᵥ u) i = s i - M.hatCurvature i) :
    ∑ i ∈ M.boundaryVerticesᶜ, s i + ∑ i ∈ M.boundaryVertices, s i =
      2 * π * (M.eulerChar : ℝ) := by
  rw [Finset.sum_compl_add_sum]
  have h0 : ∑ i, (M.laplacian *ᵥ u) i = 0 :=
    Solver.sum_eq_zero_of_mulVec_eq M.laplacian_vecMul_one rfl
  have h : ∑ i, s i = ∑ i, M.hatCurvature i := by
    have : ∑ i, (s i - M.hatCurvature i) = 0 := by
      rw [← h0]; exact Finset.sum_congr rfl fun i _ => (hstate i).symm
    rw [Finset.sum_sub_distrib] at this
    linarith
  rw [h, M.sum_hatCurvature hM hman]

/-! ## The modelling step: `ν_i = k_i + (L u)_i` is the first-order boundary curvature -/

theorem hatCurvature_of_mem {i : V} (hi : i ∈ M.boundaryVertices) :
    M.hatCurvature i = M.geodesicCurvature i := by simp [hatCurvature, hi]

theorem hatCurvature_of_notMem {i : V} (hi : i ∉ M.boundaryVertices) :
    M.hatCurvature i = M.angleDefect i := by simp [hatCurvature, hi]

/-- The conformal change of metric keeps the faces, hence the boundary. -/
theorem boundaryVertices_conformalScale (u : V → ℝ) :
    (M.conformalScale u).boundaryVertices = M.boundaryVertices := rfl

/-- **First-order curvature at every vertex.** Under `ℓ̃_ij = e^{t(u_i+u_j)/2} ℓ_ij` the curvature
`Ω̂_i` of the scaled mesh (angle defect at an interior vertex, geodesic curvature at a boundary
vertex) has derivative `(L u)_i` at `t = 0`.

[formula: boundary-optimality-system] -/
theorem hasDerivAt_hatCurvature_conformalScale (hM : M.Valid) (u : V → ℝ) (i : V) :
    HasDerivAt (fun t : ℝ => (M.conformalScale (fun j => t * u j)).hatCurvature i)
      ((M.laplacian *ᵥ u) i) 0 := by
  by_cases hi : i ∈ M.boundaryVertices
  · have h := M.hasDerivAt_geodesicCurvature_conformalScale hM u i
    refine h.congr_of_eventuallyEq (Filter.Eventually.of_forall fun t => ?_)
    simp [hatCurvature, boundaryVertices_conformalScale, hi]
  · have h := M.hasDerivAt_angleDefect_conformalScale hM u i
    refine h.congr_of_eventuallyEq (Filter.Eventually.of_forall fun t => ?_)
    simp [hatCurvature, boundaryVertices_conformalScale, hi]

/-- **The boundary rows of the constraint.** Let `(u, s)` satisfy the constraint of the problem
with boundary unknowns. At a boundary vertex `i`, the geodesic curvature
`c(t) = k̃_i(t)` of the boundary of the mesh with edge lengths `e^{t(u_i+u_j)/2} ℓ_ij` satisfies
`c(0) = k_i` and `c'(0) = (L u)_i`, and `ν_i = s_i` is exactly the first-order value
`c(0) + c'(0) = k_i + (L u)_i`.

Only the first-order statement is claimed: `c` is not affine, and the theorem says nothing about
the curvature of the actually scaled boundary at `t = 1`.

[formula: boundary-optimality-system] -/
theorem boundaryProblem_firstOrder_boundary (hM : M.Valid) {u s : V → ℝ}
    (hstate : ∀ i, (M.laplacian *ᵥ u) i = s i - M.hatCurvature i) {i : V}
    (hi : i ∈ M.boundaryVertices) :
    (M.conformalScale (fun j => (0 : ℝ) * u j)).geodesicCurvature i = M.geodesicCurvature i ∧
      HasDerivAt (fun t : ℝ => (M.conformalScale (fun j => t * u j)).geodesicCurvature i)
        (s i - M.geodesicCurvature i) 0 ∧
      s i = M.geodesicCurvature i +
        deriv (fun t : ℝ => (M.conformalScale (fun j => t * u j)).geodesicCurvature i) 0 := by
  have hs : s i - M.geodesicCurvature i = (M.laplacian *ᵥ u) i := by
    have := hstate i
    rw [M.hatCurvature_of_mem hi] at this
    linarith
  have hd := M.hasDerivAt_geodesicCurvature_conformalScale hM u i
  refine ⟨?_, hs ▸ hd, ?_⟩
  · simp only [zero_mul]
    rw [M.conformalScale_zero]
  · rw [hd.deriv]
    linarith

/-- The same at an interior vertex: `µ_i = Ω_i + (L u)_i` is the first-order angle defect.

[formula: boundary-optimality-system] -/
theorem boundaryProblem_firstOrder_interior (hM : M.Valid) {u s : V → ℝ}
    (hstate : ∀ i, (M.laplacian *ᵥ u) i = s i - M.hatCurvature i) {i : V}
    (hi : i ∉ M.boundaryVertices) :
    (M.conformalScale (fun j => (0 : ℝ) * u j)).angleDefect i = M.angleDefect i ∧
      HasDerivAt (fun t : ℝ => (M.conformalScale (fun j => t * u j)).angleDefect i)
        (s i - M.angleDefect i) 0 ∧
      s i = M.angleDefect i +
        deriv (fun t : ℝ => (M.conformalScale (fun j => t * u j)).angleDefect i) 0 := by
  have hs : s i - M.angleDefect i = (M.laplacian *ᵥ u) i := by
    have := hstate i
    rw [M.hatCurvature_of_notMem hi] at this
    linarith
  have hd := M.hasDerivAt_angleDefect_conformalScale hM u i
  refine ⟨?_, hs ▸ hd, ?_⟩
  · simp only [zero_mul]
    rw [M.conformalScale_zero]
  · rw [hd.deriv]
    linarith

end TriMesh

/-! ## The hypotheses are satisfiable -/

/-- On the single equilateral triangle (a mesh with boundary, all three vertices on it) with free
boundary sites, `u = 0`, `φ = 0`, `s = k` solves the optimality system for any `wE ≥ 0`, `λ ≥ 0`;
the theorem then gives a global minimizer, and `Σ s = 2πχ`. -/
example (wE : Fin 3 → ℝ) (hwE : ∀ v, 0 ≤ wE v) (lam : ℝ) (hlam : 0 ≤ lam) :
    exampleTriangle.Valid ∧ exampleTriangle.Manifold ∧
    (exampleTriangle.boundaryProblem wE lam fun _ => Site.free).IsMinimizer
      (fun _ => 0) exampleTriangle.hatCurvature ∧
    ∑ i, exampleTriangle.hatCurvature i = 2 * π * (exampleTriangle.eulerChar : ℝ) := by
  have hV : exampleTriangle.Valid :=
    ⟨by decide, by decide, by decide, fun f => by
      simp only [TriMesh.la, TriMesh.lb, TriMesh.lc, exampleTriangle]
      exact ⟨by norm_num, by norm_num, by norm_num⟩⟩
  have hMan : exampleTriangle.Manifold := ⟨by decide, by decide⟩
  refine ⟨hV, hMan, ?_, exampleTriangle.sum_hatCurvature hV hMan⟩
  refine (exampleTriangle.boundaryProblem_isMinimizer wE lam (fun _ => Site.free) hV hwE
    (u := fun _ => 0) (φ := fun _ => 0) (s := exampleTriangle.hatCurvature) ?_ ?_ ?_).1
  · intro i; simp [show (fun _ : Fin 3 => (0 : ℝ)) = 0 from rfl]
  · intro i; simp [show (fun _ : Fin 3 => (0 : ℝ)) = 0 from rfl]
  · intro i; exact (Site.stationary_boundary_free_iff hlam _ _).mpr rfl

end DDG
