import DdgProofs.Problem

/-!
# The optimality system characterizes the global minimizers

For the canonical discrete problem (`DdgProofs.Problem`):

* **sufficiency** (`ConeProblem.KKT.isMinimizer`): any solution `(u, φ, μ)` of the optimality
  system is a global minimizer. No assumption on `L` is needed, so the statement covers the
  Dirichlet case (`L = L_II` positive definite) and the closed case (`L` singular, constants in
  the kernel, additive constant of `u` optimized) at once; the two cases are spelled out in
  `KKT.minimizes_reduced_of_posDef` and `KKT.minimizes_closed`;
* **quadratic growth** (`KKT.objective_add_le`);
* **uniqueness** of `u` and `μ` when all `wE_i m_i > 0` (`IsMinimizer.unique`);
* **necessity** when `L` is invertible (`IsMinimizer.existsUnique_kkt`), hence an "iff" in the
  Dirichlet case (`isMinimizer_iff_exists_kkt_of_posDef`);
* the **guarantee of paper Sec. 4.5** (`IsMinimizer.distortion_le`).
-/

noncomputable section

namespace DdgProofs

open Matrix

section Abstract

variable {ι : Type*} [Fintype ι]

/-- Abstract core of sufficiency. For any matrix `A`, weights `w`, per-coordinate costs `g i`
and allowed sets `K i`: if `A u = μ − Ω`, `Aᵀφ = −w u` and `φ_i ∈ ∂(g_i + ι_{K_i})(μ_i)`, then
every feasible `(v, ν)` has objective at least that of `(u, μ)` plus `½ Σ w_i (v_i − u_i)²`. -/
theorem quadratic_growth_of_stationary (A : Matrix ι ι ℝ) (w : ι → ℝ) (g : ι → ℝ → ℝ)
    (K : ι → Set ℝ) (Ω u φ μ v ν : ι → ℝ)
    (hstate : A *ᵥ u = μ - Ω) (hadj : Aᵀ *ᵥ φ = fun i => -(w i * u i))
    (hincl : ∀ i, IsSubgradientWithinAt (g i) (K i) (φ i) (μ i))
    (hv : A *ᵥ v = ν - Ω) (hν : ∀ i, ν i ∈ K i) :
    1 / 2 * ∑ i, w i * u i ^ 2 + ∑ i, g i (μ i) + 1 / 2 * ∑ i, w i * (v i - u i) ^ 2 ≤
      1 / 2 * ∑ i, w i * v i ^ 2 + ∑ i, g i (ν i) := by
  have hd : ν - μ = A *ᵥ (v - u) := by
    rw [Matrix.mulVec_sub, hv, hstate]; ext i; simp
  have h1 : ∑ i, φ i * (ν i - μ i) = -∑ i, w i * u i * (v i - u i) := by
    have e1 : ∑ i, φ i * (ν i - μ i) = φ ⬝ᵥ (ν - μ) := by simp [dotProduct]
    rw [e1, hd, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hadj]
    simp [dotProduct, Finset.sum_neg_distrib]
  have h2 : ∑ i, w i * v i ^ 2 = ∑ i, w i * u i ^ 2 + 2 * ∑ i, w i * u i * (v i - u i) +
      ∑ i, w i * (v i - u i) ^ 2 := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun i _ => by ring
  have h3 : ∑ i, g i (μ i) + ∑ i, φ i * (ν i - μ i) ≤ ∑ i, g i (ν i) := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_le_sum fun i _ => (hincl i).2 (ν i) (hν i)
  linarith

/-- If `θ·a + θ²·C ≥ 0` for all small `θ > 0` then `a ≥ 0` (first-order optimality along a ray). -/
theorem nonneg_of_forall_small {a C : ℝ}
    (h : ∀ θ : ℝ, 0 < θ → θ ≤ 1 → 0 ≤ θ * a + θ ^ 2 * C) : 0 ≤ a := by
  by_contra ha
  push Not at ha
  have hpos : (0 : ℝ) < 2 * (|C| + 1) := by positivity
  set θ := min 1 (-a / (2 * (|C| + 1))) with hθ
  have hθ0 : 0 < θ := lt_min one_pos (div_pos (neg_pos.2 ha) hpos)
  have hθ1 : θ ≤ 1 := min_le_left _ _
  have hθ2 : θ * (2 * (|C| + 1)) ≤ -a := by
    have := min_le_right 1 (-a / (2 * (|C| + 1)))
    rwa [← hθ, le_div_iff₀ hpos] at this
  have key := h θ hθ0 hθ1
  have hC : θ * C ≤ θ * |C| := mul_le_mul_of_nonneg_left (le_abs_self C) hθ0.le
  have h3 : θ * C ≤ -a / 2 := by nlinarith [abs_nonneg C]
  have h4 : θ * a + θ ^ 2 * C = θ * (a + θ * C) := by ring
  rw [h4] at key
  have h5 : a + θ * C < 0 := by linarith
  nlinarith

/-- Abstract core of necessity (invertible `A`). If `(u, μ)` minimizes
`½ Σ w u² + Σ g_i(μ_i)` subject to `A u = μ − Ω`, `μ_i ∈ K_i`, the sets `K_i` and the costs `g_i`
are convex, and `φ` solves the adjoint equation `Aᵀφ = −w u`, then `φ_i ∈ ∂(g_i + ι_{K_i})(μ_i)`
at every `i`. -/
theorem stationary_of_minimizer [DecidableEq ι] (A : Matrix ι ι ℝ) (hA : IsUnit A.det)
    (w : ι → ℝ) (g : ι → ℝ → ℝ) (K : ι → Set ℝ) (Ω u φ μ : ι → ℝ)
    (hK : ∀ i a b θ, a ∈ K i → b ∈ K i → 0 ≤ θ → θ ≤ 1 → (1 - θ) * a + θ * b ∈ K i)
    (hg : ∀ i a b θ, a ∈ K i → b ∈ K i → 0 ≤ θ → θ ≤ 1 →
      g i ((1 - θ) * a + θ * b) ≤ (1 - θ) * g i a + θ * g i b)
    (hstate : A *ᵥ u = μ - Ω) (hμ : ∀ i, μ i ∈ K i)
    (hmin : ∀ v ν, A *ᵥ v = ν - Ω → (∀ i, ν i ∈ K i) →
      1 / 2 * ∑ i, w i * u i ^ 2 + ∑ i, g i (μ i) ≤ 1 / 2 * ∑ i, w i * v i ^ 2 + ∑ i, g i (ν i))
    (hadj : Aᵀ *ᵥ φ = fun i => -(w i * u i)) :
    ∀ i, IsSubgradientWithinAt (g i) (K i) (φ i) (μ i) := by
  intro i
  refine ⟨hμ i, fun y hy => ?_⟩
  set d := y - μ i with hd
  set z := A⁻¹ *ᵥ Pi.single i (1 : ℝ) with hz
  have hAz : A *ᵥ z = Pi.single i (1 : ℝ) := by
    rw [hz, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv A hA, Matrix.one_mulVec]
  -- the first-order term
  have hwuz : ∑ j, w j * u j * z j = -φ i := by
    have e1 : ∑ j, w j * u j * z j = -((Aᵀ *ᵥ φ) ⬝ᵥ z) := by
      rw [hadj]; simp [dotProduct, Finset.sum_neg_distrib]
    rw [e1, Matrix.mulVec_transpose, ← Matrix.dotProduct_mulVec, hAz]
    simp
  have key : ∀ θ : ℝ, 0 < θ → θ ≤ 1 →
      0 ≤ θ * (g i y - g i (μ i) - φ i * d) + θ ^ 2 * (1 / 2 * d ^ 2 * ∑ j, w j * z j ^ 2) := by
    intro θ hθ0 hθ1
    set ν := Function.update μ i ((1 - θ) * μ i + θ * y) with hν
    set v := u + (θ * d) • z with hv
    have hνi : ν i = (1 - θ) * μ i + θ * y := by simp [hν]
    have hνK : ∀ j, ν j ∈ K j := by
      intro j
      by_cases hj : j = i
      · subst hj; rw [hνi]; exact hK j _ _ θ (hμ j) hy hθ0.le hθ1
      · rw [hν, Function.update_of_ne hj]; exact hμ j
    have hfeas : A *ᵥ v = ν - Ω := by
      rw [hv, Matrix.mulVec_add, Matrix.mulVec_smul, hAz, hstate]
      ext j
      by_cases hj : j = i
      · subst hj; simp [hν, hd]; ring
      · simp [hν, Function.update_of_ne hj, Pi.single_eq_of_ne hj]
    have h1 := hmin v ν hfeas hνK
    have hgsum : ∑ j, g j (ν j) = ∑ j, g j (μ j) + (g i (ν i) - g i (μ i)) := by
      have : ∑ j, (g j (ν j) - g j (μ j)) = g i (ν i) - g i (μ i) := by
        rw [Finset.sum_eq_single i]
        · intro j _ hj; rw [hν, Function.update_of_ne hj]; ring
        · intro h; exact absurd (Finset.mem_univ i) h
      rw [← this, Finset.sum_sub_distrib]; ring
    have hconv : g i (ν i) ≤ (1 - θ) * g i (μ i) + θ * g i y := by
      rw [hνi]; exact hg i _ _ θ (hμ i) hy hθ0.le hθ1
    have hquad : ∑ j, w j * v j ^ 2 = ∑ j, w j * u j ^ 2 + 2 * (θ * d) * ∑ j, w j * u j * z j +
        (θ * d) ^ 2 * ∑ j, w j * z j ^ 2 := by
      rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun j _ => ?_
      simp only [hv, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
      ring
    rw [hgsum, hquad, hwuz] at h1
    nlinarith
  have := nonneg_of_forall_small key
  rw [hd] at this
  linarith

end Abstract

namespace ConeProblem

variable {ι : Type*} [Fintype ι] {P : ConeProblem ι} {u φ μ v ν : ι → ℝ}

/-! ## Sufficiency -/

/-- **Quadratic growth at a solution of the optimality system.** If `(u, φ, μ)` solves the
optimality system then every feasible `(v, ν)` satisfies
`J(u, μ) + ½ Σ_i wE_i m_i (v_i − u_i)² ≤ J(v, ν)`.
No assumption on `L`, the weights or the vertex data is needed.

[formula: optimality-system-sufficient] -/
theorem KKT.objective_add_le (h : P.KKT u φ μ) (hv : P.Feasible v ν) :
    P.objective u μ + 1 / 2 * ∑ i, P.wE i * P.m i * (v i - u i) ^ 2 ≤ P.objective v ν := by
  have := quadratic_growth_of_stationary P.L (fun i => P.wE i * P.m i)
    (fun i x => P.lam * (P.site i).cost x) (fun i => (P.site i).feasible) P.Ω u φ μ v ν
    h.state h.adjoint h.inclusion hv.1 hv.2
  simp only [objective, distortion, coneCost, Finset.mul_sum] at this ⊢
  linarith

/-- **Sufficiency of the optimality system.** For the canonical discrete problem with
`wE ≥ 0`, `m ≥ 0`: if `(u, φ, μ)` satisfies the state equation `L u = μ − Ω`, the adjoint equation
`Lᵀφ = −wE m u` and the per-vertex inclusions, then `(u, μ)` is a global minimizer.

The statement holds for every matrix `L`; in particular for the Dirichlet case (`L = L_II`
positive definite) and for the closed case (`L` positive semidefinite with the constants in its
kernel, where the competitors include every additive constant of `u`). Vertex weights, angle
bounds, free sites, forbidden vertices and prescribed cones are all covered by `Site`.

[formula: optimality-system-sufficient] -/
theorem KKT.isMinimizer (hwE : ∀ i, 0 ≤ P.wE i) (hm : ∀ i, 0 ≤ P.m i) (h : P.KKT u φ μ) :
    P.IsMinimizer u μ := by
  refine ⟨h.feasible, fun v ν hv => ?_⟩
  have h1 := h.objective_add_le hv
  have h2 : 0 ≤ ∑ i, P.wE i * P.m i * (v i - u i) ^ 2 :=
    Finset.sum_nonneg fun i _ => mul_nonneg (mul_nonneg (hwE i) (hm i)) (sq_nonneg _)
  linarith

/-- With a positive definite `L` (Dirichlet case) the state `u` is determined by `μ`. -/
theorem feasible_iff_of_posDef [DecidableEq ι] (hL : P.L.PosDef) :
    P.Feasible u μ ↔ u = P.L⁻¹ *ᵥ (μ - P.Ω) ∧ ∀ i, μ i ∈ (P.site i).feasible := by
  have hdet : IsUnit P.L.det := (Matrix.isUnit_iff_isUnit_det _).mp hL.isUnit
  unfold Feasible
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨?_, h2⟩
    rw [← h1, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet, Matrix.one_mulVec]
  · rintro ⟨h1, h2⟩
    refine ⟨?_, h2⟩
    rw [h1, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hdet, Matrix.one_mulVec]

/-- **Sufficiency, Dirichlet case, in reduced form.** Let the interior block `L` be positive
definite, so that `u = L⁻¹(μ − Ω)` is a function of `μ`. If `(u, φ, μ)` solves the optimality
system then `μ` minimizes the reduced objective `ν ↦ J(L⁻¹(ν − Ω), ν)` over all admissible `ν`.

[formula: optimality-system-sufficient] -/
theorem KKT.minimizes_reduced_of_posDef [DecidableEq ι] (hL : P.L.PosDef)
    (hwE : ∀ i, 0 ≤ P.wE i) (hm : ∀ i, 0 ≤ P.m i) (h : P.KKT u φ μ) :
    u = P.L⁻¹ *ᵥ (μ - P.Ω) ∧ (∀ i, μ i ∈ (P.site i).feasible) ∧
      ∀ ν : ι → ℝ, (∀ i, ν i ∈ (P.site i).feasible) →
        P.objective (P.L⁻¹ *ᵥ (μ - P.Ω)) μ ≤ P.objective (P.L⁻¹ *ᵥ (ν - P.Ω)) ν := by
  have hmin := h.isMinimizer hwE hm
  have hu := ((feasible_iff_of_posDef hL).mp hmin.1).1
  refine ⟨hu, (h.feasible).2, fun ν hν => ?_⟩
  rw [← hu]
  exact hmin.2 _ ν ((feasible_iff_of_posDef hL).mpr ⟨rfl, hν⟩)

/-- On a closed mesh (`L·1 = 0`) adding a constant to `u` keeps the constraints. -/
theorem Feasible.add_const (hL1 : P.L *ᵥ (fun _ => (1 : ℝ)) = 0) (h : P.Feasible u μ) (c : ℝ) :
    P.Feasible (fun i => u i + c) μ := by
  refine ⟨?_, h.2⟩
  have : (fun i => u i + c) = u + c • (fun _ => (1 : ℝ)) := by
    ext i; simp
  rw [this, Matrix.mulVec_add, Matrix.mulVec_smul, hL1, h.1]
  simp

/-- **Sufficiency, closed case.** Let `L·1 = 0` (closed mesh), so that `u` is determined by `μ`
only up to an additive constant. If `(u, φ, μ)` solves the optimality system then `(u, μ)` is at
least as good as `(v + c, ν)` for every admissible `ν`, every solution `v` of `L v = ν − Ω` and
every constant `c`: the cone angles *and* the additive constant are optimal.

[formula: optimality-system-sufficient] -/
theorem KKT.minimizes_closed (hL1 : P.L *ᵥ (fun _ => (1 : ℝ)) = 0)
    (hwE : ∀ i, 0 ≤ P.wE i) (hm : ∀ i, 0 ≤ P.m i) (h : P.KKT u φ μ)
    (hv : P.L *ᵥ v = ν - P.Ω) (hν : ∀ i, ν i ∈ (P.site i).feasible) (c : ℝ) :
    P.objective u μ ≤ P.objective (fun i => v i + c) ν :=
  (h.isMinimizer hwE hm).2 _ ν (Feasible.add_const hL1 ⟨hv, hν⟩ c)

/-! ## The guarantee of paper Sec. 4.5 -/

/-- **Guarantee of paper Sec. 4.5, discrete form.** Let `(u, μ)` be a minimizer and `(v, ν)` any
admissible competitor. Then `E(u) ≤ E(v) + λ(R(ν) − R(μ))`, where `E` is the distortion and `R`
the weighted cone norm `Σ wR_i |·|`.

[formula: least-distortion-guarantee] -/
theorem IsMinimizer.distortion_le_add (h : P.IsMinimizer u μ) (hv : P.Feasible v ν) :
    P.distortion u ≤ P.distortion v + P.lam * (P.coneCost ν - P.coneCost μ) := by
  have := h.2 v ν hv
  unfold objective at this
  linarith

/-- **Guarantee of paper Sec. 4.5, discrete form.** A minimizer `(u, μ)` has the least distortion
among all admissible cone configurations `ν` whose weighted norm `Σ wR_i |ν_i|` is no larger than
that of `μ`.

[formula: least-distortion-guarantee] -/
theorem IsMinimizer.distortion_le (hlam : 0 ≤ P.lam) (h : P.IsMinimizer u μ)
    (hv : P.Feasible v ν) (hnorm : P.coneCost ν ≤ P.coneCost μ) :
    P.distortion u ≤ P.distortion v := by
  have h1 := h.distortion_le_add hv
  have h2 : P.lam * (P.coneCost ν - P.coneCost μ) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos hlam (by linarith)
  linarith

/-! ## Uniqueness -/

/-- The midpoint of two feasible pairs is feasible, and the objective is midpoint convex with
the explicit defect `⅛ Σ wE m (u − v)²`. -/
private theorem midpoint_objective_le (hP : P.WellFormed) (h₁ : P.Feasible u μ)
    (h₂ : P.Feasible v ν) :
    P.Feasible (fun i => (1 - 1 / 2) * u i + 1 / 2 * v i)
        (fun i => (1 - 1 / 2) * μ i + 1 / 2 * ν i) ∧
      P.objective (fun i => (1 - 1 / 2) * u i + 1 / 2 * v i)
          (fun i => (1 - 1 / 2) * μ i + 1 / 2 * ν i) +
        1 / 8 * ∑ i, P.wE i * P.m i * (u i - v i) ^ 2 ≤
        1 / 2 * P.objective u μ + 1 / 2 * P.objective v ν := by
  constructor
  · refine ⟨?_, fun i => (P.site i).feasible_convex (h₁.2 i) (h₂.2 i) (by norm_num) (by norm_num)⟩
    have : (fun i => (1 - 1 / 2) * u i + 1 / 2 * v i) = ((1 - 1 / 2 : ℝ)) • u + (1 / 2 : ℝ) • v := by
      ext i; simp
    rw [this, Matrix.mulVec_add, Matrix.mulVec_smul, Matrix.mulVec_smul, h₁.1, h₂.1]
    ext i; simp; ring
  · have hc : P.coneCost (fun i => (1 - 1 / 2) * μ i + 1 / 2 * ν i) ≤
        1 / 2 * P.coneCost μ + 1 / 2 * P.coneCost ν := by
      unfold coneCost
      rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
      refine Finset.sum_le_sum fun i _ => ?_
      have := (P.site i).cost_convex (hP.site i) (μ i) (ν i) (θ := 1 / 2) (by norm_num)
        (by norm_num)
      linarith
    have hd : P.distortion (fun i => (1 - 1 / 2) * u i + 1 / 2 * v i) +
        1 / 8 * ∑ i, P.wE i * P.m i * (u i - v i) ^ 2 =
        1 / 2 * P.distortion u + 1 / 2 * P.distortion v := by
      unfold distortion
      simp only [Finset.mul_sum]
      rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun i _ => by ring
    have := mul_le_mul_of_nonneg_left hc hP.lam_nonneg
    unfold objective
    linarith

/-- Two minimizers have the same distortion at every vertex where distortion is measured:
`wE_i m_i (u_i − v_i)² = 0`. -/
theorem IsMinimizer.weighted_sq_sub_eq_zero (hP : P.WellFormed) (h₁ : P.IsMinimizer u μ)
    (h₂ : P.IsMinimizer v ν) (i : ι) : P.wE i * P.m i * (u i - v i) ^ 2 = 0 := by
  obtain ⟨hf, hle⟩ := midpoint_objective_le hP h₁.1 h₂.1
  have e1 := h₁.2 _ _ hf
  have e2 := h₂.2 _ _ hf
  have hnn : ∀ j ∈ Finset.univ, 0 ≤ P.wE j * P.m j * (u j - v j) ^ 2 := fun j _ =>
    mul_nonneg (mul_nonneg (hP.wE_nonneg j) (hP.m_nonneg j)) (sq_nonneg _)
  have hsum : ∑ j, P.wE j * P.m j * (u j - v j) ^ 2 = 0 :=
    le_antisymm (by linarith) (Finset.sum_nonneg hnn)
  exact (Finset.sum_eq_zero_iff_of_nonneg hnn).mp hsum i (Finset.mem_univ i)

/-- **Uniqueness.** If every vertex has `wE_i m_i > 0`, the minimizer is unique: any two
minimizers have the same `u` and the same `μ`. (Holds for every `L`, so also on closed meshes,
where it says in particular that the optimal additive constant is unique.)

[formula: minimizer-unique] -/
theorem IsMinimizer.unique (hP : P.WellFormed) (hpos : ∀ i, 0 < P.wE i * P.m i)
    (h₁ : P.IsMinimizer u μ) (h₂ : P.IsMinimizer v ν) : u = v ∧ μ = ν := by
  have huv : u = v := by
    ext i
    have h := h₁.weighted_sq_sub_eq_zero hP h₂ i
    have h' : (u i - v i) ^ 2 = 0 := by
      rcases mul_eq_zero.mp h with h0 | h0
      · exact absurd h0 (hpos i).ne'
      · exact h0
    have := pow_eq_zero_iff (two_ne_zero) |>.mp h'
    linarith
  refine ⟨huv, ?_⟩
  have e1 := h₁.1.1
  have e2 := h₂.1.1
  rw [huv, e2] at e1
  ext i
  have := congrFun e1 i
  simp only [Pi.sub_apply] at this
  linarith

/-! ## Necessity (invertible `L`) -/

/-- **Necessity of the optimality system when `L` is invertible** (Dirichlet case). If `(u, μ)`
is a global minimizer of a well-formed problem and `det L ≠ 0`, there is exactly one adjoint
state `φ` such that `(u, φ, μ)` solves the optimality system.

[formula: optimality-system-necessary] -/
theorem IsMinimizer.existsUnique_kkt [DecidableEq ι] (hP : P.WellFormed) (hL : IsUnit P.L.det)
    (h : P.IsMinimizer u μ) : ∃! φ, P.KKT u φ μ := by
  have hLt : IsUnit P.Lᵀ.det := by rw [Matrix.det_transpose]; exact hL
  set φ := P.Lᵀ⁻¹ *ᵥ (fun i => -(P.wE i * P.m i * u i)) with hφ
  have hadj : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i) := by
    rw [hφ, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hLt, Matrix.one_mulVec]
  refine ⟨φ, ⟨h.1.1, hadj, ?_⟩, fun ψ hψ => ?_⟩
  · refine stationary_of_minimizer P.L hL (fun i => P.wE i * P.m i)
      (fun i x => P.lam * (P.site i).cost x) (fun i => (P.site i).feasible) P.Ω u φ μ
      (fun i a b θ ha hb h0 h1 => (P.site i).feasible_convex ha hb h0 h1)
      (fun i a b θ _ _ h0 h1 => ?_) h.1.1 h.1.2 (fun v ν hv hν => ?_) hadj
    · have := (P.site i).cost_convex (hP.site i) a b h0 h1
      have := mul_le_mul_of_nonneg_left this hP.lam_nonneg
      linarith
    · have := h.2 v ν ⟨hv, hν⟩
      simp only [objective, distortion, coneCost, Finset.mul_sum] at this ⊢
      linarith
  · have : P.Lᵀ *ᵥ ψ = P.Lᵀ *ᵥ φ := by rw [hψ.adjoint, hadj]
    have h2 := congrArg (fun x => P.Lᵀ⁻¹ *ᵥ x) this
    simp only [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hLt, Matrix.one_mulVec] at h2
    exact h2

/-- **Dirichlet case: the optimality system characterizes the minimizers.** For a well-formed
problem whose matrix `L` (the interior block of the cotan Laplacian) is positive definite,
`(u, μ)` is a global minimizer if and only if some adjoint state `φ` completes it to a solution
of the optimality system.

[formula: optimality-system-necessary] -/
theorem isMinimizer_iff_exists_kkt_of_posDef [DecidableEq ι] (hP : P.WellFormed)
    (hL : P.L.PosDef) : P.IsMinimizer u μ ↔ ∃ φ, P.KKT u φ μ := by
  constructor
  · intro h
    exact (h.existsUnique_kkt hP ((Matrix.isUnit_iff_isUnit_det _).mp hL.isUnit)).exists
  · rintro ⟨φ, h⟩
    exact h.isMinimizer hP.wE_nonneg hP.m_nonneg

end ConeProblem

/-! ## The hypotheses are satisfiable -/

/-- A one-vertex Dirichlet instance (`L = [1]`, `Ω = 3`, `m = wE = λ = 1`, default vertex) with
a solution of the optimality system: `u = −1`, `φ = 1`, `μ = 2`. -/
example : ∃ (P : ConeProblem (Fin 1)) (u φ μ : Fin 1 → ℝ),
    P.WellFormed ∧ P.L.PosDef ∧ P.KKT u φ μ := by
  refine ⟨⟨1, fun _ => 3, fun _ => 1, fun _ => 1, 1, fun _ => Site.plain⟩,
    fun _ => -1, fun _ => 1, fun _ => 2, ?_, Matrix.PosDef.one, ?_, ?_, ?_⟩
  · exact ⟨fun _ => zero_le_one, fun _ => zero_le_one, zero_le_one,
      fun _ => ⟨zero_le_one, bot_le, le_top⟩⟩
  · ext i; simp; norm_num
  · ext i; simp
  · intro i
    rw [Site.stationary_plain_iff zero_le_one]
    rw [softThreshold_of_gt zero_le_one (by norm_num)]; norm_num

/-- A two-vertex "closed" instance (`L = [[1, −1], [−1, 1]]`, `L·1 = 0`, `Ω = (1, 1)`) with a
solution of the optimality system: `u = 0`, `φ = (1, 1)`, `μ = (1, 1)`; note `Σ μ = Σ Ω`. -/
example : ∃ (P : ConeProblem (Fin 2)) (u φ μ : Fin 2 → ℝ),
    P.WellFormed ∧ P.L *ᵥ (fun _ => (1 : ℝ)) = 0 ∧ P.KKT u φ μ := by
  refine ⟨⟨!![1, -1; -1, 1], fun _ => 1, fun _ => 1, fun _ => 1, 1, fun _ => Site.plain⟩,
    fun _ => 0, fun _ => 1, fun _ => 1, ?_, ?_, ?_, ?_, ?_⟩
  · exact ⟨fun _ => zero_le_one, fun _ => zero_le_one, zero_le_one,
      fun _ => ⟨zero_le_one, bot_le, le_top⟩⟩
  · ext i; fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  · ext i; fin_cases i <;> simp [Matrix.mulVec, dotProduct]
  · ext i; fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  · intro i
    rw [Site.stationary_plain_iff zero_le_one]
    rw [softThreshold_of_gt zero_le_one (by norm_num)]; norm_num

end DdgProofs
