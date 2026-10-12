import DdgProofs.Oracle.Problem

/-!
# The oracle's regularized system with `μ` kept as an unknown, and its piecewise linear systems

The oracle solves the regularized problems through the residual map `Fimp γ (u, φ, μ)` whose cone
row is `μ_i − clamp(c_i·P_{t_i}(φ_i + μ_i))`, `c_i = m_i/(m_i + γ)` (not in the paper, which
eliminates `μ`). Each Newton step solves the affine system `Fpiece γ σ = 0` of a choice of pieces
`σ`. This file proves, about those definitions:

* `Fimp_zero_eq_F`: `Fimp 0 = F` (the exact system is the member `γ = 0`);
* `regObjective_le_of_Fimp_eq_zero`: a zero of `Fimp γ` minimizes the regularized objective
  `J(u, μ) + Σ (γ/(2 m_i)) μ_i²`;
* `Fimp_eq_zero_of_Fpiece_eq_zero`: if a point solves the linear system of the pieces `σ` and
  lies on those pieces (`σ` equal to the state of the point, or to its canonical name), then it
  is a zero of `Fimp γ` — the stopping rule "active sets unchanged ⟹ exact solution" for the
  oracle's (possibly safeguarded) active-set iteration.
-/

noncomputable section

namespace DdgProofs

namespace Oracle

/-- When `lo < hi`, a piece and its canonical name have the same value everywhere. -/
theorem value_canon_of_lt {t lo hi : ℝ} (ht : 0 ≤ t) (h : lo < hi) (k x : ℝ)
    (σ : Ddg.ConeState) : (σ.canon t lo hi).value k t lo hi x = σ.value k t lo hi x := by
  unfold Ddg.ConeState.canon
  rw [ite_eq_left h]
  cases σ with
  | zero => rfl
  | pos => rfl
  | neg =>
    dsimp only
    simp only [zero_eq]
    split_ifs with h2
    · rfl
    · have ht0 : t = 0 := le_antisymm (not_lt.mp h2) ht
      simp [Ddg.ConeState.value, ht0]
  | lower =>
    dsimp only
    simp only [zero_eq]
    split_ifs with h2
    · rfl
    · push Not at h2
      have h0 : lo = 0 := le_antisymm h2.1 h2.2
      simp [Ddg.ConeState.value, h0]
  | upper =>
    dsimp only
    simp only [zero_eq]
    split_ifs with h2
    · rfl
    · push Not at h2
      have h0 : hi = 0 := le_antisymm h2.1 h2.2
      simp [Ddg.ConeState.value, h0]

/-- The canonical name of the piece a point lies on has the same value at that point. -/
theorem value_canon_coneState {t lo hi : ℝ} (ht : 0 ≤ t) (hlh : lo ≤ hi) (k x : ℝ) :
    ((Ddg.coneState k t lo hi x).canon t lo hi).value k t lo hi x =
      (Ddg.coneState k t lo hi x).value k t lo hi x := by
  by_cases h1 : lo < hi
  · exact value_canon_of_lt ht h1 k x _
  · have heq : lo = hi := le_antisymm hlh (not_lt.mp h1)
    unfold Ddg.ConeState.canon
    rw [ite_eq_right h1]
    unfold Ddg.coneState
    by_cases h2 : hi ≤ k * Ddg.softThreshold t x
    · rw [ite_eq_left h2]
    · have h3 : k * Ddg.softThreshold t x ≤ lo := by rw [heq]; exact le_of_not_ge h2
      rw [ite_eq_right h2, ite_eq_left h3]
      simp [Ddg.ConeState.value, heq]

/-- A vertex whose cone angle is the regularized fixed point is stationary for the regularized
problem. -/
theorem stationary_siteOf_of_eq_regularizedMeasure {t lo hi m γ : ℝ}
    (h : AdmissibleVertex t lo hi) (hm : 0 < m) (hγ : 0 ≤ γ) {μ φ : ℝ}
    (hμ : μ = Ddg.regularizedMeasure (m / (m + γ)) t lo hi (φ + μ)) :
    (siteOf t lo hi).Stationary 1 (φ - γ / m * μ) μ := by
  unfold siteOf
  split_ifs with h1
  · have h2 := (eq_regularizedMeasure_add_iff h.1 hm hγ h1.1 h1.2 μ φ).mp hμ
    rw [inSubdiff_iff h.1 h1.1 h1.2] at h2
    unfold Site.Stationary
    have hf : (fun x => (1 : ℝ) * (Site.penalized t (lo : EReal) (hi : EReal)).cost x) =
        fun x => t * |x| := by
      funext x; simp [Site.cost]
    rw [hf]
    exact h2
  · rcases h.2 with h2 | h2
    · exact absurd h2 h1
    · rw [Site.stationary_fixed_iff]
      rw [hμ]
      unfold Ddg.regularizedMeasure Ddg.clamp
      rw [← h2.2]
      exact min_eq_left (le_max_left _ _)

end Oracle

end DdgProofs

namespace Ddg.Problem

open DdgProofs DdgProofs.Oracle Matrix

variable {p : Ddg.Problem ℝ}

theorem Fimp_eq_zero_iff (p : Ddg.Problem ℝ) (γ : ℝ) (u φ μ : ℕ → ℝ) :
    (∀ r, r < 3 * p.n → p.Fimp γ u φ μ r = 0) ↔
      ∀ i, i < p.n → p.stateResidual u μ i = 0 ∧ p.adjointResidual u φ i = 0 ∧
        p.regConeResidual γ φ μ i = 0 := by
  constructor
  · intro h i hi
    refine ⟨?_, ?_, ?_⟩
    · have := h i (by omega)
      simpa [Fimp, hi] using this
    · have := h (p.n + i) (by omega)
      have h1 : ¬(p.n + i < p.n) := by omega
      have h2 : p.n + i < 2 * p.n := by omega
      simpa [Fimp, h1, h2] using this
    · have := h (2 * p.n + i) (by omega)
      have h1 : ¬(2 * p.n + i < p.n) := by omega
      have h2 : ¬(2 * p.n + i < 2 * p.n) := by omega
      simpa [Fimp, h1, h2] using this
  · intro h r hr
    unfold Fimp
    split_ifs with h1 h2
    · exact (h r h1).1
    · exact (h (r - p.n) (by omega)).2.1
    · exact (h (r - 2 * p.n) (by omega)).2.2

theorem Fpiece_eq_zero_iff (p : Ddg.Problem ℝ) (γ : ℝ) (σ : ℕ → Ddg.ConeState)
    (u φ μ : ℕ → ℝ) :
    (∀ r, r < 3 * p.n → p.Fpiece γ σ u φ μ r = 0) ↔
      ∀ i, i < p.n → p.stateResidual u μ i = 0 ∧ p.adjointResidual u φ i = 0 ∧
        p.pieceConeResidual γ σ φ μ i = 0 := by
  constructor
  · intro h i hi
    refine ⟨?_, ?_, ?_⟩
    · have := h i (by omega)
      simpa [Fpiece, hi] using this
    · have := h (p.n + i) (by omega)
      have h1 : ¬(p.n + i < p.n) := by omega
      have h2 : p.n + i < 2 * p.n := by omega
      simpa [Fpiece, h1, h2] using this
    · have := h (2 * p.n + i) (by omega)
      have h1 : ¬(2 * p.n + i < p.n) := by omega
      have h2 : ¬(2 * p.n + i < 2 * p.n) := by omega
      simpa [Fpiece, h1, h2] using this
  · intro h r hr
    unfold Fpiece
    split_ifs with h1 h2
    · exact (h r h1).1
    · exact (h (r - p.n) (by omega)).2.1
    · exact (h (r - 2 * p.n) (by omega)).2.2

/-- **The exact system is the member `γ = 0` of the oracle's regularized family**: `Fimp 0 = F`
(for nonzero masses). -/
theorem Fimp_zero_eq_F (hm : ∀ i, i < p.n → p.mass i ≠ 0) (u φ μ : ℕ → ℝ) {r : ℕ}
    (hr : r < 3 * p.n) : p.Fimp 0 u φ μ r = p.F u φ μ r := by
  unfold Fimp F
  split_ifs with h1 h2
  · rfl
  · rfl
  · have hi : r - 2 * p.n < p.n := by omega
    simp only [regConeResidual, coneResidual, coneScale, add_zero, div_self (hm _ hi),
      boundedThreshold_eq_regularizedMeasure]

/-- **Stopping rule for the oracle's piecewise linear systems.** Let the data be admissible. If
`(u, φ, μ)` solves the linear system of the pieces `σ` (`Fpiece γ σ = 0`) and every vertex lies on
its piece — `σ_i` is the state `impState γ φ μ i` of the point, or its canonical name — then
`(u, φ, μ)` is a zero of the residual map `Fimp γ`. Since a Newton step with the Jacobian
`DFpiece γ σ` lands on the zero of the affine map `Fpiece γ σ`, this is "active sets unchanged ⟹
exact solution" for the oracle's iteration, including its safeguarded choice of `σ`.

[formula: newton-stopping-rule] -/
theorem Fimp_eq_zero_of_Fpiece_eq_zero (hp : p.Admissible) {γ : ℝ} {σ : ℕ → Ddg.ConeState}
    {u φ μ : ℕ → ℝ} (hF : ∀ r, r < 3 * p.n → p.Fpiece γ σ u φ μ r = 0)
    (hσ : ∀ i, i < p.n → σ i = p.impState γ φ μ i ∨
      σ i = (p.impState γ φ μ i).canon (p.t i) (p.lo i) (p.hi i)) :
    ∀ r, r < 3 * p.n → p.Fimp γ u φ μ r = 0 := by
  rw [Fpiece_eq_zero_iff] at hF
  rw [Fimp_eq_zero_iff]
  intro i hi
  obtain ⟨h1, h2, h3⟩ := hF i hi
  refine ⟨h1, h2, ?_⟩
  have hv := hp.vertex i hi
  simp only [pieceConeResidual] at h3
  simp only [regConeResidual, regularizedMeasure_eq_value hv.1 hv.le]
  rcases hσ i hi with h | h
  · rw [h] at h3; exact h3
  · rw [h] at h3
    simp only [impState] at h3
    rw [value_canon_coneState hv.1 hv.le] at h3
    exact h3

/-- The oracle's regularized objective `J(u, μ) + Σ_i (γ/(2 m_i)) μ_i²`: the discretization of
paper Eq. 12 with the mass factor (the term is `(γ/2)‖μ‖²_{L²}` for the density `μ_i/m_i`). -/
def regObjective (p : Ddg.Problem ℝ) (γ : ℝ) (u μ : ℕ → ℝ) : ℝ :=
  p.objective u μ + ∑ i : Fin p.n, γ / p.mass i / 2 * μ i ^ 2

/-- **A zero of the oracle's regularized residual map solves the regularized problem.** For
admissible data with positive masses and `γ ≥ 0`: if `Fimp γ (u, φ, μ) = 0` then `(u, μ)` is
feasible and minimizes `J(u, μ) + Σ (γ/(2 m_i)) μ_i²` over all feasible pairs.

[formula: regularized-optimality-sufficient] -/
theorem regObjective_le_of_Fimp_eq_zero (hp : p.Admissible)
    (hm : ∀ i, i < p.n → 0 < p.mass i) {γ : ℝ} (hγ : 0 ≤ γ) {u φ μ v ν : ℕ → ℝ}
    (hF : ∀ r, r < 3 * p.n → p.Fimp γ u φ μ r = 0) (hv : p.Feasible v ν) :
    p.Feasible u μ ∧ p.regObjective γ u μ ≤ p.regObjective γ v ν := by
  rw [Fimp_eq_zero_iff] at hF
  have hP := toConeProblem_wellFormed hp
  have hstate : p.toConeProblem.L *ᵥ res p.n u = res p.n μ - p.toConeProblem.Ω := by
    ext i
    have := (hF i i.2).1
    rw [stateResidual, applyL_eq p u i] at this
    simp only [Pi.sub_apply, toConeProblem, res] at this ⊢
    linarith
  have hadj : p.toConeProblem.Lᵀ *ᵥ res p.n φ =
      fun i => -(p.toConeProblem.wE i * p.toConeProblem.m i * res p.n u i) := by
    have hsymm : p.toConeProblem.Lᵀ = p.toConeProblem.L := by
      ext i j
      simp only [toConeProblem, Matrix.transpose_apply, Matrix.of_apply]
      exact hp.L_symm j i j.2 i.2
    rw [hsymm]
    ext i
    have := (hF i i.2).2.1
    rw [adjointResidual, applyL_eq p φ i] at this
    simp only [toConeProblem, res] at this ⊢
    linarith
  have hst : ∀ i : Fin p.n, (p.toConeProblem.site i).Stationary p.toConeProblem.lam
      (res p.n φ i - γ / p.mass i * res p.n μ i) (res p.n μ i) := by
    intro i
    have h3 := (hF i i.2).2.2
    simp only [regConeResidual, coneScale] at h3
    exact stationary_siteOf_of_eq_regularizedMeasure (hp.vertex i i.2) (hm i i.2) hγ
      (by simp only [res]; linarith)
  have hmin := ConeProblem.isRegMinimizer_of_stationary hP
    (γ := fun i : Fin p.n => γ / p.mass i) (fun i => div_nonneg hγ (hm i i.2).le) hstate hadj hst
  refine ⟨(feasible_iff hp u μ).mpr hmin.1, ?_⟩
  have := hmin.2 _ _ ((feasible_iff hp v ν).mp hv)
  simp only [ConeProblem.regObjective, ← objective_eq hp] at this
  exact this

end Ddg.Problem
