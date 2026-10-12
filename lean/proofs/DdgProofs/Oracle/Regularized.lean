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

/-! ## The Newton step in the oracle's flattened form -/

namespace Ddg.Problem

open DdgProofs DdgProofs.Oracle Matrix

variable {p : Ddg.Problem ℝ}

/-- A step `(v, q, ν)` as one vector of length `3n`, in the oracle's order `(u, φ, μ)`. -/
def flat (n : ℕ) (v q ν : ℕ → ℝ) : ℕ → ℝ :=
  fun c => if c < n then v c else if c < 2 * n then q (c - n) else ν (c - 2 * n)

private theorem sumTo_three (n : ℕ) (g : ℕ → ℝ) :
    Ddg.sumTo (3 * n) g =
      (∑ c : Fin n, g c) + (∑ c : Fin n, g (n + c)) + (∑ c : Fin n, g (2 * n + c)) := by
  rw [sumTo_eq, Fin.sum_univ_eq_sum_range (fun c => g c) (3 * n),
    show 3 * n = n + (n + n) by ring, Finset.sum_range_add, Finset.sum_range_add,
    Fin.sum_univ_eq_sum_range (fun c => g c) n, Fin.sum_univ_eq_sum_range (fun c => g (n + c)) n,
    Fin.sum_univ_eq_sum_range (fun c => g (2 * n + c)) n]
  have : ∀ x, n + (n + x) = 2 * n + x := fun x => by ring
  simp only [this]
  ring

/-- **Rows of the oracle's Jacobian `DFpiece γ σ`, applied to a step.** -/
theorem DFpiece_mulVec (p : Ddg.Problem ℝ) (γ : ℝ) (σ : ℕ → Ddg.ConeState) (v q ν : ℕ → ℝ)
    {r : ℕ} (hr : r < 3 * p.n) :
    Ddg.sumTo (3 * p.n) (fun c => p.DFpiece γ σ r c * flat p.n v q ν c) =
      if r < p.n then p.applyL v r - ν r
      else if r < 2 * p.n then
        p.wE (r - p.n) * p.mass (r - p.n) * v (r - p.n) + p.applyL q (r - p.n)
      else
        -(p.coneScale γ (r - 2 * p.n) * (σ (r - 2 * p.n)).indicator) * q (r - 2 * p.n) +
          (1 - p.coneScale γ (r - 2 * p.n) * (σ (r - 2 * p.n)).indicator) * ν (r - 2 * p.n) := by
  rw [sumTo_three]
  have ind : ∀ (i : ℕ) (hi : i < p.n) (a : ℝ) (w : ℕ → ℝ),
      ∑ c : Fin p.n, (if (c : ℕ) = i then a else 0) * w c = a * w i := by
    intro i hi a w
    have : ∀ c : Fin p.n, (if (c : ℕ) = i then a else 0) * w c =
        if c = ⟨i, hi⟩ then a * w i else 0 := by
      intro c
      by_cases h : c = ⟨i, hi⟩
      · subst h; simp
      · have h' : (c : ℕ) ≠ i := fun h'' => h (Fin.ext h'')
        simp [h, h']
    simp only [this, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  by_cases h1 : r < p.n
  · rw [ite_eq_left h1]
    have s1 : ∑ c : Fin p.n, p.DFpiece γ σ r c * flat p.n v q ν c = p.applyL v r := by
      rw [applyL, sumTo_eq]
      refine Finset.sum_congr rfl fun c _ => ?_
      simp [DFpiece, flat, h1, c.2]
    have s2 : ∑ c : Fin p.n, p.DFpiece γ σ r (p.n + c) * flat p.n v q ν (p.n + c) = 0 := by
      refine Finset.sum_eq_zero fun c _ => ?_
      have hc := c.2
      have a1 : ¬(p.n + (c : ℕ) < p.n) := by omega
      have a2 : p.n + (c : ℕ) < 2 * p.n := by omega
      simp [DFpiece, h1, a1, a2]
    have s3 : ∑ c : Fin p.n, p.DFpiece γ σ r (2 * p.n + c) * flat p.n v q ν (2 * p.n + c) =
        -ν r := by
      have e : ∀ c : Fin p.n, p.DFpiece γ σ r (2 * p.n + c) * flat p.n v q ν (2 * p.n + c) =
          (if (c : ℕ) = r then (-1 : ℝ) else 0) * ν c := by
        intro c
        have a1 : ¬(2 * p.n + (c : ℕ) < p.n) := by omega
        have a2 : ¬(2 * p.n + (c : ℕ) < 2 * p.n) := by omega
        simp [DFpiece, flat, h1, a1, a2]
      rw [Finset.sum_congr rfl fun c _ => e c, ind r h1]
      ring
    rw [s1, s2, s3]
    ring
  · rw [ite_eq_right h1]
    by_cases h2 : r < 2 * p.n
    · rw [ite_eq_left h2]
      have hi : r - p.n < p.n := by omega
      have s1 : ∑ c : Fin p.n, p.DFpiece γ σ r c * flat p.n v q ν c =
          p.wE (r - p.n) * p.mass (r - p.n) * v (r - p.n) := by
        have e : ∀ c : Fin p.n, p.DFpiece γ σ r c * flat p.n v q ν c =
            (if (c : ℕ) = r - p.n then p.wE (r - p.n) * p.mass (r - p.n) else 0) * v c := by
          intro c
          simp [DFpiece, flat, h1, h2, c.2]
        rw [Finset.sum_congr rfl fun c _ => e c, ind (r - p.n) hi]
      have s2 : ∑ c : Fin p.n, p.DFpiece γ σ r (p.n + c) * flat p.n v q ν (p.n + c) =
          p.applyL q (r - p.n) := by
        rw [applyL, sumTo_eq]
        refine Finset.sum_congr rfl fun c _ => ?_
        have hc := c.2
        have a1 : ¬(p.n + (c : ℕ) < p.n) := by omega
        have a2 : p.n + (c : ℕ) < 2 * p.n := by omega
        simp [DFpiece, flat, h1, h2, a1, a2]
      have s3 : ∑ c : Fin p.n, p.DFpiece γ σ r (2 * p.n + c) * flat p.n v q ν (2 * p.n + c) =
          0 := by
        refine Finset.sum_eq_zero fun c _ => ?_
        have a1 : ¬(2 * p.n + (c : ℕ) < p.n) := by omega
        have a2 : ¬(2 * p.n + (c : ℕ) < 2 * p.n) := by omega
        simp [DFpiece, h1, h2, a1, a2]
      rw [s1, s2, s3]
      ring
    · rw [ite_eq_right h2]
      have hi : r - 2 * p.n < p.n := by omega
      have s1 : ∑ c : Fin p.n, p.DFpiece γ σ r c * flat p.n v q ν c = 0 := by
        refine Finset.sum_eq_zero fun c _ => ?_
        simp [DFpiece, h1, h2, c.2]
      have s2 : ∑ c : Fin p.n, p.DFpiece γ σ r (p.n + c) * flat p.n v q ν (p.n + c) =
          -(p.coneScale γ (r - 2 * p.n) * (σ (r - 2 * p.n)).indicator) * q (r - 2 * p.n) := by
        have e : ∀ c : Fin p.n, p.DFpiece γ σ r (p.n + c) * flat p.n v q ν (p.n + c) =
            (if (c : ℕ) = r - 2 * p.n then
              -(p.coneScale γ (r - 2 * p.n) * (σ (r - 2 * p.n)).indicator) else 0) * q c := by
          intro c
          have hc := c.2
          have a1 : ¬(p.n + (c : ℕ) < p.n) := by omega
          have a2 : p.n + (c : ℕ) < 2 * p.n := by omega
          simp [DFpiece, flat, h1, h2, a1, a2]
        rw [Finset.sum_congr rfl fun c _ => e c, ind (r - 2 * p.n) hi]
      have s3 : ∑ c : Fin p.n, p.DFpiece γ σ r (2 * p.n + c) * flat p.n v q ν (2 * p.n + c) =
          (1 - p.coneScale γ (r - 2 * p.n) * (σ (r - 2 * p.n)).indicator) * ν (r - 2 * p.n) := by
        have e : ∀ c : Fin p.n, p.DFpiece γ σ r (2 * p.n + c) * flat p.n v q ν (2 * p.n + c) =
            (if (c : ℕ) = r - 2 * p.n then
              1 - p.coneScale γ (r - 2 * p.n) * (σ (r - 2 * p.n)).indicator else 0) * ν c := by
          intro c
          have a1 : ¬(2 * p.n + (c : ℕ) < p.n) := by omega
          have a2 : ¬(2 * p.n + (c : ℕ) < 2 * p.n) := by omega
          simp [DFpiece, flat, h1, h2, a1, a2]
        rw [Finset.sum_congr rfl fun c _ => e c, ind (r - 2 * p.n) hi]
      rw [s1, s2, s3]
      ring

/-- **`DFpiece γ σ` is the Jacobian of the affine map `Fpiece γ σ`**:
`Fpiece(x + y) = Fpiece(x) + DFpiece·y`, row by row. -/
theorem Fpiece_add (p : Ddg.Problem ℝ) (γ : ℝ) (σ : ℕ → Ddg.ConeState) (u φ μ v q ν : ℕ → ℝ)
    {r : ℕ} (hr : r < 3 * p.n) :
    p.Fpiece γ σ (fun j => u j + v j) (fun j => φ j + q j) (fun j => μ j + ν j) r =
      p.Fpiece γ σ u φ μ r +
        Ddg.sumTo (3 * p.n) (fun c => p.DFpiece γ σ r c * flat p.n v q ν c) := by
  rw [DFpiece_mulVec p γ σ v q ν hr]
  unfold Fpiece
  split_ifs with h1 h2
  · simp only [stateResidual, applyL_add]
    ring
  · simp only [adjointResidual, applyL_add]
    ring
  · simp only [pieceConeResidual]
    rw [value_affine _ _ _ _ (φ (r - 2 * p.n) + μ (r - 2 * p.n))
      (φ (r - 2 * p.n) + q (r - 2 * p.n) + (μ (r - 2 * p.n) + ν (r - 2 * p.n)))]
    ring

/-- **Stopping rule of the oracle's Newton iteration, in its own matrices.** Let the data be
admissible. Let `y = (v, q, ν)` solve the linear system `DFpiece γ σ · y = −Fpiece γ σ (u, φ, μ)`
that the oracle assembles for the pieces `σ`, and let every vertex of the new iterate lie on its
piece (`σ_i` is the state of the new point, or its canonical name). Then the new iterate is a zero
of the residual map `Fimp γ`; for `γ = 0` and nonzero masses it is a zero of `F`.
(This is the unbordered system; on a closed mesh the oracle borders it.)

[formula: newton-stopping-rule] -/
theorem Fimp_eq_zero_of_newton_step (hp : p.Admissible) {γ : ℝ} {σ : ℕ → Ddg.ConeState}
    {u φ μ v q ν : ℕ → ℝ}
    (hstep : ∀ r, r < 3 * p.n →
      Ddg.sumTo (3 * p.n) (fun c => p.DFpiece γ σ r c * flat p.n v q ν c) =
        -p.Fpiece γ σ u φ μ r)
    (hσ : ∀ i, i < p.n →
      σ i = p.impState γ (fun j => φ j + q j) (fun j => μ j + ν j) i ∨
      σ i = (p.impState γ (fun j => φ j + q j) (fun j => μ j + ν j) i).canon
        (p.t i) (p.lo i) (p.hi i)) :
    ∀ r, r < 3 * p.n →
      p.Fimp γ (fun j => u j + v j) (fun j => φ j + q j) (fun j => μ j + ν j) r = 0 := by
  refine Fimp_eq_zero_of_Fpiece_eq_zero hp (fun r hr => ?_) hσ
  rw [Fpiece_add p γ σ u φ μ v q ν hr, hstep r hr]
  ring

end Ddg.Problem
