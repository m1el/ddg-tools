import Mathlib.Analysis.Complex.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-!
# Intrinsic simplification: the memoryless cost is a 1-Wasserstein distance (Eq. 5, 6)

`crates/ddg-intrinsic-simplification/SPEC.md` §4.3, §4.4, §13 (`iem-memoryless-cost`); the paper
Sec. 5.1, Eq. 5, 6, Fig. 9 ("This cost also coincides with the so-called 1-Wasserstein distance
between the old and new mass distribution", citing Peyré and Cuturi 2019, Ch. 2).

Mathlib has no optimal transport; it is modelled directly. Points form a finite type `X` with a
cost `d` that is zero on the diagonal and satisfies the triangle inequality (a pseudometric; on
the flattened star `d i j = ℓ_ij`). A coupling of two mass functions `μ ν : X → ℝ` is a
non-negative matrix `π` with row sums `μ` and column sums `ν` (`IsCoupling`); its cost is
`Σ_x Σ_y π x y · d x y`. The 1-Wasserstein distance is the least cost over all couplings.

Removing `i` sends the fraction `α_y` of `m_i` to `y` (Eq. 5): the old distribution is `m`, the new
one `m + m_i α − m_i δ_i` (`newMass`).

* `memorylessCost_eq` — with `t = 0` everywhere, the first form of Eq. 6, `Σ_j m̃_j ‖t̃_j‖`, equals
  the second, `Σ_j α_ij m_i ℓ_ij` (`‖e_ji‖ = ℓ_ij`).
* `wasserstein_eq` — that number is the least cost of a coupling of the old and the new
  distribution (`IsLeast`): the plan that keeps every `m_j` in place and moves `α_j m_i` from `i`
  to `j` attains it, and no coupling does better (the 1-Lipschitz function `d(i, ·)` bounds every
  coupling from below: Kantorovich duality, the easy direction).
* `coupling_single_source_unique` — moving the mass `m_i` alone (`m_i δ_i` to `m_i α`) has exactly
  one coupling, the forced one, of cost `m_i Σ α_j d(i, j)` (the brief's "single source").

[formula: iem-memoryless-cost]
-/

noncomputable section

namespace DDG

namespace IntrinsicSimplification

open Finset

/-! ## The two forms of Eq. 6 with `t = 0` -/

/-- Eq. 5: `m̃_j = m_j + α_j m_i`. -/
def massUpdate (mj a mi : ℝ) : ℝ := mj + a * mi

/-- Eq. 8 with the convention of SPEC §4.3 (`t̃_j = t_j` unless `m̃_j > 0`):
`t̃_j = (α_j m_i (R t_i + e_ji) + m_j t_j)/m̃_j`. -/
def vectorUpdate (mj a mi : ℝ) (R ti e tj : ℂ) : ℂ :=
  if 0 < massUpdate mj a mi then
    ((a * mi : ℝ) * (R * ti + e) + (mj : ℂ) * tj) / (massUpdate mj a mi : ℂ)
  else tj

variable {ι : Type*} [Fintype ι]

/-- **Eq. 6, memoryless.** With `t_i = 0` and `t_j = 0` (no memory), non-negative masses and
weights, and `‖e_ji‖ = ℓ_ij`: `Σ_j m̃_j ‖t̃_j‖ = Σ_j α_ij m_i ℓ_ij`, whatever the transports `R`.

[formula: iem-memoryless-cost] -/
theorem memorylessCost_eq {mi : ℝ} (hmi : 0 ≤ mi) {m a ℓ : ι → ℝ} (hm : ∀ j, 0 ≤ m j)
    (ha : ∀ j, 0 ≤ a j) (R e : ι → ℂ) (he : ∀ j, ‖e j‖ = ℓ j) :
    ∑ j, massUpdate (m j) (a j) mi * ‖vectorUpdate (m j) (a j) mi (R j) 0 (e j) 0‖ =
      ∑ j, a j * mi * ℓ j := by
  apply Finset.sum_congr rfl
  intro j _
  have hM : 0 ≤ massUpdate (m j) (a j) mi := add_nonneg (hm j) (mul_nonneg (ha j) hmi)
  unfold vectorUpdate
  split_ifs with h
  swap
  · have h0 : massUpdate (m j) (a j) mi = 0 := le_antisymm (not_lt.mp h) hM
    rw [h0, zero_mul]
    have : a j * mi = 0 := by
      unfold massUpdate at h0
      have := mul_nonneg (ha j) hmi
      linarith [hm j]
    rw [this, zero_mul]
  · simp only [mul_zero, zero_add, add_zero]
    rw [norm_div, norm_mul, he, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
      Real.norm_eq_abs, abs_of_nonneg hM, abs_of_nonneg (mul_nonneg (ha j) hmi)]
    field_simp

/-! ## Couplings -/

variable {X : Type*} [Fintype X] [DecidableEq X]

/-- `π` is a coupling of `μ` and `ν`. -/
structure IsCoupling (μ ν : X → ℝ) (π : X → X → ℝ) : Prop where
  nonneg : ∀ x y, 0 ≤ π x y
  row : ∀ x, ∑ y, π x y = μ x
  col : ∀ y, ∑ x, π x y = ν y

/-- The transport cost of a coupling. -/
def cost (d : X → X → ℝ) (π : X → X → ℝ) : ℝ := ∑ x, ∑ y, π x y * d x y

/-- The distribution after removing `i` (Eq. 5): `m_y + α_y m_i`, and `0` at `i`. -/
def newMass (m a : X → ℝ) (i : X) (y : X) : ℝ := m y + a y * m i - if y = i then m i else 0

/-- The plan of Fig. 9: keep every `m_y` in place, move `α_y m_i` from `i` to `y`. -/
def removalPlan (m a : X → ℝ) (i : X) (x y : X) : ℝ :=
  (if x = y ∧ x ≠ i then m x else 0) + (if x = i then a y * m i else 0)

theorem removalPlan_isCoupling {m a : X → ℝ} {i : X} (hm : ∀ x, 0 ≤ m x) (ha : ∀ x, 0 ≤ a x)
    (hai : a i = 0) (hsum : ∑ y, a y = 1) :
    IsCoupling m (newMass m a i) (removalPlan m a i) := by
  refine ⟨fun x y => ?_, fun x => ?_, fun y => ?_⟩
  · unfold removalPlan
    have := mul_nonneg (ha y) (hm i)
    split_ifs <;> linarith [hm x]
  · unfold removalPlan
    rw [Finset.sum_add_distrib]
    by_cases hx : x = i
    · subst hx
      simp [← Finset.sum_mul, hsum]
    · simp [hx]
  · unfold removalPlan newMass
    rw [Finset.sum_add_distrib]
    by_cases hy : y = i
    · subst hy
      simp [hai]
    · rw [Finset.sum_eq_single y (fun x _ hx => by simp [hx]) (by simp)]
      simp [hy]

theorem cost_removalPlan (d : X → X → ℝ) (hd : ∀ x, d x x = 0) (m a : X → ℝ) (i : X) :
    cost d (removalPlan m a i) = m i * ∑ y, a y * d i y := by
  unfold cost removalPlan
  simp only [add_mul, Finset.sum_add_distrib]
  have h1 : ∑ x, ∑ y, (if x = y ∧ x ≠ i then m x else 0) * d x y = 0 := by
    apply Finset.sum_eq_zero; intro x _
    apply Finset.sum_eq_zero; intro y _
    split_ifs with h
    · rw [h.1, hd, mul_zero]
    · rw [zero_mul]
  have h2 : ∑ x, ∑ y, (if x = i then a y * m i else 0) * d x y = ∑ y, a y * m i * d i y := by
    rw [Finset.sum_eq_single i (fun x _ hx => by simp [hx]) (by simp)]
    simp
  rw [h1, h2, zero_add, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro y _
  ring

/-- **Every coupling of the old and new distribution costs at least `m_i Σ α_y d(i, y)`**, the
lower bound of the 1-Lipschitz potential `d(i, ·)`. -/
theorem cost_ge (d : X → X → ℝ) (hd : ∀ x, d x x = 0) (htri : ∀ x y z, d x z ≤ d x y + d y z)
    {m a : X → ℝ} {i : X} {π : X → X → ℝ} (hπ : IsCoupling m (newMass m a i) π) :
    m i * ∑ y, a y * d i y ≤ cost d π := by
  have key : ∑ x, ∑ y, π x y * (d i y - d i x) ≤ cost d π := by
    unfold cost
    apply Finset.sum_le_sum; intro x _
    apply Finset.sum_le_sum; intro y _
    exact mul_le_mul_of_nonneg_left (by linarith [htri i x y]) (hπ.nonneg x y)
  have expand : ∑ x, ∑ y, π x y * (d i y - d i x) =
      ∑ y, newMass m a i y * d i y - ∑ x, m x * d i x := by
    simp only [mul_sub, Finset.sum_sub_distrib]
    congr 1
    · rw [Finset.sum_comm]
      apply Finset.sum_congr rfl; intro y _
      rw [← Finset.sum_mul, hπ.col]
    · apply Finset.sum_congr rfl; intro x _
      rw [← Finset.sum_mul, hπ.row]
  have hval : ∑ y, newMass m a i y * d i y - ∑ x, m x * d i x = m i * ∑ y, a y * d i y := by
    unfold newMass
    simp only [sub_mul, add_mul, Finset.sum_sub_distrib, Finset.sum_add_distrib, ite_mul,
      zero_mul, Finset.sum_ite_eq', Finset.mem_univ, ite_true, hd, mul_zero, sub_zero,
      Finset.mul_sum]
    ring_nf
    apply Finset.sum_congr rfl; intro y _; ring
  linarith

/-- **`iem-memoryless-cost`: the memoryless cost is the 1-Wasserstein distance** between the old
mass distribution and the new one: `Σ_y α_y m_i d(i, y)` is the least cost of a coupling of `m` and
`m + m_i α − m_i δ_i`. On the flattened star `d(i, j) = ℓ_ij`, so this is Eq. 6.

[formula: iem-memoryless-cost] -/
theorem wasserstein_eq (d : X → X → ℝ) (hd : ∀ x, d x x = 0)
    (htri : ∀ x y z, d x z ≤ d x y + d y z) {m a : X → ℝ} {i : X} (hm : ∀ x, 0 ≤ m x)
    (ha : ∀ x, 0 ≤ a x) (hai : a i = 0) (hsum : ∑ y, a y = 1) :
    IsLeast {c | ∃ π, IsCoupling m (newMass m a i) π ∧ cost d π = c}
      (∑ y, a y * m i * d i y) := by
  have e : ∑ y, a y * m i * d i y = m i * ∑ y, a y * d i y := by
    rw [Finset.mul_sum]; apply Finset.sum_congr rfl; intro y _; ring
  rw [e]
  refine ⟨⟨removalPlan m a i, removalPlan_isCoupling hm ha hai hsum, cost_removalPlan d hd m a i⟩,
    ?_⟩
  rintro c ⟨π, hπ, rfl⟩
  exact cost_ge d hd htri hπ

/-- **Single source: the coupling is forced.** The only coupling of `m_i δ_i` and `m_i α` is
`π x y = [x = i] α_y m_i`; its cost is `m_i Σ_y α_y d(i, y)`.

[formula: iem-memoryless-cost] -/
theorem coupling_single_source_unique {mi : ℝ} {a : X → ℝ} {i : X} {π : X → X → ℝ}
    (hπ : IsCoupling (fun x => if x = i then mi else 0) (fun y => a y * mi) π) :
    π = fun x y => if x = i then a y * mi else 0 := by
  have hrow0 : ∀ x, x ≠ i → ∀ y, π x y = 0 := by
    intro x hx y
    have hr := hπ.row x
    simp only [hx, ite_false] at hr
    exact (Finset.sum_eq_zero_iff_of_nonneg (fun y _ => hπ.nonneg x y)).mp hr y
      (Finset.mem_univ y)
  funext x y
  by_cases hx : x = i
  · subst hx
    have hc := hπ.col y
    rw [← Finset.add_sum_erase _ _ (Finset.mem_univ x)] at hc
    have : ∑ z ∈ Finset.univ.erase x, π z y = 0 :=
      Finset.sum_eq_zero (fun z hz => hrow0 z (Finset.ne_of_mem_erase hz) y)
    simp only [ite_true]
    linarith
  · simp only [hx, ite_false]
    exact hrow0 x hx y

end IntrinsicSimplification

end DDG
