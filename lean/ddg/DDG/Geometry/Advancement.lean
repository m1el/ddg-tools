import Mathlib

/-!
# Repulsive Curves: conservative advancement (SPEC §7.4, D11)

A primitive (point, segment, triangle) moving linearly is the set of convex combinations of its
corners `p_i(τ) = p_i − τ a_i`; the distance of two primitives at time `τ` is the infimum over the
weights of both (`primDist`). In any real normed space (for `ℝ³` with the Euclidean norm,
`EuclideanSpace ℝ (Fin 3)`):

* `primDist_lipschitz`: if every corner speed is at most `l_A` (resp. `l_B`), then
  `dist(τ') ≥ dist(τ) − |τ' − τ| (l_A + l_B)`: the distance is Lipschitz with the sum of the
  largest corner displacements per unit `τ` on each side.
* `advance_safe`: if `dist(τ) ≥ dist(τ_k) − l (τ − τ_k)` for `τ ≥ τ_k` (`l > 0`) and
  `dist(τ_k) > g`, then with `τ_{k+1} = τ_k + (dist(τ_k) − g)/l` the distance stays `≥ g` on
  `[τ_k, τ_{k+1}]`, and `τ_{k+1} > τ_k`; with `g = s d₀ > 0` there is no contact before the advanced
  time.
-/

noncomputable section

namespace DDG

namespace RC

open Finset

/-- Barycentric weights on the corners of a primitive. -/
def Weights (ι : Type*) [Fintype ι] := {w : ι → ℝ // (∀ i, 0 ≤ w i) ∧ ∑ i, w i = 1}

instance (ι : Type*) [Fintype ι] [DecidableEq ι] [Nonempty ι] : Nonempty (Weights ι) := by
  obtain ⟨i⟩ := ‹Nonempty ι›
  exact ⟨⟨Pi.single i 1, fun j => by
    by_cases h : j = i
    · subst h; simp
    · simp [h], by simp⟩⟩

variable {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]

/-- The point of a moving primitive with weights `w` at time `τ`: `Σ_i w_i (p_i − τ a_i)`. -/
def primPoint {ι : Type*} [Fintype ι] (p a : ι → F) (w : Weights ι) (τ : ℝ) : F :=
  ∑ i, w.1 i • (p i - τ • a i)

/-- The distance of two linearly moving primitives at time `τ`. -/
def primDist {ι κ : Type*} [Fintype ι] [Fintype κ] (p a : ι → F) (q b : κ → F) (τ : ℝ) : ℝ :=
  ⨅ x : Weights ι × Weights κ, ‖primPoint p a x.1 τ - primPoint q b x.2 τ‖

theorem primPoint_move {ι : Type*} [Fintype ι] (p a : ι → F) (w : Weights ι) (τ τ' : ℝ) :
    primPoint p a w τ' = primPoint p a w τ - (τ' - τ) • ∑ i, w.1 i • a i := by
  simp only [primPoint, Finset.smul_sum, ← Finset.sum_sub_distrib, smul_sub, smul_smul]
  refine Finset.sum_congr rfl fun i _ => ?_
  module

theorem norm_weighted_le {ι : Type*} [Fintype ι] (a : ι → F) (w : Weights ι) {l : ℝ}
    (hl : ∀ i, ‖a i‖ ≤ l) : ‖∑ i, w.1 i • a i‖ ≤ l := by
  calc ‖∑ i, w.1 i • a i‖ ≤ ∑ i, ‖w.1 i • a i‖ := norm_sum_le _ _
    _ = ∑ i, w.1 i * ‖a i‖ := by
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg (w.2.1 i)]
    _ ≤ ∑ i, w.1 i * l := Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (hl i) (w.2.1 i)
    _ = l := by rw [← Finset.sum_mul, w.2.2, one_mul]

/-- **The distance of two moving primitives is Lipschitz in `τ`** with constant `l_A + l_B`, the
largest corner speeds of the two primitives.

[formula: rc-ccd-advancement] -/
theorem primDist_lipschitz {ι κ : Type*} [Fintype ι] [Fintype κ] [DecidableEq ι] [DecidableEq κ]
    [Nonempty ι] [Nonempty κ] (p a : ι → F) (q b : κ → F) {lA lB : ℝ}
    (ha : ∀ i, ‖a i‖ ≤ lA) (hb : ∀ j, ‖b j‖ ≤ lB) (τ τ' : ℝ) :
    primDist p a q b τ - |τ' - τ| * (lA + lB) ≤ primDist p a q b τ' := by
  have hbdd : ∀ σ, BddBelow (Set.range fun x : Weights ι × Weights κ =>
      ‖primPoint p a x.1 σ - primPoint q b x.2 σ‖) :=
    fun σ => ⟨0, by rintro _ ⟨x, rfl⟩; exact norm_nonneg _⟩
  refine le_ciInf fun x => ?_
  have h1 : primDist p a q b τ ≤ ‖primPoint p a x.1 τ - primPoint q b x.2 τ‖ :=
    ciInf_le (hbdd τ) x
  have h2 : ‖primPoint p a x.1 τ - primPoint q b x.2 τ‖ ≤
      ‖primPoint p a x.1 τ' - primPoint q b x.2 τ'‖ + |τ' - τ| * (lA + lB) := by
    rw [primPoint_move p a x.1 τ τ', primPoint_move q b x.2 τ τ']
    have e : primPoint p a x.1 τ - primPoint q b x.2 τ =
        (primPoint p a x.1 τ - (τ' - τ) • ∑ i, x.1.1 i • a i
          - (primPoint q b x.2 τ - (τ' - τ) • ∑ j, x.2.1 j • b j))
        + ((τ' - τ) • ∑ i, x.1.1 i • a i - (τ' - τ) • ∑ j, x.2.1 j • b j) := by abel
    rw [e]
    refine (norm_add_le _ _).trans (add_le_add le_rfl ?_)
    refine (norm_sub_le _ _).trans ?_
    rw [norm_smul, norm_smul, Real.norm_eq_abs, mul_add]
    exact add_le_add (mul_le_mul_of_nonneg_left (norm_weighted_le a x.1 ha) (abs_nonneg _))
      (mul_le_mul_of_nonneg_left (norm_weighted_le b x.2 hb) (abs_nonneg _))
  linarith

/-- **Conservative advancement** (Additive CCD, Li, Kaufman, Jiang 2021): if the distance satisfies
the Lipschitz lower bound `dist(τ) ≥ dist(τ_k) − l (τ − τ_k)` for `τ ≥ τ_k` (`l > 0`) and
`dist(τ_k) > g`, then `τ_{k+1} = τ_k + (dist(τ_k) − g)/l > τ_k` and `dist ≥ g` on
`[τ_k, τ_{k+1}]`; for `g > 0` (the gap `s d₀`) the primitives do not touch before `τ_{k+1}`.

[formula: rc-ccd-advancement] -/
theorem advance_safe {f : ℝ → ℝ} {τk l g : ℝ} (hl : 0 < l)
    (hlip : ∀ τ, τk ≤ τ → f τk - l * (τ - τk) ≤ f τ) (hg : g < f τk) :
    τk < τk + (f τk - g) / l ∧
      ∀ τ, τk ≤ τ → τ ≤ τk + (f τk - g) / l → g ≤ f τ ∧ (0 < g → 0 < f τ) := by
  refine ⟨by have : 0 < (f τk - g) / l := div_pos (by linarith) hl; linarith, fun τ h1 h2 => ?_⟩
  have h3 : l * (τ - τk) ≤ f τk - g := by
    have : τ - τk ≤ (f τk - g) / l := by linarith
    calc l * (τ - τk) ≤ l * ((f τk - g) / l) := mul_le_mul_of_nonneg_left this hl.le
      _ = f τk - g := by field_simp
  have := hlip τ h1
  exact ⟨by linarith, fun hg0 => by linarith⟩

/-- **No contact before the advanced time**, for two linearly moving primitives: combining the two
lemmas with `l = l_A + l_B > 0`.

[formula: rc-ccd-advancement] -/
theorem primDist_advance {ι κ : Type*} [Fintype ι] [Fintype κ] [DecidableEq ι] [DecidableEq κ]
    [Nonempty ι] [Nonempty κ] (p a : ι → F) (q b : κ → F) {lA lB : ℝ}
    (ha : ∀ i, ‖a i‖ ≤ lA) (hb : ∀ j, ‖b j‖ ≤ lB) (hl : 0 < lA + lB) {τk g : ℝ}
    (hg0 : 0 < g) (hg : g < primDist p a q b τk) :
    ∀ τ, τk ≤ τ → τ ≤ τk + (primDist p a q b τk - g) / (lA + lB) →
      0 < primDist p a q b τ := by
  intro τ h1 h2
  refine ((advance_safe hl (fun τ hτ => ?_) hg).2 τ h1 h2).2 hg0
  have := primDist_lipschitz p a q b ha hb τk τ
  rwa [abs_of_nonneg (by linarith), mul_comm] at this

end RC

end DDG
