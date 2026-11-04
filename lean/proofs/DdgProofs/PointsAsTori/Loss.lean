import DdgProofs.PointsAsTori.Baselines

/-!
# Points as Tori: the training loss (Eq. 27; SPEC §10.1)

`L = mean_{q ∈ S₁} |φ(q) − φ_true(q)| + mean_{q ∈ S₂} |1 − |∇φ(q)||` over the valid
(row, query) pairs (the code drops non-finite predictions from each mean, so the two index sets
may differ). The loss is nonnegative and vanishes exactly when the prediction is the true value
and has unit gradient on the valid queries (`loss_spec`).
-/

noncomputable section

namespace DdgProofs

namespace PAT

open Finset

variable {Q : Type*}

/-- The loss of Eq. 27 as coded: two L1 means over the valid sets `S₁`, `S₂`. -/
def trainLoss (S1 S2 : Finset Q) (φ φt gnorm : Q → ℝ) : ℝ :=
  (∑ q ∈ S1, |φ q - φt q|) / S1.card + (∑ q ∈ S2, |1 - gnorm q|) / S2.card

/-- **The loss** (Eq. 27, the code's rules): `L ≥ 0`, and for nonempty valid sets `L = 0` iff
`φ = φ_true` on `S₁` and `|∇φ| = 1` on `S₂`.

[formula: pat-loss] -/
theorem loss_spec (S1 S2 : Finset Q) (φ φt gnorm : Q → ℝ) :
    0 ≤ trainLoss S1 S2 φ φt gnorm ∧
    (S1.Nonempty → S2.Nonempty →
      (trainLoss S1 S2 φ φt gnorm = 0 ↔ (∀ q ∈ S1, φ q = φt q) ∧ ∀ q ∈ S2, gnorm q = 1)) := by
  have h1 : 0 ≤ (∑ q ∈ S1, |φ q - φt q|) / S1.card :=
    div_nonneg (Finset.sum_nonneg fun _ _ => abs_nonneg _) (Nat.cast_nonneg _)
  have h2 : 0 ≤ (∑ q ∈ S2, |1 - gnorm q|) / S2.card :=
    div_nonneg (Finset.sum_nonneg fun _ _ => abs_nonneg _) (Nat.cast_nonneg _)
  refine ⟨add_nonneg h1 h2, fun hS1 hS2 => ?_⟩
  have c1 : (0 : ℝ) < S1.card := by exact_mod_cast hS1.card_pos
  have c2 : (0 : ℝ) < S2.card := by exact_mod_cast hS2.card_pos
  unfold trainLoss
  rw [add_eq_zero_iff_of_nonneg h1 h2, div_eq_zero_iff, div_eq_zero_iff,
    Finset.sum_eq_zero_iff_of_nonneg (fun _ _ => abs_nonneg _),
    Finset.sum_eq_zero_iff_of_nonneg (fun _ _ => abs_nonneg _)]
  simp only [abs_eq_zero, sub_eq_zero, c1.ne', c2.ne', or_false]
  exact ⟨fun ⟨a, b⟩ => ⟨a, fun q hq => (b q hq).symm⟩, fun ⟨a, b⟩ => ⟨a, fun q hq => (b q hq).symm⟩⟩

end PAT

end DdgProofs
