import DdgProofs.Oracle.Geometry
import DdgProofs.Cotan

/-!
# Cut flow: the clamp of the crossing parameter

*Variational Surface Cutting* (Sharp, Crane, ACM TOG 2018); the contract is
`crates/ddg-variational-cuts/SPEC.md`, §3.2. Statements about the oracle's scalar-generic
definitions (`lean/oracle/DdgOracle/Formulas/Cuts.lean`) instantiated at `ℝ`.

* **`cut-edge-crossing`** (`Ddg.clampCrossing`, `Ddg.crossingMoved`): the clamp lies in
  `[δ, 1 − δ]` for `δ ≤ ½`, is the identity exactly on that interval, moves a parameter of
  `[0, 1]` by at most `δ`, and `crossingMoved` says exactly when it moved (`clampCrossing_mem`,
  `clampCrossing_eq_self_iff`, `abs_clampCrossing_sub_le`, `crossingMoved_eq_true_iff`).

The crossing parameter itself, the barycentric lengths and the areas of sub-triangles
(`cut-edge-crossing`, `cut-barycentric-length`, `cut-subtriangle-area`) have native definitions
and their theorems in `DDG/Geometry/Barycentric.lean` (library); the bridges from the
oracle's definitions are in `DdgProofs/Oracle/Geometry.lean`. (Library step S3.)
-/

noncomputable section

namespace DDG

namespace Cuts

open scoped DDG.Oracle
open Real

/-! ## The clamp of the crossing parameter -/

/-- The oracle's clamp at `ℝ`. -/
theorem clampCrossing_eq (δ t : ℝ) : Ddg.clampCrossing δ t = max δ (min (1 - δ) t) := by
  simp only [Ddg.clampCrossing, Oracle.one_eq]

/-- **The clamped parameter lies in `[δ, 1 − δ]`** (for `δ ≤ ½`, so that the interval is not
empty).

[formula: cut-edge-crossing] -/
theorem clampCrossing_mem {δ : ℝ} (hδ : δ ≤ 1 / 2) (t : ℝ) :
    δ ≤ Ddg.clampCrossing δ t ∧ Ddg.clampCrossing δ t ≤ 1 - δ := by
  rw [clampCrossing_eq]
  exact ⟨le_max_left _ _, max_le (by linarith) (min_le_left _ _)⟩

/-- For `0 < δ ≤ ½` the clamped parameter is strictly inside the edge. -/
theorem clampCrossing_mem_Ioo {δ : ℝ} (hδ0 : 0 < δ) (hδ : δ ≤ 1 / 2) (t : ℝ) :
    0 < Ddg.clampCrossing δ t ∧ Ddg.clampCrossing δ t < 1 := by
  obtain ⟨h1, h2⟩ := clampCrossing_mem hδ t
  constructor <;> linarith

/-- **The clamp is the identity exactly on `[δ, 1 − δ]`.**

[formula: cut-edge-crossing] -/
theorem clampCrossing_eq_self_iff {δ : ℝ} (hδ : δ ≤ 1 / 2) (t : ℝ) :
    Ddg.clampCrossing δ t = t ↔ δ ≤ t ∧ t ≤ 1 - δ := by
  constructor
  · intro h
    have := clampCrossing_mem hδ t
    rw [h] at this
    exact this
  · rintro ⟨h1, h2⟩
    rw [clampCrossing_eq, min_eq_right h2, max_eq_right h1]

/-- **`crossingMoved` reports exactly the crossings that the clamp moved.**

[formula: cut-edge-crossing] -/
theorem crossingMoved_eq_true_iff {δ : ℝ} (hδ : δ ≤ 1 / 2) (t : ℝ) :
    Ddg.crossingMoved δ t = true ↔ Ddg.clampCrossing δ t ≠ t := by
  rw [Ne, clampCrossing_eq_self_iff hδ]
  simp only [Ddg.crossingMoved, Bool.or_eq_true, decide_eq_true_eq, Oracle.one_eq]
  constructor
  · rintro (h | h) ⟨h1, h2⟩ <;> linarith
  · intro h
    by_contra hcon
    push Not at hcon
    exact h ⟨hcon.1, hcon.2⟩

/-- **The clamp moves a parameter of `[0, 1]` by at most `δ`.**

[formula: cut-edge-crossing] -/
theorem abs_clampCrossing_sub_le {δ t : ℝ} (hδ0 : 0 ≤ δ) (hδ : δ ≤ 1 / 2) (ht0 : 0 ≤ t)
    (ht1 : t ≤ 1) : |Ddg.clampCrossing δ t - t| ≤ δ := by
  rw [clampCrossing_eq, abs_le]
  rcases le_total (1 - δ) t with h | h
  · rw [min_eq_left h, max_eq_right (by linarith)]
    constructor <;> linarith
  · rw [min_eq_right h]
    rcases le_total δ t with h' | h'
    · rw [max_eq_right h']
      constructor <;> linarith
    · rw [max_eq_left h']
      constructor <;> linarith

end Cuts

end DDG
