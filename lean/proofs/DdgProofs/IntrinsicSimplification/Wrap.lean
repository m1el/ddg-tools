import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Algebra.Order.ToIntervalMod

/-!
# Intrinsic simplification: wrapping angular coordinates (S4)

`crates/ddg-intrinsic-simplification/SPEC.md` S4, §13 (`iem-wrap`); the Rust function
`formulas::wrap_angle`: at an interior vertex `φ ↦ φ.rem_euclid(2π) = φ − 2π⌊φ/2π⌋`, at a boundary
vertex the identity (the coordinates lie in `[0, π]` by construction).

* `wrapAngle_eq_toIcoMod`, `wrapAngle_mem` — the interior wrap is Mathlib's `toIcoMod` onto
  `[0, 2π)`;
* `exp_wrapAngle` — it keeps the direction `e^{iφ}`;
* `wrapAngle_of_mem` — it fixes `[0, 2π)`; the boundary wrap fixes everything, in particular
  `[0, π]`;
* `wrapAngle_ccw` — it keeps the counter-clockwise order of a fan: for `φ ≤ ψ < φ + 2π` (a fan of
  less than a turn, e.g. the accumulated `φ + s θ` of SPEC §3.2 and S4), the CCW angle from the
  wrapped `φ` to the wrapped `ψ`, `wrap(wrap ψ − wrap φ)`, is `ψ − φ`.

[formula: iem-wrap]
-/

noncomputable section

namespace DdgProofs

namespace IntrinsicSimplification

open Real

/-- `wrap_angle`: `φ − 2π⌊φ/2π⌋` inside, `φ` on the boundary. -/
def wrapAngle (φ : ℝ) (boundary : Bool) : ℝ :=
  if boundary then φ else φ - 2 * π * ⌊φ / (2 * π)⌋

/-- The interior wrap is `toIcoMod` onto `[0, 2π)`.

[formula: iem-wrap] -/
theorem wrapAngle_eq_toIcoMod (φ : ℝ) : wrapAngle φ false = toIcoMod two_pi_pos 0 φ := by
  have hp := two_pi_pos
  symm
  rw [toIcoMod_eq_iff]
  refine ⟨⟨?_, ?_⟩, ⌊φ / (2 * π)⌋, ?_⟩
  · simp only [wrapAngle, Bool.false_eq_true, ite_false]
    have := Int.floor_le (φ / (2 * π))
    have : (⌊φ / (2 * π)⌋ : ℝ) * (2 * π) ≤ φ := by rwa [le_div_iff₀ hp] at this
    linarith
  · simp only [wrapAngle, zero_add, Bool.false_eq_true, ite_false]
    have := Int.lt_floor_add_one (φ / (2 * π))
    have : φ < ((⌊φ / (2 * π)⌋ : ℝ) + 1) * (2 * π) := by rwa [div_lt_iff₀ hp] at this
    linarith
  · simp only [wrapAngle, zsmul_eq_mul, Bool.false_eq_true, ite_false]; ring

/-- [formula: iem-wrap] -/
theorem wrapAngle_mem (φ : ℝ) : 0 ≤ wrapAngle φ false ∧ wrapAngle φ false < 2 * π := by
  rw [wrapAngle_eq_toIcoMod]
  have := toIcoMod_mem_Ico two_pi_pos 0 φ
  rw [zero_add] at this
  exact this

/-- **The direction is kept**: `e^{i wrap(φ)} = e^{iφ}`.

[formula: iem-wrap] -/
theorem exp_wrapAngle (φ : ℝ) (boundary : Bool) :
    Complex.exp (wrapAngle φ boundary * Complex.I) = Complex.exp (φ * Complex.I) := by
  cases boundary
  · simp only [wrapAngle, Bool.false_eq_true, ite_false]
    push_cast
    rw [show ((φ : ℂ) - 2 * π * ⌊φ / (2 * π)⌋) * Complex.I =
      φ * Complex.I + (-(⌊φ / (2 * π)⌋ : ℤ) : ℤ) * (2 * π * Complex.I) by push_cast; ring,
      Complex.exp_add, Complex.exp_int_mul_two_pi_mul_I, mul_one]
  · simp [wrapAngle]

/-- The wrap fixes `[0, 2π)` inside, everything (so `[0, π]`) on the boundary.

[formula: iem-wrap] -/
theorem wrapAngle_of_mem {φ : ℝ} (h0 : 0 ≤ φ) (h1 : φ < 2 * π) (boundary : Bool) :
    wrapAngle φ boundary = φ := by
  cases boundary
  · rw [wrapAngle_eq_toIcoMod, toIcoMod_eq_self]; exact ⟨h0, by rwa [zero_add]⟩
  · rfl

/-- **The CCW order of a fan is kept**: for `φ ≤ ψ < φ + 2π`, the counter-clockwise angle from
`wrap φ` to `wrap ψ`, wrapped, is `ψ − φ`.

[formula: iem-wrap] -/
theorem wrapAngle_ccw {φ ψ : ℝ} (h0 : φ ≤ ψ) (h1 : ψ < φ + 2 * π) :
    wrapAngle (wrapAngle ψ false - wrapAngle φ false) false = ψ - φ := by
  simp only [wrapAngle_eq_toIcoMod]
  have eψ := toIcoMod_sub_self two_pi_pos 0 ψ
  have eφ := toIcoMod_sub_self two_pi_pos 0 φ
  have : toIcoMod two_pi_pos 0 ψ - toIcoMod two_pi_pos 0 φ =
      (ψ - φ) + (toIcoDiv two_pi_pos 0 φ - toIcoDiv two_pi_pos 0 ψ) • (2 * π) := by
    simp only [neg_smul] at eψ eφ
    rw [sub_smul]; linarith
  rw [this, toIcoMod_add_zsmul, toIcoMod_eq_self]
  constructor <;> linarith

end IntrinsicSimplification

end DdgProofs
