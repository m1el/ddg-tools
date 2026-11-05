import Mathlib.Analysis.SpecialFunctions.Complex.Arg

/-!
# Input-aligned frames, the angles (`ic-input-aligned-frame`, part 1 of 4)

`crates/ddg-integer-coords/SPEC.md` §4.12 (I26); the Rust `angle_from` and `wrap` of
`crates/ddg-integer-coords/src/frames.rs`.

The layouts of input faces are in the plane `ℝ × ℝ` (as in `Strip.lean`); `toC` reads a point
as a complex number.

* `angleFrom u d` is the Rust `angle_from(u, d) = atan2(u × d, u · d)`: for real numbers
  `atan2(y, x) = arg (x + i y)` (Mathlib's `Complex.arg`, in `(−π, π]`; `atan2(0, 0) = 0 = arg 0`;
  the signed zeros of floating point are not modelled). `conj u · d = u · d + i (u × d)`
  (`conj_mul_toC`).
* **It is the angle that turns `u` into `d`** (`angleFrom_spec`): `‖u‖ d = ‖d‖ e^{iα} u`, and it
  is the only such angle in `(−π, π]` (`angleFrom_unique`). So `α` is the counter-clockwise angle
  from `u` to `d` in the layout, not just a formula.
* **On the left of `u`** (`0 ≤ u × d`), `α ∈ [0, π]` (`angleFrom_mem_Icc`).
* **Inside a corner** spanned counter-clockwise by `u` and `y` (`u × y > 0`), a direction
  `d = λu + μy`, `λ, μ ≥ 0`, `d ≠ 0`, has `0 ≤ α ≤ angleFrom u y` (`angleFrom_le_of_cone`): the
  clamp of the Rust code to `[0, corner]` does nothing; and `angleFrom u y` is
  `arccos (u·y / (‖u‖‖y‖))` (`angleFrom_eq_arccos`), the corner angle of the lengths
  (`Frames/Input.lean`).
* `wrap` is the Rust `wrap` (`rem_euclid(2π)`, whose guard `≥ 2π ↦ 0` is vacuous in exact
  arithmetic): `toIcoMod` into `[0, 2π)`; it changes an angle by a multiple of `2π`
  (`wrap_exp`), is the identity on `[0, 2π)` and sends `2π` to `0`.
-/

noncomputable section

namespace DDG.NormalCoords

namespace Frames

open Complex

/-- A point (or vector) of the layout plane as a complex number. -/
def toC (p : ℝ × ℝ) : ℂ := ⟨p.1, p.2⟩

/-- `u · d`. -/
def dot (u d : ℝ × ℝ) : ℝ := u.1 * d.1 + u.2 * d.2

/-- `u × d` (the Rust `u[0] * d[1] - u[1] * d[0]`). -/
def cross (u d : ℝ × ℝ) : ℝ := u.1 * d.2 - u.2 * d.1

/-- The Rust `angle_from(u, d) = atan2(u × d, u · d)`, the angle from `u` to `d` in `(−π, π]`. -/
def angleFrom (u d : ℝ × ℝ) : ℝ := Complex.arg ⟨dot u d, cross u d⟩

/-- The Rust `wrap`: `rem_euclid(2π)`, into `[0, 2π)`. -/
def wrap (x : ℝ) : ℝ := toIcoMod Real.two_pi_pos 0 x

lemma conj_mul_toC (u d : ℝ × ℝ) :
    (starRingEnd ℂ) (toC u) * toC d = (⟨dot u d, cross u d⟩ : ℂ) := by
  apply Complex.ext
  · simp [toC, dot]
  · simp [toC, cross]; ring

lemma toC_ne_zero {u : ℝ × ℝ} (h : u ≠ 0) : toC u ≠ 0 := by
  intro e
  apply h
  have h1 : u.1 = 0 := by simpa [toC] using congrArg Complex.re e
  have h2 : u.2 = 0 := by simpa [toC] using congrArg Complex.im e
  exact Prod.ext h1 h2

lemma norm_pair (u d : ℝ × ℝ) : ‖(⟨dot u d, cross u d⟩ : ℂ)‖ = ‖toC u‖ * ‖toC d‖ := by
  rw [← conj_mul_toC, norm_mul, Complex.norm_conj]

lemma toC_mul_conj (u : ℝ × ℝ) : toC u * (starRingEnd ℂ) (toC u) = ((‖toC u‖ ^ 2 : ℝ) : ℂ) := by
  rw [Complex.mul_conj, Complex.sq_norm]

/-- **`α` turns `u` into `d`**: `‖u‖ d = ‖d‖ e^{iα} u`. [formula: ic-input-aligned-frame] -/
theorem angleFrom_spec (u d : ℝ × ℝ) :
    (‖toC u‖ : ℂ) * toC d = ‖toC d‖ * cexp (angleFrom u d * I) * toC u := by
  have key := Complex.norm_mul_exp_arg_mul_I (⟨dot u d, cross u d⟩ : ℂ)
  rw [norm_pair, ← conj_mul_toC] at key
  unfold angleFrom
  rw [← conj_mul_toC]
  by_cases hu : toC u = 0
  · rw [hu]; simp
  have hn : (‖toC u‖ : ℂ) ≠ 0 := by exact_mod_cast (norm_ne_zero_iff.mpr hu)
  apply mul_left_cancel₀ hn
  have e2 := congrArg (fun z => toC u * z) key
  have e3 : toC u * ((starRingEnd ℂ) (toC u) * toC d) = (‖toC u‖ : ℂ) ^ 2 * toC d := by
    rw [← mul_assoc, toC_mul_conj]; push_cast; ring
  rw [e3] at e2
  push_cast at e2 ⊢
  linear_combination -e2

/-- **The only such angle**: an angle `α ∈ (−π, π]` with `‖u‖ d = ‖d‖ e^{iα} u` (`u, d ≠ 0`) is
`angleFrom u d`. [formula: ic-input-aligned-frame] -/
theorem angleFrom_unique {u d : ℝ × ℝ} (hu : u ≠ 0) (hd : d ≠ 0) {α : ℝ}
    (hα : α ∈ Set.Ioc (-Real.pi) Real.pi)
    (h : (‖toC u‖ : ℂ) * toC d = ‖toC d‖ * cexp (α * I) * toC u) : angleFrom u d = α := by
  have hu' := toC_ne_zero hu
  have hd' := toC_ne_zero hd
  have hnu : 0 < ‖toC u‖ := norm_pos_iff.mpr hu'
  have hnd : 0 < ‖toC d‖ := norm_pos_iff.mpr hd'
  have e : (⟨dot u d, cross u d⟩ : ℂ) =
      ((‖toC u‖ * ‖toC d‖ : ℝ) : ℂ) * (Complex.cos α + Complex.sin α * I) := by
    rw [← conj_mul_toC, ← Complex.exp_mul_I]
    have hn : (‖toC u‖ : ℂ) ≠ 0 := by exact_mod_cast hnu.ne'
    apply mul_left_cancel₀ hn
    have e2 := congrArg (fun z => (starRingEnd ℂ) (toC u) * z) h
    have e3 : (starRingEnd ℂ) (toC u) * toC u = ((‖toC u‖ ^ 2 : ℝ) : ℂ) := by
      rw [mul_comm, toC_mul_conj]
    push_cast at e3 ⊢
    linear_combination e2 + (‖toC d‖ : ℂ) * cexp (α * I) * e3
  unfold angleFrom
  rw [e]
  exact Complex.arg_mul_cos_add_sin_mul_I (mul_pos hnu hnd) hα

/-- **On the left of `u`** the angle lies in `[0, π]`. -/
theorem angleFrom_mem_Icc {u d : ℝ × ℝ} (h : 0 ≤ cross u d) :
    0 ≤ angleFrom u d ∧ angleFrom u d ≤ Real.pi :=
  ⟨Complex.arg_nonneg_iff.mpr h, Complex.arg_le_pi _⟩

lemma dot_comb (u y : ℝ × ℝ) (l m : ℝ) :
    dot u (l * u.1 + m * y.1, l * u.2 + m * y.2) = l * (u.1 ^ 2 + u.2 ^ 2) + m * dot u y := by
  simp only [dot]; ring

lemma cross_comb (u y : ℝ × ℝ) (l m : ℝ) :
    cross u (l * u.1 + m * y.1, l * u.2 + m * y.2) = m * cross u y := by
  simp only [cross]; ring

/-- **Inside a corner the angle is at most the corner**: for `u × y > 0` and
`d = λu + μy ≠ 0` with `λ, μ ≥ 0`, `0 ≤ angleFrom u d ≤ angleFrom u y`.
[formula: ic-input-aligned-frame] -/
theorem angleFrom_le_of_cone {u y d : ℝ × ℝ} (hy : 0 < cross u y) {l m : ℝ} (hl : 0 ≤ l)
    (hm : 0 ≤ m) (hd : d = (l * u.1 + m * y.1, l * u.2 + m * y.2)) (hne : d ≠ 0) :
    0 ≤ angleFrom u d ∧ angleFrom u d ≤ angleFrom u y := by
  have hcd : cross u d = m * cross u y := by rw [hd, cross_comb]
  have hdd : dot u d = l * (u.1 ^ 2 + u.2 ^ 2) + m * dot u y := by rw [hd, dot_comb]
  have h0 : 0 ≤ cross u d := by rw [hcd]; positivity
  refine ⟨(angleFrom_mem_Icc h0).1, ?_⟩
  have hu : u ≠ 0 := by
    rintro rfl; simp [cross] at hy
  set zd : ℂ := ⟨dot u d, cross u d⟩ with hzd
  set zy : ℂ := ⟨dot u y, cross u y⟩ with hzy
  have hy0 : y ≠ 0 := by
    rintro rfl; simp [cross] at hy
  have nzd : 0 < ‖zd‖ := by
    rw [hzd, norm_pair]
    exact mul_pos (norm_pos_iff.mpr (toC_ne_zero hu)) (norm_pos_iff.mpr (toC_ne_zero hne))
  have nzy : 0 < ‖zy‖ := by
    rw [hzy, norm_pair]
    exact mul_pos (norm_pos_iff.mpr (toC_ne_zero hu)) (norm_pos_iff.mpr (toC_ne_zero hy0))
  have zd0 : zd ≠ 0 := norm_pos_iff.mp nzd
  have zy0 : zy ≠ 0 := norm_pos_iff.mp nzy
  set A := zd.arg with hA
  set B := zy.arg with hB
  show A ≤ B
  -- `sin (B − A) ‖zd‖ ‖zy‖ = λ ‖u‖² (u × y) ≥ 0`
  have hsin : Real.sin (B - A) * (‖zd‖ * ‖zy‖) = l * (u.1 ^ 2 + u.2 ^ 2) * cross u y := by
    rw [Real.sin_sub, hA, hB, Complex.sin_arg, Complex.sin_arg, Complex.cos_arg zd0,
      Complex.cos_arg zy0]
    field_simp
    simp only [hzd, hzy, hcd, hdd]
    ring
  have hs : 0 ≤ Real.sin (B - A) := by
    have : 0 ≤ Real.sin (B - A) * (‖zd‖ * ‖zy‖) := by rw [hsin]; positivity
    exact nonneg_of_mul_nonneg_left this (mul_pos nzd nzy)
  have hB0 : 0 < B := by
    rcases (Complex.arg_nonneg_iff.mpr (show (0 : ℝ) ≤ zy.im from hy.le)).lt_or_eq with h | h
    · exact h
    · exfalso
      have := (Complex.arg_eq_zero_iff.mp h.symm).2
      simp [hzy] at this; linarith
  have hApi : A ≤ Real.pi := Complex.arg_le_pi _
  by_contra hlt
  push Not at hlt
  have : Real.sin (B - A) < 0 :=
    Real.sin_neg_of_neg_of_neg_pi_lt (by linarith) (by linarith)
  linarith

/-- **The corner angle**: for `u × y ≥ 0`, `u, y ≠ 0`,
`angleFrom u y = arccos (u·y / (‖u‖‖y‖))`. -/
theorem angleFrom_eq_arccos {u y : ℝ × ℝ} (hu : u ≠ 0) (hy0 : y ≠ 0) (hy : 0 ≤ cross u y) :
    angleFrom u y = Real.arccos (dot u y / (‖toC u‖ * ‖toC y‖)) := by
  have h0 : (⟨dot u y, cross u y⟩ : ℂ) ≠ 0 := by
    rw [← norm_pos_iff, norm_pair]
    exact mul_pos (norm_pos_iff.mpr (toC_ne_zero hu)) (norm_pos_iff.mpr (toC_ne_zero hy0))
  have h := Complex.cos_arg h0
  rw [norm_pair] at h
  unfold angleFrom
  rw [← h, Real.arccos_cos (Complex.arg_nonneg_iff.mpr hy) (Complex.arg_le_pi _)]

/-! ## `wrap` -/

theorem wrap_mem (x : ℝ) : 0 ≤ wrap x ∧ wrap x < 2 * Real.pi := by
  have := toIcoMod_mem_Ico Real.two_pi_pos 0 x
  simp only [zero_add, Set.mem_Ico] at this
  exact this

theorem wrap_eq_self {x : ℝ} (h0 : 0 ≤ x) (h1 : x < 2 * Real.pi) : wrap x = x :=
  (toIcoMod_eq_self Real.two_pi_pos).mpr ⟨h0, by simpa using h1⟩

theorem wrap_two_pi : wrap (2 * Real.pi) = 0 := by
  have h := toIcoMod_add_right Real.two_pi_pos (0 : ℝ) (0 : ℝ)
  rw [zero_add] at h
  rw [wrap, h]
  exact (toIcoMod_eq_self Real.two_pi_pos).mpr ⟨le_rfl, by linarith [Real.pi_pos]⟩

/-- `wrap` changes an angle by a multiple of `2π`: the same direction `e^{ix}`. -/
theorem wrap_exp (x : ℝ) : cexp ((wrap x : ℝ) * I) = cexp (x * I) := by
  have h := self_sub_toIcoMod Real.two_pi_pos 0 x
  have e : wrap x = x - (toIcoDiv Real.two_pi_pos 0 x : ℝ) * (2 * Real.pi) := by
    rw [wrap]; rw [zsmul_eq_mul] at h; linarith
  rw [e]
  push_cast
  rw [sub_mul, Complex.exp_sub]
  have : cexp ((toIcoDiv Real.two_pi_pos 0 x : ℂ) * (2 * Real.pi) * I) = 1 := by
    have := Complex.exp_int_mul_two_pi_mul_I (toIcoDiv Real.two_pi_pos 0 x)
    rw [← this]; congr 1; ring
  rw [this, div_one]

end Frames

end DDG.NormalCoords
