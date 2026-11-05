import DDG.Geometry.Plane
import DDG.Geometry.Stable

/-!
# Signed corner angle of three plane points

Transcription of `signed_corner_angle` in `crates/ddg-bff/src/quality.rs` (`atan2(u × v, u · v)`).

## `signed-angle-2d` — `signed_corner_angle`

* `signedCornerAngle_eq_arg` — it is `arg(conj U · V)`;
* `signedCornerAngle_mem_Ioc` — it lies in `(−π, π]`;
* `abs_signedCornerAngle_eq_angle` — its absolute value is the unsigned angle;
* `signedCornerAngle_sign` — its sign is the sign of `sda apex b c`; `0` or `π` if that is zero;
* `signedCornerAngle_sum_eq_pi` — in a counter-clockwise triangle the three angles are
  positive and sum to `π`.


Nothing is claimed about floating-point rounding.
-/

noncomputable section

namespace DDG

open Real

/-! ## Signed corner angle -/

/-- `signed_corner_angle`: `atan2(u × v, u · v)` with `u = b − apex`, `v = c − apex`, as
evaluated. -/
def signedCornerAngle (apex b c : ℝ × ℝ) : ℝ :=
  atan2 ((b.1 - apex.1) * (c.2 - apex.2) - (b.2 - apex.2) * (c.1 - apex.1))
    ((b.1 - apex.1) * (c.1 - apex.1) + (b.2 - apex.2) * (c.2 - apex.2))


open ComplexConjugate in
/-- `conj U · V = (u · v) + i (u × v)`. -/
theorem conj_toC_mul_toC (u v : ℝ × ℝ) : conj (toC u) * toC v = ⟨dot2 u v, cross2 u v⟩ := by
  apply Complex.ext <;> simp [toC, dot2, cross2]
  ring

/-- The signed corner angle is `arg` of the complex number `(u · v) + i (u × v)`, and
`u × v = sda apex b c`. -/
theorem signedCornerAngle_eq_arg_mk (apex b c : ℝ × ℝ) :
    signedCornerAngle apex b c =
      Complex.arg ⟨dot2 (b - apex) (c - apex), sda apex b c⟩ := rfl

open ComplexConjugate in
/-- **The signed corner angle is `arg(conj U · V)`** for the complex numbers `U = b − apex`,
`V = c − apex`: the angle of the rotation that takes the direction of `U` to that of `V`.

[formula: signed-angle-2d] -/
theorem signedCornerAngle_eq_arg (apex b c : ℝ × ℝ) :
    signedCornerAngle apex b c =
      Complex.arg (conj (toC b - toC apex) * (toC c - toC apex)) := by
  rw [← toC_sub, ← toC_sub, conj_toC_mul_toC]
  rfl

/-- **The signed corner angle lies in `(−π, π]`.**

[formula: signed-angle-2d] -/
theorem signedCornerAngle_mem_Ioc (apex b c : ℝ × ℝ) :
    -π < signedCornerAngle apex b c ∧ signedCornerAngle apex b c ≤ π := by
  rw [signedCornerAngle_eq_arg_mk]
  exact ⟨Complex.neg_pi_lt_arg _, Complex.arg_le_pi _⟩

/-- **The absolute value of the signed corner angle is the unsigned angle** between
`u = b − apex` and `v = c − apex` (both nonzero), measured in the Euclidean plane `ℂ`.

[formula: signed-angle-2d] -/
theorem abs_signedCornerAngle_eq_angle {apex b c : ℝ × ℝ} (hb : b ≠ apex) (hc : c ≠ apex) :
    |signedCornerAngle apex b c| =
      InnerProductGeometry.angle (toC b - toC apex) (toC c - toC apex) := by
  have hU : toC b - toC apex ≠ 0 := sub_ne_zero.mpr fun h => hb (toC_injective h)
  have hV : toC c - toC apex ≠ 0 := sub_ne_zero.mpr fun h => hc (toC_injective h)
  rw [signedCornerAngle_eq_arg, InnerProductGeometry.angle_comm, Complex.angle_eq_abs_arg hV hU]
  congr 1
  have h : (toC c - toC apex) / (toC b - toC apex) =
      (((Complex.normSq (toC b - toC apex))⁻¹ : ℝ) : ℂ) *
        ((starRingEnd ℂ) (toC b - toC apex) * (toC c - toC apex)) := by
    rw [div_eq_mul_inv, Complex.inv_def]
    ring
  rw [h, Complex.arg_real_mul _ (inv_pos.mpr (Complex.normSq_pos.mpr hU))]

/-- **The sign of the signed corner angle is the sign of `sda apex b c = u × v`**
(counter-clockwise positive):

* `0 < θ < π` iff `u × v > 0`;
* `θ < 0` iff `u × v < 0`;
* `θ = 0` iff `u × v = 0` and `u · v ≥ 0` (same direction, or one of the vectors is zero);
* `θ = π` iff `u × v = 0` and `u · v < 0` (opposite directions).

[formula: signed-angle-2d] -/
theorem signedCornerAngle_sign (apex b c : ℝ × ℝ) :
    (0 < signedCornerAngle apex b c ∧ signedCornerAngle apex b c < π ↔ 0 < sda apex b c) ∧
      (signedCornerAngle apex b c < 0 ↔ sda apex b c < 0) ∧
      (signedCornerAngle apex b c = 0 ↔
        sda apex b c = 0 ∧ 0 ≤ dot2 (b - apex) (c - apex)) ∧
      (signedCornerAngle apex b c = π ↔
        sda apex b c = 0 ∧ dot2 (b - apex) (c - apex) < 0) := by
  rw [signedCornerAngle_eq_arg_mk]
  generalize hz : (⟨dot2 (b - apex) (c - apex), sda apex b c⟩ : ℂ) = z
  have hre : dot2 (b - apex) (c - apex) = z.re := by rw [← hz]
  have him : sda apex b c = z.im := by rw [← hz]
  rw [hre, him]
  refine ⟨⟨?_, fun h => ?_⟩, Complex.arg_neg_iff, ?_, ?_⟩
  · rintro ⟨h0, hπ⟩
    rcases (Complex.arg_nonneg_iff.mp h0.le).eq_or_lt with h2 | h2
    · exfalso
      rcases le_or_gt 0 z.re with h3 | h3
      · exact h0.ne' (Complex.arg_eq_zero_iff.mpr ⟨h3, h2.symm⟩)
      · exact hπ.ne (Complex.arg_eq_pi_iff.mpr ⟨h3, h2.symm⟩)
    · exact h2
  · refine ⟨lt_of_le_of_ne (Complex.arg_nonneg_iff.mpr h.le) fun e => ?_,
      lt_of_le_of_ne (Complex.arg_le_pi z) fun e => ?_⟩
    · exact h.ne' (Complex.arg_eq_zero_iff.mp e.symm).2
    · exact h.ne' (Complex.arg_eq_pi_iff.mp e).2
  · rw [Complex.arg_eq_zero_iff, and_comm]
  · rw [Complex.arg_eq_pi_iff, and_comm]

open ComplexConjugate in
/-- **Angle sum of a counter-clockwise triangle**: if `sda a b c > 0`, the three signed corner
angles are positive and add up to `π`.

[formula: signed-angle-2d] -/
theorem signedCornerAngle_sum_eq_pi {a b c : ℝ × ℝ} (h : 0 < sda a b c) :
    0 < signedCornerAngle a b c ∧ 0 < signedCornerAngle b c a ∧ 0 < signedCornerAngle c a b ∧
      signedCornerAngle a b c + signedCornerAngle b c a + signedCornerAngle c a b = π := by
  have hb : 0 < sda b c a := by rw [sda_rotate]; exact h
  have hc : 0 < sda c a b := by rw [sda_rotate]; exact hb
  obtain ⟨h1, h1'⟩ := (signedCornerAngle_sign a b c).1.mpr h
  obtain ⟨h2, h2'⟩ := (signedCornerAngle_sign b c a).1.mpr hb
  obtain ⟨h3, h3'⟩ := (signedCornerAngle_sign c a b).1.mpr hc
  refine ⟨h1, h2, h3, ?_⟩
  -- the three complex numbers whose arguments are the angles
  have e1 := signedCornerAngle_eq_arg a b c
  have e2 := signedCornerAngle_eq_arg b c a
  have e3 := signedCornerAngle_eq_arg c a b
  set z1 : ℂ := conj (toC b - toC a) * (toC c - toC a) with hz1
  set z2 : ℂ := conj (toC c - toC b) * (toC a - toC b) with hz2
  set z3 : ℂ := conj (toC a - toC c) * (toC b - toC c) with hz3
  have n1 : z1 ≠ 0 := by
    intro h0
    rw [e1, h0, Complex.arg_zero] at h1
    exact lt_irrefl _ h1
  have n2 : z2 ≠ 0 := by
    intro h0
    rw [e2, h0, Complex.arg_zero] at h2
    exact lt_irrefl _ h2
  have n3 : z3 ≠ 0 := by
    intro h0
    rw [e3, h0, Complex.arg_zero] at h3
    exact lt_irrefl _ h3
  -- their product is a negative real number
  have hprod : z1 * z2 * z3 =
      ((-(Complex.normSq (toC b - toC a) * Complex.normSq (toC c - toC b) *
        Complex.normSq (toC a - toC c)) : ℝ) : ℂ) := by
    have m1 := Complex.mul_conj (toC b - toC a)
    have m2 := Complex.mul_conj (toC c - toC b)
    have m3 := Complex.mul_conj (toC a - toC c)
    push_cast
    rw [← m1, ← m2, ← m3, hz1, hz2, hz3]
    ring
  have hneg : -(Complex.normSq (toC b - toC a) * Complex.normSq (toC c - toC b) *
      Complex.normSq (toC a - toC c)) < 0 := by
    have p1 : 0 < Complex.normSq (toC b - toC a) :=
      Complex.normSq_pos.mpr fun h0 => n1 (by rw [hz1, h0, map_zero, zero_mul])
    have p2 : 0 < Complex.normSq (toC c - toC b) :=
      Complex.normSq_pos.mpr fun h0 => n2 (by rw [hz2, h0, map_zero, zero_mul])
    have p3 : 0 < Complex.normSq (toC a - toC c) :=
      Complex.normSq_pos.mpr fun h0 => n3 (by rw [hz3, h0, map_zero, zero_mul])
    exact neg_neg_of_pos (mul_pos (mul_pos p1 p2) p3)
  have hang : ((Complex.arg z1 + Complex.arg z2 + Complex.arg z3 : ℝ) : Real.Angle) =
      ((π : ℝ) : Real.Angle) := by
    rw [Real.Angle.coe_add, Real.Angle.coe_add, ← Complex.arg_mul_coe_angle n1 n2,
      ← Complex.arg_mul_coe_angle (mul_ne_zero n1 n2) n3, hprod, Complex.arg_ofReal_of_neg hneg]
  obtain ⟨k, hk⟩ := Real.Angle.angle_eq_iff_two_pi_dvd_sub.mp hang
  rw [e1] at h1 h1'
  rw [e2] at h2 h2'
  rw [e3] at h3 h3'
  rw [e1, e2, e3]
  have hk1 : (0 : ℝ) < 2 * k + 1 :=
    (mul_pos_iff_of_pos_left pi_pos).mp (by linarith)
  have hk2 : (0 : ℝ) < 2 - 2 * k :=
    (mul_pos_iff_of_pos_left pi_pos).mp (by linarith)
  have hk1' : (0 : ℤ) < 2 * k + 1 := by exact_mod_cast hk1
  have hk2' : (0 : ℤ) < 2 - 2 * k := by exact_mod_cast hk2
  have hk0 : k = 0 := by omega
  rw [hk0] at hk
  simp only [Int.cast_zero, mul_zero] at hk
  linarith

end DDG
