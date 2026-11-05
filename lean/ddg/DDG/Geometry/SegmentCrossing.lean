import DDG.Geometry.Plane

/-!
# Proper crossing of two plane segments

Transcription of `proper_crossing` in `crates/ddg-bff/src/quality.rs`.

## `segment-proper-crossing` — `proper_crossing`

* `strictlyOpposite_zero_iff`, `strictlyOpposite_imp` — the sign test without and with a
  tolerance;
* `properCrossing_sound` — if both tests hold, the returned point is interior to both segments;
* `properCrossing_complete` — an interior common point of non-parallel segments passes both
  tests, and it is the returned point (so it is unique);
* `properCrossing_iff` — both as one equivalence.


Nothing is claimed about floating-point rounding.
-/

noncomputable section

namespace DDG

open Real

/-! ## Proper crossing of two segments -/

/-- The test `opposite(x, y, length)` of `proper_crossing`, with the margin as a parameter:
`(x > margin ∧ y < −margin) ∨ (x < −margin ∧ y > margin)`. In the Rust code
`margin = TOUCH_TOLERANCE · max(|a₁ − a₀|, |b₁ − b₀|) · length ≥ 0`. -/
def strictlyOpposite (x y margin : ℝ) : Prop :=
  (x > margin ∧ y < -margin) ∨ (x < -margin ∧ y > margin)

/-- Without tolerance, the test says that `x` and `y` have strictly opposite signs.

[formula: segment-proper-crossing] -/
theorem strictlyOpposite_zero_iff (x y : ℝ) : strictlyOpposite x y 0 ↔ x * y < 0 := by
  unfold strictlyOpposite
  rw [mul_neg_iff, neg_zero]

/-- With any tolerance `margin ≥ 0`, the test still implies strictly opposite signs: the
tolerance only removes crossings.

[formula: segment-proper-crossing] -/
theorem strictlyOpposite_imp {x y margin : ℝ} (hm : 0 ≤ margin) (h : strictlyOpposite x y margin) :
    x * y < 0 := by
  rcases h with ⟨hx, hy⟩ | ⟨hx, hy⟩
  · exact mul_neg_iff.mpr (Or.inl ⟨by linarith, by linarith⟩)
  · exact mul_neg_iff.mpr (Or.inr ⟨by linarith, by linarith⟩)

/-- The parameter `t = d₃/(d₃ − d₄)` of the crossing point on the segment `a₀a₁`, with
`d₃ = sda b₀ b₁ a₀`, `d₄ = sda b₀ b₁ a₁`. -/
def crossingParam (a0 a1 b0 b1 : ℝ × ℝ) : ℝ :=
  sda b0 b1 a0 / (sda b0 b1 a0 - sda b0 b1 a1)

/-- The point `proper_crossing` returns: `a₀ + t (a₁ − a₀)`, coordinate by coordinate. -/
def crossingPoint (a0 a1 b0 b1 : ℝ × ℝ) : ℝ × ℝ :=
  (a0.1 + crossingParam a0 a1 b0 b1 * (a1.1 - a0.1),
    a0.2 + crossingParam a0 a1 b0 b1 * (a1.2 - a0.2))

private theorem ratio_mem_Ioo {x y : ℝ} (h : x * y < 0) : 0 < x / (x - y) ∧ x / (x - y) < 1 := by
  rcases mul_neg_iff.mp h with ⟨hx, hy⟩ | ⟨hx, hy⟩
  · have hd : 0 < x - y := by linarith
    exact ⟨div_pos hx hd, (div_lt_one hd).mpr (by linarith)⟩
  · have hd : x - y < 0 := by linarith
    exact ⟨div_pos_of_neg_of_neg hx hd, (div_lt_one_of_neg hd).mpr (by linarith)⟩

/-- `d₃ − d₄ = (a₁ − a₀) × (b₁ − b₀)` and `d₁ − d₂ = −(a₁ − a₀) × (b₁ − b₀)`. -/
theorem sda_sub_sda (a0 a1 b0 b1 : ℝ × ℝ) :
    sda b0 b1 a0 - sda b0 b1 a1 = cross2 (a1 - a0) (b1 - b0) ∧
      sda a0 a1 b0 - sda a0 a1 b1 = -cross2 (a1 - a0) (b1 - b0) := by
  simp only [sda, cross2, Prod.fst_sub, Prod.snd_sub]
  constructor <;> ring

/-- **Soundness of `proper_crossing`.** If `b₀, b₁` lie strictly on opposite sides of the line
`a₀a₁` (`d₁·d₂ < 0`) and `a₀, a₁` strictly on opposite sides of the line `b₀b₁` (`d₃·d₄ < 0`),
then with `t = d₃/(d₃ − d₄)` and `s = d₁/(d₁ − d₂)`: `0 < t < 1`, `0 < s < 1`, and the returned
point is `a₀ + t(a₁ − a₀) = b₀ + s(b₁ − b₀)`, an interior point of both segments.

The rounding tolerance `TOUCH_TOLERANCE` of the Rust code is not modelled; by
`strictlyOpposite_imp` the test with tolerance implies the hypotheses here, so this theorem
covers every point the Rust code returns (in exact arithmetic).

[formula: segment-proper-crossing] -/
theorem properCrossing_sound {a0 a1 b0 b1 : ℝ × ℝ}
    (h12 : sda a0 a1 b0 * sda a0 a1 b1 < 0) (h34 : sda b0 b1 a0 * sda b0 b1 a1 < 0) :
    (0 < crossingParam a0 a1 b0 b1 ∧ crossingParam a0 a1 b0 b1 < 1) ∧
      (0 < crossingParam b0 b1 a0 a1 ∧ crossingParam b0 b1 a0 a1 < 1) ∧
      crossingPoint a0 a1 b0 b1 = a0 + crossingParam a0 a1 b0 b1 • (a1 - a0) ∧
      crossingPoint a0 a1 b0 b1 = b0 + crossingParam b0 b1 a0 a1 • (b1 - b0) := by
  refine ⟨ratio_mem_Ioo h34, ratio_mem_Ioo h12, ?_, ?_⟩
  · ext <;> simp [crossingPoint]
  · have hd34 : sda b0 b1 a0 - sda b0 b1 a1 ≠ 0 := by
      intro h0
      have : sda b0 b1 a0 = sda b0 b1 a1 := by linarith
      rw [this] at h34
      exact absurd h34 (not_lt.mpr (mul_self_nonneg _))
    have hd12 : sda a0 a1 b0 - sda a0 a1 b1 ≠ 0 := by
      intro h0
      have : sda a0 a1 b0 = sda a0 a1 b1 := by linarith
      rw [this] at h12
      exact absurd h12 (not_lt.mpr (mul_self_nonneg _))
    have e := (sda_sub_sda a0 a1 b0 b1).1
    have e' := (sda_sub_sda a0 a1 b0 b1).2
    have hΔ : cross2 (a1 - a0) (b1 - b0) ≠ 0 := e ▸ hd34
    unfold crossingPoint crossingParam
    rw [e, e']
    ext <;> simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul,
      Prod.fst_sub, Prod.snd_sub] <;> field_simp <;>
      simp only [sda, cross2, Prod.fst_sub, Prod.snd_sub] <;> ring

/-- **Completeness of `proper_crossing` (without tolerance), and uniqueness.** If the segments
are not parallel and have a common point `a₀ + t(a₁ − a₀) = b₀ + s(b₁ − b₀)` with `0 < t < 1`,
`0 < s < 1`, then `d₁·d₂ < 0` and `d₃·d₄ < 0`, and `t = d₃/(d₃ − d₄)`, `s = d₁/(d₁ − d₂)`: the
common point is unique and is the point the function returns.

The rounding tolerance `TOUCH_TOLERANCE` of the Rust code is not modelled: with it, a crossing
at which an end point is within the tolerance of the other segment's line is not reported.

[formula: segment-proper-crossing] -/
theorem properCrossing_complete {a0 a1 b0 b1 : ℝ × ℝ} {t s : ℝ} (ht : 0 < t ∧ t < 1)
    (hs : 0 < s ∧ s < 1) (hP : a0 + t • (a1 - a0) = b0 + s • (b1 - b0))
    (hpar : cross2 (a1 - a0) (b1 - b0) ≠ 0) :
    sda a0 a1 b0 * sda a0 a1 b1 < 0 ∧ sda b0 b1 a0 * sda b0 b1 a1 < 0 ∧
      t = crossingParam a0 a1 b0 b1 ∧ s = crossingParam b0 b1 a0 a1 ∧
      crossingPoint a0 a1 b0 b1 = a0 + t • (a1 - a0) := by
  have hx : a0.1 + t * (a1.1 - a0.1) = b0.1 + s * (b1.1 - b0.1) := by
    simpa using congrArg Prod.fst hP
  have hy : a0.2 + t * (a1.2 - a0.2) = b0.2 + s * (b1.2 - b0.2) := by
    simpa using congrArg Prod.snd hP
  obtain ⟨e, e'⟩ := sda_sub_sda a0 a1 b0 b1
  have h3 : sda b0 b1 a0 = t * cross2 (a1 - a0) (b1 - b0) := by
    simp only [sda, cross2, Prod.fst_sub, Prod.snd_sub]
    linear_combination (b1.1 - b0.1) * hy - (b1.2 - b0.2) * hx
  have h1 : sda a0 a1 b0 = -(s * cross2 (a1 - a0) (b1 - b0)) := by
    simp only [sda, cross2, Prod.fst_sub, Prod.snd_sub]
    linear_combination (a1.2 - a0.2) * hx - (a1.1 - a0.1) * hy
  have h4 : sda b0 b1 a1 = -((1 - t) * cross2 (a1 - a0) (b1 - b0)) := by
    linear_combination h3 - e
  have h2 : sda a0 a1 b1 = (1 - s) * cross2 (a1 - a0) (b1 - b0) := by
    linear_combination h1 - e'
  have hΔ : 0 < cross2 (a1 - a0) (b1 - b0) ^ 2 := by positivity
  have ht' : t = crossingParam a0 a1 b0 b1 := by
    unfold crossingParam
    rw [e, h3]
    field_simp
  have hs' : s = crossingParam b0 b1 a0 a1 := by
    unfold crossingParam
    rw [e', h1]
    field_simp
  refine ⟨?_, ?_, ht', hs', ?_⟩
  · have : sda a0 a1 b0 * sda a0 a1 b1 =
        -(s * (1 - s) * cross2 (a1 - a0) (b1 - b0) ^ 2) := by
      rw [h1, h2]
      ring
    rw [this]
    exact neg_neg_of_pos (mul_pos (mul_pos hs.1 (sub_pos.mpr hs.2)) hΔ)
  · have : sda b0 b1 a0 * sda b0 b1 a1 =
        -(t * (1 - t) * cross2 (a1 - a0) (b1 - b0) ^ 2) := by
      rw [h3, h4]
      ring
    rw [this]
    exact neg_neg_of_pos (mul_pos (mul_pos ht.1 (sub_pos.mpr ht.2)) hΔ)
  · rw [ht']
    ext <;> simp [crossingPoint]

/-- **`proper_crossing` (without tolerance) decides proper crossing**: `b₀, b₁` lie strictly on
opposite sides of the line `a₀a₁` and `a₀, a₁` strictly on opposite sides of the line `b₀b₁`
iff the segments are not parallel and have a common point in the interior of both.

The rounding tolerance `TOUCH_TOLERANCE` of the Rust code is not modelled.

[formula: segment-proper-crossing] -/
theorem properCrossing_iff (a0 a1 b0 b1 : ℝ × ℝ) :
    (sda a0 a1 b0 * sda a0 a1 b1 < 0 ∧ sda b0 b1 a0 * sda b0 b1 a1 < 0) ↔
      cross2 (a1 - a0) (b1 - b0) ≠ 0 ∧
        ∃ t s : ℝ, (0 < t ∧ t < 1) ∧ (0 < s ∧ s < 1) ∧
          a0 + t • (a1 - a0) = b0 + s • (b1 - b0) := by
  constructor
  · rintro ⟨h12, h34⟩
    obtain ⟨ht, hs, hP, hP'⟩ := properCrossing_sound h12 h34
    refine ⟨?_, _, _, ht, hs, hP.symm.trans hP'⟩
    rw [← (sda_sub_sda a0 a1 b0 b1).1]
    intro h0
    have : sda b0 b1 a0 = sda b0 b1 a1 := by linarith
    rw [this] at h34
    exact absurd h34 (not_lt.mpr (mul_self_nonneg _))
  · rintro ⟨hpar, t, s, ht, hs, hP⟩
    obtain ⟨h12, h34, -⟩ := properCrossing_complete ht hs hP hpar
    exact ⟨h12, h34⟩

/-- The hypotheses of `properCrossing_sound` are satisfiable: the diagonals of the square
`[0,2]²` cross at its centre. -/
example :
    sda ((0, 0) : ℝ × ℝ) (2, 2) (0, 2) * sda ((0, 0) : ℝ × ℝ) (2, 2) (2, 0) < 0 ∧
      sda ((0, 2) : ℝ × ℝ) (2, 0) (0, 0) * sda ((0, 2) : ℝ × ℝ) (2, 0) (2, 2) < 0 ∧
      crossingPoint ((0, 0) : ℝ × ℝ) (2, 2) (0, 2) (2, 0) = (1, 1) := by
  refine ⟨by norm_num [sda], by norm_num [sda], ?_⟩
  ext <;> norm_num [crossingPoint, crossingParam, sda]

end DDG
