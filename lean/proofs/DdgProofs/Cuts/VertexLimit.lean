import DdgProofs.Cuts.Quantile
import DdgProofs.Cuts.Regions

/-!
# Cut flow: the local limit at the vertices

`crates/ddg-variational-cuts/SPEC.md` §4 ("Step size", the second form of the local limit; "On a
network"), §8.4 ("Local limit"), D26.

With the length step on the mesh the limit of the robust statistics acts on the change of the
whole step at every vertex `i`, after the length step:

```
Δ_i = φ′(i) − φ(i)
v_i = |Δ_i| / τ                       (Capped)
v_i = max(|Δ_i| / τ, ½ (ω_i − ω)₊)    (Adaptive)
φ″(i) = φ(i) + min(1, max_motion · h_i / (τ · v_i)) · Δ_i        (v_i = 0: the factor is 1)
```

and for `N` regions with `Δ_i = ½ [(φ_b′ − φ_b) − (φ_a′ − φ_a)]` (`a` the region of the vertex,
`b` the smallest of the other functions) and the same factor for every function at the vertex.

* **`cut-step-local`**, vertex form — about the oracle's `vertexChange`, `vertexSpeedCapped`,
  `vertexSpeedAdaptive`, `vertexLimitFactor` (`= localLimitFactor`) and `limitedValue`:
  for `max_motion · h_i ≥ 0`, `τ > 0` and any speed `v_i ≥ |Δ_i|/τ`, the value changes by at most
  `max_motion · h_i` (`vertexLimit_change_le`, `_capped`, `_adaptive`); for
  `max_motion · h_i > 0` the change keeps its sign, does not grow, vanishes only if `Δ_i = 0`,
  and `φ″` lies between `φ` and `φ′` (`vertexLimit_sign`, `vertexLimit_between`); a vertex with
  `τ v_i ≤ max_motion · h_i` is unchanged, `φ″ = φ′` (`vertexLimit_unchanged`, `_capped`,
  `_adaptive`). The capped rule is the clamp of `Δ_i` to `±max_motion · h_i`, whatever `τ`
  (`vertexLimit_capped_clamp`, `vertexLimit_capped_limited`). The speeds are those of the cut
  points for the rate `Δ_i/τ` (`vertexSpeedCapped_eq_abs_rate`,
  `vertexSpeedAdaptive_eq_localSpeed`); without a cut point at the vertex, or with `ω_i ≤ ω`, the
  adaptive speed is the capped one (`vertexSpeedAdaptive_of_le`).
* `N` regions: with every function at the vertex taking the factor, half the difference of the
  two functions changes by the factor times `Δ_i` (`vertexChange_limited`), hence by at most
  `max_motion · h_i`, with its sign kept (`vertexLimit_regions_change_le`, `_capped`,
  `_adaptive`, `vertexLimit_regions_sign`); every difference of two functions changes by the
  factor times its unlimited change (`limitedValue_sub_limitedValue`). On the state `[φ, −φ]`
  the change is `±(φ′ − φ)` for either region, the speeds and so the factor are those of the
  two-region rule, and the limited state is `[φ″, −φ″]` with `φ″` the limited value of the
  two-region rule (`vertexChange_two_lower`, `vertexChange_two_upper`, `abs_vertexChange_two`,
  `vertexLimit_two_regions`).

All theorems are about the oracle's definitions at `ℝ`; nothing here is a definition of this
package. Which vertices are cut-adjacent, how `h_i` and `ω_i` are collected, and the order of the
stages are not modelled.
-/

noncomputable section

namespace DdgProofs

namespace Cuts

open scoped DdgProofs.Oracle
open Real

/-! ## The oracle's definitions at `ℝ` -/

/-- The oracle's `vertexChange`: `½ [(φ_b′ − φ_b) − (φ_a′ − φ_a)]`. -/
theorem vertexChange_eq (b' b a' a : ℝ) :
    Ddg.vertexChange b' b a' a = ((b' - b) - (a' - a)) / 2 := by
  simp only [Ddg.vertexChange, Oracle.two_eq]

/-- The oracle's `vertexSpeedCapped`: `|Δ|/τ`. -/
theorem vertexSpeedCapped_eq (Δ τ : ℝ) : Ddg.vertexSpeedCapped Δ τ = |Δ| / τ := by
  simp only [Ddg.vertexSpeedCapped, Oracle.abs_eq]

/-- The oracle's `vertexSpeedAdaptive`: `max(|Δ|/τ, ½ (ω_i − ω)₊)`. -/
theorem vertexSpeedAdaptive_eq (Δ τ ωi ω : ℝ) :
    Ddg.vertexSpeedAdaptive Δ τ ωi ω = max (|Δ| / τ) (max 0 (ωi - ω) / 2) := by
  simp only [Ddg.vertexSpeedAdaptive, Oracle.abs_eq, Oracle.zero_eq, Oracle.two_eq]

/-- The oracle's `vertexLimitFactor` is its `localLimitFactor`: "the same formula". -/
theorem vertexLimitFactor_eq (mm h τ v : ℝ) :
    Ddg.vertexLimitFactor mm h τ v = Ddg.localLimitFactor mm h τ v := rfl

/-- The oracle's `limitedValue`: `φ + f (φ′ − φ)`. -/
theorem limitedValue_eq (φ φ' f : ℝ) : Ddg.limitedValue φ φ' f = φ + f * (φ' - φ) := rfl

/-- The change of a limited value is the factor times the change of the whole step. -/
theorem limitedValue_sub (φ φ' f : ℝ) : Ddg.limitedValue φ φ' f - φ = (φ' - φ) * f := by
  rw [limitedValue_eq]
  ring

/-- With the factor `1` the value of the whole step is kept. -/
theorem limitedValue_one (φ φ' : ℝ) : Ddg.limitedValue φ φ' 1 = φ' := by
  rw [limitedValue_eq]
  ring

/-- With `v_i = 0` the factor is `1` (SPEC §4). -/
theorem vertexLimitFactor_zero_speed (mm h τ : ℝ) : Ddg.vertexLimitFactor mm h τ 0 = 1 := by
  rw [vertexLimitFactor_eq, localLimitFactor_def, mul_zero, ite_eq_right (lt_irrefl 0)]

/-! ## The speeds -/

/-- The capped speed of a vertex is the size of the rate `Δ_i/τ` (`τ > 0`). -/
theorem vertexSpeedCapped_eq_abs_rate {τ : ℝ} (hτ : 0 < τ) (Δ : ℝ) :
    Ddg.vertexSpeedCapped Δ τ = |Δ / τ| := by
  rw [vertexSpeedCapped_eq, abs_div, abs_of_pos hτ]

/-- **The adaptive speed of a vertex is the oracle's `localSpeed` of the rate `Δ_i/τ`**
(`τ > 0`): the rule at the vertices is the rule at the cut points with this rate. -/
theorem vertexSpeedAdaptive_eq_localSpeed {τ : ℝ} (hτ : 0 < τ) (Δ ωi ω : ℝ) :
    Ddg.vertexSpeedAdaptive Δ τ ωi ω = Ddg.localSpeed (Δ / τ) ωi ω := by
  rw [vertexSpeedAdaptive_eq, localSpeed_eq, abs_div, abs_of_pos hτ]

/-- The adaptive speed is at least the capped one. -/
theorem vertexSpeedCapped_le_adaptive (Δ τ ωi ω : ℝ) :
    Ddg.vertexSpeedCapped Δ τ ≤ Ddg.vertexSpeedAdaptive Δ τ ωi ω := by
  rw [vertexSpeedCapped_eq, vertexSpeedAdaptive_eq]
  exact le_max_left _ _

/-- **Where `ω_i ≤ ω` the adaptive speed is the capped one** (`τ ≥ 0`); in particular at a
vertex without a cut point on its edges, where the term is absent. -/
theorem vertexSpeedAdaptive_of_le {τ : ℝ} (hτ : 0 ≤ τ) (Δ : ℝ) {ωi ω : ℝ} (h : ωi ≤ ω) :
    Ddg.vertexSpeedAdaptive Δ τ ωi ω = Ddg.vertexSpeedCapped Δ τ := by
  rw [vertexSpeedAdaptive_eq, vertexSpeedCapped_eq, max_eq_left (by linarith : ωi - ω ≤ 0),
    zero_div, max_eq_left (div_nonneg (abs_nonneg _) hτ)]

/-- Taking `max(ω, largest ω_c)` for `ω_i` (as the oracle's two-region step does) gives the same
speed as the largest `ω_c` itself. -/
theorem vertexSpeedAdaptive_max_omega (Δ τ M ω : ℝ) :
    Ddg.vertexSpeedAdaptive Δ τ (max ω M) ω = Ddg.vertexSpeedAdaptive Δ τ M ω := by
  have e : max 0 (max ω M - ω) = max 0 (M - ω) := by
    rcases le_total M ω with h | h
    · rw [max_eq_left h, sub_self, max_self, max_eq_left (by linarith : M - ω ≤ 0)]
    · rw [max_eq_right h]
  rw [vertexSpeedAdaptive_eq, vertexSpeedAdaptive_eq, e]

/-- The speeds depend on the size of the change only. -/
theorem vertexSpeedCapped_neg (Δ τ : ℝ) :
    Ddg.vertexSpeedCapped (-Δ) τ = Ddg.vertexSpeedCapped Δ τ := by
  rw [vertexSpeedCapped_eq, vertexSpeedCapped_eq, abs_neg]

theorem vertexSpeedAdaptive_neg (Δ τ ωi ω : ℝ) :
    Ddg.vertexSpeedAdaptive (-Δ) τ ωi ω = Ddg.vertexSpeedAdaptive Δ τ ωi ω := by
  rw [vertexSpeedAdaptive_eq, vertexSpeedAdaptive_eq, abs_neg]

/-! ## Two regions -/

/-- **No value changes by more than `max_motion · h_i`**: for `max_motion · h_i ≥ 0`, `τ > 0`
and a speed `v_i ≥ |φ′ − φ|/τ`, `|φ″ − φ| ≤ max_motion · h_i`.

[formula: cut-step-local] -/
theorem vertexLimit_change_le {mm h τ v φ φ' : ℝ} (ha : 0 ≤ mm * h) (hτ : 0 < τ)
    (hv : |φ' - φ| / τ ≤ v) :
    |Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ v) - φ| ≤ mm * h := by
  rw [limitedValue_sub, vertexLimitFactor_eq]
  have hr : |(φ' - φ) / τ| ≤ v := by rwa [abs_div, abs_of_pos hτ]
  have h1 := localLimit_motion_le ha hτ.le hr
  have e : τ * |(φ' - φ) / τ * Ddg.localLimitFactor mm h τ v| =
      |(φ' - φ) * Ddg.localLimitFactor mm h τ v| := by
    rw [div_mul_eq_mul_div, abs_div, abs_of_pos hτ]
    field_simp
  rwa [e] at h1

/-- The same with the capped speed `v_i = |Δ_i|/τ`.

[formula: cut-step-local] -/
theorem vertexLimit_change_le_capped {mm h τ : ℝ} (ha : 0 ≤ mm * h) (hτ : 0 < τ) (φ φ' : ℝ) :
    |Ddg.limitedValue φ φ'
        (Ddg.vertexLimitFactor mm h τ (Ddg.vertexSpeedCapped (φ' - φ) τ)) - φ| ≤ mm * h :=
  vertexLimit_change_le ha hτ (le_of_eq (vertexSpeedCapped_eq _ _).symm)

/-- The same with the adaptive speed `v_i = max(|Δ_i|/τ, ½ (ω_i − ω)₊)`.

[formula: cut-step-local] -/
theorem vertexLimit_change_le_adaptive {mm h τ : ℝ} (ha : 0 ≤ mm * h) (hτ : 0 < τ)
    (φ φ' ωi ω : ℝ) :
    |Ddg.limitedValue φ φ'
        (Ddg.vertexLimitFactor mm h τ (Ddg.vertexSpeedAdaptive (φ' - φ) τ ωi ω)) - φ| ≤
      mm * h :=
  vertexLimit_change_le ha hτ
    ((le_of_eq (vertexSpeedCapped_eq _ _).symm).trans (vertexSpeedCapped_le_adaptive _ _ _ _))

/-- **The direction of the whole step is kept at every vertex** (`max_motion · h_i > 0`): the
limited change has the sign of `φ′ − φ`, vanishes only if that does, and is not larger.

[formula: cut-step-local] -/
theorem vertexLimit_sign {mm h : ℝ} (ha : 0 < mm * h) (τ v φ φ' : ℝ) :
    (0 < Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ v) - φ ↔ 0 < φ' - φ) ∧
      (Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ v) - φ < 0 ↔ φ' - φ < 0) ∧
      (Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ v) = φ ↔ φ' = φ) ∧
      |Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ v) - φ| ≤ |φ' - φ| := by
  obtain ⟨h1, h2, h3, h4⟩ := localLimit_sign ha τ v (φ' - φ)
  rw [limitedValue_sub, vertexLimitFactor_eq]
  refine ⟨h1, h2, ?_, h4⟩
  rw [← sub_eq_zero, limitedValue_sub, h3, sub_eq_zero]

/-- The limited value lies between the value before the step and the value of the whole step
(`max_motion · h_i > 0`). -/
theorem vertexLimit_between {mm h : ℝ} (ha : 0 < mm * h) (τ v φ φ' : ℝ) :
    min φ φ' ≤ Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ v) ∧
      Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ v) ≤ max φ φ' := by
  obtain ⟨hf0, hf1⟩ := localLimitFactor_mem_Ioc ha τ v
  rw [vertexLimitFactor_eq, limitedValue_eq]
  rcases le_total φ φ' with hle | hle
  · rw [min_eq_left hle, max_eq_right hle]
    constructor
    · nlinarith [mul_nonneg hf0.le (sub_nonneg.mpr hle)]
    · nlinarith [mul_nonneg (sub_nonneg.mpr hf1) (sub_nonneg.mpr hle)]
  · rw [min_eq_right hle, max_eq_left hle]
    constructor
    · nlinarith [mul_nonneg (sub_nonneg.mpr hf1) (sub_nonneg.mpr hle)]
    · nlinarith [mul_nonneg hf0.le (sub_nonneg.mpr hle)]

/-- **A vertex within the limit is unchanged**: if `τ v_i ≤ max_motion · h_i` then `φ″ = φ′`.

[formula: cut-step-local] -/
theorem vertexLimit_unchanged {mm h τ v : ℝ} (hle : τ * v ≤ mm * h) (φ φ' : ℝ) :
    Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ v) = φ' := by
  rw [vertexLimitFactor_eq, localLimitFactor_eq_one hle, limitedValue_one]

/-- Capped: a vertex whose change is at most `max_motion · h_i` is unchanged (`τ > 0`).

[formula: cut-step-local] -/
theorem vertexLimit_unchanged_capped {mm h τ φ φ' : ℝ} (hτ : 0 < τ)
    (hle : |φ' - φ| ≤ mm * h) :
    Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ (Ddg.vertexSpeedCapped (φ' - φ) τ)) =
      φ' := by
  apply vertexLimit_unchanged
  have e : τ * (|φ' - φ| / τ) = |φ' - φ| := by field_simp
  rw [vertexSpeedCapped_eq, e]
  exact hle

/-- Adaptive: a vertex whose change is at most `max_motion · h_i`, and for which
`τ · ½ (ω_i − ω) ≤ max_motion · h_i`, is unchanged (`τ > 0`).

[formula: cut-step-local] -/
theorem vertexLimit_unchanged_adaptive {mm h τ φ φ' ωi ω : ℝ} (hτ : 0 < τ)
    (hle : |φ' - φ| ≤ mm * h) (hω : τ * ((ωi - ω) / 2) ≤ mm * h) :
    Ddg.limitedValue φ φ'
        (Ddg.vertexLimitFactor mm h τ (Ddg.vertexSpeedAdaptive (φ' - φ) τ ωi ω)) = φ' := by
  apply vertexLimit_unchanged
  have e : τ * (|φ' - φ| / τ) = |φ' - φ| := by field_simp
  have h0 : 0 ≤ mm * h := (abs_nonneg _).trans hle
  rw [vertexSpeedAdaptive_eq, mul_max_of_nonneg _ _ hτ.le, e]
  refine max_le hle ?_
  rcases le_total (ωi - ω) 0 with hn | hp
  · rw [max_eq_left hn, zero_div, mul_zero]
    exact h0
  · rw [max_eq_right hp]
    exact hω

/-- **The capped rule at a vertex is a clamp** (`max_motion · h_i ≥ 0`, `τ > 0`):
`φ″ = φ + max(−max_motion · h_i, min(max_motion · h_i, φ′ − φ))`. It does not depend on `τ`. -/
theorem vertexLimit_capped_clamp {mm h τ : ℝ} (ha : 0 ≤ mm * h) (hτ : 0 < τ) (φ φ' : ℝ) :
    Ddg.limitedValue φ φ' (Ddg.vertexLimitFactor mm h τ (Ddg.vertexSpeedCapped (φ' - φ) τ)) =
      φ + max (-(mm * h)) (min (mm * h) (φ' - φ)) := by
  by_cases hle : |φ' - φ| ≤ mm * h
  · rw [vertexLimit_unchanged_capped hτ hle]
    obtain ⟨h1, h2⟩ := abs_le.mp hle
    rw [min_eq_right h2, max_eq_right h1]
    ring
  · have hlt : mm * h < |φ' - φ| := not_le.mp hle
    have hpos : 0 < |φ' - φ| := lt_of_le_of_lt ha hlt
    have e : τ * (|φ' - φ| / τ) = |φ' - φ| := by field_simp
    rw [vertexLimitFactor_eq, vertexSpeedCapped_eq, localLimitFactor_def, e, ite_eq_left hpos,
      min_eq_right ((div_le_one hpos).mpr hlt.le), limitedValue_eq]
    rcases lt_or_gt_of_ne (abs_pos.mp hpos) with hneg | hp
    · have hlt' : mm * h < -(φ' - φ) := by rwa [abs_of_neg hneg] at hlt
      rw [abs_of_neg hneg, min_eq_right (by linarith), max_eq_left (by linarith)]
      have hne : φ' - φ ≠ 0 := hneg.ne
      field_simp
    · have hlt' : mm * h < φ' - φ := by rwa [abs_of_pos hp] at hlt
      rw [abs_of_pos hp, min_eq_left hlt'.le, max_eq_right (by linarith)]
      have hne : φ' - φ ≠ 0 := hp.ne'
      field_simp

/-- **A limited vertex of the capped rule changes by exactly `max_motion · h_i`.** -/
theorem vertexLimit_capped_limited {mm h τ φ φ' : ℝ} (ha : 0 ≤ mm * h) (hτ : 0 < τ)
    (hlt : mm * h < |φ' - φ|) :
    |Ddg.limitedValue φ φ'
        (Ddg.vertexLimitFactor mm h τ (Ddg.vertexSpeedCapped (φ' - φ) τ)) - φ| = mm * h := by
  rw [vertexLimit_capped_clamp ha hτ, add_sub_cancel_left]
  rcases le_total 0 (φ' - φ) with hp | hn
  · rw [abs_of_nonneg hp] at hlt
    rw [min_eq_left hlt.le, max_eq_right (by linarith), abs_of_nonneg ha]
  · rw [abs_of_nonpos hn] at hlt
    rw [min_eq_right (by linarith), max_eq_left (by linarith), abs_neg, abs_of_nonneg ha]

/-! ## `N` regions -/

/-- **All functions at the vertex take the factor**: half the difference of the two functions
then changes by the factor times its unlimited change. -/
theorem vertexChange_limited (b' b a' a f : ℝ) :
    Ddg.vertexChange (Ddg.limitedValue b b' f) b (Ddg.limitedValue a a' f) a =
      Ddg.vertexChange b' b a' a * f := by
  rw [vertexChange_eq, vertexChange_eq, limitedValue_eq, limitedValue_eq]
  ring

/-- Every difference of two functions at the vertex changes by the factor times its unlimited
change. -/
theorem limitedValue_sub_limitedValue (x x' y y' f : ℝ) :
    Ddg.limitedValue x x' f - Ddg.limitedValue y y' f = (x - y) + f * ((x' - x) - (y' - y)) := by
  rw [limitedValue_eq, limitedValue_eq]
  ring

/-- The change of each single function is scaled by the factor, so for `0 ≤ f ≤ 1` it does not
grow. (It is not bounded by `max_motion · h_i`; only half the difference of the two is.) -/
theorem abs_limitedValue_sub_le {f : ℝ} (h0 : 0 ≤ f) (h1 : f ≤ 1) (x x' : ℝ) :
    |Ddg.limitedValue x x' f - x| ≤ |x' - x| := by
  rw [limitedValue_sub, abs_mul, abs_of_nonneg h0]
  exact mul_le_of_le_one_right (abs_nonneg _) h1

/-- **The limited change of half the difference of the two functions is at most
`max_motion · h_i`** (`max_motion · h_i ≥ 0`, `τ > 0`, a speed `v_i ≥ |Δ_i|/τ`).

[formula: cut-step-local] -/
theorem vertexLimit_regions_change_le {mm h τ v b' b a' a : ℝ} (ha : 0 ≤ mm * h) (hτ : 0 < τ)
    (hv : |Ddg.vertexChange b' b a' a| / τ ≤ v) :
    |Ddg.vertexChange (Ddg.limitedValue b b' (Ddg.vertexLimitFactor mm h τ v)) b
        (Ddg.limitedValue a a' (Ddg.vertexLimitFactor mm h τ v)) a| ≤ mm * h := by
  rw [vertexChange_limited]
  have hv' : |Ddg.vertexChange b' b a' a - 0| / τ ≤ v := by rwa [sub_zero]
  have h1 := vertexLimit_change_le (φ := 0) (φ' := Ddg.vertexChange b' b a' a) (v := v) ha hτ hv'
  rwa [limitedValue_sub, sub_zero] at h1

/-- The same with the capped speed.

[formula: cut-step-local] -/
theorem vertexLimit_regions_change_le_capped {mm h τ : ℝ} (ha : 0 ≤ mm * h) (hτ : 0 < τ)
    (b' b a' a : ℝ) :
    |Ddg.vertexChange
        (Ddg.limitedValue b b' (Ddg.vertexLimitFactor mm h τ
          (Ddg.vertexSpeedCapped (Ddg.vertexChange b' b a' a) τ))) b
        (Ddg.limitedValue a a' (Ddg.vertexLimitFactor mm h τ
          (Ddg.vertexSpeedCapped (Ddg.vertexChange b' b a' a) τ))) a| ≤ mm * h :=
  vertexLimit_regions_change_le ha hτ (le_of_eq (vertexSpeedCapped_eq _ _).symm)

/-- The same with the adaptive speed.

[formula: cut-step-local] -/
theorem vertexLimit_regions_change_le_adaptive {mm h τ : ℝ} (ha : 0 ≤ mm * h) (hτ : 0 < τ)
    (b' b a' a ωi ω : ℝ) :
    |Ddg.vertexChange
        (Ddg.limitedValue b b' (Ddg.vertexLimitFactor mm h τ
          (Ddg.vertexSpeedAdaptive (Ddg.vertexChange b' b a' a) τ ωi ω))) b
        (Ddg.limitedValue a a' (Ddg.vertexLimitFactor mm h τ
          (Ddg.vertexSpeedAdaptive (Ddg.vertexChange b' b a' a) τ ωi ω))) a| ≤ mm * h :=
  vertexLimit_regions_change_le ha hτ
    ((le_of_eq (vertexSpeedCapped_eq _ _).symm).trans (vertexSpeedCapped_le_adaptive _ _ _ _))

/-- **The limit keeps the sign of the change of half the difference** (`max_motion · h_i > 0`),
does not increase it, and keeps a vanishing one.

[formula: cut-step-local] -/
theorem vertexLimit_regions_sign {mm h : ℝ} (ha : 0 < mm * h) (τ v b' b a' a : ℝ) :
    (0 < Ddg.vertexChange (Ddg.limitedValue b b' (Ddg.vertexLimitFactor mm h τ v)) b
        (Ddg.limitedValue a a' (Ddg.vertexLimitFactor mm h τ v)) a ↔
        0 < Ddg.vertexChange b' b a' a) ∧
      (Ddg.vertexChange (Ddg.limitedValue b b' (Ddg.vertexLimitFactor mm h τ v)) b
        (Ddg.limitedValue a a' (Ddg.vertexLimitFactor mm h τ v)) a < 0 ↔
        Ddg.vertexChange b' b a' a < 0) ∧
      (Ddg.vertexChange (Ddg.limitedValue b b' (Ddg.vertexLimitFactor mm h τ v)) b
        (Ddg.limitedValue a a' (Ddg.vertexLimitFactor mm h τ v)) a = 0 ↔
        Ddg.vertexChange b' b a' a = 0) ∧
      |Ddg.vertexChange (Ddg.limitedValue b b' (Ddg.vertexLimitFactor mm h τ v)) b
        (Ddg.limitedValue a a' (Ddg.vertexLimitFactor mm h τ v)) a| ≤
        |Ddg.vertexChange b' b a' a| := by
  rw [vertexChange_limited, vertexLimitFactor_eq]
  exact localLimit_sign ha τ v _

/-! ## The state `[φ, −φ]` -/

/-- On `[φ, −φ]` with the vertex in region 0 (`a = 0`, `b = 1`): `Δ = −(φ′ − φ)`. -/
theorem vertexChange_two_lower (φ φ' : ℝ) : Ddg.vertexChange (-φ') (-φ) φ' φ = -(φ' - φ) := by
  rw [vertexChange_eq]
  ring

/-- On `[φ, −φ]` with the vertex in region 1 (`a = 1`, `b = 0`): `Δ = φ′ − φ`. -/
theorem vertexChange_two_upper (φ φ' : ℝ) : Ddg.vertexChange φ' φ (-φ') (-φ) = φ' - φ := by
  rw [vertexChange_eq]
  ring

/-- The limited value is odd: `limitedValue (−φ) (−φ′) f = −limitedValue φ φ′ f`. -/
theorem limitedValue_neg (φ φ' f : ℝ) :
    Ddg.limitedValue (-φ) (-φ') f = -Ddg.limitedValue φ φ' f := by
  rw [limitedValue_eq, limitedValue_eq]
  ring

/-- On `[φ, −φ]` the change of half the difference has the size `|φ′ − φ|`, whichever of the two
functions is the one of the region. -/
theorem abs_vertexChange_two (φ φ' : ℝ) {a b : Fin 2} (hab : a ≠ b) :
    |Ddg.vertexChange ((![φ', -φ'] : Fin 2 → ℝ) b) ((![φ, -φ] : Fin 2 → ℝ) b)
        ((![φ', -φ'] : Fin 2 → ℝ) a) ((![φ, -φ] : Fin 2 → ℝ) a)| = |φ' - φ| := by
  fin_cases a <;> fin_cases b
  · exact absurd rfl hab
  · show |Ddg.vertexChange (-φ') (-φ) φ' φ| = |φ' - φ|
    rw [vertexChange_two_lower, abs_neg]
  · show |Ddg.vertexChange φ' φ (-φ') (-φ)| = |φ' - φ|
    rw [vertexChange_two_upper]
  · exact absurd rfl hab

/-- **For `N = 2` the rule for `N` regions is the rule for two regions, and the state stays of
the form `[φ, −φ]`.** On `[φ, −φ]` with the whole step `[φ′, −φ′]`, for `a` the region of the
vertex and `b` the other function: both speeds are those of the two-region rule for
`Δ = φ′ − φ` (so the factor is the same), and with any factor the limited state is
`[φ″, −φ″]`, `φ″ = limitedValue φ φ′ f`.

[formula: cut-step-local] -/
theorem vertexLimit_two_regions (φ φ' : ℝ) {a b : Fin 2} (hab : a ≠ b) (τ ωi ω : ℝ) :
    Ddg.vertexSpeedCapped
        (Ddg.vertexChange ((![φ', -φ'] : Fin 2 → ℝ) b) ((![φ, -φ] : Fin 2 → ℝ) b)
          ((![φ', -φ'] : Fin 2 → ℝ) a) ((![φ, -φ] : Fin 2 → ℝ) a)) τ =
        Ddg.vertexSpeedCapped (φ' - φ) τ ∧
      Ddg.vertexSpeedAdaptive
        (Ddg.vertexChange ((![φ', -φ'] : Fin 2 → ℝ) b) ((![φ, -φ] : Fin 2 → ℝ) b)
          ((![φ', -φ'] : Fin 2 → ℝ) a) ((![φ, -φ] : Fin 2 → ℝ) a)) τ ωi ω =
        Ddg.vertexSpeedAdaptive (φ' - φ) τ ωi ω ∧
      ∀ f : ℝ, (fun k : Fin 2 =>
          Ddg.limitedValue ((![φ, -φ] : Fin 2 → ℝ) k) ((![φ', -φ'] : Fin 2 → ℝ) k) f) =
        ![Ddg.limitedValue φ φ' f, -Ddg.limitedValue φ φ' f] := by
  have e := abs_vertexChange_two φ φ' hab
  refine ⟨?_, ?_, fun f => ?_⟩
  · rw [vertexSpeedCapped_eq, vertexSpeedCapped_eq, e]
  · rw [vertexSpeedAdaptive_eq, vertexSpeedAdaptive_eq, e]
  · funext k
    fin_cases k
    · rfl
    · show Ddg.limitedValue (-φ) (-φ') f = -Ddg.limitedValue φ φ' f
      exact limitedValue_neg φ φ' f

end Cuts

end DdgProofs
