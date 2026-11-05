import Mathlib

/-!
# Scalar convex facts behind the optimality system

Paper: Soliman, Slepčev, Crane, "Optimal Cone Singularities for Conformal Flattening" (2018),
Eqs. 18–19, 21 and Props. B.2–B.3.

Everything here is about one vertex: `μ` is the cone angle, `φ` the adjoint state,
`t = λ·wR ≥ 0` the threshold, `[lo, hi] ∋ 0` the allowed range of the cone angle
(`lo = ⊥`, `hi = ⊤` mean "unbounded") and `γ > 0` a regularization weight.

The subdifferential is spelled out: `IsSubgradientWithinAt f K g x` says
`g ∈ ∂(f + ι_K)(x)`, i.e. `x ∈ K` and `f x + g (y − x) ≤ f y` for every `y ∈ K`.
-/

namespace DDG

/-- The soft threshold of paper Eq. 18: `P_t(v) = max(0, v − t) + min(0, v + t)`.
The threshold `t` is `λ` in the unweighted problem and `λ·wR_i` in the weighted one. -/
def softThreshold (t v : ℝ) : ℝ := max 0 (v - t) + min 0 (v + t)

/-- The closed interval `[lo, hi] ⊆ ℝ` with extended-real end points
(`lo = ⊥`: no lower bound, `hi = ⊤`: no upper bound). -/
def boundSet (lo hi : EReal) : Set ℝ := {x | lo ≤ (x : EReal) ∧ (x : EReal) ≤ hi}

/-- Clamp of a real number to `[lo, hi]` with extended-real end points: `max(lo, min(hi, x))`.
For `lo = ⊥`, `hi = ⊤` it is the identity. (Only meaningful when `lo < ⊤` and `⊥ < hi`.) -/
noncomputable def clampE (lo hi : EReal) (x : ℝ) : ℝ := (max lo (min hi (x : EReal))).toReal

/-- `g ∈ ∂(f + ι_K)(x)`: `g` is a subgradient at `x` of `f` restricted to the set `K`.
Spelled out: `x ∈ K` and `f x + g·(y − x) ≤ f y` for every `y ∈ K`. -/
def IsSubgradientWithinAt (f : ℝ → ℝ) (K : Set ℝ) (g x : ℝ) : Prop :=
  x ∈ K ∧ ∀ y ∈ K, f x + g * (y - x) ≤ f y

/-- A subgradient inclusion is the same as a minimization: `g ∈ ∂(f + ι_K)(x)` iff
`x` minimizes `y ↦ f y − g·y` over `K`. -/
theorem isSubgradientWithinAt_iff_minimizes (f : ℝ → ℝ) (K : Set ℝ) (g x : ℝ) :
    IsSubgradientWithinAt f K g x ↔ x ∈ K ∧ ∀ y ∈ K, f x - g * x ≤ f y - g * y := by
  unfold IsSubgradientWithinAt
  constructor
  · rintro ⟨hx, h⟩
    exact ⟨hx, fun y hy => by have := h y hy; linarith⟩
  · rintro ⟨hx, h⟩
    exact ⟨hx, fun y hy => by have := h y hy; linarith⟩

/-! ## The soft threshold, piece by piece -/

theorem softThreshold_of_gt {t v : ℝ} (ht : 0 ≤ t) (h : t < v) : softThreshold t v = v - t := by
  unfold softThreshold
  rw [max_eq_right (by linarith), min_eq_left (by linarith)]; ring

theorem softThreshold_of_lt {t v : ℝ} (ht : 0 ≤ t) (h : v < -t) : softThreshold t v = v + t := by
  unfold softThreshold
  rw [max_eq_left (by linarith), min_eq_right (by linarith)]; ring

theorem softThreshold_of_abs_le {t v : ℝ} (h : |v| ≤ t) : softThreshold t v = 0 := by
  unfold softThreshold
  obtain ⟨h1, h2⟩ := abs_le.mp h
  rw [max_eq_left (by linarith), min_eq_left (by linarith)]; ring

/-- With zero threshold (`wR = 0`, a free cone site) the soft threshold is the identity. -/
theorem softThreshold_zero (v : ℝ) : softThreshold 0 v = v := by
  unfold softThreshold
  rcases le_total 0 v with h | h
  · rw [max_eq_right (by linarith), min_eq_left (by linarith)]; ring
  · rw [max_eq_left (by linarith), min_eq_right (by linarith)]; ring

/-- Closed form of the soft threshold (paper Prop. B.3, case form). -/
theorem softThreshold_eq_ite {t : ℝ} (ht : 0 ≤ t) (v : ℝ) :
    softThreshold t v = if t < v then v - t else if v < -t then v + t else 0 := by
  split_ifs with h1 h2
  · exact softThreshold_of_gt ht h1
  · exact softThreshold_of_lt ht h2
  · exact softThreshold_of_abs_le (abs_le.mpr ⟨by linarith, by linarith⟩)

/-! ## Clamping to an interval with extended-real end points -/

theorem clampE_of_mem {lo hi : EReal} {x : ℝ} (h1 : lo ≤ (x : EReal)) (h2 : (x : EReal) ≤ hi) :
    clampE lo hi x = x := by
  unfold clampE
  rw [min_eq_right h2, max_eq_right h1, EReal.toReal_coe]

theorem clampE_bot_top (x : ℝ) : clampE ⊥ ⊤ x = x := clampE_of_mem bot_le le_top

theorem clampE_of_hi_le {lo hi : EReal} {x : ℝ} (hlo : lo ≤ 0) (hhi : 0 ≤ hi)
    (h : hi ≤ (x : EReal)) : ∃ r : ℝ, hi = (r : EReal) ∧ clampE lo hi x = r := by
  induction hi using EReal.rec with
  | bot => exact absurd hhi (by simp)
  | coe r =>
    refine ⟨r, rfl, ?_⟩
    unfold clampE
    rw [min_eq_left h, max_eq_right (le_trans hlo hhi), EReal.toReal_coe]
  | top => exact absurd (top_le_iff.mp h) (EReal.coe_ne_top x)

theorem clampE_of_le_lo {lo hi : EReal} {x : ℝ} (hlo : lo ≤ 0) (hhi : 0 ≤ hi)
    (h : (x : EReal) ≤ lo) : ∃ r : ℝ, lo = (r : EReal) ∧ clampE lo hi x = r := by
  induction lo using EReal.rec with
  | bot => exact absurd (le_bot_iff.mp h) (EReal.coe_ne_bot x)
  | coe r =>
    refine ⟨r, rfl, ?_⟩
    unfold clampE
    rw [min_eq_right (le_trans h (le_trans hlo hhi)), max_eq_left h, EReal.toReal_coe]
  | top => exact absurd hlo (by simp)

/-- The clamp always lands in the interval (for `lo ≤ 0 ≤ hi`). -/
theorem clampE_mem {lo hi : EReal} (hlo : lo ≤ 0) (hhi : 0 ≤ hi) (x : ℝ) :
    clampE lo hi x ∈ boundSet lo hi := by
  by_cases h1 : lo ≤ (x : EReal)
  · by_cases h2 : (x : EReal) ≤ hi
    · rw [clampE_of_mem h1 h2]; exact ⟨h1, h2⟩
    · obtain ⟨r, hr, hc⟩ := clampE_of_hi_le hlo hhi (le_of_lt (not_le.mp h2))
      rw [hc]; exact ⟨by rw [← hr]; exact le_trans hlo hhi, by rw [hr]⟩
  · obtain ⟨r, hr, hc⟩ := clampE_of_le_lo hlo hhi (le_of_lt (not_le.mp h1))
    rw [hc]; exact ⟨by rw [hr], by rw [← hr]; exact le_trans hlo hhi⟩

private theorem exists_between_right {μ : ℝ} {hi : EReal} (h0 : 0 ≤ hi) (h : (μ : EReal) < hi) :
    ∃ y : ℝ, μ < y ∧ (y : EReal) ≤ hi ∧ (μ < 0 → y ≤ 0) := by
  by_cases hμ : μ < 0
  · exact ⟨0, hμ, by simpa using h0, fun _ => le_rfl⟩
  · induction hi using EReal.rec with
    | bot => exact absurd h not_lt_bot
    | coe r => exact ⟨r, EReal.coe_lt_coe_iff.mp h, le_rfl, fun h' => absurd h' hμ⟩
    | top => exact ⟨μ + 1, by linarith, le_top, fun h' => absurd h' hμ⟩

private theorem exists_between_left {μ : ℝ} {lo : EReal} (h0 : lo ≤ 0) (h : lo < (μ : EReal)) :
    ∃ y : ℝ, y < μ ∧ lo ≤ (y : EReal) ∧ (0 < μ → 0 ≤ y) := by
  by_cases hμ : 0 < μ
  · exact ⟨0, hμ, by simpa using h0, fun _ => le_rfl⟩
  · induction lo using EReal.rec with
    | bot => exact ⟨μ - 1, by linarith, bot_le, fun h' => absurd h' hμ⟩
    | coe r => exact ⟨r, EReal.coe_lt_coe_iff.mp h, le_rfl, fun h' => absurd h' hμ⟩
    | top => exact absurd h not_top_lt

/-! ## The subdifferential of `t·|x|` on an interval -/

/-- **Bounded, weighted subdifferential of the absolute value.**
For `t ≥ 0` and `lo ≤ 0 ≤ hi`: `g ∈ ∂(t|·| + ι_[lo,hi])(μ)` iff `μ ∈ [lo, hi]`,
`g` is at most the right slope of `t|·|` at `μ` unless `μ` sits on the upper bound, and
`g` is at least the left slope unless `μ` sits on the lower bound.

Read case by case: `μ = 0` strictly inside gives `|g| ≤ t`; `0 < μ < hi` gives `g = t`;
`μ = hi > 0` gives `g ≥ t`; `lo < μ < 0` gives `g = −t`; `μ = lo < 0` gives `g ≤ −t`.

[formula: abs-subdifferential-bounded] -/
theorem isSubgradientWithinAt_abs_iff {t : ℝ} (ht : 0 ≤ t) {lo hi : EReal} (hlo : lo ≤ 0)
    (hhi : 0 ≤ hi) (g μ : ℝ) :
    IsSubgradientWithinAt (fun x => t * |x|) (boundSet lo hi) g μ ↔
      μ ∈ boundSet lo hi ∧
      ((μ : EReal) < hi → g ≤ if 0 ≤ μ then t else -t) ∧
      (lo < (μ : EReal) → (if 0 < μ then t else -t) ≤ g) := by
  constructor
  · rintro ⟨hμ, hsub⟩
    refine ⟨hμ, fun h => ?_, fun h => ?_⟩
    · obtain ⟨y, hy1, hy2, hy3⟩ := exists_between_right hhi h
      have hyK : y ∈ boundSet lo hi := ⟨le_trans hμ.1 (EReal.coe_le_coe_iff.mpr hy1.le), hy2⟩
      have key := hsub y hyK
      simp only at key
      split_ifs with h0
      · rw [abs_of_nonneg h0, abs_of_nonneg (by linarith)] at key
        by_contra hc
        push Not at hc
        nlinarith
      · push Not at h0
        rw [abs_of_neg h0, abs_of_nonpos (hy3 h0)] at key
        by_contra hc
        push Not at hc
        nlinarith
    · obtain ⟨y, hy1, hy2, hy3⟩ := exists_between_left hlo h
      have hyK : y ∈ boundSet lo hi := ⟨hy2, le_trans (EReal.coe_le_coe_iff.mpr hy1.le) hμ.2⟩
      have key := hsub y hyK
      simp only at key
      split_ifs with h0
      · rw [abs_of_pos h0, abs_of_nonneg (hy3 h0)] at key
        by_contra hc
        push Not at hc
        nlinarith
      · push Not at h0
        rw [abs_of_nonpos h0, abs_of_nonpos (by linarith)] at key
        by_contra hc
        push Not at hc
        nlinarith
  · rintro ⟨hμ, hR, hL⟩
    refine ⟨hμ, fun y hy => ?_⟩
    simp only
    rcases lt_trichotomy y μ with hlt | heq | hgt
    · have hL' := hL (lt_of_le_of_lt hy.1 (EReal.coe_lt_coe_iff.mpr hlt))
      split_ifs at hL' with h0
      · rw [abs_of_pos h0]
        have h1 : y ≤ |y| := le_abs_self y
        nlinarith [mul_nonneg (sub_nonneg.2 hL') (sub_nonneg.2 hlt.le),
          mul_nonneg ht (sub_nonneg.2 h1)]
      · push Not at h0
        rw [abs_of_nonpos h0, abs_of_neg (by linarith)]
        nlinarith [mul_nonneg (sub_nonneg.2 hL') (sub_nonneg.2 hlt.le)]
    · subst heq; simp
    · have hR' := hR (lt_of_lt_of_le (EReal.coe_lt_coe_iff.mpr hgt) hy.2)
      split_ifs at hR' with h0
      · rw [abs_of_nonneg h0, abs_of_pos (by linarith)]
        nlinarith [mul_nonneg (sub_nonneg.2 hR') (sub_nonneg.2 hgt.le)]
      · push Not at h0
        rw [abs_of_neg h0]
        have h1 : -y ≤ |y| := neg_le_abs y
        nlinarith [mul_nonneg (sub_nonneg.2 hR') (sub_nonneg.2 hgt.le),
          mul_nonneg ht (sub_nonneg.2 h1)]


/-- **Regularized, bounded, weighted fixed-point form** (the master statement of this file).
For `t ≥ 0`, `γ > 0`, `lo ≤ 0 ≤ hi`:
`ψ − γμ ∈ ∂(t|·| + ι_[lo,hi])(μ)  ⟺  μ = clamp_[lo,hi](P_t(ψ)/γ)`.

Paper Eq. 21 / System B.4 is the case `lo = ⊥`, `hi = ⊤`; the bound handling ("a simple
projection", paper Sec. 8.2) is derived here: clamp *after* the soft threshold.

[formula: regularized-cone-update] -/
theorem isSubgradientWithinAt_abs_sub_iff_eq_clampE {t γ : ℝ} (ht : 0 ≤ t) (hγ : 0 < γ)
    {lo hi : EReal} (hlo : lo ≤ 0) (hhi : 0 ≤ hi) (ψ μ : ℝ) :
    IsSubgradientWithinAt (fun x => t * |x|) (boundSet lo hi) (ψ - γ * μ) μ ↔
      μ = clampE lo hi (softThreshold t ψ / γ) := by
  rw [isSubgradientWithinAt_abs_iff ht hlo hhi]
  have s1 : t < ψ → softThreshold t ψ = ψ - t := fun h => softThreshold_of_gt ht h
  have s2 : ψ < -t → softThreshold t ψ = ψ + t := fun h => softThreshold_of_lt ht h
  have s3 : -t ≤ ψ → ψ ≤ t → softThreshold t ψ = 0 := fun h1 h2 =>
    softThreshold_of_abs_le (abs_le.mpr ⟨h1, h2⟩)
  generalize softThreshold t ψ = s at s1 s2 s3 ⊢
  constructor
  · rintro ⟨hK, hR, hL⟩
    by_cases hhi' : (μ : EReal) < hi
    · have hR' := hR hhi'
      by_cases hlo' : lo < (μ : EReal)
      · have hL' := hL hlo'
        have hp : s / γ = μ := by
          rw [div_eq_iff hγ.ne']
          rcases lt_trichotomy μ 0 with h | h | h
          · rw [ite_eq_right (not_le.mpr h)] at hR'
            rw [ite_eq_right (not_lt.mpr h.le)] at hL'
            have hψ : ψ < -t := by nlinarith
            rw [s2 hψ]; linarith
          · subst h
            rw [ite_eq_left le_rfl] at hR'
            rw [ite_eq_right (lt_irrefl _)] at hL'
            rw [s3 (by linarith) (by linarith)]; ring
          · rw [ite_eq_left h.le] at hR'
            rw [ite_eq_left h] at hL'
            have hψ : t < ψ := by nlinarith
            rw [s1 hψ]; linarith
        rw [hp, clampE_of_mem hK.1 hK.2]
      · have hlo_eq : lo = (μ : EReal) := le_antisymm hK.1 (not_lt.mp hlo')
        have hμ0 : μ ≤ 0 := by
          have h := hlo; rw [hlo_eq] at h; exact_mod_cast h
        have hle : s / γ ≤ μ := by
          rw [div_le_iff₀ hγ]
          rcases lt_or_eq_of_le hμ0 with h | h
          · rw [ite_eq_right (not_le.mpr h)] at hR'
            have hψ : ψ < -t := by nlinarith
            rw [s2 hψ]; linarith
          · subst h
            rw [ite_eq_left le_rfl] at hR'
            rcases lt_or_ge ψ (-t) with h' | h'
            · rw [s2 h']; linarith
            · rw [s3 h' (by linarith)]; simp
        obtain ⟨r, hr, hc⟩ := clampE_of_le_lo (x := s / γ) hlo hhi
          (by rw [hlo_eq]; exact_mod_cast hle)
        rw [hc]
        rw [hlo_eq] at hr
        exact_mod_cast hr
    · have hhi_eq : hi = (μ : EReal) := le_antisymm (not_lt.mp hhi') hK.2
      have hμ0 : 0 ≤ μ := by
        have h := hhi; rw [hhi_eq] at h; exact_mod_cast h
      by_cases hlo' : lo < (μ : EReal)
      · have hL' := hL hlo'
        have hle : μ ≤ s / γ := by
          rw [le_div_iff₀ hγ]
          rcases lt_or_eq_of_le hμ0 with h | h
          · rw [ite_eq_left h] at hL'
            have hψ : t < ψ := by nlinarith
            rw [s1 hψ]; linarith
          · subst h
            rw [ite_eq_right (lt_irrefl _)] at hL'
            rcases lt_or_ge t ψ with h' | h'
            · rw [s1 h']; linarith
            · rw [s3 (by linarith) h']; simp
        obtain ⟨r, hr, hc⟩ := clampE_of_hi_le (x := s / γ) hlo hhi
          (by rw [hhi_eq]; exact_mod_cast hle)
        rw [hc]
        rw [hhi_eq] at hr
        exact_mod_cast hr
      · have hlo_eq : lo = (μ : EReal) := le_antisymm hK.1 (not_lt.mp hlo')
        have hμ0' : μ ≤ 0 := by
          have h := hlo; rw [hlo_eq] at h; exact_mod_cast h
        have hμ : μ = 0 := le_antisymm hμ0' hμ0
        subst hμ
        rcases le_total (s / γ) 0 with h | h
        · obtain ⟨r, hr, hc⟩ := clampE_of_le_lo (x := s / γ) hlo hhi
            (by rw [hlo_eq]; exact_mod_cast h)
          rw [hc]; rw [hlo_eq] at hr; exact_mod_cast hr
        · obtain ⟨r, hr, hc⟩ := clampE_of_hi_le (x := s / γ) hlo hhi
            (by rw [hhi_eq]; exact_mod_cast h)
          rw [hc]; rw [hhi_eq] at hr; exact_mod_cast hr
  · intro hμ
    by_cases h1 : lo ≤ ((s / γ : ℝ) : EReal)
    · by_cases h2 : ((s / γ : ℝ) : EReal) ≤ hi
      · rw [clampE_of_mem h1 h2] at hμ
        have hs : s = γ * μ := by rw [hμ]; field_simp
        have hpos : 0 < s → 0 < μ := fun h => by rw [hμ]; exact div_pos h hγ
        have hneg : s < 0 → μ < 0 := fun h => by rw [hμ]; exact div_neg_of_neg_of_pos h hγ
        have hzero : s = 0 → μ = 0 := fun h => by rw [hμ, h]; simp
        refine ⟨by rw [hμ]; exact ⟨h1, h2⟩, fun _ => ?_, fun _ => ?_⟩
        · rcases lt_or_ge ψ (-t) with h | h
          · have hs' := s2 h
            have := hneg (by linarith)
            rw [ite_eq_right (not_le.mpr this)]; linarith
          · rcases le_or_gt ψ t with h' | h'
            · have hs' := s3 h h'
              rw [hzero hs', ite_eq_left le_rfl]; linarith
            · have hs' := s1 h'
              have := hpos (by linarith)
              rw [ite_eq_left this.le]; linarith
        · rcases lt_or_ge ψ (-t) with h | h
          · have hs' := s2 h
            have := hneg (by linarith)
            rw [ite_eq_right (not_lt.mpr this.le)]; linarith
          · rcases le_or_gt ψ t with h' | h'
            · have hs' := s3 h h'
              rw [hzero hs', ite_eq_right (lt_irrefl _)]; linarith
            · have hs' := s1 h'
              have := hpos (by linarith)
              rw [ite_eq_left this]; linarith
      · obtain ⟨r, hr, hc⟩ := clampE_of_hi_le (x := s / γ) hlo hhi (le_of_lt (not_le.mp h2))
        rw [hc] at hμ
        subst hμ
        have hμ0 : 0 ≤ μ := by
          have h := hhi; rw [hr] at h; exact_mod_cast h
        have hlt : μ < s / γ := by
          have h := not_le.mp h2; rw [hr] at h; exact_mod_cast h
        rw [lt_div_iff₀ hγ] at hlt
        have hs0 : 0 < s := by nlinarith
        have hψ : t < ψ := by
          by_contra hc'
          push Not at hc'
          rcases lt_or_ge ψ (-t) with h | h
          · rw [s2 h] at hs0; linarith
          · rw [s3 h hc'] at hs0; exact lt_irrefl _ hs0
        rw [s1 hψ] at hlt
        refine ⟨⟨by rw [← hr]; exact le_trans hlo hhi, by rw [hr]⟩, fun h => ?_, fun _ => ?_⟩
        · rw [hr] at h; exact absurd h (lt_irrefl _)
        · split_ifs <;> nlinarith
    · obtain ⟨r, hr, hc⟩ := clampE_of_le_lo (x := s / γ) hlo hhi (le_of_lt (not_le.mp h1))
      rw [hc] at hμ
      subst hμ
      have hμ0 : μ ≤ 0 := by
        have h := hlo; rw [hr] at h; exact_mod_cast h
      have hlt : s / γ < μ := by
        have h := not_le.mp h1; rw [hr] at h; exact_mod_cast h
      rw [div_lt_iff₀ hγ] at hlt
      have hs0 : s < 0 := by nlinarith
      have hψ : ψ < -t := by
        by_contra hc'
        push Not at hc'
        rcases le_or_gt ψ t with h | h
        · rw [s3 hc' h] at hs0; exact lt_irrefl _ hs0
        · rw [s1 h] at hs0; linarith
      rw [s2 hψ] at hlt
      refine ⟨⟨by rw [hr], by rw [← hr]; exact le_trans hlo hhi⟩, fun _ => ?_, fun h => ?_⟩
      · split_ifs <;> nlinarith
      · rw [hr] at h; exact absurd h (lt_irrefl _)

/-- **Bounded, weighted fixed-point form** of the per-vertex inclusion. For `t ≥ 0`,
`lo ≤ 0 ≤ hi`: `φ ∈ ∂(t|·| + ι_[lo,hi])(μ)  ⟺  μ = clamp_[lo,hi](P_t(φ + μ))`.
This is the form of the third optimality condition with cone-angle bounds (paper Sec. 8.2
says only "a simple projection").

[formula: cone-fixed-point-bounded] -/
theorem isSubgradientWithinAt_abs_iff_eq_clampE {t : ℝ} (ht : 0 ≤ t) {lo hi : EReal}
    (hlo : lo ≤ 0) (hhi : 0 ≤ hi) (φ μ : ℝ) :
    IsSubgradientWithinAt (fun x => t * |x|) (boundSet lo hi) φ μ ↔
      μ = clampE lo hi (softThreshold t (φ + μ)) := by
  have h := isSubgradientWithinAt_abs_sub_iff_eq_clampE ht one_pos hlo hhi (φ + μ) μ
  simpa using h

theorem boundSet_bot_top : boundSet ⊥ ⊤ = Set.univ := by
  ext x; simp [boundSet]

/-- **Paper Eq. 19, weighted** (threshold `t = λ·wR ≥ 0`, including `t = 0`):
`μ = P_t(φ + μ)  ⟺  φ ∈ t·∂|μ|`, the subdifferential written out as
`t|μ| + φ·(y − μ) ≤ t|y|` for all `y`.

[formula: soft-threshold-fixed-point] -/
theorem eq_softThreshold_add_iff_isSubgradient {t : ℝ} (ht : 0 ≤ t) (φ μ : ℝ) :
    μ = softThreshold t (φ + μ) ↔
      IsSubgradientWithinAt (fun x => t * |x|) Set.univ φ μ := by
  have h := isSubgradientWithinAt_abs_iff_eq_clampE ht (lo := ⊥) (hi := ⊤) bot_le le_top φ μ
  rw [boundSet_bot_top, clampE_bot_top] at h
  exact h.symm

/-- **Paper Eq. 19, case form**: `μ = P_t(φ + μ)` iff `|φ| ≤ t`, and `φ = t` wherever `μ > 0`,
and `φ = −t` wherever `μ < 0`. (Discrete form of paper Prop. 4.1 / C.1:
`supp μ⁺ ⊆ {φ = +t}`, `supp μ⁻ ⊆ {φ = −t}`.)

[formula: soft-threshold-fixed-point] -/
theorem eq_softThreshold_add_iff_cases {t : ℝ} (ht : 0 ≤ t) (φ μ : ℝ) :
    μ = softThreshold t (φ + μ) ↔ |φ| ≤ t ∧ (0 < μ → φ = t) ∧ (μ < 0 → φ = -t) := by
  rw [eq_softThreshold_add_iff_isSubgradient ht, ← boundSet_bot_top,
    isSubgradientWithinAt_abs_iff ht bot_le le_top]
  simp only [boundSet_bot_top, Set.mem_univ, true_and, EReal.coe_lt_top, EReal.bot_lt_coe,
    forall_const]
  constructor
  · rintro ⟨hR, hL⟩
    rcases lt_trichotomy μ 0 with h | h | h
    · rw [ite_eq_right (not_le.mpr h)] at hR
      rw [ite_eq_right (not_lt.mpr h.le)] at hL
      exact ⟨abs_le.mpr ⟨by linarith, by linarith⟩, fun h' => absurd h (not_lt.mpr h'.le),
        fun _ => le_antisymm hR hL⟩
    · subst h
      rw [ite_eq_left le_rfl] at hR
      rw [ite_eq_right (lt_irrefl _)] at hL
      exact ⟨abs_le.mpr ⟨hL, hR⟩, fun h' => absurd h' (lt_irrefl _),
        fun h' => absurd h' (lt_irrefl _)⟩
    · rw [ite_eq_left h.le] at hR
      rw [ite_eq_left h] at hL
      exact ⟨abs_le.mpr ⟨by linarith, by linarith⟩, fun _ => le_antisymm hR hL,
        fun h' => absurd h (not_lt.mpr h'.le)⟩
  · rintro ⟨habs, hpos, hneg⟩
    obtain ⟨h1, h2⟩ := abs_le.mp habs
    constructor
    · split_ifs with h0
      · exact h2
      · push Not at h0; rw [hneg h0]
    · split_ifs with h0
      · rw [hpos h0]
      · exact h1

/-- **Paper Eq. 19, complementarity form**: `μ = P_t(φ + μ)` iff `|φ| ≤ t` and
`φ·μ = t·|μ|` (so wherever `μ ≠ 0`, `φ` sits on the bound with the sign of `μ`).

[formula: soft-threshold-fixed-point] -/
theorem eq_softThreshold_add_iff {t : ℝ} (ht : 0 ≤ t) (φ μ : ℝ) :
    μ = softThreshold t (φ + μ) ↔ |φ| ≤ t ∧ φ * μ = t * |μ| := by
  rw [eq_softThreshold_add_iff_cases ht]
  constructor
  · rintro ⟨habs, hpos, hneg⟩
    refine ⟨habs, ?_⟩
    rcases lt_trichotomy μ 0 with h | h | h
    · rw [hneg h, abs_of_neg h]; ring
    · subst h; simp
    · rw [hpos h, abs_of_pos h]
  · rintro ⟨habs, hc⟩
    obtain ⟨h1, h2⟩ := abs_le.mp habs
    refine ⟨habs, fun h => ?_, fun h => ?_⟩
    · rw [abs_of_pos h] at hc
      exact mul_right_cancel₀ h.ne' hc
    · rw [abs_of_neg h] at hc
      have : φ * μ = -t * μ := by rw [hc]; ring
      exact mul_right_cancel₀ h.ne this

/-- The same statement in the wording "`|φ| ≤ t`, `μ·(|φ| − t) = 0`, and `φ`, `μ` do not have
opposite signs".

[formula: soft-threshold-fixed-point] -/
theorem eq_softThreshold_add_iff_complementarity {t : ℝ} (ht : 0 ≤ t) (φ μ : ℝ) :
    μ = softThreshold t (φ + μ) ↔ |φ| ≤ t ∧ μ * (|φ| - t) = 0 ∧ 0 ≤ φ * μ := by
  rw [eq_softThreshold_add_iff ht]
  constructor
  · rintro ⟨habs, hc⟩
    refine ⟨habs, ?_, by rw [hc]; exact mul_nonneg ht (abs_nonneg _)⟩
    rcases eq_or_ne μ 0 with h | h
    · rw [h]; simp
    · have h1 : |φ| * |μ| = t * |μ| := by rw [← abs_mul, hc, abs_mul, abs_of_nonneg ht, abs_abs]
      have h2 : |φ| = t := mul_right_cancel₀ (abs_ne_zero.mpr h) h1
      rw [h2]; ring
  · rintro ⟨habs, hc, hs⟩
    refine ⟨habs, ?_⟩
    rcases eq_or_ne μ 0 with h | h
    · rw [h]; simp
    · have h2 : |φ| = t := by
        rcases mul_eq_zero.mp hc with h' | h'
        · exact absurd h' h
        · linarith
      rw [← h2, ← abs_mul, abs_of_nonneg hs]

/-- **Free cone site** (`wR = 0`): the inclusion collapses to `φ = 0`. -/
theorem eq_softThreshold_zero_add_iff (φ μ : ℝ) : μ = softThreshold 0 (φ + μ) ↔ φ = 0 := by
  rw [softThreshold_zero]
  constructor
  · intro h; linarith
  · intro h; rw [h]; ring

/-- **Regularized relation, paper Eq. 21 / System B.3–B.4** (unbounded):
`μ = P_t(φ)/γ  ⟺  φ − γμ ∈ t·∂|μ|`.

In the library `γ` here is `γ_paper / m_i`: the regularizer `(γ/2)‖µ‖²_{L²}` of paper Eq. 12 is
`(γ/2) Σ µ_i²/m_i` for integrated vertex values `µ_i`, which gives `µ_i = (m_i/γ) P_t(φ_i)`.
Paper Eq. 21 prints `µ = P_λ(φ)/γ` without the mass factor.

[formula: regularized-cone-update] -/
theorem eq_softThreshold_div_iff_isSubgradient {t γ : ℝ} (ht : 0 ≤ t) (hγ : 0 < γ) (φ μ : ℝ) :
    μ = softThreshold t φ / γ ↔
      IsSubgradientWithinAt (fun x => t * |x|) Set.univ (φ - γ * μ) μ := by
  have h := isSubgradientWithinAt_abs_sub_iff_eq_clampE ht hγ (lo := ⊥) (hi := ⊤) bot_le le_top φ μ
  rw [boundSet_bot_top, clampE_bot_top] at h
  exact h.symm

/-- **The regularized per-vertex minimization (paper Prop. B.3, with bounds).**
For `t ≥ 0`, `γ > 0`, `lo ≤ 0 ≤ hi`, the function `y ↦ t|y| + (γ/2)y² − ψ·y` has exactly one
minimizer on `[lo, hi]`, namely `clamp_[lo,hi](P_t(ψ)/γ)`.

[formula: regularized-cone-update] -/
theorem minimizes_reg_iff_eq_clampE {t γ : ℝ} (ht : 0 ≤ t) (hγ : 0 < γ)
    {lo hi : EReal} (hlo : lo ≤ 0) (hhi : 0 ≤ hi) (ψ μ : ℝ) :
    (μ ∈ boundSet lo hi ∧ ∀ y ∈ boundSet lo hi,
        t * |μ| + γ / 2 * μ ^ 2 - ψ * μ ≤ t * |y| + γ / 2 * y ^ 2 - ψ * y) ↔
      μ = clampE lo hi (softThreshold t ψ / γ) := by
  have growth : ∀ c : ℝ, c = clampE lo hi (softThreshold t ψ / γ) → ∀ y ∈ boundSet lo hi,
      t * |c| + γ / 2 * c ^ 2 - ψ * c + γ / 2 * (y - c) ^ 2 ≤ t * |y| + γ / 2 * y ^ 2 - ψ * y := by
    intro c hc y hy
    have h := (isSubgradientWithinAt_abs_sub_iff_eq_clampE ht hγ hlo hhi ψ c).mpr hc
    have := h.2 y hy
    simp only at this
    nlinarith
  constructor
  · rintro ⟨hμ, hmin⟩
    set c := clampE lo hi (softThreshold t ψ / γ) with hc
    have hcK : c ∈ boundSet lo hi := clampE_mem hlo hhi _
    have h1 := growth c hc μ hμ
    have h2 := hmin c hcK
    have h3 : γ / 2 * (μ - c) ^ 2 ≤ 0 := by linarith
    have h4 : (μ - c) ^ 2 ≤ 0 := by
      by_contra hcon
      push Not at hcon
      have : 0 < γ / 2 * (μ - c) ^ 2 := by positivity
      linarith
    have h5 : μ - c = 0 := by
      have := sq_nonneg (μ - c)
      exact pow_eq_zero_iff (two_ne_zero) |>.mp (le_antisymm h4 this)
    linarith
  · intro hμ
    refine ⟨by rw [hμ]; exact clampE_mem hlo hhi _, fun y hy => ?_⟩
    have := growth μ hμ y hy
    have h0 : 0 ≤ γ / 2 * (y - μ) ^ 2 := by positivity
    linarith

end DDG
