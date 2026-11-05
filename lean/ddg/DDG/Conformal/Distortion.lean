import DDG.Geometry.Plane

/-!
# Distortion of a triangle map: singular values, quasi-conformal error, log scale

Transcriptions of the formulas of `crates/ddg-bff/src/quality.rs` into definitions on real
numbers (same arrangement of terms), and what they compute; see `DdgProofs.Bff.Quality` for the
plane primitives (`dot2`, `cross2`, `sda`, `dist2`, the reference layout).

## `quasi-conformal-error` — `singular_values`, `quasi_conformal_error`

* `jacobian_maps_triangle` — `ℓ₀₁·J_s = q₁ − q₀` and `x·J_s + y·J_t = q₂ − q₀`;
* `sqStretchMax_add_sqStretchMin` — `S₊ + S₋ = a + c`;
* `sqStretchMax_mul_sqStretchMin` — `S₊·S₋ = ac − b² = det²`;
* `sqStretchMin_nonneg`, `sqStretchMin_le_sqStretchMax` — `0 ≤ S₋ ≤ S₊`;
* `sqStretch_bounds` — `S₋ ≤ ‖e₁J_s + e₂J_t‖² ≤ S₊` for `e₁² + e₂² = 1`;
* `sqStretch_attained` — both bounds are attained, at perpendicular unit vectors;
* `sigmaMin_eq_sqrt` — `σ_min = |det|/σ_max` equals `√S₋` when `S₊ > 0`;
* `sigmaMin_pos_iff` — the test `σ_min > 0` holds iff `det ≠ 0`;
* `one_le_qcError` — `Q ≥ 1` when `σ_min > 0`;
* `qcError_eq_one_iff` — `Q = 1 ↔ a = c ∧ b = 0` when `σ_min > 0`;
* `similarity_iff` — `a = c ∧ b = 0` iff `J` scales every vector by the same factor;
* `qcError_divPair` — `Q` and the test `σ_min > 0` do not change when `J` is divided by `k > 0`;
* `sigmaMax_mul_sigmaMin` — `σ_max·σ_min = |det|`;
* `singularValuesEval_eq` — dividing `J` by its largest entry and scaling the results back
  returns `(σ_max, σ_min, det)`, also for `J = 0`;
* `qcErrorEval_eq` — `quasi_conformal_error` returns `σ_max/σ_min` if `det ≠ 0`, else `+∞`;
* `sda_eq_jacDet_mul_area` — `sda q₀ q₁ q₂ = det·(2·area)`;
* `sda_neg_iff_jacDet_neg` — the image is flipped iff `det < 0`, degenerate iff `det = 0`.


## `face-log-scale-factor` — `log_scale_factor`

* `logScaleFactor_sq_mul` — `½ ln(s²A/A) = ln s`;
* `logScaleFactor_abs_det_mul` — `½ ln(|d|A/A) = ½ ln|d|`;
* `logScaleFactor_eq_half_log_sigma` — `… = ½(ln σ_max + ln σ_min)` when `σ_min > 0`;
* `logScaleFactor_triangle` — the same for a face as the distortion report evaluates it.


Nothing is claimed about floating-point rounding.
-/

noncomputable section

namespace DDG

open Real

/-! ## Singular values of a triangle map -/


/-- First column of the Jacobian, `J_s = (q₁ − q₀)/ℓ₀₁`. -/
def jacS (l01 : ℝ) (q0 q1 : ℝ × ℝ) : ℝ × ℝ := ((q1.1 - q0.1) / l01, (q1.2 - q0.2) / l01)

/-- Second column of the Jacobian, `J_t = ((q₂ − q₀) − x·J_s)/y`. -/
def jacT (x y : ℝ) (js q0 q2 : ℝ × ℝ) : ℝ × ℝ :=
  (((q2.1 - q0.1) - x * js.1) / y, ((q2.2 - q0.2) - x * js.2) / y)

/-- **`J` is the linear part of the triangle map**: the affine map `p ↦ q₀ + p₁·J_s + p₂·J_t`
takes `p₀ = (0,0)`, `p₁ = (ℓ₀₁, 0)`, `p₂ = (x, y)` to `q₀, q₁, q₂`. Only `ℓ₀₁ ≠ 0` and `y ≠ 0`
are needed.

[formula: quasi-conformal-error] -/
theorem jacobian_maps_triangle {l01 y : ℝ} (hl : l01 ≠ 0) (hy : y ≠ 0) (x : ℝ)
    (q0 q1 q2 : ℝ × ℝ) :
    l01 • jacS l01 q0 q1 = q1 - q0 ∧
      x • jacS l01 q0 q1 + y • jacT x y (jacS l01 q0 q1) q0 q2 = q2 - q0 := by
  constructor
  · ext <;> simp only [jacS, Prod.smul_fst, Prod.smul_snd, smul_eq_mul, Prod.fst_sub,
      Prod.snd_sub] <;> field_simp
  · ext <;> simp only [jacS, jacT, Prod.smul_fst, Prod.smul_snd, smul_eq_mul, Prod.fst_sub,
      Prod.snd_sub, Prod.fst_add, Prod.snd_add] <;> field_simp <;> ring

/-- `S₊ = ½(a + c + √((a − c)² + 4b²))` with `a = J_s·J_s`, `b = J_s·J_t`, `c = J_t·J_t`: the
larger eigenvalue of `JᵀJ`, the radicand of `σ_max` as evaluated. -/
def sqStretchMax (s t : ℝ × ℝ) : ℝ :=
  1 / 2 * (dot2 s s + dot2 t t +
    √((dot2 s s - dot2 t t) * (dot2 s s - dot2 t t) + 4 * dot2 s t * dot2 s t))

/-- `S₋ = ½(a + c − √((a − c)² + 4b²))`: the smaller eigenvalue of `JᵀJ`. (Not evaluated by the
Rust code, which uses `σ_min = |det|/σ_max`, see `sigmaMin_eq_sqrt`.) -/
def sqStretchMin (s t : ℝ × ℝ) : ℝ :=
  1 / 2 * (dot2 s s + dot2 t t -
    √((dot2 s s - dot2 t t) * (dot2 s s - dot2 t t) + 4 * dot2 s t * dot2 s t))

/-- `σ_max = √(½(a + c + √((a − c)² + 4b²)))`. -/
def sigmaMax (s t : ℝ × ℝ) : ℝ := √(sqStretchMax s t)

/-- `σ_min = |det J|/σ_max` if `σ_max > 0`, else `0`; `det J = J_s × J_t`. -/
def sigmaMin (s t : ℝ × ℝ) : ℝ :=
  if 0 < sigmaMax s t then |cross2 s t| / sigmaMax s t else 0

/-- The quasi-conformal error `Q = σ_max/σ_min`. (The Rust code returns `+∞` unless
`σ_min > 0`; that case is `cross2 s t = 0`, see `sigmaMin_pos_iff`.) -/
def qcError (s t : ℝ × ℝ) : ℝ := sigmaMax s t / sigmaMin s t

/-- Lagrange's identity: `a·c − b² = det²`. -/
theorem dot2_mul_dot2_sub_sq (s t : ℝ × ℝ) :
    dot2 s s * dot2 t t - dot2 s t ^ 2 = cross2 s t ^ 2 := by
  unfold dot2 cross2
  ring

theorem dot2_self_nonneg (s : ℝ × ℝ) : 0 ≤ dot2 s s := by
  unfold dot2
  nlinarith [mul_self_nonneg s.1, mul_self_nonneg s.2]

private theorem disc_nonneg (a b c : ℝ) : 0 ≤ (a - c) * (a - c) + 4 * b * b := by
  nlinarith [mul_self_nonneg (a - c), mul_self_nonneg b]

/-- `S₊ + S₋ = a + c` (the trace of `JᵀJ`).

[formula: quasi-conformal-error] -/
theorem sqStretchMax_add_sqStretchMin (s t : ℝ × ℝ) :
    sqStretchMax s t + sqStretchMin s t = dot2 s s + dot2 t t := by
  unfold sqStretchMax sqStretchMin
  ring

/-- `S₊ · S₋ = a·c − b² = det²` (the determinant of `JᵀJ`).

[formula: quasi-conformal-error] -/
theorem sqStretchMax_mul_sqStretchMin (s t : ℝ × ℝ) :
    sqStretchMax s t * sqStretchMin s t = dot2 s s * dot2 t t - dot2 s t ^ 2 ∧
      sqStretchMax s t * sqStretchMin s t = cross2 s t ^ 2 := by
  have h := sq_sqrt (disc_nonneg (dot2 s s) (dot2 s t) (dot2 t t))
  have h1 : sqStretchMax s t * sqStretchMin s t = dot2 s s * dot2 t t - dot2 s t ^ 2 := by
    unfold sqStretchMax sqStretchMin
    linear_combination (-1 / 4 : ℝ) * h
  exact ⟨h1, h1.trans (dot2_mul_dot2_sub_sq s t)⟩

/-- `S₋ ≤ S₊`.

[formula: quasi-conformal-error] -/
theorem sqStretchMin_le_sqStretchMax (s t : ℝ × ℝ) : sqStretchMin s t ≤ sqStretchMax s t := by
  unfold sqStretchMax sqStretchMin
  linarith [sqrt_nonneg ((dot2 s s - dot2 t t) * (dot2 s s - dot2 t t) + 4 * dot2 s t * dot2 s t)]

/-- `0 ≤ S₋`.

[formula: quasi-conformal-error] -/
theorem sqStretchMin_nonneg (s t : ℝ × ℝ) : 0 ≤ sqStretchMin s t := by
  have ha := dot2_self_nonneg s
  have hc := dot2_self_nonneg t
  have hdet : 0 ≤ dot2 s s * dot2 t t - dot2 s t ^ 2 := by
    rw [dot2_mul_dot2_sub_sq]
    positivity
  have hD : √((dot2 s s - dot2 t t) * (dot2 s s - dot2 t t) + 4 * dot2 s t * dot2 s t) ≤
      dot2 s s + dot2 t t := by
    calc _ ≤ √((dot2 s s + dot2 t t) ^ 2) := sqrt_le_sqrt (by nlinarith)
      _ = _ := sqrt_sq (by linarith)
  unfold sqStretchMin
  linarith

theorem sqStretchMax_nonneg (s t : ℝ × ℝ) : 0 ≤ sqStretchMax s t :=
  (sqStretchMin_nonneg s t).trans (sqStretchMin_le_sqStretchMax s t)

private theorem quadForm_bounds (a b c : ℝ) {e₁ e₂ : ℝ} (he : e₁ ^ 2 + e₂ ^ 2 = 1) :
    1 / 2 * (a + c - √((a - c) * (a - c) + 4 * b * b)) ≤
        a * e₁ ^ 2 + 2 * b * (e₁ * e₂) + c * e₂ ^ 2 ∧
      a * e₁ ^ 2 + 2 * b * (e₁ * e₂) + c * e₂ ^ 2 ≤
        1 / 2 * (a + c + √((a - c) * (a - c) + 4 * b * b)) := by
  have hpq : (e₁ ^ 2 - e₂ ^ 2) ^ 2 + (2 * e₁ * e₂) ^ 2 = 1 := by
    linear_combination (e₁ ^ 2 + e₂ ^ 2 + 1) * he
  have hid : (a - c) * (a - c) + 4 * b * b -
      ((a - c) * (e₁ ^ 2 - e₂ ^ 2) + 2 * b * (2 * e₁ * e₂)) ^ 2 =
      ((a - c) * (2 * e₁ * e₂) - 2 * b * (e₁ ^ 2 - e₂ ^ 2)) ^ 2 := by
    linear_combination (-((a - c) * (a - c) + 4 * b * b)) * hpq
  have hD : |(a - c) * (e₁ ^ 2 - e₂ ^ 2) + 2 * b * (2 * e₁ * e₂)| ≤
      √((a - c) * (a - c) + 4 * b * b) :=
    abs_le_sqrt (by nlinarith [sq_nonneg ((a - c) * (2 * e₁ * e₂) - 2 * b * (e₁ ^ 2 - e₂ ^ 2))])
  have hval : a * e₁ ^ 2 + 2 * b * (e₁ * e₂) + c * e₂ ^ 2 =
      1 / 2 * (a + c + ((a - c) * (e₁ ^ 2 - e₂ ^ 2) + 2 * b * (2 * e₁ * e₂))) := by
    linear_combination (1 / 2 * (a + c)) * he
  rw [hval]
  have := abs_le.mp hD
  constructor <;> linarith [this.1, this.2]

private theorem dot2_comb (s t : ℝ × ℝ) (e₁ e₂ : ℝ) :
    dot2 (e₁ • s + e₂ • t) (e₁ • s + e₂ • t) =
      dot2 s s * e₁ ^ 2 + 2 * dot2 s t * (e₁ * e₂) + dot2 t t * e₂ ^ 2 := by
  simp only [dot2, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  ring

/-- **`√S₊` and `√S₋` are the largest and the smallest stretch of `J`.** For every unit vector
`e = (e₁, e₂)`, `S₋ ≤ ‖e₁ J_s + e₂ J_t‖² ≤ S₊`.

[formula: quasi-conformal-error] -/
theorem sqStretch_bounds (s t : ℝ × ℝ) {e₁ e₂ : ℝ} (he : e₁ ^ 2 + e₂ ^ 2 = 1) :
    sqStretchMin s t ≤ dot2 (e₁ • s + e₂ • t) (e₁ • s + e₂ • t) ∧
      dot2 (e₁ • s + e₂ • t) (e₁ • s + e₂ • t) ≤ sqStretchMax s t := by
  rw [dot2_comb]
  exact quadForm_bounds (dot2 s s) (dot2 s t) (dot2 t t) he

/-- **Both bounds are attained**, at perpendicular unit vectors `e = (e₁, e₂)` and
`e⊥ = (−e₂, e₁)`: `‖J e‖² = S₊` and `‖J e⊥‖² = S₋`.

[formula: quasi-conformal-error] -/
theorem sqStretch_attained (s t : ℝ × ℝ) :
    ∃ e₁ e₂ : ℝ, e₁ ^ 2 + e₂ ^ 2 = 1 ∧
      dot2 (e₁ • s + e₂ • t) (e₁ • s + e₂ • t) = sqStretchMax s t ∧
      dot2 ((-e₂) • s + e₁ • t) ((-e₂) • s + e₁ • t) = sqStretchMin s t := by
  simp only [dot2_comb]
  unfold sqStretchMax sqStretchMin
  generalize dot2 s s = a
  generalize dot2 s t = b
  generalize dot2 t t = c
  have h4 := sq_sqrt (disc_nonneg a b c)
  rcases eq_or_ne (⟨a - c, 2 * b⟩ : ℂ) 0 with hz | hz
  · have hac : a - c = 0 := congrArg Complex.re hz
    have hb : 2 * b = 0 := congrArg Complex.im hz
    have hb' : b = 0 := by linarith
    have hac' : a = c := by linarith
    refine ⟨1, 0, by norm_num, ?_, ?_⟩ <;> simp [hb', hac'] <;> ring
  · set z : ℂ := ⟨a - c, 2 * b⟩ with hzdef
    set D := √((a - c) * (a - c) + 4 * b * b) with hDdef
    have hn : ‖z‖ = D := by
      rw [Complex.norm_eq_sqrt_sq_add_sq, hDdef]
      congr 1
      simp only [hzdef]
      ring
    have hDpos : 0 < D := by
      rw [← hn]
      exact norm_pos_iff.mpr hz
    have hc : cos (Complex.arg z) = (a - c) / D := by
      rw [Complex.cos_arg hz, hn]
    have hs : sin (Complex.arg z) = 2 * b / D := by
      rw [Complex.sin_arg, hn]
    set θ := Complex.arg z / 2 with hθ
    have h2θ : Complex.arg z = 2 * θ := by rw [hθ]; ring
    rw [h2θ, cos_two_mul] at hc
    rw [h2θ, sin_two_mul] at hs
    have h1' : (2 * cos θ ^ 2 - 1) * D = a - c := by
      rw [hc]
      field_simp
    have h2' : 2 * sin θ * cos θ * D = 2 * b := by
      rw [hs]
      field_simp
    have h3 : cos θ ^ 2 + sin θ ^ 2 = 1 := cos_sq_add_sin_sq θ
    refine ⟨cos θ, sin θ, h3, ?_, ?_⟩
    · apply mul_right_cancel₀ hDpos.ne'
      linear_combination (1 / 2 * (a - c)) * h1' + b * h2' + (c * D) * h3 - (1 / 2) * h4
    · apply mul_right_cancel₀ hDpos.ne'
      linear_combination (-(1 / 2) * (a - c)) * h1' - b * h2' + (a * D) * h3 + (1 / 2) * h4

theorem sigmaMax_nonneg (s t : ℝ × ℝ) : 0 ≤ sigmaMax s t := sqrt_nonneg _

/-- `σ_max · √S₋ = |det J|`. -/
theorem sigmaMax_mul_sqrt_sqStretchMin (s t : ℝ × ℝ) :
    sigmaMax s t * √(sqStretchMin s t) = |cross2 s t| := by
  unfold sigmaMax
  rw [← sqrt_mul (sqStretchMax_nonneg s t), (sqStretchMax_mul_sqStretchMin s t).2,
    sqrt_sq_eq_abs]

/-- **`σ_min` as evaluated is `√S₋`**: for `S₊ > 0` (i.e. `J ≠ 0`),
`|det J|/σ_max = √(½(a + c − √((a − c)² + 4b²)))`.

[formula: quasi-conformal-error] -/
theorem sigmaMin_eq_sqrt {s t : ℝ × ℝ} (h : 0 < sqStretchMax s t) :
    sigmaMin s t = √(sqStretchMin s t) := by
  have hmax : 0 < sigmaMax s t := sqrt_pos.mpr h
  unfold sigmaMin
  rw [ite_eq_left hmax, div_eq_iff hmax.ne', mul_comm, sigmaMax_mul_sqrt_sqStretchMin]

/-- The test `σ_min > 0` of `quasi_conformal_error` holds exactly when `det J ≠ 0` (the image
triangle is not degenerate).

[formula: quasi-conformal-error] -/
theorem sigmaMin_pos_iff (s t : ℝ × ℝ) : 0 < sigmaMin s t ↔ cross2 s t ≠ 0 := by
  unfold sigmaMin
  constructor
  · intro h hdet
    by_cases hmax : 0 < sigmaMax s t
    · rw [ite_eq_left hmax, hdet, abs_zero, zero_div] at h
      exact lt_irrefl _ h
    · rw [ite_eq_right hmax] at h
      exact lt_irrefl _ h
  · intro hdet
    have hprod := sigmaMax_mul_sqrt_sqStretchMin s t
    have habs : 0 < |cross2 s t| := abs_pos.mpr hdet
    have hmax : 0 < sigmaMax s t := by
      rcases (sigmaMax_nonneg s t).eq_or_lt with h0 | h0
      · rw [← h0, zero_mul] at hprod
        linarith
      · exact h0
    rw [ite_eq_left hmax]
    exact div_pos habs hmax

/-- If `σ_min > 0` then `σ_max > 0` and `σ_min = √S₋ ≤ σ_max`. -/
theorem sigmaMin_le_sigmaMax {s t : ℝ × ℝ} (h : 0 < sigmaMin s t) :
    0 < sigmaMax s t ∧ sigmaMin s t = √(sqStretchMin s t) ∧ sigmaMin s t ≤ sigmaMax s t := by
  have hmax : 0 < sigmaMax s t := by
    by_contra hmax
    unfold sigmaMin at h
    rw [ite_eq_right hmax] at h
    exact lt_irrefl _ h
  have hS : 0 < sqStretchMax s t := sqrt_pos.mp hmax
  have e := sigmaMin_eq_sqrt hS
  exact ⟨hmax, e, e ▸ sqrt_le_sqrt (sqStretchMin_le_sqStretchMax s t)⟩

/-- **`Q ≥ 1`** whenever `σ_min > 0`.

[formula: quasi-conformal-error] -/
theorem one_le_qcError {s t : ℝ × ℝ} (h : 0 < sigmaMin s t) : 1 ≤ qcError s t := by
  unfold qcError
  rw [one_le_div h]
  exact (sigmaMin_le_sigmaMax h).2.2

/-- **`Q = 1` exactly for a similarity**: for `σ_min > 0`, `σ_max/σ_min = 1` iff `a = c` and
`b = 0`, i.e. the columns of `J` have the same length and are perpendicular.

[formula: quasi-conformal-error] -/
theorem qcError_eq_one_iff {s t : ℝ × ℝ} (h : 0 < sigmaMin s t) :
    qcError s t = 1 ↔ dot2 s s = dot2 t t ∧ dot2 s t = 0 := by
  obtain ⟨-, e, -⟩ := sigmaMin_le_sigmaMax h
  have hd := disc_nonneg (dot2 s s) (dot2 s t) (dot2 t t)
  unfold qcError
  rw [div_eq_one_iff_eq h.ne', e, sigmaMax,
    sqrt_inj (sqStretchMax_nonneg s t) (sqStretchMin_nonneg s t)]
  unfold sqStretchMax sqStretchMin
  constructor
  · intro heq
    have h0 : √((dot2 s s - dot2 t t) * (dot2 s s - dot2 t t) + 4 * dot2 s t * dot2 s t) = 0 := by
      linarith
    rw [sqrt_eq_zero hd] at h0
    have h1 : dot2 s s - dot2 t t = 0 := by
      nlinarith [mul_self_nonneg (dot2 s s - dot2 t t), mul_self_nonneg (dot2 s t)]
    have h2 : dot2 s t = 0 := by
      nlinarith [mul_self_nonneg (dot2 s s - dot2 t t), mul_self_nonneg (dot2 s t)]
    exact ⟨by linarith, h2⟩
  · rintro ⟨h1, h2⟩
    rw [h1, h2]
    simp

/-- **What `a = c ∧ b = 0` means**: the columns of `J` have the same length and are
perpendicular iff `J` scales every vector by the same factor,
`‖e₁ J_s + e₂ J_t‖² = a·(e₁² + e₂²)` for all `e`; that is, the linear map is a similarity (or
zero).

[formula: quasi-conformal-error] -/
theorem similarity_iff (s t : ℝ × ℝ) :
    (dot2 s s = dot2 t t ∧ dot2 s t = 0) ↔
      ∀ e₁ e₂ : ℝ, dot2 (e₁ • s + e₂ • t) (e₁ • s + e₂ • t) = dot2 s s * (e₁ ^ 2 + e₂ ^ 2) := by
  constructor
  · rintro ⟨h1, h2⟩ e₁ e₂
    rw [dot2_comb, ← h1, h2]
    ring
  · intro h
    have h01 := h 0 1
    have h11 := h 1 1
    rw [dot2_comb] at h01 h11
    constructor <;> linarith

/-- `σ_max · σ_min = |det J|` (also when `J = 0`).

[formula: quasi-conformal-error] -/
theorem sigmaMax_mul_sigmaMin (s t : ℝ × ℝ) : sigmaMax s t * sigmaMin s t = |cross2 s t| := by
  by_cases hmax : 0 < sigmaMax s t
  · unfold sigmaMin
    rw [ite_eq_left hmax, mul_div_cancel₀ _ hmax.ne']
  · have h0 : sigmaMax s t = 0 := le_antisymm (not_lt.mp hmax) (sigmaMax_nonneg s t)
    have h := sigmaMax_mul_sqrt_sqStretchMin s t
    rw [h0, zero_mul] at h ⊢
    exact h

/-! ### Dividing the Jacobian by its largest entry -/

/-- A plane vector divided by a number, entry by entry (`v.map(|x| x / k)`). -/
def divPair (k : ℝ) (v : ℝ × ℝ) : ℝ × ℝ := (v.1 / k, v.2 / k)

/-- The largest entry of the Jacobian in absolute value, as evaluated:
`|J_s[0]|.max(|J_s[1]|).max(|J_t[0]|).max(|J_t[1]|)`. -/
def jacScale (s t : ℝ × ℝ) : ℝ := max (max (max |s.1| |s.2|) |t.1|) |t.2|

/-- `singular_values` as evaluated: with `k = jacScale J`, the triple `(0, 0, 0)` if `k = 0`, else
`(σ_max(J/k)·k, σ_min(J/k)·k, det(J/k)·k·k)`. -/
def singularValuesEval (s t : ℝ × ℝ) : ℝ × ℝ × ℝ :=
  if jacScale s t = 0 then (0, 0, 0)
  else
    (sigmaMax (divPair (jacScale s t) s) (divPair (jacScale s t) t) * jacScale s t,
      sigmaMin (divPair (jacScale s t) s) (divPair (jacScale s t) t) * jacScale s t,
      cross2 (divPair (jacScale s t) s) (divPair (jacScale s t) t) * jacScale s t * jacScale s t)

theorem dot2_divPair (k : ℝ) (u v : ℝ × ℝ) :
    dot2 (divPair k u) (divPair k v) = dot2 u v / k ^ 2 := by
  unfold dot2 divPair
  ring

theorem cross2_divPair (k : ℝ) (u v : ℝ × ℝ) :
    cross2 (divPair k u) (divPair k v) = cross2 u v / k ^ 2 := by
  unfold cross2 divPair
  ring

theorem sqStretchMax_divPair (k : ℝ) (s t : ℝ × ℝ) :
    sqStretchMax (divPair k s) (divPair k t) = sqStretchMax s t / k ^ 2 := by
  unfold sqStretchMax
  simp only [dot2_divPair]
  have h : (dot2 s s / k ^ 2 - dot2 t t / k ^ 2) * (dot2 s s / k ^ 2 - dot2 t t / k ^ 2) +
      4 * (dot2 s t / k ^ 2) * (dot2 s t / k ^ 2) =
      ((dot2 s s - dot2 t t) * (dot2 s s - dot2 t t) + 4 * dot2 s t * dot2 s t) / (k ^ 2) ^ 2 := by
    ring
  rw [h, sqrt_div (disc_nonneg _ _ _), sqrt_sq (sq_nonneg k)]
  ring

/-- `σ_max(J/k) = σ_max(J)/k` for `k > 0`. -/
theorem sigmaMax_divPair {k : ℝ} (hk : 0 < k) (s t : ℝ × ℝ) :
    sigmaMax (divPair k s) (divPair k t) = sigmaMax s t / k := by
  unfold sigmaMax
  rw [sqStretchMax_divPair, sqrt_div (sqStretchMax_nonneg s t), sqrt_sq hk.le]

/-- `σ_min(J/k) = σ_min(J)/k` for `k > 0`. -/
theorem sigmaMin_divPair {k : ℝ} (hk : 0 < k) (s t : ℝ × ℝ) :
    sigmaMin (divPair k s) (divPair k t) = sigmaMin s t / k := by
  unfold sigmaMin
  rw [sigmaMax_divPair hk, cross2_divPair]
  by_cases hmax : 0 < sigmaMax s t
  · rw [ite_eq_left hmax, ite_eq_left (div_pos hmax hk), abs_div, abs_of_pos (pow_pos hk 2)]
    field_simp
  · have h0 : sigmaMax s t = 0 := le_antisymm (not_lt.mp hmax) (sigmaMax_nonneg s t)
    simp [h0]

/-- **The quasi-conformal error does not depend on the scale of the Jacobian**: for `k > 0`,
`Q(J/k) = Q(J)`. (`quasi_conformal_error` evaluates `Q` for the Jacobian divided by its largest
entry.)

[formula: quasi-conformal-error] -/
theorem qcError_divPair {k : ℝ} (hk : 0 < k) (s t : ℝ × ℝ) :
    qcError (divPair k s) (divPair k t) = qcError s t ∧
      (0 < sigmaMin (divPair k s) (divPair k t) ↔ 0 < sigmaMin s t) := by
  unfold qcError
  rw [sigmaMax_divPair hk, sigmaMin_divPair hk]
  exact ⟨div_div_div_cancel_right₀ hk.ne' _ _, div_pos_iff_of_pos_right hk⟩

theorem jacScale_nonneg (s t : ℝ × ℝ) : 0 ≤ jacScale s t :=
  le_max_of_le_right (abs_nonneg _)

/-- The largest entry vanishes exactly for the zero matrix. -/
theorem jacScale_eq_zero_iff (s t : ℝ × ℝ) : jacScale s t = 0 ↔ s = 0 ∧ t = 0 := by
  constructor
  · intro h
    have h1 : |s.1| ≤ jacScale s t :=
      le_max_of_le_left (le_max_of_le_left (le_max_left _ _))
    have h2 : |s.2| ≤ jacScale s t :=
      le_max_of_le_left (le_max_of_le_left (le_max_right _ _))
    have h3 : |t.1| ≤ jacScale s t := le_max_of_le_left (le_max_right _ _)
    have h4 : |t.2| ≤ jacScale s t := le_max_right _ _
    rw [h] at h1 h2 h3 h4
    exact ⟨Prod.ext (abs_nonpos_iff.mp h1) (abs_nonpos_iff.mp h2),
      Prod.ext (abs_nonpos_iff.mp h3) (abs_nonpos_iff.mp h4)⟩
  · rintro ⟨rfl, rfl⟩
    simp [jacScale]

/-- **`singular_values` returns `(σ_max, σ_min, det J)`**: dividing the Jacobian by its largest
entry and scaling the results back changes nothing in exact arithmetic, and the early return
`(0, 0, 0)` for the zero matrix agrees with the formulas.

[formula: quasi-conformal-error] -/
theorem singularValuesEval_eq (s t : ℝ × ℝ) :
    singularValuesEval s t = (sigmaMax s t, sigmaMin s t, cross2 s t) := by
  unfold singularValuesEval
  by_cases hk : jacScale s t = 0
  · obtain ⟨rfl, rfl⟩ := (jacScale_eq_zero_iff s t).mp hk
    rw [ite_eq_left hk]
    simp [sigmaMax, sigmaMin, sqStretchMax, dot2, cross2]
  · have hpos : 0 < jacScale s t := lt_of_le_of_ne (jacScale_nonneg s t) (Ne.symm hk)
    rw [ite_eq_right hk, sigmaMax_divPair hpos, sigmaMin_divPair hpos, cross2_divPair]
    refine Prod.ext ?_ (Prod.ext ?_ ?_) <;> simp only <;> field_simp

/-- `quasi_conformal_error` as evaluated, with values in the extended reals: with
`k = jacScale J`, the result is `+∞` if `k = 0` or if `σ_min(J/k) > 0` fails, and
`σ_max(J/k)/σ_min(J/k)` otherwise. -/
def qcErrorEval (s t : ℝ × ℝ) : EReal :=
  if jacScale s t = 0 then ⊤
  else if 0 < sigmaMin (divPair (jacScale s t) s) (divPair (jacScale s t) t) then
    ((sigmaMax (divPair (jacScale s t) s) (divPair (jacScale s t) t) /
      sigmaMin (divPair (jacScale s t) s) (divPair (jacScale s t) t) : ℝ) : EReal)
  else ⊤

/-- **`quasi_conformal_error` returns `Q = σ_max/σ_min` if `det J ≠ 0`, and `+∞` otherwise**
(a degenerate image, including `J = 0`).

[formula: quasi-conformal-error] -/
theorem qcErrorEval_eq (s t : ℝ × ℝ) :
    qcErrorEval s t = if cross2 s t ≠ 0 then ((qcError s t : ℝ) : EReal) else ⊤ := by
  unfold qcErrorEval
  by_cases hk : jacScale s t = 0
  · obtain ⟨rfl, rfl⟩ := (jacScale_eq_zero_iff s t).mp hk
    rw [ite_eq_left hk, ite_eq_right (by simp [cross2])]
  · have hpos : 0 < jacScale s t := lt_of_le_of_ne (jacScale_nonneg s t) (Ne.symm hk)
    obtain ⟨hq, hiff⟩ := qcError_divPair hpos s t
    rw [ite_eq_right hk]
    by_cases hdet : cross2 s t ≠ 0
    · rw [ite_eq_left hdet, ite_eq_left (hiff.mpr ((sigmaMin_pos_iff s t).mpr hdet)), ← hq]
      rfl
    · rw [ite_eq_right hdet,
        ite_eq_right fun h => hdet ((sigmaMin_pos_iff s t).mp (hiff.mp h))]

/-! ### The determinant is the ratio of the signed areas -/

/-- **`det J = (signed area in the plane)/(area on the surface)`**:
`sda q₀ q₁ q₂ = det J · (2·area)`, with `J_s`, `J_t` as evaluated. Only `ℓ₀₁ ≠ 0` and
`area ≠ 0` are needed, and `x` is arbitrary (it drops out of `J_s × J_t`); in the Rust code
`x = layoutX ℓ₀₁ ℓ₁₂ ℓ₂₀`.

[formula: quasi-conformal-error] -/
theorem sda_eq_jacDet_mul_area {l01 area : ℝ} (hl : l01 ≠ 0) (hA : area ≠ 0) (x : ℝ)
    (q0 q1 q2 : ℝ × ℝ) :
    sda q0 q1 q2 =
      cross2 (jacS l01 q0 q1) (jacT x (layoutY l01 area) (jacS l01 q0 q1) q0 q2) *
        (2 * area) := by
  unfold sda cross2 jacS jacT layoutY
  field_simp
  ring

/-- **A face is flipped iff `det J < 0`**, and its image is degenerate iff `det J = 0`: for
`ℓ₀₁ ≠ 0` and `area > 0`, the signed area of the image has the sign of the determinant.

[formula: quasi-conformal-error] -/
theorem sda_neg_iff_jacDet_neg {l01 area : ℝ} (hl : l01 ≠ 0) (hA : 0 < area) (x : ℝ)
    (q0 q1 q2 : ℝ × ℝ) :
    (sda q0 q1 q2 < 0 ↔
        cross2 (jacS l01 q0 q1) (jacT x (layoutY l01 area) (jacS l01 q0 q1) q0 q2) < 0) ∧
      (sda q0 q1 q2 = 0 ↔
        cross2 (jacS l01 q0 q1) (jacT x (layoutY l01 area) (jacS l01 q0 q1) q0 q2) = 0) := by
  rw [sda_eq_jacDet_mul_area hl hA.ne' x]
  have h2 : 0 < 2 * area := by positivity
  refine ⟨⟨fun h => ?_, fun h => mul_neg_of_neg_of_pos h h2⟩, by simp [hA.ne']⟩
  by_contra hc
  exact absurd h (not_lt.mpr (mul_nonneg (not_lt.mp hc) h2.le))

/-! ## Log scale factor -/

/-- `log_scale_factor`: `u = ½ ln(B/A)` for a triangle of area `A` on the surface and `B` in the
plane. -/
def logScaleFactor (surfaceArea planeArea : ℝ) : ℝ := 1 / 2 * log (planeArea / surfaceArea)

/-- **The log scale factor of a similarity is the logarithm of its scale**: if lengths grow by
the factor `s > 0`, areas grow by `s²` and `u = ½ ln(s²A/A) = ln s`.

[formula: face-log-scale-factor] -/
theorem logScaleFactor_sq_mul {A s : ℝ} (hA : 0 < A) (hs : 0 < s) :
    logScaleFactor A (s ^ 2 * A) = log s := by
  unfold logScaleFactor
  rw [mul_div_assoc, div_self hA.ne', mul_one, sq, log_mul hs.ne' hs.ne']
  ring

/-- If the plane area is `|det|·A` then `u = ½ ln|det|`.

[formula: face-log-scale-factor] -/
theorem logScaleFactor_abs_det_mul {A : ℝ} (hA : 0 < A) (d : ℝ) :
    logScaleFactor A (|d| * A) = 1 / 2 * log |d| := by
  unfold logScaleFactor
  rw [mul_div_assoc, div_self hA.ne', mul_one]

/-- **`u = ½(ln σ_max + ln σ_min)`**: for a map with Jacobian columns `s, t` and `det ≠ 0`
(equivalently `σ_min > 0`, which implies `σ_max > 0`), the log scale factor of the areas `A`
and `|det|·A` is the mean of the logarithms of the two singular values.

[formula: face-log-scale-factor] -/
theorem logScaleFactor_eq_half_log_sigma {A : ℝ} (hA : 0 < A) {s t : ℝ × ℝ}
    (h : 0 < sigmaMin s t) :
    logScaleFactor A (|cross2 s t| * A) = 1 / 2 * (log (sigmaMax s t) + log (sigmaMin s t)) := by
  rw [logScaleFactor_abs_det_mul hA, ← sigmaMax_mul_sigmaMin,
    log_mul (sigmaMin_le_sigmaMax h).1.ne' h.ne']

/-- **The log scale factor of a face as the distortion report evaluates it**
(`log_scale_factor(area, |½·sda q₀ q₁ q₂|)`) is `½ ln|det J| = ½(ln σ_max + ln σ_min)` of the
Jacobian of the face, for `ℓ₀₁ ≠ 0`, `area > 0` and a nondegenerate image.

[formula: face-log-scale-factor] -/
theorem logScaleFactor_triangle {l01 area : ℝ} (hl : l01 ≠ 0) (hA : 0 < area) (x : ℝ)
    {q0 q1 q2 : ℝ × ℝ} (hq : sda q0 q1 q2 ≠ 0) :
    logScaleFactor area |1 / 2 * sda q0 q1 q2| =
        1 / 2 * log |cross2 (jacS l01 q0 q1) (jacT x (layoutY l01 area) (jacS l01 q0 q1) q0 q2)| ∧
      logScaleFactor area |1 / 2 * sda q0 q1 q2| =
        1 / 2 * (log (sigmaMax (jacS l01 q0 q1)
            (jacT x (layoutY l01 area) (jacS l01 q0 q1) q0 q2)) +
          log (sigmaMin (jacS l01 q0 q1) (jacT x (layoutY l01 area) (jacS l01 q0 q1) q0 q2))) := by
  have hdet := sda_eq_jacDet_mul_area hl hA.ne' x q0 q1 q2
  have hne : cross2 (jacS l01 q0 q1) (jacT x (layoutY l01 area) (jacS l01 q0 q1) q0 q2) ≠ 0 := by
    intro h0
    rw [h0, zero_mul] at hdet
    exact hq hdet
  have hB : |1 / 2 * sda q0 q1 q2| =
      |cross2 (jacS l01 q0 q1) (jacT x (layoutY l01 area) (jacS l01 q0 q1) q0 q2)| * area := by
    rw [hdet, ← abs_of_pos hA, ← abs_mul, abs_of_pos hA]
    congr 1
    ring
  rw [hB]
  exact ⟨logScaleFactor_abs_det_mul hA _,
    logScaleFactor_eq_half_log_sigma hA ((sigmaMin_pos_iff _ _).mpr hne)⟩

end DDG
