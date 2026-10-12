import DdgProofs.Evaluated

/-!
# Quality measures of a flattening, in the form the Rust code evaluates them

Transcriptions of the formulas of `crates/ddg-bff/src/quality.rs` and
`crates/ddg-bff/src/mollify.rs` into definitions on real numbers (same arrangement of terms),
and what they compute. Plane points are pairs `ℝ × ℝ`; lengths are written `√(Δx² + Δy²)`
(`dist2`), never with the metric of the type `ℝ × ℝ` (which is the sup metric). `dot2 u v` and
`cross2 u v` are `u·v` and `u×v`; `toC` turns a point into a complex number.

## `signed-area-2d` — `signed_double_area` (`sda`)

* `sda_swap_left`, `sda_swap_right`, `sda_swap_outer` — swapping two corners changes the sign;
* `sda_rotate` — cyclic rotation of the corners does not change it;
* `sda_translate` — translation does not change it;
* `sda_sq_mul_four` — `sda²·4 = 2xy + 2yz + 2zx − x² − y² − z²` in the squared side lengths;
* `abs_sda_eq_two_mul_triArea` — `|sda| = 2·triArea` of the Euclidean side lengths (no
  nondegeneracy hypothesis).

## `quasi-conformal-error` — `singular_values`, `quasi_conformal_error`

* `layout_isometric` — the reference layout `(0,0)`, `(ℓ₀₁,0)`, `(x,y)` has the given side
  lengths, and `y > 0`;
* `sda_layout` — the reference layout has signed double area `2·area`;
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

## `signed-angle-2d` — `signed_corner_angle`

* `signedCornerAngle_eq_arg` — it is `arg(conj U · V)`;
* `signedCornerAngle_mem_Ioc` — it lies in `(−π, π]`;
* `abs_signedCornerAngle_eq_angle` — its absolute value is the unsigned angle;
* `signedCornerAngle_sign` — its sign is the sign of `sda apex b c`; `0` or `π` if that is zero;
* `signedCornerAngle_sum_eq_pi` — in a counter-clockwise triangle the three angles are
  positive and sum to `π`.

## `segment-proper-crossing` — `proper_crossing`

* `strictlyOpposite_zero_iff`, `strictlyOpposite_imp` — the sign test without and with a
  tolerance;
* `properCrossing_sound` — if both tests hold, the returned point is interior to both segments;
* `properCrossing_complete` — an interior common point of non-parallel segments passes both
  tests, and it is the returned point (so it is unique);
* `properCrossing_iff` — both as one equivalence.

## `similarity-fit` — `fit_similarity`

* `simFit_residual_eq_add` — `Σ|a'p_i + b' − q_i|² = Σ|a p_i + b − q_i|² +
  Σ|(a' − a)p_i + (b' − b)|²`;
* `simFit_residual_le` — the fitted `(a, b)` minimizes the sum of squared residuals;
* `simFit_residual_eq_iff` — and is the only minimizer;
* `toC_fitMean`, `fitNorm_eq`, `fitRe_eq`, `fitIm_eq`, `toC_fitCoeff`, `toC_fitTranslation` —
  the real sums of the Rust code are the real and imaginary parts of the complex expressions;
* `fit_similarity_optimal`, `fit_similarity_unique` — hence the evaluated coefficient and
  translation are the unique least-squares optimum;
* `fitNorm_eq_zero_iff` — `norm = 0` iff all `p_i` coincide.

## `intrinsic-mollification` — `mollify_edge_lengths`

* `margin_iff_mollifyDeficit_le` — a face satisfies the inequalities with margin at `ε` iff
  its deficit is `≤ ε`;
* `mollifyNeeded_le_iff` — `needed ≤ ε` iff `0 ≤ ε` and the face satisfies them at `ε`;
* `mollifyEpsilonLoop_le_iff` — the loop returns the least `ε ≥ 0` at which all faces do;
* `mollifyEpsilonLoop_eq`, `mollifyEpsilonLoop_ofFn` — it is `ε* = max(0, max deficit)`;
* `mollifyEpsilon_le_iff`, `mollifyEpsilon_nonneg`, `mollifyEpsilon_feasible`,
  `mollifyEpsilon_le_of_feasible` — `ε* ≥ 0`, `ε*` is admissible, and every admissible
  `ε ≥ 0` is `≥ ε*`;
* `mollifyEpsilon_eq_zero_iff` — `ε* = 0` iff all faces already satisfy the inequalities;
* `isTriangle_of_margin`, `isTriangle_mollified` — for `δ > 0` the mollified sides of every
  face satisfy the strict triangle inequalities (nothing else is needed).

Nothing is claimed about floating-point rounding.
-/

noncomputable section

namespace DdgProofs

open Real

/-! ## Plane vectors -/

/-- Dot product of two plane vectors, `u·v = u₁v₁ + u₂v₂`. -/
def dot2 (u v : ℝ × ℝ) : ℝ := u.1 * v.1 + u.2 * v.2

/-- Cross product (determinant) of two plane vectors, `u×v = u₁v₂ − u₂v₁`. -/
def cross2 (u v : ℝ × ℝ) : ℝ := u.1 * v.2 - u.2 * v.1

/-- Squared Euclidean distance of two plane points, `Δx² + Δy²`. -/
def sqDist2 (p q : ℝ × ℝ) : ℝ := (q.1 - p.1) ^ 2 + (q.2 - p.2) ^ 2

/-- Euclidean distance of two plane points, `√(Δx² + Δy²)`. (The distance `dist p q` of the type
`ℝ × ℝ` is the sup metric and is not used.) -/
def dist2 (p q : ℝ × ℝ) : ℝ := √((q.1 - p.1) ^ 2 + (q.2 - p.2) ^ 2)

/-- The point `(x, y)` as the complex number `x + i·y`. -/
def toC (p : ℝ × ℝ) : ℂ := ⟨p.1, p.2⟩

theorem sqDist2_nonneg (p q : ℝ × ℝ) : 0 ≤ sqDist2 p q := by unfold sqDist2; positivity

theorem dist2_sq (p q : ℝ × ℝ) : dist2 p q ^ 2 = sqDist2 p q :=
  sq_sqrt (sqDist2_nonneg p q)

/-! ## Signed area -/

/-- `signed_double_area`: twice the signed area of the plane triangle `(a, b, c)`,
`(b₁ − a₁)(c₂ − a₂) − (b₂ − a₂)(c₁ − a₁)`. -/
def sda (a b c : ℝ × ℝ) : ℝ := (b.1 - a.1) * (c.2 - a.2) - (b.2 - a.2) * (c.1 - a.1)

/-- `sda a b c = (b − a) × (c − a)`. -/
theorem sda_eq_cross2 (a b c : ℝ × ℝ) : sda a b c = cross2 (b - a) (c - a) := rfl

/-- Swapping the first two corners changes the sign.

[formula: signed-area-2d] -/
theorem sda_swap_left (a b c : ℝ × ℝ) : sda b a c = -sda a b c := by unfold sda; ring

/-- Swapping the last two corners changes the sign.

[formula: signed-area-2d] -/
theorem sda_swap_right (a b c : ℝ × ℝ) : sda a c b = -sda a b c := by unfold sda; ring

/-- Swapping the first and the last corner changes the sign.

[formula: signed-area-2d] -/
theorem sda_swap_outer (a b c : ℝ × ℝ) : sda c b a = -sda a b c := by unfold sda; ring

/-- Cyclic rotation of the corners does not change the signed area.

[formula: signed-area-2d] -/
theorem sda_rotate (a b c : ℝ × ℝ) : sda b c a = sda a b c := by unfold sda; ring

/-- Translation does not change the signed area.

[formula: signed-area-2d] -/
theorem sda_translate (a b c t : ℝ × ℝ) : sda (a + t) (b + t) (c + t) = sda a b c := by
  simp only [sda, Prod.fst_add, Prod.snd_add]
  ring

/-- `4·(2·signed area)² = 2xy + 2yz + 2zx − x² − y² − z²` with `x, y, z` the squared lengths of
the sides `bc`, `ca`, `ab`: sixteen times the squared area, as in Heron's formula
(`sixteen_mul_triArea_sq`).

[formula: signed-area-2d] -/
theorem sda_sq_mul_four (a b c : ℝ × ℝ) :
    sda a b c ^ 2 * 4 =
      2 * sqDist2 b c * sqDist2 c a + 2 * sqDist2 c a * sqDist2 a b +
        2 * sqDist2 a b * sqDist2 b c - sqDist2 b c ^ 2 - sqDist2 c a ^ 2 - sqDist2 a b ^ 2 := by
  unfold sda sqDist2
  ring

/-- **The signed double area has the absolute value of twice the Heron area** of the Euclidean
side lengths. No nondegeneracy is needed: for collinear points both sides are zero.

[formula: signed-area-2d] -/
theorem abs_sda_eq_two_mul_triArea (a b c : ℝ × ℝ) :
    |sda a b c| = 2 * triArea (dist2 b c) (dist2 c a) (dist2 a b) := by
  have key : (dist2 b c + dist2 c a + dist2 a b) * (dist2 c a - dist2 b c + dist2 a b) *
      (dist2 b c - dist2 c a + dist2 a b) * (dist2 b c - dist2 a b + dist2 c a) =
      (2 * sda a b c) ^ 2 := by
    calc _ = 2 * dist2 b c ^ 2 * dist2 c a ^ 2 + 2 * dist2 c a ^ 2 * dist2 a b ^ 2 +
          2 * dist2 a b ^ 2 * dist2 b c ^ 2 - (dist2 b c ^ 2) ^ 2 - (dist2 c a ^ 2) ^ 2 -
          (dist2 a b ^ 2) ^ 2 := by ring
      _ = _ := by
        rw [dist2_sq, dist2_sq, dist2_sq]
        linear_combination (-1 : ℝ) * sda_sq_mul_four a b c
  unfold triArea
  rw [key, sqrt_sq_eq_abs, abs_mul, abs_two]
  ring

/-! ## Singular values of a triangle map -/

/-- `x` of the reference layout `p₀ = (0,0)`, `p₁ = (ℓ₀₁, 0)`, `p₂ = (x, y)`, as evaluated:
`(ℓ₀₁² + (ℓ₂₀ − ℓ₁₂)(ℓ₂₀ + ℓ₁₂)) / (2 ℓ₀₁)`. -/
def layoutX (l01 l12 l20 : ℝ) : ℝ := (l01 * l01 + (l20 - l12) * (l20 + l12)) / (2 * l01)

/-- `y` of the reference layout: `2·area / ℓ₀₁`. -/
def layoutY (l01 area : ℝ) : ℝ := 2 * area / l01

/-- **The reference layout is isometric.** For a nondegenerate triangle with sides
`ℓ₀₁, ℓ₁₂, ℓ₂₀` and its (Heron) area, the points `p₀ = (0,0)`, `p₁ = (ℓ₀₁, 0)`, `p₂ = (x, y)`
have `|p₁ − p₀| = ℓ₀₁`, `|p₂ − p₀| = ℓ₂₀`, `|p₂ − p₁| = ℓ₁₂`, and `y > 0` (counter-clockwise).

[formula: quasi-conformal-error] -/
theorem layout_isometric {l01 l12 l20 area : ℝ} (h : IsTriangle l01 l12 l20)
    (hA : area = triArea l01 l12 l20) :
    0 < layoutY l01 area ∧
      dist2 (0, 0) (l01, 0) = l01 ∧
      dist2 (0, 0) (layoutX l01 l12 l20, layoutY l01 area) = l20 ∧
      dist2 (l01, 0) (layoutX l01 l12 l20, layoutY l01 area) = l12 := by
  have h0 := h.pos_a
  have hT := sixteen_mul_triArea_sq h
  have hTpos := triArea_pos h
  rw [← hA] at hT hTpos
  have hx : layoutX l01 l12 l20 * (2 * l01) = l01 * l01 + (l20 - l12) * (l20 + l12) :=
    div_mul_cancel₀ _ (by positivity)
  have hy : layoutY l01 area * l01 = 2 * area := div_mul_cancel₀ _ h0.ne'
  have h4 : (2 * l01) ^ 2 ≠ 0 := by positivity
  have e20 : layoutX l01 l12 l20 ^ 2 + layoutY l01 area ^ 2 = l20 ^ 2 := by
    apply mul_right_cancel₀ h4
    linear_combination (layoutX l01 l12 l20 * (2 * l01) + (l01 * l01 + (l20 - l12) * (l20 + l12))) *
      hx + 4 * (layoutY l01 area * l01 + 2 * area) * hy + hT
  have e12 : (layoutX l01 l12 l20 - l01) ^ 2 + layoutY l01 area ^ 2 = l12 ^ 2 := by
    linear_combination e20 - hx
  refine ⟨div_pos (by positivity) h0, ?_, ?_, ?_⟩
  · unfold dist2
    simp only [sub_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, add_zero]
    exact sqrt_sq h0.le
  · unfold dist2
    simp only [sub_zero]
    rw [e20, sqrt_sq h.pos_c.le]
  · unfold dist2
    simp only [sub_zero]
    rw [e12, sqrt_sq h.pos_b.le]

/-- The hypotheses of `layout_isometric` are satisfiable: the triangle with
`[ℓ₀₁, ℓ₁₂, ℓ₂₀] = [3, 5, 4]` has area `6` and is laid out as `(0,0)`, `(3,0)`, `(0,4)`. -/
example : IsTriangle (3 : ℝ) 5 4 ∧ (6 : ℝ) = triArea 3 5 4 ∧ layoutX 3 5 4 = 0 ∧
    layoutY 3 6 = 4 := by
  refine ⟨⟨by norm_num, by norm_num, by norm_num⟩, ?_, by norm_num [layoutX],
    by norm_num [layoutY]⟩
  unfold triArea
  rw [show ((3 : ℝ) + 5 + 4) * (5 - 3 + 4) * (3 - 5 + 4) * (3 - 4 + 5) = 24 ^ 2 by norm_num,
    sqrt_sq (by norm_num)]
  norm_num

/-- The reference layout has signed double area `2·area` (for `ℓ₀₁ ≠ 0` and any `x`).

[formula: quasi-conformal-error] -/
theorem sda_layout {l01 : ℝ} (hl : l01 ≠ 0) (x area : ℝ) :
    sda (0, 0) (l01, 0) (x, layoutY l01 area) = 2 * area := by
  unfold sda layoutY
  simp only [sub_zero, zero_mul]
  field_simp

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

/-! ## Signed corner angle -/

/-- `signed_corner_angle`: `atan2(u × v, u · v)` with `u = b − apex`, `v = c − apex`, as
evaluated. -/
def signedCornerAngle (apex b c : ℝ × ℝ) : ℝ :=
  atan2 ((b.1 - apex.1) * (c.2 - apex.2) - (b.2 - apex.2) * (c.1 - apex.1))
    ((b.1 - apex.1) * (c.1 - apex.1) + (b.2 - apex.2) * (c.2 - apex.2))

theorem toC_injective : Function.Injective toC :=
  fun _ _ h => Prod.ext (congrArg Complex.re h) (congrArg Complex.im h)

theorem toC_sub (p q : ℝ × ℝ) : toC (p - q) = toC p - toC q := by
  apply Complex.ext <;> simp [toC]

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

/-! ## Similarity fit -/

section SimilarityFit

open ComplexConjugate

variable {ι : Type*} [Fintype ι]

/-- Centroid `p̄ = (Σ p_i)/n` of finitely many complex numbers. -/
def centroid (p : ι → ℂ) : ℂ := (∑ i, p i) / (Fintype.card ι : ℂ)

/-- The fitted coefficient `a = Σ (q_i − q̄)·conj(p_i − p̄) / Σ |p_i − p̄|²`. -/
def simFitA (p q : ι → ℂ) : ℂ :=
  (∑ i, (q i - centroid q) * conj (p i - centroid p)) /
    ((∑ i, ‖p i - centroid p‖ ^ 2 : ℝ) : ℂ)

/-- The fitted translation `b = q̄ − a·p̄`. -/
def simFitB (p q : ι → ℂ) : ℂ := centroid q - simFitA p q * centroid p

theorem sum_sub_centroid [Nonempty ι] (p : ι → ℂ) : ∑ i, (p i - centroid p) = 0 := by
  have hn : (Fintype.card ι : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul, centroid,
    mul_div_cancel₀ _ hn, sub_self]

/-- The residual of the fit is `a (p_i − p̄) − (q_i − q̄)`. -/
theorem simFit_residual (p q : ι → ℂ) (i : ι) :
    simFitA p q * p i + simFitB p q - q i =
      simFitA p q * (p i - centroid p) - (q i - centroid q) := by
  unfold simFitB
  ring

/-- The residuals of the fit sum to zero. -/
theorem simFit_sum_residual [Nonempty ι] (p q : ι → ℂ) :
    ∑ i, (simFitA p q * p i + simFitB p q - q i) = 0 := by
  simp only [simFit_residual]
  rw [Finset.sum_sub_distrib, ← Finset.mul_sum, sum_sub_centroid, sum_sub_centroid, mul_zero,
    sub_zero]

/-- The residuals of the fit are orthogonal to the points: `Σ r_i·conj(p_i) = 0` (the normal
equation for `a`). -/
theorem simFit_sum_residual_mul_conj [Nonempty ι] (p q : ι → ℂ)
    (hp : ∑ i, ‖p i - centroid p‖ ^ 2 ≠ 0) :
    ∑ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i) = 0 := by
  have hN : ((∑ i, ‖p i - centroid p‖ ^ 2 : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hp
  have hcentred : ∑ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i - centroid p) = 0 := by
    have hterm : ∀ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i - centroid p) =
        simFitA p q * ((p i - centroid p) * conj (p i - centroid p)) -
          (q i - centroid q) * conj (p i - centroid p) := fun i => by
      rw [simFit_residual]
      ring
    have hsq : ∑ i, (p i - centroid p) * conj (p i - centroid p) =
        ((∑ i, ‖p i - centroid p‖ ^ 2 : ℝ) : ℂ) := by
      rw [Complex.ofReal_sum]
      exact Finset.sum_congr rfl fun i _ => by rw [Complex.mul_conj', Complex.ofReal_pow]
    simp only [hterm]
    rw [Finset.sum_sub_distrib, ← Finset.mul_sum, hsq, simFitA, div_mul_cancel₀ _ hN, sub_self]
  have hsplit : ∑ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i) =
      ∑ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i - centroid p) +
        (∑ i, (simFitA p q * p i + simFitB p q - q i)) * conj (centroid p) := by
    rw [Finset.sum_mul, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [map_sub]
    ring
  rw [hsplit, hcentred, simFit_sum_residual, zero_mul, add_zero]

/-- **Pythagoras for the similarity fit.** For any similarity `z ↦ a'z + b'`, the sum of squared
residuals exceeds that of the fitted `(a, b)` by `Σ |(a' − a) p_i + (b' − b)|²`. Hypotheses: at
least one point, and not all `p_i` equal (`Σ |p_i − p̄|² ≠ 0`).

[formula: similarity-fit] -/
theorem simFit_residual_eq_add [Nonempty ι] (p q : ι → ℂ)
    (hp : ∑ i, ‖p i - centroid p‖ ^ 2 ≠ 0) (a' b' : ℂ) :
    ∑ i, ‖a' * p i + b' - q i‖ ^ 2 =
      ∑ i, ‖simFitA p q * p i + simFitB p q - q i‖ ^ 2 +
        ∑ i, ‖(a' - simFitA p q) * p i + (b' - simFitB p q)‖ ^ 2 := by
  have hr0 := simFit_sum_residual p q
  have hr1 := simFit_sum_residual_mul_conj p q hp
  generalize simFitA p q = a at *
  generalize simFitB p q = b at *
  have hsum : ∑ i, (a * p i + b - q i) * conj ((a' - a) * p i + (b' - b)) =
      conj (a' - a) * ∑ i, (a * p i + b - q i) * conj (p i) +
        conj (b' - b) * ∑ i, (a * p i + b - q i) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [map_add, map_mul]
    ring
  have hcross : ∑ i, ((a * p i + b - q i) * conj ((a' - a) * p i + (b' - b))).re = 0 := by
    rw [← Complex.re_sum, hsum, hr0, hr1, mul_zero, mul_zero, add_zero, Complex.zero_re]
  have hterm : ∀ i, ‖a' * p i + b' - q i‖ ^ 2 =
      ‖a * p i + b - q i‖ ^ 2 + ‖(a' - a) * p i + (b' - b)‖ ^ 2 +
        2 * ((a * p i + b - q i) * conj ((a' - a) * p i + (b' - b))).re := by
    intro i
    rw [Complex.sq_norm, Complex.sq_norm, Complex.sq_norm, ← Complex.normSq_add]
    congr 1
    ring
  simp only [hterm]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.mul_sum, hcross, mul_zero,
    add_zero]

/-- **The similarity fit is the least-squares optimum**: no similarity `z ↦ a'z + b'` has a
smaller sum of squared residuals `Σ |a' p_i + b' − q_i|²` than the fitted
`a = Σ (q_i − q̄)·conj(p_i − p̄) / Σ |p_i − p̄|²`, `b = q̄ − a·p̄`.

[formula: similarity-fit] -/
theorem simFit_residual_le [Nonempty ι] (p q : ι → ℂ)
    (hp : ∑ i, ‖p i - centroid p‖ ^ 2 ≠ 0) (a' b' : ℂ) :
    ∑ i, ‖simFitA p q * p i + simFitB p q - q i‖ ^ 2 ≤ ∑ i, ‖a' * p i + b' - q i‖ ^ 2 := by
  rw [simFit_residual_eq_add p q hp a' b']
  exact le_add_of_nonneg_right (Finset.sum_nonneg fun i _ => sq_nonneg _)

/-- **The optimum is unique**: a similarity attains the minimal sum of squared residuals only
if it is the fitted one.

[formula: similarity-fit] -/
theorem simFit_residual_eq_iff [Nonempty ι] (p q : ι → ℂ)
    (hp : ∑ i, ‖p i - centroid p‖ ^ 2 ≠ 0) (a' b' : ℂ) :
    ∑ i, ‖a' * p i + b' - q i‖ ^ 2 = ∑ i, ‖simFitA p q * p i + simFitB p q - q i‖ ^ 2 ↔
      a' = simFitA p q ∧ b' = simFitB p q := by
  constructor
  · intro h
    rw [simFit_residual_eq_add p q hp a' b'] at h
    have h0 : ∑ i, ‖(a' - simFitA p q) * p i + (b' - simFitB p q)‖ ^ 2 = 0 := by linarith
    rw [Finset.sum_eq_zero_iff_of_nonneg fun i _ => sq_nonneg _] at h0
    have hi : ∀ i, (a' - simFitA p q) * p i + (b' - simFitB p q) = 0 := fun i => by
      simpa using h0 i (Finset.mem_univ i)
    have hn : (Fintype.card ι : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
    have hs : ∑ i, ((a' - simFitA p q) * p i + (b' - simFitB p q)) = 0 :=
      Finset.sum_eq_zero fun i _ => hi i
    rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_const, Finset.card_univ,
      nsmul_eq_mul] at hs
    have hc : (a' - simFitA p q) * centroid p + (b' - simFitB p q) = 0 := by
      unfold centroid
      field_simp
      linear_combination hs
    have ha : a' - simFitA p q = 0 := by
      by_contra hne
      apply hp
      refine Finset.sum_eq_zero fun i _ => ?_
      have : (a' - simFitA p q) * (p i - centroid p) = 0 := by linear_combination hi i - hc
      rw [(mul_eq_zero.mp this).resolve_left hne, norm_zero]
      norm_num
    obtain ⟨i⟩ := ‹Nonempty ι›
    have hb := hi i
    rw [ha, zero_mul, zero_add] at hb
    exact ⟨sub_eq_zero.mp ha, sub_eq_zero.mp hb⟩
  · rintro ⟨rfl, rfl⟩
    rfl

/-! ### The real sums of `fit_similarity` -/

/-- `mean` in `fit_similarity`: `(Σ p_i[0] / n, Σ p_i[1] / n)`. -/
def fitMean (p : ι → ℝ × ℝ) : ℝ × ℝ :=
  ((∑ i, (p i).1) / Fintype.card ι, (∑ i, (p i).2) / Fintype.card ι)

/-- `re = Σ (qx·px + qy·py)` on centred coordinates. -/
def fitRe (p q : ι → ℝ × ℝ) : ℝ :=
  ∑ i, (((q i).1 - (fitMean q).1) * ((p i).1 - (fitMean p).1) +
    ((q i).2 - (fitMean q).2) * ((p i).2 - (fitMean p).2))

/-- `im = Σ (qy·px − qx·py)` on centred coordinates. -/
def fitIm (p q : ι → ℝ × ℝ) : ℝ :=
  ∑ i, (((q i).2 - (fitMean q).2) * ((p i).1 - (fitMean p).1) -
    ((q i).1 - (fitMean q).1) * ((p i).2 - (fitMean p).2))

/-- `norm = Σ (px·px + py·py)` on centred coordinates. -/
def fitNorm (p : ι → ℝ × ℝ) : ℝ :=
  ∑ i, (((p i).1 - (fitMean p).1) * ((p i).1 - (fitMean p).1) +
    ((p i).2 - (fitMean p).2) * ((p i).2 - (fitMean p).2))

/-- `(ar, ai) = (re / norm, im / norm)`. -/
def fitCoeff (p q : ι → ℝ × ℝ) : ℝ × ℝ := (fitRe p q / fitNorm p, fitIm p q / fitNorm p)

/-- `translation = (qm[0] − (ar·pm[0] − ai·pm[1]), qm[1] − (ai·pm[0] + ar·pm[1]))`. -/
def fitTranslation (p q : ι → ℝ × ℝ) : ℝ × ℝ :=
  ((fitMean q).1 - ((fitCoeff p q).1 * (fitMean p).1 - (fitCoeff p q).2 * (fitMean p).2),
    (fitMean q).2 - ((fitCoeff p q).2 * (fitMean p).1 + (fitCoeff p q).1 * (fitMean p).2))

/-- The residual length `d = hypot(dx, dy)` of one pair of points under the similarity with
coefficient `ar + i·ai` and translation `tr`: `dx = ar·px − ai·py + tr[0] − qx`,
`dy = ai·px + ar·py + tr[1] − qy`. -/
def fitDist (ar ai : ℝ) (tr p q : ℝ × ℝ) : ℝ :=
  √((ar * p.1 - ai * p.2 + tr.1 - q.1) ^ 2 + (ai * p.1 + ar * p.2 + tr.2 - q.2) ^ 2)

/-- `sum = Σ d·d`, the sum of squared residual lengths (`rms = √(sum / n)`). -/
def fitSqError (ar ai : ℝ) (tr : ℝ × ℝ) (p q : ι → ℝ × ℝ) : ℝ :=
  ∑ i, fitDist ar ai tr (p i) (q i) * fitDist ar ai tr (p i) (q i)

/-- The mean is the centroid.

[formula: similarity-fit] -/
theorem toC_fitMean (p : ι → ℝ × ℝ) : toC (fitMean p) = centroid fun i => toC (p i) := by
  apply Complex.ext <;> simp [toC, fitMean, centroid, Complex.div_natCast_re,
    Complex.div_natCast_im]

/-- `norm = Σ |p_i − p̄|²`.

[formula: similarity-fit] -/
theorem fitNorm_eq (p : ι → ℝ × ℝ) :
    fitNorm p = ∑ i, ‖toC (p i) - centroid fun i => toC (p i)‖ ^ 2 := by
  rw [← toC_fitMean]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Complex.sq_norm, Complex.normSq_apply]
  simp [toC]

/-- `re` is the real part of `Σ (q_i − q̄)·conj(p_i − p̄)`.

[formula: similarity-fit] -/
theorem fitRe_eq (p q : ι → ℝ × ℝ) :
    fitRe p q = (∑ i, (toC (q i) - centroid fun i => toC (q i)) *
      conj (toC (p i) - centroid fun i => toC (p i))).re := by
  rw [← toC_fitMean, ← toC_fitMean, Complex.re_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp [toC]
  ring

/-- `im` is the imaginary part of `Σ (q_i − q̄)·conj(p_i − p̄)`.

[formula: similarity-fit] -/
theorem fitIm_eq (p q : ι → ℝ × ℝ) :
    fitIm p q = (∑ i, (toC (q i) - centroid fun i => toC (q i)) *
      conj (toC (p i) - centroid fun i => toC (p i))).im := by
  rw [← toC_fitMean, ← toC_fitMean, Complex.im_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp [toC]
  ring

/-- `ar + i·ai` is the complex coefficient `a` of the fit.

[formula: similarity-fit] -/
theorem toC_fitCoeff (p q : ι → ℝ × ℝ) :
    toC (fitCoeff p q) = simFitA (fun i => toC (p i)) (fun i => toC (q i)) := by
  unfold simFitA
  rw [← fitNorm_eq]
  apply Complex.ext
  · rw [Complex.div_ofReal_re, ← fitRe_eq]
    rfl
  · rw [Complex.div_ofReal_im, ← fitIm_eq]
    rfl

/-- The translation formula is `b = q̄ − a·p̄` in real and imaginary parts.

[formula: similarity-fit] -/
theorem toC_fitTranslation (p q : ι → ℝ × ℝ) :
    toC (fitTranslation p q) = simFitB (fun i => toC (p i)) (fun i => toC (q i)) := by
  unfold simFitB
  rw [← toC_fitCoeff, ← toC_fitMean, ← toC_fitMean]
  apply Complex.ext <;> simp [toC, fitTranslation]
  ring

/-- The residual length is `|a·p + b − q|` for `a = ar + i·ai`, `b = tr`. -/
theorem fitDist_eq_norm (ar ai : ℝ) (tr p q : ℝ × ℝ) :
    fitDist ar ai tr p q = ‖(⟨ar, ai⟩ : ℂ) * toC p + toC tr - toC q‖ := by
  rw [Complex.norm_eq_sqrt_sq_add_sq]
  unfold fitDist
  congr 1
  simp [toC]
  ring

theorem fitSqError_eq (ar ai : ℝ) (tr : ℝ × ℝ) (p q : ι → ℝ × ℝ) :
    fitSqError ar ai tr p q = ∑ i, ‖(⟨ar, ai⟩ : ℂ) * toC (p i) + toC tr - toC (q i)‖ ^ 2 := by
  unfold fitSqError
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [fitDist_eq_norm, sq]

theorem fitNorm_nonneg (p : ι → ℝ × ℝ) : 0 ≤ fitNorm p :=
  Finset.sum_nonneg fun _ _ => add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)

/-- `norm = 0` (the case in which `fit_similarity` returns `None`) exactly when all points `p_i`
coincide.

[formula: similarity-fit] -/
theorem fitNorm_eq_zero_iff [Nonempty ι] (p : ι → ℝ × ℝ) :
    fitNorm p = 0 ↔ ∀ i j, p i = p j := by
  have hmean : fitNorm p = 0 ↔ ∀ i, p i = fitMean p := by
    unfold fitNorm
    rw [Finset.sum_eq_zero_iff_of_nonneg fun i _ =>
      add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)]
    constructor
    · intro h i
      have hi := h i (Finset.mem_univ i)
      have h1 : (p i).1 - (fitMean p).1 = 0 := by
        nlinarith [mul_self_nonneg ((p i).1 - (fitMean p).1),
          mul_self_nonneg ((p i).2 - (fitMean p).2)]
      have h2 : (p i).2 - (fitMean p).2 = 0 := by
        nlinarith [mul_self_nonneg ((p i).1 - (fitMean p).1),
          mul_self_nonneg ((p i).2 - (fitMean p).2)]
      exact Prod.ext (sub_eq_zero.mp h1) (sub_eq_zero.mp h2)
    · intro h i _
      rw [← h i]
      ring
  rw [hmean]
  constructor
  · intro h i j
    rw [h i, h j]
  · intro h i
    have hn : (Fintype.card ι : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
    unfold fitMean
    rw [Finset.sum_congr rfl fun j _ => congrArg Prod.fst (h j i),
      Finset.sum_congr rfl fun j _ => congrArg Prod.snd (h j i)]
    simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    exact Prod.ext (by field_simp) (by field_simp)

/-- **`fit_similarity` returns the least-squares similarity.** With at least one point and not
all `p_i` equal (`norm ≠ 0`; otherwise the function returns `None`), the evaluated coefficient
`(ar, ai)` and translation give a sum of squared residuals `sum = Σ d·d` that is at most that
of any other coefficient `(ar', ai')` and translation `tr'`. (Hence also the smallest `rms`.)

[formula: similarity-fit] -/
theorem fit_similarity_optimal [Nonempty ι] (p q : ι → ℝ × ℝ) (hp : fitNorm p ≠ 0)
    (ar' ai' : ℝ) (tr' : ℝ × ℝ) :
    fitSqError (fitCoeff p q).1 (fitCoeff p q).2 (fitTranslation p q) p q ≤
      fitSqError ar' ai' tr' p q := by
  rw [fitSqError_eq, fitSqError_eq]
  rw [fitNorm_eq] at hp
  have h := simFit_residual_le (fun i => toC (p i)) (fun i => toC (q i)) hp ⟨ar', ai'⟩ (toC tr')
  rw [← toC_fitCoeff, ← toC_fitTranslation] at h
  exact h

/-- The same with uniqueness: another coefficient and translation give the same sum of squared
residuals only if they are the evaluated ones.

[formula: similarity-fit] -/
theorem fit_similarity_unique [Nonempty ι] (p q : ι → ℝ × ℝ) (hp : fitNorm p ≠ 0)
    (ar' ai' : ℝ) (tr' : ℝ × ℝ)
    (h : fitSqError ar' ai' tr' p q =
      fitSqError (fitCoeff p q).1 (fitCoeff p q).2 (fitTranslation p q) p q) :
    (ar', ai') = fitCoeff p q ∧ tr' = fitTranslation p q := by
  rw [fitSqError_eq, fitSqError_eq] at h
  rw [fitNorm_eq] at hp
  have h' := (simFit_residual_eq_iff (fun i => toC (p i)) (fun i => toC (q i)) hp ⟨ar', ai'⟩
    (toC tr')).mp
  rw [← toC_fitCoeff, ← toC_fitTranslation] at h'
  obtain ⟨ha, hb⟩ := h' h
  exact ⟨toC_injective ha, toC_injective hb⟩

end SimilarityFit

/-! ## Intrinsic mollification -/

/-- The triangle inequalities of a face with sides `a, b, c` hold with margin `δ` after every
side is lengthened by `ε`: `(x+ε) + (y+ε) ≥ (z+ε) + δ` for each side `z`. -/
def MarginHolds (a b c δ ε : ℝ) : Prop :=
  (b + ε) + (c + ε) ≥ (a + ε) + δ ∧ (c + ε) + (a + ε) ≥ (b + ε) + δ ∧
    (a + ε) + (b + ε) ≥ (c + ε) + δ

/-- By how much a face misses the margin: the largest of `ℓ_k − ℓ_{k+1} − ℓ_{k+2} + δ`. -/
def mollifyDeficit (a b c δ : ℝ) : ℝ :=
  max (a - b - c + δ) (max (b - c - a + δ) (c - a - b + δ))

/-- `needed` as the loop over the three corners computes it: starting from `0`, the maximum
with `ℓ_k − ℓ_{k+1} − ℓ_{k+2} + δ` for `k = 0, 1, 2`. -/
def mollifyNeeded (a b c δ : ℝ) : ℝ :=
  max (max (max 0 (a - b - c + δ)) (b - c - a + δ)) (c - a - b + δ)

/-- One pass of the loop over the faces: `if needed > 0 { if needed > ε { ε = needed } }`. -/
def mollifyStep (δ ε : ℝ) (l : ℝ × ℝ × ℝ) : ℝ :=
  if 0 < mollifyNeeded l.1 l.2.1 l.2.2 δ then
    if ε < mollifyNeeded l.1 l.2.1 l.2.2 δ then mollifyNeeded l.1 l.2.1 l.2.2 δ else ε
  else ε

/-- `epsilon` as the loop over the faces computes it, starting from `0`. A face is the triple
of its side lengths. -/
def mollifyEpsilonLoop (δ : ℝ) (faces : List (ℝ × ℝ × ℝ)) : ℝ := faces.foldl (mollifyStep δ) 0

/-- **The margin inequalities of a face hold at `ε` iff its deficit is at most `ε`** (for every
real `ε`).

[formula: intrinsic-mollification] -/
theorem margin_iff_mollifyDeficit_le (a b c δ ε : ℝ) :
    MarginHolds a b c δ ε ↔ mollifyDeficit a b c δ ≤ ε := by
  unfold MarginHolds mollifyDeficit
  rw [max_le_iff, max_le_iff]
  constructor
  · rintro ⟨h1, h2, h3⟩
    exact ⟨by linarith, by linarith, by linarith⟩
  · rintro ⟨h1, h2, h3⟩
    exact ⟨by linarith, by linarith, by linarith⟩

/-- The loop over the corners starts from `0`: `needed = max(0, deficit)`. -/
theorem mollifyNeeded_eq (a b c δ : ℝ) :
    mollifyNeeded a b c δ = max 0 (mollifyDeficit a b c δ) := by
  unfold mollifyNeeded mollifyDeficit
  rw [max_assoc, max_assoc]

/-- **`needed ≤ ε` iff `ε ≥ 0` and the margin inequalities of the face hold at `ε`.**

[formula: intrinsic-mollification] -/
theorem mollifyNeeded_le_iff (a b c δ ε : ℝ) :
    mollifyNeeded a b c δ ≤ ε ↔ 0 ≤ ε ∧ MarginHolds a b c δ ε := by
  rw [mollifyNeeded_eq, max_le_iff, margin_iff_mollifyDeficit_le]

theorem mollifyNeeded_nonneg (a b c δ : ℝ) : 0 ≤ mollifyNeeded a b c δ := by
  rw [mollifyNeeded_eq]
  exact le_max_left _ _

private theorem mollifyStep_nonneg {δ ε₀ : ℝ} (h0 : 0 ≤ ε₀) (l : ℝ × ℝ × ℝ) :
    0 ≤ mollifyStep δ ε₀ l := by
  unfold mollifyStep
  split_ifs with h1 h2
  · exact h1.le
  · exact h0
  · exact h0

private theorem mollifyStep_le_iff {δ ε₀ : ℝ} (h0 : 0 ≤ ε₀) (l : ℝ × ℝ × ℝ) (ε : ℝ) :
    mollifyStep δ ε₀ l ≤ ε ↔ ε₀ ≤ ε ∧ mollifyNeeded l.1 l.2.1 l.2.2 δ ≤ ε := by
  have hn := mollifyNeeded_nonneg l.1 l.2.1 l.2.2 δ
  unfold mollifyStep
  split_ifs with h1 h2
  · exact ⟨fun h => ⟨by linarith, h⟩, fun h => h.2⟩
  · exact ⟨fun h => ⟨h, by linarith⟩, fun h => h.1⟩
  · exact ⟨fun h => ⟨h, by linarith⟩, fun h => h.1⟩

private theorem foldl_mollifyStep_le_iff (δ : ℝ) (faces : List (ℝ × ℝ × ℝ)) :
    ∀ ε₀ : ℝ, 0 ≤ ε₀ → ∀ ε : ℝ, faces.foldl (mollifyStep δ) ε₀ ≤ ε ↔
      ε₀ ≤ ε ∧ ∀ l ∈ faces, mollifyNeeded l.1 l.2.1 l.2.2 δ ≤ ε := by
  induction faces with
  | nil =>
    intro ε₀ _ ε
    simp
  | cons l ls ih =>
    intro ε₀ h0 ε
    rw [List.foldl_cons, ih _ (mollifyStep_nonneg h0 l), mollifyStep_le_iff h0,
      List.forall_mem_cons, and_assoc]

/-- **The loop computes the least admissible `ε`**: `epsilon ≤ ε` iff `ε ≥ 0` and the margin
inequalities of all faces hold at `ε`. So `epsilon` is itself `≥ 0` and admissible (take
`ε = epsilon`), and every admissible `ε ≥ 0` is at least `epsilon`.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilonLoop_le_iff (δ : ℝ) (faces : List (ℝ × ℝ × ℝ)) (ε : ℝ) :
    mollifyEpsilonLoop δ faces ≤ ε ↔
      0 ≤ ε ∧ ∀ l ∈ faces, MarginHolds l.1 l.2.1 l.2.2 δ ε := by
  unfold mollifyEpsilonLoop
  rw [foldl_mollifyStep_le_iff δ faces 0 le_rfl]
  constructor
  · rintro ⟨h0, h⟩
    exact ⟨h0, fun l hl => ((mollifyNeeded_le_iff _ _ _ _ _).mp (h l hl)).2⟩
  · rintro ⟨h0, h⟩
    exact ⟨h0, fun l hl => (mollifyNeeded_le_iff _ _ _ _ _).mpr ⟨h0, h l hl⟩⟩

section Family

variable {F : Type*} [Fintype F] [Nonempty F]

/-- `ε* = max(0, max over the faces of the deficit)` for a finite nonempty family of faces
`l f = (a, b, c)`: the closed form of what `mollify_edge_lengths` adds to every edge. -/
def mollifyEpsilon (l : F → ℝ × ℝ × ℝ) (δ : ℝ) : ℝ :=
  max 0 (Finset.univ.sup' Finset.univ_nonempty fun f =>
    mollifyDeficit (l f).1 (l f).2.1 (l f).2.2 δ)

/-- **`ε* ≤ ε` iff `ε ≥ 0` and all faces satisfy the inequalities with margin at `ε`**: `ε*` is
the least such `ε`.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_le_iff (l : F → ℝ × ℝ × ℝ) (δ ε : ℝ) :
    mollifyEpsilon l δ ≤ ε ↔ 0 ≤ ε ∧ ∀ f, MarginHolds (l f).1 (l f).2.1 (l f).2.2 δ ε := by
  unfold mollifyEpsilon
  rw [max_le_iff, Finset.sup'_le_iff]
  simp only [Finset.mem_univ, true_implies, margin_iff_mollifyDeficit_le]

/-- **The loop computes `ε*`**, in whatever order (and however often) it visits the faces.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilonLoop_eq (l : F → ℝ × ℝ × ℝ) (δ : ℝ) {faces : List (ℝ × ℝ × ℝ)}
    (hfaces : ∀ x, x ∈ faces ↔ ∃ f, l f = x) :
    mollifyEpsilonLoop δ faces = mollifyEpsilon l δ := by
  refine eq_of_forall_ge_iff fun ε => ?_
  rw [mollifyEpsilonLoop_le_iff, mollifyEpsilon_le_iff]
  constructor
  · rintro ⟨h0, h⟩
    exact ⟨h0, fun f => h _ ((hfaces _).mpr ⟨f, rfl⟩)⟩
  · rintro ⟨h0, h⟩
    refine ⟨h0, fun x hx => ?_⟩
    obtain ⟨f, rfl⟩ := (hfaces x).mp hx
    exact h f

/-- The loop over the faces `0, …, n` in this order computes `ε*`. (For no faces at all the
loop returns `0`.)

[formula: intrinsic-mollification] -/
theorem mollifyEpsilonLoop_ofFn {n : ℕ} (l : Fin (n + 1) → ℝ × ℝ × ℝ) (δ : ℝ) :
    mollifyEpsilonLoop δ (List.ofFn l) = mollifyEpsilon l δ :=
  mollifyEpsilonLoop_eq l δ fun _ => List.mem_ofFn

/-- **`ε* ≥ 0`.**

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_nonneg (l : F → ℝ × ℝ × ℝ) (δ : ℝ) : 0 ≤ mollifyEpsilon l δ :=
  ((mollifyEpsilon_le_iff l δ _).mp le_rfl).1

/-- **`ε*` is admissible**: after lengthening every side by `ε*`, every face satisfies its three
triangle inequalities with margin `δ`.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_feasible (l : F → ℝ × ℝ × ℝ) (δ : ℝ) (f : F) :
    MarginHolds (l f).1 (l f).2.1 (l f).2.2 δ (mollifyEpsilon l δ) :=
  ((mollifyEpsilon_le_iff l δ _).mp le_rfl).2 f

/-- **`ε*` is the least admissible amount**: every `ε ≥ 0` for which all faces satisfy the
inequalities with margin is at least `ε*`.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_le_of_feasible (l : F → ℝ × ℝ × ℝ) (δ : ℝ) {ε : ℝ} (hε : 0 ≤ ε)
    (h : ∀ f, MarginHolds (l f).1 (l f).2.1 (l f).2.2 δ ε) : mollifyEpsilon l δ ≤ ε :=
  (mollifyEpsilon_le_iff l δ ε).mpr ⟨hε, h⟩

/-- The lengths are left unchanged (`ε* = 0`) exactly when all faces already satisfy the
inequalities with margin.

[formula: intrinsic-mollification] -/
theorem mollifyEpsilon_eq_zero_iff (l : F → ℝ × ℝ × ℝ) (δ : ℝ) :
    mollifyEpsilon l δ = 0 ↔ ∀ f, MarginHolds (l f).1 (l f).2.1 (l f).2.2 δ 0 := by
  constructor
  · intro h
    exact ((mollifyEpsilon_le_iff l δ 0).mp h.le).2
  · intro h
    exact le_antisymm ((mollifyEpsilon_le_iff l δ 0).mpr ⟨le_rfl, h⟩) (mollifyEpsilon_nonneg l δ)

end Family

/-- **A positive margin gives a nondegenerate triangle.** If the margin inequalities hold at
`ε` with `δ > 0`, the lengthened sides satisfy the strict triangle inequalities, and each of
them is at least `δ`. Nothing else is needed: no sign condition on `a, b, c` or `ε`.

[formula: intrinsic-mollification] -/
theorem isTriangle_of_margin {a b c δ ε : ℝ} (hδ : 0 < δ) (h : MarginHolds a b c δ ε) :
    IsTriangle (a + ε) (b + ε) (c + ε) ∧ δ ≤ a + ε ∧ δ ≤ b + ε ∧ δ ≤ c + ε := by
  obtain ⟨h1, h2, h3⟩ := h
  exact ⟨⟨by linarith, by linarith, by linarith⟩, by linarith, by linarith, by linarith⟩

/-- **The mollified metric is a metric of nondegenerate triangles** whenever `δ > 0`: for every
face, the sides lengthened by `ε*` satisfy `IsTriangle` (so they are positive and the face has
positive area). For `δ = 0` only the weak inequalities hold, and a degenerate face stays
degenerate (see the example below).

[formula: intrinsic-mollification] -/
theorem isTriangle_mollified {F : Type*} [Fintype F] [Nonempty F] (l : F → ℝ × ℝ × ℝ) {δ : ℝ}
    (hδ : 0 < δ) (f : F) :
    IsTriangle ((l f).1 + mollifyEpsilon l δ) ((l f).2.1 + mollifyEpsilon l δ)
        ((l f).2.2 + mollifyEpsilon l δ) ∧
      0 < triArea ((l f).1 + mollifyEpsilon l δ) ((l f).2.1 + mollifyEpsilon l δ)
        ((l f).2.2 + mollifyEpsilon l δ) := by
  have h := (isTriangle_of_margin hδ (mollifyEpsilon_feasible l δ f)).1
  exact ⟨h, triArea_pos h⟩

/-- `δ > 0` cannot be dropped: the degenerate face `(1, 1, 2)` needs no lengthening for
`δ = 0` and is not a triangle. -/
example : mollifyNeeded 1 1 2 0 = 0 ∧ ¬IsTriangle ((1 : ℝ) + 0) (1 + 0) (2 + 0) := by
  refine ⟨by norm_num [mollifyNeeded], fun h => ?_⟩
  have := h.lt_c
  norm_num at this

end DdgProofs
