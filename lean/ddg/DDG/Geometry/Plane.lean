import DDG.Geometry.Triangle

/-!
# Plane points: dot and cross products, distances, signed area and the reference layout

Transcriptions of the formulas of `crates/ddg-bff/src/quality.rs` and
`crates/ddg-bff/src/mollify.rs` into definitions on real numbers (same arrangement of terms),
and what they compute. Plane points are pairs `ℝ × ℝ`; lengths are written `√(Δx² + Δy²)`
(`dist2`), never with the metric of the type `ℝ × ℝ` (which is the sup metric). `dot2 u v` and
`cross2 u v` are `u·v` and `u×v`; `toC` turns a point into a complex number.

(This module was the first part of the former `Bff.Quality` module; the quality measures themselves are in
`DDG.Conformal.Distortion`, `DDG.Geometry.SignedAngle`,
`DDG.Geometry.SegmentCrossing`, `DDG.Geometry.SimilarityFit` and
`DDG.Intrinsic.Mollification`.)

## `signed-area-2d` — `signed_double_area` (`sda`)

* `sda_swap_left`, `sda_swap_right`, `sda_swap_outer` — swapping two corners changes the sign;
* `sda_rotate` — cyclic rotation of the corners does not change it;
* `sda_translate` — translation does not change it;
* `sda_sq_mul_four` — `sda²·4 = 2xy + 2yz + 2zx − x² − y² − z²` in the squared side lengths;
* `abs_sda_eq_two_mul_triArea` — `|sda| = 2·triArea` of the Euclidean side lengths (no
  nondegeneracy hypothesis).

## `quasi-conformal-error` — the reference layout

* `layout_isometric` — the reference layout `(0,0)`, `(ℓ₀₁,0)`, `(x,y)` has the given side
  lengths, and `y > 0`;
* `sda_layout` — the reference layout has signed double area `2·area`;

## Complex numbers

* `toC_injective`, `toC_sub` — `toC` is injective and additive.
-/

noncomputable section

namespace DDG

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
theorem toC_injective : Function.Injective toC :=
  fun _ _ h => Prod.ext (congrArg Complex.re h) (congrArg Complex.im h)

theorem toC_sub (p q : ℝ × ℝ) : toC (p - q) = toC p - toC q := by
  apply Complex.ext <;> simp [toC]

end DDG
