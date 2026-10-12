import Mathlib

/-!
# Euclidean triangles from their edge lengths

The library is intrinsic: every triangle quantity is a function of the three edge lengths
(paper Sec. 5.1.1). Convention (the same as in `lean/oracle`): in `cornerAngle a b c`,
`cornerCot a b c` the corner is the one **opposite the first argument** `a`, so it lies between
the sides `b` and `c`.

* `cornerAngle a b c = arccos((b² + c² − a²)/(2bc))` (law of cosines), also in the numerically
  stable half-angle form `2·arctan √(((a−b+c)(a+b−c))/((a+b+c)(b+c−a)))`
  (`cornerAngle_eq_two_arctan`);
* `triArea a b c = ¼ √((a+b+c)(b+c−a)(a−b+c)(a+b−c))` (Heron);
* `cornerCot a b c = (b² + c² − a²)/(4·Area)`, the cotangent of the corner angle
  (`cornerCot_eq_cot`);
* the three corner angles sum to `π` (`cornerAngle_add_cornerAngle_add_cornerAngle`);
* for points of a Euclidean space these are the angle and the area of the triangle
  (`angle_eq_cornerAngle`, `half_mul_sin_angle_eq_triArea`);
* behaviour under uniform scaling (`cornerAngle_smul`, `cornerCot_smul`, `triArea_smul`).
-/

noncomputable section

namespace DdgProofs

open Real

/-- `a b c` are the side lengths of a nondegenerate Euclidean triangle: the three strict triangle
inequalities hold. (Positivity of `a b c` follows, see `IsTriangle.pos_a`.) -/
structure IsTriangle (a b c : ℝ) : Prop where
  lt_a : a < b + c
  lt_b : b < c + a
  lt_c : c < a + b

namespace IsTriangle

variable {a b c : ℝ}

theorem pos_a (h : IsTriangle a b c) : 0 < a := by linarith [h.lt_b, h.lt_c]
theorem pos_b (h : IsTriangle a b c) : 0 < b := by linarith [h.lt_a, h.lt_c]
theorem pos_c (h : IsTriangle a b c) : 0 < c := by linarith [h.lt_a, h.lt_b]

/-- Cyclic relabelling. -/
theorem rotate (h : IsTriangle a b c) : IsTriangle b c a := ⟨h.lt_b, h.lt_c, h.lt_a⟩

/-- The side lengths `3, 4, 5` form a triangle: the hypothesis is satisfiable. -/
example : IsTriangle 3 4 5 := ⟨by norm_num, by norm_num, by norm_num⟩

end IsTriangle

/-- Area of the triangle with side lengths `a b c`, by Heron's formula:
`¼ √((a+b+c)(b+c−a)(a−b+c)(a+b−c))`. -/
def triArea (a b c : ℝ) : ℝ := √((a + b + c) * (b - a + c) * (a - b + c) * (a - c + b)) / 4

/-- Interior angle opposite the side `a` of the triangle with side lengths `a b c`, by the law of
cosines: `arccos((b² + c² − a²)/(2bc))`. -/
def cornerAngle (a b c : ℝ) : ℝ := arccos ((b ^ 2 + c ^ 2 - a ^ 2) / (2 * b * c))

/-- Cotangent of the interior angle opposite the side `a`: `(b² + c² − a²)/(4·Area)`. -/
def cornerCot (a b c : ℝ) : ℝ := (b ^ 2 + c ^ 2 - a ^ 2) / (4 * triArea a b c)

section

variable {a b c : ℝ}

theorem IsTriangle.heron_pos (h : IsTriangle a b c) :
    0 < (a + b + c) * (b - a + c) * (a - b + c) * (a - c + b) := by
  have h1 : 0 < a + b + c := by linarith [h.pos_a, h.pos_b, h.pos_c]
  have h2 : 0 < b - a + c := by linarith [h.lt_a]
  have h3 : 0 < a - b + c := by linarith [h.lt_b]
  have h4 : 0 < a - c + b := by linarith [h.lt_c]
  positivity

theorem triArea_pos (h : IsTriangle a b c) : 0 < triArea a b c :=
  div_pos (sqrt_pos.mpr h.heron_pos) (by norm_num)

/-- `16·Area² = (a+b+c)(b+c−a)(a−b+c)(a+b−c)`. -/
theorem sq_four_mul_triArea (h : IsTriangle a b c) :
    (4 * triArea a b c) ^ 2 = (a + b + c) * (b - a + c) * (a - b + c) * (a - c + b) := by
  unfold triArea
  rw [mul_div_cancel₀ _ (by norm_num : (4 : ℝ) ≠ 0), sq_sqrt h.heron_pos.le]

/-- `16·Area² = 2a²b² + 2b²c² + 2c²a² − a⁴ − b⁴ − c⁴`. -/
theorem sixteen_mul_triArea_sq (h : IsTriangle a b c) :
    16 * triArea a b c ^ 2 =
      2 * a ^ 2 * b ^ 2 + 2 * b ^ 2 * c ^ 2 + 2 * c ^ 2 * a ^ 2 - a ^ 4 - b ^ 4 - c ^ 4 := by
  have := sq_four_mul_triArea h
  nlinarith

/-- The area does not depend on the labelling of the sides. -/
theorem triArea_rotate (a b c : ℝ) : triArea b c a = triArea a b c := by
  unfold triArea
  congr 2
  ring

theorem triArea_swap (a b c : ℝ) : triArea a c b = triArea a b c := by
  unfold triArea
  congr 2
  ring

private theorem cos_arg_bounds (h : IsTriangle a b c) :
    -1 < (b ^ 2 + c ^ 2 - a ^ 2) / (2 * b * c) ∧ (b ^ 2 + c ^ 2 - a ^ 2) / (2 * b * c) < 1 := by
  have hb := h.pos_b
  have hc := h.pos_c
  have hbc : 0 < 2 * b * c := by positivity
  constructor
  · rw [lt_div_iff₀ hbc]
    nlinarith [h.lt_a, h.pos_a, mul_pos (sub_pos.mpr h.lt_a) (by linarith [h.pos_a] : 0 < b + c + a)]
  · rw [div_lt_iff₀ hbc]
    nlinarith [mul_pos (by linarith [h.lt_b] : 0 < a - b + c) (by linarith [h.lt_c] : 0 < a - c + b)]

/-- **Law of cosines**: `cos(angle opposite a) = (b² + c² − a²)/(2bc)`.

[formula: corner-angle] -/
theorem cos_cornerAngle (h : IsTriangle a b c) :
    cos (cornerAngle a b c) = (b ^ 2 + c ^ 2 - a ^ 2) / (2 * b * c) :=
  cos_arccos (cos_arg_bounds h).1.le (cos_arg_bounds h).2.le

theorem cornerAngle_pos (h : IsTriangle a b c) : 0 < cornerAngle a b c :=
  arccos_pos.mpr (cos_arg_bounds h).2

theorem cornerAngle_lt_pi (h : IsTriangle a b c) : cornerAngle a b c < π :=
  arccos_lt_pi.mpr (cos_arg_bounds h).1

/-- `sin(angle opposite a) = 2·Area/(bc)`, i.e. `Area = ½·b·c·sin A`.

[formula: triangle-area] -/
theorem sin_cornerAngle (h : IsTriangle a b c) :
    sin (cornerAngle a b c) = 2 * triArea a b c / (b * c) := by
  have hb := h.pos_b
  have hc := h.pos_c
  have hH := h.heron_pos
  unfold cornerAngle triArea
  rw [sin_arccos]
  have key : 1 - ((b ^ 2 + c ^ 2 - a ^ 2) / (2 * b * c)) ^ 2 =
      ((a + b + c) * (b - a + c) * (a - b + c) * (a - c + b)) / (2 * b * c) ^ 2 := by
    field_simp
    ring
  rw [key, sqrt_div hH.le, sqrt_sq (by positivity)]
  field_simp
  ring

/-- **Cotangent from edge lengths**: `cornerCot a b c = (b² + c² − a²)/(4·Area)` is the cotangent
of the corner angle opposite `a`.

[formula: corner-cot] -/
theorem cornerCot_eq_cot (h : IsTriangle a b c) :
    cornerCot a b c = cot (cornerAngle a b c) := by
  have hb := h.pos_b
  have hc := h.pos_c
  have hT := triArea_pos h
  rw [cot_eq_cos_div_sin, cos_cornerAngle h, sin_cornerAngle h, cornerCot]
  field_simp
  ring

/-- **Angle sum**: the three corner angles computed from the edge lengths add up to `π`.

[formula: triangle-angle-sum] -/
theorem cornerAngle_add_cornerAngle_add_cornerAngle (h : IsTriangle a b c) :
    cornerAngle a b c + cornerAngle b c a + cornerAngle c a b = π := by
  have ha := h.pos_a
  have hb := h.pos_b
  have hc := h.pos_c
  have h' := h.rotate
  have h'' := h'.rotate
  have hT := sixteen_mul_triArea_sq h
  have hTpos := triArea_pos h
  have cA := cos_cornerAngle h
  have cB := cos_cornerAngle h'
  have cC := cos_cornerAngle h''
  have sA := sin_cornerAngle h
  have sB := sin_cornerAngle h'
  have sC := sin_cornerAngle h''
  rw [triArea_rotate] at sB
  rw [triArea_rotate, triArea_rotate] at sC
  set A := cornerAngle a b c
  set B := cornerAngle b c a
  set C := cornerAngle c a b
  set T := triArea a b c
  have hA0 : 0 < A := cornerAngle_pos h
  have hA1 : A < π := cornerAngle_lt_pi h
  have hB0 : 0 < B := cornerAngle_pos h'
  have hB1 : B < π := cornerAngle_lt_pi h'
  have hC0 : 0 < C := cornerAngle_pos h''
  have hC1 : C < π := cornerAngle_lt_pi h''
  have hcos : cos (A + B) = cos (π - C) := by
    rw [cos_add, cos_pi_sub, cA, cB, cC, sA, sB]
    field_simp
    linear_combination (-1 : ℝ) * hT
  have hsin : sin (A + B) = sin C := by
    rw [sin_add, cA, cB, sA, sB, sC]
    field_simp
    ring
  have hsinpos : 0 < sin (A + B) := by rw [hsin]; exact sin_pos_of_pos_of_lt_pi hC0 hC1
  have hAB : A + B < π := by
    by_contra hcon
    push Not at hcon
    have h1 : 0 ≤ sin (A + B - π) :=
      sin_nonneg_of_nonneg_of_le_pi (by linarith) (by linarith)
    rw [sin_sub_pi] at h1
    linarith
  have := injOn_cos ⟨by linarith, hAB.le⟩ ⟨by linarith, by linarith⟩ hcos
  linarith

/-- **Half-angle form of the corner angle** (the numerically stable form the library evaluates):
`angle opposite a = 2·arctan √(((a−b+c)(a+b−c)) / ((a+b+c)(b+c−a)))`.

[formula: corner-angle] -/
theorem cornerAngle_eq_two_arctan (h : IsTriangle a b c) :
    cornerAngle a b c =
      2 * arctan (√(((a - b + c) * (a - c + b)) / ((a + b + c) * (b - a + c)))) := by
  have hb := h.pos_b
  have hc := h.pos_c
  have h0 := cornerAngle_pos h
  have h1 := cornerAngle_lt_pi h
  have hcos := cos_cornerAngle h
  have f1 : 0 < a + b + c := by linarith [h.pos_a]
  have f2 : 0 < b - a + c := by linarith [h.lt_a]
  have f3 : 0 < a - b + c := by linarith [h.lt_b]
  have f4 : 0 < a - c + b := by linarith [h.lt_c]
  set θ := cornerAngle a b c
  have hs : sin (θ / 2) = √((1 - cos θ) / 2) := sin_half_eq_sqrt h0.le (by linarith [pi_pos])
  have hc' : cos (θ / 2) = √((1 + cos θ) / 2) := cos_half (by linarith [pi_pos]) h1.le
  have e1 : 1 - cos θ = (a - b + c) * (a - c + b) / (2 * b * c) := by
    rw [hcos]; field_simp; ring
  have e2 : 1 + cos θ = (a + b + c) * (b - a + c) / (2 * b * c) := by
    rw [hcos]; field_simp; ring
  have htan : tan (θ / 2) = √(((a - b + c) * (a - c + b)) / ((a + b + c) * (b - a + c))) := by
    rw [tan_eq_sin_div_cos, hs, hc', ← sqrt_div (by rw [e1]; positivity), e1, e2]
    congr 1
    field_simp
  rw [← htan, arctan_tan (by linarith [pi_pos]) (by linarith)]
  ring

/-! ## Scaling -/

/-- Corner angles are invariant under uniform scaling of the triangle.

[formula: scale-invariance] -/
theorem cornerAngle_smul {s : ℝ} (hs : s ≠ 0) (a b c : ℝ) :
    cornerAngle (s * a) (s * b) (s * c) = cornerAngle a b c := by
  unfold cornerAngle
  congr 1
  rcases eq_or_ne b 0 with hb | hb
  · simp [hb]
  rcases eq_or_ne c 0 with hc | hc
  · simp [hc]
  field_simp

/-- The area scales with the square of the scaling factor.

[formula: scale-invariance] -/
theorem triArea_smul {s : ℝ} (a b c : ℝ) :
    triArea (s * a) (s * b) (s * c) = s ^ 2 * triArea a b c := by
  unfold triArea
  have : (s * a + s * b + s * c) * (s * b - s * a + s * c) * (s * a - s * b + s * c) *
      (s * a - s * c + s * b) =
      (s ^ 2) ^ 2 * ((a + b + c) * (b - a + c) * (a - b + c) * (a - c + b)) := by ring
  rw [this, sqrt_mul (sq_nonneg _), sqrt_sq (sq_nonneg s)]
  ring

/-- Corner cotangents are invariant under uniform scaling of the triangle.

[formula: scale-invariance] -/
theorem cornerCot_smul {s : ℝ} (hs : s ≠ 0) (a b c : ℝ) :
    cornerCot (s * a) (s * b) (s * c) = cornerCot a b c := by
  unfold cornerCot
  rw [triArea_smul]
  rcases eq_or_ne (triArea a b c) 0 with hT | hT
  · simp [hT]
  field_simp

end

/-! ## Link with Euclidean geometry (Mathlib) -/

section Euclidean

open EuclideanGeometry

variable {V P : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [MetricSpace P]
  [NormedAddTorsor V P]

/-- **The corner angle is the Euclidean angle.** For points `p₁ p₂ p₃` of a Euclidean space with
`p₁ ≠ p₂`, `p₃ ≠ p₂`, the angle at `p₂` equals `cornerAngle` of the three distances, the first
argument being the side opposite `p₂`.

[formula: corner-angle] -/
theorem angle_eq_cornerAngle {p₁ p₂ p₃ : P} (h12 : p₁ ≠ p₂) (h32 : p₃ ≠ p₂) :
    ∠ p₁ p₂ p₃ = cornerAngle (dist p₁ p₃) (dist p₁ p₂) (dist p₃ p₂) := by
  unfold cornerAngle
  have hl := law_cos p₁ p₂ p₃
  have d1 : 0 < dist p₁ p₂ := dist_pos.mpr h12
  have d2 : 0 < dist p₃ p₂ := dist_pos.mpr h32
  have hcos : cos (∠ p₁ p₂ p₃) =
      (dist p₁ p₂ ^ 2 + dist p₃ p₂ ^ 2 - dist p₁ p₃ ^ 2) / (2 * dist p₁ p₂ * dist p₃ p₂) := by
    field_simp
    linarith
  rw [← hcos, arccos_cos (angle_nonneg _ _ _) (angle_le_pi _ _ _)]

/-- **Heron's formula.** For points `p₁ p₂ p₃` of a Euclidean space, the area
`½·|p₁p₂|·|p₃p₂|·sin(∠ p₁ p₂ p₃)` of the triangle equals `triArea` of the three distances.

[formula: triangle-area] -/
theorem half_mul_sin_angle_eq_triArea (p₁ p₂ p₃ : P) :
    1 / 2 * dist p₁ p₂ * dist p₃ p₂ * sin (∠ p₁ p₂ p₃) =
      triArea (dist p₁ p₃) (dist p₁ p₂) (dist p₃ p₂) := by
  have hl := law_cos p₁ p₂ p₃
  have hsin : 0 ≤ sin (∠ p₁ p₂ p₃) := sin_nonneg_of_nonneg_of_le_pi (angle_nonneg _ _ _)
    (angle_le_pi _ _ _)
  have hsc := sin_sq_add_cos_sq (∠ p₁ p₂ p₃)
  set x := dist p₁ p₂
  set y := dist p₃ p₂
  set z := dist p₁ p₃
  set θ := ∠ p₁ p₂ p₃
  have hx : 0 ≤ x := dist_nonneg
  have hy : 0 ≤ y := dist_nonneg
  have key : (z + x + y) * (x - z + y) * (z - x + y) * (z - y + x) =
      (2 * x * y * sin θ) ^ 2 := by
    have h1 : 2 * x * y * cos θ = x ^ 2 + y ^ 2 - z ^ 2 := by nlinarith
    have h2 : (2 * x * y * sin θ) ^ 2 = (2 * x * y) ^ 2 - (2 * x * y * cos θ) ^ 2 := by
      have : sin θ ^ 2 = 1 - cos θ ^ 2 := by linarith
      calc (2 * x * y * sin θ) ^ 2 = (2 * x * y) ^ 2 * sin θ ^ 2 := by ring
        _ = (2 * x * y) ^ 2 * (1 - cos θ ^ 2) := by rw [this]
        _ = (2 * x * y) ^ 2 - (2 * x * y * cos θ) ^ 2 := by ring
    rw [h2, h1]
    ring
  unfold triArea
  rw [key, sqrt_sq (by positivity)]
  ring

end Euclidean

end DdgProofs
