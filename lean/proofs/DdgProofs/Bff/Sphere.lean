import Mathlib

/-!
# Mapping a sphere-like surface to the unit sphere: stereographic projection and inversions

Statements about the formulas evaluated by `crates/ddg-bff/src/sphere.rs` (Sawhney & Crane,
*Boundary First Flattening*, ACM TOG 2017, Sec. 7.5; the centering is Algorithm 1 of Baden, Crane,
Kazhdan, *Möbius Registration*, CGF 37(5), SGP 2018). Every Rust expression is transcribed as a
`def` over `ℝ` (exact arithmetic, not floating point). Points of space are triples
`ℝ × ℝ × ℝ`, written out coordinate by coordinate (no metric of the product type is used).

## Formula ids and results

* `stereographic-projection` — `inverse_stereographic`, `stereographic`:
  `invStereo_mem_sphere` (the image is on the unit sphere), `invStereo_ne_south` (it is not the
  south pole), `stereo_invStereo` (projecting back gives the point of the plane),
  `invStereo_stereo` (every point of the sphere but the south pole is an image).
* `stereographic-conformal` — `inverse_stereographic`: `hasDerivAt_invStereo_X`,
  `hasDerivAt_invStereo_Y` (the partial derivatives are `stereoJacX`, `stereoJacY`),
  `stereoJac_orthogonal`, `stereoJac_norm` (they are orthogonal and both of length
  `2 / (1 + X² + Y²)`: the differential is a similarity), `stereoJac_orientation` (the inner
  product of their cross product with the image point is `4 / (1 + X² + Y²)² > 0`: orientation is
  preserved, counter-clockwise seen from outside).
* `mobius-sphere-inversion` — `sphere_inversion`: `sphereInversion_mem_sphere` (a point of the
  unit sphere goes to the unit sphere when `x + c ≠ 0`), `normSqAdd_pos` (which holds for
  `|c| < 1`), `sphereInversion_zero` (the identity for `c = 0`).
* `mobius-centering-jacobian` — `mobius_center`: `hasDerivAt_sphereInversion` (moving the center
  along `t ↦ t u`, a point `x` of the sphere moves with velocity `2 (u − (u·x) x)` at `t = 0`,
  stated for the first coordinate; the others are the same with the axes renamed; hence the
  Jacobian of `Σ A_f ψ_c(C_f)` at `c = 0` is `2 Σ A_f (I − C_f C_fᵀ)` for points `C_f` of the
  sphere).

Cited in `CITED.md`: `mobius-centering-existence` (a centering inversion exists and is unique up
to rotation), `spherical-triangle-area` (Van Oosterom–Strackee), `sphere-injectivity-criterion`.
-/

namespace DDG

/-! ## Stereographic projection -/

/-- `inverse_stereographic`: with `s = X·X + Y·Y` and `d = 1 + s`, the point
`(2X/d, 2Y/d, (1 − s)/d)`. -/
noncomputable def invStereo (X Y : ℝ) : ℝ × ℝ × ℝ :=
  (2 * X / (1 + (X * X + Y * Y)), 2 * Y / (1 + (X * X + Y * Y)),
    (1 - (X * X + Y * Y)) / (1 + (X * X + Y * Y)))

/-- `stereographic`: `(x, y, z) ↦ (x/(1 + z), y/(1 + z))`. -/
noncomputable def stereo (p : ℝ × ℝ × ℝ) : ℝ × ℝ :=
  (p.1 / (1 + p.2.2), p.2.1 / (1 + p.2.2))

theorem one_add_normSq_pos (X Y : ℝ) : 0 < 1 + (X * X + Y * Y) := by
  nlinarith [mul_self_nonneg X, mul_self_nonneg Y]

/-- The inverse stereographic projection lands on the unit sphere.

[formula: stereographic-projection] -/
theorem invStereo_mem_sphere (X Y : ℝ) :
    (invStereo X Y).1 ^ 2 + (invStereo X Y).2.1 ^ 2 + (invStereo X Y).2.2 ^ 2 = 1 := by
  have hd := (one_add_normSq_pos X Y).ne'
  simp only [invStereo]
  field_simp
  ring

theorem one_add_invStereo_z (X Y : ℝ) :
    1 + (1 - (X * X + Y * Y)) / (1 + (X * X + Y * Y)) = 2 / (1 + (X * X + Y * Y)) := by
  have hd := (one_add_normSq_pos X Y).ne'
  field_simp
  ring

/-- The image is never the south pole `(0, 0, −1)` (the image of the removed vertex).

[formula: stereographic-projection] -/
theorem invStereo_ne_south (X Y : ℝ) : 1 + (invStereo X Y).2.2 ≠ 0 := by
  have hd := one_add_normSq_pos X Y
  simp only [invStereo, one_add_invStereo_z]
  positivity

/-- Projecting back gives the point of the plane: `stereographic ∘ inverse_stereographic = id`.

[formula: stereographic-projection] -/
theorem stereo_invStereo (X Y : ℝ) : stereo (invStereo X Y) = (X, Y) := by
  have hd := (one_add_normSq_pos X Y).ne'
  simp only [stereo, invStereo, one_add_invStereo_z]
  refine Prod.ext ?_ ?_
  · field_simp
  · field_simp

/-- Every point of the unit sphere other than the south pole is the image of a point of the plane:
`inverse_stereographic ∘ stereographic = id` there.

[formula: stereographic-projection] -/
theorem invStereo_stereo (x y z : ℝ) (h : x ^ 2 + y ^ 2 + z ^ 2 = 1) (hz : 1 + z ≠ 0) :
    invStereo (stereo (x, y, z)).1 (stereo (x, y, z)).2 = (x, y, z) := by
  simp only [stereo, invStereo]
  have hs : x / (1 + z) * (x / (1 + z)) + y / (1 + z) * (y / (1 + z)) = (1 - z) / (1 + z) := by
    field_simp
    linear_combination h
  rw [hs]
  have hd : 1 + (1 - z) / (1 + z) = 2 / (1 + z) := by field_simp; ring
  rw [hd]
  refine Prod.ext ?_ (Prod.ext ?_ ?_)
  · field_simp
  · field_simp
  · field_simp
    ring

/-! ## Conformality of the inverse stereographic projection -/

/-- The partial derivative of `invStereo` in `X` (quotient rule), with `d = 1 + X² + Y²`:
`(2(d − 2X²), −4XY, −4X) / d²`. -/
noncomputable def stereoJacX (X Y : ℝ) : ℝ × ℝ × ℝ :=
  ((2 * (1 + (X * X + Y * Y)) - 4 * X * X) / (1 + (X * X + Y * Y)) ^ 2,
    -4 * X * Y / (1 + (X * X + Y * Y)) ^ 2, -4 * X / (1 + (X * X + Y * Y)) ^ 2)

/-- The partial derivative of `invStereo` in `Y`: `(−4XY, 2(d − 2Y²), −4Y) / d²`. -/
noncomputable def stereoJacY (X Y : ℝ) : ℝ × ℝ × ℝ :=
  (-4 * X * Y / (1 + (X * X + Y * Y)) ^ 2,
    (2 * (1 + (X * X + Y * Y)) - 4 * Y * Y) / (1 + (X * X + Y * Y)) ^ 2,
    -4 * Y / (1 + (X * X + Y * Y)) ^ 2)

/-- `stereoJacX` is the derivative of `invStereo` in `X`, coordinate by coordinate.

[formula: stereographic-conformal] -/
theorem hasDerivAt_invStereo_X (X Y : ℝ) :
    HasDerivAt (fun t => (invStereo t Y).1) (stereoJacX X Y).1 X ∧
    HasDerivAt (fun t => (invStereo t Y).2.1) (stereoJacX X Y).2.1 X ∧
    HasDerivAt (fun t => (invStereo t Y).2.2) (stereoJacX X Y).2.2 X := by
  have hd := (one_add_normSq_pos X Y).ne'
  have hid : HasDerivAt (fun t : ℝ => t) 1 X := hasDerivAt_id' _
  have hden : HasDerivAt (fun t : ℝ => 1 + (t * t + Y * Y)) (1 * X + X * 1) X :=
    ((hid.fun_mul hid).add_const (Y * Y)).const_add 1
  have hnum1 : HasDerivAt (fun t : ℝ => 2 * t) (2 * 1) X := hid.const_mul 2
  have hnum3 : HasDerivAt (fun t : ℝ => 1 - (t * t + Y * Y)) (-(1 * X + X * 1)) X :=
    ((hid.fun_mul hid).add_const (Y * Y)).const_sub 1
  refine ⟨?_, ?_, ?_⟩
  · show HasDerivAt (fun t => 2 * t / (1 + (t * t + Y * Y))) _ X
    refine (hnum1.fun_div hden hd).congr_deriv ?_
    simp only [stereoJacX]
    field_simp
    ring
  · show HasDerivAt (fun t => 2 * Y / (1 + (t * t + Y * Y))) _ X
    refine ((hasDerivAt_const X (2 * Y)).fun_div hden hd).congr_deriv ?_
    simp only [stereoJacX]
    field_simp
    ring
  · show HasDerivAt (fun t => (1 - (t * t + Y * Y)) / (1 + (t * t + Y * Y))) _ X
    refine (hnum3.fun_div hden hd).congr_deriv ?_
    simp only [stereoJacX]
    field_simp
    ring

/-- `stereoJacY` is the derivative of `invStereo` in `Y`, coordinate by coordinate.

[formula: stereographic-conformal] -/
theorem hasDerivAt_invStereo_Y (X Y : ℝ) :
    HasDerivAt (fun t => (invStereo X t).1) (stereoJacY X Y).1 Y ∧
    HasDerivAt (fun t => (invStereo X t).2.1) (stereoJacY X Y).2.1 Y ∧
    HasDerivAt (fun t => (invStereo X t).2.2) (stereoJacY X Y).2.2 Y := by
  have hd := (one_add_normSq_pos X Y).ne'
  have hid : HasDerivAt (fun t : ℝ => t) 1 Y := hasDerivAt_id' _
  have hden : HasDerivAt (fun t : ℝ => 1 + (X * X + t * t)) (1 * Y + Y * 1) Y :=
    ((hid.fun_mul hid).const_add (X * X)).const_add 1
  have hnum2 : HasDerivAt (fun t : ℝ => 2 * t) (2 * 1) Y := hid.const_mul 2
  have hnum3 : HasDerivAt (fun t : ℝ => 1 - (X * X + t * t)) (-(1 * Y + Y * 1)) Y :=
    ((hid.fun_mul hid).const_add (X * X)).const_sub 1
  refine ⟨?_, ?_, ?_⟩
  · show HasDerivAt (fun t => 2 * X / (1 + (X * X + t * t))) _ Y
    refine ((hasDerivAt_const Y (2 * X)).fun_div hden hd).congr_deriv ?_
    simp only [stereoJacY]
    field_simp
    ring
  · show HasDerivAt (fun t => 2 * t / (1 + (X * X + t * t))) _ Y
    refine (hnum2.fun_div hden hd).congr_deriv ?_
    simp only [stereoJacY]
    field_simp
    ring
  · show HasDerivAt (fun t => (1 - (X * X + t * t)) / (1 + (X * X + t * t))) _ Y
    refine (hnum3.fun_div hden hd).congr_deriv ?_
    simp only [stereoJacY]
    field_simp
    ring

/-- The two partial derivatives are orthogonal.

[formula: stereographic-conformal] -/
theorem stereoJac_orthogonal (X Y : ℝ) :
    (stereoJacX X Y).1 * (stereoJacY X Y).1 + (stereoJacX X Y).2.1 * (stereoJacY X Y).2.1 +
      (stereoJacX X Y).2.2 * (stereoJacY X Y).2.2 = 0 := by
  have hd := (one_add_normSq_pos X Y).ne'
  simp only [stereoJacX, stereoJacY]
  field_simp
  ring

/-- Both partial derivatives have squared length `(2 / (1 + X² + Y²))²`: the differential is a
similarity with this factor.

[formula: stereographic-conformal] -/
theorem stereoJac_norm (X Y : ℝ) :
    (stereoJacX X Y).1 ^ 2 + (stereoJacX X Y).2.1 ^ 2 + (stereoJacX X Y).2.2 ^ 2 =
        (2 / (1 + (X * X + Y * Y))) ^ 2 ∧
      (stereoJacY X Y).1 ^ 2 + (stereoJacY X Y).2.1 ^ 2 + (stereoJacY X Y).2.2 ^ 2 =
        (2 / (1 + (X * X + Y * Y))) ^ 2 := by
  have hd := (one_add_normSq_pos X Y).ne'
  constructor
  · simp only [stereoJacX]
    field_simp
    ring
  · simp only [stereoJacY]
    field_simp
    ring

/-- The cross product of the partial derivatives is `4/d²` times... more precisely its inner
product with the image point is `4 / d²` (with `d = 1 + X² + Y²`), which is positive: the cross
product points away from the origin, so counter-clockwise in the plane is counter-clockwise
seen from outside the sphere.

[formula: stereographic-conformal] -/
theorem stereoJac_orientation (X Y : ℝ) :
    ((stereoJacX X Y).2.1 * (stereoJacY X Y).2.2 - (stereoJacX X Y).2.2 * (stereoJacY X Y).2.1) *
          (invStereo X Y).1 +
        ((stereoJacX X Y).2.2 * (stereoJacY X Y).1 - (stereoJacX X Y).1 * (stereoJacY X Y).2.2) *
          (invStereo X Y).2.1 +
        ((stereoJacX X Y).1 * (stereoJacY X Y).2.1 - (stereoJacX X Y).2.1 * (stereoJacY X Y).1) *
          (invStereo X Y).2.2 =
      4 / (1 + (X * X + Y * Y)) ^ 2 := by
  have hd := (one_add_normSq_pos X Y).ne'
  simp only [stereoJacX, stereoJacY, invStereo]
  field_simp
  ring

/-! ## Inversions of the sphere -/

/-- `sphere_inversion` with center `c = (c1, c2, c3)` at `x = (x1, x2, x3)`: `y = x + c`,
`y2 = y·y`, `s = (1 − c·c)/y2`, result `s y + c`. -/
noncomputable def sphereInversion (c1 c2 c3 x1 x2 x3 : ℝ) : ℝ × ℝ × ℝ :=
  ((1 - (c1 * c1 + c2 * c2 + c3 * c3)) /
        ((x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3)) * (x1 + c1) + c1,
    (1 - (c1 * c1 + c2 * c2 + c3 * c3)) /
        ((x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3)) * (x2 + c2) + c2,
    (1 - (c1 * c1 + c2 * c2 + c3 * c3)) /
        ((x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3)) * (x3 + c3) + c3)

/-- For `|c| < 1` and `|x| = 1`, `|x + c|² > 0`: the inversion is defined on the whole sphere.

[formula: mobius-sphere-inversion] -/
theorem normSqAdd_pos (c1 c2 c3 x1 x2 x3 : ℝ) (hx : x1 ^ 2 + x2 ^ 2 + x3 ^ 2 = 1)
    (hc : c1 ^ 2 + c2 ^ 2 + c3 ^ 2 < 1) :
    0 < (x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3) := by
  by_contra hle
  have s1 := mul_self_nonneg (x1 + c1)
  have s2 := mul_self_nonneg (x2 + c2)
  have s3 := mul_self_nonneg (x3 + c3)
  have e1 : x1 + c1 = 0 := mul_self_eq_zero.mp (by linarith)
  have e2 : x2 + c2 = 0 := mul_self_eq_zero.mp (by linarith)
  have e3 : x3 + c3 = 0 := mul_self_eq_zero.mp (by linarith)
  have : x1 = -c1 := by linarith
  have : x2 = -c2 := by linarith
  have : x3 = -c3 := by linarith
  subst x1 x2 x3
  nlinarith

/-- An inversion maps the unit sphere to itself: if `|x| = 1` and `x + c ≠ 0`, then
`|ψ_c(x)| = 1`. (Explicitly, `(|ψ|² − 1)|y|⁴ = |y|² (1 − |c|²) (1 − |x|²)` with `y = x + c`.)

[formula: mobius-sphere-inversion] -/
theorem sphereInversion_mem_sphere (c1 c2 c3 x1 x2 x3 : ℝ) (hx : x1 ^ 2 + x2 ^ 2 + x3 ^ 2 = 1)
    (hy : (x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3) ≠ 0) :
    (sphereInversion c1 c2 c3 x1 x2 x3).1 ^ 2 + (sphereInversion c1 c2 c3 x1 x2 x3).2.1 ^ 2 +
      (sphereInversion c1 c2 c3 x1 x2 x3).2.2 ^ 2 = 1 := by
  show ((1 - (c1 * c1 + c2 * c2 + c3 * c3)) /
        ((x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3)) * (x1 + c1) +
          c1) ^ 2 +
      ((1 - (c1 * c1 + c2 * c2 + c3 * c3)) /
        ((x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3)) * (x2 + c2) +
          c2) ^ 2 +
      ((1 - (c1 * c1 + c2 * c2 + c3 * c3)) /
        ((x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3)) * (x3 + c3) +
          c3) ^ 2 = 1
  generalize hD : (x1 + c1) * (x1 + c1) + (x2 + c2) * (x2 + c2) + (x3 + c3) * (x3 + c3) = D
    at hy ⊢
  have hk : D * D⁻¹ = 1 := mul_inv_cancel₀ hy
  simp only [div_eq_mul_inv]
  generalize D⁻¹ = k at hk ⊢
  generalize hS : 1 - (c1 * c1 + c2 * c2 + c3 * c3) = S
  linear_combination (S ^ 2 * k ^ 2 + S * k) * hD + (S ^ 2 * k + S) * hk - (S * k) * hx -
    (S * k + 1) * hS

/-- The inversion with center `0` is the identity on the unit sphere.

[formula: mobius-sphere-inversion] -/
theorem sphereInversion_zero (x1 x2 x3 : ℝ) (hx : x1 ^ 2 + x2 ^ 2 + x3 ^ 2 = 1) :
    sphereInversion 0 0 0 x1 x2 x3 = (x1, x2, x3) := by
  have h : x1 * x1 + x2 * x2 + x3 * x3 = 1 := by linear_combination hx
  simp only [sphereInversion, add_zero, mul_zero, sub_zero, h, div_one, one_mul]

/-- Moving the center of the inversion along `t ↦ t u`, a point `x` of the unit sphere moves,
at `t = 0`, with velocity `2 (u − (u·x) x)` (Baden et al. 2018, Sec. 4: `½ η̇ = u − ⟨u, x⟩ x`).
Summed with weights `A_f` over points `C_f` of the sphere this is the Jacobian
`2 Σ A_f (I − C_f C_fᵀ)` that `mobius_center` inverts. Stated for the first coordinate; the
others follow by renaming the axes (`sphereInversion` is symmetric under permuting them).

[formula: mobius-centering-jacobian] -/
theorem hasDerivAt_sphereInversion (u1 u2 u3 x1 x2 x3 : ℝ)
    (hx : x1 ^ 2 + x2 ^ 2 + x3 ^ 2 = 1) :
    HasDerivAt (fun t : ℝ => (sphereInversion (t * u1) (t * u2) (t * u3) x1 x2 x3).1)
      (2 * (u1 - (u1 * x1 + u2 * x2 + u3 * x3) * x1)) 0 := by
  have hid : HasDerivAt (fun t : ℝ => t) 1 0 := hasDerivAt_id' _
  have m1 : HasDerivAt (fun t : ℝ => t * u1) (1 * u1) 0 := hid.mul_const u1
  have m2 : HasDerivAt (fun t : ℝ => t * u2) (1 * u2) 0 := hid.mul_const u2
  have m3 : HasDerivAt (fun t : ℝ => t * u3) (1 * u3) 0 := hid.mul_const u3
  have l1 : HasDerivAt (fun t : ℝ => x1 + t * u1) (1 * u1) 0 := m1.const_add x1
  have l2 : HasDerivAt (fun t : ℝ => x2 + t * u2) (1 * u2) 0 := m2.const_add x2
  have l3 : HasDerivAt (fun t : ℝ => x3 + t * u3) (1 * u3) 0 := m3.const_add x3
  have hS := (((m1.fun_mul m1).fun_add (m2.fun_mul m2)).fun_add (m3.fun_mul m3)).const_sub 1
  have hD := ((l1.fun_mul l1).fun_add (l2.fun_mul l2)).fun_add (l3.fun_mul l3)
  have hD0 : (x1 + 0 * u1) * (x1 + 0 * u1) + (x2 + 0 * u2) * (x2 + 0 * u2) +
      (x3 + 0 * u3) * (x3 + 0 * u3) = 1 := by
    simp only [zero_mul, add_zero]; linear_combination hx
  have hq := hS.fun_div hD (by rw [hD0]; norm_num)
  show HasDerivAt (fun t : ℝ => (1 - (t * u1 * (t * u1) + t * u2 * (t * u2) + t * u3 * (t * u3))) /
      ((x1 + t * u1) * (x1 + t * u1) + (x2 + t * u2) * (x2 + t * u2) +
        (x3 + t * u3) * (x3 + t * u3)) * (x1 + t * u1) + t * u1) _ 0
  refine ((hq.fun_mul l1).fun_add m1).congr_deriv ?_
  rw [hD0]
  ring

end DDG
