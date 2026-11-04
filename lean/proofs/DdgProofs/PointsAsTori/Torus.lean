import DdgProofs.PointsAsTori.Fitting

/-!
# Points as Tori: the torus SDF (Eq. 23–24; SPEC §7.2–7.3) and the contact of the fitted torus

For a center `c`, a unit axis `u`, a major radius `R` and a signed minor radius `r`,
`g(x) = sgn(r)(|d| − |r|)`, `d = (v − R, z)`, `v = |(x − c) × u|`, `z = (x − c)·u` (`torusSDF`).

* **Exactness** (`torus_sdf`): `|g(x)| ≤ |x − y|` for every point `y` of the torus
  `{(|(y − c) × u| − R)² + ((y − c)·u)² = r²}`, for all `R` (also spindle tori), and for `R ≥ |r|`
  the bound is attained: `|g|` is the Euclidean distance to the torus. The proof is the
  "circle twice" argument: `|x − y|` is at least the distance of the meridian points
  `(v_x, z_x)`, `(v_y, z_y)` (`dist_ge_meridian`), which is at least the distance of `(v_x, z_x)`
  to the circle of radius `|r|` around `(R, 0)` (`sqrt2_tri`).
* **Gradient** (`torus_sdf_gradient`): `∇g = sgn(r)(w₁ (x − c − z u)/v + w₂ u)`, `w = d/|d|`,
  off the axis (`v ≠ 0`) and the core circle (`d ≠ 0`).
* **Stable form** (`torus_sdf_stable`): `v − R = (v² − R²)/(v + R)` with `v² − R²` expanded
  around a point `o` with `(o − c) ⊥ u`.
* **Contact** (`torus_contact_value`, `torus_contact_normal`, `torus_contact_curvatures`): see
  the docstrings; the curvature statement is where the paper's center is wrong for saddles.
-/

noncomputable section

namespace DdgProofs

namespace PAT

open Matrix Real
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff vnorm_smul
  dot_self_nonneg)

/-! ## Inequalities in `ℝ³` and `ℝ²` -/

/-- Cauchy–Schwarz for `vnorm`. -/
theorem dot_le_vnorm_mul (A B : V3) : A ⬝ᵥ B ≤ vnorm A * vnorm B := by
  have hL := cross_dot_cross A B A B
  have h0 := dot_self_nonneg (A ⨯₃ B)
  have hsq : (A ⬝ᵥ B) ^ 2 ≤ (vnorm A * vnorm B) ^ 2 := by
    rw [mul_pow, vnorm_sq, vnorm_sq]; rw [dotProduct_comm B A] at hL; nlinarith
  exact le_of_sq_le_sq hsq (mul_nonneg (vnorm_nonneg A) (vnorm_nonneg B))

/-- `|X × u|² = X·X − (X·u)²` for a unit `u`. -/
theorem cross_unit_sq {u : V3} (hu : u ⬝ᵥ u = 1) (X : V3) :
    (X ⨯₃ u) ⬝ᵥ (X ⨯₃ u) = X ⬝ᵥ X - (X ⬝ᵥ u) ^ 2 := by
  rw [cross_dot_cross, hu, dotProduct_comm u X]; ring

/-- `(X × u)·(Y × u) = X·Y − (X·u)(Y·u)` for a unit `u`. -/
theorem cross_unit_dot {u : V3} (hu : u ⬝ᵥ u = 1) (X Y : V3) :
    (X ⨯₃ u) ⬝ᵥ (Y ⨯₃ u) = X ⬝ᵥ Y - (X ⬝ᵥ u) * (Y ⬝ᵥ u) := by
  rw [cross_dot_cross, hu, dotProduct_comm u Y]; ring

/-- The distance of two points is at least the distance of their meridian coordinates
`(|· × u|, ·⬝u)`. -/
theorem dist_ge_meridian {u : V3} (hu : u ⬝ᵥ u = 1) (X Y : V3) :
    (vnorm (X ⨯₃ u) - vnorm (Y ⨯₃ u)) ^ 2 + (X ⬝ᵥ u - Y ⬝ᵥ u) ^ 2 ≤ (X - Y) ⬝ᵥ (X - Y) := by
  have hX := cross_unit_sq hu X
  have hY := cross_unit_sq hu Y
  have hXY := cross_unit_dot hu X Y
  have hcs := dot_le_vnorm_mul (X ⨯₃ u) (Y ⨯₃ u)
  have e : (X - Y) ⬝ᵥ (X - Y) = X ⬝ᵥ X - 2 * (X ⬝ᵥ Y) + Y ⬝ᵥ Y := by
    simp only [sub_dotProduct, dotProduct_sub, dotProduct_comm Y X]; ring
  rw [e]
  have h1 := vnorm_sq (X ⨯₃ u)
  have h2 := vnorm_sq (Y ⨯₃ u)
  nlinarith

/-- Triangle inequality in the plane: `|(a, c)| ≤ |(a, c) − (b, d)| + |(b, d)|`. -/
theorem sqrt2_tri (a c b d : ℝ) :
    √(a ^ 2 + c ^ 2) ≤ √((a - b) ^ 2 + (c - d) ^ 2) + √(b ^ 2 + d ^ 2) := by
  set p := √((a - b) ^ 2 + (c - d) ^ 2)
  set q := √(b ^ 2 + d ^ 2)
  have hp : p ^ 2 = (a - b) ^ 2 + (c - d) ^ 2 := Real.sq_sqrt (by positivity)
  have hq : q ^ 2 = b ^ 2 + d ^ 2 := Real.sq_sqrt (by positivity)
  have hp0 : 0 ≤ p := Real.sqrt_nonneg _
  have hq0 : 0 ≤ q := Real.sqrt_nonneg _
  have hcs : (a - b) * b + (c - d) * d ≤ p * q := by
    have hsq : ((a - b) * b + (c - d) * d) ^ 2 ≤ (p * q) ^ 2 := by
      rw [mul_pow, hp, hq]; nlinarith [sq_nonneg ((a - b) * d - (c - d) * b)]
    exact le_of_sq_le_sq hsq (mul_nonneg hp0 hq0)
  rw [Real.sqrt_le_iff]
  refine ⟨by positivity, ?_⟩
  nlinarith

/-! ## The torus SDF and its exactness -/

/-- The torus `{y : (|(y − c) × u| − R)² + ((y − c)·u)² = r²}`. -/
def torusSet (c u : V3) (R r : ℝ) : Set V3 :=
  {y | (vnorm ((y - c) ⨯₃ u) - R) ^ 2 + ((y - c) ⬝ᵥ u) ^ 2 = r ^ 2}

/-- `|d| = √((v − R)² + z²)`, the distance to the core circle. -/
def coreDist (c u : V3) (R : ℝ) (x : V3) : ℝ :=
  √((vnorm ((x - c) ⨯₃ u) - R) ^ 2 + ((x - c) ⬝ᵥ u) ^ 2)

theorem torusSDF_eq (c u : V3) (R r : ℝ) (x : V3) :
    torusSDF c u R r x = Real.sign r * (coreDist c u R x - |r|) := rfl

theorem abs_sign_le (r : ℝ) : |Real.sign r| ≤ 1 := by
  rcases lt_trichotomy r 0 with h | h | h
  · rw [Real.sign_of_neg h]; norm_num
  · rw [h, Real.sign_zero]; norm_num
  · rw [Real.sign_of_pos h]; norm_num

theorem abs_sign_of_ne {r : ℝ} (h : r ≠ 0) : |Real.sign r| = 1 := by
  rcases lt_or_gt_of_ne h with h | h
  · rw [Real.sign_of_neg h]; norm_num
  · rw [Real.sign_of_pos h]; norm_num

/-- **The torus SDF is the distance to the torus** (Eq. 23, "a circle, twice"): for every point
`y` of the torus, `|g(x)| ≤ |x − y|` (any `R`, spindle tori included); for `R ≥ |r|` and `r ≠ 0`
some point of the torus attains it, so `|g(x)|` is the Euclidean distance from `x` to the torus.
In general `||d| − |r||` is attained (`|g| = ||d| − |r||` for `r ≠ 0`).

[formula: pat-torus-sdf] -/
theorem torus_sdf {c u : V3} (hu : u ⬝ᵥ u = 1) (R r : ℝ) (x : V3) :
    (∀ y ∈ torusSet c u R r, |torusSDF c u R r x| ≤ vnorm (x - y)) ∧
    (|r| ≤ R → ∃ y ∈ torusSet c u R r, vnorm (x - y) = abs (coreDist c u R x - |r|)) ∧
    (r ≠ 0 → |torusSDF c u R r x| = abs (coreDist c u R x - |r|)) := by
  have habs : |torusSDF c u R r x| ≤ abs (coreDist c u R x - |r|) := by
    rw [torusSDF_eq, abs_mul]
    exact mul_le_of_le_one_left (abs_nonneg _) (abs_sign_le r)
  refine ⟨?_, ?_, ?_⟩
  · intro y hy
    refine habs.trans ?_
    set X := x - c
    set Y := y - c
    have hxy : x - y = X - Y := by simp [X, Y]
    have hy' : √((vnorm (Y ⨯₃ u) - R) ^ 2 + (Y ⬝ᵥ u) ^ 2) = |r| := by
      rw [show (vnorm (Y ⨯₃ u) - R) ^ 2 + (Y ⬝ᵥ u) ^ 2 = r ^ 2 from hy, Real.sqrt_sq_eq_abs]
    set D := √((vnorm (X ⨯₃ u) - vnorm (Y ⨯₃ u)) ^ 2 + (X ⬝ᵥ u - Y ⬝ᵥ u) ^ 2)
    have hD : D ≤ vnorm (x - y) := by
      rw [hxy]; exact Real.sqrt_le_sqrt (dist_ge_meridian hu X Y)
    have t1 := sqrt2_tri (vnorm (X ⨯₃ u) - R) (X ⬝ᵥ u) (vnorm (Y ⨯₃ u) - R) (Y ⬝ᵥ u)
    have t2 := sqrt2_tri (vnorm (Y ⨯₃ u) - R) (Y ⬝ᵥ u) (vnorm (X ⨯₃ u) - R) (X ⬝ᵥ u)
    have e1 : (vnorm (X ⨯₃ u) - R - (vnorm (Y ⨯₃ u) - R)) ^ 2 + (X ⬝ᵥ u - Y ⬝ᵥ u) ^ 2 =
        (vnorm (X ⨯₃ u) - vnorm (Y ⨯₃ u)) ^ 2 + (X ⬝ᵥ u - Y ⬝ᵥ u) ^ 2 := by ring
    have e2 : (vnorm (Y ⨯₃ u) - R - (vnorm (X ⨯₃ u) - R)) ^ 2 + (Y ⬝ᵥ u - X ⬝ᵥ u) ^ 2 =
        (vnorm (X ⨯₃ u) - vnorm (Y ⨯₃ u)) ^ 2 + (X ⬝ᵥ u - Y ⬝ᵥ u) ^ 2 := by ring
    rw [e1, hy'] at t1
    rw [e2, hy'] at t2
    change abs (coreDist c u R x - |r|) ≤ vnorm (x - y)
    unfold coreDist
    rw [abs_sub_le_iff]
    constructor <;> linarith
  · intro hR
    set X := x - c
    set v := vnorm (X ⨯₃ u)
    set z := X ⬝ᵥ u
    set q := coreDist c u R x
    have hq0 : 0 ≤ q := Real.sqrt_nonneg _
    have hq2 : q ^ 2 = (v - R) ^ 2 + z ^ 2 := Real.sq_sqrt (by positivity)
    -- a unit direction `e ⊥ u` with `X = v e + z u`
    obtain ⟨e, he, heu, hX⟩ : ∃ e : V3, e ⬝ᵥ e = 1 ∧ e ⬝ᵥ u = 0 ∧ X = v • e + z • u := by
      have hPX : (X - z • u) ⬝ᵥ (X - z • u) = v ^ 2 := by
        rw [vnorm_sq, cross_unit_sq hu]
        simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, hu,
          smul_eq_mul, dotProduct_comm u X]
        simp only [z]; ring
      by_cases hv : v = 0
      · refine ⟨frameS u, (refFrame_spec u hu).2.1, (refFrame_spec u hu).2.2.1, ?_⟩
        rw [hv, zero_smul, zero_add]
        rw [hv] at hPX
        have : X - z • u = 0 := dotProduct_self_eq_zero.mp (by rw [hPX]; ring)
        exact (sub_eq_zero.mp this)
      · refine ⟨v⁻¹ • (X - z • u), ?_, ?_, ?_⟩
        · rw [smul_dotProduct, dotProduct_smul, hPX, smul_eq_mul, smul_eq_mul]; field_simp
        · simp only [smul_dotProduct, sub_dotProduct, smul_eq_mul, hu, z]; ring
        · rw [smul_smul, mul_inv_cancel₀ hv, one_smul]; abel
    -- the unit vector `w = d/|d|` (or `(1, 0)` on the core circle)
    obtain ⟨w1, w2, hw, hw1, hd⟩ : ∃ w1 w2 : ℝ, w1 ^ 2 + w2 ^ 2 = 1 ∧ |w1| ≤ 1 ∧
        v - R = q * w1 ∧ z = q * w2 := by
      by_cases hq : q = 0
      · refine ⟨1, 0, by norm_num, by norm_num, ?_⟩
        have : (v - R) ^ 2 + z ^ 2 = 0 := by rw [← hq2, hq]; ring
        have h1 : v - R = 0 := by nlinarith [sq_nonneg (v - R), sq_nonneg z]
        have h2 : z = 0 := by nlinarith [sq_nonneg (v - R), sq_nonneg z]
        rw [hq, h1, h2]; simp
      · refine ⟨(v - R) / q, z / q, ?_, ?_, by field_simp, by field_simp⟩
        · rw [div_pow, div_pow, ← add_div, ← hq2]; field_simp
        · rw [abs_div, abs_of_nonneg hq0, div_le_one (lt_of_le_of_ne hq0 (Ne.symm hq))]
          exact abs_le_of_sq_le_sq (by rw [hq2]; nlinarith [sq_nonneg z]) hq0
    have hw1' : -1 ≤ w1 := (abs_le.mp hw1).1
    have hRw : 0 ≤ R + |r| * w1 := by nlinarith [abs_nonneg r]
    refine ⟨c + (R + |r| * w1) • e + (|r| * w2) • u, ?_, ?_⟩
    · -- on the torus
      show (vnorm ((c + (R + |r| * w1) • e + (|r| * w2) • u - c) ⨯₃ u) - R) ^ 2 +
        ((c + (R + |r| * w1) • e + (|r| * w2) • u - c) ⬝ᵥ u) ^ 2 = r ^ 2
      have hY : c + (R + |r| * w1) • e + (|r| * w2) • u - c = (R + |r| * w1) • e + (|r| * w2) • u := by
        abel
      rw [hY]
      have hcr : ((R + |r| * w1) • e + (|r| * w2) • u) ⨯₃ u = (R + |r| * w1) • (e ⨯₃ u) := by
        simp [map_add, map_smul, cross_self]
      have heu1 : vnorm (e ⨯₃ u) = 1 := by
        unfold vnorm; rw [cross_unit_sq hu, he, heu]; norm_num
      rw [hcr, vnorm_smul, heu1, mul_one, abs_of_nonneg hRw]
      simp only [add_dotProduct, smul_dotProduct, heu, hu, smul_eq_mul]
      have : r ^ 2 = |r| ^ 2 := (sq_abs r).symm
      rw [this]; linear_combination |r| ^ 2 * hw
    · -- the distance
      have hdiff : x - (c + (R + |r| * w1) • e + (|r| * w2) • u) =
          (v - R - |r| * w1) • e + (z - |r| * w2) • u := by
        have : x = c + X := by simp [X]
        rw [this, hX]; simp only [sub_smul]; abel_nf; simp only [sub_smul, add_smul]; abel
      rw [hdiff]
      have hn : ((v - R - |r| * w1) • e + (z - |r| * w2) • u) ⬝ᵥ
          ((v - R - |r| * w1) • e + (z - |r| * w2) • u) = (q - |r|) ^ 2 := by
        simp only [add_dotProduct, dotProduct_add, smul_dotProduct, dotProduct_smul, he, hu, heu,
          dotProduct_comm u e, smul_eq_mul]
        rw [hd.1, hd.2]
        linear_combination (q - |r|) ^ 2 * hw
      unfold vnorm; rw [hn, Real.sqrt_sq_eq_abs]
  · intro hr
    rw [torusSDF_eq, abs_mul, abs_sign_of_ne hr, one_mul]

/-! ## The gradient -/

theorem hasFDerivAt_crossRight (c u x : V3) :
    HasFDerivAt (fun y => (y - c) ⨯₃ u) (crossRight u) x := by
  have := (crossRight u).hasFDerivAt.comp x ((hasFDerivAt_id x).sub_const c)
  simpa using this

/-- `∇|(x − c) × u| = (x − c − ((x − c)·u) u)/|(x − c) × u|` off the axis. -/
theorem hasFDerivAt_axisDist {c u : V3} (hu : u ⬝ᵥ u = 1) {x : V3} (hv : (x - c) ⨯₃ u ≠ 0) :
    HasFDerivAt (fun y => vnorm ((y - c) ⨯₃ u))
      (dualOf ((vnorm ((x - c) ⨯₃ u))⁻¹ • ((x - c) - ((x - c) ⬝ᵥ u) • u))) x := by
  have h := hasFDerivAt_vnorm_comp (hasFDerivAt_crossRight c u x) hv
  refine h.congr_fderiv ?_
  ext h
  simp only [_root_.smul_apply, ContinuousLinearMap.comp_apply, dualOf_apply, crossRight_apply,
    smul_dotProduct, smul_eq_mul, cross_unit_dot hu, sub_dotProduct, dotProduct_comm h u]

/-- The torus SDF gradient `sgn(r)(w₁ (x − c − z u)/v + w₂ u)`, `w = d/|d|`. -/
def torusGrad (c u : V3) (R r : ℝ) (x : V3) : V3 :=
  Real.sign r • (((vnorm ((x - c) ⨯₃ u) - R) / (coreDist c u R x * vnorm ((x - c) ⨯₃ u))) •
      ((x - c) - ((x - c) ⬝ᵥ u) • u) + (((x - c) ⬝ᵥ u) / coreDist c u R x) • u)

/-- **The gradient of the torus SDF** (SPEC §7.2) off the axis (`v ≠ 0`) and the core circle
(`d ≠ 0`): `∇g = sgn(r)(w₁ (I − uuᵀ)(x − c)/v + w₂ u)`, `w = d/|d|`.

[formula: pat-torus-sdf-gradient] -/
theorem torus_sdf_gradient {c u : V3} (hu : u ⬝ᵥ u = 1) (R r : ℝ) {x : V3}
    (hv : (x - c) ⨯₃ u ≠ 0)
    (hd : (vnorm ((x - c) ⨯₃ u) - R) ^ 2 + ((x - c) ⬝ᵥ u) ^ 2 ≠ 0) :
    HasFDerivAt (torusSDF c u R r) (dualOf (torusGrad c u R r x)) x := by
  have h1 := ((hasFDerivAt_axisDist hu hv).sub_const R).pow 2
  have hz : HasFDerivAt (fun y => (y - c) ⬝ᵥ u) (dualOf u) x := by
    refine ((hasFDerivAt_dot_const u x).sub_const (u ⬝ᵥ c)).congr_of_eventuallyEq
      (Filter.Eventually.of_forall fun y => ?_)
    simp [sub_dotProduct, dotProduct_comm u]
  have h2 := hz.pow 2
  have h3 := ((h1.add h2).sqrt hd).sub_const |r| |>.const_mul (Real.sign r)
  have hvpos : 0 < vnorm ((x - c) ⨯₃ u) := vnorm_pos_iff.mpr hv
  have hqpos : 0 < coreDist c u R x := Real.sqrt_pos.mpr (lt_of_le_of_ne (by positivity) (Ne.symm hd))
  refine (h3.congr_fderiv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
  · ext h
    simp only [torusGrad, ContinuousLinearMap.smul_apply, _root_.smul_apply,
      ContinuousLinearMap.add_apply, ContinuousLinearMap.comp_apply, dualOf_apply,
      smul_dotProduct, add_dotProduct, smul_eq_mul, Function.comp, nsmul_eq_mul, Pi.add_apply]
    have hq : coreDist c u R x =
        √((vnorm ((x - c) ⨯₃ u) - R) ^ 2 + ((x - c) ⬝ᵥ u) ^ 2) := rfl
    rw [← hq, dotProduct_comm u h]
    set q := coreDist c u R x
    set v := vnorm ((x - c) ⨯₃ u)
    field_simp
    ring
  · simp [torusSDF, dotProduct_comm u]

end PAT

end DdgProofs
