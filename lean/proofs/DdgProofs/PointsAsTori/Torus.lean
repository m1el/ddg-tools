import DdgProofs.PointsAsTori.Fitting
import DdgProofs.PointsAsTori.Features

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

namespace DDG

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

/-! ## The stable form (D8) -/

/-- **Stable form of `v − R`** (SPEC §7.3): with a point `o` such that `m = o − c ⊥ u` and
`δ = x − o`, `v² − R² = |δ|² + 2 δ·m − (δ·u)² + (|m|² − R²)` and `v − R = (v² − R²)/(v + R)`
where `v + R ≠ 0`; the axial coordinate is `(x − c)·u = δ·u`. For the fitted torus
(`m = (η ρ) n*`, `|n*| = 1`, `η² = 1`, `R = ρ − σ r`) the constant is
`|m|² − R² = σ r (ρ + R)`, which has no cancellation. Hence Eq. 23 evaluated with these
expressions is the same function (`torusSDF_stable_eq`).

[formula: pat-torus-sdf-stable] -/
theorem torus_sdf_stable {c u o : V3} (hu : u ⬝ᵥ u = 1) (hm : (o - c) ⬝ᵥ u = 0) (R : ℝ) (x : V3)
    (hpos : vnorm ((x - c) ⨯₃ u) + R ≠ 0) :
    vnorm ((x - c) ⨯₃ u) - R =
      ((x - o) ⬝ᵥ (x - o) + 2 * ((x - o) ⬝ᵥ (o - c)) - ((x - o) ⬝ᵥ u) ^ 2 +
        ((o - c) ⬝ᵥ (o - c) - R ^ 2)) / (vnorm ((x - c) ⨯₃ u) + R) ∧
    (x - c) ⬝ᵥ u = (x - o) ⬝ᵥ u := by
  have hsplit : x - c = (x - o) + (o - c) := by abel
  have hz : (x - c) ⬝ᵥ u = (x - o) ⬝ᵥ u := by rw [hsplit, add_dotProduct, hm, add_zero]
  refine ⟨?_, hz⟩
  rw [eq_div_iff hpos]
  have hv2 : vnorm ((x - c) ⨯₃ u) ^ 2 = (x - c) ⬝ᵥ (x - c) - ((x - c) ⬝ᵥ u) ^ 2 := by
    rw [vnorm_sq, cross_unit_sq hu]
  have hexp : (x - c) ⬝ᵥ (x - c) = (x - o) ⬝ᵥ (x - o) + 2 * ((x - o) ⬝ᵥ (o - c)) +
      (o - c) ⬝ᵥ (o - c) := by
    rw [hsplit]; simp only [add_dotProduct, dotProduct_add, dotProduct_comm (o - c) (x - o)]; ring
  rw [hz] at hv2
  nlinarith [hv2, hexp]

/-- The fitted torus's constant: `|(ηρ) n*|² − (ρ − σ r)² = σ r (ρ + R)` for `σ² = 1`. -/
theorem stable_const {ns : V3} (hn : ns ⬝ᵥ ns = 1) {η ρ σ r : ℝ} (hη : η ^ 2 = 1) (hσ : σ ^ 2 = 1) :
    ((η * ρ) • ns) ⬝ᵥ ((η * ρ) • ns) - (ρ - σ * r) ^ 2 = σ * r * (ρ + (ρ - σ * r)) := by
  rw [smul_dotProduct, dotProduct_smul, hn, smul_eq_mul, smul_eq_mul]
  linear_combination ρ ^ 2 * hη

/-- Eq. 23 with the stable `v − R` (where `v + R ≠ 0`) is the torus SDF. -/
theorem torusSDF_stable_eq {c u o : V3} (hu : u ⬝ᵥ u = 1) (hm : (o - c) ⬝ᵥ u = 0) (R r : ℝ)
    (x : V3) (hpos : vnorm ((x - c) ⨯₃ u) + R ≠ 0) :
    torusSDF c u R r x = Real.sign r *
      (√((((x - o) ⬝ᵥ (x - o) + 2 * ((x - o) ⬝ᵥ (o - c)) - ((x - o) ⬝ᵥ u) ^ 2 +
        ((o - c) ⬝ᵥ (o - c) - R ^ 2)) / (vnorm ((x - c) ⨯₃ u) + R)) ^ 2 +
        ((x - o) ⬝ᵥ u) ^ 2) - |r|) := by
  obtain ⟨h1, h2⟩ := torus_sdf_stable hu hm R x hpos
  rw [torusSDF, h1, h2]

/-- **Stable form of `v − v_o`, any expansion point** (SPEC §7.3, generalizing `torus_sdf_stable`):
for any `o` (no condition `(o − c)·u = 0`), with `m = o − c`, `δ = x − o`,
`v² − v_o² = |δ|² + 2 δ·m − (δ·u)² − 2 (δ·u)(m·u)` and `v − v_o = (v² − v_o²)/(v + v_o)` where
`v + v_o ≠ 0`, with `v = |(x − c) × u|`, `v_o = |(o − c) × u|`; the axial coordinate is
`(x − c)·u = δ·u + m·u`.

[formula: pat-torus-sdf-stable] -/
theorem torus_sdf_stable_general {c u o : V3} (hu : u ⬝ᵥ u = 1) (x : V3)
    (hpos : vnorm ((x - c) ⨯₃ u) + vnorm ((o - c) ⨯₃ u) ≠ 0) :
    vnorm ((x - c) ⨯₃ u) - vnorm ((o - c) ⨯₃ u) =
      ((x - o) ⬝ᵥ (x - o) + 2 * ((x - o) ⬝ᵥ (o - c)) - ((x - o) ⬝ᵥ u) ^ 2 -
        2 * (((x - o) ⬝ᵥ u) * ((o - c) ⬝ᵥ u))) /
        (vnorm ((x - c) ⨯₃ u) + vnorm ((o - c) ⨯₃ u)) ∧
    (x - c) ⬝ᵥ u = (x - o) ⬝ᵥ u + (o - c) ⬝ᵥ u := by
  have hsplit : x - c = (x - o) + (o - c) := by abel
  have hz : (x - c) ⬝ᵥ u = (x - o) ⬝ᵥ u + (o - c) ⬝ᵥ u := by rw [hsplit, add_dotProduct]
  refine ⟨?_, hz⟩
  rw [eq_div_iff hpos]
  have hv2 : vnorm ((x - c) ⨯₃ u) ^ 2 = (x - c) ⬝ᵥ (x - c) - ((x - c) ⬝ᵥ u) ^ 2 := by
    rw [vnorm_sq, cross_unit_sq hu]
  have hw2 : vnorm ((o - c) ⨯₃ u) ^ 2 = (o - c) ⬝ᵥ (o - c) - ((o - c) ⬝ᵥ u) ^ 2 := by
    rw [vnorm_sq, cross_unit_sq hu]
  have hexp : (x - c) ⬝ᵥ (x - c) = (x - o) ⬝ᵥ (x - o) + 2 * ((x - o) ⬝ᵥ (o - c)) +
      (o - c) ⬝ᵥ (o - c) := by
    rw [hsplit]; simp only [add_dotProduct, dotProduct_add, dotProduct_comm (o - c) (x - o)]; ring
  rw [hz] at hv2
  nlinarith [hv2, hw2, hexp]

/-- **The second cancellation** (`‖d‖ − |r|`, SPEC §7.3): for `d = (d₁ + e₁, d₂ + e₂)`,
`√((d₁ + e₁)² + (d₂ + e₂)²) − |r| = (K + 2 (d₁ e₁ + d₂ e₂) + e₁² + e₂²)/(√(…) + |r|)` with
`K = d₁² + d₂² − r²` (computed once, without `O(r)` cancellation), where the denominator is
non-zero.

[formula: pat-torus-sdf-stable] -/
theorem stable_core_dist (d₁ d₂ e₁ e₂ r : ℝ)
    (h : √((d₁ + e₁) ^ 2 + (d₂ + e₂) ^ 2) + |r| ≠ 0) :
    √((d₁ + e₁) ^ 2 + (d₂ + e₂) ^ 2) - |r| =
      ((d₁ ^ 2 + d₂ ^ 2 - r ^ 2) + 2 * (d₁ * e₁ + d₂ * e₂) + (e₁ ^ 2 + e₂ ^ 2)) /
        (√((d₁ + e₁) ^ 2 + (d₂ + e₂) ^ 2) + |r|) := by
  rw [eq_div_iff h]
  have hs := Real.sq_sqrt (show 0 ≤ (d₁ + e₁) ^ 2 + (d₂ + e₂) ^ 2 by positivity)
  linear_combination hs - sq_abs r

/-- **The torus SDF in the doubly stable form** (the F32Stable kernel, SPEC §7.3): expanded around
any `o`, with `d₁ = v_o − R`, `d₂ = (o − c)·u`, `e₁ = (v² − v_o²)/(v + v_o)`, `e₂ = (x − o)·u`
and `K = d₁² + d₂² − r²`,
`g(x) = sgn(r) (K + 2 (d₁ e₁ + d₂ e₂) + e₁² + e₂²)/(√((d₁ + e₁)² + (d₂ + e₂)²) + |r|)`.
No assumption on `o`; the denominators are assumed non-zero.

[formula: pat-torus-sdf-stable] -/
theorem torusSDF_stable_general {c u o : V3} (hu : u ⬝ᵥ u = 1) (R r : ℝ) (x : V3)
    {d₁ d₂ e₁ e₂ : ℝ}
    (hd₁ : d₁ = vnorm ((o - c) ⨯₃ u) - R) (hd₂ : d₂ = (o - c) ⬝ᵥ u)
    (he₁ : e₁ = ((x - o) ⬝ᵥ (x - o) + 2 * ((x - o) ⬝ᵥ (o - c)) - ((x - o) ⬝ᵥ u) ^ 2 -
        2 * (((x - o) ⬝ᵥ u) * ((o - c) ⬝ᵥ u))) /
        (vnorm ((x - c) ⨯₃ u) + vnorm ((o - c) ⨯₃ u)))
    (he₂ : e₂ = (x - o) ⬝ᵥ u)
    (hpos : vnorm ((x - c) ⨯₃ u) + vnorm ((o - c) ⨯₃ u) ≠ 0)
    (hq : √((d₁ + e₁) ^ 2 + (d₂ + e₂) ^ 2) + |r| ≠ 0) :
    torusSDF c u R r x = Real.sign r *
      (((d₁ ^ 2 + d₂ ^ 2 - r ^ 2) + 2 * (d₁ * e₁ + d₂ * e₂) + (e₁ ^ 2 + e₂ ^ 2)) /
        (√((d₁ + e₁) ^ 2 + (d₂ + e₂) ^ 2) + |r|)) := by
  obtain ⟨h1, h2⟩ := torus_sdf_stable_general (c := c) (u := u) (o := o) hu x hpos
  have hv : vnorm ((x - c) ⨯₃ u) - R = d₁ + e₁ := by rw [hd₁, he₁, ← h1]; ring
  have hz : (x - c) ⬝ᵥ u = d₂ + e₂ := by rw [h2, hd₂, he₂]; ring
  rw [torusSDF_eq, coreDist, hv, hz, ← stable_core_dist d₁ d₂ e₁ e₂ r hq]

/-! ## Contact of the fitted torus with the patch -/

/-- **The fitted torus passes through `o = Q*(0,0)`** (SPEC §7.1, contact): with
`o = c + (η ρ) n*` (`η = ±1` the side of the center, the paper's `η = sign(T)`), `n* ⊥ u` unit,
`ρ = 1/|κmin| > 0`, `R = ρ − σ r` with `σ = sgn(κ₊κ₋) = ±1` and the stored minor radius `τ r`
(`τ = ±1`), `g(o) = 0`.

[formula: pat-torus-contact] -/
theorem torus_contact_value {c u ns : V3} (hu : u ⬝ᵥ u = 1) (hn : ns ⬝ᵥ ns = 1)
    (hun : ns ⬝ᵥ u = 0) {ρ r σ τ η : ℝ} (hρ : 0 < ρ) (hσ : σ ^ 2 = 1) (hτ : τ ^ 2 = 1)
    (hη : η ^ 2 = 1) :
    torusSDF c u (ρ - σ * r) (τ * r) (c + (η * ρ) • ns) = 0 := by
  have hv : vnorm ((c + (η * ρ) • ns - c) ⨯₃ u) = ρ := by
    rw [add_sub_cancel_left]
    unfold vnorm
    rw [cross_unit_sq hu, smul_dotProduct, dotProduct_smul, hn, smul_dotProduct, hun]
    simp only [smul_eq_mul, mul_one, mul_zero]
    rw [show (η * ρ) * (η * ρ) - 0 ^ 2 = ρ ^ 2 by linear_combination ρ ^ 2 * hη]
    exact Real.sqrt_sq hρ.le
  have hz : (c + (η * ρ) • ns - c) ⬝ᵥ u = 0 := by
    rw [add_sub_cancel_left, smul_dotProduct, hun, smul_zero]
  rw [torusSDF, hv, hz]
  have habs : |σ * r| = |r| := by
    rw [abs_mul, show |σ| = 1 from by nlinarith [sq_abs σ, abs_nonneg σ], one_mul]
  have habs' : |τ * r| = |r| := by
    rw [abs_mul, show |τ| = 1 from by nlinarith [sq_abs τ, abs_nonneg τ], one_mul]
  rw [show (ρ - (ρ - σ * r)) ^ 2 + (0 : ℝ) ^ 2 = (σ * r) ^ 2 by ring, Real.sqrt_sq_eq_abs,
    habs, habs', sub_self, mul_zero]

/-- **The contact fails for an exact zero product** `κ₊κ₋ = 0` under the D9 rule
(`sgn(0) = 0`, so `σ = 0` and `R = ρ`): then `o` is on the core circle and
`g(o) = −sgn(τ r)|r| ≠ 0` for `r ≠ 0`. -/
theorem torus_contact_fails_sign_zero {c u ns : V3} (hu : u ⬝ᵥ u = 1) (hn : ns ⬝ᵥ ns = 1)
    (hun : ns ⬝ᵥ u = 0) {ρ r τ η : ℝ} (hρ : 0 < ρ) (hη : η ^ 2 = 1) (hτ : τ ^ 2 = 1) (hr : r ≠ 0) :
    torusSDF c u (ρ - 0 * r) (τ * r) (c + (η * ρ) • ns) = -(Real.sign (τ * r) * |τ * r|) ∧
      torusSDF c u (ρ - 0 * r) (τ * r) (c + (η * ρ) • ns) ≠ 0 := by
  have hv : vnorm ((c + (η * ρ) • ns - c) ⨯₃ u) = ρ := by
    rw [add_sub_cancel_left]
    unfold vnorm
    rw [cross_unit_sq hu, smul_dotProduct, dotProduct_smul, hn, smul_dotProduct, hun]
    simp only [smul_eq_mul, mul_one, mul_zero]
    rw [show (η * ρ) * (η * ρ) - 0 ^ 2 = ρ ^ 2 by linear_combination ρ ^ 2 * hη]
    exact Real.sqrt_sq hρ.le
  have hz : (c + (η * ρ) • ns - c) ⬝ᵥ u = 0 := by
    rw [add_sub_cancel_left, smul_dotProduct, hun, smul_zero]
  have hval : torusSDF c u (ρ - 0 * r) (τ * r) (c + (η * ρ) • ns) =
      -(Real.sign (τ * r) * |τ * r|) := by
    rw [torusSDF, hv, hz]; simp
  refine ⟨hval, ?_⟩
  have hτ0 : τ ≠ 0 := by rintro rfl; norm_num at hτ
  have hne : τ * r ≠ 0 := mul_ne_zero hτ0 hr
  rw [hval, neg_ne_zero]
  exact mul_ne_zero (by
    rcases lt_or_gt_of_ne hne with h | h
    · rw [Real.sign_of_neg h]; norm_num
    · rw [Real.sign_of_pos h]; norm_num) (abs_ne_zero.mpr hne)

/-- **The torus normal at `o` is `±n*`** (SPEC §7.1, contact): `∇g(o) = (τ σ η) n*` for `r > 0`.
With the paper's center (`η = τ = sign(T)`) this is `σ n*`: the outward normal `+n*` for
`κ₊κ₋ > 0`, but `−n*` at saddles (`σ = −1`), where the SDF is negative on the outer side.

[formula: pat-torus-contact] -/
theorem torus_contact_normal {c u ns : V3} (hu : u ⬝ᵥ u = 1) (hn : ns ⬝ᵥ ns = 1)
    (hun : ns ⬝ᵥ u = 0) {ρ r σ τ η : ℝ} (hρ : 0 < ρ) (hr : 0 < r) (hσ : σ ^ 2 = 1)
    (hτ : τ ^ 2 = 1) (hη : η ^ 2 = 1) :
    HasFDerivAt (torusSDF c u (ρ - σ * r) (τ * r)) (dualOf ((τ * σ * η) • ns))
      (c + (η * ρ) • ns) := by
  set o := c + (η * ρ) • ns
  have hX : o - c = (η * ρ) • ns := add_sub_cancel_left c _
  have hv : vnorm ((o - c) ⨯₃ u) = ρ := by
    rw [hX]; unfold vnorm
    rw [cross_unit_sq hu, smul_dotProduct, dotProduct_smul, hn, smul_dotProduct, hun]
    simp only [smul_eq_mul, mul_one, mul_zero]
    rw [show (η * ρ) * (η * ρ) - 0 ^ 2 = ρ ^ 2 by linear_combination ρ ^ 2 * hη]
    exact Real.sqrt_sq hρ.le
  have hz : (o - c) ⬝ᵥ u = 0 := by rw [hX, smul_dotProduct, hun, smul_zero]
  have hτs : Real.sign (τ * r) = τ := by
    rcases mul_self_eq_one_iff.mp (by nlinarith : τ * τ = 1) with h | h <;> subst h
    · rw [one_mul, Real.sign_of_pos hr]
    · rw [Real.sign_of_neg (by linarith)]
  have habsσ : |σ| = 1 := by nlinarith [sq_abs σ, abs_nonneg σ]
  have hq : coreDist c u (ρ - σ * r) o = r := by
    unfold coreDist
    rw [hv, hz, show (ρ - (ρ - σ * r)) ^ 2 + (0 : ℝ) ^ 2 = (σ * r) ^ 2 by ring,
      Real.sqrt_sq_eq_abs, abs_mul, habsσ, one_mul, abs_of_pos hr]
  have hvne : (o - c) ⨯₃ u ≠ 0 := by
    intro h; have := congrArg vnorm h; rw [hv, RC.vnorm_zero] at this; linarith
  have hdne : (vnorm ((o - c) ⨯₃ u) - (ρ - σ * r)) ^ 2 + ((o - c) ⬝ᵥ u) ^ 2 ≠ 0 := by
    rw [hv, hz]; ring_nf; rw [hσ]; positivity
  refine (torus_sdf_gradient hu _ _ hvne hdne).congr_fderiv ?_
  congr 1
  rw [torusGrad, hq, hv, hz, hτs, hX]
  simp only [zero_smul, sub_zero, zero_div, add_zero, smul_smul]
  congr 1
  field_simp
  ring

/-! ### Curvatures of the torus at the contact point -/

/-- The torus through `c` with axis `u` in the frame `T = (e₁, u × e₁, u)`:
`X(θ, ψ) = c + (R + ρ cos θ)(cos ψ e₁ + sin ψ e₂) + ρ sin θ u` (signed tube radius `ρ`). -/
def torusParam (T : Frame) (c : V3) (R ρ θ ψ : ℝ) : V3 :=
  c + T.comb ((R + ρ * Real.cos θ) * Real.cos ψ) ((R + ρ * Real.cos θ) * Real.sin ψ) (ρ * Real.sin θ)

/-- The parametrization lies on the torus of `torus_sdf` (for `R + ρ cos θ ≥ 0`). -/
theorem torusParam_mem (T : Frame) (c : V3) {R ρ θ ψ : ℝ} (h : 0 ≤ R + ρ * Real.cos θ) :
    torusParam T c R ρ θ ψ ∈ torusSet c T.n R ρ := by
  show (vnorm ((torusParam T c R ρ θ ψ - c) ⨯₃ T.n) - R) ^ 2 +
    ((torusParam T c R ρ θ ψ - c) ⬝ᵥ T.n) ^ 2 = ρ ^ 2
  have hn : T.n = T.comb 0 0 1 := T.comb_n.symm
  simp only [torusParam, add_sub_cancel_left]
  rw [hn, T.comb_cross, T.vnorm_comb, T.comb_dot]
  have : √(((R + ρ * Real.cos θ) * Real.sin ψ * 1 - ρ * Real.sin θ * 0) ^ 2 +
      (ρ * Real.sin θ * 0 - (R + ρ * Real.cos θ) * Real.cos ψ * 1) ^ 2 +
      ((R + ρ * Real.cos θ) * Real.cos ψ * 0 - (R + ρ * Real.cos θ) * Real.sin ψ * 0) ^ 2) =
      R + ρ * Real.cos θ := by
    rw [show ((R + ρ * Real.cos θ) * Real.sin ψ * 1 - ρ * Real.sin θ * 0) ^ 2 +
      (ρ * Real.sin θ * 0 - (R + ρ * Real.cos θ) * Real.cos ψ * 1) ^ 2 +
      ((R + ρ * Real.cos θ) * Real.cos ψ * 0 - (R + ρ * Real.cos θ) * Real.sin ψ * 0) ^ 2 =
      (R + ρ * Real.cos θ) ^ 2 * (Real.sin ψ ^ 2 + Real.cos ψ ^ 2) by ring,
      Real.sin_sq_add_cos_sq, mul_one, Real.sqrt_sq h]
  rw [this]
  linear_combination ρ ^ 2 * Real.sin_sq_add_cos_sq θ

/-- **The fundamental forms of the torus at `θ = ψ = 0`.** `X(0,0) = c + (R + ρ) e₁`; the
partials there are `X_θ = ρ u`, `X_ψ = (R + ρ) e₂`, `X_θθ = −ρ e₁`, `X_ψψ = −(R + ρ) e₁`,
`X_θψ = 0`, so with the unit normal `N = η e₁` (`η = ±1`):
`E = ρ²`, `F = 0`, `G = (R + ρ)²`, `e = −ηρ`, `f = 0`, `g = −η(R + ρ)`; the principal
curvatures are `e/E = −η/ρ` (meridian, direction `u`) and `g/G = −η/(R + ρ)` (parallel,
direction `e₂`).

[formula: pat-torus-contact] -/
theorem torus_param_forms (T : Frame) (c : V3) (R ρ η : ℝ) :
    torusParam T c R ρ 0 0 = c + (R + ρ) • T.s ∧
    HasDerivAt (fun θ => torusParam T c R ρ θ 0) (T.comb 0 0 ρ) 0 ∧
    HasDerivAt (fun ψ => torusParam T c R ρ 0 ψ) (T.comb 0 (R + ρ) 0) 0 ∧
    (∀ θ, HasDerivAt (fun θ => torusParam T c R ρ θ 0)
      (T.comb (-(ρ * Real.sin θ)) 0 (ρ * Real.cos θ)) θ) ∧
    (∀ ψ, HasDerivAt (fun ψ => torusParam T c R ρ 0 ψ)
      (T.comb (-((R + ρ) * Real.sin ψ)) ((R + ρ) * Real.cos ψ) 0) ψ) ∧
    HasDerivAt (fun θ => T.comb (-(ρ * Real.sin θ)) 0 (ρ * Real.cos θ)) (T.comb (-ρ) 0 0) 0 ∧
    HasDerivAt (fun ψ => T.comb (-((R + ρ) * Real.sin ψ)) ((R + ρ) * Real.cos ψ) 0)
      (T.comb (-(R + ρ)) 0 0) 0 ∧
    (T.comb 0 0 ρ ⬝ᵥ T.comb 0 0 ρ = ρ ^ 2 ∧ T.comb 0 0 ρ ⬝ᵥ T.comb 0 (R + ρ) 0 = 0 ∧
      T.comb 0 (R + ρ) 0 ⬝ᵥ T.comb 0 (R + ρ) 0 = (R + ρ) ^ 2) ∧
    (T.comb η 0 0 ⬝ᵥ T.comb (-ρ) 0 0 = -(η * ρ) ∧
      T.comb η 0 0 ⬝ᵥ T.comb (-(R + ρ)) 0 0 = -(η * (R + ρ))) := by
  have hc : ∀ (k x : ℝ), HasDerivAt (fun _ : ℝ => k) 0 x := fun k x => hasDerivAt_const x k
  have hθ : ∀ θ, HasDerivAt (fun θ => torusParam T c R ρ θ 0)
      (T.comb (-(ρ * Real.sin θ)) 0 (ρ * Real.cos θ)) θ := by
    intro θ
    have h1 : HasDerivAt (fun θ => (R + ρ * Real.cos θ) * Real.cos 0) (-(ρ * Real.sin θ)) θ := by
      simpa using ((Real.hasDerivAt_cos θ).const_mul ρ).const_add R
    have h2 : HasDerivAt (fun θ => (R + ρ * Real.cos θ) * Real.sin 0) 0 θ := by
      simpa using hc 0 θ
    have h3 : HasDerivAt (fun θ => ρ * Real.sin θ) (ρ * Real.cos θ) θ :=
      (Real.hasDerivAt_sin θ).const_mul ρ
    exact (hasDerivAt_comb T h1 h2 h3).const_add c
  have hψ : ∀ ψ, HasDerivAt (fun ψ => torusParam T c R ρ 0 ψ)
      (T.comb (-((R + ρ) * Real.sin ψ)) ((R + ρ) * Real.cos ψ) 0) ψ := by
    intro ψ
    have h1 : HasDerivAt (fun ψ => (R + ρ * Real.cos 0) * Real.cos ψ) (-((R + ρ) * Real.sin ψ)) ψ := by
      simpa using (Real.hasDerivAt_cos ψ).const_mul (R + ρ)
    have h2 : HasDerivAt (fun ψ => (R + ρ * Real.cos 0) * Real.sin ψ) ((R + ρ) * Real.cos ψ) ψ := by
      simpa using (Real.hasDerivAt_sin ψ).const_mul (R + ρ)
    have h3 : HasDerivAt (fun _ : ℝ => ρ * Real.sin 0) 0 ψ := hc _ ψ
    exact (hasDerivAt_comb T h1 h2 h3).const_add c
  refine ⟨?_, by simpa using hθ 0, by simpa using hψ 0, hθ, hψ, ?_, ?_, ?_, ?_⟩
  · simp [torusParam, Frame.comb]
  · have h1 : HasDerivAt (fun θ => -(ρ * Real.sin θ)) (-(ρ * Real.cos 0)) 0 :=
      ((Real.hasDerivAt_sin 0).const_mul ρ).neg
    have h3 : HasDerivAt (fun θ => ρ * Real.cos θ) (ρ * -Real.sin 0) 0 :=
      (Real.hasDerivAt_cos 0).const_mul ρ
    simpa using hasDerivAt_comb T h1 (hc 0 0) h3
  · have h1 : HasDerivAt (fun ψ => -((R + ρ) * Real.sin ψ)) (-((R + ρ) * Real.cos 0)) 0 :=
      ((Real.hasDerivAt_sin 0).const_mul (R + ρ)).neg
    have h2 : HasDerivAt (fun ψ => (R + ρ) * Real.cos ψ) ((R + ρ) * -Real.sin 0) 0 :=
      (Real.hasDerivAt_cos 0).const_mul (R + ρ)
    simpa using hasDerivAt_comb T h1 h2 (hc 0 0)
  · simp only [T.comb_dot]; refine ⟨by ring, by ring, by ring⟩
  · simp only [T.comb_dot]; refine ⟨by ring, by ring⟩

/-- **The principal curvatures of the fitted torus at `o`** (SPEC §7.1, contact). The fitted
torus has `ρ = σ r` (signed tube radius), `R + ρ = 1/|κmin|`, `r = 1/|κmax|`, `σ = sgn(κmin κmax)`,
and its normal at `o` is `n* = η e₁` where `η` is the side of the center
(`c = o − (η/|κmin|) n*`). By `torus_param_forms` its curvatures w.r.t. `n*` are
`−η/(σ r) = −ησ|κmax|` (meridian, along `u ∥ v̂max`) and `−η/(R + ρ) = −η|κmin|` (parallel, along
`v̂min`). For nonzero `κmin, κmax`:

* the paper's center (`η = sign(T)`, Eq. 22) reproduces `(κmax, κmin)` iff `κmin κmax > 0`;
* for saddles (`κmin κmax < 0`) it gives `(−κmax, −κmin)`: the mirror image of the patch;
* the corrected center `η = −sgn(κmin)`, i.e. `c = o + n*/κmin`, reproduces `(κmax, κmin)` in
  all cases, and then `∇g(o) = sign(T)·σ·η n* = +n*` (`torus_contact_normal`), while the paper's
  center gives `σ n*`.

[formula: pat-torus-contact] -/
theorem torus_contact_curvatures {kmin kmax : ℝ} (h0 : kmin ≠ 0) (h1 : kmax ≠ 0) :
    (0 < kmin * kmax →
      -(torusSign kmin kmax) * Real.sign (kmin * kmax) * |kmax| = kmax ∧
        -(torusSign kmin kmax) * |kmin| = kmin) ∧
    (kmin * kmax < 0 →
      -(torusSign kmin kmax) * Real.sign (kmin * kmax) * |kmax| = -kmax ∧
        -(torusSign kmin kmax) * |kmin| = -kmin) ∧
    (-(-Real.sign kmin) * Real.sign (kmin * kmax) * |kmax| = kmax ∧
      -(-Real.sign kmin) * |kmin| = kmin) ∧
    torusSign kmin kmax * Real.sign (kmin * kmax) * (-Real.sign kmin) = 1 ∧
    torusSign kmin kmax * Real.sign (kmin * kmax) * torusSign kmin kmax = Real.sign (kmin * kmax) := by
  rcases lt_or_gt_of_ne h0 with ha | ha <;> rcases lt_or_gt_of_ne h1 with hb | hb
  · have hp : 0 < kmin * kmax := mul_pos_of_neg_of_neg ha hb
    simp only [torusSign, hb, ha, and_self, true_or, ite_true, Real.sign_of_pos hp,
      Real.sign_of_neg ha, abs_of_neg ha, abs_of_neg hb]
    refine ⟨fun _ => ⟨by ring, by ring⟩, fun h => absurd h (not_lt.mpr hp.le), ⟨by ring, by ring⟩,
      by ring, by ring⟩
  · have hp : kmin * kmax < 0 := mul_neg_of_neg_of_pos ha hb
    have : ¬ (kmax < 0) := not_lt.mpr hb.le
    simp only [torusSign, this, false_and, or_self, ite_false, Real.sign_of_neg hp,
      Real.sign_of_neg ha, abs_of_neg ha, abs_of_pos hb]
    refine ⟨fun h => absurd h (not_lt.mpr hp.le), fun _ => ⟨by ring, by ring⟩, ⟨by ring, by ring⟩,
      by ring, by ring⟩
  · have hp : kmin * kmax < 0 := mul_neg_of_pos_of_neg ha hb
    simp only [torusSign, hb, ha, true_and, not_lt.mpr ha.le, or_true, ite_true,
      Real.sign_of_neg hp, Real.sign_of_pos ha, abs_of_pos ha, abs_of_neg hb]
    refine ⟨fun h => absurd h (not_lt.mpr hp.le), fun _ => ⟨by ring, by ring⟩, ⟨by ring, by ring⟩,
      by ring, by ring⟩
  · have hp : 0 < kmin * kmax := mul_pos ha hb
    have : ¬ (kmax < 0) := not_lt.mpr hb.le
    simp only [torusSign, this, false_and, or_self, ite_false, Real.sign_of_pos hp,
      Real.sign_of_pos ha, abs_of_pos ha, abs_of_pos hb]
    refine ⟨fun _ => ⟨by ring, by ring⟩, fun h => absurd h (not_lt.mpr hp.le), ⟨by ring, by ring⟩,
      by ring, by ring⟩

/-- `sgn(x)² = 1` for `x ≠ 0`. -/
theorem sign_sq_of_ne {x : ℝ} (h : x ≠ 0) : Real.sign x ^ 2 = 1 := by
  rcases lt_or_gt_of_ne h with h | h
  · rw [Real.sign_of_neg h]; norm_num
  · rw [Real.sign_of_pos h]; norm_num

theorem torusSign_sq (kmin kmax : ℝ) : torusSign kmin kmax ^ 2 = 1 := by
  unfold torusSign; split_ifs <;> norm_num

/-- **Contact of the fitted torus under both sign rules** (SPEC §18, resolution 1;
`FitOptions::sign_rule`). Let `κmin, κmax ≠ 0`, `n*` a unit normal at `o = Q*(0,0)` orthogonal to
the axis `u`, `r = 1/|κmax|`, `σ = sgn(κmin κmax)`, `R = 1/|κmin| − σ r`, the stored minor radius
`sign(T)·r` (Eq. 22) and `g` the torus SDF (Eq. 23–24).

* `SignRule::Geometric`, center `c = o + n*/κmin`: `g(o) = 0` and `∇g(o) = n*`, so
  `∂g/∂n*(o) = +1` (the contact condition).
* `SignRule::Reference`, center `c = o − (sign(T)/|κmin|) n*` (the paper): `g(o) = 0` and
  `∇g(o) = σ n*`, so `∂g/∂n*(o) = −1` exactly at saddles (`κmin κmax < 0`) and `+1` otherwise.

[formula: pat-fit-contact] [formula: pat-torus-contact] -/
theorem fit_contact {u ns o : V3} (hu : u ⬝ᵥ u = 1) (hn : ns ⬝ᵥ ns = 1) (hun : ns ⬝ᵥ u = 0)
    {kmin kmax : ℝ} (h0 : kmin ≠ 0) (h1 : kmax ≠ 0) :
    let r := 1 / |kmax|
    let σ := Real.sign (kmin * kmax)
    let R := 1 / |kmin| - σ * r
    let τ := torusSign kmin kmax
    (torusSDF (o + (1 / kmin) • ns) u R (τ * r) o = 0 ∧
      HasFDerivAt (torusSDF (o + (1 / kmin) • ns) u R (τ * r)) (dualOf ns) o ∧
      dualOf ns ns = 1) ∧
    (torusSDF (o - (τ / |kmin|) • ns) u R (τ * r) o = 0 ∧
      HasFDerivAt (torusSDF (o - (τ / |kmin|) • ns) u R (τ * r)) (dualOf (σ • ns)) o ∧
      dualOf (σ • ns) ns = σ ∧ (kmin * kmax < 0 → σ = -1) ∧ (0 < kmin * kmax → σ = 1)) := by
  intro r σ R τ
  have hρ : 0 < 1 / |kmin| := by have := abs_pos.mpr h0; positivity
  have hr : 0 < r := by have := abs_pos.mpr h1; simp only [r]; positivity
  have hσ : σ ^ 2 = 1 := sign_sq_of_ne (mul_ne_zero h0 h1)
  have hτ : τ ^ 2 = 1 := torusSign_sq kmin kmax
  have hη : (-Real.sign kmin) ^ 2 = 1 := by rw [neg_sq]; exact sign_sq_of_ne h0
  obtain ⟨-, -, -, hgeo, href⟩ := torus_contact_curvatures h0 h1
  -- geometric center: `o = c + (η ρ) n*` with `η = −sgn κmin`
  have hoG : o = (o + (1 / kmin) • ns) + (-Real.sign kmin * (1 / |kmin|)) • ns := by
    have : -Real.sign kmin * (1 / |kmin|) = -(1 / kmin) := by
      rcases lt_or_gt_of_ne h0 with h | h
      · rw [Real.sign_of_neg h, abs_of_neg h]; field_simp
      · rw [Real.sign_of_pos h, abs_of_pos h]; field_simp
    rw [this, neg_smul]; abel
  have hoR : o = (o - (τ / |kmin|) • ns) + (τ * (1 / |kmin|)) • ns := by
    rw [show τ * (1 / |kmin|) = τ / |kmin| by ring]; abel
  have vG := torus_contact_value (c := o + (1 / kmin) • ns) hu hn hun hρ hσ hτ hη (r := r)
  have nG := torus_contact_normal (c := o + (1 / kmin) • ns) hu hn hun hρ hr hσ hτ hη
  have vR := torus_contact_value (c := o - (τ / |kmin|) • ns) hu hn hun hρ hσ hτ hτ (r := r)
  have nR := torus_contact_normal (c := o - (τ / |kmin|) • ns) hu hn hun hρ hr hσ hτ hτ
  rw [← hoG] at vG nG
  rw [← hoR] at vR nR
  have e1 : τ * σ * -Real.sign kmin = 1 := hgeo
  have e2 : τ * σ * τ = σ := href
  rw [e1, one_smul] at nG
  rw [e2] at nR
  refine ⟨⟨vG, nG, by rw [dualOf_apply, hn]⟩, ⟨vR, nR, by rw [dualOf_apply, smul_dotProduct, hn,
    smul_eq_mul, mul_one], fun h => Real.sign_of_neg h, fun h => Real.sign_of_pos h⟩⟩

/-! ## The Laplacian (ours, D24) -/

/-- The trace of a linear map of `ℝ³` (the divergence of a field with that derivative). -/
def trace3 (D : V3 →L[ℝ] V3) : ℝ := ∑ i : Fin 3, D (Pi.single i 1) i

theorem trace3_add (A B : V3 →L[ℝ] V3) : trace3 (A + B) = trace3 A + trace3 B := by
  simp [trace3, Finset.sum_add_distrib]

theorem trace3_smul (c : ℝ) (A : V3 →L[ℝ] V3) : trace3 (c • A) = c * trace3 A := by
  simp [trace3, Finset.mul_sum]

theorem trace3_sub (A B : V3 →L[ℝ] V3) : trace3 (A - B) = trace3 A - trace3 B := by
  simp [trace3, Finset.sum_sub_distrib]

theorem trace3_id : trace3 (ContinuousLinearMap.id ℝ V3) = 3 := by
  simp [trace3]

/-- `tr(h ↦ L(h) b) = L(b)`. -/
theorem trace3_smulRight (L : V3 →L[ℝ] ℝ) (b : V3) : trace3 (L.smulRight b) = L b := by
  have hb : b = b 0 • Pi.single 0 1 + b 1 • Pi.single 1 1 + b 2 • Pi.single 2 1 := by
    ext i; fin_cases i <;> simp
  conv_rhs => rw [hb]
  simp only [trace3, Fin.sum_univ_three, ContinuousLinearMap.smulRight_apply, map_add, map_smul,
    smul_eq_mul, Pi.smul_apply]
  ring

theorem hasFDerivAt_inv_comp {f : V3 → ℝ} {f' : V3 →L[ℝ] ℝ} {x : V3} (hf : HasFDerivAt f f' x)
    (h : f x ≠ 0) : HasFDerivAt (fun y => (f y)⁻¹) ((-(f x ^ 2)⁻¹) • f') x := by
  have := (hasFDerivAt_inv h).comp x hf
  refine this.congr_fderiv ?_
  ext v; simp [mul_comm]

theorem proj_sq {u : V3} (hu : u ⬝ᵥ u = 1) (X : V3) :
    (X - (X ⬝ᵥ u) • u) ⬝ᵥ (X - (X ⬝ᵥ u) • u) = vnorm (X ⨯₃ u) ^ 2 := by
  rw [vnorm_sq, cross_unit_sq hu]
  simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, hu, smul_eq_mul,
    dotProduct_comm u X]; ring

theorem proj_dot_u {u : V3} (hu : u ⬝ᵥ u = 1) (X : V3) : (X - (X ⬝ᵥ u) • u) ⬝ᵥ u = 0 := by
  simp only [sub_dotProduct, smul_dotProduct, hu, smul_eq_mul]; ring

/-- **The Laplacian of the torus SDF** (ours, replaces the reference's placeholder, D24): off
the axis and the core circle, the gradient field `torusGrad` (the gradient of `g` there,
`torus_sdf_gradient`) has a derivative `D` with
`Δg = tr D = sgn(r)(1/|d| + d₁/(|d| v))`, `d₁ = v − R`: the 2D Laplacian `1/|d|` of the
distance to the core point in the meridian plane plus the axisymmetric term `∂_v|d|/v`.

[formula: pat-torus-sdf-laplacian] -/
theorem torus_sdf_laplacian {c u : V3} (hu : u ⬝ᵥ u = 1) (R r : ℝ) {x : V3}
    (hv : (x - c) ⨯₃ u ≠ 0)
    (hd : (vnorm ((x - c) ⨯₃ u) - R) ^ 2 + ((x - c) ⬝ᵥ u) ^ 2 ≠ 0) :
    ∃ D : V3 →L[ℝ] V3, HasFDerivAt (torusGrad c u R r) D x ∧
      trace3 D = Real.sign r * (1 / coreDist c u R x +
        (vnorm ((x - c) ⨯₃ u) - R) / (coreDist c u R x * vnorm ((x - c) ⨯₃ u))) := by
  set v := vnorm ((x - c) ⨯₃ u) with hvdef
  set q := coreDist c u R x with hqdef
  set z := (x - c) ⬝ᵥ u with hzdef
  set P := (x - c) - z • u with hPdef
  have hvpos : 0 < v := vnorm_pos_iff.mpr hv
  have hqpos : 0 < q := Real.sqrt_pos.mpr (lt_of_le_of_ne (by positivity) (Ne.symm hd))
  have hq2 : q ^ 2 = (v - R) ^ 2 + z ^ 2 := Real.sq_sqrt (by positivity)
  have hPP : P ⬝ᵥ P = v ^ 2 := proj_sq hu (x - c)
  have hPu : P ⬝ᵥ u = 0 := proj_dot_u hu (x - c)
  -- derivatives of the pieces
  have hV := hasFDerivAt_axisDist hu hv
  have hZ : HasFDerivAt (fun y => (y - c) ⬝ᵥ u) (dualOf u) x := by
    refine ((hasFDerivAt_dot_const u x).sub_const (u ⬝ᵥ c)).congr_of_eventuallyEq
      (Filter.Eventually.of_forall fun y => ?_)
    simp [sub_dotProduct, dotProduct_comm u]
  have hQ : HasFDerivAt (fun y => coreDist c u R y) (dualOf (torusGrad c u R 1 x)) x := by
    have := (torus_sdf_gradient hu R 1 hv hd).add_const 1
    refine this.congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
    simp [torusSDF, coreDist]
  have hPd : HasFDerivAt (fun y => (y - c) - ((y - c) ⬝ᵥ u) • u)
      (ContinuousLinearMap.id ℝ V3 - (dualOf u).smulRight u) x :=
    ((hasFDerivAt_id x).sub_const c).sub (hZ.smul_const u)
  have hqv : HasFDerivAt (fun y => coreDist c u R y * vnorm ((y - c) ⨯₃ u))
      (q • dualOf (v⁻¹ • P) + v • dualOf (torusGrad c u R 1 x)) x := hQ.mul hV
  have hφ1 := (hV.sub_const R).mul (hasFDerivAt_inv_comp hqv (by positivity))
  have hφ2 := hZ.mul (hasFDerivAt_inv_comp hQ hqpos.ne')
  have hG := ((hφ1.smul hPd).add (hφ2.smul_const u)).const_smul (Real.sign r)
  refine ⟨_, hG.congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_), ?_⟩
  · simp only [torusGrad, div_eq_mul_inv, mul_inv, Pi.smul_apply, Pi.add_apply, Pi.mul_apply]
    congr 1
  · rw [trace3_smul, trace3_add, trace3_add, trace3_smul, trace3_sub, trace3_id,
      trace3_smulRight, trace3_smulRight, trace3_smulRight]
    have hT : torusGrad c u R 1 x = ((v - R) / (q * v)) • P + (z / q) • u := by
      simp [torusGrad, Real.sign_one, hvdef, hqdef, hzdef, hPdef]
    simp only [Pi.mul_apply, ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply,
      _root_.smul_apply, dualOf_apply, smul_eq_mul]
    simp only [← hvdef, ← hqdef, ← hzdef, ← hPdef, hT, add_dotProduct, dotProduct_add,
      smul_dotProduct, dotProduct_smul, hu, hPP, hPu, dotProduct_comm u P, smul_eq_mul,
      ContinuousLinearMap.one_apply, (dualOf_smul), map_smul]
    have hv0 := hvpos.ne'
    have hq0 := hqpos.ne'
    field_simp
    linear_combination (Real.sign r * v) * hq2

end PAT

end DDG
