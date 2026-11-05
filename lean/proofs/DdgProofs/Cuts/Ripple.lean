import DdgProofs.Cuts.OnCut
import DdgProofs.Cuts.Dirichlet
import DdgProofs.Cuts.Analytic

/-!
# Cut flow: the roughening rate

`crates/ddg-variational-cuts/SPEC.md` §1 ("Roughening"), §4 ("Length weight"), D18. Linearizing
the flow about a straight cut gives the growth rate `ω k − w k²` of a ripple of wavenumber `k`,
with `ω = 2 (Q_0 + Q_1)` (the oracle's `roughRate`, `roughQ`); the linearization itself is smooth
theory and is **cited** (`CITED.md`, `cut-ripple-linearization`). Machine-checked here:

* **`cut-roughening-rate`** (`Ddg.lengthWeightForRipples`, `Ddg.roughRate`, `Ddg.roughQ`):
  with `α_L = ω ℓ/(2π)` the rate `ω k − α_L k²` vanishes at `k = 2π/ℓ`, is positive for longer
  ripples (`0 < k < 2π/ℓ`) and negative for shorter ones (`rippleRate_cutoff`,
  `rippleRate_pos_of_longer`, `rippleRate_neg_of_shorter`); in general it is positive exactly
  for `k < ω/w` and at most `ω²/(4w)`, attained at `k = ω/(2w)` (`rippleRate_pos_iff`,
  `rippleRate_le_max`); for the Dirichlet energy `Q_r = w_r d_r² ≥ 0`, so `ω ≥ 0`
  (`roughQ_dirichlet`, `roughRate_dirichlet_nonneg`).
* The two facts the linearization rests on, on closed forms: the density difference
  `Q_0 (1 − kδ)² − Q_1 (1 + kδ)²` has the derivative `−ω k` at `δ = 0`
  (`hasDerivAt_ripple_density`); `e^{−ky} cos(kx)` is harmonic in the half plane and its outward
  normal derivative on the boundary is `k` times its boundary value
  (`halfPlane_mode_harmonic`, `halfPlane_mode_dtn`).
* The sphere (SPEC §1, §6 item 5): the mode `tan(ψ/2)^m` satisfies `sin ψ · h′ = m h`, hence its
  normal derivative at the equator of the unit sphere is `m` times its value
  (`sphere_mode_ode`, `sphere_mode_dtn`), and the base solution has `u″ = 1 = K` there
  (`capU_second_derivative_equator`). The rate `(m − 1)(4/R³ − α_L (m + 1)/R²)` that SPEC
  derives from these (`sphereRippleRate`, the formula is cited) vanishes for `m = 1`, is the
  planar rate with `ω = 4/R²`, `k = m/R` up to the curvature correction `4/R³ − α_L/R²`, and is
  negative for every `m ≥ 2` iff `α_L > 4/(3R)` (`sphere_equator_stable_iff`); on the sphere of
  unit area the threshold is between `4.72` and `4.73` (`sphere_unit_area_threshold`).
-/

noncomputable section

namespace DDG

namespace Cuts

open scoped DDG.Oracle
open Real

/-! ## The oracle's definitions -/

theorem roughQ_eq (G β w b0 : ℝ) : Ddg.roughQ G β w b0 = G - β - w * b0 := rfl

theorem roughRate_eq (q0 q1 : ℝ) : Ddg.roughRate q0 q1 = 2 * (q0 + q1) := by
  simp only [Ddg.roughRate, Oracle.two_eq]

theorem lengthWeightForRipples_eq (ω ℓ : ℝ) :
    Ddg.lengthWeightForRipples ω ℓ = ω * ℓ / (2 * π) := by
  simp only [Ddg.lengthWeightForRipples, Oracle.two_eq, Oracle.pi_eq]

/-- For the Dirichlet energy the quadratic part of the region density is `w_r d²`. -/
theorem roughQ_dirichlet (w d β : ℝ) :
    Ddg.roughQ (Ddg.regionDensity w (Ddg.gradDirichletD d) β) β w 0 = w * d ^ 2 := by
  rw [roughQ_eq, regionDensity_eq, gradDirichletD_eq]
  ring

/-- For a pointwise energy it is `−w_r d e`. -/
theorem roughQ_pointwise (w b0 d e β : ℝ) :
    Ddg.roughQ (Ddg.regionDensity w (Ddg.gradPointwiseD b0 d e) β) β w b0 = -(w * d * e) := by
  rw [roughQ_eq, regionDensity_eq, gradPointwiseD_eq]
  ring

/-- **For the Dirichlet energy the roughening rate is nonnegative** (`w_r ≥ 0`). -/
theorem roughRate_dirichlet_nonneg {w0 w1 : ℝ} (h0 : 0 ≤ w0) (h1 : 0 ≤ w1) (d0 d1 β0 β1 : ℝ) :
    0 ≤ Ddg.roughRate (Ddg.roughQ (Ddg.regionDensity w0 (Ddg.gradDirichletD d0) β0) β0 w0 0)
      (Ddg.roughQ (Ddg.regionDensity w1 (Ddg.gradDirichletD d1) β1) β1 w1 0) := by
  rw [roughRate_eq, roughQ_dirichlet, roughQ_dirichlet]
  positivity

/-! ## The rate of a ripple -/

/-- Without length weight every ripple grows, the shorter the faster (the flow is ill-posed). -/
theorem rippleRate_zero_weight (ω k : ℝ) : rippleRate ω 0 k = ω * k := by
  unfold rippleRate
  ring

/-- **Which ripples grow**: for `w > 0` and `k > 0` the rate is positive exactly for
`k < ω/w`, i.e. for wavelengths longer than `2πw/ω`.

[formula: cut-roughening-rate] -/
theorem rippleRate_pos_iff {ω w k : ℝ} (hw : 0 < w) (hk : 0 < k) :
    0 < rippleRate ω w k ↔ k < ω / w := by
  have e : rippleRate ω w k = k * (ω - w * k) := by
    unfold rippleRate
    ring
  rw [e, lt_div_iff₀ hw]
  constructor
  · intro h
    have := (mul_pos_iff_of_pos_left hk).mp h
    linarith
  · intro h
    exact mul_pos hk (by linarith)

/-- **The fastest ripple**: the rate is at most `ω²/(4w)`, its value at `k = ω/(2w)`. -/
theorem rippleRate_le_max {w : ℝ} (hw : 0 < w) (ω k : ℝ) :
    rippleRate ω w k ≤ ω ^ 2 / (4 * w) ∧ rippleRate ω w (ω / (2 * w)) = ω ^ 2 / (4 * w) := by
  have hw' := hw.ne'
  constructor
  · rw [le_div_iff₀ (by positivity)]
    unfold rippleRate
    nlinarith [sq_nonneg (ω - 2 * w * k)]
  · unfold rippleRate
    field_simp
    ring

/-- **The cutoff wavelength**: with `α_L = lengthWeightForRipples ω ℓ = ω ℓ/(2π)` the rate
vanishes at the wavenumber `2π/ℓ` of the wavelength `ℓ`.

[formula: cut-roughening-rate] -/
theorem rippleRate_cutoff (ω : ℝ) {ℓ : ℝ} (hℓ : ℓ ≠ 0) :
    rippleRate ω (Ddg.lengthWeightForRipples ω ℓ) (2 * π / ℓ) = 0 := by
  rw [lengthWeightForRipples_eq]
  unfold rippleRate
  have := Real.pi_ne_zero
  field_simp
  ring

/-- The rate with the weight of the wavelength `ℓ`, factored. -/
theorem rippleRate_lengthWeight (ω ℓ k : ℝ) :
    rippleRate ω (Ddg.lengthWeightForRipples ω ℓ) k = ω * k * (2 * π - ℓ * k) / (2 * π) := by
  rw [lengthWeightForRipples_eq]
  unfold rippleRate
  have := Real.pi_ne_zero
  field_simp

/-- **Ripples longer than `ℓ` grow** (`ω > 0`): for `0 < k < 2π/ℓ` the rate is positive.

[formula: cut-roughening-rate] -/
theorem rippleRate_pos_of_longer {ω ℓ k : ℝ} (hω : 0 < ω) (hℓ : 0 < ℓ) (hk0 : 0 < k)
    (hk : k < 2 * π / ℓ) : 0 < rippleRate ω (Ddg.lengthWeightForRipples ω ℓ) k := by
  rw [rippleRate_lengthWeight]
  have h1 : ℓ * k < 2 * π := by
    rw [lt_div_iff₀ hℓ] at hk
    linarith
  exact div_pos (mul_pos (mul_pos hω hk0) (by linarith)) (by positivity)

/-- **Ripples shorter than `ℓ` decay** (`ω > 0`): for `k > 2π/ℓ` the rate is negative.

[formula: cut-roughening-rate] -/
theorem rippleRate_neg_of_shorter {ω ℓ k : ℝ} (hω : 0 < ω) (hℓ : 0 < ℓ)
    (hk : 2 * π / ℓ < k) : rippleRate ω (Ddg.lengthWeightForRipples ω ℓ) k < 0 := by
  rw [rippleRate_lengthWeight]
  have hk0 : 0 < k := lt_trans (by positivity) hk
  have h1 : 2 * π < ℓ * k := by
    rw [div_lt_iff₀ hℓ] at hk
    linarith
  exact div_neg_of_neg_of_pos (mul_neg_of_pos_of_neg (mul_pos hω hk0) (by linarith))
    (by positivity)

/-! ## The linearization, on closed forms -/

/-- **The density difference of a rippled cut.** If the normal derivative of region 0 falls by
the factor `1 − kδ` where the cut is displaced by `δ` into region 1, and that of region 1 rises by
`1 + kδ`, the difference of the quadratic parts `Q_0 (1 − kδ)² − Q_1 (1 + kδ)²` has the derivative
`−ω k` at `δ = 0`, `ω = roughRate Q_0 Q_1`. With `δ_t = −φ_t = −(G_0 − G_1)` (SPEC §1) this is
the growth rate `ω k` of the ripple.

[formula: cut-roughening-rate] -/
theorem hasDerivAt_ripple_density (Q0 Q1 k : ℝ) :
    HasDerivAt (fun δ : ℝ => Q0 * (1 - k * δ) ^ 2 - Q1 * (1 + k * δ) ^ 2)
      (-(Ddg.roughRate Q0 Q1 * k)) 0 := by
  have e : (fun δ : ℝ => Q0 * (1 - k * δ) ^ 2 - Q1 * (1 + k * δ) ^ 2) =
      fun δ => (Q0 - Q1) + (-(2 * (Q0 + Q1) * k)) * δ + (k ^ 2 * (Q0 - Q1)) * δ ^ 2 := by
    funext δ
    ring
  rw [e, roughRate_eq]
  have h1 : HasDerivAt (fun δ : ℝ => (Q0 - Q1) + (-(2 * (Q0 + Q1) * k)) * δ)
      (-(2 * (Q0 + Q1) * k) * 1) 0 :=
    ((hasDerivAt_id (0 : ℝ)).const_mul (-(2 * (Q0 + Q1) * k))).const_add (Q0 - Q1)
  have h2 : HasDerivAt (fun δ : ℝ => (k ^ 2 * (Q0 - Q1)) * δ ^ 2)
      ((k ^ 2 * (Q0 - Q1)) * (2 * 0)) 0 := (hasDerivAt_pow2 0).const_mul (k ^ 2 * (Q0 - Q1))
  exact (h1.add h2).congr_deriv (by ring)

private theorem hasDerivAt_exp_mul_cos_x (k y s : ℝ) :
    HasDerivAt (fun t : ℝ => Real.exp (-k * y) * Real.cos (k * t))
      (-k * Real.exp (-k * y) * Real.sin (k * s)) s := by
  have h := (((hasDerivAt_id s).const_mul k).cos).const_mul (Real.exp (-k * y))
  exact h.congr_deriv (by simp; ring)

private theorem hasDerivAt_exp_mul_cos_y (k x s : ℝ) :
    HasDerivAt (fun t : ℝ => Real.exp (-k * t) * Real.cos (k * x))
      (-k * Real.cos (k * x) * Real.exp (-k * s)) s := by
  have h := (((hasDerivAt_id s).const_mul (-k)).exp).mul_const (Real.cos (k * x))
  exact h.congr_deriv (by simp; ring)

/-- **The mode of a half plane is harmonic**: `∂²/∂x² + ∂²/∂y²` of `e^{−ky} cos(kx)` vanishes. -/
theorem halfPlane_mode_harmonic (k x y : ℝ) :
    deriv (fun s => deriv (fun t => Real.exp (-k * y) * Real.cos (k * t)) s) x +
      deriv (fun s => deriv (fun t => Real.exp (-k * t) * Real.cos (k * x)) s) y = 0 := by
  have e1 : (fun s => deriv (fun t => Real.exp (-k * y) * Real.cos (k * t)) s) =
      fun s => -k * Real.exp (-k * y) * Real.sin (k * s) :=
    funext fun s => (hasDerivAt_exp_mul_cos_x k y s).deriv
  have e2 : (fun s => deriv (fun t => Real.exp (-k * t) * Real.cos (k * x)) s) =
      fun s => -k * Real.cos (k * x) * Real.exp (-k * s) :=
    funext fun s => (hasDerivAt_exp_mul_cos_y k x s).deriv
  have hx2 : HasDerivAt (fun s : ℝ => -k * Real.exp (-k * y) * Real.sin (k * s))
      (-k * Real.exp (-k * y) * (Real.cos (k * x) * k)) x := by
    have h := (((hasDerivAt_id x).const_mul k).sin).const_mul (-k * Real.exp (-k * y))
    exact h.congr_deriv (by simp)
  have hy2 : HasDerivAt (fun s : ℝ => -k * Real.cos (k * x) * Real.exp (-k * s))
      (-k * Real.cos (k * x) * (Real.exp (-k * y) * -k)) y := by
    have h := (((hasDerivAt_id y).const_mul (-k)).exp).const_mul (-k * Real.cos (k * x))
    exact h.congr_deriv (by simp)
  rw [e1, e2, hx2.deriv, hy2.deriv]
  ring

/-- **Its Dirichlet-to-Neumann map is multiplication by `k`**: on the boundary `y = 0` of the
half plane `y ≥ 0` the outward normal derivative `−∂/∂y` of `e^{−ky} cos(kx)` is `k` times its
boundary value. (This is the factor `1 − kδ` of SPEC §1.) -/
theorem halfPlane_mode_dtn (k x : ℝ) :
    -deriv (fun t => Real.exp (-k * t) * Real.cos (k * x)) 0 =
      k * (Real.exp (-k * 0) * Real.cos (k * x)) := by
  rw [(hasDerivAt_exp_mul_cos_y k x 0).deriv]
  ring

/-! ## The sphere -/

theorem hasDerivAt_tan_half {ψ : ℝ} (h : Real.cos (ψ / 2) ≠ 0) :
    HasDerivAt (fun x : ℝ => Real.tan (x / 2)) (1 / (2 * Real.cos (ψ / 2) ^ 2)) ψ := by
  have h1 : HasDerivAt (fun x : ℝ => x / 2) (1 / 2) ψ := (hasDerivAt_id ψ).div_const 2
  have h2 := (Real.hasDerivAt_tan h).comp ψ h1
  exact h2.congr_deriv (by ring)

/-- **The base solution at the equator**: `u′ = tan(ψ/2)` (`hasDerivAt_capU`) has the derivative
`1` at `ψ = π/2`: `u″ = K = 1` on the unit sphere. -/
theorem capU_second_derivative_equator : deriv (fun x : ℝ => Real.tan (x / 2)) (π / 2) = 1 := by
  have e : π / 2 / 2 = π / 4 := by ring
  have hc : Real.cos (π / 2 / 2) ≠ 0 := by
    rw [e, Real.cos_pi_div_four]
    positivity
  rw [(hasDerivAt_tan_half hc).deriv, e, Real.cos_pi_div_four, div_pow,
    Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  norm_num

/-- **The modes of a hemisphere**: `h(ψ) = tan(ψ/2)^m` satisfies `sin ψ · h′ = m h` (so
`h(ψ) cos(mθ)` is harmonic on the sphere: `sin ψ (sin ψ h′)′ = m² h`). -/
theorem sphere_mode_ode (m : ℕ) {ψ : ℝ} (h : Real.cos (ψ / 2) ≠ 0) :
    Real.sin ψ * deriv (fun x : ℝ => Real.tan (x / 2) ^ m) ψ = m * Real.tan (ψ / 2) ^ m := by
  have hd := (hasDerivAt_tan_half h).fun_pow m
  rw [hd.deriv]
  have h1 : Real.sin ψ = 2 * Real.sin (ψ / 2) * Real.cos (ψ / 2) := by
    rw [← Real.sin_two_mul]
    ring_nf
  rcases m with _ | n
  · simp
  · rw [h1, Real.tan_eq_sin_div_cos]
    simp only [add_tsub_cancel_right, pow_succ]
    field_simp

/-- **The Dirichlet-to-Neumann map of a hemisphere of the unit sphere**: at the equator
`h′(π/2) = m · h(π/2) = m` for the mode `m` (the factor `m/R` of SPEC §1). -/
theorem sphere_mode_dtn (m : ℕ) : deriv (fun x : ℝ => Real.tan (x / 2) ^ m) (π / 2) = m := by
  have e : π / 2 / 2 = π / 4 := by ring
  have hc : Real.cos (π / 2 / 2) ≠ 0 := by
    rw [e, Real.cos_pi_div_four]
    positivity
  have h := sphere_mode_ode m hc
  rw [Real.sin_pi_div_two, one_mul, e, Real.tan_pi_div_four, one_pow, mul_one] at h
  exact h

/-- The growth rate of the ripple `ε cos(mθ)` of the equator of a sphere of radius `R` under
the Dirichlet energy with length weight `α` that SPEC §1 derives:
`(m − 1) (4/R³ − α (m + 1)/R²)`. (The formula is cited, `cut-ripple-linearization`.) -/
def sphereRippleRate (R α m : ℝ) : ℝ := (m - 1) * (4 / R ^ 3 - α * (m + 1) / R ^ 2)

/-- The mode `m = 1` (a tilt of the equator) is neutral. -/
theorem sphereRippleRate_one (R α : ℝ) : sphereRippleRate R α 1 = 0 := by
  unfold sphereRippleRate
  ring

/-- It is the planar rate `ω k − α k²` with `ω = 4/R²` (`Q_r = 1/R²` on both sides) and
`k = m/R`, minus the curvature correction `4/R³ − α/R²`. -/
theorem sphereRippleRate_eq_planar {R : ℝ} (hR : R ≠ 0) (α m : ℝ) :
    sphereRippleRate R α m = rippleRate (4 / R ^ 2) α (m / R) - (4 / R ^ 3 - α / R ^ 2) := by
  unfold sphereRippleRate rippleRate
  field_simp
  ring

/-- The mode `m > 1` decays iff `α > 4/(R (m + 1))`. -/
theorem sphereRippleRate_neg_iff {R : ℝ} (hR : 0 < R) {m : ℝ} (hm : 1 < m) (α : ℝ) :
    sphereRippleRate R α m < 0 ↔ 4 / (R * (m + 1)) < α := by
  have key : sphereRippleRate R α m = (m - 1) * (4 - α * (m + 1) * R) / R ^ 3 := by
    unfold sphereRippleRate
    have := hR.ne'
    field_simp
  have hm1 : 0 < m - 1 := by linarith
  have hR3 : 0 < R ^ 3 := by positivity
  have hpos : 0 < R * (m + 1) := by positivity
  rw [div_lt_iff₀ hpos, key]
  constructor
  · intro h
    by_contra hcon
    push Not at hcon
    have h1 : 0 ≤ 4 - α * (m + 1) * R := by nlinarith
    have : 0 ≤ (m - 1) * (4 - α * (m + 1) * R) / R ^ 3 :=
      div_nonneg (mul_nonneg hm1.le h1) hR3.le
    linarith
  · intro h
    have h1 : 4 - α * (m + 1) * R < 0 := by nlinarith
    exact div_neg_of_neg_of_pos (mul_neg_of_pos_of_neg hm1 h1) hR3

/-- **The equator is stable iff `α > 4/(3R)`**: every mode `m ≥ 2` decays exactly then; the
mode `m = 2` (the seam of a tennis ball) is the last to become stable.

[formula: cut-roughening-rate] -/
theorem sphere_equator_stable_iff {R : ℝ} (hR : 0 < R) (α : ℝ) :
    (∀ m : ℕ, 2 ≤ m → sphereRippleRate R α m < 0) ↔ 4 / (3 * R) < α := by
  constructor
  · intro h
    have h2 := h 2 le_rfl
    rw [sphereRippleRate_neg_iff hR (by norm_num)] at h2
    have e : R * (((2 : ℕ) : ℝ) + 1) = 3 * R := by
      push_cast
      ring
    rwa [e] at h2
  · intro h m hm
    have hm2 : (2 : ℝ) ≤ m := by exact_mod_cast hm
    rw [sphereRippleRate_neg_iff hR (by linarith)]
    refine lt_of_le_of_lt ?_ h
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith

/-- **On the sphere of unit area** (`R = 1/√(4π)`) the threshold `4/(3R)` lies between `4.72`
and `4.73` (SPEC §1: "`= 4.73`"). -/
theorem sphere_unit_area_threshold :
    4.72 < 4 / (3 * (1 / √(4 * π))) ∧ 4 / (3 * (1 / √(4 * π))) < 4.73 := by
  have hs : 0 < √(4 * π) := Real.sqrt_pos.mpr (by positivity)
  have e : 4 / (3 * (1 / √(4 * π))) = 4 / 3 * √(4 * π) := by
    field_simp
  have h1 : (3.54 : ℝ) < √(4 * π) := by
    rw [Real.lt_sqrt (by norm_num)]
    nlinarith [Real.pi_gt_d2]
  have h2 : √(4 * π) < 3.5475 := by
    rw [Real.sqrt_lt' (by norm_num)]
    nlinarith [Real.pi_lt_d4]
  rw [e]
  constructor <;> nlinarith

end Cuts

end DDG
