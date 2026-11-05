import DdgProofs.PointsAsTori.Apps

/-!
# Points as Tori: baselines (Sec. 5.2, Eq. 28–29; SPEC §9)

* SSPD (Eq. 29): `φ = Σ A_i ⟨x − p_i, n_i⟩ e^{−λ(ρ_i−σ)} / Σ A_i e^{−λ(ρ_i−σ)}`, the blend of the
  tangent planes with area weights; shift invariant and with gradient
  `∇φ = Σ w_i (n_i − λ(g_i − φ) r̂_i)/W` (`sspd_gradient`).
* SHC (Eq. 28): `w = Σ A_i (λρ_i + 1)⟨x − p_i, n_i⟩/(2πρ_i³) e^{−λ(ρ_i−σ)}`,
  `φ = sgn(w)(−log|w|/λ + σ)`; the gradient of one kernel term (`shc_kernel_gradient`) and
  `∇φ = −∇w/(λ|w|)` where `w ≠ 0` (`shc_gradient`).
* Fast winding numbers: the Green's function `G(y) = 1/(4π|y|)`, its gradient and Hessian (the
  order-0 and order-1 terms of Barill et al.'s expansion), `ΔG = 0`, and the dipole identity
  `⟨n, p − x⟩/(4π|p − x|³) = n·∇G(x − p)` (`fwn_expansion`).
-/

noncomputable section

namespace DDG

namespace PAT

open Matrix Real Finset
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff vnorm_smul
  dot_self_nonneg)

variable {ι : Type*} [Fintype ι]

/-! ## SSPD (Eq. 29) -/

/-- **SSPD** (Eq. 29) is the blend (Eq. 25) of the tangent-plane SDFs `⟨x − p_i, n_i⟩` with the
area weights; in particular it does not depend on the shift `σ`.

[formula: pat-sspd] -/
theorem sspd_eq (A : ι → ℝ) (p n : ι → V3) (lam σ : ℝ) (x : V3) :
    (∑ i, A i * Real.exp (-lam * (vnorm (x - p i) - σ)) * ((x - p i) ⬝ᵥ n i)) /
        ∑ i, A i * Real.exp (-lam * (vnorm (x - p i) - σ)) =
      blend A p (fun i y => (y - p i) ⬝ᵥ n i) lam 0 x :=
  blend_shift A p (fun i y => (y - p i) ⬝ᵥ n i) lam σ x

theorem hasFDerivAt_plane (p n x : V3) :
    HasFDerivAt (fun y => (y - p) ⬝ᵥ n) (dualOf n) x := by
  refine ((hasFDerivAt_dot_const n x).sub_const (n ⬝ᵥ p)).congr_of_eventuallyEq
    (Filter.Eventually.of_forall fun y => ?_)
  simp [sub_dotProduct, dotProduct_comm n]

/-- **The SSPD gradient** with frozen `λ, σ` and positive areas, away from the data points:
`∇φ = Σ w_i (n_i − λ(g_i − φ) r̂_i)/W`, `w_i = A_i e^{−λ(ρ_i−σ)}`, `g_i = ⟨x − p_i, n_i⟩`.

[formula: pat-sspd-gradient] -/
theorem sspd_gradient [Nonempty ι] {A : ι → ℝ} (hA : ∀ i, 0 < A i) (p n : ι → V3) (lam σ : ℝ)
    {x : V3} (hx : ∀ i, x ≠ p i) :
    let w := fun i => A i * Real.exp (-lam * (vnorm (x - p i) - σ))
    let φ := (∑ i, w i * ((x - p i) ⬝ᵥ n i)) / ∑ i, w i
    HasFDerivAt (fun y => (∑ i, A i * Real.exp (-lam * (vnorm (y - p i) - σ)) * ((y - p i) ⬝ᵥ n i)) /
        ∑ i, A i * Real.exp (-lam * (vnorm (y - p i) - σ)))
      (dualOf ((∑ i, w i)⁻¹ • ∑ i, w i • (n i - (lam * ((x - p i) ⬝ᵥ n i - φ)) • rhat (p i) x))) x := by
  intro w φ
  have hexp : ∀ i y, A i * Real.exp (-lam * (vnorm (y - p i) - σ)) =
      Real.exp (Real.log (A i) + -lam * (vnorm (y - p i) - σ)) := by
    intro i y; rw [Real.exp_add, Real.exp_log (hA i)]
  have he : ∀ i, HasFDerivAt (fun y => Real.log (A i) + -lam * (vnorm (y - p i) - σ))
      (dualOf ((-lam) • rhat (p i) x)) x := by
    intro i
    refine ((((hasFDerivAt_dist (hx i)).sub_const σ).const_mul (-lam)).const_add _).congr_fderiv ?_
    ext v; simp [smul_dotProduct]
  have h := blend_hasFDerivAt he (fun i => hasFDerivAt_plane (p i) (n i) x)
  refine (h.congr_fderiv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
  · congr 1
    simp only [← hexp]
    ext k
    simp only [Pi.smul_apply, Pi.add_apply, Pi.sub_apply, Finset.sum_apply, smul_eq_mul,
      Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [φ, w]
    ring
  · simp only [hexp]

/-! ## SHC (Eq. 28) -/

/-- One SHC kernel term (without the area and the constant shift factor `e^{λσ}`):
`k(x) = (λρ + 1)⟨x − p, n⟩ e^{−λρ}/(2πρ³)`, `ρ = |x − p|`. -/
def shcKernel (lam : ℝ) (p n x : V3) : ℝ :=
  (lam * vnorm (x - p) + 1) * ((x - p) ⬝ᵥ n) * Real.exp (-lam * vnorm (x - p)) /
    (2 * π * vnorm (x - p) ^ 3)

/-- **Gradient of an SHC kernel term** (`ref: signed_log_sum_exp_gradient`):
`∇k = e^{−λρ}/(2π) [(λρ + 1) n/ρ³ − (λ²ρ² + 3λρ + 3)⟨y, n⟩ y/ρ⁵]`, `y = x − p ≠ 0`.

[formula: pat-shc-gradient] -/
theorem shc_kernel_gradient (lam : ℝ) (p n : V3) {x : V3} (hx : x ≠ p) :
    HasFDerivAt (shcKernel lam p n)
      (dualOf ((Real.exp (-lam * vnorm (x - p)) / (2 * π)) •
        (((lam * vnorm (x - p) + 1) / vnorm (x - p) ^ 3) • n -
          ((lam ^ 2 * vnorm (x - p) ^ 2 + 3 * lam * vnorm (x - p) + 3) * ((x - p) ⬝ᵥ n) /
            vnorm (x - p) ^ 5) • (x - p)))) x := by
  have hρ : 0 < vnorm (x - p) := vnorm_pos_iff.mpr (sub_ne_zero.mpr hx)
  have hd := hasFDerivAt_dist hx
  have h1 := (hd.const_mul lam).add_const 1
  have h2 := hasFDerivAt_plane p n x
  have h3 := (hd.const_mul (-lam)).exp
  have h4 := (hd.pow 3).const_mul (2 * π)
  have hne : 2 * π * vnorm (x - p) ^ 3 ≠ 0 := by positivity
  have h := ((h1.mul h2).mul h3).mul (hasFDerivAt_inv_comp h4 hne)
  refine (h.congr_fderiv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
  · ext v
    simp only [rhat, _root_.add_apply, _root_.smul_apply, dualOf_apply, smul_eq_mul,
      smul_dotProduct, sub_dotProduct, nsmul_eq_mul, Pi.mul_apply]
    have hp : Real.pi ≠ 0 := Real.pi_ne_zero
    field_simp
    ring
  · simp [shcKernel, div_eq_mul_inv]

/-- **The SHC gradient** (Eq. 28) where `w ≠ 0` and `λ ≠ 0`: `φ = sgn(w)(−log|w|/λ + σ)` has
`∇φ = −∇w/(λ|w|)`.

[formula: pat-shc] -/
theorem shc_gradient {w : V3 → ℝ} {Gw : V3} {lam σ : ℝ} {x : V3}
    (hw : HasFDerivAt w (dualOf Gw) x) (h0 : w x ≠ 0) (hl : lam ≠ 0) :
    HasFDerivAt (fun y => Real.sign (w y) * (-Real.log |w y| / lam + σ))
      (dualOf ((-(lam * |w x|)⁻¹) • Gw)) x := by
  have hlog := hw.log h0
  rcases lt_or_gt_of_ne h0 with hneg | hpos
  · have ev : ∀ᶠ y in nhds x, w y < 0 := hw.continuousAt.eventually (isOpen_Iio.mem_nhds hneg)
    have hf := ((hlog.neg.const_mul lam⁻¹).add_const σ).const_mul (-1 : ℝ)
    refine (hf.congr_fderiv ?_).congr_of_eventuallyEq (ev.mono fun y hy => ?_)
    · ext v; simp [abs_of_neg hneg, smul_dotProduct]; field_simp
    · simp [Real.sign_of_neg hy, Real.log_abs, div_eq_inv_mul]
  · have ev : ∀ᶠ y in nhds x, 0 < w y := hw.continuousAt.eventually (isOpen_Ioi.mem_nhds hpos)
    have hf := ((hlog.neg.const_mul lam⁻¹).add_const σ).const_mul (1 : ℝ)
    refine (hf.congr_fderiv ?_).congr_of_eventuallyEq (ev.mono fun y hy => ?_)
    · ext v; simp [abs_of_pos hpos, smul_dotProduct]; field_simp
    · simp [Real.sign_of_pos hy, Real.log_abs, div_eq_inv_mul]

/-- SHC is shift-consistent: with `w = e^{λσ} w₀` (the shift factored out of every kernel term),
`−log|w|/λ + σ = −log|w₀|/λ`, so `φ` does not depend on `σ`.

[formula: pat-shc] -/
theorem shc_shift {w0 lam σ : ℝ} (h0 : w0 ≠ 0) (hl : lam ≠ 0) :
    Real.sign (Real.exp (lam * σ) * w0) * (-Real.log |Real.exp (lam * σ) * w0| / lam + σ) =
      Real.sign w0 * (-Real.log |w0| / lam) := by
  rw [abs_mul, abs_of_pos (Real.exp_pos _), Real.log_mul (Real.exp_pos _).ne'
    (abs_pos.mpr h0).ne', Real.log_exp]
  have hs : Real.sign (Real.exp (lam * σ) * w0) = Real.sign w0 := by
    rw [mul_comm]; exact sign_mul_pos' (Real.exp_pos _)
  rw [hs]; field_simp; ring

/-! ## Fast winding numbers -/

/-- The Green's function of the Laplacian in `ℝ³`, `G(y) = 1/(4π|y|)`. -/
def green (y : V3) : ℝ := (4 * π * vnorm y)⁻¹

/-- `∇G(y) = −y/(4π|y|³)`. -/
def greenGrad (y : V3) : V3 := (-(4 * π * vnorm y ^ 3)⁻¹) • y

/-- **The Taylor terms of the dipole kernel** (Barill et al. 2018, Sec. 4; SPEC §9): for `y ≠ 0`,
`∇G(y) = −y/(4π|y|³)`; the gradient field has derivative
`∇²G(y) h = (3(y·h) y/|y|² − h)/(4π|y|³)` (the order-1 coefficient), whose trace is `0`
(`G` is harmonic); and the GWN term is the dipole `⟨n, p − x⟩/(4π|p − x|³) = n·∇G(x − p)`.

[formula: pat-fwn-expansion] -/
theorem fwn_expansion {y : V3} (hy : y ≠ 0) :
    HasFDerivAt green (dualOf (greenGrad y)) y ∧
    HasDiv greenGrad 0 y ∧
    HasFDerivAt greenGrad
      ((4 * π * vnorm y ^ 3)⁻¹ • ((dualOf ((3 / vnorm y ^ 2) • y)).smulRight y -
        ContinuousLinearMap.id ℝ V3)) y ∧
    ∀ n p x : V3, n ⬝ᵥ (p - x) / (4 * π * vnorm (p - x) ^ 3) = n ⬝ᵥ greenGrad (x - p) := by
  have hρ : 0 < vnorm y := vnorm_pos_iff.mpr hy
  have hd : HasFDerivAt (fun z : V3 => vnorm z) (dualOf ((vnorm y)⁻¹ • y)) y := by
    have := hasFDerivAt_dist (p := 0) (x := y) (by simpa using hy)
    simpa [rhat] using this
  have hp : Real.pi ≠ 0 := Real.pi_ne_zero
  have hG : HasFDerivAt green (dualOf (greenGrad y)) y := by
    have h4 := hd.const_mul (4 * π)
    have := hasFDerivAt_inv_comp h4 (by positivity)
    refine this.congr_fderiv ?_
    ext v; simp [greenGrad, smul_dotProduct]; field_simp
  have hc : HasFDerivAt (fun z : V3 => -(4 * π * vnorm z ^ 3)⁻¹)
      (dualOf ((3 / (4 * π * vnorm y ^ 5)) • y)) y := by
    have h4 := (hd.pow 3).const_mul (4 * π)
    have := (hasFDerivAt_inv_comp h4 (by positivity)).neg
    refine this.congr_fderiv ?_
    ext v; simp [smul_dotProduct]; field_simp; try ring
  have hH := hc.smul (hasFDerivAt_id y)
  have hHess : HasFDerivAt greenGrad
      ((4 * π * vnorm y ^ 3)⁻¹ • ((dualOf ((3 / vnorm y ^ 2) • y)).smulRight y -
        ContinuousLinearMap.id ℝ V3)) y := by
    refine hH.congr_fderiv ?_
    ext v i
    simp [smul_dotProduct]
    field_simp
    ring
  refine ⟨hG, ⟨_, hHess, ?_⟩, hHess, fun n p x => ?_⟩
  · rw [trace3_smul, trace3_sub, trace3_smulRight, trace3_id, dualOf_apply, smul_dotProduct,
      ← vnorm_mul_self, smul_eq_mul]
    field_simp
    ring
  · simp only [greenGrad, dotProduct_smul, smul_eq_mul, RC.vnorm_sub_comm p x]
    rw [show p - x = -(x - p) by abel, dotProduct_neg]
    ring

end PAT

end DDG
