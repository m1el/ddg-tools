import DdgProofs.Developability.Basic

/-!
# Developability: closed forms on the regular cone (SPEC §10.3)

`crates/ddg-developability/SPEC.md` §10.3. The regular cone of valence `k ≥ 3` and height `h > 0`:
apex `p = (0, 0, h)`, ring `w_m = (cos 2πm/k, sin 2πm/k, 0)`, faces `(p, w_m, w_{m+1})`. With
`δ = π/k`, `φ_m = (2m + 1)π/k` and `ρ = √(h² + cos² δ)`:

* the face normals are `n_m = 2 sin δ (h cos φ_m, h sin φ_m, cos δ)`, the unit normals
  `N_m = (h cos φ_m, h sin φ_m, cos δ)/ρ` (`cone_faceNormal`, `cone_unitNormal`); all apex angles
  are equal (`cone_angle`);
* the objective of Eq. 3 is `Σ θ ⟨u, N_m⟩² = θ k (h²/(2ρ²) (u_x² + u_y²) + cos²δ/ρ² u_z²)`
  (`cone_covObjective`), so `E_λ = θ k min(sin²α/2, cos²α)` with `sin²α = h²/ρ²`,
  `cos²α = cos²δ/ρ²` (`cone_covariance_energy`);
* the vertex normal is `e_z`, the log-mapped normals are `Ñ_m = α (cos φ_m, sin φ_m, 0)` with
  `α = atan2(h, cos δ)`, and the intrinsic objective is `θ k α²/2` in every tangent direction
  (`cone_intrinsic_energy`): the intrinsic energy is `θ k α²/2` (a double eigenvalue).

Not proved here (checked numerically by the Rust tests and the oracle only): the closed forms of
§10.3 for `E_P`, `λ_max` and the intrinsic max energy.

The trigonometric sums `Σ_m cos φ_m = Σ_m sin φ_m = Σ_m cos 2φ_m = Σ_m sin 2φ_m = 0` (`k ≥ 3`)
come from the geometric sum of the `k`-th roots of unity.
-/

noncomputable section

namespace DDG

namespace Dev

open scoped DDG.Oracle
open Real Ddg.Dev

/-! ## Sums over the roots of unity -/

/-- `Σ_{m<k} cos(x + m y) = 0` and `Σ_{m<k} sin(x + m y) = 0` when `e^{iy} ≠ 1 = e^{iky}`. -/
theorem sum_cos_sin_arith (k : ℕ) (x y : ℝ) (h1 : Complex.exp (y * Complex.I) ≠ 1)
    (h2 : Complex.exp (k * y * Complex.I) = 1) :
    ∑ m ∈ Finset.range k, cos (x + m * y) = 0 ∧ ∑ m ∈ Finset.range k, sin (x + m * y) = 0 := by
  have hg : ∑ m ∈ Finset.range k, Complex.exp (y * Complex.I) ^ m = 0 := by
    rw [geom_sum_eq h1, ← Complex.exp_nat_mul, show (k : ℂ) * (y * Complex.I) = k * y * Complex.I by
      ring, h2, sub_self, zero_div]
  have hs : ∑ m ∈ Finset.range k, Complex.exp (((x + m * y : ℝ) : ℂ) * Complex.I) = 0 := by
    have e : ∀ m : ℕ, Complex.exp (((x + m * y : ℝ) : ℂ) * Complex.I) =
        Complex.exp (x * Complex.I) * Complex.exp (y * Complex.I) ^ m := by
      intro m
      rw [← Complex.exp_nat_mul, ← Complex.exp_add]
      congr 1
      push_cast
      ring
    simp only [e, ← Finset.mul_sum, hg, mul_zero]
  constructor
  · have := congrArg Complex.re hs
    rw [Complex.re_sum] at this
    simp only [Complex.exp_ofReal_mul_I_re, Complex.zero_re] at this
    exact this
  · have := congrArg Complex.im hs
    rw [Complex.im_sum] at this
    simp only [Complex.exp_ofReal_mul_I_im, Complex.zero_im] at this
    exact this

/-- `e^{2πi j/k} ≠ 1` for `0 < j < k`. -/
theorem exp_two_pi_div_ne_one {j k : ℕ} (hj : 0 < j) (hjk : j < k) :
    Complex.exp (((2 * π * j / k : ℝ) : ℂ) * Complex.I) ≠ 1 := by
  rw [Ne, Complex.exp_eq_one_iff]
  rintro ⟨n, hn⟩
  have hk : (0 : ℝ) < k := by exact_mod_cast (by omega : 0 < k)
  have hre : (2 * π * j / k : ℝ) = n * (2 * π) := by
    have := congrArg Complex.im hn
    simpa using this
  have hjk' : (j : ℝ) / k = n := by
    have hπ := pi_pos
    field_simp at hre ⊢
    nlinarith [hre]
  have h0 : (0 : ℝ) < j / k := div_pos (by exact_mod_cast hj) hk
  have h1 : (j : ℝ) / k < 1 := (div_lt_one hk).mpr (by exact_mod_cast hjk)
  rw [hjk'] at h0 h1
  have : (0 : ℤ) < n := by exact_mod_cast h0
  have : n < (1 : ℤ) := by exact_mod_cast h1
  omega

theorem exp_two_pi_mul_eq_one (j k : ℕ) (hk : 0 < k) :
    Complex.exp ((k : ℂ) * ((2 * π * j / k : ℝ) : ℂ) * Complex.I) = 1 := by
  rw [Complex.exp_eq_one_iff]
  refine ⟨j, ?_⟩
  have hk' : (k : ℂ) ≠ 0 := by exact_mod_cast hk.ne'
  push_cast
  field_simp

/-- The four trigonometric sums of the regular cone, `φ_m = (2m + 1)π/k`, `k ≥ 3`. -/
theorem cone_trig_sums {k : ℕ} (hk : 3 ≤ k) :
    ∑ m ∈ Finset.range k, cos ((2 * m + 1) * π / k) = 0 ∧
    ∑ m ∈ Finset.range k, sin ((2 * m + 1) * π / k) = 0 ∧
    ∑ m ∈ Finset.range k, cos (2 * ((2 * m + 1) * π / k)) = 0 ∧
    ∑ m ∈ Finset.range k, sin (2 * ((2 * m + 1) * π / k)) = 0 := by
  have hk0 : 0 < k := by omega
  have hkr : (k : ℝ) ≠ 0 := by exact_mod_cast hk0.ne'
  obtain ⟨c1, s1⟩ := sum_cos_sin_arith k (π / k) (2 * π * (1 : ℕ) / k)
    (exp_two_pi_div_ne_one (by norm_num) (by omega))
    (by have := exp_two_pi_mul_eq_one 1 k hk0; push_cast at this ⊢; exact this)
  obtain ⟨c2, s2⟩ := sum_cos_sin_arith k (2 * π / k) (2 * π * (2 : ℕ) / k)
    (exp_two_pi_div_ne_one (by norm_num) (by omega))
    (by have := exp_two_pi_mul_eq_one 2 k hk0; push_cast at this ⊢; exact this)
  have e1 : ∀ m : ℕ, (2 * m + 1) * π / k = π / k + m * (2 * π * ((1 : ℕ) : ℝ) / k) := by
    intro m; push_cast; field_simp; ring
  have e2 : ∀ m : ℕ, 2 * ((2 * m + 1) * π / k) = 2 * π / k + m * (2 * π * ((2 : ℕ) : ℝ) / k) := by
    intro m; push_cast; field_simp; ring
  simp only [e2]
  simp only [e1]
  exact ⟨c1, s1, c2, s2⟩


/-! ## The regular cone -/

/-- Ring vertex `w_m = (cos 2πm/k, sin 2πm/k, 0)`. -/
def coneRing (k m : ℕ) : V3 := ⟨cos (2 * π * m / k), sin (2 * π * m / k), 0⟩

/-- Apex `(0, 0, h)`. -/
def coneApex (h : ℝ) : V3 := ⟨0, 0, h⟩

/-- `φ_m = (2m + 1)π/k`, the azimuth of the normal of face `m`. -/
def conePhi (k m : ℕ) : ℝ := (2 * m + 1) * π / k

/-- `ρ = √(h² + cos²(π/k))`. -/
def coneRho (k : ℕ) (h : ℝ) : ℝ := √(h ^ 2 + cos (π / k) ^ 2)

/-- The common apex angle `atan2(2 sin(π/k) ρ, cos(2π/k) + h²)`. -/
def coneTheta (k : ℕ) (h : ℝ) : ℝ := Ddg.Dev.atan2 (2 * sin (π / k) * coneRho k h) (cos (2 * π / k) + h ^ 2)

/-- `(h cos φ_m, h sin φ_m, cos δ)`. -/
def coneDir (k m : ℕ) (h : ℝ) : V3 := ⟨h * cos (conePhi k m), h * sin (conePhi k m), cos (π / k)⟩

section cone

variable {k : ℕ} (hk : 3 ≤ k) (h : ℝ)
include hk

theorem cone_ring_angles (m : ℕ) :
    2 * π * (m : ℕ) / k = conePhi k m - π / k ∧ 2 * π * ((m + 1 : ℕ) : ℝ) / k = conePhi k m + π / k := by
  have hkr : (k : ℝ) ≠ 0 := by exact_mod_cast (by omega : k ≠ 0)
  unfold conePhi
  constructor <;> push_cast <;> field_simp <;> ring

theorem sin_cos_delta_pos : 0 < sin (π / k) ∧ 0 < cos (π / k) := by
  have hkr : (3 : ℝ) ≤ k := by exact_mod_cast hk
  have hk0 : (0 : ℝ) < k := by linarith
  have h1 : 0 < π / k := div_pos pi_pos hk0
  have h2 : π / k ≤ π / 3 := div_le_div_of_nonneg_left pi_pos.le (by norm_num) hkr
  constructor
  · exact sin_pos_of_pos_of_lt_pi h1 (by linarith [pi_pos])
  · exact cos_pos_of_mem_Ioo ⟨by linarith [pi_pos], by linarith [pi_pos]⟩

theorem coneRho_pos : 0 < coneRho k h := by
  unfold coneRho
  apply sqrt_pos.mpr
  have := (sin_cos_delta_pos hk).2
  positivity

/-- The face normals of the regular cone: `n_m = 2 sin δ (h cos φ_m, h sin φ_m, cos δ)`. -/
theorem cone_faceNormal (m : ℕ) :
    faceNormal (coneApex h) (coneRing k m) (coneRing k (m + 1)) =
      smul (2 * sin (π / k)) (coneDir k m h) := by
  obtain ⟨e0, e1⟩ := cone_ring_angles hk m
  unfold coneRing
  rw [e0, e1]
  have hp := sin_sq_add_cos_sq (conePhi k m)
  ext <;> simp only [faceNormal_def, coneApex, coneDir, cross_x, cross_y, cross_z, sub_x, sub_y,
    sub_z, smul_x, smul_y, smul_z] <;> simp only [sin_add, sin_sub, cos_add, cos_sub]
  · ring
  · ring
  · linear_combination (2 * sin (π / k) * cos (π / k)) * hp

omit hk in
theorem norm_coneDir (m : ℕ) : norm (coneDir k m h) = coneRho k h := by
  rw [norm_def, coneRho]
  congr 1
  simp only [norm2_def, coneDir]
  linear_combination h ^ 2 * sin_sq_add_cos_sq (conePhi k m)

/-- The unit normals of the regular cone: `N_m = (h cos φ_m, h sin φ_m, cos δ)/ρ`. -/
theorem cone_unitNormal (m : ℕ) :
    faceUnitNormal (coneApex h) (coneRing k m) (coneRing k (m + 1)) =
      smul (1 / coneRho k h) (coneDir k m h) := by
  have hs := (sin_cos_delta_pos hk).1
  have hr := coneRho_pos hk h
  rw [faceUnitNormal_def, cone_faceNormal hk, normalize_def,
    norm_smul_of_nonneg _ (by positivity), norm_coneDir]
  ext <;> simp only [smul_x, smul_y, smul_z] <;> field_simp

/-- All apex angles of the regular cone are equal, to `coneTheta k h`. -/
theorem cone_angle (m : ℕ) :
    cornerAngleAt (coneApex h) (coneRing k m) (coneRing k (m + 1)) = coneTheta k h := by
  have hs := (sin_cos_delta_pos hk).1
  unfold cornerAngleAt coneTheta
  have hn := cone_faceNormal hk h m
  rw [faceNormal_def] at hn
  rw [hn, norm_smul_of_nonneg _ (by positivity), norm_coneDir]
  congr 1
  obtain ⟨e0, e1⟩ := cone_ring_angles hk m
  have hc : cos (2 * π / k) = cos (conePhi k m + π / k) * cos (conePhi k m - π / k) +
      sin (conePhi k m + π / k) * sin (conePhi k m - π / k) := by
    rw [← cos_sub]; congr 1; ring
  rw [hc]
  unfold coneRing
  rw [e0, e1]
  simp only [dot_def, coneApex, sub_x, sub_y, sub_z]
  ring

omit hk in
/-- A quadratic form in `(cos φ, sin φ)` written with the double angle. -/
theorem quad_trig (A B C D E F φ : ℝ) :
    A * cos φ ^ 2 + B * sin φ ^ 2 + C * (cos φ * sin φ) + D * cos φ + E * sin φ + F =
      ((A + B) / 2 + F) + (A - B) / 2 * cos (2 * φ) + C / 2 * sin (2 * φ) + D * cos φ +
        E * sin φ := by
  rw [cos_two_mul, sin_two_mul]
  linear_combination B * sin_sq_add_cos_sq φ

/-- `Σ_m (K + P cos 2φ_m + Q sin 2φ_m + R cos φ_m + S sin φ_m) = k K`. -/
theorem cone_sum_trig (K P Q R S : ℝ) :
    ∑ m ∈ Finset.range k, (K + P * cos (2 * conePhi k m) + Q * sin (2 * conePhi k m) +
      R * cos (conePhi k m) + S * sin (conePhi k m)) = k * K := by
  obtain ⟨c1, s1, c2, s2⟩ := cone_trig_sums hk
  unfold conePhi
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum, c1, s1, c2, s2, Finset.sum_const,
    Finset.card_range, nsmul_eq_mul]
  ring

/-- **Eq. 3 on the regular cone**: `Σ_m θ ⟨u, N_m⟩² = θ k (h²/(2ρ²)(u_x² + u_y²) + cos²δ/ρ² u_z²)`.

[formula: dev-regular-cone-values] -/
theorem cone_covObjective (u : V3) :
    covObjective k (fun m => cornerAngleAt (coneApex h) (coneRing k m) (coneRing k (m + 1)))
        (fun m => faceUnitNormal (coneApex h) (coneRing k m) (coneRing k (m + 1))) u =
      coneTheta k h * k * (h ^ 2 / (2 * coneRho k h ^ 2) * (u.x ^ 2 + u.y ^ 2) +
        cos (π / k) ^ 2 / coneRho k h ^ 2 * u.z ^ 2) := by
  have hr := (coneRho_pos hk h).ne'
  unfold covObjective
  rw [sumTo_eq_sum_range]
  simp only [cone_angle hk, cone_unitNormal hk, dot_smul_r]
  have e : ∀ m, coneTheta k h * (1 / coneRho k h * dot u (coneDir k m h) *
      (1 / coneRho k h * dot u (coneDir k m h))) =
      coneTheta k h / coneRho k h ^ 2 * (((h ^ 2 * u.x ^ 2 + h ^ 2 * u.y ^ 2) / 2 +
        u.z ^ 2 * cos (π / k) ^ 2) + (h ^ 2 * u.x ^ 2 - h ^ 2 * u.y ^ 2) / 2 * cos (2 * conePhi k m) +
        (2 * h ^ 2 * u.x * u.y) / 2 * sin (2 * conePhi k m) +
        (2 * h * u.x * u.z * cos (π / k)) * cos (conePhi k m) +
        (2 * h * u.y * u.z * cos (π / k)) * sin (conePhi k m)) := by
    intro m
    rw [← quad_trig]
    simp only [dot_def, coneDir]
    field_simp
    ring
  rw [Finset.sum_congr rfl (fun m _ => e m), ← Finset.mul_sum, cone_sum_trig hk]
  field_simp

/-- **`E_λ` of the regular cone** (SPEC §10.3): the minimum of Eq. 3 over unit vectors is
`θ k min(sin²α/2, cos²α)` with `sin²α = h²/ρ²`, `cos²α = cos²(π/k)/ρ²`, `θ` the apex angle
(`coneTheta`); the smaller eigenvalue is double when `sin²α/2 < cos²α`.

[formula: dev-regular-cone-values] -/
theorem cone_covariance_energy :
    IsLeast {q | ∃ u : V3, norm2 u = 1 ∧
        covObjective k (fun m => cornerAngleAt (coneApex h) (coneRing k m) (coneRing k (m + 1)))
          (fun m => faceUnitNormal (coneApex h) (coneRing k m) (coneRing k (m + 1))) u = q}
      (coneTheta k h * k * min (h ^ 2 / coneRho k h ^ 2 / 2) (cos (π / k) ^ 2 / coneRho k h ^ 2)) := by
  have hθ : 0 ≤ coneTheta k h := (atan2_mem_Icc (by
    have := (sin_cos_delta_pos hk).1; have := coneRho_pos hk h; positivity) _).1
  have hkr : (0 : ℝ) ≤ k := by positivity
  set a := h ^ 2 / coneRho k h ^ 2 / 2
  set b := cos (π / k) ^ 2 / coneRho k h ^ 2
  have hval : ∀ u : V3, covObjective k (fun m => cornerAngleAt (coneApex h) (coneRing k m)
      (coneRing k (m + 1))) (fun m => faceUnitNormal (coneApex h) (coneRing k m)
      (coneRing k (m + 1))) u = coneTheta k h * k * (a * (u.x ^ 2 + u.y ^ 2) + b * u.z ^ 2) := by
    intro u
    rw [cone_covObjective hk]
    simp only [a, b]
    ring
  constructor
  · rcases le_total a b with hab | hab
    · refine ⟨⟨1, 0, 0⟩, by simp [norm2_def], ?_⟩
      rw [hval, min_eq_left hab]; ring
    · refine ⟨⟨0, 0, 1⟩, by simp [norm2_def], ?_⟩
      rw [hval, min_eq_right hab]; ring
  · rintro q ⟨u, hu, rfl⟩
    rw [hval]
    rw [norm2_def] at hu
    apply mul_le_mul_of_nonneg_left _ (mul_nonneg hθ hkr)
    have h1 : u.x ^ 2 + u.y ^ 2 = 1 - u.z ^ 2 := by nlinarith
    have hz : 0 ≤ u.z ^ 2 := sq_nonneg _
    have hz1 : u.z ^ 2 ≤ 1 := by nlinarith [sq_nonneg u.x, sq_nonneg u.y]
    rw [h1]
    have := min_le_left a b
    have := min_le_right a b
    nlinarith

/-- The vector area of the cone star is `(0, 0, k sin δ cos δ)`. -/
theorem cone_vectorArea :
    vectorArea k (fun m => (coneApex h, coneRing k m, coneRing k (m + 1))) =
      ⟨0, 0, k * (sin (π / k) * cos (π / k))⟩ := by
  obtain ⟨c1, s1, -, -⟩ := cone_trig_sums hk
  unfold vectorArea vsumTo
  simp only [cone_faceNormal hk]
  ext <;> simp only [smul_x, smul_y, smul_z, coneDir, conePhi, sumTo_eq_sum_range,
    Oracle.one_eq, Oracle.two_eq]
  · rw [show (∑ i ∈ Finset.range k, 2 * sin (π / k) * (h * cos ((2 * i + 1) * π / k))) =
      2 * sin (π / k) * h * ∑ i ∈ Finset.range k, cos ((2 * i + 1) * π / k) by
        rw [Finset.mul_sum]; congr 1; ext; ring, c1]; ring
  · rw [show (∑ i ∈ Finset.range k, 2 * sin (π / k) * (h * sin ((2 * i + 1) * π / k))) =
      2 * sin (π / k) * h * ∑ i ∈ Finset.range k, sin ((2 * i + 1) * π / k) by
        rw [Finset.mul_sum]; congr 1; ext; ring, s1]; ring
  · simp only [Finset.sum_const, Finset.card_range, nsmul_eq_mul]; ring

/-- The vertex normal of the cone apex is `e_z`. -/
theorem cone_vertexNormal :
    vertexNormal k (fun m => (coneApex h, coneRing k m, coneRing k (m + 1))) = ⟨0, 0, 1⟩ := by
  obtain ⟨hs, hc⟩ := sin_cos_delta_pos hk
  have hk0 : (0 : ℝ) < k := by exact_mod_cast (by omega : 0 < k)
  have hpos : 0 < k * (sin (π / k) * cos (π / k)) := by positivity
  unfold vertexNormal
  rw [cone_vectorArea hk, normalize_def, norm_def]
  have : norm2 (⟨0, 0, k * (sin (π / k) * cos (π / k))⟩ : V3) =
      (k * (sin (π / k) * cos (π / k))) ^ 2 := by simp [norm2_def]; ring
  rw [this, sqrt_sq hpos.le]
  ext <;> simp only [smul_x, smul_y, smul_z] <;> field_simp <;> simp

/-- **Log map on the cone**: for `h > 0`, `log_{e_z}(N_m) = α (cos φ_m, sin φ_m, 0)` with
`α = atan2(h, cos δ)`. -/
theorem cone_logMap (hh : 0 < h) (m : ℕ) :
    logMap ⟨0, 0, 1⟩ (faceUnitNormal (coneApex h) (coneRing k m) (coneRing k (m + 1))) =
      smul (Ddg.Dev.atan2 h (cos (π / k))) ⟨cos (conePhi k m), sin (conePhi k m), 0⟩ := by
  have hr := coneRho_pos hk h
  rw [cone_unitNormal hk]
  set ρ := coneRho k h
  have hcr : norm (cross ⟨0, 0, 1⟩ (smul (1 / ρ) (coneDir k m h))) = h / ρ := by
    rw [norm_def]
    have : norm2 (cross ⟨0, 0, 1⟩ (smul (1 / ρ) (coneDir k m h))) = (h / ρ) ^ 2 := by
      simp only [norm2_def, cross_x, cross_y, cross_z, smul_x, smul_y, smul_z, coneDir]
      field_simp
      linear_combination h ^ 2 * sin_sq_add_cos_sq (conePhi k m)
    rw [this, sqrt_sq (by positivity)]
  have hdot : dot ⟨0, 0, 1⟩ (smul (1 / ρ) (coneDir k m h)) = (1 / ρ) * cos (π / k) := by
    simp [dot_def, coneDir]
  have hang : logAngle ⟨0, 0, 1⟩ (smul (1 / ρ) (coneDir k m h)) = Ddg.Dev.atan2 h (cos (π / k)) := by
    unfold logAngle
    rw [hcr, hdot, show h / ρ = (1 / ρ) * h by ring, atan2_smul (by positivity)]
  unfold logMap logG
  rw [hang, hcr, hdot]
  ext <;> simp only [smul_x, smul_y, smul_z, sub_x, sub_y, sub_z, coneDir] <;> field_simp <;> ring

/-- **The intrinsic energy of the regular cone** (SPEC §10.3): for `h > 0` the vertex normal is
`e_z` and, in every unit tangent direction `u` (`u_z = 0`), the intrinsic objective
`Σ θ ⟨u, Ñ_m⟩²` equals `θ k α²/2` with `α = atan2(h, cos(π/k))`: the 2x2 matrix `Ã` is
`(θ k α²/2) I`, a double eigenvalue, and the intrinsic energy is `θ k α²/2`.

[formula: dev-regular-cone-values] -/
theorem cone_intrinsic_energy (hh : 0 < h) (u : V3) (hz : u.z = 0) (hu : u.x ^ 2 + u.y ^ 2 = 1) :
    vertexNormal k (fun m => (coneApex h, coneRing k m, coneRing k (m + 1))) = ⟨0, 0, 1⟩ ∧
    covObjective k (fun m => cornerAngleAt (coneApex h) (coneRing k m) (coneRing k (m + 1)))
        (fun m => logMap (vertexNormal k (fun m => (coneApex h, coneRing k m, coneRing k (m + 1))))
          (faceUnitNormal (coneApex h) (coneRing k m) (coneRing k (m + 1)))) u =
      coneTheta k h * k * Ddg.Dev.atan2 h (cos (π / k)) ^ 2 / 2 := by
  refine ⟨cone_vertexNormal hk h, ?_⟩
  rw [cone_vertexNormal hk h]
  unfold covObjective
  rw [sumTo_eq_sum_range]
  simp only [cone_angle hk, cone_logMap hk h hh]
  set α := Ddg.Dev.atan2 h (cos (π / k))
  have e : ∀ m, coneTheta k h * (dot u (smul α ⟨cos (conePhi k m), sin (conePhi k m), 0⟩) *
      dot u (smul α ⟨cos (conePhi k m), sin (conePhi k m), 0⟩)) =
      coneTheta k h * α ^ 2 * ((1 / 2 + 0) + (u.x ^ 2 - u.y ^ 2) / 2 * cos (2 * conePhi k m) +
        (2 * u.x * u.y) / 2 * sin (2 * conePhi k m) + 0 * cos (conePhi k m) +
        0 * sin (conePhi k m)) := by
    intro m
    have hq := quad_trig (u.x ^ 2) (u.y ^ 2) (2 * u.x * u.y) 0 0 0 (conePhi k m)
    rw [show (u.x ^ 2 + u.y ^ 2) / 2 + 0 = 1 / 2 + 0 by rw [hu]] at hq
    rw [← hq]
    simp only [dot_def, smul_x, smul_y, smul_z, hz]
    ring
  rw [Finset.sum_congr rfl (fun m _ => e m), ← Finset.mul_sum, cone_sum_trig hk]
  ring

end cone

end Dev

end DDG
