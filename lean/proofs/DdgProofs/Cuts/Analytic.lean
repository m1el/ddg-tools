import DdgProofs.Cuts.Step

/-!
# Cut flow: closed-form solutions

The smooth theory of *Variational Surface Cutting* (the Yamabe equation, the shape derivative,
Céa's method, the level-set evolution) is cited, not formalized (`CITED.md`, `cut-yamabe-equation`,
`cut-shape-derivative`, `cut-cea-method`, `cut-level-set-evolution`). This file machine-checks
the closed-form solutions on which the formulas of `crates/ddg-variational-cuts/SPEC.md` §1 can be
verified by one-variable calculus, and which the Rust tests compare against (SPEC §6, targets 2
and 4). `Δ` is the positive semidefinite Laplacian, so `Δu = −K` reads `u″ + u′/r = K` for a
radial function on the plane and `(sin ψ · u′)′ = K sin ψ` on the unit sphere.

* **`cut-disk-radial-solution`** — flat disk of radius `ρ`, constant source `K`:
  `u(r) = K (r² − ρ²)/4` (`diskU`) solves `u″ + u′/r = K`, `u(ρ) = 0` (`diskU_ode`,
  `diskU_laplacian`, `diskU_boundary`); `∫₀^ρ u′² 2πr dr = π K² ρ⁴/8` (`disk_dirichletEnergy`)
  and its derivative in `ρ` is `u′(ρ)² · 2πρ` (`hasDerivAt_disk_dirichletEnergy`: the Dirichlet
  density `G = (∂u/∂n)²` integrated over the boundary circle); with weights,
  `d/dρ [w J + β π ρ²] = (w u′(ρ)² + β) · 2πρ` (`hasDerivAt_disk_regionEnergy`).
  Hencky: `∫₀^ρ u² 2πr dr = π K² ρ⁶/48` (`disk_henckyEnergy`); the adjoint
  `p(r) = K r⁴/32 − K ρ² r²/8 + 3 K ρ⁴/32` (`diskP`) solves `p″ + p′/r = 2u = b′(u)`, `p(ρ) = 0`,
  `p′(ρ) = −K ρ³/8`; and the derivative of the energy in `ρ` is
  `(b(u(ρ)) − u′(ρ) p′(ρ)) · 2πρ` (`hasDerivAt_disk_henckyEnergy`: the pointwise density).
  Scaled Dirichlet with a constant weight `c`: `hasDerivAt_disk_scaledEnergy`.
* **`cut-sphere-cap-solution`** — cap `ψ < ψ₀` of the unit sphere (`K = 1`):
  `u(ψ) = 2 ln(cos(ψ₀/2)/cos(ψ/2))` (`capU`) has `u′ = tan(ψ/2)`, `(sin ψ · u′)′ = sin ψ`,
  `u(ψ₀) = 0`; `J(ψ₀) = ∫₀^{ψ₀} u′² 2π sin ψ dψ = −8π ln cos(ψ₀/2) − 4π sin²(ψ₀/2)`
  (`capEnergy_eq`), `J′(ψ₀) = tan²(ψ₀/2) · 2π sin ψ₀` (`hasDerivAt_capEnergy_integral`). For the
  two caps of a latitude circle, `d/dψ₀ [J(ψ₀) + J(π − ψ₀)] = (G_0 − G_1) · 2π sin ψ₀` with
  `G_0 − G_1 = cutSpeed` of the two densities (`hasDerivAt_twoCaps`): the sign convention of
  SPEC §1 on a closed-form example; the equator is stationary (`twoCaps_equator`).
* **`cut-curve-shortening-circle`** — `ρ(t) = √(ρ₀² − 2αt)` solves `ρ′ = −α/ρ`
  (`hasDerivAt_shrinkRadius`); the level-set function `φ = r − ρ(t)` satisfies
  `φ_t = α (φ_rr + φ_r/r)` on the cut (`circle_levelSet_flow`), which is (★) of SPEC §1 with
  `G_0 = G_1`; the length energy `α · 2πρ` decreases at the rate `−(α/ρ)² · 2πρ`
  (`hasDerivAt_circle_lengthEnergy`).
-/

noncomputable section

namespace DDG

namespace Cuts

open scoped DDG.Oracle
open Real

/-! ## Powers -/

theorem hasDerivAt_pow2 (r : ℝ) : HasDerivAt (fun x : ℝ => x ^ 2) (2 * r) r := by
  simpa using hasDerivAt_pow 2 r

theorem hasDerivAt_pow4 (r : ℝ) : HasDerivAt (fun x : ℝ => x ^ 4) (4 * r ^ 3) r := by
  simpa using hasDerivAt_pow 4 r

theorem hasDerivAt_pow6 (r : ℝ) : HasDerivAt (fun x : ℝ => x ^ 6) (6 * r ^ 5) r := by
  simpa using hasDerivAt_pow 6 r

/-! ## Flat disk, constant curvature source -/

/-- `u(r) = K (r² − ρ²)/4`: the solution of `Δu = −K` on the disk of radius `ρ` with `u = 0` on
the boundary circle (SPEC §1). -/
def diskU (K ρ r : ℝ) : ℝ := K * (r ^ 2 - ρ ^ 2) / 4

theorem hasDerivAt_diskU (K ρ r : ℝ) : HasDerivAt (diskU K ρ) (K * r / 2) r :=
  ((((hasDerivAt_pow2 r).sub_const (ρ ^ 2)).const_mul K).div_const 4).congr_deriv (by ring)

theorem deriv_diskU (K ρ : ℝ) : deriv (diskU K ρ) = fun r => K * r / 2 :=
  funext fun r => (hasDerivAt_diskU K ρ r).deriv

theorem deriv_deriv_diskU (K ρ r : ℝ) : deriv (deriv (diskU K ρ)) r = K / 2 := by
  rw [deriv_diskU]
  exact (((hasDerivAt_id r).const_mul K).div_const 2).congr_deriv (by ring) |>.deriv

/-- **The radial equation**: `u″ + u′/r = K` for `r ≠ 0`, i.e. `Δu = −K` with the positive
semidefinite Laplacian.

[formula: cut-disk-radial-solution] -/
theorem diskU_ode (K ρ : ℝ) {r : ℝ} (hr : r ≠ 0) :
    deriv (deriv (diskU K ρ)) r + deriv (diskU K ρ) r / r = K := by
  rw [deriv_deriv_diskU, deriv_diskU]
  field_simp
  ring

/-- **Boundary condition**: `u(ρ) = 0`.

[formula: cut-disk-radial-solution] -/
theorem diskU_boundary (K ρ : ℝ) : diskU K ρ ρ = 0 := by
  unfold diskU
  ring

/-- The same function in Cartesian coordinates. -/
theorem diskU_cartesian (K ρ x y : ℝ) :
    diskU K ρ √(x ^ 2 + y ^ 2) = K * (x ^ 2 + y ^ 2 - ρ ^ 2) / 4 := by
  unfold diskU
  rw [Real.sq_sqrt (by positivity)]

/-- **The Laplacian in the plane**: `∂²u/∂x² + ∂²u/∂y² = K` for
`u(x, y) = K (x² + y² − ρ²)/4`, everywhere (also at the centre).

[formula: cut-disk-radial-solution] -/
theorem diskU_laplacian (K ρ x y : ℝ) :
    deriv (fun s => deriv (fun t => K * (t ^ 2 + y ^ 2 - ρ ^ 2) / 4) s) x +
      deriv (fun s => deriv (fun t => K * (x ^ 2 + t ^ 2 - ρ ^ 2) / 4) s) y = K := by
  have h1 : ∀ c s : ℝ, HasDerivAt (fun t : ℝ => K * (t ^ 2 + c - ρ ^ 2) / 4) (K * s / 2) s := by
    intro c s
    exact (((((hasDerivAt_pow2 s).add_const c).sub_const (ρ ^ 2)).const_mul K).div_const
      4).congr_deriv (by ring)
  have h2 : ∀ c s : ℝ, HasDerivAt (fun t : ℝ => K * (c + t ^ 2 - ρ ^ 2) / 4) (K * s / 2) s := by
    intro c s
    exact (((((hasDerivAt_pow2 s).const_add c).sub_const (ρ ^ 2)).const_mul K).div_const
      4).congr_deriv (by ring)
  have h3 : ∀ s : ℝ, HasDerivAt (fun t : ℝ => K * t / 2) (K / 2) s := fun s =>
    (((hasDerivAt_id s).const_mul K).div_const 2).congr_deriv (by ring)
  have e1 : (fun s => deriv (fun t => K * (t ^ 2 + y ^ 2 - ρ ^ 2) / 4) s) = fun s => K * s / 2 :=
    funext fun s => (h1 (y ^ 2) s).deriv
  have e2 : (fun s => deriv (fun t => K * (x ^ 2 + t ^ 2 - ρ ^ 2) / 4) s) = fun s => K * s / 2 :=
    funext fun s => (h2 (x ^ 2) s).deriv
  rw [e1, e2, (h3 x).deriv, (h3 y).deriv]
  ring

/-- **Dirichlet energy of the disk**: `∫₀^ρ u′(r)² · 2πr dr = π K² ρ⁴ / 8`.

[formula: cut-disk-radial-solution] -/
theorem disk_dirichletEnergy (K ρ : ℝ) :
    ∫ r in (0 : ℝ)..ρ, deriv (diskU K ρ) r ^ 2 * (2 * π * r) = π * K ^ 2 * ρ ^ 4 / 8 := by
  have hF : ∀ r ∈ Set.uIcc (0 : ℝ) ρ,
      HasDerivAt (fun r : ℝ => π * K ^ 2 * r ^ 4 / 8) ((K * r / 2) ^ 2 * (2 * π * r)) r := by
    intro r _
    exact (((hasDerivAt_pow4 r).const_mul (π * K ^ 2)).div_const 8).congr_deriv (by ring)
  have hint : IntervalIntegrable (fun r : ℝ => (K * r / 2) ^ 2 * (2 * π * r))
      MeasureTheory.volume 0 ρ := (by fun_prop : Continuous _).intervalIntegrable _ _
  rw [deriv_diskU, intervalIntegral.integral_eq_sub_of_hasDerivAt hF hint]
  ring

/-- **Shape derivative of the Dirichlet energy of the disk**:
`d/dρ ∫₀^ρ u′² 2πr dr = u′(ρ)² · 2πρ`, the Dirichlet density `(∂u/∂n)²` integrated over the
boundary circle for the unit normal speed (paper Eq. 13, Table 1; SPEC §1).

[formula: cut-disk-radial-solution] [formula: cut-gradient-dirichlet] -/
theorem hasDerivAt_disk_dirichletEnergy (K ρ : ℝ) :
    HasDerivAt (fun ρ' : ℝ => ∫ r in (0 : ℝ)..ρ', deriv (diskU K ρ') r ^ 2 * (2 * π * r))
      (deriv (diskU K ρ) ρ ^ 2 * (2 * π * ρ)) ρ := by
  have e : (fun ρ' : ℝ => ∫ r in (0 : ℝ)..ρ', deriv (diskU K ρ') r ^ 2 * (2 * π * r)) =
      fun ρ' => π * K ^ 2 * ρ' ^ 4 / 8 := funext fun ρ' => disk_dirichletEnergy K ρ'
  rw [e, deriv_diskU]
  exact (((hasDerivAt_pow4 ρ).const_mul (π * K ^ 2)).div_const 8).congr_deriv (by ring)

/-- **With region weights**: `d/dρ [w · J(ρ) + β · πρ²] = (w u′(ρ)² + β) · 2πρ`, the oracle's
`regionDensity w G β` integrated over the boundary circle (SPEC §1, §3.6).

[formula: cut-disk-radial-solution] [formula: cut-speed] -/
theorem hasDerivAt_disk_regionEnergy (K w β ρ : ℝ) :
    HasDerivAt (fun ρ' : ℝ => w * (∫ r in (0 : ℝ)..ρ', deriv (diskU K ρ') r ^ 2 * (2 * π * r)) +
        β * (π * ρ' ^ 2))
      (Ddg.regionDensity w (deriv (diskU K ρ) ρ ^ 2) β * (2 * π * ρ)) ρ := by
  have h1 := (hasDerivAt_disk_dirichletEnergy K ρ).const_mul w
  have h2 := ((hasDerivAt_pow2 ρ).const_mul π).const_mul β
  exact (h1.add h2).congr_deriv (by rw [regionDensity_eq]; ring)

/-- **Flux balance on the disk**: `u′(ρ) · 2πρ = K · πρ²` (the boundary integral of `∂u/∂n` is
the total source).

[formula: cut-disk-radial-solution] -/
theorem disk_flux_balance (K ρ : ℝ) : deriv (diskU K ρ) ρ * (2 * π * ρ) = K * (π * ρ ^ 2) := by
  rw [deriv_diskU]
  ring

/-- **Scaled Dirichlet energy with a constant weight `c`** on the disk: the adjoint is
`p = −2c·u` and the derivative of `c ∫ u′² 2πr dr` in `ρ` is
`(−c u′(ρ)² − u′(ρ) p′(ρ)) · 2πρ` (paper Eq. 18, Table 1).

[formula: cut-disk-radial-solution] [formula: cut-gradient-scaled] -/
theorem hasDerivAt_disk_scaledEnergy (K c ρ : ℝ) :
    HasDerivAt
      (fun ρ' : ℝ => c * ∫ r in (0 : ℝ)..ρ', deriv (diskU K ρ') r ^ 2 * (2 * π * r))
      ((-(c * deriv (diskU K ρ) ρ ^ 2) -
        deriv (diskU K ρ) ρ * deriv (fun r => -2 * c * diskU K ρ r) ρ) * (2 * π * ρ)) ρ := by
  have hp : deriv (fun r => -2 * c * diskU K ρ r) ρ = -2 * c * (K * ρ / 2) :=
    ((hasDerivAt_diskU K ρ ρ).const_mul (-2 * c)).deriv
  refine ((hasDerivAt_disk_dirichletEnergy K ρ).const_mul c).congr_deriv ?_
  rw [hp, deriv_diskU]
  ring

/-- **Hencky energy of the disk**: `∫₀^ρ u(r)² · 2πr dr = π K² ρ⁶ / 48`.

[formula: cut-disk-radial-solution] -/
theorem disk_henckyEnergy (K ρ : ℝ) :
    ∫ r in (0 : ℝ)..ρ, diskU K ρ r ^ 2 * (2 * π * r) = π * K ^ 2 * ρ ^ 6 / 48 := by
  have hF : ∀ r ∈ Set.uIcc (0 : ℝ) ρ,
      HasDerivAt (fun r : ℝ => π * K ^ 2 / 48 *
        (r ^ 6 - 3 * ρ ^ 2 * r ^ 4 + 3 * ρ ^ 4 * r ^ 2 - ρ ^ 6))
        (diskU K ρ r ^ 2 * (2 * π * r)) r := by
    intro r _
    have h := (((((hasDerivAt_pow6 r).sub ((hasDerivAt_pow4 r).const_mul (3 * ρ ^ 2))).add
      ((hasDerivAt_pow2 r).const_mul (3 * ρ ^ 4))).sub_const (ρ ^ 6)).const_mul (π * K ^ 2 / 48))
    exact h.congr_deriv (by unfold diskU; ring)
  have hint : IntervalIntegrable (fun r : ℝ => diskU K ρ r ^ 2 * (2 * π * r))
      MeasureTheory.volume 0 ρ := by
    unfold diskU
    exact (by fun_prop : Continuous _).intervalIntegrable _ _
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hF hint]
  ring

/-- The adjoint of the Hencky energy on the disk:
`p(r) = K r⁴/32 − K ρ² r²/8 + 3 K ρ⁴/32` (SPEC §1: `Δp = −b′(u)`, `p = 0` on the boundary). -/
def diskP (K ρ r : ℝ) : ℝ := K * r ^ 4 / 32 - K * ρ ^ 2 * r ^ 2 / 8 + 3 * K * ρ ^ 4 / 32

theorem hasDerivAt_diskP (K ρ r : ℝ) :
    HasDerivAt (diskP K ρ) (K * r ^ 3 / 8 - K * ρ ^ 2 * r / 4) r := by
  have h := ((((hasDerivAt_pow4 r).const_mul K).div_const 32).sub
    (((hasDerivAt_pow2 r).const_mul (K * ρ ^ 2)).div_const 8)).add_const (3 * K * ρ ^ 4 / 32)
  exact h.congr_deriv (by ring)

theorem deriv_diskP (K ρ : ℝ) : deriv (diskP K ρ) = fun r => K * r ^ 3 / 8 - K * ρ ^ 2 * r / 4 :=
  funext fun r => (hasDerivAt_diskP K ρ r).deriv

/-- **The adjoint equation**: `p″ + p′/r = 2u = b′(u)` for `b(u) = u²` and `r ≠ 0`, i.e.
`Δp = −b′(u)`; `p(ρ) = 0`; `p′(ρ) = −K ρ³/8`.

[formula: cut-disk-radial-solution] -/
theorem diskP_ode (K ρ : ℝ) {r : ℝ} (hr : r ≠ 0) :
    deriv (deriv (diskP K ρ)) r + deriv (diskP K ρ) r / r = 2 * diskU K ρ r := by
  have h2 : deriv (deriv (diskP K ρ)) r = 3 * K * r ^ 2 / 8 - K * ρ ^ 2 / 4 := by
    rw [deriv_diskP]
    have h : HasDerivAt (fun x : ℝ => x ^ 3) (3 * r ^ 2) r := by simpa using hasDerivAt_pow 3 r
    exact ((((h.const_mul K).div_const 8).sub
      (((hasDerivAt_id r).const_mul (K * ρ ^ 2)).div_const 4)).congr_deriv (by ring)).deriv
  rw [h2, deriv_diskP]
  unfold diskU
  field_simp
  ring

theorem diskP_boundary (K ρ : ℝ) : diskP K ρ ρ = 0 := by
  unfold diskP
  ring

theorem deriv_diskP_boundary (K ρ : ℝ) : deriv (diskP K ρ) ρ = -(K * ρ ^ 3) / 8 := by
  rw [deriv_diskP]
  ring

/-- **Shape derivative of the Hencky energy of the disk**:
`d/dρ ∫₀^ρ u² 2πr dr = (b(u(ρ)) − u′(ρ) p′(ρ)) · 2πρ` with `b(u) = u²` — the pointwise density
`G = b(u) − (∂u/∂n)(∂p/∂n)` of paper Table 1 integrated over the boundary circle. Its value is
`π K² ρ⁵ / 8` (SPEC §1).

[formula: cut-disk-radial-solution] [formula: cut-gradient-pointwise] -/
theorem hasDerivAt_disk_henckyEnergy (K ρ : ℝ) :
    HasDerivAt (fun ρ' : ℝ => ∫ r in (0 : ℝ)..ρ', diskU K ρ' r ^ 2 * (2 * π * r))
      ((diskU K ρ ρ ^ 2 - deriv (diskU K ρ) ρ * deriv (diskP K ρ) ρ) * (2 * π * ρ)) ρ := by
  have e : (fun ρ' : ℝ => ∫ r in (0 : ℝ)..ρ', diskU K ρ' r ^ 2 * (2 * π * r)) =
      fun ρ' => π * K ^ 2 * ρ' ^ 6 / 48 := funext fun ρ' => disk_henckyEnergy K ρ'
  rw [e, diskU_boundary, deriv_diskU, deriv_diskP]
  exact (((hasDerivAt_pow6 ρ).const_mul (π * K ^ 2)).div_const 48).congr_deriv (by ring)

/-- The value of that derivative: `π K² ρ⁵ / 8`. -/
theorem disk_hencky_derivative_value (K ρ : ℝ) :
    (diskU K ρ ρ ^ 2 - deriv (diskU K ρ) ρ * deriv (diskP K ρ) ρ) * (2 * π * ρ) =
      π * K ^ 2 * ρ ^ 5 / 8 := by
  rw [diskU_boundary, deriv_diskU, deriv_diskP]
  ring

/-! ## Spherical cap -/

/-- `u(ψ) = 2 ln(cos(ψ₀/2)/cos(ψ/2))`: the solution of `Δu = −1` on the cap `ψ < ψ₀` of the unit
sphere (`ψ` the polar angle) with `u(ψ₀) = 0` (SPEC §6, target 2). -/
def capU (ψ₀ ψ : ℝ) : ℝ := 2 * Real.log (Real.cos (ψ₀ / 2) / Real.cos (ψ / 2))

/-- **Boundary condition**: `u(ψ₀) = 0`.

[formula: cut-sphere-cap-solution] -/
theorem capU_boundary (ψ₀ : ℝ) : capU ψ₀ ψ₀ = 0 := by
  unfold capU
  by_cases h : Real.cos (ψ₀ / 2) = 0
  · simp [h]
  · rw [div_self h]
    simp

theorem cos_half_pos {ψ : ℝ} (h : ψ ∈ Set.Ioo (-π) π) : 0 < Real.cos (ψ / 2) :=
  Real.cos_pos_of_mem_Ioo ⟨by linarith [h.1], by linarith [h.2]⟩

/-- **`u′ = tan(ψ/2)`.**

[formula: cut-sphere-cap-solution] -/
theorem hasDerivAt_capU {ψ₀ ψ : ℝ} (h0 : Real.cos (ψ₀ / 2) ≠ 0) (h : Real.cos (ψ / 2) ≠ 0) :
    HasDerivAt (capU ψ₀) (Real.tan (ψ / 2)) ψ := by
  have hc : HasDerivAt (fun x : ℝ => Real.cos (x / 2)) (-Real.sin (ψ / 2) * (1 / 2)) ψ :=
    ((hasDerivAt_id ψ).div_const 2).cos
  have hq : HasDerivAt (fun x : ℝ => Real.cos (ψ₀ / 2) / Real.cos (x / 2))
      ((0 * Real.cos (ψ / 2) - Real.cos (ψ₀ / 2) * (-Real.sin (ψ / 2) * (1 / 2))) /
        Real.cos (ψ / 2) ^ 2) ψ :=
    (hasDerivAt_const ψ (Real.cos (ψ₀ / 2))).div hc h
  have hl : HasDerivAt (fun x : ℝ => 2 * Real.log (Real.cos (ψ₀ / 2) / Real.cos (x / 2)))
      (2 * ((0 * Real.cos (ψ / 2) - Real.cos (ψ₀ / 2) * (-Real.sin (ψ / 2) * (1 / 2))) /
        Real.cos (ψ / 2) ^ 2 / (Real.cos (ψ₀ / 2) / Real.cos (ψ / 2)))) ψ :=
    (hq.log (div_ne_zero h0 h)).const_mul 2
  refine hl.congr_deriv ?_
  rw [Real.tan_eq_sin_div_cos]
  field_simp
  ring

theorem deriv_capU {ψ₀ ψ : ℝ} (h0 : Real.cos (ψ₀ / 2) ≠ 0) (h : Real.cos (ψ / 2) ≠ 0) :
    deriv (capU ψ₀) ψ = Real.tan (ψ / 2) :=
  (hasDerivAt_capU h0 h).deriv

/-- `sin ψ · tan(ψ/2) = 1 − cos ψ`. -/
theorem sin_mul_tan_half {ψ : ℝ} (h : Real.cos (ψ / 2) ≠ 0) :
    Real.sin ψ * Real.tan (ψ / 2) = 1 - Real.cos ψ := by
  have h1 : Real.sin ψ = 2 * Real.sin (ψ / 2) * Real.cos (ψ / 2) := by
    rw [← Real.sin_two_mul]
    ring_nf
  have h2 : Real.cos ψ = 2 * Real.cos (ψ / 2) ^ 2 - 1 := by
    rw [← Real.cos_two_mul]
    ring_nf
  have h3 := Real.sin_sq_add_cos_sq (ψ / 2)
  rw [Real.tan_eq_sin_div_cos, h1, h2]
  field_simp
  linear_combination (2 : ℝ) * h3

/-- **The equation on the sphere**: `(sin ψ · u′)′ = sin ψ` on `(−π, π)`, i.e. `Δu = −K` with
`K = 1` for a function of the polar angle.

[formula: cut-sphere-cap-solution] -/
theorem capU_ode {ψ₀ ψ : ℝ} (h0 : Real.cos (ψ₀ / 2) ≠ 0) (h : ψ ∈ Set.Ioo (-π) π) :
    HasDerivAt (fun x => Real.sin x * deriv (capU ψ₀) x) (Real.sin ψ) ψ := by
  have h1 : HasDerivAt (fun x : ℝ => 1 - Real.cos x) (Real.sin ψ) ψ :=
    ((Real.hasDerivAt_cos ψ).const_sub 1).congr_deriv (by ring)
  refine h1.congr_of_eventuallyEq ?_
  filter_upwards [Ioo_mem_nhds h.1 h.2] with x hx
  have hc := (cos_half_pos hx).ne'
  rw [deriv_capU h0 hc, sin_mul_tan_half hc]

/-- **Flux balance on the cap**: `u′(ψ₀) · 2π sin ψ₀ = 2π (1 − cos ψ₀)`, the area of the cap
(the total curvature for `K = 1`).

[formula: cut-sphere-cap-solution] -/
theorem cap_flux_balance {ψ₀ : ℝ} (h : Real.cos (ψ₀ / 2) ≠ 0) :
    Real.tan (ψ₀ / 2) * (2 * π * Real.sin ψ₀) = 2 * π * (1 - Real.cos ψ₀) := by
  rw [← sin_mul_tan_half h]
  ring

/-- The closed form of the Dirichlet energy of the cap:
`J(ψ₀) = −8π ln cos(ψ₀/2) − 4π sin²(ψ₀/2)`. -/
def capEnergy (ψ₀ : ℝ) : ℝ :=
  -(8 * π) * Real.log (Real.cos (ψ₀ / 2)) - 4 * π * Real.sin (ψ₀ / 2) ^ 2

/-- **`J′(ψ₀) = tan²(ψ₀/2) · 2π sin ψ₀`** for the closed form.

[formula: cut-sphere-cap-solution] -/
theorem hasDerivAt_capEnergy {ψ : ℝ} (h : Real.cos (ψ / 2) ≠ 0) :
    HasDerivAt capEnergy (Real.tan (ψ / 2) ^ 2 * (2 * π * Real.sin ψ)) ψ := by
  have hc : HasDerivAt (fun x : ℝ => Real.cos (x / 2)) (-Real.sin (ψ / 2) * (1 / 2)) ψ :=
    ((hasDerivAt_id ψ).div_const 2).cos
  have hs : HasDerivAt (fun x : ℝ => Real.sin (x / 2)) (Real.cos (ψ / 2) * (1 / 2)) ψ :=
    ((hasDerivAt_id ψ).div_const 2).sin
  have hG : HasDerivAt (fun x : ℝ => -(8 * π) * Real.log (Real.cos (x / 2)) -
      4 * π * (Real.sin (x / 2) * Real.sin (x / 2)))
      (-(8 * π) * (-Real.sin (ψ / 2) * (1 / 2) / Real.cos (ψ / 2)) -
        4 * π * (Real.cos (ψ / 2) * (1 / 2) * Real.sin (ψ / 2) +
          Real.sin (ψ / 2) * (Real.cos (ψ / 2) * (1 / 2)))) ψ :=
    ((hc.log h).const_mul (-(8 * π))).sub ((hs.mul hs).const_mul (4 * π))
  have efun : capEnergy = fun x : ℝ => -(8 * π) * Real.log (Real.cos (x / 2)) -
      4 * π * (Real.sin (x / 2) * Real.sin (x / 2)) := by
    funext x
    unfold capEnergy
    ring
  rw [efun]
  refine hG.congr_deriv ?_
  have h1 : Real.sin ψ = 2 * Real.sin (ψ / 2) * Real.cos (ψ / 2) := by
    rw [← Real.sin_two_mul]
    ring_nf
  have h3 := Real.sin_sq_add_cos_sq (ψ / 2)
  rw [Real.tan_eq_sin_div_cos, h1]
  generalize Real.sin (ψ / 2) = s at *
  generalize Real.cos (ψ / 2) = c at *
  have hs2 : 1 - c ^ 2 = s ^ 2 := by linarith
  calc -(8 * π) * (-s * (1 / 2) / c) - 4 * π * (c * (1 / 2) * s + s * (c * (1 / 2)))
      = 4 * π * s * (1 - c ^ 2) / c := by
        field_simp
        ring
    _ = 4 * π * s * s ^ 2 / c := by rw [hs2]
    _ = (s / c) ^ 2 * (2 * π * (2 * s * c)) := by
        field_simp
        ring

/-- **Dirichlet energy of the cap**, with the integrand written with `tan`:
`∫₀^{ψ₀} tan²(ψ/2) · 2π sin ψ dψ = −8π ln cos(ψ₀/2) − 4π sin²(ψ₀/2)` for `|ψ₀| < π`. -/
theorem capEnergy_eq_tan {ψ₀ : ℝ} (h : ψ₀ ∈ Set.Ioo (-π) π) :
    ∫ ψ in (0 : ℝ)..ψ₀, Real.tan (ψ / 2) ^ 2 * (2 * π * Real.sin ψ) = capEnergy ψ₀ := by
  have hsub : ∀ x ∈ Set.uIcc (0 : ℝ) ψ₀, x ∈ Set.Ioo (-π) π := by
    intro x hx
    rcases Set.mem_uIcc.mp hx with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact ⟨by linarith [Real.pi_pos], by linarith [h.2]⟩
    · exact ⟨by linarith [h.1], by linarith [Real.pi_pos]⟩
  have hF : ∀ x ∈ Set.uIcc (0 : ℝ) ψ₀,
      HasDerivAt capEnergy (Real.tan (x / 2) ^ 2 * (2 * π * Real.sin x)) x :=
    fun x hx => hasDerivAt_capEnergy (cos_half_pos (hsub x hx)).ne'
  have hcont : ContinuousOn (fun ψ : ℝ => Real.tan (ψ / 2) ^ 2 * (2 * π * Real.sin ψ))
      (Set.uIcc 0 ψ₀) := by
    have h1 : ContinuousOn (fun ψ : ℝ => Real.tan (ψ / 2)) (Set.uIcc 0 ψ₀) := by
      intro x hx
      have hc := (cos_half_pos (hsub x hx)).ne'
      have hdiv : ContinuousAt (fun y : ℝ => y / 2) x := by fun_prop
      have ht : ContinuousAt Real.tan (x / 2) := Real.continuousAt_tan.mpr hc
      exact (ContinuousAt.comp (f := fun y : ℝ => y / 2) (g := Real.tan) ht
        hdiv).continuousWithinAt
    exact (h1.pow 2).mul (by fun_prop)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hF hcont.intervalIntegrable]
  simp [capEnergy]

/-- **Dirichlet energy of the spherical cap** (SPEC §6, target 2):
`∫₀^{ψ₀} u′(ψ)² · 2π sin ψ dψ = −8π ln cos(ψ₀/2) − 4π sin²(ψ₀/2)` for `|ψ₀| < π`.

[formula: cut-sphere-cap-solution] -/
theorem capEnergy_eq {ψ₀ : ℝ} (h : ψ₀ ∈ Set.Ioo (-π) π) :
    ∫ ψ in (0 : ℝ)..ψ₀, deriv (capU ψ₀) ψ ^ 2 * (2 * π * Real.sin ψ) =
      -(8 * π) * Real.log (Real.cos (ψ₀ / 2)) - 4 * π * Real.sin (ψ₀ / 2) ^ 2 := by
  have h0 := (cos_half_pos h).ne'
  rw [← capEnergy, ← capEnergy_eq_tan h]
  refine intervalIntegral.integral_congr fun x hx => ?_
  have hx' : x ∈ Set.Ioo (-π) π := by
    rcases Set.mem_uIcc.mp hx with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact ⟨by linarith [Real.pi_pos], by linarith [h.2]⟩
    · exact ⟨by linarith [h.1], by linarith [Real.pi_pos]⟩
  simp only [deriv_capU h0 (cos_half_pos hx').ne']

/-- **Shape derivative of the Dirichlet energy of the cap**:
`d/dψ₀ ∫₀^{ψ₀} u′² 2π sin ψ dψ = u′(ψ₀)² · 2π sin ψ₀`, the Dirichlet density `(∂u/∂n)²`
integrated over the latitude circle of length `2π sin ψ₀` (unit sphere, unit normal speed).

[formula: cut-sphere-cap-solution] [formula: cut-gradient-dirichlet] -/
theorem hasDerivAt_capEnergy_integral {ψ₀ : ℝ} (h : ψ₀ ∈ Set.Ioo (-π) π) :
    HasDerivAt (fun ψ₁ : ℝ => ∫ ψ in (0 : ℝ)..ψ₁, deriv (capU ψ₁) ψ ^ 2 * (2 * π * Real.sin ψ))
      (deriv (capU ψ₀) ψ₀ ^ 2 * (2 * π * Real.sin ψ₀)) ψ₀ := by
  have h0 := (cos_half_pos h).ne'
  rw [deriv_capU h0 h0]
  refine (hasDerivAt_capEnergy h0).congr_of_eventuallyEq ?_
  filter_upwards [Ioo_mem_nhds h.1 h.2] with x hx
  rw [capEnergy_eq hx, capEnergy]

/-- In arclength on a sphere of radius `R` the integrand is the same: the Dirichlet energy does
not depend on the radius, and `∂u/∂n = tan(ψ₀/2)/R` (SPEC §6, target 2). -/
theorem cap_integrand_scale (ψ : ℝ) {R : ℝ} (hR : R ≠ 0) :
    (Real.tan (ψ / 2) / R) ^ 2 * (2 * π * (R * Real.sin ψ)) * R =
      Real.tan (ψ / 2) ^ 2 * (2 * π * Real.sin ψ) := by
  field_simp

/-- **The two caps of a latitude circle: the sign of the flow on a closed-form example.** Region
0 is the cap `ψ < ψ₀`, region 1 the cap `ψ > ψ₀`; moving the cut with unit speed along the
outward normal of region 0 is `dψ₀/dt = 1`. The total Dirichlet energy
`J(ψ₀) + J(π − ψ₀)` then changes at the rate `(G_0 − G_1) · 2π sin ψ₀` with the Dirichlet
densities `G_0 = tan²(ψ₀/2)`, `G_1 = tan²((π − ψ₀)/2)` — the formula `dE = ∫ (G_0 − G_1) σ ds`
of SPEC §1 with `σ = 1`, and `G_0 − G_1` is the oracle's `cutSpeed`.

[formula: cut-sphere-cap-solution] [formula: cut-speed] -/
theorem hasDerivAt_twoCaps {ψ₀ : ℝ} (h : ψ₀ ∈ Set.Ioo 0 π) :
    HasDerivAt (fun ψ₁ : ℝ => capEnergy ψ₁ + capEnergy (π - ψ₁))
      (Ddg.cutSpeed (Real.tan (ψ₀ / 2) ^ 2) (Real.tan ((π - ψ₀) / 2) ^ 2) *
        (2 * π * Real.sin ψ₀)) ψ₀ := by
  have h1 : ψ₀ ∈ Set.Ioo (-π) π := ⟨by linarith [h.1, Real.pi_pos], h.2⟩
  have h2 : π - ψ₀ ∈ Set.Ioo (-π) π :=
    ⟨by linarith [h.2, Real.pi_pos], by linarith [h.1]⟩
  have hN := hasDerivAt_capEnergy (cos_half_pos h1).ne'
  have hS : HasDerivAt (fun ψ₁ : ℝ => capEnergy (π - ψ₁))
      (Real.tan ((π - ψ₀) / 2) ^ 2 * (2 * π * Real.sin (π - ψ₀)) * -1) ψ₀ :=
    (hasDerivAt_capEnergy (cos_half_pos h2).ne').comp ψ₀ ((hasDerivAt_id ψ₀).const_sub π)
  refine (hN.add hS).congr_deriv ?_
  rw [Real.sin_pi_sub, cutSpeed_eq]
  ring

/-- **The equator is stationary**: at `ψ₀ = π/2` the speed `G_0 − G_1` of the two caps
vanishes (SPEC §6, target 5).

[formula: cut-sphere-cap-solution] -/
theorem twoCaps_equator :
    Ddg.cutSpeed (Real.tan (π / 2 / 2) ^ 2) (Real.tan ((π - π / 2) / 2) ^ 2) = 0 := by
  rw [cutSpeed_eq, show (π - π / 2) / 2 = π / 2 / 2 by ring]
  ring

/-! ## A shrinking circle -/

/-- `ρ(t) = √(ρ₀² − 2αt)`: the radius of a circle in the plane under the length term of the flow
(curve shortening with the factor `α`; SPEC §6, target 4). -/
def shrinkRadius (ρ₀ α t : ℝ) : ℝ := √(ρ₀ ^ 2 - 2 * α * t)

/-- **`ρ′ = −α/ρ`** as long as the circle has not vanished.

[formula: cut-curve-shortening-circle] -/
theorem hasDerivAt_shrinkRadius {ρ₀ α t : ℝ} (h : 0 < ρ₀ ^ 2 - 2 * α * t) :
    HasDerivAt (shrinkRadius ρ₀ α) (-α / shrinkRadius ρ₀ α t) t := by
  have hf : HasDerivAt (fun s : ℝ => ρ₀ ^ 2 - 2 * α * s) (-(2 * α * 1)) t :=
    ((hasDerivAt_id t).const_mul (2 * α)).const_sub (ρ₀ ^ 2)
  refine (hf.sqrt h.ne').congr_deriv ?_
  unfold shrinkRadius
  ring

/-- The initial radius. -/
theorem shrinkRadius_zero {ρ₀ : ℝ} (h : 0 ≤ ρ₀) (α : ℝ) : shrinkRadius ρ₀ α 0 = ρ₀ := by
  unfold shrinkRadius
  rw [mul_zero, sub_zero, Real.sqrt_sq h]

theorem shrinkRadius_pos {ρ₀ α t : ℝ} (h : 0 < ρ₀ ^ 2 - 2 * α * t) : 0 < shrinkRadius ρ₀ α t :=
  Real.sqrt_pos.mpr h

/-- The enclosed area decreases at the constant rate `2πα`: `π ρ(t)² = π ρ₀² − 2παt`. -/
theorem shrinkRadius_area {ρ₀ α t : ℝ} (h : 0 ≤ ρ₀ ^ 2 - 2 * α * t) :
    π * shrinkRadius ρ₀ α t ^ 2 = π * ρ₀ ^ 2 - 2 * π * α * t := by
  unfold shrinkRadius
  rw [Real.sq_sqrt h]
  ring

/-- **The level-set form.** For `φ(t, r) = r − ρ(t)` (the signed distance to the circle,
negative inside), on the cut `r = ρ(t)`: `φ_t = α (φ_rr + φ_r / r)`. This is (★) of SPEC §1,
`φ_t = (G_0 − G_1) − Δ_w φ`, with `G_0 = G_1` (flat disk, `u = 0`), `w = α` and
`−Δ_w φ = α ∇²φ = α/r`.

[formula: cut-curve-shortening-circle] -/
theorem circle_levelSet_flow {ρ₀ α t : ℝ} (h : 0 < ρ₀ ^ 2 - 2 * α * t) :
    deriv (fun s => shrinkRadius ρ₀ α t - shrinkRadius ρ₀ α s) t =
      α * (deriv (fun x => deriv (fun y => y - shrinkRadius ρ₀ α t) x) (shrinkRadius ρ₀ α t) +
        deriv (fun y => y - shrinkRadius ρ₀ α t) (shrinkRadius ρ₀ α t) / shrinkRadius ρ₀ α t) := by
  have h1 : ∀ x : ℝ, HasDerivAt (fun y : ℝ => y - shrinkRadius ρ₀ α t) 1 x := fun x =>
    (hasDerivAt_id x).sub_const _
  have e1 : (fun x => deriv (fun y => y - shrinkRadius ρ₀ α t) x) = fun _ => (1 : ℝ) :=
    funext fun x => (h1 x).deriv
  have ht : HasDerivAt (fun s => shrinkRadius ρ₀ α t - shrinkRadius ρ₀ α s)
      (-(-α / shrinkRadius ρ₀ α t)) t := (hasDerivAt_shrinkRadius h).const_sub _
  rw [ht.deriv, e1, (h1 _).deriv, deriv_const]
  ring

/-- **The length energy decreases along the flow**: `d/dt [α · 2πρ(t)] = −(α/ρ)² · 2πρ`, which
is `−Σ g² ℓ` of `firstVariation_levelSet` with `g = α κ = α/ρ` and `ℓ` the length of the circle
(for `ρ > 0`).

[formula: cut-curve-shortening-circle] -/
theorem hasDerivAt_circle_lengthEnergy {ρ₀ α t : ℝ} (h : 0 < ρ₀ ^ 2 - 2 * α * t) :
    HasDerivAt (fun s => α * (2 * π * shrinkRadius ρ₀ α s))
      (-((α / shrinkRadius ρ₀ α t) ^ 2 * (2 * π * shrinkRadius ρ₀ α t))) t := by
  have hρ := (shrinkRadius_pos h).ne'
  refine (((hasDerivAt_shrinkRadius h).const_mul (2 * π)).const_mul α).congr_deriv ?_
  field_simp

/-! ## The same shape derivatives through the oracle's densities

The Rust code evaluates the gradient densities with the oracle's `gradDirichletD`, `gradScaledD`,
`gradPointwiseD` at recovered normal derivatives `d ≈ ∂u/∂n`, `e ≈ ∂p/∂n` (SPEC §3.6). On the
closed-form solutions, with the exact derivatives in place of the recovered ones, these densities
integrated over the boundary are the derivatives of the energies. -/

/-- Disk, Dirichlet energy: `dJ/dρ = gradDirichletD(u′(ρ)) · 2πρ`.

[formula: cut-gradient-dirichlet] -/
theorem hasDerivAt_disk_dirichletEnergy_density (K ρ : ℝ) :
    HasDerivAt (fun ρ' : ℝ => ∫ r in (0 : ℝ)..ρ', deriv (diskU K ρ') r ^ 2 * (2 * π * r))
      (Ddg.gradDirichletD (deriv (diskU K ρ) ρ) * (2 * π * ρ)) ρ := by
  rw [gradDirichletD_eq]
  exact hasDerivAt_disk_dirichletEnergy K ρ

/-- Disk, scaled Dirichlet energy with a constant weight `c` (adjoint `p = −2c·u`):
`d(cJ)/dρ = gradScaledD c u′(ρ) p′(ρ) · 2πρ`.

[formula: cut-gradient-scaled] -/
theorem hasDerivAt_disk_scaledEnergy_density (K c ρ : ℝ) :
    HasDerivAt
      (fun ρ' : ℝ => c * ∫ r in (0 : ℝ)..ρ', deriv (diskU K ρ') r ^ 2 * (2 * π * r))
      (Ddg.gradScaledD c (deriv (diskU K ρ) ρ) (deriv (fun r => -2 * c * diskU K ρ r) ρ) *
        (2 * π * ρ)) ρ := by
  rw [gradScaledD_eq]
  exact hasDerivAt_disk_scaledEnergy K c ρ

/-- Disk, Hencky energy `b(u) = u²` (`b(0) = 0`, adjoint `diskP`):
`dJ/dρ = gradPointwiseD b(0) u′(ρ) p′(ρ) · 2πρ`.

[formula: cut-gradient-pointwise] -/
theorem hasDerivAt_disk_henckyEnergy_density (K ρ : ℝ) :
    HasDerivAt (fun ρ' : ℝ => ∫ r in (0 : ℝ)..ρ', diskU K ρ' r ^ 2 * (2 * π * r))
      (Ddg.gradPointwiseD ((0 : ℝ) ^ 2) (deriv (diskU K ρ) ρ) (deriv (diskP K ρ) ρ) *
        (2 * π * ρ)) ρ := by
  rw [gradPointwiseD_eq]
  have h := hasDerivAt_disk_henckyEnergy K ρ
  rw [diskU_boundary] at h
  exact h

/-- Spherical cap, Dirichlet energy: `dJ/dψ₀ = gradDirichletD(u′(ψ₀)) · 2π sin ψ₀`.

[formula: cut-gradient-dirichlet] -/
theorem hasDerivAt_capEnergy_density {ψ₀ : ℝ} (h : ψ₀ ∈ Set.Ioo (-π) π) :
    HasDerivAt (fun ψ₁ : ℝ => ∫ ψ in (0 : ℝ)..ψ₁, deriv (capU ψ₁) ψ ^ 2 * (2 * π * Real.sin ψ))
      (Ddg.gradDirichletD (deriv (capU ψ₀) ψ₀) * (2 * π * Real.sin ψ₀)) ψ₀ := by
  rw [gradDirichletD_eq]
  exact hasDerivAt_capEnergy_integral h

end Cuts

end DDG
