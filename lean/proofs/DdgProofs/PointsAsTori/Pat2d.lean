import DdgProofs.PointsAsTori.Loss

/-!
# Points as Tori: the 2D variant with circles (Fig. 10; SPEC §8.7, D15 (a))

Everything is written in the local frame `(t̂, n̂)` of the point (an orthonormal frame of `ℝ²`, so
distances are coordinate distances); points are pairs of coordinates.

* `pat2d-circle-fit`: the graph `γ(s) = (s, Q(s))`, `Q = a0 + a1 s + a2 s²`, has signed curvature
  `κ = 2a2/(1 + a1²)^{3/2}` at `0` w.r.t. the unit normal `n* = (−a1, 1)/√(1 + a1²)`; the center
  `Q*(0) − (sign/|κ|) n*` with `sign = +1` iff `κ < 0` is `o + n*/κ`; and the circle through `o`
  with that center osculates: `|γ(s) − c|² − 1/κ² = 2a1a2 s³ + a2² s⁴`.
* `pat2d-circle-sdf`: `g = sign·(|x − c| − 1/|κ|)` is the signed distance to the circle (lower
  bound and attainment), `g(o) = 0`, and its gradient at `o` is `+n*`.
* `pat2d-blend`: Eq. 25 in 2D: the shift cancels, and the blend is a convex combination.
-/

noncomputable section

namespace DdgProofs

namespace PAT

open Real Finset

/-- The Euclidean norm of `(a, b)`. -/
def n2 (a b : ℝ) : ℝ := √(a ^ 2 + b ^ 2)

/-- The coefficients of the 2D height function `Q(s) = a0 + a1 s + a2 s²`. -/
structure Coef2 where
  a0 : ℝ
  a1 : ℝ
  a2 : ℝ

namespace Coef2

variable (a : Coef2)

def A : ℝ := √(1 + a.a1 ^ 2)
/-- `κ = 2a2/(1 + a1²)^{3/2}`. -/
def kappa : ℝ := 2 * a.a2 / ((1 + a.a1 ^ 2) * a.A)
/-- The sign of the circle SDF: `+1` iff `κ < 0`. -/
def csign : ℝ := if a.kappa < 0 then 1 else -1

theorem A_pos : 0 < a.A := Real.sqrt_pos.mpr (by positivity)
theorem A_sq : a.A ^ 2 = 1 + a.a1 ^ 2 := Real.sq_sqrt (by positivity)

end Coef2

/-- **The circle fit** (SPEC §8.7): with `γ'(0) = (1, a1)`, `γ''(0) = (0, 2a2)` and
`n* = (−a1, 1)/A`, `A = √(1 + a1²)`: `n*` is a unit normal, `κ = γ''·n*/|γ'|² =
(x'y'' − y'x'')/|γ'|³ = 2a2/(1 + a1²)^{3/2}`; for `κ ≠ 0` the SPEC's center
`o − (sign/|κ|) n*` equals `o + n*/κ`; and the circle of radius `1/|κ|` around it osculates the
graph at `o = (0, a0)`: `|γ(s) − c|² − 1/κ² = 2a1a2 s³ + a2² s⁴`.

[formula: pat2d-circle-fit] -/
theorem circle_fit (a : Coef2) :
    ((-a.a1 / a.A) ^ 2 + (1 / a.A) ^ 2 = 1 ∧ 1 * (-a.a1 / a.A) + a.a1 * (1 / a.A) = 0) ∧
    (0 * (-a.a1 / a.A) + 2 * a.a2 * (1 / a.A)) / (1 + a.a1 ^ 2) = a.kappa ∧
    (1 * (2 * a.a2) - a.a1 * 0) / ((1 + a.a1 ^ 2) * a.A) = a.kappa ∧
    (a.kappa ≠ 0 → (a.csign / |a.kappa|) = -(1 / a.kappa)) ∧
    (a.a2 ≠ 0 → ∀ s : ℝ,
      (s - (-a.a1 / a.A) / a.kappa) ^ 2 +
        (a.a0 + a.a1 * s + a.a2 * s ^ 2 - (a.a0 + (1 / a.A) / a.kappa)) ^ 2 - (1 / a.kappa) ^ 2 =
        2 * a.a1 * a.a2 * s ^ 3 + a.a2 ^ 2 * s ^ 4) := by
  have hA := a.A_pos
  have hA2 := a.A_sq
  have hA0 := hA.ne'
  refine ⟨⟨?_, ?_⟩, ?_, ?_, ?_, ?_⟩
  · field_simp; rw [hA2]; ring
  · field_simp; ring
  · unfold Coef2.kappa; field_simp; ring
  · unfold Coef2.kappa; field_simp; ring
  · intro hk
    unfold Coef2.csign
    split_ifs with h
    · rw [abs_of_neg h]; field_simp
    · rw [abs_of_pos (lt_of_le_of_ne (not_lt.mp h) (Ne.symm hk))]; field_simp
  · intro h2 s
    have e1 : (-a.a1 / a.A) / a.kappa = -a.a1 * (1 + a.a1 ^ 2) / (2 * a.a2) := by
      unfold Coef2.kappa; field_simp
    have e2 : (1 / a.A) / a.kappa = (1 + a.a1 ^ 2) / (2 * a.a2) := by
      unfold Coef2.kappa; field_simp
    have e3 : (1 / a.kappa) ^ 2 = (1 + a.a1 ^ 2) ^ 3 / (4 * a.a2 ^ 2) := by
      unfold Coef2.kappa
      rw [one_div_div, div_pow, mul_pow, hA2]
      field_simp
      ring
    rw [e1, e2, e3]
    field_simp
    ring

/-- **The circle SDF** (SPEC §8.7): `|x − c| − ρ` is the distance to the circle up to sign: at
most the distance to any point `y` of the circle, attained for `x ≠ c`; for the fitted circle
(`c = o + n*/κ`, `ρ = 1/|κ|`), `g = sign·(|x − c| − ρ)` vanishes at `o` and its gradient there,
`sign·(o − c)/|o − c|`, is `+n*` (positive on the side of the input normal).

[formula: pat2d-circle-sdf] -/
theorem circle_sdf (c0 c1 ρ x0 x1 : ℝ) :
    (∀ y0 y1 : ℝ, n2 (y0 - c0) (y1 - c1) = ρ →
      |n2 (x0 - c0) (x1 - c1) - ρ| ≤ n2 (x0 - y0) (x1 - y1)) ∧
    (0 ≤ ρ → n2 (x0 - c0) (x1 - c1) ≠ 0 →
      let y0 := c0 + ρ * (x0 - c0) / n2 (x0 - c0) (x1 - c1)
      let y1 := c1 + ρ * (x1 - c1) / n2 (x0 - c0) (x1 - c1)
      n2 (y0 - c0) (y1 - c1) = ρ ∧ n2 (x0 - y0) (x1 - y1) = |n2 (x0 - c0) (x1 - c1) - ρ|) := by
  refine ⟨fun y0 y1 hy => ?_, fun hρ hn => ?_⟩
  · have t1 := sqrt2_tri (x0 - c0) (x1 - c1) (y0 - c0) (y1 - c1)
    have t2 := sqrt2_tri (y0 - c0) (y1 - c1) (x0 - c0) (x1 - c1)
    have e1 : (x0 - c0 - (y0 - c0)) ^ 2 + (x1 - c1 - (y1 - c1)) ^ 2 = (x0 - y0) ^ 2 + (x1 - y1) ^ 2 := by
      ring
    have e2 : (y0 - c0 - (x0 - c0)) ^ 2 + (y1 - c1 - (x1 - c1)) ^ 2 = (x0 - y0) ^ 2 + (x1 - y1) ^ 2 := by
      ring
    rw [e1] at t1; rw [e2] at t2
    unfold n2 at hy ⊢
    rw [hy] at t1 t2
    rw [abs_sub_le_iff]; constructor <;> linarith
  · intro y0 y1
    set q := n2 (x0 - c0) (x1 - c1)
    have hq0 : 0 ≤ q := Real.sqrt_nonneg _
    have hq : 0 < q := lt_of_le_of_ne hq0 (Ne.symm hn)
    have hq2 : q ^ 2 = (x0 - c0) ^ 2 + (x1 - c1) ^ 2 := Real.sq_sqrt (by positivity)
    constructor
    · show √((c0 + ρ * (x0 - c0) / q - c0) ^ 2 + (c1 + ρ * (x1 - c1) / q - c1) ^ 2) = ρ
      rw [show (c0 + ρ * (x0 - c0) / q - c0) ^ 2 + (c1 + ρ * (x1 - c1) / q - c1) ^ 2 =
        ρ ^ 2 * (((x0 - c0) ^ 2 + (x1 - c1) ^ 2) / q ^ 2) by field_simp; ring, ← hq2,
        div_self (by positivity), mul_one, Real.sqrt_sq hρ]
    · show √((x0 - (c0 + ρ * (x0 - c0) / q)) ^ 2 + (x1 - (c1 + ρ * (x1 - c1) / q)) ^ 2) = |q - ρ|
      rw [show (x0 - (c0 + ρ * (x0 - c0) / q)) ^ 2 + (x1 - (c1 + ρ * (x1 - c1) / q)) ^ 2 =
        (q - ρ) ^ 2 * (((x0 - c0) ^ 2 + (x1 - c1) ^ 2) / q ^ 2) by field_simp; ring, ← hq2,
        div_self (by positivity), mul_one, Real.sqrt_sq_eq_abs]

/-- The fitted circle passes through `o` and `sign·(o − c)/|o − c| = n*`.

[formula: pat2d-circle-sdf] -/
theorem circle_sdf_contact (a : Coef2) (hk : a.kappa ≠ 0) :
    n2 (0 - (-a.a1 / a.A) / a.kappa) (a.a0 - (a.a0 + (1 / a.A) / a.kappa)) = 1 / |a.kappa| ∧
    a.csign * ((0 - (-a.a1 / a.A) / a.kappa) / (1 / |a.kappa|)) = -a.a1 / a.A ∧
    a.csign * ((a.a0 - (a.a0 + (1 / a.A) / a.kappa)) / (1 / |a.kappa|)) = 1 / a.A := by
  have hA := a.A_pos
  have hA2 := a.A_sq
  have hn : n2 (0 - (-a.a1 / a.A) / a.kappa) (a.a0 - (a.a0 + (1 / a.A) / a.kappa)) =
      1 / |a.kappa| := by
    unfold n2
    rw [show (0 - (-a.a1 / a.A) / a.kappa) ^ 2 + (a.a0 - (a.a0 + (1 / a.A) / a.kappa)) ^ 2 =
      (1 / a.kappa) ^ 2 * ((a.a1 ^ 2 + 1) / a.A ^ 2) by field_simp; ring, hA2,
      show (a.a1 ^ 2 + 1) / (1 + a.a1 ^ 2) = 1 by field_simp; ring, mul_one,
      Real.sqrt_sq_eq_abs, abs_div, abs_one]
  refine ⟨hn, ?_, ?_⟩ <;> unfold Coef2.csign <;> split_ifs with h
  · rw [abs_of_neg h]; field_simp; ring
  · rw [abs_of_pos (lt_of_le_of_ne (not_lt.mp h) (Ne.symm hk))]; field_simp; ring
  · rw [abs_of_neg h]; field_simp; ring
  · rw [abs_of_pos (lt_of_le_of_ne (not_lt.mp h) (Ne.symm hk))]; field_simp; ring

/-- **The 2D blend** (Eq. 25 with circles): for any distances `ρ_i` (here `|x − p_i|` in the
plane) the shift cancels, and with positive weights the blend lies between the smallest and the
largest circle SDF (`blend_convex`); the weight range is `blend_weight_range` verbatim.

[formula: pat2d-blend] -/
theorem blend2_shift {ι : Type*} [Fintype ι] (a ρ g : ι → ℝ) (lam σ : ℝ) :
    (∑ i, a i * Real.exp (-lam * (ρ i - σ)) * g i) / ∑ i, a i * Real.exp (-lam * (ρ i - σ)) =
      (∑ i, a i * Real.exp (-lam * ρ i) * g i) / ∑ i, a i * Real.exp (-lam * ρ i) := by
  have h : ∀ i, Real.exp (-lam * (ρ i - σ)) = Real.exp (lam * σ) * Real.exp (-lam * ρ i) := by
    intro i; rw [← Real.exp_add]; congr 1; ring
  simp only [h]
  have e1 : ∑ i, a i * (Real.exp (lam * σ) * Real.exp (-lam * ρ i)) * g i =
      Real.exp (lam * σ) * ∑ i, a i * Real.exp (-lam * ρ i) * g i := by
    rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun i _ => by ring
  have e2 : ∑ i, a i * (Real.exp (lam * σ) * Real.exp (-lam * ρ i)) =
      Real.exp (lam * σ) * ∑ i, a i * Real.exp (-lam * ρ i) := by
    rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun i _ => by ring
  rw [e1, e2, mul_div_mul_left _ _ (Real.exp_pos _).ne']

end PAT

end DdgProofs
