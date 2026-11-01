import DdgProofs.RepulsiveSurfaces.Geometry
import DdgProofs.RepulsiveCurves.Kernel

/-!
# Repulsive Surfaces: the tangent-point kernel and its gradient (paper Eq. 16; SPEC §3.1, §3.2)

```
K(X, N; Y) = |N·(X − Y)|^α / |X − Y|^β
```

with real exponents (`Real.rpow`; Lean's `0 ^ y = 0` for `y ≠ 0` and `a / 0 = 0` make it a total
function, "finite" in the SPEC is `X ≠ Y` here).

* `rs-kernel`: `K = |ν|^α r^{−β}`; even in `N`; invariant under translations and orthogonal maps
  `R` (also `N ↦ −RN`, which is what an improper `R` does to a face normal); homogeneous,
  `K(cX, N; cY) = c^{α−β} K(X, N; Y)` (`c > 0`); with `P = NNᵀ` and `|N| = 1`,
  `|Pv|² = vᵀPv = (N·v)²` and `(vᵀPv)^{α/2} = |N·v|^α` (the reference's projector form).
* `rs-kernel-gradient`: with `v = X − Y`, `ν = N·v`, `r = |v|`,
  `∂K/∂v = α|ν|^{α−2}ν r^{−β} N − β|ν|^α r^{−β−2} v`, `∂K/∂N = α|ν|^{α−2}ν r^{−β} v`,
  `∂K/∂Y = −∂K/∂v`: the derivative along any differentiable curves `X(s), N(s), Y(s)` with
  `X ≠ Y` at `s = t`, for `α > 1`, **including `ν = 0`** (Mathlib's `hasDerivAt_abs_rpow`:
  `|x|^α` is differentiable with derivative `α|x|^{α−2}x` everywhere for `α > 1`). `N` need not be a
  unit vector here.
-/

noncomputable section

namespace DdgProofs

namespace RS

open Matrix Real Finset
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff
  dot_self_nonneg unit_dot_self proj hasDerivAt_dot hasDerivAt_cross hasDerivAt_vnorm
  hasDerivAt_unit orth_dot orth_vnorm vnorm_smul vnorm_neg)

/-- The tangent-point kernel `K(X, N; Y) = |N·(X − Y)|^α / |X − Y|^β` (paper Eq. 16, in the
reference's two-exponent form).

[formula: rs-kernel] -/
def kernel (α β : ℝ) (X N Y : V3) : ℝ := |N ⬝ᵥ (X - Y)| ^ α / vnorm (X - Y) ^ β

theorem kernel_nonneg (α β : ℝ) (X N Y : V3) : 0 ≤ kernel α β X N Y :=
  div_nonneg (Real.rpow_nonneg (abs_nonneg _) _) (Real.rpow_nonneg (vnorm_nonneg _) _)

/-- **The SPEC's form** `K = |ν|^α r^{−β}`.

[formula: rs-kernel] -/
theorem kernel_eq_mul_rpow_neg (α β : ℝ) (X N Y : V3) :
    kernel α β X N Y = |N ⬝ᵥ (X - Y)| ^ α * vnorm (X - Y) ^ (-β) := by
  rw [kernel, Real.rpow_neg (vnorm_nonneg _), div_eq_mul_inv]

/-- **The kernel is even in `N`** (the sign of a face normal is irrelevant).

[formula: rs-kernel] -/
theorem kernel_neg_normal (α β : ℝ) (X N Y : V3) : kernel α β X (-N) Y = kernel α β X N Y := by
  simp only [kernel, neg_dotProduct, abs_neg]

/-- **Translation invariance.**

[formula: rs-kernel] -/
theorem kernel_translate (α β : ℝ) (X N Y a : V3) :
    kernel α β (X + a) N (Y + a) = kernel α β X N Y := by
  simp only [kernel, add_sub_add_right_eq_sub]

/-- **Invariance under orthogonal maps** (`Rᵀ R = 1`) applied to both points and to the normal,
with either sign of the normal.

[formula: rs-kernel] -/
theorem kernel_orth {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : Rᵀ * R = 1) (α β : ℝ) (X N Y : V3)
    (σ : ℝ) (hσ : σ = 1 ∨ σ = -1) :
    kernel α β (R *ᵥ X) (σ • R *ᵥ N) (R *ᵥ Y) = kernel α β X N Y := by
  have hd : R *ᵥ X - R *ᵥ Y = R *ᵥ (X - Y) := (mulVec_sub R X Y).symm
  have hσ' : |σ| = 1 := by rcases hσ with h | h <;> simp [h]
  simp only [kernel, hd, smul_dotProduct, orth_dot hR, orth_vnorm hR, smul_eq_mul, abs_mul, hσ',
    one_mul]

/-- **Homogeneity**: `K(cX, N; cY) = c^{α−β} K(X, N; Y)` for `c > 0`.

[formula: rs-kernel] -/
theorem kernel_smul {c : ℝ} (hc : 0 < c) (α β : ℝ) (X N Y : V3) :
    kernel α β (c • X) N (c • Y) = c ^ (α - β) * kernel α β X N Y := by
  simp only [kernel]
  rw [← smul_sub, dotProduct_smul, smul_eq_mul, abs_mul, abs_of_pos hc, vnorm_smul, abs_of_pos hc,
    Real.mul_rpow hc.le (abs_nonneg _), Real.mul_rpow hc.le (vnorm_nonneg _), mul_div_mul_comm,
    ← Real.rpow_sub hc]

/-- **The projector form** (SPEC §3.1, §8.3): for a unit `N` and `P = NNᵀ`,
`P v = (N·v) N`, `|Pv|² = vᵀPv = (N·v)²`, and `(vᵀPv)^{α/2} = |N·v|^α` (`α` real).

[formula: rs-kernel] -/
theorem projector_form {N : V3} (hN : N ⬝ᵥ N = 1) (v : V3) (α : ℝ) :
    vecMulVec N N *ᵥ v = (N ⬝ᵥ v) • N ∧
      vnorm (vecMulVec N N *ᵥ v) ^ 2 = (N ⬝ᵥ v) ^ 2 ∧
      v ⬝ᵥ (vecMulVec N N *ᵥ v) = (N ⬝ᵥ v) ^ 2 ∧
      (v ⬝ᵥ (vecMulVec N N *ᵥ v)) ^ (α / 2) = |N ⬝ᵥ v| ^ α := by
  have hP : vecMulVec N N *ᵥ v = (N ⬝ᵥ v) • N := by
    ext i
    simp [vecMulVec, mulVec, dotProduct, Finset.mul_sum, mul_comm, mul_left_comm]
  have hq : v ⬝ᵥ (vecMulVec N N *ᵥ v) = (N ⬝ᵥ v) ^ 2 := by
    rw [hP, dotProduct_smul, smul_eq_mul, dotProduct_comm v N]; ring
  refine ⟨hP, ?_, hq, ?_⟩
  · rw [vnorm_sq, hP, smul_dotProduct, dotProduct_smul, hN, smul_eq_mul, smul_eq_mul]; ring
  · rw [hq, ← sq_abs, ← Real.rpow_natCast, ← Real.rpow_mul (abs_nonneg _)]
    congr 1; push_cast; ring

/-! ## The gradient -/

/-- `∂K/∂v = α|ν|^{α−2}ν r^{−β} N − β|ν|^α r^{−β−2} v` (`v = X − Y`, `ν = N·v`, `r = |v|`);
`∂K/∂X = ∂K/∂v`, `∂K/∂Y = −∂K/∂v`.

[formula: rs-kernel-gradient] -/
def kgradV (α β : ℝ) (N v : V3) : V3 :=
  (α * |N ⬝ᵥ v| ^ (α - 2) * (N ⬝ᵥ v) / vnorm v ^ β) • N
    - (β * |N ⬝ᵥ v| ^ α / vnorm v ^ (β + 2)) • v

/-- `∂K/∂N = α|ν|^{α−2}ν r^{−β} v`.

[formula: rs-kernel-gradient] -/
def kgradN (α β : ℝ) (N v : V3) : V3 := (α * |N ⬝ᵥ v| ^ (α - 2) * (N ⬝ᵥ v) / vnorm v ^ β) • v

/-- **The kernel's derivative along curves** (SPEC §3.2): for `α > 1` and differentiable curves
`X, N, Y` with `X(t) ≠ Y(t)`,
`d/ds K(X, N; Y) = ∂K/∂v · (X' − Y') + ∂K/∂N · N'`. No condition on `ν`: at `ν = 0` the
`|ν|^{α−2}ν` factor is `0` and is the derivative of `|ν|^α` there.

[formula: rs-kernel-gradient] -/
theorem hasDerivAt_kernel {α β : ℝ} (hα : 1 < α) {X N Y : ℝ → V3} {X' N' Y' : V3} {t : ℝ}
    (hX : HasDerivAt X X' t) (hN : HasDerivAt N N' t) (hY : HasDerivAt Y Y' t)
    (hXY : X t ≠ Y t) :
    HasDerivAt (fun s => kernel α β (X s) (N s) (Y s))
      (kgradV α β (N t) (X t - Y t) ⬝ᵥ (X' - Y') + kgradN α β (N t) (X t - Y t) ⬝ᵥ N') t := by
  have hd : HasDerivAt (fun s => X s - Y s) (X' - Y') t := hX.sub hY
  have hd0 : X t - Y t ≠ 0 := sub_ne_zero.mpr hXY
  have hρ : 0 < vnorm (X t - Y t) := vnorm_pos_iff.mpr hd0
  have hν := hasDerivAt_dot hN hd
  have hA := (hasDerivAt_abs_rpow (N t ⬝ᵥ (X t - Y t)) hα).comp t hν
  have hr := hasDerivAt_vnorm hd hd0
  have hrb := hr.rpow_const (p := β) (Or.inl hρ.ne')
  have hk := hA.div hrb (Real.rpow_pos_of_pos hρ β).ne'
  refine hk.congr_deriv ?_
  set d := X t - Y t with hd_def
  set r := vnorm d with hr_def
  set ν := N t ⬝ᵥ d with hν_def
  have hB : 0 < r ^ β := Real.rpow_pos_of_pos hρ β
  have r1 : r ^ (β - 1) = r ^ β / r := Real.rpow_sub_one hρ.ne' β
  have r2 : r ^ (β + 2) = r ^ β * r ^ 2 := by
    rw [Real.rpow_add hρ, Real.rpow_two]
  simp only [kgradV, kgradN, sub_dotProduct, smul_dotProduct, smul_eq_mul, Function.comp_apply]
  rw [r1, r2, dotProduct_comm N' d]
  simp only [← hd_def, ← hr_def, ← hν_def]
  field_simp
  ring

end RS

end DdgProofs
