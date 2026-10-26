import DdgProofs.RepulsiveCurves.Basic

/-!
# Repulsive Curves: the tangent-point kernel (paper Eq. 3; SPEC §3.1, §3.3, D2)

```
k^α_β(p, q, T) = |T × (p − q)|^α / |p − q|^β
```

with real exponents (`Real.rpow`). For a unit `T`, `|T × d|² = |d|² − (T·d)² = |P_T d|²`; the kernel
is even in `T`, symmetric in `p, q`, invariant under translations and orthogonal maps (rotations
and reflections) applied to all three arguments, and homogeneous: `k(cp, cq, T) = c^{α−β} k(p, q, T)`
for `c > 0`. `k^0_β = |p − q|^{−β}` and `k²₄ · |d|^{−(2σ+1)} = k²_{2σ+5}` (SPEC D2: the exponents of
Eq. 32 and of the trapezoid weight of Sec. 5.2.3 agree).

In Lean `x ^ y` for `x = 0` is `0` (`y ≠ 0`) or `1` (`y = 0`) and `a / 0 = 0`, so the kernel is a
total function; "finite" in the SPEC corresponds to `p ≠ q` here (`Energy.lean`).
-/

noncomputable section

namespace DdgProofs

namespace RC

open Matrix Real

/-- The tangent-point kernel `k^α_β(p, q, T) = |T × (p − q)|^α / |p − q|^β` (paper Eq. 3).

[formula: rc-kernel] -/
def kernel (α β : ℝ) (p q T : V3) : ℝ :=
  vnorm (crossProduct T (p - q)) ^ α / vnorm (p - q) ^ β

/-- Lagrange's identity `|T × d|² = |T|²|d|² − (T·d)²`. -/
theorem cross_vnorm_sq (T d : V3) :
    vnorm (crossProduct T d) ^ 2 = (T ⬝ᵥ T) * (d ⬝ᵥ d) - (T ⬝ᵥ d) ^ 2 := by
  rw [vnorm_sq, cross_dot_cross, dotProduct_comm d T]
  ring

/-- `|P_T d|² = |d|² − (T·d)²` for a unit `T`. -/
theorem proj_vnorm_sq {T : V3} (hT : T ⬝ᵥ T = 1) (d : V3) :
    vnorm (proj T d) ^ 2 = d ⬝ᵥ d - (T ⬝ᵥ d) ^ 2 := by
  rw [vnorm_sq, proj_dot, dot_proj_self hT, dot_proj, dotProduct_comm d T]
  ring

/-- **`|T × d|² = |d|² − (T·d)²`** and `|T × d| = |P_T d|` for a unit `T` (SPEC §3.1).

[formula: rc-kernel] -/
theorem cross_vnorm_sq_unit {T : V3} (hT : T ⬝ᵥ T = 1) (d : V3) :
    vnorm (crossProduct T d) ^ 2 = vnorm d ^ 2 - (T ⬝ᵥ d) ^ 2 ∧
      vnorm (crossProduct T d) = vnorm (proj T d) := by
  refine ⟨by rw [cross_vnorm_sq, hT, vnorm_sq, one_mul], ?_⟩
  refine (pow_left_inj₀ (vnorm_nonneg _) (vnorm_nonneg _) two_ne_zero).mp ?_
  rw [cross_vnorm_sq, proj_vnorm_sq hT, hT, one_mul]

theorem cross_vnorm_le {T : V3} (hT : T ⬝ᵥ T = 1) (d : V3) :
    vnorm (crossProduct T d) ≤ vnorm d := by
  refine (pow_le_pow_iff_left₀ (vnorm_nonneg _) (vnorm_nonneg _) two_ne_zero).mp ?_
  rw [(cross_vnorm_sq_unit hT d).1]
  nlinarith [sq_nonneg (T ⬝ᵥ d)]

theorem kernel_nonneg (α β : ℝ) (p q T : V3) : 0 ≤ kernel α β p q T :=
  div_nonneg (Real.rpow_nonneg (vnorm_nonneg _) _) (Real.rpow_nonneg (vnorm_nonneg _) _)

/-- **The kernel is even in `T`** (so nothing depends on the edge orientation).

[formula: rc-kernel] -/
theorem kernel_neg_tangent (α β : ℝ) (p q T : V3) : kernel α β p q (-T) = kernel α β p q T := by
  simp only [kernel, map_neg, LinearMap.neg_apply, vnorm_neg]

/-- **The kernel is symmetric in `p, q`.**

[formula: rc-kernel] -/
theorem kernel_swap (α β : ℝ) (p q T : V3) : kernel α β q p T = kernel α β p q T := by
  simp only [kernel]
  rw [← neg_sub p q, map_neg, vnorm_neg, vnorm_neg]

/-- **Translation invariance.**

[formula: rc-kernel] -/
theorem kernel_translate (α β : ℝ) (p q T a : V3) :
    kernel α β (p + a) (q + a) T = kernel α β p q T := by
  simp only [kernel, add_sub_add_right_eq_sub]

/-- **Invariance under orthogonal maps** (`Rᵀ R = 1`: rotations and reflections) applied to the two
points and the tangent.

[formula: rc-kernel] -/
theorem kernel_orth {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : Rᵀ * R = 1) (α β : ℝ) (p q T : V3) :
    kernel α β (R *ᵥ p) (R *ᵥ q) (R *ᵥ T) = kernel α β p q T := by
  have hd : R *ᵥ p - R *ᵥ q = R *ᵥ (p - q) := (mulVec_sub R p q).symm
  have hc : vnorm (crossProduct (R *ᵥ T) (R *ᵥ (p - q))) = vnorm (crossProduct T (p - q)) := by
    refine (pow_left_inj₀ (vnorm_nonneg _) (vnorm_nonneg _) two_ne_zero).mp ?_
    rw [cross_vnorm_sq, cross_vnorm_sq, orth_dot hR, orth_dot hR, orth_dot hR]
  simp only [kernel, hd, hc, orth_vnorm hR]

/-- **Homogeneity**: `k(cp, cq, T) = c^{α−β} k(p, q, T)` for `c > 0`.

[formula: rc-kernel] -/
theorem kernel_smul {c : ℝ} (hc : 0 < c) (α β : ℝ) (p q T : V3) :
    kernel α β (c • p) (c • q) T = c ^ (α - β) * kernel α β p q T := by
  simp only [kernel]
  rw [← smul_sub, map_smul, vnorm_smul, vnorm_smul, abs_of_pos hc,
    Real.mul_rpow hc.le (vnorm_nonneg _), Real.mul_rpow hc.le (vnorm_nonneg _), mul_div_mul_comm,
    ← Real.rpow_sub hc]

/-- `k^0_β(p, q, T) = |p − q|^{−β}` (the kernel `k⁰_{2σ+1}` of B.3.4 is the plain power). -/
theorem kernel_zero_alpha (β : ℝ) (p q T : V3) : kernel 0 β p q T = vnorm (p - q) ^ (-β) := by
  rw [kernel, Real.rpow_zero, Real.rpow_neg (vnorm_nonneg _), one_div]

/-- **The quadrature identity of SPEC D2**: `k²₄(p, q, T) · |p − q|^{−(2σ+1)} = k²_{2σ+5}(p, q, T)`
for `p ≠ q`. The kernels of the trapezoid weight (paper Sec. 5.2.3) and of Eq. 32 are the same;
what differs is only the quadrature (D2).

[formula: rc-kernel-matrix-form] -/
theorem kernel_sq_quadrature (σ : ℝ) {p q : V3} (hpq : p ≠ q) (T : V3) :
    kernel 2 4 p q T * vnorm (p - q) ^ (-(2 * σ + 1)) = kernel 2 (2 * σ + 5) p q T := by
  have hρ : 0 < vnorm (p - q) := vnorm_pos_iff.mpr (sub_ne_zero.mpr hpq)
  have h5 : vnorm (p - q) ^ (2 * σ + 5) = vnorm (p - q) ^ (4 : ℝ) * vnorm (p - q) ^ (2 * σ + 1) := by
    rw [← Real.rpow_add hρ]; ring_nf
  have h4 : 0 < vnorm (p - q) ^ (4 : ℝ) := Real.rpow_pos_of_pos hρ _
  have h1 : 0 < vnorm (p - q) ^ (2 * σ + 1) := Real.rpow_pos_of_pos hρ _
  simp only [kernel]
  rw [h5, Real.rpow_neg hρ.le]
  field_simp

end RC

end DdgProofs
