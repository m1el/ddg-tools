import Mathlib

/-!
# Vectors of `ℝ³`, edge quantities, derivatives along curves (Repulsive Curves, SPEC §2.1)

"Repulsive Curves", Yu, Schumacher, Crane, ACM TOG 40(2), 2021 ("the paper"); contract
`crates/ddg-repulsive-curves/SPEC.md` ("SPEC"), §2.1.

Vectors are `Fin 3 → ℝ` with Mathlib's `dotProduct` (`⬝ᵥ`) and `crossProduct` (`×₃`); `vnorm` is
the Euclidean norm `√(v·v)` (Mathlib's norm on `Fin 3 → ℝ` is the sup norm). A curve network has
vertices `V`, edges `E` and endpoints `ep : E → Fin 2 → V` (`i₁ = ep I 0`, `i₂ = ep I 1`, the
orientation the caller gives); positions are `γ : V → ℝ³`.

```
d_I = γ_{i₂} − γ_{i₁},   ℓ_I = |d_I|,   T_I = d_I / ℓ_I,   x_I = (γ_{i₁} + γ_{i₂}) / 2
m_v = ½ Σ_{I ∋ v} ℓ_I
```

**Derivatives** are taken along arbitrary differentiable curves of positions `s ↦ Γ s`
(`HasDerivAt` of each vertex position, a curve in the normed space `Fin 3 → ℝ`), as in the
developability modules: a statement "`d/ds F(Γ s) = Σ_v g_v · Γ'_v`" contains the gradient `g` and
the chain rule along any differentiable motion. The lemmas here give the derivatives of `ℓ_I`, `T_I`
and `x_I`:

```
ℓ_I' = T_I · (Γ'_{i₂} − Γ'_{i₁}),   T_I' = P_{T_I}(Γ'_{i₂} − Γ'_{i₁}) / ℓ_I,   x_I' = (Γ'_{i₁} + Γ'_{i₂}) / 2
```

with `P_T w = w − (T·w) T`.

(The dual lengths `m_v` are in `DdgProofs.Repulsive.Curves.DualLengths`.)
-/

noncomputable section

namespace DDG


open Matrix Real Filter Topology

/-- Vectors of `ℝ³`. -/
abbrev V3 := Fin 3 → ℝ

/-! ## The Euclidean norm -/

theorem dot_self_nonneg (v : V3) : 0 ≤ v ⬝ᵥ v := by
  simp only [dotProduct, Fin.sum_univ_three]
  exact add_nonneg (add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)) (mul_self_nonneg _)

/-- The Euclidean norm `|v| = √(v·v)`. -/
def vnorm (v : V3) : ℝ := √(v ⬝ᵥ v)

theorem vnorm_nonneg (v : V3) : 0 ≤ vnorm v := Real.sqrt_nonneg _

theorem vnorm_sq (v : V3) : vnorm v ^ 2 = v ⬝ᵥ v := Real.sq_sqrt (dot_self_nonneg v)

theorem vnorm_mul_self (v : V3) : vnorm v * vnorm v = v ⬝ᵥ v :=
  Real.mul_self_sqrt (dot_self_nonneg v)

theorem vnorm_eq_zero_iff {v : V3} : vnorm v = 0 ↔ v = 0 := by
  rw [vnorm, Real.sqrt_eq_zero (dot_self_nonneg v), dotProduct_self_eq_zero]

theorem vnorm_pos_iff {v : V3} : 0 < vnorm v ↔ v ≠ 0 := by
  rw [lt_iff_le_and_ne, ne_comm, Ne, vnorm_eq_zero_iff]
  exact ⟨fun h => h.2, fun h => ⟨vnorm_nonneg v, h⟩⟩

theorem vnorm_smul (c : ℝ) (v : V3) : vnorm (c • v) = |c| * vnorm v := by
  rw [vnorm, vnorm, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul, ← mul_assoc,
    Real.sqrt_mul (mul_self_nonneg c), Real.sqrt_mul_self_eq_abs]

theorem vnorm_neg (v : V3) : vnorm (-v) = vnorm v := by
  simp [vnorm]

theorem vnorm_sub_comm (p q : V3) : vnorm (p - q) = vnorm (q - p) := by
  rw [← neg_sub, vnorm_neg]

theorem vnorm_zero : vnorm 0 = 0 := by simp [vnorm]

/-- A vector divided by its length has length one. -/
theorem unit_dot_self {v : V3} (hv : v ≠ 0) : ((vnorm v)⁻¹ • v) ⬝ᵥ ((vnorm v)⁻¹ • v) = 1 := by
  have h := vnorm_pos_iff.mpr hv
  rw [smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul, ← vnorm_mul_self]
  field_simp

/-- Orthogonal matrices (`Rᵀ R = 1`, rotations and reflections) preserve dot products. -/
theorem orth_dot {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : Rᵀ * R = 1) (u v : V3) :
    (R *ᵥ u) ⬝ᵥ (R *ᵥ v) = u ⬝ᵥ v := by
  rw [dotProduct_mulVec, ← vecMul_transpose, vecMul_vecMul, hR, vecMul_one]

theorem orth_vnorm {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : Rᵀ * R = 1) (v : V3) :
    vnorm (R *ᵥ v) = vnorm v := by
  rw [vnorm, vnorm, orth_dot hR]

/-- The projection onto the plane orthogonal to `T`: `P_T w = w − (T·w) T`. -/
def proj (T w : V3) : V3 := w - (T ⬝ᵥ w) • T

theorem proj_dot (T w u : V3) : proj T w ⬝ᵥ u = w ⬝ᵥ u - (T ⬝ᵥ w) * (T ⬝ᵥ u) := by
  simp [proj, sub_dotProduct, smul_dotProduct]

theorem dot_proj (T w u : V3) : u ⬝ᵥ proj T w = u ⬝ᵥ w - (T ⬝ᵥ w) * (u ⬝ᵥ T) := by
  simp [proj, dotProduct_sub, dotProduct_smul]

/-- `P_T` is symmetric: `(P_T w)·u = w·(P_T u)`. -/
theorem proj_dot_comm (T w u : V3) : proj T w ⬝ᵥ u = w ⬝ᵥ proj T u := by
  rw [proj_dot, dot_proj, dotProduct_comm w T]
  ring

/-- For a unit `T`, `P_T T = 0` and `T · P_T w = 0`. -/
theorem proj_self {T : V3} (hT : T ⬝ᵥ T = 1) : proj T T = 0 := by
  simp [proj, hT]

theorem dot_proj_self {T : V3} (hT : T ⬝ᵥ T = 1) (w : V3) : T ⬝ᵥ proj T w = 0 := by
  rw [dot_proj, hT, dotProduct_comm]
  ring

/-! ## Edge quantities -/

section edges

variable {V E : Type*}

/-- `d_I = γ_{i₂} − γ_{i₁}`. -/
def edgeVec (γ : V → V3) (ep : E → Fin 2 → V) (I : E) : V3 := γ (ep I 1) - γ (ep I 0)

/-- The edge length `ℓ_I = |γ_{i₂} − γ_{i₁}|`. -/
def len (γ : V → V3) (ep : E → Fin 2 → V) (I : E) : ℝ := vnorm (edgeVec γ ep I)

/-- The unit tangent `T_I = (γ_{i₂} − γ_{i₁}) / ℓ_I`. -/
def tangent (γ : V → V3) (ep : E → Fin 2 → V) (I : E) : V3 := (len γ ep I)⁻¹ • edgeVec γ ep I

/-- The midpoint `x_I = (γ_{i₁} + γ_{i₂}) / 2`. -/
def mid (γ : V → V3) (ep : E → Fin 2 → V) (I : E) : V3 := (1 / 2 : ℝ) • (γ (ep I 0) + γ (ep I 1))

/-- The endpoint signs `s₁ = −1`, `s₂ = +1` (`(−1)^{a+b} = s_a s_b`). -/
def sgn : Fin 2 → ℝ := ![-1, 1]

@[simp] theorem sgn_zero : sgn 0 = -1 := rfl
@[simp] theorem sgn_one : sgn 1 = 1 := rfl

theorem sgn_mul_sgn (a b : Fin 2) : sgn a * sgn b = (-1 : ℝ) ^ ((a : ℕ) + (b : ℕ)) := by
  fin_cases a <;> fin_cases b <;> norm_num

theorem sgn_mul_self (a : Fin 2) : sgn a * sgn a = 1 := by
  fin_cases a <;> norm_num

theorem len_nonneg (γ : V → V3) (ep : E → Fin 2 → V) (I : E) : 0 ≤ len γ ep I := vnorm_nonneg _

/-- `ℓ_I > 0` exactly when the two endpoints are at different points. -/
theorem len_pos_iff {γ : V → V3} {ep : E → Fin 2 → V} {I : E} :
    0 < len γ ep I ↔ γ (ep I 0) ≠ γ (ep I 1) := by
  rw [len, vnorm_pos_iff, edgeVec, Ne, Ne, sub_eq_zero, eq_comm]

/-- **`|T_I| = 1`** for an edge of positive length.

[formula: rc-edge-quantities] -/
theorem tangent_dot_self {γ : V → V3} {ep : E → Fin 2 → V} {I : E} (h : 0 < len γ ep I) :
    tangent γ ep I ⬝ᵥ tangent γ ep I = 1 := by
  have hd : edgeVec γ ep I ≠ 0 := vnorm_pos_iff.mp h
  exact unit_dot_self hd

theorem vnorm_tangent {γ : V → V3} {ep : E → Fin 2 → V} {I : E} (h : 0 < len γ ep I) :
    vnorm (tangent γ ep I) = 1 := by
  rw [vnorm, tangent_dot_self h, Real.sqrt_one]

/-- `ℓ_I T_I = d_I` (also for `ℓ_I = 0`, where both sides vanish). -/
theorem len_smul_tangent (γ : V → V3) (ep : E → Fin 2 → V) (I : E) :
    len γ ep I • tangent γ ep I = edgeVec γ ep I := by
  by_cases h : len γ ep I = 0
  · have : edgeVec γ ep I = 0 := vnorm_eq_zero_iff.mp h
    simp [tangent, this]
  · rw [tangent, smul_smul, mul_inv_cancel₀ h, one_smul]

/-- `T_I · d_I = ℓ_I`. -/
theorem tangent_dot_edgeVec (γ : V → V3) (ep : E → Fin 2 → V) (I : E) :
    tangent γ ep I ⬝ᵥ edgeVec γ ep I = len γ ep I := by
  rw [tangent, smul_dotProduct, ← vnorm_mul_self, smul_eq_mul, ← len]
  by_cases h : len γ ep I = 0
  · simp [h]
  · field_simp

/-- **The edge quantities of SPEC §2.1**: `d_I = Σ_a s_a γ_{i_a}`, `ℓ_I T_I = d_I`, `T_I · d_I = ℓ_I`,
`x_I = γ_{i₁} + ½ d_I`, and for `ℓ_I > 0`, `|T_I| = 1`.

[formula: rc-edge-quantities] -/
theorem edge_quantities (γ : V → V3) (ep : E → Fin 2 → V) (I : E) :
    edgeVec γ ep I = ∑ a, sgn a • γ (ep I a) ∧
      len γ ep I • tangent γ ep I = edgeVec γ ep I ∧
      tangent γ ep I ⬝ᵥ edgeVec γ ep I = len γ ep I ∧
      mid γ ep I = γ (ep I 0) + (1 / 2 : ℝ) • edgeVec γ ep I ∧
      (0 < len γ ep I → tangent γ ep I ⬝ᵥ tangent γ ep I = 1) := by
  refine ⟨?_, len_smul_tangent γ ep I, tangent_dot_edgeVec γ ep I, ?_, tangent_dot_self⟩
  · simp [edgeVec, Fin.sum_univ_two, sub_eq_add_neg, add_comm]
  · ext i
    simp [mid, edgeVec]
    ring

end edges

/-! ## Derivatives along curves -/

section deriv

variable {t : ℝ}

/-- Derivative of the dot product of two curves. -/
theorem hasDerivAt_dot {f g : ℝ → V3} {f' g' : V3} (hf : HasDerivAt f f' t)
    (hg : HasDerivAt g g' t) :
    HasDerivAt (fun s => f s ⬝ᵥ g s) (f' ⬝ᵥ g t + f t ⬝ᵥ g') t := by
  have hf' := hasDerivAt_pi.mp hf
  have hg' := hasDerivAt_pi.mp hg
  have := HasDerivAt.sum (u := Finset.univ) (fun i _ => (hf' i).mul (hg' i))
  simp only [dotProduct]
  convert this using 1
  · funext s; simp
  · rw [← Finset.sum_add_distrib]

/-- Derivative of the cross product of two curves. -/
theorem hasDerivAt_cross {f g : ℝ → V3} {f' g' : V3} (hf : HasDerivAt f f' t)
    (hg : HasDerivAt g g' t) :
    HasDerivAt (fun s => crossProduct (f s) (g s)) (crossProduct f' (g t) + crossProduct (f t) g') t := by
  have hf' := hasDerivAt_pi.mp hf
  have hg' := hasDerivAt_pi.mp hg
  rw [hasDerivAt_pi]
  intro i
  fin_cases i
  · simp only [cross_apply, Fin.zero_eta, Pi.add_apply, Matrix.cons_val_zero]
    convert ((hf' 1).mul (hg' 2)).sub ((hf' 2).mul (hg' 1)) using 1; ring
  · simp only [cross_apply, Fin.mk_one, Pi.add_apply, Matrix.cons_val_one, Matrix.cons_val_zero]
    convert ((hf' 2).mul (hg' 0)).sub ((hf' 0).mul (hg' 2)) using 1; ring
  · simp only [cross_apply, Fin.reduceFinMk, Pi.add_apply, Matrix.cons_val_two,
      Matrix.tail_cons, Matrix.head_cons]
    convert ((hf' 0).mul (hg' 1)).sub ((hf' 1).mul (hg' 0)) using 1; ring

/-- Derivative of the Euclidean norm along a curve that does not pass through `0`. -/
theorem hasDerivAt_vnorm {d : ℝ → V3} {d' : V3} (hd : HasDerivAt d d' t) (h0 : d t ≠ 0) :
    HasDerivAt (fun s => vnorm (d s)) ((vnorm (d t))⁻¹ * (d t ⬝ᵥ d')) t := by
  have h2 := hasDerivAt_dot hd hd
  have hpos : 0 < vnorm (d t) := vnorm_pos_iff.mpr h0
  have hne : d t ⬝ᵥ d t ≠ 0 := by rw [← vnorm_mul_self]; positivity
  have := h2.sqrt hne
  refine this.congr_deriv ?_
  rw [dotProduct_comm d' (d t), ← vnorm]
  field_simp
  ring

/-- Derivative of the normalized vector `d/|d|` along a curve that does not pass through `0`:
`(d/|d|)' = P_{d/|d|} d' / |d|`. -/
theorem hasDerivAt_unit {d : ℝ → V3} {d' : V3} (hd : HasDerivAt d d' t) (h0 : d t ≠ 0) :
    HasDerivAt (fun s => (vnorm (d s))⁻¹ • d s)
      ((vnorm (d t))⁻¹ • proj ((vnorm (d t))⁻¹ • d t) d') t := by
  have hpos : 0 < vnorm (d t) := vnorm_pos_iff.mpr h0
  have hn := (hasDerivAt_vnorm hd h0).inv hpos.ne'
  have := hn.smul hd
  refine this.congr_deriv ?_
  ext i
  simp only [proj, Pi.smul_apply, Pi.sub_apply, smul_eq_mul, Pi.add_apply, smul_dotProduct,
    Pi.inv_apply]
  field_simp
  ring

end deriv

section curveDeriv

variable {V E : Type*} {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ} {ep : E → Fin 2 → V} {I : E}

theorem hasDerivAt_edgeVec (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) :
    HasDerivAt (fun s => edgeVec (Γ s) ep I) (Γ' (ep I 1) - Γ' (ep I 0)) t :=
  (hΓ (ep I 1)).sub (hΓ (ep I 0))

/-- **Derivative of an edge length**: `ℓ_I' = Σ_a s_a T_I · Γ'_{i_a} = T_I · (Γ'_{i₂} − Γ'_{i₁})`. -/
theorem hasDerivAt_len (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (h : 0 < len (Γ t) ep I) :
    HasDerivAt (fun s => len (Γ s) ep I) (∑ a, sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) t := by
  have hd : edgeVec (Γ t) ep I ≠ 0 := vnorm_pos_iff.mp h
  have := hasDerivAt_vnorm (hasDerivAt_edgeVec (ep := ep) (I := I) hΓ) hd
  refine this.congr_deriv ?_
  simp only [Fin.sum_univ_two, sgn_zero, sgn_one, tangent, len, smul_dotProduct, smul_eq_mul,
    dotProduct_sub]
  ring

/-- **Derivative of a unit tangent**: `T_I' = ℓ_I⁻¹ Σ_a s_a P_{T_I} Γ'_{i_a}`. -/
theorem hasDerivAt_tangent (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (h : 0 < len (Γ t) ep I) :
    HasDerivAt (fun s => tangent (Γ s) ep I)
      ((len (Γ t) ep I)⁻¹ • ∑ a, sgn a • proj (tangent (Γ t) ep I) (Γ' (ep I a))) t := by
  have hd : edgeVec (Γ t) ep I ≠ 0 := vnorm_pos_iff.mp h
  have := hasDerivAt_unit (hasDerivAt_edgeVec (ep := ep) (I := I) hΓ) hd
  refine this.congr_deriv ?_
  ext i
  simp only [Fin.sum_univ_two, sgn_zero, sgn_one, proj, tangent, len, Pi.smul_apply,
    Pi.add_apply, Pi.sub_apply, smul_eq_mul, dotProduct_sub, neg_mul, one_mul]
  ring

/-- **Derivative of a midpoint**: `x_I' = (Γ'_{i₁} + Γ'_{i₂})/2`. -/
theorem hasDerivAt_mid (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) :
    HasDerivAt (fun s => mid (Γ s) ep I) ((1 / 2 : ℝ) • (Γ' (ep I 0) + Γ' (ep I 1))) t :=
  ((hΓ (ep I 0)).add (hΓ (ep I 1))).const_smul (1 / 2 : ℝ)

/-- The derivative of a unit tangent is orthogonal to it. -/
theorem tangent_dot_deriv (h : 0 < len (Γ t) ep I) :
    tangent (Γ t) ep I ⬝ᵥ
      ((len (Γ t) ep I)⁻¹ • ∑ a, sgn a • proj (tangent (Γ t) ep I) (Γ' (ep I a))) = 0 := by
  rw [dotProduct_smul, dotProduct_sum]
  simp only [dotProduct_smul, dot_proj_self (tangent_dot_self h), smul_zero, Finset.sum_const_zero]

end curveDeriv


end DDG
