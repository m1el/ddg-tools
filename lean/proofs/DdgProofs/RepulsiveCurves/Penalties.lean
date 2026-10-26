import DdgProofs.RepulsiveCurves.Constraints

/-!
# Repulsive Curves: the penalties (paper Sec. 8.1; SPEC §5.3, D19, D20)

Values and their derivatives along a differentiable motion `Γ` of the vertices (`Basic.lean`); the
gradient blocks are read off as in `Constraints.lean`. `ℓ_I' = Σ_a s_a T_I·Γ'_{i_a}`.

| penalty | value | gradient |
|---|---|---|
| `TotalLength` | `w Σ ℓ_I` | `i₂ += w T_I`, `i₁ −= w T_I` |
| `LengthDifference` | `w Σ_v (ℓ_{I_v} − ℓ_{J_v})²` | `2w(ℓ_I − ℓ_J)(∇ℓ_I − ∇ℓ_J)` |
| `FieldAlignment` | `w Σ ℓ_I |T_I × X(x_I)|²` | `s_a|c|²T + 2s_a P_T(X × c) + ℓ J_Xᵀ(c × T)` |
| `Hencky` | `w Σ ℓ⁰_I ln²(ℓ_I/ℓ⁰_I)` | `∂/∂ℓ_I = 2w(ℓ⁰_I/ℓ_I) ln(ℓ_I/ℓ⁰_I)` |
| `SurfacePotential` | `w Σ_v m_v Σ_t A_t δ_vt^{−p}` | chain rule with `∇δ` (cited) and `∂m_v` |
-/

noncomputable section

namespace DdgProofs

namespace RC

open Matrix Real Finset

variable {V E : Type*} [Fintype E] {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
  {ep : E → Fin 2 → V}

/-- **Total-length penalty** `w Σ_I ℓ_I`: gradient `+w T_I` at `i₂`, `−w T_I` at `i₁`.

[formula: rc-penalty-total-length] -/
theorem penalty_total_length (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (h : ∀ I, 0 < len (Γ t) ep I) (w : ℝ) :
    HasDerivAt (fun s => w * ∑ I, len (Γ s) ep I)
      (∑ I, ∑ a, (w * sgn a) • tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a)) t := by
  have := (HasDerivAt.fun_sum (u := Finset.univ) fun I _ =>
    hasDerivAt_len (ep := ep) hΓ (h I)).const_mul w
  refine this.congr_deriv ?_
  simp only [Finset.mul_sum, smul_dotProduct, smul_eq_mul]
  exact Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun a _ => by ring

omit [Fintype E] in
/-- **Length-difference penalty** `w Σ_v (ℓ_{I_v} − ℓ_{J_v})²` over the degree-2 vertices `v` with
incident edges `I_v, J_v` (any index type): derivative `w Σ_v 2(ℓ_I − ℓ_J)(ℓ_I' − ℓ_J')`.

[formula: rc-penalty-length-difference] -/
theorem penalty_length_difference (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (h : ∀ I, 0 < len (Γ t) ep I) {W : Type*} [Fintype W] (Iv Jv : W → E) (w : ℝ) :
    HasDerivAt (fun s => w * ∑ v, (len (Γ s) ep (Iv v) - len (Γ s) ep (Jv v)) ^ 2)
      (w * ∑ v, 2 * (len (Γ t) ep (Iv v) - len (Γ t) ep (Jv v)) *
        (∑ a, sgn a * (tangent (Γ t) ep (Iv v) ⬝ᵥ Γ' (ep (Iv v) a))
          - ∑ a, sgn a * (tangent (Γ t) ep (Jv v) ⬝ᵥ Γ' (ep (Jv v) a)))) t := by
  have := (HasDerivAt.fun_sum (u := Finset.univ) fun v _ =>
    ((hasDerivAt_len (ep := ep) hΓ (h (Iv v))).sub
      (hasDerivAt_len (ep := ep) hΓ (h (Jv v)))).pow 2).const_mul w
  refine this.congr_deriv ?_
  congr 1
  refine Finset.sum_congr rfl fun v _ => ?_
  simp only [Pi.sub_apply]
  push_cast
  ring

/-- The Hencky term of one edge as a function of its length: `w ℓ⁰ ln²(ℓ/ℓ⁰)` has derivative
`2w (ℓ⁰/ℓ) ln(ℓ/ℓ⁰)` at `ℓ > 0` (`ℓ⁰ > 0`).

[formula: rc-penalty-hencky] -/
theorem hencky_hasDerivAt {ℓ0 ℓ : ℝ} (h0 : 0 < ℓ0) (hl : 0 < ℓ) (w : ℝ) :
    HasDerivAt (fun x => w * ℓ0 * Real.log (x / ℓ0) ^ 2)
      (2 * w * (ℓ0 / ℓ) * Real.log (ℓ / ℓ0)) ℓ := by
  have hd : HasDerivAt (fun x => x / ℓ0) (1 / ℓ0) ℓ := (hasDerivAt_id ℓ).div_const ℓ0
  have := ((hd.log (by positivity)).pow 2).const_mul (w * ℓ0)
  refine this.congr_deriv ?_
  field_simp
  push_cast
  ring

/-- **Hencky penalty** `w Σ_I ℓ⁰_I ln²(ℓ_I/ℓ⁰_I)` (Supp. 2.2): derivative
`Σ_I 2w(ℓ⁰_I/ℓ_I) ln(ℓ_I/ℓ⁰_I) ℓ_I'`.

[formula: rc-penalty-hencky] -/
theorem penalty_hencky (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (h : ∀ I, 0 < len (Γ t) ep I) {ℓ0 : E → ℝ} (h0 : ∀ I, 0 < ℓ0 I) (w : ℝ) :
    HasDerivAt (fun s => ∑ I, w * ℓ0 I * Real.log (len (Γ s) ep I / ℓ0 I) ^ 2)
      (∑ I, 2 * w * (ℓ0 I / len (Γ t) ep I) * Real.log (len (Γ t) ep I / ℓ0 I) *
        ∑ a, sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) t := by
  refine HasDerivAt.fun_sum (u := Finset.univ) fun I _ => ?_
  have := (hencky_hasDerivAt (h0 I) (h I) w).comp t (hasDerivAt_len (ep := ep) hΓ (h I))
  exact this

omit [Fintype E] in
/-- **Field-alignment penalty of one edge** `ℓ_I |c|²`, `c = T_I × X(x_I)` (`|c|² = c·c`), for a
field `X` whose derivative along the midpoint is `J (x_I)'` (`J` its Jacobian at `x_I`):
the gradient block at `i_a` is `s_a |c|² T_I + 2 s_a P_{T_I}(X × c) + ℓ_I Jᵀ(c × T_I)`.

[formula: rc-penalty-field] -/
theorem penalty_field_edge (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {I : E}
    (hI : 0 < len (Γ t) ep I) {X : V3 → V3} {J : Matrix (Fin 3) (Fin 3) ℝ}
    (hX : HasDerivAt (fun s => X (mid (Γ s) ep I))
      (J *ᵥ ((1 / 2 : ℝ) • (Γ' (ep I 0) + Γ' (ep I 1)))) t) :
    HasDerivAt (fun s => len (Γ s) ep I *
        (crossProduct (tangent (Γ s) ep I) (X (mid (Γ s) ep I)) ⬝ᵥ
          crossProduct (tangent (Γ s) ep I) (X (mid (Γ s) ep I))))
      (∑ a, ((sgn a * (crossProduct (tangent (Γ t) ep I) (X (mid (Γ t) ep I)) ⬝ᵥ
            crossProduct (tangent (Γ t) ep I) (X (mid (Γ t) ep I)))) • tangent (Γ t) ep I
          + (2 * sgn a) • proj (tangent (Γ t) ep I)
              (crossProduct (X (mid (Γ t) ep I))
                (crossProduct (tangent (Γ t) ep I) (X (mid (Γ t) ep I))))
          + len (Γ t) ep I • (Jᵀ *ᵥ crossProduct
              (crossProduct (tangent (Γ t) ep I) (X (mid (Γ t) ep I))) (tangent (Γ t) ep I)))
        ⬝ᵥ Γ' (ep I a)) t := by
  have hT := hasDerivAt_tangent (ep := ep) (I := I) hΓ hI
  have hc := hasDerivAt_cross hT hX
  have hcc := hasDerivAt_dot hc hc
  have hl := hasDerivAt_len (ep := ep) (I := I) hΓ hI
  refine (hl.mul hcc).congr_deriv ?_
  set T := tangent (Γ t) ep I
  set x := X (mid (Γ t) ep I)
  set c := crossProduct T x
  set ℓ := len (Γ t) ep I
  set cc := c ⬝ᵥ c with hccdef
  have hl0 : ℓ ≠ 0 := hI.ne'
  -- the two triple products
  have k1 : ∀ T' : V3, crossProduct T' x ⬝ᵥ c = T' ⬝ᵥ crossProduct x c := fun T' => by
    rw [dotProduct_comm, triple_product_permutation]
  have k2 : ∀ X' : V3, crossProduct T X' ⬝ᵥ c = X' ⬝ᵥ crossProduct c T := fun X' => by
    rw [dotProduct_comm, triple_product_permutation, triple_product_permutation]
  have k3 : ∀ y : V3, (J *ᵥ ((1 / 2 : ℝ) • (Γ' (ep I 0) + Γ' (ep I 1)))) ⬝ᵥ y =
      (1 / 2 : ℝ) * ((Γ' (ep I 0) ⬝ᵥ (Jᵀ *ᵥ y)) + (Γ' (ep I 1) ⬝ᵥ (Jᵀ *ᵥ y))) := fun y => by
    rw [dotProduct_comm, dotProduct_mulVec, ← mulVec_transpose, dotProduct_comm, smul_dotProduct,
      add_dotProduct, smul_eq_mul]
  have k4 : ∀ y : V3, ((ℓ)⁻¹ • ∑ a, sgn a • proj T (Γ' (ep I a))) ⬝ᵥ y =
      ℓ⁻¹ * ∑ a, sgn a * (Γ' (ep I a) ⬝ᵥ proj T y) := fun y => by
    rw [smul_dotProduct, sum_dotProduct, smul_eq_mul]
    simp only [smul_dotProduct, smul_eq_mul, proj_dot_comm]
  rw [dotProduct_comm c (crossProduct (ℓ⁻¹ • ∑ a, sgn a • proj T (Γ' (ep I a))) x
    + crossProduct T (J *ᵥ ((1 / 2 : ℝ) • (Γ' (ep I 0) + Γ' (ep I 1)))))]
  rw [add_dotProduct, k1 (ℓ⁻¹ • ∑ a, sgn a • proj T (Γ' (ep I a))),
    k2 (J *ᵥ ((1 / 2 : ℝ) • (Γ' (ep I 0) + Γ' (ep I 1)))), k3, k4]
  simp only [Fin.sum_univ_two, add_dotProduct, smul_dotProduct, smul_eq_mul, sgn_zero, sgn_one,
    dotProduct_comm T, dotProduct_comm (proj T (crossProduct x c)),
    dotProduct_comm (Jᵀ *ᵥ crossProduct c T)]
  field_simp
  ring

/-! ## Surface potential (D19) -/

omit [Fintype E] in
/-- The far-field (lumped) term `A_N |γ_v − c_N|^{−p}`: along a motion, `|g − c|^{−p}` has
derivative `−p |d|^{−p−2} (d·g')` with `d = g − c ≠ 0`.

[formula: rc-penalty-surface] -/
theorem hasDerivAt_dist_rpow {g : ℝ → V3} {g' : V3} (hg : HasDerivAt g g' t) (cN : V3)
    (hd : g t - cN ≠ 0) (p : ℝ) :
    HasDerivAt (fun s => vnorm (g s - cN) ^ (-p))
      (-p * vnorm (g t - cN) ^ (-p - 2) * ((g t - cN) ⬝ᵥ g')) t := by
  have hρ : 0 < vnorm (g t - cN) := vnorm_pos_iff.mpr hd
  have := (hasDerivAt_vnorm (hg.sub_const cN) hd).rpow_const (p := -p) (Or.inl hρ.ne')
  refine this.congr_deriv ?_
  have e : vnorm (g t - cN) ^ (-p - 1) = vnorm (g t - cN) ^ (-p - 2) * vnorm (g t - cN) := by
    rw [show -p - 1 = (-p - 2) + 1 by ring, Real.rpow_add_one hρ.ne']
  rw [e]
  field_simp

/-- **Derivative of the dual lengths** (`m_v = ½ Σ_{I∋v} ℓ_I`, the `∂m_v/∂γ` terms of the surface
potential): `m_v' = ½ Σ_{I∋v} ℓ_I'`.

[formula: rc-penalty-surface] -/
theorem hasDerivAt_dualLen [DecidableEq V] (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (h : ∀ I, 0 < len (Γ t) ep I) (v : V) :
    HasDerivAt (fun s => dualLen (Γ s) ep v)
      ((1 / 2) * ∑ I, ∑ a : Fin 2, if ep I a = v then
        ∑ b, sgn b * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I b)) else 0) t := by
  unfold dualLen
  refine HasDerivAt.const_mul _ (HasDerivAt.fun_sum fun I _ => HasDerivAt.fun_sum fun a _ => ?_)
  by_cases hv : ep I a = v
  · simp only [hv, ite_true]
    exact hasDerivAt_len hΓ (h I)
  · simp only [hv, ite_false]
    exact hasDerivAt_const _ _

/-- **The surface potential by the chain rule** (D19): with dual lengths `m_v(s)` and
point–triangle distances `δ_vk(s) > 0` differentiable along the motion (for the distance,
`δ' = ∇δ·Γ'_v` with `∇δ = (γ_v − π(γ_v))/δ`, the cited `rc-distance-gradient`), the derivative of
`w Σ_v m_v Σ_k A_k δ_vk^{−p}` is `w Σ_v (m_v' Σ_k A_k δ_vk^{−p} + m_v Σ_k A_k (−p) δ_vk^{−p−1} δ_vk')`.

[formula: rc-penalty-surface] -/
theorem penalty_surface {W K : Type*} [Fintype W] [Fintype K] {m : W → ℝ → ℝ} {m' : W → ℝ}
    {δ : W → K → ℝ → ℝ} {δ' : W → K → ℝ} (hm : ∀ v, HasDerivAt (m v) (m' v) t)
    (hδ : ∀ v k, HasDerivAt (δ v k) (δ' v k) t) (hpos : ∀ v k, 0 < δ v k t) (A : K → ℝ)
    (w p : ℝ) :
    HasDerivAt (fun s => w * ∑ v, m v s * ∑ k, A k * δ v k s ^ (-p))
      (w * ∑ v, (m' v * ∑ k, A k * δ v k t ^ (-p)
        + m v t * ∑ k, A k * (-p * δ v k t ^ (-p - 1) * δ' v k))) t := by
  refine HasDerivAt.const_mul w (HasDerivAt.fun_sum fun v _ => ?_)
  have hs : HasDerivAt (fun s => ∑ k, A k * δ v k s ^ (-p))
      (∑ k, A k * (-p * δ v k t ^ (-p - 1) * δ' v k)) t := by
    refine HasDerivAt.fun_sum fun k _ => ?_
    have := ((hδ v k).rpow_const (p := -p) (Or.inl (hpos v k).ne')).const_mul (A k)
    refine this.congr_deriv ?_
    ring
  exact ((hm v).mul hs).congr_deriv (by ring)

end RC

end DdgProofs
