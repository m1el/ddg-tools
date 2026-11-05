import DDG.Repulsive.Surfaces.Basic

/-!
# Repulsive Surfaces: derivatives of the face quantities (SPEC §3.2)

Derivatives are taken along an arbitrary differentiable motion of the vertices `s ↦ Γ s`
(`HasDerivAt` of every vertex position, velocities `Γ'`), as in the curves package: a statement
`d/ds q(Γ s) = Σ_a g_a · Γ'_{fv S a}` contains the gradient `g` and the chain rule along every
differentiable motion. For a nondegenerate face `S` with corners `a ∈ Fin 3` and opposite edges
`e_a = f_{a+2} − f_{a+1}`:

```
n_S'  = Σ_a e_a × Γ'_a
a_S'  = Σ_a (½ N_S × e_a) · Γ'_a                              (rs-face-area-gradient)
N_S'  = Σ_a (1/(2a_S)) (e_a × N_S) (N_S · Γ'_a)              (rs-face-normal-variation)
X_S'  = ⅓ Σ_a Γ'_a
```

The key step of the normal variation is `P_N(e × w) = (N·w)(e × N)` for a unit `N ⊥ e`
(`proj_cross`): the tangential part of a velocity only rotates `n` along itself.
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Real Finset

/-- `P_N(e × w) = (N·w)(e × N)` for a unit `N` with `N ⊥ e`. -/
theorem proj_cross {N e : V3} (hN : N ⬝ᵥ N = 1) (he : N ⬝ᵥ e = 0) (w : V3) :
    proj N (e ⨯₃ w) = (N ⬝ᵥ w) • (e ⨯₃ N) := by
  -- `N × ((e × w) × N) = (e × w) − (N·(e × w)) N` and `(e × w) × N = (e·N) w − (w·N) e`
  have key : ∀ i, (proj N (e ⨯₃ w)) i - ((N ⬝ᵥ w) • (e ⨯₃ N)) i
      = (1 - N ⬝ᵥ N) * (e ⨯₃ w) i + (N ⬝ᵥ e) * (N ⨯₃ w) i := by
    intro i
    fin_cases i <;> simp [proj, cross_apply, dotProduct, Fin.sum_univ_three] <;> ring
  ext i
  have := key i
  rw [hN, he] at this
  simp only [sub_self, zero_mul, add_zero] at this
  linarith

section motion

variable {V F : Type*} {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ} {fv : F → Fin 3 → V} {S : F}

theorem hasDerivAt_corner (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (fv : F → Fin 3 → V)
    (S : F) (a : Fin 3) :
    HasDerivAt (fun s => corner (Γ s) fv S a) (Γ' (fv S a)) t := hΓ (fv S a)

/-- `n_S' = Σ_a e_a × Γ'_a`. -/
theorem hasDerivAt_nvec (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) :
    HasDerivAt (fun s => nvec (Γ s) fv S) (∑ a, opp (Γ t) fv S a ⨯₃ Γ' (fv S a)) t := by
  have h1 := (hasDerivAt_corner hΓ fv S 1).sub (hasDerivAt_corner hΓ fv S 0)
  have h2 := (hasDerivAt_corner hΓ fv S 2).sub (hasDerivAt_corner hΓ fv S 0)
  have := hasDerivAt_cross h1 h2
  refine this.congr_deriv ?_
  simp only [Fin.sum_univ_three, opp, fin3_add_two_zero, fin3_add_one_zero, fin3_add_two_one,
    fin3_add_one_one, fin3_add_two_two, fin3_add_one_two]
  apply v3_ext <;> simp [cross_apply, corner] <;> ring

/-- `X_S' = ⅓ Σ_a Γ'_a`. -/
theorem hasDerivAt_bary (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) :
    HasDerivAt (fun s => bary (Γ s) fv S) ((1 / 3 : ℝ) • ∑ a, Γ' (fv S a)) t := by
  have := (((hasDerivAt_corner hΓ fv S 0).add (hasDerivAt_corner hΓ fv S 1)).add
    (hasDerivAt_corner hΓ fv S 2)).const_smul (1 / 3 : ℝ)
  refine this.congr_deriv ?_
  simp [Fin.sum_univ_three]

/-- The polynomial identity behind the area gradient: `n · n' = Σ_a (n × e_a) · w_a`. -/
theorem nvec_dot_deriv (f : V → V3) (w : Fin 3 → V3) :
    nvec f fv S ⬝ᵥ (∑ a, opp f fv S a ⨯₃ w a) = ∑ a, (nvec f fv S ⨯₃ opp f fv S a) ⬝ᵥ w a := by
  simp only [Fin.sum_univ_three, dotProduct_add]
  simp only [nvec, opp, corner, fin3_add_two_zero, fin3_add_one_zero, fin3_add_two_one,
    fin3_add_one_one, fin3_add_two_two, fin3_add_one_two, cross_apply, dotProduct,
    Fin.sum_univ_three]
  simp
  ring

/-- **The face-area gradient** (SPEC §3.2): for a nondegenerate face,
`d/ds a_S = Σ_a (½ N_S × e_a) · Γ'_a`, i.e. `∂a_S/∂f_i = ½ N_S × e_jk` (cyclically).

[formula: rs-face-area-gradient] -/
theorem hasDerivAt_area (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (hS : Nondeg (Γ t) fv S) :
    HasDerivAt (fun s => area (Γ s) fv S)
      (∑ a, ((1 / 2 : ℝ) • (normal (Γ t) fv S ⨯₃ opp (Γ t) fv S a)) ⬝ᵥ Γ' (fv S a)) t := by
  have h := (hasDerivAt_vnorm (hasDerivAt_nvec (fv := fv) (S := S) hΓ) hS).div_const 2
  refine h.congr_deriv ?_
  rw [nvec_dot_deriv]
  simp only [normal, map_smul, LinearMap.smul_apply, smul_dotProduct, smul_eq_mul,
    ← Finset.mul_sum]
  ring

/-- **The face-normal variation** (SPEC §3.2): for a nondegenerate face,
`d/ds N_S = Σ_a (1/(2a_S)) (e_a × N_S)(N_S · Γ'_a)`, i.e. `∂N_S/∂f_i = (1/(2a_S)) (e_jk × N_S) N_Sᵀ`
(cyclically). Only the normal component of a vertex velocity turns the normal.

[formula: rs-face-normal-variation] -/
theorem hasDerivAt_normal (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (hS : Nondeg (Γ t) fv S) :
    HasDerivAt (fun s => normal (Γ s) fv S)
      (∑ a, ((normal (Γ t) fv S ⬝ᵥ Γ' (fv S a)) / (2 * area (Γ t) fv S)) •
        (opp (Γ t) fv S a ⨯₃ normal (Γ t) fv S)) t := by
  have h := hasDerivAt_unit (hasDerivAt_nvec (fv := fv) (S := S) hΓ) hS
  refine h.congr_deriv ?_
  obtain ⟨hunit, -, -, horth, -⟩ := face_quantities hS
  have hpos : 0 < vnorm (nvec (Γ t) fv S) := vnorm_pos_iff.mpr hS
  -- `P_N` is additive
  have hlin : ∀ (x y : V3), proj (normal (Γ t) fv S) (x + y) =
      proj (normal (Γ t) fv S) x + proj (normal (Γ t) fv S) y := by
    intro x y; simp only [proj, dotProduct_add, add_smul]; abel
  change (vnorm (nvec (Γ t) fv S))⁻¹ • proj (normal (Γ t) fv S)
      (∑ a, opp (Γ t) fv S a ⨯₃ Γ' (fv S a)) = _
  simp only [Fin.sum_univ_three, hlin, proj_cross hunit (horth _), smul_add, smul_smul]
  have h2 : 2 * area (Γ t) fv S = vnorm (nvec (Γ t) fv S) := by rw [area]; ring
  rw [h2]
  congr 1; congr 1
  all_goals (congr 1; field_simp)

end motion

end RS

end DDG
