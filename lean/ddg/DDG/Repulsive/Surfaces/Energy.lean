import DDG.Repulsive.Surfaces.Kernel
import DDG.Repulsive.Curves.Energy

/-!
# Repulsive Surfaces: the discrete energy, its invariances and its gradient (paper Eq. 9, 10;
SPEC §3.1–§3.3)

```
Ê(f) = Σ_{S} Σ_{T ≠ S} K(X_S, N_S; X_T) a_S a_T          (ordered pairs, S = T excluded)
```

* `rs-energy`: the definition (no `2^p`, no ½, neighbours not excluded), nonnegative, invariant
  under flipping the orientation of any set of faces.
* `rs-energy-scaling`: `Ê(cf) = c^{4+α−β} Ê(f)` for `c > 0`; invariance under translations and
  orthogonal maps (rotations **and** reflections: an improper map flips every normal, which the
  kernel ignores).
* `rs-energy-gradient`: the chain rule of SPEC §3.2. Along any differentiable motion `Γ` with
  every face nondegenerate and all barycenters distinct at `Γ t`, for `α > 1`,
  `d/ds Ê(Γ s) = Σ_v g_v · Γ'_v` with `g = energyGrad`, assembled from the per-pair corner
  gradients `pairGradS` (through `X_S`, `N_S`, `a_S`) and `pairGradT` (through `X_T`, `a_T`).
  No condition on `ν = N_S·(X_S − X_T)` (coplanar neighbours are fine, `rs-kernel-gradient`).
* `rs-gradient-identities`: `Σ_v g_v = 0`, `Σ_v f_v × g_v = 0`, `Σ_v f_v·g_v = (4 + α − β) Ê`
  (through the curves package's `RC.gradient_identities_of_curves`).
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Real Finset Filter Topology
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff
  dot_self_nonneg unit_dot_self proj hasDerivAt_dot hasDerivAt_cross hasDerivAt_vnorm
  hasDerivAt_unit orth_dot orth_vnorm vnorm_smul vnorm_neg)

variable {V F : Type*} [Fintype V] [Fintype F] [DecidableEq V] [DecidableEq F]

/-- One ordered pair term `K(X_S, N_S; X_T) a_S a_T`. -/
def pairEnergy (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (S T : F) : ℝ :=
  kernel α β (bary f fv S) (normal f fv S) (bary f fv T) * area f fv S * area f fv T

/-- **The discrete tangent-point energy** (paper Eq. 10, SPEC §3.1): the sum over ordered pairs
of distinct faces, both orientations, neighbours included.

[formula: rs-energy] -/
def energy (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) : ℝ :=
  ∑ S, ∑ T, if S ≠ T then pairEnergy α β f fv S T else 0

omit [Fintype V] [DecidableEq V] in
/-- The energy is nonnegative and each ordered pair enters once.

[formula: rs-energy] -/
theorem energy_nonneg (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) : 0 ≤ energy α β f fv := by
  refine Finset.sum_nonneg fun S _ => Finset.sum_nonneg fun T _ => ?_
  split_ifs
  · exact mul_nonneg (mul_nonneg (kernel_nonneg _ _ _ _ _) (area_nonneg _ _ _)) (area_nonneg _ _ _)
  · exact le_rfl

/-! ## Invariances and scaling -/

section invariance

omit [Fintype V] [Fintype F] [DecidableEq V] [DecidableEq F]

variable (f : V → V3) (fv : F → Fin 3 → V)

theorem nvec_smul (c : ℝ) (S : F) : nvec (c • f) fv S = c ^ 2 • nvec f fv S := by
  apply v3_ext <;> simp [nvec, corner, cross_apply] <;> ring

theorem area_smul {c : ℝ} (S : F) : area (c • f) fv S = c ^ 2 * area f fv S := by
  rw [area, area, nvec_smul, vnorm_smul, abs_of_nonneg (sq_nonneg c)]; ring

theorem normal_smul {c : ℝ} (hc : c ≠ 0) (S : F) : normal (c • f) fv S = normal f fv S := by
  rw [normal, normal, nvec_smul, vnorm_smul, abs_of_nonneg (sq_nonneg c), smul_smul]
  congr 1
  field_simp

theorem bary_smul (c : ℝ) (S : F) : bary (c • f) fv S = c • bary f fv S := by
  simp only [bary, corner, Pi.smul_apply, smul_add, smul_comm c]

theorem nvec_translate (a : V3) (S : F) : nvec (fun v => f v + a) fv S = nvec f fv S := by
  simp [nvec, corner]

theorem bary_translate (a : V3) (S : F) :
    bary (fun v => f v + a) fv S = bary f fv S + a := by
  apply v3_ext <;> simp [bary, corner] <;> ring

/-- `(Ra) × (Rb) = adj(R)ᵀ (a × b)` for any `3 × 3` matrix. -/
theorem cross_mulVec (R : Matrix (Fin 3) (Fin 3) ℝ) (a b : V3) :
    (R *ᵥ a) ⨯₃ (R *ᵥ b) = (R.adjugate)ᵀ *ᵥ (a ⨯₃ b) := by
  apply v3_ext <;>
    simp [cross_apply, mulVec, dotProduct, Fin.sum_univ_three, adjugate_fin_three] <;> ring

/-- For an orthogonal `R`: `(Ra) × (Rb) = det R · R (a × b)` and `det R = ±1`. -/
theorem cross_orth {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : Rᵀ * R = 1) (a b : V3) :
    (R *ᵥ a) ⨯₃ (R *ᵥ b) = R.det • (R *ᵥ (a ⨯₃ b)) ∧ (R.det = 1 ∨ R.det = -1) := by
  have hR' : R * Rᵀ = 1 := mul_eq_one_comm.mp hR
  have hadj : R.adjugate = R.det • Rᵀ := by
    calc R.adjugate = R.adjugate * (R * Rᵀ) := by rw [hR', mul_one]
      _ = (R.adjugate * R) * Rᵀ := by rw [mul_assoc]
      _ = R.det • Rᵀ := by rw [adjugate_mul, smul_mul_assoc, one_mul]
  refine ⟨?_, ?_⟩
  · rw [cross_mulVec, hadj, transpose_smul, transpose_transpose, smul_mulVec]
  · have h := congrArg det hR
    rw [det_mul, det_transpose, det_one] at h
    have : (R.det - 1) * (R.det + 1) = 0 := by ring_nf; linarith
    rcases mul_eq_zero.mp this with h1 | h1
    · left; linarith
    · right; linarith

variable {R : Matrix (Fin 3) (Fin 3) ℝ}

theorem nvec_orth (hR : Rᵀ * R = 1) (S : F) :
    nvec (fun v => R *ᵥ f v) fv S = R.det • (R *ᵥ nvec f fv S) := by
  simp only [nvec, corner, ← mulVec_sub, (cross_orth hR _ _).1]

theorem vnorm_det_smul (hR : Rᵀ * R = 1) (x : V3) : vnorm (R.det • (R *ᵥ x)) = vnorm x := by
  rw [vnorm_smul, orth_vnorm hR]
  rcases (cross_orth hR 0 0).2 with h | h <;> simp [h]

theorem area_orth (hR : Rᵀ * R = 1) (S : F) : area (fun v => R *ᵥ f v) fv S = area f fv S := by
  rw [area, area, nvec_orth f fv hR, vnorm_det_smul hR]

theorem normal_orth (hR : Rᵀ * R = 1) (S : F) :
    normal (fun v => R *ᵥ f v) fv S = R.det • (R *ᵥ normal f fv S) := by
  rw [normal, normal, nvec_orth f fv hR, vnorm_det_smul hR, mulVec_smul, smul_comm]

theorem bary_orth (S : F) : bary (fun v => R *ᵥ f v) fv S = R *ᵥ bary f fv S := by
  simp only [bary, corner, mulVec_smul, mulVec_add]

/-- Flipping the orientation of the faces with `o S = true` (corners `1` and `2` swapped). -/
def flipFv (fv : F → Fin 3 → V) (o : F → Bool) (S : F) : Fin 3 → V :=
  if o S then ![fv S 0, fv S 2, fv S 1] else fv S

theorem flip_quantities (o : F → Bool) (S : F) :
    area f (flipFv fv o) S = area f fv S ∧ bary (f) (flipFv fv o) S = bary f fv S ∧
      (normal f (flipFv fv o) S = normal f fv S ∨ normal f (flipFv fv o) S = -normal f fv S) := by
  by_cases h : o S
  · have hn : nvec f (flipFv fv o) S = -nvec f fv S := by
      apply v3_ext <;> simp [nvec, corner, flipFv, h, cross_apply] <;> ring
    refine ⟨?_, ?_, Or.inr ?_⟩
    · rw [area, area, hn, vnorm_neg]
    · apply v3_ext <;> simp [bary, corner, flipFv, h] <;> ring
    · rw [normal, normal, hn, vnorm_neg, smul_neg]
  · simp [area, normal, nvec, bary, corner, flipFv, h]

end invariance

omit [Fintype V] [DecidableEq V] in
/-- **Homogeneity of the energy** (SPEC §3.3): `Ê(cf) = c^{4+α−β} Ê(f)` for `c > 0`.

[formula: rs-energy-scaling] -/
theorem energy_smul {c : ℝ} (hc : 0 < c) (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) :
    energy α β (c • f) fv = c ^ (4 + α - β) * energy α β f fv := by
  have hpow : c ^ (4 + α - β) = c ^ (α - β) * c ^ 2 * c ^ 2 := by
    rw [show 4 + α - β = (α - β) + ((4 : ℕ) : ℝ) by push_cast; ring, Real.rpow_add hc,
      Real.rpow_natCast]
    ring
  simp only [energy, pairEnergy, Finset.mul_sum, area_smul, normal_smul f fv hc.ne', bary_smul,
    kernel_smul hc]
  refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
  split_ifs
  · rw [hpow]; ring
  · simp

omit [Fintype V] [DecidableEq V] in
/-- **Translation invariance of the energy.**

[formula: rs-energy-scaling] -/
theorem energy_translate (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (a : V3) :
    energy α β (fun v => f v + a) fv = energy α β f fv := by
  simp only [energy, pairEnergy, area, normal, nvec_translate, bary_translate, kernel_translate]

omit [Fintype V] [DecidableEq V] in
/-- **Invariance under orthogonal maps** `Rᵀ R = 1` (rotations and reflections).

[formula: rs-energy-scaling] -/
theorem energy_orth {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : Rᵀ * R = 1) (α β : ℝ) (f : V → V3)
    (fv : F → Fin 3 → V) : energy α β (fun v => R *ᵥ f v) fv = energy α β f fv := by
  simp only [energy, pairEnergy, area_orth f fv hR, normal_orth f fv hR, bary_orth,
    kernel_orth hR α β _ _ _ _ (cross_orth hR 0 0).2]

omit [Fintype V] [DecidableEq V] in
/-- **Invariance under flipping face orientations** (SPEC §2.1: nothing in the energy depends on
the orientation).

[formula: rs-energy] -/
theorem energy_flip (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (o : F → Bool) :
    energy α β f (flipFv fv o) = energy α β f fv := by
  unfold energy pairEnergy
  refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
  obtain ⟨hS1, hS2, hS3⟩ := flip_quantities f fv o S
  obtain ⟨hT1, hT2, -⟩ := flip_quantities f fv o T
  rw [hS1, hS2, hT1, hT2]
  rcases hS3 with h | h <;> rw [h]
  rw [kernel_neg_normal]

/-! ## The differential -/

section differential

variable (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V)

/-- Regular configuration: every face nondegenerate, all barycenters distinct. -/
def Regular : Prop :=
  (∀ S, Nondeg f fv S) ∧ ∀ S T, S ≠ T → bary f fv S ≠ bary f fv T

/-- `∂e_ST/∂f_{S,a}` (corner `a` of the source face `S`): through `X_S` (`⅓ ∂K/∂v`), `N_S`
(`(1/(2a_S)) (∂K/∂N · (e_a × N_S)) N_S`, `rs-face-normal-variation`) and `a_S`
(`K a_T ½ N_S × e_a`, `rs-face-area-gradient`).

[formula: rs-energy-gradient] -/
def pairGradS (S T : F) (a : Fin 3) : V3 :=
  (area f fv S * area f fv T) •
      ((1 / 3 : ℝ) • kgradV α β (normal f fv S) (bary f fv S - bary f fv T)
        + ((kgradN α β (normal f fv S) (bary f fv S - bary f fv T) ⬝ᵥ
            (opp f fv S a ⨯₃ normal f fv S)) / (2 * area f fv S)) • normal f fv S)
    + (kernel α β (bary f fv S) (normal f fv S) (bary f fv T) * area f fv T) •
      ((1 / 2 : ℝ) • (normal f fv S ⨯₃ opp f fv S a))

/-- `∂e_ST/∂f_{T,b}` (corner `b` of the target face `T`): through `X_T` (`−⅓ ∂K/∂v`) and `a_T`.

[formula: rs-energy-gradient] -/
def pairGradT (S T : F) (b : Fin 3) : V3 :=
  -(area f fv S * area f fv T / 3) • kgradV α β (normal f fv S) (bary f fv S - bary f fv T)
    + (kernel α β (bary f fv S) (normal f fv S) (bary f fv T) * area f fv S) •
      ((1 / 2 : ℝ) • (normal f fv T ⨯₃ opp f fv T b))

variable {α β} {fv}

omit [Fintype V] [Fintype F] [DecidableEq V] [DecidableEq F] in
/-- **The differential of one pair term** (SPEC §3.2): for `α > 1`, nondegenerate `S`, `T` and
`X_S ≠ X_T` at `Γ t`, `d/ds e_ST(Γ s) = Σ_a ∂_{S,a}·Γ'_{S,a} + Σ_b ∂_{T,b}·Γ'_{T,b}`.

[formula: rs-energy-gradient] -/
theorem hasDerivAt_pairEnergy (hα : 1 < α) {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {S T : F}
    (hS : Nondeg (Γ t) fv S) (hT : Nondeg (Γ t) fv T) (hST : bary (Γ t) fv S ≠ bary (Γ t) fv T) :
    HasDerivAt (fun s => pairEnergy α β (Γ s) fv S T)
      (∑ a, pairGradS α β (Γ t) fv S T a ⬝ᵥ Γ' (fv S a)
        + ∑ b, pairGradT α β (Γ t) fv S T b ⬝ᵥ Γ' (fv T b)) t := by
  have hK := hasDerivAt_kernel (β := β) hα (hasDerivAt_bary (fv := fv) (S := S) hΓ)
    (hasDerivAt_normal hΓ hS) (hasDerivAt_bary (fv := fv) (S := T) hΓ) hST
  have h := (hK.mul (hasDerivAt_area hΓ hS)).mul (hasDerivAt_area hΓ hT)
  refine h.congr_deriv ?_
  have haS : 0 < area (Γ t) fv S := area_pos hS
  simp only [pairGradS, pairGradT, Fin.sum_univ_three, dotProduct_add, add_dotProduct,
    dotProduct_smul, smul_dotProduct, dotProduct_sub, neg_smul, neg_dotProduct, smul_eq_mul,
    Pi.mul_apply]
  field_simp
  ring

/-- **The vertex gradient** `g_v = Σ_{S ≠ T} (Σ_{a : fv S a = v} ∂_{S,a} e_ST +
Σ_{b : fv T b = v} ∂_{T,b} e_ST)` (SPEC §3.2).

[formula: rs-energy-gradient] -/
def energyGrad (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (v : V) : V3 :=
  ∑ S, ∑ T, if S ≠ T then
    (∑ a, if fv S a = v then pairGradS α β f fv S T a else 0)
      + ∑ b, if fv T b = v then pairGradT α β f fv S T b else 0
  else 0

/-- **The differential of the energy** (SPEC §3.2): for `α > 1`, along a differentiable motion
`Γ` that is regular at `t` (nondegenerate faces, distinct barycenters),
`d/ds Ê(Γ s) = Σ_v g_v · Γ'_v` with `g = energyGrad`.

[formula: rs-energy-gradient] -/
theorem hasDerivAt_energy (hα : 1 < α) {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (hreg : Regular (Γ t) fv) :
    HasDerivAt (fun s => energy α β (Γ s) fv) (∑ v, energyGrad α β (Γ t) fv v ⬝ᵥ Γ' v) t := by
  have hpair : ∀ S T, HasDerivAt
      (fun s => if S ≠ T then pairEnergy α β (Γ s) fv S T else 0)
      (if S ≠ T then
        (∑ a, pairGradS α β (Γ t) fv S T a ⬝ᵥ Γ' (fv S a)
          + ∑ b, pairGradT α β (Γ t) fv S T b ⬝ᵥ Γ' (fv T b)) else 0) t := by
    intro S T
    by_cases hd : S ≠ T
    · simp only [hd, ite_true, ne_eq, not_false_eq_true]
      exact hasDerivAt_pairEnergy hα hΓ (hreg.1 S) (hreg.1 T) (hreg.2 S T hd)
    · simp only [hd, ite_false]
      exact hasDerivAt_const _ _
  have h := HasDerivAt.fun_sum (u := Finset.univ) fun S _ =>
    HasDerivAt.fun_sum (u := Finset.univ) fun T _ => hpair S T
  have hfun : (fun s => energy α β (Γ s) fv) =
      fun s => ∑ S, ∑ T, if S ≠ T then pairEnergy α β (Γ s) fv S T else 0 := rfl
  rw [hfun]
  refine h.congr_deriv ?_
  have ite_dot : ∀ (c : Prop) [Decidable c] (x w : V3),
      (if c then x else 0) ⬝ᵥ w = if c then x ⬝ᵥ w else 0 := by
    intro c _ x w; split_ifs <;> simp
  symm
  unfold energyGrad
  simp only [sum_dotProduct]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun S _ => ?_
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun T _ => ?_
  by_cases hd : S ≠ T
  · simp only [hd, ite_true, ne_eq, not_false_eq_true, add_dotProduct, Finset.sum_add_distrib,
      sum_dotProduct, ite_dot]
    congr 1
    · rw [Finset.sum_comm]
      simp only [Finset.sum_ite_eq, Finset.mem_univ, ite_true]
    · rw [Finset.sum_comm]
      simp only [Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  · simp [hd]

/-- **The gradient identities** (SPEC §3.3): for `α > 1` and a regular configuration, the vertex
gradient satisfies `Σ_v g_v = 0`, `Σ_v f_v × g_v = 0` and `Σ_v f_v · g_v = (4 + α − β) Ê(f)`
(invariance under translations and rotations, homogeneity of degree `4 + α − β`).

[formula: rs-gradient-identities] -/
theorem energy_gradient_identities (hα : 1 < α) {f : V → V3} (hreg : Regular f fv) :
    ∑ v, energyGrad α β f fv v = 0 ∧
      ∑ v, f v ⨯₃ energyGrad α β f fv v = 0 ∧
      ∑ v, f v ⬝ᵥ energyGrad α β f fv v = (4 + α - β) * energy α β f fv :=
  RC.gradient_identities_of_curves (F := fun f => energy α β f fv)
    (fun Γ Γ' h0 hΓ => by
      subst h0
      exact hasDerivAt_energy hα hΓ hreg)
    (fun a => energy_translate α β f fv a)
    (fun _ hR => energy_orth hR α β f fv)
    (fun _ hc => energy_smul hc α β f fv)

end differential

end RS

end DDG
