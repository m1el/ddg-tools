import DdgProofs.RepulsiveCurves.KernelGradient
import DdgProofs.Developability.Invariance

/-!
# Repulsive Curves: the discrete energy and its differential (paper Eq. 17–18; SPEC §3)

```
(k̂)_IJ = ¼ Σ_{i ∈ I} Σ_{j ∈ J} k(γ_i, γ_j, T_I)             (Eq. 18)
Ê(γ)   = Σ_{(I,J) ∈ 𝒫} (k̂)_IJ ℓ_I ℓ_J                         (Eq. 17)
```

over the ordered pairs `𝒫` of edges without a common vertex (`Disj`).

* `rc-energy`: the definitions; the energy is a sum of nonnegative terms, and every kernel
  evaluation has a positive denominator iff no two vertices of a pair in `𝒫` coincide.
* `rc-energy-scaling`: `Ê(cγ) = c^{2+α−β} Ê(γ)` (`energy_smul`).
* Invariances (SPEC §3.3): translations, orthogonal maps (rotations and reflections), reversing the
  orientation of any set of edges.
* `rc-energy-differential`: along any differentiable motion, the derivative of a pair term
  `e_IJ = k̂_IJ ℓ_I ℓ_J` is `Σ_a ∂_{i_a}·Γ'_{i_a} + Σ_b ∂_{j_b}·Γ'_{j_b}` with the SPEC §3.2 formulas
  (`hasDerivAt_pairEnergy`), and that of `Ê` is `Σ_v g_v·Γ'_v` for the assembled vertex gradient
  `g = energyGrad` (`hasDerivAt_energy`). Hypothesis: the pair is regular (`ℓ > 0` and `ν > 0` in
  the four kernel evaluations); at `ν = 0` the gradient extends continuously (`tendsto_kgrad`), but
  the derivative there is not proved here.
* `rc-gradient-identities`: `Σ_v g_v = 0`, `Σ_v γ_v × g_v = 0`, `Σ_v γ_v · g_v = (2 + α − β) Ê`
  for that gradient (`energy_gradient_identities`), from invariance and homogeneity along explicit
  curves (the rotation curves and their derivatives are those of `dev-gradient-invariance`,
  `Dev.rotX`, `Dev.rot_hasDerivAt`); and the same for any Fréchet gradient on `Fin n` vertices,
  through `Dev.sum_grad_eq_zero`, `Dev.sum_cross_grad_eq_zero` (`energy_gradient_identities_fderiv`).

The evaluation order of SPEC §3.2 (rows `R_I`, per-edge endpoint derivatives `G_{I,a}`) only
regroups the finite sums of `energyGrad`.
-/

noncomputable section

namespace DdgProofs

namespace RC

open Matrix Real Filter Topology

variable {V E : Type*} [Fintype V] [Fintype E] [DecidableEq V]

/-- `(I, J) ∈ 𝒫`: the edges have no common vertex. -/
def Disj (ep : E → Fin 2 → V) (I J : E) : Prop := ∀ a b, ep I a ≠ ep J b

instance (ep : E → Fin 2 → V) : DecidableRel (Disj ep) := fun I J => by
  unfold Disj; infer_instance

/-- The discrete kernel `(k̂)_IJ = ¼ Σ_{i∈I} Σ_{j∈J} k(γ_i, γ_j, T_I)` (Eq. 18).

[formula: rc-discrete-kernel] -/
def discreteKernel (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (I J : E) : ℝ :=
  (1 / 4) * ∑ a, ∑ b, kernel α β (γ (ep I a)) (γ (ep J b)) (tangent γ ep I)

/-- The term of one ordered pair, `e_IJ = (k̂)_IJ ℓ_I ℓ_J`. -/
def pairEnergy (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (I J : E) : ℝ :=
  discreteKernel α β γ ep I J * len γ ep I * len γ ep J

/-- The discrete tangent-point energy `Ê(γ) = Σ_{(I,J) ∈ 𝒫} (k̂)_IJ ℓ_I ℓ_J` (Eq. 17), both orders
of every pair.

[formula: rc-energy] -/
def energy (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) : ℝ :=
  ∑ I, ∑ J, if Disj ep I J then pairEnergy α β γ ep I J else 0

omit [Fintype V] [Fintype E] [DecidableEq V] in
/-- The discrete kernel (Eq. 18) is the mean of the kernel over the four endpoint pairs, hence
nonnegative.

[formula: rc-discrete-kernel] -/
theorem discreteKernel_nonneg (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (I J : E) :
    0 ≤ discreteKernel α β γ ep I J :=
  mul_nonneg (by norm_num) (Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ =>
    kernel_nonneg _ _ _ _ _)

omit [Fintype V] in
/-- **The energy is a sum of nonnegative terms, and its kernel evaluations are finite exactly when
no two vertices of a pair in `𝒫` coincide** (every denominator `|γ_i − γ_j|` is then positive).

[formula: rc-energy] -/
theorem energy_nonneg_and_finite (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) :
    0 ≤ energy α β γ ep ∧
      ((∀ I J, Disj ep I J → ∀ a b, 0 < vnorm (γ (ep I a) - γ (ep J b))) ↔
        ∀ I J, Disj ep I J → ∀ a b, γ (ep I a) ≠ γ (ep J b)) := by
  refine ⟨Finset.sum_nonneg fun I _ => Finset.sum_nonneg fun J _ => ?_, ?_⟩
  · split_ifs
    · exact mul_nonneg (mul_nonneg (discreteKernel_nonneg _ _ _ _ _ _) (len_nonneg _ _ _))
        (len_nonneg _ _ _)
    · exact le_rfl
  · simp only [vnorm_pos_iff, Ne, sub_eq_zero]

/-! ## Scaling and invariances -/

section invariance

variable (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V)

omit [Fintype V] [Fintype E] [DecidableEq V] in
/-- **Eq. 18 written out**: `(k̂)_IJ` is the trapezoid rule over the four endpoint pairs.

[formula: rc-discrete-kernel] -/
theorem discreteKernel_eq (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (I J : E) :
    discreteKernel α β γ ep I J = (1 / 4) *
      (kernel α β (γ (ep I 0)) (γ (ep J 0)) (tangent γ ep I)
        + kernel α β (γ (ep I 0)) (γ (ep J 1)) (tangent γ ep I)
        + kernel α β (γ (ep I 1)) (γ (ep J 0)) (tangent γ ep I)
        + kernel α β (γ (ep I 1)) (γ (ep J 1)) (tangent γ ep I)) := by
  simp only [discreteKernel, Fin.sum_univ_two]
  ring

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem edgeVec_smul (c : ℝ) (I : E) : edgeVec (c • γ) ep I = c • edgeVec γ ep I := by
  simp [edgeVec, smul_sub]

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem len_smul {c : ℝ} (hc : 0 < c) (I : E) : len (c • γ) ep I = c * len γ ep I := by
  rw [len, edgeVec_smul, vnorm_smul, abs_of_pos hc, len]

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem tangent_smul {c : ℝ} (hc : 0 < c) (I : E) : tangent (c • γ) ep I = tangent γ ep I := by
  rw [tangent, len_smul γ ep hc, edgeVec_smul, smul_smul, tangent, mul_inv, mul_comm c⁻¹,
    mul_assoc, inv_mul_cancel₀ hc.ne', mul_one]

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem discreteKernel_smul {c : ℝ} (hc : 0 < c) (I J : E) :
    discreteKernel α β (c • γ) ep I J = c ^ (α - β) * discreteKernel α β γ ep I J := by
  simp only [discreteKernel, tangent_smul γ ep hc, Pi.smul_apply, kernel_smul hc, ← Finset.mul_sum]
  ring

omit [Fintype V] in
/-- **Homogeneity of the energy**: `Ê(cγ) = c^{2+α−β} Ê(γ)` for `c > 0` (SPEC §3.3).

[formula: rc-energy-scaling] -/
theorem energy_smul {c : ℝ} (hc : 0 < c) :
    energy α β (c • γ) ep = c ^ (2 + α - β) * energy α β γ ep := by
  have hp : c ^ (2 + α - β) = c ^ (α - β) * c * c := by
    rw [show 2 + α - β = (α - β) + 2 by ring, Real.rpow_add hc, Real.rpow_two]
    ring
  simp only [energy, Finset.mul_sum, mul_ite, mul_zero]
  refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun J _ => ?_
  simp only [pairEnergy, discreteKernel_smul α β γ ep hc, len_smul γ ep hc, hp]
  split_ifs <;> ring

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem edgeVec_translate (a : V3) (I : E) :
    edgeVec (fun v => γ v + a) ep I = edgeVec γ ep I := by
  simp [edgeVec]

omit [Fintype V] in
/-- **Translation invariance** of the energy.

[formula: rc-energy] -/
theorem energy_translate (a : V3) : energy α β (fun v => γ v + a) ep = energy α β γ ep := by
  simp only [energy, pairEnergy, discreteKernel, len, tangent, edgeVec_translate,
    kernel_translate]

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem edgeVec_orth (R : Matrix (Fin 3) (Fin 3) ℝ) (I : E) :
    edgeVec (fun v => R *ᵥ γ v) ep I = R *ᵥ edgeVec γ ep I := by
  simp [edgeVec, mulVec_sub]

omit [Fintype V] in
/-- **Invariance under orthogonal maps** (`Rᵀ R = 1`: rotations and reflections).

[formula: rc-energy] -/
theorem energy_orth {R : Matrix (Fin 3) (Fin 3) ℝ} (hR : Rᵀ * R = 1) :
    energy α β (fun v => R *ᵥ γ v) ep = energy α β γ ep := by
  have hl : ∀ I, len (fun v => R *ᵥ γ v) ep I = len γ ep I := fun I => by
    rw [len, edgeVec_orth, orth_vnorm hR, len]
  have ht : ∀ I, tangent (fun v => R *ᵥ γ v) ep I = R *ᵥ tangent γ ep I := fun I => by
    rw [tangent, hl, edgeVec_orth, tangent, mulVec_smul]
  simp only [energy, pairEnergy, discreteKernel, hl, ht, kernel_orth hR]

@[simp] theorem rev_zero2 : Fin.rev (0 : Fin 2) = 1 := rfl
@[simp] theorem rev_one2 : Fin.rev (1 : Fin 2) = 0 := rfl

/-- Reverse the orientation of the edges `I` with `o I = true`. -/
def reorient (ep : E → Fin 2 → V) (o : E → Bool) (I : E) (a : Fin 2) : V :=
  if o I then ep I (Fin.rev a) else ep I a

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem disj_reorient (o : E → Bool) (I J : E) : Disj (reorient ep o) I J ↔ Disj ep I J := by
  have key : ∀ (K : E) (x : V), (∃ a, reorient ep o K a = x) ↔ ∃ a, ep K a = x := by
    intro K x
    unfold reorient
    split_ifs
    · exact ⟨fun ⟨a, h⟩ => ⟨_, h⟩, fun ⟨a, h⟩ => ⟨Fin.rev a, by simpa using h⟩⟩
    · exact Iff.rfl
  constructor
  · intro h a b hab
    obtain ⟨a', ha'⟩ := (key I (ep I a)).mpr ⟨a, rfl⟩
    obtain ⟨b', hb'⟩ := (key J (ep J b)).mpr ⟨b, rfl⟩
    exact h a' b' (by rw [ha', hb', hab])
  · intro h a b hab
    obtain ⟨a', ha'⟩ := (key I (reorient ep o I a)).mp ⟨a, rfl⟩
    obtain ⟨b', hb'⟩ := (key J (reorient ep o J b)).mp ⟨b, rfl⟩
    exact h a' b' (by rw [ha', hb', hab])

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem len_reorient (o : E → Bool) (I : E) : len γ (reorient ep o) I = len γ ep I := by
  unfold len edgeVec reorient
  split_ifs
  · simp only [rev_one2, rev_zero2]
    exact vnorm_sub_comm _ _
  · rfl

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem tangent_reorient (o : E → Bool) (I : E) :
    tangent γ (reorient ep o) I = if o I then -tangent γ ep I else tangent γ ep I := by
  rw [tangent, len_reorient, tangent]
  unfold edgeVec reorient
  split_ifs
  · simp only [rev_one2, rev_zero2]
    rw [← neg_sub, smul_neg]
  · rfl

omit [Fintype V] in
/-- **Invariance under re-orientation** of any set of edges (`k` is even in `T`, and the double
sum runs over both endpoints of both edges).

[formula: rc-energy] -/
theorem energy_reorient (o : E → Bool) :
    energy α β γ (reorient ep o) = energy α β γ ep := by
  unfold energy
  refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun J _ => ?_
  simp only [disj_reorient, pairEnergy, len_reorient]
  congr 2
  unfold discreteKernel
  rw [tangent_reorient]
  have hT : ∀ p q, kernel α β p q (if o I then -tangent γ ep I else tangent γ ep I) =
      kernel α β p q (tangent γ ep I) := fun p q => by
    split_ifs
    · exact kernel_neg_tangent _ _ _ _ _
    · rfl
  simp only [hT]
  congr 1
  unfold reorient
  simp only [Fin.sum_univ_two, rev_zero2, rev_one2]
  split_ifs <;> ring

end invariance

/-! ## The differential -/

section differential

variable (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V)

/-- A pair `(I, J)` is regular at `γ`: both lengths are positive and `ν = |T_I × (γ_i − γ_j)| > 0`
for the four endpoint pairs. -/
def PairRegular (I J : E) : Prop :=
  0 < len γ ep I ∧ 0 < len γ ep J ∧
    ∀ a b, crossProduct (tangent γ ep I) (γ (ep I a) - γ (ep J b)) ≠ 0

/-- `g_T = ¼ Σ_i Σ_j ∇_T k(γ_i, γ_j, T_I)`. -/
def avgGradT (I J : E) : V3 :=
  (1 / 4 : ℝ) • ∑ a, ∑ b, kgradT α β (γ (ep I a)) (γ (ep J b)) (tangent γ ep I)

/-- `∂e_IJ/∂γ_{i_a} = s_a ℓ_J (k̂_IJ T_I + g_T) + ℓ_I ℓ_J ¼ Σ_{j ∈ J} ∇_p k(γ_{i_a}, γ_j, T_I)`
(SPEC §3.2).

[formula: rc-energy-differential] -/
def pairGradI (I J : E) (a : Fin 2) : V3 :=
  (sgn a * len γ ep J) • (discreteKernel α β γ ep I J • tangent γ ep I + avgGradT α β γ ep I J)
    + (len γ ep I * len γ ep J / 4) • ∑ b, kgradP α β (γ (ep I a)) (γ (ep J b)) (tangent γ ep I)

/-- `∂e_IJ/∂γ_{j_b} = s_b ℓ_I k̂_IJ T_J + ℓ_I ℓ_J ¼ Σ_{i ∈ I} ∇_q k(γ_i, γ_{j_b}, T_I)`, with
`∇_q k = −∇_p k` (SPEC §3.2).

[formula: rc-energy-differential] -/
def pairGradJ (I J : E) (b : Fin 2) : V3 :=
  (sgn b * len γ ep I * discreteKernel α β γ ep I J) • tangent γ ep J
    - (len γ ep I * len γ ep J / 4) • ∑ a, kgradP α β (γ (ep I a)) (γ (ep J b)) (tangent γ ep I)

omit [Fintype V] [Fintype E] [DecidableEq V] in
theorem kgradT_dot_proj {α β : ℝ} {p q T : V3} (hT : T ⬝ᵥ T = 1) (w : V3) :
    kgradT α β p q T ⬝ᵥ proj T w = kgradT α β p q T ⬝ᵥ w := by
  rw [dot_proj]
  have : kgradT α β p q T ⬝ᵥ T = 0 := by
    rw [kgradT, smul_dotProduct, dotProduct_comm (proj T (p - q)) T, dot_proj_self hT, smul_zero]
  rw [this, mul_zero, sub_zero]

variable {α β} {ep}

omit [Fintype V] [Fintype E] [DecidableEq V] in
/-- **The differential of one pair term** (SPEC §3.2): along a differentiable motion `Γ` with the
pair regular at `Γ t`, `d/ds e_IJ(Γ s) = Σ_a ∂_{i_a}·Γ'_{i_a} + Σ_b ∂_{j_b}·Γ'_{j_b}`.

[formula: rc-energy-differential] -/
theorem hasDerivAt_pairEnergy {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {I J : E}
    (hreg : PairRegular (Γ t) ep I J) :
    HasDerivAt (fun s => pairEnergy α β (Γ s) ep I J)
      (∑ a, pairGradI α β (Γ t) ep I J a ⬝ᵥ Γ' (ep I a)
        + ∑ b, pairGradJ α β (Γ t) ep I J b ⬝ᵥ Γ' (ep J b)) t := by
  obtain ⟨hI, hJ, hν⟩ := hreg
  have hT := hasDerivAt_tangent (ep := ep) (I := I) hΓ hI
  have hunit := tangent_dot_self hI
  have hk : ∀ a b, HasDerivAt
      (fun s => kernel α β (Γ s (ep I a)) (Γ s (ep J b)) (tangent (Γ s) ep I))
      (kgradP α β (Γ t (ep I a)) (Γ t (ep J b)) (tangent (Γ t) ep I) ⬝ᵥ (Γ' (ep I a) - Γ' (ep J b))
        + kgradT α β (Γ t (ep I a)) (Γ t (ep J b)) (tangent (Γ t) ep I) ⬝ᵥ
          ((len (Γ t) ep I)⁻¹ • ∑ c, sgn c • proj (tangent (Γ t) ep I) (Γ' (ep I c)))) t :=
    fun a b => hasDerivAt_kernel (hΓ _) (hΓ _) hT hunit (tangent_dot_deriv hI) (hν a b)
  have hkhat : HasDerivAt (fun s => discreteKernel α β (Γ s) ep I J)
      ((1 / 4) * ∑ a, ∑ b, (kgradP α β (Γ t (ep I a)) (Γ t (ep J b)) (tangent (Γ t) ep I) ⬝ᵥ
          (Γ' (ep I a) - Γ' (ep J b))
        + kgradT α β (Γ t (ep I a)) (Γ t (ep J b)) (tangent (Γ t) ep I) ⬝ᵥ
          ((len (Γ t) ep I)⁻¹ • ∑ c, sgn c • proj (tangent (Γ t) ep I) (Γ' (ep I c))))) t := by
    have := (HasDerivAt.sum (u := Finset.univ) fun a _ =>
      HasDerivAt.sum (u := Finset.univ) fun b _ => hk a b).const_mul (1 / 4 : ℝ)
    unfold discreteKernel
    convert this using 1
    funext s
    simp
  have hlI := hasDerivAt_len (ep := ep) (I := I) hΓ hI
  have hlJ := hasDerivAt_len (ep := ep) (I := J) hΓ hJ
  have := (hkhat.mul hlI).mul hlJ
  refine this.congr_deriv ?_
  have hl0 : len (Γ t) ep I ≠ 0 := hI.ne'
  simp only [dotProduct_smul, dotProduct_sum, kgradT_dot_proj hunit, smul_eq_mul, Pi.mul_apply]
  simp only [pairGradI, pairGradJ, avgGradT, add_dotProduct, sub_dotProduct, smul_dotProduct,
    dotProduct_sub, smul_eq_mul, Fin.sum_univ_two, sgn_zero, sgn_one]
  field_simp
  ring

/-- The assembled vertex gradient `g_v = ∂Ê/∂γ_v = Σ_{(I,J) ∈ 𝒫} (Σ_{a : i_a = v} ∂_{i_a} e_IJ +
Σ_{b : j_b = v} ∂_{j_b} e_IJ)`.

[formula: rc-energy-differential] -/
def energyGrad (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (v : V) : V3 :=
  ∑ I, ∑ J, if Disj ep I J then
    (∑ a, if ep I a = v then pairGradI α β γ ep I J a else 0)
      + ∑ b, if ep J b = v then pairGradJ α β γ ep I J b else 0
  else 0

/-- **The differential of the energy** (SPEC §3.2): along a differentiable motion `Γ` with every
pair of `𝒫` regular at `Γ t`, `d/ds Ê(Γ s) = Σ_v g_v · Γ'_v` with `g = energyGrad`.

[formula: rc-energy-differential] -/
theorem hasDerivAt_energy {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (hreg : ∀ I J, Disj ep I J → PairRegular (Γ t) ep I J) :
    HasDerivAt (fun s => energy α β (Γ s) ep) (∑ v, energyGrad α β (Γ t) ep v ⬝ᵥ Γ' v) t := by
  have hpair : ∀ I J, HasDerivAt
      (fun s => if Disj ep I J then pairEnergy α β (Γ s) ep I J else 0)
      (if Disj ep I J then
        (∑ a, pairGradI α β (Γ t) ep I J a ⬝ᵥ Γ' (ep I a)
          + ∑ b, pairGradJ α β (Γ t) ep I J b ⬝ᵥ Γ' (ep J b)) else 0) t := by
    intro I J
    by_cases hd : Disj ep I J
    · simp only [hd, ite_true]
      exact hasDerivAt_pairEnergy hΓ (hreg I J hd)
    · simp only [hd, ite_false]
      exact hasDerivAt_const _ _
  have h := HasDerivAt.fun_sum (u := Finset.univ) fun I _ =>
    HasDerivAt.fun_sum (u := Finset.univ) fun J _ => hpair I J
  have hfun : (fun s => energy α β (Γ s) ep) =
      fun s => ∑ I, ∑ J, if Disj ep I J then pairEnergy α β (Γ s) ep I J else 0 := rfl
  rw [hfun]
  refine h.congr_deriv ?_
  have ite_dot : ∀ (c : Prop) [Decidable c] (x w : V3),
      (if c then x else 0) ⬝ᵥ w = if c then x ⬝ᵥ w else 0 := by
    intro c _ x w; split_ifs <;> simp
  symm
  unfold energyGrad
  simp only [sum_dotProduct]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun J _ => ?_
  by_cases hd : Disj ep I J
  · simp only [hd, ite_true, add_dotProduct, Finset.sum_add_distrib, sum_dotProduct, ite_dot]
    congr 1
    · rw [Finset.sum_comm]
      simp only [Finset.sum_ite_eq, Finset.mem_univ, ite_true]
    · rw [Finset.sum_comm]
      simp only [Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  · simp [hd]

end differential

/-! ## Invariance and homogeneity ⇒ the gradient identities -/

section identities

/-- The rotation `Dev.rotX s` as a matrix. -/
def rotMatX (s : ℝ) : Matrix (Fin 3) (Fin 3) ℝ := !![1, 0, 0; 0, cos s, -sin s; 0, sin s, cos s]
/-- The rotation `Dev.rotY s` as a matrix. -/
def rotMatY (s : ℝ) : Matrix (Fin 3) (Fin 3) ℝ := !![cos s, 0, sin s; 0, 1, 0; -sin s, 0, cos s]
/-- The rotation `Dev.rotZ s` as a matrix. -/
def rotMatZ (s : ℝ) : Matrix (Fin 3) (Fin 3) ℝ := !![cos s, -sin s, 0; sin s, cos s, 0; 0, 0, 1]

theorem rotX_eq (s : ℝ) (p : V3) : Dev.rotX s p = rotMatX s *ᵥ p := by
  ext i
  fin_cases i <;> simp [Dev.rotX, rotMatX, mulVec, dotProduct, Fin.sum_univ_three]
  ring

theorem rotY_eq (s : ℝ) (p : V3) : Dev.rotY s p = rotMatY s *ᵥ p := by
  ext i
  fin_cases i <;> simp [Dev.rotY, rotMatY, mulVec, dotProduct, Fin.sum_univ_three]

theorem rotZ_eq (s : ℝ) (p : V3) : Dev.rotZ s p = rotMatZ s *ᵥ p := by
  ext i
  fin_cases i <;> simp [Dev.rotZ, rotMatZ, mulVec, dotProduct, Fin.sum_univ_three]
  ring

theorem rotMatX_orth (s : ℝ) : (rotMatX s)ᵀ * rotMatX s = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [rotMatX, Matrix.mul_apply, Fin.sum_univ_three] <;> nlinarith [sin_sq_add_cos_sq s]

theorem rotMatY_orth (s : ℝ) : (rotMatY s)ᵀ * rotMatY s = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [rotMatY, Matrix.mul_apply, Fin.sum_univ_three] <;> nlinarith [sin_sq_add_cos_sq s]

theorem rotMatZ_orth (s : ℝ) : (rotMatZ s)ᵀ * rotMatZ s = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [rotMatZ, Matrix.mul_apply, Fin.sum_univ_three] <;> nlinarith [sin_sq_add_cos_sq s]

theorem hasDerivAt_rotX (p : V3) :
    HasDerivAt (fun s => Dev.rotX s p) ![0, -p 2, p 1] 0 := by
  rw [hasDerivAt_pi]
  intro c
  fin_cases c
  · simpa [Dev.rotX] using hasDerivAt_const (0 : ℝ) (p 0)
  · simpa [Dev.rotX] using (Dev.rot_hasDerivAt (p 1) (p 2)).1
  · simpa [Dev.rotX] using (Dev.rot_hasDerivAt (p 1) (p 2)).2

theorem hasDerivAt_rotY (p : V3) :
    HasDerivAt (fun s => Dev.rotY s p) ![p 2, 0, -p 0] 0 := by
  rw [hasDerivAt_pi]
  intro c
  fin_cases c
  · simpa [Dev.rotY] using (Dev.rotY_hasDerivAt (p 0) (p 2)).1
  · simpa [Dev.rotY] using hasDerivAt_const (0 : ℝ) (p 1)
  · simpa [Dev.rotY] using (Dev.rotY_hasDerivAt (p 0) (p 2)).2

theorem hasDerivAt_rotZ (p : V3) :
    HasDerivAt (fun s => Dev.rotZ s p) ![-p 1, p 0, 0] 0 := by
  rw [hasDerivAt_pi]
  intro c
  fin_cases c
  · simpa [Dev.rotZ] using (Dev.rot_hasDerivAt (p 0) (p 1)).1
  · simpa [Dev.rotZ] using (Dev.rot_hasDerivAt (p 0) (p 1)).2
  · simpa [Dev.rotZ] using hasDerivAt_const (0 : ℝ) (p 2)

omit [Fintype E] [DecidableEq V] in
/-- **Invariance and homogeneity ⇒ the gradient identities**, for a function `F` of the positions
whose derivative along every differentiable motion through `γ` is `Σ_v g_v·Γ'_v`: invariant under
translations and orthogonal maps, homogeneous of degree `p` ⇒ `Σ g = 0`, `Σ γ × g = 0`,
`Σ γ·g = p F(γ)` (the curves are the translations, `Dev.rotX/Y/Z` and `s ↦ (1+s)γ`).

[formula: rc-gradient-identities] -/
theorem gradient_identities_of_curves {F : (V → V3) → ℝ} {γ g : V → V3} {p : ℝ}
    (hd : ∀ (Γ : ℝ → V → V3) (Γ' : V → V3), Γ 0 = γ →
      (∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) 0) →
        HasDerivAt (fun s => F (Γ s)) (∑ v, g v ⬝ᵥ Γ' v) 0)
    (htr : ∀ a : V3, F (fun v => γ v + a) = F γ)
    (hrot : ∀ R : Matrix (Fin 3) (Fin 3) ℝ, Rᵀ * R = 1 → F (fun v => R *ᵥ γ v) = F γ)
    (hhom : ∀ c : ℝ, 0 < c → F (c • γ) = c ^ p * F γ) :
    ∑ v, g v = 0 ∧ ∑ v, crossProduct (γ v) (g v) = 0 ∧ ∑ v, γ v ⬝ᵥ g v = p * F γ := by
  -- a curve along which `F` is constant has derivative `0`
  have const0 : ∀ (Γ : ℝ → V → V3) (Γ' : V → V3), Γ 0 = γ →
      (∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) 0) → (∀ s, F (Γ s) = F γ) →
        ∑ v, g v ⬝ᵥ Γ' v = 0 := by
    intro Γ Γ' h0 hΓ hc
    have h1 := hd Γ Γ' h0 hΓ
    have h2 : HasDerivAt (fun s => F (Γ s)) 0 0 := by
      simpa [hc] using hasDerivAt_const (0 : ℝ) (F γ)
    exact h1.unique h2
  refine ⟨?_, ?_, ?_⟩
  · funext c
    have := const0 (fun s v => γ v + s • Pi.single c (1 : ℝ)) (fun _ => Pi.single c 1)
      (by simp) (fun v => by simpa using ((hasDerivAt_id (0 : ℝ)).smul_const
        (Pi.single c (1 : ℝ) : V3)).const_add (γ v)) (fun s => htr _)
    simpa [dotProduct, Pi.single_apply, Finset.sum_apply] using this
  · have dX := const0 (fun s v => Dev.rotX s (γ v)) (fun v => ![0, -γ v 2, γ v 1])
      (by funext v; ext i; fin_cases i <;> simp [Dev.rotX]) (fun v => hasDerivAt_rotX (γ v))
      (fun s => by simp only [rotX_eq]; exact hrot _ (rotMatX_orth s))
    have dY := const0 (fun s v => Dev.rotY s (γ v)) (fun v => ![γ v 2, 0, -γ v 0])
      (by funext v; ext i; fin_cases i <;> simp [Dev.rotY]) (fun v => hasDerivAt_rotY (γ v))
      (fun s => by simp only [rotY_eq]; exact hrot _ (rotMatY_orth s))
    have dZ := const0 (fun s v => Dev.rotZ s (γ v)) (fun v => ![-γ v 1, γ v 0, 0])
      (by funext v; ext i; fin_cases i <;> simp [Dev.rotZ]) (fun v => hasDerivAt_rotZ (γ v))
      (fun s => by simp only [rotZ_eq]; exact hrot _ (rotMatZ_orth s))
    simp only [dotProduct, Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] at dX dY dZ
    funext c
    rw [Finset.sum_apply]
    fin_cases c <;> simp only [cross_apply, Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons, Pi.zero_apply]
    · rw [← dX]; exact Finset.sum_congr rfl fun v _ => by ring
    · rw [← dY]; exact Finset.sum_congr rfl fun v _ => by ring
    · rw [← dZ]; exact Finset.sum_congr rfl fun v _ => by ring
  · -- the scaling curve `s ↦ (1 + s) γ`
    have h1 := hd (fun s : ℝ => (1 + s) • γ) γ (by simp) (fun v => by
      simpa using ((hasDerivAt_id (0 : ℝ)).const_add 1).smul_const (γ v))
    have h2 : HasDerivAt (fun s : ℝ => F ((1 + s) • γ)) (p * F γ) 0 := by
      have hpow : HasDerivAt (fun s : ℝ => (1 + s) ^ p * F γ) (p * F γ) 0 := by
        have := (((hasDerivAt_id (0 : ℝ)).const_add 1).rpow_const (p := p)
          (Or.inl (by norm_num))).mul_const (F γ)
        simpa using this
      refine hpow.congr_of_eventuallyEq ?_
      have : ∀ᶠ s in 𝓝 (0 : ℝ), -1 < s := lt_mem_nhds (by norm_num)
      filter_upwards [this] with s hs using hhom _ (by linarith)
    have := h1.unique h2
    rw [← this]
    exact Finset.sum_congr rfl fun v _ => dotProduct_comm _ _

/-- **The gradient identities of the energy** (SPEC §3.3): if every pair of `𝒫` is regular at `γ`,
the vertex gradient `g = energyGrad` satisfies `Σ_v g_v = 0`, `Σ_v γ_v × g_v = 0` and
`Σ_v γ_v · g_v = (2 + α − β) Ê(γ)`.

[formula: rc-gradient-identities] -/
theorem energy_gradient_identities {α β : ℝ} {γ : V → V3} {ep : E → Fin 2 → V}
    (hreg : ∀ I J, Disj ep I J → PairRegular γ ep I J) :
    ∑ v, energyGrad α β γ ep v = 0 ∧
      ∑ v, crossProduct (γ v) (energyGrad α β γ ep v) = 0 ∧
      ∑ v, γ v ⬝ᵥ energyGrad α β γ ep v = (2 + α - β) * energy α β γ ep :=
  gradient_identities_of_curves (F := fun γ => energy α β γ ep)
    (fun Γ Γ' h0 hΓ => by
      subst h0
      exact hasDerivAt_energy hΓ hreg)
    (fun a => energy_translate α β γ ep a)
    (fun _ hR => energy_orth α β γ ep hR)
    (fun _ hc => energy_smul α β γ ep hc)

/-- **The gradient identities for a Fréchet gradient**, through `dev-gradient-invariance`
(`Dev.sum_grad_eq_zero`, `Dev.sum_cross_grad_eq_zero`) on `Fin n` vertices: if `Ê` has a Fréchet
derivative at `γ` represented by `g`, then `Σ g = 0` and `Σ γ × g = 0`.

[formula: rc-gradient-identities] -/
theorem energy_gradient_identities_fderiv {n : ℕ} {E : Type*} [Fintype E] {α β : ℝ}
    {ep : E → Fin 2 → Fin n} {γ : Dev.Positions n} {E' : Dev.Positions n →L[ℝ] ℝ}
    (hE : HasFDerivAt (fun γ => energy α β γ ep) E' γ) (g : Dev.Positions n)
    (hg : ∀ δ, E' δ = ∑ v, ∑ c, g v c * δ v c) :
    ∑ v, g v = 0 ∧ ∑ v, crossProduct (γ v) (g v) = 0 := by
  refine ⟨Dev.sum_grad_eq_zero hE g hg fun a s => ?_, Dev.sum_cross_grad_eq_zero hE g hg ?_ ?_ ?_⟩
  · have : γ + s • (fun _ => a) = fun v => γ v + s • a := by funext v; simp
    rw [this]; exact energy_translate α β γ ep (s • a)
  · intro s; simp only [rotX_eq]; exact energy_orth α β γ ep (rotMatX_orth s)
  · intro s; simp only [rotY_eq]; exact energy_orth α β γ ep (rotMatY_orth s)
  · intro s; simp only [rotZ_eq]; exact energy_orth α β γ ep (rotMatZ_orth s)

end identities

end RC

end DdgProofs
