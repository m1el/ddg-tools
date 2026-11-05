import DDG.Repulsive.Surfaces.Energy

/-!
# Repulsive Surfaces: the Barnes–Hut energy (paper Eq. 17, Sec. 4; SPEC §3.4, D4, D5)

The tree is abstract: faces `F`, leaves `L` with `leafOf : F → L` (every face in exactly one
leaf), nodes `N` with face sets `nodeFaces C`; per leaf `ℓ` a set of admissible nodes `adm ℓ` and
a set of inadmissible ("near") leaves `near ℓ`. Node aggregates `A_C = Σ_{T∈C} a_T`,
`X_C = Σ_{T∈C} a_T X_T / A_C`.

```
Ẽ = Σ_ℓ Σ_{S ∈ ℓ} [ Σ_{C ∈ adm ℓ} K(X_S, N_S; X_C) a_S A_C + Σ_{ℓ' ∈ near ℓ} Σ_{T ∈ ℓ'} [S ≠ T] K(X_S, N_S; X_T) a_S a_T ]
```

* `rs-bh-energy`: if, for every leaf, the admissible nodes and the near leaves cover every face
  exactly once, `Ẽ = Ê + Σ_ℓ Σ_{S∈ℓ} Σ_{C∈adm ℓ} (lumped − exact)` (`bh_energy_split`): the only
  approximation is the lumping of admissible clusters, which is exact for a one-face cluster; with
  the admissibility rule `max(r_ℓ, r_C) < θ d_box` and `θ = 0` no node is admissible
  (`no_admissible_theta_zero`), and then `Ẽ = Ê` (`bh_energy_theta_zero`).
* `rs-bh-aggregate-adjoint` (D4): with `A_C = Σ_{T∈C} a_T` and `M_C = Σ_{T∈C} a_T X_T`, the
  derivative of any function of the aggregates is
  `Σ_C (∂_{A_C} A_C' + ∂_{M_C}·M_C') = Σ_T ((Σ_{C∋T} ∂_{A_C}) a_T' + (Σ_{C∋T} ∂_{M_C})·(a_T X_T)')`:
  pushing the cluster derivatives down to the member faces gives the exact gradient of `Ẽ`.
* `rs-box-distance`: `d² = Σ_k max(0, lo_B,k − hi_A,k, lo_A,k − hi_B,k)²` is at most `|x − y|²`
  for all `x ∈ A`, `y ∈ B`; it is `0` for overlapping boxes, which are therefore never admissible.
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Real Finset
open RC (V3 vnorm vnorm_nonneg)

/-! ## The energy with the leaf rule -/

section energy

variable {V F L N : Type*} [Fintype F] [Fintype L] [Fintype N] [DecidableEq F] [DecidableEq L]

/-- The faces of leaf `ℓ`. -/
def leafFaces (leafOf : F → L) (ℓ : L) : Finset F := univ.filter (leafOf · = ℓ)

/-- The node area `A_C = Σ_{T∈C} a_T`. -/
def nodeArea (a : F → ℝ) (nodeFaces : N → Finset F) (C : N) : ℝ := ∑ T ∈ nodeFaces C, a T

/-- The node barycenter `X_C = Σ_{T∈C} a_T X_T / A_C`. -/
def nodeBary (a : F → ℝ) (X : F → V3) (nodeFaces : N → Finset F) (C : N) : V3 :=
  (nodeArea a nodeFaces C)⁻¹ • ∑ T ∈ nodeFaces C, a T • X T

/-- The lumped far-field term of face `S` and node `C`: `K(X_S, N_S; X_C) a_S A_C`. -/
def lumped (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (nodeFaces : N → Finset F) (S : F)
    (C : N) : ℝ :=
  kernel α β (bary f fv S) (normal f fv S) (nodeBary (area f fv) (bary f fv) nodeFaces C) *
    area f fv S * nodeArea (area f fv) nodeFaces C

/-- **The Barnes–Hut energy with the leaf rule** (Eq. 17, SPEC §3.4, D5).

[formula: rs-bh-energy] -/
def bhEnergy (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (leafOf : F → L)
    (nodeFaces : N → Finset F) (adm : L → Finset N) (near : L → Finset L) : ℝ :=
  ∑ ℓ, ∑ S ∈ leafFaces leafOf ℓ,
    (∑ C ∈ adm ℓ, lumped α β f fv nodeFaces S C +
      ∑ ℓ' ∈ near ℓ, ∑ T ∈ leafFaces leafOf ℓ', if S ≠ T then pairEnergy α β f fv S T else 0)

/-- Every leaf's admissible nodes and near leaves cover each face exactly once. -/
def ExactCover (leafOf : F → L) (nodeFaces : N → Finset F) (adm : L → Finset N)
    (near : L → Finset L) : Prop :=
  ∀ ℓ T, (∑ C ∈ adm ℓ, if T ∈ nodeFaces C then (1 : ℝ) else 0) +
    (∑ ℓ' ∈ near ℓ, if leafOf T = ℓ' then (1 : ℝ) else 0) = 1

omit [Fintype L] [DecidableEq F] in
theorem sum_leafFaces (leafOf : F → L) (s : Finset L) (g : F → ℝ) :
    ∑ ℓ' ∈ s, ∑ T ∈ leafFaces leafOf ℓ', g T =
      ∑ T, (∑ ℓ' ∈ s, if leafOf T = ℓ' then (1 : ℝ) else 0) * g T := by
  simp only [leafFaces, Finset.sum_filter, Finset.sum_mul, ite_mul, one_mul, zero_mul]
  exact Finset.sum_comm

omit [Fintype L] [DecidableEq L] [Fintype N] in
theorem sum_nodeFaces (nodeFaces : N → Finset F) (s : Finset N) (g : F → ℝ) :
    ∑ C ∈ s, ∑ T ∈ nodeFaces C, g T =
      ∑ T, (∑ C ∈ s, if T ∈ nodeFaces C then (1 : ℝ) else 0) * g T := by
  simp only [Finset.sum_mul, ite_mul, one_mul, zero_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun C _ => ?_
  rw [Finset.sum_ite_mem, Finset.univ_inter]

omit [Fintype N] in
/-- **The Barnes–Hut energy is the exact energy plus the lumping error** (SPEC §3.4): under an
exact cover, `Ẽ = Ê + Σ_ℓ Σ_{S∈ℓ} Σ_{C ∈ adm ℓ} (K(X_S, N_S; X_C) a_S A_C −
Σ_{T∈C} [S ≠ T] K(X_S, N_S; X_T) a_S a_T)`.

[formula: rs-bh-energy] -/
theorem bh_energy_split (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (leafOf : F → L)
    (nodeFaces : N → Finset F) (adm : L → Finset N) (near : L → Finset L)
    (hcov : ExactCover leafOf nodeFaces adm near) :
    bhEnergy α β f fv leafOf nodeFaces adm near = energy α β f fv +
      ∑ ℓ, ∑ S ∈ leafFaces leafOf ℓ, ∑ C ∈ adm ℓ, (lumped α β f fv nodeFaces S C -
        ∑ T ∈ nodeFaces C, if S ≠ T then pairEnergy α β f fv S T else 0) := by
  have key : ∀ ℓ S, ∑ C ∈ adm ℓ, lumped α β f fv nodeFaces S C +
      ∑ ℓ' ∈ near ℓ, ∑ T ∈ leafFaces leafOf ℓ', (if S ≠ T then pairEnergy α β f fv S T else 0) =
      ∑ T, (if S ≠ T then pairEnergy α β f fv S T else 0) +
        ∑ C ∈ adm ℓ, (lumped α β f fv nodeFaces S C -
          ∑ T ∈ nodeFaces C, if S ≠ T then pairEnergy α β f fv S T else 0) := by
    intro ℓ S
    rw [sum_leafFaces, Finset.sum_sub_distrib, sum_nodeFaces]
    have : ∀ T, (∑ ℓ' ∈ near ℓ, if leafOf T = ℓ' then (1 : ℝ) else 0) =
        1 - ∑ C ∈ adm ℓ, if T ∈ nodeFaces C then (1 : ℝ) else 0 := fun T => by
      linarith [hcov ℓ T]
    simp only [this, sub_mul, one_mul, Finset.sum_sub_distrib]
    ring
  unfold bhEnergy
  simp only [key, Finset.sum_add_distrib]
  congr 1
  rw [energy]
  exact Finset.sum_fiberwise univ leafOf _

omit [Fintype F] [Fintype N] in
/-- **A one-face cluster is lumped exactly**: if `C = {T}` with `T ≠ S`, the far-field term is the
exact pair term.

[formula: rs-bh-energy] -/
theorem lumped_singleton (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (nodeFaces : N → Finset F)
    {S T : F} {C : N} (hC : nodeFaces C = {T}) (hT : 0 < area f fv T) (hST : S ≠ T) :
    lumped α β f fv nodeFaces S C =
      ∑ T' ∈ nodeFaces C, if S ≠ T' then pairEnergy α β f fv S T' else 0 := by
  simp only [lumped, nodeBary, nodeArea, hC, Finset.sum_singleton, hST, ne_eq, not_false_eq_true,
    ite_true, pairEnergy, smul_smul, inv_mul_cancel₀ hT.ne', one_smul]

/-- **At `θ = 0` no node is admissible**: `max(r_ℓ, r_C) < 0 · d` is impossible for radii `≥ 0`.

[formula: rs-bh-energy] -/
theorem no_admissible_theta_zero {rℓ rC d : ℝ} (hℓ : 0 ≤ rℓ) : ¬ max rℓ rC < 0 * d := by
  rw [zero_mul, not_lt]
  exact le_max_of_le_left hℓ

omit [Fintype N] in
/-- **`θ = 0` is the exact energy** (SPEC §3.4): with no admissible node and every leaf near every
leaf, `Ẽ = Ê`.

[formula: rs-bh-energy] -/
theorem bh_energy_theta_zero (α β : ℝ) (f : V → V3) (fv : F → Fin 3 → V) (leafOf : F → L)
    (nodeFaces : N → Finset F) :
    bhEnergy α β f fv leafOf nodeFaces (fun _ => ∅) (fun _ => univ) = energy α β f fv := by
  rw [bh_energy_split α β f fv leafOf nodeFaces _ _ (fun ℓ T => by simp)]
  simp

end energy

/-! ## The aggregate adjoint (D4) -/

section adjoint

variable {F N : Type*} [Fintype F] [Fintype N] [DecidableEq F]

/-- **Pushing cluster derivatives down to the faces** (SPEC §3.4, D4): for `A_C = Σ_{T∈C} a_T`
and `M_C = Σ_{T∈C} m_T` (`m_T = a_T X_T`),
`Σ_C (g^A_C A_C' + g^M_C · M_C') = Σ_T ((Σ_{C ∋ T} g^A_C) a_T' + (Σ_{C ∋ T} g^M_C) · m_T')`.

[formula: rs-bh-aggregate-adjoint] -/
theorem aggregate_adjoint (nodeFaces : N → Finset F) (gA : N → ℝ) (gM : N → V3) (a' : F → ℝ)
    (m' : F → V3) :
    ∑ C, (gA C * ∑ T ∈ nodeFaces C, a' T + gM C ⬝ᵥ ∑ T ∈ nodeFaces C, m' T) =
      ∑ T, ((∑ C, if T ∈ nodeFaces C then gA C else 0) * a' T +
        (∑ C, if T ∈ nodeFaces C then gM C else 0) ⬝ᵥ m' T) := by
  have hs : ∀ {β : Type} [AddCommMonoid β] (s : Finset F) (g : F → β),
      ∑ T ∈ s, g T = ∑ T, if T ∈ s then g T else 0 := by
    intro β _ s g; rw [Finset.sum_ite_mem, Finset.univ_inter]
  simp only [hs (nodeFaces _)]
  simp only [Finset.mul_sum, dotProduct_sum, Finset.sum_add_distrib, Finset.sum_mul,
    sum_dotProduct]
  congr 1 <;>
  · rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun T _ => Finset.sum_congr rfl fun C _ => ?_
    split_ifs <;> simp

/-- **The chain rule through the aggregates** (D4): if `G` is differentiable in the aggregates
with partials `g^A`, `g^M` along every motion of them, then along motions `a_T(s)`, `m_T(s)` of
the face data, `d/ds G(A(s), M(s)) = Σ_T ((Σ_{C∋T} g^A_C) a_T' + (Σ_{C∋T} g^M_C)·m_T')`.

[formula: rs-bh-aggregate-adjoint] -/
theorem hasDerivAt_aggregates (nodeFaces : N → Finset F) {G : (N → ℝ) → (N → V3) → ℝ}
    {gA : N → ℝ} {gM : N → V3} {a : ℝ → F → ℝ} {m : ℝ → F → V3} {a' : F → ℝ} {m' : F → V3}
    {t : ℝ} (ha : ∀ T, HasDerivAt (fun s => a s T) (a' T) t)
    (hm : ∀ T, HasDerivAt (fun s => m s T) (m' T) t)
    (hG : ∀ (A : ℝ → N → ℝ) (M : ℝ → N → V3) (A' : N → ℝ) (M' : N → V3),
      (∀ C, HasDerivAt (fun s => A s C) (A' C) t) → (∀ C, HasDerivAt (fun s => M s C) (M' C) t) →
        HasDerivAt (fun s => G (A s) (M s)) (∑ C, (gA C * A' C + gM C ⬝ᵥ M' C)) t) :
    HasDerivAt (fun s => G (fun C => ∑ T ∈ nodeFaces C, a s T) (fun C => ∑ T ∈ nodeFaces C, m s T))
      (∑ T, ((∑ C, if T ∈ nodeFaces C then gA C else 0) * a' T +
        (∑ C, if T ∈ nodeFaces C then gM C else 0) ⬝ᵥ m' T)) t := by
  rw [← aggregate_adjoint]
  exact hG _ _ _ _ (fun C => HasDerivAt.fun_sum fun T _ => ha T)
    (fun C => HasDerivAt.fun_sum fun T _ => hm T)

end adjoint

/-! ## The box distance -/

/-- The per-axis gap between boxes `[loA, hiA]` and `[loB, hiB]`. -/
def boxGap (loA hiA loB hiB : V3) (k : Fin 3) : ℝ := max 0 (max (loB k - hiA k) (loA k - hiB k))

/-- **The box distance is a lower bound** (SPEC §3.4): for `x` in box `A` and `y` in box `B`,
`Σ_k gap_k² ≤ |x − y|²`; for overlapping boxes it is `0`, so such a pair is never admissible
(`max(r_I, r_J) < θ · 0` is false).

[formula: rs-box-distance] -/
theorem box_distance (loA hiA loB hiB x y : V3) (hx : ∀ k, loA k ≤ x k ∧ x k ≤ hiA k)
    (hy : ∀ k, loB k ≤ y k ∧ y k ≤ hiB k) :
    ∑ k, boxGap loA hiA loB hiB k ^ 2 ≤ (x - y) ⬝ᵥ (x - y) ∧
      ((∀ k, loB k ≤ hiA k ∧ loA k ≤ hiB k) → ∑ k, boxGap loA hiA loB hiB k ^ 2 = 0) := by
  constructor
  · simp only [dotProduct, Pi.sub_apply]
    refine Finset.sum_le_sum fun k _ => ?_
    obtain ⟨h1, h2⟩ := hx k
    obtain ⟨h3, h4⟩ := hy k
    have hg : 0 ≤ boxGap loA hiA loB hiB k := le_max_left _ _
    have hle : boxGap loA hiA loB hiB k ≤ |x k - y k| := by
      unfold boxGap
      refine max_le (abs_nonneg _) (max_le ?_ ?_)
      · rw [abs_sub_comm]; exact le_trans (by linarith) (le_abs_self _)
      · exact le_trans (by linarith) (le_abs_self _)
    have := mul_le_mul hle hle hg (abs_nonneg _)
    rw [abs_mul_abs_self] at this
    nlinarith
  · intro h
    refine Finset.sum_eq_zero fun k _ => ?_
    have : boxGap loA hiA loB hiB k = 0 := by
      unfold boxGap
      exact max_eq_left (max_le (by linarith [(h k).1]) (by linarith [(h k).2]))
    rw [this]; norm_num

end RS

end DDG
