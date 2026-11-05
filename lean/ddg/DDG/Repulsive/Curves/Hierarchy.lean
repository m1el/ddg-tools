import DDG.Repulsive.Curves.Gram

/-!
# Repulsive Curves: the lumped quantities of the hierarchy (paper B.1, B.2, Eq. 26, 28;
SPEC §8.1–§8.3, D15, D16)

Edges carry masses `ℓ_I`, midpoints `x_I` and tangents `T_I`; a node is a finite set `N` of edges.

* `rc-bvh-lumped`: `L_N = Σ ℓ_I`, `x̄_N = Σ ℓ_I x_I / L_N`, `T̂_N = Σ ℓ_I T_I / |Σ ℓ_I T_I|` (D15:
  normalized; it is the normalization of the paper's `T̄_N = Σ ℓ_I T_I / L_N`); the data of a parent
  is the mass-weighted combination of its children's; `x̄_N` lies in the box of the midpoints; the
  endpoint-box guard excludes every edge adjacent to `I`.
* Eq. 26: the lumped term `k(x_I, x̄_N, T_I) ℓ_I L_N` equals the sum over the node of the midpoint
  terms `k(x_I, x_J, T_I) ℓ_I ℓ_J` when all midpoints of the node coincide (exact as the cluster
  shrinks; the error bound for `r_x/d, r_T` small is cited, `rc-bh-admissibility`).
* `rc-bh-differential`: the B.1.2 increment is the derivative of `ℓ_I L_N (k(x_I, x̄_N, T_I) +
  k(x̄_N, x_I, T̂_N))` with `L_N, x̄_N, T̂_N` fixed.
* `rc-bct-lowrank`: the rank-one block `ℓ[𝒜] κ ℓ[ℬ]ᵀ` applied to `ψ` is `ℓ_I κ Σ_{J∈ℬ} ℓ_J ψ_J`
  (the 0th-order multipole pass); it equals the exact kernel matrix (midpoint weights, Eq. 32)
  when both clusters are collapsed to their lumped tangent-points.
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Real Finset

/-! ## Node data -/

section node

variable {E : Type*} [DecidableEq E] (ℓ : E → ℝ) (x T : E → V3)

/-- `L_N = Σ_{I∈N} ℓ_I`. -/
def nodeMass (N : Finset E) : ℝ := ∑ I ∈ N, ℓ I

/-- `Σ_{I∈N} ℓ_I x_I`. -/
def nodeMoment (N : Finset E) : V3 := ∑ I ∈ N, ℓ I • x I

/-- `x̄_N = Σ ℓ_I x_I / L_N`.

[formula: rc-bvh-lumped] -/
def nodeCentre (N : Finset E) : V3 := (nodeMass ℓ N)⁻¹ • nodeMoment ℓ x N

/-- `Σ_{I∈N} ℓ_I T_I`. -/
def nodeTSum (N : Finset E) : V3 := ∑ I ∈ N, ℓ I • T I

/-- `T̂_N = Σ ℓ_I T_I / |Σ ℓ_I T_I|` (normalized, D15).

[formula: rc-bvh-lumped] -/
def nodeTangent (N : Finset E) : V3 := (vnorm (nodeTSum ℓ T N))⁻¹ • nodeTSum ℓ T N

/-- **A parent's data is the mass-weighted combination of its children's** (two disjoint
children): masses, moments and tangent sums add; `x̄ = (L₁x̄₁ + L₂x̄₂)/(L₁ + L₂)`;
`T̂ = (S₁ + S₂)/|S₁ + S₂|`.

[formula: rc-bvh-lumped] -/
theorem node_union {N1 N2 : Finset E} (h : Disjoint N1 N2) :
    nodeMass ℓ (N1 ∪ N2) = nodeMass ℓ N1 + nodeMass ℓ N2 ∧
      nodeMoment ℓ x (N1 ∪ N2) = nodeMoment ℓ x N1 + nodeMoment ℓ x N2 ∧
      nodeTSum ℓ T (N1 ∪ N2) = nodeTSum ℓ T N1 + nodeTSum ℓ T N2 ∧
      (nodeMass ℓ N1 ≠ 0 → nodeMass ℓ N2 ≠ 0 →
        nodeCentre ℓ x (N1 ∪ N2) = (nodeMass ℓ N1 + nodeMass ℓ N2)⁻¹ •
          (nodeMass ℓ N1 • nodeCentre ℓ x N1 + nodeMass ℓ N2 • nodeCentre ℓ x N2)) ∧
      nodeTangent ℓ T (N1 ∪ N2) =
        (vnorm (nodeTSum ℓ T N1 + nodeTSum ℓ T N2))⁻¹ • (nodeTSum ℓ T N1 + nodeTSum ℓ T N2) := by
  have hm : nodeMass ℓ (N1 ∪ N2) = nodeMass ℓ N1 + nodeMass ℓ N2 := Finset.sum_union h
  have hx : nodeMoment ℓ x (N1 ∪ N2) = nodeMoment ℓ x N1 + nodeMoment ℓ x N2 := Finset.sum_union h
  have ht : nodeTSum ℓ T (N1 ∪ N2) = nodeTSum ℓ T N1 + nodeTSum ℓ T N2 := Finset.sum_union h
  refine ⟨hm, hx, ht, fun h1 h2 => ?_, by rw [nodeTangent, ht]⟩
  rw [nodeCentre, hm, hx, nodeCentre, nodeCentre, smul_smul, smul_smul, mul_inv_cancel₀ h1,
    mul_inv_cancel₀ h2, one_smul, one_smul]

omit [DecidableEq E] in
/-- **`T̂_N` is the normalization of the paper's `T̄_N = Σ ℓ_I T_I / L_N`** (D15), for `L_N > 0`.

[formula: rc-bvh-lumped] -/
theorem nodeTangent_eq_normalize_paper {N : Finset E} (hm : 0 < nodeMass ℓ N) :
    nodeTangent ℓ T N =
      (vnorm ((nodeMass ℓ N)⁻¹ • nodeTSum ℓ T N))⁻¹ • ((nodeMass ℓ N)⁻¹ • nodeTSum ℓ T N) := by
  rw [nodeTangent, vnorm_smul, abs_of_pos (inv_pos.mpr hm), smul_smul]
  congr 1
  by_cases h : vnorm (nodeTSum ℓ T N) = 0
  · simp [h]
  · field_simp

omit [DecidableEq E] in
/-- **The centre lies in the box of the midpoints**: for positive masses and a nonempty node, every
coordinate of `x̄_N` lies between the least and the largest coordinate of the midpoints.

[formula: rc-bvh-lumped] -/
theorem nodeCentre_mem_box {N : Finset E} (hpos : ∀ I ∈ N, 0 < ℓ I) (hN : N.Nonempty) (c : Fin 3)
    {lo hi : ℝ} (hlo : ∀ I ∈ N, lo ≤ x I c) (hhi : ∀ I ∈ N, x I c ≤ hi) :
    lo ≤ nodeCentre ℓ x N c ∧ nodeCentre ℓ x N c ≤ hi := by
  have hm : 0 < nodeMass ℓ N := Finset.sum_pos hpos hN
  have hc : nodeCentre ℓ x N c = (∑ I ∈ N, ℓ I * x I c) / nodeMass ℓ N := by
    simp [nodeCentre, nodeMoment, Finset.sum_apply, div_eq_inv_mul]
  rw [hc, le_div_iff₀ hm, div_le_iff₀ hm]
  constructor
  · rw [nodeMass, Finset.mul_sum]
    exact Finset.sum_le_sum fun I hI => by nlinarith [hlo I hI, hpos I hI]
  · rw [nodeMass, Finset.mul_sum]
    exact Finset.sum_le_sum fun I hI => by nlinarith [hhi I hI, hpos I hI]

omit [DecidableEq E] in
/-- **The endpoint-box guard excludes the adjacent edges exactly** (D15, D16): if every endpoint
of every edge of the node lies in a set `B` and neither endpoint of `I` does, then `I` has no
common vertex with any edge of the node.

[formula: rc-bvh-lumped] -/
theorem guard_disjoint {V : Type*} (γ : V → V3) (ep : E → Fin 2 → V) {N : Finset E} {B : Set V3}
    (hB : ∀ J ∈ N, ∀ b, γ (ep J b) ∈ B) {I : E} (hI : ∀ a, γ (ep I a) ∉ B) :
    ∀ J ∈ N, ∀ a b, ep I a ≠ ep J b := by
  intro J hJ a b h
  exact hI a (h ▸ hB J hJ b)

end node

/-! ## Eq. 26: the lumped energy -/

section lumped

variable {E : Type*} [DecidableEq E] (ℓ : E → ℝ) (x T : E → V3)

/-- The lumped term of Eq. 26, `k(x_I, x̄_N, T_I) ℓ_I L_N`.

[formula: rc-bvh-lumped] -/
def lumpedTerm (α β : ℝ) (xI TI : V3) (ℓI : ℝ) (N : Finset E) : ℝ :=
  kernel α β xI (nodeCentre ℓ x N) TI * ℓI * nodeMass ℓ N

omit [DecidableEq E] in
/-- **Eq. 26 is exact for a collapsed cluster**: if all midpoints of the node are `c` and its mass
is positive, then `x̄_N = c` and `Σ_{J∈N} k(x_I, x_J, T_I) ℓ_I ℓ_J = k(x_I, x̄_N, T_I) ℓ_I L_N`.

[formula: rc-bvh-lumped] -/
theorem lumped_exact (α β : ℝ) (xI TI : V3) (ℓI : ℝ) {N : Finset E} {c : V3}
    (hc : ∀ J ∈ N, x J = c) (hm : nodeMass ℓ N ≠ 0) :
    nodeCentre ℓ x N = c ∧
      ∑ J ∈ N, kernel α β xI (x J) TI * ℓI * ℓ J = lumpedTerm ℓ x α β xI TI ℓI N := by
  have hmom : nodeMoment ℓ x N = nodeMass ℓ N • c := by
    rw [nodeMoment, nodeMass, Finset.sum_smul]
    exact Finset.sum_congr rfl fun J hJ => by rw [hc J hJ]
  have hcen : nodeCentre ℓ x N = c := by
    rw [nodeCentre, hmom, smul_smul, inv_mul_cancel₀ hm, one_smul]
  refine ⟨hcen, ?_⟩
  rw [lumpedTerm, hcen, nodeMass, Finset.mul_sum]
  exact Finset.sum_congr rfl fun J hJ => by rw [hc J hJ]

end lumped

/-! ## The Barnes–Hut differential (B.1.2) -/

section bhDifferential

variable {V E : Type*} {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ} {ep : E → Fin 2 → V}

/-- **The B.1.2 increment is the derivative of the lumped term** with `L_N, x̄_N, T̂_N` held fixed:
along a motion with `ν > 0` in both kernel evaluations and a unit `T̂_N`,
`d/ds [ℓ_I L_N (k(x_I, x̄, T_I) + k(x̄, x_I, T̂))] = Σ_a g_a · Γ'_{i_a}` with
`g_a = L_N (s_a (k₁ + k₂) T_I + s_a ∇_T k₁ + (ℓ_I/2)(∇_p k₁ − ∇_p k₂))`
(`∇_q k₂ = −∇_p k₂` is the derivative of `k₂` in its second point).

[formula: rc-bh-differential] -/
theorem bh_differential {α β : ℝ} (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {I : E}
    (hI : 0 < len (Γ t) ep I) (L : ℝ) (c Th : V3) (hTh : Th ⬝ᵥ Th = 1)
    (h1 : crossProduct (tangent (Γ t) ep I) (mid (Γ t) ep I - c) ≠ 0)
    (h2 : crossProduct Th (c - mid (Γ t) ep I) ≠ 0) :
    HasDerivAt (fun s => len (Γ s) ep I * L *
        (kernel α β (mid (Γ s) ep I) c (tangent (Γ s) ep I) + kernel α β c (mid (Γ s) ep I) Th))
      (∑ a, (L • ((sgn a * (kernel α β (mid (Γ t) ep I) c (tangent (Γ t) ep I)
            + kernel α β c (mid (Γ t) ep I) Th)) • tangent (Γ t) ep I
          + sgn a • kgradT α β (mid (Γ t) ep I) c (tangent (Γ t) ep I)
          + (len (Γ t) ep I / 2) • (kgradP α β (mid (Γ t) ep I) c (tangent (Γ t) ep I)
              - kgradP α β c (mid (Γ t) ep I) Th))) ⬝ᵥ Γ' (ep I a)) t := by
  have hT := hasDerivAt_tangent (ep := ep) (I := I) hΓ hI
  have hx := hasDerivAt_mid (ep := ep) (I := I) hΓ
  have hunit := tangent_dot_self hI
  have hk1 := hasDerivAt_kernel (α := α) (β := β) hx (hasDerivAt_const t c) hT hunit
    (tangent_dot_deriv hI) h1
  have hk2 := hasDerivAt_kernel (α := α) (β := β) (hasDerivAt_const t c) hx (hasDerivAt_const t Th)
    hTh (by simp) h2
  have hl := hasDerivAt_len (ep := ep) (I := I) hΓ hI
  refine ((hl.mul_const L).mul (hk1.add hk2)).congr_deriv ?_
  have hl0 : len (Γ t) ep I ≠ 0 := hI.ne'
  simp only [dotProduct_smul, dotProduct_sum, kgradT_dot_proj hunit, smul_eq_mul, sub_zero,
    zero_sub, dotProduct_zero, add_zero, dotProduct_neg, Pi.add_apply]
  simp only [add_dotProduct, sub_dotProduct, smul_dotProduct, smul_eq_mul, Fin.sum_univ_two,
    sgn_zero, sgn_one, dotProduct_add, dotProduct_comm (tangent (Γ t) ep I)]
  field_simp
  ring

end bhDifferential

/-! ## The block cluster tree (B.2, Eq. 28) -/

section bct

variable {E : Type*} (σ : ℝ) (ℓ : E → ℝ) (x T : E → V3)

/-- The exact low-order kernel matrix with midpoint weights (Eq. 32),
`K⁰_IJ = (k²_{2σ+5}(x_I, x_J, T_I) + k²_{2σ+5}(x_J, x_I, T_J)) ℓ_I ℓ_J`. -/
def K0mid (I J : E) : ℝ :=
  (kernel 2 (2 * σ + 5) (x I) (x J) (T I) + kernel 2 (2 * σ + 5) (x J) (x I) (T J)) * ℓ I * ℓ J

/-- The exact high-order kernel matrix with midpoint weights, `K^B_IJ = 2 ℓ_I ℓ_J |x_I − x_J|^{−(2σ+1)}`. -/
def KBmid (I J : E) : ℝ := 2 * ℓ I * ℓ J * vnorm (x I - x J) ^ (-(2 * σ + 1))

/-- **The admissible block is the lumped kernel; it is exact for collapsed clusters**
(B.2.2, SPEC §8.3): applying `ℓ[𝒜] κ ℓ[ℬ]ᵀ` to `ψ` gives `ℓ_I κ Σ_{J∈ℬ} ℓ_J ψ_J` (one upward sum
per node); and if every edge of `𝒜` has the tangent-point `(T̂_𝒜, x̄_𝒜)` and every edge of `ℬ`
`(T̂_ℬ, x̄_ℬ)`, then `K⁰_IJ = ℓ_I κ⁰ ℓ_J` and `K^B_IJ = ℓ_I κ^B ℓ_J` with
`κ⁰ = k²_{2σ+5}(x̄_𝒜, x̄_ℬ, T̂_𝒜) + k²_{2σ+5}(x̄_ℬ, x̄_𝒜, T̂_ℬ)`, `κ^B = 2|x̄_𝒜 − x̄_ℬ|^{−(2σ+1)}`.

[formula: rc-bct-lowrank] -/
theorem bct_lowrank (A B : Finset E) (κ : ℝ) (ψ : E → ℝ) (I : E) {xA xB TA TB : V3}
    (hA : ∀ I ∈ A, x I = xA ∧ T I = TA) (hB : ∀ J ∈ B, x J = xB ∧ T J = TB) :
    ∑ J ∈ B, ℓ I * κ * ℓ J * ψ J = ℓ I * κ * ∑ J ∈ B, ℓ J * ψ J ∧
      ∀ I ∈ A, ∀ J ∈ B,
        K0mid σ ℓ x T I J = ℓ I * (kernel 2 (2 * σ + 5) xA xB TA
          + kernel 2 (2 * σ + 5) xB xA TB) * ℓ J ∧
        KBmid σ ℓ x I J = ℓ I * (2 * vnorm (xA - xB) ^ (-(2 * σ + 1))) * ℓ J := by
  refine ⟨by rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun J _ => by ring,
    fun I hI J hJ => ?_⟩
  obtain ⟨hx, hT⟩ := hA I hI
  obtain ⟨hx', hT'⟩ := hB J hJ
  refine ⟨?_, ?_⟩
  · rw [K0mid, hx, hT, hx', hT']; ring
  · rw [KBmid, hx, hx']; ring

/-- The kernel matrices of the hierarchy are those of the dense metric: with `x, T, ℓ` the
midpoints, tangents and lengths of a curve, `K⁰mid = w⁰ + w⁰ᵀ` and `K^Bmid = w + wᵀ` for the
midpoint weights of `Gram.lean`.

[formula: rc-bct-lowrank] -/
theorem K_mid_eq_gram {V : Type*} (γ : V → V3) (ep : E → Fin 2 → V) (I J : E) :
    K0mid σ (len γ ep) (mid γ ep) (tangent γ ep) I J = w0Mid σ γ ep I J + w0Mid σ γ ep J I ∧
      KBmid σ (len γ ep) (mid γ ep) I J = wMid σ γ ep I J + wMid σ γ ep J I := by
  obtain ⟨h0, -, hB⟩ := eq32_midpoint σ γ ep I J
  exact ⟨by rw [h0, K0mid], by rw [hB, KBmid]⟩

end bct

end RC

end DDG
