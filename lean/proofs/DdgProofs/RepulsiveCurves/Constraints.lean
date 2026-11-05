import DdgProofs.RepulsiveCurves.Basic

/-!
# Repulsive Curves: the constraint Jacobians (paper Sec. 8.1; SPEC §5.1, D6, D7)

Each Jacobian is stated as a derivative along an arbitrary differentiable motion `Γ` of the
vertices (`Basic.lean`): "`dΦ/ds = Σ_v J_v Γ'_v`" gives the nonzero blocks `J_v` of the row(s).
`s_a = sgn a` (`s₁ = −1`, `s₂ = +1`).

| constraint | `Φ` | blocks |
|---|---|---|
| `EdgeLengths` | `ℓ⁰_I − ℓ_I` | `i₁: +T_Iᵀ`, `i₂: −T_Iᵀ` |
| `TotalLength` | `L⁰ − Σ ℓ_I` | `−Σ_{I∋v} s_{I,v} T_Iᵀ` |
| `Barycenter` (D6) | `(Σ ℓ_I x_I)/L − x₀` | `(1/L) Σ_{I∋v} [(ℓ_I/2) I₃ + s_{I,v}(x_I − x̄) T_Iᵀ]` |
| `Point` | `γ_v − x_v` | `I₃` |
| `Surface` | `f(γ_v)` | `∇f(γ_v)ᵀ` |
| `Tangent` (D7) | `(e_k · s)` | `s_a e_kᵀ P_{T_I}/ℓ_I` (edge site), vertex sites likewise |

and `rc-tangent-frame`: for an orthonormal frame `(X/|X|, e₁, e₂)`, `e₁·s = e₂·s = 0 ⟺ s ∥ X`;
the frame of Duff et al. 2017 is orthonormal; the paper's three rows `T_I − X` have dependent
Jacobian rows (`Σ_k T_k P_T e_k = P_T T = 0`), which is why D7 uses two rows.
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Real Finset

variable {V E : Type*} [Fintype E] {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
  {ep : E → Fin 2 → V}

/-! ## Lengths -/

omit [Fintype E] in
/-- **Jacobian of an edge-length constraint** `Φ = ℓ⁰_I − ℓ_I`: block `+T_I` at `i₁`, `−T_I` at `i₂`.

[formula: rc-constraint-edge-length] -/
theorem jac_edge_length (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {I : E}
    (h : 0 < len (Γ t) ep I) (ℓ0 : ℝ) :
    HasDerivAt (fun s => ℓ0 - len (Γ s) ep I)
      (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I 0) - tangent (Γ t) ep I ⬝ᵥ Γ' (ep I 1)) t := by
  refine ((hasDerivAt_len hΓ h).const_sub ℓ0).congr_deriv ?_
  simp [Fin.sum_univ_two]
  ring

/-- **Jacobian of the total-length constraint** `Φ = L⁰ − Σ_I ℓ_I`: block
`−Σ_{I∋v} s_{I,v} T_Iᵀ`.

[formula: rc-constraint-total-length] -/
theorem jac_total_length (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (h : ∀ I, 0 < len (Γ t) ep I) (L0 : ℝ) :
    HasDerivAt (fun s => L0 - ∑ I, len (Γ s) ep I)
      (∑ I, ∑ a, -(sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a)))) t := by
  have := (HasDerivAt.fun_sum (u := Finset.univ) fun I _ => hasDerivAt_len (ep := ep) hΓ (h I))
  refine (this.const_sub L0).congr_deriv ?_
  simp only [← Finset.sum_neg_distrib]

/-! ## Barycenter (D6) -/

/-- The total length `L = Σ_I ℓ_I`. -/
def totalLen (γ : V → V3) (ep : E → Fin 2 → V) : ℝ := ∑ I, len γ ep I

/-- The length-weighted barycenter `x̄ = (Σ_I ℓ_I x_I)/L` (D6: the paper's `Σ ℓ_I (x_I − x₀)`
divided by `L`, same zero set, violation in length units). -/
def bary (γ : V → V3) (ep : E → Fin 2 → V) : V3 :=
  (totalLen γ ep)⁻¹ • ∑ I, len γ ep I • mid γ ep I

/-- **Jacobian of the barycenter constraint** (D6, exact): `Φ = x̄ − x₀` has
`dΦ/ds = (1/L) Σ_I Σ_a [(ℓ_I/2) Γ'_{i_a} + s_a (T_I·Γ'_{i_a})(x_I − x̄)]`, i.e. the block
`(1/L) Σ_{I∋v} [(ℓ_I/2) I₃ + s_{I,v}(x_I − x̄) T_Iᵀ]` at `v`.

[formula: rc-constraint-barycenter] -/
theorem jac_barycenter (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t)
    (h : ∀ I, 0 < len (Γ t) ep I) (hL : 0 < totalLen (Γ t) ep) (x0 : V3) :
    HasDerivAt (fun s => bary (Γ s) ep - x0)
      ((totalLen (Γ t) ep)⁻¹ • ∑ I, ∑ a,
        ((len (Γ t) ep I / 2) • Γ' (ep I a)
          + (sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) • (mid (Γ t) ep I - bary (Γ t) ep))) t := by
  set L := totalLen (Γ t) ep with hLdef
  have hl : ∀ I, HasDerivAt (fun s => len (Γ s) ep I)
      (∑ a, sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) t :=
    fun I => hasDerivAt_len hΓ (h I)
  have hLd : HasDerivAt (fun s => totalLen (Γ s) ep)
      (∑ I, ∑ a, sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) t :=
    HasDerivAt.fun_sum (u := Finset.univ) fun I _ => hl I
  have hN : HasDerivAt (fun s => ∑ I, len (Γ s) ep I • mid (Γ s) ep I)
      (∑ I, (len (Γ t) ep I • ((1 / 2 : ℝ) • (Γ' (ep I 0) + Γ' (ep I 1)))
        + (∑ a, sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) • mid (Γ t) ep I)) t :=
    HasDerivAt.fun_sum (u := Finset.univ) fun I _ => (hl I).smul (hasDerivAt_mid hΓ)
  have hinv := hLd.inv hL.ne'
  have hb := (hinv.smul hN).sub_const x0
  refine hb.congr_deriv ?_
  -- algebra
  set N := ∑ I, len (Γ t) ep I • mid (Γ t) ep I with hNdef
  set L' := ∑ I, ∑ a, sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a)) with hL'def
  have hbary : bary (Γ t) ep = L⁻¹ • N := rfl
  have hsum : ∑ I, ∑ a, ((len (Γ t) ep I / 2) • Γ' (ep I a)
      + (sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) • (mid (Γ t) ep I - bary (Γ t) ep)) =
      ∑ I, (len (Γ t) ep I • ((1 / 2 : ℝ) • (Γ' (ep I 0) + Γ' (ep I 1)))
        + (∑ a, sgn a * (tangent (Γ t) ep I ⬝ᵥ Γ' (ep I a))) • mid (Γ t) ep I)
      - L' • bary (Γ t) ep := by
    rw [hL'def, Finset.sum_smul, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun I _ => ?_
    simp only [Fin.sum_univ_two, smul_sub, add_smul]
    module
  rw [hsum, hbary, smul_sub, smul_smul]
  have : (-L' / L ^ 2) = -(L⁻¹ * L' * L⁻¹) := by field_simp
  rw [this, neg_smul, ← smul_smul, ← smul_smul]
  simp only [Pi.inv_apply, ← hLdef]
  abel

/-! ## Point and surface -/

omit [Fintype E] in
/-- **Jacobian of a point constraint** `Φ = γ_v − x_v`: the identity block.

[formula: rc-constraint-point] -/
theorem jac_point (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (v : V) (x : V3) :
    HasDerivAt (fun s => Γ s v - x) (Γ' v) t :=
  (hΓ v).sub_const x

omit [Fintype E] in
/-- **Jacobian of a surface constraint** `Φ = f(γ_v)` for a differentiable `f` with gradient `∇f`
(`df(w) = ∇f·w`): the block `∇f(γ_v)ᵀ`.

[formula: rc-constraint-surface] -/
theorem jac_surface (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (v : V) {f : V3 → ℝ}
    {f' : V3 →L[ℝ] ℝ} (hf : HasFDerivAt f f' (Γ t v)) {grad : V3} (hg : ∀ w, f' w = grad ⬝ᵥ w) :
    HasDerivAt (fun s => f (Γ s v)) (grad ⬝ᵥ Γ' v) t := by
  have := hf.comp_hasDerivAt t (hΓ v)
  rw [hg] at this
  exact this

/-! ## Tangent constraints (D7) -/

omit [Fintype E] in
/-- **Jacobian of a tangent row at an edge site** (D7, `s = T_I`): `d(e·T_I)/ds =
Σ_a (s_a/ℓ_I) (P_{T_I} e)·Γ'_{i_a}`, i.e. the blocks `±eᵀ P_{T_I}/ℓ_I`.

[formula: rc-constraint-tangent] -/
theorem jac_tangent_edge (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {I : E}
    (h : 0 < len (Γ t) ep I) (e : V3) :
    HasDerivAt (fun s => e ⬝ᵥ tangent (Γ s) ep I)
      (∑ a, (sgn a / len (Γ t) ep I) * (proj (tangent (Γ t) ep I) e ⬝ᵥ Γ' (ep I a))) t := by
  have hT := hasDerivAt_tangent (ep := ep) (I := I) hΓ h
  refine (hasDerivAt_dot (hasDerivAt_const t e) hT).congr_deriv ?_
  simp only [zero_dotProduct, zero_add, dotProduct_smul, dotProduct_sum, smul_eq_mul,
    ← proj_dot_comm, Finset.mul_sum]
  exact Finset.sum_congr rfl fun a _ => by ring

omit [Fintype E] in
/-- **Jacobian of a tangent row at a vertex site** (D7): for a vertex `v` with neighbours `a` (first
edge) and `b`, `s = (γ_v − γ_a)/|γ_v − γ_a| + (γ_b − γ_v)/|γ_b − γ_v|` and
`d(e·s)/ds = (P_{T₁} e)·(Γ'_v − Γ'_a)/|γ_v − γ_a| + (P_{T₂} e)·(Γ'_b − Γ'_v)/|γ_b − γ_v|`
(`T₁, T₂` the two unit edge directions). At an endpoint (degree 1) only the first term is present.

[formula: rc-constraint-tangent] -/
theorem jac_tangent_vertex (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {v a b : V}
    (ha : Γ t v - Γ t a ≠ 0) (hb : Γ t b - Γ t v ≠ 0) (e : V3) :
    HasDerivAt (fun s => e ⬝ᵥ ((vnorm (Γ s v - Γ s a))⁻¹ • (Γ s v - Γ s a)
        + (vnorm (Γ s b - Γ s v))⁻¹ • (Γ s b - Γ s v)))
      ((vnorm (Γ t v - Γ t a))⁻¹ *
          (proj ((vnorm (Γ t v - Γ t a))⁻¹ • (Γ t v - Γ t a)) e ⬝ᵥ (Γ' v - Γ' a))
        + (vnorm (Γ t b - Γ t v))⁻¹ *
          (proj ((vnorm (Γ t b - Γ t v))⁻¹ • (Γ t b - Γ t v)) e ⬝ᵥ (Γ' b - Γ' v))) t := by
  have h1 := hasDerivAt_unit ((hΓ v).sub (hΓ a)) ha
  have h2 := hasDerivAt_unit ((hΓ b).sub (hΓ v)) hb
  refine (hasDerivAt_dot (hasDerivAt_const t e) (h1.add h2)).congr_deriv ?_
  simp only [Pi.sub_apply, zero_dotProduct, zero_add, dotProduct_add, dotProduct_smul,
    smul_eq_mul, ← proj_dot_comm]

/-- **The paper's three tangent rows are dependent** (D7): the Jacobian rows of `T_I − X` are
`P_{T_I} e_k / ℓ_I` (`k = 1, 2, 3`) on `[−Γ'_{i₁} + Γ'_{i₂}]`, and `Σ_k (T_I)_k P_{T_I} e_k =
P_{T_I} T_I = 0` for a unit `T_I`, so `C Cᵀ` is singular; D7 keeps two rows.

[formula: rc-constraint-tangent] -/
theorem tangent_rows_dependent {T : V3} (hT : T ⬝ᵥ T = 1) :
    ∑ k, T k • proj T (Pi.single k 1) = 0 ∧ proj T T = 0 := by
  refine ⟨?_, proj_self hT⟩
  have hdot : ∀ k, T ⬝ᵥ Pi.single k (1 : ℝ) = T k := fun k => by simp [dotProduct, Pi.single_apply]
  simp only [proj, hdot, smul_sub, smul_smul, Finset.sum_sub_distrib, ← Finset.sum_smul]
  have h1 : ∑ k, T k • (Pi.single k (1 : ℝ) : V3) = T := by
    ext i; simp [Finset.sum_apply, Pi.single_apply]
  have h2 : ∑ k, T k * T k = 1 := by rw [← hT]; rfl
  rw [h1, h2, one_smul, sub_self]

/-! ## The tangent frame (D7) -/

/-- **`e₁·s = e₂·s = 0 ⟺ s ∥ u`** for an orthonormal frame `(u, e₁, e₂)` of `ℝ³` (`u = X/|X|`), and
then `s = (s·u) u`.

[formula: rc-tangent-frame] -/
theorem tangent_frame {u e1 e2 : V3} (hu : u ⬝ᵥ u = 1) (h1 : e1 ⬝ᵥ e1 = 1) (h2 : e2 ⬝ᵥ e2 = 1)
    (hu1 : u ⬝ᵥ e1 = 0) (hu2 : u ⬝ᵥ e2 = 0) (h12 : e1 ⬝ᵥ e2 = 0) (s : V3) :
    (e1 ⬝ᵥ s = 0 ∧ e2 ⬝ᵥ s = 0) ↔ s = (u ⬝ᵥ s) • u := by
  -- the matrix with rows `u, e₁, e₂` is orthogonal
  set Q : Matrix (Fin 3) (Fin 3) ℝ := Matrix.of ![u, e1, e2] with hQ
  have hQQ : Q * Qᵀ = 1 := by
    ext i j
    have : (Q * Qᵀ) i j = ![u, e1, e2] i ⬝ᵥ ![u, e1, e2] j := rfl
    rw [this]
    fin_cases i <;> fin_cases j <;>
      simp [hu, h1, h2, hu1, hu2, h12, dotProduct_comm e1 u, dotProduct_comm e2 u,
        dotProduct_comm e2 e1]
  have hQQ' : Qᵀ * Q = 1 := mul_eq_one_comm.mp hQQ
  have hs : s = (u ⬝ᵥ s) • u + (e1 ⬝ᵥ s) • e1 + (e2 ⬝ᵥ s) • e2 := by
    have : s = Qᵀ *ᵥ (Q *ᵥ s) := by rw [mulVec_mulVec, hQQ', one_mulVec]
    conv_lhs => rw [this]
    ext i
    simp [hQ, mulVec, dotProduct, Fin.sum_univ_three, transpose_apply]
    ring
  constructor
  · rintro ⟨ha, hb⟩
    rw [ha, hb, zero_smul, zero_smul, add_zero, add_zero] at hs
    exact hs
  · intro h
    rw [h, dotProduct_smul, dotProduct_smul, dotProduct_comm e1 u, dotProduct_comm e2 u, hu1, hu2]
    simp

/-- The first frame vector of Duff et al. 2017 (`sign = ±1`, `a = −1/(sign + n_z)`,
`b = n_x n_y a`): `b₁ = (1 + sign·n_x² a, sign·b, −sign·n_x)`. -/
def duffB1 (sgnz : ℝ) (n : V3) : V3 :=
  ![1 + sgnz * n 0 ^ 2 * (-1 / (sgnz + n 2)), sgnz * (n 0 * n 1 * (-1 / (sgnz + n 2))),
    -(sgnz * n 0)]

/-- The second frame vector of Duff et al. 2017: `b₂ = (b, sign + n_y² a, −n_y)`. -/
def duffB2 (sgnz : ℝ) (n : V3) : V3 :=
  ![n 0 * n 1 * (-1 / (sgnz + n 2)), sgnz + n 1 ^ 2 * (-1 / (sgnz + n 2)), -n 1]

/-- **The deterministic frame of Duff et al. 2017 is orthonormal**: for a unit `n`, `sign² = 1` and
`sign + n_z ≠ 0` (`sign = copysign(1, n_z)`), `(n, b₁, b₂)` is orthonormal.

[formula: rc-tangent-frame] -/
theorem duff_frame {n : V3} (hn : n ⬝ᵥ n = 1) {sgnz : ℝ} (hs : sgnz ^ 2 = 1)
    (hz : sgnz + n 2 ≠ 0) :
    duffB1 sgnz n ⬝ᵥ duffB1 sgnz n = 1 ∧ duffB2 sgnz n ⬝ᵥ duffB2 sgnz n = 1 ∧
      duffB1 sgnz n ⬝ᵥ duffB2 sgnz n = 0 ∧ n ⬝ᵥ duffB1 sgnz n = 0 ∧
      n ⬝ᵥ duffB2 sgnz n = 0 := by
  have hn' : n 0 ^ 2 + n 1 ^ 2 + n 2 ^ 2 = 1 := by
    rw [← hn]; simp [dotProduct, Fin.sum_univ_three]; ring
  simp only [duffB1, duffB2, dotProduct, Fin.sum_univ_three, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · field_simp
    linear_combination (n 0 ^ 2 * sgnz ^ 2) * hn' + (n 0 ^ 2 * sgnz ^ 2 + 2 * sgnz * n 2 * n 0 ^ 2) * hs
  · field_simp
    linear_combination n 1 ^ 2 * hn' + (2 * sgnz * n 2 + n 2 ^ 2 + sgnz ^ 2 - n 1 ^ 2) * hs
  · field_simp
    linear_combination (n 0 * n 1 * sgnz) * hn' + (n 0 * n 1 * n 2) * hs
  · field_simp
    linear_combination (-n 0 * sgnz) * hn' + (-n 0 * n 2) * hs
  · field_simp
    linear_combination (-n 1) * hn' + n 1 * hs

end RC

end DDG
