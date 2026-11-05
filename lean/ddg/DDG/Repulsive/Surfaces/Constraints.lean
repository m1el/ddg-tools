import DDG.Repulsive.Surfaces.Geometry

/-!
# Repulsive Surfaces: constraint values and Jacobians (paper Sec. 7.1; SPEC §5.1, D8, D9)

* `rs-area-gradient`: `A = Σ_S a_S`, `d/ds A = Σ_v (Σ_{S ∋ v} ½ N_S × e_opp) · Γ'_v` along any
  motion with nondegenerate faces.
* `rs-volume-gradient` (D9): `V = ⅙ Σ_S f_i·(f_j × f_k)`, `d/ds V = Σ_v (⅙ Σ_{S ∋ v} f_j × f_k)·Γ'_v`
  (no condition on the faces); on a **closed, consistently oriented fan** around `v`,
  `⅙ Σ_{S∋v} f_j × f_k = ⅓ Σ_{S∋v} a_S N_S`: the reference's row `Σ a_S N_S` is `3∇V`.
* `rs-barycenter` (D8): per component `c`, `Φ_c = Σ_{v∈c} w_v f_v / W_c − X₀` with frozen weights
  `w_v` (dual areas) and `W_c = Σ_{v∈c} w_v`; its derivative is `Σ_{v∈c} (w_v/W_c) Γ'_v` (the
  rows `w_v/W_c ⊗ I₃`); translating component `c` by `X₀ − X_c` makes `Φ_c = 0` exactly and leaves
  every other component's barycenter unchanged.
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Real Finset

variable {V F : Type*} [Fintype V] [Fintype F] [DecidableEq V]

/-- Regrouping a sum over the corners of all faces by vertex. -/
theorem regroup (fv : F → Fin 3 → V) (G : F → Fin 3 → V3) (w : V → V3) :
    ∑ S, ∑ a, G S a ⬝ᵥ w (fv S a) =
      ∑ v, (∑ S, ∑ a, if fv S a = v then G S a else 0) ⬝ᵥ w v := by
  symm
  calc ∑ v, (∑ S, ∑ a, if fv S a = v then G S a else 0) ⬝ᵥ w v
      = ∑ v, ∑ S, ∑ a, (if fv S a = v then G S a ⬝ᵥ w v else 0) := by
        refine Finset.sum_congr rfl fun v _ => ?_
        rw [sum_dotProduct]; refine Finset.sum_congr rfl fun S _ => ?_
        rw [sum_dotProduct]; refine Finset.sum_congr rfl fun a _ => ?_
        split_ifs <;> simp
    _ = ∑ S, ∑ a, ∑ v, (if fv S a = v then G S a ⬝ᵥ w v else 0) := by
        rw [Finset.sum_comm]; refine Finset.sum_congr rfl fun S _ => ?_; rw [Finset.sum_comm]
    _ = _ := by simp

/-! ## Area -/

/-- The total area `A = Σ_S a_S`. -/
def totalArea (f : V → V3) (fv : F → Fin 3 → V) : ℝ := ∑ S, area f fv S

/-- The area gradient `∇_v A = Σ_{S ∋ v} ½ N_S × e_opp(v)` (the Jacobian row of `Area`).

[formula: rs-area-gradient] -/
def areaGrad (f : V → V3) (fv : F → Fin 3 → V) (v : V) : V3 :=
  ∑ S, ∑ a, if fv S a = v then (1 / 2 : ℝ) • (normal f fv S ⨯₃ opp f fv S a) else 0

/-- **The area gradient** (SPEC §5.1): along any differentiable motion with all faces
nondegenerate at `t`, `d/ds A(Γ s) = Σ_v ∇_v A · Γ'_v`.

[formula: rs-area-gradient] -/
theorem hasDerivAt_totalArea {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ} {fv : F → Fin 3 → V}
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (hnd : ∀ S, Nondeg (Γ t) fv S) :
    HasDerivAt (fun s => totalArea (Γ s) fv) (∑ v, areaGrad (Γ t) fv v ⬝ᵥ Γ' v) t := by
  have h := HasDerivAt.fun_sum (u := Finset.univ) fun S _ => hasDerivAt_area hΓ (hnd S)
  refine h.congr_deriv ?_
  rw [regroup fv (fun S a => (1 / 2 : ℝ) • (normal (Γ t) fv S ⨯₃ opp (Γ t) fv S a)) Γ']
  rfl

/-! ## Volume -/

/-- The signed volume `V = ⅙ Σ_S f_i·(f_j × f_k)` (of the cone over the origin for an open
component). -/
def volume (f : V → V3) (fv : F → Fin 3 → V) : ℝ :=
  (1 / 6) * ∑ S, corner f fv S 0 ⬝ᵥ (corner f fv S 1 ⨯₃ corner f fv S 2)

/-- The volume gradient `∇_v V = ⅙ Σ_{S ∋ v} f_{a+1} × f_{a+2}` (the exact Jacobian row, D9).

[formula: rs-volume-gradient] -/
def volGrad (f : V → V3) (fv : F → Fin 3 → V) (v : V) : V3 :=
  ∑ S, ∑ a, if fv S a = v then
    (1 / 6 : ℝ) • (corner f fv S (a + 1) ⨯₃ corner f fv S (a + 2)) else 0

/-- **The volume gradient** (SPEC §5.1, D9): along any differentiable motion,
`d/ds V(Γ s) = Σ_v ∇_v V · Γ'_v`.

[formula: rs-volume-gradient] -/
theorem hasDerivAt_volume {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ} {fv : F → Fin 3 → V}
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) :
    HasDerivAt (fun s => volume (Γ s) fv) (∑ v, volGrad (Γ t) fv v ⬝ᵥ Γ' v) t := by
  have hface : ∀ S, HasDerivAt
      (fun s => corner (Γ s) fv S 0 ⬝ᵥ (corner (Γ s) fv S 1 ⨯₃ corner (Γ s) fv S 2))
      (∑ a, (corner (Γ t) fv S (a + 1) ⨯₃ corner (Γ t) fv S (a + 2)) ⬝ᵥ Γ' (fv S a)) t := by
    intro S
    have h := hasDerivAt_dot (hasDerivAt_corner hΓ fv S 0)
      (hasDerivAt_cross (hasDerivAt_corner hΓ fv S 1) (hasDerivAt_corner hΓ fv S 2))
    refine h.congr_deriv ?_
    simp only [Fin.sum_univ_three, fin3_add_two_zero, fin3_add_one_zero, fin3_add_two_one,
      fin3_add_one_one, fin3_add_two_two, fin3_add_one_two]
    simp [cross_apply, dotProduct, Fin.sum_univ_three]
    ring
  have h := (HasDerivAt.fun_sum (u := Finset.univ) fun S _ => hface S).const_mul (1 / 6 : ℝ)
  refine h.congr_deriv ?_
  have e : (1 / 6 : ℝ) * ∑ S, ∑ a,
      (corner (Γ t) fv S (a + 1) ⨯₃ corner (Γ t) fv S (a + 2)) ⬝ᵥ Γ' (fv S a) =
      ∑ S, ∑ a, ((1 / 6 : ℝ) • (corner (Γ t) fv S (a + 1) ⨯₃ corner (Γ t) fv S (a + 2))) ⬝ᵥ
        Γ' (fv S a) := by
    simp only [smul_dotProduct, smul_eq_mul, Finset.mul_sum]
  rw [e, regroup]
  rfl

omit [Fintype V] [Fintype F] [DecidableEq V] in
/-- Per corner: `f_{a+1} × f_{a+2} = n_S − f_a × (f_{a+1} − f_{a+2})`. -/
theorem cross_next_eq (f : V → V3) (fv : F → Fin 3 → V) (S : F) (a : Fin 3) :
    corner f fv S (a + 1) ⨯₃ corner f fv S (a + 2) =
      nvec f fv S - corner f fv S a ⨯₃ (corner f fv S (a + 1) - corner f fv S (a + 2)) := by
  fin_cases a <;> apply v3_ext <;> simp [nvec, cross_apply] <;> ring

omit [Fintype V] in
/-- **The factor 3** (SPEC D9): if the fan around `v` is closed and consistently oriented in the
sense `Σ_{(S,a) : fv S a = v} (f_{a+1} − f_{a+2}) = 0` (`closedFan_of_equiv`), then
`∇_v V = ⅙ Σ_{S∋v} n_S = ⅓ Σ_{S∋v} a_S N_S`; the reference's row `Σ_{S∋v} a_S N_S` is `3 ∇_v V`.

[formula: rs-volume-gradient] -/
theorem volGrad_closed_fan (f : V → V3) (fv : F → Fin 3 → V) (v : V)
    (hfan : ∑ S, ∑ a, (if fv S a = v then
      corner f fv S (a + 1) - corner f fv S (a + 2) else 0) = 0)
    (hnd : ∀ S, Nondeg f fv S) :
    volGrad f fv v = (1 / 3 : ℝ) • ∑ S, ∑ a, (if fv S a = v then
      area f fv S • normal f fv S else 0) := by
  have h1 : ∀ S a, (if fv S a = v then
      (1 / 6 : ℝ) • (corner f fv S (a + 1) ⨯₃ corner f fv S (a + 2)) else 0) =
      (1 / 6 : ℝ) • (if fv S a = v then nvec f fv S else 0) -
        (1 / 6 : ℝ) • (f v ⨯₃ (if fv S a = v then
          corner f fv S (a + 1) - corner f fv S (a + 2) else 0)) := by
    intro S a
    split_ifs with h
    · rw [cross_next_eq, smul_sub, corner, h]
    · simp
  have h2 : ∀ S a, (if fv S a = v then nvec f fv S else 0) =
      (2 : ℝ) • (if fv S a = v then area f fv S • normal f fv S else 0) := by
    intro S a
    split_ifs
    · rw [(face_quantities (hnd S)).2.1, smul_smul, mul_comm]
    · simp
  simp only [volGrad, h1, Finset.sum_sub_distrib, ← Finset.smul_sum, ← map_sum, hfan, map_zero,
    smul_zero, sub_zero, h2, smul_smul]
  norm_num

omit [Fintype V] in
/-- Corners at a vertex as a sum over a subtype. -/
theorem sum_corners_at (fv : F → Fin 3 → V) (v : V) (g : F → Fin 3 → V3) :
    ∑ S, ∑ a, (if fv S a = v then g S a else 0) =
      ∑ p : {p : F × Fin 3 // fv p.1 p.2 = v}, g p.1.1 p.1.2 := by
  rw [← Fintype.sum_prod_type' (fun S a => if fv S a = v then g S a else 0), ← Finset.sum_filter]
  exact (Finset.sum_subtype _ (by simp) _)

omit [Fintype V] in
/-- **A closed, consistently oriented fan** (SPEC D9): if the corners at `v` can be cycled by a
bijection `σ` that maps each corner to the next face around `v`, whose "next" vertex is this
corner's "previous" vertex (`fv (σ p) (+1) = fv p (+2)`: the faces traverse their common edge in
opposite directions), then `Σ_{(S,a) at v} (f_{a+1} − f_{a+2}) = 0`.

[formula: rs-volume-gradient] -/
theorem closedFan_of_equiv (f : V → V3) (fv : F → Fin 3 → V) (v : V)
    (σ : {p : F × Fin 3 // fv p.1 p.2 = v} ≃ {p : F × Fin 3 // fv p.1 p.2 = v})
    (hσ : ∀ p, fv (σ p).1.1 ((σ p).1.2 + 1) = fv p.1.1 (p.1.2 + 2)) :
    ∑ S, ∑ a, (if fv S a = v then corner f fv S (a + 1) - corner f fv S (a + 2) else 0) = 0 := by
  rw [sum_corners_at, Finset.sum_sub_distrib, sub_eq_zero]
  rw [← Equiv.sum_comp σ (fun p => corner f fv p.1.1 (p.1.2 + 1))]
  exact Finset.sum_congr rfl fun p _ => by simp only [corner, hσ]

/-! ## Barycenter (per component, D8) -/

section bary

variable {K : Type*} [DecidableEq K]

/-- The weight of component `c`: `W_c = Σ_{v ∈ c} w_v`. -/
def compWeight (comp : V → K) (w : V → ℝ) (c : K) : ℝ := ∑ v, if comp v = c then w v else 0

/-- The weighted barycenter of component `c` with frozen weights: `X_c = Σ_{v∈c} w_v f_v / W_c`.

[formula: rs-barycenter] -/
def compBary (comp : V → K) (w : V → ℝ) (f : V → V3) (c : K) : V3 :=
  (compWeight comp w c)⁻¹ • ∑ v, if comp v = c then w v • f v else 0

omit [DecidableEq V] in
/-- **The barycenter Jacobian** (SPEC §5.1, frozen weights): `d/ds X_c(Γ s) =
Σ_{v ∈ c} (w_v/W_c) Γ'_v`, the rows `w_v/W_c ⊗ I₃`.

[formula: rs-barycenter] -/
theorem hasDerivAt_compBary (comp : V → K) (w : V → ℝ) {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ}
    (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) (c : K) :
    HasDerivAt (fun s => compBary comp w (Γ s) c)
      (∑ v, if comp v = c then (w v / compWeight comp w c) • Γ' v else 0) t := by
  have h : ∀ v, HasDerivAt (fun s => if comp v = c then w v • Γ s v else 0)
      (if comp v = c then w v • Γ' v else 0) t := by
    intro v; split_ifs
    · exact (hΓ v).const_smul (w v)
    · exact hasDerivAt_const _ _
  have := (HasDerivAt.fun_sum (u := Finset.univ) fun v _ => h v).const_smul
    (compWeight comp w c)⁻¹
  refine this.congr_deriv ?_
  rw [Finset.smul_sum]
  refine Finset.sum_congr rfl fun v _ => ?_
  split_ifs
  · rw [smul_smul, div_eq_inv_mul]
  · simp

omit [DecidableEq V] in
/-- **The translation projection is exact** (SPEC §5.1, §6.2): translating the vertices of
component `c` by `X₀ − X_c` (with `W_c ≠ 0`) makes `X_c = X₀` and leaves the barycenter of every
other component unchanged.

[formula: rs-barycenter] -/
theorem compBary_translate (comp : V → K) (w : V → ℝ) (f : V → V3) {c : K}
    (hW : compWeight comp w c ≠ 0) (X₀ : V3) :
    let f' := fun v => if comp v = c then f v + (X₀ - compBary comp w f c) else f v
    compBary comp w f' c = X₀ ∧ ∀ c', c' ≠ c → compBary comp w f' c' = compBary comp w f c' := by
  intro f'
  constructor
  · have e : ∑ v, (if comp v = c then w v • f' v else 0) =
        ∑ v, (if comp v = c then w v • f v else 0) + compWeight comp w c • (X₀ - compBary comp w f c) := by
      rw [compWeight, Finset.sum_smul, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun v _ => ?_
      by_cases h : comp v = c <;> simp [f', h, smul_add]
    rw [compBary, e, smul_add, smul_smul, inv_mul_cancel₀ hW, one_smul]
    rw [show (compWeight comp w c)⁻¹ • ∑ v, (if comp v = c then w v • f v else 0) =
      compBary comp w f c from rfl]
    abel
  · intro c' hc'
    unfold compBary
    congr 1
    refine Finset.sum_congr rfl fun v _ => ?_
    by_cases h : comp v = c'
    · have : comp v ≠ c := by rw [h]; exact hc'
      simp [f', h, hc']
    · simp [h]

end bary

end RS

end DDG
