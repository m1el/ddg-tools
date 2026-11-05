import DDG.Repulsive.Surfaces.Kernel
import DDG.Repulsive.Curves.Gram

/-!
# Repulsive Surfaces: the operators of the inner product (paper Sec. 3.2, Eq. 11–13, App. A;
SPEC §4.1–§4.3, §8.3)

* `rs-df-gradient` (SPEC §4.1, D2): `D_f u(S) = (1/(2a_S)) N_S × (u_i e_jk + u_j e_ki + u_k e_ij)`
  with **un-normalized** edges is the P1 gradient: `D_f u ⊥ N_S`,
  `D_f u · (f_b − f_a) = u_b − u_a` on the triangle, the unique such vector, and `D_f 1 = 0`. The
  paper's "unit edge vectors" give a different vector: an explicit triangle where the unit-edge
  form is half the gradient (`df_unit_edge_counterexample`).
* `rs-sobolev-forms` (Eq. 11–13): the matrix assembled by increments over ordered pairs,
  `M = Σ_{S,T} w_ST Σ_c (r_{S,c} − r_{T,c})(r_{S,c} − r_{T,c})ᵀ` for feature rows `r` (the
  averaging rows of `ū`, or the three component rows of `D_f`), has the bilinear form
  `Σ_{S,T} w_ST Σ_c ((Ru)_{S,c} − (Ru)_{T,c})((Rv)_{S,c} − (Rv)_{T,c})`; it is symmetric and
  positive semidefinite for `w ≥ 0`, whatever the symmetry of `w`; constants are in its kernel
  when `R 1 = 0` per pair difference (`ū(1) = 1`, `D_f 1 = 0`): `L_τ 1 = B 1 = B⁰ 1 = 0`.
* `rs-low-order-symmetrisation`: an ordered sum with a symmetric factor `g` equals the sum with the
  symmetrised kernel `½(K_ST + K_TS)`.
* `rs-kernel-matrix-form` (App. A, SPEC §4.3): for symmetric `H` and `a_S ≠ 0`,
  `Σ_{S,T} (x_S − x_T)(y_S − y_T) a_S H_ST a_T = 2 (a⊙x)ᵀ (diag(a)⁻¹ diag(H a) − H) (a⊙y)`, with
  `a⊙ū = U u`, and the vector version (`V u = a ⊙ D_f u`, summed over the components); as a matrix
  identity `L = 2 Uᵀ (diag(a)⁻¹ diag(H a) − H) U`.
* `rs-hierarchical-null-space`: `(diag(a)⁻¹ diag(H a) − H) a = 0` for **any** `H` (also a
  compressed, non-symmetric one) when `diag(H a)` uses the same `H`; hence
  `Uᵀ(diag(a)⁻¹ diag(H a) − H) U 1 = 0`.
* `rs-projector-average` (SPEC §8.3): `P̄ = Σ_S w_S N_S N_Sᵀ` has `vᵀP̄v = Σ_S w_S |P_S v|²`
  (unit normals), and `0 ≤ vᵀP̄v ≤ |v|²` for a convex combination.
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Real Finset
open RC (cross_vnorm_sq posSemidef_of_form)

/-! ## The P1 gradient `D_f` -/

section df

variable {V F : Type*}

/-- The per-corner vector `(1/(2a_S)) N_S × e_a`: `D_f u(S) = Σ_a u_a d_a`. -/
def dvec (f : V → V3) (fv : F → Fin 3 → V) (S : F) (a : Fin 3) : V3 :=
  (1 / (2 * area f fv S)) • (normal f fv S ⨯₃ opp f fv S a)

/-- **The P1 gradient** `D_f u(S) = (1/(2a_S)) N_S × (u_i e_jk + u_j e_ki + u_k e_ij)` (SPEC §4.1,
un-normalized edges, D2).

[formula: rs-df-gradient] -/
def dfGrad (f : V → V3) (fv : F → Fin 3 → V) (u : V → ℝ) (S : F) : V3 :=
  (1 / (2 * area f fv S)) • (normal f fv S ⨯₃ ∑ a, u (fv S a) • opp f fv S a)

theorem dfGrad_eq_sum (f : V → V3) (fv : F → Fin 3 → V) (u : V → ℝ) (S : F) :
    dfGrad f fv u S = ∑ a, u (fv S a) • dvec f fv S a := by
  simp only [dfGrad, dvec, map_sum, map_smul, Finset.smul_sum, smul_comm (u _)]

/-- `D_f` is linear in `u`. -/
theorem dfGrad_add (f : V → V3) (fv : F → Fin 3 → V) (u w : V → ℝ) (S : F) :
    dfGrad f fv (u + w) S = dfGrad f fv u S + dfGrad f fv w S := by
  simp only [dfGrad_eq_sum, Pi.add_apply, add_smul, Finset.sum_add_distrib]

/-- **`D_f 1 = 0`** (and `D_f c = 0` for every constant).

[formula: rs-df-gradient] -/
theorem dfGrad_const (f : V → V3) (fv : F → Fin 3 → V) (c : ℝ) (S : F) :
    dfGrad f fv (fun _ => c) S = 0 := by
  simp only [dfGrad, ← Finset.smul_sum, sum_opp, smul_zero, map_zero]

/-- The polynomial identity behind the edge property: `(n × w)·(p₁ − p₀) = (u₁ − u₀)|n|²`. -/
theorem cross_w_dot (p0 p1 p2 : V3) (u0 u1 u2 : ℝ) :
    (((p1 - p0) ⨯₃ (p2 - p0)) ⨯₃ (u0 • (p2 - p1) + u1 • (p0 - p2) + u2 • (p1 - p0))) ⬝ᵥ (p1 - p0)
      = (u1 - u0) * (((p1 - p0) ⨯₃ (p2 - p0)) ⬝ᵥ ((p1 - p0) ⨯₃ (p2 - p0))) := by
  simp [cross_apply, dotProduct, Fin.sum_univ_three]
  ring

theorem cross_w_dot' (p0 p1 p2 : V3) (u0 u1 u2 : ℝ) :
    (((p1 - p0) ⨯₃ (p2 - p0)) ⨯₃ (u0 • (p2 - p1) + u1 • (p0 - p2) + u2 • (p1 - p0))) ⬝ᵥ (p2 - p0)
      = (u2 - u0) * (((p1 - p0) ⨯₃ (p2 - p0)) ⬝ᵥ ((p1 - p0) ⨯₃ (p2 - p0))) := by
  simp [cross_apply, dotProduct, Fin.sum_univ_three]
  ring

/-- `D_f u = |n_S|⁻² (n_S × w)` with `w = Σ_a u_a e_a`. -/
theorem dfGrad_eq_nvec {f : V → V3} {fv : F → Fin 3 → V} {S : F} (u : V → ℝ) :
    dfGrad f fv u S = ((vnorm (nvec f fv S))⁻¹ * (vnorm (nvec f fv S))⁻¹) •
      (nvec f fv S ⨯₃ ∑ a, u (fv S a) • opp f fv S a) := by
  rw [dfGrad, normal, map_smul, LinearMap.smul_apply, smul_smul, area]
  congr 1
  ring

/-- **`D_f u` is the P1 gradient** (SPEC §4.1): for a nondegenerate face, `D_f u ⊥ N_S` and
`D_f u · (f_b − f_a) = u_b − u_a` for every two corners `a`, `b`.

[formula: rs-df-gradient] -/
theorem dfGrad_props {f : V → V3} {fv : F → Fin 3 → V} {S : F} (hS : Nondeg f fv S)
    (u : V → ℝ) :
    normal f fv S ⬝ᵥ dfGrad f fv u S = 0 ∧
      ∀ a b, dfGrad f fv u S ⬝ᵥ (corner f fv S b - corner f fv S a) = u (fv S b) - u (fv S a) := by
  have hpos : 0 < vnorm (nvec f fv S) := vnorm_pos_iff.mpr hS
  have hnn : nvec f fv S ⬝ᵥ nvec f fv S = vnorm (nvec f fv S) * vnorm (nvec f fv S) :=
    (vnorm_mul_self _).symm
  have hw : ∑ a, u (fv S a) • opp f fv S a =
      u (fv S 0) • (corner f fv S 2 - corner f fv S 1) + u (fv S 1) • (corner f fv S 0 -
        corner f fv S 2) + u (fv S 2) • (corner f fv S 1 - corner f fv S 0) := by
    simp [Fin.sum_univ_three, opp]
  have e01 : dfGrad f fv u S ⬝ᵥ (corner f fv S 1 - corner f fv S 0) = u (fv S 1) - u (fv S 0) := by
    rw [dfGrad_eq_nvec, smul_dotProduct, hw, nvec, cross_w_dot, ← nvec, hnn, smul_eq_mul]
    field_simp
  have e02 : dfGrad f fv u S ⬝ᵥ (corner f fv S 2 - corner f fv S 0) = u (fv S 2) - u (fv S 0) := by
    rw [dfGrad_eq_nvec, smul_dotProduct, hw, nvec, cross_w_dot', ← nvec, hnn, smul_eq_mul]
    field_simp
  refine ⟨?_, fun a b => ?_⟩
  · rw [dfGrad, dotProduct_smul, dot_self_cross, smul_zero]
  · have hb : ∀ b : Fin 3, dfGrad f fv u S ⬝ᵥ (corner f fv S b - corner f fv S 0) =
        u (fv S b) - u (fv S 0) := by
      intro b; fin_cases b
      · simp
      · exact e01
      · exact e02
    have : corner f fv S b - corner f fv S a =
        (corner f fv S b - corner f fv S 0) - (corner f fv S a - corner f fv S 0) := by abel
    rw [this, dotProduct_sub, hb, hb]
    ring

/-- **Uniqueness**: a vector orthogonal to `N_S` with the P1 edge differences on two edges of a
nondegenerate face is `D_f u`.

[formula: rs-df-gradient] -/
theorem dfGrad_unique {f : V → V3} {fv : F → Fin 3 → V} {S : F} (hS : Nondeg f fv S)
    (u : V → ℝ) {g : V3} (hN : normal f fv S ⬝ᵥ g = 0)
    (h1 : g ⬝ᵥ (corner f fv S 1 - corner f fv S 0) = u (fv S 1) - u (fv S 0))
    (h2 : g ⬝ᵥ (corner f fv S 2 - corner f fv S 0) = u (fv S 2) - u (fv S 0)) :
    g = dfGrad f fv u S := by
  obtain ⟨hN', hE⟩ := dfGrad_props hS u
  set d := g - dfGrad f fv u S with hd
  have hd1 : d ⬝ᵥ (corner f fv S 1 - corner f fv S 0) = 0 := by
    rw [hd, sub_dotProduct, h1, hE]; ring
  have hd2 : d ⬝ᵥ (corner f fv S 2 - corner f fv S 0) = 0 := by
    rw [hd, sub_dotProduct, h2, hE]; ring
  have hdn : nvec f fv S ⬝ᵥ d = 0 := by
    obtain ⟨-, hn, -⟩ := face_quantities hS
    rw [hn, smul_dotProduct, hd, dotProduct_sub, hN, hN', sub_zero, smul_zero]
  -- `d × n = e₁ (d·e₂) − e₂ (d·e₁) = 0`
  have hx : d ⨯₃ nvec f fv S = 0 := by
    rw [nvec, cross_cross_eq_smul_sub_smul', hd2, dotProduct_comm _ d, hd1]; simp
  have hdn' : d ⬝ᵥ nvec f fv S = 0 := by rw [dotProduct_comm]; exact hdn
  have hlag := cross_vnorm_sq d (nvec f fv S)
  rw [hx, show vnorm (0 : V3) = 0 by simp [DDG.vnorm], hdn'] at hlag
  have hnn : 0 < nvec f fv S ⬝ᵥ nvec f fv S := by
    rw [← vnorm_mul_self]; have := vnorm_pos_iff.mpr hS; positivity
  have hdd : d ⬝ᵥ d = 0 := by
    have : (d ⬝ᵥ d) * (nvec f fv S ⬝ᵥ nvec f fv S) = 0 := by nlinarith
    rcases mul_eq_zero.mp this with h | h
    · exact h
    · linarith
  have : d = 0 := dotProduct_self_eq_zero.mp hdd
  rw [hd] at this
  exact sub_eq_zero.mp this

/-- The paper's printed variant with **unit** edge vectors (Sec. 3.2.1, a misprint, D2). -/
def dfGradUnit (f : V → V3) (fv : F → Fin 3 → V) (u : V → ℝ) (S : F) : V3 :=
  (1 / (2 * area f fv S)) •
    (normal f fv S ⨯₃ ∑ a, u (fv S a) • ((vnorm (opp f fv S a))⁻¹ • opp f fv S a))

theorem vnorm_eq_of_dot {v : V3} {c : ℝ} (hc : 0 ≤ c) (h : v ⬝ᵥ v = c ^ 2) : vnorm v = c := by
  rw [DDG.vnorm, h, Real.sqrt_sq hc]

/-- **The unit-edge form is not the gradient** (SPEC D2): on the triangle `(0,0,0), (2,0,0),
(0,1,0)` with `u = (0, 0, 1)`, `D_f u · (f₂ − f₀) = 1` (the P1 gradient) while the unit-edge form
gives `½`.

[formula: rs-df-gradient] -/
theorem df_unit_edge_counterexample :
    let f : Fin 3 → V3 := ![![0, 0, 0], ![2, 0, 0], ![0, 1, 0]]
    let fv : Unit → Fin 3 → Fin 3 := fun _ => id
    let u : Fin 3 → ℝ := ![0, 0, 1]
    dfGrad f fv u () ⬝ᵥ (f 2 - f 0) = 1 ∧ dfGradUnit f fv u () ⬝ᵥ (f 2 - f 0) = 1 / 2 := by
  intro f fv u
  have hn : nvec f fv () = ![0, 0, 2] := by
    apply v3_ext <;> simp [nvec, corner, f, fv, cross_apply]
  have hnn : vnorm (nvec f fv ()) = 2 := by
    rw [hn]; apply vnorm_eq_of_dot (by norm_num); simp [dotProduct, Fin.sum_univ_three]; norm_num
  have hS : Nondeg f fv () := by
    rw [Nondeg, hn]; intro h; have := congrFun h 2; simp at this
  constructor
  · have := (dfGrad_props hS u).2 0 2
    simpa [corner, f, fv, u] using this
  · have ha : area f fv () = 1 := by rw [area, hnn]; norm_num
    have hN : normal f fv () = ![0, 0, 1] := by
      rw [normal, hnn, hn]; apply v3_ext <;> simp
    have ho2 : vnorm (opp f fv () 2) = 2 := by
      apply vnorm_eq_of_dot (by norm_num)
      simp [opp, corner, f, fv, dotProduct, Fin.sum_univ_three]; norm_num
    simp only [dfGradUnit, ha, hN, Fin.sum_univ_three, ho2]
    simp [u, fv, f, opp, corner, cross_apply, dotProduct, Fin.sum_univ_three]

end df

/-! ## Assembly by increments, the forms of Eq. 11–13 -/

section forms

variable {V F ι : Type*} [Fintype V] [Fintype F] [Fintype ι] [DecidableEq V]

/-- The pair form `Σ_{S,T} w_ST Σ_c (x_{S,c} − x_{T,c})(y_{S,c} − y_{T,c})` of Eq. 11–13 (the
`S = T` terms vanish, so `w` may include them or not). -/
def pairForm (w : F → F → ℝ) (x y : F → ι → ℝ) : ℝ :=
  ∑ S, ∑ T, w S T * ∑ c, (x S c - x T c) * (y S c - y T c)

/-- The matrix assembled by increments over ordered pairs (Alg. 1): for each `(S, T)` and each
feature component `c`, `w_ST (r_{S,c} − r_{T,c})(r_{S,c} − r_{T,c})ᵀ`, where `r_{S,c} = R (S, c)`
is a row of the feature matrix (`ū` or `D_f`; at most six vertices per pair).

[formula: rs-sobolev-forms] -/
def gramOf (R : Matrix (F × ι) V ℝ) (w : F → F → ℝ) : Matrix V V ℝ :=
  ∑ S, ∑ T, w S T • ∑ c, vecMulVec (R (S, c) - R (T, c)) (R (S, c) - R (T, c))

omit [DecidableEq V] in
theorem dot_vecMulVec (p q u v : V → ℝ) : u ⬝ᵥ vecMulVec p q *ᵥ v = (u ⬝ᵥ p) * (q ⬝ᵥ v) := by
  simp only [vecMulVec, mulVec, dotProduct, Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  simp only [of_apply]; ring

/-- The feature values `(R u)(S, c)` as a function of `S` and `c`. -/
def feat (R : Matrix (F × ι) V ℝ) (u : V → ℝ) : F → ι → ℝ := fun S c => (R *ᵥ u) (S, c)

omit [DecidableEq V] in
/-- **The increments assemble the form** (Eq. 11–13): `uᵀ M v = Σ_{S,T} w_ST Σ_c
((Ru)_{S,c} − (Ru)_{T,c})((Rv)_{S,c} − (Rv)_{T,c})`.

[formula: rs-sobolev-forms] -/
theorem gramOf_form (R : Matrix (F × ι) V ℝ) (w : F → F → ℝ) (u v : V → ℝ) :
    u ⬝ᵥ gramOf R w *ᵥ v = pairForm w (feat R u) (feat R v) := by
  simp only [gramOf, pairForm, feat, sum_mulVec, dotProduct_sum, smul_mulVec, dotProduct_smul,
    dot_vecMulVec, smul_eq_mul]
  refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
  congr 1
  refine Finset.sum_congr rfl fun c _ => ?_
  simp only [mulVec, sub_dotProduct, dotProduct_comm u]

omit [Fintype V] [DecidableEq V] in
theorem pairForm_symm (w : F → F → ℝ) (x y : F → ι → ℝ) : pairForm w x y = pairForm w y x := by
  simp only [pairForm, mul_comm (x _ _ - x _ _)]

omit [Fintype V] [DecidableEq V] in
theorem pairForm_self_nonneg {w : F → F → ℝ} (hw : ∀ S T, 0 ≤ w S T) (x : F → ι → ℝ) :
    0 ≤ pairForm w x x :=
  Finset.sum_nonneg fun S _ => Finset.sum_nonneg fun T _ =>
    mul_nonneg (hw S T) (Finset.sum_nonneg fun _ _ => mul_self_nonneg _)

/-- **The assembled operator is symmetric positive semidefinite** for nonnegative weights (whether
or not `w` is symmetric: `B⁰`'s `K₂` is not).

[formula: rs-sobolev-forms] -/
theorem gramOf_posSemidef (R : Matrix (F × ι) V ℝ) {w : F → F → ℝ} (hw : ∀ S T, 0 ≤ w S T) :
    (gramOf R w).PosSemidef :=
  posSemidef_of_form (fun u v => by rw [gramOf_form, gramOf_form, pairForm_symm])
    (fun u => by rw [gramOf_form]; exact pairForm_self_nonneg hw _)

/-- **Constants are in the kernel** when every feature takes the same value on the constant
vector on all faces (`ū(1) = 1`, `D_f 1 = 0`).

[formula: rs-sobolev-forms] -/
theorem gramOf_const (R : Matrix (F × ι) V ℝ) (w : F → F → ℝ) {k : ι → ℝ}
    (hR : ∀ S c, (R *ᵥ fun _ => (1 : ℝ)) (S, c) = k c) :
    gramOf R w *ᵥ (fun _ => 1) = 0 := by
  ext x
  have := gramOf_form R w (Pi.single x 1) (fun _ => 1)
  simp only [pairForm, feat, hR, sub_self, mul_zero, Finset.sum_const_zero] at this
  simpa [dotProduct, Pi.single_apply] using this

/-- The averaging rows: `(Avg u)(S) = ū_S = (u_i + u_j + u_k)/3`. -/
def avgMat (fv : F → Fin 3 → V) : Matrix (F × Unit) V ℝ :=
  fun Sc v => ∑ a : Fin 3, if fv Sc.1 a = v then (1 / 3 : ℝ) else 0

omit [Fintype F] in
theorem avgMat_mulVec (fv : F → Fin 3 → V) (u : V → ℝ) (S : F) (c : Unit) :
    (avgMat fv *ᵥ u) (S, c) = (1 / 3) * ∑ a, u (fv S a) := by
  simp only [avgMat, mulVec, dotProduct, Finset.sum_mul, ite_mul, zero_mul]
  rw [Finset.sum_comm]
  simp [Finset.mul_sum]

/-- The component rows of `D_f`: `(D u)(S, c) = (D_f u(S))_c`. -/
def dfMat (f : V → V3) (fv : F → Fin 3 → V) : Matrix (F × Fin 3) V ℝ :=
  fun Sc v => ∑ a : Fin 3, if fv Sc.1 a = v then dvec f fv Sc.1 a Sc.2 else 0

omit [Fintype F] in
theorem dfMat_mulVec (f : V → V3) (fv : F → Fin 3 → V) (u : V → ℝ) (S : F) (c : Fin 3) :
    (dfMat f fv *ᵥ u) (S, c) = dfGrad f fv u S c := by
  rw [dfGrad_eq_sum]
  simp only [dfMat, mulVec, dotProduct, Finset.sum_mul, ite_mul, zero_mul, Finset.sum_apply,
    Pi.smul_apply, smul_eq_mul]
  rw [Finset.sum_comm]
  simp [mul_comm]

/-- The fractional operator `L_τ` of Eq. 11 (weights `a_S a_T ρ^{−(2τ+2)}`, `S ≠ T`). -/
def opL [DecidableEq F] (f : V → V3) (fv : F → Fin 3 → V) (τ : ℝ) : Matrix V V ℝ :=
  gramOf (avgMat fv) fun S T =>
    if S ≠ T then area f fv S * area f fv T * vnorm (bary f fv S - bary f fv T) ^ (-(2 * τ + 2))
    else 0

/-- The high-order operator `B` of Eq. 12 (`D_f` rows, weights `a_S a_T ρ^{−2s}`). -/
def opB [DecidableEq F] (f : V → V3) (fv : F → Fin 3 → V) (s : ℝ) : Matrix V V ℝ :=
  gramOf (dfMat f fv) fun S T =>
    if S ≠ T then area f fv S * area f fv T * vnorm (bary f fv S - bary f fv T) ^ (-(2 * s))
    else 0

/-- The low-order operator `B⁰` of Eq. 13 (`K₂(S, T) = |N_S·(X_S − X_T)|² ρ^{−4}` times `ρ^{−2s}`,
fixed `α = 2`, `β = 4`). -/
def opB0 [DecidableEq F] (f : V → V3) (fv : F → Fin 3 → V) (s : ℝ) : Matrix V V ℝ :=
  gramOf (avgMat fv) fun S T =>
    if S ≠ T then kernel 2 4 (bary f fv S) (normal f fv S) (bary f fv T) * area f fv S *
      area f fv T * vnorm (bary f fv S - bary f fv T) ^ (-(2 * s))
    else 0

/-- **Eq. 11–13 as quadratic forms, PSD, with the constants in the kernel** (SPEC §4.2):
`uᵀ L_τ v = Σ_{S≠T} (ū_S − ū_T)(v̄_S − v̄_T) a_S a_T ρ^{−(2τ+2)}`,
`uᵀ B v = Σ_{S≠T} ⟨D_f u(S) − D_f u(T), D_f v(S) − D_f v(T)⟩ a_S a_T ρ^{−2s}`,
`uᵀ B⁰ v = Σ_{S≠T} (ū_S − ū_T)(v̄_S − v̄_T) K₂(S, T) a_S a_T ρ^{−2s}`; all three are PSD and
annihilate constants.

[formula: rs-sobolev-forms] -/
theorem sobolev_forms [DecidableEq F] (f : V → V3) (fv : F → Fin 3 → V) (τ s : ℝ)
    (u v : V → ℝ) :
    let ρ := fun S T => vnorm (bary f fv S - bary f fv T)
    let ub := fun S => (1 / 3 : ℝ) * ∑ a, u (fv S a)
    let vb := fun S => (1 / 3 : ℝ) * ∑ a, v (fv S a)
    (u ⬝ᵥ opL f fv τ *ᵥ v = ∑ S, ∑ T, if S ≠ T then
        (ub S - ub T) * (vb S - vb T) * (area f fv S * area f fv T * ρ S T ^ (-(2 * τ + 2)))
        else 0) ∧
      (u ⬝ᵥ opB f fv s *ᵥ v = ∑ S, ∑ T, if S ≠ T then
        ((dfGrad f fv u S - dfGrad f fv u T) ⬝ᵥ (dfGrad f fv v S - dfGrad f fv v T)) *
          (area f fv S * area f fv T * ρ S T ^ (-(2 * s))) else 0) ∧
      (u ⬝ᵥ opB0 f fv s *ᵥ v = ∑ S, ∑ T, if S ≠ T then
        (ub S - ub T) * (vb S - vb T) * (kernel 2 4 (bary f fv S) (normal f fv S) (bary f fv T) *
          area f fv S * area f fv T * ρ S T ^ (-(2 * s))) else 0) ∧
      (opL f fv τ).PosSemidef ∧ (opB f fv s).PosSemidef ∧ (opB0 f fv s).PosSemidef ∧
      opL f fv τ *ᵥ (fun _ => 1) = 0 ∧ opB f fv s *ᵥ (fun _ => 1) = 0 ∧
      opB0 f fv s *ᵥ (fun _ => 1) = 0 := by
  intro ρ ub vb
  have hρ : ∀ S T (e : ℝ), 0 ≤ vnorm (bary f fv S - bary f fv T) ^ e :=
    fun S T e => Real.rpow_nonneg (vnorm_nonneg _) e
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [opL, gramOf_form, pairForm]
    refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
    by_cases h : S ≠ T
    · simp only [feat, avgMat_mulVec, Finset.univ_unique, Finset.sum_singleton, ub, vb, ρ, h,
        ite_true, ne_eq, not_false_eq_true]; ring
    · simp [h]
  · rw [opB, gramOf_form, pairForm]
    refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
    by_cases h : S ≠ T
    · simp only [feat, dfMat_mulVec, dotProduct, Pi.sub_apply, ρ, h, ite_true, ne_eq,
        not_false_eq_true]; ring
    · simp [h]
  · rw [opB0, gramOf_form, pairForm]
    refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
    by_cases h : S ≠ T
    · simp only [feat, avgMat_mulVec, Finset.univ_unique, Finset.sum_singleton, ub, vb, ρ, h,
        ite_true, ne_eq, not_false_eq_true]; ring
    · simp [h]
  · exact gramOf_posSemidef _ fun S T => by
      split_ifs
      · exact mul_nonneg (mul_nonneg (area_nonneg _ _ _) (area_nonneg _ _ _)) (hρ S T _)
      · exact le_rfl
  · exact gramOf_posSemidef _ fun S T => by
      split_ifs
      · exact mul_nonneg (mul_nonneg (area_nonneg _ _ _) (area_nonneg _ _ _)) (hρ S T _)
      · exact le_rfl
  · exact gramOf_posSemidef _ fun S T => by
      split_ifs
      · exact mul_nonneg (mul_nonneg (mul_nonneg (kernel_nonneg _ _ _ _ _) (area_nonneg _ _ _))
          (area_nonneg _ _ _)) (hρ S T _)
      · exact le_rfl
  · exact gramOf_const _ _ (k := fun _ => 1) fun S c => by
      rw [avgMat_mulVec]; simp
  · exact gramOf_const _ _ (k := fun _ => 0) fun S c => by
      rw [dfMat_mulVec, dfGrad_const]; rfl
  · exact gramOf_const _ _ (k := fun _ => 1) fun S c => by
      rw [avgMat_mulVec]; simp

end forms

/-! ## Symmetrisation and the kernel-matrix form (App. A) -/

section appA

variable {F : Type*} [Fintype F] [DecidableEq F]

omit [DecidableEq F] in
/-- **The low-order symmetrisation** (SPEC §4.2): for a factor `g` symmetric in `(S, T)`, the
ordered sum with the kernel `K` equals the sum with `½(K_ST + K_TS)`.

[formula: rs-low-order-symmetrisation] -/
theorem low_order_symmetrisation (g K : F → F → ℝ) (hg : ∀ S T, g S T = g T S) :
    ∑ S, ∑ T, g S T * K S T = ∑ S, ∑ T, g S T * ((K S T + K T S) / 2) := by
  have h : ∑ S, ∑ T, g S T * K T S = ∑ S, ∑ T, g S T * K S T := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => by rw [hg]
  have e : ∑ S, ∑ T, g S T * ((K S T + K T S) / 2) =
      (∑ S, ∑ T, g S T * K S T + ∑ S, ∑ T, g S T * K T S) / 2 := by
    rw [← Finset.sum_add_distrib, Finset.sum_div]
    refine Finset.sum_congr rfl fun S _ => ?_
    rw [← Finset.sum_add_distrib, Finset.sum_div]
    refine Finset.sum_congr rfl fun T _ => ?_
    ring
  rw [e, h]; ring

/-- The App. A operator `diag(a)⁻¹ diag(H a) − H`. -/
def appA (H : Matrix F F ℝ) (a : F → ℝ) : Matrix F F ℝ :=
  diagonal (fun S => (a S)⁻¹ * (H *ᵥ a) S) - H

/-- **The App. A identity, scalar version** (SPEC §4.3): for symmetric `H` and `a_S ≠ 0`,
`Σ_{S,T} (x_S − x_T)(y_S − y_T) a_S H_ST a_T = 2 (a⊙x)ᵀ (diag(a)⁻¹ diag(H a) − H) (a⊙y)`.

[formula: rs-kernel-matrix-form] -/
theorem kernel_matrix_form {H : Matrix F F ℝ} (hH : Hᵀ = H) {a : F → ℝ} (ha : ∀ S, a S ≠ 0)
    (x y : F → ℝ) :
    ∑ S, ∑ T, (x S - x T) * (y S - y T) * (a S * H S T * a T) =
      2 * ((a * x) ⬝ᵥ appA H a *ᵥ (a * y)) := by
  have hs : ∀ S T, H T S = H S T := fun S T => by rw [← transpose_apply H S T, hH]
  -- the right side
  have r1 : (a * x) ⬝ᵥ appA H a *ᵥ (a * y) =
      ∑ S, ∑ T, a S * H S T * a T * (x S * y S) - ∑ S, ∑ T, a S * H S T * a T * (x S * y T) := by
    rw [appA, sub_mulVec, dotProduct_sub, show (diagonal fun S => (a S)⁻¹ * (H *ᵥ a) S) *ᵥ (a * y)
      = fun S => (a S)⁻¹ * (H *ᵥ a) S * (a * y) S from funext fun S => mulVec_diagonal _ _ S]
    congr 1
    · simp only [dotProduct, Pi.mul_apply, mulVec, Finset.mul_sum, Finset.sum_mul]
      refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
      field_simp [ha S]
    · simp only [dotProduct, Pi.mul_apply, mulVec, Finset.mul_sum]
      refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
      ring
  have e2 : ∑ S, ∑ T, a S * H S T * a T * (x T * y T) =
      ∑ S, ∑ T, a S * H S T * a T * (x S * y S) := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => by rw [hs]; ring
  have e3 : ∑ S, ∑ T, a S * H S T * a T * (x T * y S) =
      ∑ S, ∑ T, a S * H S T * a T * (x S * y T) := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => by rw [hs]; ring
  have l : ∑ S, ∑ T, (x S - x T) * (y S - y T) * (a S * H S T * a T) =
      ∑ S, ∑ T, a S * H S T * a T * (x S * y S) + ∑ S, ∑ T, a S * H S T * a T * (x T * y T)
        - ∑ S, ∑ T, a S * H S T * a T * (x S * y T)
        - ∑ S, ∑ T, a S * H S T * a T * (x T * y S) := by
    simp only [← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun S _ => Finset.sum_congr rfl fun T _ => ?_
    ring
  rw [l, e2, e3, r1]; ring

/-- **The App. A identity, vector version** (for `V u = a ⊙ D_f u`): summed over the components.

[formula: rs-kernel-matrix-form] -/
theorem kernel_matrix_form_vec {ι : Type*} [Fintype ι] {H : Matrix F F ℝ} (hH : Hᵀ = H)
    {a : F → ℝ} (ha : ∀ S, a S ≠ 0) (x y : F → ι → ℝ) :
    ∑ S, ∑ T, (∑ c, (x S c - x T c) * (y S c - y T c)) * (a S * H S T * a T) =
      2 * ∑ c, ((a * fun S => x S c) ⬝ᵥ appA H a *ᵥ (a * fun S => y S c)) := by
  have hc : ∀ c, 2 * ((a * fun S => x S c) ⬝ᵥ appA H a *ᵥ (a * fun S => y S c)) =
      ∑ S, ∑ T, (x S c - x T c) * (y S c - y T c) * (a S * H S T * a T) :=
    fun c => (kernel_matrix_form hH ha _ _).symm
  rw [Finset.mul_sum]
  simp only [hc, Finset.sum_mul]
  symm
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun S _ => ?_
  rw [Finset.sum_comm]

/-- **The hierarchical null space** (SPEC §4.3): `(diag(a)⁻¹ diag(H a) − H) a = 0` for **any**
`H` (in particular any compressed, possibly non-symmetric, approximation), provided the diagonal
`H a` is computed with the same `H`.

[formula: rs-hierarchical-null-space] -/
theorem appA_mulVec_self (H : Matrix F F ℝ) {a : F → ℝ} (ha : ∀ S, a S ≠ 0) :
    appA H a *ᵥ a = 0 := by
  ext S
  simp only [appA, sub_mulVec, mulVec_diagonal, Pi.sub_apply, Pi.zero_apply]
  field_simp [ha S]
  ring

end appA

section appAmesh

variable {V F : Type*} [Fintype V] [Fintype F] [DecidableEq V] [DecidableEq F]

/-- The paper's `U = diag(a) Avg` (`U u = a ⊙ ū`). -/
def Umat (a : F → ℝ) (fv : F → Fin 3 → V) : Matrix F V ℝ :=
  fun S v => a S * ∑ i : Fin 3, if fv S i = v then (1 / 3 : ℝ) else 0

omit [Fintype F] [DecidableEq F] in
theorem Umat_mulVec (a : F → ℝ) (fv : F → Fin 3 → V) (u : V → ℝ) (S : F) :
    (Umat a fv *ᵥ u) S = a S * ((1 / 3) * ∑ i, u (fv S i)) := by
  have := avgMat_mulVec fv u S ()
  simp only [avgMat, mulVec, dotProduct] at this
  simp only [Umat, mulVec, dotProduct, mul_assoc, ← Finset.mul_sum, this]

/-- **The App. A operator form of `L`** (SPEC §4.3, §8.3): for symmetric `H` and nonzero `a`,
`Σ_{S,T} (ū_S − ū_T)(v̄_S − v̄_T) a_S H_ST a_T = uᵀ (2 Uᵀ (diag(a)⁻¹ diag(H a) − H) U) v`.

[formula: rs-kernel-matrix-form] -/
theorem kernel_matrix_form_U {H : Matrix F F ℝ} (hH : Hᵀ = H) {a : F → ℝ} (ha : ∀ S, a S ≠ 0)
    (fv : F → Fin 3 → V) (u v : V → ℝ) :
    ∑ S, ∑ T, ((1 / 3) * ∑ i, u (fv S i) - (1 / 3) * ∑ i, u (fv T i)) *
        ((1 / 3) * ∑ i, v (fv S i) - (1 / 3) * ∑ i, v (fv T i)) * (a S * H S T * a T) =
      u ⬝ᵥ ((2 : ℝ) • ((Umat a fv)ᵀ * appA H a * Umat a fv)) *ᵥ v := by
  rw [kernel_matrix_form hH ha (fun S => (1 / 3) * ∑ i, u (fv S i))
    (fun S => (1 / 3) * ∑ i, v (fv S i))]
  have hU : ∀ w : V → ℝ, (a * fun S => (1 / 3) * ∑ i, w (fv S i)) = Umat a fv *ᵥ w := by
    intro w; funext S; rw [Umat_mulVec]; rfl
  have hT : ∀ (M : Matrix F V ℝ) (u : V → ℝ) (w : F → ℝ), u ⬝ᵥ Mᵀ *ᵥ w = (M *ᵥ u) ⬝ᵥ w := by
    intro M u w; rw [dotProduct_mulVec, vecMul_transpose]
  rw [hU, hU, smul_mulVec, dotProduct_smul, smul_eq_mul, ← mulVec_mulVec, ← mulVec_mulVec, hT]

/-- **`L 1 = 0` for every compressed `H`** (SPEC §4.3, `rs-hierarchical-null-space`): with
`U 1 = a`, `Uᵀ (diag(a)⁻¹ diag(H a) − H) U 1 = 0`.

[formula: rs-hierarchical-null-space] -/
theorem hierarchical_null_space (H : Matrix F F ℝ) {a : F → ℝ} (ha : ∀ S, a S ≠ 0)
    (fv : F → Fin 3 → V) :
    ((Umat a fv)ᵀ * appA H a * Umat a fv) *ᵥ (fun _ => 1) = 0 := by
  have hU : Umat a fv *ᵥ (fun _ => (1 : ℝ)) = a := by
    funext S; rw [Umat_mulVec]; simp
  rw [← mulVec_mulVec, hU, ← mulVec_mulVec, appA_mulVec_self H ha, mulVec_zero]

omit [Fintype V] [Fintype F] [DecidableEq V] [DecidableEq F] in
/-- **The `V` operator annihilates constants**: `V 1 = a ⊙ D_f 1 = 0`, so `B 1 = 0` for any
compressed kernel.

[formula: rs-hierarchical-null-space] -/
theorem V_const (f : V → V3) (fv : F → Fin 3 → V) (a : F → ℝ) (S : F) (c : Fin 3) :
    a S * dfGrad f fv (fun _ => 1) S c = 0 := by
  rw [dfGrad_const]; simp

end appAmesh

/-! ## The projector average (SPEC §8.3) -/

section projector

variable {F : Type*} [Fintype F]

/-- **The projector average** (SPEC §8.3): `P̄ = Σ_S w_S N_S N_Sᵀ` with unit normals has
`vᵀ P̄ v = Σ_S w_S |P_S v|²` (`= Σ_S w_S (N_S·v)²`); for a convex combination, `0 ≤ vᵀP̄v ≤ |v|²`.

[formula: rs-projector-average] -/
theorem projector_average (w : F → ℝ) (N : F → V3) (hN : ∀ S, N S ⬝ᵥ N S = 1) (v : V3) :
    v ⬝ᵥ (∑ S, w S • vecMulVec (N S) (N S)) *ᵥ v = ∑ S, w S * vnorm (vecMulVec (N S) (N S) *ᵥ v) ^ 2 ∧
      ((∀ S, 0 ≤ w S) → ∑ S, w S = 1 →
        0 ≤ v ⬝ᵥ (∑ S, w S • vecMulVec (N S) (N S)) *ᵥ v ∧
          v ⬝ᵥ (∑ S, w S • vecMulVec (N S) (N S)) *ᵥ v ≤ v ⬝ᵥ v) := by
  have h1 : v ⬝ᵥ (∑ S, w S • vecMulVec (N S) (N S)) *ᵥ v = ∑ S, w S * (N S ⬝ᵥ v) ^ 2 := by
    simp only [sum_mulVec, dotProduct_sum, smul_mulVec, dotProduct_smul, smul_eq_mul,
      (projector_form (hN _) v 0).2.2.1]
  have hle : ∀ S, (N S ⬝ᵥ v) ^ 2 ≤ v ⬝ᵥ v := by
    intro S
    have := cross_vnorm_sq (N S) v
    rw [hN S, one_mul] at this
    nlinarith [sq_nonneg (vnorm (N S ⨯₃ v))]
  refine ⟨?_, fun hw h1w => ⟨?_, ?_⟩⟩
  · rw [h1]
    exact Finset.sum_congr rfl fun S _ => by rw [(projector_form (hN S) v 0).2.1]
  · rw [h1]; exact Finset.sum_nonneg fun S _ => mul_nonneg (hw S) (sq_nonneg _)
  · rw [h1]
    calc ∑ S, w S * (N S ⬝ᵥ v) ^ 2 ≤ ∑ S, w S * (v ⬝ᵥ v) :=
          Finset.sum_le_sum fun S _ => mul_le_mul_of_nonneg_left (hle S) (hw S)
      _ = v ⬝ᵥ v := by rw [← Finset.sum_mul, h1w, one_mul]

end projector

end RS

end DDG
