import DdgProofs.RepulsiveCurves.Energy

/-!
# Repulsive Curves: the inner products (paper Sec. 5.2, Eq. 19, 21, 31, 32; SPEC §4)

**Assembly.** A matrix assembled "by increments" (`B_{xy} += c`) is a finite sum of single-entry
matrices `incr x y c`; the theorems below are about exactly these sums, so they are statements
about the assembly loops of paper Sec. 5.2.2 and 5.2.3. Pairs are any finite set `S` of ordered
pairs of edges (in the application `𝒫`), weights `w : E → E → ℝ` are abstract (both quadratures of
SPEC D2 are instances, see below), edge data `ℓ : E → ℝ`, `T : E → ℝ³`.

* `rc-derivative-operator`: `(D u)_I = (u_{i₂} − u_{i₁}) T_I / ℓ_I`, `D_I = (1/ℓ_I)[−T_I T_I]`
  (`Dmat`), `(E u)_I = (u_{i₁} + u_{i₂})/2`; `D 𝟙 = 0`, `E 𝟙 = 𝟙`, `ℓ_I T_I·(D u)_I = u_{i₂} − u_{i₁}`.
* `rc-high-order-gram`: the increments of Sec. 5.2.2 give `uᵀ B v = Σ_{(I,J)∈S} w_IJ ⟨D_I u − D_J u,
  D_I v − D_J v⟩` (Eq. 21) for unit tangents (`gramB_form`).
* `rc-low-order-gram`: the increments of Sec. 5.2.3 give `uᵀ B⁰ v = Σ w⁰_IJ (u_I − u_J)(v_I − v_J)`.
* `rc-kernel-matrix-form`: Eq. 31: `Σ_IJ (diag(K𝟙) − K)_IJ ⟨X_I, Y_J⟩ = ½ Σ_IJ K_IJ ⟨X_I − X_J,
  Y_I − Y_J⟩` for symmetric `K`; with `K = W + Wᵀ` (`W` = the weights extended by `0` off `S`)
  `B = Dᵀ((diag(K𝟙) − K) ⊗ I₃) D` and `B⁰ = Eᵀ(diag(K⁰𝟙) − K⁰) E` as matrices; the midpoint weights
  give the kernels of Eq. 32; the D2 identity is `kernel_sq_quadrature` (`Kernel.lean`).
* `rc-gram-psd-kernel`: `B`, `B⁰` symmetric positive semidefinite for `w ≥ 0`; `B u = 0` iff
  `D_I u = D_J u` on the pairs with `w_IJ > 0`; `B⁰ u = 0` iff `u_I = u_J` on the pairs with
  `w⁰_IJ > 0`; constants along every edge are in `ker B`; `ker(B + B⁰) = ker B ∩ ker B⁰` for any two
  PSD matrices; if every edge is coupled with positive weight to an edge with an independent
  tangent and the edges are connected through `w⁰ > 0`, `ker A` is the constants and `ker Ā` the
  translations. **Not** "`ker B⁰` = constants": `ker E ⊆ ker B⁰` (`u_{i₁} = −u_{i₂}` on every
  edge: alternating values on a bipartite component, i.e. every open curve, tree or even loop), so
  `ker B⁰` is larger than the constants on such networks (`gramB0_mulVec_eq_zero_of_alternating`).
* `rc-gram-scaling`: `w(cγ) = c^{1−2σ} w(γ)` (both quadratures), `w⁰(cγ) = c^{−(2σ+1)} w⁰(γ)`
  (both), hence `B(cγ) = c^{−(2σ+1)} B(γ)`, `B⁰(cγ) = c^{−(2σ+1)} B⁰(γ)`.
* `rc-classic-metrics` (SPEC §4.6): the H¹ matrix `L` (increments `±1/ℓ`), `uᵀ L v = Σ_I (u_{i₂} −
  u_{i₁})(v_{i₂} − v_{i₁})/ℓ_I`, PSD, `ker L` = functions constant along every edge (the
  per-component constants); `H² = L M⁻¹ L` has the same kernel for positive masses.
-/

noncomputable section

namespace DdgProofs

namespace RC

open Matrix Real Finset
open scoped Kronecker

variable {V E : Type*} [Fintype V] [Fintype E] [DecidableEq V] [DecidableEq E]

/-! ## Assembly by increments -/

section incr

/-- The matrix with the single entry `c` at `(i, j)`: one increment `M_{ij} += c`. -/
def incr (i j : V) (c : ℝ) : Matrix V V ℝ := fun x y => if x = i ∧ y = j then c else 0

omit [Fintype E] [DecidableEq E] in
theorem dot_incr (u v : V → ℝ) (i j : V) (c : ℝ) : u ⬝ᵥ incr i j c *ᵥ v = u i * c * v j := by
  have : incr i j c *ᵥ v = fun x => if x = i then c * v j else 0 := by
    funext x
    simp only [incr, mulVec, dotProduct]
    split_ifs with h
    · simp [h]
    · simp [h]
  rw [this]
  simp [dotProduct, mul_assoc]

omit [Fintype E] [DecidableEq E] [DecidableEq V] in
theorem dot_sum_mulVec {ι : Type*} (s : Finset ι) (M : ι → Matrix V V ℝ) (u v : V → ℝ) :
    u ⬝ᵥ (∑ k ∈ s, M k) *ᵥ v = ∑ k ∈ s, u ⬝ᵥ M k *ᵥ v := by
  rw [sum_mulVec, dotProduct_sum]

omit [Fintype E] [DecidableEq E] [DecidableEq V] in
theorem dot_add_mulVec (M N : Matrix V V ℝ) (u v : V → ℝ) :
    u ⬝ᵥ (M + N) *ᵥ v = u ⬝ᵥ M *ᵥ v + u ⬝ᵥ N *ᵥ v := by
  rw [add_mulVec, dotProduct_add]

omit [Fintype E] [DecidableEq E] in
/-- An entry is the bilinear form on two unit vectors. -/
theorem entry_eq_dot (M : Matrix V V ℝ) (x y : V) :
    M x y = Pi.single x 1 ⬝ᵥ M *ᵥ Pi.single y 1 := by
  simp [mulVec, dotProduct, Pi.single_apply]

omit [Fintype E] [DecidableEq E] in
/-- Two matrices with the same bilinear form are equal. -/
theorem matrix_eq_of_form {M N : Matrix V V ℝ} (h : ∀ u v, u ⬝ᵥ M *ᵥ v = u ⬝ᵥ N *ᵥ v) : M = N := by
  ext x y
  rw [entry_eq_dot M, entry_eq_dot N, h]

omit [Fintype E] [DecidableEq E] in
/-- A matrix with a symmetric bilinear form is symmetric. -/
theorem transpose_eq_of_form {M : Matrix V V ℝ} (h : ∀ u v, u ⬝ᵥ M *ᵥ v = v ⬝ᵥ M *ᵥ u) :
    Mᵀ = M := by
  ext x y
  rw [transpose_apply, entry_eq_dot M y x, entry_eq_dot M x y, h]

omit [Fintype E] [DecidableEq E] in
/-- A real matrix with a symmetric, nonnegative quadratic form is positive semidefinite. -/
theorem posSemidef_of_form {M : Matrix V V ℝ} (hs : ∀ u v, u ⬝ᵥ M *ᵥ v = v ⬝ᵥ M *ᵥ u)
    (hn : ∀ u, 0 ≤ u ⬝ᵥ M *ᵥ u) : M.PosSemidef := by
  refine posSemidef_iff_dotProduct_mulVec.mpr ⟨?_, fun x => by simpa using hn x⟩
  rw [IsHermitian, conjTranspose_eq_transpose_of_trivial]
  exact transpose_eq_of_form hs

omit [Fintype E] [DecidableEq E] [DecidableEq V] in
/-- For a PSD real matrix, `M u = 0 ⟺ uᵀ M u = 0`. -/
theorem mulVec_eq_zero_iff_form {M : Matrix V V ℝ} (hM : M.PosSemidef) (u : V → ℝ) :
    M *ᵥ u = 0 ↔ u ⬝ᵥ M *ᵥ u = 0 := by
  rw [← hM.dotProduct_mulVec_zero_iff]
  simp

omit [Fintype E] [DecidableEq E] [DecidableEq V] in
theorem psd_form_nonneg {M : Matrix V V ℝ} (hM : M.PosSemidef) (u : V → ℝ) :
    0 ≤ u ⬝ᵥ M *ᵥ u := by
  simpa using (posSemidef_iff_dotProduct_mulVec.mp hM).2 u

omit [Fintype E] [DecidableEq E] [DecidableEq V] in
/-- **`ker(M + N) = ker M ∩ ker N`** for positive semidefinite `M`, `N`.

[formula: rc-gram-psd-kernel] -/
theorem ker_add_posSemidef {M N : Matrix V V ℝ} (hM : M.PosSemidef) (hN : N.PosSemidef)
    (u : V → ℝ) : (M + N) *ᵥ u = 0 ↔ M *ᵥ u = 0 ∧ N *ᵥ u = 0 := by
  rw [mulVec_eq_zero_iff_form (hM.add hN), mulVec_eq_zero_iff_form hM,
    mulVec_eq_zero_iff_form hN, dot_add_mulVec]
  have h1 := psd_form_nonneg hM u
  have h2 := psd_form_nonneg hN u
  constructor
  · intro h; constructor <;> linarith
  · rintro ⟨h1, h2⟩; rw [h1, h2, add_zero]

end incr

/-! ## The operators `D` and `E` -/

section ops

variable (ℓ : E → ℝ) (T : E → V3) (ep : E → Fin 2 → V)

/-- The discrete derivative `(D u)_I = (u_{i₂} − u_{i₁}) T_I / ℓ_I` (Sec. 5.2.1).

[formula: rc-derivative-operator] -/
def Dop (u : V → ℝ) (I : E) : V3 := ((u (ep I 1) - u (ep I 0)) / ℓ I) • T I

/-- The averaging operator `(E u)_I = (u_{i₁} + u_{i₂}) / 2` (B.3.4).

[formula: rc-derivative-operator] -/
def Eop (u : V → ℝ) (I : E) : ℝ := (u (ep I 0) + u (ep I 1)) / 2

/-- `D ∈ R^{3|E|×|V|}` as a matrix on the interleaved index `(I, c)`; its block on the two
endpoints of `I` is `D_I = (1/ℓ_I)[−T_I T_I]`.

[formula: rc-derivative-operator] -/
def Dmat : Matrix (E × Fin 3) V ℝ := fun Ic v =>
  (∑ a : Fin 2, if ep Ic.1 a = v then sgn a else 0) * T Ic.1 Ic.2 / ℓ Ic.1

/-- `E ∈ R^{|E|×|V|}` as a matrix. -/
def Emat : Matrix E V ℝ := fun I v => ∑ a : Fin 2, if ep I a = v then (1 / 2 : ℝ) else 0

omit [Fintype E] [DecidableEq E] in
/-- **The derivative operator** (Sec. 5.2.1): `D u` is the product with `Dmat`, `D_I u =
Σ_a s_a u_{i_a} T_I/ℓ_I`, constants have derivative `0`, and for a unit `T_I`,
`ℓ_I T_I · (D u)_I = u_{i₂} − u_{i₁}`; `E u` is the product with `Emat` and `E 𝟙 = 𝟙`.

[formula: rc-derivative-operator] -/
theorem derivative_operator (u : V → ℝ) (I : E) :
    (∀ c, (Dmat ℓ T ep *ᵥ u) (I, c) = Dop ℓ T ep u I c) ∧
      Dop ℓ T ep u I = ∑ a, (sgn a * u (ep I a) / ℓ I) • T I ∧
      Dop ℓ T ep (fun _ => (1 : ℝ)) I = 0 ∧
      (T I ⬝ᵥ T I = 1 → ℓ I ≠ 0 → ℓ I * (T I ⬝ᵥ Dop ℓ T ep u I) = u (ep I 1) - u (ep I 0)) ∧
      (Emat ep *ᵥ u) I = Eop ep u I ∧ Eop ep (fun _ => (1 : ℝ)) I = 1 := by
  refine ⟨fun c => ?_, ?_, by simp [Dop], fun hT hℓ => ?_, ?_, by simp [Eop]⟩
  · simp only [Dmat, mulVec, dotProduct, Dop, Pi.smul_apply, smul_eq_mul]
    simp only [Finset.sum_mul, Finset.sum_div]
    rw [Finset.sum_comm]
    simp only [ite_mul, zero_mul, ite_div, zero_div, Finset.sum_ite_eq,
      Finset.mem_univ, ite_true, Fin.sum_univ_two, sgn_zero, sgn_one]
    ring
  · simp only [Dop, Fin.sum_univ_two, sgn_zero, sgn_one, ← add_smul]
    congr 1
    ring
  · rw [Dop, dotProduct_smul, hT, smul_eq_mul]
    field_simp
  · simp only [Emat, mulVec, dotProduct, Eop]
    simp only [Finset.sum_mul]
    rw [Finset.sum_comm]
    simp only [ite_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ,
      ite_true, Fin.sum_univ_two]
    ring

end ops

/-! ## The high-order term `B` (Eq. 21) -/

section highOrder

variable (S : Finset (E × E)) (w : E → E → ℝ) (ℓ : E → ℝ) (T : E → V3) (ep : E → Fin 2 → V)

/-- The quadratic form of Eq. 21, `Σ_{(I,J) ∈ S} w_IJ ⟨D_I u − D_J u, D_I v − D_J v⟩`. -/
def formB (u v : V → ℝ) : ℝ :=
  ∑ p ∈ S, w p.1 p.2 *
    ((Dop ℓ T ep u p.1 - Dop ℓ T ep u p.2) ⬝ᵥ (Dop ℓ T ep v p.1 - Dop ℓ T ep v p.2))

/-- The Gram matrix `B` assembled by the increments of Sec. 5.2.2: for every `(I, J) ∈ S` and
`a, b ∈ {1, 2}`, `B_{i_a i_b} += (−1)^{a+b} w_IJ/ℓ_I²`, `B_{i_a j_b} −= (−1)^{a+b} w_IJ ⟨T_I,
T_J⟩/(ℓ_I ℓ_J)`, `B_{j_a j_b} += (−1)^{a+b} w_IJ/ℓ_J²`, `B_{j_a i_b} −= (−1)^{a+b} w_IJ ⟨T_J,
T_I⟩/(ℓ_J ℓ_I)`.

[formula: rc-high-order-gram] -/
def gramB : Matrix V V ℝ :=
  ∑ p ∈ S, ∑ a : Fin 2, ∑ b : Fin 2,
    (incr (ep p.1 a) (ep p.1 b) ((-1) ^ ((a : ℕ) + (b : ℕ)) * w p.1 p.2 / ℓ p.1 ^ 2)
      + incr (ep p.1 a) (ep p.2 b)
          (-((-1) ^ ((a : ℕ) + (b : ℕ)) * w p.1 p.2 * (T p.1 ⬝ᵥ T p.2) / (ℓ p.1 * ℓ p.2)))
      + incr (ep p.2 a) (ep p.2 b) ((-1) ^ ((a : ℕ) + (b : ℕ)) * w p.1 p.2 / ℓ p.2 ^ 2)
      + incr (ep p.2 a) (ep p.1 b)
          (-((-1) ^ ((a : ℕ) + (b : ℕ)) * w p.1 p.2 * (T p.2 ⬝ᵥ T p.1) / (ℓ p.2 * ℓ p.1))))

omit [DecidableEq E] [Fintype E] in
/-- **The increments of Sec. 5.2.2 assemble the quadratic form of Eq. 21**: for unit tangents,
`uᵀ B v = Σ_{(I,J) ∈ S} w_IJ ⟨D_I u − D_J u, D_I v − D_J v⟩` for all `u, v`.

[formula: rc-high-order-gram] -/
theorem gramB_form (hT : ∀ I, T I ⬝ᵥ T I = 1) (u v : V → ℝ) :
    u ⬝ᵥ gramB S w ℓ T ep *ᵥ v = formB S w ℓ T ep u v := by
  simp only [gramB, dot_sum_mulVec, dot_add_mulVec, dot_incr, formB]
  refine Finset.sum_congr rfl fun p _ => ?_
  simp only [Dop, sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, smul_eq_mul,
    hT, Fin.sum_univ_two, Fin.val_zero, Fin.val_one, dotProduct_comm (T p.2) (T p.1)]
  ring

omit [DecidableEq E] [Fintype V] [Fintype E] [DecidableEq V] in
theorem formB_symm (u v : V → ℝ) : formB S w ℓ T ep u v = formB S w ℓ T ep v u := by
  unfold formB
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [dotProduct_comm]

omit [DecidableEq E] [Fintype V] [Fintype E] [DecidableEq V] in
theorem formB_self_nonneg (hw : ∀ p ∈ S, 0 ≤ w p.1 p.2) (u : V → ℝ) :
    0 ≤ formB S w ℓ T ep u u :=
  Finset.sum_nonneg fun p hp => mul_nonneg (hw p hp) (dot_self_nonneg _)

omit [DecidableEq E] [Fintype E] in
/-- **`B` is symmetric positive semidefinite** for nonnegative weights and unit tangents.

[formula: rc-gram-psd-kernel] -/
theorem gramB_posSemidef (hT : ∀ I, T I ⬝ᵥ T I = 1) (hw : ∀ p ∈ S, 0 ≤ w p.1 p.2) :
    (gramB S w ℓ T ep).PosSemidef :=
  posSemidef_of_form (fun u v => by rw [gramB_form S w ℓ T ep hT, gramB_form S w ℓ T ep hT,
    formB_symm]) (fun u => by rw [gramB_form S w ℓ T ep hT]; exact formB_self_nonneg S w ℓ T ep hw u)

omit [DecidableEq E] [Fintype E] in
/-- **The kernel of `B`**: `B u = 0` iff `D_I u = D_J u` for every `(I, J) ∈ S` with `w_IJ > 0`.

[formula: rc-gram-psd-kernel] -/
theorem gramB_ker (hT : ∀ I, T I ⬝ᵥ T I = 1) (hw : ∀ p ∈ S, 0 ≤ w p.1 p.2) (u : V → ℝ) :
    gramB S w ℓ T ep *ᵥ u = 0 ↔
      ∀ p ∈ S, 0 < w p.1 p.2 → Dop ℓ T ep u p.1 = Dop ℓ T ep u p.2 := by
  rw [mulVec_eq_zero_iff_form (gramB_posSemidef S w ℓ T ep hT hw), gramB_form S w ℓ T ep hT,
    formB, Finset.sum_eq_zero_iff_of_nonneg fun p hp => mul_nonneg (hw p hp) (dot_self_nonneg _)]
  refine forall₂_congr fun p hp => ?_
  rw [mul_eq_zero, dotProduct_self_eq_zero, sub_eq_zero]
  constructor
  · rintro (h | h) hpos
    · exact absurd h hpos.ne'
    · exact h
  · intro h
    rcases (hw p hp).eq_or_lt with h0 | h0
    · exact Or.inl h0.symm
    · exact Or.inr (h h0)

omit [DecidableEq E] [Fintype E] in
/-- **Functions constant along every edge (constant on each component) are in `ker B`.**

[formula: rc-gram-psd-kernel] -/
theorem gramB_mulVec_eq_zero_of_edge_const (hT : ∀ I, T I ⬝ᵥ T I = 1)
    (hw : ∀ p ∈ S, 0 ≤ w p.1 p.2) {u : V → ℝ} (hu : ∀ I, u (ep I 0) = u (ep I 1)) :
    gramB S w ℓ T ep *ᵥ u = 0 := by
  rw [gramB_ker S w ℓ T ep hT hw]
  intro p _ _
  simp [Dop, hu]

end highOrder

/-! ## The low-order term `B⁰` (Sec. 5.2.3) -/

section lowOrder

variable (S : Finset (E × E)) (w0 : E → E → ℝ) (ep : E → Fin 2 → V)

/-- The quadratic form `Σ_{(I,J) ∈ S} w⁰_IJ (u_I − u_J)(v_I − v_J)`, `u_I = (E u)_I`. -/
def formB0 (u v : V → ℝ) : ℝ :=
  ∑ p ∈ S, w0 p.1 p.2 * ((Eop ep u p.1 - Eop ep u p.2) * (Eop ep v p.1 - Eop ep v p.2))

/-- The Gram matrix `B⁰` assembled by the increments of Sec. 5.2.3: `B⁰_{i_a i_b} += ¼ w⁰_IJ`,
`B⁰_{i_a j_b} −= ¼ w⁰_IJ`, `B⁰_{j_a i_b} −= ¼ w⁰_IJ`, `B⁰_{j_a j_b} += ¼ w⁰_IJ`.

[formula: rc-low-order-gram] -/
def gramB0 : Matrix V V ℝ :=
  ∑ p ∈ S, ∑ a : Fin 2, ∑ b : Fin 2,
    (incr (ep p.1 a) (ep p.1 b) (w0 p.1 p.2 / 4) + incr (ep p.1 a) (ep p.2 b) (-(w0 p.1 p.2 / 4))
      + incr (ep p.2 a) (ep p.1 b) (-(w0 p.1 p.2 / 4)) + incr (ep p.2 a) (ep p.2 b) (w0 p.1 p.2 / 4))

omit [DecidableEq E] [Fintype E] in
/-- **The increments of Sec. 5.2.3 assemble the quadratic form of `B⁰`**:
`uᵀ B⁰ v = Σ_{(I,J) ∈ S} w⁰_IJ (u_I − u_J)(v_I − v_J)`.

[formula: rc-low-order-gram] -/
theorem gramB0_form (u v : V → ℝ) : u ⬝ᵥ gramB0 S w0 ep *ᵥ v = formB0 S w0 ep u v := by
  simp only [gramB0, dot_sum_mulVec, dot_add_mulVec, dot_incr, formB0, Eop]
  refine Finset.sum_congr rfl fun p _ => ?_
  simp only [Fin.sum_univ_two]
  ring

omit [DecidableEq E] [Fintype E] in
/-- **`B⁰` is symmetric positive semidefinite** for nonnegative weights.

[formula: rc-gram-psd-kernel] -/
theorem gramB0_posSemidef (hw : ∀ p ∈ S, 0 ≤ w0 p.1 p.2) : (gramB0 S w0 ep).PosSemidef :=
  posSemidef_of_form (fun u v => by
      rw [gramB0_form, gramB0_form, formB0, formB0]
      exact Finset.sum_congr rfl fun p _ => by ring)
    (fun u => by
      rw [gramB0_form]
      exact Finset.sum_nonneg fun p hp => mul_nonneg (hw p hp) (mul_self_nonneg _))

omit [DecidableEq E] [Fintype E] in
/-- **The kernel of `B⁰`**: `B⁰ u = 0` iff `u_I = u_J` (edge averages) for every `(I, J) ∈ S` with
`w⁰_IJ > 0`.

[formula: rc-gram-psd-kernel] -/
theorem gramB0_ker (hw : ∀ p ∈ S, 0 ≤ w0 p.1 p.2) (u : V → ℝ) :
    gramB0 S w0 ep *ᵥ u = 0 ↔ ∀ p ∈ S, 0 < w0 p.1 p.2 → Eop ep u p.1 = Eop ep u p.2 := by
  rw [mulVec_eq_zero_iff_form (gramB0_posSemidef S w0 ep hw), gramB0_form, formB0,
    Finset.sum_eq_zero_iff_of_nonneg fun p hp => mul_nonneg (hw p hp) (mul_self_nonneg _)]
  refine forall₂_congr fun p hp => ?_
  rw [mul_eq_zero, mul_self_eq_zero, sub_eq_zero]
  constructor
  · rintro (h | h) hpos
    · exact absurd h hpos.ne'
    · exact h
  · intro h
    rcases (hw p hp).eq_or_lt with h0 | h0
    · exact Or.inl h0.symm
    · exact Or.inr (h h0)

omit [DecidableEq E] [Fintype E] in
/-- **Global constants are in `ker B⁰`.**

[formula: rc-gram-psd-kernel] -/
theorem gramB0_mulVec_const (hw : ∀ p ∈ S, 0 ≤ w0 p.1 p.2) (c : ℝ) :
    gramB0 S w0 ep *ᵥ (fun _ => c) = 0 :=
  (gramB0_ker S w0 ep hw _).mpr fun _ _ _ => by simp [Eop]

omit [DecidableEq E] [Fintype E] in
/-- **`ker E ⊆ ker B⁰`**: a function with `u_{i₁} = −u_{i₂}` on every edge (alternating values on a
bipartite component: an open curve, a tree, a loop with an even number of vertices) has `E u = 0`,
hence `B⁰ u = 0`. So `ker B⁰` is in general larger than the constants, even when every pair of
edges is coupled (SPEC §4.5 claims `dim ker B⁰ = 1`; `ker A` is still the constants, see
`gram_ker_eq_const`).

[formula: rc-gram-psd-kernel] -/
theorem gramB0_mulVec_eq_zero_of_alternating (hw : ∀ p ∈ S, 0 ≤ w0 p.1 p.2) {u : V → ℝ}
    (hu : ∀ I, u (ep I 0) = -u (ep I 1)) : gramB0 S w0 ep *ᵥ u = 0 :=
  (gramB0_ker S w0 ep hw _).mpr fun p _ _ => by simp [Eop, hu]

/-- Two edges are linked through the low-order kernel. -/
def Linked (I J : E) : Prop := ((I, J) ∈ S ∧ 0 < w0 I J) ∨ ((J, I) ∈ S ∧ 0 < w0 J I)

omit [DecidableEq E] [Fintype E] in
/-- **If the edges are connected through `w⁰ > 0`, `B⁰ u = 0` iff all edge averages `(E u)_I`
are equal.**

[formula: rc-gram-psd-kernel] -/
theorem gramB0_ker_connected (hw : ∀ p ∈ S, 0 ≤ w0 p.1 p.2)
    (hconn : ∀ I J, Relation.ReflTransGen (Linked S w0) I J) (u : V → ℝ) :
    gramB0 S w0 ep *ᵥ u = 0 ↔ ∀ I J, Eop ep u I = Eop ep u J := by
  rw [gramB0_ker S w0 ep hw]
  constructor
  · intro h I J
    induction hconn I J with
    | refl => rfl
    | tail _ hlink ih =>
      rw [ih]
      rcases hlink with ⟨hm, hp⟩ | ⟨hm, hp⟩
      · exact h _ hm hp
      · exact (h _ hm hp).symm
  · intro h p _ _
    exact h _ _

end lowOrder

/-! ## The metric `A = B + B⁰` and its kernel -/

section metric

variable (S : Finset (E × E)) (w w0 : E → E → ℝ) (ℓ : E → ℝ) (T : E → V3) (ep : E → Fin 2 → V)

omit [DecidableEq E] [Fintype E] in
/-- **`ker A` is the constants** (SPEC §4.5, D5) when (i) every edge has positive length and is
coupled with `w > 0` to an edge with a linearly independent tangent, (ii) the edges are connected
through `w⁰ > 0`, (iii) every vertex is an endpoint. Then `(B + B⁰) u = 0` iff `u` is constant.

[formula: rc-gram-psd-kernel] -/
theorem gram_ker_eq_const (hT : ∀ I, T I ⬝ᵥ T I = 1) (hw : ∀ p ∈ S, 0 ≤ w p.1 p.2)
    (hw0 : ∀ p ∈ S, 0 ≤ w0 p.1 p.2) (hℓ : ∀ I, ℓ I ≠ 0)
    (hgen : ∀ I, ∃ J, (I, J) ∈ S ∧ 0 < w I J ∧ ∀ a b : ℝ, a • T I = b • T J → a = 0 ∧ b = 0)
    (hconn : ∀ I J, Relation.ReflTransGen (Linked S w0) I J)
    (hcov : ∀ v, ∃ I a, ep I a = v) (u : V → ℝ) :
    (gramB S w ℓ T ep + gramB0 S w0 ep) *ᵥ u = 0 ↔ ∃ c, ∀ v, u v = c := by
  rw [ker_add_posSemidef (gramB_posSemidef S w ℓ T ep hT hw) (gramB0_posSemidef S w0 ep hw0)]
  constructor
  · rintro ⟨hB, hB0⟩
    rw [gramB_ker S w ℓ T ep hT hw] at hB
    rw [gramB0_ker_connected S w0 ep hw0 hconn] at hB0
    -- every edge has equal end values
    have hedge : ∀ I, u (ep I 0) = u (ep I 1) := by
      intro I
      obtain ⟨J, hJ, hpos, hind⟩ := hgen I
      have := hB (I, J) hJ hpos
      simp only [Dop] at this
      have h0 := (hind _ _ this).1
      rw [div_eq_zero_iff, sub_eq_zero] at h0
      rcases h0 with h0 | h0
      · exact h0.symm
      · exact absurd h0 (hℓ I)
    have hall : ∀ v w, u v = u w := by
      intro v w
      obtain ⟨I, a, rfl⟩ := hcov v
      obtain ⟨J, b, rfl⟩ := hcov w
      have hE := hB0 I J
      simp only [Eop, ← hedge] at hE
      have hI : u (ep I a) = u (ep I 0) := by
        fin_cases a
        · rfl
        · exact (hedge I).symm
      have hJ : u (ep J b) = u (ep J 0) := by
        fin_cases b
        · rfl
        · exact (hedge J).symm
      rw [hI, hJ]
      linarith
    by_cases hV : Nonempty V
    · obtain ⟨v₀⟩ := hV
      exact ⟨u v₀, fun v => hall v v₀⟩
    · exact ⟨0, fun v => absurd ⟨v⟩ hV⟩
  · rintro ⟨c, hc⟩
    have hu : u = fun _ => c := funext hc
    subst hu
    exact ⟨gramB_mulVec_eq_zero_of_edge_const S w ℓ T ep hT hw (fun _ => rfl),
      gramB0_mulVec_const S w0 ep hw0 c⟩

omit [Fintype E] [DecidableEq E] in
/-- `(A ⊗ I₃) x` at `(v, c)` is `A` applied to the `c`-th coordinates. -/
theorem kron_one_mulVec {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ)
    (x : ι × Fin 3 → ℝ) (v : ι) (c : Fin 3) :
    ((A ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)) *ᵥ x) (v, c) = ∑ w, A v w * x (w, c) := by
  simp only [mulVec, dotProduct, kroneckerMap_apply, Fintype.sum_prod_type, one_apply]
  refine Finset.sum_congr rfl fun w _ => ?_
  simp [mul_ite, ite_mul]

omit [Fintype E] [DecidableEq E] in
/-- `xᵀ (L ⊗ I₃) y = Σ_IJ L_IJ ⟨x_I, y_J⟩`. -/
theorem dot_kron_one {ι : Type*} [Fintype ι] [DecidableEq ι] (L : Matrix ι ι ℝ)
    (x y : ι × Fin 3 → ℝ) :
    x ⬝ᵥ (L ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)) *ᵥ y =
      ∑ I, ∑ J, L I J * ∑ c, x (I, c) * y (J, c) := by
  rw [dotProduct, Fintype.sum_prod_type]
  simp only [kron_one_mulVec, Finset.mul_sum]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun J _ => Finset.sum_congr rfl fun c _ => by ring

omit [DecidableEq E] in
/-- **`Ā = A ⊗ I₃` acts coordinatewise** on the interleaved vector (D29), so `Ā x = 0` iff `x` is
a translation `x_v = t` whenever `ker A` is the constants.

[formula: rc-gram-psd-kernel] -/
theorem kron_ker_translations (A : Matrix V V ℝ)
    (hA : ∀ u : V → ℝ, A *ᵥ u = 0 ↔ ∃ c, ∀ v, u v = c) (x : V × Fin 3 → ℝ) :
    (A ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)) *ᵥ x = 0 ↔ ∃ t : V3, ∀ v c, x (v, c) = t c := by
  have hk : ∀ v c, ((A ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)) *ᵥ x) (v, c) =
      (A *ᵥ fun w => x (w, c)) v := by
    intro v c
    rw [kron_one_mulVec]
    rfl
  have key : (A ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)) *ᵥ x = 0 ↔
      ∀ c, (A *ᵥ fun w => x (w, c)) = 0 := by
    constructor
    · intro h c; funext v; rw [← hk]; simp [h]
    · intro h; funext ⟨v, c⟩; rw [hk, h c]; rfl
  rw [key]
  constructor
  · intro h
    choose t ht using fun c => (hA _).mp (h c)
    exact ⟨t, fun v c => ht c v⟩
  · rintro ⟨t, ht⟩ c
    exact (hA _).mpr ⟨t c, fun v => ht v c⟩

end metric

/-! ## The kernel-matrix form (Eq. 31, 32) -/

section kernelMatrix

/-- `L_K = diag(K𝟙) − K`. -/
def lapK (K : Matrix E E ℝ) : Matrix E E ℝ := diagonal (fun I => ∑ J, K I J) - K

omit [Fintype V] [DecidableEq V] in
/-- **The identity behind Eq. 31**: for symmetric `K`,
`Σ_IJ (diag(K𝟙) − K)_IJ ⟨X_I, Y_J⟩ = ½ Σ_IJ K_IJ ⟨X_I − X_J, Y_I − Y_J⟩`.

[formula: rc-kernel-matrix-form] -/
theorem lapK_bilin {K : Matrix E E ℝ} (hK : Kᵀ = K) (X Y : E → V3) :
    ∑ I, ∑ J, lapK K I J * (X I ⬝ᵥ Y J) =
      (1 / 2) * ∑ I, ∑ J, K I J * ((X I - X J) ⬝ᵥ (Y I - Y J)) := by
  have hs : ∀ I J, K J I = K I J := fun I J => by
    rw [← transpose_apply K I J, hK]
  have e1 : ∑ I, ∑ J, lapK K I J * (X I ⬝ᵥ Y J) =
      ∑ I, ∑ J, K I J * (X I ⬝ᵥ Y I) - ∑ I, ∑ J, K I J * (X I ⬝ᵥ Y J) := by
    simp only [lapK, Matrix.sub_apply, diagonal_apply, sub_mul, Finset.sum_sub_distrib, ite_mul, zero_mul,
      Finset.sum_ite_eq, Finset.mem_univ, ite_true, Finset.sum_mul]
  have e2 : ∑ I, ∑ J, K I J * (X J ⬝ᵥ Y J) = ∑ I, ∑ J, K I J * (X I ⬝ᵥ Y I) := by
    rw [Finset.sum_comm]
    simp only [hs]
  have e3 : ∑ I, ∑ J, K I J * (X J ⬝ᵥ Y I) = ∑ I, ∑ J, K I J * (X I ⬝ᵥ Y J) := by
    rw [Finset.sum_comm]
    simp only [hs]
  rw [e1]
  simp only [sub_dotProduct, dotProduct_sub, mul_sub, Finset.sum_sub_distrib,
    e2, e3]
  have e4 : ∑ I, ∑ J, K I J * (X J ⬝ᵥ Y I) = ∑ I, ∑ J, K I J * (X I ⬝ᵥ Y J) := e3
  linarith [e2, e4]

omit [Fintype V] [DecidableEq V] in
/-- The scalar version (for `B⁰`). -/
theorem lapK_bilin_scalar {K : Matrix E E ℝ} (hK : Kᵀ = K) (X Y : E → ℝ) :
    ∑ I, ∑ J, lapK K I J * (X I * Y J) =
      (1 / 2) * ∑ I, ∑ J, K I J * ((X I - X J) * (Y I - Y J)) := by
  have h := lapK_bilin hK (fun I => fun _ => X I / 3) (fun I => fun _ => Y I)
  simp only [dotProduct, Fin.sum_univ_three, Pi.sub_apply] at h
  have l : ∀ a b : ℝ, a / 3 * b + a / 3 * b + a / 3 * b = a * b := fun a b => by ring
  have l2 : ∀ a b c d : ℝ, (a / 3 - b / 3) * (c - d) + (a / 3 - b / 3) * (c - d) +
      (a / 3 - b / 3) * (c - d) = (a - b) * (c - d) := fun a b c d => by ring
  simpa only [l, l2] using h

/-- The weights extended by `0` off the pair set: `W_IJ = w_IJ` for `(I, J) ∈ S`, else `0`. -/
def extW (S : Finset (E × E)) (w : E → E → ℝ) : Matrix E E ℝ :=
  fun I J => if (I, J) ∈ S then w I J else 0

omit [Fintype E] [Fintype V] [DecidableEq V] in
theorem symK_transpose (S : Finset (E × E)) (w : E → E → ℝ) :
    (extW S w + (extW S w)ᵀ)ᵀ = extW S w + (extW S w)ᵀ := by
  rw [transpose_add, transpose_transpose, add_comm]

omit [Fintype V] [DecidableEq V] in
/-- A sum over the pair set is the double sum against the extended weights. -/
theorem sum_pairs_eq (S : Finset (E × E)) (w : E → E → ℝ) (f : E → E → ℝ) :
    ∑ p ∈ S, w p.1 p.2 * f p.1 p.2 = ∑ I, ∑ J, extW S w I J * f I J := by
  simp only [extW, ite_mul, zero_mul]
  rw [← Fintype.sum_prod_type' (f := fun I J => if (I, J) ∈ S then w I J * f I J else 0)]
  rw [Finset.sum_ite_mem, Finset.univ_inter]

omit [Fintype V] [DecidableEq V] in
/-- For `f` symmetric, `Σ W f = ½ Σ (W + Wᵀ) f`. -/
theorem sum_extW_symm (S : Finset (E × E)) (w : E → E → ℝ) {f : E → E → ℝ}
    (hf : ∀ I J, f J I = f I J) :
    ∑ I, ∑ J, extW S w I J * f I J =
      (1 / 2) * ∑ I, ∑ J, (extW S w + (extW S w)ᵀ) I J * f I J := by
  have : ∑ I, ∑ J, (extW S w)ᵀ I J * f I J = ∑ I, ∑ J, extW S w I J * f I J := by
    simp only [transpose_apply]
    rw [Finset.sum_comm]
    simp only [hf]
  simp only [Matrix.add_apply, add_mul, Finset.sum_add_distrib, this]
  ring

variable (S : Finset (E × E)) (w w0 : E → E → ℝ) (ℓ : E → ℝ) (T : E → V3) (ep : E → Fin 2 → V)

/-- **Eq. 31 for `B`**: with `K = W + Wᵀ`,
`uᵀ B v = Σ_IJ (diag(K𝟙) − K)_IJ ⟨D_I u, D_J v⟩` for unit tangents.

[formula: rc-kernel-matrix-form] -/
theorem gramB_kernel_form (hT : ∀ I, T I ⬝ᵥ T I = 1) (u v : V → ℝ) :
    u ⬝ᵥ gramB S w ℓ T ep *ᵥ v =
      ∑ I, ∑ J, lapK (extW S w + (extW S w)ᵀ) I J * (Dop ℓ T ep u I ⬝ᵥ Dop ℓ T ep v J) := by
  rw [gramB_form S w ℓ T ep hT, lapK_bilin (symK_transpose S w), formB,
    sum_pairs_eq S w (fun I J => (Dop ℓ T ep u I - Dop ℓ T ep u J) ⬝ᵥ
      (Dop ℓ T ep v I - Dop ℓ T ep v J)),
    sum_extW_symm S w (fun I J => by rw [← neg_sub (Dop ℓ T ep u I), ← neg_sub (Dop ℓ T ep v I),
      neg_dotProduct, dotProduct_neg, neg_neg])]

/-- **Eq. 31 for `B⁰`**: with `K⁰ = W⁰ + W⁰ᵀ`, `uᵀ B⁰ v = Σ_IJ (diag(K⁰𝟙) − K⁰)_IJ u_I v_J`.

[formula: rc-kernel-matrix-form] -/
theorem gramB0_kernel_form (u v : V → ℝ) :
    u ⬝ᵥ gramB0 S w0 ep *ᵥ v =
      ∑ I, ∑ J, lapK (extW S w0 + (extW S w0)ᵀ) I J * (Eop ep u I * Eop ep v J) := by
  rw [gramB0_form, lapK_bilin_scalar (symK_transpose S w0), formB0,
    sum_pairs_eq S w0 (fun I J => (Eop ep u I - Eop ep u J) * (Eop ep v I - Eop ep v J)),
    sum_extW_symm S w0 (fun I J => by ring)]

/-- **Eq. 31 as matrices**: `B = Dᵀ ((diag(K𝟙) − K) ⊗ I₃) D` (B.3.4 with `D` for `E`).

[formula: rc-kernel-matrix-form] -/
theorem gramB_eq_kron (hT : ∀ I, T I ⬝ᵥ T I = 1) :
    gramB S w ℓ T ep = (Dmat ℓ T ep)ᵀ *
      ((lapK (extW S w + (extW S w)ᵀ)) ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)) * Dmat ℓ T ep := by
  refine matrix_eq_of_form fun u v => ?_
  rw [gramB_kernel_form S w ℓ T ep hT, ← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec,
    vecMul_transpose]
  have hD : ∀ x : V → ℝ, Dmat ℓ T ep *ᵥ x = fun Ic => Dop ℓ T ep x Ic.1 Ic.2 := fun x => by
    funext Ic; exact (derivative_operator ℓ T ep x Ic.1).1 Ic.2
  rw [hD, hD, dot_kron_one]
  rfl

/-- **Eq. 31 as matrices**: `B⁰ = Eᵀ (diag(K⁰𝟙) − K⁰) E`.

[formula: rc-kernel-matrix-form] -/
theorem gramB0_eq_kernel :
    gramB0 S w0 ep = (Emat ep)ᵀ * lapK (extW S w0 + (extW S w0)ᵀ) * Emat ep := by
  refine matrix_eq_of_form fun u v => ?_
  rw [gramB0_kernel_form, ← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec, vecMul_transpose]
  have hE : ∀ x : V → ℝ, Emat ep *ᵥ x = fun I => Eop ep x I := fun x => by
    funext I
    simp only [Emat, mulVec, dotProduct, Eop, Finset.sum_mul]
    rw [Finset.sum_comm]
    simp only [ite_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ, ite_true, Fin.sum_univ_two]
    ring
  rw [hE, hE]
  simp only [dotProduct, mulVec, Finset.mul_sum]
  refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun J _ => by ring

end kernelMatrix

/-! ## The weights of the two quadratures and the scaling -/

section weights

variable (σ : ℝ) (γ : V → V3) (ep : E → Fin 2 → V)

/-- Trapezoid weight of `B`: `w_IJ = ¼ ℓ_I ℓ_J Σ_{i∈I} Σ_{j∈J} |γ_i − γ_j|^{−(2σ+1)}` (Sec. 5.2.2).

[formula: rc-high-order-gram] -/
def wTrap (I J : E) : ℝ :=
  (1 / 4) * len γ ep I * len γ ep J *
    ∑ a, ∑ b, vnorm (γ (ep I a) - γ (ep J b)) ^ (-(2 * σ + 1))

/-- Midpoint weight of `B`: `w_IJ = ℓ_I ℓ_J |x_I − x_J|^{−(2σ+1)}` (reference, Eq. 32; D2).

[formula: rc-high-order-gram] -/
def wMid (I J : E) : ℝ :=
  len γ ep I * len γ ep J * vnorm (mid γ ep I - mid γ ep J) ^ (-(2 * σ + 1))

/-- Trapezoid weight of `B⁰`: `w⁰_IJ = ¼ ℓ_I ℓ_J Σ_{i∈I} Σ_{j∈J} k²_{2σ+5}(γ_i, γ_j, T_I)`.

[formula: rc-low-order-gram] -/
def w0Trap (I J : E) : ℝ :=
  (1 / 4) * len γ ep I * len γ ep J *
    ∑ a, ∑ b, kernel 2 (2 * σ + 5) (γ (ep I a)) (γ (ep J b)) (tangent γ ep I)

/-- Midpoint weight of `B⁰`: `w⁰_IJ = ℓ_I ℓ_J k²_{2σ+5}(x_I, x_J, T_I)`.

[formula: rc-low-order-gram] -/
def w0Mid (I J : E) : ℝ :=
  len γ ep I * len γ ep J * kernel 2 (2 * σ + 5) (mid γ ep I) (mid γ ep J) (tangent γ ep I)

omit [Fintype V] [Fintype E] [DecidableEq V] [DecidableEq E] in
/-- **The trapezoid weight of `B⁰` as printed** (Sec. 5.2.3: `k²₄ / |γ_i − γ_j|^{2σ+1}`) equals the
`k²_{2σ+5}` form when no two vertices of the pair coincide (D2).

[formula: rc-kernel-matrix-form] -/
theorem w0Trap_eq_printed {I J : E} (h : ∀ a b, γ (ep I a) ≠ γ (ep J b)) :
    w0Trap σ γ ep I J = (1 / 4) * len γ ep I * len γ ep J *
      ∑ a, ∑ b, kernel 2 4 (γ (ep I a)) (γ (ep J b)) (tangent γ ep I) *
        vnorm (γ (ep I a) - γ (ep J b)) ^ (-(2 * σ + 1)) := by
  simp only [w0Trap, kernel_sq_quadrature σ (h _ _)]

omit [Fintype V] [Fintype E] [DecidableEq V] [DecidableEq E] in
/-- **Eq. 32 is the midpoint quadrature** (D2): `K⁰_IJ = w⁰_IJ + w⁰_JI =
(k²_{2σ+5}(x_I, x_J, T_I) + k²_{2σ+5}(x_J, x_I, T_J)) ℓ_I ℓ_J`, and `K^B_IJ = w_IJ + w_JI =
(k⁰_{2σ+1}(x_I, x_J, ·) + k⁰_{2σ+1}(x_J, x_I, ·)) ℓ_I ℓ_J = 2 ℓ_I ℓ_J |x_I − x_J|^{−(2σ+1)}`.

[formula: rc-kernel-matrix-form] -/
theorem eq32_midpoint (I J : E) :
    w0Mid σ γ ep I J + w0Mid σ γ ep J I =
        (kernel 2 (2 * σ + 5) (mid γ ep I) (mid γ ep J) (tangent γ ep I)
          + kernel 2 (2 * σ + 5) (mid γ ep J) (mid γ ep I) (tangent γ ep J))
          * len γ ep I * len γ ep J ∧
      wMid σ γ ep I J + wMid σ γ ep J I =
        (kernel 0 (2 * σ + 1) (mid γ ep I) (mid γ ep J) (tangent γ ep I)
          + kernel 0 (2 * σ + 1) (mid γ ep J) (mid γ ep I) (tangent γ ep J))
          * len γ ep I * len γ ep J ∧
      wMid σ γ ep I J + wMid σ γ ep J I =
        2 * len γ ep I * len γ ep J * vnorm (mid γ ep I - mid γ ep J) ^ (-(2 * σ + 1)) := by
  refine ⟨by simp only [w0Mid]; ring, ?_, ?_⟩
  · simp only [wMid, kernel_zero_alpha, neg_add_rev, vnorm_sub_comm (mid γ ep J)]
    ring
  · simp only [wMid, vnorm_sub_comm (mid γ ep J)]
    ring

omit [Fintype V] [Fintype E] [DecidableEq V] [DecidableEq E] in
theorem rpow_scale {c : ℝ} (hc : 0 < c) (x y : ℝ) (hx : 0 ≤ x) :
    (c * x) ^ y = c ^ y * x ^ y := Real.mul_rpow hc.le hx

omit [Fintype V] [Fintype E] [DecidableEq V] [DecidableEq E] in
theorem mid_smul (c : ℝ) (I : E) : mid (c • γ) ep I = c • mid γ ep I := by
  ext i
  simp [mid]
  ring

omit [Fintype V] [Fintype E] [DecidableEq V] [DecidableEq E] in
/-- **Scaling of the weights of `B`**: `w(cγ) = c^{1−2σ} w(γ)` for both quadratures, `c > 0`.

[formula: rc-gram-scaling] -/
theorem w_smul {c : ℝ} (hc : 0 < c) (I J : E) :
    wTrap σ (c • γ) ep I J = c ^ (1 - 2 * σ) * wTrap σ γ ep I J ∧
      wMid σ (c • γ) ep I J = c ^ (1 - 2 * σ) * wMid σ γ ep I J := by
  have hp : c ^ (1 - 2 * σ) = c * c * c ^ (-(2 * σ + 1)) := by
    rw [show 1 - 2 * σ = (2 : ℝ) + (-(2 * σ + 1)) by ring, Real.rpow_add hc, Real.rpow_two]
    ring
  constructor
  · simp only [wTrap, len_smul γ ep hc, Pi.smul_apply, ← smul_sub, vnorm_smul, abs_of_pos hc,
      rpow_scale hc _ _ (vnorm_nonneg _), ← Finset.mul_sum, hp]
    ring
  · simp only [wMid, len_smul γ ep hc, mid_smul, ← smul_sub, vnorm_smul, abs_of_pos hc,
      rpow_scale hc _ _ (vnorm_nonneg _), hp]
    ring

omit [Fintype V] [Fintype E] [DecidableEq V] [DecidableEq E] in
/-- **Scaling of the weights of `B⁰`**: `w⁰(cγ) = c^{−(2σ+1)} w⁰(γ)` for both quadratures,
`c > 0` (`ℓ_I ℓ_J` gives `c²`, `k²_{2σ+5}` gives `c^{2−(2σ+5)}`). (SPEC §4.5 printed
`c^{1−2σ}` here; its conclusion about `B⁰` is right.)

[formula: rc-gram-scaling] -/
theorem w0_smul {c : ℝ} (hc : 0 < c) (I J : E) :
    w0Trap σ (c • γ) ep I J = c ^ (-(2 * σ + 1)) * w0Trap σ γ ep I J ∧
      w0Mid σ (c • γ) ep I J = c ^ (-(2 * σ + 1)) * w0Mid σ γ ep I J := by
  have hp : c ^ (-(2 * σ + 1)) = c * c * c ^ ((2 : ℝ) - (2 * σ + 5)) := by
    rw [show -(2 * σ + 1) = (2 : ℝ) + ((2 : ℝ) - (2 * σ + 5)) by ring, Real.rpow_add hc,
      Real.rpow_two]
    ring
  constructor
  · simp only [w0Trap, len_smul γ ep hc, tangent_smul γ ep hc, Pi.smul_apply, kernel_smul hc,
      ← Finset.mul_sum, hp]
    ring
  · simp only [w0Mid, len_smul γ ep hc, tangent_smul γ ep hc, mid_smul, kernel_smul hc, hp]
    ring

variable (S : Finset (E × E))

omit [Fintype V] [DecidableEq E] [Fintype E] in
/-- `B` with weights `s·w` and lengths `c·ℓ` is `(s/c²) B`. -/
theorem gramB_scale (w : E → E → ℝ) (ℓ : E → ℝ) (T : E → V3) {s c : ℝ} (hc : c ≠ 0) :
    gramB S (fun I J => s * w I J) (fun I => c * ℓ I) T ep = (s / c ^ 2) • gramB S w ℓ T ep := by
  simp only [gramB, Finset.smul_sum, smul_add]
  refine Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun a _ =>
    Finset.sum_congr rfl fun b _ => ?_
  have hs : ∀ (i j : V) (x : ℝ), (s / c ^ 2) • incr i j x = incr i j (s / c ^ 2 * x) := by
    intro i j x; ext y z; simp [incr]
  simp only [hs]
  congr 1
  · congr 1
    · congr 1
      · congr 1; field_simp
      · congr 1; field_simp
    · congr 1; field_simp
  · congr 1; field_simp

omit [DecidableEq E] [Fintype V] [Fintype E] in
/-- **Scaling of `B`**: `B(cγ) = c^{−(2σ+1)} B(γ)` for `c > 0`, for any weight with
`w(cγ) = c^{1−2σ} w(γ)` (both quadratures, `w_smul`).

[formula: rc-gram-scaling] -/
theorem gramB_smul {c : ℝ} (hc : 0 < c) (w : (V → V3) → E → E → ℝ)
    (hw : ∀ I J, w (c • γ) I J = c ^ (1 - 2 * σ) * w γ I J) :
    gramB S (w (c • γ)) (len (c • γ) ep) (tangent (c • γ) ep) ep =
      c ^ (-(2 * σ + 1)) • gramB S (w γ) (len γ ep) (tangent γ ep) ep := by
  have h1 : w (c • γ) = fun I J => c ^ (1 - 2 * σ) * w γ I J := funext fun I => funext (hw I)
  have h2 : len (c • γ) ep = fun I => c * len γ ep I := funext (len_smul γ ep hc)
  have h3 : tangent (c • γ) ep = tangent γ ep := funext (tangent_smul γ ep hc)
  rw [h1, h2, h3, gramB_scale ep S _ _ _ hc.ne']
  congr 1
  rw [← Real.rpow_natCast c 2, ← Real.rpow_sub hc]
  congr 1
  push_cast
  ring

omit [Fintype V] [Fintype E] [DecidableEq E] in
/-- **Scaling of `B⁰`**: `B⁰(cγ) = c^{−(2σ+1)} B⁰(γ)` for `c > 0`, for any weight with
`w⁰(cγ) = c^{−(2σ+1)} w⁰(γ)` (both quadratures, `w0_smul`).

[formula: rc-gram-scaling] -/
theorem gramB0_smul {c : ℝ} (w0 : (V → V3) → E → E → ℝ)
    (hw : ∀ I J, w0 (c • γ) I J = c ^ (-(2 * σ + 1)) * w0 γ I J) :
    gramB0 S (w0 (c • γ)) ep = c ^ (-(2 * σ + 1)) • gramB0 S (w0 γ) ep := by
  simp only [gramB0, hw, Finset.smul_sum, smul_add]
  have hs : ∀ (i j : V) (x : ℝ), c ^ (-(2 * σ + 1)) • incr i j x =
      incr i j (c ^ (-(2 * σ + 1)) * x) := by
    intro i j x; ext y z; simp [incr]
  simp only [hs]
  refine Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun a _ =>
    Finset.sum_congr rfl fun b _ => ?_
  congr 1
  · congr 1
    · congr 1
      · congr 1; ring
      · congr 1; ring
    · congr 1; ring
  · congr 1; ring

end weights

/-! ## The comparison metrics L², H¹, H² (SPEC §4.6) -/

section classic

variable (ℓ : E → ℝ) (ep : E → Fin 2 → V)

/-- The H¹ matrix (weak Laplacian of the curve), by increments: for every edge `I`,
`L_{i₁i₁} += 1/ℓ_I`, `L_{i₂i₂} += 1/ℓ_I`, `L_{i₁i₂} −= 1/ℓ_I`, `L_{i₂i₁} −= 1/ℓ_I`.

[formula: rc-classic-metrics] -/
def lapH1 : Matrix V V ℝ :=
  ∑ I, (incr (ep I 0) (ep I 0) (1 / ℓ I) + incr (ep I 1) (ep I 1) (1 / ℓ I)
    + incr (ep I 0) (ep I 1) (-(1 / ℓ I)) + incr (ep I 1) (ep I 0) (-(1 / ℓ I)))

/-- The H² matrix `L M⁻¹ L` with the lumped masses `M = diag(m)`.

[formula: rc-classic-metrics] -/
def metricH2 (m : V → ℝ) : Matrix V V ℝ := lapH1 ℓ ep * diagonal (fun v => (m v)⁻¹) * lapH1 ℓ ep

omit [DecidableEq E] in
/-- **The H¹ form**: `uᵀ L v = Σ_I (u_{i₂} − u_{i₁})(v_{i₂} − v_{i₁}) / ℓ_I`.

[formula: rc-classic-metrics] -/
theorem lapH1_form (u v : V → ℝ) :
    u ⬝ᵥ lapH1 ℓ ep *ᵥ v = ∑ I, (u (ep I 1) - u (ep I 0)) * (v (ep I 1) - v (ep I 0)) / ℓ I := by
  simp only [lapH1, dot_sum_mulVec, dot_add_mulVec, dot_incr]
  refine Finset.sum_congr rfl fun I _ => ?_
  ring

omit [DecidableEq E] in
theorem lapH1_symm (u v : V → ℝ) : u ⬝ᵥ lapH1 ℓ ep *ᵥ v = v ⬝ᵥ lapH1 ℓ ep *ᵥ u := by
  rw [lapH1_form, lapH1_form]
  exact Finset.sum_congr rfl fun I _ => by ring

omit [DecidableEq E] in
/-- **`L` is PSD and `ker L` = the functions constant along every edge** (the per-component
constants), for positive lengths.

[formula: rc-classic-metrics] -/
theorem lapH1_psd_ker (hℓ : ∀ I, 0 < ℓ I) :
    (lapH1 ℓ ep).PosSemidef ∧
      ∀ u, lapH1 ℓ ep *ᵥ u = 0 ↔ ∀ I, u (ep I 0) = u (ep I 1) := by
  have hn : ∀ u, 0 ≤ u ⬝ᵥ lapH1 ℓ ep *ᵥ u := fun u => by
    rw [lapH1_form]
    exact Finset.sum_nonneg fun I _ => div_nonneg (mul_self_nonneg _) (hℓ I).le
  have hpsd := posSemidef_of_form (lapH1_symm ℓ ep) hn
  refine ⟨hpsd, fun u => ?_⟩
  rw [mulVec_eq_zero_iff_form hpsd, lapH1_form,
    Finset.sum_eq_zero_iff_of_nonneg fun I _ => div_nonneg (mul_self_nonneg _) (hℓ I).le]
  refine forall_congr' fun I => ?_
  simp only [Finset.mem_univ, true_implies]
  rw [div_eq_zero_iff, or_iff_left (hℓ I).ne', mul_self_eq_zero, sub_eq_zero, eq_comm]

omit [DecidableEq E] in
/-- **`ker (L M⁻¹ L) = ker L`** for positive masses, and `uᵀ L M⁻¹ L u = Σ_v (L u)_v² / m_v`.

[formula: rc-classic-metrics] -/
theorem metricH2_ker (hℓ : ∀ I, 0 < ℓ I) {m : V → ℝ} (hm : ∀ v, 0 < m v) (u : V → ℝ) :
    u ⬝ᵥ metricH2 ℓ ep m *ᵥ u = ∑ v, (lapH1 ℓ ep *ᵥ u) v ^ 2 / m v ∧
      (metricH2 ℓ ep m *ᵥ u = 0 ↔ ∀ I, u (ep I 0) = u (ep I 1)) := by
  have hsymm : (lapH1 ℓ ep)ᵀ = lapH1 ℓ ep := transpose_eq_of_form (lapH1_symm ℓ ep)
  have hform : u ⬝ᵥ metricH2 ℓ ep m *ᵥ u = ∑ v, (lapH1 ℓ ep *ᵥ u) v ^ 2 / m v := by
    rw [metricH2, ← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec, ← mulVec_transpose, hsymm]
    simp only [dotProduct, mulVec_diagonal]
    exact Finset.sum_congr rfl fun v _ => by field_simp
  refine ⟨hform, ?_⟩
  rw [← (lapH1_psd_ker ℓ ep hℓ).2 u]
  constructor
  · intro h
    have h0 : ∑ v, (lapH1 ℓ ep *ᵥ u) v ^ 2 / m v = 0 := by rw [← hform, h, dotProduct_zero]
    rw [Finset.sum_eq_zero_iff_of_nonneg fun v _ => div_nonneg (sq_nonneg _) (hm v).le] at h0
    funext v
    have := h0 v (Finset.mem_univ v)
    rw [div_eq_zero_iff, or_iff_left (hm v).ne'] at this
    simpa using this
  · intro h
    rw [metricH2, ← mulVec_mulVec, ← mulVec_mulVec, h, mulVec_zero, mulVec_zero]

end classic

end RC

end DdgProofs
