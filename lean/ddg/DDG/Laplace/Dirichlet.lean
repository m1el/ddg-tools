import DDG.Laplace.Kernel
import DDG.Laplace.Cotan

/-!
# The Dirichlet problem and the harmonic extension for a square matrix

Linear algebra over a finite index set `V` split by a predicate `I` into the unknowns (`I v`) and
the rest (`B`). `L` is any square matrix unless a hypothesis says otherwise; for the cotan
Laplacian of a connected valid mesh with a nonempty `B` the hypotheses hold
(`TriMesh.existsUnique_harmonicExtension`, …). Used by the Variational Cuts flow
(`crates/ddg-variational-cuts/SPEC.md` §3.4 – §3.6, §3.8).

`IsDirichletSolution L I f u`: `u = 0` off `I` and `(L u)_i = f_i` on `I`, i.e. `u_B = 0`,
`L_II u_I = f_I` (`isDirichletSolution_iff_block`). `boundaryFlux L I f u` is `(L u)_i − f_i` at
`i ∉ I` and `0` on `I`.

* **`cut-yamabe-system`** (`f = −g`): `existsUnique_isDirichletSolution`: a unique solution when
  `L_II` is invertible (positive definite).
* **`cut-harmonic-extension`**: `existsUnique_dirichlet_data`: unique `σ` with `σ = s` off `O`
  and `(L σ)_i = 0` on `O` (`L_OO σ_O = −L_OC s_C`, `harmonic_iff_block`);
  `harmonicExtension_energy_le`: it minimizes `τᵀ L τ` among all `τ` with the same data, with
  equality only for `τ = σ` when `L_OO` is positive definite; the same on a mesh
  (`TriMesh.existsUnique_harmonicExtension`, `TriMesh.harmonicExtension_faceEnergy_le`).

(Split off `DdgProofs.Cuts.Dirichlet`, which keeps the oracle-facing part: the Yamabe system and
normal derivative, the gradient densities and the adjoints.)
-/

noncomputable section

namespace DDG

namespace Cuts

open Matrix

variable {V : Type*} [Fintype V]

/-! ## Block products -/

/-- If `u` vanishes off `I`, any row of `L u` only sees the columns in `I`. -/
theorem mulVec_of_zero_off {L : Matrix V V ℝ} {I : V → Prop} [DecidablePred I] {u : V → ℝ}
    (hu : ∀ v, ¬I v → u v = 0) (i : V) :
    (L *ᵥ u) i = ∑ j : {v // I v}, L i j * u j := by
  simp only [Matrix.mulVec, dotProduct]
  rw [← Finset.sum_subtype (Finset.univ.filter I) (by simp) (fun j => L i j * u j),
    Finset.sum_filter]
  refine Finset.sum_congr rfl fun j _ => ?_
  by_cases h : I j
  · simp [h]
  · simp [h, hu j h]

/-- … so at a row of `I` it is the product with the block `L_II`. -/
theorem mulVec_eq_blockII {L : Matrix V V ℝ} {I : V → Prop} [DecidablePred I] {u : V → ℝ}
    (hu : ∀ v, ¬I v → u v = 0) (i : {v // I v}) :
    (L *ᵥ u) i =
      ((L.submatrix (Subtype.val : {v // I v} → V) Subtype.val) *ᵥ fun j : {v // I v} => u j) i :=
  mulVec_of_zero_off hu i

/-- … and at a row outside `I` the product with the block `L_BI`. -/
theorem mulVec_eq_blockBI {L : Matrix V V ℝ} {I : V → Prop} [DecidablePred I] {u : V → ℝ}
    (hu : ∀ v, ¬I v → u v = 0) (i : {v // ¬I v}) :
    (L *ᵥ u) i =
      ((L.submatrix (Subtype.val : {v // ¬I v} → V) (Subtype.val : {v // I v} → V)) *ᵥ
        fun j : {v // I v} => u j) i :=
  mulVec_of_zero_off hu i

/-- Splitting a row of `L a` into the columns in `I` and the others. -/
theorem mulVec_split (L : Matrix V V ℝ) (I : V → Prop) [DecidablePred I] (a : V → ℝ) (i : V) :
    (L *ᵥ a) i = ∑ j : {v // I v}, L i j * a j + ∑ j : {v // ¬I v}, L i j * a j := by
  simp only [Matrix.mulVec, dotProduct]
  exact (Fintype.sum_subtype_add_sum_subtype I fun j => L i j * a j).symm

/-- For a symmetric matrix, `xᵀ L y = yᵀ L x`. -/
theorem dotProduct_mulVec_symm {L : Matrix V V ℝ} (hL : L.IsSymm) (x y : V → ℝ) :
    x ⬝ᵥ (L *ᵥ y) = y ⬝ᵥ (L *ᵥ x) := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hL.eq, dotProduct_comm]

/-- If the columns of `L` sum to zero, the entries of `L u` sum to zero. -/
theorem sum_mulVec_eq_zero {L : Matrix V V ℝ} (hL : (fun _ => (1 : ℝ)) ᵥ* L = 0) (u : V → ℝ) :
    ∑ v, (L *ᵥ u) v = 0 :=
  Solver.sum_eq_zero_of_mulVec_eq hL rfl

/-- The quadratic form of a vector that vanishes off `I` is the quadratic form of the block. -/
theorem dotProduct_mulVec_blockII {L : Matrix V V ℝ} {I : V → Prop} [DecidablePred I]
    {u : V → ℝ} (hu : ∀ v, ¬I v → u v = 0) :
    u ⬝ᵥ (L *ᵥ u) = (fun j : {v // I v} => u j) ⬝ᵥ
      ((L.submatrix (Subtype.val : {v // I v} → V) Subtype.val) *ᵥ fun j : {v // I v} => u j) := by
  have e1 : u ⬝ᵥ (L *ᵥ u) = ∑ i : {v // I v}, u i * (L *ᵥ u) i := by
    simp only [dotProduct]
    rw [← Finset.sum_subtype (Finset.univ.filter I) (by simp) (fun i => u i * (L *ᵥ u) i),
      Finset.sum_filter]
    refine Finset.sum_congr rfl fun i _ => ?_
    by_cases h : I i
    · simp [h]
    · simp [h, hu i h]
  rw [e1]
  simp only [dotProduct]
  exact Finset.sum_congr rfl fun i _ => by rw [mulVec_eq_blockII hu i]

/-! ## The Dirichlet problem -/

/-- `u` solves the Dirichlet problem with zero data: `u = 0` off `I` and `(L u)_i = f_i` on `I`.
The Yamabe system of SPEC §3.4 is the case `f = −g`; the adjoint systems of §3.6 are the cases
`f = r`. -/
structure IsDirichletSolution (L : Matrix V V ℝ) (I : V → Prop) (f u : V → ℝ) : Prop where
  zero_off : ∀ v, ¬I v → u v = 0
  eq_on : ∀ v, I v → (L *ᵥ u) v = f v

/-- The boundary flux of a Dirichlet solution: `(L u)_i − f_i` at `i ∉ I`, `0` on `I`. For the
Yamabe system (`f = −g`) it is the integrated outward normal derivative `v_i = (L_BI u_I)_i + g_i`
of SPEC §3.5; for an adjoint system (`f = r`) it is `q_i = (L_BI p_I)_i − r_i` of §3.6. -/
def boundaryFlux (L : Matrix V V ℝ) (I : V → Prop) [DecidablePred I] (f u : V → ℝ) : V → ℝ :=
  fun v => if I v then 0 else (L *ᵥ u) v - f v

section Dirichlet

variable {L : Matrix V V ℝ} {I : V → Prop} [DecidablePred I]

/-- **Block form**: `u_B = 0` and `L_II u_I = f_I`. -/
theorem isDirichletSolution_iff_block {f u : V → ℝ} :
    IsDirichletSolution L I f u ↔ (∀ v, ¬I v → u v = 0) ∧
      (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val) *ᵥ (fun j : {v // I v} => u j) =
        fun j : {v // I v} => f j := by
  constructor
  · rintro ⟨h0, h1⟩
    refine ⟨h0, funext fun i => ?_⟩
    rw [← mulVec_eq_blockII h0 i]
    exact h1 i i.2
  · rintro ⟨h0, h1⟩
    refine ⟨h0, fun v hv => ?_⟩
    rw [mulVec_eq_blockII h0 ⟨v, hv⟩, h1]

/-- **Existence and uniqueness.** If the block `L_II` is invertible there is exactly one solution
of the Dirichlet problem, for every right-hand side.

[formula: cut-yamabe-system] -/
theorem existsUnique_isDirichletSolution [DecidableEq V]
    (hdet : IsUnit (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val).det) (f : V → ℝ) :
    ∃! u : V → ℝ, IsDirichletSolution L I f u := by
  set A := L.submatrix (Subtype.val : {v // I v} → V) Subtype.val with hA
  set x : {v // I v} → ℝ := A⁻¹ *ᵥ (fun i : {v // I v} => f i) with hx
  set u : V → ℝ := fun v => if h : I v then x ⟨v, h⟩ else 0 with hu
  have huoff : ∀ v, ¬I v → u v = 0 := fun v hv => by simp [hu, hv]
  have hux : (fun j : {v // I v} => u j) = x := by funext j; simp [hu, j.2]
  refine ⟨u, ⟨huoff, fun v hv => ?_⟩, ?_⟩
  · rw [mulVec_eq_blockII huoff ⟨v, hv⟩, hux, hx, Matrix.mulVec_mulVec,
      Matrix.mul_nonsing_inv _ hdet, Matrix.one_mulVec]
  · rintro w ⟨hw0, hw1⟩
    have hAw : A *ᵥ (fun j : {v // I v} => w j) = fun i : {v // I v} => f i := by
      funext i
      rw [← mulVec_eq_blockII hw0 i]
      exact hw1 i i.2
    have hwx : (fun j : {v // I v} => w j) = x := by
      rw [hx, ← hAw, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet, Matrix.one_mulVec]
    funext v
    by_cases hv : I v
    · have h1 := congrFun hwx ⟨v, hv⟩
      have h2 := congrFun hux ⟨v, hv⟩
      simp only at h1 h2
      rw [h1, h2]
    · rw [hw0 v hv, huoff v hv]

/-- A positive definite block is invertible. -/
theorem isUnit_det_of_posDef [DecidableEq V]
    (hpd : (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val).PosDef) :
    IsUnit (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val).det :=
  (Matrix.isUnit_iff_isUnit_det _).mp hpd.isUnit

/-- **The partitioned system.** A Dirichlet solution satisfies `L u = f + (boundary flux)` at
every index. -/
theorem IsDirichletSolution.mulVec_eq {f u : V → ℝ} (h : IsDirichletSolution L I f u) :
    L *ᵥ u = f + boundaryFlux L I f u := by
  funext v
  by_cases hv : I v
  · simp [boundaryFlux, hv, h.eq_on v hv]
  · simp [boundaryFlux, hv]

/-- The boundary flux vanishes on `I`. -/
theorem boundaryFlux_of_mem (f u : V → ℝ) {v : V} (hv : I v) : boundaryFlux L I f u v = 0 := by
  simp [boundaryFlux, hv]

/-- **The boundary flux is determined**: it is the only vector that completes `L u = f + w`. -/
theorem IsDirichletSolution.flux_unique {f u w : V → ℝ} (h : IsDirichletSolution L I f u)
    (hw : L *ᵥ u = f + w) : w = boundaryFlux L I f u := by
  have := h.mulVec_eq
  rw [hw] at this
  exact add_left_cancel this

/-- **Flux balance.** If the columns of `L` sum to zero (a Laplacian), the boundary flux of a
Dirichlet solution sums to `−Σ f`. -/
theorem IsDirichletSolution.sum_boundaryFlux {f u : V → ℝ} (h : IsDirichletSolution L I f u)
    (hL : (fun _ => (1 : ℝ)) ᵥ* L = 0) : ∑ v, boundaryFlux L I f u v = -∑ v, f v := by
  have h1 := sum_mulVec_eq_zero hL u
  rw [h.mulVec_eq] at h1
  simp only [Pi.add_apply, Finset.sum_add_distrib] at h1
  linarith

/-- For a vector `w` that vanishes off `I`: `wᵀ L u = Σ w_i f_i`. -/
theorem IsDirichletSolution.dotProduct_mulVec {f u : V → ℝ} (h : IsDirichletSolution L I f u)
    {w : V → ℝ} (hw : ∀ v, ¬I v → w v = 0) : w ⬝ᵥ (L *ᵥ u) = ∑ v, w v * f v := by
  simp only [dotProduct]
  refine Finset.sum_congr rfl fun v _ => ?_
  by_cases hv : I v
  · rw [h.eq_on v hv]
  · rw [hw v hv, zero_mul, zero_mul]

/-- **Reciprocity** (the discrete adjoint identity). For a symmetric `L`, if `u` solves the
Dirichlet problem with right-hand side `f` and `p` the one with right-hand side `r`, then
`Σ p_i f_i = Σ u_i r_i`. -/
theorem IsDirichletSolution.reciprocity (hL : L.IsSymm) {f r u p : V → ℝ}
    (hu : IsDirichletSolution L I f u) (hp : IsDirichletSolution L I r p) :
    ∑ v, p v * f v = ∑ v, u v * r v := by
  rw [← hu.dotProduct_mulVec hp.zero_off, ← hp.dotProduct_mulVec hu.zero_off]
  exact dotProduct_mulVec_symm hL p u

omit [DecidablePred I] in
/-- Linearity: solutions add. -/
theorem IsDirichletSolution.add_smul {f g u w : V → ℝ} (hu : IsDirichletSolution L I f u)
    (hw : IsDirichletSolution L I g w) (t : ℝ) :
    IsDirichletSolution L I (f + t • g) (u + t • w) := by
  refine ⟨fun v hv => ?_, fun v hv => ?_⟩
  · simp [hu.zero_off v hv, hw.zero_off v hv]
  · rw [Matrix.mulVec_add, Matrix.mulVec_smul]
    simp [hu.eq_on v hv, hw.eq_on v hv]


/-! ## Harmonic extension -/

/-- **Dirichlet problem with data.** If `L_II` is invertible then for every source `φ` and all
prescribed values `s` there is exactly one `a` with `a = s` off `I` and `(L a)_i = φ_i` on `I`. -/
theorem existsUnique_dirichlet_data [DecidableEq V]
    (hdet : IsUnit (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val).det)
    (φ s : V → ℝ) :
    ∃! a : V → ℝ, (∀ v, ¬I v → a v = s v) ∧ ∀ v, I v → (L *ᵥ a) v = φ v := by
  set s₀ : V → ℝ := fun v => if I v then 0 else s v with hs₀
  obtain ⟨z, hz, hzu⟩ := existsUnique_isDirichletSolution (I := I) hdet
    (fun v => φ v - (L *ᵥ s₀) v)
  refine ⟨s₀ + z, ⟨fun v hv => ?_, fun v hv => ?_⟩, ?_⟩
  · simp [hs₀, hv, hz.zero_off v hv]
  · rw [Matrix.mulVec_add]
    simp only [Pi.add_apply]
    rw [hz.eq_on v hv]
    ring
  · rintro b ⟨hb1, hb2⟩
    have hbz : IsDirichletSolution L I (fun v => φ v - (L *ᵥ s₀) v) (b - s₀) := by
      refine ⟨fun v hv => by simp [hs₀, hv, hb1 v hv], fun v hv => ?_⟩
      rw [Matrix.mulVec_sub]
      simp only [Pi.sub_apply]
      rw [hb2 v hv]
    have h := hzu _ hbz
    funext v
    have hv := congrFun h v
    simp only [Pi.sub_apply] at hv
    simp only [Pi.add_apply]
    linarith

/-- **Harmonic extension** (SPEC §3.8): if `L_OO` is invertible (positive definite), the speed
`s` on the cut vertices has exactly one extension `σ` with `(L σ)_i = 0` at every unknown `i`.

[formula: cut-harmonic-extension] -/
theorem existsUnique_harmonicExtension [DecidableEq V]
    (hdet : IsUnit (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val).det) (s : V → ℝ) :
    ∃! σ : V → ℝ, (∀ v, ¬I v → σ v = s v) ∧ ∀ v, I v → (L *ᵥ σ) v = 0 :=
  existsUnique_dirichlet_data hdet (fun _ => 0) s

/-- **The system that is solved**: `(L σ)_i = 0` on `O` is `L_OO σ_O = −L_OC σ_C`.

[formula: cut-harmonic-extension] -/
theorem harmonic_iff_block {σ : V → ℝ} :
    (∀ v, I v → (L *ᵥ σ) v = 0) ↔
      (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val) *ᵥ (fun j : {v // I v} => σ j) =
        -((L.submatrix (Subtype.val : {v // I v} → V) (Subtype.val : {v // ¬I v} → V)) *ᵥ
          fun j : {v // ¬I v} => σ j) := by
  have key : ∀ i : {v // I v}, (L *ᵥ σ) i =
      ((L.submatrix (Subtype.val : {v // I v} → V) Subtype.val) *ᵥ
        (fun j : {v // I v} => σ j)) i +
      ((L.submatrix (Subtype.val : {v // I v} → V) (Subtype.val : {v // ¬I v} → V)) *ᵥ
        fun j : {v // ¬I v} => σ j) i := fun i => mulVec_split L I σ i
  constructor
  · intro h
    funext i
    have := key i
    rw [h i i.2] at this
    simp only [Pi.neg_apply]
    linarith
  · intro h v hv
    have h1 := key ⟨v, hv⟩
    have h2 := congrFun h ⟨v, hv⟩
    simp only [Pi.neg_apply] at h2
    simp only at h1
    linarith

/-- **Energy identity for a harmonic vector.** Let `L` be symmetric, `(L σ)_i = 0` on `I` and
`τ = σ` off `I`. Then `τᵀ L τ = σᵀ L σ + (τ − σ)ᵀ L (τ − σ)`. -/
theorem harmonic_energy_identity (hL : L.IsSymm) {σ τ : V → ℝ}
    (hσ : ∀ v, I v → (L *ᵥ σ) v = 0) (hτ : ∀ v, ¬I v → τ v = σ v) :
    τ ⬝ᵥ (L *ᵥ τ) = σ ⬝ᵥ (L *ᵥ σ) + (τ - σ) ⬝ᵥ (L *ᵥ (τ - σ)) := by
  have hd0 : (τ - σ) ⬝ᵥ (L *ᵥ σ) = 0 := by
    simp only [dotProduct, Pi.sub_apply]
    refine Finset.sum_eq_zero fun v _ => ?_
    by_cases hv : I v
    · rw [hσ v hv, mul_zero]
    · rw [hτ v hv, sub_self, zero_mul]
  have hsym : σ ⬝ᵥ (L *ᵥ (τ - σ)) = (τ - σ) ⬝ᵥ (L *ᵥ σ) := dotProduct_mulVec_symm hL _ _
  have hτ' : τ = σ + (τ - σ) := by
    funext v
    simp
  calc τ ⬝ᵥ (L *ᵥ τ) = (σ + (τ - σ)) ⬝ᵥ (L *ᵥ (σ + (τ - σ))) := by rw [← hτ']
    _ = σ ⬝ᵥ (L *ᵥ σ) + (τ - σ) ⬝ᵥ (L *ᵥ (τ - σ)) := by
        rw [Matrix.mulVec_add, add_dotProduct, dotProduct_add, dotProduct_add, hsym, hd0]
        ring

/-- **The harmonic extension minimizes the Dirichlet energy** among all vectors with the same
values off the unknowns (symmetric positive semidefinite `L`).

[formula: cut-harmonic-extension] -/
theorem harmonicExtension_energy_le (hpsd : L.PosSemidef) {σ τ : V → ℝ}
    (hσ : ∀ v, I v → (L *ᵥ σ) v = 0) (hτ : ∀ v, ¬I v → τ v = σ v) :
    σ ⬝ᵥ (L *ᵥ σ) ≤ τ ⬝ᵥ (L *ᵥ τ) := by
  rw [harmonic_energy_identity hpsd.1.isSymm hσ hτ]
  have := hpsd.dotProduct_mulVec_nonneg (τ - σ)
  rw [star_trivial] at this
  linarith

/-- **… and it is the only minimizer** when `L_OO` is positive definite.

[formula: cut-harmonic-extension] -/
theorem harmonicExtension_energy_lt (hL : L.IsSymm)
    (hpd : (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val).PosDef) {σ τ : V → ℝ}
    (hσ : ∀ v, I v → (L *ᵥ σ) v = 0) (hτ : ∀ v, ¬I v → τ v = σ v) (hne : τ ≠ σ) :
    σ ⬝ᵥ (L *ᵥ σ) < τ ⬝ᵥ (L *ᵥ τ) := by
  rw [harmonic_energy_identity hL hσ hτ]
  have hd0 : ∀ v, ¬I v → (τ - σ) v = 0 := fun v hv => by simp [hτ v hv]
  rw [dotProduct_mulVec_blockII hd0]
  have hx : (fun j : {v // I v} => (τ - σ) j) ≠ 0 := by
    intro h0
    apply hne
    funext v
    by_cases hv : I v
    · have := congrFun h0 ⟨v, hv⟩
      simp only [Pi.sub_apply, Pi.zero_apply] at this
      linarith
    · exact hτ v hv
  have := hpd.dotProduct_mulVec_pos hx
  rw [star_trivial] at this
  linarith

end Dirichlet

end Cuts

/-! ## On a mesh -/

namespace TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype V] [Fintype F] (M : TriMesh V F)

open Cuts Matrix

variable (I : V → Prop) [DecidablePred I]

/-- **Harmonic extension on a connected valid mesh** with at least one cut vertex: exactly one
extension, and it has the least cotan energy among all extensions of the same data.

[formula: cut-harmonic-extension] -/
theorem existsUnique_harmonicExtension (hM : M.Valid) (hconn : M.Connected) (hB : ∃ v, ¬I v)
    (s : V → ℝ) :
    ∃! σ : V → ℝ, (∀ v, ¬I v → σ v = s v) ∧ ∀ v, I v → (M.laplacian *ᵥ σ) v = 0 :=
  Cuts.existsUnique_harmonicExtension
    (isUnit_det_of_posDef (M.laplacian_interior_posDef hM hconn I hB)) s

/-- The harmonic extension on a mesh minimizes `Σ_f (face energy)`.

[formula: cut-harmonic-extension] -/
theorem harmonicExtension_faceEnergy_le (hM : M.Valid) {σ τ : V → ℝ}
    (hσ : ∀ v, I v → (M.laplacian *ᵥ σ) v = 0) (hτ : ∀ v, ¬I v → τ v = σ v) :
    ∑ f, M.faceEnergy f σ ≤ ∑ f, M.faceEnergy f τ := by
  rw [← M.dotProduct_laplacian_mulVec_eq_sum_faceEnergy hM,
    ← M.dotProduct_laplacian_mulVec_eq_sum_faceEnergy hM]
  exact harmonicExtension_energy_le (M.laplacian_posSemidef hM) hσ hτ

end TriMesh

end DDG
