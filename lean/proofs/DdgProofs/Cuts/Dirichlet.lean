import DdgProofs.Cuts.Weighted
import DdgProofs.Solver
import DdgProofs.Laplace.Cotan

/-!
# Cut flow: the Dirichlet solves in a region

`crates/ddg-variational-cuts/SPEC.md` §3.4 – §3.6, §3.8. Everything here is linear algebra over a
finite index set `V` (the vertices of a submesh) split by a predicate `I` into the unknowns
(`I v`) and the rest (`B`, the cut vertices and the vertices of `∂M`). `L` is any square matrix
unless a hypothesis says otherwise; for the cotan Laplacian of a connected valid mesh with a
nonempty `B` the hypotheses hold (`TriMesh.yamabe_existsUnique`, …).

`IsDirichletSolution L I f u`: `u = 0` off `I` and `(L u)_i = f_i` on `I`, i.e. `u_B = 0`,
`L_II u_I = f_I` (`isDirichletSolution_iff_block`). `boundaryFlux L I f u` is `(L u)_i − f_i` at
`i ∉ I` and `0` on `I`.

* **`cut-yamabe-system`** (`f = −g`): `existsUnique_isDirichletSolution`, `yamabe_existsUnique`:
  a unique solution when `L_II` is invertible (positive definite).
* **`cut-normal-derivative`** (`v = boundaryFlux L I (−g) u`, so `v_i = (L_BI u_I)_i + g_i`,
  `Ddg.normalFlux`): `normalFlux_spec`: `v` is the unique vector supported off `I` with
  `L u = −g + v`; `sum_normalFlux`: if the columns of `L` sum to zero, `Σ_B v = Σ g` (flux
  balance, the discrete divergence theorem).
* **`cut-dirichlet-energy`**: `dirichletEnergy_eq`: `uᵀ L u = −Σ_I u_i g_i`.
* **`cut-adjoint-scaled`**, **`cut-gradient-scaled`**: with `L_a = c·L` and `a ≡ c` the adjoint
  is `p = −2c·u` and its flux is `q = −2c·v` (`scaledAdjoint_const`), so the scaled density is
  `c` times the Dirichlet one (`gradScaled_const`; `c = 1`: `gradScaled_one`). For any symmetric
  `L_a` the adjoint is exact for the discrete objective: `scaledEnergy_perturb`.
* **`cut-adjoint-pointwise`**: `hasDerivAt_pointwiseEnergy`: for a perturbation `δg` of the
  source, the derivative of `Σ_i A_i b(u_i)` is `Σ_I p_i δg_i`.
* **`cut-gradient-dirichlet`**, **`cut-gradient-pointwise`**: the oracle's densities in terms of
  `∂u/∂n ≈ v/ℓ` (`gradDirichlet_eq`, `gradPointwise_eq`). They are discretizations of the shape
  derivative (cited: `cut-shape-derivative`, `cut-cea-method`), checked on closed-form solutions
  in `DdgProofs/Cuts/Analytic.lean`.
* **`cut-harmonic-extension`**: `existsUnique_dirichlet_data`: unique `σ` with `σ = s` off `O`
  and `(L σ)_i = 0` on `O` (`L_OO σ_O = −L_OC s_C`, `harmonic_iff_block`);
  `harmonicExtension_energy_le`: it minimizes `τᵀ L τ` among all `τ` with the same data, with
  equality only for `τ = σ` when `L_OO` is positive definite.
-/

noncomputable section

namespace DDG

namespace Cuts

open Matrix
open scoped DDG.Oracle

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

/-! ## The Yamabe system, its normal derivative and its energy -/

/-- **The Yamabe system has exactly one solution** when `L_II` is invertible:
`u_B = 0`, `L_II u_I = −g_I` (SPEC §3.4).

[formula: cut-yamabe-system] -/
theorem yamabe_existsUnique [DecidableEq V]
    (hdet : IsUnit (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val).det) (g : V → ℝ) :
    ∃! u : V → ℝ, IsDirichletSolution L I (-g) u :=
  existsUnique_isDirichletSolution hdet (-g)

/-- The Yamabe system in block form: `u_B = 0`, `L_II u_I = −g_I`.

[formula: cut-yamabe-system] -/
theorem yamabe_iff_block {g u : V → ℝ} :
    IsDirichletSolution L I (-g) u ↔ (∀ v, ¬I v → u v = 0) ∧
      (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val) *ᵥ (fun j : {v // I v} => u j) =
        fun j : {v // I v} => -g j :=
  isDirichletSolution_iff_block

/-- **The integrated normal derivative** at a vertex outside `I` is the oracle's
`normalFlux ((L_BI u_I)_i) g_i = (L_BI u_I)_i + g_i` (SPEC §3.5).

[formula: cut-normal-derivative] -/
theorem normalFlux_apply {g u : V → ℝ} (hu : ∀ v, ¬I v → u v = 0) (i : {v // ¬I v}) :
    boundaryFlux L I (-g) u i =
      Ddg.normalFlux (((L.submatrix (Subtype.val : {v // ¬I v} → V)
        (Subtype.val : {v // I v} → V)) *ᵥ fun j : {v // I v} => u j) i) (g i) := by
  unfold Ddg.normalFlux
  rw [← mulVec_eq_blockBI hu i]
  simp [boundaryFlux, i.2]

/-- **The partitioned Yamabe system** (paper Eq. 26, with the sign of D6): with
`v = boundaryFlux L I (−g) u`: `v` vanishes on `I`, `L u = −g + [0; v]` at every vertex, and
`v` is the only vector with this property.

[formula: cut-normal-derivative] -/
theorem normalFlux_spec {g u : V → ℝ} (h : IsDirichletSolution L I (-g) u) :
    (∀ v, I v → boundaryFlux L I (-g) u v = 0) ∧ L *ᵥ u = -g + boundaryFlux L I (-g) u ∧
      ∀ w : V → ℝ, L *ᵥ u = -g + w → w = boundaryFlux L I (-g) u :=
  ⟨fun _ hv => boundaryFlux_of_mem _ _ hv, h.mulVec_eq, fun _ hw => h.flux_unique hw⟩

/-- **Flux balance**: if the columns of `L` sum to zero then `Σ_B v = Σ_V g`. With
`g_i = Σ_t k_t A_t/3` this says that the integrated normal derivative over the boundary of the
region is its total curvature (the divergence theorem for `Δu = −K`).

[formula: cut-normal-derivative] -/
theorem sum_normalFlux {g u : V → ℝ} (h : IsDirichletSolution L I (-g) u)
    (hL : (fun _ => (1 : ℝ)) ᵥ* L = 0) : ∑ v, boundaryFlux L I (-g) u v = ∑ v, g v := by
  rw [h.sum_boundaryFlux hL]
  simp

/-- **Flux balance for a symmetric matrix with zero row sums** (the form of SPEC §3.5; for a
matrix that is not symmetric the hypothesis that matters is zero *column* sums, `sum_normalFlux`).

[formula: cut-normal-derivative] -/
theorem sum_normalFlux_of_isSymm {g u : V → ℝ} (h : IsDirichletSolution L I (-g) u)
    (hL : L.IsSymm) (hrow : L *ᵥ (fun _ => (1 : ℝ)) = 0) :
    ∑ v, boundaryFlux L I (-g) u v = ∑ v, g v := by
  refine sum_normalFlux h ?_
  rw [← Matrix.mulVec_transpose, hL.eq, hrow]

/-- **Dirichlet energy of the Yamabe solution**: `uᵀ L u = −Σ u_i g_i = −u_I · g_I`.

[formula: cut-dirichlet-energy] -/
theorem dirichletEnergy_eq {g u : V → ℝ} (h : IsDirichletSolution L I (-g) u) :
    u ⬝ᵥ (L *ᵥ u) = -∑ v, u v * g v := by
  rw [h.dotProduct_mulVec h.zero_off]
  simp

/-- The same with the sum over the unknowns only.

[formula: cut-dirichlet-energy] -/
theorem dirichletEnergy_eq_interior {g u : V → ℝ} (h : IsDirichletSolution L I (-g) u) :
    u ⬝ᵥ (L *ᵥ u) = -∑ i : {v // I v}, u i * g i := by
  rw [dirichletEnergy_eq h]
  congr 1
  rw [← Finset.sum_subtype (Finset.univ.filter I) (by simp) (fun i => u i * g i),
    Finset.sum_filter]
  refine Finset.sum_congr rfl fun v _ => ?_
  by_cases hv : I v
  · simp [hv]
  · simp [hv, h.zero_off v hv]

/-- The Dirichlet energy is nonnegative for a positive semidefinite `L`. -/
theorem dirichletEnergy_nonneg (hpsd : L.PosSemidef) {g u : V → ℝ}
    (h : IsDirichletSolution L I (-g) u) : 0 ≤ -∑ v, u v * g v := by
  rw [← dirichletEnergy_eq h]
  have := hpsd.dotProduct_mulVec_nonneg u
  rwa [star_trivial] at this

/-! ## Gradient densities -/

/-- The oracle's `normalDerivative` at `ℝ`: `v/ℓ`; multiplied by the dual length it gives back
the integrated value. -/
theorem normalDerivative_mul (v : ℝ) {ℓ : ℝ} (hℓ : ℓ ≠ 0) : Ddg.normalDerivative v ℓ * ℓ = v := by
  unfold Ddg.normalDerivative
  field_simp

/-- The oracle's `dualLength` of two boundary edges: half their sum. -/
theorem dualLength_pair (l1 l2 : ℝ) : Ddg.dualLength [l1, l2] = (l1 + l2) / 2 := by
  simp [Ddg.dualLength, Ddg.sumList]

/-- The oracle's `dualLength` is positive for two positive edges. -/
theorem dualLength_pair_pos {l1 l2 : ℝ} (h1 : 0 < l1) (h2 : 0 < l2) :
    0 < Ddg.dualLength [l1, l2] := by
  rw [dualLength_pair]
  positivity

/-- **Dirichlet density**: the oracle's `gradDirichlet v ℓ` is `(∂u/∂n)²` with
`∂u/∂n = v/ℓ` (SPEC §3.6, paper Eq. 13).

[formula: cut-gradient-dirichlet] -/
theorem gradDirichlet_eq (v ℓ : ℝ) : Ddg.gradDirichlet v ℓ = Ddg.normalDerivative v ℓ ^ 2 := by
  unfold Ddg.gradDirichlet Ddg.normalDerivative
  ring

/-- The Dirichlet density is nonnegative.

[formula: cut-gradient-dirichlet] -/
theorem gradDirichlet_nonneg (v ℓ : ℝ) : 0 ≤ Ddg.gradDirichlet v ℓ := by
  rw [gradDirichlet_eq]
  positivity

/-- **Scaled density**: `G = −a (∂u/∂n)² − (∂u/∂n)(∂p/∂n)` (paper Eq. 18, Table 1).

[formula: cut-gradient-scaled] -/
theorem gradScaled_eq (a v q ℓ : ℝ) :
    Ddg.gradScaled a v q ℓ =
      -(a * Ddg.normalDerivative v ℓ ^ 2) - Ddg.normalDerivative v ℓ * Ddg.normalDerivative q ℓ := by
  unfold Ddg.gradScaled Ddg.normalDerivative
  ring

/-- **With a constant weight `c` and the adjoint flux `q = −2c·v` the scaled density is `c`
times the Dirichlet density.**

[formula: cut-gradient-scaled] -/
theorem gradScaled_const (c v ℓ : ℝ) :
    Ddg.gradScaled c v (-2 * c * v) ℓ = c * Ddg.gradDirichlet v ℓ := by
  unfold Ddg.gradScaled Ddg.gradDirichlet
  ring

/-- **With `a ≡ 1`, `q = −2v`, the scaled density is the Dirichlet density** (SPEC §1, §3.6).

[formula: cut-gradient-scaled] -/
theorem gradScaled_one (v ℓ : ℝ) : Ddg.gradScaled 1 v (-2 * v) ℓ = Ddg.gradDirichlet v ℓ := by
  have := gradScaled_const 1 v ℓ
  simpa using this

/-- **Pointwise density**: `G = b(0) − (∂u/∂n)(∂p/∂n)` (paper Table 1; `u = 0` on the cut).

[formula: cut-gradient-pointwise] -/
theorem gradPointwise_eq (b0 v q ℓ : ℝ) :
    Ddg.gradPointwise b0 v q ℓ = b0 - Ddg.normalDerivative v ℓ * Ddg.normalDerivative q ℓ := by
  unfold Ddg.gradPointwise Ddg.normalDerivative
  ring

/-! ### The densities in terms of recovered normal derivatives

SPEC §3.6 evaluates the densities at `d ≈ ∂u/∂n`, `e ≈ ∂p/∂n` recovered from the integrated
fluxes (`cut-flux-recovery`, `DdgProofs/Cuts/Recovery.lean`); with the nodal recovery `d = v/ℓ`,
`e = q/ℓ` they are the expressions above. -/

/-- The oracle's `gradDirichletD` at `ℝ`: `G = d²`. -/
theorem gradDirichletD_eq (d : ℝ) : Ddg.gradDirichletD d = d ^ 2 := by
  unfold Ddg.gradDirichletD
  ring

/-- The oracle's `gradScaledD` at `ℝ`: `G = −a d² − d e`. -/
theorem gradScaledD_eq (a d e : ℝ) : Ddg.gradScaledD a d e = -(a * d ^ 2) - d * e := by
  unfold Ddg.gradScaledD
  ring

/-- The oracle's `gradPointwiseD` at `ℝ`: `G = b(0) − d e`. -/
theorem gradPointwiseD_eq (b0 d e : ℝ) : Ddg.gradPointwiseD b0 d e = b0 - d * e := rfl

/-- The Dirichlet density is nonnegative.

[formula: cut-gradient-dirichlet] -/
theorem gradDirichletD_nonneg (d : ℝ) : 0 ≤ Ddg.gradDirichletD d := by
  rw [gradDirichletD_eq]
  positivity

/-- **With the nodal recovery `d = v/ℓ` the Dirichlet density is the paper's `(v/ℓ)²`.**

[formula: cut-gradient-dirichlet] -/
theorem gradDirichletD_nodal (v ℓ : ℝ) :
    Ddg.gradDirichletD (Ddg.normalDerivative v ℓ) = Ddg.gradDirichlet v ℓ := by
  rw [gradDirichletD_eq, gradDirichlet_eq]

/-- **With the nodal recovery the scaled density is the paper's.**

[formula: cut-gradient-scaled] -/
theorem gradScaledD_nodal (a v q ℓ : ℝ) :
    Ddg.gradScaledD a (Ddg.normalDerivative v ℓ) (Ddg.normalDerivative q ℓ) =
      Ddg.gradScaled a v q ℓ := by
  rw [gradScaledD_eq, gradScaled_eq]

/-- **With the nodal recovery the pointwise density is the paper's.**

[formula: cut-gradient-pointwise] -/
theorem gradPointwiseD_nodal (b0 v q ℓ : ℝ) :
    Ddg.gradPointwiseD b0 (Ddg.normalDerivative v ℓ) (Ddg.normalDerivative q ℓ) =
      Ddg.gradPointwise b0 v q ℓ := by
  rw [gradPointwiseD_eq, gradPointwise_eq]

/-- **With a constant weight `c` and `e = −2c·d` the scaled density is `c` times the Dirichlet
density** (`e = −2c·d` because `q = −2c·v`, `scaledAdjoint_const`, and the recovery is linear).

[formula: cut-gradient-scaled] -/
theorem gradScaledD_const (c d : ℝ) :
    Ddg.gradScaledD c d (-2 * c * d) = c * Ddg.gradDirichletD d := by
  rw [gradScaledD_eq, gradDirichletD_eq]
  ring

/-- **With `a ≡ 1`, `e = −2d`, the scaled density is the Dirichlet density.**

[formula: cut-gradient-scaled] -/
theorem gradScaledD_one (d : ℝ) : Ddg.gradScaledD 1 d (-2 * d) = Ddg.gradDirichletD d := by
  have := gradScaledD_const 1 d
  simpa using this

/-! ## The adjoint of the scaled Dirichlet energy -/

/-- Right-hand side of the adjoint of the scaled energy (SPEC §3.6), from the oracle's entries:
`r_i = −2 (L_a u)_i` on `I`, `r_i = −2 ((L_a u)_i − a_i v_i)` outside. -/
def scaledAdjointRhs (La : Matrix V V ℝ) (I : V → Prop) [DecidablePred I] (a u v : V → ℝ) :
    V → ℝ :=
  fun i => if I i then Ddg.adjointRhsScaledInterior ((La *ᵥ u) i)
    else Ddg.adjointRhsScaledBoundary ((La *ᵥ u) i) (a i) (v i)

theorem scaledAdjointRhs_of_mem (La : Matrix V V ℝ) (a u v : V → ℝ) {i : V} (hi : I i) :
    scaledAdjointRhs La I a u v i = -2 * (La *ᵥ u) i := by
  simp [scaledAdjointRhs, hi, Ddg.adjointRhsScaledInterior]

theorem scaledAdjointRhs_of_not_mem (La : Matrix V V ℝ) (a u v : V → ℝ) {i : V} (hi : ¬I i) :
    scaledAdjointRhs La I a u v i = -2 * ((La *ᵥ u) i - a i * v i) := by
  simp [scaledAdjointRhs, hi, Ddg.adjointRhsScaledBoundary]

/-- **The adjoint flux** at a vertex outside `I` is the oracle's
`adjointFlux ((L_BI p_I)_i) r_i = (L_BI p_I)_i − r_i` (SPEC §3.6).

[formula: cut-adjoint-scaled] -/
theorem adjointFlux_apply {r p : V → ℝ} (hp : ∀ v, ¬I v → p v = 0) (i : {v // ¬I v}) :
    boundaryFlux L I r p i =
      Ddg.adjointFlux (((L.submatrix (Subtype.val : {v // ¬I v} → V)
        (Subtype.val : {v // I v} → V)) *ᵥ fun j : {v // I v} => p j) i) (r i) := by
  unfold Ddg.adjointFlux
  rw [← mulVec_eq_blockBI hp i]
  simp [boundaryFlux, i.2]

/-- **Constant weight: the scaled adjoint is `−2c·u` and its flux is `−2c·v`.** Let `u` vanish
off `I`, let `v` be any vector, `L_a = c·L` and `a ≡ c`. Then `p = −2c·u` solves the adjoint
system `L_II p_I = r_I`, `p_B = 0` with the right-hand side `r` of SPEC §3.6, and its flux
`q_i = (L_BI p_I)_i − r_i` is `−2c·v_i` at every `i ∉ I`. With `c = 1`: `p = −2u`, `q = −2v`.

[formula: cut-adjoint-scaled] -/
theorem scaledAdjoint_const (c : ℝ) {u : V → ℝ} (hu : ∀ v, ¬I v → u v = 0) (v : V → ℝ) :
    IsDirichletSolution L I (scaledAdjointRhs (c • L) I (fun _ => c) u v) ((-2 * c) • u) ∧
      ∀ i, ¬I i →
        boundaryFlux L I (scaledAdjointRhs (c • L) I (fun _ => c) u v) ((-2 * c) • u) i =
          -2 * c * v i := by
  have hLp : ∀ i, (L *ᵥ ((-2 * c) • u)) i = -2 * c * (L *ᵥ u) i := by
    intro i
    rw [Matrix.mulVec_smul]
    simp
  have hLa : ∀ i, ((c • L) *ᵥ u) i = c * (L *ᵥ u) i := by
    intro i
    rw [Matrix.smul_mulVec]
    simp
  refine ⟨⟨fun i hi => by simp [hu i hi], fun i hi => ?_⟩, fun i hi => ?_⟩
  · rw [hLp, scaledAdjointRhs_of_mem _ _ _ _ hi, hLa]
    ring
  · simp only [boundaryFlux, hi, ite_false]
    rw [hLp, scaledAdjointRhs_of_not_mem _ _ _ _ hi, hLa]
    ring

/-- **… and it is the only solution** when `L_II` is invertible.

[formula: cut-adjoint-scaled] -/
theorem scaledAdjoint_const_unique [DecidableEq V]
    (hdet : IsUnit (L.submatrix (Subtype.val : {v // I v} → V) Subtype.val).det) (c : ℝ)
    {u : V → ℝ} (hu : ∀ v, ¬I v → u v = 0) (v : V → ℝ) {p : V → ℝ}
    (hp : IsDirichletSolution L I (scaledAdjointRhs (c • L) I (fun _ => c) u v) p) :
    p = (-2 * c) • u :=
  (existsUnique_isDirichletSolution hdet _).unique hp (scaledAdjoint_const c hu v).1

/-- **The scaled adjoint is exact for the discrete objective.** Let `L` and `L_a` be symmetric,
`δu` the Yamabe solution for the source perturbation `δg` (if `u` is the solution for `g` then
`u + t δu` is the solution for `g + t δg`, `IsDirichletSolution.add_smul`), and `p` the adjoint:
`p_B = 0`, `(L p)_i = −2 (L_a u)_i` on `I`. Then for every `t`

`(u + t δu)ᵀ L_a (u + t δu) = uᵀ L_a u + t · Σ_I p_i δg_i + t² · δuᵀ L_a δu`,

so the first variation of `J = uᵀ L_a u` with respect to the source is `p_I · δg_I`.

[formula: cut-adjoint-scaled] -/
theorem scaledEnergy_perturb (hL : L.IsSymm) {La : Matrix V V ℝ} (hLa : La.IsSymm)
    {dg u du p r : V → ℝ}
    (hdu : IsDirichletSolution L I (-dg) du) (hp : IsDirichletSolution L I r p)
    (hr : ∀ i, I i → r i = Ddg.adjointRhsScaledInterior ((La *ᵥ u) i)) (t : ℝ) :
    (u + t • du) ⬝ᵥ (La *ᵥ (u + t • du)) =
      u ⬝ᵥ (La *ᵥ u) + t * ∑ v, p v * dg v + t ^ 2 * (du ⬝ᵥ (La *ᵥ du)) := by
  have hrec := hdu.reciprocity hL hp
  have h1 : du ⬝ᵥ (La *ᵥ u) = -(1 / 2) * ∑ v, du v * r v := by
    simp only [dotProduct]
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun v _ => ?_
    by_cases hv : I v
    · rw [hr v hv]
      unfold Ddg.adjointRhsScaledInterior
      simp only [Oracle.two_eq]
      ring
    · rw [hdu.zero_off v hv]
      ring
  have h2 : ∑ v, p v * (-dg) v = -∑ v, p v * dg v := by simp
  have h3 : u ⬝ᵥ (La *ᵥ du) = du ⬝ᵥ (La *ᵥ u) := dotProduct_mulVec_symm hLa u du
  rw [Matrix.mulVec_add, Matrix.mulVec_smul, add_dotProduct, dotProduct_add, dotProduct_add,
    smul_dotProduct, smul_dotProduct, dotProduct_smul, dotProduct_smul, h3, h1, ← hrec, h2]
  simp only [smul_eq_mul]
  ring

/-! ## The adjoint of a pointwise energy -/

/-- **The pointwise adjoint is exact for the discrete objective** (SPEC §3.6). Let `L` be
symmetric, `δu` the Yamabe solution for the source perturbation `δg` (`δu_B = 0`,
`L_II δu_I = −δg_I`; `u + t δu` is then the solution for `g + t δg`), and `p` the adjoint:
`p_B = 0`, `(L p)_i = r_i = −A_i b′(u_i)` on `I` (the oracle's `adjointRhsPointwise`). Then

`d/dt Σ_i A_i b(u_i + t δu_i) |_{t=0} = Σ_I p_i δg_i`.

[formula: cut-adjoint-pointwise] -/
theorem hasDerivAt_pointwiseEnergy (hL : L.IsSymm) {A u du dg p r : V → ℝ} {b b' : ℝ → ℝ}
    (hb : ∀ v, HasDerivAt b (b' (u v)) (u v))
    (hdu : IsDirichletSolution L I (-dg) du) (hp : IsDirichletSolution L I r p)
    (hr : ∀ i, I i → r i = Ddg.adjointRhsPointwise (A i) (b' (u i))) :
    HasDerivAt (fun t : ℝ => ∑ v, A v * b (u v + t * du v)) (∑ v, p v * dg v) 0 := by
  have hd : ∀ v, HasDerivAt (fun t : ℝ => A v * b (u v + t * du v))
      (A v * (b' (u v) * du v)) 0 := by
    intro v
    have h1 : HasDerivAt (fun t : ℝ => u v + t * du v) (du v) 0 := by
      simpa using ((hasDerivAt_id (0 : ℝ)).mul_const (du v)).const_add (u v)
    have h2 : HasDerivAt b (b' (u v)) (u v + 0 * du v) := by
      simpa using hb v
    exact (h2.comp 0 h1).const_mul (A v)
  have hsum := HasDerivAt.fun_sum (u := Finset.univ) (fun v _ => hd v)
  have hrec := hdu.reciprocity hL hp
  have hval : ∑ v, A v * (b' (u v) * du v) = ∑ v, p v * dg v := by
    have h1 : ∑ v, A v * (b' (u v) * du v) = -∑ v, du v * r v := by
      rw [← Finset.sum_neg_distrib]
      refine Finset.sum_congr rfl fun v _ => ?_
      by_cases hv : I v
      · rw [hr v hv]
        unfold Ddg.adjointRhsPointwise
        ring
      · rw [hdu.zero_off v hv]
        ring
    rw [h1, ← hrec]
    simp
  rw [← hval]
  exact hsum

/-- **Hencky energy** `b(u) = u²` (paper Eq. 19): with `(L p)_i = −2 A_i u_i` on `I`,
`d/dt Σ_i A_i (u_i + t δu_i)² |_{t=0} = Σ_I p_i δg_i`.

[formula: cut-adjoint-pointwise] -/
theorem hasDerivAt_henckyEnergy (hL : L.IsSymm) {A u du dg p r : V → ℝ}
    (hdu : IsDirichletSolution L I (-dg) du) (hp : IsDirichletSolution L I r p)
    (hr : ∀ i, I i → r i = Ddg.adjointRhsPointwise (A i) (2 * u i)) :
    HasDerivAt (fun t : ℝ => ∑ v, A v * (u v + t * du v) ^ 2) (∑ v, p v * dg v) 0 :=
  hasDerivAt_pointwiseEnergy (b := fun x => x ^ 2) (b' := fun x => 2 * x) hL
    (fun v => by simpa using hasDerivAt_pow 2 (u v)) hdu hp hr

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

/-- A per-face quantity distributed to the vertices in thirds: `Σ_{f ∋ v} c_f / 3`. The source
of the Yamabe system is `g = lumped (k_f A_f)` (SPEC §3.4, the paper's `W_F K`); the lumped mass
is `lumped (A_f)`. -/
def lumped (c : F → ℝ) (v : V) : ℝ :=
  (∑ f, if v = M.a f ∨ v = M.b f ∨ v = M.c f then c f else 0) / 3

/-- Distributing to the vertices keeps the total: `Σ_v lumped c v = Σ_f c_f`. -/
theorem sum_lumped (hM : M.Valid) (c : F → ℝ) : ∑ v, M.lumped c v = ∑ f, c f := by
  unfold lumped
  rw [← Finset.sum_div, Finset.sum_comm, Finset.sum_div]
  refine Finset.sum_congr rfl fun f _ => ?_
  have hab := hM.ab f
  have hbc := hM.bc f
  have hca := hM.ca f
  have h3 : ∀ v : V, (if v = M.a f ∨ v = M.b f ∨ v = M.c f then c f else 0) =
      (if v = M.a f then c f else 0) + (if v = M.b f then c f else 0) +
        (if v = M.c f then c f else 0) := by
    intro v
    by_cases h1 : v = M.a f <;> by_cases h2 : v = M.b f <;> by_cases h3 : v = M.c f <;> simp_all
  simp only [h3, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  ring

variable (I : V → Prop) [DecidablePred I]

/-- **The Yamabe system of a connected valid mesh with at least one vertex outside the unknowns
has exactly one solution** (the block `L_II` of the cotan Laplacian is positive definite).

[formula: cut-yamabe-system] -/
theorem yamabe_existsUnique (hM : M.Valid) (hconn : M.Connected) (hB : ∃ v, ¬I v) (g : V → ℝ) :
    ∃! u : V → ℝ, IsDirichletSolution M.laplacian I (-g) u :=
  Cuts.yamabe_existsUnique (isUnit_det_of_posDef (M.laplacian_interior_posDef hM hconn I hB)) g

/-- **Flux balance on a mesh**: the integrated normal derivatives sum to the total source. With
`g = lumped (k_f A_f)`: `Σ_B v = Σ_f k_f A_f`, the total curvature of the region.

[formula: cut-normal-derivative] -/
theorem sum_normalFlux_lumped (hM : M.Valid) (c : F → ℝ) {u : V → ℝ}
    (h : IsDirichletSolution M.laplacian I (-M.lumped c) u) :
    ∑ v, boundaryFlux M.laplacian I (-M.lumped c) u v = ∑ f, c f := by
  rw [sum_normalFlux h M.laplacian_vecMul_one, M.sum_lumped hM]

/-- **The Dirichlet energy on a mesh** is the sum of the cotan energies of the faces:
`−Σ_I u_i g_i = uᵀ L u = Σ_f (face energy) ≥ 0`.

[formula: cut-dirichlet-energy] -/
theorem dirichletEnergy_eq_sum_faceEnergy (hM : M.Valid) {g u : V → ℝ}
    (h : IsDirichletSolution M.laplacian I (-g) u) :
    -∑ v, u v * g v = ∑ f, M.faceEnergy f u := by
  rw [← dirichletEnergy_eq h, M.dotProduct_laplacian_mulVec_eq_sum_faceEnergy hM]

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
