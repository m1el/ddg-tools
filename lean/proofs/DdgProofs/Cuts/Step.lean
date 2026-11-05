import DdgProofs.Laplace.Dirichlet
import DdgProofs.Oracle.Scalar
import DdgProofs.Cuts.Weighted

/-!
# Cut flow: speed, explicit step, length step, symmetry, step size, feature size

`crates/ddg-variational-cuts/SPEC.md` §1 (signs), §3.6 – §3.11, §4, D2, D4, D8.

* **`cut-speed`** (`Ddg.cutSpeed`, `Ddg.regionDensity`): if the first variation of the energy for
  a normal speed `σ` of the cut (along the outward normal of region 0) is `dE = Σ_c g_c σ_c ℓ_c`,
  the level-set relation is `φ_t = −σ` and the flow sets `φ_t = g`, then
  `dE = −Σ_c g_c² ℓ_c ≤ 0` (`firstVariation_levelSet`, `firstVariation_cutSpeed`); it is the
  steepest descent direction in the `ℓ`-weighted norm (`firstVariation_ge`); the opposite sign
  (the literal reading of the paper, D2) gives `+Σ g² ℓ ≥ 0` (`firstVariation_ascent`).
* **`cut-explicit-step`** (`Ddg.explicitStep`): `φ̃ − φ = τ σ`, so by the above one explicit step
  changes the energy to first order by `−τ Σ_c s_c² ℓ_c ≤ 0` (`firstVariation_explicitStep`).
* **`cut-length-step`** (`Ddg.lengthStepEntry`): `W_V + τ L_w` is positive definite for positive
  masses, positive semidefinite `L_w` and `τ ≥ 0` (`lengthStep_posDef`); the step has exactly one
  solution (`lengthStep_existsUnique`); it does not increase `φᵀ L_w φ` nor `φᵀ W_V φ`
  (`lengthStep_energy_le`, `lengthStep_mass_le`); it keeps `Σ A_i φ_i` if the columns of `L_w`
  sum to zero (`lengthStep_sum`); it is the minimizer of
  `½ (ψ − φ̃)ᵀ W_V (ψ − φ̃) + (τ/2) ψᵀ L_w ψ` (`lengthStep_minimizes`).
* **`cut-penalty-step-paper`** (`Ddg.penaltyStepEntry`): `(W D⁻¹ + τ L) φ′ = W D⁻¹ φ̃` is
  equivalent to `(W + τ D L) φ′ = W φ̃` for a positive diagonal `D` (`penaltyStep_iff`), has
  exactly one solution (`penaltyStep_existsUnique`), and for a constant `D = d` it is the default
  step with `L_w = d L` (`penaltyStep_const_iff`).
* **`cut-symmetry-projection`** (`Ddg.classMean`): the class mean is idempotent, fixes exactly the
  class-constant vectors, is self-adjoint, and is the orthogonal projection onto the
  class-constant vectors: `Σ(φ − ψ)² = Σ(φ − Pφ)² + Σ(Pφ − ψ)²` for every class-constant `ψ`
  (`classMean_idem`, `classMean_eq_self_iff`, `classMean_selfAdjoint`, `classMean_pythagoras`).
* **`cut-step-cap`** (`Ddg.stepCap`): `τ ≤ τ_max`, `τ · max|s| ≤ max_motion · h`, and `τ` is the
  largest number with both properties (`stepCap_le_max`, `stepCap_mul_le`, `le_stepCap`).
* **`cut-feature-size`** (`Ddg.featureSizeRaw`, `Ddg.featureWeight`): `0 < r⁰ ≤ r_ε`, with
  `r⁰ = r_ε` exactly at `K = 0` and `r⁰ < 1/√|K|` otherwise (`featureSizeRaw_pos`,
  `featureSizeRaw_le`, `featureSizeRaw_lt_inv_sqrt`); the smoothing `(W_V + c L) r = W_V r⁰`
  keeps the bounds of `r⁰` **if the off-diagonal entries of `L` are nonpositive**
  (`smoothing_ge`, `smoothing_le`); with the cotan Laplacian itself that is not guaranteed on a
  mesh with negative cotangent weights, with the Laplacian `L⁺` of the edge weights
  `max(w_ij, 0)` that SPEC §3.6 uses (D17) it holds on every mesh
  (`smoothedFeatureSize_bounds_positivePart`).
* **`cut-sliver-smoothing`** (`Ddg.aspectRatio`, `Ddg.smoothedSpeed`): the smoothed speed is a
  convex combination of the two neighbours' speeds with the weights `ℓ₂/(ℓ₁+ℓ₂)`, `ℓ₁/(ℓ₁+ℓ₂)`,
  hence within their range, and it is exact for a speed that is linear in arclength
  (`smoothedSpeed_eq_convex`, `smoothedSpeed_mem`, `smoothedSpeed_linear`); the aspect ratio is
  `1` for an equilateral triangle and `≥ 1` for every triangle (`aspectRatio_equilateral`,
  `one_le_aspectRatio`).
-/

noncomputable section

namespace DDG

namespace Cuts

open Matrix
open scoped DDG.Oracle
open Real

/-! ## Speed on the cut and the sign of the flow -/

/-- The oracle's `cutSpeed` at `ℝ`: `s = G_0 − G_1`. -/
theorem cutSpeed_eq (g0 g1 : ℝ) : Ddg.cutSpeed g0 g1 = g0 - g1 := rfl

/-- The oracle's `regionDensity` at `ℝ`: `G_r = w_r G + β_r`. -/
theorem regionDensity_eq (w G β : ℝ) : Ddg.regionDensity w G β = w * G + β := rfl

/-- With the default weights `w_r = 1`, `β_r = 0` the region total is the density itself. -/
theorem regionDensity_default (G : ℝ) : Ddg.regionDensity 1 G 0 = G := by
  rw [regionDensity_eq]
  ring

/-- A region with `w_r = 0` contributes `β_r` whatever its density (it needs no solve). -/
theorem regionDensity_zero_weight (G β : ℝ) : Ddg.regionDensity 0 G β = β := by
  rw [regionDensity_eq]
  ring

/-- Exchanging the two regions changes the sign of the speed (and `φ ↦ −φ`). -/
theorem cutSpeed_swap (g0 g1 : ℝ) : Ddg.cutSpeed g1 g0 = -Ddg.cutSpeed g0 g1 := by
  simp only [cutSpeed_eq]
  ring

section FirstVariation

variable {C : Type*} [Fintype C]

/-- First variation of the energy when the cut moves with normal speed `σ` (along the outward
normal of region 0), for a gradient density `g` and dual lengths `ℓ` at the cut vertices:
`dE(σ) = Σ_c g_c σ_c ℓ_c` (the discrete `∫_γ g σ ds` of SPEC §1). -/
def firstVariation (g ℓ σ : C → ℝ) : ℝ := ∑ c, g c * σ c * ℓ c

/-- **The flow is a descent.** If the level-set function changes by `φ_t = −σ` when the cut
moves with normal speed `σ`, and the flow sets `φ_t = g`, then
`dE = −Σ_c g_c² ℓ_c`.

[formula: cut-speed] -/
theorem firstVariation_levelSet {g ℓ σ φt : C → ℝ} (hσ : ∀ c, φt c = -σ c)
    (hg : ∀ c, φt c = g c) : firstVariation g ℓ σ = -∑ c, g c ^ 2 * ℓ c := by
  unfold firstVariation
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun c _ => ?_
  have : σ c = -g c := by rw [← hg c, hσ c]; ring
  rw [this]
  ring

/-- … which is `≤ 0` for nonnegative lengths.

[formula: cut-speed] -/
theorem firstVariation_levelSet_nonpos {g ℓ σ φt : C → ℝ} (hℓ : ∀ c, 0 ≤ ℓ c)
    (hσ : ∀ c, φt c = -σ c) (hg : ∀ c, φt c = g c) : firstVariation g ℓ σ ≤ 0 := by
  rw [firstVariation_levelSet hσ hg, neg_nonpos]
  exact Finset.sum_nonneg fun c _ => mul_nonneg (sq_nonneg _) (hℓ c)

/-- With positive lengths the first variation along the flow vanishes only at a stationary cut
(`g = 0`). -/
theorem firstVariation_levelSet_eq_zero_iff {g ℓ σ φt : C → ℝ} (hℓ : ∀ c, 0 < ℓ c)
    (hσ : ∀ c, φt c = -σ c) (hg : ∀ c, φt c = g c) :
    firstVariation g ℓ σ = 0 ↔ ∀ c, g c = 0 := by
  rw [firstVariation_levelSet hσ hg, neg_eq_zero,
    Finset.sum_eq_zero_iff_of_nonneg fun c _ => mul_nonneg (sq_nonneg _) (hℓ c).le]
  constructor
  · intro h c
    have := h c (Finset.mem_univ c)
    rcases mul_eq_zero.mp this with h1 | h1
    · exact pow_eq_zero_iff (two_ne_zero) |>.mp h1
    · exact absurd h1 (hℓ c).ne'
  · intro h c _
    rw [h c]
    ring

/-- **The opposite sign is an ascent** (SPEC D2: the paper's Eq. 14–15 read literally). With
`φ_t = −σ` and `φ_t = −g`: `dE = +Σ_c g_c² ℓ_c ≥ 0`. -/
theorem firstVariation_ascent {g ℓ σ φt : C → ℝ} (hℓ : ∀ c, 0 ≤ ℓ c) (hσ : ∀ c, φt c = -σ c)
    (hg : ∀ c, φt c = -g c) :
    firstVariation g ℓ σ = ∑ c, g c ^ 2 * ℓ c ∧ 0 ≤ firstVariation g ℓ σ := by
  have h1 : firstVariation g ℓ σ = ∑ c, g c ^ 2 * ℓ c := by
    unfold firstVariation
    refine Finset.sum_congr rfl fun c _ => ?_
    have : σ c = g c := by
      have := hσ c
      rw [hg c] at this
      linarith
    rw [this]
    ring
  exact ⟨h1, h1 ▸ Finset.sum_nonneg fun c _ => mul_nonneg (sq_nonneg _) (hℓ c)⟩

/-- **Steepest descent.** Among all normal speeds `σ` that are not larger than `g` in the
`ℓ`-weighted norm, none decreases the energy faster than `σ = −g`:
`dE(σ) ≥ −Σ_c g_c² ℓ_c`. -/
theorem firstVariation_ge {g ℓ σ : C → ℝ} (hℓ : ∀ c, 0 ≤ ℓ c)
    (hnorm : ∑ c, σ c ^ 2 * ℓ c ≤ ∑ c, g c ^ 2 * ℓ c) :
    -∑ c, g c ^ 2 * ℓ c ≤ firstVariation g ℓ σ := by
  have h1 : 0 ≤ ∑ c, (g c + σ c) ^ 2 * ℓ c :=
    Finset.sum_nonneg fun c _ => mul_nonneg (sq_nonneg _) (hℓ c)
  have h2 : ∑ c, (g c + σ c) ^ 2 * ℓ c =
      ∑ c, g c ^ 2 * ℓ c + 2 * firstVariation g ℓ σ + ∑ c, σ c ^ 2 * ℓ c := by
    unfold firstVariation
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun c _ => by ring
  linarith

/-- **The speed of SPEC §3.7 is a descent direction for the distortion terms.** With the region
densities `G_0`, `G_1` at the cut vertices, the first variation of the distortion terms is
`Σ_c (G_0 − G_1)_c σ_c ℓ_c` (SPEC §1); the flow moves the cut with `σ_c = −s_c`,
`s_c = cutSpeed G_0 G_1` (that is `φ_t = s`), and then `dE = −Σ_c s_c² ℓ_c ≤ 0`.

[formula: cut-speed] -/
theorem firstVariation_cutSpeed (G0 G1 : C → ℝ) {ℓ : C → ℝ} (hℓ : ∀ c, 0 ≤ ℓ c) :
    firstVariation (fun c => G0 c - G1 c) ℓ (fun c => -Ddg.cutSpeed (G0 c) (G1 c)) =
        -∑ c, Ddg.cutSpeed (G0 c) (G1 c) ^ 2 * ℓ c ∧
      firstVariation (fun c => G0 c - G1 c) ℓ (fun c => -Ddg.cutSpeed (G0 c) (G1 c)) ≤ 0 := by
  have hσ : ∀ c, Ddg.cutSpeed (G0 c) (G1 c) = -(fun c => -Ddg.cutSpeed (G0 c) (G1 c)) c :=
    fun c => by simp
  have hg : ∀ c, Ddg.cutSpeed (G0 c) (G1 c) = (fun c => G0 c - G1 c) c := fun c => rfl
  exact ⟨firstVariation_levelSet hσ hg, firstVariation_levelSet_nonpos hℓ hσ hg⟩

/-! ## Explicit step -/

/-- The oracle's `explicitStep` at `ℝ`: `φ̃ = φ + τ σ`. -/
theorem explicitStep_eq (φ τ σ : ℝ) : Ddg.explicitStep φ τ σ = φ + τ * σ := rfl

/-- Forward Euler: `φ̃ − φ = τ σ`.

[formula: cut-explicit-step] -/
theorem explicitStep_sub (φ τ σ : ℝ) : Ddg.explicitStep φ τ σ - φ = τ * σ := by
  rw [explicitStep_eq]
  ring

/-- **The explicit step decreases the distortion terms to first order.** At the cut vertices the
step changes the level-set function by `φ̃_c − φ_c = τ s_c`, i.e. it displaces the cut by
`−τ s_c` along the outward normal of region 0 (`φ_t = −σ`). The first variation of an energy
with gradient density `s` for this displacement is `−τ Σ_c s_c² ℓ_c ≤ 0`.

[formula: cut-explicit-step] -/
theorem firstVariation_explicitStep (φ s : C → ℝ) {ℓ : C → ℝ} (hℓ : ∀ c, 0 ≤ ℓ c) {τ : ℝ}
    (hτ : 0 ≤ τ) :
    firstVariation s ℓ (fun c => -(Ddg.explicitStep (φ c) τ (s c) - φ c)) =
        -(τ * ∑ c, s c ^ 2 * ℓ c) ∧
      firstVariation s ℓ (fun c => -(Ddg.explicitStep (φ c) τ (s c) - φ c)) ≤ 0 := by
  have h1 : firstVariation s ℓ (fun c => -(Ddg.explicitStep (φ c) τ (s c) - φ c)) =
      -(τ * ∑ c, s c ^ 2 * ℓ c) := by
    unfold firstVariation
    rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun c _ => ?_
    simp only [explicitStep_sub]
    ring
  refine ⟨h1, ?_⟩
  rw [h1, neg_nonpos]
  exact mul_nonneg hτ (Finset.sum_nonneg fun c _ => mul_nonneg (sq_nonneg _) (hℓ c))

end FirstVariation

/-! ## Length step -/

section LengthStep

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- A linear system with an invertible matrix has exactly one solution. -/
theorem existsUnique_mulVec_eq {M : Matrix V V ℝ} (hdet : IsUnit M.det) (b : V → ℝ) :
    ∃! x : V → ℝ, M *ᵥ x = b := by
  refine ⟨M⁻¹ *ᵥ b, ?_, fun y hy => ?_⟩
  · show M *ᵥ (M⁻¹ *ᵥ b) = b
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hdet, Matrix.one_mulVec]
  · rw [← hy, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet, Matrix.one_mulVec]

/-- The quadratic form of `diag(A) + τ L`. -/
theorem dotProduct_diagonal_add_smul_mulVec (A : V → ℝ) (L : Matrix V V ℝ) (τ : ℝ) (x : V → ℝ) :
    x ⬝ᵥ ((Matrix.diagonal A + τ • L) *ᵥ x) = ∑ i, A i * x i ^ 2 + τ * (x ⬝ᵥ (L *ᵥ x)) := by
  rw [Matrix.add_mulVec, dotProduct_add, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul]
  congr 1
  simp only [dotProduct, Matrix.mulVec_diagonal]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- **The matrix of the length step is positive definite**: `W_V + τ L_w` for positive vertex
areas, a positive semidefinite `L_w` and `τ ≥ 0` (SPEC §3.10).

[formula: cut-length-step] -/
theorem lengthStep_posDef {A : V → ℝ} (hA : ∀ i, 0 < A i) {L : Matrix V V ℝ}
    (hL : L.PosSemidef) {τ : ℝ} (hτ : 0 ≤ τ) : (Matrix.diagonal A + τ • L).PosDef := by
  refine Matrix.PosDef.of_dotProduct_mulVec_pos ?_ fun x hx => ?_
  · exact Matrix.isHermitian_iff_isSymm.mpr
      ((Matrix.isSymm_diagonal A).add (hL.1.isSymm.smul τ))
  · rw [star_trivial, dotProduct_diagonal_add_smul_mulVec]
    have h1 : 0 ≤ x ⬝ᵥ (L *ᵥ x) := by
      have := hL.dotProduct_mulVec_nonneg x
      rwa [star_trivial] at this
    have h2 : 0 < ∑ i, A i * x i ^ 2 := by
      obtain ⟨i, hi⟩ : ∃ i, x i ≠ 0 := by
        by_contra hcon
        push Not at hcon
        exact hx (funext hcon)
      exact Finset.sum_pos' (fun j _ => mul_nonneg (hA j).le (sq_nonneg _))
        ⟨i, Finset.mem_univ i, mul_pos (hA i) (by positivity)⟩
    have := mul_nonneg hτ h1
    linarith

/-- **The length step has exactly one solution**: `(W_V + τ L_w) φ′ = W_V φ̃`.

[formula: cut-length-step] -/
theorem lengthStep_existsUnique {A : V → ℝ} (hA : ∀ i, 0 < A i) {L : Matrix V V ℝ}
    (hL : L.PosSemidef) {τ : ℝ} (hτ : 0 ≤ τ) (φt : V → ℝ) :
    ∃! φ' : V → ℝ, (Matrix.diagonal A + τ • L) *ᵥ φ' = Matrix.diagonal A *ᵥ φt :=
  existsUnique_mulVec_eq ((Matrix.isUnit_iff_isUnit_det _).mp (lengthStep_posDef hA hL hτ).isUnit) _

omit [Fintype V] in
/-- The entries of the matrix of the length step are the oracle's `lengthStepEntry`. -/
theorem lengthStep_entry (A : V → ℝ) (L : Matrix V V ℝ) (τ : ℝ) (i j : V) :
    (Matrix.diagonal A + τ • L) i j = Ddg.lengthStepEntry (decide (i = j)) τ (A i) (L i j) := by
  unfold Ddg.lengthStepEntry
  by_cases h : i = j
  · subst h
    simp
  · simp [h]

/-- Row `i` of the length step: `A_i (φ̃_i − φ′_i) = τ (L_w φ′)_i`. -/
theorem lengthStep_row {A : V → ℝ} {L : Matrix V V ℝ} {τ : ℝ} {φt φ' : V → ℝ}
    (h : (Matrix.diagonal A + τ • L) *ᵥ φ' = Matrix.diagonal A *ᵥ φt) (i : V) :
    A i * (φt i - φ' i) = τ * (L *ᵥ φ') i := by
  have := congrFun h i
  simp only [Matrix.add_mulVec, Matrix.smul_mulVec, Pi.add_apply, Pi.smul_apply,
    Matrix.mulVec_diagonal, smul_eq_mul] at this
  linarith

/-- **Backward Euler does not increase the energy**: the solution of
`(W_V + τ L_w) φ′ = W_V φ̃` has `φ′ᵀ L_w φ′ ≤ φ̃ᵀ L_w φ̃`.

[formula: cut-length-step] -/
theorem lengthStep_energy_le {A : V → ℝ} (hA : ∀ i, 0 < A i) {L : Matrix V V ℝ}
    (hL : L.PosSemidef) {τ : ℝ} (hτ : 0 ≤ τ) {φt φ' : V → ℝ}
    (h : (Matrix.diagonal A + τ • L) *ᵥ φ' = Matrix.diagonal A *ᵥ φt) :
    φ' ⬝ᵥ (L *ᵥ φ') ≤ φt ⬝ᵥ (L *ᵥ φt) := by
  have hrow := lengthStep_row h
  obtain ⟨d, hd⟩ : ∃ d : V → ℝ, d = φt - φ' := ⟨_, rfl⟩
  have hφt : φt = φ' + d := by
    funext i
    simp [hd]
  have hsym : φ' ⬝ᵥ (L *ᵥ d) = d ⬝ᵥ (L *ᵥ φ') := dotProduct_mulVec_symm hL.1.isSymm _ _
  have hexp : φt ⬝ᵥ (L *ᵥ φt) =
      φ' ⬝ᵥ (L *ᵥ φ') + 2 * (d ⬝ᵥ (L *ᵥ φ')) + d ⬝ᵥ (L *ᵥ d) := by
    rw [hφt, Matrix.mulVec_add, add_dotProduct, dotProduct_add, dotProduct_add, hsym]
    ring
  have hdd : 0 ≤ d ⬝ᵥ (L *ᵥ d) := by
    have := hL.dotProduct_mulVec_nonneg d
    rwa [star_trivial] at this
  have hcross : τ * (d ⬝ᵥ (L *ᵥ φ')) = ∑ i, A i * d i ^ 2 := by
    simp only [dotProduct, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    have hi := hrow i
    rw [hd]
    simp only [Pi.sub_apply]
    rw [mul_left_comm, ← hi]
    ring
  have hsq : 0 ≤ ∑ i, A i * d i ^ 2 :=
    Finset.sum_nonneg fun i _ => mul_nonneg (hA i).le (sq_nonneg _)
  rcases hτ.eq_or_lt with h0 | hpos
  · have hd0 : d = 0 := by
      funext i
      have hi := hrow i
      rw [← h0, zero_mul] at hi
      have := (mul_eq_zero.mp hi).resolve_left (hA i).ne'
      rw [hd]
      simpa using this
    rw [hexp, hd0]
    simp
  · have hc : 0 ≤ d ⬝ᵥ (L *ᵥ φ') :=
      nonneg_of_mul_nonneg_right (by rw [hcross]; exact hsq) hpos
    linarith

/-- **Backward Euler does not increase the mass-weighted norm**: `Σ A_i φ′_i² ≤ Σ A_i φ̃_i²`.

[formula: cut-length-step] -/
theorem lengthStep_mass_le {A : V → ℝ} (hA : ∀ i, 0 ≤ A i) {L : Matrix V V ℝ}
    (hL : L.PosSemidef) {τ : ℝ} (hτ : 0 ≤ τ) {φt φ' : V → ℝ}
    (h : (Matrix.diagonal A + τ • L) *ᵥ φ' = Matrix.diagonal A *ᵥ φt) :
    ∑ i, A i * φ' i ^ 2 ≤ ∑ i, A i * φt i ^ 2 := by
  have hrow := lengthStep_row h
  have h1 : 0 ≤ φ' ⬝ᵥ (L *ᵥ φ') := by
    have := hL.dotProduct_mulVec_nonneg φ'
    rwa [star_trivial] at this
  have h2 : ∑ i, A i * φt i ^ 2 - ∑ i, A i * φ' i ^ 2 =
      2 * (τ * (φ' ⬝ᵥ (L *ᵥ φ'))) + ∑ i, A i * (φt i - φ' i) ^ 2 := by
    simp only [dotProduct]
    rw [← Finset.sum_sub_distrib, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    have hi := hrow i
    have : φ' i * (τ * (L *ᵥ φ') i) = φ' i * (A i * (φt i - φ' i)) := by rw [hi]
    linear_combination (-2 : ℝ) * this
  have h3 : 0 ≤ ∑ i, A i * (φt i - φ' i) ^ 2 :=
    Finset.sum_nonneg fun i _ => mul_nonneg (hA i) (sq_nonneg _)
  have := mul_nonneg hτ h1
  linarith

/-- **The length step keeps the integral of `φ`**: if the columns of `L_w` sum to zero then
`Σ A_i φ′_i = Σ A_i φ̃_i`. -/
theorem lengthStep_sum {A : V → ℝ} {L : Matrix V V ℝ} (hcol : (fun _ => (1 : ℝ)) ᵥ* L = 0)
    {τ : ℝ} {φt φ' : V → ℝ}
    (h : (Matrix.diagonal A + τ • L) *ᵥ φ' = Matrix.diagonal A *ᵥ φt) :
    ∑ i, A i * φ' i = ∑ i, A i * φt i := by
  have hrow := lengthStep_row h
  have h0 := sum_mulVec_eq_zero hcol φ'
  have h1 : ∑ i, A i * (φt i - φ' i) = τ * ∑ i, (L *ᵥ φ') i := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => hrow i
  rw [h0, mul_zero] at h1
  simp only [mul_sub, Finset.sum_sub_distrib] at h1
  linarith

/-- **The length step is the minimizing movement of `½ φᵀ L_w φ`.** The solution `φ′` of
`(W_V + τ L_w) φ′ = W_V φ̃` minimizes `Φ(ψ) = ½ Σ A_i (ψ_i − φ̃_i)² + (τ/2) ψᵀ L_w ψ`:
`Φ(ψ) = Φ(φ′) + ½ Σ A_i (ψ_i − φ′_i)² + (τ/2) (ψ − φ′)ᵀ L_w (ψ − φ′)`.

[formula: cut-length-step] -/
theorem lengthStep_minimizes {A : V → ℝ} (hA : ∀ i, 0 ≤ A i) {L : Matrix V V ℝ}
    (hL : L.PosSemidef) {τ : ℝ} (hτ : 0 ≤ τ) {φt φ' : V → ℝ}
    (h : (Matrix.diagonal A + τ • L) *ᵥ φ' = Matrix.diagonal A *ᵥ φt) (ψ : V → ℝ) :
    1 / 2 * ∑ i, A i * (φ' i - φt i) ^ 2 + τ / 2 * (φ' ⬝ᵥ (L *ᵥ φ')) ≤
      1 / 2 * ∑ i, A i * (ψ i - φt i) ^ 2 + τ / 2 * (ψ ⬝ᵥ (L *ᵥ ψ)) := by
  have hrow := lengthStep_row h
  obtain ⟨e, he⟩ : ∃ e : V → ℝ, e = ψ - φ' := ⟨_, rfl⟩
  have hψ : ψ = φ' + e := by
    funext i
    simp [he]
  have hsym : φ' ⬝ᵥ (L *ᵥ e) = e ⬝ᵥ (L *ᵥ φ') := dotProduct_mulVec_symm hL.1.isSymm _ _
  have hexp : ψ ⬝ᵥ (L *ᵥ ψ) =
      φ' ⬝ᵥ (L *ᵥ φ') + 2 * (e ⬝ᵥ (L *ᵥ φ')) + e ⬝ᵥ (L *ᵥ e) := by
    rw [hψ, Matrix.mulVec_add, add_dotProduct, dotProduct_add, dotProduct_add, hsym]
    ring
  have hee : 0 ≤ e ⬝ᵥ (L *ᵥ e) := by
    have := hL.dotProduct_mulVec_nonneg e
    rwa [star_trivial] at this
  have hcross : τ * (e ⬝ᵥ (L *ᵥ φ')) = ∑ i, A i * (φt i - φ' i) * e i := by
    simp only [dotProduct, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [hrow i]
    ring
  have hsq : ∑ i, A i * (ψ i - φt i) ^ 2 =
      ∑ i, A i * (φ' i - φt i) ^ 2 - 2 * ∑ i, A i * (φt i - φ' i) * e i +
        ∑ i, A i * e i ^ 2 := by
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [hψ]
    simp only [Pi.add_apply]
    ring
  have h3 : 0 ≤ ∑ i, A i * e i ^ 2 :=
    Finset.sum_nonneg fun i _ => mul_nonneg (hA i) (sq_nonneg _)
  have h4 := mul_nonneg hτ hee
  rw [hexp, hsq]
  nlinarith

/-! ## The paper's penalty rule -/

omit [Fintype V] in
/-- The entries of the matrix of the paper's penalty rule are the oracle's `penaltyStepEntry`. -/
theorem penaltyStep_entry (A d : V → ℝ) (L : Matrix V V ℝ) (τ : ℝ) (i j : V) :
    (Matrix.diagonal (fun i => A i / d i) + τ • L) i j =
      Ddg.penaltyStepEntry (decide (i = j)) τ (A i) (d i) (L i j) := by
  unfold Ddg.penaltyStepEntry
  by_cases h : i = j
  · subst h
    simp
  · simp [h]

/-- **The paper's penalty rule in two forms** (Sec. 7.1): for a positive diagonal
`D = diag(d)`, `(W_V D⁻¹ + τ L) φ′ = W_V D⁻¹ φ̃` is equivalent to
`(W_V + τ D L) φ′ = W_V φ̃`.

[formula: cut-penalty-step-paper] -/
theorem penaltyStep_iff {A d : V → ℝ} (hd : ∀ i, 0 < d i) (L : Matrix V V ℝ) (τ : ℝ)
    (φt φ' : V → ℝ) :
    (Matrix.diagonal (fun i => A i / d i) + τ • L) *ᵥ φ' =
        Matrix.diagonal (fun i => A i / d i) *ᵥ φt ↔
      (Matrix.diagonal A + τ • (Matrix.diagonal d * L)) *ᵥ φ' = Matrix.diagonal A *ᵥ φt := by
  simp only [funext_iff, Matrix.add_mulVec, Matrix.smul_mulVec, ← Matrix.mulVec_mulVec,
    Pi.add_apply, Pi.smul_apply, Matrix.mulVec_diagonal, smul_eq_mul]
  refine forall_congr' fun i => ?_
  have hk : A i = A i / d i * d i := by
    have := (hd i).ne'
    field_simp
  constructor
  · intro h
    rw [hk]
    linear_combination d i * h
  · intro h
    rw [hk] at h
    exact mul_left_cancel₀ (hd i).ne' (by linear_combination h)

/-- **The paper's penalty rule has exactly one solution** for positive areas, positive `d`, a
positive semidefinite `L` and `τ ≥ 0`.

[formula: cut-penalty-step-paper] -/
theorem penaltyStep_existsUnique {A d : V → ℝ} (hA : ∀ i, 0 < A i) (hd : ∀ i, 0 < d i)
    {L : Matrix V V ℝ} (hL : L.PosSemidef) {τ : ℝ} (hτ : 0 ≤ τ) (φt : V → ℝ) :
    ∃! φ' : V → ℝ, (Matrix.diagonal (fun i => A i / d i) + τ • L) *ᵥ φ' =
      Matrix.diagonal (fun i => A i / d i) *ᵥ φt :=
  lengthStep_existsUnique (fun i => div_pos (hA i) (hd i)) hL hτ φt

/-- **For a constant penalty the paper's rule is the default step**: with `d_i = d > 0` it is
`(W_V + τ (d L)) φ′ = W_V φ̃`, the length step with `L_w = d·L` (the weighted Laplacian of the
constant weight `d`, `TriMesh.weightedLaplacian_const`). The two rules differ only through the
variation of `η` (SPEC D8).

[formula: cut-penalty-step-paper] -/
theorem penaltyStep_const_iff {A : V → ℝ} {d : ℝ} (hd : 0 < d) (L : Matrix V V ℝ) (τ : ℝ)
    (φt φ' : V → ℝ) :
    (Matrix.diagonal (fun i => A i / d) + τ • L) *ᵥ φ' =
        Matrix.diagonal (fun i => A i / d) *ᵥ φt ↔
      (Matrix.diagonal A + τ • (d • L)) *ᵥ φ' = Matrix.diagonal A *ᵥ φt := by
  rw [penaltyStep_iff (d := fun _ => d) (fun _ => hd) L τ φt φ']
  have : Matrix.diagonal (fun _ : V => d) * L = d • L := by
    ext i j
    simp [Matrix.diagonal_mul]
  rw [this]

end LengthStep

/-! ## Symmetry projection -/

section Symmetry

variable {V κ : Type*} [Fintype V] [DecidableEq κ]

/-- Number of vertices in the class of `v`. -/
def classSize (cls : V → κ) (v : V) : ℕ := (Finset.univ.filter fun w => cls w = cls v).card

/-- The arithmetic mean of `φ` over the class of `v` (SPEC §3.11). -/
def classMean (cls : V → κ) (φ : V → ℝ) (v : V) : ℝ :=
  (∑ w ∈ Finset.univ.filter (fun w => cls w = cls v), φ w) / classSize cls v

theorem classSize_pos (cls : V → κ) (v : V) : 0 < classSize cls v :=
  Finset.card_pos.mpr ⟨v, by simp⟩

theorem classSize_eq (cls : V → κ) {v w : V} (h : cls v = cls w) :
    classSize cls v = classSize cls w := by
  unfold classSize
  rw [h]

/-- The class mean is constant on classes. -/
theorem classMean_eq_of_class (cls : V → κ) (φ : V → ℝ) {v w : V} (h : cls v = cls w) :
    classMean cls φ v = classMean cls φ w := by
  unfold classMean classSize
  rw [h]

/-- The class mean as a sum over all vertices. -/
theorem classMean_eq_sum (cls : V → κ) (φ : V → ℝ) (v : V) :
    classMean cls φ v = ∑ w, if cls w = cls v then φ w / classSize cls v else 0 := by
  unfold classMean
  rw [Finset.sum_div, Finset.sum_filter]

/-- **Class-constant vectors are fixed.**

[formula: cut-symmetry-projection] -/
theorem classMean_of_classConstant (cls : V → κ) {ψ : V → ℝ}
    (hψ : ∀ v w, cls v = cls w → ψ v = ψ w) : classMean cls ψ = ψ := by
  funext v
  unfold classMean
  have h1 : ∑ w ∈ Finset.univ.filter (fun w => cls w = cls v), ψ w =
      (classSize cls v : ℝ) * ψ v := by
    rw [Finset.sum_congr rfl fun w hw => hψ w v (Finset.mem_filter.mp hw).2]
    simp [classSize]
  have h2 : (classSize cls v : ℝ) ≠ 0 := by exact_mod_cast (classSize_pos cls v).ne'
  rw [h1]
  field_simp

/-- **The fixed points are exactly the class-constant vectors.**

[formula: cut-symmetry-projection] -/
theorem classMean_eq_self_iff (cls : V → κ) (ψ : V → ℝ) :
    classMean cls ψ = ψ ↔ ∀ v w, cls v = cls w → ψ v = ψ w := by
  constructor
  · intro h v w hvw
    rw [← h]
    exact classMean_eq_of_class cls ψ hvw
  · exact classMean_of_classConstant cls

/-- **The class mean is idempotent.**

[formula: cut-symmetry-projection] -/
theorem classMean_idem (cls : V → κ) (φ : V → ℝ) :
    classMean cls (classMean cls φ) = classMean cls φ :=
  classMean_of_classConstant cls fun _ _ h => classMean_eq_of_class cls φ h

/-- **The class mean is self-adjoint** for the plain dot product:
`Σ_v (Pφ)_v χ_v = Σ_v φ_v (Pχ)_v`.

[formula: cut-symmetry-projection] -/
theorem classMean_selfAdjoint (cls : V → κ) (φ χ : V → ℝ) :
    ∑ v, classMean cls φ v * χ v = ∑ v, φ v * classMean cls χ v := by
  simp only [classMean_eq_sum, Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun w _ => Finset.sum_congr rfl fun v _ => ?_
  by_cases h : cls w = cls v
  · rw [ite_eq_left h, ite_eq_left h.symm, classSize_eq cls h]
    ring
  · rw [ite_eq_right h, ite_eq_right (Ne.symm h)]
    ring

/-- **The residual is orthogonal to the class-constant vectors**:
`Σ_v (φ − Pφ)_v ψ_v = 0` for every class-constant `ψ`.

[formula: cut-symmetry-projection] -/
theorem sum_sub_classMean_mul (cls : V → κ) (φ : V → ℝ) {ψ : V → ℝ}
    (hψ : ∀ v w, cls v = cls w → ψ v = ψ w) : ∑ v, (φ v - classMean cls φ v) * ψ v = 0 := by
  simp only [sub_mul, Finset.sum_sub_distrib]
  rw [classMean_selfAdjoint, classMean_of_classConstant cls hψ]
  ring

/-- **The class mean is the orthogonal projection onto the class-constant vectors**: for every
class-constant `ψ`, `Σ(φ − ψ)² = Σ(φ − Pφ)² + Σ(Pφ − ψ)²`. In particular `Pφ` is the
class-constant vector nearest to `φ`.

[formula: cut-symmetry-projection] -/
theorem classMean_pythagoras (cls : V → κ) (φ : V → ℝ) {ψ : V → ℝ}
    (hψ : ∀ v w, cls v = cls w → ψ v = ψ w) :
    ∑ v, (φ v - ψ v) ^ 2 =
      ∑ v, (φ v - classMean cls φ v) ^ 2 + ∑ v, (classMean cls φ v - ψ v) ^ 2 := by
  have hc : ∀ v w, cls v = cls w →
      (fun v => classMean cls φ v - ψ v) v = (fun v => classMean cls φ v - ψ v) w := by
    intro v w h
    simp only
    rw [classMean_eq_of_class cls φ h, hψ v w h]
  have h0 := sum_sub_classMean_mul cls φ hc
  have h1 : ∑ v, (φ v - ψ v) ^ 2 =
      ∑ v, (φ v - classMean cls φ v) ^ 2 +
        2 * ∑ v, (φ v - classMean cls φ v) * (classMean cls φ v - ψ v) +
        ∑ v, (classMean cls φ v - ψ v) ^ 2 := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun v _ => by ring
  rw [h1, h0]
  ring

/-- The class mean is the nearest class-constant vector. -/
theorem classMean_nearest (cls : V → κ) (φ : V → ℝ) {ψ : V → ℝ}
    (hψ : ∀ v w, cls v = cls w → ψ v = ψ w) :
    ∑ v, (φ v - classMean cls φ v) ^ 2 ≤ ∑ v, (φ v - ψ v) ^ 2 := by
  rw [classMean_pythagoras cls φ hψ]
  have : 0 ≤ ∑ v, (classMean cls φ v - ψ v) ^ 2 := Finset.sum_nonneg fun v _ => sq_nonneg _
  linarith

/-- The oracle's `sumList` at `ℝ` is the sum of the list. -/
theorem sumList_eq_sum (l : List ℝ) : Ddg.sumList l = l.sum := by
  induction l with
  | nil => simp [Ddg.sumList]
  | cons x xs ih => simp [Ddg.sumList, ih]

/-- **The oracle's `classMean` is the class mean**, when `members` lists the class of `v` without
repetition and `count` is its length.

[formula: cut-symmetry-projection] -/
theorem oracle_classMean_eq [DecidableEq V] (cls : V → κ) (φ : V → ℝ) (v : V) {members : List V}
    (hnd : members.Nodup) (hmem : ∀ w, w ∈ members ↔ cls w = cls v) :
    Ddg.classMean (members.map φ) (members.length : ℝ) = classMean cls φ v := by
  have hset : members.toFinset = Finset.univ.filter (fun w => cls w = cls v) := by
    ext w
    simp [hmem w]
  unfold Ddg.classMean classMean classSize
  rw [sumList_eq_sum, ← List.sum_toFinset _ hnd, ← List.toFinset_card_of_nodup hnd, hset]

end Symmetry

/-! ## Step size -/

/-- The oracle's `stepCap` at `ℝ`: `min(τ_max, max_motion · h / max|s|)`. -/
theorem stepCap_eq (τmax c h m : ℝ) : Ddg.stepCap τmax c h m = min τmax (c * h / m) := rfl

/-- **The capped step is at most `τ_max`.**

[formula: cut-step-cap] -/
theorem stepCap_le_max (τmax c h m : ℝ) : Ddg.stepCap τmax c h m ≤ τmax := by
  rw [stepCap_eq]
  exact min_le_left _ _

/-- **The capped step moves the level set by at most `max_motion · h`**:
`τ · max|s| ≤ max_motion · h` (for `max|s| > 0`).

[formula: cut-step-cap] -/
theorem stepCap_mul_le (τmax c h : ℝ) {m : ℝ} (hm : 0 < m) :
    Ddg.stepCap τmax c h m * m ≤ c * h := by
  rw [stepCap_eq]
  calc min τmax (c * h / m) * m ≤ c * h / m * m :=
        mul_le_mul_of_nonneg_right (min_le_right _ _) hm.le
    _ = c * h := by field_simp

/-- The capped step is nonnegative for nonnegative data. -/
theorem stepCap_nonneg {τmax c h m : ℝ} (hτ : 0 ≤ τmax) (hc : 0 ≤ c) (hh : 0 ≤ h) (hm : 0 < m) :
    0 ≤ Ddg.stepCap τmax c h m := by
  rw [stepCap_eq]
  exact le_min hτ (div_nonneg (mul_nonneg hc hh) hm.le)

/-- **Every cut vertex moves by at most `max_motion · h`**: if `|s_c| ≤ max|s|` then
`τ |s_c| ≤ max_motion · h`.

[formula: cut-step-cap] -/
theorem stepCap_mul_abs_le {τmax c h m : ℝ} (hτ : 0 ≤ τmax) (hc : 0 ≤ c) (hh : 0 ≤ h)
    (hm : 0 < m) {s : ℝ} (hs : |s| ≤ m) : Ddg.stepCap τmax c h m * |s| ≤ c * h :=
  (mul_le_mul_of_nonneg_left hs (stepCap_nonneg hτ hc hh hm)).trans (stepCap_mul_le τmax c h hm)

/-- **The capped step is the largest admissible one**: every `τ ≤ τ_max` with
`τ · max|s| ≤ max_motion · h` is at most `stepCap`.

[formula: cut-step-cap] -/
theorem le_stepCap {τmax c h m τ : ℝ} (hm : 0 < m) (h1 : τ ≤ τmax) (h2 : τ * m ≤ c * h) :
    τ ≤ Ddg.stepCap τmax c h m := by
  rw [stepCap_eq]
  exact le_min h1 ((le_div_iff₀ hm).mpr h2)

/-! ## Curvature-based weight -/

/-- The oracle's `vertexCurvatureDensity` at `ℝ`: `K_i = (2π − Θ_i)/A_i`. -/
theorem vertexCurvatureDensity_eq (θ A : ℝ) :
    Ddg.vertexCurvatureDensity θ A = (2 * π - θ) / A := by
  simp only [Ddg.vertexCurvatureDensity, Oracle.two_eq, Oracle.pi_eq]

/-- The oracle's `featureSizeRaw` at `ℝ`: `r⁰ = 1/(√|K| + 1/r_ε)` (SPEC D5). -/
theorem featureSizeRaw_eq (rε K : ℝ) : Ddg.featureSizeRaw rε K = 1 / (√|K| + 1 / rε) := by
  simp only [Ddg.featureSizeRaw, Oracle.sqrt_eq, Oracle.abs_eq, Oracle.one_eq]

/-- **The feature size is positive** (`r_ε > 0`).

[formula: cut-feature-size] -/
theorem featureSizeRaw_pos {rε : ℝ} (hr : 0 < rε) (K : ℝ) : 0 < Ddg.featureSizeRaw rε K := by
  rw [featureSizeRaw_eq]
  have := Real.sqrt_nonneg |K|
  positivity

/-- **The feature size is at most `r_ε`.**

[formula: cut-feature-size] -/
theorem featureSizeRaw_le {rε : ℝ} (hr : 0 < rε) (K : ℝ) : Ddg.featureSizeRaw rε K ≤ rε := by
  rw [featureSizeRaw_eq]
  have h1 := Real.sqrt_nonneg |K|
  have h2 : 0 < 1 / rε := by positivity
  rw [div_le_iff₀ (by linarith)]
  have : rε * (1 / rε) = 1 := by field_simp
  nlinarith

/-- In a flat region the feature size is `r_ε`. -/
theorem featureSizeRaw_zero (rε : ℝ) : Ddg.featureSizeRaw rε 0 = rε := by
  rw [featureSizeRaw_eq]
  simp

/-- In a curved region the feature size is below the curvature radius `1/√|K|`. -/
theorem featureSizeRaw_lt_inv_sqrt {rε : ℝ} (hr : 0 < rε) {K : ℝ} (hK : K ≠ 0) :
    Ddg.featureSizeRaw rε K < 1 / √|K| := by
  rw [featureSizeRaw_eq]
  have h1 : 0 < √|K| := Real.sqrt_pos.mpr (abs_pos.mpr hK)
  have h2 : 0 < 1 / rε := by positivity
  exact one_div_lt_one_div_of_lt h1 (by linarith)

/-- The feature size decreases with the curvature. -/
theorem featureSizeRaw_antitone {rε : ℝ} (hr : 0 < rε) {K1 K2 : ℝ} (h : |K1| ≤ |K2|) :
    Ddg.featureSizeRaw rε K2 ≤ Ddg.featureSizeRaw rε K1 := by
  rw [featureSizeRaw_eq, featureSizeRaw_eq]
  have h1 := Real.sqrt_nonneg |K1|
  have h2 : 0 < 1 / rε := by positivity
  exact one_div_le_one_div_of_le (by linarith) (by linarith [Real.sqrt_le_sqrt h])

/-- The oracle's `featureWeight` at `ℝ`: `a = r³`; positive for `r > 0`. -/
theorem featureWeight_eq (r : ℝ) : Ddg.featureWeight r = r ^ 3 := by
  unfold Ddg.featureWeight
  ring

theorem featureWeight_pos {r : ℝ} (hr : 0 < r) : 0 < Ddg.featureWeight r := by
  rw [featureWeight_eq]
  positivity

section Smoothing

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- **Discrete maximum principle for the smoothing `(W_V + c L) r = W_V r⁰`.** If the
off-diagonal entries of `L` are nonpositive (true for the cotan Laplacian when all edge weights
`cot α + cot β` are nonnegative, e.g. on a Delaunay mesh), the rows of `L` sum to zero, the
areas are positive and `c ≥ 0`, then a lower bound of `r⁰` is a lower bound of `r`.
Without the sign condition this can fail, and then `a = r³ ≥ 0` (needed for
`weightedLaplacian_faceMean_posSemidef`) is not guaranteed.

[formula: cut-feature-size] -/
theorem smoothing_ge {A : V → ℝ} (hA : ∀ i, 0 < A i) {L : Matrix V V ℝ}
    (hoff : ∀ i j, i ≠ j → L i j ≤ 0) (hrow : L *ᵥ (fun _ => (1 : ℝ)) = 0) {c : ℝ}
    (hc : 0 ≤ c) {r r0 : V → ℝ}
    (h : (Matrix.diagonal A + c • L) *ᵥ r = Matrix.diagonal A *ᵥ r0) {m : ℝ}
    (hm : ∀ i, m ≤ r0 i) (i : V) : m ≤ r i := by
  obtain ⟨i₀, -, hmin⟩ := Finset.exists_min_image Finset.univ r ⟨i, Finset.mem_univ i⟩
  have hLr : (L *ᵥ r) i₀ ≤ 0 := by
    have h1 : (L *ᵥ r) i₀ = ∑ j, L i₀ j * (r j - r i₀) := by
      have h0 := congrFun hrow i₀
      simp only [Matrix.mulVec, dotProduct, mul_one, Pi.zero_apply] at h0
      simp only [Matrix.mulVec, dotProduct, mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul,
        h0, zero_mul, sub_zero]
    rw [h1]
    refine Finset.sum_nonpos fun j _ => ?_
    by_cases hj : j = i₀
    · rw [hj, sub_self, mul_zero]
    · exact mul_nonpos_of_nonpos_of_nonneg (hoff i₀ j (Ne.symm hj))
        (sub_nonneg.mpr (hmin j (Finset.mem_univ j)))
  have hrow₀ := lengthStep_row h i₀
  have h2 : 0 ≤ A i₀ * (r i₀ - r0 i₀) := by
    have := mul_nonneg hc (neg_nonneg.mpr hLr)
    nlinarith
  have h3 : r0 i₀ ≤ r i₀ := by
    have := nonneg_of_mul_nonneg_right h2 (hA i₀)
    linarith
  exact (hm i₀).trans (h3.trans (hmin i (Finset.mem_univ i)))

/-- The same for an upper bound.

[formula: cut-feature-size] -/
theorem smoothing_le {A : V → ℝ} (hA : ∀ i, 0 < A i) {L : Matrix V V ℝ}
    (hoff : ∀ i j, i ≠ j → L i j ≤ 0) (hrow : L *ᵥ (fun _ => (1 : ℝ)) = 0) {c : ℝ}
    (hc : 0 ≤ c) {r r0 : V → ℝ}
    (h : (Matrix.diagonal A + c • L) *ᵥ r = Matrix.diagonal A *ᵥ r0) {m : ℝ}
    (hm : ∀ i, r0 i ≤ m) (i : V) : r i ≤ m := by
  have h' : (Matrix.diagonal A + c • L) *ᵥ (-r) = Matrix.diagonal A *ᵥ (-r0) := by
    rw [Matrix.mulVec_neg, Matrix.mulVec_neg, h]
  have := smoothing_ge hA hoff hrow hc h' (m := -m) (fun j => by simpa using hm j) i
  simpa using this

/-- **Smoothed feature size**: on a mesh whose cotan Laplacian has nonpositive off-diagonal
entries, the solution of `(W_V + c L) r = W_V r⁰` with `r⁰ = featureSizeRaw r_ε K` satisfies
`0 < r_i ≤ r_ε`, hence `a_i = r_i³ > 0`.

[formula: cut-feature-size] -/
theorem smoothedFeatureSize_bounds [Nonempty V] {A : V → ℝ} (hA : ∀ i, 0 < A i)
    {L : Matrix V V ℝ} (hoff : ∀ i j, i ≠ j → L i j ≤ 0)
    (hrow : L *ᵥ (fun _ => (1 : ℝ)) = 0) {c : ℝ} (hc : 0 ≤ c) {rε : ℝ} (hr : 0 < rε)
    (K : V → ℝ) {r : V → ℝ}
    (h : (Matrix.diagonal A + c • L) *ᵥ r =
      Matrix.diagonal A *ᵥ fun i => Ddg.featureSizeRaw rε (K i)) (i : V) :
    0 < r i ∧ r i ≤ rε ∧ 0 < Ddg.featureWeight (r i) := by
  obtain ⟨j₀, -, hmin⟩ := Finset.exists_min_image Finset.univ
    (fun j => Ddg.featureSizeRaw rε (K j)) ⟨i, Finset.mem_univ i⟩
  have h1 : 0 < r i :=
    (featureSizeRaw_pos hr (K j₀)).trans_le
      (smoothing_ge hA hoff hrow hc h (fun j => hmin j (Finset.mem_univ j)) i)
  exact ⟨h1, smoothing_le hA hoff hrow hc h (fun j => featureSizeRaw_le hr (K j)) i,
    featureWeight_pos h1⟩

/-- **The smoothing of SPEC §3.6 (D17) keeps the bounds on every mesh.** `L⁺` is the Laplacian
of the edge weights `max(w_ij, 0)`; with `o_ij = −w_ij` the off-diagonal entries of the cotan
Laplacian its off-diagonal entries are `min(o_ij, 0) ≤ 0`, so the hypotheses of
`smoothedFeatureSize_bounds` hold by construction: the solution of
`(W_V + c L⁺) r = W_V r⁰` has `0 < r_i ≤ r_ε` and `a_i = r_i³ > 0`.

[formula: cut-feature-size] -/
theorem smoothedFeatureSize_bounds_positivePart [Nonempty V] {A : V → ℝ} (hA : ∀ i, 0 < A i)
    (o : V → V → ℝ) {c : ℝ} (hc : 0 ≤ c) {rε : ℝ} (hr : 0 < rε) (K : V → ℝ) {r : V → ℝ}
    (h : (Matrix.diagonal A + c • lapOfOffDiag (fun i j => min (o i j) 0)) *ᵥ r =
      Matrix.diagonal A *ᵥ fun i => Ddg.featureSizeRaw rε (K i)) (i : V) :
    0 < r i ∧ r i ≤ rε ∧ 0 < Ddg.featureWeight (r i) :=
  smoothedFeatureSize_bounds hA (fun i j hij => by simp [lapOfOffDiag, hij])
    (lapOfOffDiag_mulVec_one _) hc hr K h i

end Smoothing

/-! ## Sliver smoothing -/

/-- The oracle's `smoothedSpeed` at `ℝ`: `(ℓ₂ s₁ + ℓ₁ s₂)/(ℓ₁ + ℓ₂)`. -/
theorem smoothedSpeed_eq (l1 l2 s1 s2 : ℝ) :
    Ddg.smoothedSpeed l1 l2 s1 s2 = (l2 * s1 + l1 * s2) / (l1 + l2) := rfl

/-- **The smoothed speed is a convex combination of the two neighbours' speeds**, with the
weights `ℓ₂/(ℓ₁ + ℓ₂)` and `ℓ₁/(ℓ₁ + ℓ₂)` (the nearer neighbour weighs more).

[formula: cut-sliver-smoothing] -/
theorem smoothedSpeed_eq_convex (l1 l2 s1 s2 : ℝ) :
    Ddg.smoothedSpeed l1 l2 s1 s2 = l2 / (l1 + l2) * s1 + l1 / (l1 + l2) * s2 := by
  rw [smoothedSpeed_eq]
  ring

/-- The two weights are nonnegative and sum to `1` (segments of nonnegative length, not both
zero). -/
theorem smoothedSpeed_weights {l1 l2 : ℝ} (h1 : 0 ≤ l1) (h2 : 0 ≤ l2) (h : 0 < l1 + l2) :
    0 ≤ l2 / (l1 + l2) ∧ 0 ≤ l1 / (l1 + l2) ∧ l2 / (l1 + l2) + l1 / (l1 + l2) = 1 := by
  refine ⟨div_nonneg h2 h.le, div_nonneg h1 h.le, ?_⟩
  have := h.ne'
  field_simp
  ring

/-- **The smoothed speed lies between the two neighbours' speeds.**

[formula: cut-sliver-smoothing] -/
theorem smoothedSpeed_mem {l1 l2 : ℝ} (h1 : 0 ≤ l1) (h2 : 0 ≤ l2) (h : 0 < l1 + l2)
    (s1 s2 : ℝ) :
    min s1 s2 ≤ Ddg.smoothedSpeed l1 l2 s1 s2 ∧ Ddg.smoothedSpeed l1 l2 s1 s2 ≤ max s1 s2 := by
  rw [smoothedSpeed_eq, le_div_iff₀ h, div_le_iff₀ h]
  have a1 := mul_le_mul_of_nonneg_left (min_le_left s1 s2) h2
  have a2 := mul_le_mul_of_nonneg_left (min_le_right s1 s2) h1
  have b1 := mul_le_mul_of_nonneg_left (le_max_left s1 s2) h2
  have b2 := mul_le_mul_of_nonneg_left (le_max_right s1 s2) h1
  constructor <;> linarith

/-- **The smoothing is exact for a speed that is linear in arclength**: with the neighbours at
the arclengths `−ℓ₁` and `+ℓ₂` from the vertex and `s(x) = α + β x`, the smoothed value is
`s(0) = α`.

[formula: cut-sliver-smoothing] -/
theorem smoothedSpeed_linear {l1 l2 : ℝ} (h : l1 + l2 ≠ 0) (α β : ℝ) :
    Ddg.smoothedSpeed l1 l2 (α - β * l1) (α + β * l2) = α := by
  rw [smoothedSpeed_eq, div_eq_iff h]
  ring

/-- A constant speed is not changed. -/
theorem smoothedSpeed_const {l1 l2 : ℝ} (h : l1 + l2 ≠ 0) (s : ℝ) :
    Ddg.smoothedSpeed l1 l2 s s = s := by
  have := smoothedSpeed_linear h s 0
  simpa using this

/-- The value does not depend on which neighbour is called the first. -/
theorem smoothedSpeed_swap (l1 l2 s1 s2 : ℝ) :
    Ddg.smoothedSpeed l2 l1 s2 s1 = Ddg.smoothedSpeed l1 l2 s1 s2 := by
  rw [smoothedSpeed_eq, smoothedSpeed_eq]
  ring

/-- The oracle's `aspectRatio` at `ℝ`: `(longest edge)² √3 / (4 · area)`. -/
theorem aspectRatio_eq (longest area : ℝ) :
    Ddg.aspectRatio longest area = longest * longest * √3 / (4 * area) := by
  simp only [Ddg.aspectRatio, Oracle.sqrt_eq, Oracle.three_eq, Oracle.four_eq]

/-- **The aspect ratio of an equilateral triangle is `1`.**

[formula: cut-sliver-smoothing] -/
theorem aspectRatio_equilateral {l : ℝ} (hl : 0 < l) :
    Ddg.aspectRatio l (triArea l l l) = 1 := by
  have h3 : (0 : ℝ) < √3 := Real.sqrt_pos.mpr (by norm_num)
  have hT : triArea l l l = √3 * l ^ 2 / 4 := by
    unfold triArea
    rw [show (l + l + l) * (l - l + l) * (l - l + l) * (l - l + l) = 3 * (l ^ 2) ^ 2 by ring,
      Real.sqrt_mul (by norm_num), Real.sqrt_sq (by positivity)]
  rw [aspectRatio_eq, hT]
  have := h3.ne'
  have := hl.ne'
  field_simp

/-- **The aspect ratio of every triangle is at least `1`** (`c` the longest side): the
threshold `ρ₀` of SPEC §3.7 is compared with a number that is `1` for the equilateral triangle
and larger for every other shape.

[formula: cut-sliver-smoothing] -/
theorem one_le_aspectRatio {a b c : ℝ} (h : IsTriangle a b c) (ha : a ≤ c) (hb : b ≤ c) :
    1 ≤ Ddg.aspectRatio c (triArea a b c) := by
  have hT := triArea_pos h
  have h16 := sixteen_mul_triArea_sq h
  have hc := h.pos_c
  have ha0 := h.pos_a
  have hb0 := h.pos_b
  have h3 : (4 * triArea a b c) ^ 2 ≤ (√3 * c ^ 2) ^ 2 := by
    rw [mul_pow (√3) (c ^ 2) 2, Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 3)]
    have hx : 0 ≤ c ^ 2 - a ^ 2 := by nlinarith
    have hy : 0 ≤ c ^ 2 - b ^ 2 := by nlinarith
    nlinarith [mul_nonneg (sq_nonneg c) hx, mul_nonneg (sq_nonneg c) hy, sq_nonneg (a ^ 2 - b ^ 2)]
  have h4 : 4 * triArea a b c ≤ √3 * c ^ 2 :=
    calc 4 * triArea a b c = √((4 * triArea a b c) ^ 2) := (Real.sqrt_sq (by positivity)).symm
      _ ≤ √((√3 * c ^ 2) ^ 2) := Real.sqrt_le_sqrt h3
      _ = √3 * c ^ 2 := Real.sqrt_sq (by positivity)
  rw [aspectRatio_eq, le_div_iff₀ (by positivity)]
  linarith

end Cuts

end DDG
