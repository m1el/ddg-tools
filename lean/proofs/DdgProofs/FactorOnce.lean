import DdgProofs.Solver

/-!
# The factor-once method: the reduced problem

This file is about `crates/ddg-optimal-cones/src/factor_once/` (`reduced.rs`, `dense.rs`,
`mod.rs`), the solution method that eliminates the scale factor `u` and solves a LASSO-type
problem in the cone vector `μ_S` of the free sites.

Notation (the rustdoc of `reduced.rs`). `A = P.L` is the Laplacian block of the unknowns, `ι` the
unknowns, `w_i = wE_i m_i`, `W = diag w`, `S` the *free* sites (`Site.penalized`), `F` the sites
with a fixed angle (`Site.fixed a`: forbidden `a = 0`, prescribed `a`), `a_F` the vector of
the fixed angles (zero on `S`), `Ω - a_F` the right-hand side `rhs0`.

A *linear solver* `R : (ι → ℝ) →ₗ (ι → ℝ)` stands for "solve the state equation": the columns of
the method are `z_j = R e_j`, and `w0 = R (Ω - a_F)`. Two solvers occur:

* `dirichletSolver A = A⁻¹` (`A` positive definite, no pin);
* `pinnedSolver A w p = P ∘ pad ∘ A_p⁻¹ ∘ restrict` (`A` singular with kernel the constants, `p` the
  pin, `P x = x - (Σ w x / Σ w)·1`).

Contents.

* `reduced-lasso-problem`: `reducedU_eq`, `objective_reducedU` (the objective identity),
  `isMinimizer_iff_reduced_dirichlet`, `isMinimizer_iff_reduced_pinned` (equivalence of the
  canonical problem and the reduced one).
* `reduced-gram-matrix`: `gramMatrix_posSemidef`, `gram_quadratic_pos`, `reducedGradient_*`.
* `reduced-optimality-conditions`: `kkt_of_reducedOptimal_dirichlet`,
  `kkt_of_reducedOptimal_pinned` (the direction the code relies on), converse
  `reducedOptimal_of_isMinimizer_dirichlet` (Dirichlet case only).
* `reduced-sum-constraint-pin`: `pinnedSolver_solves`, `pinnedSolver_weighted_mean`,
  `pinnedSolver_unique`, `distortion_pinned_le`.
-/

set_option linter.unusedSectionVars false

noncomputable section

namespace DdgProofs

open Matrix

namespace Site

/-- Whether the site has a fixed cone angle (forbidden or prescribed): `Reduced::new` computes
`s.fixed`. -/
def isFixed : Site → Bool
  | fixed _ => true
  | penalized .. => false

/-- The fixed angle `a` of a fixed site, `0` for a penalized one (`mu_fixed` in `Reduced::new`). -/
def angle : Site → ℝ
  | fixed a => a
  | penalized .. => 0

theorem cost_of_isFixed {s : Site} (h : s.isFixed = true) (x : ℝ) : s.cost x = 0 := by
  cases s <;> simp_all [isFixed, cost]

theorem angle_of_isFixed_eq_false {s : Site} (h : s.isFixed = false) : s.angle = 0 := by
  cases s <;> simp_all [isFixed, angle]

theorem feasible_of_isFixed {s : Site} (h : s.isFixed = true) (x : ℝ) :
    x ∈ s.feasible ↔ x = s.angle := by
  cases s <;> simp_all [isFixed, angle, feasible]

end Site

namespace ConeProblem

variable {ι : Type*} [Fintype ι] [DecidableEq ι] (P : ConeProblem ι)

/-- The distortion weight `w_i = wE_i m_i` (`Discrete::w`). -/
def weight (i : ι) : ℝ := P.wE i * P.m i

/-- The free sites `S`: the unknowns whose cone angle is a variable (`!s.fixed`). -/
def freeSites : Finset ι := Finset.univ.filter fun i => (P.site i).isFixed = false

/-- The vector `a_F` of fixed angles, zero on the free sites (`mu_fixed`). -/
def fixedAngle (i : ι) : ℝ := (P.site i).angle

/-- The right-hand side `Ω - a_F` (`rhs0`). -/
def rhs0 : ι → ℝ := P.Ω - P.fixedAngle

/-- A cone vector `μ` over all unknowns has the prescribed angle at every fixed site. -/
def FixedAt (μ : ι → ℝ) : Prop := ∀ i, (P.site i).isFixed = true → μ i = (P.site i).angle

/-- The cone vector with the fixed sites removed: `μ_S`, padded with zeros. -/
def freePart (μ : ι → ℝ) : ι → ℝ := fun j => if j ∈ P.freeSites then μ j else 0

/-- `μ` is feasible for the reduced problem without the sum constraint: it has the prescribed
angle at every fixed site and `lo_j ≤ μ_j ≤ hi_j` at every free site. -/
def RedFeasible (μ : ι → ℝ) : Prop :=
  P.FixedAt μ ∧ ∀ j ∈ P.freeSites, μ j ∈ (P.site j).feasible

theorem mem_freeSites {j : ι} : j ∈ P.freeSites ↔ (P.site j).isFixed = false := by
  simp [freeSites]

theorem freePart_add_fixedAngle {μ : ι → ℝ} (h : P.FixedAt μ) :
    P.freePart μ + P.fixedAngle = μ := by
  ext i
  by_cases hi : (P.site i).isFixed = true
  · have : i ∉ P.freeSites := by simp [mem_freeSites, hi]
    simp [freePart, this, fixedAngle, h i hi]
  · have hi' : (P.site i).isFixed = false := by simpa using hi
    have : i ∈ P.freeSites := (mem_freeSites P).2 hi'
    simp [freePart, this, fixedAngle, Site.angle_of_isFixed_eq_false hi']

theorem freePart_sub_rhs0 {μ : ι → ℝ} (h : P.FixedAt μ) : P.freePart μ - P.rhs0 = μ - P.Ω := by
  have := P.freePart_add_fixedAngle h
  ext i
  have h2 := congrFun this i
  simp only [rhs0, Pi.sub_apply, Pi.add_apply] at h2 ⊢
  linarith

theorem coneCost_eq_sum_free (μ : ι → ℝ) :
    P.coneCost μ = ∑ j ∈ P.freeSites, (P.site j).cost (μ j) := by
  unfold coneCost freeSites
  rw [Finset.sum_filter]
  refine Finset.sum_congr rfl fun i _ => ?_
  by_cases hi : (P.site i).isFixed = false
  · simp [hi]
  · have : (P.site i).isFixed = true := by simpa using hi
    simp [hi, Site.cost_of_isFixed this]

/-! ## The reduced quantities of a linear solver -/

section Solver

variable (R : (ι → ℝ) →ₗ[ℝ] (ι → ℝ))

/-- The column `z_j = R e_j`. -/
def column (j : ι) : ι → ℝ := R (Pi.single j 1)

/-- `w0 = R (Ω - a_F)`. -/
def w0 : ι → ℝ := R P.rhs0

/-- The Gram matrix `G_jk = z_jᵀ W z_k` (all pairs; the method uses the free sites). -/
def gramMatrix : Matrix ι ι ℝ := fun j k => ∑ i, P.weight i * column R j i * column R k i

/-- The linear term `q_j = -z_jᵀ W w0`. -/
def linTerm (j : ι) : ℝ := -∑ i, P.weight i * column R j i * w0 P R i

/-- The scale factor `u = Z μ_S - w0` of a cone vector. -/
def reducedU (μ : ι → ℝ) : ι → ℝ :=
  fun i => ∑ j ∈ P.freeSites, μ j * column R j i - w0 P R i

/-- The reduced objective `½ μ_Sᵀ G μ_S + qᵀ μ_S + λ Σ_{j∈S} wR_j |μ_j|`
(the cost is `Site.cost`, i.e. `wR·|·|`, so that `λ·cost = t_j |μ_j|` with `t_j = λ wR_j`). -/
def reducedObjective (μ : ι → ℝ) : ℝ :=
  1 / 2 * ∑ j ∈ P.freeSites, ∑ k ∈ P.freeSites, μ j * gramMatrix P R j k * μ k
    + ∑ j ∈ P.freeSites, linTerm P R j * μ j
    + P.lam * ∑ j ∈ P.freeSites, (P.site j).cost (μ j)

/-- **`u = Z μ_S - w0` is `R (μ_S - (Ω - a_F))`.** -/
theorem reducedU_eq (μ : ι → ℝ) : P.reducedU R μ = R (P.freePart μ - P.rhs0) := by
  have h1 : P.freePart μ = ∑ j, P.freePart μ j • (Pi.single j (1 : ℝ) : ι → ℝ) := by
    ext i; simp [Finset.sum_apply, Pi.single_apply]
  rw [map_sub, h1, map_sum]
  ext i
  simp only [reducedU, w0, column, Pi.sub_apply, Finset.sum_apply, map_smul, Pi.smul_apply,
    smul_eq_mul, freePart]
  simp only [ite_mul, zero_mul]
  rw [Finset.sum_ite_mem, Finset.univ_inter]

/-- Expansion of a weighted square of a combination of columns. -/
theorem weighted_sq_combination (s : Finset ι) (x : ι → ℝ) :
    ∑ i, P.weight i * (∑ j ∈ s, x j * column R j i) ^ 2 =
      ∑ j ∈ s, ∑ k ∈ s, x j * gramMatrix P R j k * x k := by
  have h1 : ∀ i, P.weight i * (∑ j ∈ s, x j * column R j i) ^ 2 =
      ∑ j ∈ s, ∑ k ∈ s, x j * x k * (P.weight i * column R j i * column R k i) := by
    intro i
    rw [sq, Finset.sum_mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun k _ => by ring
  rw [Finset.sum_congr rfl fun i _ => h1 i, Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun k _ => ?_
  simp only [gramMatrix, Finset.mul_sum, Finset.sum_mul]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- **The distortion is a quadratic in `μ_S`**: with `u = Z μ_S - w0`,
`½ Σ w u² = ½ μ_Sᵀ G μ_S + qᵀ μ_S + ½ Σ w w0²`. (`w_i = wE_i m_i`; `distortion` is
`½ Σ wE m u²`.) -/
theorem distortion_reducedU (μ : ι → ℝ) :
    P.distortion (P.reducedU R μ) =
      1 / 2 * ∑ j ∈ P.freeSites, ∑ k ∈ P.freeSites, μ j * gramMatrix P R j k * μ k
        + ∑ j ∈ P.freeSites, linTerm P R j * μ j
        + 1 / 2 * ∑ i, P.weight i * w0 P R i ^ 2 := by
  set a : ι → ℝ := fun i => ∑ j ∈ P.freeSites, μ j * column R j i with ha
  have hu : ∀ i, P.reducedU R μ i = a i - w0 P R i := fun i => rfl
  have e1 : P.distortion (P.reducedU R μ) =
      1 / 2 * ∑ i, P.weight i * (a i - w0 P R i) ^ 2 := by
    simp only [distortion, weight, hu]
  have e2 : ∑ i, P.weight i * (a i - w0 P R i) ^ 2 =
      ∑ i, P.weight i * a i ^ 2 - 2 * ∑ i, P.weight i * a i * w0 P R i
        + ∑ i, P.weight i * w0 P R i ^ 2 := by
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun i _ => by ring
  have e3 : ∑ i, P.weight i * a i * w0 P R i = -∑ j ∈ P.freeSites, linTerm P R j * μ j := by
    have h1 : ∀ i, P.weight i * a i * w0 P R i =
        ∑ j ∈ P.freeSites, μ j * (P.weight i * column R j i * w0 P R i) := by
      intro i
      simp only [ha, Finset.mul_sum, Finset.sum_mul]
      exact Finset.sum_congr rfl fun j _ => by ring
    rw [Finset.sum_congr rfl fun i _ => h1 i, Finset.sum_comm, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [← Finset.mul_sum]
    simp only [linTerm]
    ring
  have e4 := weighted_sq_combination P R P.freeSites μ
  rw [e1, e2, e3, ← e4]
  ring

/-- **The objective identity.** The canonical objective of `(u, μ)`, `u = Z μ_S - w0`, is the reduced objective plus the constant
`½ w0ᵀ W w0`. -/
theorem objective_reducedU (μ : ι → ℝ) :
    P.objective (P.reducedU R μ) μ =
      P.reducedObjective R μ + 1 / 2 * ∑ i, P.weight i * w0 P R i ^ 2 := by
  unfold objective reducedObjective
  rw [distortion_reducedU, coneCost_eq_sum_free P μ]
  ring

/-- **The gradient of the distortion in a free cone angle.** For a free site `j`, moving `μ_j`
by `t` moves `u` by `t z_j`, and the distortion changes exactly by
`t·(z_jᵀ W u) + (t²/2)·G_jj`; so `∂/∂μ_j (½ uᵀWu) = z_jᵀ W u`.

[formula: reduced-optimality-conditions] -/
theorem distortion_reducedU_update {j : ι} (hj : j ∈ P.freeSites) (μ : ι → ℝ) (t : ℝ) :
    P.distortion (P.reducedU R (Function.update μ j (μ j + t))) =
      P.distortion (P.reducedU R μ)
        + t * (column R j ⬝ᵥ (fun i => P.weight i * P.reducedU R μ i))
        + t ^ 2 / 2 * gramMatrix P R j j := by
  have hu : P.reducedU R (Function.update μ j (μ j + t)) =
      fun i => P.reducedU R μ i + t * column R j i := by
    have hf : P.freePart (Function.update μ j (μ j + t)) =
        P.freePart μ + t • (Pi.single j (1 : ℝ) : ι → ℝ) := by
      ext i
      by_cases hi : i = j
      · subst hi; simp [freePart, hj]
      · simp [freePart, hi]
    rw [reducedU_eq, hf, add_sub_right_comm, map_add, map_smul, ← reducedU_eq]
    ext i
    simp [column]
  rw [hu]
  unfold distortion gramMatrix dotProduct weight
  simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => by ring

/-! ### The Gram matrix -/

/-- **`G` is positive semidefinite** when the weights are nonnegative (`G = ZᵀWZ`).

[formula: reduced-gram-matrix] -/
theorem gramMatrix_posSemidef (hw : ∀ i, 0 ≤ P.weight i) : (P.gramMatrix R).PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ?_ fun x => ?_
  · rw [Matrix.isHermitian_iff_isSymm]
    ext j k
    simp only [gramMatrix, Matrix.transpose_apply]
    exact Finset.sum_congr rfl fun i _ => by ring
  · have : star x ⬝ᵥ (P.gramMatrix R *ᵥ x) =
        ∑ j, ∑ k, x j * gramMatrix P R j k * x k := by
      simp only [star_trivial, dotProduct, Matrix.mulVec, Finset.mul_sum]
      exact Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => by ring
    rw [this, ← weighted_sq_combination P R Finset.univ x]
    exact Finset.sum_nonneg fun i _ => mul_nonneg (hw i) (sq_nonneg _)

/-- **`G` is positive definite on any set of linearly independent columns**: if the columns
`z_j`, `j ∈ T`, are linearly independent and all weights are positive, then
`xᵀ G x > 0` for every nonzero `x` supported on `T`.

[formula: reduced-gram-matrix] -/
theorem gram_quadratic_pos (hw : ∀ i, 0 < P.weight i) (T : Finset ι)
    (hT : LinearIndepOn ℝ (column R) (T : Set ι)) (x : ι → ℝ) (hsupp : ∀ j, j ∉ T → x j = 0)
    (hx : x ≠ 0) : 0 < ∑ j, ∑ k, x j * gramMatrix P R j k * x k := by
  rw [← weighted_sq_combination P R Finset.univ x]
  have hsum : ∀ i, ∑ j, x j * column R j i = (∑ j ∈ T, x j • column R j) i := by
    intro i
    simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
    symm
    exact Finset.sum_subset (Finset.subset_univ T) fun j _ hj => by rw [hsupp j hj, zero_mul]
  have hv : ∑ j ∈ T, x j • column R j ≠ 0 := by
    intro h0
    rw [linearIndepOn_finset_iff] at hT
    have := hT x h0
    apply hx
    ext j
    by_cases hj : j ∈ T
    · exact this j hj
    · exact hsupp j hj
  obtain ⟨i₀, hi₀⟩ : ∃ i, (∑ j ∈ T, x j • column R j) i ≠ 0 := by
    by_contra hcon
    push Not at hcon
    exact hv (funext hcon)
  refine Finset.sum_pos' (fun i _ => mul_nonneg (hw i).le (sq_nonneg _)) ⟨i₀, Finset.mem_univ _, ?_⟩
  rw [hsum]
  exact mul_pos (hw i₀) (by positivity)

end Solver

/-! ## The reduced optimality conditions -/

/-- **The reduced optimality conditions** at a cone vector `μ` with gradient vector `y` and
multiplier `ν` (`ν = 0` without a sum constraint): the fixed sites have their prescribed
angle, and at every free site `j`, `-(y_j + ν) ∈ ∂g_j(μ_j)` with
`g_j(m) = t_j|m| + ι_[lo_j, hi_j](m)`, `t_j = λ wR_j` (`Site.Stationary`: `μ_j` is allowed and
`λ·cost_j(μ_j) + φ·(m - μ_j) ≤ λ·cost_j(m)` for all allowed `m`).

By `Site.stationary_iff_eq_prox` and `isSubgradientWithinAt_abs_iff` this reads, at `μ_j = 0`:
`|y_j + ν| ≤ t_j` (restricted to the signs the bounds allow); at `μ_j = hi_j > 0`:
`y_j + ν + t_j ≤ 0`; at `μ_j = lo_j < 0`: `y_j + ν - t_j ≥ 0`
(`Dense::exact_violation`, `violation_at_zero`). -/
structure ReducedOptimal (μ y : ι → ℝ) (ν : ℝ) : Prop where
  fixed : P.FixedAt μ
  free : ∀ j ∈ P.freeSites, (P.site j).Stationary P.lam (-(y j + ν)) (μ j)

/-- A fixed site satisfies the inclusion for every `φ` as soon as its angle is the prescribed one. -/
theorem stationary_of_fixed {s : Site} (lam φ μ : ℝ) (hs : s.isFixed = true) (h : μ = s.angle) :
    s.Stationary lam φ μ := by
  cases s with
  | penalized wR lo hi => simp [Site.isFixed] at hs
  | fixed a => exact (Site.stationary_fixed_iff lam a φ μ).2 (by simpa [Site.angle] using h)

theorem redFeasible_iff (μ : ι → ℝ) :
    P.RedFeasible μ ↔ ∀ i, μ i ∈ (P.site i).feasible := by
  constructor
  · rintro ⟨hf, hfree⟩ i
    by_cases hi : (P.site i).isFixed = true
    · rw [Site.feasible_of_isFixed hi]; exact hf i hi
    · exact hfree i ((P.mem_freeSites).2 (by simpa using hi))
  · intro h
    exact ⟨fun i hi => (Site.feasible_of_isFixed hi _).1 (h i), fun j _ => h j⟩

/-- In the reduced optimality conditions the cone vector is feasible. -/
theorem ReducedOptimal.redFeasible {μ y : ι → ℝ} {ν : ℝ} (h : P.ReducedOptimal μ y ν) :
    P.RedFeasible μ :=
  ⟨h.fixed, fun j hj => (h.free j hj).1⟩

theorem sum_freePart (μ : ι → ℝ) : ∑ i, P.freePart μ i = ∑ j ∈ P.freeSites, μ j := by
  simp only [freePart]
  rw [Finset.sum_ite_mem, Finset.univ_inter]

/-- The sum constraint `Σ_S μ_S = Σ (Ω - a_F)` says that `μ - Ω` sums to zero. -/
theorem sum_sub_defect_iff {μ : ι → ℝ} (h : P.FixedAt μ) :
    ∑ i, (μ i - P.Ω i) = 0 ↔ ∑ j ∈ P.freeSites, μ j = ∑ i, P.rhs0 i := by
  have h1 := congrArg (fun f : ι → ℝ => ∑ i, f i) (P.freePart_add_fixedAngle h)
  simp only [Pi.add_apply, Finset.sum_add_distrib, sum_freePart] at h1
  simp only [rhs0, Pi.sub_apply, Finset.sum_sub_distrib]
  constructor <;> intro h2 <;> linarith

/-! ## The Dirichlet case: `A` positive definite, no pin -/

section Dirichlet

/-- The solver `A⁻¹`. -/
def dirichletSolver (A : Matrix ι ι ℝ) : (ι → ℝ) →ₗ[ℝ] (ι → ℝ) := Matrix.mulVecLin A⁻¹

theorem dirichletSolver_apply (A : Matrix ι ι ℝ) (b : ι → ℝ) :
    dirichletSolver A b = A⁻¹ *ᵥ b := rfl

theorem mulVec_dirichletSolver {A : Matrix ι ι ℝ} (hA : IsUnit A.det) (b : ι → ℝ) :
    A *ᵥ dirichletSolver A b = b := by
  rw [dirichletSolver_apply, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv A hA, Matrix.one_mulVec]

/-- **The gradient identity** (Dirichlet): `z_jᵀ W u = y_j` with `y = A⁻¹ (W u)`, for symmetric
`A`. (`z_j = A⁻¹ e_j`.)

[formula: reduced-optimality-conditions] -/
theorem column_dot_dirichlet {A : Matrix ι ι ℝ} (hs : Aᵀ = A) (j : ι)
    (v : ι → ℝ) :
    column (dirichletSolver A) j ⬝ᵥ v = dirichletSolver A v j := by
  have hM : (A⁻¹)ᵀ = A⁻¹ := by rw [Matrix.transpose_nonsing_inv, hs]
  unfold column
  rw [dirichletSolver_apply, dirichletSolver_apply, dotProduct_comm, Matrix.dotProduct_mulVec,
    ← Matrix.mulVec_transpose, hM, dotProduct_comm]
  simp

variable {P}

/-- **The Dirichlet case of the link to the optimality system**
(`reduced-optimality-conditions`; the direction the solver relies on).
Let `A = L` be positive definite, `μ` a cone vector with the reduced optimality conditions for
`ν = 0`, `u = Z μ_S - w0` and `y = A⁻¹(W u)`. Then `(u, φ, μ)` with `φ = -y` solves the optimality
system `KKT` of the canonical problem.

[formula: reduced-optimality-conditions] -/
theorem kkt_of_reducedOptimal_dirichlet (hL : P.L.PosDef) {μ y : ι → ℝ}
    (h : P.ReducedOptimal μ y 0)
    (hy : y = dirichletSolver P.L
      (fun i => P.weight i * P.reducedU (dirichletSolver P.L) μ i)) :
    P.KKT (P.reducedU (dirichletSolver P.L) μ) (fun i => -y i) μ := by
  have hdet : IsUnit P.L.det := (Matrix.isUnit_iff_isUnit_det _).mp hL.isUnit
  have hs : P.Lᵀ = P.L := by simpa using hL.isHermitian.eq
  refine ⟨?_, ?_, fun i => ?_⟩
  · rw [reducedU_eq, mulVec_dirichletSolver hdet, P.freePart_sub_rhs0 h.fixed]
  · rw [hs]
    ext i
    have : P.L *ᵥ (fun i => -y i) = -(P.L *ᵥ y) := by
      rw [← Matrix.mulVec_neg]; rfl
    rw [this, hy, mulVec_dirichletSolver hdet]
    simp [weight]
  · by_cases hi : (P.site i).isFixed = true
    · exact stationary_of_fixed _ _ _ hi (h.fixed i hi)
    · have := h.free i ((P.mem_freeSites).2 (by simpa using hi))
      simpa using this

/-- … hence `(u, μ)` is a global minimizer of the canonical problem
(`optimality-system-sufficient`). -/
theorem isMinimizer_of_reducedOptimal_dirichlet (hP : P.WellFormed) (hL : P.L.PosDef)
    {μ y : ι → ℝ} (h : P.ReducedOptimal μ y 0)
    (hy : y = dirichletSolver P.L
      (fun i => P.weight i * P.reducedU (dirichletSolver P.L) μ i)) :
    P.IsMinimizer (P.reducedU (dirichletSolver P.L) μ) μ :=
  (kkt_of_reducedOptimal_dirichlet hL h hy).isMinimizer hP.wE_nonneg hP.m_nonneg

/-- **Converse (Dirichlet case only).** For a well-formed problem with `L` positive definite,
a minimizer `(u, μ)` has `u = Z μ_S - w0` and satisfies the reduced optimality conditions with
`y = A⁻¹(W u)`, `ν = 0`.

[formula: reduced-optimality-conditions] -/
theorem reducedOptimal_of_isMinimizer_dirichlet (hP : P.WellFormed)
    (hL : P.L.PosDef) {u μ : ι → ℝ} (h : P.IsMinimizer u μ) :
    u = P.reducedU (dirichletSolver P.L) μ ∧
      P.ReducedOptimal μ (dirichletSolver P.L (fun i => P.weight i * u i)) 0 := by
  have hdet : IsUnit P.L.det := (Matrix.isUnit_iff_isUnit_det _).mp hL.isUnit
  have hs : P.Lᵀ = P.L := by simpa using hL.isHermitian.eq
  obtain ⟨φ, hk⟩ := (isMinimizer_iff_exists_kkt_of_posDef hP hL).1 h
  have hfix : P.FixedAt μ := fun i hi =>
    (Site.feasible_of_isFixed hi _).1 (hk.feasible.2 i)
  have hu : u = P.reducedU (dirichletSolver P.L) μ := by
    have := (feasible_iff_of_posDef (P := P) hL).1 hk.feasible
    rw [reducedU_eq, P.freePart_sub_rhs0 hfix, dirichletSolver_apply]
    exact this.1
  refine ⟨hu, hfix, fun j _ => ?_⟩
  have hφ : φ = fun i => -(dirichletSolver P.L (fun i => P.weight i * u i) i + 0) := by
    have h1 : P.L *ᵥ φ = fun i => -(P.weight i * u i) := by
      have := hk.adjoint; rw [hs] at this; exact this
    have h2 : φ = P.L⁻¹ *ᵥ (P.L *ᵥ φ) := by
      rw [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet, Matrix.one_mulVec]
    rw [h2, h1]
    ext i
    simp only [dirichletSolver_apply, add_zero]
    rw [show (fun i => -(P.weight i * u i)) = -(fun i => P.weight i * u i) by ext; simp,
      Matrix.mulVec_neg]
    rfl
  have := hk.inclusion j
  rwa [hφ] at this

/-- **`reduced-lasso-problem`, Dirichlet case.** With `A = L` positive definite,
`(u, μ)` is a minimizer of the canonical problem iff `μ` is feasible (prescribed angles at the
fixed sites, bounds at the free ones), `u = Z μ_S - w0`, and `μ` minimizes the reduced objective
`½ μ_Sᵀ G μ_S + qᵀ μ_S + λ Σ_S wR_j |μ_j|` over all feasible cone vectors.

[formula: reduced-lasso-problem] -/
theorem isMinimizer_iff_reduced_dirichlet (hL : P.L.PosDef) {u μ : ι → ℝ} :
    P.IsMinimizer u μ ↔
      P.RedFeasible μ ∧ u = P.reducedU (dirichletSolver P.L) μ ∧
        ∀ μ', P.RedFeasible μ' →
          P.reducedObjective (dirichletSolver P.L) μ ≤ P.reducedObjective (dirichletSolver P.L) μ' := by
  have hdet : IsUnit P.L.det := (Matrix.isUnit_iff_isUnit_det _).mp hL.isUnit
  have hfe : ∀ μ', P.RedFeasible μ' → P.Feasible (P.reducedU (dirichletSolver P.L) μ') μ' := by
    intro μ' h'
    refine ⟨?_, (P.redFeasible_iff μ').1 h'⟩
    rw [reducedU_eq, mulVec_dirichletSolver hdet, P.freePart_sub_rhs0 h'.1]
  constructor
  · intro h
    have hf := h.1
    have hr : P.RedFeasible μ := (P.redFeasible_iff μ).2 hf.2
    have hu : u = P.reducedU (dirichletSolver P.L) μ := by
      have := (feasible_iff_of_posDef (P := P) hL).1 hf
      rw [reducedU_eq, P.freePart_sub_rhs0 hr.1, dirichletSolver_apply]
      exact this.1
    refine ⟨hr, hu, fun μ' h' => ?_⟩
    have := h.2 _ μ' (hfe μ' h')
    rw [hu, objective_reducedU, objective_reducedU] at this
    linarith
  · rintro ⟨hr, hu, hmin⟩
    subst hu
    refine ⟨hfe μ hr, fun v ν hv => ?_⟩
    have hν : P.RedFeasible ν := (P.redFeasible_iff ν).2 hv.2
    have hvu : v = P.reducedU (dirichletSolver P.L) ν := by
      have := (feasible_iff_of_posDef (P := P) hL).1 hv
      rw [reducedU_eq, P.freePart_sub_rhs0 hν.1, dirichletSolver_apply]
      exact this.1
    rw [hvu, objective_reducedU, objective_reducedU]
    linarith [hmin ν hν]

/-- In the Dirichlet case the columns `z_j = A⁻¹ e_j` are linearly independent, so by
`gram_quadratic_pos` `G` is positive definite on all of `S` when the weights are positive.

[formula: reduced-gram-matrix] -/
theorem linearIndepOn_column_dirichlet {A : Matrix ι ι ℝ} (hA : IsUnit A.det) (T : Finset ι) :
    LinearIndepOn ℝ (column (dirichletSolver A)) (T : Set ι) := by
  rw [linearIndepOn_finset_iff]
  intro f hf i hi
  have h1 : A⁻¹ *ᵥ (∑ j ∈ T, f j • (Pi.single j (1 : ℝ) : ι → ℝ)) = 0 := by
    rw [Matrix.mulVec_sum]
    simpa [column, dirichletSolver_apply, Matrix.mulVec_smul] using hf
  have h2 : (∑ j ∈ T, f j • (Pi.single j (1 : ℝ) : ι → ℝ)) = 0 := by
    have := congrArg (fun x => A *ᵥ x) h1
    simpa [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv A hA] using this
  have := congrFun h2 i
  simpa [Finset.sum_apply, Pi.single_apply, hi] using this

end Dirichlet


/-! ## The closed case: a pin -/

section Pinned

variable (A : Matrix ι ι ℝ) (p : ι)

/-- `A_p`: `A` without row and column `p`. -/
def pinBlock : Matrix {i // i ≠ p} {i // i ≠ p} ℝ := A.submatrix Subtype.val Subtype.val

/-- `pad`: a vector on the unknowns without `p`, extended by `0` at the pin. -/
def padLin : ({i // i ≠ p} → ℝ) →ₗ[ℝ] (ι → ℝ) where
  toFun x i := if h : i = p then 0 else x ⟨i, h⟩
  map_add' x y := by ext i; by_cases h : i = p <;> simp [h]
  map_smul' c x := by ext i; by_cases h : i = p <;> simp [h]

/-- The entries of a vector on the unknowns without `p` (`gather`). -/
def restrictLin : (ι → ℝ) →ₗ[ℝ] ({i // i ≠ p} → ℝ) where
  toFun b k := b k.1
  map_add' x y := by ext; simp
  map_smul' c x := by ext; simp

/-- The pinned solve `b ↦ pad (A_p⁻¹ b_p̂)` (`solve_rows`, `scatter`). -/
def pinSolve : (ι → ℝ) →ₗ[ℝ] (ι → ℝ) :=
  padLin p ∘ₗ Matrix.mulVecLin (pinBlock A p)⁻¹ ∘ₗ restrictLin p

/-- The `w`-weighted mean removal `P x = x - (Σ w x / Σ w)·1` (`Reduced::project`). -/
def meanProj (w : ι → ℝ) : (ι → ℝ) →ₗ[ℝ] (ι → ℝ) where
  toFun x i := x i - (∑ k, w k * x k) / ∑ k, w k
  map_add' x y := by
    ext i
    simp only [Pi.add_apply, mul_add, Finset.sum_add_distrib]
    ring
  map_smul' c x := by
    ext i
    have h : ∑ k, w k * (c * x k) = c * ∑ k, w k * x k := by
      rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun k _ => by ring
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply, h]
    ring

/-- The solver of the closed case, `x = P·pad(A_p⁻¹ b_p̂)`. -/
def pinnedSolver (w : ι → ℝ) : (ι → ℝ) →ₗ[ℝ] (ι → ℝ) := meanProj w ∘ₗ pinSolve A p

theorem padLin_apply (x : {i // i ≠ p} → ℝ) (i : ι) :
    padLin p x i = if h : i = p then 0 else x ⟨i, h⟩ := rfl

theorem restrictLin_apply (b : ι → ℝ) (k : {i // i ≠ p}) : restrictLin p b k = b k.1 := rfl

theorem meanProj_apply (w x : ι → ℝ) (i : ι) :
    meanProj w x i = x i - (∑ k, w k * x k) / ∑ k, w k := rfl

theorem sum_dite_ne (F : {i // i ≠ p} → ℝ) :
    ∑ i, (if h : i = p then (0 : ℝ) else F ⟨i, h⟩) = ∑ k : {i // i ≠ p}, F k := by
  rw [Fintype.sum_eq_add_sum_subtype_ne _ p]
  simp only [dite_true, zero_add]
  exact Finset.sum_congr rfl fun x _ => by simp [x.2]

/-- `padLin` is adjoint to `restrictLin`. -/
theorem padLin_dot (x : {i // i ≠ p} → ℝ) (v : ι → ℝ) :
    padLin p x ⬝ᵥ v = x ⬝ᵥ restrictLin p v := by
  simp only [dotProduct, padLin_apply, restrictLin_apply]
  rw [← sum_dite_ne p (fun k => x k * v k.1)]
  exact Finset.sum_congr rfl fun i _ => by by_cases h : i = p <;> simp [h]

theorem restrictLin_dot (a : ι → ℝ) (x : {i // i ≠ p} → ℝ) :
    restrictLin p a ⬝ᵥ x = a ⬝ᵥ padLin p x := by
  rw [dotProduct_comm, ← padLin_dot, dotProduct_comm]

theorem mulVec_padLin (x : {i // i ≠ p} → ℝ) (i : ι) (hi : i ≠ p) :
    (A *ᵥ padLin p x) i = (pinBlock A p *ᵥ x) ⟨i, hi⟩ := by
  simp only [Matrix.mulVec, dotProduct, padLin_apply, pinBlock, Matrix.submatrix_apply]
  rw [← sum_dite_ne p (fun k => A i k.1 * x k)]
  exact Finset.sum_congr rfl fun j _ => by by_cases h : j = p <;> simp [h]

/-- The hypotheses of a closed mesh: `A` is positive semidefinite, and its kernel is exactly the
constants (`cotan-laplacian-psd`, `cotan-laplacian-kernel`, `laplacian_mulVec_one`). -/
structure KernelIsConstants (A : Matrix ι ι ℝ) : Prop where
  posSemidef : A.PosSemidef
  mulVec_one : A *ᵥ (fun _ => (1 : ℝ)) = 0
  ker : ∀ x, A *ᵥ x = 0 → ∃ c : ℝ, x = fun _ => c

variable {A p}

theorem KernelIsConstants.transpose_eq (h : KernelIsConstants A) : Aᵀ = A := by
  simpa using h.posSemidef.isHermitian.eq

theorem KernelIsConstants.vecMul_one (h : KernelIsConstants A) :
    (fun _ => (1 : ℝ)) ᵥ* A = 0 := by
  rw [← Matrix.mulVec_transpose, h.transpose_eq, h.mulVec_one]

/-- `A_p` is invertible: a positive semidefinite `A` whose kernel is the constants stays definite
once one row and column are removed. -/
theorem KernelIsConstants.pinBlock_isUnit (h : KernelIsConstants A) (p : ι) :
    IsUnit (pinBlock A p).det := by
  rw [← Matrix.isUnit_iff_isUnit_det, ← Matrix.mulVec_injective_iff_isUnit]
  have hinj : Function.Injective (Matrix.mulVecLin (pinBlock A p)) := by
    rw [injective_iff_map_eq_zero]
    intro x hx
    have hx' : pinBlock A p *ᵥ x = 0 := hx
    set X := padLin p x with hX
    have hXp : X p = 0 := by simp [hX, padLin_apply]
    have hAX : ∀ i, i ≠ p → (A *ᵥ X) i = 0 := by
      intro i hi
      rw [hX, mulVec_padLin A p x i hi, hx']; rfl
    have h0 : star X ⬝ᵥ (A *ᵥ X) = 0 := by
      simp only [star_trivial, dotProduct]
      refine Finset.sum_eq_zero fun i _ => ?_
      by_cases hi : i = p
      · subst hi; rw [hXp, zero_mul]
      · rw [hAX i hi, mul_zero]
    have h1 : A *ᵥ X = 0 := h.posSemidef.dotProduct_mulVec_zero_iff.1 h0
    obtain ⟨c, hc⟩ := h.ker X h1
    have hc0 : c = 0 := by rw [← hXp, hc]
    ext k
    have := congrFun hc k.1
    simp only [hc0] at this
    have e : X k.1 = x k := by simp [hX, padLin_apply, k.2]
    rw [← e, this]; rfl
  exact hinj


theorem pinSolve_apply_pin (b : ι → ℝ) : pinSolve A p b p = 0 := by
  simp [pinSolve, padLin_apply]

theorem mulVec_pinSolve_of_ne (h : KernelIsConstants A) (b : ι → ℝ) (i : ι) (hi : i ≠ p) :
    (A *ᵥ pinSolve A p b) i = b i := by
  unfold pinSolve
  simp only [LinearMap.comp_apply]
  rw [mulVec_padLin A p _ i hi, Matrix.mulVecLin_apply, Matrix.mulVec_mulVec,
    Matrix.mul_nonsing_inv _ (h.pinBlock_isUnit p), Matrix.one_mulVec]
  rfl

/-- **The pinned solve solves the singular system** when the right-hand side sums to zero: the
equations of all rows but the pin's hold by construction, the pin's by `1ᵀA = 0` and `Σ b = 0`.

[formula: reduced-sum-constraint-pin] -/
theorem mulVec_pinSolve (h : KernelIsConstants A) {b : ι → ℝ} (hb : ∑ i, b i = 0) :
    A *ᵥ pinSolve A p b = b := by
  funext i
  by_cases hi : i = p
  · subst hi
    exact Solver.mulVec_eq_of_sum_eq h.vecMul_one hb i (fun j hj => mulVec_pinSolve_of_ne h b j hj)
  · exact mulVec_pinSolve_of_ne h b i hi

theorem weighted_meanProj {w : ι → ℝ} (hW : ∑ k, w k ≠ 0) (x : ι → ℝ) :
    ∑ k, w k * meanProj w x k = 0 := by
  simp only [meanProj_apply, mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul]
  field_simp
  ring

theorem mulVec_meanProj (h : KernelIsConstants A) (w x : ι → ℝ) :
    A *ᵥ meanProj w x = A *ᵥ x := by
  have : meanProj w x = x - ((∑ k, w k * x k) / ∑ k, w k) • (fun _ => (1 : ℝ)) := by
    ext i; simp [meanProj_apply]
  rw [this, Matrix.mulVec_sub, Matrix.mulVec_smul, h.mulVec_one]
  simp

/-- **`x = P·pad(A_p⁻¹ b_p̂)` solves `A x = b`** when `Σ b = 0` and `Σ w ≠ 0`.

[formula: reduced-sum-constraint-pin] -/
theorem mulVec_pinnedSolver (h : KernelIsConstants A) (w : ι → ℝ) {b : ι → ℝ}
    (hb : ∑ i, b i = 0) : A *ᵥ pinnedSolver A p w b = b := by
  unfold pinnedSolver
  rw [LinearMap.comp_apply, mulVec_meanProj h, mulVec_pinSolve h hb]

/-- **… and has weighted mean zero**: `Σ w x = 0`.

[formula: reduced-sum-constraint-pin] -/
theorem weighted_pinnedSolver {w : ι → ℝ} (hW : ∑ k, w k ≠ 0) (b : ι → ℝ) :
    ∑ k, w k * pinnedSolver A p w b k = 0 :=
  weighted_meanProj hW _

/-- **The weighted-mean-zero representative is the only one**: if `A v = b` and `Σ w v = 0`
then `v = P·pad(A_p⁻¹ b_p̂)`. Every solution of `A v = b` is this one plus a constant.

[formula: reduced-sum-constraint-pin] -/
theorem pinnedSolver_unique (h : KernelIsConstants A) {w : ι → ℝ} (hW : ∑ k, w k ≠ 0)
    {v b : ι → ℝ} (hv : A *ᵥ v = b) (hmean : ∑ k, w k * v k = 0) :
    v = pinnedSolver A p w b := by
  have hb : ∑ i, b i = 0 := Solver.sum_eq_zero_of_mulVec_eq h.vecMul_one hv
  have hx := mulVec_pinnedSolver (p := p) h w hb
  obtain ⟨c, hc⟩ := h.ker (v - pinnedSolver A p w b) (by rw [Matrix.mulVec_sub, hv, hx, sub_self])
  have hv' : v = pinnedSolver A p w b + fun _ => c := by
    ext i; have := congrFun hc i; simp only [Pi.sub_apply] at this; simp only [Pi.add_apply]; linarith
  have h2 : ∑ k, w k * v k = c * ∑ k, w k := by
    have h3 : ∀ k, w k * v k = w k * pinnedSolver A p w b k + c * w k := by
      intro k; rw [hv']; simp only [Pi.add_apply]; ring
    rw [Finset.sum_congr rfl fun k _ => h3 k, Finset.sum_add_distrib, weighted_pinnedSolver hW,
      zero_add, Finset.mul_sum]
  have hc0 : c = 0 := by
    rw [hmean] at h2
    exact (mul_eq_zero.1 h2.symm).resolve_right hW
  rw [hv', hc0]; ext i; simp

/-- **The gradient identity** (closed case): for symmetric `A` with `A_p` invertible and
`Σ v = 0` (here `v = W u` and `Σ w u = 0`): `z_jᵀ v = (pad (A_p⁻¹ v_p̂))_j`, with
`z_j = P·pad(A_p⁻¹ (e_j)_p̂)`; in particular `0` at the pin.

[formula: reduced-optimality-conditions] -/
theorem column_dot_pinned (hs : Aᵀ = A) (w : ι → ℝ) (j : ι) {v : ι → ℝ} (hv : ∑ i, v i = 0) :
    column (pinnedSolver A p w) j ⬝ᵥ v = pinSolve A p v j := by
  have hM : ((pinBlock A p)⁻¹)ᵀ = (pinBlock A p)⁻¹ := by
    rw [Matrix.transpose_nonsing_inv]
    congr 1
    ext a b
    simpa [pinBlock] using congrFun (congrFun hs a.1) b.1
  have e1 : column (pinnedSolver A p w) j ⬝ᵥ v = pinSolve A p (Pi.single j 1) ⬝ᵥ v := by
    simp only [column, pinnedSolver, LinearMap.comp_apply, dotProduct, meanProj_apply, sub_mul,
      Finset.sum_sub_distrib, ← Finset.mul_sum, hv, mul_zero, sub_zero]
  rw [e1]
  unfold pinSolve
  simp only [LinearMap.comp_apply, Matrix.mulVecLin_apply]
  rw [padLin_dot, dotProduct_comm, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hM,
    dotProduct_comm, restrictLin_dot]
  simp [padLin_apply]

end Pinned

/-! ### The closed case of the reduced problem -/

section Closed

variable {P}

/-- The solver of the closed case for the problem `P` with pin `p`: `x = P·pad(A_p⁻¹ b_p̂)` with
the weights `w_i = wE_i m_i`. -/
abbrev closedSolver (P : ConeProblem ι) (p : ι) : (ι → ℝ) →ₗ[ℝ] (ι → ℝ) :=
  pinnedSolver P.L p P.weight

theorem sum_weight_ne_zero (hw : ∀ i, 0 < P.weight i) (p : ι) : ∑ k, P.weight k ≠ 0 :=
  (Finset.sum_pos (fun i _ => hw i) ⟨p, Finset.mem_univ p⟩).ne'

theorem sum_freePart_sub_rhs0 {μ : ι → ℝ} (hsum : ∑ j ∈ P.freeSites, μ j = ∑ i, P.rhs0 i) :
    ∑ i, (P.freePart μ - P.rhs0) i = 0 := by
  simp only [Pi.sub_apply, Finset.sum_sub_distrib, sum_freePart, hsum, sub_self]

/-- **The state equation of the closed case**: for a feasible cone vector with the sum constraint
`Σ_S μ_S = Σ (Ω - a_F)`, `u = Z μ_S - w0` solves `L u = μ - Ω` and has `Σ w u = 0`.

[formula: reduced-sum-constraint-pin] -/
theorem reducedU_closed_spec (hK : KernelIsConstants P.L) (hw : ∀ i, 0 < P.weight i) (p : ι)
    {μ : ι → ℝ} (hfix : P.FixedAt μ) (hsum : ∑ j ∈ P.freeSites, μ j = ∑ i, P.rhs0 i) :
    P.L *ᵥ P.reducedU (closedSolver P p) μ = μ - P.Ω ∧
      ∑ k, P.weight k * P.reducedU (closedSolver P p) μ k = 0 := by
  rw [reducedU_eq]
  refine ⟨?_, weighted_pinnedSolver (sum_weight_ne_zero hw p) _⟩
  rw [mulVec_pinnedSolver hK _ (sum_freePart_sub_rhs0 hsum), P.freePart_sub_rhs0 hfix]

/-- **The closed case of the link to the optimality system**
(`reduced-optimality-conditions`; the direction the solver relies on).
Let `A = L` be positive semidefinite with kernel the constants, `w > 0`, `p` any index, `μ` a cone
vector with the sum constraint `Σ_S μ_S = Σ (Ω - a_F)` and the reduced optimality conditions for
some multiplier `ν`, where `u = Z μ_S - w0` (with `Z, w0` from `P·pad(A_p⁻¹ ·)`) and
`y = pad(A_p⁻¹ (W u)_p̂)` (so `y_p = 0`). Then `(u, φ, μ)` with `φ = -(y + ν)` solves the
optimality system `KKT` of the canonical problem. The pin need not be a free site for this
direction (the solver chooses a free one).

[formula: reduced-optimality-conditions] -/
theorem kkt_of_reducedOptimal_pinned (hK : KernelIsConstants P.L) (hw : ∀ i, 0 < P.weight i)
    (p : ι) {μ y : ι → ℝ} {ν : ℝ} (hsum : ∑ j ∈ P.freeSites, μ j = ∑ i, P.rhs0 i)
    (h : P.ReducedOptimal μ y ν)
    (hy : y = pinSolve P.L p (fun i => P.weight i * P.reducedU (closedSolver P p) μ i)) :
    P.KKT (P.reducedU (closedSolver P p) μ) (fun i => -(y i + ν)) μ := by
  obtain ⟨hst, hmean⟩ := reducedU_closed_spec hK hw p h.fixed hsum
  refine ⟨hst, ?_, fun i => ?_⟩
  · rw [hK.transpose_eq]
    have e : (fun i => -(y i + ν)) = -y - ν • (fun _ => (1 : ℝ)) := by
      ext i; simp; ring
    rw [e, Matrix.mulVec_sub, Matrix.mulVec_neg, Matrix.mulVec_smul, hK.mulVec_one, hy,
      mulVec_pinSolve hK (b := fun i => P.weight i * P.reducedU (closedSolver P p) μ i) hmean]
    ext i; simp [weight]
  · by_cases hi : (P.site i).isFixed = true
    · exact stationary_of_fixed _ _ _ hi (h.fixed i hi)
    · exact h.free i ((P.mem_freeSites).2 (by simpa using hi))

/-- … hence `(u, μ)` is a global minimizer of the canonical problem. -/
theorem isMinimizer_of_reducedOptimal_pinned (hP : P.WellFormed) (hK : KernelIsConstants P.L)
    (hw : ∀ i, 0 < P.weight i) (p : ι) {μ y : ι → ℝ} {ν : ℝ}
    (hsum : ∑ j ∈ P.freeSites, μ j = ∑ i, P.rhs0 i) (h : P.ReducedOptimal μ y ν)
    (hy : y = pinSolve P.L p (fun i => P.weight i * P.reducedU (closedSolver P p) μ i)) :
    P.IsMinimizer (P.reducedU (closedSolver P p) μ) μ :=
  (kkt_of_reducedOptimal_pinned hK hw p hsum h hy).isMinimizer hP.wE_nonneg hP.m_nonneg

/-- **The gradient of the distortion** (closed case): `∂/∂μ_j (½ uᵀWu) = z_jᵀ W u = y_j`
with `y = pad(A_p⁻¹ (W u)_p̂)`, for every `u` with `Σ w u = 0` (every `u = Z μ_S - w0`).

[formula: reduced-optimality-conditions] -/
theorem column_dot_weighted_pinned (hK : KernelIsConstants P.L) (p : ι) (j : ι) {u : ι → ℝ}
    (hu : ∑ k, P.weight k * u k = 0) :
    column (closedSolver P p) j ⬝ᵥ (fun i => P.weight i * u i) =
      pinSolve P.L p (fun i => P.weight i * u i) j :=
  column_dot_pinned hK.transpose_eq P.weight j hu

/-- **The distortion as a function of the additive constant, on the solutions of `L v = b`.**
Every solution is `x + c·1` with `x = P·pad(A_p⁻¹ b_p̂)`, and
`E(x + c) = E(x) + ½ (Σ w) c²`: the weighted-mean-zero representative minimizes the distortion,
uniquely.

[formula: reduced-sum-constraint-pin] -/
theorem distortion_pinned (hK : KernelIsConstants P.L) (hw : ∀ i, 0 < P.weight i) (p : ι)
    {v b : ι → ℝ} (hv : P.L *ᵥ v = b) :
    ∃ c : ℝ, v = (fun i => closedSolver P p b i + c) ∧
      P.distortion v = P.distortion (closedSolver P p b) + 1 / 2 * (∑ k, P.weight k) * c ^ 2 := by
  have hW := sum_weight_ne_zero hw p
  have hb : ∑ i, b i = 0 := Solver.sum_eq_zero_of_mulVec_eq hK.vecMul_one hv
  obtain ⟨c, hc⟩ := hK.ker (v - closedSolver P p b)
    (by rw [Matrix.mulVec_sub, hv, mulVec_pinnedSolver hK _ hb, sub_self])
  have hv' : v = fun i => closedSolver P p b i + c := by
    ext i; have := congrFun hc i; simp only [Pi.sub_apply] at this; linarith
  refine ⟨c, hv', ?_⟩
  rw [hv', distortion_add_const]
  have : ∑ i, P.wE i * P.m i * closedSolver P p b i = 0 := weighted_pinnedSolver hW b
  rw [this]
  simp only [weight]
  ring

/-- **`reduced-lasso-problem`, closed case.** With `A = L` positive semidefinite with kernel the
constants and `w > 0`, `(u, μ)` is a minimizer of the canonical problem iff `μ` is feasible, has
`Σ_S μ_S = Σ (Ω - a_F)`, `u = Z μ_S - w0`, and `μ` minimizes the reduced objective over all
feasible cone vectors with that sum.

[formula: reduced-lasso-problem] -/
theorem isMinimizer_iff_reduced_pinned (hK : KernelIsConstants P.L) (hw : ∀ i, 0 < P.weight i)
    (p : ι) {u μ : ι → ℝ} :
    P.IsMinimizer u μ ↔
      (P.RedFeasible μ ∧ ∑ j ∈ P.freeSites, μ j = ∑ i, P.rhs0 i) ∧
        u = P.reducedU (closedSolver P p) μ ∧
        ∀ μ', P.RedFeasible μ' → ∑ j ∈ P.freeSites, μ' j = ∑ i, P.rhs0 i →
          P.reducedObjective (closedSolver P p) μ ≤ P.reducedObjective (closedSolver P p) μ' := by
  have hW := sum_weight_ne_zero hw p
  have hfe : ∀ μ', P.RedFeasible μ' → ∑ j ∈ P.freeSites, μ' j = ∑ i, P.rhs0 i →
      P.Feasible (P.reducedU (closedSolver P p) μ') μ' := fun μ' h' hs' =>
    ⟨(reducedU_closed_spec hK hw p h'.1 hs').1, (P.redFeasible_iff μ').1 h'⟩
  have hsum_of : ∀ ν, P.Feasible (P.reducedU (closedSolver P p) ν) ν →
      ∑ j ∈ P.freeSites, ν j = ∑ i, P.rhs0 i := by
    intro ν hf
    have hr : P.FixedAt ν := ((P.redFeasible_iff ν).2 hf.2).1
    rw [← sum_sub_defect_iff P hr]
    have := Solver.sum_eq_zero_of_mulVec_eq hK.vecMul_one hf.1
    exact this
  constructor
  · intro h
    have hf := h.1
    have hr : P.RedFeasible μ := (P.redFeasible_iff μ).2 hf.2
    have hs : ∑ j ∈ P.freeSites, μ j = ∑ i, P.rhs0 i := by
      rw [← sum_sub_defect_iff P hr.1]
      exact Solver.sum_eq_zero_of_mulVec_eq hK.vecMul_one hf.1
    have hmean : ∑ k, P.weight k * u k = 0 := h.weighted_mean_eq_zero hK.mulVec_one
    have hu : u = P.reducedU (closedSolver P p) μ := by
      rw [reducedU_eq, P.freePart_sub_rhs0 hr.1]
      exact pinnedSolver_unique hK hW hf.1 hmean
    refine ⟨⟨hr, hs⟩, hu, fun μ' h' hs' => ?_⟩
    have := h.2 _ μ' (hfe μ' h' hs')
    rw [hu, objective_reducedU, objective_reducedU] at this
    linarith
  · rintro ⟨⟨hr, hs⟩, hu, hmin⟩
    subst hu
    refine ⟨hfe μ hr hs, fun v ν hv => ?_⟩
    have hν : P.RedFeasible ν := (P.redFeasible_iff ν).2 hv.2
    have hsν : ∑ j ∈ P.freeSites, ν j = ∑ i, P.rhs0 i := by
      rw [← sum_sub_defect_iff P hν.1]
      exact Solver.sum_eq_zero_of_mulVec_eq hK.vecMul_one hv.1
    obtain ⟨c, hvc, hD⟩ := distortion_pinned hK hw p hv.1
    have hxν : closedSolver P p (ν - P.Ω) = P.reducedU (closedSolver P p) ν := by
      rw [reducedU_eq, P.freePart_sub_rhs0 hν.1]
    have h1 : P.reducedU (closedSolver P p) ν = closedSolver P p (ν - P.Ω) := hxν.symm
    have hge : P.distortion (P.reducedU (closedSolver P p) ν) ≤ P.distortion v := by
      rw [hD, h1]
      have : 0 ≤ 1 / 2 * (∑ k, P.weight k) * c ^ 2 :=
        mul_nonneg (mul_nonneg (by norm_num) (Finset.sum_nonneg fun i _ => (hw i).le))
          (sq_nonneg _)
      linarith
    have h2 := hmin ν hν hsν
    have h3 := objective_reducedU P (closedSolver P p) μ
    have h4 := objective_reducedU P (closedSolver P p) ν
    have h5 : P.objective (P.reducedU (closedSolver P p) ν) ν ≤ P.objective v ν := by
      unfold objective
      linarith
    linarith

end Closed

end ConeProblem

end DdgProofs
