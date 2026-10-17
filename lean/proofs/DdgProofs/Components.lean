import DdgProofs.Closed

/-!
# The cone problem separates over the connected components of a mesh

On a mesh with several connected components the cotan Laplacian has no entry between two
components: ordering the unknowns component by component, `L` is block diagonal. The objective
`½ Σ wE m u² + λ Σ wR |μ|` is a sum over the vertices, and the constraints (`L u = μ − Ω`, the
per-vertex bounds) are per block. So the problem of the whole mesh is the **sum** of the problems
of its components, with one common `λ` and the masses normalized once, by the total area of the
whole mesh.

This file is about `crates/ddg-optimal-cones/src/problem.rs` (`Discrete::components`, which splits
the discrete problem by connected component) and `crates/ddg-optimal-cones/src/solve.rs` (which
solves each component on its own and assembles the result).

* `ConeProblem.sum P Q`: the problem on `ι ⊕ κ` with `L = fromBlocks P.L 0 0 Q.L`, the data of
  `P` on `ι` and of `Q` on `κ`, and the parameter `λ` of `P`. It is the problem of the mesh made
  of the two pieces when `P.lam = Q.lam` (the hypothesis of every theorem below that involves
  `λ`). Several components are obtained by iterating (`ι ⊕ (κ ⊕ …)`).
* `objective_sum`, `feasible_sum_iff`, `wellFormed_sum`: objective, constraints and admissibility
  split.
* `isMinimizer_sum_iff` (main theorem), `kkt_sum_iff`: a pair is a global minimizer of the whole
  problem iff each restriction minimizes its component's problem; a triple solves the
  optimality system iff each restriction does. On a closed component the competitors are all
  solutions of that component's state equation, so the additive constant of `u` on each closed
  component is a separate free variable.
* `Feasible.sum_cone_eq_sum_defect_inl`, `…_inr`: Gauss–Bonnet per closed component.

The masses. The whole mesh has total area `A`, a component has area `A_c`, and the global
problem uses masses `a_i / A`. Solving the component **alone** with its masses normalized by its
own area, `a_i / A_c = (A / A_c)·(a_i / A)`, is the global component problem with all masses
multiplied by `c = A / A_c`; with `λ_c = c·λ` it has the same minimizers, and its adjoint is
`c` times the global adjoint (`φ_global = (A_c / A)·φ_alone`):

* `Site.stationary_mul_iff`, `objective_scaleMass_withLam`, `isMinimizer_scaleMass_withLam_iff`,
  `kkt_scaleMass_withLam_iff`, `kkt_scaleMass_withLam_iff'`;
* `isMinimizer_normalizeMass_alone_iff`, `kkt_normalizeMass_alone_iff`: the same, written with the
  two normalizations `normalizeMass A_c` and `normalizeMass A` of the raw component problem;
* `isMinimizer_sum_iff_alone`: the whole problem against the two alone-solves.
-/

noncomputable section

namespace DdgProofs

open Matrix

namespace Site

/-- **The inclusion is homogeneous in `(λ, φ)`.** For `c > 0`:
`c·φ ∈ ∂(c·λ·cost + ι)(μ)  ⟺  φ ∈ ∂(λ·cost + ι)(μ)`.

[formula: component-lambda-rescaling] -/
theorem stationary_mul_iff (s : Site) {c : ℝ} (hc : 0 < c) (lam φ μ : ℝ) :
    s.Stationary (c * lam) (c * φ) μ ↔ s.Stationary lam φ μ := by
  unfold Stationary IsSubgradientWithinAt
  refine and_congr Iff.rfl (forall_congr' fun y => imp_congr_right fun _ => ?_)
  dsimp only
  have e1 : c * lam * s.cost μ + c * φ * (y - μ) = c * (lam * s.cost μ + φ * (y - μ)) := by ring
  have e2 : c * lam * s.cost y = c * (lam * s.cost y) := by ring
  rw [e1, e2]
  exact mul_le_mul_iff_of_pos_left hc

end Site

namespace ConeProblem

variable {ι κ : Type*}

/-! ## The sum of two problems -/

/-- The problem on the disjoint union `ι ⊕ κ`: block diagonal `L = fromBlocks P.L 0 0 Q.L`, and
`Ω`, `m`, `wE`, `site` of `P` on `ι` and of `Q` on `κ`. The parameter is `P.lam`; the sum is the
problem of a mesh with the two pieces as its components only when `P.lam = Q.lam`, which the
theorems below assume where `λ` enters. -/
def sum (P : ConeProblem ι) (Q : ConeProblem κ) : ConeProblem (ι ⊕ κ) where
  L := Matrix.fromBlocks P.L 0 0 Q.L
  Ω := Sum.elim P.Ω Q.Ω
  m := Sum.elim P.m Q.m
  wE := Sum.elim P.wE Q.wE
  lam := P.lam
  site := Sum.elim P.site Q.site

/-- Every vector on `ι ⊕ κ` is `Sum.elim` of its two restrictions, so the statements below in
terms of `Sum.elim u v` cover all vectors. -/
theorem eq_sum_elim_comp {α : Type*} (U : ι ⊕ κ → α) :
    U = Sum.elim (U ∘ Sum.inl) (U ∘ Sum.inr) :=
  (Sum.elim_comp_inl_inr U).symm

theorem sum_elim_eq_sum_elim_iff {α : Type*} {f f' : ι → α} {g g' : κ → α} :
    Sum.elim f g = Sum.elim f' g' ↔ f = f' ∧ g = g' := by
  constructor
  · intro h
    exact ⟨funext fun i => congrFun h (Sum.inl i), funext fun j => congrFun h (Sum.inr j)⟩
  · rintro ⟨rfl, rfl⟩; rfl

variable [Fintype ι] [Fintype κ] (P : ConeProblem ι) (Q : ConeProblem κ)

/-- The block diagonal `L` acts per block. -/
theorem sum_L_mulVec (u : ι → ℝ) (v : κ → ℝ) :
    (P.sum Q).L *ᵥ Sum.elim u v = Sum.elim (P.L *ᵥ u) (Q.L *ᵥ v) := by
  simp only [sum, Matrix.fromBlocks_mulVec, Sum.elim_comp_inl, Sum.elim_comp_inr,
    Matrix.zero_mulVec, add_zero, zero_add]

/-- The transpose of the block diagonal `L` acts per block. -/
theorem sum_L_transpose_mulVec (u : ι → ℝ) (v : κ → ℝ) :
    (P.sum Q).Lᵀ *ᵥ Sum.elim u v = Sum.elim (P.Lᵀ *ᵥ u) (Q.Lᵀ *ᵥ v) := by
  simp only [sum, Matrix.fromBlocks_transpose, Matrix.transpose_zero, Matrix.fromBlocks_mulVec,
    Sum.elim_comp_inl, Sum.elim_comp_inr, Matrix.zero_mulVec, add_zero, zero_add]

theorem distortion_sum (u : ι → ℝ) (v : κ → ℝ) :
    (P.sum Q).distortion (Sum.elim u v) = P.distortion u + Q.distortion v := by
  simp only [distortion, sum, Fintype.sum_sum_type, Sum.elim_inl, Sum.elim_inr]
  ring

theorem coneCost_sum (μ : ι → ℝ) (ν : κ → ℝ) :
    (P.sum Q).coneCost (Sum.elim μ ν) = P.coneCost μ + Q.coneCost ν := by
  simp only [coneCost, sum, Fintype.sum_sum_type, Sum.elim_inl, Sum.elim_inr]

/-- **The objective of the sum is the sum of the objectives** (common `λ`).

[formula: component-separation] -/
theorem objective_sum (hlam : P.lam = Q.lam) (u : ι → ℝ) (v : κ → ℝ) (μ : ι → ℝ) (ν : κ → ℝ) :
    (P.sum Q).objective (Sum.elim u v) (Sum.elim μ ν) = P.objective u μ + Q.objective v ν := by
  simp only [objective, distortion_sum, coneCost_sum]
  have : (P.sum Q).lam = P.lam := rfl
  rw [this, mul_add, ← hlam]
  ring

/-- **The constraints are per block**: `(u, v; μ, ν)` is feasible for the sum iff `(u, μ)` is
feasible for `P` and `(v, ν)` for `Q`. No hypothesis on `λ`.

[formula: component-separation] -/
theorem feasible_sum_iff (u : ι → ℝ) (v : κ → ℝ) (μ : ι → ℝ) (ν : κ → ℝ) :
    (P.sum Q).Feasible (Sum.elim u v) (Sum.elim μ ν) ↔ P.Feasible u μ ∧ Q.Feasible v ν := by
  have hrhs : Sum.elim μ ν - (P.sum Q).Ω = Sum.elim (μ - P.Ω) (ν - Q.Ω) := by
    funext i; cases i <;> rfl
  unfold Feasible
  rw [sum_L_mulVec, hrhs, sum_elim_eq_sum_elim_iff, Sum.forall]
  simp only [sum, Sum.elim_inl, Sum.elim_inr]
  tauto

omit [Fintype ι] [Fintype κ] in
/-- The sum of two well-formed problems is well formed (its `λ` is `P.lam`). -/
theorem wellFormed_sum (hP : P.WellFormed) (hQ : Q.WellFormed) :
    (P.sum Q).WellFormed where
  wE_nonneg := Sum.forall.mpr ⟨hP.wE_nonneg, hQ.wE_nonneg⟩
  m_nonneg := Sum.forall.mpr ⟨hP.m_nonneg, hQ.m_nonneg⟩
  lam_nonneg := hP.lam_nonneg
  site := Sum.forall.mpr ⟨hP.site, hQ.site⟩

/-- **The cone problem separates over the connected components.** On a mesh with several
connected components, ordering the unknowns by component makes `L` block diagonal; the objective
is the sum of the objectives of the components (`objective_sum`), and the constraints are per
block (`feasible_sum_iff`). With a common `λ` (and the masses of all components normalized by
the same total area, as in the problem of the whole mesh), a pair is a global minimizer of the
whole problem iff each restriction is a global minimizer of its component's problem. So the
components may be solved separately and the results assembled.

In particular, on every closed component the competitors of that component are all solutions of
its own state equation: the additive constant of `u` on each closed component is a separate free
variable of the minimization (the competitors `(v₁, v₂)` are arbitrary, not tied by one common
constant). No hypothesis besides the common `λ` is needed (no well-formedness, no assumption on
`L`).

[formula: component-separation] -/
theorem isMinimizer_sum_iff (hlam : P.lam = Q.lam) (u : ι → ℝ) (v : κ → ℝ) (μ : ι → ℝ)
    (ν : κ → ℝ) :
    (P.sum Q).IsMinimizer (Sum.elim u v) (Sum.elim μ ν) ↔
      P.IsMinimizer u μ ∧ Q.IsMinimizer v ν := by
  constructor
  · rintro ⟨hf, hmin⟩
    obtain ⟨hfP, hfQ⟩ := (feasible_sum_iff P Q u v μ ν).mp hf
    refine ⟨⟨hfP, fun v' ν' hv' => ?_⟩, ⟨hfQ, fun v' ν' hv' => ?_⟩⟩
    · have h := hmin (Sum.elim v' v) (Sum.elim ν' ν)
        ((feasible_sum_iff P Q v' v ν' ν).mpr ⟨hv', hfQ⟩)
      rw [objective_sum P Q hlam, objective_sum P Q hlam] at h
      linarith
    · have h := hmin (Sum.elim u v') (Sum.elim μ ν')
        ((feasible_sum_iff P Q u v' μ ν').mpr ⟨hfP, hv'⟩)
      rw [objective_sum P Q hlam, objective_sum P Q hlam] at h
      linarith
  · rintro ⟨⟨hfP, hminP⟩, ⟨hfQ, hminQ⟩⟩
    refine ⟨(feasible_sum_iff P Q u v μ ν).mpr ⟨hfP, hfQ⟩, fun W N hW => ?_⟩
    rw [eq_sum_elim_comp W, eq_sum_elim_comp N] at hW ⊢
    obtain ⟨hW1, hW2⟩ := (feasible_sum_iff P Q _ _ _ _).mp hW
    rw [objective_sum P Q hlam, objective_sum P Q hlam]
    linarith [hminP _ _ hW1, hminQ _ _ hW2]

/-- **The optimality system separates over the components** (common `λ`): `(u, v; φ, ψ; μ, ν)`
solves the optimality system of the whole problem iff `(u, φ, μ)` solves that of `P` and
`(v, ψ, ν)` that of `Q`. In particular the adjoint of the whole problem is assembled from the
adjoints of the components.

[formula: component-separation] -/
theorem kkt_sum_iff (hlam : P.lam = Q.lam) (u : ι → ℝ) (v : κ → ℝ) (φ : ι → ℝ) (ψ : κ → ℝ)
    (μ : ι → ℝ) (ν : κ → ℝ) :
    (P.sum Q).KKT (Sum.elim u v) (Sum.elim φ ψ) (Sum.elim μ ν) ↔
      P.KKT u φ μ ∧ Q.KKT v ψ ν := by
  have hrhs : Sum.elim μ ν - (P.sum Q).Ω = Sum.elim (μ - P.Ω) (ν - Q.Ω) := by
    funext i; cases i <;> rfl
  have hadj : (fun i => -((P.sum Q).wE i * (P.sum Q).m i * Sum.elim u v i)) =
      Sum.elim (fun i => -(P.wE i * P.m i * u i)) (fun j => -(Q.wE j * Q.m j * v j)) := by
    funext i; cases i <;> rfl
  constructor
  · intro h
    have hs := h.state
    have ha := h.adjoint
    have hi := h.inclusion
    rw [sum_L_mulVec, hrhs, sum_elim_eq_sum_elim_iff] at hs
    rw [sum_L_transpose_mulVec, hadj, sum_elim_eq_sum_elim_iff] at ha
    rw [Sum.forall] at hi
    refine ⟨⟨hs.1, ha.1, hi.1⟩, ⟨hs.2, ha.2, fun j => ?_⟩⟩
    have := hi.2 j
    rwa [show (P.sum Q).lam = Q.lam from hlam] at this
  · rintro ⟨hP, hQ⟩
    refine ⟨?_, ?_, ?_⟩
    · rw [sum_L_mulVec, hrhs, sum_elim_eq_sum_elim_iff]; exact ⟨hP.state, hQ.state⟩
    · rw [sum_L_transpose_mulVec, hadj, sum_elim_eq_sum_elim_iff]; exact ⟨hP.adjoint, hQ.adjoint⟩
    · rw [Sum.forall]
      refine ⟨hP.inclusion, fun j => ?_⟩
      rw [show (P.sum Q).lam = Q.lam from hlam]
      exact hQ.inclusion j

/-! ## Gauss–Bonnet per closed component -/

/-- **Gauss–Bonnet per component.** If the first component is closed (`1ᵀ P.L = 0`), every
feasible pair of the whole problem has `Σ_{i ∈ ι} μ_i = Σ_{i ∈ ι} Ω_i` on that component
(`2πχ` of the component, by discrete Gauss–Bonnet).

[formula: closed-mesh-cone-sum] -/
theorem Feasible.sum_cone_eq_sum_defect_inl {u : ι → ℝ} {v : κ → ℝ} {μ : ι → ℝ} {ν : κ → ℝ}
    (hcol : (fun _ => (1 : ℝ)) ᵥ* P.L = 0)
    (h : (P.sum Q).Feasible (Sum.elim u v) (Sum.elim μ ν)) : ∑ i, μ i = ∑ i, P.Ω i :=
  ConeProblem.Feasible.sum_cone_eq_sum_defect hcol ((feasible_sum_iff P Q u v μ ν).mp h).1

/-- Gauss–Bonnet on the second component, if it is closed (`1ᵀ Q.L = 0`).

[formula: closed-mesh-cone-sum] -/
theorem Feasible.sum_cone_eq_sum_defect_inr {u : ι → ℝ} {v : κ → ℝ} {μ : ι → ℝ} {ν : κ → ℝ}
    (hcol : (fun _ => (1 : ℝ)) ᵥ* Q.L = 0)
    (h : (P.sum Q).Feasible (Sum.elim u v) (Sum.elim μ ν)) : ∑ j, ν j = ∑ j, Q.Ω j :=
  ConeProblem.Feasible.sum_cone_eq_sum_defect hcol ((feasible_sum_iff P Q u v μ ν).mp h).2

/-! ## A component solved alone: masses and `λ` -/

variable {P}

/-- Multiplying the masses by `c` and `λ` by `c` multiplies the objective by `c`. -/
theorem objective_scaleMass_withLam (c : ℝ) (u μ : ι → ℝ) :
    ((P.scaleMass c).withLam (c * P.lam)).objective u μ = c * P.objective u μ := by
  simp only [objective, distortion, coneCost, scaleMass, withLam, Finset.mul_sum, mul_add]
  congr 1 <;> exact Finset.sum_congr rfl fun i _ => by ring

/-- **A component solved with its own mass normalization.** For `c > 0`, multiplying all masses
by `c` and `λ` by `c` does not change the minimizers. With `c = A / A_c` (`A` the total area of
the mesh, `A_c` the area of the component) this says: the component solved alone, with masses
`a_i / A_c` and `λ_c = (A / A_c)·λ`, has the minimizers of the component's part of the whole
problem (masses `a_i / A`, parameter `λ`).

[formula: component-lambda-rescaling] -/
theorem isMinimizer_scaleMass_withLam_iff {c : ℝ} (hc : 0 < c) (u μ : ι → ℝ) :
    ((P.scaleMass c).withLam (c * P.lam)).IsMinimizer u μ ↔ P.IsMinimizer u μ := by
  have hfeas : ∀ v ν, ((P.scaleMass c).withLam (c * P.lam)).Feasible v ν ↔ P.Feasible v ν :=
    fun _ _ => Iff.rfl
  simp only [IsMinimizer, hfeas, objective_scaleMass_withLam]
  refine and_congr Iff.rfl (forall_congr' fun v => forall_congr' fun ν =>
    imp_congr_right fun _ => mul_le_mul_iff_of_pos_left hc)

/-- **The adjoint of a component solved alone.** For `c > 0`: `(u, φ, μ)` solves the optimality
system of `P` iff `(u, c·φ, μ)` solves that of the problem with masses `c·m` and parameter
`c·λ`. With `c = A / A_c`: `φ_alone = (A / A_c)·φ_global`, i.e. `φ_global = (A_c / A)·φ_alone`.

[formula: component-lambda-rescaling] -/
theorem kkt_scaleMass_withLam_iff {c : ℝ} (hc : 0 < c) (u φ μ : ι → ℝ) :
    ((P.scaleMass c).withLam (c * P.lam)).KKT u (c • φ) μ ↔ P.KKT u φ μ := by
  have hadj : ((P.scaleMass c).withLam (c * P.lam)).Lᵀ *ᵥ (c • φ) =
        (fun i => -(((P.scaleMass c).withLam (c * P.lam)).wE i *
          ((P.scaleMass c).withLam (c * P.lam)).m i * u i)) ↔
      P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i) := by
    simp only [scaleMass, withLam, Matrix.mulVec_smul]
    have e : (fun i => -(P.wE i * (c * P.m i) * u i)) = c • fun i => -(P.wE i * P.m i * u i) := by
      funext i; simp only [Pi.smul_apply, smul_eq_mul]; ring
    rw [e]
    exact (smul_right_injective _ hc.ne').eq_iff
  constructor
  · intro h
    refine ⟨h.state, hadj.mp h.adjoint, fun i => ?_⟩
    have := h.inclusion i
    simp only [scaleMass, withLam, Pi.smul_apply, smul_eq_mul] at this
    exact (Site.stationary_mul_iff _ hc _ _ _).mp this
  · intro h
    refine ⟨h.state, hadj.mpr h.adjoint, fun i => ?_⟩
    simp only [scaleMass, withLam, Pi.smul_apply, smul_eq_mul]
    exact (Site.stationary_mul_iff _ hc _ _ _).mpr (h.inclusion i)

/-- The same, solved for the adjoint of the rescaled problem: `φ` is an adjoint of the problem
with masses `c·m` and parameter `c·λ` iff `c⁻¹·φ` is an adjoint of `P`.

[formula: component-lambda-rescaling] -/
theorem kkt_scaleMass_withLam_iff' {c : ℝ} (hc : 0 < c) (u φ μ : ι → ℝ) :
    ((P.scaleMass c).withLam (c * P.lam)).KKT u φ μ ↔ P.KKT u (c⁻¹ • φ) μ := by
  rw [← kkt_scaleMass_withLam_iff (P := P) hc u (c⁻¹ • φ) μ, smul_smul, mul_inv_cancel₀ hc.ne', one_smul]

omit [Fintype ι] in
/-- The component normalized by its own area `A_c` with `λ_c = (A / A_c)·λ` is the component
normalized by the total area `A` with `λ`, its masses then multiplied by `A / A_c`. -/
theorem normalizeMass_alone_eq (R : ConeProblem ι) {A Ac : ℝ} (hA : A ≠ 0) (lam : ℝ) :
    (R.normalizeMass Ac).withLam (A / Ac * lam) =
      (((R.normalizeMass A).withLam lam).scaleMass (A / Ac)).withLam
        (A / Ac * ((R.normalizeMass A).withLam lam).lam) := by
  simp only [normalizeMass, withLam, scaleMass]
  congr 1
  funext i
  field_simp

/-- **Solving a component alone** (the form `solve.rs` uses). `R` is the component with its
raw lumped masses `a_i`; `A > 0` is the total area of the whole mesh and `A_c > 0` the area of the
component. The component normalized by its own area with `λ_c = (A / A_c)·λ` has the minimizers
of the component normalized by the total area with `λ` (its part of the whole problem).

[formula: component-lambda-rescaling] -/
theorem isMinimizer_normalizeMass_alone_iff (R : ConeProblem ι) {A Ac : ℝ} (hA : 0 < A)
    (hAc : 0 < Ac) (lam : ℝ) (u μ : ι → ℝ) :
    ((R.normalizeMass Ac).withLam (A / Ac * lam)).IsMinimizer u μ ↔
      ((R.normalizeMass A).withLam lam).IsMinimizer u μ := by
  rw [normalizeMass_alone_eq R hA.ne' lam]
  exact isMinimizer_scaleMass_withLam_iff (div_pos hA hAc) u μ

/-- **The adjoint of a component solved alone**: with the normalizations of
`isMinimizer_normalizeMass_alone_iff`, `(u, φ_alone, μ)` solves the optimality system of the
alone-solve iff `(u, (A_c / A)·φ_alone, μ)` solves that of the component's part of the whole
problem: `φ_global = (A_c / A)·φ_alone`.

[formula: component-lambda-rescaling] -/
theorem kkt_normalizeMass_alone_iff (R : ConeProblem ι) {A Ac : ℝ} (hA : 0 < A) (hAc : 0 < Ac)
    (lam : ℝ) (u φ μ : ι → ℝ) :
    ((R.normalizeMass Ac).withLam (A / Ac * lam)).KKT u φ μ ↔
      ((R.normalizeMass A).withLam lam).KKT u ((Ac / A) • φ) μ := by
  rw [normalizeMass_alone_eq R hA.ne' lam, kkt_scaleMass_withLam_iff' (div_pos hA hAc), inv_div]

variable (P)

/-- **The whole problem against the alone-solves.** For a common `λ` and any `a, b > 0`, a pair
minimizes the problem of the whole mesh iff each restriction minimizes its component's problem
with masses multiplied by its factor (`a = A / A_P`, `b = A / A_Q`: each component normalized by
its own area) and `λ` multiplied by the same factor.

[formula: component-separation] -/
theorem isMinimizer_sum_iff_alone (hlam : P.lam = Q.lam) {a b : ℝ} (ha : 0 < a) (hb : 0 < b)
    (u : ι → ℝ) (v : κ → ℝ) (μ : ι → ℝ) (ν : κ → ℝ) :
    (P.sum Q).IsMinimizer (Sum.elim u v) (Sum.elim μ ν) ↔
      ((P.scaleMass a).withLam (a * P.lam)).IsMinimizer u μ ∧
        ((Q.scaleMass b).withLam (b * Q.lam)).IsMinimizer v ν := by
  rw [isMinimizer_sum_iff P Q hlam, isMinimizer_scaleMass_withLam_iff ha,
    isMinimizer_scaleMass_withLam_iff hb]

end ConeProblem

/-! ## The hypotheses are satisfiable -/

open ConeProblem in
/-- Two one-vertex components with a common `λ`: the assembled triple solves the optimality
system of the whole problem, hence is a global minimizer of it. -/
example : ∃ (P Q : ConeProblem (Fin 1)) (u v φ ψ μ ν : Fin 1 → ℝ),
    P.lam = Q.lam ∧ (P.sum Q).WellFormed ∧
      (P.sum Q).KKT (Sum.elim u v) (Sum.elim φ ψ) (Sum.elim μ ν) ∧
      (P.sum Q).IsMinimizer (Sum.elim u v) (Sum.elim μ ν) := by
  let P : ConeProblem (Fin 1) := ⟨1, fun _ => 3, fun _ => 1, fun _ => 1, 1, fun _ => Site.plain⟩
  have hP : P.WellFormed := ⟨fun _ => zero_le_one, fun _ => zero_le_one, zero_le_one,
    fun _ => ⟨zero_le_one, bot_le, le_top⟩⟩
  have hK : P.KKT (fun _ => -1) (fun _ => 1) (fun _ => 2) := by
    refine ⟨?_, ?_, ?_⟩
    · ext i; simp [P]; norm_num
    · ext i; simp [P]
    · intro i
      rw [Site.stationary_plain_iff zero_le_one]
      simp [softThreshold]; norm_num
  have hS : (P.sum P).KKT (Sum.elim (fun _ => -1) (fun _ => -1)) (Sum.elim (fun _ => 1) fun _ => 1)
      (Sum.elim (fun _ => 2) fun _ => 2) := (kkt_sum_iff P P rfl _ _ _ _ _ _).mpr ⟨hK, hK⟩
  have hW := wellFormed_sum P P hP hP
  exact ⟨P, P, _, _, _, _, _, _, rfl, hW, hS,
    hS.isMinimizer hW.wE_nonneg hW.m_nonneg⟩

open ConeProblem in
/-- A closed two-vertex component (`1ᵀL = 0`) next to a one-vertex component: the Gauss–Bonnet
corollary applies to a feasible pair of the whole problem. -/
example : ∃ (P : ConeProblem (Fin 2)) (Q : ConeProblem (Fin 1)) (u μ : Fin 2 → ℝ)
    (v ν : Fin 1 → ℝ), (fun _ => (1 : ℝ)) ᵥ* P.L = 0 ∧
      (P.sum Q).Feasible (Sum.elim u v) (Sum.elim μ ν) ∧ ∑ i, μ i = ∑ i, P.Ω i := by
  let P : ConeProblem (Fin 2) :=
    ⟨!![1, -1; -1, 1], fun _ => 1, fun _ => 1, fun _ => 1, 1, fun _ => Site.plain⟩
  let Q : ConeProblem (Fin 1) := ⟨1, fun _ => 3, fun _ => 1, fun _ => 1, 1, fun _ => Site.plain⟩
  have hcol : (fun _ => (1 : ℝ)) ᵥ* P.L = 0 := by
    ext i; fin_cases i <;> simp [P, Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  have hf : (P.sum Q).Feasible (Sum.elim (fun _ => 0) (fun _ => -1))
      (Sum.elim (fun _ => 1) (fun _ => 2)) := by
    refine (feasible_sum_iff P Q _ _ _ _).mpr ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩⟩
    · ext i; fin_cases i <;> simp [P, Matrix.mulVec, dotProduct]
    · intro i; exact ⟨bot_le, le_top⟩
    · ext i; simp [Q]; norm_num
    · intro i; exact ⟨bot_le, le_top⟩
  exact ⟨P, Q, _, _, _, _, hcol, hf, Feasible.sum_cone_eq_sum_defect_inl P Q hcol hf⟩

open ConeProblem in
/-- The rescaling of a component: with `c = 2`, the adjoint of the rescaled problem is twice the
adjoint of the original one. -/
example : ∃ (P : ConeProblem (Fin 1)) (u φ μ : Fin 1 → ℝ),
    P.KKT u φ μ ∧ ((P.scaleMass 2).withLam (2 * P.lam)).KKT u ((2 : ℝ) • φ) μ := by
  let P : ConeProblem (Fin 1) := ⟨1, fun _ => 3, fun _ => 1, fun _ => 1, 1, fun _ => Site.plain⟩
  have hK : P.KKT (fun _ => -1) (fun _ => 1) (fun _ => 2) := by
    refine ⟨?_, ?_, ?_⟩
    · ext i; simp [P]; norm_num
    · ext i; simp [P]
    · intro i
      rw [Site.stationary_plain_iff zero_le_one]
      simp [softThreshold]; norm_num
  exact ⟨P, _, _, _, hK, (kkt_scaleMass_withLam_iff two_pos _ _ _).mpr hK⟩

end DdgProofs
