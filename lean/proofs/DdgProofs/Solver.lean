import DdgProofs.Laplace.Kernel
import DdgProofs.Certificate
import DdgProofs.Newton
import DdgProofs.Closed

/-!
# The linear systems the cone solver solves, its closed-mesh treatment, and cone merging

This file is about `crates/ddg-optimal-cones`: `linsys.rs` (the Newton systems), `newton.rs` (the
two Newton iterations, the shift of the adjoint constant, the line search), `evaluate.rs` (the
scale factor of a given cone configuration) and `extract.rs` (merging clusters of cones).

* **`saddle-point-newton-system`.** Every Newton step solves a system in `(u, φ)` described by one
  `Row` per vertex (`flat`, `sloped`, `pinned`). `isSolution_iff`: its solutions are exactly the
  pairs for which the adjoint equation holds and the induced cone angles `μ = A u + Ω` satisfy
  the relation of each row. `isAssembledSolution_iff`: the matrix the code assembles (unknowns
  `(u, ψ)` with `ψ_i = φ_i − offset_i`, known values moved to the right-hand side, pinned rows
  replaced by the identity) has the same solutions. `newton_step_iff_isSolution`: the
  semismooth Newton step of the residual map of `DDG.Newton` is the row system of the
  current pieces. `kkt_of_isSolution`, `isMinimizer_of_isSolution`:
  a solution all of whose vertices lie on their pieces is a solution of the optimality system,
  hence a global minimizer; `isProxMinimizer_of_isSolution`: the same for the regularized
  (proximal) stage; `Site.stationary_of_pos` … `stationary_of_lower`: the test "on its piece" for
  each kind of row.
* **Closed meshes.** `mulVec_eq_of_sum_eq`: with zero column sums, an equation of `A u = b` is
  implied by the others when `Σ b = 0` (the dropped state equation of the pinned vertex; the
  dropped row of `LaplaceSolver`). `IsSolution.sum_weighted_eq_zero`: the adjoint equations fix
  the additive constant of `u` (`Σ w u = 0`). `IsSolution.add_const`: without a sloped or pinned
  row, `φ` is determined only up to a constant.
* **`adjoint-constant-interval`.** `kkt_add_const_iff`, `Site.stationary_iff_mem_subdiff`: the
  constants `c` for which `φ + c` completes `(u, μ)` to a solution of the optimality system are
  those in the intersection of the intervals `∂g_i(μ_i) − φ_i`.
* **`gauss-bonnet-adjoint-shift`.** `sum_proxReg_add_monotone`: the total cone angle that
  `φ + c` asks for is nondecreasing in `c`; `isSolution_of_pinned_of_sum_eq`: once it equals `Σ Ω`
  the pinned system is the full one.
* **`regularized-dual-objective`.** `dualSlope_eq`: the quantity the line search evaluates is
  minus the pairing of the direction `dφ` with the state residual at the trial point;
  `dualSlope_monotone`: it is nondecreasing in the step length. (The dual objective `Φ` itself,
  an antiderivative, is not formalized.)
* **`proximal-fixed-point`.** `IsProxMinimizer.isMinimizer_of_fixedPoint`: a minimizer of the
  proximal problem whose centre is its own cone vector minimizes the objective itself;
  `IsProxMinimizer.objective_add_le`: a proximal step does not increase the objective.
* **`componentwise-backward-error`.** `oettli_prager`: the accepted-solve criterion of the linear
  solver is exactly "solves a system with entries perturbed by at most the relative amount `ε`".
* **`cone-cluster-frechet-mean`.** `sum_mergeCones`, `sum_abs_mergeCones_le`,
  `mergeCones_mem_feasible`, `exists_frechetMean`: merging a cluster into one cone at *any* of its
  vertices keeps the total angle, does not increase `Σ|μ|`, and is admissible iff the summed angle
  is admissible at the chosen vertex; a vertex of least Fréchet variance exists.
-/

noncomputable section

namespace DDG

open Matrix

/-! ## Monotonicity of the proximal maps -/

theorem softThreshold_monotone (t : ℝ) : Monotone (softThreshold t) := by
  intro x y h
  unfold softThreshold
  exact add_le_add (max_le_max le_rfl (by linarith)) (min_le_min le_rfl (by linarith))

/-- The clamp to `[lo, hi]` is nondecreasing (`lo ≤ 0 ≤ hi`). -/
theorem clampE_monotone {lo hi : EReal} (hlo : lo ≤ 0) (hhi : 0 ≤ hi) :
    Monotone (clampE lo hi) := by
  intro x y h
  by_cases h1 : lo ≤ (y : EReal)
  · by_cases h2 : (y : EReal) ≤ hi
    · rw [clampE_of_mem h1 h2]
      by_cases h3 : lo ≤ (x : EReal)
      · rw [clampE_of_mem h3 (le_trans (by exact_mod_cast h) h2)]; exact h
      · obtain ⟨r, hr, hc⟩ := clampE_of_le_lo hlo hhi (le_of_lt (not_le.mp h3))
        rw [hc]
        rw [hr] at h1
        exact_mod_cast h1
    · obtain ⟨r, hr, hc⟩ := clampE_of_hi_le hlo hhi (le_of_lt (not_le.mp h2))
      rw [hc]
      have hmem : ((clampE lo hi x : ℝ) : EReal) ≤ (r : EReal) :=
        le_of_le_of_eq (clampE_mem hlo hhi x).2 hr
      exact_mod_cast hmem
  · have h3 : (x : EReal) ≤ lo :=
      le_trans (by exact_mod_cast h) (le_of_lt (not_le.mp h1))
    obtain ⟨r, hr, hc⟩ := clampE_of_le_lo hlo hhi (le_of_lt (not_le.mp h1))
    obtain ⟨r', hr', hc'⟩ := clampE_of_le_lo hlo hhi h3
    rw [hc, hc']
    have : (r' : EReal) = r := by rw [← hr, ← hr']
    exact le_of_eq (by exact_mod_cast this)

/-- The (regularized) proximal map of every vertex kind is nondecreasing. -/
theorem Site.proxReg_monotone {lam γ : ℝ} (hγ : 0 < γ) {s : Site} (hs : s.WellFormed) :
    Monotone (s.proxReg lam γ) := by
  cases s with
  | penalized wR lo hi =>
    intro x y h
    simp only [Site.proxReg]
    exact clampE_monotone hs.2.1 hs.2.2
      (div_le_div_of_nonneg_right (softThreshold_monotone _ h) hγ.le)
  | fixed a => intro x y _; simp [Site.proxReg]

namespace Solver

variable {ι : Type*} [Fintype ι]

/-! ## The Newton systems (`linsys.rs`) -/

/-- How one vertex enters a Newton system (transcription of `linsys::Row`). -/
inductive Row where
  /-- The cone angle is the given number `mu`. -/
  | flat (mu : ℝ)
  /-- The cone angle is the affine function `base + k·(φ_i − offset)` of the adjoint. -/
  | sloped (k offset base : ℝ)
  /-- The adjoint is the given number `phi`; the cone angle is free. -/
  | pinned (phi : ℝ)

namespace Row

/-- The equation of the row of vertex `i` (module documentation of `linsys.rs`):
`(A u)_i = μ_i − Ω_i`, `(A u)_i − k (φ_i − o) = b − Ω_i`, or `φ_i = p`. -/
def Holds (A : Matrix ι ι ℝ) (Ω u φ : ι → ℝ) (i : ι) : Row → Prop
  | flat mu => (A *ᵥ u) i = mu - Ω i
  | sloped k o b => (A *ᵥ u) i - k * (φ i - o) = b - Ω i
  | pinned p => φ i = p

/-- The relation between the cone angle `μ_i` and the adjoint `φ_i` that the row stands for. -/
def Relation (μ φ : ℝ) : Row → Prop
  | flat mu => μ = mu
  | sloped k o b => μ = b + k * (φ - o)
  | pinned p => φ = p

/-- The value that is known in advance and moved to the right-hand side: the pinned adjoint
value, or the offset of a sloped row (whose unknown is `ψ_i = φ_i − offset`). -/
def known : Row → ℝ
  | flat _ => 0
  | sloped _ o _ => o
  | pinned p => p

/-- Whether the row is pinned. -/
def IsPinned : Row → Prop
  | pinned _ => True
  | _ => False

instance : DecidablePred IsPinned := fun r => by cases r <;> unfold IsPinned <;> infer_instance

/-- The adjoint value returned for a row, from the unknown `ψ_i` of the assembled system
(`SaddleSolver::solve`: a pinned entry is returned as given, the offset is added back on a
sloped row). -/
def recover (ψ : ℝ) : Row → ℝ
  | flat _ => ψ
  | sloped _ o _ => ψ + o
  | pinned p => p

end Row

/-- `(u, φ)` solves the Newton system of the rows `rows`: the adjoint equation
`W u + A φ = 0` at every vertex and the equation of every row. -/
structure IsSolution (A : Matrix ι ι ℝ) (w Ω : ι → ℝ) (rows : ι → Row) (u φ : ι → ℝ) : Prop where
  adjoint : ∀ i, w i * u i + (A *ᵥ φ) i = 0
  row : ∀ i, (rows i).Holds A Ω u φ i

/-- **What a solution of a Newton system is.** `(u, φ)` solves the system of `rows` iff the
adjoint equation `A φ = −w∘u` holds and the cone angles `μ = A u + Ω` that `u` induces (so the
state equation `A u = μ − Ω` holds by definition) satisfy the relation of every row: `μ_i` is the
given number (flat), `μ_i = b_i + k_i (φ_i − o_i)` (sloped), `φ_i = p_i` (pinned).

The code reads `μ_i` off as the given number on a flat row and as `(A u)_i + Ω_i` on the other
rows; for an exact solution these are the same vector `A u + Ω`.

[formula: saddle-point-newton-system] -/
theorem isSolution_iff (A : Matrix ι ι ℝ) (w Ω : ι → ℝ) (rows : ι → Row) (u φ : ι → ℝ) :
    IsSolution A w Ω rows u φ ↔
      (A *ᵥ φ = fun i => -(w i * u i)) ∧
        ∀ i, (rows i).Relation ((A *ᵥ u) i + Ω i) (φ i) := by
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨funext fun i => by have := h1 i; linarith, fun i => ?_⟩
    have := h2 i
    cases hr : rows i with
    | flat mu => rw [hr] at this; simp only [Row.Holds] at this; simp only [Row.Relation]; linarith
    | sloped k o b =>
      rw [hr] at this; simp only [Row.Holds] at this; simp only [Row.Relation]; linarith
    | pinned p => rw [hr] at this; exact this
  · rintro ⟨h1, h2⟩
    refine ⟨fun i => by
      have : (A *ᵥ φ) i = -(w i * u i) := congrFun h1 i
      linarith, fun i => ?_⟩
    have := h2 i
    cases hr : rows i with
    | flat mu => rw [hr] at this; simp only [Row.Relation] at this; simp only [Row.Holds]; linarith
    | sloped k o b =>
      rw [hr] at this; simp only [Row.Relation] at this; simp only [Row.Holds]; linarith
    | pinned p => rw [hr] at this; exact this

/-- `(u, ψ)` solves the system **as `SaddleSolver::assemble` writes it**: the unknown of a sloped
row is `ψ_i = φ_i − offset_i`, the columns of pinned rows are removed and their known values,
like the offsets, are moved to the right-hand side:

* row `j` of the first block: `w_j u_j + Σ_{i not pinned} A_ji ψ_i = −Σ_i A_ji·known_i`;
* row `i` of the second block: `(A u)_i = μ_i − Ω_i` (flat), `(A u)_i − k_i ψ_i = b_i − Ω_i`
  (sloped), `ψ_i = p_i` (pinned).

For symmetric `A` its matrix is the symmetric `[[W, A], [A, −K]]` with pinned rows and columns
replaced by the identity. -/
structure IsAssembledSolution (A : Matrix ι ι ℝ) (w Ω : ι → ℝ) (rows : ι → Row)
    (u ψ : ι → ℝ) : Prop where
  first : ∀ j, w j * u j + ∑ i, (if (rows i).IsPinned then 0 else A j i * ψ i) =
    -∑ i, A j i * (rows i).known
  second : ∀ i, match rows i with
    | Row.flat mu => (A *ᵥ u) i = mu - Ω i
    | Row.sloped k _ b => (A *ᵥ u) i - k * ψ i = b - Ω i
    | Row.pinned p => ψ i = p

private theorem adjoint_split (A : Matrix ι ι ℝ) (rows : ι → Row) (ψ : ι → ℝ) (j : ι) :
    (A *ᵥ fun i => (rows i).recover (ψ i)) j =
      ∑ i, (if (rows i).IsPinned then 0 else A j i * ψ i) + ∑ i, A j i * (rows i).known := by
  simp only [Matrix.mulVec, dotProduct]
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  cases rows i <;> simp [Row.recover, Row.IsPinned, Row.known, mul_add]

/-- **The assembled system is the Newton system.** `(u, ψ)` solves the system as assembled iff
`(u, φ)` with `φ_i = recover(ψ_i)` (the vector the solver returns) solves the Newton system of the
rows and `ψ` has the given values on the pinned rows. Conversely every solution `(u, φ)` of the
Newton system comes from the assembled solution `ψ_i = φ_i − known_i` on the rows that are not
pinned.

[formula: saddle-point-newton-system] -/
theorem isAssembledSolution_iff (A : Matrix ι ι ℝ) (w Ω : ι → ℝ) (rows : ι → Row)
    (u ψ : ι → ℝ) :
    IsAssembledSolution A w Ω rows u ψ ↔
      IsSolution A w Ω rows u (fun i => (rows i).recover (ψ i)) ∧
        ∀ i p, rows i = Row.pinned p → ψ i = p := by
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨⟨fun j => ?_, fun i => ?_⟩, fun i p hp => ?_⟩
    · rw [adjoint_split]; have := h1 j; linarith
    · have := h2 i
      cases hr : rows i with
      | flat mu => rw [hr] at this; simpa [Row.Holds] using this
      | sloped k o b =>
        rw [hr] at this
        simp only [Row.Holds, Row.recover, hr]
        simp only at this
        linarith
      | pinned p => simp [Row.Holds, Row.recover, hr]
    · have := h2 i; rw [hp] at this; exact this
  · rintro ⟨⟨h1, h2⟩, h3⟩
    refine ⟨fun j => ?_, fun i => ?_⟩
    · have := h1 j; rw [adjoint_split] at this; linarith
    · have := h2 i
      cases hr : rows i with
      | flat mu => rw [hr] at this; simpa [Row.Holds] using this
      | sloped k o b =>
        rw [hr] at this
        simp only [Row.Holds, Row.recover, hr] at this
        simp only
        linarith
      | pinned p => exact h3 i p hr

/-- Every solution of the Newton system is returned by the assembled one. -/
theorem IsSolution.isAssembledSolution {A : Matrix ι ι ℝ} {w Ω : ι → ℝ} {rows : ι → Row}
    {u φ : ι → ℝ} (h : IsSolution A w Ω rows u φ) :
    IsAssembledSolution A w Ω rows u (fun i => φ i - (rows i).known +
        (match rows i with | Row.pinned p => p | _ => 0)) ∧
      (fun i => (rows i).recover (φ i - (rows i).known +
        (match rows i with | Row.pinned p => p | _ => 0))) = φ := by
  have e : (fun i => (rows i).recover (φ i - (rows i).known +
      (match rows i with | Row.pinned p => p | _ => 0))) = φ := by
    funext i
    have := h.row i
    cases hr : rows i with
    | flat mu => simp [Row.recover, Row.known]
    | sloped k o b => simp [Row.recover, Row.known]
    | pinned p => rw [hr] at this; simp [Row.recover]; exact this.symm
  refine ⟨(isAssembledSolution_iff A w Ω rows u _).mpr ⟨by rw [e]; exact h, fun i p hp => ?_⟩, e⟩
  have := h.row i
  rw [hp] at this
  simp only [hp, Row.known, Row.Holds] at this ⊢
  linarith

/-! ## Closed meshes: the constants of `u` and `φ`, and the dropped equation -/

/-- **The adjoint equations fix the additive constant of `u`** on a closed mesh: if `1ᵀA = 0`,
every solution of a Newton system has `Σ_i w_i u_i = 0`.

[formula: closed-mesh-optimal-constant] -/
theorem IsSolution.sum_weighted_eq_zero {A : Matrix ι ι ℝ} {w Ω : ι → ℝ} {rows : ι → Row}
    {u φ : ι → ℝ} (hcol : (fun _ => (1 : ℝ)) ᵥ* A = 0) (h : IsSolution A w Ω rows u φ) :
    ∑ i, w i * u i = 0 := by
  have h1 : ∑ i, (A *ᵥ φ) i = 0 := by
    have : (fun _ => (1 : ℝ)) ⬝ᵥ (A *ᵥ φ) = 0 := by
      rw [Matrix.dotProduct_mulVec, hcol]; simp
    simpa [dotProduct] using this
  have h2 : ∑ i, (w i * u i + (A *ᵥ φ) i) = 0 := Finset.sum_eq_zero fun i _ => h.adjoint i
  rw [Finset.sum_add_distrib, h1] at h2
  linarith

/-- **Without an active vertex the adjoint is determined only up to a constant** on a closed
mesh: if `A·1 = 0` and every row is flat, `φ ↦ φ + c` maps solutions to solutions. -/
theorem IsSolution.add_const {A : Matrix ι ι ℝ} {w Ω : ι → ℝ} {rows : ι → Row} {u φ : ι → ℝ}
    (hrow : A *ᵥ (fun _ => (1 : ℝ)) = 0) (hflat : ∀ i, ∃ mu, rows i = Row.flat mu)
    (h : IsSolution A w Ω rows u φ) (c : ℝ) :
    IsSolution A w Ω rows u (fun i => φ i + c) := by
  have e : (A *ᵥ fun i => φ i + c) = A *ᵥ φ := by
    have : (fun i => φ i + c) = φ + c • (fun _ => (1 : ℝ)) := by ext i; simp
    rw [this, Matrix.mulVec_add, Matrix.mulVec_smul, hrow]; simp
  refine ⟨fun i => by rw [e]; exact h.adjoint i, fun i => ?_⟩
  obtain ⟨mu, hmu⟩ := hflat i
  have := h.row i
  rw [hmu] at this ⊢
  exact this

/-- **The pinned system of a closed mesh is the full one once Gauss–Bonnet holds.** Let
`1ᵀA = 0`. Let every row but that of `i₀` be flat with cone angle `mu i`, let the row of `i₀` be
pinned (its state equation dropped), and let `Σ_i mu_i = Σ_i Ω_i`. Then a solution of the pinned
system satisfies the dropped state equation `(A u)_{i₀} = mu_{i₀} − Ω_{i₀}` as well, so it solves
the system in which every row is flat.

[formula: gauss-bonnet-adjoint-shift] -/
theorem isSolution_of_pinned_of_sum_eq [DecidableEq ι] {A : Matrix ι ι ℝ} {w Ω mu : ι → ℝ}
    {u φ : ι → ℝ} (hcol : (fun _ => (1 : ℝ)) ᵥ* A = 0) (i₀ : ι) (p : ℝ)
    (hsum : ∑ i, mu i = ∑ i, Ω i)
    (h : IsSolution A w Ω (Function.update (fun i => Row.flat (mu i)) i₀ (Row.pinned p)) u φ) :
    IsSolution A w Ω (fun i => Row.flat (mu i)) u φ := by
  refine ⟨h.adjoint, fun i => ?_⟩
  simp only [Row.Holds]
  have hne : ∀ j, j ≠ i₀ → (A *ᵥ u) j = mu j - Ω j := by
    intro j hj
    have := h.row j
    rw [Function.update_of_ne hj] at this
    exact this
  by_cases hi : i = i₀
  · subst hi
    exact mulVec_eq_of_sum_eq hcol (b := fun j => mu j - Ω j)
      (by rw [Finset.sum_sub_distrib, hsum]; ring) i hne
  · exact hne i hi

end Solver

/-! ## From a solution of the solver's system to a minimizer -/

namespace ConeProblem

variable {ι : Type*} [Fintype ι] {P : ConeProblem ι} {u φ μ : ι → ℝ}

open Solver

/-- **A solution of a Newton system whose vertices satisfy the inclusion solves the optimality
system.** Let `L` be symmetric and `(u, φ)` a solution of the Newton system of any rows. Put
`μ = L u + Ω`. If `φ_i ∈ ∂g_i(μ_i)` at every vertex, then `(u, φ, μ)` solves the optimality system.

The rows enter only through the lemmas below, which give the inclusion from the relation of a
row together with the test "the vertex is still on the piece its row was built from".

[formula: saddle-point-newton-system] -/
theorem kkt_of_isSolution (hsymm : P.Lᵀ = P.L) {rows : ι → Row}
    (h : IsSolution P.L (fun i => P.wE i * P.m i) P.Ω rows u φ)
    (hincl : ∀ i, (P.site i).Stationary P.lam (φ i) (P.inducedCone u i)) :
    P.KKT u φ (P.inducedCone u) := by
  rw [isSolution_iff] at h
  refine ⟨?_, by rw [hsymm]; exact h.1, hincl⟩
  ext i; simp [inducedCone]

/-- **… and is a global minimizer.** This is the stopping test of the exact stage
(`solve_exact`): the step's system was solved and every vertex lies on its piece.

[formula: newton-stopping-rule] -/
theorem isMinimizer_of_isSolution (hP : P.WellFormed) (hsymm : P.Lᵀ = P.L) {rows : ι → Row}
    (h : IsSolution P.L (fun i => P.wE i * P.m i) P.Ω rows u φ)
    (hincl : ∀ i, (P.site i).Stationary P.lam (φ i) (P.inducedCone u i)) :
    P.IsMinimizer u (P.inducedCone u) :=
  (kkt_of_isSolution hsymm h hincl).isMinimizer hP.wE_nonneg hP.m_nonneg

/-! ### The Newton step of the residual map is the solver's linear system -/

/-- The slope of the proximal map on a piece is `0` or `1` (for `γ = 1`). -/
theorem slopeReg_one_eq_zero_or_one (lam : ℝ) (s : Site) (x : ℝ) :
    s.slopeReg lam 1 x = 0 ∨ s.slopeReg lam 1 x = 1 := by
  unfold Site.slopeReg
  cases s.pieceReg lam 1 x <;> simp [Site.pieceSlope]

open Classical in
/-- The row of a vertex in the exact stage, from the current point `x = φ_i + μ_i`: where the
proximal map has slope `0` the cone angle is fixed to its current value `prox(x)` (`0`, a bound,
or the fixed angle); where it has slope `1` the adjoint is pinned to `x − prox(x)` (which is
`λ wR` on the positive and `−λ wR` on the negative active piece). -/
def newtonRow (lam : ℝ) (s : Site) (x : ℝ) : Row :=
  if s.slopeReg lam 1 x = 0 then Row.flat (s.prox lam x) else Row.pinned (x - s.prox lam x)

/-- **The semismooth Newton step of Stage II is the solver's linear system.** For symmetric `L`:
`(v, q, ν)` solves the Newton system `DF(u, φ, μ)(v, q, ν) = −F(u, φ, μ)` of the residual map of
`DDG.Newton` (paper Eq. 26) if and only if the new iterate `(u + v, φ + q)` solves the
solver's row system for the rows `newtonRow` of the current point, and `μ + ν` is the cone vector
`L(u + v) + Ω` that the new `u` induces.

Together with `ConeProblem.isMinimizer_of_newton_step` this is the chain: solve the row system;
if every vertex is still on its piece, the result is a global minimizer.

[formula: saddle-point-newton-system] -/
theorem newton_step_iff_isSolution (hsymm : P.Lᵀ = P.L) (u φ μ v q ν : ι → ℝ) :
    P.residualDeriv u φ μ (v, q, ν) = -P.residual u φ μ ↔
      IsSolution P.L (fun i => P.wE i * P.m i) P.Ω
          (fun i => newtonRow P.lam (P.site i) (φ i + μ i)) (u + v) (φ + q) ∧
        μ + ν = P.inducedCone (u + v) := by
  simp only [residualDeriv, residual, Prod.neg_mk, Prod.mk.injEq, hsymm]
  have hLv : ∀ i, (P.L *ᵥ (u + v)) i = (P.L *ᵥ u) i + (P.L *ᵥ v) i := fun i => by
    rw [Matrix.mulVec_add]; rfl
  have hLq : ∀ i, (P.L *ᵥ (φ + q)) i = (P.L *ᵥ φ) i + (P.L *ᵥ q) i := fun i => by
    rw [Matrix.mulVec_add]; rfl
  constructor
  · rintro ⟨h1, h2, h3⟩
    have hμ : ∀ i, μ i + ν i = (P.L *ᵥ (u + v)) i + P.Ω i := by
      intro i
      have := congrFun h1 i
      simp only [Pi.sub_apply, Pi.neg_apply, Pi.add_apply] at this
      rw [hLv]; linarith
    refine ⟨⟨fun i => ?_, fun i => ?_⟩, funext fun i => by simpa [inducedCone] using hμ i⟩
    · have := congrFun h2 i
      simp only [Pi.add_apply, Pi.neg_apply] at this
      rw [hLq]; simp only [Pi.add_apply]; linarith
    · have h3i := congrFun h3 i
      simp only [Pi.neg_apply] at h3i
      unfold newtonRow
      rcases slopeReg_one_eq_zero_or_one P.lam (P.site i) (φ i + μ i) with hs | hs
      · rw [ite_eq_left hs]
        rw [hs] at h3i
        simp only [Row.Holds]
        have := hμ i
        linarith
      · rw [ite_eq_right (by rw [hs]; norm_num)]
        rw [hs] at h3i
        simp only [Row.Holds, Pi.add_apply]
        linarith
  · rintro ⟨⟨h1, h2⟩, h3⟩
    have hμ : ∀ i, μ i + ν i = (P.L *ᵥ (u + v)) i + P.Ω i := fun i => by
      have := congrFun h3 i; simpa [inducedCone] using this
    refine ⟨funext fun i => ?_, funext fun i => ?_, funext fun i => ?_⟩
    · simp only [Pi.sub_apply, Pi.neg_apply, Pi.add_apply]
      have := hμ i; rw [hLv] at this; linarith
    · simp only [Pi.add_apply, Pi.neg_apply]
      have := h1 i; rw [hLq] at this; simp only [Pi.add_apply] at this; linarith
    · simp only [Pi.neg_apply]
      have hrow := h2 i
      unfold newtonRow at hrow
      rcases slopeReg_one_eq_zero_or_one P.lam (P.site i) (φ i + μ i) with hs | hs
      · rw [ite_eq_left hs] at hrow
        rw [hs]
        simp only [Row.Holds] at hrow
        have := hμ i
        linarith
      · rw [ite_eq_right (by rw [hs]; norm_num)] at hrow
        rw [hs]
        simp only [Row.Holds, Pi.add_apply] at hrow
        linarith

/-! ### The regularized (proximal) stage -/

/-- The objective of the regularized problem the solver's Stage I and its safeguard solve
(`newton.rs`, `solve_regularized`): `J(u, μ) + Σ_i (μ_i − μ̄_i)²/(2 k_i)` with stiffness
`k_i > 0` and centre `μ̄` (`μ̄ = 0`, `k_i = m_i/γ` is the regularization of the paper's Stage I). -/
def proxObjective (P : ConeProblem ι) (k c : ι → ℝ) (u μ : ι → ℝ) : ℝ :=
  P.objective u μ + ∑ i, (μ i - c i) ^ 2 / (2 * k i)

/-- **Sufficiency for the proximal problem.** For a well-formed problem and `k_i > 0`: if
`L u = μ − Ω`, `Lᵀφ = −wE m u` and `φ_i − (μ_i − μ̄_i)/k_i ∈ ∂g_i(μ_i)` at every vertex, then
`(u, μ)` is feasible and minimizes `J(u, μ) + Σ (μ_i − μ̄_i)²/(2 k_i)` over all feasible pairs.

[formula: regularized-optimality-sufficient] -/
theorem isProxMinimizer_of_stationary (hP : P.WellFormed) {k c : ι → ℝ} (hk : ∀ i, 0 < k i)
    (hstate : P.L *ᵥ u = μ - P.Ω) (hadj : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i))
    (hst : ∀ i, (P.site i).Stationary P.lam (φ i - (μ i - c i) / k i) (μ i)) :
    P.Feasible u μ ∧ ∀ v ν, P.Feasible v ν → P.proxObjective k c u μ ≤ P.proxObjective k c v ν := by
  have hincl : ∀ i, IsSubgradientWithinAt
      (fun x => P.lam * (P.site i).cost x + (x - c i) ^ 2 / (2 * k i)) (P.site i).feasible
      (φ i) (μ i) := by
    intro i
    refine ⟨(hst i).1, fun y hy => ?_⟩
    have h1 := (hst i).2 y hy
    simp only at h1 ⊢
    have hk' := hk i
    have h2 : (y - c i) ^ 2 / (2 * k i) - (μ i - c i) ^ 2 / (2 * k i) -
        (μ i - c i) / k i * (y - μ i) = (y - μ i) ^ 2 / (2 * k i) := by
      field_simp; ring
    have h3 : 0 ≤ (y - μ i) ^ 2 / (2 * k i) := by positivity
    linarith
  refine ⟨⟨hstate, fun i => (hst i).1⟩, fun v ν hv => ?_⟩
  have := quadratic_growth_of_stationary P.L (fun i => P.wE i * P.m i)
    (fun i x => P.lam * (P.site i).cost x + (x - c i) ^ 2 / (2 * k i))
    (fun i => (P.site i).feasible) P.Ω u φ μ v ν hstate hadj hincl hv.1 hv.2
  have hnn : 0 ≤ ∑ i, 1 / 2 * (P.wE i * P.m i * (v i - u i) ^ 2) :=
    Finset.sum_nonneg fun i _ => mul_nonneg (by norm_num)
      (mul_nonneg (mul_nonneg (hP.wE_nonneg i) (hP.m_nonneg i)) (sq_nonneg _))
  simp only [proxObjective, objective, distortion, coneCost, Finset.mul_sum,
    Finset.sum_add_distrib] at this ⊢
  linarith

/-- The proximal inclusion in the form the code states it: for `k > 0`,
`φ − (μ − μ̄)/k ∈ ∂g(μ)  ⟺  μ = proxReg_{1/k}(φ + μ̄/k)`; for a penalized vertex the right side
is `clamp(k·P_t(φ + μ̄/k), lo, hi)`, the code's `Site::measure(k, φ + μ̄/k)`.

[formula: regularized-cone-update] -/
theorem stationary_prox_iff {lam k : ℝ} (hlam : 0 ≤ lam) (hk : 0 < k) {s : Site}
    (hs : s.WellFormed) (φ μ c : ℝ) :
    s.Stationary lam (φ - (μ - c) / k) μ ↔ μ = s.proxReg lam (1 / k) (φ + c / k) := by
  rw [← Site.stationary_sub_iff_eq_proxReg hlam (one_div_pos.mpr hk) hs]
  have e : φ + c / k - 1 / k * μ = φ - (μ - c) / k := by
    have := hk.ne'; field_simp; ring
  rw [e]

/-- **A fixed point of the proximal iteration solves the exact problem**: if the state and
adjoint equations hold and `μ` satisfies the proximal inclusion with centre `μ̄ = μ` itself, then
`(u, φ, μ)` solves the optimality system. -/
theorem kkt_of_prox_fixedPoint {k : ι → ℝ}
    (hstate : P.L *ᵥ u = μ - P.Ω) (hadj : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i))
    (hst : ∀ i, (P.site i).Stationary P.lam (φ i - (μ i - μ i) / k i) (μ i)) : P.KKT u φ μ := by
  refine ⟨hstate, hadj, fun i => ?_⟩
  have := hst i
  simpa using this

/-- `(u, μ)` is a global minimizer of the proximal problem with stiffness `k` and centre `c`:
feasible, and no feasible pair has a smaller `J + Σ (μ_i − c_i)²/(2 k_i)`. -/
def IsProxMinimizer (P : ConeProblem ι) (k c : ι → ℝ) (u μ : ι → ℝ) : Prop :=
  P.Feasible u μ ∧ ∀ v ν, P.Feasible v ν → P.proxObjective k c u μ ≤ P.proxObjective k c v ν

/-- The objective is convex along the segment between two pairs. -/
private theorem objective_convex_le (hP : P.WellFormed) (u μ v ν : ι → ℝ) {θ : ℝ} (h0 : 0 ≤ θ)
    (h1 : θ ≤ 1) :
    P.objective (fun i => (1 - θ) * u i + θ * v i) (fun i => (1 - θ) * μ i + θ * ν i) ≤
      (1 - θ) * P.objective u μ + θ * P.objective v ν := by
  have hc : P.coneCost (fun i => (1 - θ) * μ i + θ * ν i) ≤
      (1 - θ) * P.coneCost μ + θ * P.coneCost ν := by
    unfold coneCost
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    exact Finset.sum_le_sum fun i _ => (P.site i).cost_convex (hP.site i) (μ i) (ν i) h0 h1
  have hd : P.distortion (fun i => (1 - θ) * u i + θ * v i) ≤
      (1 - θ) * P.distortion u + θ * P.distortion v := by
    unfold distortion
    have : ∑ i, P.wE i * P.m i * ((1 - θ) * u i + θ * v i) ^ 2 ≤
        (1 - θ) * ∑ i, P.wE i * P.m i * u i ^ 2 + θ * ∑ i, P.wE i * P.m i * v i ^ 2 := by
      rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
      refine Finset.sum_le_sum fun i _ => ?_
      have hw : 0 ≤ P.wE i * P.m i := mul_nonneg (hP.wE_nonneg i) (hP.m_nonneg i)
      have key : (1 - θ) * u i ^ 2 + θ * v i ^ 2 - ((1 - θ) * u i + θ * v i) ^ 2 =
          θ * (1 - θ) * (u i - v i) ^ 2 := by ring
      have hnn : 0 ≤ θ * (1 - θ) * (u i - v i) ^ 2 :=
        mul_nonneg (mul_nonneg h0 (by linarith)) (sq_nonneg _)
      nlinarith [mul_nonneg hw hnn]
    linarith
  have := mul_le_mul_of_nonneg_left hc hP.lam_nonneg
  unfold objective
  nlinarith

/-- **A fixed point of the proximal iteration is a minimizer of the objective itself.** For a
well-formed problem and `k_i > 0`: if `(u, μ)` minimizes `J + Σ (ν_i − μ_i)²/(2 k_i)` — the
proximal problem whose centre is `μ` itself — then `(u, μ)` is a global minimizer of `J`.
(No optimality system is involved: the penalty and its first-order term vanish at the centre.)

[formula: proximal-fixed-point] -/
theorem IsProxMinimizer.isMinimizer_of_fixedPoint (hP : P.WellFormed) {k : ι → ℝ}
    (hk : ∀ i, 0 < k i) (h : P.IsProxMinimizer k μ u μ) : P.IsMinimizer u μ := by
  refine ⟨h.1, fun v ν hv => ?_⟩
  set C : ℝ := ∑ i, (ν i - μ i) ^ 2 / (2 * k i) with hC
  have hCnn : 0 ≤ C := Finset.sum_nonneg fun i _ => by
    have := hk i; positivity
  -- along the segment from (u, μ) to (v, ν)
  have key : ∀ θ : ℝ, 0 < θ → θ ≤ 1 →
      0 ≤ θ * (P.objective v ν - P.objective u μ) + θ ^ 2 * C := by
    intro θ h0 h1
    have hfeas : P.Feasible (fun i => (1 - θ) * u i + θ * v i)
        (fun i => (1 - θ) * μ i + θ * ν i) := by
      refine ⟨?_, fun i => (P.site i).feasible_convex (h.1.2 i) (hv.2 i) h0.le h1⟩
      have e : (fun i => (1 - θ) * u i + θ * v i) = (1 - θ) • u + θ • v := by ext i; simp
      rw [e, Matrix.mulVec_add, Matrix.mulVec_smul, Matrix.mulVec_smul, h.1.1, hv.1]
      ext i; simp; ring
    have h2 := h.2 _ _ hfeas
    have h3 := objective_convex_le hP u μ v ν h0.le h1
    have e1 : P.proxObjective k μ u μ = P.objective u μ := by
      simp [proxObjective]
    have e2 : P.proxObjective k μ (fun i => (1 - θ) * u i + θ * v i)
        (fun i => (1 - θ) * μ i + θ * ν i) =
        P.objective (fun i => (1 - θ) * u i + θ * v i) (fun i => (1 - θ) * μ i + θ * ν i) +
          θ ^ 2 * C := by
      simp only [proxObjective, hC, Finset.mul_sum]
      congr 1
      refine Finset.sum_congr rfl fun i _ => ?_
      have := (hk i).ne'
      field_simp; ring
    rw [e1, e2] at h2
    nlinarith
  have := nonneg_of_forall_small key
  linarith

/-- **A proximal step does not increase the objective**: if `(u, μ)` minimizes the proximal
problem with centre `c` and `(u₀, c)` is feasible (the previous iterate), then
`J(u, μ) + Σ (μ_i − c_i)²/(2 k_i) ≤ J(u₀, c)`; in particular `J(u, μ) ≤ J(u₀, c)` for `k_i > 0`.

[formula: proximal-fixed-point] -/
theorem IsProxMinimizer.objective_add_le {k c : ι → ℝ} {u₀ : ι → ℝ}
    (h : P.IsProxMinimizer k c u μ) (h0 : P.Feasible u₀ c) :
    P.objective u μ + ∑ i, (μ i - c i) ^ 2 / (2 * k i) ≤ P.objective u₀ c := by
  have := h.2 _ _ h0
  simpa [proxObjective] using this

/-- **A solution of the regularized stage's linear system, with every vertex on its piece,
solves the proximal problem.** Let `L` be symmetric, `k_i > 0`, and `(u, φ)` a solution of the
Newton system of any rows. Put `μ = L u + Ω`. If `φ_i − (μ_i − μ̄_i)/k_i ∈ ∂g_i(μ_i)` at every
vertex, then `(u, μ)` is feasible and minimizes `J + Σ (μ_i − μ̄_i)²/(2 k_i)`.

[formula: saddle-point-newton-system] -/
theorem isProxMinimizer_of_isSolution (hP : P.WellFormed) (hsymm : P.Lᵀ = P.L) {k c : ι → ℝ}
    (hk : ∀ i, 0 < k i) {rows : ι → Row}
    (h : IsSolution P.L (fun i => P.wE i * P.m i) P.Ω rows u φ)
    (hst : ∀ i, (P.site i).Stationary P.lam
      (φ i - (P.inducedCone u i - c i) / k i) (P.inducedCone u i)) :
    P.Feasible u (P.inducedCone u) ∧ ∀ v ν, P.Feasible v ν →
      P.proxObjective k c u (P.inducedCone u) ≤ P.proxObjective k c v ν := by
  rw [isSolution_iff] at h
  refine isProxMinimizer_of_stationary hP hk (φ := φ) ?_ (by rw [hsymm]; exact h.1) hst
  ext i; simp [inducedCone]

end ConeProblem

/-! ### When a vertex lies on the piece of its row

The rows the code builds, and the test under which the vertex "has not left its piece"
(`Site::next_piece`, pieces closed at their kinks). Each lemma gives the inclusion `g ∈ ∂g_i(μ)`
that `kkt_of_isSolution` (with `g = φ_i`) and `isProxMinimizer_of_isSolution` (with
`g = φ_i − (μ_i − μ̄_i)/k_i`) ask for. `t = λ·wR`.

| piece | exact stage | regularized stage | relation gives | test |
|---|---|---|---|---|
| positive cone | `pinned t` | `sloped k t μ̄` | `g = t` | `0 ≤ μ ≤ hi` |
| negative cone | `pinned (−t)` | `sloped k (−t) μ̄` | `g = −t` | `lo ≤ μ ≤ 0` |
| no cone | `flat 0` | `flat 0` | `μ = 0` | `g ≤ t` unless `hi = 0`, `−t ≤ g` unless `lo = 0` |
| upper bound | `flat hi` | `flat hi` | `μ = hi` | `t ≤ g` (or `−t ≤ g` if `hi = 0`) |
| lower bound | `flat lo` | `flat lo` | `μ = lo` | `g ≤ −t` (or `g ≤ t` if `lo = 0`) |
| fixed | `flat a` | `flat a` | `μ = a` | none |

For a sloped row `μ = μ̄ + k (φ − o)` the relation gives `φ − (μ − μ̄)/k = o`
(`Solver.Row.sloped_relation_iff`). -/

namespace Solver.Row

/-- The relation of a sloped row, solved for the argument of the inclusion:
`μ = b + k (φ − o)  ⟺  φ − (μ − b)/k = o` (`k ≠ 0`). -/
theorem sloped_relation_iff {k : ℝ} (hk : k ≠ 0) (o b μ φ : ℝ) :
    (Row.sloped k o b).Relation μ φ ↔ φ - (μ - b) / k = o := by
  simp only [Relation]
  constructor
  · intro h; rw [h]; field_simp; ring
  · intro h; rw [← h]; field_simp; ring

end Solver.Row

namespace Site

private theorem stationary_penalized_iff {lam wR : ℝ} {lo hi : EReal} (hlam : 0 ≤ lam)
    (hs : (penalized wR lo hi).WellFormed) (g μ : ℝ) :
    (penalized wR lo hi).Stationary lam g μ ↔
      μ ∈ boundSet lo hi ∧
      ((μ : EReal) < hi → g ≤ if 0 ≤ μ then lam * wR else -(lam * wR)) ∧
      (lo < (μ : EReal) → (if 0 < μ then lam * wR else -(lam * wR)) ≤ g) := by
  have hf : (fun x => lam * (penalized wR lo hi).cost x) = fun x => lam * wR * |x| := by
    funext x; simp [cost, mul_assoc]
  unfold Stationary
  rw [hf, feasible, isSubgradientWithinAt_abs_iff (mul_nonneg hlam hs.1) hs.2.1 hs.2.2]

/-- **Positive cone**: `g = λ wR` and `0 ≤ μ ≤ hi` give the inclusion.

[formula: soft-threshold-active-sets] -/
theorem stationary_of_pos {lam wR : ℝ} {lo hi : EReal} (hlam : 0 ≤ lam)
    (hs : (penalized wR lo hi).WellFormed) {μ : ℝ} (h0 : 0 ≤ μ) (hhi : (μ : EReal) ≤ hi) :
    (penalized wR lo hi).Stationary lam (lam * wR) μ := by
  have ht : 0 ≤ lam * wR := mul_nonneg hlam hs.1
  rw [stationary_penalized_iff hlam hs]
  refine ⟨⟨le_trans hs.2.1 (by exact_mod_cast h0), hhi⟩, fun _ => ?_, fun _ => ?_⟩
  · rw [ite_eq_left h0]
  · split_ifs <;> linarith

/-- **Negative cone**: `g = −λ wR` and `lo ≤ μ ≤ 0` give the inclusion.

[formula: soft-threshold-active-sets] -/
theorem stationary_of_neg {lam wR : ℝ} {lo hi : EReal} (hlam : 0 ≤ lam)
    (hs : (penalized wR lo hi).WellFormed) {μ : ℝ} (h0 : μ ≤ 0) (hlo : lo ≤ (μ : EReal)) :
    (penalized wR lo hi).Stationary lam (-(lam * wR)) μ := by
  have ht : 0 ≤ lam * wR := mul_nonneg hlam hs.1
  rw [stationary_penalized_iff hlam hs]
  refine ⟨⟨hlo, le_trans (by exact_mod_cast h0) hs.2.2⟩, fun _ => ?_, fun _ => ?_⟩
  · split_ifs <;> linarith
  · rw [ite_eq_right (not_lt.mpr h0)]

/-- **No cone**: `μ = 0` with `g ≤ λ wR` (not needed if `hi = 0`) and `−λ wR ≤ g` (not needed if
`lo = 0`) gives the inclusion.

[formula: soft-threshold-active-sets] -/
theorem stationary_of_zero {lam wR : ℝ} {lo hi : EReal} (hlam : 0 ≤ lam)
    (hs : (penalized wR lo hi).WellFormed) {g : ℝ} (h1 : 0 < hi → g ≤ lam * wR)
    (h2 : lo < 0 → -(lam * wR) ≤ g) : (penalized wR lo hi).Stationary lam g 0 := by
  rw [stationary_penalized_iff hlam hs]
  refine ⟨⟨by simpa using hs.2.1, by simpa using hs.2.2⟩, fun h => ?_, fun h => ?_⟩
  · rw [ite_eq_left le_rfl]; exact h1 (by simpa using h)
  · rw [ite_eq_right (lt_irrefl _)]; exact h2 (by simpa using h)

/-- **Upper bound**: `μ = hi` (finite) with `λ wR ≤ g` gives the inclusion; if `hi = 0`,
`−λ wR ≤ g` is enough.

[formula: soft-threshold-active-sets] -/
theorem stationary_of_upper {lam wR : ℝ} {lo : EReal} {h : ℝ} (hlam : 0 ≤ lam)
    (hs : (penalized wR lo (h : EReal)).WellFormed) {g : ℝ}
    (hg : (if 0 < h then lam * wR else -(lam * wR)) ≤ g) :
    (penalized wR lo (h : EReal)).Stationary lam g h := by
  rw [stationary_penalized_iff hlam hs]
  exact ⟨⟨le_trans hs.2.1 hs.2.2, le_rfl⟩, fun c => absurd c (lt_irrefl _), fun _ => hg⟩

/-- **Lower bound**: `μ = lo` (finite) with `g ≤ −λ wR` gives the inclusion; if `lo = 0`,
`g ≤ λ wR` is enough.

[formula: soft-threshold-active-sets] -/
theorem stationary_of_lower {lam wR : ℝ} {hi : EReal} {l : ℝ} (hlam : 0 ≤ lam)
    (hs : (penalized wR (l : EReal) hi).WellFormed) {g : ℝ}
    (hg : g ≤ if 0 ≤ l then lam * wR else -(lam * wR)) :
    (penalized wR (l : EReal) hi).Stationary lam g l := by
  rw [stationary_penalized_iff hlam hs]
  exact ⟨⟨le_rfl, le_trans hs.2.1 hs.2.2⟩, fun _ => hg, fun c => absurd c (lt_irrefl _)⟩

end Site

/-! ## The constant of the adjoint on a closed mesh (`adjoint_constant_interval`) -/

namespace Site

open Classical in
/-- Lower end of the subdifferential `∂(λ·cost + ι_feasible)(μ)` as an extended real
(`Site::subdifferential`): `−∞` at the lower bound and at a fixed vertex, otherwise `λ wR` if
`μ > 0` and `−λ wR` if not. -/
def subdiffLo (lam : ℝ) : Site → ℝ → EReal
  | penalized wR lo _, μ =>
    if lo < (μ : EReal) then ((if 0 < μ then lam * wR else -(lam * wR) : ℝ) : EReal) else ⊥
  | fixed _, _ => ⊥

open Classical in
/-- Upper end of the subdifferential: `+∞` at the upper bound and at a fixed vertex, otherwise
`λ wR` if `μ ≥ 0` and `−λ wR` if not. -/
def subdiffHi (lam : ℝ) : Site → ℝ → EReal
  | penalized wR _ hi, μ =>
    if (μ : EReal) < hi then ((if 0 ≤ μ then lam * wR else -(lam * wR) : ℝ) : EReal) else ⊤
  | fixed _, _ => ⊤

/-- **The per-vertex inclusion as an interval**: for an admissible `μ`,
`g ∈ ∂(λ·cost + ι_feasible)(μ)` iff `subdiffLo ≤ g ≤ subdiffHi`.

[formula: abs-subdifferential-bounded] -/
theorem stationary_iff_mem_subdiff {lam : ℝ} (hlam : 0 ≤ lam) {s : Site} (hs : s.WellFormed)
    {μ : ℝ} (hμ : μ ∈ s.feasible) (g : ℝ) :
    s.Stationary lam g μ ↔ s.subdiffLo lam μ ≤ (g : EReal) ∧ (g : EReal) ≤ s.subdiffHi lam μ := by
  cases s with
  | penalized wR lo hi =>
    have hf : (fun x => lam * (penalized wR lo hi).cost x) = fun x => lam * wR * |x| := by
      funext x; simp [cost, mul_assoc]
    unfold Stationary
    rw [hf, feasible, isSubgradientWithinAt_abs_iff (mul_nonneg hlam hs.1) hs.2.1 hs.2.2]
    simp only [subdiffLo, subdiffHi]
    constructor
    · rintro ⟨_, h1, h2⟩
      constructor
      · by_cases c : lo < (μ : EReal)
        · rw [ite_eq_left c]; exact_mod_cast h2 c
        · rw [ite_eq_right c]; exact bot_le
      · by_cases c : (μ : EReal) < hi
        · rw [ite_eq_left c]; exact_mod_cast h1 c
        · rw [ite_eq_right c]; exact le_top
    · rintro ⟨h1, h2⟩
      refine ⟨hμ, fun c => ?_, fun c => ?_⟩
      · rw [ite_eq_left c] at h2; exact_mod_cast h2
      · rw [ite_eq_left c] at h1; exact_mod_cast h1
  | fixed a =>
    simp only [subdiffLo, subdiffHi, bot_le, le_top, and_self, iff_true]
    exact (stationary_fixed_iff lam a g μ).mpr hμ

end Site

namespace ConeProblem

variable {ι : Type*} [Fintype ι] {P : ConeProblem ι} {u φ μ : ι → ℝ}

/-- **Shifting the adjoint by a constant on a closed mesh.** If `Lᵀ·1 = 0` and `(u, φ, μ)`
satisfies the state and the adjoint equation, then for every constant `c` the triple
`(u, φ + c, μ)` solves the optimality system iff `φ_i + c` satisfies the inclusion at every vertex.

[formula: adjoint-constant-interval] -/
theorem kkt_add_const_iff (hcol : P.Lᵀ *ᵥ (fun _ => (1 : ℝ)) = 0)
    (hstate : P.L *ᵥ u = μ - P.Ω) (hadj : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i))
    (c : ℝ) :
    P.KKT u (fun i => φ i + c) μ ↔ ∀ i, (P.site i).Stationary P.lam (φ i + c) (μ i) := by
  have e : (P.Lᵀ *ᵥ fun i => φ i + c) = P.Lᵀ *ᵥ φ := by
    have : (fun i => φ i + c) = φ + c • (fun _ => (1 : ℝ)) := by ext i; simp
    rw [this, Matrix.mulVec_add, Matrix.mulVec_smul, hcol]; simp
  constructor
  · intro h; exact h.inclusion
  · intro h; exact ⟨hstate, by rw [e]; exact hadj, h⟩

/-- **The admissible constants form an interval** (`adjoint_constant_interval`): under the
hypotheses of `kkt_add_const_iff`, for a well-formed problem and admissible `μ`, the triple
`(u, φ + c, μ)` solves the optimality system iff `a_i − φ_i ≤ c ≤ b_i − φ_i` for all `i`, where
`[a_i, b_i] = ∂g_i(μ_i)` (`Site.subdiffLo`, `Site.subdiffHi`); i.e. iff
`max_i (a_i − φ_i) ≤ c ≤ min_i (b_i − φ_i)`. With an active vertex (`a_i = b_i`) at most one
constant qualifies.

[formula: adjoint-constant-interval] -/
theorem kkt_add_const_iff_interval (hP : P.WellFormed) (hcol : P.Lᵀ *ᵥ (fun _ => (1 : ℝ)) = 0)
    (hstate : P.L *ᵥ u = μ - P.Ω) (hadj : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i))
    (hμ : ∀ i, μ i ∈ (P.site i).feasible) (c : ℝ) :
    P.KKT u (fun i => φ i + c) μ ↔
      ∀ i, (P.site i).subdiffLo P.lam (μ i) ≤ ((φ i + c : ℝ) : EReal) ∧
        ((φ i + c : ℝ) : EReal) ≤ (P.site i).subdiffHi P.lam (μ i) := by
  rw [kkt_add_const_iff hcol hstate hadj]
  exact forall_congr' fun i =>
    Site.stationary_iff_mem_subdiff hP.lam_nonneg (hP.site i) (hμ i) _

/-! ## The shift of the adjoint constant that restores Gauss–Bonnet (`gauss_bonnet_shift`) -/

/-- **The total cone angle asked for by `φ + c` is nondecreasing in `c`**: for a well-formed
problem and `γ_i > 0`, `c ↦ Σ_i proxReg_{γ_i}(x_i + c)` is monotone. (In the code
`proxReg_{1/k}(x) = clamp(k·P_t(x))` is `Site::measure(k, x)`.) Hence the constants `c` with
`Σ_i measure_i(x_i + c) = Σ Ω` form an interval, and the one of smallest magnitude that
`gauss_bonnet_shift` looks for is well defined when the interval is not empty. The search itself
(breakpoints, bisection) is not formalized.

[formula: gauss-bonnet-adjoint-shift] -/
theorem sum_proxReg_add_monotone (hP : P.WellFormed) {γ : ι → ℝ} (hγ : ∀ i, 0 < γ i)
    (x : ι → ℝ) : Monotone fun c : ℝ => ∑ i, (P.site i).proxReg P.lam (γ i) (x i + c) := by
  intro c d h
  exact Finset.sum_le_sum fun i _ =>
    Site.proxReg_monotone (hγ i) (hP.site i) (by linarith)

/-! ## The slope of the regularized dual objective along a Newton direction (`DualSlope`) -/

/-- The quantity `g(α)` the line search of the regularized stage evaluates (`DualSlope::at`):
`Σ w_i du_i (u_i + α du_i) − Σ Ω_i dφ_i + Σ m_i(x_i + α dφ_i)·dφ_i`, where `m_i` is the cone
angle the adjoint asks for at vertex `i` (`Site::measure`) and `x_i = φ_i + μ̄_i/k_i`. -/
def dualSlope (w Ω : ι → ℝ) (m : ι → ℝ → ℝ) (u du x dφ : ι → ℝ) (α : ℝ) : ℝ :=
  ∑ i, w i * du i * (u i + α * du i) - ∑ i, Ω i * dφ i + ∑ i, m i (x i + α * dφ i) * dφ i

/-- **The line-search slope is minus the pairing of the direction with the state residual.**
Let `A` be symmetric and let the direction satisfy the adjoint relation `w∘du + A dφ = 0`. Then

`g(α) = −Σ_i dφ_i · ((A (u + α du))_i − μ^α_i + Ω_i)`,  `μ^α_i = m_i(x_i + α dφ_i)`:

the derivative, along the Newton direction, of a function whose gradient with respect to `φ` is
minus the residual of the state equation. A zero of `g` is a point where that residual is
orthogonal to the direction.

[formula: regularized-dual-objective] -/
theorem dualSlope_eq {A : Matrix ι ι ℝ} (hA : Aᵀ = A) (w Ω : ι → ℝ) (m : ι → ℝ → ℝ)
    (u du x dφ : ι → ℝ) (hdir : ∀ i, w i * du i + (A *ᵥ dφ) i = 0) (α : ℝ) :
    dualSlope w Ω m u du x dφ α =
      -∑ i, dφ i * ((A *ᵥ fun j => u j + α * du j) i - m i (x i + α * dφ i) + Ω i) := by
  have h1 : ∑ i, w i * du i * (u i + α * du i) =
      -∑ i, dφ i * (A *ᵥ fun j => u j + α * du j) i := by
    have e1 : ∑ i, w i * du i * (u i + α * du i) =
        -((A *ᵥ dφ) ⬝ᵥ fun j => u j + α * du j) := by
      simp only [dotProduct]
      rw [← Finset.sum_neg_distrib]
      refine Finset.sum_congr rfl fun i _ => ?_
      have := hdir i
      have e : w i * du i = -(A *ᵥ dφ) i := by linarith
      rw [e]; ring
    rw [e1]
    have e2 : (A *ᵥ dφ) ⬝ᵥ (fun j => u j + α * du j) =
        dφ ⬝ᵥ (A *ᵥ fun j => u j + α * du j) := by
      rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hA]
    rw [e2]; simp [dotProduct]
  unfold dualSlope
  rw [h1]
  simp only [mul_add, mul_sub, Finset.sum_add_distrib, Finset.sum_sub_distrib]
  have e3 : ∑ i, Ω i * dφ i = ∑ i, dφ i * Ω i := Finset.sum_congr rfl fun i _ => mul_comm _ _
  have e4 : ∑ i, m i (x i + α * dφ i) * dφ i = ∑ i, dφ i * m i (x i + α * dφ i) :=
    Finset.sum_congr rfl fun i _ => mul_comm _ _
  rw [e3, e4]; ring

/-- **The line-search slope is nondecreasing in the step length** when `w ≥ 0` and every `m_i` is
nondecreasing (which `Site::measure` is, `Site.proxReg_monotone`). So the function it is the
derivative of is convex along the direction, and the root the line search computes is its
minimizer on the segment.

[formula: regularized-dual-objective] -/
theorem dualSlope_monotone {w Ω : ι → ℝ} (hw : ∀ i, 0 ≤ w i) {m : ι → ℝ → ℝ}
    (hm : ∀ i, Monotone (m i)) (u du x dφ : ι → ℝ) :
    Monotone (dualSlope w Ω m u du x dφ) := by
  intro α β h
  unfold dualSlope
  have h1 : ∑ i, w i * du i * (u i + α * du i) ≤ ∑ i, w i * du i * (u i + β * du i) := by
    refine Finset.sum_le_sum fun i _ => ?_
    have : 0 ≤ w i * du i ^ 2 * (β - α) :=
      mul_nonneg (mul_nonneg (hw i) (sq_nonneg _)) (by linarith)
    nlinarith
  have h2 : ∑ i, m i (x i + α * dφ i) * dφ i ≤ ∑ i, m i (x i + β * dφ i) * dφ i := by
    refine Finset.sum_le_sum fun i _ => ?_
    rcases le_total 0 (dφ i) with hd | hd
    · exact mul_le_mul_of_nonneg_right (hm i (by nlinarith)) hd
    · exact mul_le_mul_of_nonpos_right (hm i (by nlinarith)) hd
  linarith

end ConeProblem

/-! ## The componentwise backward error of a linear solve (`linsys.rs`, `backward_error`) -/

section BackwardError

variable {κ ι : Type*} [Fintype ι]

/-- **Oettli–Prager.** For a linear system `K x = b` and `ε ≥ 0`, the residual satisfies
`|b − K x|_i ≤ ε·(|K||x| + |b|)_i` in every row iff `x` solves exactly a system whose entries
differ from those of `K` and `b` by at most the relative amount `ε`:
`(K + δK) x = b + δb` with `|δK_ij| ≤ ε|K_ij|`, `|δb_i| ≤ ε|b_i|`.

So the number `max_i |b − K x|_i / (|K||x| + |b|)_i` that `backward_error` returns (rows with a
zero denominator have zero residual, `residual_eq_zero_of_magnitude_eq_zero`) is the smallest such
`ε`; the solver accepts a solution when it is at most `1e-11`.

[formula: componentwise-backward-error] -/
theorem oettli_prager (K : Matrix κ ι ℝ) (x : ι → ℝ) (b : κ → ℝ) {ε : ℝ} (hε : 0 ≤ ε) :
    (∀ i, |b i - (K *ᵥ x) i| ≤ ε * (∑ j, |K i j| * |x j| + |b i|)) ↔
      ∃ (δK : Matrix κ ι ℝ) (δb : κ → ℝ), (∀ i j, |δK i j| ≤ ε * |K i j|) ∧
        (∀ i, |δb i| ≤ ε * |b i|) ∧ (K + δK) *ᵥ x = b + δb := by
  constructor
  · intro h
    -- θ_i = r_i / m_i, with |θ_i| ≤ ε
    set m : κ → ℝ := fun i => ∑ j, |K i j| * |x j| + |b i| with hm
    have hm0 : ∀ i, 0 ≤ m i := fun i =>
      add_nonneg (Finset.sum_nonneg fun j _ => mul_nonneg (abs_nonneg _) (abs_nonneg _))
        (abs_nonneg _)
    set θ : κ → ℝ := fun i => (b i - (K *ᵥ x) i) / m i with hθ
    have hθε : ∀ i, |θ i| ≤ ε := by
      intro i
      rcases (hm0 i).eq_or_lt with h0 | h0
      · simp [hθ, ← h0, hε]
      · rw [hθ]; simp only
        rw [abs_div, abs_of_pos h0, div_le_iff₀ h0]
        exact h i
    have hθm : ∀ i, θ i * m i = b i - (K *ᵥ x) i := by
      intro i
      rcases (hm0 i).eq_or_lt with h0 | h0
      · have hi : |b i - (K *ᵥ x) i| ≤ ε * m i := h i
        rw [← h0, mul_zero] at hi ⊢
        exact (abs_eq_zero.mp (le_antisymm hi (abs_nonneg _))).symm
      · rw [hθ]; simp only
        field_simp
    refine ⟨Matrix.of fun i j => θ i * |K i j| * SignType.sign (x j),
      fun i => -(θ i * |b i|), fun i j => ?_, fun i => ?_, ?_⟩
    · simp only [Matrix.of_apply, abs_mul, abs_abs]
      have h1 : |((SignType.sign (x j) : SignType) : ℝ)| ≤ 1 := by
        rcases lt_trichotomy (x j) 0 with hx | hx | hx <;> simp [hx, sign_neg, sign_pos]
      calc |θ i| * |K i j| * |((SignType.sign (x j) : SignType) : ℝ)|
          ≤ |θ i| * |K i j| * 1 :=
            mul_le_mul_of_nonneg_left h1 (mul_nonneg (abs_nonneg _) (abs_nonneg _))
        _ ≤ ε * |K i j| := by
            rw [mul_one]; exact mul_le_mul_of_nonneg_right (hθε i) (abs_nonneg _)
    · rw [abs_neg, abs_mul, abs_abs]
      exact mul_le_mul_of_nonneg_right (hθε i) (abs_nonneg _)
    · funext i
      have hsx : ∀ j, ((SignType.sign (x j) : SignType) : ℝ) * x j = |x j| := by
        intro j
        rcases lt_trichotomy (x j) 0 with hx | hx | hx
        · simp [sign_neg hx, abs_of_neg hx]
        · simp [hx]
        · simp [sign_pos hx, abs_of_pos hx]
      have e : ((Matrix.of fun i j => θ i * |K i j| * SignType.sign (x j)) *ᵥ x) i =
          θ i * ∑ j, |K i j| * |x j| := by
        simp only [Matrix.mulVec, dotProduct, Matrix.of_apply, Finset.mul_sum]
        refine Finset.sum_congr rfl fun j _ => ?_
        rw [mul_assoc, mul_assoc, hsx j]
      rw [Matrix.add_mulVec]
      simp only [Pi.add_apply, e]
      have := hθm i
      simp only [hm, mul_add] at this
      linarith
  · rintro ⟨δK, δb, h1, h2, h3⟩ i
    have e : b i - (K *ᵥ x) i = (δK *ᵥ x) i - δb i := by
      have := congrFun h3 i
      rw [Matrix.add_mulVec] at this
      simp only [Pi.add_apply] at this
      linarith
    rw [e]
    have hK : |(δK *ᵥ x) i| ≤ ε * ∑ j, |K i j| * |x j| := by
      simp only [Matrix.mulVec, dotProduct]
      calc |∑ j, δK i j * x j| ≤ ∑ j, |δK i j * x j| := Finset.abs_sum_le_sum_abs _ _
        _ ≤ ∑ j, ε * (|K i j| * |x j|) := Finset.sum_le_sum fun j _ => by
            rw [abs_mul, ← mul_assoc]
            exact mul_le_mul_of_nonneg_right (h1 i j) (abs_nonneg _)
        _ = ε * ∑ j, |K i j| * |x j| := by rw [Finset.mul_sum]
    calc |(δK *ᵥ x) i - δb i| ≤ |(δK *ᵥ x) i| + |δb i| := abs_sub _ _
      _ ≤ ε * ∑ j, |K i j| * |x j| + ε * |b i| := add_le_add hK (h2 i)
      _ = ε * (∑ j, |K i j| * |x j| + |b i|) := by ring

/-- A row whose magnitude `(|K||x| + |b|)_i` is zero has zero residual, so skipping it in the
maximum (as `backward_error` does) loses nothing. -/
theorem residual_eq_zero_of_magnitude_eq_zero (K : Matrix κ ι ℝ) (x : ι → ℝ) (b : κ → ℝ) (i : κ)
    (h : ∑ j, |K i j| * |x j| + |b i| = 0) : b i - (K *ᵥ x) i = 0 := by
  have h1 : 0 ≤ ∑ j, |K i j| * |x j| :=
    Finset.sum_nonneg fun j _ => mul_nonneg (abs_nonneg _) (abs_nonneg _)
  have hb : |b i| = 0 := by linarith [abs_nonneg (b i)]
  have hs : ∑ j, |K i j| * |x j| = 0 := by linarith [abs_nonneg (b i)]
  have hterm : ∀ j, K i j * x j = 0 := by
    intro j
    have := (Finset.sum_eq_zero_iff_of_nonneg
      (fun j _ => mul_nonneg (abs_nonneg (K i j)) (abs_nonneg (x j)))).mp hs j (Finset.mem_univ j)
    rw [← abs_mul] at this
    exact abs_eq_zero.mp this
  simp only [Matrix.mulVec, dotProduct, hterm, Finset.sum_const_zero, abs_eq_zero.mp hb, sub_zero]

end BackwardError

/-! ## Merging a cluster of cones (`extract.rs`) -/

section Merge

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The cone angles after merging the cluster `C` into one cone at `v`: the sum of the cluster's
angles at `v`, zero at the other vertices of the cluster, unchanged elsewhere. -/
def mergeCones (C : Finset ι) (v : ι) (μ : ι → ℝ) : ι → ℝ :=
  fun i => if i = v then ∑ j ∈ C, μ j else if i ∈ C then 0 else μ i

/-- **Merging keeps the total cone angle**, whichever vertex `v` of the cluster receives the cone.
So on a closed mesh `Σ μ = Σ Ω = 2πχ` survives merging.

[formula: cone-cluster-frechet-mean] -/
theorem sum_mergeCones {C : Finset ι} {v : ι} (hv : v ∈ C) (μ : ι → ℝ) :
    ∑ i, mergeCones C v μ i = ∑ i, μ i := by
  have h1 : ∑ i, mergeCones C v μ i =
      ∑ i ∈ C, mergeCones C v μ i + ∑ i ∈ Cᶜ, mergeCones C v μ i :=
    (Finset.sum_add_sum_compl C _).symm
  have h2 : ∑ i, μ i = ∑ i ∈ C, μ i + ∑ i ∈ Cᶜ, μ i := (Finset.sum_add_sum_compl C _).symm
  have h3 : ∑ i ∈ Cᶜ, mergeCones C v μ i = ∑ i ∈ Cᶜ, μ i := by
    refine Finset.sum_congr rfl fun i hi => ?_
    have hiC : i ∉ C := Finset.mem_compl.mp hi
    have hiv : i ≠ v := fun h => hiC (h ▸ hv)
    simp [mergeCones, hiC, hiv]
  have h4 : ∑ i ∈ C, mergeCones C v μ i = ∑ j ∈ C, μ j := by
    rw [Finset.sum_eq_single_of_mem v hv]
    · simp [mergeCones]
    · intro i hi hiv; simp [mergeCones, hi, hiv]
  rw [h1, h2, h3, h4]

/-- **Merging does not increase the total cone angle `Σ|μ|`** (it lowers it when the cluster
has angles of both signs). -/
theorem sum_abs_mergeCones_le {C : Finset ι} {v : ι} (hv : v ∈ C) (μ : ι → ℝ) :
    ∑ i, |mergeCones C v μ i| ≤ ∑ i, |μ i| := by
  have h1 : ∑ i, |mergeCones C v μ i| =
      ∑ i ∈ C, |mergeCones C v μ i| + ∑ i ∈ Cᶜ, |mergeCones C v μ i| :=
    (Finset.sum_add_sum_compl C _).symm
  have h2 : ∑ i, |μ i| = ∑ i ∈ C, |μ i| + ∑ i ∈ Cᶜ, |μ i| :=
    (Finset.sum_add_sum_compl C _).symm
  have h3 : ∑ i ∈ Cᶜ, |mergeCones C v μ i| = ∑ i ∈ Cᶜ, |μ i| := by
    refine Finset.sum_congr rfl fun i hi => ?_
    have hiC : i ∉ C := Finset.mem_compl.mp hi
    have hiv : i ≠ v := fun h => hiC (h ▸ hv)
    simp [mergeCones, hiC, hiv]
  have h4 : ∑ i ∈ C, |mergeCones C v μ i| = |∑ j ∈ C, μ j| := by
    rw [Finset.sum_eq_single_of_mem v hv]
    · simp [mergeCones]
    · intro i hi hiv; simp [mergeCones, hi, hiv]
  rw [h1, h2, h3, h4]
  have := Finset.abs_sum_le_sum_abs (fun j => μ j) C
  linarith

omit [Fintype ι] in
/-- **What merging does to admissibility.** If `μ` is admissible, the merged angles are
admissible iff the summed angle is admissible at the receiving vertex `v` and a zero angle is
admissible at the other vertices of the cluster (it always is at a penalized vertex of a
well-formed problem; it is not at a prescribed cone with nonzero angle, which is why prescribed
cones are kept out of clusters).

[formula: cone-cluster-frechet-mean] -/
theorem mergeCones_mem_feasible_iff (P : ConeProblem ι) {C : Finset ι} {v : ι} (μ : ι → ℝ)
    (hμ : ∀ i, μ i ∈ (P.site i).feasible) :
    (∀ i, mergeCones C v μ i ∈ (P.site i).feasible) ↔
      (∑ j ∈ C, μ j) ∈ (P.site v).feasible ∧ ∀ i ∈ C, i ≠ v → (0 : ℝ) ∈ (P.site i).feasible := by
  constructor
  · intro h
    refine ⟨by have := h v; simpa [mergeCones] using this, fun i hi hiv => ?_⟩
    have := h i
    simpa [mergeCones, hi, hiv] using this
  · rintro ⟨h1, h2⟩ i
    by_cases hiv : i = v
    · subst hiv; simpa [mergeCones] using h1
    · by_cases hi : i ∈ C
      · simpa [mergeCones, hi, hiv] using h2 i hi hiv
      · simpa [mergeCones, hi, hiv] using hμ i

/-- A zero cone angle is admissible at every penalized vertex of a well-formed problem. -/
theorem zero_mem_feasible_penalized {wR : ℝ} {lo hi : EReal}
    (hs : (Site.penalized wR lo hi).WellFormed) : (0 : ℝ) ∈ (Site.penalized wR lo hi).feasible :=
  ⟨by simpa using hs.2.1, by simpa using hs.2.2⟩

/-- **The merged configuration is feasible for the scale factor recomputed from it.** If the
merged angles are admissible and `u'` solves `L u' = merged − Ω`, then `(u', merged)` is a
feasible pair of the canonical problem (not, in general, a minimizer: its objective is at least
the optimal one). On a closed mesh such a `u'` can exist because merging keeps `Σ μ = Σ Ω`
(`sum_mergeCones`, `ConeProblem.sum_cone_eq_sum_defect`). -/
theorem feasible_mergeCones (P : ConeProblem ι) {C : Finset ι} {v : ι} (μ u' : ι → ℝ)
    (hadm : ∀ i, mergeCones C v μ i ∈ (P.site i).feasible)
    (hu : P.L *ᵥ u' = mergeCones C v μ - P.Ω) : P.Feasible u' (mergeCones C v μ) :=
  ⟨hu, hadm⟩

omit [Fintype ι] [DecidableEq ι] in
/-- **A Fréchet mean of a cluster exists**: for a nonempty cluster and any variance function
(in the code `v ↦ Σ_k |μ_k|·dist(v, c_k)²`) some vertex of the cluster has the least variance.
Which vertex is chosen affects none of the statements above.

[formula: cone-cluster-frechet-mean] -/
theorem exists_frechetMean {C : Finset ι} (hC : C.Nonempty) (variance : ι → ℝ) :
    ∃ v ∈ C, ∀ w ∈ C, variance v ≤ variance w :=
  Finset.exists_min_image C variance hC

omit [DecidableEq ι] in
/-- Thresholding: the kept angles sum to the total minus the dropped ones
(`ExtractionReport::dropped_angle`). -/
theorem sum_filter_abs_gt (μ : ι → ℝ) (threshold : ℝ) :
    ∑ i ∈ Finset.univ.filter (fun i => threshold < |μ i|), μ i =
      ∑ i, μ i - ∑ i ∈ Finset.univ.filter (fun i => ¬threshold < |μ i|), μ i := by
  have := Finset.sum_filter_add_sum_filter_not Finset.univ (fun i => threshold < |μ i|) μ
  linarith

end Merge

/-! ## The scale factor of a given cone configuration on a closed mesh (`evaluate.rs`) -/

section Evaluate

variable {ι : Type*} [Fintype ι]

/-- **Spreading the Gauss–Bonnet defect by area makes the right-hand side admissible**
(`log_scale_of_measure`): with `defect = Σ μ − Σ Ω` and `Σ m ≠ 0`, the right-hand side
`μ_i − Ω_i − defect·m_i/Σ m` sums to zero, which is the condition under which `A u = rhs` can be
solved on a closed mesh (`Solver.sum_eq_zero_of_mulVec_eq`; given it, the equation dropped by the
pinned solve holds, `Solver.mulVec_eq_of_sum_eq`).

[formula: closed-mesh-cone-sum] -/
theorem sum_sub_defect_eq_zero (μ Ω m : ι → ℝ) (hm : ∑ i, m i ≠ 0) :
    ∑ i, (μ i - Ω i - (∑ j, μ j - ∑ j, Ω j) * m i / ∑ j, m j) = 0 := by
  have h : ∑ i, (∑ j, μ j - ∑ j, Ω j) * m i / ∑ j, m j = ∑ j, μ j - ∑ j, Ω j := by
    rw [← Finset.sum_div, ← Finset.mul_sum]
    field_simp
  rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib, h]
  ring

end Evaluate

end DDG
