import DdgProofs.Oracle.Scalar
import DdgProofs.Duality
import DdgProofs.Closed

/-!
# The optimality theorems, stated about the oracle's `Ddg.Problem`

`Ddg.Problem α` (in `lean/oracle`) is the data of the canonical problem on `n` unknowns, with the
vectors as functions `Nat → α` of which only the indices `< n` matter, the thresholds
`t_i = λ·wR_i`, and **finite** bounds `lo_i`, `hi_i` (in the `Float` executable an absent bound is
`±∞`; at the real numbers a bound is always finite, so the statements below are about problems
with finite bounds — the unbounded problem is covered by the math-level theorems on
`DdgProofs.ConeProblem`).

A vertex is admissible if `t_i ≥ 0` and either `lo_i ≤ 0 ≤ hi_i`, or `t_i = 0` and
`lo_i = hi_i` (the oracle's encoding of a prescribed cone; `lo = hi = 0` is a forbidden vertex).

Results, all about the oracle's own definitions `F`, `stateResidual`, `adjointResidual`,
`coneResidual`, `state`, `objective`:

* `objective_le_of_F_eq_zero`: a zero of the residual map `F` is a global minimizer;
* `unique_of_minimizers`: uniqueness when all `wE_i m_i > 0`;
* `F_eq_zero_of_newton_step`: the stopping rule of the Newton method (Stage II), and
  `Freg_eq_zero_of_newton_step` (Stage I);
* `objective_le_add_gap`: the duality-gap certificate.
-/

noncomputable section

namespace DdgProofs

namespace Oracle

open Matrix

/-! ## Sums and vectors -/

/-- The oracle's `sumTo n f` is `Σ_{i < n} f i`. -/
theorem sumTo_eq (n : ℕ) (f : ℕ → ℝ) : Ddg.sumTo n f = ∑ i : Fin n, f i := by
  induction n with
  | zero => simp [Ddg.sumTo]
  | succ n ih => rw [Ddg.sumTo, ih, Fin.sum_univ_castSucc]; simp

/-- Restriction of an oracle vector to its first `n` entries. -/
def res (n : ℕ) (u : ℕ → ℝ) : Fin n → ℝ := fun i => u i

/-! ## Vertices -/

open Classical in
/-- The vertex data `(t, lo, hi)` of the oracle as a `Site`: a penalized vertex with weight `t`
(`λ` is absorbed, so `λ = 1`) and bounds `[lo, hi]` when `lo ≤ 0 ≤ hi`, otherwise the fixed
value `lo` (the oracle encodes a prescribed cone `c` as `t = 0`, `lo = hi = c`). -/
def siteOf (t lo hi : ℝ) : Site :=
  if lo ≤ 0 ∧ 0 ≤ hi then Site.penalized t (lo : EReal) (hi : EReal) else Site.fixed lo

/-- Admissible vertex data: `t ≥ 0`, and `lo ≤ 0 ≤ hi` or (`t = 0` and `lo = hi`). -/
def AdmissibleVertex (t lo hi : ℝ) : Prop :=
  0 ≤ t ∧ ((lo ≤ 0 ∧ 0 ≤ hi) ∨ (t = 0 ∧ lo = hi))

theorem AdmissibleVertex.le {t lo hi : ℝ} (h : AdmissibleVertex t lo hi) : lo ≤ hi := by
  rcases h.2 with h1 | h1
  · exact le_trans h1.1 h1.2
  · exact h1.2.le

theorem siteOf_wellFormed {t lo hi : ℝ} (h : AdmissibleVertex t lo hi) :
    (siteOf t lo hi).WellFormed := by
  unfold siteOf
  split_ifs with h1
  · exact ⟨h.1, by exact_mod_cast h1.1, by exact_mod_cast h1.2⟩
  · trivial

theorem mem_feasible_siteOf {t lo hi : ℝ} (h : AdmissibleVertex t lo hi) (μ : ℝ) :
    μ ∈ (siteOf t lo hi).feasible ↔ lo ≤ μ ∧ μ ≤ hi := by
  unfold siteOf
  split_ifs with h1
  · simp [Site.feasible, boundSet]
  · rcases h.2 with h2 | h2
    · exact absurd h2 h1
    · simp only [Site.feasible, Set.mem_singleton_iff]
      rw [← h2.2]
      constructor
      · rintro rfl; exact ⟨le_rfl, le_rfl⟩
      · rintro ⟨h3, h4⟩; exact le_antisymm h4 h3

theorem cost_siteOf {t lo hi : ℝ} (h : AdmissibleVertex t lo hi) (μ : ℝ) :
    (siteOf t lo hi).cost μ = t * |μ| := by
  unfold siteOf
  split_ifs with h1
  · rfl
  · rcases h.2 with h2 | h2
    · exact absurd h2 h1
    · simp [Site.cost, h2.1]

theorem prox_siteOf {t lo hi : ℝ} (h : AdmissibleVertex t lo hi) (x : ℝ) :
    (siteOf t lo hi).prox 1 x = Ddg.boundedThreshold t lo hi x := by
  unfold siteOf Ddg.boundedThreshold
  split_ifs with h1
  · simp only [Site.prox, one_mul]
    rw [clamp_eq (le_trans h1.1 h1.2), softThreshold_eq]
  · rcases h.2 with h2 | h2
    · exact absurd h2 h1
    · simp only [Site.prox, Ddg.clamp]
      rw [← h2.2]
      exact (min_eq_left (le_max_left _ _)).symm

/-- The conjugate of the per-vertex cost for finite bounds:
`hi·max(0, p − t) + lo·min(0, p + t)`; for a prescribed cone (`t = 0`, `lo = hi = c`) this is
`c·p`. -/
def vertexConj (t lo hi p : ℝ) : ℝ := hi * max 0 (p - t) + lo * min 0 (p + t)

theorem conj_siteOf {t lo hi : ℝ} (h : AdmissibleVertex t lo hi) (p : ℝ) :
    (siteOf t lo hi).conj 1 p = vertexConj t lo hi p ∧ (siteOf t lo hi).DualFeasible 1 p := by
  unfold siteOf vertexConj
  split_ifs with h1
  · simp [Site.conj, Site.DualFeasible]
  · rcases h.2 with h2 | h2
    · exact absurd h2 h1
    · refine ⟨?_, trivial⟩
      simp only [Site.conj, h2.1, sub_zero, add_zero, ← h2.2]
      rcases le_total 0 p with hp | hp
      · rw [max_eq_right hp, min_eq_left hp]; ring
      · rw [max_eq_left hp, min_eq_right hp]; ring

end Oracle

end DdgProofs

namespace Ddg.Problem

open DdgProofs DdgProofs.Oracle Matrix

/-- The oracle's problem data are admissible: `L` symmetric on the unknowns, `wE ≥ 0`,
`mass ≥ 0`, and every vertex admissible (`t ≥ 0`; `lo ≤ 0 ≤ hi`, or `t = 0` and `lo = hi`). -/
structure Admissible (p : Ddg.Problem ℝ) : Prop where
  L_symm : ∀ i j, i < p.n → j < p.n → p.L i j = p.L j i
  wE_nonneg : ∀ i, i < p.n → 0 ≤ p.wE i
  mass_nonneg : ∀ i, i < p.n → 0 ≤ p.mass i
  vertex : ∀ i, i < p.n → AdmissibleVertex (p.t i) (p.lo i) (p.hi i)

/-- Feasibility in the oracle's terms: the state residual vanishes and the cone angles are
within their bounds, at every unknown. -/
def Feasible (p : Ddg.Problem ℝ) (u μ : ℕ → ℝ) : Prop :=
  (∀ i, i < p.n → p.stateResidual u μ i = 0) ∧ ∀ i, i < p.n → p.lo i ≤ μ i ∧ μ i ≤ p.hi i

/-- The oracle's problem as a `ConeProblem` on `Fin p.n` (with `λ = 1`, thresholds as weights). -/
def toConeProblem (p : Ddg.Problem ℝ) : ConeProblem (Fin p.n) where
  L := Matrix.of fun i j => p.L i j
  Ω := fun i => p.omega i
  m := fun i => p.mass i
  wE := fun i => p.wE i
  lam := 1
  site := fun i => siteOf (p.t i) (p.lo i) (p.hi i)

variable {p : Ddg.Problem ℝ}

theorem applyL_eq (p : Ddg.Problem ℝ) (v : ℕ → ℝ) (i : Fin p.n) :
    p.applyL v i = (p.toConeProblem.L *ᵥ res p.n v) i := by
  simp [applyL, sumTo_eq, toConeProblem, Matrix.mulVec, dotProduct, res]

theorem applyL_add (p : Ddg.Problem ℝ) (v w : ℕ → ℝ) (i : ℕ) :
    p.applyL (fun j => v j + w j) i = p.applyL v i + p.applyL w i := by
  simp only [applyL, sumTo_eq, mul_add, Finset.sum_add_distrib]

theorem toConeProblem_wellFormed (hp : p.Admissible) : p.toConeProblem.WellFormed where
  wE_nonneg := fun i => hp.wE_nonneg i i.2
  m_nonneg := fun i => hp.mass_nonneg i i.2
  lam_nonneg := zero_le_one
  site := fun i => siteOf_wellFormed (hp.vertex i i.2)

theorem objective_eq (hp : p.Admissible) (u μ : ℕ → ℝ) :
    p.objective u μ = p.toConeProblem.objective (res p.n u) (res p.n μ) := by
  simp only [objective, distortion, coneCost, sumTo_eq, ConeProblem.objective,
    ConeProblem.distortion, ConeProblem.coneCost, toConeProblem, res, one_mul, two_eq, abs_eq]
  have h1 : ∑ i : Fin p.n, (siteOf (p.t i) (p.lo i) (p.hi i)).cost (μ i) =
      ∑ i : Fin p.n, p.t i * |μ i| :=
    Finset.sum_congr rfl fun i _ => cost_siteOf (hp.vertex i i.2) _
  rw [h1]
  have h2 : ∑ i : Fin p.n, p.wE i * p.mass i * (u i * u i) =
      ∑ i : Fin p.n, p.wE i * p.mass i * u i ^ 2 :=
    Finset.sum_congr rfl fun i _ => by ring
  rw [h2]
  ring

theorem feasible_iff (hp : p.Admissible) (u μ : ℕ → ℝ) :
    p.Feasible u μ ↔ p.toConeProblem.Feasible (res p.n u) (res p.n μ) := by
  unfold Feasible ConeProblem.Feasible
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨?_, fun i => ?_⟩
    · ext i
      have := h1 i i.2
      rw [stateResidual, applyL_eq p u i] at this
      simp only [Pi.sub_apply, toConeProblem, res] at this ⊢
      linarith
    · exact (mem_feasible_siteOf (hp.vertex i i.2) _).mpr (h2 i i.2)
  · rintro ⟨h1, h2⟩
    refine ⟨fun i hi => ?_, fun i hi => ?_⟩
    · have := congrFun h1 ⟨i, hi⟩
      rw [stateResidual, applyL_eq p u ⟨i, hi⟩]
      simp only [Pi.sub_apply, toConeProblem, res] at this ⊢
      linarith
    · exact (mem_feasible_siteOf (hp.vertex i hi) _).mp (h2 ⟨i, hi⟩)

/-- The flattened residual `F` vanishes iff its three blocks vanish. -/
theorem F_eq_zero_iff (p : Ddg.Problem ℝ) (u φ μ : ℕ → ℝ) :
    (∀ r, r < 3 * p.n → p.F u φ μ r = 0) ↔
      ∀ i, i < p.n → p.stateResidual u μ i = 0 ∧ p.adjointResidual u φ i = 0 ∧
        p.coneResidual φ μ i = 0 := by
  constructor
  · intro h i hi
    refine ⟨?_, ?_, ?_⟩
    · have := h i (by omega)
      simpa [F, hi] using this
    · have := h (p.n + i) (by omega)
      have h1 : ¬(p.n + i < p.n) := by omega
      have h2 : p.n + i < 2 * p.n := by omega
      simpa [F, h1, h2] using this
    · have := h (2 * p.n + i) (by omega)
      have h1 : ¬(2 * p.n + i < p.n) := by omega
      have h2 : ¬(2 * p.n + i < 2 * p.n) := by omega
      simpa [F, h1, h2] using this
  · intro h r hr
    unfold F
    split_ifs with h1 h2
    · exact (h r h1).1
    · exact (h (r - p.n) (by omega)).2.1
    · exact (h (r - 2 * p.n) (by omega)).2.2

/-- A zero of the oracle's residual map is a solution of the optimality system of the
corresponding `ConeProblem`. -/
theorem kkt_of_F_eq_zero (hp : p.Admissible) {u φ μ : ℕ → ℝ}
    (hF : ∀ r, r < 3 * p.n → p.F u φ μ r = 0) :
    p.toConeProblem.KKT (res p.n u) (res p.n φ) (res p.n μ) := by
  rw [F_eq_zero_iff] at hF
  rw [ConeProblem.kkt_iff_fixedPoint _ (toConeProblem_wellFormed hp)]
  refine ⟨?_, ?_, fun i => ?_⟩
  · ext i
    have := (hF i i.2).1
    rw [stateResidual, applyL_eq p u i] at this
    simp only [Pi.sub_apply, toConeProblem, res] at this ⊢
    linarith
  · have hsymm : p.toConeProblem.Lᵀ = p.toConeProblem.L := by
      ext i j
      simp only [toConeProblem, Matrix.transpose_apply, Matrix.of_apply]
      exact hp.L_symm j i j.2 i.2
    rw [hsymm]
    ext i
    have := (hF i i.2).2.1
    rw [adjointResidual, applyL_eq p φ i] at this
    simp only [toConeProblem, res] at this ⊢
    linarith
  · have := (hF i i.2).2.2
    rw [coneResidual] at this
    simp only [toConeProblem, res]
    rw [prox_siteOf (hp.vertex i i.2)]
    linarith

/-- **Sufficiency, on the oracle's definitions.** For admissible data, if the residual map
`F(u, φ, μ)` (state, adjoint and cone residuals of all `n` unknowns) vanishes, then `(u, μ)` is
feasible and its objective is at most that of every feasible `(v, ν)`.

[formula: optimality-system-sufficient] -/
theorem objective_le_of_F_eq_zero (hp : p.Admissible) {u φ μ v ν : ℕ → ℝ}
    (hF : ∀ r, r < 3 * p.n → p.F u φ μ r = 0) (hv : p.Feasible v ν) :
    p.Feasible u μ ∧ p.objective u μ ≤ p.objective v ν := by
  have hk := kkt_of_F_eq_zero hp hF
  have hmin := hk.isMinimizer (toConeProblem_wellFormed hp).wE_nonneg
    (toConeProblem_wellFormed hp).m_nonneg
  refine ⟨(feasible_iff hp u μ).mpr hmin.1, ?_⟩
  rw [objective_eq hp, objective_eq hp]
  exact hmin.2 _ _ ((feasible_iff hp v ν).mp hv)

/-- **Uniqueness, on the oracle's definitions.** If all `wE_i·mass_i > 0`, two feasible pairs
that both minimize the objective agree at every unknown.

[formula: minimizer-unique] -/
theorem unique_of_minimizers (hp : p.Admissible) (hpos : ∀ i, i < p.n → 0 < p.wE i * p.mass i)
    {u μ v ν : ℕ → ℝ} (h₁ : p.Feasible u μ)
    (hmin₁ : ∀ w ω, p.Feasible w ω → p.objective u μ ≤ p.objective w ω)
    (h₂ : p.Feasible v ν)
    (hmin₂ : ∀ w ω, p.Feasible w ω → p.objective v ν ≤ p.objective w ω) :
    ∀ i, i < p.n → u i = v i ∧ μ i = ν i := by
  have ext : ∀ (w ω : Fin p.n → ℝ), p.toConeProblem.Feasible w ω →
      ∃ w' ω' : ℕ → ℝ, res p.n w' = w ∧ res p.n ω' = ω ∧ p.Feasible w' ω' := by
    intro w ω hw
    refine ⟨fun i => if h : i < p.n then w ⟨i, h⟩ else 0,
      fun i => if h : i < p.n then ω ⟨i, h⟩ else 0, ?_, ?_, ?_⟩
    · funext i; simp [res]
    · funext i; simp [res]
    · rw [feasible_iff hp]
      have e1 : res p.n (fun i => if h : i < p.n then w ⟨i, h⟩ else 0) = w := by
        funext i; simp [res]
      have e2 : res p.n (fun i => if h : i < p.n then ω ⟨i, h⟩ else 0) = ω := by
        funext i; simp [res]
      rw [e1, e2]; exact hw
  have m₁ : p.toConeProblem.IsMinimizer (res p.n u) (res p.n μ) := by
    refine ⟨(feasible_iff hp u μ).mp h₁, fun w ω hw => ?_⟩
    obtain ⟨w', ω', rfl, rfl, hf⟩ := ext w ω hw
    rw [← objective_eq hp, ← objective_eq hp]
    exact hmin₁ w' ω' hf
  have m₂ : p.toConeProblem.IsMinimizer (res p.n v) (res p.n ν) := by
    refine ⟨(feasible_iff hp v ν).mp h₂, fun w ω hw => ?_⟩
    obtain ⟨w', ω', rfl, rfl, hf⟩ := ext w ω hw
    rw [← objective_eq hp, ← objective_eq hp]
    exact hmin₂ w' ω' hf
  obtain ⟨e1, e2⟩ := m₁.unique (toConeProblem_wellFormed hp) (fun i => hpos i i.2) m₂
  intro i hi
  exact ⟨congrFun e1 ⟨i, hi⟩, congrFun e2 ⟨i, hi⟩⟩

/-- **Stopping rule of the Newton method, Stage II, on the oracle's definitions.** Let
`(v, q, ν)` solve the linearized system at `(u, φ, μ)` block by block,

* `(L v)_i − ν_i = −stateResidual_i`,
* `(L q)_i + wE_i m_i v_i = −adjointResidual_i`,
* `ν_i − D_i (q_i + ν_i) = −coneResidual_i`, `D_i` the indicator of the active set
  `state φ μ i` (rows of the oracle's `DF`),

and suppose the step does not change the state (active set) of any vertex. Then the residual
map `F` vanishes at the new iterate.

[formula: newton-stopping-rule] -/
theorem F_eq_zero_of_newton_step (hp : p.Admissible) {u φ μ v q ν : ℕ → ℝ}
    (h1 : ∀ i, i < p.n → p.applyL v i - ν i = -p.stateResidual u μ i)
    (h2 : ∀ i, i < p.n → p.applyL q i + p.wE i * p.mass i * v i = -p.adjointResidual u φ i)
    (h3 : ∀ i, i < p.n →
      ν i - (p.state φ μ i).indicator * (q i + ν i) = -p.coneResidual φ μ i)
    (hact : ∀ i, i < p.n →
      p.state (fun j => φ j + q j) (fun j => μ j + ν j) i = p.state φ μ i) :
    ∀ r, r < 3 * p.n →
      p.F (fun j => u j + v j) (fun j => φ j + q j) (fun j => μ j + ν j) r = 0 := by
  rw [F_eq_zero_iff]
  intro i hi
  refine ⟨?_, ?_, ?_⟩
  · have := h1 i hi
    simp only [stateResidual, applyL_add] at this ⊢
    linarith
  · have := h2 i hi
    simp only [adjointResidual, applyL_add] at this ⊢
    linarith
  · have h := h3 i hi
    have ha := hact i hi
    have hv := hp.vertex i hi
    simp only [state, one_eq] at ha h
    simp only [coneResidual, boundedThreshold_eq_regularizedMeasure] at h ⊢
    have hq := regularizedMeasure_isPiecewiseAffine hv.1 hv.le 1 (φ i + μ i)
      ((φ i + q i) + (μ i + ν i)) ha
    rw [hq]
    simp only [one_mul] at h ⊢
    linarith

/-- **Stopping rule, Stage II: the result is optimal** (oracle's definitions): under the
hypotheses of `F_eq_zero_of_newton_step` the new iterate minimizes the objective.

[formula: newton-stopping-rule] -/
theorem objective_le_of_newton_step (hp : p.Admissible) {u φ μ v q ν w ω : ℕ → ℝ}
    (h1 : ∀ i, i < p.n → p.applyL v i - ν i = -p.stateResidual u μ i)
    (h2 : ∀ i, i < p.n → p.applyL q i + p.wE i * p.mass i * v i = -p.adjointResidual u φ i)
    (h3 : ∀ i, i < p.n →
      ν i - (p.state φ μ i).indicator * (q i + ν i) = -p.coneResidual φ μ i)
    (hact : ∀ i, i < p.n →
      p.state (fun j => φ j + q j) (fun j => μ j + ν j) i = p.state φ μ i)
    (hw : p.Feasible w ω) :
    p.objective (fun j => u j + v j) (fun j => μ j + ν j) ≤ p.objective w ω :=
  (objective_le_of_F_eq_zero hp (F_eq_zero_of_newton_step hp h1 h2 h3 hact) hw).2

/-- **Stopping rule of the Newton method, Stage I, on the oracle's definitions** (regularized
residual `Freg`, with the mass factor `m_i/γ`). Let `(v, q)` solve the linearized system at
`(u, φ)`,

* `(L v)_i − (m_i/γ) D_i q_i = −regStateResidual_i`, `D_i` the indicator of `regState γ φ i`,
* `(L q)_i + wE_i m_i v_i = −adjointResidual_i`,

and suppose the step does not change `regState` at any vertex. Then `Freg` vanishes at the new
iterate.

[formula: newton-stopping-rule] -/
theorem Freg_eq_zero_of_newton_step (hp : p.Admissible) (γ : ℝ) {u φ v q : ℕ → ℝ}
    (h1 : ∀ i, i < p.n →
      p.applyL v i - p.mass i / γ * (p.regState γ φ i).indicator * q i =
        -p.regStateResidual γ u φ i)
    (h2 : ∀ i, i < p.n → p.applyL q i + p.wE i * p.mass i * v i = -p.adjointResidual u φ i)
    (hact : ∀ i, i < p.n → p.regState γ (fun j => φ j + q j) i = p.regState γ φ i) :
    ∀ r, r < 2 * p.n → p.Freg γ (fun j => u j + v j) (fun j => φ j + q j) r = 0 := by
  intro r hr
  unfold Freg
  split_ifs with hrn
  · have h := h1 r hrn
    have ha := hact r hrn
    have hv := hp.vertex r hrn
    simp only [regState] at ha h
    simp only [regStateResidual, regMeasure, applyL_add] at h ⊢
    have hq := regularizedMeasure_isPiecewiseAffine hv.1 hv.le (p.mass r / γ) (φ r)
      (φ r + q r) ha
    rw [hq]
    linarith
  · have := h2 (r - p.n) (by omega)
    simp only [adjointResidual, applyL_add] at this ⊢
    linarith

/-! ## The certificate -/

/-- The dual value of the oracle's problem at `(ũ, φ)`:
`−½ Σ wE m ũ² + Σ φ Ω − Σ (hi·max(0, φ − t) + lo·min(0, φ + t))`. With finite bounds every `φ` is
admissible; `ũ` must satisfy the adjoint equation. -/
def dualValue (p : Ddg.Problem ℝ) (ũ φ : ℕ → ℝ) : ℝ :=
  -(1 / 2 * ∑ i : Fin p.n, p.wE i * p.mass i * ũ i ^ 2) + ∑ i : Fin p.n, φ i * p.omega i -
    ∑ i : Fin p.n, vertexConj (p.t i) (p.lo i) (p.hi i) (φ i)

theorem dualValue_eq (hp : p.Admissible) (ũ φ : ℕ → ℝ) :
    p.dualValue ũ φ = p.toConeProblem.dualObjective (res p.n ũ) (res p.n φ) := by
  simp only [dualValue, ConeProblem.dualObjective, toConeProblem, res]
  congr 1
  exact Finset.sum_congr rfl fun i _ => ((conj_siteOf (hp.vertex i i.2) _).1).symm

/-- **Duality-gap certificate, on the oracle's definitions.** Let the data be admissible, let
`(u, μ)` be feasible, and let `(ũ, φ)` satisfy the adjoint equation
`(L φ)_i + wE_i m_i ũ_i = 0` at every unknown. Then `dualValue ũ φ` is a lower bound for the
objective of every feasible point; so if `objective u μ − dualValue ũ φ ≤ gap`, the objective
of `(u, μ)` is within `gap` of the optimum.

[formula: duality-gap-certificate] -/
theorem objective_le_add_gap (hp : p.Admissible) {u μ ũ φ : ℕ → ℝ} (hprimal : p.Feasible u μ)
    (hadj : ∀ i, i < p.n → p.adjointResidual ũ φ i = 0) {gap : ℝ}
    (hgap : p.objective u μ - p.dualValue ũ φ ≤ gap) :
    0 ≤ gap ∧ ∀ v ν, p.Feasible v ν → p.objective u μ ≤ p.objective v ν + gap := by
  have hdual : p.toConeProblem.DualFeasible (res p.n ũ) (res p.n φ) := by
    refine ⟨?_, fun i => (conj_siteOf (hp.vertex i i.2) _).2⟩
    have hsymm : p.toConeProblem.Lᵀ = p.toConeProblem.L := by
      ext i j
      simp only [toConeProblem, Matrix.transpose_apply, Matrix.of_apply]
      exact hp.L_symm j i j.2 i.2
    rw [hsymm]
    ext i
    have := hadj i i.2
    rw [adjointResidual, applyL_eq p φ i] at this
    simp only [toConeProblem, res] at this ⊢
    linarith
  rw [objective_eq hp, dualValue_eq hp] at hgap
  have h := ConeProblem.duality_gap_certificate (toConeProblem_wellFormed hp)
    ((feasible_iff hp u μ).mp hprimal) hdual hgap
  refine ⟨h.1, fun v ν hv => ?_⟩
  rw [objective_eq hp, objective_eq hp]
  exact h.2 _ _ ((feasible_iff hp v ν).mp hv)

end Ddg.Problem
