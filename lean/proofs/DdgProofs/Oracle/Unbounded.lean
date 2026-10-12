import DdgProofs.Oracle.Regularized

/-!
# The oracle's theorems for vertices without bounds

At the real numbers the bounds `lo_i`, `hi_i` of `Ddg.Problem ℝ` are finite numbers, so the
theorems of `DdgProofs.Oracle.Problem` and `DdgProofs.Oracle.Regularized` are about problems in
which every cone angle is bounded. The `Float` executable's default is `lo_i = −∞`, `hi_i = +∞`
(and one-sided choices such as `lo_i = 0`, `hi_i = +∞` are allowed). This file extends the
theorems to that case.

An absent bound is described by two predicates on the vertices, `Lo i` ("the lower bound of
vertex `i` is `−∞`") and `Hi i` ("the upper bound is `+∞`"). Where a bound is absent the stored
number `p.lo i` resp. `p.hi i` is ignored. The cone map of such a vertex is
`clampOff (Lo i) (Hi i) lo hi ∘ P_t`, with

`clampOff loInf hiInf lo hi y = y`, `min hi y`, `max lo y` or `min hi (max lo y)`

according to which bounds are absent. This is what the oracle's `clamp lo hi y = min hi (max lo y)`
evaluates to in IEEE arithmetic when `lo = −∞` and/or `hi = +∞` and `y` is a number (that
remark about `Float` is not machine-checked; nothing here is about `Float`).

The reduction that is machine-checked: **the points that occur in a statement lie in a finite
box.** `Ddg.Problem.boxed` replaces every absent bound by a finite one that is out of the way of
two given vectors; on those vectors the oracle's own definitions for the boxed problem
(`coneResidual`, `regConeResidual`, `Feasible`, `objective`, `dualValue`) agree with the
unbounded ones, and the theorems for finite bounds apply.

Results (`Lo`, `Hi` arbitrary; `Lo = Hi = fun _ => False` gives back the finite theorems):

* `objective_le_of_zero_off`: a zero of the residual map is a global minimizer;
* `unique_of_minimizers_off`: uniqueness when all `wE_i·mass_i > 0`;
* `objective_le_add_gap_off`: the duality-gap certificate (at a vertex with an absent bound the
  dual point must satisfy `φ_i ≤ t_i` resp. `−t_i ≤ φ_i`, and the corresponding term of the
  conjugate is absent);
* `regObjective_le_of_zero_off`: a zero of the regularized residual map minimizes the regularized
  objective;
* `zero_off_of_Fpiece_eq_zero`: the stopping rule of the oracle's iteration.
-/

noncomputable section

namespace DdgProofs

namespace Oracle

open Classical in
/-- The clamp to an interval whose lower bound is absent if `loInf` and whose upper bound is
absent if `hiInf`. With both bounds present it is the oracle's `clamp lo hi`. -/
def clampOff (loInf hiInf : Prop) (lo hi y : ℝ) : ℝ :=
  if hiInf then (if loInf then y else max lo y)
  else (if loInf then min hi y else Ddg.clamp lo hi y)

theorem clampOff_false (lo hi y : ℝ) : clampOff False False lo hi y = Ddg.clamp lo hi y := by
  simp [clampOff]

theorem clampOff_true (lo hi y : ℝ) : clampOff True True lo hi y = y := by
  simp [clampOff]

/-- If `μ` is the clamped value of `y`, then `μ ≤ y` when the lower bound is absent and `y ≤ μ`
when the upper bound is absent. -/
theorem le_of_eq_clampOff {loInf hiInf : Prop} {lo hi y μ : ℝ}
    (h : μ = clampOff loInf hiInf lo hi y) : (loInf → μ ≤ y) ∧ (hiInf → y ≤ μ) := by
  unfold clampOff at h
  by_cases h1 : loInf <;> by_cases h2 : hiInf
  · rw [ite_eq_left h2, ite_eq_left h1] at h
    exact ⟨fun _ => h.le, fun _ => h.ge⟩
  · rw [ite_eq_right h2, ite_eq_left h1] at h
    exact ⟨fun _ => by rw [h]; exact min_le_right _ _, fun h' => absurd h' h2⟩
  · rw [ite_eq_left h2, ite_eq_right h1] at h
    exact ⟨fun h' => absurd h' h1, fun _ => by rw [h]; exact le_max_right _ _⟩
  · exact ⟨fun h' => absurd h' h1, fun h' => absurd h' h2⟩

open Classical in
/-- **A finite box that is out of the way acts like an absent bound.** Replace an absent lower
bound by a number `bl ≤ y` and an absent upper bound by a number `bh ≥ y` (and `bh ≥ lo` if the
lower bound is present). Then the oracle's `clamp` with the finite bounds takes the value of
`clampOff` at `y`. -/
theorem clamp_box_eq {loInf hiInf : Prop} {lo hi bl bh y : ℝ} (h1 : loInf → bl ≤ y)
    (h2 : hiInf → y ≤ bh) (h3 : hiInf → ¬loInf → lo ≤ bh) :
    Ddg.clamp (if loInf then bl else lo) (if hiInf then bh else hi) y =
      clampOff loInf hiInf lo hi y := by
  unfold clampOff Ddg.clamp
  by_cases a1 : loInf <;> by_cases a2 : hiInf
  · rw [ite_eq_left a1, ite_eq_left a2, ite_eq_left a2, ite_eq_left a1, max_eq_right (h1 a1),
      min_eq_right (h2 a2)]
  · rw [ite_eq_left a1, ite_eq_right a2, ite_eq_right a2, ite_eq_left a1, max_eq_right (h1 a1)]
  · rw [ite_eq_right a1, ite_eq_left a2, ite_eq_left a2, ite_eq_right a1,
      min_eq_right (max_le (h3 a2 a1) (h2 a2))]
  · rw [ite_eq_right a1, ite_eq_right a2, ite_eq_right a2, ite_eq_right a1]

open Classical in
/-- The conjugate of the per-vertex cost when bounds may be absent:
`hi·max(0, p − t) + lo·min(0, p + t)` without the term of an absent bound. -/
def conjOff (loInf hiInf : Prop) (t lo hi p : ℝ) : ℝ :=
  (if hiInf then 0 else hi * max 0 (p - t)) + (if loInf then 0 else lo * min 0 (p + t))

open Classical in
/-- On the domain of the conjugate (`p ≤ t` if the upper bound is absent, `−t ≤ p` if the lower
one is) the conjugate for any finite replacement of the absent bounds is `conjOff`. -/
theorem vertexConj_box_eq {loInf hiInf : Prop} {t lo hi bl bh p : ℝ} (h1 : hiInf → p ≤ t)
    (h2 : loInf → -t ≤ p) :
    vertexConj t (if loInf then bl else lo) (if hiInf then bh else hi) p =
      conjOff loInf hiInf t lo hi p := by
  unfold vertexConj conjOff
  by_cases a1 : loInf <;> by_cases a2 : hiInf
  · rw [ite_eq_left a1, ite_eq_left a2, ite_eq_left a2, ite_eq_left a1,
      max_eq_left (by linarith [h1 a2]), min_eq_left (by linarith [h2 a1])]
    ring
  · rw [ite_eq_left a1, ite_eq_right a2, ite_eq_right a2, ite_eq_left a1,
      min_eq_left (by linarith [h2 a1])]
    ring
  · rw [ite_eq_right a1, ite_eq_left a2, ite_eq_left a2, ite_eq_right a1,
      max_eq_left (by linarith [h1 a2])]
    ring
  · rw [ite_eq_right a1, ite_eq_right a2, ite_eq_right a2, ite_eq_right a1]

open Classical in
/-- The piece of `x ↦ clampOff loInf hiInf lo hi (k·P_t(x))` at `x`: the oracle's `coneState`
with the test of an absent bound left out. -/
def stateOff (loInf hiInf : Prop) (k t lo hi x : ℝ) : Ddg.ConeState :=
  if ¬hiInf ∧ hi ≤ k * Ddg.softThreshold t x then .upper
  else if ¬loInf ∧ k * Ddg.softThreshold t x ≤ lo then .lower
  else if Ddg.ActivePos t x then .pos
  else if Ddg.ActiveNeg t x then .neg
  else .zero

/-- **The pieces are correct when bounds are absent**: for `t ≥ 0` (and `lo ≤ hi` if both bounds
are present) the map `x ↦ clampOff (k·P_t(x))` equals, at every `x`, the value of the piece
`stateOff … x` at `x`.

[formula: soft-threshold-active-sets] -/
theorem clampOff_eq_value {loInf hiInf : Prop} {t lo hi : ℝ} (ht : 0 ≤ t)
    (hlh : ¬loInf → ¬hiInf → lo ≤ hi) (k x : ℝ) :
    clampOff loInf hiInf lo hi (k * Ddg.softThreshold t x) =
      (stateOff loInf hiInf k t lo hi x).value k t lo hi x := by
  have hs : ∀ s : Ddg.ConeState, s = (if Ddg.ActivePos t x then Ddg.ConeState.pos
      else if Ddg.ActiveNeg t x then Ddg.ConeState.neg else Ddg.ConeState.zero) →
      s.value k t lo hi x = k * Ddg.softThreshold t x := by
    intro s hs
    rw [softThreshold_eq]
    split_ifs at hs with h3 h4
    · subst hs
      simp only [Ddg.ConeState.value]
      rcases h3 with h3 | ⟨h3, h3'⟩
      · rw [softThreshold_of_gt ht h3]
      · simp only [zero_eq] at h3 h3'
        have : t = 0 := le_antisymm h3 ht
        subst this
        rw [softThreshold_zero]; ring
    · subst hs
      simp only [Ddg.ConeState.value]
      rw [softThreshold_of_lt ht h4]
    · subst hs
      unfold Ddg.ActivePos at h3
      unfold Ddg.ActiveNeg at h4
      push Not at h3 h4
      simp only [Ddg.ConeState.value, zero_eq]
      rw [softThreshold_of_abs_le (abs_le.mpr ⟨h4, h3.1⟩), mul_zero]
  unfold stateOff clampOff Ddg.clamp
  by_cases a1 : ¬hiInf ∧ hi ≤ k * Ddg.softThreshold t x
  · rw [ite_eq_left a1, ite_eq_right a1.1]
    simp only [Ddg.ConeState.value]
    by_cases b : loInf
    · rw [ite_eq_left b]; exact min_eq_left a1.2
    · rw [ite_eq_right b]; exact min_eq_left (le_trans a1.2 (le_max_right _ _))
  · rw [ite_eq_right a1]
    by_cases a2 : ¬loInf ∧ k * Ddg.softThreshold t x ≤ lo
    · rw [ite_eq_left a2]
      simp only [Ddg.ConeState.value]
      by_cases b : hiInf
      · rw [ite_eq_left b, ite_eq_right a2.1]; exact max_eq_left a2.2
      · rw [ite_eq_right b, ite_eq_right a2.1, max_eq_left a2.2]; exact min_eq_right (hlh a2.1 b)
    · rw [ite_eq_right a2, hs _ rfl]
      have c1 : ¬hiInf → k * Ddg.softThreshold t x ≤ hi := fun h =>
        le_of_lt (not_le.mp fun h' => a1 ⟨h, h'⟩)
      have c2 : ¬loInf → lo ≤ k * Ddg.softThreshold t x := fun h =>
        le_of_lt (not_le.mp fun h' => a2 ⟨h, h'⟩)
      by_cases b1 : hiInf <;> by_cases b2 : loInf
      · rw [ite_eq_left b1, ite_eq_left b2]
      · rw [ite_eq_left b1, ite_eq_right b2, max_eq_right (c2 b2)]
      · rw [ite_eq_right b1, ite_eq_left b2, min_eq_right (c1 b1)]
      · rw [ite_eq_right b1, ite_eq_right b2, max_eq_right (c2 b2), min_eq_right (c1 b1)]

end Oracle

end DdgProofs

namespace Ddg.Problem

open DdgProofs DdgProofs.Oracle Matrix

/-- Admissible data when the lower bound of the vertices in `Lo` and the upper bound of the
vertices in `Hi` are absent: `L` symmetric, `wE, mass, t ≥ 0`, and at every vertex either each
*present* bound has the right sign (`lo ≤ 0`, `0 ≤ hi`), or both bounds are present and the vertex
is the oracle's encoding of a forbidden or prescribed one (`t = 0`, `lo = hi`). -/
structure AdmissibleOff (p : Ddg.Problem ℝ) (Lo Hi : ℕ → Prop) : Prop where
  L_symm : ∀ i j, i < p.n → j < p.n → p.L i j = p.L j i
  wE_nonneg : ∀ i, i < p.n → 0 ≤ p.wE i
  mass_nonneg : ∀ i, i < p.n → 0 ≤ p.mass i
  t_nonneg : ∀ i, i < p.n → 0 ≤ p.t i
  vertex : ∀ i, i < p.n →
    ((Lo i ∨ p.lo i ≤ 0) ∧ (Hi i ∨ 0 ≤ p.hi i)) ∨
      (¬Lo i ∧ ¬Hi i ∧ p.t i = 0 ∧ p.lo i = p.hi i)

/-- Feasibility when some bounds are absent: the state residual vanishes and every *present*
bound holds. -/
def FeasibleOff (p : Ddg.Problem ℝ) (Lo Hi : ℕ → Prop) (u μ : ℕ → ℝ) : Prop :=
  (∀ i, i < p.n → p.stateResidual u μ i = 0) ∧
    ∀ i, i < p.n → (¬Lo i → p.lo i ≤ μ i) ∧ (¬Hi i → μ i ≤ p.hi i)

/-- The cone residual when some bounds are absent: `μ_i − clampOff(P_{t_i}(φ_i + μ_i))`; for a
vertex without bounds, `μ_i − P_{t_i}(φ_i + μ_i)` (paper Eq. 20, third row). -/
def coneResidualOff (p : Ddg.Problem ℝ) (Lo Hi : ℕ → Prop) (φ μ : ℕ → ℝ) (i : ℕ) : ℝ :=
  μ i - clampOff (Lo i) (Hi i) (p.lo i) (p.hi i) (Ddg.softThreshold (p.t i) (φ i + μ i))

/-- The regularized cone residual when some bounds are absent:
`μ_i − clampOff(c_i·P_{t_i}(φ_i + μ_i))`, `c_i = m_i/(m_i + γ)`. -/
def regConeResidualOff (p : Ddg.Problem ℝ) (Lo Hi : ℕ → Prop) (γ : ℝ) (φ μ : ℕ → ℝ)
    (i : ℕ) : ℝ :=
  μ i - clampOff (Lo i) (Hi i) (p.lo i) (p.hi i)
    (p.coneScale γ i * Ddg.softThreshold (p.t i) (φ i + μ i))

/-- The dual value when some bounds are absent (the terms of absent bounds are left out of the
conjugate). -/
def dualValueOff (p : Ddg.Problem ℝ) (Lo Hi : ℕ → Prop) (ũ φ : ℕ → ℝ) : ℝ :=
  -(1 / 2 * ∑ i : Fin p.n, p.wE i * p.mass i * ũ i ^ 2) + ∑ i : Fin p.n, φ i * p.omega i -
    ∑ i : Fin p.n, conjOff (Lo i) (Hi i) (p.t i) (p.lo i) (p.hi i) (φ i)

open Classical in
/-- The problem with every absent bound replaced by a finite one that is out of the way of the
vectors `a`, `b`: lower bound `min(0, a_i, b_i) − 1`, upper bound `max(0, a_i, b_i) + 1`. -/
def boxed (p : Ddg.Problem ℝ) (Lo Hi : ℕ → Prop) (a b : ℕ → ℝ) : Ddg.Problem ℝ :=
  { p with
    lo := fun i => if Lo i then min 0 (min (a i) (b i)) - 1 else p.lo i
    hi := fun i => if Hi i then max 0 (max (a i) (b i)) + 1 else p.hi i }

variable {p : Ddg.Problem ℝ} {Lo Hi : ℕ → Prop}

/-- With no bound absent, `AdmissibleOff` is `Admissible`. -/
theorem admissibleOff_false_iff (p : Ddg.Problem ℝ) :
    p.AdmissibleOff (fun _ => False) (fun _ => False) ↔ p.Admissible := by
  constructor
  · intro h
    refine ⟨h.L_symm, h.wE_nonneg, h.mass_nonneg, fun i hi => ⟨h.t_nonneg i hi, ?_⟩⟩
    rcases h.vertex i hi with h1 | h1
    · exact Or.inl ⟨h1.1.resolve_left id, h1.2.resolve_left id⟩
    · exact Or.inr ⟨h1.2.2.1, h1.2.2.2⟩
  · intro h
    refine ⟨h.L_symm, h.wE_nonneg, h.mass_nonneg, fun i hi => (h.vertex i hi).1, fun i hi => ?_⟩
    rcases (h.vertex i hi).2 with h1 | h1
    · exact Or.inl ⟨Or.inr h1.1, Or.inr h1.2⟩
    · exact Or.inr ⟨id, id, h1.1, h1.2⟩

private theorem box_lo_le (a b : ℝ) : min 0 (min a b) - 1 ≤ a ∧ min 0 (min a b) - 1 ≤ b ∧
    min 0 (min a b) - 1 ≤ 0 := by
  have h1 : min 0 (min a b) ≤ 0 := min_le_left _ _
  have h2 : min 0 (min a b) ≤ a := le_trans (min_le_right _ _) (min_le_left _ _)
  have h3 : min 0 (min a b) ≤ b := le_trans (min_le_right _ _) (min_le_right _ _)
  exact ⟨by linarith, by linarith, by linarith⟩

private theorem box_hi_ge (a b : ℝ) : a ≤ max 0 (max a b) + 1 ∧ b ≤ max 0 (max a b) + 1 ∧
    0 ≤ max 0 (max a b) + 1 := by
  have h1 : 0 ≤ max 0 (max a b) := le_max_left _ _
  have h2 : a ≤ max 0 (max a b) := le_trans (le_max_left _ _) (le_max_right _ _)
  have h3 : b ≤ max 0 (max a b) := le_trans (le_max_right _ _) (le_max_right _ _)
  exact ⟨by linarith, by linarith, by linarith⟩

/-- The boxed problem has admissible data in the sense of the finite theory. -/
theorem boxed_admissible (hp : p.AdmissibleOff Lo Hi) (a b : ℕ → ℝ) :
    (p.boxed Lo Hi a b).Admissible := by
  refine ⟨hp.L_symm, hp.wE_nonneg, hp.mass_nonneg, fun i hi => ⟨hp.t_nonneg i hi, ?_⟩⟩
  simp only [boxed]
  rcases hp.vertex i hi with ⟨h1, h2⟩ | ⟨h1, h2, h3, h4⟩
  · refine Or.inl ⟨?_, ?_⟩
    · by_cases c : Lo i
      · rw [ite_eq_left c]; exact (box_lo_le _ _).2.2
      · rw [ite_eq_right c]; exact h1.resolve_left c
    · by_cases c : Hi i
      · rw [ite_eq_left c]; exact (box_hi_ge _ _).2.2
      · rw [ite_eq_right c]; exact h2.resolve_left c
  · rw [ite_eq_right h1, ite_eq_right h2]
    exact Or.inr ⟨h3, h4⟩

/-- A pair that is feasible with the bounds absent is feasible for every box that contains its
cone angles. -/
theorem boxed_feasible {v ν a b : ℕ → ℝ} (hv : p.FeasibleOff Lo Hi v ν)
    (hbox : ∀ i, i < p.n → ν i = a i ∨ ν i = b i) : (p.boxed Lo Hi a b).Feasible v ν := by
  refine ⟨hv.1, fun i hi => ?_⟩
  simp only [boxed]
  constructor
  · by_cases c : Lo i
    · rw [ite_eq_left c]
      rcases hbox i hi with h | h <;> rw [h]
      · exact (box_lo_le _ _).1
      · exact (box_lo_le _ _).2.1
    · rw [ite_eq_right c]; exact (hv.2 i hi).1 c
  · by_cases c : Hi i
    · rw [ite_eq_left c]
      rcases hbox i hi with h | h <;> rw [h]
      · exact (box_hi_ge _ _).1
      · exact (box_hi_ge _ _).2.1
    · rw [ite_eq_right c]; exact (hv.2 i hi).2 c

/-- A pair that is feasible for a box is feasible with the bounds absent. -/
theorem feasibleOff_of_boxed {v ν a b : ℕ → ℝ} (hv : (p.boxed Lo Hi a b).Feasible v ν) :
    p.FeasibleOff Lo Hi v ν := by
  refine ⟨hv.1, fun i hi => ⟨fun c => ?_, fun c => ?_⟩⟩
  · have := (hv.2 i hi).1
    simp only [boxed] at this
    rwa [ite_eq_right c] at this
  · have := (hv.2 i hi).2
    simp only [boxed] at this
    rwa [ite_eq_right c] at this

/-- **The oracle's clamp for the boxed problem agrees with the unbounded one** at a value `y`
whose clamped value `μ_i` lies in the box. -/
theorem boxed_clamp_eq (hp : p.AdmissibleOff Lo Hi) {a b μ : ℕ → ℝ} {i : ℕ} (hi : i < p.n)
    (hbox : μ i = a i ∨ μ i = b i) {y : ℝ}
    (hμ : μ i = clampOff (Lo i) (Hi i) (p.lo i) (p.hi i) y) :
    Ddg.clamp ((p.boxed Lo Hi a b).lo i) ((p.boxed Lo Hi a b).hi i) y = μ i := by
  obtain ⟨h1, h2⟩ := le_of_eq_clampOff hμ
  have hl : min 0 (min (a i) (b i)) - 1 ≤ μ i := by
    rcases hbox with h | h <;> rw [h]
    · exact (box_lo_le _ _).1
    · exact (box_lo_le _ _).2.1
  have hh : μ i ≤ max 0 (max (a i) (b i)) + 1 := by
    rcases hbox with h | h <;> rw [h]
    · exact (box_hi_ge _ _).1
    · exact (box_hi_ge _ _).2.1
  simp only [boxed]
  rw [hμ]
  refine clamp_box_eq (fun c => ?_) (fun c => ?_) (fun c c' => ?_)
  · have := h1 c; linarith
  · have := h2 c; linarith
  · rcases hp.vertex i hi with ⟨g1, _⟩ | ⟨_, g2, _⟩
    · have := g1.resolve_left c'
      have := (box_hi_ge (a i) (b i)).2.2
      linarith
    · exact absurd c g2

/-- **Sufficiency on the oracle's definitions, with absent bounds.** Let the data be admissible
(`AdmissibleOff`). If at every unknown the state residual, the adjoint residual and the cone
residual `μ_i − clampOff(P_{t_i}(φ_i + μ_i))` vanish, then `(u, μ)` is feasible and its objective
is at most that of every pair that is feasible with the same bounds absent. With `Lo`, `Hi`
everywhere true this is the paper's unbounded problem: `μ_i = P_{t_i}(φ_i + μ_i)`.

[formula: optimality-system-sufficient] -/
theorem objective_le_of_zero_off (hp : p.AdmissibleOff Lo Hi) {u φ μ v ν : ℕ → ℝ}
    (hstate : ∀ i, i < p.n → p.stateResidual u μ i = 0)
    (hadj : ∀ i, i < p.n → p.adjointResidual u φ i = 0)
    (hcone : ∀ i, i < p.n → p.coneResidualOff Lo Hi φ μ i = 0)
    (hv : p.FeasibleOff Lo Hi v ν) :
    p.FeasibleOff Lo Hi u μ ∧ p.objective u μ ≤ p.objective v ν := by
  have hp' := boxed_admissible hp μ ν
  have hF : ∀ r, r < 3 * (p.boxed Lo Hi μ ν).n → (p.boxed Lo Hi μ ν).F u φ μ r = 0 := by
    rw [F_eq_zero_iff]
    intro i hi
    refine ⟨hstate i hi, hadj i hi, ?_⟩
    have h := hcone i hi
    simp only [coneResidualOff] at h
    simp only [coneResidual, Ddg.boundedThreshold]
    rw [show (p.boxed Lo Hi μ ν).t i = p.t i from rfl,
      boxed_clamp_eq hp hi (Or.inl rfl) (by linarith)]
    ring
  have h := objective_le_of_F_eq_zero hp' hF
    (boxed_feasible hv fun _ _ => Or.inr rfl)
  exact ⟨feasibleOff_of_boxed h.1, h.2⟩

/-- **Uniqueness on the oracle's definitions, with absent bounds.** If all `wE_i·mass_i > 0`, two
pairs that are feasible and minimize the objective (among the pairs feasible with the same
bounds absent) agree at every unknown.

[formula: minimizer-unique] -/
theorem unique_of_minimizers_off (hp : p.AdmissibleOff Lo Hi)
    (hpos : ∀ i, i < p.n → 0 < p.wE i * p.mass i) {u μ v ν : ℕ → ℝ}
    (h₁ : p.FeasibleOff Lo Hi u μ)
    (hmin₁ : ∀ w ω, p.FeasibleOff Lo Hi w ω → p.objective u μ ≤ p.objective w ω)
    (h₂ : p.FeasibleOff Lo Hi v ν)
    (hmin₂ : ∀ w ω, p.FeasibleOff Lo Hi w ω → p.objective v ν ≤ p.objective w ω) :
    ∀ i, i < p.n → u i = v i ∧ μ i = ν i :=
  unique_of_minimizers (boxed_admissible hp μ ν) hpos
    (boxed_feasible h₁ fun _ _ => Or.inl rfl)
    (fun w ω hw => hmin₁ w ω (feasibleOff_of_boxed hw))
    (boxed_feasible h₂ fun _ _ => Or.inr rfl)
    (fun w ω hw => hmin₂ w ω (feasibleOff_of_boxed hw))

/-- **Duality-gap certificate on the oracle's definitions, with absent bounds.** Let the data be
admissible, `(u, μ)` feasible, and let `(ũ, φ)` satisfy the adjoint equation at every unknown and
the dual constraints of the absent bounds: `φ_i ≤ t_i` where the upper bound is absent,
`−t_i ≤ φ_i` where the lower bound is absent (`|φ_i| ≤ t_i` at a vertex without bounds). Then
`dualValueOff ũ φ` is a lower bound for the objective of every feasible point; so if
`objective u μ − dualValueOff ũ φ ≤ gap`, the objective of `(u, μ)` is within `gap` of the optimum.
The hypotheses are exact: see `DdgProofs.Certificate` for what holds for inexact data.

[formula: duality-gap-certificate] -/
theorem objective_le_add_gap_off (hp : p.AdmissibleOff Lo Hi) {u μ ũ φ : ℕ → ℝ}
    (hprimal : p.FeasibleOff Lo Hi u μ)
    (hadj : ∀ i, i < p.n → p.adjointResidual ũ φ i = 0)
    (hφ : ∀ i, i < p.n → (Hi i → φ i ≤ p.t i) ∧ (Lo i → -p.t i ≤ φ i)) {gap : ℝ}
    (hgap : p.objective u μ - p.dualValueOff Lo Hi ũ φ ≤ gap) :
    0 ≤ gap ∧ ∀ v ν, p.FeasibleOff Lo Hi v ν → p.objective u μ ≤ p.objective v ν + gap := by
  have key : ∀ v ν, p.FeasibleOff Lo Hi v ν →
      0 ≤ gap ∧ p.objective u μ ≤ p.objective v ν + gap := by
    intro v ν hv
    have hd : (p.boxed Lo Hi μ ν).dualValue ũ φ = p.dualValueOff Lo Hi ũ φ := by
      simp only [dualValue, dualValueOff]
      congr 1
      refine Finset.sum_congr rfl fun i _ => ?_
      exact vertexConj_box_eq (hφ i i.2).1 (hφ i i.2).2
    have h := objective_le_add_gap (boxed_admissible hp μ ν)
      (boxed_feasible hprimal fun _ _ => Or.inl rfl) hadj (gap := gap)
      (by rw [hd]; exact hgap)
    exact ⟨h.1, h.2 v ν (boxed_feasible hv fun _ _ => Or.inr rfl)⟩
  exact ⟨(key u μ hprimal).1, fun v ν hv => (key v ν hv).2⟩

/-- **A zero of the regularized residual map solves the regularized problem, with absent
bounds.** For admissible data with positive masses and `γ ≥ 0`: if the state residual, the
adjoint residual and the regularized cone residual `μ_i − clampOff(c_i·P_{t_i}(φ_i + μ_i))`,
`c_i = m_i/(m_i + γ)`, vanish at every unknown, then `(u, μ)` is feasible and minimizes
`J(u, μ) + Σ (γ/(2 m_i)) μ_i²` over all pairs feasible with the same bounds absent.

[formula: regularized-optimality-sufficient] -/
theorem regObjective_le_of_zero_off (hp : p.AdmissibleOff Lo Hi)
    (hm : ∀ i, i < p.n → 0 < p.mass i) {γ : ℝ} (hγ : 0 ≤ γ) {u φ μ v ν : ℕ → ℝ}
    (hstate : ∀ i, i < p.n → p.stateResidual u μ i = 0)
    (hadj : ∀ i, i < p.n → p.adjointResidual u φ i = 0)
    (hcone : ∀ i, i < p.n → p.regConeResidualOff Lo Hi γ φ μ i = 0)
    (hv : p.FeasibleOff Lo Hi v ν) :
    p.FeasibleOff Lo Hi u μ ∧ p.regObjective γ u μ ≤ p.regObjective γ v ν := by
  have hp' := boxed_admissible hp μ ν
  have hF : ∀ r, r < 3 * (p.boxed Lo Hi μ ν).n → (p.boxed Lo Hi μ ν).Fimp γ u φ μ r = 0 := by
    rw [Fimp_eq_zero_iff]
    intro i hi
    refine ⟨hstate i hi, hadj i hi, ?_⟩
    have h := hcone i hi
    simp only [regConeResidualOff] at h
    simp only [regConeResidual, Ddg.regularizedMeasure]
    rw [show (p.boxed Lo Hi μ ν).t i = p.t i from rfl,
      show (p.boxed Lo Hi μ ν).coneScale γ i = p.coneScale γ i from rfl,
      boxed_clamp_eq hp hi (Or.inl rfl) (by linarith)]
    ring
  have h := regObjective_le_of_Fimp_eq_zero hp' hm hγ hF
    (boxed_feasible hv fun _ _ => Or.inr rfl)
  exact ⟨feasibleOff_of_boxed h.1, h.2⟩

/-- **Stopping rule for the oracle's piecewise linear systems, with absent bounds.** Let
`t_i ≥ 0` (and `lo_i ≤ hi_i` where both bounds are present). If `(u, φ, μ)` solves the linear
system of the pieces `σ` (`Fpiece γ σ = 0`) and every vertex lies on its piece — the affine
function of `σ_i` takes, at `φ_i + μ_i`, the value of the piece `stateOff` the point is on — then
the state residual, the adjoint residual and the regularized cone residual (`γ = 0`: the exact
cone residual) vanish at every unknown.

The hypothesis `hσ` holds in particular when `σ_i` *is* the piece of the point
(`σ_i = stateOff …`), and when `σ_i` is another name for the same affine function (the pieces
`pos` and `neg` of a free site; a clamped piece whose bound is `0` and the piece `zero`).

[formula: newton-stopping-rule] -/
theorem zero_off_of_Fpiece_eq_zero (ht : ∀ i, i < p.n → 0 ≤ p.t i)
    (hlh : ∀ i, i < p.n → ¬Lo i → ¬Hi i → p.lo i ≤ p.hi i) {γ : ℝ}
    {σ : ℕ → Ddg.ConeState} {u φ μ : ℕ → ℝ}
    (hF : ∀ r, r < 3 * p.n → p.Fpiece γ σ u φ μ r = 0)
    (hσ : ∀ i, i < p.n →
      (σ i).value (p.coneScale γ i) (p.t i) (p.lo i) (p.hi i) (φ i + μ i) =
        (stateOff (Lo i) (Hi i) (p.coneScale γ i) (p.t i) (p.lo i) (p.hi i) (φ i + μ i)).value
          (p.coneScale γ i) (p.t i) (p.lo i) (p.hi i) (φ i + μ i)) :
    ∀ i, i < p.n → p.stateResidual u μ i = 0 ∧ p.adjointResidual u φ i = 0 ∧
      p.regConeResidualOff Lo Hi γ φ μ i = 0 := by
  rw [Fpiece_eq_zero_iff] at hF
  intro i hi
  obtain ⟨h1, h2, h3⟩ := hF i hi
  refine ⟨h1, h2, ?_⟩
  simp only [pieceConeResidual] at h3
  simp only [regConeResidualOff]
  rw [clampOff_eq_value (ht i hi) (hlh i hi), ← hσ i hi]
  exact h3

/-- The regularized cone residual with `γ = 0` is the exact one (for nonzero masses). -/
theorem regConeResidualOff_zero {i : ℕ} (hm : p.mass i ≠ 0) (φ μ : ℕ → ℝ) :
    p.regConeResidualOff Lo Hi 0 φ μ i = p.coneResidualOff Lo Hi φ μ i := by
  simp only [regConeResidualOff, coneResidualOff, coneScale, add_zero, div_self hm, one_mul]

end Ddg.Problem
