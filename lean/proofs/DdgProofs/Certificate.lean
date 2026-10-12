import DdgProofs.Duality

/-!
# The certificate when the data are inexact

`DdgProofs.Duality` proves the duality-gap certificate for a pair `(u, μ)` that satisfies the
state equation `L u = μ − Ω` **exactly** and a dual point that satisfies the adjoint equation and
the dual constraints **exactly**. No computed triple `(u, φ, μ)` does. This file says what remains
true for an **arbitrary** triple, in exact real arithmetic (rounding in *evaluating* the formulas
below is not covered anywhere in this package).

Notation: `w_i = wE_i m_i`, `r = L u − μ + Ω` (state residual, `stateResidual`),
`s = Lᵀφ + w∘ũ` (adjoint residual, `adjointResidual`), `J` the objective, `D` the dual objective.

1. **Gap identity for arbitrary points** (`gap_identity_general`):
   `J(u, μ) − D(ũ, φ) = ½ Σ w_i (u_i − ũ_i)² + Σ_i (λ cost_i(μ_i) + conj_i(φ_i) − φ_i μ_i)
                         − Σ_i φ_i r_i + Σ_i s_i u_i`.
   The last two terms are what `gap_identity` assumes to vanish; they have no sign, so the
   computed "gap" `J − D` can be negative (`example` at the end).
2. **Lower bound from any `φ`.** `dualClip φ` clips `φ_i` to the domain of the conjugate of
   vertex `i`; `dualState φ̂ = −Lᵀφ̂/w` makes the adjoint equation hold by construction. Then
   `dualObjectiveOfAdjoint φ̂` is a lower bound for the objective of every feasible point
   (`dual_lower_bound_of_clip`) provided `(Lᵀφ̂)_i = 0` at every vertex with `w_i = 0`. Without that
   proviso the bound holds only up to the term `Σ_{w_i = 0} (Lᵀφ̂)_i v_i`, which depends on the
   competitor `v` (`dualObjectiveOfAdjoint_add_le`) and can be arbitrarily bad (`example`).
3. **Upper bound from any `(u, μ)`.** The only cone vector that satisfies the state equation
   with the given `u` is `μ + r = L u + Ω` (`feasible_iff_eq_inducedCone`). If `μ_i + r_i` is
   admissible at every vertex, `(u, μ + r)` is feasible and
   `J(u, μ + r) ≤ J(u, μ) + λ Σ wR_i |r_i|` (`objective_inducedCone_le`). Admissible means:
   inside `[lo_i, hi_i]` at a penalized vertex, and `r_i = a_i − μ_i` at a vertex fixed to `a_i`
   (so `r_i = 0` exactly at a forbidden or prescribed vertex whose `μ_i` is exact). If it fails
   at one vertex, no feasible point has this `u`, and the upper bound is false in general
   (`example`s: a forbidden vertex, a vertex at its bound).
4. **Enclosure** (`certificate_enclosure`, `certificate_enclosure_of_clip`): under 2 and 3,
   `dualObjectiveOfAdjoint φ̂ ≤ optimum ≤ J(u, μ + r) ≤ J(u, μ) + λ Σ wR_i |r_i|`, and `(u, μ + r)`
   is suboptimal by at most the width `J(u, μ) − dualObjectiveOfAdjoint φ̂ + λ Σ wR_i |r_i|`.
5. **What is left when `μ + r` is not admissible.**
   * `repaired_certificate`: with `μ̄_i` the admissible angle nearest to `μ_i + r_i` and
     `ρ = μ + r − μ̄` the part of the residual that cannot be absorbed (nonzero only at the
     offending vertices), `(u, μ̄)` is feasible, and certified with an explicit gap, for the
     problem whose angle defects are `Ω − ρ`.
   * `backward_certificate`: `(u, μ)` itself is feasible for the problem with defects `Ω − r`
     and is suboptimal there by at most the complementarity gap
     `½ Σ w (u − ũ)² + Σ (λ cost(μ) + conj(φ̂) − φ̂ μ)`, a sum of nonnegative terms.
   * `KKT.objective_add_le_of_residual`: for the original problem,
     `optimum ≤ J(u, μ) + Σ φ*_i r_i` with `φ*` the adjoint state of an exact solution. At an
     unbounded penalized vertex `|φ*_i| ≤ λ wR_i` is known in advance (that is item 3); at a
     fixed vertex or an active bound `φ*_i` is not bounded by the data, and that is where the
     rigorous statement stops.
-/

noncomputable section

namespace DdgProofs

open Matrix

namespace Site

open Classical in
/-- `p` clipped to the domain of the conjugate of the per-vertex cost (`Site.DualFeasible`):
`max(a, min(b, p))` with `a = −λ·wR` if the lower bound is infinite (else `a = −∞`) and
`b = λ·wR` if the upper bound is infinite (else `b = +∞`). At a fixed vertex and at a vertex
with two finite bounds nothing is clipped; at an unbounded free site (`wR = 0`) the result is `0`. -/
def dualClip (lam : ℝ) : Site → ℝ → ℝ
  | penalized wR lo hi, p =>
    if lo = ⊥ then max (-(lam * wR)) (if hi = ⊤ then min (lam * wR) p else p)
    else if hi = ⊤ then min (lam * wR) p else p
  | fixed _, p => p

/-- The clipped value lies in the domain of the conjugate. -/
theorem dualClip_dualFeasible {lam : ℝ} (hlam : 0 ≤ lam) {s : Site} (hs : s.WellFormed) (p : ℝ) :
    s.DualFeasible lam (s.dualClip lam p) := by
  cases s with
  | penalized wR lo hi =>
    have ht : 0 ≤ lam * wR := mul_nonneg hlam hs.1
    simp only [dualClip, DualFeasible]
    refine ⟨fun h2 => ?_, fun h1 => ?_⟩
    · rw [ite_eq_left h2]
      by_cases h1 : lo = ⊥
      · rw [ite_eq_left h1]; exact max_le (by linarith) (min_le_left _ _)
      · rw [ite_eq_right h1]; exact min_le_left _ _
    · rw [ite_eq_left h1]; exact le_max_left _ _
  | fixed a => trivial

/-- Clipping does nothing to a value that is already in the domain of the conjugate. -/
theorem dualClip_of_dualFeasible {lam : ℝ} {s : Site} {p : ℝ} (h : s.DualFeasible lam p) :
    s.dualClip lam p = p := by
  cases s with
  | penalized wR lo hi =>
    simp only [dualClip]
    by_cases h1 : lo = ⊥
    · rw [ite_eq_left h1]
      by_cases h2 : hi = ⊤
      · rw [ite_eq_left h2, min_eq_right (h.1 h2), max_eq_right (h.2 h1)]
      · rw [ite_eq_right h2, max_eq_right (h.2 h1)]
    · rw [ite_eq_right h1]
      by_cases h2 : hi = ⊤
      · rw [ite_eq_left h2, min_eq_right (h.1 h2)]
      · rw [ite_eq_right h2]
  | fixed a => rfl

/-- The admissible cone angle nearest to `x`: the clamp to `[lo, hi]` at a penalized vertex,
the fixed angle at a fixed vertex. -/
def proj : Site → ℝ → ℝ
  | penalized _ lo hi, x => clampE lo hi x
  | fixed a, _ => a

theorem proj_mem {s : Site} (hs : s.WellFormed) (x : ℝ) : s.proj x ∈ s.feasible := by
  cases s with
  | penalized wR lo hi => exact clampE_mem hs.2.1 hs.2.2 x
  | fixed a => rfl

theorem proj_of_mem {s : Site} {x : ℝ} (h : x ∈ s.feasible) : s.proj x = x := by
  cases s with
  | penalized wR lo hi => exact clampE_of_mem h.1 h.2
  | fixed a => exact h.symm

/-- The cone cost is subadditive: `cost(x + y) ≤ cost(x) + cost(y)`. -/
theorem cost_add_le {s : Site} (hs : s.WellFormed) (x y : ℝ) :
    s.cost (x + y) ≤ s.cost x + s.cost y := by
  cases s with
  | penalized wR lo hi =>
    simp only [cost]
    have := mul_le_mul_of_nonneg_left (abs_add_le x y) hs.1
    linarith
  | fixed a => simp [cost]

/-- The cone cost is monotone in `|x|`. -/
theorem cost_le_of_abs_le {s : Site} (hs : s.WellFormed) {x y : ℝ} (h : |x| ≤ |y|) :
    s.cost x ≤ s.cost y := by
  cases s with
  | penalized wR lo hi => exact mul_le_mul_of_nonneg_left h hs.1
  | fixed a => simp [cost]

/-- The projection to the admissible set does not increase the distance to an admissible
point. -/
theorem abs_proj_sub_le {s : Site} (hs : s.WellFormed) {x : ℝ} (hx : x ∈ s.feasible) (y : ℝ) :
    |s.proj y - x| ≤ |y - x| := by
  cases s with
  | penalized wR lo hi =>
    obtain ⟨_, hlo, hhi⟩ := hs
    obtain ⟨hx1, hx2⟩ := hx
    simp only [proj]
    by_cases h1 : lo ≤ (y : EReal)
    · by_cases h2 : (y : EReal) ≤ hi
      · rw [clampE_of_mem h1 h2]
      · obtain ⟨r, hr, hc⟩ := clampE_of_hi_le hlo hhi (le_of_lt (not_le.mp h2))
        rw [hc]
        have hxr : x ≤ r := by rw [hr] at hx2; exact_mod_cast hx2
        have hry : r < y := by have := not_le.mp h2; rw [hr] at this; exact_mod_cast this
        rw [abs_of_nonneg (by linarith), abs_of_nonneg (by linarith)]
        linarith
    · obtain ⟨r, hr, hc⟩ := clampE_of_le_lo hlo hhi (le_of_lt (not_le.mp h1))
      rw [hc]
      have hxr : r ≤ x := by rw [hr] at hx1; exact_mod_cast hx1
      have hry : y < r := by have := not_le.mp h1; rw [hr] at this; exact_mod_cast this
      rw [abs_of_nonpos (by linarith), abs_of_nonpos (by linarith)]
      linarith
  | fixed a =>
    simp only [feasible, Set.mem_singleton_iff] at hx
    simp [proj, hx]

end Site

namespace ConeProblem

variable {ι : Type*} [Fintype ι] (P : ConeProblem ι)

/-! ## Residuals and the quantities of the certificate -/

/-- The state residual `r = L u − μ + Ω` of an arbitrary pair `(u, μ)`. -/
def stateResidual (u μ : ι → ℝ) : ι → ℝ := fun i => (P.L *ᵥ u) i - μ i + P.Ω i

/-- The adjoint residual `s = Lᵀφ + wE m ũ` of an arbitrary pair `(ũ, φ)`. -/
def adjointResidual (ũ φ : ι → ℝ) : ι → ℝ := fun i => (P.Lᵀ *ᵥ φ) i + P.wE i * P.m i * ũ i

/-- The cone angles that a scale factor `u` produces: `L u + Ω`. For every `μ` this is
`μ + stateResidual u μ` (`inducedCone_eq`), and `(u, ν)` satisfies the state equation iff
`ν = inducedCone u`. -/
def inducedCone (u : ι → ℝ) : ι → ℝ := fun i => (P.L *ᵥ u) i + P.Ω i

/-- The dual state built from `φ` alone, `ũ_i = −(Lᵀφ)_i/(wE_i m_i)`, so that the adjoint
equation holds by construction wherever `wE_i m_i ≠ 0`. Where `wE_i m_i = 0` it is `0` (division
by zero), and the adjoint residual there is `(Lᵀφ)_i` (`adjointResidual_dualState`). -/
def dualState (φ : ι → ℝ) : ι → ℝ := fun i => -(P.Lᵀ *ᵥ φ) i / (P.wE i * P.m i)

/-- `φ` clipped, vertex by vertex, to the domain of the conjugates (`Site.dualClip`). -/
def dualClip (φ : ι → ℝ) : ι → ℝ := fun i => (P.site i).dualClip P.lam (φ i)

/-- The same problem with other angle defects. -/
def withDefect (Ω' : ι → ℝ) : ConeProblem ι := { P with Ω := Ω' }

/-- The complementarity gap of `(u, φ, μ)`:
`½ Σ w_i (u_i − ũ_i)² + Σ_i (λ cost_i(μ_i) + conj_i(φ_i) − φ_i μ_i)` with `ũ = dualState φ`.
Every term is nonnegative when `μ` is admissible and `φ` lies in the domains of the conjugates
(`complementarityGap_nonneg`); it does not involve `Ω`. -/
def complementarityGap (u φ μ : ι → ℝ) : ℝ :=
  1 / 2 * ∑ i, P.wE i * P.m i * (u i - P.dualState φ i) ^ 2 +
    ∑ i, (P.lam * (P.site i).cost (μ i) + (P.site i).conj P.lam (φ i) - φ i * μ i)

/-- The part of the adjoint equation that `dualState` cannot enforce, paired with a state `v`:
`Σ_{i : wE_i m_i = 0} (Lᵀφ)_i v_i`. It vanishes if all `wE_i m_i > 0`, or if `(Lᵀφ)_i = 0`
wherever `wE_i m_i = 0`. -/
def adjointDefect (φ v : ι → ℝ) : ℝ :=
  ∑ i, if P.wE i * P.m i = 0 then (P.Lᵀ *ᵥ φ) i * v i else 0

/-- The admissible cone vector nearest to `L u + Ω`, vertex by vertex. -/
def repairedCone (u : ι → ℝ) : ι → ℝ := fun i => (P.site i).proj (P.inducedCone u i)

/-- The part `ρ = (L u + Ω) − repairedCone u` of the state residual that cannot be absorbed
into the cone angles. It is zero at every vertex where `L u + Ω` is admissible. -/
def unabsorbed (u : ι → ℝ) : ι → ℝ := fun i => P.inducedCone u i - P.repairedCone u i

variable {P} {u ũ φ μ v ν : ι → ℝ}

theorem inducedCone_eq (P : ConeProblem ι) (u μ : ι → ℝ) :
    P.inducedCone u = fun i => μ i + P.stateResidual u μ i := by
  funext i; simp only [inducedCone, stateResidual]; ring

/-- `(u, μ)` satisfies the state equation iff its state residual vanishes. -/
theorem stateResidual_eq_zero_iff : P.stateResidual u μ = 0 ↔ P.L *ᵥ u = μ - P.Ω := by
  constructor
  · intro h; ext i
    have := congrFun h i
    simp only [stateResidual, Pi.zero_apply] at this
    simp only [Pi.sub_apply]; linarith
  · intro h; funext i
    simp only [stateResidual, Pi.zero_apply, h, Pi.sub_apply]; ring

/-- **There is exactly one candidate repair that keeps `u`.** `(u, ν)` is feasible iff
`ν = L u + Ω` and `L u + Ω` is admissible at every vertex. So if `μ + r` is not admissible at
some vertex, no cone vector makes `u` feasible. -/
theorem feasible_iff_eq_inducedCone :
    P.Feasible u ν ↔ ν = P.inducedCone u ∧ ∀ i, P.inducedCone u i ∈ (P.site i).feasible := by
  constructor
  · rintro ⟨h1, h2⟩
    have e : ν = P.inducedCone u := by
      funext i
      have := congrFun h1 i
      simp only [Pi.sub_apply] at this
      simp only [inducedCone]; linarith
    exact ⟨e, fun i => by rw [← e]; exact h2 i⟩
  · rintro ⟨rfl, h2⟩
    refine ⟨?_, h2⟩
    ext i; simp [inducedCone]

/-! ## 1. The gap identity for arbitrary points -/

/-- **Gap identity for arbitrary points.** For *any* `u, ũ, φ, μ` (no equation assumed), with the
state residual `r = L u − μ + Ω` and the adjoint residual `s = Lᵀφ + wE m ũ`:

`J(u, μ) − D(ũ, φ) = ½ Σ wE_i m_i (u_i − ũ_i)² + Σ_i (λ·cost_i(μ_i) + conj_i(φ_i) − φ_i μ_i)
                      − Σ_i φ_i r_i + Σ_i s_i u_i`.

`gap_identity` is the case `r = 0`, `s = 0`. The two residual terms have no sign.

[formula: duality-gap-certificate] -/
theorem gap_identity_general (P : ConeProblem ι) (u ũ φ μ : ι → ℝ) :
    P.objective u μ - P.dualObjective ũ φ =
      1 / 2 * ∑ i, P.wE i * P.m i * (u i - ũ i) ^ 2 +
        ∑ i, (P.lam * (P.site i).cost (μ i) + (P.site i).conj P.lam (φ i) - φ i * μ i) -
        ∑ i, φ i * P.stateResidual u μ i + ∑ i, P.adjointResidual ũ φ i * u i := by
  have h : ∑ i, φ i * (P.L *ᵥ u) i = ∑ i, (P.Lᵀ *ᵥ φ) i * u i := by
    have := Matrix.dotProduct_mulVec φ P.L u
    rw [← Matrix.mulVec_transpose] at this
    simpa [dotProduct] using this
  have e1 : ∑ i, P.wE i * P.m i * (u i - ũ i) ^ 2 =
      ∑ i, P.wE i * P.m i * u i ^ 2 + ∑ i, P.wE i * P.m i * ũ i ^ 2 -
        2 * ∑ i, P.wE i * P.m i * ũ i * u i := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun i _ => by ring
  have e2 : ∑ i, (P.lam * (P.site i).cost (μ i) + (P.site i).conj P.lam (φ i) - φ i * μ i) =
      P.lam * ∑ i, (P.site i).cost (μ i) + ∑ i, (P.site i).conj P.lam (φ i) - ∑ i, φ i * μ i := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
  have e3 : ∑ i, φ i * P.stateResidual u μ i =
      ∑ i, φ i * (P.L *ᵥ u) i - ∑ i, φ i * μ i + ∑ i, φ i * P.Ω i := by
    rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun i _ => by simp only [stateResidual]; ring
  have e4 : ∑ i, P.adjointResidual ũ φ i * u i =
      ∑ i, (P.Lᵀ *ᵥ φ) i * u i + ∑ i, P.wE i * P.m i * ũ i * u i := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun i _ => by simp only [adjointResidual]; ring
  simp only [objective, dualObjective, distortion, coneCost]
  rw [e1, e2, e3, e4, h]
  ring

/-! ## 2. The dual point built from `φ`, and the lower bound -/

/-- The adjoint residual of the constructed dual point: zero where `wE_i m_i ≠ 0`, and
`(Lᵀφ)_i` where `wE_i m_i = 0`. -/
theorem adjointResidual_dualState (P : ConeProblem ι) (φ : ι → ℝ) (i : ι) :
    P.adjointResidual (P.dualState φ) φ i =
      if P.wE i * P.m i = 0 then (P.Lᵀ *ᵥ φ) i else 0 := by
  unfold adjointResidual dualState
  split_ifs with h
  · rw [h]; simp
  · generalize P.wE i * P.m i = w at h
    field_simp
    ring

/-- The dual objective at the constructed dual point is `dualObjectiveOfAdjoint φ` (a term with
`wE_i m_i = 0` contributes nothing to the quadratic part). -/
theorem dualObjective_dualState (P : ConeProblem ι) (φ : ι → ℝ) :
    P.dualObjective (P.dualState φ) φ = P.dualObjectiveOfAdjoint φ := by
  have e : ∀ i, P.wE i * P.m i * P.dualState φ i ^ 2 =
      (P.Lᵀ *ᵥ φ) i ^ 2 / (P.wE i * P.m i) := by
    intro i
    unfold dualState
    by_cases h : P.wE i * P.m i = 0
    · rw [h]; simp
    · field_simp
  simp only [dualObjective, dualObjectiveOfAdjoint, e]

/-- **The computed gap, exactly.** For arbitrary `u, φ, μ`:

`J(u, μ) − dualObjectiveOfAdjoint φ
   = complementarityGap u φ μ − Σ_i φ_i r_i + Σ_{i : wE_i m_i = 0} (Lᵀφ)_i u_i`,

`r = L u − μ + Ω`. This is the relation between the two forms of the gap that the Rust
certificate reports (`duality_gap` on the left, `complementarity_gap` on the right).

[formula: duality-gap-certificate] -/
theorem objective_sub_dualObjectiveOfAdjoint (P : ConeProblem ι) (u φ μ : ι → ℝ) :
    P.objective u μ - P.dualObjectiveOfAdjoint φ =
      P.complementarityGap u φ μ - ∑ i, φ i * P.stateResidual u μ i + P.adjointDefect φ u := by
  rw [← dualObjective_dualState, gap_identity_general]
  unfold complementarityGap adjointDefect
  congr 1
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [adjointResidual_dualState]
  split_ifs <;> simp

/-- The complementarity gap is a sum of nonnegative terms when `μ` is admissible and `φ` lies in
the domains of the conjugates (no equation is needed). -/
theorem complementarityGap_nonneg (hP : P.WellFormed)
    (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i)) (hμ : ∀ i, μ i ∈ (P.site i).feasible)
    (u : ι → ℝ) : 0 ≤ P.complementarityGap u φ μ := by
  have h1 : 0 ≤ ∑ i, P.wE i * P.m i * (u i - P.dualState φ i) ^ 2 :=
    Finset.sum_nonneg fun i _ =>
      mul_nonneg (mul_nonneg (hP.wE_nonneg i) (hP.m_nonneg i)) (sq_nonneg _)
  have h2 : 0 ≤ ∑ i, (P.lam * (P.site i).cost (μ i) + (P.site i).conj P.lam (φ i) - φ i * μ i) :=
    Finset.sum_nonneg fun i _ => by
      have := Site.fenchel_young hP.lam_nonneg (hP.site i) (hφ i) (hμ i)
      linarith
  unfold complementarityGap
  linarith

theorem adjointDefect_eq_zero (hzero : ∀ i, P.wE i * P.m i = 0 → (P.Lᵀ *ᵥ φ) i = 0)
    (v : ι → ℝ) : P.adjointDefect φ v = 0 := by
  unfold adjointDefect
  refine Finset.sum_eq_zero fun i _ => ?_
  split_ifs with h
  · rw [hzero i h, zero_mul]
  · rfl

/-- **Lower bound, in general.** For a well-formed problem, any `φ` in the domains of the
conjugates and any feasible `(v, ν)`:

`dualObjectiveOfAdjoint φ + Σ_{i : wE_i m_i = 0} (Lᵀφ)_i v_i ≤ J(v, ν)`.

The second term is the price of vertices excluded from the distortion measure (`wE_i = 0`): it
depends on the competitor and cannot be evaluated from the data.

[formula: duality-gap-certificate] -/
theorem dualObjectiveOfAdjoint_add_le (hP : P.WellFormed)
    (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i)) (hv : P.Feasible v ν) :
    P.dualObjectiveOfAdjoint φ + P.adjointDefect φ v ≤ P.objective v ν := by
  have h := P.objective_sub_dualObjectiveOfAdjoint v φ ν
  have hr : ∑ i, φ i * P.stateResidual v ν i = 0 := by
    rw [stateResidual_eq_zero_iff.mpr hv.1]; simp
  have hg := complementarityGap_nonneg hP hφ hv.2 v
  linarith

/-- **Lower bound.** For a well-formed problem and any `φ` in the domains of the conjugates with
`(Lᵀφ)_i = 0` wherever `wE_i m_i = 0` (no condition if all `wE_i m_i > 0`),
`dualObjectiveOfAdjoint φ = −½ Σ (Lᵀφ)_i²/(wE_i m_i) + Σ φ_i Ω_i − Σ conj_i(φ_i)` is at most the
objective of every feasible point.

[formula: duality-gap-certificate] -/
theorem dual_lower_bound (hP : P.WellFormed) (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i))
    (hzero : ∀ i, P.wE i * P.m i = 0 → (P.Lᵀ *ᵥ φ) i = 0) (hv : P.Feasible v ν) :
    P.dualObjectiveOfAdjoint φ ≤ P.objective v ν := by
  have := dualObjectiveOfAdjoint_add_le hP hφ hv
  rw [adjointDefect_eq_zero hzero] at this
  linarith

omit [Fintype ι] in
/-- The clipped vector lies in the domain of the conjugate at every vertex. -/
theorem dualClip_dualFeasible (hP : P.WellFormed) (φ : ι → ℝ) (i : ι) :
    (P.site i).DualFeasible P.lam (P.dualClip φ i) :=
  Site.dualClip_dualFeasible hP.lam_nonneg (hP.site i) _

/-- **Lower bound from an arbitrary `φ`.** Clip `φ` to the domains of the conjugates
(`φ̂ = dualClip φ`). If `(Lᵀφ̂)_i = 0` wherever `wE_i m_i = 0` (no condition if all
`wE_i m_i > 0`), then `dualObjectiveOfAdjoint φ̂` is at most the objective of every feasible
point. Nothing is assumed about `φ`.

[formula: duality-gap-certificate] -/
theorem dual_lower_bound_of_clip (hP : P.WellFormed) (φ : ι → ℝ)
    (hzero : ∀ i, P.wE i * P.m i = 0 → (P.Lᵀ *ᵥ P.dualClip φ) i = 0) (hv : P.Feasible v ν) :
    P.dualObjectiveOfAdjoint (P.dualClip φ) ≤ P.objective v ν :=
  dual_lower_bound hP (dualClip_dualFeasible hP φ) hzero hv

/-! ## 3. The upper bound: absorbing the state residual into the cone angles -/

/-- If `μ + r` is admissible at every vertex, `(u, μ + r)` is feasible. -/
theorem feasible_add_stateResidual
    (hadm : ∀ i, μ i + P.stateResidual u μ i ∈ (P.site i).feasible) :
    P.Feasible u (fun i => μ i + P.stateResidual u μ i) := by
  rw [← P.inducedCone_eq u μ, feasible_iff_eq_inducedCone]
  refine ⟨rfl, fun i => ?_⟩
  rw [P.inducedCone_eq u μ]; exact hadm i

/-- **Cost of absorbing the residual**: `J(u, μ + r) ≤ J(u, μ) + λ Σ_i cost_i(r_i)`, i.e.
`+ λ Σ wR_i |r_i|` over the penalized vertices (a fixed vertex contributes nothing). No
admissibility is needed for this inequality.

[formula: duality-gap-certificate] -/
theorem objective_add_stateResidual_le (hP : P.WellFormed) (u μ : ι → ℝ) :
    P.objective u (fun i => μ i + P.stateResidual u μ i) ≤
      P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) := by
  have h : P.coneCost (fun i => μ i + P.stateResidual u μ i) ≤
      P.coneCost μ + ∑ i, (P.site i).cost (P.stateResidual u μ i) := by
    unfold coneCost
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_le_sum fun i _ => Site.cost_add_le (hP.site i) _ _
  have := mul_le_mul_of_nonneg_left h hP.lam_nonneg
  unfold objective
  linarith

/-- The cost of absorbing the residual in terms of its maximum norm: if `|r_i| ≤ ε` at every
vertex then `λ Σ_i cost_i(r_i) ≤ ε · λ Σ_i cost_i(1)`, and `Σ_i cost_i(1)` is the sum of the
weights `wR_i` over the penalized vertices. -/
theorem repairCost_le_of_abs_le (hP : P.WellFormed) {r : ι → ℝ} {ε : ℝ} (hr : ∀ i, |r i| ≤ ε) :
    P.lam * ∑ i, (P.site i).cost (r i) ≤ ε * (P.lam * ∑ i, (P.site i).cost 1) := by
  have h : ∑ i, (P.site i).cost (r i) ≤ ε * ∑ i, (P.site i).cost 1 := by
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum fun i _ => ?_
    have hs := hP.site i
    cases hsite : P.site i with
    | penalized wR lo hi =>
      rw [hsite] at hs
      simp only [Site.cost, abs_one, mul_one]
      have := mul_le_mul_of_nonneg_left (hr i) hs.1
      linarith
    | fixed a => simp [Site.cost]
  have := mul_le_mul_of_nonneg_left h hP.lam_nonneg
  linarith

/-! ## 4. The enclosure -/

/-- **Certificate for inexact data: enclosure of the optimal value.** Let the problem be
well-formed. Let `u`, `μ` be arbitrary with state residual `r = L u − μ + Ω`, and let `φ` be any
vector in the domains of the conjugates with `(Lᵀφ)_i = 0` wherever `wE_i m_i = 0`. Assume
`μ_i + r_i` is admissible at every vertex. Write

* `lower = dualObjectiveOfAdjoint φ`,
* `upper = J(u, μ) + λ Σ_i cost_i(r_i)`  (`cost_i(r_i) = wR_i |r_i|`, `0` at a fixed vertex).

Then

1. `(u, μ + r)` is feasible (exactly);
2. `lower ≤ J(v, ν)` for every feasible `(v, ν)`;
3. `J(u, μ + r) ≤ upper`;
4. `J(u, μ + r) ≤ J(v, ν) + (upper − lower)` for every feasible `(v, ν)`: the repaired point is
   suboptimal by at most the width of the enclosure;
5. `0 ≤ upper − lower`.

In particular the optimal value lies in `[lower, upper]`.

[formula: duality-gap-certificate] -/
theorem certificate_enclosure (hP : P.WellFormed) (u μ φ : ι → ℝ)
    (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i))
    (hzero : ∀ i, P.wE i * P.m i = 0 → (P.Lᵀ *ᵥ φ) i = 0)
    (hadm : ∀ i, μ i + P.stateResidual u μ i ∈ (P.site i).feasible) :
    P.Feasible u (fun i => μ i + P.stateResidual u μ i) ∧
      (∀ v ν, P.Feasible v ν → P.dualObjectiveOfAdjoint φ ≤ P.objective v ν) ∧
      P.objective u (fun i => μ i + P.stateResidual u μ i) ≤
        P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) ∧
      (∀ v ν, P.Feasible v ν →
        P.objective u (fun i => μ i + P.stateResidual u μ i) ≤
          P.objective v ν +
            (P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) -
              P.dualObjectiveOfAdjoint φ)) ∧
      0 ≤ P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) -
        P.dualObjectiveOfAdjoint φ := by
  have hf := feasible_add_stateResidual hadm
  have hup := objective_add_stateResidual_le hP u μ
  have hlow : ∀ v ν, P.Feasible v ν → P.dualObjectiveOfAdjoint φ ≤ P.objective v ν :=
    fun v ν hv => dual_lower_bound hP hφ hzero hv
  refine ⟨hf, hlow, hup, fun v ν hv => ?_, ?_⟩
  · have := hlow v ν hv; linarith
  · have := hlow _ _ hf; linarith

/-- **Certificate for inexact data, from a computed triple.** Let the problem be well-formed
with all `wE_i m_i > 0`. Let `(u, φ, μ)` be *arbitrary*; put `r = L u − μ + Ω` and
`φ̂ = dualClip φ`. If `μ_i + r_i` is admissible at every vertex (inside `[lo_i, hi_i]`; equal to
the fixed angle at a forbidden or prescribed vertex), and

* `J(u, μ) − dualObjectiveOfAdjoint φ̂ ≤ gap`   (the computed duality gap; it may be negative),
* `λ Σ_i cost_i(r_i) ≤ ε`                       (`λ Σ wR_i |r_i|` over the penalized vertices),

then `(u, μ + r)` is feasible, `0 ≤ gap + ε`, and `J(u, μ + r) ≤ J(v, ν) + gap + ε` for every
feasible `(v, ν)`.

[formula: duality-gap-certificate] -/
theorem certificate_of_computed (hP : P.WellFormed) (hpos : ∀ i, 0 < P.wE i * P.m i)
    (u φ μ : ι → ℝ) (hadm : ∀ i, μ i + P.stateResidual u μ i ∈ (P.site i).feasible)
    {gap ε : ℝ} (hgap : P.objective u μ - P.dualObjectiveOfAdjoint (P.dualClip φ) ≤ gap)
    (hres : P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) ≤ ε) :
    P.Feasible u (fun i => μ i + P.stateResidual u μ i) ∧ 0 ≤ gap + ε ∧
      ∀ v ν, P.Feasible v ν →
        P.objective u (fun i => μ i + P.stateResidual u μ i) ≤ P.objective v ν + (gap + ε) := by
  obtain ⟨h1, _, _, h4, h5⟩ := certificate_enclosure hP u μ (P.dualClip φ)
    (dualClip_dualFeasible hP φ) (fun i h => absurd h (hpos i).ne') hadm
  refine ⟨h1, by linarith, fun v ν hv => ?_⟩
  have := h4 v ν hv
  linarith

/-- The enclosure for a minimizer: under the hypotheses of `certificate_enclosure`, the optimal
value `J(u*, μ*)` of any minimizer lies between `dualObjectiveOfAdjoint φ` and
`J(u, μ) + λ Σ_i cost_i(r_i)`. -/
theorem IsMinimizer.objective_mem_enclosure (hP : P.WellFormed) (u μ φ : ι → ℝ)
    (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i))
    (hzero : ∀ i, P.wE i * P.m i = 0 → (P.Lᵀ *ᵥ φ) i = 0)
    (hadm : ∀ i, μ i + P.stateResidual u μ i ∈ (P.site i).feasible)
    {u' μ' : ι → ℝ} (hmin : P.IsMinimizer u' μ') :
    P.dualObjectiveOfAdjoint φ ≤ P.objective u' μ' ∧
      P.objective u' μ' ≤
        P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) := by
  obtain ⟨h1, h2, h3, _, _⟩ := certificate_enclosure hP u μ φ hφ hzero hadm
  exact ⟨h2 _ _ hmin.1, le_trans (hmin.2 _ _ h1) h3⟩

/-- **Distance of the computed scale factor from the optimal one.** Under the hypotheses of
`certificate_enclosure`, every minimizer `(u*, μ*)` satisfies
`¼ Σ_i wE_i m_i (u_i − u*_i)² ≤ upper − lower`, the width of the enclosure. (The scale factor `u`
is the same in the computed pair `(u, μ)` and in the repaired pair `(u, μ + r)`.)

[formula: duality-gap-certificate] -/
theorem IsMinimizer.weighted_sq_sub_le_width (hP : P.WellFormed) (u μ φ : ι → ℝ)
    (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i))
    (hzero : ∀ i, P.wE i * P.m i = 0 → (P.Lᵀ *ᵥ φ) i = 0)
    (hadm : ∀ i, μ i + P.stateResidual u μ i ∈ (P.site i).feasible)
    {u' μ' : ι → ℝ} (hmin : P.IsMinimizer u' μ') :
    1 / 4 * ∑ i, P.wE i * P.m i * (u i - u' i) ^ 2 ≤
      P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) -
        P.dualObjectiveOfAdjoint φ := by
  obtain ⟨h1, h2, h3, _, _⟩ := certificate_enclosure hP u μ φ hφ hzero hadm
  have h4 := hmin.weighted_sq_sub_le hP h1
  have h5 := h2 _ _ hmin.1
  linarith

/-! ## 5. When `μ + r` is not admissible: backward statements -/

omit [Fintype ι] in
theorem withDefect_wellFormed (hP : P.WellFormed) (Ω' : ι → ℝ) : (P.withDefect Ω').WellFormed :=
  ⟨hP.wE_nonneg, hP.m_nonneg, hP.lam_nonneg, hP.site⟩

/-- The dual objective depends on the angle defects only through the term `Σ φ_i Ω_i`. -/
theorem dualObjectiveOfAdjoint_withDefect (P : ConeProblem ι) (ρ φ : ι → ℝ) :
    (P.withDefect fun i => P.Ω i - ρ i).dualObjectiveOfAdjoint φ =
      P.dualObjectiveOfAdjoint φ - ∑ i, φ i * ρ i := by
  simp only [dualObjectiveOfAdjoint, withDefect, mul_sub, Finset.sum_sub_distrib]
  ring

/-- `(u, repairedCone u)` satisfies the constraints of the problem whose defects are `Ω − ρ`,
`ρ = unabsorbed u`. -/
theorem feasible_repairedCone (hP : P.WellFormed) (u : ι → ℝ) :
    (P.withDefect fun i => P.Ω i - P.unabsorbed u i).Feasible u (P.repairedCone u) := by
  refine ⟨?_, fun i => Site.proj_mem (hP.site i) _⟩
  ext i
  simp only [withDefect, unabsorbed, inducedCone, Pi.sub_apply]
  ring

/-- Where `L u + Ω` is admissible nothing is left over: `ρ_i = 0`. -/
theorem unabsorbed_eq_zero_of_mem {i : ι} (h : P.inducedCone u i ∈ (P.site i).feasible) :
    P.unabsorbed u i = 0 := by
  simp only [unabsorbed, repairedCone, Site.proj_of_mem h, sub_self]

/-- The part of the residual that cannot be absorbed is at most the residual, when `μ` itself
is admissible: `|ρ_i| ≤ |r_i|`. -/
theorem abs_unabsorbed_le (hP : P.WellFormed) (u μ : ι → ℝ) {i : ι}
    (hμ : μ i ∈ (P.site i).feasible) : |P.unabsorbed u i| ≤ |P.stateResidual u μ i| := by
  have e : P.inducedCone u i = μ i + P.stateResidual u μ i := congrFun (P.inducedCone_eq u μ) i
  simp only [unabsorbed, repairedCone]
  set y := P.inducedCone u i with hy
  -- |y − proj y| ≤ |y − μ_i| because proj y is the nearest admissible point
  have hr : P.stateResidual u μ i = y - μ i := by rw [e]; ring
  rw [hr]
  cases hsite : P.site i with
  | penalized wR lo hi =>
    have hs := hP.site i
    rw [hsite] at hs hμ
    obtain ⟨_, hlo, hhi⟩ := hs
    obtain ⟨hx1, hx2⟩ := hμ
    simp only [Site.proj]
    by_cases h1 : lo ≤ (y : EReal)
    · by_cases h2 : (y : EReal) ≤ hi
      · rw [clampE_of_mem h1 h2]; simp
      · obtain ⟨c, hc, hcl⟩ := clampE_of_hi_le hlo hhi (le_of_lt (not_le.mp h2))
        rw [hcl]
        have hxc : μ i ≤ c := by rw [hc] at hx2; exact_mod_cast hx2
        have hcy : c < y := by have := not_le.mp h2; rw [hc] at this; exact_mod_cast this
        rw [abs_of_nonneg (by linarith), abs_of_nonneg (by linarith)]
        linarith
    · obtain ⟨c, hc, hcl⟩ := clampE_of_le_lo hlo hhi (le_of_lt (not_le.mp h1))
      rw [hcl]
      have hxc : c ≤ μ i := by rw [hc] at hx1; exact_mod_cast hx1
      have hcy : y < c := by have := not_le.mp h1; rw [hc] at this; exact_mod_cast this
      rw [abs_of_nonpos (by linarith), abs_of_nonpos (by linarith)]
      linarith
  | fixed a =>
    rw [hsite] at hμ
    simp only [Site.feasible, Set.mem_singleton_iff] at hμ
    simp [Site.proj, hμ]

/-- **Certificate for inexact data when the residual cannot be absorbed everywhere.** Let the
problem be well-formed, `u`, `μ` arbitrary with every `μ_i` admissible, `r = L u − μ + Ω`, and
`φ` in the domains of the conjugates with `(Lᵀφ)_i = 0` wherever `wE_i m_i = 0`. Let
`μ̄ = repairedCone u` (the admissible vector nearest to `μ + r`) and `ρ = unabsorbed u = μ + r − μ̄`
(zero wherever `μ_i + r_i` is admissible, and `|ρ_i| ≤ |r_i|` everywhere). Then, for the problem
`P'` that differs from `P` only in the angle defects, `Ω' = Ω − ρ`:

1. `(u, μ̄)` is feasible for `P'`;
2. with `width = J(u, μ) + λ Σ_i cost_i(r_i) − dualObjectiveOfAdjoint φ + Σ_i φ_i ρ_i`:
   `0 ≤ width` and `J(u, μ̄) ≤ J(v, ν) + width` for every `(v, ν)` feasible for `P'`.

If `ρ = 0` this is `certificate_enclosure`. If `ρ ≠ 0` (a forbidden or prescribed vertex with a
nonzero residual, or a vertex at a bound whose residual points outward) the statement is about
the perturbed problem `P'` only.

[formula: duality-gap-certificate] -/
theorem repaired_certificate (hP : P.WellFormed) (u μ φ : ι → ℝ)
    (hμ : ∀ i, μ i ∈ (P.site i).feasible)
    (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i))
    (hzero : ∀ i, P.wE i * P.m i = 0 → (P.Lᵀ *ᵥ φ) i = 0) :
    (P.withDefect fun i => P.Ω i - P.unabsorbed u i).Feasible u (P.repairedCone u) ∧
      0 ≤ P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) -
        P.dualObjectiveOfAdjoint φ + ∑ i, φ i * P.unabsorbed u i ∧
      ∀ v ν, (P.withDefect fun i => P.Ω i - P.unabsorbed u i).Feasible v ν →
        P.objective u (P.repairedCone u) ≤
          P.objective v ν +
            (P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) -
              P.dualObjectiveOfAdjoint φ + ∑ i, φ i * P.unabsorbed u i) := by
  set P' := P.withDefect fun i => P.Ω i - P.unabsorbed u i with hP'
  have hf := feasible_repairedCone hP u
  have hlow : ∀ v ν, P'.Feasible v ν →
      P.dualObjectiveOfAdjoint φ - ∑ i, φ i * P.unabsorbed u i ≤ P.objective v ν := by
    intro v ν hv
    have := dual_lower_bound (P := P') (withDefect_wellFormed hP _) hφ hzero hv
    rw [hP', dualObjectiveOfAdjoint_withDefect] at this
    exact this
  -- cost of moving μ to the repaired cone vector
  have hcost : P.objective u (P.repairedCone u) ≤
      P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) := by
    have h : P.coneCost (P.repairedCone u) ≤
        P.coneCost μ + ∑ i, (P.site i).cost (P.stateResidual u μ i) := by
      unfold coneCost
      rw [← Finset.sum_add_distrib]
      refine Finset.sum_le_sum fun i _ => ?_
      have h1 : P.repairedCone u i = μ i + (P.repairedCone u i - μ i) := by ring
      have h2 := Site.cost_add_le (hP.site i) (μ i) (P.repairedCone u i - μ i)
      have h3 : |P.repairedCone u i - μ i| ≤ |P.stateResidual u μ i| := by
        have e : P.inducedCone u i = μ i + P.stateResidual u μ i :=
          congrFun (P.inducedCone_eq u μ) i
        have := Site.abs_proj_sub_le (hP.site i) (hμ i) (P.inducedCone u i)
        rw [e] at this
        simpa [repairedCone, e] using this
      have h4 := Site.cost_le_of_abs_le (hP.site i) h3
      rw [← h1] at h2
      linarith
    have := mul_le_mul_of_nonneg_left h hP.lam_nonneg
    unfold objective
    linarith
  refine ⟨hf, ?_, fun v ν hv => ?_⟩
  · have := hlow _ _ hf
    have e : P'.objective u (P.repairedCone u) = P.objective u (P.repairedCone u) := rfl
    linarith
  · have := hlow v ν hv
    linarith

/-- **Backward certificate.** Let the problem be well-formed, `u`, `μ` arbitrary with every
`μ_i` admissible, `r = L u − μ + Ω`, and `φ` in the domains of the conjugates with
`(Lᵀφ)_i = 0` wherever `wE_i m_i = 0`. Then `(u, μ)` itself is feasible for the problem `P'` with
angle defects `Ω − r`, the complementarity gap is nonnegative, and

`J(u, μ) ≤ J(v, ν) + complementarityGap u φ μ`   for every `(v, ν)` feasible for `P'`.

No condition on `r`: this covers forbidden and prescribed vertices and active bounds. The
statement is about the perturbed problem; `J` does not depend on `Ω`.

[formula: duality-gap-certificate] -/
theorem backward_certificate (hP : P.WellFormed) (u μ φ : ι → ℝ)
    (hμ : ∀ i, μ i ∈ (P.site i).feasible)
    (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i))
    (hzero : ∀ i, P.wE i * P.m i = 0 → (P.Lᵀ *ᵥ φ) i = 0) :
    (P.withDefect fun i => P.Ω i - P.stateResidual u μ i).Feasible u μ ∧
      0 ≤ P.complementarityGap u φ μ ∧
      ∀ v ν, (P.withDefect fun i => P.Ω i - P.stateResidual u μ i).Feasible v ν →
        P.objective u μ ≤ P.objective v ν + P.complementarityGap u φ μ := by
  set P' := P.withDefect fun i => P.Ω i - P.stateResidual u μ i with hP'
  have hf : P'.Feasible u μ := by
    refine ⟨?_, hμ⟩
    ext i
    simp only [hP', withDefect, stateResidual, Pi.sub_apply]
    ring
  have hgap := P.objective_sub_dualObjectiveOfAdjoint u φ μ
  rw [adjointDefect_eq_zero hzero] at hgap
  refine ⟨hf, complementarityGap_nonneg hP hφ hμ u, fun v ν hv => ?_⟩
  have := dual_lower_bound (P := P') (withDefect_wellFormed hP _) hφ hzero hv
  rw [hP', dualObjectiveOfAdjoint_withDefect] at this
  have e : P'.objective v ν = P.objective v ν := rfl
  linarith

/-- **Sensitivity to the state residual, with the exact multiplier.** Let `(u*, φ*, μ*)` solve
the optimality system of a well-formed problem, and let `(u, μ)` be arbitrary with every `μ_i`
admissible and state residual `r = L u − μ + Ω`. Then

`J(u*, μ*) + ½ Σ wE_i m_i (u_i − u*_i)² ≤ J(u, μ) + Σ_i φ*_i r_i`.

So the optimal value is at most `J(u, μ) + Σ |φ*_i| |r_i|`. At a penalized vertex without bounds
`|φ*_i| ≤ λ wR_i` (`Site.DualFeasible`), which gives back `certificate_enclosure`; at a fixed
vertex and at an active bound `φ*_i` is not bounded by the data of the problem.

[formula: duality-gap-certificate] -/
theorem KKT.objective_add_le_of_residual (hP : P.WellFormed) {u' φ' μ' : ι → ℝ}
    (h : P.KKT u' φ' μ') (u μ : ι → ℝ) (hμ : ∀ i, μ i ∈ (P.site i).feasible) :
    P.objective u' μ' + 1 / 2 * ∑ i, P.wE i * P.m i * (u i - u' i) ^ 2 ≤
      P.objective u μ + ∑ i, φ' i * P.stateResidual u μ i := by
  have hg := P.gap_identity_general u u' φ' μ
  have hs : ∑ i, P.adjointResidual u' φ' i * u i = 0 := by
    refine Finset.sum_eq_zero fun i _ => ?_
    have := congrFun h.adjoint i
    simp only [adjointResidual, this]
    ring
  have hfy : 0 ≤ ∑ i, (P.lam * (P.site i).cost (μ i) + (P.site i).conj P.lam (φ' i) -
      φ' i * μ i) :=
    Finset.sum_nonneg fun i _ => by
      have := Site.fenchel_young hP.lam_nonneg (hP.site i)
        ((h.inclusion i).dualFeasible hP.lam_nonneg (hP.site i)) (hμ i)
      linarith
  have hD := h.dualObjective_eq hP
  linarith

end ConeProblem

/-! ## The hypotheses are satisfiable, and none of them can be dropped

All four instances have one vertex, `L = [1]`, `Ω = 3`, `m = λ = 1`. -/

section Examples

open ConeProblem

private theorem sum_fin_one (f : Fin 1 → ℝ) : ∑ i, f i = f 0 := by simp

/-- **`certificate_of_computed` applies with nonzero residuals.** Default vertex,
`u = −9/10`, `φ = 1`, `μ = 2`: the state residual is `r = 1/10`, `μ + r = 21/10` is admissible,
the computed gap `J − D` is `−19/200 < 0`, the repair cost is `1/10`, and the certified bound
`gap + ε = 1/200` is exactly the suboptimality of `(u, μ + r)`: `J(−9/10, 21/10) = 501/200`,
optimum `5/2`. -/
example : ∃ (P : ConeProblem (Fin 1)) (u φ μ : Fin 1 → ℝ),
    P.WellFormed ∧ (∀ i, 0 < P.wE i * P.m i) ∧
    (∀ i, μ i + P.stateResidual u μ i ∈ (P.site i).feasible) ∧
    P.stateResidual u μ 0 = 1 / 10 ∧
    P.objective u μ - P.dualObjectiveOfAdjoint (P.dualClip φ) = -(19 / 200) ∧
    P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) = 1 / 10 := by
  refine ⟨⟨1, fun _ => 3, fun _ => 1, fun _ => 1, 1, fun _ => Site.plain⟩,
    fun _ => -9 / 10, fun _ => 1, fun _ => 2, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact ⟨fun _ => zero_le_one, fun _ => zero_le_one, zero_le_one,
      fun _ => ⟨zero_le_one, bot_le, le_top⟩⟩
  · intro i; norm_num
  · intro i; simp [Site.plain, Site.feasible, boundSet]
  · simp [stateResidual]; norm_num
  · simp [objective, distortion, coneCost, dualObjectiveOfAdjoint, ConeProblem.dualClip,
      Site.dualClip, Site.plain, Site.cost, Site.conj]
    norm_num
  · simp [stateResidual, Site.plain, Site.cost]; norm_num

/-- **The computed gap can be negative, and `J(u, μ)` can lie below the optimum**, although `μ`
is admissible and `φ` is dual feasible: the same instance, `J(u, μ) = 481/200 < 5/2 =` optimum
(attained at `u = −1`, `μ = 2`). So `J(u, μ) − D` is not a bound on anything by itself. -/
example : ∃ (P : ConeProblem (Fin 1)) (u φ μ u' μ' : Fin 1 → ℝ),
    P.WellFormed ∧ (∀ i, μ i ∈ (P.site i).feasible) ∧
    (∀ i, (P.site i).DualFeasible P.lam (φ i)) ∧ P.IsMinimizer u' μ' ∧
    P.objective u μ < P.objective u' μ' ∧
    P.objective u μ - P.dualObjectiveOfAdjoint φ < 0 := by
  refine ⟨⟨1, fun _ => 3, fun _ => 1, fun _ => 1, 1, fun _ => Site.plain⟩,
    fun _ => -9 / 10, fun _ => 1, fun _ => 2, fun _ => -1, fun _ => 2, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact ⟨fun _ => zero_le_one, fun _ => zero_le_one, zero_le_one,
      fun _ => ⟨zero_le_one, bot_le, le_top⟩⟩
  · intro i; simp [Site.plain, Site.feasible, boundSet]
  · intro i; simp [Site.plain, Site.DualFeasible]
  · refine KKT.isMinimizer (φ := fun _ => 1) (fun _ => zero_le_one) (fun _ => zero_le_one)
      ⟨?_, ?_, fun i => ?_⟩
    · ext i; simp; norm_num
    · ext i; simp
    · rw [Site.stationary_plain_iff zero_le_one]
      rw [softThreshold_of_gt zero_le_one (by norm_num)]; norm_num
  · simp [objective, distortion, coneCost, Site.plain, Site.cost]; norm_num
  · simp [objective, distortion, coneCost, dualObjectiveOfAdjoint, Site.plain, Site.cost,
      Site.conj]
    norm_num

/-- **Admissibility of `μ + r` cannot be dropped: a forbidden vertex.** Vertex fixed to `0`.
The only feasible point is `u = −3`, `μ = 0` with `J = 9/2`. For `u = −29/10`, `μ = 0` the
residual is `r = 1/10`, `μ` is admissible, the repair cost `λ·cost(r)` is `0`, and
`J(u, μ) + λ·cost(r) = 841/200 < 9/2`: the would-be upper bound is below the optimum. -/
example : ∃ (P : ConeProblem (Fin 1)) (u μ u' μ' : Fin 1 → ℝ),
    P.WellFormed ∧ (∀ i, μ i ∈ (P.site i).feasible) ∧ P.IsMinimizer u' μ' ∧
    P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) <
      P.objective u' μ' := by
  refine ⟨⟨1, fun _ => 3, fun _ => 1, fun _ => 1, 1, fun _ => Site.forbidden⟩,
    fun _ => -29 / 10, fun _ => 0, fun _ => -3, fun _ => 0, ?_, ?_, ?_, ?_⟩
  · exact ⟨fun _ => zero_le_one, fun _ => zero_le_one, zero_le_one, fun _ => trivial⟩
  · intro i; simp [Site.forbidden, Site.feasible]
  · refine KKT.isMinimizer (φ := fun _ => 3) (fun _ => zero_le_one) (fun _ => zero_le_one)
      ⟨?_, ?_, fun i => ?_⟩
    · ext i; simp
    · ext i; simp
    · exact (Site.stationary_fixed_iff _ _ _ _).mpr rfl
  · simp [objective, distortion, coneCost, Site.forbidden, Site.cost]; norm_num

/-- **Admissibility of `μ + r` cannot be dropped: a vertex at its bound.** Default cost, upper
bound `hi = 1`. The minimizer is `u = −2`, `μ = 1` with `J = 3`. For `u = −19/10`, `μ = 1` the
residual is `r = 1/10`, `μ + r` exceeds the bound, and `J(u, μ) + λ·wR·|r| = 581/200 < 3`. -/
example : ∃ (P : ConeProblem (Fin 1)) (u μ u' μ' : Fin 1 → ℝ),
    P.WellFormed ∧ (∀ i, μ i ∈ (P.site i).feasible) ∧ P.IsMinimizer u' μ' ∧
    P.objective u μ + P.lam * ∑ i, (P.site i).cost (P.stateResidual u μ i) <
      P.objective u' μ' := by
  refine ⟨⟨1, fun _ => 3, fun _ => 1, fun _ => 1, 1, fun _ => Site.penalized 1 ⊥ (1 : ℝ)⟩,
    fun _ => -19 / 10, fun _ => 1, fun _ => -2, fun _ => 1, ?_, ?_, ?_, ?_⟩
  · exact ⟨fun _ => zero_le_one, fun _ => zero_le_one, zero_le_one,
      fun _ => ⟨zero_le_one, bot_le, by exact_mod_cast zero_le_one⟩⟩
  · intro i; simp [Site.feasible, boundSet]
  · refine ⟨⟨by ext i; simp; norm_num, fun i => by simp [Site.feasible, boundSet]⟩,
      fun v ν hv => ?_⟩
    obtain ⟨h1, h2⟩ := hv
    have e1 := congrFun h1 0
    have e2 : ν 0 ≤ 1 := by
      have := (h2 0).2
      exact_mod_cast this
    simp only [one_mulVec, Pi.sub_apply] at e1
    simp only [objective, distortion, coneCost, sum_fin_one, Site.cost, one_mul, mul_one]
    have h3 : ν 0 ≤ |ν 0| := le_abs_self _
    rw [e1]
    norm_num
    nlinarith [sq_nonneg (ν 0 - 1)]
  · simp [objective, distortion, coneCost, stateResidual, Site.cost]; norm_num

/-- **The proviso on vertices with `wE_i m_i = 0` cannot be dropped.** `wE = 0`, default vertex:
the minimizer is `u = −3`, `μ = 0` with `J = 0`, while `φ = 1` is dual feasible and
`dualObjectiveOfAdjoint φ = 3`. (Here `(Lᵀφ)_0 = 1 ≠ 0`.) -/
example : ∃ (P : ConeProblem (Fin 1)) (φ u' μ' : Fin 1 → ℝ),
    P.WellFormed ∧ (∀ i, (P.site i).DualFeasible P.lam (φ i)) ∧ P.IsMinimizer u' μ' ∧
    P.objective u' μ' < P.dualObjectiveOfAdjoint φ := by
  refine ⟨⟨1, fun _ => 3, fun _ => 1, fun _ => 0, 1, fun _ => Site.plain⟩,
    fun _ => 1, fun _ => -3, fun _ => 0, ?_, ?_, ?_, ?_⟩
  · exact ⟨fun _ => le_rfl, fun _ => zero_le_one, zero_le_one,
      fun _ => ⟨zero_le_one, bot_le, le_top⟩⟩
  · intro i; simp [Site.plain, Site.DualFeasible]
  · refine KKT.isMinimizer (φ := fun _ => 0) (fun _ => le_rfl) (fun _ => zero_le_one)
      ⟨?_, ?_, fun i => ?_⟩
    · ext i; simp
    · ext i; simp
    · rw [Site.stationary_plain_iff zero_le_one]
      rw [softThreshold_of_abs_le (by norm_num)]
  · simp [objective, distortion, coneCost, dualObjectiveOfAdjoint, Site.plain, Site.cost,
      Site.conj]

end Examples

end DdgProofs
