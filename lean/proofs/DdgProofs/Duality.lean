import DdgProofs.Optimality

/-!
# The discrete dual problem and the duality-gap certificate

The discrete analogue of the paper's pre-dual (Problem 10 / A.1) for the canonical discrete
problem is

```
maximize    D(ũ, φ) = −½ Σ_i wE_i m_i ũ_i²  +  Σ_i φ_i Ω_i  −  Σ_i g_i*(φ_i)
subject to  Lᵀφ = −wE m ũ,   φ_i ∈ dom g_i*
```

where `g_i*` is the convex conjugate of the per-vertex cost `λ·wR_i|·| + ι_[lo_i,hi_i]`
(`Site.conj`) and `dom g_i*` its domain (`Site.DualFeasible`). For the paper's basic problem
(unit weights, no bounds) `g_i* = 0` on `|φ_i| ≤ λ`, and `−D` is exactly the objective
`½ ∫|u|² − ∫ Kφ` of Problem 10 with the constraint `|φ| ≤ λ`.

Main results: the gap identity (`gap_identity`), weak duality (`weak_duality`), the certificate
(`duality_gap_certificate`), strong duality at a solution of the optimality system
(`KKT.dualObjective_eq`), and a version in which the dual point is built from `φ` alone
(`adjoint_gap_certificate`).
-/

noncomputable section

namespace DDG

open Matrix

namespace Site

/-- The convex conjugate of the per-vertex cost, `sup_{x allowed} (p·x − λ·cost(x))`, on its
domain `DualFeasible`:

* penalized vertex with threshold `t = λ·wR` and bounds `[lo, hi]`:
  `hi·max(0, p − t) + lo·min(0, p + t)` (a term with an infinite bound is absent);
* fixed vertex `μ = a`: `a·p`. -/
def conj (lam : ℝ) : Site → ℝ → ℝ
  | penalized wR lo hi, p => hi.toReal * max 0 (p - lam * wR) + lo.toReal * min 0 (p + lam * wR)
  | fixed a, p => a * p

/-- The domain of the conjugate: where a bound is infinite, `p` must respect the threshold on
that side (`p ≤ λ·wR` if `hi = ⊤`, `−λ·wR ≤ p` if `lo = ⊥`). For an unbounded vertex this is
`|p| ≤ λ·wR`; for a vertex with two finite bounds, and for a fixed vertex, there is no
condition. -/
def DualFeasible (lam : ℝ) : Site → ℝ → Prop
  | penalized wR lo hi, p => (hi = ⊤ → p ≤ lam * wR) ∧ (lo = ⊥ → -(lam * wR) ≤ p)
  | fixed _, _ => True

theorem dualFeasible_plain_iff (lam p : ℝ) : plain.DualFeasible lam p ↔ |p| ≤ lam := by
  simp only [plain, DualFeasible, mul_one, forall_const]
  rw [abs_le]; tauto

theorem conj_plain (lam p : ℝ) : plain.conj lam p = 0 := by
  simp [plain, conj]

/-- **Fenchel–Young inequality** for the per-vertex cost:
`p·x ≤ λ·cost(x) + conj(p)` for every allowed `x` and every `p` in the domain of the conjugate. -/
theorem fenchel_young {lam : ℝ} (hlam : 0 ≤ lam) {s : Site} (hs : s.WellFormed) {p x : ℝ}
    (hp : s.DualFeasible lam p) (hx : x ∈ s.feasible) : p * x ≤ lam * s.cost x + s.conj lam p := by
  cases s with
  | penalized wR lo hi =>
    obtain ⟨hw, hlo, hhi⟩ := hs
    obtain ⟨hp1, hp2⟩ := hp
    obtain ⟨hx1, hx2⟩ := hx
    simp only [cost, conj]
    have ht : 0 ≤ lam * wR := mul_nonneg hlam hw
    set t := lam * wR
    have hhi' : 0 ≤ hi.toReal := EReal.toReal_nonneg hhi
    have hlo' : lo.toReal ≤ 0 := EReal.toReal_nonpos hlo
    have hA : 0 ≤ hi.toReal * max 0 (p - t) := mul_nonneg hhi' (le_max_left _ _)
    have hB : 0 ≤ lo.toReal * min 0 (p + t) :=
      mul_nonneg_of_nonpos_of_nonpos hlo' (min_le_left _ _)
    rcases le_or_gt 0 x with h0 | h0
    · rw [abs_of_nonneg h0]
      rcases le_or_gt p t with hpt | hpt
      · nlinarith [mul_nonneg (sub_nonneg.2 hpt) h0]
      · have hne : hi ≠ ⊤ := fun h => absurd (hp1 h) (not_le.mpr hpt)
        have hne' : hi ≠ ⊥ := fun h => by rw [h] at hhi; simp at hhi
        have hxle : x ≤ hi.toReal := by
          have := hx2; rw [← EReal.coe_toReal hne hne'] at this; exact_mod_cast this
        rw [max_eq_right (by linarith : 0 ≤ p - t)]
        nlinarith [mul_nonneg (sub_nonneg.2 hpt.le) (sub_nonneg.2 hxle)]
    · rw [abs_of_neg h0]
      rcases le_or_gt (-t) p with hpt | hpt
      · nlinarith [mul_nonneg (by linarith : 0 ≤ p + t) (by linarith : 0 ≤ -x)]
      · have hne : lo ≠ ⊥ := fun h => absurd (hp2 h) (not_le.mpr hpt)
        have hne' : lo ≠ ⊤ := fun h => by rw [h] at hlo; simp at hlo
        have hxle : lo.toReal ≤ x := by
          have := hx1; rw [← EReal.coe_toReal hne' hne] at this; exact_mod_cast this
        rw [min_eq_right (by linarith : p + t ≤ 0)]
        nlinarith [mul_nonneg (by linarith : 0 ≤ -(p + t)) (sub_nonneg.2 hxle)]
  | fixed a =>
    simp only [feasible, Set.mem_singleton_iff] at hx
    simp [cost, conj, hx, mul_comm]

/-- The conjugate value is attained: some allowed `y` has `p·y − λ·cost(y) = conj(p)`.
(So `conj` is the conjugate, not merely an upper bound.) -/
theorem exists_conj_eq {lam : ℝ} (hlam : 0 ≤ lam) {s : Site} (hs : s.WellFormed) {p : ℝ}
    (hp : s.DualFeasible lam p) : ∃ y ∈ s.feasible, p * y - lam * s.cost y = s.conj lam p := by
  cases s with
  | penalized wR lo hi =>
    obtain ⟨hw, hlo, hhi⟩ := hs
    obtain ⟨hp1, hp2⟩ := hp
    simp only [cost, conj, feasible]
    have ht : 0 ≤ lam * wR := mul_nonneg hlam hw
    set t := lam * wR
    have hhi' : 0 ≤ hi.toReal := EReal.toReal_nonneg hhi
    have hlo' : lo.toReal ≤ 0 := EReal.toReal_nonpos hlo
    rcases lt_or_ge t p with hpt | hpt
    · have hne : hi ≠ ⊤ := fun h => absurd (hp1 h) (not_le.mpr hpt)
      have hne' : hi ≠ ⊥ := fun h => by rw [h] at hhi; simp at hhi
      refine ⟨hi.toReal, ⟨?_, ?_⟩, ?_⟩
      · rw [EReal.coe_toReal hne hne']; exact le_trans hlo hhi
      · rw [EReal.coe_toReal hne hne']
      · rw [abs_of_nonneg hhi', max_eq_right (by linarith), min_eq_left (by linarith)]; ring
    · rcases lt_or_ge p (-t) with hpt' | hpt'
      · have hne : lo ≠ ⊥ := fun h => absurd (hp2 h) (not_le.mpr hpt')
        have hne' : lo ≠ ⊤ := fun h => by rw [h] at hlo; simp at hlo
        refine ⟨lo.toReal, ⟨?_, ?_⟩, ?_⟩
        · rw [EReal.coe_toReal hne' hne]
        · rw [EReal.coe_toReal hne' hne]; exact le_trans hlo hhi
        · rw [abs_of_nonpos hlo', max_eq_left (by linarith), min_eq_right (by linarith)]; ring
      · refine ⟨0, ⟨by simpa using hlo, by simpa using hhi⟩, ?_⟩
        rw [max_eq_left (by linarith), min_eq_left (by linarith)]; simp
  | fixed a =>
    exact ⟨a, rfl, by simp [cost, conj, mul_comm]⟩

/-- A stationary pair is dual feasible. -/
theorem Stationary.dualFeasible {lam : ℝ} (hlam : 0 ≤ lam) {s : Site} (hs : s.WellFormed)
    {p x : ℝ} (h : s.Stationary lam p x) : s.DualFeasible lam p := by
  cases s with
  | penalized wR lo hi =>
    obtain ⟨hw, hlo, hhi⟩ := hs
    have ht : 0 ≤ lam * wR := mul_nonneg hlam hw
    have hf : (fun x => lam * (penalized wR lo hi).cost x) = fun x => lam * wR * |x| := by
      funext x; simp [cost, mul_assoc]
    unfold Stationary at h
    rw [hf, feasible, isSubgradientWithinAt_abs_iff ht hlo hhi] at h
    obtain ⟨_, hR, hL⟩ := h
    constructor
    · intro htop
      have := hR (by rw [htop]; exact EReal.coe_lt_top x)
      split_ifs at this <;> linarith
    · intro hbot
      have := hL (by rw [hbot]; exact EReal.bot_lt_coe x)
      split_ifs at this <;> linarith
  | fixed a => trivial

/-- **Equality in Fenchel–Young** at a stationary pair: `p·x = λ·cost(x) + conj(p)`. -/
theorem Stationary.fenchel_young_eq {lam : ℝ} (hlam : 0 ≤ lam) {s : Site} (hs : s.WellFormed)
    {p x : ℝ} (h : s.Stationary lam p x) : p * x = lam * s.cost x + s.conj lam p := by
  have hd := h.dualFeasible hlam hs
  have h1 := fenchel_young hlam hs hd h.1
  obtain ⟨y, hy, hyeq⟩ := exists_conj_eq hlam hs hd
  have h2 := h.2 y hy
  simp only at h2
  linarith

end Site

namespace ConeProblem

variable {ι : Type*} [Fintype ι] (P : ConeProblem ι)

/-- The dual objective `D(ũ, φ) = −½ Σ wE m ũ² + Σ φ_i Ω_i − Σ conj_i(φ_i)`. -/
def dualObjective (ũ φ : ι → ℝ) : ℝ :=
  -(1 / 2 * ∑ i, P.wE i * P.m i * ũ i ^ 2) + ∑ i, φ i * P.Ω i - ∑ i, (P.site i).conj P.lam (φ i)

/-- Dual feasibility: the adjoint equation `Lᵀφ = −wE m ũ` and `φ_i` in the domain of the
conjugate at every vertex (`|φ_i| ≤ λ·wR_i` at an unbounded vertex). -/
def DualFeasible (ũ φ : ι → ℝ) : Prop :=
  (P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * ũ i)) ∧ ∀ i, (P.site i).DualFeasible P.lam (φ i)

/-- The dual objective as a function of `φ` alone, available when all `wE_i m_i > 0`
(`ũ = −(Lᵀφ)/(wE m)` eliminated): `−½ Σ (Lᵀφ)_i²/(wE_i m_i) + Σ φ_i Ω_i − Σ conj_i(φ_i)`. -/
def dualObjectiveOfAdjoint (φ : ι → ℝ) : ℝ :=
  -(1 / 2 * ∑ i, (P.Lᵀ *ᵥ φ) i ^ 2 / (P.wE i * P.m i)) + ∑ i, φ i * P.Ω i -
    ∑ i, (P.site i).conj P.lam (φ i)

variable {P} {u ũ φ μ v ν : ι → ℝ}

/-- **Gap identity.** If `L u = μ − Ω` and `Lᵀφ = −wE m ũ`, then
`J(u, μ) − D(ũ, φ) = ½ Σ wE_i m_i (u_i − ũ_i)² + Σ_i (λ·cost_i(μ_i) + conj_i(φ_i) − φ_i μ_i)`.
For feasible points every term on the right is nonnegative; this is the cancellation-free way
to evaluate the gap.

[formula: duality-gap-certificate] -/
theorem gap_identity (hstate : P.L *ᵥ u = μ - P.Ω)
    (hadj : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * ũ i)) :
    P.objective u μ - P.dualObjective ũ φ =
      1 / 2 * ∑ i, P.wE i * P.m i * (u i - ũ i) ^ 2 +
        ∑ i, (P.lam * (P.site i).cost (μ i) + (P.site i).conj P.lam (φ i) - φ i * μ i) := by
  have hΩ : ∑ i, φ i * P.Ω i = ∑ i, φ i * μ i + ∑ i, P.wE i * P.m i * ũ i * u i := by
    have e1 : ∑ i, φ i * (μ i - P.Ω i) = φ ⬝ᵥ (P.L *ᵥ u) := by
      rw [hstate]; simp [dotProduct]
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hadj] at e1
    simp only [dotProduct, neg_mul, Finset.sum_neg_distrib, mul_sub, Finset.sum_sub_distrib] at e1
    linarith
  simp only [objective, dualObjective, distortion, coneCost]
  rw [hΩ]
  simp only [Finset.mul_sum, Finset.sum_add_distrib, Finset.sum_sub_distrib]
  have : ∑ i, 1 / 2 * (P.wE i * P.m i * (u i - ũ i) ^ 2) =
      ∑ i, 1 / 2 * (P.wE i * P.m i * u i ^ 2) + ∑ i, 1 / 2 * (P.wE i * P.m i * ũ i ^ 2) -
        ∑ i, P.wE i * P.m i * ũ i * u i := by
    rw [← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun i _ => by ring
  rw [this]; ring

/-- **Weak duality.** Every dual-feasible value is a lower bound for every primal-feasible
value: `D(ũ, φ) ≤ J(v, ν)`.

[formula: duality-gap-certificate] -/
theorem weak_duality (hP : P.WellFormed) (hv : P.Feasible v ν) (hd : P.DualFeasible ũ φ) :
    P.dualObjective ũ φ ≤ P.objective v ν := by
  have h := gap_identity hv.1 hd.1
  have h1 : 0 ≤ ∑ i, P.wE i * P.m i * (v i - ũ i) ^ 2 :=
    Finset.sum_nonneg fun i _ =>
      mul_nonneg (mul_nonneg (hP.wE_nonneg i) (hP.m_nonneg i)) (sq_nonneg _)
  have h2 : 0 ≤ ∑ i, (P.lam * (P.site i).cost (ν i) + (P.site i).conj P.lam (φ i) - φ i * ν i) :=
    Finset.sum_nonneg fun i _ => by
      have := Site.fenchel_young hP.lam_nonneg (hP.site i) (hd.2 i) (hv.2 i)
      linarith
  linarith

/-- **Certificate soundness.** Given a primal-feasible point `(u, μ)` and a dual-feasible point
`(ũ, φ)` of a well-formed problem, the duality gap `J(u, μ) − D(ũ, φ)` is nonnegative and bounds
the suboptimality of `(u, μ)`: for every number `gap` with `J(u, μ) − D(ũ, φ) ≤ gap`, the objective
of `(u, μ)` is within `gap` of the objective of every feasible point, hence of the optimum.

[formula: duality-gap-certificate] -/
theorem duality_gap_certificate (hP : P.WellFormed) (hprimal : P.Feasible u μ)
    (hdual : P.DualFeasible ũ φ) {gap : ℝ} (hgap : P.objective u μ - P.dualObjective ũ φ ≤ gap) :
    0 ≤ gap ∧ ∀ v ν, P.Feasible v ν → P.objective u μ ≤ P.objective v ν + gap := by
  refine ⟨?_, fun v ν hv => ?_⟩
  · have := weak_duality hP hprimal hdual; linarith
  · have := weak_duality hP hv hdual; linarith

/-- A solution of the optimality system is dual feasible (with `ũ = u`). -/
theorem KKT.dualFeasible (hP : P.WellFormed) (h : P.KKT u φ μ) : P.DualFeasible u φ :=
  ⟨h.adjoint, fun i => (h.inclusion i).dualFeasible hP.lam_nonneg (hP.site i)⟩

/-- **Strong duality at a solution of the optimality system**: the dual value at `(u, φ)` equals
the primal value at `(u, μ)`; the gap of the certificate is zero at the optimum, so the
certificate is tight.

[formula: duality-gap-certificate] -/
theorem KKT.dualObjective_eq (hP : P.WellFormed) (h : P.KKT u φ μ) :
    P.dualObjective u φ = P.objective u μ := by
  have hg := gap_identity h.state h.adjoint
  have h1 : ∑ i, P.wE i * P.m i * (u i - u i) ^ 2 = 0 := by simp
  have h2 : ∑ i, (P.lam * (P.site i).cost (μ i) + (P.site i).conj P.lam (φ i) - φ i * μ i) = 0 :=
    Finset.sum_eq_zero fun i _ => by
      have := (h.inclusion i).fenchel_young_eq hP.lam_nonneg (hP.site i)
      linarith
  rw [h1, h2] at hg
  linarith

/-- When all `wE_i m_i > 0`, any `φ` becomes a dual point by *defining*
`ũ_i = −(Lᵀφ)_i/(wE_i m_i)`; the adjoint equation then holds exactly and the dual objective is
`dualObjectiveOfAdjoint φ`. -/
theorem dualObjective_of_adjoint (hpos : ∀ i, 0 < P.wE i * P.m i) (φ : ι → ℝ) :
    (P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * (-(P.Lᵀ *ᵥ φ) i / (P.wE i * P.m i)))) ∧
      P.dualObjective (fun i => -(P.Lᵀ *ᵥ φ) i / (P.wE i * P.m i)) φ =
        P.dualObjectiveOfAdjoint φ := by
  constructor
  · ext i
    have hw := (hpos i).ne'
    generalize P.wE i * P.m i = w at hw
    generalize (P.Lᵀ *ᵥ φ) i = a
    field_simp
  · have e : ∑ i, P.wE i * P.m i * (-(P.Lᵀ *ᵥ φ) i / (P.wE i * P.m i)) ^ 2 =
        ∑ i, (P.Lᵀ *ᵥ φ) i ^ 2 / (P.wE i * P.m i) := by
      refine Finset.sum_congr rfl fun i _ => ?_
      have hw := (hpos i).ne'
      generalize P.wE i * P.m i = w at hw
      generalize (P.Lᵀ *ᵥ φ) i = a
      field_simp
    simp only [dualObjective, dualObjectiveOfAdjoint, e]

/-- **Certificate from the adjoint state alone** (all `wE_i m_i > 0`). Let `(u, μ)` be primal
feasible and let `φ` be *any* vector that lies in the domain of the conjugate at every vertex
(`|φ_i| ≤ λ·wR_i` at an unbounded vertex). Then
`gap := J(u, μ) − dualObjectiveOfAdjoint φ` is nonnegative and `J(u, μ)` is within `gap` of the
optimum. No linear solve is involved in evaluating the bound: only `Lᵀφ`.

[formula: duality-gap-certificate] -/
theorem adjoint_gap_certificate (hP : P.WellFormed) (hpos : ∀ i, 0 < P.wE i * P.m i)
    (hprimal : P.Feasible u μ) (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i)) {gap : ℝ}
    (hgap : P.objective u μ - P.dualObjectiveOfAdjoint φ ≤ gap) :
    0 ≤ gap ∧ ∀ v ν, P.Feasible v ν → P.objective u μ ≤ P.objective v ν + gap := by
  obtain ⟨hadj, hobj⟩ := dualObjective_of_adjoint hpos φ
  rw [← hobj] at hgap
  exact duality_gap_certificate hP hprimal ⟨hadj, hφ⟩ hgap

/-- The lower bound alone: for any `φ` in the domain of the conjugates,
`dualObjectiveOfAdjoint φ` is a lower bound for the optimal value (all `wE_i m_i > 0`). -/
theorem dualObjectiveOfAdjoint_le (hP : P.WellFormed) (hpos : ∀ i, 0 < P.wE i * P.m i)
    (hφ : ∀ i, (P.site i).DualFeasible P.lam (φ i)) (hv : P.Feasible v ν) :
    P.dualObjectiveOfAdjoint φ ≤ P.objective v ν := by
  obtain ⟨hadj, hobj⟩ := dualObjective_of_adjoint hpos φ
  rw [← hobj]
  exact weak_duality hP hv ⟨hadj, hφ⟩

end ConeProblem

end DDG
