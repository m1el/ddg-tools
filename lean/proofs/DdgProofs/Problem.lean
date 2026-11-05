import DdgProofs.SoftThreshold

/-!
# The canonical discrete problem

This file transcribes the "Canonical discrete problem" of the repository's `CLAUDE.md`:

```
minimize    ½ Σ_i wE_i m_i u_i²  +  λ Σ_i wR_i |μ_i|
subject to  (L u)_i = μ_i − Ω_i  (i ∈ I),   u_i = 0  (i ∈ B),   lo_i ≤ μ_i ≤ hi_i
```

The index type `ι` is the set of **unknowns**: the interior vertices `I` on a mesh with boundary
(`L` is then the interior block `L_II`; the constraint `u_i = 0` on `B` has been eliminated, see
`DDG.Cotan` for the block reduction), or all vertices on a closed mesh (`L` is the whole
cotan Laplacian; the additive constant of `u` is then a free variable of the minimization).

Sign conventions are those of `CLAUDE.md`: `L` is positive semidefinite and `u` is the true log
scale factor, so the state equation is `L u = μ − Ω`. Paper Eq. 17 prints `L u = Ω − μ`; with the
paper's own PSD `L` its `u` is the negative of ours.
-/

noncomputable section

namespace DDG

open Matrix

/-- What the problem allows at one vertex.

* `penalized wR lo hi`: the cone angle `μ_i` may be anything in `[lo, hi]` and costs `λ·wR·|μ_i|`.
  `wR = 0` is a free cone site; `lo = ⊥`, `hi = ⊤` mean "unbounded".
* `fixed a`: `μ_i = a` is forced, at no cost. `fixed 0` is a forbidden vertex (`wR = +∞` in
  `CLAUDE.md`); `fixed a` with `a ≠ 0` is a prescribed cone. -/
inductive Site where
  | penalized (wR : ℝ) (lo hi : EReal)
  | fixed (a : ℝ)

namespace Site

/-- The default vertex: unit cone cost, no bounds. -/
def plain : Site := penalized 1 ⊥ ⊤
/-- A free cone site (`wR = 0`): any angle, no cost. -/
def free : Site := penalized 0 ⊥ ⊤
/-- A forbidden vertex (`wR = +∞`): no cone allowed. -/
def forbidden : Site := fixed 0
/-- A prescribed cone of angle `a`. -/
def prescribed (a : ℝ) : Site := fixed a

/-- The set of cone angles allowed at the vertex. -/
def feasible : Site → Set ℝ
  | penalized _ lo hi => boundSet lo hi
  | fixed a => {a}

/-- The cone cost at the vertex, without the factor `λ`: `wR·|μ|`. -/
def cost : Site → ℝ → ℝ
  | penalized wR _ _, x => wR * |x|
  | fixed _, _ => 0

/-- Admissible data: `wR ≥ 0` and `lo ≤ 0 ≤ hi`. -/
def WellFormed : Site → Prop
  | penalized wR lo hi => 0 ≤ wR ∧ lo ≤ 0 ∧ 0 ≤ hi
  | fixed _ => True

/-- The per-vertex optimality inclusion `φ ∈ ∂(λ·cost + ι_feasible)(μ)`: `μ` is allowed and
`λ·cost(μ) + φ·(y − μ) ≤ λ·cost(y)` for every allowed `y`. -/
def Stationary (lam : ℝ) (s : Site) (φ μ : ℝ) : Prop :=
  IsSubgradientWithinAt (fun x => lam * s.cost x) s.feasible φ μ

/-- The per-vertex proximal map: `clamp_[lo,hi] ∘ P_{λ wR}` for a penalized vertex, the constant
`a` for a fixed one. -/
noncomputable def prox (lam : ℝ) : Site → ℝ → ℝ
  | penalized wR lo hi, x => clampE lo hi (softThreshold (lam * wR) x)
  | fixed a, _ => a

theorem cost_nonneg {s : Site} (hs : s.WellFormed) (x : ℝ) : 0 ≤ s.cost x := by
  cases s with
  | penalized wR lo hi => exact mul_nonneg hs.1 (abs_nonneg x)
  | fixed a => exact le_rfl

/-- The allowed set is convex (an interval or a point). -/
theorem feasible_convex (s : Site) {a b : ℝ} (ha : a ∈ s.feasible) (hb : b ∈ s.feasible)
    {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ 1) : (1 - θ) * a + θ * b ∈ s.feasible := by
  cases s with
  | penalized wR lo hi =>
    have hmin : min a b ≤ (1 - θ) * a + θ * b := by
      rcases le_total a b with h | h
      · rw [min_eq_left h]; nlinarith
      · rw [min_eq_right h]; nlinarith
    have hmax : (1 - θ) * a + θ * b ≤ max a b := by
      rcases le_total a b with h | h
      · rw [max_eq_right h]; nlinarith
      · rw [max_eq_left h]; nlinarith
    constructor
    · refine le_trans ?_ (EReal.coe_le_coe_iff.mpr hmin)
      rcases le_total a b with h | h
      · rw [min_eq_left h]; exact ha.1
      · rw [min_eq_right h]; exact hb.1
    · refine le_trans (EReal.coe_le_coe_iff.mpr hmax) ?_
      rcases le_total a b with h | h
      · rw [max_eq_right h]; exact hb.2
      · rw [max_eq_left h]; exact ha.2
  | fixed c =>
    simp only [feasible, Set.mem_singleton_iff] at ha hb ⊢
    rw [ha, hb]; ring

/-- The cone cost is convex. -/
theorem cost_convex {s : Site} (hs : s.WellFormed) (a b : ℝ) {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ 1) :
    s.cost ((1 - θ) * a + θ * b) ≤ (1 - θ) * s.cost a + θ * s.cost b := by
  cases s with
  | penalized wR lo hi =>
    simp only [cost]
    have h : |(1 - θ) * a + θ * b| ≤ (1 - θ) * |a| + θ * |b| := by
      calc |(1 - θ) * a + θ * b| ≤ |(1 - θ) * a| + |θ * b| := abs_add_le _ _
        _ = (1 - θ) * |a| + θ * |b| := by
          rw [abs_mul, abs_mul, abs_of_nonneg h0, abs_of_nonneg (by linarith : 0 ≤ 1 - θ)]
    have := mul_le_mul_of_nonneg_left h hs.1
    linarith
  | fixed c => simp [cost]

/-- **Per-vertex inclusion ⟺ fixed point of the proximal map**, for every kind of vertex:
weighted, bounded, free (`wR = 0`), forbidden and prescribed.
`φ ∈ ∂(λ·cost + ι_feasible)(μ)  ⟺  μ = prox(φ + μ)`.
For the default vertex this is paper Eq. 19, `μ = P_λ(φ + μ)`.

[formula: cone-fixed-point-bounded] -/
theorem stationary_iff_eq_prox {lam : ℝ} (hlam : 0 ≤ lam) {s : Site} (hs : s.WellFormed)
    (φ μ : ℝ) : s.Stationary lam φ μ ↔ μ = s.prox lam (φ + μ) := by
  cases s with
  | penalized wR lo hi =>
    have h := isSubgradientWithinAt_abs_iff_eq_clampE (mul_nonneg hlam hs.1) hs.2.1 hs.2.2 φ μ
    have hf : (fun x => lam * (penalized wR lo hi).cost x) = fun x => lam * wR * |x| := by
      funext x; simp [cost, mul_assoc]
    unfold Stationary
    rw [hf]
    exact h
  | fixed a =>
    simp only [Stationary, IsSubgradientWithinAt, feasible, cost, prox, Set.mem_singleton_iff]
    constructor
    · rintro ⟨h, _⟩; exact h
    · intro h; exact ⟨h, fun y hy => by rw [hy, h]; simp⟩

/-- At a fixed vertex (forbidden or prescribed) the inclusion says only `μ = a`; the adjoint
value `φ` is unconstrained there. -/
theorem stationary_fixed_iff (lam a φ μ : ℝ) : (fixed a).Stationary lam φ μ ↔ μ = a := by
  simp only [Stationary, IsSubgradientWithinAt, feasible, cost, Set.mem_singleton_iff]
  constructor
  · rintro ⟨h, _⟩; exact h
  · intro h; exact ⟨h, fun y hy => by rw [hy, h]; simp⟩

/-- At the default vertex the inclusion is paper Eq. 19: `μ = P_λ(φ + μ)`. -/
theorem stationary_plain_iff {lam : ℝ} (hlam : 0 ≤ lam) (φ μ : ℝ) :
    plain.Stationary lam φ μ ↔ μ = softThreshold lam (φ + μ) := by
  rw [stationary_iff_eq_prox hlam (s := plain) ⟨zero_le_one, bot_le, le_top⟩]
  simp [plain, prox, clampE_bot_top]

/-- At a free cone site the inclusion is `φ = 0`. -/
theorem stationary_free_iff {lam : ℝ} (hlam : 0 ≤ lam) (φ μ : ℝ) :
    free.Stationary lam φ μ ↔ φ = 0 := by
  rw [stationary_iff_eq_prox hlam (s := free) ⟨le_rfl, bot_le, le_top⟩]
  simp only [free, prox, clampE_bot_top, mul_zero]
  exact eq_softThreshold_zero_add_iff φ μ

end Site

/-- Data of the canonical discrete problem on the unknowns `ι`.

* `L`: the cotan Laplacian restricted to the unknowns (interior block `L_II`, or all of `L` on a
  closed mesh);
* `Ω`: angle defects; `m`: lumped vertex masses; `wE`: distortion weights; `lam`: the parameter `λ`;
* `site i`: cone-cost weight and angle bounds at vertex `i`, or a fixed (forbidden / prescribed)
  cone angle. -/
structure ConeProblem (ι : Type*) where
  L : Matrix ι ι ℝ
  Ω : ι → ℝ
  m : ι → ℝ
  wE : ι → ℝ
  lam : ℝ
  site : ι → Site

namespace ConeProblem

variable {ι : Type*} [Fintype ι] (P : ConeProblem ι)

/-- The distortion term `½ Σ_i wE_i m_i u_i²`. -/
def distortion (u : ι → ℝ) : ℝ := 1 / 2 * ∑ i, P.wE i * P.m i * u i ^ 2

/-- The weighted cone norm `Σ_i wR_i |μ_i|` (fixed vertices contribute nothing). -/
def coneCost (μ : ι → ℝ) : ℝ := ∑ i, (P.site i).cost (μ i)

/-- The objective `½ Σ_i wE_i m_i u_i² + λ Σ_i wR_i |μ_i|`. -/
def objective (u μ : ι → ℝ) : ℝ := P.distortion u + P.lam * P.coneCost μ

/-- `(u, μ)` satisfies the constraints: `L u = μ − Ω` and every `μ_i` is allowed at its vertex
(inside its bounds; zero at a forbidden vertex; the given angle at a prescribed cone). -/
def Feasible (u μ : ι → ℝ) : Prop :=
  P.L *ᵥ u = μ - P.Ω ∧ ∀ i, μ i ∈ (P.site i).feasible

/-- `(u, μ)` is a global minimizer: feasible, and no feasible pair has a smaller objective.
On a closed mesh the competitors `v` range over *all* solutions of `L v = ν − Ω`, so the
additive constant is part of the minimization. -/
def IsMinimizer (u μ : ι → ℝ) : Prop :=
  P.Feasible u μ ∧ ∀ v ν, P.Feasible v ν → P.objective u μ ≤ P.objective v ν

/-- Admissible data: `wE ≥ 0`, `m ≥ 0`, `λ ≥ 0`, `wR ≥ 0`, `lo ≤ 0 ≤ hi`. -/
structure WellFormed : Prop where
  wE_nonneg : ∀ i, 0 ≤ P.wE i
  m_nonneg : ∀ i, 0 ≤ P.m i
  lam_nonneg : 0 ≤ P.lam
  site : ∀ i, (P.site i).WellFormed

/-- **The discrete optimality system** for `(u, φ, μ)`, in the sign convention of `CLAUDE.md`:

* state equation `L u = μ − Ω`;
* adjoint equation `Lᵀ φ = −(wE ∘ m ∘ u)`, i.e. `(Lᵀφ)_i = −wE_i m_i u_i`;
* per-vertex inclusion `φ_i ∈ ∂(λ·cost_i + ι_feasible_i)(μ_i)`.

Paper Eq. 17 / System B.2 reads `L u = Ω − μ`, `Lᵀφ = M u`, `μ ∈ ∂1_λ(φ)`; substituting the
paper's `u = −u_ours` gives exactly this system with the same `φ` and `μ`, and with the weights
`wE` of the paper's generalized system (App. A.3) included. -/
structure KKT (u φ μ : ι → ℝ) : Prop where
  state : P.L *ᵥ u = μ - P.Ω
  adjoint : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i)
  inclusion : ∀ i, (P.site i).Stationary P.lam (φ i) (μ i)

/-- The inclusions of a KKT triple contain feasibility of `μ`. -/
theorem KKT.feasible {P : ConeProblem ι} {u φ μ : ι → ℝ} (h : P.KKT u φ μ) : P.Feasible u μ :=
  ⟨h.state, fun i => (h.inclusion i).1⟩

/-- The optimality system in fixed-point form (the form the solver's residual uses):
`L u = μ − Ω`, `Lᵀφ = −wE m u`, `μ_i = prox_i(φ_i + μ_i)`. -/
theorem kkt_iff_fixedPoint (hP : P.WellFormed) (u φ μ : ι → ℝ) :
    P.KKT u φ μ ↔
      P.L *ᵥ u = μ - P.Ω ∧ (P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i)) ∧
        ∀ i, μ i = (P.site i).prox P.lam (φ i + μ i) := by
  constructor
  · intro h
    exact ⟨h.state, h.adjoint, fun i =>
      (Site.stationary_iff_eq_prox hP.lam_nonneg (hP.site i) _ _).mp (h.inclusion i)⟩
  · rintro ⟨h1, h2, h3⟩
    exact ⟨h1, h2, fun i =>
      (Site.stationary_iff_eq_prox hP.lam_nonneg (hP.site i) _ _).mpr (h3 i)⟩

end ConeProblem

end DDG
