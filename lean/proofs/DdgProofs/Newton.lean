import DdgProofs.Optimality
import DDG.Optimization.PiecewiseAffine

/-!
# The semismooth Newton method and its stopping rule (paper Sec. 5.3)

The residual maps of the two stages are piecewise affine: the only nonlinearity is a per-vertex
function that is affine on each of finitely many "pieces" (the active sets of paper Sec. 5.3.1).
A semismooth Newton step solves the linearization on the current pieces. Hence:

> if a Newton step does not change the active sets, the new iterate is an **exact** zero of the
> residual.

This is the convergence criterion of paper Sec. 5.3.3. It is proved here

* abstractly, for any per-vertex piecewise-affine function
  (`newton_stageI_exact`, `newton_stageII_exact`);
* in the paper's form, with the soft threshold and the active sets `A⁺`, `A⁻` of Eq. 24
  (`paper_stageI_stopping_rule`, `paper_stageII_stopping_rule`);
* for the canonical problem with all vertex kinds (bounds, free, forbidden, prescribed), where
  the conclusion is upgraded to "the new iterate is a global minimizer"
  (`ConeProblem.isMinimizer_of_newton_step`, `ConeProblem.isRegMinimizer_of_newton_step`).

Sign convention: `CLAUDE.md` (`L u = μ − Ω`). The residuals below are the paper's Eqs. 20 and 22
after the substitution `u_paper = −u`.
-/

noncomputable section

namespace DDG

open Matrix

/-! ## The paper's form: soft threshold and the active sets `A⁺`, `A⁻` -/

/-- Diagonal entry of the active-set matrix `A⁺` of paper Eq. 24: `1` if `x > t`, else `0`. -/
def activePos (t x : ℝ) : ℝ := if t < x then 1 else 0

/-- Diagonal entry of the active-set matrix `A⁻` of paper Eq. 24: `1` if `x < −t`, else `0`. -/
def activeNeg (t x : ℝ) : ℝ := if x < -t then 1 else 0

/-- `P_t(x) = A⁺(x)·(x − t) + A⁻(x)·(x + t)`. -/
theorem softThreshold_eq_active {t : ℝ} (ht : 0 ≤ t) (x : ℝ) :
    softThreshold t x = activePos t x * (x - t) + activeNeg t x * (x + t) := by
  unfold activePos activeNeg
  rcases lt_or_ge t x with h | h
  · rw [softThreshold_of_gt ht h, ite_eq_left h, ite_eq_right (by linarith)]; ring
  · rcases lt_or_ge x (-t) with h' | h'
    · rw [softThreshold_of_lt ht h', ite_eq_right (not_lt.mpr h), ite_eq_left h']; ring
    · rw [softThreshold_of_abs_le (abs_le.mpr ⟨h', h⟩), ite_eq_right (not_lt.mpr h),
        ite_eq_right (not_lt.mpr h')]; ring

/-- The soft threshold is affine on each configuration of the active sets, with generalized
derivative `D = A⁺ + A⁻` (paper Sec. 5.3.2, citing Hintermüller et al. 2002, Lemma 3.1).

[formula: soft-threshold-active-sets] -/
theorem softThreshold_isPiecewiseAffine {t : ℝ} (ht : 0 ≤ t) :
    IsPiecewiseAffine (softThreshold t) (fun x => (activePos t x, activeNeg t x))
      (fun x => activePos t x + activeNeg t x) := by
  intro x x' h
  simp only [Prod.mk.injEq] at h
  rw [softThreshold_eq_active ht x', softThreshold_eq_active ht x, h.1, h.2]
  ring

section Paper

variable {ι : Type*} [Fintype ι]

/-- **Stopping rule of paper Sec. 5.3.3, Stage II** (weighted thresholds `t_i = λ·wR_i ≥ 0`,
distortion weights `w_i = wE_i m_i`). With

* residual `F(u, φ, μ) = (L u − μ + Ω,  Lᵀφ + w u,  μ − P_t(φ + μ))` (paper Eq. 20),
* differential `DF(v, q, ν) = (L v − ν,  Lᵀq + w v,  −D q + (I − D) ν)`, `D = A⁺ + A⁻` evaluated at
  `φ + μ` (paper Eq. 26),

if `(v, q, ν)` solves `DF(u, φ, μ)(v, q, ν) = −F(u, φ, μ)` and the active sets `A⁺(φ + μ)`,
`A⁻(φ + μ)` are unchanged by the step, then the new iterate is an exact zero of `F`.

[formula: newton-stopping-rule] -/
theorem paper_stageII_stopping_rule (L : Matrix ι ι ℝ) (w Ω t : ι → ℝ) (ht : ∀ i, 0 ≤ t i)
    (u φ μ v q ν : ι → ℝ)
    (h1 : L *ᵥ v - ν = -(L *ᵥ u - μ + Ω))
    (h2 : Lᵀ *ᵥ q + (fun i => w i * v i) = -(Lᵀ *ᵥ φ + fun i => w i * u i))
    (h3 : ∀ i, -(activePos (t i) (φ i + μ i) + activeNeg (t i) (φ i + μ i)) * q i +
        (1 - (activePos (t i) (φ i + μ i) + activeNeg (t i) (φ i + μ i))) * ν i =
        -(μ i - softThreshold (t i) (φ i + μ i)))
    (hpos : ∀ i, activePos (t i) ((φ i + q i) + (μ i + ν i)) = activePos (t i) (φ i + μ i))
    (hneg : ∀ i, activeNeg (t i) ((φ i + q i) + (μ i + ν i)) = activeNeg (t i) (φ i + μ i)) :
    L *ᵥ (u + v) - (μ + ν) + Ω = 0 ∧
      Lᵀ *ᵥ (φ + q) + (fun i => w i * (u i + v i)) = 0 ∧
      ∀ i, μ i + ν i = softThreshold (t i) ((φ i + q i) + (μ i + ν i)) := by
  have := newton_stageII_exact L w Ω (fun i => softThreshold (t i))
    (fun i x => (activePos (t i) x, activeNeg (t i) x))
    (fun i x => activePos (t i) x + activeNeg (t i) x)
    (fun i => softThreshold_isPiecewiseAffine (ht i)) u φ μ v q ν h1 h2
    (fun i => by have := h3 i; linarith)
    (fun i => by simp only [hpos i, hneg i])
  exact ⟨this.1, this.2.1, fun i => by have := this.2.2 i; linarith⟩

/-- **Stopping rule of paper Sec. 5.3.3, Stage I** (thresholds `t_i ≥ 0`, per-vertex
regularization weights `γ_i > 0`). With

* residual `F_γ(u, φ) = (L u + Ω − P_t(φ)/γ,  Lᵀφ + w u)` (paper Eq. 22),
* differential `DF_γ(v, q) = (L v − (D/γ) q,  Lᵀq + w v)`, `D = A⁺ + A⁻` evaluated at `φ`
  (paper Eq. 25),

if `(v, q)` solves `DF_γ(u, φ)(v, q) = −F_γ(u, φ)` and the active sets `A⁺(φ)`, `A⁻(φ)` are
unchanged by the step, then the new iterate is an exact zero of `F_γ`.

The mass-consistent regularization is `γ_i = γ/m_i`, i.e. `μ_i = (m_i/γ)·P_t(φ_i)`; paper Eq. 22
prints the case `γ_i = γ`. The statement covers both.

[formula: newton-stopping-rule] -/
theorem paper_stageI_stopping_rule (L : Matrix ι ι ℝ) (w Ω t γ : ι → ℝ) (ht : ∀ i, 0 ≤ t i)
    (u φ v q : ι → ℝ)
    (h1 : ∀ i, (L *ᵥ v) i - (activePos (t i) (φ i) + activeNeg (t i) (φ i)) / γ i * q i =
        -((L *ᵥ u) i + Ω i - softThreshold (t i) (φ i) / γ i))
    (h2 : Lᵀ *ᵥ q + (fun i => w i * v i) = -(Lᵀ *ᵥ φ + fun i => w i * u i))
    (hpos : ∀ i, activePos (t i) (φ i + q i) = activePos (t i) (φ i))
    (hneg : ∀ i, activeNeg (t i) (φ i + q i) = activeNeg (t i) (φ i)) :
    (∀ i, (L *ᵥ (u + v)) i + Ω i - softThreshold (t i) (φ i + q i) / γ i = 0) ∧
      Lᵀ *ᵥ (φ + q) + (fun i => w i * (u i + v i)) = 0 := by
  refine newton_stageI_exact L w Ω (fun i x => softThreshold (t i) x / γ i)
    (fun i x => (activePos (t i) x, activeNeg (t i) x))
    (fun i x => (activePos (t i) x + activeNeg (t i) x) / γ i)
    (fun i x x' h => ?_) u φ v q h1 h2 (fun i => by simp only [hpos i, hneg i])
  show softThreshold (t i) x' / γ i = softThreshold (t i) x / γ i +
    (activePos (t i) x + activeNeg (t i) x) / γ i * (x' - x)
  rw [softThreshold_isPiecewiseAffine (ht i) x x' h]
  ring

end Paper

/-! ## All vertex kinds: pieces of the (regularized) proximal map -/

/-- The five pieces of `x ↦ clamp_[lo,hi](P_t(x)/γ)`: clamped at the lower bound, negative
active (`A⁻`), inactive, positive active (`A⁺`), clamped at the upper bound. -/
inductive Piece
  | lower | neg | zero | pos | upper
  deriving DecidableEq

namespace Site

/-- The regularized per-vertex proximal map `x ↦ clamp_[lo,hi](P_{λ wR}(x)/γ)` (a constant at a
fixed vertex). `γ = 1` is `Site.prox`. -/
def proxReg (lam γ : ℝ) : Site → ℝ → ℝ
  | penalized wR lo hi, x => clampE lo hi (softThreshold (lam * wR) x / γ)
  | fixed a, _ => a

theorem prox_eq_proxReg (lam : ℝ) (s : Site) (x : ℝ) : s.prox lam x = s.proxReg lam 1 x := by
  cases s <;> simp [prox, proxReg]

open Classical in
/-- The active set (piece) of `x` for `proxReg`. For an unbounded vertex only `neg`, `zero`,
`pos` occur and they are the paper's `A⁻`, inactive, `A⁺`. -/
def pieceReg (lam γ : ℝ) : Site → ℝ → Piece
  | penalized wR lo hi, x =>
    if hi < (((x - lam * wR) / γ : ℝ) : EReal) then Piece.upper
    else if lam * wR < x then Piece.pos
    else if (((x + lam * wR) / γ : ℝ) : EReal) < lo then Piece.lower
    else if x < -(lam * wR) then Piece.neg
    else Piece.zero
  | fixed _, _ => Piece.zero

/-- Slope of `proxReg` on a piece: `1/γ` on the two active pieces, `0` elsewhere. -/
def pieceSlope (γ : ℝ) : Piece → ℝ
  | Piece.pos => 1 / γ
  | Piece.neg => 1 / γ
  | _ => 0

/-- Intercept of `proxReg` on a piece. -/
def pieceIntercept (lam γ : ℝ) : Site → Piece → ℝ
  | penalized _ _ hi, Piece.upper => hi.toReal
  | penalized wR _ _, Piece.pos => -(lam * wR) / γ
  | penalized _ _ _, Piece.zero => 0
  | penalized wR _ _, Piece.neg => lam * wR / γ
  | penalized _ lo _, Piece.lower => lo.toReal
  | fixed a, _ => a

/-- The generalized derivative of `proxReg` used by the semismooth Newton method. -/
def slopeReg (lam γ : ℝ) (s : Site) (x : ℝ) : ℝ := pieceSlope γ (s.pieceReg lam γ x)

/-- On each piece `proxReg` is the affine function `slope·x + intercept` of that piece. -/
theorem proxReg_eq_affine {lam γ : ℝ} (hlam : 0 ≤ lam) (hγ : 0 < γ) {s : Site}
    (hs : s.WellFormed) (x : ℝ) :
    s.proxReg lam γ x =
      pieceSlope γ (s.pieceReg lam γ x) * x + s.pieceIntercept lam γ (s.pieceReg lam γ x) := by
  cases s with
  | penalized wR lo hi =>
    obtain ⟨hw, hlo, hhi⟩ := hs
    have ht : 0 ≤ lam * wR := mul_nonneg hlam hw
    simp only [proxReg, pieceReg]
    set t := lam * wR
    split_ifs with h1 h2 h3 h4
    · -- clamped at the upper bound
      have hpos : (0 : EReal) < (((x - t) / γ : ℝ) : EReal) := lt_of_le_of_lt hhi h1
      have hpos' : 0 < (x - t) / γ := by exact_mod_cast hpos
      have hx : t < x := by
        by_contra hc
        push Not at hc
        have : (x - t) / γ ≤ 0 := div_nonpos_of_nonpos_of_nonneg (by linarith) hγ.le
        linarith
      rw [softThreshold_of_gt ht hx]
      obtain ⟨r, hr, hc⟩ := clampE_of_hi_le (x := (x - t) / γ) hlo hhi h1.le
      rw [hc]
      subst hr
      simp [pieceSlope, pieceIntercept]
    · -- positive active
      rw [softThreshold_of_gt ht h2]
      have h0 : 0 ≤ (x - t) / γ := div_nonneg (by linarith) hγ.le
      rw [clampE_of_mem (le_trans hlo (by exact_mod_cast h0)) (not_lt.mp h1)]
      simp only [pieceSlope, pieceIntercept]
      ring
    · -- clamped at the lower bound
      have hneg : (((x + t) / γ : ℝ) : EReal) < 0 := lt_of_lt_of_le h3 hlo
      have hneg' : (x + t) / γ < 0 := by exact_mod_cast hneg
      have hx : x < -t := by
        by_contra hc
        push Not at hc
        have : 0 ≤ (x + t) / γ := div_nonneg (by linarith) hγ.le
        linarith
      rw [softThreshold_of_lt ht hx]
      obtain ⟨r, hr, hc⟩ := clampE_of_le_lo (x := (x + t) / γ) hlo hhi h3.le
      rw [hc]
      subst hr
      simp [pieceSlope, pieceIntercept]
    · -- negative active
      rw [softThreshold_of_lt ht h4]
      have h0 : (x + t) / γ ≤ 0 := div_nonpos_of_nonpos_of_nonneg (by linarith) hγ.le
      rw [clampE_of_mem (not_lt.mp h3) (le_trans (by exact_mod_cast h0) hhi)]
      simp only [pieceSlope, pieceIntercept]
      ring
    · -- inactive
      rw [softThreshold_of_abs_le (abs_le.mpr ⟨not_lt.mp h4, not_lt.mp h2⟩), zero_div,
        clampE_of_mem (by simpa using hlo) (by simpa using hhi)]
      simp [pieceSlope, pieceIntercept]
  | fixed a =>
    simp [proxReg, pieceReg, pieceSlope, pieceIntercept]

/-- The (regularized) proximal map of every vertex kind is piecewise affine with respect to its
active sets, with generalized derivative `slopeReg`.

[formula: soft-threshold-active-sets] -/
theorem proxReg_isPiecewiseAffine {lam γ : ℝ} (hlam : 0 ≤ lam) (hγ : 0 < γ) {s : Site}
    (hs : s.WellFormed) :
    IsPiecewiseAffine (s.proxReg lam γ) (s.pieceReg lam γ) (s.slopeReg lam γ) := by
  intro x x' h
  rw [proxReg_eq_affine hlam hγ hs x', proxReg_eq_affine hlam hγ hs x, h, slopeReg]
  ring

/-- For an unbounded vertex the generalized derivative is the paper's `D = A⁺ + A⁻`, divided
by `γ`. -/
theorem slopeReg_unbounded {lam γ wR : ℝ} (ht : 0 ≤ lam * wR) (x : ℝ) :
    (penalized wR ⊥ ⊤).slopeReg lam γ x =
      (activePos (lam * wR) x + activeNeg (lam * wR) x) / γ := by
  simp only [slopeReg, pieceReg, activePos, activeNeg, not_top_lt, not_lt_bot, ite_false]
  by_cases h1 : lam * wR < x
  · have h2 : ¬x < -(lam * wR) := by linarith
    simp [h1, h2, pieceSlope]
  · by_cases h2 : x < -(lam * wR)
    · simp [h1, h2, pieceSlope]
    · simp [h1, h2, pieceSlope]

/-- The regularized per-vertex inclusion in fixed-point form:
`φ − γμ ∈ ∂(λ·cost + ι_feasible)(μ)  ⟺  μ = proxReg(φ)`.

[formula: regularized-cone-update] -/
theorem stationary_sub_iff_eq_proxReg {lam γ : ℝ} (hlam : 0 ≤ lam) (hγ : 0 < γ) {s : Site}
    (hs : s.WellFormed) (φ μ : ℝ) :
    s.Stationary lam (φ - γ * μ) μ ↔ μ = s.proxReg lam γ φ := by
  cases s with
  | penalized wR lo hi =>
    have h := isSubgradientWithinAt_abs_sub_iff_eq_clampE (mul_nonneg hlam hs.1) hγ hs.2.1 hs.2.2
      φ μ
    have hf : (fun x => lam * (penalized wR lo hi).cost x) = fun x => lam * wR * |x| := by
      funext x; simp [cost, mul_assoc]
    unfold Stationary
    rw [hf]
    exact h
  | fixed a => exact stationary_fixed_iff lam a _ μ

end Site

/-! ## The canonical problem: residuals, Newton steps, and optimality of the result -/

namespace ConeProblem

variable {ι : Type*} [Fintype ι] (P : ConeProblem ι)

/-- **Stage II residual** `F(u, φ, μ) = (L u − μ + Ω,  Lᵀφ + wE m u,  μ − prox(φ + μ))`
(paper Eq. 20 in the sign convention of `CLAUDE.md`, with weights, bounds and fixed vertices). -/
def residual (u φ μ : ι → ℝ) : (ι → ℝ) × (ι → ℝ) × (ι → ℝ) :=
  (P.L *ᵥ u - μ + P.Ω,
   P.Lᵀ *ᵥ φ + (fun i => P.wE i * P.m i * u i),
   fun i => μ i - (P.site i).prox P.lam (φ i + μ i))

/-- **Stage II semismooth differential** `DF(u, φ, μ)` applied to `(v, q, ν)`:
`(L v − ν,  Lᵀq + wE m v,  ν − D·(q + ν))` with `D` the generalized derivative of the proximal map
at `φ + μ` (paper Eq. 26: third row `−D q + (I − D) ν`). -/
def residualDeriv (_u φ μ : ι → ℝ) (y : (ι → ℝ) × (ι → ℝ) × (ι → ℝ)) :
    (ι → ℝ) × (ι → ℝ) × (ι → ℝ) :=
  (P.L *ᵥ y.1 - y.2.2,
   P.Lᵀ *ᵥ y.2.1 + (fun i => P.wE i * P.m i * y.1 i),
   fun i => y.2.2 i - (P.site i).slopeReg P.lam 1 (φ i + μ i) * (y.2.1 i + y.2.2 i))

/-- **Stage I residual** `F_γ(u, φ) = (L u + Ω − proxReg_γ(φ),  Lᵀφ + wE m u)` (paper Eq. 22 in
the sign convention of `CLAUDE.md`). `γ i > 0` is the per-vertex regularization weight; the
mass-consistent choice is `γ i = γ / m i`, the paper prints `γ i = γ`. -/
def residualReg (γ : ι → ℝ) (u φ : ι → ℝ) : (ι → ℝ) × (ι → ℝ) :=
  (fun i => (P.L *ᵥ u) i + P.Ω i - (P.site i).proxReg P.lam (γ i) (φ i),
   P.Lᵀ *ᵥ φ + fun i => P.wE i * P.m i * u i)

/-- **Stage I semismooth differential** `DF_γ(u, φ)` applied to `(v, q)`:
`(L v − D_γ q,  Lᵀq + wE m v)` (paper Eq. 25). -/
def residualRegDeriv (γ : ι → ℝ) (_u φ : ι → ℝ) (y : (ι → ℝ) × (ι → ℝ)) :
    (ι → ℝ) × (ι → ℝ) :=
  (fun i => (P.L *ᵥ y.1) i - (P.site i).slopeReg P.lam (γ i) (φ i) * y.2 i,
   P.Lᵀ *ᵥ y.2 + fun i => P.wE i * P.m i * y.1 i)

/-- The objective of the regularized problem (paper Eq. 12, discretized):
`J(u, μ) + Σ_i (γ_i/2) μ_i²`. With `γ_i = γ/m_i` the last term is `(γ/2)‖μ‖²_{L²}` for the
density `μ_i/m_i`. -/
def regObjective (γ : ι → ℝ) (u μ : ι → ℝ) : ℝ := P.objective u μ + ∑ i, γ i / 2 * μ i ^ 2

/-- `(u, μ)` is a global minimizer of the regularized problem. -/
def IsRegMinimizer (γ : ι → ℝ) (u μ : ι → ℝ) : Prop :=
  P.Feasible u μ ∧ ∀ v ν, P.Feasible v ν → P.regObjective γ u μ ≤ P.regObjective γ v ν

variable {P} {u φ μ v q ν : ι → ℝ}

/-- The zeros of the Stage II residual are exactly the solutions of the optimality system. -/
theorem residual_eq_zero_iff (hP : P.WellFormed) : P.residual u φ μ = 0 ↔ P.KKT u φ μ := by
  rw [kkt_iff_fixedPoint P hP]
  simp only [residual, Prod.mk_eq_zero]
  constructor
  · rintro ⟨h1, h2, h3⟩
    refine ⟨?_, ?_, fun i => ?_⟩
    · ext i
      have := congrFun h1 i
      simp only [Pi.add_apply, Pi.sub_apply, Pi.zero_apply] at this ⊢
      linarith
    · ext i
      have := congrFun h2 i
      simp only [Pi.add_apply, Pi.zero_apply] at this ⊢
      linarith
    · have := congrFun h3 i
      simp only [Pi.zero_apply] at this
      linarith
  · rintro ⟨h1, h2, h3⟩
    refine ⟨?_, ?_, ?_⟩
    · rw [h1]; ext i; simp
    · rw [h2]; ext i; simp
    · ext i
      have := h3 i
      simp only [Pi.zero_apply]
      linarith

/-- **Stopping rule, Stage II, canonical problem.** If `(v, q, ν)` solves the Newton system
`DF(u, φ, μ)(v, q, ν) = −F(u, φ, μ)` and the active set of every vertex is unchanged by the step,
then the new iterate is an exact zero of the residual.

[formula: newton-stopping-rule] -/
theorem residual_eq_zero_of_newton_step (hP : P.WellFormed)
    (hstep : P.residualDeriv u φ μ (v, q, ν) = -P.residual u φ μ)
    (hact : ∀ i, (P.site i).pieceReg P.lam 1 ((φ i + q i) + (μ i + ν i)) =
      (P.site i).pieceReg P.lam 1 (φ i + μ i)) :
    P.residual (u + v) (φ + q) (μ + ν) = 0 := by
  simp only [residualDeriv, residual, Prod.neg_mk, Prod.mk.injEq] at hstep
  obtain ⟨h1, h2, h3⟩ := hstep
  have := newton_stageII_exact P.L (fun i => P.wE i * P.m i) P.Ω
    (fun i => (P.site i).proxReg P.lam 1) (fun i => (P.site i).pieceReg P.lam 1)
    (fun i => (P.site i).slopeReg P.lam 1)
    (fun i => Site.proxReg_isPiecewiseAffine hP.lam_nonneg one_pos (hP.site i))
    u φ μ v q ν h1 h2
    (fun i => by
      have := congrFun h3 i
      simp only [Pi.neg_apply, Site.prox_eq_proxReg] at this
      exact this)
    hact
  obtain ⟨g1, g2, g3⟩ := this
  simp only [residual, Prod.mk_eq_zero]
  refine ⟨g1, g2, ?_⟩
  ext i
  simp only [Pi.add_apply, Pi.zero_apply, Site.prox_eq_proxReg]
  exact g3 i

/-- **Stopping rule, Stage II: the result is optimal.** If a semismooth Newton step for the
optimality system leaves every active set unchanged, the new iterate `(u + v, μ + ν)` is a global
minimizer of the canonical problem (with adjoint state `φ + q`).

[formula: newton-stopping-rule] -/
theorem isMinimizer_of_newton_step (hP : P.WellFormed)
    (hstep : P.residualDeriv u φ μ (v, q, ν) = -P.residual u φ μ)
    (hact : ∀ i, (P.site i).pieceReg P.lam 1 ((φ i + q i) + (μ i + ν i)) =
      (P.site i).pieceReg P.lam 1 (φ i + μ i)) :
    P.KKT (u + v) (φ + q) (μ + ν) ∧ P.IsMinimizer (u + v) (μ + ν) := by
  have h := (residual_eq_zero_iff hP).mp (residual_eq_zero_of_newton_step hP hstep hact)
  exact ⟨h, h.isMinimizer hP.wE_nonneg hP.m_nonneg⟩

/-- **Stopping rule, Stage I, canonical problem.** If `(v, q)` solves the Newton system
`DF_γ(u, φ)(v, q) = −F_γ(u, φ)` and the active set of every vertex is unchanged by the step, then
the new iterate is an exact zero of the regularized residual.

[formula: newton-stopping-rule] -/
theorem residualReg_eq_zero_of_newton_step (hP : P.WellFormed) {γ : ι → ℝ}
    (hγ : ∀ i, 0 < γ i)
    (hstep : P.residualRegDeriv γ u φ (v, q) = -P.residualReg γ u φ)
    (hact : ∀ i, (P.site i).pieceReg P.lam (γ i) (φ i + q i) =
      (P.site i).pieceReg P.lam (γ i) (φ i)) :
    P.residualReg γ (u + v) (φ + q) = 0 := by
  simp only [residualRegDeriv, residualReg, Prod.neg_mk, Prod.mk.injEq] at hstep
  obtain ⟨h1, h2⟩ := hstep
  have := newton_stageI_exact P.L (fun i => P.wE i * P.m i) P.Ω
    (fun i => (P.site i).proxReg P.lam (γ i)) (fun i => (P.site i).pieceReg P.lam (γ i))
    (fun i => (P.site i).slopeReg P.lam (γ i))
    (fun i => Site.proxReg_isPiecewiseAffine hP.lam_nonneg (hγ i) (hP.site i))
    u φ v q
    (fun i => by
      have := congrFun h1 i
      simp only [Pi.neg_apply] at this
      exact this)
    h2 hact
  obtain ⟨g1, g2⟩ := this
  simp only [residualReg, Prod.mk_eq_zero]
  refine ⟨?_, g2⟩
  ext i
  simp only [Pi.add_apply, Pi.zero_apply]
  exact g1 i

/-- **Sufficiency for the regularized problem.** If `L u = μ − Ω`, `Lᵀφ = −wE m u` and
`φ_i − γ_i μ_i ∈ ∂(λ·cost_i + ι_feasible_i)(μ_i)` at every vertex (`γ_i ≥ 0`), then `(u, μ)` is a
global minimizer of the regularized objective `J(u, μ) + Σ (γ_i/2) μ_i²`. For `γ = 0` this is the
sufficiency theorem of the unregularized problem.

[formula: regularized-optimality-sufficient] -/
theorem isRegMinimizer_of_stationary (hP : P.WellFormed) {γ : ι → ℝ} (hγ : ∀ i, 0 ≤ γ i)
    (hstate : P.L *ᵥ u = μ - P.Ω) (hadj : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i))
    (hst : ∀ i, (P.site i).Stationary P.lam (φ i - γ i * μ i) (μ i)) :
    P.IsRegMinimizer γ u μ := by
  have hincl : ∀ i, IsSubgradientWithinAt
      (fun x => P.lam * (P.site i).cost x + γ i / 2 * x ^ 2) (P.site i).feasible (φ i) (μ i) := by
    intro i
    refine ⟨(hst i).1, fun y hy => ?_⟩
    have := (hst i).2 y hy
    simp only at this ⊢
    nlinarith [sq_nonneg (y - μ i), hγ i, mul_nonneg (hγ i) (sq_nonneg (y - μ i))]
  refine ⟨⟨hstate, fun i => (hst i).1⟩, fun v ν hv => ?_⟩
  have := quadratic_growth_of_stationary P.L (fun i => P.wE i * P.m i)
    (fun i x => P.lam * (P.site i).cost x + γ i / 2 * x ^ 2) (fun i => (P.site i).feasible)
    P.Ω u φ μ v ν hstate hadj hincl hv.1 hv.2
  have hnn : 0 ≤ ∑ i, 1 / 2 * (P.wE i * P.m i * (v i - u i) ^ 2) :=
    Finset.sum_nonneg fun i _ => mul_nonneg (by norm_num)
      (mul_nonneg (mul_nonneg (hP.wE_nonneg i) (hP.m_nonneg i)) (sq_nonneg _))
  simp only [regObjective, objective, distortion, coneCost, Finset.mul_sum,
    Finset.sum_add_distrib] at this ⊢
  linarith

/-- **A zero of the Stage I residual solves the regularized problem.** If `F_γ(u, φ) = 0` then,
with `μ_i := proxReg_γ(φ_i)` (for a default vertex `μ_i = P_λ(φ_i)/γ_i`), the pair `(u, μ)` is a
global minimizer of the regularized objective `J(u, μ) + Σ (γ_i/2) μ_i²`.

[formula: regularized-optimality-sufficient] -/
theorem isRegMinimizer_of_residualReg_eq_zero (hP : P.WellFormed) {γ : ι → ℝ}
    (hγ : ∀ i, 0 < γ i) (h : P.residualReg γ u φ = 0) :
    P.IsRegMinimizer γ u (fun i => (P.site i).proxReg P.lam (γ i) (φ i)) := by
  simp only [residualReg, Prod.mk_eq_zero] at h
  obtain ⟨h1, h2⟩ := h
  set μ : ι → ℝ := fun i => (P.site i).proxReg P.lam (γ i) (φ i) with hμ
  have hstate : P.L *ᵥ u = μ - P.Ω := by
    ext i
    have := congrFun h1 i
    simp only [Pi.zero_apply, Pi.sub_apply] at this ⊢
    linarith
  have hadj : P.Lᵀ *ᵥ φ = fun i => -(P.wE i * P.m i * u i) := by
    ext i
    have := congrFun h2 i
    simp only [Pi.add_apply, Pi.zero_apply] at this ⊢
    linarith
  exact isRegMinimizer_of_stationary hP (fun i => (hγ i).le) hstate hadj fun i =>
    (Site.stationary_sub_iff_eq_proxReg hP.lam_nonneg (hγ i) (hP.site i) _ _).mpr rfl

/-- **Stopping rule, Stage I: the result is optimal for the regularized problem.**

[formula: newton-stopping-rule] -/
theorem isRegMinimizer_of_newton_step (hP : P.WellFormed) {γ : ι → ℝ} (hγ : ∀ i, 0 < γ i)
    (hstep : P.residualRegDeriv γ u φ (v, q) = -P.residualReg γ u φ)
    (hact : ∀ i, (P.site i).pieceReg P.lam (γ i) (φ i + q i) =
      (P.site i).pieceReg P.lam (γ i) (φ i)) :
    P.IsRegMinimizer γ (u + v) (fun i => (P.site i).proxReg P.lam (γ i) ((φ + q) i)) :=
  isRegMinimizer_of_residualReg_eq_zero hP hγ
    (residualReg_eq_zero_of_newton_step hP hγ hstep hact)

/-- **The regularized solution is nearly optimal for the true problem.** If `(u_γ, μ_γ)`
minimizes the regularized objective and `(u, μ)` is any feasible pair (e.g. the true minimizer),
then `J(u_γ, μ_γ) ≤ J(u, μ) + Σ_i (γ_i/2) μ_i²`.

[formula: regularized-suboptimality] -/
theorem IsRegMinimizer.objective_le {γ : ι → ℝ} (hγ : ∀ i, 0 ≤ γ i) {uγ μγ : ι → ℝ}
    (h : P.IsRegMinimizer γ uγ μγ) (hf : P.Feasible u μ) :
    P.objective uγ μγ ≤ P.objective u μ + ∑ i, γ i / 2 * μ i ^ 2 := by
  have h1 := h.2 u μ hf
  have h2 : 0 ≤ ∑ i, γ i / 2 * μγ i ^ 2 :=
    Finset.sum_nonneg fun i _ => mul_nonneg (by linarith [hγ i]) (sq_nonneg _)
  unfold regObjective at h1
  linarith

end ConeProblem

end DDG
