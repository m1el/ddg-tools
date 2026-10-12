import DdgOracle.Scalar

/-!
# Soft thresholding, its bounded variant, active sets

Per vertex the cone term of the canonical problem is the convex function

  `g(µ) = t·|µ| + indicator of [lo, hi]`,   `t = λ·wR ≥ 0`,

and all the nonsmoothness of the algorithm is in the proximal map of `g`. The paper treats
`lo = −∞`, `hi = +∞` (Eq. 18); bounds get one sentence there (Sec. 8.2), so the bounded variant
below is derived from the canonical problem: in one dimension the proximal map of
`f + indicator of [lo, hi]` is the proximal map of `f` followed by the projection to `[lo, hi]`.

In the executable an absent bound is the IEEE infinity `±∞`, for which `clamp` is the identity;
a forbidden cone is `t = 0, lo = hi = 0`, and a prescribed cone `c` is `t = 0, lo = hi = c`.
-/

namespace Ddg

variable {α : Type} [Scalar α]

/-- The soft threshold `P_t(x) = max(0, x − t) + min(0, x + t)` (paper Eq. 18), the proximal map
of `t·|·|`. -/
def softThreshold (t x : α) : α := max 0 (x - t) + min 0 (x + t)

/-- Projection of `x` to the interval `[lo, hi]`. -/
def clamp (lo hi x : α) : α := min hi (max lo x)

/-- Bounded variant of Eq. 18: the proximal map of `t·|·| + indicator of [lo, hi]`.
The optimality condition `φ ∈ ∂g(µ)` is equivalent to the fixed-point equation
`µ = boundedThreshold t lo hi (φ + µ)` (paper Eq. 19 for the unbounded case). -/
def boundedThreshold (t lo hi x : α) : α := clamp lo hi (softThreshold t x)

/-- The regularized measure `µ = clamp(k · P_t(φ))` with `k = m/γ` (Stage I). It is the unique
solution of `φ − µ/k ∈ ∂g(µ)`, the optimality condition of the Moreau–Yosida regularized
problem with the term `(γ/2)·µ²/m` added per vertex. -/
def regularizedMeasure (k t lo hi x : α) : α := clamp lo hi (k * softThreshold t x)

/-! ## Active sets (paper Eq. 24) -/

/-- `x` is in the positive active set `A⁺`: `x > t`. A free site (`t = 0`) is active for every
`x`, split as `x ≥ 0` here and `x < 0` in `ActiveNeg`: for `t = 0` the map `P_0` is the identity,
whose derivative is `1` also at `x = 0`. -/
def ActivePos (t x : α) : Prop := t < x ∨ (t ≤ 0 ∧ 0 ≤ x)

/-- `x` is in the negative active set `A⁻`: `x < −t`. -/
def ActiveNeg (t x : α) : Prop := x < -t

instance (t x : α) : Decidable (ActivePos t x) := by unfold ActivePos; infer_instance
instance (t x : α) : Decidable (ActiveNeg t x) := by unfold ActiveNeg; infer_instance

/-- The linear piece of `x ↦ clamp(k · P_t(x))` that a point lies on. -/
inductive ConeState where
  /-- `|x| ≤ t`: the value is `0` -/
  | zero
  /-- `x > t`: the value is `k (x − t)` -/
  | pos
  /-- `x < −t`: the value is `k (x + t)` -/
  | neg
  /-- clamped at the lower bound: the value is `lo` -/
  | lower
  /-- clamped at the upper bound: the value is `hi` -/
  | upper
  deriving Repr, BEq, DecidableEq, Inhabited

/-- The piece of `x ↦ clamp lo hi (k · P_t(x))` at `x`. With `k = 1` and `x = φ + µ` these are
the active sets of Stage II, with `k = m/γ` and `x = φ` those of Stage I. -/
def coneState (k t lo hi x : α) : ConeState :=
  if hi ≤ k * softThreshold t x then .upper
  else if k * softThreshold t x ≤ lo then .lower
  else if ActivePos t x then .pos
  else if ActiveNeg t x then .neg
  else .zero

namespace ConeState

/-- The entry `D_ii ∈ {0, 1}` of the semismooth differential of `P` (paper Sec. 5.3.2):
`1` on `A⁺ ∪ A⁻`, `0` elsewhere (including the clamped states). -/
def indicator : ConeState → α
  | .pos | .neg => 1
  | _ => 0

/-- `true` on `A⁺ ∪ A⁻`. -/
def isActive : ConeState → Bool
  | .pos | .neg => true
  | _ => false

/-- Value of `clamp lo hi (k · P_t(x))` on a given piece. -/
def value (k t lo hi x : α) : ConeState → α
  | .zero => 0
  | .pos => k * (x - t)
  | .neg => k * (x + t)
  | .lower => lo
  | .upper => hi

end ConeState

/-! ## The subdifferential inclusion itself -/

/-- `φ ∈ ∂g(µ)` for `g = t·|·| + indicator of [lo, hi]`: `µ` is feasible, and `φ` lies in
`∂(t|·|)(µ)` plus the normal cone of `[lo, hi]` at `µ`. Spelled out:
the upper end of the interval is `t` (or `−t` if `µ < 0`) unless `µ = hi`, where it is `+∞`;
the lower end is `−t` (or `t` if `µ > 0`) unless `µ = lo`, where it is `−∞`. -/
def InSubdiff (t lo hi mu phi : α) : Prop :=
  lo ≤ mu ∧ mu ≤ hi ∧
  (mu < hi → phi ≤ (if mu < 0 then -t else t)) ∧
  (lo < mu → (if 0 < mu then t else -t) ≤ phi)

/-- Distance from `φ` to the interval `∂g(µ)` (for feasible `µ`); `0` iff `InSubdiff`. -/
def subdiffViolation (t lo hi mu phi : α) : α :=
  max (if mu < hi then max 0 (phi - (if mu < 0 then -t else t)) else 0)
      (if lo < mu then max 0 ((if 0 < mu then t else -t) - phi) else 0)

/-- Distance from `µ` to the interval `[lo, hi]`. -/
def boundViolation (lo hi mu : α) : α := max (max 0 (lo - mu)) (mu - hi)

end Ddg
