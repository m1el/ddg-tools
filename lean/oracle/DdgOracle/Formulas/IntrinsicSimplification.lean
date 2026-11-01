import DdgOracle.Formulas.Triangle

/-!
# Scalar-generic formulas of intrinsic simplification (`crates/ddg-intrinsic-simplification/SPEC.md`)

"Surface Simplification using Intrinsic Error Metrics", Liu, Gillespie, Chislett, Sharp, Jacobson,
Crane 2023. Every definition carries the id of the formula it belongs to as `[formula: iem-…]`;
`lean/proofs/DdgProofs/IntrinsicSimplification/Oracle.lean` identifies them at `ℝ` with the
definitions the theorems are about. The executable (`DdgOracle.IntrinsicSimplification.*`) uses
them for every number it prints; exponentials and trigonometric functions (not in `Ddg.Scalar`)
are evaluated by the executable and passed in.

Complex numbers are pairs `(re, im)`.
-/

namespace Ddg

variable {α : Type} [Scalar α]

/-- The extended corner angle of CETM Sec. 3.1, opposite the side `a`: `π` if `a ≥ b + c`, `0` if
`b ≥ a + c` or `c ≥ a + b`, else the angle of the triangle. [formula: iem-flatten-monotone] -/
def iemExtAngle (a b c : α) : α :=
  if b + c ≤ a then Scalar.pi
  else if a + c ≤ b then 0
  else if a + b ≤ c then 0
  else cornerAngle a b c

/-- Eq. 4: `u − (Θ̃ − Θ)/H` with `H = ½ Σ_f hessTerm f`. [formula: iem-flatten-newton] -/
def iemNewtonStep (u target theta hess : α) : α := u - (target - theta) / hess

/-- Eq. 9 for one neighbour, with S11's flat-vertex rule: `|ΔK_j| / Σ_l |ΔK_l|`, or `1/|N_i|` if
`|K_i| ≤ tol`. [formula: iem-redistribution-weights] -/
def iemWeight (dKj sumAbs Ki tol count : α) : α :=
  if Scalar.abs Ki ≤ tol then 1 / count else Scalar.abs dKj / sumAbs

/-- Eq. 5: `m̃_j = m_j + α_j m_i`. [formula: iem-memoryless-cost] -/
def iemMassUpdate (mj a mi : α) : α := mj + a * mi

/-- Product of two complex numbers given as pairs. -/
def cmul (x y : α × α) : α × α := (x.1 * y.1 - x.2 * y.2, x.1 * y.2 + x.2 * y.1)

/-- Eq. 8 (SPEC §4.3): `t̃_j = (α m_i (R t_i + e_ji) + m_j t_j)/m̃_j`, and `t̃_j = t_j` unless
`m̃_j > 0`. `R` and `e` are given as pairs. [formula: iem-center-of-mass-update] -/
def iemVectorUpdate (mj a mi : α) (R ti e tj : α × α) : α × α :=
  let m := iemMassUpdate mj a mi
  if 0 < m then
    let v := cmul R ti
    let w := a * mi
    ((w * (v.1 + e.1) + mj * tj.1) / m, (w * (v.2 + e.2) + mj * tj.2) / m)
  else tj

/-- One channel's term of the cost of Eq. 6 (first form): `m̃_j ‖t̃_j‖`.
[formula: iem-memoryless-cost] -/
def iemCostTerm (m : α) (t : α × α) : α := m * Scalar.sqrt (t.1 * t.1 + t.2 * t.2)

/-- Eq. 10 with corner weights `w` (`e^{u}` at the corners equal to `i`, `1` elsewhere).
[formula: iem-projective-interpolation] -/
def iemProjInterp (w0 w1 w2 b0 b1 b2 : α) : α × α × α :=
  let s := w0 * b0 + w1 * b1 + w2 * b2
  (w0 * b0 / s, w1 * b1 / s, w2 * b2 / s)

end Ddg
