import DdgOracle.Formulas
import DdgOracle.Hex

/-!
# The input of one flow step

See `README.md`, section "Cut flow", for the text format; this file is the data model.
-/

namespace Ddg.Cuts

/-- Largest number of vertices the cut-flow oracle accepts. -/
def maxCutVertices : Nat := 1500

/-- The distortion energy `J` (SPEC §1). -/
inductive Energy where
  /-- `J = ∫ |∇u|²` -/
  | dirichlet
  /-- pointwise `b(u) = u²` -/
  | hencky
  /-- `J = ∫ a |∇u|²` with `a` given per vertex -/
  | scaledA (a : FloatArray)
  /-- `J = ∫ a |∇u|²` with the curvature-based `a` of SPEC §3.6 (`r_ε`, `c`) -/
  | scaledCurv (rEps c : Float)

/-- Source of the curvature density that drives the Yamabe equation. -/
inductive CurvSource where
  /-- `k_f = Ω_f / A_f` -/
  | gaussian
  /-- the same density `c` on every face -/
  | uniform (c : Float)
  /-- one density per face -/
  | perFace (k : FloatArray)

/-- Penalty rule for `η` (SPEC §3.10). -/
inductive Penalty where
  /-- `(W_V + τ L_w) φ′ = W_V φ̃`, `w = α_L + η` -/
  | gradient
  /-- `(W_V D⁻¹ + τ L) φ′ = W_V D⁻¹ φ̃`, `D = diag(α_L + η)` -/
  | paperDiagonal

/-- Step size rule (SPEC §4). -/
inductive StepRule where
  /-- `τ` given -/
  | fixed (tau : Float)
  /-- `τ = min(τ_max, c · h / max|s|)` -/
  | capped (tauMax maxMotion : Float)

/-- One case: a mesh given by faces and edge lengths (caller's unit) and everything else of one
flow step. All parameters refer to the surface scaled to unit total area, `φ` is in the caller's
unit. -/
structure CutCase where
  /-- identifier -/
  name : String
  /-- number of vertices -/
  nv : Nat
  /-- oriented triangles -/
  faces : Array Face
  /-- edges `(a, b, length)` with `a < b` -/
  edges : Array (Nat × Nat × Float)
  /-- level set function, caller's unit -/
  phi : FloatArray
  /-- distortion energy -/
  energy : Energy
  /-- curvature source -/
  curv : CurvSource
  /-- length weight `α_L` -/
  alphaL : Float
  /-- penalty `η` per vertex -/
  eta : Option FloatArray
  /-- penalty rule -/
  penalty : Penalty
  /-- distortion weights `w_0, w_1` -/
  w : Float × Float
  /-- area weights `β_0, β_1` -/
  beta : Float × Float
  /-- step size rule -/
  step : StepRule
  /-- `min_crossing_offset` -/
  delta : Float
  /-- sliver threshold `ρ₀` -/
  rho0 : Float
  /-- number of smoothing passes -/
  passes : Nat
  /-- symmetry class label of every vertex -/
  classes : Option (Array Nat)

end Ddg.Cuts
