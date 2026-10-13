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
  deriving BEq

/-- How the integrated fluxes are turned into normal derivatives (SPEC §3.5). -/
inductive FluxRule where
  /-- `d_p = v_p / ℓ_p` (the paper) -/
  | nodal
  /-- flux recovery with `c` -/
  | smoothed (c : Float)

/-- Step size rule (SPEC §4). -/
inductive StepRule where
  /-- `τ` given -/
  | fixed (tau : Float)
  /-- `τ = min(τ_max, c · h / max|ρ|)` -/
  | capped (tauMax maxMotion : Float)
  /-- `τ = min(τ_max, time_scale / ω, max_motion · h / max|ρ|)` (format 2) -/
  | adaptive (timeScale maxMotion tauMax : Float)

/-- Where the length term acts (SPEC §3.10). -/
inductive LengthStep where
  /-- `(W_V + τ L_w) φ′ = W_V φ̃` after the explicit step (the paper) -/
  | onMesh
  /-- length force and linearly implicit rate on the cut (default) -/
  | onCut
  deriving BEq

/-- How the rate on the cut becomes a rate at the vertices (SPEC §3.8). -/
inductive Extension where
  /-- harmonic in each submesh (the paper) -/
  | harmonic
  /-- from the cut: transpose of the interpolation, then ring by ring (default) -/
  | local
  deriving BEq

/-- Scaling of the rate by the gradient of `φ` (SPEC §3.9). -/
inductive SpeedScaling where
  /-- `ḡ = 1` (the paper) -/
  | unit
  /-- `ḡ` = mean `|∇φ|` over the piece of the cut (default) -/
  | meanGradient
  deriving BEq

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
  /-- flux recovery rule -/
  flux : FluxRule
  /-- sliver threshold `ρ₀` -/
  rho0 : Float
  /-- number of smoothing passes -/
  passes : Nat
  /-- symmetry class label of every vertex -/
  classes : Option (Array Nat)
  /-- text format version: 1 (the stages of the paper only) or 2 -/
  format : Nat := 1
  /-- where the length term acts -/
  lengthStep : LengthStep := .onMesh
  /-- extension of the rate -/
  extension : Extension := .harmonic
  /-- scaling by `|∇φ|` -/
  speed : SpeedScaling := .unit
  /-- feature size smoothing with the positive edge weights `L⁺` (SPEC §3.6, D17) -/
  featPos : Bool := false

end Ddg.Cuts
