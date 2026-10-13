import DdgOracle.Cuts

/-!
# N regions (SPEC §8): the input of one step
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- Which stages (SPEC §8.4 or §8.5). -/
inductive Stages where
  /-- the paper's: harmonic extension of `±s`, `φ̇_j = −φ̇_k`, `(W_V + τ L_w)` on every function -/
  | paper
  /-- the default stages on networks of SPEC §8.5 -/
  | dflt
  deriving BEq

/-- One case. All parameters refer to the surface scaled to unit area; the functions `phis` are in
the caller's unit. -/
structure RegionCase where
  /-- identifier -/
  name : String
  /-- number of vertices -/
  nv : Nat
  /-- oriented triangles -/
  faces : Array Face
  /-- edges `(a, b, length)`, `a < b` -/
  edges : Array (Nat × Nat × Float)
  /-- the `N` functions `φ_k`, caller's unit -/
  phis : Array FloatArray
  /-- distortion energy -/
  energy : Energy
  /-- curvature source -/
  curv : CurvSource
  /-- length weight -/
  alphaL : Float
  /-- penalty -/
  eta : Option FloatArray
  /-- penalty rule -/
  penalty : Penalty
  /-- distortion weights `w_k` -/
  ws : Array Float
  /-- area weights `β_k` -/
  betas : Array Float
  /-- step size rule -/
  step : StepRule
  /-- `δ` -/
  delta : Float
  /-- flux recovery -/
  flux : FluxRule
  /-- stages -/
  stages : Stages
  /-- positive feature smoothing -/
  featPos : Bool
  /-- symmetry classes -/
  classes : Option (Array Nat)
  /-- text format version: 1, or 2 (statistics of the step rule) -/
  format : Nat := 1
  /-- statistics of the step rule -/
  stats : StepStats := .extremes

/-- The two-region case of the same data, for the base mesh. -/
def RegionCase.toCut (c : RegionCase) : CutCase :=
  { name := c.name, nv := c.nv, faces := c.faces, edges := c.edges, phi := c.phis[0]!
    energy := c.energy, curv := c.curv, alphaL := c.alphaL, eta := c.eta, penalty := c.penalty
    w := (c.ws[0]!, c.ws[1]!), beta := (c.betas[0]!, c.betas[1]!), step := c.step
    delta := c.delta, flux := c.flux, rho0 := 5.0, passes := 0, classes := c.classes
    featPos := c.featPos }

end Ddg.Regions
