import DdgOracle.Regions.Step
import DdgOracle.Cuts.Gen

/-!
# Test states for N regions
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- The N-region case that holds the same data as a two-region case: `φ_0 = φ`, `φ_1 = −φ`. -/
def ofCut (c : CutCase) (stages : Stages) : RegionCase :=
  { name := c.name, nv := c.nv, faces := c.faces, edges := c.edges
    phis := #[c.phi, FloatArray.ofFn' c.phi.size fun i => -(c.phi.get! i)]
    energy := c.energy, curv := c.curv, alphaL := c.alphaL, eta := c.eta, penalty := c.penalty
    ws := #[c.w.1, c.w.2], betas := #[c.beta.1, c.beta.2], step := c.step, delta := c.delta
    flux := c.flux, stages, featPos := c.featPos, classes := c.classes }

/-- A region case with default parameters (Dirichlet, Gaussian curvature, `α_L = 0.05`, flux
recovery `c = 0.5`, `w_k = 1`, `β_k = 0`, `δ = 10⁻³`, the paper's stages, fixed `τ = 0.001`). -/
def mkRegionCase (name : String) (m : RawMesh) (phis : Array FloatArray) : RegionCase :=
  let c := mkCase name m (phis[0]!)
  { name, nv := c.nv, faces := c.faces, edges := c.edges, phis
    energy := .dirichlet, curv := .gaussian, alphaL := 0.05, eta := none, penalty := .gradient
    ws := Array.replicate phis.size 1.0, betas := Array.replicate phis.size 0.0
    step := .fixed 0.001, delta := 1e-3, flux := .smoothed 0.5, stages := .paper, featPos := false
    classes := none }

/-- The states of the Voronoi partition of the given centres (Euclidean distance in `ℝ³`):
`φ_k = (d_k − min_{j≠k} d_j)/2`, negative in the cell of `k`. -/
def voronoiState (m : RawMesh) (centers : Array (Point Float)) : Array FloatArray :=
  Array.ofFn (n := centers.size) fun k =>
    FloatArray.ofFn' m.pos.size fun i => Id.run do
      let d := fun (j : Nat) => dist m.pos[i]! centers[j]!
      let mut best := posInf
      for j in [0:centers.size] do
        if j != k.val && d j < best then best := d j
      return (d k.val - best) / 2

/-- Bands of the sphere between the polar angles `psis` (increasing): region `k` is
`ψ_{k-1} < ψ ≤ ψ_k`. `φ_k = max(ψ_{k-1} − ψ, ψ − ψ_k)` (a signed distance up to the Lipschitz
constant 1). -/
def bandState (m : RawMesh) (psis : Array Float) : Array FloatArray :=
  Array.ofFn (n := psis.size + 1) fun k =>
    FloatArray.ofFn' m.pos.size fun i =>
      let psi := Float.acos m.pos[i]!.z
      let lo := if k.val == 0 then negInf else psis[k.val - 1]!
      let hi := if k.val == psis.size then posInf else psis[k.val]!
      max (lo - psi) (psi - hi)

/-- A state from vertex labels `lab : vertex → region`: function `k` is `±` the distance in
`ℝ³` to the nearest vertex of the other side (`+` outside, `−` inside), halved like the
two-region helper (midpoint crossings). A region without vertices gets `+1` everywhere. -/
def labelState (m : RawMesh) (n : Nat) (lab : Array Nat) : Array FloatArray :=
  Array.ofFn (n := n) fun k =>
    FloatArray.ofFn' m.pos.size fun i => Id.run do
      let inK := lab[i]! == k.val
      let mut best := posInf
      for j in [0:m.pos.size] do
        if (lab[j]! == k.val) != inK then
          let d := dist m.pos[i]! m.pos[j]!
          if d < best then best := d
      return (if best == posInf then 1.0 else if inK then -(best / 2) else best / 2)

/-- Redistance a state (SPEC §8.4 step 2, with the projection) so that it is a signed distance near
the interfaces: the functions of `voronoiState` are exact only near their own interface. -/
def redistancedState (m : RawMesh) (phis : Array FloatArray) : Array FloatArray :=
  let c := mkRegionCase "tmp" m phis
  match buildBase c.toCut with
  | .error _ => phis
  | .ok b =>
    let sc := phis.map fun p => FloatArray.ofFn' p.size fun i => b.scale * p.get! i
    let r := redistanceN b c.faces sc
    r.phi.map fun p => FloatArray.ofFn' p.size fun i => p.get! i / b.scale

end Ddg.Regions
