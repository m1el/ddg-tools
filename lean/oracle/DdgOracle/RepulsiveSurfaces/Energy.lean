import DdgOracle.RepulsiveSurfaces.Basic

/-!
# The exact tangent-point energy, its gradient, and `D_f` (SPEC §3.1, §3.2, §4.1)
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

/-- `|ν|^{α−2} ν = sign(ν) |ν|^{α−1}`, continuous at `ν = 0` for `α > 1`. -/
def signPow (nu alpha : Float) : Float :=
  if nu == 0.0 then 0.0 else (if nu > 0.0 then 1.0 else -1.0) * Float.pow nu.abs (alpha - 1.0)

/-- The kernel `K(X, N; Y) = |N·(X−Y)|^α |X−Y|^{−β}` (SPEC §3.1). -/
def kernel (ex : Ex) (x n y : V3) : Float :=
  let v := x - y
  Float.pow (n.dot v).abs ex.alpha * Float.pow v.norm (-ex.beta)

/-- Exact energy `Σ_{S≠T} K(X_S,N_S;X_T) a_S a_T` (ordered pairs, neighbours included). Error
`NonFiniteEnergy` when two barycenters coincide. -/
def energy (faces : Array Tri) (g : Array FaceGeo) (ex : Ex) : Except String Float := do
  let mut e := 0.0
  for s in [0:faces.size] do
    for t in [0:faces.size] do
      if s != t then
        let v := g[s]!.x - g[t]!.x
        if v.norm2 == 0.0 then throw s!"NonFiniteEnergy: faces ({s}, {t}) have the same barycenter"
        e := e + kernel ex g[s]!.x g[s]!.n g[t]!.x * g[s]!.area * g[t]!.area
  if !e.isFinite then throw "NonFiniteEnergy: the energy is not finite"
  return e

/-- Energy and gradient `∂ℰ̂/∂f_v` (SPEC §3.2). -/
def energyGrad (faces : Array Tri) (pos : Array V3) (ex : Ex) : Except String (Float × Array V3) := do
  let g := faceGeo faces pos
  let nf := faces.size
  let mut e := 0.0
  let mut dX : Array V3 := Array.replicate nf V3.zero
  let mut dN : Array V3 := Array.replicate nf V3.zero
  let mut dA : Array Float := Array.replicate nf 0.0
  for s in [0:nf] do
    for t in [0:nf] do
      if s != t then
        let v := g[s]!.x - g[t]!.x
        let r2 := v.norm2
        if r2 == 0.0 then throw s!"NonFiniteEnergy: faces ({s}, {t}) have the same barycenter"
        let r := Float.sqrt r2
        let n := g[s]!.n
        let nu := n.dot v
        let rb := Float.pow r (-ex.beta)
        let k := Float.pow nu.abs ex.alpha * rb
        let aa := g[s]!.area * g[t]!.area
        e := e + k * aa
        let sp := signPow nu ex.alpha
        -- ∂K/∂v = α|ν|^{α−2}ν r^{−β} N − β|ν|^α r^{−β−2} v
        let dKdv := (ex.alpha * sp * rb) * n - (ex.beta * Float.pow nu.abs ex.alpha * rb / r2) * v
        let dKdN := (ex.alpha * sp * rb) * v
        dX := dX.modify s (· + aa * dKdv)
        dX := dX.modify t (· - aa * dKdv)
        dN := dN.modify s (· + aa * dKdN)
        dA := dA.modify s (· + k * g[t]!.area)
        dA := dA.modify t (· + k * g[s]!.area)
  if !e.isFinite then throw "NonFiniteEnergy: the energy is not finite"
  let mut grad : Array V3 := Array.replicate pos.size V3.zero
  for fi in [0:nf] do
    let f := faces[fi]!
    let n := g[fi]!.n
    let a := g[fi]!.area
    for t in [0:3] do
      let ee := oppEdge pos f t
      let w := ee.cross n
      let contrib := (1.0 / 3.0 : Float) * dX[fi]! + (dA[fi]! * 0.5) * n.cross ee
        + ((w.dot dN[fi]!) / (2.0 * a)) * n
      grad := grad.modify (f.get t) (· + contrib)
  return (e, grad)

/-- The P1 gradient of a vertex function on one face: `D_f u = (1/2a) N × (u_i e_jk + u_j e_ki + u_k e_ij)`
(un-normalized edges, SPEC §4.1). -/
def dfFace (pos : Array V3) (f : Tri) (g : FaceGeo) (u : Array Float) : V3 := Id.run do
  let mut s := V3.zero
  for t in [0:3] do
    s := s + u[f.get t]! * oppEdge pos f t
  return (1.0 / (2.0 * g.area)) * g.n.cross s

/-- The coefficient vector of `u_{f.get t}` in `D_f u`: `(1/2a) N × e_opp`. -/
def dfCoef (pos : Array V3) (f : Tri) (g : FaceGeo) (t : Nat) : V3 :=
  (1.0 / (2.0 * g.area)) * g.n.cross (oppEdge pos f t)

end Ddg.RepulsiveSurfaces
