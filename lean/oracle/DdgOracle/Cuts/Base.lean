import DdgOracle.Topology
import DdgOracle.Cuts.Case
import DdgOracle.Cuts.Mat

/-!
# The mesh `M` of the cut flow, scaled to unit area

Everything is intrinsic: the input is faces and edge lengths. Lengths are multiplied by
`s = 1/√(total area)` (SPEC §2); all quantities below live in the scaled unit.
-/

namespace Ddg.Cuts

/-- The scaled mesh and its operators. -/
structure Base where
  /-- number of vertices -/
  nv : Nat
  /-- connectivity (the mesh may be disconnected) -/
  topo : Topology
  /-- `s = 1/√(total area)`, multiplies the caller's lengths -/
  scale : Float
  /-- total area in the caller's unit -/
  totalArea : Float
  /-- scaled edge lengths, row-major `nv × nv`, `0` where there is no edge -/
  len : FloatArray
  /-- area of every face -/
  faceArea : Array Float
  /-- cotangents of the corners at `a`, `b`, `c` of every face -/
  cots : Array (Float × Float × Float)
  /-- `Θ_i`, the sum of the corner angles at vertex `i` -/
  theta : FloatArray
  /-- `Ω_f` -/
  omega : Array Float
  /-- `Ω_f / A_f` -/
  kGauss : Array Float
  /-- the density that drives the Yamabe equation (after a curvature override) -/
  kUsed : Array Float
  /-- lumped vertex area `A_i` -/
  vertArea : FloatArray
  /-- cotan Laplacian of `M` -/
  lap : Mat
  /-- feature size `r` (zeros unless the energy is curvature-scaled) -/
  feat : FloatArray
  /-- the weight `a` per vertex (zeros unless the energy is scaled) -/
  aVert : FloatArray

namespace Base

/-- Length of the edge `{i, j}` (only meaningful on edges). -/
def l (b : Base) (i j : Nat) : Float := b.len.get! (i * b.nv + j)

end Base

/-- Is the energy one of the scaled kind. -/
def Energy.isScaled : Energy → Bool
  | .scaledA _ => true
  | .scaledCurv _ _ => true
  | _ => false

/-- Validate the case and build the scaled mesh, its curvature, vertex areas, Laplacian and (when
needed) the weight `a`. -/
def buildBase (c : CutCase) : Except String Base := do
  let nv := c.nv
  let topo ← buildTopology nv c.faces (connected := false)
  -- edge table
  let mut raw : FloatArray := zeros (nv * nv)
  let mut nEdges := 0
  for (a, b, l) in c.edges do
    if !(a < b && b < nv) then throw s!"edge ({a}, {b}): expected a < b < {nv}"
    if !(l.isFinite && l > 0) then throw s!"edge ({a}, {b}) has a length that is not finite and positive"
    if raw.get! (a * nv + b) != 0 then throw s!"edge ({a}, {b}) is listed twice"
    raw := (raw.set! (a * nv + b) l).set! (b * nv + a) l
    nEdges := nEdges + 1
  if nEdges != topo.numEdges then
    throw s!"{nEdges} edge lengths given but the faces have {topo.numEdges} edges"
  for fi in [0:c.faces.size] do
    let f := c.faces[fi]!
    for (p, q) in [(f.a, f.b), (f.b, f.c), (f.c, f.a)] do
      if raw.get! (p * nv + q) == 0 then
        throw s!"face {fi}: the length of edge ({min p q}, {max p q}) is missing"
  -- total area in the caller's unit
  let rawLen := fun (i j : Nat) => raw.get! (i * nv + j)
  let mut total := fzero
  for fi in [0:c.faces.size] do
    let f := c.faces[fi]!
    let a := rawLen f.b f.c
    let b := rawLen f.c f.a
    let cc := rawLen f.a f.b
    if !(b - a + cc > 0 && a - b + cc > 0 && a - cc + b > 0) then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) is degenerate: its edge lengths violate the strict triangle inequality"
    total := total + f.area rawLen
  if !(total > 0 && total.isFinite) then throw "total area is not finite and positive"
  let s := 1 / Float.sqrt total
  let mut len := zeros (nv * nv)
  for i in [0:nv * nv] do
    len := len.set! i (s * raw.get! i)
  let lenF := fun (i j : Nat) => len.get! (i * nv + j)
  -- faces
  let nf := c.faces.size
  let mut faceArea : Array Float := Array.mkEmpty nf
  let mut cots : Array (Float × Float × Float) := Array.mkEmpty nf
  for fi in [0:nf] do
    let f := c.faces[fi]!
    let ar := f.area lenF
    if !(ar > 0 && ar.isFinite && (f.cotA lenF).isFinite && (f.cotB lenF).isFinite && (f.cotC lenF).isFinite) then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) is numerically degenerate"
    faceArea := faceArea.push ar
    cots := cots.push (f.cotA lenF, f.cotB lenF, f.cotC lenF)
  -- vertices
  let mut theta := zeros nv
  let mut vertArea := zeros nv
  for v in [0:nv] do
    let star := topo.facesAt v
    theta := theta.set! v (angleSum star lenF v)
    vertArea := vertArea.set! v (areaAt star lenF v / 3)
  let isB := fun (v : Nat) => topo.isBoundary[v]!
  let thetaF := fun (v : Nat) => theta.get! v
  let mut omega : Array Float := Array.mkEmpty nf
  let mut kGauss : Array Float := Array.mkEmpty nf
  for fi in [0:nf] do
    let o := (c.faces[fi]!).omega lenF thetaF isB
    omega := omega.push o
    kGauss := kGauss.push (faceCurvatureDensity o faceArea[fi]!)
  let kUsed ← match c.curv with
    | .gaussian => pure kGauss
    | .uniform k =>
      if !k.isFinite then throw "uniform curvature is not finite"
      pure (Array.replicate nf k)
    | .perFace ks =>
      if ks.size != nf then throw s!"per-face curvature has {ks.size} values for {nf} faces"
      let mut out : Array Float := Array.mkEmpty nf
      for i in [0:nf] do
        if !(ks.get! i).isFinite then throw s!"per-face curvature of face {i} is not finite"
        out := out.push (ks.get! i)
      pure out
  -- cotan Laplacian of M
  let mut lap := Mat.zero nv
  for fi in [0:nf] do
    let f := c.faces[fi]!
    let (c0, c1, c2) := cots[fi]!
    lap := lap.addTri f.a f.b f.c c0 c1 c2 1
  -- weight a
  let mut feat := zeros nv
  let mut aVert := zeros nv
  match c.energy with
  | .scaledA a =>
    if a.size != nv then throw s!"scaled energy: {a.size} weights for {nv} vertices"
    for i in [0:nv] do
      if !((a.get! i).isFinite && a.get! i > 0) then throw s!"scaled energy: a[{i}] is not finite and positive"
    aVert := a
  | .scaledCurv rEps cc =>
    if !(rEps.isFinite && rEps > 0) then throw "r_eps must be finite and positive"
    if !(cc.isFinite && cc > 0) then throw "c must be finite and positive"
    let mut rhs := zeros nv
    let mut sys := zeros (nv * nv)
    for i in [0:nv] do
      let kd := if isB i then fzero else vertexCurvatureDensity (theta.get! i) (vertArea.get! i)
      let r0 := featureSizeRaw rEps kd
      rhs := rhs.set! i (vertArea.get! i * r0)
      for j in [0:nv] do
        sys := sys.set! (i * nv + j) (cc * lap.get i j + (if i == j then vertArea.get! i else fzero))
    match luFactor nv sys with
    | none => throw "curvature-based weight: singular system"
    | some lu =>
      feat := lu.solve rhs
      for i in [0:nv] do
        aVert := aVert.set! i (featureWeight (feat.get! i))
  | _ => pure ()
  return { nv, topo, scale := s, totalArea := total, len, faceArea, cots, theta, omega, kGauss
           kUsed, vertArea, lap, feat, aVert }

end Ddg.Cuts
