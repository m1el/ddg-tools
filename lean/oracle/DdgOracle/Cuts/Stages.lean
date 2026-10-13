import Std.Data.HashMap
import DdgOracle.Cuts.Region

/-!
# The default stages of one flow step (SPEC §3.8 – §3.10, §4)

The length step on the cut, the extension from the cut, the scaling by the mean gradient and the
roughening rate. The functions work on the cut mesh of the two-region code.
-/

namespace Ddg.Cuts

/-- Solve the chain system `(diag ℓ + Σ_σ κ_σ [[1,−1],[−1,1]]) x = rhs` over `n` nodes with the
couplings `(p, q, κ)`. -/
def chainSolve (n : Nat) (ell : FloatArray) (edges : Array (Nat × Nat × Float)) (rhs : FloatArray) :
    Except String FloatArray := do
  if n == 0 then return zeros 0
  let mut m := zeros (n * n)
  for k in [0:n] do m := m.set! (k * n + k) (fluxDiagEntry (ell.get! k) 0)
  for (p, q, kap) in edges do
    m := m.set! (p * n + p) (m.get! (p * n + p) + kap)
    m := m.set! (q * n + q) (m.get! (q * n + q) + kap)
    m := m.set! (p * n + q) (m.get! (p * n + q) + fluxOffEntry kap)
    m := m.set! (q * n + p) (m.get! (q * n + p) + fluxOffEntry kap)
  match luFactor n m with
  | none => throw "chain system is singular"
  | some lu => return lu.solve rhs

/-- `w` at a node of the cut mesh (original vertex or cut vertex). -/
def wNode (cm : CutMesh) (wv : FloatArray) (node : Nat) : Float :=
  if node < cm.nv then wv.get! node
  else
    let cr := cm.crossings[node - cm.nv]!
    (1 - cr.t) * wv.get! cr.a + cr.t * wv.get! cr.b

/-- The coupling constant `c` of a flux rule. -/
def fluxC (f : FluxRule) : Float := match f with | .nodal => 0.0 | .smoothed c => c

/-- Quantities of the length force at the cut vertices. -/
structure CutForce where
  /-- `w` at the cut vertex -/
  wc : FloatArray
  /-- dual length on the cut: half the lengths of its segments -/
  ellC : FloatArray
  /-- number of segments at the cut vertex (1 on `∂M`) -/
  deg : Array Nat
  /-- interior angle of region 0 at the cut vertex -/
  beta : FloatArray
  /-- `(∂w/∂n)` out of region 0 across each segment, in the order of `cm.segs` -/
  dwdn : FloatArray
  /-- integrated length force `Θ` -/
  theta : FloatArray
  deriving Inhabited

/-- Length force at the cut vertices (SPEC §3.10). `paperRule` drops the `∂w/∂n` term. -/
def cutForce (cm : CutMesh) (wv : FloatArray) (paperRule : Bool) : CutForce := Id.run do
  let nv := cm.nv
  let nc := cm.nc
  let mut ellC := zeros nc
  let mut deg : Array Nat := Array.replicate nc 0
  for sg in cm.segs do
    ellC := (ellC.set! sg.p (ellC.get! sg.p + sg.len / 2)).set! sg.q (ellC.get! sg.q + sg.len / 2)
    deg := (deg.set! sg.p (deg[sg.p]! + 1)).set! sg.q (deg[sg.q]! + 1)
  let wc := FloatArray.ofFn' nc fun k => wNode cm wv (nv + k)
  let mut beta := zeros nc
  for t in cm.tris do
    if t.region == 0 then
      for (n, c) in [(t.n0, t.c0), (t.n1, t.c1), (t.n2, t.c2)] do
        if n ≥ nv then beta := beta.set! (n - nv) (beta.get! (n - nv) + angleFromCot c)
  let mut dwdn := zeros cm.segs.size
  for si in [0:cm.segs.size] do
    let sg := cm.segs[si]!
    let t := cm.tris[sg.tri]!
    let lam := footFraction t.c1 t.c2
    let h := heightOver t.area sg.len
    let d := normalDerivativeOfLinear (wNode cm wv (nv + sg.p)) (wNode cm wv (nv + sg.q))
      (wv.get! sg.lone) lam h
    dwdn := dwdn.set! si (if cm.region[sg.lone]! == 0 then d else -d)
  let mut theta := zeros nc
  for k in [0:nc] do
    let gamma := if deg[k]! == 1 then Scalar.pi / 2 else Scalar.pi
    let mut segs : List (Float × Float) := []
    for si in [0:cm.segs.size] do
      let sg := cm.segs[si]!
      if sg.p == k || sg.q == k then
        segs := segs ++ [(sg.len, if paperRule then fzero else dwdn.get! si)]
    theta := theta.set! k (lengthForce (wc.get! k) gamma (beta.get! k) segs)
  return { wc, ellC, deg, beta, dwdn, theta }

/-- Connected pieces of the cut: the smallest cut-vertex index of the piece of every cut vertex. -/
def cutPieces (cm : CutMesh) : Array Nat := Id.run do
  let nc := cm.nc
  let mut lab : Array Nat := Array.range nc
  let mut changed := true
  while changed do
    changed := false
    for sg in cm.segs do
      let m := min lab[sg.p]! lab[sg.q]!
      if lab[sg.p]! != m then
        lab := lab.set! sg.p m
        changed := true
      if lab[sg.q]! != m then
        lab := lab.set! sg.q m
        changed := true
  return lab

/-- `ḡ` of every cut vertex: the mean `|∇φ|` over its piece (SPEC §3.9). -/
def meanGradients (cm : CutMesh) (phi : FloatArray) : FloatArray × Array Nat := Id.run do
  let piece := cutPieces cm
  let nc := cm.nc
  let mut sumL := zeros nc
  let mut sumLG := zeros nc
  for sg in cm.segs do
    let t := cm.tris[sg.tri]!
    let g := gradientNormLone (phi.get! sg.lone) sg.len t.area
    let r := piece[sg.p]!
    sumL := sumL.set! r (sumL.get! r + sg.len)
    sumLG := sumLG.set! r (sumLG.get! r + sg.len / g)
  return (FloatArray.ofFn' nc fun k => meanGradient (sumL.get! piece[k]!) (sumLG.get! piece[k]!), piece)

/-- The extension from the cut (SPEC §3.8, default): `σ` at every original vertex and the ring of
every vertex (`none` where unreachable). -/
def localExtension (b : Base) (cm : CutMesh) (ellC rg : FloatArray) :
    FloatArray × Array (Option Nat) := Id.run do
  let nv := cm.nv
  let mut sigma := zeros nv
  -- first line: transpose of the interpolation
  let mut ws : Array (List Float) := Array.replicate nv []
  let mut vs : Array (List Float) := Array.replicate nv []
  for k in [0:cm.nc] do
    let cr := cm.crossings[k]!
    for (i, beta) in [(cr.a, 1 - cr.t), (cr.b, cr.t)] do
      ws := ws.set! i (ws[i]! ++ [ellC.get! k * beta])
      vs := vs.set! i (vs[i]! ++ [rg.get! k])
  let mut ring : Array (Option Nat) := Array.replicate nv none
  for i in [0:nv] do
    if cm.adj[i]! then
      ring := ring.set! i (some 0)
      sigma := sigma.set! i (weightedMean ws[i]! vs[i]!)
  -- second line: ring by ring
  let mut level := 0
  let mut go := true
  while go do
    let mut next : Array Nat := #[]
    for i in [0:nv] do
      if ring[i]! == some level then
        for j in b.topo.neighbors[i]! do
          if ring[j]!.isNone && cm.region[j]! == cm.region[i]! then
            ring := ring.set! j (some (level + 1))
            next := next.push j
    if next.size == 0 then go := false
    else
      for i in next do
        let mut nb : List (Float × Float) := []
        for j in b.topo.neighbors[i]! do
          if ring[j]! == some level && cm.region[j]! == cm.region[i]! then
            nb := nb ++ [(sigma.get! j, b.l i j)]
        sigma := sigma.set! i (ringMean nb)
      level := level + 1
  return (sigma, ring)

/-- The harmonic extension (SPEC §3.8, the paper's): rates `r` at the cut vertices. -/
def harmonicExtension (cm : CutMesh) (regs : Array RegionData) (r : FloatArray) :
    Except String FloatArray := do
  let nv := cm.nv
  let N := cm.nodes
  let mut sigma := zeros nv
  for rd in regs do
    let mut groups : Std.HashMap Nat (Array Nat) := {}
    for i in [0:N] do
      if rd.inS[i]! then
        groups := groups.insert rd.comp[i]! ((groups.getD rd.comp[i]! #[]).push i)
    for (_, members) in groups.toList do
      let O := members.filter (· < nv)
      let C := members.filter (· ≥ nv)
      if C.size > 0 && O.size > 0 then
        let rhs := FloatArray.ofFn' O.size fun k => Id.run do
          let mut sum := fzero
          for cn in C do
            sum := sum + rd.lap.get O[k]! cn * r.get! (cn - nv)
          return -sum
        let some sg := rd.lap.solveSub O rhs | throw s!"region {rd.r}: singular extension system"
        for k in [0:O.size] do sigma := sigma.set! O[k]! (sg.get! k)
  return sigma

end Ddg.Cuts
