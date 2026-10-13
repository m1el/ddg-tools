import DdgOracle.Regions.Cut

/-!
# Region solves, speed and energy for N regions (SPEC §8.3, §6)

Every region is solved by the code of the two-region case on its own view of the cut mesh: the
nodes are the original vertices and the copies of cut and triple points in that region, renumbered
`nv + j`.
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- A region's view of the cut mesh and its solve. -/
structure RView where
  /-- the region -/
  X : Nat
  /-- local cut mesh -/
  cm : CutMesh
  /-- local cut index `j` to the copy in the global cut mesh -/
  toCopy : Array Nat
  /-- copy to local cut index -/
  localOf : Std.HashMap Nat Nat
  /-- the solve -/
  rd : RegionData
  /-- weight `a` at the local nodes (empty unless scaled) -/
  aNode : FloatArray
  deriving Inhabited

/-- The geometric id of a copy and, for a triple point, its corner. -/
def geoOfCopy (net : Net) (c : Nat) : Nat × Nat :=
  if c < net.nv + 2 * net.nP then (net.nv + (c - net.nv) / 2, 0)
  else
    let r := c - net.nv - 2 * net.nP
    (net.nv + net.nP + r / 3, r % 3)

/-- `x` at a geometric node from the values at the vertices (cut point: linear on the edge with the
moved `t`; triple point: barycentric). -/
def valueAtGeo (net : Net) (faces : Array Face) (vals : FloatArray) (g : Nat) : Float :=
  if g < net.nv then vals.get! g
  else if g < net.nv + net.nP then
    let cr := net.cps[g - net.nv]!
    (1 - cr.t) * vals.get! cr.a + cr.t * vals.get! cr.b
  else
    let t := net.triples[g - net.nv - net.nP]!
    let f := faces[t.face]!
    t.b0 * vals.get! f.a + t.b1 * vals.get! f.b + t.b2 * vals.get! f.c

/-- Build the view of region `X` and solve it. -/
def solveRegion (b : Base) (c : RegionCase) (net : Net) (X : Nat) : Except String RView := do
  let nv := net.nv
  let trisX := net.tris.filter (·.region == X)
  let mut copies : Array Nat := #[]
  for t in trisX do
    for n in [t.n0, t.n1, t.n2] do
      if n ≥ nv && !copies.contains n then copies := copies.push n
  copies := copies.qsort (· < ·)
  let mut localOf : Std.HashMap Nat Nat := {}
  for j in [0:copies.size] do localOf := localOf.insert copies[j]! j
  let loc := fun (n : Nat) => if n < nv then n else nv + localOf.getD n 0
  let tris := trisX.map fun t => { t with n0 := loc t.n0, n1 := loc t.n1, n2 := loc t.n2 }
  let mut segs : Array Seg := #[]
  for s in net.segs do
    if s.regA == X || s.regB == X then
      segs := segs.push { p := localOf.getD (net.copyIdx s.g1 X) 0, q := localOf.getD (net.copyIdx s.g2 X) 0
                          face := s.face, len := s.len, tri := 0, lone := 0 }
  let cm : CutMesh := { nv, region := net.region, crossings := Array.replicate copies.size default
                        tris, segs, adj := net.adj, nearTies := 0 }
  let aGeo := fun (g : Nat) => valueAtGeo net c.faces b.aVert g
  let aNode := if c.energy.isScaled then
      FloatArray.ofFn' (nv + copies.size) fun i =>
        if i < nv then b.aVert.get! i else aGeo (geoOfCopy net copies[i - nv]!).1
    else zeros 0
  let rd ← computeRegion b cm c.energy aNode b.kUsed c.flux X
  return { X, cm, toCopy := copies, localOf, rd, aNode }

/-- Evaluation of the network. -/
structure NEval where
  /-- the network -/
  net : Net
  /-- `w = α_L + η` at the vertices -/
  wv : FloatArray
  /-- the region solves -/
  views : Array RView
  /-- speed `s_p = G_i − G_j` at every cut point (`i < j` the regions of its edge) -/
  sCut : FloatArray
  /-- lower and higher region of every cut point -/
  lowHigh : Array (Nat × Nat)
  /-- roughening rate at every cut point -/
  omegaP : FloatArray
  /-- rate `s_v` at the copies of every region, by local cut index -/
  sNode : Array FloatArray
  /-- distortion energy per region -/
  Js : Array Float
  /-- area of the solved part per region -/
  As : Array Float
  /-- length term -/
  lengthTerm : Float
  /-- total energy -/
  total : Float
  deriving Inhabited

/-- Evaluate everything that does not need `τ` at the scaled functions `phis`. -/
def evaluateNet (c : RegionCase) (b : Base) (phis : Array FloatArray) : Except String NEval := do
  let net ← buildNet b c.faces phis c.delta
  let nv := net.nv
  let N := net.nreg
  let wv := FloatArray.ofFn' nv fun i =>
    c.alphaL + (match c.eta with | some e => e.get! i | none => fzero)
  let mut views : Array RView := #[]
  for X in [0:N] do views := views.push (← solveRegion b c net X)
  let nP := net.nP
  -- G of region X at the copy of a cut point
  let gAt := fun (X p : Nat) =>
    let v := views[X]!
    let j := v.localOf.getD (net.copyIdx (nv + p) X) 0
    regionDensity c.ws[X]! (v.rd.G.get! (nv + j)) c.betas[X]!
  let qAt := fun (X p : Nat) =>
    let v := views[X]!
    let j := v.localOf.getD (net.copyIdx (nv + p) X) 0
    roughQ (regionDensity c.ws[X]! (v.rd.G.get! (nv + j)) c.betas[X]!) c.betas[X]! c.ws[X]! 0
  let mut sCut := zeros nP
  let mut omegaP := zeros nP
  let mut lowHigh : Array (Nat × Nat) := #[]
  for p in [0:nP] do
    let cr := net.cps[p]!
    let ra := net.region[cr.a]!
    let rb := net.region[cr.b]!
    let i := min ra rb
    let j := max ra rb
    lowHigh := lowHigh.push (i, j)
    sCut := sCut.set! p (cutSpeed (gAt i p) (gAt j p))
    omegaP := omegaP.set! p (roughRate (qAt i p) (qAt j p))
  -- rates at the copies
  let sInX := fun (X p : Nat) =>
    let (i, _) := lowHigh[p]!
    if X == i then sCut.get! p else -(sCut.get! p)
  let segLen := fun (p m : Nat) => Id.run do
    let mut l := fzero
    for s in net.segs do
      if s.toTriple && s.g1 == nv + p && s.g2 == nv + nP + m then l := s.len
    return l
  let mut sNode : Array FloatArray := #[]
  for X in [0:N] do
    let v := views[X]!
    let mut arr := zeros v.toCopy.size
    for j in [0:v.toCopy.size] do
      let (g, k) := geoOfCopy net v.toCopy[j]!
      if g < nv + nP then arr := arr.set! j (sInX X (g - nv))
      else
        let m := g - nv - nP
        let t := net.triples[m]!
        let pk := #[t.p0, t.p1, t.p2]
        let p := pk[k]!
        let p' := pk[(k + 2) % 3]!
        arr := arr.set! j (tripleRate (segLen p m) (segLen p' m) (sInX X p) (sInX X p'))
    sNode := sNode.push arr
  -- energy
  let mut lengthTerm := fzero
  for s in net.segs do
    lengthTerm := lengthTerm + (valueAtGeo net c.faces wv s.g1 + valueAtGeo net c.faces wv s.g2) / 2 * s.len
  let mut total := lengthTerm
  let mut Js : Array Float := #[]
  let mut As : Array Float := #[]
  for X in [0:N] do
    let rd := views[X]!.rd
    Js := Js.push rd.J
    As := As.push rd.areaSolved
    total := total + c.ws[X]! * rd.J + c.betas[X]! * rd.areaSolved
  return { net, wv, views, sCut, lowHigh, omegaP, sNode, Js, As, lengthTerm, total }

end Ddg.Regions
