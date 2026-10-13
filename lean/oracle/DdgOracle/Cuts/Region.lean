import DdgOracle.Cuts.CutMesh

/-!
# One region: scale factor, adjoint, gradient density (SPEC §3.4 – §3.6)
-/

namespace Ddg.Cuts

/-- Everything computed on the submesh `S_r` of one region. Arrays are indexed by node. -/
structure RegionData where
  /-- the region -/
  r : Nat
  /-- node belongs to `S` -/
  inS : Array Bool
  /-- node in `B_S` (a cut vertex or a vertex of `∂M`) -/
  isB : Array Bool
  /-- component (a representative node); meaningful where `inS` -/
  comp : Array Nat
  /-- node lies in a component that has a boundary vertex -/
  solved : Array Bool
  /-- cotan Laplacian of `S` -/
  lap : Mat
  /-- scale factor -/
  u : FloatArray
  /-- adjoint -/
  p : FloatArray
  /-- integrated outward normal derivative of `u` -/
  v : FloatArray
  /-- lumped boundary length -/
  ell : FloatArray
  /-- integrated outward normal derivative of the adjoint -/
  q : FloatArray
  /-- recovered normal derivative of `u` (`v/ℓ` away from the cut) -/
  d : FloatArray
  /-- recovered normal derivative of `p` -/
  e : FloatArray
  /-- gradient density `G` (before `w_r`, `β_r`) -/
  G : FloatArray
  /-- distortion energy `J` -/
  J : Float
  /-- area of the whole submesh -/
  area : Float
  /-- area of the solved components -/
  areaSolved : Float
  /-- number of triangles -/
  nTri : Nat
  /-- number of components -/
  nComp : Nat
  /-- number of components without boundary vertex -/
  nUnsolvedComp : Nat
  /-- number of nodes in them -/
  nUnsolvedNodes : Nat
  deriving Inhabited

private partial def findRoot (par : Array Nat) (x : Nat) : Nat :=
  let p := par[x]!
  if p == x then x else findRoot par p

/-- Flux recovery of SPEC §3.5: the normal derivatives from the integrated fluxes `vals` (`v` or
`q`) at the nodes of `B_S`. At the cut vertices, `smoothed c` solves
`(diag ℓ + Σ_σ κ_σ [[1,-1],[-1,1]]) d = vals`; everywhere else, and for `nodal`, `d = vals/ℓ`. -/
def recoverFlux (b : Base) (cm : CutMesh) (flux : FluxRule) (isB : Array Bool) (ell vals : FloatArray) :
    Except String FloatArray := do
  let N := cm.nodes
  let nv := cm.nv
  let nc := cm.nc
  let mut out := zeros N
  for i in [0:N] do
    if isB[i]! then out := out.set! i (normalDerivative (vals.get! i) (ell.get! i))
  match flux with
  | .nodal => return out
  | .smoothed c =>
    if nc == 0 then return out
    let mut m := zeros (nc * nc)
    let mut rhs := zeros nc
    for k in [0:nc] do
      m := m.set! (k * nc + k) (ell.get! (nv + k))
      rhs := rhs.set! k (fluxRhs (vals.get! (nv + k)))
    for sg in cm.segs do
      let f := b.topo.faces[sg.face]!
      let h2 := faceMeanSquareEdge (b.l f.a f.b) (b.l f.b f.c) (b.l f.c f.a)
      let kap := fluxCoupling c h2 sg.len
      m := m.set! (sg.p * nc + sg.p) (m.get! (sg.p * nc + sg.p) + kap)
      m := m.set! (sg.q * nc + sg.q) (m.get! (sg.q * nc + sg.q) + kap)
      m := m.set! (sg.p * nc + sg.q) (m.get! (sg.p * nc + sg.q) + fluxOffEntry kap)
      m := m.set! (sg.q * nc + sg.p) (m.get! (sg.q * nc + sg.p) + fluxOffEntry kap)
    match luFactor nc m with
    | none => throw "flux recovery: singular system"
    | some lu =>
      let d := lu.solve rhs
      for k in [0:nc] do out := out.set! (nv + k) (d.get! k)
      return out

/-- Solve region `r` for the case's energy. `aNode` is the weight `a` at every node (empty for
unscaled energies), `kf` the curvature density per face of `M`. -/
def computeRegion (b : Base) (cm : CutMesh) (en : Energy) (aNode : FloatArray)
    (kf : Array Float) (flux : FluxRule) (r : Nat) : Except String RegionData := do
  let N := cm.nodes
  let nv := cm.nv
  let tris := cm.tris.filter (·.region == r)
  let mut inS : Array Bool := Array.replicate N false
  let mut lap := Mat.zero N
  let mut g := zeros N
  let mut areaS := zeros N
  let mut par : Array Nat := Array.range N
  let mut area := fzero
  for t in tris do
    for n in [t.n0, t.n1, t.n2] do
      inS := inS.set! n true
      g := g.set! n (g.get! n + kf[t.parent]! * t.area / 3)
      areaS := areaS.set! n (areaS.get! n + t.area / 3)
    lap := lap.addTri t.n0 t.n1 t.n2 t.c0 t.c1 t.c2 1
    area := area + t.area
    let r0 := findRoot par t.n0
    let r1 := findRoot par t.n1
    par := par.set! r1 r0
    let r2 := findRoot par t.n2
    let r0' := findRoot par t.n0
    par := par.set! r2 r0'
  let isB : Array Bool := Array.ofFn (n := N) fun i =>
    inS[i.val]! && (i.val ≥ nv || b.topo.isBoundary[i.val]!)
  let comp : Array Nat := Array.ofFn (n := N) fun i => findRoot par i.val
  let mut hasB : Array Bool := Array.replicate N false
  for i in [0:N] do
    if isB[i]! then hasB := hasB.set! comp[i]! true
  let solved : Array Bool := Array.ofFn (n := N) fun i => inS[i.val]! && hasB[comp[i.val]!]!
  let mut roots : Array Nat := #[]
  let mut nUnsolvedComp := 0
  let mut nUnsolvedNodes := 0
  for i in [0:N] do
    if inS[i]! && comp[i]! == i then
      roots := roots.push i
      if !hasB[i]! then nUnsolvedComp := nUnsolvedComp + 1
    if inS[i]! && !solved[i]! then nUnsolvedNodes := nUnsolvedNodes + 1
  let mut areaSolved := fzero
  for t in tris do
    if solved[t.n0]! then areaSolved := areaSolved + t.area
  -- index sets
  let mut I : Array Nat := #[]
  let mut B : Array Nat := #[]
  for i in [0:N] do
    if solved[i]! then
      if isB[i]! then B := B.push i else I := I.push i
  -- u
  let gI := FloatArray.ofFn' I.size fun k => -(g.get! I[k]!)
  let some uI := lap.solveSub I gI | throw s!"region {r}: singular Yamabe system"
  let mut u := zeros N
  for k in [0:I.size] do u := u.set! I[k]! (uI.get! k)
  let lu := lap.mulVec u
  -- v
  let mut v := zeros N
  for i in B do v := v.set! i (normalFlux (lu.get! i) (g.get! i))
  -- boundary lengths
  let mut edgeCount : Std.HashMap (Nat × Nat) (Nat × Float) := {}
  for t in tris do
    for (x, y, l) in [(t.n0, t.n1, t.l01), (t.n1, t.n2, t.l12), (t.n2, t.n0, t.l20)] do
      let key := (min x y, max x y)
      let (cnt, _) := edgeCount.getD key (0, 0.0)
      edgeCount := edgeCount.insert key (cnt + 1, l)
  let mut bl : Array (List Float) := Array.replicate N []
  for ((x, y), (cnt, l)) in edgeCount.toList do
    if cnt == 1 then
      bl := (bl.set! x (l :: bl[x]!)).set! y (l :: bl[y]!)
  let mut ell := zeros N
  for i in B do ell := ell.set! i (dualLength bl[i]!)
  for i in B do
    if !(ell.get! i > 0) then throw s!"region {r}: node {i} on the boundary has no boundary edge"
  -- energy specific part
  let mut p := zeros N
  let mut q := zeros N
  let mut d := zeros N
  let mut e := zeros N
  let mut G := zeros N
  let mut J := fzero
  match en with
  | .dirichlet =>
    for k in [0:I.size] do J := J - uI.get! k * g.get! I[k]!
    d ← recoverFlux b cm flux isB ell v
    for i in B do G := G.set! i (gradDirichletD (d.get! i))
  | .scaledA _ | .scaledCurv _ _ =>
    let mut la := Mat.zero N
    for t in tris do
      la := la.addTri t.n0 t.n1 t.n2 t.c0 t.c1 t.c2
        (triangleMean (aNode.get! t.n0) (aNode.get! t.n1) (aNode.get! t.n2))
    let lau := la.mulVec u
    for i in [0:N] do J := J + u.get! i * lau.get! i
    let rI := FloatArray.ofFn' I.size fun k => adjointRhsScaledInterior (lau.get! I[k]!)
    let some pI := lap.solveSub I rI | throw s!"region {r}: singular adjoint system"
    for k in [0:I.size] do p := p.set! I[k]! (pI.get! k)
    let lp := lap.mulVec p
    for i in B do
      let ri := adjointRhsScaledBoundary (lau.get! i) (aNode.get! i) (v.get! i)
      q := q.set! i (adjointFlux (lp.get! i) ri)
    d ← recoverFlux b cm flux isB ell v
    e ← recoverFlux b cm flux isB ell q
    for i in B do G := G.set! i (gradScaledD (aNode.get! i) (d.get! i) (e.get! i))
  | .hencky =>
    -- b(u) = u², b'(u) = 2u, b(0) = 0
    for i in [0:N] do
      if solved[i]! then J := J + areaS.get! i * (u.get! i * u.get! i)
    let rI := FloatArray.ofFn' I.size fun k =>
      adjointRhsPointwise (areaS.get! I[k]!) (2 * u.get! I[k]!)
    let some pI := lap.solveSub I rI | throw s!"region {r}: singular adjoint system"
    for k in [0:I.size] do p := p.set! I[k]! (pI.get! k)
    let lp := lap.mulVec p
    for i in B do
      let ri := adjointRhsPointwise (areaS.get! i) (2 * u.get! i)
      q := q.set! i (adjointFlux (lp.get! i) ri)
    d ← recoverFlux b cm flux isB ell v
    e ← recoverFlux b cm flux isB ell q
    for i in B do G := G.set! i (gradPointwiseD 0 (d.get! i) (e.get! i))
  return { r, inS, isB, comp, solved, lap, u, p, v, ell, q, d, e, G, J, area, areaSolved
           nTri := tris.size, nComp := roots.size, nUnsolvedComp, nUnsolvedNodes }

end Ddg.Cuts
