import DdgOracle.Regions.Case

/-!
# Crossings, triple points and the cut mesh of N regions (SPEC §8.2)

*Geometric ids*: an original vertex `v < nv`; a cut point `p`: `nv + p`; a triple point `m`:
`nv + nP + m`. *Cut-mesh nodes* (copies): a vertex `v`; a cut point on `[a, b]` has `nv + 2p`
(region of `a`) and `nv + 2p + 1` (region of `b`); a triple point has `nv + 2 nP + 3m + k` in the
region of the corner `c_k`.
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- A triple point. -/
structure Triple where
  /-- the face -/
  face : Nat
  /-- corners `c_0, c_1, c_2` (the face's `a, b, c`) -/
  c0 : Nat
  /-- see `c0` -/
  c1 : Nat
  /-- see `c0` -/
  c2 : Nat
  /-- cut points on the edges `c_0 c_1`, `c_1 c_2`, `c_2 c_0` -/
  p0 : Nat
  /-- see `p0` -/
  p1 : Nat
  /-- see `p0` -/
  p2 : Nat
  /-- barycentric coordinates in the face (corners `a, b, c`) -/
  b0 : Float
  /-- see `b0` -/
  b1 : Float
  /-- see `b0` -/
  b2 : Float
  deriving Inhabited

/-- A segment of the network between two geometric nodes. -/
structure NSeg where
  /-- first geometric id (a cut point) -/
  g1 : Nat
  /-- second geometric id (a cut point or a triple point) -/
  g2 : Nat
  /-- the face -/
  face : Nat
  /-- length -/
  len : Float
  /-- region of the smaller end of the edge of `g1` -/
  regA : Nat
  /-- region of the larger end of the edge of `g1` -/
  regB : Nat
  /-- the second node is a triple point -/
  toTriple : Bool
  /-- index of the triangle `(a, p, q)` for a face with two regions -/
  tri : Nat
  /-- the lone vertex of a face with two regions -/
  lone : Nat
  deriving Inhabited

/-- The network: regions, cut points, triple points, cut mesh. -/
structure Net where
  /-- vertices -/
  nv : Nat
  /-- regions -/
  nreg : Nat
  /-- region of every vertex -/
  region : Array Nat
  /-- cut points in the order of their edges -/
  cps : Array Crossing
  /-- triple points in the order of their faces -/
  triples : Array Triple
  /-- triangles of the cut mesh -/
  tris : Array Tri
  /-- segments in the order of the faces -/
  segs : Array NSeg
  /-- vertices of crossed edges -/
  adj : Array Bool
  /-- quads whose diagonals tie -/
  nearTies : Nat
  deriving Inhabited

namespace Net

/-- Number of cut points. -/
def nP (n : Net) : Nat := n.cps.size

/-- Number of triple points. -/
def nM (n : Net) : Nat := n.triples.size

/-- Number of nodes of the cut mesh. -/
def nodes (n : Net) : Nat := n.nv + 2 * n.cps.size + 3 * n.triples.size

/-- The copy of a geometric node in region `X`. -/
def copyIdx (n : Net) (g : Nat) (X : Nat) : Nat :=
  if g < n.nv then g
  else if g < n.nv + n.nP then
    let p := g - n.nv
    n.nv + 2 * p + (if X == n.region[n.cps[p]!.b]! then 1 else 0)
  else
    let m := g - n.nv - n.nP
    let t := n.triples[m]!
    let k := if X == n.region[t.c0]! then 0 else if X == n.region[t.c1]! then 1 else 2
    n.nv + 2 * n.nP + 3 * m + k

end Net

/-- Region of a vertex: the index of the smallest `φ_k` (ties: the lowest `k`). -/
def argminRegion (phis : Array FloatArray) (i : Nat) : Nat := Id.run do
  let mut best := 0
  let mut bv := (phis[0]!).get! i
  for k in [1:phis.size] do
    if (phis[k]!).get! i < bv then
      bv := (phis[k]!).get! i
      best := k
  return best

/-- Build the network of the (scaled) functions `phis`. -/
def buildNet (b : Base) (faces : Array Face) (phis : Array FloatArray) (delta : Float) :
    Except String Net := do
  let nv := b.nv
  let N := phis.size
  let region : Array Nat := Array.ofFn (n := nv) fun i => argminRegion phis i.val
  -- crossings (Eq. 31)
  let mut cps : Array Crossing := #[]
  let mut idx : Std.HashMap (Nat × Nat) Nat := {}
  let mut adj : Array Bool := Array.replicate nv false
  for a in [0:nv] do
    for c in b.topo.neighbors[a]! do
      if a < c && region[a]! != region[c]! then
        let i := region[a]!
        let j := region[c]!
        let tRaw := crossingParamRegions ((phis[i]!).get! a) ((phis[j]!).get! a)
          ((phis[i]!).get! c) ((phis[j]!).get! c)
        let t := clampCrossing delta tRaw
        idx := idx.insert (a, c) cps.size
        cps := cps.push { a, b := c, tRaw, t, moved := crossingMoved delta tRaw, len := b.l a c }
        adj := (adj.set! a true).set! c true
  let nP := cps.size
  -- triple points: first pass, to know their number
  let mut triples : Array Triple := #[]
  for fi in [0:faces.size] do
    let f := faces[fi]!
    let r0 := region[f.a]!
    let r1 := region[f.b]!
    let r2 := region[f.c]!
    if r0 != r1 && r1 != r2 && r2 != r0 then
      let cs := #[f.a, f.b, f.c]
      let pk := fun (k : Nat) => idx.getD (min cs[k % 3]! cs[(k + 1) % 3]!, max cs[k % 3]! cs[(k + 1) % 3]!) 0
      -- parameter of the cut point on the edge c_k c_{k+1}, from c_k
      let sk := fun (k : Nat) =>
        let cr := cps[pk k]!
        if cs[k % 3]! < cs[(k + 1) % 3]! then cr.t else 1 - cr.t
      let coord := fun (k : Nat) => triplePointCoord (sk k) (sk (k + 2))
      triples := triples.push { face := fi, c0 := f.a, c1 := f.b, c2 := f.c
                                p0 := pk 0, p1 := pk 1, p2 := pk 2
                                b0 := coord 0, b1 := coord 1, b2 := coord 2 }
  let nM := triples.size
  let net0 : Net := { nv, nreg := N, region, cps, triples, tris := #[], segs := #[], adj, nearTies := 0 }
  -- barycentric coordinates of a geometric node in a face
  let bary := fun (f : Face) (g : Nat) =>
    let unit := fun (v : Nat) =>
      if v == f.a then (1.0, 0.0, 0.0) else if v == f.b then (0.0, 1.0, 0.0) else (0.0, 0.0, 1.0)
    if g < nv then unit g
    else if g < nv + nP then
      let cr := cps[g - nv]!
      let (xa, xb, xc) := unit cr.a
      let (ya, yb, yc) := unit cr.b
      ((1 - cr.t) * xa + cr.t * ya, (1 - cr.t) * xb + cr.t * yb, (1 - cr.t) * xc + cr.t * yc)
    else
      let t := triples[g - nv - nP]!
      (t.b0, t.b1, t.b2)
  let mut tris : Array Tri := #[]
  let mut segs : Array NSeg := #[]
  let mut nearTies := 0
  let mut mi := 0
  for fi in [0:faces.size] do
    let f := faces[fi]!
    let lab := b.l f.a f.b
    let lbc := b.l f.b f.c
    let lca := b.l f.c f.a
    let af := b.faceArea[fi]!
    let mk := fun (g0 g1 g2 reg : Nat) =>
      let p0 := bary f g0
      let p1 := bary f g1
      let p2 := bary f g2
      let dist := fun (p q : Float × Float × Float) =>
        baryLength lab lbc lca (p.1 - q.1) (p.2.1 - q.2.1) (p.2.2 - q.2.2)
      let l01 := dist p0 p1
      let l12 := dist p1 p2
      let l20 := dist p2 p0
      let d := det3 p0.1 p0.2.1 p0.2.2 p1.1 p1.2.1 p1.2.2 p2.1 p2.2.1 p2.2.2
      let area := subTriangleArea d af
      let sub := fun (p q : Float × Float × Float) => (p.1 - q.1, p.2.1 - q.2.1, p.2.2 - q.2.2)
      let inner := fun (o x y : Float × Float × Float) =>
        let u := sub x o
        let w := sub y o
        baryInner lab lbc lca u.1 u.2.1 u.2.2 w.1 w.2.1 w.2.2
      let uncut := g0 == f.a && g1 == f.b && g2 == f.c
      let (k0, k1, k2) := b.cots[fi]!
      ({ n0 := net0.copyIdx g0 reg, n1 := net0.copyIdx g1 reg, n2 := net0.copyIdx g2 reg
         region := reg, parent := fi, l01, l12, l20
         area := if uncut then af else area
         c0 := if uncut then k0 else cotFromInner (inner p0 p1 p2) area
         c1 := if uncut then k1 else cotFromInner (inner p1 p2 p0) area
         c2 := if uncut then k2 else cotFromInner (inner p2 p0 p1) area } : Tri)
    let vs := #[f.a, f.b, f.c]
    let rg := fun (k : Nat) => region[vs[k % 3]!]!
    let cpOf := fun (x y : Nat) => nv + idx.getD (min x y, max x y) 0
    if rg 0 == rg 1 && rg 1 == rg 2 then
      tris := tris.push (mk f.a f.b f.c (rg 0))
    else if rg 0 != rg 1 && rg 1 != rg 2 && rg 2 != rg 0 then
      -- three regions: triple point
      let g := nv + nP + mi
      mi := mi + 1
      let p := fun (k : Nat) => cpOf vs[k % 3]! vs[(k + 1) % 3]!
      for k in [0:3] do
        let reg := rg k
        tris := (tris.push (mk vs[k]! (p k) g reg)).push (mk vs[k]! g (p (k + 2)) reg)
      for k in [0:3] do
        let cr := cps[p k - nv]!
        let pb := bary f (p k)
        let mb := bary f g
        let len := baryLength lab lbc lca (pb.1 - mb.1) (pb.2.1 - mb.2.1) (pb.2.2 - mb.2.2)
        segs := segs.push { g1 := p k, g2 := g, face := fi, len, regA := region[cr.a]!, regB := region[cr.b]!
                            toTriple := true, tri := 0, lone := nv }
    else
      let mut k := 0
      for kk in [0:3] do
        if rg kk != rg (kk + 1) && rg kk != rg (kk + 2) then k := kk
      let a := vs[k]!
      let bb := vs[(k + 1) % 3]!
      let cc := vs[(k + 2) % 3]!
      let p := cpOf a bb
      let q := cpOf a cc
      let t1 := mk a p q (rg k)
      tris := tris.push t1
      let cr := cps[p - nv]!
      segs := segs.push { g1 := p, g2 := q, face := fi, len := t1.l12, regA := region[cr.a]!, regB := region[cr.b]!
                          toTriple := false, tri := tris.size - 1, lone := a }
      let reg2 := rg (k + 1)
      let len2 := fun (x y : Nat) =>
        let px := bary f x
        let py := bary f y
        baryLength2 lab lbc lca (px.1 - py.1) (px.2.1 - py.2.1) (px.2.2 - py.2.2)
      let dpc2 := len2 p cc
      let dbq2 := len2 bb q
      if (dpc2 - dbq2).abs ≤ 1e-9 * (dpc2 + dbq2) then nearTies := nearTies + 1
      if dpc2 ≤ dbq2 then
        tris := (tris.push (mk p bb cc reg2)).push (mk p cc q reg2)
      else
        tris := (tris.push (mk p bb q reg2)).push (mk bb cc q reg2)
  for t in tris do
    if !(t.area > 0) then
      throw s!"sub-triangle of face {t.parent} has no area (a crossing on a vertex; delta = 0?)"
  return { net0 with tris, segs, nearTies }

end Ddg.Regions
