import Std.Data.HashMap
import DdgOracle.Cuts.Base

/-!
# Regions, crossings and the cut mesh (SPEC §3.2, §3.3)

Nodes of the cut mesh: the original vertices `0 … nv-1`, then one *cut vertex* per crossed edge,
node `nv + k` for the `k`-th crossed edge in increasing lexicographic order of `(a, b)`.
-/

namespace Ddg.Cuts

/-- A crossed edge `[a, b]`, `a < b`. -/
structure Crossing where
  /-- smaller endpoint -/
  a : Nat
  /-- larger endpoint -/
  b : Nat
  /-- `φ_a/(φ_a − φ_b)` before clamping -/
  tRaw : Float
  /-- the parameter used, from `a` -/
  t : Float
  /-- `t` was moved into `[δ, 1 − δ]` -/
  moved : Bool
  /-- length of the edge (scaled) -/
  len : Float
  deriving Inhabited

/-- A triangle of the cut mesh. Corner `k` is at node `nk`; `l01` is the length between `n0` and
`n1`, and so on. -/
structure Tri where
  /-- node 0 -/
  n0 : Nat
  /-- node 1 -/
  n1 : Nat
  /-- node 2 -/
  n2 : Nat
  /-- region (0 or 1) -/
  region : Nat
  /-- face of `M` that contains it -/
  parent : Nat
  /-- length `|n0 n1|` -/
  l01 : Float
  /-- length `|n1 n2|` -/
  l12 : Float
  /-- length `|n2 n0|` -/
  l20 : Float
  /-- area -/
  area : Float
  /-- cotangent of the corner at `n0` -/
  c0 : Float
  /-- cotangent of the corner at `n1` -/
  c1 : Float
  /-- cotangent of the corner at `n2` -/
  c2 : Float
  deriving Inhabited

/-- A cut segment: the two cut vertices (indices `k`, not nodes) of a crossed face. -/
structure Seg where
  /-- first cut vertex (on the edge `a b` of the face, see `Tri` order `a p q`) -/
  p : Nat
  /-- second cut vertex (on the edge `a c`) -/
  q : Nat
  /-- the face -/
  face : Nat
  /-- length -/
  len : Float
  /-- index of the triangle `(a, p, q)` -/
  tri : Nat
  /-- the lone vertex `a` -/
  lone : Nat
  deriving Inhabited

/-- The cut mesh of a level set function. -/
structure CutMesh where
  /-- number of original vertices -/
  nv : Nat
  /-- region of every original vertex -/
  region : Array Nat
  /-- crossed edges -/
  crossings : Array Crossing
  /-- triangles -/
  tris : Array Tri
  /-- segments of the cut -/
  segs : Array Seg
  /-- endpoints of crossed edges -/
  adj : Array Bool
  /-- quads whose two diagonals agree to a relative `10⁻⁹`: the choice of diagonal is fragile -/
  nearTies : Nat
  deriving Inhabited

namespace CutMesh

/-- Number of cut vertices. -/
def nc (cm : CutMesh) : Nat := cm.crossings.size

/-- Number of nodes. -/
def nodes (cm : CutMesh) : Nat := cm.nv + cm.crossings.size

end CutMesh

/-- Region of a value of `φ`: `0` iff `φ ≤ 0`. -/
def regionOf (x : Float) : Nat := if x ≤ 0 then 0 else 1

/-- Barycentric coordinates of a node within a face. -/
private def bary (cm : Array Crossing) (nv : Nat) (f : Face) (node : Nat) : Float × Float × Float :=
  let unit := fun (v : Nat) =>
    if v == f.a then (1.0, 0.0, 0.0) else if v == f.b then (0.0, 1.0, 0.0) else (0.0, 0.0, 1.0)
  if node < nv then unit node
  else
    let cr := cm[node - nv]!
    let (xa, xb, xc) := unit cr.a
    let (ya, yb, yc) := unit cr.b
    ((1 - cr.t) * xa + cr.t * ya, (1 - cr.t) * xb + cr.t * yb, (1 - cr.t) * xc + cr.t * yc)

/-- Build the regions, crossings and the cut mesh for the (scaled) level set `phi`. -/
def buildCutMesh (b : Base) (faces : Array Face) (phi : FloatArray) (delta : Float) : CutMesh :=
  Id.run do
  let nv := b.nv
  let region : Array Nat := Array.ofFn (n := nv) fun i => regionOf (phi.get! i.val)
  let mut crossings : Array Crossing := #[]
  let mut idx : Std.HashMap (Nat × Nat) Nat := {}
  let mut adj : Array Bool := Array.replicate nv false
  for a in [0:nv] do
    for c in b.topo.neighbors[a]! do
      if a < c && region[a]! != region[c]! then
        let tRaw := crossingParam (phi.get! a) (phi.get! c)
        let t := clampCrossing delta tRaw
        idx := idx.insert (a, c) crossings.size
        crossings := crossings.push
          { a, b := c, tRaw, t, moved := crossingMoved delta tRaw, len := b.l a c }
        adj := (adj.set! a true).set! c true
  let mut tris : Array Tri := #[]
  let mut segs : Array Seg := #[]
  let mut nearTies := 0
  for fi in [0:faces.size] do
    let f := faces[fi]!
    let lab := b.l f.a f.b
    let lbc := b.l f.b f.c
    let lca := b.l f.c f.a
    let af := b.faceArea[fi]!
    let mk := fun (n0 n1 n2 reg : Nat) =>
      let p0 := bary crossings nv f n0
      let p1 := bary crossings nv f n1
      let p2 := bary crossings nv f n2
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
      -- uncut faces keep the cotangents and the area of `M`
      let uncut := n0 == f.a && n1 == f.b && n2 == f.c
      let (k0, k1, k2) := b.cots[fi]!
      ({ n0, n1, n2, region := reg, parent := fi, l01, l12, l20
         area := if uncut then af else area
         c0 := if uncut then k0 else cotFromInner (inner p0 p1 p2) area
         c1 := if uncut then k1 else cotFromInner (inner p1 p2 p0) area
         c2 := if uncut then k2 else cotFromInner (inner p2 p0 p1) area } : Tri)
    let vs := #[f.a, f.b, f.c]
    let rg := fun (k : Nat) => region[vs[k % 3]!]!
    if rg 0 == rg 1 && rg 1 == rg 2 then
      tris := tris.push (mk f.a f.b f.c (rg 0))
    else
      -- the alone vertex: its region differs from the other two
      let mut k := 0
      for kk in [0:3] do
        if rg kk != rg (kk + 1) && rg kk != rg (kk + 2) then k := kk
      let a := vs[k]!
      let bb := vs[(k + 1) % 3]!
      let cc := vs[(k + 2) % 3]!
      let edgeNode := fun (x y : Nat) => nv + (idx.getD (min x y, max x y) 0)
      let p := edgeNode a bb
      let q := edgeNode a cc
      let t1 := mk a p q (rg k)
      tris := tris.push t1
      segs := segs.push { p := p - nv, q := q - nv, face := fi, len := t1.l12, tri := tris.size - 1, lone := a }
      let reg2 := rg (k + 1)
      -- diagonals p–c and b–q of the quad (p, b, c, q)
      let len2 := fun (x y : Nat) =>
        let px := bary crossings nv f x
        let py := bary crossings nv f y
        baryLength2 lab lbc lca (px.1 - py.1) (px.2.1 - py.2.1) (px.2.2 - py.2.2)
      let dpc2 := len2 p cc
      let dbq2 := len2 bb q
      if (dpc2 - dbq2).abs ≤ 1e-9 * (dpc2 + dbq2) then nearTies := nearTies + 1
      if dpc2 ≤ dbq2 then
        tris := (tris.push (mk p bb cc reg2)).push (mk p cc q reg2)
      else
        tris := (tris.push (mk p bb q reg2)).push (mk bb cc q reg2)
  return { nv, region, crossings, tris, segs, adj, nearTies }

end Ddg.Cuts
