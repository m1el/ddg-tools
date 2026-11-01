import DdgOracle.Dense
import DdgOracle.Hex
import DdgOracle.Generators
import DdgOracle.RepulsiveCurves.Basic

/-!
# Repulsive surfaces: triangle meshes, topology, validation

Plain `Float` code written from `crates/ddg-repulsive-surfaces/SPEC.md` §2 for the `rs` oracle. It
shares nothing with the Rust crate. Only the 3-vector type of the curves oracle is reused.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

/-- A triangle `(a, b, c)` in the caller's vertex order. -/
structure Tri where
  /-- first vertex -/
  a : Nat
  /-- second vertex -/
  b : Nat
  /-- third vertex -/
  c : Nat
  deriving Inhabited, BEq

/-- Vertex `t ∈ {0,1,2}` of the triangle. -/
def Tri.get (f : Tri) (t : Nat) : Nat := if t == 0 then f.a else if t == 1 then f.b else f.c

/-- Triangle mesh: positions and faces. -/
structure Mesh where
  /-- positions -/
  pos : Array V3
  /-- faces -/
  faces : Array Tri
  deriving Inhabited

/-- Number of vertices. -/
def Mesh.nv (m : Mesh) : Nat := m.pos.size
/-- Number of faces. -/
def Mesh.nf (m : Mesh) : Nat := m.faces.size

/-- The exponent pair `(α, β)` (SPEC §2.2). -/
structure Ex where
  /-- `α` -/
  alpha : Float
  /-- `β` -/
  beta : Float
  deriving Inhabited

/-- `s = (β − 2)/α` (SPEC §2.2). -/
def Ex.s (e : Ex) : Float := (e.beta - 2.0) / e.alpha

/-- `α` finite and `> 1`, `α + 2 < β < 2α + 2` (SPEC §2.2). -/
def Ex.check (e : Ex) : Except String Unit := do
  if !(e.alpha.isFinite && e.beta.isFinite) then throw "InvalidExponents: both exponents must be finite"
  if !(e.alpha > 1.0) then throw "InvalidExponents: alpha must be greater than 1"
  if !(e.beta > e.alpha + 2.0) then throw "InvalidExponents: beta must be greater than alpha + 2"
  if !(e.beta < 2.0 * e.alpha + 2.0) then throw "InvalidExponents: beta must be smaller than 2 alpha + 2"

/-- Per-face quantities (SPEC §2.1). -/
structure FaceGeo where
  /-- area `a_S = ½|n_S|` -/
  area : Float
  /-- unit normal `N_S` -/
  n : V3
  /-- barycenter `X_S` -/
  x : V3
  deriving Inhabited

/-- Face quantities at positions `pos`. -/
def faceGeo (faces : Array Tri) (pos : Array V3) : Array FaceGeo :=
  faces.map fun f =>
    let pi := pos[f.a]!
    let pj := pos[f.b]!
    let pk := pos[f.c]!
    let nn := (pj - pi).cross (pk - pi)
    let l := nn.norm
    { area := 0.5 * l, n := (1.0 / l) * nn, x := (1.0 / 3.0 : Float) * (pi + pj + pk) }

/-- The edge vector opposite vertex slot `t` of a face: `e_jk, e_ki, e_ij` for `t = 0, 1, 2`. -/
def oppEdge (pos : Array V3) (f : Tri) (t : Nat) : V3 :=
  pos[f.get ((t + 2) % 3)]! - pos[f.get ((t + 1) % 3)]!

/-! ## Topology -/

/-- Edge with its incident faces: `(lo, hi, [(face, traversed lo→hi)])`, sorted by `(lo, hi)`. -/
def edgeTable (faces : Array Tri) : Array (Nat × Nat × Array (Nat × Bool)) := Id.run do
  let mut raw : Array (Nat × Nat × Nat × Bool) := #[]
  for fi in [0:faces.size] do
    let f := faces[fi]!
    for t in [0:3] do
      let u := f.get t
      let v := f.get ((t + 1) % 3)
      raw := raw.push (min u v, max u v, fi, u < v)
  let sorted := raw.qsort fun x y => x.1 < y.1 || (x.1 == y.1 && (x.2.1 < y.2.1 || (x.2.1 == y.2.1 && x.2.2.1 < y.2.2.1)))
  let mut out : Array (Nat × Nat × Array (Nat × Bool)) := #[]
  for (lo, hi, fi, fw) in sorted do
    if out.size > 0 && out.back!.1 == lo && out.back!.2.1 == hi then
      let (a, b, l) := out.back!
      out := out.pop.push (a, b, l.push (fi, fw))
    else out := out.push (lo, hi, #[(fi, fw)])
  return out

/-- Validated topology of a mesh. -/
structure Topo where
  /-- number of vertices -/
  nv : Nat
  /-- faces -/
  faces : Array Tri
  /-- edge table -/
  edges : Array (Nat × Nat × Array (Nat × Bool))
  /-- component of each face (components numbered by smallest face index) -/
  faceComp : Array Nat
  /-- component of each vertex -/
  vertComp : Array Nat
  /-- number of components -/
  ncomp : Nat
  /-- boundary vertices -/
  isBoundary : Array Bool

/-- Validate the connectivity (SPEC §2.1). Errors are named like the Rust variants. -/
def buildTopo (nv : Nat) (faces : Array Tri) : Except String Topo := do
  for fi in [0:faces.size] do
    let f := faces[fi]!
    for t in [0:3] do
      if f.get t ≥ nv then throw s!"VertexOutOfRange: face {fi} = ({f.a}, {f.b}, {f.c}), {nv} vertices"
    if f.a == f.b || f.b == f.c || f.a == f.c then
      throw s!"DegenerateFace: face {fi} = ({f.a}, {f.b}, {f.c}) repeats a vertex"
  let mut used : Array Bool := Array.replicate nv false
  for f in faces do
    for t in [0:3] do used := used.set! (f.get t) true
  for v in [0:nv] do
    if !used[v]! then throw s!"UnreferencedVertex: vertex {v}"
  let edges := edgeTable faces
  for (a, b, l) in edges do
    if l.size > 2 then
      let ids : String := ", ".intercalate (l.toList.map fun x => toString x.1)
      throw s!"NonManifoldEdge: edge ({a}, {b}) is shared by {l.size} faces: {ids}"
  -- face adjacency
  let mut adj : Array (Array Nat) := Array.replicate faces.size #[]
  for (_, _, l) in edges do
    if l.size == 2 then
      adj := adj.modify l[0]!.1 (·.push l[1]!.1)
      adj := adj.modify l[1]!.1 (·.push l[0]!.1)
  -- components by flood fill
  let mut comp : Array Nat := Array.replicate faces.size faces.size
  let mut ncomp := 0
  for s in [0:faces.size] do
    if comp[s]! == faces.size then
      comp := comp.set! s ncomp
      let mut stack := #[s]
      while stack.size > 0 do
        let x := stack.back!
        stack := stack.pop
        for y in adj[x]! do
          if comp[y]! == faces.size then
            comp := comp.set! y ncomp
            stack := stack.push y
      ncomp := ncomp + 1
  -- vertex stars must be one fan
  let mut star : Array (Array Nat) := Array.replicate nv #[]
  for fi in [0:faces.size] do
    for t in [0:3] do star := star.modify (faces[fi]!.get t) (·.push fi)
  for v in [0:nv] do
    let fs := star[v]!
    let mut seen : Array Bool := Array.replicate fs.size false
    let mut count := 0
    for s in [0:fs.size] do
      if !seen[s]! then
        count := count + 1
        seen := seen.set! s true
        let mut stack := #[s]
        while stack.size > 0 do
          let x := stack.back!
          stack := stack.pop
          for y in [0:fs.size] do
            if !seen[y]! then
              -- faces share an edge through v iff they share a second vertex
              let fx := faces[fs[x]!]!
              let fy := faces[fs[y]!]!
              let mut shared := 0
              for p in [0:3] do
                for q in [0:3] do
                  if fx.get p == fy.get q && fx.get p != v then shared := shared + 1
              if shared > 0 then
                seen := seen.set! y true
                stack := stack.push y
    if count > 1 then throw s!"NonManifoldVertex: the faces around vertex {v} form {count} fans"
  let mut vcomp : Array Nat := Array.replicate nv 0
  for fi in [0:faces.size] do
    for t in [0:3] do vcomp := vcomp.set! (faces[fi]!.get t) comp[fi]!
  let mut bnd : Array Bool := Array.replicate nv false
  for (a, b, l) in edges do
    if l.size == 1 then
      bnd := (bnd.set! a true).set! b true
  return { nv, faces, edges, faceComp := comp, vertComp := vcomp, ncomp, isBoundary := bnd }

/-- Positions: finite and no zero-area face (SPEC §2.3). -/
def checkPositions (t : Topo) (pos : Array V3) : Except String Unit := do
  if pos.size != t.nv then throw s!"PositionCount: {pos.size} positions for {t.nv} vertices"
  for v in [0:pos.size] do
    if !pos[v]!.finite then throw s!"NonFinitePosition: vertex {v}"
  let g := faceGeo t.faces pos
  for fi in [0:g.size] do
    if !(g[fi]!.area > 0.0 && g[fi]!.area.isFinite) then
      throw s!"DegenerateFace: face {fi} has zero area"

/-- Topology and positions together. -/
def Mesh.validate (m : Mesh) : Except String Topo := do
  let t ← buildTopo m.nv m.faces
  checkPositions t m.pos
  return t

/-- Mean edge length. -/
def meanEdge (t : Topo) (pos : Array V3) : Float := Id.run do
  let mut s := 0.0
  for (a, b, _) in t.edges do s := s + (pos[a]! - pos[b]!).norm
  return s / t.edges.size.toFloat

/-- Vertex (barycentric dual) areas `a_v = ⅓ Σ a_S` (SPEC §2.1). -/
def vertexAreas (t : Topo) (g : Array FaceGeo) : Array Float := Id.run do
  let mut a : Array Float := Array.replicate t.nv 0.0
  for fi in [0:t.faces.size] do
    for k in [0:3] do
      a := a.modify (t.faces[fi]!.get k) (· + g[fi]!.area / 3.0)
  return a

/-- First interior edge of component `c` whose two faces traverse it in the same direction. -/
def inconsistentEdge (t : Topo) (c : Nat) : Option (Nat × Nat × Nat × Nat) := Id.run do
  for (a, b, l) in t.edges do
    if l.size == 2 && t.faceComp[l[0]!.1]! == c && l[0]!.2 == l[1]!.2 then
      return some (a, b, l[0]!.1, l[1]!.1)
  return none

end Ddg.RepulsiveSurfaces
