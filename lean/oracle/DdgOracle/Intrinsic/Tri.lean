import Std.Data.HashMap
import DdgOracle.Topology
import DdgOracle.Dense
import DdgOracle.Hex
import DdgOracle.Formulas

/-!
# A signpost triangulation as a halfedge mesh (SPEC `ddg-intrinsic` §3, §4)

Halfedges `h = 3 f + k` of the input mesh; `next`, `twin` and `tail` are stored, so that flips can
produce self-edges, double edges and triangles that use a vertex twice. A halfedge without a twin is
a boundary halfedge (the surface lies on its left). Every halfedge has a length (the same on both
halfedges of an edge) and a signpost: the normalized angle of its direction at its tail.
-/

namespace Ddg.Intrinsic

/-- "No twin". -/
def noTwin : Nat := 4294967295

/-- The triangulation with its signposts. -/
structure Tri where
  /-- number of vertices -/
  nv : Nat
  /-- number of halfedges (three per face) -/
  nh : Nat
  /-- tail vertex of every halfedge -/
  tail : Array Nat
  /-- next halfedge in the face -/
  next : Array Nat
  /-- twin halfedge or `noTwin` -/
  twin : Array Nat
  /-- face of the halfedge -/
  face : Array Nat
  /-- one halfedge of every face (the first corner) -/
  fh : Array Nat
  /-- length of the edge -/
  len : Array Float
  /-- signpost at the tail -/
  sign : Array Float
  /-- edge id of the halfedge (a halfedge id of the edge, kept under flips) -/
  eid : Array Nat
  /-- angle sum of every vertex of the input mesh (never changes) -/
  theta : FloatArray
  /-- vertex lies on the boundary -/
  isB : Array Bool
  deriving Inhabited

namespace Tri

/-- Number of faces. -/
def nf (t : Tri) : Nat := t.nh / 3

/-- The previous halfedge in the face. -/
def prev (t : Tri) (h : Nat) : Nat := t.next[t.next[h]!]!

/-- Head vertex of a halfedge. -/
def head (t : Tri) (h : Nat) : Nat := t.tail[t.next[h]!]!

/-- Has a twin. -/
def hasTwin (t : Tri) (h : Nat) : Bool := t.twin[h]! != noTwin

/-- Corner angle at the tail of `h` in its face, from the lengths. -/
def corner (t : Tri) (h : Nat) : Float :=
  cornerAngle t.len[t.next[h]!]! t.len[t.prev h]! t.len[h]!

/-- The scale of the signposts at a vertex. -/
def scaleAt (t : Tri) (v : Nat) : Float := angleScale t.isB[v]! (t.theta.get! v)

/-- Number of edges (an edge is counted at its representative halfedge). -/
def numEdges (t : Tri) : Nat := Id.run do
  let mut n := 0
  for h in [0:t.nh] do
    if t.eid[h]! == h then n := n + 1
  return n

/-- The three halfedges of a face, starting at its first corner. -/
def faceHalfedges (t : Tri) (f : Nat) : Array Nat :=
  let h0 := t.fh[f]!
  #[h0, t.next[h0]!, t.next[t.next[h0]!]!]

/-- Total area. -/
def totalArea (t : Tri) : Float := Id.run do
  let mut a := fzero
  for f in [0:t.nf] do
    let hs := t.faceHalfedges f
    a := a + triArea t.len[hs[0]!]! t.len[hs[1]!]! t.len[hs[2]!]!
  return a

/-- Angle sum at every vertex, computed from the current lengths. -/
def angleSums (t : Tri) : FloatArray := Id.run do
  let mut s := zeros t.nv
  for h in [0:t.nh] do
    s := s.set! t.tail[h]! (s.get! t.tail[h]! + t.corner h)
  return s

end Tri

/-- Build the triangulation of an input mesh (SPEC Alg. 4): lengths from the edge table `len` (row-major,
`nv × nv`), signposts set. The first halfedge of an interior vertex is its outgoing halfedge with the
smallest id; the first halfedge of a boundary vertex is its outgoing halfedge without a twin (N2). -/
def initTri (nv : Nat) (faces : Array Face) (lenTab : FloatArray) : Except String Tri := do
  let nh := 3 * faces.size
  let mut tail : Array Nat := Array.replicate nh 0
  let mut next : Array Nat := Array.replicate nh 0
  let mut face : Array Nat := Array.replicate nh 0
  let mut len : Array Float := Array.replicate nh 0.0
  let mut dir : Std.HashMap (Nat × Nat) Nat := {}
  for f in [0:faces.size] do
    let vs := #[faces[f]!.a, faces[f]!.b, faces[f]!.c]
    for k in [0:3] do
      let h := 3 * f + k
      tail := tail.set! h vs[k]!
      next := next.set! h (3 * f + (k + 1) % 3)
      face := face.set! h f
      len := len.set! h (lenTab.get! (vs[k]! * nv + vs[(k + 1) % 3]!))
      dir := dir.insert (vs[k]!, vs[(k + 1) % 3]!) h
  let mut twin : Array Nat := Array.replicate nh noTwin
  let mut eid : Array Nat := Array.range nh
  for h in [0:nh] do
    let a := tail[h]!
    let b := tail[next[h]!]!
    match dir.get? (b, a) with
    | some t =>
      twin := twin.set! h t
      eid := eid.set! h (min h t)
    | none => pure ()
  let mut isB : Array Bool := Array.replicate nv false
  for h in [0:nh] do
    if twin[h]! == noTwin then isB := isB.set! tail[h]! true
  let t0 : Tri := { nv, nh, tail, next, twin, face, fh := Array.ofFn (n := faces.size) fun f => 3 * f.val
                    len, sign := Array.replicate nh 0.0, eid, theta := zeros nv, isB }
  -- angle sums
  let theta := t0.angleSums
  let mut t := { t0 with theta }
  -- signposts
  let mut sign : Array Float := Array.replicate nh 0.0
  for v in [0:nv] do
    -- first halfedge
    let mut first := nh
    if t.isB[v]! then
      for h in [0:nh] do
        if first == nh && t.tail[h]! == v && t.twin[h]! == noTwin then first := h
    else
      for h in [0:nh] do
        if first == nh && t.tail[h]! == v then first := h
    if first == nh then throw s!"vertex {v} is not used by any face"
    let mut cur := first
    let mut phi := fzero
    sign := sign.set! cur phi
    let mut go := true
    let mut guard := 0
    while go && guard ≤ nh do
      guard := guard + 1
      let p := t.prev cur
      let nxt := t.twin[p]!
      if nxt == noTwin || nxt == first then go := false
      else
        phi := signpostUpdate t.isB[v]! phi (t.corner cur) (theta.get! v)
        sign := sign.set! nxt phi
        cur := nxt
  t := { t with sign }
  return t

end Ddg.Intrinsic
