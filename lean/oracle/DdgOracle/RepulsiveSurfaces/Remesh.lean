import DdgOracle.RepulsiveSurfaces.Step

/-!
# Remeshing operations (SPEC §10): one split, one collapse, one flip

Operation semantics only. Numbering conventions (so that the output is canonical):

* **split `a b`**: the midpoint is vertex `nv` (appended). For each incident face, in increasing face
  index, with the directed edge `p → q` and opposite vertex `c` (cyclic order): the face keeps its
  index and becomes `(p, m, c)`, and `(m, q, c)` is appended. Parents: the appended face's parent is
  the face it was cut from.
* **collapse `a b`**: the vertex `min(a,b)` survives, `max(a,b)` is removed, vertices above it shift
  down by one. Faces containing both ends are deleted, the others keep their relative order. The
  surviving vertex goes to the pinned endpoint if exactly one is pinned, else to the boundary
  endpoint if exactly one is on the boundary, else to the midpoint. Both pinned: `BothPinned`.
* **flip `a b`**: faces `F1 < F2` on the edge, `F1` traverses `p → q`, `F2` traverses `q → p`;
  `c` opposite in `F1`, `d` opposite in `F2`. `F1 := (p, d, c)`, `F2 := (d, q, c)`.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

/-- Output of an operation. -/
structure RemeshOut where
  /-- positions -/
  pos : Array V3
  /-- faces -/
  faces : Array Tri
  /-- old face of every new face -/
  parent : Array Nat
  /-- old vertices of every new vertex (for collapse: the survivor maps to both) -/
  vmap : Array (Array Nat)

/-- Angle between two vectors in `[0, π]`. -/
def angleBetween (u v : V3) : Float := Float.atan2 (u.cross v).norm (u.dot v)

/-- The foldover threshold `π − 2` (SPEC §10). -/
def foldoverLimit : Float := 3.141592653589793 - 2.0

/-- Edge `a b`: the table entry. -/
def findEdge (t : Topo) (a b : Nat) : Except String (Array (Nat × Bool)) := do
  let lo := min a b
  let hi := max a b
  match t.edges.find? (fun e => e.1 == lo && e.2.1 == hi) with
  | some (_, _, l) => pure l
  | none => throw s!"NoSuchEdge: ({a}, {b}) is not an edge"

private def orient (f : Tri) (a b : Nat) : Option (Nat × Nat × Nat) := Id.run do
  -- directed edge p→q of f among {a,b} and the opposite vertex c
  for k in [0:3] do
    let p := f.get k
    let q := f.get ((k + 1) % 3)
    let c := f.get ((k + 2) % 3)
    if (p == a && q == b) || (p == b && q == a) then return some (p, q, c)
  return none

/-- Edge split at the midpoint. -/
def splitEdge (m : Mesh) (t : Topo) (a b : Nat) : Except String RemeshOut := do
  let inc ← findEdge t a b
  let mid := (0.5 : Float) * (m.pos[a]! + m.pos[b]!)
  let nv := m.nv
  let mut faces := m.faces
  let mut parent := Array.range m.nf
  let sorted := inc.qsort fun x y => x.1 < y.1
  for (fi, _) in sorted do
    let some (p, q, c) := orient m.faces[fi]! a b | throw "internal: edge not in face"
    faces := faces.set! fi ⟨p, nv, c⟩
    faces := faces.push ⟨nv, q, c⟩
    parent := parent.push fi
  let vmap := (Array.range nv).map (fun v => #[v]) |>.push #[a, b]
  return { pos := m.pos.push mid, faces, parent, vmap }

/-- Edge collapse (link condition always; foldover and empty-interior checks only if asked). -/
def collapseEdge (m : Mesh) (t : Topo) (a b : Nat) (pinned : Array Nat) (foldover : Bool) :
    Except String RemeshOut := do
  let inc ← findEdge t a b
  let pa := pinned.contains a
  let pb := pinned.contains b
  if pa && pb then throw s!"BothPinned: both ends of ({a}, {b}) are pinned"
  -- link condition
  let nbrs := fun (v : Nat) => Id.run do
    let mut r : Array Nat := #[]
    for f in m.faces do
      for k in [0:3] do
        if f.get k == v then
          for j in [0:3] do
            if f.get j != v && !r.contains (f.get j) then r := r.push (f.get j)
    return r
  let na := nbrs a
  let nb := nbrs b
  let common := na.filter (nb.contains ·)
  if common.size != inc.size then
    throw s!"LinkCondition: the ends of ({a}, {b}) have {common.size} common neighbours but the edge has {inc.size} faces"
  if t.isBoundary[a]! && t.isBoundary[b]! && inc.size == 2 then
    throw s!"LinkCondition: ({a}, {b}) is an interior edge between two boundary vertices"
  let keep := min a b
  let gone := max a b
  let newPos :=
    if pa then m.pos[a]! else if pb then m.pos[b]!
    else if t.isBoundary[a]! && !t.isBoundary[b]! then m.pos[a]!
    else if t.isBoundary[b]! && !t.isBoundary[a]! then m.pos[b]!
    else (0.5 : Float) * (m.pos[a]! + m.pos[b]!)
  let relabel := fun (v : Nat) => if v == gone then keep else if v > gone then v - 1 else v
  let mut faces : Array Tri := #[]
  let mut parent : Array Nat := #[]
  for fi in [0:m.nf] do
    let f := m.faces[fi]!
    let hasA := f.a == a || f.b == a || f.c == a
    let hasB := f.a == b || f.b == b || f.c == b
    if !(hasA && hasB) then
      if foldover && (hasA || hasB) then
        let moved := m.pos.set! keep newPos |>.set! gone newPos
        let g0 := faceGeo #[f] m.pos
        let g1 := faceGeo #[f] moved
        if !(g1[0]!.area > 0.0) then throw s!"Foldover: face {fi} becomes degenerate"
        if angleBetween g0[0]!.n g1[0]!.n > foldoverLimit then
          throw s!"Foldover: face {fi} turns its normal by more than pi - 2"
      faces := faces.push ⟨relabel f.a, relabel f.b, relabel f.c⟩
      parent := parent.push fi
  let key := fun (f : Tri) => (#[f.a, f.b, f.c].qsort (· < ·))
  for i in [0:faces.size] do
    for j in [i+1:faces.size] do
      if key faces[i]! == key faces[j]! then
        throw s!"NonManifoldResult: faces {parent[i]!} and {parent[j]!} would have the same three vertices"
  let mut pos : Array V3 := #[]
  let mut vmap : Array (Array Nat) := #[]
  for v in [0:m.nv] do
    if v != gone then
      pos := pos.push (if v == keep then newPos else m.pos[v]!)
      vmap := vmap.push (if v == keep then #[a, b] else #[v])
  match buildTopo pos.size faces with
  | .error e => throw s!"NonManifoldResult: {e}"
  | .ok _ => pure ()
  return { pos, faces, parent, vmap }

/-- Edge flip (interior edge, consistent orientation, no existing edge `c d`; Delaunay and foldover
checks only if asked). -/
def flipEdge (m : Mesh) (t : Topo) (a b : Nat) (delaunay foldover : Bool) : Except String RemeshOut := do
  let inc ← findEdge t a b
  if inc.size != 2 then throw s!"BoundaryEdge: ({a}, {b}) has {inc.size} face(s); only interior edges flip"
  let sorted := inc.qsort fun x y => x.1 < y.1
  let f1 := sorted[0]!.1
  let f2 := sorted[1]!.1
  if sorted[0]!.2 == sorted[1]!.2 then
    throw s!"InconsistentOrientation: faces {f1} and {f2} traverse ({a}, {b}) in the same direction"
  let some (p, q, c) := orient m.faces[f1]! a b | throw "internal"
  let some (_, _, d) := orient m.faces[f2]! a b | throw "internal"
  if (t.edges.find? fun e => e.1 == min c d && e.2.1 == max c d).isSome then
    throw s!"FlipNonManifold: vertices {c} and {d} already share an edge"
  let ang := fun (v w x : Nat) => angleBetween (m.pos[w]! - m.pos[v]!) (m.pos[x]! - m.pos[v]!)
  if delaunay && !(ang c p q + ang d p q > 3.141592653589793) then
    throw s!"NotDelaunay: opposite angles {ang c p q} + {ang d p q} do not exceed pi"
  let n1 : Tri := ⟨p, d, c⟩
  let n2 : Tri := ⟨d, q, c⟩
  if foldover then
    let g := faceGeo #[m.faces[f1]!, m.faces[f2]!, n1, n2] m.pos
    if !(g[2]!.area > 0.0) || !(g[3]!.area > 0.0) then throw "Foldover: a new face has zero area"
    if angleBetween g[0]!.n g[2]!.n > foldoverLimit || angleBetween g[1]!.n g[3]!.n > foldoverLimit then
      throw "Foldover: a new face turns its normal by more than pi - 2"
  let faces := (m.faces.set! f1 n1).set! f2 n2
  match buildTopo m.nv faces with
  | .error e => throw s!"NonManifoldResult: {e}"
  | .ok _ => pure ()
  return { pos := m.pos, faces, parent := Array.range m.nf, vmap := (Array.range m.nv).map (fun v => #[v]) }

end Ddg.RepulsiveSurfaces
