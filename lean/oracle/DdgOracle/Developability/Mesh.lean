import DdgOracle.Developability.Vec
import DdgOracle.Formulas.Mesh

/-!
# Meshes for the developability oracle: validation, stars, per-face data

A mesh is positions and CCW triangles. Validation (manifold, consistently oriented, no isolated
vertex) and the vertex stars of SPEC section 1 are computed from scratch with a dense table of
directed edges (the oracle is meant for meshes of up to a few hundred vertices).
-/

namespace Ddg.Developability

/-- Largest number of vertices accepted. -/
def maxVertices : Nat := 600

/-- Positions and triangles. -/
structure DMesh where
  /-- positions -/
  pos : Array V3
  /-- triangles, counter-clockwise -/
  faces : Array (Nat × Nat × Nat)
  deriving Inhabited

/-- The vertex of slot `s ∈ {0,1,2}` of a face. -/
def faceVertex (f : Nat × Nat × Nat) (s : Nat) : Nat :=
  if s % 3 == 0 then f.1 else if s % 3 == 1 then f.2.1 else f.2.2

/-- The star of a vertex: its faces in counter-clockwise order, with the slot of the vertex in
each. Interior: a cycle starting at the face of smallest index. Boundary: a path starting at the
boundary halfedge leaving the vertex. -/
structure Star where
  /-- the vertex -/
  v : Nat
  /-- faces in order -/
  faces : Array Nat
  /-- slot of `v` in each face -/
  slots : Array Nat
  /-- `true` iff interior -/
  closed : Bool
  deriving Inhabited

/-- A validated mesh with stars. -/
structure Topo where
  /-- number of vertices -/
  nv : Nat
  /-- the stars -/
  stars : Array Star
  /-- edges `(a, b)`, `a < b`, sorted lexicographically -/
  edges : Array (Nat × Nat)

/-- Validate the connectivity and build the stars. -/
def buildTopo (m : DMesh) : Except String Topo := do
  let nv := m.pos.size
  if nv < 3 then throw s!"mesh has {nv} vertices; at least 3 are required"
  if nv > maxVertices then throw s!"mesh has {nv} vertices; the oracle accepts at most {maxVertices}"
  if m.faces.size == 0 then throw "mesh has no faces"
  let mut he : Array Nat := Array.replicate (nv * nv) 0
  let mut incident : Array (Array (Nat × Nat)) := Array.replicate nv #[]
  for fi in [0:m.faces.size] do
    let f := m.faces[fi]!
    for s in [0:3] do
      let a := faceVertex f s
      let b := faceVertex f (s + 1)
      if a ≥ nv || b ≥ nv then throw s!"face {fi} has a vertex index out of range (the mesh has {nv} vertices)"
      if a == b then throw s!"face {fi} is degenerate: repeated vertex {a}"
      if he[a * nv + b]! != 0 then throw s!"directed edge ({a}, {b}) is used twice (non-manifold or inconsistent orientation)"
      he := he.set! (a * nv + b) (fi + 1)
      incident := incident.modify a (·.push (fi, s))
  let mut stars : Array Star := #[]
  for v in [0:nv] do
    let inc := incident[v]!
    if inc.size == 0 then throw s!"vertex {v} is isolated"
    -- boundary start: a face whose outgoing edge v -> a has no opposite a -> v
    let mut start := inc[0]!
    let mut starts := 0
    for (fi, s) in inc do
      let a := faceVertex m.faces[fi]! (s + 1)
      if he[a * nv + v]! == 0 then
        start := (fi, s)
        starts := starts + 1
    if starts > 1 then throw s!"vertex {v} is a non-manifold vertex (several boundary fans)"
    let closed := starts == 0
    let mut faces : Array Nat := #[]
    let mut slots : Array Nat := #[]
    let mut cur := start
    let mut fuel := inc.size + 1
    let mut done := false
    while !done && fuel > 0 do
      fuel := fuel - 1
      faces := faces.push cur.1
      slots := slots.push cur.2
      let b := faceVertex m.faces[cur.1]! (cur.2 + 2)
      let nxt := he[v * nv + b]!
      if nxt == 0 then done := true
      else
        let nf := nxt - 1
        if nf == start.1 && closed then done := true
        else
          let ns := (List.range 3).find? (fun s => faceVertex m.faces[nf]! s == v)
          cur := (nf, ns.getD 0)
    if faces.size != inc.size then throw s!"vertex {v} is a non-manifold vertex (its faces do not form one fan)"
    stars := stars.push { v, faces, slots, closed }
  let mut edges : Array (Nat × Nat) := #[]
  for a in [0:nv] do
    for b in [a + 1:nv] do
      if he[a * nv + b]! != 0 || he[b * nv + a]! != 0 then edges := edges.push (a, b)
  return { nv, stars, edges }

/-- Per-face geometry. -/
structure FaceData where
  /-- `n = (f_j - f_i) × (f_k - f_i)` -/
  n : V3
  /-- unit normal -/
  N : V3
  /-- `|n|` -/
  nn : Float
  /-- corner angles by slot -/
  theta : Array Float
  deriving Inhabited

/-- Corner angle at `p` between the edges to `q` and `r`: `atan2 (|a × b|, a · b)`. -/
def cornerAngle (p q r : V3) : Float :=
  let a := q - p
  let b := r - p
  Float.atan2 (V3.cross a b).norm (V3.dot a b)

/-- Geometry of face `fi`; `none` when the face is degenerate (zero or non-finite normal). -/
def faceData (m : DMesh) (fi : Nat) : Option FaceData :=
  let f := m.faces[fi]!
  let pi := m.pos[f.1]!
  let pj := m.pos[f.2.1]!
  let pk := m.pos[f.2.2]!
  let n := V3.cross (pj - pi) (pk - pi)
  let nn := n.norm
  if !(nn.isFinite) || nn == 0 then none
  else some { n, N := V3.smul (1 / nn) n, nn,
              theta := #[cornerAngle pi pj pk, cornerAngle pj pk pi, cornerAngle pk pi pj] }

end Ddg.Developability
