import DdgOracle.Formulas.Mesh

/-!
# Mesh connectivity and validation

Input meshes must be connected, manifold, consistently oriented triangle meshes, with or without
boundary. Everything else is rejected with a message that says what is wrong and where.
-/

namespace Ddg

/-- Largest number of vertices the dense oracle accepts. -/
def maxVertices : Nat := 2000

/-- Validated connectivity of a triangle mesh. -/
structure Topology where
  /-- number of vertices -/
  nv : Nat
  /-- the triangles -/
  faces : Array Face
  /-- `isBoundary[v]` iff `v` lies on a boundary edge -/
  isBoundary : Array Bool
  /-- indices of the faces at each vertex, increasing -/
  vertexFaces : Array (Array Nat)
  /-- vertices adjacent to each vertex, increasing -/
  neighbors : Array (Array Nat)
  /-- number of edges -/
  numEdges : Nat
  /-- number of boundary vertices -/
  numBoundary : Nat
  /-- Euler characteristic `V − E + F` -/
  chi : Int

namespace Topology

/-- The faces at `v`, in input order. -/
def facesAt (t : Topology) (v : Nat) : List Face :=
  (t.vertexFaces[v]!.map fun f => t.faces[f]!).toList

/-- `true` iff the mesh has no boundary. -/
def closed (t : Topology) : Bool := t.numBoundary == 0

end Topology

private def insertSorted (xs : Array Nat) (x : Nat) : Array Nat :=
  if xs.contains x then xs else (xs.push x).qsort (· < ·)

/-- Validate connectivity and build the adjacency tables. `connected := false` accepts meshes
with several connected components (the cut flow). -/
def buildTopology (nv : Nat) (faces : Array Face) (connected : Bool := true) :
    Except String Topology := do
  if nv < 3 then throw s!"mesh has {nv} vertices; at least 3 are required"
  if nv > maxVertices then
    throw s!"mesh has {nv} vertices; the dense oracle accepts at most {maxVertices}"
  if faces.size == 0 then throw "mesh has no faces"
  -- directed edges: `he[a * nv + b] = f + 1` if face `f` contains the directed edge `a → b`
  let mut he : Array Nat := Array.replicate (nv * nv) 0
  let mut vertexFaces : Array (Array Nat) := Array.replicate nv #[]
  let mut neighbors : Array (Array Nat) := Array.replicate nv #[]
  for fi in [0:faces.size] do
    let f := faces[fi]!
    if f.a ≥ nv || f.b ≥ nv || f.c ≥ nv then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) has a vertex index outside 0..{nv - 1}"
    if f.a == f.b || f.b == f.c || f.c == f.a then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) has a repeated vertex"
    for (p, q) in [(f.a, f.b), (f.b, f.c), (f.c, f.a)] do
      let g := he[p * nv + q]!
      if g != 0 then
        throw s!"directed edge ({p}, {q}) occurs in faces {g - 1} and {fi}: the mesh is non-manifold or not consistently oriented"
      he := he.set! (p * nv + q) (fi + 1)
      neighbors := neighbors.set! p (insertSorted neighbors[p]! q)
      neighbors := neighbors.set! q (insertSorted neighbors[q]! p)
    for v in [f.a, f.b, f.c] do
      vertexFaces := vertexFaces.set! v (vertexFaces[v]!.push fi)
  -- boundary edges and vertices
  let mut isBoundary : Array Bool := Array.replicate nv false
  let mut outB : Array Nat := Array.replicate nv 0
  let mut inB : Array Nat := Array.replicate nv 0
  let mut directed := 0
  let mut boundaryEdges := 0
  for p in [0:nv] do
    for q in neighbors[p]! do
      if he[p * nv + q]! != 0 then
        directed := directed + 1
        if he[q * nv + p]! == 0 then
          boundaryEdges := boundaryEdges + 1
          isBoundary := (isBoundary.set! p true).set! q true
          outB := outB.set! p (outB[p]! + 1)
          inB := inB.set! q (inB[q]! + 1)
  for v in [0:nv] do
    if vertexFaces[v]!.size == 0 then throw s!"vertex {v} is not used by any face"
    if outB[v]! > 1 || inB[v]! > 1 then
      throw s!"vertex {v} lies on {outB[v]! + inB[v]!} boundary edges: the mesh is non-manifold at this vertex"
  -- every vertex star must be a single fan (boundary) or a single cycle (interior)
  for v in [0:nv] do
    let star := vertexFaces[v]!
    -- next face around `v`: across the edge `q → v` of the face `(v, p, q)`
    let step := fun (fi : Nat) =>
      let f := faces[fi]!
      let q := if f.a == v then f.c else if f.b == v then f.a else f.b
      he[v * nv + q]!
    -- first face: the one whose edge leaving `v` is a boundary edge, or any face
    let mut start := star[0]!
    for fi in star do
      let f := faces[fi]!
      let p := if f.a == v then f.b else if f.b == v then f.c else f.a
      if he[p * nv + v]! == 0 then start := fi
    let mut count := 1
    let mut cur := start
    for _ in [0:star.size] do
      let nxt := step cur
      if nxt == 0 || nxt - 1 == start then break
      cur := nxt - 1
      count := count + 1
    if count != star.size then
      throw s!"vertex {v} is non-manifold: its {star.size} faces do not form a single fan"
  -- connectedness
  let mut seen : Array Bool := Array.replicate nv false
  let mut stack : Array Nat := #[0]
  seen := seen.set! 0 true
  let mut reached := 1
  while stack.size > 0 do
    let v := stack.back!
    stack := stack.pop
    for w in neighbors[v]! do
      if !seen[w]! then
        seen := seen.set! w true
        reached := reached + 1
        stack := stack.push w
  if connected && reached != nv then
    throw s!"mesh is not connected: only {reached} of {nv} vertices are reachable from vertex 0"
  let numEdges := (directed + boundaryEdges) / 2
  let numBoundary := (isBoundary.filter id).size
  return {
    nv, faces, isBoundary, vertexFaces, neighbors, numEdges, numBoundary
    chi := (nv : Int) - (numEdges : Int) + (faces.size : Int)
  }

end Ddg
