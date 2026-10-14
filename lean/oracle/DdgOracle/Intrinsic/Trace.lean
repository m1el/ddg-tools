import DdgOracle.Intrinsic.Flip

/-!
# Tracing (SPEC §4, N7)

The triangles along the path are unfolded one after the other into one plane in which the path is the
positive x-axis; the edge through which the path leaves a triangle is decided by the signs of the
y-coordinates (`it-unfold-apex`), the crossing is `it-ray-edge-intersection`.
-/

namespace Ddg.Intrinsic

/-- Result of a tracing query. -/
structure TraceResult where
  /-- `end` (length used up inside a face), `edge` (along an edge, inside it), `vertex` (a vertex was
  reached with length left, or the end of an edge), `boundary` (the boundary was reached) -/
  stop : String
  /-- the face of the end point -/
  face : Nat
  /-- barycentric coordinates in the corner order of the face (tails of `fh`, `next`, `next next`) -/
  bary : Float × Float × Float
  /-- length left unused -/
  left : Float
  /-- the vertex reached (`noTwin` if none) -/
  vertex : Nat
  /-- crossed halfedges with the parameter from their tail -/
  crossings : Array (Nat × Float)
  /-- the direction was moved into the corner of the vertex (N3) -/
  clamp : Bool
  deriving Inhabited

/-- The three halfedges of the face of `h0` starting at `h0`. -/
def Tri.hs3 (t : Tri) (h0 : Nat) : Array Nat := #[h0, t.next[h0]!, t.next[t.next[h0]!]!]

/-- The face laid out in the plane: corner 0 at the origin, corner 1 on the x-axis. -/
def layoutFace (t : Tri) (h0 : Nat) : Array (Float × Float) :=
  let hs := t.hs3 h0
  let l0 := t.len[hs[0]!]!
  let apex := unfoldApex 0.0 0.0 l0 0.0 l0 t.len[hs[1]!]! t.len[hs[2]!]!
  #[(0.0, 0.0), (l0, 0.0), apex]

/-- Barycentric coordinates in the canonical corner order of the face from those in the order of `hs`. -/
def canonBary (t : Tri) (f : Nat) (hs : Array Nat) (b : Float × Float × Float) : Float × Float × Float :=
  let arr := #[b.1, b.2.1, b.2.2]
  let r := if hs[0]! == t.fh[f]! then 0 else if hs[1]! == t.fh[f]! then 1 else 2
  (arr[r % 3]!, arr[(r + 1) % 3]!, arr[(r + 2) % 3]!)

/-- Core of tracing: the face of `h0` laid out as `layoutFace`; the ray starts at `start` in that layout
with the unit direction `dir`, length `dist`. -/
def traceCore (t : Tri) (h0 : Nat) (start dir : Float × Float) (dist : Float) (maxCross : Nat) :
    Except String TraceResult := do
  let P := layoutFace t h0
  let (dx, dy) := dir
  let toRay := fun (p : Float × Float) =>
    let x := p.1 - start.1
    let y := p.2 - start.2
    (dx * x + dy * y, -dy * x + dx * y)
  let mut hs := t.hs3 h0
  let mut Q : Array (Float × Float) := P.map toRay
  let mut xcur := fzero
  let mut entering : Nat := 3      -- no entering edge in the first triangle
  let mut crossings : Array (Nat × Float) := #[]
  let mut result : Option TraceResult := none
  let mut count := 0
  while result.isNone do
    -- the exit edge
    let mut exitK := 3
    let mut tt := fzero
    let mut xc := fzero
    for k in [0:3] do
      if exitK == 3 && k != entering then
        let a := Q[k]!
        let b := Q[(k + 1) % 3]!
        if (a.2 < 0 && b.2 > 0) || (a.2 > 0 && b.2 < 0) then
          let t' := axisCrossingT a.2 b.2
          let x := axisCrossingX a.1 b.1 t'
          if entering != 3 || x > xcur then
            exitK := k
            tt := t'
            xc := x
    let f := t.face[hs[0]!]!
    let bary := fun (x : Float) =>
      canonBary t f hs (barycentric Q[0]!.1 Q[0]!.2 Q[1]!.1 Q[1]!.2 Q[2]!.1 Q[2]!.2 x 0.0)
    if exitK == 3 then
      -- through a vertex on the axis?
      let mut vm := 3
      for m in [0:3] do
        if vm == 3 && Q[m]!.2 == 0 && Q[m]!.1 > xcur then vm := m
      if vm == 3 then
        throw s!"the ray leaves face {f} through no edge and no vertex"
      let xv := Q[vm]!.1
      if dist ≤ xv then
        result := some { stop := "end", face := f, bary := bary dist, left := 0.0, vertex := noTwin, crossings, clamp := false }
      else
        let bv := if vm == 0 then (1.0, 0.0, 0.0) else if vm == 1 then (0.0, 1.0, 0.0) else (0.0, 0.0, 1.0)
        result := some { stop := "vertex", face := f, bary := canonBary t f hs bv, left := dist - xv
                         vertex := t.tail[hs[vm]!]!, crossings, clamp := false }
    else if dist ≤ xc then
      result := some { stop := "end", face := f, bary := bary dist, left := 0.0, vertex := noTwin, crossings, clamp := false }
    else
      let hk := hs[exitK]!
      crossings := crossings.push (hk, tt)
      let tw := t.twin[hk]!
      if tw == noTwin then
        result := some { stop := "boundary", face := f, bary := bary xc, left := dist - xc, vertex := noTwin
                         crossings, clamp := false }
      else
        count := count + 1
        if count > maxCross then throw s!"more than {maxCross} crossings: the trace does not end"
        let A := Q[exitK]!
        let B := Q[(exitK + 1) % 3]!
        let n1 := t.next[tw]!
        let n2 := t.next[n1]!
        let C := unfoldApex B.1 B.2 (A.1 - B.1) (A.2 - B.2) t.len[tw]! t.len[n1]! t.len[n2]!
        hs := #[tw, n1, n2]
        Q := #[B, A, C]
        xcur := xc
        entering := 0
  match result with
  | some r => return r
  | none => throw "trace: no result"

/-- Trace from a vertex in a normalized direction (SPEC Alg. 2). -/
def traceFromVertex (t : Tri) (v : Nat) (phi dist : Float) (maxCross : Nat := 100000) :
    Except String TraceResult := do
  if v ≥ t.nv then throw s!"vertex {v} does not exist"
  if !(dist ≥ 0 && dist.isFinite) then throw "the distance must be finite and non-negative"
  let full := if t.isB[v]! then floatPi else 2 * floatPi
  if !(phi ≥ 0 && phi ≤ full && phi.isFinite) then
    throw s!"the direction {phi} is outside [0, {full}] at vertex {v}"
  -- the wedge: the outgoing halfedge with the largest signpost not above the direction
  let mut best := t.nh
  for h in [0:t.nh] do
    if t.tail[h]! == v && t.sign[h]! ≤ phi then
      if best == t.nh || t.sign[h]! > t.sign[best]! then best := h
  if best == t.nh then throw s!"vertex {v} has no outgoing halfedge"
  let delta := phi - t.sign[best]!
  let f := t.face[best]!
  if delta == 0 then
    -- along the edge
    let l := t.len[best]!
    let bv := fun (s : Float) => canonBary t f (t.hs3 best) (1 - s, s, 0.0)
    if dist < l then
      return { stop := "edge", face := f, bary := bv (dist / l), left := 0.0, vertex := noTwin, crossings := #[], clamp := false }
    else
      return { stop := "vertex", face := f, bary := bv 1.0, left := dist - l, vertex := t.head best, crossings := #[], clamp := false }
  let corner := t.corner best
  let mut phit := delta / t.scaleAt v
  let mut clamp := false
  if phit ≥ corner then
    phit := corner * (1 - 2.220446049250313e-16)
    clamp := true
  let r ← traceCore t best (0.0, 0.0) (Float.cos phit, Float.sin phit) dist maxCross
  return { r with clamp }

/-- Trace from a point of a face (layout of `layoutFace t (fh f)`), direction a unit vector. -/
def traceFromFace (t : Tri) (f : Nat) (start dir : Float × Float) (dist : Float) (maxCross : Nat := 100000) :
    Except String TraceResult :=
  traceCore t t.fh[f]! start dir dist maxCross

/-- Position in the layout of a point with barycentric coordinates `b` in the face `f`. -/
def facePoint (t : Tri) (f : Nat) (b : Float × Float × Float) : Float × Float :=
  let P := layoutFace t t.fh[f]!
  (b.1 * P[0]!.1 + b.2.1 * P[1]!.1 + b.2.2 * P[2]!.1, b.1 * P[0]!.2 + b.2.1 * P[1]!.2 + b.2.2 * P[2]!.2)

/-- The point of a trace result in the layout of its face. -/
def endPoint (t : Tri) (r : TraceResult) : Float × Float := facePoint t r.face r.bary

/-- Distance of the end of a trace from the vertex `w`, if `w` is a corner of its face (`posInf`
otherwise): the check that a trace of an edge arrives at its far vertex (SPEC N8). -/
def missDistance (t : Tri) (r : TraceResult) (w : Nat) : Float := Id.run do
  let hs := t.faceHalfedges r.face
  let P := layoutFace t t.fh[r.face]!
  let e := endPoint t r
  let mut best := posInf
  for k in [0:3] do
    if t.tail[hs[k]!]! == w then
      let d := planeDistance e.1 e.2 P[k]!.1 P[k]!.2
      if d < best then best := d
  return best

end Ddg.Intrinsic
