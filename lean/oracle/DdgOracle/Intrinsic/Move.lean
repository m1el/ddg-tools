import DdgOracle.Intrinsic.Trace

/-!
# Vertex insertion, edge split and moving an inserted vertex (SPEC §4, N10, N19, N21, N24)

Written from the paper (Alg. 7, 8, 10, Sec. 3.2.3, 3.3.3) and `crates/ddg-intrinsic/SPEC.md`; the
crate's code was not read. Insertion and split act on the triangulation fresh from the input mesh
(the task `move` inserts one vertex first), so the face and the edge are input elements and the place
of the new vertex is known exactly. Moving uses `movedDoubleArea` (`it-move-area`) for the star test,
`movedLength` (`it-move-length`) for the new lengths, a trace over the input mesh for the new place,
and the directions of the neighbours from the new place in the unfolded star for the new signposts.
-/

namespace Ddg.Intrinsic

open Ddg

/-- The place of a vertex on the input mesh (SPEC §3). -/
inductive Loc where
  /-- an input vertex -/
  | vertex (v : Nat)
  /-- a point of the input edge `{lo, hi}`, `lo < hi`, at the parameter `t` from `lo` -/
  | edge (lo hi : Nat) (t : Float)
  /-- a point of the input face `f`, barycentric coordinates in the corner order of its `f` line -/
  | face (f : Nat) (b : Float × Float × Float)
  deriving Inhabited

/-- `vertex V`, `edge LO HI T`, `face F B0 B1 B2`. -/
def Loc.print : Loc → String
  | .vertex v => s!"vertex {v}"
  | .edge lo hi t => s!"edge {lo} {hi} {floatToHex t}"
  | .face f (b0, b1, b2) => s!"face {f} {floatToHex b0} {floatToHex b1} {floatToHex b2}"

/-- `atan2` into `[0, 2π)`. -/
def angleOf (x y : Float) : Float :=
  let a := Float.atan2 y x
  if a < 0 then a + 2 * floatPi else a

/-- An angle into `[0, 2π)` (a few turns at most). -/
def wrapTurn (a : Float) : Float := Id.run do
  let mut x := a
  for _ in [0:4] do
    if x < 0 then x := x + 2 * floatPi
    if x ≥ 2 * floatPi then x := x - 2 * floatPi
  return x

/-- Rotate a vector by the angle `a`. -/
def rotateBy (a : Float) (p : Float × Float) : Float × Float :=
  (p.1 * Float.cos a - p.2 * Float.sin a, p.1 * Float.sin a + p.2 * Float.cos a)

/-- The triangulation with one more vertex (no halfedges yet); `Θ = π` on the boundary, `2π` inside. -/
def Tri.addVertex (t : Tri) (onBoundary : Bool) : Tri :=
  { t with nv := t.nv + 1, theta := t.theta.push (if onBoundary then floatPi else 2 * floatPi)
           isB := t.isB.push onBoundary }

/-- Room for `k` more halfedges (`k / 3` more faces); the new entries are placeholders. -/
def Tri.grow (t : Tri) (k : Nat) : Tri :=
  { t with nh := t.nh + k
           tail := t.tail ++ Array.replicate k 0, next := t.next ++ Array.replicate k 0
           twin := t.twin ++ Array.replicate k noTwin, face := t.face ++ Array.replicate k 0
           fh := t.fh ++ Array.replicate (k / 3) 0, len := t.len ++ Array.replicate k 0.0
           sign := t.sign ++ Array.replicate k 0.0, eid := t.eid ++ Array.replicate k 0 }

/-- Set the three halfedges of face `f` (in order) and their `next`, `face`, `fh`. -/
def Tri.setFace (t : Tri) (f a b c : Nat) : Tri :=
  { t with next := ((t.next.set! a b).set! b c).set! c a
           face := ((t.face.set! a f).set! b f).set! c f
           fh := t.fh.set! f a }

/-- Make `a`, `b` twins with the edge id `min a b` and the length `l`. -/
def Tri.pair (t : Tri) (a b : Nat) (l : Float) : Tri :=
  { t with twin := (t.twin.set! a b).set! b a, eid := (t.eid.set! a (min a b)).set! b (min a b)
           len := (t.len.set! a l).set! b l }

/-- A halfedge without twin with the length `l`. -/
def Tri.single (t : Tri) (a : Nat) (l : Float) : Tri :=
  { t with twin := t.twin.set! a noTwin, eid := t.eid.set! a a, len := t.len.set! a l }

/-- Insert a vertex into face `f` of a triangulation that is fresh from the input mesh (Alg. 7), at the
barycentric coordinates `b` in the corner order of the face. Returns the triangulation and the place.
The new vertex has the id `t.nv`; its signposts are angles in the layout of `f` (N10). -/
def insertInFace (t : Tri) (f : Nat) (b : Float × Float × Float) : Except String (Tri × Loc) := do
  if f ≥ t.nf then throw s!"face {f} does not exist"
  let (b0, b1, b2) := b
  if !(b0 > 0 && b1 > 0 && b2 > 0) then throw "insert: the barycentric coordinates must be positive"
  if (b0 + b1 + b2 - 1).abs > 1e-12 then throw "insert: the barycentric coordinates must sum to 1"
  let hs := t.faceHalfedges f
  let h0 := hs[0]!
  let h1 := hs[1]!
  let h2 := hs[2]!
  let i := t.tail[h0]!
  let j := t.tail[h1]!
  let k := t.tail[h2]!
  let (di, dj, dk) := cornerDistances t.len[h0]! t.len[h1]! t.len[h2]! b0 b1 b2
  let P := layoutFace t h0
  let px := b0 * P[0]!.1 + b1 * P[1]!.1 + b2 * P[2]!.1
  let py := b0 * P[0]!.2 + b1 * P[1]!.2 + b2 * P[2]!.2
  let p := t.nv
  let nh := t.nh
  let g := t.nf
  let a := nh       -- j → p
  let bb := nh + 1  -- p → i
  let c := nh + 2   -- k → p
  let d := nh + 3   -- p → j
  let e := nh + 4   -- i → p
  let ff := nh + 5  -- p → k
  let mut u := (t.addVertex false).grow 6
  u := { u with tail := (((((u.tail.set! a j).set! bb p).set! c k).set! d p).set! e i).set! ff p }
  u := u.setFace f h0 a bb
  u := u.setFace g h1 c d
  u := u.setFace (g + 1) h2 e ff
  u := ((u.pair a d dj).pair bb e di).pair c ff dk
  -- outgoing signposts of p: directions in the layout of f
  let dirTo := fun (q : Float × Float) => angleOf (q.1 - px) (q.2 - py)
  let mut sg := u.sign
  sg := ((sg.set! bb (dirTo P[0]!)).set! d (dirTo P[1]!)).set! ff (dirTo P[2]!)
  -- incoming: the next halfedge counter-clockwise at the corner, by the new corner angle
  sg := sg.set! e (signpostUpdate t.isB[i]! t.sign[h0]! (u.corner h0) (t.theta.get! i))
  sg := sg.set! a (signpostUpdate t.isB[j]! t.sign[h1]! (u.corner h1) (t.theta.get! j))
  sg := sg.set! c (signpostUpdate t.isB[k]! t.sign[h2]! (u.corner h2) (t.theta.get! k))
  u := { u with sign := sg }
  return (u, .face f b)

/-- Split the edge `{A, B}` of a triangulation fresh from the input mesh at the parameter `T` from `A`
(Alg. 8, SPEC N10, N17): the half from the tail keeps the halfedge, lengths `T ℓ` and `(1 − T) ℓ`
(`edgeSplitLengths`). On a boundary edge the new vertex is a boundary vertex (signpost 0 along the
boundary halfedge that leaves it); inside, the angle 0 is along the edge from its lower-numbered end. -/
def splitEdge (t : Tri) (A B : Nat) (T : Float) : Except String (Tri × Loc) := do
  if !(T > 0 && T < 1) then throw "split: the parameter must lie in (0, 1)"
  -- the halfedge A → B, or B → A (then the parameter from its tail is 1 − T)
  let (h, tt) ← match findHalfedge t A B with
    | .ok h => pure (h, T)
    | .error _ => do
      let h ← findHalfedge t B A
      pure (h, 1 - T)
  let tw := t.twin[h]!
  let onB := tw == noTwin
  let h1 := t.next[h]!
  let h2 := t.next[h1]!
  let a := t.tail[h]!
  let b := t.tail[h1]!
  let c := t.tail[h2]!
  let l := t.len[h]!
  let (lap, lpb) := edgeSplitLengths l (1 - tt) tt
  let (_, _, lpc) := cornerDistances l t.len[h1]! t.len[h2]! (1 - tt) tt 0.0
  let p := t.nv
  let nh := t.nh
  let g := t.nf
  let f1 := t.face[h]!
  let x1 := nh      -- p → c
  let y1 := nh + 1  -- p → b
  let z1 := nh + 2  -- c → p
  let mut u := (t.addVertex onB).grow (if onB then 3 else 6)
  u := { u with tail := ((u.tail.set! x1 p).set! y1 p).set! z1 c }
  u := u.setFace f1 h x1 h2
  u := u.setFace g y1 h1 z1
  u := u.pair x1 z1 lpc
  let lo := min a b
  let hi := max a b
  let loc : Loc := .edge lo hi (if a < b then tt else 1 - tt)
  if onB then
    u := (u.single h lap).single y1 lpb
    let mut sg := u.sign
    sg := (sg.set! y1 0.0).set! x1 (u.corner y1)
    sg := sg.set! z1 (signpostUpdate t.isB[c]! t.sign[h2]! (u.corner h2) (t.theta.get! c))
    u := { u with sign := sg }
    return (u, loc)
  else
    let t1 := t.next[tw]!
    let t2 := t.next[t1]!
    let dd := t.tail[t2]!
    let f2 := t.face[tw]!
    let (_, _, lpd) := cornerDistances l t.len[t1]! t.len[t2]! tt (1 - tt) 0.0
    let x2 := nh + 3  -- p → d
    let y2 := nh + 4  -- p → a
    let z2 := nh + 5  -- d → p
    u := { u with tail := ((u.tail.set! x2 p).set! y2 p).set! z2 dd }
    u := u.setFace f2 tw x2 t2
    u := u.setFace (g + 1) y2 t1 z2
    u := ((u.pair h y2 lap).pair tw y1 lpb).pair x2 z2 lpd
    let mut sg := u.sign
    if a < b then
      sg := (((sg.set! y1 0.0).set! x1 (u.corner y1)).set! y2 floatPi).set! x2 (floatPi + u.corner y2)
    else
      sg := (((sg.set! y2 0.0).set! x2 (u.corner y2)).set! y1 floatPi).set! x1 (floatPi + u.corner y1)
    sg := sg.set! z1 (signpostUpdate t.isB[c]! t.sign[h2]! (u.corner h2) (t.theta.get! c))
    sg := sg.set! z2 (signpostUpdate t.isB[dd]! t.sign[t2]! (u.corner t2) (t.theta.get! dd))
    u := { u with sign := sg }
    return (u, loc)

/-- Do two triangulations have the same halfedges `(tail, head, length, signpost)` (halfedge ids, face
ids and rotations are ignored)? Returns the largest deviation, or `posInf` if they do not match. -/
def sameTriangulation (a b : Tri) : Float := Id.run do
  if a.nh != b.nh then return posInf
  let mut used : Array Bool := Array.replicate b.nh false
  let mut worst := fzero
  for h in [0:a.nh] do
    let mut best := posInf
    let mut bi := b.nh
    for g in [0:b.nh] do
      if !used[g]! && a.tail[h]! == b.tail[g]! && a.head h == b.head g then
        let ds := (a.sign[h]! - b.sign[g]!).abs
        let d := max (a.len[h]! - b.len[g]!).abs (min ds (2 * floatPi - ds))
        if d < best then
          best := d
          bi := g
    if bi == b.nh then return posInf
    used := used.set! bi true
    worst := max worst best
  return worst

/-- Outgoing halfedges of `v`. -/
def Tri.outgoing (t : Tri) (v : Nat) : Array Nat := Id.run do
  let mut out := #[]
  for h in [0:t.nh] do
    if t.tail[h]! == v then out := out.push h
  return out

/-- The triangle of the star of a flat vertex `v` at its outgoing halfedge `h`, unfolded with `v` at the
origin: the far end `j` of `h` and the third corner `k`, and twice the area. -/
def starTriangle (t : Tri) (h : Nat) : (Float × Float) × (Float × Float) × Float :=
  let phi := t.sign[h]!
  let th := t.corner h
  let lj := t.len[h]!
  let lk := t.len[t.prev h]!
  ((lj * Float.cos phi, lj * Float.sin phi), (lk * Float.cos (phi + th), lk * Float.sin (phi + th)),
   2 * triArea lj t.len[t.next[h]!]! lk)

/-- Largest admissible fraction of the displacement `d` (SPEC N19): `min 2A/(d × (k − j))` over the
triangles of the star with `d × (k − j) > 0`; `posInf` if there is none. Also whether every triangle
keeps a positive area at the full displacement (`movedDoubleArea`). -/
def starKernel (t : Tri) (v : Nat) (d : Float × Float) : Float × Bool := Id.run do
  let mut frac := posInf
  let mut inside := true
  for h in t.outgoing v do
    let (j, k, a2) := starTriangle t h
    let side := (k.1 - j.1, k.2 - j.2)
    let cross := d.1 * side.2 - d.2 * side.1
    if cross > 0 then frac := min frac (a2 / cross)
    if !(movedDoubleArea a2 d side > 0) then inside := false
  return (frac, inside)

/-- The input halfedge `a → b`, if any. -/
def inputHalfedge (t0 : Tri) (a b : Nat) : Option Nat :=
  match findHalfedge t0 a b with
  | .ok h => some h
  | .error _ => none

/-- The snapping threshold of a moved vertex (SPEC N24): relative to the longest edge of the input face
and to the shortest new edge at the vertex. -/
def moveSnap : Float := 1e-12

/-- The shape limit of N5, N19. -/
def moveShapeLimit : Float := 1e-6

/-- Result of a move. -/
inductive MoveOutcome where
  /-- moved: the new triangulation, the new place, whether it was snapped onto an input edge, and the
  angle from the old frame of the vertex to the new one -/
  | ok (t : Tri) (loc : Loc) (snapped : Bool) (rho : Float)
  /-- refused, with the reason and (for `not_in_star`) the admissible fraction -/
  | refused (reason : String) (frac : Option Float)

/-- Where the end of a trace over the input mesh lies, with snapping onto an input edge (N10, N24:
within `moveSnap` of the longest edge of the face and of `minEdge`, the shortest new edge at the
vertex): the place, the angle of the new frame in the layout of the end face (0 for a face, the
direction of the edge from its lower-numbered end for an edge), whether it was snapped. -/
def placeOfEnd (t0 : Tri) (f : Nat) (b : Float × Float × Float) (minEdge : Float) :
    Loc × Float × Bool := Id.run do
  let P := layoutFace t0 t0.fh[f]!
  let hs := t0.faceHalfedges f
  let qx := b.1 * P[0]!.1 + b.2.1 * P[1]!.1 + b.2.2 * P[2]!.1
  let qy := b.1 * P[0]!.2 + b.2.1 * P[1]!.2 + b.2.2 * P[2]!.2
  let mut longest := fzero
  for k in [0:3] do longest := max longest t0.len[hs[k]!]!
  let mut best := 3
  let mut bestD := posInf
  for k in [0:3] do
    let A := P[k]!
    let B := P[(k + 1) % 3]!
    let ex := B.1 - A.1
    let ey := B.2 - A.2
    let n := Float.sqrt (ex * ex + ey * ey)
    let dist := (ex * (qy - A.2) - ey * (qx - A.1)).abs / n
    if dist < bestD then
      bestD := dist
      best := k
  if bestD ≤ moveSnap * longest && bestD ≤ moveSnap * minEdge then
    let A := P[best]!
    let B := P[(best + 1) % 3]!
    let ex := B.1 - A.1
    let ey := B.2 - A.2
    let s := (ex * (qx - A.1) + ey * (qy - A.2)) / (ex * ex + ey * ey)
    let u := t0.tail[hs[best]!]!
    let w := t0.tail[hs[(best + 1) % 3]!]!
    if u < w then return (.edge u w s, angleOf ex ey, true)
    else return (.edge w u (1 - s), angleOf (-ex) (-ey), true)
  return (.face f b, 0.0, false)

/-- Move the inserted vertex `v` (placed at `loc` on the input mesh `t0`) by the displacement `d` in
its frame (SPEC N19). `nInput` is the number of input vertices. -/
def moveVertex (t0 t : Tri) (nInput v : Nat) (loc : Loc) (d : Float × Float) :
    Except String MoveOutcome := do
  if v < nInput then return .refused "input_vertex" none
  if v ≥ t.nv then throw s!"vertex {v} does not exist"
  let (dx, dy) := d
  if !(dx.isFinite && dy.isFinite) then return .refused "not_finite" none
  let hs := t.outgoing v
  for h in hs do
    if t.head h == v then return .refused "self_edge" none
  if t.isB[v]! && dy != 0 then return .refused "leaves_boundary" none
  if dx == 0 && dy == 0 then return .ok t loc false 0.0
  -- the star
  let (frac, inside) := starKernel t v d
  if !inside then return .refused "not_in_star" (some frac)
  -- new lengths (it-move-length); `α` is the angle from the motion to the edge
  let r := planeDistance 0.0 0.0 dx dy
  let psi := Float.atan2 dy dx
  let newLen := fun (l phi : Float) => movedLength l r (Float.cos (phi - psi)) (Float.sin (phi - psi))
  -- shape of the new triangles
  for h in hs do
    let (j, k, a2) := starTriangle t h
    let side := (k.1 - j.1, k.2 - j.2)
    let na2 := movedDoubleArea a2 d side
    let lvj := newLen t.len[h]! t.sign[h]!
    let lvk := newLen t.len[t.prev h]! (t.sign[h]! + t.corner h)
    let longer := max lvj lvk
    let longest := max longer t.len[t.next[h]!]!
    if na2 < moveShapeLimit * longer * longest then return .refused "too_thin" none
  -- the shortest new edge at the vertex (N24)
  let mut minEdge := posInf
  for h in hs do minEdge := min minEdge (newLen t.len[h]! t.sign[h]!)
  for g in [0:t.nh] do
    if t.head g == v && !t.hasTwin g then minEdge := min minEdge (newLen t.len[g]! floatPi)
  -- the new place and the rotation `rho` from the old frame to the new one
  let mut newLoc := loc
  let mut rho := fzero
  let mut snapped := false
  match loc with
  | .vertex _ => throw s!"vertex {v} is an input vertex"
  | .edge lo hi tt =>
    let hl := inputHalfedge t0 lo hi
    let hr := inputHalfedge t0 hi lo
    let l := match hl, hr with
      | some x, _ => t0.len[x]!
      | none, some x => t0.len[x]!
      | none, none => 0.0
    if l == 0 then throw s!"edge {lo} {hi} is not an input edge"
    if dy == 0 then
      -- along the edge (interior: the frame is lo → hi; boundary: the boundary halfedge)
      let sgn := if t.isB[v]! then (if hl.isSome && t0.twin[hl.get!]! == noTwin then 1.0 else -1.0) else 1.0
      let t' := tt + sgn * dx / l
      if !(t' > 0 && t' < 1) then return .refused "not_located" none
      newLoc := .edge lo hi t'
    else
      let (h0, start, dir) := if dy > 0 then (hl, (tt * l, 0.0), (dx / r, dy / r))
        else (hr, ((1 - tt) * l, 0.0), (-dx / r, -dy / r))
      match h0 with
      | none => return .refused "not_located" none
      | some h0 =>
        let res ← traceCore t0 h0 start dir r 100000 (enter0 := 0)
        if res.stop == "vertex" then return .refused "on_input_vertex" none
        if res.stop != "end" then return .refused "not_located" none
        let (pl, frameAngle, sn) := placeOfEnd t0 res.face res.bary minEdge
        newLoc := pl
        snapped := sn
        rho := angleOf res.dirEnd.1 res.dirEnd.2 - frameAngle - psi
  | .face f b =>
    let res ← traceFromFace t0 f (facePoint t0 f b) (dx / r, dy / r) r
    if res.stop == "vertex" then return .refused "on_input_vertex" none
    if res.stop != "end" then return .refused "not_located" none
    let (pl, frameAngle, sn) := placeOfEnd t0 res.face res.bary minEdge
    newLoc := pl
    snapped := sn
    rho := angleOf res.dirEnd.1 res.dirEnd.2 - frameAngle - psi
  -- lengths
  let mut u := t
  for h in hs do
    let l := newLen t.len[h]! t.sign[h]!
    u := { u with len := u.len.set! h l }
    if t.hasTwin h then u := { u with len := u.len.set! t.twin[h]! l }
  -- the boundary edge that arrives at a boundary vertex (direction π, no halfedge leaving `v`)
  for g in [0:t.nh] do
    if t.head g == v && !t.hasTwin g then
      u := { u with len := u.len.set! g (newLen t.len[g]! floatPi) }
  -- outgoing signposts: the directions of the neighbours from the new place, in the new frame
  let mut sg := u.sign
  for h in hs do
    let phi := t.sign[h]!
    let jx := t.len[h]! * Float.cos phi - dx
    let jy := t.len[h]! * Float.sin phi - dy
    let beta := Float.atan2 jy jx
    let s := if t.isB[v]! then (if beta < 0 then beta + 2 * floatPi else beta) else wrapTurn (beta + rho)
    -- at a boundary vertex the boundary halfedge that leaves it stays at 0 exactly
    sg := sg.set! h (if t.isB[v]! && !t.hasTwin h then 0.0 else s)
  u := { u with sign := sg }
  -- incoming: the side of the link that leaves the neighbour, by the new corner there
  for h in hs do
    if t.hasTwin h then
      let s := t.next[h]!
      let j := t.tail[s]!
      sg := sg.set! t.twin[h]! (signpostUpdate t.isB[j]! t.sign[s]! (u.corner s) (t.theta.get! j))
  u := { u with sign := sg }
  return .ok u newLoc snapped rho

/-- The update of the optimal Delaunay triangulation at `v` (SPEC N21, N23) in the frame of `v`:
`Σ_T |T| (c_T − v) / Σ_T |T|` (`weightedCircumcenter`, rotated by the signpost of the triangle's
halfedge leaving `v`, `odtDisplacement` with the step 1), and the same with the barycentre
(`weightedBarycenter`) for the triangles that have a vertex on the boundary. Valid for a flat vertex
(an inserted one, whose signposts are unscaled angles). -/
def odtTarget (t : Tri) (v : Nat) : (Float × Float) × (Float × Float) := Id.run do
  let mut sc : Float × Float := (0.0, 0.0)
  let mut sb : Float × Float := (0.0, 0.0)
  let mut total := fzero
  for h in t.outgoing v do
    let lab := t.len[h]!
    let lac := t.len[t.prev h]!
    let lbc := t.len[t.next[h]!]!
    let A := triArea lab lbc lac
    let wc := rotateBy t.sign[h]! (weightedCircumcenter lab lac lbc A)
    let atB := t.isB[v]! || t.isB[t.head h]! || t.isB[t.head t.next[h]!]!
    let wb := if atB then rotateBy t.sign[h]! (weightedBarycenter lab lac lbc A) else wc
    sc := (sc.1 + wc.1, sc.2 + wc.2)
    sb := (sb.1 + wb.1, sb.2 + wb.2)
    total := total + A
  return (odtDisplacement sc total 1.0, odtDisplacement sb total 1.0)

/-- The point of the place `loc` in the layout of the input face `f` (`layoutFace t0 (fh f)`), if `loc`
lies in the closed face `f`. -/
def locPointIn (t0 : Tri) (loc : Loc) (f : Nat) : Option (Float × Float) := Id.run do
  let P := layoutFace t0 t0.fh[f]!
  let hs := t0.faceHalfedges f
  match loc with
  | .face g b => if g == f then return some (facePoint t0 f b) else return none
  | .vertex w =>
    for k in [0:3] do
      if t0.tail[hs[k]!]! == w then return some P[k]!
    return none
  | .edge lo hi tt =>
    for k in [0:3] do
      let u := t0.tail[hs[k]!]!
      let w := t0.tail[hs[(k + 1) % 3]!]!
      if (u == lo && w == hi) || (u == hi && w == lo) then
        let s := if u == lo then tt else 1 - tt
        let A := P[k]!
        let B := P[(k + 1) % 3]!
        return some (A.1 + s * (B.1 - A.1), A.2 + s * (B.2 - A.2))
    return none

/-- Trace over the input mesh from a place `loc` (a face or an edge) in the direction `phi` of its frame
for `dist`. -/
def traceFromLoc (t0 : Tri) (loc : Loc) (phi dist : Float) : Except String TraceResult := do
  let dir := (Float.cos phi, Float.sin phi)
  match loc with
  | .face f b => traceFromFace t0 f (facePoint t0 f b) dir dist
  | .edge lo hi tt =>
    let hl := inputHalfedge t0 lo hi
    let hr := inputHalfedge t0 hi lo
    if dir.2 ≥ 0 then
      match hl with
      | some h => traceCore t0 h (tt * t0.len[h]!, 0.0) dir dist 100000 (enter0 := 0)
      | none => throw "no face on the left of the edge"
    else
      match hr with
      | some h => traceCore t0 h ((1 - tt) * t0.len[h]!, 0.0) (-dir.1, -dir.2) dist 100000 (enter0 := 0)
      | none => throw "no face on the right of the edge"
  | .vertex _ => throw "traceFromLoc: start at a vertex"

end Ddg.Intrinsic
