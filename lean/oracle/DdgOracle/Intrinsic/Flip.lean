import DdgOracle.Intrinsic.Tri

/-!
# Edge flip and Delaunay flipping (SPEC §4, N5, N6, N9)
-/

namespace Ddg.Intrinsic

/-- Why a flip is refused. -/
inductive Refusal where
  /-- the edge is on the boundary -/
  | boundary
  /-- an end point has degree one -/
  | degreeOne
  /-- the two triangles do not form a strictly convex quadrilateral -/
  | notConvex
  /-- a new triangle would be thinner than the shape limit (N5) -/
  | tooThin
  deriving BEq

/-- Name of a refusal. -/
def Refusal.name : Refusal → String
  | .boundary => "boundary"
  | .degreeOne => "degree_one"
  | .notConvex => "not_convex"
  | .tooThin => "too_thin"

/-- The shape limit of N5. -/
def flipShapeLimit : Float := 1e-6

/-- The data of a flip. -/
structure FlipInfo where
  /-- vertices `i, j` of the edge and the apexes `k` (left of `i → j`), `l` -/
  i : Nat
  /-- see `i` -/
  j : Nat
  /-- see `i` -/
  k : Nat
  /-- see `i` -/
  l : Nat
  /-- the five lengths `a, b, c, d, e` -/
  lens : Array Float
  /-- the unfolded diamond: `x_k, h_k, x_l, h_l` -/
  diamond : Array Float
  /-- the new length -/
  lnew : Float

/-- Flip the edge of halfedge `h` (SPEC Alg. 5): the halfedge `h` becomes `l → k` in the face that
keeps the id of its face, its twin `k → l`; see the module documentation. -/
def flipHalfedge (t : Tri) (h : Nat) : Except Refusal (Tri × FlipInfo) := do
  let tw := t.twin[h]!
  if tw == noTwin then throw .boundary
  let h1 := t.next[h]!
  let h2 := t.next[h1]!
  let t1 := t.next[tw]!
  let t2 := t.next[t1]!
  if t.twin[h2]! == t1 || t.twin[h1]! == t2 then throw .degreeOne
  let i := t.tail[h]!
  let j := t.tail[h1]!
  let k := t.tail[h2]!
  let l := t.tail[t2]!
  let a := t.len[h]!
  let b := t.len[h2]!
  let c := t.len[h1]!
  let d := t.len[t1]!
  let e := t.len[t2]!
  let xk := diamondX a b c
  let hk := diamondH a b c
  let xl := diamondX a d e
  let hl := diamondH a d e
  if !flipConvex a xk hk xl hl then throw .notConvex
  let lnew := flippedLength xk hk xl hl
  let a1 := flipArea1 xk hk xl hl
  let a2 := flipArea2 a xk hk xl hl
  if a1 < flipShapeLimit * lnew * max d (max b lnew) || a2 < flipShapeLimit * lnew * max c (max e lnew) then
    throw .tooThin
  let f := t.face[h]!
  let g := t.face[tw]!
  -- connectivity: f = (i, l, k): t1, h, h2 ; g = (j, k, l): h1, tw, t2
  let mut tl := t
  tl := { tl with tail := (tl.tail.set! h l).set! tw k }
  tl := { tl with next := ((((((tl.next.set! t1 h).set! h h2).set! h2 t1).set! h1 tw).set! tw t2).set! t2 h1) }
  tl := { tl with face := ((((((tl.face.set! t1 f).set! h f).set! h2 f).set! h1 g).set! tw g).set! t2 g) }
  tl := { tl with fh := (tl.fh.set! f h).set! g tw }
  tl := { tl with len := (tl.len.set! h lnew).set! tw lnew }
  -- signposts of the two new halfedges (Alg. 5, UpdateSignpost)
  let thetaL := cornerAngle c lnew e      -- angle at l in the triangle (j, k, l) between l → j and l → k
  let thetaK := cornerAngle d b lnew      -- angle at k in the triangle (i, l, k) between k → i and k → l
  let signH := signpostUpdate t.isB[l]! (t.sign[t2]!) thetaL (t.theta.get! l)
  let signT := signpostUpdate t.isB[k]! (t.sign[h2]!) thetaK (t.theta.get! k)
  tl := { tl with sign := (tl.sign.set! h signH).set! tw signT }
  return (tl, { i, j, k, l, lens := #[a, b, c, d, e], diamond := #[xk, hk, xl, hl], lnew })

/-- The halfedge `a → b` of the triangulation if exactly one exists. -/
def findHalfedge (t : Tri) (a b : Nat) : Except String Nat := do
  let mut found : Array Nat := #[]
  for h in [0:t.nh] do
    if t.tail[h]! == a && t.head h == b then found := found.push h
  if found.size == 0 then throw s!"there is no edge from {a} to {b}"
  if found.size > 1 then throw s!"{found.size} edges from {a} to {b}: the edge is not determined by its end points"
  return found[0]!

/-- The Delaunay criterion. -/
inductive Criterion where
  /-- flip when `θ_k + θ_l > π + tol` (the default, N6) -/
  | angleSum
  /-- flip when the cotan weight `½ (cot θ_k + cot θ_l) < −1e-6` (the reference's test) -/
  | cotan
  deriving BEq

/-- Result of Delaunay flipping. -/
structure DelaunayReport where
  /-- flips performed -/
  nFlips : Nat
  /-- edges that failed the test and whose flip was refused -/
  nRefused : Nat
  /-- interior edges within the tolerance of the test (ties) -/
  nTies : Nat
  /-- `done`, `flip_limit` or `blocked` -/
  status : String

/-- Does the interior edge `h` fail the test? Returns `(fails, excess-or-weight)`. -/
def failsTest (t : Tri) (h : Nat) (crit : Criterion) (tol : Float) : Bool × Float :=
  let tw := t.twin[h]!
  let a := t.len[h]!
  let b := t.len[t.prev h]!
  let c := t.len[t.next[h]!]!
  let d := t.len[t.next[tw]!]!
  let e := t.len[t.prev tw]!
  match crit with
  | .angleSum =>
    let (thK, thL) := oppositeAngles a b c d e
    let ex := delaunayExcess thK thL
    (!delaunayHolds tol thK thL, ex)
  | .cotan =>
    let w := cotanWeightFromLengths a b c d e / 2
    (w < -1e-6, w)

/-- Delaunay flipping with the plain queue algorithm (SPEC §4.1, N6): all edges in ascending order of
their id, first in first out, the four neighbours of a flipped edge appended unless queued. -/
def delaunayFlip (t0 : Tri) (crit : Criterion := .angleSum) (tol : Float := 1e-10)
    (limitPerEdge : Nat := 64) : Tri × DelaunayReport := Id.run do
  let mut t := t0
  let mut queue : Array Nat := #[]
  let mut inQ : Array Bool := Array.replicate t.nh false
  for h in [0:t.nh] do
    if t.eid[h]! == h then
      queue := queue.push h
      inQ := inQ.set! h true
  let limit := limitPerEdge * t.numEdges
  let mut nFlips := 0
  let mut nRef := 0
  let mut status := "done"
  let mut qi := 0
  while qi < queue.size do
    let e := queue[qi]!
    qi := qi + 1
    inQ := inQ.set! e false
    if t.hasTwin e then
      let (fails, _) := failsTest t e crit tol
      if fails then
        if nFlips ≥ limit then
          status := "flip_limit"
          qi := queue.size
        else
          let h1 := t.next[e]!
          let h2 := t.next[h1]!
          let tw := t.twin[e]!
          let t1 := t.next[tw]!
          let t2 := t.next[t1]!
          match flipHalfedge t e with
          | .ok (t', _) =>
            t := t'
            nFlips := nFlips + 1
            for n in [h1, h2, t1, t2] do
              let id := t.eid[n]!
              if !inQ[id]! then
                queue := queue.push id
                inQ := inQ.set! id true
          | .error _ => nRef := nRef + 1
  if status == "done" && nRef > 0 then status := "blocked"
  -- ties
  let mut nTies := 0
  for h in [0:t.nh] do
    if t.eid[h]! == h && t.hasTwin h then
      let (_, ex) := failsTest t h crit tol
      if ex.abs ≤ (match crit with | .angleSum => tol | .cotan => 1e-6) then nTies := nTies + 1
  return (t, { nFlips, nRefused := nRef, nTies, status })

end Ddg.Intrinsic
