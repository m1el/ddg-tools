import DdgOracle.Selftest
import DdgOracle.Intrinsic.MoveCases

/-!
# `intrinsic-selftest`
-/

namespace Ddg.Intrinsic

open Ddg Ddg.Cuts

def sheared (nx ny : Nat) : RawMesh :=
  ((gridDisk nx ny fun x y => 0.1 * Float.sin (4 * x) * y).mapPos fun p => ⟨p.x + 2.5 * p.y, p.y * 0.4, p.z⟩).jitter 0.01 5

/-- A torus with `nu × nv` vertices whose rows are staggered by a third of a step: the surface is a
torus, the triangles are thin parallelograms (not Delaunay), and flipping them produces edges from a
vertex to itself and double edges. Major radius `R`, minor radius `r`. -/
def twistedTorus (nu nv : Nat) (R r : Float) : RawMesh := Id.run do
  let idx := fun (i j : Nat) => (i % nu) * nv + j
  let mut pos : Array (Point Float) := Array.replicate (nu * nv) default
  for i in [0:nu] do
    for j in [0:nv] do
      let a := 2 * floatPi * (Float.ofNat i + Float.ofNat j / Float.ofNat nv) / Float.ofNat nu
      let b := 2 * floatPi * Float.ofNat j / Float.ofNat nv
      let d := R + r * Float.cos b
      pos := pos.set! (idx i j) ⟨d * Float.cos a, d * Float.sin a, r * Float.sin b⟩
  let mut faces : Array Face := #[]
  for i in [0:nu] do
    for j in [0:nv] do
      let p00 := idx i j
      let p10 := idx (i + 1) j
      let p01 := if j + 1 < nv then idx i (j + 1) else idx (i + 1) 0
      let p11 := if j + 1 < nv then idx (i + 1) (j + 1) else idx (i + 2) 0
      faces := (faces.push ⟨p00, p10, p11⟩).push ⟨p00, p11, p01⟩
  return { pos, faces }

/-- The meshes of the tests: name and mesh. -/
def testMeshes : List (String × RawMesh) :=
  [("icosphere 2, noisy", (icosphere 2).jitter 0.05 3),
   ("sheared grid with boundary", sheared 9 8),
   ("long thin torus", torus 36 3 2.0 0.06 fun _ _ => 0.0),
   ("twisted thin torus", (twistedTorus 36 3 2.0 0.06).jitter 0.002 17),
   ("thin tube torus (Delaunay has self-edges and double edges)", (twistedTorus 12 3 1.5 0.03).jitter 0.001 3),
   ("annulus", (polar 4 14 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.2 * y * y).jitter 0.03 7),
   ("bumpy grid disk", (gridDisk 9 8 fun x y => 0.5 * Float.sin (5 * x) * Float.cos (4 * y)).jitter 0.02 11)]

/-- Run all checks; returns the number of failures. -/
def selftest : IO Nat := do
  let tally ← IO.mkRef ({} : Tally)
  let pi := floatPi

  -- generic formulas on random triangles
  let mut rng : Rng := ⟨2024⟩
  let mut worstPt := fzero
  let mut worstApex := fzero
  let mut worstCirc := fzero
  let mut worstBary := fzero
  for _ in [0:300] do
    let (x1, r1) := rng.range 0.5 2.0
    let (x2, r2) := r1.range (-1.0) 2.0
    let (y2, r3) := r2.range 0.2 2.0
    let (u, r4) := r3.range 0.1 0.8
    let (v, r5) := r4.range 0.1 0.8
    rng := r5
    -- triangle (0,0), (x1,0), (x2,y2)
    let l01 := x1
    let l12 := Float.sqrt ((x2 - x1) * (x2 - x1) + y2 * y2)
    let l20 := Float.sqrt (x2 * x2 + y2 * y2)
    -- layout reproduces the lengths
    let (ax, ay) := unfoldApex 0.0 0.0 l01 0.0 l01 l12 l20
    worstApex := max worstApex (max (ax - x2).abs (ay - y2).abs)
    -- point from barycentric coordinates -> distances -> position -> barycentric
    let b0 := u
    let b1 := v
    let b2 := 1 - u - v
    let (d0, d1, d2) := cornerDistances l01 l12 l20 b0 b1 b2
    let (px, py) := pointFromDistances l01 l12 l20 d0 d1 d2
    let (c0, c1, c2) := barycentric 0.0 0.0 l01 0.0 ax ay px py
    worstPt := max worstPt (max (c0 - b0).abs (max (c1 - b1).abs (c2 - b2).abs))
    -- the baryLength2 formula of the cut flow gives the same distances
    let dd := baryLength2 l01 l12 l20 (b0 - 1) b1 b2
    worstBary := max worstBary ((Float.sqrt dd - d0).abs)
    -- circumcenter equidistant
    let (w0, w1, w2) := circumcenterBary l12 l20 l01
    let (r0, rr1, rr2) := cornerDistances l01 l12 l20 w0 w1 w2
    worstCirc := max worstCirc (max (r0 - rr1).abs (max (r0 - rr2).abs (w0 + w1 + w2 - 1).abs))
  check tally "formulas on 300 random triangles: unfold apex, barycentric coordinates from distances, point distances against baryLength2, circumcenter equidistant"
    (worstApex ≤ 1e-12 && worstPt ≤ 1e-11 && worstBary ≤ 1e-12 && worstCirc ≤ 1e-11)
    s!"apex {sci worstApex} point round trip {sci worstPt} distances {sci worstBary} circumcenter {sci worstCirc}"

  for (nm, m) in testMeshes do
    let c := mkICase nm "flip" m
    match buildInput c with
    | .error e => check tally s!"{nm}: input builds" false e
    | .ok (t0, topo) =>
      -- angle sums and signposts at the start
      let sums := t0.angleSums
      let mut worstNorm := fzero
      for v in [0:t0.nv] do
        -- the scaled corner angles at a vertex sum to the full angle: follow the halfedges
        let mut total := fzero
        for h in [0:t0.nh] do
          if t0.tail[h]! == v then total := total + t0.corner h * t0.scaleAt v
        let full := if t0.isB[v]! then pi else 2 * pi
        worstNorm := max worstNorm ((total - full).abs)
      check tally s!"{nm}: angle normalization (scaled corners sum to 2 pi / pi at every vertex), {t0.nv} vertices, {topo.numBoundary} on the boundary"
        (worstNorm ≤ 1e-12) s!"max deviation {sci worstNorm}"
      -- flips: twice restores; area and angle sums; flipped length against the unfolded diamond
      let mut tried := 0
      let mut refused := 0
      let mut worstTwice := fzero
      let mut worstArea := fzero
      let mut worstSum := fzero
      let mut worstLen := fzero
      for h in [0:t0.nh] do
        if t0.eid[h]! == h then
          match flipHalfedge t0 h with
          | .error _ => refused := refused + 1
          | .ok (t1, info) =>
            tried := tried + 1
            worstArea := max worstArea ((t1.totalArea - t0.totalArea).abs / t0.totalArea)
            let s1 := t1.angleSums
            for v in [0:t0.nv] do
              worstSum := max worstSum ((s1.get! v - sums.get! v).abs)
            -- unfolded diamond: i = (0,0), j = (a,0), k above, l below
            let a := info.lens[0]!
            let (kx, ky) := unfoldApex 0.0 0.0 a 0.0 a info.lens[2]! info.lens[1]!
            let (lx, ly) := unfoldApex a 0.0 (-a) 0.0 a info.lens[3]! info.lens[4]!
            let direct := planeDistance kx ky lx ly
            worstLen := max worstLen ((direct - info.lnew).abs / info.lnew)
            match flipHalfedge t1 h with
            | .ok (t2, _) =>
              worstTwice := max worstTwice (sameTriangulation t0 t2)
            | .error _ => worstTwice := posInf
      check tally s!"{nm}: {tried} edges flipped ({refused} refused); flipping twice restores lengths and signposts; area and angle sums unchanged; flipped length = distance in the unfolded diamond"
        (worstTwice ≤ 1e-12 && worstArea ≤ 1e-13 && worstSum ≤ 1e-12 && worstLen ≤ 1e-12)
        s!"twice {sci worstTwice} area {sci worstArea} angle sums {sci worstSum} diamond {sci worstLen}"
      -- Delaunay
      let (td, rep) := delaunayFlip t0
      let (ex, w) := delaunayMeasures td
      let sd := td.angleSums
      let mut dsum := fzero
      for v in [0:t0.nv] do dsum := max dsum ((sd.get! v - sums.get! v).abs)
      let mut nonSimplicial := 0
      for h in [0:td.nh] do
        if td.tail[h]! == td.head h then nonSimplicial := nonSimplicial + 1
        for g in [h+1:td.nh] do
          if td.tail[h]! == td.tail[g]! && td.head h == td.head g then nonSimplicial := nonSimplicial + 1
      check tally s!"{nm}: Delaunay flipping ({rep.nFlips} flips, status {rep.status}, {rep.nTies} ties, {nonSimplicial} self-edge halfedges and double-edge pairs): every interior edge satisfies the condition and has a non-negative cotan weight; area and angle sums unchanged"
        (rep.status == "done" && ex ≤ 1e-10 && w ≥ -1e-9 && (td.totalArea - t0.totalArea).abs / t0.totalArea ≤ 1e-12 && dsum ≤ 1e-11)
        s!"max excess {sci ex} min weight {sci w} angle sums {sci dsum}"
      -- tracing the signpost of every edge for its length over the input mesh
      let mut worstMiss := fzero
      let mut nOk := 0
      let mut nBad := 0
      for h in [0:td.nh] do
        match traceFromVertex t0 td.tail[h]! td.sign[h]! td.len[h]! with
        | .error _ => nBad := nBad + 1
        | .ok r =>
          let miss := if r.stop == "vertex" && r.vertex == td.head h then 0.0 else missDistance t0 r (td.head h)
          if miss ≤ 1e-9 * td.len[h]! then
            nOk := nOk + 1
            worstMiss := max worstMiss (miss / td.len[h]!)
          else nBad := nBad + 1
      check tally s!"{nm}: tracing the signpost of every edge of the Delaunay triangulation for its length over the input mesh arrives at the other end ({nOk} of {td.nh} halfedges)"
        (nBad == 0) s!"largest miss / length {sci worstMiss}"
      -- circumcenters
      let mut worstR := fzero
      for f in [0:t0.nf] do
        let hs := t0.faceHalfedges f
        let (w0, w1, w2) := circumcenterBary t0.len[hs[1]!]! t0.len[hs[2]!]! t0.len[hs[0]!]!
        let (r0, r1, r2) := cornerDistances t0.len[hs[0]!]! t0.len[hs[1]!]! t0.len[hs[2]!]! w0 w1 w2
        worstR := max worstR (max (r0 - r1).abs (r0 - r2).abs / r0)
      check tally s!"{nm}: circumcenter of every face equidistant from the three corners" (worstR ≤ 1e-9) s!"max relative difference {sci worstR}"

  -- signpost update wraps
  let w := signpostUpdate false (2 * pi - 0.1) 0.5 (2 * pi)
  check tally "signpost update reduces to [0, 2 pi) at an interior vertex and not at a boundary vertex"
    ((w - 0.4).abs ≤ 1e-12 && (signpostUpdate true 0.5 0.5 pi - 1.0).abs ≤ 1e-12)

  -- text format: parse, print, parse
  let c := { mkICase "rt" "trace" (icosphere 1) with traces := #[(3, 1.5, 0.2)] }
  let text := printCase c
  let ok := match (Id.run do
      let lines := (text.splitOn "\n").filterMap fun l =>
        let t := tokens l
        if t.size == 0 || t[0]! == "case" || t[0]! == "end" then none else some t
      pure (parseBody "rt" lines.toArray)) with
    | .ok c2 => printCase c2 == text
    | .error _ => false
  check tally "text format: print, parse, print is the identity" ok

  -- moving a vertex (task move)
  moveChecks tally

  let t ← tally.get
  IO.println s!"intrinsic-selftest: {t.passed} passed, {t.failed} failed"
  return t.failed

end Ddg.Intrinsic
