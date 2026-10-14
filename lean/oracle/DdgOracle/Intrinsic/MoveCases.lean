import DdgOracle.Intrinsic.Cases
import DdgOracle.Selftest

/-!
# Cases and checks of the task `move` (SPEC N19, N21, N24)

Fixture cases are generated here: the flips around the new vertex are chosen by trying the edges
opposite it, and on flat meshes the displacements are computed from targets in the plane of the mesh
(the frame of the vertex is the layout of its input face or the direction of its input edge, which on a
flat mesh is a rotation of the plane). `moveChecks` are the self-tests.
-/

namespace Ddg.Intrinsic

open Ddg Ddg.Cuts

/-- The state after the insertion (or split) and the flips of a case: input, triangulation, place. -/
def moveSetup (c : ICase) : Except String (Tri × Tri × Loc) := do
  let (t0, _) ← buildInput c
  let (t1, loc) ← match c.insert, c.split with
    | some (f, b0, b1, b2), _ => insertInFace t0 f (b0, b1, b2)
    | none, some (a, b, tt) => splitEdge t0 a b tt
    | none, none => throw "no insertion"
  let t ← applyFlips t1 c.preflips
  return (t0, t, loc)

/-- Flip the edges opposite `v` (each at most once, up to `maxFlips`): raises the degree of `v`.
Returns the pairs of end points of the flipped edges, each of which identifies one edge when it is
flipped. -/
def flipsAround (t0 : Tri) (v : Nat) (maxFlips : Nat) : Array (Nat × Nat) := Id.run do
  let mut t := t0
  let mut out : Array (Nat × Nat) := #[]
  let mut tried : Array (Nat × Nat) := #[]
  for _ in [0:maxFlips] do
    let mut done := false
    for h in t.outgoing v do
      if !done then
        let s := t.next[h]!
        let a := t.tail[s]!
        let b := t.head s
        if t.hasTwin s && !tried.contains (a, b) && !tried.contains (b, a) && a != v && b != v then
          tried := tried.push (a, b)
          match applyFlips t #[(a, b)] with
          | .ok t' =>
            t := t'
            out := out.push (a, b)
            done := true
          | .error _ => pure ()
    if !done then return out
  return out

/-- Position of a place on the input mesh (positions of the input vertices `pos`). -/
def locPos (m : RawMesh) (c : ICase) (loc : Loc) : Point Float :=
  let p := m.pos
  match loc with
  | .vertex v => p[v]!
  | .edge lo hi t => ⟨p[lo]!.x + t * (p[hi]!.x - p[lo]!.x), p[lo]!.y + t * (p[hi]!.y - p[lo]!.y),
      p[lo]!.z + t * (p[hi]!.z - p[lo]!.z)⟩
  | .face f (b0, b1, b2) =>
    let fc := c.faces[f]!
    ⟨b0 * p[fc.a]!.x + b1 * p[fc.b]!.x + b2 * p[fc.c]!.x, b0 * p[fc.a]!.y + b1 * p[fc.b]!.y + b2 * p[fc.c]!.y,
     b0 * p[fc.a]!.z + b1 * p[fc.b]!.z + b2 * p[fc.c]!.z⟩

/-- On a mesh in the plane `z = 0` (counter-clockwise faces): the angle of the frame of a place in
the plane. Face: its corner 0 → corner 1; interior edge: lower → higher end; boundary place: along the
boundary halfedge of the input (the surface on its left). -/
def frameAngle (m : RawMesh) (c : ICase) (t0 : Tri) (onB : Bool) (loc : Loc) : Float :=
  let p := m.pos
  let ang := fun (a b : Nat) => Float.atan2 (p[b]!.y - p[a]!.y) (p[b]!.x - p[a]!.x)
  match loc with
  | .vertex _ => 0.0
  | .edge lo hi _ =>
    if onB then
      match inputHalfedge t0 lo hi with
      | some h => if t0.twin[h]! == noTwin then ang lo hi else ang hi lo
      | none => ang hi lo
    else ang lo hi
  | .face f _ => ang c.faces[f]!.a c.faces[f]!.b

/-- A displacement in the frame of the vertex that moves it to the plane point `(x, y)`. -/
def towards (m : RawMesh) (c : ICase) (t0 : Tri) (onB : Bool) (loc : Loc) (x y : Float) : Float × Float :=
  let q := locPos m c loc
  let a := frameAngle m c t0 onB loc
  rotateBy (-a) (x - q.x, y - q.y)

/-- Plan a sequence of moves of the new vertex of a case on a flat mesh: each target is a plane
point, or `none` for a raw displacement given in the frame. Simulates the moves. -/
def planFlat (m : RawMesh) (c : ICase) (targets : Array (Sum (Float × Float) (Nat × Float × Float))) :
    Array (Nat × Float × Float) := Id.run do
  match moveSetup c with
  | .error _ => return #[]
  | .ok (t0, t1, loc0) =>
    let v := t0.nv
    let mut t := t1
    let mut loc := loc0
    let mut out := #[]
    for tg in targets do
      let (w, d) := match tg with
        | .inl (x, y) => (v, towards m c t0 t.isB[v]! loc x y)
        | .inr (w, dx, dy) => (w, (dx, dy))
      out := out.push (w, d.1, d.2)
      let wl := if w == v then loc else Loc.vertex w
      match moveVertex t0 t t0.nv w wl d with
      | .ok (.ok t' l' _ _) =>
        t := t'
        loc := l'
      | _ => pure ()
    return out

/-- Plan moves on any mesh: each move goes in the direction `a` (radians, in the frame of the vertex)
for the fraction `q` of the distance to the border of the star. -/
def planFractions (c : ICase) (moves : Array (Float × Float)) : Array (Nat × Float × Float) := Id.run do
  match moveSetup c with
  | .error _ => return #[]
  | .ok (t0, t1, loc0) =>
    let v := t0.nv
    let mut t := t1
    let mut loc := loc0
    let mut out := #[]
    for (a, q) in moves do
      let e := (Float.cos a, Float.sin a)
      let (frac, _) := starKernel t v e
      let d := (q * frac * e.1, q * frac * e.2)
      out := out.push (v, d.1, d.2)
      match moveVertex t0 t t0.nv v loc d with
      | .ok (.ok t' l' _ _) =>
        t := t'
        loc := l'
      | _ => pure ()
    return out

/-- A flat grid of `n × n` vertices on the unit square, the diagonals alternating. -/
def flatGrid (n : Nat) : RawMesh := gridDisk n n fun _ _ => 0.0

/-- The fixtures of the task `move`: case and description. -/
def moveFixtures : Array (ICase × String) := Id.run do
  let mut fs : Array (ICase × String) := #[]
  -- a flat 6 × 6 grid (spacing 0.2); face 26 has the corners (0.4,0.4), (0.6,0.4), (0.6,0.6)
  let g := flatGrid 6
  let cg0 : ICase := { mkICase "it-move-grid-face" "move" g with insert := some (24, 0.5, 0.3, 0.2) }
  let cg := match moveSetup cg0 with
    | .ok (t0, t, _) => { cg0 with preflips := flipsAround t t0.nv 3 }
    | .error _ => cg0
  let (x0, y0) : Float × Float := match moveSetup cg with
    | .ok (_, _, loc) =>
      let q : Point Float := locPos g cg loc
      (q.x, q.y)
    | .error _ => (0.0, 0.0)
  let tg := fun (x y : Float) => (Sum.inl (x, y) : Sum (Float × Float) (Nat × Float × Float))
  -- the star after the flips is the polygon 2, 9, 15, 21, 14, 8; its kernel is
  -- 0.4 ≤ x ≤ 0.6, x − 0.4 ≤ y ≤ x, which meets the input faces 24, 15, 14 and 5
  fs := fs.push ({ cg with moves := planFlat g cg #[tg 0.52 0.46, tg 0.55 0.35, tg 0.45 0.25, tg 0.48 0.15,
      tg x0 y0] },
    "flat grid: a vertex inserted in face 24, three flips (degree 6), moves inside the face, into the input faces 15, 14, 5 and back to the start")
  -- onto an input edge, along it, off it
  let ce := { cg with name := "it-move-grid-edge" }
  let v0 := g.pos.size
  let raw := fun (dx dy : Float) => (Sum.inr (v0, dx, dy) : Sum (Float × Float) (Nat × Float × Float))
  fs := fs.push ({ ce with moves := planFlat g ce #[tg 0.52 0.4, raw 0.04 0.0, tg 0.53 0.47, tg 0.45 0.35,
      raw (0.03 * Float.sqrt 2) 0.0, tg 0.52 0.33] },
    "flat grid: onto the input edge 14-15 (snapped), along it, off it into face 24, onto the input edge 9-14 (snapped), along it, off it into face 15")
  -- refusals and a sequence after them
  let cr := { cg with name := "it-move-refused" }
  let infF := posInf
  let cr0 := match moveSetup cr with
    | .ok (t0, t, _) =>
      let v := t0.nv
      let e := (Float.cos 0.7, Float.sin 0.7)
      let (frac, _) := starKernel t v e
      #[Sum.inr (v, 0.5, 0.2), Sum.inr (0, 0.01, 0.0), Sum.inr (v, infF, 0.0),
        Sum.inr (v, frac * (1 - 1e-9) * e.1, frac * (1 - 1e-9) * e.2),
        Sum.inr (v, 0.0, 0.0), tg (x0 + 0.01) (y0 + 0.01)]
    | .error _ => #[]
  fs := fs.push ({ cr with moves := planFlat g cr cr0 },
    "flat grid: refused moves: outside the star (not_in_star with the fraction), an input vertex, an infinite component, a move that leaves a triangle too thin; then a zero move and an accepted one")
  -- curved: a noisy icosphere, degree 6 after flips
  let ico := (icosphere 2).jitter 0.05 3
  let ci0 : ICase := { mkICase "it-move-ico" "move" ico with insert := some (17, 0.2, 0.5, 0.3) }
  let ci := match moveSetup ci0 with
    | .ok (t0, t, _) => { ci0 with preflips := flipsAround t t0.nv 4 }
    | .error _ => ci0
  fs := fs.push ({ ci with moves := planFractions ci #[(0.3, 0.9), (2.0, 0.85), (4.0, 0.9), (5.5, 0.8), (1.2, 0.5)] },
    "noisy icosphere: a vertex inserted in a face, flips to degree >= 5, five moves in sequence across input faces")
  -- an interior input edge split, then moved off it (curved grid)
  let bump := (gridDisk 6 5 fun x y => 0.3 * Float.sin (3 * x) * Float.cos (2 * y)).jitter 0.01 5
  let cs0 : ICase := { mkICase "it-move-split" "move" bump with split := some (8, 14, 0.4) }
  let cs := match moveSetup cs0 with
    | .ok (t0, t, _) => { cs0 with preflips := flipsAround t t0.nv 2 }
    | .error _ => cs0
  fs := fs.push ({ cs with moves := planFractions cs #[(0.0, 0.5), (4.4, 0.5), (1.3, 0.5), (3.6, 0.6), (5.0, 0.3)] },
    "bumpy grid: a vertex on an interior input edge (split at 0.4), moved along it, off it to the right (into the face on the right of the edge from its lower end), and on")
  -- boundary
  let cb0 : ICase := { mkICase "it-move-boundary" "move" bump with split := some (2, 3, 0.3) }
  fs := fs.push ({ cb0 with moves := #[(30, 0.05, 0.0), (30, -0.08, 0.0), (30, 0.0, 0.01), (30, 0.3, 0.0),
      (30, 0.04, 0.0)] },
    "bumpy grid: a vertex on a boundary edge (split at 0.3) moved along the boundary both ways; refused: off the boundary, beyond its neighbour")
  return fs

/-- Distance of two places on the input mesh, measured in the layout of an input face that contains
both (`posInf` if there is none among the faces of the two places). -/
def locDistance (t0 : Tri) (a b : Loc) : Float := Id.run do
  let faces := fun (l : Loc) => match l with
    | .face f _ => #[f]
    | .edge lo hi _ => Id.run do
      let mut out := #[]
      for x in [inputHalfedge t0 lo hi, inputHalfedge t0 hi lo] do
        if let some h := x then out := out.push t0.face[h]!
      return out
    | .vertex _ => #[]
  for f in faces a ++ faces b do
    match locPointIn t0 a f, locPointIn t0 b f with
    | some p, some q => return planeDistance p.1 p.2 q.1 q.2
    | _, _ => pure ()
  return posInf

/-- The checks of `move` (selftest). -/
def moveChecks (tally : IO.Ref Tally) : IO Unit := do
  let pi := floatPi
  -- text format of the task: print, parse, print
  let mut rt := true
  for (c, _) in moveFixtures do
    let text := printCase c
    let lines := (text.splitOn "\n").filterMap fun l =>
      let t := tokens l
      if t.size == 0 || t[0]! == "case" || t[0]! == "end" then none else some t
    match parseBody c.name lines.toArray with
    | .ok c2 => if printCase c2 != text then rt := false
    | .error _ => rt := false
  check tally s!"move: text format of the {moveFixtures.size} fixtures: print, parse, print is the identity" rt
  let meshes : List (String × RawMesh × Bool) :=
    [("flat grid", flatGrid 7, true),
     ("sheared flat grid", (gridDisk 7 6 fun _ _ => 0.0).mapPos fun p => ⟨p.x + 0.4 * p.y, 0.8 * p.y, 0.0⟩, true),
     ("noisy icosphere", (icosphere 2).jitter 0.05 3, false),
     ("bumpy grid", (gridDisk 7 6 fun x y => 0.3 * Float.sin (3 * x) * Float.cos (2 * y)).jitter 0.01 5, false)]
  for (nm, m, flat) in meshes do
    let c0 := mkICase nm "move" m
    let some t00 := triOf c0 | check tally s!"move {nm}: input builds" false
    let mut rng : Rng := ⟨77⟩
    let mut nOk := 0
    let mut nRef := 0
    let mut wArea := fzero
    let mut wSum := fzero
    let mut wSumV := fzero
    let mut wTrace := fzero
    let mut nTraceBad := 0
    let mut wBack := fzero
    let mut wBackLoc := fzero
    let mut wFlat := fzero
    let mut wFlatPos := fzero
    let mut wFrac := fzero
    let mut nFracBad := 0
    let mut wOdt := fzero
    let mut wFd := fzero
    let mut nSnapped := 0
    let mut nCases := 0
    let sums0 := t00.angleSums
    let area0 := t00.totalArea
    for trial in [0:24] do
      -- an insertion in a random face, or a split of a random interior edge (every fourth trial)
      let (f, r1) := rng.below t00.nf
      let (u, r2) := r1.range 0.15 0.7
      let (w, r3) := r2.range 0.15 0.7
      rng := r3
      let s := u + w
      let (b1, b2) := if s < 0.9 then (u, w) else (u * 0.9 / s, w * 0.9 / s)
      let c1 : ICase := if trial % 4 == 3 then
          let h := 3 * f
          let a := t00.tail[h]!
          let b := t00.head h
          { c0 with split := some (a, b, b1 / (b1 + b2)) }
        else { c0 with insert := some (f, 1 - b1 - b2, b1, b2) }
      let c := match moveSetup c1 with
        | .ok (t0, t, _) => { c1 with preflips := flipsAround t t0.nv 4 }
        | .error _ => c1
      match moveSetup c with
      | .error _ => pure ()
      | .ok (t0, tStart, locStart) =>
        nCases := nCases + 1
        let v := t0.nv
        let mut t := tStart
        let mut loc := locStart
        for _ in [0:6] do
          let (a, r4) := rng.range 0.0 (2 * pi)
          let (q, r5) := r4.range 0.1 0.9
          rng := r5
          let e := if t.isB[v]! then (if a < pi then (1.0, 0.0) else (-1.0, 0.0)) else (Float.cos a, Float.sin a)
          let (frac, _) := starKernel t v e
          if frac.isFinite then
            -- kernel fraction: just inside and just outside
            match moveVertex t0 t t0.nv v loc (0.999 * frac * e.1, 0.999 * frac * e.2) with
            | .ok (.refused "not_in_star" _) => nFracBad := nFracBad + 1
            | _ => pure ()
            match moveVertex t0 t t0.nv v loc (1.001 * frac * e.1, 1.001 * frac * e.2) with
            | .ok (.refused "not_in_star" (some g)) => wFrac := max wFrac (g - 1 / 1.001).abs
            | _ => nFracBad := nFracBad + 1
          let d := if frac.isFinite then (q * frac * e.1, q * frac * e.2) else (0.1 * q * e.1, 0.1 * q * e.2)
          match moveVertex t0 t t0.nv v loc d with
          | .error _ => nRef := nRef + 1
          | .ok (.refused _ _) => nRef := nRef + 1
          | .ok (.ok t' loc' sn rho) =>
            nOk := nOk + 1
            if sn then nSnapped := nSnapped + 1
            wArea := max wArea ((t'.totalArea - area0).abs / area0)
            let s' := t'.angleSums
            for i in [0:t0.nv] do wSum := max wSum (s'.get! i - sums0.get! i).abs
            wSumV := max wSumV (s'.get! v - (if t'.isB[v]! then pi else 2 * pi)).abs
            -- traces of the edges at v over the input mesh, both ways (SPEC §5 item 3)
            for h in t'.outgoing v do
              let j := t'.head h
              let l := t'.len[h]!
              if j < t0.nv then
                match traceFromLoc t0 loc' t'.sign[h]! l with
                | .ok r =>
                  let miss := if r.stop == "vertex" && r.vertex == j then 0.0 else missDistance t0 r j
                  if miss ≤ 1e-9 * l then wTrace := max wTrace (miss / l) else nTraceBad := nTraceBad + 1
                | .error _ => nTraceBad := nTraceBad + 1
                if t'.hasTwin h then
                  let g := t'.twin[h]!
                  match traceFromVertex t0 j t'.sign[g]! l with
                  | .ok r =>
                    match locPointIn t0 loc' r.face with
                    | some pnt =>
                      let en := endPoint t0 r
                      let miss := planeDistance en.1 en.2 pnt.1 pnt.2
                      if miss ≤ 1e-9 * l then wTrace := max wTrace (miss / l) else nTraceBad := nTraceBad + 1
                    | none => nTraceBad := nTraceBad + 1
                  | .error _ => nTraceBad := nTraceBad + 1
            -- flat: the new place is the old one plus the displacement, the lengths are straight
            if flat then
              let q0 := locPos m c loc
              let q1 := locPos m c loc'
              let dg := rotateBy (frameAngle m c t0 t.isB[v]! loc) d
              wFlatPos := max wFlatPos (max (q1.x - q0.x - dg.1).abs (q1.y - q0.y - dg.2).abs)
              for h in t'.outgoing v do
                let j := t'.head h
                if j < t0.nv then
                  let pj := m.pos[j]!
                  let dd := Float.sqrt ((pj.x - q1.x) ^ 2 + (pj.y - q1.y) ^ 2)
                  wFlat := max wFlat ((dd - t'.len[h]!).abs / dd)
              -- the ODT update against circumcentres computed in the plane, and against the
              -- gradient of the energy of the star by central differences
              let (oc, _) := odtTarget t' v
              let mut acc : Float × Float := (0.0, 0.0)
              let mut tot := fzero
              let mut star : Array ((Float × Float) × (Float × Float)) := #[]
              let mut allInput := true
              for h in t'.outgoing v do
                let j := t'.head h
                let k := t'.head t'.next[h]!
                if j < t0.nv && k < t0.nv then
                  let pj := (m.pos[j]!.x, m.pos[j]!.y)
                  let pk := (m.pos[k]!.x, m.pos[k]!.y)
                  star := star.push (pj, pk)
                  let A := signedArea2 q1.x q1.y pj.1 pj.2 pk.1 pk.2 / 2
                  let ax := pj.1 - q1.x
                  let ay := pj.2 - q1.y
                  let bx := pk.1 - q1.x
                  let byy := pk.2 - q1.y
                  let dn := 2 * (ax * byy - ay * bx)
                  let ccx := (byy * (ax * ax + ay * ay) - ay * (bx * bx + byy * byy)) / dn
                  let ccy := (ax * (bx * bx + byy * byy) - bx * (ax * ax + ay * ay)) / dn
                  acc := (acc.1 + A * ccx, acc.2 + A * ccy)
                  tot := tot + A
                else allInput := false
              if allInput && !t'.isB[v]! then
                let og := rotateBy (frameAngle m c t0 false loc') oc
                wOdt := max wOdt (max (og.1 - acc.1 / tot).abs (og.2 - acc.2 / tot).abs)
                let energy := fun (x y : Float) => Id.run do
                  let mut e := fzero
                  for (pj, pk) in star do
                    let A := signedArea2 x y pj.1 pj.2 pk.1 pk.2 / 2
                    e := e + (x * x + y * y + pj.1 * pj.1 + pj.2 * pj.2 + pk.1 * pk.1 + pk.2 * pk.2) * A / 3
                  return e
                let hh := 1e-5
                let gx := (energy (q1.x + hh) q1.y - energy (q1.x - hh) q1.y) / (2 * hh)
                let gy := (energy q1.x (q1.y + hh) - energy q1.x (q1.y - hh)) / (2 * hh)
                let pred := (-3 / (2 * tot) * gx, -3 / (2 * tot) * gy)
                wFd := max wFd (max (pred.1 - acc.1 / tot).abs (pred.2 - acc.2 / tot).abs)
            -- move and move back
            let r := planeDistance 0.0 0.0 d.1 d.2
            let back := if t'.isB[v]! then (-d.1, -d.2) else
              let ang := Float.atan2 d.2 d.1 + pi + rho
              (r * Float.cos ang, r * Float.sin ang)
            match moveVertex t0 t' t0.nv v loc' back with
            | .ok (.ok t2 loc2 _ _) =>
              wBack := max wBack (sameTriangulation t t2)
              wBackLoc := max wBackLoc (locDistance t0 loc loc2)
            | _ => wBack := posInf
            t := t'
            loc := loc'
    check tally s!"move {nm}: {nCases} insertions and splits with flips, {nOk} moves accepted ({nRef} refused, {nSnapped} snapped): area and angle sums unchanged, angle sum 2 pi (pi) at the moved vertex"
      (nOk > 20 && wArea ≤ 1e-12 && wSum ≤ 1e-12 && wSumV ≤ 1e-12)
      s!"area {sci wArea} angle sums {sci wSum} at v {sci wSumV}"
    check tally s!"move {nm}: every edge at the moved vertex traced over the input mesh from both ends arrives at the other end"
      (nTraceBad == 0) s!"largest miss / length {sci wTrace}"
    check tally s!"move {nm}: moving back by the same distance in the direction of the old place restores lengths, signposts and the place"
      (wBack ≤ 1e-10 && wBackLoc ≤ 1e-12) s!"triangulation {sci wBack} place {sci wBackLoc}"
    check tally s!"move {nm}: the admissible fraction of N19: 0.999 of it accepted, 1.001 of it refused with the fraction 1/1.001"
      (nFracBad == 0 && wFrac ≤ 1e-12) s!"fraction error {sci wFrac}"
    if flat then
      check tally s!"move {nm}: flat: the new place is the old one plus the displacement, the new lengths are the straight distances, the ODT update is the area-weighted mean of the circumcentres and -3/(2|w|) times the gradient of the energy (central differences)"
        (wFlatPos ≤ 1e-12 && wFlat ≤ 1e-12 && wOdt ≤ 1e-12 && wFd ≤ 1e-7)
        s!"place {sci wFlatPos} lengths {sci wFlat} odt {sci wOdt} gradient {sci wFd}"

end Ddg.Intrinsic
