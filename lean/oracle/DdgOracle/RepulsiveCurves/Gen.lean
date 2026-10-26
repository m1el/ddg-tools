import DdgOracle.Generators
import DdgOracle.RepulsiveCurves.Format

/-!
# Test curves and the seeded case generator of the `rc` oracle

Curves: circle, trefoil, open helix, tripod (3-junction), Hopf link, planar loop, random loops.
Case `i` of seed `s` depends on `(s, i)` only.
-/

namespace Ddg.RepulsiveCurves

private abbrev G := StateM Rng

private def uni (a b : Float) : G Float := modifyGet fun r => r.range a b
private def below (n : Nat) : G Nat := modifyGet fun r => r.below n

private def twoPi : Float := 6.283185307179586

/-- Edges `i → i+1` of a closed chain of `n` vertices starting at `off`. -/
def loopEdges (n off : Nat) : Array (Nat × Nat) :=
  (Array.range n).map fun i => (off + i, off + (i + 1) % n)

/-- Regular `n`-gon of radius `r` in the plane `z = 0`. -/
def circle (n : Nat) (r : Float) : Net :=
  { pos := (Array.range n).map fun i =>
      let a := twoPi * i.toFloat / n.toFloat
      ⟨r * Float.cos a, r * Float.sin a, 0.0⟩
    edges := loopEdges n 0 }

/-- The trefoil `(sin t + 2 sin 2t, cos t − 2 cos 2t, −sin 3t)` scaled by `s`, `n` vertices. -/
def trefoil (n : Nat) (s : Float) : Net :=
  { pos := (Array.range n).map fun i =>
      let t := twoPi * i.toFloat / n.toFloat
      ⟨s * (Float.sin t + 2.0 * Float.sin (2.0 * t)), s * (Float.cos t - 2.0 * Float.cos (2.0 * t)),
       -(s * Float.sin (3.0 * t))⟩
    edges := loopEdges n 0 }

/-- An open helical arc with `n` vertices (`n - 1` edges), made irregular by `noise`. -/
def helix (n : Nat) (noise : Float) : G Net := do
  let mut pos : Array V3 := #[]
  for i in [0:n] do
    let t := 1.7 * i.toFloat / n.toFloat * twoPi
    let dx ← uni (-noise) noise
    let dy ← uni (-noise) noise
    let dz ← uni (-noise) noise
    pos := pos.push ⟨Float.cos t + dx, Float.sin t + dy, 0.35 * t + dz⟩
  return { pos, edges := (Array.range (n - 1)).map fun i => (i, i + 1) }

/-- A tripod: vertex 0 joined to three arms of `len` vertices going out in three directions. -/
def tripod (len : Nat) (noise : Float) : G Net := do
  let mut pos : Array V3 := #[V3.zero]
  let mut edges : Array (Nat × Nat) := #[]
  let dirs : Array V3 := #[⟨1.0, 0.0, 0.1⟩, ⟨-0.5, 0.866, -0.1⟩, ⟨-0.5, -0.866, 0.15⟩]
  for a in [0:3] do
    for j in [0:len] do
      let dx ← uni (-noise) noise
      let dy ← uni (-noise) noise
      let dz ← uni (-noise) noise
      let r := 0.4 * (j + 1).toFloat
      let bend := 0.08 * (j + 1).toFloat * (j + 1).toFloat
      pos := pos.push (r * dirs[a]! + ⟨dx, dy, dz + bend * (a.toFloat - 1.0)⟩)
      edges := edges.push (if j == 0 then 0 else 1 + a * len + j - 1, 1 + a * len + j)
  return { pos, edges }

/-- Two linked circles of `n` vertices each: one in the plane `z = 0` around the origin, one in
the plane `y = 0` around `(1, 0, 0)`. -/
def hopf (n : Nat) (noise : Float) : G Net := do
  let mut pos : Array V3 := #[]
  for i in [0:n] do
    let a := twoPi * i.toFloat / n.toFloat
    let d ← uni (-noise) noise
    pos := pos.push ⟨(1.0 + d) * Float.cos a, (1.0 + d) * Float.sin a, 0.0⟩
  for i in [0:n] do
    let a := twoPi * i.toFloat / n.toFloat
    let d ← uni (-noise) noise
    pos := pos.push ⟨1.0 + (1.0 + d) * Float.cos a, 0.0, (1.0 + d) * Float.sin a⟩
  return { pos, edges := loopEdges n 0 ++ loopEdges n n }

/-- A closed curve `r(θ)(cos θ, sin θ, h(θ))`, `r = 1 + a sin(kθ) + noise`, heights `h = zamp sin(2θ) + noise`. -/
def wobbleLoop (n : Nat) (a : Float) (k : Nat) (zamp noise : Float) : G Net := do
  let mut pos : Array V3 := #[]
  for i in [0:n] do
    let th := twoPi * (i.toFloat + 0.0) / n.toFloat
    let dr ← uni (-noise) noise
    let dz ← uni (-noise) noise
    let r := 1.0 + a * Float.sin (k.toFloat * th) + dr
    pos := pos.push ⟨r * Float.cos th, r * Float.sin th, zamp * Float.sin (2.0 * th) + dz⟩
  return { pos, edges := loopEdges n 0 }

/-- Move a network by a rigid-looking offset (not needed for correctness; breaks symmetric zeros). -/
def shifted (n : Net) (d : V3) : Net := { n with pos := n.pos.map (· + d) }

/-- Centroid and mean distance to it. -/
def centroidRadius (n : Net) : V3 × Float :=
  let c := (1.0 / n.nv.toFloat) * n.pos.foldl (· + ·) V3.zero
  (c, (n.pos.foldl (fun acc p => acc + (p - c).norm) 0.0) / n.nv.toFloat)

private def perturb (s : V3) (amt : Float) : G V3 := do
  let dx ← uni (-amt) amt
  let dy ← uni (-amt) amt
  let dz ← uni (-amt) amt
  return s + ⟨dx, dy, dz⟩

/-- A random non-redundant constraint list for `net` (SPEC §5.1; the menu keeps the rows independent). -/
def randomConstraints (net0 : Net) (planar : Bool) : G (Net × Array Cons) := do
  let net := net0
  let n := net.nv
  let (c, rad) := centroidRadius net
  let menu ← below (if planar then 6 else 9)
  let sub (k : Nat) : G (Array Nat) := do
    -- k distinct vertices
    let mut vs : Array Nat := #[]
    let mut tries := 0
    while vs.size < k && tries < 100 do
      let v ← below n
      if !vs.contains v then vs := vs.push v
      tries := tries + 1
    return vs.qsort (· < ·)
  match menu with
  | 0 => return (net, #[.barycenter none, .edgeLengths none])
  | 1 => return (net, #[.barycenter none, .totalLength none])
  | 2 =>
    let es ← (do
      let mut es : Array Nat := #[]
      for e in [0:net.ne] do
        let r ← below 2
        if r == 0 then es := es.push e
      pure es)
    return (net, #[.barycenter none, .edgeLengths (some es), .totalLength none])
  | 3 =>
    let v ← below n
    return (net, #[.point v none, .totalLength none])
  | 4 =>
    let vs ← sub (min 4 n)
    -- put the constrained vertices on the sphere, so that the constraint holds at the base curve
    let mut pos := net.pos
    for v in vs do
      let d := pos[v]! - c
      pos := pos.set! v (c + ((rad / d.norm) * d))
    return ({ net with pos }, #[.barycenter none, .sphere c rad (some vs)])
  | 5 =>
    let v ← below n
    let w ← below n
    -- two points on non-adjacent vertices (no edge between them, so no redundant length row)
    let w := if v == w then (v + 2) % n else w
    let adjacent := net.edges.any fun (a, b) => (a == v && b == w) || (a == w && b == v)
    if adjacent then return (net, #[.point v none, .totalLength none])
    return (net, #[.point v none, .point w none])
  | 6 =>
    let e ← below net.ne
    let (a, b) := net.edges[e]!
    let s := (net.pos[b]! - net.pos[a]!)
    let x ← perturb ((1.0 / s.norm) * s) 0.4
    return (net, #[.barycenter none, .totalLength none, .tangent (.edge e) (some x)])
  | 7 =>
    -- a vertex of degree 1 or 2
    let mut cand : Array Nat := #[]
    for v in [0:n] do
      let d := (net.incident v).size
      if d == 1 || d == 2 then cand := cand.push v
    if cand.isEmpty then return (net, #[.barycenter none, .edgeLengths none])
    let v := cand[(← below cand.size)]!
    let exact ← below 2
    let dir ← match tangentOf net net.pos (.vertex v) with
      | .ok (s, _) => if exact == 0 then pure none else do pure (some (← perturb s 0.3))
      | .error _ => pure none
    return (net, #[.barycenter none, .edgeLengths none, .tangent (.vertex v) dir])
  | _ =>
    let e ← below net.ne
    return (net, #[.barycenter none, .totalLength none, .tangent (.edge e) none])

/-- The kinds of curve of the generator. -/
private def randomNet (kind : Nat) : G (Net × Option Nat) := do
  match kind with
  | 0 => let n ← below 10; return (← wobbleLoop (8 + n) 0.25 3 0.2 0.04, none)
  | 1 => let n ← below 12; return (← helix (7 + n) 0.04, none)
  | 2 => let l ← below 3; return (← tripod (3 + l) 0.03, none)
  | 3 => let n ← below 5; return (← hopf (8 + n) 0.05, none)
  | 4 =>
    let n ← below 10
    let net ← wobbleLoop (9 + n) 0.3 3 0.0 0.03
    return ({ net with pos := net.pos.map fun p => ⟨p.x, p.y, 0.0⟩ }, some 2)
  | 5 => let n ← below 12; return (trefoil (14 + n) 0.3, none)
  | 6 =>
    -- planar network in the plane y = 0.5 (held axis 1)
    let n ← below 8
    let net ← wobbleLoop (9 + n) 0.2 2 0.0 0.03
    return ({ net with pos := net.pos.map fun p => ⟨p.x, 0.5, p.y⟩ }, some 1)
  | _ => let n ← below 10; return (← wobbleLoop (10 + n) 0.35 5 0.3 0.05, none)

private def expoChoices : Array Expo := #[⟨3.0, 6.0⟩, ⟨2.0, 4.5⟩, ⟨2.5, 5.5⟩, ⟨3.0, 5.5⟩, ⟨4.0, 7.5⟩, ⟨1.5, 3.6⟩]

/-- Case `i` of seed `s`. -/
def genCase (seed i : Nat) : RCCase :=
  let rng : Rng := ⟨(seed.toUInt64 * 0x100000001B3) ^^^ (i.toUInt64 * 0x9E3779B97F4A7C15 + 12345)⟩
  let act : G RCCase := do
    let kind ← below 8
    let (net, planar) ← randomNet kind
    let ex := expoChoices[(← below expoChoices.size)]!
    let quad := if (← below 2) == 0 then Quad.midpoint else Quad.trapezoid
    let tk ← below 8
    let task := if tk == 0 then "energy" else if tk == 1 then "gram" else if tk < 5 then "saddle" else "project"
    let (net, cons) ← if task == "saddle" || task == "project" then randomConstraints net planar.isSome else pure (net, #[])
    let f ← uni 0.02 0.6
    let tolExp ← uni 2.0 5.0
    let maxit ← below 4
    return { name := s!"rc-{seed}-{i}", task, ex, quad, planar, allMatrices := true, net, cons,
             tau := f, tauRel := true, tolMean := true, tolVal := Float.pow 10.0 (-tolExp), maxit := 1 + maxit }
  (act.run rng).1

end Ddg.RepulsiveCurves
