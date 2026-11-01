import DdgOracle.RepulsiveSurfaces.Format

/-!
# Test meshes and the seeded case generator of the `rs` oracle

Meshes: tetrahedron, octahedron, icosahedron, torus grids, open grid patches, Möbius bands (the
non-orientable case), two disjoint components. Case `i` of seed `s` depends on `(s, i)` only.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

private abbrev G := StateM Rng

private def uni (a b : Float) : G Float := modifyGet fun r => r.range a b
private def below (n : Nat) : G Nat := modifyGet fun r => r.below n

private def twoPi : Float := 6.283185307179586

private def tri (a b c : Nat) : Tri := ⟨a, b, c⟩

/-- Regular tetrahedron (outward faces). -/
def tetra : Mesh :=
  { pos := #[⟨1.0, 1.0, 1.0⟩, ⟨1.0, -1.0, -1.0⟩, ⟨-1.0, 1.0, -1.0⟩, ⟨-1.0, -1.0, 1.0⟩]
    faces := #[tri 0 1 2, tri 0 3 1, tri 0 2 3, tri 1 3 2] }

/-- Octahedron (outward faces). -/
def octa : Mesh :=
  { pos := #[⟨1.0, 0.0, 0.0⟩, ⟨-1.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩, ⟨0.0, -1.0, 0.0⟩, ⟨0.0, 0.0, 1.0⟩, ⟨0.0, 0.0, -1.0⟩]
    faces := #[tri 0 2 4, tri 2 1 4, tri 1 3 4, tri 3 0 4, tri 2 0 5, tri 1 2 5, tri 3 1 5, tri 0 3 5] }

/-- Icosahedron (outward faces), circumradius about 1.9. -/
def icosa : Mesh :=
  let p := 1.618033988749895
  { pos := #[⟨-1.0, p, 0.0⟩, ⟨1.0, p, 0.0⟩, ⟨-1.0, -p, 0.0⟩, ⟨1.0, -p, 0.0⟩, ⟨0.0, -1.0, p⟩, ⟨0.0, 1.0, p⟩,
             ⟨0.0, -1.0, -p⟩, ⟨0.0, 1.0, -p⟩, ⟨p, 0.0, -1.0⟩, ⟨p, 0.0, 1.0⟩, ⟨-p, 0.0, -1.0⟩, ⟨-p, 0.0, 1.0⟩]
    faces := #[tri 0 11 5, tri 0 5 1, tri 0 1 7, tri 0 7 10, tri 0 10 11, tri 1 5 9, tri 5 11 4, tri 11 10 2,
               tri 10 7 6, tri 7 1 8, tri 3 9 4, tri 3 4 2, tri 3 2 6, tri 3 6 8, tri 3 8 9, tri 4 9 5,
               tri 2 4 11, tri 6 2 10, tri 8 6 7, tri 9 8 1] }

/-- Torus grid `nu × nv` (outward faces), radii `R`, `r`. -/
def torus (nu nv : Nat) (rr r : Float) : Mesh := Id.run do
  let mut pos : Array V3 := #[]
  for i in [0:nu] do
    for j in [0:nv] do
      let th := twoPi * i.toFloat / nu.toFloat
      let ph := twoPi * j.toFloat / nv.toFloat
      pos := pos.push ⟨(rr + r * Float.cos ph) * Float.cos th, (rr + r * Float.cos ph) * Float.sin th, r * Float.sin ph⟩
  let mut faces : Array Tri := #[]
  for i in [0:nu] do
    for j in [0:nv] do
      let a := i * nv + j
      let b := ((i + 1) % nu) * nv + j
      let c := ((i + 1) % nu) * nv + (j + 1) % nv
      let d := i * nv + (j + 1) % nv
      faces := faces.push (tri a b c) |>.push (tri a c d)
  return { pos, faces }

/-- Open grid patch `nx × ny` quads in the plane `z = bump · sin`, normal `+z`. -/
def grid (nx ny : Nat) (bump : Float) : Mesh := Id.run do
  let mut pos : Array V3 := #[]
  for j in [0:ny + 1] do
    for i in [0:nx + 1] do
      let x := i.toFloat
      let y := j.toFloat
      pos := pos.push ⟨x, y, bump * Float.sin (0.9 * x) * Float.cos (0.7 * y)⟩
  let mut faces : Array Tri := #[]
  for j in [0:ny] do
    for i in [0:nx] do
      let a := j * (nx + 1) + i
      faces := faces.push (tri a (a + 1) (a + nx + 2)) |>.push (tri a (a + nx + 2) (a + nx + 1))
  return { pos, faces }

/-- Möbius band with `n` segments (`2n` faces), non-orientable. -/
def mobius (n : Nat) : Mesh := Id.run do
  let mut pos : Array V3 := #[]
  for i in [0:n] do
    let t := twoPi * i.toFloat / n.toFloat
    for w in [1.0, -1.0] do
      let r := 1.5 + 0.5 * w * Float.cos (t / 2.0)
      pos := pos.push ⟨r * Float.cos t, r * Float.sin t, 0.5 * w * Float.sin (t / 2.0)⟩
  -- vertex 2i = upper u_i, 2i+1 = lower l_i
  let mut faces : Array Tri := #[]
  for i in [0:n] do
    let u0 := 2 * i
    let l0 := 2 * i + 1
    let wrap := i + 1 == n
    let u1 := if wrap then 1 else 2 * (i + 1)
    let l1 := if wrap then 0 else 2 * (i + 1) + 1
    faces := faces.push (tri u0 l0 l1) |>.push (tri u0 l1 u1)
  return { pos, faces }

/-- The union of two meshes, the second moved by `off`. -/
def disjoint (m1 m2 : Mesh) (off : V3) : Mesh :=
  { pos := m1.pos ++ m2.pos.map (· + off)
    faces := m1.faces ++ m2.faces.map fun f => tri (f.a + m1.nv) (f.b + m1.nv) (f.c + m1.nv) }

/-- Scale by `s`. -/
def Mesh.scale (m : Mesh) (s : Float) : Mesh := { m with pos := m.pos.map fun p => s * p }

/-- Independent uniform noise of amplitude `amp` on every coordinate. -/
def Mesh.jitter (m : Mesh) (amp : Float) (seed : Nat) : Mesh := Id.run do
  let mut rng : Rng := ⟨seed.toUInt64 * 0x9E3779B97F4A7C15 + 7⟩
  let mut pos : Array V3 := #[]
  for p in m.pos do
    let (dx, r1) := rng.range (-amp) amp
    let (dy, r2) := r1.range (-amp) amp
    let (dz, r3) := r2.range (-amp) amp
    rng := r3
    pos := pos.push (p + ⟨dx, dy, dz⟩)
  return { m with pos }

/-- Does the mesh, as given, have a consistently oriented component structure? (for the generator) -/
def Mesh.consistent (m : Mesh) : Bool :=
  match m.validate with
  | .ok t => (Array.range t.ncomp).all fun c => (inconsistentEdge t c).isNone
  | .error _ => false

/-- A random mesh of one of the kinds: `(mesh, kind name)`. -/
def randomMesh (kind : Nat) : G (Mesh × String) := do
  let seed ← below 1000000
  match kind with
  | 0 => return (tetra.jitter 0.12 seed, "tetra")
  | 1 => return (octa.jitter 0.1 seed, "octa")
  | 2 => return (icosa.jitter 0.12 seed, "icosa")
  | 3 =>
    let (nu, nv) := #[(5, 3), (6, 3), (6, 4), (5, 4)][(← below 4)]!
    return ((torus nu nv 2.0 0.8).jitter 0.05 seed, "torus")
  | 4 =>
    let (nx, ny) := #[(2, 2), (3, 2), (4, 3), (5, 3)][(← below 4)]!
    return ((grid nx ny 0.3).jitter 0.05 seed, "patch")
  | 5 => let n ← below 8; return ((mobius (6 + n)).jitter 0.03 seed, "mobius")
  | 6 => return ((disjoint (Mesh.scale tetra 0.6) octa ⟨4.0, 0.3, 0.2⟩).jitter 0.05 seed, "two-components")
  | _ =>
    let n ← below 2
    return ((disjoint (octa.jitter 0.05 (seed + 1)) (if n == 0 then octa else Mesh.scale tetra 0.9) ⟨3.5, 0.0, 0.4⟩).jitter 0.03 seed, "two-components-b")

private def exChoices : Array Ex := #[⟨6.0, 12.0⟩, ⟨4.0, 8.0⟩, ⟨3.0, 6.0⟩, ⟨2.0, 5.0⟩, ⟨2.5, 6.0⟩, ⟨4.0, 9.0⟩, ⟨3.0, 7.5⟩]

/-- Random positional constraints that fix every component, plus optional nonlinear ones. -/
private def randomCons (m : Mesh) (nonlinear : Bool) : G (Array Cons) := do
  let some t := (m.validate).toOption | return #[]
  let mut cons : Array Cons := #[]
  let how ← below 3
  if how == 0 then cons := cons.push (.barycenter none none)
  else if how == 1 then
    for c in [0:t.ncomp] do cons := cons.push (.barycenter (some c) none)
  else
    -- one pin per component, a random vertex of it, target moved a little; plus a barycenter on component 0
    for c in [0:t.ncomp] do
      let vs := (Array.range t.nv).filter (t.vertComp[·]! == c)
      let v := vs[(← below vs.size)]!
      let dx ← uni (-0.05) 0.05
      let dy ← uni (-0.05) 0.05
      let dz ← uni (-0.05) 0.05
      cons := cons.push (.pins #[v] (some #[m.pos[v]! + ⟨dx, dy, dz⟩]))
  if nonlinear then
    let k ← below 4
    let g := faceGeo t.faces m.pos
    let a0 := totalArea t g none
    let d1 ← uni (-0.04) 0.04
    let d2 ← uni (-0.04) 0.04
    if k == 1 || k == 3 then cons := cons.push (.area none (some (a0 * (1.0 + d1))))
    if k == 2 || k == 3 then
      let v0 := totalVolume t m.pos none
      if m.consistent && v0.abs > 1e-3 then cons := cons.push (.volume none (some (v0 * (1.0 + d2))))
  return cons

/-- Case `i` of seed `s`. -/
def genCase (seed i : Nat) : RSCase :=
  let rng : Rng := ⟨(seed.toUInt64 * 0x100000001B3) ^^^ (i.toUInt64 * 0x9E3779B97F4A7C15 + 4242)⟩
  let act : G RSCase := do
    let kind ← below 8
    let (mesh0, _) ← randomMesh kind
    let ex := exChoices[(← below exChoices.size)]!
    let tk ← below 16
    let task := if tk < 3 then "energy" else if tk < 5 then "df" else if tk < 7 then "metric" else if tk < 9 then "constraints"
                else if tk < 12 then "saddle" else if tk < 13 then "precond" else if tk < 15 then "step" else "remesh-ops"
    -- volume constraints need orientation: only on consistent meshes (the generator checks)
    let nonlinear := task == "saddle" || task == "step" || task == "constraints"
    let cons ← if task == "saddle" || task == "step" || task == "precond" || task == "constraints"
      then randomCons mesh0 (nonlinear && task != "precond") else pure #[]
    let name := s!"rs-{seed}-{i}"
    let nv := mesh0.nv
    let mut u : Array (Array Float) := #[]
    for _ in [0:2] do
      let mut a : Array Float := #[]
      for _ in [0:nv] do a := a.push (← uni (-1.0) 1.0)
      u := u.push a
    let mut r : Array V3 := #[]
    for _ in [0:nv] do r := r.push ⟨← uni (-1.0) 1.0, ← uni (-1.0) 1.0, ← uni (-1.0) 1.0⟩
    let tau ← uni 0.2 1.4
    let wtp := #[1.0, 0.5, 2.0][(← below 3)]!
    let clo := #[1.0, 0.0][(← below 2)]!
    let ww := #[0.0, 0.0, 0.3][(← below 3)]!
    let armijo := #[1e-4, 1e-2, 0.3][(← below 3)]!
    let startK ← below 3
    let lengths ← uni 0.02 0.3
    let start := if startK == 0 then Start.maxdisp lengths else if startK == 1 then Start.reference else Start.abs (lengths * 0.05)
    let fixedN ← below 2
    -- remesh op on a random edge
    let (op, delaunay, foldover, pinned) ← do
      match mesh0.validate.toOption with
      | none => pure (none, false, false, #[])
      | some t =>
        let e := t.edges[(← below t.edges.size)]!
        let kinds := #["split", "collapse", "flip"]
        let k := kinds[(← below 3)]!
        let d ← below 2
        let fo ← below 2
        let pv ← below 3
        let pinned := if pv == 0 then #[e.1] else #[]
        pure (some (k, e.1, e.2.1), d == 0, fo == 0, pinned)
    return { name, task, ex, mesh := mesh0, u, tau, mo := { wtp, clo, ww }, allMatrices := true, epsRel := 1e-10,
             cons, r, so := { armijo, start, maxHalv := 60, newton := 5, newtonFixed := fixedN == 0 }, op,
             pinned, delaunay, foldover }
  (act.run rng).1

end Ddg.RepulsiveSurfaces
