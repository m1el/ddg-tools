import DdgOracle.Cuts.Gen
import DdgOracle.Developability.Format

/-!
# Seeded case generator for the developability oracle

Stars of valence 3-10 (interior and boundary fans), near-hinge and near-flat stars, regular cones,
and small closed and open meshes (jittered icospheres, cube, octahedron, tetrahedron, a bumpy torus,
crumpled grids, annuli). Case `i` of seed `s` depends on `(s, i)` only.
-/

namespace Ddg.Developability

open Ddg

/-- Convert a raw mesh of the other tasks. -/
def ofRaw (m : RawMesh) : DMesh :=
  { pos := m.pos.map fun p => ⟨p.x, p.y, p.z⟩, faces := m.faces.map fun f => (f.a, f.b, f.c) }

private abbrev G := StateM Rng

private def uni (a b : Float) : G Float := modifyGet fun r => r.range a b
private def below (n : Nat) : G Nat := modifyGet fun r => r.below n

private def twoPi : Float := 6.283185307179586

/-- Interior star: apex `0` at height `h`, a ring of `k` vertices with random angles, radii and
heights of amplitude `amp`. -/
def interiorStar (k : Nat) (h amp : Float) : G DMesh := do
  let mut pos : Array V3 := #[⟨0, 0, h⟩]
  for m in [0:k] do
    let jit ← uni (-0.35) 0.35
    let a := twoPi * (m.toFloat + jit) / k.toFloat
    let r ← uni 0.6 1.4
    let z ← uni (-amp) amp
    pos := pos.push ⟨r * Float.cos a, r * Float.sin a, z⟩
  let faces := (Array.range k).map fun m => (0, 1 + m, 1 + (m + 1) % k)
  return { pos, faces }

/-- Boundary fan: `k` faces over an angular span `span`. -/
def boundaryFan (k : Nat) (span h amp : Float) : G DMesh := do
  let mut pos : Array V3 := #[⟨0, 0, h⟩]
  for m in [0:k + 1] do
    let jit ← uni (-0.3) 0.3
    let a := span * (m.toFloat + jit * (if m == 0 || m == k then 0 else 1)) / k.toFloat
    let r ← uni 0.6 1.4
    let z ← uni (-amp) amp
    pos := pos.push ⟨r * Float.cos a, r * Float.sin a, z⟩
  let faces := (Array.range k).map fun m => (0, 1 + m, 2 + m)
  return { pos, faces }

/-- Near-hinge star: ring points on the planes `z = ± a x`, with two ring points on the crease,
plus noise of amplitude `eps`. `kr + kl + 2` faces. -/
def hingeStar (kr kl : Nat) (a eps : Float) : G DMesh := do
  let mut pos : Array V3 := #[⟨0, 0, 0⟩]
  let ring : Array Float := Id.run do
    let mut angs : Array Float := #[]
    angs := angs.push (twoPi / 4)
    for i in [0:kl] do angs := angs.push (twoPi / 4 + twoPi / 2 * (i + 1).toFloat / (kl + 1).toFloat)
    angs := angs.push (3 * twoPi / 4)
    for i in [0:kr] do angs := angs.push (3 * twoPi / 4 + twoPi / 2 * (i + 1).toFloat / (kr + 1).toFloat)
    return angs
  for ang in ring do
    let r ← uni 0.6 1.4
    let x := r * Float.cos ang
    let y := r * Float.sin ang
    let z := a * x.abs
    pos := pos.push ⟨x, y, z⟩
  let mut noisy : Array V3 := #[]
  for p in pos do
    let dx ← uni (-eps) eps
    let dy ← uni (-eps) eps
    let dz ← uni (-eps) eps
    noisy := noisy.push ⟨p.x + dx, p.y + dy, p.z + dz⟩
  let k := ring.size
  let faces := (Array.range k).map fun m => (0, 1 + m, 1 + (m + 1) % k)
  return { pos := noisy, faces }

/-- Regular cone of valence `k` and height `h` (SPEC 10.3). -/
def regularCone (k : Nat) (h : Float) : DMesh :=
  ofRaw (fan k h)

private def pick {α : Type} [Inhabited α] (xs : Array α) : G α := do
  let i ← below xs.size
  return xs[i]!

/-- Names of the generator kinds, indexed by `i % kinds.size`. -/
def kinds : Array String :=
  #["star", "fan", "hinge", "flat", "cone", "icosphere", "polyhedron", "torus", "crumpled",
    "annulus", "folded"]

/-- The mesh of kind `kind`. -/
def genMesh (kind : Nat) : G (DMesh × String) := do
  match kind % kinds.size with
  | 0 =>
    let k := 3 + (← below 8)
    let h ← uni 0.1 1.5
    let amp ← uni 0.02 0.8
    return (← interiorStar k h amp, s!"star{k}")
  | 1 =>
    let k := 1 + (← below 9)
    let span ← uni 1.0 5.5
    let h ← uni 0.1 1.5
    let amp ← uni 0.02 0.8
    return (← boundaryFan k span h amp, s!"fan{k}")
  | 2 =>
    let kr := 1 + (← below 4)
    let kl := 1 + (← below 4)
    let a ← uni 0.3 2.0
    let eps ← pick #[0.0, 1e-9, 1e-6, 1e-3, 1e-2]
    return (← hingeStar kr kl a eps, s!"hinge{kr}{kl}")
  | 3 =>
    let k := 3 + (← below 8)
    let eps ← pick #[1e-9, 1e-6, 1e-4, 1e-3]
    return (← interiorStar k 0.0 eps, s!"flat{k}")
  | 4 =>
    let k := 3 + (← below 8)
    let h ← uni 0.05 3.0
    return (regularCone k h, s!"cone{k}")
  | 5 =>
    let lvl := (← below 2)
    let amp ← uni 0.01 0.15
    let seed ← below 1000000
    let m := ((Ddg.Cuts.icosphere lvl).jitter amp seed)
    return (ofRaw m, s!"icosphere{lvl}")
  | 6 =>
    let which ← below 3
    let amp ← uni 0.005 0.3
    let seed ← below 1000000
    let base := if which == 0 then cube else if which == 1 then tetrahedron else bipyramid 4 1.0
    return (ofRaw (base.jitter amp seed), s!"polyhedron{which}")
  | 7 =>
    let nu := 4 + (← below 3)
    let nv := 3 + (← below 3)
    let amp ← uni 0.0 0.1
    let seed ← below 1000000
    let m := (torus nu nv 2.0 0.7 fun a b => 0.1 * Float.sin (a + 2 * b)).jitter amp seed
    return (ofRaw m, s!"torus{nu}x{nv}")
  | 8 =>
    let n := 4 + (← below 5)
    let amp ← uni 0.05 0.5
    let seed ← below 1000000
    let m := (gridDisk n n fun x y => 0.3 * Float.sin (5 * x) * y).jitter amp seed
    return (ofRaw m, s!"crumpled{n}")
  | 9 =>
    let rings := 2 + (← below 3)
    let segs := 5 + (← below 6)
    let amp ← uni 0.0 0.2
    let seed ← below 1000000
    let m := (polar rings segs 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.2 * y * y).jitter amp seed
    return (ofRaw m, s!"annulus{rings}x{segs}")
  | _ =>
    let k := 3 + (← below 8)
    let h ← uni (-0.5) 1.0
    return (← interiorStar k h 1.8, s!"folded{k}")

/-- Case `i` of seed `seed`. -/
def genCase (seed i : Nat) : DevCase :=
  let rng : Rng := ⟨(seed * 1000003 + i * 7919 + 12345).toUInt64⟩
  let (c, _) := (do
    let (mesh, tag) ← genMesh i
    let r ← uni 0 1
    let taskName := if r < 0.5 then "energy" else if r < 0.9 then "step" else "subdivide"
    let incl ← pick #[true, true, true, false]
    let omitV ← pick #[true, true, true, false]
    let e ← pick EnergyKind.all
    let init ← pick #[0.1, 0.1, 0.01, 1.0, 3.0]
    let maxd ← pick #[0.5, 0.5, 2.0, 0.2]
    let c1 ← pick #[0.0, 0.0, 1e-4]
    let c2 ← pick #[0.5, 0.5, 0.9, 0.1]
    let np ← below 4
    let mut pins : Array Nat := #[]
    if np == 3 then
      for _ in [0:2] do pins := pins.push (← below mesh.pos.size)
    return ({ name := s!"{tag}-s{seed}-{i}", task := taskName,
              energies := if taskName == "energy" then EnergyKind.all else #[e],
              includeBoundary := incl, omit3 := omitV,
              stepOpts := { initialStep := init, maxDisplacement := maxd, c1, c2, pinned := pins },
              mesh } : DevCase) : G DevCase).run rng
  c

end Ddg.Developability
