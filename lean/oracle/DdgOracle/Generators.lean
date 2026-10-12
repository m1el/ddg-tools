import DdgOracle.Operators
import DdgOracle.Hex

/-!
# Deterministic test meshes

Small meshes of every topology the oracle supports, and a reproducible random number generator.
Used by `selftest`, the fixture writer and the throughput benchmark.
-/

namespace Ddg

/-- SplitMix64. -/
structure Rng where
  /-- generator state -/
  state : UInt64

namespace Rng

/-- Next 64 random bits. -/
def next (r : Rng) : UInt64 × Rng :=
  let s := r.state + 0x9E3779B97F4A7C15
  let z := (s ^^^ (s >>> 30)) * 0xBF58476D1CE4E5B9
  let z := (z ^^^ (z >>> 27)) * 0x94D049BB133111EB
  (z ^^^ (z >>> 31), ⟨s⟩)

/-- Uniform in `[0, 1)`. -/
def float (r : Rng) : Float × Rng :=
  let (x, r) := r.next
  (Float.ofNat (x >>> 11).toNat / 9007199254740992.0, r)

/-- Uniform in `[a, b)`. -/
def range (r : Rng) (a b : Float) : Float × Rng :=
  let (x, r) := r.float
  (a + (b - a) * x, r)

/-- Uniform in `{0, …, n-1}`. -/
def below (r : Rng) (n : Nat) : Nat × Rng :=
  let (x, r) := r.next
  (x.toNat % n, r)

end Rng

/-- Vertex positions and triangles. -/
structure RawMesh where
  /-- positions -/
  pos : Array (Point Float)
  /-- triangles -/
  faces : Array Face
  deriving Inhabited

namespace RawMesh

/-- Add independent uniform noise of amplitude `amp` to every coordinate. -/
def jitter (m : RawMesh) (amp : Float) (seed : Nat) : RawMesh := Id.run do
  let mut rng : Rng := ⟨seed.toUInt64⟩
  let mut pos : Array (Point Float) := Array.mkEmpty m.pos.size
  for p in m.pos do
    let (dx, r1) := rng.range (-amp) amp
    let (dy, r2) := r1.range (-amp) amp
    let (dz, r3) := r2.range (-amp) amp
    rng := r3
    pos := pos.push ⟨p.x + dx, p.y + dy, p.z + dz⟩
  return { m with pos }

/-- Apply a map to every vertex position. -/
def mapPos (m : RawMesh) (f : Point Float → Point Float) : RawMesh :=
  { m with pos := m.pos.map f }

/-- Relabel the vertices: old vertex `v` becomes `perm[v]`. -/
def relabel (m : RawMesh) (perm : Array Nat) : RawMesh := Id.run do
  let mut pos : Array (Point Float) := Array.replicate m.pos.size default
  for v in [0:m.pos.size] do
    pos := pos.set! perm[v]! m.pos[v]!
  return { pos, faces := m.faces.map fun f => ⟨perm[f.a]!, perm[f.b]!, perm[f.c]!⟩ }

end RawMesh

/-- A random permutation of `0, …, n-1` (Fisher–Yates). -/
def randomPerm (n seed : Nat) : Array Nat := Id.run do
  let mut rng : Rng := ⟨seed.toUInt64⟩
  let mut a := Array.range n
  for ii in [0:n] do
    let i := n - 1 - ii
    let (j, r) := rng.below (i + 1)
    rng := r
    let x := a[i]!
    a := (a.set! i a[j]!).set! j x
  return a

/-- Rectangular grid `nx × ny` over `[0,1]²` lifted by `height`, a disk with
`(nx-2)(ny-2)` interior vertices. -/
def gridDisk (nx ny : Nat) (height : Float → Float → Float) : RawMesh := Id.run do
  let mut pos : Array (Point Float) := #[]
  let mut faces : Array Face := #[]
  for j in [0:ny] do
    for i in [0:nx] do
      let x := Float.ofNat i / Float.ofNat (nx - 1)
      let y := Float.ofNat j / Float.ofNat (ny - 1)
      pos := pos.push ⟨x, y, height x y⟩
  for j in [0:ny - 1] do
    for i in [0:nx - 1] do
      let a := j * nx + i
      let b := a + 1
      let c := a + nx
      let d := c + 1
      if (i + j) % 2 == 0 then
        faces := (faces.push ⟨a, b, d⟩).push ⟨a, d, c⟩
      else
        faces := (faces.push ⟨a, b, c⟩).push ⟨b, d, c⟩
  return { pos, faces }

/-- Rings of `segs` vertices between the radii `rIn` and `rOut`, lifted by `height`.
With `rIn = 0` the inner ring collapses to a single centre vertex (a disk); otherwise the mesh
is an annulus. -/
def polar (rings segs : Nat) (rIn rOut : Float) (height : Float → Float → Float) : RawMesh :=
  Id.run do
  let disk := rIn == 0.0
  let mut pos : Array (Point Float) := #[]
  let mut faces : Array Face := #[]
  if disk then pos := pos.push ⟨0.0, 0.0, height 0.0 0.0⟩
  let first := if disk then 1 else 0
  let offset := if disk then 1 else 0
  for r in [first:rings + 1] do
    let rad := rIn + (rOut - rIn) * Float.ofNat r / Float.ofNat rings
    for k in [0:segs] do
      -- stagger alternate rings by half a segment
      let ang := 2 * floatPi * (Float.ofNat k + (if r % 2 == 1 then 0.5 else 0.0)) / Float.ofNat segs
      let x := rad * Float.cos ang
      let y := rad * Float.sin ang
      pos := pos.push ⟨x, y, height x y⟩
  let ring := fun (r k : Nat) => offset + (r - first) * segs + k % segs
  if disk then
    for k in [0:segs] do
      faces := faces.push ⟨0, ring 1 k, ring 1 (k + 1)⟩
  for r in [first:rings] do
    for k in [0:segs] do
      let a := ring r k
      let b := ring r (k + 1)
      let c := ring (r + 1) k
      let d := ring (r + 1) (k + 1)
      if r % 2 == 1 then
        -- the inner ring is the staggered one
        faces := (faces.push ⟨a, c, d⟩).push ⟨a, d, b⟩
      else
        faces := (faces.push ⟨a, c, b⟩).push ⟨b, c, d⟩
  return { pos, faces }

/-- Latitude–longitude sphere with `stacks ≥ 2` bands and `slices ≥ 3` meridians, mapped
through `shape` (unit sphere point ↦ surface point). `2 + (stacks-1)·slices` vertices. -/
def uvSphere (stacks slices : Nat) (shape : Point Float → Point Float) : RawMesh := Id.run do
  let mut pos : Array (Point Float) := #[shape ⟨0.0, 0.0, 1.0⟩]
  let mut faces : Array Face := #[]
  for s in [1:stacks] do
    let th := floatPi * Float.ofNat s / Float.ofNat stacks
    for k in [0:slices] do
      let ph := 2 * floatPi * (Float.ofNat k + (if s % 2 == 1 then 0.5 else 0.0)) / Float.ofNat slices
      pos := pos.push (shape ⟨Float.sin th * Float.cos ph, Float.sin th * Float.sin ph, Float.cos th⟩)
  pos := pos.push (shape ⟨0.0, 0.0, -1.0⟩)
  let south := pos.size - 1
  let ring := fun (s k : Nat) => 1 + (s - 1) * slices + k % slices
  for k in [0:slices] do
    faces := faces.push ⟨0, ring 1 k, ring 1 (k + 1)⟩
    faces := faces.push ⟨south, ring (stacks - 1) (k + 1), ring (stacks - 1) k⟩
  for s in [1:stacks - 1] do
    for k in [0:slices] do
      let a := ring s k
      let b := ring s (k + 1)
      let c := ring (s + 1) k
      let d := ring (s + 1) (k + 1)
      if s % 2 == 1 then
        faces := (faces.push ⟨a, d, b⟩).push ⟨a, c, d⟩
      else
        faces := (faces.push ⟨a, c, b⟩).push ⟨b, c, d⟩
  return { pos, faces }

/-- Torus with `nu × nv` vertices, major radius `R`, minor radius `r`, with the tube radius
modulated by `bump` (a function of the two angles). -/
def torus (nu nv : Nat) (R r : Float) (bump : Float → Float → Float) : RawMesh := Id.run do
  let mut pos : Array (Point Float) := #[]
  let mut faces : Array Face := #[]
  for i in [0:nu] do
    let a := 2 * floatPi * Float.ofNat i / Float.ofNat nu
    for j in [0:nv] do
      let b := 2 * floatPi * Float.ofNat j / Float.ofNat nv
      let rr := r * (1 + bump a b)
      let d := R + rr * Float.cos b
      pos := pos.push ⟨d * Float.cos a, d * Float.sin a, rr * Float.sin b⟩
  let at' := fun (i j : Nat) => (i % nu) * nv + j % nv
  for i in [0:nu] do
    for j in [0:nv] do
      let a := at' i j
      let b := at' (i + 1) j
      let c := at' i (j + 1)
      let d := at' (i + 1) (j + 1)
      faces := (faces.push ⟨a, b, d⟩).push ⟨a, d, c⟩
  return { pos, faces }

/-- Regular tetrahedron. -/
def tetrahedron : RawMesh :=
  { pos := #[⟨1, 1, 1⟩, ⟨1, -1, -1⟩, ⟨-1, 1, -1⟩, ⟨-1, -1, 1⟩]
    faces := #[⟨0, 1, 2⟩, ⟨0, 3, 1⟩, ⟨0, 2, 3⟩, ⟨1, 3, 2⟩] }

/-- Bipyramid over a regular `k`-gon with apex heights `±h`; `k = 4`, `h = 1` is the regular
octahedron. `k + 2` vertices. -/
def bipyramid (k : Nat) (h : Float) : RawMesh := Id.run do
  let mut pos : Array (Point Float) := #[⟨0.0, 0.0, h⟩, ⟨0.0, 0.0, -h⟩]
  let mut faces : Array Face := #[]
  for i in [0:k] do
    let a := 2 * floatPi * Float.ofNat i / Float.ofNat k
    pos := pos.push ⟨Float.cos a, Float.sin a, 0.0⟩
  for i in [0:k] do
    let p := 2 + i
    let q := 2 + (i + 1) % k
    faces := (faces.push ⟨0, p, q⟩).push ⟨1, q, p⟩
  return { pos, faces }

/-- The cube `[-1,1]³` with each square split along a diagonal. -/
def cube : RawMesh :=
  { pos := #[⟨-1, -1, -1⟩, ⟨1, -1, -1⟩, ⟨1, 1, -1⟩, ⟨-1, 1, -1⟩,
             ⟨-1, -1, 1⟩, ⟨1, -1, 1⟩, ⟨1, 1, 1⟩, ⟨-1, 1, 1⟩]
    faces := #[⟨0, 2, 1⟩, ⟨0, 3, 2⟩, ⟨4, 5, 6⟩, ⟨4, 6, 7⟩, ⟨0, 1, 5⟩, ⟨0, 5, 4⟩,
               ⟨1, 2, 6⟩, ⟨1, 6, 5⟩, ⟨2, 3, 7⟩, ⟨2, 7, 6⟩, ⟨3, 0, 4⟩, ⟨3, 4, 7⟩] }

/-- A fan of `k` triangles around one interior vertex (index 0) raised to height `h`. -/
def fan (k : Nat) (h : Float) : RawMesh := Id.run do
  let mut pos : Array (Point Float) := #[⟨0.0, 0.0, h⟩]
  let mut faces : Array Face := #[]
  for i in [0:k] do
    let a := 2 * floatPi * Float.ofNat i / Float.ofNat k
    pos := pos.push ⟨Float.cos a, Float.sin a, 0.0⟩
  for i in [0:k] do
    faces := faces.push ⟨0, 1 + i, 1 + (i + 1) % k⟩
  return { pos, faces }

/-- Topology and operators of a raw mesh. -/
def RawMesh.operators (m : RawMesh) : Except String Operators := do
  let topo ← buildTopology m.pos.size m.faces
  buildOperators topo m.pos

end Ddg
