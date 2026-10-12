import DdgOracle.Generators
import DdgOracle.Certificate

/-!
# Test cases: a mesh together with parameters

Deterministic builders for the parameter features of the canonical problem, used by `selftest`,
the fixtures and the benchmark.
-/

namespace Ddg

/-- Uniform parameters: `wE = wR = 1`, no bounds, no prescribed cones. -/
def uniformParams (nv : Nat) (lambda : Float) : Params :=
  { lambda
    wE := FloatArray.ofFn' nv fun _ => 1.0
    wR := FloatArray.ofFn' nv fun _ => 1.0
    lo := FloatArray.ofFn' nv fun _ => negInf
    hi := FloatArray.ofFn' nv fun _ => posInf
    cones := #[] }

namespace Params

/-- Set `wR = 0` at the given vertices (free cone sites). -/
def withFree (p : Params) (vs : Array Nat) : Params :=
  { p with wR := vs.foldl (fun a v => a.set! v 0.0) p.wR }

/-- Set `wR = +∞` at the vertices selected by `sel` (forbidden region). -/
def withForbidden (p : Params) (sel : Nat → Bool) : Params :=
  { p with wR := FloatArray.ofFn' p.wR.size fun v => if sel v then posInf else p.wR.get! v }

/-- Set `wE = 0` at the vertices selected by `sel`. -/
def withZeroDistortion (p : Params) (sel : Nat → Bool) : Params :=
  { p with wE := FloatArray.ofFn' p.wE.size fun v => if sel v then 0.0 else p.wE.get! v }

/-- The same bounds `[lo, hi]` at every vertex. -/
def withBounds (p : Params) (lo hi : Float) : Params :=
  { p with
    lo := FloatArray.ofFn' p.lo.size fun _ => lo
    hi := FloatArray.ofFn' p.hi.size fun _ => hi }

/-- Prescribe cones. -/
def withCones (p : Params) (cones : Array (Nat × Float)) : Params :=
  { p with cones := p.cones ++ cones }

/-- Independent random weights, log-uniform in `[1/r, r]`. -/
def withRandomWeights (p : Params) (r : Float) (seed : Nat) : Params := Id.run do
  let mut rng : Rng := ⟨seed.toUInt64⟩
  let mut wE := FloatArray.emptyWithCapacity p.wE.size
  let mut wR := FloatArray.emptyWithCapacity p.wR.size
  let l := Float.log r
  for v in [0:p.wE.size] do
    let (a, r1) := rng.range (-l) l
    let (b, r2) := r1.range (-l) l
    rng := r2
    wE := wE.push (p.wE.get! v * Float.exp a)
    wR := wR.push (p.wR.get! v * Float.exp b)
  return { p with wE, wR }

end Params

/-- A named test case. -/
structure Case where
  /-- identifier (no whitespace) -/
  name : String
  /-- the mesh -/
  mesh : RawMesh
  /-- the parameters -/
  params : Params
  deriving Inhabited

/-- The interior vertices of a mesh, in increasing order (empty if the mesh is invalid). -/
def interiorVertices (m : RawMesh) : Array Nat :=
  match buildTopology m.pos.size m.faces with
  | .ok t => (Array.range t.nv).filter fun v => !t.isBoundary[v]!
  | .error _ => #[]

/-- A random mesh of one of the supported topologies with roughly `size` vertices. -/
def randomMesh (size seed : Nat) : String × RawMesh := Id.run do
  let rng : Rng := ⟨(seed * 7919 + 13).toUInt64⟩
  let (kind, rng) := rng.below 6
  let (a, rng) := rng.range 0.1 0.6
  let (b, rng) := rng.range 0.2 0.8
  let (c, rng) := rng.range 0.2 0.8
  let (amp, _) := rng.range 0.0 0.04
  let side := (Float.sqrt (Float.ofNat size)).toUInt64.toNat
  let side := if side < 3 then 3 else side
  match kind with
  | 0 =>
    let m := gridDisk (side + 1) side fun x y => a * Float.exp (-9 * ((x - b) * (x - b) + (y - c) * (y - c)))
    ("grid", m.jitter (amp / Float.ofNat side) seed)
  | 1 =>
    let segs := if side < 5 then 5 else side + 2
    let rings := if size / segs < 2 then 2 else size / segs
    let m := polar rings segs 0.0 1.0 fun x y => a * (1 - x * x - y * y) + 0.3 * b * x * y
    ("disk", m.jitter (amp / Float.ofNat rings) seed)
  | 2 =>
    let segs := if side < 6 then 6 else side + 3
    let rings := if size / segs < 2 then 2 else size / segs - 1
    let m := polar rings segs 0.35 1.0 fun x y => a * Float.sin (3 * x + b) * Float.cos (2 * y + c)
    ("annulus", m.jitter (amp / Float.ofNat (rings + 2)) seed)
  | 3 =>
    let slices := if side < 4 then 4 else side
    let stacks := if size / slices < 3 then 3 else size / slices
    let m := uvSphere stacks slices fun p => ⟨p.x * (1 + a), p.y * (0.6 + b), p.z * (0.6 + c)⟩
    ("sphere", m.jitter (amp / Float.ofNat stacks) seed)
  | 4 =>
    let nu := if side < 4 then 4 else side + 1
    let nv := if size / nu < 3 then 3 else size / nu
    let m := torus nu nv 1.0 (0.25 + 0.3 * a) fun s t => b * 0.3 * Float.sin (2 * s) * Float.cos t
    ("torus", m.jitter (amp / Float.ofNat nu) seed)
  | _ =>
    let k := if size < 5 then 3 else size - 2
    ("bipyramid", (bipyramid k (0.4 + a)).jitter (amp / Float.ofNat k) seed)

/-- Random parameters exercising every feature of the canonical problem. -/
def randomParams (m : RawMesh) (seed : Nat) : Params := Id.run do
  let nv := m.pos.size
  let interior := interiorVertices m
  let rng : Rng := ⟨(seed * 104729 + 71).toUInt64⟩
  let (ll, rng) := rng.range (-3.0) 0.5
  let mut p := uniformParams nv (Float.exp (ll * Float.log 10.0))
  let (weighted, rng) := rng.below 3
  if weighted == 0 then p := p.withRandomWeights 3.0 (seed + 1)
  let (bk, rng) := rng.below 6
  let mut rng := rng
  p := match bk with
    | 0 => p.withBounds 0.0 posInf
    | 1 => p.withBounds (-(floatPi / 2)) (floatPi / 2)
    | 2 => p.withBounds 0.0 (floatPi / 2)
    | 3 => p.withBounds negInf 0.0
    | _ => p
  if interior.size > 0 then
    let (nfree, r) := rng.below 4
    rng := r
    if nfree == 0 then
      let (i, r) := rng.below interior.size
      rng := r
      p := p.withFree #[interior[i]!]
    let (forb, r) := rng.below 4
    rng := r
    if forb == 0 then
      let (phase, r) := rng.below 3
      rng := r
      p := p.withForbidden fun v => v % 3 == phase
    let (zeroE, r) := rng.below 6
    rng := r
    if zeroE == 0 then
      let (phase, r) := rng.below 4
      rng := r
      p := p.withZeroDistortion fun v => v % 4 == phase
    let (pres, r) := rng.below 4
    rng := r
    if pres == 0 then
      let (i, r) := rng.below interior.size
      let (ang, r) := r.range (-1.0) 1.5
      rng := r
      p := p.withCones #[(interior[i]!, ang)]
  return p

end Ddg
