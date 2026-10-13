import DdgOracle.Fixtures
import DdgOracle.Cuts.Selftest

/-!
# Golden fixtures of the cut flow

`cuts-fixtures write DIR` writes, for every case `NAME`, `NAME.in` (the case in the input format)
and `NAME.out` (the output of `ddg-oracle cuts NAME.in`), and `MANIFEST.tsv`. `cuts-fixtures
check DIR` recomputes them and compares with the files (floats to `10⁻⁹`).
-/

namespace Ddg.Cuts

/-- A fixture: a case and what it is for. -/
structure Fixture where
  /-- the case; its name is the file stem -/
  case : CutCase
  /-- one-line description for the manifest -/
  what : String

def polarAngleFrom' (axis : Point Float) (p : Point Float) : Float :=
  let n := Float.sqrt (axis.x * axis.x + axis.y * axis.y + axis.z * axis.z)
  Float.acos ((axis.x * p.x + axis.y * p.y + axis.z * p.z) / n)

/-- Class label of every vertex of `m` under the reflection `x ↦ -x`: the smaller index of the
vertex and its mirror image. -/
def mirrorClasses (m : RawMesh) : Array Nat := Id.run do
  let n := m.pos.size
  let mut out : Array Nat := Array.replicate n 0
  for i in [0:n] do
    let p := m.pos[i]!
    let mut best := i
    let mut bd := posInf
    for j in [0:n] do
      let q := m.pos[j]!
      let d := (q.x + p.x) * (q.x + p.x) + (q.y - p.y) * (q.y - p.y) + (q.z - p.z) * (q.z - p.z)
      if d < bd then
        bd := d
        best := j
    out := out.set! i (min i best)
  return out

/-- The fixture list. -/
def fixtures : Array Fixture := Id.run do
  -- x and y are warped (z is kept, so that "z = 0" stays exact): the symmetric icosphere gives
  -- exact ties between the two diagonals of a cut quad, which would make the fixtures fragile
  let warp := fun (m : RawMesh) => m.mapPos fun p =>
    ⟨p.x * (1 + 0.05 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.04 * Float.cos (4 * p.x - 2 * p.z)), p.z⟩
  let ico1 := warp (icosphere 1)
  let ico2 := warp (icosphere 2)
  let ico2s := icosphere 2
  let cap := fun (m : RawMesh) (psi0 : Float) => phiOf m fun p => Float.acos p.z - psi0
  let capped := StepRule.capped 0.01 0.25
  let n2 := ico2.pos.size
  let n1 := ico1.pos.size
  let torusM := (torus 12 8 1.0 0.42 fun a b => 0.1 * Float.sin (2 * a) * Float.cos b).jitter 0.005 21
  let grid := (gridDisk 12 11 fun x y => 0.3 * Float.sin (5 * x) * Float.cos (4 * y) + 0.2 * x * x).jitter 0.004 22
  let ann := (polar 5 16 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.25 * y * y).jitter 0.005 23
  let flat := (polar 10 28 0.0 1.0 fun _ _ => 0.0).jitter 0.004 24
  let two := RawMesh.union ico1 ico1 ⟨3.0, 0.0, 0.0⟩
  let mut fs : Array Fixture := #[]
  let add := fun (fs : Array Fixture) (c : CutCase) (what : String) => fs.push { case := c, what }
  fs := add fs { mkCase "ico1-tilted-dirichlet" ico1 (tilted ico1 0.3 0.2 0.93 0.1) with step := capped }
    "icosphere level 1 (42 vertices), tilted-plane level set, Dirichlet, capped step"
  fs := add fs { mkCase "ico2-tilted-dirichlet" ico2 (tilted ico2 0.3 0.2 0.93 0.15) with step := capped, flux := .nodal, passes := 1 }
    "icosphere level 2 (162 vertices), tilted-plane level set, Dirichlet, capped step; the paper's procedure (nodal v/l, one smoothing pass)"
  fs := add fs { mkCase "ico2-tilted-dirichlet-recovery" ico2 (tilted ico2 0.3 0.2 0.93 0.15) with step := capped }
    "the same with the default flux recovery (c = 0.5, no sliver smoothing): sane speeds where the nodal quotient blows up"
  fs := add fs { mkCase "ico2-cap-dirichlet" ico2 (cap ico2 1.0) with step := .fixed 0.001 }
    "icosphere level 2, geodesic cap of polar angle 1, Dirichlet, fixed tau"
  fs := add fs { mkCase "ico2-cap-hencky" ico2 (cap ico2 1.0) with energy := .hencky, step := .fixed 0.001 }
    "icosphere level 2, geodesic cap, Hencky energy"
  fs := add fs { mkCase "ico2-cap-scaled-curv" ico2 (cap ico2 1.0) with energy := .scaledCurv 1.0 0.02, step := .fixed 0.001 }
    "icosphere level 2, geodesic cap, curvature-scaled Dirichlet (r_eps = 1, c = 0.02)"
  fs := add fs { mkCase "torus-band-dirichlet" torusM (phiOf torusM fun p => Float.cos (3 * Float.atan2 p.y p.x) - 0.2) with step := capped }
    "torus (96 vertices), two-band level set, Dirichlet; cut components are loops"
  fs := add fs { mkCase "grid-boundary-scaled-a" grid (phiOf grid fun p => p.x + 0.3 * p.y - 0.65) with
      energy := .scaledA (FloatArray.ofFn' grid.pos.size fun i => 1.0 + 0.5 * Float.sin (7.0 * Float.ofNat i)),
      step := capped }
    "grid disk (132 vertices), straight cut that meets the boundary twice, scaled Dirichlet with an explicit a"
  fs := add fs { mkCase "annulus-ring-hencky" ann (phiOf ann fun p => Float.sqrt (p.x * p.x + p.y * p.y) - 0.7) with
      energy := .hencky, step := capped }
    "annulus, circular cut between the two boundary circles, Hencky"
  fs := add fs { mkCase "annulus-radial-dirichlet" ann (phiOf ann fun p => p.x - 0.1) with step := capped }
    "annulus, straight cut across the ring (meets both boundary circles), Dirichlet"
  let etaZ := FloatArray.ofFn' n2 fun i => 0.2 * (1 + ico2.pos[i]!.z)
  fs := add fs { mkCase "ico2-eta-gradient" ico2 (cap ico2 1.0) with eta := some etaZ, step := .fixed 0.001 }
    "geodesic cap with a penalty eta = 0.2 (1 + z), default (gradient) rule (W_V + tau L_w) phi' = W_V phi~"
  fs := add fs { mkCase "ico2-eta-paper" ico2 (cap ico2 1.0) with eta := some etaZ, penalty := .paperDiagonal, step := .fixed 0.001 }
    "the same, the paper's diagonal rule (W_V D^-1 + tau L) phi' = W_V D^-1 phi~"
  fs := add fs { mkCase "ico2-region-weights" ico2 (cap ico2 1.0) with w := (2.0, 0.5), beta := (0.7, 0.2), step := .fixed 0.001 }
    "geodesic cap, w = (2, 0.5), beta = (0.7, 0.2)"
  let tiltedCap := phiOf ico2s fun p => polarAngleFrom' ⟨0.15, 0.05, 1.0⟩ p - 1.0
  fs := add fs { mkCase "ico2-symmetry" ico2s tiltedCap with classes := some (mirrorClasses ico2s), step := .fixed 0.001 }
    "geodesic cap about a tilted axis (not mirror symmetric), classes = mirror x -> -x: both projections are active"
  fs := add fs { mkCase "ico1-phi-zero-vertex" ico1 (phiOf ico1 fun p => p.z) with step := capped }
    "icosphere level 1, phi = z: phi is exactly 0 at 8 vertices, so crossings have t = 0 and are moved to delta; default flux recovery"
  fs := add fs { mkCase "ico1-phi-zero-vertex-nodal" ico1 (phiOf ico1 fun p => p.z) with flux := .nodal, passes := 1, step := capped }
    "the same with the paper's procedure (nodal, one pass): max |s| about 1e6"
  fs := add fs { mkCase "ico1-delta-clamp" ico1 (phiOf ico1 fun p => p.z - 0.02) with delta := 0.2, flux := .nodal, rho0 := 2.0, passes := 2, step := capped }
    "icosphere level 1, delta = 0.2 so several crossings hit the clamp; nodal, rho0 = 2, two smoothing passes"
  fs := add fs { mkCase "ico2-sliver-smoothing" ico2 (tilted ico2 0.3 0.2 0.93 0.15) with flux := .nodal, rho0 := 2.0, passes := 3, step := capped }
    "level 2, tilted plane, nodal, rho0 = 2 with three smoothing passes"
  fs := add fs { mkCase "ico2-recovery-sliver-passes" ico2 (tilted ico2 0.3 0.2 0.93 0.15) with rho0 := 2.0, passes := 2, step := capped }
    "level 2, tilted plane, flux recovery c = 0.5 followed by two sliver-smoothing passes (rho0 = 2)"
  fs := add fs { mkCase "ico2-cap-flux-c1" ico2 (cap ico2 1.0) with flux := .smoothed 1.0, step := .fixed 0.001 }
    "geodesic cap, flux recovery with c = 1"
  fs := add fs { mkCase "ico2-cap-delta-zero" ico2 (cap ico2 1.0) with delta := 0.0, step := .fixed 0.001 }
    "geodesic cap, delta = 0: no crossing is moved"
  fs := add fs { mkCase "flat-disk-uniform-curvature" flat (phiOf flat fun p => Float.sqrt (p.x * p.x + p.y * p.y) - 0.3 * Float.sqrt Ddg.floatPi) with
      curv := .uniform 5.0, step := .fixed 0.001 }
    "flat disk with a curvature override (uniform 5) and a circular cut of scaled radius 0.3: closed forms of SPEC section 1"
  fs := add fs { mkCase "flat-disk-gaussian" flat (phiOf flat fun p => Float.sqrt (p.x * p.x + p.y * p.y) - 0.3 * Float.sqrt Ddg.floatPi) with
      alphaL := 0.1, step := .fixed 0.001 }
    "flat disk, Gaussian curvature (zero): u = 0, s = 0, pure curve shortening"
  let perFace := FloatArray.ofFn' ico1.faces.size fun i => 3.0 + 2.0 * Float.sin (1.7 * Float.ofNat i)
  fs := add fs { mkCase "ico1-curvature-per-face" ico1 (cap ico1 1.0) with curv := .perFace perFace, energy := .hencky, step := .fixed 0.001 }
    "icosphere level 1, per-face curvature override, Hencky"
  let phiTwoA := FloatArray.ofFn' two.pos.size fun i => if i < n1 then Float.acos two.pos[i]!.z - 1.0 else 1.0
  fs := add fs { mkCase "two-spheres-untouched-1" two phiTwoA with step := .fixed 0.001 }
    "two disjoint icospheres; the cut lies on the first, the second is entirely region 1: a closed component the cut does not touch"
  let phiTwoB := FloatArray.ofFn' two.pos.size fun i => if i < n1 then Float.acos two.pos[i]!.z - 1.0 else -1.0
  fs := add fs { mkCase "two-spheres-untouched-0" two phiTwoB with energy := .hencky, beta := (0.5, 0.5), step := .fixed 0.001 }
    "the same with the second sphere entirely region 0, Hencky, beta = (0.5, 0.5)"
  return fs

/-- The two texts of a fixture. -/
def Fixture.render (f : Fixture) : String × String :=
  (printCase f.case, runCase f.case)

/-- Write a list of fixtures to `dir`. -/
def writeList (fs : Array Fixture) (dir : System.FilePath) : IO Unit := do
  IO.FS.createDirAll dir
  let mut manifest := ""
  for f in fs do
    let (input, output) := f.render
    IO.FS.writeFile (dir / s!"{f.case.name}.in") input
    IO.FS.writeFile (dir / s!"{f.case.name}.out") output
    manifest := manifest ++ s!"{f.case.name}\t{statusOf output}\t{f.what}\n"
    IO.println s!"{f.case.name}: {statusOf output}"
  IO.FS.writeFile (dir / "MANIFEST.tsv") manifest
  IO.println s!"wrote {fs.size} fixtures to {dir}"

/-- Recompute a list of fixtures and compare with the files in `dir`; returns the number of
differences. -/
def checkList (fs : Array Fixture) (dir : System.FilePath) : IO Nat := do
  let mut bad := 0
  for f in fs do
    let (input, output) := f.render
    let name := f.case.name
    let mut problems : Array String := #[]
    for (suffix, text, exact) in [(".in", input, true), (".out", output, false)] do
      let path := dir / s!"{name}{suffix}"
      if !(← path.pathExists) then
        problems := problems.push s!"{name}{suffix} is missing"
      else
        let stored ← IO.FS.readFile path
        if exact then
          if stored != text then problems := problems.push s!"{name}{suffix} differs from the generated input"
        else if let some d := compareOutputs stored text then
          problems := problems.push s!"{name}{suffix}: {d}"
    if problems.size == 0 then IO.println s!"ok    {name}"
    else
      bad := bad + 1
      for p in problems do IO.println s!"FAIL  {p}"
  IO.println s!"cuts fixtures: {fs.size - bad} of {fs.size} match"
  return bad

/-- Write the fixtures of the paper's stages. -/
def writeFixtures (dir : System.FilePath) : IO Unit := writeList fixtures dir

/-- Check the fixtures of the paper's stages. -/
def checkFixtures (dir : System.FilePath) : IO Nat := checkList fixtures dir

end Ddg.Cuts
