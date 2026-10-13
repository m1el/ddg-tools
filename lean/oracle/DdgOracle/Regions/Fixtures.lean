import DdgOracle.Cuts.FixturesDefault
import DdgOracle.Regions.Format

/-!
# Golden fixtures of the N-region flow

`fixtures/regions/NAME.in|out`, `MANIFEST.tsv`; written and checked like the others.
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- A fixture. -/
structure RFixture where
  /-- the case -/
  case : RegionCase
  /-- description for the manifest -/
  what : String

/-- The fixtures. -/
def fixtures : Array RFixture := Id.run do
  let warp := fun (m : RawMesh) => m.mapPos fun p =>
    ⟨p.x * (1 + 0.05 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.04 * Float.cos (4 * p.x - 2 * p.z)), p.z⟩
  let ico2 := warp (icosphere 2)
  -- a warp that is even in x: the mesh stays mirror symmetric (x -> -x) but is irregular
  let ico2s := (icosphere 2).mapPos fun p =>
    ⟨p.x * (1 + 0.11 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.09 * Float.cos (4 * p.x * p.x - 2 * p.z)), p.z⟩
  let torusM := (torus 14 8 1.0 0.42 fun a b => 0.1 * Float.sin (2 * a) * Float.cos b).jitter 0.005 21
  let grid := (gridDisk 14 13 fun x y => 0.3 * Float.sin (5 * x) * Float.cos (4 * y) + 0.2 * x * x).jitter 0.004 22
  let pt := fun (x y z : Float) => (⟨x, y, z⟩ : Point Float)
  let nrm := fun (p : Point Float) =>
    let n := Float.sqrt (p.x * p.x + p.y * p.y + p.z * p.z)
    (⟨p.x / n, p.y / n, p.z / n⟩ : Point Float)
  let tetra : Array (Point Float) := #[nrm (pt 1.0 1.13 0.87), nrm (pt 1.07 (-0.93) (-1.1)), nrm (pt (-1.1) 0.97 (-0.9)), nrm (pt (-0.9) (-1.05) 1.13)]
  let cube : Array (Point Float) := #[pt 1.0 0.05 0.0, pt (-1.0) 0.0 0.05, pt 0.05 1.0 0.0, pt 0.0 (-1.0) (-0.05), pt 0.0 0.05 1.0, pt (-0.05) 0.0 (-1.0)]
  let three : Array (Point Float) := #[nrm (pt 1.0 0.1 0.2), nrm (pt (-0.6) 0.9 (-0.1)), nrm (pt (-0.5) (-0.8) 0.3)]
  let mut fs : Array RFixture := #[]
  let add := fun (fs : Array RFixture) (c : RegionCase) (what : String) => fs.push { case := c, what }
  let dflt := fun (c : RegionCase) =>
    { c with stages := .dflt, featPos := true, step := .adaptive 0.05 0.25 1.0 }
  let V := fun (nm : String) (m : RawMesh) (cs : Array (Point Float)) => mkRegionCase nm m (redistancedState m (voronoiState m cs))
  -- the paper's stages first
  fs := add fs { V "r-sphere3-paper" ico2 three with alphaL := 0.5, step := .fixed 0.0005 }
    "paper's stages: 3 regions on a sphere (Voronoi of 3 centres, 2 triple points), Dirichlet"
  fs := add fs { V "r-sphere4-paper" ico2 tetra with alphaL := 1.0, step := .fixed 0.0003 }
    "paper's stages: 4 regions (tetrahedral partition), Dirichlet"
  fs := add fs { V "r-sphere6-paper" ico2 cube with alphaL := 1.0, energy := .hencky, step := .capped 0.01 0.25 }
    "paper's stages: 6 regions (cubical partition), Hencky, capped step"
  let tc : Array (Point Float) := #[pt 1.3 0.0 0.0, pt 0.0 1.3 0.1, pt (-1.3) 0.1 0.0, pt 0.1 (-1.3) 0.0, pt 1.0 0.0 0.45]
  fs := add fs { V "r-torus5-paper" torusM tc with alphaL := 0.5, energy := .scaledCurv 1.0 0.02, step := .capped 0.01 0.25 }
    "paper's stages: 5 regions on a torus (96... vertices), curvature-scaled Dirichlet"
  let dc : Array (Point Float) := #[pt 0.2 0.25 0.0, pt 0.8 0.3 0.0, pt 0.3 0.8 0.0, pt 0.85 0.85 0.0]
  fs := add fs { V "r-disk4-paper" grid dc with alphaL := 0.5, energy := .scaledA (FloatArray.ofFn' grid.pos.size fun i => 1.0 + 0.5 * Float.sin (7.0 * Float.ofNat i)), step := .capped 0.01 0.25 }
    "paper's stages: 4 regions on a disk with boundary (cuts end on the boundary, triple points), scaled Dirichlet with an explicit a"
  let bands := bandState ico2 #[0.9, 1.9]
  fs := add fs { mkRegionCase "r-bands3-weights-paper" ico2 bands with ws := #[2.0, 0.5, 1.5], betas := #[0.4, 0.0, 0.1], alphaL := 1.0, step := .fixed 0.0005 }
    "paper's stages: 3 latitude bands (no triple point), w = (2, 0.5, 1.5), beta = (0.4, 0, 0.1)"
  let empty := (bandState ico2 #[1.2]).push (FloatArray.ofFn' ico2.pos.size fun _ => 1.0)
  fs := add fs { mkRegionCase "r-empty-region-paper" ico2 empty with alphaL := 1.0, step := .fixed 0.0005 }
    "paper's stages: 3 functions of which the third is positive everywhere (empty region 2); two caps"
  let etaZ := FloatArray.ofFn' ico2.pos.size fun i => 0.2 * (1 + ico2.pos[i]!.z)
  fs := add fs { V "r-sphere3-eta-paper" ico2 three with alphaL := 0.5, eta := some etaZ, step := .fixed 0.0005 }
    "paper's stages: 3 regions, penalty eta = 0.2 (1 + z), gradient rule"
  fs := add fs { V "r-sphere3-eta-paperrule" ico2 three with alphaL := 0.5, eta := some etaZ, penalty := .paperDiagonal, step := .fixed 0.0005 }
    "paper's stages: 3 regions, penalty eta, the paper's diagonal rule"
  let symC : Array (Point Float) := #[pt 0.0 0.9 0.4, pt 0.0 (-0.8) 0.5, pt 0.0 0.05 (-1.0)]
  fs := add fs { mkRegionCase "r-sphere3-symmetry-paper" ico2s (redistancedState ico2s (voronoiState ico2s symC)) with alphaL := 0.5, step := .fixed 0.0005, classes := some (mirrorClasses ico2s) }
    "paper's stages: 3 regions invariant under x -> -x on the symmetric icosphere, symmetry classes"
  fs := add fs { V "r-sphere4-nodal-paper" ico2 tetra with alphaL := 1.0, flux := .nodal, step := .fixed 0.0003 }
    "paper's stages with nodal flux (the paper's evaluation), 4 regions"
  -- the default stages of section 8.5
  fs := add fs (dflt { V "r-sphere3-default" ico2 three with alphaL := 0.5 })
    "default stages (8.5): 3 regions on a sphere"
  fs := add fs (dflt { V "r-sphere4-default" ico2 tetra with alphaL := 1.0, energy := .hencky })
    "default stages: 4 regions, Hencky"
  fs := add fs (dflt { V "r-sphere6-default" ico2 cube with alphaL := 1.0 })
    "default stages: 6 regions"
  fs := add fs (dflt { V "r-torus5-default" torusM tc with alphaL := 0.5 })
    "default stages: 5 regions on a torus"
  fs := add fs (dflt { V "r-disk4-default" grid dc with alphaL := 0.5 })
    "default stages: 4 regions on a disk with boundary (boundary cut points, triple points)"
  fs := add fs (dflt { mkRegionCase "r-bands3-weights-default" ico2 bands with ws := #[2.0, 0.5, 1.5], betas := #[0.4, 0.0, 0.1], alphaL := 1.0 })
    "default stages: 3 bands, per-region weights"
  fs := add fs (dflt { V "r-sphere3-eta-default" ico2 three with alphaL := 0.5, eta := some etaZ })
    "default stages: 3 regions, penalty eta (gradient rule: dw/dn term)"
  fs := add fs (dflt { mkRegionCase "r-empty-region-default" ico2 empty with alphaL := 1.0 })
    "default stages: an empty region"
  return fs

/-- Render a fixture. -/
def RFixture.render (f : RFixture) : String × String := (printCase f.case, runCase f.case)

/-- Write the fixtures. -/
def writeFixtures (dir : System.FilePath) : IO Unit := do
  IO.FS.createDirAll dir
  let mut manifest := ""
  for f in fixtures do
    let (input, output) := f.render
    IO.FS.writeFile (dir / s!"{f.case.name}.in") input
    IO.FS.writeFile (dir / s!"{f.case.name}.out") output
    manifest := manifest ++ s!"{f.case.name}\t{statusOf output}\t{f.what}\n"
    IO.println s!"{f.case.name}: {statusOf output}"
  IO.FS.writeFile (dir / "MANIFEST.tsv") manifest
  IO.println s!"wrote {fixtures.size} fixtures to {dir}"

/-- Check the fixtures; returns the number of differences. -/
def checkFixtures (dir : System.FilePath) : IO Nat := do
  let mut bad := 0
  for f in fixtures do
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
  IO.println s!"region fixtures: {fixtures.size - bad} of {fixtures.size} match"
  return bad

end Ddg.Regions
