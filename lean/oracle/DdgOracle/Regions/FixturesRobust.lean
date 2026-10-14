import DdgOracle.Regions.Fixtures
import DdgOracle.Cuts.FixturesRobust

/-!
# Golden fixtures of the robust statistics for N regions (format 2)

`fixtures/regions-robust/`.
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- The fixtures. -/
def robustFixtures : Array RFixture := Id.run do
  let warp := fun (m : RawMesh) => m.mapPos fun p =>
    ⟨p.x * (1 + 0.05 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.04 * Float.cos (4 * p.x - 2 * p.z)), p.z⟩
  let ico2 := warp (icosphere 2)
  let grid := (gridDisk 14 13 fun x y => 0.3 * Float.sin (5 * x) * Float.cos (4 * y) + 0.2 * x * x).jitter 0.004 22
  let pt := fun (x y z : Float) => (⟨x, y, z⟩ : Point Float)
  let nrm := fun (p : Point Float) =>
    let n := Float.sqrt (p.x * p.x + p.y * p.y + p.z * p.z)
    (⟨p.x / n, p.y / n, p.z / n⟩ : Point Float)
  let tetra : Array (Point Float) := #[nrm (pt 1.0 1.13 0.87), nrm (pt 1.07 (-0.93) (-1.1)), nrm (pt (-1.1) 0.97 (-0.9)), nrm (pt (-0.9) (-1.05) 1.13)]
  let three : Array (Point Float) := #[nrm (pt 1.0 0.1 0.2), nrm (pt (-0.6) 0.9 (-0.1)), nrm (pt (-0.5) (-0.8) 0.3)]
  let dc : Array (Point Float) := #[pt 0.2 0.25 0.0, pt 0.8 0.3 0.0, pt 0.3 0.8 0.0, pt 0.85 0.85 0.0]
  let mid := nrm (pt 0.4 1.0 0.1)
  let spike := spikeAt ico2 (nearestVertex ico2 mid) 2.4
  let V := fun (nm : String) (m : RawMesh) (cs : Array (Point Float)) =>
    { mkRegionCase nm m (redistancedState m (voronoiState m cs)) with format := 2, stats := .robust }
  let dflt := fun (c : RegionCase) =>
    { c with stages := .dflt, featPos := true, step := .adaptive 0.05 0.25 1.0 }
  let paperAd := fun (c : RegionCase) => { c with step := .adaptive 0.05 0.25 1.0 }
  let mut fs : Array RFixture := #[]
  let add := fun (fs : Array RFixture) (c : RegionCase) (what : String) => fs.push { case := c, what }
  fs := add fs (paperAd { V "rr-sphere4-paper" ico2 tetra with alphaL := 1.0 })
    "robust statistics, paper's stages, 4 regions (m-hat for s over the cut points and for s + lambda0 over the boundary nodes)"
  fs := add fs (dflt { V "rr-sphere4-default" ico2 tetra with alphaL := 1.0 })
    "robust statistics, default stages, 4 regions"
  fs := add fs (paperAd { V "rr-disk4-paper" grid dc with alphaL := 0.5 })
    "robust, paper's stages, 4 regions on a disk with boundary"
  fs := add fs (dflt { mkRegionCase "rr-sphere3-spike-default" spike (redistancedState spike (voronoiState spike three)) with format := 2, stats := .robust, alphaL := 1.0 })
    "robust, default stages, 3 regions with a spike of curvature at one vertex"
  fs := add fs (paperAd { mkRegionCase "rr-sphere3-spike-paper" spike (redistancedState spike (voronoiState spike three)) with format := 2, stats := .robust, alphaL := 1.0 })
    "robust, paper's stages, 3 regions with the spike (local limit on s_v)"
  -- D30: vertices on an interface to rounding (phi_i(a) - phi_j(a) = -1e-15): their crossings are moved by delta
  let st4 := redistancedState ico2 (voronoiState ico2 tetra)
  let moved := Id.run do
    let nvv := ico2.pos.size
    let mut st := st4
    let mut count := 0
    match triOfMesh ico2 with
    | none => return st
    | some topo =>
      for v in [0:nvv] do
        if count < 6 then
          let mut order : Array Nat := (Array.range st.size).qsort fun a b => (st[a]!).get! v < (st[b]!).get! v
          let i := order[0]!
          let j := order[1]!
          -- the vertex has a neighbour in another region
          let mut other := false
          for w in topo.neighbors[v]! do
            let mut best := 0
            for k in [1:st.size] do
              if (st[k]!).get! w < (st[best]!).get! w then best := k
            if best != i then other := true
          if other && v % 3 == 0 then
            st := st.set! i ((st[i]!).set! v ((st[j]!).get! v - 1e-15))
            count := count + 1
    return st
  fs := add fs (dflt { mkRegionCase "rr-moved-default" ico2 moved with format := 3, stats := .robust, alphaL := 1.0 })
    "format 3: six vertices on an interface to rounding, crossings moved by delta, default stages: those faces do not count in the mean slope"
  fs := add fs (paperAd { mkRegionCase "rr-moved-paper" ico2 moved with format := 3, stats := .robust, alphaL := 1.0 })
    "the same with the paper's stages"
  return fs

/-- Fixtures that are adopted as files (inputs built elsewhere): name and description. -/
def adoptedFixtures : Array (String × String) :=
  #[("triple-chain-paper", "a region (one vertex whose neighbours alternate between two other regions) whose boundary is a closed chain through triple-point copies only; paper's stages, robust"),
    ("triple-chain-default", "the same, default stages")]

/-- Run an adopted input file. -/
def runAdopted (dir : System.FilePath) (name : String) : IO (Option String) := do
  let path := dir / s!"{name}.in"
  if !(← path.pathExists) then return none
  let text ← IO.FS.readFile path
  let lines := (text.splitOn "\n").filterMap fun l =>
    let t := tokens l
    if t.size == 0 || t[0]! == "case" || t[0]! == "end" then none else some t
  match parseBody name lines.toArray with
  | .error e => return some (printError name e)
  | .ok c => return some (runCase c)

/-- Write the fixtures (the adopted ones are only listed in the manifest and checked by `checkRobust`). -/
def writeRobust (dir : System.FilePath) : IO Unit := do
  writeListR robustFixtures dir
  let mut manifest ← IO.FS.readFile (dir / "MANIFEST.tsv")
  for (n, what) in adoptedFixtures do
    match ← runAdopted dir n with
    | some out =>
      IO.FS.writeFile (dir / s!"{n}.out") out
      manifest := manifest ++ s!"{n}\t{statusOf out}\t{what}\n"
    | none => IO.println s!"{n}: input file missing"
  IO.FS.writeFile (dir / "MANIFEST.tsv") manifest

/-- Check the fixtures, the adopted ones included; returns the number of differences. -/
def checkRobust (dir : System.FilePath) : IO Nat := do
  let mut bad ← checkListR robustFixtures dir
  for (n, _) in adoptedFixtures do
    match ← runAdopted dir n with
    | none =>
      IO.println s!"FAIL  {n}.in is missing"
      bad := bad + 1
    | some out =>
      let stored ← IO.FS.readFile (dir / s!"{n}.out")
      match compareOutputs stored out with
      | none => IO.println s!"ok    {n}"
      | some d =>
        IO.println s!"FAIL  {n}: {d}"
        bad := bad + 1
  return bad

end Ddg.Regions
