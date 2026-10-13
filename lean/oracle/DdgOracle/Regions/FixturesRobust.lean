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
  return fs

/-- Write the fixtures. -/
def writeRobust (dir : System.FilePath) : IO Unit := writeListR robustFixtures dir

/-- Check the fixtures. -/
def checkRobust (dir : System.FilePath) : IO Nat := checkListR robustFixtures dir

end Ddg.Regions
