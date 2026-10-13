import DdgOracle.Cuts.Fixtures

/-!
# Golden fixtures of the default stages (format 2)

`fixtures/cuts-default/`: the same kinds of cases as `fixtures/cuts/`, with the default stages of
SPEC §4 (length step on the cut, extension from the cut, mean-gradient scaling, adaptive step,
positive feature-size smoothing), and mixed variants that change one choice at a time.
-/

namespace Ddg.Cuts

/-- The default stages of SPEC §4 (`FlowOptions::new`) for a case. -/
def withDefaultStages (c : CutCase) : CutCase :=
  { c with format := 2, lengthStep := .onCut, extension := .local, speed := .meanGradient
           featPos := true, step := .adaptive 0.05 0.25 1.0, flux := .smoothed 0.5, passes := 0 }

/-- The fixtures of the default stages. -/
def defaultFixtures : Array Fixture := Id.run do
  let warp := fun (m : RawMesh) => m.mapPos fun p =>
    ⟨p.x * (1 + 0.05 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.04 * Float.cos (4 * p.x - 2 * p.z)), p.z⟩
  let ico1 := warp (icosphere 1)
  let ico2 := warp (icosphere 2)
  let ico2s := icosphere 2
  let cap := fun (m : RawMesh) (psi0 : Float) => phiOf m fun p => Float.acos p.z - psi0
  let torusM := (torus 12 8 1.0 0.42 fun a b => 0.1 * Float.sin (2 * a) * Float.cos b).jitter 0.005 21
  let grid := (gridDisk 12 11 fun x y => 0.3 * Float.sin (5 * x) * Float.cos (4 * y) + 0.2 * x * x).jitter 0.004 22
  let ann := (polar 5 16 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.25 * y * y).jitter 0.005 23
  let flat := (polar 10 28 0.0 1.0 fun _ _ => 0.0).jitter 0.004 24
  let two := RawMesh.union ico1 ico1 ⟨3.0, 0.0, 0.0⟩
  let n1 := ico1.pos.size
  let n2 := ico2.pos.size
  let D := fun (name : String) (m : RawMesh) (phi : FloatArray) => withDefaultStages (mkCase name m phi)
  let mut fs : Array Fixture := #[]
  let add := fun (fs : Array Fixture) (c : CutCase) (what : String) => fs.push { case := c, what }
  let latPhi := cap ico2 1.0
  let tiltPhi := tilted ico2 0.3 0.2 0.93 0.15
  -- default column
  fs := add fs { D "d-ico1-tilted-dirichlet" ico1 (tilted ico1 0.3 0.2 0.93 0.1) with alphaL := 0.5 }
    "default stages: icosphere level 1, tilted plane, Dirichlet, alpha_L = 0.5"
  fs := add fs { D "d-ico2-cap-dirichlet" ico2 latPhi with alphaL := 2.0 }
    "default stages: geodesic cap on icosphere level 2, Dirichlet, alpha_L = 2"
  fs := add fs { D "d-ico2-tilted-hencky" ico2 tiltPhi with energy := .hencky, alphaL := 1.0 }
    "default stages: tilted plane, Hencky"
  fs := add fs { D "d-ico2-cap-scaled-curv" ico2 latPhi with energy := .scaledCurv 1.0 0.02, alphaL := 0.2 }
    "default stages: geodesic cap, curvature-scaled Dirichlet (featurelap positive)"
  fs := add fs { D "d-torus-band" torusM (phiOf torusM fun p => Float.cos (3 * Float.atan2 p.y p.x) - 0.2) with alphaL := 0.5 }
    "default stages: torus, two-band level set"
  fs := add fs { D "d-grid-boundary-scaled-a" grid (phiOf grid fun p => p.x + 0.3 * p.y - 0.65) with
      energy := .scaledA (FloatArray.ofFn' grid.pos.size fun i => 1.0 + 0.5 * Float.sin (7.0 * Float.ofNat i)), alphaL := 0.3 }
    "default stages: grid disk, a straight cut ending on the boundary at both ends (boundary cut points), scaled Dirichlet"
  fs := add fs { D "d-annulus-radial" ann (phiOf ann fun p => p.x - 0.1) with alphaL := 0.3 }
    "default stages: annulus, a cut across the ring ending on both boundary circles"
  fs := add fs { D "d-annulus-ring-hencky" ann (phiOf ann fun p => Float.sqrt (p.x * p.x + p.y * p.y) - 0.7) with
      energy := .hencky, alphaL := 0.3 }
    "default stages: annulus, closed circular cut, Hencky"
  let etaZ := FloatArray.ofFn' n2 fun i => 0.2 * (1 + ico2.pos[i]!.z)
  fs := add fs { D "d-ico2-eta-gradient" ico2 latPhi with eta := some etaZ, alphaL := 1.0 }
    "default stages with a penalty eta = 0.2 (1 + z), gradient rule (includes the dw/dn term)"
  fs := add fs { D "d-ico2-eta-paper" ico2 latPhi with eta := some etaZ, penalty := .paperDiagonal, alphaL := 1.0 }
    "the same with the paper's rule (no dw/dn term)"
  let tiltedCap := phiOf ico2s fun p => polarAngleFrom' ⟨0.15, 0.05, 1.0⟩ p - 1.0
  fs := add fs { D "d-ico2-symmetry" ico2s tiltedCap with classes := some (mirrorClasses ico2s), alphaL := 1.0 }
    "default stages, symmetry classes (mirror x -> -x), cap about a tilted axis"
  fs := add fs { D "d-ico2-region-weights" ico2 latPhi with w := (2.0, 0.5), beta := (0.7, 0.2), alphaL := 1.0 }
    "default stages, w = (2, 0.5), beta = (0.7, 0.2)"
  fs := add fs { D "d-flat-disk-circle" flat (phiOf flat fun p => Float.sqrt (p.x * p.x + p.y * p.y) - 0.3 * Float.sqrt Ddg.floatPi) with
      alphaL := 1.0 }
    "default stages: flat disk, circular cut, Gaussian curvature: pure curve shortening"
  let phiTwo := FloatArray.ofFn' two.pos.size fun i => if i < n1 then Float.acos two.pos[i]!.z - 1.0 else -1.0
  fs := add fs { D "d-two-spheres-untouched" two phiTwo with energy := .hencky, alphaL := 1.0 }
    "default stages: two disjoint spheres, one untouched (region 0) component"
  fs := add fs { { D "d-ico1-phi-zero-unit-speed" ico1 (phiOf ico1 fun p => p.z) with alphaL := 1.0 } with speed := .unit }
    "default stages except speed unit: phi = z exactly 0 at vertices (mean gradient would be 0 there)"
  -- mixed variants, one change from the default at a time
  let base := D "m" ico2 latPhi
  fs := add fs { base with name := "m-length-onmesh", alphaL := 1.0, lengthStep := .onMesh }
    "default except the length step on the mesh (extension local, mean gradient, adaptive)"
  fs := add fs { base with name := "m-extension-harmonic", alphaL := 1.0, extension := .harmonic }
    "default except the harmonic extension"
  fs := add fs { base with name := "m-speed-unit", alphaL := 1.0, speed := .unit }
    "default except ḡ = 1"
  fs := add fs { base with name := "m-step-capped", alphaL := 1.0, step := .capped 0.01 0.25 }
    "default except the capped step rule"
  fs := add fs { base with name := "m-step-fixed", alphaL := 1.0, step := .fixed 0.0005 }
    "default except a fixed step"
  fs := add fs { base with name := "m-featurelap-cotan", alphaL := 0.2, energy := .scaledCurv 1.0 0.02, featPos := false }
    "default except the cotan Laplacian in the feature size smoothing (scaled_curv)"
  fs := add fs { base with name := "m-flux-nodal", alphaL := 1.0, flux := .nodal, passes := 1 }
    "default except nodal flux and one sliver-smoothing pass (the paper's evaluation)"
  fs := add fs { { mkCase "m-paper-stages-format2" ico2 latPhi with format := 2, alphaL := 1.0, step := .fixed 0.001 } with
      featPos := false }
    "format 2 with all of the paper's stages (equals the same case in format 1 up to the extra lines)"
  return fs

end Ddg.Cuts
