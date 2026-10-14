import DdgOracle.Cuts.FixturesDefault

/-!
# Golden fixtures of the robust step statistics (format 3)

`fixtures/cuts-robust/`: two regions, `statistics robust` (and one `extremes` case that pins the
format-3 output of the old statistics).
-/

namespace Ddg.Cuts

/-- A copy of `m` in which vertex `idx` is pushed radially away from the origin by `factor`
(a spike of curvature). -/
def spikeAt (m : RawMesh) (idx : Nat) (factor : Float) : RawMesh :=
  { m with pos := m.pos.mapIdx fun i p => if i == idx then ⟨p.x * factor, p.y * factor, p.z * factor⟩ else p }

/-- The vertex of the mesh whose polar angle is closest to `psi` (and whose `x` is positive). -/
def vertexAtPolarAngle (m : RawMesh) (psi : Float) : Nat := Id.run do
  let mut best := 0
  let mut bd := posInf
  for i in [0:m.pos.size] do
    let p := m.pos[i]!
    let d := (Float.acos (p.z / Float.sqrt (p.x * p.x + p.y * p.y + p.z * p.z)) - psi).abs + (if p.x > 0.2 then 0.0 else 1.0)
    if d < bd then
      bd := d
      best := i
  return best

/-- Topology of a mesh (for the neighbours of a vertex). -/
def triOfMesh (m : RawMesh) : Option Topology :=
  match buildTopology m.pos.size m.faces (connected := false) with
  | .ok t => some t
  | .error _ => none

/-- The vertex nearest to a point. -/
def nearestVertex (m : RawMesh) (q : Point Float) : Nat := Id.run do
  let mut best := 0
  let mut bd := posInf
  for i in [0:m.pos.size] do
    let d := dist m.pos[i]! q
    if d < bd then
      bd := d
      best := i
  return best

/-- The default stages with the robust statistics (format 3). -/
def withRobust (c : CutCase) : CutCase :=
  { withDefaultStages c with format := 3, stats := .robust }

/-- The fixtures of the robust statistics. -/
def robustFixtures : Array Fixture := Id.run do
  let warp := fun (m : RawMesh) => m.mapPos fun p =>
    ⟨p.x * (1 + 0.05 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.04 * Float.cos (4 * p.x - 2 * p.z)), p.z⟩
  let ico2 := warp (icosphere 2)
  let cap := fun (m : RawMesh) (psi0 : Float) => phiOf m fun p => Float.acos p.z - psi0
  -- a face with a much larger curvature density next to the cut (an outlier of omega_c)
  let hot := vertexAtPolarAngle ico2 1.0
  let kf := FloatArray.ofFn' ico2.faces.size fun fi =>
    let f := ico2.faces[fi]!
    if f.a == hot || f.b == hot || f.c == hot then 40.0 * 4 * Ddg.floatPi else 4 * Ddg.floatPi
  let spike := ico2
  let torusM := (torus 12 8 1.0 0.42 fun a b => 0.1 * Float.sin (2 * a) * Float.cos b).jitter 0.005 21
  let grid := (gridDisk 12 11 fun x y => 0.3 * Float.sin (5 * x) * Float.cos (4 * y) + 0.2 * x * x).jitter 0.004 22
  let mut fs : Array Fixture := #[]
  let add := fun (fs : Array Fixture) (c : CutCase) (what : String) => fs.push { case := c, what }
  let capPhi := cap ico2 1.0
  fs := add fs { withRobust (mkCase "rb-ico2-cap-adaptive" ico2 capPhi) with alphaL := 2.0 }
    "robust statistics, default stages, adaptive step (uniform-ish cut: robust and extremes nearly agree)"
  fs := add fs { withRobust (mkCase "rb-ico2-cap-capped" ico2 capPhi) with alphaL := 2.0, step := .capped 0.01 0.25 }
    "robust statistics, default stages, capped step"
  fs := add fs { { mkCase "rb-ico2-cap-paper-adaptive" ico2 capPhi with format := 3, stats := .robust, alphaL := 2.0 } with
      step := .adaptive 0.05 0.25 1.0 }
    "robust statistics with the paper's stages (length step on the mesh: the larger of the two m-hats)"
  fs := add fs { withRobust (mkCase "rb-spike-adaptive" spike capPhi) with curv := .perFace kf, alphaL := 2.0 }
    "robust, default stages: one face next to the cut with a 40 times larger curvature density (a spike of curvature): outlier omega_c, cut points limited"
  fs := add fs { withRobust (mkCase "rb-spike-capped" spike capPhi) with curv := .perFace kf, alphaL := 2.0, step := .capped 0.01 0.25 }
    "the same with the capped step (limit |r_c|)"
  fs := add fs { { withRobust (mkCase "rb-spike-paper" spike capPhi) with lengthStep := .onMesh, extension := .harmonic, speed := .unit, featPos := false } with curv := .perFace kf, alphaL := 2.0 }
    "robust with the paper's stages on the spike (both m-hats, limit on r = s)"
  fs := add fs { withRobust (mkCase "rb-spike-extremes" spike capPhi) with stats := .extremes, curv := .perFace kf, alphaL := 2.0 }
    "the same spike with the extremes in format 3: pins the old statistics and the extra lines"
  fs := add fs { withRobust (mkCase "rb-grid-boundary-adaptive" grid (phiOf grid fun p => p.x + 0.3 * p.y - 0.65)) with
      energy := .scaledA (FloatArray.ofFn' grid.pos.size fun i => 1.0 + 0.5 * Float.sin (7.0 * Float.ofNat i)), alphaL := 0.3 }
    "robust, default stages, cut ending on the boundary, scaled energy"
  fs := add fs { withRobust (mkCase "rb-torus-hencky" torusM (phiOf torusM fun p => Float.cos (3 * Float.atan2 p.y p.x) - 0.2)) with
      energy := .hencky, alphaL := 0.5 }
    "robust, default stages, torus, Hencky"
  -- D30: vertices on the cut to rounding, so that crossings are moved by delta
  let flat := (polar 10 28 0.0 1.0 fun _ _ => 0.0).jitter 0.004 24
  let rho := 0.3 * Float.sqrt Ddg.floatPi
  let radial := phiOf flat fun p => Float.sqrt (p.x * p.x + p.y * p.y) - rho
  -- the six vertices nearest to the circle lie on it to rounding (alternating sides)
  let order := (Array.range flat.pos.size).qsort fun i j => (radial.get! i).abs < (radial.get! j).abs
  let snapped := Id.run do
    let mut phi := radial
    for k in [0:6] do
      phi := phi.set! order[k]! (if k % 2 == 0 then 1e-15 else -1e-15)
    return phi
  -- every vertex of the inside that has an outside neighbour lies on the cut to rounding: every crossing is moved
  let allMoved := Id.run do
    let mut phi := radial
    match triOfMesh flat with
    | none => return phi
    | some topo =>
      for v in [0:flat.pos.size] do
        if radial.get! v ≤ 0 then
          let mut touches := false
          for w in topo.neighbors[v]! do
            if radial.get! w > 0 then touches := true
          if touches then phi := phi.set! v (-1e-15)
    return phi
  fs := add fs { { withRobust (mkCase "rb-moved-some" flat snapped) with format := 4 } with alphaL := 1.0 }
    "format 4: six vertices on the start circle to rounding, their crossings are moved by delta: the faces with a moved crossing do not count in the mean slope (gr lines)"
  fs := add fs { { withRobust (mkCase "rb-moved-all" flat allMoved) with format := 4 } with alphaL := 1.0 }
    "format 4: every crossing of the closed cut is moved: no face counts and the mean slope is 1"
  return fs

end Ddg.Cuts
