import DdgOracle.Cuts.Selftest
import DdgOracle.Regions.Gen

/-!
# `regions-selftest`
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- A step that must succeed. -/
def runOkR (c : RegionCase) : Option RStep :=
  match runRegionStep c with
  | .ok r => some r
  | .error _ => none

/-- Largest absolute difference between two lists of functions. -/
def maxDiffs (a b : Array FloatArray) : Float := Id.run do
  let mut d := fzero
  for k in [0:a.size] do
    d := max d (maxDiff a[k]! b[k]!)
  return d

/-- Run all checks; returns the number of failures. -/
def selftest : IO Nat := do
  let tally ← IO.mkRef ({} : Tally)
  let pi := floatPi
  let ico2 := icosphere 2
  let ico3 := icosphere 3

  -- 1. N = 2 through the N-region code equals the two-region code, for both stage sets
  let torusM := (torus 12 8 1.0 0.42 fun a b => 0.1 * Float.sin (2 * a) * Float.cos b).jitter 0.005 21
  let grid := (gridDisk 12 11 fun x y => 0.3 * Float.sin (5 * x) * Float.cos (4 * y) + 0.2 * x * x).jitter 0.004 22
  let ann := (polar 5 16 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.25 * y * y).jitter 0.005 23
  let etaZ := FloatArray.ofFn' ico2.pos.size fun i => 0.2 * (1 + ico2.pos[i]!.z)
  let capPhi := phiOf ico2 fun p => Float.acos p.z - 1.0
  let cases : List (String × CutCase) :=
    [("sphere cap, Dirichlet", mkCase "a" ico2 capPhi),
     ("sphere cap, Hencky, w = (2, .5), beta = (.7, .2)", { mkCase "b" ico2 capPhi with energy := .hencky, w := (2.0, 0.5), beta := (0.7, 0.2) }),
     ("sphere cap, eta, paper rule", { mkCase "c" ico2 capPhi with eta := some etaZ, penalty := .paperDiagonal, alphaL := 0.5 }),
     ("sphere cap, eta, gradient rule, curvature-scaled", { mkCase "d" ico2 capPhi with eta := some etaZ, energy := .scaledCurv 1.0 0.02, alphaL := 0.5 }),
     ("torus band", mkCase "e" torusM (phiOf torusM fun p => Float.cos (3 * Float.atan2 p.y p.x) - 0.2)),
     ("grid disk, cut ends on the boundary, scaled a", { mkCase "f" grid (phiOf grid fun p => p.x + 0.3 * p.y - 0.65) with
        energy := .scaledA (FloatArray.ofFn' grid.pos.size fun i => 1.0 + 0.5 * Float.sin (7.0 * Float.ofNat i)) }),
     ("annulus, radial cut, classes none", mkCase "g" ann (phiOf ann fun p => p.x - 0.1))]
  for (nm, c0) in cases do
    for (stNm, st) in [("paper", Stages.paper), ("default", Stages.dflt)] do
      let c1 : CutCase :=
        if st == .paper then { c0 with step := .fixed 0.001 }
        else { c0 with format := 2, lengthStep := .onCut, extension := .local, speed := .meanGradient
                       featPos := true, step := .adaptive 0.05 0.25 1.0 }
      match runOk c1, runOkR (ofCut c1 st) with
      | some a, some b =>
        let dPhi := maxDiff a.phiNew (b.phiNew[0]!)
        let dPhi1 := maxDiff (FloatArray.ofFn' a.phiNew.size fun i => -(a.phiNew.get! i)) (b.phiNew[1]!)
        let dTau := (a.tau - b.tau).abs / a.tau
        let dE := (a.ev.total - b.ev.total).abs / (1 + a.ev.total.abs)
        let dEnew := (a.evNew.total - b.evNew.total).abs / (1 + a.evNew.total.abs)
        let mut dS := fzero
        for k in [0:a.ev.cm.nc] do dS := max dS ((a.ev.sRaw.get! k - b.ev.sCut.get! k).abs / (1 + (a.ev.sRaw.get! k).abs))
        check tally s!"N = 2 through the N-region code equals the two-region code ({stNm} stages): {nm}"
          (dPhi ≤ 1e-12 && dPhi1 ≤ 1e-12 && dTau ≤ 1e-12 && dE ≤ 1e-12 && dEnew ≤ 1e-12 && dS ≤ 1e-12)
          s!"phi {sci dPhi} (-phi: {sci dPhi1}) tau {sci dTau} E {sci dE} E_new {sci dEnew} s {sci dS}"
      | _, _ => check tally s!"N = 2 ({stNm}): {nm}: steps run" false

  -- 2. cut mesh invariants for networks; projection; junction loads on a flat mesh
  let warp := fun (m : RawMesh) => m.mapPos fun p =>
    ⟨p.x * (1 + 0.05 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.04 * Float.cos (4 * p.x - 2 * p.z)), p.z⟩
  let wico := warp ico2
  let pt := fun (x y z : Float) => (⟨x, y, z⟩ : Point Float)
  let nrm := fun (p : Point Float) =>
    let n := Float.sqrt (p.x * p.x + p.y * p.y + p.z * p.z)
    (⟨p.x / n, p.y / n, p.z / n⟩ : Point Float)
  let tetra : Array (Point Float) := #[nrm (pt 1.0 1.13 0.87), nrm (pt 1.07 (-0.93) (-1.1)), nrm (pt (-1.1) 0.97 (-0.9)), nrm (pt (-0.9) (-1.05) 1.13)]
  let cube : Array (Point Float) := #[pt 1.0 0.05 0.0, pt (-1.0) 0.0 0.05, pt 0.05 1.0 0.0, pt 0.0 (-1.0) (-0.05), pt 0.0 0.05 1.0, pt (-0.05) 0.0 (-1.0)]
  let three : Array (Point Float) := #[nrm (pt 1.0 0.1 0.2), nrm (pt (-0.6) 0.9 (-0.1)), nrm (pt (-0.5) (-0.8) 0.3)]
  let diskM := (gridDisk 14 13 fun x y => 0.3 * Float.sin (5 * x) * Float.cos (4 * y) + 0.2 * x * x).jitter 0.004 22
  let dc : Array (Point Float) := #[pt 0.2 0.25 0.0, pt 0.8 0.3 0.0, pt 0.3 0.8 0.0, pt 0.85 0.85 0.0]
  let netCases : List (String × RawMesh × Array (Point Float)) :=
    [("sphere, 3 regions", wico, three), ("sphere, 4 regions", wico, tetra), ("sphere, 6 regions", wico, cube),
     ("disk with boundary, 4 regions", diskM, dc)]
  for (nm, m, cs) in netCases do
    let c := mkRegionCase "net" m (redistancedState m (voronoiState m cs))
    match runOkR c with
    | none => check tally s!"network, {nm}: step runs" false
    | some r =>
      let net := r.ev.net
      let mut area := fzero
      for t in net.tris do area := area + t.area
      let mut chiS : Int := 0
      for rv in r.ev.views do
        chiS := chiS + regionChi rv.cm rv.X
      let chiM : Int := (match buildTopology m.pos.size m.faces with | .ok t => t.chi | .error _ => 0)
      let lhs := chiS - (net.nP : Int) - 2 * (net.nM : Int) + (net.segs.size : Int)
      -- projection: smallest = -second, differences kept
      let mut projErr := fzero
      let mut diffErr := fzero
      for i in [0:r.base.nv] do
        if r.redist.adj[i]! then
          let vals := (List.range cs.size).map fun k => (r.redist.projected[k]!).get! i
          let sorted := vals.toArray.qsort (· < ·)
          projErr := max projErr ((sorted[0]! + sorted[1]!).abs)
          for k in [1:cs.size] do
            let d0 := (r.phiSym[k]!).get! i - (r.phiSym[0]!).get! i
            let d1 := (r.redist.projected[k]!).get! i - (r.redist.projected[0]!).get! i
            diffErr := max diffErr ((d0 - d1).abs)
      check tally s!"network, {nm}: {net.nP} cut points, {net.nM} triple points; areas sum to 1; Euler characteristics; projection"
        ((area - 1).abs ≤ 1e-12 && lhs == chiM && projErr ≤ 1e-15 && diffErr ≤ 1e-15)
        s!"area-1={sci (area - 1)} chi {lhs} vs {chiM}; smallest+second={sci projErr}; differences kept {sci diffErr}"

  -- 3. flat mesh: junction loads of a balanced junction sum to zero, opposite turning along a curve
  let flatM := (gridDisk 16 15 fun _ _ => 0.0).jitter 0.004 5
  let flatC : Array (Point Float) := #[pt 0.3 0.3 0.0, pt 0.75 0.4 0.0, pt 0.45 0.8 0.0]
  match runOkR { mkRegionCase "flatnet" flatM (redistancedState flatM (voronoiState flatM flatC)) with
                 stages := .dflt, alphaL := 1.0, step := .fixed 0.0001 } with
  | none => check tally "flat network: step runs" false
  | some r =>
    let net := r.ev.net
    let nv := net.nv
    let mut loadTriple := fzero
    let mut loadCut := fzero
    -- accumulate per geometric node
    let mut tripleSum : Array Float := Array.replicate net.nM fzero
    let mut cutSum : Array Float := Array.replicate net.nP fzero
    let mut cutIsBd : Array Bool := Array.replicate net.nP false
    for rv in r.ev.views do
      for j in [0:rv.toCopy.size] do
        let d := (r.dflt[rv.X]!)[j]!
        let (g, _) := geoOfCopy net rv.toCopy[j]!
        let load := d[3]! - d[1]!
        if g ≥ nv + net.nP then tripleSum := tripleSum.set! (g - nv - net.nP) (tripleSum[g - nv - net.nP]! + load)
        else
          cutSum := cutSum.set! (g - nv) (cutSum[g - nv]! + load)
          cutIsBd := cutIsBd.set! (g - nv) (d[3]! < 3.0)
    for x in tripleSum do loadTriple := max loadTriple x.abs
    for p in [0:net.nP] do
      if !cutIsBd[p]! then loadCut := max loadCut (cutSum[p]!).abs
    check tally "flat network: the loads gamma - beta of the regions sum to zero at every triple point and every interior cut point"
      (net.nM > 0 && loadTriple ≤ 1e-12 && loadCut ≤ 1e-12)
      s!"{net.nM} triple points: max |sum| {sci loadTriple}; cut points: {sci loadCut}"

  -- 4. chain sums per region, opposite rates on the two sides of a curve without junctions
  let bands3 := bandState wico #[0.9, 1.9]
  match runOkR { mkRegionCase "bands" wico bands3 with stages := .dflt, alphaL := 1.0, step := .adaptive 0.05 0.25 1.0 } with
  | none => check tally "bands, default stages: step runs" false
  | some r =>
    let net := r.ev.net
    let mut sumErr := fzero
    for rv in r.ev.views do
      let mut st := fzero
      let mut sl := fzero
      let mut sr := fzero
      let mut srhs := fzero
      for j in [0:rv.toCopy.size] do
        let d := (r.dflt[rv.X]!)[j]!
        st := st + d[4]!
        sl := sl + d[0]! * d[5]!
        sr := sr + d[0]! * d[6]!
        srhs := srhs + d[0]! * ((r.ev.sNode[rv.X]!).get! j + d[5]!)
      sumErr := max sumErr (max ((sl - st).abs) ((sr - srhs).abs))
    -- opposite rates: the two copies of every cut point
    let mut oppErr := fzero
    for p in [0:net.nP] do
      let cr := net.cps[p]!
      let ra := net.region[cr.a]!
      let rb := net.region[cr.b]!
      let ia := (r.ev.views[ra]!).localOf.getD (net.copyIdx (net.nv + p) ra) 0
      let ib := (r.ev.views[rb]!).localOf.getD (net.copyIdx (net.nv + p) rb) 0
      let da := (r.dflt[ra]!)[ia]!
      let db := (r.dflt[rb]!)[ib]!
      oppErr := max oppErr ((da[6]! + db[6]!).abs)
    check tally "bands, default stages: sum l lam0 = sum Theta, sum l r = sum l (s + lam0) per region; rates opposite on the two sides of a curve"
      (sumErr ≤ 1e-12 && oppErr ≤ 1e-12) s!"sums {sci sumErr}; opposite {sci oppErr}"

  -- 5. energy decreases over one small step for N = 3 (both stage sets)
  for (stNm, st) in [("paper", Stages.paper), ("default", Stages.dflt)] do
    let c := { mkRegionCase "dec3" wico (redistancedState wico (voronoiState wico three)) with
      stages := st, alphaL := 1.0, step := .fixed 0.0002 }
    match runOkR c with
    | none => check tally s!"energy decrease N = 3 ({stNm}): step runs" false
    | some r =>
      check tally s!"energy decreases along one small step, 3 regions with 2 triple points ({stNm} stages)"
        (r.evNew.total < r.ev.total) s!"E {sci r.ev.total} -> {sci r.evNew.total} (change {sci (r.evNew.total - r.ev.total)})"

  -- 6. finite differences for an N-region state (3 bands): shifting phi_0 by +-e
  let c6 := { mkRegionCase "fd3" ico3 (bandState ico3 #[0.9, 1.9]) with alphaL := 0.0 }
  match buildBase c6.toCut with
  | .error e => check tally "finite differences, 3 bands" false e
  | .ok b =>
    let p0 := c6.phis.map fun p => FloatArray.ofFn' p.size fun i => b.scale * p.get! i
    let shifted := fun (e : Float) => p0.mapIdx fun k p => if k == 0 then FloatArray.ofFn' p.size fun i => p.get! i + e else p
    let e := 0.0005
    match evaluateNet c6 b p0, evaluateNet c6 b (shifted e), evaluateNet c6 b (shifted (-e)) with
    | .ok ev, .ok evP, .ok evM =>
      let net := ev.net
      let mut dual := zeros net.nP
      for sg in net.segs do
        dual := dual.set! (sg.g1 - net.nv) (dual.get! (sg.g1 - net.nv) + sg.len / 2)
        dual := dual.set! (sg.g2 - net.nv) (dual.get! (sg.g2 - net.nv) + sg.len / 2)
      let mut pred := fzero
      for p in [0:net.nP] do
        if ev.lowHigh[p]! == (0, 1) then pred := pred - e * ev.sCut.get! p * dual.get! p
      let fd := evP.total - evM.total
      let rel := (fd - pred).abs / pred.abs
      check tally "finite differences of E for 3 latitude bands (shift of phi_0 moves the interface 0|1 by e/2), icosphere 3 (10%)"
        (rel ≤ 0.10) s!"E(+e)-E(-e)={sci fd} predicted {sci pred} (rel {sci rel})"
    | _, _, _ => check tally "finite differences, 3 bands: evaluations run" false

  let t ← tally.get
  IO.println s!"regions-selftest: {t.passed} passed, {t.failed} failed"
  return t.failed

end Ddg.Regions
