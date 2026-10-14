import DdgOracle.Regions.FixturesRobust
import DdgOracle.Cuts.SelftestRobust
import DdgOracle.Regions.Selftest

/-!
# Self-tests of the robust step statistics for N regions
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- Run the checks; returns the number of failures. -/
def selftestRobust : IO Nat := do
  let tally ← IO.mkRef ({} : Tally)
  let ico2 := icosphere 2
  let capPhi := phiOf ico2 fun p => Float.acos p.z - 1.0
  let kf := FloatArray.ofFn' ico2.faces.size fun fi =>
    let f := ico2.faces[fi]!
    let hot := vertexAtPolarAngle ico2 1.0
    if f.a == hot || f.b == hot || f.c == hot then 40.0 * 4 * floatPi else 4 * floatPi
  -- N = 2 with the robust statistics against the two-region code, both stage sets
  for (stNm, st) in [("paper", Stages.paper), ("default", Stages.dflt)] do
    for (cn, step) in [("adaptive", StepRule.adaptive 0.05 0.25 1.0), ("capped", StepRule.capped 0.01 0.25)] do
      let c1 : CutCase :=
        if st == .paper then
          { mkCase "r2" ico2 capPhi with format := 3, stats := .robust, curv := .perFace kf, alphaL := 2.0, step }
        else { withRobust (mkCase "r2" ico2 capPhi) with curv := .perFace kf, alphaL := 2.0, step }
      let c2 := { ofCut c1 st with format := 2, stats := .robust }
      match runOk c1, runOkR c2 with
      | some a, some b =>
        let dPhi := maxDiff a.phiNew (b.phiNew[0]!)
        let dTau := (a.tau - b.tau).abs / a.tau
        let st2 := a.stage2.get!
        let dOm := (st2.omega - b.omega).abs / (1 + st2.omega.abs)
        let dM := (st2.mused - b.mhat.2.2).abs / (1 + st2.mused.abs)
        check tally s!"N = 2, robust statistics, {stNm} stages, {cn} step: N-region code equals the two-region code"
          (dPhi ≤ 1e-12 && dTau ≤ 1e-12 && dOm ≤ 1e-12 && dM ≤ 1e-9)
          s!"phi {sci dPhi} tau {sci dTau} omega {sci dOm} mhat {sci dM}; limited {(st2.factor.toList.filter (· < 1)).length}"
      | _, _ => check tally s!"N = 2 robust ({stNm}, {cn}): steps run" false
  -- D30 through the N-region code: N = 2 with moved crossings equals the two-region code
  for f in Ddg.Cuts.robustFixtures do
    if f.case.name == "rb-moved-some" || f.case.name == "rb-moved-all" then
      let c2 := { ofCut { f.case with format := 3 } .dflt with format := 3, stats := .robust }
      match runOk { f.case with format := 3 }, runOkR c2 with
      | some a, some b =>
        let d := maxDiff a.phiNew (b.phiNew[0]!)
        check tally s!"{f.case.name} through the N-region code (N = 2, default stages, moved crossings): same new state as the two-region code"
          (d ≤ 1e-12) s!"difference {sci d}"
      | _, _ => check tally s!"{f.case.name}: N = 2 steps run" false
  -- local limit on networks: motion of every boundary copy <= max_motion h_c
  for f in robustFixtures.filter (·.case.stages == .dflt) do
    match runOkR f.case with
    | none => check tally s!"{f.case.name}: step runs" false
    | some r =>
      let mm := match f.case.step with | .capped _ m => m | .adaptive _ m _ => m | .fixed _ => 1.0
      let net := r.ev.net
      let nv := net.nv
      let mut worst := fzero
      let mut nlim := 0
      for rv in r.ev.views do
        for j in [0:rv.toCopy.size] do
          let (g, _) := geoOfCopy net rv.toCopy[j]!
          let h := if g < nv + net.nP then r.hP.get! (g - nv)
            else
              let t := net.triples[g - nv - net.nP]!
              min (r.hP.get! t.p0) (min (r.hP.get! t.p1) (r.hP.get! t.p2))
          worst := max worst (r.tau * ((r.rAfter[rv.X]!).get! j).abs / (mm * h))
          if (r.factor[rv.X]!).get! j < 1 then nlim := nlim + 1
      check tally s!"{f.case.name}: local limit bounds the motion of every boundary copy by max_motion h_c ({nlim} limited)"
        (worst ≤ 1 + 1e-12) s!"largest motion / (max_motion h_c) = {sci worst}; omega {sci r.omega} max {sci r.omegaMax}"
  -- the paper-stage fixtures: the limit acts at the vertices
  for f in robustFixtures.filter (·.case.stages == .paper) do
    match runOkR f.case with
    | none => check tally s!"{f.case.name}: step runs" false
    | some r =>
      let mm := match f.case.step with | .capped _ m => m | .adaptive _ m _ => m | .fixed _ => 1.0
      let net := r.ev.net
      let mut nlim := 0
      if let some (_, _, fv) := r.vlimit then
        for i in [0:net.nv] do
          if net.adj[i]! && fv.get! i < 1 then nlim := nlim + 1
      match r.vlimit with
      | none => check tally s!"{f.case.name}: vertex limit present" false
      | some (hv, _, fv) =>
        let mut worst := fzero
        for i in [0:net.nv] do
          let a := net.region[i]!
          let mut bb := if a == 0 then 1 else 0
          for k in [0:f.case.phis.size] do
            if k != a && (r.phis[k]!).get! i < (r.phis[bb]!).get! i then bb := k
          let delta := vertexChange ((r.phiLen[bb]!).get! i) ((r.phis[bb]!).get! i) ((r.phiLen[a]!).get! i) ((r.phis[a]!).get! i)
          worst := max worst (fv.get! i * delta.abs / (mm * hv.get! i))
        check tally s!"{f.case.name}: D26 vertex limit bounds the half difference of the two smallest functions by max_motion h_i ({nlim} cut-adjacent vertices limited)"
          (worst ≤ 1 + 1e-12) s!"largest change / (max_motion h_i) = {sci worst}"
  -- D26 for N regions (paper's stages): change of the whole step bounded at every vertex; energy
  let three : Array (Point Float) := #[⟨1.0, 0.1, 0.2⟩, ⟨-0.6, 0.9, -0.1⟩, ⟨-0.5, -0.8, 0.3⟩]
  let mid : Point Float := ⟨0.4, 1.0, 0.1⟩
  let spike := spikeAt ico2 (nearestVertex ico2 mid) 2.4
  let st3 := redistancedState spike (voronoiState spike three)
  let mut worstV := fzero
  let mut anyLim := false
  let mut okE := true
  for al in [0.5, 5.0, 20.0, 40.0] do
    for (step, mm) in [(StepRule.capped 0.5 0.25, 0.25), (StepRule.adaptive 0.05 0.25 1.0, 0.25)] do
      let c := { mkRegionCase "d26n" spike st3 with format := 2, stats := .robust, alphaL := al, step }
      match runOkR c with
      | none => okE := false
      | some r =>
        match r.vlimit with
        | none => okE := false
        | some (hv, _, fv) =>
          let net := r.ev.net
          for i in [0:net.nv] do
            let a := net.region[i]!
            let mut bb := if a == 0 then 1 else 0
            for k in [0:3] do
              if k != a && (r.phis[k]!).get! i < (r.phis[bb]!).get! i then bb := k
            let delta := vertexChange ((r.phiLen[bb]!).get! i) ((r.phis[bb]!).get! i) ((r.phiLen[a]!).get! i) ((r.phis[a]!).get! i)
            worstV := max worstV (fv.get! i * delta.abs / (mm * hv.get! i))
            if fv.get! i < 1 then anyLim := true
          if (al == 0.5 || al == 40.0) && !(r.evNew.total ≤ r.ev.total) then okE := false
  check tally "D26, 3 regions (paper's stages, spike): the half difference of the two smallest functions changes by at most max_motion h_i, for the paper fixtures too; the energy does not rise for alpha_L = 0.5 and 40 (it does for 5 and 20, with or without the limit: the junctions of a Voronoi start are off 120 degrees)"
    (worstV ≤ 1 + 1e-12 && okE) s!"largest change / (max_motion h_i) = {sci worstV}; vertices limited somewhere: {anyLim}"
  -- a region whose boundary is a closed chain through triple-point copies only
  for n in ["triple-chain-paper", "triple-chain-default"] do
    let path : System.FilePath := "fixtures/regions-robust" / s!"{n}.in"
    if !(← path.pathExists) then
      check tally s!"{n}: input file present" false
    else
      let text ← IO.FS.readFile path
      let lines := (text.splitOn "\n").filterMap fun l =>
        let t := tokens l
        if t.size == 0 || t[0]! == "case" || t[0]! == "end" then none else some t
      match parseBody n lines.toArray with
      | .error e => check tally s!"{n}: parses" false e
      | .ok c =>
        match runOkR c with
        | none => check tally s!"{n}: step runs" false
        | some r =>
          let net := r.ev.net
          -- regions every segment of whose boundary ends in a triple point: the chain alternates
          -- cut-point and triple-point copies and no two cut points are joined
          let mut tripleOnly := 0
          for X in [0:c.phis.size] do
            let segsX := net.segs.filter fun sg => sg.regA == X || sg.regB == X
            if segsX.size > 0 && segsX.all (·.toTriple) then tripleOnly := tripleOnly + 1
          let mm := match c.step with | .capped _ m => m | .adaptive _ m _ => m | .fixed _ => 1.0
          let mut worst := fzero
          for rv in r.ev.views do
            for j in [0:rv.toCopy.size] do
              let (g, _) := geoOfCopy net rv.toCopy[j]!
              let h := if g < net.nv + net.nP then r.hP.get! (g - net.nv)
                else
                  let t := net.triples[g - net.nv - net.nP]!
                  min (r.hP.get! t.p0) (min (r.hP.get! t.p1) (r.hP.get! t.p2))
              worst := max worst (r.tau * ((r.rAfter[rv.X]!).get! j).abs / (mm * h))
          check tally s!"{n}: {tripleOnly} region(s) whose boundary chain consists of segments to triple points only; motion of the copies bounded by max_motion h_c (default stages)"
            (tripleOnly > 0 && (c.stages == .paper || worst ≤ 1 + 1e-12) && r.tau > 0 && r.tau.isFinite)
            s!"largest motion / (max_motion h_c) = {sci worst}; tau {sci r.tau}"
  let t ← tally.get
  IO.println s!"regions-selftest (robust): {t.passed} passed, {t.failed} failed"
  return t.failed

end Ddg.Regions
