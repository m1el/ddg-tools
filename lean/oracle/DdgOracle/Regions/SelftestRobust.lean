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
  -- local limit on networks: motion of every boundary copy <= max_motion h_c
  for f in robustFixtures do
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
  let t ← tally.get
  IO.println s!"regions-selftest (robust): {t.passed} passed, {t.failed} failed"
  return t.failed

end Ddg.Regions
