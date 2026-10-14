import DdgOracle.Cuts.FixturesRobust

/-!
# Self-tests of the robust step statistics (SPEC §1, §3.9, §4, D22, D23)
-/

namespace Ddg.Cuts

/-- Brute-force weighted quantile: among the entries `(x_i, i)` ordered by value then index, the
smallest `x_i` whose cumulative weight (all entries not after it in that order) reaches `q · Σ ℓ`. -/
def bruteQuantile (q : Float) (xs ls : Array Float) : Float := Id.run do
  let n := xs.size
  if n == 0 then return 0
  let before := fun (j i : Nat) => xs[j]! < xs[i]! || (xs[j]! == xs[i]! && j ≤ i)
  let mut total := fzero
  for j in [0:n] do total := total + ls[j]!
  let mut best := posInf
  let mut bestIdx := n
  for i in [0:n] do
    let mut cum := fzero
    for j in [0:n] do
      if before j i then cum := cum + ls[j]!
    if cum ≥ q * total then
      if xs[i]! < best || (xs[i]! == best && i < bestIdx) then
        best := xs[i]!
        bestIdx := i
  return best

/-- Run the checks; returns the number of failures. -/
def selftestRobust : IO Nat := do
  let tally ← IO.mkRef ({} : Tally)
  -- quantile against the brute force (dyadic weights: the sums are exact)
  let mut rng : Rng := ⟨4242⟩
  let mut worst := fzero
  let mut worstMax := fzero
  let mut n := 0
  for trial in [0:300] do
    let (cnt, r0) := rng.below 12
    rng := r0
    let m := cnt + 1
    let mut xs : Array Float := #[]
    let mut ls : Array Float := #[]
    for _ in [0:m] do
      let (x, r1) := rng.below 6      -- many ties
      let (l, r2) := rng.below 9
      rng := r2
      let _ := r1
      xs := xs.push (Float.ofNat x - 2.5)
      ls := ls.push (Float.ofNat l / 8 + 0.125)
    let q := Float.ofNat (trial % 11) / 10
    let a := weightedQuantile q xs.toList ls.toList
    let b := bruteQuantile q xs ls
    worst := max worst (a - b).abs
    let mx := xs.foldl (fun acc x => max acc x) negInf
    worstMax := max worstMax ((weightedQuantile 1.0 xs.toList ls.toList) - mx).abs
    n := n + 1
  check tally s!"weighted quantile against the brute-force definition ({n} random lists with ties); q = 1 is the maximum"
    (worst == 0 && worstMax == 0) s!"max difference {sci worst}; q = 1 vs maximum {sci worstMax}"
  check tally "weighted quantile: empty list gives 0, one value gives it"
    (weightedQuantile 0.9 [] [] == 0 && weightedQuantile 0.9 [3.5] [2.0] == 3.5)
  -- uniform cut: the robust step equals the step of the extremes
  let tauR := robustStep 1.0 0.05 100.0 0.25 (weightedQuantile 0.9 [40.0, 40.0, 40.0] [1.0, 2.0, 3.0])
  let tauE := adaptiveStep 1.0 0.05 100.0 0.25 0.1 4.0
  check tally "uniform cut: robust and extremes give the same adaptive step (speed 40 in edges of 0.1, omega 100)"
    ((tauR - tauE).abs ≤ 1e-15) s!"{sci tauR} vs {sci tauE}"
  -- equator of a regular icosphere: nearly uniform, so the two agree to a few percent
  let ico3 := icosphere 3
  let eq := phiOf ico3 fun p => Float.acos p.z - floatPi / 2
  let base := { mkCase "eq" ico3 eq with alphaL := 2.0, format := 3, step := .adaptive 0.05 0.25 1.0, flux := .smoothed 0.5 }
  match runOk { base with stats := .robust }, runOk { base with stats := .extremes } with
  | some a, some b =>
    check tally "equator of icosphere 3 (nearly uniform cut): robust and extremes steps agree within 5 percent, robust >= extremes"
      (a.tau ≥ b.tau && a.tau ≤ 1.05 * b.tau) s!"tau robust {sci a.tau} extremes {sci b.tau}"
  | _, _ => check tally "equator: steps run" false
  -- local limit: every point moves at most max_motion * h_c (two regions)
  let ico2 := icosphere 2
  let kf := FloatArray.ofFn' ico2.faces.size fun fi =>
    let f := ico2.faces[fi]!
    let hot := vertexAtPolarAngle ico2 1.0
    if f.a == hot || f.b == hot || f.c == hot then 40.0 * 4 * floatPi else 4 * floatPi
  let capPhi := phiOf ico2 fun p => Float.acos p.z - 1.0
  for (nm, step, mm) in [("capped", StepRule.capped 0.01 0.25, 0.25), ("adaptive", StepRule.adaptive 0.05 0.25 1.0, 0.25)] do
    match runOk { withRobust (mkCase "lim" ico2 capPhi) with curv := .perFace kf, alphaL := 2.0, step } with
    | none => check tally s!"local limit ({nm}): step runs" false
    | some r =>
      let st := r.stage2.get!
      let mut worstMove := fzero
      let mut nlim := 0
      for k in [0:r.ev.cm.nc] do
        let move := r.tau * (st.rlim.get! k).abs / (mm * st.hc.get! k)
        worstMove := max worstMove move
        if st.factor.get! k < 1 then nlim := nlim + 1
      check tally s!"local limit ({nm}, curvature spike): motion of every cut point <= max_motion h_c; {nlim} points limited"
        (worstMove ≤ 1 + 1e-12 && nlim > 0) s!"largest motion / (max_motion h_c) = {sci worstMove}"
  -- D26: robust statistics and the length step on the mesh: the limit acts at the vertices on the whole step
  let mut worstV := fzero
  let mut anyLimited := false
  let mut energyOk := true
  let mut worstDelta := fzero
  for al in [0.5, 5.0, 20.0, 40.0] do
    for (step, mm) in [(StepRule.capped 0.5 0.25, 0.25), (StepRule.adaptive 0.05 0.25 1.0, 0.25)] do
      let c := { mkCase "d26" ico2 capPhi with format := 3, stats := .robust, curv := .perFace kf, alphaL := al, step }
      match runOk c with
      | none => energyOk := false
      | some r =>
        match r.vlimit with
        | none => energyOk := false
        | some (hv, _, fv, pl) =>
          for i in [0:r.base.nv] do
            let ch := (pl.get! i - r.phi.get! i).abs
            worstV := max worstV (ch / (mm * hv.get! i))
            if fv.get! i < 1 then anyLimited := true
          if !(r.evNew.total ≤ r.ev.total) then energyOk := false
          worstDelta := max worstDelta (r.ev.total - r.evNew.total)
  check tally "D26 (robust, length step on the mesh, curvature spike, alpha_L 0.5 .. 40, capped and adaptive): no vertex changes by more than max_motion h_i; vertices are limited; the energy after one step is not above the energy before"
    (worstV ≤ 1 + 1e-12 && anyLimited && energyOk) s!"largest change / (max_motion h_i) = {sci worstV}"
  -- with the limit not active the step equals the unlimited one
  let free := { mkCase "free" ico2 capPhi with format := 3, stats := .robust, curv := .perFace kf, alphaL := 2.0 }
  match runOk { free with step := .capped 0.0005 1.0e6 }, runOk { free with step := .fixed 0.0005 } with
  | some a, some b =>
    let mut d := fzero
    for i in [0:a.base.nv] do d := max d ((a.phiNew.get! i - b.phiNew.get! i).abs)
    check tally "limit not active (huge max_motion): the robust step equals the step with a fixed tau" (d ≤ 1e-15 && a.tau == b.tau) s!"difference {sci d}"
  | _, _ => check tally "limit not active: steps run" false
  -- zero slope: a face with |grad phi| = 0 is left out of the mean (D23)
  let ico1 := icosphere 1
  let zero := withDefaultStages (mkCase "zs" (ico1.mapPos fun p => ⟨p.x * (1 + 0.05 * Float.sin (5 * p.y + p.z)), p.y * (1 + 0.04 * Float.cos (4 * p.x - 2 * p.z)), p.z⟩) (phiOf ico1 fun p => p.z))
  match runOk { zero with alphaL := 1.0 } with
  | none => check tally "zero slope: step runs" false
  | some r =>
    let st := r.stage2.get!
    let mut allPos := true
    for k in [0:r.ev.cm.nc] do
      if !(st.gbar.get! k > 0 && (st.gbar.get! k).isFinite) then allPos := false
    -- count faces with a zero slope
    let mut nz := 0
    for sg in r.ev.cm.segs do
      let t := r.ev.cm.tris[sg.tri]!
      if !(gradientNormLone (r.phi.get! sg.lone) sg.len t.area > 0) then nz := nz + 1
    check tally s!"zero slope (phi = z, exactly 0 at vertices): {nz} crossed faces with |grad phi| = 0 are left out, gbar is positive and finite"
      (nz > 0 && allPos) s!"gbar in {sci (st.gbar.foldl min posInf)} .. {sci (st.gbar.foldl max negInf)}"
  -- D30: a face with a moved crossing does not count in the mean slope
  let oldMean := fun (cm : CutMesh) (phi : FloatArray) (countMoved : Bool) => Id.run do
    -- one piece assumed per call: sums over all segments (the closed circle / the single cut)
    let mut sL := fzero
    let mut sLG := fzero
    for sg in cm.segs do
      let t := cm.tris[sg.tri]!
      let g := gradientNormLone (phi.get! sg.lone) sg.len t.area
      let moved := cm.crossings[sg.p]!.moved || cm.crossings[sg.q]!.moved
      if g > 0 && (countMoved || !moved) then
        sL := sL + sg.len
        sLG := sLG + sg.len / g
    return if sL > 0 then sL / sLG else 1.0
  for f in robustFixtures do
    if f.case.name == "rb-moved-some" || f.case.name == "rb-moved-all" then
      match runOk f.case with
      | none => check tally s!"{f.case.name}: step runs" false
      | some r =>
        let st := r.stage2.get!
        let cm := r.ev.cm
        let mut nMoved := 0
        for sg in cm.segs do
          if cm.crossings[sg.p]!.moved || cm.crossings[sg.q]!.moved then nMoved := nMoved + 1
        let expected := oldMean cm r.phi false
        let withMoved := oldMean cm r.phi true
        let mut maxDiff := fzero
        for k in [0:cm.nc] do maxDiff := max maxDiff ((st.gbar.get! k - expected).abs)
        if f.case.name == "rb-moved-some" then
          check tally s!"{f.case.name}: {nMoved} of {cm.segs.size} faces have a moved crossing; gbar of the piece = mean over the other faces ({sci expected}; counting them would give {sci withMoved})"
            (nMoved > 0 && nMoved < cm.segs.size && maxDiff ≤ 1e-15 && withMoved < 0.5 * expected) s!"difference {sci maxDiff}"
        else
          check tally s!"{f.case.name}: every one of the {cm.segs.size} faces has a moved crossing; gbar = 1"
            (nMoved == cm.segs.size && maxDiff == 0 && expected == 1.0) s!"gbar {sci (st.gbar.get! 0)}"
  -- a state without moved crossings: gbar equals the mean over all faces with a positive slope, to the last bit
  let mut sameBits := true
  for (nm, m, phi) in [("cap 1.0", ico2, capPhi),
                       ("cap 0.7", ico2, phiOf ico2 fun p => Float.acos p.z - 0.7),
                       ("cap 1.3", ico2, phiOf ico2 fun p => Float.acos p.z - 1.3)] do
    match runOk (withRobust (mkCase nm m phi)) with
    | none => sameBits := false
    | some r =>
      let cm := r.ev.cm
      let mut anyMoved := false
      for cr in cm.crossings do if cr.moved then anyMoved := true
      let st := r.stage2.get!
      let a := (st.gbar.get! 0)
      let b := oldMean cm r.phi true
      if anyMoved then sameBits := false
      if a.toBits != b.toBits then sameBits := false
  check tally "no moved crossing: gbar is unchanged to the last bit (three caps on icosphere 2)" sameBits

  let t ← tally.get
  IO.println s!"cuts-selftest (robust): {t.passed} passed, {t.failed} failed"
  return t.failed

end Ddg.Cuts
