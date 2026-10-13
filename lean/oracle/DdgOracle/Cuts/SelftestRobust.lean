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
  let t ← tally.get
  IO.println s!"cuts-selftest (robust): {t.passed} passed, {t.failed} failed"
  return t.failed

end Ddg.Cuts
