import DdgOracle.Selftest
import DdgOracle.Developability.Gen

/-!
# Internal checks of the developability oracle

Closed forms of SPEC 10.3, brute-force searches that do not use the candidate enumerations,
finite differences of the oracle's own gradient away from ties, invariances, subdivision, and the
line search.
-/

namespace Ddg.Developability

open Ddg

private def say (s : String) : IO Unit := do
  IO.println s
  (← IO.getStdout).flush

private structure Tally where
  failures : IO.Ref Nat
  checks : IO.Ref Nat

private def Tally.check (t : Tally) (ok : Bool) (what : String) : IO Unit := do
  t.checks.modify (· + 1)
  if !ok then
    t.failures.modify (· + 1)
    say s!"FAIL  {what}"


private def relErr (a b : Float) : Float := (a - b).abs / max (max a.abs b.abs) 1e-300

/-- Energy values at vertex `v` for all six energies. -/
private def energiesAt (m : DMesh) (topo : Topo) (v : Nat) : Array Float :=
  EnergyKind.all.map fun e =>
    match evaluate m topo { energy := e, omitValenceThree := false } with
    | .ok ev => ev.e[v]!
    | .error _ => nan

/-- Brute-force `E_P` and `π_max` of a star: all subsets, connectedness tested directly. -/
private def bruteCombinatorial (ns : Array V3) (closed : Bool) : Float × Float := Id.run do
  let k := ns.size
  let mut bestP := posInf
  let mut bestM := posInf
  for mask in [1:2 ^ k - 1] do
    let inS := fun (i : Nat) => (mask >>> i) % 2 == 1
    let runs := (Array.range k).foldl (fun acc i =>
      let nxt := if closed then (i + 1) % k else i + 1
      let cur := inS i
      if nxt < k then (if cur && !inS nxt then acc + 1 else acc)
      else (if cur then acc + 1 else acc)) 0
    -- for a path: number of runs of S (a run ends where S stops or at the last index);
    -- both S and its complement must be a single run
    let runsC := (Array.range k).foldl (fun acc i =>
      let nxt := if closed then (i + 1) % k else i + 1
      let cur := !inS i
      if nxt < k then (if cur && inS nxt then acc + 1 else acc)
      else (if cur then acc + 1 else acc)) 0
    if runs == 1 && runsC == 1 then
      let f1 := (Array.range k).filter inS
      let f2 := (Array.range k).filter (!inS ·)
      let width := fun (f : Array Nat) =>
        let mean := V3.smul (1 / f.size.toFloat) (f.foldl (fun a i => a + ns[i]!) V3.zero)
        (f.foldl (fun a i => a + (ns[i]! - mean).norm2) 0) / f.size.toFloat
      let wmax := fun (f : Array Nat) => Id.run do
        let mut b := 0.0
        for a in f do
          for c in f do
            b := max b (ns[a]! - ns[c]!).norm2
        return b
      bestP := min bestP (width f1 + width f2)
      bestM := min bestM (max (wmax f1) (wmax f2))
  if bestP == posInf then return (0, 0)
  return (bestP, bestM)

/-- Dense search for `min_u max_σ ⟨u, N_σ⟩²` over the sphere (upper bound on the true minimum). -/
private def denseCovMax (ns : Array V3) (res : Nat) : Float := Id.run do
  let mut best := posInf
  for i in [0:res + 1] do
    let th := 3.141592653589793 * i.toFloat / res.toFloat
    for j in [0:2 * res] do
      let ph := 3.141592653589793 * j.toFloat / res.toFloat
      let u : V3 := ⟨Float.sin th * Float.cos ph, Float.sin th * Float.sin ph, Float.cos th⟩
      let phi := ns.foldl (fun a n => max a ((V3.dot u n) ^ 2)) 0
      best := min best phi
  return best

/-- Total energy of `m` (domain: everything). -/
private def totalOf (m : DMesh) (topo : Topo) (e : EnergyKind) : Float :=
  match evaluate m topo { energy := e, omitValenceThree := false } with
  | .ok ev => ev.total
  | .error _ => nan

/-- Run all checks; returns the number of failures. -/
def selftest : IO Nat := do
  let t : Tally := { failures := ← IO.mkRef 0, checks := ← IO.mkRef 0 }
  -- 1. closed forms on regular cones
  let pi := 3.141592653589793
  let mut worst := #[0.0, 0.0, 0.0, 0.0, 0.0]
  for k in [3:10] do
    for h in [0.1, 0.6, 2.0] do
      let m := regularCone k h
      match buildTopo m with
      | .error e => t.check false s!"cone {k} {h}: {e}"
      | .ok topo =>
        let es := energiesAt m topo 0
        let kf := k.toFloat
        let c := Float.cos (pi / kf)
        let s2 := h * h / (h * h + c * c)
        let c2 := 1 - s2
        let theta := 2 * Float.asin (Float.sin (pi / kf) / Float.sqrt (1 + h * h))
        let alpha := Float.atan2 h c
        let eLam := theta * kf * min (s2 / 2) c2
        let rho := fun (mm : Nat) => Float.sin (mm.toFloat * pi / kf) / (mm.toFloat * Float.sin (pi / kf))
        let mut eP := posInf
        for mm in [1:k] do
          eP := min eP (2 - rho mm ^ 2 - rho (k - mm) ^ 2)
        eP := s2 * eP
        let lmax := min c2 (s2 * c * c)
        let eInt := theta * kf * alpha * alpha / 2
        let eIntMax := if k % 2 == 0 then (2 * alpha * c) ^ 2 else (alpha * (1 + c)) ^ 2
        let errs := #[relErr es[2]! eLam, relErr es[0]! eP, if k % 2 == 0 then relErr es[3]! lmax else 0.0,
                      relErr es[4]! eInt, relErr es[5]! eIntMax]
        for i in [0:5] do worst := worst.set! i (max worst[i]! errs[i]!)
        t.check (errs.all (· < 1e-9)) s!"cone k={k} h={h}: relative errors E_lambda {sci errs[0]!}, E_P {sci errs[1]!}, lmax {sci errs[2]!}, intrinsic {sci errs[3]!}, intrinsic max {sci errs[4]!}"
  say s!"cones: worst relative error E_lambda {sci worst[0]!}, E_P {sci worst[1]!}, lambda_max {sci worst[2]!}, intrinsic {sci worst[3]!}, intrinsic max {sci worst[4]!}"
  -- 2. exact hinges and flat stars vanish; valence-3 flag; two-face boundary stars
  let mut hingeMax := 0.0
  for seed in [0:20] do
    let (m, _) := (hingeStar (seed % 4 + 1) (seed % 3 + 1) 0.7 0.0).run ⟨seed.toUInt64 + 1⟩
    match buildTopo m with
    | .error e => t.check false s!"hinge: {e}"
    | .ok topo =>
      for x in energiesAt m topo 0 do hingeMax := max hingeMax x.abs
  say s!"exact hinge stars: largest |energy| {sci hingeMax}"
  t.check (hingeMax < 1e-14) s!"hinge energies vanish ({sci hingeMax})"
  -- 3. brute force: E_P, pi_max, lambda_max, dense search
  let mut maxEP := 0.0
  let mut maxPM := 0.0
  let mut maxLM := 0.0
  let mut worstDense := 0.0
  let mut nstars := 0
  for i in [0:60] do
    let c := genCase 7 i
    match buildTopo c.mesh with
    | .error _ => pure ()
    | .ok topo =>
      match evaluate c.mesh topo { energy := .combinatorial, omitValenceThree := false },
            evaluate c.mesh topo { energy := .combinatorialMax, omitValenceThree := false },
            evaluate c.mesh topo { energy := .covarianceMax, omitValenceThree := false } with
      | .ok eP, .ok eM, .ok eL =>
        for v in [0:topo.nv] do
          let st := topo.stars[v]!
          if st.faces.size ≤ 10 then
            let ns := st.faces.map fun f => (faceData c.mesh f).get!.N
            let (bp, bm) := bruteCombinatorial ns st.closed
            maxEP := max maxEP (relErr eP.e[v]! bp)
            maxPM := max maxPM (relErr eM.e[v]! bm)
            let dense := denseCovMax ns 90
            -- the candidate minimum must not exceed the dense upper bound and must be close to it
            if eL.e[v]! > dense + 1e-12 then maxLM := max maxLM (eL.e[v]! - dense)
            worstDense := max worstDense (dense - eL.e[v]!)
            nstars := nstars + 1
      | _, _, _ => pure ()
  say s!"brute force over {nstars} stars: E_P {sci maxEP}, pi_max {sci maxPM}, lambda_max above dense search by {sci maxLM}; dense search above by up to {sci worstDense}"
  t.check (maxEP < 1e-12 && maxPM < 1e-12) "combinatorial energies equal brute force over all subsets"
  t.check (maxLM == 0) "lambda_max never exceeds the dense search"
  t.check (worstDense < 5e-2) "lambda_max is not far below the dense search (sanity)"
  -- 4. finite differences of the gradient away from ties
  let mut fdWorst := EnergyKind.all.map fun _ => 0.0
  let mut fdCount := EnergyKind.all.map fun _ => 0
  for i in [0:80] do
    let c := genCase 11 i
    match buildTopo c.mesh with
    | .error _ => pure ()
    | .ok topo =>
      for ei in [0:6] do
        let e := EnergyKind.all[ei]!
        let o : Opts := { energy := e, omitValenceThree := false }
        match evaluate c.mesh topo o with
        | .error _ => pure ()
        | .ok ev0 =>
          let h := 1e-6
          for p in [0:min topo.nv 6] do
            -- stars touching p: p and its neighbours
            let touching := (Array.range topo.nv).filter fun u =>
              u == p || topo.edges.any fun (a, b) => (a == u && b == p) || (a == p && b == u)
            let okMargin := touching.all fun u => ev0.margin[u]! > 1e-5 * max 1 ev0.e[u]!.abs
            if okMargin then
              for coord in [0:3] do
                let shift := fun (s : Float) => ({ c.mesh with pos := c.mesh.pos.modify p fun q =>
                  if coord == 0 then ⟨q.x + s, q.y, q.z⟩ else if coord == 1 then ⟨q.x, q.y + s, q.z⟩ else ⟨q.x, q.y, q.z + s⟩ } : DMesh)
                match evaluate (shift h) topo o, evaluate (shift (-h)) topo o with
                | .ok a, .ok b =>
                  let same := touching.all fun u => a.choice[u]! == ev0.choice[u]! && b.choice[u]! == ev0.choice[u]!
                  if same then
                    let fd := (a.total - b.total) / (2 * h)
                    let g := ev0.grad[p]!.get coord
                    let err := (fd - g).abs / max (ev0.grad[p]!.maxAbs) 1e-2
                    fdWorst := fdWorst.modify ei (max · err)
                    fdCount := fdCount.modify ei (· + 1)
                | _, _ => pure ()
  for ei in [0:6] do
    say s!"finite differences, {EnergyKind.all[ei]!.name}: {fdCount[ei]!} components, worst scaled error {sci fdWorst[ei]!}"
    t.check (fdCount[ei]! > 100) s!"enough finite-difference samples for {EnergyKind.all[ei]!.name}"
    t.check (fdWorst[ei]! < 1e-5) s!"gradient of {EnergyKind.all[ei]!.name} matches finite differences"
  -- 5. invariances: translation, rotation (about z), scaling; sum of gradients; moment
  let mut invWorst := 0.0
  let mut sumWorst := 0.0
  for i in [0:40] do
    let c := genCase 5 i
    match buildTopo c.mesh with
    | .error _ => pure ()
    | .ok topo =>
      for e in EnergyKind.all do
        let o : Opts := { energy := e, omitValenceThree := false }
        match evaluate c.mesh topo o with
        | .error _ => pure ()
        | .ok ev0 =>
          let cs := Float.cos 0.7
          let sn := Float.sin 0.7
          let rigid := fun (s : Float) (q : V3) => V3.smul s (⟨cs * q.x - sn * q.y + 3, sn * q.x + cs * q.y - 2, q.z + 5⟩ : V3)
          let moved := fun (s : Float) => ({ pos := c.mesh.pos.map (rigid s), faces := c.mesh.faces } : DMesh)
          match evaluate (moved 1000) topo o with
          | .ok ev1 => invWorst := max invWorst ((ev0.total - ev1.total).abs / (ev0.total.abs + 1e-4))
          | .error _ => pure ()
          let gs := ev0.grad.foldl (· + ·) V3.zero
          let gmax := ev0.grad.foldl (fun a g => max a g.maxAbs) 1e-300
          let s1 := gs.maxAbs / gmax / topo.nv.toFloat
          if s1 > 1e-9 then say s!"  sum of gradient {c.name} {e.name}: {sci s1}"
          sumWorst := max sumWorst s1
          let mom := (Array.range topo.nv).foldl (fun a v => a + V3.cross c.mesh.pos[v]! ev0.grad[v]!) V3.zero
          let scl := (Array.range topo.nv).foldl (fun a v => a + (V3.cross c.mesh.pos[v]! ev0.grad[v]!).norm + c.mesh.pos[v]!.norm * ev0.grad[v]!.norm) 1e-300
          if ev0.total > 1e-3 && mom.maxAbs / scl > 1e-9 then say s!"  moment {c.name} {e.name}: {sci (mom.maxAbs / scl)}"
          if ev0.total > 1e-3 then sumWorst := max sumWorst (mom.maxAbs / scl)
          let rad := (Array.range topo.nv).foldl (fun a v => a + V3.dot c.mesh.pos[v]! ev0.grad[v]!) 0
          if rad.abs / scl > 1e-9 then say s!"  radial {c.name} {e.name}: {sci (rad.abs / scl)}"
          if ev0.total > 1e-3 then sumWorst := max sumWorst (rad.abs / scl)
  say s!"invariances: rigid motion and scaling x1000 {sci invWorst}; sum, moment, radial moment of the gradient {sci sumWorst}"
  t.check (invWorst < 1e-9) "energy invariant under rigid motion and scaling"
  t.check (sumWorst < 1e-9) "gradient has zero sum, zero moment, zero radial moment"
  -- 6. subdivision preserves every energy; counts
  let mut subWorst := 0.0
  for i in [0:30] do
    let c := genCase 3 i
    match buildTopo c.mesh with
    | .error _ => pure ()
    | .ok topo =>
      let (m2, mids) := subdivide c.mesh topo
      t.check (m2.faces.size == 4 * c.mesh.faces.size && m2.pos.size == c.mesh.pos.size + topo.edges.size && mids.size == topo.edges.size)
        s!"subdivision counts of {c.name}"
      match buildTopo m2 with
      | .error err => t.check false s!"subdivided {c.name} is not a valid mesh: {err}"
      | .ok topo2 =>
        for e in EnergyKind.all do
          let a := totalOf c.mesh topo e
          let b := totalOf m2 topo2 e
          if a.isFinite && b.isFinite then
            subWorst := max subWorst ((a - b).abs / max a.abs 1e-3)
  say s!"subdivision: largest relative change of any energy {sci subWorst}"
  t.check (subWorst < 1e-9) "4-1 subdivision preserves the energies"
  -- 7. line search: accepted steps satisfy the conditions that define them
  let mut stepChecked := 0
  let mut stepBad := 0
  for i in [0:60] do
    let c := genCase 9 i
    match buildTopo c.mesh with
    | .error _ => pure ()
    | .ok topo =>
      let o : Opts := { energy := EnergyKind.all[i % 6]!, includeBoundary := c.includeBoundary, omitValenceThree := c.omit3 }
      match step c.mesh topo o c.stepOpts with
      | .error _ => pure ()
      | .ok r =>
        if r.kind == .wolfe || r.kind == .capped || r.kind == .weak then
          stepChecked := stepChecked + 1
          let h := r.e1 - r.e0
          if !(h < c.stepOpts.c1 * r.slope * r.t || (c.stepOpts.c1 == 0 && h < 0)) then stepBad := stepBad + 1
          if !(r.t ≤ r.tmax * (1 + 1e-15)) then stepBad := stepBad + 1
  say s!"line search: {stepChecked} steps checked, {stepBad} violations"
  t.check (stepBad == 0 && stepChecked > 10) "accepted steps satisfy the Armijo condition and the cap"
  -- 8. bad input
  let tri : DMesh := { pos := #[⟨0, 0, 0⟩, ⟨1, 0, 0⟩, ⟨0, 1, 0⟩, ⟨1, 1, 0⟩], faces := #[(0, 1, 2), (1, 3, 2)] }
  t.check (buildTopo tri |>.isOk) "a two-triangle disk is accepted"
  t.check (!(buildTopo { tri with faces := #[(0, 1, 2), (0, 1, 3)] } |>.isOk)) "inconsistent orientation is rejected"
  t.check (!(buildTopo { tri with faces := #[(0, 1, 2), (1, 3, 5)] } |>.isOk)) "out-of-range vertex is rejected"
  match buildTopo { tri with pos := #[⟨0, 0, 0⟩, ⟨1, 0, 0⟩, ⟨2, 0, 0⟩, ⟨1, 1, 0⟩] } with
  | .ok topo =>
    t.check (match evaluate { tri with pos := #[⟨0, 0, 0⟩, ⟨1, 0, 0⟩, ⟨2, 0, 0⟩, ⟨1, 1, 0⟩] } topo {} with
      | .error (.degenerate _) => true | _ => false) "a degenerate face is reported"
  | .error _ => t.check false "degenerate-face mesh topology"
  -- 9. D33: boundary stars with at most two faces contribute exactly nothing
  let hingeB := (hingeStar 3 3 0.8 1e-9).run ⟨5⟩ |>.1
  match buildTopo hingeB with
  | .error _ => t.check false "hinge mesh topology"
  | .ok topo =>
    for e in EnergyKind.all do
      match evaluate hingeB topo { energy := e, includeBoundary := true }, evaluate hingeB topo { energy := e, includeBoundary := false } with
      | .ok a, .ok b =>
        let same := (Array.range topo.nv).all fun v => (a.grad[v]! - b.grad[v]!).norm == 0 && (!topo.stars[v]!.closed → a.e[v]! == 0)
        t.check (same && a.total == b.total) s!"D33: two-face boundary stars add nothing ({e.name})"
      | _, _ => t.check false s!"D33 evaluation ({e.name})"
  -- 10. SPEC 4.9: floor, blocking, start error
  let sq : DMesh := { pos := #[⟨0, 0, 0⟩, ⟨1, 0, 0⟩, ⟨1, 1, 0.3⟩, ⟨0, 1, 0⟩, ⟨0.5, 0.5, 0.6⟩], faces := #[(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)] }
  match buildTopo sq with
  | .error _ => t.check false "floor test mesh topology"
  | .ok topo =>
    let o : Opts := { energy := .covariance }
    t.check (match step sq topo o { floor := 1.2 } with | .error (.nearDegenerate _) => true | _ => false) "a start below the floor is an error"
    t.check (match step sq topo o { floor := 0.0 } with | .ok _ => true | _ => false) "floor 0 is off"
    t.check (match step sq topo o { floor := 1e-4 } with | .ok r => r.inadmissible == 0 && r.blockedCount == 0 | _ => false) "the default floor blocks nothing here"
    let minA := (Array.range 4).foldl (fun a f => min a (minAngle sq f)) 10.0
    t.check (match step sq topo o { floor := minA * 0.75 } with | .ok r => r.blockedCount == 5 | _ => false) "faces within twice the floor block their vertices"
  let failures ← t.failures.get
  let checks ← t.checks.get
  say s!"dev-selftest: {checks - failures} of {checks} checks passed"
  return failures

end Ddg.Developability
