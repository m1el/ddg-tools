import DdgOracle.RepulsiveSurfaces.Saddle

/-!
# One dense projected-gradient step (SPEC §7.1 steps 2-5, §7.2, §6.2)

`g = dℰ̂`, `x` the nested-saddle direction, backtracking (halving) with the standard Armijo rule
`F(f − τx) ≤ F(f) − c τ ⟨g, x⟩`, then the projection: positional constraints exactly, nonlinear ones by
the Newton correction `h` of the step's own `Z`, `S` (old Jacobian), repeated while a row exceeds
`1e-8 |target|` (or a fixed count). No remeshing, no CCD, no penalties.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

/-- Start rule of the line search (SPEC §7.2). -/
inductive Start where
  /-- `τ₀` given -/
  | abs (tau : Float)
  /-- `τ₀` moves the farthest vertex by `lengths · ℓ̄₀` -/
  | maxdisp (lengths : Float)
  /-- `τ₀ = 2/‖x‖_F` -/
  | reference

/-- Options of one step. -/
structure StepOpts where
  /-- Armijo constant `c` -/
  armijo : Float := 1.0e-2
  /-- start rule -/
  start : Start := .maxdisp 0.1
  /-- maximal number of halvings -/
  maxHalv : Nat := 60
  /-- Newton iterations: maximal count -/
  newton : Nat := 5
  /-- `true`: exactly `newton` iterations (the reference's `Fixed`), else until tolerance -/
  newtonFixed : Bool := false

/-- One line-search trial. -/
structure Trial where
  /-- step size -/
  tau : Float
  /-- `F(f − τx)` (NaN if not finite) -/
  energy : Float
  /-- `(F0 − c τ slope − F)/|F0|`, positive iff Armijo holds -/
  margin : Float
  /-- `accepted`, `armijo-failed` or `non-finite` -/
  outcome : String
  deriving Inhabited

/-- Result of a step. -/
structure StepResult where
  /-- `F(f)` -/
  f0 : Float
  /-- `⟨g, x⟩` -/
  slope : Float
  /-- `‖x‖_F` (Frobenius over the `3|V|` entries) -/
  xnorm : Float
  /-- mean edge length `ℓ̄₀` -/
  lbar : Float
  /-- first step size -/
  tau0 : Float
  /-- the trials -/
  trials : Array Trial
  /-- accepted? (else stuck) -/
  accepted : Bool
  /-- accepted `τ` (0 when stuck) -/
  tau : Float
  /-- number of halvings -/
  halvings : Nat
  /-- `F` after the line search (before projection) -/
  e1 : Float
  /-- Newton iterations performed -/
  newtonIters : Nat
  /-- final positions -/
  pos : Array V3
  /-- final energy -/
  eFinal : Float
  /-- `Φ` of every row after the line search, before the projection -/
  phiBefore : Array Float
  /-- `Φ` of every row at the final positions -/
  phi : Array Float
  /-- pivot ratio of `A′` -/
  pivot : Float

/-- Exact positional projection: translate each constrained component, reset pins. -/
def projectPositional (t : Topo) (rc : Array RCons) (pos : Array V3) : Array V3 := Id.run do
  let mut p := pos
  for c in rc do
    match c with
    | .bary comp target =>
      let b := barycenterOf t p comp
      let sh := target - b
      for v in [0:t.nv] do
        if t.vertComp[v]! == comp then p := p.modify v (· + sh)
    | .pins vs targets =>
      for k in [0:vs.size] do p := p.set! vs[k]! targets[k]!
    | _ => pure ()
  return p

/-- Targets of the nonlinear rows (for the tolerance). -/
private def nonlinearTargets (rc : Array RCons) : Array Float :=
  rc.filterMap fun
    | .area _ tg => some tg
    | .volume _ tg => some tg
    | _ => none

/-- One step from `pos0`. `NonFiniteEnergy` at the start is an error (SPEC §2.3). -/
def projectedStep (t : Topo) (ex : Ex) (mo : MetricOpts) (rc : Array RCons) (so : StepOpts)
    (pos0 : Array V3) : Except String StepResult := do
  checkTranslation t rc
  let (f0, grad) ← energyGrad t.faces pos0 ex
  let ops ← buildOps t pos0 ex mo
  let rows := evalAll t pos0 rc
  let nst ← buildNested t.nv ops.a rows
  let gf := flatten grad
  let x := nst.direction gf
  let slope := dotF gf x
  let xnorm := Float.sqrt (dotF x x)
  let mut xinf := 0.0
  for i in [0:x.size] do xinf := max xinf (x.get! i).abs
  let lbar := meanEdge t pos0
  let tau0 := match so.start with
    | .abs tau => tau
    | .maxdisp l => Id.run do
      let mut dmax := 0.0
      for v in [0:t.nv] do
        let d := Float.sqrt (x.get! (3 * v) * x.get! (3 * v) + x.get! (3 * v + 1) * x.get! (3 * v + 1) + x.get! (3 * v + 2) * x.get! (3 * v + 2))
        dmax := max dmax d
      return l * lbar / dmax
    | .reference => 2.0 / xnorm
  let tauMin := 1.0e-12 * lbar / xinf
  let mut trials : Array Trial := #[]
  let mut accepted := false
  let mut tau := tau0
  let mut halv := 0
  let mut e1 := f0
  let mut p1 := pos0
  if slope > 0.0 && tau0.isFinite && tau0 > 0.0 then
    let mut go := true
    while go do
      if tau < tauMin || halv > so.maxHalv then
        go := false
      else
        let trial := (Array.range t.nv).map fun v =>
          pos0[v]! - tau * (⟨x.get! (3 * v), x.get! (3 * v + 1), x.get! (3 * v + 2)⟩ : V3)
        let ft := match energy t.faces (faceGeo t.faces trial) ex with
          | .ok e => e
          | .error _ => 0.0 / 0.0
        let margin := (f0 - so.armijo * tau * slope - ft) / f0.abs
        if !ft.isFinite then
          trials := trials.push { tau, energy := ft, margin := 0.0 / 0.0, outcome := "non-finite" }
          tau := tau / 2.0
          halv := halv + 1
        else if ft ≤ f0 - so.armijo * tau * slope then
          trials := trials.push { tau, energy := ft, margin, outcome := "accepted" }
          accepted := true
          e1 := ft
          p1 := trial
          go := false
        else
          trials := trials.push { tau, energy := ft, margin, outcome := "armijo-failed" }
          tau := tau / 2.0
          halv := halv + 1
  -- projection (runs also when the step is stuck)
  let phiBefore := (evalAll t p1 rc).map (·.2.2.phi)
  let mut p := projectPositional t rc p1
  let tg := nonlinearTargets rc
  let mut iters := 0
  let mut done := false
  while !done && iters < so.newton do
    let nr := (evalAll t p rc).filter (fun r => !r.2.1)
    let phi := FloatArray.ofFn' nr.size fun i => nr[i]!.2.2.phi
    let mut ok := true
    if !so.newtonFixed then
      for i in [0:nr.size] do
        let tol := 1.0e-8 * (if tg[i]! != 0.0 then tg[i]!.abs else 1.0)
        if (phi.get! i).abs > tol then ok := false
    if nr.size == 0 || (ok && !so.newtonFixed) then
      done := true
    else
      let h := nst.correction phi
      p := (Array.range t.nv).map fun v =>
        p[v]! - (⟨h.get! (3 * v), h.get! (3 * v + 1), h.get! (3 * v + 2)⟩ : V3)
      p := projectPositional t rc p
      iters := iters + 1
  let newtonDone := iters
  let eFinal := match energy t.faces (faceGeo t.faces p) ex with | .ok e => e | .error _ => 0.0 / 0.0
  let phiAll := ((evalAll t p rc).map (·.2.2.phi))
  return { f0, slope, xnorm, lbar, tau0, trials, accepted, tau := if accepted then tau else 0.0, halvings := halv,
           e1, newtonIters := newtonDone, pos := p, eFinal, phiBefore, phi := phiAll, pivot := nst.ap.pivot }

end Ddg.RepulsiveSurfaces
