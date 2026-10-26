import DdgOracle.RepulsiveCurves.Saddle

/-!
# One (or a few chained) normalized-Armijo steps (SPEC §7.1, §7.3)

`[formula: rc-armijo]`. Step `k`: energy and differential at the current curve, `A` and `C` there,
`g̃` from the saddle system (§6.1), then the line search of §7.3: start
`τ₀ = 1/n` if `n = ‖g̃‖₂ > 1` else `1/√n`, and if the previous accepted step size `τ_prev > 1e-5`,
`τ₀ = min(1.5 τ_prev, 4 τ₀)`; then, for `τ = τ₀, ŝτ₀, ŝ²τ₀, …`: `γ_τ = γ − τ g̃`, chord
projection (§6.2, `Ā` and `C` frozen at the current curve, tolerance `τ_len` of the session's
`ℓ̄₀`); accept when the projection converged and `Ê(γ̂_τ) ≤ Ê(γ) − c τ (dÊ·g̃)`; give up after more
than `maxHalvings` shrinks or when `τ < 1e-12 τ₀` (the step is *stuck*, the curve does not move).
Constraint targets are those of the first curve, `ℓ̄₀` is the mean edge length of the first curve.
-/

namespace Ddg.RepulsiveCurves

/-- Parameters of the line search. -/
structure SearchOpts where
  /-- Armijo constant `c` -/
  armijo : Float := 1e-4
  /-- shrink factor -/
  shrink : Float := 0.5
  /-- maximal number of shrinks -/
  maxHalvings : Nat := 40
  /-- number of chained steps -/
  nSteps : Nat := 1
  deriving Inhabited

/-- One trial of the line search. -/
structure TrialRec where
  /-- step size -/
  tau : Float
  /-- `accepted`, `armijo-failed`, `non-finite`, or `projection-failed:STATUS` -/
  outcome : String
  /-- chord iterations -/
  iters : Nat
  /-- `Ê` at the projected curve (NaN if not evaluated) -/
  energy : Float
  /-- `max|Φ|` after the projection -/
  viol : Float
  /-- `(F₀ − c τ slope − F)/|F₀|` (NaN if not evaluated): positive iff Armijo holds -/
  margin : Float
  deriving Inhabited

/-- One step. -/
structure StepRec where
  /-- first trial size -/
  tau0 : Float
  /-- `dÊ · g̃` -/
  slope : Float
  /-- `Ê` at the base curve -/
  f0 : Float
  /-- `‖g̃‖₂` -/
  gnorm : Float
  /-- all trials -/
  trials : Array TrialRec
  /-- accepted? -/
  accepted : Bool
  /-- accepted step size (0 if stuck) -/
  tau : Float
  /-- number of shrinks -/
  halvings : Nat
  /-- `Ê` after the step -/
  e1 : Float
  /-- curve after the step -/
  pos : Array V3
  /-- worst conditioning of the saddle -/
  cond : Float
  /-- smallest equilibrated pivot ratio -/
  pivot : Float
  deriving Inhabited

private def nanF : Float := Float.ofBits 0x7ff8000000000000

private def statusToken : ProjStatus → String
  | .converged => "converged"
  | .maxIterations => "max-iterations"
  | .diverged => "diverged"
  | .tangentFlipped c => s!"tangent-flipped-{c}"

/-- One normalized-Armijo step from `pos`. `lbar0` and `tol` are the session's `ℓ̄₀` and `τ_len`. -/
def normalizedStep (net : Net) (ex : Expo) (q : Quad) (coords : Array Nat) (rc : Array RCons)
    (opts : SearchOpts) (tol lbar0 : Float) (maxit : Nat) (pos : Array V3) (prevTau : Option Float) :
    Except String StepRec := do
  let cur : Net := { net with pos }
  let (f0, d0) := energyAndDifferential cur ex
  -- no component along a held axis
  let d := d0.map fun p => Id.run do
    let mut r := p
    for ax in [0:3] do
      if !coords.contains ax then r := r.set ax 0.0
    return r
  let rows := evalRows cur pos coords rc
  let G := gram cur ex q
  let sad ← Saddle.build net.nv coords G.a rows
  let gr := projectedGradient sad d
  let n := Float.sqrt (gr.g.foldl (fun acc p => acc + p.norm2) 0.0)
  let slope := gr.dfg
  let maxSpeed := gr.g.foldl (fun acc p => max acc p.norm) 0.0
  let mut res : StepRec := { tau0 := 0.0, slope, f0, gnorm := n, trials := #[], accepted := false, tau := 0.0,
                             halvings := 0, e1 := f0, pos, cond := sad.cond, pivot := sad.pivotRatio }
  if !(slope > 0.0 && slope.isFinite && maxSpeed > 0.0 && maxSpeed.isFinite) then return res
  let mut tau0 := if n > 1.0 then 1.0 / n else 1.0 / Float.sqrt n
  match prevTau with
  | some p => if p > 1e-5 then tau0 := min (1.5 * p) (4.0 * tau0)
  | none => pure ()
  res := { res with tau0 }
  if !(tau0 > 0.0 && tau0.isFinite) then return res
  let mut tau := tau0
  let mut halv := 0
  let mut trials : Array TrialRec := #[]
  let mut done := false
  while !done do
    let trial := (Array.range net.nv).map fun v => pos[v]! - V3.smul tau gr.g[v]!
    let pr := chordProject net rc sad trial tol lbar0 maxit
    let mut outcome := "armijo-failed"
    let mut ok := true
    match pr.status with
    | .converged => pure ()
    | st =>
      outcome := s!"projection-failed:{statusToken st}"
      ok := false
    let mut f := nanF
    let mut margin := nanF
    if ok then
      let (fe, _) := energyAndDifferential { net with pos := pr.pos } ex
      f := fe
      margin := (f0 - opts.armijo * tau * slope - f) / (if f0.abs > 1e-300 then f0.abs else 1.0)
      if !f.isFinite then outcome := "non-finite"
      else if f ≤ f0 - opts.armijo * tau * slope then
        outcome := "accepted"
        trials := trials.push { tau, outcome, iters := pr.iterations, energy := f, viol := pr.maxViolation, margin }
        res := { res with trials, accepted := true, tau, halvings := halv, e1 := f, pos := pr.pos }
        done := true
    if !done then
      trials := trials.push { tau, outcome, iters := pr.iterations, energy := f, viol := pr.maxViolation, margin }
      tau := tau * opts.shrink
      halv := halv + 1
      if halv > opts.maxHalvings || tau < 1e-12 * tau0 then
        res := { res with trials, halvings := halv }
        done := true
  return res

end Ddg.RepulsiveCurves
