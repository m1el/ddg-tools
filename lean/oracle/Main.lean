import DdgOracle
open Ddg

def describe (p : Params) : String := Id.run do
  let nv := p.wE.size
  let mut free := 0
  let mut forb := 0
  let mut zeroE := 0
  for v in [0:nv] do
    if p.wR.get! v == 0.0 then free := free + 1
    if p.wR.get! v == posInf then forb := forb + 1
    if p.wE.get! v == 0.0 then zeroE := zeroE + 1
  return s!"lambda={p.lambda} lo={p.lo.get! 0} hi={p.hi.get! 0} free={free} forb={forb} zeroE={zeroE} cones={p.cones.size} wR0={p.wR.get! 0}"

def report (name : String) (m : RawMesh) (p : Params) (verbose := true) : IO Bool := do
  match m.operators with
  | .error e => IO.println s!"{name}: ERROR {e}"; return false
  | .ok ops =>
    let t0 ← IO.monoMsNow
    match solve ops p with
    | .error e => IO.println s!"{name}: PARAM ERROR {e}"; return false
    | .ok r =>
      let t1 ← IO.monoMsNow
      let cones := (List.range ops.nv).filter (fun v => r.mu.get! v != 0.0) |>.length
      let good := match r.status with | .converged => true | .infeasible _ => true | _ => false
      let slow := r.itersOne.any (· ≥ 60)
      if verbose || !good || slow then
        IO.println s!"{name}: nv={ops.nv} n={r.nInterior} chi={ops.topo.chi} {describe p} status={repr r.status} stop={repr r.stop} it1={r.itersOne} it2={r.itersTwo} cones={cones} obj={r.cert.objective} resP={r.cert.resPrimal} resA={r.cert.resAdjoint} resI={r.cert.resInclusion} resB={r.cert.resBounds} gap={r.cert.gap} shift=[{r.cert.shiftLo},{r.cert.shiftHi}] ms={t1 - t0}"
      return good

def main (args : List String) : IO UInt32 := do
  let count := (args.head? >>= String.toNat?).getD 50
  let size := (args[1]? >>= String.toNat?).getD 30
  let verbose := args.length > 2
  let mut bad := 0
  for seed in [0:count] do
    let (kind, m) := randomMesh size seed
    let p := randomParams m seed
    let ok ← report s!"{kind}#{seed}" m p verbose
    if !ok then bad := bad + 1
  IO.println s!"stress: {count} cases, {bad} bad"
  return 0
