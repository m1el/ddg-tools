import DdgOracle
open Ddg

def stateStr (st : Array ConeState) : String :=
  String.ofList (st.toList.map fun s => match s with
    | .zero => '.' | .pos => '+' | .neg => '-' | .lower => 'L' | .upper => 'U')

def main (args : List String) : IO UInt32 := do
  let seed := (args.head? >>= String.toNat?).getD 19
  let size := (args[1]? >>= String.toNat?).getD 30
  let mode := (args[2]?).getD "trace"
  let (_, m) := randomMesh size seed
  let p := randomParams m seed
  match m.operators with
  | .error e => IO.println e; return 1
  | .ok ops =>
    let s := buildSetup ops p
    let pr := s.problem
    let n := s.n
    let opts : Options := {}
    if mode == "sol" then
      match solve ops p with
      | .error e => IO.println e
      | .ok r =>
        IO.println s!"status {repr r.status} stop {repr r.stop} it1 {r.itersOne} it2 {r.itersTwo}"
        for i in [0:n] do
          let v := s.interior[i]!
          let f := r.phi.get! v
          let mm := r.mu.get! v
          IO.println s!"{v}: t={s.t.get! i} lo={s.lo.get! i} hi={s.hi.get! i} wR={p.wR.get! v} lo0={p.lo.get! v} hi0={p.hi.get! v} u={r.u.get! v} phi={f} mu={mm} F3={pr.coneResidual (vecFn (FloatArray.ofFn' n fun j => r.phi.get! s.interior[j]!)) (vecFn (FloatArray.ofFn' n fun j => r.mu.get! s.interior[j]!)) i} viol={subdiffViolation (p.lambda * p.wR.get! v) (p.lo.get! v) (p.hi.get! v) mm f}"
      return 0
    let mut u := zeros n
    let mut phi := zeros n
    let mut mu := zeros n
    let mut sigma : Option (Array ConeState) := none
    for gamma in opts.gammas.push 0.0 do
      IO.println s!"gamma={gamma}"
      for it in [0:14] do
        let states := arrayOf n (pr.impState gamma (vecFn phi) (vecFn mu))
        let residual := FloatArray.ofFn' (3 * n) (pr.Fimp gamma (vecFn u) (vecFn phi) (vecFn mu))
        IO.println s!"  it={it} target={stateStr states} used={(sigma.map stateStr).getD "-"} |F|={maxAbs residual}"
        let r := newtonStage s gamma 1 0 u phi mu sigma
        u := r.u; phi := r.phi; mu := r.mu; sigma := r.sigma
        if r.iters == 0 then
          IO.println s!"  stop {repr r.stop}"
          break
    return 0
