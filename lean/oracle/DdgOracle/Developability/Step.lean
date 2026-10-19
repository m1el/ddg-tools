import DdgOracle.Developability.Energy

/-!
# One L-BFGS step, and 4-1 subdivision

`step` is the first iteration of L-BFGS in SPEC 4.2 (the memory is empty, so the direction is scaled
steepest descent, SPEC 4.1) with the Lewis-Overton weak Wolfe line search of SPEC 4.3.
`subdivide` is SPEC 5.2.
-/

namespace Ddg.Developability

/-- Options of one step (SPEC 4.1, 4.3). -/
structure StepOpts where
  /-- `initial_step` -/
  initialStep : Float := 0.1
  /-- `max_displacement` -/
  maxDisplacement : Float := 0.5
  /-- `c1` -/
  c1 : Float := 0.0
  /-- `c2` -/
  c2 : Float := 0.5
  /-- `max_doublings` -/
  maxDoublings : Nat := 30
  /-- `max_bisections` -/
  maxBisections : Nat := 30
  /-- pinned vertices -/
  pinned : Array Nat := #[]
  deriving Inhabited

/-- One trial of the line search. -/
structure Trial where
  /-- step -/
  t : Float
  /-- energy at the trial point (NaN if undefined) -/
  e : Float
  /-- Armijo condition `S` -/
  s : Bool
  /-- curvature condition `C`, evaluated only if `S` holds -/
  c : Option Bool
  /-- margin of the Armijo test, `(c1 s t - h) / |s t|`; positive iff `S` holds -/
  sm : Float := 0
  /-- margin of the curvature test, `(g(t)·d - c2 s) / |s|`; `0` if not evaluated -/
  cm : Float := 0
  deriving Inhabited

/-- How the line search ended. -/
inductive StepKind
  | wolfe | capped | weak | failed | stationary
  deriving BEq

/-- Name in the output. -/
def StepKind.name : StepKind → String
  | .wolfe => "wolfe" | .capped => "capped" | .weak => "weak" | .failed => "failed"
  | .stationary => "stationary"

/-- Outcome of a step. -/
structure StepResult where
  /-- mean edge length -/
  lbar : Float
  /-- largest per-vertex gradient norm -/
  gmax : Float
  /-- energy at the start -/
  e0 : Float
  /-- `g·d` -/
  slope : Float
  /-- cap on the step -/
  tmax : Float
  /-- trials -/
  trials : Array Trial
  /-- how it ended -/
  kind : StepKind
  /-- accepted step (0 if none) -/
  t : Float
  /-- energy afterwards -/
  e1 : Float
  /-- new positions -/
  pos : Array V3
  /-- smallest margin over the counted vertices, relative to `max 1 |e_v|` -/
  margin : Float

/-- Mean length of the unique edges. -/
def meanEdgeLength (m : DMesh) (topo : Topo) : Float :=
  (topo.edges.foldl (fun acc (a, b) => acc + (m.pos[a]! - m.pos[b]!).norm) 0) / topo.edges.size.toFloat

private def inf : Float := posInf

/-- One step from the positions of `m`. -/
def step (m : DMesh) (topo : Topo) (o : Opts) (so : StepOpts) : Except EvalError StepResult := do
  let ev0 ← evaluate m topo o
  let margin := (Array.range topo.nv).foldl (fun acc v => if ev0.counted[v]! then min acc (ev0.margin[v]! / max 1 ev0.e[v]!.abs) else acc) posInf
  let pinned := fun v => so.pinned.contains v
  let g := Array.ofFn (n := topo.nv) fun v => if pinned v.val then V3.zero else ev0.grad[v]!
  let lbar := meanEdgeLength m topo
  let gmax := g.foldl (fun acc x => max acc x.norm) 0
  if !(gmax > 0) then
    return { lbar, gmax, e0 := ev0.total, slope := 0, tmax := posInf, trials := #[], kind := .stationary,
             t := 0, e1 := ev0.total, pos := m.pos, margin }
  let scale := so.initialStep * lbar / gmax
  let d := g.map (V3.smul (-scale))
  let dinf := d.foldl (fun acc x => max acc x.norm) 0
  let tmax := so.maxDisplacement * lbar / dinf
  let slope := (Array.range topo.nv).foldl (fun acc v => acc + V3.dot g[v]! d[v]!) 0
  let moved := fun (t : Float) => ({ m with pos := Array.ofFn (n := topo.nv) fun v => m.pos[v]! + V3.smul t d[v]! } : DMesh)
  let mut alpha := 0.0
  let mut beta := posInf
  let mut t := min 1.0 tmax
  let mut nd := 0
  let mut nb := 0
  let mut trials : Array Trial := #[]
  let mut kind := StepKind.failed
  let mut accepted := 0.0
  let mut eAcc := nan
  let mut eAlpha := nan
  let mut done := false
  let mut fuel := so.maxDoublings + so.maxBisections + 4
  while !done && fuel > 0 do
    fuel := fuel - 1
    let mt := moved t
    let (e, sOk, cOk, sm, cm) :=
      match evaluate mt topo o with
      | .error _ => (nan, false, none, -inf, 0.0)
      | .ok ev =>
        let h := ev.total - ev0.total
        let smv := (so.c1 * slope * t - h) / (slope.abs * t)
        if h < so.c1 * slope * t then
          let gt := Array.ofFn (n := topo.nv) fun v => if pinned v.val then V3.zero else ev.grad[v]!
          let gd := (Array.range topo.nv).foldl (fun acc v => acc + V3.dot gt[v]! d[v]!) 0
          (ev.total, true, some (decide (gd > so.c2 * slope)), smv, (gd - so.c2 * slope) / slope.abs)
        else (ev.total, false, none, smv, 0.0)
    trials := trials.push { t, e, s := sOk, c := cOk, sm, cm }
    if !sOk then beta := t
    else if cOk == some false then
      alpha := t
      eAlpha := e
    else
      kind := .wolfe
      accepted := t
      eAcc := e
      done := true
    if !done then
      if beta < posInf then
        if nb ≥ so.maxBisections then
          done := true
          if alpha > 0 then
            kind := .weak
            accepted := alpha
            eAcc := eAlpha
        else
          nb := nb + 1
          t := (alpha + beta) / 2
      else if alpha ≥ tmax then
        kind := .capped
        accepted := alpha
        eAcc := eAlpha
        done := true
      else if nd ≥ so.maxDoublings then
        done := true
        if alpha > 0 then
          kind := .weak
          accepted := alpha
          eAcc := eAlpha
      else
        nd := nd + 1
        t := min (2 * alpha) tmax
  let pos := if kind == .failed then m.pos else (moved accepted).pos
  return { lbar, gmax, e0 := ev0.total, slope, tmax, trials, kind, t := accepted,
           e1 := if kind == .failed then ev0.total else eAcc, pos, margin }

/-- 4-1 subdivision: returns the new mesh and the edge midpoints `(a, b, id)`. -/
def subdivide (m : DMesh) (topo : Topo) : DMesh × Array (Nat × Nat × Nat) := Id.run do
  let nv := topo.nv
  let mut mids : Array (Nat × Nat × Nat) := #[]
  let mut table : Array Nat := Array.replicate (nv * nv) 0
  let mut pos := m.pos
  for i in [0:topo.edges.size] do
    let (a, b) := topo.edges[i]!
    mids := mids.push (a, b, nv + i)
    table := (table.set! (a * nv + b) (nv + i + 1)).set! (b * nv + a) (nv + i + 1)
    pos := pos.push (V3.smul 0.5 (m.pos[a]! + m.pos[b]!))
  let mut faces : Array (Nat × Nat × Nat) := #[]
  for f in m.faces do
    let (a, b, c) := f
    let mab := table[a * nv + b]! - 1
    let mbc := table[b * nv + c]! - 1
    let mca := table[c * nv + a]! - 1
    faces := faces.push (a, mab, mca)
    faces := faces.push (b, mbc, mab)
    faces := faces.push (c, mca, mbc)
    faces := faces.push (mab, mbc, mca)
  return ({ pos, faces }, mids)

end Ddg.Developability
