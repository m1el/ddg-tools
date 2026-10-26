import DdgOracle.Cuts.Format
import DdgOracle.RepulsiveCurves.Step

/-!
# Text format of the `rc` command

Same conventions as the other tasks (README, "Text format"): version 1, line oriented, floats as 16
hex digits, 0-based indices, one `case NAME ... end` block per case. Specified in the README section
"Repulsive curves".
-/

namespace Ddg.RepulsiveCurves

/-- One case. -/
structure RCCase where
  /-- identifier -/
  name : String
  /-- `energy`, `gram`, `saddle` or `project` -/
  task : String
  /-- exponents -/
  ex : Expo := ⟨3.0, 6.0⟩
  /-- quadrature of `saddle` and `project` (`gram` prints both) -/
  quad : Quad := .midpoint
  /-- held axis of a planar network -/
  planar : Option Nat := none
  /-- `gram`: also print `B`, `B⁰`, `K^B`, `K⁰` (otherwise only `A`) -/
  allMatrices : Bool := true
  /-- the network -/
  net : Net
  /-- constraints (`saddle`, `project`) -/
  cons : Array Cons := #[]
  /-- step size `τ` of `project`: trial `γ − τ g̃` -/
  tau : Float := 0.0
  /-- `step` is a fraction: `τ = tau · ℓ̄ / max|g̃|` (so the trial moves the curve by a fraction of
  an edge length whatever the scale of the gradient) -/
  tauRel : Bool := false
  /-- line search of `step` -/
  search : SearchOpts := {}
  /-- tolerance relative to the mean edge length (else absolute) -/
  tolMean : Bool := true
  /-- tolerance value -/
  tolVal : Float := 1e-4
  /-- maximal chord iterations -/
  maxit : Nat := 3
  deriving Inhabited

private def hx (x : Float) : String := floatToHex x
private def joinTokens (l : Array String) : String := " ".intercalate l.toList

private def natTok (ctx tok : String) : Except String Nat :=
  match tok.toNat? with
  | some n => pure n
  | none => throw s!"{ctx}: expected a decimal integer, got '{tok}'"

private def floatTok (ctx tok : String) : Except String Float :=
  match hexToFloat? tok with
  | some x => pure x
  | none => throw s!"{ctx}: expected 16 hex digits (IEEE-754 bit pattern), got '{tok}'"

private def v3Toks (ctx : String) (l : Array String) (i : Nat) : Except String V3 := do
  return ⟨← floatTok ctx l[i]!, ← floatTok ctx l[i+1]!, ← floatTok ctx l[i+2]!⟩

private def natList (ctx : String) (l : Array String) (i : Nat) : Except String (Array Nat) := do
  let mut r : Array Nat := #[]
  for k in [i:l.size] do r := r.push (← natTok ctx l[k]!)
  return r

/-- Largest network the oracle accepts. -/
def maxVertices : Nat := 400

private def parseCons (l : Array String) : Except String Cons := do
  if l.size < 2 then throw "c: expected a constraint kind"
  match l[1]! with
  | "barycenter" =>
    if l.size == 2 then return .barycenter none
    if l.size == 5 then return .barycenter (some (← v3Toks "c barycenter" l 2))
    throw "c barycenter: expected no target or three coordinates"
  | "total-length" =>
    if l.size == 2 then return .totalLength none
    if l.size == 3 then return .totalLength (some (← floatTok "c total-length" l[2]!))
    throw "c total-length: expected at most one target"
  | "edge-lengths" =>
    if l.size == 3 && l[2]! == "all" then return .edgeLengths none
    return .edgeLengths (some (← natList "c edge-lengths" l 2))
  | "point" =>
    if l.size == 3 then return .point (← natTok "c point" l[2]!) none
    if l.size == 6 then return .point (← natTok "c point" l[2]!) (some (← v3Toks "c point" l 3))
    throw "c point: expected 'c point V [X Y Z]'"
  | "sphere" =>
    if l.size < 7 then throw "c sphere: expected 'c sphere CX CY CZ R all|V V ...'"
    let c ← v3Toks "c sphere" l 2
    let r ← floatTok "c sphere" l[5]!
    if l.size == 7 && l[6]! == "all" then return .sphere c r none
    return .sphere c r (some (← natList "c sphere" l 6))
  | "tangent" =>
    if l.size != 4 && l.size != 7 then throw "c tangent: expected 'c tangent edge|vertex N [X Y Z]'"
    let n ← natTok "c tangent" l[3]!
    let site ← match l[2]! with
      | "edge" => pure (Site.edge n)
      | "vertex" => pure (Site.vertex n)
      | other => throw s!"c tangent: unknown site '{other}'"
    let dir ← if l.size == 7 then (do let d ← v3Toks "c tangent" l 4; pure (some d)) else pure none
    return .tangent site dir
  | other => throw s!"c: unknown constraint '{other}'"

/-- Parse the lines of one case (between `case` and `end`). -/
def parseBody (name : String) (lines : Array (Array String)) : Except String RCCase := do
  let mut task := ""
  let mut ex : Expo := ⟨3.0, 6.0⟩
  let mut quad : Quad := .midpoint
  let mut planar : Option Nat := none
  let mut allM := true
  let mut cons : Array Cons := #[]
  let mut tau : Float := 0.0
  let mut tauRel := false
  let mut search : SearchOpts := {}
  let mut tolMean := true
  let mut tolVal : Float := 1e-4
  let mut maxit := 3
  let mut net : Option Net := none
  let mut i := 0
  while i < lines.size do
    let l := lines[i]!
    match l[0]! with
    | "task" =>
      if l.size != 2 || !["energy", "gram", "saddle", "project", "step"].contains l[1]! then
        throw "task: expected 'task energy|gram|saddle|project|step'"
      task := l[1]!
    | "exponents" =>
      if l.size != 3 then throw "exponents: expected 'exponents ALPHA BETA'"
      ex := ⟨← floatTok "exponents" l[1]!, ← floatTok "exponents" l[2]!⟩
    | "quadrature" =>
      if l.size != 2 then throw "quadrature: expected midpoint|trapezoid"
      quad ← match l[1]! with
        | "midpoint" => pure Quad.midpoint
        | "trapezoid" => pure Quad.trapezoid
        | other => throw s!"unknown quadrature '{other}'"
    | "planar" =>
      if l.size != 2 then throw "planar: expected an axis 0, 1 or 2"
      let a ← natTok "planar" l[1]!
      if a > 2 then throw "planar: axis must be 0, 1 or 2"
      planar := some a
    | "matrices" =>
      if l.size != 2 then throw "matrices: expected all|a"
      allM := l[1]! == "all"
    | "step" =>
      if l.size != 2 then throw "step: expected 'step TAU'"
      tau ← floatTok "step" l[1]!
    | "steprel" =>
      if l.size != 2 then throw "steprel: expected 'steprel FRACTION'"
      tau ← floatTok "steprel" l[1]!
      tauRel := true
    | "search" =>
      if l.size != 5 then throw "search: expected 'search ARMIJO SHRINK MAXHALVINGS NSTEPS'"
      search := { armijo := ← floatTok "search" l[1]!, shrink := ← floatTok "search" l[2]!,
                  maxHalvings := ← natTok "search" l[3]!, nSteps := ← natTok "search" l[4]! }
    | "tol" =>
      if l.size != 3 then throw "tol: expected 'tol mean|abs VALUE'"
      tolMean := l[1]! == "mean"
      tolVal ← floatTok "tol" l[2]!
    | "maxit" =>
      if l.size != 2 then throw "maxit: expected a count"
      maxit ← natTok "maxit" l[1]!
    | "c" => cons := cons.push (← parseCons l)
    | "net" =>
      if l.size != 3 then throw "net: expected 'net NV NE'"
      let nv ← natTok "net" l[1]!
      let ne ← natTok "net" l[2]!
      if nv > maxVertices then throw s!"net has {nv} vertices; the oracle accepts at most {maxVertices}"
      if i + nv + ne ≥ lines.size then throw "net: missing lines"
      let mut pos : Array V3 := #[]
      let mut edges : Array (Nat × Nat) := #[]
      for k in [0:nv] do
        let p := lines[i + 1 + k]!
        if p[0]! != "p" || p.size != 4 then throw s!"malformed line '{joinTokens p}': expected 'p <x> <y> <z>'"
        pos := pos.push (← v3Toks "p" p 1)
      for k in [0:ne] do
        let e := lines[i + 1 + nv + k]!
        if e[0]! != "e" || e.size != 3 then throw s!"malformed line '{joinTokens e}': expected 'e <a> <b>'"
        edges := edges.push (← natTok "e" e[1]!, ← natTok "e" e[2]!)
      net := some { pos, edges }
      i := i + nv + ne
    | other => throw s!"unknown line '{other}'"
    i := i + 1
  if task == "" then throw "missing 'task' line"
  let some nt := net | throw "missing 'net' block"
  return { name, task, ex, quad, planar, allMatrices := allM, net := nt, cons, tau, tauRel, search, tolMean, tolVal, maxit }

private def consLine : Cons → String
  | .barycenter none => "c barycenter\n"
  | .barycenter (some t) => s!"c barycenter {hx t.x} {hx t.y} {hx t.z}\n"
  | .totalLength none => "c total-length\n"
  | .totalLength (some t) => s!"c total-length {hx t}\n"
  | .edgeLengths none => "c edge-lengths all\n"
  | .edgeLengths (some es) => "c edge-lengths " ++ " ".intercalate (es.map toString).toList ++ "\n"
  | .point v none => s!"c point {v}\n"
  | .point v (some t) => s!"c point {v} {hx t.x} {hx t.y} {hx t.z}\n"
  | .sphere c r none => s!"c sphere {hx c.x} {hx c.y} {hx c.z} {hx r} all\n"
  | .sphere c r (some vs) => s!"c sphere {hx c.x} {hx c.y} {hx c.z} {hx r} " ++ " ".intercalate (vs.map toString).toList ++ "\n"
  | .tangent site dir =>
    let (kind, n) := match site with | .edge e => ("edge", e) | .vertex v => ("vertex", v)
    s!"c tangent {kind} {n}" ++ (match dir with | none => "" | some d => s!" {hx d.x} {hx d.y} {hx d.z}") ++ "\n"

/-- The case in the input format. -/
def printCase (c : RCCase) : String := Id.run do
  let mut s := s!"case {c.name}\ntask {c.task}\nexponents {hx c.ex.alpha} {hx c.ex.beta}\n"
  if c.task == "saddle" || c.task == "project" || c.task == "step" then s := s ++ s!"quadrature {c.quad.name}\n"
  if c.task == "gram" && !c.allMatrices then s := s ++ "matrices a\n"
  if let some a := c.planar then s := s ++ s!"planar {a}\n"
  if c.task == "project" then
    s := s ++ s!"{if c.tauRel then "steprel" else "step"} {hx c.tau}\ntol {if c.tolMean then "mean" else "abs"} {hx c.tolVal}\nmaxit {c.maxit}\n"
  if c.task == "step" then
    s := s ++ s!"search {hx c.search.armijo} {hx c.search.shrink} {c.search.maxHalvings} {c.search.nSteps}\ntol {if c.tolMean then "mean" else "abs"} {hx c.tolVal}\nmaxit {c.maxit}\n"
  for k in c.cons do s := s ++ consLine k
  s := s ++ s!"net {c.net.nv} {c.net.ne}\n"
  for p in c.net.pos do s := s ++ s!"p {hx p.x} {hx p.y} {hx p.z}\n"
  for (a, b) in c.net.edges do s := s ++ s!"e {a} {b}\n"
  return s ++ "end\n"

/-- Output block of a case that could not be processed; `msg` starts with the error kind. -/
def printError (name task msg : String) : String :=
  let kind := (msg.splitOn ":").head!
  let rest := (msg.drop (kind.length + 1)).trimAscii.toString
  s!"case {name}\ntask {task}\nstatus error {kind} {rest}\nend\n"

private def row (tag : String) (i : Nat) (a : FloatArray) (off len : Nat) : String := Id.run do
  let mut s := s!"{tag} {i}"
  for j in [0:len] do s := s ++ " " ++ hx (a.get! (off + j))
  return s ++ "\n"

private def matrix (tag : String) (a : FloatArray) (rows cols : Nat) : String := Id.run do
  let mut s := ""
  for i in [0:rows] do s := s ++ row tag i a (i * cols) cols
  return s

/-- Validate and evaluate the energy of a case (shared by all tasks). -/
def common (c : RCCase) : Except String (Array V3 × Float × Array V3) := do
  c.ex.check
  c.net.validate
  if let some ax := c.planar then
    let h := (c.net.pos[0]!).get ax
    for v in [0:c.net.nv] do
      if (c.net.pos[v]!).get ax != h then throw s!"NotPlanar: vertex {v} is off the plane"
  let (e, d) := energyAndDifferential c.net c.ex
  let coords := coordsOf c.planar
  -- no component along the held axis
  let d := d.map fun p => Id.run do
    let mut q := p
    for ax in [0:3] do
      if !coords.contains ax then q := q.set ax 0.0
    return q
  return (c.net.pos, e, d)

private def printEnergy (c : RCCase) (e : Float) (d : Array V3) : String := Id.run do
  let mut s := s!"case {c.name}\ntask energy\nstatus ok\nnet {c.net.nv} {c.net.ne} {c.net.nPairs}\nenergy {hx e}\n"
  for v in [0:d.size] do s := s ++ s!"d {v} {hx d[v]!.x} {hx d[v]!.y} {hx d[v]!.z}\n"
  return s ++ "end\n"

private def printGram (c : RCCase) : String := Id.run do
  let n := c.net.nv
  let m := c.net.ne
  let mut s := s!"case {c.name}\ntask gram\nstatus ok\nnet {n} {m} {c.net.nPairs}\n"
  for q in [Quad.midpoint, Quad.trapezoid] do
    let G := gram c.net c.ex q
    s := s ++ s!"section {q.name}\n"
    if c.allMatrices then
      s := s ++ matrix "B" G.b n n ++ matrix "B0" G.b0 n n
    s := s ++ matrix "A" G.a n n
    if c.allMatrices then
      s := s ++ matrix "KB" G.kb m m ++ matrix "K0" G.k0 m m
  return s ++ "end\n"

private def printSaddle (c : RCCase) (d : Array V3) : Except String String := do
  let coords := coordsOf c.planar
  let rc ← resolve c.net c.net.pos c.planar c.cons
  let rows := evalRows c.net c.net.pos coords rc
  let G := gram c.net c.ex c.quad
  let sad ← Saddle.build c.net.nv coords G.a rows
  let gr := projectedGradient sad d
  let mut s := s!"case {c.name}\ntask saddle\nstatus ok\nnet {c.net.nv} {c.net.ne} {c.net.nPairs}\nquadrature {c.quad.name}\n"
  s := s ++ s!"rows {sad.k} {sad.nd}\npivot {hx sad.pivotRatio}\ncond {hx sad.cond}\n"
  s := s ++ "phi" ++ (rows.foldl (fun acc r => acc ++ " " ++ hx r.phi) "") ++ "\n"
  s := s ++ matrix "c" sad.c sad.k sad.nd
  for v in [0:gr.g.size] do s := s ++ s!"g {v} {hx gr.g[v]!.x} {hx gr.g[v]!.y} {hx gr.g[v]!.z}\n"
  s := s ++ "lambda" ++ ((Array.range sad.k).foldl (fun acc i => acc ++ " " ++ hx (gr.lam.get! i)) "") ++ "\n"
  s := s ++ s!"dfg {hx gr.dfg}\ngnorm {hx gr.gnorm}\n"
  return s ++ "end\n"

private def printStep (c : RCCase) : Except String String := do
  let coords := coordsOf c.planar
  let rc ← resolve c.net c.net.pos c.planar c.cons
  let lbar := (c.net.geo.len.foldl (· + ·) 0.0) / c.net.ne.toFloat
  let tol := if c.tolMean then c.tolVal * lbar else c.tolVal
  let mut pos := c.net.pos
  let mut prev : Option Float := none
  let mut recs : Array StepRec := #[]
  for _ in [0:c.search.nSteps] do
    let r ← normalizedStep c.net c.ex c.quad coords rc c.search tol lbar c.maxit pos prev
    recs := recs.push r
    if !r.accepted then break
    pos := r.pos
    prev := some r.tau
  let mut s := s!"case {c.name}\ntask step\nstatus ok\nnet {c.net.nv} {c.net.ne} {c.net.nPairs}\nquadrature {c.quad.name}\n"
  s := s ++ s!"rows {rc.size} {coords.size * c.net.nv}\nlbar {hx lbar}\nsteps {recs.size}\n"
  for k in [0:recs.size] do
    let r := recs[k]!
    s := s ++ s!"step {k} {hx r.tau0} {hx r.slope} {hx r.f0} {hx r.gnorm} {hx r.pivot} {hx r.cond}\n"
    for i in [0:r.trials.size] do
      let t := r.trials[i]!
      s := s ++ s!"trial {k} {i} {hx t.tau} {t.outcome} {t.iters} {hx t.energy} {hx t.viol} {hx t.margin}\n"
    s := s ++ s!"result {k} {if r.accepted then "accepted" else "stuck"} {hx r.tau} {r.halvings} {hx r.e1}\n"
    for v in [0:r.pos.size] do s := s ++ s!"x {k} {v} {hx r.pos[v]!.x} {hx r.pos[v]!.y} {hx r.pos[v]!.z}\n"
  return s ++ "end\n"

private def printProject (c : RCCase) (d : Array V3) : Except String String := do
  let coords := coordsOf c.planar
  let rc ← resolve c.net c.net.pos c.planar c.cons
  let rows := evalRows c.net c.net.pos coords rc
  let G := gram c.net c.ex c.quad
  let sad ← Saddle.build c.net.nv coords G.a rows
  let gr := projectedGradient sad d
  let lbar := (c.net.geo.len.foldl (· + ·) 0.0) / c.net.ne.toFloat
  let tol := if c.tolMean then c.tolVal * lbar else c.tolVal
  let gmax := gr.g.foldl (fun acc p => max acc (max p.x.abs (max p.y.abs p.z.abs))) 0.0
  let step := if c.tauRel then (if gmax > 0.0 then c.tau * lbar / gmax else 0.0) else c.tau
  let trial := (Array.range c.net.nv).map fun v => c.net.pos[v]! - V3.smul step gr.g[v]!
  let r := chordProject c.net rc sad trial tol lbar c.maxit
  let mut s := s!"case {c.name}\ntask project\nstatus ok\nnet {c.net.nv} {c.net.ne} {c.net.nPairs}\nquadrature {c.quad.name}\n"
  s := s ++ s!"rows {sad.k} {sad.nd}\npivot {hx sad.pivotRatio}\ncond {hx sad.cond}\nlbar {hx lbar}\ntau {hx step}\n"
  s := s ++ "phi" ++ (rows.foldl (fun acc r => acc ++ " " ++ hx r.phi) "") ++ "\n"
  for v in [0:gr.g.size] do s := s ++ s!"g {v} {hx gr.g[v]!.x} {hx gr.g[v]!.y} {hx gr.g[v]!.z}\n"
  s := s ++ s!"result {r.status.name} {r.iterations}\nmaxviol {hx r.maxViolation}\n"
  s := s ++ "hist" ++ (r.history.foldl (fun acc h => acc ++ " " ++ hx h) "") ++ "\n"
  for v in [0:r.pos.size] do s := s ++ s!"x {v} {hx r.pos[v]!.x} {hx r.pos[v]!.y} {hx r.pos[v]!.z}\n"
  return s ++ "end\n"

/-- Run one case and print its output block. -/
def runCase (c : RCCase) : String :=
  match common c with
  | .error e => printError c.name c.task e
  | .ok (_, e, d) =>
    match c.task with
    | "energy" => printEnergy c e d
    | "gram" => printGram c
    | "saddle" => match printSaddle c d with | .error err => printError c.name c.task err | .ok s => s
    | "step" => match printStep c with | .error err => printError c.name c.task err | .ok s => s
    | _ => match printProject c d with | .error err => printError c.name c.task err | .ok s => s

/-- Read cases until the end of `input`, write one block per case, flush after each. -/
partial def runStream (input output : IO.FS.Stream) : IO Nat := do
  let reader ← LineReader.new input
  let rec loop (count : Nat) : IO Nat := do
    match ← Cuts.readCaseLines reader with
    | none => return count
    | some (name, parsed) =>
      let block := match parsed with
        | .error e => printError name "?" e
        | .ok lines =>
          match parseBody name lines with
          | .error e => printError name "?" e
          | .ok c => runCase c
      output.putStr block
      output.flush
      loop (count + 1)
  loop 0

end Ddg.RepulsiveCurves
