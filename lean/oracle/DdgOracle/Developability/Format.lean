import DdgOracle.Cuts.Format
import DdgOracle.Developability.Step

/-!
# Text format of the `dev` command

Same conventions as the other tasks (README, "Text format"): version 1, line oriented, floats as 16
hex digits, 0-based indices, one `case NAME ... end` block per case. Specified in the README section
"Developability".
-/

namespace Ddg.Developability

/-- One case. -/
structure DevCase where
  /-- identifier -/
  name : String
  /-- `energy`, `step` or `subdivide` -/
  task : String
  /-- energies to evaluate (`step` uses the first) -/
  energies : Array EnergyKind := #[.covariance]
  /-- count boundary vertices -/
  includeBoundary : Bool := true
  /-- leave out interior vertices of valence 3 -/
  omit3 : Bool := true
  /-- options of the step -/
  stepOpts : StepOpts := {}
  /-- the mesh -/
  mesh : DMesh
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

private def boolTok (ctx tok : String) : Except String Bool :=
  if tok == "0" then pure false else if tok == "1" then pure true
  else throw s!"{ctx}: expected 0 or 1, got '{tok}'"

/-- Parse the lines of one case (between `case` and `end`). -/
def parseBody (name : String) (lines : Array (Array String)) : Except String DevCase := do
  let mut task := ""
  let mut energies : Array EnergyKind := #[]
  let mut incl := true
  let mut omit3 := true
  let mut so : StepOpts := {}
  let mut mesh : Option DMesh := none
  let mut i := 0
  while i < lines.size do
    let l := lines[i]!
    match l[0]! with
    | "task" =>
      if l.size != 2 then throw "task: expected 'task <energy|step|subdivide>'"
      if !["energy", "step", "subdivide"].contains l[1]! then throw s!"unknown task '{l[1]!}'"
      task := l[1]!
    | "energy" =>
      if l.size != 2 then throw "energy: expected one name or 'all'"
      if l[1]! == "all" then energies := energies ++ EnergyKind.all
      else match EnergyKind.ofName? l[1]! with
        | some e => energies := energies.push e
        | none => throw s!"unknown energy '{l[1]!}'"
    | "boundary" =>
      if l.size != 2 then throw "boundary: expected 0 or 1"
      incl ← boolTok "boundary" l[1]!
    | "omit3" =>
      if l.size != 2 then throw "omit3: expected 0 or 1"
      omit3 ← boolTok "omit3" l[1]!
    | "step" =>
      if l.size != 7 then throw "step: expected 'step <initial> <max_displacement> <c1> <c2> <max_doublings> <max_bisections>'"
      so := { so with initialStep := ← floatTok "step" l[1]!, maxDisplacement := ← floatTok "step" l[2]!,
                      c1 := ← floatTok "step" l[3]!, c2 := ← floatTok "step" l[4]!,
                      maxDoublings := ← natTok "step" l[5]!, maxBisections := ← natTok "step" l[6]! }
    | "pinned" =>
      let mut ps := so.pinned
      for k in [1:l.size] do ps := ps.push (← natTok "pinned" l[k]!)
      so := { so with pinned := ps }
    | "mesh" =>
      if l.size != 3 then throw "mesh: expected 'mesh <nv> <nf>'"
      let nv ← natTok "mesh" l[1]!
      let nf ← natTok "mesh" l[2]!
      if nv > maxVertices then throw s!"mesh has {nv} vertices; the oracle accepts at most {maxVertices}"
      if i + nv + nf ≥ lines.size then throw "mesh: missing lines"
      let mut pos : Array V3 := #[]
      let mut faces : Array (Nat × Nat × Nat) := #[]
      for k in [0:nv] do
        let p := lines[i + 1 + k]!
        if p[0]! != "p" || p.size != 4 then throw s!"malformed line '{joinTokens p}': expected 'p <x> <y> <z>'"
        pos := pos.push ⟨← floatTok "p" p[1]!, ← floatTok "p" p[2]!, ← floatTok "p" p[3]!⟩
      for k in [0:nf] do
        let f := lines[i + 1 + nv + k]!
        if f[0]! != "f" || f.size != 4 then throw s!"malformed line '{joinTokens f}': expected 'f <a> <b> <c>'"
        faces := faces.push (← natTok "f" f[1]!, ← natTok "f" f[2]!, ← natTok "f" f[3]!)
      mesh := some { pos, faces }
      i := i + nv + nf
    | other => throw s!"unknown line '{other}'"
    i := i + 1
  if task == "" then throw "missing 'task' line"
  let some msh := mesh | throw "missing 'mesh' block"
  if energies.isEmpty then energies := #[.covariance]
  if task == "step" && energies.size != 1 then throw "task step needs exactly one 'energy' line"
  return { name, task, energies, includeBoundary := incl, omit3, stepOpts := so, mesh := msh }

/-- The case in the input format. -/
def printCase (c : DevCase) : String := Id.run do
  let mut s := s!"case {c.name}\ntask {c.task}\n"
  if c.task != "subdivide" then
    for e in c.energies do s := s ++ s!"energy {e.name}\n"
    s := s ++ s!"boundary {if c.includeBoundary then 1 else 0}\nomit3 {if c.omit3 then 1 else 0}\n"
  if c.task == "step" then
    let o := c.stepOpts
    s := s ++ s!"step {hx o.initialStep} {hx o.maxDisplacement} {hx o.c1} {hx o.c2} {o.maxDoublings} {o.maxBisections}\n"
    if o.pinned.size > 0 then s := s ++ "pinned " ++ " ".intercalate (o.pinned.map toString).toList ++ "\n"
  s := s ++ s!"mesh {c.mesh.pos.size} {c.mesh.faces.size}\n"
  for p in c.mesh.pos do s := s ++ s!"p {hx p.x} {hx p.y} {hx p.z}\n"
  for (a, b, cc) in c.mesh.faces do s := s ++ s!"f {a} {b} {cc}\n"
  return s ++ "end\n"

/-- Output block of a case that could not be processed. -/
def printError (name task msg : String) : String := s!"case {name}\ntask {task}\nstatus error {msg}\nend\n"

private def b2n (b : Bool) : String := if b then "1" else "0"

private def printEnergy (c : DevCase) (topo : Topo) : String := Id.run do
  let m := c.mesh
  let nb := (topo.stars.filter (!·.closed)).size
  let mut s := s!"case {c.name}\ntask energy\nstatus ok\nmesh {topo.nv} {m.faces.size} {topo.edges.size} {nb}\n"
  for e in c.energies do
    s := s ++ s!"section {e.name}\n"
    match evaluate m topo { energy := e, includeBoundary := c.includeBoundary, omitValenceThree := c.omit3 } with
    | .error err => s := s ++ s!"result error {err.text}\n"
    | .ok ev =>
      s := s ++ s!"result ok {hx ev.total}\n"
      for v in [0:topo.nv] do
        s := s ++ s!"v {v} {b2n ev.counted[v]!} {ev.valence[v]!} {b2n ev.boundary[v]!} {hx ev.e[v]!} {hx ev.margin[v]!}\n"
        let g := ev.grad[v]!
        s := s ++ s!"g {v} {hx g.x} {hx g.y} {hx g.z}\n"
        s := s ++ s!"ch {v}" ++ (ev.choice[v]!.foldl (fun acc k => acc ++ s!" {k}") "") ++ "\n"
  return s ++ "end\n"

private def printStep (c : DevCase) (topo : Topo) : String := Id.run do
  let o : Opts := { energy := c.energies[0]!, includeBoundary := c.includeBoundary, omitValenceThree := c.omit3 }
  match step c.mesh topo o c.stepOpts with
  | .error err => return printError c.name "step" err.text
  | .ok r =>
    let mut s := s!"case {c.name}\ntask step\nstatus ok\n"
    s := s ++ s!"lbar {hx r.lbar}\ngmax {hx r.gmax}\nenergy0 {hx r.e0}\nslope {hx r.slope}\ntmax {hx r.tmax}\nmargin {hx r.margin}\n"
    for i in [0:r.trials.size] do
      let t := r.trials[i]!
      let cc := match t.c with | none => "-1" | some b => b2n b
      s := s ++ s!"trial {i} {hx t.t} {hx t.e} {b2n t.s} {cc}\n"
    s := s ++ s!"result {r.kind.name} {hx r.t} {r.trials.size} {hx r.e1}\n"
    for v in [0:topo.nv] do
      let p := r.pos[v]!
      s := s ++ s!"x {v} {hx p.x} {hx p.y} {hx p.z}\n"
    return s ++ "end\n"

private def printSubdivide (c : DevCase) (topo : Topo) : String := Id.run do
  let (m2, mids) := subdivide c.mesh topo
  let mut s := s!"case {c.name}\ntask subdivide\nstatus ok\nmesh {m2.pos.size} {m2.faces.size} {mids.size}\n"
  for (a, b, id) in mids do s := s ++ s!"em {a} {b} {id}\n"
  for i in [0:m2.pos.size] do
    let p := m2.pos[i]!
    s := s ++ s!"x {i} {hx p.x} {hx p.y} {hx p.z}\n"
  for (a, b, cc) in m2.faces do s := s ++ s!"f {a} {b} {cc}\n"
  return s ++ "end\n"

/-- Run one case and print its output block. -/
def runCase (c : DevCase) : String :=
  match buildTopo c.mesh with
  | .error e => printError c.name c.task e
  | .ok topo =>
    match c.task with
    | "energy" => printEnergy c topo
    | "step" => printStep c topo
    | _ => printSubdivide c topo

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

end Ddg.Developability
