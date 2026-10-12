import DdgOracle.Certificate
import DdgOracle.Cases

/-!
# The text format

Line oriented; every line is a tag followed by whitespace-separated tokens. Integers are decimal,
floats are the 16 hex digits of their IEEE-754 bit pattern (`DdgOracle.Hex`). The grammar is
specified in `README.md`; this file is its implementation.
-/

namespace Ddg

/-- What to compute for a case. -/
inductive Task where
  /-- curvatures, masses, cotan Laplacian -/
  | operators
  /-- solve the canonical problem -/
  | solve
  deriving BEq, Repr, Inhabited

/-- Name of a task in the text format. -/
def Task.name : Task → String
  | .operators => "operators"
  | .solve => "solve"

/-- A parsed input case. -/
structure CaseInput where
  /-- identifier -/
  name : String
  /-- task requested by a `task` line, if any -/
  task : Option Task
  /-- vertex positions -/
  pos : Array (Point Float)
  /-- triangles -/
  faces : Array Face
  /-- parameters, if the case has a parameter section -/
  params : Option Params
  deriving Inhabited

/-! ## Parsing -/

private def joinTokens (l : Array String) : String := " ".intercalate l.toList

private def natTok (ctx tok : String) : Except String Nat :=
  match tok.toNat? with
  | some n => pure n
  | none => throw s!"{ctx}: expected a decimal integer, got '{tok}'"

private def floatTok (ctx tok : String) : Except String Float :=
  match hexToFloat? tok with
  | some x => pure x
  | none => throw s!"{ctx}: expected 16 hex digits (IEEE-754 bit pattern), got '{tok}'"

/-- Parse the lines of one case (everything between the `case` line and the `end` line). -/
def parseBody (name : String) (lines : Array (Array String)) : Except String CaseInput := do
  let mut i := 0
  let line := fun (k : Nat) (what : String) =>
    if h : k < lines.size then pure (f := Except String) lines[k]
    else throw s!"missing line: expected {what}"
  -- optional task
  let mut task : Option Task := none
  let l ← line i "'mesh <nv> <nf>'"
  if l[0]! == "task" then
    if l.size != 2 then throw s!"malformed line '{joinTokens l}': expected 'task <operators|solve>'"
    task ← match l[1]! with
      | "operators" => pure (some Task.operators)
      | "solve" => pure (some Task.solve)
      | other => throw s!"unknown task '{other}'"
    i := i + 1
  let l ← line i "'mesh <nv> <nf>'"
  if l[0]! != "mesh" || l.size != 3 then
    throw s!"malformed line '{joinTokens l}': expected 'mesh <nv> <nf>'"
  let nv ← natTok "mesh" l[1]!
  let nf ← natTok "mesh" l[2]!
  if nv > maxVertices then
    throw s!"mesh has {nv} vertices; the dense oracle accepts at most {maxVertices}"
  if nf > 4 * maxVertices then throw s!"mesh has {nf} faces; too many for {maxVertices} vertices"
  i := i + 1
  let mut pos : Array (Point Float) := Array.mkEmpty nv
  for v in [0:nv] do
    let l ← line i s!"'v <x> <y> <z>' for vertex {v}"
    if l[0]! != "v" || l.size != 4 then
      throw s!"malformed line '{joinTokens l}': expected 'v <x> <y> <z>' for vertex {v}"
    let ctx := s!"vertex {v}"
    pos := pos.push ⟨← floatTok ctx l[1]!, ← floatTok ctx l[2]!, ← floatTok ctx l[3]!⟩
    i := i + 1
  let mut faces : Array Face := Array.mkEmpty nf
  for f in [0:nf] do
    let l ← line i s!"'f <a> <b> <c>' for face {f}"
    if l[0]! != "f" || l.size != 4 then
      throw s!"malformed line '{joinTokens l}': expected 'f <a> <b> <c>' for face {f}"
    let ctx := s!"face {f}"
    faces := faces.push ⟨← natTok ctx l[1]!, ← natTok ctx l[2]!, ← natTok ctx l[3]!⟩
    i := i + 1
  if i == lines.size then return { name, task, pos, faces, params := none }
  -- parameter section
  let l ← line i "'lambda <value>'"
  if l[0]! != "lambda" || l.size != 2 then
    throw s!"malformed line '{joinTokens l}': expected 'lambda <value>' or 'end'"
  let lambda ← floatTok "lambda" l[1]!
  i := i + 1
  let mut wE := FloatArray.emptyWithCapacity nv
  let mut wR := FloatArray.emptyWithCapacity nv
  let mut lo := FloatArray.emptyWithCapacity nv
  let mut hi := FloatArray.emptyWithCapacity nv
  for v in [0:nv] do
    let l ← line i s!"'w <wE> <wR> <lo> <hi>' for vertex {v}"
    if l[0]! != "w" || l.size != 5 then
      throw s!"malformed line '{joinTokens l}': expected 'w <wE> <wR> <lo> <hi>' for vertex {v}"
    let ctx := s!"parameters of vertex {v}"
    wE := wE.push (← floatTok ctx l[1]!)
    wR := wR.push (← floatTok ctx l[2]!)
    lo := lo.push (← floatTok ctx l[3]!)
    hi := hi.push (← floatTok ctx l[4]!)
    i := i + 1
  let l ← line i "'cones <count>'"
  if l[0]! != "cones" || l.size != 2 then
    throw s!"malformed line '{joinTokens l}': expected 'cones <count>'"
  let k ← natTok "cones" l[1]!
  if k > nv then throw s!"{k} prescribed cones for {nv} vertices"
  i := i + 1
  let mut cones : Array (Nat × Float) := Array.mkEmpty k
  for j in [0:k] do
    let l ← line i s!"'c <vertex> <angle>' for prescribed cone {j}"
    if l[0]! != "c" || l.size != 3 then
      throw s!"malformed line '{joinTokens l}': expected 'c <vertex> <angle>' for prescribed cone {j}"
    cones := cones.push (← natTok "prescribed cone" l[1]!, ← floatTok "prescribed cone" l[2]!)
    i := i + 1
  if i != lines.size then
    throw s!"unexpected line '{joinTokens lines[i]!}': expected 'end'"
  return { name, task, pos, faces, params := some { lambda, wE, wR, lo, hi, cones } }

/-- Reader of token lines with one line of lookahead. Blank lines and lines starting with `#`
are skipped. -/
structure LineReader where
  /-- source -/
  stream : IO.FS.Stream
  /-- a line read but not yet consumed -/
  pending : IO.Ref (Option (Array String))

/-- A reader on a stream. -/
def LineReader.new (stream : IO.FS.Stream) : IO LineReader := do
  return { stream, pending := ← IO.mkRef none }

/-- Next non-blank, non-comment line as tokens; `none` at end of input. -/
partial def LineReader.next (r : LineReader) : IO (Option (Array String)) := do
  match ← r.pending.get with
  | some l =>
    r.pending.set none
    return some l
  | none =>
    let raw ← r.stream.getLine
    if raw.isEmpty then return none
    let toks := tokens raw
    if toks.size == 0 || toks[0]!.startsWith "#" then r.next else return some toks

/-- Read the next case. `none` at end of input; otherwise the case name and either the parsed
case or the reason it is malformed. The reader is left after the `end` line of the case, so a
malformed case never affects the cases after it. -/
partial def readCase (r : LineReader) : IO (Option (String × Except String CaseInput)) := do
  match ← r.next with
  | none => return none
  | some first =>
    if first[0]! != "case" || first.size != 2 then
      -- resynchronize: drop everything up to and including the next `end`
      let rec skip : IO Unit := do
        match ← r.next with
        | none => pure ()
        | some l => if l[0]! == "end" then pure () else skip
      if first[0]! != "end" then skip
      return some ("?", .error s!"expected 'case <name>', got '{joinTokens first}'")
    let name := first[1]!
    let mut lines : Array (Array String) := #[]
    let mut closed := false
    let mut problem : Option String := none
    while !closed && problem.isNone do
      match ← r.next with
      | none => problem := some "unexpected end of input: missing 'end'"
      | some l =>
        if l[0]! == "end" then closed := true
        else if l[0]! == "case" then
          r.pending.set (some l)
          problem := some "missing 'end' before the next 'case'"
        else lines := lines.push l
    match problem with
    | some e => return some (name, .error e)
    | none => return some (name, parseBody name lines)

/-! ## Printing -/

/-- A case in the input format. -/
def printCase (c : Case) (task : Option Task := none) (withParams : Bool := true) : String :=
  Id.run do
  let nv := c.mesh.pos.size
  let mut s := s!"case {c.name}\n"
  if let some t := task then s := s ++ s!"task {t.name}\n"
  s := s ++ s!"mesh {nv} {c.mesh.faces.size}\n"
  for p in c.mesh.pos do
    s := s ++ s!"v {floatToHex p.x} {floatToHex p.y} {floatToHex p.z}\n"
  for f in c.mesh.faces do
    s := s ++ s!"f {f.a} {f.b} {f.c}\n"
  if withParams then
    let p := c.params
    s := s ++ s!"lambda {floatToHex p.lambda}\n"
    for v in [0:nv] do
      s := s ++ s!"w {floatToHex (p.wE.get! v)} {floatToHex (p.wR.get! v)} {floatToHex (p.lo.get! v)} {floatToHex (p.hi.get! v)}\n"
    s := s ++ s!"cones {p.cones.size}\n"
    for (v, a) in p.cones do
      s := s ++ s!"c {v} {floatToHex a}\n"
  return s ++ "end\n"

private def header (name : String) (task : Task) (status : String) : String :=
  s!"case {name}\ntask {task.name}\nstatus {status}\n"

/-- Output block of a case that could not be processed. -/
def printError (name : String) (task : Task) (message : String) : String :=
  header name task s!"error {message}" ++ "end\n"

private def meshLine (t : Topology) : String :=
  s!"mesh {t.nv} {t.faces.size} {t.numEdges} {t.numBoundary} {t.chi}\n"

/-- Output block of the `operators` task. -/
def printOperators (name : String) (ops : Operators) : String := Id.run do
  let nv := ops.nv
  let mut s := header name .operators "ok" ++ meshLine ops.topo
  s := s ++ s!"area {floatToHex ops.totalArea}\n"
  for v in [0:nv] do
    let b := if ops.topo.isBoundary[v]! then 1 else 0
    s := s ++ s!"o {v} {b} {floatToHex (ops.curvature.get! v)} {floatToHex (ops.mass.get! v)}\n"
  for f in [0:ops.faceGeom.size] do
    let g := ops.faceGeom[f]!
    s := s ++ s!"t {f} {floatToHex g.area} {floatToHex g.angA} {floatToHex g.angB} {floatToHex g.angC} {floatToHex g.cotA} {floatToHex g.cotB} {floatToHex g.cotC}\n"
  let mut count := 0
  let mut body := ""
  for i in [0:nv] do
    body := body ++ s!"l {i} {i} {floatToHex (ops.L i i)}\n"
    count := count + 1
    for j in ops.topo.neighbors[i]! do
      if i < j then
        body := body ++ s!"l {i} {j} {floatToHex (ops.L i j)}\n"
        count := count + 1
  return s ++ s!"laplacian {count}\n" ++ body ++ "end\n"

/-- Output block of the `solve` task. -/
def printSolve (name : String) (ops : Operators) (r : SolveResult) : String := Id.run do
  match r.status with
  | .infeasible msg => return header name .solve s!"infeasible {msg}" ++ "end\n"
  | _ => pure ()
  let (status, full) := match r.status with
    | .converged => ("converged", true)
    | .notConverged msg => (s!"not_converged {msg}", false)
    | .infeasible msg => (s!"infeasible {msg}", false)
  let c := r.cert
  let mut s := header name .solve status ++ meshLine ops.topo
  s := s ++ s!"interior {r.nInterior}\n"
  s := s ++ s!"closed {if r.closed then 1 else 0}\n"
  s := s ++ s!"iterations {r.itersOne.foldl (· + ·) 0} {r.itersTwo} {r.stop.name}\n"
  s := s ++ s!"objective {floatToHex c.objective} {floatToHex c.distortion} {floatToHex c.coneCost}\n"
  s := s ++ s!"residuals {floatToHex c.resPrimal} {floatToHex c.resAdjoint} {floatToHex c.resInclusion} {floatToHex c.resBounds} {floatToHex c.gap}\n"
  s := s ++ s!"gauss_bonnet {floatToHex c.sumMu} {floatToHex c.sumOmega}\n"
  s := s ++ s!"phi_shift {floatToHex c.shiftLo} {floatToHex c.shiftHi}\n"
  if full then
    for v in [0:ops.nv] do
      s := s ++ s!"s {v} {floatToHex (r.u.get! v)} {floatToHex (r.phi.get! v)} {floatToHex (r.mu.get! v)}\n"
  return s ++ "end\n"

/-! ## Running a case -/

/-- Process one parsed case and return its output block. -/
def runCase (c : CaseInput) (task : Task) (opts : Options := {}) : String :=
  let result : Except String String := do
    let topo ← buildTopology c.pos.size c.faces
    let ops ← buildOperators topo c.pos
    match task with
    | .operators => return printOperators c.name ops
    | .solve =>
      match c.params with
      | none => throw "task solve needs the parameter section ('lambda', 'w' lines, 'cones')"
      | some p =>
        let r ← solve ops p opts
        return printSolve c.name ops r
  match result with
  | .ok s => s
  | .error e => printError c.name task e

/-- Read cases from `input` until end of input, write one output block per case to `output`,
flushing after each. `fixed` is the task of the `operators`/`solve` commands; in `batch` mode
(`fixed = none`) every case must carry its own `task` line. Returns the number of cases. -/
partial def runStream (input output : IO.FS.Stream) (fixed : Option Task) : IO Nat := do
  let reader ← LineReader.new input
  let rec loop (count : Nat) : IO Nat := do
    match ← readCase reader with
    | none => return count
    | some (name, parsed) =>
      let block := match parsed with
        | .error e => printError name (fixed.getD .solve) e
        | .ok c =>
          match fixed, c.task with
          | some t, none => runCase c t
          | some t, some t' =>
            if t == t' then runCase c t
            else printError name t s!"the case asks for task {t'.name} but the command is {t.name}"
          | none, some t => runCase c t
          | none, none => printError name .solve "batch mode needs a 'task <operators|solve>' line after 'case'"
      output.putStr block
      output.flush
      loop (count + 1)
  loop 0

end Ddg
