import DdgOracle.IntrinsicSimplification.Star
import DdgOracle.Cuts.Format

/-!
# Text format of the `iem` task (intrinsic simplification, seventh paper)

Specified in `README.md` ("Intrinsic simplification"). Floats are 16-digit hexadecimal bit patterns
(text format 1 of this package), labels and counts are decimal. Every number printed is computed in
`Star.lean` from the formulas of `DdgOracle/Formulas/IntrinsicSimplification.lean`; this file only
parses, loops and prints.
-/

namespace Ddg.IEM

open Ddg

/-- One case as read: name, task, body lines (tokens). -/
structure ICase where
  name : String
  task : String
  lines : Array (Array String)
  deriving Inhabited

/-- The tasks. -/
def tasks : List String := ["step", "flatten", "weights", "removal"]

private def joinTokens (l : Array String) : String := " ".intercalate l.toList

/-- A case in the input format. -/
def printCase (c : ICase) : String := Id.run do
  let mut s := s!"case {c.name}\ntask {c.task}\n"
  for l in c.lines do s := s ++ joinTokens l ++ "\n"
  return s ++ "end\n"

/-- An error block. -/
def printError (name task msg : String) : String :=
  s!"case {name}\ntask {task}\nstatus error {msg}\nend\n"

private def hexTok (l : Array String) (i : Nat) : Except String Float :=
  match l[i]? with
  | none => throw s!"line '{joinTokens l}': missing token {i}"
  | some t => match hexToFloat? t with
    | some x => pure x
    | none => throw s!"line '{joinTokens l}': expected a 16-digit hex float, got '{t}'"

private def natTok (l : Array String) (i : Nat) : Except String Nat :=
  match l[i]? with
  | none => throw s!"line '{joinTokens l}': missing token {i}"
  | some t => match t.toNat? with
    | some x => pure x
    | none => throw s!"line '{joinTokens l}': expected a nonnegative integer, got '{t}'"

private def expectSize (l : Array String) (n : Nat) : Except String Unit :=
  if l.size != n then throw s!"line '{joinTokens l}': expected {n - 1} tokens after '{l[0]!}'"
  else pure ()

private def h (x : Float) : String := floatToHex x

/-- The parsed body of a star task. -/
structure Body where
  target : Float := 2.0 * floatPi
  faces : Array SFace := #[]
  u : Option Float := none
  nbrs : Array Nbr := #[]
  chans : Array Channel := #[]
  memory : Bool := true
  ki : Option Float := none
  tol : Float := flatTol
  dks : Array (Nat × Float) := #[]

/-- Parse the body lines. -/
def parseBody (c : ICase) : Except String Body := do
  let mut b : Body := {}
  for l in c.lines do
    match l[0]! with
    | "target" =>
      expectSize l 2
      match l[1]! with
      | "interior" => b := { b with target := 2.0 * floatPi }
      | "boundary" => b := { b with target := floatPi }
      | _ => b := { b with target := ← hexTok l 1 }
    | "face" =>
      expectSize l 7
      let v := #[← natTok l 1, ← natTok l 2, ← natTok l 3]
      let s := #[← hexTok l 4, ← hexTok l 5, ← hexTok l 6]
      if !(s.all (· > 0.0)) then throw s!"line '{joinTokens l}': sides must be positive"
      if !(v.any (· == 0)) then throw s!"line '{joinTokens l}': no corner is the vertex 0"
      b := { b with faces := b.faces.push { v, l := s } }
    | "u" => expectSize l 2; b := { b with u := some (← hexTok l 1) }
    | "nbr" =>
      expectSize l 5
      let j ← natTok l 1
      if j == 0 then throw s!"line '{joinTokens l}': a neighbour label must be ≥ 1"
      let pij ← hexTok l 2
      let pji ← hexTok l 3
      let len ← hexTok l 4
      let n : Nbr := { j := j, phiIJ := pij, phiJI := pji, len := len }
      b := { b with nbrs := b.nbrs.push n }
    | "channel" =>
      expectSize l 5
      let w ← hexTok l 1
      let mi ← hexTok l 2
      let tx ← hexTok l 3
      let ty ← hexTok l 4
      let ch : Channel := { w := w, mi := mi, ti := (tx, ty), nb := #[] }
      b := { b with chans := b.chans.push ch }
    | "mass" =>
      expectSize l 5
      if b.chans.size == 0 then throw s!"line '{joinTokens l}': 'mass' before any 'channel'"
      let j ← natTok l 1
      let k := b.chans.size - 1
      let ch := b.chans[k]!
      let pos := ch.nb.size
      match b.nbrs[pos]? with
      | none => throw s!"line '{joinTokens l}': more 'mass' lines than neighbours"
      | some n =>
        if n.j != j then
          throw s!"line '{joinTokens l}': expected the neighbour {n.j} (the order of the 'nbr' lines)"
      let mj ← hexTok l 2
      let tx ← hexTok l 3
      let ty ← hexTok l 4
      let ch' := { ch with nb := ch.nb.push (mj, (tx, ty)) }
      b := { b with chans := b.chans.set! k ch' }
    | "memory" =>
      expectSize l 2
      match l[1]! with
      | "0" => b := { b with memory := false }
      | "1" => b := { b with memory := true }
      | _ => throw s!"line '{joinTokens l}': expected 'memory 0' or 'memory 1'"
    | "ki" => expectSize l 2; b := { b with ki := some (← hexTok l 1) }
    | "tol" => expectSize l 2; b := { b with tol := ← hexTok l 1 }
    | "dk" => expectSize l 3; b := { b with dks := b.dks.push (← natTok l 1, ← hexTok l 2) }
    | t => throw s!"unknown line tag '{t}'"
  return b

private def needFaces (b : Body) : Except String Unit :=
  if b.faces.size == 0 then throw "no 'face' lines" else pure ()

/-- Output lines of `step`. -/
def runStep (b : Body) : Except String String := do
  needFaces b
  let some u := b.u | throw "missing 'u U'"
  let th := thetaI b.faces u
  let hs := hessian b.faces u
  let next := iemNewtonStep u b.target th hs
  let mut s := s!"step {h u} {h th} {h hs} {h next}\n"
  s := s ++ s!"paperhess {h (paperHessian b.faces u)}\n"
  for k in [0:b.faces.size] do
    let L := b.faces[k]!.sides u
    s := s ++ s!"face {k} {h L[0]!} {h L[1]!} {h L[2]!} {if isValid L then 1 else 0}\n"
  return s

/-- Output lines of `flatten`. -/
def runFlatten (b : Body) : Except String String := do
  needFaces b
  let r := flatten b.faces b.target
  let mut s := s!"theta0 {h (thetaI b.faces 0.0)}\n"
  s := s ++ s!"solve {r.status.toStr} {r.iters} {h r.u} {h r.theta}\n"
  s := s ++ s!"hess {h (hessian b.faces r.u)}\n"
  for k in [0:b.faces.size] do
    let L := b.faces[k]!.sides r.u
    s := s ++ s!"face {k} {h L[0]!} {h L[1]!} {h L[2]!} {if isValid L then 1 else 0}\n"
  return s

/-- Output lines of `weights`. -/
def runWeights (b : Body) : Except String String := do
  let some ki := b.ki | throw "missing 'ki K'"
  if b.dks.size == 0 then throw "no 'dk' lines"
  let sumAbs := b.dks.foldl (fun s x => s + Float.abs x.2) 0.0
  let count := Float.ofNat b.dks.size
  let mut s := ""
  for (j, d) in b.dks do
    s := s ++ s!"w {j} {h (iemWeight d sumAbs ki b.tol count)}\n"
  return s

/-- Output lines of `removal`. -/
def runRemoval (b : Body) : Except String String := do
  needFaces b
  if b.nbrs.size == 0 then throw "no 'nbr' lines"
  -- every label ≠ 0 of the faces must be a neighbour, each once
  for f in b.faces do
    for v in f.v do
      if v != 0 && !(b.nbrs.any (·.j == v)) then throw s!"label {v} has no 'nbr' line"
  for k in [0:b.nbrs.size] do
    for k' in [0:k] do
      if b.nbrs[k]!.j == b.nbrs[k']!.j then throw s!"neighbour {b.nbrs[k]!.j} listed twice"
  for c in b.chans do
    if c.nb.size != b.nbrs.size then throw "every channel needs one 'mass' line per neighbour"
  let r := removal b.faces b.target b.nbrs b.chans b.memory
  let mut s := s!"flatten {r.flat.status.toStr} {r.flat.iters} {h r.flat.u}\n"
  if r.flat.status != .ok then
    return s ++ s!"cost {h posInf}\n"
  s := s ++ s!"ki {h r.Ki}\n"
  for k in [0:b.nbrs.size] do
    let (d, a, len) := r.nbr[k]!
    s := s ++ s!"nbr {b.nbrs[k]!.j} {h d} {h a} {h len}\n"
  for p in [0:r.chan.size] do
    for k in [0:b.nbrs.size] do
      let (m, t) := r.chan[p]![k]!
      s := s ++ s!"new {p} {b.nbrs[k]!.j} {h m} {h t.1} {h t.2}\n"
  s := s ++ s!"cost {h r.cost}\nmemoryless {h r.memoryless}\n"
  return s

/-- The output lines of one case. -/
def runBody (c : ICase) : Except String String := do
  let b ← parseBody c
  match c.task with
  | "step" => runStep b
  | "flatten" => runFlatten b
  | "weights" => runWeights b
  | "removal" => runRemoval b
  | t => throw s!"unknown task '{t}' (expected one of {", ".intercalate tasks})"

/-- Run one case and print its output block. -/
def runCase (c : ICase) : String :=
  match runBody c with
  | .ok body => s!"case {c.name}\ntask {c.task}\nstatus ok\n" ++ body ++ "end\n"
  | .error e => printError c.name c.task e

/-- Parse the lines of one case (after `case NAME`). -/
def parseCase (name : String) (lines : Array (Array String)) : Except String ICase := do
  match lines[0]? with
  | some l =>
    if l[0]! != "task" || l.size != 2 then throw s!"expected 'task T', got '{joinTokens l}'"
    return { name, task := l[1]!, lines := lines.extract 1 lines.size }
  | none => throw "missing 'task T'"

/-- Read cases until the end of `input`, write one block per case, flush after each. -/
partial def runStream (input output : IO.FS.Stream) : IO Nat := do
  let reader ← LineReader.new input
  let rec loop (count : Nat) : IO Nat := do
    match ← Ddg.Cuts.readCaseLines reader with
    | none => return count
    | some (name, parsed) =>
      let block := match parsed with
        | .error e => printError name "?" e
        | .ok lines =>
          match parseCase name lines with
          | .error e => printError name "?" e
          | .ok c => runCase c
      output.putStr block
      output.flush
      loop (count + 1)
  loop 0

end Ddg.IEM
