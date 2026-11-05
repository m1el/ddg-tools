import DdgOracle.Cuts.Format
import DdgOracle.IntegerCoords.Formulas

/-!
# Text format of the `ic` task (integer coordinates, sixth paper)

Specified in `README.md` ("Integer coordinates (sixth paper)"). Everything is a decimal integer; the
outputs are compared with exact equality. Every number printed is computed by the definitions of
`DdgOracle/IntegerCoords/Formulas.lean` (the ones `lean/ddg/DDG/NormalCoords` and the bridge
`lean/proofs/DdgProofs/Oracle/NormalCoords.lean` prove theorems about); this file only parses, loops and prints.

A case is `case NAME`, `task T`, task-specific lines, `end`. The output repeats `case`/`task`, gives
`status ok` (or `status error MESSAGE` for a malformed case), the task's lines, and `end`. A line of
the input that is well formed but describes something invalid (an invalid face, a crossing index
out of range) produces an output line that says so; it does not make the case fail.
-/

namespace Ddg.IntegerCoords

/-- One case: its name, its task and its body lines (tokens), kept as read. -/
structure ICCase where
  name : String
  task : String
  lines : Array (Array String)

/-- The tasks. -/
def tasks : List String :=
  ["corner", "flip", "splitedge", "trace", "tracecurve", "identify", "subdivide", "regions", "split"]

private def joinTokens (l : Array String) : String := " ".intercalate l.toList

/-- A case in the input format. -/
def printCase (c : ICCase) : String := Id.run do
  let mut s := s!"case {c.name}\ntask {c.task}\n"
  for l in c.lines do s := s ++ joinTokens l ++ "\n"
  return s ++ "end\n"

/-- An error block. -/
def printError (name task msg : String) : String :=
  s!"case {name}\ntask {task}\nstatus error {msg}\nend\n"

private def intTok (l : Array String) (i : Nat) : Except String Int :=
  match l[i]? with
  | none => throw s!"line '{joinTokens l}': missing token {i}"
  | some t => match t.toInt? with
    | some x => pure x
    | none => throw s!"line '{joinTokens l}': expected an integer, got '{t}'"

private def natTok (l : Array String) (i : Nat) : Except String Nat := do
  let x ← intTok l i
  if x < 0 then throw s!"line '{joinTokens l}': token {i} must be nonnegative"
  return x.toNat

private def expectSize (l : Array String) (n : Nat) : Except String Unit :=
  if l.size != n then throw s!"line '{joinTokens l}': expected {n - 1} numbers after '{l[0]!}'"
  else pure ()

/-! ## Printing helpers -/

/-- A region as tokens. -/
def Region.toStr : Region → String
  | .corner k t => s!"corner {k} {t}"
  | .middle s => s!"middle {s}"

/-- A boundary point as tokens. -/
def Pt.toStr : Pt → String
  | .corner k => s!"corner {k}"
  | .side k p => s!"side {k} {p}"

private def b2s (b : Bool) : String := if b then "1" else "0"

/-- Parse a region from tokens starting at `i`. -/
def parseRegion (l : Array String) (i : Nat) : Except String Region := do
  match l[i]? with
  | some "corner" => return .corner (← natTok l (i + 1)) (← natTok l (i + 2))
  | some "middle" => return .middle (← natTok l (i + 1))
  | _ => throw s!"line '{joinTokens l}': expected 'corner K T' or 'middle S'"

/-! ## Tracing a curve through a small Δ-complex (Alg. 1 repeated) -/

/-- A small Δ-complex: vertex labels per corner, twin per halfedge (`-1`: boundary), `n` per
halfedge (equal on twins). Halfedge `3f + k` runs from corner `k` to corner `k + 1` of face `f`. -/
structure Complex where
  vert : Array Nat
  twin : Array Int
  n : Array Int

/-- Number of faces. -/
def Complex.nf (c : Complex) : Nat := c.vert.size / 3

/-- Check the gluing and the faces. -/
def Complex.check (c : Complex) : Except String Unit := do
  let nh := c.vert.size
  if c.twin.size != nh || c.n.size != nh then throw "complex: every face needs 'f', 'tw' and 'nf'"
  for h in [0:nh] do
    let t := c.twin[h]!
    if t >= 0 then
      let tn := t.toNat
      if tn >= nh then throw s!"halfedge {h}: twin {t} does not exist"
      if c.twin[tn]! != (h : Int) then throw s!"halfedge {h}: twin {t} does not point back"
      if tn == h then throw s!"halfedge {h} is its own twin"
      if c.n[tn]! != c.n[h]! then throw s!"halfedge {h}: n = {c.n[h]!} but its twin {t} has {c.n[tn]!}"
      let f := h / 3
      let g := tn / 3
      let tailH := c.vert[h]!
      let headH := c.vert[3 * f + (h % 3 + 1) % 3]!
      let tailT := c.vert[tn]!
      let headT := c.vert[3 * g + (tn % 3 + 1) % 3]!
      if tailH != headT || headH != tailT then
        throw s!"halfedge {h} ({tailH} -> {headH}) and its twin {t} ({tailT} -> {headT}) disagree"
    else if t != -1 then throw s!"halfedge {h}: twin must be a halfedge or -1"
  for f in [0:c.nf] do
    if !faceValid c.n[3 * f]! c.n[3 * f + 1]! c.n[3 * f + 2]! then
      throw s!"face {f} is not valid (an odd numerator of Eq. 1)"

/-- Trace from the crossing `(h, p)` (Alg. 1, TraceFrom): the crossings visited, starting with
`(h, p)`, and the end: `ok (face, corner, vertex)` or an error message. A crossing visited twice (in
either direction) is an error (SPEC I18), as is a walk longer than `bound` steps. -/
def traceCurve (c : Complex) (h0 : Nat) (p0 : Int) (bound : Nat) :
    Array (Nat × Int) × Except String (Nat × Nat × Nat) := Id.run do
  let nh := c.vert.size
  let mut out : Array (Nat × Int) := #[]
  let mut seen : Array (Nat × Int) := #[]
  let mut h := h0
  let mut p := p0
  for _ in [0:bound + 1] do
    if h >= nh then return (out, .error s!"halfedge {h} does not exist")
    let nij := c.n[h]!
    if p < 0 || p >= nPlus nij then
      return (out, .error s!"crossing ({h}, {p}) out of range: n = {nij}")
    out := out.push (h, p)
    let t := c.twin[h]!
    if t < 0 then return (out, .error s!"crossing ({h}, {p}) is on a boundary halfedge")
    let tn := t.toNat
    -- the crossing, independent of the direction
    let key := if h < tn then (h, p) else (tn, nPlus nij - p - 1)
    if seen.contains key then
      return (out, .error s!"crossing ({h}, {p}) visited twice: the curve is closed")
    seen := seen.push key
    let f := tn / 3
    let s := tn % 3
    let hik := 3 * f + (s + 1) % 3
    let hkj := 3 * f + (s + 2) % 3
    match traceStep nij c.n[hik]! c.n[hkj]! p with
    | .toIK q => h := hik; p := q
    | .toKJ q => h := hkj; p := q
    | .vertex =>
      let corner := (s + 2) % 3
      return (out, .ok (f, corner, c.vert[3 * f + corner]!))
  return (out, .error s!"no vertex reached after {bound} steps")

/-! ## Running one case -/

private def faceTriple (l : Array String) : Except String (Int × Int × Int) := do
  return (← intTok l 1, ← intTok l 2, ← intTok l 3)

private def runCorner (c : ICCase) : Except String String := do
  let mut s := ""
  for l in c.lines do
    if l[0]! != "face" then throw s!"line '{joinTokens l}': expected 'face N0 N1 N2'"
    expectSize l 4
    let (n0, n1, n2) ← faceTriple l
    let sd := side n0 n1 n2
    let nums := [0, 1, 2].map fun k => cornerNumerator (sd (next k)) (sd k) (sd (prev k))
    let cs := [0, 1, 2].map fun k => match cornerCoord? (sd (next k)) (sd k) (sd (prev k)) with
      | some x => toString x
      | none => "odd"
    let es := [0, 1, 2].map fun k => faceEmanating n0 n1 n2 k
    s := s ++ s!"face {n0} {n1} {n2} {b2s (faceValid n0 n1 n2)} {" ".intercalate (nums.map toString)} {" ".intercalate cs} {" ".intercalate (es.map toString)}\n"
  return s

private def flipOne (q : Int × Int × Int × Int × Int) (rb : Option (Int × Nat × Int × Nat))
    (again : Bool) : String := Id.run do
  let (nij, njk, nki, nil, nlj) := q
  let v1 := faceValid nij njk nki
  let v2 := faceValid nij nil nlj
  let mut s := s!"quad {nij} {njk} {nki} {nil} {nlj}\n"
  if !(v1 && v2) then return s ++ s!"faces {b2s v1} {b2s v2} invalid\n"
  s := s ++ s!"faces {b2s v1} {b2s v2}\n"
  let d := flipDoubled nij njk nki nil nlj
  let nkl := flipNormal nij njk nki nil nlj
  s := s ++ s!"flip {d} {nkl}\n"
  s := s ++ s!"newfaces {b2s (faceValid nki nil nkl)} {b2s (faceValid nlj njk nkl)}\n"
  let ek := emanating nil nki nkl
  let el := emanating njk nlj nkl
  s := s ++ s!"emanating {ek} {el}\n"
  if let some (rki, dk, rlj, dl) := rb then
    s := s ++ s!"roundabout {roundaboutUpdate rki ek (nMinus nki) dk} {roundaboutUpdate rlj el (nMinus nlj) dl}\n"
  if again then
    s := s ++ s!"again {flipNormal nkl nki nil nlj njk}\n"
  return s

private def runFlip (c : ICCase) : Except String String := do
  -- each `quad` line starts a block; `roundabout` and `again` belong to the last `quad`
  let mut blocks : Array ((Int × Int × Int × Int × Int) × Option (Int × Nat × Int × Nat) × Bool) := #[]
  for l in c.lines do
    if l[0]! == "quad" then
      expectSize l 6
      blocks := blocks.push ((← intTok l 1, ← intTok l 2, ← intTok l 3, ← intTok l 4, ← intTok l 5), none, false)
    else if blocks.size == 0 then throw s!"line '{joinTokens l}': a 'quad' line must come first"
    else if l[0]! == "roundabout" then
      expectSize l 5
      let rb := (← intTok l 1, ← natTok l 2, ← intTok l 3, ← natTok l 4)
      blocks := blocks.modify (blocks.size - 1) fun (q, _, a) => (q, some rb, a)
    else if l[0]! == "again" then
      expectSize l 1
      blocks := blocks.modify (blocks.size - 1) fun (q, r, _) => (q, r, true)
    else throw s!"line '{joinTokens l}': expected 'quad', 'roundabout' or 'again'"
  return String.join (blocks.toList.map fun (q, r, a) => flipOne q r a)

private def runSplitEdge (c : ICCase) : Except String String := do
  let mut s := ""
  for l in c.lines do
    if l[0]! != "edge" then throw s!"line '{joinTokens l}': expected 'edge NKI NJK'"
    expectSize l 3
    let a ← intTok l 1
    let b ← intTok l 2
    s := s ++ s!"edge {a} {b} {splitEdgeNormal a b}\n"
  return s

private def runTrace (c : ICCase) : Except String String := do
  let mut s := ""
  for l in c.lines do
    if l[0]! != "cross" then throw s!"line '{joinTokens l}': expected 'cross NIJ NIK NKJ P'"
    expectSize l 5
    let nij ← intTok l 1
    let nik ← intTok l 2
    let nkj ← intTok l 3
    let p ← intTok l 4
    let pre := s!"cross {nij} {nik} {nkj} {p}"
    if !faceValid nij nik nkj then s := s ++ s!"{pre} error face not valid\n"
    else if p < 0 || p >= nPlus nij then s := s ++ s!"{pre} error crossing out of range\n"
    else
      let r := match traceStep nij nik nkj p with
        | .toIK q => s!"ik {q}"
        | .toKJ q => s!"kj {q}"
        | .vertex => "vertex"
      s := s ++ s!"{pre} step {r}\n"
  return s

private def set3 {α : Type} (a : Array α) (f : Nat) (x y z : α) : Array α :=
  ((a.set! (3 * f) x).set! (3 * f + 1) y).set! (3 * f + 2) z

private def runTraceCurve (c : ICCase) : Except String String := do
  let mut nf := 0
  let mut vert : Array Nat := #[]
  let mut twin : Array Int := #[]
  let mut nn : Array Int := #[]
  let mut starts : Array (Nat × Int) := #[]
  let mut bound : Option Nat := none
  for l in c.lines do
    let tag := l[0]!
    if tag == "faces" then
      expectSize l 2
      nf ← natTok l 1
      vert := Array.replicate (3 * nf) 0
      twin := Array.replicate (3 * nf) (-2)
      nn := Array.replicate (3 * nf) 0
    else if tag == "f" || tag == "tw" || tag == "nf" then
      expectSize l 5
      let f ← natTok l 1
      if f >= nf then throw s!"line '{joinTokens l}': face out of range (declare 'faces NF' first)"
      if tag == "f" then vert := set3 vert f (← natTok l 2) (← natTok l 3) (← natTok l 4)
      else if tag == "tw" then twin := set3 twin f (← intTok l 2) (← intTok l 3) (← intTok l 4)
      else nn := set3 nn f (← intTok l 2) (← intTok l 3) (← intTok l 4)
    else if tag == "start" then
      expectSize l 3
      starts := starts.push (← natTok l 1, ← intTok l 2)
    else if tag == "bound" then
      expectSize l 2
      bound := some (← natTok l 1)
    else throw s!"line '{joinTokens l}': expected faces, f, tw, nf, start or bound"
  if twin.contains (-2) then throw "every face needs a 'tw' line"
  let cx : Complex := { vert, twin, n := nn }
  cx.check
  let b := bound.getD (nn.foldl (fun a x => a + (nPlus x).toNat) 0 + nn.size)
  let mut s := ""
  for i in [0:starts.size] do
    let (h, p) := starts[i]!
    s := s ++ s!"curve {i} {h} {p}\n"
    let (xs, r) := traceCurve cx h p b
    for (hh, pp) in xs do s := s ++ s!"x {hh} {pp}\n"
    match r with
    | .ok (f, k, v) => s := s ++ s!"stop vertex {v} {f} {k} {xs.size}\n"
    | .error e => s := s ++ s!"stop error {e}\n"
  return s

private def runIdentify (c : ICCase) : Except String String := do
  let mut s := ""
  for l in c.lines do
    if l[0]! != "walk" then throw s!"line '{joinTokens l}': expected 'walk RAI NAI NIJ NJA P DEG'"
    expectSize l 7
    let r ← intTok l 1
    let nai ← intTok l 2
    let nij ← intTok l 3
    let nja ← intTok l 4
    let p ← intTok l 5
    let deg ← natTok l 6
    s := s ++ s!"walk {r} {nai} {nij} {nja} {p} {deg} index {identifyIndex r nai nij nja p deg}\n"
  return s

private def runFaces (c : ICCase) : Except String String := do
  let mut s := ""
  for l in c.lines do
    if l[0]! != "face" then throw s!"line '{joinTokens l}': expected 'face N0 N1 N2 ...'"
    if c.task != "split" then expectSize l 4
    else if l.size != 4 && l.size != 6 && l.size != 7 then
      throw s!"line '{joinTokens l}': expected 'face N0 N1 N2 [corner K T | middle S]'"
    let (n0, n1, n2) ← faceTriple l
    let pre := s!"face {n0} {n1} {n2}"
    match decode? n0 n1 n2 with
    | none => s := s ++ s!"{pre} invalid\n"
    | some d =>
      if c.task == "subdivide" then
        s := s ++ s!"{pre} decoded {d.c 0} {d.c 1} {d.c 2} {d.kappa} {d.e}\n"
        for k in [0:3] do
          for p in [0:d.n k] do
            s := s ++ s!"arc {k} {p} {(partner d.n d.c k p).toStr}\n"
      else if c.task == "regions" then
        let rs := regions d
        s := s ++ s!"{pre} regions {rs.length}\n"
        for r in rs do s := s ++ s!"r {r.toStr}\n"
      else
        let rs ← if l.size == 4 then pure (regions d) else do pure [← parseRegion l 4]
        for r in rs do
          match faceSplit? n0 n1 n2 r with
          | some (a, b, e) => s := s ++ s!"split {n0} {n1} {n2} {r.toStr} {a} {b} {e}\n"
          | none => s := s ++ s!"split {n0} {n1} {n2} {r.toStr} none\n"
  return s

/-- The output lines of one task. -/
def runBody (c : ICCase) : Except String String :=
  match c.task with
  | "corner" => runCorner c
  | "flip" => runFlip c
  | "splitedge" => runSplitEdge c
  | "trace" => runTrace c
  | "tracecurve" => runTraceCurve c
  | "identify" => runIdentify c
  | "subdivide" | "regions" | "split" => runFaces c
  | t => throw s!"unknown task '{t}' (expected one of {", ".intercalate tasks})"

/-- Run one case and print its output block. -/
def runCase (c : ICCase) : String :=
  match runBody c with
  | .ok body => s!"case {c.name}\ntask {c.task}\nstatus ok\n" ++ body ++ "end\n"
  | .error e => printError c.name c.task e

/-- Parse the lines of one case (after `case NAME`). -/
def parseBody (name : String) (lines : Array (Array String)) : Except String ICCase := do
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
          match parseBody name lines with
          | .error e => printError name "?" e
          | .ok c => runCase c
      output.putStr block
      output.flush
      loop (count + 1)
  loop 0

end Ddg.IntegerCoords
