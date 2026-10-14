import DdgOracle.Cuts.Format
import DdgOracle.Intrinsic.Trace

/-!
# Text format of the `intrinsic` task (single operations)

Specified in `README.md` ("Intrinsic triangulations"). Same conventions as the other tasks.
-/

namespace Ddg.Intrinsic

open Ddg Ddg.Cuts

/-- One case. -/
structure ICase where
  /-- identifier -/
  name : String
  /-- `flip`, `delaunay`, `trace`, `edgetrace` or `circumcenter` -/
  task : String
  /-- number of vertices -/
  nv : Nat
  /-- faces -/
  faces : Array Face
  /-- edges `(a, b, length)`, `a < b` -/
  edges : Array (Nat × Nat × Float)
  /-- `flip A B` (task flip) -/
  flip : Option (Nat × Nat) := none
  /-- flip the new edge a second time -/
  again : Bool := false
  /-- `preflip A B` lines (task edgetrace): flips applied first -/
  preflips : Array (Nat × Nat) := #[]
  /-- Delaunay flipping first (task edgetrace) -/
  predelaunay : Bool := false
  /-- tolerance of the Delaunay test -/
  tol : Float := 1e-10
  /-- flips per edge allowed -/
  limit : Nat := 64
  /-- criterion -/
  crit : Criterion := .angleSum
  /-- `trace V PHI DIST` lines -/
  traces : Array (Nat × Float × Float) := #[]
  /-- faces for `circumcenter` -/
  circ : Array Nat := #[]
  /-- `circumcenter all` -/
  circAll : Bool := false

private def natTok (ctx tok : String) : Except String Nat :=
  match tok.toNat? with
  | some n => pure n
  | none => throw s!"{ctx}: expected a decimal integer, got '{tok}'"

private def floatTok (ctx tok : String) : Except String Float :=
  match hexToFloat? tok with
  | some x => pure x
  | none => throw s!"{ctx}: expected 16 hex digits (IEEE-754 bit pattern), got '{tok}'"

private def joinTokens (l : Array String) : String := " ".intercalate l.toList

/-- Parse the lines of one case. -/
def parseBody (name : String) (lines : Array (Array String)) : Except String ICase := do
  let mut i := 0
  let line := fun (k : Nat) (what : String) =>
    if h : k < lines.size then pure (f := Except String) lines[k]
    else throw s!"missing line: expected {what}"
  let l ← line i "'task <flip|delaunay|trace|edgetrace|circumcenter>'"
  if l[0]! != "task" || l.size != 2 then throw s!"malformed line '{joinTokens l}': expected 'task <name>'"
  let task := l[1]!
  if !["flip", "delaunay", "trace", "edgetrace", "circumcenter"].contains task then
    throw s!"unknown task '{task}'"
  i := i + 1
  let l ← line i "'mesh <nv> <nf> <ne>'"
  if l[0]! != "mesh" || l.size != 4 then throw s!"malformed line '{joinTokens l}': expected 'mesh <nv> <nf> <ne>'"
  let nv ← natTok "mesh" l[1]!
  let nf ← natTok "mesh" l[2]!
  let ne ← natTok "mesh" l[3]!
  if nv > maxCutVertices then throw s!"mesh has {nv} vertices; the dense oracle accepts at most {maxCutVertices}"
  i := i + 1
  let mut faces : Array Face := Array.mkEmpty nf
  for f in [0:nf] do
    let l ← line i s!"'f <a> <b> <c>' for face {f}"
    if l[0]! != "f" || l.size != 4 then throw s!"malformed line '{joinTokens l}': expected 'f <a> <b> <c>'"
    faces := faces.push ⟨← natTok "face" l[1]!, ← natTok "face" l[2]!, ← natTok "face" l[3]!⟩
    i := i + 1
  let mut edges : Array (Nat × Nat × Float) := Array.mkEmpty ne
  for e in [0:ne] do
    let l ← line i s!"'e <a> <b> <length>' for edge {e}"
    if l[0]! != "e" || l.size != 4 then throw s!"malformed line '{joinTokens l}': expected 'e <a> <b> <length>'"
    edges := edges.push (← natTok "edge" l[1]!, ← natTok "edge" l[2]!, ← floatTok "edge" l[3]!)
    i := i + 1
  let mut c : ICase := { name, task, nv, faces, edges }
  let mut seenSingle : Array String := #[]
  while i < lines.size do
    let l := lines[i]!
    let key := l[0]!
    let single := fun (k : String) => do
      if seenSingle.contains k then throw s!"line '{k}' occurs twice"
    match key with
    | "flip" =>
      single "flip"
      seenSingle := seenSingle.push "flip"
      if l.size != 3 then throw "flip: expected 'flip <a> <b>'"
      c := { c with flip := some (← natTok "flip" l[1]!, ← natTok "flip" l[2]!) }
    | "again" =>
      single "again"
      seenSingle := seenSingle.push "again"
      if l.size != 2 then throw "again: expected 'again 0|1'"
      c := { c with again := (← natTok "again" l[1]!) != 0 }
    | "preflip" =>
      if l.size != 3 then throw "preflip: expected 'preflip <a> <b>'"
      c := { c with preflips := c.preflips.push (← natTok "preflip" l[1]!, ← natTok "preflip" l[2]!) }
    | "predelaunay" =>
      single "predelaunay"
      seenSingle := seenSingle.push "predelaunay"
      if l.size != 2 then throw "predelaunay: expected 'predelaunay 0|1'"
      c := { c with predelaunay := (← natTok "predelaunay" l[1]!) != 0 }
    | "tolerance" =>
      single "tolerance"
      seenSingle := seenSingle.push "tolerance"
      if l.size != 2 then throw "tolerance: expected one value"
      c := { c with tol := ← floatTok "tolerance" l[1]! }
    | "limit" =>
      single "limit"
      seenSingle := seenSingle.push "limit"
      if l.size != 2 then throw "limit: expected one integer"
      c := { c with limit := ← natTok "limit" l[1]! }
    | "criterion" =>
      single "criterion"
      seenSingle := seenSingle.push "criterion"
      if l.size != 2 || (l[1]! != "angle" && l[1]! != "cotan") then throw "criterion: expected angle | cotan"
      c := { c with crit := if l[1]! == "cotan" then .cotan else .angleSum }
    | "trace" =>
      if l.size != 4 then throw "trace: expected 'trace <vertex> <phi> <distance>'"
      c := { c with traces := c.traces.push (← natTok "trace" l[1]!, ← floatTok "trace" l[2]!, ← floatTok "trace" l[3]!) }
    | "circumcenter" =>
      if l.size != 2 then throw "circumcenter: expected 'circumcenter <face>' or 'circumcenter all'"
      if l[1]! == "all" then c := { c with circAll := true }
      else c := { c with circ := c.circ.push (← natTok "circumcenter" l[1]!) }
    | other => throw s!"unknown line '{other}'"
    i := i + 1
  match task with
  | "flip" => if c.flip.isNone then throw "task flip needs a line 'flip <a> <b>'"
  | "trace" => if c.traces.isEmpty then throw "task trace needs at least one line 'trace <vertex> <phi> <distance>'"
  | "circumcenter" => if c.circ.isEmpty && !c.circAll then throw "task circumcenter needs 'circumcenter <face>' or 'circumcenter all'"
  | _ => pure ()
  return c

/-- The case in the input format. -/
def printCase (c : ICase) : String := Id.run do
  let hx := floatToHex
  let mut s := s!"case {c.name}\ntask {c.task}\nmesh {c.nv} {c.faces.size} {c.edges.size}\n"
  for f in c.faces do s := s ++ s!"f {f.a} {f.b} {f.c}\n"
  for (a, b, l) in c.edges do s := s ++ s!"e {a} {b} {hx l}\n"
  if let some (a, b) := c.flip then s := s ++ s!"flip {a} {b}\n"
  if c.again then s := s ++ "again 1\n"
  for (a, b) in c.preflips do s := s ++ s!"preflip {a} {b}\n"
  if c.predelaunay then s := s ++ "predelaunay 1\n"
  if c.task == "delaunay" || c.predelaunay then
    s := s ++ s!"tolerance {hx c.tol}\nlimit {c.limit}\ncriterion {if c.crit == .cotan then "cotan" else "angle"}\n"
  for (v, phi, d) in c.traces do s := s ++ s!"trace {v} {hx phi} {hx d}\n"
  if c.circAll then s := s ++ "circumcenter all\n"
  for f in c.circ do s := s ++ s!"circumcenter {f}\n"
  return s ++ "end\n"

/-- Build the input triangulation of a case; checks the mesh. -/
def buildInput (c : ICase) : Except String (Tri × Topology) := do
  let topo ← buildTopology c.nv c.faces (connected := false)
  let nv := c.nv
  let mut tab := zeros (nv * nv)
  let mut n := 0
  for (a, b, l) in c.edges do
    if !(a < b && b < nv) then throw s!"edge ({a}, {b}): expected a < b < {nv}"
    if !(l.isFinite && l > 0) then throw s!"edge ({a}, {b}) has a length that is not finite and positive"
    if tab.get! (a * nv + b) != 0 then throw s!"edge ({a}, {b}) is listed twice"
    tab := (tab.set! (a * nv + b) l).set! (b * nv + a) l
    n := n + 1
  if n != topo.numEdges then throw s!"{n} edge lengths given but the faces have {topo.numEdges} edges"
  for fi in [0:c.faces.size] do
    let f := c.faces[fi]!
    for (p, q) in [(f.a, f.b), (f.b, f.c), (f.c, f.a)] do
      if tab.get! (p * nv + q) == 0 then throw s!"face {fi}: the length of edge ({min p q}, {max p q}) is missing"
    let a := tab.get! (f.b * nv + f.c)
    let b := tab.get! (f.c * nv + f.a)
    let cc := tab.get! (f.a * nv + f.b)
    if !(b - a + cc > 0 && a - b + cc > 0 && a - cc + b > 0) then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) violates the strict triangle inequality"
  let t ← initTri nv c.faces tab
  return (t, topo)

private def hx (x : Float) : String := floatToHex x

/-- The corners of a face in the canonical rotation: the one that starts with the smallest
`(vertex, next vertex, third vertex, length, signpost)`. -/
def canonicalFace (t : Tri) (f : Nat) : Array (Nat × Float × Float) := Id.run do
  let hs := t.faceHalfedges f
  let key := fun (r : Nat) => (t.tail[hs[r % 3]!]!, t.tail[hs[(r + 1) % 3]!]!, t.tail[hs[(r + 2) % 3]!]!)
  let mut best := 0
  for r in [1:3] do
    let a := key r
    let b := key best
    let less := a.1 < b.1 || (a.1 == b.1 && (a.2.1 < b.2.1 || (a.2.1 == b.2.1 && a.2.2 < b.2.2))) ||
      (a == b && (t.len[hs[r]!]! < t.len[hs[best]!]! || (t.len[hs[r]!]! == t.len[hs[best]!]! && t.sign[hs[r]!]! < t.sign[hs[best]!]!)))
    if less then best := r
  return Array.ofFn (n := 3) fun k => (t.tail[hs[(best + k.val) % 3]!]!, t.sign[hs[(best + k.val) % 3]!]!, t.len[hs[(best + k.val) % 3]!]!)

/-- Print a state: `state LABEL NF NE` and one `cf` line per face. -/
def printState (label : String) (t : Tri) : String := Id.run do
  let mut s := s!"state {label} {t.nf} {t.numEdges}\n"
  for f in [0:t.nf] do
    let cs := canonicalFace t f
    s := s ++ s!"cf {f} {cs[0]!.1} {cs[1]!.1} {cs[2]!.1} {hx cs[0]!.2.1} {hx cs[1]!.2.1} {hx cs[2]!.2.1} {hx cs[0]!.2.2} {hx cs[1]!.2.2} {hx cs[2]!.2.2}\n"
  return s

/-- A crossing as `U V T` with `U < V` and `T` from `U`. -/
def crossingKey (t : Tri) (c : Nat × Float) : Nat × Nat × Float :=
  let u := t.tail[c.1]!
  let v := t.head c.1
  if u < v then (u, v, c.2) else (v, u, 1 - c.2)

private def printTrace (tag : String) (t : Tri) (r : TraceResult) : String := Id.run do
  let (b0, b1, b2) := r.bary
  let vtx : String := if r.vertex == noTwin then "-1" else toString r.vertex
  let mut s := s!"{tag} {r.stop} {r.face} {hx b0} {hx b1} {hx b2} {hx r.left} {vtx} {if r.clamp then 1 else 0} {r.crossings.size}\n"
  for c in r.crossings do
    let (u, v, tt) := crossingKey t c
    s := s ++ s!"{tag}x {u} {v} {hx tt}\n"
  return s

/-- Largest Delaunay excess (angle sum minus pi) and smallest cotan weight over the interior edges. -/
def delaunayMeasures (t : Tri) : Float × Float := Id.run do
  let mut ex := negInf
  let mut w := posInf
  for h in [0:t.nh] do
    if t.eid[h]! == h && t.hasTwin h then
      let tw := t.twin[h]!
      let a := t.len[h]!
      let b := t.len[t.prev h]!
      let c := t.len[t.next[h]!]!
      let d := t.len[t.next[tw]!]!
      let e := t.len[t.prev tw]!
      let (thK, thL) := oppositeAngles a b c d e
      ex := max ex (delaunayExcess thK thL)
      w := min w (cotanWeightFromLengths a b c d e)
  return (ex, w)

/-- Header lines common to every successful case. -/
private def header (c : ICase) (t : Tri) (topo : Topology) (status : String) : String := Id.run do
  let mut s := s!"case {c.name}\ntask {c.task}\nstatus {status}\n"
  s := s ++ s!"mesh {t.nv} {t.nf} {topo.numEdges} {topo.numBoundary} {topo.chi}\n"
  for v in [0:t.nv] do
    s := s ++ s!"vt {v} {hx (t.theta.get! v)} {if t.isB[v]! then 1 else 0}\n"
  return s

/-- Output block of a case that could not be processed. -/
def printError (name task msg : String) : String := s!"case {name}\ntask {task}\nstatus error {msg}\nend\n"

/-- Apply `preflip` lines: each pair must identify one edge. -/
private def applyFlips (t0 : Tri) (pairs : Array (Nat × Nat)) : Except String Tri := do
  let mut t := t0
  for (a, b) in pairs do
    let h ← match findHalfedge t a b with
      | .ok h => pure h
      | .error _ => findHalfedge t b a
    match flipHalfedge t h with
    | .ok (t', _) => t := t'
    | .error r => throw s!"flip of {a} {b} refused: {r.name}"
  return t

/-- Run one case and print its output block. -/
def runCase (c : ICase) : String :=
  match buildInput c with
  | .error e => printError c.name c.task e
  | .ok (t0, topo) =>
    let res : Except String String := do
      match c.task with
      | "flip" =>
        let (a, b) := c.flip.get!
        let h ← match findHalfedge t0 a b with
          | .ok h => pure h
          | .error _ => findHalfedge t0 b a
        match flipHalfedge t0 h with
        | .error r =>
          return header c t0 topo s!"refused {r.name}" ++ printState "before" t0 ++ "end\n"
        | .ok (t1, info) =>
          let mut s := header c t0 topo "ok" ++ printState "before" t0
          s := s ++ s!"flipinfo {info.i} {info.j} {info.k} {info.l} {" ".intercalate (info.lens.toList.map hx)} {" ".intercalate (info.diamond.toList.map hx)} {hx info.lnew}\n"
          s := s ++ printState "after" t1
          if c.again then
            match flipHalfedge t1 h with
            | .error r => s := s ++ s!"again refused {r.name}\n"
            | .ok (t2, _) => s := s ++ printState "again" t2
          return s ++ "end\n"
      | "delaunay" =>
        let (t1, rep) := delaunayFlip t0 c.crit c.tol c.limit
        let (ex, w) := delaunayMeasures t1
        let mut s := header c t0 topo "ok" ++ printState "input" t0
        s := s ++ s!"report {rep.nFlips} {rep.nRefused} {rep.nTies} {rep.status}\n"
        s := s ++ s!"measures {hx ex} {hx w} {hx t1.totalArea}\n"
        s := s ++ printState "final" t1
        return s ++ "end\n"
      | "trace" =>
        let mut s := header c t0 topo "ok"
        for ti in [0:c.traces.size] do
          let (v, phi, d) := c.traces[ti]!
          s := s ++ s!"query {ti} {v} {hx phi} {hx d}\n"
          match traceFromVertex t0 v phi d with
          | .ok r => s := s ++ printTrace s!"tr{ti}" t0 r
          | .error e => s := s ++ s!"tr{ti} error {e}\n"
        return s ++ "end\n"
      | "edgetrace" =>
        let t1 ← applyFlips t0 c.preflips
        let t2 := if c.predelaunay then (delaunayFlip t1 c.crit c.tol c.limit).1 else t1
        let mut s := header c t0 topo "ok" ++ printState "state" t2
        for h in [0:t2.nh] do
          let u := t2.tail[h]!
          let w := t2.head h
          let phi := t2.sign[h]!
          let l := t2.len[h]!
          match traceFromVertex t0 u phi l with
          | .error e => s := s ++ s!"te {u} {w} {hx phi} {hx l} error {e}\n"
          | .ok r =>
            let miss := if r.stop == "vertex" && r.vertex == w then 0.0 else missDistance t0 r w
            let (b0, b1, b2) := r.bary
            s := s ++ s!"te {u} {w} {hx phi} {hx l} {r.stop} {r.face} {hx b0} {hx b1} {hx b2} {hx miss} {r.crossings.size}\n"
            for cr in r.crossings do
              let (p, q, tt) := crossingKey t0 cr
              s := s ++ s!"tx {u} {w} {hx phi} {p} {q} {hx tt}\n"
        return s ++ "end\n"
      | _ =>
        let faces := if c.circAll then (List.range t0.nf).toArray else c.circ
        let mut s := header c t0 topo "ok"
        for f in faces do
          if f ≥ t0.nf then throw s!"face {f} does not exist"
          let hs := t0.faceHalfedges f
          let l01 := t0.len[hs[0]!]!
          let l12 := t0.len[hs[1]!]!
          let l20 := t0.len[hs[2]!]!
          let (w0, w1, w2) := circumcenterBary l12 l20 l01
          let (r0, r1, r2) := cornerDistances l01 l12 l20 w0 w1 w2
          let inside := w0 > 0 && w1 > 0 && w2 > 0
          s := s ++ s!"cc {f} {hx w0} {hx w1} {hx w2} {hx r0} {hx r1} {hx r2} {if inside then 1 else 0}\n"
          let g := facePoint t0 f (1.0 / 3, 1.0 / 3, 1.0 / 3)
          let p := facePoint t0 f (w0, w1, w2)
          let dx := p.1 - g.1
          let dy := p.2 - g.2
          let dist := Float.sqrt (dx * dx + dy * dy)
          if dist == 0 then
            s := s ++ s!"cr{f} end {f} {hx (1.0 / 3)} {hx (1.0 / 3)} {hx (1.0 / 3)} {hx 0.0} -1 0 0\n"
          else
            match traceFromFace t0 f g (dx / dist, dy / dist) dist with
            | .ok r => s := s ++ printTrace s!"cr{f}" t0 r
            | .error e => s := s ++ s!"cr{f} error {e}\n"
        return s ++ "end\n"
    match res with
    | .ok s => s
    | .error e => printError c.name c.task e

/-- Read cases until the end of `input`, write one block per case, flush after each. -/
partial def runStream (input output : IO.FS.Stream) : IO Nat := do
  let reader ← LineReader.new input
  let rec loop (count : Nat) : IO Nat := do
    match ← readCaseLines reader with
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

end Ddg.Intrinsic
