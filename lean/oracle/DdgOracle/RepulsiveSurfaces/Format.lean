import DdgOracle.Cuts.Format
import DdgOracle.RepulsiveSurfaces.Remesh

/-!
# Text format of the `rs` command

Same conventions as `rc` (README, "Repulsive surfaces"): line oriented, floats as 16 hex digits,
0-based indices, one `case NAME ... end` block per case.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

/-- One case. -/
structure RSCase where
  /-- identifier -/
  name : String
  /-- `energy`, `df`, `metric`, `constraints`, `saddle`, `precond`, `step` or `remesh-ops` -/
  task : String
  /-- exponents -/
  ex : Ex := ⟨6.0, 12.0⟩
  /-- the mesh -/
  mesh : Mesh
  /-- `df`: vertex functions -/
  u : Array (Array Float) := #[]
  /-- `metric`: the order `τ` of `L_τ` -/
  tau : Float := 0.75
  /-- `metric`, `saddle`, `precond`, `step`: weights of `A` -/
  mo : MetricOpts := {}
  /-- `metric`, `precond`: print every matrix (else only `A`, `Lfr`) -/
  allMatrices : Bool := true
  /-- `precond`: `ε_rel` -/
  epsRel : Float := 1.0e-10
  /-- constraints -/
  cons : Array Cons := #[]
  /-- `precond`: the vector `r` -/
  r : Array V3 := #[]
  /-- `step` options -/
  so : StepOpts := {}
  /-- `remesh-ops`: operation, edge -/
  op : Option (String × Nat × Nat) := none
  /-- `remesh-ops`: pinned vertices -/
  pinned : Array Nat := #[]
  /-- `remesh-ops`: refuse flips that are not Delaunay -/
  delaunay : Bool := false
  /-- `remesh-ops`: refuse foldovers (normal turned by more than `π − 2`) -/
  foldover : Bool := false
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
  if i + 2 ≥ l.size then throw s!"{ctx}: expected three coordinates"
  return ⟨← floatTok ctx l[i]!, ← floatTok ctx l[i+1]!, ← floatTok ctx l[i+2]!⟩

/-- Largest mesh the oracle accepts (vertices). -/
def maxVertices : Nat := 200

private def compTok (ctx tok : String) : Except String (Option Nat) :=
  if tok == "all" then pure none else do pure (some (← natTok ctx tok))

private def parseCons (l : Array String) : Except String Cons := do
  if l.size < 2 then throw "c: expected a constraint kind"
  match l[1]! with
  | "barycenter" =>
    if l.size == 3 then return .barycenter (← compTok "c barycenter" l[2]!) none
    if l.size == 6 then return .barycenter (← compTok "c barycenter" l[2]!) (some (← v3Toks "c barycenter" l 3))
    throw "c barycenter: expected 'c barycenter all|C [X Y Z]'"
  | "pins" =>
    if l.size < 3 then throw "c pins: expected 'c pins N V ... [X Y Z ...]'"
    let n ← natTok "c pins" l[2]!
    if l.size != 3 + n && l.size != 3 + n + 3 * n then throw "c pins: expected N vertices and optionally 3N coordinates"
    let mut vs : Array Nat := #[]
    for k in [0:n] do vs := vs.push (← natTok "c pins" l[3 + k]!)
    if l.size == 3 + n then return .pins vs none
    let mut tg : Array V3 := #[]
    for k in [0:n] do tg := tg.push (← v3Toks "c pins" l (3 + n + 3 * k))
    return .pins vs (some tg)
  | "area" =>
    if l.size == 3 then return .area (← compTok "c area" l[2]!) none
    if l.size == 4 then return .area (← compTok "c area" l[2]!) (some (← floatTok "c area" l[3]!))
    throw "c area: expected 'c area all|C [TARGET]'"
  | "volume" =>
    if l.size == 3 then return .volume (← compTok "c volume" l[2]!) none
    if l.size == 4 then return .volume (← compTok "c volume" l[2]!) (some (← floatTok "c volume" l[3]!))
    throw "c volume: expected 'c volume all|C [TARGET]'"
  | other => throw s!"c: unknown constraint '{other}'"

/-- Parse the lines of one case (between `case` and `end`). -/
def parseBody (name : String) (lines : Array (Array String)) : Except String RSCase := do
  let mut c : RSCase := { name, task := "", mesh := { pos := #[], faces := #[] } }
  let mut mesh : Option Mesh := none
  let mut i := 0
  while i < lines.size do
    let l := lines[i]!
    match l[0]! with
    | "task" =>
      if l.size != 2 || !["energy", "df", "metric", "constraints", "saddle", "precond", "step", "remesh-ops"].contains l[1]! then
        throw "task: expected 'task energy|df|metric|constraints|saddle|precond|step|remesh-ops'"
      c := { c with task := l[1]! }
    | "exponents" =>
      if l.size != 3 then throw "exponents: expected 'exponents ALPHA BETA'"
      c := { c with ex := ⟨← floatTok "exponents" l[1]!, ← floatTok "exponents" l[2]!⟩ }
    | "u" =>
      let mut vals : Array Float := #[]
      for k in [1:l.size] do vals := vals.push (← floatTok "u" l[k]!)
      c := { c with u := c.u.push vals }
    | "tau" =>
      if l.size != 2 then throw "tau: expected 'tau T'"
      c := { c with tau := ← floatTok "tau" l[1]! }
    | "weights" =>
      if l.size != 4 then throw "weights: expected 'weights W_TP C_LO W_W'"
      c := { c with mo := { wtp := ← floatTok "weights" l[1]!, clo := ← floatTok "weights" l[2]!, ww := ← floatTok "weights" l[3]! } }
    | "matrices" =>
      if l.size != 2 then throw "matrices: expected all|a"
      c := { c with allMatrices := l[1]! == "all" }
    | "eps-rel" =>
      if l.size != 2 then throw "eps-rel: expected 'eps-rel E'"
      c := { c with epsRel := ← floatTok "eps-rel" l[1]! }
    | "c" => c := { c with cons := c.cons.push (← parseCons l) }
    | "r" =>
      if l.size != 5 then throw "r: expected 'r V X Y Z'"
      let v ← natTok "r" l[1]!
      if v != c.r.size then throw "r: lines must come in vertex order"
      c := { c with r := c.r.push (← v3Toks "r" l 2) }
    | "armijo" =>
      if l.size != 2 then throw "armijo: expected 'armijo C'"
      c := { c with so := { c.so with armijo := ← floatTok "armijo" l[1]! } }
    | "start" =>
      if l.size == 3 && l[1]! == "abs" then c := { c with so := { c.so with start := .abs (← floatTok "start" l[2]!) } }
      else if l.size == 3 && l[1]! == "maxdisp" then c := { c with so := { c.so with start := .maxdisp (← floatTok "start" l[2]!) } }
      else if l.size == 2 && l[1]! == "reference" then c := { c with so := { c.so with start := .reference } }
      else throw "start: expected 'start abs T | maxdisp LENGTHS | reference'"
    | "maxhalv" =>
      if l.size != 2 then throw "maxhalv: expected a count"
      c := { c with so := { c.so with maxHalv := ← natTok "maxhalv" l[1]! } }
    | "newton" =>
      if l.size != 3 || (l[1]! != "fixed" && l[1]! != "tol") then throw "newton: expected 'newton fixed|tol N'"
      c := { c with so := { c.so with newton := ← natTok "newton" l[2]!, newtonFixed := l[1]! == "fixed" } }
    | "op" =>
      if l.size != 4 || !["split", "collapse", "flip"].contains l[1]! then throw "op: expected 'op split|collapse|flip A B'"
      c := { c with op := some (l[1]!, ← natTok "op" l[2]!, ← natTok "op" l[3]!) }
    | "pinned" =>
      let mut vs : Array Nat := #[]
      for k in [1:l.size] do vs := vs.push (← natTok "pinned" l[k]!)
      c := { c with pinned := vs }
    | "legality" =>
      for k in [1:l.size] do
        match l[k]! with
        | "delaunay" => c := { c with delaunay := true }
        | "foldover" => c := { c with foldover := true }
        | other => throw s!"legality: unknown check '{other}'"
    | "mesh" =>
      if l.size != 3 then throw "mesh: expected 'mesh NV NF'"
      let nv ← natTok "mesh" l[1]!
      let nf ← natTok "mesh" l[2]!
      if nv > maxVertices then throw s!"mesh has {nv} vertices; the oracle accepts at most {maxVertices}"
      if i + nv + nf ≥ lines.size then throw "mesh: missing lines"
      let mut pos : Array V3 := #[]
      let mut faces : Array Tri := #[]
      for k in [0:nv] do
        let p := lines[i + 1 + k]!
        if p[0]! != "p" || p.size != 4 then throw s!"malformed line '{joinTokens p}': expected 'p <x> <y> <z>'"
        pos := pos.push (← v3Toks "p" p 1)
      for k in [0:nf] do
        let f := lines[i + 1 + nv + k]!
        if f[0]! != "f" || f.size != 4 then throw s!"malformed line '{joinTokens f}': expected 'f <a> <b> <c>'"
        faces := faces.push ⟨← natTok "f" f[1]!, ← natTok "f" f[2]!, ← natTok "f" f[3]!⟩
      mesh := some { pos, faces }
      i := i + nv + nf
    | other => throw s!"unknown line '{other}'"
    i := i + 1
  if c.task == "" then throw "missing 'task' line"
  let some m := mesh | throw "missing 'mesh' block"
  return { c with mesh := m }

private def consLine : Cons → String
  | .barycenter comp tg =>
    let cs := match comp with | none => "all" | some k => toString k
    s!"c barycenter {cs}" ++ (match tg with | none => "" | some t => s!" {hx t.x} {hx t.y} {hx t.z}") ++ "\n"
  | .pins vs tg =>
    s!"c pins {vs.size} " ++ " ".intercalate (vs.map toString).toList
      ++ (match tg with | none => "" | some ts => ts.foldl (fun acc t => acc ++ s!" {hx t.x} {hx t.y} {hx t.z}") "") ++ "\n"
  | .area comp tg =>
    s!"c area {match comp with | none => "all" | some k => toString k}" ++ (match tg with | none => "" | some t => s!" {hx t}") ++ "\n"
  | .volume comp tg =>
    s!"c volume {match comp with | none => "all" | some k => toString k}" ++ (match tg with | none => "" | some t => s!" {hx t}") ++ "\n"

private def startLine : Start → String
  | .abs t => s!"start abs {hx t}\n"
  | .maxdisp l => s!"start maxdisp {hx l}\n"
  | .reference => "start reference\n"

/-- The case in the input format. -/
def printCase (c : RSCase) : String := Id.run do
  let mut s := s!"case {c.name}\ntask {c.task}\n"
  if c.task != "remesh-ops" && c.task != "constraints" && c.task != "df" then
    s := s ++ s!"exponents {hx c.ex.alpha} {hx c.ex.beta}\n"
  if c.task == "df" then
    for u in c.u do s := s ++ "u" ++ u.foldl (fun acc x => acc ++ " " ++ hx x) "" ++ "\n"
  if c.task == "metric" then
    s := s ++ s!"tau {hx c.tau}\nweights {hx c.mo.wtp} {hx c.mo.clo} {hx c.mo.ww}\n"
    if !c.allMatrices then s := s ++ "matrices a\n"
  if c.task == "saddle" || c.task == "step" then
    s := s ++ s!"weights {hx c.mo.wtp} {hx c.mo.clo} {hx c.mo.ww}\n"
  if c.task == "precond" then
    s := s ++ s!"eps-rel {hx c.epsRel}\n"
    if !c.allMatrices then s := s ++ "matrices a\n"
    for k in [0:c.r.size] do s := s ++ s!"r {k} {hx c.r[k]!.x} {hx c.r[k]!.y} {hx c.r[k]!.z}\n"
  if c.task == "step" then
    s := s ++ s!"armijo {hx c.so.armijo}\n" ++ startLine c.so.start ++ s!"maxhalv {c.so.maxHalv}\n"
    s := s ++ s!"newton {if c.so.newtonFixed then "fixed" else "tol"} {c.so.newton}\n"
  if c.task == "remesh-ops" then
    if let some (k, a, b) := c.op then s := s ++ s!"op {k} {a} {b}\n"
    if c.pinned.size > 0 then s := s ++ "pinned " ++ " ".intercalate (c.pinned.map toString).toList ++ "\n"
    if c.delaunay || c.foldover then
      s := s ++ "legality" ++ (if c.delaunay then " delaunay" else "") ++ (if c.foldover then " foldover" else "") ++ "\n"
  for k in c.cons do s := s ++ consLine k
  s := s ++ s!"mesh {c.mesh.nv} {c.mesh.nf}\n"
  for p in c.mesh.pos do s := s ++ s!"p {hx p.x} {hx p.y} {hx p.z}\n"
  for f in c.mesh.faces do s := s ++ s!"f {f.a} {f.b} {f.c}\n"
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

private def v3Lines (tag : String) (a : Array V3) : String :=
  a.foldl (fun acc p => acc) "" ++ (Id.run do
    let mut s := ""
    for v in [0:a.size] do s := s ++ s!"{tag} {v} {hx a[v]!.x} {hx a[v]!.y} {hx a[v]!.z}\n"
    return s)

private def flatV3Lines (tag : String) (x : FloatArray) : String := v3Lines tag (unflatten x)

private def head (c : RSCase) (t : Topo) : String :=
  s!"case {c.name}\ntask {c.task}\nstatus ok\nmesh {c.mesh.nv} {c.mesh.nf} {t.ncomp}\n"

private def maxAbsDiff (a b : FloatArray) : Float := Id.run do
  let mut m := 0.0
  for i in [0:a.size] do m := max m (a.get! i - b.get! i).abs
  return m

private def phiLine (rows : Array (String × Bool × CRow)) : String :=
  "phi" ++ rows.foldl (fun acc r => acc ++ " " ++ hx r.2.2.phi) "" ++ "\n"

private def runEnergy (c : RSCase) (t : Topo) : Except String String := do
  c.ex.check
  let (e, g) ← energyGrad t.faces c.mesh.pos c.ex
  return head c t ++ s!"energy {hx e}\n" ++ v3Lines "d" g ++ "end\n"

private def runDf (c : RSCase) (t : Topo) : Except String String := do
  let g := faceGeo t.faces c.mesh.pos
  let mut s := head c t
  for k in [0:c.u.size] do
    let u := c.u[k]!
    if u.size != t.nv then throw s!"ValueCount: u line {k} has {u.size} values for {t.nv} vertices"
    for fi in [0:t.faces.size] do
      let f := t.faces[fi]!
      let ub := (u[f.a]! + u[f.b]! + u[f.c]!) / 3.0
      let d := dfFace c.mesh.pos f g[fi]! u
      s := s ++ s!"df {k} {fi} {hx ub} {hx d.x} {hx d.y} {hx d.z}\n"
  return s ++ "end\n"

private def runMetric (c : RSCase) (t : Topo) : Except String String := do
  c.ex.check
  let p ← assemblePairs t c.mesh.pos c.tau c.ex.s
  let ops ← buildOps t c.mesh.pos c.ex c.mo
  let n := t.nv
  let mut s := head c t ++ s!"tau {hx c.tau}\ns {hx c.ex.s}\n"
  if c.allMatrices then
    s := s ++ matrix "L" p.l n n ++ matrix "B" ops.b n n ++ matrix "B0" ops.b0 n n ++ matrix "W" ops.w n n
      ++ matrix "WMW" ops.wmw n n
  s := s ++ matrix "A" ops.a n n
  s := s ++ "M" ++ ops.m.foldl (fun acc x => acc ++ " " ++ hx x) "" ++ "\n"
  return s ++ "end\n"

private def runConstraints (c : RSCase) (t : Topo) : Except String String := do
  let rc ← resolveCons t c.mesh.pos c.cons
  let rows := evalAll t c.mesh.pos rc
  let n3 := 3 * t.nv
  let mut s := head c t ++ s!"rows {rows.size} {n3}\n"
  s := s ++ "vcomp" ++ t.vertComp.foldl (fun acc k => acc ++ s!" {k}") "" ++ "\n"
  for i in [0:rows.size] do
    let (k, q, r) := rows[i]!
    s := s ++ s!"row {i} {k} {if q then "positional" else "nonlinear"} {hx r.phi}\n"
  for i in [0:rows.size] do
    s := s ++ row "J" i (denseRow n3 rows[i]!.2.2) 0 n3
  return s ++ "end\n"

private def runSaddle (c : RSCase) (t : Topo) : Except String String := do
  c.ex.check
  let rc ← resolveCons t c.mesh.pos c.cons
  checkTranslation t rc
  let (_, grad) ← energyGrad t.faces c.mesh.pos c.ex
  let ops ← buildOps t c.mesh.pos c.ex c.mo
  let rows := evalAll t c.mesh.pos rc
  let nst ← buildNested t.nv ops.a rows
  let gf := flatten grad
  let x := nst.direction gf
  let nrows := (rows.filter (fun r => !r.2.1)).map (·.2.2)
  let phiN := FloatArray.ofFn' nrows.size fun i => nrows[i]!.phi
  let h := nst.correction phiN
  let cnx := nst.applyCn x
  let mut cq := 0.0
  for r in [0:nst.ap.kq] do
    let mut sum := 0.0
    for j in [0:nst.ap.n3] do sum := sum + nst.ap.cq.get! (r * nst.ap.n3 + j) * x.get! j
    cq := max cq sum.abs
  let mut s := head c t ++ s!"rows {nst.ap.kq} {nst.kn}\npivot {hx nst.ap.pivot}\n"
  s := s ++ "phin" ++ (Array.range nst.kn).foldl (fun acc i => acc ++ " " ++ hx (phiN.get! i)) "" ++ "\n"
  s := s ++ v3Lines "g" grad ++ flatV3Lines "x" x
  s := s ++ matrix "S" nst.s nst.kn nst.kn
  s := s ++ flatV3Lines "h" h
  s := s ++ s!"gx {hx (dotF gf x)}\ncnx {hx (maxAbs cnx)}\ncqx {hx cq}\n"
  return s ++ "end\n"

private def runPrecond (c : RSCase) (t : Topo) : Except String String := do
  c.ex.check
  let rc ← resolveCons t c.mesh.pos c.cons
  for k in rc do
    if !k.positional then throw "Unsupported: precond takes positional constraints (barycenter, pins) only"
  if c.r.size != t.nv then throw s!"ValueCount: {c.r.size} r lines for {t.nv} vertices"
  let rows := evalAll t c.mesh.pos rc
  let w := stiffness t c.mesh.pos
  let m := vertexAreas t (faceGeo t.faces c.mesh.pos)
  let (eps, lmat) := buildPrecondMatrix t.nv w m c.epsRel
  let ap ← buildAp t.nv lmat (rows.map (·.2.2))
  let lfr := (← assemblePairs t c.mesh.pos (2.0 - c.ex.s) c.ex.s).l
  let r := flatten c.r
  let z := applyPrecond t.nv ap lfr r
  let n := t.nv
  let n3 := 3 * n
  let nn := n3 + ap.kq
  let mut s := head c t ++ s!"rows {ap.kq}\npivot {hx ap.pivot}\neps {hx eps}\n"
  s := s ++ matrix "Lfr" lfr n n
  if c.allMatrices then
    -- L′ rebuilt densely for printing
    let mut lp := zeros (nn * nn)
    for v in [0:n] do
      for u in [0:n] do
        for d in [0:3] do
          lp := lp.set! ((3 * v + d) * nn + 3 * u + d) (lmat.get! (v * n + u))
    for k in [0:ap.kq] do
      for j in [0:n3] do
        lp := lp.set! ((n3 + k) * nn + j) (ap.cq.get! (k * n3 + j))
        lp := lp.set! (j * nn + n3 + k) (ap.cq.get! (k * n3 + j))
    s := s ++ matrix "Lp" lp nn nn
  s := s ++ flatV3Lines "z" z
  return s ++ "end\n"

private def runStep (c : RSCase) (t : Topo) : Except String String := do
  c.ex.check
  let rc ← resolveCons t c.mesh.pos c.cons
  let r ← projectedStep t c.ex c.mo rc c.so c.mesh.pos
  let mut s := head c t ++ s!"pivot {hx r.pivot}\nenergy0 {hx r.f0}\nslope {hx r.slope}\nxnorm {hx r.xnorm}\nlbar {hx r.lbar}\ntau0 {hx r.tau0}\n"
  for i in [0:r.trials.size] do
    let tr := r.trials[i]!
    s := s ++ s!"trial {i} {hx tr.tau} {tr.outcome} {hx tr.energy} {hx tr.margin}\n"
  s := s ++ s!"result {if r.accepted then "accepted" else "stuck"} {hx r.tau} {r.halvings}\n"
  s := s ++ s!"energy1 {hx r.e1}\nnewton {r.newtonIters}\n"
  s := s ++ "phi-before" ++ r.phiBefore.foldl (fun acc x => acc ++ " " ++ hx x) "" ++ "\n"
  s := s ++ "phi" ++ r.phi.foldl (fun acc x => acc ++ " " ++ hx x) "" ++ "\n"
  s := s ++ s!"energy-final {hx r.eFinal}\n" ++ v3Lines "x" r.pos
  return s ++ "end\n"

private def runRemesh (c : RSCase) (t : Topo) : Except String String := do
  let some (kind, a, b) := c.op | throw "MissingOp: no 'op' line"
  for v in c.pinned do
    if v ≥ t.nv then throw s!"ConstraintOutOfRange: pinned vertex {v}"
  let out ← match kind with
    | "split" => splitEdge c.mesh t a b
    | "collapse" => collapseEdge c.mesh t a b c.pinned c.foldover
    | _ => flipEdge c.mesh t a b c.delaunay c.foldover
  let mut s := s!"case {c.name}\ntask {c.task}\nstatus ok\nmesh {out.pos.size} {out.faces.size}\n"
  for v in [0:out.pos.size] do s := s ++ s!"p {hx out.pos[v]!.x} {hx out.pos[v]!.y} {hx out.pos[v]!.z}\n"
  for f in out.faces do s := s ++ s!"f {f.a} {f.b} {f.c}\n"
  for i in [0:out.parent.size] do s := s ++ s!"parent {i} {out.parent[i]!}\n"
  for v in [0:out.vmap.size] do s := s ++ s!"vmap {v} " ++ " ".intercalate (out.vmap[v]!.map toString).toList ++ "\n"
  return s ++ "end\n"

/-- Run one case and print its output block. -/
def runCase (c : RSCase) : String :=
  let res : Except String String := do
    let t ← c.mesh.validate
    match c.task with
    | "energy" => runEnergy c t
    | "df" => runDf c t
    | "metric" => runMetric c t
    | "constraints" => runConstraints c t
    | "saddle" => runSaddle c t
    | "precond" => runPrecond c t
    | "step" => runStep c t
    | _ => runRemesh c t
  match res with
  | .ok s => s
  | .error e => printError c.name c.task e

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

end Ddg.RepulsiveSurfaces
