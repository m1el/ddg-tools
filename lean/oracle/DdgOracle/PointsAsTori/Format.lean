import DdgOracle.Cuts.Format
import DdgOracle.PointsAsTori.Eval

/-!
# Text format of the `pat` command

Same conventions as the other papers (README, "Points as Tori"): line oriented, floats as 16 hex
digits, 0-based indices, one `case NAME ... end` block per case.
-/

namespace Ddg.PointsAsTori

/-- One case. -/
structure PatCase where
  name : String
  /-- `features`, `fit`, `sdf`, `blend` or `baselines` -/
  task : String
  pts : Array V3 := #[]
  nrm : Array V3 := #[]
  areas : Array Float := #[]
  center : Nat := 0
  /-- `nbhd K`: fitting neighborhood `[center] ++ (K-1 nearest others)` -/
  nbhdK : Option Nat := none
  /-- `nbr j0 j1 ...`: explicit neighborhood -/
  nbr : Option (Array Nat) := none
  y : Option (Array Float) := none
  coefs : Array (Nat × Array Float) := #[]
  fitOpts : FitOpts := {}
  tori : Array Torus := #[]
  queries : Array V3 := #[]
  opts : EvalOpts := {}
  deriving Inhabited

private def hx (x : Float) : String := floatToHex x
private def join (l : Array String) : String := " ".intercalate l.toList

private def natTok (ctx tok : String) : Except String Nat :=
  match tok.toNat? with
  | some n => pure n
  | none => throw s!"Syntax: {ctx}: expected a decimal integer, got '{tok}'"

private def floatTok (ctx tok : String) : Except String Float :=
  match hexToFloat? tok with
  | some x => pure x
  | none => throw s!"Syntax: {ctx}: expected 16 hex digits (IEEE-754 bit pattern), got '{tok}'"

private def floats (ctx : String) (l : Array String) (from_ cnt : Nat) : Except String (Array Float) := do
  if l.size != from_ + cnt then throw s!"Syntax: {ctx}: expected {cnt} values"
  let mut out := #[]
  for i in [from_:from_ + cnt] do out := out.push (← floatTok ctx l[i]!)
  return out

private def v3Of (a : Array Float) (i : Nat) : V3 := ⟨a[i]!, a[i + 1]!, a[i + 2]!⟩

/-- Largest cloud / query set accepted. -/
def maxPoints : Nat := 400
def maxQueries : Nat := 5000

/-- Parse the body lines of a case. -/
def parseBody (name : String) (lines : Array (Array String)) : Except String PatCase := do
  let mut c : PatCase := { name, task := "" }
  let mut declared : Option Nat := none
  for l in lines do
    match l[0]! with
    | "task" =>
      if l.size != 2 then throw "Syntax: task: expected one name"
      c := { c with task := l[1]! }
    | "cloud" =>
      if l.size != 2 then throw "Syntax: cloud: expected 'cloud N'"
      declared := some (← natTok "cloud" l[1]!)
    | "p" =>
      if l.size != 7 && l.size != 8 then throw "Syntax: p: expected 'p X Y Z NX NY NZ [AREA]'"
      let a ← floats "p" (l.extract 0 7) 1 6
      let ar ← if l.size == 8 then floatTok "p" l[7]! else pure 1.0
      c := { c with pts := c.pts.push (v3Of a 0), nrm := c.nrm.push (v3Of a 3), areas := c.areas.push ar }
    | "center" =>
      if l.size != 2 then throw "Syntax: center: expected 'center I'"
      c := { c with center := ← natTok "center" l[1]! }
    | "nbhd" =>
      if l.size != 2 then throw "Syntax: nbhd: expected 'nbhd K'"
      c := { c with nbhdK := some (← natTok "nbhd" l[1]!) }
    | "nbr" =>
      let mut a := #[]
      for i in [1:l.size] do a := a.push (← natTok "nbr" l[i]!)
      c := { c with nbr := some a }
    | "y" => c := { c with y := some (← floats "y" l 1 6) }
    | "coef" =>
      if l.size != 8 then throw "Syntax: coef: expected 'coef I a00 a01 a10 a11 a02 a20'"
      c := { c with coefs := c.coefs.push (← natTok "coef" l[1]!, ← floats "coef" l 2 6) }
    | "fitopts" =>
      let a ← floats "fitopts" l 1 2
      c := { c with fitOpts := { eps := a[0]!, floor := a[1]! } }
    | "torus" =>
      let a ← floats "torus" l 1 8
      c := { c with tori := c.tori.push { c := v3Of a 0, u := v3Of a 3, R := a[6]!, r := a[7]! } }
    | "x" =>
      let a ← floats "x" l 1 3
      c := { c with queries := c.queries.push (v3Of a 0) }
    | "neighborhood" =>
      match l.toList with
      | ["neighborhood", "knn", k] => c := { c with opts := { c.opts with nbhd := .knn (← natTok "knn" k) } }
      | ["neighborhood", "all"] => c := { c with opts := { c.opts with nbhd := .all } }
      | ["neighborhood", "radius", "auto", fb] =>
        c := { c with opts := { c.opts with nbhd := .radius none (← natTok "radius" fb) } }
      | ["neighborhood", "radius", r, fb] =>
        c := { c with opts := { c.opts with nbhd := .radius (some (← floatTok "radius" r)) (← natTok "radius" fb) } }
      | _ => throw "Syntax: neighborhood: expected 'knn K', 'all', 'radius R FALLBACK' or 'radius auto FALLBACK'"
    | "lambda" =>
      match l.toList with
      | ["lambda", "perquery", cc] => c := { c with opts := { c.opts with lam := .perQuery (← floatTok "lambda" cc) } }
      | ["lambda", "constant", v] => c := { c with opts := { c.opts with lam := .constant (← floatTok "lambda" v) } }
      | ["lambda", "paper", cc, f, km] =>
        c := { c with opts := { c.opts with lam := .paper (← floatTok "lambda" cc) (← floatTok "lambda" f) (← natTok "lambda" km) } }
      | _ => throw "Syntax: lambda: expected 'perquery C', 'constant L' or 'paper C FACTOR KMEAN'"
    | "areas" =>
      if l.size != 2 then throw "Syntax: areas: expected 'areas 0|1'"
      c := { c with opts := { c.opts with areas := l[1]! == "1" } }
    | "iso" => c := { c with opts := { c.opts with iso := (← floats "iso" l 1 1)[0]! } }
    | t => throw s!"Syntax: unknown line '{t}'"
  if let some d := declared then
    if d != c.pts.size then throw s!"Syntax: cloud declares {d} points, {c.pts.size} 'p' lines"
  if c.task == "" then throw "Syntax: missing 'task' line"
  if c.pts.size > maxPoints then throw s!"TooLarge: {c.pts.size} points, at most {maxPoints}"
  if c.queries.size > maxQueries then throw s!"TooLarge: {c.queries.size} queries, at most {maxQueries}"
  return c

/-- Print the input of a case (canonical form of the format). -/
def printCase (c : PatCase) : String := Id.run do
  let mut s := s!"case {c.name}\ntask {c.task}\n"
  if !c.pts.isEmpty then
    s := s ++ s!"cloud {c.pts.size}\n"
    for i in [0:c.pts.size] do
      let p := c.pts[i]!
      let n := c.nrm[i]!
      s := s ++ join #["p", hx p.x, hx p.y, hx p.z, hx n.x, hx n.y, hx n.z, hx c.areas[i]!] ++ "\n"
  if c.task == "features" then
    s := s ++ s!"center {c.center}\n"
    if let some k := c.nbhdK then s := s ++ s!"nbhd {k}\n"
    if let some a := c.nbr then s := s ++ join (#["nbr"] ++ a.map toString) ++ "\n"
    if let some y := c.y then s := s ++ join (#["y"] ++ y.map hx) ++ "\n"
  if c.task == "fit" || (c.task == "blend" && !c.coefs.isEmpty) then
    if c.task == "fit" then s := s ++ join #["fitopts", hx c.fitOpts.eps, hx c.fitOpts.floor] ++ "\n"
    for (i, a) in c.coefs do s := s ++ join (#["coef", toString i] ++ a.map hx) ++ "\n"
  for t in c.tori do
    s := s ++ join #["torus", hx t.c.x, hx t.c.y, hx t.c.z, hx t.u.x, hx t.u.y, hx t.u.z, hx t.R, hx t.r] ++ "\n"
  if c.task == "blend" || c.task == "baselines" then
    s := s ++ (match c.opts.nbhd with
      | .knn k => s!"neighborhood knn {k}\n"
      | .all => "neighborhood all\n"
      | .radius none fb => s!"neighborhood radius auto {fb}\n"
      | .radius (some r) fb => s!"neighborhood radius {hx r} {fb}\n")
    s := s ++ (match c.opts.lam with
      | .perQuery cc => s!"lambda perquery {hx cc}\n"
      | .constant l => s!"lambda constant {hx l}\n"
      | .paper cc f km => s!"lambda paper {hx cc} {hx f} {km}\n")
    s := s ++ s!"areas {if c.opts.areas then 1 else 0}\niso {hx c.opts.iso}\n"
  for q in c.queries do s := s ++ join #["x", hx q.x, hx q.y, hx q.z] ++ "\n"
  return s ++ "end\n"

def printError (name task msg : String) : String :=
  let kind := (msg.splitOn ":").head!
  let rest := (msg.drop (kind.length + 1)).trimAscii.toString
  s!"case {name}\ntask {task}\nstatus error {kind} {rest}\nend\n"

private def row (tag : String) (ints : Array Nat) (fs : Array Float) : String :=
  join (#[tag] ++ ints.map toString ++ fs.map hx) ++ "\n"

private def rowT (tag : String) (ints : Array Nat) (fs : Array Float) (tail : Array Nat) : String :=
  join (#[tag] ++ ints.map toString ++ fs.map hx ++ tail.map toString) ++ "\n"

private def v3a (v : V3) : Array Float := #[v.x, v.y, v.z]
private def b2s (b : Bool) : Nat := if b then 1 else 0

private def checkNormals (c : PatCase) : Except String Unit := do
  for i in [0:c.nrm.size] do
    let l := c.nrm[i]!.norm
    if !((l - 1.0).abs ≤ 1.0e-6) then throw s!"NonUnitNormal: normal {i} has norm {l.toString}"
  for p in c.pts do
    if !p.isFinite then throw "NonFinite: a point is not finite"

/-- The fitting neighborhood `[i] ++ (k-1 nearest j ≠ i by (ρ², j))`. -/
def fitNeighborhood (pts : Array V3) (i k : Nat) : Array Nat :=
  let s := sortedDist2 pts pts[i]!
  let others := (s.filter fun u => u.2 != i).map (·.2)
  #[i] ++ others.extract 0 (if k - 1 < others.size then k - 1 else others.size)

private def runFeatures (c : PatCase) : Except String String := do
  checkNormals c
  if c.pts.isEmpty then throw "Syntax: features needs a cloud"
  if c.center ≥ c.pts.size then throw s!"IndexOutOfRange: center {c.center}, {c.pts.size} points"
  let idx ← match c.nbr, c.nbhdK with
    | some a, _ => pure a
    | none, some k =>
      if k == 0 then throw "Syntax: nbhd 0" else pure (fitNeighborhood c.pts c.center k)
    | none, none => throw "Syntax: features needs 'nbhd K' or 'nbr ...'"
  if idx.isEmpty then throw "Syntax: empty neighborhood"
  for j in idx do
    if j ≥ c.pts.size then throw s!"IndexOutOfRange: neighbor {j}, {c.pts.size} points"
  if idx[0]! != c.center then throw "Syntax: the first neighbor must be the center"
  let f := features c.pts c.nrm idx
  if !(f.sigma > 0.0) then throw "DegenerateNeighborhood: median distance sigma is 0"
  let mut s := s!"status ok\ncenter {c.center} {idx.size}\n"
  s := s ++ row "sigma" #[] #[f.sigma]
  s := s ++ row "frame" #[] (v3a f.s ++ v3a f.t ++ v3a f.c)
  s := s ++ row "nbr" idx #[]
  for m in [0:f.rows.size] do s := s ++ row "f" #[m] f.rows[m]!
  if let some y := c.y then s := s ++ row "coef" #[] (outputScaling f.sigma y)
  return s

private def fitAll (c : PatCase) : Except String (Array FitResult) := do
  checkNormals c
  let mut out := #[]
  for (i, a) in c.coefs do
    if i ≥ c.pts.size then throw s!"IndexOutOfRange: coef point {i}, {c.pts.size} points"
    let r := fitTorus c.fitOpts c.pts[i]! c.nrm[i]! a
    if !r.finite then throw s!"NonFiniteTorus: the torus at point {i} is not finite"
    out := out.push r
  return out

private def runFit (c : PatCase) : Except String String := do
  let rs ← fitAll c
  let mut s := s!"status ok\nfitopts {hx c.fitOpts.eps} {hx c.fitOpts.floor}\ntori {rs.size}\n"
  for m in [0:rs.size] do
    let i := c.coefs[m]!.1
    let r := rs[m]!
    let t := r.torus
    s := s ++ row "t" #[i] (v3a t.c ++ v3a t.u ++ #[t.R, t.r, r.sign])
    s := s ++ row "d" #[i] #[r.A, r.H, r.K, r.kPaperPlus, r.kPaperMinus, r.kPlus, r.kMinus, r.kMin, r.kMax]
    s := s ++ row "dir" #[i] (v3a r.vMin ++ v3a r.vMax ++ v3a r.nStar)
    s := s ++ row "flags" #[i, b2s r.clampDisc, b2s r.clampEig, b2s r.clampDenom, b2s r.clampCurv, b2s r.umbilic] #[]
  return s

private def checkTori (ts : Array Torus) : Except String Unit := do
  for i in [0:ts.size] do
    let t := ts[i]!
    if !(t.c.isFinite && t.R.isFinite && t.r.isFinite && t.u.isFinite) then
      throw s!"NonFinite: torus {i}"
    if !((t.u.norm - 1.0).abs ≤ 1.0e-6) then throw s!"NonUnitAxis: torus {i} has axis norm {t.u.norm.toString}"
    if t.r == 0.0 then throw s!"ZeroMinorRadius: torus {i}"

private def runSdf (c : PatCase) : Except String String := do
  checkTori c.tori
  let mut s := s!"status ok\ntori {c.tori.size} queries {c.queries.size}\n"
  for q in [0:c.queries.size] do
    for t in [0:c.tori.size] do
      let o := c.tori[t]!.sdf c.queries[q]!
      s := s ++ rowT "g" #[q, t] (#[o.g] ++ v3a o.grad ++ #[o.lap]) #[o.kind]
  return s

private def toriOfCase (c : PatCase) : Except String (Array Torus) := do
  if !c.tori.isEmpty then
    if c.tori.size != c.pts.size then throw s!"ValueCount: {c.tori.size} tori for {c.pts.size} points"
    checkTori c.tori
    return c.tori
  if c.coefs.size != c.pts.size then throw s!"ValueCount: {c.coefs.size} coef lines for {c.pts.size} points"
  for i in [0:c.coefs.size] do
    if c.coefs[i]!.1 != i then throw s!"ValueCount: coef line {i} is for point {c.coefs[i]!.1}"
  let rs ← fitAll c
  return rs.map (·.torus)

private def paperLine (c : PatCase) : String :=
  match c.opts.lam with
  | .paper cc f km => let (D, l, R) := paperData c.pts cc f km; row "paper" #[] #[D, l, R]
  | _ => row "paper" #[] #[0.0, 0.0, 0.0]

private def runBlend (c : PatCase) : Except String String := do
  checkNormals c
  if c.pts.isEmpty then throw "Syntax: blend needs a cloud"
  let tori ← toriOfCase c
  let mut s := s!"status ok\ncloud {c.pts.size} queries {c.queries.size}\n" ++ paperLine c
  for q in [0:c.queries.size] do
    let b ← blend c.pts tori c.areas c.opts c.queries[q]!
    s := s ++ rowT "b" #[q] (#[b.phi] ++ v3a b.grad ++ v3a b.gradExact ++ #[b.lap, b.nset.sigma, b.nset.lambda, b.nset.shift]) #[b.nset.idx.size, b.singG, b.singR]
    s := s ++ row "nb" (#[q] ++ b.nset.idx) #[]
  return s

private def runBaselines (c : PatCase) : Except String String := do
  checkNormals c
  if c.pts.isEmpty then throw "Syntax: baselines needs a cloud"
  let mut s := s!"status ok\ncloud {c.pts.size} queries {c.queries.size}\n" ++ paperLine c
  for q in [0:c.queries.size] do
    let b ← baselines c.pts c.nrm c.areas c.opts c.queries[q]!
    s := s ++ rowT "s" #[q] (#[b.shc] ++ v3a b.shcGrad ++ #[b.sspd] ++ v3a b.sspdGrad ++ #[b.nset.sigma, b.nset.lambda, b.nset.shift]) #[b.nset.idx.size, b.skipped]
  return s

/-- Run a case; the output block. -/
def runCase (c : PatCase) : String :=
  let r : Except String String := match c.task with
    | "features" => runFeatures c
    | "fit" => runFit c
    | "sdf" => runSdf c
    | "blend" => runBlend c
    | "baselines" => runBaselines c
    | t => throw s!"Unsupported: unknown task '{t}'"
  match r with
  | .ok body => s!"case {c.name}\ntask {c.task}\n{body}end\n"
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

end Ddg.PointsAsTori
