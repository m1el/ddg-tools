import DdgOracle.Format
import DdgOracle.Cuts.Step

/-!
# Text format of the `cuts` task

Specified in `README.md` ("Cut flow"). Same conventions as the first task: version 1, line
oriented, floats as 16 hex digits, 0-based indices.
-/

namespace Ddg.Cuts

private def joinTokens (l : Array String) : String := " ".intercalate l.toList

private def natTok (ctx tok : String) : Except String Nat :=
  match tok.toNat? with
  | some n => pure n
  | none => throw s!"{ctx}: expected a decimal integer, got '{tok}'"

private def floatTok (ctx tok : String) : Except String Float :=
  match hexToFloat? tok with
  | some x => pure x
  | none => throw s!"{ctx}: expected 16 hex digits (IEEE-754 bit pattern), got '{tok}'"

private def floatsFrom (ctx : String) (l : Array String) (start count : Nat) :
    Except String FloatArray := do
  if l.size != start + count then
    throw s!"{ctx}: expected {count} values, got {l.size - start}"
  let mut out := zeros 0
  for i in [start:start + count] do
    out := out.push (← floatTok ctx l[i]!)
  return out

/-- Parse the lines of one case (between `case` and `end`). -/
def parseBody (name : String) (lines : Array (Array String)) : Except String CutCase := do
  let mut i := 0
  let line := fun (k : Nat) (what : String) =>
    if h : k < lines.size then pure (f := Except String) lines[k]
    else throw s!"missing line: expected {what}"
  let l ← line i "'task cuts' or 'mesh <nv> <nf> <ne>'"
  if l[0]! == "task" then
    if l.size != 2 || l[1]! != "cuts" then
      throw s!"malformed line '{joinTokens l}': expected 'task cuts'"
    i := i + 1
  let l ← line i "'mesh <nv> <nf> <ne>'"
  if l[0]! != "mesh" || l.size != 4 then
    throw s!"malformed line '{joinTokens l}': expected 'mesh <nv> <nf> <ne>'"
  let nv ← natTok "mesh" l[1]!
  let nf ← natTok "mesh" l[2]!
  let ne ← natTok "mesh" l[3]!
  if nv > maxCutVertices then throw s!"mesh has {nv} vertices; the dense oracle accepts at most {maxCutVertices}"
  if nf > 4 * maxCutVertices || ne > 6 * maxCutVertices then throw "mesh too large"
  i := i + 1
  let mut faces : Array Face := Array.mkEmpty nf
  for f in [0:nf] do
    let l ← line i s!"'f <a> <b> <c>' for face {f}"
    if l[0]! != "f" || l.size != 4 then
      throw s!"malformed line '{joinTokens l}': expected 'f <a> <b> <c>' for face {f}"
    let ctx := s!"face {f}"
    faces := faces.push ⟨← natTok ctx l[1]!, ← natTok ctx l[2]!, ← natTok ctx l[3]!⟩
    i := i + 1
  let mut edges : Array (Nat × Nat × Float) := Array.mkEmpty ne
  for e in [0:ne] do
    let l ← line i s!"'e <a> <b> <length>' for edge {e}"
    if l[0]! != "e" || l.size != 4 then
      throw s!"malformed line '{joinTokens l}': expected 'e <a> <b> <length>' for edge {e}"
    let ctx := s!"edge {e}"
    edges := edges.push (← natTok ctx l[1]!, ← natTok ctx l[2]!, ← floatTok ctx l[3]!)
    i := i + 1
  -- keyed parameter lines, in any order
  let mut kv : Std.HashMap String (Array String) := {}
  while i < lines.size do
    let l := lines[i]!
    if kv.contains l[0]! then throw s!"line '{l[0]!}' occurs twice"
    kv := kv.insert l[0]! l
    i := i + 1
  let allowed := ["phi", "energy", "curvature", "alpha", "eta", "penalty", "w", "beta", "step",
    "delta", "flux", "rho0", "passes", "classes", "format", "stages", "featurelap", "statistics"]
  for (k, _) in kv.toList do
    if !allowed.contains k then throw s!"unknown line '{k}'"
  let format ← match kv.get? "format" with
    | none => pure 1
    | some l =>
      if l.size != 2 || (l[1]! != "2" && l[1]! != "3") then throw "format: expected 'format 2' or 'format 3'"
      pure (if l[1]! == "3" then 3 else 2)
  if format == 1 then
    for k in ["stages", "featurelap"] do
      if kv.contains k then throw s!"line '{k}' needs 'format 2'"
  if format < 3 && kv.contains "statistics" then throw "line 'statistics' needs 'format 3'"
  let need := fun (k : String) =>
    match kv.get? k with
    | some l => pure (f := Except String) l
    | none => throw s!"missing line '{k}'"
  let phi ← do let l ← need "phi"; floatsFrom "phi" l 1 nv
  let energy ← do
    let l ← need "energy"
    match l[1]? with
    | some "dirichlet" =>
      if l.size != 2 then throw "energy dirichlet takes no arguments"
      pure Energy.dirichlet
    | some "hencky" =>
      if l.size != 2 then throw "energy hencky takes no arguments"
      pure Energy.hencky
    | some "scaled_a" => pure (Energy.scaledA (← floatsFrom "energy scaled_a" l 2 nv))
    | some "scaled_curv" =>
      if l.size != 4 then throw "energy scaled_curv: expected 'energy scaled_curv <r_eps> <c>'"
      pure (Energy.scaledCurv (← floatTok "r_eps" l[2]!) (← floatTok "c" l[3]!))
    | _ => throw s!"malformed line '{joinTokens l}': energy is dirichlet | hencky | scaled_a | scaled_curv"
  let curv ← do
    let l ← need "curvature"
    match l[1]? with
    | some "gaussian" =>
      if l.size != 2 then throw "curvature gaussian takes no arguments"
      pure CurvSource.gaussian
    | some "uniform" =>
      if l.size != 3 then throw "curvature uniform: expected one value"
      pure (CurvSource.uniform (← floatTok "curvature" l[2]!))
    | some "per_face" => pure (CurvSource.perFace (← floatsFrom "curvature per_face" l 2 nf))
    | _ => throw s!"malformed line '{joinTokens l}': curvature is gaussian | uniform | per_face"
  let alphaL ← do
    let l ← need "alpha"
    if l.size != 2 then throw "alpha: expected one value"
    floatTok "alpha" l[1]!
  let eta ← match kv.get? "eta" with
    | some l => pure (some (← floatsFrom "eta" l 1 nv))
    | none => pure none
  let penalty ← do
    let l ← need "penalty"
    match l[1]? with
    | some "gradient" => pure Penalty.gradient
    | some "paper" => pure Penalty.paperDiagonal
    | _ => throw s!"malformed line '{joinTokens l}': penalty is gradient | paper"
  let two := fun (k : String) => do
    let l ← need k
    if l.size != 3 then throw s!"{k}: expected two values"
    pure (← floatTok k l[1]!, ← floatTok k l[2]!)
  let w ← two "w"
  let beta ← two "beta"
  let step ← do
    let l ← need "step"
    match l[1]? with
    | some "fixed" =>
      if l.size != 3 then throw "step fixed: expected one value"
      pure (StepRule.fixed (← floatTok "step" l[2]!))
    | some "capped" =>
      if l.size != 4 then throw "step capped: expected two values"
      pure (StepRule.capped (← floatTok "step" l[2]!) (← floatTok "step" l[3]!))
    | some "adaptive" =>
      if format < 2 then throw "step adaptive needs 'format 2'"
      if l.size != 5 then throw "step adaptive: expected 'step adaptive <time_scale> <max_motion> <tau_max>'"
      pure (StepRule.adaptive (← floatTok "step" l[2]!) (← floatTok "step" l[3]!) (← floatTok "step" l[4]!))
    | _ => throw s!"malformed line '{joinTokens l}': step is fixed <tau> | capped <tau_max> <max_motion> | adaptive <time_scale> <max_motion> <tau_max>"
  let delta ← do
    let l ← need "delta"
    if l.size != 2 then throw "delta: expected one value"
    floatTok "delta" l[1]!
  let flux ← do
    let l ← need "flux"
    match l[1]? with
    | some "nodal" =>
      if l.size != 2 then throw "flux nodal takes no arguments"
      pure FluxRule.nodal
    | some "smoothed" =>
      if l.size != 3 then throw "flux smoothed: expected one value"
      pure (FluxRule.smoothed (← floatTok "flux" l[2]!))
    | _ => throw s!"malformed line '{joinTokens l}': flux is nodal | smoothed <c>"
  let rho0 ← do
    let l ← need "rho0"
    if l.size != 2 then throw "rho0: expected one value"
    floatTok "rho0" l[1]!
  let passes ← do
    let l ← need "passes"
    if l.size != 2 then throw "passes: expected one integer"
    natTok "passes" l[1]!
  let classes ← match kv.get? "classes" with
    | some l =>
      if l.size != nv + 1 then throw s!"classes: expected {nv} values"
      let mut out : Array Nat := Array.mkEmpty nv
      for j in [1:nv + 1] do out := out.push (← natTok "classes" l[j]!)
      pure (some out)
    | none => pure none
  let (lengthStep, extension, speed) ← do
    if format == 1 then pure (LengthStep.onMesh, Extension.harmonic, SpeedScaling.unit)
    else
      let l ← need "stages"
      if l.size != 4 then throw "stages: expected 'stages <onmesh|oncut> <harmonic|local> <unit|mean_gradient>'"
      let ls ← match l[1]! with
        | "onmesh" => pure LengthStep.onMesh
        | "oncut" => pure LengthStep.onCut
        | o => throw s!"stages: unknown length step '{o}'"
      let ex ← match l[2]! with
        | "harmonic" => pure Extension.harmonic
        | "local" => pure Extension.local
        | o => throw s!"stages: unknown extension '{o}'"
      let sp ← match l[3]! with
        | "unit" => pure SpeedScaling.unit
        | "mean_gradient" => pure SpeedScaling.meanGradient
        | o => throw s!"stages: unknown speed scaling '{o}'"
      pure (ls, ex, sp)
  let featPos ← match kv.get? "featurelap" with
    | none => pure false
    | some l =>
      match l[1]? with
      | some "cotan" => pure false
      | some "positive" => pure true
      | _ => throw "featurelap: expected cotan | positive"
  let stats ← match kv.get? "statistics" with
    | none => pure StepStats.extremes
    | some l =>
      match l[1]? with
      | some "extremes" => pure StepStats.extremes
      | some "robust" => pure StepStats.robust
      | _ => throw "statistics: expected extremes | robust"
  return { name, nv, faces, edges, phi, energy, curv, alphaL, eta, penalty, w, beta, step
           delta, flux, rho0, passes, classes, format, lengthStep, extension, speed, featPos, stats }

/-- Read the next case as raw lines: `(name, lines)`, or an error message. -/
partial def readCaseLines (r : LineReader) :
    IO (Option (String × Except String (Array (Array String)))) := do
  match ← r.next with
  | none => return none
  | some first =>
    if first[0]! != "case" || first.size != 2 then
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
    | none => return some (name, .ok lines)

/-! ## Printing -/

private def hx (x : Float) : String := floatToHex x

private def b2n (b : Bool) : String := if b then "1" else "0"

/-- A case in the input format. -/
def printCase (c : CutCase) : String := Id.run do
  let mut s := s!"case {c.name}\ntask cuts\nmesh {c.nv} {c.faces.size} {c.edges.size}\n"
  for f in c.faces do s := s ++ s!"f {f.a} {f.b} {f.c}\n"
  for (a, b, l) in c.edges do s := s ++ s!"e {a} {b} {hx l}\n"
  let arr := fun (x : FloatArray) => " ".intercalate (x.toList.map hx)
  s := s ++ s!"phi {arr c.phi}\n"
  match c.energy with
  | .dirichlet => s := s ++ "energy dirichlet\n"
  | .hencky => s := s ++ "energy hencky\n"
  | .scaledA a => s := s ++ s!"energy scaled_a {arr a}\n"
  | .scaledCurv r k => s := s ++ s!"energy scaled_curv {hx r} {hx k}\n"
  match c.curv with
  | .gaussian => s := s ++ "curvature gaussian\n"
  | .uniform k => s := s ++ s!"curvature uniform {hx k}\n"
  | .perFace k => s := s ++ s!"curvature per_face {arr k}\n"
  s := s ++ s!"alpha {hx c.alphaL}\n"
  if let some e := c.eta then s := s ++ s!"eta {arr e}\n"
  s := s ++ (match c.penalty with | .gradient => "penalty gradient\n" | .paperDiagonal => "penalty paper\n")
  s := s ++ s!"w {hx c.w.1} {hx c.w.2}\nbeta {hx c.beta.1} {hx c.beta.2}\n"
  match c.step with
  | .fixed t => s := s ++ s!"step fixed {hx t}\n"
  | .capped t m => s := s ++ s!"step capped {hx t} {hx m}\n"
  | .adaptive ts m tm => s := s ++ s!"step adaptive {hx ts} {hx m} {hx tm}\n"
  s := s ++ (match c.flux with | .nodal => "flux nodal\n" | .smoothed k => s!"flux smoothed {hx k}\n")
  s := s ++ s!"delta {hx c.delta}\nrho0 {hx c.rho0}\npasses {c.passes}\n"
  if c.format ≥ 2 then
    s := s ++ s!"format {c.format}\n"
    if c.format == 3 then
      s := s ++ s!"statistics {match c.stats with | .extremes => "extremes" | .robust => "robust"}\n"
    s := s ++ s!"stages {match c.lengthStep with | .onMesh => "onmesh" | .onCut => "oncut"} {match c.extension with | .harmonic => "harmonic" | .local => "local"} {match c.speed with | .unit => "unit" | .meanGradient => "mean_gradient"}\n"
    s := s ++ s!"featurelap {if c.featPos then "positive" else "cotan"}\n"
  if let some cl := c.classes then
    s := s ++ s!"classes {" ".intercalate (cl.toList.map toString)}\n"
  return s ++ "end\n"

/-- Output block of a case that could not be processed. -/
def printError (name msg : String) : String :=
  s!"case {name}\ntask cuts\nstatus error {msg}\nend\n"

private def energyLine (tag : String) (ev : Eval) : String :=
  s!"{tag} {hx ev.total} {hx ev.J0} {hx ev.J1} {hx ev.A0} {hx ev.A1} {hx ev.lengthTerm}\n"

/-- The lines that describe the evaluation at the input. -/
private def printEval (c : CutCase) (b : Base) (phi : FloatArray) (ev : Eval) : String := Id.run do
  let nv := b.nv
  let cm := ev.cm
  let mut movedCount := 0
  for cr in cm.crossings do if cr.moved then movedCount := movedCount + 1
  let mut sliverVerts := 0
  for x in ev.sliver do if x then sliverVerts := sliverVerts + 1
  let mut sumOmega := fzero
  for o in b.omega do sumOmega := sumOmega + o
  let mut s := s!"scale {hx b.scale} {hx b.totalArea}\n"
  s := s ++ s!"mesh {nv} {c.faces.size} {b.topo.numEdges} {b.topo.numBoundary} {b.topo.chi}\n"
  s := s ++ s!"sum_omega {hx sumOmega}\n"
  s := s ++ s!"counts {cm.nc} {cm.segs.size} {cm.tris.size} {movedCount} {sliverVerts} {ev.nSliverTris} {cm.nearTies}\n"
  for fi in [0:c.faces.size] do
    s := s ++ s!"f {fi} {hx b.omega[fi]!} {hx b.kGauss[fi]!} {hx b.kUsed[fi]!} {hx b.faceArea[fi]!}\n"
  for i in [0:nv] do
    s := s ++ s!"v {i} {cm.region[i]!} {hx (b.theta.get! i)} {hx (b.vertArea.get! i)} {hx (b.feat.get! i)} {hx (b.aVert.get! i)} {hx (ev.wv.get! i)} {hx (phi.get! i)} {hx (ev.sigma.get! i)}\n"
  for k in [0:cm.nc] do
    let cr := cm.crossings[k]!
    s := s ++ s!"x {k} {cr.a} {cr.b} {hx cr.tRaw} {hx cr.t} {b2n cr.moved} {hx cr.len}\n"
  for ti in [0:cm.tris.size] do
    let t := cm.tris[ti]!
    s := s ++ s!"t {ti} {t.parent} {t.region} {t.n0} {t.n1} {t.n2} {hx t.l01} {hx t.l12} {hx t.l20} {hx t.area}\n"
  for sg in cm.segs do
    s := s ++ s!"g {sg.p} {sg.q} {sg.face} {hx sg.len}\n"
  for rd in ev.regs do
    s := s ++ s!"r {rd.r} {rd.nTri} {rd.inS.foldl (fun n x => if x then n + 1 else n) 0} {rd.nComp} {rd.nUnsolvedComp} {rd.nUnsolvedNodes} {hx rd.J} {hx rd.area} {hx rd.areaSolved}\n"
  for rd in ev.regs do
    for n in [0:cm.nodes] do
      if rd.inS[n]! then
        s := s ++ s!"u {rd.r} {n} {b2n rd.solved[n]!} {hx (rd.u.get! n)} {hx (rd.p.get! n)}\n"
    for n in [0:cm.nodes] do
      if rd.isB[n]! then
        s := s ++ s!"b {rd.r} {n} {hx (rd.v.get! n)} {hx (rd.ell.get! n)} {hx (rd.d.get! n)} {hx (rd.q.get! n)} {hx (rd.e.get! n)} {hx (rd.G.get! n)}\n"
  for k in [0:cm.nc] do
    s := s ++ s!"s {k} {hx (ev.sRaw.get! k)} {hx (ev.sSmooth.get! k)} {b2n ev.sliver[k]!}\n"
  s := s ++ energyLine "energy" ev
  return s

/-- Output block of the `cuts` task. -/
def printOutcome (c : CutCase) (o : Outcome) : String := Id.run do
  match o with
  | .noCut b phi ev =>
    return s!"case {c.name}\ntask cuts\nstatus no_cut\n" ++ printEval c b phi ev ++ "end\n"
  | .step r =>
    let evP := match r.stage2 with | some st => { r.ev with sigma := st.sigma } | none => r.ev
    let mut s := s!"case {c.name}\ntask cuts\nstatus ok\n" ++ printEval c r.base r.phi evP
    s := s ++ s!"tau {hx r.tau} {hx r.hMin} {hx r.sMax}\n"
    if let some st := r.stage2 then
      s := s ++ s!"omega {hx st.omega}\nlimits {hx st.limits.1} {hx st.limits.2.1} {hx st.limits.2.2}\nlammax {hx st.lamMax}\n"
      for k in [0:r.ev.cm.nc] do
        s := s ++ s!"c {k} {hx (st.force.wc.get! k)} {hx (st.force.ellC.get! k)} {st.force.deg[k]!} {hx (st.force.beta.get! k)} {hx (st.force.theta.get! k)} {hx (st.lam0.get! k)} {hx (st.rate.get! k)} {hx (st.gbar.get! k)} {hx (st.rg.get! k)} {st.piece[k]!} {hx (st.omegaP.get! k)}\n"
      for si in [0:r.ev.cm.segs.size] do
        s := s ++ s!"dw {si} {hx (st.force.dwdn.get! si)}\n"
      if c.format == 3 then
        let mut nlim := 0
        match r.vlimit with
        | some (_, _, fv, _) =>
          for i in [0:r.base.nv] do
            if r.ev.cm.adj[i]! && fv.get! i < 1 then nlim := nlim + 1
        | none =>
          for k in [0:r.ev.cm.nc] do
            if st.factor.get! k < 1 then nlim := nlim + 1
        s := s ++ s!"stats {match c.stats with | .extremes => "extremes" | .robust => "robust"}\n"
        s := s ++ s!"omegamax {hx st.omegaMax}\nmhat {hx st.mhat.1} {hx st.mhat.2} {hx st.mused}\nlimited {nlim}\n"
        for k in [0:r.ev.cm.nc] do
          s := s ++ s!"h {k} {hx (st.hc.get! k)} {hx (st.mc1.get! k)} {hx (st.mc2.get! k)} {hx (st.factor.get! k)} {hx (st.rate.get! k)} {hx (st.rlim.get! k)}\n"
        if let some (hv, omv, fv, pl) := r.vlimit then
          s := s ++ "vlimit 1\n"
          for i in [0:r.base.nv] do
            s := s ++ s!"hv {i} {hx (hv.get! i)} {hx (omv.get! i)} {hx (fv.get! i)} {hx (pl.get! i)}\n"
      for i in [0:r.base.nv] do
        let rn := match st.ring[i]! with | some n => toString n | none => "-1"
        s := s ++ s!"m {i} {rn} {hx (st.sigma.get! i)}\n"
    let nAdj := r.redist.adj.foldl (fun n x => if x then n + 1 else n) 0
    for i in [0:r.base.nv] do
      s := s ++ s!"n {i} {hx (r.phiTilde.get! i)} {hx (r.phiLen.get! i)} {hx (r.phiSym.get! i)} {r.redist.region[i]!} {b2n r.redist.adj[i]!} {hx (r.redist.T.get! i)} {hx (r.phiNew.get! i)} {hx (r.phiNew.get! i / r.base.scale)}\n"
    s := s ++ s!"redist {r.redist.sweeps} {nAdj} {r.redist.unreached}\n"
    s := s ++ energyLine "energy_new" r.evNew
    s := s ++ s!"cuts_new {r.evNew.cm.nc}\n"
    return s ++ "end\n"

/-- Run one case. -/
def runCase (c : CutCase) : String :=
  match runStep c with
  | .ok o => printOutcome c o
  | .error e => printError c.name e

/-- Read cases until the end of `input`, write one block per case, flush after each. -/
partial def runStream (input output : IO.FS.Stream) : IO Nat := do
  let reader ← LineReader.new input
  let rec loop (count : Nat) : IO Nat := do
    match ← readCaseLines reader with
    | none => return count
    | some (name, parsed) =>
      let block := match parsed with
        | .error e => printError name e
        | .ok lines =>
          match parseBody name lines with
          | .error e => printError name e
          | .ok c => runCase c
      output.putStr block
      output.flush
      loop (count + 1)
  loop 0

end Ddg.Cuts
