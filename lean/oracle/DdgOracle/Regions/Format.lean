import DdgOracle.Regions.Gen
import DdgOracle.Cuts.Format

/-!
# Text format of the `regions` task

Specified in `README.md` ("N regions"). Same conventions as the other tasks.
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

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

/-- Parse the lines of one case. -/
def parseBody (name : String) (lines : Array (Array String)) : Except String RegionCase := do
  let mut i := 0
  let line := fun (k : Nat) (what : String) =>
    if h : k < lines.size then pure (f := Except String) lines[k]
    else throw s!"missing line: expected {what}"
  let l ← line i "'task regions' or 'mesh <nv> <nf> <ne>'"
  if l[0]! == "task" then
    if l.size != 2 || l[1]! != "regions" then
      throw s!"malformed line '{joinTokens l}': expected 'task regions'"
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
  let l ← line i "'regions <n>'"
  if l[0]! != "regions" || l.size != 2 then throw s!"malformed line '{joinTokens l}': expected 'regions <n>'"
  let N ← natTok "regions" l[1]!
  if N < 2 || N > 12 then throw "regions: expected 2 <= n <= 12"
  i := i + 1
  let mut phis : Array FloatArray := Array.replicate N (zeros 0)
  let mut seen : Array Bool := Array.replicate N false
  let mut kv : Std.HashMap String (Array String) := {}
  while i < lines.size do
    let l := lines[i]!
    if l[0]! == "phi" then
      if l.size < 2 then throw "phi: expected 'phi <k> <values>'"
      let k ← natTok "phi" l[1]!
      if k ≥ N then throw s!"phi {k}: only {N} regions"
      if seen[k]! then throw s!"phi {k} occurs twice"
      seen := seen.set! k true
      phis := phis.set! k (← floatsFrom s!"phi {k}" l 2 nv)
    else
      if kv.contains l[0]! then throw s!"line '{l[0]!}' occurs twice"
      kv := kv.insert l[0]! l
    i := i + 1
  for k in [0:N] do
    if !seen[k]! then throw s!"missing line 'phi {k}'"
  let allowed := ["energy", "curvature", "alpha", "eta", "penalty", "w", "beta", "step",
    "delta", "flux", "stages", "featurelap", "classes", "format", "statistics"]
  for (k, _) in kv.toList do
    if !allowed.contains k then throw s!"unknown line '{k}'"
  let format ← match kv.get? "format" with
    | none => pure 1
    | some l =>
      if l.size != 2 || l[1]! != "2" then throw "format: expected 'format 2'"
      pure 2
  if format == 1 && kv.contains "statistics" then throw "line 'statistics' needs 'format 2'"
  let need := fun (k : String) =>
    match kv.get? k with
    | some l => pure (f := Except String) l
    | none => throw s!"missing line '{k}'"
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
  let ws ← do let l ← need "w"; pure (← floatsFrom "w" l 1 N).toList.toArray
  let betas ← do let l ← need "beta"; pure (← floatsFrom "beta" l 1 N).toList.toArray
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
      if l.size != 5 then throw "step adaptive: expected three values"
      pure (StepRule.adaptive (← floatTok "step" l[2]!) (← floatTok "step" l[3]!) (← floatTok "step" l[4]!))
    | _ => throw s!"malformed line '{joinTokens l}': step is fixed | capped | adaptive"
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
  let stages ← do
    let l ← need "stages"
    match l[1]? with
    | some "paper" => pure Stages.paper
    | some "default" => pure Stages.dflt
    | _ => throw s!"malformed line '{joinTokens l}': stages is paper | default"
  let featPos ← do
    let l ← need "featurelap"
    match l[1]? with
    | some "cotan" => pure false
    | some "positive" => pure true
    | _ => throw "featurelap: expected cotan | positive"
  let classes ← match kv.get? "classes" with
    | some l =>
      if l.size != nv + 1 then throw s!"classes: expected {nv} values"
      let mut out : Array Nat := Array.mkEmpty nv
      for j in [1:nv + 1] do out := out.push (← natTok "classes" l[j]!)
      pure (some out)
    | none => pure none
  let stats ← match kv.get? "statistics" with
    | none => pure StepStats.extremes
    | some l =>
      match l[1]? with
      | some "extremes" => pure StepStats.extremes
      | some "robust" => pure StepStats.robust
      | _ => throw "statistics: expected extremes | robust"
  return { name, nv, faces, edges, phis, energy, curv, alphaL, eta, penalty, ws, betas, step, delta,
           flux, stages, featPos, classes, format, stats }

/-! ## Printing -/

private def hx (x : Float) : String := floatToHex x

private def b2n (b : Bool) : String := if b then "1" else "0"

/-- A case in the input format. -/
def printCase (c : RegionCase) : String := Id.run do
  let mut s := s!"case {c.name}\ntask regions\nmesh {c.nv} {c.faces.size} {c.edges.size}\n"
  for f in c.faces do s := s ++ s!"f {f.a} {f.b} {f.c}\n"
  for (a, b, l) in c.edges do s := s ++ s!"e {a} {b} {hx l}\n"
  let arr := fun (x : FloatArray) => " ".intercalate (x.toList.map hx)
  s := s ++ s!"regions {c.phis.size}\n"
  for k in [0:c.phis.size] do s := s ++ s!"phi {k} {arr c.phis[k]!}\n"
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
  s := s ++ s!"w {" ".intercalate (c.ws.toList.map hx)}\nbeta {" ".intercalate (c.betas.toList.map hx)}\n"
  match c.step with
  | .fixed t => s := s ++ s!"step fixed {hx t}\n"
  | .capped t m => s := s ++ s!"step capped {hx t} {hx m}\n"
  | .adaptive ts m tm => s := s ++ s!"step adaptive {hx ts} {hx m} {hx tm}\n"
  s := s ++ (match c.flux with | .nodal => "flux nodal\n" | .smoothed k => s!"flux smoothed {hx k}\n")
  s := s ++ s!"delta {hx c.delta}\n"
  s := s ++ s!"stages {match c.stages with | .paper => "paper" | .dflt => "default"}\n"
  s := s ++ s!"featurelap {if c.featPos then "positive" else "cotan"}\n"
  if c.format == 2 then
    s := s ++ s!"format 2\nstatistics {match c.stats with | .extremes => "extremes" | .robust => "robust"}\n"
  if let some cl := c.classes then
    s := s ++ s!"classes {" ".intercalate (cl.toList.map toString)}\n"
  return s ++ "end\n"

/-- Output block of a case that could not be processed. -/
def printError (name msg : String) : String :=
  s!"case {name}\ntask regions\nstatus error {msg}\nend\n"

private def energyLine (tag : String) (ev : NEval) : String :=
  s!"{tag} {hx ev.total} {hx ev.lengthTerm}\n"

/-- Output block of a step. -/
def printStep (c : RegionCase) (r : RStep) : String := Id.run do
  let b := r.base
  let nv := b.nv
  let N := c.phis.size
  let ev := r.ev
  let net := ev.net
  let nP := net.nP
  let mut moved := 0
  for cr in net.cps do if cr.moved then moved := moved + 1
  let mut sumOmega := fzero
  for o in b.omega do sumOmega := sumOmega + o
  let globalNode := fun (rv : RView) (node : Nat) => if node < nv then node else rv.toCopy[node - nv]!
  let mut s := s!"case {c.name}\ntask regions\nstatus ok\n"
  s := s ++ s!"scale {hx b.scale} {hx b.totalArea}\n"
  s := s ++ s!"mesh {nv} {c.faces.size} {b.topo.numEdges} {b.topo.numBoundary} {b.topo.chi}\n"
  s := s ++ s!"regions {N}\n"
  s := s ++ s!"counts {nP} {net.nM} {net.segs.size} {net.tris.size} {moved} {net.nearTies}\n"
  s := s ++ s!"sum_omega {hx sumOmega}\n"
  for fi in [0:c.faces.size] do
    s := s ++ s!"f {fi} {hx b.omega[fi]!} {hx b.kGauss[fi]!} {hx b.kUsed[fi]!} {hx b.faceArea[fi]!}\n"
  for i in [0:nv] do
    s := s ++ s!"v {i} {net.region[i]!} {hx (b.theta.get! i)} {hx (b.vertArea.get! i)} {hx (b.feat.get! i)} {hx (b.aVert.get! i)} {hx (ev.wv.get! i)}\n"
  for i in [0:nv] do
    s := s ++ s!"p {i} {" ".intercalate ((r.phis.map fun p => hx (p.get! i)).toList)}\n"
  for p in [0:nP] do
    let cr := net.cps[p]!
    s := s ++ s!"x {p} {cr.a} {cr.b} {hx cr.tRaw} {hx cr.t} {b2n cr.moved} {hx cr.len} {ev.lowHigh[p]!.1} {ev.lowHigh[p]!.2}\n"
  for m in [0:net.nM] do
    let t := net.triples[m]!
    s := s ++ s!"q {m} {t.face} {t.c0} {t.c1} {t.c2} {t.p0} {t.p1} {t.p2} {hx t.b0} {hx t.b1} {hx t.b2}\n"
  let gid := fun (g : Nat) => g - nv
  for sg in net.segs do
    s := s ++ s!"g {gid sg.g1} {gid sg.g2} {sg.face} {hx sg.len} {sg.regA} {sg.regB} {b2n sg.toTriple}\n"
  for ti in [0:net.tris.size] do
    let t := net.tris[ti]!
    s := s ++ s!"t {ti} {t.parent} {t.region} {t.n0} {t.n1} {t.n2} {hx t.l01} {hx t.l12} {hx t.l20} {hx t.area}\n"
  for rv in ev.views do
    let rd := rv.rd
    s := s ++ s!"r {rv.X} {rd.nTri} {rd.inS.foldl (fun n x => if x then n + 1 else n) 0} {rd.nComp} {rd.nUnsolvedComp} {rd.nUnsolvedNodes} {hx rd.J} {hx rd.area} {hx rd.areaSolved}\n"
  for rv in ev.views do
    let rd := rv.rd
    for n in [0:rv.cm.nodes] do
      if rd.inS[n]! then
        s := s ++ s!"u {rv.X} {globalNode rv n} {b2n rd.solved[n]!} {hx (rd.u.get! n)} {hx (rd.p.get! n)}\n"
    for n in [0:rv.cm.nodes] do
      if rd.isB[n]! then
        s := s ++ s!"b {rv.X} {globalNode rv n} {hx (rd.v.get! n)} {hx (rd.ell.get! n)} {hx (rd.d.get! n)} {hx (rd.q.get! n)} {hx (rd.e.get! n)} {hx (rd.G.get! n)}\n"
  for p in [0:nP] do
    s := s ++ s!"s {p} {hx (ev.sCut.get! p)} {hx (ev.omegaP.get! p)}\n"
  for rv in ev.views do
    for j in [0:rv.toCopy.size] do
      s := s ++ s!"sn {rv.X} {rv.toCopy[j]!} {hx ((ev.sNode[rv.X]!).get! j)}\n"
  s := s ++ energyLine "energy" ev
  for X in [0:N] do
    s := s ++ s!"er {X} {hx ev.Js[X]!} {hx ev.As[X]!}\n"
  s := s ++ s!"tau {hx r.tau} {hx r.hMin} {hx r.rhoMax}\n"
  s := s ++ s!"omega {hx r.omega}\nlimits {hx r.limits.1} {hx r.limits.2.1} {hx r.limits.2.2}\n"
  if c.format == 2 then
    let mut nlim := 0
    match r.vlimit with
    | some (_, _, fv) =>
      for i in [0:nv] do
        if net.adj[i]! && fv.get! i < 1 then nlim := nlim + 1
    | none =>
      for X in [0:N] do
        for j in [0:(r.factor[X]!).size] do
          if (r.factor[X]!).get! j < 1 then nlim := nlim + 1
    s := s ++ s!"stats {match c.stats with | .extremes => "extremes" | .robust => "robust"}\n"
    s := s ++ s!"omegamax {hx r.omegaMax}\nmhat {hx r.mhat.1} {hx r.mhat.2.1} {hx r.mhat.2.2}\nlimited {nlim}\n"
    for p in [0:nP] do
      s := s ++ s!"hp {p} {hx (r.hP.get! p)} {hx (r.ellP.get! p)} {hx (ev.omegaP.get! p)}\n"
    for rv in ev.views do
      for j in [0:rv.toCopy.size] do
        let (g, _) := geoOfCopy net rv.toCopy[j]!
        let h := if g < nv + nP then r.hP.get! (g - nv)
          else
            let t := net.triples[g - nv - nP]!
            min (r.hP.get! t.p0) (min (r.hP.get! t.p1) (r.hP.get! t.p2))
        s := s ++ s!"hc {rv.X} {rv.toCopy[j]!} {hx h} {hx ((r.factor[rv.X]!).get! j)} {hx ((r.rBefore[rv.X]!).get! j)} {hx ((r.rAfter[rv.X]!).get! j)}\n"
  if c.format == 2 then
    if let some (hv, omv, fv) := r.vlimit then
      s := s ++ "vlimit 1\n"
      for i in [0:nv] do
        s := s ++ s!"hv {i} {hx (hv.get! i)} {hx (omv.get! i)} {hx (fv.get! i)}\n"
  for k in [0:N] do
    for i in [0:nv] do
      s := s ++ s!"pd {k} {i} {hx ((r.rates[k]!).get! i)}\n"
  if c.stages == .dflt then
    for rv in ev.views do
      for j in [0:rv.toCopy.size] do
        let d := (r.dflt[rv.X]!)[j]!
        s := s ++ s!"d {rv.X} {rv.toCopy[j]!} {" ".intercalate (d.toList.map hx)}\n"
      let dw := r.dwdn[rv.X]!
      for si in [0:rv.cm.segs.size] do
        s := s ++ s!"dw {rv.X} {si} {hx (dw.get! si)}\n"
    for i in [0:nv] do
      let rn := match r.ring[i]! with | some n => toString n | none => "-1"
      s := s ++ s!"ev {i} {hx (r.eVert.get! i)} {rn}\n"
      for (j, d) in r.deltas[i]! do
        s := s ++ s!"ed {i} {j} {hx d}\n"
  for k in [0:N] do
    for i in [0:nv] do
      s := s ++ s!"n {k} {i} {hx ((r.phiTilde[k]!).get! i)} {hx ((r.phiLen[k]!).get! i)} {hx ((r.phiSym[k]!).get! i)} {hx ((r.redist.projected[k]!).get! i)} {hx ((r.redist.T[k]!).get! i)} {hx ((r.phiNew[k]!).get! i)} {hx ((r.phiNew[k]!).get! i / b.scale)}\n"
  for i in [0:nv] do
    let bits := String.join ((List.range N).map fun k => b2n (r.redist.adjK[k]!)[i]!)
    s := s ++ s!"a {i} {r.redist.region[i]!} {b2n r.redist.adj[i]!} {bits}\n"
  s := s ++ s!"redist {r.redist.accepted} {r.redist.unreached}\n"
  s := s ++ s!"motion {hx r.motion}\n"
  s := s ++ energyLine "energy_new" r.evNew
  s := s ++ s!"cuts_new {r.evNew.net.nP} {r.evNew.net.nM}\n"
  return s ++ "end\n"

/-- Run one case. -/
def runCase (c : RegionCase) : String :=
  match runRegionStep c with
  | .ok r => printStep c r
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

end Ddg.Regions
