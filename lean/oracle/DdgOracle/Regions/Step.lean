import DdgOracle.Regions.Eval

/-!
# One step for N regions (SPEC §8.4, §8.5)
-/

namespace Ddg.Regions

open Ddg Ddg.Cuts

/-- Redistancing of N functions (SPEC §8.4). -/
structure RRedist where
  /-- regions of `φ′` -/
  region : Array Nat
  /-- cut-adjacent vertices of every region -/
  adjK : Array (Array Bool)
  /-- union -/
  adj : Array Bool
  /-- after the projection of Losasso et al. -/
  projected : Array FloatArray
  /-- fast-marching values per function (`+∞` where not reached) -/
  T : Array FloatArray
  /-- the new functions -/
  phi : Array FloatArray
  /-- accepted non-seed vertices, summed over the functions -/
  accepted : Nat
  /-- vertices that some front did not reach -/
  unreached : Nat
  deriving Inhabited

/-- Single-pass fast marching on the classes `cls` with the given seeds (SPEC §3.12). -/
def marchClass (b : Base) (faces : Array Face) (cls : Array Nat) (seed : Array Bool)
    (seedVal : FloatArray) : FloatArray × Nat := Id.run do
  let nv := b.nv
  let mut T := FloatArray.ofFn' nv fun i => if seed[i]! then seedVal.get! i else posInf
  let mut accepted := seed
  let mut count := 0
  let mut go := true
  while go do
    let mut bestI := nv
    let mut bestV := posInf
    for i in [0:nv] do
      if !accepted[i]! then
        let mut v := posInf
        for fi in b.topo.vertexFaces[i]! do
          match hopfLaxCandidate b cls T faces[fi]! i with
          | some cand => if cand < v then v := cand
          | none => pure ()
        if v < bestV then
          bestV := v
          bestI := i
    if bestI < nv then
      T := T.set! bestI bestV
      accepted := accepted.set! bestI true
      count := count + 1
    else go := false
  return (T, count)

/-- Projection and redistancing of the functions `phiP` (already symmetrized). -/
def redistanceN (b : Base) (faces : Array Face) (phiP : Array FloatArray) : RRedist := Id.run do
  let nv := b.nv
  let N := phiP.size
  let region : Array Nat := Array.ofFn (n := nv) fun i => argminRegion phiP i.val
  let mut adjK : Array (Array Bool) := Array.replicate N (Array.replicate nv false)
  let mut adj : Array Bool := Array.replicate nv false
  for a in [0:nv] do
    for c in b.topo.neighbors[a]! do
      if region[a]! != region[c]! then
        adj := adj.set! a true
        adjK := adjK.set! region[a]! (adjK[region[a]!]!.set! a true)
        adjK := adjK.set! region[c]! (adjK[region[c]!]!.set! a true)
  -- projection
  let mut proj := phiP
  for i in [0:nv] do
    if adj[i]! then
      let mut lo := posInf
      let mut lo2 := posInf
      for k in [0:N] do
        let v := (phiP[k]!).get! i
        if v < lo then
          lo2 := lo
          lo := v
        else if v < lo2 then lo2 := v
      let m := projectionShift lo lo2
      for k in [0:N] do
        proj := proj.set! k ((proj[k]!).set! i (projectedValue ((phiP[k]!).get! i) m))
  -- fast marching per function
  let mut Ts : Array FloatArray := #[]
  let mut out := proj
  let mut accepted := 0
  let mut reached : Array Bool := Array.replicate nv true
  for k in [0:N] do
    let anyAdj := adjK[k]!.any id
    if !anyAdj then
      Ts := Ts.push (FloatArray.ofFn' nv fun _ => posInf)
    else
      let cls : Array Nat := Array.ofFn (n := nv) fun i => if region[i.val]! == k then 0 else 1
      let seedVal := FloatArray.ofFn' nv fun i => ((proj[k]!).get! i).abs
      let (T, cnt) := marchClass b faces cls adjK[k]! seedVal
      accepted := accepted + cnt
      Ts := Ts.push T
      let mut phik := proj[k]!
      for i in [0:nv] do
        if !adjK[k]![i]! then
          if T.get! i < posInf then
            phik := phik.set! i (if region[i]! == k then -(T.get! i) else T.get! i)
          else reached := reached.set! i false
      out := out.set! k phik
  let unreached := (reached.filter (!·)).size
  return { region, adjK, adj, projected := proj, T := Ts, phi := out, accepted, unreached }

/-- Everything of one step for N regions. -/
structure RStep where
  /-- base mesh -/
  base : Base
  /-- scaled input -/
  phis : Array FloatArray
  /-- evaluation at the input -/
  ev : NEval
  /-- step size -/
  tau : Float
  /-- shortest crossed edge -/
  hMin : Float
  /-- largest `|ρ|` (`|s|` for the paper's stages, `|s + λ⁰|` for the default stages) -/
  rhoMax : Float
  /-- roughening rate (largest over the cut points) -/
  omega : Float
  /-- limits of the step rule -/
  limits : Float × Float × Float
  /-- `φ̇_k` (paper) or the explicit rates `τ`-free per function (default) -/
  rates : Array FloatArray
  /-- `φ̃_k` -/
  phiTilde : Array FloatArray
  /-- after the length step (paper) or equal to `phiTilde` (default) -/
  phiLen : Array FloatArray
  /-- after the symmetry projection -/
  phiSym : Array FloatArray
  /-- redistancing -/
  redist : RRedist
  /-- new functions (after the second projection) -/
  phiNew : Array FloatArray
  /-- evaluation at the new functions -/
  evNew : NEval
  /-- cut motion -/
  motion : Float
  /-- default-stage data per region: local cut copies' `(ellC, beta, w, gamma, theta, lam0, rate, gbar, rho)` -/
  dflt : Array (Array (Array Float))
  /-- `dw/dn` per region per segment of its view -/
  dwdn : Array FloatArray
  /-- D26 (robust statistics, paper's stages, `Capped`/`Adaptive`): `h_i`, `ω_i` and the factor of the
  limit at every vertex -/
  vlimit : Option (FloatArray × FloatArray × FloatArray) := none
  /-- largest local roughening rate -/
  omegaMax : Float
  /-- `m̂` for `s` over the cut points, for `s + λ⁰` over the boundary nodes, the one used -/
  mhat : Float × Float × Float
  /-- `h_c` of every cut point -/
  hP : FloatArray
  /-- weight `ℓ_c` of every cut point -/
  ellP : FloatArray
  /-- local-limit factor of every boundary copy of every region -/
  factor : Array FloatArray
  /-- rate before and after the local limit -/
  rBefore : Array FloatArray
  /-- see `rBefore` -/
  rAfter : Array FloatArray
  /-- `e_i` of the default extension (at every vertex) and its ring -/
  eVert : FloatArray
  /-- ring -/
  ring : Array (Option Nat)
  /-- `δ_i^j`: list of `(j, value)` per vertex -/
  deltas : Array (List (Nat × Float))

/-- Ring-by-ring extension of values given at the cut-adjacent vertices (SPEC §3.8, second line). -/
def ringExtend (b : Base) (region : Array Nat) (adj : Array Bool) (sigma0 : FloatArray) :
    FloatArray × Array (Option Nat) := Id.run do
  let nv := b.nv
  let mut sigma := sigma0
  let mut ring : Array (Option Nat) := Array.replicate nv none
  for i in [0:nv] do
    if adj[i]! then ring := ring.set! i (some 0)
  let mut level := 0
  let mut go := true
  while go do
    let mut next : Array Nat := #[]
    for i in [0:nv] do
      if ring[i]! == some level then
        for j in b.topo.neighbors[i]! do
          if ring[j]!.isNone && region[j]! == region[i]! then
            ring := ring.set! j (some (level + 1))
            next := next.push j
    if next.size == 0 then go := false
    else
      for i in next do
        let mut nb : List (Float × Float) := []
        for j in b.topo.neighbors[i]! do
          if ring[j]! == some level && region[j]! == region[i]! then
            nb := nb ++ [(sigma.get! j, b.l i j)]
        sigma := sigma.set! i (ringMean nb)
      level := level + 1
  return (sigma, ring)

/-- Number of faces at the edge `{a, c}` is one: the edge is on `∂M`. -/
def isBoundaryEdge (b : Base) (a c : Nat) : Bool := Id.run do
  let mut n := 0
  for fi in b.topo.vertexFaces[a]! do
    let f := b.topo.faces[fi]!
    if f.a == c || f.b == c || f.c == c then n := n + 1
  return n == 1

/-- Validate the parameters of a region case. -/
def validateRegions (c : RegionCase) : Except String Unit := do
  let N := c.phis.size
  if N < 2 then throw "at least two regions are required"
  if N > 12 then throw "at most 12 regions"
  if c.ws.size != N || c.betas.size != N then throw s!"w and beta need {N} values"
  for k in [0:N] do
    if (c.phis[k]!).size != c.nv then throw s!"phi_{k} has {(c.phis[k]!).size} values for {c.nv} vertices"
    for i in [0:c.nv] do
      if !((c.phis[k]!).get! i).isFinite then throw s!"phi_{k}[{i}] is not finite"
    if !(c.ws[k]!.isFinite && c.ws[k]! ≥ 0) then throw s!"w_{k} must be finite and >= 0"
    if !(c.betas[k]!).isFinite then throw s!"beta_{k} must be finite"
  let c2 : CutCase := { c.toCut with phi := zeros c.nv }
  Ddg.Cuts.validateParams c2

/-- One step. -/
def runRegionStep (c : RegionCase) : Except String RStep := do
  let N := c.phis.size
  if N < 2 then throw "at least two regions are required"
  validateRegions c
  let b ← buildBase c.toCut
  let nv := b.nv
  let phis := c.phis.map fun p => FloatArray.ofFn' nv fun i => b.scale * p.get! i
  let ev ← evaluateNet c b phis
  let net := ev.net
  let nP := net.nP
  if nP == 0 then throw "no edge is crossed: the state has a single region"
  let mut hMin := posInf
  for cr in net.cps do
    if cr.len < hMin then hMin := cr.len
  let mut omegaMaxV := negInf
  for p in [0:nP] do
    if ev.omegaP.get! p > omegaMaxV then omegaMaxV := ev.omegaP.get! p
  -- weights of the cut points (half a segment to a cut point, a whole one to a triple point) and h_c
  let mut ellP := zeros nP
  for sg in net.segs do
    let p := sg.g1 - nv
    if sg.toTriple then ellP := ellP.set! p (ellP.get! p + sg.len)
    else
      ellP := ellP.set! p (ellP.get! p + sg.len / 2)
      let q := sg.g2 - nv
      ellP := ellP.set! q (ellP.get! q + sg.len / 2)
  let hP := FloatArray.ofFn' nP fun p => Id.run do
    let cr := net.cps[p]!
    let mut h := posInf
    for p2 in [0:nP] do
      let c2 := net.cps[p2]!
      if c2.a == cr.a || c2.a == cr.b || c2.b == cr.a || c2.b == cr.b then
        if c2.len < h then h := c2.len
    return h
  let robust := c.stats == .robust
  let omegaMax := omegaMaxV
  let omega := if robust then weightedQuantile 0.9 ev.omegaP.toList ellP.toList else omegaMaxV
  let paperRule := c.penalty == .paperDiagonal
  let cc := fluxC c.flux
  let wGeoNode := fun (rv : RView) (node : Nat) =>
    if node < nv then ev.wv.get! node
    else valueAtGeo net c.faces ev.wv (geoOfCopy net rv.toCopy[node - nv]!).1
  -- stage-dependent part
  let mut rates : Array FloatArray := Array.replicate N (zeros nv)
  let mut tau := fzero
  let mut rhoMax := fzero
  let mut limits : Float × Float × Float := (posInf, posInf, posInf)
  let mut dfl : Array (Array (Array Float)) := #[]
  let mut dwdnAll : Array FloatArray := #[]
  let mut eVert := zeros nv
  let mut ring : Array (Option Nat) := Array.replicate nv none
  let mut deltas : Array (List (Nat × Float)) := Array.replicate nv []
  let stepTau := fun (rm : Float) =>
    match c.step with
    | .fixed t => t
    | .capped tm mm => if rm == 0 then tm else stepCap tm mm hMin rm
    | .adaptive ts mm tm => adaptiveStep tm ts omega mm hMin rm
  let stepLimits := fun (rm : Float) =>
    match c.step with
    | .fixed t => (t, posInf, posInf)
    | .capped tm mm => (tm, posInf, if rm > 0 then mm * hMin / rm else posInf)
    | .adaptive ts mm tm =>
      (tm, if omega > 0 then ts / omega else posInf, if rm > 0 then mm * hMin / rm else posInf)
  -- limits and rates of §8.4 / §8.5 per region
  let needForce := c.stages == .dflt || robust
  let mut lam0s : Array FloatArray := #[]
  let mut thetas : Array FloatArray := #[]
  let mut ellCs : Array FloatArray := #[]
  let mut betas : Array FloatArray := #[]
  let mut wns : Array FloatArray := #[]
  let mut gammas : Array FloatArray := #[]
  if needForce then
    for X in [0:N] do
      let rv := ev.views[X]!
      let nl := rv.toCopy.size
      let mut ellC := zeros nl
      for sg in rv.cm.segs do
        ellC := (ellC.set! sg.p (ellC.get! sg.p + sg.len / 2)).set! sg.q (ellC.get! sg.q + sg.len / 2)
      let mut beta := zeros nl
      for t in rv.cm.tris do
        for (n, cot) in [(t.n0, t.c0), (t.n1, t.c1), (t.n2, t.c2)] do
          if n ≥ nv then beta := beta.set! (n - nv) (beta.get! (n - nv) + angleFromCot cot)
      let wn := FloatArray.ofFn' nl fun j => wGeoNode rv (nv + j)
      let gamma := FloatArray.ofFn' nl fun j =>
        let (g, _) := geoOfCopy net rv.toCopy[j]!
        if g ≥ nv + nP then 2 * Scalar.pi / 3
        else if isBoundaryEdge b net.cps[g - nv]!.a net.cps[g - nv]!.b then Scalar.pi / 2
        else Scalar.pi
      -- dw/dn out of X across every segment
      let mut dwdn := zeros rv.cm.segs.size
      for si in [0:rv.cm.segs.size] do
        let sg := rv.cm.segs[si]!
        let mut found := false
        for t in rv.cm.tris do
          if !found then
            let ns := #[t.n0, t.n1, t.n2]
            let cs := #[t.c0, t.c1, t.c2]
            let iu := ns.findIdx? (· == nv + sg.p)
            let iv := ns.findIdx? (· == nv + sg.q)
            match iu, iv with
            | some a, some bq =>
              found := true
              let iz := 3 - a - bq
              let lam := footFraction cs[a]! cs[bq]!
              let h := heightOver t.area sg.len
              dwdn := dwdn.set! si (normalDerivativeOfLinear (wn.get! sg.p) (wn.get! sg.q)
                (wGeoNode rv ns[iz]!) lam h)
            | _, _ => pure ()
      let mut theta := zeros nl
      for j in [0:nl] do
        let mut segs : List (Float × Float) := []
        for si in [0:rv.cm.segs.size] do
          let sg := rv.cm.segs[si]!
          if sg.p == j || sg.q == j then
            segs := segs ++ [(sg.len, if paperRule then fzero else dwdn.get! si)]
        theta := theta.set! j (lengthForce (wn.get! j) (gamma.get! j) (beta.get! j) segs)
      let mut edges : Array (Nat × Nat × Float) := #[]
      for sg in rv.cm.segs do
        let f := b.topo.faces[sg.face]!
        let h2 := faceMeanSquareEdge (b.l f.a f.b) (b.l f.b f.c) (b.l f.c f.a)
        edges := edges.push (sg.p, sg.q, fluxCoupling cc h2 sg.len)
      let lam0 ← if c.alphaL == 0 && c.eta.isNone then pure (zeros nl) else chainSolve nl ellC edges theta
      lam0s := lam0s.push lam0
      thetas := thetas.push theta
      ellCs := ellCs.push ellC
      betas := betas.push beta
      wns := wns.push wn
      gammas := gammas.push gamma
      dwdnAll := dwdnAll.push dwdn
  -- the local data of every boundary copy: h, local roughening rate, weight
  let copyH := fun (X j : Nat) =>
    let rv := ev.views[X]!
    let (g, _) := geoOfCopy net rv.toCopy[j]!
    if g < nv + nP then hP.get! (g - nv)
    else
      let t := net.triples[g - nv - nP]!
      min (hP.get! t.p0) (min (hP.get! t.p1) (hP.get! t.p2))
  let copyOmega := fun (X j : Nat) =>
    let rv := ev.views[X]!
    let (g, _) := geoOfCopy net rv.toCopy[j]!
    if g < nv + nP then ev.omegaP.get! (g - nv) else omega
  let nodesOf := fun (X : Nat) => (ev.views[X]!).toCopy.size
  -- speeds for the step rule
  let mut rmPaper := fzero
  for p in [0:nP] do
    if (ev.sCut.get! p).abs > rmPaper then rmPaper := (ev.sCut.get! p).abs
  let mut rmDflt := fzero
  if needForce then
    for X in [0:N] do
      for j in [0:nodesOf X] do
        let a := ((ev.sNode[X]!).get! j + (lam0s[X]!).get! j).abs
        if a > rmDflt then rmDflt := a
  let rm := if c.stages == .paper then rmPaper else rmDflt
  rhoMax := rm
  -- robust statistics m-hat
  let mut m1 := fzero
  let mut m2 := fzero
  let mut mhat := fzero
  if robust then
    m1 := weightedQuantile 0.9 ((List.range nP).map fun p => speedInEdges (ev.sCut.get! p) (hP.get! p)) ellP.toList
    -- pooled over all boundary copies; equal values in the order of the node numbers of §8.2
    let mut ents : Array (Nat × Nat × Nat) := #[]
    for X in [0:N] do
      for j in [0:nodesOf X] do
        ents := ents.push ((ev.views[X]!).toCopy[j]!, X, j)
    ents := ents.qsort (fun a b => a.1 < b.1)
    let mut xs : List Float := []
    let mut ls : List Float := []
    for (_, X, j) in ents do
      xs := xs ++ [speedInEdges ((ev.sNode[X]!).get! j + (lam0s[X]!).get! j) (copyH X j)]
      ls := ls ++ [(ellCs[X]!).get! j]
    m2 := weightedQuantile 0.9 xs ls
    mhat := if c.stages == .paper then max m1 m2 else m2
    tau := match c.step with
      | .fixed t => t
      | .capped tm mm => robustStep tm posInf 0 mm mhat
      | .adaptive ts mm tm => robustStep tm ts omega mm mhat
    limits := match c.step with
      | .fixed t => (t, posInf, posInf)
      | .capped tm mm => (tm, posInf, if mhat > 0 then mm / mhat else posInf)
      | .adaptive ts mm tm => (tm, if omega > 0 then ts / omega else posInf, if mhat > 0 then mm / mhat else posInf)
  else
    tau := stepTau rm
    limits := stepLimits rm
  -- the local limit
  let limitFactor := fun (X j : Nat) (r : Float) =>
    if !robust || c.stages == .paper then 1.0
    else match c.step with
      | .fixed _ => 1.0
      | .capped _ mm => if r.abs / copyH X j ≤ mhat then 1.0 else localLimitFactor mm (copyH X j) tau r.abs
      | .adaptive _ mm _ =>
        let v := localSpeed r (copyOmega X j) omega
        if v / copyH X j ≤ mhat then 1.0 else localLimitFactor mm (copyH X j) tau v
  let sLim : Array FloatArray := Array.ofFn (n := N) fun X =>
    FloatArray.ofFn' (nodesOf X.val) fun j => (ev.sNode[X.val]!).get! j * limitFactor X.val j ((ev.sNode[X.val]!).get! j)
  let mut factorAll : Array FloatArray := Array.ofFn (n := N) fun X =>
    FloatArray.ofFn' (nodesOf X.val) fun j => limitFactor X.val j ((ev.sNode[X.val]!).get! j)
  let mut rBefore : Array FloatArray := ev.sNode
  let mut rAfter : Array FloatArray := sLim
  if c.stages == .paper then
    -- §8.4
    for X in [0:N] do
      let rv := ev.views[X]!
      let e ← harmonicExtension rv.cm #[rv.rd] (sLim[X]!)
      for i in [0:nv] do
        if net.region[i]! == X then
          for k in [0:N] do
            rates := rates.set! k ((rates[k]!).set! i (if k == X then e.get! i else -(e.get! i)))
  else
    -- §8.5, per region
    -- rates per region
    let allZero := Id.run do
      let mut z := true
      for i in [0:nv] do
        if ev.wv.get! i != 0 then z := false
      return z
    let mut rhoV : Array FloatArray := #[]
    let mut rateV : Array FloatArray := #[]
    for X in [0:N] do
      let rv := ev.views[X]!
      let nl := rv.toCopy.size
      let rate ←
        if allZero then pure (ev.sNode[X]!)
        else do
          let mut edges : Array (Nat × Nat × Float) := #[]
          for sg in rv.cm.segs do
            let wbar := segmentMeanWeight ((wns[X]!).get! sg.p) ((wns[X]!).get! sg.q)
            edges := edges.push (sg.p, sg.q, rateCoupling tau wbar sg.len)
          let rhs := FloatArray.ofFn' nl fun j => rateRhs ((ellCs[X]!).get! j) ((ev.sNode[X]!).get! j) ((lam0s[X]!).get! j)
          chainSolve nl (ellCs[X]!) edges rhs
      rateV := rateV.push rate
    -- pieces of the network over the geometric ids (cut points, then triple points)
    let nG := nP + net.nM
    let gid := fun (g : Nat) => g - nv
    let mut lab : Array Nat := Array.range nG
    let mut changed := true
    while changed do
      changed := false
      for s in net.segs do
        let a := gid s.g1
        let bb := gid s.g2
        let m := min lab[a]! lab[bb]!
        if lab[a]! != m then
          lab := lab.set! a m
          changed := true
        if lab[bb]! != m then
          lab := lab.set! bb m
          changed := true
    let mut sumL := zeros nG
    let mut sumLG := zeros nG
    for s in net.segs do
      if !s.toTriple then
        let t := net.tris[s.tri]!
        let i := net.region[s.lone]!
        let j := if s.regA == i then s.regB else s.regA
        let dphi := (((phis[i]!).get! s.lone) - ((phis[j]!).get! s.lone)).abs / 2
        let g := gradientNormLone dphi s.len t.area
        let r := lab[gid s.g1]!
        -- a face with |grad phi| = 0 is left out of both sums (SPEC §3.9, D23)
        if g > 0 then
          sumL := sumL.set! r (sumL.get! r + s.len)
          sumLG := sumLG.set! r (sumLG.get! r + s.len / g)
    let gbarOf := fun (g : Nat) =>
      let r := lab[gid g]!
      if sumL.get! r > 0 then meanGradient (sumL.get! r) (sumLG.get! r) else 1.0
    -- rho at the copies, scaled by gbar
    for X in [0:N] do
      let rv := ev.views[X]!
      let nl := rv.toCopy.size
      let gb := FloatArray.ofFn' nl fun j => gbarOf (geoOfCopy net rv.toCopy[j]!).1
      let fac := FloatArray.ofFn' nl fun j => limitFactor X j ((rateV[X]!).get! j)
      factorAll := factorAll.set! X fac
      rBefore := rBefore.set! X (rateV[X]!)
      let rlim := FloatArray.ofFn' nl fun j => (rateV[X]!).get! j * fac.get! j
      rAfter := rAfter.set! X rlim
      let rho := FloatArray.ofFn' nl fun j => gb.get! j * rlim.get! j
      rhoV := rhoV.push rho
      dfl := dfl.push (Array.ofFn (n := nl) fun j =>
        #[(ellCs[X]!).get! j.val, (betas[X]!).get! j.val, (wns[X]!).get! j.val, (gammas[X]!).get! j.val,
          (thetas[X]!).get! j.val, (lam0s[X]!).get! j.val, (rateV[X]!).get! j.val, gb.get! j.val, rho.get! j.val])
    let mut eW : Array (List Float) := Array.replicate nv []
    let mut eV : Array (List Float) := Array.replicate nv []
    let mut dW : Array (List (Nat × Float × Float)) := Array.replicate nv []
    for p in [0:nP] do
      let cr := net.cps[p]!
      let ra := net.region[cr.a]!
      let rb := net.region[cr.b]!
      let rhoA := (rhoV[ra]!).get! ((ev.views[ra]!).localOf.getD (net.copyIdx (nv + p) ra) 0)
      let rhoB := (rhoV[rb]!).get! ((ev.views[rb]!).localOf.getD (net.copyIdx (nv + p) rb) 0)
      let wa := ellP.get! p * (1 - cr.t)
      let wb := ellP.get! p * cr.t
      eW := eW.set! cr.a (eW[cr.a]! ++ [wa])
      eV := eV.set! cr.a (eV[cr.a]! ++ [rhoA])
      dW := dW.set! cr.a (dW[cr.a]! ++ [(rb, wa, rhoA - rhoB)])
      eW := eW.set! cr.b (eW[cr.b]! ++ [wb])
      eV := eV.set! cr.b (eV[cr.b]! ++ [rhoB])
      dW := dW.set! cr.b (dW[cr.b]! ++ [(ra, wb, rhoB - rhoA)])
    let mut e0 := zeros nv
    for i in [0:nv] do
      if net.adj[i]! then
        e0 := e0.set! i (weightedMean eW[i]! eV[i]!)
        -- delta per region
        let mut regs : List Nat := []
        for (j, _, _) in dW[i]! do
          if !regs.contains j then regs := regs ++ [j]
        let mut dl : List (Nat × Float) := []
        for j in regs.toArray.qsort (· < ·) do
          let items := dW[i]!.filter (·.1 == j)
          dl := dl ++ [(j, weightedMean (items.map (·.2.1)) (items.map (·.2.2)))]
        deltas := deltas.set! i dl
    let (ex, rn) := ringExtend b net.region net.adj e0
    eVert := ex
    ring := rn
    for i in [0:nv] do
      let k := net.region[i]!
      for j in [0:N] do
        let v :=
          if j == k then ex.get! i
          else match (deltas[i]!).find? (·.1 == j) with
            | some (_, d) => ex.get! i - d
            | none => -(ex.get! i)
        rates := rates.set! j ((rates[j]!).set! i v)
  -- explicit step
  let phiTilde := Array.ofFn (n := N) fun k =>
    FloatArray.ofFn' nv fun i => explicitStep ((phis[k]!).get! i) tau ((rates[k]!).get! i)
  -- length step
  let mut phiLen : Array FloatArray := #[]
  for k in [0:N] do
    if c.stages == .paper then
      phiLen := phiLen.push (← lengthStepMesh c.toCut b ev.wv (phiTilde[k]!) tau)
    else phiLen := phiLen.push (phiTilde[k]!)
  -- D26: the local limit of the paper's stages acts at the vertices on the whole step
  let mut vlimit : Option (FloatArray × FloatArray × FloatArray) := none
  let mut phiLim := phiLen
  if robust && c.stages == .paper then
    match c.step with
    | .fixed _ => pure ()
    | step =>
      let mm := match step with | .capped _ m => m | .adaptive _ m _ => m | .fixed _ => 1.0
      let mut hv := zeros nv
      let mut omv := zeros nv
      let mut fv := zeros nv
      for i in [0:nv] do
        let mut h := posInf
        let mut om := omega
        let mut anyCross := false
        for p in [0:nP] do
          let cr := net.cps[p]!
          if cr.a == i || cr.b == i then
            anyCross := true
            if cr.len < h then h := cr.len
            if ev.omegaP.get! p > om then om := ev.omegaP.get! p
        if !anyCross then
          for j in b.topo.neighbors[i]! do
            if b.l i j < h then h := b.l i j
        -- a = region of i, b = the smallest of the other functions (ties: lowest index), before the step
        let a := net.region[i]!
        let mut bb := if a == 0 then 1 else 0
        for k in [0:N] do
          if k != a && (phis[k]!).get! i < (phis[bb]!).get! i then bb := k
        let delta := vertexChange ((phiLen[bb]!).get! i) ((phis[bb]!).get! i) ((phiLen[a]!).get! i) ((phis[a]!).get! i)
        let v := match step with
          | .capped _ _ => vertexSpeedCapped delta tau
          | _ => vertexSpeedAdaptive delta tau om omega
        let f := vertexLimitFactor mm h tau v
        hv := hv.set! i h
        omv := omv.set! i om
        fv := fv.set! i f
        for k in [0:N] do
          phiLim := phiLim.set! k ((phiLim[k]!).set! i (limitedValue ((phis[k]!).get! i) ((phiLen[k]!).get! i) f))
      vlimit := some (hv, omv, fv)
  let phiSym := phiLim.map (symmetrize c.classes)
  let redist := redistanceN b c.faces phiSym
  let phiNew := redist.phi.map (symmetrize c.classes)
  let evNew ← evaluateNet c b phiNew
  -- cut motion
  let mut motion := fzero
  for i in [0:nv] do
    if net.adj[i]! || evNew.net.adj[i]! then
      for k in [0:N] do
        let d := ((phiNew[k]!).get! i - (phis[k]!).get! i).abs
        if d > motion then motion := d
  return { base := b, phis, ev, tau, hMin, rhoMax, omega, limits, rates, phiTilde, phiLen, phiSym
           redist, phiNew, evNew, motion, dflt := dfl, dwdn := dwdnAll, eVert, ring, deltas
           vlimit, omegaMax, mhat := (m1, m2, mhat), hP, ellP, factor := factorAll, rBefore, rAfter }

end Ddg.Regions
