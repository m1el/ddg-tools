import DdgOracle.Cuts.Stages

/-!
# One flow step (SPEC §3.7 – §3.12, §4, §6)
-/

namespace Ddg.Cuts

/-- Evaluation at a level set function: cut mesh, both regions, speed, extension, energy. -/
structure Eval where
  /-- the cut mesh -/
  cm : CutMesh
  /-- weight `a` at every node (empty unless the energy is scaled) -/
  aNode : FloatArray
  /-- `w_i = α_L + η_i` at every original vertex -/
  wv : FloatArray
  /-- the two regions -/
  regs : Array RegionData
  /-- speed at every cut vertex before smoothing -/
  sRaw : FloatArray
  /-- sliver vertex flags -/
  sliver : Array Bool
  /-- speed after smoothing -/
  sSmooth : FloatArray
  /-- harmonic extension on the original vertices -/
  sigma : FloatArray
  /-- energy terms: `J_0`, `J_1`, areas (solved parts), length term, total -/
  J0 : Float
  /-- see `J0` -/
  J1 : Float
  /-- see `J0` -/
  A0 : Float
  /-- see `J0` -/
  A1 : Float
  /-- see `J0` -/
  lengthTerm : Float
  /-- see `J0` -/
  total : Float
  /-- number of triangles whose aspect ratio is above `ρ₀` -/
  nSliverTris : Nat

private def wAt (c : CutCase) (b : Base) (cm : CutMesh) (wv : FloatArray) (node : Nat) : Float :=
  if node < cm.nv then wv.get! node
  else
    let cr := cm.crossings[node - cm.nv]!
    (1 - cr.t) * wv.get! cr.a + cr.t * wv.get! cr.b

/-- Evaluate everything that does not need `τ` at the (scaled) level set function `phi`. -/
def evaluate (c : CutCase) (b : Base) (phi : FloatArray) : Except String Eval := do
  let cm := buildCutMesh b c.faces phi c.delta
  let nv := b.nv
  let N := cm.nodes
  let nc := cm.nc
  for t in cm.tris do
    if !(t.area > 0) then
      throw s!"sub-triangle of face {t.parent} has no area (a crossing on a vertex; delta = 0?)"
  let wv := FloatArray.ofFn' nv fun i =>
    c.alphaL + (match c.eta with | some e => e.get! i | none => fzero)
  -- weight a at nodes
  let mut aNode := zeros 0
  if c.energy.isScaled then
    aNode := zeros N
    for i in [0:nv] do aNode := aNode.set! i (b.aVert.get! i)
    for k in [0:nc] do
      let cr := cm.crossings[k]!
      aNode := aNode.set! (nv + k) ((1 - cr.t) * b.aVert.get! cr.a + cr.t * b.aVert.get! cr.b)
  let r0 ← computeRegion b cm c.energy aNode b.kUsed c.flux 0
  let r1 ← computeRegion b cm c.energy aNode b.kUsed c.flux 1
  -- speed
  let mut sRaw := zeros nc
  for k in [0:nc] do
    let g0 := regionDensity c.w.1 (r0.G.get! (nv + k)) c.beta.1
    let g1 := regionDensity c.w.2 (r1.G.get! (nv + k)) c.beta.2
    sRaw := sRaw.set! k (cutSpeed g0 g1)
  -- sliver vertices
  let mut sliver : Array Bool := Array.replicate nc false
  let mut nSliverTris := 0
  for t in cm.tris do
    let longest := max t.l01 (max t.l12 t.l20)
    if aspectRatio longest t.area > c.rho0 then
      nSliverTris := nSliverTris + 1
      for n in [t.n0, t.n1, t.n2] do
        if n ≥ nv then sliver := sliver.set! (n - nv) true
  -- neighbours along the cut
  let mut nbrs : Array (List (Nat × Float)) := Array.replicate nc []
  for s in cm.segs do
    nbrs := (nbrs.set! s.p ((s.q, s.len) :: nbrs[s.p]!)).set! s.q ((s.p, s.len) :: nbrs[s.q]!)
  let mut sSmooth := sRaw
  for _ in [0:c.passes] do
    let prev := sSmooth
    for k in [0:nc] do
      if sliver[k]! then
        match nbrs[k]! with
        | [(n1, l1), (n2, l2)] =>
          sSmooth := sSmooth.set! k (smoothedSpeed l1 l2 (prev.get! n1) (prev.get! n2))
        | _ => pure ()
  -- harmonic extension
  let mut sigma := zeros nv
  for rd in [r0, r1] do
    let mut groups : Std.HashMap Nat (Array Nat) := {}
    for i in [0:N] do
      if rd.inS[i]! then
        groups := groups.insert rd.comp[i]! ((groups.getD rd.comp[i]! #[]).push i)
    for (_, members) in groups.toList do
      let O := members.filter (· < nv)
      let C := members.filter (· ≥ nv)
      if C.size > 0 && O.size > 0 then
        let rhs := FloatArray.ofFn' O.size fun k => Id.run do
          let mut sum := fzero
          for cn in C do
            sum := sum + rd.lap.get O[k]! cn * sSmooth.get! (cn - nv)
          return -sum
        let some sg := rd.lap.solveSub O rhs | throw s!"region {rd.r}: singular extension system"
        for k in [0:O.size] do sigma := sigma.set! O[k]! (sg.get! k)
  -- energy
  let mut lengthTerm := fzero
  for s in cm.segs do
    lengthTerm := lengthTerm + (wAt c b cm wv (nv + s.p) + wAt c b cm wv (nv + s.q)) / 2 * s.len
  let total := c.w.1 * r0.J + c.beta.1 * r0.areaSolved + c.w.2 * r1.J + c.beta.2 * r1.areaSolved
    + lengthTerm
  return { cm, aNode, wv, regs := #[r0, r1], sRaw, sliver, sSmooth, sigma
           J0 := r0.J, J1 := r1.J, A0 := r0.areaSolved, A1 := r1.areaSolved
           lengthTerm, total, nSliverTris }

/-- Replace every value by the mean over its class (the members in index order). -/
def symmetrize (classes : Option (Array Nat)) (phi : FloatArray) : FloatArray :=
  match classes with
  | none => phi
  | some cl => Id.run do
    let n := phi.size
    let mut out := phi
    for i in [0:n] do
      let mut vals : List Float := []
      let mut cnt := 0
      for j in [0:n] do
        if cl[j]! == cl[i]! then
          vals := phi.get! j :: vals
          cnt := cnt + 1
      -- `vals` is in reverse index order; sum it in index order
      out := out.set! i (classMean vals.reverse (Float.ofNat cnt))
    return out

/-- Result of the redistancing. -/
structure Redist where
  /-- region of every vertex of `φ′` -/
  region : Array Nat
  /-- cut-adjacent vertices of `φ′` -/
  adj : Array Bool
  /-- the Hopf–Lax values (`+∞` where no front arrives; `|φ′|` on `adj`) -/
  T : FloatArray
  /-- the new values -/
  phi : FloatArray
  /-- sweeps until nothing changed -/
  sweeps : Nat
  /-- vertices outside `adj` that no front reached -/
  unreached : Nat

/-- The two other vertices of the face `f` besides `i`, in the cyclic order of `f`. -/
def otherTwo (f : Face) (i : Nat) : Nat × Nat :=
  if f.a == i then (f.b, f.c) else if f.b == i then (f.c, f.a) else (f.a, f.b)

/-- The Hopf–Lax candidate for vertex `i` from the face `f`, given the current values `T` and
the regions; `none` if neither of the two other vertices is in `i`'s sign class with a finite
value. -/
def hopfLaxCandidate (b : Base) (region : Array Nat) (T : FloatArray) (f : Face) (i : Nat) :
    Option Float :=
  let (j, k) := otherTwo f i
  let okJ := region[j]! == region[i]! && T.get! j < posInf
  let okK := region[k]! == region[i]! && T.get! k < posInf
  if okJ && okK then some (hopfLaxTriangle (b.l i j) (b.l i k) (b.l j k) (T.get! j) (T.get! k))
  else if okJ then some (hopfLaxEdge (T.get! j) (b.l i j))
  else if okK then some (hopfLaxEdge (T.get! k) (b.l i k))
  else none

/-- Seeds, regions and the initial `T` of the redistancing. -/
def redistSetup (b : Base) (phiP : FloatArray) : Array Nat × Array Bool × FloatArray := Id.run do
  let nv := b.nv
  let region : Array Nat := Array.ofFn (n := nv) fun i => regionOf (phiP.get! i.val)
  let mut adj : Array Bool := Array.replicate nv false
  for a in [0:nv] do
    for c in b.topo.neighbors[a]! do
      if region[a]! != region[c]! then adj := adj.set! a true
  let T := FloatArray.ofFn' nv fun i => if adj[i]! then (phiP.get! i).abs else posInf
  return (region, adj, T)

/-- The result from the final `T`. -/
def redistFinish (nv : Nat) (phiP : FloatArray) (region : Array Nat) (adj : Array Bool)
    (T : FloatArray) (sweeps : Nat) : Redist := Id.run do
  let mut phi := phiP
  let mut unreached := 0
  for i in [0:nv] do
    if !adj[i]! then
      if T.get! i < posInf then
        phi := phi.set! i (if region[i]! == 0 then -(T.get! i) else T.get! i)
      else unreached := unreached + 1
  return { region, adj, T, phi, sweeps, unreached }

/-- Redistancing of SPEC §3.12: fast marching in a single pass. All cut-adjacent vertices are
accepted seeds; repeatedly the not yet accepted vertex with the smallest tentative value
(ties: lower index) is accepted. `sweeps` reports the number of accepted vertices. -/
def redistance (b : Base) (faces : Array Face) (phiP : FloatArray) : Redist := Id.run do
  let nv := b.nv
  let (region, adj, T0) := redistSetup b phiP
  let mut T := T0
  let mut accepted := adj
  let mut count := 0
  let mut go := true
  while go do
    let mut bestI := nv
    let mut bestV := posInf
    for i in [0:nv] do
      if !accepted[i]! then
        let mut v := posInf
        for fi in b.topo.vertexFaces[i]! do
          match hopfLaxCandidate b region T faces[fi]! i with
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
  return redistFinish nv phiP region adj T count

/-- The fixed point of the upwind scheme by Gauss–Seidel sweeps from `+∞` (the earlier
definition of SPEC §3.12). It coincides with `redistance` on meshes without obtuse angles; kept
as a cross-check. -/
def redistanceFixedPoint (b : Base) (faces : Array Face) (phiP : FloatArray) : Redist := Id.run do
  let nv := b.nv
  let (region, adj, T0) := redistSetup b phiP
  let mut T := T0
  let mut sweeps := 0
  let mut changed := true
  while changed do
    changed := false
    sweeps := sweeps + 1
    for i in [0:nv] do
      if !adj[i]! then
        let mut best := T.get! i
        for fi in b.topo.vertexFaces[i]! do
          match hopfLaxCandidate b region T faces[fi]! i with
          | some cand => if cand < best then best := cand
          | none => pure ()
        if best < T.get! i then
          T := T.set! i best
          changed := true
  return redistFinish nv phiP region adj T sweeps

/-- Quantities of the default stages (format 2). -/
structure Stage2Out where
  /-- length force data at the cut vertices -/
  force : CutForce
  /-- recovered length force per unit length `λ⁰` (zero for the length step on the mesh) -/
  lam0 : FloatArray
  /-- the rate `r` (after the length step on the cut, else `s`) -/
  rate : FloatArray
  /-- `ḡ` of the piece of every cut vertex and the piece label -/
  gbar : FloatArray
  /-- piece label (smallest cut-vertex index of the piece) -/
  piece : Array Nat
  /-- `ḡ · r` -/
  rg : FloatArray
  /-- rings of the local extension -/
  ring : Array (Option Nat)
  /-- the extension actually used -/
  sigma : FloatArray
  /-- roughening rate `ω` of every cut vertex -/
  omegaP : FloatArray
  /-- `ω` (largest over the cut vertices) -/
  omega : Float
  /-- the three limits of the step rule: `τ_max`, `time_scale/ω`, `max_motion h/max|ρ|` -/
  limits : Float × Float × Float
  /-- `max |λ⁰|` -/
  lamMax : Float
  /-- largest local roughening rate -/
  omegaMax : Float
  /-- local crossed-edge length `h_c` -/
  hc : FloatArray
  /-- `|s|/h_c` -/
  mc1 : FloatArray
  /-- `|s + λ⁰|/h_c` -/
  mc2 : FloatArray
  /-- the two statistics `m̂` (for `s`, for `s + λ⁰`; the extremes put the one used in both) -/
  mhat : Float × Float
  /-- the `m̂` that the step rule used -/
  mused : Float
  /-- factor of the local limit -/
  factor : FloatArray
  /-- rate after the local limit -/
  rlim : FloatArray
  deriving Inhabited

/-- Everything of one step. -/
structure StepResult where
  /-- scale `s` and total area of the caller's mesh -/
  base : Base
  /-- the scaled input -/
  phi : FloatArray
  /-- evaluation at the input -/
  ev : Eval
  /-- step size used -/
  tau : Float
  /-- shortest crossed edge -/
  hMin : Float
  /-- largest `|s|` after smoothing -/
  sMax : Float
  /-- `φ̃` -/
  phiTilde : FloatArray
  /-- `φ′` after the length / penalty step, before the symmetry projection -/
  phiLen : FloatArray
  /-- `φ′` after the symmetry projection -/
  phiSym : FloatArray
  /-- redistancing of `phiSym` -/
  redist : Redist
  /-- new level set function (after the second symmetry projection), scaled unit -/
  phiNew : FloatArray
  /-- evaluation at `phiNew` -/
  evNew : Eval
  /-- default-stage quantities (format 2) -/
  stage2 : Option Stage2Out := none
  /-- the local limit acted at the vertices (robust statistics, length step on the mesh, D26):
  `h_i`, `ω_i`, factor and limited value of every vertex -/
  vlimit : Option (FloatArray × FloatArray × FloatArray × FloatArray) := none

/-- Outcome: either a step, or the cut is absent (one region is empty). -/
inductive Outcome where
  /-- a full step -/
  | step (r : StepResult)
  /-- no crossing: only the evaluation at the input -/
  | noCut (base : Base) (phi : FloatArray) (ev : Eval)

/-- Validate the parameters of the case that do not concern the mesh. -/
def validateParams (c : CutCase) : Except String Unit := do
  let nv := c.nv
  if c.phi.size != nv then throw s!"phi has {c.phi.size} values for {nv} vertices"
  for i in [0:nv] do
    if !(c.phi.get! i).isFinite then throw s!"phi[{i}] is not finite"
  let pos := fun (x : Float) => x.isFinite && x ≥ 0
  if !pos c.alphaL then throw "alpha_L must be finite and >= 0"
  if let some e := c.eta then
    if e.size != nv then throw s!"eta has {e.size} values for {nv} vertices"
    for i in [0:nv] do
      if !pos (e.get! i) then throw s!"eta[{i}] must be finite and >= 0"
  if !(pos c.w.1 && pos c.w.2) then throw "w_0 and w_1 must be finite and >= 0"
  if !(c.beta.1.isFinite && c.beta.2.isFinite) then throw "beta_0 and beta_1 must be finite"
  if !(c.delta ≥ 0 && c.delta < 0.5) then throw "delta must be in [0, 1/2)"
  if !(c.rho0.isFinite && c.rho0 > 0) then throw "rho0 must be finite and positive"
  match c.step with
  | .fixed tau => if !(tau.isFinite && tau > 0) then throw "tau must be finite and positive"
  | .capped tm mm =>
    if !(tm.isFinite && tm > 0 && mm.isFinite && mm > 0) then
      throw "tau_max and max_motion must be finite and positive"
  | .adaptive ts mm tm =>
    if !(ts.isFinite && ts > 0 && tm.isFinite && tm > 0 && mm.isFinite && mm > 0) then
      throw "time_scale, max_motion and tau_max must be finite and positive"
  match c.flux with
  | .nodal => pure ()
  | .smoothed k => if !(k.isFinite && k ≥ 0) then throw "flux: c must be finite and >= 0"
  if let some cl := c.classes then
    if cl.size != nv then throw s!"classes has {cl.size} values for {nv} vertices"

/-- The length and penalty step on the mesh (SPEC §3.10, the paper's). -/
def lengthStepMesh (c : CutCase) (b : Base) (wv phiTilde : FloatArray) (tau : Float) :
    Except String FloatArray := do
  let nv := b.nv
  let allZero := Id.run do
    let mut z := true
    for i in [0:nv] do
      if wv.get! i != 0 then z := false
    return z
  if allZero then return phiTilde
  let mut mat := zeros (nv * nv)
  let mut rhs := zeros nv
  match c.penalty with
  | .gradient =>
    let mut lw := Mat.zero nv
    for fi in [0:c.faces.size] do
      let f := c.faces[fi]!
      let (c0, c1, c2) := b.cots[fi]!
      lw := lw.addTri f.a f.b f.c c0 c1 c2
        (triangleMean (wv.get! f.a) (wv.get! f.b) (wv.get! f.c))
    for i in [0:nv] do
      for j in [0:nv] do
        mat := mat.set! (i * nv + j) (lengthStepEntry (i == j) tau (b.vertArea.get! i) (lw.get i j))
      rhs := rhs.set! i (b.vertArea.get! i * phiTilde.get! i)
  | .paperDiagonal =>
    for i in [0:nv] do
      if !(wv.get! i > 0) then
        throw s!"the paper's penalty rule needs alpha_L + eta_i > 0; vertex {i} has {wv.get! i}"
    for i in [0:nv] do
      for j in [0:nv] do
        mat := mat.set! (i * nv + j)
          (penaltyStepEntry (i == j) tau (b.vertArea.get! i) (wv.get! i) (b.lap.get i j))
      rhs := rhs.set! i (b.vertArea.get! i / wv.get! i * phiTilde.get! i)
  match luFactor nv mat with
  | none => throw "length step: singular system"
  | some lu => return lu.solve rhs

/-- One flow step of SPEC §4. -/
def runStep (c : CutCase) : Except String Outcome := do
  validateParams c
  let b ← buildBase c
  let nv := b.nv
  let phi := FloatArray.ofFn' nv fun i => b.scale * c.phi.get! i
  let ev ← evaluate c b phi
  if ev.cm.nc == 0 then return .noCut b phi ev
  let cm := ev.cm
  let nc := cm.nc
  -- shortest crossed edge
  let mut hMin := posInf
  for cr in cm.crossings do
    if cr.len < hMin then hMin := cr.len
  let mut stage2 : Option Stage2Out := none
  let mut sMax := fzero
  let mut tau := fzero
  let mut sigma := ev.sigma
  if c.format == 1 then
    for k in [0:nc] do
      let a := (ev.sSmooth.get! k).abs
      if a > sMax then sMax := a
    tau := match c.step with
      | .fixed t => t
      | .capped tm mm => if sMax == 0 then tm else stepCap tm mm hMin sMax
      | .adaptive _ _ tm => tm
  else
    -- default stages / mixed variants
    let force := cutForce cm ev.wv (c.penalty == .paperDiagonal)
    let cc := fluxC c.flux
    let allZero := Id.run do
      let mut z := true
      for i in [0:nv] do
        if ev.wv.get! i != 0 then z := false
      return z
    -- length force per unit length on the cut (used by the rate, or for the step rule only)
    let mut lam0 := zeros nc
    if (c.lengthStep == .onCut || c.stats == .robust) && !allZero then
      let mut edges : Array (Nat × Nat × Float) := #[]
      for sg in cm.segs do
        let f := b.topo.faces[sg.face]!
        let h2 := faceMeanSquareEdge (b.l f.a f.b) (b.l f.b f.c) (b.l f.c f.a)
        edges := edges.push (sg.p, sg.q, fluxCoupling cc h2 sg.len)
      lam0 ← chainSolve nc force.ellC edges force.theta
    let lamUsed := if c.lengthStep == .onCut then lam0 else zeros nc
    let rho := FloatArray.ofFn' nc fun k => ev.sSmooth.get! k + lamUsed.get! k
    let rhoMaxV := maxAbs rho
    let lamMax := maxAbs lamUsed
    -- roughening rates
    let r0 := ev.regs[0]!
    let r1 := ev.regs[1]!
    let omegaP := FloatArray.ofFn' nc fun k =>
      roughRate (roughQ (regionDensity c.w.1 (r0.G.get! (nv + k)) c.beta.1) c.beta.1 c.w.1 0)
                (roughQ (regionDensity c.w.2 (r1.G.get! (nv + k)) c.beta.2) c.beta.2 c.w.2 0)
    let mut omegaMax := negInf
    for k in [0:nc] do
      if omegaP.get! k > omegaMax then omegaMax := omegaP.get! k
    -- local crossed-edge lengths h_c
    let hc := FloatArray.ofFn' nc fun k => Id.run do
      let cr := cm.crossings[k]!
      let mut h := posInf
      for k2 in [0:nc] do
        let c2 := cm.crossings[k2]!
        if c2.a == cr.a || c2.a == cr.b || c2.b == cr.a || c2.b == cr.b then
          if c2.len < h then h := c2.len
      return h
    let robust := c.stats == .robust
    let omegaList := omegaP.toList
    let ellList := force.ellC.toList
    let omega := if robust then weightedQuantile 0.9 omegaList ellList else omegaMax
    let mc1 := FloatArray.ofFn' nc fun k => speedInEdges (ev.sSmooth.get! k) (hc.get! k)
    let mc2 := FloatArray.ofFn' nc fun k => speedInEdges (ev.sSmooth.get! k + lam0.get! k) (hc.get! k)
    let q1 := weightedQuantile 0.9 mc1.toList ellList
    let q2 := weightedQuantile 0.9 mc2.toList ellList
    -- m-hat of the step rule
    let mhat : Float :=
      if robust then
        (if c.lengthStep == .onCut then q2 else max q1 q2)
      else if rhoMaxV > 0 then rhoMaxV / hMin else 0
    sMax := rhoMaxV
    let limits : Float × Float × Float := match c.step with
      | .fixed t => (t, posInf, posInf)
      | .capped tm mm =>
        (tm, posInf, if robust then (if mhat > 0 then mm / mhat else posInf)
                     else if rhoMaxV > 0 then mm * hMin / rhoMaxV else posInf)
      | .adaptive ts mm tm =>
        (tm, if omega > 0 then ts / omega else posInf,
         if robust then (if mhat > 0 then mm / mhat else posInf)
         else if rhoMaxV > 0 then mm * hMin / rhoMaxV else posInf)
    tau := match c.step with
      | .fixed t => t
      | .capped tm mm =>
        if robust then robustStep tm posInf 0 mm mhat
        else if rhoMaxV == 0 then tm else stepCap tm mm hMin rhoMaxV
      | .adaptive ts mm tm =>
        if robust then robustStep tm ts omega mm mhat
        else adaptiveStep tm ts omega mm hMin rhoMaxV
    -- the rate
    let mut rate := ev.sSmooth
    if c.lengthStep == .onCut && !allZero then
      let mut edges : Array (Nat × Nat × Float) := #[]
      for sg in cm.segs do
        let wbar := segmentMeanWeight (force.wc.get! sg.p) (force.wc.get! sg.q)
        edges := edges.push (sg.p, sg.q, rateCoupling tau wbar sg.len)
      let rhs := FloatArray.ofFn' nc fun k => rateRhs (force.ellC.get! k) (ev.sSmooth.get! k) (lam0.get! k)
      rate ← chainSolve nc force.ellC edges rhs
    -- local limit (robust statistics, capped and adaptive steps)
    let mut factor := FloatArray.ofFn' nc fun _ => 1.0
    if robust && c.lengthStep == .onCut then
      match c.step with
      | .fixed _ => pure ()
      | .capped _ mm =>
        for k in [0:nc] do
          factor := factor.set! k (if (rate.get! k).abs / hc.get! k ≤ mhat then 1.0 else localLimitFactor mm (hc.get! k) tau (rate.get! k).abs)
      | .adaptive _ mm _ =>
        for k in [0:nc] do
          let v := localSpeed (rate.get! k) (omegaP.get! k) omega
          factor := factor.set! k (if v / hc.get! k ≤ mhat then 1.0 else localLimitFactor mm (hc.get! k) tau v)
    let rlim := FloatArray.ofFn' nc fun k => rate.get! k * factor.get! k
    let (gbar, piece) :=
      if c.speed == .meanGradient then meanGradients cm phi
      else (FloatArray.ofFn' nc fun _ => 1.0, cutPieces cm)
    let rg := FloatArray.ofFn' nc fun k => gbar.get! k * rlim.get! k
    let mut ring : Array (Option Nat) := Array.replicate nv none
    match c.extension with
    | .local =>
      let (sg, rn) := localExtension b cm force.ellC rg
      sigma := sg
      ring := rn
    | .harmonic => sigma ← harmonicExtension cm ev.regs rg
    stage2 := some { force, lam0 := lamUsed, rate, gbar, piece, rg, ring, sigma, omegaP, omega, limits
                     lamMax, omegaMax, hc, mc1, mc2, mhat := (q1, q2), mused := mhat, factor, rlim }
  -- explicit step
  let phiTilde := FloatArray.ofFn' nv fun i => explicitStep (phi.get! i) tau (sigma.get! i)
  -- length and penalty step
  let phiLen ←
    if c.format ≥ 2 && c.lengthStep == .onCut then pure phiTilde
    else lengthStepMesh c b ev.wv phiTilde tau
  -- D26: with the length step on the mesh the local limit acts at the vertices on the whole step
  let mut vlimit : Option (FloatArray × FloatArray × FloatArray × FloatArray) := none
  let mut phiLim := phiLen
  if c.format ≥ 2 && c.stats == .robust && c.lengthStep == .onMesh then
    match c.step, stage2 with
    | .fixed _, _ => pure ()
    | _, none => pure ()
    | step, some st =>
      let mm := match step with | .capped _ m => m | .adaptive _ m _ => m | .fixed _ => 1.0
      let mut hv := zeros nv
      let mut omv := zeros nv
      let mut fv := zeros nv
      for i in [0:nv] do
        let mut h := posInf
        let mut om := st.omega
        let mut anyCross := false
        for k in [0:nc] do
          let cr := cm.crossings[k]!
          if cr.a == i || cr.b == i then
            anyCross := true
            if cr.len < h then h := cr.len
            if st.omegaP.get! k > om || om == st.omega then
              if st.omegaP.get! k > om then om := st.omegaP.get! k
        if !anyCross then
          for j in b.topo.neighbors[i]! do
            if b.l i j < h then h := b.l i j
        let delta := phiLen.get! i - phi.get! i
        let v := match step with
          | .capped _ _ => vertexSpeedCapped delta tau
          | _ => vertexSpeedAdaptive delta tau om st.omega
        let f := vertexLimitFactor mm h tau v
        hv := hv.set! i h
        omv := omv.set! i om
        fv := fv.set! i f
        phiLim := phiLim.set! i (limitedValue (phi.get! i) (phiLen.get! i) f)
      vlimit := some (hv, omv, fv, phiLim)
  let phiSym := symmetrize c.classes phiLim
  let redist := redistance b c.faces phiSym
  let phiNew := symmetrize c.classes redist.phi
  let evNew ← evaluate c b phiNew
  return .step { base := b, phi, ev, tau, hMin, sMax, phiTilde, phiLen, phiSym, redist, phiNew, evNew, stage2, vlimit }

end Ddg.Cuts
