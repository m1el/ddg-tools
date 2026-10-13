import DdgOracle.Selftest
import DdgOracle.Cuts.Gen

/-!
# `cuts-selftest`

Checks of the cut-flow oracle that need no external reference: identities, closed forms, analytic
solutions at loose tolerances, invariances, the redistancing fixed point, the text format.
-/

namespace Ddg.Cuts

/-- A step that must succeed. -/
def runOk (c : CutCase) : Option StepResult :=
  match runStep c with
  | .ok (.step r) => some r
  | _ => none

/-- Euler characteristic of the submesh of region `r`. -/
def regionChi (cm : CutMesh) (r : Nat) : Int := Id.run do
  let tris := cm.tris.filter (·.region == r)
  let mut nodes : Std.HashMap Nat Unit := {}
  let mut edges : Std.HashMap (Nat × Nat) Unit := {}
  for t in tris do
    for n in [t.n0, t.n1, t.n2] do nodes := nodes.insert n ()
    for (x, y) in [(t.n0, t.n1), (t.n1, t.n2), (t.n2, t.n0)] do
      edges := edges.insert (min x y, max x y) ()
  return (nodes.size : Int) - (edges.size : Int) + (tris.size : Int)

/-- Largest `|x|`. -/
def maxAbsF (a : FloatArray) : Float := maxAbs a

/-- Tilted-plane level set on a mesh. -/
def tilted (m : RawMesh) (nx ny nz d : Float) : FloatArray :=
  phiOf m fun p => nx * p.x + ny * p.y + nz * p.z - d

private def pointRadius (p : Point Float) : Float := Float.sqrt (p.x * p.x + p.y * p.y)

/-- Verify a fast-marching result against its own rule: the value of every non-seed reached
vertex is the minimum of the updates that use only vertices accepted before it (value `<=` its
own; ties are resolved by the index, so equal values are allowed). Returns the largest
violation. -/
def redistResidual (b : Base) (faces : Array Face) (rd : Redist) : Float := Id.run do
  let mut worst := fzero
  for i in [0:b.nv] do
    if !rd.adj[i]! && rd.T.get! i < posInf then
      let ti := rd.T.get! i
      let masked := FloatArray.ofFn' b.nv fun j =>
        if rd.T.get! j < ti || (rd.T.get! j == ti && j < i) then rd.T.get! j else posInf
      let mut best := posInf
      for fi in b.topo.vertexFaces[i]! do
        match hopfLaxCandidate b rd.region masked faces[fi]! i with
        | some cand => if cand < best then best := cand
        | none => pure ()
      let d := (best - ti).abs
      if d > worst then worst := d
  return worst

/-- Run all checks; returns the number of failures. -/
def selftest : IO Nat := do
  let tally ← IO.mkRef ({} : Tally)
  let pi := floatPi

  let ico3 := icosphere 3
  let psi0 := 1.0
  let capPhi := phiOf ico3 fun p => Float.acos p.z - psi0
  -- 1. meshes: icosphere, orientation, Euler characteristic
  for lvl in [0, 1, 2, 3] do
    let m := icosphere lvl
    let nv := m.pos.size
    let mut inward := 0
    for f in m.faces do
      let p := m.pos[f.a]!
      let q := m.pos[f.b]!
      let r := m.pos[f.c]!
      let cx := (q.y - p.y) * (r.z - p.z) - (q.z - p.z) * (r.y - p.y)
      let cy := (q.z - p.z) * (r.x - p.x) - (q.x - p.x) * (r.z - p.z)
      let cz := (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x)
      if cx * (p.x + q.x + r.x) + cy * (p.y + q.y + r.y) + cz * (p.z + q.z + r.z) ≤ 0 then
        inward := inward + 1
    let chi : Int := match buildTopology nv m.faces with | .ok t => t.chi | .error _ => 0
    check tally s!"icosphere level {lvl}: {nv} vertices, {m.faces.size} faces, all outward, chi = 2"
      (inward == 0 && chi == 2 && nv == 10 * 4 ^ lvl + 2)

  -- 2. sum of face curvatures = sum of interior angle defects (and = 2 pi chi on closed meshes)
  let meshes : List (String × RawMesh) :=
    [("icosphere 2", icosphere 2),
     ("torus", (torus 9 6 1.0 0.42 fun a b => 0.22 * Float.sin (2 * a) * Float.cos b).jitter 0.01 106),
     ("bumpy grid disk", (gridDisk 9 8 fun x y => 0.35 * Float.sin (4 * x) * Float.sin (3 * y)).jitter 0.01 101),
     ("annulus", (polar 4 12 0.4 1.0 fun x y => 0.3 * Float.sin (3 * x) + 0.25 * y * y).jitter 0.01 103),
     ("lumpy sphere", (uvSphere 8 10 fun p =>
        let r := 1 + 0.45 * Float.sin (3 * p.x) * Float.sin (3 * p.y)
        ⟨r * p.x, r * p.y, r * p.z⟩).jitter 0.004 105)]
  for (nm, m) in meshes do
    match buildBase (mkCase nm m (phiOf m fun p => p.x)) with
    | .error e => check tally s!"sum of face curvatures on {nm}" false e
    | .ok b =>
      let mut sumO := fzero
      for o in b.omega do sumO := sumO + o
      let mut sumD := fzero
      for v in [0:b.nv] do
        if !b.topo.isBoundary[v]! then sumD := sumD + (2 * pi - b.theta.get! v)
      let expected := if b.topo.numBoundary == 0 then 2 * pi * Float.ofInt b.topo.chi else sumD
      check tally s!"{nm}: sum_f Omega_f = sum over interior vertices of (2 pi - Theta)"
        ((sumO - sumD).abs ≤ 1e-12 && (sumO - expected).abs ≤ 1e-12)
        s!"sum={sci sumO} defects={sci sumD} diff={sci (sumO - sumD)}"

  -- 3. cut mesh invariants on a sphere cut by a tilted plane and on a torus
  let ico2 := icosphere 2
  let tCases : List (String × RawMesh × FloatArray) :=
    [("icosphere 2 / tilted plane", ico2, tilted ico2 0.3 0.2 0.93 0.15),
     ("torus / band", (torus 12 8 1.0 0.42 fun _ _ => 0.0),
       phiOf (torus 12 8 1.0 0.42 fun _ _ => 0.0) fun p => Float.cos (3 * Float.atan2 p.y p.x) - 0.2),
     ("grid disk / plane", gridDisk 10 9 (fun x y => 0.3 * x * y), tilted (gridDisk 10 9 (fun x y => 0.3 * x * y)) 1 0.3 0 0.47)]
  for (nm, m, phi) in tCases do
    match runOk (mkCase nm m phi) with
    | none => check tally s!"{nm}: step runs" false
    | some r =>
      let cm := r.ev.cm
      let mut area := fzero
      for t in cm.tris do area := area + t.area
      let chiM : Int := (match buildTopology m.pos.size m.faces with | .ok t => t.chi | .error _ => 0)
      let chiG : Int := (cm.nc : Int) - (cm.segs.size : Int)
      let lhs := regionChi cm 0 + regionChi cm 1 - chiG
      check tally s!"{nm}: cut mesh area = 1, chi(S0)+chi(S1)-chi(cut) = chi(M)"
        ((area - 1).abs ≤ 1e-12 && lhs == chiM)
        s!"area-1={sci (area - 1)} chi {lhs} vs {chiM}, {cm.nc} cut vertices, {cm.tris.size} triangles"

  -- 4. flat disk: u = 0, s = 0
  let flat := polar 8 24 0.0 1.0 fun _ _ => 0.0
  let flatPhi := phiOf flat fun p => pointRadius p - 0.55
  match runOk (mkCase "flat" flat flatPhi) with
  | none => check tally "flat disk: step runs" false
  | some r =>
    let mut uMax := fzero
    let mut vMax := fzero
    for rd in r.ev.regs do
      uMax := max uMax (maxAbs rd.u)
      vMax := max vMax (maxAbs rd.v)
    let sMax := maxAbs r.ev.sRaw
    check tally "flat disk: u = 0, v = 0, s = 0 (to rounding)"
      (uMax ≤ 1e-13 && vMax ≤ 1e-13 && sMax ≤ 1e-12)
      s!"max|u|={sci uMax} max|v|={sci vMax} max|s|={sci sMax}"

  -- 5. scaled Dirichlet with a = 1 equals Dirichlet; Hencky / Dirichlet differ
  let bumpy := (gridDisk 11 10 fun x y => 0.4 * Float.sin (4 * x) * Float.sin (3 * y)).jitter 0.004 7
  let bphi := tilted bumpy 1 0.4 0 0.55
  let cD := mkCase "dir" bumpy bphi
  let cS := { cD with energy := .scaledA (FloatArray.ofFn' bumpy.pos.size fun _ => 1.0) }
  match runOk cD, runOk cS with
  | some a, some b =>
    let mut dJ := (a.ev.J0 - b.ev.J0).abs + (a.ev.J1 - b.ev.J1).abs
    let mut dG := fzero
    for k in [0:a.ev.cm.nc] do
      dG := max dG ((a.ev.sRaw.get! k - b.ev.sRaw.get! k).abs)
    let mut dP := fzero
    for rr in [0:2] do
      let ra := a.ev.regs[rr]!
      let rb := b.ev.regs[rr]!
      for n in [0:a.ev.cm.nodes] do
        dP := max dP ((rb.p.get! n + 2 * ra.u.get! n).abs)
        dP := max dP ((rb.q.get! n + 2 * ra.v.get! n).abs)
    dJ := dJ + dG
    check tally "scaled Dirichlet with a = 1: J, s equal Dirichlet; p = -2u, q = -2v"
      (dJ ≤ 1e-11 && dP ≤ 1e-11) s!"J,s diff={sci dJ} p,q diff={sci dP}"
  | _, _ => check tally "scaled a = 1: steps run" false

  -- 6. flat disk with a constant source K: closed forms of SPEC §1
  let disk := polar 16 48 0.0 1.0 fun _ _ => 0.0
  let area0 := pi  -- area of the disk of radius 1
  let sc := 1 / Float.sqrt area0
  let rho := 0.3          -- scaled units
  let kconst := 5.0
  let diskPhi := phiOf disk fun p => pointRadius p - rho / sc
  let cF := { mkCase "flatK" disk diskPhi with curv := .uniform kconst }
  let jDir := pi * kconst * kconst * rho ^ 4 / 8
  let gDir := (kconst * rho / 2) * (kconst * rho / 2)
  let jHen := pi * kconst * kconst * rho ^ 6 / 48
  let gHen := kconst * kconst * rho ^ 4 / 16
  for (nm, en) in [("Dirichlet", Energy.dirichlet), ("Hencky", Energy.hencky)] do
    match runOk { cF with energy := en } with
    | none => check tally s!"flat disk + constant source, {nm}: step runs" false
    | some r =>
      let rd := r.ev.regs[0]!
      let nv := r.ev.cm.nv
      let nc := r.ev.cm.nc
      let mut gSum := fzero
      let mut lSum := fzero
      for k in [0:nc] do
        gSum := gSum + rd.G.get! (nv + k) * rd.ell.get! (nv + k)
        lSum := lSum + rd.ell.get! (nv + k)
      let gMean := gSum / lSum
      let jExp := if nm == "Dirichlet" then jDir else jHen
      let gExp := if nm == "Dirichlet" then gDir else gHen
      -- u against K (r^2 - rho^2)/4 at the vertices of the disk
      let mut uErr := fzero
      for i in [0:nv] do
        if r.ev.cm.region[i]! == 0 then
          let rr := pointRadius disk.pos[i]! * sc
          uErr := max uErr ((rd.u.get! i - kconst * (rr * rr - rho * rho) / 4).abs)
      let relJ := (rd.J - jExp).abs / jExp
      let relG := (gMean - gExp).abs / gExp
      check tally s!"flat disk + constant source K = 5, {nm}: J and mean G against the closed form (5%, 10%)"
        (relJ ≤ 0.05 && relG ≤ 0.10 && uErr ≤ 0.01)
        s!"J={sci rd.J} exact {sci jExp} (rel {sci relJ}); mean G={sci gMean} exact {sci gExp} (rel {sci relG}); max|u-exact|={sci uErr}"

  -- 7. sphere cap (SPEC §6 item 2) on the icosphere of level 3 (642 vertices)
  match runOk (mkCase "cap" ico3 capPhi) with
  | none => check tally "sphere cap: step runs" false
  | some r =>
    let rd := r.ev.regs[0]!
    let nv := r.ev.cm.nv
    let nc := r.ev.cm.nc
    let Rs := 1 / Float.sqrt (4 * pi)
    let jExp := -8 * pi * Float.log (Float.cos (psi0 / 2)) - 4 * pi * (Float.sin (psi0 / 2)) ^ 2
    let dnExp := Float.tan (psi0 / 2) / Rs
    let mut vSum := fzero
    let mut lSum := fzero
    for k in [0:nc] do
      vSum := vSum + rd.v.get! (nv + k)
      lSum := lSum + rd.ell.get! (nv + k)
    let mut gSumAll := fzero
    for n in [0:r.ev.cm.nodes] do
      if rd.isB[n]! then gSumAll := gSumAll + rd.v.get! n
    let mut sumg := fzero
    for t in r.ev.cm.tris do
      if t.region == 0 then sumg := sumg + r.base.kUsed[t.parent]! * t.area
    let mut uErr := fzero
    for i in [0:nv] do
      if r.ev.cm.region[i]! == 0 then
        let psi := Float.acos ico3.pos[i]!.z
        let ue := 2 * Float.log (Float.cos (psi0 / 2) / Float.cos (psi / 2))
        uErr := max uErr ((rd.u.get! i - ue).abs)
    let relJ := (rd.J - jExp).abs / jExp
    let dn := vSum / lSum
    let relD := (dn - dnExp).abs / dnExp
    check tally "sphere cap, icosphere 3: sum of boundary fluxes v = total curvature of the cap (rounding)"
      ((gSumAll - sumg).abs ≤ 1e-11) s!"sum v={sci gSumAll} sum k A={sci sumg}"
    check tally "sphere cap, icosphere 3: J, du/dn and u against SPEC section 6 item 2 (3%, 5%, 0.03)"
      (relJ ≤ 0.03 && relD ≤ 0.05 && uErr ≤ 0.03)
      s!"J={sci rd.J} exact {sci jExp} (rel {sci relJ}); du/dn={sci dn} exact {sci dnExp} (rel {sci relD}); max|u-exact|={sci uErr}"

  -- 8. redistancing output satisfies its own fixed-point equation; Adj keep phi'
  for (nm, m, phi) in tCases do
    match runOk { mkCase nm m phi with step := .capped 0.01 0.2 } with
    | none => check tally s!"{nm}: step runs" false
    | some r =>
      let res := redistResidual r.base (mkCase nm m phi).faces r.redist
      let mut keep := fzero
      for i in [0:r.base.nv] do
        if r.redist.adj[i]! then keep := max keep ((r.phiSym.get! i - r.redist.phi.get! i).abs)
      check tally s!"{nm}: fast marching satisfies its own acceptance rule; cut-adjacent values kept"
        (res ≤ 1e-13 && keep == 0 && r.redist.unreached == 0)
        s!"residual={sci res} sweeps={r.redist.sweeps} unreached={r.redist.unreached} tau={sci r.tau} max|s|={sci r.sMax} nc={r.ev.cm.nc} max|phi|={sci (maxAbs r.phi)} new cut vertices={r.evNew.cm.nc}"

  -- 9. closed form of the local Hopf–Lax update against a brute-force scan
  let mut rng : Rng := ⟨12345⟩
  let mut worstBelow := fzero
  let mut worstAbove := fzero
  for _ in [0:300] do
    let (xj, r1) := rng.range (-2.0) 2.0
    let (yj, r2) := r1.range 0.1 2.0
    let (xk, r3) := r2.range (-2.0) 2.0
    let (yk, r4) := r3.range (-2.0) 0.0
    let (tj, r5) := r4.range 0.0 3.0
    let (tk, r6) := r5.range 0.0 3.0
    rng := r6
    -- x_i = origin; x_j above, x_k below the x axis, so the triangle is not degenerate generically
    let lij := Float.sqrt (xj * xj + yj * yj)
    let lik := Float.sqrt (xk * xk + yk * yk)
    let ljk := Float.sqrt ((xj - xk) * (xj - xk) + (yj - yk) * (yj - yk))
    let closed := hopfLaxTriangle lij lik ljk tj tk
    let f := fun (lam : Float) =>
      let px := lam * xj + (1 - lam) * xk
      let py := lam * yj + (1 - lam) * yk
      lam * tj + (1 - lam) * tk + Float.sqrt (px * px + py * py)
    -- the function is convex: ternary search
    let mut lo := 0.0
    let mut hi := 1.0
    for _ in [0:200] do
      let m1 := lo + (hi - lo) / 3
      let m2 := hi - (hi - lo) / 3
      if f m1 < f m2 then hi := m2 else lo := m1
    let scan := min (f ((lo + hi) / 2)) (min (f 0) (f 1))
    worstBelow := max worstBelow (closed - scan)   -- closed form may not exceed the scan
    worstAbove := max worstAbove (scan - closed)   -- the scan is coarse
  check tally "Hopf-Lax closed form against a ternary search over lambda (300 random triangles)"
    (worstBelow ≤ 1e-12 && worstAbove ≤ 1e-12)
    s!"closed - search <= {sci worstBelow}; search - closed <= {sci worstAbove}"

  -- 10. invariance under a change of unit: lengths and phi times 3
  let c1 := { mkCase "inv" ico2 (tilted ico2 0.3 0.2 0.93 0.15) with energy := .hencky, step := .capped 0.01 0.2 }
  let c3 := { c1 with edges := c1.edges.map (fun (a, b, l) => (a, b, 3 * l)), phi := FloatArray.ofFn' c1.phi.size fun i => 3 * c1.phi.get! i }
  match runOk c1, runOk c3 with
  | some a, some b =>
    let dE := (a.ev.total - b.ev.total).abs
    let mut dPhi := fzero
    let mut dSig := fzero
    for i in [0:a.base.nv] do
      dPhi := max dPhi ((a.phiNew.get! i - b.phiNew.get! i).abs)
      dSig := max dSig ((a.ev.sigma.get! i - b.ev.sigma.get! i).abs)
    let mut dS := fzero
    for k in [0:a.ev.cm.nc] do
      dS := max dS ((a.ev.sSmooth.get! k - b.ev.sSmooth.get! k).abs)
    check tally "uniform change of unit (lengths and phi times 3): scaled quantities unchanged"
      (dE ≤ 1e-11 * (1 + a.ev.total.abs) && dPhi ≤ 1e-11 && dSig ≤ 1e-11 * (1 + maxAbs a.ev.sigma)) s!"E diff={sci dE} sigma diff={sci dSig} s diff={sci dS} max|s|={sci a.sMax} phi_new diff={sci dPhi} nc {a.ev.cm.nc} {b.ev.cm.nc} sweeps {a.redist.sweeps} {b.redist.sweeps}"
  | _, _ => check tally "scale invariance: steps run" false

  -- 11. energy decreases for a small step (informational check on one mesh)
  let c11 := { mkCase "dec" ico2 (tilted ico2 0.3 0.2 0.93 0.15) with step := .capped 0.01 0.1 }
  match runOk c11 with
  | some r =>
    check tally "energy decreases along one small step on icosphere 2 (tilted plane)"
      (r.evNew.total < r.ev.total) s!"E {sci r.ev.total} -> {sci r.evNew.total}"
  | none => check tally "energy decrease: step runs" false

  -- 11b. finite differences (SPEC §6 item 3): moving the whole level set by a constant e changes
  -- the energy by -2 e * sum_c s_c l_c (no length term: alpha_L = 0)
  for (nm, en, w, beta) in [("Dirichlet", Energy.dirichlet, (1.0, 1.0), (0.0, 0.0)),
                            ("Hencky, w = (1, 2), beta = (0.5, 0.1)", Energy.hencky, (1.0, 2.0), (0.5, 0.1)),
                            ("curvature-scaled", Energy.scaledCurv 1.0 0.02, (1.5, 1.0), (0.0, 0.0))] do
    let cFD := { mkCase "fd" ico3 capPhi with alphaL := 0.0, energy := en, w, beta }
    match buildBase cFD with
    | .error e => check tally s!"finite differences, {nm}" false e
    | .ok b =>
      let phi0 := FloatArray.ofFn' b.nv fun i => b.scale * capPhi.get! i
      let e := 0.0005
      match evaluate cFD b phi0, evaluate cFD b (FloatArray.ofFn' b.nv fun i => phi0.get! i + e),
            evaluate cFD b (FloatArray.ofFn' b.nv fun i => phi0.get! i - e) with
      | .ok ev, .ok evP, .ok evM =>
        let mut dual := zeros ev.cm.nc
        for sg in ev.cm.segs do
          dual := dual.set! sg.p (dual.get! sg.p + sg.len / 2)
          dual := dual.set! sg.q (dual.get! sg.q + sg.len / 2)
        let mut pred := fzero
        for k in [0:ev.cm.nc] do pred := pred - 2 * e * ev.sRaw.get! k * dual.get! k
        let fd := evP.total - evM.total
        let rel := (fd - pred).abs / pred.abs
        check tally s!"finite differences of E under a shift of the level set, {nm} (icosphere 3, 10%)"
          (rel ≤ 0.10) s!"E(+e)-E(-e)={sci fd} predicted {sci pred} (rel {sci rel})"
      | _, _, _ => check tally s!"finite differences, {nm}: evaluations run" false

  -- 8b. redistancing on a flat grid against the exact distance to a straight cut
  let gridF := gridDisk 21 21 (fun _ _ => 0.0)
  let phiP := phiOf gridF fun p => p.x - 0.4871
  match buildBase (mkCase "rd" gridF phiP) with
  | .error e => check tally "redistancing on a flat grid" false e
  | .ok b =>
    let faces := (mkCase "rd" gridF phiP).faces
    let rd := redistance b faces phiP
    let rf := redistanceFixedPoint b faces phiP
    let mut dFix := fzero
    for i in [0:b.nv] do dFix := max dFix ((rd.phi.get! i - rf.phi.get! i).abs)
    check tally "single-pass fast marching equals the fixed point of the scheme on a mesh without obtuse angles"
      (dFix ≤ 1e-13) s!"max diff={sci dFix}"
    let mut err := fzero
    for i in [0:b.nv] do
      err := max err ((rd.phi.get! i - phiP.get! i).abs)
    check tally "redistancing of a flat grid with a straight cut reproduces the exact signed distance"
      (err ≤ 1e-12) s!"max error={sci err} sweeps={rd.sweeps}"

  -- 11c. equations of the stages, checked on the output: Yamabe system, extension, (re)assembled
  -- from the cut-mesh triangles
  match runOk { mkCase "eq" ico2 (phiOf ico2 fun p => Float.acos p.z - 1.0) with energy := .hencky, step := .fixed 0.001 } with
  | none => check tally "stage equations: step runs" false
  | some r =>
    let cm := r.ev.cm
    let mut yam := fzero
    let mut ext := fzero
    let mut adj := fzero
    for rd in r.ev.regs do
      let lu := rd.lap.mulVec rd.u
      let mut g := zeros cm.nodes
      for t in cm.tris do
        if t.region == rd.r then
          for n in [t.n0, t.n1, t.n2] do
            g := g.set! n (g.get! n + r.base.kUsed[t.parent]! * t.area / 3)
      let mut sg := zeros cm.nodes
      for i in [0:cm.nv] do sg := sg.set! i (r.ev.sigma.get! i)
      for k in [0:cm.nc] do sg := sg.set! (cm.nv + k) (r.ev.sSmooth.get! k)
      let ls := rd.lap.mulVec sg
      let lp := rd.lap.mulVec rd.p
      for n in [0:cm.nodes] do
        if rd.solved[n]! && !rd.isB[n]! then
          yam := max yam ((lu.get! n + g.get! n).abs)
          -- Hencky adjoint: L p = -A_S b'(u) = -2 A_S u
          let mut aS := fzero
          for t in cm.tris do
            if t.region == rd.r && (t.n0 == n || t.n1 == n || t.n2 == n) then aS := aS + t.area / 3
          adj := max adj ((lp.get! n + 2 * aS * rd.u.get! n).abs)
        if rd.inS[n]! && n < cm.nv then
          ext := max ext (ls.get! n).abs
    check tally "stage equations: (L u)_I = -g_I, L_II p = -A b'(u), harmonic extension (L sigma)_O = 0 (to rounding)"
      (yam ≤ 1e-11 && adj ≤ 1e-10 && ext ≤ 1e-10) s!"yamabe={sci yam} adjoint={sci adj} extension={sci ext}"

  -- 11d. flux recovery (SPEC §3.5): conservation, maximum principle, c = 0 is nodal, blow-up cases
  let blow : List (String × RawMesh × FloatArray) :=
    [("icosphere 2, tilted plane (nodal max|s| was 4.7e6)", ico2, tilted ico2 0.3 0.2 0.93 0.15),
     ("icosphere 1, phi = z, zeros at vertices (nodal max|s| was 1.1e6)", icosphere 1, phiOf (icosphere 1) fun p => p.z),
     ("icosphere 2, geodesic cap", ico2, phiOf ico2 fun p => Float.acos p.z - 1.0)]
  for (nm, m, phi) in blow do
    let cS := { mkCase "rec" m phi with energy := .hencky }
    let cN := { cS with flux := .nodal }
    let c0 := { cS with flux := .smoothed 0.0 }
    match runOk cS, runOk cN, runOk c0 with
    | some a, some bN, some z =>
      let cm := a.ev.cm
      let nv := cm.nv
      -- pieces of the cut
      let mut par : Array Nat := Array.range cm.nc
      let rec find (par : Array Nat) (x : Nat) (fuel : Nat) : Nat :=
        match fuel with
        | 0 => x
        | f + 1 => if par[x]! == x then x else find par par[x]! f
      for sg in cm.segs do
        let rp := find par sg.p cm.nc
        let rq := find par sg.q cm.nc
        par := par.set! rp rq
      let mut worstCons := fzero
      let mut worstRange := fzero
      for rd in a.ev.regs do
        for root in [0:cm.nc] do
          if find par root cm.nc == root then
            let mut sv := fzero
            let mut sld := fzero
            let mut lo := posInf
            let mut hi := negInf
            let mut scale := fzero
            for k in [0:cm.nc] do
              if find par k cm.nc == root then
                let n := nv + k
                let q := rd.v.get! n / rd.ell.get! n
                sv := sv + rd.v.get! n
                sld := sld + rd.ell.get! n * rd.d.get! n
                scale := scale + (rd.v.get! n).abs
                lo := min lo q
                hi := max hi q
            worstCons := max worstCons ((sv - sld).abs / (scale + 1e-300))
            for k in [0:cm.nc] do
              if find par k cm.nc == root then
                let dk := rd.d.get! (nv + k)
                let viol := max (lo - dk) (dk - hi)
                worstRange := max worstRange (viol / (lo.abs + hi.abs + 1e-300))
      let sMaxS := maxAbs a.ev.sRaw
      let sMaxN := maxAbs bN.ev.sRaw
      let mut dZero := fzero
      for k in [0:cm.nc] do
        dZero := max dZero ((z.ev.sRaw.get! k - bN.ev.sRaw.get! k).abs / (1 + sMaxN))
      check tally s!"flux recovery, {nm}: sum l d = sum v per piece, d in the range of v/l, c = 0 is nodal"
        (worstCons ≤ 1e-12 && worstRange ≤ 1e-12 && dZero ≤ 1e-12)
        s!"conservation {sci worstCons} range {sci worstRange} c=0 vs nodal {sci dZero}; max|s| nodal {sci sMaxN} -> smoothed {sci sMaxS}"
    | _, _, _ => check tally s!"flux recovery, {nm}: steps run" false

  -- 11e. delta = 0: allowed; a crossing on a vertex is then an error
  let d0 := { mkCase "d0" ico2 (phiOf ico2 fun p => Float.acos p.z - 1.0) with delta := 0.0 }
  let d0z := { mkCase "d0z" (icosphere 1) (phiOf (icosphere 1) fun p => p.z) with delta := 0.0 }
  check tally "delta = 0: runs when no crossing is on a vertex, no crossing is moved; error when one is"
    ((match runOk d0 with | some r => r.ev.cm.crossings.all (fun c => !c.moved) | none => false)
      && (runOk d0z).isNone && ((runCase d0z).splitOn "\n")[2]!.startsWith "status error")

  -- 11f. default stages (SPEC §3.8 - §3.10, §4)
  let dflt := fun (c : CutCase) =>
    { c with format := 2, lengthStep := .onCut, extension := .local, speed := .meanGradient,
             featPos := true, step := .adaptive 0.05 0.25 1.0 }
  -- flat disk: the length forces of a closed cut sum to 2 pi w; circle law for one small step
  let flatD := polar 12 36 0.0 1.0 fun _ _ => 0.0
  let sc2 := 1 / Float.sqrt floatPi
  let rho2 := 0.3
  let cF2 := dflt { mkCase "circ" flatD (phiOf flatD fun p => pointRadius p - rho2 / sc2) with
    alphaL := 1.0 }
  match runOk cF2 with
  | none => check tally "default stages, flat disk: step runs" false
  | some r =>
    match r.stage2 with
    | none => check tally "default stages: stage data present" false
    | some st =>
      let mut sumTheta := fzero
      for k in [0:r.ev.cm.nc] do sumTheta := sumTheta + st.force.theta.get! k
      let mut sumEll := fzero
      let mut sumLam := fzero
      let mut sumEllR := fzero
      let mut sumEllSL := fzero
      for k in [0:r.ev.cm.nc] do
        sumLam := sumLam + st.force.ellC.get! k * st.lam0.get! k
        sumEllR := sumEllR + st.force.ellC.get! k * st.rate.get! k
        sumEllSL := sumEllSL + st.force.ellC.get! k * (r.ev.sSmooth.get! k + st.lam0.get! k)
        sumEll := sumEll + st.force.ellC.get! k
      -- mean radius of the new cut
      let nvv := r.base.nv
      let mut rad := fzero
      let mut cnt := fzero
      for cr in r.evNew.cm.crossings do
        let pa := flatD.pos[cr.a]!
        let pb := flatD.pos[cr.b]!
        let x := (1 - cr.t) * pa.x + cr.t * pb.x
        let y := (1 - cr.t) * pa.y + cr.t * pb.y
        rad := rad + Float.sqrt (x * x + y * y) * sc2
        cnt := cnt + 1
      let radius := rad / cnt
      let mut rad0 := fzero
      for cr in r.ev.cm.crossings do
        let pa := flatD.pos[cr.a]!
        let pb := flatD.pos[cr.b]!
        let x := (1 - cr.t) * pa.x + cr.t * pb.x
        let y := (1 - cr.t) * pa.y + cr.t * pb.y
        rad0 := rad0 + Float.sqrt (x * x + y * y) * sc2
      let radius0 := rad0 / Float.ofNat r.ev.cm.nc
      let exact := radius0 - (rho2 - Float.sqrt (rho2 * rho2 - 2 * 1.0 * r.tau))
      check tally "default stages, flat disk: sum of length forces = 2 pi w; sums kept by both chain systems; circle law after one step"
        ((sumTheta - 2 * floatPi * 1.0).abs ≤ 1e-12 && (sumLam - sumTheta).abs ≤ 1e-12
          && (sumEllR - sumEllSL).abs ≤ 1e-12 && (radius - exact).abs / exact ≤ 0.01)
        s!"sum Theta={sci sumTheta} (2 pi={sci (2 * floatPi)}); sum l lam0 - sum Theta={sci (sumLam - sumTheta)}; sum l r - sum l (s+lam0)={sci (sumEllR - sumEllSL)}; radius {sci radius} (was {sci radius0}) vs exact {sci exact} (rel {sci ((radius - exact).abs / exact)}); nv={nvv}"
  -- chain systems: row sums, range, rest state
  let mut rng2 : Rng := ⟨777⟩
  let mut worstSum := fzero
  let mut worstRange := fzero
  let mut worstRest := fzero
  for _ in [0:40] do
    let n := 7
    let mut ell := zeros n
    let mut rhs := zeros n
    for k in [0:n] do
      let (x, r1) := rng2.range 0.1 2.0
      let (y, r2) := r1.range (-3.0) 3.0
      rng2 := r2
      ell := ell.set! k x
      rhs := rhs.set! k (x * y)
    let mut edges : Array (Nat × Nat × Float) := #[]
    for k in [0:n] do
      let (kap, r1) := rng2.range 0.0 5.0
      rng2 := r1
      edges := edges.push (k, (k + 1) % n, kap)
    match chainSolve n ell edges rhs, chainSolve n ell edges (zeros n) with
    | .ok x, .ok z =>
      let mut sl := fzero
      let mut sr := fzero
      let mut lo := posInf
      let mut hi := negInf
      for k in [0:n] do
        sl := sl + ell.get! k * x.get! k
        sr := sr + rhs.get! k
        lo := min lo (rhs.get! k / ell.get! k)
        hi := max hi (rhs.get! k / ell.get! k)
      worstSum := max worstSum ((sl - sr).abs)
      for k in [0:n] do
        worstRange := max worstRange (max (lo - x.get! k) (x.get! k - hi))
      worstRest := max worstRest (maxAbs z)
    | _, _ => worstSum := 1.0
  check tally "chain systems (flux recovery, lambda0, rate): sum l x = sum rhs, x in the range of rhs/l, zero data stays zero"
    (worstSum ≤ 1e-12 && worstRange ≤ 1e-12 && worstRest == 0)
    s!"sum {sci worstSum} range {sci worstRange} rest {sci worstRest}"
  -- local extension: range and mean rate
  for (nm, m, phi) in [("icosphere 2 geodesic cap", ico2, phiOf ico2 fun p => Float.acos p.z - 1.0),
                       ("torus band", (torus 12 8 1.0 0.42 fun _ _ => 0.0),
                        phiOf (torus 12 8 1.0 0.42 fun _ _ => 0.0) fun p => Float.cos (3 * Float.atan2 p.y p.x) - 0.2)] do
    match runOk (dflt { mkCase "ext" m phi with energy := .hencky, alphaL := 0.5 }) with
    | none => check tally s!"local extension, {nm}: step runs" false
    | some r =>
      let st := r.stage2.get!
      let cm := r.ev.cm
      let mut lo := posInf
      let mut hi := negInf
      for k in [0:cm.nc] do
        lo := min lo (st.rg.get! k)
        hi := max hi (st.rg.get! k)
      let mut viol := fzero
      for i in [0:cm.nv] do
        if st.ring[i]!.isSome then
          viol := max viol (max (lo - st.sigma.get! i) (st.sigma.get! i - hi))
      let mut lhs := fzero
      let mut rhs := fzero
      for k in [0:cm.nc] do
        let cr := cm.crossings[k]!
        lhs := lhs + st.force.ellC.get! k * ((1 - cr.t) * st.sigma.get! cr.a + cr.t * st.sigma.get! cr.b)
        rhs := rhs + st.force.ellC.get! k * st.rg.get! k
      check tally s!"extension from the cut, {nm}: sigma within the range of the rates; mean rate of the cut kept"
        (viol ≤ 1e-12 * (1 + hi.abs + lo.abs) && (lhs - rhs).abs ≤ 1e-11 * (1 + rhs.abs))
        s!"range violation {sci viol}; sum l interp(sigma) - sum l r = {sci (lhs - rhs)}"
  -- energy decrease for a small step, default stages; format 2 with the paper's stages equals format 1
  match runOk (dflt { mkCase "dec2" ico2 (phiOf ico2 fun p => Float.acos p.z - 1.0) with alphaL := 2.0, step := .fixed 0.0003 }) with
  | some r =>
    check tally "default stages: energy decreases along one small step (icosphere 2, cap, alpha_L = 2)"
      (r.evNew.total < r.ev.total) s!"E {sci r.ev.total} -> {sci r.evNew.total}"
  | none => check tally "default stages: energy decrease: step runs" false
  let cP1 := { mkCase "p1" ico2 (phiOf ico2 fun p => Float.acos p.z - 1.0) with step := .fixed 0.001 }
  let cP2 := { cP1 with format := 2 }
  match runOk cP1, runOk cP2 with
  | some a, some b =>
    let mut d := fzero
    for i in [0:a.base.nv] do d := max d ((a.phiNew.get! i - b.phiNew.get! i).abs)
    check tally "format 2 with the paper's stages gives the same new state as format 1" (d == 0) s!"max diff={sci d}"
  | _, _ => check tally "format 2 vs format 1: steps run" false

  -- 12. text format: print / parse round trip and error blocks
  let c12 := { mkCase "rt" ico2 (tilted ico2 0.3 0.2 0.93 0.15) with
    energy := .scaledCurv 0.5 0.02, eta := some (FloatArray.ofFn' ico2.pos.size fun i => 0.01 * Float.ofNat (i % 3)),
    classes := some (Array.ofFn (n := ico2.pos.size) fun i => i.val % 5), step := .capped 0.1 0.3 }
  let text := printCase c12
  let ok12 := match (Id.run do
      let lines := (text.splitOn "\n").filterMap fun l =>
        let t := tokens l
        if t.size == 0 || t[0]! == "case" || t[0]! == "end" then none else some t
      pure (parseBody "rt" lines.toArray)) with
    | .ok c => printCase c == text
    | .error _ => false
  check tally "text format: print, parse, print is the identity" ok12
  let bad := runCase { c12 with delta := 0.7 }
  check tally "invalid delta gives an error block" ((bad.splitOn "\n")[2]!.startsWith "status error")
  let bad2 := runCase { c12 with phi := zeros 3 }
  check tally "wrong phi length gives an error block" ((bad2.splitOn "\n")[2]!.startsWith "status error")

  let t ← tally.get
  IO.println s!"cuts-selftest: {t.passed} passed, {t.failed} failed"
  return t.failed

end Ddg.Cuts
