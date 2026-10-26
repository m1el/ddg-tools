import DdgOracle.RepulsiveCurves.Gen

/-!
# Internal consistency checks of the `rc` oracle

Each check compares two independent routes through the SPEC: central differences against the
analytic differential and Jacobians, the quadratic forms of §4.2-4.3 against the increments and the
kernel-matrix form of §4.4, the saddle equations against the solution, scaling laws, the symmetric
circle. Run on seeded random networks.
-/

namespace Ddg.RepulsiveCurves

private def relErr (a b scale : Float) : Float := (a - b).abs / (if scale > 1e-300 then scale else 1.0)

/-- Quadratic form `uᵀ M v` of a dense `n × n` matrix. -/
private def bilinear (n : Nat) (m : FloatArray) (u v : Array Float) : Float := Id.run do
  let mut s : Float := 0.0
  for i in [0:n] do
    for j in [0:n] do
      s := s + u[i]! * m.get! (i * n + j) * v[j]!
  return s

/-- `uᵀ B v` straight from the definition of §4.2: `Σ_{(I,J)} w_IJ ⟨D_I u − D_J u, D_I v − D_J v⟩` with
`(D u)_I = (u_{i₂} − u_{i₁}) T_I / ℓ_I`, and `uᵀ B⁰ v` of §4.3. -/
private def formsDirect (net : Net) (ex : Expo) (q : Quad) (u v : Array Float) : Float × Float := Id.run do
  let g := net.geo
  let mut sb : Float := 0.0
  let mut s0 : Float := 0.0
  let d := fun (w : Array Float) (I : Nat) =>
    let (a, b) := net.edges[I]!
    ((w[b]! - w[a]!) / g.len[I]!) * g.tan[I]!
  let e := fun (w : Array Float) (I : Nat) =>
    let (a, b) := net.edges[I]!
    0.5 * (w[a]! + w[b]!)
  for I in [0:net.ne] do
    for J in [0:net.ne] do
      if net.adj I J then continue
      let (w, w0) := pairWeights net g ex q I J
      sb := sb + w * (d u I - d u J).dot (d v I - d v J)
      s0 := s0 + w0 * (e u I - e u J) * (e v I - e v J)
  return (sb, s0)

/-- Quadratic forms through Eq. 31: `uᵀ Dᵀ((diag(K𝟙) − K) ⊗ I₃) D v` and `uᵀ Eᵀ(diag(K⁰𝟙) − K⁰) E v`. -/
private def formsKernel (net : Net) (G : Gram) (u v : Array Float) : Float × Float := Id.run do
  let g := net.geo
  let m := net.ne
  let d := fun (w : Array Float) (I : Nat) =>
    let (a, b) := net.edges[I]!
    ((w[b]! - w[a]!) / g.len[I]!) * g.tan[I]!
  let e := fun (w : Array Float) (I : Nat) =>
    let (a, b) := net.edges[I]!
    0.5 * (w[a]! + w[b]!)
  let mut sb : Float := 0.0
  let mut s0 : Float := 0.0
  for I in [0:m] do
    let mut rowB : Float := 0.0
    let mut row0 : Float := 0.0
    for J in [0:m] do
      rowB := rowB + G.kb.get! (I * m + J)
      row0 := row0 + G.k0.get! (I * m + J)
      sb := sb - G.kb.get! (I * m + J) * (d u I).dot (d v J)
      s0 := s0 - G.k0.get! (I * m + J) * e u I * e v J
    sb := sb + rowB * (d u I).dot (d v I)
    s0 := s0 + row0 * e u I * e v I
  return (sb, s0)

private def perturbPos (net : Net) (v c : Nat) (h : Float) : Net :=
  { net with pos := net.pos.set! v ((net.pos[v]!).set c ((net.pos[v]!).get c + h)) }

/-- Max relative disagreement of the differential with central differences. -/
private def fdEnergy (net : Net) (ex : Expo) : Float := Id.run do
  let (_, d) := energyAndDifferential net ex
  let mut scale : Float := 0.0
  for p in d do scale := max scale (max p.x.abs (max p.y.abs p.z.abs))
  let h := 1e-6 * (net.geo.len.foldl (· + ·) 0.0) / net.ne.toFloat
  let mut worst : Float := 0.0
  for v in [0:net.nv] do
    for c in [0:3] do
      let (ep, _) := energyAndDifferential (perturbPos net v c h) ex
      let (em, _) := energyAndDifferential (perturbPos net v c (-h)) ex
      worst := max worst (relErr ((ep - em) / (2.0 * h)) ((d[v]!).get c) scale)
  return worst

/-- Max disagreement of the constraint Jacobian with central differences, relative to its largest entry. -/
private def fdJacobian (net : Net) (planar : Option Nat) (rc : Array RCons) : Float := Id.run do
  let coords := coordsOf planar
  let rows := evalRows net net.pos coords rc
  let C := denseJacobian net.nv coords rows
  let nd := net.nv * coords.size
  let h := 1e-6 * (net.geo.len.foldl (· + ·) 0.0) / net.ne.toFloat
  let mut worst : Float := 0.0
  let scale := maxAbs C
  for v in [0:net.nv] do
    for s in [0:coords.size] do
      let rp := evalRows (perturbPos net v coords[s]! h) (perturbPos net v coords[s]! h).pos coords rc
      let rm := evalRows (perturbPos net v coords[s]! (-h)) (perturbPos net v coords[s]! (-h)).pos coords rc
      for k in [0:rows.size] do
        let fd := (rp[k]!.phi - rm[k]!.phi) / (2.0 * h)
        worst := max worst (relErr fd (C.get! (k * nd + v * coords.size + s)) scale)
  return worst

/-- Run the selftest on `count` seeded cases; returns the number of failures. -/
def selftest (count : Nat := 120) : IO Nat := do
  let mut bad := 0
  let mut checks := 0
  let fail := fun (bad : Nat) (what : String) => do
    IO.println s!"FAIL  {what}"
    return bad + 1
  -- the regular circle: the gradient of the energy is radial, equal at all vertices
  for (n, ex) in [(12, (⟨3.0, 6.0⟩ : Expo)), (17, ⟨2.0, 4.5⟩)] do
    let net := circle n 1.0
    let (_, d) := energyAndDifferential net ex
    let r0 := d[0]!.norm
    let mut worst : Float := 0.0
    for v in [0:n] do
      let p := net.pos[v]!
      let radial := (d[v]!.dot p)
      worst := max worst (relErr (d[v]!.norm) r0 r0)
      worst := max worst (relErr radial (r0 * (if radial < 0.0 then -1.0 else 1.0)) r0)
    checks := checks + 1
    if worst > 1e-12 then bad ← fail bad s!"circle {n}: gradient is not radial and equal ({worst})"
    -- fixed edge lengths and barycenter: the projected gradient vanishes
    let rc := (resolve net net.pos none #[.barycenter none, .edgeLengths none]).toOption.get!
    let rows := evalRows net net.pos #[0, 1, 2] rc
    match Saddle.build net.nv #[0, 1, 2] (gram net ex .midpoint).a rows with
    | .ok sad =>
      let gr := projectedGradient sad d
      checks := checks + 1
      if gr.gnorm > 1e-12 * r0 * 10.0 then bad ← fail bad s!"circle {n}: projected gradient {gr.gnorm} does not vanish"
    | .error e => bad ← fail bad s!"circle {n}: {e}"
  for i in [0:count] do
    let c := genCase 7 i
    let tag := c.name
    match common c with
    | .error e => bad ← fail bad s!"{tag}: generator produced an invalid case: {e}"
    | .ok (_, e0, d) =>
      let net := c.net
      let ex := c.ex
      -- differential against central differences (full R³)
      if c.planar.isNone then
        let w := fdEnergy net ex
        checks := checks + 1
        if w > 2e-6 then bad ← fail bad s!"{tag}: differential vs central differences {w}"
      -- scaling laws (SPEC §3.3, §4.5)
      let sc : Float := 1.7
      let net2 : Net := { net with pos := net.pos.map (V3.smul sc) }
      let (e2, _) := energyAndDifferential net2 ex
      checks := checks + 1
      if relErr e2 (e0 * sc.pow (2.0 + ex.alpha - ex.beta)) e0 > 1e-12 then bad ← fail bad s!"{tag}: energy scaling"
      let G1 := gram net ex c.quad
      let G2 := gram net2 ex c.quad
      let f := sc.pow (-(2.0 * ex.sigma + 1.0))
      let mut ws : Float := 0.0
      for k in [0:G1.b.size] do
        ws := max ws (relErr (G2.b.get! k) (f * G1.b.get! k) (maxAbs G1.b * f))
        ws := max ws (relErr (G2.b0.get! k) (f * G1.b0.get! k) (maxAbs G1.b0 * f))
      checks := checks + 1
      if ws > 1e-12 then bad ← fail bad s!"{tag}: Gram scaling {ws}"
      -- the increments against the definition of the quadratic form, and against Eq. 31
      let n := net.nv
      let mut u : Array Float := #[]
      let mut v : Array Float := #[]
      for k in [0:n] do
        u := u.push (Float.sin (1.3 * k.toFloat + 0.4))
        v := v.push (Float.cos (0.7 * k.toFloat * k.toFloat + 0.1))
      for q in [Quad.midpoint, Quad.trapezoid] do
        let G := gram net ex q
        let (sb, s0) := formsDirect net ex q u v
        let ub := bilinear n G.b u v
        let u0 := bilinear n G.b0 u v
        let (kb, k0) := formsKernel net G u v
        checks := checks + 1
        -- uᵀ B v with u,v scalar functions: the vector version is Ā = A ⊗ I₃, uᵀ Ā u = Σ uᵀ A u per coordinate
        if relErr ub sb (sb.abs + ub.abs) > 1e-10 then
          bad ← fail bad s!"{tag} {q.name}: B increments vs quadratic form ({ub} vs {sb})"
        if relErr u0 s0 (s0.abs + u0.abs) > 1e-10 then bad ← fail bad s!"{tag} {q.name}: B0 increments vs quadratic form ({u0} vs {s0})"
        if relErr kb sb (sb.abs + kb.abs) > 1e-10 then bad ← fail bad s!"{tag} {q.name}: Eq. 31 B vs quadratic form ({kb} vs {sb})"
        if relErr k0 s0 (s0.abs + k0.abs) > 1e-10 then bad ← fail bad s!"{tag} {q.name}: Eq. 31 B0 vs quadratic form ({k0} vs {s0})"
        -- symmetry
        let mut asym : Float := 0.0
        for a in [0:n] do
          for b in [0:n] do
            asym := max asym ((G.b.get! (a * n + b) - G.b.get! (b * n + a)).abs)
            asym := max asym ((G.b0.get! (a * n + b) - G.b0.get! (b * n + a)).abs)
        if asym > 1e-12 * (maxAbs G.b + maxAbs G.b0) then bad ← fail bad s!"{tag} {q.name}: Gram matrix not symmetric ({asym})"
      -- constraints
      if c.cons.size > 0 then
        match resolve net net.pos c.planar c.cons with
        | .error err => bad ← fail bad s!"{tag}: {err}"
        | .ok rc =>
          let jd := fdJacobian net c.planar rc
          checks := checks + 1
          if jd > 2e-6 then bad ← fail bad s!"{tag}: constraint Jacobian vs central differences {jd}"
          -- saddle equations
          let coords := coordsOf c.planar
          let rows := evalRows net net.pos coords rc
          match Saddle.build n coords (gram net ex c.quad).a rows with
          | .error err => IO.println s!"note  {tag}: {err}"
          | .ok sad =>
            let gr := projectedGradient sad d
            let gx := gather coords gr.g
            let a := gather coords d
            -- C g = 0 and Ā g + Cᵀ λ = dF
            let mut r1 : Float := 0.0
            for k in [0:sad.k] do
              let mut s : Float := 0.0
              for j in [0:sad.nd] do s := s + sad.c.get! (k * sad.nd + j) * gx.get! j
              r1 := max r1 s.abs
            let A := (gram net ex c.quad).a
            let nc := coords.size
            let mut r2 : Float := 0.0
            for vv in [0:n] do
              for s in [0:nc] do
                let mut acc : Float := 0.0
                for w in [0:n] do acc := acc + A.get! (vv * n + w) * gx.get! (w * nc + s)
                for k in [0:sad.k] do acc := acc + sad.c.get! (k * sad.nd + vv * nc + s) * gr.lam.get! k
                r2 := max r2 (acc - a.get! (vv * nc + s)).abs
            checks := checks + 2
            let gscale := maxAbs a + 1e-300
            if r2 > 1e-8 * gscale then bad ← fail bad s!"{tag}: saddle residual (stationarity) {r2} of {gscale}"
            if r1 > 1e-8 * (gr.gnorm + 1e-300) + 1e-14 then bad ← fail bad s!"{tag}: C g̃ = {r1} not zero (|g̃| = {gr.gnorm})"
            if gr.dfg < -1e-9 * gscale * gr.gnorm then bad ← fail bad s!"{tag}: dF·g̃ = {gr.dfg} < 0, not a descent direction"
            -- the projection of the unmoved curve changes nothing, and a small step projects back
            let lbar := (net.geo.len.foldl (· + ·) 0.0) / net.ne.toFloat
            let trial := (Array.range n).map fun vv => net.pos[vv]! - V3.smul (0.01 * lbar / (maxAbs gx + 1e-300)) gr.g[vv]!
            let pr := chordProject net rc sad trial (1e-9 * lbar) lbar 8
            -- only from a feasible base curve (a perturbed tangent direction starts infeasible)
            let feasible := rows.all fun r => r.phi.abs ≤ 1e-9 * (if r.unit == 1 then 1.0 else lbar) * (if r.unit == 2 then r.gscale else 1.0)
            if feasible then checks := checks + 1
            match pr.status with
            | _ => if !feasible then pure () else match pr.status with
              | .converged => pure ()
              | .tangentFlipped _ => pure ()
              | _ => bad ← fail bad s!"{tag}: projection of a 1% step did not converge to 1e-9 in 8 iterations ({pr.history})"
  IO.println s!"rc selftest: {checks} checks, {bad} failures"
  return bad

end Ddg.RepulsiveCurves
