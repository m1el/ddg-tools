import DdgOracle.RepulsiveSurfaces.Gen

/-!
# Internal consistency checks of the `rs` oracle

Each check compares two independent routes through the SPEC: central differences against the exact
gradient and the Jacobians, the pair-sum matrices against the kernel-matrix form of §4.3, the nested
saddle against the full saddle solved by Gauss-Jordan elimination, scaling laws, remeshing invariants.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

private def relErr (a b scale : Float) : Float := (a - b).abs / (if scale > 1e-300 then scale else 1.0)

private structure Tally where
  bad : Nat := 0
  total : Nat := 0

private def report (tally : IO.Ref Tally) (name : String) (ok : Bool) (detail : String) : IO Unit := do
  tally.modify fun t => { bad := t.bad + (if ok then 0 else 1), total := t.total + 1 }
  if !ok then IO.println s!"FAIL  {name}: {detail}"

private def randVec (rng : Rng) (n : Nat) : Array Float × Rng := Id.run do
  let mut r := rng
  let mut a : Array Float := #[]
  for _ in [0:n] do
    let (x, r2) := r.range (-1.0) 1.0
    a := a.push x
    r := r2
  return (a, r)

private def bilinear (n : Nat) (m : Mat) (u v : Array Float) : Float := Id.run do
  let mut s := 0.0
  for i in [0:n] do
    for j in [0:n] do
      s := s + u[i]! * m.get! (i * n + j) * v[j]!
  return s

/-- Kernel-matrix form (SPEC §4.3): `2 xᵀ (diag(a)⁻¹ diag(H a) − H) y` per component, with `X u = a ⊙ (face values)`. -/
private def kernelForm (a : Array Float) (h : Array (Array Float)) (xs ys : Array (Array Float)) : Float := Id.run do
  let nf := a.size
  let mut total := 0.0
  for c in [0:xs.size] do
    let x := (Array.range nf).map fun s => a[s]! * xs[c]![s]!
    let y := (Array.range nf).map fun s => a[s]! * ys[c]![s]!
    for s in [0:nf] do
      let mut ha := 0.0
      for t in [0:nf] do ha := ha + h[s]![t]! * a[t]!
      total := total + 2.0 * x[s]! * (ha / a[s]!) * y[s]!
      for t in [0:nf] do total := total - 2.0 * x[s]! * h[s]![t]! * y[t]!
  return total

private def faceScalar (t : Topo) (u : Array Float) : Array Float :=
  t.faces.map fun f => (u[f.a]! + u[f.b]! + u[f.c]!) / 3.0

private def faceGrad (t : Topo) (pos : Array V3) (g : Array FaceGeo) (u : Array Float) : Array (Array Float) :=
  let ds := (Array.range t.faces.size).map fun s => dfFace pos t.faces[s]! g[s]! u
  #[ds.map (·.x), ds.map (·.y), ds.map (·.z)]

private def checkMesh (tally : IO.Ref Tally) (label : String) (m : Mesh) (ex : Ex) (seed : Nat) : IO Unit := do
  let rpt := report tally
  match m.validate with
  | .error e => rpt s!"{label} valid" false e
  | .ok t =>
    let g := faceGeo t.faces m.pos
    let nv := t.nv
    let nf := t.faces.size
    let lbar := meanEdge t m.pos
    let mut rng : Rng := ⟨seed.toUInt64 * 7919 + 17⟩
    -- energy gradient by central differences
    match energyGrad t.faces m.pos ex with
    | .error e => rpt s!"{label} energy" false e
    | .ok (e, grad) =>
      let mut gmax := 0.0
      for p in grad do gmax := max gmax (max p.x.abs (max p.y.abs p.z.abs))
      let h := 1.0e-6 * lbar
      let mut worst := 0.0
      for _ in [0:8] do
        let (vi, r1) := rng.below nv
        let (ax, r2) := r1.below 3
        rng := r2
        let bump := fun (s : Float) => m.pos.modify vi fun p => p.set ax (p.get ax + s)
        let ep := (energyGrad t.faces (bump h) ex).toOption.map (·.1) |>.getD (0.0 / 0.0)
        let em := (energyGrad t.faces (bump (-h)) ex).toOption.map (·.1) |>.getD (0.0 / 0.0)
        let fd := (ep - em) / (2.0 * h)
        worst := max worst (relErr fd ((grad[vi]!).get ax) gmax)
      rpt s!"{label} energy FD" (worst < 2e-5) s!"relative error {worst}"
      -- identities
      let mut sg := V3.zero
      let mut sm := V3.zero
      let mut sd := 0.0
      for v in [0:nv] do
        sg := sg + grad[v]!
        sm := sm + m.pos[v]!.cross grad[v]!
        sd := sd + m.pos[v]!.dot grad[v]!
      let es := (4.0 + ex.alpha - ex.beta) * e
      rpt s!"{label} sum g = 0" (sg.norm ≤ 1e-9 * (gmax * nv.toFloat + 1e-300)) s!"{sg.norm}"
      rpt s!"{label} sum f x g = 0" (sm.norm ≤ 1e-9 * (gmax * nv.toFloat * (1.0 + lbar * 10.0) + 1e-300)) s!"{sm.norm}"
      rpt s!"{label} Euler identity" (relErr sd es (e.abs) < 1e-9) s!"{sd} vs {es}"
      let e2 := (energy t.faces (faceGeo t.faces (m.pos.map fun p => (2.0 : Float) * p)) ex).toOption.getD 0.0
      let expect := Float.pow 2.0 (4.0 + ex.alpha - ex.beta) * e
      rpt s!"{label} scaling" (relErr e2 expect e.abs < 1e-10) s!"{e2} vs {expect}"
    -- D_f
    let (cv, r1) := randVec rng 3
    rng := r1
    let cvec : V3 := ⟨cv[0]!, cv[1]!, cv[2]!⟩
    let lin := m.pos.map fun p => cvec.dot p
    let mut dfErr := 0.0
    for s in [0:nf] do
      let d := dfFace m.pos t.faces[s]! g[s]! lin
      let want := cvec - (cvec.dot g[s]!.n) * g[s]!.n
      dfErr := max dfErr (d - want).norm
    rpt s!"{label} D_f of a linear function" (dfErr < 1e-10 * (1.0 + cvec.norm)) s!"{dfErr}"
    let ones := Array.replicate nv 1.0
    let mut d1 := 0.0
    for s in [0:nf] do d1 := max d1 (dfFace m.pos t.faces[s]! g[s]! ones).norm
    rpt s!"{label} D_f 1 = 0" (d1 < 1e-12) s!"{d1}"
    -- metric
    match buildOps t m.pos ex { wtp := 1.0, clo := 1.0, ww := 0.3 } with
    | .error e => rpt s!"{label} metric" false e
    | .ok ops =>
      match assemblePairs t m.pos 0.6 ex.s with
      | .error e => rpt s!"{label} pairs" false e
      | .ok p =>
        let (u, r2) := randVec rng nv
        let (v, r3) := randVec r2 nv
        rng := r3
        let a := g.map (·.area)
        let rho := fun (s t' : Nat) => (g[s]!.x - g[t']!.x).norm
        let hL := (Array.range nf).map fun s => (Array.range nf).map fun t' => if s == t' then 0.0 else Float.pow (rho s t') (-(2.0 * 0.6 + 2.0))
        let hB := (Array.range nf).map fun s => (Array.range nf).map fun t' => if s == t' then 0.0 else Float.pow (rho s t') (-2.0 * ex.s)
        let k2 := fun (s t' : Nat) => let nu := g[s]!.n.dot (g[s]!.x - g[t']!.x); nu * nu / Float.pow (rho s t') 4.0
        let h0 := (Array.range nf).map fun s => (Array.range nf).map fun t' =>
          if s == t' then 0.0 else 0.5 * (k2 s t' + k2 t' s) * Float.pow (rho s t') (-2.0 * ex.s)
        let su := faceScalar t u
        let sv := faceScalar t v
        let qL := kernelForm a hL #[su] #[sv]
        let qB := kernelForm a hB (faceGrad t m.pos g u) (faceGrad t m.pos g v)
        let q0 := kernelForm a h0 #[su] #[sv]
        rpt s!"{label} L form" (relErr (bilinear nv p.l u v) qL (qL.abs + 1e-300) < 1e-9) s!"{bilinear nv p.l u v} vs {qL}"
        rpt s!"{label} B form" (relErr (bilinear nv ops.b u v) qB (qB.abs + 1e-300) < 1e-9) s!"{bilinear nv ops.b u v} vs {qB}"
        rpt s!"{label} B0 form" (relErr (bilinear nv ops.b0 u v) q0 (q0.abs + 1e-300) < 1e-9) s!"{bilinear nv ops.b0 u v} vs {q0}"
        -- symmetry, null space, PSD
        let mut asym := 0.0
        let mut rowsum := 0.0
        let mut amax := 0.0
        for i in [0:nv] do
          let mut rs := 0.0
          for j in [0:nv] do
            asym := max asym ((ops.a.get! (i * nv + j) - ops.a.get! (j * nv + i)).abs)
            amax := max amax (ops.a.get! (i * nv + j)).abs
            rs := rs + ops.a.get! (i * nv + j)
          rowsum := max rowsum rs.abs
        rpt s!"{label} A symmetric" (asym ≤ 1e-10 * amax) s!"{asym} (max {amax})"
        rpt s!"{label} A 1 = 0" (rowsum ≤ 1e-9 * amax) s!"{rowsum} (max {amax})"
        let mut minq := 0.0
        for mm in [p.l, ops.b, ops.b0, ops.a] do
          let mut sc := 0.0
          for i in [0:nv * nv] do sc := max sc (mm.get! i).abs
          for _ in [0:4] do
            let (w, r4) := randVec rng nv
            rng := r4
            minq := min minq (bilinear nv mm w w / (sc + 1e-300))
        rpt s!"{label} PSD" (minq > -1e-10) s!"min uᵀMu/scale = {minq}"
        let mut msum := 0.0
        for x in ops.m do msum := msum + x
        let mut tot := 0.0
        for x in a do tot := tot + x
        rpt s!"{label} masses sum to the area" (relErr msum tot tot < 1e-12) s!"{msum} vs {tot}"
        let mut wrow := 0.0
        for i in [0:nv] do
          let mut rs := 0.0
          for j in [0:nv] do rs := rs + ops.w.get! (i * nv + j)
          wrow := max wrow rs.abs
        rpt s!"{label} W 1 = 0" (wrow < 1e-10) s!"{wrow}"
    -- constraints and saddle
    let (dxs, r5) := randVec rng 3
    rng := r5
    let shift : V3 := ⟨dxs[0]!, dxs[1]!, dxs[2]!⟩
    let rcAll := (resolveCons t m.pos #[.barycenter none none, .area none none]).toOption.getD #[]
    -- translation covariance of the barycenter, frozen-weight Jacobian
    for k in rcAll do
      match k with
      | .bary comp tgt =>
        let rows := evalCons t m.pos k
        let moved := evalCons t (m.pos.map (· + shift)) k
        let mut err := 0.0
        for ax in [0:3] do err := max err ((moved[ax]!.phi - rows[ax]!.phi) - shift.get ax).abs
        -- Jacobian times the translation of the component = the shift
        let mut jerr := 0.0
        for ax in [0:3] do
          let mut s := 0.0
          for (j, w) in rows[ax]!.jac do if j % 3 == ax then s := s + w * shift.get ax
          jerr := max jerr (s - shift.get ax).abs
        rpt s!"{label} barycenter {comp} translation" (err < 1e-12 && jerr < 1e-12) s!"{err} {jerr}"
        let _ := tgt
      | _ => pure ()
    -- area and volume Jacobians by central differences
    let mut kinds : Array Cons := #[.area none none]
    for c in [0:t.ncomp] do
      kinds := kinds.push (.area (some c) none)
      if (inconsistentEdge t c).isNone then kinds := kinds.push (.volume (some c) none)
    if (Array.range t.ncomp).all (fun c => (inconsistentEdge t c).isNone) then kinds := kinds.push (.volume none none)
    for kc in kinds do
      match resolveCons t m.pos #[kc] with
      | .error e => rpt s!"{label} resolve" false e
      | .ok rc =>
        let rows := evalCons t m.pos rc[0]!
        let h := 1.0e-6 * lbar
        let mut worst := 0.0
        let mut jmax := 0.0
        for (_, w) in rows[0]!.jac do jmax := max jmax w.abs
        for _ in [0:6] do
          let (vi, r1) := rng.below nv
          let (ax, r2) := r1.below 3
          rng := r2
          let bump := fun (s : Float) => m.pos.modify vi fun p => p.set ax (p.get ax + s)
          let pp := (evalCons t (bump h) rc[0]!)[0]!.phi
          let pm := (evalCons t (bump (-h)) rc[0]!)[0]!.phi
          let mut jv := 0.0
          for (j, w) in rows[0]!.jac do if j == 3 * vi + ax then jv := jv + w
          worst := max worst (relErr ((pp - pm) / (2.0 * h)) jv jmax)
        rpt s!"{label} {(rc[0]!).name} Jacobian FD" (worst < 1e-6) s!"{worst}"
    -- the nested saddle against the full saddle
    let (posTargets, r6) := randVec rng 3
    rng := r6
    let mut cons : Array Cons := #[]
    for c in [0:t.ncomp] do
      let vs := (Array.range nv).filter (t.vertComp[·]! == c)
      cons := cons.push (.pins #[vs[0]!] (some #[m.pos[vs[0]!]! + (0.01 : Float) * (⟨posTargets[0]!, posTargets[1]!, posTargets[2]!⟩ : V3)]))
    if t.ncomp == 1 then cons := cons.push (.barycenter none none)
    let g0 := faceGeo t.faces m.pos
    cons := cons.push (.area none (some (totalArea t g0 none * 1.03)))
    if (Array.range t.ncomp).all (fun c => (inconsistentEdge t c).isNone) && (totalVolume t m.pos none).abs > 1e-3 then
      cons := cons.push (.volume none (some (totalVolume t m.pos none * 0.97)))
    match resolveCons t m.pos cons with
    | .error e => rpt s!"{label} saddle constraints" false e
    | .ok rc =>
      match energyGrad t.faces m.pos ex, buildOps t m.pos ex {} with
      | .ok (_, grad), .ok ops =>
        let rows := evalAll t m.pos rc
        match buildNested nv ops.a rows with
        | .error e => rpt s!"{label} nested saddle" false e
        | .ok nst =>
          let gf := flatten grad
          let x := nst.direction gf
          let n3 := 3 * nv
          let nrows := (rows.filter (fun r => !r.2.1)).map (·.2.2)
          let phiN := FloatArray.ofFn' nrows.size fun i => nrows[i]!.phi
          let h := nst.correction phiN
          let cx := maxAbs (nst.applyCn x)
          let mut cqx := 0.0
          for r in [0:nst.ap.kq] do
            let mut s := 0.0
            for j in [0:n3] do s := s + nst.ap.cq.get! (r * n3 + j) * x.get! j
            cqx := max cqx s.abs
          let xs := Float.sqrt (dotF x x) + 1e-300
          rpt s!"{label} C_N x = 0" (cx ≤ 1e-8 * xs * (1.0 + maxAbs nst.cn)) s!"{cx} (|x| {xs})"
          rpt s!"{label} C_Q x = 0" (cqx ≤ 1e-8 * xs * (1.0 + maxAbs nst.ap.cq)) s!"{cqx}"
          let ch := nst.applyCn h
          let mut herr := 0.0
          for i in [0:nrows.size] do herr := max herr (ch.get! i - phiN.get! i).abs
          rpt s!"{label} C_N h = Phi" (herr ≤ 1e-8 * (1.0 + maxAbs phiN)) s!"{herr}"
          -- full saddle by Gauss-Jordan
          let kall := rows.size
          let nn := n3 + kall
          let mut mfull := zeros (nn * nn)
          for v in [0:nv] do
            for w in [0:nv] do
              for d in [0:3] do
                mfull := mfull.set! ((3 * v + d) * nn + 3 * w + d) (ops.a.get! (v * nv + w))
          for r in [0:kall] do
            let row := denseRow n3 rows[r]!.2.2
            for j in [0:n3] do
              mfull := mfull.set! ((n3 + r) * nn + j) (row.get! j)
              mfull := mfull.set! (j * nn + n3 + r) (row.get! j)
          let rhs := FloatArray.ofFn' nn fun i => if i < n3 then gf.get! i else 0.0
          match gaussJordanSolve nn mfull rhs with
          | none => rpt s!"{label} full saddle" false "Gauss-Jordan failed"
          | some sol =>
            let mut dmax := 0.0
            for j in [0:n3] do dmax := max dmax (sol.get! j - x.get! j).abs
            rpt s!"{label} nested = full saddle" (dmax ≤ 1e-7 * xs) s!"{dmax} (|x| {xs})"
          -- preconditioner
          let w := stiffness t m.pos
          let (_, lmat) := buildPrecondMatrix nv w ops.m 1e-10
          let qrows := (rows.filter (·.2.1)).map (·.2.2)
          match buildAp nv lmat qrows with
          | .error e => rpt s!"{label} L'" false e
          | .ok ap =>
            match assemblePairs t m.pos (2.0 - ex.s) ex.s with
            | .error e => rpt s!"{label} L_fr" false e
            | .ok pr =>
              let (r, r7) := randVec rng n3
              let (r', r8) := randVec r7 n3
              rng := r8
              let fa := FloatArray.mk r
              let fb := FloatArray.mk r'
              let za := applyPrecond nv ap pr.l fa
              let zb := applyPrecond nv ap pr.l fb
              let s1 := dotF fb za
              let s2 := dotF fa zb
              rpt s!"{label} P^-1 symmetric" (relErr s1 s2 (s1.abs + s2.abs + 1e-300) < 1e-6) s!"{s1} vs {s2}"
              let mut cz := 0.0
              for k in [0:ap.kq] do
                let mut s := 0.0
                for j in [0:n3] do s := s + ap.cq.get! (k * n3 + j) * za.get! j
                cz := max cz s.abs
              rpt s!"{label} C_Q P^-1 r = 0" (cz ≤ 1e-6 * (Float.sqrt (dotF za za) + 1e-300) * (1.0 + maxAbs ap.cq)) s!"{cz}"
          -- one step, with area alone and with volume alone (together they are nearly dependent on spheres)
          let posRc := rc.filter (·.positional)
          let nlSets : Array (Array Cons) :=
            #[#[.area none none]] ++
            (if (Array.range t.ncomp).all (fun c => (inconsistentEdge t c).isNone) && (totalVolume t m.pos none).abs > 1e-3
             then #[#[.volume none none]] else #[])
          for nl in nlSets do
            match resolveCons t m.pos nl with
            | .error e => rpt s!"{label} step constraints" false e
            | .ok nrc =>
              let rcs := posRc ++ nrc
              match projectedStep t ex {} rcs {} m.pos with
              | .error e => rpt s!"{label} step" false e
              | .ok sr =>
                rpt s!"{label} step slope > 0" (sr.slope > 0.0) s!"{sr.slope}"
                if sr.accepted then rpt s!"{label} step decreases the energy" (sr.e1 < sr.f0) s!"{sr.e1} vs {sr.f0}"
                -- the projection reduces the nonlinear violation (old Jacobian: linear convergence)
                let before := sr.phiBefore.foldl (fun acc x => max acc x.abs) 0.0
                let after := sr.phi.foldl (fun acc x => max acc x.abs) 0.0
                rpt s!"{label} step projection" (after ≤ 1e-7 * (1.0 + before) || after ≤ 0.1 * before) s!"{before} -> {after}"
      | _, _ => rpt s!"{label} saddle setup" false "energy or ops failed"
    -- remeshing
    for e in t.edges do
      let (a, b, _) := e
      let area0 := totalArea t g none
      match splitEdge m t a b with
      | .error err => rpt s!"{label} split {a} {b}" false err
      | .ok o =>
        match buildTopo o.pos.size o.faces with
        | .error err => rpt s!"{label} split {a} {b} valid" false err
        | .ok t2 =>
          let a2 := totalArea t2 (faceGeo t2.faces o.pos) none
          rpt s!"{label} split {a} {b}" (o.pos.size == nv + 1 && o.faces.size == nf + e.2.2.size && relErr a2 area0 area0 < 1e-12 && t2.ncomp == t.ncomp) s!"{o.pos.size} {o.faces.size} {a2} {area0}"
          let cons1 := (Array.range t.ncomp).all fun c => (inconsistentEdge t c).isNone
          let cons2 := (Array.range t2.ncomp).all fun c => (inconsistentEdge t2 c).isNone
          rpt s!"{label} split {a} {b} orientation" (cons1 == cons2) ""
      match flipEdge m t a b false false with
      | .error _ => pure ()
      | .ok o =>
        match buildTopo o.pos.size o.faces with
        | .error err => rpt s!"{label} flip {a} {b} valid" false err
        | .ok t2 =>
          let a2 := totalArea t2 (faceGeo t2.faces o.pos) none
          let _ := (a2, area0)
          rpt s!"{label} flip {a} {b}" (o.faces.size == nf && t2.edges.size == t.edges.size) s!"{a2} {area0}"
      match collapseEdge m t a b #[] false with
      | .error _ => pure ()
      | .ok o =>
        match buildTopo o.pos.size o.faces with
        | .error err => rpt s!"{label} collapse {a} {b} valid" false err
        | .ok t2 =>
          rpt s!"{label} collapse {a} {b}" (o.pos.size == nv - 1 && o.faces.size + e.2.2.size == nf && t2.ncomp == t.ncomp) s!"{o.pos.size} {o.faces.size}"

/-- Analytic checks on fixed meshes. -/
private def analytic (tally : IO.Ref Tally) : IO Unit := do
  let rpt := report tally
  -- a flat patch: N·(X−Y) = 0 for every pair, so the energy and its gradient vanish
  let flat := grid 3 3 0.0
  match flat.validate with
  | .error e => rpt "flat patch valid" false e
  | .ok t =>
    match energyGrad t.faces flat.pos ⟨6.0, 12.0⟩ with
    | .error e => rpt "flat patch energy" false e
    | .ok (e, gr) =>
      let mut gm := 0.0
      for p in gr do gm := max gm p.norm
      rpt "flat patch energy = 0" (e == 0.0 && gm == 0.0) s!"{e} {gm}"
  -- the outward faces of the closed meshes are consistent and enclose positive volume
  for (nm, m) in [("tetra", tetra), ("octa", octa), ("icosa", icosa), ("torus", torus 5 3 2.0 0.8)] do
    match m.validate with
    | .error e => rpt s!"{nm} valid" false e
    | .ok t =>
      let ok := (Array.range t.ncomp).all fun c => (inconsistentEdge t c).isNone
      rpt s!"{nm} consistent" ok ""
      rpt s!"{nm} volume > 0" (totalVolume t m.pos none > 0.0) s!"{totalVolume t m.pos none}"
  match (mobius 8).validate with
  | .error e => rpt "mobius valid" false e
  | .ok t => rpt "mobius is not consistent" (!((inconsistentEdge t 0).isNone)) ""
  -- tetrahedron: the area of the unit-edge-scaled regular tetrahedron with vertices (±1, ±1, ±1) is 4 * sqrt 3 * (2 sqrt 2)^2 / 4... check numerically: edge 2 sqrt 2
  match tetra.validate with
  | .error _ => pure ()
  | .ok t =>
    let area := totalArea t (faceGeo t.faces tetra.pos) none
    let edge := 2.0 * Float.sqrt 2.0
    rpt "tetra area" (relErr area (Float.sqrt 3.0 * edge * edge) area < 1e-14) s!"{area}"
    rpt "tetra volume" (relErr (totalVolume t tetra.pos none) (8.0 / 3.0) 3.0 < 1e-14) s!"{totalVolume t tetra.pos none}"
  -- error cases
  let degenerate : Mesh := { pos := #[⟨0.0, 0.0, 0.0⟩, ⟨1.0, 0.0, 0.0⟩, ⟨2.0, 0.0, 0.0⟩], faces := #[⟨0, 1, 2⟩] }
  rpt "degenerate face is an error" (match degenerate.validate with | .error e => e.startsWith "DegenerateFace" | .ok _ => false) ""
  let nonman : Mesh := { pos := #[⟨0.0, 0.0, 0.0⟩, ⟨1.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩, ⟨0.0, 0.0, 1.0⟩, ⟨0.0, -1.0, 0.0⟩],
                         faces := #[⟨0, 1, 2⟩, ⟨0, 1, 3⟩, ⟨0, 1, 4⟩] }
  rpt "non-manifold edge is an error" (match nonman.validate with | .error e => e.startsWith "NonManifoldEdge" | .ok _ => false) ""
  rpt "exponent range" ((Ex.check ⟨6.0, 8.0⟩).isOk == false && (Ex.check ⟨6.0, 12.0⟩).isOk) ""

/-- Run all checks on `n` seeded meshes; returns the number of failures. -/
def selftest (n : Nat := 120) : IO Nat := do
  let tally ← IO.mkRef ({} : Tally)
  analytic tally
  -- cases are spread over the meshes of the generator; every mesh kind is hit
  for i in [0:n] do
    let c := genCase 11 i
    checkMesh tally s!"[{i} {c.name}]" c.mesh c.ex i
  for i in [0:200] do
    let c := genCase 5 i
    let lines := ((printCase c).splitOn "\n").toArray.filterMap fun l =>
      let tk := tokens l
      if tk.isEmpty || tk[0]! == "case" || tk[0]! == "end" then none else some tk
    let ok := match parseBody c.name lines with
      | .error _ => false
      | .ok c2 => printCase c2 == printCase c && runCase c2 == runCase c
    report tally s!"generated case {i} round trip" ok ""
  let t ← tally.get
  IO.println s!"rs selftest: {t.total} checks, {t.bad} failures"
  return t.bad

end Ddg.RepulsiveSurfaces
