import DdgOracle.PointsAsTori.Gen

/-!
# Internal consistency checks of the `pat` oracle

Two independent routes through the SPEC: central differences against the closed-form gradients and
Laplacians, Rayleigh quotients against Eq. 19, the torus contact conditions (`pat-torus-contact`),
translation / scale invariances, and print / parse round trips.
-/

namespace Ddg.PointsAsTori

private structure Tally where
  bad : Nat := 0
  total : Nat := 0

private def report (tally : IO.Ref Tally) (name : String) (ok : Bool) (detail : String) : IO Unit := do
  tally.modify fun t => { bad := t.bad + (if ok then 0 else 1), total := t.total + 1 }
  if !ok then IO.println s!"FAIL  {name}: {detail}"

private def rel (a b : Float) : Float := (a - b).abs / (if a.abs > b.abs then (if a.abs > 1.0 then a.abs else 1.0) else (if b.abs > 1.0 then b.abs else 1.0))

private def relV (a b : V3) : Float :=
  let m := (a - b).norm
  m / (let s := if a.norm > b.norm then a.norm else b.norm; if s > 1.0 then s else 1.0)

private def ex : V3 := ⟨1.0, 0.0, 0.0⟩
private def ey : V3 := ⟨0.0, 1.0, 0.0⟩
private def ez : V3 := ⟨0.0, 0.0, 1.0⟩

private def fdGrad (f : V3 → Float) (x : V3) (h : Float) : V3 :=
  let d := fun (e : V3) => (f (x + e.smul h) - f (x - e.smul h)) / (2.0 * h)
  ⟨d ex, d ey, d ez⟩

private def fdLap (f : V3 → Float) (x : V3) (h : Float) : Float :=
  let d2 := fun (e : V3) => (f (x + e.smul h) - 2.0 * f x + f (x - e.smul h)) / (h * h)
  d2 ex + d2 ey + d2 ez

/-- Run the checks. Returns the number of failures. -/
def selftest (n : Nat := 60) : IO Nat := do
  let tally ← IO.mkRef ({} : Tally)
  -- frame
  for k in [0:n] do
    let (c, _) := (unitVec.run (mkRng k))
    let (s, t) := frame c
    let ok := (s.norm - 1.0).abs < 1e-12 && (t.norm - 1.0).abs < 1e-12 && (s.dot t).abs < 1e-12 &&
      (s.dot c).abs < 1e-12 && (s.cross t - c).norm < 1e-12
    report tally s!"frame {k}" ok s!"{c.x} {c.y} {c.z}"
  for c in [ex, ey, ez, (⟨0.95, 0.0, Float.sqrt (1.0 - 0.9025)⟩ : V3), (⟨-0.95, 0.0, Float.sqrt (1.0 - 0.9025)⟩ : V3)] do
    let (s, t) := frame c
    report tally s!"frame axis" ((s.cross t - c).norm < 1e-12 && (s.dot t).abs < 1e-12) s!"{c.x}"
  -- features: invariances under translation and uniform scaling
  for k in [0:n] do
    let act : G (Array V3 × Array V3 × V3 × Float) := do
      let mut ps := #[]
      let mut ns := #[]
      for _ in [0:9] do
        ps := ps.push (← vecIn (-1.0) 1.0)
        ns := ns.push (← unitVec)
      return (ps, ns, ← vecIn (-3.0) 3.0, ← uni 0.2 5.0)
    let ((ps, ns, tr, sc), _) := act.run (mkRng (1000 + k))
    let idx := (Array.range 9)
    let f0 := features ps ns idx
    let ps2 := ps.map fun p => (p + tr).smul sc
    let f1 := features ps2 ns idx
    let mut worst := rel f1.sigma (sc * f0.sigma)
    for m in [0:9] do
      for j in [0:6] do
        let e := (f0.rows[m]![j]! - f1.rows[m]![j]!).abs
        if e > worst then worst := e
    report tally s!"features invariance {k}" (worst < 1e-8) s!"{worst}"
  -- median
  report tally "median even" (median #[3.0, 1.0, 2.0, 10.0] == 2.5) ""
  report tally "median odd" (median #[3.0, 1.0, 2.0] == 2.0) ""
  -- fitting
  for k in [0:n * 4] do
    let act : G (V3 × V3 × Array Float) := do return (← vecIn (-1.0) 1.0, ← unitVec, ← randCoef)
    let ((p, nn, a), _) := act.run (mkRng (2000 + k))
    let r := fitTorus {} p nn a
    let o := p + nn.smul a[0]!
    let s := r.torus.sdf o
    -- passes through o, tube normal is ±n*, Rayleigh = Eq. 19 as a set
    let scale := 1.0 / r.kMin.abs
    let c1 := s.g.abs / (if scale > 1.0 then scale else 1.0)
    let dir := if s.kind == 0 then (s.grad.smul (sgn0 (s.grad.dot r.nStar)) - r.nStar).norm else 0.0
    let ks := (r.kPlus - r.kPaperPlus).abs + (r.kMinus - r.kPaperMinus).abs
    let ks2 := (r.kPlus - r.kPaperMinus).abs + (r.kMinus - r.kPaperPlus).abs
    let kk := (if ks < ks2 then ks else ks2) / (r.kMax.abs + r.kMin.abs)
    report tally s!"fit contact {k}" (c1 < 1e-9 && dir < 1e-7 && kk < 1e-9 && !r.clampDisc)
      s!"g(o)={s.g} scale={scale} dir={dir} kk={kk}"
    report tally s!"fit unit axis {k}" ((r.torus.u.norm - 1.0).abs < 1e-12 && (r.torus.u.dot r.nStar).abs < 1e-12) ""
    -- tangent plane at o is spanned by v±; v± orthogonal to n*
    report tally s!"fit v perp n* {k}" ((r.vMin.dot r.nStar).abs < 1e-12 && (r.vMax.dot r.nStar).abs < 1e-12) ""
  -- sphere init
  let (sh, th) := frame ez
  report tally "sphere" (let r := fitTorus {} ⟨0.0, 0.0, 0.0⟩ ez #[0.0, 0.0, 0.0, 0.0, -0.5, -0.5]
      (r.torus.c - (⟨0.0, 0.0, -1.0⟩ : V3)).norm < 1e-12 && r.torus.R == 0.0 && (r.torus.r - 1.0).abs < 1e-12 && r.clampDisc) s!"{sh.x}{th.x}"
  -- SDF gradient and Laplacian against differences
  for k in [0:n * 3] do
    let act : G (Torus × V3) := do
      let c ← vecIn (-1.0) 1.0
      let u ← unitVec
      let R ← uni 0.0 2.0
      let r0 ← uni 0.1 1.2
      let sg ← below 2
      let x ← vecIn (-3.0) 3.0
      return ({ c, u, R, r := if sg == 0 then r0 else -r0 }, x)
    let ((t, x), _) := act.run (mkRng (3000 + k))
    let s := t.sdf x
    if s.kind == 0 then
      let fg := fdGrad (fun y => (t.sdf y).g) x 1e-6
      let fl := fdLap (fun y => (t.sdf y).g) x 1e-4
      -- the Laplacian check needs a point away from the singular sets
      let xc := x - t.c
      let v := (xc.cross t.u).norm
      let z := xc.dot t.u
      let dn := Float.sqrt ((v - t.R) * (v - t.R) + z * z)
      report tally s!"sdf gradient {k}" (relV fg s.grad < 1e-7) s!"{relV fg s.grad}"
      if dn > 0.15 && v > 0.15 then
        report tally s!"sdf laplacian {k}" (rel fl s.lap < 1e-5) s!"fd {fl} exact {s.lap}"
  -- blend and baselines
  for k in [0:n / 2] do
    let act : G (Array V3 × Array V3 × Array (Array Float) × V3 × Array Float) := do
      let mut ps := #[]
      let mut ns := #[]
      let mut cs := #[]
      let mut ar := #[]
      for _ in [0:12] do
        ps := ps.push (← vecIn (-1.0) 1.0)
        ns := ns.push (← unitVec)
        cs := cs.push (← randCoef)
        ar := ar.push (← uni 0.5 1.5)
      return (ps, ns, cs, ← vecIn (-1.5) 1.5, ar)
    let ((ps, ns, cs, x, ar), _) := act.run (mkRng (4000 + k))
    let tori := (Array.range 12).map fun i => (fitTorus {} ps[i]! ns[i]! cs[i]!).torus
    -- constant lambda over all points: phi is smooth, frozen gradient = exact gradient
    let o1 : EvalOpts := { nbhd := .all, lam := .constant 3.0, areas := k % 2 == 0 }
    match blend ps tori ar o1 x with
    | .error e => report tally s!"blend {k}" false e
    | .ok b =>
      let phi := fun y => match blend ps tori ar o1 y with | .ok r => r.phi | .error _ => 0.0
      let fg := fdGrad phi x 1e-6
      let fl := fdLap phi x 1e-4
      report tally s!"blend gradient (constant) {k}" (relV fg b.grad < 1e-6 && relV b.grad b.gradExact < 1e-15) s!"{relV fg b.grad}"
      report tally s!"blend laplacian (constant) {k}" (rel fl b.lap < 1e-4) s!"fd {fl} exact {b.lap}"
    -- per-query: the exact gradient matches differences while the neighbor set is locally constant
    let o2 : EvalOpts := { nbhd := .knn 8, lam := .perQuery 8.0 }
    match blend ps tori ar o2 x with
    | .error e => report tally s!"blend pq {k}" false e
    | .ok b =>
      let phi := fun y => match blend ps tori ar o2 y with | .ok r => r.phi | .error _ => 0.0
      let fg := fdGrad phi x 1e-7
      report tally s!"blend gradient (exact) {k}" (relV fg b.gradExact < 1e-5) s!"{relV fg b.gradExact} frozen {relV fg b.grad}"
    -- convexity
    match blend ps tori ar o2 x with
    | .error _ => pure ()
    | .ok b =>
      let gs := b.nset.idx.map fun i => (tori[i]!.sdf x).g
      let lo := gs.foldl min gs[0]!
      let hi := gs.foldl max gs[0]!
      report tally s!"blend convex {k}" (lo - 1e-12 ≤ b.phi && b.phi ≤ hi + 1e-12) s!"{lo} {b.phi} {hi}"
    -- baselines
    match baselines ps ns ar o1 x with
    | .error e => report tally s!"baselines {k}" false e
    | .ok b =>
      let shc := fun y => match baselines ps ns ar o1 y with | .ok r => r.shc | .error _ => 0.0
      let sspd := fun y => match baselines ps ns ar o1 y with | .ok r => r.sspd | .error _ => 0.0
      report tally s!"shc gradient {k}" (relV (fdGrad shc x 1e-6) b.shcGrad < 1e-6) s!"{relV (fdGrad shc x 1e-6) b.shcGrad}"
      report tally s!"sspd gradient {k}" (relV (fdGrad sspd x 1e-6) b.sspdGrad < 1e-6) s!"{relV (fdGrad sspd x 1e-6) b.sspdGrad}"
  -- the blend of identical tori is the torus
  do
    let t : Torus := { c := ⟨0.1, 0.2, 0.3⟩, u := ez, R := 1.0, r := 0.3 }
    let ps := #[(⟨0.0, 0.0, 0.0⟩ : V3), ⟨1.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩, ⟨0.5, 0.5, 0.5⟩]
    let x : V3 := ⟨0.7, -0.2, 0.4⟩
    match blend ps (ps.map fun _ => t) (ps.map fun _ => 1.0) {} x with
    | .ok b => report tally "blend identical" (rel b.phi (t.sdf x).g < 1e-14 && relV b.grad (t.sdf x).grad < 1e-13) ""
    | .error e => report tally "blend identical" false e
  -- print / parse round trips of generated cases
  for k in [0:6] do
    let act : G PatCase := do
      let (ps, ns) := fibonacci 10 1.0 V3.zero
      let mut c := cloudCase s!"rt{k}" "blend" ps ns
      for i in [0:10] do c := { c with coefs := c.coefs.push (i, ← randCoef) }
      let qs ← randQueries 5 (-1.5) 1.5
      return { c with queries := qs, opts := { nbhd := .radius none 4, lam := .paper 64.0 1000.0 5 } }
    let (c, _) := act.run (mkRng (5000 + k))
    let text := printCase c
    let lines := ((text.splitOn "\n").toArray.filterMap fun l =>
      let t := tokens l
      if t.isEmpty || t[0]! == "case" || t[0]! == "end" then none else some t)
    match parseBody c.name lines with
    | .error e => report tally s!"round trip {k}" false e
    | .ok c2 => report tally s!"round trip {k}" (printCase c2 == text && runCase c2 == runCase c) ""
  let t ← tally.get
  IO.println s!"pat selftest: {t.total} checks, {t.bad} failures"
  return t.bad

end Ddg.PointsAsTori
