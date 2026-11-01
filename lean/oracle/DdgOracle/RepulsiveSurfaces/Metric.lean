import DdgOracle.RepulsiveSurfaces.Energy

/-!
# Dense operators: `L_τ`, `B`, `B⁰`, `A`, `W`, `M` (SPEC §4.2, §4.4)

All assembled densely from the ordered face pairs `S ≠ T`; vertex functions enter through the
sparse rows `ū_S = (u_i+u_j+u_k)/3` and `D_f u(S) = Σ u_t g_t(S)`.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

/-- Dense square matrix, row-major. -/
abbrev Mat := FloatArray

/-- `a * b` for `n × n` matrices. -/
def matMul (n : Nat) (a b : Mat) : Mat := Id.run do
  let mut c := zeros (n * n)
  for i in [0:n] do
    for k in [0:n] do
      let x := a.get! (i * n + k)
      if x != 0.0 then
        for j in [0:n] do
          c := c.set! (i * n + j) (c.get! (i * n + j) + x * b.get! (k * n + j))
  return c

/-- Cotan stiffness `W` (PSD) and lumped mass `a_v` (SPEC §2.1). -/
def stiffness (t : Topo) (pos : Array V3) : Mat := Id.run do
  let n := t.nv
  let mut w := zeros (n * n)
  for f in t.faces do
    for k in [0:3] do
      let p := f.get k
      let b := f.get ((k + 1) % 3)
      let c := f.get ((k + 2) % 3)
      let u := pos[b]! - pos[p]!
      let v := pos[c]! - pos[p]!
      let cot := u.dot v / (u.cross v).norm
      w := w.set! (b * n + c) (w.get! (b * n + c) - 0.5 * cot)
      w := w.set! (c * n + b) (w.get! (c * n + b) - 0.5 * cot)
      w := w.set! (b * n + b) (w.get! (b * n + b) + 0.5 * cot)
      w := w.set! (c * n + c) (w.get! (c * n + c) + 0.5 * cot)
  return w

/-- `W M⁻¹ W`. -/
def willmoreMat (t : Topo) (w : Mat) (m : Array Float) : Mat := Id.run do
  let n := t.nv
  let mut mw := zeros (n * n)
  for i in [0:n] do
    for j in [0:n] do
      mw := mw.set! (i * n + j) (w.get! (i * n + j) / m[i]!)
  return matMul n w mw

/-- Sparse row entries of one face for the pair form: `(vertex, scalar coefficient, vector coefficient)`. -/
private def faceEntries (pos : Array V3) (f : Tri) (g : FaceGeo) (sign : Float) : Array (Nat × Float × V3) :=
  (#[0, 1, 2] : Array Nat).map fun k => (f.get k, sign / 3.0, sign * dfCoef pos f g k)

/-- Merge entries with the same vertex. -/
private def merge (es : Array (Nat × Float × V3)) : Array (Nat × Float × V3) := Id.run do
  let mut out : Array (Nat × Float × V3) := #[]
  for (v, c, d) in es do
    match out.findIdx? (·.1 == v) with
    | some i => out := out.modify i fun (v0, c0, d0) => (v0, c0 + c, d0 + d)
    | none => out := out.push (v, c, d)
  return out

/-- The three pair-sum matrices of one pass: `L_τ`, `B`, `B⁰`. -/
structure Pairs where
  /-- `L_τ` -/
  l : Mat
  /-- `B` -/
  b : Mat
  /-- `B⁰` -/
  b0 : Mat

/-- Assemble `L_τ` (kernel `ρ^{−(2τ+2)}`), `B` and `B⁰` (kernel `ρ^{−2s}`) over all ordered face pairs.
`B⁰` uses `K₂(S,T) = |N_S·(X_S − X_T)|² ρ^{−4}` whatever the energy's exponents. -/
def assemblePairs (t : Topo) (pos : Array V3) (tau s : Float) : Except String Pairs := do
  let n := t.nv
  let nf := t.faces.size
  let g := faceGeo t.faces pos
  let mut l := zeros (n * n)
  let mut b := zeros (n * n)
  let mut b0 := zeros (n * n)
  for si in [0:nf] do
    for ti in [0:nf] do
      if si != ti then
        let d := g[si]!.x - g[ti]!.x
        let rho2 := d.norm2
        if rho2 == 0.0 then throw s!"NonFiniteEnergy: faces ({si}, {ti}) have the same barycenter"
        let rho := Float.sqrt rho2
        let aa := g[si]!.area * g[ti]!.area
        let wl := aa * Float.pow rho (-(2.0 * tau + 2.0))
        let wb := aa * Float.pow rho (-2.0 * s)
        let nu := g[si]!.n.dot d
        let k2 := nu * nu / (rho2 * rho2)
        let es := merge (faceEntries pos t.faces[si]! g[si]! 1.0 ++ faceEntries pos t.faces[ti]! g[ti]! (-1.0))
        for (p, cp, dp) in es do
          for (q, cq, dq) in es do
            let i := p * n + q
            l := l.set! i (l.get! i + wl * cp * cq)
            b := b.set! i (b.get! i + wb * dp.dot dq)
            b0 := b0.set! i (b0.get! i + wb * k2 * cp * cq)
  return { l, b, b0 }

/-- Options of `A` (SPEC §4.4). -/
structure MetricOpts where
  /-- tangent-point weight -/
  wtp : Float := 1.0
  /-- `c_lo`, the weight of `B⁰` (0 or 1) -/
  clo : Float := 1.0
  /-- Willmore weight -/
  ww : Float := 0.0

/-- `A = w_TP (B + c_lo B⁰) + w_W W M⁻¹ W`. -/
def metricA (n : Nat) (b b0 wmw : Mat) (o : MetricOpts) : Mat := Id.run do
  let mut a := zeros (n * n)
  for i in [0:n * n] do
    a := a.set! i (o.wtp * (b.get! i + o.clo * b0.get! i) + o.ww * wmw.get! i)
  return a

/-- Everything the `metric`, `saddle`, `precond` and `step` tasks need at a position. -/
structure Ops where
  /-- stiffness -/
  w : Mat
  /-- vertex areas -/
  m : Array Float
  /-- `B` -/
  b : Mat
  /-- `B⁰` -/
  b0 : Mat
  /-- `W M⁻¹ W` -/
  wmw : Mat
  /-- `A` -/
  a : Mat

/-- Operators of the metric at `pos`. -/
def buildOps (t : Topo) (pos : Array V3) (ex : Ex) (o : MetricOpts) : Except String Ops := do
  let p ← assemblePairs t pos 0.5 ex.s
  let w := stiffness t pos
  let m := vertexAreas t (faceGeo t.faces pos)
  let wmw := willmoreMat t w m
  return { w, m, b := p.b, b0 := p.b0, wmw, a := metricA t.nv p.b p.b0 wmw o }

end Ddg.RepulsiveSurfaces
