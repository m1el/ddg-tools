import DdgOracle.RepulsiveCurves.Energy

/-!
# Near-pair quadrature (SPEC §3.4, WP16; not in the paper)

`[formula: rc-near-h-constant]`, `[formula: rc-near-inner]`, `[formula: rc-near-cos-series]`,
`[formula: rc-near-outer]`, `[formula: rc-near-blend]`, written from the SPEC text (§3.4), not
from the crate's code: the inner integral `F(y) = h^{1+α−β}(S(θ₁) − S(θ₀))` with the two
binomial series of `S`, the graded outer quadrature over `J` (breakpoints `t*` and the
projections of `I`'s ends, sinh grading with `δ = d/ℓ_J`, `M` sub-panels of 8-point
Gauss–Legendre per half-panel), and the C² blend `e = e_p + w(d/ℓ_IJ)(e_n − e_p)`. All edge pairs
are visited (no broad phase). The differential is not derived here: `nearEnergyAndDifferential`
adds the paper's analytic differential (`Energy.lean`) and a Richardson-extrapolated central
difference of the correction (an independent check of the crate's analytic differential).
-/

namespace Ddg.RepulsiveCurves

/-- The 8-point Gauss–Legendre rule on `[−1, 1]`. -/
def gl8 : Array (Float × Float) := #[
  (-0.9602898564975362, 0.10122853629037669),
  (-0.7966664774136267, 0.22238103445337434),
  (-0.525532409916329, 0.31370664587788705),
  (-0.18343464249564978, 0.36268378337836177),
  (0.18343464249564978, 0.36268378337836177),
  (0.525532409916329, 0.31370664587788705),
  (0.7966664774136267, 0.22238103445337434),
  (0.9602898564975362, 0.10122853629037669)]

/-- Options of the near-pair quadrature: `c₁`, `c₂`, `M`. -/
structure NearOpts where
  c1 : Float := 0.25
  c2 : Float := 0.75
  m : Nat := 4
  deriving Inhabited

/-- Series tables for `m = β − 2`. -/
structure Series where
  m : Float
  a : Array Float
  r : Array Float
  half : Float

private def nTerms : Nat := 60

private def evalLow (t : Series) (s : Float) : Float := Id.run do
  -- Σ_k C(a,k)(−1)^k s^{2k+1}/(2k+1)
  let mut acc := 0.0
  let mut p := s
  for k in [0:nTerms] do
    acc := acc + t.a[k]! * p
    p := p * s * s
  return acc

private def evalHigh (t : Series) (c : Float) : Float := Id.run do
  -- c^{m+1} Σ_k C(−½,k)(−1)^k c^{2k}/(m+1+2k)
  let mut acc := 0.0
  let mut p := 1.0
  for k in [0:nTerms] do
    acc := acc + t.r[k]! * p
    p := p * c * c
  return acc * Float.pow c (t.m + 1.0)

/-- The tables (SPEC §3.4). -/
def Series.make (ex : Expo) : Series := Id.run do
  let m := ex.beta - 2.0
  let a0 := 0.5 * (m - 1.0)
  let mut ca := 1.0
  let mut cr := 1.0
  let mut aa : Array Float := #[]
  let mut rr : Array Float := #[]
  for k in [0:nTerms] do
    let kf := k.toFloat
    aa := aa.push (ca / (2.0 * kf + 1.0))
    rr := rr.push (cr / (m + 1.0 + 2.0 * kf))
    ca := ca * (-(a0 - kf) / (kf + 1.0))
    cr := cr * ((kf + 0.5) / (kf + 1.0))
  let t : Series := { m, a := aa, r := rr, half := 0.0 }
  let c := Float.sqrt 0.5
  return { t with half := evalLow t c + evalHigh t c }

/-- `S(θ) = ∫_0^θ cos^m` from `sin θ`, `cos θ ≥ 0`. -/
def Series.S (t : Series) (sn cs : Float) : Float :=
  if sn.abs ≤ cs then evalLow t sn
  else if sn > 0.0 then t.half - evalHigh t cs
  else evalHigh t cs - t.half

/-- `F(y) = ∫_I k(x, y, T_I) dx` for `I = (a, a1)`. -/
def innerF (t : Series) (ex : Expo) (y a a1 : V3) : Float := Id.run do
  let e := a1 - a
  let l := e.norm
  let T := (1.0 / l) * e
  let r := y - a
  let u0 := T.dot r
  let h := (r - u0 * T).norm
  if h == 0.0 then return 0.0
  let x1 := l - u0
  let x0 := -u0
  let r1 := Float.sqrt (x1 * x1 + h * h)
  let r0 := Float.sqrt (x0 * x0 + h * h)
  let c :=
    if x0 ≥ h then evalHigh t (h / r0) - evalHigh t (h / r1)
    else if x1 ≤ -h then evalHigh t (h / r1) - evalHigh t (h / r0)
    else t.S (x1 / r1) (h / r1) - t.S (x0 / r0) (h / r0)
  return Float.pow h (1.0 + ex.alpha - ex.beta) * c

/-- Closest points `(d, s, t)` of the segments `p0p1`, `q0q1` (Ericson §5.1.9; parallel: `s = 0`). -/
def closestST (p0 p1 q0 q1 : V3) : Float × Float × Float := Id.run do
  let A := p1 - p0
  let B := q1 - q0
  let r := p0 - q0
  let aa := A.dot A
  let ee := B.dot B
  let f := B.dot r
  let c := A.dot r
  let b := A.dot B
  let den := aa * ee - b * b
  let clamp := fun (x : Float) => if x < 0.0 then 0.0 else if x > 1.0 then 1.0 else x
  let mut s := if den > 1e-14 * aa * ee then clamp ((b * f - c * ee) / den) else 0.0
  let mut t := (b * s + f) / ee
  if t ≤ 0.0 || t ≥ 1.0 then
    t := clamp t
    s := clamp ((b * t - c) / aa)
  let d := (r + s * A - t * B).norm
  return (d, s, t)

private def projParam (x q0 q1 : V3) : Float :=
  let B := q1 - q0
  let q := (x - q0).dot B / B.dot B
  if q < 0.0 then 0.0 else if q > 1.0 then 1.0 else q

/-- `e^near_IJ`: the graded quadrature of `∫_I ∫_J k(x, y, T_I)` (SPEC §3.4). `tStar` and `d`
are the closest-point parameter on `J` and the segment distance. -/
def nearOrdered (t : Series) (ex : Expo) (o : NearOpts) (p0 p1 q0 q1 : V3) (tStar d : Float) :
    Float := Id.run do
  let B := q1 - q0
  let lj := B.norm
  let delta := d / lj
  let bps := #[tStar, projParam p0 q0 q1, projParam p1 q0 q1].qsort (· < ·)
  let zs := #[0.0, bps[0]!, bps[1]!, bps[2]!, 1.0]
  let mut sum := 0.0
  let mf := o.m.toFloat
  for k in [0:4] do
    let lo := zs[k]!
    let hi := zs[k+1]!
    let half := 0.5 * (hi - lo)
    if half ≤ 0.0 then continue
    -- asinh(x) = ln(x + √(x² + 1))
    let x := half / delta
    let wmax := Float.log (x + Float.sqrt (x * x + 1.0))
    for (c, sg) in #[(lo, 1.0), (hi, -1.0)] do
      for j in [0:o.m] do
        for (xi, wk) in gl8 do
          let eta := (j.toFloat + 0.5 * (1.0 + xi)) / mf
          let w := wmax * eta
          let sh := 0.5 * (Float.exp w - Float.exp (-w))
          let ch := 0.5 * (Float.exp w + Float.exp (-w))
          let tt := c + sg * delta * sh
          let wt := wk / (2.0 * mf) * delta * wmax * ch
          sum := sum + wt * innerF t ex (q0 + tt * B) p0 p1
  return lj * sum

/-- The paper's ordered term `k̂_IJ ℓ_I ℓ_J` (Eq. 17–18). -/
def paperOrdered (ex : Expo) (p0 p1 q0 q1 : V3) : Float := Id.run do
  let li := (p1 - p0).norm
  let lj := (q1 - q0).norm
  let T := (1.0 / li) * (p1 - p0)
  let mut k := 0.0
  for a in #[p0, p1] do
    for b in #[q0, q1] do
      let (v, _, _) := kernelParts ex T (a - b)
      k := k + v
  return 0.25 * k * li * lj

/-- `w(r)`: 1 for `r ≤ c₁`, 0 for `r ≥ c₂`, `τ³(10 − 15τ + 6τ²)` between. -/
def blendW (o : NearOpts) (r : Float) : Float :=
  let tau := (o.c2 - r) / (o.c2 - o.c1)
  if tau ≥ 1.0 then 1.0 else if tau ≤ 0.0 then 0.0
  else tau * tau * tau * (10.0 - 15.0 * tau + 6.0 * tau * tau)

/-- Some edge joins an endpoint of `I` to an endpoint of `J` (second neighbours along the
network; excluded from the correction, SPEC §3.4). -/
def joined (net : Net) (I J : Nat) : Bool :=
  let (a, b) := net.edges[I]!
  let (c, d) := net.edges[J]!
  net.edges.any fun (x, y) =>
    ((x == a || x == b) && (y == c || y == d)) || ((y == a || y == b) && (x == c || x == d))

/-- One near pair of the output. -/
structure NearPairRec where
  i : Nat
  j : Nat
  d : Float
  w : Float
  corr : Float

/-- The correction `Σ_{I<J} w (Σ_orders e_n − e_p)` at `pos`, with the pairs in the band. -/
def nearCorrection (net : Net) (pos : Array V3) (ex : Expo) (o : NearOpts) :
    Float × Array NearPairRec := Id.run do
  let t := Series.make ex
  let mut total := 0.0
  let mut recs : Array NearPairRec := #[]
  for I in [0:net.ne] do
    for J in [I+1:net.ne] do
      if net.adj I J || joined net I J then continue
      let (i1, i2) := net.edges[I]!
      let (j1, j2) := net.edges[J]!
      let p0 := pos[i1]!
      let p1 := pos[i2]!
      let q0 := pos[j1]!
      let q1 := pos[j2]!
      let li := (p1 - p0).norm
      let lj := (q1 - q0).norm
      let (d, s, tt) := closestST p0 p1 q0 q1
      let lm := Float.sqrt (0.5 * (li * li + lj * lj))
      let r := d / lm
      if r ≥ o.c2 then continue
      let w := blendW o r
      let eIJ := nearOrdered t ex o p0 p1 q0 q1 tt d - paperOrdered ex p0 p1 q0 q1
      let eJI := nearOrdered t ex o q0 q1 p0 p1 s d - paperOrdered ex q0 q1 p0 p1
      let c := w * (eIJ + eJI)
      total := total + c
      recs := recs.push { i := I, j := J, d, w, corr := c }
  return (total, recs)

/-- Energy with the near-pair quadrature and its differential: the paper's analytic differential
plus central differences of the correction (steps `h` and `h/2`, Richardson),
`h = 1e-3 min(ℓ̄, min d_IJ)` over the near pairs. -/
def nearEnergyAndDifferential (net : Net) (ex : Expo) (o : NearOpts) :
    Float × Float × Array NearPairRec × Array V3 := Id.run do
  let (e, d) := energyAndDifferential net ex
  let (c, recs) := nearCorrection net net.pos ex o
  let mut grad := d
  if recs.size > 0 then
    let lbar := (net.geo.len.foldl (· + ·) 0.0) / net.ne.toFloat
    let h := 1e-3 * recs.foldl (fun acc r => min acc r.d) lbar
    -- vertices of near pairs only (the correction does not depend on the others)
    let mut verts : Array Nat := #[]
    for rec in recs do
      let (a, b) := net.edges[rec.i]!
      let (c', d') := net.edges[rec.j]!
      for v in #[a, b, c', d'] do
        if !verts.contains v then verts := verts.push v
    for v in verts do
      let mut g := V3.zero
      for ax in [0:3] do
        let f := fun (step : Float) =>
          (nearCorrection net (net.pos.set! v ((net.pos[v]!).set ax ((net.pos[v]!).get ax + step))) ex o).1
        let d1 := (f h - f (-h)) / (2.0 * h)
        let d2 := (f (0.5 * h) - f (-0.5 * h)) / h
        g := g.set ax ((4.0 * d2 - d1) / 3.0)
      grad := grad.set! v (grad[v]! + g)
  return (e, e + c, recs, grad)

end Ddg.RepulsiveCurves
