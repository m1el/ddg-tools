import DdgOracle.RepulsiveCurves.Energy

/-!
# The Sobolev–Slobodeckij Gram matrices `B`, `B⁰`, `A = B + B⁰` (SPEC §4)

By their increments (§4.2, §4.3), for both quadratures of §4.4, and the kernel matrices `K^B`, `K⁰`.
`[formula: rc-high-order-gram]`, `[formula: rc-low-order-gram]`, `[formula: rc-kernel-matrix-form]`.
-/

namespace Ddg.RepulsiveCurves

/-- Quadrature of the weights (SPEC §4.4, D2). -/
inductive Quad where
  /-- one point at the edge midpoints (Eq. 32) -/
  | midpoint
  /-- the four endpoint pairs (Eq. 20-21) -/
  | trapezoid
  deriving BEq, Inhabited

/-- Name in the text format. -/
def Quad.name : Quad → String
  | .midpoint => "midpoint"
  | .trapezoid => "trapezoid"

/-- The weights `(w_IJ, w⁰_IJ)` of the ordered pair `(I, J)` (SPEC §4.2-4.3). -/
def pairWeights (net : Net) (g : Geo) (ex : Expo) (q : Quad) (I J : Nat) : Float × Float :=
  let sig := ex.sigma
  let lI := g.len[I]!
  let lJ := g.len[J]!
  let T := g.tan[I]!
  match q with
  | .midpoint =>
    let d := g.mid[I]! - g.mid[J]!
    let r := d.norm
    let w := lI * lJ * r.pow (-(2.0 * sig + 1.0))
    let w0 := lI * lJ * (V3.cross T d).norm2 / r.pow (2.0 * sig + 5.0)
    (w, w0)
  | .trapezoid => Id.run do
    let (i1, i2) := net.edges[I]!
    let (j1, j2) := net.edges[J]!
    let mut s : Float := 0.0
    let mut s0 : Float := 0.0
    for i in [i1, i2] do
      for j in [j1, j2] do
        let d := net.pos[i]! - net.pos[j]!
        let r := d.norm
        s := s + r.pow (-(2.0 * sig + 1.0))
        s0 := s0 + (V3.cross T d).norm2 / r.pow (2.0 * sig + 5.0)
    return (0.25 * lI * lJ * s, 0.25 * lI * lJ * s0)

/-- Dense `|V| × |V|` matrices, row-major. -/
structure Gram where
  /-- `B` -/
  b : FloatArray
  /-- `B⁰` -/
  b0 : FloatArray
  /-- `K^B_IJ = w_IJ + w_JI` (`|E| × |E|`) -/
  kb : FloatArray
  /-- `K⁰_IJ = w⁰_IJ + w⁰_JI` (`|E| × |E|`) -/
  k0 : FloatArray

/-- `B`, `B⁰` by the increments of SPEC §4.2-4.3 over all ordered pairs `(I, J) ∈ 𝒫`; `K^B`, `K⁰`. -/
def gram (net : Net) (ex : Expo) (q : Quad) : Gram := Id.run do
  let g := net.geo
  let n := net.nv
  let m := net.ne
  let mut b := zeros (n * n)
  let mut b0 := zeros (n * n)
  let mut kb := zeros (m * m)
  let mut k0 := zeros (m * m)
  for I in [0:m] do
    for J in [0:m] do
      if net.adj I J then continue
      let (w, w0) := pairWeights net g ex q I J
      let (wt, w0t) := pairWeights net g ex q J I
      kb := kb.set! (I * m + J) (w + wt)
      k0 := k0.set! (I * m + J) (w0 + w0t)
      let (i1, i2) := net.edges[I]!
      let (j1, j2) := net.edges[J]!
      let iv := #[i1, i2]
      let jv := #[j1, j2]
      let lI := g.len[I]!
      let lJ := g.len[J]!
      let tt := (g.tan[I]!).dot (g.tan[J]!)
      for a in [0:2] do
        for c in [0:2] do
          let sgn : Float := if (a + c) % 2 == 0 then 1.0 else -1.0
          let add := fun (arr : FloatArray) (r col : Nat) (v : Float) => arr.set! (r * n + col) (arr.get! (r * n + col) + v)
          b := add b iv[a]! iv[c]! (sgn * w / (lI * lI))
          b := add b iv[a]! jv[c]! (-(sgn * w * tt / (lI * lJ)))
          b := add b jv[a]! jv[c]! (sgn * w / (lJ * lJ))
          b := add b jv[a]! iv[c]! (-(sgn * w * tt / (lJ * lI)))
          b0 := add b0 iv[a]! iv[c]! (0.25 * w0)
          b0 := add b0 iv[a]! jv[c]! (-(0.25 * w0))
          b0 := add b0 jv[a]! iv[c]! (-(0.25 * w0))
          b0 := add b0 jv[a]! jv[c]! (0.25 * w0)
  return { b, b0, kb, k0 }

/-- `A = B + B⁰`. -/
def Gram.a (G : Gram) : FloatArray := Id.run do
  let mut a := G.b
  for i in [0:a.size] do
    a := a.set! i (a.get! i + G.b0.get! i)
  return a

end Ddg.RepulsiveCurves
