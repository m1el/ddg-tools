import DdgOracle.RepulsiveCurves.Basic

/-!
# The discrete tangent-point energy and its differential (SPEC §3)

`[formula: rc-kernel]`, `[formula: rc-discrete-kernel]`, `[formula: rc-energy]`,
`[formula: rc-kernel-gradient]`, `[formula: rc-energy-differential]`, written from the SPEC text.
-/

namespace Ddg.RepulsiveCurves

/-- `k(p, q, T) = |T × d|^α / |d|^β` for `d = p − q`, with `∇_p k` and `∇_T k` (SPEC §3.1–3.2).
`ν^{α−1} n̂` is `0` at `ν = 0`. -/
def kernelParts (ex : Expo) (T d : V3) : Float × V3 × V3 :=
  let rho := d.norm
  let k := (V3.cross T d).norm.pow ex.alpha / rho.pow ex.beta
  let td := T.dot d
  let n := d - td * T
  let nu := n.norm
  let w : V3 := if nu == 0.0 then V3.zero else (nu.pow (ex.alpha - 1.0) / nu) * n
  let invRhoBeta := 1.0 / rho.pow ex.beta
  let gp := (ex.alpha * invRhoBeta) * w - (ex.beta * k / (rho * rho)) * d
  let gT := (-(ex.alpha * td) * invRhoBeta) * w
  (k, gp, gT)

/-- `Ê` and the differential `∂Ê/∂γ_v` of every vertex. Sums over ordered pairs `(I, J)` of
edges without a common vertex. -/
def energyAndDifferential (net : Net) (ex : Expo) : Float × Array V3 := Id.run do
  let g := net.geo
  let m := net.ne
  let mut total : Float := 0.0
  let mut grad : Array V3 := Array.replicate net.nv V3.zero
  for I in [0:m] do
    for J in [0:m] do
      if net.adj I J then continue
      let (i1, i2) := net.edges[I]!
      let (j1, j2) := net.edges[J]!
      let ivs := #[i1, i2]
      let jvs := #[j1, j2]
      let T := g.tan[I]!
      let lI := g.len[I]!
      let lJ := g.len[J]!
      let mut khat : Float := 0.0
      let mut gT := V3.zero
      let mut gpI : Array V3 := #[V3.zero, V3.zero]   -- ¼ Σ_j ∇_p k(γ_{i_a}, γ_j, T_I)
      let mut gqJ : Array V3 := #[V3.zero, V3.zero]   -- ¼ Σ_i ∇_q k(γ_i, γ_{j_b}, T_I)
      for a in [0:2] do
        for b in [0:2] do
          let d := net.pos[ivs[a]!]! - net.pos[jvs[b]!]!
          let (k, gp, gt) := kernelParts ex T d
          khat := khat + 0.25 * k
          gT := gT + (0.25 : Float) * gt
          gpI := gpI.set! a (gpI[a]! + (0.25 : Float) * gp)
          gqJ := gqJ.set! b (gqJ[b]! - (0.25 : Float) * gp)
      total := total + khat * lI * lJ
      let TJ := g.tan[J]!
      for a in [0:2] do
        let s : Float := if a == 0 then -1.0 else 1.0
        let c1 := (s * lJ) * ((khat * T) + gT)
        let c2 := (lI * lJ) * gpI[a]!
        grad := grad.modify ivs[a]! (· + c1 + c2)
      for b in [0:2] do
        let s : Float := if b == 0 then -1.0 else 1.0
        let c1 := (s * lI * khat) * TJ
        let c2 := (lI * lJ) * gqJ[b]!
        grad := grad.modify jvs[b]! (· + c1 + c2)
  return (total, grad)

end Ddg.RepulsiveCurves
