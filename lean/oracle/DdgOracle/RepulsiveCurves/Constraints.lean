import DdgOracle.RepulsiveCurves.Gram

/-!
# Constraints and their Jacobians (SPEC §5.1)

`[formula: rc-constraint-barycenter]`, `rc-constraint-total-length`, `rc-constraint-edge-length`,
`rc-constraint-point`, `rc-constraint-surface`, `rc-constraint-tangent`, `rc-tangent-frame`.
Targets default to the values at the positions the constraints are resolved at.
-/

namespace Ddg.RepulsiveCurves

/-- Where a tangent constraint applies (D7). -/
inductive Site where
  /-- the tangent of an edge -/
  | edge (e : Nat)
  /-- the tangent at a vertex of degree 1 or 2 -/
  | vertex (v : Nat)
  deriving Inhabited

/-- A constraint as given (`none`: the current value). -/
inductive Cons where
  /-- `Σ ℓ_I x_I / Σ ℓ_I − x₀` -/
  | barycenter (target : Option V3)
  /-- `L⁰ − Σ ℓ_I` -/
  | totalLength (target : Option Float)
  /-- `ℓ⁰_I − ℓ_I`; `none`: every edge -/
  | edgeLengths (edges : Option (Array Nat))
  /-- `γ_v − x_v` -/
  | point (v : Nat) (target : Option V3)
  /-- `|γ_v − c|² − R²` on the listed vertices (`none`: all) -/
  | sphere (c : V3) (r : Float) (vs : Option (Array Nat))
  /-- `(e₁·s, e₂·s)` for the tangent `s` of the site; `none`: `X = s` -/
  | tangent (site : Site) (dir : Option V3)

/-- A constraint with its targets filled in. -/
inductive RCons where
  /-- barycenter target -/
  | barycenter (x0 : V3)
  /-- total length target -/
  | totalLength (l0 : Float)
  /-- listed edges and their target lengths -/
  | edgeLengths (es : Array Nat) (ts : Array Float)
  /-- vertex and target -/
  | point (v : Nat) (x0 : V3)
  /-- sphere centre, radius, vertices -/
  | sphere (c : V3) (r : Float) (vs : Array Nat)
  /-- site, direction `X`, frame `e₁, e₂` -/
  | tangent (site : Site) (dir e1 e2 : V3)
  deriving Inhabited

/-- A scalar row `Φ_r` with its Jacobian entries `(vertex, coordinate, value)`. -/
structure Row where
  /-- `Φ_r` -/
  phi : Float
  /-- 0 length, 1 dimensionless (tangent), 2 surface (SPEC §5.1 table) -/
  unit : Nat
  /-- `|∇f|` over the active coordinates for surface rows, else 1 -/
  gscale : Float
  /-- index of the constraint -/
  cons : Nat
  /-- Jacobian entries (repeated positions add) -/
  ents : Array (Nat × Nat × Float)
  deriving Inhabited

/-- Deterministic orthonormal pair completing `X/|X|` (Duff et al. 2017; SPEC §5.1, D7). -/
def tangentFrame (x : V3) : V3 × V3 :=
  let n := (1.0 / x.norm) * x
  let neg := n.z < 0.0 || (n.z == 0.0 && 1.0 / n.z < 0.0)
  let sign : Float := if neg then -1.0 else 1.0
  let a := -1.0 / (sign + n.z)
  let b := n.x * n.y * a
  (⟨1.0 + sign * n.x * n.x * a, sign * b, -sign * n.x⟩, ⟨b, sign + n.y * n.y * a, -n.y⟩)

/-- The other end of edge `e` seen from `v`. -/
def otherEnd (net : Net) (e v : Nat) : Nat :=
  let (a, b) := net.edges[e]!
  if a == v then b else a

/-- `s` of a site and the pieces `(vertex, sign, ℓ, T)` of its derivative: `∂s/∂γ_vertex` has the
block `sign · P_T / ℓ` (`P_T = I − T Tᵀ`, symmetric). Degree 2 (D7), `a`, `b` the other ends of the
first and second incident edge: `s = (γ_v − γ_a)/|·| + (γ_b − γ_v)/|·|`. -/
def tangentOf (net : Net) (pos : Array V3) (site : Site) :
    Except String (V3 × Array (Nat × Float × Float × V3)) := do
  let unit := fun (p q : Nat) =>
    let d := pos[q]! - pos[p]!
    let l := d.norm
    ((1.0 / l) * d, l)
  match site with
  | .edge e =>
    if e ≥ net.ne then throw s!"ConstraintOutOfRange: edge {e}"
    let (a, b) := net.edges[e]!
    let (t, l) := unit a b
    return (t, #[(a, -1.0, l, t), (b, 1.0, l, t)])
  | .vertex v =>
    if v ≥ net.nv then throw s!"ConstraintOutOfRange: vertex {v}"
    let inc := net.incident v
    if inc.size == 1 then
      let w := otherEnd net inc[0]! v
      let (t, l) := unit v w
      return (t, #[(v, -1.0, l, t), (w, 1.0, l, t)])
    else if inc.size == 2 then
      let a := otherEnd net inc[0]! v
      let b := otherEnd net inc[1]! v
      let (u1, l1) := unit a v
      let (u2, l2) := unit v b
      return (u1 + u2, #[(a, -1.0, l1, u1), (v, 1.0, l1, u1), (v, -1.0, l2, u2), (b, 1.0, l2, u2)])
    else throw s!"InvalidTangentSite: vertex {v} has degree {inc.size}"

/-- Fill in the targets at `pos` (SPEC §5.1); errors named like the Rust variants. -/
def resolve (net : Net) (pos : Array V3) (planar : Option Nat) (cs : Array Cons) :
    Except String (Array RCons) := do
  let geo := edgeGeo net.edges pos
  let mut out : Array RCons := #[]
  for c in cs do
    match c with
    | .barycenter t =>
      let mut num := V3.zero
      let mut L : Float := 0.0
      for i in [0:net.ne] do
        num := num + (geo.len[i]! * geo.mid[i]!)
        L := L + geo.len[i]!
      out := out.push (RCons.barycenter (t.getD ((1.0 / L) * num)))
    | .totalLength t =>
      out := out.push (.totalLength (t.getD (geo.len.foldl (· + ·) 0.0)))
    | .edgeLengths es =>
      let es := es.getD (Array.range net.ne)
      for e in es do
        if e ≥ net.ne then throw s!"ConstraintOutOfRange: edge {e}"
      out := out.push (.edgeLengths es (es.map fun e => geo.len[e]!))
    | .point v t =>
      if v ≥ net.nv then throw s!"ConstraintOutOfRange: vertex {v}"
      out := out.push (.point v (t.getD pos[v]!))
    | .sphere c r vs =>
      let vs := vs.getD (Array.range net.nv)
      for v in vs do
        if v ≥ net.nv then throw s!"ConstraintOutOfRange: vertex {v}"
      out := out.push (.sphere c r vs)
    | .tangent site dir =>
      if planar.isSome then throw "Unsupported: tangent constraints in planar mode are not modelled by the oracle"
      let (s, _) ← tangentOf net pos site
      let x := dir.getD s
      if x.norm == 0.0 || !x.finite then throw "InvalidOption: zero or non-finite tangent direction"
      if !(s.dot x > 0.0) then throw "InfeasibleConstraint: tangent opposite to its direction"
      let (e1, e2) := tangentFrame x
      out := out.push (.tangent site x e1 e2)
  return out

/-- `Φ` and the Jacobian entries of every row at `pos`, in the order the constraints are listed
(`coords`: the active coordinates, SPEC §2.2). -/
def evalRows (net : Net) (pos : Array V3) (coords : Array Nat) (rc : Array RCons) : Array Row := Id.run do
  let geo := edgeGeo net.edges pos
  let mut rows : Array Row := #[]
  for ci in [0:rc.size] do
    match rc[ci]! with
    | .barycenter x0 =>
      let mut num := V3.zero
      let mut L : Float := 0.0
      for i in [0:net.ne] do
        num := num + (geo.len[i]! * geo.mid[i]!)
        L := L + geo.len[i]!
      let xbar := (1.0 / L) * num
      for r in coords do
        let mut ents : Array (Nat × Nat × Float) := #[]
        for i in [0:net.ne] do
          let (a, b) := net.edges[i]!
          let lever := geo.mid[i]! - xbar
          let T := geo.tan[i]!
          for (v, s) in [(a, (-1.0 : Float)), (b, 1.0)] do
            for c in [0:3] do
              let diag : Float := if c == r then 0.5 * geo.len[i]! else 0.0
              ents := ents.push (v, c, (diag + s * lever.get r * T.get c) / L)
        rows := rows.push { phi := xbar.get r - x0.get r, unit := 0, gscale := 1.0, cons := ci, ents }
    | .totalLength l0 =>
      let mut ents : Array (Nat × Nat × Float) := #[]
      for i in [0:net.ne] do
        let (a, b) := net.edges[i]!
        for (v, s) in [(a, (-1.0 : Float)), (b, 1.0)] do
          for c in [0:3] do
            ents := ents.push (v, c, -(s * (geo.tan[i]!).get c))
      rows := rows.push { phi := l0 - geo.len.foldl (· + ·) 0.0, unit := 0, gscale := 1.0, cons := ci, ents }
    | .edgeLengths es ts =>
      for k in [0:es.size] do
        let e := es[k]!
        let (a, b) := net.edges[e]!
        let mut ents : Array (Nat × Nat × Float) := #[]
        for c in [0:3] do
          ents := ents.push (a, c, (geo.tan[e]!).get c)
          ents := ents.push (b, c, -((geo.tan[e]!).get c))
        rows := rows.push { phi := ts[k]! - geo.len[e]!, unit := 0, gscale := 1.0, cons := ci, ents }
    | .point v x0 =>
      for r in coords do
        rows := rows.push { phi := (pos[v]!).get r - x0.get r, unit := 0, gscale := 1.0, cons := ci,
                            ents := #[(v, r, 1.0)] }
    | .sphere c rad vs =>
      for v in vs do
        let d := pos[v]! - c
        let grad := (2.0 : Float) * d
        let mut gs : Float := 0.0
        for r in coords do gs := gs + (grad.get r) * (grad.get r)
        rows := rows.push { phi := d.norm2 - rad * rad, unit := 2, gscale := Float.sqrt gs, cons := ci,
                            ents := #[(v, 0, grad.x), (v, 1, grad.y), (v, 2, grad.z)] }
    | .tangent site _ e1 e2 =>
      match tangentOf net pos site with
      | .error _ => pure ()
      | .ok (s, parts) =>
        for e in [e1, e2] do
          let mut ents : Array (Nat × Nat × Float) := #[]
          for (w, sign, l, T) in parts do
            let pe := e - (T.dot e) * T      -- P_T e, then the block is sign · P_T / ℓ
            for c in [0:3] do
              ents := ents.push (w, c, sign * pe.get c / l)
          rows := rows.push { phi := e.dot s, unit := 1, gscale := 1.0, cons := ci, ents }
  return rows

/-- The dense Jacobian `C` (`k × d|V|`, row-major; column `v · |coords| + slot`). -/
def denseJacobian (nv : Nat) (coords : Array Nat) (rows : Array Row) : FloatArray := Id.run do
  let nc := coords.size
  let nd := nv * nc
  let mut c := zeros (rows.size * nd)
  for k in [0:rows.size] do
    for (v, co, val) in rows[k]!.ents do
      match coords.findIdx? (· == co) with
      | some slot => c := c.set! (k * nd + v * nc + slot) (c.get! (k * nd + v * nc + slot) + val)
      | none => pure ()
  return c

end Ddg.RepulsiveCurves
