import DdgOracle.RepulsiveSurfaces.Metric

/-!
# Constraints: values and Jacobians (SPEC §5.1)

Barycenter (per component, frozen area weights), pins, total area and volume. Unknowns are ordered
`3 v + axis`.
-/

namespace Ddg.RepulsiveSurfaces

open Ddg.RepulsiveCurves (V3)

/-- A constraint as the case states it (targets default to the initial value). -/
inductive Cons where
  /-- barycenter of component `c` (`none`: every component separately) -/
  | barycenter (comp : Option Nat) (target : Option V3)
  /-- pinned vertices -/
  | pins (vs : Array Nat) (targets : Option (Array V3))
  /-- total area of a component (`none`: the whole mesh) -/
  | area (comp : Option Nat) (target : Option Float)
  /-- signed volume of a component (`none`: the whole mesh) -/
  | volume (comp : Option Nat) (target : Option Float)

/-- A constraint with its targets resolved. A barycenter is one component. -/
inductive RCons where
  /-- barycenter of component, target -/
  | bary (comp : Nat) (target : V3)
  /-- pins and targets -/
  | pins (vs : Array Nat) (targets : Array V3)
  /-- area of component (or all), target -/
  | area (comp : Option Nat) (target : Float)
  /-- volume of component (or all), target -/
  | volume (comp : Option Nat) (target : Float)
  deriving Inhabited

/-- Positional constraints live in `A′`, the others go through the Schur complement. -/
def RCons.positional : RCons → Bool
  | .bary .. => true
  | .pins .. => true
  | _ => false

/-- One scalar row: value and sparse Jacobian (column `3 v + axis`). -/
structure CRow where
  /-- `Φ` -/
  phi : Float
  /-- Jacobian entries -/
  jac : Array (Nat × Float)
  deriving Inhabited

/-- Sum of face areas of a component (or all). -/
def totalArea (t : Topo) (g : Array FaceGeo) (comp : Option Nat) : Float := Id.run do
  let mut s := 0.0
  for fi in [0:t.faces.size] do
    if comp.all (· == t.faceComp[fi]!) then s := s + g[fi]!.area
  return s

/-- Signed volume `⅙ Σ f_i·(f_j × f_k)` of a component (or all). -/
def totalVolume (t : Topo) (pos : Array V3) (comp : Option Nat) : Float := Id.run do
  let mut s := 0.0
  for fi in [0:t.faces.size] do
    if comp.all (· == t.faceComp[fi]!) then
      let f := t.faces[fi]!
      s := s + pos[f.a]!.dot (pos[f.b]!.cross pos[f.c]!) / 6.0
  return s

/-- Barycenter `Σ a_v f_v / Σ a_v` of a component. -/
def barycenterOf (t : Topo) (pos : Array V3) (comp : Nat) : V3 := Id.run do
  let g := faceGeo t.faces pos
  let av := vertexAreas t g
  let mut s := V3.zero
  let mut tot := 0.0
  for v in [0:t.nv] do
    if t.vertComp[v]! == comp then
      s := s + av[v]! * pos[v]!
      tot := tot + av[v]!
  return (1.0 / tot) * s

/-- The volume constraint needs consistently oriented components (SPEC §2.1). -/
def checkVolumeComps (t : Topo) (comp : Option Nat) : Except String Unit := do
  for c in [0:t.ncomp] do
    if comp.all (· == c) then
      match inconsistentEdge t c with
      | some (a, b, f1, f2) =>
        throw s!"NonOrientable: component {c} is not consistently oriented: edge ({a}, {b}) is traversed in the same direction by faces {f1} and {f2}"
      | none => pure ()

/-- Resolve defaults and component ranges at the base positions. -/
def resolveCons (t : Topo) (pos : Array V3) (cons : Array Cons) : Except String (Array RCons) := do
  let g := faceGeo t.faces pos
  let mut out : Array RCons := #[]
  for c in cons do
    match c with
    | .barycenter comp target =>
      match comp with
      | some k =>
        if k ≥ t.ncomp then throw s!"ConstraintOutOfRange: component {k}, {t.ncomp} components"
        out := out.push (.bary k (target.getD (barycenterOf t pos k)))
      | none =>
        for k in [0:t.ncomp] do
          out := out.push (.bary k (target.getD (barycenterOf t pos k)))
    | .pins vs targets =>
      for v in vs do
        if v ≥ t.nv then throw s!"ConstraintOutOfRange: pinned vertex {v}, {t.nv} vertices"
      match targets with
      | some tg =>
        if tg.size != vs.size then throw "ConstraintOutOfRange: pins and targets differ in number"
        out := out.push (.pins vs tg)
      | none => out := out.push (.pins vs (vs.map (pos[·]!)))
    | .area comp target =>
      if let some k := comp then
        if k ≥ t.ncomp then throw s!"ConstraintOutOfRange: component {k}, {t.ncomp} components"
      out := out.push (.area comp (target.getD (totalArea t g comp)))
    | .volume comp target =>
      if let some k := comp then
        if k ≥ t.ncomp then throw s!"ConstraintOutOfRange: component {k}, {t.ncomp} components"
      checkVolumeComps t comp
      out := out.push (.volume comp (target.getD (totalVolume t pos comp)))
  return out

/-- Rows (values and Jacobians) of a resolved constraint at `pos`. A barycenter has 3 rows (frozen
weights `a_v/A_c` of `pos`), a pin 3 rows per vertex, area and volume one row. -/
def evalCons (t : Topo) (pos : Array V3) (c : RCons) : Array CRow := Id.run do
  let g := faceGeo t.faces pos
  match c with
  | .bary comp target =>
    let av := vertexAreas t g
    let mut tot := 0.0
    for v in [0:t.nv] do
      if t.vertComp[v]! == comp then tot := tot + av[v]!
    let b := barycenterOf t pos comp
    let mut rows : Array CRow := #[]
    for ax in [0:3] do
      let mut jac : Array (Nat × Float) := #[]
      for v in [0:t.nv] do
        if t.vertComp[v]! == comp then jac := jac.push (3 * v + ax, av[v]! / tot)
      rows := rows.push { phi := b.get ax - target.get ax, jac }
    return rows
  | .pins vs targets =>
    let mut rows : Array CRow := #[]
    for k in [0:vs.size] do
      for ax in [0:3] do
        rows := rows.push { phi := (pos[vs[k]!]!).get ax - (targets[k]!).get ax, jac := #[(3 * vs[k]! + ax, 1.0)] }
    return rows
  | .area comp target =>
    let mut gr : Array V3 := Array.replicate t.nv V3.zero
    for fi in [0:t.faces.size] do
      if comp.all (· == t.faceComp[fi]!) then
        let f := t.faces[fi]!
        for k in [0:3] do
          gr := gr.modify (f.get k) (· + (0.5 : Float) * g[fi]!.n.cross (oppEdge pos f k))
    let mut jac : Array (Nat × Float) := #[]
    for v in [0:t.nv] do
      for ax in [0:3] do
        let x := (gr[v]!).get ax
        if x != 0.0 then jac := jac.push (3 * v + ax, x)
    return #[{ phi := totalArea t g comp - target, jac }]
  | .volume comp target =>
    let mut gr : Array V3 := Array.replicate t.nv V3.zero
    for fi in [0:t.faces.size] do
      if comp.all (· == t.faceComp[fi]!) then
        let f := t.faces[fi]!
        for k in [0:3] do
          let pj := pos[f.get ((k + 1) % 3)]!
          let pk := pos[f.get ((k + 2) % 3)]!
          gr := gr.modify (f.get k) (· + (1.0 / 6.0 : Float) * pj.cross pk)
    let mut jac : Array (Nat × Float) := #[]
    for v in [0:t.nv] do
      for ax in [0:3] do
        let x := (gr[v]!).get ax
        if x != 0.0 then jac := jac.push (3 * v + ax, x)
    return #[{ phi := totalVolume t pos comp - target, jac }]

/-- Row kind names of a resolved constraint. -/
def RCons.name : RCons → String
  | .bary .. => "barycenter"
  | .pins .. => "pins"
  | .area .. => "area"
  | .volume .. => "volume"

/-- All rows of the resolved constraints, tagged `(kind, positional, row)`, in order. -/
def evalAll (t : Topo) (pos : Array V3) (rc : Array RCons) : Array (String × Bool × CRow) := Id.run do
  let mut out := #[]
  for c in rc do
    for r in evalCons t pos c do out := out.push (c.name, c.positional, r)
  return out

/-- Every component must be pinned by a barycenter row or at least one pin (SPEC §5.1). -/
def checkTranslation (t : Topo) (rc : Array RCons) : Except String Unit := do
  for c in [0:t.ncomp] do
    let mut fixed := false
    for k in rc do
      match k with
      | .bary comp _ => if comp == c then fixed := true
      | .pins vs _ => for v in vs do if t.vertComp[v]! == c then fixed := true
      | _ => pure ()
    if !fixed then throw s!"TranslationNotFixed: component {c}"

end Ddg.RepulsiveSurfaces
