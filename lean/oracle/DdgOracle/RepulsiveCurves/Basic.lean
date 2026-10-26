import DdgOracle.Dense
import DdgOracle.Hex

/-!
# Repulsive curves: vectors, curve networks, validation

Plain `Float` code written from `crates/ddg-repulsive-curves/SPEC.md` §2 for the `rc` oracle. It
shares nothing with the Rust crate (and does not use the scalar-generic formulas of the other
papers).
-/

namespace Ddg.RepulsiveCurves

/-- A vector of `R³`. -/
structure V3 where
  /-- first coordinate -/
  x : Float
  /-- second coordinate -/
  y : Float
  /-- third coordinate -/
  z : Float
  deriving Inhabited

namespace V3

/-- zero vector -/
def zero : V3 := ⟨0.0, 0.0, 0.0⟩
/-- sum -/
def add (a b : V3) : V3 := ⟨a.x + b.x, a.y + b.y, a.z + b.z⟩
/-- difference -/
def sub (a b : V3) : V3 := ⟨a.x - b.x, a.y - b.y, a.z - b.z⟩
/-- negation -/
def neg (a : V3) : V3 := ⟨-a.x, -a.y, -a.z⟩
/-- scalar multiple -/
def smul (s : Float) (a : V3) : V3 := ⟨s * a.x, s * a.y, s * a.z⟩
/-- dot product -/
def dot (a b : V3) : Float := a.x * b.x + a.y * b.y + a.z * b.z
/-- cross product -/
def cross (a b : V3) : V3 :=
  ⟨a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x⟩
/-- squared norm -/
def norm2 (a : V3) : Float := dot a a
/-- Euclidean norm -/
def norm (a : V3) : Float := Float.sqrt (norm2 a)
/-- all coordinates finite -/
def finite (a : V3) : Bool := a.x.isFinite && a.y.isFinite && a.z.isFinite
/-- coordinate `i ∈ {0,1,2}` -/
def get (a : V3) (i : Nat) : Float := if i == 0 then a.x else if i == 1 then a.y else a.z
/-- `a` with coordinate `i` replaced -/
def set (a : V3) (i : Nat) (v : Float) : V3 :=
  if i == 0 then ⟨v, a.y, a.z⟩ else if i == 1 then ⟨a.x, v, a.z⟩ else ⟨a.x, a.y, v⟩

instance : Add V3 := ⟨add⟩
instance : Sub V3 := ⟨sub⟩
instance : Neg V3 := ⟨neg⟩
instance : HMul Float V3 V3 := ⟨smul⟩

end V3

/-- A curve network: positions and edges `(i₁, i₂)` (SPEC §2.1). -/
structure Net where
  /-- vertex positions -/
  pos : Array V3
  /-- oriented edges -/
  edges : Array (Nat × Nat)
  deriving Inhabited

/-- The exponent pair `(α, β)` of the energy (SPEC §1). -/
structure Expo where
  /-- `α` -/
  alpha : Float
  /-- `β` -/
  beta : Float
  deriving Inhabited

/-- `σ = s − 1 = (β − 1)/α − 1` (SPEC §1). -/
def Expo.sigma (e : Expo) : Float := (e.beta - 1.0) / e.alpha - 1.0

/-- `α > 1`, `α + 2 ≤ β < 2α + 1`, both finite (SPEC §2.3, D24). -/
def Expo.check (e : Expo) : Except String Unit := do
  if !(e.alpha.isFinite && e.beta.isFinite) then throw "InvalidExponents: both exponents must be finite"
  if !(e.alpha > 1.0) then throw "InvalidExponents: alpha must be greater than 1"
  if e.beta < e.alpha + 2.0 then throw "InvalidExponents: beta must be at least alpha + 2"
  if !(e.beta < 2.0 * e.alpha + 1.0) then throw "InvalidExponents: beta must be smaller than 2 alpha + 1"

/-- `ℓ_I`, `T_I`, `x_I` of every edge (SPEC §2.1). -/
structure Geo where
  /-- lengths -/
  len : Array Float
  /-- unit tangents -/
  tan : Array V3
  /-- midpoints -/
  mid : Array V3

/-- Edge quantities at positions `pos`. -/
def edgeGeo (edges : Array (Nat × Nat)) (pos : Array V3) : Geo := Id.run do
  let mut len : Array Float := #[]
  let mut tan : Array V3 := #[]
  let mut mid : Array V3 := #[]
  for (a, b) in edges do
    let d := pos[b]! - pos[a]!
    let l := d.norm
    len := len.push l
    tan := tan.push ((1.0 / l) * d)
    mid := mid.push ((0.5 : Float) * (pos[a]! + pos[b]!))
  return { len, tan, mid }

/-- The geometry of the network. -/
def Net.geo (n : Net) : Geo := edgeGeo n.edges n.pos

/-- Number of vertices. -/
def Net.nv (n : Net) : Nat := n.pos.size
/-- Number of edges. -/
def Net.ne (n : Net) : Nat := n.edges.size

/-- Do edges `i` and `j` share a vertex (also `i = j`)? -/
def Net.adj (n : Net) (i j : Nat) : Bool :=
  let (a, b) := n.edges[i]!
  let (c, d) := n.edges[j]!
  a == c || a == d || b == c || b == d

/-- Edges incident to `v`, in increasing edge index. -/
def Net.incident (n : Net) (v : Nat) : Array Nat := Id.run do
  let mut r : Array Nat := #[]
  for i in [0:n.edges.size] do
    let (a, b) := n.edges[i]!
    if a == v || b == v then r := r.push i
  return r

/-- Number of ordered pairs of edges without a common vertex (`|𝒫|`). -/
def Net.nPairs (n : Net) : Nat := Id.run do
  let mut c := 0
  for i in [0:n.ne] do
    for j in [0:n.ne] do
      if !n.adj i j then c := c + 1
  return c

/-- The checks of SPEC §2.3 (all but `EdgesIntersect`), named like the Rust error variants. -/
def Net.validate (n : Net) : Except String Unit := do
  for i in [0:n.edges.size] do
    let (a, b) := n.edges[i]!
    if a ≥ n.nv || b ≥ n.nv then throw s!"VertexOutOfRange: edge {i} = ({a}, {b}), {n.nv} vertices"
    if a == b then throw s!"SelfLoop: edge {i} at vertex {a}"
  for i in [0:n.edges.size] do
    for j in [i+1:n.edges.size] do
      let (a, b) := n.edges[i]!
      let (c, d) := n.edges[j]!
      if (a == c && b == d) || (a == d && b == c) then throw s!"DuplicateEdge: edges {i} and {j}"
  for v in [0:n.nv] do
    if (n.incident v).isEmpty then throw s!"IsolatedVertex: vertex {v}"
  for v in [0:n.nv] do
    if !n.pos[v]!.finite then throw s!"NonFinitePosition: vertex {v}"
  for i in [0:n.edges.size] do
    let (a, b) := n.edges[i]!
    let d := n.pos[a]! - n.pos[b]!
    if d.x == 0.0 && d.y == 0.0 && d.z == 0.0 then throw s!"ZeroLengthEdge: edge {i}"
  for a in [0:n.nv] do
    for b in [a+1:n.nv] do
      let d := n.pos[a]! - n.pos[b]!
      if d.x == 0.0 && d.y == 0.0 && d.z == 0.0 then throw s!"CoincidentVertices: {a} and {b}"
  if n.nPairs == 0 then throw "NoInteractingPairs: no two edges are disjoint"

/-- The coordinates carried by the unknowns: all three, or two in planar mode (SPEC §2.2). -/
def coordsOf (planar : Option Nat) : Array Nat :=
  match planar with
  | none => #[0, 1, 2]
  | some ax => (#[0, 1, 2] : Array Nat).filter (· != ax)

end Ddg.RepulsiveCurves
