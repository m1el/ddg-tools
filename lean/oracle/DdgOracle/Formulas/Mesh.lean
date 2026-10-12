import DdgOracle.Formulas.Sum
import DdgOracle.Formulas.Triangle

/-!
# Discrete curvature, cotan Laplacian and lumped mass

Conventions are those of `CLAUDE.md` (the contract between the Rust crates, this oracle and the
proofs). A mesh is given by a list of oriented triangles `faces` and a symmetric edge-length
function `len : Nat → Nat → α` (`len i j` is the length of the edge between vertices `i` and `j`;
it is only ever evaluated on pairs of vertices of a common face).

Every per-vertex or per-edge quantity below is a sum over `faces` in which the faces that do not
contain the vertex (or edge) contribute an exact `0`. The executable therefore passes only the
faces incident to the vertex, in their original order, which gives the same value.
-/

namespace Ddg

variable {α : Type} [Scalar α]

/-- An oriented triangle with vertex indices `a b c` (counter-clockwise). -/
structure Face where
  /-- first vertex -/
  a : Nat
  /-- second vertex -/
  b : Nat
  /-- third vertex -/
  c : Nat
  deriving Repr, BEq, DecidableEq, Inhabited

namespace Face

/-- `v` is a vertex of `f`. -/
def Has (f : Face) (v : Nat) : Prop := v = f.a ∨ v = f.b ∨ v = f.c

instance (f : Face) (v : Nat) : Decidable (f.Has v) := by unfold Has; infer_instance

/-- `{i, j}` is the (undirected) edge of `f` opposite its corner at `f.c`. -/
def EdgeAB (f : Face) (i j : Nat) : Prop := (i = f.a ∧ j = f.b) ∨ (i = f.b ∧ j = f.a)

/-- `{i, j}` is the edge of `f` opposite its corner at `f.a`. -/
def EdgeBC (f : Face) (i j : Nat) : Prop := (i = f.b ∧ j = f.c) ∨ (i = f.c ∧ j = f.b)

/-- `{i, j}` is the edge of `f` opposite its corner at `f.b`. -/
def EdgeCA (f : Face) (i j : Nat) : Prop := (i = f.c ∧ j = f.a) ∨ (i = f.a ∧ j = f.c)

instance (f : Face) (i j : Nat) : Decidable (f.EdgeAB i j) := by unfold EdgeAB; infer_instance
instance (f : Face) (i j : Nat) : Decidable (f.EdgeBC i j) := by unfold EdgeBC; infer_instance
instance (f : Face) (i j : Nat) : Decidable (f.EdgeCA i j) := by unfold EdgeCA; infer_instance

/-- Interior angle of `f` at its corner `a`. -/
def angleA (len : Nat → Nat → α) (f : Face) : α :=
  cornerAngle (len f.b f.c) (len f.c f.a) (len f.a f.b)

/-- Interior angle of `f` at its corner `b`. -/
def angleB (len : Nat → Nat → α) (f : Face) : α :=
  cornerAngle (len f.c f.a) (len f.a f.b) (len f.b f.c)

/-- Interior angle of `f` at its corner `c`. -/
def angleC (len : Nat → Nat → α) (f : Face) : α :=
  cornerAngle (len f.a f.b) (len f.b f.c) (len f.c f.a)

/-- Cotangent of the interior angle of `f` at its corner `a`. -/
def cotA (len : Nat → Nat → α) (f : Face) : α :=
  cornerCot (len f.b f.c) (len f.c f.a) (len f.a f.b)

/-- Cotangent of the interior angle of `f` at its corner `b`. -/
def cotB (len : Nat → Nat → α) (f : Face) : α :=
  cornerCot (len f.c f.a) (len f.a f.b) (len f.b f.c)

/-- Cotangent of the interior angle of `f` at its corner `c`. -/
def cotC (len : Nat → Nat → α) (f : Face) : α :=
  cornerCot (len f.a f.b) (len f.b f.c) (len f.c f.a)

/-- Area of `f`. -/
def area (len : Nat → Nat → α) (f : Face) : α :=
  triArea (len f.a f.b) (len f.b f.c) (len f.c f.a)

/-- Interior angle of `f` at the vertex `v`; `0` if `v` is not a vertex of `f`. -/
def angleAt (len : Nat → Nat → α) (f : Face) (v : Nat) : α :=
  if v = f.a then f.angleA len
  else if v = f.b then f.angleB len
  else if v = f.c then f.angleC len
  else 0

/-- Cotangent of the corner of `f` opposite the edge `{i, j}`; `0` if `{i, j}` is not an edge
of `f`. -/
def cotOpposite (len : Nat → Nat → α) (f : Face) (i j : Nat) : α :=
  if f.EdgeAB i j then f.cotC len
  else if f.EdgeBC i j then f.cotA len
  else if f.EdgeCA i j then f.cotB len
  else 0

/-- Area of `f` if `v` is one of its vertices, `0` otherwise. -/
def areaAt (len : Nat → Nat → α) (f : Face) (v : Nat) : α :=
  if f.Has v then f.area len else 0

end Face

/-! ## Curvature (paper Eqs. 14, 15) -/

/-- Sum of the corner angles at the vertex `v`. -/
def angleSum (faces : List Face) (len : Nat → Nat → α) (v : Nat) : α :=
  sumList (faces.map fun f => f.angleAt len v)

/-- Angle defect `Ω_v = 2π − Σ (corner angles at v)` of an interior vertex (Eq. 14). -/
def angleDefect (faces : List Face) (len : Nat → Nat → α) (v : Nat) : α :=
  2 * Scalar.pi - angleSum faces len v

/-- Geodesic curvature `k_v = π − Σ (corner angles at v)` of a boundary vertex (Eq. 15). -/
def geodesicCurvature (faces : List Face) (len : Nat → Nat → α) (v : Nat) : α :=
  Scalar.pi - angleSum faces len v

/-! ## Cotan Laplacian (paper Sec. 5.1.3), positive semidefinite -/

/-- Off-diagonal entry `L_ij = −½ Σ cot(angle opposite edge ij)`, the sum being over the one or
two faces that contain the edge `{i, j}`. -/
def laplacianOffDiag (faces : List Face) (len : Nat → Nat → α) (i j : Nat) : α :=
  -(sumList (faces.map fun f => f.cotOpposite len i j)) / 2

/-- Diagonal entry `L_ii = −Σ_{j ≠ i} L_ij` of a matrix on `nv` vertices whose off-diagonal
entries are `offDiag`. -/
def laplacianDiag (nv : Nat) (offDiag : Nat → Nat → α) (i : Nat) : α :=
  -(sumTo nv fun j => if j = i then 0 else offDiag i j)

/-- Entry `L_ij` of the cotan Laplacian of a mesh with `nv` vertices. -/
def laplacian (faces : List Face) (len : Nat → Nat → α) (nv : Nat) (i j : Nat) : α :=
  if i = j then laplacianDiag nv (laplacianOffDiag faces len) i
  else laplacianOffDiag faces len i j

/-! ## Lumped mass (paper Eq. 16), normalized by the total area -/

/-- Sum of the areas of the triangles at `v`. -/
def areaAt (faces : List Face) (len : Nat → Nat → α) (v : Nat) : α :=
  sumList (faces.map fun f => f.areaAt len v)

/-- Total area of the mesh. -/
def totalArea (faces : List Face) (len : Nat → Nat → α) : α :=
  sumList (faces.map fun f => f.area len)

/-- `m = (⅓ · area of the triangles at the vertex) / (total area)`. -/
def massOf (areaAtVertex total : α) : α := areaAtVertex / 3 / total

/-- Lumped mass `m_v = ⅓ Σ (area of triangles at v)`, divided by the total area so that
`Σ_v m_v = 1`. -/
def lumpedMass (faces : List Face) (len : Nat → Nat → α) (v : Nat) : α :=
  massOf (areaAt faces len v) (totalArea faces len)

end Ddg
