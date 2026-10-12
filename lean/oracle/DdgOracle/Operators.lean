import DdgOracle.Formulas
import DdgOracle.Topology
import DdgOracle.Dense

/-!
# The discrete operators of a mesh, at `Float`

This file only gathers data: every number is produced by one of the scalar-generic formulas of
`DdgOracle.Formulas.Mesh`, instantiated at `Float`.
-/

namespace Ddg

instance : Inhabited (Point Float) := ⟨⟨0.0, 0.0, 0.0⟩⟩

/-- Corner data of one face. -/
structure FaceGeom where
  /-- area -/
  area : Float
  /-- interior angles at the corners `a`, `b`, `c` -/
  angA : Float
  /-- see `angA` -/
  angB : Float
  /-- see `angA` -/
  angC : Float
  /-- cotangents of the interior angles at the corners `a`, `b`, `c` -/
  cotA : Float
  /-- see `cotA` -/
  cotB : Float
  /-- see `cotA` -/
  cotC : Float
  deriving Inhabited

/-- All operators of a mesh, on all of its vertices (boundary included). -/
structure Operators where
  /-- connectivity -/
  topo : Topology
  /-- vertex positions -/
  pos : Array (Point Float)
  /-- per-face corner data -/
  faceGeom : Array FaceGeom
  /-- total area -/
  totalArea : Float
  /-- angle defect `Ω_v` at interior vertices, geodesic curvature `k_v` at boundary vertices -/
  curvature : FloatArray
  /-- lumped masses, `Σ_v mass[v] = 1` -/
  mass : FloatArray
  /-- dense cotan Laplacian, row-major `nv × nv`, positive semidefinite -/
  lap : FloatArray

namespace Operators

/-- Number of vertices. -/
def nv (o : Operators) : Nat := o.topo.nv

/-- Entry `L_ij`. -/
def L (o : Operators) (i j : Nat) : Float := o.lap.get! (i * o.topo.nv + j)

end Operators

/-- Edge-length function of a mesh with the given vertex positions. -/
def edgeLength (pos : Array (Point Float)) (i j : Nat) : Float := dist pos[i]! pos[j]!

/-- Compute all operators; fails on non-finite positions and degenerate triangles. -/
def buildOperators (topo : Topology) (pos : Array (Point Float)) : Except String Operators := do
  let nv := topo.nv
  if pos.size != nv then throw s!"expected {nv} vertex positions, got {pos.size}"
  for v in [0:nv] do
    let p := pos[v]!
    if !(p.x.isFinite && p.y.isFinite && p.z.isFinite) then
      throw s!"vertex {v} has a non-finite coordinate"
  let len := edgeLength pos
  -- faces
  let mut faceGeom : Array FaceGeom := Array.mkEmpty topo.faces.size
  for fi in [0:topo.faces.size] do
    let f := topo.faces[fi]!
    let a := len f.b f.c
    let b := len f.c f.a
    let c := len f.a f.b
    if !(a.isFinite && b.isFinite && c.isFinite) then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) has an edge length that overflows"
    if !(a > 0 && b > 0 && c > 0) then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) has a zero-length edge"
    if !(b - a + c > 0 && a - b + c > 0 && a - c + b > 0) then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) is degenerate: its edge lengths violate the strict triangle inequality"
    let g : FaceGeom := {
      area := f.area len
      angA := f.angleA len, angB := f.angleB len, angC := f.angleC len
      cotA := f.cotA len, cotB := f.cotB len, cotC := f.cotC len }
    if !(g.area > 0 && g.area.isFinite && g.cotA.isFinite && g.cotB.isFinite && g.cotC.isFinite) then
      throw s!"face {fi} = ({f.a}, {f.b}, {f.c}) is numerically degenerate: area or cotangent not finite and positive"
    faceGeom := faceGeom.push g
  let allFaces := topo.faces.toList
  let total : Float := totalArea allFaces len
  if !(total > 0 && total.isFinite) then throw "total area is not finite and positive"
  -- vertices
  let mut curvature : FloatArray := FloatArray.emptyWithCapacity nv
  let mut mass : FloatArray := FloatArray.emptyWithCapacity nv
  for v in [0:nv] do
    let star := topo.facesAt v
    let k : Float :=
      if topo.isBoundary[v]! then geodesicCurvature star len v else angleDefect star len v
    curvature := curvature.push k
    mass := mass.push (massOf (areaAt star len v) total)
  -- Laplacian: off-diagonal entries on edges, then the diagonal from the rows
  let mut lap : FloatArray := FloatArray.emptyWithCapacity (nv * nv)
  for _ in [0:nv * nv] do
    lap := lap.push fzero
  for i in [0:nv] do
    let star := topo.facesAt i
    for j in topo.neighbors[i]! do
      lap := lap.set! (i * nv + j) (laplacianOffDiag star len i j)
  let off := lap
  for i in [0:nv] do
    lap := lap.set! (i * nv + i) (laplacianDiag nv (fun r c => off.get! (r * nv + c)) i)
  return { topo, pos, faceGeom, totalArea := total, curvature, mass, lap }

end Ddg
