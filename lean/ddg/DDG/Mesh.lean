import DDG.Mesh.GaussBonnet
import DDG.Mesh.Gluing
import DDG.Mesh.GluingFlip
import DDG.Mesh.Orientation
import DDG.Mesh.TriMesh

/-!
# DDG.Mesh: triangle meshes and triangulated surfaces

**Scope.** Combinatorial and metric structure of triangle meshes, in two models: `TriMesh`, a
simplicial mesh with faces as vertex triples and edge lengths, and `Gluing`, a triangulated surface
as gluing data (a Δ-complex, in which two edges may join the same vertices and a face may be glued
to itself, as intrinsic triangulations require).

**Main results.**
* `TriMesh`: per-face lengths, angles and cotangents, validity, scaling.
* `GaussBonnet`: angle defect, boundary curvature, Euler characteristic and discrete
  Gauss–Bonnet (`TriMesh.gauss_bonnet`, `gauss_bonnet_closed`).
* `Orientation`: consistent orientation, halfedges and twins (`TriMesh.Oriented`, `TriMesh.twin`),
  boundary halfedges are boundary edges, telescoping of antisymmetric sums (`sum_halfedge_antisymm`).
* `Gluing`: halfedges, twins, edges and vertices of gluing data, `χ = V − E + F`, labellings.
* `GluingFlip`: the combinatorial edge flip (`Gluing.flip`, `Gluing.Flippable`) keeps the vertices
  (`flipVertexEquiv`), the Euler characteristic and the boundary; flipping twice is an isomorphism.

**Sources.** A. Hatcher, *Algebraic Topology*, 2002, Sec. 2.1 (Δ-complexes); N. Sharp, Y. Soliman,
K. Crane, *Navigating Intrinsic Triangulations*, ACM TOG 2019, Sec. 2.1 and Alg. 5;
A. Bobenko, B. Springborn, *A discrete Laplace–Beltrami operator for simplicial surfaces*, 2007,
Sec. 2; C. Indermitte, Th. Liebling, M. Troyanov, H. Clémençon, *Voronoi diagrams on piecewise flat
surfaces and an application to biological growth*, TCS 2001, Sec. 4.

**Depends on:** `DDG.Geometry`.
-/
