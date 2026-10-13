# Intrinsic triangulations — contract

The contract of the crate `ddg-intrinsic` and of what is built on it, for

> Nicholas Sharp, Yousuf Soliman, Keenan Crane. *Navigating Intrinsic Triangulations.*
> ACM Transactions on Graphics 38(4), 2019. <https://doi.org/10.1145/3306346.3322979>
> (`../NavigatingIntrinsicTriangulations/paper.pdf`, `supplemental.pdf`)

The authors' demo (`../NavigatingIntrinsicTriangulations/navigating-intrinsic-triangulations-demo`,
MIT, built on geometry-central, MIT) may be read and used as a reference and as a target for
differential tests. Code taken from it keeps its licence notice (`LICENSES/`), as in `ddg-bff`.
The paper's Appendix B has pseudocode for every operation; the paper is the specification, this
file says what it does not: our conventions, what is in scope, where we depart, and what must hold.
Where the paper is silent or the reference does something the paper does not say, decide, and
record it in §8 with an id `N<n>`. Never drop an entry.

Workspace conventions of `CLAUDE.md` apply: `f64`, `u32` stored indices, explicit option structs,
error enums that name the element, a report from every operation that can do a variable amount of
work, nothing repaired or normalized silently, no panics on bad input.

## 1. Terms

The *input mesh* is a `ddg_core::TriMesh` with a `ddg_core::Geometry` (edge lengths). It is what
the paper calls the extrinsic mesh; here it is intrinsic data too (no positions are needed
anywhere in this crate). An *intrinsic triangulation* is another triangulation of the same
polyhedral surface: its vertices are points of the input surface (all input vertices, plus
inserted points), its edges are geodesics, its triangles are flat. Angles, areas, curvature and
distances of the surface do not depend on the triangulation; that is what every test leans on.

## 2. Scope

In: the signpost data structure (paper Sec. 3.1), signpost update, tracing queries, vertex update
(3.2), edge flip, vertex insertion, edge split (3.3), point queries and the overlay / common
subdivision (3.4), Delaunay flipping (4.1), Delaunay refinement (4.2) including surfaces with
boundary (N1), and what the cut flow needs to run on an intrinsic triangulation (§6).

Later: vertex repositioning and the optimal Delaunay triangulation (3.3.3, 4.3). Leave room for
`move_vertex`; do not implement it now.

Out: the applications of the paper's Sec. 5.

## 3. Data model

* A halfedge mesh that allows what flips produce (paper Sec. 2.1, Fig. 6): two edges between the
  same vertices, an edge from a vertex to itself, a triangle that uses a vertex twice, vertices of
  degree one. It stays an oriented manifold, with or without boundary. `ddg_core::TriMesh` rejects
  these, so this is a type of its own; conversion from a `TriMesh` is exact, conversion to one is
  §6.
* Per edge a length; per halfedge the signpost angle at its tail (paper Sec. 2.3, 3.1).
* Angle normalization: at an interior vertex angles are scaled by `2π/Θ_i`; at a boundary vertex
  by `π/Θ_i`, measured from the boundary halfedge that has the surface on its left (N2). `Θ_i` of
  an input vertex never changes; an inserted interior vertex has `Θ = 2π`, one inserted on the
  boundary `π`.
* Every vertex knows where it is on the input mesh: an input vertex, a point of an input edge
  (edge, parameter from the lower-numbered end), or a point of an input face (face, barycentric
  coordinates).
* The triangulation never changes the input mesh and holds a reference to it or a copy of what it
  needs; say which and why.

## 4. Operations and what each must guarantee

Each follows the paper's algorithm of the same name; the reference shows the details.

* **Init** from the input mesh (Alg. 4): the same triangulation, signposts set.
* **Flip** (Alg. 5): refused, with the reason, for a boundary edge, for an edge whose two
  triangles do not form a strictly convex quadrilateral when unfolded, and for an edge with an
  endpoint of degree one. Otherwise lengths and signposts are updated. Flipping twice restores
  lengths and angles to rounding.
* **Delaunay flipping** (4.1): the queue algorithm, until every interior edge has opposite angles
  summing to at most `π` within a stated tolerance; a report with the number of flips; a limit on
  the number of flips as an option, reported when hit (the paper gives no bound).
* **Tracing** (Alg. 2, 3.2.2): from a vertex in a signpost direction, and from a point inside a
  face, over the intrinsic triangulation and over the input mesh; returns the end point and, on
  request, every edge crossed with its parameter. Paths through vertices need not be supported
  except where an operation can produce one; then say what happens.
* **Insert a vertex** in a face, **split an edge** including a boundary edge (Alg. 7, 8); the new
  vertex's place on the input mesh comes from a tracing query (3.2.3).
* **Trace an intrinsic edge** over the input mesh: the ordered crossings with input edges. Every
  such trace is checked: it must end at the edge's far vertex within a stated tolerance (the
  paper's Appendix A describes the checks and repairs it uses; here a repair is allowed only if
  it is one of those and is counted in the report). A trace that fails the check is an error that
  names the intrinsic edge and the mismatch. Nothing is perturbed (N3).
* **Point query** (Alg. 11) both ways: a point of the input mesh on the intrinsic triangulation
  and back.
* **Overlay** (3.4.2): the common subdivision as a mesh, every vertex of it located on both
  triangulations, and the pieces of every intrinsic triangle. Also, because the cut flow needs
  exactly this: a straight segment inside an intrinsic triangle, given by two points in
  barycentric coordinates, as the sequence of its pieces on input faces.
* **Delaunay refinement** (4.2): Chew's second algorithm as the paper describes it: flip to
  Delaunay, find a triangle below the angle bound, insert its circumcenter (traced from the vertex
  of the obtuse angle when it lies outside the triangle), repeat. Options: the angle bound (the
  paper reaches 30°; say which default you choose and why), an optional bound on the circumradius
  or edge length, a limit on inserted vertices. Vertices whose cone angle makes the bound
  unreachable are left alone and listed in the report (the paper's rule, and whatever the
  reference does beyond it, N4). With boundary (N1): a circumcenter that would leave the surface
  or encroach on a boundary edge is not inserted; that boundary edge is split instead, in the
  manner of the planar algorithm. The report gives the inserted vertices, the smallest and largest
  angle reached, and why it stopped.

## 5. What must hold (tests)

After any sequence of operations, on closed surfaces and surfaces with boundary, on well-shaped
and on badly shaped inputs (use the generators of `ddg_core`, sheared and needle meshes as in
`crates/ddg-variational-cuts/examples/bad_mesh_study.rs`, the models in `data/models/`, the
deliberately bad meshes of `data/optimal-cones/…/data/`):

1. Every triangle satisfies the triangle inequalities strictly; the total area is the input's.
2. The angle sum at every input vertex is the input's; `2π` (or `π`) at inserted vertices; the
   Gauss–Bonnet sum is unchanged.
3. Signposts are consistent: tracing from a vertex along a halfedge's signpost for the edge's
   length arrives at the other end, on the intrinsic triangulation trivially and on the input
   mesh to tolerance.
4. After Delaunay flipping every interior edge satisfies the Delaunay condition and has a
   non-negative cotan weight; the result does not depend on the order of the queue except at
   cocircular ties, which are to be characterized.
5. Geodesic distances between input vertices computed on the input mesh and on the intrinsic one
   agree to the accuracy of the distance solver, and better after refinement.
6. The overlay has the input's area, is a valid mesh, and each of its faces lies in one input face
   and one intrinsic face.
7. After refinement the angle bound holds everywhere except at the reported vertices.
8. Against the reference program, where it can be built: the same Delaunay triangulation (edge
   lengths per edge, compared as sets) on the same inputs; for refinement, comparable counts and
   angle statistics.
9. Bad input and failed traces are errors with a message that names the element; no NaN anywhere.

Performance: be aware of it (the paper's Fig. 10 and its timing tables are the reference points);
benches for flipping, edge tracing and overlay extraction on the same meshes as the other crates.

## 6. For the cut flow

The flow (`ddg-variational-cuts`) takes a `TriMesh` with a `Geometry` and uses nothing else, so it
runs on an intrinsic triangulation unchanged once that is handed to it as such a pair:

* **Simplicial export.** An intrinsic triangulation that is not a simplicial complex (a self-edge,
  a double edge, a triangle using a vertex twice) is made one by splitting offending edges (an
  intrinsic operation: the surface does not change); then a `TriMesh`, a `Geometry` from the
  intrinsic edge lengths, and for every vertex its place on the input mesh. The number of splits
  is reported.
* **The way back.** The flow's result is a cut: points on edges of the exported mesh and straight
  segments inside its faces (and triple points inside faces). With the overlay operation of §4
  each becomes points on input edges and inside input faces, joined by straight pieces that each
  lie in one input face. That is what the insertion of `ddg_uv::cutting` must be given.
* **State.** The flow's level-set functions live on the vertices of the exported mesh, which has
  more vertices than the input when refinement or the export inserted some. Refinement and export
  are deterministic, so the same input and options give the same vertices again; say what
  identifies a state as belonging to a triangulation, so that a stale one is refused.

The pipeline (`ddg-uv`), the bindings and the Blender operators use this internally, with an
option to turn it off; there is no visible change in Blender other than the result.

## 7. Verification layers

As for the other papers (`CLAUDE.md`): `[formula: <id>]` tags with ids starting `it-`, Lean
theorems for the closed-form formulas (the length of a flipped edge, the signpost update, the
Delaunay condition and the cotan weight, the circumcenter in barycentric coordinates, the
ray–edge intersection of a tracing step, the angle normalization), citations for what Mathlib
cannot reach (termination of Delaunay flipping, existence and uniqueness of the intrinsic Delaunay
triangulation, Rippa's theorem, the guarantees of Chew's algorithm), and the Lean oracle for
differential tests of single operations.

## 8. Decisions and departures

* **N1** The paper's refinement does not consider surfaces with boundary (Sec. 4.2, Sec. 7). Here
  boundary edges are split when encroached, as in the planar algorithm. The exact rule and why it
  terminates, or that termination is not proved, goes here.
* **N2** The paper describes the angle normalization at interior vertices only. The boundary rule
  of §3 is the reference's, to be confirmed against it.
* **N3** The paper handles the meshes on which tracing fails in double precision by perturbing
  vertex positions or by quad precision (Sec. 6.2). Here such a failure is an error; nothing is
  perturbed.
* **N4** Vertices with small cone angles in refinement: the paper's one-sentence rule against what
  the reference does.
