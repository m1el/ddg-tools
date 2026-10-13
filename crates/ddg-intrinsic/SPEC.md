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

Entries N5 to N9 were added with the first half of the crate (data model, flips, Delaunay
flipping, tracing); N2 and N3 were completed there. The numbers are those of
`crates/ddg-intrinsic/README.md`.

* **N2, confirmed.** geometry-central 76da12c, `signpost_intrinsic_triangulation.cpp`
  (`updateAngleFromCWNeighor`, `vertexAngleScaling`): at a boundary vertex the first halfedge is
  the one whose twin is exterior, its angle is 0, the boundary edge that arrives has the angle
  `Θ`, and directions are scaled by `π/Θ` where a vector in the tangent space of the vertex is
  needed. The reference stores the angles unscaled, in `[0, Θ)`, and scales on use. Here the
  stored signpost is the scaled angle (paper Sec. 3.1): in `[0, 2π)` at an interior vertex, in
  `[0, π]` at a boundary vertex. The boundary edge that arrives at a boundary vertex has no
  halfedge leaving it (numbering of `ddg_core`, N9); its direction there is the constant `π`.
* **N3, completed.** What the check of a trace is: N8. Appendix A's measures, one by one:
  intersections are computed only with edges that do not contain the origin of the ray (here by
  construction, N7); edge parameters are clamped to `[0, 1]` (not needed, N7); an end point with
  a barycentric coordinate near zero is snapped to the edge (not done; it belongs to vertex
  insertion); crossings after the one into the star of the target are removed and a path that
  stops short is connected through a neighbouring triangle (neither is done, N8). What is moved
  into its range, and reported as `Trace::clamp` instead of silently: a barycentric coordinate
  of an end point that comes out negative, and an angle inside the corner of a vertex that comes
  out beyond the corner. Both were zero in every run recorded in the README.
* **N5** *Strictly convex, in floating point.* A flip is refused (§4) unless the unfolded
  quadrilateral is strictly convex. Decided here: it is also refused when the quadrilateral is
  convex so barely that a triangle `T` the flip would create has
  `2 A_T < ρ · ℓ_new · ℓ_max(T)`, `ρ = 1e-6` (`formulas::FLIP_SHAPE_LIMIT`; refusal `TooThin`).
  Reason: the new length is rounded to a double, and for a triangle with sides `a, b, c` and
  area `A`, `∂γ/∂c = c/2A` and `|∂α/∂c| ≤ a/2A`: a relative error `ε` of the new side moves the
  angles of `T`, as its three lengths give them, by up to `ε / shape` with
  `shape = 2A / (ℓ_new ℓ_max)`. Below `shape ≈ √ε = 1e-8` the lengths say nothing about `T`.
  Found by test: on a grid, flipping an edge whose quadrilateral has three collinear vertices
  passed the sign test by rounding and left a triangle of height `1e-8 ℓ`; angle sums were off
  by `1e-8` and edge traces missed their vertex by `3e-8 ℓ`. With the limit one flip moves an
  angle by `2e-10` at most. Delaunay flipping is not affected on any mesh tested: the smallest
  shape number of a Delaunay flip was `4e-4` (camelshotnoise; `1e-2` or more on the others).
  The paper's remedy for the same loss is a Taylor expansion of the length formula for nearly
  collinear triples (Appendix A); the reference refuses `flipEdgeIfPossible` when a new signed
  area is below `1e-6` of the two together and tests nothing in `flipEdgeIfNotDelaunay`.
  The length of the flipped edge is computed from the unfolded quadrilateral with heights from
  Kahan's area formula (`it-flip-length`), not by the paper's law of cosines (Fig. 9).
* **N6** *The Delaunay test, its tolerance, ties, the limit.* The paper's Eq. 4 has no
  tolerance; in floating point both diagonals of a cocircular quadrilateral can fail it in turn.
  Default (`DelaunayCriterion::AngleSum`): flip when `α + β > π + 1e-10`. The reference flips
  when the cotan weight `½ (cot α + cot β) < −1e-6`; that test is
  `DelaunayCriterion::REFERENCE` and is what the comparison with the reference uses. The two
  differ only at edges whose weight lies in `(−1e-6, 0)`, which the default flips and the
  reference's test leaves: 6 edges on `bob`, none on the other 14 inputs of the comparison.
  The queue is the reference's: all edges in ascending order, first in first out, the four
  neighbours of a flipped edge appended unless queued. An edge that fails the test but whose
  flip is refused (N5) is left and counted (`n_refused`, status `Blocked`); none occurred.
  Limit: `FlipLimit`, default 64 flips per edge of the triangulation, reported as status
  `FlipLimit` (the reference has none).
  *Ties* (§5 item 4): an interior edge is counted as a tie when `|α + β − π|` (or `|weight|`)
  is within the tolerance. The result is the Delaunay triangulation up to the diagonals of
  ties: two orders of flipping differ in at most `2 · n_ties` edges and in none where
  `n_ties = 0`, and the Dirichlet energy of every vertex function is the same for both.
  *Boundary.* Boundary edges are never flipped, so a boundary edge keeps the weight `½ cot` of
  its one opposite angle, negative if that angle is obtuse: on a surface with boundary the
  Laplacian after flipping has no positive off-diagonal entry between vertices that are not
  both on the boundary, and can have one between two boundary vertices.
* **N7** *Tracing.* The triangles along the path are unfolded one after the other into one
  plane in which the path is the positive x-axis (paper Sec. 3.2.2); the edge through which the
  path leaves a triangle is decided by the sign of the y-coordinate of the triangle's third
  vertex (`it-unfold-apex`), and the crossing is `it-ray-edge-intersection`, whose parameter
  lies in `[0, 1]` without clamping. The reference instead carries a direction from face to
  face in barycentric coordinates and projects it into range (constants `1e-12` and `1e-9`,
  not counted); nothing of that kind is done here.
  Starts: a vertex with a normalized direction (the wedge is the outgoing halfedge with the
  largest signpost not above the direction, Alg. 2), a point of an edge, a point of a face. A
  query that starts at the end of an edge or the corner of a face is refused: start at the
  vertex.
  Paths through vertices: a path exactly along an edge is followed to a point of that edge or
  to its end vertex; a path that reaches a vertex with length left stops there
  (`TraceStop::Vertex`) and is not continued; a path that reaches the boundary stops there
  (`TraceStop::Boundary`). No operation of §4 needs more.
  Frames (paper Sec. 2.3): a vertex: N2; a point of an edge: angle zero along the edge, on the
  input mesh from its lower-numbered vertex; a point of face `f`: angle zero along halfedge
  `3 f`. A trace returns the direction that points back along the path in the frame of its end
  (`back_direction`): the first signpost of a vertex inserted there (Sec. 3.2.3).
  Limit on the crossings of one query: `100 · n_faces + 100` unless given; reaching it is an
  error (`CrossingLimit`).
* **N8** *The check of an edge trace.* The place of the edge's far vertex on the input mesh is
  unfolded with every input triangle along the path that contains it. The trace succeeds when
  that place comes within `tolerance · ℓ` of the end of the path (`EdgeTraceOptions::tolerance`,
  default `1e-6`; `ℓ` the edge's length), and the path ends there. It fails when the path has
  been followed for `ℓ (1 + tolerance)` without that, or left the surface, or ran along an
  input edge into another vertex, or reached the crossing limit; the error names the edge, its
  end points and the smallest distance reached. This is in place of Appendix A's rule (the
  last crossing must enter the star of the target, later crossings are removed, a path that
  stops short is connected through a neighbouring triangle), which has no metric bound; here no
  crossing is removed or added. An edge that has not been flipped since the triangulation was
  made lies along its input edge and is not traced (`edge_on_input`; the reference's
  `edgeIsOriginal`). An edge flipped away and back is traced and recognized by its trace.
* **N9** *Data model* (§3). Halfedges are numbered as in `ddg_core::TriMesh` (`h = 3 f + k`,
  none outside the surface), so a triangulation fresh from the input mesh has its ids; the twin
  of a halfedge may lie in the same face. Edges have ids of their own, kept under flips; the
  direction of an edge is that of `edge_halfedge(e)`, which a flip turns counter-clockwise (as
  the reference). A flip moves two halfedges between the two faces, so halfedge ids are not
  stable; vertex, edge and face ids are. The input mesh is held as an immutable copy behind an
  `Arc` (`InputSurface`: the `TriMesh`, its `Geometry`, its signposts and angle sums): no
  lifetime parameter, and clones of a triangulation share it. `to_tri_mesh` hands out a
  triangulation that is a simplicial complex and refuses one that is not; splitting edges to
  make it one is §6. Open for the second half: how elements are numbered when inserted vertices
  are removed again (the reference's refinement does that), and what identifies a state (§6).

Formula ids in use: `it-angle-normalization`, `it-signpost-update`, `it-unfold-diamond`,
`it-flip-length`, `it-flip-convexity`, `it-delaunay-condition`, `it-cotan-weight`,
`it-unfold-apex`, `it-ray-edge-intersection`, `it-barycentric-coordinates`; from `ddg-core`:
`cotan-laplacian`, `lumped-mass`, `angle-defect`, `boundary-geodesic-curvature`. The Lean
theorems and citations of §7 are not written yet, so `scripts/check-formulas.sh` reports the
`it-` ids as missing.
