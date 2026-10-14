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

Entries N10 to N16 were added with the second half of the crate (insertion, point queries,
overlay, refinement, export); N1 and N4 are written out there. N17 was added with the replay of
the Lean oracle's fixtures. Numbers are from the tests named
and from `crates/ddg-intrinsic/README.md`.

* **N10** *Insertion and edge split* (Alg. 7, 8). The lengths of the new edges are the
  distances of the point from the corners in the triangle laid out in the plane
  (`it-point-distances`; the paper's closed form cancels in thin triangles). The place of the
  new vertex on the input mesh is found before anything is changed, by a trace over the input
  mesh from the nearest corner of the face; a split of an edge that lies along an input edge
  (always the case on the boundary) is placed by interpolating along that edge instead, which
  is exact. The signposts of an inserted vertex are angles in the frame of the input element
  it lies in (N7), not scaled, since the vertex is flat; a vertex on the boundary has, like
  every boundary vertex, the angle 0 on the boundary halfedge that leaves it and `π` on the
  one that arrives.
  *Snapping (Appendix A).* A traced place whose distance from an edge of its input face is at
  most `1e-12` of the longest edge of that face (`insert::LOCATION_SNAP`) is stored as a point
  of that input edge, and `Insertion::snapped` says so. The measure was a barycentric
  coordinate at first; found by test: in a needle of the input (0.5 long, 1e-4 wide,
  camelshotnoise) a vertex one rounding error from a long edge has the coordinate `2e-12`,
  stayed inside the face, and made an overlay piece with an area of `1e-32`.
  *Refused* (`InsertionRefused`, nothing changed): a point that would make a triangle below
  the shape limit of N5 (`TooThin`, with the nearest edge: split that edge instead); a split of
  an edge with the same face on both sides (`SameFaceOnBothSides`: the edge from a vertex of
  degree one; flip or split a neighbour first); a place that cannot be found (`NotLocated`).
  Inserted vertices are never removed, so vertex, edge and face ids only grow (this answers the
  question left open in N9; the reference's refinement removes vertices again).
* **N11** *Point queries* (Alg. 11). A point of the triangulation is found on the input mesh by
  a trace over the input mesh from a vertex of its face; a point of the input mesh is found on
  the triangulation by a trace over the triangulation from an input vertex of its face. Both
  are exact inverses to rounding (`point_queries_are_inverse_to_each_other`: `1e-9` of the
  mean edge after flips, insertions and refinement). A path over the triangulation that runs
  into an inserted vertex continues straight through it (the vertex is flat); into an input
  vertex it cannot run, because the query starts at one and ends inside the star of one face.
  A path to a point of a boundary edge may reach the boundary a rounding error early; that is
  accepted when the length left is at most `1e-9` of the path, otherwise the query is an error.
* **N12** *Overlay* (3.4.2). Built per face of the triangulation from the edge traces alone
  (no positions): the numbers of crossings on the three sides decide which segments of input
  edges cross the face (`overlay` module documentation: corner-cutting segments, or a fan from
  one corner), and the input edges of the crossings that a segment joins must agree. Traces
  that fit neither case, or whose labels disagree, are refused (`InconsistentOverlay`, naming
  the face) instead of being paired wrongly; traces of another triangulation are refused. The
  paper pairs crossings the same way (its cases 1 and 2) without the check.
  Result: the vertices of the triangulation, then the crossings, each as a point of both
  triangulations; faces as convex polygons with 3 to 6 corners, each with its input face and
  its face of the triangulation. Invariants tested (`tests/overlay.rs`, §5 item 6): the pieces
  in an input face have its area and those in a face of the triangulation have its area, each
  within `1e-9` of the total area, and all together the total within `1e-11`; every corner of
  a piece lies in the closed input face and in the closed face of the triangulation that the
  piece names; every piece is convex and counter-clockwise (every turn strictly positive);
  fan-triangulated, the overlay is a manifold `TriMesh` with the Euler characteristic, the
  boundary loops and the components of the input, and no triangle of it is degenerate. Not
  excluded in general: a piece of the size of a rounding error where a vertex lies that close
  to an input edge without being on it (N10 removes the case that occurred).
* **N13** *Segments, and the way back for a cut* (§4, §6). `segments_on_input` takes nodes
  (points of the triangulation) and segments (a face and two nodes). Every node is located on
  the input mesh once, so the paths of two segments that share a node meet in the same point,
  bit for bit. Every segment is then traced over the input mesh from its first node towards
  its second and must arrive there within `tolerance` times its length (the tolerance and the
  check of N8, default `1e-6`); a segment that does not is an error naming the face
  (`SegmentTrace`), nothing is connected or dropped. The result per segment: its points (the
  nodes at the ends, crossings of input edges between) and the input face of every piece.
  `segment_on_input` is the same for one segment given by barycentric coordinates, which also
  say which corner is meant in a face that uses a vertex twice. Measured on cuts of exported
  meshes (`tests/export.rs`): nodes located within `3e-15`, lengths of the pieces equal to the
  segment within `7e-13`, mismatch at the far node below `7e-13`.
* **N1, written out.** *Boundary.* Ruppert's rule with diametral circles. (a) A boundary edge
  is *encroached* when the angle opposite it exceeds `90° + 1e-9` (in a Delaunay triangulation
  that is the case exactly when some vertex lies inside the circle on the edge as diameter).
  Encroached boundary edges are split before any triangle is looked at. (b) The circumcentre
  of a triangle is not inserted when the path to it leaves the surface (the boundary edge it
  leaves through is split), or when it lies inside the circle on a boundary edge of one of the
  triangles its insertion would remove (those whose circumcircle contains it, found by
  unfolding around the point; that edge is split). The triangle is then looked at again. (c) A
  boundary edge is split at its midpoint, except an edge at a boundary vertex whose interior
  angle is below 60°: there the split is at the distance `u · 2^k` from that vertex that is
  nearest half the edge (`u` the shorter of the two boundary edges at the vertex when
  refinement starts), kept within the middle half of the edge, so that the two boundary edges
  at the vertex carry vertices at the same distances and stop encroaching on each other
  (Ruppert's concentric shells). (d) A triangle whose circumcentre is still not insertable
  after 64 such splits is given up and reported (`NotInsertable`); this did not occur.
  The reference instead inserts the midpoint of the boundary edge that the path hits, tests no
  encroachment, and removes inserted vertices near a split again; on the sheared grid it stops
  with angles of 9.46° (25°) and 2.51° (30°) and largest angles of 153° and 164°, where this
  rule ends with every angle at the bound except the two corners of 9.46° that the input
  boundary has, and a largest angle of 102° and 115°.
  *Termination.* **Closed surface, every vertex with an angle sum of 60° or more, bound up to
  30°: argued, not proved.** The planar argument carries over: the circumdisk of a Delaunay
  triangle is an immersed disk without vertices; a triangle with an angle below 30° has a
  circumradius above its shortest edge; an inserted circumcentre is therefore farther from
  every vertex than the shortest edge of its triangle, and a geodesic loop from it around one
  vertex of angle sum `Θ ≥ 60°` at distance `d ≥ R` has the length `2 d sin(Θ/2) ≥ R`; so no
  edge shorter than the shortest present appears, unless a loop around several cone points or
  a handle is shorter, and the length of those is bounded below by a constant of the surface;
  with lengths bounded below, the area bounds the number of vertices. Not proved: the bound on
  those loops, and that the triangulation after the flips has no shorter edge than the
  distances say. **With boundary: observed.** Ruppert's rule is proved to end in the plane for
  bounds up to 20.7° (boundary angles of 60° or more; with shells below that); at 25° and 30°
  its end is observed, as it is in the plane: status `Done` at 25° and 30° on every mesh tested
  (17 generated shapes, 16 models and bad meshes; tetrahedra with a needle apex also at 20°,
  28° and 33°). Chew's second
  algorithm, which is proved to 26.5°, removes vertices near a split boundary edge; that is not
  done here (N10). **The limit on insertions is the net**: reaching it is reported as status
  `InsertionLimit` with every corner still below the bound listed.
* **N4, written out.** *Small angle sums.* (a) A corner at a vertex whose angle sum (interior
  angle, on the boundary) is below the bound cannot meet the bound in any triangulation. Such
  corners do not make a triangle one to refine; the other two corners of the triangle do. They
  are listed in the report with the reason `AngleSum` and the angle sum, and
  `min_angle_reachable` is the smallest angle over all other corners. (b) An interior vertex
  with an angle sum `Θ` between the bound and 60° can meet the bound only with few triangles
  around it, in the limit one triangle with both other corners at the same vertex (a
  self-edge); refinement ends there with such triangles (tetrahedra with apex sums from 26° to
  59°: `tests/refinement.rs`). (c) *Needle guard.* A circumcentre at distance `d` from an
  interior vertex with `Θ < 60°` would get an edge to itself around that vertex of length
  `2 d sin(Θ/2)` (tested against a traced loop). If that is shorter than the shortest edge of
  the triangle to be refined, the point is not inserted and the triangle is set aside with the
  reason `NeedleVertex`; it is looked at again when it changes. This is what keeps (N1) the
  shortest edge from shrinking around a needle. In the final reports it occurred once: at a
  vertex whose angle sum equals the bound to rounding. The paper's rule ("vertices of degree
  one are not refined further", 4.2) is the special case; the reference skips, without
  reporting, every triangle with a vertex of angle sum below 60°, with a vertex of degree one,
  or between two boundary edges. With (a) to (c) the bound is reached at every corner where a
  triangulation can reach it on every mesh tested, the reference's included
  (camelshotnoise: 22 corners of 0.042° at spikes with an angle sum of 0.042°, everything
  else at the bound). **Termination with such vertices: observed**, not argued: the guard sees
  only vertices of the triangles the insertion would remove.
* **N14** *Refinement, the rest.* Departure from §4: the circumcentre is found by a trace over
  the triangulation from the barycentre of the triangle, as the reference does, not from the
  vertex of the obtuse angle: a trace cannot start at a corner into the face without choosing
  a wedge, and from the barycentre the path leaves the triangle through the side it should.
  A circumcentre within the shape limit of an edge (N5, N10) is inserted by splitting that
  edge at the foot of the point. Order: the triangle of largest area first (a heap; a face is
  looked at again whenever it changes), so the result is deterministic: the same input and
  options give the same triangulation bit for bit. Default bound **25°**, not the paper's 30°:
  both are reached on every mesh tested, 30° costs two to three times the inserted vertices of
  25° (README table) for a largest angle of 120° instead of 130°, and the user's criterion is
  robustness and speed. Allowed: 0 to 35°; 0 splits encroached boundary edges only. Option
  `max_circumradius`: triangles above it are refined too. Limit: `InsertionLimit`, default 10
  vertices per vertex and at least 1000 (a small mesh with thin triangles needs many times
  its own vertices: needles and caps 6: 405 for 37).
  Against the reference (§5 item 8, `tests/refinement.rs`): the inserted counts agree within
  1 % at both bounds on 12 of the 15 inputs (identical on spot, bob, bunny, bump, the
  icospheres and the long torus; blub 3651 against 3646 at 30°); camelshotnoise has 2.6 %
  fewer at 30°, the thin cylinder 66 against 68 and 150 against 146, and the sheared grid 256
  against 170 and 395 against 326, where the boundary rules differ (N1). The smallest angle is
  never below the reference's; the largest is within 0.2° of it except on the sheared grid.
* **N15** *Simplicial export and the identity of a state* (§6). The triangulation is made a
  simplicial complex by splitting, at their midpoints, every edge from a vertex to itself and
  all but one of every group of edges between the same two vertices, repeated until none is
  left; with `ExportOptions::delaunay` (default) it is flipped to Delaunay before and after,
  and splits and flips alternate until both hold. The export is the triangulation in that
  state as a `TriMesh` (same vertex, face and corner ids) with a `Geometry` from the edge
  lengths, the place of every vertex on the input mesh, the map from mesh edges to edges of
  the triangulation, and a report: splits, flips, whether it is Delaunay (tested on the
  result), its smallest and largest angle.
  **Splitting does not keep an angle bound** (found by test,
  `the_export_keeps_the_angles_only_where_nothing_is_split`): 34 of 50 refined triangulations
  need no split and are handed out unchanged; the others lose the bound where they are split
  (long torus at 25°: 9.6° after, 65 corners below; camel_noise2: 7.7°, 100 corners;
  camelshotnoise: 0.014° and a largest angle of 179.96°). Two reasons. Three vertices on a
  geodesic loop around a thin tube or a spike are almost collinear, and their Delaunay
  triangulation has slivers. And at an interior vertex with an angle sum below three times
  the bound no simplicial complex has the bound, because it has at least three corners there,
  where refinement had one (N4 b). Option `ExportOptions::refine`: refine between the rounds
  of splits. That restores the bound on thin tubes and handles in one to three rounds (long
  torus: 25.7°, 35 more vertices). At spikes it does not converge, so the refinements stop
  when two in a row leave no fewer edges to split, or when one leaves more than twice as many
  (that one is undone; without this camelshotnoise grows from 17 000 to 580 000 vertices in six
  rounds), and the export ends with splits and flips. `SimplicialExport::refinement` carries
  the report of the last refinement exactly when the export ended with it; otherwise
  `min_angle` and `max_angle` say what the splits left. Open: a refinement that knows the
  result must be simplicial (corners at vertices below three times the bound exempt, new
  vertices around a spike placed at equal distance from it instead of on the loop).
  *State.* `fingerprint()` is an FNV-1a hash (64 bits) of the number of vertices and faces,
  the vertices of every face in order with the bit patterns of the lengths of its sides, and
  the place of every vertex on the input mesh (kind, element, bit patterns of the
  coordinates). Refinement and export are deterministic (N14), so the same input mesh and the
  same options give the same fingerprint, and any flip, insertion, split or other geometry
  gives another (collisions of the hash are not excluded). A state of the flow is kept with
  the fingerprint of the export it was made for; `SimplicialExport::check_state` refuses
  another fingerprint (`StaleState`, both values in the message) and a wrong number of values.
* **N16** *§5 item 5, "better after refinement": measured, not found as stated.* On sheared
  lattices with known distances (`distances_after_refinement`, errors of the first-order solver
  at the input vertices in lattice widths): on the input mesh the mean error is 1 to 7 widths,
  on the intrinsic Delaunay triangulation 0.05 to 0.26, after refinement to 25° or 30° 0.07
  to 0.26, with the circumradius bounded by half a width 0.05 to 0.29. Refinement is as good as
  flipping, far better than the input, and not better than flipping: the solver's error is of
  first order in the mesh width, which an angle bound does not reduce, and refinement leaves
  angles up to 130°, which the Delaunay triangulation of a lattice does not have. What
  refinement improves for the cut flow is measured in `tests/export.rs` (README, end to end).

* **N17** *The Lean oracle's fixtures, and what they made exact.* `lean/oracle` has an
  implementation of single operations of this paper written from the paper and this file
  without this crate's code; `tests/oracle.rs` replays its 18 fixtures, written from the
  oracle's README without the oracle's code, and compares everything they give (angle sums;
  triangulations before and after a flip and after Delaunay flipping as sets of faces with
  lengths and signposts; the unfolded diamond; refusals; report and measures; traces from
  vertices with every crossing; the trace of every halfedge of three flipped triangulations;
  circumcentres and their traces): 23 000 comparisons, largest difference `4e-14`, tolerance
  `1e-9`. Sentences this file did not have, now exact:
  - *First halfedge (N9).* The halfedge of a vertex with the signpost 0 is, at a boundary
    vertex, its outgoing halfedge without a twin, and at an interior vertex its outgoing
    halfedge with the smallest id (`InputSurface::new` chooses it so, whatever
    `TriMesh::outgoing_halfedge` returns). The signpost of every other halfedge is the scale
    of N2 times the sum of the corner angles between the first halfedge and it,
    counter-clockwise. Any other choice shifts all signposts of a vertex by a constant; the
    oracle chose the same (its reading 50), so signposts are compared as they are.
  - *The weight and its tag (N6).* The cotan weight of an edge is `w = ½ (cot α + cot β)`:
    that is what `edge_cotan_weight`, `DelaunayReport::min_cotan_weight` and the tolerance of
    `DelaunayCriterion::CotanWeight` mean. The function tagged `it-cotan-weight`
    (`formulas::opposite_cotan_sum`) returns `2w = cot α + cot β`. The oracle's `MINWEIGHT` is
    `2w`; that is the one number in which the two sides differ, and the oracle's to change or
    to name (the replay accepts `w` and `2w`).
  - *The number of flips (N6).* "All edges in ascending order" means the edge ids of the
    triangulation (for a fresh one those of `ddg_core::TriMesh`, sorted vertex pairs). The
    number of flips depends on that order and is not part of the contract: on the oracle's
    sheared grid 119 here, 121 in the oracle (edges ordered by halfedge), 119 to 124 over 600
    random orders with no test within `1e-6` of a tie. The triangulation that results does
    not depend on it (identical in all six fixtures, among them the tube torus with 36
    self-edges).
  - *A direction equal to a signpost (N7).* It runs along that edge: the path ends on the
    edge (`Completed`, the end a point of the edge) if the length is shorter than the edge,
    and at the far vertex otherwise (`Completed` if the two lengths are equal,
    `TraceStop::Vertex` with the rest not travelled if length is left). At a boundary vertex
    the direction `π` exactly is the boundary edge that arrives (N2), which has no halfedge
    leaving the vertex; a direction above `π` leaves the surface at once (`Boundary`, length
    0). Found by the test of this sentence: for `π` the angle inside the last wedge came out a
    rounding error short of the corner in some cases and the path ended beside the edge,
    inside the face; `π` is now taken as the edge.
  - *Crossings (N7).* The boundary edge at which a path leaves the surface is its end point
    and is not among its crossings (the oracle lists it as the last crossing). A path that
    ends at a vertex may or may not cross an edge at that vertex within rounding; the check
    of N8 ends the path at the target instead, so `EdgeTraces::crossings` has no such
    crossing (46 of 216 untargeted traces on the oracle's tube torus have one on one side
    only).
  - *Edge split lengths (N10).* `split_edge(e, t)` puts the vertex at `(1 − t)·tail +
    t·head`; the half from the tail keeps the id `e` and has the length `t ℓ`, the half to
    the head `(1 − t) ℓ`. With `b_i` the barycentric weight of end `i` the distance from `i`
    is `ℓ_ij b_j`. The paper's Alg. 8 prints `ℓ_ip = ℓ_ij b_i`; that is the same only if its
    `b_i` is the weight of the other end, and is taken as a slip (the oracle's reading 51 is
    the same). `it-point-distances` with the point on a side gives these two lengths and the
    lengths of the new edges to the opposite vertices.
  - *Areas (N5).* Heights and the unfolded diamond use Kahan's formula for the area here and
    Heron's in the oracle; in the fixtures the two differ by `2e-16` of the edge at most. No
    fixture has a flip near the shape limit.
  - *The circumcentre's trace (N14).* From the barycentre of the face, in the direction and
    over the distance that the layout of the face gives (corner 0 at the origin, corner 1 on
    the x-axis); a circumcentre inside its face is not traced by refinement. Same in the
    oracle (reading 56): 464 traces, end points within `2e-14`.
  Not covered by the oracle: insertion and split as operations, point queries, the overlay,
  refinement, the export.

Formula ids in use: `it-angle-normalization`, `it-signpost-update`, `it-unfold-diamond`,
`it-flip-length`, `it-flip-convexity`, `it-delaunay-condition`, `it-cotan-weight`,
`it-unfold-apex`, `it-ray-edge-intersection`, `it-barycentric-coordinates`,
`it-point-distances`, `it-circumcenter`; from `ddg-core`: `cotan-laplacian`, `lumped-mass`,
`angle-defect`, `boundary-geodesic-curvature`. The lengths of an edge split are
`it-point-distances` with the point on a side. The Lean theorems and citations of §7 are not
written yet, so `scripts/check-formulas.sh` reports the `it-` ids as missing.
