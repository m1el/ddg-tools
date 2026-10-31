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
(3.2), edge flip, vertex insertion, edge split, vertex repositioning (3.3), point queries and the
overlay / common subdivision (3.4), Delaunay flipping (4.1), Delaunay refinement (4.2) including
surfaces with boundary (N1), the intrinsic optimal Delaunay triangulation (4.3; N19 to N23), and
what the cut flow needs to run on an intrinsic triangulation (§6).

(Vertex repositioning and the optimal Delaunay triangulation were "later" in the first version
of this file and were added on 2026-10-14.)

Out: the applications of the paper's Sec. 5; they live in the crate `ddg-intrinsic-apps` (its
`SPEC.md`, decisions N30 to N49). For its adaptive refinement (Sec. 5.4) this crate gained
`refine_with(options, extra)` on 2026-10-17: Delaunay refinement with an extra per-face
criterion; `refine` is `refine_with` with a criterion that is always false, bit for bit.

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

* **Move a vertex** (Alg. 10): an inserted vertex by a vector in its own frame, inside its flat
  star: new lengths by the law of cosines, signposts by the vertex update, the place on the
  input mesh by a trace. Input vertices never move. What happens to a larger motion, and to
  vertices on the boundary: N19, N20.
* **Optimal Delaunay triangulation** (4.3): split edges longer than a target length, move every
  movable inserted vertex towards the area-weighted mean of the vectors to the circumcentres of
  its triangles (Chen and Holst 2011, Eq. 4.13), flip to Delaunay, repeat until a limit on the
  iterations or a tolerance on the motion. Options and report: N21 to N23.

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
10. After moves of inserted vertices (N19): items 1 to 3, 6 and the point queries hold; on a flat
    input the moved vertex is at its old place plus the displacement and its edges have the
    lengths of the straight segments; a move and the move back restore lengths and signposts; a
    refused move changes nothing; input vertices never move.
11. The optimal Delaunay triangulation (N21): the result is Delaunay after every iteration,
    items 1 to 3 and 6 hold, input vertices are where they were, the run is deterministic, and
    its update at a vertex is the area-weighted mean of the circumcentres in coordinates and
    `−3/(2|ω|)` times the gradient of the energy of Chen and Holst (finite differences).

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
  *Chosen by the caller since N52:* this test is `FlipShapeTest::N5`, the default of every flip
  and of every operation of this crate; a caller may ask a single flip for the reference's
  relative test instead (FlipOut's Alg. 1 does by default).
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
  A place on the input mesh within `1e-12` of an input edge is returned as a point of that edge
  (N29).
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
  segment within `7e-13`, mismatch at the far node below `7e-13`. A segment along a side that
  lies along an input edge is not traced, and one along an input edge otherwise meets its end
  as a point of that edge (N29).
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
    (`formulas::opposite_cotan_sum`) returns `2w = cot α + cot β`. The oracle's `MINWEIGHT` was
    `2w`, the one number in which the two sides differed; it now prints `w` (the six
    `it-delaunay-*` fixtures were regenerated, only their `measures` line changed). The replay
    still accepts `w` and `2w`.
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

* **N18** *Alg. 1 has no reduction modulo 2π; the crate and the oracle subtract 2π once.* The
  paper keeps every signpost `φ` in `[0, 2π)` (Sec. 2.3, 3.1) and updates it by `φ_ik = φ_ij +
  (2π/Θ_i) θ_i^{jk}` (Eq. 3, Alg. 1, `UpdateSignpost`), but Alg. 1 as printed does not reduce
  the result, so as printed the update can leave `[0, 2π)`. `formulas::signpost_update` (and the
  oracle's `signpostUpdate`) subtract `2π` once if the sum reaches `2π`, at an interior vertex;
  at a boundary vertex directions do not wrap and nothing is subtracted (N2). One subtraction is
  exact, not a rounding convention, for `0 ≤ φ < 2π` and `0 ≤ θ ≤ Θ`: then `(2π/Θ) θ ∈ [0, 2π]`,
  so the sum is in `[0, 4π)` and a single subtraction brings it into `[0, 2π)`. The Lean
  theorems `signpostUpdate_mem` (the result is in `[0, 2π)` and equals the sum or the sum minus
  `2π`) and `signpostUpdate_direction` (it has the same cosine and sine as the sum, so it is the
  direction of Eq. 3) in `lean/proofs/DdgProofs/Intrinsic/Signpost.lean` state this for the
  oracle's function, at an interior vertex. Outside the hypotheses (a negative `φ`, a corner
  larger than `Θ`, a signpost of `2π` or more) one subtraction is not enough, and the crate does
  not guard against them.

Entries N19 to N27 were added with vertex repositioning and the optimal Delaunay triangulation
(2026-10-14). The reference implementation has neither (geometry-central 76da12c has
`removeInsertedVertex`, which flips the vertex to degree three and deletes it, and an extrinsic
`smoothByCircumcenter` in `remeshing.cpp` that is not used by the demo); everything here is
written from the paper and the cited papers.

* **N19** *Moving a vertex* (Alg. 10). `move_vertex(v, displacement)`: `displacement` is a
  vector in the frame of `v`, the frame its signposts refer to (N10: the frame of the input face
  or edge it lies in; at a boundary vertex the boundary frame, `0` along the boundary halfedge
  that leaves it, `π` along the one that arrives). Only inserted vertices move; an input vertex
  is refused (`MoveRefusal::InputVertex`).
  *The star.* `v` is flat (angle sum `2π`, or `π` on the boundary), so its star unfolds with `v`
  at the origin and the far end of the halfedge `h` at `ℓ_h (cos φ_h, sin φ_h)`, `φ_h` the
  signpost; in each triangle the end of the arriving side is at the direction `φ_h + θ` (`θ` the
  corner angle from the lengths). A vertex with an edge to itself is refused (`SelfEdge`): its
  star is not a disk around it.
  *When the local update applies.* Exactly when the new place lies in the kernel of the star
  polygon, i.e. every triangle of the star keeps a positive area. For the triangle on the link
  side `d = k − j` twice the new area is `2A − v × d` (`it-move-area`), affine in `v`, so the
  test is exact, and the largest admissible fraction of `v` is `min 2A/(v × d)` over the
  triangles with `v × d > 0`. A move that turns a triangle over or makes it degenerate is
  refused with that fraction (`NotInStar { face, fraction }`); a move that keeps the areas
  positive but leaves a triangle below the shape limit of N5 (twice the area over the product
  of the longer new edge and the longest side, `1e-6`) is refused (`TooThin`). The paper's test
  (`|v|` below the smallest height of a triangle of the star over its link side) is sufficient
  and stricter; it is not required, and `VertexMove::within_height_bound` says whether it held
  (in `tests/reposition.rs` about half of the accepted random moves are beyond it).
  *Lengths* (`it-move-length`): `ℓ' = √(r² + ℓ² − 2rℓ cos α)` with `r = |v|` and `α` the angle
  between `v` and the edge, evaluated as `√((ℓ − r cos α)² + (r sin α)²)`. Only the edges at
  `v` change.
  *The place on the input mesh* is found before anything changes, by a trace over the input
  mesh from the old place of `v` in the direction of `v` for `|v|`. Departure from Alg. 10, whose
  vertex update traces from the neighbour `j₀` to `v`: the segment from the old to the new place
  lies in the star, which has no other vertex, so in exact arithmetic both give the same point;
  this trace is shorter and needs no other frame. Snapping as for an insertion (N10,
  `LOCATION_SNAP`, reported in `VertexMove::snapped`). A trace that does not complete although
  the move is inside the star is an inconsistency and an error (`MoveNotLocated`), not
  repaired; an end at an input vertex is refused (`OnInputVertex`; it cannot happen in exact
  arithmetic, because the star contains no vertex).
  *Signposts* (vertex update, Sec. 3.2.3). Outgoing: the first halfedge of `v` gets the
  direction of the trace back to the old place plus the angle, in the unfolded star, from the
  direction back to the old place to the first neighbour; the others follow corner by corner
  from the new lengths (`it-signpost-update`). Incoming: the halfedge from each neighbour `j`
  gets the signpost of the side of the link that leaves `j`, updated by the new corner at `j`.
  Nothing else changes.
  *Edges along input edges.* After a move off the boundary the edges at `v` are no longer
  marked as lying along an input edge (`edge_on_input`); they are traced (N8) like any other.
  *Boundary.* A vertex inserted on the boundary moves along the boundary only (the paper does
  not consider boundaries, cf. N1): in its frame the boundary is the x-axis, so the second
  component of the displacement must be zero (`LeavesBoundary` otherwise). Its place is
  interpolated on the input edge it lies on, from the parameter of its neighbour along the
  boundary, as for an edge split (N10); the boundary edges keep their input edge. Moving it off
  the boundary would change the surface.
  *Order of the tests* (the first that applies is the refusal): an index out of range, an input
  vertex, a displacement that is not finite, a boundary vertex with a second component, an edge
  to itself; then a zero displacement is accepted and changes nothing; then the star
  (`NotInStar`), the shapes (`TooThin`), the place (`OnInputVertex`, `MoveNotLocated`). The
  oracle tests the edge to itself before the boundary (its reading 62); no fixture has a move to
  which both apply.
  *Refusals change nothing*. `VertexMove::back_direction`
  is the direction of the old place in the new frame: moving by the same distance in it undoes
  the move (lengths to `1e-12`, signposts to `1e-10` in the tests).
* **N20** *A larger motion.* `move_vertex` never shortens a step and never moves a vertex out
  of its star; a move outside the kernel is refused with the admissible fraction, and the caller
  decides. The paper's general route (remove the vertex, retriangulate the flat polygon, insert
  it again at the target, Sec. 3.3.3) is **not implemented**: it would remove an inserted vertex,
  which nothing in this crate does (N10: ids only grow), and the optimal Delaunay iteration does
  not need it (N21: a shortened step followed by flips reaches the same region of the surface in
  later iterations). Open: a target outside the star, reached in one operation.

* **N21** *The optimal Delaunay triangulation* (Sec. 4.3). `optimal_delaunay(options)`. The
  update at a movable vertex `x` (`odt_target`) is Chen and Holst's Eq. 4.13 with uniform density
  and step `α`: `x_new − x = α Σ_T |T| (c_T − x) / Σ_T |T|` over the triangles `T` of the star,
  `c_T` the circumcentre. As the paper says, the vectors to the circumcentres are averaged, each
  computed in the planar layout of its triangle around `x` and rotated into the frame of `x` by
  the signpost of the triangle's halfedge leaving `x`; for an inserted vertex (flat) that is
  the unfolded star. `|T| (c_T − x)` is evaluated without division as
  `(|T| ℓ_ab / 2, ℓ_ab (ℓ_ac² + ℓ_bc² − ℓ_ab²) / 8)` in the frame of the side `ab` leaving `x`
  (`it-circumcenter-vector`); the mean is `it-odt-update`. For a planar star with a fixed link
  this is `−3/(2|ω|)` times the gradient of the energy of Chen and Holst's Lemma 4.2, Eq. 4.5
  (checked by finite differences in `tests/reposition.rs` and in the oracle's selftest;
  machine-checked as a derivative, `hasFDerivAt_odtEnergy` in
  `lean/proofs/DdgProofs/Intrinsic/Move.lean`, which also proves that this planar energy is
  exactly quadratic in `x`, so that the step changes it by `−3|G|²/(4|ω|) α (2 − α)`: a descent
  for `0 < α < 2`, the minimum at `α = 1`; that the optimal Delaunay triangulation minimizes the
  interpolation error is cited, `it-odt-energy`). Chen and Holst halve `α` to keep the elements
  valid, not for descent; so does N21's line search.
  *One iteration*: (1) split every edge longer than the target length at its midpoint, longest
  first (ties by id; edges made by these splits wait for the next iteration); a refused split
  (N10, too thin) is counted, not an error; (2) flip to Delaunay if anything was split;
  (3) move every movable inserted vertex in ascending order of its id, each with the update
  computed from the current triangulation (Gauss–Seidel; the paper does not say), with the step
  `α` of the options; a move that would leave the star or make a triangle too thin (N19) is
  retried with `α/2, α/4, …` at most `max_halvings` times (Chen and Holst's line search, right
  after their Eq. 4.13: "we use α = 1, 1/2, 1/4, … as a simple line search strategy"), then
  the vertex stays for this iteration (counted, `n_not_moved`); a vertex with an edge to itself
  is skipped (counted); (4) flip to Delaunay. The triangulation is Delaunay (the report of the
  flips says so) after every iteration.
  *Stopping*: converged when an iteration split nothing and moved no vertex by more than
  `motion_tolerance` times the mean length of its edges (before the move); otherwise at
  `max_iterations` (status `IterationLimit`) or when the splits reach `max_insertions` (status
  `InsertionLimit`; that iteration's moves and flips are still made). Deterministic: the same
  input and options give the same triangulation bit for bit.
  *Target length*: `TargetLength::Absolute(h)`, `MeanInputEdge(k)` (`k` times the mean edge
  length of the input mesh) or `NoSplits` (only vertices already inserted move, e.g. after
  refinement). The paper names "a user-defined target length" (after Tournois et al. 2009) and
  no default; ours is `MeanInputEdge(1.0)`.
  *Report*: status, iterations, vertices added, every iteration (splits, refused splits, moves,
  halved moves, vertices not moved or skipped, flips, largest motion absolute and relative),
  the largest motion of the last iteration, the last Delaunay report, and the angle and area
  distribution of the result (`ShapeStatistics`: angle quantiles, fractions below 10°/20°/30°
  and above 90°/120°/150°, area max/min and coefficient of variation).
  *Defaults*: 50 iterations, motion tolerance `1e-3`, `α = 1` (Chen and Holst's initial step) with
  up to 8 halvings, the shape bound `1e-2` of N25, at most 10 insertions per vertex (1000 at least,
  as refinement). With them it converged on 10 of 12 meshes of the study (N27), in at most 50
  iterations (spot 25, camelbad 44).
* **N22** *Moves beyond the star in the optimal Delaunay iteration* (N20): the step is halved,
  never clamped to the border of the star, and every halved or abandoned move is counted in the
  report. The paper's general route (removal and re-insertion) is not used.
* **N23** *Boundary in the optimal Delaunay triangulation* (the paper does not consider it; cf.
  N1). `BoundaryVertices::Fixed` (default): inserted boundary vertices (from splits of long
  boundary edges) stay where the split put them, as Chen and Holst keep boundary vertices
  fixed. `AlongBoundary`: they move along the boundary by the component of their update along
  it (N19). `OdtCenters::BarycentersAtBoundary` (default): a triangle with a vertex on the
  boundary contributes its barycentre instead of its circumcentre (Chen and Holst 2011, Sec. 4.2.1
  of the preprint, "for boundary elements, use the barycenter": a stretched triangle at a boundary edge cannot
  be flipped away); `Circumcenters` is the paper's rule everywhere. On a closed surface the two
  are the same.

* **N24** *Snapping a moved vertex onto an input edge* (Appendix A, as N10). A moved vertex
  whose traced place is within `MOVE_SNAP = 1e-12` of the longest edge of its input face, and
  within `1e-12` of the shortest new edge at the vertex, from an edge of that input face is stored
  as a point of that input edge (`VertexMove::snapped`, counted per iteration in
  `OdtIteration::n_snapped`). The threshold is that of an insertion; the second bound keeps the
  relative mismatch it causes on every edge at the vertex below `1e-12` where the input faces are
  much longer than the edges of the triangulation. It matters for the optimal Delaunay iteration:
  a vertex made by splitting an edge along an input edge, whose star is symmetric about that
  edge, is moved along it again and again and stays within rounding of it (1055 moves were
  snapped on blub in 35 iterations).
* **N25** *The shape of the triangles a move makes, in the optimal Delaunay iteration* (found by
  the study, blub). The line search of N21 is also used when a move would leave a triangle of
  the star with a shape number (N5, N19) below `OptimalDelaunayOptions::min_shape` (default
  `1e-2`) and below the smallest shape number the star had before the move; `move_vertex`
  itself only refuses the limit of N5 (`1e-6`). Reason: the signposts at both ends of the moved
  edges are updated through corner angles, whose rounding error is about `ε / shape`; with
  shapes down to `1e-6` they are good to `1e-10` only. Without this bound (and with the shape
  limit of N5 alone) the converged iODT of blub had two of 53 649 edge traces that failed N8: the
  edge from an input vertex to a vertex that had stayed `1e-10` from an input edge (N24) had a
  stored direction `3e-10` rad and a stored length `2e-10` (relative) away from what the
  vertex's traced place says, and its trace ended in the face on the other side of that input
  edge. Raising the snapping threshold to `1e-10` hid it on every mesh tested, at the price of
  a mismatch of up to `1e-7` of an edge on cells 1:16; the shape bound removes the cause: with
  it no trace failed on any mesh of the study and the largest relative mismatch of a trace
  after iODT was `1.3e-11` on blub, `1.3e-10` on camelshotnoise (spikes, where most moves are
  stopped by the bound: 1640 vertex-iterations not moved), `8e-11` on cells 1:16, otherwise
  `4e-12` or less. The iteration also converges in fewer iterations (blub: 35 instead of 50).

* **N26** *The oracle's moves* (as N17). `lean/oracle` has a task `move` written from the paper and
  N10, N19, N21, N23, N24 without this crate's code: an insertion or a split on the fresh
  triangulation, flips, then moves of the new vertex, printing after every move its place, the
  triangulation and the update of N21 (both centre rules), and every refusal with its reason and
  the admissible fraction. Its new signposts are computed differently (directions in the
  unfolded star turned into the new frame, reading 67) and its place by its own tracing.
  `tests/oracle.rs` replays the six fixtures `it-move-*` (flat grid inside a face, across faces,
  onto and along input edges with snapping; refusals; a noisy icosphere with a vertex of degree
  7; a vertex on an interior input edge; a boundary vertex) from the oracle's README. What the
  oracle reads differently is in its readings 59 to 69 (the vertex is named in `move V DX DY`; the
  order of the tests, see N19).

* **N27** *Sec. 4.3 measured* (README, "The optimal Delaunay triangulation"; 12 meshes). At the
  vertex count of refinement the paper's claim reproduces: the areas are more even (coefficient
  of variation lower on every mesh but one, by up to a factor of 5) at the cost of skinnier
  angles (the smallest angle stays that of the flips, since input vertices cannot move; 1 to 8 %
  of the corners below 30° where refinement has none). The condition number of the cotan
  Laplacian with lumped masses is lower than refinement's in 13 of 18 comparisons. iODT takes 3
  to 40 times the time of refinement. In the cut flow (`flow_on_bad_meshes.rs`, TREAT=odt, odt15)
  it brings no gain in robustness or correctness over flips that refinement does not bring, and
  costs time; the default of `ddg_uv::cutting` stays Delaunay flips, iODT is an explicit option
  (`IntrinsicOptions::optimal_delaunay`).

* **N28** *Remeshing for the cut flow measured* (cut paper Sec. 6.1; README, "Remeshing for the
  cut flow"; `crates/ddg-uv/examples/remesh_study.rs`, 1 260 runs on 10 meshes). Uniform
  refinement by midpoint splits and flips, the paper's curvature-adaptive refinement at the same
  vertex count, and iODT, each with the curvature of the triangulation (at the input vertices
  only) or with the curvature density of the input faces carried over per face ("spread"),
  against the flips, on a common yardstick (a fine subdivision of the input). Refinement does
  not make the flow converge (1.2 to 3 times the steps, no run ends by motion), costs 11 to 70
  times per step, and is even or worse than the flips on meshes of a few thousand vertices and
  more; only with the spread and only on meshes of a few hundred vertices is it better (−2 to
  −3 % in the median); without the spread it solves another problem (cones at the input
  vertices; outliers to +1 160 %). Decision: no refinement option for the cut flow; the default
  of `ddg_uv::cutting` stays Delaunay flips. Nothing in this crate changed.
  *A segment along an input edge* (found by the study, fixed afterwards: N29).
  `segments_on_input` failed on 2 of 900 cuts from refined triangulations (`remesh_study repro
  bob_control ka4 dirichlet 1 1`, and `ka2 dirichlet 0.5 1`; deterministic): two segments of the
  cut have both ends on input edge (300, 298) of bob's control mesh (its mirror plane, which is
  also the start plane), next to a vertex the refinement inserted on that edge, so they run
  along the edge; the trace from the first end leaves into one input face and the second end is
  placed in the other ("the place of its end is in none of the input faces the path went
  through"), at every tolerance up to 0.1.

* **N29** *The way back along an input edge* (N11, N13, N28). *Cause.* A point of the
  triangulation that lies on an input edge (on bob: where edges of the triangulation cross the
  mirror plane, on which the symmetric start function of the flow vanishes) is found on the
  input mesh by a trace (N11) that ends a rounding error from that edge, and was returned as a
  point of whichever of the two input faces rounding chose (barycentric coordinate `0` or
  `2e-16` opposite the edge). A segment between two such points runs along the edge; its trace
  stays in one face, and the check of N8 meets a target given as a point of a face only in that
  face (`trace.rs`, `meet_target`), so an end placed in the other face is never met, at any
  tolerance. *Fix.* (a) A point query (`point_on_input`, and the places of the nodes of
  `segments_on_input` and `segment_on_input`) whose end lies within `LOCATION_SNAP = 1e-12` of
  the longest edge of its input face, and within `1e-12` of the edge of the triangulation the
  point lies on (of the shortest side of its face, for a point inside a face), from an edge of
  that input face returns the point of that input edge: the snapping of N10 and N24
  (`snap_onto_input_edge`, shared with insertions and moves), applied to point queries. A point
  of an input edge belongs to both faces of the edge: the check of N8 considers a target on an
  edge in every face that has the edge, and a trace that starts on an edge enters the face its
  direction points into, so the trace of a segment along the edge meets its end in the first
  face and reports that it runs along the edge. A point exactly at a corner (two coordinates
  exactly zero) stays a point of the face. (b) A segment whose two ends lie on one side of its
  face, where that side lies along an input edge (`edge_on_input`), is not traced: its path is
  the piece of the input edge between the places of its ends (both are on that edge: a vertex
  of the side by its stored place, a point of the side by interpolation, N11), and its mismatch
  is the difference between the length of that piece and the length of the segment, checked
  against the tolerance as for a trace (as an edge along an input edge is not traced, N8).
  Found by the new test on a boundary edge: the direction of such a segment from a boundary
  vertex or a point of a boundary edge comes out up to `4e-16` beyond `π`, which points out of
  the surface, and the trace stopped at once. Nothing is loosened: the tolerance of N8 and N13
  is unchanged, and a mismatch beyond it is still an error.
  *Tests* (`tests/way_back_along_input_edges.rs`): strips along a straight line of input edges,
  flat and folded along the line (`z = fold · |y|`, mirror-symmetric as on bob), in both
  orientations, the line inside the strip or on its boundary; vertices inserted on the line
  edges and inside the faces, then flips; segments between two points where sides of a face
  cross the same line edge, from a corner on the line (an inserted vertex on it, or an input
  vertex) to such a point and back to the corner, and along sides that lie on the line (between
  two points of the side, to and from its corners). Before the fix 18 of 142 segments of the
  first interior strip failed with the error seen on bob, and 3 of 78 on the first boundary
  strip; after it all pass, every end is a point of the line edge, every piece lies in the
  input face it names, the mismatch is below `1e-12` of the strip's cells (below `5.2e-10` of
  the segment for segments of `4e-4` in thin faces), and point queries of crossings of input
  edges return points of those edges. Both reproducers come back (ka4: 588 segments, largest
  relative mismatch `1.4e-11`; ka2: 283 segments, `7.4e-12`; the study's other runs on bob's
  control mesh went to `8e-10`), on the triangulations they failed on (checked with the
  snapping switched off in `point_on_input` alone, for the build). For the callers: a node a
  rounding error from an input edge used to come back as a point of one face and now comes back
  as a point of the edge; every piece still lies in the face it names. `point_on_input` changed
  the same way, and the remeshing study reads the curvature of the input face under the centroid
  of a face from it: on bob's control mesh centroids on the mirror plane now get the mean of
  both faces, and the `ka4` build has 12 452 vertices instead of 12 454 (the 216 cuts of the
  control meshes rerun with ka2, ka4, u2, u4, odt2, odt4 all come back, largest relative
  mismatch `8e-10` as before, 1 to 10 fewer pieces per cut where nodes lay on input edges).
  *Still possible:* a point within `1e-12` of an input vertex without being at it is moved onto
  the nearest edge at that vertex only; a path that arrives at it from a face around the vertex
  that does not have that edge does not meet it.

Entry N50 was added with the second representation of intrinsic triangulations, integer
coordinates (`crates/ddg-integer-coords/SPEC.md`, decisions I1, I9, I13; work package WP0), and
FlipOut (`crates/ddg-flip-geodesics/SPEC.md`, requests R1, R3, R4). The numbers N30 to N49 are
those of `crates/ddg-intrinsic-apps/SPEC.md`.

* **N50** *Connectivity and representation, and the shared interface* (2026-10-31). The
  triangulation is split in two. `complex::Complex` is the connectivity alone: the halfedge
  Δ-complex of N9 (faces, twins, the edge of every halfedge, one halfedge per edge and per
  vertex) and a mask of *constrained* edges, with its topological operations: `flip_edge`,
  `split_face`, `split_edge`, and `remove_vertex` for an interior vertex of degree three (new:
  integer coordinates need removal, IC §4.6). Each returns its index bookkeeping (`ComplexFlip`,
  `FaceSplit`, `EdgeSplit`, `VertexRemoval`): the new vertex, edges and faces, and every
  halfedge whose data moved (a halfedge *moved* when it leaves the same vertex in the same
  direction under another id), with helpers that carry a representation's per-halfedge data
  along. A flip and a split keep every vertex, edge and face id and the direction of every edge
  they do not flip or split (FlipOut F8); removal is the only renumbering operation: the
  freed ids are taken by the last vertex, the last three edges and the last two faces, and
  `VertexRemoval` lists every move (`apply_*_data`). `IntrinsicTriangulation` embeds a
  `Complex` and keeps only its own data (lengths, signposts, angle sums, places, `edge_on_input`);
  its flip, insertion and split call the topological operations and then update that data in
  the same order and with the same arithmetic as before. **Bitwise identical**:
  `tests/bitwise_pins.rs` pins hashes of every stored bit (connectivity, lengths, signposts,
  angle sums, places), the fingerprint of N15 and edge lengths as bit patterns after random
  flips, insertions, splits, Delaunay flipping, edge traces, the overlay, moves, refinement,
  export and the optimal Delaunay iteration on 8 generated shapes and on spot and camelbad,
  recorded on the code before the split (commit `3b4ad07`); they are unchanged, and so are the
  test suites of this crate, `ddg-intrinsic-apps`, `ddg-uv` and `ddg-py` (one match arm for the
  new refusal was added to `tests/oracle.rs`).
  *The interface.* `interface::Triangulation` is the trait both representations implement and
  their consumers program against (IC §2.1, FlipOut R1): `complex()`, `input()`,
  `edge_lengths()`, `vertex_angle_sum`, `vertex_location`, `can_flip` / `flip_edge` (the
  refusals of §4 and N5, decided from the lengths alone, plus `FlipRefusal::Constrained`),
  `insert_vertex`, `split_edge` (the lengths of N10; `Insertion::halves` gives both halves in
  order from the old tail, directed like the edge), `set_constrained`, `trace_edge_on_input`
  (`EdgeTrace`), and provided methods for lengths, corner angles, areas, the Delaunay test and
  Delaunay flipping. Delaunay flipping (N6) is written once over the trait and
  `IntrinsicTriangulation::flip_to_delaunay` is that function, so both representations flip the
  same edges in the same order given the same lengths (IC §4.3).
  *Constrained edges* (FlipOut R3, R4; IC §2.1). A flip refuses them; Delaunay flipping leaves
  them, counts them as Delaunay and reports their number (`DelaunayReport::n_constrained`; the
  other numbers of the report are over the interior edges that are not constrained), as it does
  boundary edges; a split of a constrained edge gives two constrained halves and new edges are
  not constrained; `remove_vertex` refuses a vertex with a constrained edge. Refinement treats
  boundary edges specially (encroachment, N1; a circumcentre trace stops at the boundary);
  doing the same for constrained edges is FlipOut R5, not done yet, so `refine`, `refine_with`
  and `optimal_delaunay` refuse to run while an edge is constrained
  (`IntrinsicError::ConstrainedEdges`, nothing changed).
  *Removal in this crate.* `IntrinsicTriangulation::remove_vertex` removes an inserted interior
  vertex of degree three: the outer halfedges keep their signposts and the outer edges their
  lengths, so nothing is computed. None of this crate's operations calls it (N10 stands: they
  never remove vertices). Inserting a vertex and removing it again restores the triangulation
  exactly (`tests/complex.rs`).
  *Tracing without signposts* (IC I9). `frames::AngleFrames` is a view of a connectivity with
  lengths whose signposts are computed from the corner angles when the tracing kernel asks: at a
  vertex, direction `0` is its outgoing halfedge (the boundary halfedge at a boundary vertex, as
  N2), and the `k`-th halfedge counter-clockwise has `s (θ_0 + … + θ_{k−1})` with `s` of
  `it-angle-normalization` and `Θ` the computed angle sum (a running sum scaled once, as for the
  input mesh, so no reduction of N18 is needed). A frame belongs to the current connectivity (a
  flip at an interior vertex can change its outgoing halfedge). Measured on the 17 standard
  shapes after random flips, insertions and splits (`tests/angle_frames.rs`): computed
  signposts equal the stored ones, in the stored frame turned to put the outgoing halfedge at
  `0`, within `6.2e-12` rad (largest on the needle tetrahedron, whose stored input angle sums
  and the corner angles of its flipped needles differ by `1.2e-11` rad; below `4e-12` on the
  others, `0` before any operation); 6800 traces from vertices and points of faces end within
  `6.2e-12` of the mean edge length of the stored-signpost traces, about half of them bitwise
  at the same point, all with the same stop reason.
  *Mollification* (IC I13) moved to `ddg_core::mollify`; `ddg_bff::mollify_edge_lengths` wraps
  it with `BffError`.
  *Errors of other representations.* The trait's methods return `IntrinsicError`; another
  implementation passes its own failures as `IntrinsicError::Representation { representation,
  operation, message }` and an operation it does not provide (yet) as `Unsupported {
  representation, operation, reason }`, so the provided methods (Delaunay flipping) run over any
  implementation.
  *Cost.* `benches/retriangulation.rs` (iDT, iDR 25°, iODT, a sweep of moves on four meshes),
  run single-threaded on one pinned core, the code before the split and after it alternately
  twice (load average 3 to 18 from other jobs): the differences are within ±6 % in both
  directions and change sign between the two rounds (e.g. spot iDT +2.1 %, then −6.0 %;
  camelbad iDT +3.2 %, then −0.6 %); no slowdown is measurable above that noise.

* **N51** *Refinement over the shared interface, and the rules of the integer-coordinates paper*
  (2026-10-31, IC-SPEC §5, I11, I23). The driver of N1, N4, N14 is written once over
  `Triangulation`: `refine::refine(&mut t, &RefineOptions, RefineRules, extra)`;
  `IntrinsicTriangulation::refine` and `refine_with` are it with `RefineRules::Ruppert`, and the
  provided `Triangulation::refine_with_rules` chooses the rules. `RefineOptions` is unchanged (its
  callers list every field); the rules are a separate argument.
  *What the trait gained* (all provided, so other implementations compile unchanged):
  `trace_geodesic(start, direction, length, options)` (default: `AngleFrames`; this crate
  overrides it with its stored signposts; from a point of a face the path, its crossings and its
  end depend on the lengths alone, so both agree bit for bit), `remove_inserted_vertex(v) ->
  InsertedVertexRemoval` (default `Unsupported`; here: flips by `removal_flip_choice` until the
  degree is three, then `remove_vertex`; a star where no flip lowers the degree is
  `RemovalRefusal::Stuck { degree }`, the triangulation restored from a copy), and
  `refine_with_rules`. `remove::removal_flip_choice(&Complex, &[f64], v)` is IC I24 (self-edges
  first, then the largest `min(A₁, A₂)/(A₁ + A₂)`, ties by edge id), decided from the lengths alone
  and used by both representations, so they remove a vertex by the same flips.
  `EdgeQueue::apply_removal` keeps the Delaunay queue valid across the renumbering of a removal;
  the refinement heap skips face ids past the end and re-stamps faces that a removal moved.
  *The two rules.* `RefineRules::Ruppert` (default; N1, N4): encroached boundary edges (opposite
  angle `> 90°`) are split first; a circumcentre that would encroach a boundary edge of the
  triangles its insertion removes splits that edge instead; boundary edges next to a boundary
  angle below 60° are split at power-of-two distances, others at the midpoint; corners at vertices
  whose angle sum is below the bound are exempt; the needle guard of N4 (a circumcentre whose loop
  around an interior vertex of angle sum below 60° would be shorter than the triangle's shortest
  edge is not inserted); vertices are never removed. `RefineRules::Chew` (Gillespie, Sharp, Crane
  2021, Supp. Alg. 10, as their reference's `delaunayRefine`): no edge is split for being
  encroached; when the circumcentre path from the barycentre leaves the surface through a
  boundary edge `lm`, `lm` is split at its midpoint, the triangulation is flipped to Delaunay,
  every inserted interior vertex at Dijkstra distance (over the edges) `< ℓ_lm` from the new
  vertex is removed in decreasing id order (`ic-dijkstra-ball`; the reference uses `≤`), and it is
  flipped to Delaunay again; a refused removal is counted, not fatal; a circumcentre on an
  interior edge splits it at the foot as under N14, without removals (the reference removes after
  every split); the narrow-cone rule of IC I23: a face with exactly one corner at a vertex of angle
  sum below 60°, or whose vertices' places all lie in one input face with such a vertex (the
  reference's `getParentFace`, from the places alone), is not refined; a corner at a vertex of
  degree one or between two boundary edges does not count as small. Under both rules a
  circumcentre path that ends elsewhere than in a face or at the boundary skips the face and is
  counted (under `Ruppert` it was skipped silently; an error of the trace itself stays an error
  there, and is counted under `Chew`). Order, tie-breaking and the insertion limit are N14's under
  both.
  *Report.* `RefineReport` gained `rules`, `n_removals`, `n_removal_flips` (flips made by
  removals, not in `n_flips`), `n_removals_refused`, `n_trace_failures`, `n_narrow_cone_faces` (the
  faces of the result the I23 rule leaves, whichever rules ran) and `min_angle_valid` (the
  smallest corner angle over the other faces: the paper's "minimum corner angle of the output
  mesh"; the reference's statistic leaves out faces with *any* narrow vertex). Under `Chew` a
  corner below the bound in a face the rule leaves is listed as `BelowReason::NeedleVertex` with
  the narrow vertex (of the triangulation or of the input face). The existing behaviour that a
  run whose insertion count reaches the limit exactly reports `InsertionLimit`, not `Done`, is
  unchanged.
  *Bitwise.* With `Ruppert` the results are those before the change: `tests/bitwise_pins.rs` and
  the suites of this crate, `ddg-intrinsic-apps`, `ddg-uv` and `ddg-py` pass unchanged.
  *Measured* (`ddg-integer-coords/tests/refine.rs`, 256 runs on the generator shapes, spot,
  camelbad, cones with small angle sums, meshes with boundary and Δ-complex starts, both rules, 25°
  and 30°): this crate and the integer coordinates make the same operations and reach the same
  connectivity with bitwise the same lengths in 255 runs; the other one parts where a
  circumcentre lies on an interior edge exactly at an input-edge crossing (the integer split
  refuses it, IC I20; this crate snaps the place onto the input edge). Under `Chew` 2 307 vertices
  were removed and 175 removals refused, 0 trace failures.
  *Cost* (`benches/retriangulation.rs`, `idr25`, the code before the change and after it
  alternately, one pinned core, load average 0.2 to 3.4 from other jobs, five rounds): within
  −2.8 % to +3.1 % on the noisy icosphere, the sheared grid and spot, +1.8 % to +3 % on camelbad
  (one outlier +6 %). Instruction counts (`perf stat`, 20 refinements): +0.5 % on spot, +0.4 % on
  camelbad, of which the end-of-run narrow-cone statistic was about two thirds; it is now skipped
  when no vertex is narrow (spot and camelbad have none; 97.2 % of Thingi10k, paper Sec. 4.2) and
  otherwise evaluated only for faces with a corner it can involve. The remaining ~0.25 M
  instructions per refinement (spot) are spread over the generic driver; the last two rounds
  measure −0.9 % to +1.2 % on spot and +2 % on camelbad.

Entry N52 was added with FlipOut's requests R9 and R10 (`crates/ddg-flip-geodesics/SPEC.md` §11,
§7 test 9).

* **N52** *The shape test of a flip, chosen by the caller, and the corner angles a flip returns on request*
  (2026-10-31, user decision on R10 and R9).
  *The test.* `FlipShapeTest` is the test a flip's two new triangles must pass after strict
  convexity (the signs of `it-flip-convexity`) and before the triangle inequalities of both new
  triangles, which stay under every test. With `A_tail`, `A_head` twice the areas of the new
  triangles `ilk` and `jkl` (`Diamond::flipped_areas`): `Absolute { limit }` is N5, refused
  unless `A_T ≥ limit · ℓ_new · ℓ_max(T)` for both; `FlipShapeTest::N5 = Absolute { limit:
  FLIP_SHAPE_LIMIT }` is the default. `Relative { ratio }` is the reference's
  (`flipEdgeIfPossible`, geometry-central `signpost_intrinsic_triangulation.cpp:333-338`), refused
  unless `min(A_tail, A_head) ≥ ratio · (A_tail + A_head)`; `FlipShapeTest::REFERENCE` has its
  `ratio = 1e-6`. Both are invariant under scaling. A flip the test refuses is
  `FlipRefusal::TooThin` with N5's shape number, whichever test refused it. The two differ where
  a new triangle is a needle whose short side is not the new edge and the two areas are
  comparable (the relative test accepts it; on `tree` shapes `1.2e-7` to `6.9e-7` with ratios
  `5e-5` to `1e-3`), and where one new triangle is far smaller than the other but not thin (the
  relative test refuses it: camelbad 8294 → 9415 of FlipOut's 9 × 500 pairs has one refused flip
  under it and none under N5). By N5's argument the accepted
  needles have angles that are uncertain by `ε / shape ≈ 1e-9` rad, which is what FlipOut's
  certificate tolerates (`1e-4`).
  *Interface.* `plan_flip_with(complex, lengths, e, test)`; `plan_flip` is it with
  `FlipShapeTest::N5`. `can_flip_with(e, test)` takes the test; `flip_edge_with(e, request)` takes
  a `FlipRequest { shape: FlipShapeTest, corner_angles: bool }` (`FlipRequest::DEFAULT` = N5, no
  angles; it is `Default`), on `IntrinsicTriangulation` and in `ddg-integer-coords`;
  `Triangulation::{flip_edge_with, can_flip_with}` are provided (the default test, and for
  `flip_edge_with` the default request, go to `can_flip` / `flip_edge`, anything else is
  `IntrinsicError::Unsupported`, so wrappers that implement the trait compile unchanged) and
  overridden by both representations. `flip_edge` and `can_flip` keep N5 and return no angles. *Who uses which:* every operation of this crate and of the
  integer coordinates (Delaunay flipping, refinement, removal, insertion, the optimal Delaunay
  iteration, the export) keeps N5; only FlipOut's Alg. 1 asks for `REFERENCE`, by default
  (`ShortenOptions::flip_shape`, which can be set back to `N5`); FlipOut's iCDT is Delaunay
  flipping and keeps N5.
  *Corner angles* (R9), on request. `Flip::new_corner_angles` is an `Option<[f64; 2]>`: `None`
  unless `FlipRequest::corner_angles` was set, then `Some([θ_l, θ_k])`, in the order of
  `new_vertices = [l, k]`: the angle at `l` in the new face `(j, k, l)` and the angle at `k` in
  the new face `(i, l, k)`. Both representations fill them if and only if they are asked, so the
  same request gives the same `Flip` from both (`ddg-integer-coords/tests/ops.rs` compares them
  for plain `flip_edge` and for requests with N5 and with the reference's test).
  `Flip::new_corners` (their halfedges after the flip, `l → j` and `k → i`, the two that stay in
  place, `ComplexFlip`) is always given. Each angle is bitwise `Triangulation::corner_angle` of
  that corner (`ddg_core::geometry::corner_angle` orders its two sides, so the argument order
  does not matter). The signposts compute exactly these two for their signpost update, so asking
  costs them nothing; the integer coordinates compute them from the plan
  (`FlipPlan::new_corner_angles`, two `corner_angle`), which costs about 60 ns per flip, so only
  a caller that uses them asks. *Who asks:* FlipOut's Alg. 1 (its memo takes them); Delaunay
  flipping, refinement, removal, insertion, the optimal Delaunay iteration and iCDT do not. The
  other two corners at `l` and `k` are not computed.
  *Bitwise.* With the default test nothing changes: `tests/bitwise_pins.rs` and the suites of
  this crate, `ddg-intrinsic-apps`, `ddg-uv`, `ddg-py` and `ddg-integer-coords` pass unchanged.
  FlipOut's fingerprints (9 meshes × 500 pairs × 2 representations, `examples/fo_pin`) with
  `flip_shape = N5` are those from before; with `REFERENCE` 99 pairs change on both
  representations alike, all 4500 certify; using the returned angles in FlipOut's memo changes
  no fingerprint.
  *Cost* (one pinned core, load average 1.7 to 8.7 from other jobs). Computing the angles in every
  flip, as the first version of N52 did, cost the integer coordinates +5.7 % in camelbad's Delaunay
  flipping (8110 flips, 8.84 → 9.35 ms) and +3 % to +7 % in FlipOut. With the request: Delaunay
  flipping from the input mesh (best of 3 interleaved rounds, integers; signposts unchanged within
  1 %): camelbad 9.47 → 8.89 ms (−6 %; before the angles existed 8.99), spot (281 flips) 1.58 →
  1.54 (1.57), tree (594 flips) 8.91 → 8.71 (8.69). FlipOut's shortening (median per pair, 9
  meshes × 500 pairs, best of 4 interleaved rounds, N5 on both sides, against the code before the
  angles existed, `adf9d80^`): signposts −5.8 % to −7.5 % (the memo takes angles the signpost update
  computed anyway); integers, which compute them for the memo, +2.4 % to +8.2 % (+0.7 % to +3.2 %
  when FlipOut does not ask), so the two angles cost an integer flip more than the memo saves.
  Against the first version of N52 the integers are 2 % to 4 % faster.

Formula ids in use: `it-angle-normalization`, `it-signpost-update`, `it-unfold-diamond`,
`it-flip-length`, `it-flip-convexity`, `it-delaunay-condition`, `it-cotan-weight`,
`it-unfold-apex`, `it-ray-edge-intersection`, `it-barycentric-coordinates`,
`it-point-distances`, `it-circumcenter`, `it-move-length`, `it-move-area`,
`it-circumcenter-vector`, `it-odt-update`; from `ddg-core`: `cotan-laplacian`, `lumped-mass`,
`angle-defect`, `boundary-geodesic-curvature`. The lengths of an edge split are
`it-point-distances` with the point on a side. The Lean theorems of §7 are in
`lean/proofs/DdgProofs/Intrinsic/`, the citations in `lean/proofs/CITED.md`, and the oracle's
single operations are replayed by `tests/oracle.rs`.
