# FlipOut geodesics — contract

The contract between the Rust crate `crates/ddg-flip-geodesics`, the Lean oracle (`lean/oracle`),
the Lean proofs (`lean/proofs`), the Python binding (`crates/ddg-py`) and the Blender add-on for

> Nicholas Sharp, Keenan Crane. *You Can Find Geodesic Paths in Triangle Meshes by Just Flipping
> Edges.* ACM Transactions on Graphics 39(6), Article 249, 2020.
> <https://doi.org/10.1145/3414685.3417839>
> (`../FlipOut/FlipOut.pdf`, "the paper", 15 pages with Appendices A and B; equation, algorithm,
> theorem and figure numbers are the paper's),

run on either intrinsic triangulation of this workspace, and for Sec. 5.2 and Fig. 17 of the
integer-coordinates paper (Gillespie, Sharp, Crane 2021, `../IntegerCoordinates/integer_coordinates.pdf`,
"IntCoords"), which runs FlipOut on the integer representation.

**Reference.** geometry-central at commit `ab116f1` (MIT, Copyright (c) 2017-2019 Nicholas Sharp
and the geometry-central contributors) in
`../IntegerCoordinates/intrinsic-triangulations-demo/deps/geometry-central/`, cited as
`gc: <file>:<line>` with `fg.cpp` = `src/surface/flip_geodesics.cpp`, `fg.h` =
`include/geometrycentral/surface/flip_geodesics.h`, `ic.cpp` =
`src/surface/integer_coordinates_intrinsic_triangulation.cpp`, `it.cpp`/`it.ipp` =
`src/surface/intrinsic_triangulation.cpp`/`include/…/intrinsic_triangulation.ipp`, `mga.cpp` =
`src/surface/mesh_graph_algorithms.cpp`. May be read and ported; ported code keeps its notice in
`LICENSES/`; the README credits both authors. Two facts about the reference shape this contract:
`FlipEdgeNetwork` is hard-wired to `SignpostIntrinsicTriangulation` (gc: fg.h:127, fg.cpp:224), and
`IntegerCoordinatesIntrinsicTriangulation::splitEdge(Halfedge, double)` throws "not implemented"
(gc: ic.cpp:1246-1250). So the IntCoords Sec. 5.2 runs (paths, loops, Bézier curves on integer
coordinates) used code that is not at this commit; the integer half of this crate is written from
the papers and the two contracts below.

**Other contracts.** Signposts: `crates/ddg-intrinsic/SPEC.md` ("IT-SPEC", decisions `N<n>`).
Integer coordinates: `crates/ddg-integer-coords/SPEC.md` ("IC-SPEC", decisions `I<n>`), whose
shared Δ-complex connectivity (I1, WP0) this crate builds on. Everything they say about the input
mesh, halfedge conventions (IT N9: `h = 3f + k`, stable vertex/edge/face ids), flips (IT §4, N5)
and edge traces (IT N7, N8) holds here. Decisions of this file have ids `F<n>` (§10); what this
file needs changed in the other two is §11.

Conventions of `CLAUDE.md` apply: `f64`, stored indices `u32`, explicit option structs, a report
from every operation that does a variable amount of work, bad input is an error naming the
element and never a panic, nothing repaired or normalized silently, deterministic results at any
thread count, formula tags `[formula: fo-…]` matched by Lean theorems or `lean/proofs/CITED.md`.

## 1. The problem and scope

### 1.1 Problem

`T⁰` is the input mesh (a manifold `TriMesh` with edge lengths), `T` an intrinsic triangulation of
the same polyhedral surface (IT §1; initially `T = T⁰`). An **edge network** `Γ` on `T` is a set of
edge paths (open paths and closed loops) with **pinned** vertices (paper Sec. 3.2.1, 6.1). FlipOut
changes `T` by intrinsic edge flips (never moving the surface, never inserting a vertex) and
replaces sub-paths by shorter edge paths of the new `T` until every joint is locally shortest
(Eq. 2) or blocked by a pinned vertex, keeping `Γ` in its isotopy class (non-crossing curves stay
non-crossing, Sec. 3.2.2). The result is extracted as polylines on `T⁰`.

What must be true, in one sentence: the returned paths are edge paths of the returned `T` that are
isotopic to the input relative to the pinned vertices, no longer than the input, and geodesic
(every interior joint not at a pinned vertex has, on each side, a wedge angle `≥ π − tol` or a
wedge blocked by a path end; F27), and the polyline extracted on `T⁰` is that same curve: connected
exactly (integer coordinates) or within the trace check of IT N8 (signposts).

### 1.2 Scope

In: paths, loops and networks (Sec. 4.2, 4.3, 6.1), multiply covered curves (Fig. 8, 11), initial
paths by Dijkstra on the edge graph (Sec. 5.1; gc: mga.cpp:14-84), through waypoints
(gc: fg.cpp:283-312) and from an edge set with junctions (gc: fg.cpp:315-396), curve shortening
with early stopping (Sec. 4.2 end, Fig. 20, 21), geodesic Bézier curves (Sec. 6.2, Fig. 23;
IntCoords Fig. 17 right), the intrinsic constrained Delaunay triangulation (Sec. 6.3), extraction
on `T⁰` for both representations, the Sec. 5 evaluation and the IntCoords Sec. 5.2 benchmark.
Later (WP F9, §12): Delaunay refinement conforming to the paths (Sec. 6.4). Out (F21): the
single-source geodesic tree and log map (Sec. 6.5), the Lagrangian baselines of Fig. 17/18,
PDE boundary conditions and aligned cross fields (Sec. 6.4 examples, Fig. 26).

### 1.3 Two representations, one algorithm

FlipOut touches the triangulation only through the operations below (F1). Both
`ddg_intrinsic::IntrinsicTriangulation` and `ddg_integer_coords::IntegerTriangulation` provide them
through the representation-neutral interface of the shared layer (IC I1/WP0, extended by §11 R1):

| operation | used for | signposts (`ddg-intrinsic`) | integer coordinates (`ddg-integer-coords`) |
|---|---|---|---|
| connectivity: `next`, `twin`, `edge`, `face`, `from_vertex`, `edge_halfedge`, boundary | wedges, outer arcs | shared Δ-complex (IC I1) | same type |
| `edge_length(e)`, corner angle from lengths (`corner-angle`) | Eq. 1, `β_i` | stored lengths | stored lengths |
| `can_flip(e)` / `flip_edge(e)` with IT §4 + N5 refusals | Alg. 1 | Alg. 5 + signpost update | IC §4.2 (Eq. 3, 4) |
| `split_edge(e, t)` returning both halves in `e`'s direction | Bézier (§6.1), refinement (later) | IT N10 | IC §4.5 (§11 R4) |
| constrained-edge mask honoured by Delaunay flipping | iCDT (§6.3) | §11 R3 | §11 R3 |
| `trace_edge_on_input(e) -> ddg_intrinsic::EdgeTrace` | extraction (§5) | `trace_edge` (IT N7, N8) | `ExtractIntrinsicEdge` (§11 R2) |
| `vertex_location(v) -> SurfacePoint on T⁰` | ends of polylines | IT §3 | `V⁰`: itself; `V*`: `q⁰` (IC §2.2) |

Every lengths-producing operation uses the same formula in both (IT `it-flip-length`,
`it-point-distances`; IC §4.3 already requires bitwise-identical Delaunay flipping). Hence, on the
same input, options and queue, **FlipOut makes the same flips and returns the same edge paths with
bitwise-identical lengths on both representations** (§7 item 8). They differ only in extraction:

* **Signposts.** Per segment, `trace_edge` of its edge (reversed for a backward segment): the
  floating-point walk over `T⁰` from the tail's signpost, checked by IT N8 (it must reach the far
  vertex within `1e-6 ℓ`). A failed check is an error naming the edge (F18); the polyline is then
  not connected. This is the failure of IntCoords Fig. 17 left.
* **Integer coordinates.** Per segment, `ExtractIntrinsicEdge` (§11 R2, the dual of IC Alg. 4):
  the `n⁺_e` crossings of `e` with `T⁰` in order `p = 0 … n_e − 1`; the input edge of every
  crossing identified from the integers alone (TraceFrom + IdentifyTracedEdge, IC Alg. 1, 3);
  consecutive crossed input edges checked to bound a common input face; positions `t` from the
  straight segment through the unfolded strip of those input faces. Connectivity is exact; floating
  point decides only positions (clamped and counted, IC I5). `n_e = −1`: the segment runs along an
  input edge (its part between the segment's ends).

In both, the path's polyline is the concatenation of its segments' traces, the shared vertex at a
joint written once (it is the same `SurfacePoint`, so the reference's `onSameElement` test,
gc: fg.cpp:1588-1591, is not needed).

## 2. Data

### 2.1 Network

* **Segment** `(edge: u32, forward: bool)`: a directed edge of `T`, `forward` relative to
  `edge_halfedge(edge)`. Stored by edge id, not halfedge id (halfedge ids move under flips, IT N9);
  the orientation is stable because path edges are never flipped and an unflipped edge keeps its
  direction when a neighbour is flipped (ddg-intrinsic `flip.rs:225-231`; required of the shared
  type, F8).
* **Path**: `segments: Vec<Segment>` (consecutive segments share a vertex), `closed: bool`, explicit
  (F9). Joints: `γ_abc = (ab, bc)` for consecutive segments, and last→first if closed. A path of
  one segment that is closed is a single-segment loop (§4.3).
* **Pinned vertices** `P ⊆ V`: no joint at a pinned vertex is ever shortened. The two ends of an
  open path are pinned implicitly. (The paper's "marked vertices", Sec. 6.1; F10.)
* **Layering** (Fig. 11): per edge, the segments on it ordered across the edge, from the face of
  `edge_halfedge(e)` to the face of its twin (gc: fg.h:133, a deque). "Outermost on side `h`" is the
  first entry if `h = edge_halfedge(e)`, the last otherwise (gc: fg.cpp:1319-1351).
* **Validity** (checked on construction and, with `check_invariants`, after every iteration): paths
  connected; open-path ends pinned; every segment appears exactly once in the layering of its edge;
  **non-crossing**: at every vertex `b`, list the segment ends at `b` in the cyclic order of the
  edges around `b`, each edge's segments in their layer order; every joint at `b` pairs two of these
  slots, and no two pairs interleave (`fo-non-crossing`). Path ends at `b` pair nothing and
  constrain nothing. A violation is `FlipError::Crossing { vertex, paths }`.

### 2.2 Initial paths

* `dijkstra_path(tri, from, to)`: shortest path in the edge graph of the current `T`, weights `ℓ`,
  ties by `(distance, halfedge id)` (F15; gc: mga.cpp:29-30 ties by halfedge). Same vertex:
  `SameVertex`; not connected: `NotConnected { from, to }` (gc returns an empty path,
  fg.cpp:275-278).
* `waypoint_path(tri, points, closed, pin_waypoints)`: Dijkstra between consecutive waypoints,
  concatenated (gc: fg.cpp:283-312, `markInterior`).
* `network_from_edges(tri, edge_mask, extra_pinned)`: the reference's walk (gc: fg.cpp:315-396):
  vertices of path degree `≠ 2` and `extra_pinned` are ends (pinned); paths between ends; the rest
  closed loops. Used for seams and cut graphs (Fig. 19, 22).
* Explicit: `GeodesicNetwork::new(tri, paths, pinned, layering: Option<…>)`. Without a layering,
  each segment is pushed outermost on the side of its own halfedge in path order (gc:
  fg.cpp:29-32, F11), then the non-crossing check runs; inputs that cross are refused, never
  repaired. Input paths need not lie on `T⁰`'s edges: they are edge paths of `T`.

### 2.3 Output

`GeodesicNetwork<T>` owns `T` and `Γ`; `paths()`, `length()` (running sum, F14),
`certify(tol) -> GeodesicCertificate` (§7 layer 1), `extract(&ExtractOptions) -> Vec<PathOnInput>`
(`ddg_intrinsic::PathOnInput`: points on `T⁰`, the input face of every piece, mismatch),
`polylines_3d(positions)` (barycentric interpolation, for Python/Blender), `triangulation()` /
`into_triangulation()` (the conforming `T`, Fig. 1 bottom), and the iCDT (§6.3).

## 3. Geometry of a joint (paper Sec. 3.2.3–3.2.5)

For a joint `γ_abc` at `b`, the two **wedges** are the corners at `b` swept counter-clockwise from
`bc` to `ba` (`W_abc`, left of the path) and from `ba` to `bc` (right), with counter-clockwise
faces (the reference's `leftAngle = angleIn − angleOut`, gc: fg.cpp:489-521; tested geometrically).
Corner angles are summed counter-clockwise from the start of the sweep, independent of
`outgoing_halfedge`, so repeated evaluations are bitwise equal (the oracle sums in the same order).
**Out-and-back joint** (F26): when `ab` and `bc` lie on the same edge, the layering decides: the
side between the two segments has no corner and angle `0`, the other side has `Θ_b`. With corner angles `θ` from
the edge lengths (`ddg_core::geometry::corner_angle`):

```
α(W) = Σ_{corners jk_b ∈ W} θ^{jk}_b                                       [formula: fo-wedge-angle]   (Eq. 1)
α(W) = +∞  if the sweep crosses the boundary gap at b                      [formula: fo-boundary-wedge] (Sec. 3.2.3)
locally shortest:  α^min_b ≥ π                                             [formula: fo-locally-shortest] (Eq. 2)
Θ_b = α^min_b + α^max_b,  so α^min_b ≥ π ⇒ Θ_b ≥ 2π ⇒ Ω_b ≤ 0               [formula: fo-geodesic-saddle] (Eq. 3)
```

Both sides are summed corner by corner (F2), so `α_L + α_R` is the current angle sum to rounding;
the reference subtracts signposts instead (gc: fg.cpp:489-522), which in `ddg-intrinsic` are stored
scaled by `2π/Θ` (IT N2) and would need unscaling.

**Classification** (gc: fg.cpp:524-587): `Straight` if `α^min > π − straight_tolerance` (default
`1e-5`, gc: fg.h:161; F3), otherwise `Left` if `α_L < α_R` strictly, else `Right` (the reference's
tie rule); the queue takes sides with `α ≤ π − tol`, the certificate accepts `α ≥ π − tol`; the larger side is also
queued when it is `< π − tol` too (`Θ_b < 2π`, gc: fg.cpp:1005-1007).

**Flexible** (Sec. 3.2.3; gc: fg.cpp:589-658): the joint's own segments `ab`, `bc` are outermost on
the wedge side of their edges, and no edge strictly inside the wedge carries a segment. **Blocked
by an end** (gc: fg.cpp:662-730): some segment in the wedge, or outermost on a bounding edge, is
the first or last segment of an open path; such a joint is locally shortest in its isotopy class
(Sec. 4.2, case (b)) and the certificate exempts it. A path end at `b` itself pins `b`, so the
exemption matters only for end segments whose far vertex is the end. **Flexible** in code: no
segment end strictly inside the sweep, which also covers the degree-0 wedge and self-edges.

**Vertex cases** (what the definitions imply; tested in §7):

| vertex `b` | consequence |
|---|---|
| interior, `Θ_b < π` (incl. degree 1, a single corner) | both wedges `< π`: every joint through `b` is shortened, smaller side first |
| interior, `π ≤ Θ_b < 2π` | `α^min ≤ Θ_b/2 < π`: no geodesic passes through `b` (cone point) |
| interior, `Θ_b = 2π` (flat; every inserted vertex) | passes only if both sides are exactly `π` (non-generic, Sec. 3.2.5) |
| interior, `Θ_b > 2π` (saddle) | a joint with both sides `≥ π` stays: geodesics may bend there |
| boundary, interior angle `Θ_b < π` | the boundary side is `∞`, the other `< π`: the path is pushed into the interior |
| boundary, `Θ_b ≥ π` (reflex corner) | the path may wrap around it and stay (geodesic touching the boundary) |
| pinned | never shortened; the path may bend there (waypoints, control points, network junctions) |

## 4. The algorithm

### 4.1 Local shortening, Alg. 1 `FlipOut(T, γ_abc)`

Precondition: `γ_abc` flexible, `α_abc < π`, `ab ≠ bc` as segments. Unfold the wedge (only
conceptually): outer arc `a = n₀, n₁, …, n_k = c`, `k = deg(W)`, and at `n_i` (`0 < i < k`)

```
β_i = θ^{b n_{i+1}}_{n_i} + θ^{n_{i−1} b}_{n_i}   (the two corners at n_i of the wedge's triangles)   [formula: fo-outer-arc-angle]
```

1. While some `β_i < π`: `j = min{i : β_i < π}`; flip `b n_j` (Alg. 1 line 2–3; App. A explains why
   the *first* such edge: it avoids degree-1 vertices, Thm A.5, and cycling, Thm A.8). After a flip
   only `β_{j−1}` can newly be `< π`, so the scan resumes at `j − 1` (Sec. 4.1, gc: fg.cpp:814).
2. `γ_shorter = (a, n₁, …, n_{k−1}, c)`; replace the joint by it (gc: fg.cpp:79-183): the new
   segments go innermost (wedge side) on the outer-arc edges; joints at `a`, `c` and the outermost
   segments that shared `ab`, `bc` are re-queued (gc: fg.cpp:162-182).
3. `deg(W) = 1`: `γ_shorter = (a, c)` along the single opposite edge. `deg(W) = 0` (out and back
   along one edge, `α = 0`): the joint contracts, both segments are removed (Sec. 4.1, F12).

Why every flip is valid: the diamond of `b n_j` has angles `β_j < π`, two triangle corners `< π`,
and the inner corner `< α < π`, so it is strictly convex (`fo-wedge-flip-convex`, Thm A.6); the
degree condition holds by Lemmas A.1–A.4 and Thm A.5 (`fo-wedge-flip-valid`). Termination: `k`
never increases, and each flip lowers `k` or the index `ρ` of the first twice-incident edge
(Thm A.7, A.8; `fo-flipout-terminates`; at most `k₀²` flips). Result: `|γ_shorter| < |γ_abc|`
(Thm 4.1; `fo-outer-arc-shorter`).

**Floating point** (F6). A flip refused by IT N5 (`TooThin`) although `β_j < π − 0` is skipped and
the scan continues past it, as the reference's `flipEdgeIfPossible` loop does (gc: fg.cpp:785-819);
counted (`refused_flips`). Before replacing, the new length must be strictly smaller than
`ℓ_ab + ℓ_bc` (gc: fg.cpp:836-838); otherwise the path is left as it was, the flips are kept, and
the joint is counted `stalled` and not re-queued until a neighbour changes. An edge carrying a
segment is never flipped (F7; the reference's skip at gc: fg.cpp:787-790 is here an invariant
check returning `FlipError::Internal`).

### 4.2 Geodesic paths, Alg. 2 `MakePathGeodesic`

A min-priority queue of joints keyed `(α, side, path index, segment id)` with `f64::total_cmp`
(F4; the reference ties by pointer, gc: fg.h:147-148, fg.cpp:1813-1830). Initially every joint of
every path (both sides when both are `< π − tol`). Loop (gc: fg.cpp:919-975):

1. Pop. **Stale** if the segment is gone or its recomputed angle differs bitwise from the key
   (gc: fg.cpp:937-945): skip. A joint at a pinned vertex: skip.
2. Not flexible: drop (gc: fg.cpp:950-952). It is re-queued when a neighbour's replacement makes
   it outermost (step 2 of §4.1).
3. Alg. 1, then update. Iteration counted.
4. Stop rules (F14): `max_iterations` (status `IterationLimit`), `stop_length_ratio` (current/initial
   `<` ratio, the paper's 0.9 for segmentations, Sec. 6.1; gc: fg.cpp:967-973), `stop_min_angle`
   (all queued `α ≥` threshold, Sec. 4.2 end).
5. Queue empty: a **final sweep** recomputes every joint and re-queues each flexible one that is
   not straight (F5, a departure: the reference stops at an empty queue); at most `max_sweeps`
   (default 4) sweeps that find work, reported. Status `Geodesic` when a sweep finds none.

The queue is rebuilt when stale entries outnumber live ones (gc purges every 1000 iterations,
fg.cpp:962-965); this never changes which joint is processed next.

### 4.3 Loops (Sec. 4.3, App. B)

Closed paths use Alg. 2 unchanged, except the **single-segment loop** `(aa)` with `α < π` on a side:
replace it by the two other edges of the triangle on that side (`ab`, `ba`), an isotopy inside one
face (gc: fg.cpp:855-917). The length grows; App. B shows the next step does not undo it
(`φ₁ + φ₂ < ω`, Eq. 4–6; `fo-loop-perturbation`). A loop that contracts to a vertex (§4.1 step 3
leaves no segment) is removed and reported in `contracted` (Fig. 20 "contracts to point"; F12).
Termination for loops is a conjecture (App. B end; Sec. 7), so loops always run under
`max_iterations` (default `100 (|E| + m₀)`, `m₀` the initial number of segments; F13).

### 4.4 Networks (Sec. 6.1)

Several paths in one `T`: one queue, one layering; nothing else changes. Junctions are pinned
(they are path ends, §2.2), so the topology of the network is kept (Fig. 6, 19).

### 4.5 Guarantees: what is proved, what is cited

| statement | source | here |
|---|---|---|
| a joint with `α^min < π` is shortened by Alg. 1, every flip valid, finitely many flips | Thm 4.1, A.5–A.8 | Lean (§8, tiers B, C) |
| Alg. 2 on paths and networks terminates | Thm 4.2, via finiteness of geodesics of bounded length (Indermitte et al. 2001, Prop. 1) | the reduction proved (A), finiteness cited |
| Alg. 2 never runs out of flexible joints before the curve is geodesic in its isotopy class | Sec. 4.2 (argument in prose) | tier D attempt, else cited |
| isotopy class / non-crossing preserved | Sec. 4.2 | tier C (`fo-non-crossing`, `fo-isotopy`) |
| the result is a polyhedral geodesic (Mitchell et al. 1987, Lemma 3.5) | Eq. 2, 3 | certified per run (§7 layer 1) |
| loops terminate | conjecture (App. B, Sec. 7) | cited as open; `max_iterations` |
| minimality | not claimed (Sec. 6.5) | measured against exact distance (§7) |

In floating point, Alg. 2 stops with every joint straight to `straight_tolerance` or blocked, or
with a reported status; the certificate (§7) is computed on the final state, never assumed.

## 5. Extraction on the input mesh (paper Sec. 5.1; IntCoords Sec. 5.2)

`extract(options)` returns one `PathOnInput` per path, segments in path order, as described in
§1.3. Per edge the work is independent, so edges are extracted in parallel and assembled in order
(F24). Options: the representation's trace options (IT `EdgeTraceOptions`; IC `IntegerOptions`),
`certify_input: bool` (default on in tests and the evaluation).

**Certificate on `T⁰`** (independent of `T`, §7 layer 1b): along the polyline, at every crossing
of an input edge the two adjacent pieces, unfolded across that edge (`it-unfold-diamond`), make an
angle `π` within `input_tolerance` (`fo-straight-crossing`); at every polyline vertex that is an
input vertex, the two sides measured on `T⁰` (input corner angles plus the partial angles of the
first and last pieces) satisfy Eq. 2 within `input_tolerance`; the length of the polyline equals
the path's length on `T` within `1e-9` relative. Integer extraction also reports, per edge, the
distance between the unfolded ends and `ℓ_e` (`length_mismatch`, a consistency check between
`T` and `T⁰`).

## 6. Applications in scope

### 6.1 Geodesic Bézier curves (Sec. 6.2, Fig. 23; Morera et al. 2008; IntCoords Fig. 17 right)

Input: control vertices `P₀ … P_d` (`d ≥ 1`), rounds `r`. Initial path by `waypoint_path` with the
control points pinned. The reference's recursion (gc: fg.cpp:1068-1248), which is the paper's four
steps applied to each half:

1. Shorten (Alg. 2): every piece between consecutive pinned points is a geodesic.
2. In the active range, for each piece of length `L` insert its midpoint: walk its segments to the
   one containing arc length `L/2` and `split_edge` at `t = (L/2 − s)/ℓ` (`s` the length before it;
   `fo-bezier-midpoint-split`), pin the new vertex. If the point is within `1e-4 ℓ` of either end
   of that edge, pin that vertex instead (F19; the reference checks only the far end,
   gc: fg.cpp:1195-1229). Unpin the old interior control points of the active range except its
   first and last (gc: fg.cpp:1164-1183); shrink the range to the new points; repeat until one
   point remains: `B(½)`, the de Casteljau midpoint.
3. Recurse on `[P₀, B(½)]` and `[B(½), P_d]` with `r − 1`.

The pinned points left after `r` rounds are the control polygons of the `2^r` pieces; the path
through them is the output, converging to the geodesic Bézier curve (Morera et al. 2008, cited:
`fo-geodesic-bezier`). On a flat domain the pinned points are exactly the planar de Casteljau
subdivision points (`fo-de-casteljau-midpoint`, §7). Splits go through the representation
(§1.3); a split edge's segments are replaced by segments on both halves in the same layer order
(the reference supports one segment only, gc: fg.cpp:1744-1747). Bézier works on one open path;
the integer split is a floating-point decision only in which region of the face the new vertex
lands (IC §4.4, I6), its connectivity stays exact.

### 6.2 Loops, networks, curve shortening

Loops from a closed edge path or closed waypoints; networks from edge sets (seams, cuts,
segmentation boundaries, Fig. 19, 21, 22); early stopping (Fig. 20, 21) by `stop_length_ratio` /
`stop_min_angle`.

### 6.3 Intrinsic constrained Delaunay triangulation (Sec. 6.3, Fig. 25)

After shortening, `make_delaunay(&DelaunayOptions)`: Delaunay flipping (IT N6, IC §4.3) with every
edge that carries a segment constrained (gc: fg.cpp:1035-1045). Same vertex set, every path an
edge path. Report as `DelaunayReport`, plus the number of constrained edges. Needs §11 R3.

### 6.4 Later: refinement conforming to the paths (Sec. 6.4, Fig. 24, 25 right, 29)

Delaunay refinement treating constrained edges like boundary edges (Chew 1993, Shewchuk 1997; the
reference: marked edges are fixed and stop the circumcenter trace, gc: it.ipp:8-12,
it.cpp:244-246, 359), with path segments following the splits (gc: fg.cpp:1047-1066, 1737-1810).
Needs §11 R5 in both crates; the reference's integer split of such edges throws (IC reference map
B5). WP F9, after everything else.

## 7. Verification

**Layers.**

1. **Certificates** of every run (in tests and the evaluation): (a) on `T`: §2.1 validity, every
   interior joint not pinned has, on each side, `α ≥ π − certify_tolerance` (default `1e-4` rad,
   the paper's Sec. 5.3 criterion; the reference only declares `isStraight(1e-4)`, gc: fg.h:188) or
   a wedge blocked by an end (F27); reports the minimum `α^min` and the number of blocked joints; (b) on `T⁰`: §5; (c) for
   integer coordinates, IC `certify()` once at the end of a run (F25).
2. **Differential, signposts vs integers** (§1.3): same flip sequence, same final segments,
   bitwise-identical lengths; polylines within `1e-9 ·` total length where the signpost trace
   succeeds.
3. **Lean oracle** (`lean/oracle`, no Mathlib, Float): Alg. 1 and Alg. 2 written from the paper
   and §4 of this file (queue order included) on Δ-complexes given by gluing data and lengths;
   output: flip sequence, final segments, lengths. Rust replays committed fixtures: flip sequences
   and segments exactly, lengths to `1e-12` relative. Fixtures include a wedge with a twice-incident
   edge and a degree-1 vertex (App. A cases; the report counts `twice_incident_wedges` so the
   tests can assert they were exercised).
4. **Lean proofs** (§8).
5. **Against the reference** (Q2): a headless driver of gc's `FlipEdgeNetwork` (signposts) on the
   same trials: final lengths to `1e-9` relative where both reach the same geodesic, success rates,
   flip counts as statistics (the angle formula and ties differ, F2, F4).

**Tests** (closed surfaces and surfaces with boundary; `ddg_core` generators, `data/models/`,
the bad meshes of IT §5):

1. *Flat convex domain* (planar triangulated polygon, flat interior vertices): the geodesic between
   two vertices is the straight segment; length = Euclidean distance to `1e-12` relative; the
   extracted polyline is collinear.
2. *Flat non-convex domain* (L-shape, a polygon with holes): the result is the shortest path in its
   homotopy class; for an L-shape the visibility-graph distance (the generator of
   `ddg_intrinsic_apps::exact`'s tests), to `1e-12`.
3. *Closed forms*: the lateral surface of a pyramid (a polyhedral cone; unfold, `Θ < 2π` apex never
   on the path), a cube (corner to opposite corner: `√5 a`), a loop around a prism tube (the
   cross-section's perimeter), a path through a saddle fan with `Θ = 3π` and both sides `≥ π`
   (unchanged, zero flips), a convex boundary corner (pushed off) and a reflex one (kept).
4. *Exact distance*: `ddg_intrinsic_apps::exact::exact_distance` (MMP-type, IT-apps N37) between
   the path's ends: `|γ| ≥ d_exact (1 − 1e-9)` always; the fraction with `|γ| ≤ d_exact (1 + 1e-6)`
   is reported (paper Sec. 6.5: > 95 % minimal). Not gc's `exact_polyhedral_geodesics`, which at
   this commit is a radius-limited greedy unfolding with its MMP code commented out
   (gc: `include/…/exact_polyhedral_geodesics.h:14-21`; F23).
5. *Monotone length*: each Alg. 1 replacement strictly shortens (paths and networks); for loops the
   App. B moves are the only increases and are counted.
6. *Isotopy*: §2.1 non-crossing after every iteration (`check_invariants`); multiply covered input
   (Fig. 8 top as an open path spiralling twice around a tube, and as two parallel loops given as
   two paths; a path pulled around its own end, Fig. 8 bottom). A single closed path winding twice
   around a tube crosses itself under every layering (a curve in the class `a²` of an annulus is not
   simple), so §2.1 refuses it; the reference does not check.
7. *Bézier*: on a plane, the pinned points after `r` rounds equal the planar de Casteljau
   subdivision points to `1e-12`; on curved inputs, symmetric under reversing the control points
   (same point set) to `1e-12`.
8. *Representation identity* (layer 2) on every input of items 1–7.
9. *Robustness* (IntCoords Fig. 17 left): the evaluation's trials where a signpost trace fails are
   committed as fixtures (mesh id, seed); on them the integer extraction must succeed and pass both
   certificates.
10. *Determinism*: results bitwise identical at 1 and 64 threads; a refused input changes nothing.

## 8. Lean: model and theorems

**Model.** (B) The Euclidean plane (Mathlib): a wedge laid out as a fan of triangles around `b`
with apex angles summing to `α < π` (an embedded fan, so planar statements suffice). (C) A
combinatorial wedge in a Δ-complex: the corners `n_i b n_{i+1}` of `W` (distinct) with the
identification of outer corners and twice-incident edges of App. A, a flip acting on it; and a
vertex model for non-crossing: the cyclic slot order at `b` and joints as chords.

| id | statement | tier |
|---|---|---|
| `fo-wedge-angle`, `fo-boundary-wedge` | definitions (Eq. 1) and `α_L + α_R = Θ_b` at an interior vertex | A |
| `fo-locally-shortest`, `fo-geodesic-saddle` | Eq. 2 ⇒ Eq. 3: `Θ_b ≥ 2π`, `Ω_b ≤ 0` | A |
| `fo-outer-arc-angle` | `β_i` is the interior angle of the unfolded fan region at `n_i` | A |
| `fo-wedge-flip-convex` | `β_j < π ∧ α < π` ⇒ the diamond of `b n_j` is strictly convex (uses `it-flip-convexity`) | B |
| `fo-outer-arc-shorter` | Thm 4.1: in the fan, if all `β_i ≥ π` then `|a n₁ … c| < |ab| + |bc|`. Proof by cutting off: extend `a n₁` to `p ∈ bc`, induct on `k`, triangle inequality (no Crofton formula) | B |
| `fo-straight-crossing` | two pieces unfolded across an edge are collinear iff the angle is `π` (uses `it-unfold-diamond`) | A |
| `fo-loop-perturbation` | App. B: Eq. 4, 5 ⇒ Eq. 6 | A |
| `fo-de-casteljau-midpoint`, `fo-bezier-midpoint-split` | planar midpoint subdivision gives `B(½)` and the halves' control polygons; the split parameter | A |
| `fo-wedge-flip-valid` | Lemmas A.1–A.4, Thm A.5, A.6 (topological part) in the wedge model | C |
| `fo-flipout-terminates` | Thm A.7, A.8: `(k, ρ)` decreases lexicographically; at most `k₀²` flips | C |
| `fo-non-crossing`, `fo-isotopy` | replacing a flexible joint by its outer arc, innermost on the outer-arc edges, keeps every vertex's chord diagram non-interleaving | C |
| `fo-path-terminates` | Thm 4.2 given the finiteness of geodesics of length `≤ L` as a hypothesis | A |
| `fo-geodesics-finite` | Indermitte, Liebling, Troyanov, Clémençon 2001, Prop. 1 | E, CITED |
| `fo-flexible-progress` | Sec. 4.2: a non-geodesic network has a flexible non-straight joint or is blocked by ends | D attempt, else CITED (paper's argument) |
| `fo-loop-termination` | App. B conjecture | CITED as open |
| `fo-geodesic-bezier` | Morera, Carvalho, Velho 2008 | CITED |
| `fo-polyhedral-geodesic` | locally shortest ⇔ Eq. 2 at vertices, straight in unfolded strips (Mitchell, Mount, Papadimitriou 1987, Lemma 3.5) | CITED |

Tiers: A formulas, B one planar wedge, C the combinatorial wedge/vertex model, D whole network,
E cited. WP F5 does A, B; WP F6 does C, D and can be cut without blocking anything else.

## 9. Measurements

### 9.1 What the paper reports (the floor to match or beat)

* Machine: one core of an i7-4790K, 16 GB (Sec. 5.1). Timings include extracting the path.
* Datasets (Sec. 5.1): manifold, non-degenerate models of Princeton Shape Benchmark (~600), MPZ
  (Myles et al. 2014, ~100), Thingi10k (~6200). Task: per mesh 10 trials, two random connected
  distinct vertices, Dijkstra path, Alg. 2, extraction; trials already geodesic omitted from the
  timing plots.
* Speed (Sec. 5.2): typically a few ms; 2M faces 9 ms; worst 0.5 s (1500 segments to one).
  Figures: Fig. 4 < 1 ms; Fig. 14: 12, < 1, 5 ms; Fig. 15 (loops): 5, < 1, 3 ms; Fig. 17: 17 ms,
  3938 iterations, 14k flips (iterative relaxation 3.3 s, 4.5 M iterations); Fig. 19: 26, 50 ms;
  Fig. 20: 12 ms; Fig. 21: 4, 14 ms; Fig. 22: < 1 ms; Fig. 29: 27 ms, 361 inserted vertices
  (Lagrangian + extrinsic refinement 742, + intrinsic 898).
* Fig. 18 (MPZ): runtime normalized by Dijkstra's ≈ 1.1 (FlipOut adds ~10 %), ~100× faster than
  Martínez et al. 2005, comparable to Xin and Wang 2007.
* Fig. 16 (MPZ): flips `∝ m^1.5`, `m` the initial number of segments; `O(V^0.75)` per shortest path.
* Sec. 5.3: geodesic to `1e-4` rad on 99.4 % of the Sec. 5.1 trials (Xin–Wang 96.2 %, Martínez
  84.4 %); 100 % on MPZ; failures attributed to extremely skinny triangles.
* Sec. 6.5, Fig. 27 (gear): single-source tree 348 ms, 99.6 % minimal (Xin–Wang per path 6240 ms,
  84.9 %); > 95 % minimal to `1e-6`, max error 1.04 % (out of scope here, F21).

### 9.2 IntCoords Sec. 5.2

FlipOut on integer coordinates, the Sec. 5.1 task repeated on Thingi10k: ~62 000 trials (≈ 6200
meshes × 10, so FlipOut's subset), every one a connected polyline. No timings. Fig. 17 left: a
model where the signpost extraction breaks and the integer one does not.

### 9.3 Our plan

* **Thingi10k** (`data/thingi10k/Thingi10K.tar.gz`, downloaded for IC-SPEC Q1): every mesh
  `TriMesh` accepts with valid triangles (counts reported against ~6200 and IC's 7696); 10 trials per
  mesh with seeds from `(thingi id, trial)`; the same trials on both representations, unmollified
  (the paper) and mollified `1e-5 h̄` (IC §7, F22). Per trial: success on `T` (§7 1a), extraction
  success (signposts: every trace passes N8; integers: connected and certified), the `T⁰`
  certificate, length, `|γ| / d_exact` on meshes up to 20 000 vertices, times of Dijkstra, FlipOut
  and extraction, flips, iterations, `m`, refused flips, stalled joints, sweeps. One job per mesh,
  64 in parallel, 10-minute cap, one log per job (`runs_tmp/fo-thingi/`). Targets: integers 100 %
  connected (IntCoords), `≥ 99.4 %` straight (paper), signposts' extraction failures counted.
* **Timings** single thread, alone, on a 50-mesh subset and `data/models`: median and worst per
  trial; FlipOut time / Dijkstra time (floor ≈ 0.1 extra); fit of flips vs `m` (floor exponent
  1.5). Benches (`benches/`): one path on spot, bob, a 2M-face mesh; a loop; a Bézier with 4 rounds;
  the seam network of a `ddg-uv` cut.
* **MPZ, PSB** if Q1 is answered yes (Fig. 16, 18 and the 100 % on MPZ).
* **Figures** reproduced qualitatively on `data/models` (paths, loops, Bézier, networks, iCDT) as
  Blender/Python smoke tests, with their times.

## 10. Decisions

* **F1** One generic algorithm (`GeodesicNetwork<T: Triangulation>`), the trait in the shared layer
  (§11 R1); no representation-specific code outside extraction.
* **F2** Wedge angles by summing corner angles (Eq. 1), not signpost differences (gc:
  fg.cpp:489-522): identical on both representations; ddg-intrinsic's signposts are scaled (IT N2).
  Cost `O(deg b)` per test; measured, cached per vertex only if the timings ask for it.
* **F3** Queue tolerance `1e-5` (gc: fg.h:161, fg.cpp:534-543); certificate `1e-4` (paper Sec. 5.3).
* **F4** Total queue order `(α, side Left < Right, path, segment id)`.
* **F5** Stale and blocked entries as the reference (gc: fg.cpp:937-952, 1009-1013), plus the final
  sweep of §4.2 (the reference ends at an empty queue).
* **F6** Refused flips skipped, the length guard of gc: fg.cpp:836-838, `stalled` counted.
* **F7** Edges carrying segments are never flipped; a violation is `Internal`, not a panic.
* **F8** Segments by `(edge, forward)`; the shared connectivity must keep an unflipped edge's
  direction (§11 R1).
* **F9** `closed` is explicit (gc infers it from first = last vertex, fg.cpp:234, 265, 1303), so an
  open path from `x` to `x` (a loop through a pinned point, Sec. 2.2 inset) is expressible.
* **F10** Pinning is explicit and always honoured. The reference marks vertices but straightens
  through them by default (`straightenAroundMarkedVertices = true`, gc: fg.h:160), honouring marks
  only during Bézier subdivision (fg.cpp:1070-1081).
* **F11** Default layering as the reference (gc: fg.cpp:29-32), then the non-crossing check;
  crossing inputs refused.
* **F12** `deg(W) = 0` contracts the joint; a loop that contracts is removed and reported (the
  reference has no explicit case).
* **F13** Loops run under `max_iterations`, default `100 (|E| + m₀)`; paths have the same default
  limit (the reference has none), reported as `IterationLimit`.
* **F14** Stop rules `max_iterations`, `stop_length_ratio`, `stop_min_angle`; the length is a running
  sum (the reference recomputes it every iteration, gc: fg.cpp:968-969).
* **F15** Dijkstra on the current `T`, ties by `(distance, halfedge id)`; errors for same vertex and
  disconnected ends. A boundary edge walked against its only halfedge gets the key
  `n_halfedges + edge` (the reference has exterior halfedges; we don't). `path_from_vertices(tri,
  vertices, closed)`: where two consecutive vertices share several edges the shortest is taken; a
  closed path through one vertex is a loop along a self-edge.
* **F16** No rewind (gc: fg.cpp:1250-1313): trials clone the triangulation (input behind an `Arc`,
  IT N9); clone time reported, excluded from FlipOut time.
* **F17** Integer extraction per edge (§11 R2), not via the whole common subdivision (gc:
  ic.cpp:69-88; IC reference map §3.5).
* **F18** A failed signpost trace is an error naming the edge, as IT N8; never repaired here.
* **F19** Bézier: the reference's recursion; midpoint snapped to either end within `1e-4 ℓ`; splits
  of edges with several segments supported.
* **F20** iCDT in scope; conforming refinement later (§6.4).
* **F21** Out of scope: single-source tree and log map (not in the reference; the paper's
  description is a sketch), Lagrangian baselines, PDE/cross-field demos.
* **F22** Mollification off in the library; both settings in the evaluation.
* **F23** Ground truth for minimality: `ddg_intrinsic_apps::exact`.
* **F24** Single-threaded algorithm; extraction parallel per edge, assembled in order.
* **F25** Integer `certify` once per run, not per flip.
* **F26** *Out-and-back joint.* The reference measures `Θ` on both sides when `ab` and `bc` share an
  edge (`measureSideAngles` with equal signposts), so it never sees `α = 0`; here the layering puts
  angle `0` on the side between the two segments (§3), which Alg. 1 then shortens.
* **F27** *Certificate per side.* "The small wedge is blocked by an end" (the earlier wording)
  accepts a non-geodesic at `Θ_b < 2π`, where both sides can be `< π` and the large side may still
  shorten; every side with `α < π − tol` must be blocked by an end.

## 11. Requests to the other contracts

To IC-SPEC (and, for R1, R3, R4, R5, IT-SPEC's crate through IC WP0):

* **R1** (I1/WP0) The shared layer exposes a representation-neutral trait used by FlipOut and by
  IC §4.3: the connectivity, lengths, corner angles, `can_flip`/`flip_edge` with IT §4 + N5
  refusals, `split_edge(e, t)` returning both halves oriented like `e`, the constrained-edge mask
  (R3), `trace_edge_on_input(e) -> ddg_intrinsic::EdgeTrace`, `vertex_location(v)`. Guarantee: a
  flip keeps the direction of every edge it does not flip (ddg-intrinsic does, `flip.rs:225-231`).
  **Done (IC WP0, IT N50):** the trait is `ddg_intrinsic::Triangulation` over the connectivity
  `ddg_intrinsic::Complex`; `trace_edge_on_input(e, &EdgeTraceOptions)`; `split_edge` returns
  `Insertion` with `halves: Some([e, e2])`; flips of constrained edges are refused with
  `FlipRefusal::Constrained`. The direction guarantee (F8) is documented on `flip_edge` and tested
  in `crates/ddg-intrinsic/tests/complex.rs`.
* **R2** (§3.1, new operation) `ExtractIntrinsicEdge(e)`: an edge of `T¹` as a path over `T⁰`,
  returning `ddg_intrinsic::EdgeTrace`; each crossing's input edge identified exactly (Alg. 1 + 3,
  memoized along each walk), consecutive input edges checked to share an input face (else
  `InvalidCoordinates` naming the edge), positions from the straight segment in the unfolded strip
  of input faces (input lengths), clamps per I5; `n = −1` handled as an input-edge piece; ends in
  `V*` from `q⁰`. Test: positions agree with those from IC Alg. 2's strip (the other direction) to
  `1e-9`; the common subdivision's edges agree with it.
* **R3** (§4.3) Delaunay flipping with a set of constrained edges never flipped (iCDT), identical
  on both representations. **Done:** the provided `Triangulation::flip_to_delaunay` (one
  implementation for both); `DelaunayReport::n_constrained`, and the other statistics leave
  constrained edges out.
* **R4** (§4.5) Edge split: accepts constrained edges (the flip inside "face split then flip" is part
  of the split, fixing reference bug B5); returns both halves; assigns lengths `tℓ`, `(1 − t)ℓ` and
  `it-point-distances`, exactly as ddg-intrinsic's `split_edge` (IT N10), so both representations
  stay bitwise identical through Bézier subdivision. **Done:** `ddg_intrinsic::plan_edge_split`
  decides refusal and lengths from connectivity and lengths alone; both crates call it; splitting a
  constrained edge gives two constrained halves.
* **R5** (§5, later) Refinement treats constrained interior edges as boundary edges (IT N1's
  encroachment and splitting rules), in both crates. **Open:** until then `refine`, `refine_with`
  and `optimal_delaunay` refuse a triangulation with constrained edges
  (`IntrinsicError::ConstrainedEdges`) and change nothing.
* **R6** (§2.4, I4) `certify` is not run per flip by default; FlipOut makes thousands of flips per
  path and certifies once at the end.
* **R7** (§10.1 last bullet, §13 Q4) Answered: FlipOut is in scope; the Sec. 5.2 benchmark lives in
  this file (§9.2, §9.3). (§13 Q1: Thingi10k is already in `data/thingi10k/`.)
* **R8** (§13 Q5) Answered by the user: a visible Blender operator for the refinement / common
  subdivision is wanted; its drawing of `T¹` edges on the input mesh is reused here (§13).

## 12. Work packages

| WP | content | files | depends on | model |
|---|---|---|---|---|
| F0 | **done** with IC WP0 (2814e8a–aca460e): R1, R3, R4 | `crates/ddg-intrinsic/src/…` | — | Opus |
| F1 | crate skeleton: network, segments, layering, pinning, validity and non-crossing check, wedge angles, errors, reports | **done** (2144fed, d3c11e3): `src/{lib,error,network,wedge,certify,init,report}.rs` (`init`: Dijkstra and vertex-sequence paths, F3 adds waypoints and edge sets; `report`: options and reports for F2) | F0 (R1 for signposts only) | Opus |
| F2 | Alg. 1, Alg. 2, loops, contraction, stop rules, final sweep; tests 1–6 on signposts | `src/{flipout,shorten}.rs`, `tests/{flat,closed_form,loops,isotopy}.rs` | F1 | Opus |
| F3 | initial paths (Dijkstra, waypoints, edge sets); extraction both ways; `T⁰` certificate; 3D polylines | `src/{init,extract}.rs`, `tests/extract.rs` | F1; IC WP2 + R2 for the integer half | Sonnet |
| F4 | Bézier, split bookkeeping, iCDT; tests 7, 8 | `src/{bezier,cdt}.rs`, `tests/{bezier,identity}.rs` | F2, IC WP3/WP4 + R4 | Opus |
| F5 | Lean tiers A, B; oracle Alg. 1/2 and fixtures; `CITED.md` section | `lean/proofs/DdgProofs/FlipGeodesics/…`, `lean/oracle/DdgOracle/FlipGeodesics…` | this file | Opus |
| F6 | Lean tiers C, D | `lean/proofs/DdgProofs/FlipGeodesics/Wedge/…` | F5 | Opus |
| F7 | evaluation: Thingi10k both representations, exact distance, timings, benches, reference driver (Q2), MPZ/PSB (Q1) | `examples/fo_thingi.rs`, `benches/`, `scripts/fo-eval.sh` | F3, F4 | Sonnet |
| F8 | Python and Blender (§13) | `crates/ddg-py/src/flip_geodesics.rs`, `blender/ddg_uv_unwrap/geodesic_operator.py`, `blender/tests/test_geodesic.py` | F3, F4 | Sonnet |
| F9 | later: conforming refinement (§6.4, R5) | both intrinsic crates + `src/refine.rs` | F4, IC WP6 | Opus |

F1–F3 run on signposts before the integer crate exists; the integer half of F3/F4 waits for IC
WP2–WP4. At most three agents at a time; timings alone.

## 13. Python and Blender

**Python** (`ddg_tools`): `geodesic_path(positions, triangles, vertices, *, closed=False,
pin_waypoints=True, representation=…, mollify=None, tolerance=1e-5, max_iterations=None,
stop_length_ratio=None)`, `geodesic_network(positions, triangles, edges, *, pinned=None, …)`,
`geodesic_bezier(positions, triangles, control_vertices, *, rounds=4, …)`. Each returns
`polylines` (list of `(n, 3)` arrays), `surface_points` (per point: input vertex / edge + `t` /
face + barycentric), `lengths`, and the report and certificates as a dict. `representation`
defaults as IC-SPEC Q3 decides for all consumers.

**Blender** (proposal, Q3). Operator "Geodesic Path" (`mesh.ddg_geodesic_path`), Edit Mode, in the
add-on's sidebar panel:

* Input: selected vertices in selection order (`bmesh.select_history`; with more than two vertices
  and no complete history, an error asking to shift-click them in order). Mode: *Path* (geodesic
  through the points, points pinned), *Smooth* (points only seed the path: one geodesic between the
  first and last), *Loop* (closed through the points, pinned or not), *Bézier* (points are control
  points; rounds 1–8, default 4). Options: representation, mollify, stop at length ratio.
* Output: a new Curve object `<mesh>_geodesic` (poly splines from `polylines_3d`, one per path),
  parented to the mesh with its transform, bevel depth 0.2 % of the bounding-box diagonal, custom
  properties: length, iterations, flips, representation, certificate. The mesh is not changed.
  Re-running replaces the object of the same name.
* "Straighten Seams" (`mesh.ddg_straighten_seams`): the edges marked as seams as a network
  (`network_from_edges`, selected vertices extra-pinned), stop ratio option (default none; 0.9 as a
  preset for segmentations), output as above.
* Option "Show triangulation": also draws the edges of the final `T` on the input mesh as a second
  curve object, using the drawing of the IC operator (R8).
* Tests: `blender/tests/test_geodesic.py` in the bpy module (as the other operators).

## 14. The user's decisions (2026-10-31)

* **Q1** Download MPZ (Myles et al. 2014) and the Princeton Shape Benchmark into `data/` (with
  Thingi10k, already in `data/thingi10k/`) and repeat Fig. 16, Fig. 18 and the 99.4 % / 100 % figures
  over the paper's three datasets.
* **Q2** Build the headless driver of gc's `FlipEdgeNetwork` (`ab116f1`, signposts, geometry-central
  + Eigen only) outside the repository, for §7 layer 5.
* **Q3** The Blender proposal of §13 as written: "Geodesic Path" (Path / Smooth / Loop / Bézier,
  Curve object output) and "Straighten Seams", with "Show triangulation".
