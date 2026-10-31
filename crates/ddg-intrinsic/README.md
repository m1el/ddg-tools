# ddg-intrinsic

Intrinsic triangulations of a triangle mesh in the signpost representation of

> Nicholas Sharp, Yousuf Soliman, Keenan Crane. *Navigating Intrinsic Triangulations.*
> ACM Transactions on Graphics 38(4), 2019. <https://doi.org/10.1145/3306346.3322979>

An intrinsic triangulation is another triangulation of the same polyhedral surface: the same
points, the same distances, angles and areas, other edges. Its edges are straight paths over
the input mesh, its triangles are flat, and it need not be a simplicial complex. Flipping
edges until it is Delaunay removes the negative cotan weights of the input without moving a
single point of the surface.

The contract of the crate is [`SPEC.md`](SPEC.md): scope, conventions, what must hold, and
every decision the paper leaves open (entries `N1` to `N25`).

## What is implemented

- The data structure: a halfedge mesh that allows what flips produce (two edges between the
  same vertices, an edge from a vertex to itself, a triangle that uses a vertex twice,
  vertices of degree one), a length per edge, a signpost per halfedge (normalized to `2π` at
  interior vertices and to `π` at boundary vertices), the angle sum of every vertex, and the
  place of every vertex on the input mesh (paper Sec. 3.1).
- Building it from a `ddg_core::TriMesh` with a `Geometry` (Alg. 4).
- The edge flip (Alg. 5) with every refusal: boundary edge, end point of degree one, not
  strictly convex, too thin (SPEC N5). The caller can choose the shape test per flip: N5's
  (the default, used by every operation here) or the reference's relative one, which FlipOut
  uses (`FlipShapeTest`, SPEC N52). A flip returns the two corner angles it computes at the new
  edge's ends (`Flip::new_corner_angles`).
- Delaunay flipping (Sec. 4.1) with a report, a tolerance, a flip limit and the count of ties.
- Tracing queries (Alg. 2, Sec. 3.2.2) over the input mesh and over the intrinsic
  triangulation, from a vertex, from a point of an edge and from a point of a face, with the
  crossings on request.
- The path of every intrinsic edge over the input mesh, checked: a trace that does not arrive
  at the edge's far vertex is an error that names the edge (SPEC N3, N8). Nothing is perturbed
  or repaired.
- Vertex insertion in a face and edge split, boundary edges included (Alg. 7, 8); the new
  vertex is placed on the input mesh by a tracing query (SPEC N10).
- Point queries both ways (Alg. 11), and a straight segment inside a face of the
  triangulation as the sequence of its pieces on input faces (SPEC N11, N13).
- The overlay (Sec. 3.4.2): the common subdivision as a mesh of convex polygons, every vertex
  located on both triangulations (SPEC N12).
- Delaunay refinement (Sec. 4.2) to an angle bound, on closed surfaces and with boundary
  (Ruppert's rule, SPEC N1), with vertices of small angle sum (SPEC N4). Everything left below
  the bound is listed with its reason. The driver is written once over `Triangulation`
  (`refine::refine`, SPEC N51), so the integer coordinates run the same code; it also has the
  rules of Gillespie, Sharp and Crane 2021 (`RefineRules::Chew`: midpoint splits of the boundary
  edge a circumcentre path leaves through, removal of the inserted vertices in its Dijkstra
  ball, the reference's narrow-cone rule), opt-in here, the default of `ddg-integer-coords`.
- The simplicial export for code that takes a `TriMesh` with a `Geometry`, the way back for a
  cut found there, and a fingerprint that identifies the state (SPEC §6, N15).
- The cotan Laplacian (positive semidefinite) and the lumped vertex areas of the triangulation,
  angle defects, and quality numbers.

- Vertex repositioning (Sec. 3.3.3, Alg. 10): an inserted vertex moved by a vector in its own
  frame inside its flat star, with an exact test of the star, the place on the input mesh by a
  trace, vertices on the boundary along the boundary (SPEC N19, N20, N24).
- The intrinsic optimal Delaunay triangulation (Sec. 4.3, after Chen and Holst 2011): splits to
  a target length, the area-weighted circumcentre update with a line search, Delaunay flips,
  with a report and the angle and area distribution (SPEC N21 to N25).

- The connectivity on its own (`Complex`: flips, face and edge splits, removal of a vertex of
  degree three, a constrained-edge mask, each operation with its index bookkeeping), the trait
  `Triangulation` that this representation and the integer coordinates of Gillespie, Sharp and
  Crane (2021) share, Delaunay flipping written once over it, and tracing over a triangulation
  that stores no signposts (`AngleFrames`), SPEC N50.

Not implemented: the paper's general route for a move beyond the star of a vertex (SPEC N20;
moves do not remove vertices; removal by flips exists since N51); refinement with constrained edges (N50); a
refinement that keeps
the angle bound through the export at spikes (SPEC N15: on thin tubes and handles the export
can restore it, at vertices with a small angle sum it cannot). Benches: `benches/retriangulation.rs`
(criterion: iDT, iDR, iODT and a sweep of moves); the tables' timings are those of the example
programs. The sixteen `it-` ids below have machine-checked theorems in
`lean/proofs/DdgProofs/Intrinsic/` (about the oracle's definitions in
`lean/oracle/DdgOracle/Formulas/Intrinsic.lean`); what is cited instead of proved (existence,
uniqueness and termination for the Delaunay triangulation, Rippa's theorem, the planar
refinement guarantees, straightest geodesics) is in `lean/proofs/CITED.md`. A whole tracing
query, the overlay, refinement and the export have no theorem.

## Use

```rust
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{DelaunayOptions, DelaunayStatus, EdgeTraceOptions, IntrinsicTriangulation};

let mesh = TriMesh::new(positions.len(), faces)?;
let geometry = Geometry::from_positions(&mesh, &positions)?;

let mut triangulation = IntrinsicTriangulation::new(&mesh, &geometry)?;
let report = triangulation.flip_to_delaunay(&DelaunayOptions::default())?;
assert_eq!(report.status, DelaunayStatus::Delaunay); // or FlipLimit, Blocked: see the report

let laplacian = triangulation.cotan_laplacian();      // same conventions as ddg_core
let mass = triangulation.lumped_mass_normalized();

// Where the edges run on the input mesh. Fails if a trace misses its vertex.
let traces = triangulation.trace_edges(&EdgeTraceOptions::default())?;
for e in 0..triangulation.n_edges() {
    for crossing in traces.crossings(e) {
        // crossing.edge, crossing.t: a point of an input edge; crossing.distance: along e
    }
}
```

Refinement, and the cut flow on the result (`tests/export.rs` runs exactly this):

```rust
use ddg_intrinsic::{EdgeTraceOptions, ExportOptions, RefineOptions, RefineStatus, Segment};
use ddg_variational_cuts::{Cut, Flow};

let mut triangulation = IntrinsicTriangulation::new(&mesh, &geometry)?;
let options = RefineOptions::default();                 // 25 degrees
let report = triangulation.refine(&options)?;
assert_eq!(report.status, RefineStatus::Done);          // or InsertionLimit
for corner in &report.below_bound { /* face, vertex, angle, reason */ }

// A simplicial mesh of the same surface. `refine` restores the bound after splits where it can.
let export = triangulation.export_simplicial(&ExportOptions {
    refine: Some(options),
    ..ExportOptions::default()
})?;
// export.mesh, export.geometry: what the flow takes. export.locations[v]: where vertex v is
// on the input mesh. export.min_angle, export.max_angle: what the export really has.
let mut flow = Flow::new(&export.mesh, &export.geometry, &phi)?;
// ... flow.step(...) ...; keep (export.fingerprint, flow.level_set()) as the state.
export.check_state(saved_fingerprint, saved_phi.len())?; // refuses a state of another export

// The cut, back on the input mesh.
let cut = Cut::from_level_set(&export.mesh, &flow.level_set(), 1e-3)?;
let nodes: Vec<_> = cut.points.iter().map(|p| export.edge_point(p.edge as usize, p.t)).collect();
let segments: Vec<_> =
    cut.segments.iter().map(|s| Segment { face: s.face, nodes: s.points }).collect();
let paths = triangulation.segments_on_input(&nodes, &segments, &EdgeTraceOptions::default())?;
// paths.nodes[i]: node i as a point of the input mesh (vertex, edge or inside a face).
// paths.paths[k].points / .faces: segment k as straight pieces, each in one input face.
```

`examples/delaunay_report.rs` reads OBJ files and prints what the tables below are made of.

## API

| item | what it is |
|------|------------|
| `IntrinsicTriangulation::new(&TriMesh, &Geometry)` | the input mesh as an intrinsic triangulation, with the same ids |
| `n_vertices`, `n_edges`, `n_faces`, `faces`, `next`, `prev`, `twin`, `from_vertex`, `to_vertex`, `edge`, `edge_halfedge`, `edge_vertices`, `outgoing_halfedges`, `vertex_degree`, `is_boundary_edge`, `is_boundary_vertex` | connectivity; halfedge `h = 3 f + k` as in `ddg_core` |
| `edge_length`, `signpost`, `vertex_angle_sum`, `vertex_location`, `edge_on_input`, `corner_angle`, `face_area`, `total_area` | geometry and the correspondence with the input |
| `flip_edge(e) -> Result<Flip, IntrinsicError>`, `can_flip(e)` | one flip, or the reason it is refused (`FlipRefusal`); `Flip::new_corner_angles` at `new_corners` (SPEC N52) |
| `flip_edge_with(e, FlipShapeTest)`, `can_flip_with(e, _)`, `plan_flip_with(&Complex, &lengths, e, _)` | the same with the shape test chosen by the caller: `FlipShapeTest::N5` (default, `Absolute { limit: 1e-6 }`) or `FlipShapeTest::REFERENCE` (`Relative { ratio: 1e-6 }`, the reference's); also on the trait and the integer coordinates (SPEC N52) |
| `flip_to_delaunay(&DelaunayOptions) -> Result<DelaunayReport, _>`, `is_delaunay`, `edge_angle_excess`, `edge_cotan_weight` | Delaunay flipping; `DelaunayCriterion::{AngleSum, CotanWeight, REFERENCE}`, `FlipLimit` |
| `InputSurface::trace`, `IntrinsicTriangulation::trace` `(start, direction, length, &TraceOptions) -> Result<Trace, _>` | a straight path over the input mesh / over the triangulation from a `SurfacePoint` |
| `trace_edge(e, &EdgeTraceOptions)`, `trace_edges(&EdgeTraceOptions) -> Result<EdgeTraces, _>` | the checked path of an edge / of all edges over the input mesh |
| `cotan_laplacian`, `lumped_mass`, `lumped_mass_normalized`, `corner_angle_sums`, `angle_defects`, `boundary_geodesic_curvatures`, `quality` | operators and statistics |
| `insert_vertex(f, bary)`, `split_edge(e, t) -> Result<Insertion, _>` | a new vertex in a face / on an edge, placed on the input mesh; or the reason it is refused (`InsertRefusal`) |
| `point_on_input(point)`, `point_on_triangulation(point)` | a point of one triangulation as a point of the other |
| `complex() -> &Complex`, `set_constrained(e, bool)`, `is_constrained(e)`, `remove_vertex(v)` | the shared connectivity, constrained edges, removal of an inserted vertex of degree three (SPEC N50) |
| `Triangulation` (trait), `AngleFrames::of(&t).trace(...)` | the interface shared with integer coordinates; tracing with signposts computed from corner angles (SPEC N50) |
| `segments_on_input(&nodes, &segments, &EdgeTraceOptions) -> Result<PathsOnInput, _>`, `segment_on_input(f, from, to, _)` | straight segments inside faces of the triangulation as pieces on input faces, checked like edge traces |
| `overlay(&EdgeTraces) -> Result<Overlay, _>` | the common subdivision: `vertices` (on both triangulations), polygons `face(f)` with `face_input`, `face_triangulation`; `triangles()` |
| `refine(&RefineOptions) -> Result<RefineReport, _>` | Delaunay refinement; `min_angle` (default 25°), `max_circumradius`, `InsertionLimit`; the report lists `below_bound` with `BelowReason` |
| `refine_with_rules(&RefineOptions, RefineRules)` (trait), `refine::refine(&mut t, _, rules, extra)` | the same driver on any representation; `RefineRules::Ruppert` (default, N1, N4) or `RefineRules::Chew` (IC Alg. 10, SPEC N51); the report also counts removals, refused removals, trace failures and the faces left by the narrow-cone rule (`min_angle_valid`) |
| `remove_inserted_vertex(v) -> Result<InsertedVertexRemoval, _>`, `removal_flip_choice` | removal of an inserted interior vertex of any degree by flips (IC I24), shared by both representations; `RemovalRefusal::Stuck` (SPEC N51) |
| `move_vertex(v, [dx, dy]) -> Result<VertexMove, _>`, `odt_target(v, barycenters_at_boundary)`, `star_min_shape(v)` | move an inserted vertex inside its star (or the reason it is refused, `MoveRefusal`); the update of the optimal Delaunay triangulation at a vertex, in its frame |
| `optimal_delaunay(&OptimalDelaunayOptions) -> Result<OptimalDelaunayReport, _>` | iODT; `TargetLength`, `max_iterations`, `motion_tolerance`, `step`, `max_halvings`, `min_shape`, `BoundaryVertices`, `OdtCenters`, `InsertionLimit`; the report has every `OdtIteration` and `ShapeStatistics` |
| `shape_statistics() -> ShapeStatistics` | distribution of the corner angles and the face areas (`Summary`: min, max, mean, deviation, quantiles) |
| `export_simplicial(&ExportOptions) -> Result<SimplicialExport, _>` | splits what is not simplicial; `mesh`, `geometry`, `locations`, `n_splits`, `min_angle`, `max_angle`, `is_delaunay`, `refinement`, `fingerprint`, `edge_point`, `triangulation_edge`, `check_state` |
| `fingerprint()` | a hash of connectivity, lengths and places |
| `irregularities`, `to_tri_mesh`, `check` | what is not simplicial; the triangulation as a `TriMesh` if it is; structure and identities |
| `formulas::*` | the closed-form formulas as plain functions |
| `IntrinsicError`, `FlipRefusal`, `InsertRefusal`, `MoveRefusal`, `EdgeTraceFailure` | every error, with the element concerned |

What the signatures do not say:

- **Ids.** A fresh triangulation has the vertex, edge, face and halfedge ids of the input mesh.
  Vertices, edges and faces keep their ids under flips; halfedges do not (a flip moves two of
  them to the other face). An edge has a direction, that of `edge_halfedge(e)`; a flip turns
  it counter-clockwise.
- **Signposts** are normalized angles: actual angles times `2π/Θ` at an interior vertex and
  `π/Θ` at a boundary vertex, where they run from the boundary halfedge that leaves the vertex
  (0) to the boundary edge that arrives (`π`). The frame of a vertex is that of the input mesh,
  so a signpost can be handed to `InputSurface::trace`.
- **The input** is copied once and shared by clones of the triangulation (`Arc`).
- **A point** (`SurfacePoint`) is a vertex, a point of an edge (`t` from the lower-numbered
  vertex on the input mesh, from the tail of `edge_halfedge` on the triangulation) or a point
  of a face in barycentric coordinates. Directions at a point of an edge are measured from the
  edge, at a point of face `f` from halfedge `3 f`.
- **A trace** ends in a face, on an edge it runs along, at the vertex such an edge leads to
  (`TraceStop::Vertex` if length is left over), or at the boundary (`TraceStop::Boundary`).
  `back_direction` is the direction back along the path in the frame of the end point.
- **Boundary edges are never flipped**, so on a surface with boundary a boundary edge can keep
  a negative cotan weight after Delaunay flipping (SPEC N6).
- **Inserted vertices** are appended and never removed; they are flat. Their place on the
  input mesh is a vertex-free point of an input face or of an input edge; a traced place
  within `1e-12` of the face's longest edge from an input edge is moved onto it and reported
  (`Insertion::snapped`).
- **Moving a vertex** (`move_vertex`) never shortens the step and never leaves the star: a move
  that would turn a triangle over is refused with the fraction that would still be inside
  (`MoveRefusal::NotInStar`). The vector is in the frame of the vertex, which is that of the
  input face or edge it lies in, so the frame changes when the vertex moves into another input
  face; `VertexMove::back_direction` is the direction of the old place in the new frame. Input
  vertices are never moved; a vertex on the boundary only along it (`displacement[1] == 0`).
- **The optimal Delaunay triangulation** is deterministic, Delaunay after every iteration, and
  moves only inserted vertices. With `TargetLength::MeanInputEdge(1.0)` (the default) it adds
  1.5 to 3 times the input's vertices on most meshes; the target is the knob for the vertex
  count. It does not improve angles at input vertices.
- **Refinement is deterministic** and so is the export: the same input and options give the
  same triangulation bit for bit, hence the same fingerprint.
- **`export_simplicial` changes the triangulation** (splits, and flips unless turned off); the
  way back must be asked of the triangulation in that state. **Splitting does not keep an
  angle bound**: read `min_angle`, `max_angle` and `refinement` of the export, not the report
  of an earlier refinement.
- **Three numbers decide what is accepted** (all in SPEC §8): an edge is flipped by Delaunay
  flipping when its opposite angles sum to more than `π + 1e-10`; a flip is refused when it
  would create a triangle with `2 A < 1e-6 · ℓ_new · ℓ_max`; an edge trace must arrive within
  `1e-6` of the edge's length of its vertex. The first and the third are options; the second
  can be chosen per flip with `flip_edge_with` (SPEC N52), and no operation of this crate does.

## What was measured

All numbers are from this repository's tests and example on the machine the crate was written
on (one thread). "Bad meshes" are `camelbad`, `camelbadnonuniform`, `camel_noise1`,
`camel_noise2` and `camelshotnoise` of the optimal-cones data set.

### Delaunay flipping: quality and time

`examples/delaunay_report.rs`, release build, one thread, fastest of 7 runs. Angles in degrees;
"neg." is the share of interior edges with a negative cotan weight (below `−1e-12`); "ref." is
the time of the authors' program for the same operation (`data/intrinsic-reference`,
`stats.txt`, measured by another process on the same machine).

| mesh | faces | flips | min angle before → after | largest angle before → after | faces > 120° before → after | neg. before → after | init ms | flip ms | ref. ms | trace ms |
|---|---|---|---|---|---|---|---|---|---|---|
| spot | 5 856 | 281 | 10.21 → 11.43 | 131.72 → 132.56 | 0.48 % → 0.10 % | 3.06 % → 0 | 0.28 | 1.56 | 1.68 | 0.20 |
| bob | 10 688 | 657 | 11.06 → 11.06 | 119.83 → 103.58 | 0 → 0 | 4.10 % → 0 | 0.50 | 2.86 | 2.75 | 0.41 |
| blub | 14 208 | 1 865 | 4.49 → 6.92 | 147.12 → 158.97 | 4.11 % → 0.58 % | 6.98 % → 0 | 0.73 | 4.23 | 4.56 | 0.94 |
| camelbad | 22 725 | 8 110 | 0.0012 → 1.34 | 180.00 → 164.84 | 15.00 % → 2.23 % | 15.78 % → 0 | 1.40 | 8.83 | 9.34 | 3.44 |
| camelbadnonuniform | 18 125 | 6 966 | 0.0022 → 1.34 | 180.00 → 170.00 | 16.72 % → 2.84 % | 16.72 % → 0 | 0.94 | 7.24 | 7.25 | 2.87 |
| camel_noise1 | 22 704 | 7 584 | 0.112 → 0.807 | 179.22 → 160.49 | 16.34 % → 3.33 % | 16.19 % → 0 | 1.18 | 8.61 | 8.73 | 3.11 |
| camel_noise2 | 22 704 | 7 406 | 0.322 → 1.18 | 179.07 → 161.36 | 20.79 % → 3.86 % | 16.57 % → 0 | 1.21 | 8.61 | | 3.01 |
| camelshotnoise | 22 968 | 6 681 | 0.0030 → 0.026 | 177.30 → 179.71 | 7.78 % → 2.66 % | 14.46 % → 0 | 1.19 | 8.20 | 8.10 | 2.75 |

Every run ends with status `Delaunay`, no refused flip, and the largest angle excess below zero
(`1.7e-15` on blub, which has 43 ties). The flip time is that of testing every edge at least once
(about 5 to 6 million edge tests per second); a flip alone, measured by flipping every edge that
can be flipped there and back, runs at 4.6 to 5.8 million flips per second; tracing at 1.5 to 4
million crossings per second. After flipping, camel_noise1 has 1 self-edge and 5 double edges, camel_noise2 4 and 75,
camelshotnoise 22 and 24; the others are simplicial complexes.

### Edge traces

After Delaunay flipping, with the default tolerance of `1e-6` of the edge length, **no trace
fails its check on any mesh**: 0 of the 281 to 7331 flipped edges per mesh (`tests/models.rs`),
0 on the 17 generated shapes of `tests/common`, 0 on the 15 inputs of the reference comparison.
The largest mismatch between the unfolded far vertex and the end of the path is `2.1e-11` of
the edge length on `camelshotnoise` (24 of its 6095 traced edges are above `1e-12`, none above
`1e-10`) and `7e-14` or less on every other mesh; a tolerance of `1e-9` gives the same result.
The length of each traced path, measured in space from its crossing points, equals the length
of the edge to `2e-12` relative. No value was moved into range (`clamp` is 0 throughout).

Stress beyond what Delaunay flipping does (not a kept test; the first row of each mesh is
`after_random_flips_every_edge_still_traces`): edges flipped at random, `k · n_edges` attempts,
which makes triangles with angles down to `1e-4` degrees and edges that cross up to 46 input
edges.

| attempts | flips made (camelbad) | failed at `1e-6` | largest mismatch | failed at `1e-9` |
|---|---|---|---|---|
| `1 · n_edges` | 22 006 | 0 on all 8 meshes | `7e-11` | 0 |
| `4 · n_edges` | 74 398 | 0 | `1.8e-10` | 0 |
| `16 · n_edges` | 275 781 | 0 | `4.7e-10` | 0 |
| `64 · n_edges` | 1 082 993 | 0 | `2.8e-9` | bob 4 of 15 939, camelbad 7 of 33 958, camelbadnonuniform 4 of 27 066; 0 on the other five |

The failures at `1e-9` are paths that end `1e-9` to `3e-9` of their length from the vertex:
the accumulated effect of flips at the shape limit (SPEC N5), each of which may move an angle
by `2e-10`.

### Against the authors' program

`tests/reference.rs`, on the 15 inputs of `data/intrinsic-reference` (spot, bob, blub, bunny,
five camels, bump, two icospheres, a thin cylinder, a long torus, a sheared grid), with the
reference's Delaunay test: the same triangulation on 14 inputs (every halfedge matched by its
vertex pair and length), the same numbers of self-edges and double edges (long torus 6 and 16
pairs, camel_noise1 1 and 5, camelshotnoise 22 and 24, camel 0 and 2); on `blub` 68 halfedges
differ, within its 43 cocircular ties. Flip counts are equal on spot, bob, bunny, bump, the
icospheres, the cylinder and the grid and differ by at most 27 on the others (the edges are
numbered differently, so the queues differ). Lengths agree to `4e-12` relative or better,
except on `camelshotnoise`, where 12 edges differ by up to `1.9e-9`: there this crate's
lengths agree with the traced paths to `1.8e-12` and the reference's do not. The Laplacian
agrees with the reference's file to `3e-12` of its largest entry (`2e-7` on `camelshotnoise`,
which is the accuracy of that file by the reference's own check).

### Distances

`tests/distance.rs`: the first-order upwind solver of `ddg_core` on the input mesh and on the
intrinsic Delaunay triangulation of sheared lattices (all triangles obtuse), against the exact
distance, in units of the lattice width `h`:

| | mean error | largest error |
|---|---|---|
| intrinsic Delaunay triangulation | 0.03 to 0.30 | 0.06 to 0.96 |
| input mesh, seed at the corner the shear leans away from | 0.004 to 0.07 | 0.02 to 0.15 |
| input mesh, seed at another corner or in the middle | 0.4 to 6.1 (12.6 at 48 cells) | 1.3 to 16 (32) |

On the intrinsic triangulation the error does not grow with the resolution; on the input mesh
it does, except from the one corner where its long edges happen to point along the paths. On
curved shapes of good quality the two solutions differ by 0.02 to 0.09 `h` in the mean.

### Refinement at 25° and 30°

`refine` with the default options and with 30°, against the authors' program on the same
inputs (`data/intrinsic-reference/*/refine25|refine30`; it was not run on camel_noise2). Every
run ends with status `Done`. "Where reachable" is the smallest angle over the corners at
vertices whose angle sum is not below the bound. Times are one thread, fastest of five, on
the same machine for both programs (the reference's own timer around the operation).

| mesh | V | bound | inserted (reference) | smallest angle (reference) | where reachable | largest (reference) | left below the bound | ms (reference) |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| spot | 2930 | 25° | 255 (255) | 25.047° (25.047°) | 25.047° | 126.82° (126.82°) | none | 4.62 (13.48) |
| spot | 2930 | 30° | 827 (827) | 30.031° (30.031°) | 30.031° | 119.02° (119.02°) | none | 6.88 (18.16) |
| bob | 5344 | 25° | 93 (93) | 25.044° (25.044°) | 25.044° | 123.11° (123.11°) | none | 6.31 (15.93) |
| bob | 5344 | 30° | 1108 (1108) | 30.017° (30.017°) | 30.017° | 119.18° (119.18°) | none | 11.10 (25.65) |
| blub | 7106 | 25° | 1771 (1771) | 25.005° (25.005°) | 25.005° | 129.86° (129.86°) | none | 18.32 (49.20) |
| blub | 7106 | 30° | 3651 (3646) | 30.004° (30.004°) | 30.004° | 119.84° (119.84°) | none | 27.01 (72.46) |
| bunny | 14290 | 25° | 0 (0) | 35.383° (35.383°) | 35.383° | 92.04° (92.04°) | none | 15.65 (33.93) |
| bunny | 14290 | 30° | 0 (0) | 35.383° (35.383°) | 35.383° | 92.04° (92.04°) | none | 15.62 (38.18) |
| camel | 11381 | 25° | 3544 (3548) | 25.002° (25.002°) | 25.002° | 129.73° (129.73°) | none | 33.22 (97.58) |
| camel | 11381 | 30° | 8963 (8892) | 30.000° (30.000°) | 30.000° | 119.92° (119.92°) | none | 56.81 (152.88) |
| camelbad | 11392 | 25° | 4169 (4171) | 25.000° (25.000°) | 25.000° | 129.95° (129.95°) | none | 35.64 (109.43) |
| camelbad | 11392 | 30° | 10922 (11004) | 30.001° (30.001°) | 30.001° | 119.97° (119.97°) | none | 64.09 (169.06) |
| camelbadnonuniform | 9092 | 25° | 3756 (3759) | 25.001° (25.001°) | 25.001° | 129.94° (129.94°) | none | 30.69 (94.58) |
| camelbadnonuniform | 9092 | 30° | 9791 (9733) | 30.001° (30.001°) | 30.001° | 119.93° (119.93°) | none | 56.69 (160.10) |
| camel_noise1 | 11381 | 25° | 6698 (6679) | 25.000° (25.000°) | 25.000° | 129.98° (129.97°) | none | 47.70 (124.60) |
| camel_noise1 | 11381 | 30° | 15936 (15833) | 30.001° (30.000°) | 30.001° | 119.88° (119.73°) | none | 89.22 (231.40) |
| camel_noise2 | 11381 | 25° | 8368 (not run) | 25.001° | 25.001° | 129.87° | none | 51.89 |
| camel_noise2 | 11381 | 30° | 20232 (not run) | 30.000° | 30.000° | 119.88° | none | 108.01 |
| camelshotnoise | 11513 | 25° | 5660 (5630) | 0.042° (0.042°) | 25.002° | 129.73° (129.73°) | 22 corners, angle sum of their vertex below the bound | 42.33 (133.01) |
| camelshotnoise | 11513 | 30° | 15752 (16166) | 0.042° (0.042°) | 30.001° | 119.93° (119.93°) | 22 corners, angle sum of their vertex below the bound | 84.04 (226.26) |
| bump | 3329 | 25° | 0 (0) | 25.370° (25.370°) | 25.370° | 95.94° (95.94°) | none | 3.13 (6.38) |
| bump | 3329 | 30° | 13 (13) | 30.085° (30.085°) | 30.085° | 114.87° (114.87°) | none | 3.21 (8.92) |
| icosphere | 42 | 25° | 0 (0) | 55.569° (55.569°) | 55.569° | 68.86° (68.86°) | none | 0.04 (0.07) |
| icosphere | 42 | 30° | 0 (0) | 55.569° (55.569°) | 55.569° | 68.86° (68.86°) | none | 0.04 (0.07) |
| icosphere_jitter | 42 | 25° | 0 (0) | 26.098° (26.098°) | 26.098° | 102.06° (102.06°) | none | 0.04 (0.11) |
| icosphere_jitter | 42 | 30° | 2 (2) | 30.252° (30.252°) | 30.252° | 102.49° (102.49°) | none | 0.05 (0.17) |
| thin_cylinder | 52 | 25° | 66 (68) | 26.180° (26.180°) | 26.180° | 114.76° (114.76°) | none | 0.39 (1.08) |
| thin_cylinder | 52 | 30° | 150 (146) | 31.212° (30.230°) | 31.212° | 110.56° (110.56°) | none | 0.72 (1.76) |
| long_torus | 24 | 25° | 54 (54) | 25.112° (25.112°) | 25.112° | 107.90° (107.90°) | none | 0.26 (0.81) |
| long_torus | 24 | 30° | 72 (72) | 33.949° (33.949°) | 33.949° | 106.25° (106.25°) | none | 0.33 (0.99) |
| sheared_grid | 121 | 25° | 256 (170) | 9.462° (9.462°) | 25.058° | 102.34° (153.43°) | 2 corners, angle sum of their vertex below the bound | 0.80 (2.00) |
| sheared_grid | 121 | 30° | 395 (326) | 9.462° (2.514°) | 30.213° | 114.97° (163.51°) | 2 corners, angle sum of their vertex below the bound | 1.28 (4.91) |

- The bound is reached at every corner where a triangulation can reach it. What is left: the
  22 corners at the tips of the spikes of camelshotnoise (angle sum 0.042°) and the two acute
  corners of the sheared grid's boundary (9.46°), each reported with the reason. No triangle
  was set aside for another reason (`NeedleVertex`, `NotInsertable`) and no run reached the
  insertion limit.
- 30° costs 2 to 3 times the inserted vertices of 25° on the meshes that need refinement
  (spot 3.2, blub 2.1, camel 2.5, camelbad 2.6, camelshotnoise 2.8; bob 12, which needs almost
  none at 25°), for a largest angle of 120° instead of 130°. That is why the default is 25°.
- The counts are the reference's within 1 % on 12 of the 15 inputs; where there is a boundary
  with bad triangles the rules differ (SPEC N1): the sheared grid ends at 25.06° and 102°
  here, at 9.46° and 153° there.
- 2.5 to 3 times faster than the reference on the larger meshes.
- Termination: argued for closed surfaces up to 30°, observed with boundary and with spikes
  (SPEC N1, N4); the insertion limit is the net.

### Traces and overlay after refinement

| mesh | bound | edges traced | crossings | largest mismatch | trace ms | overlay vertices | overlay faces | overlay ms |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| spot | 25° | 1658 | 1460 | 3.5e-11 | 0.56 | 4645 | 8719 | 1.13 |
| spot | 30° | 4191 | 2812 | 3.3e-13 | 1.19 | 6569 | 12032 | 1.45 |
| blub | 25° | 10400 | 13041 | 1.8e-13 | 3.32 | 21918 | 35884 | 4.18 |
| blub | 30° | 17528 | 17165 | 2.3e-12 | 5.06 | 27922 | 45270 | 4.99 |
| camel | 25° | 21086 | 28433 | 4.7e-12 | 6.81 | 43358 | 68679 | 7.26 |
| camel | 30° | 40665 | 42046 | 3.3e-12 | 11.93 | 62390 | 96452 | 10.68 |
| camelbad | 25° | 24536 | 32360 | 1.1e-12 | 7.96 | 47921 | 75452 | 9.22 |
| camelbad | 30° | 48831 | 48891 | 7.9e-12 | 13.83 | 71205 | 109531 | 11.45 |
| camel_noise1 | 25° | 34479 | 30396 | 1.1e-12 | 10.55 | 48475 | 80881 | 9.09 |
| camel_noise1 | 30° | 67306 | 49369 | 5.1e-12 | 18.29 | 76686 | 123445 | 13.53 |
| camelshotnoise | 25° | 28005 | 33961 | 3.8e-11 | 8.58 | 51134 | 79275 | 7.78 |
| camelshotnoise | 30° | 62197 | 53027 | 2.7e-11 | 15.61 | 80292 | 122443 | 12.71 |
| sheared_grid | 25° | 638 | 3188 | 5.7e-14 | 0.25 | 3565 | 3956 | 0.16 |
| sheared_grid | 30° | 973 | 4170 | 9.5e-14 | 0.34 | 4686 | 5176 | 0.22 |

Every edge trace arrives (largest mismatch `8e-12` of the edge length, tolerance `1e-6`).
The overlay passes the checks of SPEC N12 on all 17 generated shapes and on the eight models
and bad meshes, for the Delaunay triangulation, after random flips and insertions, after
refinement and after the export.

### The export

34 of 50 refined triangulations (25 meshes, two bounds) are simplicial complexes already and
are handed out unchanged, among them spot, bob, blub, camelbad and camelbadnonuniform. The
others, with what splitting does to their angles, and what refinement between the rounds of
splits (`ExportOptions::refine`) recovers (`tests/export.rs`,
`the_export_keeps_the_angles_only_where_nothing_is_split`; angles over the corners at vertices
whose angle sum is not below the bound):

| mesh | bound | self-edges, edges in groups | splits and flips only: splits | smallest reachable, largest | corners below | with refinement between: rounds | ended refined | splits, inserted | smallest reachable, largest | corners below |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| torus 5x3 | 25° | 5, 35 | 30 | 25.938°, 121.24° | 0 | 1 | yes | 30, 0 | 25.938°, 121.24° | 0 |
| torus 5x3 | 30° | 0, 10 | 5 | 25.021°, 103.51° | 10 | 1 | yes | 5, 15 | 30.404°, 112.53° | 0 |
| long torus 8x3 R 4 r 0.2 twist 0.37 | 25° | 13, 88 | 95 | 9.564°, 144.21° | 65 | 3 | yes | 87, 35 | 25.743°, 115.12° | 0 |
| long torus 8x3 R 4 r 0.2 twist 0.37 | 30° | 4, 116 | 78 | 11.463°, 128.25° | 123 | 2 | yes | 69, 75 | 30.060°, 116.63° | 0 |
| needle tetrahedron h 8 lift 4 | 25° | 2, 4 | 8 | 6.418°, 130.41° | 8 | 4 | no | 24, 18 | 6.418°, 130.41° | 10 |
| needle tetrahedron h 8 lift 4 | 30° | 1, 4 | 5 | 6.418°, 132.77° | 10 | 2 | no | 13, 17 | 6.418°, 130.41° | 9 |
| tube 4x12 r 0.1 twist 0.37 stagger 0 | 25° | 0, 98 | 49 | 24.792°, 114.33° | 3 | 1 | yes | 49, 3 | 25.414°, 112.96° | 0 |
| tube 4x12 r 0.1 twist 0.37 stagger 0 | 30° | 0, 4 | 2 | 29.089°, 110.56° | 1 | 1 | yes | 2, 1 | 31.270°, 110.56° | 0 |
| tube 3x6 r 0.02 twist 0.37 stagger 1 | 25° | 65, 132 | 466 | 14.442°, 129.60° | 242 | 3 | yes | 402, 327 | 25.078°, 121.34° | 0 |
| tube 3x6 r 0.02 twist 0.37 stagger 1 | 30° | 65, 132 | 466 | 14.442°, 129.60° | 321 | 2 | yes | 273, 580 | 30.301°, 115.30° | 0 |
| camel_noise1 | 25° | 1, 12 | 11 | 6.377°, 133.73° | 17 | 4 | no | 22, 19 | 11.507°, 129.98° | 5 |
| camel_noise1 | 30° | 1, 10 | 8 | 8.471°, 119.88° | 16 | 4 | no | 21, 112 | 7.842°, 119.88° | 9 |
| camel_noise2 | 25° | 1, 128 | 89 | 7.662°, 130.88° | 100 | 6 | no | 148, 172 | 8.678°, 129.87° | 10 |
| camel_noise2 | 30° | 1, 88 | 53 | 7.659°, 124.20° | 88 | 6 | no | 127, 824 | 7.149°, 122.18° | 30 |
| camelshotnoise | 25° | 42, 134 | 299 | 0.014°, 179.96° | 195 | 0 | no | 299, 0 | 0.014°, 179.96° | 195 |
| camelshotnoise | 30° | 40, 158 | 242 | 0.014°, 179.96° | 294 | 0 | no | 242, 0 | 0.014°, 179.96° | 294 |

Thin tubes and handles: the bound comes back in one to three rounds. Spikes (a vertex with a
small angle sum: the needle tetrahedron, the noise camels): it does not, and cannot at a
vertex whose angle sum is below three times the bound; the damage is local (5 to 30 corners)
except on camelshotnoise, where the loops around its 22 spikes become slivers of 0.005° and
179.96° (SPEC N15). The export says so in `min_angle`, `max_angle` and `refinement`.

### End to end: the cut flow on a refined export

The two flat bad meshes of `ddg-variational-cuts/examples/bad_mesh_study.rs` at its size
(`FLAT_N = 32`), and its two measurements: redistancing from the exact distance to a circle
(mean and largest error in input mesh widths `h = √(area / vertices)`), and the circle law
(a circle of radius 0.3 shrinking under the length term, 96 steps to `t = 0.03`: relative
error of the mean radius, and standard deviation of the radius over its mean), with the
flow's default stages and with the paper's. Then the cut of the default stages carried back
to the input mesh. `tests/export.rs`, `the_cut_flow_runs_on_a_refined_export_and_its_cut_comes_back`.

| mesh | treatment | V | redistancing, circle: mean, max | circle law, default stages: radius error, roundness | circle law, paper's stages | cut on the input mesh: nodes, pieces | nodes on input edges / inside input faces |
|---|---|---:|---:|---:|---:|---:|---:|
| regular grid 32² | input mesh | 1089 | 0.057, 0.422 | +0.24%, 0.54% | -15.70%, 2.60% | | |
| sheared grid 153° | input mesh | 1085 | 1.772, 8.588 | +1.87%, 4.10% | the cut vanished | | |
| sheared grid 153° | flipped | 1085 | 0.059, 0.425 | +0.25%, 0.56% | -16.10%, 2.70% | 74, 188 | 22 / 52 |
| sheared grid 153° | refined 25° | 1213 | 0.059, 0.425 | +0.25%, 0.56% | -16.10%, 2.70% | 74, 188 | 22 / 52 |
| sheared grid 153° | refined 30° | 1213 | 0.059, 0.425 | +0.25%, 0.56% | -16.10%, 2.70% | 74, 188 | 22 / 52 |
| needles and caps | input mesh | 703 | 0.100, 0.365 | -0.47%, 1.30% | -39.10%, 4.50% | | |
| needles and caps | flipped | 703 | 0.098, 0.436 | -0.67%, 0.93% | -35.60%, 3.40% | 100, 112 | 54 / 46 |
| needles and caps | refined 25° | 2395 | 0.032, 0.146 | +0.21%, 0.16% | -11.70%, 0.44% | 85, 145 | 6 / 79 |
| needles and caps | refined 30° | 4812 | 0.025, 0.114 | +0.24%, 0.12% | -10.00%, 0.47% | 116, 176 | 6 / 110 |

- Both bad meshes go through refine, export, flow and back, and behave like the regular grid
  afterwards: the circle law holds to 0.2 % with the default stages; the paper's stages have
  their own error of 16 % on the regular grid too and reach 10 to 12 % on the refined needles.
- The sheared grid is a lattice: flipping alone makes it the regular grid's triangulation, and
  refinement adds only vertices on its ragged boundary (128), away from the circle.
- On the input meshes themselves: redistancing is off by up to 8.6 mesh widths on the sheared
  grid and the paper's stages lose the cut; on needles and caps the circle is 8 times less
  round than after refinement.
- The way back: every node of the cut is found on the input mesh within `3e-16` of where the
  positions of the exported vertices put it, the pieces of a segment add up to its length
  within `5e-13`, and the cut on the input mesh has the same radius and roundness as on the
  export. Most nodes lie *inside* input faces (79 of 85 on the refined needles).

### Distances after refinement

SPEC §5 item 5 expects distances to be better after refinement. Measured
(`tests/refinement.rs`, `distances_after_refinement`; first-order solver, sheared lattices,
errors at the input vertices in lattice widths): 1 to 7 on the input mesh, 0.05 to 0.26 on
the Delaunay triangulation, 0.07 to 0.26 after refinement to 25° or 30°. Refinement is as
good as flipping and not better (SPEC N16).

### The optimal Delaunay triangulation (paper Sec. 4.3, Fig. 14, 17)

`examples/odt_study.rs` on 12 meshes (logs in `runs_tmp/odt/`, not in the repository): flips
alone (iDT), refinement to 25° and 30° (iDR), the optimal Delaunay triangulation with its
target length chosen by bisection so that its vertex count is that of iDR at 25° and at 30°
(≈, within 1 %), and with the default options (target = the mean input edge). `s`: the time of
the operation, one thread, run alone (the bisection runs are not counted). Angles: the smallest,
the 1 % and 99 % quantiles and the largest corner, percent of corners below 30° and above 90°.
Areas: largest over smallest, coefficient of variation. `λ_max/λ_2`: condition number of the
pencil of the cotan Laplacian and the lumped masses (power iteration, inverse iteration with a
Cholesky factorization). Trace mismatch: the largest relative mismatch of an edge trace over the
input mesh (SPEC N8; the tolerance is `1e-6`). Not shown: camelbadnonuniform and camel_noise1
(like camelbad and camel_noise2), the sheared grid at 135° (all six are the regular grid of the
flips: nothing is split or moved).

| mesh | triangulation | V | s | min° | 1 % | 99 % | max° | < 30° | > 90° | area max/min | area cv | λ_max/λ_2 | trace mismatch |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| spot | iDT | 2930 | 0.002 | 11.428 | 19.98 | 97.66 | 132.56 | 5.13 | 4.23 | 1.840e2 | 0.644 | 1.188e5 | 7.9e-16 |
| spot | iDR 25° | 3185 | 0.005 | 25.047 | 27.94 | 106.28 | 126.82 | 2.29 | 4.83 | 3.065e2 | 0.718 | 1.737e5 | 3.5e-11 |
| spot | iODT ≈ iDR 25° | 3202 | 0.039 | 11.428 | 20.29 | 98.90 | 132.56 | 4.49 | 4.21 | 1.174e2 | 0.529 | 1.188e5 | 2.0e-14 |
| spot | iDR 30° | 3757 | 0.007 | 30.031 | 31.53 | 105.35 | 119.02 | 0.00 | 4.45 | 3.838e2 | 0.817 | 1.469e5 | 3.3e-13 |
| spot | iODT ≈ iDR 30° | 3792 | 0.064 | 11.428 | 21.28 | 99.06 | 132.56 | 3.28 | 4.34 | 8.511e1 | 0.458 | 1.189e5 | 2.1e-14 |
| spot | iODT default | 7688 | 0.483 | 12.446 | 27.90 | 98.13 | 118.51 | 1.49 | 4.29 | 4.257e1 | 0.326 | 1.190e5 | 2.3e-12 |
| bob | iDT | 5344 | 0.003 | 11.061 | 26.78 | 91.85 | 103.58 | 5.41 | 11.94 | 5.623e1 | 0.475 | 5.617e4 | 8.8e-16 |
| bob | iDR 25° | 5437 | 0.007 | 25.044 | 28.32 | 91.87 | 123.11 | 4.62 | 12.21 | 8.100e1 | 0.494 | 6.644e4 | 8.0e-14 |
| bob | iODT ≈ iDR 25° | 5402 | 0.024 | 11.061 | 26.81 | 91.86 | 111.00 | 5.35 | 11.87 | 5.087e1 | 0.446 | 5.617e4 | 5.7e-15 |
| bob | iDR 30° | 6452 | 0.012 | 30.017 | 31.27 | 91.87 | 119.18 | 0.00 | 8.75 | 1.067e2 | 0.624 | 6.324e4 | 8.1e-13 |
| bob | iODT ≈ iDR 30° | 6453 | 0.116 | 11.061 | 28.11 | 94.51 | 116.36 | 4.19 | 10.29 | 3.500e1 | 0.344 | 5.618e4 | 2.8e-14 |
| bob | iODT default | 12186 | 0.731 | 14.183 | 28.34 | 95.57 | 117.31 | 2.56 | 6.25 | 1.698e1 | 0.293 | 5.618e4 | 9.3e-12 |
| blub | iDT | 7106 | 0.004 | 6.921 | 12.10 | 97.98 | 158.97 | 11.82 | 4.64 | 8.205e2 | 0.714 | 3.783e5 | 2.3e-15 |
| blub | iDR 25° | 8877 | 0.019 | 25.005 | 26.86 | 111.97 | 129.86 | 3.77 | 6.88 | 4.213e2 | 0.904 | 4.023e5 | 1.8e-13 |
| blub | iODT ≈ iDR 25° | 8836 | 0.345 | 6.921 | 13.46 | 99.24 | 158.97 | 8.46 | 3.92 | 4.722e2 | 0.559 | 3.784e5 | 2.7e-13 |
| blub | iDR 30° | 10757 | 0.027 | 30.004 | 31.17 | 107.82 | 119.84 | 0.00 | 5.49 | 4.382e2 | 1.048 | 5.256e5 | 2.3e-12 |
| blub | iODT ≈ iDR 30° | 10768 | 0.783 | 7.682 | 15.83 | 99.16 | 158.97 | 6.49 | 3.86 | 3.914e2 | 0.483 | 3.784e5 | 3.5e-12 |
| blub | iODT default | 17885 | 1.499 | 7.682 | 20.17 | 97.31 | 158.97 | 3.49 | 3.46 | 1.943e2 | 0.344 | 3.785e5 | 1.3e-11 |
| camelbad | iDT | 11392 | 0.009 | 1.342 | 12.38 | 116.31 | 164.84 | 12.32 | 8.46 | 3.837e5 | 1.047 | 2.159e7 | 6.9e-14 |
| camelbad | iDR 25° | 15561 | 0.040 | 25.000 | 26.81 | 112.23 | 129.95 | 4.04 | 7.94 | 1.437e5 | 1.236 | 5.420e7 | 1.2e-12 |
| camelbad | iODT ≈ iDR 25° | 15692 | 0.801 | 1.497 | 15.01 | 109.54 | 164.84 | 7.32 | 6.42 | 7.093e4 | 0.535 | 2.159e7 | 2.1e-13 |
| camelbad | iDR 30° | 22314 | 0.072 | 30.001 | 31.18 | 106.60 | 119.97 | 0.00 | 5.81 | 9.990e4 | 1.332 | 5.124e7 | 7.9e-12 |
| camelbad | iODT ≈ iDR 30° | 22232 | 2.074 | 1.942 | 19.32 | 104.13 | 164.27 | 4.40 | 5.07 | 4.133e4 | 0.401 | 2.159e7 | 1.9e-12 |
| camelbad | iODT default | 26712 | 3.059 | 2.313 | 22.09 | 101.88 | 162.50 | 3.20 | 4.30 | 3.146e4 | 0.357 | 2.159e7 | 1.3e-12 |
| camel_noise2 | iDT | 11381 | 0.008 | 1.185 | 12.06 | 122.16 | 161.36 | 11.52 | 12.95 | 2.570e3 | 0.687 | 4.761e5 | 3.5e-15 |
| camel_noise2 | iDR 25° | 19749 | 0.055 | 25.001 | 26.52 | 110.96 | 129.87 | 4.22 | 8.07 | 4.115e3 | 0.783 | 2.302e6 | 7.5e-13 |
| camel_noise2 | iODT ≈ iDR 25° | 19589 | 0.840 | 1.558 | 19.22 | 106.88 | 158.70 | 3.72 | 5.82 | 5.263e2 | 0.352 | 4.419e5 | 1.4e-13 |
| camel_noise2 | iDR 30° | 31613 | 0.112 | 30.000 | 31.07 | 105.61 | 119.88 | 0.00 | 5.88 | 3.737e3 | 0.857 | 2.858e6 | 2.8e-12 |
| camel_noise2 | iODT ≈ iDR 30° | 31709 | 3.157 | 2.621 | 28.44 | 97.98 | 158.70 | 1.29 | 3.01 | 2.783e2 | 0.252 | 3.493e6 | 3.3e-12 |
| camel_noise2 | iODT default | 25951 | 2.779 | 2.421 | 24.73 | 101.11 | 158.70 | 2.03 | 3.96 | 3.676e2 | 0.285 | 4.333e5 | 4.3e-12 |
| camelshotnoise | iDT | 11513 | 0.008 | 0.026 | 10.40 | 117.04 | 179.71 | 13.46 | 6.23 | 1.025e7 | 0.998 | 1.792e10 | 2.1e-11 |
| camelshotnoise | iDR 25° | 17173 | 0.044 | 0.042 | 26.83 | 111.96 | 129.73 | 4.00 | 7.31 | 5.718e5 | 1.307 | 1.879e8 | 3.8e-11 |
| camelshotnoise | iODT ≈ iDR 25° | 17072 | 1.261 | 0.026 | 14.06 | 105.85 | 179.71 | 7.56 | 5.37 | 3.302e6 | 0.527 | 2.105e8 | 2.1e-11 |
| camelshotnoise | iDR 30° | 27265 | 0.090 | 0.042 | 31.09 | 106.55 | 119.93 | 0.01 | 5.66 | 6.234e5 | 1.562 | 2.558e8 | 2.7e-11 |
| camelshotnoise | iODT ≈ iDR 30° | 27257 | 3.428 | 0.026 | 20.83 | 101.49 | 179.71 | 3.96 | 4.27 | 3.302e6 | 0.399 | 9.643e8 | 1.0e-10 |
| camelshotnoise | iODT default | 23130 | 2.593 | 0.026 | 18.17 | 102.86 | 179.71 | 5.14 | 4.77 | 3.302e6 | 0.426 | 9.684e8 | 1.3e-10 |
| jitter 0.33 | iDT | 1089 | 0.001 | 14.730 | 22.85 | 107.73 | 125.34 | 4.52 | 7.50 | 9.204e0 | 0.342 | 1.507e3 | 7.7e-16 |
| jitter 0.33 | iDR 25° | 1214 | 0.002 | 25.028 | 27.06 | 106.99 | 123.77 | 3.12 | 7.46 | 1.281e1 | 0.379 | 3.157e3 | 1.4e-13 |
| jitter 0.33 | iODT ≈ iDR 25° | 1205 | 0.011 | 14.730 | 24.40 | 105.85 | 125.34 | 3.29 | 6.73 | 8.200e0 | 0.314 | 2.002e3 | 1.1e-14 |
| jitter 0.33 | iDR 30° | 1514 | 0.003 | 30.011 | 31.45 | 105.00 | 119.48 | 0.00 | 6.19 | 1.726e1 | 0.432 | 3.395e3 | 2.5e-13 |
| jitter 0.33 | iODT ≈ iDR 30° | 1501 | 0.032 | 17.794 | 27.35 | 103.58 | 125.34 | 1.80 | 5.57 | 7.374e0 | 0.315 | 2.369e3 | 2.0e-14 |
| jitter 0.33 | iODT default | 2469 | 0.093 | 20.917 | 32.24 | 98.20 | 118.68 | 0.46 | 3.54 | 5.554e0 | 0.258 | 3.515e3 | 4.2e-14 |
| sheared 165° | iDT | 1089 | 0.002 | 5.104 | 5.10 | 159.90 | 159.90 | 5.97 | 3.12 | 1.000e0 | 0.000 | 8.170e4 | 9.6e-15 |
| sheared 165° | iDR 25° | 1409 | 0.003 | 15.000 | 26.20 | 90.00 | 104.83 | 4.45 | 0.99 | 8.000e0 | 0.306 | 1.261e5 | 9.6e-15 |
| sheared 165° | iODT ≈ iDR 25° | 1409 | 0.010 | 15.000 | 15.00 | 95.76 | 114.91 | 3.25 | 3.45 | 4.500e0 | 0.323 | 1.098e5 | 1.2e-13 |
| sheared 165° | iDR 30° | 2223 | 0.006 | 15.000 | 31.15 | 107.02 | 110.83 | 0.02 | 3.28 | 1.600e1 | 0.716 | 1.740e5 | 3.3e-14 |
| sheared 165° | iODT ≈ iDR 30° | 2433 | 0.037 | 15.000 | 26.77 | 93.32 | 120.02 | 1.80 | 3.95 | 3.653e0 | 0.228 | 1.098e5 | 3.8e-12 |
| sheared 165° | iODT default | 1153 | 0.003 | 15.000 | 15.00 | 150.00 | 150.00 | 3.80 | 3.03 | 3.000e0 | 0.177 | 5.840e4 | 9.6e-15 |
| cells 1:16 | iDT | 16929 | 0.009 | 3.576 | 3.58 | 90.00 | 90.00 | 33.33 | 0.00 | 1.000e0 | 0.000 | 3.066e7 | 0.0e0 |
| cells 1:16 | iDR 25° | 56865 | 0.283 | 26.096 | 26.54 | 113.96 | 115.06 | 8.10 | 12.06 | 3.202e1 | 1.152 | 5.233e7 | 3.8e-14 |
| cells 1:16 | iODT ≈ iDR 25° | 50213 | 2.371 | 3.381 | 3.58 | 90.02 | 143.53 | 33.25 | 6.26 | 7.760e0 | 0.354 | 1.251e8 | 7.6e-11 |
| cells 1:16 | iDR 30° | 106529 | 0.494 | 30.373 | 32.01 | 103.36 | 113.65 | 0.00 | 6.37 | 2.397e1 | 0.861 | 4.770e7 | 6.0e-14 |
| cells 1:16 | iODT ≈ iDR 30° | 113573 | 16.310 | 9.070 | 22.90 | 110.54 | 159.57 | 3.04 | 3.49 | 1.270e1 | 0.250 | 5.782e7 | 9.6e-12 |
| cells 1:16 | iODT default | 50209 | 3.475 | 3.381 | 3.58 | 90.02 | 141.18 | 33.25 | 6.26 | 9.502e0 | 0.354 | 1.251e8 | 8.0e-11 |

**Sec. 4.3 reproduces.** At the vertex count of refinement the areas are far more even: the
coefficient of variation is lower on every mesh but the sheared 165° grid at 25° (0.32 against
0.31), by a factor of 1.1 (bob) to 5 (camelbadnonuniform: 0.58 against 2.9), and the largest
area ratio is lower on 6 of the 8 curved meshes at 25° and on 7 at 30° (not on camelshotnoise,
whose spikes keep their tiny triangles). The price is the skinny
angles: the smallest angle stays that of the flips (an input vertex cannot move, and iODT does
not split an edge for an angle), the 1 % quantile is 13 to 28° (3.6° on cells 1:16) where
refinement has 26 to 32°, 1 to 8 % of the corners are below 30° (33 % on cells 1:16 at 25°), and the largest angle stays near that of the flips (159°
on blub, 164° on the camels, where refinement has at most 130°). The 99 % quantile of the
largest angle is lower than refinement's on most meshes. The paper's Fig. 14 shows the same:
fewer outliers in area, a wider angle histogram than refinement.
**Fig. 17 in part.** `λ_max/λ_2` after iODT at the count of 25° is within 10 % of that of the flips
on 7 of the 8 curved meshes (the smallest elements are those at input vertices, which neither changes) and
lower than refinement's at the same vertex count in 13 of 18 comparisons (camelbad 2.2e7 against
5.4e7 and 5.1e7); it is higher on camelshotnoise (both counts), camel_noise2 at 30° and
cells 1:16 (both): spikes and thin input faces that iODT leaves. Fig. 17 compares with extrinsic
ODT, which is not implemented here.
**Time.** iODT costs 3 to 40 times the time of refinement at the same vertex count (bob 25°: 0.024
against 0.007 s; camelbad 30°: 2.1 against 0.07 s; cells 1:16 at 30°: 16 against 0.5 s), because
every iteration flips the whole triangulation twice and it needs 6 to 50 iterations to converge
(motion below `1e-3` of an edge). With the default target it converged on 10 of the 12 meshes;
camelbadnonuniform and camelshotnoise keep moving near their spikes (on camelshotnoise the
shape bound of SPEC N25 stopped 1640 vertex moves) and reached the limit of 50 iterations, as
did two of the matched runs (camelbadnonuniform, cells 1:16 at 30°).
**The default target length** (the mean input edge) adds 1 to 2 times the input's vertex count
on the curved meshes and nothing on the sheared lattices; on a regular mesh it splits every
diagonal and lowers the smallest angle (regular grid 45° to 39°, icosphere 54° to 29°,
`flow_on_bad_meshes.rs quality`). It is the knob for the vertex count, not a quality setting.

### The optimal Delaunay triangulation under the cut flow

`examples/flow_on_bad_meshes.rs` with `TREAT=odt` (default options) and `TREAT=odt15` (target 1.5
times the mean input edge), all six treatments rerun with this code (logs in
`runs_tmp/flow/`, not in the repository), default stages. The user's criteria are robustness,
correctness and speed.

| measure | flips (default) | refined 25° | iODT, target 1.5 | iODT, default |
|---|---|---|---|---|
| vertices added, real meshes (9) | 0 to 133 | 93 to 8 688 | 964 to 5 155 | 4 758 to 16 686 |
| building the triangulation, camelbad (s, alone) | 0.027 | 0.072 | 0.83 | 2.9 |
| flow step, camelbad / camel_noise2 / bob (ms, alone) | 15 / 15 / 8 | 22 / 28 / 9 | 19 / 19 / 10 | 38 / 33 / 15 |
| evaluation, largest relative error of the speed: sheared 165° / needles / regular grid | 0.028 / 0.061 / 0.031 | 0.028 / 0.029 / 0.031 | 0.028 / 0.020 / 0.031 | 0.028 / 0.020 / 0.022 |
| redistancing, largest error in widths: sheared 165° / needles / cells 1:16 | 0.48 / 0.44 / 1.02 | 0.59 / 0.27 / 0.36 | 0.48 / 0.19 / 1.02 | 0.50 / 0.19 / 0.98 |
| circle law, needles (radius error, roundness) | −0.67 %, 0.93 % | +0.20 %, 0.16 % | +0.23 %, 0.14 % | +0.22 %, 0.14 % |
| real meshes: final energy on the input mesh, bob / blub / camelbadnonuniform | 0.942 / 0.922 / 1.482 | 0.941 / 0.937 / 1.484 | 1.022 / 0.990 / 1.479 | 0.944 / 0.921 / 1.427 |
| way back, largest mismatch of a segment | `3e-13` (`3e-11` camelshotnoise) | `5e-12` (`2e-10`) | `5e-11` | `1.2e-10` |
| failures | none | none | none | latlon 41×82, `α_L = 3`: the run stalls at the start energy (10.17 against 9.77) |

- **Correctness.** Where flips leave a mesh well shaped (lattices, sheared grids, the
  icospheres) the target 1.5 splits nothing and iODT is the flip triangulation; where the input
  has long needles (needles and caps) both targets make the flow as accurate as refinement does
  (circle law 0.2 % instead of 0.7 %, evaluation and redistancing better than refinement). The
  default target, with 2 to 3 times the vertices, makes the evaluation on flat meshes more
  accurate by the finer mesh alone, and on the latitude–longitude spheres less accurate than
  flips (`s` up to 4.5e-3 against 1.1e-3), though more than refinement.
- **Robustness.** No error, no failed trace or way back on any mesh. One stall of the default
  target on a latitude–longitude sphere (above) that flips and refinement do not have; on the
  real meshes the target 1.5 stopped the flow earlier, at a higher energy, on bob and blub (the
  stop rule after a rise of the energy, SPEC D27), and the default target did better than all
  others on camelbadnonuniform. The runs on these meshes are sensitive to the triangulation
  (ROADMAP: "chaotic"), so single energies are not a ranking.
- **Speed.** iODT is the slowest treatment to build (up to 4.8 s alone with the default target
  and 1.2 s with the target 1.5, where flips take at most 0.07 s) and, with the default target,
  makes every flow step 1.7 to 4 times as slow as on the flips; the target 1.5 adds 0 to 70 %.
- **Recommendation: keep Delaunay flips as the default** of `ddg_uv::cutting`. iODT does not
  improve robustness or correctness over flips where refinement does not already, and costs time;
  it is offered as `IntrinsicOptions::optimal_delaunay` (`IntrinsicOptions::optimized()`) for
  callers who want even triangle areas.

### Remeshing for the cut flow (cut paper Sec. 6.1): does more resolution help?

The paper of the cut flow (Sharp, Crane 2018, Sec. 6.1, Fig. 6) optionally remeshes the input,
finer where `|K|` is large, because the cut can be no finer than the mesh. On an intrinsic
triangulation vertices can be added without changing the surface or the caller's mesh, and the
cut is carried back. `crates/ddg-uv/examples/remesh_study.rs` (SPEC N28) runs the flow on

| treatment | triangulation the flow runs on |
|---|---|
| `input`, `flip` | the caller's triangles; Delaunay flips (the default of `ddg_uv::cutting`, through its `Session`) |
| `u2`, `u4` | every edge longer than 1/2 (1/4) of a base length split at its midpoint, flips, repeated (9 to 40 times the vertices) |
| `ka2`, `ka4` | curvature-adaptive as in the paper: the longest edge of a face whose `\|K\| A` exceeds a bound split, flips, repeated; the bound chosen so that as many vertices are added as by `u2` (`u4`) |
| `odt2`, `odt4` | the optimal Delaunay triangulation with target 1/2 (1/4) of the base length (through the `Session`) |
| suffix `S` | the flow uses per face the curvature density of the input faces under it ("spread", `CurvatureSource::PerFace`) instead of the Gaussian curvature of the triangulation, which on an intrinsic triangulation sits only at the input vertices (an inserted vertex is flat) |

Meshes: the control meshes of spot, bob, blub (112 to 334 vertices, "coarse"), spot, bob, blub
("fine", 2 930 to 7 106; blub at the 1/2 level only), camel, camelbad, camelbadnonuniform,
camelshotnoise (9 092 to 11 513, 1/2 level only; `u4` of camel has 418 000 vertices). Base
length: the mean input edge; for the camels the side of the equilateral triangle of the mean
face area (`TARGET=area`), because on the non-uniform camels the median edge is 0.77 to 0.91 of
the mean and either edge statistic follows how the area is cut into slivers. Cases: Hencky and
Dirichlet, the automatic length weight of the bindings (ripple wavelength 0.25, from the flip
triangulation) and 1/2, 1/4 of it, two planar starts: 12 per mesh, 1 260 runs, default stop
rule (at most 1 000 steps). Every final cut is carried back to the input faces and evaluated on
one yardstick per mesh, the regular 1:n² subdivision of every input face (`n` odd, 9 to 33,
up to 2 million faces): with the curvature density of the input face under each face
("spread": the problem the caller's mesh poses, at fine resolution) and with the
subdivision's own Gaussian curvature ("polyhedral": cones at the input vertices). A third
measure is the energy of the final level-set function at the input vertices, evaluated on the
input mesh (the measure of the section above; it cannot see detail between input vertices).
Runs are deterministic (a rerun of the control meshes agreed to the last digit).

**Convergence**: median steps; stops EnergyStalled / EnergyRose / MaxSteps of runs; median
motion of the cut over the last 10 steps (input edges per step). No run of any treatment ends
`Converged` (SPEC D31): the runs at rest end by the energy rule.

| treatment | coarse | fine | camels |
|---|---|---|---|
| input | 70; 9/26/0 of 36; 0.066 | 216; 30/6/0 of 36; 0.089 | 259; 25/23/0 of 48; 0.124 |
| flip | 82; 5/30/0 of 36; 0.070 | 208; 28/8/0 of 36; 0.087 | 284; 25/23/0 of 48; 0.130 |
| flipS | 75; 8/28/0 of 36; 0.067 | 183; 29/7/0 of 36; 0.103 | 223; 34/14/0 of 48; 0.109 |
| u2 | 95; 12/24/0 of 36; 0.071 | 476; 28/4/4 of 36; 0.069 | 579; 17/12/19 of 48; 0.084 |
| u2S | 122; 14/21/0 of 36; 0.071 | 582; 27/3/6 of 36; 0.060 | 584; 21/10/17 of 48; 0.083 |
| u4 | 176; 22/14/0 of 36; 0.040 | 553; 15/1/8 of 24; 0.040 | – |
| u4S | 204; 14/22/0 of 36; 0.039 | 402; 12/5/7 of 24; 0.033 | – |
| ka2 | 135; 10/25/0 of 35; 0.050 | 373; 26/3/7 of 36; 0.034 | 682; 19/17/12 of 48; 0.055 |
| ka2S | 138; 16/20/0 of 36; 0.052 | 384; 25/4/7 of 36; 0.037 | 628; 15/12/21 of 48; 0.049 |
| ka4 | 171; 22/13/0 of 35; 0.034 | 410; 17/1/6 of 24; 0.022 | – |
| ka4S | 224; 19/17/0 of 36; 0.031 | 382; 18/1/5 of 24; 0.019 | – |
| odt2 | 90; 13/23/0 of 36; 0.076 | 550; 28/4/4 of 36; 0.063 | 550; 21/12/15 of 48; 0.085 |
| odt4 | 110; 15/21/0 of 36; 0.067 | 361; 19/2/3 of 24; 0.062 | – |

A finer triangulation takes 1.2 to 3 times the steps, stops less often on a rising energy and
moves less at the end (per input edge), but reaches the step limit in up to 44 % of the camel
runs and never converges by motion: it does not make the flow converge.

**Quality** against `flip` of the same case: median relative difference in %, cases lower than
`flip` / cases, cases more than 10 % higher (spread / polyhedral yardstick / on the input mesh).

| treatment | coarse | fine | camels |
|---|---|---|---|
| input | +0.4, 13/35, 6 / +0.1, 16/35, 7 / −0.1, 21/35, 3 | +0.0, 16/36, 0 / −0.0, 19/36, 0 / −0.0, 19/36, 0 | +0.2, 14/48, 4 / +0.8, 10/48, 7 / +0.0, 20/48, 4 |
| flipS | +0.0, 16/35, 10 / −0.0, 18/35, 8 / +0.0, 15/35, 8 | +0.0, 13/36, 3 / +0.0, 13/36, 3 / +0.0, 10/36, 3 | +1.4, 16/48, 11 / +1.9, 13/48, 11 / +1.2, 16/48, 12 |
| u2 | +1.1, 12/35, 14 / −6.9, 27/35, 7 / +9.1, 9/35, 17 | +0.0, 18/36, 5 / −0.0, 22/36, 5 / +0.0, 18/36, 5 | +0.2, 11/48, 7 / +0.4, 14/48, 7 / +0.2, 11/48, 8 |
| u2S | **−2.5, 25/34, 3** / −1.8, 25/34, 3 / +0.2, 15/34, 11 | −0.1, 21/36, 4 / −0.0, 23/36, 4 / −0.0, 22/36, 4 | +0.3, 14/48, 5 / +0.2, 19/48, 5 / +0.3, 11/48, 5 |
| u4 | +4.5, 10/35, 16 / −1.8, 25/35, 10 / +9.3, 11/35, 17 | +0.6, 4/24, 7 / +0.4, 5/24, 7 / +0.6, 3/24, 7 | – |
| u4S | **−2.1, 29/35, 1** / −3.1, 30/35, 1 / +0.0, 16/35, 9 | +0.4, 6/24, 7 / +0.3, 4/24, 7 / +0.3, 4/24, 7 | – |
| ka2 | +1.7, 10/34, 12 / −7.1, 27/34, 5 / +7.3, 11/34, 17 | +0.2, 10/36, 5 / +0.0, 17/36, 4 / +0.2, 12/36, 4 | +1.0, 4/48, 9 / +0.8, 9/48, 9 / +1.5, 4/48, 12 |
| ka2S | **−2.9, 26/35, 3** / −2.5, 25/35, 5 / +0.1, 17/35, 11 | +0.0, 16/36, 4 / +0.0, 14/36, 4 / +0.0, 13/36, 5 | +0.6, 13/48, 10 / +0.6, 14/48, 9 / +0.7, 9/48, 11 |
| ka4 | +4.9, 8/34, 16 / −6.7, 26/34, 6 / +10.7, 9/34, 17 | +0.6, 2/24, 8 / +0.5, 3/24, 8 / +0.6, 2/24, 8 | – |
| ka4S | **−2.1, 27/35, 1** / −2.7, 26/35, 1 / +0.6, 17/35, 12 | +0.3, 4/24, 8 / +0.2, 5/24, 8 / +0.2, 3/24, 7 | – |
| odt2 | +1.0, 9/35, 14 / −6.4, 27/35, 6 / +12.4, 9/35, 18 | −0.0, 22/36, 4 / −0.1, 22/36, 4 / −0.0, 21/36, 4 | +0.2, 15/48, 5 / +0.2, 18/48, 6 / +0.1, 20/48, 8 |
| odt4 | +6.1, 11/35, 16 / −4.6, 24/35, 7 / +6.9, 9/35, 17 | +0.1, 7/24, 4 / +0.0, 10/24, 4 / +0.1, 8/24, 4 | – |

- **Native curvature on a refined triangulation** (`u`, `ka`, `odt` without `S`) poses another
  problem: the curvature of the input sits at its vertices as cones. It wins on the polyhedral
  yardstick of the coarse meshes (−2 to −7 %), which is that problem, and loses on the spread one
  with outliers up to +1 160 % (blub_control). On the camels'
  dirichlet/α×0.25/start 0 case every refined run, with or without `S`, ends in another basin
  (energy 55 against 14.4, a shorter cut with more distortion).
- **Spread curvature** (`S`) removes those outliers and is the only refinement that is better
  than `flip` on the spread yardstick of the coarse meshes (−2 to −3 %, 25 to 29 of 35 cases,
  up to −40 % on spot_control). On the fine meshes it is even (median −0.1 to +0.4 %; the 1/4
  level is above `flip` in 18 of 24 cases, by more than 10 % in 7, 5 of them at α×1/4 and up to
  +171 % on bob: the finer mesh resolves the ripples of a small length weight); on the camels it is worse (+0.3 to +0.6 %,
  13 to 14 of 48 below). `flipS`, the spread without refinement, is −10.7 % on camelshotnoise
  (whose input has spikes) and +11 % on camelbad. Curvature-adaptive refinement is not better
  than uniform refinement at the same vertex count.
- On the input mesh (third measure) no refinement is clearly better than `flip` (medians −0.1
  to +12 %).

**Exact laws** (`remesh_study laws`; flat circle: radius error after 96 fixed steps, whose
time-step error alone is about +2e-3, SPEC D30, and roundness; latitude on icosphere(4): the ODE
after 50 steps, and the rest latitude from above and below; `ka` is left out, the curvature is
uniform):

| treatment | V (grid / ico4) | circle radius, roundness | latitude ODE | rest latitude ψ0 = 0.5, 1.1 |
|---|---|---|---|---|
| input | 1 089 / 2 562 | +1.8e-3, 3.7e-3 | +2.19e-3 | −1.05e-2, +1.11e-2 |
| flip | 1 089 / 2 562 | +1.4e-3, 2.0e-3 | +2.19e-3 | −1.05e-2, +1.11e-2 |
| u2 / u2S | 9 659 / 25 602 | +2.6e-3, 3.4e-4 | +2.13e-3 / +1.66e-3 | −8.7e-3, +8.4e-3 / −1.14e-2, +8.2e-3 |
| u4 / u4S | 38 142 / 102 402 | +2.8e-3, 1.6e-4 | +2.11e-3 / +1.91e-3 | −5.2e-3, +6.8e-3 / −1.10e-2, +1.08e-2 |
| odt2 | 8 884 / 25 422 | +2.6e-3, 1.8e-4 | +2.13e-3 | −8.8e-3, +8.2e-3 |
| odt4 | 11 979 / 28 182 | +2.7e-3, 2.5e-4 | +2.15e-3 | −8.2e-3, +8.1e-3 |

Every treatment meets the laws; refinement makes the circle 6 to 12 times rounder (its radius
error is the time step's) and, with native curvature, halves the known bias of the rest latitude
(ROADMAP); the spread keeps the bias of the input (its density is the coarse polyhedron's).

**Cost** (one thread, run alone, `remesh_study speed`; build = triangulation and export, the
whole `Session` for `flip` and `odt`; step = median of steps 5 to 25, Hencky):

| treatment | spot: V, build s, ms per step | camel: V, build s, ms per step |
|---|---|---|
| input | 2 930, 0, 3.1 | 11 381, 0, 14.1 |
| flip | 2 930, 0.016, 3.1 | 11 383, 0.076, 19.2 |
| u2 | 29 635, 0.12, 42 | 118 232, 0.73, 214 |
| u4 | 118 174, 0.56, 209 | – |
| ka2 | 30 690, 0.62, 45 | 123 107, 5.1, 233 |
| ka4 | 124 098, 3.2, 199 | – |
| odt2 | 29 056, 5.3, 43 | 112 440, 25.7, 220 |
| odt4 | 32 230, 0.32, 51 | – |

A step costs 11 to 70 times the step on the flips, and the runs take 2 to 3 times the steps;
`S` costs nothing per step. A first attempt at the 1/4 level on camel (418 000 vertices, 16
runs at a time) was stopped at 14 GB.

**Robustness.** 2 of the 900 way backs from a refined triangulation failed (`bob_control`,
`ka4`, Dirichlet, α×1, start 1; and `ka2`, α×1/2, start 1), none on any other mesh or treatment;
the largest relative mismatch of a segment in the other direct runs was `8e-8`. Reproducer
(`remesh_study repro bob_control ka4 dirichlet 1 1`, deterministic): the start plane z = 0 is
bob's mirror plane and contains input vertices and edges, the flow returns its start state, and
two segments of the cut run *along* input edge (300, 298) (between input faces 575 and 646),
near a vertex that the refinement inserted on that edge: both ends lie on the edge (barycentric
coordinates 0 to 2e-16). The trace of such a segment ends on the other side of the edge than
the place of its end; a larger tolerance (up to 0.1) does not help, because the test is which
input face the end lies in, not a distance (ka4: segments 524 and 528 in faces 22094, 22104;
ka2: 147 and 209 in faces 2862, 4391). **Fixed since** (SPEC N29): the ends were returned as
points of whichever of the two input faces rounding chose; a point query now returns a place
within `1e-12` of an input edge as a point of that edge, which belongs to both faces, and a
segment along a side of the triangulation that lies along an input edge is the piece of that
edge between its ends instead of a trace. Both reproducers now come back (largest relative
mismatch `1.4e-11` and `7.4e-12`); `tests/way_back_along_input_edges.rs` covers the class on
strips (interior and boundary input edges, segments between two crossings of the edge, from a
vertex on it, along a side on it). The `ka` builds read the curvature under face centroids
through `point_on_input`, so on bob's control mesh, whose mirror plane holds such centroids, the
`ka4` build now has 12 452 vertices instead of 12 454; the tables above are from before. It is
not the class of the open flipped
UV triangles of the iODT spot atlas (`blender/README.md`): there the way back succeeds, and the
flipped triangles come after the insertion; what the two share is a finer cut on the input
faces (from iODT with target 1/2 the cut comes back in 3 times the pieces of the flips' cut, median 1 142
against 361, and the bend merge at 0.1 removes about 570 bends where it removes 20 from the
flips' cut).

**Conclusion.** Remeshing does not improve the convergence of the cut flow (more steps, no run
converges by motion) and does not improve its quality on meshes that already resolve the
surface (spot, bob, blub, the camels): there it is even or worse, at 11 to 70 times the cost
per step. Only on very coarse meshes (a few hundred vertices) does refinement with the curvature
spread from the input faces lower the energy, by 2 to 3 % in the median. The paper's
curvature-adaptive refinement is no better than uniform refinement at the same vertex count,
and refinement with the curvature at the input vertices only (what an intrinsic or a planar
subdivision gives) solves a different problem, with large outliers. **Recommendation: no
remeshing in `ddg_uv::cutting`** (no option was added; the default stays Delaunay flips). A
caller with a very coarse mesh is better served by subdividing it before the flow.

## Verification

- `tests/invariants.rs`: a fresh triangulation is the input mesh, bit for bit; after random
  flips and after Delaunay flipping: strict triangle inequalities, total area (`1e-12`
  relative after Delaunay, `1e-9` after random flips), angle sums and signposts (`1e-11` rad,
  `1e-9`), Gauss–Bonnet (`1e-10`); every signpost leads to its vertex on the triangulation
  (exactly) and on the input mesh (`1e-9 ℓ`), and the traced path has the edge's length in
  space (`1e-10`); on flat shapes every edge is the straight segment between its end points
  (`1e-12`) and its crossings lie on it (`1e-11`).
- `tests/flip.rs`: what a flip changes; twice restores length and signposts to
  `1e-14 · c · ℓ` and `1e-13 · c²` rad, `c` the condition `ℓ_max / h_min` of the two triangles;
  every refusal, with its message; a refused flip changes nothing.
- `tests/delaunay.rs`: Delaunay condition and cotan weights afterwards; independence of the
  order of flips up to ties (three random orders per shape against the queue); empty
  circumcircles on convex flat shapes; the flip limit; the reference's test; self-edges, double
  edges and vertices of degree one on a long torus, a helical tube and a tilted needle
  tetrahedron.
- `tests/tracing.rs`: paths against straight segments on flat convex shapes, over the input
  and over the triangulation, from vertices, edge points and face points (`1e-12`); boundary
  stops; paths along edges; the way back on curved shapes; invalid queries.
- `tests/operators.rs`: the Laplacian of a fresh triangulation is `ddg_core`'s bit for bit;
  entry by entry against the edge weights; double edges and self-edges; Rippa's theorem (no
  Delaunay flip increases a Dirichlet energy); known quality numbers.
- `tests/distance.rs`, `tests/models.rs`, `tests/reference.rs`: the tables above.
- `tests/bad_input.rs` and the unit tests of `src/edge_trace.rs`: errors name the element; a
  damaged signpost or length is a failed trace; no NaN.
- A mutation check (24 one-line changes to the crate, each run against the test suite):
  all 24 make at least one test fail (one survived at first, the rule of SPEC N8 that a path
  is followed for the tolerance beyond its length; a unit test for it was added).

- `tests/insertion.rs`: insertions and splits keep area, angle sums, signposts and
  Gauss–Bonnet; a vertex is where its coordinates say (`1e-15` on flat shapes); snapping,
  also in a needle; every refusal; point queries are inverse to each other (`1e-9` of the
  mean edge) after flips, insertions and refinement; invalid points are errors.
- `tests/overlay.rs`: the invariants of SPEC N12 on fresh, Delaunay, randomly flipped,
  inserted-into, refined and exported triangulations, on the generated shapes and on the
  models and bad meshes; segments as paths on the input mesh, in faces that use a vertex
  twice too; segments that share a node meet in one point.
- `tests/refinement.rs`: SPEC §5 item 7 (every corner below the bound is listed with a reason
  that holds; Delaunay; no boundary edge encroached, so every cotan weight is non-negative,
  the boundary's included; consistency; every edge traces); determinism; the boundary rule;
  angle sums below, at and above the bound; limits and invalid options; the reference's counts
  and angles.
- Unit tests of `src/refine.rs`: the three rules at their thresholds (a point just inside and
  just outside the circle on a boundary edge, also from a neighbouring face; the needle guard
  against a traced loop; split distances at a boundary vertex of 20°).
- `tests/oracle.rs`: the 24 fixtures of the Lean oracle (`lean/oracle/fixtures/intrinsic`),
  an implementation of single operations written from the paper and `SPEC.md` without this
  crate's code: flips with the triangulation before and after, Delaunay flipping (also with
  self-edges and double edges), traces from vertices, the trace of every halfedge, and
  circumcentres with their traces; 23 000 comparisons, largest difference `4e-14` (SPEC N17);
  and the six `it-move-*` fixtures (an insertion or a split, flips, moves of the new vertex: its
  place, the triangulation and the update of the optimal Delaunay triangulation after every
  move, every refusal with the admissible fraction; SPEC N26): 7 294 comparisons, largest
  difference `9e-14`.
  And the oracle's readings where no fixture exists: a direction equal to a signpost, the
  first halfedge of a vertex, the lengths of an edge split.
- `tests/export.rs`: the export is a simplicial mesh of the same surface with the same ids,
  lengths and places; reproducible bit for bit; a stale state is refused; the angles before
  and after; a cut (a closed curve, and a network with a triple point) comes back as pieces
  on input faces; the end-to-end run above.

- `tests/reposition.rs` (SPEC §5 item 10): random moves of every inserted vertex on the 17
  generated shapes, 3 rounds with Delaunay flips between: items 1 to 3 (area `1e-11`, angle sums
  and signposts `1e-9`, Gauss–Bonnet, every edge traces to its vertex with `1e-8` and has its
  length in space), the overlay has the input's area, point queries are inverse (`1e-9` of the
  mean edge), input vertices stay; about half of the accepted moves are beyond the paper's
  height bound; on a flat grid the moved vertex is at its old place plus the displacement
  (`1e-14`) and its signposts are the world angles; a move and the move back restore lengths
  (`1e-12`) and signposts (`1e-10`); the ODT update against coordinates (`1e-14`) and against a
  central difference of Chen and Holst's energy (`1e-7`); every refusal changes nothing.
- `tests/odt.rs` (§5 item 11): iODT on the 17 shapes with 1, 3 and 12 iterations, both boundary
  rules: Delaunay, items 1 to 3, overlay, point queries, input vertices stay; on the models and
  bad meshes (8 iterations): Delaunay, consistency, traces; determinism and the prefix property
  of the iterations; convergence and even areas on a jittered grid (area cv 0.17 against 0.29
  for flips and 0.31 for refinement); fixed and sliding boundary vertices; limits; invalid
  options change nothing.

Not verified: termination of refinement beyond what SPEC N1 and N4 argue (it is observed);
the export's refinement rounds on meshes other than those in the table; the way back on
curved shapes is checked by the trace's own end-point check and by round trips, not against
positions (the position checks are on flat shapes); cuts whose nodes fall on vertices of the
export (the flow keeps crossings off vertices); the Lean proofs of SPEC §7; the oracle has no
fixtures for insertion, split (except as the start of a move), point queries, overlay, refinement,
export and the optimal Delaunay iteration as a whole; no mutation check was run on the second
half or on the moves.

## Formula ids

| id | function | formula |
|---|---|---|
| `it-angle-normalization` | `formulas::angle_scale` | `2π/Θ` at an interior vertex, `π/Θ` at a boundary vertex |
| `it-signpost-update` | `formulas::signpost_update` | `φ_ik = φ_ij + s θ_i^{jk}`, reduced to `[0, 2π)` at an interior vertex |
| `it-unfold-diamond` | `formulas::unfold_diamond` | the two triangles on an edge in the plane: `x = (a² + (b − c)(b + c)) / 2a`, `h = 2A/a` |
| `it-flip-length` | `Diamond::flipped_length` | `ℓ_kl² = (x_k − x_l)² + (h_k + h_l)²` |
| `it-flip-convexity` | `Diamond::flipped_areas` | `x_k h_l + x_l h_k > 0` and `(a − x_k) h_l + (a − x_l) h_k > 0` |
| `it-delaunay-condition` | `formulas::delaunay_angle_excess` | `θ_k + θ_l − π ≤ 0` |
| `it-cotan-weight` | `formulas::opposite_cotan_sum` | `cot θ_k + cot θ_l`, with the sign of `π − (θ_k + θ_l)`; this is twice the cotan weight `w = ½ (cot θ_k + cot θ_l)` of the edge, which is what every reported weight is |
| `it-unfold-apex` | `formulas::unfold_apex` | `c = a + α d + β d⊥`, `α = (ℓ_ab² + ℓ_ca² − ℓ_bc²)/2ℓ_ab²`, `β = 2A/ℓ_ab²` |
| `it-ray-edge-intersection` | `formulas::axis_crossing` | `t = y_tail / (y_tail − y_head)`, `x = x_tail + t (x_head − x_tail)`: the intersection of the line through the edge with the x-axis, the ray's line in the unfolded layout. The point is on the edge (`t ∈ [0, 1]`; `t = 0` exactly if the tail is on the axis, `t = 1` exactly if the head is) if the ends are on opposite sides of the axis or on it (`y_tail ≤ 0 ≤ y_head` or the reverse, not both zero). It is ahead of the entry crossing on the ray if the unfolded triangle `(up, down, apex)` of the entry edge, with `up` above and `down` below the axis, is counter-clockwise (signed double area `D > 0`): leaving through `down–apex` (apex above the axis) `x_exit − x_entry = t D / (y_up − y_down) > 0`, through `apex–up` (apex below) `(1 − t) D / (y_up − y_down) > 0` |
| `it-barycentric-coordinates` | `formulas::barycentric` | ratios of signed areas |
| `it-point-distances` | `formulas::corner_distances` | `d_k = |Σ_j b_j v_j − v_k|` with the triangle laid out in the plane, for barycentric coordinates that sum to one (the function does not read `b_0`: it places the point at `b_1 v_1 + b_2 v_2`, with `v_0` at the origin): the lengths of the edges of an inserted vertex, and of an edge split with the point on a side (on side `ij`: `d_i = ℓ_ij b_j`, `d_j = ℓ_ij b_i`) |
| `it-circumcenter` | `formulas::circumcenter_barycentric` | weights `a_k² (a_{k+1}² + a_{k+2}² − a_k²)`, `a_k` the side opposite corner `k`, normalized |
| `it-move-length` | `formulas::moved_length` | `ℓ' = √((ℓ − r cos α)² + (r sin α)²) = √(r² + ℓ² − 2rℓ cos α)` (Alg. 10) |
| `it-move-area` | `formulas::moved_double_area` | `2A' = 2A − v × (k − j)`: the area of a triangle of the star after its vertex moved by `v`; the areas of a closed star add up to the same total |
| `it-circumcenter-vector` | `formulas::weighted_circumcenter_vector` | `A (o − a) = (A ℓ_ab / 2, ℓ_ab (ℓ_ac² + ℓ_bc² − ℓ_ab²) / 8)` in the frame of side `ab` |
| `it-odt-update` | `formulas::odt_displacement`, `formulas::weighted_barycenter_vector` | `α Σ |T| (c_T − x) / Σ |T|` (Chen and Holst 2011, Eq. 4.13); proved: it is `−3/(2|ω|)` times the gradient of the planar ODT energy of the star (their Eq. 4.5) |

From `ddg-core`, on the triangulation: `cotan-laplacian`, `lumped-mass`, `angle-defect`,
`boundary-geodesic-curvature`.

## The reference implementation

The authors' demo (`navigating-intrinsic-triangulations-demo`, © 2019 Nick Sharp) and
geometry-central (© 2017–2019 Nicholas Sharp and the geometry-central contributors), both MIT,
were read for what the paper leaves out; the licence notices are in [`LICENSES/`](LICENSES).
No code was translated. Taken from the reference: the normalization of angles at boundary
vertices, the cotan-weight form of the Delaunay test with its tolerance (as an option), the
order of the flip queue, the flag for edges that still lie along an input edge. Different by
design: the tracing (one unfolding, no projections; SPEC N7), the check of an edge trace
(SPEC N8), the length of a flipped edge and the refusal of flips that are too thin (SPEC N5),
the default Delaunay test (SPEC N6). For the second half, taken from the reference: tracing to
the circumcentre from the barycentre of the triangle, and splitting an edge when the
circumcentre lies on it. Different by design: the boundary rule (encroachment is tested,
shells at small angles, no vertex is removed; SPEC N1), what is left below the bound is
reported instead of skipped (SPEC N4), the overlay is built from the traces with a
consistency check and as a mesh (the reference draws it; SPEC N12). The reference has no vertex
repositioning and no optimal Delaunay triangulation (only `removeInsertedVertex`, and an
extrinsic `smoothByCircumcenter` in geometry-central's `remeshing.cpp`); those are written from
the paper, Chen and Holst 2011 and SPEC N19 to N26.

## Reproducing the tables

```sh
export CARGO_TARGET_DIR=$PWD/target-intrinsic
cargo run --release -p ddg-intrinsic --example delaunay_report -- data/models/spot/spot/spot_triangulated.obj ...
cargo test --release -p ddg-intrinsic --test models --test reference --test distance -- --nocapture
cargo test --release -p ddg-intrinsic --test refinement --test export --test overlay -- --nocapture
```

For the optimal Delaunay triangulation (one process per mesh in parallel for the table, then
one after the other with `--no-cond` for the times) and the cut flow (one process per
subcommand and treatment, `TREAT=input|flip|r25|r30|odt|odt15`; `speed` alone):

```sh
cargo run --release -p ddg-intrinsic --example odt_study -- data/models/spot/spot/spot_triangulated.obj sheared165 ...
TREAT=odt cargo run --release -p ddg-intrinsic --example flow_on_bad_meshes -- laws
cargo bench -p ddg-intrinsic --bench retriangulation
```

For the remeshing study (one process per mesh in parallel, one log each; `blub` and the camels
with `TREATS=input,flip,flipS,u2,u2S,ka2,ka2S,odt2`, the camels also with `TARGET=area` and
`REF_FACES=2000000`; `speed` alone; `SKIP_FROM=<log>` resumes an interrupted mesh):

```sh
STARTS=0,1 cargo run --release -p ddg-uv --example remesh_study -- real spot_control
cargo run --release -p ddg-uv --example remesh_study -- laws
cargo run --release -p ddg-uv --example remesh_study -- speed spot
cargo run --release -p ddg-uv --example remesh_study -- repro bob_control ka4 dirichlet 1 1
uv run scripts/remesh_tables.py runs_tmp/remesh/real2/*.log
```

The example prints, per mesh and angle bound, the refinement report, traces, overlay and
export, and a line starting with `refine |` for the table; run it alone when the timings
matter. The tests run one job per mesh in parallel and print the other tables.

The meshes are not part of the repository (`data/` is ignored); tests that need them print
`SKIPPED` and pass when a file is absent.
