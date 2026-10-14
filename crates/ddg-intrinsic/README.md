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
every decision the paper leaves open (entries `N1` to `N17`).

## What is implemented

- The data structure: a halfedge mesh that allows what flips produce (two edges between the
  same vertices, an edge from a vertex to itself, a triangle that uses a vertex twice,
  vertices of degree one), a length per edge, a signpost per halfedge (normalized to `2π` at
  interior vertices and to `π` at boundary vertices), the angle sum of every vertex, and the
  place of every vertex on the input mesh (paper Sec. 3.1).
- Building it from a `ddg_core::TriMesh` with a `Geometry` (Alg. 4).
- The edge flip (Alg. 5) with every refusal: boundary edge, end point of degree one, not
  strictly convex, too thin (SPEC N5).
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
  the bound is listed with its reason.
- The simplicial export for code that takes a `TriMesh` with a `Geometry`, the way back for a
  cut found there, and a fingerprint that identifies the state (SPEC §6, N15).
- The cotan Laplacian (positive semidefinite) and the lumped vertex areas of the triangulation,
  angle defects, and quality numbers.

Not implemented: vertex repositioning (Alg. 10) and vertex removal; a refinement that keeps
the angle bound through the export at spikes (SPEC N15: on thin tubes and handles the export
can restore it, at vertices with a small angle sum it cannot); benches (timings are those of
the example program). The Lean theorems for the formula ids below are not written:
`scripts/check-formulas.sh` reports the twelve `it-` ids as missing (the four ids shared with
`ddg-core` resolve to their proofs).

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
| `flip_edge(e) -> Result<Flip, IntrinsicError>`, `can_flip(e)` | one flip, or the reason it is refused (`FlipRefusal`) |
| `flip_to_delaunay(&DelaunayOptions) -> Result<DelaunayReport, _>`, `is_delaunay`, `edge_angle_excess`, `edge_cotan_weight` | Delaunay flipping; `DelaunayCriterion::{AngleSum, CotanWeight, REFERENCE}`, `FlipLimit` |
| `InputSurface::trace`, `IntrinsicTriangulation::trace` `(start, direction, length, &TraceOptions) -> Result<Trace, _>` | a straight path over the input mesh / over the triangulation from a `SurfacePoint` |
| `trace_edge(e, &EdgeTraceOptions)`, `trace_edges(&EdgeTraceOptions) -> Result<EdgeTraces, _>` | the checked path of an edge / of all edges over the input mesh |
| `cotan_laplacian`, `lumped_mass`, `lumped_mass_normalized`, `corner_angle_sums`, `angle_defects`, `boundary_geodesic_curvatures`, `quality` | operators and statistics |
| `insert_vertex(f, bary)`, `split_edge(e, t) -> Result<Insertion, _>` | a new vertex in a face / on an edge, placed on the input mesh; or the reason it is refused (`InsertRefusal`) |
| `point_on_input(point)`, `point_on_triangulation(point)` | a point of one triangulation as a point of the other |
| `segments_on_input(&nodes, &segments, &EdgeTraceOptions) -> Result<PathsOnInput, _>`, `segment_on_input(f, from, to, _)` | straight segments inside faces of the triangulation as pieces on input faces, checked like edge traces |
| `overlay(&EdgeTraces) -> Result<Overlay, _>` | the common subdivision: `vertices` (on both triangulations), polygons `face(f)` with `face_input`, `face_triangulation`; `triangles()` |
| `refine(&RefineOptions) -> Result<RefineReport, _>` | Delaunay refinement; `min_angle` (default 25°), `max_circumradius`, `InsertionLimit`; the report lists `below_bound` with `BelowReason` |
| `export_simplicial(&ExportOptions) -> Result<SimplicialExport, _>` | splits what is not simplicial; `mesh`, `geometry`, `locations`, `n_splits`, `min_angle`, `max_angle`, `is_delaunay`, `refinement`, `fingerprint`, `edge_point`, `triangulation_edge`, `check_state` |
| `fingerprint()` | a hash of connectivity, lengths and places |
| `irregularities`, `to_tri_mesh`, `check` | what is not simplicial; the triangulation as a `TriMesh` if it is; structure and identities |
| `formulas::*` | the closed-form formulas as plain functions |
| `IntrinsicError`, `FlipRefusal`, `EdgeTraceFailure` | every error, with the element concerned |

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
- **Refinement is deterministic** and so is the export: the same input and options give the
  same triangulation bit for bit, hence the same fingerprint.
- **`export_simplicial` changes the triangulation** (splits, and flips unless turned off); the
  way back must be asked of the triangulation in that state. **Splitting does not keep an
  angle bound**: read `min_angle`, `max_angle` and `refinement` of the export, not the report
  of an earlier refinement.
- **Three numbers decide what is accepted** (all in SPEC §8): an edge is flipped by Delaunay
  flipping when its opposite angles sum to more than `π + 1e-10`; a flip is refused when it
  would create a triangle with `2 A < 1e-6 · ℓ_new · ℓ_max`; an edge trace must arrive within
  `1e-6` of the edge's length of its vertex. The first and the third are options.

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
- `tests/oracle.rs`: the 18 fixtures of the Lean oracle (`lean/oracle/fixtures/intrinsic`),
  an implementation of single operations written from the paper and `SPEC.md` without this
  crate's code: flips with the triangulation before and after, Delaunay flipping (also with
  self-edges and double edges), traces from vertices, the trace of every halfedge, and
  circumcentres with their traces; 23 000 comparisons, largest difference `4e-14` (SPEC N17).
  And the oracle's readings where no fixture exists: a direction equal to a signpost, the
  first halfedge of a vertex, the lengths of an edge split.
- `tests/export.rs`: the export is a simplicial mesh of the same surface with the same ids,
  lengths and places; reproducible bit for bit; a stale state is refused; the angles before
  and after; a cut (a closed curve, and a network with a triple point) comes back as pieces
  on input faces; the end-to-end run above.

Not verified: termination of refinement beyond what SPEC N1 and N4 argue (it is observed);
the export's refinement rounds on meshes other than those in the table; the way back on
curved shapes is checked by the trace's own end-point check and by round trips, not against
positions (the position checks are on flat shapes); cuts whose nodes fall on vertices of the
export (the flow keeps crossings off vertices); the Lean proofs of SPEC §7; the oracle has no
fixtures for insertion, split, point queries, overlay, refinement and export; no mutation check
was run on the second half.

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
consistency check and as a mesh (the reference draws it; SPEC N12).

## Reproducing the tables

```sh
export CARGO_TARGET_DIR=$PWD/target-intrinsic
cargo run --release -p ddg-intrinsic --example delaunay_report -- data/models/spot/spot/spot_triangulated.obj ...
cargo test --release -p ddg-intrinsic --test models --test reference --test distance -- --nocapture
cargo test --release -p ddg-intrinsic --test refinement --test export --test overlay -- --nocapture
```

The example prints, per mesh and angle bound, the refinement report, traces, overlay and
export, and a line starting with `refine |` for the table; run it alone when the timings
matter. The tests run one job per mesh in parallel and print the other tables.

The meshes are not part of the repository (`data/` is ignored); tests that need them print
`SKIPPED` and pass when a file is absent.
