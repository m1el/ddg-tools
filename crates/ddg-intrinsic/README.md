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
every decision the paper leaves open (entries `N1` to `N9`).

## What is implemented

This is the first half of the crate.

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
- The cotan Laplacian (positive semidefinite) and the lumped vertex areas of the triangulation,
  angle defects, and quality numbers.

Not yet: vertex insertion and edge split (Alg. 7, 8), vertex repositioning (Alg. 10), point
queries (Alg. 11), the overlay (Sec. 3.4.2), Delaunay refinement (Sec. 4.2), the export of a
triangulation that is not simplicial, and what the cut flow needs on top (SPEC §6). The Lean
theorems for the formula ids below are not written either: `scripts/check-formulas.sh` reports
the ten `it-` ids as missing.

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
(about 5.8 million edge tests per second); tracing runs at 1.5 to 4 million crossings per
second. After flipping, camel_noise1 has 1 self-edge and 5 double edges, camel_noise2 4 and 75,
camelshotnoise 22 and 24; the others are simplicial complexes.

### Edge traces

After Delaunay flipping, with the default tolerance of `1e-6` of the edge length, **no trace
fails its check on any mesh**: 0 of the 281 to 7331 flipped edges per mesh (`tests/models.rs`),
0 on the 25 generated shapes of `tests/common`, 0 on the 15 inputs of the reference comparison.
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

Not verified: anything about vertices that are not input vertices (the tracing kernel accepts
start points and targets on edges and in faces, and the tests start queries there, but no
operation creates such vertices yet); the Lean proofs and the oracle of SPEC §7.

## Formula ids

| id | function | formula |
|---|---|---|
| `it-angle-normalization` | `formulas::angle_scale` | `2π/Θ` at an interior vertex, `π/Θ` at a boundary vertex |
| `it-signpost-update` | `formulas::signpost_update` | `φ_ik = φ_ij + s θ_i^{jk}`, reduced to `[0, 2π)` at an interior vertex |
| `it-unfold-diamond` | `formulas::unfold_diamond` | the two triangles on an edge in the plane: `x = (a² + (b − c)(b + c)) / 2a`, `h = 2A/a` |
| `it-flip-length` | `Diamond::flipped_length` | `ℓ_kl² = (x_k − x_l)² + (h_k + h_l)²` |
| `it-flip-convexity` | `Diamond::flipped_areas` | `x_k h_l + x_l h_k > 0` and `(a − x_k) h_l + (a − x_l) h_k > 0` |
| `it-delaunay-condition` | `formulas::delaunay_angle_excess` | `θ_k + θ_l − π ≤ 0` |
| `it-cotan-weight` | `formulas::opposite_cotan_sum` | `cot θ_k + cot θ_l`, with the sign of `π − (θ_k + θ_l)` |
| `it-unfold-apex` | `formulas::unfold_apex` | `c = a + α d + β d⊥`, `α = (ℓ_ab² + ℓ_ca² − ℓ_bc²)/2ℓ_ab²`, `β = 2A/ℓ_ab²` |
| `it-ray-edge-intersection` | `formulas::axis_crossing` | `t = y_tail / (y_tail − y_head)`, `x = x_tail + t (x_head − x_tail)` |
| `it-barycentric-coordinates` | `formulas::barycentric` | ratios of signed areas |

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
the default Delaunay test (SPEC N6).

## Reproducing the tables

```sh
export CARGO_TARGET_DIR=$PWD/target-intrinsic
cargo run --release -p ddg-intrinsic --example delaunay_report -- data/models/spot/spot/spot_triangulated.obj ...
cargo test --release -p ddg-intrinsic --test models --test reference --test distance -- --nocapture
```

The meshes are not part of the repository (`data/` is ignored); tests that need them print
`SKIPPED` and pass when a file is absent.
