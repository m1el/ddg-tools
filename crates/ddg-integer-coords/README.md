# ddg-integer-coords

Intrinsic triangulations stored by integer coordinates, after

> Mark Gillespie, Nicholas Sharp, Keenan Crane. *Integer Coordinates for Intrinsic Geometry
> Processing.* ACM Transactions on Graphics 40(6), Article 252, 2021.
> <https://doi.org/10.1145/3478513.3480522>

The contract is [`SPEC.md`](SPEC.md) (decisions `I1`–`I27`, work packages in §12). This crate is
the second representation of intrinsic triangulations in the workspace; the first, signposts
(Sharp, Soliman, Crane 2019), is `crates/ddg-intrinsic`, whose `Complex` (the connectivity) it
shares.

## What is here (WP1)

| | |
|---|---|
| `IntegerTriangulation` | connectivity (`ddg_intrinsic::Complex`), `Arc<InputSurface>`, lengths `f64`, normal coordinates `n: i64` per edge (SPEC I2), roundabouts `u32` per halfedge (`NO_ROUNDABOUT` for tails in `V*`), places `q⁰` of inserted vertices. `new(mesh, geometry)` is `T¹ = T⁰`; `with_options(mesh, lengths, &InitOptions)` mollifies first if asked (`ddg_core::mollify`, off by default, `InitReport::mollified`). `from_parts` / `into_parts` for stored data. Read-only accessors. |
| `LocalIndices` | the local indices of Supp. B with the start of SPEC I3 (`TriMesh::outgoing_halfedge`; at a boundary vertex the last index is the arriving boundary edge, no wrap). |
| `formulas` | Eq. 1 (`corner_coordinate`, never a floor division, I14), Eq. 2 (`emanating`), face validity, `reverse_crossing`, Alg. 1's step (`trace_step`), Alg. 3's index (`identify_index`), Eq. 6 (`split_edge_normal`), `FaceCoordinates`; all equal to the Lean oracle's. |
| `trace` | combinatorial Supp. Alg. 1–4: `trace_from` / `trace_with` (Alg. 1), `extract_curve` (Alg. 2 without geometry), `identify_traced_edge` (Alg. 3, walking back past inserted vertices on the edge, halfedge `a→i`, SPEC §3.1), `extract_edge` (Alg. 4 with the continuation through inserted vertices, I17), `input_directions` (the window of Alg. 4 at a vertex, all directions at once), `curve_ends`. Generic over `NormalGluing` so that hand-built Δ-complexes can be traced. |
| `certify()` | SPEC §2.4, all four items (below). |
| `from_intrinsic` | integer coordinates read off a signpost triangulation by tracing its edges (below). |

**Not yet** (later work packages): refinement (WP6), transfer (WP7).

## Operations (WP3)

| | |
|---|---|
| `flip_edge`, `can_flip` | refusals and new length from `ddg_intrinsic::plan_flip` (bitwise as signposts); Eq. 3 with the corrected first terms `c^{ji}_l + c^{ij}_k` in its doubled form (`i128`, an odd value is an error, I14; outside `i64` → `Overflow { edge }`; `n_kl < −1` → `InvalidCoordinates` naming the edge); Eq. 4 for both new halfedges with tails in `V⁰` (`LocalIndices::modulus`, no reduction at the boundary). All from the pre-flip state, then `Complex::flip_edge` + `move_halfedge_data`. |
| Delaunay flipping | the shared `Triangulation::flip_to_delaunay`; `flip_to_delaunay_certified` certifies after it (I4). |
| `split_edge` | refusals and lengths from `plan_edge_split`. `n = −1` (all boundary edges): Eq. 6 on both sides, `q⁰ = (1 − t) q⁰_i + t q⁰_j` on the input edge (the signposts' arithmetic). `n ≥ 0` (I20, explicit, also for constrained edges): crossings before the point from `extract_intrinsic_edge` (`distance < t ℓ`), `q⁰` interpolated in the input face between the two neighbouring crossings (refused within rounding of a crossing), the integers by `formulas::split_side` (the arcs separated by the new chord; Eq. 6 is its `n = −1` case). `split_edge_given` takes both decisions from the caller. |
| `remove_vertex` | inserted interior vertices: flips by I24 (degree must drop; loops first, then `min(A₁, A₂)/(A₁ + A₂)`, ties by edge id; bound `10 deg`), then `Complex::remove_vertex` with `VertexRemoval` moving `n`, `r`, lengths, `q⁰`. Input vertex / constrained spoke: refused before anything changes; boundary: `Unsupported` (I7); no flip left: `RemovalStuck`, with the state restored exactly (connectivity from a copy, the rest from a journal). Returns `Removal` (the flips, the renumbering). |
| `impl Triangulation` | lengths, angle sums, places, flips, splits, the constrained mask; `trace_edge_on_input` = `extract_intrinsic_edge` (WP2); errors through `to_intrinsic_error` (`IntrinsicError::Representation`, `"integer coordinates"`). |

**I27.** A vertex made by an edge split has two spokes at an angle of `π`: their flips would make
a triangle without area (`plan_flip` refuses them), and of the two halves the one with the angle
below `π` at the vertex can be flipped only if its quadrilateral is convex at the far end. When
it cannot, the removal is `RemovalStuck` (state restored); merging the halves would need an
inverse of `Complex::split_edge`, which `ddg_intrinsic` does not have. Measured below: about a
third of fresh split vertices, half after Delaunay flipping; face insertions almost never (12 of
1 160).

**Timings** (release, one pinned core, load average 3.3, 2026-10-31, best of 5 / 3;
`cargo test --release --test ops -- --ignored --nocapture timing`):

| mesh | \|E\| | Delaunay flipping (flips) | signposts | integer | certify after | random flips/s, signposts | integer | removal |
|---|---|---|---|---|---|---|---|---|
| spot | 8 784 | 281 | 1.57 ms | 1.63 ms | 1.48 ms | 3.74 M | 3.83 M | 7.8 µs |
| camelbad | 34 116 | 8 110 | 8.87 ms | 9.29 ms | 6.94 ms | 3.27 M | 3.23 M | 25 µs |

Integer flips cost the same as signpost flips (Eq. 3, 4 are a few dozen integer operations; the
time is in `plan_flip`, shared). A removal copies the connectivity once (O(|E|), most of its
cost on camelbad) so that a stuck removal can be undone exactly.

## Face split, App. A, moves (WP4)

| | |
|---|---|
| `split_face(f, bary, &FaceSplitOptions)`, `insert_vertex(f, bary)` | Supp. Alg. 5, 6. Lengths and refusals from `plan_face_split` (I25, bitwise as signposts). The arcs of the face come from the integers (`insert::FaceArrangement`: `c` arcs around each corner, the fan of `e` arcs); their end points from `extract_intrinsic_edge` of the three sides (`distance / ℓ`, clamped into `[0, 1]` and counted). The regions are enumerated in the oracle's labels (`Region::Corner { corner, index }`, `Region::Middle { index }`); floating point only **picks** one: the largest margin of side tests against its bounding chords, evaluated in barycentric form in the frame of each chord's corner (`ab − b v₁ − a v₂`, accurate however close a chord is to its corner); the counts are then `FaceArrangement::split_counts` of that region (I15; Lean `faceSplit?_gapRegion`, valid whichever region). New roundabouts at input-vertex corners: `r + n⁻ + the arcs the corner emanates before m`. `FaceSplitReport` has the `Insertion`, the region choice (`margin`, `runner_up`, `curve_distance`: distance to the nearest arc in units of the shortest side, I6) and the recovery. `too_close > distance` → `TooCloseToCurve`, nothing changed (default 0: never). |
| `locate_in_face`, `point_on_input` | the same without inserting (SPEC §4.9, `T¹ → T⁰`; same signature as the signposts'); edge points between neighbouring crossings, points along input edges interpolated; results within `1e-12` of a side of the input face are returned as points of the input edge (vertex), as WP5's `LOCATE_SNAP`. |
| `recover::recover_barycentric` | App. A: minimum-norm `ξ` of Eq. 15 by QR (MGS with reorthogonalization) over **every** corner of the region (3 to 6, I16), then Eq. 16. The input face of the region is read off a gap of the region on a side (the input face `ExtractIntrinsicEdge` says the side runs through there), or the roundabout when all three sides run along input edges. An insertion place within `1e-12` outside `(0, 1)` is moved towards the region's centroid by the smallest fraction (`RecoverReport::moved`), farther is `RecoverFailed`. |
| `move_vertex(v, displacement, &options)` | Sec. 3.6: `AngleFrames::of(&t).trace` from `v` (I9; frame 0 = `outgoing_halfedge`), face split (or edge split if the path ends on an edge), `remove_vertex(v)`; the moved vertex keeps id `v`. Any refusal (incl. `RemovalStuck`, I27) restores the whole state (the reference ignores a failed removal). Boundary vertices: `Unsupported` (I7). |

Results (`tests/insert.rs`, release; the 17 shapes, spot, camelbad, 7 scenarios, 2 seeds; points
inside faces, near sides down to `1e-7`, near corners, in the most crowded faces (up to 70
crossings), and at `1e-3 … 1e-16` from the middle of an arc; 834 insertions in faces with a
repeated vertex):

* **26 099 insertions in both representations: the integer state equals `from_intrinsic` of the
  signposts exactly in all but 9**, all 9 within `1.2e-10` shortest sides of an arc (7 below
  `3.2e-14`): there the two representations put the point on different sides of the arc by
  rounding, and both states are valid. 816 further insertions at `≤ 1e-12` from an arc are snapped
  by the signposts onto the input edge (a vertex on the input edge, by design different). Refusals
  identical (3 581); neither side refused what the other inserted. Places `q⁰` within `5.6e-10`
  mean edge lengths of the signposts' `vertex_location` (worst on a jittered icosphere face with
  20 crossings). App. A used 3 / 4 / 5 / 6 corners 16 975 / 7 520 / 1 592 / 12 times, always
  exactly the corners of the region; 12 places moved inside (smallest coordinate `−2.2e-13`);
  residual of Eq. 15 `≤ 5.6e-16`.
* Before the corner-frame side tests, one insertion in a face whose arc around a corner had both
  ends within `1e-12` of it (a vertex `4e-12` from an input edge) got the wrong region by `0.04`
  shortest sides from side tests in the plane, and App. A refused it (`RecoverFailed`): the
  integers stay valid whichever region, the place shows the error.
* App. A on flat meshes (the place is known): within `8.6e-14` mean edge lengths over 2 652
  insertions; on `T¹ = T⁰` the place is the point itself (error 0). Random affine maps, 3–6
  corners: `4e-14`.
* Many insertions (8 691) in 3 batches with Delaunay flipping (54 939 flips) and certification
  per batch, then removal of all inserted vertices: 8 639 removed, 52 `RemovalStuck`.
* Point queries: 10 640 agree with the signposts' `point_on_input` within `6e-12` mean edge
  lengths; a vertex inserted at a point is where the query says (`1e-14` barycentric).
* I6: the reported distance refuses with a limit above it (561, state unchanged) and not below.
* Moves: 1 603 inside the star against the signposts' local `move_vertex`: new places within
  `7.6e-13` mean edge lengths, after Delaunay flipping both the lengths agree within `1.4e-13`
  and the normal coordinates as multisets agree (2 more differ where Delaunay ties were broken
  differently); 951 moves beyond the star (the signposts cannot), all certify; 71 refused as
  `RemovalStuck`, 606 because the path leaves the surface, every refused move restores the state.
* Oracle (`tests/oracle.rs`): `regions` (514 faces, 2 469 regions) and `split` (516 faces, 2 457
  region counts) exactly.

**Timings** (release, one pinned core, load average 4.7 from other jobs, best of 5, 2 000 random
points; `cargo test --release --test insert -- --ignored --nocapture timing`):

| mesh | state | crossings per face | signposts' insertion | integer insertion | point query |
|---|---|---|---|---|---|
| spot | Delaunay | 0.1 | 0.90 µs | 1.17 µs | 0.72 µs |
| spot | 5 \|E\| random flips | 3.0 | 0.96 µs | 2.63 µs | 2.15 µs |
| camelbad | Delaunay | 1.1 | 1.67 µs | 2.64 µs | 1.59 µs |
| camelbad | 5 \|E\| random flips | 2.8 | 1.56 µs | 3.39 µs | 2.47 µs |

The integer insertion pays for `ExtractIntrinsicEdge` of the three sides (walks along the curves
through them); the signposts trace one path over the input mesh.

## Common subdivision, P0, P1, point queries (WP5)

* `IntegerTriangulation::common_subdivision()` -> `CommonSubdivision { overlay, fingerprint, report }`:
  the signpost crate's `Overlay` type (same vertex order: `V^1`, then the crossings edge by edge from
  the tail of the stored halfedge), built per face of `T^1` from the integers (Supp. Alg. 9; `face_arcs`
  is the per-face pairing, tested exactly against the oracle's `subdivide`) and the crossings of
  `extract_all_edges` + Alg. 2. `fingerprint()` is computed, never cached (SPEC I19).
* `interpolation_matrices(&subdivision)` -> `(P0, P1)` as `ddg_core::SparseMatrix`; refuses a
  subdivision of another state.
* `point_on_triangulation(point)` / `point_locator()`: `T^0 -> T^1`, pieces of an input face found
  lazily by a search over the faces of `T^1` (see `subdivision.rs`, deviation from SPEC I10).
* `extract_all_edges()`: `extract_edge` of every input edge with `input_directions` once per vertex
  (2.5-3x faster than one by one on spot and camelbad).
* Results (tests/subdivision.rs, 266 states): same connectivity as labelled polygons as the
  signposts' overlay (0 mismatches over 863 301 pieces, 272 651 crossings); crossing positions
  `|d| sin(theta) <= 3.5e-10` (raw up to 6.5e-8 at grazing crossings); the signpost overlay checks
  pass; `P0 x0` reproduces the positions to 5e-15 edge lengths. Timings (release, one thread, load
  about 6): spot after Delaunay 4.0 ms (signposts trace + overlay 1.15), camelbad 24 ms (8.7).

## Attribute transfer, Fig. 14 (WP7)

`Transfer::new(&tri, &common_subdivision)` (or `Transfer::from_state(&tri)`) builds `M_S`, `P₀`,
`P₁`, `M₀ = P₀ᵀ M_S P₀`, `M₁ = P₁ᵀ M_S P₁` and their Cholesky factorizations once; it refuses a
stale subdivision (fingerprint, SPEC I19). `to_input(&f)` solves Eq. 14 (`T¹ → T⁰`),
`to_intrinsic(&g)` the reverse (`M₁ ĝ = P₁ᵀ M_S P₀ g`), the `_report` variants return the relative
residual of the normal equations and the `L²` distance (Eq. 13); `copy_to_input` /
`copy_to_intrinsic` are plain copying. `M_S` is the **unnormalized** consistent P1 mass matrix of
the pieces (`Σ M_S` = total area; the area normalization of `CLAUDE.md` belongs to the cones
problem), with piece areas from barycentric coordinates in the input face and the face areas, no
positions needed. The crate does not depend on `ddg-intrinsic-apps`.

Measured (`tests/transfer.rs`, 133 states: the standard shapes, spot and camelbad, 7 scenarios):

* `P₀ᵀ M_S P₀` against the input mesh's own `A/12 (1 + δ)` matrix (`ic-transfer-mass-identity`):
  worst entry error `2.1e-15`, relative to the larger diagonal entry of the pair; total area of
  the pieces within `1.1e-14`;
* constants reproduced to `5.3e-15` (both directions), affine functions on the flat shapes to
  rounding; normal-equation residual at most `2.6e-16` (relative);
* the transfer is never farther from `P₁ f` (resp. `P₀ g`) than copying; on camelbad after
  Delaunay flipping `2.3e-4` against `3.5e-4`;
* stale subdivisions, wrong lengths and NaN/infinite values are errors.

`examples/fig14.rs` (`cargo run --release -p ddg-integer-coords --example fig14`) refines by
inserting circumcenters of the largest faces with Delaunay flipping (WP6's refinement was not
committed yet), samples a smooth function at the vertices of `T¹` and prints relative `L²` errors
against it on the input surface (PLY with colours in `runs_tmp/wp7-fig14/`):

| mesh | extra vertices | copy | `L²` transfer | `T¹` itself | copy / transfer |
|---|---|---|---|---|---|
| spot | 0 | 1.069e-2 | 1.074e-2 | 1.081e-2 | 0.99 |
| spot | 1465 | 1.069e-2 | 7.248e-3 | 6.572e-3 | 1.47 |
| spot | 2930 | 1.069e-2 | 6.073e-3 | 4.582e-3 | 1.76 |
| camelbad | 0 | 5.275e-3 | 5.020e-3 | 4.946e-3 | 1.05 |
| camelbad | 5696 | 5.275e-3 | 3.071e-3 | 2.158e-3 | 1.72 |
| camelbad | 11392 | 5.275e-3 | 2.811e-3 | 1.603e-3 | 1.88 |

The transfer beats copying and the gap grows with the refinement, but the factors are far from the
paper's 28x: the paper's function and refinement size are not reproducible from the text, and our
sampled function is smooth and low-frequency. The transfer cannot beat the refined triangulation's
own function here (it projects it to the coarser space).

Timings (release, single thread, pinned core, load average about 3): spot (final level, 5 860
vertices, `|S|` = 12 720): common subdivision 8.6 ms, `M_S` 7.9 ms, `P₀ P₁` 1.6 ms, `M₀ M₁` 19.5 ms,
factorize both 11.1 ms, one transfer 0.5 ms. Camelbad (22 784 vertices, `|S|` = 62 825): 38.8,
42.5, 7.5, 113, 46.9 ms, one transfer 2.8 ms.

## Geometry of the crossings (WP2)

Connectivity comes from the integers; floating point decides only where the crossings lie. Both
directions unfold a strip of triangles (`it-unfold-apex`) and cut one straight segment with its
interior edges (`it-ray-edge-intersection`); a parameter outside `[0, 1]` by rounding is clamped
and counted (`StripReport::clamped_crossings`, SPEC I5).

| | |
|---|---|
| `extract_geometric_curve(ζ)` | Supp. Alg. 2: the segment through a crossing, its input edge (Alg. 3), `u_start`/`u_end` from the places of its ends (I17), `GeometricCrossing { crossing, input_edge, u, v, sin_angle }` (`u` in the input edge's direction, `v` along the crossed halfedge, `sin_angle` the conditioning). `extract_geometric_crossing(ζ)`: one crossing of it. `segment_geometry`: the strip of a given `CurveSegment`. |
| `extract_edge_geometry(h⁰)` | Alg. 4 with positions: WP1's `InputEdgePath`, the parameter of every vertex on it, the crossings of every piece. |
| `extract_intrinsic_edge(e, &EdgeTraceOptions)` | `ExtractIntrinsicEdge` (FlipOut R2), returned as the signposts' `EdgeTrace`: each crossing's input edge identified from the integers (Alg. 1 + 3, memoized along the walks), the input faces followed exactly from the tail's roundabout or `q⁰` to the head's (consecutive input edges must bound a common face), positions from the straight segment through those input faces. `mismatch` and `clamp` are reported, never a failure. `_with_report` adds the `StripReport`. It is `trace_edge_on_input` of the shared interface. |

Tests (`tests/geometry.rs`, the 266 states of `tests/states.rs`, 37 of them with self-edges):
`extract_intrinsic_edge` of all 678 806 edges crosses exactly the input halfedges the signposts'
`trace_edge` crosses, in the same order (272 651 crossings, first faces and `along` equal); the
two directions of every crossing (Alg. 2 over `T¹`, `extract_intrinsic_edge` over `T⁰`) and the
signposts' overlay agree to `|Δ| sin θ ≤ 3.5e-10`; raw `|Δt|` exceeds `1e-9` at 2 grazing
crossings only (in 3 states; sin θ = 1.4e-4 and 9.3e-4, `|Δt|` up to 5.1e-8); 0 clamps; mutated integers give
results or errors, never panics.

Cost (`examples/extract_cost.rs`, after Delaunay flipping, release, one pinned core, load
average 10.7 from other jobs, best of five): all edges by `extract_intrinsic_edge` spot 0.33 ms
(37 ns per edge + crossing), camelbad 5.3 ms (112 ns); the signposts' `trace_edges` 0.20 / 3.3 ms.
All input edges by `extract_edge_geometry` 7.0 / 31 ms, of which WP1's combinatorial
`extract_edge` is 6.4 / 27 ms (`input_directions` at both ends of every edge).

## The certificate

`certify()` returns a `CertifyReport` or `IntegerError::InvalidCoordinates` with the item of
SPEC §2.4 (`ValidityItem`), the place (face, edge, halfedge, vertex, crossing or input edge) and
the values. In order:

0. ranges: array sizes, `n ≥ −1`, `n = −1` on boundary edges, roundabouts in `[0, deg₀)` exactly
   on halfedges leaving input vertices, `q⁰` of the right kind, boundary status of input vertices;
1. item 1: even numerators of Eq. 1 in every face (item 2 follows from Eq. 2);
2. item 3: every edge with `n = −1` joins two distinct vertices along one input edge, named by the
   roundabouts at its input-vertex ends and by the places of its inserted ends;
3. item 4: at every input vertex the roundabouts agree with the emanating counts (going ccw, the
   next roundabout is `r + n⁻ + e`, the counts sum to `deg₀`); every input edge followed from its
   tail arrives at its head by the end the roundabouts there give, through inserted vertices that
   lie on it in order (two curve ends each); every crossing and every `n = −1` edge is used by
   exactly one input edge (a bitmap of `Σ n⁺` bits).

**Cost (SPEC I4).** Linear: every crossing is visited once. `examples/certify_cost.rs`, release,
one pinned core, load average 0.6 (2026-10-31), best of five, after `k·|E|` random flips:

| mesh | flips | \|E\| | Σn⁺ | certify | ns / (\|E\| + Σn⁺) | ≈ flips | signposts' `trace_edges` |
|---|---|---|---|---|---|---|---|
| icosphere 4 | 0 | 7 680 | 0 | 1.2 ms | 159 | – | 0.1 ms |
| icosphere 4 | 50 E | 7 680 | 11 704 | 2.0 ms | 101 | 13 400 | 2.6 ms |
| torus 48×16 | 50 E | 2 304 | 4 708 | 0.59 ms | 84 | 4 300 | 0.8 ms |
| bumpy grid 40×40 | 50 E | 4 880 | 6 477 | 1.1 ms | 101 | 8 700 | 1.5 ms |
| spot | 10 E | 8 784 | 9 876 | 2.1 ms | 114 | 14 600 | 2.5 ms |
| spot | 50 E | 8 784 | 13 353 | 2.3 ms | 103 | 15 800 | 3.0 ms |
| camelbad | 50 E | 34 116 | 55 244 | 9.8 ms | 110 | 60 000 | 13.1 ms |

"≈ flips" is the certificate's time over the time of one signpost flip on the same mesh (0.13 to
0.16 µs): one certificate costs 1.7 to 1.9 `|E|` flips, a little less than tracing every edge with
signposts. Recommended default: certify at the end
of every public operation that does a variable amount of work (Delaunay flipping, refinement,
before the common subdivision, a FlipOut path), never per flip or split; in debug builds and tests,
after every operation on small meshes and per batch on large ones.

## Conversion from signposts

`IntegerTriangulation::from_intrinsic(&IntrinsicTriangulation, &EdgeTraceOptions)` traces every
edge over the input mesh with the signposts: `n` is the number of input edges crossed (`−1` along
an input edge), the roundabout of an end at an input vertex is the local index after the input
corner where the path leaves it, `q⁰` is the signposts' `vertex_location`. Floating point decides
this, so the result is certified before it is returned; a trace that fails or integers that do not
certify are an error. Two exact readings are applied before: a path without crossings between two
points of one side of its input face runs along that side, and a last crossing with the input edge
on which the far vertex lies is the arrival. On the 17 standard shapes of `ddg-intrinsic`'s tests,
spot and camelbad, in 7 scenarios (initial, random flips, Delaunay, insertions, edge splits,
splits of input edges then flips, mixed) and 2 seeds, **266 of 266 conversions certified**; neither
exact reading was needed.

## Tests

* `tests/states.rs`: the 266 conversions; every crossing of every edge (272 651), extracted in both
  directions and identified by Alg. 3, names the input edge the signposts' trace crosses there
  (1 328 of them by walking back past inserted vertices); `extract_edge` of every input halfedge
  arrives where it should and the paths use `Σ n⁺` crossings; **15 941 single mutations** (`n ± 1`,
  a roundabout `± 1`, two roundabouts at a vertex swapped) all rejected — by item 4 (extraction):
  6 167, item 3: 7 414, ranges: 2 240, item 1: 120. Run with `--nocapture` for the tables.
* `tests/oracle.rs`: the Lean oracle's committed fixtures (`lean/oracle/fixtures/ic`), exact
  equality: `corner` (519 faces), `trace` (192 steps), `tracecurve` (14 curves on the Fig. 5 complex
  and two hand-built Δ-complexes, 2 closed curves), `identify` (168 indices; 34 arguments of the
  random fixture have an invalid face `ija`, where the oracle rounds Eq. 1 and this crate refuses,
  I14), `splitedge` (52), `flip` (305 quads, 302 roundabout updates, 304 flips back).
* `tests/ops.rs` (WP3), against the signposts on the 17 shapes, spot and camelbad, 7 scenarios, 2
  seeds; after every batch the integer state **equals** `from_intrinsic` of the flipped signposts
  (connectivity, lengths bitwise, every `n`, roundabout and place) and certifies:
  * flips: 266 states, 1 008 batches, 912 744 attempts, 437 595 flips, 475 149 refusals (all the
    same in both), then Delaunay flipping with the identical flip sequence and report (183 002
    flips); Δ-complexes up to 23 self-edges and 46 faces with a repeated vertex in one state
    (I21); 0 mismatches;
  * flip then flip back: 10 490 edges, integers restored exactly;
  * edge splits: 4 256 along input edges known to the signposts (exact, places bitwise), 1 060
    along input edges the signposts trace (all but the place exact, places within 2e-11), 5 030
    across input edges with the signposts' crossings and place (exact) and the same 5 030 with this
    crate's own crossings (all but the place exact, places within 1.5e-10 barycentric); refusals
    identical;
  * removals of inserted vertices, the same flips and removal replayed in the signposts, compared
    exactly: 3 884 removed (9 416 flips, 2.4 per removal, at most 14), 1 301 stuck (1 287 of them at
    a straight star from an edge split), states restored exactly; refused removals (input vertex,
    boundary, constrained spoke) change nothing;
  * I27: 2 588 vertices from splits of interior edges (913 across input edges) removed at once:
    1 728 removed, 860 stuck; after Delaunay flipping 1 692 removed, 1 898 stuck (1 873 straight);
  * overflow: hand-built states near `i64::MAX` give `Overflow`, never wrap; Eq. 3 below `−1`
    (quad `2,0,0,0,0`) is refused as invalid.

## The reference implementation

geometry-central (MIT, © 2017–2019 Nicholas Sharp and the geometry-central contributors), commit
`ab116f1`, pinned by the authors' demo, was read for the details the paper leaves out (map:
`runs_tmp/ic-reference-map.md`); nothing was translated from it. Its notice is in
`LICENSES/geometry-central-MIT.txt`.
