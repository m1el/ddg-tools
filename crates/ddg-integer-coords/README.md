# ddg-integer-coords

Intrinsic triangulations stored by integer coordinates, after

> Mark Gillespie, Nicholas Sharp, Keenan Crane. *Integer Coordinates for Intrinsic Geometry
> Processing.* ACM Transactions on Graphics 40(6), Article 252, 2021.
> <https://doi.org/10.1145/3478513.3480522>

The contract is [`SPEC.md`](SPEC.md) (decisions `I1`–`I30`, work packages in §12). This crate is
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

Delaunay refinement (WP6) and attribute transfer (WP7) are described below.

## Operations (WP3)

| | |
|---|---|
| `flip_edge`, `can_flip` | refusals and new length from `ddg_intrinsic::plan_flip` (bitwise as signposts); Eq. 3 with the corrected first terms `c^{ji}_l + c^{ij}_k` in its doubled form (`i128`, an odd value is an error, I14; outside `i64` → `Overflow { edge }`; `n_kl < −1` → `InvalidCoordinates` naming the edge); Eq. 4 for both new halfedges with tails in `V⁰` (`LocalIndices::modulus`, no reduction at the boundary). All from the pre-flip state, then `Complex::flip_edge` + `move_halfedge_data`. |
| `flip_edge_with(e, FlipRequest)`, `can_flip_with(e, FlipShapeTest)` | the same with the caller's shape test (`plan_flip_with`; FlipOut uses the reference's relative one) and, when `corner_angles` is requested, `Flip::new_corner_angles` from `FlipPlan::new_corner_angles` (two `corner_angle`), bitwise the signposts'; `flip_edge` and every operation of this crate (Delaunay flipping, refinement, removal, insertion) use N5 and request no angles, so they pay nothing for them (IT-SPEC N52) |
| Delaunay flipping | the shared `Triangulation::flip_to_delaunay`; `flip_to_delaunay_certified` certifies after it (I4). |
| `split_edge` | refusals and lengths from `plan_edge_split`. `n = −1` (all boundary edges): Eq. 6 on both sides, `q⁰ = (1 − t) q⁰_i + t q⁰_j` on the input edge (the signposts' arithmetic). `n ≥ 0` (I20, explicit, also for constrained edges): crossings before the point from `extract_intrinsic_edge` (`distance < t ℓ`), `q⁰` interpolated in the input face between the two neighbouring crossings (refused within rounding of a crossing), the integers by `formulas::split_side` (the arcs separated by the new chord; Eq. 6 is its `n = −1` case). `split_edge_given` takes both decisions from the caller. |
| `remove_vertex` | inserted vertices: flips by I24 (degree must drop; loops first, then `min(A₁, A₂)/(A₁ + A₂)`, ties by edge id; bound `10 deg`), then `Complex::remove_vertex` with `VertexRemoval` moving `n`, `r`, lengths, `q⁰`. Where no flip is left at a straight star, and for an inserted boundary vertex: the halves are merged (I27, `merge_edge` below). Input vertex / constrained spoke: refused before anything changes; no flip and no merge: `RemovalStuck`, with the state restored exactly (from a journal of the flips). Returns `Removal` (the flips, the last step `RemovalStep::{Star, Merge}` with its renumbering). |
| `merge_edge(m, e)` | the inverse of `split_edge` (I27, SPEC §4.6): `plan_edge_merge` (degree 4, or 3 on the boundary; the halves straight within `1e-9` rad; `ℓ₁ + ℓ₂`), `Complex::merge_edge` (`EdgeMerge`: smaller ids kept, the last elements fill the freed ones); `n` of the merged edge by `ic-merge-edge-normal` (`−1` along an input edge, `n₁ + n₂ + 1` through a vertex on an input edge, `n₁ + n₂` through a vertex in a face), roundabouts unchanged; refused where the integers see a bend (`CurvesBend`). Split then merge restores the state, the merged length within an ulp. |
| `impl Triangulation` | lengths, angle sums, places, flips, splits, the constrained mask; `trace_edge_on_input` = `extract_intrinsic_edge` (WP2); errors through `to_intrinsic_error` (`IntrinsicError::Representation`, `"integer coordinates"`). |

**I27 (resolved).** A vertex made by an edge split has two spokes at an angle of `π`: their flips
would make a triangle without area (`plan_flip` refuses them), and of the two halves the one with
the angle below `π` at the vertex can be flipped only if its quadrilateral is convex at the far
end. When it cannot, Alg. 8 is stuck; the removal then merges the halves back into one edge
(`ddg_intrinsic::Complex::merge_edge`, the inverse of `split_edge`, after flipping the other edges
down to degree four), in both representations alike. Before: about a third of fresh split
vertices and half after Delaunay flipping were `RemovalStuck`; now none of the fresh ones, and
after Delaunay flipping only stars with an edge from the vertex to itself (74 of 1 574).

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
  per batch, then removal of all inserted vertices: 8 639 removed, 52 `RemovalStuck` (with the
  merge of I27: 8 640, 51; 6 straight-by-angle stars within rounding of an input edge are refused
  by the integers, `CurvesBend`).
* Point queries: 10 640 agree with the signposts' `point_on_input` within `6e-12` mean edge
  lengths; a vertex inserted at a point is where the query says (`1e-14` barycentric).
* I6: the reported distance refuses with a limit above it (561, state unchanged) and not below.
* Moves: 1 603 inside the star against the signposts' local `move_vertex`: new places within
  `7.6e-13` mean edge lengths, after Delaunay flipping both the lengths agree within `1.4e-13`
  and the normal coordinates as multisets agree (2 more differ where Delaunay ties were broken
  differently); 951 moves beyond the star (the signposts cannot), all certify; 71 refused as
  `RemovalStuck`, 606 because the path leaves the surface, every refused move restores the state.
  With the merge of I27 no move is `RemovalStuck` (1 021 beyond the star); of moved edge-split
  vertices 463 done, 0 stuck (before 453, 10).
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

## Input-aligned tangent frames (I26)

* `IntegerTriangulation::frame_directions()` -> `(Vec<f64>, FrameReport)`: the direction of every
  halfedge at its tail in the frame the signposts use (input mesh at `V^0`, input face or edge at
  `V*`), normalized as `IntrinsicTriangulation::signpost`, from the roundabouts and the straight
  segment of each edge in its unfolded strip of input faces (SPEC §4.12, `ic-input-aligned-frame`);
  one `extract_all_edges` plus one strip layout per edge. What `ddg-intrinsic-apps`' tangent fields
  (`vector::TangentFrames`) read on integers.
* Results (`tests/frames.rs`): within `1.14e-10` rad of the signposts on 266 converted states
  (1 352 388 halfedges), within `3.5e-11` after the same operations on both representations where
  the places agree; a flip leaves every other edge's directions bitwise unchanged; mutated integers
  are errors. 2.7 ms on spot (iDT) to 38 ms on camelbad refined (one core).

## Delaunay refinement (WP6)

| | |
|---|---|
| `refine(&IntegerRefineOptions) -> IntegerRefineReport`, `refine_with(_, extra)` | Paper Sec. 4.2, Supp. Alg. 10 (SPEC §5). The driver is `ddg_intrinsic::refine::refine`, written once over `Triangulation` (IT-SPEC N51), so both representations run the same queue, circumcentres, boundary rule and flips (I11). This crate supplies the operations: flips (Eq. 3, 4), the circumcentre traced from the barycentre over `T¹` with angle frames (I9), face splits (Supp. Alg. 5, 6, App. A), edge splits (Alg. 7, I20), removals (Alg. 8 with the shared I24 choice `ddg_intrinsic::removal_flip_choice`). The certificate runs at the end (I4; `certify: false` turns it off). |
| `IntegerRefineOptions` | `refine: RefineOptions` (25°, no circumradius bound, 10 insertions per vertex and at least 1000, the Delaunay test), `rules` (default `RefineRules::Chew`, the paper's: a boundary edge is split at its midpoint when a circumcentre path leaves through it, then the inserted interior vertices at Dijkstra distance `< ℓ` are removed and the triangulation is flipped to Delaunay; the reference's narrow-cone rule, I23. `RefineRules::Ruppert`: the signpost crate's default rules N1, N4), `face_split` (I6; default never refuses). |
| `IntegerRefineReport` | the shared `RefineReport` (inserted, boundary and edge splits, flips, removals, refused removals, trace failures, `min_angle_valid`, `n_narrow_cone_faces`, `below_bound` with reasons) plus face splits, those within `NEAR_ARC = 1e-9` of an arc and the nearest distance (I6), face splits refused because App. A put the place outside its input face (`n_recover_failed`: the face is left, as a refused insertion), places moved inside (I16), clamped crossings (I5), the certificate. |

Results (`tests/refine.rs`, release; run with `--nocapture` for the tables):

* **Against the signposts, 256 runs**: the 17 generator shapes, spot, camelbad, bob, blub,
  camelbadnonuniform, camel_noise1, bump, hand, a bumpy grid 20×20, two closed double cones and
  two open cones with apex angle sums 20°–59° (narrow) and one at 70°; starts from the input,
  from random flips (Δ-complexes), from random insertions and from random edge splits; both rules;
  25° and 30°. **255 runs make the same operations in the same order and end with the same
  connectivity and bitwise the same lengths**; 242 of them equal `from_intrinsic` of the refined
  signposts exactly, the other 13 differ only around signpost places snapped onto input edges
  (129 snaps, by design, WP4). Places within `4.2e-11` mean edge lengths. 33 face splits were
  within `1e-9` shortest sides of an arc (nearest `9.5e-17`), 3 crossings clamped, 4 places
  moved inside, 0 App. A refusals, 0 trace failures. The one run that parts (sheared
  parallelogram 10, σ = 6, Chew, 30°): a circumcentre lies on an interior edge exactly where an
  input edge crosses it; the integer split refuses a point within rounding of a crossing (I20),
  the signposts split and snap the place onto the input edge.
* **Removals** (Chew): 2 307 vertices removed and 175 refused (`RemovalStuck`, mostly I27
  straight stars) over the differential runs; from starts full of edge-split vertices 450
  removed, 80 refused. Refused removals are counted, never fatal. Since straight stars are
  merged (I27, below): 2 304 removed, 0 refused; from split starts 441 removed, 0 refused; 254
  runs identical, the new divergence (needles and caps 6, Chew, 30°) is I28 again.
* **Angles** (Chew, from the input): every run ends `Done`; every corner below the bound is at a
  vertex whose angle sum is below it or in a face the narrow-cone rule leaves. 25° everywhere else
  (camelbad 25.000° with 4 169 insertions, spot 25.047°, hand 25.269°, bumpy grid 25.084°).
  **Thm. C.1**: the closed shapes whose angle sums are all at least 60° reach 30° (spot 30.031°,
  jittered icosphere 30.013°, the tori, the double cone 70°/300°). Faces left by the narrow-cone
  rule: needle tetrahedron 4, tube 3×6 8, sheared parallelograms 2, double cones 6 / 6 / 16, open
  cones 6 / 8; the smallest angle in the other faces there: e.g. double cone 20°/200° 33.3°, open
  cone 30° 25.55° (the double cone 59°/59° has no other face).
* **Termination**: insertion limits 0, 1, 7, 40 are honoured on every shape (`InsertionLimit`
  when the full run needs more; a run that reaches the limit exactly also says
  `InsertionLimit`, as the signposts); every result certifies. Bad options and constrained edges
  are refused with nothing changed.

**Timings** (release, one thread pinned, load average 3.0, best of 5; `taskset -c N cargo test
--release --test refine -- --ignored --nocapture timing`). Refinement to 25° from the input,
including the first Delaunay flipping:

| mesh | \|V⁰\| | rules | inserted | removed | signposts | integer | integer + certificate |
|---|---|---|---|---|---|---|---|
| spot | 2 930 | Ruppert | 255 | – | 4.7 ms | 4.9 ms | 6.5 ms |
| spot | | Chew | 255 | 0 | 4.5 ms | 4.8 ms | 6.3 ms |
| camelbad | 11 392 | Ruppert | 4 169 | – | 39.0 ms | 48.7 ms | 57.3 ms |
| camelbad | | Chew | 4 169 | 0 | 36.8 ms | 45.9 ms | 54.6 ms |
| bunny | 14 290 | Ruppert | 0 | – | 16.8 ms | 16.4 ms | 30.6 ms |
| bunny | | Chew | 0 | 0 | 18.1 ms | 18.5 ms | 31.9 ms |
| bump (boundary) | 3 329 | Ruppert | 0 | – | 3.5 ms | 3.5 ms | 5.3 ms |
| hand (boundary) | 3 372 | Ruppert | 50 | – | 4.6 ms | 4.5 ms | 6.9 ms |
| hand | | Chew | 61 | 11 | 4.9 ms | 4.8 ms | 7.0 ms |
| face (boundary) | 13 746 | Ruppert | 104 | – | 18.1 ms | 18.3 ms | 28.0 ms |
| face | | Chew | 66 | 2 | 18.7 ms | 17.4 ms | 24.5 ms |

The integers cost 0.93 to 1.25 times the signposts (the paper: about 2.7 times, 4 min against
1.5 min on Thingi10k 719791, SPEC §10.1); the difference is in the face splits, which extract
the three sides (WP4). Without insertions the time is the Delaunay flipping and the driver's
passes over the faces.

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

`examples/fig14.rs` (`cargo run --release -p ddg-integer-coords --example fig14 [trials [generator]]`)
is the paper's Fig. 14: `Δf = g`, `g = sin(πx) sin(πy)` on the unit square, `f = 0` on the boundary,
exact `f = −g / (2π²)`, solved as `L f = −M g` (positive semidefinite cotan `L`, Galerkin `M`, `g`
sampled at the vertices) on random low-quality meshes of the square, 100 trials per size, `2⁶ … 2¹¹`
vertices in half octaves like the figure's axis (the paper marks about `2^7.5`). The generator is
ours (the paper gives none): the diagonal-cut square, then repeated splitting of the edge of a
random point's face that is nearest to it, at the projection of the point (seeds `1000 k + t`).
The IDT (same vertices) is `flip_to_delaunay` of the integer triangulation, with `L` and `M`
assembled in the example from its edge lengths. All errors are absolute `‖f − f̃‖_{L²}` over the
unit square (`‖f‖ = 0.0253`), by quadrature on the pieces of the common subdivision, where `f̃` is
`P₀` of the nodal vector for `original` (the P1 solution on `T⁰`), `IDT - copy values` (the IDT's
nodal values on the same vertices) and `IDT - L²` (`Transfer::to_input`), and `P₁` of it for
`IDT itself` (not in the paper). Means over the 100 trials (the paper plots means), and the ratio
of the means against `original`:

| `\|V\|` | original | IDT copy | IDT `L²` | IDT itself | original / copy | original / `L²` |
|---|---|---|---|---|---|---|
| 64 | 9.36e-3 | 5.18e-3 | 3.89e-3 | 3.13e-3 | 1.81 | 2.41 |
| 128 | 7.91e-3 | 3.28e-3 | 2.26e-3 | 1.62e-3 | 2.41 | 3.50 |
| 181 (`2^7.5`) | 6.93e-3 | 2.54e-3 | 1.70e-3 | 1.15e-3 | 2.73 | 4.09 |
| 256 | 6.57e-3 | 2.02e-3 | 1.32e-3 | 8.13e-4 | 3.26 | 4.99 |
| 512 | 5.90e-3 | 1.34e-3 | 8.36e-4 | 4.13e-4 | 4.40 | 7.06 |
| 1024 | 5.45e-3 | 8.90e-4 | 5.45e-4 | 2.06e-4 | 6.13 | 10.01 |
| 2048 | 5.51e-3 | 6.30e-4 | 3.89e-4 | 1.02e-4 | 8.76 | 14.18 |

The medians per trial ratio are a little lower (181: 2.27 and 3.41; 2048: 7.20 and 12.10); the
normal-equation residual is at most `5e-16`. The qualitative picture of the paper is reproduced:
the error on the original mesh stays almost flat, the IDT converges at rate `~1/|V|`, and the
`L²` transfer beats copying at every size (by 1.3 to 1.6; the paper's marked size shows 2.4). The
numbers are not: at the paper's marked size the IDT is 11.5x and the transfer 28.1x more accurate
than the original, ours 2.7x and 4.1x, and the factors grow with `|V|` (8.8x and 14.2x at 2048) but do not reach the paper's. Our meshes are less bad than the paper's (original error
`≈ 7e-3` against its `3e-2`, whose figure shows long fans from the corners), and the ratio depends on
that: the generator is the free parameter, not the transfer. The transfer is never more accurate
than the IDT's own function, as it must be.

`examples/refined_transfer.rs` is the other experiment (on spot and camelbad): refine by inserting
circumcenters of the largest faces with Delaunay flipping (WP6's refinement was not committed
then), sample a smooth function at the vertices of `T¹`, and print relative `L²` errors against it
on the input surface (PLY with colours in `runs_tmp/wp7-refined/`):

| mesh | extra vertices | copy | `L²` transfer | `T¹` itself | copy / transfer |
|---|---|---|---|---|---|
| spot | 0 | 1.069e-2 | 1.074e-2 | 1.081e-2 | 0.99 |
| spot | 1465 | 1.069e-2 | 7.248e-3 | 6.572e-3 | 1.47 |
| spot | 2930 | 1.069e-2 | 6.073e-3 | 4.582e-3 | 1.76 |
| camelbad | 0 | 5.275e-3 | 5.020e-3 | 4.946e-3 | 1.05 |
| camelbad | 5696 | 5.275e-3 | 3.071e-3 | 2.158e-3 | 1.72 |
| camelbad | 11392 | 5.275e-3 | 2.811e-3 | 1.603e-3 | 1.88 |

The transfer beats copying and the gap grows with the refinement; it cannot beat the refined
triangulation's own function (it projects it to the coarser space).

Timings (release, single thread, pinned core, load average about 3): spot (final level, 5 860
vertices, `|S|` = 12 720): common subdivision 8.6 ms, `M_S` 7.9 ms, `P₀ P₁` 1.6 ms, `M₀ M₁` 19.5 ms,
factorize both 11.1 ms, one transfer 0.5 ms. Camelbad (22 784 vertices, `|S|` = 62 825): 38.8,
42.5, 7.5, 113, 46.9 ms, one transfer 2.8 ms.

## Simplicial export, the way back of segments, the optimal Delaunay iteration (I29)

What `ddg-uv`'s cut flow needs, so that it runs on integer coordinates by default (SPEC §4.10,
§4.11):

| | |
|---|---|
| `export_simplicial(&ExportOptions)` | `ddg_intrinsic::export::export_simplicial` on this representation (`impl Exportable`): the signposts' splits and flips; a midpoint within rounding of a crossing (I28) is split at `½ ± 2⁻¹⁶` instead (`n_off_midpoint_splits`); refinement between rounds by the signposts' rules, certified; places are the stored `q⁰`; the fingerprint hashes `n` and `r` too. |
| `segments_on_input(nodes, segments, &EdgeTraceOptions)` | the signposts' signature, output and errors. The input edges a segment crosses are the arcs on the path between the regions of its ends in the tree of a face's arcs (from the integers); positions by intersecting two lines in the face's layout. Checked: every piece lies in its input face, the path's length equals the segment's within the tolerance (`mismatch`). |
| `optimal_delaunay(&IntegerOdtOptions)` | `ddg_intrinsic::odt::optimal_delaunay` (`impl OdtMoves`) in `AngleFrames` frames, certified. |
| `move_vertex_in_star(v, displacement, min_shape)` | the signposts' local move on integers: the star keeps its triangles; new `n` of the spokes from the chords of the star polygon that separate the new place from each link vertex, roundabouts by `ic-roundabout-update`, `q⁰` by App. A. Never needs a removal, so the straight stars of I27 do not stop it. |
| `move_vertex` | unchanged route (Sec. 3.6), but a refused removal after a face split is undone by removing the inserted vertex again instead of restoring a copy of the whole state. |

Against the signposts (`tests/export.rs`, 27 meshes, flips and refined to 25°): 52 of 54 exports
identical in connectivity, bitwise lengths and counts (places within `3.3e-13` mean edge lengths;
the other 2 split off the midpoint on the torus 5×3); 30 041 of 30 042 segments cross the same
input edges in the same order and all lie within `1e-9` mean edge lengths of the signposts' traces;
the optimal Delaunay iteration makes the same iterations, insertions and halvings on all 27 meshes,
places within `3.3e-13` (the bumpy grid takes another trajectory from a last-bit tie). Cost, in
`ddg-uv` session builds (its README): flips and export 1.1×, with refinement 1.2–1.4×, with the
optimal Delaunay iteration 1.2–1.3× the signposts' (SPEC §4.11: a move that stays in its input
face changes no integer and skips the arcs; 3.4–4.2× before); the way back was not timed
separately.

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
  * removals of inserted vertices, the same flips and removal (or merge) replayed in the
    signposts, compared exactly: 3 884 removed (9 416 flips, 2.4 per removal, at most 14), 1 301
    stuck (1 287 of them at a straight star from an edge split), states restored exactly; refused
    removals (input vertex, constrained spoke) change nothing. With the merge (I27): 4 539 removed
    (595 by a merge), 70 stuck (56 at straight stars, all with an edge from the vertex to itself),
    0 mismatches; inserted boundary vertices are merged and certify;
  * I27: 2 588 vertices from splits of interior edges (913 across input edges) removed at once:
    1 728 removed, 860 stuck; after Delaunay flipping 1 692 removed, 1 898 stuck (1 873 straight).
    With the merge: 2 461 of 2 461 removed at once (657 merged); after Delaunay flipping 1 500
    removed (388 merged), 74 stuck (all with a self-edge);
  * merge (`tests/merge.rs`): 18 672 splits merged at once (along interior input edges, on the
    boundary, across input edges) give back the state exactly, the merged length bitwise in 16 330
    and one ulp off otherwise; 2 340 splits of the iDT, Delaunay flipping, every split vertex
    removed (700 by a merge), each state equal to the signposts' and certified, and flipped to
    Delaunay again the iDT itself in the 15 runs without ties (lengths within `9.2e-13`);
  * overflow: hand-built states near `i64::MAX` give `Overflow`, never wrap; Eq. 3 below `−1`
    (quad `2,0,0,0,0`) is refused as invalid.

## Evaluation on Thingi10k and against the authors' demo (WP10)

Thingi10K: Qingnan Zhou, Alec Jacobson, *Thingi10K: A Dataset of 10,000 3D-Printing Models*,
arXiv:1605.04797, 2016. Each model is its Thingiverse author's, under the licence they chose; the
models named below are credited by author and licence in [`DATASETS.md`](../../DATASETS.md).

Everything here is reproduced by `scripts/ic-eval.sh` (conversion, the per-mesh runs, the summary),
`scripts/ic-demo-build.sh` and `scripts/ic-demo-eval.sh` (the authors' demo, built headless outside
the repository), `scripts/ic-timings.sh` (single-thread timings) and the examples
`thingi_convert.rs`, `thingi_eval.rs`; the tables below are `runs_tmp/wp10/summary.md` and
`timings.md`. One process per mesh and representation, 56 at a time, a 10-minute cap and an
address-space limit per job (4 GB up to 100 000 faces, 12 GB up to 1 M, 24 GB beyond); one log and
one JSON per job under `runs_tmp/wp10/jobs/`.

### Conversion (SPEC §13 Q1)

`thingi_convert` reads the 10 000 raw files (9 956 STL, binary and ASCII; 42 OBJ, 1 OFF, 1 PLY),
merges vertices whose coordinates are exactly equal, fan-triangulates polygons, drops unused
vertices and keeps the meshes that `TriMesh::new` accepts (no repair, no removal of degenerate or
duplicate faces). 39 s on 56 threads. **7 745 kept**, 2 255 excluded:

| excluded | files |
|---|---|
| vertex shared by several fans (non-manifold vertex) | 649 |
| two faces with the same three vertices | 616 |
| face with a repeated vertex after the merge | 482 |
| edge on 3 or more faces | 343 |
| inconsistent orientation (95 of them would be repaired by `orient_faces`, not done) | 160 |
| no faces (empty or points/lines only) | 3 |
| unreadable (truncated ASCII STL: 49911, 77942) | 2 |

Against the dataset's own metadata (`geometry_data.csv`: 7 754 vertex- and edge-manifold, all
oriented): 7 744 of those are kept; the 10 that are not: 5 repeated vertex after the merge
(285440, 298321, 298322, 298323, 376253), 3 without faces (74463, 81313, 286163), 1 duplicate face
(527631), 1 non-manifold vertex (82541). One mesh is kept that the metadata calls non-manifold
(73178). The paper's 7 696 (MeshLab conversion to PLY, rule not published) is 49 fewer than ours; the
difference cannot be attributed without MeshLab's rules.

### Table 1 as reproduced

Success = the common subdivision's connectivity is recovered (integers: certificate, then
`common_subdivision`; signposts: traces of all edges, then `overlay`, IT-SPEC N12's refusals), and
the pieces, fan-triangulated, form a valid `TriMesh` with the Euler characteristic, components and
boundary loops of the input. Both after the intrinsic Delaunay triangulation (IDT) and after
refinement to 25° (no insertion limit, as the paper's benchmark; mollified `1e-5`).

| | paper (7 696) IDT | paper refinement | ours (7 745) IDT | ours refinement |
|---|---|---|---|---|
| explicit overlay | 100 % | – | – | – |
| signposts (paper) | 96.0 % | 69.1 % | | |
| signposts, Ruppert rules (this workspace) | | | 100.00 % (7745) | 99.81 % (7730) |
| signposts, Chew rules (this workspace) | | | 100.00 % (7745) | 99.81 % (7730) |
| integer coordinates | 100 % | 100 % | **100.00 % (7745)** | **100.00 % (7745)** |
| integer coordinates, the authors' Delaunay test | | | 100.00 % (7745) | 100.00 % (7745) |
| the authors' demo, integer backend (7 444 meshes up to 200 000 faces) | | | 100.00 % (7444) | 99.97 % (7442) |
| the authors' demo, signpost backend (same) | | | 87.86 % (6540) | 20.71 % (1542) |

Signposts were run on every mesh, with their boundary if any (the paper's signposts were
evaluated without boundary): on the 7 339 closed meshes the numbers are the same to within a mesh
(integers 100.00 %/100.00 %; in the WP10 run, before I30's fix, signposts Ruppert 100.00 %/99.77 %, demo
signposts 87.76 %/20.68 %). The
demo's signpost backend does not mollify its input (the integer backend does, `1e-5`); ours use
the same mollified lengths as our integers, which is why they do so much better than the paper's and
the demo's signposts. The demo's failures are its own assertions in the common-subdivision recovery
(over both operations: `could not identify source face on mesh A` 5 336, `normal coordinates ... must
sum to an even number` 1 213, `non finite edge length` 8, 40 segmentation faults and 200 time-outs).

**Integer failures: none** (after I30's fix, a589919; the WP10 run had two). No certificate failure
(7 745 IDTs, 7 745 refinements), no invalid subdivision, no trace failure, no overflow, no
out-of-memory. The signposts' 15 failures under Ruppert's rules are trace failures (the same meshes
as in WP10).

The two meshes of the WP10 run, 176193 and 190253 (about 16 000 faces, closed), did not finish their
refinement, with integers and signposts alike. Cause (SPEC I30): mollified near-duplicate input
vertices make Delaunay flips whose shape number (about 2e-7) is below IT-SPEC N5's limit `1e-6`;
refused, they left edges that are not Delaunay (4 in 176193's IDT, 35 flips during 190253's
refinement), and the circumcentres of their triangles are not in empty circles: 4 057 of the first
16 244 insertions on 176193 landed closer than the circumradius to an existing vertex (the reference:
0 of 1 592), each making short edges to refine in turn. Refinement now retries such a flip with only
strict convexity and valid triangles required (geometry-central's `flipEdgeIfNotDelaunay` has no
shape test; by Lawson, a Delaunay flip never lowers the pair's smallest angle). 176193: 1 595
insertions, 25.04° (the demo 1 592); 190253: 21 695, 25.01° (the demo 21 704). `RefineReport::
n_thin_flips` counts these flips: 756 on 105 meshes, and results changed only on those (`|V¹|` by
tens of vertices). `tests/thingi_regressions.rs` runs both (needs `data/thingi10k/converted/`).
Reaching 25° but not on every face: with the default Delaunay test mesh 551019 has 6 face splits
refused by App. A's recovery (`RecoverFailed`, minimum angle 21.17°); with the authors' test it reaches
25°. The other 169 meshes with corners below 25° at the end are narrow cones (angle sum below 60°, the
narrow-cone rule of Chew's rules: 3 000 corners at needle vertices, 871 at angle sums below the
bound). Chew's rules leave the boundary edges encroached by the diametral circles (16 059 edges on
276 of the 406 meshes with boundary: not a failure, the minimum angle is reached).

Counters over the integer runs: 140.07 M vertices inserted, 184 706 removed, 1 266 removals
refused (I8, I27; on 110 meshes), 756 thin flips (I30; 105 meshes), 2 840 crossings clamped in the
refinements (I5; 385 meshes) and 10 338 more in the refined subdivisions, 0 clamps in the IDT subdivisions, 791 places
moved into their face (I16), 7 face splits refused (551019: 6, 97660: 1; the I28 refusal of an edge
split is not counted separately by the driver), 0 trace failures. The Delaunay test of the
reference (`DelaunayCriterion::REFERENCE`) gives the same picture (1 refusal).

Integers against signposts under the same rules: the IDTs have the same flips and the same
`Σ n⁺` on all 7 743 meshes where both completed; after refinement `|V¹|` is equal on 7 712 of 7 728
(the differences are the signposts' snaps onto input edges, I28/N10), and `|V¹|` and `|V^S|` both on 6 734
(WP10 run).

### Sizes

| | `|V^S|/|V⁰|` mean | 95th percentile | `|V¹|/|V⁰|` mean | 95th percentile |
|---|---|---|---|---|
| paper, after refinement | 20 | 45 | 3.7 | 7.8 |
| integer coordinates, after refinement (7 745 meshes) | 19.82 | 46.35 | 3.80 | 8.00 |
| after the IDT | 5.17 | 15.04 | 1 | 1 |

### Against the authors' demo

Built headless (`--noGUI`) with `scripts/ic-demo-build.sh` in `runs_tmp/ic-demo-build/`, geometry-central
`ab116f1`, polyscope `56b82cd`, happly `88f7972`; nothing is written into the reference tree. Needed:
`GC_EIGEN_LOCATION` (Eigen 3.4.0), `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` (polyscope's glad), the demo's
`deps/happly` on the include path (geometry-central's own is an empty submodule). `setMarkedEdges`
is not needed to link; no system package was missing. Run as the paper's `benchmark/run_benchmark.py`
does (`--refineMaxInsertions=0`), on the 7 444 converted meshes up to 200 000 faces.

* **IDT** (integers, 6 004 meshes up to 20 000 faces; ours with the demo's Delaunay test `cotan weight <
  −1e-6`): `|V¹|` equal on all; the subdivision's `|V^S| = |V¹| + Σ n⁺` and its triangulated V/E/F equal
  on 4 899 (81.6 %), and on 1 392 of 1 393 meshes without ties; the rest differ by the direction in
  which a near-tie was flipped (the two flip sequences are not the same, I23). The sets of distinct
  edge lengths agree to 1e-9 relative on 4 581 meshes; 9 516 of 22.5 M distinct lengths (0.04 %) have
  no partner. With the default test (tolerance `1e-10` on the angle sum instead of the demo's `1e-6`
  on the cotan weight) the agreement is lower, which is why the table uses the demo's test.
* **Refinement** (7 440 meshes, statistics, I23): `|V¹|` ours/demo mean 0.9998, median 1.0000, equal on 805,
  within 1 % on 6 531, within 5 % on 7 389; `Σ|V¹|` −0.03 %; `|V^S|` ours/demo mean 0.994, median
  0.9992; `Σ|V^S|` −0.14 %. Minimum angle over valid faces: the demo ends below 25° on 9 meshes
  (56527, 78735, 106830, 106839, 106840, 107202, 264708, 521601, 521603), ours on none (with the
  authors' test).
* Success disagrees on 2 meshes: the demo times out on 65399 and 106831, where ours finishes (in the WP10
  run also on 176193 and 190253, where ours did not finish before I30's fix).
* **Reference tests** ported in `tests/reference_ports.rs` on `fox.ply` and `cat_head.obj`
  (`tests/assets/`, MIT notice): `IntegerFlip`, `DelaunayTriangulationsAgree` (lengths equal to 0, the
  reference's bound is 1e-5), `IntegerTrace`, `IntegerRefine` (fox: 313 → 516 vertices, 25.04°;
  cat_head: 131 → 147, 25.13°, with boundary), `IntegerCommonSubdivision`,
  `CommonSubdivisionCompareIntegerSignpost`, `FunctionTransfer` (L² error 4.3e2 against 9.1e2 for
  pointwise transfer; normal equations to 1e-6). Not ported: the two trace-versus-bulk comparisons
  and `CommonSubdivisionGeometry`.

### Timings

One core (`taskset -c 7`), one job after the other, nothing else of ours running; 1-minute load
average median 1.03, 90th percentile 1.34 (the other agent's work). Seconds; the paper's machine
was an i9-9980XE. 'refine' includes the certificate (integers).

| mesh | paper | ours integers refine | ours signposts refine | demo integers refine | demo signposts refine |
|---|---|---|---|---|---|
| 44395 (Fig. 15, 2 948 vertices) | 0.2 s | 0.029 | 0.021 | 0.084 | 0.031 |
| 54674 (Fig. 16, 9 389) | 0.6 s | 0.096 | 0.064 | 0.305 | 0.101 |
| 719791 (Fig. 19, 871 434) | 4 min / 1.5 min | 58.6 s | 36.1 s | 181.1 s | 50.9 s |
| 176193 (I30, after the fix) | – | 0.048 | 0.028 | 0.17 | – |
| 190253 (I30, after the fix) | – | 0.177 | 0.136 | 0.52 | – |
| 54621 (after the undo journal) | – | 1.16 | 0.604 | 2.14 | – |

719791: `|V¹|` 3 349 214, **`|V^S|` 33 873 704** (the paper: about 30 M and not extractable for memory);
the common subdivision of the integers took 16.1 s and 12.1 GB peak resident (24 GB limit), the demo's 229 s.
The 50-mesh subset (every 155th kept mesh by id; the ids are in `runs_tmp/wp10/timings.md`): totals over the
subset ours integers IDT 1.65 s, refine 8.79 s, subdivision 4.29 s; ours signposts 0.83 / 5.24 / 2.56 s;
demo integers 1.68 / 24.10 / 28.71 s; demo signposts 0.94 / 7.66 s. The demo's integer refinement takes
a median 2.8 times ours (1.3 to 6.7). Integer refinement (with certificate) costs about 1.7 times the signposts' here, the IDT 2 times (the paper: about 2.7 times).

Criterion benches (`cargo bench -p ddg-integer-coords --bench representation`, single thread pinned to a
core, nothing else running; medians in ms; the signposts' own benches are `-p ddg-intrinsic --bench
retriangulation`: spot IDT 1.50, refinement 4.26, camelbad IDT 8.50, refinement 37.2). 1 000 edges flipped and flipped back
cost (2 000 flips) 0.33 ms on every mesh, integers and signposts alike (a flip is O(1) in both).

| mesh (`|V0|`) | IDT ints, certified | IDT ints | IDT signposts | refine ints (Chew, certified) | refine signposts (Ruppert) | subdivision ints | subdivision signposts (traces + overlay) | transfer operators (ints) | transfer solve |
|---|---|---|---|---|---|---|---|---|---|
| spot | 2.95 | 1.56 | 1.50 | 5.88 | 4.22 | 4.11 | 1.63 | 18.6 | 0.27 |
| bob | 3.09 | 1.54 | 1.49 | 4.71 | 3.10 | 3.92 | 1.15 | 16.0 | 0.24 |
| blub | 7.81 | 4.23 | 4.08 | 23.5 | 17.4 | 18.4 | 6.90 | 79.7 | 0.99 |
| camelbad | 15.5 | 8.92 | 8.58 | 50.4 | 35.5 | 28.5 | 15.6 | 190.5 | 2.15 |
| 44395 (Fig. 15, 2 948) | 3.95 | 2.29 | 2.23 | 28.5 | 19.2 | 9.30 | 6.32 | 58.7 | 0.77 |
| 54674 (Fig. 16, 9 389) | 12.6 | 7.51 | 7.21 | 93.7 | 59.3 | 30.9 | 19.9 | 231 | 2.40 |

The certificate roughly doubles the IDT's cost and the integers' subdivision 1.5 to 2.5 times the signposts'
(it also recovers the crossings' geometry from the integers); the integer IDT itself costs the signposts'.


**A slow case, fixed.** Mesh 54621 (159 858 vertices, 1 450 boundary edges, 12 951 removals under
Chew's rules) took 28.9 s with integers, 64.9 s with signposts under Chew's rules, 0.6 s under Ruppert's
(no removals) and 2.1 s in the authors' demo; 439 s in the 56-job run. `perf`: 87 % of the time was
`memcpy`. A removal that might need flips first copied the whole state to restore it on a refusal
(about 20 MB with integers, 50 MB with signposts at 226 000 vertices), O(n) per removal. Now a refused
removal undoes its flips from a journal (`Complex::flip_undo`/`undo_flip` plus the representation's
own data, 16210ca), with bit-identical results: 1.16 s with integers, 0.80 s with signposts under
Chew's rules. The removal schedule (Alg. 10 line 20, Dijkstra ball `< ℓ` after the split's flips;
the reference removes in its split callback with `≤ ℓ`) was not the cost and is kept. Over the whole
Thingi10k run the integer refinement time went from 5 999 s to 3 488 s; the slowest mesh is now
719791.

### Verdict on SPEC §13 Q3

**The gate passes.** Q3 asks for integer success at least the paper's 100 %/100 % minus explained
exceptions such as time-outs or memory, and no certificate failure.

* Integer coordinates: 7 745 of 7 745 meshes complete the IDT and the refinement with a valid common
  subdivision (100.00 % in both columns).
* Certificate failures: **none** (7 745 IDT, 7 745 refinement certificates), no invalid subdivision, no
  trace failure, no overflow, no memory failure; the 30 M-vertex subdivision of 719791 is extracted
  (33.9 M vertices, 12 GB, 58.6 s + 16.1 s, against the paper's 4 min and "not extractable").
* The WP10 run's two exceptions (176193, 190253) were a defect of the shared refinement driver,
  fixed in a589919 (SPEC I30); their tests in `tests/thingi_regressions.rs` pass.
* The integers do better than the paper's signposts (96.0 %/69.1 %) and than our own signposts
  (100.00 %/99.81 % under Ruppert's rules, whose 15 failures are trace failures); the demo's
  signposts reach 87.9 %/20.7 % unmollified.

With the gate passed, the consumers' default is now integer coordinates (Q3):
`Representation::DEFAULT`, the Python `representation=` defaults and Blender's
`REPRESENTATION_DEFAULT`. Signposts stay an option.

## The reference implementation

geometry-central (MIT, © 2017–2019 Nicholas Sharp and the geometry-central contributors), commit
`ab116f1`, pinned by the authors' demo, was read for the details the paper leaves out (map:
`runs_tmp/ic-reference-map.md`); nothing was translated from it. Its notice is in
`LICENSES/geometry-central-MIT.txt`.
