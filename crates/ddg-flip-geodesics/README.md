 both extractions certify and agree (1378 paths, 6210 points compared, `|Δ| sin θ ≤ 9.5e-14`) |# ddg-flip-geodesics

Geodesic paths, loops and networks on triangle meshes by intrinsic edge flips, after

> Nicholas Sharp, Keenan Crane. *You Can Find Geodesic Paths in Triangle Meshes by Just
> Flipping Edges.* ACM Transactions on Graphics 39(6), Article 249, 2020.
> <https://doi.org/10.1145/3414685.3417839>

run on either intrinsic triangulation of this workspace: the signposts of `ddg-intrinsic`
(Sharp, Soliman, Crane 2019) and the integer coordinates of `ddg-integer-coords` (Gillespie,
Sharp, Crane 2021). Everything is generic over `ddg_intrinsic::Triangulation`.

The contract of the crate is [`SPEC.md`](SPEC.md): scope, data, the geometry of a joint, the
algorithm, verification, and every decision the paper leaves open (`F1`–`F30`).

## What is implemented (work package F1)

- Edge networks: paths and loops as sequences of directed edges (`Segment`, `EdgePath`),
  `closed` explicit, pinned vertices (path ends pinned implicitly), multiply covered edges with
  the layering of the segments that share an edge (paper Fig. 11; default as the reference,
  or given explicitly), stable segment ids.
- Validity (`Network::check`): connected paths, closed loops, pinned ends, a complete layering
  and non-crossing at every vertex (the chords of the joints, in the cyclic order of the
  segment ends, do not interleave). Crossing inputs are refused, never repaired.
- Initial paths: through a vertex sequence, and Dijkstra on the edge graph of the current
  triangulation with ties by `(distance, halfedge id)`.
- The geometry of a joint: both wedge angles summed corner by corner (Eq. 1), boundary wedges
  (`∞`), the Δ-complex cases (a self-edge incident on a wedge twice, the out-and-back joint at a
  vertex of degree one with a wedge of degree zero), the classification straight / left /
  right with the reference's tolerance, flexible wedges and wedges blocked by a path end.
- The geodesic certificate on the triangulation (SPEC §7 layer 1a).
- The options and report types of the shortening.

## What is implemented (work package F2)

- **Alg. 1** (`FlipOut(T, γ_abc)`, SPEC §4.1, `src/flipout.rs`): the wedge's outer arc is
  made an edge path by flipping the first edge `b n_j` with `β_j < π`, resuming at `j − 1`
  after a flip; refused flips are skipped and counted, the replacement is made only if
  strictly shorter (else the joint *stalls*), new segments go innermost on the outer-arc edges.
  The Δ-complex cases of App. A need no special code (the wedge is recomputed from the
  connectivity after every flip); an out-and-back joint (`deg W = 0`) contracts; a
  single-segment loop is moved across the face on its wedge's side (App. B).
- **Alg. 2** (`MakePathGeodesic`, SPEC §4.2–§4.4, `src/shorten.rs`):
  `GeodesicNetwork::shorten(&ShortenOptions)` on paths, loops and networks: the queue
  ordered `(α, side, path, segment)`, stale and blocked entries, stalls, a final sweep after
  the queue empties, the stop rules `max_iterations` (default `100 (|E| + m₀)`),
  `stop_length_ratio`, `stop_min_angle`, the running length, contracted loops reported.
  `shorten_with` calls an observer after every iteration; `shorten_joint` runs Alg. 1 once.
- **Report** (`ShortenReport`): status, iterations, flips, refused flips, stalls, sweeps,
  App. B moves, contractions, the App. A cases met, stale and blocked entries.

Every test run ends with the certificate (SPEC §7 layer 1a) and, with `check_invariants`,
checks the network (non-crossing included) after every iteration.

| test (SPEC §7) | file | result |
|---|---|---|
| 1 flat convex domains: straight segment | `tests/flat.rs` | 250 paths on grids, sheared and jittered flat domains (input and Delaunay): length = Euclidean distance to 8.0e-15 relative, every vertex on the segment |
| 2 flat L: visibility distance, reflex corner kept | `tests/flat.rs` | 300 paths: visibility-graph distance to 9.3e-16 relative; 22 bend at the reflex corner and keep it |
| 3 pyramid cone, cube `√5 a`, prism loop, `3π` saddle, boundary corners | `tests/closed_form.rs`, `tests/flat.rs` | pyramid 1.6e-16, cube 5.7e-16, prism loop 1.5e-16 (relative, `≤ 1e-12` asserted); the `3π` saddle path unchanged with 0 flips, 0 iterations; the convex corner is left |
| 4 exact distance (`ddg_intrinsic_apps::exact`) | `tests/random_paths.rs` | 360 paths on 6 surfaces: never below `d_exact (1 − 1e-9)`; 80.6 % minimal to `1e-6`, 89.4 % to `1e-3` (the rest are other certified local geodesics, e.g. around a torus the long way) |
| 5 monotone length | every run (`tests/common`) | each replacement strictly shorter; only App. B moves lengthen |
| 6 isotopy: Fig. 8 top (spiral, two loops), Fig. 8 bottom, networks | `tests/isotopy.rs` | valid after every iteration |
| App. A cases | `tests/isotopy.rs`, `tests/loops.rs` | over 11 136 runs: 18 wedges with a twice-incident edge, 18 outer arcs through a degree-one vertex, 18 out-and-back contractions, 44 App. B moves |
| 8 representation identity (paths, loops, networks) | `tests/identity.rs` | identical (see below) |
| 10 determinism, refused inputs | `tests/random_paths.rs` | bitwise identical, serial and 16 runs on parallel threads |

Evaluation on `data/` meshes (`tests/random_paths.rs`, ignored test, release, one thread):

| mesh | V | geodesic (`1e-4` rad) | FlipOut ms median / p99 / max | FlipOut / Dijkstra | flips median / max | m median | refused flips | minimal to `1e-6` |
|---|---|---|---|---|---|---|---|---|
| spot | 2930 | 200/200 | 0.27 / 1.11 / 1.79 | 1.06 | 101 / 759 | 23 | 4 | 160/200 |
| bob | 5344 | 200/200 | 0.56 / 1.96 / 2.50 | 1.03 | 237 / 1131 | 32 | 8 | 180/200 |
| blub | 7106 | 200/200 | 0.67 / 5.69 / 5.85 | 1.12 | 279 / 2544 | 36 | 6 | 177/200 |
| CrumpledDevelopable | 289 | 200/200 | 0.06 / 0.39 / 0.56 | 2.66 | 24 / 220 | 9 | 0 | 200/200 |
| camelbad | 11392 | 200/200 | 0.91 / 5.55 / 6.81 | 0.76 | 331 / 1914 | 47 | 2 | 130/200 |
| bunny | 14290 | 200/200 | 0.83 / 3.75 / 3.99 | 0.48 | 318 / 1549 | 49 | 1 | 122/200 |
| cow1 | 15049 | 200/200 | 0.76 / 2.87 / 3.62 | 0.51 | 299 / 1473 | 50 | 0 | 183/200 |
| rhino | 20205 | 200/200 | 1.09 / 3.66 / 6.29 | 0.50 | 435 / 2562 | 57 | 7 | – |
| tree | 16917 | 200/200 | 0.97 / 3.35 / 4.03 | 0.54 | 405 / 1769 | 54 | 21 | 191/200 |

1800 of 1800 trials certified geodesic (paper Sec. 5.3: 99.4 % on its datasets, 100 % on
MPZ); one stall (blub); no final sweep found work; flips grow as `m^1.59` (paper Fig. 16:
`m^1.5`). Times exclude the clone of the triangulation (F16) and extraction (F3; the paper's
times include it: "a few ms"). Single thread; the machine was shared with other agents' jobs
(load average 4.6–10 on 64 cores). Run with `cargo test --release -p ddg-flip-geodesics
--test random_paths -- --ignored --nocapture --test-threads=1`.

On both representations (`tests/identity.rs`, SPEC §7 test 8): the same iterations, flips,
segments and bitwise-identical lengths on signposts and on `ddg_integer_coords`, from the
input and from the Delaunay triangulation (self-edges, degree-one vertices), with the integer
certificate after every run.

## What is implemented (work package F3)

- **Initial paths** (`src/init.rs`, SPEC §2.2): `waypoint_path` (Dijkstra between consecutive
  waypoints, concatenated, waypoints pinned or not; gc: fg.cpp:283-312) and `network_from_edges`
  (the reference's walk, gc: fg.cpp:315-396: vertices of path degree `≠ 2` and `extra_pinned` are
  ends and are pinned, paths between ends, the rest closed loops; boundary edges are walked in
  both directions, which the reference's halfedge loop cannot). Both return `InitialPaths { paths,
  pinned }` for `GeodesicNetwork::new`. Errors name the leg (`WaypointLeg`), the mask size or the
  vertex.
- **Extraction** (`src/extract.rs`, SPEC §5): `GeodesicNetwork::extract(&ExtractOptions)` is one
  `PathOnInput` per path, segments in path order, each the `trace_edge_on_input` of its edge
  (reversed for a backward segment, the joint written once; a closed path repeats its first
  point at the end). The distinct edges are traced in parallel (rayon) and assembled in order:
  bitwise identical at any thread count. A failed trace is `FlipError::Extraction` naming path,
  segment and edge. `extract_with_report` adds the per-edge counts, `segment_starts` and the
  certificate; `polylines_3d(positions)` and `polyline_3d` interpolate barycentrically;
  `crossing_sines` gives the conditioning `sin θ` of each crossing position.
- **Certificate on `T⁰`** (`InputCertificate`, SPEC §7 layer 1b), from the input mesh's lengths
  and corner angles alone: the angle at every crossing of an input edge is `π` (second face unfolded
  across the edge; tolerance `1e-6`); at every joint both sides measured on `T⁰` are `≥ π − 1e-4`
  unless the joint is pinned or blocked by an end on `T`; the sides are also compared with the wedge
  angles `T` computed (`max_wedge_mismatch`); the polyline's length equals the edge path's within
  `1e-9` relative. The per-edge `length_mismatch` of the integer extraction is
  `ExtractReport::edges[..].mismatch`.

### Results (`tests/extract.rs`)

| test | result |
|---|---|
| waypoints, edge sets | legs equal the concatenated Dijkstra paths; tree, tree-cotree cut graphs, random edge sets, boundaries: every masked edge once, ends exactly the vertices of degree `≠ 2`, valid networks; a Dijkstra path or a tube ring round-trips |
| certificates on `T⁰`, both representations | 223 paths on 17 generated and data shapes (random pairs, three at a time, waypoint loops, self-edge loops, tube loops): all pass; crossing deviation `≤ 2.6e-11` (signposts), `1.3e-12` (integers); length error `≤ 5e-15` |
| signposts vs integers (§7 layer 2) | same crossed input edges at all 3401 compared points; `|Δ| sin θ ≤ 5.8e-13` (raw `|Δ|` `≤ 3.9e-11`) |
| seam networks | tree-cotree cut graph of a torus, BFS trees of a sphere and of spot, shortened with junctions pinned: both extractions certify and agree |
| determinism (§7 test 10) | `Extraction` equal, `t` bitwise, at 1, 3 and 16 threads |

Hunt for signpost trace failures (SPEC §7 test 9): 9 data meshes x 4 seeds x 1000 random pairs
(36 000 trials, `FO_MESH`, `FO_SEED`, `FO_TRIALS` of the ignored `evaluation` test): **no
signpost trace failed**, so there is no committed fixture yet (the failures of IntCoords Fig. 17
need Thingi10k, F7). Two trials (tree, seed 2, pairs 372 and 934) ended with FlipOut not
geodesic to `1e-4` (the paper's 99.4 %); the certificate on `T⁰` measured the same angle as `T`
(`max_wedge_mismatch` 7e-11). Both were IT N5's shape test refusing a needle flip; with the
reference's relative test, now the default, they certify (see "Runs that stall" below).

Timings (`evaluation`, release, one thread, 200 random pairs per mesh, median per path; the
machine was shared: load average 1.4 before, 4.3 after the run):

| mesh | V | Dijkstra ms | FlipOut ms | extraction ms, signposts | extraction ms, integers | (FlipOut + extraction) / Dijkstra |
|---|---|---|---|---|---|---|
| spot | 2930 | 0.30 | 0.28 | 0.010 | 0.021 | 1.00 |
| bob | 5344 | 0.79 | 0.62 | 0.010 | 0.027 | 0.82 |
| blub | 7106 | 0.79 | 0.62 | 0.012 | 0.030 | 0.82 |
| CrumpledDevelopable | 289 | 0.026 | 0.052 | 0.002 | 0.005 | 2.2 |
| camelbad | 11392 | 1.40 | 0.87 | 0.028 | 0.057 | 0.67 |
| bunny | 14290 | 2.03 | 0.79 | 0.047 | 0.082 | 0.43 |
| cow1 | 15049 | 1.86 | 0.84 | 0.043 | 0.086 | 0.50 |
| rhino | 20205 | 2.84 | 1.17 | 0.037 | 0.078 | 0.44 |
| tree | 16917 | 2.52 | 1.10 | 0.046 | 0.089 | 0.47 |

Extraction costs up to 10 % of FlipOut on the integers and half of that on signposts (the trace of
an edge is the cheap part; the integer walk identifies each crossing). The paper's total
(Dijkstra + FlipOut + extraction) is about 1.1 times Dijkstra's time (Fig. 18); ours is
`1 + (0.4 .. 1.0)` on the large meshes, i.e. 1.4-2 times: FlipOut costs about as much as our
Dijkstra on small meshes. In absolute terms every path takes under 2 ms on one thread (paper:
"a few ms" on an i7-4790K).

## What is implemented (work package F4)

- **Geodesic Bézier curves** (`src/bezier.rs`, SPEC §6.1, paper Sec. 6.2, Fig. 23, after
  Morera, Carvalho, Velho 2008): `geodesic_bezier(tri, control_vertices, rounds, &options)`
  returns the network (one open path, the curve) and a `BezierReport`;
  `GeodesicNetwork::bezier_subdivide` works on an existing open path whose pinned vertices
  are the control points. It follows the reference's recursion (`gc: fg.cpp:1068-1248`):
  - shorten;
  - split every part between consecutive control points at arc length `L/2`, snapping to
    either end of the edge within `1e-4 ℓ` (F19);
  - unpin the old interior control points, shorten, and repeat on the new points until one
    part is left;
  - recurse on both halves.

  After `r` rounds the pinned points (`BezierReport::control_points`, `d 2^r + 1` of them)
  are the control polygons of the `2^r` pieces. The control points are kept as a list and
  checked against the pinned vertices along the path after every change. A curve that meets
  a control point twice is refused (`FlipError::Bezier`); the reference assumes the curve is
  simple. Two retries the reference does not need, both counted in the report:
  - **thin split**: when the triangulation refuses a split as too thin (IT N5), the
    triangulation is flipped to Delaunay with the curve's edges held fixed, and the split is
    made again. On a flat square a straight geodesic is one long edge between needles.
  - **point on an input edge**: the integer split refuses a point on an input edge
    (`NotLocated`), and lattice-like meshes give many. A midpoint within `1e-9 ℓ` of a
    crossing of its edge with an input edge (`trace_edge_on_input`, on both
    representations) is moved by `2e-9 ℓ` towards the middle of the edge.
- **Edge splits that carry the network** (`GeodesicNetwork::split_edge`, `network.rs`):
  every segment on the split edge, however many share it, becomes two segments on the halves
  (`Insertion::halves`), each in its old place across the edge. The first part keeps the
  segment's id.
- **The iCDT** (`src/cdt.rs`, SPEC §6.3, paper Sec. 6.3, Fig. 25): `make_delaunay(&options)`
  marks every edge that carries a segment as constrained and runs the shared
  `flip_to_delaunay` (`gc: fg.cpp:1035-1045`). It returns a `CdtReport`: the
  `DelaunayReport`, the number of network edges and the number newly constrained. The marks
  stay, as in the reference; `clear_constraints` removes them before the network changes
  again.
- **Certificate on `T⁰` at inserted vertices** (`src/extract.rs`): a joint at an inserted
  vertex (inside an input face, or on an input edge) is measured in that face, and the face
  after it is unfolded across the edge. The joint is checked like the others; before, such a
  bend on an input edge was taken for a crossing. A crossing closer than `1e-8` of the face's
  size to such a joint belongs to the joint. The integer split can put a vertex within
  rounding of an input edge (barycentric `7e-16`), and the angle between pieces that short
  is rounding.

### Results (`tests/bezier.rs`, `tests/cdt.rs`, `tests/identity.rs`)

| test | result |
|---|---|
| 7, plane: pinned points = planar de Casteljau subdivision | jittered flat squares 12 × 12 and 24 × 24; `d = 1 … 4`; 0–4 rounds; input and Delaunay start; 75 cases. Places on `T⁰`: integers `≤ 2.2e-14`, signposts `≤ 1.2e-11`. Part lengths: `≤ 1.2e-14`. The bound is `1e-12` for integers and part lengths, `1e-10` for signpost places (the trace of a new vertex carries the rounding of its edge's signpost after many flips) |
| 7, curved: reversing the control points | icosphere, jittered icosphere, torus, bumpy grid, spot, camelbad; `d = 3`; 3 rounds. Signposts: 22 cases, `≤ 1.3e-13` of the bounding-box diagonal. Integers without a nudge: 20 cases, `≤ 4.6e-15`. Two integer cases had a midpoint exactly on an input edge: the two directions split different edges there and the nudges differ, `2.5e-10` (bound `1e-8`) |
| Bézier on both representations | both certificates after every curve (`d = 2`, 4 rounds); reports, paths and extraction identical; deterministic; extraction equal at 1 and 8 threads |
| 8, identity, Bézier + iCDT | 53 curves on 15 shapes (input and Delaunay start): same splits, flips, segments, iCDT flips, faces and lengths bit for bit; integer `certify()` after each |
| splits of shared edges | three paths side by side on the same edges, two running one way and one the other: after each split both halves list the segments in the old order; valid network; joints straight; layering rebuilt equal. A single-segment loop on a self-edge is split into a two-segment loop. Refused splits change nothing |
| snapping (F19) | a grid row, 3 rounds: all 7 midpoints snap to grid vertices, nothing is split |
| iCDT | 59 networks on 16 shapes, plus Bézier curves on spot, camelbad and an icosphere, both representations: paths and lengths unchanged bitwise; every unconstrained interior edge Delaunay (status `Delaunay`); both certificates; integer `certify()`; representations identical; deterministic |
| refused inputs | too few or repeated control points, too many rounds, bad snap tolerance, a loop or two paths given to `bezier_subdivide`: errors; nothing changed |

Timings (`tests/bezier.rs::timings`, ignored; release, one pinned core, run alone; three runs
within 1 %; load average 3.2 from other agents' jobs). Ten curves per mesh with `d = 3` and
4 rounds (90 splits each), then the iCDT of each. Medians:

| mesh | Bézier ms, signposts / integers | iCDT ms, signposts / integers | flips per curve: Bézier / iCDT | plain Delaunay of the input (with building it) ms, sp / ic |
|---|---|---|---|---|
| spot | 0.80 / 0.88 | 1.81 / 1.82 | 149 / 567 | 1.85 / 2.07 |
| camelbad | 3.23 / 3.43 | 9.6 / 9.7 | 451 / 8559 | 11.8 / 13.2 |

The iCDT costs what Delaunay flipping of the whole mesh costs; the curve's 50 edges are a
small part of it.

## Speed of the shortening, and runs that stall

**Speed.** `shorten` keeps its working memory in `src/scratch.rs` and no longer allocates per
call. Each item below keeps every result bitwise:

- **Corner angles** are memoized per face (F2's cache). Alg. 1 invalidates the two faces of
  every flip, and debug builds check every hit against the face's lengths.
- **Interior joints** whose two segment ends lie on different edge ends are walked around their
  vertex (`h → twin(prev(h))`). Their flexibility is tested as `gc: fg.cpp:589-658` does. The
  star stays the fallback, and debug builds assert both agree.
- **β tests:** Alg. 1 decides `β_j < π` from Kahan's factors without computing angles, since
  `θ₁ + θ₂ − π = 2 atan((P − 1)/(T₁ + T₂))` with `T = tan(θ/2)` and `P = T₁T₂`. The angles are
  summed only within `1e-13` of `π`, and debug builds assert the decision.
- **After a flip** Alg. 1 merges the two corners of the flipped edge instead of rebuilding the
  wedge, and checks the Δ-complex cases (App. A) only where they can change.
- **Flips** call `flip_edge_with` once, with a `FlipRequest` for the shape test and the corner
  angles; it refuses exactly when `can_flip_with` would.
- **The flip's corner angles** (`Flip::new_corner_angles`, IT N52): Alg. 1 requests the two
  corner angles at the new edge's ends (`FlipRequest::corner_angles`; no other caller does) where
  they are free (`Triangulation::flip_corner_angles_are_free`: the signposts), and they go into
  the memo instead of being recomputed. The integer coordinates are not asked: there the two
  angles cost more than the memo saves.

`examples/fo_pin` fingerprints every result: lengths, counts, every `ShortenStep`, the final
triangulation and the layering. The fingerprints of 9 meshes × 500 pairs × 2 representations are
unchanged by all of this.

`start_shorten` + `ShortenSession::run` is `shorten` split after the initial queue.
`examples/fo_pairs` prints the queue's share as a last column, `queue_ms`. The reference builds
its queue when the network is made, so its `iterativeShorten` time compares with
`shorten − queue`.

Same pairs as the reference driver of SPEC §7 layer 5 (geometry-central at `ab116f1`). Release
build, one pinned core, everything run one after another; load average 2.8–3.6 from other jobs
on the machine. Columns: "sp" = signposts, "ic" = integers; "before" = commit `e8d9682`.

| mesh | Dijkstra ms gc / sp / ic | shorten ms gc / sp before→after / ic before→after | loop alone ms sp / ic | extract ms gc / sp / ic | shorten ratio sp/gc per pair (median, p90) | loop ratio sp/gc (median) | same length as gc (1e-9) sp / ic | geodesic gc / sp / ic |
|---|---|---|---|---|---|---|---|---|
| spot | 0.489 / 0.306 / 0.299 | 0.0434 / 0.2662→0.0683 / 0.2454→0.0669 | 0.0592 / 0.0571 | 0.0289 / 0.0149 / 0.0237 | 1.57, 1.81 | 1.34 | 500 / 500 of 500 | 500 / 500 / 500 |
| bob | 0.845 / 0.517 / 0.518 | 0.0859 / 0.5005→0.1276 / 0.4979→0.1262 | 0.1142 / 0.1153 | 0.0412 / 0.0101 / 0.0237 | 1.48, 1.65 | 1.34 | 500 / 500 of 500 | 500 / 500 / 500 |
| blub | 1.216 / 0.735 / 0.729 | 0.1129 / 0.6407→0.1595 / 0.6268→0.1590 | 0.1424 / 0.1427 | 0.0565 / 0.0129 / 0.0277 | 1.40, 1.51 | 1.27 | 500 / 500 of 500 | 500 / 500 / 500 |
| crumpled | 0.040 / 0.026 / 0.026 | 0.0104 / 0.0633→0.0175 / 0.0618→0.0172 | 0.0138 / 0.0137 | 0.0108 / 0.0021 / 0.0045 | 1.69, 2.14 | 1.31 | 498 / 498 of 498 | 498 / 498 / 498 |
| camelbad | 2.163 / 1.312 / 1.337 | 0.1581 / 0.8091→0.2028 / 0.8052→0.2018 | 0.1829 / 0.1843 | 0.0794 / 0.0199 / 0.0403 | 1.29, 1.36 | 1.18 | 500 / 500 of 500 | 500 / 500 / 500 |
| bunny | 2.955 / 1.649 / 1.673 | 0.1632 / 0.7283→0.1845 / 0.7184→0.1846 | 0.1651 / 0.1658 | 0.1006 / 0.0267 / 0.0534 | 1.13, 1.20 | 1.01 | 500 / 500 of 500 | 500 / 500 / 500 |
| cow1 | 2.869 / 1.717 / 1.699 | 0.1740 / 0.8151→0.2061 / 0.8004→0.2040 | 0.1837 / 0.1816 | 0.0976 / 0.0303 / 0.0498 | 1.18, 1.26 | 1.05 | 500 / 500 of 500 | 500 / 500 / 500 |
| rhino | 3.822 / 2.435 / 2.408 | 0.2261 / 1.1484→0.2835 / 1.1205→0.2788 | 0.2607 / 0.2590 | 0.0972 / 0.0308 / 0.0570 | 1.25, 1.33 | 1.14 | 500 / 500 of 500 | 500 / 500 / 500 |
| tree | 3.536 / 2.116 / 2.109 | 0.2259 / 1.0467→0.2577 / 1.0297→0.2564 | 0.2369 / 0.2355 | 0.1036 / 0.0263 / 0.0539 | 1.14, 1.21 | 1.04 | 499 / 499 of 500 | 500 / 499 / 499 |

The shortening is 3.6–4.1 times faster than before on both representations. It is still
1.13–1.69 times the reference per pair, or 1.01–1.34 times for the loop alone.

- **Corner angles.** The reference reads wedge angles from signposts in `O(1)`. Here they are
  sums of corner angles (F2), which need the corners at every vertex a replacement puts on the
  path; those corners are about a third of the remaining time on spot.
- **The flips themselves.** The triangulation's flips (`plan_flip`, the flip, two corner angles
  for the signposts) are about a quarter of it.

**The flip test and the flip's corner angles** (IT N52). Alg. 1's flips now use the reference's
relative shape test (`ShortenOptions::flip_shape`, below), and the memo takes the two corner
angles each flip returns. Same harness and pairs, one pinned core, three interleaved rounds
(median per pair over the rounds, then over the pairs), load average 2.8–7.7 from other jobs;
"before" = commit `3afbe2c`, "after" = both changes; the gc column is the reference's from the
table above.

| mesh | gc shorten | sp shorten before → after | sp loop before → after | ic shorten before → after | ic loop before → after | loop / gc sp, ic (median per pair) |
|---|---|---|---|---|---|---|
| spot | 0.0434 | 0.0679 → 0.0660 | 0.0578 → 0.0569 | 0.0669 → 0.0692 | 0.0577 → 0.0599 | 1.28, 1.38 |
| bob | 0.0859 | 0.1297 → 0.1234 | 0.1155 → 0.1101 | 0.1268 → 0.1338 | 0.1150 → 0.1221 | 1.28, 1.42 |
| blub | 0.1129 | 0.1607 → 0.1533 | 0.1436 → 0.1362 | 0.1574 → 0.1684 | 0.1417 → 0.1503 | 1.22, 1.35 |
| crumpled | 0.0104 | 0.0180 → 0.0173 | 0.0141 → 0.0135 | 0.0174 → 0.0185 | 0.0139 → 0.0148 | 1.29, 1.41 |
| camelbad | 0.1581 | 0.2063 → 0.1968 | 0.1849 → 0.1740 | 0.2018 → 0.2167 | 0.1826 → 0.1960 | 1.13, 1.26 |
| bunny | 0.1632 | 0.1867 → 0.1796 | 0.1681 → 0.1594 | 0.1848 → 0.1941 | 0.1666 → 0.1748 | 0.97, 1.07 |
| cow1 | 0.1740 | 0.2073 → 0.1973 | 0.1848 → 0.1739 | 0.2036 → 0.2147 | 0.1824 → 0.1923 | 1.00, 1.10 |
| rhino | 0.2261 | 0.2895 → 0.2778 | 0.2653 → 0.2536 | 0.2828 → 0.2999 | 0.2618 → 0.2783 | 1.10, 1.22 |
| tree | 0.2259 | 0.2591 → 0.2479 | 0.2369 → 0.2268 | 0.2571 → 0.2719 | 0.2357 → 0.2525 | 1.00, 1.11 |

- Signposts: shortening −2.7 % to −4.9 %; the loop alone is 0.97–1.29 times the reference.
- Integers: +3.4 % to +7.4 %. Their flip computes the two corner angles only for `Flip` (the
  signposts compute them anyway for their signpost update), and that costs more than the memo
  saves. Since then a flip returns the angles only when asked (`FlipRequest::corner_angles`,
  IT N52), so Delaunay flipping, refinement and the other operations no longer pay for them
  (integer Delaunay flipping on camelbad 9.47 → 8.89 ms). Alg. 1 asks only where the angles are
  free (the signposts); against the code before the angles (`adf9d80^`, N5 on both sides, best
  of 4 interleaved rounds, load 2.2–3.0), shortening is −4.2 % to −6.7 % on the signposts and
  +0.8 % to +2.6 % on the integers (+2.4 % to +8.2 % if Alg. 1 asked there too).
- The relative test alone moves the times by −2.2 % to +1.3 %.
- All 4498 pairs now have the reference's length to `1e-9` on both representations, and all
  are geodesic (tree 12804 → 1144 was the one exception before).

**Runs that stall** (`tests/stalled.rs`).
- On the `tree` mesh, three pairs ended with a joint bent by `6e-4`–`4e-3` rad under IT N5's
  shape test: 12804 → 1144, and pairs 372 and 934 of seed 2 (F3's open trials).
- Each Dijkstra path passes two vertices `1e-4`–`6e-4` apart. Straightening the joint next to
  them needs a last flip that makes a needle triangle of shape `1.2e-7`–`6.9e-7`.
- N5 refuses it (`FLIP_SHAPE_LIMIT = 1e-6`). The outer arc is then longer than the joint, and
  the joint stalls (F6). Such a run ends with status `Stalled`, not `Geodesic`.
- The reference's test is relative, `A ≥ 1e-6 (A₁ + A₂)` (ratios `5e-5`–`1e-3` here), so it
  makes the flip and its paths skip both close vertices.
- **Now the default** (`ShortenOptions::flip_shape = FlipShapeTest::REFERENCE`): the three
  pairs certify, end `Geodesic`, and have the reference's length to `1e-12` on both
  representations (`tree_pairs_certify_with_the_reference_flip_test`).
- With `flip_shape = FlipShapeTest::N5` the results are bitwise those from before (9 × 500 pairs
  × 2 representations), and the three pairs still stall
  (`tree_pairs_with_a_refused_needle_flip_stall_under_n5`).
- What the default changes, against N5: 99 of the 4500 pairs, identically on both
  representations; every one certifies. Mostly a flip or two more and a length within `1e-10`
  (tree 1938 → 8644 now ends on the reference's geodesic, `4.6e-8` longer); two statuses change
  from `Stalled` to `Geodesic` (tree 12804 → 1144 and 10499 → 1703); on camelbad
  8294 → 9415 the relative test refuses one flip that N5 makes, with the same length. Only
  Alg. 1's flips use the relative test; the triangulation's own operations (and the iCDT, which
  is Delaunay flipping) keep N5.

## Credits

The algorithm is the paper's. The authors' implementation in geometry-central (MIT,
© 2017–2019 Nicholas Sharp and the geometry-central contributors; notice in
[`LICENSES/geometry-central-MIT.txt`](LICENSES/geometry-central-MIT.txt)) was read for the
details the paper leaves out; the places that follow it are cited `gc: <file>:<line>` in the
docs and in `SPEC.md`.
