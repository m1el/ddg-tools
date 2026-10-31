 both extractions certify and agree (1378 paths, 6210 points compared, `|Δ| sin θ ≤ 9.5e-14`) |# ddg-flip-geodesics

Geodesic paths, loops and networks on triangle meshes by intrinsic edge flips, after

> Nicholas Sharp, Keenan Crane. *You Can Find Geodesic Paths in Triangle Meshes by Just
> Flipping Edges.* ACM Transactions on Graphics 39(6), Article 249, 2020.
> <https://doi.org/10.1145/3414685.3417839>

run on either intrinsic triangulation of this workspace: the signposts of `ddg-intrinsic`
(Sharp, Soliman, Crane 2019) and the integer coordinates of `ddg-integer-coords` (Gillespie,
Sharp, Crane 2021). Everything is generic over `ddg_intrinsic::Triangulation`.

The contract of the crate is [`SPEC.md`](SPEC.md): scope, data, the geometry of a joint, the
algorithm, verification, and every decision the paper leaves open (`F1`–`F27`).

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
(`max_wedge_mismatch` 7e-11).

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

Next: geodesic Bézier curves and the intrinsic constrained Delaunay triangulation (F4).

## Credits

The algorithm is the paper's. The authors' implementation in geometry-central (MIT,
© 2017–2019 Nicholas Sharp and the geometry-central contributors; notice in
[`LICENSES/geometry-central-MIT.txt`](LICENSES/geometry-central-MIT.txt)) was read for the
details the paper leaves out; the places that follow it are cited `gc: <file>:<line>` in the
docs and in `SPEC.md`.
