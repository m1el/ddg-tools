# ddg-flip-geodesics

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

Next: waypoints and edge sets, extraction on the input mesh for both representations (F3);
geodesic Bézier curves and the intrinsic constrained Delaunay triangulation (F4).

## Credits

The algorithm is the paper's. The authors' implementation in geometry-central (MIT,
© 2017–2019 Nicholas Sharp and the geometry-central contributors; notice in
[`LICENSES/geometry-central-MIT.txt`](LICENSES/geometry-central-MIT.txt)) was read for the
details the paper leaves out; the places that follow it are cited `gc: <file>:<line>` in the
docs and in `SPEC.md`.
