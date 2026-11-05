# WP8 notes (continuous collision detection)

Build/test: `CARGO_TARGET_DIR=$PWD/target-wp8 cargo test --release -p ddg-repulsive-surfaces --test ccd`
(16 tests, ~15 s release, random scenes run in parallel). Bench (run alone): `cargo bench -p
ddg-repulsive-surfaces --bench ccd [-- quick] > runs_tmp/bench_full.csv`. clippy `--all-targets` clean.

## API (src/ccd.rs) for WP5
```rust
let pool = energy::exact::make_pool(threads)?;
let obstacles: Vec<CcdObstacle> = /* CcdObstacle::new(&Obstacle, leaf_size) once per run */;
let opts = CcdOptions { reference_length: Some(l0), ..Default::default() };   // gap 0.05, 100 iterations
let rep = first_contact(&pool, &surface, &positions, &velocity, tau0, &obstacles, &opts)?;
let tau = tau0.min(2.0 / 3.0 * rep.tau_max);          // SPEC 7.3; rep.tau_max == tau0 without contact
```
`velocity` is the direction of motion (`-x` for the flow). `CcdReport { tau_max, contact:
Option<Contact{tau, primitive}>, gap_abs, iterations, pairs_tested, windows, close_pairs,
touching_pairs, seconds }`. `Primitive`: `VertexTriangle`, `EdgeEdge`, `VertexObstacleTriangle`,
`TriangleObstacleVertex`, `EdgeObstacleEdge`. Pass `reference_length = Some(l0)` (mean edge
length of the INITIAL mesh); `None` takes the current mean edge length.
Also: `min_distance`, `intersects`, `intersecting_pairs`, `triangles_intersect`,
`triangle_distance`, `vertex_triangle_advance`, `edge_edge_advance` (tagged `rs-ccd-advancement`),
`CcdCollapseGuard::new(&opts, l0)` (implements `remesh::CollapseGuard`; one per `remesh` call;
`refused()`, `pairs_tested()`).

## Decisions
- Gap = `opts.gap * reference_length` absolute. A pair with `d0 < gap` is measured against its own
  distance: target `close_factor * d0` (0.5), so such pairs never freeze the flow (counted in
  `close_pairs`). `d0 = 0` pairs (vertex exactly on a face) are not constrained (`touching_pairs`).
- Lipschitz bound `l`: cheap `max speed(A) + max speed(B)` (precomputed per-vertex speeds) first;
  only if the pair is not cleared, the better of {rest frame, mean-velocity frame} (as the curves
  crate). Both are the SPEC bound (relative motion `max|u_p| + max|u_tri|`).
- Windows `[0,t1],[t1,2t1],[2t1,4t1],..`, `t1` = time for the fastest vertex to move one mean edge
  length; the first window with a contact gives `tau_max` (min over its pairs, tie-break by
  primitive order). Per query, later pairs are advanced only up to the best time found so far.
- Broad phase is NOT the curves `Bvh<3>`: its build is sequential (11 ms of a 46 ms CCD on 10k
  faces, no gain from threads) and it has no per-element box test. `BoxTree`: Morton sort of the
  box centres (parallel), blocks of 8, complete binary tree of unions; per-element box test at the
  leaves. This cut the candidate pairs 5x. Leaf size 2..16 makes no difference.
- Adjacent = shares a vertex (vertex v vs face f with v in f; edges sharing a vertex). A vertex
  opposite an adjacent face is tested (as in the SPEC); on a sharp crease it can cost steps.
- `intersects` = no face pair without a common vertex touches, and for faces sharing exactly one
  vertex the opposite edge does not meet the other face; faces sharing an edge are not tested.
- Collapse guard: for each new star triangle, all triangles (found through a vertex hash grid kept
  up to date by watching which collapses were applied) that share no vertex with it and are not
  incident to the collapsed edge: refuse if they intersect or are closer than `gap_abs`. A pair
  that was already closer than the gap before the collapse therefore also blocks it (stacked
  sheets closer than the gap: all collapses between them are refused).

## Measurements (release; machine load 10-18 from other agents, so +-20 %)
Bench `benches/ccd.rs`, CSV in `runs_tmp/bench_full.csv`. `t_energy` = Barnes-Hut (theta 0.5,
leaf 8) energy + differential via `energy::energy`; `t_ccd` = `first_contact`, `tau_end = 1`,
velocity scaled so the farthest vertex moves one mean edge length (approach_x25: 2.5).
Ratio = t_ccd / t_energy, median over the 21 (mesh, motion) cases, 5k-50k faces:

| set | 1 thread | 8 threads |
|---|---|---|
| all (21) | 0.29 | 0.19 |
| gradient motion (10) | 0.21 | 0.15 |
| near contact: sheet pairs 3 edge lengths apart approaching, torus squeezed (11) | 0.36 | 0.25 |
| range | 0.14 - 0.79 | 0.09 - 0.57 |

Largest: tight torus (R=1, r=0.8) squeezed to its axis (0.67-0.79 single thread), sheets x2.5
approach on 5k faces (0.66). At 50k faces the gradient motion costs 0.14 (1 thread) / 0.09 (8).
Where the time goes: ~15 % tree build, the rest the per-primitive queries (mostly rejecting
adjacent edges that share a vertex) and 1-3 distance evaluations per surviving pair.

## Recommendation for the default (SPEC 7.3, 10 %-of-a-step rule)
Against the energy + gradient alone the median is 0.19-0.29, above 10 %. A full step is more than
the energy + gradient: the fractional solve (GMRES with Barnes-Hut/hierarchical products, tens of
operator applications) and the line-search trials; the paper's cost is dominated by the solve. A
step is therefore plausibly 5-20 times the energy + gradient, which puts CCD at 1-5 % of a step.
I could not measure whole steps (WP5/WP6 not merged). Recommendation: default ON with
`gap = 0.05`, and WP5 confirms with `seconds.ccd / seconds.total` per step (StepReport); turn it
off by default only if the median there exceeds 10 %.

## Remaining
- `AdvancementPairs.lean` is the proof of `rs-ccd-advancement`; the code tags `advance`, `vertex_triangle_advance`,
  `edge_edge_advance`, `first_contact`.
- Not done: obstacle point clouds (only triangle meshes), moving obstacles, a CCD-aware `Session`
  (WP5), the "number of runs where the energy alone lets surfaces cross" study (needs the flow).
