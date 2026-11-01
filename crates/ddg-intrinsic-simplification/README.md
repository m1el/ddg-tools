# ddg-intrinsic-simplification

Surface simplification using intrinsic error metrics:

> Hsueh-Ti Derek Liu, Mark Gillespie, Benjamin Chislett, Nicholas Sharp, Alec Jacobson, Keenan Crane.
> *Surface Simplification using Intrinsic Error Metrics.* ACM Transactions on Graphics 42(4), 2023.
> <https://doi.org/10.1145/3592403>

The authors' reference code (<https://github.com/HTDerekLiu/intrinsic-simplification>, MIT,
Copyright (c) 2023 Hsueh-Ti Derek Liu) is used as a reference; its notice is in `LICENSES/`.

Work in progress. The contract is [`SPEC.md`](SPEC.md).

## The coarse mesh and vertex removal (WP0–WP2)

* `CoarseMesh` (SPEC §2): a Δ-complex (`ddg_intrinsic::Complex`) with one length per edge, the
  angular coordinate of every halfedge with a reference halfedge per vertex (S4), the tracked
  points (every input vertex, plus any `InputPoint`s given at construction) with a per-face index,
  fixed vertices, and stable vertex ids: every public vertex id is the input vertex id.
  `CoarseMesh::from_positions(positions, faces, points, options)` validates the input and names
  the bad element; `certify(delaunay)` checks the invariants of SPEC §2.4.
* `flip_edge`, `flip_to_delaunay`, `flip_to_delaunay_from`: intrinsic flips with point tracking
  through the planar layout of the diamond, and greedy Delaunay flipping (SPEC §3.2, §3.3).
* `remove_vertex(v)`: flattening (bracketed Newton on Springborn's extended angles, with his flips
  as the fallback), flips to degree three (two on the boundary; an ear's opposite edge first),
  removal with the points of the star re-expressed, Delaunay flips among the neighbours
  (SPEC §3.4, §3.5). In two phases for the error metric: `begin_removal` returns a
  `PreparedRemoval` (neighbours, their curvature before and after, edges and frames after
  flattening, every flip), then `finish_removal` or `revert_removal` (bitwise);
  `tentative_removal` is the same code without point tracking, reverted (SPEC §3.6, S5).
* `Complex::remove_boundary_vertex` in `ddg-intrinsic` (IT-SPEC N54) merges the two faces of a
  boundary vertex.

## The ICE metric and intrinsic coarsening (WP3)

* `metric` (SPEC §4): channels of mass (`ChannelWeights::curvature()` by default,
  `::mixture(w)`, area, per-vertex user arrays split by sign) with weights; `ChannelState` holds
  the masses and error vectors per vertex and evaluates or commits a removal's cost from its
  `PreparedRemoval` (Eq. 9 weights, Eq. 5/8 transport or the memoryless Eq. 6).
* `simplify(&mut mesh, &SimplifyOptions) -> SimplifyReport` (paper Alg. 1, SPEC §5): initial
  Delaunay, a queue keyed by `(cost, vertex)` with lazy deletion, targets (`Vertices`, `Ratio`,
  `PerComponent`, `MaxCost`, `AllRemovable`), fixed vertices, snapshots at given vertex counts
  for a hierarchy, and a report (removals by kind, skips and commit failures by reason, the
  check that a fresh key equals the commit cost, timings, the optional removal sequence).
* `SimplifyOptions::reference_compat(w, target)` reproduces the reference's non-bug choices for
  differential tests (SPEC §12.2, §12.3).
* `cargo run --release -p ddg-intrinsic-simplification --example simplify -- MESH.obj --target N`
  prints the report and quality statistics; `--compat --ref DIR --replay` compares with a run of
  the reference driver.

Tests: `cargo test --release -p ddg-intrinsic-simplification` (certificate after every flip and
removal on generated shapes, Keenan Crane's models and the reference's meshes, the planar-grid
test, curvature conservation, bitwise reverts, ears); `-- --ignored` adds every removable vertex
of the reference's larger meshes. The meshes live in the gitignored `data/` (`data/models/`,
`data/intrinsic-simplification/`, see its `README.md`); missing files are skipped with a notice.
Timing: `cargo run --release -p ddg-intrinsic-simplification --example remove_all -- MESH.obj`
removes every removable vertex by ascending id (no error metric) and prints the time.
