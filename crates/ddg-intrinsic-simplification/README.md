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

## The map and the prolongation (WP4; SPEC §7)

* `CoarseMesh::compact()` returns a `CoarseOutput`: the coarse faces in compacted vertex ids
  (with the map to the fine ids), the side lengths and corner angles of every face, the angular
  coordinates, and every tracked point as `CoarsePoint::Vertex(v̄)` or `Face { face, bary }`.
* `scalar_prolongation()` is `P` (`|V| × |V̄|`, rows sum to one); `point_prolongation()` the same
  for every tracked point (texels, centroids).
* `vector_prolongation(&FineFrames, n)` is `P^vec` (complex, `ComplexSparse`) with a
  `VectorProlongationReport` of how each row was found. It is the paper's construction with
  the unnormalization (`unnormalize_field` / `normalize_field`, `unnormalize` / `normalize`),
  the coarse interpolation `β`, `w`, `z` and the fine sampling `ψ`, with fallbacks that never
  leave a row empty (same face, edge-adjacent face, average of the neighbours' rows). `FineFrames`
  holds the fine mesh's normalized frames (those of `IntrinsicTriangulation::new` on the input);
  `n` takes `n`-direction fields. Note that the vector code of `ddg-intrinsic-apps` stores the
  vectors of a boundary vertex in the unnormalized frame, `P^vec` takes and gives the normalized
  one (`π/Θ`); `unnormalize`/`normalize` convert.
* `colour_faces()` is a deterministic greedy colouring of the coarse faces, distinct across faces
  that share a vertex, with no palette limit; `face_of_points()` and `point_colours()` give the
  coarse face and colour of every tracked point (paper Sec. 8.1).

Tests: `cargo test --release -p ddg-intrinsic-simplification --test prolongation` (flat grids:
linear functions and constant vector/line/cross fields reproduced to 1e-13; row sums; rotation
consistency; no empty rows on the models; the smoothest field of `ddg-intrinsic-apps` on the
coarse mesh against the fine one; the colouring).

## Anisotropic coarsening (WP8; SPEC §8, paper Sec. 8.6.2, Fig. 27)

`anisotropy`: edge lengths scaled by `(1 - tau) + (tau/2)((u_i.e_ij)^2 + (u_j.e_ji)^2)` for per-vertex
direction lines `u` (angles in the normalized frames, from `curvature::direction_angles` or
`principal_direction_field(positions, faces, Max | Min, options)`), the strict triangle inequality
checked per face (an error listing the faces, or `clamp_tau`: the exact largest valid `tau`,
reported), `anisotropic_coarse_mesh` (= `CoarseMesh::from_lengths` of the scaled lengths, tracked
points still on the input faces), and `coarse_face_alignment`, the measurement of Fig. 27. The
coarse triangles are isotropic in the scaled metric, hence **elongated across `u`**: `u = d_max`
gives triangles long along the minimum direction (a cylinder's axis). Study:
`cargo run --release -p ddg-intrinsic-simplification --example anisotropy_study -- MESH.obj
max|min 300 0 0.2 0.4 0.6`; tests: `tests/anisotropy.rs`. Numbers: SPEC §8, "As built (WP8)".
