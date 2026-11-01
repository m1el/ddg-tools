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

## Applications on the coarse mesh (WP6; SPEC §10.1)

`applications::CoarseApplications::new(&CoarseMesh)` compacts the mesh and builds `P`, the coarse
cotan Laplacian `L̄` and lumped mass `M̄` (from the lengths, per halfedge, safe for Δ-complexes:
`CoarseOutput::cotan_laplacian` / `lumped_mass`). On it:

* **Poisson** (Fig. 21): `PoissonSolver` solves `L x = b` with Dirichlet values or, on a closed
  surface, the mean-zero problem (compatible right-hand side, `M`-mean-zero solution); it takes any
  `L`, `M`, so the fine problem uses the same code. `poisson(f̄, fixed, values)` solves
  `L̄ x̄ = M̄ f̄` and returns `x = P x̄`. Fine data are restricted pointwise (`sample`) or as a load
  (`restrict_load`, `b̄ = Pᵀ M f`).
* **Geodesic distance** (Fig. 22, 23): the source (a tracked point: a coarse vertex or a point of
  a coarse face) goes to `ddg_intrinsic_apps::exact::exact_distance_complex` on the coarse
  Δ-complex (`geodesic_exact`), or through `CoarseHeat`, the heat method reimplemented for a
  Δ-complex with lengths (the heat solver of `ddg-intrinsic-apps` needs positions and an
  `IntrinsicTriangulation`); the result is prolongated with `P`.
* **All pairs** (Fig. 24): `all_pairs(pool, ..)` is `D̄` (one exact run per coarse vertex, over a
  `ThreadPool`, bitwise independent of the thread count); `prolongate_rows` / `prolongate_dense`
  evaluate `D̃ = P D̄ Pᵀ`.
* **Voronoi** (Fig. 25): `voronoi(sources, pool, ..)` labels the fine vertices by the nearest source
  in the prolongated exact distances (and the coarse vertices by the coarse ones).
* **Exponential map** (Sec. 8.5.4): `exp_map_face` / `exp_map_vertex` trace the straightest geodesic
  on the coarse mesh by unfolding faces; at a vertex the ray continues by the equal-angle rule. The
  result is a coarse point (the map is fine → coarse, S21), reversible.
* `greedy_coarsening(&mut CoarseMesh, target)` is a stand-in removal order (smallest memoryless
  cost) for tests and the example until `simplify` is merged; **numbers measured with it are not
  those of the ICE order**.

Tests: `cargo test --release -p ddg-intrinsic-simplification --test applications`. Example (errors
and speedups against the fine solutions, one row per coarse size):
`cargo run --release -p ddg-intrinsic-simplification --example applications -- MESH.obj
[--fractions 0.5,0.1,0.02]`; to rerun with ICE, replace `make_coarse` in the example.
