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

## Pre- and post-refinement (WP5; SPEC §6, paper Sec. 6.3)

`refine` module (the workspace's own Delaunay refinement, integer coordinates by default, signposts
as an option; the reference's refinement is slow and has bugs and is not ported):

* `CoarseMesh::from_refined(mesh, geometry, points, mesh_options, &RefinementOptions)`:
  pre-refinement (Fig. 12). Refines the input, makes the coarse mesh of the refined `T¹` with the
  input vertices and `points` tracked (located by the `T⁰ → T¹` query). Inserted vertices have fine
  ids `>= |V|`, are ordinary removable vertices (loop over `CoarseMesh::vertices()` to remove them
  too), are not rows of `P`, and are columns if they stay. `CoarseMesh::from_intrinsic(&T, ...)` is the
  same from a triangulation you refined yourself.
* `CoarseMesh::post_refine(&options)`: refines a coarse mesh to a minimum angle (default 25°) as a
  new intrinsic triangulation and composes the maps (Fig. 13); the result is a `CoarseMesh` with the
  same tracked points, rebuilt frames and fresh ids for the new vertices. The report has the minimum
  angle reached and `guarantee_applies` (closed and all input angle sums >= 60°, Gillespie et al.
  2021a).
* A coarse mesh must be simplicial to be refined (the refinement's input is a `TriMesh`): a
  self-edge, a bigon or two faces on the same three vertices is refused with a message naming them
  (`RefineError::NotSimplicial`); `RefinementOptions::make_simplicial_flips` tries tracked edge flips
  first (`CoarseMesh::make_simplicial`).

Measured with `cargo run --release -p ddg-intrinsic-simplification --example refine_study -- MESH
--points [--signposts] [--repair]` (single thread; random removal order down to the fraction, no error
metric; other agents' load 3 to 5 on 64 cores):

| mesh | pre-refinement | post-refinement of 50 / 25 / 10 / 5 / 2 % of the pre-refined vertices |
|---|---|---|
| blub (7106 v) | +1771 v in 0.029 s (reference: 33 s) | +1278 / +1360 / +863 / +497 / refused (4-face edge), repaired by 1 flip |
| spot (2930 v) | +255 v in 0.008 s | +428 / +519 / +402 / refused (duplicate faces), repaired by 1 flip / +106 |
| bob (5344 v) | +93 v in 0.012 s | all five fractions accepted |
| capsule (10658 v) | +38 v in 0.022 s | all five accepted |
| boxpart (451 v) | +1029 v in 0.009 s | refused at all five without repair; with repair: 50 % (3 flips), 25 % (11 flips) accepted, 10 / 5 / 2 % still refused |

After every post-refinement the minimum angle is 25.00 to 25.10 degrees on these (closed, all angle
sums >= 60°). Signposts give the same vertex counts and angles; their point query fails for points
exactly on a boundary edge of a coarse mesh (see SPEC §6.2, "As built").

Tests: `cargo test --release -p ddg-intrinsic-simplification --test refine`: a planar grid through
pre-refinement, removal of all but a scattering of vertices and post-refinement (every vertex and
extra point reproduced to 1e-13, 3 grids x 2 representations), `P` rows summing to 1 after each
step, certificates, small shapes and the models (spot, bob, blub, boxpart, capsule, ...), the
refusals, the flip repair, both rules on a closed torus.

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

## Hierarchy, multigrid and conformalized mean curvature flow (WP7; SPEC §9, paper Sec. 8.7, Fig. 29-31)

`hierarchy::build_hierarchy(&CoarseMesh, &HierarchyOptions)`: one `simplify` run with snapshots at
`n0, n0/4, n0/16, …` down to at most 500 vertices (or the reference's fresh run per level,
`fresh_levels`), giving per level the faces, side lengths, cotan Laplacian and lumped masses
(Δ-complex safe) and the scalar prolongations `P_k` (`n_k x n_{k+1}`, rows sum to 1). The Galerkin
multigrid of Liu et al. 2021 (`multigrid`) runs on it: `Hierarchy::multigrid(A, known, options)`; with
2+2 symmetric Gauss-Seidel sweeps the V-cycle contracts by 0.23-0.30 (spot, bob) and 0.6 (blub,
spot_subdiv), and PCG with the V-cycle as the preconditioner needs 10-22 iterations for `1e-10`
(SPEC §9.3 "As built (WP7)" has the table). `cmcf::CmcfFlow`: the conformalized mean curvature flow of
Kazhdan, Solomon, Ben-Chen 2012, `(M_t + h L0) X_{t+1} = M_t X_t`, by a direct Cholesky solve
(`step_direct`) or PCG + multigrid on the input's intrinsic hierarchy (`step_multigrid`); the two agree to
`1e-11` per step. Against this crate's supernodal Cholesky the multigrid step is not faster up to 187k
vertices (0.3x at 3k, 0.8x at 187k); the paper's 20x is against a slower direct solver. Studies:
`cargo run --release -p ddg-intrinsic-simplification --example multigrid_study -- MESH.obj` and
`--example cmcf -- MESH.obj --steps 100 --out runs_tmp/cmcf_x` (OBJ per step); tests:
`tests/hierarchy.rs`.
