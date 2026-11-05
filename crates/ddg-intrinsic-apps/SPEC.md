# Applications of intrinsic triangulations — contract

The contract of the crate `ddg-intrinsic-apps`, for Sec. 5 of

> Nicholas Sharp, Yousuf Soliman, Keenan Crane. *Navigating Intrinsic Triangulations.*
> ACM Transactions on Graphics 38(4), 2019. <https://doi.org/10.1145/3306346.3322979>
> (`../NavigatingIntrinsicTriangulations/paper.pdf`, `supplemental.pdf`)

The data structure (signposts, flips, Delaunay flipping and refinement, the optimal Delaunay
triangulation, tracing, point queries) is the crate `ddg-intrinsic`; its contract is
`crates/ddg-intrinsic/SPEC.md`, whose decisions run from `N1` to `N29`. This file continues the
numbering: `N30` to `N39` for Sec. 5.1 to 5.3, `N40` to `N49` for Sec. 5.4 and 5.5. Where the
paper is silent or the authors' demo and geometry-central (both MIT, may be read and ported
with the licence notice) do something the paper does not say, the decision is recorded in §6.
Never drop an entry.

Workspace conventions of `CLAUDE.md` apply: `f64`, `u32` stored indices, explicit option
structs, error enums that name the element, a report from every operation that can do a
variable amount of work, nothing repaired or normalized silently, no panics on bad input. `L`
is the positive semidefinite cotan Laplacian; `Δ` the positive semidefinite Laplace–Beltrami
operator.

## 1. Scope

In:

* 5.1 Steiner tree approximation (`steiner`): the MST of the intrinsic Delaunay triangulation,
  shortened by Fermat points inserted as intrinsic vertices; the cut it defines and the check
  that the cut leaves a flat disk; the extrinsic comparison (MST of the input edges).
* 5.2 Finite elements (`fem`): stiffness, Galerkin and lumped mass matrices on an intrinsic
  triangulation; the regular 4-1 subdivision as the extrinsic comparison; extreme eigenvalues
  and condition numbers.
* 5.3 Geodesic distance (`heat`): the heat method on an intrinsic triangulation, factored once;
  values at vertices and at any point of the input mesh. Its ground truth, the exact polyhedral
  distance from vertices (`exact`).
* 5.4 Adaptive mesh refinement (`amr`): residual a posteriori estimator, bulk marking,
  refinement of the marked faces through `IntrinsicTriangulation::refine_with`, for the
  harmonic Green's function, the short-time heat kernel and a harmonic map of a disk.
* 5.5 Tangent vector fields (`vector`): discrete connection from the signposts, finite
  difference and finite element connection Laplacians, harmonic fields with boundary values
  and the vector maximum principle, smoothest `n`-direction fields and their singularities.

Out (open; `ROADMAP.md`): the extrinsic ODT of Fig. 17; cut graphs on surfaces of genus above
zero (5.1); inserting the 5.1 cut into the mesh as seams for Blender; Fig. 16 (the `L²` error
of a Poisson problem); the Blender side of every application.

## 2. Operations and what each must guarantee

* `steiner::steiner_tree(tri, &SteinerOptions) -> SteinerTree` (edges of a spanning tree of all
  vertices on the triangulation it returns, its length, a report with the MST length and the
  Fermat points). With the default candidates the tree is never longer than the MST of the same
  triangulation. `cut_to_disk(tri, edges) -> CutReport` (N33).
* `fem::{stiffness_matrix, galerkin_mass_matrix, lumped_mass, regular_subdivision,
  condition_number, PinnedLaplacian}` (N34).
* `heat::HeatSolver::new(tri, &HeatOptions)`, `distance(&[HeatSource]) -> HeatDistance`,
  `value_at(values, point)`; `heat::geodesic_distance` for one call (N35, N36).
* `exact::exact_distance(mesh, geometry, sources, &ExactOptions) -> ExactDistance` (N37).
* `amr::{solve, estimate, mark, refine_marked, adaptive_solve, l2_error, input_vertex_error,
  map_quality, disk_boundary, side_flux}` (N46 to N49).
* `vector::{connection_matrices, smoothest_field, singularities, harmonic_field,
  check_maximum_principle, transport_angle, face_holonomy, TangentField::input_vectors_3d}`
  (N40 to N45).

Every option struct has a `Default` that is the paper's choice where it has one; every solve
returns a report; bad input is an error naming the element (`SteinerError`, `FemError`,
`HeatError`, `ExactError`, `AmrError`, `VectorError`).

## 3. What must hold (tests)

The tests are listed in the README's *Verification*. In short: exactness on linear functions
and flat domains (stiffness matrix, heat method, exact distance, harmonic fields); invariance of
surface quantities under retriangulation (the exact distance on the input, the iDT and the iDR;
the vector index sum `nχ`); the analytic cases (Fermat point, cube distances, round sphere
eigenvalues `λ ≈ n`, Knöppel's closed forms); agreement of two independent algorithms where
there are two (Lanczos against dense Jacobi; the finite element connection Laplacian against a
quadrature of its model; the exact distance against Dijkstra on a visibility graph in the
plane); and that every adaptive step reduces the error of a problem with a known solution.

## 4. Verification layers

1. Reports and self-checks: `cut_to_disk` (Euler characteristic, defect, consistent layout),
   `HeatReport::inconsistency`, `ExactReport` counters, eigen-residuals, linear-solve residuals,
   `Singularities::max_rounding`.
2. Independent algorithms in the tests (§3).
3. Lean proofs (`lean/ddg/DDG/FEM/Heat.lean`,
   `lean/ddg/DDG/VectorFields/Holonomy.lean`) and cited results
   (`lean/proofs/CITED.md`), tagged `[formula: app-…]`.

## 5. Formula ids

Machine-checked: `app-heat-gradient`, `app-heat-divergence`, `app-fem-galerkin-mass`,
`app-fem-midpoint-subdivision`, `app-vector-face-holonomy`,
`app-vector-fd-connection-laplacian`, `app-amr-map-distortion`. Cited: `app-steiner-fermat-point`,
`app-vector-transport`, `app-vector-fem-connection-laplacian`, `app-vector-singularity-index`,
`app-vector-maximum-principle`, `app-amr-residual-estimator`,
`app-amr-bulk-marking`. Machine-checked as integrals over a triangle
(`lean/proofs/DdgProofs/TriangleIntegral/`): `fem-area-coordinate-rule` (`monomial_integral`),
`app-amr-element-norms` (the element norms of `estimate`), `app-amr-edge-flux` (`side_flux` is
the outward flux, and twin fluxes add to the jump), and the integrals behind
`app-fem-galerkin-mass`. Shared with the workspace: `cotan-fem-stiffness` (proved),
`cotan-laplacian`.

## 6. Decisions and departures

### Sec. 5.1 to 5.3

* **N30** *The Fermat point.* The paper gives the Fermat point of triangle `ijk` as
  "barycentric coordinates `b_i := csc(θ_i + π/3)`". Those are its trilinear coordinates
  (distances to the sides up to a common factor); the barycentric coordinates are
  `b_k ∝ a_k csc(θ_k + π/3)`, `a_k` the side opposite corner `k`. The two agree only on an
  equilateral triangle. Here: the barycentric form, normalized (`steiner::fermat_point`,
  `[formula: app-steiner-fermat-point]`; the test checks that the point sees every side under
  `120°` and that the tree has the length `S² = (a²+b²+c²)/2 + 2√3 A`). A triangle with an
  angle of `120°` or more has its obtuse corner as Fermat point and is not a candidate.
* **N31** *The greedy order and the candidates.* The paper's summary of Smith, Lee and Liebman
  ("enqueue all edges and per-triangle Steiner trees; greedily add the shortest element") would
  never take a triangle: its tree is longer than each of its sides, whose two shorter ones join
  its corners first. The original ranks a triangle by `|Steiner tree| / |MST of its corners|`
  (Zachariasen and Winter 1997, Sec. 1, 6.2); here all candidate triangles come first in that
  order, then all edges by length (Kruskal over both). Candidates (`TriangleCandidates`): three
  different corners, all angles below `120°`, and by default (`MstAdjacent`, as the original)
  two sides in the MST; `All` takes every such triangle. On a surface with boundary the
  vertices of every boundary loop count as joined from the start (the boundary is a cut
  already). The triangulation is not flipped again after the insertions: the tree's edges
  must stay. On a surface of genus `g > 0` a tree does not cut a disk (Euler characteristic
  `1 − 2g` after the cut); a cut graph needs `2g` loops (Erickson and Har-Peled 2004), not done.
* **N32** *A refused insertion.* `insert_vertex` refuses a point too close to a side of its
  face. Then the triangle is dropped from the candidates and the greedy pass and the
  insertions are repeated without it (`SteinerReport::n_refused`, `n_rounds`).
* **N33** *What `cut_to_disk` checks.* Cut along the given edges (and the boundary): the
  number of pieces, the Euler characteristic of the cut surface, its boundary length (twice the
  tree plus the boundary), the largest angle defect of a vertex not on the cut (zero when the
  disk is flat inside), the turning of the boundary (`2π` for a flat disk), and a layout of
  every triangle by unfolding across the uncut edges, with the largest disagreement between two
  places of one vertex of the cut surface. `is_disk`: one piece of Euler characteristic 1.
* **N34** *Finite elements and the condition number.* Stiffness: the cotan matrix from the
  intrinsic lengths (`cotan-fem-stiffness`). Mass: Galerkin (`A/6`, `A/12`;
  `app-fem-galerkin-mass`) or lumped (`A/3` per corner). The paper does not say which matrix
  Fig. 17 conditions; all three of `L x = λ M x` with `M` Galerkin, lumped and the identity are
  reported. `λ_max`: Lanczos with full reorthogonalization in the `M` inner product on
  `M⁻¹ L`. `λ_min` (smallest non-zero): Lanczos on `L⁺ M`, shift-and-invert at zero through a
  Cholesky factorization of `L` with one vertex per component pinned (`PinnedLaplacian`; exact
  for right-hand sides with zero sum per component), the constants of every component projected
  out. The 4-1 subdivision is computed from lengths alone (`app-fem-midpoint-subdivision`): the
  same surface and angles as the input.
* **N35** *The heat method.* Mass lumped; time step `t = m h²`, `h` the mean edge length of
  the triangulation that is solved on (not of the input), `m = 1` (`HeatOptions::time_factor`;
  or `time` directly); right-hand side the indicator of the source vertices, a source point in a
  face spread over its corners by its barycentric coordinates (as geometry-central). A face in
  which the heat has no gradient (it underflowed to a constant) gets `X = 0` and is counted
  (`HeatReport::n_flat_faces`). The Poisson step is solved exactly with the pinned Laplacian
  (N34) instead of geometry-central's shift of `L` by `1e-6 M`; the result is shifted by the
  mean, weighted by the barycentric coordinates of the sources, of the value at each corner of
  a source face minus the distance of the source from that corner (geometry-central's rule; a
  vertex source gets exactly zero). Boundary: Neumann (default), zero Dirichlet on the heat step, or
  the mean of both (Crane et al. 2017, Sec. 3.4); the Poisson step always has natural
  conditions. Under Dirichlet conditions a source with weight at a boundary vertex is refused
  (`HeatError::SourceOnDirichletBoundary`). Formulas `app-heat-gradient`,
  `app-heat-divergence`.
* **N36** *Where values live.* `distance` returns a value per vertex of the triangulation,
  input vertices first with their ids. An inserted vertex has its place on the input mesh in
  `IntrinsicTriangulation::vertex_location`; `HeatSolver::value_at(values, point)` evaluates at
  any point of the input mesh by a point query (paper Alg. 11) and linear interpolation in the
  intrinsic face. Values at inserted vertices agree with `value_at` at their places (test).
* **N37** *The exact distance.* Written for this crate from Mitchell, Mount, Papadimitriou
  (1987), Surazhsky et al. (2005) and Xin and Wang (2009); no code taken. Vertices whose angle
  sum is at least `2π` (saddles, and flat vertices too) and boundary vertices become
  pseudo-sources. Child windows narrower than `1e-12` of their edge are not created (counted);
  end points that rounding moves off the edge are moved back (largest move reported). The
  vertex filter and the "cannot improve" rule never drop a shortest path. Validated on flat
  convex domains (Euclidean distance, `1e-10`), an annulus against Dijkstra on the visibility
  graph (`1e-10`), a cube (`1e-12`), and by invariance: the same distances at the input
  vertices on the input, the exported iDT and the exported iDR (`1e-9` of the largest).
* **N37b** *The exact distance on a Δ-complex* (Intrinsic Simplification S20).
  `exact::exact_distance_complex(&Complex, &[edge lengths], &[ExactSource], &ExactOptions)`
  shares the propagation of N37 through a private adjacency trait, so the `TriMesh` entry point
  is bit-for-bit unchanged (`tests/exact_pins.rs`: hashes of distances and report on seven
  meshes, recorded before the refactoring) and the `Complex` entry point on `from_tri_mesh`
  gives the same bits (`tests/exact_complex.rs`). Lengths must be positive and form triangles
  of positive area (errors `Lengths`, `InvalidPoint`). Sources: `Vertex(v)` or
  `Point { face, bary }` (corners `faces()[f][0..3]`); a point seeds the three corners at
  their in-face distances and one window per edge across it; on an edge it is moved off by a
  relative `1e-13`, at a corner it is the vertex. Validated by: brute-force unfolding with
  Dijkstra over vertices on complexes coarsened by flips and removals (self-edges, doubled
  edges, repeated corners present; `1e-7`), the lattice distance on a flat 5x4 torus
  coarsened to 8, 4, 3, 2 vertices (`1e-9`), and a face-point source against the vertex source
  of `Complex::split_face` (`1e-11`). The heat method has no Δ-complex form, so there is no
  heat comparison on Δ-complexes.
* **N38** *Measurement of 5.3.* Errors at the input vertices against `exact`, three sources
  per mesh (vertices `0`, `n/3`, `2n/3`): mean relative error over vertices with `d* > 0`, mean
  and largest absolute error over the largest distance. The paper's "iDR" is refinement to an
  angle bound; on well-shaped meshes that adds few vertices and does not reduce the mesh width
  `h`, which is what bounds the heat method's error, so a second column "iDR r/2" also bounds
  every circumradius by half the mean edge of the iDT. Finding: on a flat parallelogram of
  obtuse triangles (the sheared square) the iDT barely helps, because boundary edges cannot be
  flipped and a row of slivers stays along the slanted sides; refinement splits the
  encroached boundary edges and fixes it.
* **N39** *Bindings.* Python (`ddg-py`): `geodesic_distance(positions, triangles, sources, *,
  triangulation, time_factor, boundary)` and `steiner_tree(positions, triangles, *, intrinsic,
  steiner_points)`, results per input vertex and as polylines on the surface. Blender: nothing.
  The natural use, the 5.1 cut as seams for a one-chart unwrap, needs its intrinsic edges
  inserted into the mesh first (as `ddg_uv::cutting` does for the cut flow).

### Sec. 5.4 and 5.5

* **N40** *Frames and storage of fields.* An `n`-direction field is stored by its `n`-th
  power `u_v ∈ C`, the angle measured in the tangent frame of `v`: its signpost frame, so an
  input vertex has the frame of the input mesh on every triangulation and a field computed on
  any triangulation is read at the input vertices as it is. The transport along an edge is
  `ρ_ij = (φ_ji + π) − φ_ij` (VHM Eq. 10; `app-vector-transport`), `r_ij = e^{inρ_ij}`. Complex
  matrices are stored in real form (`a + ib` as `[[a, −b], [b, a]]`) so that the workspace's
  real sparse Cholesky applies. `TangentField::input_vectors_3d` puts one of the `n` vectors in
  the plane of the input face its direction points into. *Representations* (Integer
  Coordinates I26): the frames are those of the signposts on either representation; integer
  coordinates compute them (`vector::TangentFrames`, Integer Coordinates SPEC §4.12,
  `ic-input-aligned-frame`), equal to the stored signposts to rounding on the same
  triangulation; the signposts' results are bitwise those of the signpost-only code
  (`tests/vector_pins.rs`).
* **N41** *Boundary frames.* At an interior vertex the signposts are normalized to `2π`. At a
  boundary vertex the tangent angles are **not** normalized (`s = 1`, GODF Eq. 11): angle 0 is
  the boundary halfedge leaving the vertex, `Θ` the boundary edge arriving. Normalizing by
  `π/Θ` (as the signposts of `ddg-intrinsic`, N2 there) would put curvature into every face at
  the boundary of a flat domain, where a parallel field must have zero energy. The face
  holonomy is `Ω_f = Σ_c (s_c − 1) θ_c` (`app-vector-face-holonomy`).
* **N42** *Sign of the finite difference Laplacian.* `L_ii = Σ_j w_ij`, `L_ij = −w_ij r_ji`,
  so `u* L u = Σ_edges w_ij |u_j − r_ij u_i|²` (`app-vector-fd-connection-laplacian`) and
  `(L u)_i = Σ_j w_ij (u_i − r_ji u_j)`; lumped mass. Positive semidefinite exactly when the
  cotan weights allow it (all non-negative on an iDT); `smoothest_field` refuses an indefinite
  one (`VectorError::NotPositiveDefinite`). Edges between the same two vertices add up; a
  self-edge adds `2w(1 − cos nρ)` to its diagonal.
* **N43** *The finite element connection Laplacian.* Knöppel et al.'s closed forms (GODF Eq.
  17, App. D.3) are partly garbled in the text layer; the model is derived here
  (`vector/knoeppel.rs` module docs): basis sections `ψ_j = b_j ξ_j`, `ξ_j` parallel along
  rays from corner `j`, curvature `Ω/A` constant on the triangle, and both matrices integrated
  exactly monomial by monomial. `Ω` is the unwrapped face holonomy of N41 (not reduced to
  `(−π, π]`): the factors `e^{iΩ b_i}` then see the curvature of the face itself, and the
  indices of `singularities` add up to `nχ` because `Σ_f Ω_f = 2πχ` on a closed surface. Tests: Knöppel's
  mass entry and diagonal, the flat limit (the cotan matrix), and all nine entries against a
  quadrature of the model.
* **N44** *The eigen solver.* The smallest eigenpair of `A u = λ M u` by block inverse
  iteration with Rayleigh–Ritz (GODF Alg. 2 is the block of one), a Cholesky factorization of
  the real form of `A + εM`, `ε = 1e-10 · tr A / tr M` (GODF Alg. 1's shift: makes `A` definite
  when a parallel field exists, does not change the eigenvectors). Block of 4 by default: on a
  round sphere the smallest eigenvalues of 4-direction fields nearly coincide and one vector
  converges slowly. Stop at a relative eigen-residual of `1e-10` or after 1000 iterations
  (reported, not an error). The energy reported is `λ · area` (independent of scale).
* **N45** *Harmonic fields and the maximum principle.* `(L∇ u)_i = 0` at interior vertices,
  `u` prescribed at boundary vertices: the unit tangent, the unit inward normal, or values at
  the input boundary vertices; a vertex inserted on the boundary gets the linear interpolation
  (by length along the boundary) of its neighbouring input values, both transported to it
  along the boundary. The interior block is factored by Cholesky, or LU when it is not positive
  definite (reported). `check_maximum_principle` measures how far each interior value lies
  inside the convex cone of its transported neighbours: margin `π − (a⁺ − a⁻)` (radians);
  outside when below `−tolerance` (`app-vector-maximum-principle`).
* **N46** *The model problems of 5.4.* All of the form `(c + Δ) u = f` with piecewise linear
  elements (cotan stiffness, Galerkin mass): the Green's function of an input vertex,
  `Δu = δ_s − 1/|M|`, `∫u = 0`, natural boundary conditions; the short-time heat kernel as one
  backward Euler step `(1 + tΔ) u = δ_s` (the study and the Python binding take `t` the
  squared mean input edge length); a harmonic map of a disk onto the unit disk (Fig. 20), the boundary by arc length
  from the input boundary vertex with the smallest id. The paper does not say how its kernels
  are discretized beyond "FEM"; these are the standard choices.
* **N47** *The estimator.* `η_T² = h_T² ‖R‖²_{L²(T)} + Σ_e ℓ_e ‖J_e‖²_{L²(e)}`,
  `osc_T = h_T ‖R − R̄_T‖_{L²(T)}`, `R = f − c u` (Mekchay–Nochetto Eqs. 1.3.3 to 1.3.5, 1.3.8;
  `app-amr-residual-estimator`). The supplement prints the area `A_T` where the papers have
  `h_T`: that has the wrong scaling (`A² ~ h⁴`) and would underweight the element residual on
  small faces; `h_T` is the diameter (longest side). The point source is a vertex load of the
  Galerkin system and is not in `R`, which is in `L²`; on the boundary the natural conditions
  contribute the normal derivative itself, Dirichlet conditions nothing. The jump is the flux
  difference through a side (`app-amr-edge-flux`).
* **N48** *Marking and refinement.* The supplement's criterion as printed:
  `Σ_marked η_T > θ Σ η_T` and the same for `osc`, `θ = 1/2` (`BulkSum::Indicators`, default),
  or Dörfler's on squares (`BulkSum::Squares`, as the papers state it): the fewest faces of
  largest indicator, then faces of largest oscillation (`app-amr-bulk-marking`). The marked
  faces are refined through the new hook of `ddg-intrinsic`, `refine_with(options, extra)`: a
  face counts as marked while it is unchanged since the estimate (same id, same three edges),
  so the criterion is false everywhere after finitely many insertions; its circumcentre is
  inserted by the same rules as for a face below the angle bound. `refine` is `refine_with`
  with a criterion that is always false, bit for bit. `adaptive_solve` first refines to the
  angle bound, then repeats solve, estimate, mark, refine until a tolerance, a vertex budget,
  an iteration limit, or a step that inserts nothing (`AmrStatus`).
* **N49** *Measurement of 5.4.* The paper compares adaptive and uniform refinement "at the same
  accuracy" (Fig. 19) without saying against what. Here: a reference solution with several
  times the vertices of any run (an adaptive refinement, then refined further with a
  circumradius bound), checked against a second reference of half the size; the relative `L²`
  error of a run by the reference's lumped quadrature with point queries (`l2_error`), and at
  the input vertices (`input_vertex_error`, the source left out for the Green's function).
  The uniform runs are iDR with decreasing circumradius bounds; for every adaptive step the
  uniform vertex count and time at the same error are interpolated log–log between uniform
  runs. All on one thread; timings from runs made one at a time. The disk map is judged by
  flipped faces and the conformal distortion `σ₁/σ₂` (`map_quality`,
  `app-amr-map-distortion`).
