# Intrinsic Simplification — contract

The contract between the Rust crate `crates/ddg-intrinsic-simplification`, the Lean oracle
(`lean/oracle`), the Lean proofs (`lean/proofs`), the Python binding (`crates/ddg-py`) and the
Blender add-on for

> Hsueh-Ti Derek Liu, Mark Gillespie, Benjamin Chislett, Nicholas Sharp, Alec Jacobson, Keenan
> Crane. *Surface Simplification using Intrinsic Error Metrics.* ACM Transactions on Graphics
> 42(4), Article 118, 2023. <https://doi.org/10.1145/3592403>
> (`../IntrinsicErrorMetric/IntrinsicErrorMetric.pdf`, "the paper", 16 pages; equation, figure
> and algorithm numbers are the paper's).

**Reference.** The authors' code, <https://github.com/HTDerekLiu/intrinsic-simplification> at
commit `708bc7c` (MIT, Copyright (c) 2023 Hsueh-Ti Derek Liu), cloned at
`../IntrinsicErrorMetric/ref/intrinsic-simplification/`, cited as `ref: <file>:<line>` (relative to
its `src/` unless an example is named). It may be read and ported; ported code keeps its notice in
`LICENSES/`, and the crate's README credits all six authors. Its bundled libigl (MPL-2.0) is not
copied from, and **`externals/libigl/include/igl/copyleft/` (GPL) is never opened** (the reference
does not use it). The map of the reference is `runs_tmp/iem/reference-map.md`; its differences from
the paper are numbered `R1…R54` there and are decided below (§14). A headless driver of the
reference (`../IntrinsicErrorMetric/ref-driver/`, built in `../IntrinsicErrorMetric/ref-build/`,
documented in `runs_tmp/iem/reference-driver.md`) dumps coarse meshes, tracked barycentrics, the
removal order and the prolongations for differential tests.

**Papers it builds on** (PDFs in `../IntrinsicErrorMetric/papers/`): Springborn, Schröder, Pinkall
2008, *Conformal Equivalence of Triangle Meshes* (CETM: Eq. 3, 4 and 10 of the paper); Liu, Zhang,
Ben-Chen, Jacobson 2021, *Surface Multigrid via Intrinsic Prolongation* (§9); Kazhdan, Solomon,
Ben-Chen 2012, *Can Mean-Curvature Flow be Modified to be Non-singular?* (§9.3); Garland, Heckbert
1997, *Surface Simplification Using Quadric Error Metrics* (§10.4); Panozzo, Puppo, Rocca 2010,
*Efficient Multi-scale Curvature and Crease Estimation* (§8.6); Khodakovsky, Litke, Schröder 2003,
*Globally Smooth Parameterizations with Low Distortion* (§10.3). Gillespie, Sharp, Crane 2021
(integer coordinates, `crates/ddg-integer-coords`, "IC-SPEC") and Sharp, Soliman, Crane 2019
(signposts, `crates/ddg-intrinsic`, "IT-SPEC") provide the connectivity, Delaunay flipping,
refinement and the input-surface point queries this crate uses.

Conventions of `CLAUDE.md` apply: `f64` lengths, stored indices `u32`, explicit option structs, a
report from every operation that does a variable amount of work, bad input is an error naming the
element and never a panic, nothing repaired or normalized silently, deterministic results at any
thread count, formula tags `[formula: iem-…]` (§13). Decisions of this file have ids `S<n>`.

---------------------------------------------------------------------------------------------------

## 1. The problem

Input: a manifold, oriented triangle mesh `M = (V, E, F)` with edge lengths `ℓ` (from vertex
positions, or given), possibly with boundary and several components. Output: a coarse intrinsic
triangulation `M̄ = (V̄, Ē, F̄)` with edge lengths `ℓ̄`, `V̄ ⊆ V` (plus inserted vertices when
pre-refinement is on, §6), together with a **bijective map** `φ: M → M̄` that sends every fine
vertex (and any other tracked point) to a coarse face and barycentric coordinates (§7.1), the
**prolongation** operators `P` (scalars) and `P^vec` (tangent vectors) from `M̄` to `M` (§7.2,
§7.3), and a report.

Unlike the intrinsic triangulations of IT-SPEC and IC-SPEC, the coarse mesh is **not** a
retriangulation of the input metric: every removal first *flattens* the vertex (§3.4), which
changes the lengths of its star conformally and moves its curvature to its neighbours. The coarse
metric approximates the fine one; the ICE error metric (§4) decides which vertex to remove next so
that the curvature (and other "mass") moves as little as possible. Hence this crate has its own
data structure (S1) and does not implement `ddg_intrinsic::Triangulation` (which promises a fixed
input metric).

Non-goals: extrinsic positions of the coarse mesh (it has none; it is drawn on the fine mesh,
§7.4); topology changes (a vertex whose removal would change the topology is never removed, §3.5);
the "replay" of operations for points not known in advance (paper Sec. 7.1, last sentence: all
points to track are given before simplification; S17).

---------------------------------------------------------------------------------------------------

## 2. Data model

### 2.1 Connectivity

`ddg_intrinsic::Complex` (IT-SPEC N50): a halfedge Δ-complex with its conventions (halfedge
`3f + k` from `faces[f][k]` to `faces[f][(k+1)%3]`, one halfedge per boundary edge, CCW fans), its
operations (`flip_edge`, `remove_vertex`, `merge_edge`, `split_face`, `split_edge`) and the id moves
they report. A boundary vertex of degree 2 is removed by a new operation
`Complex::remove_boundary_vertex` (WP0, IT-SPEC N54): for a boundary vertex `i` with boundary
neighbours `b` (before it on the boundary) and `c` (after it) and one interior spoke `ia`, the two
faces `iab`, `ica` become one face `abc` (paper Sec. 4.2, inset right): the spoke `ia` is deleted,
the two boundary edges `bi`, `ic` become one boundary edge `bc` (its length, `ℓ_bi + ℓ_ic` for a
flat `i`, is the caller's business), and the outer edges `ab`, `ca` are kept. Id rules as
`remove_vertex`: the merged face takes the smaller face id, freed ids are filled by the last ones,
every move is reported.

**Stable vertex ids (S2).** `Complex` renumbers vertices on removal (the last vertex takes the
freed id). The crate hides this: every public vertex id is the **fine id** (`0..|V|`, then
`|V|..` for inserted vertices in creation order), and an internal map `complex ↔ fine` is updated
from each removal's report. Edge and face ids are internal; the output mesh is compacted
(`CoarseMesh::faces()` lists fine ids).

### 2.2 Stored quantities

* `ℓ_e > 0` per edge.
* **Angular coordinates** `φ_h` per halfedge, at its tail `i`: the angle of the halfedge from the
  vertex's reference direction, normalized (paper Sec. 3.2.2): `φ = s_i θ`, `θ` the unnormalized
  angle measured CCW, `s_i = 2π/Θ_i` at an interior vertex and `s_i = π/Θ_i` at a boundary vertex
  (S4; IT-SPEC N41, the signposts' convention, so that `φ ∈ [0, 2π)` inside and `φ ∈ [0, π]` on the
  boundary with the boundary halfedge at `0` and the incoming boundary edge's virtual direction at
  `π`). Each vertex has a **reference halfedge** `r_i` (any outgoing halfedge; on the boundary the
  outgoing boundary halfedge) whose `φ` is held fixed when the lengths of the star change (S4).
* **Channels** `p = 0..k` (§4.1): per vertex a mass `m_i^p ≥ 0` and an error vector
  `t_i^p ∈ ℂ` in the normalized frame of `i`; weights `w^p ≥ 0`.
* **Tracked points**: per point `b` either `Vertex(v)` (a fine vertex that is still a vertex) or
  `Face(f, [b₀, b₁, b₂])` with barycentric coordinates w.r.t. the corners of face `f`
  (`faces[f][k]`), `b_k ≥ 0`, `Σ b_k = 1`; and per face the list of points in it (an index kept
  in sync by every operation). Points `0..|V|` are the fine vertices; further points (texels,
  sample points, fine face centroids for the Blender view) are added by the caller before
  simplification (S17).
* **Fixed vertices** (paper Sec. 8.6.3): a mask; a fixed vertex is never removed (cost `+∞`).

### 2.3 Derived quantities

Corner angles from the lengths by the law of cosines (`ddg_core` / `ddg_intrinsic::formulas`),
angle sums `Θ_i`, `K_i = 2π − Θ_i` (interior, Eq. 1), `κ_i = π − Θ_i` (boundary, Eq. 2), face areas
(Heron, stable form), the cotan Laplacian and lumped/Galerkin mass matrices (the formulas of
IT-SPEC / `ddg_core`, here over a Δ-complex: per halfedge, so repeated vertices and self-edges are
handled as in `ddg_intrinsic::generic_ops`).

*As built (WP1).* The angle sum `Θ_i` is cached per vertex and refreshed (from the lengths,
summed counter-clockwise from the reference halfedge) whenever the lengths at `i` change
(flattening) or its corners change (the corners of a merged face after a removal); flips keep it
(they do not change the metric). `K_i`, `κ_i` and the scale `s_i` of the angular coordinates come
from the cache, so every operation sees bitwise the same values.

### 2.4 Validity (`certify()`)

1. `Complex::check()`; the fine/complex vertex maps are inverse bijections on live vertices.
2. Every face satisfies the strict triangle inequality.
3. Angular coordinates: per vertex, the `φ` of the fan are increasing CCW from the reference
   modulo `2π` (inside) / `[0, π]` (boundary), consecutive differences equal `s_i ×` the corner
   angles to `1e-9` rad, the boundary halfedge has `φ = 0`.
4. Total curvature: `Σ K + Σ κ = 2πχ` to `1e-9 · (1 + Σ|K| + Σ|κ|)` (Gauss–Bonnet; flattening
   preserves it, `iem-flatten-conserves-curvature`).
5. Every tracked point is in exactly one place; `Face` points have `b_k ≥ 0`, `|Σ b − 1| ≤ 1e-12`,
   and appear in their face's list; `Vertex(v)` points name a live vertex.
6. Masses `m ≥ 0`, all values finite.
7. After `simplify` (unless Delaunay is off): the intrinsic Delaunay condition at every interior
   edge (the test of §3.3).

`certify()` is cheap enough to run after every operation in tests (`O(|F| + #points)`).
*As built:* `certify(delaunay: bool)` checks item 7 only when asked and always measures it
(`Certificate::{max_delaunay_excess, n_not_delaunay}`); item 3 recomputes the angle sums from the
lengths (not the cache); item 6 covers lengths, angular coordinates and point coordinates (masses
come with WP3; they live in the run's `ChannelState`, not in `CoarseMesh`, and are checked
by the tests of `simplify`).

---------------------------------------------------------------------------------------------------

## 3. Operations

### 3.1 Construction

* `from_mesh(mesh, positions | lengths, options)`: validates (manifold, oriented, finite positive
  lengths, strict triangle inequality: errors name the face), builds the `Complex`, the angular
  coordinates (reference = `Complex::outgoing_halfedge`), the channels (§4.1), the tracked fine
  vertices and any extra points (`Face(f, b)` on the input mesh).
* `from_intrinsic(triangulation, points)` (pre-refinement, §6.1): the `Complex` and lengths of an
  integer-coordinate (default, Q3 of IC-SPEC) or signpost triangulation `T¹` of the input, the
  fine vertices as `Vertex` points, other input-surface points located in `T¹` by the
  representation's `T⁰ → T¹` point query.

### 3.2 Edge flip

An intrinsic flip (paper Sec. 3.3.1): flippable iff both endpoints have degree ≥ 2 and the two
triangles form a strictly convex quadrilateral (`Complex::flip_refusal` plus the convexity test of
`ddg_intrinsic::formulas::unfold_diamond`). New length from the planar unfolding
(`it-flip-length`). Angular coordinates: the four outer halfedges keep their `φ`; each new
halfedge gets `φ` from its CW neighbour plus `s ×` the corner angle between them
(`it-signpost-update`, the signposts' rule), wrapped (S4). If a flipped halfedge was a reference,
the reference passes to its CW neighbour (whose `φ` is kept). Tracked points of the two faces are
re-expressed in the new faces through the planar layout of the diamond (S16).

*As built (WP1).* Flippable is `ddg_intrinsic::plan_flip_with` with the option `flip_test`,
default IT-SPEC N5's shape test (strict convexity plus the shape limit `1e-6`, so that no flip
creates a triangle its lengths do not describe). Point relocation: the diamond is laid out by
`it-unfold-diamond`; the end point of the old edge that a new face lacks is expressed in that
face (`it-barycentric-coordinates`), and a point's coordinates in both new faces follow by linear
substitution. The inside test of S16 uses the signed distances `b_k h_k` to the sides (`h_k` the
heights of the new face from its lengths) against `−point_slack · L`, `L` the mean of the five
diamond lengths, `point_slack = 1e-12`.

### 3.3 Delaunay flipping

The intrinsic Delaunay test of IT-SPEC (the opposite-angle sum `θ^k_ij + θ^l_ji ≤ π`, evaluated as
`ddg_intrinsic::formulas::delaunay_angle_excess` with IT-SPEC's tolerance), greedy flipping from a
queue (Bobenko–Springborn): `flip_to_delaunay` from all interior edges at construction (S3), and
after each removal from the edges whose two endpoints are in `N_i` (paper Sec. 6.2, S15), pushing
the four outer edges of every flip. Edges that fail the test but are not flippable are skipped and
counted.

### 3.4 Flattening (paper Sec. 4.1, Eq. 3–4)

Find `u = u_i` such that scaling every edge of the star of `i` by Eq. 3 with `u_j = 0` for all
`j ≠ i` (`ℓ̃_e = e^{(u_i + u_j)/2} ℓ_e`, i.e. `e^{u/2}` per endpoint equal to `i`: a self-edge `ii`
gets `e^{u}`, outer edges are unchanged) makes `Θ̃_i = 2π` (interior) or `π` (boundary).

**Solver (S9).** With the extended angles of Springborn et al. 2008 Sec. 3.1 (a face violating the
triangle inequality has angles `π, 0, 0`), `Θ_i(u)` is continuous and non-increasing in `u`, with
`Θ_i(u) → deg(i)·π` (`≥ 2π`; boundary: `≥ π`) as `u → −∞` and `→ 0` as `u → +∞`
(`iem-flatten-monotone`), so a root exists. Newton's step is Eq. 4,
`u ← u − (Θ̃_i − Θ_i)/(½ Σ_{ijk ∋ i} (cot θ^k_ij + cot θ^j_ki))` (sum over the faces of the star,
each face once per corner at `i`; in a face with a self-edge `ii` the cotangent opposite the
self-edge is dropped, since that corner's angle does not depend on `u` — the per-corner sum of Eq. 4,
which the reference uses (R55), is not the derivative there: a gap of 0.42 on a self-edge star in
the oracle's self-test; only Newton's speed is affected, not the root), safeguarded by a bracket `[u_lo, u_hi]` with
`Θ_i(u_lo) ≥ Θ̃_i ≥ Θ_i(u_hi)` that is expanded geometrically from `u = 0` and replaced by
bisection whenever the Newton step leaves it or the Hessian is not positive (cot terms of
degenerate corners count 0). Stop when `|Θ_i − Θ̃_i| ≤ 1e-12 · Θ̃_i`; iteration cap 100 (the
report counts iterations; the reference's 5-typical holds for Newton). The reference's NaN-halving
(R7, R8) is not used.

**Validity and Springborn's fallback.** At the root the scaled star must satisfy the strict
triangle inequality. If it does not (an angle `π` in some face, `iem-flatten-monotone`'s extended
case), the flip fallback of Springborn et al. 2008 Sec. 3.2: in the *pre-flattening* metric, flip
an edge of the star opposite a corner whose scaled angle is `π` (if flippable; spokes are tried
before outer edges, outer edges are not boundary), and solve again; at most `deg(i) + 3` attempts,
then the flattening fails (the vertex gets cost `∞`, §4.4). The fallback flips are part of the
removal and are reverted with it if the removal fails.

**Commit.** Lengths of the star scaled; `K_j` of the neighbours updated (their angle sums change,
`iem-flatten-conserves-curvature`: `Σ_{j ∈ N_i} K̃_j = Σ_{j ∈ N_i} K_j + K_i`); the angular
coordinates of `i` and of every neighbour recomputed from their references (S4); tracked points
of the star faces moved by the projective interpolation of Springborn et al. 2008 Sec. 3.4,
`b̃ = (e^{−u_i} b_i, b_j, b_k)/(e^{−u_i} b_i + b_j + b_k)` with the factor applied to every corner
equal to `i` (the paper cites this map but prints Eq. 10 with `e^{+u_i}`, its inverse; S25:
`ProjectiveSign::Negative` is the default, `Positive` an option; `iem-projective-interpolation`: either is continuous
across edges and a bijection of each face onto itself).

*As built (WP2).* The bracket starts open on both sides at `u = 0`; while one side is open a
step leaving the bracket moves `1, 2, 4, …` past the known end; `|u| ≤ 60`
(`FLATTEN_MAX_FACTOR`, a factor `e^{30}` on the lengths) else the solve fails. For a Δ-face with
`i` at several slots the derivative is the sum of Springborn's `∂θ_a/∂u_s` over the slots `a, s`
of `i` (`∂θ_a/∂u_a = −½(cot θ_b + cot θ_c)`, `∂θ_a/∂u_b = ½ cot θ_c`), which is Eq. 4 for a
simplicial face. A bracket shrunk to adjacent doubles is accepted if the residual is within
`1000 ×` the tolerance, reported. Fallback candidates: in every scaled face that violates the
triangle inequality, its longest side; spokes first, then non-boundary outer edges, each by edge
id; the first flippable one is flipped; at most `vertex_degree(i) + 3` solves. The scaled
lengths of the solve and of the commit are bitwise equal, so the validity found by the solve
holds after the commit.

### 3.5 Vertex removal (paper Sec. 4.2)

After flattening, the star of `i` is flat (angle sum `2π`, or `π` on the boundary), so the
following steps do not change the geometry.

* **Interior.** Flip spokes until `deg(i) = 3`: first any self-edge at `i` (paper), else the
  flippable spoke `ij` with the largest opposite-angle sum `θ^k_ij + θ^l_ji`; non-flippable spokes
  are skipped in that round and re-scored every round (S14, R23). By Gillespie et al. 2021a
  Thm. D.1 (`ic-vertex-removal`, proved, tier B′) a flat simplicial star of degree > 3 always has
  a flippable spoke; in a Δ-complex it may not (`ic-vertex-removal-delta`, open): if no spoke is
  flippable at degree > 3, the removal fails and is reverted (paper "skip this vertex removal and
  revert", S7; never a panic, R6). Then `Complex::remove_vertex`; `i` becomes a tracked point of
  the merged face with barycentric coordinates from the layout of the three faces
  (`iem-degree-three-barycentric`: the normalized signed sub-triangle areas, exact for a flat
  star), and the points of the three faces are re-expressed in the merged face (linear
  substitution, S16).
* **Regular boundary vertex** (degree ≥ 2, `κ̃ = 0` after flattening): flip spokes (not the two
  boundary edges) until degree 2, same rule; then `Complex::remove_boundary_vertex`; `i` lies on
  the new boundary edge `bc` (paper Sec. 4.2, S13; the reference's flip to degree 1 and ear
  deletion, R9, is not used).
* **Ear** (boundary vertex of degree 1): flip the edge opposite `i` first (paper Sec. 4,
  "Special cases"), which makes `i` a regular boundary vertex of degree 2; if that edge is a
  boundary edge or not flippable, the removal fails.
* **Skipped vertices** (cost `∞`): an endpoint of a boundary self-edge (removal would delete the
  last vertex of a boundary loop), a vertex of a self-face (a face with one distinct vertex),
  a fixed vertex, an interior vertex that cannot reach degree 3 (e.g. degree 2), and a vertex
  whose flattening or flips fail. There is no separate minimum vertex count: the operations never
  change the topology, and a component stops when none of its vertices is removable (S12; high-genus
  components stop early, paper Fig. 17).

Every removal reports what it did (flattening iterations, fallback flips, flips to degree 3,
Delaunay flips after it) and can be undone exactly (S5).

*As built (WP2).* Order: an ear's flip, then the flattening with its fallback flips, then the flips
to degree 3/2, then the checks that the removal can complete (the complex accepts it; the merged
triangle — the outer lengths, or `ℓ_bi + ℓ_ic, ℓ_ca, ℓ_ab` on the boundary — is strictly valid,
else `Degenerate`), all on the undo log. `N_i` is the set of distinct vertices `≠ i` of the star
at flattening (after any fallback flips); its curvature before and after, and for each neighbour
the shortest connecting edge after flattening with `φ_ij`, `φ_ji` (`π` for the virtual direction
of an arriving boundary edge, S4), are in the `PreparedRemoval` for the metric. Flip choice (S14):
scores first, then the flip test in order (self-edges, then the larger opposite-angle sum, then the
lower edge id); the first flippable spoke is flipped; at most `4 vertex_degree(i) + 8` flips, else
`Stuck`. Degree 3 (2) reached with a self-edge at `i` is the refusal `SelfEdge`. Removal from the
complex: `Complex::remove_vertex` (interior) or `Complex::remove_boundary_vertex` (IT-SPEC N54);
a reference halfedge that dies passes to its clockwise neighbour, all per-halfedge, per-edge,
per-face and per-vertex data and the point lists follow the reported moves; the removed fine
vertex becomes a point of the merged face (degree three: `iem-degree-three-barycentric`;
boundary: `(1 − t, t)` on `bc` with `t = ℓ_bi/(ℓ_bi + ℓ_ic)`, `iem-boundary-merge`); the points of
the star are re-expressed by linear substitution (exact for a flat star, no clamping); the angle
sums and angular coordinates of the merged face's corners are refreshed from their references.
Delaunay seeding (S15): the interior edges with both end points in `N_i`.

### 3.6 Tentative removal

The cost of removing `i` (§4) is evaluated by performing flattening, the flips to degree 3/2 and
the mass/vector transport **on an undo log** and reverting (paper Sec. 6.1). The tentative and the
committed removal run the same code, so the queued key equals the cost of the removal if the
neighbourhood is unchanged (S5, R16).

*As built (WP2).* Two phases: `begin_removal(v) -> PreparedRemoval` (everything up to the removal
from the complex, on the undo log) and then `finish_removal()` (the removal, the points, local
Delaunay; not undoable) or `revert_removal()` (bitwise restore, including the order of the
per-face point lists and the clamp counter). `tentative_removal(v)` is the first phase without
point tracking followed by the revert: points do not enter the cost, so the `PreparedRemoval` is
the committed one's (asserted, `tests/removal.rs`). While a removal is pending every other
mutation is refused (`RemovalPending`).

---------------------------------------------------------------------------------------------------

## 4. The ICE error metric (paper Sec. 5)

### 4.1 Channels (Sec. 5.4, 5.4.1, 8.6.1)

Each channel is a non-negative mass function. Built-in sources:

* `CurvaturePositive`, `CurvatureNegative`: `K⁺ = max(K, 0)`, `K⁻ = −min(K, 0)` at interior
  vertices, `κ⁺`, `κ⁻` at boundary vertices (Sec. 5.4);
* `Area`: the lumped vertex area `⅓ Σ A_f`, divided by the total area of the mesh (as the
  reference, R2; so the channel totals 1);
* `PerVertex(values)`: any user array; a signed array is split into two channels `max(x, 0)`,
  `−min(x, 0)` (Sec. 5.4.1), e.g. the spatial adaptivity masses `1/d_x²(i)` of Sec. 8.6.1.

with weights `w^p ≥ 0`. **Default (S10):** the two curvature channels with weight 1 and nothing
else (the reference's default `w = 0`). `ChannelWeights::mixture(w)` reproduces the reference's
single knob: curvature channels `1 − w`, area `w` (Fig. 11).

### 4.2 Redistribution weights (Eq. 9)

`α_ij = |K̃_j − K_j| / Σ_{l ∈ N_i} |K̃_l − K_l|` with `K` replaced by `κ` at boundary neighbours.
`N_i` are the distinct neighbours `j ≠ i` (a neighbour joined by several edges appears once; S12).
If `|K_i| ≤ 1e-12` (resp. `|κ_i|`) the vertex is already flat and `α_ij = 1/|N_i|` (paper; S11
replaces the reference's exact-zero test and `+1e-6`, R12). Otherwise `Σ |K̃_l − K_l| ≥ |K_i| > 0`
(`iem-flatten-conserves-curvature`), so Eq. 9 is defined. The same `α` is used for every channel
(paper).

### 4.3 Masses and error vectors (Eq. 5, 7, 8)

For each channel and each `j ∈ N_i`:

```
m̃_j = m_j + α_ij m_i                                                   (Eq. 5)
t̃_j = (α_ij m_i (R_ij t_i + e_ji) + m_j t_j) / m̃_j     (m̃_j > 0; else t̃_j = t_j)   (Eq. 8)
```

with `e_ji` the vector of the edge from `j` to `i` in the normalized frame of `j`
(`ℓ_ij e^{i φ_ji}`), `R_ij = e^{i(φ_ji + π − φ_ij)}` the transport from `i` to `j` along that edge
(Sec. 3.2.2), both **after** flattening (the frames of the flattened star). Which edge: the
shortest edge between `i` and `j`, ties by the lower edge id (S12; R13 used the first in CCW
order). The memoryless metric (Sec. 5.1) is the option `memory = false`: `t_i = 0` always and the
cost `Σ α_ij m_i ℓ_ij` (Eq. 6, second form; `iem-memoryless-cost`: it equals the 1-Wasserstein
distance of moving `m_i` from `i` to its neighbours on the flat star). Memory-based is the default
(paper: "we use the memory-based metric for all results").

At a boundary neighbour, the edge `ij` may be a boundary edge; `φ_ji` and `φ_ij` are then `0` or
`π` by S4, and the formulas apply unchanged (the reference's three boundary cases, R14, are not
needed).

### 4.4 The cost (Eq. 6, Sec. 5.4.1)

`C_i = Σ_{j ∈ N_i} Σ_p w^p m̃_j^p ‖t̃_j^p‖`, over **all** neighbours and channels, including those
that receive no mass (their existing `m_j‖t_j‖` counts; paper Eq. 6, S11 against R11).
`C_i = ∞` when `i` is skipped (§3.5) or the tentative removal fails. After a removal the vertex's
own masses and vectors are dropped (they now live at the neighbours).

*As built (WP3, `src/metric.rs`).* `ChannelWeights { channels: Vec<Channel { source, weight }> }`
with `ChannelWeights::curvature()` (the default, S10), `::mixture(w)` and `.with(source, w)`;
sources `CurvaturePositive`, `CurvatureNegative`, `Area`, `PerVertex(values)` (one value per fine
vertex; split into `user<k>+`, `user<k>−` only if a value is negative). The state of a run is a
`ChannelState` (masses and vectors per fine id and channel, outside `CoarseMesh`, so `certify`
does not see them; the tests check `m ≥ 0`, finite, zero at removed vertices, conserved totals).
Channels of weight 0 are not tracked (they cannot change a cost). `fold_weights` multiplies the
weights into the masses (R2; `t̃` is invariant under a common factor, so the cost is the same up to
rounding). The plain formulas are public: `metric::redistribution_weights` (Eq. 9,
`iem-redistribution-weights`), `metric::transport_update` (Eq. 5 and 8, `iem-center-of-mass-update`;
`m̃_j = 0` keeps `t_j`), `metric::memoryless_cost` (Eq. 6 second form, `iem-memoryless-cost`). If
`|K_i| > 1e-12` but `Σ|K̃_l − K_l|` is not a positive finite number (rounding only), the shares are
equal too. A NaN cost is queued as `∞`. `MetricOptions { memory, flat: FlatRule, receiving_only,
boundary_frames }`; `FlatRule::{Threshold(1e-12), ReferenceEpsilon}` (R12),
`BoundaryFrames::{Paper, Reference}` (R3), `Frames::{AfterFlattening, BeforeFlattening,
ReferenceCommit}` (R16, R19) select the angular coordinates of the transport: `NeighbourChange`
carries `φ_ij`, `φ_ji` of the connecting edge both after and before flattening.

---------------------------------------------------------------------------------------------------

## 5. Simplification (paper Sec. 6, Alg. 1)

`simplify(options) -> SimplifyReport`:

1. Flip to intrinsic Delaunay (S3; option `initial_delaunay`, default on; the reference skips it,
   R1).
2. Masses from the channels, `t = 0`, the tentative cost of every vertex, a priority queue keyed
   by `(cost, fine id)` (deterministic ties) with lazy deletion (a version counter per vertex).
3. Loop: pop the minimum valid entry. If its cost is `∞`, stop (paper; R17). Otherwise commit the
   removal (flatten, flip, remove, transport, Delaunay flipping of §3.3). If the commit fails
   (it cannot when the neighbourhood is unchanged since the key was computed, but a stale key may
   meet a changed star), revert, mark the vertex, and do not re-queue it until a vertex of its
   star changes (paper; S7). Re-evaluate the cost of every vertex of `N_i` and every endpoint of
   an edge flipped by the Delaunay pass (S6) and push them.
4. Stop at the target: `TargetVertices(n)` (total live vertices), `Ratio(r)` (of the input,
   rounded), `PerComponent(r)` (each component to `max(min_c, round(r |V_c|))`; Sec. 8.2 normalizes
   by components), `MaxCost(c)`, or when nothing is removable (all `∞`).
5. Optional post-processing: Delaunay refinement (§6.2).

The report: removals by kind (interior, boundary, ear), skips by reason, commit failures, stale
pops, flattening iterations (max, mean), fallback flips, flips to degree 3/2, Delaunay flips,
points clamped (S16), final counts per component, time per phase, and the sequence of
`(vertex, cost)` if asked (for differential tests).

**Snapshots (S24).** `simplify` can record snapshots at given vertex counts (a descending list):
each snapshot is a coarse mesh plus the tracked points of the previous snapshot's vertices,
giving the hierarchy of §9 from one run (the masses and vectors carry across levels).

**Performance.** The paper decimates ≈ 10 000 vertices/s (i7-8750H), the reference ≈ 12–15 000/s
on this machine (single thread). Floor: ≥ the reference's rate on spot, bob, blub, spot_subdiv
with tracking on; near-linear scaling with `|V|` (Fig. 16).

*As built (WP3, `src/simplify.rs`).* `simplify(&mut CoarseMesh, &SimplifyOptions) ->
Result<SimplifyReport, IemError>`. Options: `target: Target::{Vertices(n), Ratio(r),
PerComponent(r), MaxCost(c), AllRemovable}` (`Vertices` is the `TargetVertices` above; `Ratio` and
`PerComponent` round `r ·` the count when `simplify` starts; a popped vertex of a component at its
target is discarded; `MaxCost` stops when the minimum valid key exceeds `c`), `channels`,
`metric`, `initial_delaunay`, `requeue: Requeue::{StarsChanged, OneRing}`, `fixed` (added to the
mesh's mask and left set), `snapshots` (strictly descending counts), `record_sequence`,
`check_keys`, `reference: Option<ReferenceCompat>` (§12.2). Bad options are
`IemError::InvalidOptions` and change nothing. The queue holds `(cost, stamp, fine id, version)`
ordered by `f64::total_cmp`, then stamp (always 0 unless the reference's ties are asked for), then
id; an entry is valid iff its version is the vertex's current one and the vertex is alive. Every
live vertex is evaluated once at the start (ascending fine id). A popped valid vertex whose
removal fails at commit is counted by reason and not re-queued (S7: it comes back when a
re-queue reaches it). **Re-queue (S6, as built):** `N_i` plus the four vertices (old and new end
points) of every flip of the removal — the ear flip, Springborn's fallback flips, the flips to
degree 3/2 and the Delaunay flips after it. **Key check (S5, as built):** each vertex records the
removal count of its last change (it is the removed vertex, in `N_i`, or a vertex of one of those
flips) and of its last key; a key is *clean* if neither the vertex nor any of its neighbours
changed since the key was computed, and a clean key must equal the commit cost bitwise
(`KeyCheck::n_clean_mismatch = 0`, asserted on all test shapes in default and `OneRing` modes;
keys that are not clean differ when a 2-ring vertex's masses changed, the paper's approximation).
Report (`SimplifyReport`): stop reason, counts before/after, initial Delaunay flips and status,
removals by kind, `skips` (live vertices at the end whose last evaluation was `∞`, by reason),
`infinite_evaluations` and `commit_failures` by reason (`refusal_name`), pops, stale pops,
component-done pops, evaluations, flattening steps (max, mean), fallback/degree/Delaunay flips,
refused Delaunay flips, clamped points, per-component counts, `KeyCheck`, `SimplifyTimings`
(initial Delaunay, channels, initial costs, removals, re-queue, snapshots, total), the optional
`(vertex, kind, key, cost)` sequence, the snapshots, the final `ChannelState`.
**Snapshots (S24, as built):** `Snapshot { n_removed, vertices (fine ids; index = compact id),
faces (compact ids), lengths (per face side k → k+1), previous_vertices, previous:
Vec<PointLocation> }`, taken when the vertex count first drops to each requested count (also at
the start if it is already there; counts never reached are not taken). The compact ids are the
complex's (the complex keeps its ids dense), i.e. the ids of WP4's `CoarseMesh::compact()`.
`previous` locates the previous level's vertices (the vertices at the start for the first
snapshot) in this one: `P_k` of §9.1 is built from it.

*Performance (WP3).* `examples/simplify.rs --repeat 3` (best of 3, default options, fine vertices
tracked, `simplify` from the initial Delaunay to the end) against `iem_ref --unlogged`
(`coarsen_mesh` alone, best of 3), both pinned to one core, one run at a time, load average
2.7–2.8 on 64 cores from other work (`runs_tmp/iem/wp3/bench1/`); "compat" is
`reference_compat(0, …)`:

| mesh, target | removals | ours s | ours removals/s | reference s | reference removals/s | ratio | compat ratio |
|---|---|---|---|---|---|---|---|
| spot → 500 | 2430 | 0.168 | 14 500 | 0.190 | 12 800 | 1.14 | 1.17 |
| spot → min | 2927 | 0.190 | 15 400 | 0.219 | 13 400 | 1.15 | 1.19 |
| bob → 500 | 4844 | 0.297 | 16 300 | 0.326 | 14 800 | 1.10 | 1.12 |
| bob → min | 5343 | 0.326 | 16 400 | 0.366 | 14 600 | 1.12 | 1.15 |
| blub → 500 | 6606 | 0.413 | 16 000 | 0.474 | 13 900 | 1.15 | 1.17 |
| blub → min | 7103 | 0.444 | 16 000 | 0.517 | 13 800 | 1.16 | 1.19 |
| spot_subdiv → 500 | 24255 | 1.500 | 16 200 | 1.735 | 14 000 | 1.16 | 1.18 |
| spot_subdiv → min | 24752 | 1.540 | 16 100 | 1.798 | 13 800 | 1.17 | 1.19 |

Two thirds of the time is the re-queue (≈ 6 tentative removals per removal, ≈ 7 µs each) and
13 % the initial costs; the
profile is dominated by corner angles (Kahan's `atan2` form: ≈ 35 %), computed again when a
neighbour's fan is refreshed after the tentative flattening and when the spokes are scored for
the flips to degree 3 (the spoke scores are cached between the rounds of one removal: 2–3 %). A
per-halfedge cache of corner angles (kept through flips, flattening and
the undo log) would remove most of them; it changes the core's cache contract (§2.3) and is left
for later.

---------------------------------------------------------------------------------------------------

## 6. Intrinsic retriangulation (paper Sec. 6.3)

### 6.1 Pre-refinement

Intrinsic Delaunay refinement of the input (IC-SPEC §5 / IT-SPEC refinement, default angle bound
25°, the representation `Representation::DEFAULT` = integers) inserts candidate vertices; the
coarsening then starts from the refined `T¹` (`from_intrinsic`, §3.1). Inserted vertices are
ordinary vertices (they may stay as coarse vertices: `V̄ ⊄ V` then); they are not rows of `P`
(rows are the fine vertices); they are columns if they stay. The reference recipe pairs
pre-refinement with an area weight `w = 0.1` (else all the zero-curvature inserted vertices go
first); the default weights are kept and the recipe is documented, not imposed.

### 6.2 Post-refinement

Delaunay refinement of the coarse mesh for element quality (corner angles ≥ 25° when `M̄` is closed
and every input angle sum `Θ_i ≥ 60°`, Gillespie et al. 2021a; otherwise best effort, reported).
Done by treating the coarse mesh as the input of a new intrinsic triangulation (its `Complex` and
lengths, **if it is simplicial**: the input of `ddg-integer-coords` / `ddg-intrinsic` is a `TriMesh`,
not a Δ-complex, see "As built" below for what is refused), refining it, and
composing the maps: a tracked point `Face(f, b)` of the coarse mesh is a point of that input and is
located in the refined triangulation by its `T⁰ → T¹` query. The result is again a `CoarseMesh`
(lengths, angular coordinates rebuilt, inserted vertices get new ids `≥` all previous ids).

*As built (WP5; `src/refine.rs`, `tests/refine.rs`, `examples/refine_study.rs`).*

* **API.** `CoarseMesh::from_intrinsic<T: LocatedTriangulation>(&T, &[InputPoint], options)`
  (§3.1; `LocatedTriangulation: ddg_intrinsic_apps::representation::InputPoints` with a batch
  query, implemented for `IntegerTriangulation`, which caches the common subdivision per input
  face, and `IntrinsicTriangulation`); `CoarseMesh::from_refined(mesh, geometry, points,
  mesh_options, &RefinementOptions)` = build `T⁰`, refine, `from_intrinsic`; `CoarseMesh::post_refine(&self,
  &RefinementOptions) -> (CoarseMesh, RefinementReport)` (`self` unchanged). `RefinementOptions`:
  `representation` (default `Representation::DEFAULT`), the shared driver's `RefineOptions`
  (25°), `rules` (`None`: the representation's default, `Chew` for integers, `Ruppert` for
  signposts), `certify`, `make_simplicial_flips`. Errors: `RefineError`
  (`Coarse(IemError)`, `NotSimplicial`, `Integer`, `Intrinsic`, `Locate { point, what }`).
  Dependencies added: `ddg-integer-coords` and `ddg-intrinsic-apps` (the `InputPoints` trait and
  `Representation`; moved from dev-dependencies), both workspace crates.
* **Ids (S2).** Pre-refinement: fine ids are the ids of `T¹` (the input vertices `0..|V|`, then the
  inserted ones in the order of `T¹`'s vertices; Chew's removals renumber, so "creation order"
  means the final order). `n_fine_vertices` stays `|V|`; the tracked points are the input vertices
  (as `Vertex`) and the extra points (as `Face`); inserted vertices are neither tracked nor fixed
  (`removal.rs` already only re-expresses a removed vertex that is `< n_input`). Post-refinement:
  old vertices keep their fine ids, inserted ones get `complex_of.len() + (v − m)` for the vertex
  `v` of `T¹` (`m` = old vertex count), all tracked points keep their ids, the fixed mask is
  carried over (new vertices free).
* **Located points** are put in a face (S16 style): a result on a vertex goes to the face of the
  vertex's outgoing halfedge with a unit barycentric vector, one on an edge to the face of the edge's
  halfedge; negatives down to `−1e-9` are clamped to 0 and counted (`n_clamped`), further out is
  `RefineError::Locate`; coordinates are renormalized to sum 1.
* **Frames.** The angular coordinates are rebuilt from scratch (S4: the reference halfedge is the
  new complex's outgoing halfedge, `φ = 0` there). The new frames are not the old ones and no
  rotation is kept; `P` and `P^vec` of the new mesh are computed from its own frames and the fine
  frames, so they are consistent without it (a coarse vector field expressed in the old frames
  would need re-expressing; not provided).
* **What is refused (post-refinement).** The refinement crates take a `TriMesh` (distinct vertices
  in a face, no two faces on the same three vertices, one edge per vertex pair), a coarse mesh is a
  Δ-complex. `as_simplicial_input()` returns `NotSimplicial { what }` naming the face or edge (fine
  ids) for: a face with a repeated vertex (self-edge/loop), two faces on the same three vertices,
  two different edges between the same two vertices (bigon), any other `TriMesh` error. Nothing is
  repaired silently. Option `make_simplicial_flips: Some(n)`: on a working copy,
  `CoarseMesh::make_simplicial(n)` greedily flips (tracked, metric and map unchanged) an edge at
  or next to a defect when that lowers the number of defects (self-faces + surplus parallel edges +
  surplus faces) and the flip is geometrically possible; the report has the flips and defects; if
  defects remain the refusal stands. Pre-refinement has no such restriction (its input is a
  `TriMesh`). Measured (random removal order down to a fraction of the pre-refined mesh, README):
  refused without the option on blub only at 2 %, on spot at 5 %, on boxpart at every fraction
  (up to 50 %); with the option all but boxpart at ≤ 10 %.
* **Guarantee.** `RefinementReport::guarantee_applies` = the input of the refinement is closed and
  every angle sum ≥ 60°; then the bound is reached (tested: 25.00–25.10° on every closed model).
  `min_angle` is reported either way, as are `refine.min_angle_reachable` (the smallest angle at
  vertices that are not narrow) and the shared driver's report.
* **Measured** (single thread, release, machine load ≈ 3–5 from other agents): blub (7106 vertices):
  pre-refinement +1771 vertices (the same count as the reference; reference 33 s, here 0.029 s
  refine + 0.015 s build/locate, signposts 0.027 s), 14208 extra points located in 0.07 s
  (integers); post-refinement of 2219 coarse vertices +1360 in 0.043 s.
* **Convenience.** "Refine, then coarsen" in one call is `from_refined` followed by `simplify`
  (WP3); it is to be wired when `simplify` lands. The reference recipe's area weight 0.1 is not
  imposed.
* **Signposts limitation found.** The signposts' `T⁰ → T¹` query fails ("ran into a vertex") for a
  point exactly on an input boundary edge of a coarse mesh whose boundary vertices are straight only
  to rounding (it continues through a vertex only if the angle sum `== π`): the planar-grid
  tests locate such points (removed boundary vertices) and so run post-refinement on integer
  coordinates only, signposts report `RefineError::Locate` naming the point. The integer
  representation has no such failure.

---------------------------------------------------------------------------------------------------

## 7. Mapping and prolongation (paper Sec. 7)

### 7.1 Pointwise map

`point_location(b) -> CoarsePoint` (`Vertex(v̄)` or `Face(f̄, b̄)`) for every tracked point, with
the coarse face given in the compacted output ids. Bijectivity is by construction: flips and
removals re-express points in a planar layout of a flat region, flattening is Eq. 10.

### 7.2 Scalar prolongation `P` (Sec. 7.2)

`P ∈ ℝ^{|V| × |V̄|}`: row `i` has `1` at `i`'s column if `i` is a coarse vertex, else `b_k` at the
three corners of its coarse face (summed if a corner repeats in a Δ-face). Rows sum to 1.

### 7.3 Vector prolongation `P^vec` (Sec. 7.3)

Literal paper (S18): coarse vectors are given in the **unnormalized** frames of the coarse vertices
(`u_i = ‖u⁰_i‖ e^{i Θ_i φ_i / 2π}` from the normalized `u⁰_i`; helpers `normalize`/`unnormalize`
are provided; the reference skips this step, R28). Coarse interpolation in each face with the
`β`, `w`, `z` of Sec. 7.3.1; fine sampling with `ψ_i = w_ij / w̄_ij` averaged over the fine
neighbours `j` whose images lie in the same coarse face as `x_i` (paper), then — if there are none
— over fine neighbours in an edge-adjacent coarse face (laid out across the shared edge; the
reference's addition R29 with its sign bug R31 fixed), then — if there are none — the average of
the `ψ` of the fine neighbours already known (breadth-first, paper's "average of known interpolated
values"; never an empty row, R29). `ψ` is normalized to unit modulus. A fine vertex that is a coarse
vertex gets `ψ` from its fine neighbours' images in the coarse faces around it (not `1`; R29). The
output is in the fine mesh's normalized frames, the frames of `ddg_intrinsic::IntrinsicTriangulation`
on the input (so the result plugs into `ddg-intrinsic-apps`' vector code). n-direction fields
(Knöppel encoding) use `z^n`, `ψ^n`.

*As built (WP4).* `CoarseMesh::compact() -> CoarseOutput` (`src/output.rs`; the complex's ids
are already dense, so the compacted ids are its vertex and face ids); `src/prolongation.rs`
holds `P`, `P^vec`, `FineFrames` and the helpers. Decisions where the text left room:

* **Fine frames.** The fine side of `ψ` is `e^{iφ_ij}` with `φ_ij` the *normalized* signpost of
  the fine edge in the frame of `i` (the paper's "unit vector along the fine edge in the
  normalized coordinate system at `i`"), so the output is in the normalized frames
  (`2π/Θ` inside, `π/Θ` on the boundary), and the matrix acts on *unnormalized* coarse vectors.
  `ddg-intrinsic-apps` stores a boundary vertex's vectors in the unnormalized frame; the
  conversion is `unnormalize`/`normalize` (boundary scale `π/Θ`).
* **Boundary unnormalization** generalizes the paper's `Θ/2π` to `Θ/π` at a boundary vertex:
  `w = e^{iφ/s}`, `s` the scale of S4.
* **Neighbours that are coarse vertices** count as lying in every face around them (at the
  corner); the reference ignores them (R29). A coarse vertex's row has one entry, the average
  over the faces around it and the neighbours in them of `ψ z`, normalized.
* **Fallback 3 differs from S18's text:** instead of averaging the `ψ` of the known neighbours
  (which live in the frames of different coarse faces, so their average has no meaning) the row
  is the average of the known neighbours' *rows*, each carried to the frame of `i` by the fine
  transport `e^{i n (φ_ij − φ_ji − π)}` (the paper's "average of known interpolated values"),
  breadth-first in rounds. A tier whose `Σ ψ_ij` cancels (`< 1e-8` per neighbour) counts as
  failed. A last resort (no row known in a component) takes the face's own corners with
  `ψ = 1` and is counted in `n_default` (never happened).
* Edge-adjacent faces are laid out across the shared side from the two faces' lengths (no
  obtuse-sign bug, R31); a face adjacent to itself (Δ-complex) is unfolded the same way.
* Fallback 3 uses the normalized angles as they are, which is exact only where the scale is one
  (boundary vertices of cone angle `≠ π`); it fires for a handful of rows only.
* The `β` of the paper are the angles of the sides in the face frame (`0`, `π − θ_j`,
  `π + θ_i`), the printed subscripts of the paper being garbled.

### 7.4 Drawing the coarse mesh on the fine mesh (Sec. 8.1)

* `face_of_points(points)`: coarse face per tracked point; a greedy colouring of the coarse faces
  (adjacent through a vertex → different colours; deterministic; no fixed palette size: the
  colouring returns colour indices, the caller maps them) — the reference's 24-colour overflow
  (R40) cannot occur.
* *As built (WP4).* `CoarseOutput::{face_of_point, face_of_points, colour_faces, point_colours}`:
  a point that is a coarse vertex belongs to its incident face of lowest id; the colouring visits
  the faces in id order and gives each the smallest colour absent around its corners.
* The Blender view (§11) tracks one point per fine face (its centroid) and colours each fine face
  by its coarse face; with a UV map it can instead bake texels (Sec. 8.1, Fig. 15).

---------------------------------------------------------------------------------------------------

## 8. Adaptive and anisotropic coarsening, constraints (paper Sec. 8.6)

* **Spatial adaptivity** (8.6.1): `PerVertex` channels (§4.1).
* **Anisotropy** (8.6.2): before simplification, every edge length is scaled by
  `(1 − τ) + (τ/2)((u_i · ê_ij)² + (u_j · ê_ji)²)`, `τ ∈ [0, 1]`, `u_i` unit tangent vectors at the
  vertices (given in the normalized frames of the input, or as 3D directions projected into the
  vertex's tangent plane and converted to an angle: S22), `ê_ij` the unit edge vector at `i`.
  Lengths that violate the triangle inequality are an error listing the faces (the paper notes this
  limits the strength; option `clamp_tau` lowers `τ` globally to the largest valid value and
  reports it). The fine points still refer to the unscaled input; the metric being coarsened is the
  scaled one.
* **Principal curvature directions** (for anisotropy, Fig. 27): Panozzo, Puppo, Rocca 2010
  (quadric fitting in a local frame over a geodesic-ish neighbourhood at scale `k ×` mean edge
  length; implemented from the paper; libigl's `principal_curvature` is not copied, S22). Output:
  `k_max`, `k_min`, unit directions per vertex.
* **Boundary preservation / constraints** (8.6.3, Fig. 28): fixed vertices (§2.2).

*As built (WP8; `src/anisotropy.rs`).*
* **Frames.** `u_i` is a normalized frame angle (`curvature::direction_angles`). The edge angle at
  `i` is rebuilt from the corner angles of the input lengths, counter-clockwise from
  `outgoing_halfedge(i)` (the unnormalized angle `θ`; the closing edge of a boundary fan is at
  `Θ`); `u` is unnormalized by `Θ/2π` (interior) or `Θ/π` (boundary), so the scalar product is
  taken in the real tangent plane. At a cone vertex the angle between `u` and the edge is the
  smaller of the two ways round, `min(δ, Θ − δ)`; on the boundary `|δ|`. `(u·ê)² = cos²` is
  invariant under `u → −u` at flat vertices. A direction pointing out of a boundary vertex's fan is
  snapped to the nearest boundary edge by `direction_angles`; this is its rule, not ours.
* **Direction of the effect.** The factor is `1` for edges parallel to `u` and `1 − τ` for edges
  orthogonal to it, so coarse triangles (isotropic in the scaled metric) are **elongated across
  `u`**: `u = d_max` → long along `d_min`.
* **Validity is exact, not bisected.** Each scaled length is affine in `τ`, so every triangle
  inequality is affine in `τ` and holds on `[0, τ*)`; `max_valid_tau` is the minimum of the roots
  (a test compares it with bisection to `1e-9`). `clamp_tau` uses `τ_max (1 − clamp_safety)`
  (default `1e-6`) and reports `tau_used`, `tau_max`, the violating-face count and the factor
  range in `AnisotropyReport`; without it, `AnisotropyError::TriangleInequality` names all faces.
  `τ = 0` returns the input lengths bitwise.
* **Tracking** needs no change: `CoarseMesh::from_lengths` takes the scaled lengths and the
  points as (input face, barycentric coordinates), which are metric independent. Tests: the
  `τ = 0` mesh is `bitwise_eq` to the plain one; `certify` passes after coarsening; every point is
  located; the total area equals that of the scaled mesh.
* **Measured** (`examples/anisotropy_study.rs`, 300 coarse vertices, stand-in order
  `coarsen_by_area`: lightest lumped area first, in rounds; mean angle in degrees between the
  long axis of the coarse faces, from the covariance of their tracked fine vertices, and the
  direction line; 45 = no relation, 90 = across). Open cylinder 96x60 (`Max`): 31 at τ = 0 →
  74 at τ = 0.6, and 16 to `d_min`. Spot `Max`: `τ_max = 0.275`, 46 → 55 (τ = 0.2) → 57 (clamped);
  `Min`: `τ_max = 0.360`, 44 → 61. Bob `Max`: `τ_max = 0.463`, 48 → 69 (τ = 0.4) → 71;
  `Min`: `τ_max = 0.498`, 42 → 74. Capsule (reference): `τ_max = 3e-4` for both (the input has
  nearly degenerate triangles), no effect possible. Mean aspect of the faces 1.7 → 2.1-2.4 (cylinder
  3-4.5); mean minimum corner angle of the coarse faces in the scaled metric stays 41 degrees, in
  3D (chord triangles) it falls from 41 to 34.

---------------------------------------------------------------------------------------------------

## 9. Hierarchies, multigrid, mean curvature flow (paper Sec. 8.7)

### 9.1 Hierarchy

Levels `|V| = n₀ > n₁ > … > n_L ≥ n_min` (default ratio ¼, `n_min = 500`, as the reference), from
one `simplify` run with snapshots (S24; the reference re-runs a fresh coarsening per level, R33,
available as an option for comparison). `P_k ∈ ℝ^{n_{k} × n_{k+1}}` from the tracked points of
level-`k` vertices in level `k + 1`.

### 9.2 Multigrid solver (Liu et al. 2021)

For an SPD system `A x = b` on the fine mesh (e.g. `L + εM`, or `M − τL` for flows), with
Dirichlet rows: Galerkin coarse operators `A_{k+1} = P_kᵀ A_k P_k` (`iem-galerkin-spd`: SPD if `A_k`
SPD and `P_k` has full column rank, which `P` rows of the coarse vertices guarantee), V-cycles with
`ν₁ = ν₂ = 2` forward/backward Gauss–Seidel sweeps (symmetric smoothing so the cycle can also
precondition CG), an exact sparse Cholesky solve at the coarsest level, stopping at a relative
residual (default `1e-10`) or a cycle cap, with a report (residual per cycle, convergence factor).
Dirichlet values are eliminated as in the reference (`min_quad_with_fixed_mg`), the dropped
columns cascading through the levels. Compared against a direct solve.

### 9.3 Conformalized mean curvature flow (Kazhdan et al. 2012, Fig. 30)

`(M_t − h L₀) X_{t+1} = M_t X_t` with the fixed cotan Laplacian `L₀` of the input and the mass
matrix `M_t` of the current positions (Kazhdan et al. Sec. 5), solved per step by the multigrid of
§9.2 on the **intrinsic** hierarchy of the input (built once). Reported: time per step against a
direct solve, and the difference of the results.

*As built (WP7, `src/hierarchy.rs`, `src/cmcf.rs`; study `examples/multigrid_study.rs`, flow
`examples/cmcf.rs`; tests `tests/hierarchy.rs`).*

**Hierarchy.** `build_hierarchy(&CoarseMesh, &HierarchyOptions) -> Result<Hierarchy, HierarchyError>`;
options `ratio` (0.25), `n_min` (500), `fresh_levels` (false), `simplify` (channels, metric, fixed
vertices; `target`, `snapshots` and `record_sequence` are overwritten). The mesh is cloned, not changed.
Level counts as the reference (`mg_precompute.cpp`): `n_{k+1} = round(ρ n_k)` while `n_k > n_min`, so the
coarsest level is the first with at most `n_min` vertices (spot 2930, 733, 183; bob 5344, 1336, 334;
blub 7106, 1777, 444; spot_subdiv 24755, 6189, 1547, 387). One run: `Target::Vertices(last)` with
`snapshots = [n_1, …]` (S24); `fresh_levels`: one `simplify` call per level on the previous result with
one snapshot at its target (each call recomputes masses and error vectors from the current curvature, which
is R33; the tracked fine vertices carry the locations). `Hierarchy`: `levels()` (`Level { vertices (fine
ids), faces, lengths, laplacian, mass }`), `prolongations()` (`P_k`, `n_k × n_{k+1}`), `report()`
(targets, counts, the `SimplifyReport`s, nnz of the `P_k`, seconds), `multigrid(A, known, options)`
(= `Multigrid::new(A, P_k, …)`). Level 0 is the mesh as given, **before** the initial Delaunay flips of the
run (its Laplacian is that of the input connectivity); `P_0` rows are indexed by the mesh's vertices.
Level Laplacians and lumped masses are computed from each snapshot's faces and side lengths (self-edges
contribute nothing, so they are safe on Δ-complexes); the test compares them with
`CoarseMesh::cotan_laplacian`/`lumped_mass` of the same run stopped at that count. `P_k` rows: 1 at the own
column for a kept vertex, else the barycentric coordinates of the coarse face (summed on a repeated
corner); `HierarchyError` for a ratio outside (0,1), a refused/failed `simplify`, a level whose lengths are
not a triangulation, an inconsistent snapshot.

**Multigrid on it** (`examples/multigrid_study.rs`; V-cycle `ν₁ = ν₂ = 2` symmetric Gauss–Seidel, tolerance
`1e-10` relative; `A = L + M/Σm` or the reference's Dirichlet example `L z = M 1`, `z = 0` at 0, 1, 2,
random start; 1 thread; load average 2–3 from other work). Cycles to `1e-10`, convergence factor (geometric
mean), and PCG iterations (the V-cycle as preconditioner), one-run snapshots:

| mesh (levels) | problem | V-cycles / factor | PCG its (factor) | V 1+1 | V 3+3 |
|---|---|---|---|---|---|
| spot (2930, 733, 183) | `L+M` | 19 / 0.29 | 11 (0.11) | 45 / 0.59 | 13 / 0.17 |
| spot | Dirichlet | 26 / 0.30 | 14 (0.10) | 65 / 0.62 | 20 / 0.20 |
| bob (5344, 1336, 334) | `L+M` | 16 / 0.23 | 10 (0.08) | 43 / 0.58 | 11 / 0.12 |
| bob | Dirichlet | 22 / 0.23 | 12 (0.07) | 58 / 0.59 | 19 / 0.19 |
| blub (7106, 1777, 444) | `L+M` | 49 / 0.62 | 17 (0.25) | 138 / 0.85 | 31 / 0.47 |
| blub | Dirichlet | 66 / 0.63 | 22 (0.23) | 186 / 0.85 | 42 / 0.48 |
| spot_subdiv (24755 … 387) | `L+M` | 46 / 0.61 | 16 (0.24) | 102 / 0.80 | 31 / 0.47 |
| spot_subdiv | Dirichlet | 60 / 0.58 | 21 (0.21) | 138 / 0.79 | 37 / 0.42 |

The solutions agree with the direct solve to `1e-11` relative or better on every row.
The reference's example (spot, Dirichlet): residual `157 → 4.2e-11` in 19 cycles (factor 0.22); here
`149 → 4.2e-11` in 23 cycles (0.29; the reference uses forward sweeps twice, this solver forward then
backward so that the cycle can precondition CG). The earlier subdivision hierarchy of spot gave 0.58 per
cycle with 2+2 sweeps: on the intrinsic hierarchy it is 0.29. Blub (low-quality triangles) is the hardest,
as expected. **One run (S24) against fresh levels (R33): no measurable difference** (cycle counts equal or
within one on all rows; e.g. spot Dirichlet 26 / 26, bob 22 / 23, blub 66 / 66, spot_subdiv 60 / 60), also
with `reference_compat(0, …)` on spot (23 / 23 cycles to `4.2e-11`). The one-run hierarchy is kept as the
default (no repeated curvature computation). **Defaults:** 2+2 sweeps (the paper's), and PCG with one
V-cycle as the preconditioner as the recommended driver (a third to a half of the iterations of stand-alone
V-cycles at the same cost per iteration; 3+3 sweeps cut iterations further, but the PCG solve time is
44-48 ms for 1+1, 2+2 and 3+3 on spot_subdiv, so 2+2 stays). Timing of one solve (spot_subdiv, 24.7k, PCG 2+2): setup 16 ms, solve
44 ms, against 38 ms for submatrix + sparse Cholesky + solve of `ddg_core` (supernodal, fill-reducing):
this multigrid is not faster than the direct solve for a *single* system at these sizes; it pays when the
system changes per solve and the factorization cannot be reused (cMCF below), and at sizes where
`O(n^{3/2})` factorization dominates.

**cMCF** (`cmcf::CmcfFlow`): `(M_t + h L₀) X_{t+1} = M_t X_t` with this crate's PSD `L₀` (§9.3 prints
`M_t − h L₀` for a negative semidefinite stiffness matrix; same system), `M_t` the **unnormalized** lumped
mass of the current positions, so `h` is an area (the unit sphere obeys `r² = 1 − 4t`; tested within 3 %
for 10 backward-Euler steps, `h = 0.002`). `step_direct` (symbolic Cholesky once, numeric every step) and
`step_multigrid(&Hierarchy, &CmcfMultigridOptions)` (PCG + V-cycle per coordinate, previous positions as
the initial guess, Galerkin products and coarsest factorization rebuilt per step with the hierarchy's
`P_k`; tolerance `1e-10` relative, cap 100); `normalize()` re-centres and rescales to the input's area
(the flow otherwise reaches a point in finite time, a unit-area sphere at `t = 1/16π`, after which `M_t` is
singular; the example applies it after every step). Reports per step: assembly, setup, solve seconds, PCG
iterations per coordinate, convergence, the relative residual. Results (unit-area input, `h = 1e-3`,
rescaled each step, 2+2, load average 2–3):

| mesh | steps | direct s/step | multigrid s/step (assemble + Galerkin + PCG) | PCG its | one-step diff | 100-step trajectory diff |
|---|---|---|---|---|---|---|
| spot 2930 | 100 | 0.0028 | 0.0092 (0.0008 + 0.0019 + 0.0065) | 8 | `3e-11` | `6e-10` |
| bob 5344 | 100 | 0.0097 | 0.0155 (0.0014 + 0.0035 + 0.0106) | 7 | `4e-11` | `3e-10` |
| blub 7106 | 100 | 0.0110 | 0.0295 (0.0020 + 0.0045 + 0.0230) | 11-12 | `2e-11` | `2e-10` |
| spot_subdiv 24755 | 50 | 0.033 | 0.12 (0.007 + 0.015 + 0.10) | 11-12 | `1e-11` | `4e-11` |
| spot ×16 (46850, tol `1e-10`) | 3 | 0.063 | 0.165 | 10 | `1e-12` | `1e-12` |
| spot ×16 (tol `1e-5`) | 3 | 0.061 | 0.09 | 4 | `1.7e-7` | `2e-7` |
| spot ×64 (187394, tol `1e-5`) | 3 | 0.41 | 0.51 (0.05 + 0.19 + 0.27) | 4-5 | `1e-7` | `1.5e-7` |

Reusing the input hierarchy's `P_k` with fresh Galerkin products is enough: the iteration count is
constant over 100 steps while the shape goes from the input to a sphere (spot 9 → 8, bob 7 → 7, blub
15 → 11), and the multigrid steps agree with the direct steps to `1e-11` (one step from the same state) and
`≤ 6e-10` after 100 steps (the test asserts `1e-8`). The Galerkin rebuild is 20-25 % of a multigrid step
(0.19 s of 0.51 s at 187k; the same order as one numeric factorization, 0.29 s). **The paper's 20× speedup
is not reproduced**: against this crate's sequential supernodal Cholesky (0.41 s per step at 187k vertices,
which is much faster than the Eigen `SimplicialLDLT` of the authors' stack) the multigrid step is 0.3-0.6×
at 3-25k vertices and 0.8× at 187k with `1e-5`; the trend (0.27 → 0.5 → 0.8 from 3k to 187k) says a
crossover somewhere near a million vertices. What would help, not done: the split
`P^T(M_t + hL₀)P = h·P^T L₀ P (once) + P^T M_t P` (about a third of the setup), a threaded or colour-ordered
Gauss–Seidel, solving the three coordinates together. Meshes: `runs_tmp/cmcf_*/step_XXXX.obj` (every 10th
step) from `examples/cmcf.rs`.

---------------------------------------------------------------------------------------------------

## 10. Applications and evaluation (paper Sec. 8.2–8.5)

### 10.1 Applications (Sec. 8.5)

On the coarse mesh, prolongated with `P`:
* **Poisson** (Fig. 21): `L̄ x̄ = M̄ f̄` (`f̄ = Pᵀ`-restricted or user data), `x = P x̄`.
* **Geodesic distance** (Fig. 22, 23): single source (a fine vertex mapped to `M̄`; if it is not a
  coarse vertex, the coarse distance from a point in a face) by the heat method (`ddg-intrinsic-apps`
  heat on the coarse lengths) or exactly (window propagation, `ddg-intrinsic-apps::exact`, which
  WP6 generalizes from a `TriMesh` to a `Complex` with lengths: Δ-complexes, self-edges, face
  points as sources; S20).
* **All-pairs** (Fig. 24): `D̃ = P D̄ Pᵀ` (dense, `|V̄|²` exact runs).
* **Geodesic Voronoi** (Fig. 25): nearest-source labels from the coarse distances.
* **Exponential map** (Sec. 8.5.4): straightest geodesics traced on the coarse mesh from a mapped
  point (the tracing of `ddg-intrinsic`), mapped back to the fine mesh through the tracked map is
  *not* available (the map is fine → coarse): the result is reported as coarse points and compared
  on coarse points (S21).

*As built (WP6).* `src/applications.rs`: `CoarseApplications::new(&CoarseMesh)` (output, Δ-complex,
`P`, `L̄`, `M̄`); everything takes a `CoarseMesh`, independent of `simplify`.
* `CoarseOutput::{cotan_laplacian, lumped_mass, face_areas}`: per halfedge, self-edges skipped,
  as `ddg_intrinsic::generic_ops`; equal to `CoarseMesh::cotan_laplacian` (tested).
* **Poisson**: `PoissonSolver::new(L, M, fixed)` (Cholesky, any `L`, `M`: the fine problem uses it
  too). No fixed vertices: Neumann/mean-zero (the load is projected to sum zero, `rhs_shift`
  reported; `x` has `M`-mean zero; internally vertex 0 is pinned). `poisson(f̄, fixed, values)` is
  `L̄ x̄ = M̄ f̄`, `x = P x̄`. Fine data: `sample` (`f̄_i` = `f` at the fine vertex of `i`) or
  `restrict_load` (`b̄ = Pᵀ M f`, sums conserved); both measured, similar errors; the example uses
  `sample`.
* **Geodesics**: `geodesic_exact` = `source_of` (vertex or face point, S20) →
  `exact_distance_complex` → `P`. Heat: **reimplemented** as `CoarseHeat` for a Δ-complex with
  lengths (the heat solver of `ddg-intrinsic-apps` needs a `Triangulation`): `(M̄ + tL̄)u = δ`,
  `X = −∇u/|∇u|` from the planar layout of each face, `L̄φ = −∇·X`, `t = (factor)·h̄²` (factor 1),
  `φ` shifted to zero at the sources (mean of `φ` at their barycentric vertices). Factorizations
  once per `CoarseHeat`.
* **All pairs**: `all_pairs(pool, opts)` (a `ddg_core::sparse::ThreadPool`, one exact run per row,
  bitwise equal for 1 and 4 threads, tested); `prolongate_rows` (rows of `P D̄ Pᵀ`),
  `prolongate_dense` (`8|V|²` bytes). **Voronoi**: `voronoi(sources, pool, opts)`, one exact run per
  source, prolongated, argmin, ties to the lower index; coarse labels returned too. (The paper
  computes the diagram on the coarse mesh; the fine labels here come from the prolongated
  distances, as the brief of WP6 asks.)
* **Exponential map**: `exp_map_face(face, bary, vector in the face's planar frame)` and
  `exp_map_vertex(v, angle, length)` (the *unnormalized* angle `φ/s` of the angular coordinates, SPEC
  §2.2). Unfolding across faces from their lengths; a ray that hits a vertex (relative `1e-10`)
  continues at `θ_back + Θ/2` (Polthier–Schmies straightest geodesics); stops at the boundary
  (`ExpStop::Boundary`). Returns a coarse point, traced length, faces crossed, vertices hit and the
  final unit direction in the last face's frame (reversible, tested on a Δ-complex and on a flat
  grid against the straight line). It is written on the output's `twin`/`face_lengths`, not on
  `ddg-intrinsic`'s tracer (which needs a `Triangulation` with positions or signposts).
* `greedy_coarsening(&mut CoarseMesh, target)`: stand-in order (smallest memoryless cost with
  `α_ij = (K̃_j − K_j)/K_i`); not the ICE order.
* Measured numbers (greedy order, not ICE; the example `applications`) are in the WP6 report; the
  paper's ≈1 % at three orders of magnitude (20k-vertex mesh) is not reproduced by the stand-in
  order. Rerun with the ICE order (`simplify`): §10.5 (WP9a); the example's `--method` picks the
  order (`ice` default).
* The exponential-map comparison of the example is weak: the coarse and fine start directions are
  matched by the normalized angular coordinate (valid while the reference halfedge's direction is
  kept), and the endpoint distance is an upper bound (distance on the coarse mesh to the corners of
  the tracked endpoint's face plus the planar distance).

### 10.2 Benchmark (Sec. 8.2; Fig. 16, 17)

Thingi10k (`data/thingi10k/`, the 7745 manifold meshes of IC-SPEC Q1), preprocessed by intrinsic
Delaunay refinement at 25°, coarsened as far as possible (target 0, stop at `∞`), the ratio
`|V̄|/|V|` normalized by the number of components. Report: the share reaching ≤ 10 % and ≤ 1 %
(paper: 98 % and 84 %), failures by reason (paper: 0.1 % from floating point), time per removal,
peak memory. Fig. 16: time against `|V|` on a subdivided cow/spot (`10³…10⁶`) and against the
percent reduction on a 750k-vertex subdivision. Runs in parallel per mesh (64 cores; timings alone).

**As built (WP9b).** `refine::refine_then_simplify` (pre-refinement then `simplify`, one combined
report); `examples/thingi_simplify.rs` (one mesh, one JSON line: modes `refine`, `--no-refine`,
`--signposts`; edge lengths intrinsically mollified with `ε = 1e-5 h̄` as in the integer-coordinates
evaluation, because Thingi10k has zero-area faces which `Geometry` refuses; target `AllRemovable`);
`scripts/iem-thingi.sh run|summary` (phases by size, ≤ 24 jobs, `ulimit -v`, `timeout`) and
`scripts/iem_thingi_summary.py`. The ratio is `|V̄|/|V|` (the paper's normalization of both counts by
the number of components cancels); `ratio_mean_component` (mean of per-component final/input) is
reported too. A failure of a mesh is counted when the job gives no result; the vertices left
unremovable are "skips" by reason, commit failures (reverted) separately.

**Results of record** (7745 meshes, default ICE metric, 24 parallel jobs on a shared machine; load
about 24 from this run, plus other agents'; timings indicative only; JSON per mesh and
`summary.md` in `runs_tmp/iem/thingi/`):

| | pre-refinement (integer, Chew, 25°) | none |
|---|---|---|
| completed, no crash/panic/timeout | 7745 | 7745 |
| certificate failures | 0 | 0 |
| `|V̄|/|V|` ≤ 10 % (paper 98 %) | 97.29 % | 97.24 % |
| `|V̄|/|V|` ≤ 1 % (paper 84 %) | 76.48 % | 76.50 % |
| same, meshes with `|V|` ≥ 300 (6366) | 99.62 %, 89.38 % | 99.59 %, 89.32 % |
| same, meshes with `|V|` ≥ 10⁴ (1745) | 99.66 %, 98.45 % | 99.66 %, 98.45 % |
| mean per-component ratio ≤ 10 % / ≤ 1 % | 94.47 % / 72.47 % | 94.42 % / 72.54 % |
| median / p90 / p99 ratio | 0.0027 / 0.027 / 0.214 | 0.0026 / 0.026 / 0.214 |
| meshes with a skip by flattening failure (`flatten_failed`) | 836 (10.8 %) | 814 (10.5 %) |
| meshes with reverted commit failures | 679 (8.8 %) | 619 (8.0 %) |
| removals per second, simplify phase (median / p10) | 12.0k / 9.7k | 12.4k / 10.3k |
| total time per mesh, median / p99 / max | 0.59 s / 51 s / 671 s | 0.15 s / 33 s / 340 s |
| peak RSS, median / max | 8 MB / 2.1 GB | 5 MB / 1.4 GB |

Findings. (1) All meshes finish; every coarse mesh passes `certify`. (2) The shares below the paper's
are the floor of 3 vertices per closed component (a tetrahedron cannot go below 0.75, a mesh of
100 vertices not below 3 %), not failures: the sets that can reach the bound do (99.6 % ≤ 10 %,
89 %-98 % ≤ 1 %). The remaining 10 %-outliers are tiny meshes, many-component meshes (3 vertices
each), and a few high-genus ones (13 closed meshes of total genus ≥ 10 above 10 %), as the paper
says for Fig. 17. (3) Numerical trouble is not at the paper's 0.1 %: 7-11 % of the meshes
leave at least one vertex unremovable because its flattening failed (`flatten_failed`, a
non-convergent Newton iteration) or is `not_flattenable` (a triangle inequality on the star
violated after the flips); most of them are meshes that still reach ≤ 10 %; only 14-19 of the
meshes above 1 % are dominated by these reasons. We did not separate floating point from true
topological obstruction for `not_flattenable` (the paper does not say how it counted). (4)
Pre-refinement changes the ratio on 2186 of 7745 meshes (smaller on 1051, larger on 1135), the
share ≤ 10 % by 0.05 %: it gives candidate locations but does not help the ratio on average; it
costs 2.4× the total CPU time (8.7 h against 3.7 h; median per mesh 4×; 9 % of it the refinement
itself, the rest simplifying the inserted vertices) and inserts a median 2.2× `|V|`. (5) Speed matches the paper's "about 10 000 vertices per
second". Slowest: 719791 (871k vertices, 671 s with pre-refinement).

### 10.3 Distortion (Sec. 8.3, Fig. 18)

For every fine edge `ij`: `ℓ̂_ij` = the length of the image path, approximated (paper) by the
distance between `φ(i)` and `φ(j)` on `M̄`: in the layout of the coarse face if both are in one
face, else the exact geodesic distance on `M̄` (§10.1). For each fine face, the triangle with
lengths `ℓ̂` against the input triangle: with `γ_max ≥ γ_min` the eigenvalues of the metric tensor
`JᵀJ` of the affine map (`γ = σ²`), the anisotropic distortion `γ_max/γ_min` and the area distortion
`log(γ_max γ_min)`/2 `= log(σ₁σ₂)` as Khodakovsky et al. 2003 Sec. 2 define them (the paper says
"as defined by Khodakovsky et al."); `σ₁/σ₂` is reported too. Mean and max (`iem-distortion`).

### 10.4 Extrinsic baselines (Sec. 8.4; Fig. 18, 19, 21, 22)

* **QEM** (Garland, Heckbert 1997): edge collapse by quadrics with optimal placement (the 3×3
  system, else the optimum on the segment, else the best of the endpoints and the midpoint; Sec. 4),
  boundary planes for open meshes with a large penalty (Sec. 6), rejecting collapses that break
  the link condition or flip a face normal; to a target vertex count. Defaults as in the paper:
  unweighted plane quadrics (area weighting is an option), unbounded optimal placement, candidate
  pairs are the edges only (`t = 0`, Sec. 3.2). Area error (Fig. 19).
* **Map for the extrinsic result**, for prolongation-based comparisons: closest-point projection of
  the fine vertices onto the coarse extrinsic mesh (simple, documented as ours), and the successive
  self-parameterization of Liu et al. 2021 (the paper's Fig. 18 and 21 baseline) built during the
  QEM collapses.
* Comparisons of record (Fig. 18, 19, 21, 22): distortion, area error, Poisson error, geodesic
  distance error, at equal vertex counts, ICE against QEM.

*As built (WP9a).* `src/evaluation.rs`, `src/baseline/successive.rs`, examples `evaluate`,
`applications --method`, `scaling`.

* **Distortion** (`evaluation::distortion(positions, faces, &CoarseApplications, pool, opts)`,
  `[formula: iem-distortion]`): the map is the tracked points `output.points[..n_fine]`.
  `ℓ̂_ij` per unique fine edge by `PointDistance`, a **local unfolding search** instead of one
  exact run per fine vertex (too slow, and `exact_distance_complex` gives distances to
  vertices, not between face points): from the image `s` faces are unfolded across their sides
  carrying the visibility cone (an exact-window propagation without merging), the target image
  `t` is tested inside each unfolded copy of its face (straight paths touching no vertex);
  vertices seen inside a cone give `d_s(c)`, the same search from `t` gives `d_t(c)`, and
  `min_c d_s(c) + d_t(c)` is the best path with one bend at a vertex; branches farther than the
  best length are pruned; a pair in one face is searched too (on a sharp cone the way round the
  apex can be shorter than the segment). It is **exact** for shortest paths crossing at most
  `max_depth = 16` faces with at most one bend and an **upper bound** otherwise; if nothing is
  found the vertex route `|s − a| + d(a, b) + |b − t|` through exact vertex distances is used
  (`n_fallback`; zero in every run of record). Points on a side or corner are nudged `1e-9` into
  the face. Tested against the closed forms on a plane and on cones of angle sum below and above
  `2π` (400 random pairs, error `< 1e-6`). Per fine face: `J` from the planar layout of the input
  lengths to that of `ℓ̂` (both counter-clockwise), `σ₁ ≥ σ₂` of `J`, `γ = σ²`; reported:
  `γ_max/γ_min`, `log(σ₁σ₂)` (signed and absolute) and `σ₁/σ₂`; mean weighted by the input face
  area, maximum, median and 99th percentile in the example; faces with `ℓ̂` violating the strict
  triangle inequality are counted (`n_invalid_faces`, `invalid_area_fraction`) and excluded from
  the statistics. The **mean of `γ_max/γ_min` is dominated by a few near-degenerate images**
  (QEM has faces of `10^15`); read it with the median and the 99th percentile.
* **Eq. 10 sign (S25).** `CoarseMeshOptions::projective_sign: ProjectiveSign { Positive, Negative }`
  (default `Negative` since the user's decision of 2026-11-01, S25), `formulas::projective_interpolation_signed`; the coarse mesh and the removal
  order do not depend on it (tested), only the map.
* **Successive self-parameterization (S23)** in `baseline::successive`, switched on by
  `QemOptions::successive_map` (`QemResult::successive`: face and barycentric coordinates of each
  fine vertex, `QemResult::successive_stats`). Per collapse `(u, v) → k` of the QEM run: the faces
  around `u` and `v` before and around `k` after are flattened jointly by the **least-squares
  conformal energy** (the variant Liu et al. found best; the boundary vertices are shared
  unknowns, two far apart ones pinned), every tracked point is lifted to the shared plane and
  dropped into the after-face that contains it, so the maps compose over the run. Gaps against
  Liu et al. 2021 (the validity test is in their earlier paper, App. B, which we did not have):
  *closed meshes only* (the boundary cases of Sec. 4.3 are not built; an open mesh is an error);
  a layout with a non-positive triangle falls back to the **Tutte embedding** (circle, chord-length
  parameter, uniform weights), and if that fails too the collapse is rejected; no ARAP. In
  the runs of record (114.8k collapses) the Tutte fallback fired 5 times, no collapse was rejected
  and no point was clamped. Cost: the relocation
  of all points in the patch per collapse (the mean fine-to-image distance equals the
  closest-point projection's to `5 %` on the sphere test).
* **Extrinsic result as applications:** `evaluation::extrinsic_applications(&QemResult, n_fine,
  ExtrinsicMap::{ClosestPoint, Successive})` builds a `CoarseMesh` from the QEM mesh (extrinsic
  lengths), replaces the tracked points by the images and `QemResult::origin` (the fine vertex
  each coarse vertex is, for `sample`), and wraps it with `CoarseApplications::from_output`
  (new, a split of `CoarseApplications::new`). `qem_to_target` runs QEM to a target with a
  ladder of relaxations if the defaults stop short.
* **QEM changes made on the way (WP9a):** (1) the link test rejected every collapse of a vertex of
  degree three (its face `(u, a, b)` between the two opposite vertices): it now rejects only if the
  other end already has the face `(x, a, b)` (the tetrahedron); before, runs stalled far above the
  target (spot_subdiv: at 318 of 200). (2) When the queue runs dry before the target, the edges are
  queued again while the last sweep collapsed something (a rejected edge only came back when a
  neighbour changed). (3) `QemResult::origin`.

### 10.5 Results of record (WP9a)

Raw outputs: `runs_tmp/iem/wp9/` (gitignored; `dist_*.log` distortion and area, `perface/` the
per-face distortions, `app_*.log` applications, `area_*.log`, `scaling_*.tsv`; the driver scripts
are next to them). Commands: `examples/evaluate`, `examples/applications --method …`,
`examples/scaling`. Meshes: `data/models/{spot,bob,blub}` (triangulated), `spot_subdiv` (24 755
vertices) and `bs_rest`, `CrumpledDevelopable` as named. Not here: the Thingi10k benchmark (§10.2,
WP9b, waits for the pre-refinement of WP5). **Substitutions** (the bunny and the paper's meshes
are not in the tree): Fig. 18's 28k bunny at 200 vertices is replaced by `spot_subdiv` (24.7k) and
`blub` (7.1k) at 200; Fig. 19's developable of Rabinovich et al. 2018 by `CrumpledDevelopable`
(289 vertices, Crane's repository) and by `bs_rest` (an open sheet of 19 985 vertices from the
reference's meshes; it may be the reference's cloth, not verified); Fig. 21's cloth mesh by
`bs_rest` at 500 vertices and by closed meshes. ICE is `simplify` with the defaults; QEM is the
unweighted Garland–Heckbert of §10.4 (reached every target here with its defaults after the fixes
listed there); "QEM-s" is the successive self-parameterization (closed meshes), "QEM-c" the
closest-point projection. Equal vertex counts (ICE's target is `ceil`, so ICE may have one more
vertex). The machine was loaded by another job (load average 14–33) for all of this; no number
below is a timing except in the scaling tables.

**Distortion and area error (Fig. 18, 19; `evaluate`).** "aniso" is `γ_max/γ_min` (area-weighted
mean, median, 99th percentile and maximum over the fine faces), "abs log-area" is the mean of
`|log(σ₁σ₂)|`, "area err" is `100 |A_coarse − A_fine|/A_fine` with the coarse area from the coarse
lengths (ICE: intrinsic; QEM: extrinsic). "invalid" are faces whose `ℓ̂` violate the triangle
inequality (not in the statistics). The `ice-negative` rows are S25's other sign (same coarse mesh).
Edge distances: no pair needed the exact-vertex fallback or was unresolved in any run.

| mesh | n | method | area err % | aniso mean | median | p99 | max | abs log-area mean | max | invalid faces (area) |
|---|---|---|---|---|---|---|---|---|---|---|
| blub_triangulated | 142 | ice-positive | 2.26 | 1.138 | 1.117 | 1.942 | 7.489 | 0.055 | 0.917 | 0 (0.0e+00) |
| blub_triangulated | 142 | ice-negative | 2.26 | 1.052 | 1.036 | 1.465 | 4.485 | 0.059 | 0.781 | 0 (0.0e+00) |
| blub_triangulated | 142 | qem-closest | 0.11 | 6.6e+10 | 1.038 | 6.375 | 2.5e+15 | 0.045 | 17.752 | 6 (2.4e-04) |
| blub_triangulated | 142 | qem-successive | 0.11 | 1.382 | 1.107 | 6.917 | 945 | 0.056 | 2.554 | 0 (0.0e+00) |
| blub_triangulated | 200 | ice-positive | 1.56 | 1.091 | 1.078 | 1.582 | 11.016 | 0.039 | 1.239 | 0 (0.0e+00) |
| blub_triangulated | 200 | qem-closest | 0.24 | 6.6e+10 | 1.027 | 3.436 | 2.5e+15 | 0.031 | 17.752 | 4 (1.7e-04) |
| blub_triangulated | 200 | qem-successive | 0.24 | 1.136 | 1.070 | 2.970 | 178 | 0.033 | 2.489 | 0 (0.0e+00) |
| blub_triangulated | 711 | ice-positive | 0.42 | 1.018 | 1.017 | 1.107 | 1.213 | 0.009 | 0.144 | 0 (0.0e+00) |
| blub_triangulated | 711 | ice-negative | 0.42 | 1.010 | 1.008 | 1.087 | 1.187 | 0.010 | 0.115 | 0 (0.0e+00) |
| blub_triangulated | 711 | qem-closest | 0.03 | 1.019 | 1.009 | 1.428 | 3.569 | 0.008 | 0.636 | 0 (0.0e+00) |
| blub_triangulated | 711 | qem-successive | 0.03 | 1.022 | 1.015 | 1.307 | 2.968 | 0.008 | 0.469 | 0 (0.0e+00) |
| bob_tri | 107 | ice-positive | 0.83 | 1.183 | 1.157 | 2.051 | 10.491 | 0.082 | 1.420 | 0 (0.0e+00) |
| bob_tri | 107 | ice-negative | 0.83 | 1.060 | 1.045 | 1.521 | 6.184 | 0.089 | 1.624 | 0 (0.0e+00) |
| bob_tri | 107 | qem-closest | 0.62 | 1.197 | 1.080 | 3.530 | 27.289 | 0.063 | 1.639 | 0 (0.0e+00) |
| bob_tri | 107 | qem-successive | 0.62 | 1.208 | 1.152 | 3.216 | 28.724 | 0.056 | 1.306 | 0 (0.0e+00) |
| bob_tri | 200 | ice-positive | 0.47 | 1.094 | 1.084 | 1.429 | 3.242 | 0.047 | 0.632 | 0 (0.0e+00) |
| bob_tri | 200 | qem-closest | 0.44 | 1.086 | 1.046 | 1.790 | 4.971 | 0.032 | 0.702 | 0 (0.0e+00) |
| bob_tri | 200 | qem-successive | 0.44 | 1.094 | 1.075 | 1.654 | 4.309 | 0.028 | 0.478 | 0 (0.0e+00) |
| bob_tri | 534 | ice-positive | 0.18 | 1.031 | 1.028 | 1.121 | 1.778 | 0.017 | 0.268 | 0 (0.0e+00) |
| bob_tri | 534 | ice-negative | 0.18 | 1.016 | 1.011 | 1.078 | 1.354 | 0.017 | 0.250 | 0 (0.0e+00) |
| bob_tri | 534 | qem-closest | 0.13 | 1.026 | 1.018 | 1.181 | 1.669 | 0.010 | 0.248 | 0 (0.0e+00) |
| bob_tri | 534 | qem-successive | 0.13 | 1.030 | 1.025 | 1.190 | 1.784 | 0.010 | 0.202 | 0 (0.0e+00) |
| spot_subdiv | 200 | ice-positive | 2.01 | 1.195 | 1.134 | 2.135 | 9.786 | 0.075 | 1.248 | 0 (0.0e+00) |
| spot_subdiv | 200 | qem-closest | 1.10 | 1.663 | 1.032 | 5.561 | 23125 | 0.055 | 4.026 | 6 (1.0e-04) |
| spot_subdiv | 200 | qem-successive | 1.10 | 1.307 | 1.111 | 4.220 | 426 | 0.058 | 2.528 | 0 (0.0e+00) |
| spot_subdiv | 495 | ice-positive | 0.93 | 1.071 | 1.050 | 1.342 | 2.074 | 0.033 | 0.419 | 0 (0.0e+00) |
| spot_subdiv | 495 | ice-negative | 0.93 | 1.025 | 1.015 | 1.167 | 1.669 | 0.035 | 0.336 | 0 (0.0e+00) |
| spot_subdiv | 495 | qem-closest | 0.27 | 1.063 | 1.015 | 1.830 | 349 | 0.022 | 2.929 | 0 (0.0e+00) |
| spot_subdiv | 495 | qem-successive | 0.27 | 1.076 | 1.038 | 1.670 | 7.816 | 0.021 | 0.957 | 0 (0.0e+00) |
| spot_subdiv | 2476 | ice-positive | 0.22 | 1.011 | 1.008 | 1.047 | 1.113 | 0.006 | 0.060 | 0 (0.0e+00) |
| spot_subdiv | 2476 | ice-negative | 0.22 | 1.006 | 1.003 | 1.033 | 1.078 | 0.006 | 0.052 | 0 (0.0e+00) |
| spot_subdiv | 2476 | qem-closest | 0.03 | 1.009 | 1.004 | 1.085 | 3.535 | 0.004 | 0.629 | 0 (0.0e+00) |
| spot_subdiv | 2476 | qem-successive | 0.03 | 1.011 | 1.006 | 1.081 | 2.263 | 0.004 | 0.392 | 0 (0.0e+00) |
| spot_triangulated | 59 | ice-positive | 7.58 | 1.734 | 1.413 | 10.209 | 239 | 0.187 | 2.806 | 0 (0.0e+00) |
| spot_triangulated | 59 | ice-negative | 7.58 | 1.250 | 1.139 | 3.677 | 23.598 | 0.223 | 2.307 | 0 (0.0e+00) |
| spot_triangulated | 59 | qem-closest | 1.26 | 9.5e+11 | 1.163 | 73.476 | 1.1e+16 | 0.168 | 18.485 | 50 (2.8e-03) |
| spot_triangulated | 59 | qem-successive | 1.26 | 17.734 | 1.652 | 667 | 2.9e+05 | 0.205 | 5.654 | 0 (0.0e+00) |
| spot_triangulated | 200 | ice-positive | 1.85 | 1.155 | 1.134 | 2.027 | 6.860 | 0.060 | 0.990 | 0 (0.0e+00) |
| spot_triangulated | 200 | qem-closest | 0.53 | 1.113 | 1.054 | 3.979 | 18.696 | 0.038 | 1.440 | 0 (0.0e+00) |
| spot_triangulated | 200 | qem-successive | 0.53 | 1.170 | 1.116 | 3.275 | 19.023 | 0.048 | 1.497 | 0 (0.0e+00) |
| spot_triangulated | 293 | ice-positive | 1.26 | 1.101 | 1.090 | 1.544 | 2.245 | 0.042 | 0.431 | 2 (3.2e-04) |
| spot_triangulated | 293 | ice-negative | 1.26 | 1.056 | 1.041 | 1.378 | 1.969 | 0.047 | 0.594 | 2 (3.2e-04) |
| spot_triangulated | 293 | qem-closest | 0.36 | 1.073 | 1.037 | 2.513 | 6.209 | 0.027 | 0.909 | 0 (0.0e+00) |
| spot_triangulated | 293 | qem-successive | 0.36 | 1.102 | 1.073 | 2.287 | 5.612 | 0.033 | 0.745 | 0 (0.0e+00) |

Reading (all from the table; mean aniso of QEM-c is dominated by a few degenerate images, use the
median/p99/max):
* **Anisotropic distortion:** against QEM-s, ICE has the smaller or equal mean at the smallest
  target of every mesh (spot_subdiv 200: 1.195 against 1.307; blub 200: 1.091 against 1.136; blub
  142: 1.138 against 1.382; spot 59: 1.73 against 17.7; bob 200 equal), and a smaller or equal p99
  and maximum in every 2 %/200-vertex row (blub 142 max 7.5 against 945; spot_subdiv 200: 9.8
  against 426; QEM-c has maxima of 10⁴–10¹⁵ and invalid faces). But ICE's *median* is slightly
  larger at 200 vertices (spot 1.134 against 1.054 QEM-c / 1.116 QEM-s; bob 1.084 against 1.046 /
  1.075; blub 1.078 against 1.027 / 1.070), and at 10 % the means are comparable (bob 534:
  1.031 / 1.026 / 1.030). The paper's "modest reduction" holds for the tails and for the mean at
  strong reductions, not for the median.
* **Area distortion (`|log σ₁σ₂|`):** *not reproduced*: ICE's mean is larger than QEM's in every row
  except blub 142 and spot 59 against QEM-s (where QEM-s has degenerate faces; ICE's maximum is
  smaller there: spot 59 `2.8` against `5.7`; blub 200 `1.24` against `2.5`). Example: spot_subdiv
  200 `0.075` ICE against `0.055` / `0.058`; blub 200 `0.039` against `0.031` / `0.033`. ICE's
  `log σ₁σ₂` is biased low (signed mean −0.006 … −0.026 at 200 vertices, QEM ≈ 0): the intrinsic
  coarse mesh is smaller than the fine surface, by the area error (0.5–2.0 % at 200 vertices;
  QEM's 0.2–1.1 %). On closed smooth surfaces QEM's area error is the smaller (spot 293: 1.26 %
  ICE against 0.36 % QEM; blub 711: 0.42 against 0.03).
* **Developable (Fig. 19):** the result of the paper is reproduced. `CrumpledDevelopable` (289):
  area error ICE / QEM at 100 vertices `0.0001 %` / `5.05 %`, at 50 `0.023 %` / `9.46 %`, at 25
  `0.46 %` / `15.2 %`; aniso at 100: `1.0000` (an isometric coarse mesh) against 3.0 (max 674).
  `bs_rest`: 2000 `0.15 %` / `0.44 %`, 500 `0.46 %` / `2.10 %`, 200 `0.075 %` / `1.38 %`.

```
CrumpledDevelopable  vertices  method       area_err_%  aniso_mean  aniso_max  |logarea|_mean
                     100       ice          0.0001      1.0000      1.0000     0.0000
                     100       qem-closest  5.0473      3.0074      673.6      0.1350
                     50        ice          0.0227      1.0176      1.166      0.0050
                     50        qem-closest  9.4590      4280        2.0e6      0.3208
                     25        ice          0.4583      1.2244      3.348      0.0510
                     25        qem-closest  15.181      435.3       2.0e5      0.4226
```
(`area_crumpled.log`, `area_bs_rest.log`.)

**Eq. 10 sign (S25).** Same coarse mesh, same area error and lengths; only the map differs
(`ice-positive` against `ice-negative` rows above, at 10 % and 2 % of spot, bob, blub and spot_subdiv):
the `e^{−u}` map (`Negative`) has the smaller anisotropic mean *and* maximum in all 8 of 8
configurations (spot 293: 1.101 → 1.056, max 2.25 → 1.97; spot 59: 1.734 → 1.250, max 239 → 23.6;
bob 107: 1.183 → 1.060, max 10.5 → 6.2; blub 142: 1.138 → 1.052, max 7.5 → 4.5; spot_subdiv 495:
1.071 → 1.025), while the mean `|log σ₁σ₂|` is 3–19 % *larger* in 8 of 8 (spot 59: 0.187 → 0.223).
Decision (the user's, 2026-11-01): `Negative` is the default, `Positive` an option (S25). The
applications of §10.1 were not rerun with `Negative`.

**Applications (Fig. 21, 22, 23, 24, 25; `applications --method`).** Errors are against the fine
solution: exact geodesic and heat method = mean relative L1 error of the prolongated distance over
the fine vertices (5 sources); Poisson = relative L2 error (`f` smooth in the coordinates, `sample`
restriction); all pairs = relative L1 on 8 sampled rows of `P D̄ Pᵀ`; Voronoi = share of fine
vertices with a different nearest source (20 sources). Speedups are ratios of wall times, single
thread, on a loaded machine (indicative). Rows with the same size are the same row of the example
(ICE's target rounds up). `bs_rest` is open: QEM-s is not available, the QEM row is QEM-c.

### blub
# fine: exact distance 528.3 ms/source, Poisson 13.9 ms, Voronoi (20 sources) 10.81 s, all-pairs (extrapolated) 4319.4 s, exp map 0.00 ms (100 traces)
| n̄ | method | exact geod. L1 % | heat geod. L1 % | Poisson L2 % | all-pairs L1 % | Voronoi wrong % | speedup exact / Poisson / all-pairs / Voronoi |
|---|---|---|---|---|---|---|---|
| 143 | ice | 2.508 | 6.302 | 10.533 | 2.333 | 4.10 | 719 / 110.1 / 5596 / 613 |
| 143 | qem-closest | 1.579 | 8.679 | 9.427 | 1.843 | 4.78 | 754 / 99.7 / 5917 / 741 |
| 143 | qem-successive | 1.643 | 8.598 | 9.359 | 1.777 | 4.12 | 727 / 95.4 / 6787 / 755 |
| 200 | ice | 1.884 | 5.259 | 6.935 | 1.861 | 3.39 | 448 / 78.9 / 4595 / 442 |
| 200 | qem-closest | 1.314 | 9.344 | 7.532 | 1.414 | 2.69 | 466 / 72.7 / 5130 / 472 |
| 200 | qem-successive | 1.330 | 9.277 | 7.523 | 1.417 | 2.84 | 445 / 69.1 / 4683 / 468 |
| 711 | ice | 0.537 | 3.135 | 2.384 | 0.519 | 0.90 | 65 / 20.3 / 688 / 68 |
| 711 | qem-closest | 0.449 | 15.588 | 2.925 | 0.466 | 1.13 | 63 / 19.7 / 677 / 62 |
| 711 | qem-successive | 0.453 | 15.588 | 2.925 | 0.464 | 1.07 | 62 / 18.4 / 681 / 68 |

### bob
# fine: exact distance 508.3 ms/source, Poisson 44.7 ms, Voronoi (20 sources) 9.95 s, all-pairs (extrapolated) 2665.6 s, exp map 0.00 ms (100 traces)
| n̄ | method | exact geod. L1 % | heat geod. L1 % | Poisson L2 % | all-pairs L1 % | Voronoi wrong % | speedup exact / Poisson / all-pairs / Voronoi |
|---|---|---|---|---|---|---|---|
| 107 | ice | 5.045 | 11.207 | 6.169 | 5.368 | 6.12 | 965 / 347.2 / 5572 / 894 |
| 107 | qem-closest | 3.703 | 12.238 | 17.011 | 4.796 | 5.88 | 644 / 288.5 / 6648 / 611 |
| 107 | qem-successive | 3.635 | 12.228 | 17.048 | 4.735 | 6.34 | 794 / 402.8 / 7285 / 743 |
| 200 | ice | 3.039 | 7.525 | 2.850 | 3.497 | 2.86 | 418 / 175.8 / 3990 / 392 |
| 200 | qem-closest | 2.088 | 8.093 | 9.488 | 3.130 | 3.61 | 341 / 216.0 / 3790 / 265 |
| 200 | qem-successive | 2.099 | 8.038 | 9.505 | 3.111 | 3.65 | 346 / 217.9 / 4173 / 327 |
| 535 | ice | 1.213 | 4.526 | 1.105 | 1.455 | 1.35 | 88 / 65.1 / 751 / 70 |
| 535 | qem-closest | 0.898 | 6.623 | 4.267 | 1.111 | 1.50 | 70 / 76.4 / 677 / 67 |
| 535 | qem-successive | 0.888 | 6.608 | 4.268 | 1.105 | 1.48 | 70 / 75.9 / 662 / 56 |

### bs_rest
# fine: exact distance 1822.3 ms/source, Poisson 40.5 ms, Voronoi (20 sources) 36.85 s, all-pairs (extrapolated) 34838.2 s, exp map 0.01 ms (99 traces)
| n̄ | method | exact geod. L1 % | heat geod. L1 % | Poisson L2 % | all-pairs L1 % | Voronoi wrong % | speedup exact / Poisson / all-pairs / Voronoi |
|---|---|---|---|---|---|---|---|
| 200 | ice | 2.742 | 9.556 | 7.557 | 3.293 | 4.95 | 1808 / 177.3 / 6683 / 1399 |
| 200 | qem-closest | 2.108 | 5.943 | 5.440 | 2.014 | 4.62 | 1883 / 100.3 / 7322 / 1800 |
| 500 | ice | 1.590 | 4.733 | 2.102 | 1.889 | 2.34 | 531 / 80.1 / 5219 / 523 |
| 500 | qem-closest | 1.037 | 4.177 | 2.602 | 2.172 | 2.36 | 587 / 87.4 / 5204 / 565 |

### spot_subdiv
# fine: exact distance 6637.7 ms/source, Poisson 72.4 ms, Voronoi (20 sources) 126.85 s, all-pairs (extrapolated) 156045.7 s, exp map 0.01 ms (100 traces)
| n̄ | method | exact geod. L1 % | heat geod. L1 % | Poisson L2 % | all-pairs L1 % | Voronoi wrong % | speedup exact / Poisson / all-pairs / Voronoi |
|---|---|---|---|---|---|---|---|
| 201 | ice | 1.699 | 4.806 | 11.838 | 2.887 | 3.55 | 5430 / 307.0 / 22935 / 4530 |
| 201 | qem-closest | 1.423 | 6.006 | 11.930 | 2.547 | 2.84 | 5458 / 325.3 / 24995 / 5372 |
| 201 | qem-successive | 1.442 | 6.004 | 11.895 | 2.586 | 3.11 | 6038 / 318.7 / 24228 / 5256 |
| 496 | ice | 0.786 | 4.716 | 6.775 | 1.455 | 1.46 | 1564 / 132.4 / 18411 / 1593 |
| 496 | qem-closest | 0.600 | 4.473 | 15.636 | 1.228 | 1.28 | 1488 / 140.5 / 19377 / 1494 |
| 496 | qem-successive | 0.596 | 4.465 | 15.626 | 1.230 | 1.34 | 1631 / 137.0 / 17505 / 1429 |
| 2476 | ice | 0.172 | 2.366 | 1.275 | 0.292 | 0.31 | 110 / 25.2 / 1012 / 117 |
| 2476 | qem-closest | 0.140 | 9.842 | 2.709 | 0.278 | 0.29 | 97 / 23.7 / 956 / 101 |
| 2476 | qem-successive | 0.140 | 9.840 | 2.709 | 0.278 | 0.28 | 103 / 23.4 / 973 / 103 |

Reading: **Poisson (Fig. 21): ICE is better** in 8 of 11 rows, equal in 1 (spot_subdiv 201), QEM
better in 2 (blub 143: `10.5` against `9.4`; bs_rest 200: `7.6` against `5.4`); where the size is
not tiny by a factor 2–4 (bob 535 `1.1 %` against `4.3 %`, bob 200 `2.9` against `9.5`, spot_subdiv
496 `6.8` against `15.6`, 2476 `1.3` against `2.7`; bs_rest 500 `2.1` against `2.6`).
**Heat-method geodesics: ICE is better** in 8 of 11 rows by 1.1–5× (blub 711 `3.1` against `15.6`,
spot_subdiv 2476 `2.4` against `9.8`; QEM's extrinsic coarse mesh has the wrong lengths for the
diffusion), QEM in bs_rest (both rows) and spot_subdiv 496. **Exact geodesics (Fig. 22): QEM is
better in all 11 rows, by 16–37 %** (blub 711 `0.54` ICE / `0.45` QEM; spot_subdiv 2476 `0.17` /
`0.14`; bob 535 `1.2` / `0.9`; bob 200 `3.0` / `2.1`): the paper's 4× advantage for geodesics
(against the earlier method of Lee et al., not against QEM) is *not reproduced against QEM*.
All pairs: QEM lower in 10 of 11 rows by 5–40 % (bs_rest 500 the exception); Voronoi: mixed,
within ±1 point. The absolute ICE numbers agree with Fig. 23/25 (spot_subdiv 496 vertices:
`0.79 %` exact error at a 1564× speedup; Voronoi 1.5 % wrong at a 1593× speedup; all pairs `1.5 %`
at 18 000×). The exponential-map column is only produced for ICE (median endpoint distance
`0.10–0.33` of the traced length: a loose upper bound, §10.1).

**Fig. 16 scaling (`scaling`).** `simplify` to 10 % of the vertices with a tracked point per
vertex, plus the compaction and `P`, spot 4-1-subdivided (`regular_subdivision`; the first row is
`spot_control_mesh` subdivided once, the others `spot_triangulated` subdivided 0…5 times), one
process per size, single thread, pinned (`taskset -c 63`), sequential. **The machine was loaded**
(another job, up to 24 threads, load average 14–33 during the runs; the early rows overlapped
with this agent's own application jobs): the numbers are *indicative*, to be rerun on a quiet
machine.

| mesh | levels | \|V\| | \|V̄\| | build s | simplify s | P s | total s | s per 10⁶ vertices | load avg (1/5/15 min) at start |
|---|---|---|---|---|---|---|---|---|---|
| spot_control_mesh | 1 | 746 | 75 | 0.001 | 0.077 | 0.000 | 0.078 | 104.6 | 32.39/23.93/13.12 |
| spot_triangulated | 0 | 2930 | 293 | 0.003 | 0.286 | 0.000 | 0.290 | 99.0 | 32.39/23.93/13.12 |
| spot_triangulated | 1 | 11714 | 1171 | 0.012 | 1.132 | 0.001 | 1.145 | 97.7 | 32.39/23.93/13.12 |
| spot_triangulated | 2 | 46850 | 4685 | 0.055 | 4.248 | 0.005 | 4.308 | 92.0 | 32.39/23.93/13.12 |
| spot_triangulated | 3 | 187394 | 18739 | 0.172 | 12.801 | 0.019 | 12.992 | 69.3 | 32.43/24.08/13.22 |
| spot_triangulated | 4 | 749570 | 74957 | 0.735 | 55.935 | 0.131 | 56.801 | 75.8 | 32.60/24.53/13.55 |
| spot_triangulated | 5 | 2998274 | 299827 | 2.849 | 220.828 | 0.591 | 224.269 | 74.8 | 32.09/25.63/14.55 |

Log–log least-squares fit of the total time against `|V|` over the seven sizes (7.5·10² … 3.0·10⁶):
**exponent 0.95** (`t ≈ 1.5·10⁻⁴ |V|^0.95` s; 0.95 for `simplify` alone; 0.947 from 11 714 up;
0.935 up to 187 394): linear to within the noise (the per-vertex time falls from 105 to 75 µs as
the caches stop mattering, as the removal count is `0.9 |V|` here). Throughput at 750k vertices:
674 613 removals in 55.9 s = 12.1·10³ per second (paper: about 10 000 per second on its machine).
Against the percent reduction (spot subdivided 4 times, 749 570 vertices):

| target | \|V̄\| | simplify s | P s | total s | load avg at start |
|---|---|---|---|---|---|
| 90.0 % of \|V\| | 674613 | 9.857 | 0.389 | 10.974 | 21.89/24.70/16.78 |
| 50.0 % of \|V\| | 374785 | 27.805 | 0.243 | 28.744 | 22.26/24.67/16.86 |
| 10.0 % of \|V\| | 74957 | 51.452 | 0.115 | 52.254 | 23.86/24.81/17.15 |
| 1.0 % of \|V\| | 7496 | 55.844 | 0.085 | 56.611 | 21.00/23.94/17.30 |
| 0.1 % of \|V\| | 750 | 56.746 | 0.076 | 57.532 | 16.25/22.12/17.06 |

The time grows with the number of removals (10 s for 10 %, 28 s for 50 %, 51 s for 90 %) and
saturates (57 s at 99.9 %), the queue work being done up front.

---------------------------------------------------------------------------------------------------

## 11. Python and Blender

* **Python** (`ddg_tools`): `simplify(positions, triangles, *, target=…, ratio=…, weights=…,
  memory=True, fixed=…, masses=…, points=…, refine_before=…, refine_after=…, anisotropy=…)` → an
  object with the coarse faces (fine ids), lengths, point locations, `P` and `P^vec` (scipy-free:
  COO triplets, with a helper to build scipy matrices if scipy is present), the report; the
  multigrid solver; the applications of §10.1. Stubs, README, tests.
* **Blender** (the user's choice, Q3): an operator "Intrinsic Simplify" (coarse-mesh view) on the
  active mesh: target (count or percent), curvature/area weight, memory on/off, selected vertices
  fixed, pre/post refinement toggles; result: a face-corner colour attribute colouring every fine
  face by its coarse face (Fig. 1 style), plus the coarse vertices as a vertex group, and the
  report in the info area. No change to the mesh geometry. Tests in `blender/tests/`.

**As built (WP11).** Python (`crates/ddg-py/src/simplification.rs`, glue `src/python/simplification.rs`,
stubs, README, `python-tests/test_intrinsic_simplification.py`): `simplify(positions, triangles, *,
target | ratio, weights=0.0 (w of mixture(w)), memory, fixed, masses + mass_weight, points=(faces,
bary), refine_before, refine_after, repair_simplicial=True, anisotropy=tau + anisotropy_direction
("max" | "min" | 3 reals per vertex) + anisotropy_clamp, representation="integer")` returns a
`Simplification`: `vertices`/`faces` (fine ids), `faces_compact`, `lengths`, `point_face`/`point_bary`/
`point_vertex` for every tracked point (input vertices, then the extra points), `face_colours`/
`point_colours`/`n_colours`, `report` (dict), `prolongation()` and `point_prolongation()` as COO,
`vector_prolongation(n)` as COO with `re`/`im` (refused with refinement or anisotropy: the coarse
frames are not those of the input mesh), `*_scipy()` helpers (CSR; ImportError without scipy and
numpy), and the applications `poisson`, `geodesic(source, "exact" | "heat")`, `voronoi`.
`ratio` is relative to the input's vertex count (also with `refine_before`). `build_hierarchy`
returns a `Hierarchy` (levels, `P_k`, `L_k`, masses, report) whose `multigrid(matrix | shift,
known, ...)` gives a `MultigridSolver` (`solve`, `pcg`, `info`). `qem_simplify(target, ...)` is the
baseline with its closest-point projection. Errors: `InputError` (message starts with the
parameter), `MeshError`, `SolveError`; the GIL is released in every long call. Blender
(`blender/ddg_uv_unwrap/simplify_operator.py`): `object.ddg_intrinsic_simplify`, panel DDG >
Intrinsic Simplify, Object and Mesh menu entries; one centroid per face of Blender's tessellation
is tracked; writes the corner colour attribute `IntrinsicSimplifyFaces` (12-colour palette cycled
over the greedy colouring), the face attribute `IntrinsicSimplifyFaceId` and the vertex group
`IntrinsicSimplifyVertices` (vertices of the input that stay, weight 1); like Show Developability it
does not change the viewport shading. Post-refinement with signposts can fail to locate a point
that lies exactly on an edge (§6.2); the integer representation is the default.

---------------------------------------------------------------------------------------------------

## 12. Numerics and verification

### 12.1 Tolerances

All tolerances are relative to the local scale (mean edge length of the face or star) unless they
are angles: Delaunay and convexity as IT-SPEC; flattening `1e-12 Θ̃`; barycentric inside test by
the signs of `det(p − x_i, p − x_j)` (paper Sec. 6.4) with a relative slack `1e-12` before
clamping; clamps are counted (S16). No absolute `1e-5` as in the reference (R26).

### 12.2 Verification layers

1. **Certificates**: `certify()` (§2.4) after every operation in tests; the bijectivity checks of
   the map (every point in one face, barycentrics valid); curvature conservation; the planar-grid
   test (a flat grid with boundary: the tracked map reproduces positions to `1e-9`; reference
   1.9e-8).
2. **Differential against the reference driver**: a *fixed* build of the driver
   (`../IntrinsicErrorMetric/ref-driver/`, with R8, R10, R18 fixed by copying those files into the
   driver directory and documented there) and a `ReferenceCompat` option set in this crate that
   reproduces the reference's non-bug choices (no initial Delaunay R1, `mixture(w)` R2, boundary
   normalization `2π/(Θ+π)` R3, tentative keys before the frame update R16, 2-ring Delaunay
   seeding R21, the sum over receiving neighbours R11, the `+1e-6` R12, degree-1 boundary removal
   R9). On spot, bob, blub, capsule, crumpled developable: the same removal order and coarse
   lengths to `1e-9` in compat mode; statistics compared in default mode. *(The same order is
   not reachable; see §12.3 for what is reached and why.)*

   *As built (WP3).* `SimplifyOptions::reference_compat(w, target)`: no initial Delaunay (R1);
   `mixture(w)` with the weights folded into the masses (R2); the reference's boundary frames
   `2π/(Θ + π)` in the transport (R3: feasible without touching the mesh, as `A = φ · 2Θ/(Θ + π)`
   from the crate's `φ = (π/Θ) θ`, the virtual direction `π` becoming `2πΘ/(Θ + π)`; the
   reference's ear range `π` is not reproduced); the sum over receiving neighbours only (R11);
   `ReferenceEpsilon` (R12); keys in the frames before flattening and the commit before flattening
   at `i`, after at `j` (R16, R19); `Requeue::OneRing`; and in `ReferenceCompat`: 2-ring Delaunay
   seeding (R21: every edge of the faces within two edge steps of the merged face), the Delaunay
   tolerance `1e-5` (R22, the angle-sum test only; the reference's looser convexity test is not
   reproduced), the flattening tolerance `|K| < 1e-7` (`1e-7/2π` relative), ties by push time then
   vertex (the reference's `(cost, stamp, v)`). **R9 needs no option:** after the flips to degree 2
   the reference flips the last spoke `ia` of the flat star (the quad `i, c, a, b` with the angle
   `π` at `i`), giving the faces `(i, c, b)` (degenerate, deleted) and `(c, a, b)` with
   `ℓ_bc = ℓ_bi + ℓ_ic`; our merge gives the same face `abc`, length and point layout, and the
   flips before it follow the same greedy rule (R9 vs S13 differ only in rounding: confirmed by
   the replays of §12.3, which end in the same connectivity with 48 boundary and 13 ear
   removals). Not reproduced: R4 (the reference halfedge is the face side of least `A`,
   re-chosen after flips and removals; ours is held fixed and passed clockwise), R13 (first CCW
   edge instead of the shortest; only multi-edges), R23 (permanent zero scores), the convexity
   tolerance of R22, R5/R6 (the tentative removal tests the flips; never differs in practice).
3. **Lean oracle** (`lean/oracle`): Float implementations of the flattening step, the transport
   update and the cost of one removal on small stars; golden fixtures checked against the crate.
4. **Lean proofs** (`lean/proofs`, §13).

### 12.3 Against the reference (WP3)

Runs (2026-11-01; logs and tables in `runs_tmp/iem/wp3/`, driver in `runs_tmp/iem/reference-driver.md`):
the fixed reference `iem_ref_fixed` (R8, R10, R18 fixed in copies) against
`examples/simplify.rs --compat --ref DIR --replay`. Three comparisons: the removal order of both
runs; a **replay** that removes the vertices in the reference's order with the compat metric
(public API only) and compares every commit cost and the final lengths, which does not depend on
how near-ties are broken; and default-mode statistics.

| run | removals | first divergence (step) | prefix sets equal (steps) / max diff / at end | key rel. diff. before it | replay: commit cost rel. diff. median / 90% / 99% | replay: final edges matched, max rel. length diff. |
|---|---|---|---|---|---|---|
| spot_500 | 2430 | 4 | 1526 / 4 / 2 | 1.0e-9 | 1.3e-4 / 1.6e-3 / 9.3e-3 | 1494/1494, 1.4e-6 |
| spot_min | 2927 | 4 | 1526 / 36 / 3 | 1.0e-9 | 2.1e-4 / 3.7e-3 / 2.1e-2 | 3/3, 1.1e-8 |
| bob_500 | 4844 | 5 | 1013 / 105 / 95 | 1.5e-9 | stops at step 1501 (Stuck, see below) | – |
| blub_500 | 6606 | 3 | 66 / 214 / 187 | 5.4e-9 | 5.9e-5 / 6.0e-4 / 7.0e-3 | 1493/1494, 1.2e-4 |
| blub_min | 7103 | 3 | 66 / 214 / 3 | 5.4e-9 | 7.0e-5 / 9.5e-4 / 1.8e-2 | 3/3, 2.9e-7 |
| capsule_500 (w 0) | 10158 | 23 | 174 / 154 / 128 | 1.6e-7 | 1.0e-5 / 3.4e-4 / 1.7e-2 | 1494/1494, 1.8e-15 |
| capsule_w0.8 | 10158 | 3 | 1031 / 227 / 193 | 1.8e-11 | 3.4e-5 / 4.9e-4 / 5.5e-2 | 1494/1494, 1.1e-15 |
| crumpled_50 | 239 | 91 | 236 / 1 / 0 | 3.6e-3 | 1.4e-7 / 9.7e-5 / 1.5e-1 | 97/97, 5.6e-9 |
| crumpled_min | 286 | 91 | 258 / 4 / 3 | 3.6e-3 | 3.7e-7 / 8.1e-2 / 4.3e-1 | 3/3, 2.1e-10 |
| spot_subdiv_500 | 24255 | 1 | 619 / 326 / 226 | – | 2.8e-5 / 4.0e-4 / 2.8e-3 | 1494/1494, 3.1e-6 |

Why they differ:
* **The order diverges at near-ties.** Up to the first divergence the keys agree to `1e-9`
  relative (the first removals are almost flat vertices whose masses `K = 2π − Θ ≈ 1e-7` carry
  the rounding of `Θ`: Kahan's angle formula here, `acos` of normalized lengths there). The
  first divergence is always a pair of keys within that noise (spot: `3.07849016166e-7` against
  `3.07849016176e-7`, a mirrored pair; spot_subdiv at step 1, `2.49266761e-10` against
  `2.49266748e-10`); the swapped pair is removed one step later. The prefix sets then stay close
  (spot: never more than 4 vertices apart) but drift on larger meshes as keys diverge.
* **Commit costs in the replay** (same order, same connectivity, same masses): the transport
  itself is exact — on a planar grid with rigid frames the error vectors point to the mass
  centroids to `1e-15` (`tests/simplify.rs`), for interior and boundary vertices. On the same
  grid the reference's costs differ at most steps because R12 makes the shares of a flat vertex
  rounding noise (`Σ|ΔK| ≈ 1e-16` is not exactly 0); with a diagnostic build `iem_ref_flat`
  (S11's `|K_i| ≤ 1e-12` rule in the reference) the grid replays agree to a median of
  `1e-13`–`1e-15` until a corner is removed and curvature appears. On curved meshes the
  remaining differences (median `1e-5`–`2e-4`) come from how the error vectors are carried when
  a neighbour's frame changes: the frames are anchored at a reference halfedge whose `φ` is held
  fixed, and the reference's anchor (R4: the face side of least `A`, re-chosen after flips) is
  not ours (fixed, passed clockwise). The flattening tolerance is not the cause (a replay with
  `1e-14` changes nothing); masses agree.
* **Lengths** agree to `1e-15` where no flattening is inexact (capsule), else to the reference's
  flattening tolerance `|K| < 1e-7` (`1e-6`). blub_500: one edge differs (a Delaunay
  co-circular tie). bob: at replay step 1501 vertex 1160 cannot be flipped to degree 3 here (our
  flip shape test, IT-SPEC N5, refuses a nearly degenerate flip that the reference's convexity
  test `< π + 1e-5` allows, R22).
* crumpled: the keys before the divergence differ by up to `3.6e-3` at four interior vertices
  of the developable sheet, where `K` is close to rounding size (R12's shares) and the vectors
  already carry frames (R4); the order still agrees for 90 steps and the prefix sets end equal.

Default mode against the reference (`mixture(w)` with `w` as given, else the curvature channels;
initial Delaunay; S6; S11), at the same counts — the quality is the reference's:

| run | ours: min angle °, \|K\| median / 90% / max, Σ\|K\| | reference (fixed) |
|---|---|---|
| spot_500 | 9.9, 0.187 / 0.415 / 0.925, 108.27 | 6.9, 0.188 / 0.396 / 0.802, 107.79 |
| bob_500 | 16.8, 0.085 / 0.146 / 0.473, 46.82 | 15.0, 0.085 / 0.151 / 0.364, 46.92 |
| blub_500 | 8.2, 0.090 / 0.309 / 1.306, 70.67 | 9.7, 0.089 / 0.309 / 1.296, 70.74 |
| capsule_500 | 6.6, 0.028 / 0.035 / 0.048, 12.57 | 7.9, 0.028 / 0.035 / 0.045, 12.57 |
| capsule_w0.8 | 17.2, 0.032 / 0.045 / 0.058, 12.57 | 14.1, 0.032 / 0.045 / 0.058, 12.57 |
| crumpled_50 | 4.3, 0.323 / 1.234 / 2.357, 24.67 | 4.3, 0.318 / 1.234 / 2.357, 24.67 |
| spot_subdiv_500 | 7.0, 0.182 / 0.401 / 0.771, 103.60 | 11.3, 0.181 / 0.396 / 0.807, 103.96 |

Distortion (§10.3) is not compared yet (WP4/WP9).

---------------------------------------------------------------------------------------------------

## 13. Formulas and Lean

Formula ids `iem-…`; tiers as IC-SPEC §9: A (algebraic/analytic identities), B (one star/face), C
(whole mesh), cited.

| id | statement | tier |
|---|---|---|
| `iem-flatten-newton` | Eq. 4 is Newton's step for `Θ_i(u) = Θ̃_i` (Springborn et al. Eq. 9, 10 restricted to one vertex): `∂Θ_i/∂u_i = −½ Σ_f Σ_{corners a at i}` (cotangents of the angles opposite the sides of `a` that are not self-edges); for a face with one corner at `i` this is `cot θ^k_ij + cot θ^j_ki`; with a self-edge `ii` the per-corner sum of §3.4 / Eq. 4 is **not** the derivative (the cotangent opposite the self-edge is dropped, as the Rust `face_angle_sum_and_derivative` does); positive for valid faces without a self-face | A (proved) |
| `iem-flatten-monotone` | with extended angles, `Θ_i(u)` is continuous, non-increasing, `→ deg·π` at `−∞` (`deg` = faces of the star), `→ 0` at `+∞` (no self-face); a bracket `Θ(lo) ≥ Θ̃ ≥ Θ(hi)` contains a root; a root exists for `0 < Θ̃ < deg·π`; strictly decreasing between `u₁ < u₂` when some face (not a self-face) is valid at `u₁` or `u₂`; the root is unique if some face is valid at it | A/B (proved) |
| `iem-flatten-conserves-curvature` | flattening `i` changes only `Θ` of `i` and `N_i`, and `Σ_{N_i} K̃ = Σ_{N_i} K + K_i` (local Gauss–Bonnet: each face's angles sum to π) | A (proved) |
| `iem-projective-interpolation` | Eq. 10 maps each face bijectively onto itself, fixes the corners, and agrees on a shared edge for the two faces (C⁰), for any positive per-vertex weights. Note: Eq. 10's `e^{+u}` is the *inverse* of CETM's circumcircle-preserving map (`e^{−u}`, CETM Sec. 3.4, App. D); both are valid bijective tracking maps | A (proved) |
| `iem-degree-three-barycentric` | in a flat degree-3 star, the barycentric coordinates of `i` in the merged face are the normalized sub-triangle areas, and the points of the three faces map affinely | A (proved) |
| `iem-boundary-merge` | a flat boundary vertex of degree 2 lies on the segment `bc` of the merged triangle, at the coordinate given by the layout | A (proved) |
| `iem-flip-to-degree-three` | = `ic-vertex-removal` (Thm. D.1) for simplicial stars (corners `≤ π`; the strict version is false: the cross star); the boundary analogue: a flat half-star of degree `> 2` has an interior spoke with a **strictly** convex quadrilateral; Δ-complexes: `ic-vertex-removal-delta` (open) | B′ (proved) |
| `iem-wrap` | S4: the wrap of an angular coordinate to `[0, 2π)` inside / `[0, π]` on the boundary preserves the direction (`e^{iφ}` inside) and the CCW order of the fan | A (proved) |
| `iem-redistribution-weights` | Eq. 9: `α_ij ≥ 0`, `Σ α = 1`, and the denominator is ≥ `|K_i| > 0` for a non-flat vertex | A (proved) |
| `iem-memoryless-cost` | Eq. 6: `Σ m̃_j‖t̃_j‖ = Σ α_ij m_i ℓ_ij` when `t = 0`, and this is the 1-Wasserstein distance of moving the mass `m_i` at `i` to the neighbours (a single source: every coupling is the forced one; distances along the flat star are `ℓ_ij`) | A (proved) |
| `iem-center-of-mass-update` | Eq. 7/8 in the plane: `t̃_j` points from `x_j` to the combined centre of mass when `t_i`, `t_j` point to the centres of mass of `m_i`, `m_j` | A (proved) |
| `iem-transport-rotation` | `R_ij = e^{i(φ_ji + π − φ_ij)}` maps the frame of `i` to the frame of `j` along the edge, in a flat region (on a curved surface it is the definition, `app-vector-transport`, cited) | A (proved) |
| `iem-vector-interpolation` | Sec. 7.3.1: the `β` angles of the three edges in the face frame and the rotations `z` | A (proved) |
| `iem-anisotropic-scale` | the length factor is in `[1 − τ, 1]` for unit `u`, `ê` | A (proved) |
| `iem-galerkin-spd` | `Pᵀ A P` is SPD for SPD `A` and full-column-rank `P` (also on the free unknowns after Dirichlet elimination); unit rows of the coarse vertices give full column rank | A (proved) |
| `iem-distortion` | the singular values of the affine map between two triangles from their lengths (`‖J‖²_F = Σ cot θ_k ℓ̂_k²/(2A)`, `|det J| =` area ratio, `σ² = (s ± √(s² − 4d²))/2`); `σ₁σ₂ =` area ratio. Khodakovsky et al. Sec. 2 define the measures from the eigenvalues `γ = σ²` of the metric tensor (`γ_max/γ_min`, `γ_max γ_min`); §10.3's `σ₁/σ₂`, `log(σ₁σ₂)` are monotone in them | A (proved) |
| `iem-qem-quadric` | the quadric of a plane gives the squared distance (weighted: `w (n·(x − p))²`); quadrics add; a solution of `A v = −b` (the cofactor formula of `optimal_point`) minimizes a quadric with symmetric PSD `A` (uniquely for PD); along a segment the clamped vertex `clamp(−β/α, 0, 1)` minimizes it for `α > 0` | A (proved) |
| `iem-wasserstein-memory` | the memory-based cost approximates the transport cost (paper Sec. 5.2: "only approximate") | cited |
| `iem-cetm-existence` | Springborn et al. 2008 Thm. (convex energy, uniqueness up to constants) | cited |
| `iem-delaunay-quality` | 30°/120° corner bound after post-refinement (Gillespie et al. 2021a) | cited (`ic-refinement-termination`) |

Lean (WP10): `lean/proofs/DdgProofs/IntrinsicSimplification/*.lean`; the exact statements and theorem
names are in `lean/proofs/README.md`, section "Intrinsic Simplification"; the cited ids in
`lean/proofs/CITED.md`. Oracle: `ddg-oracle iem` (tasks `step`, `flatten`, `weights`, `removal` on one
star), fixtures `lean/oracle/fixtures/iem/` (format in `lean/oracle/README.md`); its side `k` is
**opposite** corner `k`, i.e. Rust's `l[(k + 1) % 3]`.

---------------------------------------------------------------------------------------------------

## 14. Decisions

* **S1 Own data structure.** `CoarseMesh` = `Complex` + lengths + angular coordinates + channels +
  tracked points; not a `Triangulation` (the metric changes). The cotan Laplacian, masses and the
  vector operators are computed from the lengths with the shared formulas.
* **S2 Stable fine ids** (§2.1).
* **S3 Initial Delaunay** on (paper Alg. 1 line 3; R1 off in the reference).
* **S4 Angular coordinates**: normalization `2π/Θ` inside, `π/Θ` on the boundary (IT-SPEC N41;
  R3 used `2π/(Θ+π)`); a reference halfedge per vertex whose `φ` is held fixed under length
  changes, passed to its CW neighbour when flipped or removed; every value wrapped to `[0, 2π)`
  (inside) and kept in `[0, π]` (boundary) by one function (`iem-wrap`), which fixes R18; after a
  flattening the vertex's own coordinates are recomputed too (R19).
* **S5 Tentative = commit**: one code path on an undo log; the key is the commit cost for an
  unchanged neighbourhood; tests assert it (R16, R5).
* **S6 Re-queue set**: `N_i` (paper) plus the endpoints of every edge flipped by the Delaunay pass
  after the removal (their stars changed). Option `requeue = OneRing` for the paper's literal set.
  As built: the four vertices (old and new end points) of every flip of the removal (§5).
* **S7 Failures revert and skip** (paper; R6, R17): never a panic or an abort.
* **S8 Stop at the first `∞`** (paper, R17).
* **S9 Flattening**: bracketed Newton on the extended angles (§3.4), Springborn's flip fallback
  (R7, R8).
* **S10 Default channels**: curvature `K⁺`, `K⁻` weight 1; `mixture(w)` as the reference (R2).
* **S11 Flat vertex** `|K_i| ≤ 1e-12` → equal shares (R12); the cost sums over all of `N_i` (R11).
* **S12 Neighbours** are distinct vertices `≠ i`; transport along the shortest connecting edge
  (R13). No minimum vertex count beyond what the operations allow (§3.5).
* **S13 Boundary removal** by degree 2 and a merge (paper; R9); ears by a flip first.
* **S14 Flip order to degree 3/2**: self-edges first, then the largest opposite-angle sum among
  flippable spokes, re-scored every round (R23).
* **S15 Local Delaunay seeding**: the edges with both endpoints in `N_i` (paper; R21).
* **S16 Point relocation**: planar layout of the flat region, inside test by orientation signs with
  relative slack, tiny negatives clamped and counted, a point outside by more than the slack is an
  internal error reported with its id (never silently placed; R26); points on a new edge go to the
  face with the lower id (deterministic).
* **S17 Points are known in advance** (paper Sec. 7.1); no replay.
* **S18 `P^vec` literal** (unnormalization; R28); `ψ` fallbacks never leave an empty row (R29, R31).
* **S19 Pre-refinement on integer coordinates** by default (Q3 of IC-SPEC), signposts an option.
* **S20 Exact geodesics on Δ-complexes**: `exact_distance` generalized (WP6).
* **S21 Exponential map**: on the coarse mesh only (the map is fine → coarse).
* **S22 Principal curvatures from Panozzo et al. 2010**, implemented from the paper.
* **S23 Extrinsic baseline**: QEM (Garland–Heckbert) with closest-point and Liu et al. 2021
  successive maps (§10.4); not in the reference (R43).
* **S24 Hierarchy from one run with snapshots** (masses carry across levels); the reference's fresh
  per-level runs (R33) as an option.
* **S25 The map through a flattening: Springborn et al.'s `e^{−u}` (default), Eq. 10 as printed
  an option.** The paper cites the circumcircle-preserving projective interpolation of Springborn
  et al. 2008 (Sec. 3.4, App. D: linear interpolation of `(x, 1)` and `e^{−u}(s, 1)`), but prints
  Eq. 10 with `e^{+u}` (and the reference implements that): its inverse, under the same length
  convention (Eq. 3). Found by WP10. Both are continuous bijections of each face
  (`iem-projective-interpolation`), so tracking is correct either way; the coarse mesh and the
  removal order do not depend on it. WP9a measured the distortion (§10.5): `Negative` has the lower
  anisotropic distortion in 8 of 8 configurations (mean and maximum), `Positive` a 3–19 % lower area
  distortion in 8 of 8. **Decision (the user's, 2026-11-01): `ProjectiveSign::Negative` (`e^{−u}`) is the
  default, `Positive` (`e^{+u}`) an option** (`CoarseMeshOptions::projective_sign`; Python
  `projective_sign="negative"|"positive"`; Blender "Point Map"); `ReferenceCompat` comparisons use
  `Positive`. The circumcircle property of
  `e^{−u}` is checked numerically (oracle self-test), not proved.
* **R55** (new, found by WP10): the reference's Newton denominator sums `cot` per corner also in
  faces with a self-edge, which is not the derivative there; S9 drops the cotangent opposite the
  self-edge. The root is the same.
* **Reference bugs not reproduced**: R8, R10, R14, R18, R25, R26, R30, R31, R35, R36, R39, R40,
  R41, R42, R44 (each fixed in the fixed driver build where it affects a comparison, §12.2).

---------------------------------------------------------------------------------------------------

## 15. Work packages

| WP | content | depends on | model |
|---|---|---|---|
| WP0 | `Complex::remove_boundary_vertex` (degree 2), IT-SPEC N54, tests (ids, undo, Δ-cases) | — | with WP1 |
| WP1 | crate skeleton: `CoarseMesh` (§2), construction (§3.1 `from_mesh`), angular coordinates (S4), flips with point tracking (§3.2, S16), Delaunay flipping (§3.3), `certify` (§2.4), errors, report types; copy the reference meshes into `data/intrinsic-simplification/` | WP0 | Opus |
| WP2 | flattening (§3.4, S9), removals (§3.5, interior, boundary, ear; S13, S14), point tracking through them, the undo log and tentative removal (§3.6, S5) | WP1 | Opus |
| WP3 | the ICE metric (§4) and Alg. 1 (§5): channels, weights, memory on/off, fixed vertices, targets, the queue, snapshots, the report; `ReferenceCompat`; the fixed reference driver and the differential tests (§12.2); performance | WP2 | Opus |
| WP4 | the map and prolongation (§7): point locations, `P`, `P^vec` (S18), normalize/unnormalize, colouring | WP3 | Sonnet |
| WP5 | pre/post refinement (§6) on integer coordinates and signposts | WP3 | Sonnet |
| WP6 | applications (§10.1): Poisson, heat and exact distance (exact generalized to `Complex`, S20; this part can start any time), all-pairs, Voronoi, exp map | WP4 (exact: —) | Sonnet |
| WP7 | hierarchy, multigrid solver, cMCF (§9) | WP4 | Sonnet |
| WP8 | anisotropic coarsening and principal curvatures (§8) | WP3 | Sonnet |
| WP9 | evaluation (§10.2–10.4): Thingi10k, timings, distortion, QEM, the extrinsic maps, figures of record | WP4–WP7 (QEM: —) | Sonnet |
| WP10 | Lean (§13, tiers A/B) and the oracle with fixtures | this file | Opus |
| WP11 | Python and Blender (§11) | WP4–WP8 | Sonnet |

At most three agents at a time; timings alone. Contracts are edited by the agent that changes the
behaviour, in its own entries; `ROADMAP.md` by the orchestrator.

---------------------------------------------------------------------------------------------------

## 16. The user's decisions (2026-11-01)

* **Q1 Reference**: fetch, read and build the authors' MIT code (done: clone, map, headless driver).
* **Q2 Scope**: the core (Alg. 1, the ICE metric with channels and weights, fixed vertices and
  masses, pre/post refinement, the map, scalar and vector prolongation) plus the applications
  (Sec. 8.5), the multigrid hierarchy and solver (Sec. 8.7, with cMCF), anisotropic coarsening
  (Sec. 8.6.2) and the evaluation with baselines (Sec. 8.2–8.4, QEM).
* **Q3 Blender**: a coarse-mesh view (§11); no acceleration of other operators.
