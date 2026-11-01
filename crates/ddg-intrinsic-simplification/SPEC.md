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
come with WP3).

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
of the star faces moved by Eq. 10, `b̃ = (e^{u_i} b_i, b_j, b_k)/(e^{u_i} b_i + b_j + b_k)` with the
factor applied to every corner equal to `i` (`iem-projective-interpolation`: continuous across
edges and a bijection of each face onto itself).

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
lengths are a valid Δ-complex input of `ddg-integer-coords` / `ddg-intrinsic`), refining it, and
composing the maps: a tracked point `Face(f, b)` of the coarse mesh is a point of that input and is
located in the refined triangulation by its `T⁰ → T¹` query. The result is again a `CoarseMesh`
(lengths, angular coordinates rebuilt, inserted vertices get new ids `≥` all previous ids).

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

### 10.2 Benchmark (Sec. 8.2; Fig. 16, 17)

Thingi10k (`data/thingi10k/`, the 7745 manifold meshes of IC-SPEC Q1), preprocessed by intrinsic
Delaunay refinement at 25°, coarsened as far as possible (target 0, stop at `∞`), the ratio
`|V̄|/|V|` normalized by the number of components. Report: the share reaching ≤ 10 % and ≤ 1 %
(paper: 98 % and 84 %), failures by reason (paper: 0.1 % from floating point), time per removal,
peak memory. Fig. 16: time against `|V|` on a subdivided cow/spot (`10³…10⁶`) and against the
percent reduction on a 750k-vertex subdivision. Runs in parallel per mesh (64 cores; timings alone).

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
   lengths to `1e-9` in compat mode; statistics compared in default mode.
3. **Lean oracle** (`lean/oracle`): Float implementations of the flattening step, the transport
   update and the cost of one removal on small stars; golden fixtures checked against the crate.
4. **Lean proofs** (`lean/proofs`, §13).

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
* **S25 Eq. 10 as printed** (`b̃_i ∝ e^{+u} b_i`): the paper's and the reference's map, which is
  the *inverse* of Springborn et al.'s circumcircle-preserving projective interpolation (`e^{−u}`,
  CETM Sec. 3.4, App. D; found by WP10). Both are continuous bijections of each face
  (`iem-projective-interpolation`), so tracking is correct either way; the default follows the
  paper, and WP9 measures the distortion (§10.3) with both (an option `projective_sign`).
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
