# ddg-uv

One operation: a polygon mesh goes in, a UV map comes out — one UV per polygon corner, the
chart of every polygon, the cones, the edges that were cut, and a report.

The crate ties the workspace together. It contains no flattening or optimization method of
its own; it places cone singularities with

> Yousuf Soliman, Dejan Slepčev, Keenan Crane.
> *Optimal Cone Singularities for Conformal Flattening.*
> ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201367>

(crate `ddg-optimal-cones`) and flattens through those cones with

> Rohan Sawhney, Keenan Crane. *Boundary First Flattening.*
> ACM Transactions on Graphics 37(1), 2017. <https://doi.org/10.1145/3132705>

(crate `ddg-bff`, a port of the authors' MIT-licensed reference implementation, which also
plans and makes the cuts). The intrinsic mollification used for degenerate triangles is from
N. Sharp and K. Crane, *A Laplacian for Nonmanifold Triangle Meshes*, SGP 2020, through
`ddg-bff`. What is new here: validation in the caller's numbering, polygons with a supplied
triangulation, seams, steering, policies for what the papers leave open, packing, the report.

`ddg-core` has no triangulator, which is why the triangulation of the polygons is an input and
not computed here (Blender triangulates; `PolygonMesh::from_polygons` fans convex polygons).

Contents: [example](#example) · [API](#api) · [stages and policies](#stages-and-policies) ·
[the four decisions](#the-four-decisions) · [the dataset](#the-dataset) · [the boundary](#the-boundary) ·
[the cut flow on polygon meshes](#the-cut-flow-on-polygon-meshes) · [benchmarks](#benchmarks) · [tests](#tests) · [not done, not verified](#not-done-not-verified)

## Example

```rust
use ddg_uv::primitives::grid_cube;
use ddg_uv::{MeshInput, UnwrapOptions, unwrap};

let mesh = grid_cube(4, 2.0);                 // any polygon mesh with a triangulation
let quad = mesh.polygon(0);
let seams = [[quad[0], quad[1]]];             // seams are vertex pairs
let input = MeshInput { seam_edges: &seams, ..mesh.input() };

let map = unwrap(&input, &UnwrapOptions::default())?;
// map.uvs[c]            UV of polygon corner c, packed into the unit square
// map.polygon_charts[p] chart of polygon p
// map.cones             vertex and angle of every cone
// map.cut_edges         [min, max] vertex pairs, seams included
// map.report            timings, cone statistics, certificate, distortion, overlap, warnings
```

## API

| item | what it is |
|------|------------|
| `unwrap(&MeshInput, &UnwrapOptions) -> Result<UvMap, UvError>` | the whole pipeline |
| `MeshInput { positions, corner_vertices, polygon_offsets, triangles, seam_edges, pinned_vertices, prescribed_cones, forbidden_vertices }` | borrowed input; `triangles` are triples of polygon-**corner** indices; `MeshInput::new`, `validate` |
| `UvMap { uvs, polygon_charts, cones, cut_edges, report }` | the result |
| `UnwrapOptions { lambda, cone_signs, max_cone_angle, closed_surfaces, merge, orientation, degenerate, packing, overlap_samples_per_triangle, reoptimize_angles, boundary, require_certified_cones, spike_check }` | every policy, see below; `validate` |
| `ConeSigns`, `ClosedSurfaces`, `SlitPlacement`, `MergeClusters`, `Orientation`, `Degenerate`, `PackOptions`, `SpikeCheck`, `SpikeAction`, `BoundaryValues` | the policy types |
| `prepare`, `solve_cones`, `flatten`, `assemble` (`Prepared`, `Piece`, `ConeStage`, `PieceCones`, `FlatStage`, `PieceFlat`) | the stages as separate functions |
| `pack_charts(&mut uvs, polygon_offsets, polygon_charts, &PackOptions) -> Result<PackReport, UvError>` | the packer on its own |
| `area_distortion(positions, corner_vertices, triangles, uvs) -> AreaDistortion` | the dataset README's distortion measure of any UV map |
| `Report` (`ValidationReport`, `SeamReport`, `DegenerateReport`, `SteeringReport`, `ChartReport`, `ConeSolveReport`, `FlattenReport`, `ConeStats`, `SeamLengths`, `ClosedHandling`, `MergeOutcome`, `Timings`, `PackReport`, `SpikeReport`), `Report::summary() -> Vec<(&str, Value)>` | the report, typed and as a flat key–value list for bindings |
| `UvError`, `UvError::kind() -> ErrorKind { Input, Mesh, Solve }` | every refusal, in the caller's vertex and polygon numbers |
| `PolygonMesh` (`from_triangles`, `from_polygons`, `input`, `edges`), `fan_triangulation` | an owned mesh, for tests and callers without their own arrays |
| `primitives::{uv_sphere, capsule, cylinder, cone, torus, cube, grid_cube, holed_slab, voxel_surface, revolution, jittered}` | Blender-like polygon meshes for tests and benchmarks |

## Stages and policies

**1. Validation.** Array shapes, index ranges, finite numbers, that the triangles triangulate
the polygons, that no vertex is steered twice (`ErrorKind::Input`). Then the surface
(`ErrorKind::Mesh`): a polygon that repeats a vertex, duplicate polygons, an edge with more
than two polygons (`edge (0, 1) is shared by 3 faces: 0, 1, 2`), a diagonal of the supplied
triangulation that coincides with another edge, a vertex whose polygons form several fans, a
surface that cannot be oriented.

- *Loose vertices* are ignored and counted (`validation.n_loose_vertices`).
- *Inconsistent winding*: `Orientation::Repair` (default) flips the minority of every
  connected component, counts the flipped polygons and warns; their UVs are mirrored with
  respect to their own winding. `Orientation::Reject` is an error naming the edge and the two
  polygons.
- *Triangles without area*: `Degenerate::MollifyIfNeeded { relative_delta: 1e-6 }` (default)
  leaves a mesh without such triangles untouched; otherwise all edge lengths grow by the
  smallest common `ε` that restores the triangle inequalities with a margin
  (`ddg_bff::mollify_edge_lengths`). `ε`, the margin, the number of triangles and the worst
  polygon are in `validation.degenerate`, with a warning. `MollifyAlways` is what the reference
  implementation of Boundary First Flattening does; `Reject` is an error naming the polygon.

**2. Seams first.** Seams are cut before anything else (`ddg_bff::cut_mesh`), so the cone
solver and the flattener see the same seam-cut surface: seams are boundary, `u = 0` along
them, both sides keep their length (with the default `BoundaryValues::Zero`; with
`Optimized` they do not, see [the boundary](#the-boundary)). The connected pieces are the charts, numbered by their
lowest polygon.

- A seam that is not an edge of a polygon (a diagonal, two unrelated vertices, a loose vertex)
  is ignored, counted and warned about. A seam on the boundary of the mesh needs no cut.
- A seam of a single edge between two interior vertices that no other seam touches has no two
  sides in an indexed mesh. It is extended by the shortest neighbouring polygon edge, counted
  (`n_extension_edges`) and warned about.

**3. Spike check, cones, then flattening, per piece.** Before the cone solve, the
[spike check](#the-spike-check) looks at the piece without cones. Steering is translated at the vertices that are
interior to their piece: pinned → cone cost 0 (free cone site), forbidden → cost `+∞`,
prescribed → fixed angle. On a boundary, seam or slit there are no cones: pinned and forbidden
vertices have no effect there (counted), a prescribed cone is listed in
`steering.prescribed_not_applied` with a warning. `λ` applies to every piece with its own area
normalized to 1.

| option | default | alternatives |
|--------|---------|--------------|
| `lambda` | `0.1` | any `λ ≥ 0` |
| `cone_signs` | `PositiveWhereFeasible`: `µ ≥ 0` on pieces with boundary and closed pieces of genus 0 and 1; any sign, with a warning, on closed pieces of genus ≥ 2 | `PositiveOnly` (genus ≥ 2 closed is an error), `Any` |
| `max_cone_angle` | `0.95 · 2π`, the upper bound of every cone angle | up to `2π` |
| `closed_surfaces` (closed pieces of genus 0) | `Auto { max_cones: 32, max_cone_fraction: 0.01, fallback_sites: 8 }`: the raw optimum, unless it is dense | `RawOptimum`, `FreeSites { count, max_angle, forbid_elsewhere }`, `Slit { placement, relative_length }`, `Refuse` (any closed piece is an error) |
| `merge` | `Guarded { max_cluster_size: ∞, max_factor: 10, max_increase: 0.05 }` | `Never`, `Always { max_cluster_size }` |
| `reoptimize_angles` | `false` | `true`: a second solve for the best angles at the chosen positions |
| `boundary` | `Zero`: `u = 0` on the boundary of every piece, boundary lengths (and so the two sides of a seam) preserved | `Optimized`: the solver chooses `u` on the boundary, [see below](#the-boundary) |
| `spike_check.action` | `ExcludeFromDistortion`: spike tips get zero distortion weight, see [the spike check](#the-spike-check) | `WarnOnly` |
| `require_certified_cones` | `false`: uncertified cones are used with a warning | `true`: an error |

Then the piece is cut through its cones by `ddg-bff`'s planner with a mask that allows
**polygon edges only**, and flattened with the solver's scale factors. Every piece comes out
as one chart. At the end the pipeline checks that all triangle corners of one polygon corner
received the same UV, that every UV is finite, that every cut edge is a polygon edge and that
every seam is in the cut; a violation is `UvError::Internal`.

**4. Packing.** `PackOptions { margin: 1/256, rotate: true }`: every chart is rotated into its
smallest bounding box (rigid, never mirrored), all charts get **one** scale, boxes are shelved
tallest first, the scale is found by bisection. `packing: None` leaves the charts where the
flattener put them.

**Warnings** are raised, besides the cases above, when thin spikes dominate the distortion
(next section), when cones sit at `max_cone_angle`, when
more than 1 % of the triangles of a chart are mirrored, and when more than 1 % of a chart's
footprint is covered twice.

## The four decisions

Everything in this section is measured with `examples/evaluate.rs`; the logs of the runs
quoted are `runs_tmp/uv-eval-*.log`. "Distortion" is the measure of the dataset's README,
`Σ A u²` with `u = ½ ln(B/A)` per triangle, with the areas normalized to sum 1 and minimized
over the scale of the whole map (`AreaDistortion::scale_free`); see `src/measure.rs` for why
the scale has to be taken out. It is about twice the solver's `½ Σ m u²`.

### Closed surfaces without seams

`evaluate closed`: positive cones, merging without size limit and with the stricter guard the
example hard-codes (factor 2, +0.01; because of it two cells, cube and cone primitive with the
slit at `λ = 0.02`, show unmerged optima). Cells are
cones / distortion / overlap % / mirrored triangles / cut edges / fill of the unit square.

| mesh | (a) raw optimum, merged | (b) 1 free site, `µ ≤ 0.95·2π` | (b) 8 free sites, none elsewhere, `µ ≤ π` | (c) slit, full diameter, λ 0.5 | (c) same, λ 0.1 | (c) same, λ 0.02 | (c) slit of 2 edges at the largest curvature, λ 0.5 |
|---|---|---|---|---|---|---|---|
| UV sphere 32×16 | 482 / 0.0000 / 0 / 0 / 481 / 0.18 | 110 / 0.310 / 11.4 / 101 / 119 / 0.18 | 8 / 0.0107 / 0 / 0 / 54 / 0.34 | 0 / 0.216 / 0 / 0 / 16 / 0.53 | 1 / 0.0665 / 0 / 0 / 24 / 0.43 | 89 / 0.0129 / 0 / 0 / 110 / 0.32 | 21 / 0.138 / 0 / 0 / 35 / 0.48 |
| icosphere, 320 faces | 162 / 0.0000 / 0 / 0 / 161 / 0.29 | 53 / 0.172 / 12.0 / 24 / 58 / 0.27 | 8 / 0.0110 / 0 / 0 / 25 / 0.29 | 0 / 0.207 / 0 / 0 / 10 / 0.52 | 1 / 0.0659 / 0 / 0 / 15 / 0.43 | 37 / 0.0117 / 0 / 0 / 50 / 0.29 | 7 / 0.120 / 0 / 0 / 16 / 0.74 |
| cube | 8 / 0.0000 / 0 / 0 / 7 / 0.37 | 5 / 0.134 / 0 / 0 / 5 / 0.30 | 8 / 0.0000 / 0 / 0 / 7 / 0.37 | 0 / 0.148 / 0 / 0 / 3 / 0.72 | 1 / 0.132 / 0 / 0 / 4 / 0.66 | 4 / 0.0057 / 0 / 0 / 7 / 0.36 | 0 / 0.279 / 0 / 0 / 2 / 0.90 |
| cube, 8×8 quads per side | 8 / 0.0000 / 0 / 0 / 56 / 0.37 | 3 / 0.267 / 0.18 / 10 / 24 / 0.17 | 8 / 0.0248 / 0 / 0 / 53 / 0.39 | 0 / 0.197 / 0 / 0 / 24 / 0.42 | 1 / 0.123 / 0 / 0 / 36 / 0.37 | 5 / 0.0139 / 0 / 0 / 54 / 0.28 | 23 / 0.170 / 0 / 0 / 42 / 0.49 |
| cylinder 32, n-gon caps | 64 / 0.0000 / 0 / 0 / 63 / 0.36 | 22 / 0.212 / 0.61 / 6 / 29 / 0.30 | 8 / 0.0035 / 0 / 0 / 47 / 0.44 | 0 / 0.149 / 0 / 0 / 17 / 0.44 | 1 / 0.0901 / 0 / 0 / 19 / 0.49 | 3 / 0.0205 / 0 / 0 / 32 / 0.42 | 1 / 0.216 / 0 / 1 / 18 / 0.34 |
| cone 32, n-gon base | 33 / 0.0000 / 0 / 0 / 32 / 0.26 | 17 / 0.072 / 0 / 0 / 24 / 0.40 | 8 / 0.0054 / 0 / 0 / 27 / 0.28 | 1 / 0.627 / 0 / 0 / 17 / 0.74 | 1 / 0.0759 / 0 / 0 / 17 / 0.39 | 3 / 0.0177 / 0 / 0 / 25 / 0.30 | 0 / 0.159 / 0 / 0 / 2 / 0.11 |
| capsule 32×8 | 514 / 0.0000 / 0 / 0 / 513 / 0.24 | 172 / 0.266 / 15.1 / 106 / 183 / 0.22 | 8 / 0.0072 / 0 / 0 / 63 / 0.55 | 0 / 0.270 / 0 / 0 / 17 / 0.60 | 1 / 0.0703 / 0 / 0 / 29 / 0.53 | 1 / 0.0436 / 0 / 0 / 28 / 0.47 | 0 / 0.428 / 0 / 0 / 2 / 0.11 |
| axes (7 466 triangles) | 4 / 0.807 / 0 / 0 / 180 / 0.29 | 5 / 1.241 / 3.15 / 5216 / 144 / 0.04 | 4 / 0.810 / 0 / 0 / 177 / 0.30 | 2 / 1.399 / 0 / 0 / 142 / 0.32 | 2 / 0.408 / 0 / 0 / 142 / 0.39 | 2 / 0.343 / 0 / 0 / 142 / 0.46 | 4 / 0.779 / 0 / 0 / 145 / 0.55 |
| bunny (28 576) | 15 / 0.224 / 0 / 0 / 300 / 0.69 | 4 / 0.461 / 0 / 14 / 191 / 0.32 | 6 / 0.258 / 0 / 0 / 279 / 0.44 | 0 / 1.153 / 0 / 0 / 116 / 0.59 | 2 / 0.188 / 0 / 0 / 216 / 0.45 | 8 / 0.0949 / 0 / 0 / 264 / 0.61 | 4 / 0.351 / 0 / 0 / 224 / 0.38 |
| chair (11 448) | 7 / 0.0245 / 0 / 0 / 165 / 0.55 | 5 / 0.108 / 0 / 0 / 118 / 0.32 | 5 / 0.0629 / 0 / 0 / 136 / 0.49 | 0 / 0.249 / 0 / 0 / 66 / 0.59 | 2 / 0.127 / 0 / 0 / 116 / 0.56 | 4 / 0.0395 / 0 / 0 / 172 / 0.31 | 11 / 0.262 / 0 / 0 / 129 / 0.49 |

Also measured: one free site **without** a bound fails on all ten meshes at `λ = 0.5` (the
site takes `4π`; cone primitive 9.04, bunny 11.25) and on six of ten at `λ = 0.05` (6.8 to
8.2), because no cone has `2π` or more. The torus 48×12 (genus 1, `µ ≥ 0`, so no cone at all) is the same under every
policy: 0 cones, distortion 0.0293, 60 cut edges, fill 0.25. (d), refusing, is
`ClosedSurfaces::Refuse`.

What the table says:

- **(a)** is exact where the cones are the corners (cube, subdivided cube: 8 cones of `π/2`)
  and on developable polyhedra (the cylinder and cone primitives: a cone at every vertex is
  their unfolding). It is good on the bumpy closed meshes once merged (axes 14 → 4, bunny
  254 → 15, chair 73 → 7 cones). On smooth round shapes it is a cone at every vertex: 482, 162
  and 514 cones and as many cuts.
- **(b)** with one site is not usable: unbounded it asks for a cone of `2π` or more; bounded
  it gives 17 to 172 cones on the smooth primitives, with up to 15 % overlap and mirrored
  triangles. With **eight sites and cones forbidden elsewhere** it gives at most 8 cones, no
  overlap, no mirrored triangle, and distortion 0.0035 to 0.025 on the primitives. The sites
  come from farthest-point sampling; their positions are not optimized.
- **(c)** With `λ = 0.5` a slit gives 0 or 1 cone and distortion 0.15 to 0.63 on the
  primitives. Reaching 0.01 to 0.02 takes `λ = 0.02`, and then 37 to 89 cones on the spheres.
  The reason is visible in the solver's own number: with `u = 0` on a short slit its objective
  is dominated by a common scale offset of the rest of the surface (UV sphere, slit of half
  the diameter, `λ = 0.5`: `½ Σ m u² = 1.76`, distortion after removing the scale 0.34), which
  packing removes anyway. Cones are spent on something that does not matter for a UV map.

**Default: `Auto`.** The raw optimum is used unless it is dense — more than 32 cones *and*
more than 1 % of the vertices after merging — in which case the piece is solved again with 8
free sites and cones forbidden elsewhere. What it costs:

- a second solve when the fallback triggers;
- fallback cones whose positions are not optimized, always 8 or fewer;
- the exact unfolding of a developable primitive with more than 32 vertices is replaced by an
  approximate map (cylinder 0.0035, cone 0.0054 instead of 0) unless seams are marked;
- the criterion is not calibrated between the cases measured: dense optima had a cone at
  every vertex, sparse ones at most 0.1 % of the vertices (brain: 55 of 152 893).

Closed pieces of genus 1 get no cone with positive cones (a flat torus cut along two loops).
Closed pieces of genus ≥ 2 are the mirror image of the sphere: the optimum is a diffuse
negative measure and `λ` steers nothing down to 0.01 (a slab with two holes, 800 quads: 62
negative cones in one cluster, merged into a single cone of `−4π`, distortion 0.308, 2 mirrored
triangles; three holes: 132 negative cones, merge rejected by the guard, distortion 0.098,
0.44 % overlap). The pipeline returns a map with a warning; for such surfaces seams are the
answer.

### Negative cones

`evaluate signs`: cones (negative ones) / distortion / overlap %.

| mesh | `µ ≥ 0`, λ 0.5 | any, λ 0.5 | `µ ≥ 0`, λ 0.1 | any, λ 0.1 | `µ ≥ 0`, λ 0.02 | any, λ 0.02 |
|---|---|---|---|---|---|---|
| airplane | 4 (0) / 2.383 / 0 | same | 4 (0) / 1.664 / 0 | same | 4 (0) / 1.572 / 0 | 19 (13) / 0.273 / 34.6 |
| body | 4 (0) / 2.142 / 0 | same | 8 (0) / 1.324 / 0 | same | 8 (0) / 1.170 / 0 | 15 (6) / 0.476 / 11.5 |
| hand | 3 (0) / 2.292 / 0 | same | 5 (0) / 1.465 / 0 | same | 5 (0) / 1.314 / 0 | 7 (2) / 0.417 / 0.36 |
| cow | 2 (0) / 1.598 / 0 | same | 5 (0) / 0.655 / 0 | same | 5 (0) / 0.537 / 0 | 6 (1) / 0.531 / 0 |
| bump, camel, chair, mannequin, rhino | | same as `µ ≥ 0` | | same | | same |
| torus 48×12, tube radius 0.25 | 0 / 0.0293 / 0 | 0 / 0.0293 / 0 | same | same | same | same |
| torus 48×12, tube radius 0.6 | 0 / 0.1571 / 0 | 0 / 0.1571 / 0 | same | same | same | same |

At `λ ≥ 0.1` the optimum has no negative cone on any mesh measured, so the bound costs
nothing. At `λ = 0.02` negative cones appear on four meshes: they lower the distortion by a
factor 2.5 to 5.8 (cow: 1 %) and bring up to 35 % overlap. On the two tori the unrestricted
optimum has no cone either.

**Default: `PositiveWhereFeasible`.** Pieces with boundary and closed pieces of genus 0 and 1
get `µ ≥ 0`; closed pieces of genus ≥ 2, where the angles must sum to `2πχ < 0`, get any sign
and a warning. `PositiveOnly` refuses those; `Any` lifts the bound everywhere.

### Cluster merging

`evaluate merge`: 21 dataset meshes × `λ ∈ {0.5, 0.02}` × both sign policies = 84 solves of
the raw problem, merging without size limit. The solver's `½ Σ m u²` before and after:

- 59 solves merged something. The ratio after/before is 1.000 to 1.022 on 54 of them.
- airplane, any sign, `λ = 0.02`: 83 → 19 cones, ratio 1.080.
  camelshotnoise, any sign: 452 → 51 cones, 3.493 → 3.847 (1.101) at `λ = 0.5`;
  618 → 64 cones, 0.0898 → 0.1449 (1.614) at `λ = 0.02`.
- mannequin, `λ = 0.02`, both sign policies: 519 cones in two clusters, 0.01505 → 0.03946
  (2.622), the largest ratio measured.
- The closed meshes: axes 14 → 4 (1.001), bunny 254 → 15 (1.005), chair 73 → 7 (1.022).

The failure the guard was asked for — `camelshotnoise` going to `10⁴` — no longer occurs: the
solver now merges only cones of one sign. The guard stays, as a comparison of the solver's
distortion before and after merging.

**Default: `Guarded { max_factor: 10, max_increase: 0.05 }`, no size limit.** It accepts
everything measured and would have rejected the old failure (a factor 3000). When it rejects,
the unmerged optimum is used and the report warns. A stricter guard is a worse default: with
factor 2 the mannequin at `λ = 0.02` keeps 519 cones instead of 2.

A merged cone can be large: on the mannequin at `λ = 0.1` a ring of 232 small cones becomes
one cone of 5.84 (the bound is 5.97). The solver's `max_merge_cluster_diameter` would prevent
that; the pipeline does not set it, for lack of data on a good value.

### λ

`evaluate lambda`: positive cones; cones / distortion.

| mesh | λ 2 | 1 | 0.5 | 0.2 | 0.1 | 0.05 | 0.02 | 0.01 | 0.002 |
|---|---|---|---|---|---|---|---|---|---|
| airplane | 2 / 4.256 | 2 / 2.945 | 4 / 2.383 | 4 / 1.804 | 4 / 1.664 | 4 / 1.605 | 4 / 1.572 | 4 / 1.556 | 4 / 1.548 |
| body | 4 / 7.628 | 4 / 3.253 | 4 / 2.142 | 8 / 1.555 | 8 / 1.324 | 8 / 1.225 | 8 / 1.170 | 8 / 1.152 | 8 / 1.139 |
| bump | 0 / 0.125 | 0 / 0.125 | 0 / 0.125 | 0 / 0.125 | 0 / 0.125 | 1 / 0.0399 | 1 / 0.0258 | 4 / 0.0236 | 5 / 0.0196 |
| camel | 0 / 1.724 | 1 / 1.621 | 1 / 0.734 | 4 / 0.313 | 4 / 0.169 | 4 / 0.121 | 4 / 0.101 | 4 / 0.0954 | 4 / 0.0916 |
| camel_noise1 | 0 / 1.731 | 1 / 1.549 | 1 / 0.764 | 10 / 0.334 | 12 / 0.187 | 13 / 0.139 | 15 / 0.118 | 15 / 0.112 | 15 / 0.108 |
| camel_noise2 | 0 / 1.555 | 2 / 1.439 | 6 / 0.829 | 13 / 0.498 | 16 / 0.305 | 22 / 0.243 | 26 / 0.216 | 29 / 0.213 | 28 / 0.206 |
| camelbad | 0 / 1.727 | 1 / 1.624 | 2 / 0.735 | 4 / 0.313 | 4 / 0.169 | 6 / 0.121 | 7 / 0.101 | 6 / 0.0956 | 5 / 0.0917 |
| camelbadnonuniform | 0 / 1.723 | 1 / 1.614 | 1 / 0.713 | 4 / 0.310 | 4 / 0.164 | 5 / 0.117 | 5 / 0.0971 | 5 / 0.0918 | 5 / 0.0880 |
| cow (four files, one mesh) | 0 / 2.908 | 0 / 2.908 | 2 / 1.598 | 5 / 0.863 | 5 / 0.655 | 5 / 0.576 | 5 / 0.537 | 5 / 0.525 | 5 / 0.516 |
| face | 0 / 0.0203 | 0 / 0.0203 | 0 / 0.0203 | 0 / 0.0203 | 0 / 0.0203 | 0 / 0.0203 | 1 / 0.0137 | 1 / 0.0063 | 3 / 0.0040 |
| hand | 1 / 6.654 | 3 / 3.351 | 3 / 2.292 | 5 / 1.731 | 5 / 1.465 | 5 / 1.364 | 5 / 1.314 | 5 / 1.299 | 5 / 1.288 |
| rhino | 0 / 2.228 | 0 / 2.228 | 1 / 1.158 | 4 / 0.850 | 7 / 0.587 | 9 / 0.381 | 9 / 0.307 | 9 / 0.290 | 9 / 0.280 |
| tree | 0 / 1.033 | 0 / 1.033 | 4 / 0.929 | 4 / 0.278 | 4 / 0.159 | 4 / 0.118 | 4 / 0.0998 | 4 / 0.0949 | 4 / 0.0912 |
| closed, genus 0: 7 primitives, axes, bunny, chair | the same at every λ | | | | | | | | |

The `lambda` and `signs` runs were made before the merge guard got its final default (they
used factor 2, +0.01). The only rows this changes are the mannequin's below `λ = 0.1` and the
voxel slabs, which are quoted from a rerun.

Total cone angle of the camel: 0, 0.22, 3.48, 6.66, 8.09, 8.80, 9.23, 9.37, 9.49 for
`λ = 2 … 0.002`. (The mannequin and camelshotnoise rows are in the log: the first has a ring of
76 to 611 raw cones below `λ = 0.5` that merging turns into one or two, the second has 22
cones at the angle bound at every `λ`.)

- With positive cones the **number** of cones saturates near `λ = 0.2` on the clean meshes
  (camel 4, cow 5, body 8, airplane 4, hand 5, tree 4; rhino 7 of its 9 at 0.1). Below that
  only the angles grow, because the penalty shrinks them less.
- From 0.5 to 0.1 the distortion falls by a factor 1.4 (airplane) to 5.8 (tree) on the meshes
  that get cones. From 0.1 to 0.02 it falls by 1.06 to 1.7 with the same cones on the clean
  meshes (rhino goes from 7 to 9 cones and gains 1.9), while the two noisy camels get 25 % and
  63 % more cones.
- Nearly developable disks (bump, face) get no cone at 0.1; their first appears at 0.05 and
  0.02.
- On closed pieces of genus 0 `λ` does nothing, as the theory says.

**Default: `λ = 0.1`**, the largest value at which the cone set of the clean meshes with
boundary is saturated. This is the solver's convention (area of every piece normalized to 1);
its README has what is known about the relation to the paper's values.

**Re-optimized angles** (`evaluate reopt`, option `reoptimize_angles`): a second solve in
which cones cost nothing at the chosen vertices and are forbidden elsewhere, the step the
paper's authors evidently applied to their published cones. Distortion before → after:

| mesh | λ 0.5 | λ 0.2 | λ 0.1 | λ 0.02 |
|---|---|---|---|---|
| airplane | 2.383 → 1.546 | 1.804 → 1.546 | 1.664 → 1.546 | 1.572 → 1.546 |
| body | 2.142 → 1.748 | 1.555 → 1.135 | 1.324 → 1.135 | 1.170 → 1.135 |
| camel | 0.734 → 0.759 (one cone, driven to the angle bound) | 0.313 → 0.0909 | 0.169 → 0.0908 | 0.101 → 0.0907 |
| cow | 1.598 → 0.936 | 0.863 → 0.515 | 0.655 → 0.514 | 0.536 → 0.514 |
| hand | 2.292 → 1.961 | 1.731 → 1.285 | 1.465 → 1.285 | 1.314 → 1.285 |
| mannequin | 0.377 → 0.107 (cone at the bound) | 0.132 → 0.107 | 0.111 → 0.107 | 0.0662 → 0.0549 |
| bump | no cone | no cone | no cone | 0.0258 → 0.0336 |
| bunny, chair (closed) | 0.22405 → 0.22319, 0.02450 → 0.02450 | same | same | same |

It gives, for the cones of `λ = 0.2` or 0.1, the distortion that `λ → 0` approaches, without
the extra cones a small `λ` puts on noisy meshes, for 40 to 90 % more time. It does nothing on
closed pieces of genus 0. With a single cone it can drive the angle to the bound and make the
map worse by the scale-free measure (camel at 0.5, bump at 0.02): the solver's energy, with
`u = 0` on the boundary, contains a scale offset that the measure ignores. Hence **off by
default**.

## The spike check

**Why.** A vertex pulled far out of the surface dominates the distortion measure of this
method (`ddg-optimal-cones` README, "Needle-like spikes"): a spike of height `h` over a base of
radius `ρ` adds about `h³/ρ` to `½ Σ m u²` without cones, and the optimum then spends its cones
on the spike tips. On `camelshotnoise` the result is unusable (22 cones at the angle bound,
distortion 53, 17.9 % overlap, 14 742 of 22 968 triangles mirrored). Mirroring and overlap show
that afterwards; the check says what the cause is, before the solve, with the vertices.

**Criterion** (`src/spikes.rs`, options `UnwrapOptions::spike_check = Some(SpikeCheck { .. })`).
Per piece, evaluate the state with no cones (`ddg_optimal_cones::evaluate_cones` with an empty
list; cones the caller prescribed are kept, since the solver cannot change them): the log scale
`u_i` and `E0 = ½ Σ m_i u_i²`, masses normalized to sum 1. Sort `c_i = m_i u_i²` in decreasing
order and let `k` be the smallest number whose `c_i` add up to at least a share `s = 0.99` of
`Σ c_i`. The piece **has spikes** if `E0 ≥ 1e-9` (below that the `u_i` are rounding noise: a
flat grid has `E0 = 2·10⁻³¹`), `k ≤ 64`, and `k ≤ 1 %` of the piece's vertices (so a piece of
fewer than 100 vertices never has any). The spikes are those `k` vertices.

**Output.** `ChartReport::spikes: Option<SpikeReport { cone_free_distortion, n_carrying,
vertices, share, excluded, not_excluded_boundary }>` with the vertices in the caller's numbering, ascending (empty unless the
criterion holds); in `Report::summary()` the keys `cone_free_distortion` (largest over the
charts), `n_spike_vertices`, `spike_vertices` (text, at most 64 numbers; only if there are
spikes) and `seconds_spike_check`; and a warning per piece with spikes that names up to 16 of
the vertices, says that the map is likely to be poor and suggests repairing the mesh (moving the
vertices back onto the surface or removing them). `spike_check: None` switches it off; with `action: WarnOnly` the map does not depend on it
(tested). If the evaluation itself fails (the Laplacian cannot be
factorized) the check is skipped with a warning.

**Measured** (`cargo run --release -p ddg-uv --example spikes`; `tests/dataset.rs`,
`tests/spikes.rs`). `k` for 99 % against the number of vertices of the piece:

| mesh | vertices | `E0` | `k` (99 %) | `k / vertices` |
|---|---:|---:|---:|---:|
| **camelshotnoise** | 11 513 | 8.48·10⁴ | **22** | **0.19 %** |
| the 22 other dataset meshes | 2 349 to 152 893 | 0.0196 to 21.1 | 1 599 (airplane) to 108 942 (brain) | 67 % to 88 % |
| generated: UV sphere 32×16, icosphere, cube, cube 8×8, cylinder, cone, capsule, quad torus, triangle torus, genus-2 slab, bumpy disk | 8 to 625 | 8·10⁻⁸ to 0.22 | 8 to 467 | 63 % to 100 % |

The flat disk and the flat open cylinder (no curvature, `E0` at most 10⁻³¹) have no spikes by the
`E0` rule. The camel's 22 spike tips are the 22 vertices reported, and every one of them lies on an edge
more than ten times longer than the median edge (the test checks this). Between the smallest
fraction of any other mesh with 100 or more vertices (63 %) and the largest of the spiked one (0.19 %) there is a factor
over 300, against the limit of 1 %; the limits are not tuned further. Cost, from the report's
own timer in the release run of `tests/dataset.rs`: 0.004 to 0.061 s per mesh, 0.72 s on the
brain (152 893 vertices), 1.3 s of the 70 s for all 23 meshes (1.9 %); 1 % to 6 % of one mesh's
time (6 % on `bump` and 5 % on `face`, which have 0 cones and so no long cone solve).

**What it does not do.** It fires only when the spikes carry 99 % of the distortion: a
mesh with a few spikes on top of a rough surface whose own cone-free distortion is more than 1 %
of theirs is not flagged. It does not repair anything, and (with `WarnOnly`) does not set the distortion weight of
the spike tips to zero.
A thin spike that is a real feature (a very thin horn) is flagged too: for this method it is the
same thing.

**Exclusion** (`SpikeCheck::action`, default `SpikeAction::ExcludeFromDistortion`). The solver's
energy gives a one-vertex needle a log scale factor of about -15 000 on `camelshotnoise`, while in
the realized map the needle triangles come out at about -4 (the flattener places the tip by
interpolating between its neighbours) on 0.26 % of the area: the energy overstates their cost by
orders of magnitude and the optimum spends cones on them. With the default the detected tips
get `distortion_weights = 0` in the cone problem (the remedy of the solver's README) and the
warning says how many were excluded. `WarnOnly` keeps the old behaviour. Only vertices that are
interior to their piece are excluded; a tip on the boundary is listed in
`SpikeReport::not_excluded_boundary` (with optimized boundary values a zero weight at a free
boundary vertex makes the problem ill-posed). If the weights would all be zero the piece is left
alone (the detector's limits make that impossible: at most 1 % of the vertices). The detection
itself always uses the state with `u = 0` on the boundary: with a free boundary the cone-free
state spreads a needle's curvature over the piece (on a disk with two needles 317 of 625
vertices then carry 99 %), which is no longer a local signature.

`camelshotnoise` (release, all other settings default; `evaluate spikeaction`):

| boundary | spikes | cones | distortion | published | overlap % | mirrored | seconds |
|---|---|---:|---:|---:|---:|---:|---:|
| `Zero` | `WarnOnly` | 22 | 52.967 | 0.0953 | 17.883 | 14 742 | 0.63 |
| `Zero` | `ExcludeFromDistortion` (default) | 4 | 0.2595 | 0.0953 | 0.001 | 116 | 4.5 |
| `Optimized` | `WarnOnly` | 22 | 29.287 | 0.0953 | 19.832 | 10 908 | 0.50 |
| `Optimized` | `ExcludeFromDistortion` | 2 | 0.9612 | 0.0953 | 13.123 | 20 804 | 5.95 |

The other 22 dataset meshes and the generated meshes are unchanged (the detector does not fire
on them; `tests/dataset.rs` still asserts their warning-free results). The cone solve is slower
with the exclusion on `camelshotnoise` (4.5 s instead of 0.5 s). Excluded tips still count in the
distortion of the realized map, so the pipeline's measured distortion (0.26) stays above the
published map's 0.095.

## The dataset

`cargo run --release -p ddg-uv --example evaluate -- dataset` (`tests/dataset.rs` runs the same
23 meshes and asserts what is said below). All defaults, one thread, nothing else running;
AMD Ryzen Threadripper 3970X. "Distortion" is the scale-free README measure of our map,
"published" the same measure on the UV map stored in the OBJ file; times in seconds (cutting
seams takes at most 1 ms here and is left out; there are none).

| mesh | triangles | cones | Σ\|µ\| | distortion | published | overlap % | mirrored | charts | validate | cones | flatten | pack + verify | total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| airplane | 4656 | 4 | 9.955 | 1.66445 | 0.12701 | 0.000 | 0 | 1 | 0.002 | 0.205 | 0.008 | 0.001 | 0.215 |
| axes | 7466 | 4 | 12.566 | 0.80698 | 0.80892 | 0.000 | 0 | 1 | 0.003 | 0.742 | 0.012 | 0.001 | 0.759 |
| body | 23892 | 8 | 11.754 | 1.32411 | 0.39335 | 0.000 | 0 | 1 | 0.011 | 2.036 | 0.036 | 0.004 | 2.087 |
| brain | 305782 | 55 | 12.566 | 0.04700 | 0.03300 | 0.000 | 38 | 1 | 0.150 | 33.154 | 0.748 | 0.058 | 34.112 |
| bump | 6528 | 0 | 0.000 | 0.12501 | 0.01969 | 0.000 | 0 | 1 | 0.002 | 0.061 | 0.010 | 0.001 | 0.075 |
| bunny | 28576 | 15 | 12.566 | 0.22405 | 0.23229 | 0.000 | 0 | 1 | 0.012 | 2.916 | 0.055 | 0.005 | 2.988 |
| camel | 22704 | 4 | 8.088 | 0.16874 | 0.03740 | 0.000 | 0 | 1 | 0.009 | 1.411 | 0.042 | 0.004 | 1.466 |
| camel_noise1 | 22704 | 12 | 8.406 | 0.18690 | 0.05123 | 0.001 | 19 | 1 | 0.009 | 1.569 | 0.043 | 0.004 | 1.624 |
| camel_noise2 | 22704 | 16 | 9.715 | 0.30481 | 0.18262 | 0.028 | 108 | 1 | 0.009 | 1.576 | 0.041 | 0.004 | 1.630 |
| camelbad | 22725 | 4 | 8.085 | 0.16894 | 0.03681 | 0.000 | 0 | 1 | 0.009 | 1.231 | 0.041 | 0.004 | 1.285 |
| camelbadnonuniform | 18125 | 4 | 8.048 | 0.16391 | 0.03593 | 0.000 | 0 | 1 | 0.007 | 0.991 | 0.033 | 0.003 | 1.034 |
| camelshotnoise | 22968 | 22 | 131.319 | 52.96730 | 0.09531 | 17.883 | 14742 | 1 | 0.009 | 0.487 | 0.047 | 0.004 | 0.548 |
| chair | 11448 | 7 | 12.566 | 0.02450 | 0.01275 | 0.000 | 0 | 1 | 0.005 | 0.492 | 0.018 | 0.001 | 0.517 |
| cow00001 | 29999 | 5 | 8.753 | 0.65469 | 0.00000 | 0.000 | 0 | 1 | 0.012 | 1.463 | 0.053 | 0.005 | 1.533 |
| cow01 | 29999 | 5 | 8.753 | 0.65469 | 0.03008 | 0.000 | 0 | 1 | 0.013 | 1.536 | 0.059 | 0.005 | 1.613 |
| cow1 | 29999 | 5 | 8.753 | 0.65469 | 0.90993 | 0.000 | 0 | 1 | 0.014 | 1.547 | 0.057 | 0.005 | 1.623 |
| cow5e-1 | 29999 | 5 | 8.753 | 0.65469 | 0.42601 | 0.000 | 0 | 1 | 0.013 | 1.501 | 0.054 | 0.005 | 1.573 |
| cow75e-2 | 29999 | 5 | 8.753 | 0.65469 | 0.72038 | 0.000 | 0 | 1 | 0.013 | 1.534 | 0.053 | 0.006 | 1.605 |
| face | 27150 | 0 | 0.000 | 0.02030 | 0.00381 | 0.000 | 0 | 1 | 0.012 | 0.440 | 0.044 | 0.005 | 0.501 |
| hand | 6660 | 5 | 10.193 | 1.46453 | 0.01050 | 0.000 | 0 | 1 | 0.003 | 0.485 | 0.010 | 0.001 | 0.499 |
| mannequin | 12195 | 1 | 5.844 | 0.11059 | 0.01757 | 0.000 | 2 | 1 | 0.005 | 0.428 | 0.020 | 0.002 | 0.456 |
| rhino | 40304 | 7 | 6.888 | 0.58652 | 0.03754 | 0.000 | 0 | 1 | 0.016 | 2.736 | 0.070 | 0.006 | 2.829 |
| tree | 33592 | 4 | 6.854 | 0.15879 | 0.04316 | 0.000 | 0 | 1 | 0.014 | 2.028 | 0.065 | 0.006 | 2.113 |

- **Every mesh is unwrapped**, as one chart, with positive cones only; no overlap above
  0.03 % of a footprint except on `camelshotnoise`.
- **Closed meshes** come out close to the published maps: axes 0.807 against 0.809 (4 cones
  each), bunny 0.224 with 15 cones against 0.232 with 4, chair 0.0245 with 7 positive cones
  against 0.0128 with 14 cones of both signs, brain 0.0470 with 55 positive cones against
  0.0330 with 45 of both signs.
- **Meshes with boundary have more distortion than the published maps**, by a factor 1.5 to
  140: cow5e-1 1.5, camel_noise2 1.7, body 3.4, tree 3.7, camel 4.5, face 5.3, bump and
  mannequin 6.3, airplane 13, rhino 16, cow01 22, hand 140. Against the two cow maps with the
  fewest published cones ours is lower (0.655 against 0.910 and 0.720). The published maps
  have more cones (hand 29, rhino 36, airplane 13 against our 5, 7, 4), negative ones among
  them, and half of them were computed with optimized boundary values (paper Sec. 8.3), which
  the solver does not implement. What smaller `λ`, negative cones and re-optimized angles each buy is in the
  tables above; none of them closes this gap.
- **Mirrored triangles** occur on four more meshes: brain 38, camel_noise1 19, camel_noise2
  108, mannequin 2. The flattener's README has 38, 19, 113 and 2 on these meshes with the
  published cones, and traces them to the triangulation.
- **`camelshotnoise`** (22 one-vertex spikes): with `SpikeAction::WarnOnly` the result is not
  usable (22 cones, all at the angle bound, distortion 53, 17.9 % overlap, 14 742 of 22 968
  triangles mirrored); the table row above is that run. The default excludes the spike tips from
  the distortion measure of the cone problem ([the spike check](#the-spike-check)): 4 cones
  (Σ|µ| 7.669), distortion 0.2595 against the published 0.0953, 0.001 % overlap, 116 mirrored
  triangles, 4.7 s; the report keeps one warning (the spike warning).
- **Time**: the cone stage is 96.6 % of the total over the 23 meshes (60.6 of 62.7 s); the
  brain takes 33.2 of its 34.1 s there.

## The boundary

`UnwrapOptions::boundary` chooses the scale factor `u` on the boundary of every piece (the
mesh's own boundary and the seams; closed pieces are untouched).

* `BoundaryValues::Zero` (default): `u = 0`, the cone solver's `Fixed` treatment (paper Fig. 10 b).
  Boundary lengths are preserved, so the two sides of a seam, which belong to different pieces,
  have the same length.
* `BoundaryValues::Optimized`: `u` on the boundary is a variable of the cone problem
  (`BoundaryTreatment::Free`, paper Sec. 8.3, Fig. 10 c); the geodesic curvature of the flattened
  boundary is free. The solver's `log_scale` (nonzero on the boundary) is handed to
  `ddg_bff::flatten_with_cones` as `scale_factors`; the flattener checks that it is flat away from
  the cut, takes the boundary lengths `e^{(u_i+u_j)/2} l_ij` and exterior angles from it, and
  closes the boundary curve. If the flatness check fails the flattener solves for its own scale
  factors again (with `u = 0` on the boundary); the warning says that the optimized values of the
  piece were dropped and `ChartReport::boundary_optimized` is `false` for it.

What changes in the pipeline with `Optimized` (the rest was checked to be indifferent to `u` on
the boundary): the solver's `measure` holds the boundary curvature at boundary vertices, so cones
are taken from interior vertices only (also when a cluster merge is rejected); the angle
re-optimization solves with the same treatment; the solver's distortion includes the boundary
vertices (cluster merging and its guard compare two such numbers, which is consistent); the
closed-surface logic, including the automatic slit, runs with `Fixed` (the slit's two sides keep
their lengths), so a closed piece gives exactly the same map as before (tested); the spike check
measures the state with `u = 0` (see its section).

**The report.**

* `Report::seam_lengths: SeamLengths`: for every cut edge that has a polygon on each side, the
  relative difference `|a - b| / max(a, b)` of the lengths of its two images in the layout (taken
  before packing, which scales all charts alike). `max_seam_mismatch`, `mean_seam_mismatch` and
  `worst_seam_edge` are over the seam edges (given ones and extension edges),
  `max_other_mismatch` over the other cut edges (slits, cuts through cones). In `summary()`:
  `max_seam_length_mismatch`, `mean_seam_length_mismatch`, `max_other_cut_length_mismatch`. A
  texture painted across a seam is stretched by this much differently on its two sides. With
  `Zero` it is what closing the boundary curve changed (about 1e-13 on benign meshes, 5 % to 28 %
  where a cone sits at the angle bound).
* `FlattenReport::boundary_length_error`: the largest `|l' - l*| / l*` over the boundary edges of
  a piece, `l'` the length in the layout, `l* = e^{(u_i+u_j)/2} l` the length the scale factors
  ask for. In `summary()`: `max_boundary_length_error`, and `boundary_values_optimized`.
* `ChartReport::boundary_optimized`.

**Does the flattening realize the solver's values?** Measured independently of the flattener's
own bookkeeping (`evaluate verify`, `tests/boundary.rs`: the solver's `log_scale` and the 3D
lengths against the corner UVs of the layout):

| mesh | boundary edges | solver `u` on the boundary (min to max) | largest relative error of a boundary length |
|---|---:|---|---:|
| tree | 240 | 1.150 to 41.3 | 1.9e-5 (but every triangle is mirrored) |
| bump | 128 | 0.420 to 0.634 | 7.4e-5 |
| body | 48 | -1.145 to -0.154 | 9.6e-4 |
| face | 340 | 0 to 0.511 | 4.2e-3 |
| cow (all five) | 97 | 0.199 to 2.673 | 8.2e-3 |
| airplane | 40 | -0.968 to -0.293 | 1.0e-2 |
| hand | 82 | -0.737 to 2.845 | 1.1e-2 |
| rhino | 104 | 0.268 to 3.607 | 9.5e-3 |
| mannequin | 63 | 0 to 5.536 | 0.135 |
| camel_noise2 | 56 | -2.697 to 8.152 | 0.166 |
| camel, camelbad, camelbadnonuniform, camel_noise1 | 56 to 57 | -0.13 to 6.2 | 0.48 to 0.55 |
| camelshotnoise (before the spike exclusion) | 56 | -122 to 1.9 | 1.9 |
| generated: vase with seams, sphere caps, torus | 64 to 120 | up to 1.25 | at most 1.7e-4 |

Where the solver's boundary values are large (a boundary vertex with `u` of 3 to 8 is common: the
objective weights a vertex by its mass, and a boundary vertex has little), the flattening does
not follow them. The likely cause, read from `ddg-bff` and not tested by changing it: the exterior
angles come from the *linearized* relation between `u` and the boundary curvature, so the curve
built from the target lengths does not close; the closing step changes lengths where it can
(reported by the flattener as `closure_length_change`, 0 on all of these meshes) and spreads what
is left of the gap linearly along the curve, which changes edge lengths without any report. The
figure above is the only place where that shows. `ddg-bff` would have to report the length change
after that spreading, or close the curve with exact angles; that is not done here.

**Evidence for the default.** All 19 dataset meshes with boundary, default `lambda = 0.1`,
positive cones, default spike handling, release build (`evaluate boundary`). "u=0" is the
default, "optimized" the option, "reopt" is `reoptimize_angles`. Distortion is the dataset
README's measure on our UVs, "published" the same measure on the UV map in the OBJ file (whose
`lambda` is that of the file name, e.g. `cow01` 0.1); overlap and mirrored triangles are of the
final map; seconds are for the whole `unwrap`.

| mesh | config | cones | distortion | published | overlap % | mirrored | seconds | boundary length error |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| airplane | u=0 | 4 | 1.66445 | 0.12701 | 0.000 | 0 | 0.25 | 9.8e-3 |
| airplane | u=0 + reopt | 4 | 1.54561 | 0.12701 | 0.000 | 0 | 0.36 | 3.4e-2 |
| airplane | optimized | 5 | 1.42606 | 0.12701 | 0.000 | 0 | 0.36 | 1.0e-2 |
| airplane | optimized + reopt | 5 | 1.08638 | 0.12701 | 0.000 | 0 | 0.56 | 3.0e-1 |
| body | u=0 | 8 | 1.32411 | 0.39335 | 0.000 | 0 | 2.20 | 5.6e-4 |
| body | u=0 + reopt | 8 | 1.13532 | 0.39335 | 0.000 | 0 | 3.62 | 6.7e-4 |
| body | optimized | 8 | 0.96929 | 0.39335 | 0.000 | 0 | 2.81 | 9.6e-4 |
| body | optimized + reopt | 8 | 0.66832 | 0.39335 | 0.000 | 0 | 4.60 | 5.0e-2 |
| bump | u=0 | 0 | 0.12501 | 0.01969 | 0.000 | 0 | 0.09 | 3.4e-5 |
| bump | u=0 + reopt | 0 | 0.12501 | 0.01969 | 0.000 | 0 | 0.09 | 3.4e-5 |
| bump | optimized | 0 | 0.12496 | 0.01969 | 0.000 | 0 | 0.15 | 7.4e-5 |
| bump | optimized + reopt | 0 | 0.12496 | 0.01969 | 0.000 | 0 | 0.14 | 7.4e-5 |
| camel | u=0 | 4 | 0.16874 | 0.03740 | 0.000 | 0 | 1.56 | 1.1e-3 |
| camel | u=0 + reopt | 4 | 0.09082 | 0.03740 | 0.000 | 0 | 2.75 | 1.3e-3 |
| camel | optimized | 2 | 1.04995 | 0.03740 | 21.438 | 21388 | 1.82 | 4.8e-1 |
| camel | optimized + reopt | 2 | 0.59655 | 0.03740 | 1.465 | 10 | 2.64 | 2.2e-1 |
| camel_noise1 | u=0 | 12 | 0.18690 | 0.05123 | 0.001 | 19 | 1.69 | 6.6e-3 |
| camel_noise1 | u=0 + reopt | 11 | 0.10786 | 0.05123 | 0.001 | 20 | 2.97 | 1.0e-2 |
| camel_noise1 | optimized | 2 | 0.93399 | 0.05123 | 25.637 | 20379 | 2.78 | 5.5e-1 |
| camel_noise1 | optimized + reopt | 2 | 0.65170 | 0.05123 | 2.376 | 40 | 3.33 | 3.1e-1 |
| camel_noise2 | u=0 | 16 | 0.30481 | 0.18262 | 0.028 | 108 | 1.71 | 6.2e-2 |
| camel_noise2 | u=0 + reopt | 14 | 0.20711 | 0.18262 | 0.030 | 108 | 3.16 | 8.9e-2 |
| camel_noise2 | optimized | 1 | 1.18063 | 0.18262 | 60.189 | 6466 | 1.39 | 1.7e-1 |
| camel_noise2 | optimized + reopt | 1 | 0.83843 | 0.18262 | 1.362 | 147 | 2.02 | 2.5e-1 |
| camelbad | u=0 | 4 | 0.16894 | 0.03681 | 0.000 | 0 | 1.29 | 1.6e-3 |
| camelbad | u=0 + reopt | 4 | 0.09101 | 0.03681 | 0.000 | 1 | 2.20 | 5.1e-4 |
| camelbad | optimized | 2 | 0.92531 | 0.03681 | 10.717 | 20401 | 1.57 | 5.2e-1 |
| camelbad | optimized + reopt | 2 | 0.60986 | 0.03681 | 3.093 | 39 | 2.25 | 2.4e-1 |
| camelbadnonuniform | u=0 | 4 | 0.16391 | 0.03593 | 0.000 | 0 | 0.96 | 2.5e-3 |
| camelbadnonuniform | u=0 + reopt | 4 | 0.08723 | 0.03593 | 0.000 | 1 | 1.54 | 6.1e-3 |
| camelbadnonuniform | optimized | 2 | 0.86746 | 0.03593 | 14.575 | 15783 | 1.04 | 5.0e-1 |
| camelbadnonuniform | optimized + reopt | 2 | 0.57271 | 0.03593 | 0.711 | 15 | 1.49 | 1.8e-1 |
| camelshotnoise | u=0 | 4 | 0.25945 | 0.09531 | 0.001 | 116 | 4.39 | 5.6e-2 |
| camelshotnoise | u=0 + reopt | 4 | 0.16974 | 0.09531 | 0.000 | 73 | 7.68 | 5.8e-2 |
| camelshotnoise | optimized | 2 | 0.96115 | 0.09531 | 13.123 | 20804 | 6.18 | 4.8e-1 |
| camelshotnoise | optimized + reopt | 2 | 0.64656 | 0.09531 | 2.237 | 85 | 8.11 | 2.7e-1 |
| cow00001 | u=0 | 5 | 0.65469 | 0.00000 | 0.000 | 0 | 1.58 | 3.6e-4 |
| cow00001 | u=0 + reopt | 5 | 0.51448 | 0.00000 | 0.000 | 0 | 2.62 | 6.7e-4 |
| cow00001 | optimized | 2 | 1.26973 | 0.00000 | 0.000 | 0 | 1.72 | 8.2e-3 |
| cow00001 | optimized + reopt | 2 | 0.89933 | 0.00000 | 0.000 | 0 | 2.35 | 2.4e-2 |
| cow01 | u=0 | 5 | 0.65469 | 0.03008 | 0.000 | 0 | 1.56 | 3.6e-4 |
| cow01 | u=0 + reopt | 5 | 0.51448 | 0.03008 | 0.000 | 0 | 2.48 | 6.7e-4 |
| cow01 | optimized | 2 | 1.26973 | 0.03008 | 0.000 | 0 | 1.60 | 8.2e-3 |
| cow01 | optimized + reopt | 2 | 0.89933 | 0.03008 | 0.000 | 0 | 2.25 | 2.4e-2 |
| cow1 | u=0 | 5 | 0.65469 | 0.90993 | 0.000 | 0 | 1.56 | 3.6e-4 |
| cow1 | u=0 + reopt | 5 | 0.51448 | 0.90993 | 0.000 | 0 | 2.47 | 6.7e-4 |
| cow1 | optimized | 2 | 1.26973 | 0.90993 | 0.000 | 0 | 1.63 | 8.2e-3 |
| cow1 | optimized + reopt | 2 | 0.89933 | 0.90993 | 0.000 | 0 | 2.27 | 2.4e-2 |
| cow5e-1 | u=0 | 5 | 0.65469 | 0.42601 | 0.000 | 0 | 1.62 | 3.6e-4 |
| cow5e-1 | u=0 + reopt | 5 | 0.51448 | 0.42601 | 0.000 | 0 | 2.49 | 6.7e-4 |
| cow5e-1 | optimized | 2 | 1.26973 | 0.42601 | 0.000 | 0 | 1.69 | 8.2e-3 |
| cow5e-1 | optimized + reopt | 2 | 0.89933 | 0.42601 | 0.000 | 0 | 2.26 | 2.4e-2 |
| cow75e-2 | u=0 | 5 | 0.65469 | 0.72038 | 0.000 | 0 | 1.59 | 3.6e-4 |
| cow75e-2 | u=0 + reopt | 5 | 0.51448 | 0.72038 | 0.000 | 0 | 2.50 | 6.7e-4 |
| cow75e-2 | optimized | 2 | 1.26973 | 0.72038 | 0.000 | 0 | 1.61 | 8.2e-3 |
| cow75e-2 | optimized + reopt | 2 | 0.89933 | 0.72038 | 0.000 | 0 | 2.27 | 2.4e-2 |
| face | u=0 | 0 | 0.02030 | 0.00381 | 0.000 | 0 | 0.49 | 4.2e-3 |
| face | u=0 + reopt | 0 | 0.02030 | 0.00381 | 0.000 | 0 | 0.49 | 4.2e-3 |
| face | optimized | 0 | 0.01886 | 0.00381 | 0.000 | 0 | 0.95 | 4.2e-3 |
| face | optimized + reopt | 0 | 0.01886 | 0.00381 | 0.000 | 0 | 0.94 | 4.2e-3 |
| hand | u=0 | 5 | 1.46453 | 0.01050 | 0.000 | 0 | 0.50 | 1.1e-2 |
| hand | u=0 + reopt | 5 | 1.28487 | 0.01050 | 0.000 | 0 | 0.76 | 1.1e-2 |
| hand | optimized | 5 | 1.60215 | 0.01050 | 0.011 | 35 | 0.80 | 1.1e-2 |
| hand | optimized + reopt | 5 | 0.98333 | 0.01050 | 0.004 | 6 | 1.14 | 2.8e-1 |
| mannequin | u=0 | 1 | 0.11059 | 0.01757 | 0.000 | 2 | 0.46 | 1.8e-1 |
| mannequin | u=0 + reopt | 1 | 0.10722 | 0.01757 | 0.000 | 2 | 0.64 | 2.0e-1 |
| mannequin | optimized | 1 | 0.42475 | 0.01757 | 2.261 | 19 | 0.72 | 1.3e-1 |
| mannequin | optimized + reopt | 1 | 0.10243 | 0.01757 | 0.008 | 4 | 1.03 | 1.7e-1 |
| rhino | u=0 | 7 | 0.58652 | 0.03754 | 0.000 | 0 | 3.07 | 3.4e-3 |
| rhino | u=0 + reopt | 7 | 0.27787 | 0.03754 | 0.000 | 0 | 4.34 | 3.7e-3 |
| rhino | optimized | 3 | 0.87693 | 0.03754 | 0.000 | 0 | 3.70 | 9.5e-3 |
| rhino | optimized + reopt | 3 | 0.65078 | 0.03754 | 0.000 | 0 | 5.39 | 1.5e-1 |
| tree | u=0 | 4 | 0.15879 | 0.04316 | 0.000 | 0 | 2.17 | 9.8e-13 |
| tree | u=0 + reopt | 4 | 0.09308 | 0.04316 | 0.000 | 0 | 3.11 | 1.1e-12 |
| tree | optimized | 0 | 1.03477 | 0.04316 | 0.000 | 33296 | 1.80 | 1.9e-5 |
| tree | optimized + reopt | 0 | 1.03477 | 0.04316 | 0.000 | 33296 | 1.90 | 1.9e-5 |

Reading the table:

* Without angle re-optimization `Optimized` has the lower distortion on `airplane` (1.426 against
  1.664), `body` (0.969 against 1.324), and marginally on `bump` and `face`; with re-optimization
  also on `hand` (0.983 against 1.285) and `mannequin`. On the five `cow` files (2 cones against 5)
  and `rhino` (3 against 7) the distortion is higher: it uses fewer cones, which is the trade the
  solver's objective makes (a lower objective, not a lower distortion).
* It breaks the map on the meshes where the boundary values are large: `tree` (0 cones, 33 296 of
  33 592 triangles mirrored), the camels (10 % to 60 % overlap and 6 000 to 21 000 mirrored
  triangles without re-optimization; 0.7 % to 3.1 % overlap with it) and `mannequin` (2.3 %
  overlap, 19 mirrored). The pipeline's warnings fire on them.
* The gap to the published maps is not closed: on every mesh with boundary the distortion stays
  several times the published one (`cow01`: 0.030 published, 0.65 ours with `u=0`, 1.27 with
  `Optimized`; `cow1` and `cow75e-2`, published with 2 and 3 cones: 0.91 and 0.72 against 1.27
  with 2 cones).
* Time: `Optimized` costs 1.0 to 1.9 times as much as `u=0` on these meshes.

Generated meshes with seams (`evaluate boundary`): a UV sphere cut by a latitude ring into two
unequal caps, an open vase (profile `1 + 0.6 sin 3t + 0.3 t`) without seam, with a straight seam
(one piece cut open) and cut by a ring into two tubes, and a torus cut along a meridian and a
longitude.

| mesh | config | cones | distortion | overlap % | mirrored | boundary length error | seam edges | max seam mismatch |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| uv_sphere 32x16, ring 3 seam (two unequal caps) | u=0 | 1 | 0.06831 | 0.000 | 0 | 4.9e-2 | 32 | 4.95e-2 |
| uv_sphere 32x16, ring 3 seam (two unequal caps) | u=0 + reopt | 1 | 0.07024 | 0.000 | 0 | 2.8e-1 | 32 | 2.82e-1 |
| uv_sphere 32x16, ring 3 seam (two unequal caps) | optimized | 0 | 0.21599 | 0.000 | 0 | 3.3e-14 | 32 | 6.89e-1 |
| uv_sphere 32x16, ring 3 seam (two unequal caps) | optimized + reopt | 0 | 0.21599 | 0.000 | 0 | 3.3e-14 | 32 | 6.89e-1 |
| uv_sphere 32x16, ring 5 seam (two unequal caps) | u=0 | 1 | 0.07789 | 0.000 | 0 | 6.1e-4 | 32 | 6.13e-4 |
| uv_sphere 32x16, ring 5 seam (two unequal caps) | u=0 + reopt | 1 | 0.02682 | 0.000 | 0 | 1.0e-1 | 32 | 1.05e-1 |
| uv_sphere 32x16, ring 5 seam (two unequal caps) | optimized | 0 | 0.07612 | 0.000 | 0 | 2.5e-14 | 32 | 3.96e-1 |
| uv_sphere 32x16, ring 5 seam (two unequal caps) | optimized + reopt | 0 | 0.07612 | 0.000 | 0 | 2.5e-14 | 32 | 3.96e-1 |
| open vase 32x16, no seam | u=0 | 0 | 0.01165 | 0.000 | 0 | 7.5e-5 | 0 | 0.00e0 |
| open vase 32x16, no seam | u=0 + reopt | 0 | 0.01165 | 0.000 | 0 | 7.5e-5 | 0 | 0.00e0 |
| open vase 32x16, no seam | optimized | 0 | 0.01163 | 0.000 | 0 | 5.7e-5 | 0 | 0.00e0 |
| open vase 32x16, no seam | optimized + reopt | 0 | 0.01163 | 0.000 | 0 | 5.7e-5 | 0 | 0.00e0 |
| open vase 32x16, straight seam | u=0 | 0 | 0.01234 | 0.000 | 0 | 1.2e-4 | 15 | 1.21e-13 |
| open vase 32x16, straight seam | u=0 + reopt | 0 | 0.01234 | 0.000 | 0 | 1.2e-4 | 15 | 1.21e-13 |
| open vase 32x16, straight seam | optimized | 0 | 0.00964 | 0.000 | 0 | 1.7e-4 | 15 | 2.20e-3 |
| open vase 32x16, straight seam | optimized + reopt | 0 | 0.00964 | 0.000 | 0 | 1.7e-4 | 15 | 2.20e-3 |
| open vase 32x16, ring 4 seam (two tubes) | u=0 | 0 | 0.00528 | 0.000 | 0 | 2.7e-6 | 32 | 3.73e-6 |
| open vase 32x16, ring 4 seam (two tubes) | u=0 + reopt | 0 | 0.00528 | 0.000 | 0 | 2.7e-6 | 32 | 3.73e-6 |
| open vase 32x16, ring 4 seam (two tubes) | optimized | 0 | 0.00241 | 0.000 | 0 | 1.4e-5 | 32 | 1.23e-1 |
| open vase 32x16, ring 4 seam (two tubes) | optimized + reopt | 0 | 0.00241 | 0.000 | 0 | 1.4e-5 | 32 | 1.23e-1 |
| torus 48x12, meridian + longitude seams | u=0 | 0 | 0.02568 | 0.000 | 0 | 1.7e-13 | 60 | 3.39e-13 |
| torus 48x12, meridian + longitude seams | u=0 + reopt | 0 | 0.02568 | 0.000 | 0 | 1.7e-13 | 60 | 3.39e-13 |
| torus 48x12, meridian + longitude seams | optimized | 0 | 0.02288 | 0.000 | 0 | 2.3e-13 | 60 | 8.79e-3 |
| torus 48x12, meridian + longitude seams | optimized + reopt | 0 | 0.02288 | 0.000 | 0 | 2.3e-13 | 60 | 8.79e-3 |

The seam mismatch with `Optimized` is 0.2 % (straight seam on the vase), 0.9 % (torus), 12 %
(vase cut into two tubes), 40 % and 69 % (sphere caps). The distortion falls on the vase
(0.0053 to 0.0024 for the two tubes) and the torus (0.0257 to 0.0229) and, on the sphere caps,
the optimized map has no cone where `u=0` has one, with higher distortion at ring 3 (0.216 against
0.068) and about equal at ring 5. With `u=0` the mismatch is at most 6e-4 here, except where a cone
sits at the angle bound (sphere with `reopt`: 10 % to 28 %, an effect of the closing step that
exists without the option and is visible now that it is measured).

**Recommendation.** Keep `Zero` as the default of the pipeline and of the add-on, which is how
they ship. `Optimized` lowers the distortion only on some meshes (`airplane`, `body`, and `hand`
and `mannequin` with re-optimization, and simple open surfaces), produces overlap or mirrored
triangles on eight of the nineteen dataset meshes, does not match the published maps either, and
costs equal seam lengths. What keeping `Zero` costs: on meshes with boundary the distortion stays
at 1.5 to 140 times the published one (the table in [the dataset](#the-dataset)), and the lower
objective the solver finds for the free boundary stays unused. It would become a reasonable
default for the add-on if `ddg-bff` closed the boundary curve with exact angles (so that large
`u` is realized) or if the solver bounded `u` on the boundary; both are changes outside this
crate. Until then the option is for pieces whose boundary values are small and for users who
accept unequal seams; the report says which of the two happened (`max_boundary_length_error`,
`max_seam_length_mismatch`, the warnings).

## The cut flow on polygon meshes

`ddg_uv::cutting` makes the cut flow of

> Nicholas Sharp, Keenan Crane. *Variational Surface Cutting.*
> ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201356>

(crate `ddg-variational-cuts`, implemented clean-room from the text of the paper) usable from
polygon meshes as an editor provides them (polygons, a triangulation, seams), end to end to UVs.
It is meant to sit under a thin binding. The flow itself and its options are described in that
crate's README and `SPEC.md`.

```text
Session::new(&MeshInput, n_regions, &Start) -> Session   the flow on the triangulation cut open at the seams
  .step / .run / .roughening_rate / .evaluate            the flow's own problem, options, stop rule and reports
  .state() / .set_state(&regions)                        region functions per polygon corner
  .restart(&Start)  .corner_regions()
  .preview(delta) -> Vec<PreviewPolyline>                the cut as polylines in space
  .describe_cut(delta) -> CutOnPolygons                  edge points and the paths through each polygon
insert(&PolygonMesh, &CutOnPolygons) -> Inserted         the crossed polygons split along the paths
Session::unwrap / unwrap_with_cut(&MeshInput, &CutOnPolygons, &CutUnwrapOptions) -> CutUnwrap
                                                         insert, then `unwrap` with seams = old seams + cut
```

* **Session.** Built from the caller's `MeshInput` (positions, corners, triangulation, seams).
  The seams are handled by `prepare` (a seam that is not an edge is ignored with a warning, a
  seam of one edge is extended, `Session::warnings`); the triangle mesh the flow runs on is the
  triangulation cut along them, so the seams are boundary for the flow. The session fixes two
  policies: `Orientation::Reject` and `Degenerate::Reject`. The flow works on true geometry, and
  a flipped or mollified mesh would no longer be the caller's. The pieces that the seams
  separate are one mesh of several components. `flow_vertex_to_vertex`, `corner_to_flow_vertex`,
  `flow_face_to_triangle` give the maps.
* **State.** `state()` returns `n_regions` arrays with one `f64` per polygon corner in the
  caller's length unit; a corner belongs to the region whose function is smallest (ties: lowest).
  For two regions the functions are `phi` and `-phi`, and `set_state` uses `phi = (f0 - f1) / 2`.
  Corners of one vertex that no seam separates must carry the same value (exact equality):
  `CuttingError::CornersDisagree` names the vertex, the region and the two corners. Across a
  seam the corners may differ; each side is its own vertex of the flow.
* **Starts** (`Start`): `GeodesicCircle { vertex, radius }` (distance on the surface cut open
  along the seams, from every copy of the vertex; it does not reach across a seam),
  `PolygonLabels` (one label per polygon, below `n_regions`; the flow's `init::from_face_labels`
  on the triangles), `Corners` (the state as above).
* **Controls.** `problem.penalty` (`eta`) and `ScaleWeight::PerVertex` are per vertex of the
  polygon mesh, `CurvatureSource::PerFace` per triangle of the supplied triangulation,
  `options.symmetry` one class per polygon-mesh vertex. A vertex that a seam splits into several
  flow vertices is left out of the symmetry (each copy is its own class). Reports are the flow's
  (`StepReport`, `RunReport`).
* **The cut on the polygons** (`CutOnPolygons`): `edge_points` (the polygon edge as a vertex pair
  `[a, b]`, `a < b`, the parameter `t` from `a`, the position, and `side: Some(polygon)` on an
  edge that is a seam, where the two sides carry separate points) and `paths` (per polygon, an
  ordered list from one edge point to another, the crossed diagonals of the triangulation
  in between as `InteriorPoint { corners, t, position }`, and the regions on its left and
  right). `PathEnd::Junction` and `CutOnPolygons::junctions` are there for the day the flow has
  three or more regions; `insert` refuses a junction with `CuttingError::Unsupported`.
* **Insertion.** `insert` needs only the polygon mesh and the description, which is what the
  Blender operator will do with bmesh. Every edge point becomes a vertex on its edge and **both**
  polygons of the edge get it; on a seam the two sides' points are separate vertices, and both
  polygons get both (so a polygon can gain collinear corners; two points at exactly the same
  parameter are one vertex). Polygons crossed by paths are cut along them; the new polygons are
  unions of pieces of the old triangles, and their triangulation is made of those pieces (a
  quadrilateral piece takes the diagonal with the better worse triangle; a piece with extra
  collinear corners is fanned). Old vertices keep their numbers; new vertices follow. The result
  (`Inserted`) holds the new mesh with its triangulation, `new_seams` (exactly the pieces between
  consecutive points of the paths), `polygon_parent`, `vertex_source`, `path_vertices`, and
  `seams_after(old_seams)` (an old seam edge that received points becomes the pieces of its chain).
* **End to end.** `Session::unwrap` / `unwrap_with_cut` insert and call `unwrap` with
  `seams_after(old seams)`. **Cones.** The existing options have no switch for cones as such:
  they are steered per vertex (`forbidden_vertices` gives infinite cost, `pinned_vertices` and
  `prescribed_cones` ask for cones). `ConePlacement::Optimal` (default) leaves `unwrap` as it is
  (the steering of the old mesh stays valid, the vertices keep their numbers);
  `ConePlacement::Off` forbids every vertex of the new mesh and refuses pinned or prescribed
  ones. Nothing in `unwrap` changed.
* **The triangulation of split polygons.** `Inserted::mesh.triangles` is the triangulation the
  flattening should see: the pieces are sub-triangles of the supplied triangulation, ear-clipped
  by shape. An editor that inserts the cut itself should hand `unwrap` these triangles for the
  polygons that were split (the Blender add-on stores them at Apply and the unwrap reads them
  back); letting the editor triangulate the split polygons again flipped 2 to 5 thin triangles on
  the spot atlases (8 clusters) that the library's own triangulation does not (measurements:
  `blender/README.md`). No change in this crate was needed.
* **Several regions.** `n_regions > 2` runs the core crate's `RegionFlow` (SPEC §8); `Session::new_with(.., FlowKind::Regions)` runs it for two regions too (the two-region session
  stays on `Flow`, bit for bit the core's flow, tested). State: `n_regions` functions per corner;
  `corner_regions` is the argmin. Starts: `PolygonLabels` (N-valued: the core's
  `regions_from_face_labels`), `NormalClusters { max_iterations }` (the core's `normal_clusters`
  into `n_regions` clusters, then `regions_from_vertex_labels`), `Corners`; a geodesic circle
  starts two regions only. **Options:** `Session::recommended_options(step_size)` is
  `FlowOptions::new` for two regions and `FlowOptions::for_regions` for more (the stages of
  the paper with the flux recovery; the core's SPEC D20 says why: on networks the default
  stages come to rest with more patches and a higher energy from a speckled start);
  `Session::new_problem(energy, alpha)` is `CutProblem::for_regions`. Options built with
  `FlowOptions::new` are accepted for any N and are not overridden. `Session::flow()` /
  `flow_mut()` are the two-region flow; in a session that runs `RegionFlow`
  (`uses_region_flow()`) they are an unused flow on the same mesh, and `region_flow()` is the
  flow.
  **The cut is a network.** `describe_cut` gives a `JunctionPoint` for every triple point
  (`position`, `polygon`, `triangle` of the supplied triangulation, its `corners`, and the
  `barycentric` weights of those corners: the junction is `sum(w_k * x[corners[k]])`), and three
  paths end there, one through each side of the triangle. A path runs between edge points and
  junctions; it is directed with the lower region id on its left, `regions` is `[lower, higher]`.
  `preview` gives the curves between junctions (`Cut::curves`) with their two regions.
  **`insert` with junctions** builds a vertex `J` at the weights inside the junction's triangle
  and replaces the triangle by three quadrilaterals `(c_k, q_k, J, q_{k-1})`, one per corner,
  where `q_k` is the vertex of the path on side `k` of the triangle; the whole construction
  (edge vertices, pieces, parts, triangulation) is written out in the rustdoc of
  `ddg_uv::cutting::insert`. Vertex numbers: old, edge points, junctions
  (`VertexSource::Junction`), interior points. Every region component is a chart in
  `unwrap_with_cut` (the network is in the seams).

### Tests (`tests/cutting.rs`, 36 tests in all; the intrinsic ones are below)

Insertion invariants on `uv_sphere`, `grid_cube`, `capsule`, `torus`, `cylinder` (cuts from the
flow) and on 5156 paths from random states of three meshes with and without seams
(`random_states_insert_cleanly`): total triangulated area kept to `1e-12`, no triangle without
area, the new mesh passes `prepare` (manifold, no flipped polygon), old vertices keep position
and number, every inserted point lies on its old edge or diagonal to `1e-12`, every cut segment
is an edge of the new mesh and `new_seams` is exactly those, old seams are kept (as they were or
as their chains), each new polygon lies in one region of the state, a closed cut on the sphere
gives two charts. State round trip through corners (bit for bit), the vertex-naming error across
no seam and its absence across a seam, every error variant, descriptions that do not fit, a cut
that ends on a seam (edge points on the two sides of a seam), a ring of seams that splits the
surface into two flow components, a one-edge seam, end to end with `ConePlacement::Off` and
`Optimal` on the generated meshes, the quad model `spot_quadrangulated.obj` (`DDG_MODELS`,
default `data/models/spot/spot`) and the dataset (`DDG_DATASET`; `bump` in a debug build, `camel`,
`bunny`, `hand` in a release build); a missing file is skipped with a notice. The UV checks are
the report's: no flipped triangle (asserted), every seam of the cut is a cut edge (asserted),
the seam length mismatch (asserted below 0.5, reported), the area distortion (reported).
Measured in a release build (30 flow steps on spot, 5 on the others, Hencky energy, default
options; `ConePlacement::Off` / `Optimal`): no flipped triangle in any run.

| mesh | polygons -> after insertion | max seam mismatch | scale-free distortion | cones |
|---|---|---|---|---|
| uv_sphere 24x12 | 288 -> 316 | 2.7e-3 / 5.7e-3 | 0.209 / 0.071 | 0 / 1 |
| grid_cube 6 | 216 -> 239 | 9.2e-3 / 0.13 | 0.488 / 0.098 | 0 / 1 |
| capsule 20x6 | 260 -> 290 | 0.45 / 1.4e-2 | 1.109 / 0.012 | 0 / 65 |
| torus 24x10 | 240 -> 256 | 1.2e-3 / 1.2e-3 | 0.044 / 0.044 | 0 / 0 |
| spot (quads), 30 steps | 2928 | 7.6e-3 / 4.3e-3 | 3.285 / 0.324 | 0 / 7 |
| bump | 6528 -> 6661 | 9.4e-5 / 9.4e-5 | 0.093 / 0.093 | 0 / 0 |
| camel | 22704 -> 22850 | 2.0e-3 / 4.7e-3 | 1.453 / 0.244 | 0 / 3 |
| bunny | 28576 -> 28885 | 6.3e-4 / 2.7e-3 | 4.424 / 0.286 | 0 / 2 |
| hand | 6660 -> 6759 | 2.5e-2 / 5.3e-2 | 4.942 / 0.663 | 0 / 5 |

What this says: the flow was run for a handful of steps from a geodesic circle, which is not a
cut that flattens well; with cones off the charts are flattened as they are and the distortion
is high on everything that is not nearly developable. The numbers show that the chain works and
what it costs, not how good the cuts are. The seam mismatch of the cones-off runs is what
closing the boundary of a chart changes (see [the boundary](#the-boundary)); the large value on
the capsule (0.45) is that.

**Networks (N regions).** `networks_insert_cleanly`: random states of 3, 4 and 6 regions on
sphere, cube, torus and a cylinder (its caps are 12-gons, fanned), with and without seams: 6749
junctions in 23491 paths, 1716 polygons with two or more junctions, 68 junctions in a polygon
that also has an ordinary path, 943 junctions in a polygon next to a seam; all insertion
invariants of the two-region tests hold, and a junction is a vertex of exactly three new
polygons and three seam edges. `one_junction_by_hand` (a fan of six triangles, three regions),
`n_equals_two_through_the_region_flow` (two regions on `RegionFlow` against `Flow`: energies to
`1e-9`, same description), `two_regions_are_the_core_flow_bit_for_bit` (the two-region session
against `ddg_variational_cuts::Flow` on the same triangle mesh: the bits of the energy, the
step and the state are equal), `n_region_state_and_starts`. Atlas end to end
(`atlas_end_to_end_*`: `NormalClusters` with 6 regions (8 on spot), `recommended_options`,
rescaled Dirichlet energy with `r_eps = 1`, `alpha_L` from the ripple rule with wavelength 0.5,
60 steps (100 on spot), release build; charts are the pieces the seams leave, equal to the
number of charts in the report):

| mesh | `alpha_L` | steps (stop) | energy | junctions, paths | charts | cones (off / optimal) | flipped | seam mismatch | distortion (scale-free) |
|---|---|---|---|---|---|---|---|---|---|
| uv_sphere 32x16 | 0.01 | 28 (stalled) | 0.071 to 0.064 | 8 to 8, 120 to 104 | 6 | 0 / 0 | 0 | 9.8e-6 | 0.003 |
| grid_cube 6 | 6.48 | 60 (max) | 32.15 to 29.98 | 8 to 8, 86 to 86 | 6 | 0 / 0 | 1 | 1.3e-2 | 0.003 |
| capsule 24x8 | 0.04 | 36 (converged) | 0.218 to 0.211 | 8 to 8, 100 to 100 | 6 | 0 / 0 | 0 | 1.3e-6 | 0.009 |
| torus 32x12 | 0.05 | 34 (stalled) | 0.405 to 0.195 | 16 to 0, 176 to 80 | 2 | 0 / 0 | 0 | 1.0e-4 | 0.003 |
| spot (quads) | 1.94 | 100 (max) | 30.40 to 4.91 | 120 to 2, 1152 to 182 | 3 | 0 / 2 | 0 | 3.9e-3 / 4.3e-3 | 0.348 / 0.343 |

The one flipped triangle on the cube is a sliver (3D area 5.8e-6, a cell is 0.11) where the cut
passes within `delta = 1e-3` of a vertex: its image has an area of `-6e-9`. The test accepts
one flipped triangle per thousand. These patches are nearly developable (the cube's faces, the
cylinder of the capsule), which is why cones off and optimal agree here; the numbers say that the
chain works for networks, not that these atlases are good ones.

### Benchmarks of the cut flow (`benches/cutting.rs`)

`cargo bench -p ddg-uv --bench cutting` (`-- cut/steps` for one group). Defaults everywhere, one
thread, `DDG_MODELS` / `DDG_DATASET` as for the tests, a missing mesh is skipped with a notice.
Every run starts as the geodesic circle of radius `0.2 sqrt(area)` around the vertex at a third
of the vertex list and takes exactly 10 flow steps (Hencky energy, `alpha_L` from the roughening
rate of the start at ripple length 0.25). One run on the development machine (criterion median
of 10 samples; two runs differed by about 10 %, more for `start` on the largest meshes). The cut
has 64 to 304 edge points after the steps, so `preview`, `describe` and `insert` are small by
construction.

| mesh (polygons / triangles) | start | 10 steps | preview | describe | insert | unwrap, cones off | unwrap, cones optimal | whole chain, cones off |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| uv_sphere 64x32 (2048 / 3968) | 4.7 ms | 23.3 ms | 15 us | 24 us | 1.0 ms | 15.5 ms | 22.4 ms | 41.7 ms |
| grid_cube 16 (1536 / 3072) | 3.7 ms | 16.2 ms | 13 us | 21 us | 0.81 ms | 11.9 ms | 15.4 ms | 32.5 ms |
| torus 96x24 (2304 / 4608) | 5.7 ms | 23.4 ms | 16 us | 26 us | 1.1 ms | 16.7 ms | 16.8 ms | 44.8 ms |
| spot quads (2928 / 5856) | 8.2 ms | 31.1 ms | 21 us | 31 us | 1.4 ms | 21.6 ms | 26.5 ms | 52.8 ms |
| spot triangles (5856 / 5856) | 8.7 ms | 39.8 ms | 21 us | 29 us | 1.9 ms | 24.5 ms | 31.1 ms | 62.8 ms |
| camel (22704 / 22704) | 60.5 ms | 161.9 ms | 72 us | 81 us | 6.1 ms | 119.5 ms | 160.0 ms | 305.7 ms |
| bunny (28576 / 28576) | 81.8 ms | 222.2 ms | 103 us | 121 us | 8.1 ms | 167.3 ms | 214.9 ms | 429.9 ms |

`start` is `Session::new`: validation, the seam cut, the maps, the flow's symbolic analysis and
the distance. "unwrap" is `ddg_uv::unwrap` of the inserted mesh with the cut as seams (the
insertion is its own column). "Whole chain" is `Session::new`, 10 steps, `describe_cut`, `insert`
and `unwrap` with cones off.

### Benchmarks of the atlas (`benches/cutting.rs`, groups `regions/`)

`cargo bench -p ddg-uv --bench cutting -- regions/`. 6 regions from `NormalClusters`,
`recommended_options`, rescaled Dirichlet energy, `alpha_L` from the ripple rule (wavelength 0.5),
10 flow steps, one thread; one run (criterion median of 10).

| mesh | cut after 10 steps | start | 10 steps | preview | describe | insert | unwrap, cones off | unwrap, optimal | whole chain, cones off |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| uv_sphere 64x32 | 196 edge points, 8 junctions, 208 paths | 7.6 ms | 62.6 ms | 47 us | 128 us | 1.6 ms | 14.5 ms | 14.8 ms | 93.6 ms |
| grid_cube 16 | 179, 6, 188 | 5.8 ms | 45.5 ms | 39 us | 113 us | 1.3 ms | 11.0 ms | 11.5 ms | 66.7 ms |
| spot quads | 447, 26, 486 | 13.5 ms | 82.3 ms | 90 us | 267 us | 3.1 ms | 26.8 ms | 26.7 ms | 126.4 ms |

`start` includes the clustering of the normals and the distance computations of the functions.

### The flow on an intrinsic triangulation

The cut flow is sensitive to badly shaped triangles. `Session::with_triangulation(.., &Triangulation)`
chooses the triangulation it runs on:

* `Triangulation::Caller`: the seam-cut triangle mesh of the caller's triangulation. What
  `Session::new` and `Session::new_with` do, bit for bit as before (tested: the energies, the
  steps and the state agree to the bit with and without the option).
* `Triangulation::Intrinsic(IntrinsicOptions)`: from that mesh (so the seams stay boundary) an
  intrinsic triangulation of the same surface (`ddg-intrinsic`, after Sharp, Soliman, Crane,
  *Navigating Intrinsic Triangulations*, 2019), flipped to Delaunay, optionally refined, exported
  as a simplicial mesh. The flow runs on the export. `Triangulation::default()` and
  `IntrinsicOptions::default()` are **Delaunay flips only** (`IntrinsicOptions::flips_only()`):
  the measurements of the intrinsic crate and of the flow on it say that flips are what makes
  badly shaped meshes behave, and that refinement (`IntrinsicOptions::refined()`: 25 degrees,
  an explicit choice, with its bound in `refine`) adds 30 % to 140 % vertices for little gain on
  most meshes. `Session::new` keeps the caller's triangulation so that existing callers do not
  change.

`Session::intrinsic_report()` gives `IntrinsicReport`: vertices of the mesh and of the
triangulation, **vertices added** (`n_added`), flips (first pass, export), refined, splits, the
refinement rounds of the export, whether it is Delaunay, the smallest and largest angle reached,
the angle bound asked for, the number of corners of the export below it, the corners the last
refinement left below the bound with the reason (`BelowReason`), seconds, the fingerprint.

**When is the per-corner view the whole state?** The vertices of the export are the vertices of
the seam-cut mesh first, then the added ones. With flips only nothing is added unless the export
had to split an edge to make the triangulation a simplicial complex (on the meshes of the table
below: 0 added on all of them). `Session::state_is_per_corner()` is true exactly then
(`IntrinsicReport::n_added == 0`, or the caller's triangulation), and `state()` / `set_state`
carry everything. Otherwise `state()` has the values at the caller's vertices only, and the
whole state is a `SessionState`.

**`SessionState`** (opaque; `Session::snapshot()`, `restore(&state)`,
`Session::new_from_state(&input, &state)`): the level-set functions on the vertices of the
triangulation the flow runs on (for the two-region flow one `φ`, otherwise one per region), the
fingerprint of that triangulation (a hash of its connectivity, lengths and places; for the
caller's triangulation a hash of the seam-cut mesh), the number of regions, the kind of flow and
the `Triangulation` option that reproduces the triangulation (refinement and export are
deterministic). `to_bytes()` / `from_bytes()` is a versioned little-endian format; `with_values()`
edits. A state of another triangulation is refused with `CuttingError::StaleState` (both
fingerprints are in the message); another number of regions or kind of flow, a wrong number of
values or bytes that are not a state is `CuttingError::BadState`.

**Starts on an intrinsic triangulation.** A geodesic circle is computed on the triangulation
itself, from every copy of the vertex. Corner values, polygon labels and normal clusters are
computed on the seam-cut mesh first; the value at a vertex the triangulation added is the linear
interpolation of the mesh's piecewise linear function at the vertex's place on the mesh. Corner
values are left as they are (at the caller's vertices exactly as given). Labels and clusters are
then made signed distances by the flow's redistancing (cut-adjacent values kept, fast marching
on the triangulation elsewhere). Per-vertex controls (`penalty`, the scaled-Dirichlet weight) are
interpolated the same way; symmetry classes leave out the added vertices (and the seam copies);
a curvature per triangle (`CurvatureSource::PerFace`) is refused on an intrinsic triangulation.
`Session::evaluate` returns the evaluation in the numbering of `Session::engine_mesh()`.

**The cut back on the caller's polygons** (`describe_cut`, `preview`). A cut on the export
(points on its edges, segments and triple points in its faces) is carried back by
`segments_on_input`: every point becomes a point of the seam-cut mesh (inside a face, or on an
edge) and every segment a chain of straight pieces that each lie in one face. In the types of
`CutOnPolygons`: a point on a polygon edge is an `EdgePoint` (an edge can carry several), a point
on a diagonal an `InteriorPoint`, a point inside a triangle a `JunctionPoint` **of any degree**
(a bend has two paths, a triple point three); a path runs between edge points and junctions and
passes through diagonal crossings; a triangle can hold several paths; consecutive points of a
path can be in the same triangle; a path can enter and leave a triangle through the same side.
The regions of the parts come from the sides of the oriented pieces (lower region id on the
left of a path), not from the regions of the vertices: `Inserted::polygon_region`. The preview
comes from the same description (`CutOnPolygons::lines`).
The rule for near-degenerate parts: every point keeps at least `IntrinsicOptions::min_offset`
(`1e-6`) of the side or of the barycentric coordinate from the vertices and sides of its
triangle (moved points are counted in `n_clamped_crossings`); a point exactly at a vertex of the
mesh is `CuttingError::Unsupported`. Refused as unsupported: a piece along a diagonal, a closed
piece of the cut inside one triangle (it would be a hole in a polygon; this happens on the
refined atlas of spot at a ripple wavelength of 0.25).

**Simplification of the cut (`IntrinsicOptions::merge_tolerance`, default 0.1).** The cut that
comes back has a bend wherever one of the flow's segments crosses an edge of the intrinsic
triangulation inside a triangle of the caller's. Such a bend is an artefact of the
triangulation, and near a vertex of the mesh the unwrap of the inserted cut flipped triangles
because of it. A bend inside a triangle (two pieces, never a junction of three) is removed
when its distance from the chord of the two pieces, added to what the pieces already stand
for, is at most `merge_tolerance` times the longest edge of the triangle, the chord is not along
a side, no other point of the triangle is in the swept triangle and no piece crosses the
chord; the swept triangle has no vertex of the mesh inside, so the cut keeps passing on the same
side of every vertex. `Session::describe_cut_with_report` returns a `DescribeReport` (bends
removed, largest displacement as a length and as a fraction of an edge, points moved by
`min_offset`). `0` switches it off.

**`insert`, generally.** The construction (written out in the rustdoc of
`ddg_uv::cutting::insert`; a bmesh implementation must mirror it): vertices on polygon edges and
on diagonals; every triangle of a crossed polygon is subdivided by the pieces it holds into
cells (planar graph in the barycentric coordinates of the triangle, faces walked
counter-clockwise); cells are joined across the parts of diagonals that are not cut into parts;
the boundary of a part is its polygon; every cell is triangulated by ear clipping in barycentric
coordinates, the best-shaped ear first. For the old cases (a segment from side to side, a triple
point) it makes the polygons of before; the triangulation can differ in a diagonal.

#### Tests and measurements

New tests (`tests/cutting.rs`, 36 in all): the general insertion by hand (a V with a bend
inside one triangle: parts, regions, area 0.04; two parallel cuts with several points on an
edge, one 1e-7 from a corner); 144 random level-set functions on the vertices of the
intrinsic triangulations (flips and refined, 2, 3 and 5 regions, sphere, cube, torus and a
cylinder whose caps are 12-gons, with and without seams) through the whole chain: 2495 bends,
904 edges with several points, 15412 paths at junctions, 1423 points on the two sides of seams,
none refused, every invariant of the earlier insertion tests (area kept, manifold, points on
their edge or diagonal, every cut segment an edge, old seams kept, junctions of degree d are
vertices of d polygons and d seams) and the regions of the parts (equal across an uncut edge,
different across a cut one); the caller's triangulation unchanged by the option (bits); the
intrinsic session against a direct run of `ddg-intrinsic` and the flow on the same mesh (same
fingerprint, same energies and state to the bit); the state round trip, `restore`,
`new_from_state`, stale state, bad bytes; end to end on the generated meshes, spot, the sheared
grid, the needles-and-caps grid (both built as in `bad_mesh_study`) and `camelbad`.

Release build, circle start of radius `0.2 sqrt(area)` (flat meshes: `alpha_L = 1`, 8 steps;
others: ripple rule, 40 steps, camelbad 20), Hencky energy, cones optimal; "energy" is the first
and the last step's:

| mesh | triangulation | vertices | energy | charts | cones | flipped | distortion |
|---|---|---:|---|---:|---:|---:|---:|
| sheared 153, 32x32 | caller / flips / refined | 1089 / 1089 / 1155 | 0.992 to 0.744 / 1.216 to 0.908 / 1.216 to 0.908 | 2 | 0 | 0 | 0 |
| needles and caps, 16 | caller / flips / refined | 192 / 192 / 830 | 1.151 to 0.797 / 1.167 to 1.167 (the stop rule ends the run after 1 step) / 1.200 to 0.821 | 2 | 0 | 0 | 0 |
| uv_sphere 32x16 | caller / flips / refined | 482 / 482 / 690 | 4.78 to 4.17 / same / 4.89 to 4.32 | 2 | 1 | 0 / 0 / 0 (5 without the simplification) | 0.074 / 0.074 / 0.073 |
| grid_cube 6 | all three | 218 | 4.30 to 4.04 | 2 | 1 | 0 | 0.100 |
| spot (quads) | caller / flips / refined | 2930 / 2930 / 3185 | 4.81 to 4.25 / 4.90 to 4.36 / 4.90 to 4.36 | 2 | 2 | 0 | 0.434 / 0.430 / 0.432 |
| spot atlas, 8 clusters, 60 steps | caller / flips / refined | 2930 / 2930 / 3185 | 0.0618 to 0.0379 / 0.0628 to 0.0401 / 0.0626 to 0.0405 | 17 / 21 / 23 | 0 | 0 / 0 / 0 (2 / 12 without) | 0.0047 / 0.0044 / 0.0042 |
| camelbad, 22725 faces | caller / flips / refined | 11392 / 11392 / 15561 | 2.456 to 2.341 / 2.475 to 2.355 / 2.458 to 2.350 | 2 / 3 / 3 | 2 | 0 / 0 / 0 (1 / 1 without) | 0.501 / 0.506 / 0.506 |

(A flat mesh has no distortion: its flattening is isometric.) The tests assert **zero flipped
triangles** for the end-to-end cases above, on all three triangulations.

**Why the intrinsic runs flipped triangles, and what removed them** (`examples/flip_study.rs`:
for every flipped triangle of the unwrap its area against the mean, its smallest angle, its
construction and how close its cut points lie to a vertex or side of the old triangle, then the
same flow state unwrapped with `min_offset` of 1e-6, 1e-4, 1e-3, 1e-2; run in parallel over the
cases). Without the simplification (`merge_tolerance = 0`) the flipped triangles are:

| case | flipped | area / mean area | smallest angle | corners |
|---|---:|---|---|---|
| uv_sphere, refined | 5 | 5e-3 to 3e-2 | 1.9 to 39 degrees | all: a vertex of the mesh, a point of an edge, a bend |
| camelbad, flips / refined | 1 / 1 | 8e-5 / 6e-5 | 0.1 / 0.3 degrees | the same |
| spot atlas (wavelength 0.5), flips / refined | 1 / 12 | 2e-3 / 1e-5 to 4e-2 | 1.3 / 0.01 to 50 degrees | the same, a few with two bends or a diagonal point |
| spot atlas (0.25), flips | 2 | 5e-4 to 1e-2 | 0.5 degrees | the same |
| spot atlas (1.0), flips / refined | 4 / 16 | 2e-3 to 0.25 | 0.3 to 15 / 0.1 to 47 degrees | the same |
| sheared, needles, cube, spot quads | 0 | | | |

Every flipped triangle is a triangle of the fan around a vertex of the mesh whose other two
corners are cut vertices, at least one of them a bend. The distance of those cut vertices from
the vertices and sides of their old triangle is not small (`8e-3` to `0.3` of an edge, most above
`0.05`), and raising `min_offset` to `1e-4`, `1e-3` or `1e-2` leaves the counts as they are (and
at `1e-2` makes them a little worse): it is not a point close to a corner. The cut that comes
back wraps around the vertex in short pieces with bends that are only there because of the
triangulation the flow ran on, and the flattening of the chart does not keep that boundary
embedded. With `merge_tolerance` of 0.01 the counts fall (sphere 5 to 0, camelbad 1 to 0,
spot atlas 0.5 flips 1 to 0, refined 12 to 4); with 0.1 they are 0 in every case measured
(atlas wavelengths 0.25, 0.5, 1.0 with flips and refined, camelbad, sphere, spot, cube, the bad
meshes). The largest displacement it reached: 0.07 of an edge on the atlas with flips, 0.10 on
the refined atlas and camelbad refined, 0.01 to 0.04 on the single circles; the bound is the
option. The caller's triangulation never needed it (no bends). The cost is in `describe` (a
few tens of microseconds more, below).

**The needles-and-caps run with flips only** (the row above, energy 1.167 over 8 steps) was
measured before a defect of the flow was fixed. Six vertices lay on the start circle to
rounding; their crossings were moved away by the smallest offset, the slope of the level-set
function in those faces came out as `2·10⁻¹²`, and the mean slope that scales the step became
`2·10⁻⁸`, so the first step moved the cut by `6·10⁻¹¹` and the stop rule called the run
converged. Faces with a moved crossing no longer count in the mean slope, and a run converges
by motion only if the predicted motion is small as well (`ddg-variational-cuts` SPEC D30,
D31). With that the same start moves from the first step: 1.167 to 0.853 in 8 steps, measured
in the flow crate's test on the same exported mesh; the row above was not re-measured through
the session.

#### Benchmarks (`-- intrinsic/`)

`cargo bench -p ddg-uv --bench cutting -- intrinsic/`, one thread, 10 steps, circle start, one
run alone, criterion medians. "added" is the number of vertices the intrinsic triangulation adds.

| mesh | triangulation | vertices (added) | angles | start | 10 steps | describe | insert | unwrap, cones off | chain |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| uv_sphere 64x32 | caller | 1986 | | 4.1 ms | 21.9 ms | 25 us | 1.4 ms | 17.0 ms | 42.7 ms |
| | flips | 1986 (0) | 5.6 to 92.8 | 8.0 ms | 22.2 ms | 90 us | 1.4 ms | 16.9 ms | 44.8 ms |
| | refined | 2633 (647) | 25.0 to 119.6 | 17.8 ms | 29.6 ms | 149 us | 1.3 ms | 16.6 ms | 66.8 ms |
| sheared 165, 64x64 | caller | 4225 | | 10.2 ms | 38.4 ms | 36 us | 2.6 ms | 38.5 ms | 87.8 ms |
| | flips | 4225 (0) | 5.1 to 159.9 | 22.8 ms | 24.3 ms | 296 us | 4.0 ms | 39.2 ms | 83.9 ms |
| | refined | 4865 (640) | 15.0 to 104.8 | 32.0 ms | 23.7 ms | 322 us | 3.8 ms | 36.2 ms | 95.4 ms |
| spot quads | caller | 2930 | | 7.1 ms | 22.6 ms | 30 us | 1.8 ms | 23.0 ms | 58.3 ms |
| | flips | 2930 (0) | 11.4 to 133.6 | 13.4 ms | 29.3 ms | 158 us | 1.8 ms | 23.6 ms | 66.5 ms |
| | refined | 3185 (255) | 25.0 to 126.8 | 21.2 ms | 32.0 ms | 181 us | 1.8 ms | 23.5 ms | 75.7 ms |

(The flow steps on the sheared grid are faster on the Delaunay triangulation, 22.6 against 39.2
ms: a better conditioned system.) The construction of the triangulation (the start column
minus the caller's) costs 3 to 13 ms on these meshes with flips and 14 to 22 with refinement. The
simplification of the cut is in `describe` (150 to 300 us on the intrinsic rows against 90 to 160
before it); the chain is within the noise of the earlier run (about 10 %).

#### What the Python binding and the Blender operators will have to change

* The state can no longer live in per-corner attributes alone: with refinement, or when the
  export had to split edges, it has values at vertices that are in no polygon. Keep the
  `SessionState` bytes (`snapshot().to_bytes()`) with the session (an object property, a custom
  data layer on the mesh datablock, a file next to the .blend) and hand them back (`restore`,
  `new_from_state`); `state_is_per_corner()` is the check that per-corner attributes suffice. A
  state is tied to the mesh and seams (fingerprint): when the mesh or the seams change the
  state is stale and must be restarted from corner values or labels.
* The description has new shapes: several edge points per edge, junctions of degree two and
  more inside triangles (`triangle`, `corners`, `barycentric`), paths that start and end at
  junctions, paths with both ends on the same polygon edge, `VertexSource::Junction`,
  `Inserted::polygon_region`.
* The bmesh insertion must mirror the general construction of `insert` (rustdoc), not the
  per-triangle cases of the earlier version: subdivide every crossed triangle by its pieces
  (cells), join cells across uncut parts of diagonals, triangulate by ear clipping. The
  triangulation it produces need not equal ours; the polygons, the vertices and the seams must.
* The preview is `CutOnPolygons::lines()` / `Session::preview` (with `regions`).
* `Session::new` is unchanged; the intrinsic triangulation is `with_triangulation`.
  `Session::flow_mesh()` is the seam-cut mesh, `Session::engine_mesh()` the triangulation the
  flow runs on.

## Benchmarks

`cargo bench -p ddg-uv` (Criterion; `benches/pipeline.rs`). End to end is `unwrap` from the
arrays in memory to packed UVs with all defaults, one thread; low / median / high of 10
samples; same machine, nothing else running.

| benchmark | input | time |
|---|---|---|
| `unwrap/bumpy_disk_64` | generated, 8 192 triangles, boundary | 83.4 / 84.1 / 85.2 ms |
| `unwrap/uv_sphere_64x32` | generated, 2 048 polygons, closed; free-site fallback, two solves | 366 / 367 / 368 ms |
| `unwrap/grid_cube_16` | generated, 1 536 quads, closed; raw optimum | 125 / 128 / 131 ms |
| `unwrap/torus_96x24` | generated, 2 304 quads, genus 1; no cone solve | 10.7 / 11.1 / 11.4 ms |
| `unwrap/hand` | dataset, 6 660 triangles | 504 / 511 / 519 ms |
| `unwrap/camel` | dataset, 22 704 triangles | 1.465 / 1.496 / 1.527 s |
| `unwrap/bunny` | dataset, 28 576 triangles, closed | 3.002 / 3.058 / 3.118 s |
| `unwrap/brain` | dataset, 305 782 triangles, closed | not run by Criterion here; 34.1 s in the table above |

Per stage (medians):

| mesh | 1 prepare (validation, seams) | 2 cones | 3 flatten | 4 verify and pack |
|---|---|---|---|---|
| `stages/bumpy_disk_64` | 2.82 ms | 68.3 ms | 11.9 ms | 1.19 ms |
| `stages/uv_sphere_64x32` | 1.29 ms | 374 ms | 5.14 ms | 0.43 ms |
| `stages/camel` | 8.91 ms | 1.414 s | 44.8 ms | 3.52 ms |

The cone solver is where the time goes; these benchmarks are what a faster solver will be
measured with:

```sh
cargo bench -p ddg-uv -- --save-baseline before      # now
cargo bench -p ddg-uv -- --baseline before           # after the change
cargo bench -p ddg-uv -- stages/                     # per-stage only
DDG_BENCH_LARGE=1 cargo bench -p ddg-uv -- brain     # 306k triangles, ten samples of about 34 s each; not run here
```

The generated meshes need no files; the dataset meshes are skipped with a notice when
`data/optimal-cones/optimal-cones-data/data` (or `DDG_DATASET`) is absent. The same per-stage
times are in every `Report` (`timings`), so `evaluate dataset` gives the breakdown for all 23
meshes in about a minute.

## Tests

```sh
export CARGO_TARGET_DIR=$PWD/target-uv
cargo test -p ddg-uv                                         # unit, pipeline, error, property tests; one dataset mesh
cargo test -p ddg-uv --release --test dataset -- --nocapture # all 23 dataset meshes
cargo test -p ddg-uv --release --test properties census -- --ignored --nocapture
cargo clippy -p ddg-uv --all-targets -- -D warnings
```

| file | what it checks |
|------|----------------|
| unit tests in `src/` | the distortion measure, the packer, the primitives |
| `tests/pipeline.rs` | end to end: a bumpy disk (quality, determinism), generated triangle meshes, the polygonal primitives, a cylinder cut along its rims and one side edge (three charts, exactly the seams cut, an isometry to `1e-12`), ignored and extended seams, steering, loose vertices, the merge guard, the free-site fallback, re-optimization |
| `tests/errors.rs` | every refusal with its exact message; inconsistent winding repaired and reported; degenerate triangles mollified and reported |
| `tests/properties.rs` | random meshes of ten kinds (disk, sphere, cube, capped cylinder, cone, torus, genus 2, disk with n-gon, annulus, two components), jittered, quads split along either diagonal, random seams, steering and options: finite UVs, cuts on polygon edges only, UVs continuous across every edge that is not cut, seams cut, prescribed cones exact or reported, counts consistent, positive cones where asked, input unchanged, identical second run |
| `tests/dataset.rs` | the 23 dataset meshes with the defaults (release build; a debug build runs `bump` only, because the unoptimized cone solver needs minutes for the others); skipped with a notice without the data |

The census of the property test's generator (2000 fixed cases): 1982 unwrapped and passed
every check, 9 refused because positive cones were demanded on a closed genus-2 surface, 9
because the random steering was infeasible; 1694 had a seam cut, 1471 had cones, 211 several
charts, 1941 steering.

Four error variants are not triggered by any test input, only their wording is checked:
`ConeSolver`, `ConesNotCertified`, `Flatten`, `Internal`.

## Not done, not verified

- The spike check flags only meshes where a handful of vertices carry 99 % of the cone-free
  distortion; weaker spikes are not flagged (see its section). Spike tips on the boundary are
  not excluded.
- Optimized boundary values (`boundary: Optimized`) are not recommended as a default, see [the
  boundary](#the-boundary): the flattening does not realize them on several dataset meshes,
  and the result then has overlap and mirrored triangles. Closed pieces, also those with an
  automatic slit, never use them; with a slit the optimized values were not tried.
- Closed surfaces of genus ≥ 2 without seams: a map comes back, its quality is not
  controlled, and it was tried on two voxel slabs only.
- The `Auto` criterion (32 cones and 1 % of the vertices) separates the cases measured, which
  lie at 100 % and at 0.1 % or less; nothing in between was measured.
- `λ` is per piece. A surface split by seams into pieces of very different size gets the same
  `λ` on each; weighting it by the area share of a piece (which would make the split objective
  equal to the unsplit one) is not implemented.
- The packer is a shelf packer of bounding boxes; fill is 0.2 to 0.7 for single charts.
- Degenerate triangles: mollification is tested on a quad with three collinear corners; no
  Blender mesh with degenerate triangles was run.
- The supplied triangulation is trusted to be a triangulation in the combinatorial sense
  (checked); a fan across a non-convex polygon gives triangles with flipped or zero area,
  which the pipeline treats like any other geometry.
- Polygons that were flipped to repair the orientation have mirrored UVs with respect to
  their own winding; the report counts them, the UVs are not post-processed.
- The Python binding and the Blender add-on expose the boundary option (`optimize_boundary`);
  their tests (`scripts/py-test.sh`, `scripts/blender-test.sh`) only check that it reaches the
  pipeline and is reported, not its quality on meshes from Blender.
- The intrinsic triangulation of the cut session: a circle law was not checked beyond the end-to-end rows; a cut through a vertex of the mesh exactly, a piece along a diagonal and a closed piece inside one triangle are refused; the `3.5` times difference of energy the study saw on curved bad meshes was not reproduced (only the flat case, see above); nothing was run on Blender meshes; the flow on flipped sheared grids was run for 8 steps only.
- The cut session (`ddg_uv::cutting`), see [the section](#the-cut-flow-on-polygon-meshes):
  no Blender mesh was used (the quad model is `spot_quadrangulated.obj`, whose polygons are
  fanned from their first corner by `PolygonMesh::from_polygons`, not Blender's tessellation),
  and `insert` was not compared with an editor's bmesh implementation. The new polygons'
  diagonals are not checked against the other edges of the mesh before `unwrap` does it
  (`DiagonalOnEdge`): a diagonal that joins two vertices that another polygon also joins is
  an error of `unwrap`, never produced in any test. A geodesic-circle start does not reach across
  a seam. Symmetry classes ignore vertices that a seam splits. The N-region flow was run for tens of steps on the atlas tests and benches, not to convergence, and its atlases were not looked at (the core README says the same of its own); `Session::flow()` of an N-region session is a placeholder. Adding `VertexSource::Junction` is the one change to an existing public enum. The cuts in
  the tests and benchmarks come from a few steps of the flow from a circle; no claim is made on
  how well they flatten, and the flow's convergence on these meshes was not studied here.

## License

MIT.
