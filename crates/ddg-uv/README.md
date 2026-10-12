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

`ddg-core` is used with `default-features = false`: nothing GPL-derived is linked, which is
why the triangulation of the polygons is an input and not computed here.

Contents: [example](#example) · [API](#api) · [stages and policies](#stages-and-policies) ·
[the four decisions](#the-four-decisions) · [the dataset](#the-dataset) ·
[benchmarks](#benchmarks) · [tests](#tests) · [not done, not verified](#not-done-not-verified)

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
| `UnwrapOptions { lambda, cone_signs, max_cone_angle, closed_surfaces, merge, orientation, degenerate, packing, overlap_samples_per_triangle, reoptimize_angles, require_certified_cones }` | every policy, see below; `validate` |
| `ConeSigns`, `ClosedSurfaces`, `SlitPlacement`, `MergeClusters`, `Orientation`, `Degenerate`, `PackOptions` | the policy types |
| `prepare`, `solve_cones`, `flatten`, `assemble` (`Prepared`, `Piece`, `ConeStage`, `PieceCones`, `FlatStage`, `PieceFlat`) | the stages as separate functions |
| `pack_charts(&mut uvs, polygon_offsets, polygon_charts, &PackOptions) -> Result<PackReport, UvError>` | the packer on its own |
| `area_distortion(positions, corner_vertices, triangles, uvs) -> AreaDistortion` | the dataset README's distortion measure of any UV map |
| `Report` (`ValidationReport`, `SeamReport`, `DegenerateReport`, `SteeringReport`, `ChartReport`, `ConeSolveReport`, `FlattenReport`, `ConeStats`, `ClosedHandling`, `MergeOutcome`, `Timings`, `PackReport`), `Report::summary() -> Vec<(&str, Value)>` | the report, typed and as a flat key–value list for bindings |
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
them, both sides keep their length. The connected pieces are the charts, numbered by their
lowest polygon.

- A seam that is not an edge of a polygon (a diagonal, two unrelated vertices, a loose vertex)
  is ignored, counted and warned about. A seam on the boundary of the mesh needs no cut.
- A seam of a single edge between two interior vertices that no other seam touches has no two
  sides in an indexed mesh. It is extended by the shortest neighbouring polygon edge, counted
  (`n_extension_edges`) and warned about.

**3. Cones, then flattening, per piece.** Steering is translated at the vertices that are
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

**Warnings** are raised, besides the cases above, when cones sit at `max_cone_angle`, when
more than 1 % of the triangles of a chart are mirrored, and when more than 1 % of a chart's
footprint is covered twice.

## The four decisions

Everything in this section is measured with `examples/evaluate.rs`; the logs of the runs
quoted are `runs_tmp/uv-eval-*.log`. "Distortion" is the measure of the dataset's README,
`Σ A u²` with `u = ½ ln(B/A)` per triangle, with the areas normalized to sum 1 and minimized
over the scale of the whole map (`AreaDistortion::scale_free`); see `src/measure.rs` for why
the scale has to be taken out. It is about twice the solver's `½ Σ m u²`.

### Closed surfaces without seams

`evaluate closed`: positive cones, merging without size limit. Cells are
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
- **`camelshotnoise`** (132 one-vertex spikes) is the one result that is not usable: 22 cones,
  all at the angle bound, distortion 53, 17.9 % overlap, 14 742 of 22 968 triangles mirrored.
  The report carries three warnings saying so. The solver's README explains the cause (the
  distortion of a spike in the linearized model) and a remedy the pipeline does not apply
  (zero distortion weight at the spike tips).
- **Time**: the cone stage is 96.6 % of the total over the 23 meshes (60.6 of 62.7 s); the
  brain takes 33.2 of its 34.1 s there.

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

- No warning yet for "a few vertices carry almost all of the distortion of the cone-free
  map" (spikes). What exists instead: the warnings for cones at the angle bound, mirrored
  triangles and overlap, which fire on `camelshotnoise`; the result there is not usable.
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
- Nothing here was run inside Blender or through the Python binding.

## License

MIT.
