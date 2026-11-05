# ddg-bff

Boundary First Flattening with cone singularities, and the mesh cutting it needs: given a
triangle mesh and cones (vertex, cone angle), cut the surface through the cones into disks and
flatten every disk conformally. Output: one UV per triangle corner, the cut edges, the chart of
every face, and a report per chart. Charts are not packed. A closed genus-0 mesh can also be
mapped conformally to the unit sphere (`map_to_sphere`, see *Mapping to the sphere*).

This crate is a Rust port of the method and the reference implementation of

> Rohan Sawhney and Keenan Crane. *Boundary First Flattening.*
> ACM Transactions on Graphics 37(1), 2017. <https://doi.org/10.1145/3132705>

Reference C++ code: <https://github.com/GeometryCollective/boundary-first-flattening>, MIT,
© 2017 Rohan Sawhney and Keenan Crane. Its licence notice is kept in
[`LICENSES/BFF-MIT.txt`](LICENSES/BFF-MIT.txt). `disk.rs`, `cones.rs`, `curve.rs` and
`quality.rs` follow `Bff.cpp` and `Distortion.cpp`; the cutting (`cut.rs`, `plan.rs`) is new.

Also used: intrinsic mollification from N. Sharp and K. Crane, *A Laplacian for Nonmanifold
Triangle Meshes*, SGP 2020 (as in the reference code); for the cut planner K. Mehlhorn, *A
faster approximation algorithm for the Steiner problem in graphs*, IPL 27 (1988), D. Eppstein,
*Dynamic generators of topologically embedded graphs*, SODA 2003, and J. Erickson and
K. Whittlesey, *Greedy optimal homotopy and homology generators*, SODA 2005.

It depends on `ddg-core` and on nothing else; it does not depend on the cone solver.

## Example

```rust
use ddg_bff::{Cone, ConeFlattenInput, ConeFlattenOptions, CutSource, flatten_with_cones};
use ddg_core::generators::icosphere;
use ddg_core::{Geometry, TriMesh};

let sphere = icosphere(2)?;
let mesh = TriMesh::new(sphere.positions.len(), sphere.faces.clone())?;
let geometry = Geometry::from_positions(&mesh, &sphere.positions)?;
// The 12 vertices of the icosahedron, each with a cone angle of π/3: they sum to 4π.
let cones: Vec<Cone> = (0..12)
    .map(|vertex| Cone { vertex, angle: std::f64::consts::PI / 3.0 })
    .collect();
let flat = flatten_with_cones(
    &ConeFlattenInput {
        mesh: &mesh,
        geometry: &geometry,
        cones: &cones,
        scale_factors: None,                                   // or Some(u) from a cone solver
        cut: CutSource::Plan { may_cut: None, must_cut: &[] }, // masks and seams go here
        anchors: &[],
    },
    &ConeFlattenOptions::default(),
)?;
let uv_of_corner = &flat.corner_uvs; // [3 f + k] is the corner of face f at faces[f][k]
let cut = flat.cut_edges();          // edge ids of `mesh`, ascending
let chart_of_face = &flat.face_chart;
for chart in &flat.charts {
    assert_eq!(chart.distortion.flipped_triangles, 0);
    assert!(chart.is_injective());
}
assert!(flat.seams.max_length_mismatch < 1e-12);
```

## API

| item | what it is |
|------|------------|
| `flatten_with_cones(&ConeFlattenInput, &ConeFlattenOptions) -> Result<ConeFlattening, BffError>` | the whole thing: scale factors, cut, boundary curves, harmonic extension, reports |
| `solve_scale_factors(&TriMesh, &Geometry, &[Cone], &ScaleFactorOptions) -> Result<ScaleFactors, BffError>` | `(L u)_i = µ_i − Ω_i`, `u = 0` on the boundary; Gauss–Bonnet check on closed components |
| `plan_cut(&CutPlanInput) -> Result<CutPlan, BffError>` | chooses the cut: disks, cones on the boundary, forbidden and required edges |
| `cut_mesh(&TriMesh, &[u32]) -> Result<CutMesh, BffError>` | cuts along edges; faces and corners keep their ids, `vertex_to_original` and `edge_to_original` map back |
| `component_topology(&TriMesh) -> Vec<ComponentTopology>`, `discontinuous_edges(&TriMesh, &[T]) -> Result<Vec<u32>, BffError>` | topology of every component; seams of an existing per-corner attribute |
| `DiskBff::new(&TriMesh, &Geometry, DiskOptions)` with `flatten_automatic`, `flatten_with_scale_factors`, `flatten_with_exterior_angles`, `flatten_to_disk`, `flatten_with_boundary_data` | the flattenings of one disk; also `dirichlet_to_neumann`, `neumann_to_dirichlet`, `compatible_exterior_angles`, `compatible_scale_factors` |
| `mollify_edge_lengths(&TriMesh, &[f64], f64) -> Result<Mollified, BffError>` | the reference's mollification as an explicit step |
| `ChartReport { closure, distortion, overlap, .. }`, `ConeReport`, `SeamReport`, `ClosedComponentReport`, `Timings` | what a result says about itself |
| `quality::{coverage, polygon_self_crossings, quasi_conformal_error, fit_similarity, ..}` | the measures, usable on any layout |
| `map_to_sphere(&TriMesh, &Geometry, &SphereOptions) -> Result<SphereMap, BffError>` | closed connected genus-0 mesh to the unit sphere: `SphereMap { positions, report }`; see *Mapping to the sphere* |
| `mobius_center`, `sphere_distortion`, `most_regular_star`, `stereographic`, `inverse_stereographic`, `sphere_inversion` | the steps of the sphere map, usable on their own |
| `BffError` | every refusal, with the vertex, edge, face or chart concerned |

Conventions are those of the repository's `CLAUDE.md`: `L` is the positive semidefinite cotan
Laplacian, `u` the true log scale factor (plane lengths ≈ `e^u` × surface lengths), a cone
angle is `µ = 2π − (angle sum in the plane)`, positive for a cone tip.

What a caller should know that the signatures do not say:

- **Nothing is normalized.** Every chart starts with the first vertex of its boundary loop (or
  its anchor) at the origin and its first boundary edge along `+x`, at the scale `e^u`. Charts
  therefore lie on top of each other and their relative sizes are meaningful. Packing is the
  caller's job.
- **`u = 0` on the original boundary** and nowhere else; on a closed component `u` has
  area-weighted mean zero.
- **Cones must be interior vertices**, at most one per vertex, with angle below `2π`.
- **Closed components need `Σµ = 2πχ`.** The default refuses a discrepancy above `1e-8`
  (`GaussBonnetViolated`, with the numbers). `GaussBonnetPolicy::DistributeOverCones` changes
  the cone angles in proportion to their absolute values; `SpreadOverVertices` keeps them and
  gives every vertex the same small curvature, which is then reported
  (`max_interior_curvature`, `closure.turning_error`) and lands as an angle error at the
  vertex where the chart's boundary curve starts.
- **Edge lengths are never mollified implicitly.** A mesh with a degenerate triangle is refused
  when its `Geometry` is built. To do what the reference does, call `mollify_edge_lengths` and
  build the geometry from the result; `Mollified` says how much was added to every edge and how
  many faces made it necessary — put that in your own report.
- **Seamless** means: the two sides of every cut edge have the same length in the plane
  (`seams.max_length_mismatch`, rounding level) and the angles at the copies of a vertex add up
  to `2π − µ` exactly; the two sides differ by a rigid motion, not necessarily by a grid
  symmetry.
- All solves run on the calling thread.

## The cut planner

`plan_cut` receives the mesh, a cost per edge (its length), the cone vertices, an optional
per-edge mask `may_cut` and a list `must_cut` of edges that have to be cut (seams).

1. The mesh is cut along the required edges. They may split it into several pieces; every
   piece becomes one chart and is not split further.
2. In each piece, every boundary loop and every cone that is still interior is a terminal. A
   closed piece of positive genus without cones gets one vertex as base point.
3. One Dijkstra run from all terminals over the cuttable edges, then Kruskal over the regions:
   Mehlhorn's approximation of a Steiner forest joining all terminals.
4. A maximum spanning tree of the dual graph that avoids the forest; the cuttable edges in
   neither tree close the handle loops (tree–cotree, `2g` loop edges per piece of genus `g`).
5. Cut edges that end at an interior vertex which is not a cone are pruned.
6. A cut that would be a single edge between two interior vertices cannot be represented by an
   indexed mesh (its two sides would be the same vertex pair); it is extended by one edge.
7. The mesh is cut and the result is checked.

**Guarantees** when `plan_cut` returns `Ok`:

- every chart is a disk (connected, Euler characteristic 1, one boundary loop) and every copy
  of every cone is on the boundary of its chart — checked at run time on the mesh that was
  actually cut;
- one chart per piece left by the required cuts — checked at run time;
- no edge with `may_cut == false` and no boundary edge is cut, every required interior edge is
  cut — by construction, checked by the tests;
- the result is deterministic (ties are broken by ids).

**Refusals** name their cause, and the tests verify each stated cause independently:
`ConeEnclosedByUncuttableEdges`, `UncuttableRegionNotDisk`, `IsolatedCutEdge` (a one-edge slit
that cannot be extended), `RequiredEdgeForbidden`, and `SphereNeedsTwoCones`. The last one is
a design decision, not a topological fact: a closed surface of genus 0 can be cut along any
path, but the planner does not invent the ends of a cut.

What the tests establish (`tests/planner.rs`): sixteen surfaces — closed of genus 0, 1, 2 and
4; the same with one to four holes; disks with holes; an annulus; a mesh with two components —
each with 0 to 6 cones and 12 seeds; 640 instances with random forbidden and required edges;
and 1440 instances on nine meshes small enough to try every subset of the cuttable edges. On
those 1440 the planner found a cut exactly when one exists (1063 planned, 265 refused where
none exists), apart from 112 refusals of the `SphereNeedsTwoCones` kind. That it always does
so is not proved. The cut is short, not shortest: on the dataset below its length is 0.885 to
1.000 of the cut stored with the reference's output.

**Negative cones.** The planner gets vertices, not angles, and treats all cones alike. A cone
with `µ < 0` has more than a full turn around it; see the next section for what that costs.

## Flipped triangles and overlap

A chart report has two exact signs of overlap — `overlap.wedges_over_2pi` (boundary vertices
with more than a full turn around them) and `overlap.boundary_crossings` — and an estimate of
how much is concerned: `overlap.coverage`, the area covered at all and the area covered more
than once, by point sampling (16 points per average triangle by default; a triangulation
without overlap gives exactly zero). `ChartReport::is_injective()` combines the exact signs
with the absence of flipped and degenerate triangles.

Measured on the 23 meshes of the optimal-cones dataset with their published cones and the
planner's cut (`examples/dataset_report.rs`; `ovl` is the share of the layout's footprint that
is covered more than once):

| mesh | cones | negative | most negative | flipped | wedges > 2π | boundary crossings | ovl |
|------|------:|---------:|--------------:|--------:|------------:|-------------------:|----:|
| axes, bump, bunny, cow1, cow5e-1, cow75e-2, face | 2–5 | 0 | – | 0 | 0 | 0 | 0.00 % |
| camel | 6 | 2 | −2.72 | 0 | 0 | 0 | 0.00 % |
| camel_noise1 | 7 | 2 | −2.84 | 19 | 0 | 0 | 0.00 % |
| camelshotnoise | 6 | 2 | −2.65 | 72 | 0 | 0 | 0.00 % |
| camel_noise2 | 7 | 1 | −0.06 | 113 | 0 | 0 | 0.03 % |
| camelbad | 6 | 2 | −2.70 | 1 | 1 | 1 | 0.00 % |
| camelbadnonuniform | 6 | 2 | −2.82 | 1 | 1 | 1 | 0.01 % |
| brain | 45 | 11 | −3.16 | 38 | 1 | 1 | 0.01 % |
| chair | 14 | 2 | −2.07 | 0 | 1 | 3 | 1.14 % |
| airplane | 13 | 5 | −3.49 | 4 | 0 | 4 | 2.17 % |
| mannequin | 8 | 1 | −1.73 | 2 | 1 | 1 | 5.55 % |
| tree | 8 | 4 | −1.68 | 0 | 2 | 2 | 6.61 % |
| cow01 | 439 | 244 | −0.33 | 0 | 48 | 134 | 7.77 % |
| hand | 29 | 11 | −5.51 | 0 | 4 | 6 | 9.28 % |
| rhino | 36 | 19 | −3.50 | 0 | 5 | 7 | 10.05 % |
| body | 10 | 2 | −4.74 | 0 | 1 | 7 | 14.92 % |
| cow00001 | 13836 | 7821 | −0.20 | 0 | 2362 | 15470 | 15.42 % |

What drives it:

- **Overlap comes from negative cones.** The seven meshes whose cones are all positive have no
  flipped triangle, no crossing and no sample covered twice. Of the sixteen with negative
  cones, nine have 1 % to 15 % of their footprint covered more than once. Every wedge over
  `2π` sits at a negative cone: at all 2417 negative cones that end a cut, and at 10 through
  which the cut passes unevenly.
- The UV maps stored in the dataset (output of the reference implementation) overlap just as
  much: 15.65 % for body, 17.93 % for cow00001, 9.57 % for hand, 18.87 % for airplane with its
  cut (ours: 2.17 %). This is a property of the method with these cones, not of the port.
- **Flipped triangles come from the triangulation.** They occur in 8 meshes, 1 to 113 of them,
  with at most `3.4e-4` of the layout area; of the 250 flips, 106 are next to an edge with a
  negative cotan weight, 3 at a cone, the rest elsewhere in meshes where 6 % to 17 % of the
  edges have negative weight. The stored maps have the same flips (brain 35 against our 38,
  camel_noise2 113 and 113, camelshotnoise 70 and 72).
- With positive cones the boundary can still cross itself in slivers: where a zig-zag cut from
  a small positive cone runs back towards the cone, the slit opens with negative width. Seen
  with cones from the cone solver (bunny 6 crossings, mannequin 30): 0.000 % to 0.024 % of the
  footprint. Decide by `coverage`, not by `is_injective()`.

What a caller can do about negative cones:

1. **Do not ask for them**: bound the cone angles from below by zero in the cone solver.
   Then no wedge exceeds `2π`.
2. **Pass seams through them** with `must_cut`, from boundary to boundary or cone to cone, so
   that no wedge holds more than `2π`. This makes more charts. (Extending the cut beyond a
   negative cone to an ordinary vertex does not help: the two sides of a slit that ends at a
   flat vertex coincide in the plane.)
3. **Accept and report**: `ConeReport::max_wedge_angle` and `overlap.coverage` say where and
   how much.

A scratch experiment with the cone solver of this repository (commit `b25c588`, cluster
merging on; the experiment is not part of this crate, which must not depend on the solver)
gave, with our `λ`:

| `λ` | negative cones without the bound | overlap without the bound | overlap with `µ ≥ 0` | cost of the bound |
|-----|----------------------------------|---------------------------|----------------------|-------------------|
| 0.5, 0.1 | none on 22 of 23 meshes | ≤ 0.03 % | the same cones | none |
| 0.02 | airplane 54 of 60, axes 4 of 50, body 116 of 125, hand 32 of 47, the five cow files 1 of 6 | 1.24 %, 5.35 %, 12.85 %, 2.62 %, 0.00 % | 0.00 % on these, ≤ 0.03 % on the others | solver's distortion 0.138 → 0.826, 0.143 → 0.399, 0.238 → 0.635, 0.210 → 0.679, 0.2678 → 0.2704 |

(The 23rd mesh, camelshotnoise, gets 420 to 589 cones and a distortion of `1e4` from the
solver at all three values; that is a matter for the solver.)

## Verification

### Against the reference implementation

`reference/bff_ref.cpp` is a driver around the reference library that writes, in full
precision, the positions it used, the UV of every corner, the cut it chose and where its
boundary walk starts. `reference/build.sh` builds it (the reference sources are compiled where
they are; nothing of them is in this repository), `reference/run.sh` runs it and compares
(`examples/reference_diff.rs`): the same mesh, cones and cut go through this crate, and the
two layouts are compared after the best similarity per chart, relative to the diagonal of the
layout's bounding box.

| input | cases | this crate with the reference's shift `1e-8` (rms / max) | this crate with exact solves (rms / max) |
|-------|------:|-----------------------------------------------------------|------------------------------------------|
| dataset meshes, published cones | 20 | 3.1e-15 … 8.7e-14 / ≤ 9.1e-13 | 2.9e-8 … 9.3e-6 / ≤ 2.1e-5 |
| camelbad, camelbadnonuniform (mollified by the reference) | 2 | 5.0e-12, 6.1e-12 / 2.2e-10 | 1.1e-6 / ≤ 3.2e-6 |
| camelshotnoise | 1 | 4.6e-10 / 1.3e-9 | 1.5e-4 / 3.9e-4 |
| generated spheres and disks with 0 to 25 cones of both signs | 7 | 1.2e-15 … 1.6e-14 / ≤ 4.5e-14 | 5.0e-8 … 1.2e-6 / ≤ 2.4e-6 |
| five disks of the reference, unit disk | 5 | 6.2e-15 … 1.0e-11 / ≤ 4.4e-10 | 2.3e-7 … 3.1e-6 / ≤ 3.2e-5 |
| the same, automatic flattening (holomorphic extension) | 5 | 4.2e-12 … 5.5e-8 / ≤ 4.9e-7 | 2.3e-7 … 3.1e-4 / ≤ 6.3e-4 |

With the reference's diagonal shift the two programs compute the same thing and agree to
rounding. The numbers in the last column are the effect of that shift on the result. In the
last row the agreement stops at `5e-8`: the holomorphic extension depends on the direction of
the boundary curve in the plane (see *Deviations*), and the shift makes the exterior angles
miss `2π` by `ε Σu ≈ 1e-4`, which turns the reference's curve relative to ours.

### The dataset

`data/optimal-cones/optimal-cones-data/data` holds 23 meshes with the UV map the reference
produced (six significant digits) and the cones "used to generate" it. Two comparisons are
tested (`tests/dataset.rs`) and tabulated (`examples/dataset_report.rs`,
`examples/dataset_investigate.rs`): with the cut that is visible in the stored map, our layout
against the stored one; with our own cut, the lengths of all triangle edges in the plane
against the stored ones after one global scale.

*Reference, run today* is the driver above with the published cones; it agrees with this crate
to the precision of the previous table, so that column says how well the files reproduce at
all. *Unaccounted turning* is measured on the stored map alone: the sum of the cone angles it
realizes at the vertices of its cut that are not listed as cones, plus, for airplane and
camelshotnoise, the amount by which it misses the published angle at one and two cones. A map
computed by Boundary First Flattening from the published cones has none.

| mesh | reference, run today, against the stored map (rms) | this crate, stored cut, against the stored map (rms) | unaccounted turning of the stored map | own cut: edge length ratio to the stored map, 1 % / 99 % | the same between our two cuts |
|------|---------:|---------:|---------:|-----------------|-----------------|
| face | 1.8e-6 | 1.9e-6 | −1.9e-4 | 0.9997 / 1.0003 | 0.9997 / 1.0002 |
| bump | 2.3e-6 | 2.3e-6 | 4.3e-4 | 0.9998 / 1.0001 | 1.0000 / 1.0000 |
| axes | 4.0e-6 | 3.8e-6 | 9.8e-5 | 0.9990 / 1.0009 | 0.9992 / 1.0008 |
| camel_noise2 | 4.6e-6 | 4.7e-6 | 1.3e-4 | 0.9916 / 1.0084 | 0.9916 / 1.0085 |
| cow1 | 1.2e-5 | 1.3e-5 | 3.3e-4 | 0.9938 / 1.0058 | 0.9949 / 1.0044 |
| hand | 1.7e-5 | 1.7e-5 | 2.4e-5 | 0.9846 / 1.0108 | 0.9845 / 1.0107 |
| cow75e-2 | 1.8e-5 | 1.8e-5 | −2.1e-4 | 0.9937 / 1.0064 | 0.9954 / 1.0041 |
| cow5e-1 | 1.9e-5 | 1.4e-5 | −4.7e-4 | 0.9984 / 1.0012 | 0.9994 / 1.0006 |
| mannequin | 2.3e-5 | 2.3e-5 | −7.5e-4 | 0.9961 / 1.0060 | 0.9960 / 1.0060 |
| cow01 | 2.4e-5 | 2.4e-5 | −7.2e-6 | 0.9992 / 1.0008 | 0.9992 / 1.0008 |
| bunny | 4.1e-5 | 4.2e-5 | 7.7e-4 | 0.9981 / 1.0011 | 0.9990 / 1.0011 |
| camelbad | 4.1e-5 | 4.1e-5 | −1.6e-3 | 0.9654 / 1.0222 | 0.9648 / 1.0216 |
| chair | 4.2e-5 | 4.3e-5 | −7.9e-4 | 0.9960 / 1.0048 | 0.9960 / 1.0049 |
| camel | 4.3e-5 | 4.3e-5 | −1.6e-3 | 0.9739 / 1.0267 | 0.9731 / 1.0262 |
| camelbadnonuniform | 4.6e-5 | 4.6e-5 | −3.1e-3 | 0.9738 / 1.0419 | 0.9730 / 1.0408 |
| tree | 5.8e-5 | 5.8e-5 | −8.8e-4 | 0.9993 / 1.0008 | 1.0000 / 1.0000 |
| brain | 1.3e-4 | 1.3e-4 | 2.7e-3 | 0.9970 / 1.0059 | 0.9975 / 1.0051 |
| cow00001 | other cut | 2.1e-4 | 9.1e-3 | 0.9983 / 1.0017 | 1.0000 / 1.0000 |
| rhino | 7.9e-4 | 7.9e-4 | −1.4e-2 | 0.9744 / 1.0229 | 0.9936 / 1.0110 |
| camel_noise1 | 2.5e-3 | 2.5e-3 | −3.4e-4 | 0.9590 / 1.0587 | 0.9921 / 1.0139 |
| camelshotnoise | 5.5e-3 | 5.4e-3 | −0.32 | 0.9505 / 1.1429 | 0.9609 / 1.0326 |
| airplane | 6.7e-3 | 6.7e-3 | −0.11 | 0.9102 / 1.0670 | 0.9375 / 1.0541 |
| body | other cut | 9.6e-3 | −0.105 | 0.9635 / 1.3143 | 0.9995 / 1.0006 |

How to read it:

- **This crate and the reference miss every stored map by the same amount.** What the table
  shows is how well the files reproduce from the published data, not how good the port is.
- **Positions, angles and UVs are stored with six digits.** Rounding the input like that
  changes our result by `1e-7` to `2e-6` (`7e-5` for camelshotnoise). Only face, bump, axes and
  camel_noise2 reproduce that well.
- **The stored maps turn by more than the published cones say**, and the residual grows with
  it: in 20 of the 23 meshes it is between 0.005 and 0.1 times the unaccounted turning. (On
  the four closed meshes — axes, bunny, chair, brain — that turning necessarily equals minus
  the amount by which the published angles miss Gauss–Bonnet.) This fits maps that were
  computed from cone data with small angles at further vertices which the published lists
  leave out: when the exterior angles do not sum to `2π`, the reference puts the whole
  deficit at the vertex where its boundary walk starts, and that vertex is on the cut. It is
  an inference from the files; what the maps were computed from is not recorded.
  For airplane and camelshotnoise the turning sits at listed cones (the stored map has −3.497
  instead of −3.387 at one cone, and −2.720 and −2.716 instead of −2.461 and −2.653 at two);
  with those angles the layouts agree to `1.4e-4` and `1.1e-4`
  (`stored_maps_with_other_cone_angles_than_published`).
- **Not explained:** camel_noise1 (and hand, cow01 at a much lower level) do not fit that
  pattern. In the stored maps of camel_noise1 and rhino the scale along the original boundary
  varies by a factor 3.04 and 1.43 (in ours and the reference's: 1.19 and 1.006), so their
  boundary curves were closed with a far larger correction than these cones need.
- **The flat metric does not depend on the cut, but a layout does**, by its discretization
  error: the last column compares our own two layouts. Where the fifth column is as wide as
  the sixth (camel, hand, the two camelbad), the difference from the stored map with another
  cut is that dependence and nothing else.
- *body* and *cow00001* are the two meshes for which the reference, run today, does not choose
  the cut of the stored map. Vertex numbers below count from 0.
- **cow00001** (13836 cones on 15049 vertices). The cuts differ in 4 edges of 14232: one leads
  to a cone whose published angle is `0.0`, so its two sides coincide in the stored UVs and
  cannot be seen; three are exchanged for edges of exactly the same length (0.018824,
  0.014801, 0.018437): ties. The maps agree: all edge length ratios lie between 0.9972 and
  1.0022, which is the rounding of the file for edges this short, for the reference and for
  this crate. (The `2.1e-4` in the table is for the stored cut completed by an edge of the
  planner's choice.)
- **body** is mirror symmetric with its cones in mirror pairs. The stored cut contains a path
  of 51 edges between the cones at vertices 2685 and 11646; the reference today takes its
  mirror image between 5617 and 11824 instead. The two paths have the same length (0.738661):
  a tie, which the reference breaks by memory address and this crate by vertex and edge ids.
  The stored cut is a valid cut through the published cones. The stored *map* is not
  reproduced on either cut: on the stored cut our layout misses it by `9.6e-3`, and the
  reference on its mirror cut misses its edge lengths by 0.958 / 1.306, while our own layouts
  on different cuts agree to 0.9995 / 1.0006. The stored map has an unaccounted turning of
  −0.105 (summed over 430 vertices, where single values are lost in the rounding of UV edges
  as short as `1e-4`), the largest after the two meshes with other cone angles, and it
  contains edges of length zero.

### Other tests

- `tests/flatten.rs`: isometric unfoldings (flat tetrahedron, octahedron and icosahedron with
  their own angle defects as cones; a flat annulus, a cylinder, a flat torus, a flat disk with
  three holes); angle sums `2π − µ` at all vertices; closed genus 2 with four cones of `−π`;
  surfaces of genus 1 and 2 with holes; two components with anchors; seams between charts;
  Gauss–Bonnet policies; overlap at negative cones; every error; no NaN.
- `tests/properties.rs`: property tests of cutting, planning and flattening on generated
  meshes.
- `tests/planner.rs`: see *The cut planner*.

## Timings

`flatten_with_cones` with the planner's cut, one thread, best of three, milliseconds (AMD
Threadripper 3970X; `examples/dataset_report.rs`). Between runs of the table the totals vary
by about 10 % (brain: 1.00 to 1.24 s over five runs). `report` includes about 100 ms per
300 000 faces for the overlap estimate, which `overlap_samples_per_triangle: 0.0` switches off.

| mesh | faces | cones | total | scale factors | cut | curves | extension | report |
|------|------:|------:|------:|--------------:|----:|-------:|----------:|-------:|
| airplane | 4656 | 13 | 9.5 | 2.9 | 1.6 | 0.1 | 2.4 | 2.6 |
| bump | 6528 | 5 | 14 | 5.1 | 2.1 | 0.1 | 4.0 | 2.7 |
| hand | 6660 | 29 | 15 | 5.3 | 2.5 | 0.1 | 4.3 | 3.2 |
| axes | 7466 | 4 | 18 | 6.4 | 2.6 | 0.1 | 5.3 | 3.4 |
| chair | 11448 | 14 | 29 | 12 | 3.9 | 0.2 | 8.0 | 5.0 |
| mannequin | 12195 | 8 | 32 | 13 | 4.3 | 0.2 | 8.6 | 5.9 |
| camelbadnonuniform | 18125 | 6 | 45 | 16 | 6.3 | 0.2 | 14 | 8.4 |
| camel, camel_noise1, camel_noise2, camelbad, camelshotnoise | 22704–22968 | 6–7 | 59–61 | 19–22 | 7.7–8.1 | 0.3 | 19–22 | 11–12 |
| body | 23892 | 10 | 59 | 22 | 8.7 | 0.4 | 17 | 11 |
| face | 27150 | 2 | 70 | 23 | 9.4 | 0.4 | 26 | 11 |
| bunny | 28576 | 4 | 79 | 29 | 11 | 0.5 | 27 | 12 |
| cow01, cow1, cow5e-1, cow75e-2 | 29999 | 2–439 | 82–97 | 31–37 | 11–14 | 0.4–0.5 | 23–31 | 13–15 |
| cow00001 | 29999 | 13836 | 79 | 31 | 15 | 2.6 | 1.7 | 29 |
| tree | 33592 | 8 | 90 | 32 | 13 | 0.6 | 30 | 15 |
| rhino | 40304 | 36 | 111 | 46 | 15 | 0.6 | 31 | 19 |
| brain | 305782 | 45 | 1139 | 433 | 134 | 4.8 | 412 | 155 |

Two factorizations dominate: the Laplacian of the surface (scale factors) and the Dirichlet
Laplacian of the cut surface (extension). Passing `scale_factors` from a cone solver saves the
first.

## Mapping to the sphere

`map_to_sphere` maps a closed, connected surface of genus 0 conformally to the unit sphere. It
follows `BFF::mapToSphere` of the reference (`Bff.cpp`) and its Möbius centering `BFF::centerMobius`,
whose algorithm is Baden, Crane and Kazhdan, *Möbius Registration*, CGF 37(5) (SGP 2018),
Algorithm 1. The code is in `src/sphere.rs`. The method:

1. Remove the star of one vertex (the pole): what is left is a disk (a submesh; the caller's mesh
   is not touched). `SphereOptions::pole`: `MostRegularStar` (default) or `Vertex(i)`.
2. `DiskBff::flatten_to_disk` (`disk_iterations = 10`), then normalize the disk (mean boundary
   radius 1).
3. Choose the radius `R` of the disk by ternary search of `ln R` in `[1, 1e4]` (tolerance
   `1e-3`), maximizing the area-weighted mean of `ln(a_f / A_f)` (the reference's objective).
4. Inverse stereographic projection of the disk of radius `R`, from the south pole; the pole
   vertex is placed at `(0, 0, -1)`.
5. Möbius centering: inversions of the sphere (Gauss–Newton step `c = -J⁻¹µ`) until the center of
   mass `µ` of the surface's area measure on the sphere has `|µ| <= 1e-10`. A step longer than
   `0.5` is shortened, a step that does not decrease `|µ|` is halved.

`SphereMap::report` (`SphereReport`) carries the radius search and centering reports
(`CenteringReport::converged()`), timings, and `SphereDistortion`: quasi-conformal error (mean,
max), area distortion (rms of `u = ½ ln(a/A)`), flipped triangles (negative orientation seen from
outside), the covering degree (signed area of the image over `4π`) and `is_injective()` (no
flipped triangle and degree 1).

What a caller should know:

- **Errors, not repairs.** `NotASphere` for a mesh that is open, disconnected or of another genus
  (map components one by one with `TriMesh::split_components`); `PoleOutOfRange`; options out of
  their domain are `NonFiniteInput`.
- **Orientation:** the image is oriented like the input, seen from outside.
- **The map is determined up to a rotation** (centering fixes the inversion, not the rotation);
  the tests compare after a best rotation fit. Rigid motions and scaling of the input do not
  change it.
- **The pole star is where quality is worst.** Next to the removed star the quasi-conformal error
  is larger (icospheres: max about 1.24, mean about 1.00); the star is not refined.
- **Flipped triangles are reported, not corrected** (brain, 153k vertices: 38, as in the
  reference's map).
- `radius_search: None` projects the normalized disk (`R = 1`, northern hemisphere);
  `centering: None` leaves the pole at `(0, 0, -1)`.

### Measurements

`examples/sphere_report.rs`, one thread, run alone, best of 3, milliseconds (AMD Threadripper
3970X). QC: quasi-conformal error `σ_max/σ_min` per triangle; area rms: rms of the log area scale
factor `u`. The last column is the same with the reference's pole (vertex 0), where it differs
from ours.

| mesh | V | QC mean / max | area rms `u` | flipped | ms | pole 0: QC mean / max, flipped |
|---|---:|---|---:|---:|---:|---|
| icosphere 3 | 642 | 1.0096 / 1.240 | 0.006 | 0 | 4 | same |
| icosphere 5 | 10242 | 1.0008 / 1.244 | 0.002 | 0 | 80 | same |
| ellipsoid 2:1:0.5, noisy | 2562 | 1.0746 / 1.332 | 0.500 | 0 | 19 | 2.0658 / 28.0, 4 |
| spot | 2930 | 1.1332 / 1.780 | 1.022 | 0 | 18 | 1.3028 / 2.451, 0 |
| blub | 7106 | 1.1039 / 1.843 | 2.173 | 0 | 51 | 1.2364 / 2.088, 0 |
| bob | – | refused: genus 1 | | | | |
| bunny | 14290 | 1.0410 / 1.289 | 1.631 | 0 | 121 | 1.0500 / 1.308, 0 |
| axes | 3735 | 1.1219 / 1.243 | 1.827 | 0 | 26 | |
| chair | 5726 | 1.0544 / 1.653 | 0.560 | 0 | 43 | 1.2777 / 1.960, 0 |
| brain | 152893 | 1.0529 / 71.9 | 0.459 | 38 | 1554 | 38 flips too |

### Against the reference

`reference/bff_ref.cpp` has a mode `sphere` (`reference/run_sphere.sh`); the driver is
`examples/reference_sphere_diff.rs`. Not run by `cargo test`.

- Our disk on the punctured mesh equals the reference's own disk mode to `2e-15` … `6e-13` rms.
- The reference's in-place puncture inside `mapToSphere` differs from that by `2e-7` … `1e-4`
  rms (unexplained; a property of its code).
- With both maps centered the same way (ours), the sphere maps agree to `3e-6` (icosphere),
  `8e-4` (bunny), `5e-3` (axes), `5e-2` (brain, where both have 38 flips).
- The reference's own centering (image areas) leaves `|µ|` between `0.03` and `0.6` by the
  surface measure; its area distortion is higher (axes: `u` rms 2.41 against our 1.83).

### Tests and proofs

`tests/sphere.rs` (9 tests): an icosphere maps to itself up to a rotation; without centering the
pole is at the south pole; centering puts the area center at the origin; genus-0 models map
injectively; invariance under rigid motion and scale; different poles give the same map up to a
rotation; non-spheres and bad options are refused; the formulas below hold numerically. Lean:
`lean/ddg/DDG/Conformal/Sphere.lean` (stereographic projection on the sphere, its inverse,
conformality and orientation; inversions keep the sphere; the derivative `2(u - (u·x)x)`).

## Deviations from the reference implementation

1. **Exact solves.** The reference adds `1e-8` to the diagonal of every Laplacian to make it
   definite. Here the Dirichlet problems are solved as they are, and singular ones (pure
   Neumann; closed components) with a compatible right-hand side and one pinned vertex.
   *Reason:* the shift perturbs results by `1e-7` to `1e-4` of the layout size (table above)
   and, on a closed component whose cone angles miss Gauss–Bonnet by `d`, adds the constant
   `d / (1e-8 · n)` to `u`, which rescales the layout by its exponential.
   `diagonal_shift: REFERENCE_DIAGONAL_SHIFT` reproduces the reference.
2. **No implicit mollification.** The reference lengthens all edges whenever a triangle misses
   the triangle inequality by less than `1e-6` of the mean edge length, which also happens for
   valid thin triangles (camelbad, camelbadnonuniform: 5 and 4 faces, `2e-8` of the mesh
   radius, layout changed by `1.5e-7` rms, `6e-6` max). *Reason:* inputs are not changed
   silently. `mollify_edge_lengths` does it on request and reports it.
3. **Gauss–Bonnet is checked.** The reference solves whatever it is given. Here the default
   refuses; the two policies are explicit. *Reason:* see 1.
4. **Free constants.** Where a solution is determined up to a constant (`u` on a closed
   component, `u` from prescribed angles, the conjugate coordinate) the reference takes what
   its shifted system returns, roughly mean zero over the vertices. Here `u` has area-weighted
   mean zero, the choice of least area distortion. Only the size of the layout depends on it
   (unit-disk maps differ in scale by 0.99 to 1.53 from the reference's).
5. **Several boundary loops are cut, not filled.** The reference fills all holes but the
   longest boundary loop with fans of triangles (all of them if the longest loop is shorter
   than half the diameter) and flattens the filled surface. Here loops are joined by cuts, and
   `u = 0` on every original boundary loop. *Reason:* no geometry is invented, and the
   boundary conditions are the same on every loop.
6. **Handles are part of the cut.** The reference computes generators from breadth-first
   spanning trees, turns them into a boundary, and treats that boundary like an original one.
   Here the handle loops come from the same shortest-path and cotree computation as the cuts
   through cones and are cut edges like any other: equal lengths on both sides, exact angles.
7. **The cut planner** (new): forbidden and required edges, several charts, deterministic
   tie-breaking, verification of the result, errors with causes. The reference's `Cutter` is a
   similar approximate Steiner tree without these; its ties are broken by memory address.
8. **Seams between charts.** Lengths are shared across cut edges between different charts as
   well, by closing the boundary curves of coupled charts together (the reference has one
   boundary curve per mesh).
9. **Closing the curve.** Same least-squares problem (paper Eq. 19, 20). Added: constraints
   that are redundant (the two sides of all cut edges parallel, as on a flat torus) are
   detected instead of inverted; a change of any length by more than `max_length_adjustment`
   (default one half) is not made — the weakest closing direction is dropped, the remaining
   gap is spread along the curve and reported (`ClosureReport`). The reference changes lengths
   by whatever closes the curve, including to negative values.
10. **One-edge slits.** A cut of a single edge between two interior vertices is extended by an
    edge, because an indexed triangle mesh cannot represent it; the reference's wedge-based
    mesh can.
11. **Orientation and placement.** The boundary is walked counter-clockwise (the reference:
    clockwise, then mirrored). The layout is left where the curve was built — first boundary
    vertex at the origin, first edge along `+x` — instead of being recentred at the centroid
    of its corners, rotated by a quarter turn, aligned to principal axes and packed.
12. **Where the boundary starts** is an option (`anchors`, `boundary_start`); in the reference
    it follows from the order of its data. With harmonic extension this is a rigid motion.
    With the holomorphic extension (`Extension::Conjugate`, used by `flatten_automatic` as in
    the reference) it is not: the first coordinate is extended and the second conjugated, and
    turning the curve changes the map by `1e-4` (rms) to `1e-3` (max) of its size on the
    meshes tried. This is a property of the method; the reference has it too.
13. **Non-finite results are errors**, as is a layout all of whose faces collapse. The
    reference replaces NaN coordinates by zero.
14. **Validation.** Cones on boundary vertices, duplicate cones, angles of `2π` or more,
    options that are not numbers: refused.
15. **Area distortion** is reported as `u = ½ ln(B/A)` per face, the log scale factor; the
    reference reports `ln(B/A) = 2u`.
16. **Not ported:** flattening to a given target curve, cone placement (the reference's own
    greedy cone picker), direct editing of boundary lengths, bin packing, the viewer.
17. **Mapping to the sphere** (`map_to_sphere`) differs from `BFF::mapToSphere` in:
    - *Pole.* The reference removes vertex 0. Here the default is the vertex with the most
      regular star (`most_regular_star`, edge lengths only); `PoleChoice::Vertex(0)` is the
      reference's choice. *Reason:* a poor pole gives flipped triangles (noisy ellipsoid: 4 and a
      maximum error of 28 against none and 1.33).
    - *Puncture.* The reference punctures its mesh in place; here a submesh is built and the
      input is untouched. The two disks differ by `2e-7` to `1e-4` (see above).
    - *Centering.* The reference weights by the areas of the image triangles, with tolerance
      `1e-3` and 1000 iterations. Here the weights are the surface areas (as in Baden et al.),
      the centroids are normalized, the tolerance is `1e-10`, the step is capped at `0.5` and
      halved when `|µ|` does not decrease. Reported as `CenteringReport`.
    - *Radius search.* Same objective; the disk is normalized first and the search is over `ln R`,
      so the result does not depend on the scale of the input.
    - *Projection from the south pole*, the pole vertex at `(0, 0, -1)`.
    - *Errors* for meshes that are not closed connected genus 0 (the reference does not check).
    - `diagonal_shift` defaults to `0` as everywhere else in this crate.

The reference code itself deviates from the paper in two places, and so does this port: the
unit-disk iteration does not average successive iterates (paper Sec. 6.5), and the closing
problem weights lengths by the boundary edge lengths `ℓ_ij` (Eq. 19) rather than by the dual
lengths (Alg. 5, line 7).

## Formula tags

Functions that implement a formula carry `[formula: <id>]` in their documentation
(`scripts/check-formulas.sh`). Ids that already exist in `lean/proofs` are reused:
`angle-defect`, `boundary-geodesic-curvature`, `cotan-laplacian`, `lumped-mass`,
`euler-characteristic`, `genus-from-euler`, `closed-mesh-cone-sum`,
`closed-mesh-optimal-constant`, and `discrete-yamabe-linearization` for `(L u)_i = µ_i − Ω_i`.

New ids, with exact statements. `B`: boundary vertices in loop order with the surface on the
left; edge `i` runs from boundary vertex `i` to `i + 1`. BFF = Sawhney and Crane 2017.

| id | statement | source | function |
|----|-----------|--------|----------|
| `bff-dirichlet-to-neumann` | For a source `φ` and boundary values `g`: solve `L_II a_I = φ_I − L_IB g`; the Neumann data are `h = φ_B − L_BI a_I − L_BB g`. | BFF Eq. 14, Alg. 3 | `DiskBff::dirichlet_to_neumann` |
| `bff-neumann-to-dirichlet` | Solve `L a = φ − [0; h]` (solvable iff `Σφ = Σh`, unique up to a constant); the boundary values are `a_B`. | BFF Eq. 5, Alg. 4 | `DiskBff::neumann_to_dirichlet` |
| `bff-cherrier-boundary-curvature` | With `(L u)_I = −Ω_I` and Neumann data `h` of `u`: target exterior angles `k̃ = k − h`. | BFF Eq. 11, 12 | `compatible_exterior_angles`, `compatible_scale_factors`, `cone_loop_data` |
| `bff-exterior-angle-sum` | For a disk, the angles `k̃ = k − h` of the previous line sum to `2π`. | BFF Proposition 1 (App. A) | `compatible_exterior_angles` |
| `bff-cone-angle-sum` | Cut through cones, with `(L u)_i = µ_i − Ω_i` on the uncut surface and `h` evaluated with the Laplacian of the cut surface: the interior angles `π − k̃` at the copies of vertex `i` sum to `2π − µ_i`. | BFF Proposition 2 (App. A), Sec. 6.4 | `cone_loop_data` |
| `bff-target-boundary-length` | `ℓ*_i = e^{(u_i + u_{i+1})/2} ℓ_i`. | BFF Sec. 5 step II, Alg. 1 line 9 | `target_boundary_lengths`, `cone_loop_data` |
| `bff-dual-boundary-length` | `ℓ_j = ½ (ℓ_{j−1} + ℓ_j)` at boundary vertex `j` (half of its two boundary edges). | BFF Sec. 3.3 | `dual_boundary_lengths` |
| `bff-boundary-tangents` | `φ_0 = 0`, `φ_i = φ_{i−1} + k̃_i`, `T_i = (cos φ_i, sin φ_i)`. | BFF Sec. 4.5.1, Alg. 5 lines 1–5 | `boundary_tangents` |
| `bff-closing-lengths` | Minimize `½ Σ_d N_d (ℓ̃_d − ℓ*_d)²` subject to `Σ_{i ∈ loop} ℓ̃_{d(i)} T_i = 0` for every loop, `d(i)` the length unknown of edge `i` (shared by the two sides of a cut edge), `N_d = Σ 1/ℓ_i` over its sides: `ℓ̃ = ℓ* − N⁻¹ Tᵀ (T N⁻¹ Tᵀ)⁺ T ℓ*`. | BFF Eq. 19, 20; Sec. 6.4 | `best_fit_curves` |
| `bff-boundary-curve` | `γ_0 = 0`, `γ_p = Σ_{i<p} ℓ̃_i T_i`. | BFF Sec. 4.5.1, Alg. 5 lines 9–11 | `integrate_curve` |
| `bff-harmonic-extension` | `L_II a_I = −L_IB γ` for each coordinate of the curve. | BFF Sec. 4.4, Eq. 6, Alg. 6 line 1 | `PoissonSystem::solve` |
| `bff-hilbert-transform` | The conjugate `b` of `a` solves `(L b)_j = ½ (a_{j−1} − a_{j+1})` at boundary vertex `j` and `(L b)_i = 0` inside. (The paper prints `½ (a_k − a_i)` for a clockwise walk.) | BFF Eq. 17, Alg. 6 lines 2–3 | `DiskBff::hilbert_transform` |
| `bff-uniformization-angles` | `k̃_j = 2π ℓ̃_j / Σ ℓ̃` with the dual target lengths of the last iterate. | BFF Sec. 6.5 | `DiskBff::flatten_to_disk` |
| `bff-close-exterior-angles` | `k̃_j ← k̃_j + (ℓ_j / Σ ℓ)(2π − Σ k̃)` with dual lengths `ℓ_j`. | reference `BFF::closeCurvatures` (BFF Sec. 6.2: angles "normalized to 2π") | `curve::close_exterior_angles` |
| `intrinsic-mollification` | `ε = max(0, max over faces and sides c of (ℓ_c − ℓ_a − ℓ_b) + δ)`, `δ = 10⁻⁶ ·` mean edge length; every length becomes `ℓ + ε`. | reference `BFFData::computeEdgeLengths`; Sharp and Crane 2020, Sec. 4.5 | `mollify_edge_lengths` |
| `cut-euler-characteristic` | Cutting `c` interior edges and thereby adding `a` vertex copies: `χ' = χ + a − c`. | counting | `cut_euler_characteristic` |
| `tree-cotree-disk` | On a connected surface with every boundary loop contracted to a point: if `T` is a set of interior edges without cycles and `T*` a spanning tree of the dual graph that crosses no edge of `T`, then cutting every interior edge not crossed by `T*` gives a disk, and the number of edges in neither tree is `2g` when `T` spans. | Eppstein 2003 | `plan_pieces` |
| `cut-leaf-pruning` | Gluing a cut edge back at an interior vertex without another cut edge turns a disk into a disk, unless these two boundary edges are its whole boundary; conversely for adding one. | elementary | `plan_pieces`, `extend_isolated_slits` |
| `quasi-conformal-error` | For the affine map of a triangle with Jacobian `J = [J_s J_t]`, `a = J_s·J_s`, `b = J_s·J_t`, `c = J_t·J_t`: `σ² = ½ (a + c ± √((a − c)² + 4b²))`, `Q = σ_max / σ_min`. | BFF Sec. 7 (Fig. 11); Sander et al. 2001 | `singular_values`, `quasi_conformal_error` |
| `face-log-scale-factor` | `u = ½ ln(B / A)` for a face of area `A` on the surface and `B` in the plane. | dataset README; reference `computeAreaScaling` (which returns `2u`) | `log_scale_factor` |
| `signed-area-2d` | `2A = (b − a) × (c − a)`, positive for counter-clockwise `a, b, c`. | – | `signed_double_area` |
| `signed-angle-2d` | The angle at `p` from `b` to `c`: `atan2((b − p) × (c − p), (b − p)·(c − p)) ∈ (−π, π]`. | – | `signed_corner_angle` |
| `segment-proper-crossing` | Open segments `a_0a_1`, `b_0b_1` cross in one interior point iff `b_0`, `b_1` are strictly on opposite sides of line `a_0a_1` and `a_0`, `a_1` strictly on opposite sides of line `b_0b_1`. | – | `proper_crossing` |
| `disk-injectivity-criterion` | An orientation-preserving piecewise linear map of a triangulated disk to the plane covers a point as often as the image of the boundary winds around it; a boundary polygon without proper crossings and without a vertex of more than a full turn winds at most once, so the map is injective up to boundary contact. | degree of a map; for convex boundaries Floater 2003 (BFF Sec. 7.5) | `ChartReport::is_injective` |
| `similarity-fit` | Least squares `q ≈ a p + b` over complex `a`, `b`: `a = Σ (q_i − q̄) conj(p_i − p̄) / Σ |p_i − p̄|²`, `b = q̄ − a p̄`. | – | `fit_similarity` |
| `stereographic-projection` | `p = (2w, \|w\|² − 1) / (\|w\|² + 1)` for `w ∈ R²` (from the south pole `(0, 0, −1)`); inverse `w = (x, y) / (1 + z)`. | elementary | `inverse_stereographic`, `stereographic` |
| `stereographic-conformal` | The projection is conformal and orientation preserving (proved in Lean). | elementary | `inverse_stereographic` |
| `mobius-sphere-inversion` | `η_c(x) = (1 − \|c\|²)(x + c)/\|x + c\|² + c`, `\|c\| < 1`, maps the unit sphere to itself. | Baden, Crane, Kazhdan 2018, Eq. 2 | `sphere_inversion`, `mobius_center` |
| `mobius-centering-jacobian` | The derivative of `∫ η_c(x) dA` at `c = 0` is `J = 2 Σ A_f (I − C_f C_fᵀ)`, `C_f` unit face centroids. | Baden, Crane, Kazhdan 2018, Alg. 1 | `mobius_center` |
| `mobius-centering-existence` | A centering inversion exists, unique up to rotation. Cited, not proved. | Baden, Crane, Kazhdan 2018, Sec. 4.1 | `mobius_center` |
| `spherical-triangle-area` | `tan(E/2) = det(a, b, c) / (1 + a·b + b·c + c·a)`. Cited. | van Oosterom, Strackee 1983 | `sphere_distortion` |
| `sphere-injectivity-criterion` | No flipped triangle and degree 1 imply a homeomorphism. Cited. | Milnor 1965, §5 | `SphereDistortion::is_injective` |

None of the BFF ids above (up to `similarity-fit`) has a theorem in `lean/proofs` or an entry in
`lean/proofs/CITED.md` yet, so `scripts/check-formulas.sh` lists them as missing. The sphere ids
have theorems in `DDG/Conformal/Sphere.lean` or entries in `CITED.md`.

## Not done, not verified

- No chart packing, no rescaling, no rotation of charts.
- Nothing is done about negative cones beyond reporting (see above).
- That the planner finds a cut whenever one exists is tested, not proved.
- `flatten_to_disk` and prescribed data on very elongated disks: the scale factors span more
  than double precision resolves in one layout (a strip 60 times as long as wide: `|u|` up to
  47, half of the faces collapse; counted in `degenerate_triangles`).
- Target curves, cone placement and direct length editing of the reference are not ported.
- Mapping to the sphere: injectivity is reported, not guaranteed (the brain has 38 flipped
  triangles, as the reference's map); the cited facts (existence of the centering, the
  injectivity criterion) are not proved; the `2e-7` to `1e-4` difference between the reference's
  in-place puncture and our submesh is not explained; the reference was run only through the
  `sphere` mode of our `reference/bff_ref.cpp`, on the meshes of the table.
- Why the stored maps of the dataset differ from what the public reference code computes
  from the published cones is inferred, not known; camel_noise1 is not explained at all.
- The driver of the differential test needs the reference sources, SuiteSparse and OpenBLAS;
  it is not run by `cargo test`.

## Reproducing

```sh
export CARGO_TARGET_DIR=$PWD/target-bff
cargo test -p ddg-bff                                             # 64 tests
cargo test --release -p ddg-bff --test dataset -- --nocapture     # all 23 dataset meshes
cargo clippy -p ddg-bff --all-targets -- -D warnings
cargo run --release -p ddg-bff --example dataset_report           # tables of this README
cargo run --release -p ddg-bff --example dataset_investigate      # what the stored maps contain
cargo test --release -p ddg-bff --test sphere                     # mapping to the sphere
cargo run --release -p ddg-bff --example sphere_report            # the sphere table
crates/ddg-bff/reference/run_sphere.sh /tmp/bff-ref               # sphere maps against the reference
crates/ddg-bff/reference/build.sh /tmp/bff-ref                    # needs cmake, ninja, OpenBLAS
crates/ddg-bff/reference/run.sh /tmp/bff-ref                      # comparison with the reference
```

The dataset is expected at `data/optimal-cones/optimal-cones-data/data` and the reference
sources at `../ddg-course/boundary-first-flattening` (or `BFF_REF_SRC`). Without the dataset the
dataset tests print a notice and pass. A debug build checks four of the dataset meshes, a
release build all of them.

## License

MIT. The flattening code is a port of the reference implementation, © 2017 Rohan Sawhney and
Keenan Crane, MIT: see [`LICENSES/BFF-MIT.txt`](LICENSES/BFF-MIT.txt).
