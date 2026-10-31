# ddg-py

The Python module `ddg_tools`: a mesh goes in as flat buffers, a UV per polygon corner comes out.
It is the bridge between Blender (`blender/` in this repository) and the Rust crates of the
workspace, built with [PyO3](https://pyo3.rs) for the stable ABI, so one binary loads in every
Python from 3.11 on. It never needs numpy.

> `ddg_tools.unwrap` runs the real pipeline, `crates/ddg-uv`: marked seams are cut, optimal cone
> singularities are computed per piece (`crates/ddg-optimal-cones`), the surface is cut through
> the cones and flattened conformally (`crates/ddg-bff`). Defaults: `lambda_=0.1`, positive cones
> only. `ddg_tools.PIPELINE` is `"ddg-uv"`. `ddg_tools.IS_PLACEHOLDER` is `False`: the flag lets
> a build that links a stand-in pipeline say so (the add-on then labels itself), and no build
> does. The pipeline's behaviour, policies and measurements are in `crates/ddg-uv/README.md`.

## Python API

```python
ddg_tools.unwrap(
    positions, corner_vertices, polygon_starts, triangles, *,
    seam_edges=None,
    pinned_vertices=None,
    prescribed_vertices=None, prescribed_angles=None,
    forbidden_vertices=None,
    lambda_=0.1,
    positive_cones_only=True,
    optimize_boundary=False,
    split_overlaps=False, split_max=32, split_threshold=0.001,
    cone_feedback=True, feedback_threshold=20.0, feedback_max_rounds=4,
) -> ddg_tools.UnwrapResult

ddg_tools.pack_charts(uvs, polygon_starts, polygon_charts, *, margin=0.0)
    -> (array('f'), dict)
```

`lambda` is a Python keyword, hence `lambda_`. Options are keyword-only.

### Inputs

Every array argument is any object with the buffer protocol: `array.array`, `bytes`,
`bytearray`, `memoryview`, a numpy array. It must be C-contiguous, and either flat or
two-dimensional with the row width given below. Nothing is modified.

| argument | per | width | element type | meaning |
|----------|-----|-------|--------------|---------|
| `positions` | vertex | 3 | real | vertex coordinates |
| `corner_vertices` | polygon corner ("loop") | 1 | index | the vertex of the corner |
| `polygon_starts` | polygon | 1 | index | the polygon's first corner: Blender's `MeshPolygon.loop_start` |
| `triangles` | triangle | 3 | index | **corner** indices of Blender's tessellation: `MeshLoopTriangle.loops` |
| `seam_edges` | edge | 2 | index | vertex pairs; the surface is cut along these edges first |
| `pinned_vertices` | vertex | 1 | index | vertices where a cone costs nothing |
| `prescribed_vertices` | vertex | 1 | index | vertices with a fixed cone angle |
| `prescribed_angles` | vertex | 1 | real | those angles, radians, each below 2π |
| `forbidden_vertices` | vertex | 1 | index | vertices where no cone may be placed |

Element types:

- **real**: float32 (`'f'`) or float64 (`'d'`).
- **index**: a signed or unsigned integer of 16, 32 or 64 bits (`'h' 'H' 'i' 'I' 'l' 'L' 'q' 'Q'`).
  Values must be non-negative and below 2³² − 1.
- **raw bytes** (`bytes`, `bytearray`, any buffer of one-byte items): read as native-endian
  float32 where a real is expected and as native-endian int32 where an index is expected. A
  one-byte buffer is never read as a list of small numbers: do not pass a `uint8` array of
  indices, its bytes would be taken as int32 data. Use a 16-bit or wider type.
- A byte-order prefix `@`, `=` or the machine's own (`<` on x86-64 and ARM64) is accepted; data in
  the other byte order is rejected.

Rules for the mesh arrays (violations raise `InputError`):

- Polygons own consecutive corners: polygon `p` has the corners
  `polygon_starts[p] .. polygon_starts[p+1]` (the last one up to `len(corner_vertices)`), in
  winding order. `polygon_starts[0]` is 0 and every polygon has at least 3 corners.
- There are exactly `corners − 2 · polygons` triangles, and the triangles of each polygon are a
  triangulation of it: wound like the polygon, every side of the polygon in exactly one
  triangle, every diagonal in exactly two. Blender's `mesh.loop_triangles` satisfies this.
  The Rust side does not triangulate (`ddg-core` has no triangulator). Passing Blender's own
  triangles means polygons are split exactly as Blender splits them, and the UVs come back for
  the original polygon corners; any other triangulation of the polygons works as well.
- Positions and angles are finite. Indices are in range.
- `seam_edges`, `pinned_vertices` and `forbidden_vertices` are sets: order, orientation and
  duplicates do not matter. A vertex may be in at most one of the pinned, prescribed and
  forbidden sets, and at most once in `prescribed_vertices`.
- `lambda_` is finite and not negative. It is the trade-off `λ` of the paper between area
  distortion and total cone angle: dimensionless, applied to every connected piece with its own
  area normalized to 1; larger means fewer cones. Default 0.1, the largest value at which the
  cone set of the clean dataset meshes with boundary no longer grows (the measurements are in
  the `ddg-uv` README: 0.002 to 2 were tried). On closed surfaces of genus 0 it has no effect.
- `positive_cones_only=True` (default) maps to `ConeSigns::PositiveWhereFeasible`: cone angles
  are `>= 0` on pieces with a boundary (this includes every seam-cut piece) and on closed
  pieces of genus 0 and 1. A closed piece of genus 2 or more cannot have only positive cones
  (the angles must sum to `2πχ < 0`); it gets cones of either sign and a warning.
  `False` allows both signs everywhere; negative cones lower the distortion and can make the
  map overlap itself.
- `optimize_boundary=False` (default): the scale factor is zero along the boundary of every piece
  (the mesh's boundary and the seams), so boundary lengths, and the two sides of a seam, keep
  their lengths. `True` lets the cone solver optimize it (paper Sec. 8.3): lower distortion on
  some meshes, but the two sides of a seam get different lengths (`max_seam_length_mismatch`),
  and where the optimized values are large the flattening does not realize them
  (`max_boundary_length_error`) and the map can overlap itself or have mirrored triangles; on
  8 of 19 dataset meshes with boundary it does (`crates/ddg-uv/README.md`, "The boundary").
  Closed pieces are not affected.
- `split_overlaps=False` (default): charts are left as the cones and seams made them. `True`
  cuts every chart that overlaps itself by more than `split_threshold` of its footprint across
  the layers and flattens again, keeping a cut only if it lowers the overlap without mirroring a
  triangle, at most `split_max` cuts (`crates/ddg-uv/README.md`, "Splitting overlapping charts").
  The report gains `n_splits`, `n_splits_rejected`, `n_split_no_cut`, `split_runs`,
  `split_exhausted`, `split_cut_length`, `overlapped_fraction_before_split`; the added cuts are in
  `cut_edges`.
- `cone_feedback=True` (default): on every piece the cone solver's linearized log scale factor
  is compared with the one the flattening realizes (on a piece with a boundary first without
  cones), the distortion weight of the vertices where they differ by more than
  `feedback_threshold` (log units; the tips of needles, where the linearized value is off by
  thousands) goes to zero, and the piece is solved and flattened again, at most
  `feedback_max_rounds` solves in all (`crates/ddg-uv/README.md`, "The cone feedback loop").
  On every measured mesh without needles nothing changes; it costs about a quarter more time on
  meshes with a boundary and nothing on closed ones. `False` solves once per piece. The report gains `cone_feedback_rounds`, `n_cone_feedback_reweighted`,
  `n_cone_feedback_pending`, `cone_feedback_max_mismatch`, `cone_feedback_final_mismatch` and
  `cone_feedback_vertices`.
- Steering acts on vertices that are interior to their piece: a pinned vertex is a cone site
  that costs nothing, a forbidden vertex gets no cone, a prescribed vertex gets exactly the given
  angle. On a mesh boundary or along a seam there are no cones: pinned and forbidden vertices
  have no effect there (`n_pinned_applied` and `n_forbidden_applied` count the ones that did),
  and a prescribed cone there is not applied and raises a warning.

What is not checked at this level, because it is the pipeline's business: manifoldness,
orientation, whether a seam edge is an edge of the mesh, degenerate triangles. The pipeline
reacts as follows (`MeshError` where it refuses):

- non-manifold edges and vertices, repeated vertices in a polygon, duplicate polygons:
  `MeshError` with the offending indices;
- inconsistent winding is repaired (the minority of each connected component is flipped,
  `n_flipped_faces` counts them, a warning says so); the UVs of a flipped polygon are mirrored
  with respect to its own winding;
- vertices that no polygon uses are ignored and counted (`n_unreferenced_vertices`);
- a seam that is not a polygon edge (a diagonal, two unrelated vertices) is ignored, counted
  (`n_seam_edges_ignored`) and warned about; a seam on the mesh boundary needs no cut;
- triangles without area: if there are any, all edge lengths are grown by the smallest common
  amount that restores the triangle inequalities (mollification: `mollified`,
  `mollify_epsilon`, `n_zero_area_triangles`, a warning).

### Outputs

`UnwrapResult` has read-only attributes. All arrays are `array.array`; they work with Blender's
`foreach_set`, with `memoryview` and with `numpy.frombuffer` / `numpy.asarray`.

| attribute | type | contents |
|-----------|------|----------|
| `uvs` | `array('f')`, 2 per polygon corner | the UV of corner `c` is `uvs[2c], uvs[2c+1]`, in corner (loop) order. Computed in float64, rounded to float32. **Not packed**: every chart has its own origin and the scale the flattening gave it, so charts overlap until `pack_charts` lays them out. |
| `polygon_charts` | `array('i')`, 1 per polygon | the chart of the polygon. A chart is a connected piece of the surface after the seams are cut; it is flattened as one piece, and the cuts through its cones stay inside it. Numbered by the lowest polygon of each chart. |
| `cone_vertices` | `array('i')` | vertices that carry a cone, ascending. Cones are interior to a chart; none is on a boundary or seam. |
| `cone_angles` | `array('d')` | their cone angles in radians: `2π − (angle sum in the flat metric)`, positive for a cone tip |
| `cut_edges` | `array('i')`, 2 per edge | the edges the surface was cut along, as vertex pairs `min, max`, sorted. Includes the seams that were honoured and the cuts the pipeline made through the cones. Cuts follow polygon edges only, never a diagonal of the triangulation. |
| `report` | `dict` | plain values, see below |

### The report

Values are `bool`, `int`, `float` or `str`; `warnings` is a `list[str]`. The dict is ordered:
the binding's own keys, the pipeline's keys (`ddg-uv`'s `Report::summary`), the binding's
measurements, then `warnings`. The keys of a run on a small grid are listed here in order;
`overlapped_fraction` and `max_chart_overlapped_fraction` are present only when the overlap was
estimated (it is, by default).

| key | meaning |
|-----|---------|
| `pipeline`, `placeholder` | which pipeline ran (`"ddg-uv"`), and whether it is a stand-in (`False`) |
| `n_vertices`, `n_polygons`, `n_corners`, `n_triangles` | size of the input |
| `n_seam_edges`, `n_pinned`, `n_prescribed`, `n_forbidden` | the steering as the Rust side received it (after removing duplicates) |
| `lambda`, `positive_cones_only`, `optimize_boundary`, `cone_feedback` | the options as received |
| `n_charts`, `n_cones`, `n_cut_edges` | size of the result |
| `n_unreferenced_vertices` | vertices that no polygon uses (ignored) |
| `n_flipped_faces` | polygons whose winding was repaired |
| `n_zero_area_triangles`, `mollified`, `mollify_epsilon` | triangles without area, whether edge lengths were grown to remove them, and by how much |
| `n_seam_edges_cut`, `n_seam_edges_on_boundary`, `n_seam_edges_ignored`, `n_seam_extension_edges` | seams that were cut, seams on the mesh boundary (nothing to cut), seams that are not polygon edges (ignored), and edges added to extend a one-edge seam |
| `n_pinned_applied`, `n_prescribed_applied`, `n_prescribed_not_applied`, `n_forbidden_applied` | steering that took effect |
| `n_closed_pieces`, `n_slits` | closed pieces after the seams, and automatic slits made in them |
| `n_negative_cones`, `total_cone_angle`, `max_cone_angle`, `min_cone_angle` | cone statistics (radians) |
| `cones_certified`, `cone_solver_max_residual`, `cone_solver_max_gap`, `cone_solver_factorizations` | the solver's optimality certificate and effort |
| `boundary_values_optimized`, `max_boundary_length_error` | whether any chart used optimized boundary values, and the largest relative difference between a boundary edge's length in the layout and the length the scale factors ask for |
| `max_seam_length_mismatch`, `mean_seam_length_mismatch`, `max_other_cut_length_mismatch` | relative difference of the lengths of the two images of a seam edge (largest, mean), and of the other cut edges (largest); about 1e-13 for `optimize_boundary=False` away from cones at the angle bound |
| `n_merges_rejected` | cone clusters whose merge the guard rejected |
| `area_distortion`, `area_distortion_unscaled` | the dataset README's distortion measure `Σ A u²` of the map, with and without taking the scale out |
| `flatten_flipped_triangles`, `flatten_collapsed_triangles`, `max_quasi_conformal_error` | the flattener's own checks |
| `cone_free_distortion`, `n_spike_vertices`, `spike_vertices` | the spike check: the distortion the surface would have without cones (largest over the charts, area normalized to 1), the number of spike vertices, and their numbers in your mesh, as text (`"11491, 11492, ..."`, at most 64; only present if there are spikes). A thin spike (a vertex pulled far out of the surface) dominates the distortion measure of this method and the map is then poor; the check fires when a handful of vertices (at most 64 and 1 % of the piece) carry 99 % of the distortion without cones, and adds a warning that suggests repairing the mesh. Measured: it fires on the dataset's `camelshotnoise` and on none of 22 other meshes (`crates/ddg-uv/README.md`). |
| `overlapped_fraction`, `max_chart_overlapped_fraction` | estimated share of the footprint that is covered twice, over all charts and for the worst chart |
| `packed` | `False`: the binding switches the pipeline's packer off; use `pack_charts` |
| `cone_feedback_rounds`, `n_cone_feedback_reweighted`, `n_cone_feedback_pending`, `cone_feedback_max_mismatch`, `cone_feedback_final_mismatch`, `cone_feedback_vertices` | only with `cone_feedback=True`: the most solves on a chart, the number of vertices whose distortion weight was lowered and of those the loop found but left (round limit, or more than 1 % of a piece), the largest mismatch `\|δ\|` (log units) among the vertices looked at in the first and in the last round, and the reweighted vertices in your mesh as text (at most 64; only present if there are any) |
| `seconds_stage_validate`, `seconds_stage_seams`, `seconds_stage_cones`, `seconds_cone_solver`, `seconds_spike_check`, `seconds_stage_flatten`, `seconds_cone_feedback`, `seconds_stage_pack`, `seconds_stage_verify`, `seconds_unwrap` | the pipeline's wall-clock times per stage (`seconds_cone_solver` and `seconds_spike_check` are parts of `seconds_stage_cones`; `seconds_cone_feedback` is the feedback loop, 0 when it is off) |
| `n_degenerate_triangles` | triangles with no area in 3D (left out of the measurements) |
| `n_collapsed_uv_triangles` | triangles with area in 3D but none in UV |
| `n_flipped_uv_triangles` | triangles whose UV image is mirrored |
| `measured_triangles` | triangles the following figures are computed over |
| `measured_angle_error_max`, `measured_angle_error_mean` | corner angle in UV minus corner angle in 3D, absolute, radians; the mean is weighted by 3D area. Zero when angles are preserved exactly (an isometry or a similarity); a conformal map with cones keeps them only approximately. |
| `measured_log_scale_std`, `measured_log_scale_min`, `measured_log_scale_max` | per triangle, `u = ½ ln(area in UV / area in 3D)`, the log of the local length scale; reported relative to its area-weighted mean (standard deviation, smallest, largest). Zero for an isometry. `exp()` of min and max is the range of local scale relative to the average. |
| `seconds_pipeline`, `seconds_measure`, `seconds_decode`, `seconds_validate`, `seconds_encode`, `seconds_total` | wall-clock time: the pipeline, the measurements, copying the buffers, validating them, building the result, and all of it |
| `warnings` | things the caller should know that did not stop the unwrap, e.g. `1 seam edges are not edges of a face and were ignored (the first is (0, 6))` |

The `measured_*` figures and the three triangle counts are computed by the binding from the
final UVs, independently of the pipeline. They are a check on the whole chain, including the
corner order.

### Errors

| exception | raised when |
|-----------|-------------|
| `ddg_tools.InputError` | an argument is malformed. The message starts with the parameter name: `corner_vertices[3] = 8 is out of range: there are 8 vertices (valid indices 0..=7)` |
| `ddg_tools.MeshError` | the arrays are fine but the mesh cannot be unwrapped. The message uses the caller's vertex and polygon ("face") indices: `edge (12, 57) is shared by 3 faces: 4, 9, 31` |
| `ddg_tools.SolveError` | the mesh was accepted but no valid map came out: the cone solver or the flattening failed, a cone could not be realized, or a result failed the pipeline's own final checks. Also raised if the pipeline returns malformed data. |
| `ddg_tools.NotAvailableError` | a valid request the library does not implement yet (`CutSession` with more than two regions) |
| `TypeError` | an array argument has no buffer protocol, or an option has the wrong type |

The module exceptions derive from `ddg_tools.DdgError`.

### `pack_charts`

Lays the charts out in the unit square and returns a new `array('f')` plus a dict with
`n_charts`, `scale`, `margin` and `coverage`.

1. Every chart gets the axis-aligned bounding box of its UVs. Charts are translated and scaled by
   **one common factor**; they are never rotated or mirrored.
2. Every box is padded by `margin` (a fraction of the unit square) on all sides.
3. The padded boxes go on shelves: tallest first (ties: lower chart index), left to right from
   the lower left corner, a new shelf above when the next box does not fit the remaining width.
4. The scale is the largest for which the shelves fit, found by bisection.

A single chart ends up with its longer side spanning exactly `[margin, 1 − margin]`. With `n`
charts the margin is capped at `1 / (4 ⌈√n⌉)` so that padding cannot eat the whole square; the
dict says which margin was used. The layout is deterministic and wastes space compared with a
real packer.

## `reoptimize_cones`: the best angles for given cone positions

```python
ddg_tools.reoptimize_cones(
    positions, triangles, cone_vertices, *,
    boundary="fixed",
    distortion_weights=None, lower_bounds=None, upper_bounds=None,
    prescribed_vertices=None, prescribed_angles=None,
    method="auto",
    lambda_=0.1,
) -> ddg_tools.ReoptimizeResult
```

Given a triangle mesh and the vertices where cones are to sit, it returns the cone angles that
give the least distortion of a conformal flattening. This is the canonical problem of the
optimal-cones solver with the cone cost `0` at the given vertices and `+∞` at every other
interior vertex, so there is no `λ` and no sparsity term: only the positions are fixed. The
method: if no angle bound is active, the optimality conditions are a linear system, which is
solved directly with one sparse factorization (`method="direct"`); for a few cones the reduced
system of the free sites is cheaper (`"auto"` picks); the generic solver (`"solve"`) takes over
when a bound is active. Every result carries the optimality certificate of the problem
(`crates/ddg-optimal-cones`). On a closed component the angles sum to `2πχ`.

Arguments (arrays are read like those of `unwrap`):

| argument | meaning |
|----------|---------|
| `positions` | 3 reals per vertex |
| `triangles` | 3 **vertex** indices per triangle: a manifold, consistently oriented mesh in which every vertex is used. Several connected components are allowed; each is solved as part of the whole problem |
| `cone_vertices` | the vertices where cones are placed; their angles are optimized (a result angle may be 0). Interior vertices only, each at most once, none with a prescribed cone. May be empty on a mesh with boundary |
| `boundary` | `"fixed"` (default): `u = 0` on the boundary. `"free"`: the boundary scale factor is optimized too. `"convex"`: optimized, with a convex boundary. `"penalized"`: optimized, with a cost on the boundary curvature (weight `lambda_`) |
| `distortion_weights` | one real `>= 0` per vertex (default 1) |
| `lower_bounds`, `upper_bounds` | one real per vertex, `<= 0` and `>= 0` (`-inf` / `inf` allowed); they bound the cone angles |
| `prescribed_vertices`, `prescribed_angles` | cones with a fixed angle (radians), parallel arrays |
| `method` | `"auto"` (default), `"direct"`, `"solve"` |
| `lambda_` | with `boundary="penalized"` the boundary curvature `ν` costs `lambda_·\|ν\|` (finite, `>= 0`, default `0.1` like `unwrap`); no effect with the other treatments |

`ReoptimizeResult` (read-only; arrays are `array.array`):

| attribute | contents |
|-----------|----------|
| `cone_vertices` | `array('i')`, as given |
| `cone_angles` | `array('d')`, parallel to `cone_vertices`, radians (positive for a cone tip, may be 0) |
| `log_scale` | `array('d')`, the log conformal scale factor `u`, one per vertex |
| `boundary_curvature` | `array('d')`, one per vertex: the curvature of the flattened boundary with an optimized boundary, 0 elsewhere |
| `report` | `dict`: `status` (`"converged"`, `"not_converged"`, `"linear_solver_failed"`), `method` (`"direct_kkt"` or `"solve"`), `distortion`, `objective`, `total_cone_angle` (sum of `|angle|`), `n_components`, `certificate_holds` (at tolerance `1e-9`), `state_residual`, `duality_gap`, `seconds` (the solver), `message` (`str` or `None`), `seconds_total` (everything, including decoding) |

Errors: `InputError` (message starts with the parameter name) for malformed buffers, wrong
lengths, an unknown `boundary` or `method`, a cone vertex out of range, repeated, on the boundary
or carrying a prescribed cone, invalid weights or bounds, and a closed component whose cones
cannot carry `2πχ` (no cones, or bounds that forbid the total angle). `MeshError` for a mesh
that is not a manifold, consistently oriented triangle mesh with positive-area triangles and
no unused vertices. A result that did not converge is not an exception: read `report["status"]`.

## `map_to_sphere`: a closed genus-0 mesh to the unit sphere

```python
ddg_tools.map_to_sphere(positions, triangles, *, pole=None, center=True) -> ddg_tools.SphereResult
```

A conformal map of a connected closed triangle mesh of genus 0 to the unit sphere, by Boundary
First Flattening (Sawhney, Crane 2017, Sec. 7.5): one vertex star is removed, the rest is
flattened to a round disk, projected stereographically and Möbius-centered (Baden, Crane,
Kazhdan 2018) so that the surface's area is balanced around the origin. The result is
determined up to a rotation of the sphere. The binding of `ddg_bff::map_to_sphere`; the crate
README of `ddg-bff` has the method, the deviations from the reference and measurements.

| argument | meaning |
|----------|---------|
| `positions` | 3 reals per vertex |
| `triangles` | 3 vertex indices per triangle: a connected, closed, consistently oriented manifold mesh of genus 0 in which every vertex is used |
| `pole` | the vertex that is removed (it lands on the sphere like every other vertex); `None` (default) picks the vertex with the most regular star |
| `center` | Möbius-center the result (default `True`); without it the pole is at `(0, 0, −1)` |

`SphereResult` (read-only):

| attribute | contents |
|-----------|----------|
| `positions` | `array('d')`, 3 per vertex: a point of the unit sphere (up to rounding; `max_radius_error`) |
| `report` | `dict`: `pole`, `is_injective` (no flipped or degenerate triangle and covering degree 1), `flipped_triangles`, `degenerate_triangles`, `covering_degree` (signed area of the image / `4π`), `quasi_conformal_mean`, `quasi_conformal_max` (`σ_max/σ_min` per triangle, 1 is conformal), `area_distortion` (area-weighted rms of the log area ratio about its mean, 0 when all areas scale alike), `log_scale_min`, `log_scale_max`, `max_radius_error`, `centering_status` (`"converged"`, `"max_iterations"`, `"stalled"`, `"singular_jacobian"`, or `None` without centering), `centering_iterations`, `center_of_mass` (`|µ|` at the end, or `None`), `seconds`, `seconds_total` |

Errors: `InputError` for malformed buffers and a `pole` that is not a vertex (message starts with
the parameter name); `MeshError` for a mesh that is not a valid manifold triangle mesh, has a
triangle of zero area, has a boundary, several components or a genus other than 0 (the message
says which; split a mesh into its components first); `SolveError` if the computation fails;
`TypeError` for an argument without the buffer protocol.

## The cut flow: `CutSession`

`ddg_tools.CutSession` binds `ddg_uv::cutting` (`crates/ddg-uv`, README section on the cut flow):
the cut flow of "Variational Surface Cutting" (Sharp, Crane 2018, `crates/ddg-variational-cuts`,
implemented clean-room from the text of the paper) on a **polygon mesh with seams**. A session
holds the state of the flow, advances it, and describes the cut on the caller's polygons. Same
conventions as `unwrap`: buffers in (any C-contiguous buffer, flat or `(n, width)`), `array.array`
out, reports as plain dictionaries, errors as exceptions carrying the Rust message. numpy is
never required (the state is returned as a `memoryview` that `numpy.asarray` turns into a 2-D
array).

```python
s = ddg_tools.CutSession(
    positions, corner_vertices, polygon_starts, triangles, *,
    n_regions=None, triangulation="intrinsic", refine_bound=None, optimal_delaunay=False, state=None, seam_edges=None, pinned_vertices=None, forbidden_vertices=None,
    prescribed_vertices=None, prescribed_angles=None,
    # exactly one start:
    circle_vertex=None, circle_radius=None,     # geodesic disk around a vertex (region 0)
    polygon_labels=None,                        # one label per polygon
    corner_values=None,                         # (n_regions, n_corners), as `state`
    cluster_iterations=None,                    # normal clusters, any n_regions (the atlas)
)
s.state                       # memoryview, float64, shape (n_regions, n_corners); assignable
s.step(**controls) -> dict    # one step
s.run(**controls, max_steps=1000, motion_tolerance=1e-5, energy_tolerance=1e-4,
      energy_window=25) -> dict
s.preview() -> (positions, starts, closed, regions) # the cut as polylines, 2 regions each
s.describe_cut() -> CutDescription                  # the cut on the polygons, as arrays
s.insert(cut=None) -> InsertedMesh                  # the Rust reference insertion
s.unwrap(cones=False, ...) -> CutUnwrapResult       # insert, then unwrap with the cut as seams
s.restart(...), s.roughening_rate(**controls), s.corner_regions(), s.state_row(k)
ddg_tools.visibility_penalty(positions, corner_vertices, polygon_starts, triangles, camera,
                             visibility=None)       # paper Sec. 7.1, no rays cast
ddg_tools.length_weight_for_ripples(roughening_rate, wavelength)
```

**Units.** The state and the positions are in the unit of the mesh. The weights, the penalty, the
step sizes and the ripple wavelength refer, as in the library, to the surface scaled to unit total
area; `s.length_scale` converts a length of the mesh to that unit (`s.total_area` is the area).
Reports give the cut length in both (`cut_length`, `cut_length_unit_area`).

**The state.** One function per region and per polygon **corner**, so that the two sides of a seam
can differ. With two regions the rows are `phi` and `-phi`; region 0 is where `phi <= 0`. Setting
the state accepts any float buffer of that shape (or flat, row-major), resets time and step count,
and refuses corners of one vertex that no seam separates but that carry different values (the
error names the vertex).

**The triangulation the flow runs on.** `triangulation="intrinsic"` (the default) runs the flow on
an intrinsic triangulation of the surface flipped to Delaunay (`ddg-intrinsic`, after Sharp,
Soliman, Crane, *Navigating Intrinsic Triangulations*, 2019): the flow's stages degrade on badly
shaped triangles (on a 165 degree sheared grid the speed is off by 110 % with the paper's stages
and the cut vanishes), flips alone restore regular-grid behaviour, and usually no vertex is added.
`"refined"` adds refinement to `refine_bound` degrees (default 25, at most 35) and adds
vertices; `"caller"` is the caller's triangulation, as the versions before this one.
`optimal_delaunay=True` adds the intrinsic optimal Delaunay triangulation (iODT, default options)
on top of `"intrinsic"` or `"refined"`: more even triangle areas, much slower (setup 4.5 s at 11.7k
vertices and 20 s at 46.9k against 0.15 s and 0.55 s without, flow steps about 2.5x); `ValueError` with `"caller"` or with `state=` (a
state carries its flag; read it back from `CutSession.optimal_delaunay`). The report gains
`n_optimized` and, with the flag, `optimal_delaunay_status`/`_iterations`/`_last_motion`. The cut
comes back onto the caller's polygons whatever the choice. `intrinsic_report` is a dictionary
(`n_added`, `n_initial_flips`, `n_splits`, `min_angle_degrees`, `max_angle_degrees`,
`is_delaunay`, `seconds`, `fingerprint`, ...; `None` for `"caller"`). **State.** `state` (the
`(n_regions, n_corners)` array) is the whole state exactly when `state_is_per_corner` is true
(no vertex added). Otherwise keep `snapshot()` (bytes: the levels on all vertices of the
triangulation, its fingerprint, the number of regions and the option): `CutSession(..., state=bytes)`
rebuilds the session on the same mesh (the triangulation and `n_regions` come from the bytes) and
`restore(bytes)` sets a live one; a state of another mesh or other seams raises
`StaleStateError`. A curvature per triangle is refused on an intrinsic triangulation
(`InputError`); use `"caller"` for it. `merge_tolerance` (default 0.1, 0 keeps every bend; also settable when a session is rebuilt from
state bytes) removes bends of the returned cut that deviate from the straight line between their
neighbours by less than that fraction of the longest edge of their triangle:
`CutDescription.report` has `n_merged_nodes`, `max_displacement` and `max_displacement_relative`
(and `n_moved_points`, `n_removed_excursions`). `CutSession.n_patches` counts the connected pieces of the surface in one
region; `InsertedMesh.changed_polygons` (with `changed_vertex_starts`, `changed_vertices`,
`changed_triangle_starts`, `changed_triangles`) is the library's triangulation of the polygons that
were split or got extra vertices, which an editor should keep. The description has several points per edge, junctions of
any degree (a bend is a junction with two paths), paths between junctions and paths with both
ends on one polygon edge; `CutDescription.lines()` / `CutSession.lines()` give the paths joined
for display; `InsertedMesh.polygon_region` is the region of each new polygon where known. The
Blender operators use all of this (`blender/README.md`).

**More than two regions.** `n_regions` can be any number from 2. The starts: `cluster_iterations`
(the vertex normals clustered into `n_regions` clusters, at most that many rounds: the start of
the atlas of the paper), `polygon_labels` and `corner_values` with N values, and the geodesic circle
for two regions only (more raise `InputError`). `stages` defaults to the recommended set: the
default stages for two regions, the paper's stages with the flux recovery of the core crate
(`FlowOptions::for_regions`, SPEC D20) for more; `stages="default"` and `"paper"` choose the
others, and every individual stage option overrides as before. `distortion_weights` and
`area_weights` have one entry per region. `step`/`run` report `distortion_k` and `area_k` for
every region. `preview()` returns the curves of the network between its junctions, each with the
two regions it separates. `describe_cut()` has the junctions (a triple point inside a triangle of
Blender's triangulation): `junction_polygon`, `junction_triangle`, `junction_corners` (3 per
junction), `junction_weights` (3 per junction; the junction is the weighted sum of the corner
positions), `junction_positions`; three paths end at each junction (`path_start`/`path_end` of
`-(j+1)`). `insert()` returns `junction_vertex` and marks junction vertices as kind 3 in
`vertex_source`. The construction is in the rustdoc of `ddg_uv::cutting::insert`.
`NotAvailableError` stays as a class for compatibility but is not raised for this.

**Controls** of `step` / `run` (keyword-only; `None` is the library's default; an unknown keyword
is a `TypeError`):

| keyword | meaning | default |
|---------|---------|---------|
| `energy` | `"dirichlet"`, `"scaled_dirichlet"`, `"hencky"` | `"hencky"` (the example program's) |
| `scale_weights` | per-vertex weight `a > 0` of the scaled Dirichlet energy | the curvature-based weight |
| `feature_r_eps`, `feature_smoothing` | cap `r_eps` and smoothing `c` of the curvature-based weight | `1.0` (example program), `0.02` (paper) |
| `length_weight` / `ripple_wavelength` | `alpha_L` explicit, or set from the roughening rate of the state so that ripples longer than the wavelength grow (`length_weight_for_ripples`); in `run` once at the start | `ripple_wavelength=0.25` (example program) |
| `penalty`, `penalty_rule` | `eta` per vertex; `"gradient"` or `"paper_diagonal"` | none, `"gradient"` |
| `distortion_weights`, `area_weights` | `w_k`, `beta_k`, one per region | all 1, all 0 |
| `symmetry` | class per vertex (equal classes are kept at equal values; a vertex split by a seam is left out) | none |
| `curvature_uniform` / `curvature_per_triangle` | curvature override (density per unit area of the unit-area surface) | Gaussian curvature |
| `stages` | `"default"`, `"paper"` (stages as printed) or `"recommended"` | recommended: `"default"` for 2 regions, the paper's stages with flux recovery for more |
| `flux_recovery`, `flux_smoothing`, `sliver_aspect_ratio`, `smoothing_passes`, `length_step`, `extension`, `speed_scaling`, `redistance` | the individual stages | those of `stages` |
| `step_rule`, `time_scale`, `max_motion`, `tau_max`, `tau` | `"adaptive"`, `"capped"` or `"fixed"` (`tau`) | adaptive, 0.05, 0.25, 1 |
| `min_crossing_offset` | `delta`, smallest distance of a crossing from a vertex | `1e-3` |
| `max_steps`, `motion_tolerance`, `energy_tolerance`, `energy_window` | the stop rule of `run`; a tolerance of 0 never stops for it | 1000, `1e-5`, `1e-4`, 25 |

A pointwise energy given as a Python callable (`Energy::Pointwise` of the library) is **not
offered**: the energies above cover the paper's examples, and a callable would need a call back
into Python at every vertex of every step.

**Reports.** `step` returns the energy terms of the state the step started from
(`energy_total`, `distortion_k`, `area_k`, `cut_length`, `penalty`), `tau`, speeds,
`roughening_rate`, `cut_motion`, counters, residuals, stage timings and the `length_weight` used.
`run` returns `steps_taken`, `stop_reason` (`"max_steps"`, `"converged"`, `"energy_stalled"`,
`"energy_rose"`, `"cut_vanished"`), `returned_step`, `lowest_energy`, `lowest_energy_step`,
`last_energy` (after `"energy_rose"` the session is back at the state of lowest energy, which
has `returned_step` steps behind it), `energy_before`, `energy_after`, `cut_length_before` / `cut_length_after`, the
final `distortion_k`, `area_k`, `length_weight`, `time`, and `steps`: the list of the `step` dicts.

**The cut on polygons** (`CutDescription`, all arrays flat): *edge points* (`edge_vertices` 2 per
point `[min, max]`, `edge_t` the parameter from `edge_vertices[0]`, `edge_positions`, `edge_side`:
on a seam the polygon whose side the point is on, else -1) and *paths*, one piece of the cut inside
one polygon (`path_polygon`, ascending; `path_start`, `path_end`: an edge point, or `-(j+1)` for a
junction; `path_regions` left and right; the diagonals a path crosses are `interior_corners`
(2 polygon corners), `interior_t`, `interior_positions`, delimited by `path_interior_starts`).
`insert` returns the reference result (`InsertedMesh`: positions, corners, polygon starts,
triangles, `new_seams`, `polygon_parent`, `vertex_source`, `edge_point_vertex`).

`unwrap(cones=False)` unwraps the inserted mesh with the old seams and the cut as seams and no
cone besides; `cones=True` lets the cone solver place cones as well. The UVs are not packed.

Errors are `InputError`, `MeshError` and `SolveError` as for `unwrap`, plus `NotAvailableError`
(a request the library does not implement yet, such as `n_regions > 2`).

**Keep the library's triangulation of the polygons you split.** `insert()`/`unwrap()` triangulate
the polygons the cut splits (`InsertedMesh.triangles`, `polygon_starts`, `vertex_source`). An editor
that inserts the cut itself and then lets its own triangulation run (Blender) gets a few flipped
slivers after the flattening (spot, eight clusters: 2 to 5 at the crossing offset `1e-3`); passing
the library's triangles for the split polygons removes them at every offset (the Blender add-on
does this, `blender/README.md`, "stored triangulation"). `min_crossing_offset` of `describe_cut`,
`preview`, `insert` and `unwrap` (default `1e-3`) only moves the cut away from vertices.

## Applications of intrinsic triangulations

Four functions over the crate `ddg-intrinsic-apps` (Sharp, Soliman, Crane, *Navigating Intrinsic
Triangulations*, 2019, Sec. 5; measurements and decisions in `crates/ddg-intrinsic-apps/README.md`
and `SPEC.md`). Each takes `positions` (3 floats per vertex) and `triangles` (3 indices per
triangle, a manifold oriented mesh), works on an intrinsic triangulation of that surface, and
returns per-vertex `array.array` results with a `report` dictionary. The triangulation is one of
`"input"` (the mesh), `"delaunay"` (its intrinsic Delaunay triangulation) and `"refined"`
(intrinsic Delaunay refinement, 25 degrees by default, as `CutSession`'s `"refined"`).

```python
geodesic_distance(positions, triangles, sources, *, triangulation="refined", time_factor=1.0,
                  boundary="neumann", representation="signposts")
                                                 # -> .distances, .report         (Sec. 5.3)
steiner_tree(positions, triangles, *, intrinsic=True, steiner_points=True,
             representation="signposts")
                                                 # -> .points, .offsets, .report  (Sec. 5.1)
smoothest_direction_field(positions, triangles, *, n=1, triangulation="delaunay",
                          laplacian="fem", min_angle_degrees=25.0)
                  # -> .vectors, .singularity_points, .singularity_indices, .report (Sec. 5.5)
adaptive_solve(positions, triangles, problem, *, source=0, time=None, max_vertices=20000,
               bulk=0.5, min_angle_degrees=25.0, representation="signposts")
                                                  # -> .values, .n_components, .report (Sec. 5.4)
```

`representation=` (`geodesic_distance`, `steiner_tree`, `adaptive_solve`) is `"signposts"` (the
default, Sharp, Soliman, Crane 2019) or `"integer"` (integer coordinates, Gillespie, Sharp, Crane
2021; the results agree with the signposts' to rounding, bitwise on the tested meshes, see
`crates/ddg-intrinsic-apps/README.md`). The default is signposts until the Thingi10k evaluation
(Integer Coordinates WP10) passes; one later commit then changes it (Q3).
`smoothest_direction_field` stays on signposts (it needs input-aligned tangent frames, contract
I26).

### `IntrinsicTriangulation` (integer coordinates)

```python
t = IntrinsicTriangulation(positions, triangles, representation="integer", mollify=None)
t.flip_to_delaunay()                                   # -> dict (status, n_flips, ...)
t.refine(min_angle_degrees=25.0, rules="chew", max_insertions=None, max_circumradius=None)
t.counts(); t.triangles(); t.vertex_places()
t.common_subdivision()   # positions, on_input, on_triangulation, polygon_starts, corners,
                         # input_face, intrinsic_face   (Fig. 12 of the integer paper)
t.edge_polylines(edges=None)     # intrinsic edges as polylines on the input mesh (for drawing)
t.certify()                      # SPEC 2.4; integer only
t.to_input(values); t.to_intrinsic(values)   # L2 transfer between T1 and the input; integer only
```

A session starts as the input mesh (`T1 = T0`, optionally mollified by the relative margin
`mollify`; the paper uses `1e-4`), is flipped and refined in place and never changes the input.
`rules="chew"` are the integer paper's refinement rules, `"ruppert"` the signpost paper's;
`representation="signposts"` runs the same session on signposts (no certificate, no transfer).
Vertex ids below `n_input_vertices` are the input's. The common subdivision lists the
triangulation's vertices first and then the crossings of intrinsic edges with input edges; every
polygon is convex with its input face and intrinsic face. Types: `ddg_tools.pyi`; tests:
`python-tests/test_integer_coords.py` (run by `scripts/py-test.sh`).

* `geodesic_distance`: the heat method, factored once; distance 0 at the source vertices.
* `steiner_tree`: a short cut through all vertices (the MST of the iDT shortened by Fermat
  points), as polylines on the surface (`offsets[k]` is where the polyline of tree edge `k`
  starts in `points`); `report["is_disk"]` says whether cutting along it leaves a disk.
* `smoothest_direction_field`: the smoothest `n`-direction field (1 vectors, 2 lines,
  4 crosses), one unit vector per vertex in the plane of an input face, and its singularities
  (index in units of `1/n`; on a closed surface they add up to `n` times the Euler
  characteristic, `report["index_sum"]`).
* `adaptive_solve`: `"green"` (harmonic Green's function of `source`), `"heat"` (one backward
  Euler step of length `time` from a unit point source; default the squared mean edge length)
  or `"disk"` (a harmonic map of a disk-shaped mesh onto the unit disk, 2 components per
  vertex), solved with adaptive refinement up to `max_vertices`; values at the mesh's vertices.

Errors: `InputError` (bad arguments, the message starts with the parameter), `MeshError` (not a
manifold triangle mesh; not a disk for `"disk"`), `SolveError`. Glue: `src/python/intrinsic_apps.rs`,
`src/python/intrinsic_vector.rs`; logic (pure Rust, tested by `cargo test`):
`src/intrinsic_apps.rs`, `src/intrinsic_vector.rs`; tests of the wheel:
`python-tests/test_intrinsic_apps.py`, `python-tests/test_intrinsic_vector.py`.

## Developability of triangle meshes

Four functions over the crate `ddg-developability` (Stein, Grinspun, Crane, *Developability of
Triangle Meshes*, ACM TOG 2018, implemented clean-room from the text; contract
`crates/ddg-developability/SPEC.md`, section 7.2). Input is the polygon mesh of `unwrap`
(`positions`, `corner_vertices`, `polygon_starts`, Blender's `triangles`; the triangulation of
the caller is used, the polygons are never changed). The mesh must be a manifold, consistently
oriented triangle mesh; vertices no polygon uses are allowed. Errors name polygons and vertices of
your mesh (`MeshError`: "the mesh is not manifold: edge (0, 2) is shared by 3 faces: polygons 0, 4,
6"), `InputError` starts with the parameter name. Results are objects with attributes
(`array.array` and a `report` dict), as for `map_to_sphere`.

```python
developability(positions, corner_vertices, polygon_starts, triangles, *, energy="covariance",
               include_boundary=True, omit_valence_three=True, threads=None)
        # -> .vertex_energy ('d'), .counted ('B'), .energy, .report
make_developable(..., energy=, include_boundary=, omit_valence_three=, method="lbfgs", memory=10,
                 step=1e-3, line_search="wolfe", first_trial="unit", first_trial_growth=2.0,
                 first_trial_max=100.0, iterations=500, subdivisions=2,
                 coarse_energy="combinatorial", max_displacement=0.5, gradient_tolerance=None, energy_tolerance=None,
                 stall_iterations=None, initial_step=None, min_face_angle=None,
                 pinned_vertices=None, remesh=True,
                 small_angle=radians(1), collapse_boundary_edges=True, fallback_operations=True,
                 return_faces=False, threads=None, progress=None, progress_interval=0.1)
        # -> .positions ('d'), .faces ('i' | None), .face_polygon ('i' | None),
        #    .vertex_source ('i'), .report
developable_seams(..., seam_edges=None, extra_terminals=None, energy=None, epsilon=1e-3,
                  defect_tolerance=0.1, length_weight=1.0, alignment_weight=0.0, threads=None)
        # -> .seam_edges ('i', 2 per edge), .terminals, .vertex_energy, .vertex_defect,
        #    .star_embeddedness ('d'), .report
developable_unwrap(..., <the same options>, extension="harmonic",
                   overlap_samples_per_triangle=None, margin=1/256)
        # -> .corner_uvs ('f', packed), .seam_edges, .polygon_chart, .corner_log_scale, .terminals,
        #    .vertex_energy, .vertex_defect, .star_embeddedness, .report
ENERGIES, DEFAULT_DEFECT_TOLERANCE                # module constants
```

* **Energies** (`ENERGIES`): `covariance` (E_lambda, the default), `combinatorial` (E_P),
  `combinatorial_max`, `covariance_max`, `intrinsic`, `intrinsic_max`. `omit_valence_three`
  leaves out interior vertices with three faces; `include_boundary` counts boundary vertices. The
  results do not depend on `threads`.
* **`make_developable`** runs the flow in one library call (`subdivisions + 1` rounds of `iterations`
  L-BFGS iterations, a 4-1 subdivision before every round but the first, small-angle remeshing
  during the rounds). `coarse_energy` (default `"combinatorial"`, D36) is the energy of round 0
  when `subdivisions > 0`, as in the paper (Sec. 5.2: `E_P` on the coarse mesh, then `energy` after
  each subdivision); `None` uses `energy` in every round (the flow of earlier versions, bitwise).
  Start from a coarse, uniform mesh (about 1k triangles); `E_λ` alone, or a subdivision cage,
  flattens thin parts (`crates/ddg-developability/QUALITY.md`). `fallback_operations` (D35) lets
  the remeshing try the other flips/collapses of a triangle when the paper's is refused. `method` is `lbfgs`, `gradient_descent` or `fixed_step` (`step` in units
  of the squared mean edge length); `line_search` is `wolfe` or `backtracking` (c1 = 1e-4, shrink
  1/2); `first_trial` is the first trial step of the `wolfe` search (D37): `unit` (1 at every
  iteration) or `grow` (`min(first_trial_growth` × the previous accepted step`, first_trial_max)`,
  defaults 2 and 100), lower final energy in 28 of 36 study runs but it can get stuck (spot,
  `intrinsic`); `max_displacement` caps a step in mean edge lengths; `min_face_angle`
  is the admissibility floor of the library on the smallest corner angle (D31; `None` keeps its value, 0 switches the
  check off). `pinned_vertices` do not move. **D24:** `faces` is `None` exactly when no
  subdivision ran and no remeshing operation was applied; then `positions` has one entry per
  vertex of your mesh and your polygons stay. The report then holds the largest
  non-planarity of polygons with more than three corners (`nonplanarity_before/after`,
  `nonplanarity_after_relative` = relative to the mean side of that polygon,
  `nonplanarity_polygon`). Otherwise `positions`/`faces` are the new triangle mesh (vertices no
  polygon used are dropped), `face_polygon` gives for every triangle the polygon it descends from
  (exact, also after several rounds with flips and collapses: the library tracks the parent
  triangle of every output triangle, `FlowResult.face_parent`) and `vertex_source`
  the vertex of your mesh each vertex descends from (-1 for a subdivision midpoint).
  `return_faces=True` returns `faces`/`face_polygon` also when nothing changed (the input's
  triangles). `progress(round, iteration, energy, gradient_norm)` is called from the calling
  thread at most every `progress_interval` seconds (0: every iteration) with the GIL; returning
  `False` stops the flow and the state at that moment is returned (`report["cancelled"]`); an
  exception, and Ctrl-C (checked at the same moments), stop the flow and propagate. Report:
  `energy_initial/final`, `topology_changed`, `n_*`, `iterations`, `min_angle`, `rounds` (a list of
  dicts: `energy` (the energy of that round), `status` (also `stationary`, D34), `iterations`,
  `flips`, `collapses`, `fallback_flips`, `fallback_collapses`, `remesh_rejected`,
  `inadmissible_trials`, ...), `bad_faces` (triangles below `small_angle` that nothing could repair),
  `inadmissible_trials` (steps refused by `min_face_angle`), `stopped` (`None`, or a dict
  `round`, `after_subdivision`, `vertex`, `reason`: the flow ended before that round because the
  energy is undefined at its starting state; the result is the state before it; also a warning),
  `warnings`, timings. `min_face_angle` also sets the floor of the remeshing's collapses and
  must be below `small_angle`; an input with a triangle below it is a `MeshError` (naming the
  polygon) unless remeshing is on (then it is repaired first, and the message names a face of
  the current mesh).
* **`developable_seams`** is the cut of the paper (Sec. 4.4) without the flattening: terminals are
  the vertices whose energy exceeds `epsilon` (evaluated with boundary vertices included and
  valence 3 counted) and, since D30, interior vertices whose angle defect exceeds
  `defect_tolerance` (`None`: the paper's rule), plus `extra_terminals`; a cut graph through
  them (Steiner tree approximation, handle loops, `length_weight`/`alignment_weight` in the edge
  cost) along **polygon edges only**; `seam_edges` are cut as well. It returns the edges as vertex
  pairs `(min, max)` and the per-vertex energy, angle defect and star embeddedness (the library's
  diagnostic of D28/D30) that explain the choice. Feed the edges to `unwrap(seam_edges=...)`: the
  atlas pipeline then makes the UV map with its cones and packing (this is what the Blender
  operator does).
* **`developable_unwrap`** is the paper's own flattening: each chart by Boundary First Flattening
  with the boundary lengths kept, packed by `pack_charts`. It is isometric where the surface is
  developable, but a chart can overlap itself (spot: one chart, about half of it covered twice),
  so the report counts flips and overlap per chart (`charts`). `corner_log_scale` is
  `0.5 ln(UV area / area)` of the triangle that owns the corner.

Glue: `src/python/developability.rs`; logic (pure Rust, tested by `cargo test`):
`src/developability.rs`; tests of the wheel: `python-tests/test_developability.py`.

## Repulsive curves: `RepulsiveSession`, `repulsive_flow`, `curve_energy`

The flow of Yu, Schumacher, Crane, *Repulsive Curves*, ACM TOG 40(2), 2021, over the crate
`ddg-repulsive-curves` (contract `crates/ddg-repulsive-curves/SPEC.md`, section 9.2): the
tangent-point energy of a curve network (closed curves, open curves, junctions) is decreased by
fractional Sobolev descent, with constraints (edge lengths, pins, tangents, surfaces),
penalties, obstacles that are never crossed, and optional growth, shrinking and subdivision.
`positions` is 3 reals per vertex, `edges` 2 vertex indices per edge; both take any buffer
(`array.array`, numpy) and, for convenience, plain Python sequences (nested or flat). Results are
`array.array` objects (`numpy.frombuffer` makes a view) and dicts.

```python
import ddg_tools
s = ddg_tools.RepulsiveSession(positions, edges, length="edges", pinned=[0], max_subdivisions=1)
report = s.step()                      # one step, GIL released; a dict
s.run(max_steps=50, progress=lambda r: True)   # until final / 50 steps / progress -> False
pos, edges, iteration, energy, status = s.snapshot()   # never waits for a running step
s.cancel()                             # from any thread; reset() starts again
positions, edges, report = ddg_tools.repulsive_flow(positions, edges, steps=100)
energy, differential = ddg_tools.curve_energy(positions, edges, alpha=3, beta=6)
```

* **Live use (Blender).** A worker thread calls `step()` (or `run`), a timer reads `snapshot()`
  or the properties `positions`, `edges`, `energy`, `iteration`, `status`; they come from the last
  completed step and never wait. `cancel()` ends a running step between line-search trials; the
  curve is that of the last completed step. After `cancel()` further `step()` calls compute nothing
  and return `status == "cancelled"`; `reset()` returns to the input curve. A `RepulsiveSession`
  stepped *n* times gives bitwise the same curve as `repulsive_flow(..., steps=n)`, whatever
  `threads` is.
* **Status.** `"running"`, `"converged"`, `"max_steps"`, `"time_limit"`, `"stuck"` (the line search
  found no step), `"cancelled"`, `"error"` (then `error` holds the message). A final session does
  not step in `run()`, but `step()` always tries.
* **Subdivision** (`max_subdivisions=k`): when the mean edge length exceeds twice its start value,
  every edge is split. The vertex count changes (`n_vertices`, report key `subdivided`). Old
  vertices keep their indices, the midpoint of edge `I` is vertex `|V| + I`, and edge `I` becomes
  edges `I` and `|E| + I`; `vertex_source` (input vertex or `-1`) and `edge_source` (input edge)
  follow that through any number of subdivisions.
* **Options** (keywords; `repulsive_defaults()` lists them all with the library's defaults, and
  `ddg_tools.pyi` has them typed; an unknown name is an error): `alpha=3`, `beta=6`;
  `length="edges"|"total"|"free"`, `grow_factor=1` and `grow_rate=0.01` (the kept length moves to
  `factor` times itself, `rate` of the start length per step; `session.set_schedule(factor, rate)`
  changes it on the fly); `pinned`, `tangent_vertices`, `tangent_edges` (index lists), `surface`
  (a dict, `{"type": "sphere", "radius": 1.0}`, also torus, plane, double_torus, rounded_box,
  metaballs, union, smooth_union; key `"vertices"` restricts it), `barycenter` (default: fixed
  unless a vertex is pinned), `planar_axis="x"|"y"|"z"`; `obstacles=[(vertices, triangles)]`
  (never crossed with the default line search; `obstacle_weight > 0` adds a repulsive potential);
  penalties `total_length_weight`, `length_difference_weight`, `field` (`{"type": "constant",
  "vector": (0, 0, 1)}`, `"circular"`, `"sphere_tangent"`) with `field_weight`;
  `line_search="collision"|"armijo"` (`"armijo"` does not check collisions), `line_search_start`,
  `armijo`, `shrink`, `max_halvings`, `gap`; `backend="auto"|"dense"|"barnes_hut"|"multigrid"|
  "accelerated"` (Barnes-Hut energy and/or the multigrid solver; `"auto"` is the library's choice,
  at present exact and dense); `threads`; `projection_tolerance` (a fraction of the start mean
  edge length, or an absolute length with `projection_absolute=True`),
  `projection_iterations`; `stop="gradient"|"relative"|"stall"` with `stop_threshold` (default
  1e-4 / 1e-3; required for `"stall"`) and `stall_window`; `max_steps=1000`, `time_limit`
  (seconds, for `run`).
* **Soft methods** (`method="lbfgs"` with `lbfgs_memory=30`, or `"ncg"`; default
  `"projected_gradient"`): the paper's soft problem. The barycenter and length constraints are
  dropped and a Hencky length penalty of weight `hencky_weight=10` at the start lengths is added.
  `pinned`, `tangent_vertices`, `tangent_edges`, `surface` and `grow_factor != 1` need the hard
  constraints and raise `RepulsiveInputError` naming the option.
* **Library options beyond the paper** (defaults are the library's): `near_pair=True` (accurate
  quadrature for close edge pairs; `False` is the paper's energy) with `near_pair_blend_start`,
  `near_pair_blend_end`, `near_pair_subpanels`, `second_neighbours`, `second_neighbour_rule=
  "fold_gated"|"one_sided"|"two_sided"|"gated"|"one_sided_gated"|"one_sided_gated_either"` (the
  last three need `second_neighbour_width`), `second_neighbour_contact_band`,
  `second_neighbour_cos_shut`, `second_neighbour_cos_open`, `second_neighbour_band`;
  `back_projection`, `second_order_correction`, `newton_polish` (`None`: chord method),
  `projection_early_abort`, `stuck_projected_reference`, `stuck_retry_chord`; `mg_max_threads`
  (`None`: no cap on the multigrid pool). `time_limit` is honoured by `run` only, not by `step`.
* **Step report** keys: `iteration`, `status` (`"accepted"`, `"stuck"`, `"converged"`,
  `"retarget_stalled"`, `"cancelled"`), `session_status`, `n_vertices`, `n_edges`,
  `energy_before`, `energy_after`, `tangent_point_after`, `penalties_after`, `gradient_norm`,
  `sobolev_norm_sq`, `tau`, `tau_max`, `halvings`, `trials`, `max_violation`,
  `projection_iterations`, `linear` (solver, iterations, residual, levels), `retarget`,
  `subdivided`, `energy_backend`, `solver_backend`, `seconds` (per phase).
* **Errors.** Anything wrong in the input or the options (a bad edge, a pin out of range, an
  infeasible constraint, an unknown option, a degenerate obstacle triangle) raises
  `ddg_tools.RepulsiveInputError`, which is both a `ValueError` and a `ddg_tools.InputError`, with
  the library's message; a failed solve raises `SolveError`. Nothing panics.

Glue: `src/python/repulsive.rs`; logic (pure Rust, tested by `cargo test`, including that the
defaults equal the library's): `src/repulsive.rs`; tests of the wheel:
`python-tests/test_repulsive.py` (which also checks the `.pyi` defaults against
`repulsive_defaults()`).

## What `unwrap` does

`src/pipeline.rs`, function `run_pipeline`, is the only call into a pipeline: it hands the
validated input to `ddg_uv::unwrap` and converts the result and the errors. Besides the
steering and `lambda_`, it sets `cone_signs` from `positive_cones_only` and switches packing
off (`packing: None`: the contract of `unwrap` is "not packed", `pack_charts` packs); every
other policy is `ddg-uv`'s default. In short (details and the measurements behind every default:
`crates/ddg-uv/README.md`):

1. Validation and repair as listed under *Inputs*.
2. Seams are cut first; the connected pieces are the charts.
3. Per piece: the spike check (see the report keys), optimal cones (`ddg-optimal-cones`), then the piece is cut through its cones along
   polygon edges and flattened conformally (`ddg-bff`). A closed piece of genus 0 uses the raw
   optimum unless it is dense (more than 32 cones and more than 1 % of the vertices); then it
   is solved again with 8 free cone sites. With `cone_feedback=True` the piece is then solved
   and flattened again with lowered distortion weights where the solver's scale factor and the
   flattening's disagree, until nothing changes.
4. Final checks: every UV is finite, all triangle corners of one polygon corner got the same
   UV, every cut edge is a polygon edge, every seam is cut.

The driver `pipeline::unwrap` validates the options, runs the pipeline, checks its output
against the input (lengths, index ranges, finite numbers), measures the distortion and
assembles the report. New options become fields of `UnwrapOptions` and keyword arguments in
`src/python.rs`.

## Layout

| file | contents |
|------|----------|
| `src/buffers.rs` | decoding buffer-protocol memory into `Vec<f64>` / `Vec<u32>` |
| `src/contract.rs` | `MeshBuffers`, `UnwrapInput` (validation), `UnwrapOptions`, `UnwrapOutput`, `Report` |
| `src/error.rs` | `InputError`, `UnwrapError` |
| `src/pipeline.rs` | `run_pipeline` (the call into `ddg-uv`) and the driver `unwrap` |
| `src/cones.rs` | `reoptimize_cones`: validation and the call into `ddg-optimal-cones` (pure Rust) |
| `src/sphere.rs` | `map_to_sphere`: validation and the call into `ddg-bff` (pure Rust) |
| `src/developability.rs` | `developability`, `make_developable`, `developable_unwrap` over `ddg-developability` (pure Rust) |
| `src/cutting.rs` | the cut session: controls, reports, arrays (pure Rust, tested by `cargo test`) |
| `src/measure.rs` | distortion measured on the result |
| `src/pack.rs` | packing charts into the unit square |
| `src/python.rs`, `src/python/cut.rs`, `src/python/sphere.rs`, `src/python/developability.rs` | the PyO3 glue (cargo feature `python`) |
| `ddg_tools.pyi` | type stubs, shipped in the wheel |
| `python-tests/` | tests of the built wheel |

### The crate inside the workspace

PyO3 is an optional dependency behind the cargo feature `python`, which is **off by default**.
A plain `cargo build`, `cargo test` or `cargo clippy`, also at the workspace root, therefore
compiles only the pure-Rust part of this crate: no PyO3, no Python interpreter, no libpython to
link against. `maturin` turns the feature on through `pyproject.toml`
(`[tool.maturin] features = ["python"]`). All logic lives in the pure-Rust part; `src/python.rs`
only converts.

The wheel contains MIT-licensed code only: ddg-core has no triangulator (Blender triangulates).

## Build and test

From the repository root (Python tooling goes through `uv`):

```sh
export CARGO_TARGET_DIR=$PWD/target-py

cargo test -p ddg-py                                              # pure-Rust tests
cargo clippy -p ddg-py --all-targets -- -D warnings
cargo clippy -p ddg-py --all-targets --features python -- -D warnings   # includes the PyO3 glue

scripts/py-build-wheel.sh       # dist/wheels/ddg_tools-<version>-cp311-abi3-<platform>.whl
scripts/py-test.sh              # builds the wheel, tests it under Python 3.11+numpy, 3.12, 3.14
```

`scripts/py-test.sh` installs the wheel into throw-away environments with `uv run`; the test
files `test_ddg_tools.py` and `test_cut_session.py` check that the module it imported is the wheel that was just built.

The Blender side is built and tested with `scripts/blender-test.sh`; see `blender/README.md`.

## References

This crate contains no algorithm of its own. The pipeline it calls (`ddg-uv`) implements

- Yousuf Soliman, Dejan Slepčev, Keenan Crane. *Optimal Cone Singularities for Conformal
  Flattening.* ACM Transactions on Graphics 37(4), 2018.
  <https://doi.org/10.1145/3197517.3201367>
- Rohan Sawhney, Keenan Crane. *Boundary First Flattening.* ACM Transactions on Graphics 37(1),
  2017. <https://doi.org/10.1145/3132705>
- Oded Stein, Eitan Grinspun, Keenan Crane. *Developability of Triangle Meshes.* ACM Transactions
  on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201303> (the developability functions,
  through `ddg-developability`)

## License

MIT.
