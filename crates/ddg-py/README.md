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
| `lambda`, `positive_cones_only`, `optimize_boundary` | the options as received |
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
| `seconds_stage_validate`, `seconds_stage_seams`, `seconds_stage_cones`, `seconds_cone_solver`, `seconds_spike_check`, `seconds_stage_flatten`, `seconds_stage_pack`, `seconds_stage_verify`, `seconds_unwrap` | the pipeline's wall-clock times per stage (`seconds_cone_solver` and `seconds_spike_check` are parts of `seconds_stage_cones`) |
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
    n_regions=2, seam_edges=None, pinned_vertices=None, forbidden_vertices=None,
    prescribed_vertices=None, prescribed_angles=None,
    # exactly one start:
    circle_vertex=None, circle_radius=None,     # geodesic disk around a vertex (region 0)
    polygon_labels=None,                        # one label per polygon
    corner_values=None,                         # (n_regions, n_corners), as `state`
)
s.state                       # memoryview, float64, shape (n_regions, n_corners); assignable
s.step(**controls) -> dict    # one step
s.run(**controls, max_steps=1000, motion_tolerance=1e-5, energy_tolerance=1e-4,
      energy_window=25) -> dict
s.preview() -> (positions, starts, closed)          # the cut as polylines
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
error names the vertex). `n_regions > 2` raises `ddg_tools.NotAvailableError` ("... N regions are
not available yet"); every array is already shaped for N regions.

**Controls** of `step` / `run` (keyword-only; `None` is the library's default; an unknown keyword
is a `TypeError`):

| keyword | meaning | default |
|---------|---------|---------|
| `energy` | `"dirichlet"`, `"scaled_dirichlet"`, `"hencky"` | `"hencky"` (the example program's) |
| `scale_weights` | per-vertex weight `a > 0` of the scaled Dirichlet energy | the curvature-based weight |
| `feature_r_eps`, `feature_smoothing` | cap `r_eps` and smoothing `c` of the curvature-based weight | `1.0` (example program), `0.02` (paper) |
| `length_weight` / `ripple_wavelength` | `alpha_L` explicit, or set from the roughening rate of the state so that ripples longer than the wavelength grow (`length_weight_for_ripples`); in `run` once at the start | `ripple_wavelength=0.25` (example program) |
| `penalty`, `penalty_rule` | `eta` per vertex; `"gradient"` or `"paper_diagonal"` | none, `"gradient"` |
| `distortion_weights`, `area_weights` | `w_k`, `beta_k` per region | all 1, all 0 |
| `symmetry` | class per vertex (equal classes are kept at equal values; a vertex split by a seam is left out) | none |
| `curvature_uniform` / `curvature_per_triangle` | curvature override (density per unit area of the unit-area surface) | Gaussian curvature |
| `stages` | `"default"` or `"paper"` (stages as printed) | `"default"` |
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
`"cut_vanished"`), `energy_before`, `energy_after`, `cut_length_before` / `cut_length_after`, the
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
   is solved again with 8 free cone sites.
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
| `src/cutting.rs` | the cut session: controls, reports, arrays (pure Rust, tested by `cargo test`) |
| `src/measure.rs` | distortion measured on the result |
| `src/pack.rs` | packing charts into the unit square |
| `src/python.rs`, `src/python/cut.rs` | the PyO3 glue (cargo feature `python`) |
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

## License

MIT.
