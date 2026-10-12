# ddg-py

The Python module `ddg_tools`: a mesh goes in as flat buffers, a UV per polygon corner comes out.
It is the bridge between Blender (`blender/` in this repository) and the Rust crates of the
workspace, built with [PyO3](https://pyo3.rs) for the stable ABI, so one binary loads in every
Python from 3.11 on. It never needs numpy.

> `ddg_tools.unwrap` runs the real pipeline, `crates/ddg-uv`: marked seams are cut, optimal cone
> singularities are computed per piece (`crates/ddg-optimal-cones`), the surface is cut through
> the cones and flattened conformally (`crates/ddg-bff`). `ddg_tools.PIPELINE` is `"ddg-uv"` and
> `ddg_tools.IS_PLACEHOLDER` is `False`. Sections below that still describe a stand-in planar
> projection are from before the pipeline was wired in and need rewriting; the signature, buffer
> formats and error types they describe are current.

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
  Passing Blender's own triangles is how the Rust side triangulates quads and n-gons exactly
  like Blender.
- Positions and angles are finite. Indices are in range.
- `seam_edges`, `pinned_vertices` and `forbidden_vertices` are sets: order, orientation and
  duplicates do not matter. A vertex may be in at most one of the pinned, prescribed and
  forbidden sets, and at most once in `prescribed_vertices`.
- `lambda_` is finite and not negative. It is the trade-off `λ` of the paper between area
  distortion and total cone angle (dimensionless, default 0.5, useful range 0.1 to 1; larger
  means fewer cones).

What is not checked at this level, because it is the pipeline's business: manifoldness,
orientation, whether a seam edge is an edge of the mesh, degenerate triangles.

### Outputs

`UnwrapResult` has read-only attributes. All arrays are `array.array`; they work with Blender's
`foreach_set`, with `memoryview` and with `numpy.frombuffer` / `numpy.asarray`.

| attribute | type | contents |
|-----------|------|----------|
| `uvs` | `array('f')`, 2 per polygon corner | the UV of corner `c` is `uvs[2c], uvs[2c+1]`, in corner (loop) order. Computed in float64, rounded to float32. **Not packed**: every chart has its own origin and the scale the flattening gave it, so charts overlap until `pack_charts` lays them out. |
| `polygon_charts` | `array('i')`, 1 per polygon | the chart of the polygon. A chart is a piece of surface flattened as one. |
| `cone_vertices` | `array('i')` | vertices that carry a cone, ascending |
| `cone_angles` | `array('d')` | their cone angles in radians: `2π − (angle sum in the flat metric)`, positive for a cone tip |
| `cut_edges` | `array('i')`, 2 per edge | the edges the surface was cut along, as vertex pairs `min, max`, sorted. Includes the seams that were honoured. Cuts follow polygon edges only, never a diagonal of the triangulation. |
| `report` | `dict` | plain values, see below |

### The report

Values are `bool`, `int`, `float` or `str`; `warnings` is a `list[str]`. The dict is ordered.

| key | meaning |
|-----|---------|
| `pipeline`, `placeholder` | which pipeline ran, and whether it is the stand-in |
| `n_vertices`, `n_polygons`, `n_corners`, `n_triangles` | size of the input |
| `n_seam_edges`, `n_pinned`, `n_prescribed`, `n_forbidden` | the steering as the Rust side received it (after removing duplicates) |
| `lambda`, `positive_cones_only` | the options as received |
| `n_charts`, `n_cones`, `n_cut_edges` | size of the result |
| *(pipeline keys)* | whatever the pipeline reports. The stand-in adds `n_unreferenced_vertices`, `n_seam_edges_ignored`, `n_fold_edges`. |
| `n_degenerate_triangles` | triangles with no area in 3D (left out of the measurements) |
| `n_collapsed_uv_triangles` | triangles with area in 3D but none in UV |
| `n_flipped_uv_triangles` | triangles whose UV image is mirrored |
| `measured_triangles` | triangles the following figures are computed over |
| `measured_angle_error_max`, `measured_angle_error_mean` | corner angle in UV minus corner angle in 3D, absolute, radians; the mean is weighted by 3D area. Zero for a conformal map. |
| `measured_log_scale_std`, `measured_log_scale_min`, `measured_log_scale_max` | per triangle, `u = ½ ln(area in UV / area in 3D)`, the log of the local length scale; reported relative to its area-weighted mean (standard deviation, smallest, largest). Zero for an isometry. `exp()` of min and max is the range of local scale relative to the average. |
| `seconds_decode`, `seconds_validate`, `seconds_pipeline`, `seconds_measure`, `seconds_encode`, `seconds_total` | wall-clock time: copying the buffers, validating them, the pipeline, the measurements, building the result, and all of it |
| `warnings` | things the caller should know that did not stop the unwrap |

The `measured_*` figures and the three triangle counts are computed by the binding from the
final UVs, independently of the pipeline. They are a check on the whole chain, including the
corner order.

### Errors

| exception | raised when |
|-----------|-------------|
| `ddg_tools.InputError` | an argument is malformed. The message starts with the parameter name: `corner_vertices[3] = 8 is out of range: there are 8 vertices (valid indices 0..=7)` |
| `ddg_tools.MeshError` | the arrays are fine but the mesh cannot be unwrapped. The message uses the caller's vertex and polygon ("face") indices: `edge (12, 57) is shared by 3 faces: 4, 9, 31` |
| `ddg_tools.SolveError` | the mesh was accepted but no valid map came out (reserved for the real pipeline; also raised if a pipeline returns malformed data) |
| `TypeError` | an array argument has no buffer protocol, or an option has the wrong type |

The three module exceptions derive from `ddg_tools.DdgError`.

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

## What the stand-in does

`src/placeholder.rs`, one function, `planar_projection`:

1. Builds the triangle mesh from Blender's triangles, drops vertices that no face uses (counted
   in `n_unreferenced_vertices`), and validates it with `ddg_core::TriMesh::new`. A non-manifold
   or inconsistently oriented mesh is rejected with ddg-core's message, rewritten in the caller's
   vertex and polygon indices.
2. Groups polygons into charts: joined across every shared edge that is not a seam. Seam edges
   that are not on any face are ignored with a warning (`n_seam_edges_ignored`).
3. Projects each chart onto the plane perpendicular to its area-weighted normal. For a chart
   whose normals cancel (a closed surface) it uses the coordinate plane perpendicular to the
   thinnest side of the bounding box.
4. Returns as "cones" the steering echoed back (pinned vertices with angle 0, prescribed cones
   with their angle) and as "cut edges" the honoured seams plus the edges where the projection
   folds over (`n_fold_edges`).

`lambda_` and `positive_cones_only` have no effect.

## Wiring in the real pipeline

`src/pipeline.rs` contains the only call into a pipeline:

```rust
pub fn run_pipeline(input: &UnwrapInput, options: &UnwrapOptions)
    -> Result<UnwrapOutput, UnwrapError>
{
    // STAND-IN. Replace this line by the call into `ddg-uv`.
    crate::placeholder::planar_projection(input, options)
}
```

To switch:

1. Add `ddg-uv` to `[dependencies]` in `Cargo.toml` (mind its features: the wheel is MIT only as
   long as ddg-core's `blender-tessellation` feature stays off).
2. Replace the body of `run_pipeline`. The input offers `positions() -> &[[f64; 3]]`,
   `corner_vertices() -> &[u32]`, `polygon_offsets() -> &[u32]` (`n_polygons + 1` entries),
   `triangles() -> &[[u32; 3]]` (corner indices), `seam_edges() -> &[[u32; 2]]`,
   `pinned_vertices()`, `prescribed_cones() -> &[Cone]`, `forbidden_vertices()`; all validated
   and canonical as described in `src/contract.rs`. The output to build is
   `UnwrapOutput { uvs: Vec<[f64; 2]>, polygon_charts: Vec<u32>, cones: Vec<Cone>,
   cut_edges: Vec<[u32; 2]>, report: Report }`. Map the pipeline's error to
   `UnwrapError::Mesh(message)` or `UnwrapError::Solve(message)`.
3. Set `PIPELINE_NAME` and `IS_PLACEHOLDER = false` in `src/pipeline.rs`.
4. Delete `src/placeholder.rs` and its `mod` line in `src/lib.rs`.
5. In `python-tests/test_ddg_tools.py` and `blender/tests/test_addon.py`, the expectations that
   only hold for the stand-in are guarded by `ddg_tools.IS_PLACEHOLDER`; review them.

Everything else stays as it is: the driver `pipeline::unwrap` validates the options, runs the
pipeline, checks its output against the input (lengths, index ranges, finite numbers), measures
the distortion and assembles the report. New options become fields of `UnwrapOptions` and
keyword arguments in `src/python.rs`.

## Layout

| file | contents |
|------|----------|
| `src/buffers.rs` | decoding buffer-protocol memory into `Vec<f64>` / `Vec<u32>` |
| `src/contract.rs` | `MeshBuffers`, `UnwrapInput` (validation), `UnwrapOptions`, `UnwrapOutput`, `Report` |
| `src/error.rs` | `InputError`, `UnwrapError` |
| `src/pipeline.rs` | the seam `run_pipeline` and the driver `unwrap` |
| `src/placeholder.rs` | the stand-in pipeline |
| `src/measure.rs` | distortion measured on the result |
| `src/pack.rs` | packing charts into the unit square |
| `src/python.rs` | the PyO3 glue (cargo feature `python`) |
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
file checks that the module it imported is the wheel that was just built.

The Blender side is built and tested with `scripts/blender-test.sh`; see `blender/README.md`.

## References

This crate contains no algorithm of its own. The pipeline it is meant for implements

- Yousuf Soliman, Dejan Slepčev, Keenan Crane. *Optimal Cone Singularities for Conformal
  Flattening.* ACM Transactions on Graphics 37(4), 2018.
  <https://doi.org/10.1145/3197517.3201367>
- Rohan Sawhney, Keenan Crane. *Boundary First Flattening.* ACM Transactions on Graphics 37(1),
  2017. <https://doi.org/10.1145/3132705>

## License

MIT.
