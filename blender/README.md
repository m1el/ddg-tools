# DDG Conformal UV Unwrap (Blender extension)

One operator takes the active mesh and generates a UV map for it; three more cut a mesh open
where it flattens best (*The variational cut*, below). The work is done in Rust by the
native module `ddg_tools` (`crates/ddg-py`), bundled with the extension as a Python wheel.

> The native module runs the real pipeline (`crates/ddg-uv`): marked seams are cut, optimal cone
> singularities are computed, the surface is cut through the cones and flattened conformally.
> Defaults: lambda 0.1 and positive cones only. What the pipeline does with a mesh, and why
> those defaults, is in `crates/ddg-uv/README.md`; the Python interface is in
> `crates/ddg-py/README.md`.

## The operator

**Conformal Unwrap (Optimal Cones)**, `bpy.ops.uv.ddg_conformal_unwrap`

- 3D viewport, Edit mode: **UV** menu. Object mode: **Object** menu. Also in operator search.
- Works on the active object if it is a mesh, from Object mode and from Edit mode.
- Unwraps the **whole mesh** (not just the selection, hidden faces included), in object space.
  As with Blender's own unwrap, the object's scale is not applied.
- Writes the **active UV map**, or creates one named "UVMap" if the mesh has none. Other UV maps
  are left alone.
- Does not change vertices, edges or faces, their order, or the selection. The only other thing
  it can change is seam flags, and only when asked to.
- Supports undo, and has a panel for **Adjust Last Operation**.

| property | Python name | default | meaning |
|----------|-------------|---------|---------|
| Lambda | `lambda_` | 0.1 | trade-off between area distortion and cones; larger gives fewer cones. The default is the largest value at which the cone set of the clean test meshes stops growing; the soft range of the slider is 0.01 to 1 |
| Positive Cones Only | `positive_cones_only` | on | allow only cone tips, no saddle cones. Negative cones lower the distortion but make the map overlap itself. A closed surface of genus 2 or more cannot be flattened with positive cones only: it gets cones of either sign and a warning |
| Optimize Boundary | `optimize_boundary` | off | let the solver choose how much to stretch the border of the mesh and of every seam: lower distortion on some meshes (measured in `crates/ddg-uv/README.md`, "The boundary"), but the two sides of a seam get different lengths, and on some meshes the result overlaps itself or has mirrored triangles. Off keeps every border edge at its length. The status line and the warnings report the largest seam length difference |
| Cut Marked Seams | `use_seams` | on | cut the surface along edges marked as seams before placing cones |
| Pinned Cones | `pinned_group` | (none) | vertex group; its vertices are cone sites that cost nothing |
| No Cones | `forbidden_group` | (none) | vertex group; no cone is placed on its vertices |
| Mark Cuts as Seams | `mark_seams` | off | afterwards, mark the edges the surface was cut along as seams (existing seams are kept) |

A vertex belongs to a vertex group for this purpose when its weight in the group is above 0. A
vertex in both groups is an error.

In Edit mode the operator leaves Edit mode, does its work on the mesh and re-enters Edit mode.
That is the only way to read an edited mesh and write its UVs in bulk; on a 300k-triangle mesh
the round trip costs Blender about 0.15 s (see below).

### Packing

The native module returns every chart at the scale the flattening gave it, all charts on top of
each other. The add-on lays them out in the unit square with `ddg_tools.pack_charts`:

1. Every chart gets the axis-aligned bounding box of its UVs. Charts are moved and scaled by one
   common factor; they are not rotated or mirrored, so relative chart sizes (texel density) are
   preserved.
2. Every box is padded by 1/256 of the UV square on each side: charts stay 1/128 apart and 1/256
   away from the border.
3. Boxes are placed on shelves: tallest first, left to right from the lower left corner, a new
   shelf above when the next box does not fit.
4. The scale is the largest for which everything fits.

A single chart ends up with its longer side spanning the square (minus the padding), anchored at
the lower left. The layout is deterministic and not space-efficient; run Blender's
*UV > Pack Islands* afterwards for a tight one. `crates/ddg-uv` has a packer of its own (it
rotates each chart into its smallest bounding box); the add-on does not use it, because
`ddg_tools.unwrap` is called with packing off and `pack_charts` does the layout.

### Feedback

After a run the status bar and the Info editor show one line, for example (Suzanne, defaults):

```
3 charts, 3 cones, 48 cut edges | angle error max 75.67°, mean 7.35° | scale 0.06x to 2.09x
of the mean | 2 flipped, 0 collapsed triangles | 0.024 s (native 0.024 s, Blender side 0.000 s)
```

(one line in Blender; the line breaks are added here). For a cube it is `1 chart, 8 cones, 7 cut
edges | angle error max 0.00°, mean 0.00° | scale 1.00x to 1.00x of the mean | ...`: eight cones
of 90° unfold a cube without distortion.

Charts, cones and cut edges are counts (a chart is a connected piece of the surface after the seams are cut; cuts through the cones stay inside a chart). The distortion figures are measured on the UVs that were
written: the largest and the area-weighted mean difference between a triangle corner's angle in
3D and in UV (0 when angles are preserved exactly, as for an isometry; a conformal map with cones keeps them only approximately), and the range of local length scale relative to its mean
(1.00x to 1.00x for an isometry). Flipped and collapsed triangles are mentioned when there are
any. Warnings from the native module are shown as separate warning lines. One of them tells you that a handful of vertices (thin spikes: vertices pulled far out of the surface) dominate the distortion, names them by their Blender vertex index, and suggests repairing the mesh first; on such a mesh the map is likely to be poor.

### Errors

A mesh that cannot be unwrapped ends in an error report, not in a traceback, and the mesh is
left as it was. The messages come from the Rust side and use Blender's vertex and face indices:

```
edge (0, 1) is shared by 3 faces: 0, 3, 6
vertex 1 is non-manifold: its incident faces form 2 separate fans
faces 2 and 3 traverse their shared edge (4, 5) in the same direction: the mesh is not
    consistently oriented (one of them is flipped)
vertex 12 is both pinned and forbidden; a vertex can be in only one of the pinned, prescribed
    and forbidden sets
Vertex group 'nope' does not exist on object 'Grid'
Mesh 'wire' has no faces
```

(Turn on *Preferences > Interface > Developer Extras* and the *Indices* overlay to see indices
in the viewport.)

## Build, install, test

Everything goes through scripts in `scripts/`; they need `cargo`, `uv` and `blender`.
All of them run Blender with `--factory-startup` and with its user directories and temporary
directory pointed at a scratch directory under `runs_tmp/`. Your Blender configuration is never
read or written.

```sh
scripts/blender-build-extension.sh          # wheel + dist/blender/ddg_uv_unwrap-<version>.zip
DDG_WHEEL_PORTABLE=1 scripts/blender-build-extension.sh   # the same, wheel built for glibc 2.28

scripts/blender-test.sh                     # build, install into a scratch Blender, run all tests
scripts/blender-test.sh --no-build          # ... with the zip that is already in dist/
scripts/blender-test.sh --bench             # ... and measure the overhead of the Blender side
scripts/blender-test.sh -k seams            # only the headless tests whose name contains "seams"

scripts/blender-test-official.sh            # the headless tests in official Blender 5.0.1
scripts/blender-test-official.sh 5.2.2 3.13 # ... in official Blender 5.2.2 (Python 3.13)
```

To install the zip in your own Blender: *Edit > Preferences > Get Extensions*, the menu at the
top right, *Install from Disk...*; or drag the zip into the Blender window; or

```sh
blender --command extension install-file -r user_default -e dist/blender/ddg_uv_unwrap-0.1.0.zip
```

Use a zip built with `DDG_WHEEL_PORTABLE=1` for anything you hand to someone else: the default
build links against the C library of the build machine (here it needs glibc 2.34), the portable
one needs glibc 2.28, which is what official Blender builds for Linux require.

### What the build does

1. `scripts/py-build-wheel.sh` builds `ddg_tools-<version>-cp311-abi3-<platform>.whl` with
   maturin. `abi3` means one binary for every Python from 3.11 on.
2. `scripts/blender_stage_extension.py` copies `blender/ddg_uv_unwrap/` to
   `dist/blender/stage/ddg_uv_unwrap/`, adds the wheel under `wheels/`, the repository's
   `LICENSE`, and completes `blender_manifest.toml` with the `wheels` and `platforms` lines for
   exactly the wheels it bundled. If the wheel reports `ddg_tools.IS_PLACEHOLDER` (the current one never does: it is `False`),
   the extension name and tagline are marked as a placeholder build.
3. `blender --command extension validate` and `... extension build` make the zip.

### Tests

`scripts/blender-test.sh` installs the zip into a fresh scratch user directory with
`blender --command extension install-file`, so the tests exercise the installed extension and
its bundled wheel, the way a user gets them.

- `blender/tests/test_addon.py` (headless, 24 tests): cube (quads), n-gons (cylinder caps, a
  concave hexagon, a non-planar pentagon), a face with a zero-area triangle (mollified, with a warning), Suzanne (three
  components that must not overlap in UV), torus, UV sphere, marked seams, marking cuts as
  seams, Edit mode (on a mesh edited after entering Edit mode), non-manifold meshes (three faces
  on an edge, a face inside a cube, a bow-tie vertex), a flipped face, a mesh without faces, a
  linked mesh, a native module that fails to load, vertex groups (including which pinned and forbidden vertices took effect), loose vertices and a wire edge marked as a seam (both ignored and counted),
  the status line, and thin spikes named by vertex index, a flat sheet that must come out as a similarity transform of its own x and y at every corner (which also checks the corner order). Every unwrap runs twice and asserts that vertex, edge, face and corner
  counts, vertex positions, edge vertices, the corner order of every face and the seam flags are
  bit-for-bit unchanged, that the active UV map has one finite UV per corner inside the unit
  square, and that both runs give identical UVs. One test feeds the checks deliberately broken
  data to show that they fail when they should.
- `blender/tests/ui_smoke.py` (windowed Blender on a private Xvfb display, 27 checks): the
  operator invoked from the viewport as a menu click would, registration for Adjust Last
  Operation, the operator panel and both menu entries actually drawn, undo and redo of UVs and
  of marked seams in Object mode and in Edit mode. Skipped with a notice if `xvfb-run` is
  missing.
- `blender/tests/bench_overhead.py`: the measurements below.
- `scripts/blender-test-official.sh` runs `test_addon.py` unchanged in Blender-as-a-Python-module
  (the `bpy` package the Blender Foundation publishes on PyPI, fetched by `uv`), with numpy
  installed, which is what official builds look like from the inside.

Observed on the development machine on 2026-10-12, all with the same zip and the same wheel (built without
`DDG_WHEEL_PORTABLE`):

| Blender | Python | numpy | headless tests | UI checks |
|---------|--------|-------|----------------|-----------|
| 5.0.1, Ubuntu's build | 3.14.4 (system) | absent | 24 of 24 | 27 of 27 |
| 5.0.1, official (`bpy` from PyPI) | 3.11.15 | 1.26.4 | 24 of 24 | not run |
| 5.2.2 LTS, official (`bpy` from PyPI) | 3.13.14 | 2.5.4 | 24 of 24 | not run |

## Cost of the Blender side

`scripts/blender-test.sh --bench`, Blender 5.0.1, Python 3.14.4, AMD Threadripper 3970X, best of
7 runs. Two meshes: a 388 x 388 grid of quads (151 321 vertices, 301 088 triangles, 602 176
face corners) and an icosphere with 8 subdivisions (163 842 vertices, 327 680 triangles,
983 040 face corners).

| step | quad grid | icosphere |
|------|----------:|----------:|
| read the mesh: positions, corners, face starts, Blender's triangles | 1.9 ms | 2.5 ms |
| the same when Blender still has to tessellate (first call) | 4.7 ms | 7.7 ms |
| read seams: none marked / 1510 resp. 2458 marked | 0.01 / 0.7 ms | 0.01 / 1.6 ms |
| write the UVs | 0.12 ms | 0.24 ms |
| **Blender side with default options** | **2.0 ms** | **2.7 ms** |
| read two vertex groups (a Python loop over all vertices) | 167 ms | 181 ms |
| mark 1503 resp. 3224 cut edges as seams (a Python loop over all edges) | 40 ms | 67 ms |
| Blender side with both vertex groups and Mark Cuts as Seams | 209 ms | 252 ms |
| leaving and re-entering Edit mode (Blender itself; only when run from Edit mode) | 124 ms | 182 ms |
| native module, without the pipeline: copy and validate the buffers, build the result | 9 ms | 17 ms |
| native module: measure the distortion of the result | 63 ms | 73 ms |
| native module: pack the charts | 9 ms | 18 ms |

With default options the Blender side is about 2 to 3 ms and the binding about 80 to 110 ms
(most of it the distortion measurement), on a mesh for which the solver is expected to take
seconds. The two optional inputs that cost real time are vertex groups (Blender has no bulk
access to group weights, so it is one Python loop over the vertices, about 1 µs each) and
marking seams.

## Things about Blender 5.0 this code depends on

- **`foreach_get` / `foreach_set` only copy memory when the buffer's element type is exactly the
  property's**, including signedness, and the property is stored as a plain array. Otherwise
  Blender converts item by item. On the meshes above: `loop_triangles.foreach_get("loops", ...)`
  takes 0.09 ms into `array('I')` and 39 ms into `array('i')`, because the property is unsigned;
  `loops.foreach_get("vertex_index", ...)` takes 45 to 75 ms whatever the buffer, while the
  attribute `".corner_vert"` gives the same values in 0.05 ms; `edges.foreach_get("use_seam")`
  takes 42 to 70 ms and the attribute `"uv_seam"` 0.01 ms. `core.py` therefore uses unsigned
  buffers for face starts, loop triangles and edge vertices, and attributes for positions,
  corner vertices and seam flags. The attribute names are verified once per session against the
  element API on a two-triangle mesh; if they ever stop matching, the element API is used.
- **Writing UVs through `uv_layer.uv` breaks undo.** That collection writes into the layer's
  array in place even while the array is shared with an undo step, so Undo "restores" the new
  UVs. Writing through the attribute with the UV map's name (or through `uv_layer.data`, which
  is 250 times slower) is correct. Found by the UI smoke test.
- **In Edit mode the UV arrays of the mesh datablock read as empty**, even after
  `Object.update_from_editmode()`. Hence the mode switch.
- **`array.array` has no bool type**, and bool properties need a bool buffer: seam flags go
  through `memoryview(bytearray).cast('?')`.
- **A windowed Blender writes `quit.blend` into the temporary directory on exit.** The scripts
  point `TMPDIR` at the scratch directory so that the user's own "recover last session" file in
  `/tmp` is not replaced.
- **`--factory-startup` ignores saved preferences**, including which extensions are enabled:
  after `extension install-file -e` the tests still have to enable the add-on in the session.
- Blender picks bundled wheels by file name. Blender 5.0 accepts a `cp311-abi3` wheel on any
  newer Python (it runs on 3.14 here); it only looks at the machine part of the platform tag,
  not at the glibc version, so a wheel that needs a newer glibc than the system has is installed
  and then fails to import. The operator reports that as "The native module 'ddg_tools' could
  not be loaded".
- The system Python of this machine's Blender (3.14) has no numpy; official builds have it. The
  add-on and the native module use neither.

## The variational cut: Start, Flow, Apply

Three more operators, offered where the unwrap operator is (Edit mode: **UV** menu; Object mode:
**Object** menu): `bpy.ops.uv.ddg_cut_start`, `uv.ddg_cut_flow`, `uv.ddg_cut_apply`. They cut a
mesh open where it flattens with the least distortion, with the flow of Sharp and Crane,
*Variational Surface Cutting* (ACM TOG 2018), implemented in `crates/ddg-variational-cuts`
**clean-room, from the text of the paper** (the authors' code was never consulted) and bound
through `ddg_tools.CutSession` (`crates/ddg-py/README.md`). A cut is advanced in several runs; an
applied cut becomes seams; further cuts can be added; **existing seams are fixed cuts** (the flow
cannot cross them and they stay seams).

Typical use: select some faces (or one vertex) in Edit mode, **UV > Cut: Start**, then
**Cut: Flow** one or more times (the cut in the wire object "<name> cut" moves), then
**Cut: Apply**, then **Conformal Unwrap**.

1. **Cut: Start** (`mode`, `radius`, `vertex`, `regions`, `cluster_iterations`, `attribute`). Creates the state as one float
   attribute per region on the face-corner domain, `ddg_cut_region_0`, `ddg_cut_region_1`, ... (a
   corner attribute, so that the two sides of a seam can differ; the second is the negative of
   the first), and the preview. `SELECTED`: the selected faces are region 0, the rest region 1.
   `CIRCLE`: the geodesic disk of `radius` (object units, scale not applied; 0: a tenth of the
   bounding box diagonal) around `vertex` (-1: the active vertex of the edit mesh) is region 0.
   `CLUSTERS`: the vertex normals are clustered into `regions` regions (2 to 16, at most
   `cluster_iterations` rounds): the start of an atlas (paper Sec. 7). `MATERIALS`: the
   material slot of each face is its region (one region per slot in use). `ATTRIBUTE`: the
   integer face attribute named by `attribute` labels the regions (one per distinct value,
   ascending). `SELECTED` and `CIRCLE` have two regions, the other modes up to 16.
2. **Cut: Flow** reads the attributes, runs `steps` steps (fewer with `use_stop_rule`: the
   library's stop rule), writes them back, replaces the preview, and reports the energy before
   and after, the cut length, the steps taken and the reason it stopped (status line; also in the
   mesh's custom property `ddg_cut_last_flow`). Running it again continues from the stored state.

   | property | meaning | default |
   |----------|---------|---------|
   | `steps`, `use_stop_rule` | number of steps; stop early by the rule | 10, on |
   | `energy` | Dirichlet, Scaled Dirichlet (curvature-based weight), Hencky | Hencky |
   | `length_mode`, `ripple_wavelength` / `length_weight` | the length weight from a ripple wavelength (a fraction of the square root of the area; ripples longer than it grow), or explicit (on the surface scaled to unit area). Computed at the start of each run | ripple, 0.25 |
   | `stages` | recommended (the default stages for two regions, the paper's stages with the library's flux recovery for more: the core crate's SPEC D20 says why), the default stages for any N, or the paper's | recommended |
   | `penalty_group`, `penalty_strength` | vertex group as penalty `eta` (paper Sec. 7.1): `eta = strength * length weight * group weight`: **the strength is multiplied by the length weight in force** (so it is in units of that weight), and the tooltips say so | none, 1 |
   | `hide_from_camera`, `hide_strength` | penalty where the active camera sees the surface: ray casts from the camera (anything in the scene can occlude) give the visibility, the library multiplies it by the facing term `max(n.e, 0)`; `eta += strength * length weight * that`. Perspective position of the camera, object-space normals | off, 1 |
   | `distortion_weights`, `area_weights` | `w_k`, `beta_k` per region (16 entries, the first N are used; the panel shows one row per region); `w = 0` for a region makes it a one-sided patch (paper Fig. 14) | all 1, all 0 |
   | `symmetry`, `symmetry_tolerance` | mirror symmetry about the object's X, Y or Z plane; classes from mirrored vertices matched within the tolerance (fraction of the bounding box diagonal); unmatched vertices are reported | none |
   | `use_uniform_curvature`, `uniform_curvature` | replace the Gaussian curvature by a constant density (paper Sec. 7.4) | off |
   | *Advanced Options* (a toggle in the panel; `layout.panel` crashed Blender 5.0.1 in the redo popup) | cap and smoothing of the feature size, step size rule and its numbers, crossing offset, stop tolerances and window, the individual stages (normal derivative, length step, extension, speed scaling, redistancing, penalty rule) | the library's defaults |

   For more than two regions the report adds the number of patches (connected pieces in one
   region) and of triple points; the preview shows the network as curves between junctions.
   Not offered: a pointwise energy given by a function.
3. **Cut: Apply** (`unwrap_after`). Inserts the cut into the mesh with bmesh, marks the new
   edges as seams, removes the preview and the attributes, and optionally runs **Conformal
   Unwrap** with its defaults.

   *How the insertion works.* The description of the cut that the session gives (points on
   polygon edges with a parameter, per polygon the paths with their points on diagonals of
   Blender's triangulation) is inserted in two steps. The points on edges are made with
   `bmesh.utils.edge_split` (the two sides of a seam can carry different points; both are put on
   the edge), which interpolates vertex data (groups, attributes) and corner data (UVs, corner
   attributes). Then every polygon is split along its paths with `bmesh.utils.face_split`,
   which gives each part the data of the polygon (material, smooth flag, face attributes); the
   points inside a polygon are created by `face_split` and get vertex-group weights, vertex
   data and corner data interpolated here from the two corners of their diagonal. Polygons that
   the cut does not cross (and whose edges got no point) are untouched, nothing is triangulated.
   *Junctions* (more than two regions). The description has, per junction, the polygon, the
   triangle of Blender's triangulation, its three corners and the barycentric weights; three
   paths end at it. Inside a polygon the network is a forest whose leaves are edge points and
   whose inner nodes are junctions. It is embedded by splitting the polygon with
   `bmesh.utils.face_split` along chains of paths that run from a node that is a vertex already
   (at first the edge points) to another one, through the junctions in between: `coords` of
   the call are the interior points of the paths and the junction positions
   (`sum(w_k x[corner_k])`); the junctions a chain passed through are vertices afterwards, and the
   remaining arms are split the same way from them. This gives the polygons the Rust
   construction gives (its pieces `(c_k, q_k, J, q_{k-1})` merged into parts), with the vertex
   sets compared in the tests. A junction gets vertex-group weights, vertex data and corner data
   blended with its three weights.
   *Difference from `CutSession.insert`* (the Rust reference): the result has the same
   vertices (positions), polygons (as vertex sets) and seams, but Blender numbers them in its own
   order, and the reference also hands out a triangulation of the new polygons, which Blender
   computes itself. Not supported: meshes with shape keys (refused with a message).

Tests: `blender/tests/test_cut.py` (see *Tests*): start, flow and apply on a UV sphere, a
subdivided cube, a torus, the monkey and, if present, `data/models/spot/spot/spot_quadrangulated.obj`;
attributes and preview exist and are removed by Apply; the result equals the Rust insertion; the
seams are exactly the inserted edges plus the old ones; Euler characteristic and manifoldness are
kept; a UV map, a vertex group, vertex and corner attributes, materials and a face attribute that
are linear functions of position or inherited are still those functions afterwards; no flipped
UV triangle after unwrap; `k` steps twice equal `2k` once (float32 attributes, within 1e-4 of the
size); a second cut with the first as seams; every error path has a readable message.
`scripts/blender-test.sh --bench` also times Start, Flow and Apply against the library time
(`blender/tests/bench_cut.py`).

### Credits

The cut operators are a front end for Sharp, Crane, *Variational Surface Cutting*, ACM
Transactions on Graphics 37(4), 2018, implemented from the text of the paper. No code of the
authors is used or was read. The flow runs by default on an intrinsic triangulation:
Sharp, Soliman, Crane, *Navigating Intrinsic Triangulations*, ACM TOG 38(4), 2019.

## Layout

| path | contents |
|------|----------|
| `ddg_uv_unwrap/blender_manifest.toml` | extension manifest; `wheels` and `platforms` are filled in at build time |
| `ddg_uv_unwrap/__init__.py` | registration and the menu entries |
| `ddg_uv_unwrap/operator.py` | the operator: properties, panel, mode handling, reports |
| `ddg_uv_unwrap/cut_operators.py`, `cut.py` | the three cut operators; attributes, preview, flow settings, bmesh insertion |
| `ddg_uv_unwrap/core.py` | mesh to buffers, the call into `ddg_tools`, UVs and seams back to the mesh |
| `tests/` | the tests described above |

## License

MIT, like the rest of the repository. Blender asks add-ons to be licensed compatibly with the
GNU GPL; MIT is. The bundled native module is built from MIT-licensed code only; no part of the workspace is
derived from Blender. It statically links third-party Rust crates under
permissive licences (MIT, Apache-2.0, BSD-2-Clause and similar). maturin puts a CycloneDX list
of the whole dependency tree into the wheel (`ddg_tools-<version>.dist-info/sboms/`).

## Credits

The add-on is a front end for two papers:

- Yousuf Soliman, Dejan Slepčev, Keenan Crane. *Optimal Cone Singularities for Conformal Flattening.*
  ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201367>
  (where the cones are placed)
- Rohan Sawhney, Keenan Crane. *Boundary First Flattening.*
  ACM Transactions on Graphics 37(1), 2017. <https://doi.org/10.1145/3132705>
  (how the surface is flattened through them)

The flattening code (`crates/ddg-bff`) is a port of the authors' MIT-licensed reference
implementation; its licence notice is in `crates/ddg-bff/LICENSES/BFF-MIT.txt`. The operator's
tooltip in Blender names both papers.

### Flipped triangles after the unwrap: the stored triangulation

On the spot quads with eight clusters the unwrap of an applied atlas had a few flipped
triangles by the library's count (2 to 5 at the cut offset 1e-3; the two-region cuts of the UV
sphere, cube, torus and spot and the cube and sphere atlases had none). **Diagnosis.** (1) The
hypothesis "Blender's retriangulated polygons leave the surface the flow saw" is mostly wrong: of
17 flipped triangles on five atlas runs (quads and triangulated spot) 15 lie inside one triangle of
the original triangulation, on its plane to 1e-7 of an edge (normal angle 0 to 0.02 degrees);
only 2 span an old diagonal (distance from the surface 0.013 of an edge, normal angle 11 and 43
degrees). On the triangulated spot, which has no diagonal to span, Blender's triangulation flips
2, 0 and 7 triangles. (2) They are thin triangles of split polygons (3D area 6e-4 to 1.8e-1 of
the mean, mostly under 1e-2, smallest angle mostly 0.1 to 0.2 degrees) next to a cut that passes
near a vertex. (3) Float32 positions are not the cause: unwrapping the library's own insertion
with its positions rounded to float32 flips nothing. (4) What decides is the choice of diagonals
in the split polygons: the library's triangulation of the same polygons flips nothing at any
offset (`CutSession.unwrap`, and this add-on with the library's triangles, below), Blender's
flips 1 to 5; choosing better diagonals ourselves (best smallest angle) flipped 5, 7 and 5 where
Blender's own flipped 4, 2 and 2. **Remedy.** Apply computes the Rust insertion of the same cut
(`CutSession.insert`), matches its polygons to the new Blender faces by vertex set and stores its
triangles for the polygons that were split or got extra vertices in the mesh property
`ddg_cut_triangulation` (a flat integer list: counts, then per polygon its ascending vertex
indices and its triangles as vertex indices). `core.unwrap_object` (so the Conformal Unwrap
operator, also after Apply with `unwrap_after`) uses these triangles for exactly those faces when
the number of vertices and faces is still what it was at Apply and the face has the same vertex
set; every other face, and every mesh without the property, keeps Blender's triangulation, so
nothing changes for meshes that never had a cut applied. Nothing is triangulated in the mesh. The
reports say what happened: Apply says how many faces have a stored triangulation, the unwrap's
status line says "library triangulation for N faces of the applied cut", and when the stored
triangulation does not fit any more (faces or vertices added or deleted since Apply) the operator
warns "Used Blender's triangulation for ... (the mesh was edited since the cut was applied ...)".
Moving vertices keeps it valid. A second Apply keeps the entries of faces it did not touch.
**Result.** Flipped triangles by the library's count after Apply and unwrap, with the stored
triangulation, at the offsets 1e-3, 0.01, 0.05, 0.1 and 0.2: zero in all 45 combinations of nine
cases (spot atlases with 8 clusters at weights 2, 1 and 6 clusters at 0.3, cube atlas, sphere
atlas, and the sphere, cube, torus and spot circles), and zero on two atlases of the triangulated spot at 1e-3. The cut
offset therefore returned to the library's default **1e-3** (`cut_offset`, on Start, Flow and
Apply; the flow keeps its own `min_crossing_offset`); a larger value only avoids very thin faces
at the price of moving the cut (largest displacement of a cut point on the spot atlases 0.0004 /
0.0024 / 0.006 / 0.014 at 0.01 / 0.05 / 0.1 / 0.2, cut length +0.003 % / +0.1 % / +0.4 % / +2 %).
What is not covered: the viewport shows UVs in Blender's own triangulation, which can differ from
the one the unwrap used, and a mesh edited after Apply gets Blender's triangulation back.

### The triangulation the flow runs on

Start has a **Triangulation** property (`triangulation`, `refine_bound`): **Intrinsic
(Delaunay)** (the default), **Intrinsic, Refined** (flips and refinement to a minimum angle,
default 25 degrees, at most 35: vertices are added) and **Mesh's Own** (Blender's triangulation,
as the versions before). The flow's stages degrade on badly shaped triangles (measured in
`crates/ddg-variational-cuts`: on a 165 degree sheared grid the speed is off by 110 % with the
paper's stages and the cut vanishes); flipping the triangulation of the same surface to Delaunay
(Sharp, Soliman, Crane, *Navigating Intrinsic Triangulations*, ACM TOG 2019, `crates/ddg-intrinsic`)
restores regular-grid behaviour at 0.01 to 0.06 s on the real meshes and usually adds no vertex;
refinement only helps where the input has too few vertices for thin geometry. The cut is carried
back onto the mesh's own faces whatever the choice (bends inside faces, several points per
edge, junctions of any degree), so nothing else changes for the user.

* **Optimal Delaunay** (`optimal_delaunay`, a checkbox next to Triangulation on Start, off by
  default) also builds the intrinsic optimal Delaunay triangulation (iODT, `crates/ddg-intrinsic`)
  on top of Intrinsic or Intrinsic, Refined: the hidden intrinsic triangulation gets more even
  triangle areas. It costs time: measured on spot (4 regions, 4 runs in parallel), Start takes
  4.5 s at 11.7k vertices and 20 s at 46.9k with it against 0.15 s and 0.55 s without, and ten flow
  steps 1.4 s and 6.8 s against 0.57 s and 2.5 s (about 2.5 times as long). It changes nothing on the mesh except through the cut
  the flow finds; the cut is carried back to the mesh's own faces as always. It is part of the
  stored option (`ddg_cut_option` reads `intrinsic+odt` or `refined:<degrees>+odt`) and so of the
  key of the kept session: a session built without it is never reused with it. Flow and Apply have
  no property of their own: they rebuild what Start stored. The iODT inserts vertices, so the
  state then usually lives in the bytes of `ddg_cut_state` too. With Mesh's Own the checkbox is
  greyed out and ignored, with the warning "Optimal Delaunay needs an intrinsic triangulation and
  is ignored with 'Mesh's Own'" (nothing is stored). The advanced Flow line reports the iODT
  status, vertices inserted and iterations.
* **State.** The per-corner attributes `ddg_cut_region_k` are the whole state when the
  triangulation added no vertex (the usual case of the default). When it did (Refined, or an
  export that had to split an edge) the whole state is the bytes of `CutSession.snapshot()` in the
  mesh property `ddg_cut_state` (base64: levels on every vertex of the triangulation, its
  fingerprint, the number of regions and the option), and the attributes stay as the view of the
  mesh's own vertices. The property `ddg_cut_option` (`intrinsic`, `refined:<degrees>`, `caller`)
  says which triangulation to rebuild. A stored state that does not belong to the mesh any more
  (its fingerprint differs: the mesh or its seams changed since Start) is refused with "The mesh
  or its seams changed since 'Cut: Start', so the stored cut no longer fits it ... Run 'Cut:
  Start' again." and nothing is overwritten. Editing the attributes by hand only works while
  there are no state bytes.
* **Cut Offset** (`cut_offset`, default 1e-3) is the smallest distance of a point of the cut from
  a vertex of the mesh's faces, as a fraction of an edge: for the intrinsic triangulations it is
  the library's `min_offset` (set when the session is built from the stored state, so Apply can
  use another one), for the mesh's own triangulation it is the crossing offset of the shown and
  inserted cut.
* **A curvature per triangle** is not offered by the operators; the binding refuses it on an
  intrinsic triangulation (`triangulation="caller"` takes it).
* **Apply** inserts the general cut. The description gives, per polygon, a graph whose nodes are
  edge points and junctions (a bend is a junction with two paths; any number of points per edge)
  and whose edges are paths through the polygon's diagonals. The graph is embedded by splitting the
  polygon with `bmesh.utils.face_split` along chains of paths: a shortest chain (breadth first)
  from a node that is a vertex already to another one, through the junctions in between; the
  junctions it passes become vertices and the next chains start from them. Of the pieces that
  contain both ends the one the chain lies in is chosen (`bmesh.geometry.intersect_face_point`);
  the new vertices are ordered along the chain by the edges of the split. A piece along a polygon
  edge (both ends on the same edge) only marks the sides between the two points as seams. Data
  is transferred as before (edge points by `edge_split`, interior points, junctions and bends
  by the barycentric weights the description gives). The Rust `insert` (rustdoc: cells in the
  barycentric coordinates of each triangle, parts across uncut diagonals) is the reference; the
  tests compare vertices (one to one by position), polygons (as vertex sets) and seams (exactly).
* **The flow's report** gains a line about the triangulation (vertices added, angle range) in the
  advanced view only (the `show_advanced` toggle); `ddg_cut_last_flow` always has the numbers.
* **Tests.** `scripts/blender-test.sh` runs `test_cut.py` five times in parallel, once per
  variant (`--triangulation intrinsic|caller|refined` changes the default of Start for that
  run, so that every expectation of the earlier tests is checked on each; `--optimal-delaunay`
  also turns Optimal Delaunay on by default: the variants `intrinsic-odt` and `refined-odt`):
  49, 48 (+1 skipped: the state bytes do not exist for the caller's triangulation), 49, 49 and
  49 tests, with the same assertions on every variant. (Until 2026-10-17 the `-odt` variants
  skipped `the_cut_offset_moves_points_away_from_vertices` and only printed the flipped UV
  triangles of the spot atlas: with iODT the 8-cluster atlas, weight 2, 60 steps, cut offset
  0.2, had 4 of 7354 flipped triangles. The cause was the order of the merge of bends and the
  rule of the cut offset in `ddg_uv`, see "Bends and flipped UV triangles" below; now 0 of 7262.)
  `the_cut_offset_moves_points_away_from_vertices` measures the smallest triangle on the stored
  triangulation (the one the unwrap uses; Blender's own triangulation of the split polygons is
  printed: it does not follow the cut and gave a ratio of 13 instead of 20 with iODT on the
  refined triangulation; the stored one 169, and 351 on the intrinsic one). New tests: the
  option string, session key and stored option of the flag; an iODT cut on a sphere (intrinsic) and a
  subdivided cube (refined) through Start, Flow, Apply and unwrap with no flipped UV triangle;
  the flag ignored, with a message, on Mesh's Own. New: the 153 degree and
  165 degree sheared grids and the needles-and-caps grid of `bad_mesh_study.rs` (built in the test)
  through Start, Flow, Apply and unwrap; a refined session with state bytes, the stale-state
  message and the unchanged bytes after a refused run; the triangulation line.
* **Bends and flipped UV triangles.** The cut that comes back from an intrinsic triangulation
  has bends inside faces; a bend that deviates from the straight line between its neighbours by
  less than **Merge Bends** (`merge_tolerance`, default 0.1, 0 keeps every bend; on Start, Flow
  (advanced) and Apply) of the longest edge of its triangle is removed before the cut is inserted
  (`ddg_uv`'s `IntrinsicOptions::merge_tolerance`). The inserted cut can then differ from the one
  the flow found by up to that fraction of an edge; the advanced Flow line and the Apply message
  say how many bends were merged and the largest displacement. Flipped UV triangles after Apply
  and unwrap, by the library's count (our own count on the triangulation the unwrap used agrees except
  in two refined cells before merging: 4 against 1 and 2 against 1), **before** (merge 0) and **after** (0.1), `test_cut.py`, spot atlas with 8 clusters at the
  length weights 2 and 10, cut offsets 1e-3 and 0.2:

  | case | intrinsic | refined | caller | intrinsic + iODT | refined + iODT |
  |------|----------:|--------:|-------:|-----------------:|---------------:|
  | spot w2, offset 1e-3 | 1 / 0 | 2 / 0 | 0 / 0 | 39 / 0 | 12 / 0 |
  | spot w2, offset 0.2 | 1 / 0 | 2 / 0 | 0 / 0 | 50 / 0 | 52 / 0 |
  | spot w10, offset 1e-3 | 0 / 0 | 0 / 0 | 0 / 0 | 15 / 0 | 25 / 0 |
  | spot w10, offset 0.2 | 0 / 0 | 4 / 0 | 0 / 0 | 51 / 0 | 36 / 0 |
  | sphere atlas 4 regions | 0 / 0 | 4 / 0 | 0 / 0 | 11 / 0 | 10 / 0 |
  | sphere circle | 0 / 0 | 0 / 0 | 0 / 0 | 17 / 0 | 19 / 0 |
  | the three bad grids (153, 165 degrees, needles and caps) | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

  Zero is asserted for every triangulation and case. Until 2026-10-17 the refined cases at offset
  0.2 kept one flipped triangle each (a bend "0.103 of an edge" from its chord) and iODT at 0.2
  kept a few: the cut offset moved every point of the cut first, and the merge measured its
  displacement at the moved points, so that a bend the offset had pushed off a side counted as
  displaced and stayed. The merge now measures on the cut the flow found, and removes the
  excursions of the cut across a side whose points the offset moved (`crates/ddg-uv/README.md`,
  "A large cut offset, and the optimal Delaunay triangulation").

### Cost of the cut operators

`scripts/blender-test.sh --bench` (`blender/tests/bench_cut.py`, run alone), UV sphere, 32 514
vertices, 32 768 faces, 130 560 corners; "library" is the time inside `ddg_tools.CutSession` (the
intrinsic triangulation included), "add-on" the rest (mesh buffers, attributes, preview object,
bmesh):

| triangulation | operator | total ms | library ms | add-on ms |
|---------------|----------|---------:|-----------:|----------:|
| caller | Start | 111.8 | 92.9 | 18.9 |
| caller | Flow (10 steps) | 641.6 | 577.9 | 63.7 |
| caller | Apply | 122.2 | 22.6 | 99.6 |
| intrinsic | Start | 188.1 | 170.5 | 17.6 |
| intrinsic | Flow (10 steps) | 595.0 | 527.6 | 67.4 |
| intrinsic | Apply | 128.1 | 22.6 | 105.6 |
| refined | Start | 386.8 | 367.4 | 19.5 |
| refined | Flow (10 steps) | 823.0 | 788.8 | 34.2 |
| refined | Apply | 81.4 | 21.6 | 59.9 |

The add-on's share was 394 to 415 ms per Flow and 230 to 240 ms per Apply before this work
(first measurements: 36 and 72 ms). Where it went (cProfile inside Blender, 32.5k-vertex
sphere, the profiler inflates Python): **Flow**: the patch count, a Python union-find over all
130k corners with 261k `join` calls (0.67 s profiled, about 330 ms real); then the session built
anew at every call (0.17 s: the intrinsic triangulation, the symbolic analysis, the factorization)
for a state that had just been written; the mesh buffers, the attributes and the preview about
50 ms. **Apply**: the session built anew (0.17 s), the table of the library's triangulation of the
split faces (0.16 s: a Python decoding of every vertex of the insertion and of every triangle
of every polygon), the `sorted` of every face's vertex tuple to validate the stored entries (62
ms), the per-face loop for sizes and starts (15 ms), then the bmesh work itself (about 90 ms: the
chains, the data transfer, `from_mesh`, `to_mesh`). What was changed: the patch count is
`CutSession.n_patches` (Rust, the same union-find; its cost is inside the library time); the library
returns the triangulation of the polygons that changed (`InsertedMesh.changed_*`) so the add-on
decodes only the vertices of those; the validation of stored entries runs only when an older
triangulation exists; the sizes and starts come from `foreach_get`; and **a session is kept**
between operator calls (`cut._CACHE`, four entries): the mesh buffers are extracted every time
and compared by hash (positions, corners, polygons, the triangulation, seams), together with the
option and the offsets; the state is compared with what the last call wrote (attributes as float32
bytes, or the state string). If only the attributes differ they are handed to the session
(`session.state = ...`; a refusal rebuilds it and so gives the readable message), if anything else
differs the session is rebuilt. Start, Flow and Apply therefore share one session and Apply's
library time falls from 100-190 ms to 22 ms. What is left in the add-on: Flow about 65 ms (buffers
and hash 15, reading the attributes for the comparison 16, writing them 16, the preview, the
mesh update), Apply about 100 ms (the bmesh chains and data transfer 60 to 90, `from_mesh` and
`to_mesh` 20, buffers and attributes 30); the refined rows are smaller in Flow because the
library's time grows. Not removed: the chains loop of Apply is Python (a bmesh call per chain).
