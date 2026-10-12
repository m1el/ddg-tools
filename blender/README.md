# DDG Conformal UV Unwrap (Blender extension)

One operator: take the active mesh and generate a UV map for it. The work is done in Rust by the
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
| Optimize Boundary | `optimize_boundary` | off | let the solver choose how much to stretch the border of the mesh and of every seam: lower distortion (measured below), but the two sides of a seam get different lengths, and on some meshes the result overlaps itself or has mirrored triangles. Off keeps every border edge at its length. The status line and the warnings report the largest seam length difference |
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

## Layout

| path | contents |
|------|----------|
| `ddg_uv_unwrap/blender_manifest.toml` | extension manifest; `wheels` and `platforms` are filled in at build time |
| `ddg_uv_unwrap/__init__.py` | registration and the two menu entries |
| `ddg_uv_unwrap/operator.py` | the operator: properties, panel, mode handling, reports |
| `ddg_uv_unwrap/core.py` | mesh to buffers, the call into `ddg_tools`, UVs and seams back to the mesh |
| `tests/` | the tests described above |

## License

MIT, like the rest of the repository. Blender asks add-ons to be licensed compatibly with the
GNU GPL; MIT is. The bundled native module is built from MIT-licensed code only; no part of the workspace is
derived from Blender. It statically links third-party Rust crates under
permissive licences (MIT, Apache-2.0, BSD-2-Clause and similar). maturin puts a CycloneDX list
of the whole dependency tree into the wheel (`ddg_tools-<version>.dist-info/sboms/`).
