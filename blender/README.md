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

## Developability: Make Developable, Show Developability, Developable Seams + UV

Three operators after Stein, Grinspun, Crane, *Developability of Triangle Meshes* (ACM TOG
2018; implemented clean-room from the text; crates `ddg-developability`, `ddg-py`; contract
`crates/ddg-developability/SPEC.md` section 7.2). All work on the active mesh from Object mode
and from Edit mode (the mesh is read in Object mode and Edit mode is re-entered), use Blender's own
triangulation of the polygons (`mesh.loop_triangles`), support undo and Adjust Last Operation, and
refuse non-manifold or inconsistently oriented meshes with a message naming polygons and
vertices. Menus: **Object** (all three), **Mesh** in Edit mode (Make Developable, Show
Developability), **UV** in Edit mode (Developable Seams + UV).

**Make Developable**, `bpy.ops.object.ddg_make_developable`: a flow that moves the vertices until
the faces around each vertex fold along straight lines (rounds of L-BFGS, a 4-1 subdivision
between rounds, small-angle remeshing). Runs in place.

| property | Python name | default | meaning |
|----------|-------------|---------|---------|
| Energy | `energy` | Eigenvalue | `COVARIANCE` (E-lambda), `COMBINATORIAL` (E_P), `COMBINATORIAL_MAX`, `COVARIANCE_MAX`, `INTRINSIC`, `INTRINSIC_MAX` |
| Method | `method` | L-BFGS | `LBFGS`, `GRADIENT_DESCENT`, `FIXED_STEP` |
| Iterations | `iterations` | 500 | optimizer iterations per round |
| Subdivision Rounds | `subdivisions` | 2 | 4-1 subdivisions; the flow has this + 1 rounds. Each one quadruples the triangles |
| Remesh | `remesh` | on | repair triangles with a corner angle below the small angle by edge flips and collapses |
| Small Angle | `small_angle` | 1 degree | the threshold of the remeshing |
| Step Cap | `max_displacement` | 0.5 | no vertex moves farther than this many mean edge lengths per step |
| Include Boundary | `include_boundary` | on | count boundary vertices in the energy |
| Omit Valence 3 | `omit_valence_three` | on | leave out interior vertices with three faces |
| Pin Selected | `pin_selected` | off | vertices selected (in Edit mode) do not move |
| Pinned Group | `pinned_group` | (none) | vertex group (weight above 0) of vertices that do not move |
| Triangulate | `triangulate` | off | replace the mesh by the algorithm's triangle mesh even if the topology is unchanged |
| Write Energy Preview | `write_energy` | off | afterwards write the preview attributes of Show Developability |

**What happens to the mesh (D24).** If no subdivision ran and no remeshing operation acted, and
Triangulate is off, the new positions are written into the existing mesh: polygons (quads), UV maps,
attributes, materials are untouched, and the report says the largest non-planarity of polygons with
more than three corners (a warning above 1 % of a side: the flow makes the *triangulation*
developable; n-gons follow only as far as they are planar). Otherwise the mesh data is replaced by
the triangle mesh the flow returns: the material index and smooth flag of every triangle come from
the polygon it descends from; **UV maps, vertex groups' weights and other attributes are dropped**
and a warning lists them. Shape keys are refused (the flow would move only the basis). Vertices
that no polygon used are dropped by a replacement.

**Progress and cancelling.** The operator is synchronous (no existing operator of this add-on is
modal). While it runs it reports to Blender's progress cursor (`window_manager.progress_*`, up to
rounds x iterations), and Ctrl-C in the terminal stops it, leaving the mesh unchanged. There is no
Escape key cancel: a flow of 500 iterations x 3 rounds on 20k triangles takes seconds; measured
times are below. (The Python function accepts a `progress` callable that can stop the flow.)

**Show Developability**, `bpy.ops.object.ddg_show_developability` (a third, small operator; the
Make Developable option `write_energy` does the same after the flow): writes the energy of every
vertex as the point color attribute `developability` (white below 1e-6, through yellow to red at
1e-1 and above, logarithmic) and the float values as `developability_energy`, and makes the color
attribute active (set the viewport's solid color to *Attribute*). Properties: `energy`,
`include_boundary`, `omit_valence_three`, `attribute_name`. The mesh is not changed.

**Developable Seams + UV**, `bpy.ops.uv.ddg_developable_seams`: (1) the library chooses the cut
(terminals: vertices whose energy exceeds the threshold, and interior vertices whose angle defect
exceeds the tolerance; a Steiner-tree-like cut graph through them along **polygon edges**); (2) the
edges are **marked as seams** (existing seams are kept); (3) the UV map is made by the same
pipeline as Conformal Unwrap with *Cut Marked Seams* on, so cones, extra cuts, packing and the
overlap estimate come from it. Vertices are not moved. The paper's own per-chart flattening is
`ddg_tools.developable_unwrap` (Python only): it keeps lengths along the cut but a chart can
overlap itself, which is why the operator uses the atlas pipeline.

| property | Python name | default | meaning |
|----------|-------------|---------|---------|
| Energy | `energy` | Eigenvalue | the energy that chooses the vertices to cut through |
| Threshold | `epsilon` | 0.001 | a vertex above this is cut through (paper: 0.01 to 0.001) |
| Cut Through Cone Vertices | `use_defect_rule` | on | also cut through interior vertices with a large angle defect (D30); off is the paper's rule |
| Defect Tolerance | `defect_tolerance` | 0.1 rad | the angle defect above which that happens |
| Length Weight | `length_weight` | 1 | weight of the edge length in the cost of a seam edge |
| Alignment Weight | `alignment_weight` | 0 | weight of the edge's misalignment with the creases at its ends (0: length only) |
| Keep Marked Seams | `use_existing_seams` | on | marked seams are cut too and the cut is planned around them |
| UV Map | `uv_name` | Developable | name of the UV map written (created if missing, then made active) |
| Lambda | `lambda_` | 0.1 | as in Conformal Unwrap |
| Positive Cones Only | `positive_cones_only` | on | as in Conformal Unwrap |
| Mark All Cuts | `mark_all_cuts` | off | also mark the extra cuts the unwrap made |

The status line reports the seam edges (and how many were newly marked), terminals, charts, cones,
the angle error and scale range, the **overlap** estimate, and the mirrored and collapsed UV
triangles.

Use it on a mesh that was made developable: on an untouched mesh a large part of the vertices has
energy above the threshold and the seams are everywhere.

### Developability tests and measurements

`blender/tests/test_developable.py` (headless, 18 tests, one skipped without an undo stack, run by `scripts/blender-test.sh` in
parallel with the cut tests; `--source DIR --site DIR` for development without building the
extension): subdivided cube (quads kept, energy down), noisy plane, UV sphere (all energies and
methods), cylinder, the replacement branch (parent material, dropped UV map listed), Triangulate,
pinned vertices (selected in Edit mode and a group), the preview attributes, readable errors
(non-manifold, shape keys, empty name), seams marked and a UV map written with no mirrored
triangle on cube, cylinder, sphere and spot, Edit mode, existing seams kept. Undo, redo and
Adjust Last Operation of both branches and of Seams + UV are checked in `ui_smoke.py` (windowed:
headless Blender has no undo stack, that test skips with a notice).

### Measured

End to end through the Blender operators (`bpy.ops`-level functions of
`developable_operator.py`), headless, Blender 5.0.1, the wheel in the installed extension:
`scripts/e2e-developable.sh` runs **Make Developable** with the defaults (2 subdivisions, 500
iterations per round, `E_lambda`, remeshing on) and with 1 subdivision, and then **Developable
Seams + UV** with its defaults (threshold 1e-3, defect rule 0.1 rad), on the models of
`data/models` (spot, bob, blub: the quadrangulated, triangulated and control meshes, imported as
OBJ) and on Blender primitives (the independent cases run in parallel, one Blender each;
`blender/tests/e2e_developable.py`; `scripts/e2e_table.py` makes the tables from the
`.jsonl` of the runs). The rows "0 x N" and "N x 5000" are longer flows run to see what a
converged flow does. Columns: *iterations done* (all rounds), *rounds end* (the optimizer status
of each round: `max` = ran out of iterations, `grad` = gradient tolerance, `LS` = the line search
could no longer find a decrease, at an energy of 1e-7 to 1e-13 that is convergence, nothing was
`STOPPED` by the D32 rule in any run), *energy* before and after (`E_lambda` of the paper, a sum over
vertices, so it grows with the mesh), *thin faces* (the library's `bad_faces`), seam edges and
terminals of the cut, charts, cones of the UV pipeline, *overlap* (the pipeline's sampled
estimate of the UV area covered more than once), *mirrored* UV triangles (clockwise in UV), and
seconds. **The seconds of the rows "2 sub x 500", "1 sub x 500" and "none" were re-measured on
2026-10-19 alone on a quiet machine** (AMD Ryzen Threadripper 3970X, 32 cores / 64 threads, load
0.1 to 1.5; one Blender at a time, one rayon thread, `scripts/e2e-developable.sh --no-sweep` with
`E2E_JOBS=1 RAYON_NUM_THREADS=1`), single-thread: spot (92k triangles after two subdivisions) 68
to 70 s, bob 102 to 110 s (170k), blub 167 to 169 s (226k); the control meshes and primitives 0 to
8 s. The seam and terminal counts of those rows agree with the earlier runs except for
`prim_sphere` (new numbers in its two rows; the code has changed since the earlier runs). **The
seconds of all other rows (the longer flows "N x 5000", "0 sub", "3000", the converged control
meshes) are from the earlier runs under load (40 to 350) and are not benchmarks**; they have not
been re-measured. The paper reports "seconds to a few minutes" single-threaded for its examples
(control-mesh-sized inputs, 1k to 5k triangles): our control meshes take 2 to 6 s for the
defaults, the dense models 1 to 3 minutes.

| mesh | flow (subdivisions x iterations per round) | iterations done | rounds end | energy | thin faces | seam edges | terminals | charts | cones | overlap | mirrored | flow s | seams+UV s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| spot_quad | 2 sub x 500 | 1500 | max/max/max | 111 to 5.47 | 2 | 2870 | 1222 | 1 | 0 | 2.9 % | 0 | 70 | 0.73 |
| spot_quad | 1 sub x 500 | 1000 | max/max | 111 to 8.84 | 2 | 2067 | 1538 | 1 | 0 | 15.2 % | 0 | 16 | 0.16 |
| spot_quad | 1 sub x 3000 | 6000 | max/max | 111 to 3.77 | 1 | 1303 | 907 | 1 | 0 | 10.6 % | 0 | 312 | 0.4 |
| spot_quad | 0 sub x 5000 | 5000 | max | 111 to 5.24 | 0 | 902 | 828 | 1 | 0 | 4.9 % | 0 | 61 | 0.1 |
| spot_quad | 0 sub x 20000 | 20000 | max | 111 to 1.96 | 2 | 511 | 456 | 1 | 0 | 15.9 % | 0 | 435 | 0.8 |
| spot_quad | none | - | - | - |  | 2617 | 2593 | 1 | 0 | 11.0 % | 0 | - | 0.04 |
| spot_tri | 2 sub x 500 | 1500 | max/max/max | 112 to 5.95 | 0 | 2522 | 1310 | 1 | 0 | 20.0 % | 0 | 68 | 0.81 |
| spot_tri | 1 sub x 500 | 1000 | max/max | 112 to 10.2 | 0 | 2049 | 1606 | 1 | 0 | 9.9 % | 0 | 15 | 0.16 |
| spot_tri | none | - | - | - |  | 2674 | 2663 | 1 | 0 | 16.0 % | 0 | - | 0.05 |
| spot_control | 2 sub x 500 | 1500 | max/max/max | 101 to 0.245 | 0 | 355 | 136 | 1 | 0 | 37.6 % | 0 | 3.0 | 0.04 |
| spot_control | 2 sub x 5000 | 15000 | max/max/max | 101 to 3.88e-07 | 9 | 296 | 65 | 1 | 0 | 34.2 % | 0 | 131 | 0.1 |
| spot_control | 1 sub x 500 | 1000 | max/max | 101 to 0.386 | 0 | 211 | 142 | 1 | 0 | 33.6 % | 0 | 0.9 | 0.01 |
| spot_control | 1 sub x 5000 | 10000 | max/max | 101 to 1.69e-07 | 3 | 148 | 65 | 1 | 0 | 30.8 % | 0 | 9 | 0.0 |
| spot_control | 0 sub x 5000 | 5000 | max | 101 to 8.02e-08 | 2 | 74 | 65 | 1 | 0 | 30.7 % | 0 | 2 | 0.0 |
| spot_control | none | - | - | - |  | 182 | 183 | 1 | 0 | 40.0 % | 0 | - | 0.01 |
| bob_quad | 2 sub x 500 | 1500 | max/max/max | 42.1 to 3.52 | 0 | 1392 | 284 | 1 | 0 | 1.0 % | 0 | 110 | 2.8 |
| bob_quad | 1 sub x 500 | 1000 | max/max | 42.1 to 6.63 | 0 | 3110 | 2610 | 1 | 0 | 6.0 % | 0 | 24 | 0.30 |
| bob_quad | 0 sub x 5000 | 5000 | max | 42.1 to 4.34 | 0 | 1127 | 1030 | 1 | 0 | 13.6 % | 0 | 223 | 0.1 |
| bob_quad | none | - | - | - |  | 4340 | 4331 | 1 | 0 | 11.5 % | 0 | - | 0.07 |
| bob_tri | 2 sub x 500 | 1500 | max/max/max | 42.9 to 3.24 | 0 | 1557 | 392 | 1 | 0 | 8.2 % | 0 | 102 | 1.8 |
| bob_tri | 1 sub x 500 | 1000 | max/max | 42.9 to 5.72 | 0 | 2235 | 1514 | 1 | 0 | 4.0 % | 0 | 22 | 0.31 |
| bob_tri | none | - | - | - |  | 4352 | 4342 | 1 | 0 | 14.5 % | 0 | - | 0.08 |
| bob_control | 2 sub x 500 | 1500 | max/max/max | 47.4 to 1.09 | 0 | 472 | 284 | 1 | 0 | 6.9 % | 0 | 5.9 | 0.07 |
| bob_control | 2 sub x 5000 | 15000 | max/max/max | 47.4 to 0.0179 | 1 | 191 | 65 | 1 | 0 | 8.0 % | 0 | 248 | 0.4 |
| bob_control | 1 sub x 500 | 1000 | max/max | 47.4 to 1.97 | 0 | 372 | 306 | 1 | 0 | 8.4 % | 0 | 1.6 | 0.02 |
| bob_control | 1 sub x 5000 | 10000 | max/max | 47.4 to 0.0347 | 3 | 102 | 58 | 1 | 0 | 7.3 % | 0 | 24 | 0.0 |
| bob_control | none | - | - | - |  | 335 | 334 | 1 | 0 | 14.7 % | 0 | - | 0.01 |
| blub_quad | 2 sub x 500 | 1500 | max/max/max | 46.3 to 2.93 | 0 | 2201 | 433 | 1 | 0 | 10.8 % | 0 | 169 | 2.0 |
| blub_quad | 1 sub x 500 | 1000 | max/max | 46.3 to 4.7 | 0 | 2263 | 1061 | 1 | 0 | 22.4 % | 0 | 36 | 0.42 |
| blub_quad | 0 sub x 5000 | 5000 | max | 46.3 to 4.07 | 0 | 1181 | 912 | 1 | 0 | 22.2 % | 0 | 232 | 0.2 |
| blub_quad | none | - | - | - |  | 5073 | 5024 | 1 | 0 | 11.8 % | 0 | - | 0.10 |
| blub_tri | 2 sub x 500 | 1500 | max/max/max | 45.6 to 2.88 | 0 | 2567 | 406 | 1 | 0 | 16.4 % | 0 | 167 | 2.2 |
| blub_tri | 1 sub x 500 | 1000 | max/max | 45.6 to 4.61 | 0 | 2278 | 1087 | 1 | 0 | 17.2 % | 0 | 34 | 0.41 |
| blub_tri | none | - | - | - |  | 5028 | 4982 | 1 | 0 | 18.9 % | 0 | - | 0.11 |
| blub_control | 2 sub x 500 | 1500 | max/max/max | 38 to 0.0142 | 2 | 207 | 70 | 1 | 0 | 9.3 % | 0 | 1.5 | 0.02 |
| blub_control | 2 sub x 5000 | 3131 | LS/LS/LS | 38 to 5.02e-13 | 2 | 210 | 53 | 1 | 0 | 24.7 % | 0 | 2 | 0.0 |
| blub_control | 1 sub x 500 | 1000 | max/max | 38 to 0.0232 | 3 | 104 | 65 | 1 | 0 | 18.7 % | 0 | 0.5 | 0.00 |
| blub_control | 1 sub x 5000 | 3100 | LS/LS | 38 to 6.33e-13 | 1 | 105 | 53 | 1 | 0 | 24.6 % | 0 | 1 | 0.0 |
| blub_control | none | - | - | - |  | 110 | 111 | 1 | 0 | 25.2 % | 0 | - | 0.01 |
| prim_cube | 2 sub x 500 | 1500 | max/max/max | 12.6 to 0.0267 | 0 | 198 | 71 | 1 | 0 | 20.7 % | 0 | 2.1 | 0.03 |
| prim_cube | 1 sub x 500 | 1000 | max/max | 12.6 to 0.0414 | 0 | 102 | 75 | 1 | 0 | 19.1 % | 0 | 0.6 | 0.01 |
| prim_cube | none | - | - | - |  | 35 | 8 | 1 | 0 | 0.0 % | 0 | - | 0.01 |
| prim_sphere | 2 sub x 500 | 1500 | max/max/max | 15.5 to 0.391 | 6 | 177 | 76 | 1 | 0 | 0.0 % | 0 | 6.8 | 0.12 |
| prim_sphere | 1 sub x 500 | 1000 | max/max | 15.5 to 0.851 | 2 | 302 | 281 | 1 | 0 | 3.4 % | 0 | 1.8 | 0.02 |
| prim_sphere | none | - | - | - |  | 481 | 482 | 1 | 0 | 0.0 % | 0 | - | 0.02 |
| prim_torus | 2 sub x 500 | 1500 | max/max/max | 15.5 to 0.0123 | 0 | 302 | 102 | 1 | 0 | 0.0 % | 0 | 7.4 | 0.12 |
| prim_torus | 1 sub x 500 | 1000 | max/max | 15.5 to 0.0225 | 0 | 159 | 128 | 1 | 0 | 0.0 % | 0 | 1.8 | 0.03 |
| prim_torus | none | - | - | - |  | 513 | 512 | 1 | 0 | 0.2 % | 0 | - | 0.01 |
| prim_monkey | 2 sub x 500 | 1500 | max/max/max | 148 to 0.878 | 0 | 837 | 330 | 3 | 0 | 9.3 % | 0 | 8.4 | 0.09 |
| prim_monkey | 1 sub x 500 | 1000 | max/max | 148 to 1.56 | 0 | 502 | 344 | 3 | 0 | 7.5 % | 0 | 2.2 | 0.02 |
| prim_monkey | none | - | - | - |  | 464 | 483 | 3 | 0 | 12.0 % | 0 | - | 0.02 |
| prim_cylinder | 2 sub x 500 | 161 | grad/LS/LS | 1.93 to -4.17e-15 | 0 | 200 | 44 | 1 | 0 | 0.0 % | 0 | 0.2 | 0.01 |
| prim_cylinder | 1 sub x 500 | 157 | grad/LS | 1.93 to -6.69e-16 | 0 | 100 | 44 | 1 | 0 | 0.0 % | 0 | 0.0 | 0.00 |
| prim_cylinder | none | - | - | - |  | 63 | 64 | 1 | 0 | 0.0 % | 0 | - | 0.01 |
| prim_plane | 2 sub x 500 | 1500 | max/max/max | 149 to 0.00168 | 0 | 116 | 25 | 1 | 0 | 6.3 % | 0 | 4.3 | 0.06 |
| prim_plane | 1 sub x 500 | 1000 | max/max | 149 to 0.00248 | 0 | 58 | 25 | 1 | 0 | 6.3 % | 0 | 1.2 | 0.01 |
| prim_plane | none | - | - | - |  | 225 | 269 | 1 | 0 | 17.6 % | 0 | - | 0.01 |

What the table says.

* **Nothing mirrors.** 0 mirrored UV triangles in every row with the default settings, on all
  models, control meshes included; no flow or default Seams + UV run raised an error (the only errors
  in the logs are the study's own: its settings sweep made more than the 8 UV maps Blender allows,
  fixed in the script since).
* **The flow of the defaults is not converged on the dense meshes.** Every round of the 500
  iterations ends by the budget (`max`). Spot (5856 triangles) needs more than 20 000 iterations
  to get to an energy of 2 (from 111); with 3 x 500 it is at 5.5. The control meshes (the paper's
  inputs) converge: spot's to 4e-7 in 15 000 iterations, blub's to 5e-13 in 3100 (`LS`), bob's to
  0.018.
* **A longer flow does give a much shorter cut.** Spot, quadrangulated: the untouched mesh has
  2617 seam edges (2593 terminals); after the defaults (3 x 500) 2870; after 1 x 3000 + 3000: 1303;
  without subdivision, 5000 iterations: 902; 20 000 iterations: **511** (the 2472 of the earlier
  run, which ended early, were those of a mesh that was not yet developable). The control meshes
  have 182 edges (spot), 335 (bob), 110 (blub) untouched and 74 to 355, 102 to 472, 105 to 210
  after the flows (65 terminals for spot, all of them cone-like vertices).
  On bob and blub the same: 4340 / 5073 edges untouched, 1127 / 1181 after 5000 iterations
  without subdivision.
* **The overlap does not follow.** Spot quadrangulated: 11 % untouched, 2.9 % after the defaults,
  15.2 % (1 x 500 + 500), 10.6 % (1 x 3000), 4.9 % (0 x 5000) and 15.9 % (0 x 20 000) for the longer
  flows: the estimate jumps around by a factor of three between neighbouring runs without a trend, and it stays at 31 % to 38 % on spot's control
  mesh whatever the flow does (40 % untouched). Bob and blub are between 1 % and 25 %. So a
  converged flow shortens the cut by 2 to 5 times, but the overlap of the unwrap is **not
  systematically lower**. One chart with a tree cut folds over itself where the cut has long
  branches; the paper deals with that by hand (SPEC D23). The cylinder, the sphere and the torus
  have 0 % overlap with the defaults (the sphere 4.6 % with 1 subdivision), and the cube goes from 0 % (8 terminals, 35
  seam edges) to 20 % after the flow (198 edges): the flow rounds the cube's corners, which then
  have to be cut through.
* **No cones in any run**: the cut is a tree, one chart (3 for the monkey, whose eyes are separate
  components), the pipeline needs no extra cones. That is also why `Lambda` has no effect (the
  seam sets and the overlap are identical for 0.1 and 1).
* The cylinder is made developable to rounding in 160 iterations (the status `grad`; its final
  energy is noise and slightly negative, SPEC D14); the flow rounds the sphere and the torus as far
  as the budget goes (energies 0.4 and 0.012).

The levers that were tried on the converged flows, **not applied to the defaults** (each row: the
same mesh after the same flow, only the settings of Developable Seams + UV change; "paper's rule"
is *Cut Through Cone Vertices* off):

| mesh (flow) | setting | seam edges | terminals | cones | overlap | mirrored |
|---|---|---|---|---|---|---|
| spot_quad (0 x 20000) | default (threshold 1e-3, defect 0.1) | 511 | 456 | 0 | 15.9 % | 0 |
| | threshold 1e-2 | 289 | 198 | 0 | 10.4 % | 0 |
| | threshold 1e-4 | 1046 | 992 | 0 | 51.6 % | 0 |
| | defect tolerance 0.3 | 481 | 422 | 0 | 12.7 % | 1 |
| | paper's rule (no defect terminals) | 447 | 374 | 0 | 14.4 % | 44 |
| | threshold 1e-2, defect 0.3 | 228 | 129 | 0 | 29.1 % | 0 |
| spot_quad (0 x 5000) | default | 902 | 828 | 0 | 4.9 % | 0 |
| | threshold 1e-2 | 362 | 252 | 0 | 17.3 % | 0 |
| | paper's rule | 881 | 801 | 0 | 4.6 % | 1 |
| bob_quad (0 x 5000) | default | 1127 | 1030 | 0 | 13.6 % | 0 |
| | threshold 1e-2 | 268 | 84 | 0 | 5.8 % | 0 |
| | threshold 1e-2, defect 0.3 | 251 | 45 | 0 | 0.7 % | 0 |
| blub_quad (0 x 5000) | default | 1181 | 912 | 0 | 22.2 % | 0 |
| | threshold 1e-2 | 335 | 196 | 0 | 18.0 % | 0 |
| | threshold 1e-2, defect 0.3 | 249 | 98 | 0 | 9.3 % | 0 |
| spot_control (2 x 5000) | default | 296 | 65 | 0 | 34.2 % | 0 |
| | defect tolerance 0.3 | 284 | 61 | 0 | 26.2 % | 0 |
| | paper's rule | 2 | 2 | 9 | 0.0 % | 42 |
| bob_control (2 x 5000) | default | 191 | 65 | 0 | 8.0 % | 0 |
| | defect tolerance 0.3 | 166 | 39 | 0 | 2.8 % | 0 |
| | paper's rule | 65 | 2 | 3 | 0.0 % | 17 |
| blub_control (2 x 5000) | default | 210 | 53 | 0 | 24.7 % | 0 |
| | defect tolerance 0.3 | 178 | 42 | 0 | 11.7 % | 1 |
| | paper's rule | 25 | 2 | 7 | 8.3 % | 68 |

* The threshold is the strong lever for the length of the cut: 1e-2 instead of 1e-3 gives
  1.4 to 4 times fewer seam edges on the dense meshes (spot 511 to 289, bob 1127 to 268, blub
  1181 to 335) and usually the same or lower overlap (but spot after 5000 iterations 4.9 % to 17 %);
  1e-4 gives 2 to 3 times more edges and more overlap. On the control meshes the threshold is
  irrelevant (the defect rule decides, all terminals are already defect terminals).
* The defect tolerance matters little on the dense meshes (481 against 511 edges for 0.3) and a
  little on the control meshes (spot 296 to 284 edges, overlap 34 % to 26 %; bob 8.0 % to 2.8 %;
  blub 24.7 % to 11.7 %, for one mirrored triangle on blub).
* **The defect rule is what keeps the control meshes from mirroring:** the paper's rule alone gives
  2 terminals, 7 to 9 cones and 17 to 68 mirrored UV triangles on the three control meshes (and 44
  on spot at 20 000 iterations), exactly the failure of the library's cut study (SPEC D30, §12.3),
  now seen in the Blender path.
* 1e-2 with a defect tolerance of 0.3 gives the shortest cuts (spot 228, bob 251, blub 249) and
  the lowest overlap on bob (0.7 %) and blub (9.3 %), but 29 % on spot at 20 000 iterations.
  There is no setting that is best on all models; the defaults are kept.

Limits of this study: one run per cell (the flow is deterministic; the overlap is a sampled
estimate); the dense meshes were converged by the flow only in the sense of the table; 1 or 2
subdivisions of the dense meshes were not run to convergence (hours). The earlier table of this
README (spot, 134 iterations then a stop) is replaced by these rows.

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
| `ddg_uv_unwrap/developable_operator.py` | Make Developable, Show Developability, Developable Seams + UV |
| `ddg_uv_unwrap/core.py` | mesh to buffers, the call into `ddg_tools`, UVs and seams back to the mesh |
| `tests/` | the tests described above |

## License

MIT, like the rest of the repository. Blender asks add-ons to be licensed compatibly with the
GNU GPL; MIT is. The bundled native module is built from MIT-licensed code only; no part of the workspace is
derived from Blender. It statically links third-party Rust crates under
permissive licences (MIT, Apache-2.0, BSD-2-Clause and similar). maturin puts a CycloneDX list
of the whole dependency tree into the wheel (`ddg_tools-<version>.dist-info/sboms/`).

## Credits

The add-on is a front end for these papers:

- Yousuf Soliman, Dejan Slepčev, Keenan Crane. *Optimal Cone Singularities for Conformal Flattening.*
  ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201367>
  (where the cones are placed)
- Rohan Sawhney, Keenan Crane. *Boundary First Flattening.*
  ACM Transactions on Graphics 37(1), 2017. <https://doi.org/10.1145/3132705>
  (how the surface is flattened through them)

- Oded Stein, Eitan Grinspun, Keenan Crane. *Developability of Triangle Meshes.*
  ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201303>
  (the developability operators; implemented clean-room from the text of the paper)

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
