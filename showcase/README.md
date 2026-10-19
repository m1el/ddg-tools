# Showcase

Stills and animations of every capability of the repository, rendered with Blender (Cycles) from
the results of the libraries themselves. The numbers in the captions are what these runs reported.
The media is in `media/`; everything in `work/` (data, frames, test renders) is regenerable and
ignored by git.

## Media

| Capability | Files | What it shows |
|---|---|---|
| Optimal cones (Soliman, Slepčev, Crane 2018) | `cones_hero.png` | The dataset airplane colored by the log scale factor `u`, with isolines. Red spheres are cone tips, blue are saddles. λ = 0.019 gives 90 cones. 2,349 vertices, certified optimum in 0.01 s. |
| | `cones_sweep.mp4/.gif` | λ sweeps from 0.002 to 0.3 and back: from 194 cones and distortion 0.017 down to 6 cones and distortion 1.07. All 72 solves take 1.05 s together. |
| UV unwrap (cones + BFF + cut) | `unwrap_spot.png` | Spot as a quad mesh. 9 cones, cut through them, one chart, 0 flipped triangles, max quasi-conformal error 1.78. The whole pipeline takes 0.06 s. |
| | `unwrap_morph.mp4/.gif` | Spot unfolds into its flat layout and folds back. The checker stays fixed. |
| Conformal map to the sphere (BFF) | `sphere_morph.mp4/.gif` | Spot keeps its texture while it morphs to its sphere image and back. |
| | `sphere_conformal.png` | The same checkerboard on the cow and on the sphere: the squares stay square. Quasi-conformal error mean 1.13, max 1.78, no flips. |
| Variational surface cutting (Sharp, Crane 2018) | `cuts_flow.mp4/.gif` | A cut grows from a geodesic disk to its rest shape on spot. 120 steps, energy 3.96 → 3.21, 0.59 s. |
| | `cuts_atlas.mp4/.gif`, `cuts_atlas.png`, `cuts_atlas_start_end.png`, `cuts_atlas_bob.png` | 8 regions start from normal clusters and smooth into a clean cut network. On spot: cut length 34.9 → 19.4 over 200 steps in 4.4 s. Bob is shown with 6 regions. |
| | `cuts_unwrap.png` | Spot cut along the atlas and flattened: 21 charts, 0 flipped triangles. |
| | `cuts_visibility.png` | Without and with the visibility penalty, seen from the camera and from the other side. The effect is partial: part of the cut stays on the visible side. |
| Intrinsic triangulations (Sharp, Soliman, Crane 2019) | `intrinsic_torus_stages.png`, `intrinsic_torus_delaunay.png` | A sliver torus with 880 faces. The input has min angle 0.6°. Intrinsic Delaunay takes 1,760 flips and reaches 16.4°. Refinement adds 138 vertices and reaches 25.8°. Intrinsic edges (orange) are traced over the input edges (grey). |
| | `intrinsic_torus.mp4/.gif` | Flips, then refinement, then the optimal Delaunay triangulation (+880 vertices). |
| Geodesic distance (heat method, Sec. 5.3) | `geodesic_camel.png`, `geodesic_camel.mp4/.gif` | camelbadnonuniform, source vertex 6061. Mean relative error against the exact polyhedral distance is 33 % on the input mesh and 2.1 % on its intrinsic Delaunay refinement (+3,756 vertices). Measured by `compute/geodesic_error.py`. |
| Direction fields (Sec. 5.5) | `field_spot_cross.png`, `field_bob_cross.png`, `field_blub_line.png` | Smoothest cross and line fields. Singularities are colored by index, and their indices sum to `n·χ`. |
| Steiner trees (Sec. 5.1) | `steiner_spot.png`, `steiner_spot_tree.png` | Spot's control cage: the mesh-edge MST has length 29.12; the intrinsic Steiner tree is 28.41 (−2.5 %) with 48 Steiner points. |
| Adaptive refinement (Sec. 5.4) | `adaptive_green_spot.png`, `adaptive_disk_face.png` | A Green's function on spot (2,930 → 17,927 vertices in 0.5 s). A harmonic map of the face mesh onto the disk. |
| Developability (Stein, Grinspun, Crane 2018) | `developable_hero.png` | The paper's recipe on a uniform input of about 1k triangles (MeshLab isotropic remeshing of the dense model): `E_P` on the input, then `E_λ` after each of 2 4-1 subdivisions, 500 iterations per round (the library's `FlowOptions::coarse_to_fine(2, 500)`). Flat shading; the seam curves (long connected components of the edges above 20°) are the red tubes; the third row is the angle defect, which sits on the seams. Spot: 1,102 → 16,854 faces, energy 32.3 → 6.34 (E_P), 16.1 → 6.27, 6.27 → 3.53, Hausdorff 0.045 of the diagonal, volume 0.96, 24 folds, 19 seam curves, 15 s. Bob: 934 → 14,346, 15.4 → 4.1, 12.7 → 5.33, 5.33 → 3.09, Hausdorff 0.058, volume 0.97, 0 folds, 19 curves, 12 s. Blub: 1,058 → 16,430, 22.5 → 1.59, 3.02 → 1.25, 1.25 → 0.691, Hausdorff 0.053, volume 0.89, 185 folds (its fins fold flat), 7 curves, 12 s. Timings are one thread on a machine with load 16 to 24, not quiet. |
| | `developable_flow_spot.mp4/.gif` | The same recipe on spot, 18 states (round, iteration) of the one-call flow, each from a run stopped at that iteration; the seam curves are found on every frame. Ears and legs start as spikes and settle into wedges and sheets. |
| | `developable_inputs.png` | The same recipe on spot's subdivision cage (372 faces → 5,138) and on the uniform 1k mesh. The cage: Hausdorff 0.221, volume 1.15, 168 folds; the head and ears fold into a flat double sheet. The uniform input: 0.045, 0.96, 24 folds. |
| | `developable_energies.png` | The six energies, each in all three rounds (3 × 500 iterations, subdivision between) on the uniform spot. `E_λ` alone flattens the legs into spikes (volume 0.91, 101 folds); the combinatorial energies keep more volume. Values are in each energy's own measure. |
| | `developable_crumpled.png` | `crumpleddevelopable` (512 faces, intrinsically flat): with the recipe 8,160 faces, 5 s, energy 33.9 → 5.2e-4; 80 of 12,368 edges above 30°, 0 folds, 6 seam curves at the rim. |
| | `developable_cut_spot.png`, `developable_cut_bob.png`, `developable_cut_crumpled.png` | Cut and flatten of the results: seams from `developable_seams` (epsilon 1e-3), unwrapped with the seams, without and with Split Overlapping Charts (`unwrap(..., split_overlaps=True)`). Spot: 1,117 seam edges, 1 chart with 2.4 % of its area covered twice; with splitting 10 charts (9 splits, 47 cut edges), overlap 0.027 %, 0 flipped triangles, 0.5 s. Bob: 917 seam edges, 1 chart, 9.0 % overlap; 21 charts (20 splits, 115 cut edges), 0.047 %. The crumpled sheet needs no seams: 1 chart, area distortion 1e-7, no overlap. Charts in colors, seams red, splits black. |
| | `blender_developable_*.png` | Make Developable and its panel with the Coarse Energy option (spot, uniform 1k input, 2 subdivisions × 500 iterations), Show Developability as vertex colors, and Developable Seams + UV (Split Overlapping Charts on: 7 charts, 6 splits, overlap 0.1 %). |
| Blender extension | `blender_*.png`, `blender_*.mp4/.gif` | Screenshots of the add-on's operators: Conformal Unwrap, Cut: Start/Flow/Apply and Conformal Map to Sphere, with their option panels. |

## How it is made

- `run.sh` runs a compute script with the `ddg_tools` wheel from `dist/wheels/`, through `uv`.
- `compute/*.py` write `work/<scene>.npz`.
- `rust/` is a standalone helper. It depends on the crates by path and is not a workspace member. It provides what the Python module does not expose: the raw cones optimum for a λ sweep, snapshots of intrinsic triangulations traced over the input mesh, and the exact polyhedral distance. Build it with `CARGO_TARGET_DIR=work/target cargo build --release` inside `rust/`.
- `render/*.py` are Blender scenes. Run them as `blender -b --python render/<scene>.py -- <args>`; the arguments are in each script's docstring. All of them use `render/blib.py` (studio lights, shadow catcher, colormaps, tubes, spheres, camera fitting). Blender's Python finds numpy in `showcase/.venv`. Create that venv with `uv venv .venv --python /usr/bin/python3.14 && uv pip install --python .venv numpy scipy matplotlib <wheel>`.
- `post.py` composites renders onto the background, adds labels, and encodes the mp4 and gif. The `render/*_compose.sh` and `render/*_anim.sh` scripts put multi-panel images together.
- Developability: the flows are the library's `develop` example (`cargo build --release -p ddg-developability --example develop`; the Python module has no per-round energies in the CLI sense, and the snapshots need one process per state). `uv run crates/ddg-developability/examples/coarse_inputs.py showcase/work/inputs 1000` makes the uniform inputs; `compute/dev_paper_runs.sh` runs the flows in parallel (one thread each) and `compute/dev_paper.py` measures them with `examples/quality.py` and extracts the seam curves; `compute/dev_paper_snap.sh` makes the animation states; `compute/dev_paper_cut.py` the seams and charts. `render/developable_paper.py` renders them (flat shading, seams as tubes; `developable_paper_stills.sh`, `developable_paper_cuts.sh`, `developable_paper_anim.sh`), `developable_paper_compose.py` assembles the figures. `ui/shoot.sh developable <display>` with `DDG_SPOT_OBJ=<uniform input>` then `DEV_DISPLAY=<display> ui/finish_developable.sh` make the add-on screenshots. The old `compute/developable.py` / `render/developable.py` (in-place flow on the dense models, edges darkened by dihedral angle) remain only as helpers.
- `ui/` drives a windowed Blender under Xvfb in a scratch user directory for the add-on screenshots: `install.sh`, then `shoot.sh <scene> <display>`, then `finish.sh`.

Several Blender renders sharing one GPU can run out of memory. `render/cuts_blender.sh` and `render/retry_frames.sh` retry, and with `SHOW_KEEP=1` they keep the frames already rendered.

## Developability: what the pictures showed

- The paper's look needs a uniform coarse input and `E_P` on it (`crates/ddg-developability/QUALITY.md`). Spot's subdivision cage under the same recipe folds its head into a double sheet (`developable_inputs.png`); blub's thin fins fold flat on any input.
- Flat shading and seam tubes show the pieces; the earlier pictures darkened every edge above 30°, and a seam that crosses the triangles diagonally then looked like hatching. The early frames of the animation still show many short seams: the input is not developable at all.
- On closed meshes the seams from the energy threshold 1e-3 leave one chart that overlaps itself (spot 2.4 %, bob 9.0 % of the area); the automatic split removes it with 9 and 20 extra cuts.
- The flow leaves two-face folds (blub 185, spot 24 edges above 170°) and keeps the angle defect on the seams: it does not preserve the metric.
- `show_developability` leaves the float attribute `developability_energy` as the active color attribute, which Blender's Attribute shading draws black; `ui/drive.py` sets the color attribute `developability` active before the screenshot.
