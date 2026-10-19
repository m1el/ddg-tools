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
| Developability (Stein, Grinspun, Crane 2018) | `developable_hero.png`, `developable_energy.png` | Spot (2,930 vertices) before and after 6,000 L-BFGS iterations of E_λ, no subdivision: energy 112 → 6.22 in 38 s (one thread). Edges are shaded by their dihedral angle (white = flat, dark ≥ 30°): 567 of 8,094 edges are above 30° afterwards. The energy picture uses one fixed log color scale. |
| | `developable_bob.png`, `developable_bob_energy.png`, `developable_blub.png` | The same on bob (4,000 iterations, 42.9 → 6.17, 46 s) and blub (3,000 iterations, 45.6 → 4.82, 50 s). |
| | `developable_flow_spot.mp4/.gif`, `developable_flow_crumpled.mp4/.gif` | The flow as it runs: 1,500 iterations on spot (energy 112 → 11.4); on the crumpled sheet 2 subdivisions × 500 iterations (61.9 → 2.5e-4, 289 → 4,123 vertices). Each frame is the state of the one-call flow stopped at that iteration. Vertex energy on a fixed scale. |
| | `developable_energies.png` | The six energies on spot, 1,500 iterations each, no subdivision. Energy values are in each energy's own measure. |
| | `developable_crumpled.png` | `crumpleddevelopable` (289 vertices, intrinsically flat): 61.9 → 2.5e-4 in 4.6 s; edges above 30° fall from 354 of 736 to 60 of 11,862; straight rulings appear. The angle defect sum grows from 0.0015 to 0.58 rad: the flow does not keep the metric. |
| | `developable_defect.png` | Spot result colored by angle defect (±0.1 rad) with the cut (threshold 0.01: 381 seam edges, 280 terminals): the defect and the seams sit on the creases. |
| | `developable_cut_bob.png`, `developable_cut_sweep.png`, `developable_cut_crumpled.png` | The paper's own flattening (`developable_unwrap`) colored by the log scale factor. Bob at threshold 0.1: 224 seam edges, 51 terminals, 1 chart, mean \|u\| 0.176, 0 flipped triangles, no overlap; at 0.001 mean \|u\| falls to 0.0014 but the chart overlaps itself (14.4 %). The crumpled sheet flattens with mean \|u\| 1.5e-4, 0 flips, 0.2 % overlap. |
| | `blender_developable_*.png` | Make Developable and its Adjust Last Operation panel (spot quad mesh, 1,500 iterations, no subdivision), Show Developability as vertex colors, and Developable Seams + UV (seams in Edit mode, the UV editor, the options panel). |
| Blender extension | `blender_*.png`, `blender_*.mp4/.gif` | Screenshots of the add-on's operators: Conformal Unwrap, Cut: Start/Flow/Apply and Conformal Map to Sphere, with their option panels. |

## How it is made

- `run.sh` runs a compute script with the `ddg_tools` wheel from `dist/wheels/`, through `uv`.
- `compute/*.py` write `work/<scene>.npz`.
- `rust/` is a standalone helper. It depends on the crates by path and is not a workspace member. It provides what the Python module does not expose: the raw cones optimum for a λ sweep, snapshots of intrinsic triangulations traced over the input mesh, and the exact polyhedral distance. Build it with `CARGO_TARGET_DIR=work/target cargo build --release` inside `rust/`.
- `render/*.py` are Blender scenes. Run them as `blender -b --python render/<scene>.py -- <args>`; the arguments are in each script's docstring. All of them use `render/blib.py` (studio lights, shadow catcher, colormaps, tubes, spheres, camera fitting). Blender's Python finds numpy in `showcase/.venv`. Create that venv with `uv venv .venv --python /usr/bin/python3.14 && uv pip install --python .venv numpy scipy matplotlib <wheel>`.
- `post.py` composites renders onto the background, adds labels, and encodes the mp4 and gif. The `render/*_compose.sh` and `render/*_anim.sh` scripts put multi-panel images together.
- Developability: `compute/developable.py` (one flow, the input and the result with per-vertex energy, angle defect, dihedral angles, seams), `compute/developable_snap.py` (states of the flow at chosen iterations, one process each, by stopping the one-call flow through its progress callback), `compute/developable_cut.py` (seams and `developable_unwrap` for several thresholds). `render/developable.py` renders them (edges shaded by dihedral angle in the shader from per-face attributes, so the line width is constant in world space); `developable_stills.sh`, `developable_compose.py` (needs numpy: `uv run --with numpy`) and `developable_anim.sh` make the media. The scene-wide fixed scales are at the top of `render/developable.py`. `ui/shoot.sh developable 81` then `ui/finish_developable.sh` make the add-on screenshots.
- `ui/` drives a windowed Blender under Xvfb in a scratch user directory for the add-on screenshots: `install.sh`, then `shoot.sh <scene> <display>`, then `finish.sh`.

Several Blender renders sharing one GPU can run out of memory. `render/cuts_blender.sh` and `render/retry_frames.sh` retry, and with `SHOW_KEEP=1` they keep the frames already rendered.

## Developability: what the pictures showed

- Control meshes (spot, bob, blub, 100 to 330 vertices) are too far from developable: with 2 subdivisions and 100 to 500 iterations the flow turns them into spiky, folded shapes (spot grows spikes at the legs and tail). The stills therefore use the fine meshes with the flow run in place; control meshes appear only in the animation-free tests.
- On closed genus 0 and torus meshes the cut and flatten gives one chart with distortion concentrated where the surface is not developable (bob's hole). A lower threshold cuts through more vertices: mean \|u\| falls, but the chart then overlaps itself.
- `show_developability` leaves the float attribute `developability_energy` as the active color attribute, which Blender's Attribute shading draws black; `ui/drive.py` sets the color attribute `developability` active before the screenshot.
- The flow leaves two-face folds (maximum dihedral angle of 180° in every result) and, on the control meshes, triangles with a corner angle near 1e-4 rad.
