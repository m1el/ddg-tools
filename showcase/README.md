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
| Blender extension | `blender_*.png`, `blender_*.mp4/.gif` | Screenshots of the add-on's operators: Conformal Unwrap, Cut: Start/Flow/Apply and Conformal Map to Sphere, with their option panels. |

## How it is made

- `run.sh` runs a compute script with the `ddg_tools` wheel from `dist/wheels/`, through `uv`.
- `compute/*.py` write `work/<scene>.npz`.
- `rust/` is a standalone helper. It depends on the crates by path and is not a workspace member. It provides what the Python module does not expose: the raw cones optimum for a λ sweep, snapshots of intrinsic triangulations traced over the input mesh, and the exact polyhedral distance. Build it with `CARGO_TARGET_DIR=work/target cargo build --release` inside `rust/`.
- `render/*.py` are Blender scenes. Run them as `blender -b --python render/<scene>.py -- <args>`; the arguments are in each script's docstring. All of them use `render/blib.py` (studio lights, shadow catcher, colormaps, tubes, spheres, camera fitting). Blender's Python finds numpy in `showcase/.venv`. Create that venv with `uv venv .venv --python /usr/bin/python3.14 && uv pip install --python .venv numpy scipy matplotlib <wheel>`.
- `post.py` composites renders onto the background, adds labels, and encodes the mp4 and gif. The `render/*_compose.sh` and `render/*_anim.sh` scripts put multi-panel images together.
- `ui/` drives a windowed Blender under Xvfb in a scratch user directory for the add-on screenshots: `install.sh`, then `shoot.sh <scene> <display>`, then `finish.sh`.

Several Blender renders sharing one GPU can run out of memory. `render/cuts_blender.sh` and `render/retry_frames.sh` retry, and with `SHOW_KEEP=1` they keep the frames already rendered.
