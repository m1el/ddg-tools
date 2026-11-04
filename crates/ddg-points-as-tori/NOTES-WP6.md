# WP6 notes: marching cubes, CSG, level-set sampling, morphology

Files: `src/{csg,mc}.rs`, `src/apps/{mod,morph}.rs`, `tests/{mc,apps}.rs`, `examples/reconstruct.rs`.
Run tests with `--release` (`cargo test --release -p ddg-points-as-tori --test mc --test apps`, about 10 s each).
Logs: `runs_tmp/pat-wp6/`.

## What exists

- `csg::ImplicitField` (`distance`, `distance_gradient`, `bbox`, batch `distances` / `distance_gradients`, `point_count_hint`),
  implemented for `TorusField` (batch via `evaluate`, own pool) and `Csg`. `Csg` = `Leaf(Arc<dyn ImplicitField>)`, `Union`, `Intersect`,
  `Subtract`, `Offset`, plus `SmoothUnion/Intersect/Subtract` (Quilez polynomial smin, off unless built). `[formula: pat-boolean]`,
  `[formula: pat-offset]`. Gradient of the active argument, ties to A.
- `mc::marching_cubes(&dyn ImplicitField, &McOptions) -> (McMesh, McReport)` (`[formula: pat-mc-interpolation]`).
- `apps::sample_level_set` (MC, area-uniform stratified sampling by a counter-based RNG, Newton projections, residual/move checks,
  replacement rounds) and `apps::morphology` (`Morph::{Erode,Dilate,Open,Close}`, refit by network + fit, report per step).

## Marching cubes design

- No 256-entry table. Per cell, each face gives directed segments (2 crossings: one; 4 crossings: the asymptotic decider), the
  segments chain into loops over the cube edges, loops are triangulated. Orientation comes from the segment direction (inside on the
  right of each face boundary traversal), outward normals, verified by the volume test.
- **Decider** (Nielson–Hamann): with inside pair P and outside pair Q of an ambiguous face, the inside diagonal is connected iff
  `f_P1 f_P2 > f_Q1 f_Q2`. It depends only on the four values through products and so is bitwise the same from both cells.
- **Triangulation of loops with 4 or more crossings.** A fan whose diagonals join two vertices of one cube face is forbidden (found
  by the random-lattice test: the neighbour across that face draws the same triangle or diagonal, giving duplicate faces and
  4-face edges). Among the allowed fans the one with the best worst-triangle quality wins; if none exists, a star around the loop
  centroid (a vertex private to the cell, key above all edge keys).
- Vertices are keyed by grid edge, numbered by sorted key, so output is independent of threads. Triangles follow block/cell order.
- **Narrow band.** Blocks of 8³ cells; corner values at block corners; a block is skipped if the balls `B(c, |v_c|/L)` around its
  corners cover it (single ball or 2 levels of octants). **Growth rule (ours):** a kept block whose face towards a skipped
  neighbour has a sign change keeps the neighbour (repeat). So `L` too small loses only isolated pieces, never leaves cracks.
  Test: torus with `L = 0.7` equals the dense mesh. The bound is stated in the report with the measured `max_edge_slope`.
- Edge refinement `refine_steps`: safeguarded Newton along each edge (bracket kept, bisection fallback); vertices stay on their edge.
- Grid values are stored densely (`(n+1)³` f64: 256³ = 136 MB, 512³ = 1.1 GB; limit 2^29 points, error above).

## Tests (`tests/mc.rs`, `tests/apps.rs`)

Sphere (res 16/33/64: watertight, chi 2, vertex error <= h/2, volume within 5 %), torus (chi 0, refinement error 7e-6 vs 6e-4
without), union of boxes (chi 2, volume within 3 %), box union minus a through-hole (chi 0), narrow band equal to dense, the
growth rule, thread determinism (1, 3, 8, 16; bits), 30,000 random lattices (3 to 5 points per axis, tens of thousands of
ambiguous faces; every mesh a valid `TriMesh`, boundary edges only on the box), decider on a checkerboard face, CSG
gradients against central differences (7 trees incl. smooth), ties to A, batch = single, bad input (7 option errors, no panic).
Apps: sampling residuals <= 1e-12 and unit normals on 3 levels, octant balance, thread independence, empty level set is an error;
MC and offset on a fitted sphere (512 Fibonacci points); morphology on the same: erode/dilate shift the radius within 0.02,
open/close keep the 512 input points within 0.03 of the surface and the volume within 15 %; golden `synth_torus` close/open within 0.04.

## Example (`examples/reconstruct.rs`)

`reconstruct <cloud.ply|.obj> <out-dir> [res,res,..] [threads] [refine] [bound|0=dense]`. Machine load from other agents was
about 30 (`runs_tmp/pat-wp6/uptime-before.txt`: load 33.6/27.9/19.8 at start, 30.9 after), 64 threads, the three runs one after the other.
Box `[-1.1, 1.1]³`, cloud normalized to `[-1,1]³`, default band (L = 2, block 8), no refinement.

| cloud (n) | build | 64³ | 128³ | 256³ | 256³ triangles | 256³ grid values |
|---|---|---|---|---|---|---|
| snail (4096) | 1.32 s | 0.055 s | 0.144 s | 0.478 s | 203,712 | 14.2 % |
| holey-bunny (12917) | 2.90 s | 0.069 s | 0.194 s | 0.688 s | 299,118 | 18.0 % |
| femur (3897) | 0.95 s | 0.032 s | 0.069 s | 0.231 s | 46,772 | 3.4 % |

(Time = marching cubes only, the field is built once; 90 % is field evaluation. Dense bunny 128³ takes 0.85 s: the band is 2.5x
faster, and the mesh is identical with `L = 0.5` thanks to growth.) OBJ files are in `runs_tmp/pat-wp6/out/`.

## Findings and SPEC proposals

1. **The fitted blend is far from 1-Lipschitz on the reference clouds.** `max_edge_slope` (largest |Δφ|/h along grid edges inside the
   kept blocks) grows like 1/h: bunny 24 / 48 / 99 at 64 / 128 / 256, snail 8 / 16 / 33, femur 4.8 / 9.8 / 22. That is a jump (a
   discontinuity of the blend, I suspect where the radius neighborhood `R_eval` changes or falls back to 32-NN), not a gradient. The
   §8.6 statement "L = 1 from ‖∇φ‖ ≈ 1" does not hold with the paper's radius options; with the `Knn` neighborhood I did not measure.
   WP14 should look at where the jumps are (size of the jump, whether it matters for the surface). Proposal: §8.6 says "L is a
   heuristic accelerator (default 2); the growth rule makes cracks impossible".
2. **Spurious closed components.** Dense and banded meshes agree, so these are in the field: femur 3 / 6 / 24 components at 64 / 128 / 256
   (chi 6 / 12 / 44), bunny 18 components (it has holes, chi -4 at 128), snail 3 / 1 / 4. Tiny blobs appear where φ crosses zero away
   from the surface. Worth a look in WP14 (eikonality map, which neighborhood option avoids them).
3. `Csg` has `Leaf(Arc<dyn ImplicitField>)` instead of `Field(TorusField)`, so baselines can be leaves; `Csg::field(f)` wraps a `TorusField`.
4. `ImplicitField` gained `point_count_hint` (default size of resampling) and the batch methods; the default batch is parallel on the
   current rayon pool, MC calls it inside its own pool.
5. `McOptions::isovalue` is on top of the field's own `EvalOptions::isovalue` (the latter is already subtracted by `distance`).
6. §8.6 "one Newton step" became `refine_steps` of safeguarded Newton along the edge (vertices stay in their cell; a projection along
   ∇φ can leave the cell and fold the mesh at coarse grids).
7. **`Cargo.toml` needs `ddg-core = { path = "../ddg-core" }`** (I added the line in the working tree; the file also holds WP4's
   uncommitted edits, so I did not commit it. Whoever commits Cargo.toml must include the line, otherwise the crate does not build).
8. The first part of this session built into the default `target/` by mistake (a missing `export`); no source effect.

## Open

- `sample_level_set` returns fewer points than requested only when rejections persist over 4 rounds; the count is reported.
- Morphology runs with `paper_eval` options for each refit (radius `2C/λ`); other settings via `MorphOptions::build`.
- Python/Blender/render (WP11-13) can use `marching_cubes`, `Csg`, `morphology` as they are; `McMesh::to_trimesh()` gives a `ddg_core::TriMesh`.
