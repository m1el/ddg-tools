# NOTES-WP8: analytic shapes and data generators

Code: `src/analytic/*`, `src/generate/*`, tests `tests/analytic.rs`, `tests/generate.rs`. Data record:
`data/points-as-tori/analytic.json` (gitignored; regenerate with
`cargo test --release --test analytic -- --ignored --nocapture write_data`, which also prints saddle counts).

## Analytic
- 16 Quilez SDFs (`analytic::Shape2d`), closed-form gradients, polygons share one exact routine
  (edge normal on edges, `s (p-v)/|p-v|` at vertices, so gradients are exact on the boundary where the
  reference's `q/|q|` is 0/0). Ellipse: exact closest point (Eberly bisection), not 4 Newton steps.
  Heart straight edges closed form. Cross: swaps arms if `by > bx` (same set).
- Vesica with `d < 0` (about half the draws, from the reference's range `d = u 2r - r`) is the union of the
  two discs, not a lens (the Quilez formula). Kept as in the reference; documented in `Kind::Vesica`.
- `Solid` = revolution (about y, offset `o`) or extrusion, implements `ImplicitField`; gradients per
  pat-revolution / pat-extrusion / pat-extrusion-gradient. Revolution keeps only the profile part with
  `s >= -o` (the part left of the axis is not on the surface); the distance is exact outside.
- Sampling: 2D boundary traced by marching squares (2048 cells), vertices Newton-projected on `phi=0`;
  extrusion sides by arclength x z, caps by rejection in the profile; revolution by `length x mean radius`
  with exact inverse CDF inside a segment. Surface samples have `|phi| < 1e-10`; area relative error ~1e-5.
- 70 instances: `seed/5` picks `PAPER_KINDS`, `Rng::new(seed)` draws params, coin, `o`/`h`. 
- Saddle-rich set (`analytic::saddle_set`, 9 shapes + sphere and cylinder controls), curvature by central
  differences of the exact gradient. Saddle = `kmin < -1e-4 and kmax > 1e-4`; 20,000 area-uniform points each:

| shape | saddle fraction |
|---|---|
| torus (r .25, o .5) | 0.347 (exact 1/2 - r/(pi o) = 0.341) |
| thin-torus (.12, .6) | 0.444 |
| fat-torus (.35, .4) | 0.226 |
| tall-ellipse-ring | 0.411 |
| lens-ring | 0.298 |
| arc-ring | 0.297 |
| moon-ring | 0.688 |
| heart-ring | 0.142 |
| washer | 0.072 (0.57 flat) |
| sphere, cylinder (controls) | 0 |

  In the 70-instance set, saddles occur only in revolutions of circle-like/curved profiles: pie-9 0.20,
  arc-10/11/13 0.12-0.14, vesica-22/23 0.11/0.04, ellipse-45..49 0.04-0.24, moon-55/59 0.32/0.29; every
  other instance has at most 4 of 20,000 points (finite-difference noise at creases), i.e. none.
  So the 70-set has 13 instances with real saddle regions; the saddle set above supplies WP14 with more.

## Generate
- `Rng` (xoshiro256**/SplitMix64), per-item streams (`Rng::for_item`): no shared RNG, no data races.
- `Perlin`: improved noise (fade, 12 gradients) written from the paper; the permutation is a shuffle
  seeded by `base % 256`, NOT Perlin's table, so fields differ from the `noise` library's (same distribution).
- Blobs: `F = |x|^2 - (1 + 2 alpha n)^2` on 64^3 over `1.2(1+alpha)`, our MC (dense), closed-manifold check, then
  **our QEM simplifier** (`simplify.rs`, replaces pyfqmr) to <= 2048 faces, `center_and_scale`. Persistence 0.5,
  alpha (0.1,0.7), scale (0.5,2), octaves 1..4, as the reference's final call. Multi-component blobs are accepted
  (as the reference) unless `require_single_component`.
- Scanner: `bvh.rs` (own BVH: `truth.rs` has no ray casting), 4 cameras Fibonacci x random rotation, radius sqrt 3,
  grid side 2 sqrt 3, spacing U(.01,.1), face normals facing the ray, FPS/random subsample. FPS and ray order are
  thread-count independent (tested).
- `dataset.rs`: `write_training_npz` (own stored-zip writer `npz.rs`, no ZIP64: files must stay < 4 GiB, i.e.
  <= ~500 shapes) in the authors' layout; neighbors from our kd-tree; unit normals (D23). Read back by WP7's
  `train::data::load_file` (test with `--features train`) and `npy::parse_npz_stored`.
- `d2.rs` for WP9b: `Curve2d` (polygon / traced noise blob, exact polyline SDF, O(segments) per query), random polygons,
  2D noise blobs, `simulate_scan_2d`, `Instance2d` (analytic / polygon / blob), `random_instance_2d`, ground truth,
  oriented clouds with exact normals, 2D FPS.

## SPEC proposals / open items
1. Mode-0 clouds and blob meshes need a simplifier; SPEC §17 does not list one. Added `generate::simplify`.
2. Blob `F` is not an SDF; `BlobField` implements `ImplicitField` with `Band::Dense` only.
3. Vesica `d < 0` shapes are discs unions (see above); if WP14 wants lens-only, restrict `d > 0`.
4. Our noise differs from the `noise` library's permutation; exact regeneration of the authors' blobs is impossible.
5. Not done: 500-shape dataset runs (not requested); simulated files with n < 2048 are padded by `write_training_npz`.
