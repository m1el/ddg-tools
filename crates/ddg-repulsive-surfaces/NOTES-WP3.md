# WP3 notes (Barnes–Hut energy and its exact gradient)

Build: `CARGO_TARGET_DIR=$PWD/target-wp3 cargo test -p ddg-repulsive-surfaces --test barnes_hut`
(release recommended: `--release`; the 40k timing is `-- --ignored --nocapture`, run alone).

## Done
- `src/spatial/mod.rs`: `FaceElement` (`Element<3>`: barycenter, triangle vertex box, area),
  `Cluster` (A_I, X_I, P̄_I as `[xx,xy,xz,yy,yz,zz]`, vertex AABB, radius to the farthest box
  corner), `FaceTree` (`Bvh<3>` + clusters + leaf bookkeeping; `build`, `refit`,
  `refit_faces`), `box_distance_sq` (rs-box-distance), `farthest_corner_distance`,
  `bh_admissible` (strict, squared form as the reference).
- `src/energy/barnes_hut.rs`: `BhTree` (`build(pool, surface, positions, opts)`,
  `build_with_geometry`, `refit(positions)`, `energy(pool, e)`, `energy_and_differential(pool,
  e) -> Result<EnergyEvaluation>`, `face_partials`, `lists()`, `face_tree()`, `geometry()`),
  `InteractionLists` (CSR per leaf + transposes), `build_lists`, the WP1 hook (unchanged
  signature).
- `tests/barnes_hut.rs`: θ = 0 vs exact, FD with frozen lists, list partition/admissibility,
  cluster data, refit vs naive evaluation of Ẽ with the same lists, thread determinism, errors,
  θ decay table, 40k timing (ignored).

## Decisions
- Leaves can be lumped (reference does it, SPEC §3.4 "node C").
- `FaceTree::refit` calls `Bvh::refit` (so the BVH's own node boxes, masses and moments are
  current), but `Cluster::bbox`/`radius` stay those of the build: the admissibility used by the
  lists is frozen. `A_I`, `X_I` are the BVH node mass / centroid; `P̄_I` from `fold_up`.
- `BhTree::refit` errors only on count mismatch / non-finite positions; a degenerate face gives
  a non-finite energy (rejected trial, SPEC §2.3).
- Gradient: per leaf (parallel) own partials, node terms `(b, u)` per lumped entry, partner
  partials per exact block (|C| entries; memory ~ 32 B × Σ|C| over exact blocks); node sums and
  partner sums over transposed lists (increasing entry); downward pass re-centred as the curves
  crate's `rc-bh-adjoint-shift`; WP1's `face_pull_back`. Bitwise independent of threads.
- `energy()` does exactly the value operations of `energy_and_differential()` (bitwise equal).
- `NonFiniteEnergy` culprit: first non-finite exact pair of the lists, else largest exact term.

## Results (release, 2026-11-01)
- θ = 0 vs exact (6 meshes × 2 exponent pairs × leaf 1/3/8): energy rel 2.0e-15, gradient rel
  3.1e-15.
- FD (h = 1e-6 ℓ̄, all coordinates, frozen lists), θ = 0.5 (6,12) and θ = 1 (4.5,9.5):
  icosphere 4.5e-9 / 1.0e-8, torus 3.6e-9 / 7.2e-9, disk 9.1e-10 / 1.6e-9, Möbius 2.0e-10 /
  3.5e-10, two components 2.3e-9 / 4.7e-9, reference Klein bottle 2.0e-8 / 3.3e-8.
- Error vs θ, 10k-face torus (R 1.5, r 0.6, 0.02 ℓ̄ noise), (6, 12), leaf 8:

  | θ | energy rel err | gradient rel err | lumped | exact pairs |
  |---|---|---|---|---|
  | 1 | 8.4e-2 | 3.5e-2 | 169628 | 1.48e6 |
  | 0.5 | 2.3e-2 | 9.2e-3 | 363937 | 3.23e6 |
  | 0.25 | 4.1e-3 | 2.5e-3 | 781290 | 9.26e6 |
  | 0.125 | 2.1e-4 | 2.0e-4 | 1158331 | 3.25e7 |
  | 0.0625 | 1.3e-6 | 3.4e-6 | 651638 | 8.12e7 |

  Rates per halving (energy): 1.85, 2.52, 4.29, 7.37 — at least O(θ²).
- Timing: see below.

## Remaining
Nothing planned. Possible speedups (not needed so far): share `r^{-β}` between the two
directions of exact blocks present in both lists; SIMD-friendly SoA leaf data.

## For WP6 (block cluster tree on the same BVH)
- `bh.face_tree()` (or `FaceTree::build(surface, positions, &geom, leaf_size)`) gives
  `bvh(): &Bvh<3>` (the curves crate's; `BlockClusterTree::build` takes it) and
  `clusters()[n]` with `area`, `centroid`, `projector` (P̄, upper triangle; 
  `projector_quadratic_form(v)` = vᵀP̄v), `bbox`, `radius` (frozen at build). Separation for the
  BCT is `max(r_I, r_J)² ≤ χ² box_distance_sq(...)` (`≤`, unlike BH's `<`).
- `FaceTree::refit` refreshes `area`, `centroid`, `projector` only.
