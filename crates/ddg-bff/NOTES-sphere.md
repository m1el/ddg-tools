# Map to the sphere: working notes (checkpoint)

Branch `bff-sphere` in the clone
`/var/data/discrete-geometry/NavigatingIntrinsicTriangulations/navigating-intrinsic-triangulations-demo/.claude/worktrees/agent-af604adf9c2b21fae/ddg-tools-sphere`.
The clone was made from `/var/data/discrete-geometry/ddg-tools` because the assigned worktree was
of the wrong repository. Base: `75b6e2c`. To merge: `git fetch <clone path> bff-sphere`.

## Done (committed)

- `src/sphere.rs`: `map_to_sphere(&TriMesh, &Geometry, &SphereOptions) -> Result<SphereMap, BffError>`;
  `mobius_center`, `sphere_distortion`, `inverse_stereographic`, `stereographic`,
  `sphere_inversion`, `most_regular_star`; errors `NotASphere`, `PoleOutOfRange`.
  The method: remove the pole star (submesh), `DiskBff::flatten_to_disk(10)`, normalize the disk
  (boundary mean radius 1), ternary search of `ln R` in `[1, 1e4]` maximizing the area-weighted
  mean of `ln(a/A)` (reference objective), stereographic projection from the south pole, Möbius
  centering after Baden–Crane–Kazhdan 2018 Alg. 1 with *surface* areas, normalized centroids,
  tolerance 1e-10, step cap 0.5 and halving line search.
- `tests/sphere.rs` (9 tests, pass), `examples/sphere_report.rs` (table below).
- Differential test: `reference/bff_ref.cpp` mode `sphere`, `reference/run_sphere.sh`,
  `examples/reference_sphere_diff.rs`. Results (`runs_tmp/sphere-refdiff.txt`): our disk on the
  punctured mesh = reference's own disk mode to 2e-15..6e-13 rms; the reference's in-place
  puncture inside `mapToSphere` differs from that by 2e-7..1e-4 (unexplained reference artifact);
  after both are centered the same way the sphere maps agree to 3e-6 (icosphere) .. 8e-4 (bunny),
  5e-3 (axes), 5e-2 (brain, 38 flips in both). The reference's own centering (image areas) leaves
  |µ| 0.03–0.6 by the surface measure; its area distortion is higher (axes 2.41 vs 1.83).
- Lean `DdgProofs/Bff/Sphere.lean` (stereographic projection on the sphere / inverse / conformal /
  orientation; inversion keeps the sphere; derivative 2(u − (u·x)x)); CITED: centering existence,
  spherical triangle area, sphere injectivity. `lake build` OK, `check-formulas.sh --build` OK.
- ddg-py: `ddg_tools.map_to_sphere` (`src/sphere.rs`, `src/python/sphere.rs`), stub, README,
  `python-tests/test_sphere.py`; `scripts/py-test.sh` PASS (3 Pythons).
- Blender: `Conformal Map to Sphere` operator (`object.ddg_map_to_sphere`, point attribute
  "sphere", optional shape key), 2 headless tests; full `scripts/blender-test.sh` passes
  (27 headless, cut 49/48+1 skip, UI 39 checks).

## Measurements (single thread, run alone, best of 3, ms)

| mesh | V | QC mean / max | area rms u | flipped | ms | pole 0 (reference's): QC mean / max, flipped |
|---|---:|---|---:|---:|---:|---|
| icosphere 3 | 642 | 1.0096 / 1.240 | 0.006 | 0 | 4 | same |
| icosphere 5 | 10242 | 1.0008 / 1.244 | 0.002 | 0 | 80 | same |
| ellipsoid 2:1:0.5 noisy | 2562 | 1.0746 / 1.332 | 0.500 | 0 | 19 | 2.0658 / 28.0, 4 |
| spot | 2930 | 1.1332 / 1.780 | 1.022 | 0 | 18 | 1.3028 / 2.451, 0 |
| blub | 7106 | 1.1039 / 1.843 | 2.173 | 0 | 51 | 1.2364 / 2.088, 0 |
| bob | – | error: genus 1 | | | | |
| bunny | 14290 | 1.0410 / 1.289 | 1.631 | 0 | 121 | 1.0500 / 1.308, 0 |
| axes | 3735 | 1.1219 / 1.243 | 1.827 | 0 | 26 | |
| chair | 5726 | 1.0544 / 1.653 | 0.560 | 0 | 43 | 1.2777 / 1.960, 0 |
| brain | 152893 | 1.0529 / 71.9 | 0.459 | 38 | 1554 | 38 flips too |

QC max ≈ 1.24 on icospheres sits next to the removed star (not refined away).

## Remaining

1. Crate README `crates/ddg-bff/README.md`: section "Mapping to the sphere" (API row, method,
   the table above, the differential table), deviations list (pole rule; centering weights =
   surface areas + normalized centroids + tol 1e-10 + step control vs reference image areas,
   1e-3, 1000; scale-invariant radius search; disk normalization; south pole; submesh instead of
   in-place star; errors for non-genus-0 / disconnected), remove "mapping to the sphere" from item
   16 and from "Not done", formula tags table rows, "Other tests"/"Reproducing" lines; crate doc
   in `src/lib.rs` (mention `sphere`).
2. `ROADMAP.md`: add a "Landed" entry (it had no sphere entry).
3. `cargo test --workspace --release > runs_tmp/sphere-ws-test.log 2>&1` and
   `cargo clippy -p ddg-bff -p ddg-py --all-targets -- -D warnings` (only ddg-bff clippy was run;
   ddg-py with `--features python` not clippy-checked).
4. Delete this file once the README carries the content.

## Commands

```sh
export CARGO_TARGET_DIR=$PWD/target-sphere
cargo test --release -p ddg-bff --test sphere
cargo run --release -p ddg-bff --example sphere_report
crates/ddg-bff/reference/build.sh /tmp/.../bffref    # SuiteSparse can be copied from /tmp/ddg-bff-takeover/ref3/ss-install
CARGO_TARGET_DIR=... crates/ddg-bff/reference/run_sphere.sh /tmp/.../bffref
scripts/check-formulas.sh --build
CARGO_TARGET_DIR=$PWD/target-sphere-py scripts/py-test.sh
CARGO_TARGET_DIR=$PWD/target-sphere-py scripts/blender-test.sh --no-build
```
