# WP15 and WP4b: oracle replay, determinism audit, formula tags, API items, crate-owned `exp`

Logs: `runs_tmp/pat-wp15/`. Everything below was run with `CARGO_TARGET_DIR=target-pat15`, `nice -n 19`.

## 1. Oracle replay (`tests/oracle.rs`)

All 35 fixtures of `lean/oracle/fixtures/pat/` (71,426 tokens) are replayed through the library:
features (`neighborhood_features_f64`, `fitting_neighborhoods`, `scale_outputs_f64`), fit
(`fit_torus`, f64 coefficients as in the fixtures), SDF (`torus_sdf_laplacian`), blend
(`TorusField::blend_at`, every neighborhood and λ rule), SHC/SSPD (`KernelField`), and the error cases.
Floats at 1e-12 relative to `max(|x|, 1)`, integers and tags exactly, neighbor sets exactly (as sets
where Rust sums in index order).

First run: 11 of 35 disagreed. Causes, all found, none a tolerance change:

* **Test harness** (6 fixtures): the fixture's Rayleigh pair is `(κ₊, κ₋)` in eigenvector order, the library
  reports it as `(κmin, κmax)` (compared as an unordered pair); a one-point cloud is not a `PointCloud`
  (n ≥ 2), so the unit-normal check and the single-point blend run on the raw data / `blend()`.
* **Operation order** (library, 5 fixtures, Δφ and the gradients in the last bits): the oracle computes the
  torus gradient as `((x−c−zu)/v)·w₁ + w₂u`, the library `w₁·(x−c−zu)/v + …`; the oracle's
  `r̂·(λg)` against the library's `(λ r̂)·g`; and the Laplacian as `(ΔN − φΔW)/W + 2B·∇φ`, the library as
  `A·B + ∇φ·B + φ(…)`. These are algebraically equal, but Δφ sums terms of size λ² (10³–10⁴) to a
  result of order 1, so the split form lost up to 1e-9. **Changed**: `sdf.rs`, `blend.rs` and
  `gpu/kernels.cu` now follow the oracle's order (the Laplacian is also cheaper: two sums fewer).
  Golden parity (`tests/parity.rs`) still passes.
* **Index-order sums** (`All`, kNN with k ≥ n): the README's note. The harness replays the library's `blend`
  on the oracle's `(ρ², index)`-sorted set and checks that the index-order evaluation agrees with it to
  `1e-10·λ²`.
* **Fit round-off in Δφ** (remaining, 3 fixtures): the blend fixtures that give `coef` lines fit their tori
  in the oracle and in the library with different (equally valid) operation orders inside the fit, so
  the tori differ by 1e-16 relative, and Δφ amplifies that by λ². Observed 0.7–1.2 ε λ² absolute, e.g.
  1.5e-12 at λ = 76. The Δφ token alone is compared at `1e-12 + 16 ε λ²`; φ, ∇φ (frozen, exact), σ, λ,
  shift and all counters stay at 1e-12 / exact. Lean's `Float.exp` and glibc's `exp` agree bit for bit
  on 10⁵ arguments (checked), so `exp` is not the cause.

## 2. Determinism audit (`tests/determinism_audit.rs`, plus `determinism.rs`)

Each result is compared through its `Debug` text (shortest round-trip digits: equal text is equal bits)
at 1 thread twice, then 3, 8, 32: build + evaluation (value, gradient, Laplacian, the tori), `sample_level_set`,
`marching_cubes`, Booleans/offsets of two fitted fields + MC, `morphology`, `orient_normals`,
`pca_normals_random_sign`, point areas, SHC/SSPD/LSE, winding numbers, FWN, `MeshTruth` (exact and fast sign),
`generate_blob`, `sample_cloud`, `generate_blob_records`, `simulate_scan`, mesh sampling, FPS, simplification,
2D generators, analytic sampling. The existing 2D (`pat2d.rs`) and MC tests were moved to 1/3/8/32.
All pass. (Nothing found.) GPU entry points are covered by `tests/gpu.rs` (run to run) when the GPU is free.

## 3. Formulas

`scripts/check-formulas.sh`: OK, 0 missing. Added CITED entries `pat-truth-closest-point`, `pat-solid-angle`,
`pat-triangle-sampling`, `pat-query-sampling` (the last two were used by `src/train/queries.rs`), tags on
`truth::closest_point_on_triangle` and `baselines::fwn::triangle_winding` (it was tagged `pat-winding-dipole`,
which is the point-element dipole; the triangle element is the solid angle).
Lean (`DDG/SDF/Torus.lean`, builds without Mathlib recompiling): `torus_sdf_stable_general`
(any expansion point `o`, no `(o−c)·u = 0`: `v − v_o = (v² − v_o²)/(v + v_o)` with
`v² − v_o² = |δ|² + 2δ·m − (δ·u)² − 2(δ·u)(m·u)`), `stable_core_dist` (`√((d₁+e₁)²+(d₂+e₂)²) − |r| =
(K + 2(d₁e₁+d₂e₂) + e₁²+e₂²)/(√… + |r|)`, `K = d₁²+d₂²−r²`) and `torusSDF_stable_general` (the SDF in the doubly
stable form of the F32Stable kernel). The old `torus_sdf_stable` stays (it is the special case used by the
fitted torus's constant).

## 4. API items (SPEC §18, WP12/WP13)

* `EigenRule`, `GpuMatmul`, `Progress` re-exported at the crate root; `Morph::levels` is public.
* `morphology` returns `(TorusField, PointCloud, MorphReport)`. **Breaking** for the three callers (tests, ddg-py).
* `TorusField::{cloud, cloud_hash, sign_rule, with_sign_rule, to_bytes, totals}`; the field now keeps its
  `PointCloud` (normals included).
* Tori file: the reserved `u32` at offset 12 is the sign-rule marker (0 not recorded, 1 Reference,
  2 Geometric). `Tori::to_bytes_with` / `from_bytes_with`; old files and the old functions are unchanged
  (old readers ignore the word). An unknown marker is a `Format` error. `PointsAsTori::build` records the rule.
* `reconstruct` fills an eval report: `ImplicitField::eval_totals` (default `None`; `TorusField` and `Csg`
  implement it) and `McReport::eval` = difference of two snapshots; batch evaluations only (single-point
  queries, used by the refinement, are not counted). Python: `eval_report` after `reconstruct`.
* `smooth = 0` is the sharp Boolean (`Csg::Smooth*` with `k = 0`); negative/non-finite is still an error.
* Progress/cancel: `Progress` (atomic counter + flag, clones share state) in `BuildOptions::progress`
  (stages `neighborhoods`, `network`, `fit`, `tree`) and `EvalOptions::progress` (`evaluate`);
  `PatError::Cancelled`; checked per CHUNK of the network, per 512 queries. Does not change results (tested
  bitwise). A CUDA network run is checked before it starts only. Orientation, level-set sampling and
  the other native calls have no hook yet.
* `ddg-py`: `PatProgress`, `progress=`, `sign_rule` property, tori files marked, `reconstruct` sets
  `eval_report`, `with_options` never inherits a progress handle. Blender: Esc cancels the build and the MC
  evaluations; `smooth` is passed through; new test `cancel_stops_the_native_build_and_caches_nothing`.

## 5. WP4b: `exp` (`src/exp.rs`, `gpu/kernels.cu: pat_exp`)

Cody–Waite reduction (`ln2_hi` with 21 trailing zero bits), degree-13 Taylor polynomial in Horner form on
`|r| ≤ ½ ln 2`, `1 + (r + r²P)`, scaling by `2^k` by bits. No FMA anywhere; the CUDA text is the same
sequence, compiled with `--fmad=false` (NVRTC compile of the kernel file checked: OK).

* Against glibc over 10⁸ uniform arguments in `[−800, 0]` and in `[−745, 709.78]`: max 1 ulp difference,
  8.8 % / 9.6 % of arguments differ. Against mpmath (200 bits) on 3·10⁵ arguments: **max error 0.94 ulp**
  (glibc ≈ 0.5). Edge cases (±0, ±∞, NaN, overflow at 709.78, underflow at −745.13, subnormals) match glibc.
* Cost, single thread, 10⁶ blend-range arguments: libm 4.9 ns, crate 11.0 ns per `exp` (machine at load ≈ 13).
  Blend (`examples/blend_bench.rs`, 20k points, 200k queries, k = 32, one thread, load 13): distance only
  5.0–5.5 µs/query before, 5.5–6.7 after (≈ +10–15 %); distance + gradient + Laplacian 5.7–6.0 before, 6.1–7.6
  after (≈ +5–10 %). Noisy (shared machine); the exp itself is +6 ns of 32 terms ≈ 0.2 µs/query.
* CPU blend bits changed by round-off; `tests/parity.rs` (golden), `tests/oracle.rs`, `tests/sdf.rs` pass.
* **GPU bitwise test: written, not run** (the training process holds the GPU: `nvidia-smi` showed
  `target-pat10/release/examples/train`). Run later: `cargo test --release --features cuda -p ddg-points-as-tori
  --test gpu pat_exp_bitwise_host_device blend_bitwise_cpu_gpu` (needs the GPU free). Until then the blend is
  level 2 by construction (identical operation sequences, `--fmad=false`), not by measurement. The old
  `blend_matches_cpu_f64` tolerance test stays. NOTES-GPU.md's "level 1 blend" wording is superseded once the test passes.

## 6. Clippy

`cargo clippy -p ddg-points-as-tori --lib --tests -- -D warnings` clean for default, `--features train` and
`--features cuda` (and this WP's examples). **Not clean, not mine**: WP14's `examples/eval.rs`, `examples/density.rs`,
`examples/eval_util/util.rs` (`neg_cmp_op_on_partial_ord`, `type_complexity`, `needless_range_loop`,
`collapsible_if`); `cargo clippy --all-targets` fails on them until WP14 fixes them.

## SPEC changes proposed

1. §7.3 / §18 WP4 item 4: replace "needs generalizing" by the three new Lean theorems (section 3).
2. §18 WP4 item 3: "the blend is level 2" after the GPU test above has run; record the exp accuracy (0.94 ulp).
3. §7.4: the Laplacian is `(ΔN − φΔW)/W + 2B·∇φ` (the form the code and the oracle use).
4. §15.3: Δφ in oracle replays carries a `λ²`-scaled round-off allowance (section 1).
