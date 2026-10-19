# Energies: progress notes (SPEC §1–§3, §7.1 evaluation, §10.2–§10.4, §10.8 energy bench)

## Done
- `src/error.rs`: `DevError` as in SPEC §7.1 (Display names element and reason).
- `src/geometry.rs`: `face_frame`, `corner_angle(s)`, `area_gradient` (Eq. 7), `normal_jacobian`
  (Eq. 8, D2), `angle_gradients` (Eq. 9, D3), `vector_area`, `backprop_face` (§2 adjoint).
- `src/energy/`: `Energy`, `EnergyOptions` (+`Default`, `validate`), `Evaluation`, `evaluate`,
  `evaluate_energy`, `star_energy`, `domain`; star functions `combinatorial`, `combinatorial_max`,
  `covariance`, `covariance_max`, `intrinsic`, `intrinsic_max`, log map + adjoint.
- Pipeline §3.9: faces → stars (fixed 256-vertex chunks) → per-face backprop → per-vertex gather
  in star order. Bitwise identical for 1/2/7/64 threads (test).
- Tests: `energy_fd` (FD on random interior stars k=3..9, boundary fans k=1..7, icosphere, bumpy
  grid, torus, crumpled, spot control), `energy_values` (§10.3), `energy_invariance` (§10.4:
  rigid/scale, Σg, Σf×g, Σf·g, retessellation, determinism), `energy_errors`, `energy_blocks`
  (Eq. 7–9, backprop, log map adjoint, circumcenter, altitude).
- Bench `benches/energy.rs`.

## Clarifications (reported)
- Margins: `π_max` uses min(partition gap, active-pair gap within the chosen partition);
  intrinsic width divides the eigen-gap by `Σθ` (= tr A_v), not tr Ã (scale-free noise when
  nearly flat); boundary stars with ≤ 2 faces get `+∞` (energy identically 0).
- `H(φ)` series used below `1e-2` (three terms) instead of `1e-3` (closed form loses ~ε/φ²).
- `threads: None` runs on the current rayon pool (caller's or rayon's default); `Some(n)` builds
  a local pool unless already inside an `n`-thread pool.
- `dev-regular-cone-values` not tagged in Rust (tests compare closed forms; Lean side decides).

## Bench (machine loaded by fuzzers, load 37–88 on 64 cores: indicative only)
Single thread, ns/vertex (value + gradient): spot E_P 1530, π_max 1310, E_λ 1020, λ_max 3850,
intrinsic 1340, intrinsic max 1820 (load ~58); bob/blub similar or lower (790–1290, λ_max 3.9–4.1 µs)
at load ~48. 16 threads at load ~48: bob/blub 290–620 ns/vertex. Multi-thread scaling not measurable.

## Remains
- Re-run `cargo bench -p ddg-developability --bench energy` on an idle machine (1T and 64T).
- The star layout (`Stars`) is rebuilt on every call; an evaluator that keeps it across optimizer
  iterations would save part of each call (unmeasured; not in the SPEC API).
- Lean ids (`check-formulas.sh` MISSING) wait for the Lean agent; all Rust ids match SPEC §9.
- 4-1 subdivision invariance is tested by the optimizer agent (`tests/subdivide.rs`).

## Two-face boundary stars (D33)

`star_energy` zeroes value, margin-tie, `y`, `s`, `z` of a boundary star with at most two faces after
computing it (so the intrinsic energies still report D12 undefined states). Before, `E_λ` took an
arbitrary eigenvector of the rank-deficient matrix and fed sensitivities of the size of the normals'
angle into the gradient (found by the oracle fuzz; `tests/robust.rs`).
