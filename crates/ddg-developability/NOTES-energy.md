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

## Remains
- (see final report)
