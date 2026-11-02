# WP11R notes (Rust side of the Lean oracle test, `tests/oracle.rs`)

Run: `cd lean/oracle && lake build`, then
`DDG_RS_FUZZ_CASES=10000 DDG_RS_FUZZ_SEED=7 cargo test --release -p ddg-repulsive-surfaces --test oracle -- --ignored fuzz --nocapture`
(env: `DDG_RS_FUZZ_CASES/SEED/THREADS/DIR`, `DDG_ORACLE`). Tolerances and skips are in the module doc.

## Fixtures (`fixtures_are_reproduced`, always runs)

70 fixtures: 44 answered ok, 19 oracle errors (all matched to the `RsError` variant or the remesh
`RefusalKind`), 7 `step` pending. Per task: energy 14 (9 ok), df 4, metric 7, constraints 8, saddle 8,
precond 3, remesh-ops 16. Worst relative errors: `x` 2e-12, `gx` 2e-12 (ratio to tolerance <= 0.003).
Special readings asserted: `rs-error-singular` (Rust rank check drops the duplicate pin, D10),
`rs-mobius-flip` (Rust reorders and succeeds, SPEC §10).

## Fuzz (`rs-gen` -> `rs`, replayed on 64 threads)

| seed | cases | answered | oracle errors | step (pending WP5) | comparisons | disagreements |
|---|---|---|---|---|---|---|
| 1 | 2400 | 2062 | 48 | 290 | 7274 | 0 |
| 7 | 10000 | 8612 | 231 | 1157 | 30081 | 0 |

Max relative error (seed 7, normwise): energy 4e-15, gradient 3e-14, df 6e-16, A/B/B0/L/Lfr 2-3e-15,
W/M/WMW/Phi/Lp/eps exactly 0, J 7e-16, `g` 3e-14, `x` 1.4e-11 (worst ratio 0.14 of the 1e-10
tolerance), `S` 1e-12, `h` 3e-12, `gx` 7e-12, `z` 1.5e-14, remesh positions exactly equal; remesh
connectivity, parents and vertex maps equal in all 641 cases (including the reordered flips).
No tolerance was loosened; no disagreement between Rust and oracle was found.

## Pending

`step` (1157 of 10000 generated cases, 7 fixtures): `step_fixtures_are_reproduced` is `#[ignore]`d
until WP5 (`Session`) can run one dense projected-gradient step with the standard Armijo rule and the
fixture's start rule (`abs`, `maxdisp`, `reference`) and `newton fixed|tol N`.

## Refactor

`tests/wp2_fixture/mod.rs` gained text parsers (`parse_case`, `parse_output`, `blocks`, `manifest`) and
the `remesh-ops` fields; `tests/remesh.rs` now reads its fixtures through it (no behaviour change).
