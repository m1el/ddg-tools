# WP16B notes (method comparison: Fig. 5 and Fig. 12)

Reproduce: `cargo build --release -p ddg-repulsive-surfaces --example compare_methods` (with `CARGO_TARGET_DIR=$PWD/target-wp16b`),
`scripts/rs-eval.sh gen && scripts/rs-eval.sh run` (56 runs, 1 thread each, ~90 min), `uv run --with matplotlib scripts/rs_eval_table.py`
(tables and `fig12.png`, `fig5.png` in `runs_tmp/wp16b/`, not committed). Per-run CSV + log: `runs_tmp/wp16b/{fig12_<mesh>,fig5_<faces>}_<method>.*`.

## Setup (as run)
- Methods: Hs, L2 (mass), H1, H2 projected gradient with barycenter + area constraints; AQP (kappa 100), H1 L-BFGS (memory 20),
  BCQN (memory 20, Standard) on `soft_problem` (area penalty weight 100). Defaults otherwise (remeshing, CCD on, Barnes-Hut theta 0.5), `StopRule::Off`.
- Fig. 12: `generate::fig12(0..5)` (CLI `fig12:1..5`), 25 076 / 26 348 / 21 640 / 44 848 / 23 344 faces. Budget 2400 s for meshes 1-3
  and 3600 s for the two tunnels (4, 5), wall clock of the steps only (exact-energy checkpoints excluded). Exact energy (1 thread) at the first step
  after 60/300/1200 s (in the CSV) and at the end; table values at fixed times are Barnes-Hut energies, "end (exact)" is exact.
- Fig. 5: `genus_surface(2, n)` at 10 952 / 39 844 / 134 832 faces, 300 steps, wall cap 5400 s per run; energy is Barnes-Hut (exact at steps
  100, 300 for 11k/40k is in the CSV `E_exact` column).
- **Normalisation.** The unweighted energy scales like 1/area (alpha, beta = 6, 12), and the soft methods let the area grow (up to 1.1-2.2x on
  the Fig. 12 meshes), which lowers the raw E. All values are reported as `E * A / A0`. This assumes exact scale covariance; it is
  not an independent check, and the soft methods still solve a slightly different problem (area is not held).
- **Machine load (fact).** `uptime` at start: load average 79.05 / 64.52 / 36.69 (before the runs started; 56 of ours were added);
  at end 51.23 / 53.53 / 59.83. The 64-core machine was oversubscribed by processes outside this work, and 56 runs shared memory
  bandwidth, so absolute times are inflated by an unknown factor (not measured); Hs (memory-bound hierarchical products) is probably
  affected more than the sparse methods.
- No run ended Stuck or Error. All Fig. 12 runs and the 135k Fig. 5 runs (except L2) ended on the wall budget.

## Fig. 12: E*A/A0 (exact at end)
| mesh | Hs 60/300/1200/end | best other at end | steps (Hs / others) |
|---|---|---|---|
| 1 linked | 126 / 50.9 / 14.4 / 7.7 | 26.0 (H1-LBFGS, BCQN) | 261 / 490-716 |
| 2 doubly-wound | 167 / 134 / 59.6 / 7.2 | 60.0 (H1-LBFGS) | 217 / 448-626 |
| 3 6x-wound | 2853 / 1308 / 1197 / 902 | 1179 (BCQN), 1210 (H1); AQP diverged (area x1.34, E up to 6e6 at 300 s) | 134 / 543-772 |
| 4 figure-8 tunnel | 1426 / 561 / 425 / 353 | 341 (BCQN), 367 (H1-LBFGS), 391 (H1) | 100 / 285-515 |
| 5 trefoil tunnel | 566 / 201 / 164 / 158 | 154 (H1-LBFGS, BCQN), 156 (AQP), 158 (H1) | 244 / 802-1198 |
Full per-method tables: `runs_tmp/wp16b/tables.md` (regenerate with the script). Hs GMRES iterations per solve: 9.5 / 12.1 / 30.1 / 28.9 / 20.5
(meshes 1-5); the sparse methods need 1-2.

## Fig. 5: E*A/A0 after n steps (n = 10 / 100)
| method | 11k | 40k | 135k |
|---|---|---|---|
| Hs | 64.3 / 41.1 | 42.3 / 39.5 | 42.6 / 39.1 |
| L2 | 328 / 139 | 406 / 227 | 437 / 365 |
| H1 | 116 / 49.7 | 185 / 49.8 | 141 / 49.3 |
| H2 | 220 / 42.0 | 50.5 / 40.2 | 52.4 / 39.8 |
| AQP | 96 / 49.8 | 287 / 49.5 | 223 / 49.7 |
| H1-LBFGS | 118 / 44.4 | 125 / 45.5 | 237 / 50.4 |
| BCQN | 111 / 45.7 | 276 / 48.4 | 236 / 52.2 |
Seconds per step 11k -> 40k -> 135k: Hs 3.5 / 14.2 / 54; H1 1.3 / 6.2 / 21.6. Steps reached by the 5400 s cap at 135k: Hs 100, H1 251, H2 169, AQP 276,
H1-LBFGS 172, BCQN 174 (L2 300).

## Comparison with the paper's claims
- *Hs reaches lower energy sooner on every mesh*: holds on meshes 1 and 2 (Hs best from about 300 s on, a factor 2-8 lower at the end), and on
  mesh 3 at the end (902 vs 1179) and from 1200 s; **not on meshes 4 and 5**: Hs is best only at its end value against the projected methods
  (L2, H2 far behind); BCQN/H1-LBFGS are at equal or lower energy at 300 s and at the end (341 vs 353; 154 vs 158), and H1/AQP/L-BFGS
  win at 60 s on mesh 4 and 5. Hs took 36 s/step on the 54k-face mesh 4 (100 steps in 3600 s) against 7-8 s for H1-type steps.
  Per iteration Hs is the strongest on every Fig. 12 mesh (lowest energy at equal step count); per second it is not, in our single-thread, loaded run.
- *Per-iteration progress of Hs is resolution independent while the others slow down*: Hs at step 100: 41.1 / 39.5 / 39.1 (flat). L2 slows
  strongly (139 -> 227 -> 365). The H1 family at step 100 is also nearly flat at the energy plateau (H1 49.7 / 49.8 / 49.3) on this surface, and
  H2 is flat from 40k on; L-BFGS/BCQN degrade mildly (44 -> 50, 46 -> 52). At step 10 the others are noisy (initial remeshing transient).
  So L2 and, mildly, the quasi-Newton methods reproduce the claim; H1/H2 do not in this experiment. Unresolved here: the genus-2 slab
  is a smooth, easy surface (all methods reach 39-50 within 100 steps), and the 135k runs are capped before 300 steps, so the step-300 column is
  missing for 5 of 7 methods.
- The Hs step is 2.7-3x (per step) more costly than H1-type steps at equal resolution and 10+ GMRES iterations/solve (up to 30 on the 6x handcuffs).

## Caveats
- Single seed/mesh analogues (D31); meshes differ from the paper's. Wall times under external load as above; no repeat to measure run-to-run noise.
- Remeshing changes the face counts differently per method (e.g. mesh 3: 23k for Hs vs 30k for the L-BFGS/BCQN runs).
- AQP on mesh 3 diverged in energy (peaks 6e6) and recovered only partially (E 1.9e5 at the end); not investigated (soft area penalty, kappa 100).

## Proposed SPEC §14.2 "Results of record"
> Results of record (WP16B, 1 thread per run, 56 runs in parallel on a loaded 64-core machine; energy `E*A/A0`, exact at the end).
> Fig. 12 end values Hs vs best other: 7.7 vs 26.0 (linked), 7.2 vs 60.0 (doubly wound), 902 vs 1179 (6x), 353 vs 341 (figure-8), 158 vs 154 (trefoil).
> Hs is lowest per iteration on all five meshes and lowest per second on the first three; on the two tunnels the H1 L-BFGS / BCQN methods are
> as good or slightly better per second (Hs 20-36 s/step on 24-54k faces, 20-30 GMRES iterations/solve on meshes 3-5). Fig. 5: Hs progress at step 100 is
> resolution independent (41.1 / 39.5 / 39.1 at 11k / 40k / 135k faces); L2 slows by 2.6x; H1, H2, AQP are flat on this surface, L-BFGS/BCQN slow mildly.
> The paper's wall-clock claim is therefore only partly reproduced; see NOTES-WP16B.md for conditions.
