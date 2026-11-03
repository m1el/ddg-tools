# GPU-2 notes (the GMRES preconditioner on the GPU)

Status: stage 1 done (2026-11-03); stage 2 not done (the numbers below say it would not pay); GMRES
vectors stay on the CPU (below). Sections: done / remaining, design, determinism, tests, measurements,
regression set, what is left on the CPU, stage 2 and GMRES decisions, proposed SPEC text, how to resume.

## Done / remaining

* **Done.**
  * `src/solve/ldlt.rs`: `SimplicialLdlt` — faer's sparse `LDLᵀ` of `L′` (same elimination order as
    before: AMD of the Laplacian block, then the multipliers; refined by a postorder of the elimination
    tree), read out as a simplicial factor (`L` by rows and by columns, `D`; relaxed supernodes' explicit
    zeros kept). Triangular solves in a **fixed chunked order** (below). `LaplacianSaddle` uses it by
    default (`TriangularSolve::Chunked`; `TriangularSolve::Faer` keeps the old faer solves, option
    `PreconditionerOptions::triangular`). The chunked CPU solve is *faster* than faer's (45k: 2.3 vs
    3.3 ms for 3 columns; 135k: 7.9 vs 16.8 ms), so it is the CPU default and determinism level 2 holds.
  * The GPU schedule `LdltSchedule` (leaf subtrees, chains, levels, chunk slots) and its CPU emulation
    `solve_scheduled` (bitwise `solve_rows`; unit tests on a 1-D saddle and a 2-D grid with dense rows).
  * CUDA (`src/gpu/kernels.cu` `ldl_*`, `prec_aug`, `vec_*`; host `src/gpu/cuda_precond.rs`,
    `gpu::GpuPreconditioner`): the whole fine preconditioner `P⁻¹` on the device — gather, `L′` solve,
    vertex-pair fractional product on device vectors (GPU-1's tree, `apply_device`), augmented middle
    factor, second `L′` solve, `z += y`, scatter; one upload and one download per application. The coarse
    correction (several components, global barycenter) stays on the host around it.
    `PreconditionerOptions::gpu_solves` (default `true`) switches it; with `refine = true`, face pairs or
    `TriangularSolve::Faer` the GPU-1 split is used.
  * Bench flags `--no-gpu-solves`, `--faer` (`benches/step.rs`); `RS_DEVICE=L1|L2`,
    `RS_TRIANGULAR=faer` for `examples/apps`; `examples/gpu2_ldlt.rs` (factor structure, schedule, CPU
    solve times), `examples/gpu2_precond.rs` (one application: CPU, GPU-1 split, GPU-2; `L′` solve
    profile; vertex-pair product alone).
* **Not done (decisions with numbers below):** stage 2 (FSAI / AMG); GMRES vectors on the GPU.
* **Remaining / next:** the vertex-pair product kernel (GPU-1's `bct_rows`; now ≈ 40 % of a GPU
  application, latency bound: 2.0 ms for 1 vector, 3.6 ms for 8 at 45k), CUDA graphs for the ≈ 100
  launches of an `L′` solve, then GMRES vectors (needs the CPU dot order and the constraint-row products
  changed to GPU-friendly fixed orders, below).

## Design

**Operation order (the contract between CPU and GPU).** `X` (`n × m`, row-major in the factor's
numbering; column `c` = one axis of one vector) is solved by

```text
forward    x_j = b_j − Σ_{k ∈ row j of L} L_jk x_k          (terms k increasing)
diagonal   x_j = x_j / d_j
backward   x_j = x_j − Σ_{i ∈ column j of L} L_ij x_i        (terms i decreasing)
```

where every sum is cut into consecutive chunks of 32 terms (`ldlt::CHUNK`), each chunk summed from
`0.0` in order (`p = 0; p += l·x`), the chunk sums subtracted in order (`s = b; s −= p`). Every value is a
fixed function of the factor and the right-hand side whatever the schedule, so any schedule that
respects the dependencies and these orders reproduces it bit for bit. (A plain sequential sum would
force one thread per row of the multiplier rows, which are dense over a component: 22k–67k dependent
additions.)

**GPU schedule** (`SimplicialLdlt::schedule`). Postorder makes every subtree an index range and a node
with one child have it just before itself.
* *Leaf subtrees*: maximal subtrees of ≤ 32 nodes (`LEAF_MAX`), one thread per (subtree, column) in index
  order (`ldl_leaves`), before the chains (forward) or after them (backward). 45k: 756 subtrees with
  15,729 of the 22,425 nodes.
* *Chains* of the remaining tree (`j, j+1, …` with `parent(j) = j+1` and no other non-leaf child; ≤ 1024
  long), grouped into levels (forward: by height; backward: by depth). Per level, `ldl_partials` sums the
  chunks of terms from earlier levels (one thread per (chunk, column), all SMs), including the partial
  chunk that continues into the chain; then `ldl_chains` runs every chain with one block per (chain,
  column), one row per thread: start from the partials, then for each chain node in order one
  `__syncthreads`, and every later row of the chain adds its in-chain term (two terms prefetched in
  registers). 45k: 623 chains in 24 forward/backward levels; the sum of the longest chain per level (the
  critical path) is 1,187 nodes. The schedule is validated (every entry of `L` respects the levels) or
  refused (the GPU path then stays off).
* *Preconditioner pipeline* (`DevicePreconditioner::apply_many`): exactly the CPU's `apply_fine`
  (`v = s·L y − K y` per `(j, axis)` over `K`'s column in order, `z = L′⁻¹ v + y`).

## Determinism

* Gathers only, one thread per output in the fixed orders above; the only synchronization is
  `__syncthreads` inside a chain block (no atomics, no reductions whose shape depends on the launch).
* **Level 1** (`RunToRun`): bitwise reproducible (test). Against the CPU: FMA contraction, preconditioner
  outputs within `2e-7` relative (the `ε = 1e-10` near-null space of `L′` amplifies round-off), flows with
  the same GMRES iterations to ±1–4 per step (table below).
* **Level 2** (`MatchCpu`): the CPU default runs the same chunked solves, so the preconditioner, and with
  it every flow, is the CPU's bit for bit (tests; the 45k and 135k benches and the whole regression set
  below: identical iterations and energies every step).
* The CPU default changed by round-off: chunked solves instead of faer's (and the factor's postordered
  elimination order). Iterations against the faer solves: identical on 9 of 14 Fig. 13 cases, ±0.2–0.6
  per step on the others, same statuses, no fallbacks; Klein bottle (self-intersecting, chaotic) drifts
  after step 8 (table).

## Tests

* `solve::ldlt` unit tests: the chunked solve solves the saddle; the schedule's CPU emulation is bitwise
  `solve_rows` for leaf sizes 0, 8, 32 and 1, 3, 6 columns (dense constraint rows included).
* `tests/gpu.rs` `preconditioner_matches_cpu`: torus, torus with 12 pins, two tori (per-component and
  global barycenter) × augmented / congruence × batches 1, 2, 4: level 2 bitwise the CPU, level 1
  reproducible and within `1e-6` (measured `2.0e-7`), a vector's result independent of the batch.
* `tests/gpu.rs` `flow_trajectories_match_cpu` (GPU-1's) now runs with the GPU solves: level 2 equal to
  the CPU bit for bit (energies, GMRES iterations, faces, positions), level 1 reproducible.
* All crate tests pass with and without `--features cuda`; clippy clean on the touched files.

## Measurements (RTX 4080, Threadripper 3970X; logs `runs_tmp/gpu2/`)

Outside load during the timed passes: another session's jobs kept the 1-minute load average at 2.4–10
(recorded at the start and end of every run in the logs); the CPU-at-4-threads baselines agree with
GPU-1's quiet-machine numbers to a few percent.

**One preconditioner application** (`examples/gpu2_precond.rs`, median of 11; 4 CPU threads; GPU times
include the upload and download):

| | 1 vector | 2 vectors | 4 vectors |
|---|---|---|---|
| 45k, CPU chunked / CPU faer | 8.5 / 10.7 ms | 9.9 / 12.6 | 12.6 / 15.2 |
| 45k, GPU-1 split (fractional on GPU, faer solves on CPU) | 7.6 | 9.1 | 11.2 |
| 45k, GPU-2 L1 / L2 | 5.5 / 6.0 | 5.2 / 5.7 | 5.6 / 6.4 |
| 45k, one `L′` solve on the GPU (L1) | 1.1 ms (GPU-2 first version 1.8, before leaf subtrees and finished partial chunks) | | |
| 135k, CPU chunked / CPU faer | 26.3 / 41.1 | 34.6 / 50.2 | 43.3 / 59.7 |
| 135k, GPU-1 split | 20.5 | 27.9 | 36.1 |
| 135k, GPU-2 L1 / L2 | 7.7 / 8.3 | 7.9 / 8.6 | 9.5 / 10.9 |
| 135k, one `L′` solve (L1) | 3.0 ms | 2.7 | 3.4 |

Of a 45k GPU application (≈ 5.5 ms): two `L′` solves ≈ 2.2 ms (critical path through the elimination
tree, latency bound), vertex-pair product ≈ 2 ms (GPU-1's kernel), the rest transfers and small kernels.

**Whole steps** (`benches/step.rs`, warm-up 2, mean of 5, 4 threads unless said; `timedB/`):

| 45k figure-8 tunnel (`fig12:3`) | step | precond | solve | factor | GMRES its/step | F after 7 steps |
|---|---|---|---|---|---|---|
| CPU 4 threads (chunked) | 2827 / 2826 ms | 418 | 1525 | 412 | 59.2 | 560.282535507900548 |
| CPU 4 threads, faer solves (= before GPU-2) | 2892 / 2964 | 548 / 563 | 1644 | 399 | 59.6 | 561.208135632529547 |
| CPU 16 threads | 1968 / 1949 | 393 / 385 | 1140 | 286 | 59.2 | = CPU 4 |
| GPU-1 split, L1 | 1369 / 1337 | 485 / 467 | 677 | 212 | 59.0 | 562.084822098171003 |
| **GPU-2 L1** | **1082 / 1100** | **176 / 182** | 385 | 241 | 59.4 | 560.660808090795854 |
| GPU-2 L2 | 1398 / 1392 | 198 / 198 | 447 | 484 | 59.2 | = CPU, bit for bit |

| 135k (`tunnel8:135000`, 134,980 → 162k faces) | step | precond | solve | GMRES its/step |
|---|---|---|---|---|
| CPU 4 threads (chunked) | 10140 ms | 1644 | 5989 | 64.2 |
| CPU 16 threads | 7364 | 1440 | 4711 | 64.2 |
| GPU-1 split, L1 (GPU-1 notes: 5.59 s) | 5071 | 2215 | 2854 | 64.2 |
| **GPU-2 L1** | **3432** | **563** | 1154 | 64.4 |
| GPU-2 L2 (= CPU bit for bit) | 4446 | 547 | 1364 | 64.2 |

Against GPU-1: the preconditioner 2.7× (45k) and 3.9× (135k) faster per step, whole steps 1.24× (45k,
1.35 → 1.09 s) and 1.48× (135k, 5.07 → 3.43 s). Against the CPU: 45k 2.6× (4 threads) / 1.8× (16), 135k
3.0× / 2.1×. The chunked solves also speed up the CPU itself (45k precond 548 → 418 ms, 135k step
11.26 s in GPU-1's notes → 10.14 s).

## Regression set (iterations, fallbacks, energies; `runs_tmp/gpu2/reg/`, 20 steps unless said)

Per-step GMRES iterations summed (mean per step); `=` means identical every step to the CPU default.

| case | CPU (chunked) | CPU faer | GPU-2 L1 | GPU-2 L2 | fallbacks |
|---|---|---|---|---|---|
| Hopf link (scene, 9.2k faces) | 392 (19.6) | = | = | = (energies =) | 0 |
| Borromean rings (scene, 11.9k) | 762 (38.1) | 764 (±1 on 3 steps) | 759 (±1 on 3 steps) | = | 0 |
| Walnut (scene, obstacles, 1.3k) | 149 (7.5) | = | = | = | 0 |
| Lung growing (scene, obstacle, 5.1k) | 346 (17.3) | = | = | = | 0 |
| Klein bottle (scene, non-orientable, self-intersecting) | 977 (48.9) | 865 (43.2; drifts after step 8) | 971 (48.5) | = | 0 |
| H with 24 pins + volume (`apps h_sparse_volume`, 3.8k) | 2302 (115.1) | 2299 | 2301 | = | 0 |
| cobblestones (`apps cobblestones`, 14 components, 4.5k) | 1853 (92.7) | = | = | = | 0 |
| Fig. 13 suite (14 cases, mean its/step) | | 9 of 14 equal, others ±0.2–0.6 | 10 of 14 equal, others ±0.2–0.8 | all = (energies =) | 0 |
| 135k | 64.2 | — | 64.4 | 64.2 (= ) | 0 |

## What is left on the CPU (GPU-2 L1, share of the step)

| part | 45k (1082 ms) | 135k (3432 ms) | note |
|---|---|---|---|
| factor: block structure + far coefficients, vertex BVH, `L′` factor + schedule + upload, constraint set | 233 (22 %) | 820 (24 %) | metric tree 95 / 306, fractional 42 / 142, `L′` 55 / 223 |
| remeshing | 181 (17 %) | 586 (17 %) | sequential |
| Barnes–Hut tree, lists, upload | 101 (9 %) | 350 (10 %) | |
| CCD | 75 (7 %) | 203 (6 %) | |
| GMRES vectors (solve − products − precond) | 55 (5 %) | 187 (5 %) | |
| projection (CPU part) | 53 (5 %) | 150 (4 %) | |
| on the GPU: precond 176 / 563, products 149 / 404, energy 41 / 115, line search 19 / 53 | 385 (36 %) | 1135 (33 %) | |

## Stage 2 and GMRES vectors: not done, why

* **Stage 2** (a preconditioner cheap by construction) was conditional on stage 1 leaving the
  preconditioner a large share. It is now 16 % of a step at both sizes, and ≈ 40 % of that is the
  vertex-pair product, which FSAI/AMG would not touch; the two `L′` solves are ≈ 2.2 ms (45k) / 5–6 ms
  (135k) per application, ≈ 8–10 % of the step. An AMG V-cycle on the GPU would have its own ≈ 30–50
  launches per application and would likely need ≥ 2 cycles to keep the iteration counts; FSAI of
  `W + εM` cannot represent the global `(−Δ)⁻¹` that the preconditioner's order bookkeeping
  (`rs-preconditioner-order`) needs, so iterations would grow with the mesh. Expected gain ≤ 5 % of a step
  with a real risk on the regression set: not worth it now.
* **GMRES vectors on the GPU** would remove the 55 / 187 ms of host vector work and the PCIe copies of
  every product and application (≈ 0.3 ms per iteration): ≤ 6–7 % of a step. Level 2 then needs the CPU
  to compute the GMRES dot products and the dense constraint rows of `A′` (`C_Q x`, barycenter rows
  dense over a component) in GPU-friendly fixed orders (e.g. the pairwise tree of `gpu::pairwise_dot`,
  chunked rows) — a round-off change of every CPU flow — plus a device-resident `A′` (metric product,
  Willmore term, `C_Q`, `C_Qᵀ`) and a second GMRES driver mirroring `GmresRun`'s control flow. Deferred
  behind the vertex-pair product kernel, which costs as much per step and is a local change.

## Proposed SPEC text (not applied; SPEC.md is the orchestrator's)

§6.3, after "(the multiplier part is zeroed between the two solves; `L′` factored once per step by sparse
LDLT)":

> `L′` is factored by faer's sparse `LDLᵀ` in the order "AMD of the Laplacian block, then the
> multipliers", postordered, and solved by `solve::ldlt` in a fixed chunked order (each sum of a
> triangular solve in chunks of 32 terms summed from zero, chunk sums subtracted in order; forward terms
> increasing, backward decreasing), the same on the CPU and the GPU (`TriangularSolve::Chunked`, default;
> `Faer` keeps faer's supernodal solves).

§8.5, replace "On the CPU: GMRES, `L′` and its solves, …" by:

> On the GPU also (GPU-2): the whole fine preconditioner `P⁻¹` — the two `L′` solves by a schedule of
> the elimination tree (leaf subtrees one thread each, chains one block each, levels by height/depth;
> chunk sums from earlier levels in parallel), the fractional product and the augmented factor, one
> upload and one download per application (`PreconditionerOptions::gpu_solves`, default on). On the CPU:
> GMRES, the factorization of `L′` and its schedule, the coarse correction, tree structures and far
> coefficients, constraints, CCD, remeshing; the face-pair `L_{2−s}` and the obstacles' evaluators.
> Measured (RTX 4080, 4 CPU threads): preconditioner per step 2.7× (45k) / 3.9× (135k) faster than GPU-1,
> steps 1.08 s (45k figure-8 tunnel; GPU-1 1.35 s, CPU 2.83 s) and 3.43 s (135k; GPU-1 5.07 s, CPU 10.1 s);
> GMRES iterations unchanged (level 2 bitwise the CPU on the regression set).

Decision row:

> | D35 | `L′` solves | — | faer supernodal solves | chunked fixed-order solves on CPU and GPU, whole preconditioner on the GPU (§6.3, §8.5) |

## How to resume

* Build/test: `CARGO_TARGET_DIR=$PWD/target-gpu2 cargo test --release -p ddg-repulsive-surfaces
  --features cuda` (`--test gpu`, `--lib ldlt`); `DDG_REQUIRE_GPU=1` makes a missing device a failure.
* Measure: `cargo run --release --features cuda -p ddg-repulsive-surfaces --example gpu2_precond --
  fig12:3` (one application; `--pins K`); `--example gpu2_ldlt -- fig12:3` (factor, schedule);
  `cargo bench ... --bench step -- --mesh fig12:3 --gpu L1 [--no-gpu-solves --faer]`.
* Scripts of the passes: `runs_tmp/gpu2/` (`timedB/`, `reg/`, logs with `uptime` at start and end).
* Code: `src/solve/ldlt.rs` (factor, solves, schedule, emulation), `src/solve/precond.rs`
  (`TriangularSolve`, `gpu_solves`, `attach_gpu_solves`), `src/gpu/cuda_precond.rs`,
  `src/gpu/kernels.cu` (`ldl_*`, `prec_aug`, `vec_*`), `src/gpu/cuda.rs` (`apply_device`).
