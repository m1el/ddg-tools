# GPU-1 notes (CUDA back end: hierarchical products and Barnes–Hut on the GPU)

Status: done (2026-11-03). Sections: binding choice, design, determinism, tests, measurements, what
stays on the CPU, GMRES vectors, proposed SPEC text, how to resume.

## Done / remaining

* Done: cargo feature `cuda`; `Device` selection; GPU products of `B̃`, `B̃⁰` (+ obstacle diagonal)
  and of the vertex-pair `L̃^V_{2−s}`; near fields and `H̃a` built on the device at level 1;
  Barnes–Hut energy and exact gradient (and value-only trials) on the device; determinism levels 1
  and 2 with tests (level 2 = the CPU flow bit for bit); benches `--gpu`; measurements below.
* Not done (by decision, below): GMRES vectors on the GPU; the `L′` back substitutions (a numerical
  change, out of scope); face-pair `L_{2−s}` (reference mode) and the obstacles' cross evaluators on
  the GPU; Python/Blender exposure of `Device` was done afterwards in `ddg-py` and `blender/` (it added `gpu::devices()` and `Session::gpu_device_name()` here).

## Binding: `cudarc` + NVRTC (not NVIDIA's CUDA Rust), checked 2026-11-03

The user pointed at the official route. What exists today (read from the repository
`github.com/NVIDIA/cuda-rust` = `NVlabs/cuda-oxide`, its `rust-toolchain.toml` files, the cuda-oxide
book's installation page, cutile-rs's README and `cutile-book/reference/compatibility.md`, and
crates.io, on 2026-11-03):

* **NVIDIA "CUDA Rust"** (Apache-2.0, announced 2026-09, "an early-stage research project. Expect bugs,
  incomplete features, and API breakage"; last commit 2026-11-03), three parts:
  * **cuda-oxide** (SIMT kernels in Rust → PTX ahead of time): a custom rustc codegen backend driven by
    `cargo oxide`; needs the pinned Rust **1.99.1** with `rustc-dev`, `rust-src`, `llvm-tools`, **LLVM
    21+ with NVPTX and clang 21+** (bindgen for the host bindings), and the **CUDA Toolkit 13.0+**
    (`cuda.h` of 13.x; driver R580+). Ubuntu 24.04 tested, Linux only, `sm_80+`. Not on crates.io
    (the crate name `cuda-oxide` there is an unrelated 2021 GPL crate).
  * **cutile-rs** (tile model, proc-macro, JIT to Tile IR → cubin by `tileiras`): stable Rust 1.89+,
    but needs the **CUDA Toolkit 13.2+** (13.3 recommended) for `tileiras`; `sm_89` is supported from
    Tile IR 13.2. Crates `cutile`, `cuda-core`, `cuda-bindings` 0.4.0 on crates.io (MSRV 1.89).
  * Host runtime crates `cuda-bindings` / `cuda-core` / `cuda-async`: generated with bindgen against
    the CUDA 13 headers.
* **Rust-CUDA** (`rust-gpu/Rust-CUDA`, `rustc_codegen_nvvm` + `cust`): rebooted, nightly-only codegen
  on libNVVM; last crates.io releases 2022.
* **cudarc** 0.19.10 (MIT/Apache-2.0, 2026-09-24, 11M downloads): safe driver API + NVRTC bindings,
  `dynamic-loading` (libcuda / libnvrtc opened at run time, nothing linked or needed at build time),
  feature `cuda-12040` for this machine's toolkit; kernels are CUDA C compiled at run time by NVRTC.

This machine: driver 595.99 (CUDA 13.x capable), **toolkit 12.4 only** (`nvcc`, `libnvrtc.so.12`),
RTX 4080 `sm_89`, workspace toolchain `stable` = rustc 1.93.1 (1.99.1 is installed by someone else
but not the workspace's). The task allows no system installs.

| criterion | cuda-oxide | cutile-rs | cudarc + NVRTC |
|---|---|---|---|
| (1) Rust stable rule | pinned 1.99.1 + `rustc-dev` + custom codegen backend: would have to live in a separate kernel crate built by `cargo oxide`; the default build stays stable only if that crate is outside the workspace | stable 1.89+, OK | stable, OK (no MSRV declared; builds on 1.93) |
| (2) FMA / float flags per kernel | not documented; Rust semantics forbid contraction, so plain `a*b+c` stays unfused through LLVM, but there is no per-kernel switch to *allow* contraction (level 1's speed) | per-op `fma`, rounding and `ftz` arguments; contraction of the `+ *` operators by `tileiras` not documented | NVRTC `--fmad=true/false` per compiled module (documented PTX semantics: `mul.rn`/`add.rn` are never fused) |
| (3) maintenance, CUDA 12.4 / `sm_89` | active, research grade; **needs CUDA 13.0+ toolkit + LLVM 21 + clang 21** | active; **needs CUDA 13.2+ `tileiras`**; `sm_89` OK | maintained; CUDA 11.4–13.4 feature flags, 12.4 present here, `sm_89` via NVRTC `-arch` |
| (4) cost for CPU-only users | none if the kernel crate is optional, but anyone enabling the feature needs the pinned toolchain, LLVM 21, CUDA 13 | none (optional dep), CUDA 13.2+ to use it | none: optional dependency, nothing CUDA at build time even with the feature on (dynamic loading); at run time only the driver and `libnvrtc.so.12` |
| fit for this workload | SIMT, would fit | tile model for dense tensors; the BCT near field is sparse/irregular (the README itself recommends SIMT for "sparse, pointer chasing") | SIMT, fits |

**Choice: cudarc + NVRTC.** The official tracks both require a CUDA 13 toolkit (not installed; no
system installs allowed), cuda-oxide additionally a pinned compiler with `rustc-dev` and LLVM/clang 21,
and both are labelled research-grade with API breakage expected. cudarc builds on the workspace's
stable toolchain, loads the driver at run time (so `--features cuda` builds on machines without CUDA),
and gives an explicit, documented FMA switch per module (`--fmad=false`), which level 2 needs. The
kernels are ~550 lines of CUDA C in `src/gpu/kernels.cu`; moving them to cuda-oxide later is a
port of those kernels only (the host side is a thin wrapper). Revisit when CUDA 13 is on the machine
and cuda-oxide is on crates.io / stable.

## Design

* **Feature** `cuda` (off by default) adds the optional dependency `cudarc` (driver + NVRTC,
  `dynamic-loading`, `cuda-12040`). Without it nothing CUDA is compiled; the API is the same (the GPU
  types are uninhabited and `GpuContext::new` returns `RsError::Gpu`).
* **Selection** (`src/gpu/mod.rs`): `HierarchicalOptions::device: Device` with `Device::Cpu`
  (default), `Device::Cuda(CudaOptions { ordinal, determinism })` (an error if not usable) and
  `Device::Auto(..)` (GPU if usable, else CPU). `Device::resolve` opens a `GpuContext` (driver
  context, one stream, kernels compiled by NVRTC for the device's `sm_XY`: ~1 s). `Session::new`
  resolves once and keeps the `Arc<GpuContext>`; it is passed explicitly to
  `hierarchical_saddle_with`, `StepSolver::build`, `Objective::build`. `hierarchical_saddle`
  (without a context) resolves `options.device` itself. No global state.
* **What runs on the GPU** (per step, trees built per step as before):
  * `B̃`, `B̃⁰` (+ obstacle diagonal) products of `HierarchicalOperator` (`GpuFaceMetric`), and the
    vertex-pair `L̃^V_{2−s}` products inside the preconditioner (`GpuVertexFractional`). Pipeline per
    product: `face_m` (M x, a M x per face, BVH order) → upward sums (leaves, then one launch per
    depth) → far field per node → `bct_rows` (one thread per (row, channel): far part via the chain
    root → leaf, then the near blocks of the chain) → `face_z` → `face_gather` (vertex gather over
    the CSR incidence) → download. The face-pair `L_{2−s}` (reference mode) stays on the CPU.
  * The near fields: level 1 builds the block structure and far coefficients on the CPU
    (`BlockTree::build_structure`) and evaluates the near entries and `H̃a` on the device
    (`near_fill_face`, `near_fill_vertex`; nothing large crosses PCIe); level 2 builds them on the CPU
    (as before) and uploads them (≈ 290 MB at 45k faces, ≈ 33 ms), because the device `pow` is not
    the CPU's.
  * The Barnes–Hut energy and exact gradient (`BhTree::attach_gpu`; lists uploaded once per step,
    face data and node aggregates per evaluation): `bh_own`, `bh_node`, `bh_partner`,
    `bh_leaf_energy`, `bh_node_sums`, `bh_final`; the downward pass over the nodes (11k at 45k) and
    the pull-back to the vertices stay on the CPU. The value-only `try_energy` (line-search trials,
    the evaluation after the projection) too. Obstacle trees (cross evaluators) stay on the CPU.
* **What stays on the CPU**: GMRES (vectors, Gram–Schmidt, Givens), the sparse `L′` LDLᵀ factor and
  its back substitutions in the preconditioner, the BVH/face tree and list builds, the block
  structure and far coefficients, constraints and projection, CCD, remeshing, the Willmore and
  constraint rows of `A′`.

## Determinism

* Every kernel computes one output per thread in the operation order of the CPU loop it replaces;
  every sum that the CPU does sequentially is sequential in one thread (rows of the block tree,
  far lists, upward sums node by node, gathers over the CSR incidence, the Barnes–Hut lists and
  their transposes); the downward accumulations are gathers along the chain root → leaf. No atomics,
  no block-shape-dependent reductions: results do not depend on the device or the launch
  configuration, so **level 1** (`Determinism::RunToRun`, default) is bitwise reproducible.
* **Level 2** (`Determinism::MatchCpu`) compiles with `--fmad=false` (also `--prec-div=true`,
  `--prec-sqrt=true`, `--ftz=false`), uses the CPU-built near fields, and runs the Barnes–Hut
  kernels only for integer exponents (`ipow` by the same binary powering; for non-integer exponents
  the CPU's `powf` cannot be reproduced and the CPU path is kept). The CPU loops are compiled
  without FMA (`block_tree::simd` uses AVX without FMA). Result: every product, energy and gradient
  is the CPU's bit for bit, so the whole flow is the CPU flow bit for bit (`tests/gpu.rs`
  `flow_trajectories_match_cpu`: 6 steps with remeshing, energies/GMRES iterations/face
  counts/positions identical; the 45k bench's final energy after 7 steps is the CPU's to 17 digits).
* The global sums of the flow stay on the CPU (GMRES dots and norms in fixed 4096-chunks, energy
  totals over the leaves summed on the host in leaf order). The fixed-shape pairwise tree
  reduction asked for the GPU (`gpu::pairwise_sum` / `pairwise_dot` on the CPU, `GpuContext::dot` with
  the `pairwise_level` kernel: pad with `+0.0` to a power of two, add adjacent pairs level by level;
  shape depends only on `n`) is implemented and tested bitwise CPU = GPU
  (`pairwise_dot_matches_cpu`) but not used by the flow, since GMRES stays on the CPU (below); the
  CPU default is unchanged.
* Level 1 against the CPU: products within `3e-16` relative (FMA), device-built near fields
  within `1e-11` (pow), flows: same GMRES iterations on the test flow, energies within `1e-12`
  after 6 steps; on the 45k bench the final energy after 7 steps differs by 0.16 % (remeshing
  amplifies round-off; the same spread as between equivalent settings in
  WP17's tolerance study).

## Tests (`tests/gpu.rs`, `--features cuda`; skip with a note when no device unless `DDG_REQUIRE_GPU`)

* `face_metric_products_match_cpu`: 5 meshes × χ ∈ {0, 0.25, 0.5, 1} × leaf sizes × 5 weight/diagonal
  cases × batches 1–6: level 2 bitwise, level 1 `< 1e-12` and reproducible.
* `vertex_fractional_products_match_cpu`: same for `L̃^V_{2−s}`.
* `gpu_built_near_fields_match_cpu`: level 1 device-built `H̃a` (`< 1e-13`) and products
  (`< 1e-11`), reproducible; the structure-only CPU trees refuse CPU products.
* `barnes_hut_matches_cpu`: energy and exact differential, θ ∈ {0, 0.25, 0.5}, exponents (6, 12) and
  (4, 9), before and after a refit: level 2 bitwise, level 1 `< 1e-10`; (5.5, 11) stays on the CPU at
  level 2.
* `flow_trajectories_match_cpu`: whole flows (above).
* `pairwise_dot_matches_cpu`, `device_selection` (bad ordinal: `Cuda` is an error, `Auto` falls back).

## Measurements (timed pass 2026-11-03, `runs_tmp/gpu/timed/`)

RTX 4080 (16 GB, `sm_89`, f64 at 1/64 of f32), driver 595.99, CUDA 12.4 NVRTC; Threadripper 3970X.
One run at a time; outside load (another session's `trialrule-*` jobs, which had pushed the
1-minute load to 55 an hour earlier and were winding down): 1-minute load average 7.6 → 2.0 over the
pass (recorded at the start and end of every run in the logs). The CPU baselines of the same pass
agree with the quiet-machine baseline taken earlier (45k step 3.04 s at load 3.4; products 23.4 ms).

**Products** (`benches/product.rs --gpu`, perturbed torus, median of 7; GPU times include the
host→device copy of the vectors and the copy back):

| | CPU 1 thr | CPU 4 thr | CPU 16 thr | GPU L1 | GPU L2 | L1 vs 4 / 16 thr | L2 vs 4 / 16 thr |
|---|---|---|---|---|---|---|---|
| 45k faces, `B̃+B̃⁰`, 1 vector | 61.8 ms | 17.4 | 11.6 | 1.81 | 2.48 | 9.6× / 6.4× | 7.0× / 4.7× |
| 45k, 2 vectors | 80.2 | 25.0 | 16.7 | 2.41 | 3.97 | 10.4× / 6.9× | 6.3× / 4.2× |
| 45k, 4 vectors | 188.0 | 51.5 | 27.2 | 4.18 | 7.06 | 12.3× / 6.5× | 7.3× / 3.9× |
| 45k, vertex `L̃^V`, 2 vectors | – | 2.40 | – | 1.17 | 1.41 | 2.1× | 1.7× |
| 135k faces, 2 vectors | – | 79.2 | 44.9 | 6.69 | 11.73 | 11.8× / 6.7× | 6.8× / 3.8× |
| 135k, 4 vectors | – | 246.4 | 117.2 | 13.06 | 21.00 | 18.9× / 9.0× | 11.7× / 5.6× |

The near field (17.9M entries at 45k, 55.5M at 135k) is read once per product: the kernels are
bound by the f64 rate of this consumer card, not by memory (level 2 without FMA is 1.6× slower than
level 1). Uploading the CPU-built near field (level 2): 32–37 ms at 45k, 105–110 ms at 135k;
level 1 builds it on the device instead (metric tree 108 ms vs 280 ms per step at 45k, including the
CPU block structure). Against 16 threads level 1 is 6.4–9×, level 2 3.8–5.6× (below the 5×
reassessment mark for 1–2 vectors; level 2 is the verification mode).

**Whole steps, 45k figure-8 tunnel** (`fig12:3`, 44,848 → 53,544 faces, 4 threads, barycenter +
area, warm-up 2, mean of 5; interleaved CPU, L1, L2, CPU, L1, L2):

| | step | BVH+lists | energy | factor (metric tree / fractional / `L′`) | solve (products / precond) | line search | CCD | projection | remesh | GMRES its/step | F after 7 steps |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CPU | 3218 / 3159 ms | 101 | 361 | 450 (280 / 80 / 53) | 1789 (1057 / 644) | 100 | 74 | 134 | 210 | 59.6 | 561.208135632529547 |
| GPU L1 | 1522 / 1583 | 119 | 42 | 238 (108 / 50 / 47) | 781 (143 / 559) | 18 | 75 | 52 | 197 | 59.0 | 562.084822098171003 |
| GPU L2 | 1865 / 1886 | 119 | 50 | 483 (270 / 79 / 50) | 854 (212 / 563) | 20 | 74 | 53 | 212 | 59.6 | 561.208135632529547 (= CPU) |

Speed-up 2.1× (L1), 1.7× (L2). The paper's figure-8 tunnel: ≈ 2.7 s/step (floor ≈ 2.5 s); ours now
1.44–1.58 s at level 1. Level 2's trajectory is the CPU's bit for bit (energies of every step
identical, `f8_*_steps.csv`); level 1's iterations are within round-off effects (59.0 vs 59.6) and
its energy after 7 steps differs by 0.16 % (remeshing amplifies round-off; WP17's tolerance study
shows the same spread between equivalent settings).

**135k faces** (`tunnel8:135000`, 134,980 → 162,277 faces, same settings): CPU 11.26 s/step, GPU L1
5.59 s (2.0×); GMRES 64.4 / 64.2 iterations per step. L1 breakdown: precond 2409 ms (43 %),
solve rest 228, products 443, factor 782, remesh 743, BVH 425, CCD 203, projection 162, energy 131.

**Fig. 13 suite** (4 threads, warm-up 2, mean of 5; CPU and GPU L1 back to back):

| # | model | faces (→ remeshed) | paper | CPU | GPU L1 | speed-up | GMRES its CPU / L1 |
|---|---|---|---|---|---|---|---|
| 1 | triple torus | 6,200 → 6,758 | ≈0.25 | 0.211 | 0.108 | 1.96× | 17.4 / 17.4 |
| 2 | moving bars (HandcuffsPole) | 10,688 | ≈0.45 | 0.339 | 0.200 | 1.70× | 17.2 / 17.2 |
| 3 | broken figure-8 | 11,740 → 11,4xx | ≈0.37 | 0.490 | 0.268 | 1.83× | 38.6 / 38.4 |
| 4 | punctured torus | 9,886 → 10,046 | ≈0.65 | 0.390 | 0.208 | 1.88× | 36.6 / 36.6 |
| 5 | handcuffs | 10,960 → 12,338 | ≈0.45 | 0.394 | 0.209 | 1.88× | 22.2 / 22.2 |
| 6 | 2×-wound handcuffs | 16,020 → 18,270 | ≈0.75 | 0.703 | 0.369 | 1.90× | 32.0 / 31.8 |
| 7 | point cloud | 20,480 | ≈1.6 | 0.677 | 0.440 | 1.54× | 14.6 / 14.6 |
| 8 | nested cages (mesh obstacle) | 20,480 | ≈1.75 | 0.750 | 0.502 | 1.49× | 13.2 / 13.2 |
| 9 | walnut (genus-0 slab) | 23,328 → 25,634 | ≈1.4 | 0.734 | 0.438 | 1.67× | 12.6 / 12.6 |
| 10 | trefoil tunnel | 23,344 → 26,584 | ≈1.4 | 1.077 | 0.559 | 1.93× | 39.0 / 39.2 |
| 11 | 6×-wound handcuffs | 21,640 → 24.4k | ≈1.35 | 1.340 | 0.626 | 2.14× | 61.6 / 62.2 |
| 12 | double torus | 35,454 → 35,530 | ≈1.8 | 1.644 | 0.860 | 1.91× | 31.2 / 31.2 |
| 13 | keys (genus-4 slab) | 37,312 → 39,772 | ≈2.05 | 1.405 | 0.704 | 2.00× | 18.6 / 18.6 |
| 14 | figure-8 tunnel | 44,848 → 53.5k | ≈2.7 | 2.956 | 1.441 | 2.05× | 59.6 / 59.0 |

(seconds per step). All 14 at or below the paper with the GPU (CPU: 12 of 14, #3 and #14 above).
#7–#8 gain least: their obstacles' cross evaluators (point cloud, mesh obstacle) stay on the CPU.

## What stays on the CPU (45k, level 1, share of the 1.52 s step)

| part | ms | share | note |
|---|---|---|---|
| preconditioner `L′` back substitutions | ≈ 490 | 32 % | sequential triangular solves per vector (2 vectors in parallel); the vertex-pair products inside the preconditioner (≈ 70 ms of its 559) are on the GPU |
| remeshing | 197 | 13 % | sequential |
| factor (block structure + far coefficients, vertex BVH, `L′` LDLᵀ, constraint set) | 238 | 16 % | GPU fill included |
| Barnes–Hut tree, lists, GPU list upload | 119 | 8 % | |
| GMRES vectors (solve − products − preconditioner) | 79 | 5 % | see below |
| CCD | 75 | 5 % | |
| projection (CPU part; its energy evaluation is on the GPU) | 52 | 3 % | |
| on the GPU: products (incl. PCIe) 143, energy 42, line search 18 | 203 | 13 % | |

At 135k the `L′` solves are 43 % of the step. **Next step:** a GPU-friendly preconditioner (the
`L′` solves are the largest remaining share; a sparse triangular solve on the GPU or a different
inner solver for `L′` — e.g. a multigrid or a few CG iterations on `W + εM` — would change the
numerics, so it needs its own study of GMRES iterations), then remeshing in parallel and the
structure/far-field build of the block trees on the device.

## GMRES vector operations: not moved (plan step 3), why

With right preconditioning every Arnoldi step needs `P⁻¹ v` on the CPU (the `L′` back
substitutions), so each iteration crosses PCIe twice whatever happens to the Gram–Schmidt; the
transfers are ≈ 0.65 MB per vector at 45k faces (≈ 30 µs each, included in the product times
below). What GMRES itself costs on the CPU (solve − products − preconditioner) is ≈ 65 ms of a
step at 45k faces, 4 threads (≈ 4 % of the GPU step). On the GPU, modified Gram–Schmidt with the
fixed-shape tree reduction is one reduction (≈ 17 dependent launches) per basis vector and
iteration: ≈ 60 iterations × 30 basis vectors × 17 launches ≈ 30k launches per step, more than
the CPU's 65 ms. Classical Gram–Schmidt (one batched pass) would change the arithmetic (a
numerical change, not round-off). So the vectors stay on the CPU; the pairwise-tree reduction is in
place (`GpuContext::dot`) for when the preconditioner moves.

## Proposed SPEC text (not applied; SPEC.md is the orchestrator's)

New **§8.5 GPU back end (optional; GPU-1)**:

> Cargo feature `cuda` (off by default; `cudarc` with the CUDA driver and NVRTC opened at run time,
> kernels in `src/gpu/kernels.cu` compiled for the device on first use). `HierarchicalOptions::device:
> Device` = `Cpu` (default) | `Cuda(CudaOptions { ordinal, determinism })` | `Auto(..)`; `Cuda` is
> `RsError::Gpu` when not usable (no feature, driver, device or NVRTC), `Auto` falls back to the CPU.
> A `Session` opens one `GpuContext` and passes it explicitly; no global state. On the GPU: the
> products of `B̃`, `B̃⁰` with `O_diag` and of the vertex-pair `L̃^V_{2−s}` (§8.3), the near fields
> and `H̃a` (level 1), the Barnes–Hut energy and exact gradient (§3.4, §8.2). On the CPU: GMRES,
> `L′` and its solves, tree structures and far coefficients, constraints, CCD, remeshing; the face-pair
> `L_{2−s}` and the obstacles' evaluators. Determinism: every kernel is one thread per output in the
> CPU's operation order, gathers only; level 1 (`RunToRun`, default) is bitwise reproducible on any
> device; level 2 (`MatchCpu`: `--fmad=false`, CPU-built near fields, Barnes–Hut on the device only
> for integer `α`, `β`) is bitwise the CPU flow. Measured (RTX 4080, 4 CPU threads): products
> 10–12× the CPU at 4 threads (6.4–9× at 16), whole steps 1.5–2.1× (Fig. 13 suite; 45k figure-8
> tunnel 3.2 → 1.5 s, 135k 11.3 → 5.6 s), GMRES iterations unchanged; the `L′` solves are then
> 32–43 % of a step.

Decision row:

> | D34 | GPU acceleration | — | — (CPU only) | optional `cuda` feature, explicit `Device`; products, near fields, Barnes–Hut on the device; level 1 bitwise reproducible, level 2 bitwise the CPU (§8.5) |

§14.5 addition: "With `Device::Cuda` (RTX 4080) all 14 Fig. 13 models are below the paper's times
(0.11–1.44 s/step), the 45k figure-8 tunnel 1.44–1.58 s against the paper's ≈ 2.7 s."

## How to resume

* Build/test: `CARGO_TARGET_DIR=$PWD/target-gpu cargo test --release -p ddg-repulsive-surfaces
  --features cuda --test gpu` (`DDG_REQUIRE_GPU=1` makes a missing device a failure).
* Benches: `cargo bench -p ddg-repulsive-surfaces --features cuda --bench product -- --threads 1,4,16
  --gpu` (CPU rows, then `gpu-L1` / `gpu-L2` rows); `--bench step -- --mesh fig12:3 --gpu L1`
  (`--suite fig13`; the scene cases need `scripts/rs-fetch-scenes.sh`). Scripts and logs of the runs
  below: `runs_tmp/gpu/` (not committed).
* Code: `src/gpu/{mod.rs, cuda.rs, kernels.cu}`; hooks in `metric/hierarchical.rs` (`upload`,
  `build_on_gpu`), `spatial/block_tree.rs` (`build_structure`, `parts`), `solve/iterative.rs`
  (`HierarchicalOptions::device`, `hierarchical_saddle_with`, `HierarchicalOperator::with_gpu`),
  `solve/precond.rs` (`new_on`, `with_gpu`), `energy/barnes_hut.rs` (`attach_gpu`, `try_energy`,
  `try_face_partials`), `step/{session.rs, objective.rs, direction.rs}`.
