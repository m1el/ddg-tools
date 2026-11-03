# GPU-1 notes (CUDA back end for the hierarchical products)

Status: in progress. Sections: binding choice, design, determinism, done / remaining, measurements, how
to resume.

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
  after 6 steps; on the 45k bench the final energy after 7 steps differs by 0.1 % (remeshing
  amplifies round-off, as between CPU thread counts with different GMRES chunking — same as
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
