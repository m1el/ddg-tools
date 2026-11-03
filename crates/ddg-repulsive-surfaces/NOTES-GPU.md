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
kernels are ~300 lines of CUDA C in `src/gpu/kernels.cu`; moving them to cuda-oxide later is a
port of those kernels only (the host side is a thin wrapper). Revisit when CUDA 13 is on the machine
and cuda-oxide is on crates.io / stable.
