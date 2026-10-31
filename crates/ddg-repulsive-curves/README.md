# ddg-repulsive-curves

A Rust implementation of

> Chris Yu, Henrik Schumacher, Keenan Crane. *Repulsive Curves.* ACM Transactions on Graphics
> 40(2), Article 10, 2021. <https://doi.org/10.1145/3439429>

(with the paper's supplemental material): the tangent-point energy of curve networks in R³,
its descent in the fractional Sobolev–Slobodeckij metric of the paper, hard constraints and soft
penalties for design, collision-guarded time steps, and the Barnes–Hut, hierarchical-matrix and
multigrid acceleration of the paper's Appendix B.

**Status: implemented.** The contract is [`SPEC.md`](SPEC.md). It lists:
* every discrete formula with its equation number and `rc-` formula id;
* where the paper or the reference code is wrong or silent, and what we do (decisions D1–D30);
* the API and the defaults;
* the measurements, in §17 (§17s holds the evaluation's results of record);
* the user's decisions, in §18.

Python bindings: `ddg_tools.RepulsiveSession` and `repulsive_flow` (`crates/ddg-py`). Blender:
**Repulsive Flow (Live)** and **Repulsive Flow (Steps)** ([`blender/README.md`](../../blender/README.md)).

## What is in it

| part | paper | module |
|---|---|---|
| Curve networks (closed and open curves, junctions), OBJ `l` lines; planar mode | Sec. 3, 8 | `network`, `io` |
| The tangent-point kernel and the edge-pair energy, exact and parallel; differential | Eq. 3, 17, 18 | `kernel`, `energy` |
| The fractional Sobolev metric `Ā = B + B⁰` (midpoint or trapezoid rule); L², H¹, H² for comparisons | Eq. 19–21, Supp. 2.1 | `metric` |
| Constraints: barycenter, edge lengths, total length, points, tangents, implicit surfaces; projected gradient and back-projection | Sec. 5.3, 8.1, Eq. 22–23 | `constraints`, `solve` |
| Penalties: surface and obstacle potentials, total length, length difference, field alignment, Hencky strain | Sec. 8.1, Supp. 2.2 | `penalties`, `obstacle`, `surfaces` |
| Target-length schedules and subdivision (packing) | Sec. 8.1 | `step::schedule`, `subdivide` |
| Time stepping: a line search guarded by continuous collision detection (the curves never cross), or the paper's normalized Armijo search | Sec. 5.4 | `step`, `ccd` |
| Soft problems: L-BFGS (memory 30) and Polak–Ribière NCG in the Sobolev metric | Supp. 3.2–3.3 | `optimize` |
| Barnes–Hut energy and differential, block cluster tree, hierarchical metric product | App. B.1–B.2 | `spatial`, `energy::barnes_hut`, `metric::hierarchical` |
| Multigrid saddle solver on curve networks (coarsening with constraints, V-cycle, FMG, flexible CG) | App. B.3 | `multigrid` |
| The authors' scenes (reader for their scene format) | Sec. 8 | `scene`, `examples/scenes.rs` |

All runs are deterministic: the same results at any thread count (D25). Bad input returns an
`RcError` that names what is wrong and where; nothing panics.

### Additions to the paper (all documented, all can be switched off)

* **Near-pair quadrature** (`FlowOptions::near_pair`, §3.4, §17o, §17r). Edge pairs that come
  closer than `0.75 ℓ` are evaluated with an accurate pair integral, blended into Eq. 18. The
  endpoint rule cannot see two edges that are close in their interiors. Without the correction,
  20 of the 228 Hs runs drove edges into each other until the collision guard stopped them.
  Second neighbours use the fold-gated rule: never below Eq. 18, so the paper's corner barrier
  stays. Sec. 9 of the paper suggests this ("adaptive quadrature for edge pairs that are close in
  space"), but did not test it. `None` gives the paper's energy exactly. When no pair is that
  close, the flow is bitwise the paper's.
* **Projection.** The constraint projection finishes with Newton steps once the violation is
  below 100× the tolerance (`newton_polish`). It also stops trials that cannot converge early
  (`projection_early_abort`). Measured separately on 24 knots, the polish takes 0.69× the time
  to 1.1 E* of the chord method, and the early abort 0.76× (§17l).
  A search that gets stuck with the polish is repeated with the paper's chord steps
  (`stuck_retry_chord`, §17p).
* **Edge-collapse guard** (`min_edge_length_ratio`, §17j). A step that would shrink an edge
  below `1e-3` of the mean edge length fails with a named error, before the metric becomes
  singular.
* **Dense solve.** A Schur complement replaces the LU of the whole saddle matrix (about 30× fewer
  flops). `Backend::Auto` picks exact + dense up to 2500 vertices, and Barnes–Hut + multigrid
  above.

## Use

```rust
use ddg_repulsive_curves::{FlowOptions, Method, flow, read_curve_obj, soft_options};

let (network, positions) = read_curve_obj("trefoil0002.obj")?;

// Hard constraints (barycenter, edge lengths), the fractional Sobolev metric and the
// collision-guarded line search: the defaults of SPEC.md §10.
let options = FlowOptions { max_iterations: 200, ..FlowOptions::default() };
let result = flow(network.clone(), positions.clone(), options, None)?;
println!("{:?}: {} steps, energy {:.4}", result.report.status,
         result.report.iterations, result.report.final_energy);

// The soft problem of the paper's Sec. 7: L-BFGS with a Hencky penalty on the lengths.
let soft = soft_options(Method::Lbfgs { memory: 30 }, 10.0, &network, &positions)?;
let result = flow(network, positions, FlowOptions { max_iterations: 200, ..soft }, None)?;
```

`Session` gives the same flow one step at a time (`Session::step`), as the live Blender
operator uses it. Every step returns a `StepReport`: energy, gradient norm, line-search trials,
projection iterations, solver residuals, timings and status. Note that only `Session::run`
honours `time_limit`; `Session::step` does not check it.

## The paper's Sec. 7 comparison

The comparison is reproduced on all 228 knots (Knot128 perturbed + Trefoil100, 1494–3351
vertices) with exponents (2, 4.5). It has two phases:
* phase A: every method, one thread, the paper's 24-minute budget;
* phase B: a 24-knot subset, 64 threads, the machine alone.

The driver is `examples/eval_knots.rs` with `scripts/rc-eval.sh`. Setup, departures and every
number are in `SPEC.md` §17n (first run) and §17s (results of record, 2026-10-31).

Failure rates, phase A, against the paper's Fig. 19:

| method | this crate | paper |
|---|---:|---:|
| Hs projected gradient, ended stuck | **0 %** | 6.5 % |
| soft Hs NCG, ended stuck | **0 %** | 3.3 % |
| soft Hs L-BFGS, ended stuck | **0 %** | 14.8 % |
| L² failed | 100 % | ≥ 97.7 % |
| H¹ failed | 68.4 % | 64.5 % |
| H² failed | 92.5 % | 92.1 % |

All 2052 runs finished with 0 crossings in the audit of every accepted step, 0 invalid curves
and 0 errors. One soft Hs L-BFGS run (trefoil0016) used up its budget 2.5 % above its target.

Phase B, median seconds to 1.1× the reference energy at 64 threads (24 knots):

| method | seconds | reached |
|---|---:|---:|
| Hs projected gradient | 24.2 | 24/24 |
| soft Hs L-BFGS | 24.4 | 23/24 |
| soft Hs NCG | 50.1 | 24/24 |
| H¹ projected gradient | 60.5 | 16/24 |
| Hs normalized Armijo | 97.1 | 2/24 |
| L², H² | – | 0/24 |

The videos in `runs_tmp/rc-eval/runB3/videos` follow the authors' video: one column per method,
frames at equal wall-clock intervals.

## Verification

* **Lean** (`lean/proofs/DdgProofs/RepulsiveCurves/`): 73 `rc-` formula ids, 55 proved with
  Mathlib and 18 cited in `lean/proofs/CITED.md` (the paper's smooth theory: Blatt's energy space,
  the order of the differential, and so on). `scripts/check-formulas.sh` cross-checks every
  `[formula: rc-…]` tag in the Rust code.
* **Oracle** (`lean/oracle/DdgOracle/RepulsiveCurves*`): an independent Float implementation,
  with golden fixtures (`tests/oracle.rs`) and a differential fuzz test
  (`fuzz_against_the_oracle`, ignored by default). The fuzz test reports 47 disagreements in the
  chord-iteration count of one step. These predate the near-pair work (same 47 at c657176); they
  are not explained yet.
* **The authors' program**: Gram matrices, constraint Jacobians and one coarsening level match
  exactly (see "Comparison with the authors' program" below).
* **Tests**: 259 tests, in the modules and 25 integration-test files. They cover:
  * central finite differences of every energy, penalty and constraint;
  * invariances and scaling;
  * the metric against its kernel-matrix form;
  * Barnes–Hut and multigrid against exact and dense;
  * CCD on constructed crossings;
  * every stop status;
  * the application scenes, each with a CCD audit of every step.

## Benchmarks

`cargo bench --bench fig20|threads|components|projdiag|trialrule` (`benches/`, results and method
in `SPEC.md` §17k and §17l). Mean seconds per iteration of a running flow (Freedman–He–Wang
resampled to the paper's sizes, (2, 4.5), barycenter + edge lengths), 1 thread, machine alone,
against the values read from the paper's Fig. 20:

| configuration | n | paper | this crate |
|---|---:|---:|---:|
| Exact | 1000 / 4000 | 6 / 145 s | 0.76 / 13.9 s |
| Barnes–Hut + dense | 2000 / 4000 | 10 / 60 s | 1.8 / 5.9 s |
| Multigrid | 2000 / 8000 / 16000 | 6 / 25 / 57 s | 1.7 / 3.7 / 5.2 s |
| Multigrid without projection | 2000 / 16000 | 3 / 38 s | 1.3 / 4.1 s |
| Multigrid, knot 96 ((3, 6)) | 2048 / 8192 / 16384 | 6 / 25 / 57 s (the paper's curve) | 2.5 / 18.5 / 33.0 s (20 steps) |

The dense solve is a Schur complement of `A + ρ𝟙𝟙ᵀ` (about 30× fewer flops than the LU of the
whole saddle matrix); the line search's projection finishes with Newton steps below 100 times the
tolerance and gives up early on trials whose chord iteration cannot converge (both measured on 24
knots, `SPEC.md` §17l; the paper's projection is `newton_polish: None`,
`projection_early_abort: false`). `Backend::Auto` switches from exact + dense to Barnes–Hut +
multigrid above 2500 vertices (measured crossover 2000 at 8 threads, 3000 at 64).

These timings predate the near-pair default. On these curves no pair is near, and the pass costs
0.3–1.7 % of an iteration (`benches/nearpair.rs`, `SPEC.md` §17o).

## Credits

* The method, the discretization, the metric, the constraints and penalties, the acceleration
  scheme, the knot datasets (Knot128, Trefoil100) and the Freedman–He–Wang meshes are by
  **Chris Yu, Henrik Schumacher and Keenan Crane**.
* The authors' reference implementation (`repulsive-curves`, MIT, Copyright (c) 2019 Christopher
  Yu) was read. The implicit surfaces (`surfaces`) are ported from it, and the block cluster
  tree (`metric::hierarchical`) follows its structure. Its licence notice is in
  [`LICENSES/repulsive-curves-MIT.txt`](LICENSES/repulsive-curves-MIT.txt).
* The V-cycle follows the authors' libgmultigrid (MIT, Copyright (c) 2020 Christopher Yu,
  [`LICENSES/libgmultigrid-MIT.txt`](LICENSES/libgmultigrid-MIT.txt)).
* L-BFGS: D. C. Liu, J. Nocedal, *On the limited memory BFGS method for large scale
  optimization*, Math. Programming 45 (1989). The rate test that stops a chord projection early
  is from E. Hairer, G. Wanner, *Solving Ordinary Differential Equations II*, §IV.8.
* The tangent-point energy and its energy space: S. Blatt, *The energy spaces of the tangent
  point energies*, J. Topol. Anal. 5 (2013); S. Blatt, Ph. Reiter, *Regularity theory for
  tangent-point energies: the non-degenerate sub-critical case*, Adv. Calc. Var. 8 (2015).
* Conservative continuous collision detection: M. Li, D. M. Kaufman, C. Jiang, *Codimensional
  Incremental Potential Contact*, ACM TOG 40(4), 2021 (additive CCD); closest-point routines
  after C. Ericson, *Real-Time Collision Detection*, 2004.

## Licence

MIT, like the rest of the workspace.

## Comparison with the authors' program

The authors' program (`repulsive-curves`, MIT) is built outside this repository with a headless
batch mode and a fixture driver (`RepulsiveCurves-build/batch-mode.patch`,
`dump-fixtures.patch`; `SPEC.md` §17c). Our energy is the paper's edge-pair discretization and the
reference's is a vertex-pair one (D1), so **energies are never matched**. What is comparable is
compared.

**Matched exactly** (`tests/reference_fixtures.rs`, fixtures in `tests/fixtures/reference/`; trefoil
with 40 vertices, a 99-vertex subsample of the Freedman curve, the full 592-vertex Freedman 600
(digest), an open curve with pins, a junction network, small loops; exponents (3, 6) and (2, 4.5)):

| quantity | agreement |
|---|---|
| dense Gram matrix `A = B + B⁰`, and `B`, `B⁰` separately (midpoint rule, edge tangents) | relative to the largest entry 4e-16 to 5e-15 (tolerance 1e-12); the trapezoid rule differs by 20 % to 200 % |
| edge-length, total-length and point constraint Jacobians | 5e-16 (tolerance 1e-13); our rows are the negatives of the reference's (we differentiate `ℓ⁰ − ℓ`) |
| one coarsening level: surviving vertices, neighbours of the removed ones, prolongation weights, coarse edges (as sets of fine indices; pins as `Point` and `Tangent(Vertex)`) | identical |

**Not comparable**, with reasons in the test file's header: energies and differentials (D1);
barycenter Jacobian (D6; the reference uses constant dual-length weights, checked against its
formula); tangent and surface constraints and penalties (D7, D8, D19, D20); coarsening of closed
loops with fewer than 6 vertices (D17; the test asserts the documented difference), coarsening beyond
level 1 and coarse positions (D17); block-cluster-tree and BVH statistics (D15, D16).

**Compared by final shape** (`examples/compare_reference.rs`, `runs_tmp/rc-wp15/`): the reference's
batch mode (`--no-bh`: exact energy, dense solver; Armijo with the Sobolev–L² start; stop when the
Sobolev dot is below 1e-4) against `flow` with the exact energy and the dense solver, constraints
edge lengths + barycenter, absolute projection tolerance 1e-4, `threads = 5`
(the reference with `OMP_NUM_THREADS = 5`), the guarded line search of the time (commit
`6af096e` plus working tree; the later multigrid, Barnes–Hut adjoint and line-search-default
changes are not in these runs), stop `‖g̃‖₂ < 1e-4`, 5 h wall-time cap. All 12 jobs ran at the
same time on a machine that other jobs were also using (load about 100 to 200), so **times are
approximate**. Curves: Freedman 600 with (2, 4.5) and the first five `perturbed` Knot128 knots
with (3, 6) (about 2000 vertices; dense LU of 6000 unknowns).
Energy: our Eq. 17 evaluated on the output curve. The reference's output is the last OBJ it wrote
(every 10th step).

| curve | vertices | input energy | reference: steps, time, our energy | ours: steps, time, our energy | RMS distance of the curves after arc-length resampling and Procrustes |
|---|---|---|---|---|---|
| Freedman 600 | 592 | 2281 | 504 (converged), 2493 s, 6.0821 | 647 (converged), 2159 s, 6.0815 | 0.0045 (curve length 21.5, bounding box 2.7) |
| Knot128 perturbed 0001 | 1998 | 1652 | 446, cap, 34.54 | 785, cap, 3.258 | not converged |
| 0002 | 1611 | 1638 | 721, cap, 21.86 | 1485, cap (plateau, 1.1935), 1.1935 | not converged |
| 0003 | 1653 | 1569 | 688, cap, 4.324 | 1439, cap (plateau), 1.9958 | not converged |
| 0004 | 2175 | 2315 | 359, cap, 52.54 | 805, cap, 9.906 | not converged |
| 0005 | 2160 | 1982 | 358, cap, 48.96 | 760, cap (plateau), 2.2565 | not converged |

On Freedman 600 both programs reach the same shape and the same energy of our discretization to
1e-4 relative; the remaining distance is a small fraction of the curve size. The five knots did
not converge within 5 hours in the reference, which makes about 1 step/min at this size under
load; ours takes about 2.5 steps/min and larger steps (on 0002, 0003 and 0005 it reaches an energy
plateau with gradient norm below 1e-2 and no decrease in the 8th digit, but the stop threshold
`1e-4` is not met), so the final curves are not the same shape (the knots expand with the fixed edge
lengths, bounding box 130 for the reference and 350 for ours on 0001) and a Procrustes distance is
not meaningful; only energies, steps and times are reported. The reference was still descending when the cap ended; how much longer it would need was not measured.
