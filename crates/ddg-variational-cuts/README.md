# ddg-variational-cuts

A clean-room Rust implementation, from the text of the paper only, of

> Nicholas Sharp, Keenan Crane. *Variational Surface Cutting.*
> ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201356>

The paper moves a cut on a triangle mesh by a gradient flow so that the surface, cut open, can be
flattened conformally with little area distortion while the cut stays short. The discrete method
as implemented here, with its sign conventions, is specified in [`SPEC.md`](SPEC.md); every place
where it differs from the paper as printed is listed in **SPEC §9** (D1 to D13). The authors' code
was not consulted (it is GPL; this workspace is MIT).

## What is implemented so far

The **static evaluation** of a cut (SPEC §2, §3.1–§3.7, §6, §7), for two regions described by one
level-set function `φ` (region 0 where `φ ≤ 0`, region 1 where `φ > 0`):

| stage | SPEC | what |
|---|---|---|
| curvature per face | §3.1 | Gaussian, or an override: uniform, per face |
| cut | §3.2 | crossings `(edge, t)`, segments, moved crossings counted |
| cut mesh | §3.3 | split faces, the two sides separated, geometry from barycentric coordinates |
| scale factor | §3.4 | `L u = −g` in every region, `u = 0` on the cut and on `∂M` |
| normal derivative | §3.5 | integrated flux `v`, dual length `ℓ`, **smoothed flux recovery** |
| energy and adjoint | §3.6 | Dirichlet, scaled Dirichlet (`a` per vertex or curvature based), Hencky, pointwise `b, b′` |
| speed | §3.7 | `s = G_0 − G_1` with `w_r`, `β_r`; optional sliver smoothing of the paper |
| energy | §6 | `E = Σ_r [w_r J_r + β_r Area_r] + ∫_γ (α_L + η) ds`, every term |

Not here yet: the flow itself (harmonic extension, implicit length step, symmetry, redistancing;
SPEC §3.8–§3.12, §4) and `N` regions (SPEC §8). The types are shaped for both: region ids are
`u32`, per-region data is per vertex of the cut mesh, and the factorization of an evaluation is
kept for the harmonic extension.

## Use

```rust
use ddg_core::generators::icosphere;
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{CutProblem, Energy, EvalOptions, Evaluator};

let sphere = icosphere(3)?;
let mesh = TriMesh::new(sphere.positions.len(), sphere.faces.clone())?;
let geometry = Geometry::from_positions(&mesh, &sphere.positions)?;

let mut evaluator = Evaluator::new(&mesh, &geometry)?;      // once per mesh: symbolic analysis
let phi: Vec<f64> = sphere.positions.iter().map(|p| 0.3 - p[2]).collect();
let problem = CutProblem::new(Energy::Dirichlet, 0.1);       // energy, α_L; weights are fields
let e = evaluator.evaluate(&phi, &problem, &EvalOptions::default())?;

e.energy.total;            // E, and its terms next to it
e.speed;                   // s at every cut point: dφ/dt on the cut, before the length term
e.log_scale;               // u on the vertices of the cut mesh
e.cut_mesh()?;             // TriMesh + Geometry + maps: what a flattening consumes
```

* `Evaluator` owns the symbolic factorization (the pattern of the system does not depend on the
  cut) and, after `evaluate`, the numeric factorization of that cut
  (`Evaluator::factorization`). On a closed mesh it is the matrix of the harmonic extension.
* `CutProblem` is plain data and may change between evaluations: `energy`, `curvature`,
  `regions` (`w_r`, `β_r`), `length_weight` (`α_L`), `penalty` (`η` per vertex).
* `EvalOptions`: `min_crossing_offset` (`δ`), `flux_recovery`, `sliver_aspect_ratio`,
  `smoothing_passes`, `parallelism`. `EvalOptions::paper()` is the procedure of the paper.
* Bad input returns a `CutError` whose message names the vertex, edge or face; nothing panics.

**Units.** The method runs on the surface scaled to unit total area (`UnitAreaScale`). Parameters
and every number of an `Evaluation` refer to that surface; only `CutMesh::geometry` is in the
caller's unit. Only signs and ratios of `φ` matter for an evaluation.

**The cut mesh** is the cut *surface*: a cut point has one vertex per side, so every vertex lies
in one region and `u`, `v`, `ℓ`, `∂u/∂n`, `p`, `q`, `∂p/∂n`, `G_r` are plain per-vertex arrays.
Vertices of the input keep their index; `CutTriangles::point_vertices` gives the two copies of a
cut point.

## The normal derivative on the cut (SPEC §3.5, D12)

The paper divides the integrated normal derivative of each cut vertex by its dual length,
`∂u/∂n ≈ v/ℓ`, and smooths the result at sliver vertices. Where the cut passes close to a mesh
vertex this fails by orders of magnitude, and sliver smoothing does not repair it. The default here
smooths the fluxes along the cut, implicitly and with the mass `ℓ`, before dividing:

```text
ℓ_p d_p + Σ_{segments σ = pq at p} κ_σ (d_p − d_q) = v_p,     κ_σ = (c · h_f)² / |pq|,   c = 0.5
```

(`h_f²` the mean squared edge length of the face of the segment). It conserves the total flux,
stays within the range of the nodal values, and is the same for the adjoint. Measured against exact
solutions, largest / mean error in percent over all cut vertices and both sides
(`examples/estimator_study.rs`; sphere caps cut by 20 random planes, and a flat annulus on a
perturbed grid with a source that varies along the cut with 12 periods):

| `∂u/∂n` | paper `v/ℓ` | + 1 sliver pass | interpolate over slivers | segment gradient | clusters `θ = 0.25` | smoothed `c = 0.5` |
|---|---|---|---|---|---|---|
| sphere, 642 vertices | 5·10⁴ / 35 | 5·10⁴ / 69 | 24 / 2.1 | 39 / 3.8 | 12 / 1.4 | **3.1 / 0.52** |
| sphere, 2562 | 3·10⁴ / 14 | 3·10⁴ / 23 | 10 / 1.0 | 72 / 2.0 | 9.2 / 0.75 | **1.5 / 0.24** |
| sphere, 10242 | 5·10⁴ / 28 | 5·10⁴ / 54 | 5.4 / 0.56 | 66 / 1.0 | 3.9 / 0.38 | **0.75 / 0.13** |
| sphere, 40962 | 5·10⁴ / 7.7 | 5·10⁴ / 15 | 1.7 / 0.30 | 61 / 0.51 | 2.6 / 0.18 | **0.39 / 0.06** |
| annulus, 540 | 4·10³ / 22 | 4·10³ / 30 | 104 / 30 | 60 / 16 | 87 / 9.0 | **22 / 9.6** |
| annulus, 2040 | 5·10³ / 12 | 5·10³ / 18 | 78 / 26 | 38 / 9.0 | 49 / 4.9 | **9.1 / 3.3** |
| annulus, 7920 | 4·10⁴ / 27 | 4·10⁴ / 39 | 77 / 21 | 73 / 4.6 | 37 / 2.3 | **4.8 / 1.1** |
| annulus, 31200 | 5·10⁴ / 13 | 5·10⁴ / 19 | 70 / 11 | 79 / 2.4 | 25 / 1.2 | **2.3 / 0.44** |

The price is a bias where `∂u/∂n` varies along the cut over a few edges (annulus with 540
vertices: 5 edges per period). `c` is the option `FluxRecovery::Smoothed { length }`;
`FluxRecovery::Nodal` is the paper's quotient.

## Verification

`cargo test -p ddg-variational-cuts` (about half a minute in a debug build; the tables below come from
the ignored tests, see the end).

**Finite differences decide the signs** (`tests/finite_differences.rs`, SPEC §6 item 3). For a
perturbation `φ + εψ` the cut point `c` sweeps the area `dA_c` into region 0, and SPEC §1 predicts
`dE = Σ_c (G_0 − G_1)_c dA_c`. The central difference of the discrete energy agrees with it in
sign and size for every energy, for each region alone, and with `w_r`, `β_r`. Largest difference
relative to `Σ_c |s_c dA_c|` over 9 problems and 2 directions `ψ`, default options | paper options:

| sphere, non-symmetric cut | 162 | 642 | 2562 | 10242 | 40962 vertices |
|---|---|---|---|---|---|
| default | 5.1·10⁻² | 1.8·10⁻² | 4.4·10⁻³ | 1.1·10⁻³ | 2.7·10⁻⁴ |
| paper | 4.8·10⁻² | 2.8·10⁻² | 2.8·10⁻² | 7.9·10⁻³ | 1.2·10⁻² |

| bumpy disk, cut ends on `∂M` | 289 | 1089 | 4225 | 16641 vertices |
|---|---|---|---|---|
| default | 2.0·10⁻¹ | 5.0·10⁻² | 1.7·10⁻² | 3.9·10⁻³ |
| paper | 2.9·10⁻¹ | 7.3·10⁻² | 2.2·10⁻² | 2.9·10⁻² |

A wrong sign of one term would give an error of 1 or 2. This confirms SPEC D2: the region where
`φ ≤ 0` enters the speed with a plus, and a region whose own gradient density is large shrinks.

**Sphere cap** (`tests/analytic.rs`, SPEC §6 item 2), polar angle 1, relative errors:

| level | vertices | max \|u − u*\| | `J_0` | `∂u/∂n` max / mean, region 0 | region 1 | `s` max / mean |
|---|---|---|---|---|---|---|
| 3 | 642 | 2.7·10⁻³ | 5.5·10⁻³ | 1.7·10⁻² / 5.1·10⁻³ | 7.3·10⁻³ / 2.2·10⁻³ | 1.7·10⁻² / 5.3·10⁻³ |
| 4 | 2562 | 6.7·10⁻⁴ | 1.4·10⁻³ | 5.9·10⁻³ / 2.4·10⁻³ | 4.8·10⁻³ / 1.3·10⁻³ | 1.0·10⁻² / 2.7·10⁻³ |
| 5 | 10242 | 1.7·10⁻⁴ | 3.8·10⁻⁴ | 3.8·10⁻³ / 1.3·10⁻³ | 4.1·10⁻³ / 8.2·10⁻⁴ | 9.5·10⁻³ / 1.7·10⁻³ |
| 6 | 40962 | 4.2·10⁻⁵ | 9.6·10⁻⁵ | 1.7·10⁻³ / 6.4·10⁻⁴ | 1.3·10⁻³ / 3.5·10⁻⁴ | 2.6·10⁻³ / 7.3·10⁻⁴ |

With the paper's `v/ℓ` the largest error of `s` on the same four meshes is 13 %, 930 %, 150 %,
6700 % (and the same with one sliver pass).

**Lean oracle** (`tests/oracle.rs`). The 21 committed fixtures of `lean/oracle/fixtures/cuts`, the
output of an independent dense implementation written from `SPEC.md` alone, are replayed: the
scaling, face curvatures, crossings, every triangle of the cut mesh (nodes, lengths, areas), the
segments, `u`, `p`, `v`, `ℓ`, `q` at every node, `J` and the areas per region, the energy, and
(with the nodal recovery the fixtures were generated with) `G` and `s` before and after sliver
smoothing: 59 954 numbers, largest difference `3·10⁻¹²` relative to `max(1, |value|)`.

Further tests: the flat disk with a uniform source (SPEC §1 "Checks": `u`, `J`, `∂u/∂n`, the Hencky
adjoint `∂p/∂n = −Kρ³/8`); length and penalty terms and their derivative on a sphere; the curvature
identity `Σ_f Ω_f = Σ (2π − Θ_i)`; areas and Euler characteristics of the cut mesh; the mesh that
`cut_mesh()` returns reproduces `u` through the cotan Laplacian of `ddg-core` (an independent
assembly); `Σ_B v = Σ g` in every region; conservation and range of the flux recovery; a cut
through vertices; closed components without a cut; invariance under relabeling and under a change
of unit; equivalent descriptions of one energy; bad input, and a property test that any finite
level-set function evaluates.

## Formula ids

Functions that implement a formula carry `[formula: <id>]`; the ids are those of SPEC §10.

| id | function |
|---|---|
| `cut-face-curvature`, `cut-face-curvature-sum` | `face_curvatures`, `face_curvature_sum_residual` |
| `cut-edge-crossing` | `edge_crossing` |
| `cut-barycentric-length`, `cut-subtriangle-area` | `barycentric_dot`, `subtriangle_area` |
| `cut-yamabe-system` | `LaplaceSystem::assemble` (internal) |
| `cut-weighted-laplacian` | `CutTriangles::apply_laplacian` |
| `cut-normal-derivative` | `normal_flux`, `dual_lengths` |
| `cut-flux-recovery` (new, SPEC D12) | `flux_smoothing_stiffness`, `FluxSolver` (internal) |
| `cut-dirichlet-energy`, `cut-gradient-dirichlet` | `dirichlet_energy_term`, `gradient_dirichlet` |
| `cut-adjoint-scaled`, `cut-gradient-scaled` | `adjoint_source_scaled`, `gradient_scaled` |
| `cut-adjoint-pointwise`, `cut-gradient-pointwise` | `adjoint_source_pointwise`, `gradient_pointwise` |
| `cut-feature-size` | `Evaluator::feature_size_weight` |
| `cut-speed`, `cut-sliver-smoothing` | `cut_speed`, `smooth_slivers` |

## Reproducing the tables

```text
cargo test --release -p ddg-variational-cuts --test finite_differences -- --ignored --nocapture
cargo test --release -p ddg-variational-cuts --test analytic -- --ignored --nocapture
cargo run  --release -p ddg-variational-cuts --example estimator_study -- 3 6
```
