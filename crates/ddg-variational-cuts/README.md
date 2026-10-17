# ddg-variational-cuts

A clean-room Rust implementation, from the text of the paper only, of

> Nicholas Sharp, Keenan Crane. *Variational Surface Cutting.*
> ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201356>

The paper moves a cut on a triangle mesh by a gradient flow so that the surface, cut open, can be
flattened conformally with little area distortion while the cut stays short. The discrete method
as implemented here, with its sign conventions, is specified in [`SPEC.md`](SPEC.md); every place
where it differs from the paper as printed is listed in **SPEC §9** (D1 to D31). The authors' code
was not consulted (it is GPL; this workspace is MIT).

## What is implemented

The **flow** for two regions (SPEC §3.8–§3.12, §4, §5; sections "The flow" to "Timings" below), on
top of the **static evaluation** of a cut (SPEC §2, §3.1–§3.7, §6, §7), for two regions described by one
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

The same for **`N` regions** (SPEC §8; section "Several regions" below): one level-set function
per region, cut networks with triple points, [`RegionFlow`]. The two-region API is unchanged; the
structs `Cut`, `CutTriangles` and `EvalReport` gained fields for triple points, and a copy of a
triple point is reported as `CutVertex::Point(n_points + m)`.

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

## The flow

```rust
use ddg_variational_cuts::{CutProblem, Energy, Flow, FlowOptions, StepSize, StopRule, init};

let phi0 = init::geodesic_circle(&mesh, &geometry, 0, 0.3)?;        // or any numbers per vertex
let mut flow = Flow::new(&mesh, &geometry, &phi0)?;                 // owns evaluator and state
let mut problem = CutProblem::new(Energy::Hencky, 0.0);
let omega = flow.roughening_rate(&problem, &Default::default())?;   // scale of α_L and of τ
problem.length_weight = ddg_variational_cuts::length_weight_for_ripples(omega, 0.25);
let options = FlowOptions::new(StepSize::default());                // FlowOptions::paper(..): as printed
let report = flow.run(&problem, &options, &StopRule::default())?;   // or flow.step(..) with new controls
let cut = flow.evaluate(&problem, &options.eval)?.cut.polylines(&mesh, &positions)?;
```

One step (SPEC §4). Each line is an option; the right column is `FlowOptions::paper`.

| stage | default | paper, selectable |
|---|---|---|
| normal derivative | smoothed flux recovery (D12) | `v/ℓ` and sliver smoothing |
| length and penalty term | on the cut: turning angle of the cut polyline, implicit along the cut (D14) | on the mesh: `(W_V + τ L_w) φ′ = W_V φ̃` |
| extension of the rate | from the cut, conservative, no solve (D15) | harmonic on the cut mesh |
| rate of `φ` | speed times the mean `|∇φ|` of the piece of the cut (D16) | speed |
| penalty `η` | full gradient `∂η/∂n + ηκ` (D8) | `ηκ` |
| step | `Adaptive`: `min(0.05/ω, 0.25/m̂)`, `ω` and the speed in edges `m̂` as 0.9-quantiles along the cut, and a limit per cut point (D22) | `Fixed(τ)`; `Capped` (Sec. 7.4) |
| redistancing | cut-adjacent values kept, fast marching elsewhere | same |

The three changes of the default were forced by the first law below; each cause and each
alternative that was measured and rejected is in SPEC D14–D16.

### Laws the flow is tested against

`examples/flow_study.rs` prints the tables; `tests/flow.rs` asserts them on small meshes. Full
tables: SPEC §6.

**Curve shortening** (flat square, circle of radius 0.3, `α_L = 1`, to `t = 0.03`, where the exact
radius is 0.17321). Relative error of the radius for 6 / 24 / 96 / 384 steps:

| cells | default stages | stages of the paper |
|---|---|---|
| 32² | `+3.7e-2 +9.6e-3 +2.0e-3 −4e-5` | `−16 % −17 % −17 % −20 %` |
| 64² | `+3.9e-2 +1.1e-2 +2.8e-3 +8.4e-4` | `−12 % −11 % −11 % −15 %` |
| 128² | `+3.9e-2 +1.1e-2 +2.9e-3 +8.8e-4` | `−9.4 % −7.5 % −6.4 % −9.1 %` |
| 256² | `+4.0e-2 +1.1e-2 +2.8e-3` | `−8.3 % −5.8 % −4.0 %` |

The default stages are first order in the step with an error that does not depend on the mesh (it
is the error of the explicit Euler method for the radius). The paper's step does not converge in
the step and converges slowly in the mesh width.

**Latitude circle on the sphere**: `ψ′ = cos ψ/(R² sin ψ) · (4/(R sin ψ) − α_L)`. Error of `ψ` with
50 / 200 steps on icosphere(5): `α_L = 8`, from 0.8 to 1.115: `+2.1e-3 / +5.3e-4` (paper stages:
`−2.2e-3 / −7.5e-3`); `α_L = 20`, from 0.5 to 0.732: `+1.7e-3 / +3.7e-4` (`−4.5e-2 / −3.0e-2`). The
balanced latitude `sin ψ* = 4/(α_L R)` (0.788 for `α_L = 20`) is reached to `±0.005` from both
sides; the paper's step stops at 0.739 and 0.741. Without a length weight the flow is ill-posed
and the law can be followed on a coarse mesh over a short time only (icosphere(3), `t = 10⁻³`:
`10⁻⁴`–`6·10⁻⁴`).

**Ripples of the equator**: `ε cos(mθ)` grows at the rate `(m − 1)[4/R³ − α_L (m + 1)/R²]`
(SPEC §1, D18): the distortion term roughens the cut, the length term smooths it, and the equator
is stable iff `α_L > 4/(3R) = 4.73`. Measured on icosphere(5) (exact in brackets):

| `α_L` | `m = 2` | `m = 3` | `m = 4` | `m = 6` |
|---|---|---|---|---|
| 2 | `+101 (+103)` | `+150 (+155)` | `+155 (+158)` | `+10 (+11)` |
| 4 | `+19 (+27)` | `−49 (−46)` | `−207 (−219)` | `−833 (−868)` |
| 6 | `−52 (−48)` | `−228 (−247)` | `−577 (−596)` | `−1690 (−1748)` |
| 10 | `−174 (−199)` | `−612 (−649)` | `−1317 (−1350)` | `−3405 (−3507)` |

With `α_L = 0`: `+529 (+535)` and `+871 (+891)` for `m = 4, 6`. This law is what `α_L` and the time
step are scaled by: every step reports the roughening rate `ω`; ripples longer than `2π α_L/ω`
grow, and `length_weight_for_ripples(ω, ℓ)` returns the weight for a given `ℓ`.

**Stopping** (SPEC §4, D27, D31). A run ends after `max_steps`, when the cut is at rest, when it
has vanished, or by the energy over a window of steps (default 25). At rest means that a step
moved the level-set function at the cut by at most `motion_tolerance` *and* that the rate of the
cut predicted no more for that step (`StepReport::predicted_motion`: the motion without the
local limit and with unit slope). The motion of one step alone can be small because something
held the step back: a cut through vertices of an intrinsic Delaunay export moved by `1.5·10⁻¹⁰`
in its first step, with a rate of 7, and the run was reported as converged (the cause was the
mean slope, which counted faces whose crossing had been moved off a vertex, SPEC D30). The
energy criterion: if the energy
changed by at most `10⁻⁴` of its size the run is *stalled* and stays in its last state; if it
*rose* by more than that, the run ends with `StopReason::EnergyRose` and the flow is set back to
the state of lowest energy of the run. Every `RunReport` has the lowest energy, the step at which
it was reached, the step of the state returned and the energy after the last step. A run that
only fluctuates at its noise floor can end either way; how far the last energy lies above the
lowest tells it from a run that went wrong (0.4 % in the third run of the table of "Several
regions", 64 % in the same run before SPEC D26, which was the reason for the rule). A tolerance
of zero turns a criterion off (also the motion criterion after a step that moved nothing), and
the energy criterion needs more than a window of steps. A number of a step that is not finite (a
term of the energy, a speed, the step size, a level-set value) is an error that names it and the
step, `CutError::NonFiniteStep`, never a stop reason; the flow stays in the state before.

Also tested: energy decrease in each of the four energies for small steps; symmetry (exact with
classes; to rounding without, unless the diagonals of cut quads tie, SPEC §3.11); independence of
numbering (`10⁻¹⁶`) and scaling (`10⁻¹⁵`); a straight cut across a linear penalty moves at the
speed of the slope with the gradient rule and not at all with the paper's rule (D8); the three
ways a run ends; resuming from a stored `φ`; changing every control between steps; bad input.

### Lean oracle

`tests/oracle.rs` replays the five fixture sets of `lean/oracle` (an independent implementation
written from `SPEC.md`), reading them by the oracle's README alone. Every stage of one full step
is compared (`τ` and what it was computed from, the length force and both rates at every cut
point, `ḡ`, the extension, `φ̃`, `φ′`, symmetry, the projection of Losasso et al., adjacency,
the redistanced functions, counters, the new energy and cut), with the statistics of the step
rule that the fixture names:

| set | fixtures | numbers compared | largest difference relative to `max(1, \|value\|)` |
|---|---|---|---|
| `cuts`: the paper's stages with the flux recovery of SPEC §3.5 | 26 | evaluation 70907; `d, e, G, s` 6218; step 24444 | `1.0·10⁻¹²`; `3.4·10⁻¹²`; `1.6·10⁻¹³` |
| `cuts-default`: the default stages (D14–D16), each stage choice, the three step rules | 23 | 67349; 5691; 32391 | `1.0·10⁻¹²`; `8.1·10⁻¹⁴`; `7.2·10⁻¹³` |
| `regions`: `N` regions, both stage sets, 15 fixtures with triple points | 19 | 125667 | `3.9·10⁻¹¹` |
| `cuts-robust`: the robust step rule (D22, D26), both stage sets, 8 robust and 1 with the extremes | 9 | 28211; 2573; 17065 | `8.2·10⁻¹⁴`; `1.7·10⁻¹⁴`; `4.5·10⁻¹⁴` |
| `regions-robust`: the same for `N` regions | 7 | 48514 | `8.2·10⁻¹⁴` |

The first three sets have the step rule with the extremes (`StepStatistics::Extremes`, the rule
before D22); with the robust statistics, which are the default, `τ` or the limited cut points
differ in every one of their fixtures whose step is not fixed (9 of `cuts`, 21 of `cuts-default`,
11 of `regions`; the test prints their names). The last two sets add the quantile `ω`, both
speeds `m̂` of the step rule, `h_c` and `ℓ_c`, and at every cut point or boundary node the factor
of the local limit and the rate before and after it.

Six of the fixtures of the last two sets have the robust statistics with the length step on the
mesh, where the local limit acts at the vertices on the whole step (SPEC D26). The oracle
implemented that rule after the crate; they replay in full, with `h_i`, `ω_i` and the factor at
every vertex and the limited function, and no difference was found.

The oracle does not have the rule of SPEC D30 yet (a face with a moved crossing does not count
in the mean slope); none of its fixtures has such a face with the mean-slope scaling, so the five
sets agree as they are.

One difference was found earlier, in `rr-sphere3-spike-paper`: the two copies of the cut point of the
quantile were limited by the factor `1 − 2·10⁻¹⁶` here and not at all in the oracle. That cut
point lies on the limit exactly, and the rounding of `τ` decided. SPEC §4 now says that a cut
point at or below the quantile is not limited (D24), the crate does that, and the five sets
agree. The oracle can also smooth the feature size with the cotan Laplacian, which the crate
does not have (D17): `m-featurelap-cotan` has no negative edge weight and agrees as it is, for
`r-torus5-paper` the replay takes the oracle's weight `a` per vertex. The helpers for initial
cuts are not in the oracle.

One case that the oracle's fixtures did not have is written by an ignored test, answered by the
oracle's executable and replayed like a fixture (the oracle has since taken the two files into
`regions-robust`): a region that is a single vertex whose
neighbours alternate between two other regions, so that its boundary runs through triple points
only (robust statistics, 4 boundary nodes limited). Default stages: 5241 numbers, largest
difference `1.6·10⁻¹⁴`; paper's stages: 4737, `1.5·10⁻¹⁴`. `tests/regions.rs` checks the same
state without the oracle (the weights of the quantiles add up to the length of the network, no
boundary node moves by more than `max_motion · h`).

### Figures of the paper, with numbers

`flow_study fig1`, `fig5`, `balance`; lengths on the unit-area surface.

**Fig. 1** (sphere, equator with a perturbation of `0.01` in the modes 2–8, Dirichlet energy,
default options, icosphere(5), in brackets icosphere(4)):

| `α_L` | result | cut length | distortion | steps |
|---|---|---|---|---|
| 6 | back to the equator (largest deviation 0.017 rad) | 1.7724 (1.7723); equator 1.7725 | 4.836 → 4.854 | 105, converged |
| 3 | seam with two lobes, amplitude 0.71 rad | 2.2552 (2.2552) | 4.836 → 3.003 (2.999) | 196, converged |
| 1.5 | seam with two lobes, amplitude 1.49 rad | 2.7593 (2.7582) | 4.836 → 1.883 (1.881) | 418, converged |

Three cuts of three lengths, the longer the less distortion, as in the figure; both regions keep
half of the area (`0.5000 ± 0.0001`). The linear theory predicts which: at 6 every mode decays, at
3 the modes 2 and 3 grow, at 1.5 the modes up to 8, and the fastest is `m ≈ 5`; the end state has
two lobes in both cases, so the ripples coarsen. The energy rose in 2 to 10 of the steps, by at
most `5·10⁻⁴` of the total drop. (With the extremes of the step rule, before D22, the same end
states took 155, 357 and 771 steps. With the motion criterion before D31 the third run ended
after 344 steps, at the same energy to five digits.)

**Fig. 5** (flat unit square, 64² cells, with a bump `0.15 exp(−r²/0.02)`; the cut starts as a
straight line beside the bump):

| energy | start | `α_L` | result | cut length (file units) | distortion |
|---|---|---|---|---|---|
| Hencky | `x = 0.2` | 0.0080, 0.0040, 0.0020 | bends into the bump, passes 0.034, 0.027, 0.023 from its top | 1.000 → 1.081, 1.090, 1.121 | 0.00446 → 0.00099, 0.00093, 0.00085 |
| Dirichlet | `x = 0.36` | 0.74 | passes 0.004 from the top | 1.014 → 1.090 | 0.555 → 0.144 |
| Dirichlet | `x = 0.36` | 0.18 | passes 0.008 from the top, in two pieces | 1.014 → 1.190 | 0.555 → 0.112 |
| Dirichlet | `x = 0.2` or `0.36` | 3 to 12 | retreats to the flat part (`x ≈ 0.07`) | → 1.000 | unchanged |

"A longer cut through the region of high curvature", as in the figure, for the Hencky energy from
far away and for the Dirichlet energy from nearby. Honest remarks: the Dirichlet rates near the
bump are `10⁵` times those at `x = 0.2`, so a length weight that controls the cut at the bump
(`α_L = ω ℓ/2π` with the largest local rate, 93, of a line through its top; the rate of that
state, the 0.9-quantile, is 25) pulls a distant cut away from the bump before the distortion
term is felt (the line across the foot of the bump is not a geodesic). The paper says the same of
this energy and rescales it (Sec. 3.5.1). Two of the three Hencky runs and the Dirichlet run with
0.18 did not come to rest within 3000 steps (they were run without the energy criterion; the
third Hencky run took 256 steps). The Dirichlet run with 0.74 does not end by the motion
criterion either: its energy stands at 0.9367 from step 381 on, where the motion of one step was
below the tolerance, but the cut points at the top of the bump are held back by the local limit
and their rates are not small (before D31 that step ended the run as converged). The Dirichlet run with 0.18 lets ripples of three mesh edges
grow; its energy rose in half of the steps, by up to 1.5 % of the total drop, and the cut ended in
two pieces (one piece, 0.004 from the top, length 1.228 with the step rule before D22).

**Fig. 14** (distortion on one side, area reward): on the sphere a cap with `w = (1, 0)`,
`β = (−8.638, 0)`, `α_L = 2` has its stable size at `ψ* = 1.2`. Caps started at 1.0 and 1.4 end at
1.203 and 1.205 (icosphere(4)), 1.195 and 1.200 (icosphere(5)): it grows, and stops.

**Spot** (`data/models/spot/spot/spot_triangulated.obj`, 2930 vertices; start: the plane across
the longest axis through the centroid, cut length 1.177; default step and stop rule):

| energy | `α_L` | steps, stop, time | distortion | cut length |
|---|---|---|---|---|
| Hencky | 2 | 144, energy rose, 0.5 s | 3.688 → 1.030 | 1.658 |
| Hencky | 0.5 | 200, energy stalled, 0.7 s | 3.688 → 0.339 | 2.431 |
| Hencky | 0.1 | 234, energy stalled, 0.9 s | 3.688 → 0.061 | 3.477 |
| Hencky | auto, `ℓ = 0.25`: 2.17 | 151, energy stalled, 0.5 s | 3.688 → 1.116 | 1.618 |
| scaled, `r_ε = 0.1` | 0.003 | 63, energy stalled, 0.2 s | 0.02045 → 0.01774 | 1.541 |
| scaled, `r_ε = 0.1` | 0.03 | 132, energy stalled, 0.4 s | 0.02045 → 0.02154 | 1.012 |
| scaled, `r_ε = 1` | auto: 0.171 | 35, energy stalled, 0.1 s | 0.4070 → 0.4086 | 1.151 |
| Dirichlet | auto: 35.3 | 32, energy stalled, 0.1 s | 106.7 → 109.0 | 1.089 |
| Dirichlet | 1 | 285, energy rose, 1.1 s | 106.7 → 5.31 | 8.851 |

The energy fell over every run. What the numbers say: the smaller the length weight, the longer
the cut and the lower the distortion, as it must be. No run ended by the motion criterion: on this
mesh the rate of the cut has a noise floor (SPEC §4), the cut keeps moving by a few percent of an
edge per step, and the energy criterion ends the run. With the automatic weight the Hencky run
goes where `α_L = 2` goes; for the scaled and the Dirichlet energy the cut still mostly shortens
and the run ends after 32 to 35 steps with the energy lower by half a percent. Single steps raise
the energy: in 23 of 144 steps by at most 0.6 % of the energy (Hencky, 2), in 27 of 200 and 17 of
234 by at most 0.08 % and 0.03 % (0.5, 0.1), never for the scaled energy with 0.003. ("Energy
rose": the run ended because the energy rose over the last window, and is in its state of
lowest energy, SPEC D27.) The
Dirichlet energy without rescaling is not usable on this model with a weight that lets the cut
move (the last line: speeds up to `9·10⁶`, the cut grows to 7.5 times its length, and the run
ends in a cycle of period two in which the energy alternates by 0.1 % and a fifth of the 1100 cut
points are limited), which is the paper's reason for the rescaled and Hencky energies. The cuts
were not inspected visually.

The same runs with the extremes of the step rule (before D22) took 350, 708 and 1000 steps (the
limit) for the three Hencky weights and ended at the distortion 1.034, 0.350 and 0.112; the
automatic weights were 5.87 (Hencky), 0.666 (scaled) and 135 (Dirichlet), and the Dirichlet run
with `α_L = 1` was at the distortion 27.9 and the length 6.75 after 1000 steps, still growing.

### Timings

One step, one thread, release build, seconds (median of five steps after the first; the paper's
implicit step factorizes once: first step 0.063 s and 0.376 s more):

| mesh | vertices | evaluate | length | extension | redistance | step |
|---|---|---|---|---|---|---|
| icosphere(6), default | 40 962 | 0.0565 | 0.0006 | 0.0024 | 0.0098 | 0.069 |
| icosphere(6), paper stages | | 0.0519 | 0.0066 | 0.0065 | 0.0094 | 0.075 |
| icosphere(7), default | 163 842 | 0.308 | 0.0015 | 0.0092 | 0.043 | 0.362 |
| icosphere(7), paper stages | | 0.315 | 0.033 | 0.034 | 0.055 | 0.438 |
| spot, default | 2 930 | 0.0021 | 0.0002 | 0.0001 | 0.0006 | 0.0030 |
| spot, paper stages | | 0.0021 | 0.0002 | 0.0002 | 0.0007 | 0.0031 |

The evaluation (cut mesh, one numeric factorization, the solves) is 70–85 % of a step. The default
stages need no linear solve after it except two tridiagonal-like systems on the cut. On a closed
mesh the harmonic extension reuses the factorization of the evaluation; on a mesh with boundary it
factorizes a second matrix. One symbolic analysis per pattern, at construction.

### Example program

```text
cargo run --release -p ddg-variational-cuts --example cut_flow -- MESH.obj --energy hencky --alpha 0.5 --out DIR
```

reads an OBJ (polygons are split into fans), takes every control from the command line (energy,
`α_L` or the ripple rule, `w_r`, `β_r`, curvature override, visibility penalty and its rule,
initial cut, symmetry classes, step and stop rule, every stage choice, `--paper`), and writes
`steps.csv` (energy terms, cut length, `τ`, speeds, `ω`, cut motion, counters, `|∇φ|`, stage
times per step), `phi.txt`, `cut.obj` (polylines) and `cut_mesh.obj` (one group per region). The
options are listed at the top of `examples/cut_flow.rs`.

## Several regions (SPEC §8)

```rust
use ddg_variational_cuts::{CutProblem, Energy, FlowOptions, RegionFlow, StepSize, StopRule, init};

let labels = init::normal_clusters(&mesh, &positions, 8, 100)?;              // or labels of your own
let phis = init::regions_from_vertex_labels(&mesh, &geometry, &labels, 8)?;  // one function per region
let mut flow = RegionFlow::new(&mesh, &geometry, &phis)?;
let mut problem = CutProblem::for_regions(Energy::Hencky, 0.0, 8);            // w_k, β_k per region
let options = FlowOptions::for_regions(StepSize::default());                 // see below
problem.length_weight =
    ddg_variational_cuts::length_weight_for_ripples(flow.roughening_rate(&problem, &options.eval)?, 0.5);
let report = flow.run(&problem, &options, &StopRule::default())?;
let network = flow.evaluate(&problem, &options.eval)?;
let curves = network.cut.curves(&mesh, &positions)?;   // polylines between junctions, with their two regions
let triples = &network.cut.triple_points;              // face and barycentric coordinates
let atlas = network.cut_mesh()?;                       // one component per patch, a region per face
```

`examples/cut_flow.rs` does this from the command line (`--init clusters:8`, or a label file with
more than two labels) and writes `network.obj`, `triples.txt`, `regions.txt`, `cut_mesh.obj`.

**`N = 2` reproduces the two-region flow**: bit for bit with the paper's stages, to `3·10⁻¹⁵`
with the default ones (sphere and disk with boundary, three energies, weights, penalty, symmetry
classes, fixed and adaptive steps).

**Which stages.** Both stage sets run on networks. The paper's (every `φ_k` diffused on the mesh,
harmonic extension per region with the substitution of App. A.3, the projection of Losasso et
al. 2006 and fast marching per function) are what `FlowOptions::for_regions` selects, with the
flux recovery of SPEC §3.5. The default stages carry over as "every region moves its own
boundary through its triple points" (SPEC §8.5); two junction forces tried before that failed.
What was measured:

| test | default stages | paper's stages |
|---|---|---|
| three bands, no triple point: excess length over `t = 0.004` (exact factor 0.579) | 0.598 | 0.580 |
| honeycomb on a flat torus given by edge lengths (42 × 24), at rest, 20 steps: junctions drift by | 3.4° | 3.6° |
| perturbed honeycomb, 160 steps: excess length `+2.74 %` → | `+0.71 %` | `+0.78 %` |
| junction angles off 120°: `27.7°` → | `8.7°` | `5.2°` |
| steps in which the length rose | 22 of 161, ≤ 0.18 % | 50 of 161, ≤ 0.19 % |
| the same at 84 × 48, 320 steps: excess, angles | `+0.39 %`, 4.2° | `+0.37 %`, 2.5°–3.8° |
| tetrahedral partition of the sphere, 30 steps: triple points move by (mesh width 0.15) | 0.013 | 0.075 |
| cubical partition | 0.010 | 0.073 |
| energy from a start with five regions, four energies, small steps | falls in every step | falls in every step |
| renumbering, permuting regions, scaling | `10⁻¹⁶` | `10⁻¹⁶` |

Stated plainly: neither stage set keeps the length of a relaxing honeycomb monotone, and neither
brings the junctions to 120° better than a few degrees at these mesh widths (the triple point is
the centroid of three cut points, and the discrete network is longer than the honeycomb by a
multiple of the mesh width per junction). The default stages are the more accurate ones where
there is distortion. `FlowOptions::for_regions` selects the paper's stages; with them, the defect
of the paper's length step on single curves (D14) is back.

**Which to use** (`network_study compare`; SPEC §8.5 "After D26"). Both stage sets from 8 clusters
of the vertex normals on four shapes, with two energies (Hencky; rescaled Dirichlet with
`r_ε = 1`) and three length weights each: `auto`, the ripple rule with `ℓ = 0.5`; `fixed`,
`α_L = 7.4·10⁻⁴`; `strong`, the ripple rule with the largest local rate of the start in place of
`ω`. Default step and stop rule, at most 4000 steps. "Final" is the state the run returns,
"lowest" the lowest energy of the run and the number of steps after which it was reached, "rose"
the number of steps that raised the energy and the largest such rise relative to the energy.
Times: 48 runs at once, one thread each, on a 64-core machine that was otherwise idle (load
average 1 before the start).

| energy | model | weight: `α_L` | stages | final energy | lowest (after steps) | steps that rose, largest | patches | cut length | steps | time, s | stop |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Hencky | spot | auto: 0.00307 | paper's | 2.6675e-02 | 2.6675e-02 (159) | 31, 0.79 % | 59 → 11 | 14.67 → 7.42 | 168 | 2.9 | rose |
| Hencky | spot | auto: 0.00307 | default | 2.7951e-02 | 2.7944e-02 (173) | 43, 1.9 % | 59 → 18 | 14.67 → 8.14 | 181 | 2.7 | stalled |
| Hencky | spot | fixed: 0.00074 | paper's | 8.7118e-03 | 8.7118e-03 (157) | 51, 1.5 % | 59 → 32 | 14.67 → 9.60 | 159 | 2.5 | rose |
| Hencky | spot | fixed: 0.00074 | default | 9.2388e-03 | 9.2320e-03 (102) | 27, 0.2 % | 59 → 39 | 14.67 → 11.10 | 104 | 1.5 | stalled |
| Hencky | spot | strong: 0.0332 | paper's | 1.7165e-01 | 1.7165e-01 (144) | 30, 57 % | 59 → 8 | 14.67 → 4.58 | 149 | 2.5 | rose |
| Hencky | spot | strong: 0.0332 | default | 2.0516e-01 | 2.0500e-01 (175) | 27, 0.24 % | 59 → 7 | 14.67 → 5.76 | 188 | 2.9 | stalled |
| Hencky | bob | auto: 0.00308 | paper's | 2.0069e-02 | 2.0069e-02 (376) | 61, 0.67 % | 23 → 12 | 11.79 → 5.15 | 389 | 13.1 | rose |
| Hencky | bob | auto: 0.00308 | default | 2.0604e-02 | 2.0604e-02 (323) | 41, 0.13 % | 23 → 9 | 11.79 → 5.63 | 337 | 8.0 | rose |
| Hencky | bob | fixed: 0.00074 | paper's | 7.0839e-03 | 7.0761e-03 (232) | 88, 0.24 % | 23 → 16 | 11.79 → 7.52 | 239 | 8.5 | stalled |
| Hencky | bob | fixed: 0.00074 | default | 7.3579e-03 | 7.3579e-03 (190) | 32, 0.098 % | 23 → 18 | 11.79 → 8.53 | 195 | 5.5 | rose |
| Hencky | bob | strong: 0.00946 | paper's | 4.8427e-02 | 4.8427e-02 (281) | 64, 0.73 % | 23 → 9 | 11.79 → 4.29 | 287 | 9.8 | rose |
| Hencky | bob | strong: 0.00946 | default | 5.0758e-02 | 5.0758e-02 (315) | 19, 0.76 % | 23 → 5 | 11.79 → 4.66 | 315 | 7.3 | rose |
| Hencky | blub | auto: 0.00183 | paper's | 1.3682e-02 | 1.3682e-02 (151) | 32, 1.7 % | 44 → 12 | 15.42 → 6.63 | 153 | 6.9 | rose |
| Hencky | blub | auto: 0.00183 | default | 1.3707e-02 | 1.3707e-02 (313) | 68, 0.23 % | 44 → 12 | 15.42 → 6.80 | 315 | 9.9 | rose |
| Hencky | blub | fixed: 0.00074 | paper's | 6.2927e-03 | 6.2927e-03 (183) | 77, 0.92 % | 44 → 12 | 15.42 → 7.06 | 188 | 8.0 | rose |
| Hencky | blub | fixed: 0.00074 | default | 6.9826e-03 | 6.9809e-03 (319) | 93, 0.37 % | 44 → 13 | 15.42 → 8.43 | 322 | 10.0 | stalled |
| Hencky | blub | strong: 0.0239 | paper's | 1.3628e-01 | 1.3628e-01 (170) | 16, 2.8 % | 44 → 11 | 15.42 → 5.40 | 170 | 6.8 | rose |
| Hencky | blub | strong: 0.0239 | default | 1.4278e-01 | 1.4278e-01 (263) | 41, 0.19 % | 44 → 9 | 15.42 → 5.77 | 275 | 8.2 | rose |
| Hencky | icosphere(4) | auto: 0.00786 | paper's | 4.4660e-02 | 4.4645e-02 (27) | 26, 0.027 % | 8 → 8 | 5.36 → 5.00 | 48 | 1.3 | stalled |
| Hencky | icosphere(4) | auto: 0.00786 | default | 4.4689e-02 | 4.4689e-02 (12) | 14, 0.02 % | 8 → 8 | 5.36 → 5.00 | 36 | 0.7 | rose |
| Hencky | icosphere(4) | fixed: 0.00074 | paper's | 6.0384e-03 | 6.0384e-03 (229) | 19, 0.58 % | 8 → 9 | 5.36 → 6.55 | 231 | 4.3 | rose |
| Hencky | icosphere(4) | fixed: 0.00074 | default | 6.0647e-03 | 6.0612e-03 (182) | 40, 0.034 % | 8 → 8 | 5.36 → 6.84 | 199 | 3.3 | stalled |
| Hencky | icosphere(4) | strong: 0.0115 | paper's | 6.2687e-02 | 6.2687e-02 (89) | 18, 0.037 % | 8 → 8 | 5.36 → 4.98 | 89 | 2.2 | stalled |
| Hencky | icosphere(4) | strong: 0.0115 | default | 6.2719e-02 | 6.2713e-02 (74) | 34, 0.037 % | 8 → 8 | 5.36 → 4.98 | 91 | 1.5 | stalled |
| scaled | spot | auto: 0.0486 | paper's | 3.1627e-01 | 3.1627e-01 (176) | 41, 0.3 % | 59 → 11 | 14.67 → 5.46 | 178 | 2.5 | rose |
| scaled | spot | auto: 0.0486 | default | 3.3247e-01 | 3.3247e-01 (251) | 40, 0.11 % | 59 → 8 | 14.67 → 5.92 | 257 | 3.0 | rose |
| scaled | spot | fixed: 0.00074 | paper's | 1.8439e-02 | 1.8439e-02 (62) | 24, 1.6 % | 59 → 92 | 14.67 → 16.91 | 62 | 1.8 | rose |
| scaled | spot | fixed: 0.00074 | default | 1.7958e-02 | 1.7958e-02 (50) | 16, 0.7 % | 59 → 129 | 14.67 → 18.78 | 54 | 1.1 | rose |
| scaled | spot | strong: 1.49 | paper's | 3.2057e+00 | 3.2057e+00 (353) | 25, 1.7 % | 59 → 4 | 14.67 → 1.91 | 373 | 4.5 | rose |
| scaled | spot | strong: 1.49 | default | 2.3699e+00 | 2.3699e+00 (443) | 59, 2.2 % | 59 → 3 | 14.67 → 1.26 | 443 | 3.6 | stalled |
| scaled | bob | auto: 0.0153 | paper's | 8.4213e-02 | 8.4213e-02 (263) | 56, 0.55 % | 23 → 10 | 11.79 → 4.19 | 268 | 9.6 | rose |
| scaled | bob | auto: 0.0153 | default | 8.7910e-02 | 8.7910e-02 (292) | 22, 0.69 % | 23 → 8 | 11.79 → 4.54 | 292 | 7.0 | rose |
| scaled | bob | fixed: 0.00074 | paper's | 1.2271e-02 | 1.2271e-02 (83) | 22, 0.43 % | 23 → 27 | 11.79 → 10.83 | 102 | 4.0 | rose |
| scaled | bob | fixed: 0.00074 | default | 1.2232e-02 | 1.2232e-02 (146) | 34, 0.25 % | 23 → 25 | 11.79 → 11.85 | 148 | 4.8 | rose |
| scaled | bob | strong: 0.859 | paper's | 9.2510e-01 | 9.2510e-01 (267) | 15, 0.91 % | 23 → 3 | 11.79 → 0.94 | 283 | 8.7 | rose |
| scaled | bob | strong: 0.859 | default | 2.2731e+00 | 2.2720e+00 (433) | 35, 0.49 % | 23 → 4 | 11.79 → 2.54 | 447 | 9.4 | stalled |
| scaled | blub | auto: 0.11 | paper's | 6.7054e-01 | 6.7054e-01 (244) | 53, 0.84 % | 44 → 9 | 15.42 → 5.66 | 264 | 10.4 | rose |
| scaled | blub | auto: 0.11 | default | 7.6729e-01 | 7.6729e-01 (252) | 47, 0.23 % | 44 → 10 | 15.42 → 6.67 | 258 | 7.6 | rose |
| scaled | blub | fixed: 0.00074 | paper's | 1.7879e-02 | 1.7879e-02 (47) | 13, 6.7 % | 44 → 83 | 15.42 → 15.64 | 47 | 3.3 | rose |
| scaled | blub | fixed: 0.00074 | default | 1.7348e-02 | 1.7223e-02 (76) | 31, 1.1 % | 44 → 121 | 15.42 → 18.01 | 89 | 3.8 | stalled |
| scaled | blub | strong: 1.39 | paper's | 4.2662e+00 | 4.2662e+00 (250) | 11, 2.9 % | 44 → 6 | 15.42 → 2.25 | 253 | 9.8 | rose |
| scaled | blub | strong: 1.39 | default | 3.4268e+00 | 3.4268e+00 (603) | 128, 5.7 % | 44 → 6 | 15.42 → 1.82 | 617 | 15.1 | rose |
| scaled | icosphere(4) | auto: 0.00786 | paper's | 4.7704e-02 | 4.7699e-02 (28) | 9, 0.034 % | 8 → 8 | 5.36 → 4.98 | 30 | 0.8 | stalled |
| scaled | icosphere(4) | auto: 0.00786 | default | 4.7689e-02 | 4.7683e-02 (48) | 26, 0.028 % | 8 → 8 | 5.36 → 4.98 | 64 | 1.1 | stalled |
| scaled | icosphere(4) | fixed: 0.00074 | paper's | 8.5237e-03 | 8.5135e-03 (177) | 14, 0.041 % | 8 → 8 | 5.36 → 7.90 | 184 | 3.9 | stalled |
| scaled | icosphere(4) | fixed: 0.00074 | default | 8.4662e-03 | 8.4635e-03 (206) | 4, 0.037 % | 8 → 8 | 5.36 → 8.13 | 208 | 3.3 | stalled |
| scaled | icosphere(4) | strong: 0.0104 | paper's | 6.0541e-02 | 6.0484e-02 (24) | 18, 0.054 % | 8 → 8 | 5.36 → 4.95 | 45 | 1.2 | stalled |
| scaled | icosphere(4) | strong: 0.0104 | default | 6.0419e-02 | 6.0419e-02 (38) | 23, 0.024 % | 8 → 8 | 5.36 → 4.96 | 56 | 1.0 | rose |

All 48 runs end below the energy they started with. The paper's stages end at the lower energy
in 16 of the 24 pairs: in all 12 with the Hencky energy (by 0.05 % to 16 %) and in 4 with the
scaled energy (by 4 % to 59 %); the default stages in 8, all with the scaled energy (by 0.03 %
to 3 % in six of them, by 20 % and 26 % on blub and spot with the strong weight). What the
paper's stages cost:

* Single steps raise the energy more: by at most 6.7 % in 23 of their runs and by 57 % in one
  (spot, Hencky, strong: the step before the stop, in which a patch was lost), against at most
  5.7 % with the default stages. 18 of their 24 runs end because the energy rose over the window
  and are set back to their lowest state (SPEC D27), against 13: they rely on that rule.
* For one curve between two regions their length step has the defect of D14: two regions on
  spot (Hencky, `α_L = 2` and `0.5`) end at 4.431 and 1.733 against 4.345 and 1.554.
* On a mesh with spikes of curvature they do not descend (`camelshotnoise.obj`, two regions,
  Hencky: the energy stays above the 5028 of the start, where the default stages reach 4525).

They take fewer steps (4359 against 5491 in total) in about the same time (129 s against 122 s). So
for several regions on meshes like these the paper's stages are the better choice by this table,
and `FlowOptions::for_regions` keeps selecting them; the default stages are the choice for two
regions, for noisy meshes, and where a run must not depend on being set back.

The run that led to this table (spot, Hencky, `α_L = 7.4·10⁻⁴`, the paper's stages) used to end
*above* its start: the energy fell from `1.1955·10⁻²` to `1.104·10⁻²` in 8 steps, rose in every
step after that to `1.8155·10⁻²`, and the run was reported as stalled. The cause was the local
limit of the step rule, which cut back the explicit rate and not the length step on the mesh
(SPEC D26; `examples/rise_study.rs` prints the change of the energy by every stage). The limit
now acts on the whole step, and the run is the third line of the table. The atlases were not
inspected visually and not flattened.

**Different starts (Fig. 10)**, `network_study starts`, 6 regions, one length weight `α_L = 1`
for every start, paper's stages:

| start | sphere (icosphere(4), Dirichlet): patches, length, energy | spot (rescaled Dirichlet) |
|---|---|---|
| clusters of normals | 6, 4.179, 5.300 | 3, 1.13, 1.43 |
| cells around six seed vertices | 6, 4.179, 5.304 | 5, 2.38, 2.71 |
| the same from another first seed | 6, 4.181, 5.302 | 5, 2.56, 2.92 |
| a random label at every vertex | 4, 3.52, 4.93 | 6, 3.23, 3.49 |

The three structured starts end with the same numbers on the sphere (length and energy within
0.1 %; 6 patches, 8 triple points and the length of the cubical partition), which is what the
plots of Fig. 10 show. On spot the two starts from cells end alike (5 patches, length within 8 %),
the clusters of normals (43 patches at the start) coarsen to 3 patches, and the cuts are not in
the same place: 12 % to 47 % of the cut of one of these three results lies within two edges of the
cut of another. The random start (1440 and 1665 patches) ends with fewer patches and a lower
energy than the structured ones on the sphere, and with more patches and a higher energy on
spot: different starts do not all reach the same minimum here.

**Timings**, one step, one thread, seconds (median of five; `network_study timings`):

| mesh | `N` | evaluate | length | extension | redistance | step (default) | step (paper's) |
|---|---|---|---|---|---|---|---|
| icosphere(6), 40 962 vertices | 2 | 0.058 | 0.0013 | 0.0025 | 0.019 | 0.082 | 0.098 |
| | 4 | 0.059 | 0.0021 | 0.0025 | 0.037 | 0.103 | 0.129 |
| | 8 | 0.060 | 0.0030 | 0.0027 | 0.073 | 0.143 | 0.193 |
| spot | 2 | 0.0025 | 0.0007 | 0.0002 | 0.0013 | 0.0048 | 0.0044 |
| | 4 | 0.0029 | 0.0014 | 0.0002 | 0.0026 | 0.0075 | 0.0066 |
| | 8 | 0.0037 | 0.0025 | 0.0002 | 0.0054 | 0.0125 | 0.0106 |

The evaluation does not depend on `N` (one factorization for all regions). What grows with `N`:
the redistancing, one fast marching over the whole mesh per function (the largest part of a step
from `N = 8` on), and the paper's implicit step, one solve per function (0.013, 0.026, 0.051 s on
icosphere(6)). The length systems on the cut grow with the number of cut points.

## Spikes of curvature (SPEC §1, §4, D22)

The roughening rate `ω` of a state, which sets the automatic length weight
(`length_weight_for_ripples`) and the adaptive step, is the 0.9-quantile of the local rates along
the cut, weighted by length; the step rule takes the same quantile of the speed, measured in local
crossed edges, and slows down the cut points beyond it individually
(`StepStatistics::Robust`, the default; `StepStatistics::Extremes` is the rule before D22). With
the length step on the mesh the speed is taken with and without the length force, whichever
gives the smaller step.

| start: a plane, Hencky energy | local rate: 50 % / 90 % / 99 % / largest | `ω` | automatic `α_L` (`ℓ = 0.25`) |
|---|---|---|---|
| `camelshotnoise.obj`, 11 513 vertices | 11 / 44 / 82 / `4.5·10⁶` | 44 | 1.75 (extremes: `1.8·10⁵`) |
| spot | 24 / 55 / 147 / 147 | 55 | 2.17 (extremes: 5.87) |
| bob | 1.8 / 4.5 / 4.6 / 4.6 | 4.5 | |
| icosphere(5) | 3.09 / 3.09 / 3.09 / 3.10 | 3.09 | |

On the camel with shot noise one cut point next to a noise spike has the rate `4.5·10⁶`. With the
extremes the run freezes (611 steps cover the flow time `1.2·10⁻⁷`). With the robust statistics
600 steps cover `5·10⁻³` and the energy falls from 5028 to 3475 while the cut grows from 325 to
775 cut points; up to 75 of them are limited in a step. The run is not smooth: the energy rises in
236 of those steps (by up to 6.7 % of the energy in one step), and with the default stop rule it
ends after 27 steps, set back to the state after 16 (5028 → 4525), because the total energy is
dominated by the noise that no cut removes; at least a tenth of the cut of the camel has a
negative local rate. On
spot the same `α_L` now gives the same or a lower end energy in fewer steps (Hencky, `α_L = 2`:
144 steps instead of 350; `0.5`: 200 instead of 708). Tests: `tests/flow.rs` (a quantile test,
and a sphere with a spike of curvature density on the cut: `ω` 148 against a largest local rate
of 1724, the step `2.8·10⁻⁴` against `2.9·10⁻⁵` with the extremes).

With the length step on the mesh (the stages of the paper) the limit does not act on the rate
at the cut points but at the vertices, on the change of the whole step, after the length step
(SPEC D26): the rate on the cut is only the explicit part of that step, and cutting it back alone
made the cut move uphill (section "Several regions").

What the larger step costs: single steps raise the energy more often and by more than before
(spot, Hencky, `α_L = 2`: 23 of 144 steps, at most `5.8·10⁻³` of the energy, against 13 of 350 and
`1.2·10⁻⁴`), and where the weight lets ripples of a few mesh edges grow the run does not settle
(the Dirichlet runs with small weights under Fig. 5 and Spot; the latter ends in a cycle of period
two instead of crawling). The tests of the laws use fixed steps and are unchanged; what this file
reports from runs with the adaptive step (Fig. 1, 5 and 14, spot, the atlases, the starts) was
measured again with the new rule. On a sphere all local rates agree (they differ by 0.8 % along
the equator for the Dirichlet energy, 0.4 % for Hencky), so the ripple rates that `ω` is verified
against are untouched.

## Mesh quality

`examples/bad_mesh_study.rs` measures every stage and the laws on meshes with the same geometry
and worse triangles (`redist`, `eval`, `laws`, `laws` with `DELAUNAY=1`, `attr`). The numbers
below are from those five runs, repeated with the present code (flat unit square, 32² cells,
circle of radius 0.3; the circle law is the one of "Laws" above, 96 steps):

| mesh | redistancing: mean relative error of the distance | circle law, default stages | circle law, paper's stages |
|---|---|---|---|
| regular grid | 1.2 % | `+2.3e-3` | −16 % |
| vertices jittered by 0.3 cells | 1.4 % | `+2.0e-3` | −17 % |
| sheared, largest angle 117° | 3.3 % | `+5.9e-3` | −41 % |
| sheared, 135° | 11 % | `+1.1e-2` | the cut vanishes |
| sheared, 153° | 34 % | `+1.8e-2` | the cut vanishes |
| sheared, 165° | 91 % | `+2.5e-2` | the cut vanishes |
| cells 1 : 16 | 0.9 % | `+2.7e-3` | −13 % |
| sheared 135° and 153°, edges flipped to Delaunay | | `+2.5e-3`, `+2.5e-3` | −16 % |

Random irregularity is harmless. Where every triangle is obtuse the fast marching degrades first
(iterating it to its fixed point does not repair it: 16 % and 30 % largest relative error at 135°
and 153° against 38 % and 110 %), the evaluation second (largest error of `∂u/∂n` 9.4 % at 135°
against 1.5 % on the regular grid). The default stages still follow the circle law to 1 %–2.5 %
there. The mean-gradient scaling is what takes the wrong slope of the redistanced function out
(mean `|∇φ|` in the crossed faces 1.03 at 153° and 1.81 at 165°): with the unit speed in its
place the errors at 135°, 153° and 165° are −4.2 %, −5.5 % and +18 %. The stages of the paper do
not survive, and the length step on the mesh is the cause: the default stages with that one
option changed lose the cut on the same three meshes. Flipping the same vertices to a Delaunay
triangulation removes the damage. The crate does not flip edges: a caller with obtuse triangles
should.

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

## Comparison with the authors' program (black box)

The authors' program was run as a black box, twice: on 2026-10-13 (before D22) and again on
2026-10-17 with the code of D32 (SPEC §9). Its sources were never opened; it was only built,
launched and driven. Outputs, scripts and pictures are in `runs_tmp/blackbox/` (not in the
repository); the pictures of 2026-10-17 are `runs_tmp/blackbox/pics_2026-10-17/` (one
`compare_<shape>.png` and `compare_<shape>_dirichlet.png` per shape, and `contact_sheet.png`).

**Method.**

* The program has a GUI only and no batch mode. It was driven under Xvfb by synthetic mouse
  clicks and typed values, on the same OBJ files as our runs (spot, bob, blub from Keenan Crane's
  model repository; an icosphere of 2 562 vertices). Every run: initialize by its "Normal
  clustering" (40 patches on spot), set the weights, "Take many steps" with its default fixed
  step size 0.01, then save "obj with injective texcoords".
* What was used: the energies it prints after every step (length, Dirichlet, Hencky, ...), the
  saved mesh (the surface with the cut inserted, and a UV map), screenshots, and wall-clock time.
* Its cut is read from the saved mesh: a seam is an edge whose two sides have different UVs, a
  chart is a UV-connected component. To evaluate it in our energy the charts are carried to the
  input mesh as face labels (nearest face centroid) and our flow is started from those labels with
  0 steps. That turns the cut into a staircase along mesh edges, which makes it longer (spot:
  11.08 instead of 10.55) and adds distortion (`1.42·10⁻³` instead of `1.28·10⁻³`); where small
  charts are lost on the way it doubles the distortion (blub: 16 charts become 14 patches,
  `2.06·10⁻³` instead of `0.82·10⁻³`; the icosphere: 14 become 10). This evaluation favours us.
  A second yardstick avoids the transfer: the area distortion of a
  flattening with `u = 0` on the cut, `∫u²` with `u = ½ log(area ratio)` per triangle, on the
  unit-area surface, plus `α_L` times the cut length; for the reference that is its own UV map
  and its own seams, for us our unwrap (`ddg_uv`, no cones) of our cut.
* **Length-weight correspondence**, unchanged since 2026-10-13: the reference's length term is the
  sum of the perimeters of the regions, twice the cut length, on the surface scaled to unit area
  (it prints `Length scale = √area`); its Hencky value is about 900 times our `J` on the same
  cut (below). Its objective `w_L · 2ℓ + w_H · H` is then ours with
  `α_L = 2 w_L / (900 w_H)`: `7.4·10⁻⁴` for `w_H = 3`, `2.22·10⁻³` for 1, `2.22·10⁻⁴` for 10
  (`w_L = 1`). Its Dirichlet value is our Dirichlet `J` to within 0.86–1.05 on the same cut, so
  `w_D = 3` corresponds to `α_L = 2/3`.
* Ours: 8 clusters of normals, as `Session` does it (`Start::NormalClusters`), `StopRule` default
  with at most 4000 steps, two stage sets: the regions stages (`FlowOptions::for_regions`, the
  paper's stages with the flux recovery of §3.5, the recommended set for `N > 2`) and the default
  stages (length step on the cut, local extension, mean-gradient speed); each on the caller's
  triangles and on the intrinsic Delaunay triangulation (`Triangulation::default()`, flips only:
  spot 281 flips, bob 657, blub 1 865, the icosphere none; no vertex added). Times are one
  thread, release build, one run at a time; the reference's are wall clock, GUI included, all
  cores available (its solver is multi-threaded).

**Hencky energy, our energy of the final state** (`E = J + α_L ℓ`, ×10⁻³; patches in
parentheses; "regions / default" are the two stage sets; the reference's cut transferred to the
input mesh as above):

| shape, `w_H` | ours 10-13, before D22 | ours 10-13, after D22 | ours now, caller | ours now, intrinsic | reference, 300 steps | reference, run to rest |
|---|---|---|---|---|---|---|
| spot, 3 | 9.675 (42) / – | 18.16 (33) / 9.231 (39) | **8.712** (32) / 9.239 (39) | 8.866 (34) / 9.281 (40) | 9.626 (29) | 8.963 (30) |
| sphere, 3 | 6.072 (8) / – | 6.073 (8) / 6.045 (8) | 6.049 (8) / **6.045** (8) | 6.049 (8) / **6.045** (8) | 7.545 (10) | 7.428 (7) |
| bob, 3 | 7.706 (20) / – | 7.140 (15) / 7.363 (18) | 7.084 (16) / 7.358 (18) | **7.073** (16) / 7.358 (18) | 8.630 (19) | 7.766 (17) |
| blub, 3 | 7.848 (18) / – | 6.663 (10) / 6.985 (14) | **6.293** (12) / 6.983 (13) | **6.293** (13) / 6.974 (14) | 8.376 (14) | 7.022 (11) |
| spot, 1 | | | 21.38 (19) / 21.89 (22) | **19.68** (13) / 22.09 (23) | 24.53 (21) | |
| spot, 10 | | | 3.472 (46) / **3.424** (51) | 3.434 (48) / 3.446 (52) | 17.59 (34) | |

**The same, measured by flattening** (`∫u²` of the flattening with `u = 0` on the cut, plus
`α_L` times the cut length; ×10⁻³; nothing transferred):

| shape, `w_H` | ours, caller: regions / default | ours, intrinsic: regions / default | reference, 300 steps | reference, run to rest |
|---|---|---|---|---|
| spot, 3 | 8.714 / 9.197 | 8.873 / 9.236 | 9.091 | **8.518** |
| sphere, 3 | 6.037 / **6.035** | 6.037 / **6.035** | 6.535 | 6.245 |
| bob, 3 | 7.085 / 7.360 | **7.073** / 7.360 | 7.962 | 7.232 |
| blub, 3 | 6.177 / 6.937 | **6.170** / 6.929 | 6.880 | 6.321 |
| spot, 1 | 21.35 / 21.86 | **19.70** / 22.05 | 22.90 | |
| spot, 10 | 3.488 / **3.434** | 3.451 / 3.444 | 18.26 | |

**Steps and times** (Hencky 3; ours: steps, stop reason, seconds):

| shape | ours caller, regions | ours caller, default | ours intrinsic, regions | ours intrinsic, default | reference |
|---|---|---|---|---|---|
| spot | 159, rose, 2.14 | 104, stalled, 1.14 | 113, rose, 1.52 | 104, stalled, 1.15 | 300 steps in 28 s; at rest after about 5 300 |
| sphere | 318, stalled, 4.11 | 232, stalled, 2.22 | 318, stalled, 4.17 | 232, stalled, 2.20 | 300 in 23 s; about 3 300 |
| bob | 239, stalled, 6.96 | 195, rose, 4.19 | 268, rose, 7.63 | 199, rose, 4.17 | 300 in 37 s; about 4 300 |
| blub | 188, rose, 6.57 | 322, stalled, 8.66 | 220, stalled, 7.77 | 238, rose, 6.19 | 300 in 45 s; about 8 300 |

On 2026-10-13 the same runs took, regions stages before D22: 1161, 484, 1113, 936 steps (13.6,
6.4, 33.0, 33.3 s); after D22, regions / default: 26 / 102 steps on spot (0.38 / 1.13 s), 300 /
232 on the sphere, 293 / 195 on bob, 247 / 407 on blub (3.8 to 10.4 s). "rose" is `EnergyRose`
(D27: the run returns its state of lowest energy).

What the tables say:

* **Against the reference after 300 steps** (the comparison of 2026-10-13) we are below it on
  every shape and weight: by flattening by 4 % (spot) to 14 % (spot, `w_H = 1`), in our energy
  by 9 % to 25 % (inflated by the transfer); at `w_H = 10` its cut on spot has a chart with distortion up to
  `e^1.2` and is 5 times our energy.
* **Against the reference run to rest** (below), by flattening, the fair yardstick, we are below
  it on three shapes by 2–3 % (sphere 3 %, bob 2 %, blub 2 %, our best stage set), and on
  **spot the reference at rest is 2.2 % below our best run** (8.518 against 8.714): its cut is
  longer (9.79 against our 9.60) with less distortion (`1.27·10⁻³` against `1.61·10⁻³`) and
  30 charts against our 32. In our energy we are below it on all four (3 % on spot, 9–19 %
  elsewhere), but there the transfer of its cut costs it more than these differences.
* **D22–D31.** Spot with the regions stages is the run of D26: after D22 it rose to
  `1.82·10⁻²` in 26 steps; now it descends to `8.71·10⁻³` (−52 %) and stops by D27. The other
  regions-stage runs end 0.4–6 % lower than after D22 (sphere 0.4 %, bob 0.8 %, blub 5.6 %); the default
  stages moved in the fourth digit (D30). Against before D22 every run ends lower and needs
  1.5 (sphere) to 7 (spot) times fewer steps and less time.
* **Intrinsic triangulation.** On the icosphere it changes nothing (already Delaunay). On spot,
  bob and blub the final energies on the two triangulations differ by −8 % (spot, `w_H = 1`) to
  +2 % with the regions stages and by less than 1 % with the default stages, in both directions; the intrinsic run evaluated on the caller's triangles
  (its state per input vertex, no vertex was added) differs from its own value by at most 0.3 %.
  It is a different discretization of the same problem, not an improvement on these meshes,
  whose triangles are reasonable (smallest angle after the flips 7–11°). Its time is the same.
* **Speed.** Ours comes to rest in 1–9 s on one thread; the reference takes 23–45 s for 300
  steps and 3 300 to 8 300 steps (4–17 minutes of wall clock, four instances at a time) to come
  to rest.

**Convergence of the reference.** The reference was driven in batches of 1 000 more steps until
the mean of its objective over the last 100 steps changed by less than 0.02 % between batches
(its fixed step makes the objective oscillate from step to step by up to 1 %). Hencky 3:

| shape | objective after 300 steps | at rest | steps | charts 300 → rest | cut length 300 → rest |
|---|---|---|---|---|---|
| spot | 24.89 | 23.23 (−6.6 %) | about 5 300 | 29 → 30 | 10.55 → 9.79 |
| sphere | 17.10 | 16.32 (−4.6 %) | about 3 300 | 14 → 13 | 7.43 → 7.17 |
| bob | 21.45 | 19.33 (−9.9 %) | about 4 300 | 21 → 20 | 9.35 → 8.02 |
| blub | 19.00 | 17.55 (−7.6 %) | about 8 300 | 16 → 12 | 8.19 → 7.49 |

The final cut does change: shorter, with fewer charts on blub, and with the borders moved
(pictures). The first 300 steps of these runs reproduce the runs of 2026-10-13 exactly; with 4
instead of 64 threads they differ from step 43 on in the fifth digit (final objective
24.90 against 24.89, same 29 charts): its parallel reductions are not deterministic in the
thread count.

**The factor between the Hencky values.** Black-box probes, Hencky 3:

* Spot and the icosphere scaled by 2 and by 0.5: every printed value of every step and the saved
  cut are identical to the unscaled runs (300 steps on spot, 1 on the icosphere). The program
  normalizes the scale; the factor is not a unit.
* Spot and the icosphere subdivided 1:4 (same surface, four times the triangles): the Hencky
  value of the start changes by +0.2 % (spot) and −2.7 % (the icosphere, whose clustering start
  differs), after 300 steps on spot by −0.7 %. The value is an integral, not a sum over elements.
* The printed Hencky value divided by `∫u²` of the program's own UV map of the same state
  (12 states: four shapes, three weights, 1 and 300 steps, subdivided): 933 to 1114, median about
  1 010. Divided by our `J` of the transferred cut it is 441 to 1052, low exactly where the
  transfer loses charts (blub 441, the icosphere 531): that was the spread "440 to 1010"
  reported on 2026-10-13, an artifact of the transfer, not of the program.

So its Hencky value is about 10³ times `∫u² dA` on the unit-area surface, the same for every
shape, scale and resolution tried; the spread is the difference between its discretization and
the flattening. Where the factor of about 10³ comes from (a fixed weight or a convention of its
energy) cannot be told from outside. The correspondence above uses 900 (calibrated on one shape
on 2026-10-13); with 1 000 the length weights would be 10 % smaller. Its Dirichlet value has no
such factor.

**Dirichlet energy** (`w_D = 3`, `w_L = 1`, `α_L = 2/3`; our energy of the final state, patches
in parentheses):

| shape | ours, caller: regions / default | ours, intrinsic: regions / default | reference, 300 steps | reference, after 8 300 steps |
|---|---|---|---|---|
| spot | 11.13 (55) / 11.04 (64) | **10.88** (57) / 11.03 (69) | 19.68 (35) | 21.95 (34) |
| sphere | 4.127 (8) / **4.124** (8) | 4.127 (8) / **4.124** (8) | 4.312 (6) | |
| bob | 6.702 (18) / 6.693 (18) | **6.546** (16) / 6.631 (18) | 8.045 (19) | 7.450 (14) |
| blub | **7.735** (38) / 8.501 (40) | 7.853 (39) / 8.060 (38) | 10.32 (18) | |

The reference does not come to rest with the Dirichlet energy and its fixed step: its objective
oscillates over the whole run (spot: between 35 and 88 over 8 300 steps, 55.9 after 300; bob
between 19 and 30), so "after 8 300 steps" is an arbitrary state. Our runs end in 0.4–11 s
(sphere 44–46 steps, bob up to 539). We are below it on every shape, by 4 % (sphere) to 45 %
(spot); the transfer of its cut favours us here too, and with this energy our cuts have many
more patches on spot and blub (the reference's cut on the sphere has 6 charts against our 8).

**Flattening on the same cut.** The program saves a UV map. On its own saved mesh with its
seams, our flattening with `u = 0` on the cut gives the same distortion to within 3.2 % on
every run (14 runs; e.g. spot `1.281·10⁻³` against `1.302·10⁻³`, bob `1.047·10⁻³` against
`1.046·10⁻³`, blub `0.816·10⁻³` against `0.842·10⁻³`, with one flipped triangle of ours on
blub); its UV map is evidently the flattening with the boundary lengths kept, the map whose
distortion the energy measures. Our flattening with optimized boundary values (`BoundaryValues::
Optimized`) has 2 to 6 times less area distortion on the same cut but flips triangles
per atlas (11 to 87 on these 14).

**Not explained or not done.** Where the factor of about 10³ comes from. Why the reference at rest
beats us on spot by flattening; we did not start our flow from its converged cut. The reference's
Dirichlet runs never come to rest with its default step; its line search option and other step
sizes were not tried. Its "use local scale" options, symmetry, area and visibility terms were not
compared. Its times are with a GUI that renders every frame in software.

## Verification

`cargo test -p ddg-variational-cuts` (about half a minute in a debug build; the tables below come from
the ignored tests, see the end).

**Finite differences decide the signs** (`tests/finite_differences.rs`, SPEC §6 item 3). For a
perturbation `φ + εψ` the cut point `c` sweeps the area `dA_c` into region 0, and SPEC §1 predicts
`dE = Σ_c (G_0 − G_1)_c dA_c`. The central difference of the discrete energy agrees with it in
sign and size for every energy, for each region alone, and with `w_r`, `β_r`. Largest difference
relative to `Σ_c |s_c dA_c|` over 10 problems and 2 directions `ψ`, default options | paper options:

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

**Lean oracle** (`tests/oracle.rs`). The committed fixtures of `lean/oracle/fixtures` (`cuts`,
`cuts-default`, `regions`), the output of an independent dense implementation written from
`SPEC.md` alone, are replayed. Of the evaluation: the scaling, face curvatures, crossings, every
triangle of the cut mesh (nodes, lengths, areas), the segments, `u`, `p`, `v`, `ℓ`, `q` at every
node, `J` and the areas per region, the energy, `d`, `e`, `G` and `s` (before and after sliver
smoothing where the fixture has it). Of the step: every stage. The counts and the largest
differences are in the section "Lean oracle" above.

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
| `cut-length-step-on-cut` (new, SPEC D14) | `rate_on_cut` (internal) |
| `cut-extension-local` (new, D15), `cut-harmonic-extension` | `Flow::extend_from_cut`, `Flow::extend_harmonically` (internal) |
| `cut-gradient-norm` (new, D16), `cut-explicit-step` | `gradient_norms` (internal), `explicit_step` |
| `cut-length-step`, `cut-penalty-step-paper` | `Flow::length_step` (internal) |
| `cut-symmetry-projection`, `cut-redistance` | `symmetry_projection`, `init::redistanced` |
| `cut-step-cap` | `StepSize::step`, `rule_step` (internal) |
| `cut-stop-energy`, `cut-run-lowest` (new, D27) | `energy_criterion`, `run_outcome` (internal) |
| `cut-stop-motion` (new, D31) | `motion_criterion` |
| `cut-step-local` (new, D22) | `StepSize::local_factor` |
| `cut-weighted-quantile` (new, D22) | `weighted_quantile` |
| `cut-roughening-rate` (new, D18; a quantile since D22) | `Flow::roughening_rate`, `length_weight_for_ripples` |
| `cut-edge-crossing-regions`, `cut-triple-point` | `Cut::from_regions`, `TriplePoint` |
| `cut-region-extension`, `cut-region-projection` | `boundary_rates`, `redistance_regions` (internal) |
| `cut-length-step-regions` (new, D20) | `region_rates` (internal) |

## Reproducing the tables

```text
cargo test --release -p ddg-variational-cuts --test finite_differences -- --ignored --nocapture
cargo test --release -p ddg-variational-cuts --test analytic -- --ignored --nocapture
cargo run  --release -p ddg-variational-cuts --example estimator_study -- 3 6
cargo run  --release -p ddg-variational-cuts --example flow_study -- circle        # also: latitude-ode,
        # modes, balance, fig1 [LEVEL], fig5 [CELLS HEIGHT X0 MOTION], speeds MESH.obj, timings [MESH.obj],
        # omega MESH.obj [ENERGY] (quantiles of the local roughening rate of a planar cut)
cargo run  --release -p ddg-variational-cuts --example redistance_study
cargo run  --release -p ddg-variational-cuts --example network_study -- honeycomb 84 paper   # also:
        # stripes, atlas MESH.obj [N RIPPLE ENERGY STEPS], starts MESH.obj ..., timings MESH.obj ...
        # (STUDY_ALPHA=A sets the length weight instead of the ripple rule),
        # compare MESH.obj ENERGY auto|strong|ALPHA paper|default   (one line of the table of stage sets)
cargo run  --release -p ddg-variational-cuts --example rise_study -- MESH.obj hencky 0.00074 stages=paper
        # the change of the energy by every stage of a step; see its header for the variants
cargo run  --release -p ddg-variational-cuts --example cut_flow -- MESH.obj --energy hencky --alpha 2
cargo run  --release -p ddg-variational-cuts --example bad_mesh_study      # see its header for the parts
cargo test --release -p ddg-variational-cuts --test regions -- --nocapture
cargo test --release -p ddg-variational-cuts --test flow -- --nocapture
cargo test --release -p ddg-variational-cuts --test oracle -- --nocapture
cargo test --release -p ddg-variational-cuts --test oracle -- --ignored --nocapture   # needs lean/oracle built
```
