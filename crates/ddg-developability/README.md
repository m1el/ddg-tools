# ddg-developability

A clean-room Rust implementation, from the text of the paper only, of

> Oded Stein, Eitan Grinspun, Keenan Crane. *Developability of Triangle Meshes.*
> ACM Transactions on Graphics 37(4), Article 77, 2018. <https://doi.org/10.1145/3197517.3201303>

A vertex of a triangle mesh is *discrete developable* (a hinge) when its triangles split into two
flat, edge-connected pieces; equivalently, when its triangle normals lie on one great circle of the
Gauss sphere. The paper measures how far a vertex is from that state by six energies, minimizes
their sum over the vertex positions with L-BFGS, and cuts the result open along the seams that
remain so that it flattens with little distortion. The discrete method as implemented here, with
every formula, its derivation and its sign conventions, is specified in [`SPEC.md`](SPEC.md); every
place where it differs from the paper as printed is listed in **SPEC §8** (D1 to D33). The
authors' code was not consulted (it is GPL; this workspace is MIT).

## Credits and clean room

* The method, the six energies, the flow and the cut are the paper's (Stein, Grinspun, Crane 2018).
* The line search is Algorithm 2.6 of A. S. Lewis, M. L. Overton, *Nonsmooth optimization via
  quasi-Newton methods*, Math. Program. 141 (2013) (weak Wolfe conditions by bracketing and
  bisection, Sec. 5.1 for the setup), and L-BFGS is the two-loop recursion of J. Nocedal, S. J.
  Wright, *Numerical Optimization* (2nd ed., 2006), Alg. 7.4 and 7.5; Alg. 3.1 for backtracking.
* The 3x3 eigenproblem of `E_λ` is solved by cyclic Jacobi rotations (`eig.rs`, a textbook
  method, written from Golub & Van Loan). The paper uses J. Kopp, *Efficient numerical
  diagonalization of hermitian 3x3 matrices*, Int. J. Mod. Phys. C 19 (2008) 523 (arXiv
  physics/0610206); **his paper is cited and nothing else of his is used: his C code, and the
  authors' implementation which contains it, were never opened or copied.**
* The cut graph is the cut planner of `ddg-bff` (`plan_cut`): K. Mehlhorn's approximation of the
  minimum Steiner tree (*A faster approximation algorithm for the Steiner problem in graphs*,
  1988) through the terminals, tree–cotree handle loops by the greedy loops of J. Erickson and
  K. Whittlesey (SODA 2005), and cuts between boundary loops. The paper itself uses the cut of
  J. Erickson and S. Har-Peled (*Optimally cutting a surface into a disk*, 2004); both
  approximate an NP-hard problem and for genus 0 both reduce to a tree through the terminals
  (SPEC D22).
* Each chart is flattened by **Boundary First Flattening**, R. Sawhney and K. Crane, ACM TOG 36(4),
  2017 (crate `ddg-bff`), with zero boundary scale factors, which keeps boundary lengths.
* The graph-theoretic facts used by the remeshing legality rules are cited in SPEC §9: the flip of
  an interior edge (Edelsbrunner 2001) and the link condition of a collapse (Dey, Edelsbrunner,
  Guha, Nekhayev 1999).

**Clean room.** `SPEC.md` was written from the text of the paper, with the formulas transcribed
from the page images; all formulas were re-derived and checked against finite differences, and the
places where the printed ones are wrong are recorded. The authors' implementation (GPL, with Kopp's
LGPL eigensolver) is quarantined outside this repository and was never opened, searched, fetched
or paraphrased by anyone who worked on this crate. Papers the text cites may be read.

## What is implemented

| stage | SPEC | what |
|---|---|---|
| energies | §1–§3 | `Energy::{Combinatorial (E_P, Eq. 1–2), CombinatorialMax (π_max, Eq. 5), Covariance (E_λ, Eq. 3–4), CovarianceMax (λ_max, Eq. 6), Intrinsic (App. B.5), IntrinsicMax (App. B.5.1)}`, per vertex, with analytic gradients through one adjoint pass (Eq. 7–9), parallel and bitwise independent of the thread count |
| optimizer | §4 | L-BFGS, gradient descent (weak Wolfe or backtracking), fixed step; units of the mean edge length; step cap; pinned vertices; admissibility floor on the smallest corner angle (D31) with blocked vertices |
| subdivision | §5.2 | 4-1 midpoint subdivision; energy-preserving; pinned ids and parent faces tracked |
| remeshing | §5.3 | small-angle repair by edge flips (caps) and collapses (needles), all legality rules, every refusal named |
| flow | §5.1 | rounds of subdivide, remesh, optimize with the remeshing checked every iteration; `vertex_source` and **`face_parent`** (the exact input face of every output face); stops cleanly before a round whose energy is undefined (D32) |
| cut and flatten | §6 | terminals from the energy (`ε`) and from the angle defect (D30), `plan_cut`, BFF per chart; per-chart distortion, flips and overlap |

Not in the paper and not here: gluing small cracks after the cut and avoiding global overlap of a
chart with itself (the paper does both by hand; SPEC D23). The overlap is measured and reported.

## Use

```rust
use ddg_core::TriMesh;
use ddg_developability::{CutOptions, EnergyOptions, FlowOptions, cut_and_flatten, evaluate, flow};

// `mesh: TriMesh` (manifold, consistently oriented) and `positions: Vec<[f64; 3]>`.
let before = evaluate(&mesh, &positions, &EnergyOptions::default())?;   // energy, vertex_energy, gradient
// before.vertex_energy is what Show Developability paints.

// Two subdivisions (three rounds) of 500 L-BFGS iterations of E_λ, remeshing on.
let options = FlowOptions::uniform(EnergyOptions::default(), 2, 500);
let result = flow(&mesh, &positions, &options, None)?;
result.mesh;                       // the final TriMesh
result.positions;                  // its positions
result.vertex_source;              // input vertex of every output vertex (NONE: created by subdivision)
result.face_parent;                // input face of every output face (exact through subdivision, flips, collapses)
result.report.rounds[0].optimize;  // status, iterations, energies, histories, inadmissible_trials, ...
result.report.stopped;             // Some(FlowStop { round, after_subdivision, vertex, reason }): D32
result.report.bad_faces;           // faces left below the small angle
result.report.min_angle;           // smallest corner angle of the result (>= min_face_angle)

// Cut through the vertices of high energy / large angle defect and flatten every chart with BFF.
let flat = cut_and_flatten(&result.mesh, &result.positions, &CutOptions::default())?;
flat.cut_edges;                    // edge ids, ascending
flat.corner_uvs;                   // 3 per face, charts not packed
flat.charts;                       // per chart: flips, overlap, quasi-conformal and log-scale distortion
flat.vertex_energy; flat.vertex_defect; flat.star_embeddedness;
```

The API in brief (SPEC §7.1 has every signature):

* `evaluate(mesh, positions, &EnergyOptions) -> Evaluation` (value, per-vertex values, gradient,
  per-vertex margins of the active piece, near-tie count) and `evaluate_energy` without the gradient.
  `EnergyOptions { energy, include_boundary, omit_valence_three, threads }`.
* `optimize(mesh, positions, &EnergyOptions, &OptimizeOptions, progress)`, `subdivide`, `remesh`,
  `flow(mesh, positions, &FlowOptions, progress)`, `cut_and_flatten(mesh, positions, &CutOptions)`.
  `progress` is a callback that may return `Control::Break`; the state at that moment is returned
  with `cancelled` set.
* Option structs are plain data with `Default` and `validate()`; nothing is clamped or rescaled
  (positions are never rescaled; tolerances are in units of the mean edge length).
* Every result has a report. Bad input is a `DevError` whose message names the face or vertex
  (`NearDegenerateFace { face, min_angle, limit }`, `UndefinedEnergy { vertex, reason }`, ...);
  nothing panics.
* Nothing mutates its input. A polygon mesh is the caller's business: `ddg-py` and the Blender
  add-on triangulate with the caller's triangulation and map `face_parent` back to polygons.

Command line: `cargo run --release -p ddg-developability --example develop -- IN.obj [options]`
(OBJ in; flow; OBJ and per-vertex energy out; optionally the cut and the flattening with UVs).
Python: `ddg_tools.developability`, `make_developable`, `developable_seams`, `developable_unwrap`
(`crates/ddg-py/README.md`). Blender: Make Developable, Show Developability, Developable Seams + UV
(`blender/README.md`, with end-to-end numbers on spot, bob and blub).

## Where this departs from the paper

Full text, with the evidence, in **SPEC §8**. The ones that matter:

* **Misprints in the paper's formulas** (all found by finite differences, each with the printed
  and the corrected form): D1 (Eq. 1, ordered pairs count twice), **D2** (Eq. 8: `1/(2A)`, not
  `1/A`), **D3** (Eq. 9: divide by `|f_i − f_j|²`), D4, **D5** (App. B.4: sign of antipodal sites;
  the printed gradient has relative error about 2 there), **D6** (Eq. 11: `μ_i μ_iᵀ`, and `v_i^jk`
  is the binormal), **D7** (App. B.5.1: the altitude gradient is the negated transpose; the `±`
  enumeration is four times too large).
* **D8, D9–D12, D14**: which vertices count (boundary included, valence 3 omitted, both options),
  bipartitions, ties and repeated eigenvalues (a deterministic active piece, near ties reported),
  `λ_max` candidates, and the intrinsic energies being undefined (an error, `+∞` in a line search).
* **D15–D19, D26, D27**: optimizer defaults the paper leaves open (L-BFGS memory 10, `c1 = 0`,
  `c2 = 0.5`, units of the mean edge length, step cap 0.5, stopping rules).
* **D21**: the remeshing thresholds and legality rules (the paper gives neither); 1 degree.
* **D22, D23**: the cut. Terminals from energies without the valence-3 omission, `plan_cut`
  instead of Erickson–Har-Peled, terminals added on closed genus-0 pieces, the edge weight; the
  charts are not glued and not made overlap free. The UV map of the bindings comes from the atlas
  pipeline of `ddg-uv` with the developable cut as seams, because the plain charts overlap
  themselves (spot after the default flow: one chart, half of it covered twice).
* **D30: defect terminals.** The paper cuts through the vertices with energy above `ε`. The flow
  leaves stars that are not embedded (spikes, double covers; D28) whose normals are coplanar, so
  their energy is about 0, but whose angle defect is large; flattening absorbs that as area
  distortion and flips triangles (1.9 % of the UV triangles on spot's control mesh, 5.7 % when the
  flow is run long). Every interior vertex with `|Ω| > 0.1` rad is a terminal as well
  (`CutOptions::defect_tolerance`, `None` restores the paper's rule): 0 flips in all 6 control-mesh
  cases, free on dense meshes, for 1.1 to 5 times the cut length of the paper's rule (SPEC §12.3).
* **D31: admissibility floor.** The energies are scale free, so squeezing a face to zero area costs
  almost nothing and the weak Wolfe search accepts it: 40 of 48 small bare optimizer runs ended
  with a face below `3e-9` rad. A step to a face below `min_face_angle` (1e-4 rad) is an infinite
  energy, an input below it is an error, a collapse that would create one is refused, and the
  vertices of faces within a factor 2 of the floor are blocked so that one squeezed face does not
  stop the rest.
* **D32: the flow does not walk into an undefined energy.** The state at the end of a round can be
  a flat pancake whose subdivision has a vertex where the intrinsic energy is undefined. The flow
  checks the next round's starting state first and, if it is undefined, stops, returns the state
  before and says so (`FlowReport::stopped`). A hard wall on the fold angle was tried at four
  values and rejected (SPEC §12.4). A remeshing result that is undefined or below the floor is not
  applied (`RemeshReport::rejected`).
* **D33: two-face boundary stars** contribute exactly nothing (the eigenvector of a rank-deficient
  matrix is arbitrary and had fed noise into the gradient; found by the oracle fuzz).
* **D24**: the polygon rule of the bindings (topology unchanged: positions only; otherwise the
  triangle mesh with the material of the parent polygon, which `face_parent` makes exact).

## Measured

All numbers are from this tree; the studies and their tables are in SPEC §12 (the logs under
`runs_tmp/` are not committed; the scripts in `examples/` regenerate them). Timings were taken on
a 64-core machine loaded by other jobs and are upper bounds.

**Gradients (SPEC §10.2).** Central differences, `h = 1e-6 ℓ̄`, against the analytic gradient, for
all six energies on random stars of valence 3 to 9, boundary fans, an icosphere, a bumpy grid, a
torus, the crumpled model and spot's control mesh, only where the active piece has a margin:
relative error 1e-10 to 2e-8. The building blocks (Eq. 7–9, the adjoint, the log map, the
circumcenter and altitude derivatives) have their own FD tests.

**Exact values (SPEC §10.3).** Hinge and flat stars: all six energies below 1e-28. Closed forms of
`E_λ`, `E_P`, `λ_max` (even valence) and both intrinsic energies on the regular cone of valence 4 to 8
and three heights, agreeing to 1e-12.

**Lean oracle** (`lean/oracle`, an independent Float implementation written from `SPEC.md`
without this crate's code; `tests/oracle.rs` replays golden fixtures and fresh seeded cases).
Seed 3, 400 000 random cases, 209.7 s on 64 threads: about 3 million vertices per energy compared,
energies agree to `4e-16` (`E_P`), `4e-15` (`E_λ`) up to `2.2e-13` (intrinsic max) absolutely, and to
`1e-12` relative where the energy exceeds 1e-3; gradients to `1e-9` of the largest entry
(worst `5.9e-11` for the intrinsic energy, `1.5e-11` for `λ_max`); 108 609 first L-BFGS steps
(scaled steepest descent plus the Lewis–Overton search) agree in the trial sequence and the ending
(largest position difference `6e-13` mean edge lengths); 39 976 subdivisions identical.
0 disagreements. The fixtures also pin the admissibility floor, the blocked vertices and the start
error (D31) and the two-face boundary stars (D33).

**Flow and cut.** The numbers behind D30–D32 (SPEC §12):
spot, bob and blub control meshes with two subdivisions and 300 iterations per round, paper rule
against the default: spot 99 of 5112 UV triangles flipped (1.9 %), mean |u| 0.90, cut 2.1 % of the
edge length, against 0 flipped, mean |u| 4.0e-3, cut 4.1 %; blub 4.9 % against 0 flipped (cut 0.8 %
against 3.7 %); over-run to 3 x 3000 iterations the paper's rule flips up to 5.7 %. On the dense
meshes (23k to 57k faces) the defect criterion changes nothing. Over-runs of all six energies on
the three control meshes (3 x 3000 iterations) complete without error and without a face below the
floor (18 of 18; two flows stop before a later round by D32). The study of every model in
`data/models` with six energies and two budgets (132 jobs) finished without error, non-finite
value, line-search failure or `NonDescent`. On the crumpled developable (the paper's Fig. 4 input),
`E_λ` falls from 61.9 to 9.1e-4 in 500 iterations and the dihedral angles above 30 degrees from 354 to
35 of 800 edges. The Blender end-to-end table on spot, bob, blub and primitives is in
`blender/README.md`.

**Speed** (taken under load 36 to 150, single thread; indicative only): value and gradient cost
about 1 to 1.8 µs per vertex for `E_P`, `π_max`, `E_λ` and the intrinsic energies and about 4 µs
for `λ_max` (spot, SPEC §12.7). The flow on spot's control mesh with two subdivisions (5112 faces,
900 iterations, `E_λ`) takes 1.7 s; one iteration on a 10.7k-face mesh about 8 ms. Dense meshes of
23k to 57k faces take 60 to 190 s per 800 iterations. The paper reports "seconds to a few minutes"
on a 2018 laptop. Multi-threading does not help at these sizes under load; the 94k-face case gains
3 to 4 times with all threads.

## Verification layers

There is no certifiable minimizer. The layers (SPEC §10):

1. Finite differences of every energy and of every building block.
2. Exact values: zeros on hinges, closed forms on regular cones.
3. Invariances: rigid motions, scaling, retessellation, 4-1 subdivision (all six energies equal to
   1e-12), `Σ g = 0`, `Σ f×g = 0`, `Σ f·g = 0`, bitwise determinism over 1, 2, 7, 64 threads.
4. The optimizer: monotone energy histories, weak Wolfe conditions on every accepted step, no stored
   pair with `sᵀy <= 0`, `|x|` behaving as Lewis–Overton describe, superlinear convergence on a
   quadratic, pinned vertices fixed, cancellation.
5. Topology: Euler characteristic and boundary preserved by subdivision, every remeshing refusal
   hit by a dedicated test, `face_parent` and `vertex_source` consistent.
6. The Lean oracle: differential fuzzing and golden fixtures (above).
7. The Lean proofs (`lean/proofs`): 43 `dev-` formula ids, 30 machine-checked and 13 cited with the
   reason (the paper's smooth theory, Lewis–Overton's termination theorems, combinatorial
   topology); `scripts/check-formulas.sh` cross-checks the tags.
8. End to end on `data/models` (crumpled developable, spot, bob and blub, control and dense meshes):
   energy decrease, charts are disks, handle loops are cut, no flipped UV triangle, and the
   robustness tests of D31–D33 (`tests/robust.rs`).

```text
cargo test -p ddg-developability                                       # fast layers
cargo test --release -p ddg-developability -- --ignored --nocapture     # models, over-runs, long studies
cargo test --release -p ddg-developability --test oracle -- --ignored --nocapture   # needs `lake build` in lean/oracle
cargo bench -p ddg-developability                                       # eig, energy, flow
```

## Known gaps

* The cut is an approximation of the shortest cut and the charts are not glued or made overlap
  free (D22, D23): the paper's own flattening of spot overlaps itself in half of the chart. The
  bindings take the cut as seams and unwrap with the `ddg-uv` atlas pipeline; whether a converged
  flow changes that is measured in `blender/README.md`.
* The minimizers of the energies include stars that are not embedded (spikes, double covers; D28,
  paper Fig. 9); the energies cannot see them. The cut finds them by the angle defect (D30); the
  flow does not avoid them, and a very long flow produces them.
* A flow can end before a round (D32) or leave faces below the small angle (`bad_faces`): both are
  reported, neither is repaired.
* `E_λ` is nonsmooth at ties and repeated eigenvalues; a deterministic active piece is used and
  near ties are counted. L-BFGS on such an energy has no convergence theorem (Lewis–Overton
  observe it works).
* The intrinsic energies are undefined at a vertex whose face normal is antipodal to its vertex
  normal; with `omit_valence_three = false` the cut stage propagates that error.
* Timings on a quiet machine, and the speed-up with threads, are not measured yet.
* Not implemented (not in the paper's method as we read it): remeshing by anything but flips and
  collapses of small-angle triangles, anisotropic remeshing, constraints other than pinned
  vertices.

MIT licensed, like the rest of the workspace.
