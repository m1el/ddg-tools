# ddg-repulsive-surfaces

A Rust implementation of

> Chris Yu, Caleb Brakensiek, Henrik Schumacher, Keenan Crane. *Repulsive Surfaces.* ACM
> Transactions on Graphics 40(6), Article 268, 2021. <https://doi.org/10.1145/3478513.3480521>

the tangent-point energy of triangle meshes in R³ (orientation not required: Klein bottles and
Möbius bands are fine), its descent in the fractional Sobolev metric of the paper, constraints and
penalties, obstacles, Barnes–Hut and hierarchical-matrix acceleration with preconditioned GMRES, and
dynamic remeshing. The contract is [`SPEC.md`](SPEC.md): every discrete formula with its `rs-` id,
every place where the paper or the authors' code is wrong or silent (decisions D1–D33, §11), the
API, the measurements (§14) and the open items (§16). The `NOTES-*.md` files are the work-package
logs; [`COMPARISON.md`](COMPARISON.md) is the black-box comparison with the authors' program.

## What is in it

| part | paper | module |
|---|---|---|
| `Surface`: validated manifold triangle mesh, orientation optional, boundary loops, components; OBJ I/O | Sec. 2 | `surface`, `io` |
| Face geometry and derivatives, cotan stiffness `W`, lumped mass `M` | Sec. 3 | `geometry` |
| The kernel `K(X, N; Y)` and the exact face-pair energy with its exact differential (parallel, bitwise independent of the thread count); Coulomb and Willmore energies | Eq. 9, 10, 16 | `kernel`, `energy` |
| Barnes–Hut energy with its exact adjoint gradient on a BVH | Sec. 4, Eq. 17 | `spatial`, `energy::barnes_hut` |
| The fractional Sobolev metric `A = B + B⁰` (dense), and the L², H¹, H² comparison metrics | Eq. 11–13, App. A | `metric` |
| Hierarchical matrices: block cluster tree, `A′` product, preconditioner | Sec. 5, App. B | `metric::hierarchical`, `solve::iterative` |
| Dense saddle with Schur complement, GMRES (right or left), fallbacks | Sec. 5, Eq. 14, 18–20 | `solve` |
| Constraints: barycenter, area, volume (growth schedules), pins and moving pins; projection | Sec. 7, Eq. 15, 21 | `constraints` |
| Penalties and obstacles: mesh, point cloud, implicit (incl. `Cylinder`), attractors | Sec. 7, 9 | `penalties`, `obstacle`, `implicit` |
| Remeshing: split, collapse, flip, tangential smoothing, pins untouched, vertex and face maps | Sec. 6 | `remesh` |
| Continuous collision detection (suggested in Sec. 10), also for inputs that already intersect | Sec. 10 | `ccd` |
| The flow: `Session` (step by step, cancellable), `flow`, line searches, schedules, stop rules | Sec. 4–7 | `step` |
| Comparison methods: AQP, H¹ L-BFGS, BCQN (reference variants and ours), soft problems | Sec. 8.2 | `optimize` |
| Generators (genus-g, tubes around curve networks, text, punctured torus, cobblestones, soups) and a reader of the authors' scene files | Sec. 9 | `generate`, `scene` |

### Departures from the paper and from the reference

Each is an option; `reference_compat: true` switches the reference's choices back on (used for the
black-box comparison only). Full table: SPEC §11.

* **D2, `D_f` misprint.** The paper says "unit edge vectors"; plain edge vectors are correct.
* **D4, exact Barnes–Hut gradient.** The paper differentiates the forward terms only; we
  differentiate the approximated energy exactly, so gradient and line-search energy agree.
* **D9, volume gradient.** The reference's volume row is `Σ aN`, three times the gradient; ours is
  exact (`⅓ Σ aN`). The parity tests show exactly this 3×.
* **D20, CCD on by default** (0.5–2 % of a step). The paper's algorithm without CCD tunnels on
  EntangledDoubleTorus: in step 6 it accepts `τ = 3.72` and an exact swept test finds 2226 crossings
  (`NOTES-CCDINV.md`). With CCD all 200 steps are crossing-free. Inputs that already intersect
  (Klein bottle) get a band rule.
* **D32, vertex-pair preconditioner.** The printed `L_{2−s}` is the face form, whose averaging nearly
  annihilates high frequencies: condition numbers 12.6 / 68 / 224 against `M Δ^{2−s}` (80 / 320 /
  1280 faces); GMRES to 1e-4 needs 38 / 74 / > 200 iterations and stalls at 4e-2 on a bumpy disk. The
  reference's left GMRES stops after 13–19 iterations with a true residual of 0.1–0.2. The vertex-pair
  form has condition 1.3–2.1 and needs 5–27 iterations. The face form stays as a mode.
* **D33, augmented preconditioner.** As printed, `L′⁻¹ blockdiag(L_{2−s}, 0) L′⁻¹` is singular and
  cannot represent the multiplier part: with pins or per-component barycenters right GMRES stalls
  (Hopf link 8.6e-4 after 200 iterations; 3 pins on a torus 0.12). We use
  `L′⁻¹ [[L_{2−s}, C_Qᵀ],[C_Q, 0]] L′⁻¹` at the same cost (proved, `RS.aug_precond_props`). Also from
  `NOTES-SOLVERFIX.md`: the preconditioner scale `L_{2−s}/(w_TP · 4 I(σ) I(2−s))` (pinned runs stalled
  at 0.965, now converge) and a coarse correction on per-component translations (cobblestones 2800 →
  134 iterations per step). A diverged solve falls back to augmented right GMRES, then a dense LU.
* **D8, barycenter `Auto`.** Global with the low-order term on; per component with a mean-motion
  shift when `low_order: false`; weights `a_v/A_c`.
* **Remeshing (D16, D17).** Order flip → smooth; collapse, flip and smoothing checks in f64; a growth
  cap; pins and per-vertex data transferred, pinned vertices never touched.
* **D24, penalty gradients.** The reference has several inconsistent ones; ours are exact gradients
  of the printed values.
* **Other:** Willmore `H²` term in `A` only (D7); Newton correction until tolerance (D12); mass-matrix
  L² metric (D25); standard H¹ L-BFGS two-loop (D26); BCQN with a spectral-norm estimate (D27); an
  energy-stall stop rule (D29). Defaults moved by measurement (WP17): GMRES tolerance 1e-3 (reference
  1e-4), restart 60, `dense_max_faces = 1000`.
* **The reference's Klein bottle limitation.** Its geometry-central cannot load non-orientable meshes
  ("duplicate edge in list"), so Fig. 26 (bottom) cannot be made with the published `master`. The
  authors' unmerged branch `nonorientable` (commit `ce89c2c`) uses the general `SurfaceMesh` and
  switches off splits and collapses. Our `Surface` accepts them, with full remeshing.

## Use

### Rust

```rust
use ddg_repulsive_surfaces::{FlowOptions, flow, read_obj};

let (surface, positions) = read_obj("torus.obj")?;
let options = FlowOptions { max_iterations: 200, ..FlowOptions::default() };
let result = flow(surface, positions, options, None)?;
println!("{:?}: {} steps", result.report.status, result.report.iterations);
// result.surface and result.positions: the mesh after remeshing.
```

`Session::new(surface, positions, options)` and `Session::step` give one step at a time (cancellable;
the surface changes with remeshing, `face_origin` maps faces to the input's). `run_scene(path,
&RunSceneOptions::default())` runs the authors' scene files (`examples/run_scene.rs`; the scenes are
copied into `data/repulsive-surfaces/` by `scripts/rs-fetch-scenes.sh`). Each application of Sec. 9
is an example in `examples/apps/`.

### Python

```python
import ddg_tools
s = ddg_tools.RepulsiveSurfaceSession(positions, faces, pins=[0, 1], remesh=True)
report = s.step()
snap = s.snapshot()           # positions, faces (change with remeshing), pinned, face_origin
positions, faces, report = ddg_tools.repulsive_surface_flow(positions, faces, steps=100)
energy, differential, report = ddg_tools.surface_energy(positions, faces)
```

Also `repulsive_surface_defaults()` and `repulsive_surface_scene(path)`; see `crates/ddg-py/README.md`.

### Blender

"Repulsive Surface Flow (Live)" (`object.ddg_repulsive_surface_live`: modal, worker thread,
Accept/Cancel) and "Repulsive Surface Flow (Steps)" (`object.ddg_repulsive_surface_flow`: N steps, one
undo step), panel in the DDG sidebar tab. Pins from a vertex group, obstacles from a collection,
implicit obstacles and attractors from empties, barycenter/area/volume constraints, penalties. With
remeshing the mesh data is replaced (UVs and colours are dropped with a warning). The panel defaults
to an area constraint, because a free flow inflates the surface. See `blender/README.md`.

## Results of record

* **Consistency (Fig. 10; `NOTES-WP16A.md`, `examples/consistency.rs`).** Perturbed torus
  (R = 1, r = 1/3), smooth energy by our own cubature, ℰ(6,12) = 460.71622. Discrete energy on 648 to
  85,554 faces: fitted rate 1.97 (1.86 over the last four points), the `O(h²)` end of the paper's
  "between `O(h)` and `O(h²)`". Barnes–Hut error `~ θ²` for `θ ≥ 0.35`; for `θ ≤ 0.25` the signed
  error changes sign (`|error| ≤ 2e-4` at `θ = 0.125`), not investigated further.
* **Performance (Fig. 13; `NOTES-WP17.md`, `benches/`).** At 4 threads 12 of the 14 models are at or
  below the paper's per-step time; above it: #3 broken figure-8 0.50 s against ≈ 0.37 s, #14 figure-8
  tunnel 3.15 s against ≈ 2.7 s. **The paper's 2.5 s at 45k faces is not met** (outside machine load
  6–16 during the timings). Thread scaling on a 40k torus: 3.32 / 2.08 / 1.44 / 1.14 / 0.97 / 0.93 /
  1.03 s at 1–64 threads. Against the reference at 1 thread: faster on 5 of 8 scenes (walnut 140
  against 259 ms, bunny 3522 against 4205, reconstruction 7638 against 9764, lung, blocky Klein
  bottle), slower on EDT 1.09×, handcuffs 1.18×, Hopf link 1.45× (30/27/24 GMRES iterations against
  its 21/13/9, which stops on the preconditioned residual).
* **Black-box comparison with the authors' program (`COMPARISON.md`).** Parity: energy 2e-16 to 7e-16,
  differential ≤ 1.1e-14, dense constrained metric ≤ 1.1e-14, Barnes–Hut at `θ = 0.5` 3e-6 to 3e-5
  (different splits); the volume row differs by exactly 3× (D9). With `reference_compat` our
  trajectories track its energy to ≈ 1 % and the final shape to ≈ 0.5 % on seven scenes; our defaults
  reach the plateaus sooner (EDT `F` 8.3 against 19.9 at 300 steps). Borromean rings: the reference
  stalls at 7302 and self-intersects in steps 1–9, ours reaches ≈ 4540 without intersection. Klein
  bottles: the published `master` loads none; its non-orientable variant flows the blocky one (ours
  within 0.2 % with splits and collapses off) and stalls on `klein_bottle.obj`.
* **Flow fidelity (`NOTES-FIDELITY.md`).** A torus of revolution stays one within 0.6 % (reference
  0.78 %). Genus-g: 800 steps are too few; g = 3 converges in about 2080 steps to E 39.91 (reference
  at 1000 steps 39.54), g = 2 to 38.59 (reference 38.40). Handcuffs unlink in both programs without
  the plane attractor; a strong plane keeps them linked (use weight ≤ 0.5). Remeshing every step sets
  an energy noise floor; polishing without remeshing lowers `E` by 2–6 %.
* **Method comparison (Fig. 5, 12; Sec. 8.2).** (results pending: WP16 method comparison; Fig. 16 transcription)
* **Knotted genus-2 table (Fig. 16, Ishii et al. 2012).** (results pending: WP16 method comparison; Fig. 16 transcription)

## Verification

* **Finite differences** of the exact energy, the Barnes–Hut energy (against its own gradient with
  frozen lists), every penalty and constraint, the Willmore term; relative error < 1e-6.
* **Invariances and operator checks:** scaling, kernel of `A`, hierarchical against dense products
  (1e-12 at `χ = 0`), dense Schur against the full saddle LU, quadratic Newton correction.
* **Lean proofs** (`lean/proofs/DdgProofs/RepulsiveSurfaces/`): 35 `rs-` ids proved and 8 cited at
  WP10, then the augmented preconditioner (`rs-fractional-preconditioner`, WP10b); the cited ones are
  in `lean/proofs/CITED.md`. `scripts/check-formulas.sh` cross-checks the `[formula: rs-…]` tags.
* **Oracle** (`lean/oracle/DdgOracle/RepulsiveSurfaces*`): an independent Float implementation of
  energy, `D_f`, metric, constraints, saddle, preconditioner, one dense step and the remeshing
  operations; 70 golden fixtures replayed by `tests/oracle.rs`, plus a differential fuzz test (ignored
  by default).
* **Parity fixtures with the authors' program** (`tests/reference_fixtures.rs`,
  `tests/fixtures/reference/`): Gram matrices, Jacobians and one direction on small meshes.
* **Applications:** `tests/applications.rs` runs each application at reduced size.

## Open items

SPEC §16 is the list. The 2.5 s per step at 45k faces is not met; the obstacle part of the
preconditioner does not help (option, off), so envelopes converge slowly; Klein bottle remeshing
inside the CCD band fails at step ≈ 184; the punctured torus diverges after 38–80 steps; open surfaces
without remeshing need 225–678 GMRES iterations; needle inputs of the subcritical spindle; eversion
(Fig. 20) is not reproduced (it needs a gradually raised volume weight and a symmetry search);
`reference_compat` cannot reproduce D9, the absolute `ε` and the penalty value conventions; the final
self-energy of the reference's reconstruction scene (4267 against 14.9) is unexplained; Fig. 16 is
being transcribed.

## Credits and licence

The paper: Chris Yu, Caleb Brakensiek, Henrik Schumacher, Keenan Crane. The authors' reference
implementation (`repulsive-surfaces`, MIT, "Copyright (c) 2019 Nick Sharp") may be read and ported;
its licence notice is kept in [`LICENSES/repulsive-surfaces-MIT.txt`](LICENSES/repulsive-surfaces-MIT.txt).
The crate depends on `ddg-repulsive-curves` (Yu, Schumacher, Crane, *Repulsive Curves*, ACM TOG 40(2),
2021) for the bounding volume hierarchy, block cluster tree, L-BFGS, implicit surfaces and CCD
primitives. Further works: Ishii et al., 2012, the table of knotted genus-2 surfaces used for Fig. 16
(the PDF was supplied by the user and is not in the repository); Kovalsky et al., 2016 (the AQP
method); Zhu et al., 2018 (BCQN); Chen and Holst, 2011 (mesh optimization). This crate is MIT
licensed, as the workspace.
