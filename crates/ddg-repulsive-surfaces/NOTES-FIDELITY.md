# Flow fidelity: do we reach the paper's minimizers? (Sec. 9.1, Fig. 1/14, 15, 20)

Tool: `examples/fidelity.rs` (`gen`, `flow`, `measure`, `scan`, `symmetrize`; see its header). Build:
`CARGO_TARGET_DIR=$PWD/target-fidelity cargo build --release -p ddg-repulsive-surfaces --example fidelity`.
Every `flow` run: area + barycenter constraints (as `apps`), `StopRule::Off`, defaults unless stated
(Barnes–Hut θ 0.5, `Backend::Auto` = hierarchical above 1000 faces, remeshing, CCD, `GrowPrevious ×2` start).
Indicators: exact all-pairs energy, eigenvalues of the area-weighted second-moment tensor about the area
centroid (`aniso_top = (e2 − e1)/e2`, 0 for a torus of revolution), angular Fourier amplitudes of ρ and z about
the least-moment axis. Reference: `RepulsiveSurfaces-build/build/bin/rsurfaces_headless` on the same OBJ with a
scene `constrain barycenter / constrain area` (+ `implicit plane attract 2 5 0 0 0 0 0 1`), `--threads 1`.
Logs, CSVs, OBJs: `runs_tmp/fidelity/` (not committed). Machine load 55–120 from other jobs throughout;
no timings are quoted.

## Bugs found and fixed

1. **`apps::eig3` (examples/apps/apps.rs) did not diagonalize.** The Jacobi angle had its `atan2`
   arguments swapped, so for non-diagonal input the returned "eigenvalues" were the diagonal of a
   partially rotated matrix (random 3×3: off-diagonal 0.1–2 left). The symmetry table of `NOTES-WP15.md`
   and the control-torus "drift" (moments 0.278/0.422/0.554) came from it. Corrected values for the WP15
   finals are below. Fix: `0.5 * (2 a_pq).atan2(a_qq − a_pp)`; test `tests/applications.rs::
   eig3_and_inertia_of_rotated_input` (Q diag(1,2,5) Qᵀ to 1e-12; a tilted torus has multiplicity 2).
2. **`ccd::intersecting_pairs` false positive on coplanar neighbours** (`src/ccd.rs`,
   `segment_meets_triangle`). Points in the plane of a generic triangle give `orient = ±1e-20`, not 0;
   the coplanar branch needed exact zeros, so the sign test ran on rounding noise and reported two
   coplanar faces with one common vertex as intersecting. The `apps genus02` start mesh (flat sides of
   the rounded slab) had one such pair (faces 1795/1988), so the CCD intersection band was active on an
   embedded mesh (113 exempt pairs in the first steps of every genus-2 run). Fix: treat
   `|orient| ≤ 1e-12 |n| ℓ` for both endpoints as coplanar. Test `tests/ccd.rs::
   coplanar_neighbours_do_not_intersect` (genus-2 start mesh, 50 flat grids in random planes: 0 pairs; a
   genuine coplanar overlap is still found). `tests/ccd.rs`, `flow`, `applications` pass.

## 1. Torus (g = 1): no drift; ≤ 1 % anisotropy, the same as the reference

With the corrected eigen-solver the WP15 control torus (64×32, R/r = 2.5, area 17.547, 800 steps) has moments
0.1431 / 0.5539 / 0.5575 (`aniso_top` 0.6 %, ρ Fourier modes ≤ 2.4e-3): a torus of revolution to 0.6 %, exact
energy 30.276. The WP15 genus-1 run (slab start) is one too (0.140 / 0.556 / 0.564, 1.6 %), exact energy 31.58
(not yet at the plateau). The reference binary on the same OBJ, 800 steps: exact energy 30.144, `aniso_top`
0.78 %; ours repeated with the convergence rule (B2, 500 steps): 30.274, 1.1 %.

Radius ratio scan (circular profile, area 17.547, ~4.6k faces, exact energy): R/r = 1.3 / 1.5 / 1.6 / 1.65 /
2.0 / 2.5 / 3.0 → 56.6 / 31.58 / **30.30** / 30.54 / 38.5 / 58.1 (start) / 84.0. Flowed tori end below the best
circular one (non-circular profile): 29.67–30.27 at 4.1–4.6k faces, so a torus of revolution is where the
energy goes; there is no lower non-symmetric state.

Ablation, torus 48×24 (2304 faces, C₄₈-symmetric tessellation), 400–600 steps; `aniso_top` at the end:

| run | start rule | aniso_top | exact E |
|---|---|---|---|
| exact energy, dense metric, no remesh, no CCD (A1h) | τ₀ = 0.5 fixed | **6.6e-15** (exactly symmetric) | 30.440 (400) |
| + hierarchical metric (A4h) | τ₀ = 0.5 | 1.8e-6 | 30.440 |
| + Barnes–Hut θ = 0.25 (A3h) | τ₀ = 0.5 | 2.1e-4 | 30.441 |
| + Barnes–Hut θ = 0.5 (A2h) | τ₀ = 0.5 | 1.6e-3 | 30.484 |
| + remeshing (A5h) | τ₀ = 0.5 | 4.3e-3 | 31.457 |
| + CCD (A6 = A1 bit for bit: no step capped) | default | 1.3e-3 | 30.253 (600) |
| defaults (A7h / A7) | τ₀ = 0.5 / default | 1.9e-2 / 4.6e-3 | 31.87 / 31.99 |
| exact, dense, no remesh, no CCD (A1) | **default** (GrowPrevious ×2) | 1.3e-3 | 30.253 (600) |

* **Step overshoot amplifies round-off** (A1 vs A1h, torus 32×16 C1–C4). With the default start (τ doubles,
  then one halving: τ ≈ 0.8–1.3) the k = 2 Fourier mode of the exactly symmetric flow grows from 1e-15 by
  ≈ 1.66× per step (steps 20–64) and saturates at ~1e-3; with τ₀ = 0.5 fixed it stays at 1e-16 for 200
  steps, with τ₀ = 1.0 it grows. A stiff mode with τλ > 2 oscillates and grows while it is too small for
  Armijo to see. It is not an instability of the energy: the C₃₂-symmetrized state is **lower** (exact
  energy 32.32825 vs 32.32856, `symmetrize`), and the energy along X_sym + t(X − X_sym) is an even parabola
  (t = ±0.5 equal to 1e-9) — the torus of revolution is a local minimum in that direction.
* Barnes–Hut θ = 0.5 adds 1.6e-3 (θ = 0.25: 2e-4), the hierarchical metric 2e-6.
* Remeshing on an irregular mesh: ~0.5–1 %. Without remeshing, flowing a remeshed (irregular) torus with the
  exact energy and τ₀ = 0.5 (P3) *raises* `aniso_top` from 0.46 % to 1.24 % while the energy falls: the
  discrete minimizer on an irregular mesh is ~1 % anisotropic. A discretization effect, as in the reference.
* **Remeshing every step holds the energy up** (the stall of WP15). Energy changes within steps vs between
  steps (2nd half of the run): exact + dense + remesh (A5) −9.59 / +9.54; defaults on 2304 faces (A7) −11.64
  / +11.72; WP15 control torus −6.05 / +5.88, genus 1–3 and handcuffs the same pattern. The "between" part
  is almost all in steps without split/collapse/flip, i.e. the 5 tangential smoothing rounds. Without
  remeshing it is exactly 0 (exact energy) or +0.09 (θ = 0.5). Polishing the remeshed states without
  remeshing: A7 31.99 → 30.00 (P2, exact, dense, 600 steps), B2 30.27 → 29.67 (P1x, exact, hierarchical);
  the structured 48×24 grid without remeshing reaches 30.25 and 64×32 29.69 (800 steps, still falling slowly). So the remeshed plateau sits
  2–6 % above the minimum the same mesh can reach. The paper prescribes remeshing (N = 5, ρ = 0.5) after
  every step (Sec. 6) and the reference does the same (its torus plateau 30.14 has the same anisotropy).

## 2. Genus g (Fig. 15): 800 steps are too few; at convergence ours and the reference agree

Corrected moments of the WP15 finals (800 steps; exact energy):

| g | moments | E | | g | moments | E |
|---|---|---|---|---|---|---|
| 0 | 0.368 0.370 0.371 (triple) | 2.234 | | 6 | 1.476 1.504 1.637 (double, 1.9 %) | 55.87 |
| 1 | 0.140 0.556 0.564 (double, 1.6 %) | 31.58 | | 7 | 1.394 1.472 2.540 | 65.60 |
| 2 | 0.517 0.562 0.951 | 38.59 | | 8 | 1.430 1.559 3.107 | 71.91 |
| 3 | 0.501 0.946 1.581 | 52.21 | | 9 | 1.531 1.579 3.667 | 72.44 |
| 4 | 0.728 1.067 1.947 | 55.92 | | 10 | 1.407 1.639 5.118 | 93.72 |
| 5 | 1.092 1.192 1.840 | 54.31 | | 11 | 1.586 1.711 5.372 | 90.76 |

g = 0, 1, 6 are symmetric to 2 %; g = 7..11 are still elongated (one moment 2–3× the others).

Runs to convergence (stop: mean F over 200 steps changes < 1e-5 relative), same start OBJ, and the
reference for 4000 steps (threads 1):

| | ours (defaults) | ours θ = 0.25 | reference |
|---|---|---|---|
| g = 2 | converged at 1010: E 38.589, moments 0.517 0.563 0.949 | 1150: 38.604, 0.517 0.562 0.952 | step 1000: E 38.405, 0.522 0.563 0.943; BH F 38.50 at 1810 |
| g = 3 | converged at 2080: E 39.906, 0.874 0.923 0.943 | 1490: 40.224, 0.887 0.916 0.933 | step 1000: E 39.539, 0.893 0.918 0.931 |

* **g = 3 leaves the linear arrangement late**: F (BH) 58.3 (250) → 57.2 (500) → 52.5 (750) → 39.9 (1000)
  ours; the reference 58.1 (250) → 40.2 (500) → 39.4 (750). The linear arrangement is a slow saddle in both
  programs; the escape time depends on the noise (θ = 0.25 escaped at ≈ 700, the defaults at ≈ 900, the
  reference at ≈ 450). The WP15 run stopped at 800, during the escape (E 52.2). After it, ours and the
  reference reach the same near-tetrahedral shape (moments within 4–7 %, picture: four handles on the
  faces of a tetrahedron); our exact energy is 0.9 % above the reference's.
* **g = 2**: ours and the reference reach the same shape (two holes at an angle, a bent saddle; moments
  equal to 1 %, energies 0.5 % apart). It is not the paper's triangular-prism shape (a 3-fold axis would make
  two moments equal; here 0.517 vs 0.563). It is a local minimizer of the discrete flow: area-preserving
  random perturbations of 0.1 and 0.2 ℓ̄₀ return to it (F 38.64 / 38.69 after 750 steps, same moments).
  Neither program finds the prism from the linear start; the paper's start mesh and step count are not
  given.
* The remaining asymmetry (moment spreads of 2–7 % where the paper shows exact symmetry) and the 1 %
  energy gap to the reference are of the size of the remeshing/Barnes–Hut noise floor of §1; they are not
  a stall in the step sizes (τ 0.5–2, one halving per step, no CCD caps in any genus run) and not CCD.
* g ≥ 7: not rerun (each needs thousands of steps at 7–11k faces); given g = 3, the WP15 elongated shapes
  are very likely unfinished escapes from the linear arrangement.

## 3. Handcuffs (Fig. 1/14): the plane attractor, not CCD, keeps them linked

Same OBJs as `apps handcuffs` (linked: `shapes::handcuffs(1, 9000)`, 9084 faces; unlinked: two rings and a
bar, same area 152.44). Area + barycenter; `--plane W` = `ImplicitAttractor` plane z = 0, power 2,
`DualArea` weights (the `apps` setting is W = 5).

| run | F (BH) by step | outcome |
|---|---|---|
| reference, linked, **no plane** (RHf) | 186.8 → 27.1 (100) → 8.28 (300) → 7.94 (1100) | unlinks: one compact genus-2 blob (as EDT) |
| ours, linked, no plane, defaults incl. CCD (HF1) | → 26.9 (100) → 7.97 (300) → 7.96 (400–800) | **same as the reference**, curve for curve; 3 CCD-capped steps |
| ours, linked, no plane, CCD off (HF3) | → 24.5 (100) → 7.96 (400) | same end |
| ours, linked, no plane, CCD off, no remesh, swept check (HF2) | → 34.2 (200; 20 s per step, 296 GMRES iterations; stopped) | no swept contact in 200 steps; not yet unlinked |
| ours, linked, plane 5, CCD **off** (H1) | 375.6 (700–950) | stays linked, flattened (as with CCD on in WP15: 375.4) |
| ours, linked, plane 5, CCD off, no remesh, swept check (H2) | 374.3 (800) | stays linked; no swept contact |
| ours, unlinked blob (HF1 step 250) + plane 5 (HS) | 8.0 → 268 (300) → 264.4 (550) | flat two-hole plate, but not WP15's unlinked-start value 242.5: several local minimizers at weight 5 |
| ours, linked, plane 0.5 (HP05) | 126 (100) → 122.7 (200) → 65.2 (300) → 53.28 (700) | **unlinks** (flat plate, two holes) |
| ours, unlinked, plane 0.5 (HU05) | 54.2 (50) → 53.79 (250) | same plate; HP05's moments moving towards it (4.02/9.39 vs 3.22/11.25) |
| ours, linked, plane 0.05 (HP005) | → 15.04 (200–400) → 15.07 (850) | unlinks |
| ours, unlinked, plane 0.05 (HU005) | 15.08 (50) → 15.02 (350) | **same minimizer** as HP005: moments 0.741/3.558/7.353 vs 0.741/3.577/7.326 |
| reference, linked / unlinked, `implicit plane attract 2 5` (RHp / RHu) | 3281 / 2930 (600–700), TP 1138 / 983 | different values: the reference stays linked too |

* The reference's plane term is not area-weighted (start value 15 697 against our ≈ 500 for the same weight 5,
  ratio ≈ 30 ≈ 1 / mean dual area): its weight 5 is about our weight 150. With that plane it also keeps
  the linked handcuffs linked, flattened onto the plane (picture), at a higher energy than its unlinked
  start.
* CCD does not block anything here: no CCD cap in H1/WP15 linked runs, 3 caps in HF1, and with CCD off
  (H1) the result is the same 375.6. With the plane, the flattened linked state is a local minimizer of
  TP + plane: unlinking needs the rings to leave the plane (an IH-like move), which the attractor opposes.
* So "both flow to the same minimizer" holds for our flow without the plane and with a weak plane
  (0.05: same F and moments to 0.5 %; 0.5: F within 1 %, still converging): the paper's statement is
  reproduced with a weak attractor. At weight 5 (the `apps` value) even unlinked starts end in different
  local minimizers (242.5 vs 264.4), so switching a strong plane on after unlinking is not enough either.
* Unlinking with CCD on (HF1) is legitimate by construction (every line-search motion is collision-checked
  with the 0.05 ℓ̄₀ gap) and follows the reference's energy curve; the CCD-off runs end at the same F.

## 4. Eversion (Fig. 20): the setup crushes the torus flat

`apps eversion`: punctured torus 64×32 (R 2, r 0.8, hole radius 0.55 on the outer equator), hard area +
barycenter, `TotalVolume { weight 1e4, power 2 }`, boundary length (weight 30, target = initial length 5.38)
and curvature (weight 1), no remeshing, CCD on, 800 steps.

* The volume term dominates: F = 1e4 · 24.16² = 5.8e6 at the start; in 6 steps the tangent-point term goes
  22 → 1668 while V → 0.01: the torus is **pressed flat** into a two-sheet elliptic ring (moments 0.13 /
  2.44 / 10.99; the hole a flare at one end) and then relaxes slowly (TP 610 at step 800, τ ≈ 6e-5, still
  falling). Zero volume is reached by flattening, not by the paper's midsurface.
* Symmetry search with principal axes, all 48 signed permutations and ICP (not only coordinate
  reflections × quarter turns about z, which is all `apps eversion` tries): the WP15 midsurface has a C₂
  rotation (chamfer 7.7e-5 of the RMS radius) and a mirror through its mid-plane (chamfer 8.4e-3) that maps
  the top sheet to the bottom one — **with normals reversed** (agreement −0.997). That symmetry maps the
  oriented surface to its opposite, so reflecting and reversing the flow cannot evert it; no
  orientation-reversing isometry with agreeing normals exists (best +1 candidate: none below chamfer 4e-2).
  Same for the 1500-step rerun (E0) and for a hard volume constraint scheduled V₀ → 0 over 400 steps with
  remeshing (E3: V = 1.7e-4 at step 650, a dumbbell with the handle shrunk to a small ring; best
  orientation-reversing candidate chamfer 4.2e-2, normals −0.97).
* Without remeshing, two further variants **stall in the CCD**: volume weight 0.1 (E2) and a larger hole
  with boundary target 3× (E6) end `Stuck` after 46 / 35 steps. The limiting pair is a boundary vertex and a
  boundary triangle of its 2-ring (ring distances 1, 1, 2), approaching the gap (0.05 ℓ̄₀) with τ_max halving
  every step (4.3e-1 → 6e-11 → 0): a fold forming at the boundary, which the boundary penalties drive and
  the CCD freezes. Without remeshing nothing removes it. The hard volume schedule without remeshing (E1) is
  also slow: 225 GMRES iterations per step at step 100 (10 s per step), 74 of 150 steps CCD-capped.

No variant tried here produced the paper's midsurface: penalty 1e4 (flat ring), hard volume schedule with
remeshing (dumbbell), penalty 0.1 without remeshing (CCD stall). Two more (penalty 1, boundary target 6×
the initial length; a larger hole 1.3 with target 3×, both with remeshing) were stopped after 20 steps: the
stiff boundary-length term made GMRES fail every step and fall back to the dense LU (40 s per step,
τ ≈ 1e-7–1e-10). Suggested corrections, untested beyond the above: (1) a volume weight comparable to the
tangent-point energy (w V₀² ≈ ℰ₀, w ≈ 0.03, raised gradually over hundreds of steps) instead of 1e4;
(2) remeshing on, as the paper does every step (without it boundary folds stall the CCD); (3) a
boundary-length target that lets the hole open, approached gradually (a schedule on the target, not a
stiff jump); (4) the symmetry test over principal axes, all 48 signed axis permutations and ICP (a scratch
script did this here; `apps eversion` only tries coordinate reflections × quarter turns about z), keeping
only candidates whose image has agreeing normals.

## Other observations (for the solver / CCD owners)

* Open surfaces without remeshing (punctured torus): GMRES iterations grow from 34–58 at step 1 to 225–678
  per step by step 100–180 (10–27 s per step); with a strong boundary-length penalty (target 3–6× the
  initial length) GMRES fails and every step falls back to the dense LU (40 s, τ ≈ 1e-7–1e-10). The
  linked handcuffs without remeshing (HF2) go from 34 to 296 iterations in 200 steps.
* CCD stall on a boundary fold (E2, E6 above): a boundary vertex approaching a triangle of its own 2-ring
  is treated as a collision and τ_max halves every step until `Stuck`. A fold is a mesh-quality event, not
  a collision of distinct sheets; remeshing is what removes it.

## Proposals

SPEC §14.4 / NOTES-WP15 (text):

> Fig. 15: run each genus to convergence, not a fixed 800 steps: the linear arrangement is a slow saddle
> (g = 3 leaves it after 450–900 steps, in the reference too); use `--converge`-style windowed means
> (mean F over 200 steps changing < 1e-5) with a cap of 5000 steps. Judge symmetry with a working
> eigen-solver and expect 1–7 % moment spreads from the remeshing / Barnes–Hut noise floor.
> Fig. 14: the attracting plane (`apps` weight 5, `DualArea`) traps the linked handcuffs; use weight
> ≤ 0.5 (linked and unlinked then reach the same minimizer: F 53.3 / 53.8 at 0.5, 15.07 / 15.02 at 0.05,
> 1000 steps). The reference's `implicit plane attract 2 5` is about our weight 150 (its plane term is not
> area-weighted).
> Fig. 20: see §4 of NOTES-FIDELITY (volume weight ≈ ℰ₀/V₀², remeshing on, a symmetry test over
> principal axes + ICP with agreeing normals).

SPEC §10 / §7 (text, new paragraph):

> Remeshing every step (paper Sec. 6, N = 5, ρ = 0.5; the reference does the same) sets a noise floor:
> near a minimizer the tangential smoothing raises the energy by as much as the step lowers it (torus,
> exact energy, dense metric: −9.6 / +9.5 per 300 steps), holding F 2–6 % above what the same mesh reaches
> without remeshing. For converged shapes, end with a polish phase without remeshing (or remesh only when
> an edge leaves [L₀/2, 3L₀/2]), and evaluate final energies exactly. Not changed in the library.

SPEC §7.2 (text):

> `GrowPrevious ×2` accepts τ above the stability limit of stiff modes (τλ > 2); such a mode grows from
> round-off by ≈ 1.7× per step until Armijo sees it (exact symmetric torus: k = 2 mode 1e-15 → 1e-3), which
> breaks exact symmetries at the 1e-3 level. A fixed τ₀ = 0.5 keeps them to round-off but converges more
> slowly. Harmless for minimizers (the symmetric state is lower); documented, not changed.
