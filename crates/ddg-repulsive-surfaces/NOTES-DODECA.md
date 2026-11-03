# NOTES-DODECA: genus 11 — row start vs the dodecahedron (Fig. 15)

Question: our g = 11 flow (`NOTES-GENUS.md`) ends in an irregular, elongated shape (moments 1.73 / 2.03 / 4.12,
E 81.9). The paper's Fig. 15 shows a g = 11 minimizer with the symmetries of the dodecahedron, "initialized by
a linear arrangement of handles". How does the authors' program get there, and what is different in ours?

**Answer: neither program reaches the dodecahedron from a linear start, and the two do the same thing.**
From the same row start, the authors' program goes to the same kind of irregular elongated shape as ours, with
about the same energy and moments. A dodecahedral start relaxes in both programs to a near-dodecahedral shape
whose energy is 8 % below the row results (4.5 % after polishing, see below). So Fig. 15's g = 11 figure is very
likely not the end of a flow from a linear row like ours. The paper gives no start mesh, face count or step count
for it: it may have been flowed from a symmetric start or a different arrangement, or picked from several
starts. There is no setting of ours that decides the outcome (experiment 4 was not needed: the reference does
not reach it either), and no library change was made.

## Method

Tool: `examples/dodeca.rs` (`gen`, `ladder`, `perturb`, `flow`, `measure`; see its header). All runs: area +
barycenter constraints, (α, β) = (6, 12), no attractor. Ours: `FlowOptions::default()` (Barnes–Hut θ 0.5,
hierarchical solver, CCD, remeshing) unless stated, stopping when the mean of F over the last 200 steps is within
1e-5 (relative) of the mean over the 200 before (the `ishii genus` rule), cap 8000; 4 or 8 threads. Reference:
`RepulsiveSurfaces-build/build/bin/rsurfaces_headless --threads 1`, scene `repel_mesh mesh.obj 6 12 /
constrain barycenter / constrain area`, its default remeshing, OBJ every 100 steps; it has no stop rule, so it
was run for 2300–2500 steps (the same windowed rule, applied to its `steps.csv`, is met at step 1298 for the row
and 414 for the dodecahedron) and stopped.

Starts, all at the same area 173.4343:
* **row**: `shapes::genus_surface(11, 9600)` (10,136 faces; the `ishii genus 11` start, `frame_00000.obj` of
  `runs_tmp/genus2/g11/`), a rounded slab with 11 holes in a row. Exact E 7956 (sharp slab edges).
* **ladder**: a tube (radius 0.22 of the spacing) around the ladder graph with 12 rungs (24 vertices, 34 edges,
  genus 11), 10,300 faces: a chain of 11 rings, closer to the row drawn in Fig. 15's top row.
* **dodecahedral**: `generate::graph_tube` around the 30 edges of a regular dodecahedron (20 vertices,
  30 − 20 + 1 = 11 handles), tube radius 0.30 of the edge length (10,116 faces, E 302.5); a second one with radius
  0.38 (10,220 faces). Both embedded, genus 11, inertia isotropic to 1e-3.
* **perturbed**: each vertex moved by ε mean edge lengths in a random direction (seeded xorshift), rescaled to
  the area. Note: on the thin row slab these starts **self-intersect** (local folds; ε = 0.2 and 0.5). Ours
  untangles them (CCD band rule, every final is embedded); the reference blew up on the ε = 0.5 ones (F 1e6–1e9
  after 200–300 steps, stopped) and flows the ε = 0.2 ones.

Indicators: exact all-pairs energy (`EnergyEval::Exact`), normalized to the start area (E ∝ area⁻¹ at (6, 12);
all final areas are within 2e-5, so the correction is below 1e-4); eigenvalues of the area-weighted second-moment
tensor about the area centroid (fixed `eig3`); `Q_l = Σ_{f,g} a_f a_g P_l(u_f·u_g) / A²`, the rotation-invariant
degree-l power of the area distribution over the directions u from the centroid. Icosahedral symmetry has no
invariant harmonics of degree 2 and 4, so a dodecahedral shape has Q2 = Q4 = 0 and Q6 > 0 (the dodecahedral
start: Q2 9e-8, Q4 2e-7, Q6 1.1e-2); a cube has Q4 > 0, an elongated shape Q2 > 0.

## Results (remeshed final states, exact energy at area 173.4343)

| start | program | steps | E | moments | Q2 | Q4 | Q6 |
|---|---|---|---|---|---|---|---|
| row | ours (`NOTES-GENUS`) | 1711 | 82.63 | 1.733 2.034 4.123 | 3.2e-2 | 1.2e-3 | 1.4e-3 |
| row | ours `reference_compat` | 2563 | 81.62 | 1.852 2.016 3.947 | 2.3e-2 | 8.9e-4 | 8.2e-4 |
| row | **reference** | 2500 | **80.37** | 1.787 1.913 4.261 | 3.4e-2 | 1.1e-3 | 8.5e-4 |
| ladder | ours | 1472 | 81.34 | 1.735 2.040 4.137 | 3.0e-2 | 1.1e-3 | 1.0e-3 |
| ladder | reference | 2000 (rule 1583) | 81.49 | 1.764 2.136 3.922 | 2.4e-2 | 1.6e-3 | 1.1e-3 |
| dodecahedral | ours | 435 (rule) | 75.94 | 2.408 2.427 2.639 | 2.3e-4 | 3.5e-5 | 3.1e-3 |
| dodecahedral | ours, 3000 steps | 3000 | 76.24 | 2.414 2.427 2.632 | 2.1e-4 | 3.2e-5 | 3.1e-3 |
| dodecahedral | ours θ = 0.25 | 2000 | 76.71 | 2.416 2.424 2.638 | 2.2e-4 | 3.6e-5 | 3.2e-3 |
| dodecahedral | ours `reference_compat` | 696 | 75.83 | 2.423 2.482 2.567 | 7.5e-5 | 1.5e-5 | 3.1e-3 |
| dodecahedral r 0.38 | ours | 2000 | 77.10 | 2.419 2.487 2.530 | 5.3e-5 | 7.1e-6 | 1.8e-3 |
| dodecahedral | **reference** | 2300 | **75.15** | 2.417 2.457 2.597 | 1.3e-4 | 1.1e-5 | 3.0e-3 |

Trajectories from the row start (exact E; moments): ours 126.6 (step 500; 0.91 1.55 10.49, still a row), 83.3
(1000), 82.8 (1500), 82.6 (1711); the reference 119.3 (500; 1.17 1.36 9.67), 80.4 (1000; 1.72 1.89 4.45), 80.6
(1500), 80.5 (2000), 80.4 (2500; 1.79 1.91 4.26). Both leave the row between steps 500 and 1000 and fold it into
the same kind of lumpy shape (Q2 ≈ 0.03, one moment 2.1–2.4× the others), not a polyhedral cage. The
dodecahedral start: ours 76.1 at step 100 and flat after; the reference 75.04 (100), 75.08 (500), 75.09 (1000),
75.15 (2000).

Perturbed row and ladder starts, ours (none reaches a symmetric shape; Q2 1.4e-2 – 5.4e-2):

| start | E | moments | | start | E | moments |
|---|---|---|---|---|---|---|
| row ε 0.2 seed 1 | 87.12 | 1.836 1.909 4.178 | | row ε 0.5 seed 4 | 85.01 | 1.772 2.309 3.629 |
| row ε 0.2 seed 2 | 83.64 | 1.577 2.595 3.466 | | row ε 0.5 seed 5 | 86.98 | 1.405 2.886 3.463 |
| row ε 0.2 seed 3 | 83.88 | 1.682 1.792 4.688 | | row ε 0.5 seed 6 | 85.29 | 1.563 2.700 3.355 |
| row ε 0.2 seed 7 | 82.44 | 1.691 1.693 4.914 | | ladder ε 0.2 seed 1 | 85.42 | 1.766 1.869 4.443 |
| ladder ε 0.5 seed 2 | 84.52 | 1.750 1.949 4.263 | | | | |

Reference from perturbed rows (ε = 0.2, 2000 steps; the rule is met at 1610 / 1272): seed 1 84.51 (1.858 1.959
4.128, Q2 2.8e-2), seed 2 82.38 (1.941 2.142 3.659, Q2 1.4e-2). Same shape class; ours from the same starts 87.12 /
83.64 (our remeshed values run 1–3 % higher, see the polishing below).

Perturbed dodecahedral starts, ours, 2000 steps: ε 0.2 seed 1 77.05 (2.445 2.497 2.521, Q2 2.6e-5), seed 2 78.08
(2.429 2.474 2.569, Q2 6.5e-5), ε 0.5 seed 3 79.69 (2.431 2.500 2.536, Q2 4.7e-5): all return to the
near-dodecahedral shape (Q4 ≤ 2.4e-5, Q6 3.0–3.3e-3); the dodecahedral state is a stable (local) minimizer.

### Polishing: equal-footing energies and the residual distortion

Remeshing every step sets an energy noise floor (`NOTES-FIDELITY.md` §1). The final states above were flowed
300 more steps with ours, remeshing off (`--no-remesh`, everything else default):

| final state | E before → after 300 steps | moments after | Q2 after |
|---|---|---|---|
| row, ours | 82.63 → 71.48 | 1.790 2.048 3.998 | 2.6e-2 |
| row, ours `reference_compat` | 81.62 → 71.68 | 1.837 2.073 3.873 | 2.1e-2 |
| row, reference (step 2300) | 80.37 → 71.59 | 1.802 2.026 4.021 | 2.6e-2 |
| ladder, ours | 81.34 → 71.37 | 1.796 2.034 4.012 | 2.6e-2 |
| dodecahedral, ours (3000 steps) | 76.24 → 68.24 | 2.461 2.475 2.510 | 9.3e-6 |
| dodecahedral, ours `reference_compat` | 75.83 → 68.33 | 2.460 2.482 2.505 | 8.0e-6 |
| dodecahedral, reference (step 2100) | 75.15 → 68.27 | 2.473 2.483 2.492 | 1.2e-6 |

* After polishing, ours and the reference agree to 0.15 % on both shapes (row 71.48 / 71.59, dodecahedron 68.24 /
  68.27): the 2–3 % gap of the remeshed states (82.6 vs 80.4) is the remeshing noise floor, not a different
  minimizer. The row shape stays at 71.4–71.7 in all four variants; the dodecahedron is 4.5 % lower.
* The polished energies are still falling slowly (68.23 → 68.01 F over the last 50 steps) and part of the drop
  is tangential mesh drift: edge lengths spread from 0.71–1.28 to 0.54–1.74 of the mean (1st/99th percentile),
  minimum angle 25° → 10°. So 68.2 / 71.5 are values on the frozen connectivity, not the continuum energy; they
  compare shapes at equal footing, the ranking is the same as for the remeshed states.
* **The remeshed dodecahedral states are stretched by 7–9 % in one moment** (ours along a 3-fold axis at step
  435, drifting towards a 2-fold axis by step 3000; the reference 6.6 % at step 800). This is not a symmetry
  breaking of the energy: it appears in both programs, its axis wanders, other starts (r 0.38, perturbed,
  `reference_compat`) show 2–4 % in other directions, and polishing without remeshing removes most of it
  (spread 2.0 % ours, 0.8 % the reference; Q2 2e-4 → 1e-5 / 1e-6). It is the same remeshing noise floor.

## Verdict

1. **The reference does not reach the dodecahedron from our row start.** It ends at E 80.4 (remeshed; 71.6
   polished) in a lumpy, elongated shape like ours (82.6; 71.5 polished). Both leave the row between steps 500 and
   1000 and, at equal footing (polished), agree in energy to 0.15 %.
2. **A dodecahedral start is a stable near-dodecahedral minimizer in both programs**, 8 % (remeshed) / 4.5 %
   (polished) below the row results: ours 75.9–76.2, the reference 75.15; polished 68.24 / 68.27. Perturbed
   dodecahedral starts return to it.
3. No start in the "linear" family reached it: 8 perturbed rows, the ladder (a chain of rings) and 2 perturbed
   ladders with ours; the row, the ladder and 2 perturbed rows with the reference. The linear start falls into a different, higher local
   minimum (the many hole arrangements of genus 11 are separate basins). Fig. 15's g = 11 image was most likely
   flowed from a different start than a plain row (for g = 7 the paper says it also tried a "symmetric
   configuration" as a start), or the row was perturbed in a way not described.
4. No setting decides it: `reference_compat` (the reference's barycenter, Armijo rule and start, preconditioner
   form and left GMRES, Newton step, remeshing order) gives the same shape class from both starts, and so do
   θ = 0.25 and long runs. The ablations of experiment 4 were not needed. No defect was found; `src/` is unchanged.

## Files (`runs_tmp/dodeca/`, not committed)

Starts `row_start.obj`, `ladder_start.obj`, `dodeca_start.obj`, `dodeca_start_r038.obj`, `row_p0*_s*.obj`; ours
`ours_*/{steps.csv,frame_*.obj,final.obj}` + `ours_*.log`; reference `ref_*/{scene.txt,mesh.obj,out/}`; polishing
`polish_*`; measurements `measure_*.log`; the comparison sheet `dodeca_comparison.png` (row and dodecahedral
results of both programs, and the polished reference shapes; `render/` holds the panels); the paper's page
`fig15-13.png`. Uptime: `uptime_start.txt` (load 0.26 at the start, 3.1 at the end, 15–55 during the runs, all from
these jobs).
