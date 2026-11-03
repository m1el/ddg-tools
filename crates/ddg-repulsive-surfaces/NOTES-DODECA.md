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

**Update ("The video's start", at the end).** The paper video shows the start: a chain of 11 fat, fused tori with
small holes, not our row slab or ladder. From a replica of it, both programs follow the video's path and reach the
dodecahedron, at the same energy as from the dodecahedral start. This holds for hole 0.25 of the spacing (the best
match to the video) but not for 0.15, 0.2 or 0.3; 3 of 18 chain runs reach it. The linear start does reach the
dodecahedron, but only from the right chain.

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

## The video's start

The paper video (2:22–2:30) shows a genus-11 flow from a chain of rings ending in the dodecahedral cage ("the method
does not easily get stuck in local minima"). Its start is not our row slab or ladder: it is a chain of 11 fat round
tori with small holes, fused by necks. This section replicates that start and flows it with both programs.

**Answer: from the replicated start, both programs follow the video's path (chain → flat bar with small holes →
bend and fold → blob → cage) and reach the dodecahedron, but only for some of the starts.** The hole size 0.25
(the best match to the video) reaches it in both programs. Hole sizes 0.15, 0.2 and 0.3 do not. A perturbed
hole-0.3 start reaches it in the reference only. In total, 3 of 18 chain runs reached it: ours 1 of 11, the
reference 2 of 7. The 16 row and ladder runs above reached it 0 times. So the video's start is far better than a
row, but the outcome hangs on details of the start: the dodecahedron is one basin among several.
Nothing was found to fix; `src/` is unchanged.

### Measuring the start

Frame 2:22.5 (1920×1080, `runs_tmp/dodeca2/video/full_142.5.png`; the genus 3 and genus 5 starts are in the same
frame and come from the same family), column profile of the green pixels:
* chain 928 px long, ring spacing s = 86.2 px. Projected height at the rings 48 px = 0.557 s, at the necks
  20 px = 0.232 s, length / height 19.3.
* The shadow is close to a top view. Its holes are 24–27 px wide (0.29–0.33 s, a lower bound because the
  penumbra shrinks them). It is compressed 0.54× vertically: the round holes are 13 px tall. Corrected for
  that, the chain is 0.83 s wide and the neck 0.31 s wide.
* Viewing elevation: the genus-3 hole opening (107 × 49 px) gives about 27°. The height fit gives 26°.

The start (`dodeca chain`): 11 round tori (axis z) at unit spacing, outer diameter 0.86, hole (throat) diameter H,
so major radius (0.86 + H)/4 and tube radius (0.86 − H)/4. Neighbours are fused by a polynomial smooth union whose
blend width (0.368) makes the neck 0.30 wide. The surface is meshed by marching tetrahedra (about 12k faces;
the flow's own remeshing runs from step 1) and scaled to the row's area 173.4343. All starts are genus 11,
embedded, one component. Rendered from the video's view (azimuth 0, elevation 26°; `render/start_h*.png`,
`start_compare.png`):

| H | height (rings) | height (necks) | length / height |
|---|---|---|---|
| video | 0.557 s | 0.232 s | 19.3 |
| 0.15 | 0.579 | 0.243 | 18.8 |
| 0.20 | 0.573 | 0.243 | 19.1 |
| **0.25** | **0.550** | **0.236** | **19.8** |
| 0.30 | 0.545 | 0.229 | 20.1 |

H = 0.25 matches best. The silhouette alone fixes H only to about ±0.05, which is why 0.15–0.30 were all flowed.
The video's holes look a little smaller than ours in the top rim, and its necks are a little fuller. Neighbouring
tori do not overlap (outer radius 0.43 < 0.5): the necks are fillets.

### Flows (area + barycenter, as in the method section above)

Ours: `dodeca flow`, defaults, `ishii genus` stop rule. Reference: `rsurfaces_headless --threads 1`,
3000 iterations, OBJ every 100. Exact energy at area 173.4343; moments; Q2 Q4 Q6. For comparison, the dodecahedral
start gave ours 75.94 / reference 75.15 (Q2 ≈ 1–2e-4).

| start | program | steps | E | moments | Q2 | Q4 | Q6 | shape |
|---|---|---|---|---|---|---|---|---|
| H 0.15 | ours | 2096 | 79.23 | 1.814 2.052 3.956 | 2.4e-2 | 1.0e-3 | 8.8e-4 | elongated |
| H 0.15 | reference | 3000 (rule 1962) | 78.82 | 1.731 2.098 3.975 | 2.6e-2 | 9.9e-4 | 1.4e-3 | elongated |
| H 0.20 | ours | 1877 | 88.16 | 1.548 1.722 5.533 | 7.7e-2 | 2.3e-3 | 8.0e-4 | elongated |
| H 0.20 | reference | 3000 | 81.76 | 1.870 2.071 3.880 | 2.1e-2 | 1.6e-3 | 9.4e-4 | elongated |
| **H 0.25** | **ours** | 2601 | **76.48** | 2.354 2.464 2.650 | 3.8e-4 | 2.2e-5 | 2.4e-3 | **dodecahedral** |
| **H 0.25** | **reference** | 3000 (rule 2516) | **75.81** | 2.334 2.509 2.615 | 3.2e-4 | 3.1e-5 | 2.3e-3 | **dodecahedral** |
| H 0.30 | ours | 1912 | 80.69 | 1.569 2.789 3.218 | 1.8e-2 | 1.1e-4 | 8.9e-4 | flattened cage |
| H 0.30 | reference | 3000 (rule 1948) | 79.11 | 1.805 1.857 4.350 | 3.6e-2 | 1.1e-3 | 8.8e-4 | elongated |

The fine mesh, the small perturbations and the no-area runs are in "Variants" below.

Polished (300 steps more with ours, no remeshing, the footing of the polishing table above): **H 0.25 ours
76.48 → 68.51, reference 75.81 → 68.21** (moments 2.455 2.480 2.515 / 2.452 2.497 2.499, Q2 1.5e-5 / 1.1e-5, Q4
6e-7 / 5e-7). This is the dodecahedral state of the polishing table above (68.24 ours, 68.27 the reference). The
H 0.30 finals polish to 71.56 (ours) and 71.05 (reference), as the row's 71.4–71.7 do. **The dodecahedron reached
from the video's start is the same minimizer as the one from the dodecahedral start**: same energy within 0.4 %,
same symmetry, both programs.

**Trajectory** (`render/traj_ours_h25.png`, `render/traj_ref_h25.png`, `render/traj_ours_h30.png`; panels at
steps 0, 100, 300, 500, 700, 1000, 1500, final). Both programs follow the video's sequence:
* the chain first shortens and fattens into a straight flat bar with 11 small holes (ours step 100, reference
  step 300);
* the bar bends and folds (300–500 ours, 500–700 reference);
* it crumples into a compact blob (700 ours, 700–1000 reference);
* the blob opens into a cage (1000–1500).

The H 0.30 run (rendered) takes the same first steps (H 0.15 and H 0.20 by their moments); they differ in how the bar folds. Two details of the
successful runs are not visible in the video, whose playback speed is not the step count:
* Ours from H 0.25 passes through the same elongated lumpy state as the failed runs: at step 1000, E 80.1, moments
  1.91 2.12 3.68, Q2 1.5e-2. It leaves that state between steps 1000 and 1500 (Q2 1.3e-3 at 1500, 3.8e-4 at 2000).
* The reference from H 0.25 reaches the cage between steps 1500 (79.0) and 2000 (75.6).

### Variants (experiment 3)

* **Small perturbations** (each vertex moved by ε mean edge lengths, as above). The reference from H 0.30, ε 0.05,
  seed 3 **reaches the dodecahedron**: 75.04 at step 3000, moments 2.439 2.489 2.531, Q2 3.5e-5, Q4 1.9e-5, Q6
  2.4e-3. It sat in an elongated state first: 79.06 at step 2000 (Q2 1.4e-2). Its F falls from 78.7 to 74.8
  between steps 2000 and 2800, and the stop rule is not yet met at 3000, so a longer run may go a little lower.
  Ours from the **identical** start (the same `perturb`) ends elongated at 80.81 (3455 steps, Q2 1.8e-2), as do
  ours from H 0.30 with ε 0.05 seeds 4 and 5 (79.92, 83.67), ε 0.2 seeds 1 and 2 (79.96, 80.33), and H 0.20
  ε 0.05 seed 3 (82.52). The reference from H 0.30 ε 0.2 seed 1: 79.29, elongated. On the same start, the two
  programs can fold the bar differently.
* **Finer mesh** (H 0.30, 30k faces): ours 1830 steps, E 71.88 at 32k faces; the reference was stopped at step
  2224 (the rule was met at 2193), E 71.45. Both are elongated, with moments about 1.8 1.9–2.0 4.1–4.2 and Q2
  2.8–3.0e-2. The energies are lower because the finer mesh lowers the remeshing floor, so they are not comparable
  to the 12k values. The fine mesh does not help.
* **No area constraint** (barycenter only, H 0.30). The (6, 12) energy scales as area⁻¹, so the flow mainly
  inflates the surface. The reference grows the area 6× by step 600 and 17× by step 1500 (122k faces, since
  remeshing keeps the edge length). It stays a straight chain: moment ratio 81 at step 1500, no fold. It was
  stopped there. Ours grew to 268k faces by step 75 and was stopped. This is not the video's path, which folds.
  The growth of the blob on screen after 2:26 is then camera zoom, not a missing area constraint.
* **Are the elongated finals minima?** Four of our finals got 1700–1800 more steps without the stop rule:
  H 0.30 80.20, H 0.15 78.29, H 0.20 87.75, H 0.30 ε 0.2 seed 1 79.48. The reference's H 0.30 and H 0.15
  finals got 800 more iterations: 78.84 and 78.53. Moments and Q2 did not change in any of them. Two runs above
  did leave an elongated plateau after 1000 steps (ours H 0.25 at steps 1000–1500, the reference ε 0.05 at
  2000–2800). So some plateaus are slow saddles, but these six held for 2300–3700 steps in total (counting from
  where the plateau began) without moving.

### Verdict

1. **The video's claim is reproduced, with a qualification.** From the replicated start (11 fused fat tori, outer
   0.86, hole 0.25, neck 0.30 of the spacing), both programs follow the video's path chain → bar → fold → blob →
   cage. Both reach the dodecahedral minimizer at the same energy as the dodecahedral start: 76.5 / 75.8 remeshed,
   68.5 / 68.2 polished, against 75.9 / 75.15 and 68.24 / 68.27.
2. **The start decides it, not a setting.** The row slab and the ladder (thin tubes, large rectangular holes)
   never reached it in 16 runs. The video's chain reached it in 3 of 18. Hole 0.15, 0.2 and 0.3 and most small
   perturbations fold the bar into an elongated or flattened cage that is 4–6 % higher (polished 71.0–71.6).
   That state is stable for thousands of steps in both programs. "Does not easily get stuck" holds for this start
   family relative to the row starts. It does not mean that every chain reaches the dodecahedron: the authors'
   run very likely used a start like H ≈ 0.25 that happens to fold the right way.
3. Neither the finer mesh nor removing the area constraint is needed or helps. Without it, the chain only inflates.

### Files (`runs_tmp/dodeca2/`, not committed)

* Starts `chain_h{15,20,25,30}.obj`, `chain_h30_fine.obj`, `chain_h30_p0*_s*.obj`; generation logs `gen_*.log`.
* Ours `ours_*/` + `ours_*.log`; reference `ref_*/{scene.txt,mesh.obj,out/}`; continuations `cont_*`; polishing
  `polish_*`.
* Measurements `measure_*.log`; job lists and runners in `bin/`.
* Video frames and crops `video/`.
* Renders `render/`: the starts in the video's view `start_h*.png`, trajectories `traj_*.png`, finals `final_*.png`,
  principal-axis views `views1.png`.
* `start_compare.png`: the video start above the H 0.15 / 0.2 / 0.3 starts.
* **`dodeca2_comparison.png`**: the video start, our H 0.25 start, both trajectories, and the video's result at
  2:30 next to the four finals.
* The frame-measurement and sheet scripts were run from the session's scratch directory and are not kept.

Uptime: 3.5 at the start (00:14). 63 at the end of the main batch (04:39) and 56 at the last stop (05:49). From
about 04:30, another agent's `ishii` runs (about 24 jobs × 3 threads) shared the machine. The step counts are
not affected; the seconds in the logs are.
