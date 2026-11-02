# WP15 notes (applications of Sec. 9, SPEC §14.4)

Code: `examples/apps/{main,apps}.rs` (one example, one subcommand per application; `apps list`),
`tests/applications.rs` (smoke tests at reduced size, 8 tests + 1 `#[ignore]` heavy one, 5 s),
`showcase/repulsive-surfaces/{run_all.sh,obj_stats.py,README.md}`.
Build/run: `CARGO_TARGET_DIR=$PWD/target-wp15 cargo build --release -p ddg-repulsive-surfaces --example apps`,
`scripts/rs-fetch-scenes.sh`, then `showcase/repulsive-surfaces/run_all.sh [--smoke] [NAME...]` (16 jobs x 4 threads;
logs `runs_tmp/wp15/full/<name>.log`, `RESULT` lines; frames `showcase/media/repulsive-surfaces/<name>/`, untracked).
All runs: `StopRule::Off` (see "Stop rule" below), CCD on (off for the subcritical runs), remeshing on unless noted,
exponents (6, 12) unless noted. Numbers below are from the full runs of 2026-11-02, 4 threads each, 16 jobs at a time
on a machine with load 30 to 50 (timings are not clean).

## Results

| application | start faces -> end | steps | time | outcome |
|---|---|---|---|---|
| genus 0..11 (Fig. 15) | 3.0k+0.6k g -> 3.3k..11.1k | 800 each | 6 to 46 min | no step error, no self-intersection, area kept to 1e-5. Sphere (g = 0) is round: moment eigenvalues equal to 0.5 %. Other genera did **not** reach the paper's symmetric shapes in 800 steps (table below) |
| genus 2..6 + plane, weight 5 and 100 (Fig. 17) | same | 800 | 14 to 37 min | flatter: the smallest moment eigenvalue is 0.15 to 0.31 of the middle one at weight 5 (g = 2..6), 0.01 to 0.8 at weight 100 (see the logs); no self-intersection. A hard plane constraint is not expressible (no constraint of that kind), weight 100 is the stand-in |
| handcuffs, linked and unlinked (Fig. 1, 14) | 9.1k -> 10.2k | 1500 | 45 and 20 min | **fails to reproduce Fig. 14**: final energies 375.4 (linked start) against 242.5 (unlinked start); the linked pair stays linked at a local minimizer (energy flat to 4e-4 relative over the last 900 steps). Both genus 2, no self-intersection |
| HandcuffsPole stage 1 + stage 2 (Fig. 18) | 10.7k -> 14.1k -> 15.1k | 1500 + 600 | 47 + 12 min | handle comes off the pole: smallest distance of any vertex to the pole axis 0.105 (start) -> 0.36 (step 500) and the surface flattens onto the attractor plane (z range +-0.77 -> +-0.16); stage 2 (cylinder weight 0.1 -> 0.000125, new session) settles at distance 0.115 from the axis. Area kept, no self-intersection. "Handle off the pole" is judged by distances only (no link test) |
| SlideOffBar stage 1 + 2 | 5.2k -> 6.0k -> 6.8k | 1500 + 600 | 31 + 6 min | runs, energy 2.1e5 -> 7.4e3 -> 9.5e2, no self-intersection. `end.obj` is missing in the reference, stage 2 continues from stage 1 |
| punctured torus eversion (Fig. 20) | 4.0k (no remeshing) | 800 | 33 min | volume 24.16 -> 0.012 (5e-4 of the start), area and boundary length kept, no self-intersection; **but the midsurface has no orientation-preserving symmetry** among reflection x quarter-turn candidates (best chamfer 0.063 of radius 5.8, normals disagree: -0.98), so the reversed flow does not end inside-out (V_end = +24.2 with the winding that keeps the path continuous). The OBJ sequence is written but is not an eversion |
| H shape, 7 variants (Fig. 21/22) | 3840 (no remeshing) | 300 | 8 to 38 min | all run; pins stay put (drift 0), volume x2 (or area x1.5) reached to 1e-8, no self-intersection; whether the box variants stay inside the box was not measured. Hard pins with the **default hierarchical solver fail** at step 0 (below); the results use `Backend::Dense` |
| shrink wrap bunny points / spot points / bunny soup (Fig. 23) | sphere20k 20480 -> 15.3k / 13.4k / 13.3k | 220 | 15 / 15 / 13 min | volume 44.6 -> 4.90 (bunny, soup) and 4.01 (spot), no self-intersection. Soup: 18560 of 28576 bunny faces kept, as a mesh obstacle. Not compared to the data quantitatively (no ground-truth distance computed) |
| nested envelopes, 4 levels (Fig. 24) | spot 16.8k -> 8.7k -> 4.0k -> 1.9k -> 0.9k | 150 each | 10 + 2 + 0.4 + 0.1 min | volume x1.150 per level (target 1.15), no self-intersection; previous level's vertices inside the next: 1.0000, 1.0000, 0.9975, 0.957. Level 3 with the default solver fails after 2 steps (below); levels 2..4 use dense |
| walnut (Fig. 25) | 1280 -> 7142 | 220 | 17 min | area x8.00 exactly, no self-intersection |
| cobblestones (Fig. 25) | 14 spheres 4480 -> 7192 | 400 | 14 min | each area x2.5000, no self-intersection between spheres, **but 218 vertices (max z 1.42) are outside the slab** (half thickness 0.9) at the end; also without remeshing (152 vertices, step 320 `Stuck`) and with the slab weight 30 and 300 (see failures). Default solver needs 2800 GMRES iterations per step (30 to 60 s per step); dense takes 2 s |
| subcritical, spindle torus (Fig. 26 top) | 4096 | 600 / 2000 | 9 min / see below | does **not** resolve the intersection at 600 steps: intersecting triangle pairs 13096 -> 791 (3.5, 7) and -> 856 (3, 6). The start is "slightly too fat" (R = 1, r = 1.02), not a small patch |
| Klein bottle, (3.5, 7), no low-order term (Fig. 26 bottom) | 1800 | 700 | 2.3 min | keeps intersecting as expected: 70 -> 49 pairs; area kept (6e-7) |

Genus symmetry indicators (area-weighted second-moment tensor about the area centroid, eigenvalues ascending; multiplicity at 3 %):

| g | eigenvalues | ratios | mult. | | g | eigenvalues | ratios | mult. |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.3684 0.3701 0.3703 | 1.005 1.001 | 3 | | 6 | 1.476 1.507 1.633 | 1.022 1.083 | 2 |
| 1 | 0.2945 0.4060 0.5595 | 1.378 1.378 | 1 | | 7 | 1.471 1.950 1.985 | 1.326 1.018 | 2 |
| 2 | 0.5171 0.5759 0.9370 | 1.114 1.627 | 1 | | 8 | 1.430 1.790 2.876 | 1.252 1.607 | 1 |
| 3 | 0.5010 0.9473 1.5798 | 1.891 1.668 | 1 | | 9 | 1.533 1.884 3.360 | 1.228 1.784 | 1 |
| 4 | 1.066 1.210 1.465 | 1.134 1.211 | 1 | | 10 | 1.421 1.751 4.993 | 1.232 2.852 | 1 |
| 5 | 1.093 1.351 1.681 | 1.236 1.244 | 1 | | 11 | 1.614 1.684 5.372 | 1.043 3.191 | 1 |

Reading: nothing is at the paper's symmetry (torus of revolution g = 1, tetrahedron g = 3, cube g = 5 would have
multiplicity 2 or 3 for these moments; g = 4, 8, 9 prism-like shapes would have two equal ones). Only g = 6 and g = 7 show
a double eigenvalue. g = 9..11 are still elongated: they are not converged in 800 steps (the energy still falls,
non-monotonically, at the end). Control: a torus of revolution (R/r = 2.5, the area of genus01) flowed with the same
options (`apps control_torus`, 800 steps) ends at energy 30.69 with moments 0.278 0.422 0.554 (not round either) and lower
than the linear-arrangement run (31.9 to 32.0, plateau): so the g = 1 run is not at the lowest energy found, and the flow
itself breaks the rotational symmetry of a torus (remeshing and Barnes-Hut noise are the suspects; not investigated).

## Failures and reproducers

1. **Hierarchical GMRES diverges** (`RsError` "the linear solver diverged: relative residual 0.96 (above 1e-1)") at step 0
   whenever hard `Pins` are present on a mesh above the dense limit: `apps h_sparse_volume`, `h_sparse_area`,
   `h_sparse_volume_box` with `RS_AUTO_BACKEND=1` (3840 faces; the pin-free `h_nopins_volume` is fine). Also on
   self-intersecting inputs: `RS_AUTO_BACKEND=1 apps subcritical` (1024 faces smoke, step 0) and `RS_AUTO_BACKEND=1 apps klein`
   (step 74, residual 0.29); and `envelope_level3` after 2 steps (residual 0.1285, `RS_AUTO_BACKEND=1 apps envelopes`).
   Workaround used: `Backend::Dense` (3k to 7k faces). Same family as the open defect of SPEC §16 (compat-noedge divergence).
2. **Many components: GMRES needs about 2800 iterations per step** (`RS_AUTO_BACKEND=1 apps cobblestones`, 14 spheres, one
   `Area` constraint per component, Global barycenter): 30 to 60 s per step against 2 s dense.
3. **Obstacle breached in cobblestones**: `apps cobblestones` (also `RS_SLAB_W=30`, `300`, `RS_NO_REMESH=1`): the area growth
   (hard constraint x2.5 in 300 steps) goes on after the spheres jam against the slab; the step size goes to zero
   (`Stuck` at step 320 without remeshing) and vertices end up outside the slab (up to |z| = 1.42 against 0.9), the system
   drifting towards +z (no-remesh run: z range [-0.63, 0.80] at step 20, [-0.91, 1.22] at step 200). CCD guards the line search only;
   the Newton projection onto the area target and the remeshing are not guarded, which is the likely cause (not isolated).
   No sphere crosses another (`intersects` false).
4. **Subcritical flow with remeshing explodes** (`RS_REMESH=1 apps subcritical`): faces 4096 -> 121 760 in 94 steps, energy
   1.5e8, then a 118 GB dense allocation. Without remeshing (used) the intersection count falls from 13 096 to 791 in 600 steps
   (3.5, 7) but does not reach zero.
5. **Linked handcuffs do not unlink** (`apps handcuffs`): see the table; the IH-move of Fig. 18/19 is not found from our
   linked start (a Hopf-type torus-link tube), which has a different geometry from the paper's.
6. **Eversion midsurface is not orientation-symmetric** (`apps eversion`), see the table.
7. `RS_SPINDLE_R=1.005 apps subcritical` stops with "face 31 is degenerate: its area is 0e0" (input generation near r = 1).

## Ishii et al. 2012 (Fig. 16): not done

Approved: arXiv preprint of Ishii, Kishimoto, Moriuchi, Suzuki, "A table of genus two handlebody-knots up to six crossings"
into `data/repulsive-surfaces/ishii2012/`. Searches (web search for the title and authors, arXiv API and search page queries
by title) found **no arXiv version** of this paper: it appears only as a journal reference (J. Knot Theory Ramifications 21(4),
1250035, 2012) in other papers (Lee and Lee, AGT 12, 2012; arXiv 2003.03748; 2511.12194, a seven-crossing extension which is a
different paper and was not downloaded). Nothing was downloaded and no graph was transcribed or invented. `generate/ishii.rs`
was not created. Skipped: all of 4_1, 5_1..5_4, 6_1..6_16. If the user supplies the paper's PDF, `graph_tube(curves, radius,
faces)` takes polylines (over/under crossings as z offsets), as noted in NOTES-WP12.

## Stop rule

The default `EnergyStall { relative: 1e-7, window: 10 }` ended the smoke genus runs after 26 to 80 steps as `Converged`
while the energy was still falling (e.g. genus 11: F 4880 -> 320 in 38 steps, then 2.9e2 at step 36 still decreasing by 2 %
per 3 steps). F is not monotone along accepted steps (new Barnes-Hut lists and remeshing each step), so a window in which F
rose satisfies `F_{k-10} - F_k <= 1e-7 |F_k|`. All applications here use `StopRule::Off` with a fixed step count.

## SPEC proposals

- §14.4: say that Fig. 16 is blocked because the preprint does not exist on arXiv (the user's permission covered only the
  preprint); keep the item open until the PDF is supplied.
- §7.4 / D29: `EnergyStall` should compare windowed minima (or means), not `F_{k-w}` against `F_k`; as written it stops early
  on non-monotone F.
- §6.3: the hierarchical back end diverges with hard pins, on self-intersecting inputs and with many components; document and
  consider a fallback to dense below about 10 000 faces, or fix the preconditioner (failures 1 and 2).
- §7.3 (D20): state that CCD does not guard the constraint projection or the remeshing, so a growing hard area/volume
  constraint can push the surface through an obstacle (failure 3).
