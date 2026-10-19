# Quality against the paper's figures

Why our results looked worse than Figs. 1, 22–24 of Stein, Grinspun, Crane 2018, measured. Short
answer: **the input mesh and the energy of the coarse phase**, not the optimizer or the remeshing.
Run the paper's recipe (a coarse mesh of uniform triangles, `E_P` coarse then `E_λ` fine, 4-1
subdivision between rounds) and the results look like the paper's: smooth developable pieces
separated by regular seam curves. Run it on a subdivision cage (spot's, bob's, blub's
`*_control_mesh.obj`) and thin parts of the cage fold flat into double sheets and spikes. The
"short creases" of the in-place flows on the dense models are mostly how a seam that crosses
triangles looks when every edge above 30 degrees is drawn dark (the showcase's shading).

Decisions that came out of this: SPEC D34 (status `Stationary`), D35 (remeshing fallback), D36
(default flow = `E_P` coarse then `E_λ` fine). The recipe is in the README ("Reproducing the
paper's look").

## Measures (`examples/quality.py`)

Every run is measured against the dense model (`spot_triangulated`, `bob_tri`,
`blub_triangulated`; the crumpled sheet against itself):

* `e30`, `e30_frac`: interior edges with a dihedral angle above 30 degrees (count, fraction).
* `folds`: edges above 170 degrees (a face folded back onto its neighbour; also the rim of a fin
  that has folded flat into two coincident sheets).
* seams: the edges above 20 degrees as a graph; `seam_comp` connected components, `seam_short`
  those of at most 2 edges (isolated short creases), `seam_long_frac` the fraction of the crease
  length in components longer than 5 % of the diagonal (clean curves).
* `d_mean`, `d_haus`: mean and symmetric Hausdorff distance to the dense model in units of its
  bounding-box diagonal (60k area-uniform samples each way; approximate).
* `vol_ratio`: enclosed volume / the dense model's (squashing shows here).
* `min_angle` (degrees), `faces_lt1`, the flow's `bad_faces`, energy, statuses, final `ℓ̄‖g‖∞`.
* a flat-shaded orthographic render (`--png`, `--crease-deg` draws the crease edges).

Scripts: `examples/coarse_inputs.py` (uniform coarse inputs, see below),
`examples/quality_study.sh JOBS OUT` (one `develop` + `quality.py` per job line, in parallel),
`examples/quality_summary.py OUT` (table). `develop` has `--rounds` (`[+]energy:iterations`, `+`
= subdivide first), `--history`, `--memory`, `--gradient-tol`, `--initial-step`, `--no-fallback`.
Logs, meshes and renders of this study: `runs_tmp/dev-quality/s1` … `s6` (not committed).

**Inputs.** The paper's coarse inputs come from Instant Meshes (uniform, triangle field,
extrinsically aligned; about 1k triangles for Fig. 1). We have no Instant Meshes; MeshLab's
isotropic explicit remeshing (`pymeshlab`, only to make the inputs) of the dense models to
about 1k triangles gives `spot_iso1000` (1102 triangles), `bob_iso1000` (934), `blub_iso1000`
(1058); also 400 and 2500. The control meshes are 372 (spot), 668 (bob) and 220 (blub) triangles
after fan triangulation of the quads.

## Headline: old default vs the paper's recipe (500 iterations per round, 2 subdivisions)

`old` = the previous default (`E_λ` in all three rounds, no remeshing fallback), `new` = the new
default (`E_P`, then `E_λ` after each subdivision, fallback on), `fig1` = `E_P` in all rounds.
Dense = the flow in place on the dense model, `E_λ`, 4000 iterations, no subdivision.

| run | faces | e30 frac | folds | seams (short) | long frac | d_mean | d_haus | vol | bad |
|---|---|---|---|---|---|---|---|---|---|
| spot control, old | 4882 | 0.075 | 282 | 5 (2) | 0.998 | 0.055 | 0.328 | 0.51 | 12 |
| spot control, new | 5138 | 0.056 | 168 | 16 (4) | 0.989 | 0.050 | 0.221 | 1.15 | 0 |
| spot uniform 1k, old | 16678 | 0.039 | 99 | 16 (2) | 0.994 | 0.011 | 0.096 | 0.91 | 0 |
| spot uniform 1k, **new** | 16854 | 0.035 | 24 | 31 (4) | 0.968 | 0.007 | **0.045** | **0.96** | 0 |
| spot uniform 1k, fig1 | 17066 | 0.033 | 19 | 28 (0) | 0.993 | 0.007 | 0.045 | 0.97 | 0 |
| spot dense in place | 5530 | 0.070 | 8 | 35 (20) | 0.970 | 0.005 | 0.028 | 0.98 | 0 |
| bob control, old | 9842 | 0.029 | 57 | 16 (7) | 0.981 | 0.024 | 0.108 | 0.87 | 0 |
| bob control, new | 10210 | 0.034 | 4 | 5 (0) | 1.000 | 0.020 | 0.084 | 1.06 | 0 |
| bob uniform 1k, old | 14110 | 0.033 | 31 | 8 (1) | 0.993 | 0.017 | 0.085 | 0.83 | 0 |
| bob uniform 1k, **new** | 14346 | 0.035 | 0 | 31 (7) | 0.974 | 0.009 | **0.058** | **0.97** | 0 |
| bob uniform 1k, fig1 | 14512 | 0.035 | 0 | 13 (0) | 0.997 | 0.009 | 0.058 | 0.97 | 0 |
| bob dense in place | 10518 | 0.026 | 3 | 41 (6) | 0.985 | 0.003 | 0.022 | 0.99 | 0 |
| blub control, old | 2866 | 0.086 | 287 | 1 (0) | 1.000 | 0.047 | 0.324 | 0.08 | 2 |
| blub control, new | 2966 | 0.067 | 134 | 6 (0) | 0.995 | 0.024 | 0.142 | 0.70 | 0 |
| blub uniform 1k, old | 16440 | 0.022 | 264 | 10 (3) | 0.975 | 0.009 | 0.050 | 0.84 | 0 |
| blub uniform 1k, **new** | 16430 | 0.022 | 185 | 12 (3) | 0.979 | 0.008 | **0.053** | **0.89** | 0 |
| blub uniform 1k, fig1 | 16508 | 0.021 | 169 | 14 (1) | 0.988 | 0.008 | 0.049 | 0.89 | 0 |
| blub dense in place | 13912 | 0.024 | 6 | 59 (20) | 0.910 | 0.002 | 0.016 | 0.99 | 0 |
| crumpled, `E_λ` 500, old / new | 501 / 502 | 0.023 / 0.022 | 1 / 3 | 10 / 9 | 1.0 | 0.009 | 0.050 / 0.054 | – | 0 |

Reading: the control meshes are the wrong input (Hausdorff 0.14 to 0.33 of the diagonal even with
the paper's energies); on uniform 1k inputs the paper's recipe keeps the volume (0.89 to 0.97)
and halves the Hausdorff distance of the old default on spot and bob. Blub's folds are its fins:
they are thin in the model and fold flat into two coincident sheets (a hinge of 180 degrees is a
zero-energy star for every energy). Renders: `runs_tmp/dev-quality/m6.png` (control old/new, uniform
new), `m1.png`, `m2.png`, `m3.png` (crease overlays).

## Hypotheses, ranked by the evidence

1. **Input mesh (strong).** Same recipe, control mesh vs uniform 1k: spot Hausdorff 0.22 -> 0.045,
   volume 1.15 -> 0.96, folds 168 -> 24; bob 0.084 -> 0.058; blub 0.142 -> 0.053. On the cages,
   thin parts (spot's ears and head, blub's fins, bob's tail) become flat double sheets and
   spikes for every energy and budget tried (50, 100, 300, 1000, 3000 it/round; step caps 0.05 to
   0.5; initial step 0.01; `runs_tmp/dev-quality/s3`): fewer iterations only make it less
   pronounced (spot `E_P` 50 it/round: Hausdorff 0.115, 300: 0.20). Coarse uniform inputs of 400
   and 2500 triangles behave like 1k (Hausdorff spot 0.051 / 0.052, blub 0.051 / 0.025), except
   that 400 is too coarse for bob (0.10, volume 0.71; 2500: 0.052).
2. **Energy of the coarse phase (strong).** `E_λ` keeps flattening coarse meshes as it converges:
   one round of 5000 iterations on the uniform 1k inputs ends at volume 0.70 / 0.30 / 0.50 (spot,
   bob, blub) with `E_λ` and 0.99 / 0.92 / 0.96 with `E_P`; on the blub control mesh 3 x 1000
   iterations of `E_λ` end at volume 1e-4 (a flat pancake, E = 1.5e-7). The flat pancake is a
   minimizer (every star a hinge of 180 degrees, D32). `E_P` does not do this as fast. The paper
   uses `E_P` for Figs. 1, 2, 25 and `E_P` coarse, `E_λ` fine for Figs. 23, 24. -> D36.
   The intrinsic energy (Figs. 20–22) gets stuck at stiff states: `ℓ̄‖g‖∞` 4e2 to 7e7, the energy
   rises like `1e15 |Δx|²` down to displacements of 1e-9 (a star whose vector area nearly
   vanishes); bob uniform 1k with intrinsic made 0 iterations in rounds 1 and 2.
3. **The dense in-place flow lacks the coarse phase (inherent).** Sec. 4.2: the energy of a seam
   shrinks with the edge length, so the coarse mesh is what regularizes the seams. In place on the
   dense models: 35 to 59 crease components of which 6 to 20 are isolated short creases, against
   12 to 31 (3 to 7 short) for the coarse-to-fine flow at more faces; but 91 to 98 % of the crease
   length is in long curves either way. The showcase's look ("covered in short creases") comes
   mostly from its shading: it darkens every edge above 30 degrees, and a seam that crosses the
   triangles diagonally is a zigzag of such edges, drawn as hatching. Flat-shaded, the same mesh
   shows pieces and seam curves (`runs_tmp/dev-quality/z1.png`, `z2.png`: spot, 6000 iterations,
   the showcase run). 180-degree folds on the dense flows: 3 to 18.
4. **Remeshing threshold (no).** 5, 10, 15 degrees instead of 1: on the coarse-to-fine flows they
   collapse many edges (spot 16.7k -> 9.2k faces at 15 degrees) and the Hausdorff distance doubles
   (0.049 -> 0.10 to 0.13; bob at 15 degrees: volume 0.60); folds do not go down. On the dense
   in-place flows 10 degrees halves spot's crease components (46 -> 21) but not bob's (48 -> 47).
   Default stays 1 degree. Remeshing off changes the coarse round little (folds 44 vs 48 on blub)
   but leaves 14 to 28 faces below 1 degree.
5. **Remeshing legality (minor; D35).** The faces left below the small angle (`bad_faces`) are
   refused collapses (link condition, most of them) and flips (normal turn, opposite vertices
   adjacent). Trying the other operation of the triangle halves them over the 45 robustness flows
   (62 -> 35). On the uniform inputs there are none either way, and the admissibility floor of D31
   never engages (0 inadmissible trials, smallest angle above 1 degree), so blocked vertices do
   not freeze folds there.
6. **Optimizer (no).** The step cap of 0.5 ℓ̄ never acts on the uniform inputs (0.2 gives bitwise
   the same flows); 0.1 and 0.05 change the results only chaotically (Hausdorff 0.049 / 0.051 /
   0.039 on spot, folds 27 / 24 / 12, blub 196 / 181 / 164). L-BFGS memory 3, 10, 30: no
   systematic effect (spot Hausdorff 0.044 / 0.045 / 0.052). The gradient never becomes "small":
   `ℓ̄‖g‖∞` stays 0.02 to 3 for `E_P` and `E_λ` on the uniform inputs after 1000 to 5000 iterations (nonsmooth energies;
   seam vertices sit on kinks), so "minimize until the gradient is sufficiently small" is not
   operational; the energy drops most in the first 50 to 200 iterations of a round (spot `E_P`
   32 -> 9.8 in 50, 5.6 at 1000, 2.5 at 5000), and 300, 1000, 3000 iterations per round give the
   same look (Hausdorff grows slowly with the budget: 0.043, 0.049, 0.069 on spot).
7. **Subdivision (fine).** Midpoint 4-1 preserves the energy (tested, SPEC §5.2), and every
   subdivided round keeps decreasing it (spot `E_P`: 5.6 -> 3.6 -> 2.3).
8. **Valence-3 omission, spikes under `E_λ`.** Not separately varied. On the uniform inputs `E_λ`
   alone makes spot's legs into spikes (old rows above); `E_P` first keeps them as wedges.

## The converged status (D34)

A flow run to convergence used to end with `LineSearchFailed` at energies of 1e-7 to 1e-13. The
line search from steepest descent fails there because the energy is flat to rounding: on the
over-runs of the control meshes (3 x 3000 iterations) the `E_λ` failures have `|ΔE| = 4.7e-15` at
the smallest trial against a rounding noise of `3e-11` (E 5.3e-9, `ℓ̄‖g‖∞` 1.3e-6), while the
intrinsic failures have `|ΔE|` of 4e-3 to 1.1 at gradients of 4e5 to 7e7 (a stiff point, not a
minimum) or hit a wall. Those are now `Stationary` and `LineSearchFailed`, respectively.

## Open

* Thin parts of a cage become flat double sheets for every energy; a cage is not a good input
  unless it is regular and has no thin parts (the paper's Fig. 23 cages). An Instant-Meshes-like
  remesher in the workspace would let the Blender add-on run the paper's workflow on any model.
* The intrinsic energy's stiff states (vector area nearly zero) stop its flows (D12 threshold
  `π − φ < 1e-6`); not changed here.
* (Done) The Python and Blender layers offer the recipe (`coarse_energy`, default `E_P` in round
  0). Blender's Remesh (Voxel) + Decimate to about 1k triangles gives inputs as good as the
  MeshLab ones (`blender/README.md`, "Starting mesh").
