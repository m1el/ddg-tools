# Comparison with the authors' program (black box)

This file compares `ddg-developability` with the reference program of Stein, Grinspun & Crane,
*Developability of Triangle Meshes* (ACM TOG 2018). It was written on 2026-10-19.

**Clean room.** The reference (GPL-2.0, with libigl and Kopp's LGPL eigensolver) was built and
run as a black box only. None of its sources, headers, CMake files, scripts or run files were
opened. Every statement below about "their" energies, defaults or behaviour is an
**observation**: a value shown in its window, a file it wrote, its exit status, or its stdout.
None of it describes their code. Nothing from the reference was copied into this repository.

- Versions:
  - Reference: commit `67848c1` (2018-11-15), with libigl `d4b6726` (v1.2.1, the pinned
    submodule; fetched by the user).
  - Ours: `ddg-developability` at `01eb71f` (includes D34/D35). An earlier pass at `6a259bf` differs only in
    10 remeshed runs (D35 changed the remeshing), with the same overall picture.
- Machine: 64 cores, shared with other jobs. `uptime` load was 1–37 during the parallel
  phase and 0.9–3.7 during the timing runs.
- Work directory: `runs_tmp/dev-blackbox/` (gitignored). It holds all scripts, logs, meshes,
  screenshots and the CSV of every measurement (`an/results.csv`).

## Summary

| question | answer |
|---|---|
| Their E_λ (energy 2), intrinsic L2 (0) and intrinsic L∞ (1) on identical geometry | **Equal to ours** to every displayed digit (6 significant), once their boundary rule is applied (next row). |
| Boundary vertices | **They leave every boundary vertex out**; we count them by default (D8, `include_boundary`). With `include_boundary = false` our three energies above reproduce theirs. |
| Interior valence-3 vertices | Left out by both programs. |
| Their L∞ covariance (energy 3) | **Exactly ½ × our `covariance-max`.** |
| Their L∞ combinatorial (energy 5) | **½ × our π_max, with arcs of one face not allowed** (see §2). |
| Their L2 combinatorial (energy 4) | **A different normalization from ours (D1).** Per star it is `min over 2-arc splits with both arcs ≥ 2 faces of Σ_{pairs inside the arcs} ‖N_i − N_j‖² / (2 · #pairs)`. Ours is `Σ_p (1/n_p) Σ_{σ∈F_p} ‖N_σ − N̄_p‖²`, single-face arcs allowed. |
| Flow, L-BFGS, 200 steps, same energy, measured with our evaluator | Ours is lower at step 10 in all 36 E_λ / intrinsic runs. At step 200 ours is lower in 25 of 30 runs where both finished: 15/16 intrinsic, 10/15 E_λ. On the L∞ covariance energy the reference is lower in every finished run (5 of 5). |
| Robustness, 72 reference runs | **The reference failed in 8 runs**: NaN positions (4), a degenerate face (2), duplicate faces after its remeshing (1), a crash (1). Its `f` never returned on any mesh with boundary (12 of 12). Ours: 72 of 72 runs completed with finite results and no invalid mesh. |
| Time per L-BFGS iteration, dense meshes (5.9k–14k faces) | Reference through its GUI: 40–160 ms per step, at most 60 steps/s on small meshes (frame bound). Ours: 4–14 ms on 1 thread, 2–6 ms on 8 threads. |
| Cut and flatten (`f`) vs ours (cut + BFF, zero boundary scale) | Their layouts have much larger area distortion: median mean ‖u − ū‖ 1.1 on the coarse meshes and 0.24 on the dense ones; ours 5e-5 and 7e-4. On the coarse meshes their layouts flip 61 triangles in 17 of 28 layouts; ours flips none with the defect rule (D30). Their cuts are longer: median 0.26 / 0.19 of the total edge length, ours 0.20 / 0.14. Their printed "Flattening error" is the max and unweighted mean of the quasi-conformal error. |
| Suspected problems in ours | **Our L∞ covariance flow makes slow progress** (§3.3, reproducer). The `develop` example's `RESULT energy=` field is wrong under `--rounds`. No incorrect values were found. |

## 1. Building and driving the reference

### Build

The reference was built out of tree in `/var/data/ddg-quarantine/build-dev`. The source
directory was opened (`chmod 755`) only inside trap'd scripts and is left at mode 000. No
source was edited. Steps:

1. `cmake` 3.31.8 (micromamba prefix) and gcc 15.2. The configure command:
   `cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-include /var/data/ddg-quarantine/prefix/shim/glew_shim.h" -DCMAKE_CXX_STANDARD_LIBRARIES=-lGLU <ref>`.
   It finds libigl at `../libigl`. `cmake -LH` already lists `LIBIGL_WITH_MATLAB=OFF`, and the
   program builds and links without MATLAB.
2. Two Linux-only problems, both fixed by flags:
   - The compile error "gl.h included before glew.h".
   - An undefined `gluErrorString` at link time, fixed by `-lGLU`.
3. Shim (approved by the user). The program first built, but crashed at start-up with
   `pc = 0`, a call through a null GL function pointer. The force-included shim
   `glew_shim.h` sets `GL_GLEXT_PROTOTYPES`, includes the system `GL/gl.h` and `glext.h`,
   and pre-sets glew's include guard. GL calls then link directly to libGL.
4. `make -j16` then gives `applications/main/main`.

Scripts: `runs_tmp/dev-blackbox/scripts/{configure,build}.sh`.

### Driving the program

- Each job runs under its own `Xvfb` (its PID recorded and killed). Keys are sent with
  XTEST (`scripts/xa.py`).
- After every batch of keys the driver presses `w` and waits for `out.obj`. Keys are handled
  in order, so the file proves that the preceding steps are done. Screenshots of the text
  panel give the displayed values.
- Driver: `scripts/probe.sh DISPLAY MESH OUTDIR KEY[*N][:s|o] ...`. A matrix job:
  `scripts/ref_job.sh DISPLAY MESH EIDX def|n`, i.e. `m`, then `e`×EIDX, optional `n`, then
  TAB×10, ×40, ×150 with a sync after each batch, then `f`.

Behaviour observed from the window and the outputs:

- Default state:
  - energy 0. Its displayed value equals our intrinsic energy on the written mesh.
  - method: a fixed step. "Timestep: 1e-05" stays constant, and the energy can rise.
  - remeshing on. Without `n` the vertex count changes and stdout prints "... edges
    collapsed, ... flipped. A structural change happened to the mesh."; with `n` it never
    changes.
- `m` once gives L-BFGS. Its timesteps double from 1e-5 by powers of two, as an Overton-type
  weak Wolfe search does. `m` twice gives a non-dyadic step, consistent with backtracking;
  `m` three times returns to the fixed step. All comparisons use `m` once.
- The panel shows:
  - "Energy (E)" and "Scaled L2 gradient norm" of the current state;
  - "Timestep", "total time" (sum of step sizes) and "total steps";
  - "Smallest triangle area" and "smallest triangle angle". The angle equals our minimum
    corner angle on the same file (0.00294921 rad on `spot_ctrl`, ours 2.949e-3).
  - "Hausdorff dist to original" (not compared).
- `f` writes `flattened_out.obj` and `mesh_with_flattened_texcoord.obj`. It prints
  "Flattening error, max: A, avg: B". **A and B are the maximum and the unweighted mean of
  the quasi-conformal error σ_max/σ_min of the layout**: equal to all printed digits on 4 of
  4 checked layouts (`scripts/qc_unweighted.py`). The layout is globally rescaled (mean
  log scale ≈ −0.7 … −3.6).

## 2. Energy values on identical geometry

Method:

1. Run their program for one or more L-BFGS steps and read the displayed energy.
2. Sync with `w` and evaluate our energies on the written `out.obj`, i.e. the same state.
3. Find a definition, built from our per-vertex terms or from a small NumPy re-evaluation
   (`scripts/variants*.py`, `scripts/hyp_ep.py`), that reproduces the displayed value.

Single-star test meshes isolate the per-vertex term. Each is a fan with one interior vertex
of valence k = 3…10 and random heights (`meshes/fan*.obj`). Since boundary vertices do not
count, the displayed value is that one vertex's energy.

Data: spot control mesh, all six energies; noisy open cylinder (boundary), energies 0, 2, 4;
14 fans (k = 4…10), all six energies on 6 fans and energies 4, 5 on 8 more; fan k = 3; final
states of 18 runs of the matrix (dense meshes, boundary meshes, remeshing on and off, all six
energies). With the rules below **every one of these values is reproduced to all displayed
digits**.

| their `e` index (README name) | our energy | observed relation |
|---|---|---|
| 0 L2 covariance "modified to reduce spikes" | `Intrinsic` | equal (interior vertices only, valence 3 omitted) |
| 1 L∞ of the same | `IntrinsicMax` | equal (same domain) |
| 2 L2 covariance | `Covariance` (E_λ) | equal (same domain) |
| 3 L∞ covariance | `CovarianceMax` | **their value = ½ × ours** (same domain) |
| 4 L2 combinatorial | `Combinatorial` (E_P) | **different normalization**, see below |
| 5 L∞ combinatorial | `CombinatorialMax` (π_max) | **½ × ours, and both arcs must have ≥ 2 faces** |

Domain:

- **Boundary vertices: left out by them.** Cylinder, intrinsic: shown 22.0069; our sum over
  interior vertices 22.006939; ours with boundary 22.835868. Covariance: shown 20.9509, our
  interior sum 20.950908.
- **Interior valence-3 vertices: left out by both.** The k = 3 fan shows 0 for all six
  energies. On spot the matches hold only with valence 3 omitted.

Energy 4, observed per-star definition. Over the contiguous 2-arc splits `{F_1, F_2}` of the
star with `|F_1|, |F_2| ≥ 2`:

```
e_v = min  ( Σ_{p=1,2} Σ_{σ<τ ∈ F_p} ‖N_σ − N_τ‖² ) / ( 2 · Σ_p C(|F_p|, 2) )
```

This is half the mean over within-arc pairs, pooled over both arcs. Ours (Eq. 1 first line,
D1) is `Σ_p (1/n_p) Σ ‖N − N̄_p‖² = Σ_p (1/n_p²) Σ_{pairs} ‖N_σ − N_τ‖²`, with single-face
arcs allowed. The two agree only for 2 | 2 splits (k = 4). For a k/2 | k/2 split the ratio
theirs/ours is `m / (2(m−1))`: 0.75 for m = 3, ⅔ for m = 4, as observed. For other valences
the minimizing split can differ.

Checks:

- spot after one step: shown 46.0857; this formula gives 46.085704; our E_P 43.301227.
- cylinder: shown 6.35269; formula 6.352689; our E_P (interior only) 7.106113.
- With single-face arcs allowed instead, the fans still match but spot does not (43.13).

Energy 5: `e_v = ½ · min over splits with both arcs ≥ 2 faces of max_{σ,τ in one arc} ‖N_σ − N_τ‖²`.

- spot: shown 104.761; this formula 104.760563; ½ × our π_max 92.59 (ours allows single-face
  arcs).
- On every fan whose optimal split has no single-face arc the ratio to ours is exactly
  0.500000 (13 fans). On the k = 9 fan, whose optimum has one, it is 0.645, and the
  restricted formula reproduces the shown 0.433508.

Notes:

- The factor ½ on energies 3 and 5 changes no minimizer.
- The arc-size restriction and the pooled normalization of energy 4 do change the energy
  landscape. Our E_P and theirs are different objectives. That is visible in §3: each flow
  wins on its own definition of E_P.
- Reproduce:
  - energies 0–3: `devbb eval MESH PREFIX`, then sum the per-vertex table over
    `boundary == 0 && valence != 3`, with ½ for `covariance-max`;
  - energies 4–5: `scripts/obs_defs.py MESH`.

## 3. Flow outcomes

Setup. Both programs run L-BFGS from the same triangle OBJ.

- Meshes: the spot, bob and blub control meshes, fan-triangulated: `*_ctrl`, 220–668 faces.
  Their dense triangulations: `*_tri`, 5.9k–14k faces. `CrumpledDevelopable` (a disk). A
  noisy open cylinder (2 boundary loops). A level-2 icosphere.
  Script: `scripts/meshes.py`.
- Reference: `m`, energy chosen with `e`, remeshing on (default) or off (`n`), 200 TAB steps,
  each one L-BFGS iteration ("total steps" counts them).
- Ours: one round without subdivision, 200 iterations, all early stops disabled.
  `--rounds ENERGY:200 --gradient-tol 0 --stall-iters 100000000 --threads 2 [--no-remesh]`.
  Script: `scripts/ours_job.sh`. All 72 of our runs ended with `max_iter`.

Columns:

- Energy values are measured with **our** evaluator and **our** definition of the same energy:
  `E` columns; "their def." columns use the rules of §2.
- Displacement is in mean edge lengths of the input: per vertex if the topology is unchanged,
  else nearest vertex.
- `V` = vertex count after 200 steps (it changes with remeshing).
- "NaN pos." = the reference wrote non-finite positions.
- Remeshing "on" = their default and our `RemeshOptions` default.

Pooled over these tables and the L∞ ones (§3.2), counting runs where both finished:

| energy | step 10: ours lower / ref lower | step 50 | step 200 | reference failed by step 200 |
|---|---|---|---|---|
| intrinsic (0) | 18 / 0 | 15 / 3 | 15 / 1 | 2 |
| E_λ (2) | 18 / 0 | 14 / 3 | 10 / 5 | 3 |
| E_P, our definition (4) | 16 / 2 | 11 / 7 | 11 / 6 | 1 (crash) |
| E_P, their definition (4) | 10 / 8 | 4 / 14 | 2 / 15 | |
| intrinsic L∞ (1) | 1 / 5 | 1 / 5 | 2 / 4 | 0 |
| covariance L∞ (3) | 6 / 0 | 0 / 6 | 0 / 5 | 1 |
| combinatorial L∞ (5), ours / theirs | 6 / 0 · 2 / 4 | 3 / 3 · 0 / 6 | 4 / 2 · 0 / 6 | 0 |

Observations:

- Ours descends much faster in the first iterations. Its first trial moves the fastest vertex
  by 0.1 ℓ̄ (D17). Their L-BFGS starts from a step of 1e-5 and doubles.
- Late in the run the two are comparable on intrinsic and E_λ. The dense meshes all end
  within about 15 %; ours is lower in 11 of 12 intrinsic/E_λ dense runs, the exception being
  spot_tri E_λ without remeshing, 17.5 vs 20.1.
- For E_P each program wins on its own definition. This is expected, since they minimize
  different functions (§2).
- On the boundary meshes both reach near-zero energy. Ours is lower at step 200 in every one
  of these 12 runs (our definition).

### 3.1 Intrinsic, E_λ and E_P (their energies 0, 2, 4)

**Energy 0 of theirs = our `intrinsic`**

| mesh | remesh | E input (ours) | ref E@10 | ours E@10 | ref E@50 | ours E@50 | ref E@200 | ours E@200 | ref E@200 (their def.) | ours E@200 (their def.) | ref disp mean/max | ours disp mean/max | ref V | ours V |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| spot_ctrl | on | 175.9 | 169.9 | 58.98 | 24.39 | 25.48 | NaN pos. | 12.61 | NaN pos. | 12.61 | — | 0.39 / 1.19 | NaN pos. | 175 |
| spot_ctrl | off | 175.9 | 169.9 | 121.4 | 24.39 | 22.44 | NaN pos. | 10.92 | NaN pos. | 10.92 | — | 0.52 / 2.41 | NaN pos. | 188 |
| bob_ctrl | on | 67.72 | 41.06 | 29.43 | 21.77 | 14.36 | 11.39 | 7.499 | 11.39 | 7.499 | 0.28 / 0.93 | 0.41 / 1.38 | 333 | 325 |
| bob_ctrl | off | 67.72 | 41.06 | 29.43 | 23.47 | 14.12 | 8.466 | 9.504 | 8.466 | 9.504 | 0.44 / 1.84 | 0.43 / 1.34 | 334 | 334 |
| blub_ctrl | on | 89.01 | 58.76 | 26.32 | 17.55 | 7.018 | 8.801 | 4.488 | 8.801 | 4.488 | 0.17 / 0.47 | 0.24 / 0.61 | 107 | 98 |
| blub_ctrl | off | 89.01 | 58.76 | 30.73 | 16.28 | 5.665 | 6.689 | 4.605 | 6.689 | 4.605 | 0.30 / 0.96 | 0.37 / 0.91 | 112 | 112 |
| spot_tri | on | 117.8 | 77.56 | 72.1 | 37.69 | 41.92 | 32.37 | 30.4 | 32.37 | 30.4 | 0.11 / 0.67 | 0.12 / 0.67 | 2923 | 2917 |
| spot_tri | off | 117.8 | 77.56 | 72.1 | 37.7 | 41.53 | 31.7 | 27.91 | 31.7 | 27.91 | 0.13 / 0.77 | 0.16 / 0.88 | 2930 | 2930 |
| bob_tri | on | 43.82 | 36.33 | 33.35 | 23.45 | 21.89 | 15.99 | 14.28 | 15.99 | 14.28 | 0.08 / 0.62 | 0.09 / 0.84 | 5343 | 5341 |
| bob_tri | off | 43.82 | 36.33 | 33.35 | 23.45 | 21.89 | 15.43 | 13.46 | 15.43 | 13.46 | 0.09 / 1.01 | 0.13 / 1.22 | 5344 | 5344 |
| blub_tri | on | 48.36 | 41.6 | 36.31 | 24.26 | 18.57 | 12.78 | 12.1 | 12.78 | 12.1 | 0.06 / 0.76 | 0.07 / 0.82 | 7100 | 7092 |
| blub_tri | off | 48.36 | 41.6 | 36.31 | 24.26 | 18.57 | 12.4 | 11.71 | 12.4 | 11.71 | 0.08 / 0.93 | 0.08 / 1.01 | 7106 | 7106 |
| crumpled | on | 79.59 | 8.601 | 2.041 | 0.2446 | 0.08936 | 0.01441 | 0.002251 | 0.01245 | 0.002207 | 0.34 / 0.74 | 0.35 / 0.79 | 289 | 287 |
| crumpled | off | 79.59 | 8.601 | 2.041 | 0.2446 | 0.07555 | 0.006747 | 0.005209 | 0.004875 | 0.005134 | 0.37 / 1.02 | 0.37 / 0.93 | 289 | 289 |
| cylinder_noisy | on | 23.41 | 1.697 | 0.6488 | 0.007664 | 0.002681 | 1.548e-05 | 9.305e-06 | 1.487e-05 | 9.149e-06 | 0.08 / 0.21 | 0.08 / 0.24 | 416 | 416 |
| cylinder_noisy | off | 23.41 | 1.697 | 0.6488 | 0.007664 | 0.002681 | 1.548e-05 | 9.305e-06 | 1.487e-05 | 9.149e-06 | 0.08 / 0.21 | 0.08 / 0.24 | 416 | 416 |
| icosphere2 | on | 14.3 | 6.377 | 4.972 | 3.855 | 2.783 | 3.855 | 2.27 | 3.855 | 2.27 | 0.31 / 0.53 | 0.63 / 0.76 | 162 | 151 |
| icosphere2 | off | 14.3 | 6.377 | 4.972 | 3.855 | 2.783 | 3.855 | 2.218 | 3.855 | 2.218 | 0.31 / 0.53 | 1.15 / 2.23 | 162 | 162 |

**Energy 2 of theirs = our `covariance`**

| mesh | remesh | E input (ours) | ref E@10 | ours E@10 | ref E@50 | ours E@50 | ref E@200 | ours E@200 | ref E@200 (their def.) | ours E@200 (their def.) | ref disp mean/max | ours disp mean/max | ref V | ours V |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| spot_ctrl | on | 99.74 | 71.26 | 27.24 | 4.893 | 7.89 | 0.8642 | 2.724 | 0.8642 | 2.724 | 0.74 / 2.32 | 0.52 / 2.19 | 183 | 171 |
| spot_ctrl | off | 99.74 | 87.12 | 31 | 3.345 | 7.195 | 1.486 | 4.106 | 1.486 | 4.106 | 1.62 / 5.03 | 0.74 / 4.04 | 188 | 188 |
| bob_ctrl | on | 47.2 | 33.56 | 23.64 | 9.744 | 9.619 | 5.295 | 4.374 | 5.295 | 4.374 | 0.47 / 1.64 | 0.44 / 1.45 | 331 | 325 |
| bob_ctrl | off | 47.2 | 33.56 | 23.64 | 9.744 | 9.684 | NaN pos. | 5.493 | NaN pos. | 5.493 | — | 0.52 / 1.76 | NaN pos. | 334 |
| blub_ctrl | on | 37.58 | 17.34 | 7.432 | NaN pos. | 1.738 | NaN pos. | 0.2446 | NaN pos. | 0.2446 | — | 0.37 / 1.96 | NaN pos. | 95 |
| blub_ctrl | off | 37.58 | 17.34 | 5.396 | 0.5613 | 0.8227 | NaN pos. | 0.2104 | NaN pos. | 0.2104 | — | 0.59 / 2.25 | NaN pos. | 112 |
| spot_tri | on | 112 | 81.95 | 66.01 | 33.53 | 32.45 | 24.17 | 22.03 | 24.17 | 22.03 | 0.14 / 0.82 | 0.16 / 0.85 | 2928 | 2913 |
| spot_tri | off | 112 | 81.95 | 66.01 | 31.05 | 30.7 | 17.46 | 20.13 | 17.46 | 20.13 | 0.25 / 1.54 | 0.21 / 1.29 | 2930 | 2930 |
| bob_tri | on | 42.91 | 36.99 | 32.45 | 23.26 | 21.58 | 15.26 | 14.64 | 15.26 | 14.64 | 0.08 / 0.78 | 0.08 / 0.78 | 5342 | 5339 |
| bob_tri | off | 42.91 | 36.99 | 32.45 | 23.26 | 21.58 | 14.09 | 12.46 | 14.09 | 12.46 | 0.12 / 1.25 | 0.14 / 1.42 | 5344 | 5344 |
| blub_tri | on | 45.62 | 38.45 | 34.71 | 19.65 | 17.47 | 11.5 | 11.18 | 11.5 | 11.18 | 0.06 / 0.81 | 0.06 / 0.82 | 7096 | 7091 |
| blub_tri | off | 45.62 | 38.45 | 34.71 | 19.65 | 17.47 | 8.654 | 8.538 | 8.654 | 8.538 | 0.10 / 1.11 | 0.10 / 1.39 | 7106 | 7106 |
| crumpled | on | 61.94 | 7.754 | 2.661 | 0.3773 | 0.0932 | 0.056 | 0.005389 | 0.04476 | 0.005338 | 0.32 / 0.73 | 0.33 / 0.75 | 289 | 284 |
| crumpled | off | 61.94 | 7.754 | 1.988 | 0.198 | 0.09272 | 0.01455 | 0.006862 | 0.01256 | 0.006706 | 0.37 / 0.95 | 0.39 / 1.08 | 289 | 289 |
| cylinder_noisy | on | 22.84 | 1.676 | 0.7137 | 0.003474 | 0.002674 | 1.06e-05 | 1.008e-05 | 1.029e-05 | 9.981e-06 | 0.07 / 0.20 | 0.08 / 0.22 | 416 | 416 |
| cylinder_noisy | off | 22.84 | 1.676 | 0.7137 | 0.003474 | 0.002674 | 1.06e-05 | 1.008e-05 | 1.029e-05 | 9.981e-06 | 0.07 / 0.20 | 0.08 / 0.22 | 416 | 416 |
| icosphere2 | on | 14.12 | 6.52 | 4.754 | 2.559 | 2.516 | 1.806 | 2.139 | 1.806 | 2.139 | 0.74 / 2.74 | 0.52 / 0.86 | 159 | 161 |
| icosphere2 | off | 14.12 | 6.52 | 4.754 | 2.559 | 2.516 | 1.916 | 2.21 | 1.916 | 2.21 | 0.80 / 2.74 | 0.83 / 2.10 | 162 | 162 |

**Energy 4 of theirs = our `combinatorial`**

| mesh | remesh | E input (ours) | ref E@10 | ours E@10 | ref E@50 | ours E@50 | ref E@200 | ours E@200 | ref E@200 (their def.) | ours E@200 (their def.) | ref disp mean/max | ours disp mean/max | ref V | ours V |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| spot_ctrl | on | 43.49 | 33.74 | 19.18 | 9.478 | 6.95 | 7.85 | 3.588 | 6.376 | 8.708 | 0.35 / 1.07 | 0.44 / 1.90 | 184 | 180 |
| spot_ctrl | off | 43.49 | 36.71 | 22.56 | 8.687 | 9.036 | 7.954 | 6.766 | 7.299 | 11.96 | 0.40 / 1.13 | 0.37 / 1.14 | 188 | 188 |
| bob_ctrl | on | 16.29 | 12.67 | 10.37 | 5.285 | 5.984 | 3.031 | 4.229 | 2.304 | 5.427 | 0.46 / 2.64 | 0.29 / 1.13 | 333 | 328 |
| bob_ctrl | off | 16.29 | 12.67 | 10.37 | 5.285 | 5.966 | 3.208 | 4.155 | 2.424 | 5.702 | 0.51 / 2.71 | 0.35 / 1.46 | 334 | 334 |
| blub_ctrl | on | 23.74 | 16.93 | 10.22 | 4.3 | 4.049 | — | 3.283 | — | 9.465 | — | 0.31 / 0.92 | — | 109 |
| blub_ctrl | off | 23.74 | 16.93 | 10.22 | 4.3 | 4.049 | 3.999 | 3.398 | 3.004 | 11.25 | 0.34 / 0.82 | 0.38 / 1.10 | 112 | 112 |
| spot_tri | on | 39.08 | 32.01 | 28.87 | 17.63 | 18.51 | 12.61 | 13.76 | 8.633 | 15.85 | 0.15 / 0.82 | 0.11 / 0.62 | 2930 | 2927 |
| spot_tri | off | 39.08 | 32.01 | 28.87 | 17.63 | 18.51 | 12.61 | 13.49 | 8.633 | 16.08 | 0.15 / 0.82 | 0.13 / 0.91 | 2930 | 2930 |
| bob_tri | on | 13.57 | 12.49 | 11.56 | 9.306 | 8.275 | 6.237 | 5.824 | 4.598 | 5.151 | 0.08 / 1.00 | 0.07 / 0.97 | 5344 | 5344 |
| bob_tri | off | 13.57 | 12.49 | 11.56 | 9.306 | 8.275 | 6.237 | 5.824 | 4.598 | 5.151 | 0.08 / 1.00 | 0.07 / 0.97 | 5344 | 5344 |
| blub_tri | on | 19.61 | 17.85 | 16.94 | 13.4 | 11.99 | 8.687 | 6.51 | 6.598 | 7.315 | 0.05 / 0.61 | 0.07 / 1.07 | 7106 | 7106 |
| blub_tri | off | 19.61 | 17.85 | 16.94 | 13.4 | 11.99 | 8.687 | 6.51 | 6.598 | 7.315 | 0.05 / 0.61 | 0.07 / 1.07 | 7106 | 7106 |
| crumpled | on | 33.9 | 5.923 | 2.379 | 0.3423 | 0.1023 | 0.05751 | 0.01407 | 0.00334 | 0.03034 | 0.40 / 0.96 | 0.37 / 1.07 | 289 | 289 |
| crumpled | off | 33.9 | 5.923 | 2.379 | 0.3423 | 0.1023 | 0.05751 | 0.01407 | 0.00334 | 0.03034 | 0.40 / 0.96 | 0.37 / 1.07 | 289 | 289 |
| cylinder_noisy | on | 8.038 | 1.171 | 0.2017 | 0.03893 | 0.000352 | 0.0002535 | 1.668e-06 | 0.000159 | 1e-06 | 0.10 / 0.69 | 0.07 / 0.19 | 416 | 416 |
| cylinder_noisy | off | 8.038 | 1.171 | 0.2017 | 0.03893 | 0.000352 | 0.0002535 | 1.668e-06 | 0.000159 | 1e-06 | 0.10 / 0.69 | 0.07 / 0.19 | 416 | 416 |
| icosphere2 | on | 4.547 | 2.084 | 3.583 | 1.454 | 2.07 | 1.395 | 1.554 | 1.011 | 1.49 | 0.40 / 0.85 | 0.35 / 0.66 | 160 | 160 |
| icosphere2 | off | 4.547 | 2.084 | 3.583 | 1.454 | 2.07 | 1.404 | 1.732 | 1.034 | 1.855 | 0.47 / 1.35 | 0.42 / 1.13 | 162 | 162 |

### 3.2 L∞ energies (their 1, 3, 5), control meshes

<details><summary>Tables</summary>

**Energy 1 of theirs = our `intrinsic-max`**

| mesh | remesh | E input (ours) | ref E@10 | ours E@10 | ref E@50 | ours E@50 | ref E@200 | ours E@200 | ref E@200 (their def.) | ours E@200 (their def.) | ref disp mean/max | ours disp mean/max | ref V | ours V |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| spot_ctrl | on | 179.3 | 159.9 | 97.38 | 40.11 | 46.81 | 22.23 | 35.15 | 22.23 | 35.15 | 0.28 / 0.78 | 0.24 / 0.68 | 184 | 183 |
| spot_ctrl | off | 179.3 | 144.2 | 152.7 | 37.96 | 60.83 | 27.49 | 23.97 | 27.49 | 23.97 | 0.33 / 1.02 | 0.34 / 1.03 | 188 | 188 |
| bob_ctrl | on | 58.49 | 33.7 | 36.13 | 18.56 | 19.17 | 11.51 | 18.24 | 11.51 | 18.24 | 0.22 / 1.24 | 0.13 / 0.94 | 331 | 330 |
| bob_ctrl | off | 58.49 | 33.7 | 36.13 | 15.87 | 19.17 | 15.75 | 17.38 | 15.75 | 17.38 | 0.17 / 1.13 | 0.14 / 0.99 | 334 | 334 |
| blub_ctrl | on | 76.03 | 30.81 | 32.78 | 9.855 | 13.04 | 5.148 | 6.361 | 5.148 | 6.361 | 0.20 / 0.47 | 0.17 / 0.45 | 110 | 108 |
| blub_ctrl | off | 76.03 | 30.81 | 32.78 | 16.87 | 10.17 | 8.076 | 6.222 | 8.076 | 6.222 | 0.20 / 0.61 | 0.22 / 0.51 | 112 | 112 |

**Energy 3 of theirs = our `covariance-max`**

| mesh | remesh | E input (ours) | ref E@10 | ours E@10 | ref E@50 | ours E@50 | ref E@200 | ours E@200 | ref E@200 (their def.) | ours E@200 (their def.) | ref disp mean/max | ours disp mean/max | ref V | ours V |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| spot_ctrl | on | 25.58 | 19.26 | 10.35 | 1.743 | 3.194 | 1.175 | 1.623 | 0.5876 | 0.8116 | 0.56 / 2.10 | 0.41 / 1.78 | 184 | 180 |
| spot_ctrl | off | 25.58 | 24.87 | 12.28 | 1.488 | 2.751 | 0.3458 | 1.704 | 0.1729 | 0.852 | 1.90 / 8.78 | 0.62 / 2.83 | 188 | 188 |
| bob_ctrl | on | 11.11 | 10.72 | 7.577 | 3.538 | 3.613 | 2.065 | 2.473 | 1.032 | 1.236 | 0.32 / 1.51 | 0.23 / 1.29 | 334 | 333 |
| bob_ctrl | off | 11.11 | 10.72 | 7.577 | 3.538 | 3.613 | 1.73 | 2.498 | 0.865 | 1.249 | 0.44 / 1.84 | 0.27 / 1.55 | 334 | 334 |
| blub_ctrl | on | 10.33 | 3.587 | 2.639 | 0.4394 | 0.7565 | degenerate face | 0.388 | 0 | 0.194 | 0.36 / 1.83 | 0.24 / 0.60 | 112 | 104 |
| blub_ctrl | off | 10.33 | 3.587 | 2.759 | 0.3211 | 0.785 | 0.1273 | 0.3875 | 0.06364 | 0.1938 | 0.62 / 3.14 | 0.39 / 0.89 | 112 | 112 |

**Energy 5 of theirs = our `combinatorial-max`**

| mesh | remesh | E input (ours) | ref E@10 | ours E@10 | ref E@50 | ours E@50 | ref E@200 | ours E@200 | ref E@200 (their def.) | ours E@200 (their def.) | ref disp mean/max | ours disp mean/max | ref V | ours V |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| spot_ctrl | on | 187.9 | 133.3 | 107.8 | 35.05 | 37.22 | 28.12 | 27.7 | 14.06 | 28.85 | 0.35 / 1.06 | 0.31 / 1.16 | 185 | 181 |
| spot_ctrl | off | 187.9 | 134.6 | 105.5 | 42.25 | 37.29 | 33.83 | 31.68 | 17.1 | 38.94 | 0.34 / 0.96 | 0.33 / 0.97 | 188 | 188 |
| bob_ctrl | on | 57.98 | 45.68 | 44.8 | 24.77 | 25.07 | 16.8 | 19.75 | 8.402 | 14.57 | 0.23 / 1.66 | 0.16 / 0.88 | 333 | 333 |
| bob_ctrl | off | 57.98 | 45.68 | 44.8 | 24.77 | 25.07 | 18.94 | 19.73 | 9.47 | 14.99 | 0.21 / 1.47 | 0.16 / 0.89 | 334 | 334 |
| blub_ctrl | on | 98.2 | 45.32 | 39.8 | 26.84 | 15.53 | 20.59 | 6.547 | 11.02 | 22.43 | 0.23 / 0.51 | 0.24 / 0.67 | 112 | 109 |
| blub_ctrl | off | 98.2 | 45.32 | 39.8 | 26.85 | 15.39 | 26.73 | 8.948 | 14.5 | 20.34 | 0.25 / 0.63 | 0.29 / 0.77 | 112 | 112 |

</details>

### 3.3 Suspected weakness in ours: slow L∞ covariance flow

On `covariance-max` (their energy 3) the reference is ahead from step 50 on in every finished
run. On spot_ctrl without remeshing, measured with our evaluator:

| steps | reference | ours |
|---|---|---|
| 50 | 1.49 | 2.75 |
| 200 | 0.346 | 1.70 |
| 1000 | (not run) | 1.21 |

Changing our options does not close the gap. Max displacement 0.1 or 2, initial step 0.01,
0.1 or 1, and L-BFGS memory 3–100 give 1.62–2.20 at 200 iterations.

What we know:

- The energy values agree with theirs up to the factor ½ (§2), and the factor cannot slow an
  L-BFGS run.
- So the difference is in the gradient/subgradient that is returned at the nonsmooth points,
  or in how our line search handles the max structure. Neither is identified.

**Follow-up (2026-10-20, two-team clean room).** A separate agent read the reference and returned
a prose-only list of differences (no code crossed; the clean side never saw its transcript). Its
findings: the energy, the candidate set and the gradient of this energy are in effect the same as
ours (full derivative through the circumcenter, no smoothing, Euclidean metric, same L-BFGS
two-loop and `H_0`); the optimizers differ in the first trial step (theirs: twice the previous
accepted step, up to 100, every iteration; ours: 1), in the constants (`c1 = 1e-4`, `c2 = 0.99`;
ours `0`, `0.5`), in the step cap (none), and on a failed line search (they keep the last trial
point even when it is higher, reset the memory and go on; ours retries once, then stops). They
have no face-angle floor.

Measured on our side (spot_ctrl, `covariance-max`, 200 iterations, no remeshing):

| change to ours | E@200 |
|---|---|
| none | 1.70 |
| `--min-face-angle 0` (no floor) | 1.70 (bitwise the same run) |
| `--wolfe 1e-4:0.99` | 1.83 |
| `--first-trial grow` (D37) | 1.10 |
| `--first-trial grow --wolfe 1e-4:0.99` | 1.12 |
| reference | 0.35 |

With `grow`, the 12 combinations initial step 0.003–1 × memory 8/10 give 1.05–1.82 (unit:
1.48–2.45), and 1000 iterations give 0.56 (unit: 1.21). Over six energies and the three control
meshes `grow` ends lower in 28 of 36 runs but traps spot `intrinsic` (SPEC D37), so it is an
option, not the default.

**The rest of the gap is the step cap.** A per-step black-box measurement of the reference (same
setup, every step synced and evaluated, by the same separate agent) found no failed line search,
no memory reset and no rise of the energy in 200 steps. Steps 1–14 barely move it (25.6 -> 24.3),
steps 15–24 take it to 8.6 with single-vertex moves of 0.6 to 3.5 mean edge lengths per step, then
it settles (1.49 at 50, 0.61 at 100, 0.35 at 200). The mesh inflates: at step 200 the mean edge
length is 3.1 × the input, the bounding-box diagonal 2.0 ×, the vertex displacement mean 1.90 and
max 8.78 input edge lengths, smallest angle 1.7e-4. Our cap (D19, 0.5 ℓ̄ per step) forbids those
steps. With `grow` and a larger cap (spot, 200 iterations, `E@200`):

| max displacement | initial step 0.001 | 0.01 | 0.1 |
|---|---|---|---|
| 0.5 (default) | 1.31 | 1.52 | 1.10 |
| 2 | 0.26 | **0.062** | 1.48 |
| 4 | 0.63 | 0.062 | 1.46 |
| none (1e6) | 0.015 | 0.66 | 0.92 |

The 0.062 run is as distorted as the reference's (edge length × 3.08, diagonal × 1.80,
displacement mean 2.02 / max 7.41, smallest angle 1.5e-3) at a sixth of its energy. Without a cap
one run reaches 0.015 by blowing the mesh up 375 × (the reason for D19). The results jump with the
options (0.062 to 1.48 at cap 2): the L∞ energy rewards a few large early moves, and whether the
flow finds them depends on the path. No error in our energy or gradient is implicated.

```
develop runs_tmp/dev-blackbox/meshes/spot_ctrl.obj --rounds covariance-max:200 --gradient-tol 0 \
  --stall-iters 100000000 --no-remesh --no-cut --first-trial grow [--wolfe 1e-4:0.99]
```

Reproducer (energy history in the file):

```
cargo run --release -p ddg-developability --example develop -- runs_tmp/dev-blackbox/meshes/spot_ctrl.obj \
  --rounds covariance-max:200 --gradient-tol 0 --stall-iters 100000000 --no-remesh --no-cut --history h.txt
# E 25.58 -> 1.70 after 200 iterations; the reference reaches 0.35
# (their displayed value x 2, or our evaluator on their out.obj)
```

The intrinsic L∞ energy (1) shows a milder gap of the same kind: the reference is lower in 4 of
6 runs at step 200.

## 4. Robustness (all 72 reference runs, 72 of ours)

| reference run (mesh, their energy, remeshing) | what happened |
|---|---|
| blub_ctrl, 2, on | NaN positions in `out.obj` by step 50 |
| blub_ctrl, 2, off | NaN positions by step 200 |
| bob_ctrl, 2, off | NaN positions by step 200 |
| spot_ctrl, 0, off | NaN positions by step 200 |
| spot_ctrl, 0, on | step 200: two identical faces (vertices 83, 84, 87) after its remeshing; then `f` did not return |
| blub_ctrl, 3, on | step 200: a zero-area face (our evaluator refuses it); `f` printed `nan` |
| blub_ctrl, 4, off | step 200: a face that violates the triangle inequality by rounding (our cut refuses it); `f` printed `nan` |
| blub_ctrl, 4, on | the process died between step 50 and 200, right after "A structural change happened to the mesh" |
| all 12 runs on crumpled and cylinder_noisy | `f` never returned (25 min; main thread at 100 % CPU, frame frozen, no files) |

- Their remeshing was on in 3 of the 8 failing runs.
- `f` printed `nan` as its error in 6 runs: the 4 with NaN positions and the 2 with a
  degenerate face.
- Ours: all 72 runs finished their 200 iterations with finite positions and a valid manifold
  mesh (`stopped=none` in every run). Our cut and flatten succeeded on all 72 results. It
  also succeeded on all 64 reference results that were valid meshes, including the 12
  boundary-mesh results where their `f` hangs. The other 7 have NaN positions (4), duplicate
  faces (1) or a degenerate face (2); the crashed run left no final mesh.

## 5. Time per step

Reference:

- Measured alone (load 0.9–3.7). Continuous optimization (Space) for 20 s with remeshing off,
  then "total steps" read from the panel (`scripts/timing_space.sh`). The process used about
  1.1 cores (`/proc` CPU time 21–25 s over 20 s on the dense meshes).
- Its GUI never does more than ~60 steps per second (1168–1207 steps in 20 s on the control
  meshes), so it does at most one step per frame. On the small meshes its time per step is the
  frame time.
- An attempt with TABs sent back to back lost most of the keys (4–5 of 100 were logged). It is
  not used.

Ours: the `develop` example, 100 iterations, remeshing off, `seconds / 100`. It includes the
line-search evaluations, 1.1–1.8 per iteration.

| mesh | faces | their energy / ours | reference s/step (GUI, 1 frame per step) | ours s/iteration, 1 thread | ours, 8 threads |
|---|---|---|---|---|---|
| spot_ctrl | 372 | 0 / intrinsic | 0.0168 (frame bound) | 0.00034 | 0.00047 |
| spot_ctrl | 372 | 2 / covariance | 0.0167 (frame bound) | 0.00031 | 0.00046 |
| spot_ctrl | 372 | 4 / combinatorial | 0.0170 (frame bound) | 0.00049 | 0.00064 |
| bob_ctrl | 668 | 0 | 0.0173 | 0.00061 | 0.00075 |
| bob_ctrl | 668 | 2 | 0.0168 | 0.00046 | 0.00055 |
| bob_ctrl | 668 | 4 | 0.0170 | 0.00080 | 0.00087 |
| spot_tri | 5856 | 0 | 0.068 | 0.0052 | 0.0030 |
| spot_tri | 5856 | 2 | 0.048 | 0.0042 | 0.0025 |
| spot_tri | 5856 | 4 | 0.040 | 0.0061 | 0.0030 |
| bob_tri | 10688 | 0 | 0.127 | 0.0081 | 0.0041 |
| bob_tri | 10688 | 2 | 0.074 | 0.0076 | 0.0041 |
| bob_tri | 10688 | 4 | 0.063 | 0.0101 | 0.0043 |
| blub_tri | 14208 | 0 | 0.159 | 0.0102 | 0.0053 |
| blub_tri | 14208 | 2 | 0.085 | 0.0100 | 0.0056 |
| blub_tri | 14208 | 4 | 0.081 | 0.0129 | 0.0056 |

The GUI numbers include whatever it does per step besides the optimization (rendering with
llvmpipe, the panel, the plot). They are upper bounds on its optimizer's cost, not a
measurement of it. With that caveat, on the dense meshes ours is 6–16× faster per iteration
on one thread.

## 6. Cut and flatten

All layouts were measured with our tools (`devbb uv`, `scripts/cutlen.py`).

Layouts:

- Their `f`: the file `mesh_with_flattened_texcoord.obj`.
- Ours on their final mesh: `develop --cut-only`, with the defect rule D30 (0.1 rad) or with
  the paper rule (`--defect-tol` unset).
- Ours on our own final mesh.

Quantities:

- `u` = per-face log scale factor ½ ln(A_uv/A); ū = its area-weighted mean. Both programs'
  layouts have arbitrary global scale.
- Q = quasi-conformal error.
- Flipped = triangles of minority orientation.
- Overlap = sampled self-overlap per chart.
- Cut length is a fraction of the total edge length.

| meshes | layout | layouts measured | flipped triangles (sum) | layouts with flips | median mean abs(u−ū) | median max abs(u−ū) | median max Q | median cut length |
|---|---|---|---|---|---|---|---|---|
| coarse (spot/bob/blub control) | reference `f` | 28 of 35 | 61 | 17 | 1.1 | 4.7 | 21.3 | 0.258 |
| coarse (spot/bob/blub control) | ours on their mesh, defect rule 0.1 | 28 of 35 | 0 | 0 | 4.8e-05 | 0.0026 | 1.01 | 0.200 |
| coarse (spot/bob/blub control) | ours on their mesh, paper rule | 28 of 35 | 352 | 16 | 0.068 | 1.4 | 47.2 | 0.191 |
| coarse (spot/bob/blub control) | ours on our mesh, defect rule 0.1 | 35 of 35 | 0 | 0 | 5e-05 | 0.0031 | 1.01 | 0.198 |
| dense (spot/bob/blub tri) | reference `f` | 18 of 18 | 1 | 1 | 0.24 | 1.5 | 1.17 | 0.187 |
| dense (spot/bob/blub tri) | ours on their mesh, defect rule 0.1 | 18 of 18 | 2 | 1 | 0.00072 | 0.022 | 1.05 | 0.141 |
| dense (spot/bob/blub tri) | ours on their mesh, paper rule | 18 of 18 | 2 | 2 | 0.00072 | 0.066 | 1.13 | 0.141 |
| dense (spot/bob/blub tri) | ours on our mesh, defect rule 0.1 | 18 of 18 | 1 | 1 | 0.00074 | 0.053 | 1.11 | 0.136 |
| boundary (crumpled, cylinder) | reference `f` | 0 of 12 (no output) | | | | | | |
| boundary (crumpled, cylinder) | ours on their mesh, defect rule 0.1 | 12 of 12 | 0 | 0 | 0.0004 | 0.0018 | 1 | 0.008 |
| boundary (crumpled, cylinder) | ours on their mesh, paper rule | 12 of 12 | 0 | 0 | 0.0004 | 0.0018 | 1 | 0.008 |
| boundary (crumpled, cylinder) | ours on our mesh, defect rule 0.1 | 12 of 12 | 0 | 0 | 0.00029 | 0.0018 | 1 | 0.008 |
| icosphere2 | reference `f` | 6 of 6 | 0 | 0 | 0.00017 | 0.0019 | 1 | 0.230 |
| icosphere2 | ours on their mesh, defect rule 0.1 | 6 of 6 | 0 | 0 | 0.00017 | 0.0023 | 1 | 0.182 |
| icosphere2 | ours on their mesh, paper rule | 6 of 6 | 0 | 0 | 0.00017 | 0.0023 | 1 | 0.182 |
| icosphere2 | ours on our mesh, defect rule 0.1 | 6 of 6 | 0 | 0 | 0.0007 | 0.004 | 1.01 | 0.136 |


Observations:

- **Area distortion.** Their layouts carry far more area distortion than ours, and their
  coarse-mesh layouts flip triangles (61 triangles in 17 of 28 layouts).
  - Their max Q on the dense meshes (median 1.17) is close to ours (1.05–1.13). Their
    printed error is Q.
  - The pattern fits a flattening that keeps angles but not lengths along the cut. Ours
    prescribes zero scale on the boundary (D23). That reading is inferred from the numbers
    only.
- **Cut length.** Their cuts are longer than ours: median 0.26 vs 0.20 of the total edge
  length on the coarse meshes, 0.19 vs 0.14 on the dense ones. Their terminal rule cannot be
  inferred from the files.
- **Paper rule vs defect rule on their meshes.** The paper rule (energy threshold only) gives
  shorter cuts (0.191). But it flips triangles in 16 of 28 coarse layouts and leaves median
  mean |u − ū| at 0.068. The defect rule D30 flips none, at mean |u − ū| ≈ 5e-5. This repeats
  the finding of SPEC D30 on a second program's output.
- **Meshes with boundary.** Their `f` produced nothing on these meshes (§4). Ours flattens all 12 of them without flips. The medians are a max |u − ū| of 0.0018
  and a cut of 0.8 % of the edge length. The worst case is our cut of their crumpled state
  (energy 0, remeshed), with max |u − ū| 0.34.

<details><summary>All 71 runs (cells: flipped / mean |u−ū| / max |u−ū| / mean Q / max Q / max chart overlap, then the cut length)</summary>

| run (mesh_energy_remesh) | reference `f` | reference printed error max / avg | ours on ref mesh, defect 0.1 | ours on ref mesh, paper rule | ours on our mesh, defect 0.1 |
|---|---|---|---|---|---|
| spot_ctrl_e0_def | — | no output (f did not return) | — | — | 0 / 7.93e-05 / 0.00769 / 1 / 1.02 / 0.25 (cut 0.213) |
| spot_ctrl_e0_n | err (cut nan) | -nan / -nan | — | — | 0 / 8.46e-05 / 0.00513 / 1 / 1.01 / 0.33 (cut 0.203) |
| spot_ctrl_e1_def | 2 / 1.23 / 4.71 / 1.273 / 11.6 / 0.14 (cut 0.270) | 11.5614 / 1.23982 | 0 / 3.91e-05 / 0.00134 / 1 / 1 / 0.44 (cut 0.203) | 6 / 0.114 / 2.6 / 2.376 / 305 / 0.43 (cut 0.200) | 0 / 6.61e-06 / 0.000349 / 1 / 1 / 0.27 (cut 0.211) |
| spot_ctrl_e1_n | 2 / 1.26 / 4.57 / 1.474 / 47.6 / 0.19 (cut 0.279) | 47.5882 / 1.36752 | 0 / 4.47e-06 / 0.000221 / 1 / 1 / 0.43 (cut 0.208) | 8 / 0.129 / 1.75 / 1.683 / 86.3 / 0.42 (cut 0.201) | 0 / 1.33e-05 / 0.0057 / 1 / 1.01 / 0.41 (cut 0.206) |
| spot_ctrl_e2_def | 3 / 1.35 / 4.75 / 1.288 / 31.2 / 0.17 (cut 0.226) | 31.1735 / 1.60409 | 0 / 0.000409 / 0.0125 / 1 / 1.03 / 0.32 (cut 0.179) | 60 / 0.459 / 2.15 / 5.137 / 152 / 0.42 (cut 0.134) | 0 / 4.92e-05 / 0.00606 / 1 / 1.01 / 0.44 (cut 0.178) |
| spot_ctrl_e2_n | 4 / 1.04 / 5.17 / 1.39 / 27.2 / 0.18 (cut 0.234) | 27.2445 / 1.31358 | 0 / 0.000182 / 0.00608 / 1 / 1.01 / 0.35 (cut 0.184) | 31 / 0.286 / 2.89 / 7.228 / 1.37e+03 / 0.39 (cut 0.160) | 0 / 0.000216 / 0.0224 / 1 / 1.05 / 0.20 (cut 0.185) |
| spot_ctrl_e3_def | 2 / 2.06 / 5.27 / 1.261 / 13.3 / 0.28 (cut 0.259) | 13.2969 / 1.23449 | 0 / 7.22e-05 / 0.00424 / 1 / 1.01 / 0.36 (cut 0.200) | 24 / 0.327 / 2.32 / 4.812 / 160 / 0.28 (cut 0.177) | 0 / 8.68e-05 / 0.012 / 1 / 1.02 / 0.34 (cut 0.198) |
| spot_ctrl_e3_n | 6 / 1.1 / 3.23 / 1.36 / 16.5 / 0.20 (cut 0.220) | 16.4819 / 1.38421 | 0 / 0.000136 / 0.059 / 1 / 1.13 / 0.24 (cut 0.186) | 47 / 0.403 / 3.18 / 15.05 / 4.4e+03 / 0.30 (cut 0.134) | 0 / 6.59e-05 / 0.00747 / 1 / 1.02 / 0.36 (cut 0.194) |
| spot_ctrl_e4_def | 2 / 1.92 / 6.47 / 1.915 / 57.2 / 0.10 (cut 0.274) | 57.184 / 1.38404 | 0 / 3.45e-05 / 0.02 / 1 / 1.04 / 0.44 (cut 0.206) | 0 / 3.45e-05 / 0.02 / 1 / 1.04 / 0.44 (cut 0.206) | 0 / 5.02e-05 / 0.0127 / 1 / 1.04 / 0.23 (cut 0.192) |
| spot_ctrl_e4_n | 4 / 2.68 / 6.07 / 2.105 / 74.7 / 0.11 (cut 0.268) | 74.7209 / 1.45185 | 0 / 6.12e-05 / 0.00157 / 1 / 1 / 0.30 (cut 0.207) | 0 / 0.0256 / 1.43 / 1.058 / 6.32 / 0.30 (cut 0.204) | 0 / 1.21e-05 / 0.00101 / 1 / 1 / 0.52 (cut 0.208) |
| spot_ctrl_e5_def | 3 / 2.17 / 4.4 / 2.552 / 45 / 0.07 (cut 0.257) | 44.9732 / 1.63393 | 0 / 1.48e-05 / 0.000687 / 1 / 1 / 0.29 (cut 0.207) | 0 / 1.48e-05 / 0.000687 / 1 / 1 / 0.29 (cut 0.207) | 0 / 8.81e-06 / 0.000355 / 1 / 1 / 0.38 (cut 0.212) |
| spot_ctrl_e5_n | 2 / 2.32 / 4.22 / 1.674 / 18.8 / 0.15 (cut 0.274) | 18.7645 / 1.30649 | 0 / 4.83e-05 / 0.00101 / 1 / 1 / 0.56 (cut 0.209) | 6 / 0.0525 / 1.45 / 1.229 / 19.6 / 0.54 (cut 0.206) | 0 / 1.75e-05 / 0.000346 / 1 / 1 / 0.33 (cut 0.210) |
| bob_ctrl_e0_def | 0 / 0.537 / 1.54 / 1.025 / 2.07 / 0.11 (cut 0.258) | 2.07105 / 1.02376 | 0 / 3.3e-05 / 0.000845 / 1 / 1 / 0.05 (cut 0.204) | 0 / 0.00521 / 0.678 / 1.016 / 3.9 / 0.05 (cut 0.204) | 0 / 4.98e-05 / 0.000868 / 1 / 1 / 0.05 (cut 0.186) |
| bob_ctrl_e0_n | 2 / 0.604 / 5.01 / 1.088 / 33.6 / 0.09 (cut 0.249) | 33.589 / 1.14913 | 0 / 4.52e-05 / 0.0408 / 1 / 1.09 / 0.04 (cut 0.189) | 5 / 0.0836 / 1.14 / 1.401 / 385 / 0.04 (cut 0.187) | 0 / 3.24e-05 / 0.00109 / 1 / 1 / 0.03 (cut 0.187) |
| bob_ctrl_e1_def | 2 / 0.624 / 5.82 / 1.079 / 19.5 / 0.44 (cut 0.252) | 19.5326 / 1.08705 | 0 / 7.62e-05 / 0.000963 / 1 / 1 / 0.03 (cut 0.215) | 23 / 0.14 / 2.07 / 2.954 / 476 / 0.06 (cut 0.213) | 0 / 6.89e-07 / 5.2e-05 / 1 / 1 / 0.14 (cut 0.233) |
| bob_ctrl_e1_n | 0 / 0.497 / 9.3 / 1.02 / 2.88e+04 / 0.00 (cut 0.265) | 28790.8 / 44.1164 | 0 / 1.69e-05 / 0.00141 / 1 / 1 / 0.06 (cut 0.227) | 0 / 1.69e-05 / 0.00141 / 1 / 1 / 0.06 (cut 0.227) | 0 / 1.41e-06 / 0.00149 / 1 / 1 / 0.15 (cut 0.233) |
| bob_ctrl_e2_def | 0 / 1.2 / 3.5 / 1.102 / 5.93 / 0.16 (cut 0.211) | 5.93027 / 1.08808 | 0 / 7.71e-05 / 0.00418 / 1 / 1.01 / 0.05 (cut 0.171) | 1 / 0.0402 / 1.29 / 1.151 / 12.6 / 0.07 (cut 0.171) | 0 / 0.00014 / 0.0126 / 1 / 1.03 / 0.06 (cut 0.175) |
| bob_ctrl_e2_n | err (cut nan) | nan / -nan | — | — | 0 / 8.44e-05 / 0.00205 / 1 / 1 / 0.07 (cut 0.175) |
| bob_ctrl_e3_def | 0 / 0.639 / 2.48 / 1.021 / 1.91 / 0.08 (cut 0.256) | 1.91258 / 1.03032 | 0 / 7.42e-05 / 0.000921 / 1 / 1 / 0.04 (cut 0.197) | 0 / 0.000861 / 0.116 / 1.002 / 1.24 / 0.04 (cut 0.195) | 0 / 0.000263 / 0.00746 / 1 / 1.02 / 0.14 (cut 0.203) |
| bob_ctrl_e3_n | 2 / 0.699 / 4.96 / 1.058 / 28.8 / 0.06 (cut 0.241) | 28.8247 / 1.10422 | 0 / 4.16e-05 / 0.00669 / 1 / 1.01 / 0.10 (cut 0.187) | 14 / 0.0998 / 3.01 / 2.892 / 431 / 0.11 (cut 0.184) | 0 / 3.02e-05 / 0.000852 / 1 / 1 / 0.09 (cut 0.208) |
| bob_ctrl_e4_def | 0 / 0.774 / 3.71 / 1.054 / 7.57 / 0.18 (cut 0.206) | 7.57207 / 1.09753 | 0 / 7.48e-05 / 0.00253 / 1 / 1.01 / 0.07 (cut 0.166) | 0 / 7.48e-05 / 0.00253 / 1 / 1.01 / 0.07 (cut 0.166) | 0 / 0.000175 / 0.00275 / 1 / 1 / 0.09 (cut 0.189) |
| bob_ctrl_e4_n | 0 / 0.788 / 11.2 / 1.049 / 2.29e+04 / 0.00 (cut 0.218) | 22908.3 / 35.3527 | 0 / 6.96e-05 / 0.00275 / 1 / 1.01 / 0.09 (cut 0.175) | 0 / 6.96e-05 / 0.00275 / 1 / 1.01 / 0.09 (cut 0.175) | 0 / 6.45e-05 / 0.000978 / 1 / 1 / 0.10 (cut 0.184) |
| bob_ctrl_e5_def | 0 / 0.401 / 4.24 / 1.047 / 8.91 / 0.12 (cut 0.265) | 8.90594 / 1.08193 | 0 / 4.78e-05 / 0.00126 / 1 / 1 / 0.06 (cut 0.229) | 0 / 4.78e-05 / 0.00126 / 1 / 1 / 0.06 (cut 0.229) | 0 / 9.19e-06 / 0.00227 / 1 / 1.01 / 0.11 (cut 0.232) |
| bob_ctrl_e5_n | 0 / 0.7 / 2.74 / 1.037 / 2.25 / 0.05 (cut 0.265) | 2.25433 / 1.0525 | 0 / 4.69e-05 / 0.00159 / 1 / 1 / 0.04 (cut 0.222) | 0 / 4.69e-05 / 0.00159 / 1 / 1 / 0.04 (cut 0.222) | 0 / 1.42e-05 / 0.00313 / 1 / 1.01 / 0.12 (cut 0.232) |
| blub_ctrl_e0_def | 0 / 1.2 / 5.36 / 1.322 / 23.2 / 0.27 (cut 0.275) | 23.1531 / 1.36936 | 0 / 0.000139 / 0.00344 / 1 / 1.01 / 0.29 (cut 0.208) | 35 / 0.339 / 2.46 / 5.801 / 1.4e+03 / 0.38 (cut 0.186) | 0 / 0.000163 / 0.00523 / 1 / 1.01 / 0.34 (cut 0.202) |
| blub_ctrl_e0_n | 2 / 1.71 / 6 / 1.932 / 17.9 / 0.07 (cut 0.271) | 17.913 / 1.68005 | 0 / 1.69e-05 / 0.000656 / 1 / 1 / 0.19 (cut 0.199) | 22 / 0.323 / 2.88 / 5.376 / 289 / 0.25 (cut 0.177) | 0 / 3.58e-05 / 0.00118 / 1 / 1 / 0.19 (cut 0.203) |
| blub_ctrl_e1_def | 21 / 2.58 / 6.64 / 10.45 / 1.17e+03 / 0.06 (cut 0.205) | 1172.47 / 15.8631 | 0 / 2.79e-05 / 0.00282 / 1 / 1.01 / 0.28 (cut 0.195) | 34 / 1.47 / 5.19 / 55.92 / 6.52e+03 / 0.03 (cut 0.137) | 0 / 0.000276 / 0.0386 / 1 / 1.08 / 0.15 (cut 0.196) |
| blub_ctrl_e1_n | 1 / 1.68 / 6.71 / 2.693 / 154 / 0.27 (cut 0.255) | 154.483 / 2.01548 | 0 / 1.26e-05 / 0.00334 / 1 / 1.01 / 0.17 (cut 0.198) | 6 / 0.126 / 2.01 / 1.852 / 74.7 / 0.12 (cut 0.185) | 0 / 3.01e-05 / 0.00135 / 1 / 1 / 0.11 (cut 0.195) |
| blub_ctrl_e2_def | err (cut nan) | -nan / nan | — | — | 0 / 0.000708 / 0.0195 / 1.001 / 1.04 / 0.22 (cut 0.173) |
| blub_ctrl_e2_n | err (cut nan) | -nan / nan | — | — | 0 / 0.00121 / 0.166 / 1.002 / 1.4 / 0.21 (cut 0.163) |
| blub_ctrl_e3_def | err (cut 0.252) | -nan / -nan | — | — | 0 / 0.000247 / 0.0344 / 1 / 1.07 / 0.13 (cut 0.196) |
| blub_ctrl_e3_n | 1 / 1.13 / 3.15 / 1.317 / 41.4 / 0.28 (cut 0.223) | 41.3802 / 1.39857 | 0 / 0.0015 / 0.0388 / 1.003 / 1.07 / 0.24 (cut 0.177) | 30 / 0.757 / 3.38 / 17.97 / 427 / 0.17 (cut 0.100) | 0 / 0.000254 / 0.324 / 1.001 / 1.91 / 0.32 (cut 0.182) |
| blub_ctrl_e4_n | err (cut 0.237) | nan / nan | — | — | 0 / 0.000291 / 0.0119 / 1 / 1.02 / 0.28 (cut 0.187) |
| blub_ctrl_e5_def | 0 / 0.209 / 2.42 / 1.019 / 2.46 / 0.09 (cut 0.266) | 2.45665 / 1.03866 | 0 / 0.000344 / 0.0142 / 1 / 1.03 / 0.25 (cut 0.201) | 0 / 0.0239 / 0.86 / 1.059 / 7.22 / 0.24 (cut 0.198) | 0 / 9.01e-06 / 0.000158 / 1 / 1 / 0.14 (cut 0.210) |
| blub_ctrl_e5_n | 0 / 8.9e-05 / 0.00179 / 1 / 1 / 0.40 (cut 0.267) | 1.00347 / 1.00012 | 0 / 2.83e-05 / 0.00111 / 1 / 1 / 0.12 (cut 0.206) | 0 / 2.83e-05 / 0.00111 / 1 / 1 / 0.12 (cut 0.206) | 0 / 1.18e-05 / 0.000596 / 1 / 1 / 0.17 (cut 0.200) |
| spot_tri_e0_def | 0 / 0.114 / 0.354 / 1 / 1.22 / 0.54 (cut 0.237) | 1.22269 / 1.00048 | 0 / 0.000235 / 0.0312 / 1 / 1.06 / 0.07 (cut 0.171) | 0 / 0.000243 / 0.104 / 1 / 1.27 / 0.07 (cut 0.171) | 0 / 0.000269 / 0.0143 / 1 / 1.04 / 0.13 (cut 0.168) |
| spot_tri_e0_n | 0 / 0.202 / 1.02 / 1.001 / 1.4 / 0.61 (cut 0.243) | 1.39995 / 1.00092 | 0 / 0.000197 / 0.0149 / 1 / 1.04 / 0.16 (cut 0.172) | 0 / 0.000216 / 0.226 / 1 / 1.6 / 0.16 (cut 0.172) | 0 / 0.000279 / 0.214 / 1 / 1.54 / 0.19 (cut 0.166) |
| spot_tri_e2_def | 0 / 0.0551 / 0.15 / 1 / 1.07 / 0.49 (cut 0.229) | 1.06939 / 1.00042 | 0 / 0.000325 / 0.00846 / 1 / 1.03 / 0.14 (cut 0.157) | 0 / 0.000327 / 0.0644 / 1 / 1.16 / 0.14 (cut 0.157) | 0 / 0.000367 / 0.071 / 1 / 1.17 / 0.20 (cut 0.153) |
| spot_tri_e2_n | 0 / 0.233 / 2.27 / 1.001 / 2.67 / 0.15 (cut 0.206) | 2.67152 / 1.00412 | 0 / 0.000341 / 0.31 / 1 / 1.86 / 0.13 (cut 0.139) | 1 / 0.000657 / 0.751 / 1.001 / 4.45 / 0.13 (cut 0.138) | 0 / 0.000481 / 0.0508 / 1 / 1.11 / 0.28 (cut 0.148) |
| spot_tri_e4_def | 0 / 0.0375 / 0.0973 / 1 / 1 / 0.41 (cut 0.207) | 1.0048 / 1.00033 | 0 / 0.000454 / 0.00799 / 1 / 1.01 / 0.13 (cut 0.155) | 0 / 0.000454 / 0.00799 / 1 / 1.01 / 0.13 (cut 0.155) | 0 / 0.000244 / 0.00822 / 1 / 1.02 / 0.22 (cut 0.166) |
| spot_tri_e4_n | 0 / 0.0375 / 0.0973 / 1 / 1 / 0.41 (cut 0.207) | 1.0048 / 1.00033 | 0 / 0.000454 / 0.00799 / 1 / 1.01 / 0.13 (cut 0.155) | 0 / 0.000454 / 0.00799 / 1 / 1.01 / 0.13 (cut 0.155) | 0 / 0.00025 / 0.164 / 1 / 1.38 / 0.15 (cut 0.164) |
| bob_tri_e0_def | 0 / 0.219 / 1.75 / 1.001 / 1.09 / 0.12 (cut 0.175) | 1.09437 / 1.00177 | 0 / 0.000718 / 0.013 / 1 / 1.03 / 0.19 (cut 0.144) | 0 / 0.000718 / 0.013 / 1 / 1.03 / 0.19 (cut 0.144) | 0 / 0.00079 / 0.0454 / 1.001 / 1.09 / 0.15 (cut 0.134) |
| bob_tri_e0_n | 1 / 1.14 / 5.96 / 1.008 / 9.2 / 0.05 (cut 0.172) | 9.19905 / 1.00798 | 0 / 0.000718 / 0.248 / 1.001 / 1.66 / 0.13 (cut 0.141) | 0 / 0.00072 / 0.121 / 1.001 / 1.25 / 0.13 (cut 0.141) | 0 / 0.000739 / 0.0547 / 1.001 / 1.12 / 0.10 (cut 0.131) |
| bob_tri_e2_def | 0 / 0.184 / 0.697 / 1.001 / 1.13 / 0.24 (cut 0.174) | 1.12547 / 1.00143 | 0 / 0.000776 / 0.275 / 1.001 / 1.73 / 0.11 (cut 0.142) | 0 / 0.000792 / 0.296 / 1.001 / 1.86 / 0.11 (cut 0.142) | 0 / 0.00082 / 0.00937 / 1.001 / 1.01 / 0.14 (cut 0.137) |
| bob_tri_e2_n | 0 / 0.245 / 2.55 / 1.002 / 1.59 / 0.07 (cut 0.169) | 1.58817 / 1.00339 | 0 / 0.000718 / 0.0257 / 1.001 / 1.06 / 0.10 (cut 0.135) | 0 / 0.000718 / 0.283 / 1.001 / 1.79 / 0.10 (cut 0.135) | 0 / 0.000743 / 0.0955 / 1.001 / 1.21 / 0.07 (cut 0.127) |
| bob_tri_e4_def | 0 / 0.823 / 2.36 / 1.002 / 1.08 / 0.14 (cut 0.188) | 1.08448 / 1.00246 | 0 / 0.000703 / 0.00661 / 1 / 1.01 / 0.06 (cut 0.152) | 0 / 0.000703 / 0.00661 / 1 / 1.01 / 0.06 (cut 0.152) | 0 / 0.000735 / 0.0094 / 1 / 1.01 / 0.14 (cut 0.149) |
| bob_tri_e4_n | 0 / 0.823 / 2.36 / 1.002 / 1.08 / 0.14 (cut 0.188) | 1.08448 / 1.00246 | 0 / 0.000703 / 0.00661 / 1 / 1.01 / 0.06 (cut 0.152) | 0 / 0.000703 / 0.00661 / 1 / 1.01 / 0.06 (cut 0.152) | 0 / 0.000735 / 0.0094 / 1 / 1.01 / 0.14 (cut 0.149) |
| blub_tri_e0_def | 0 / 0.213 / 0.651 / 1.001 / 1.06 / 0.13 (cut 0.145) | 1.06168 / 1.00126 | 0 / 0.00127 / 0.0604 / 1.001 / 1.1 / 0.12 (cut 0.100) | 0 / 0.00127 / 0.0613 / 1.001 / 1.1 / 0.12 (cut 0.100) | 0 / 0.00122 / 0.0794 / 1.001 / 1.17 / 0.22 (cut 0.098) |
| blub_tri_e0_n | 0 / 0.457 / 1.25 / 1.001 / 9.97 / 0.12 (cut 0.144) | 9.97243 / 1.0019 | 0 / 0.00143 / 0.419 / 1.001 / 2.32 / 0.11 (cut 0.100) | 0 / 0.00143 / 0.419 / 1.001 / 2.32 / 0.11 (cut 0.100) | 0 / 0.00139 / 0.29 / 1.001 / 1.89 / 0.15 (cut 0.098) |
| blub_tri_e2_def | 0 / 0.478 / 1.2 / 1.001 / 1.07 / 0.15 (cut 0.144) | 1.06548 / 1.00112 | 0 / 0.00136 / 0.0667 / 1.001 / 1.07 / 0.23 (cut 0.099) | 0 / 0.00139 / 0.0668 / 1.001 / 1.07 / 0.23 (cut 0.099) | 0 / 0.00114 / 0.0668 / 1.001 / 1.14 / 0.21 (cut 0.099) |
| blub_tri_e2_n | 0 / 0.478 / 3.45 / 1.001 / 1.52 / 0.16 (cut 0.124) | 1.5175 / 1.00246 | 2 / 0.00167 / 1.51 / 1.001 / 20.7 / 0.19 (cut 0.084) | 1 / 0.00204 / 1.53 / 1.002 / 21.5 / 0.18 (cut 0.084) | 1 / 0.00159 / 1.44 / 1.001 / 17.9 / 0.25 (cut 0.087) |
| blub_tri_e4_def | 0 / 0.362 / 2.06 / 1.001 / 1.3 / 0.13 (cut 0.187) | 1.29526 / 1.0012 | 0 / 0.000965 / 0.019 / 1.001 / 1.02 / 0.22 (cut 0.123) | 0 / 0.000965 / 0.019 / 1.001 / 1.02 / 0.22 (cut 0.123) | 0 / 0.000964 / 0.0294 / 1.001 / 1.03 / 0.28 (cut 0.110) |
| blub_tri_e4_n | 0 / 0.362 / 2.06 / 1.001 / 1.3 / 0.13 (cut 0.187) | 1.29526 / 1.0012 | 0 / 0.000965 / 0.019 / 1.001 / 1.02 / 0.22 (cut 0.123) | 0 / 0.000965 / 0.019 / 1.001 / 1.02 / 0.22 (cut 0.123) | 0 / 0.000964 / 0.0294 / 1.001 / 1.03 / 0.28 (cut 0.110) |
| crumpled_e0_def | — | no output (f did not return) | 0 / 0.00179 / 0.34 / 1.001 / 1.97 / 0.00 (cut 0.000) | 0 / 0.00179 / 0.34 / 1.001 / 1.97 / 0.00 (cut 0.000) | 0 / 0.000548 / 0.00404 / 1 / 1 / 0.00 (cut 0.000) |
| crumpled_e0_n | — | no output (f did not return) | 0 / 0.00142 / 0.00603 / 1.001 / 1.01 / 0.00 (cut 0.001) | 0 / 0.0062 / 0.0822 / 1.004 / 1.12 / 0.00 (cut 0.000) | 0 / 0.000899 / 0.0553 / 1.001 / 1.11 / 0.00 (cut 0.000) |
| crumpled_e2_def | — | no output (f did not return) | 0 / 0.00255 / 0.0371 / 1.003 / 1.09 / 0.01 (cut 0.018) | 0 / 0.0027 / 0.0574 / 1.003 / 1.08 / 0.00 (cut 0.018) | 0 / 0.000943 / 0.00347 / 1 / 1.01 / 0.00 (cut 0.003) |
| crumpled_e2_n | — | no output (f did not return) | 0 / 0.00153 / 0.0507 / 1.001 / 1.1 / 0.01 (cut 0.017) | 0 / 0.0112 / 0.806 / 1.02 / 4.72 / 0.00 (cut 0.000) | 0 / 0.000809 / 0.12 / 1.001 / 1.27 / 0.01 (cut 0.015) |
| crumpled_e4_def | — | no output (f did not return) | 0 / 0.000535 / 0.0024 / 1 / 1 / 0.00 (cut 0.000) | 0 / 0.000535 / 0.0024 / 1 / 1 / 0.00 (cut 0.000) | 0 / 0.000648 / 0.00435 / 1.001 / 1.01 / 0.00 (cut 0.016) |
| crumpled_e4_n | — | no output (f did not return) | 0 / 0.000535 / 0.0024 / 1 / 1 / 0.00 (cut 0.000) | 0 / 0.000535 / 0.0024 / 1 / 1 / 0.00 (cut 0.000) | 0 / 0.000648 / 0.00435 / 1.001 / 1.01 / 0.00 (cut 0.016) |
| cylinder_noisy_e0_def | — | no output (f did not return) | 0 / 3.31e-05 / 0.00011 / 1 / 1 / 0.00 (cut 0.008) | 0 / 3.31e-05 / 0.00011 / 1 / 1 / 0.00 (cut 0.008) | 0 / 3.32e-05 / 0.000133 / 1 / 1 / 0.00 (cut 0.008) |
| cylinder_noisy_e0_n | — | no output (f did not return) | 0 / 3.31e-05 / 0.00011 / 1 / 1 / 0.00 (cut 0.008) | 0 / 3.31e-05 / 0.00011 / 1 / 1 / 0.00 (cut 0.008) | 0 / 3.32e-05 / 0.000133 / 1 / 1 / 0.00 (cut 0.008) |
| cylinder_noisy_e2_def | — | no output (f did not return) | 0 / 2.53e-05 / 0.000104 / 1 / 1 / 0.00 (cut 0.008) | 0 / 2.53e-05 / 0.000104 / 1 / 1 / 0.00 (cut 0.008) | 0 / 3.14e-05 / 0.000125 / 1 / 1 / 0.00 (cut 0.008) |
| cylinder_noisy_e2_n | — | no output (f did not return) | 0 / 2.53e-05 / 0.000104 / 1 / 1 / 0.00 (cut 0.008) | 0 / 2.53e-05 / 0.000104 / 1 / 1 / 0.00 (cut 0.008) | 0 / 3.14e-05 / 0.000125 / 1 / 1 / 0.00 (cut 0.008) |
| cylinder_noisy_e4_def | — | no output (f did not return) | 0 / 0.00027 / 0.00128 / 1 / 1 / 0.00 (cut 0.008) | 0 / 0.00027 / 0.00128 / 1 / 1 / 0.00 (cut 0.008) | 0 / 1.99e-05 / 7.17e-05 / 1 / 1 / 0.00 (cut 0.008) |
| cylinder_noisy_e4_n | — | no output (f did not return) | 0 / 0.00027 / 0.00128 / 1 / 1 / 0.00 (cut 0.008) | 0 / 0.00027 / 0.00128 / 1 / 1 / 0.00 (cut 0.008) | 0 / 1.99e-05 / 7.17e-05 / 1 / 1 / 0.00 (cut 0.008) |
| icosphere2_e0_def | 0 / 7.16e-05 / 0.000687 / 1 / 1 / 0.00 (cut 0.248) | 1.0011 / 1.00015 | 0 / 0.0001 / 0.000874 / 1 / 1 / 0.00 (cut 0.208) | 0 / 0.0001 / 0.000874 / 1 / 1 / 0.00 (cut 0.208) | 0 / 0.00177 / 0.00392 / 1.003 / 1.01 / 0.00 (cut 0.089) |
| icosphere2_e0_n | 0 / 7.16e-05 / 0.000687 / 1 / 1 / 0.00 (cut 0.248) | 1.0011 / 1.00015 | 0 / 0.0001 / 0.000874 / 1 / 1 / 0.00 (cut 0.208) | 0 / 0.0001 / 0.000874 / 1 / 1 / 0.00 (cut 0.208) | 0 / 0.0011 / 0.621 / 1.002 / 3.46 / 0.00 (cut 0.094) |
| icosphere2_e2_def | 0 / 0.000211 / 0.00159 / 1 / 1 / 0.01 (cut 0.221) | 1.00307 / 1.00016 | 0 / 0.000164 / 0.00271 / 1 / 1.01 / 0.00 (cut 0.166) | 0 / 0.000164 / 0.00271 / 1 / 1.01 / 0.00 (cut 0.166) | 0 / 0.000953 / 0.0041 / 1.001 / 1.01 / 0.00 (cut 0.120) |
| icosphere2_e2_n | 0 / 0.000124 / 0.0392 / 1 / 1.59 / 0.01 (cut 0.238) | 1.59383 / 1.00201 | 0 / 0.000173 / 0.0138 / 1 / 1.03 / 0.01 (cut 0.178) | 0 / 0.000173 / 0.0138 / 1 / 1.03 / 0.01 (cut 0.178) | 0 / 0.000444 / 0.0106 / 1.001 / 1.02 / 0.00 (cut 0.152) |
| icosphere2_e4_def | 0 / 0.00113 / 0.00722 / 1 / 1.01 / 0.00 (cut 0.221) | 1.01388 / 1.00051 | 0 / 0.000317 / 0.0118 / 1.001 / 1.02 / 0.00 (cut 0.185) | 0 / 0.000317 / 0.0118 / 1.001 / 1.02 / 0.00 (cut 0.185) | 0 / 0.000195 / 0.00338 / 1 / 1.01 / 0.00 (cut 0.185) |
| icosphere2_e4_n | 0 / 0.00038 / 0.00228 / 1 / 1 / 0.00 (cut 0.201) | 1.00301 / 1.00034 | 0 / 0.000224 / 0.00186 / 1 / 1 / 0.00 (cut 0.180) | 0 / 0.000224 / 0.00186 / 1 / 1 / 0.00 (cut 0.180) | 0 / 9.5e-05 / 0.00264 / 1 / 1.01 / 0.00 (cut 0.185) |

</details>

## 7. Pictures

- `runs_tmp/dev-blackbox/an/sheet.png`, a contact sheet of six cases:
  - spot_ctrl E_λ, bob_ctrl intrinsic, blub_ctrl E_P, spot_tri E_λ, blub_tri intrinsic,
    crumpled E_P, all without remeshing;
  - columns: input, reference after 200 steps, ours after 200 iterations, their `f` layout,
    our cut + BFF layout;
  - layouts are coloured by |u − ū| on a common 0…1 scale, with flipped triangles green.
- Screenshots of every reference state used in §2 are under `ref/probe_e/`, `ref/probe_cyl/`,
  `ref/fan/`, `ref/fan2/` and `ref/main/*/NNN_Tabx*.png` (text panel). `an/panels_probe_e.png`
  stacks the spot panels for all six energies.

## 8. Not compared, and why

- Remeshing details:
  - Their thresholds and operations cannot be observed beyond stdout's counts of
    collapses/flips.
  - Both programs run with remeshing on (their default, our default `RemeshOptions`) and off.
    The vertex counts differ (spot_ctrl E_λ: theirs 183, ours 171).
- Gradient descent with backtracking and the fixed step (their `m` ×2, ×3) were identified but
  not run in the matrix. Only L-BFGS was compared.
- 4-1 subdivision (`s`), multi-round coarse-to-fine schedules (our D36 default flow), runs
  longer than 200 steps, and the `p` plot / `c`, `d` colourings were not exercised.
- Their "Hausdorff dist to original" was not compared; we report per-vertex or nearest-vertex
  displacement instead.
- Energy before the first step: the panel shows no value until a step is taken, so every
  identical-geometry check uses the state after a step (the file written by `w`).
- Timing:
  - Not apples to apples. Their numbers come from the GUI loop, which does at most one step
    per frame (≈ 60 Hz).
  - Back-to-back TABs are dropped, so a per-key measurement was impossible.
  - The reference's thread use was not controlled beyond `OMP_NUM_THREADS=1` and
    `LP_NUM_THREADS=1`; it used about 1.1 cores.
- Their `f` on meshes with boundary: never produced output (§4), so there is no comparison on
  those meshes.
- Overlap between charts and packing: neither program packs; overlap is per chart only.

## 9. Findings about ours

1. **L∞ covariance flow progresses slowly** compared with the reference (§3.3, reproducer
   there). The energy definition agrees up to the factor ½. Partly explained (§3.3 follow-up): the
   first trial step of the line search; `FirstTrial::GrowPrevious` (SPEC D37) brings 1.70 to 1.10
   at 200 iterations against their 0.35. The rest is our step cap (D19): with `grow` and a cap of
   2 ℓ̄ ours reaches 0.062 at the same distortion as theirs. Not a correctness error.
2. **`examples/develop.rs` reports the wrong energy name in `RESULT energy=`** when `--rounds`
   is given. It prints `--energy` (default `Covariance`) instead of the energy of the rounds.
   Repro: `develop spot.obj --rounds intrinsic:10 --no-cut | tail -1` shows
   `energy=Covariance`. Cosmetic; the numbers are the intrinsic ones.
3. Definition differences from the reference, all observed (§2), none an error on our side.
   Our choices are documented in SPEC D1/D8:
   - boundary vertices are counted by us and not by them;
   - E_P is normalized per arc by us and pooled over within-arc pairs by them, and they forbid
     single-face arcs in E_P and π_max;
   - their L∞ covariance and L∞ combinatorial energies are half of ours.
4. Robustness: on the same inputs and step budgets ours never produced NaN, a degenerate or
   duplicate face, a crash or a hang. The reference did so in 8 of 72 flows, and its `f` hung
   in 13 of 72 runs: the 12 boundary-mesh runs and the duplicate-face state.

## 10. Commands

All paths are relative to `runs_tmp/dev-blackbox/`. Our binaries live in
`CARGO_TARGET_DIR=$PWD/target-dev-bb` (repository root).

```
# reference build (out of tree; the source dir is chmod 755 only inside these scripts, then 000)
scripts/configure.sh -DCMAKE_CXX_FLAGS="-include /var/data/ddg-quarantine/prefix/shim/glew_shim.h" \
                     -DCMAKE_CXX_STANDARD_LIBRARIES="-lGLU"
scripts/build.sh

# test meshes
python3 -I scripts/meshes.py; python3 -I scripts/fans.py; python3 -I scripts/fans2.py

# our tools
cargo build --release -p ddg-developability --example develop
(cd driver && cargo build --release)          # devbb: eval / disp / uv

# energy identification (one display per job)
ZOOM=-80 scripts/probe.sh 92 meshes/spot_ctrl.obj ref/probe_e m Tab:so e Tab:so e Tab:so e Tab:so e Tab:so e Tab:so e Tab:so
devbb eval ref/probe_e/002_Tab.obj an/x                         # then sum interior, valence != 3
uv run --with numpy python -I scripts/obs_defs.py ref/probe_e/010_Tab.obj   # their E4 / E5

# flow matrix (72 reference jobs, 16 in parallel; 216 of ours)
xargs -P 16 -L 1 scripts/ref_job.sh  < an/ref_jobs.txt
xargs -P 10 -L 1 scripts/ours_job.sh < an/ours_jobs.txt
ls ref/main | sed 's#^#ref/main/#' | xargs -P 16 -L 1 scripts/cut_jobs.sh
uv run python -I scripts/analyze.py     # -> an/results.csv
python3 -I scripts/tables.py flow 0,2,4; python3 -I scripts/tables.py flow 1,3,5; python3 -I scripts/cut_table.py

# timing (alone)
scripts/timing_space.sh 251 MESH EIDX 20
develop meshes/MESH.obj --rounds ENERGY:100 --gradient-tol 0 --stall-iters 100000000 --threads 1 --no-remesh --no-cut

# pictures
uv run scripts/sheet.py an/sheet.png "spot_ctrl 2 30 -70" "bob_ctrl 0 20 -70" "blub_ctrl 4 30 -70" \
    "spot_tri 2 30 -70" "blub_tri 0 30 -70" "crumpled 4 0 -40"
```
