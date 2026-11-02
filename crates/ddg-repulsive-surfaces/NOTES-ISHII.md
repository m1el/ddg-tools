# NOTES-ISHII: Fig. 16 (genus-two handlebody-knots of Ishii, Kishimoto, Moriuchi, Suzuki 2012)

Code: `src/generate/ishii.rs` (+ generated data `src/generate/ishii_data.rs`), `tools/ishii_extract.py` (the tracer),
`examples/ishii.rs` (`list`, `check`, `run`), `tests/ishii.rs` (run with `--release`, 11 s).
Source: `/var/data/discrete-geometry/ishii2012.pdf` (not copied into the repository), page 3 = Table 1.
Build: `CARGO_TARGET_DIR=$PWD/target-ishii cargo build --release -p ddg-repulsive-surfaces --example ishii`.

## Transcription

All 21 entries that Fig. 16 shows are transcribed: 4_1, 5_1..5_4, 6_1..6_16; none skipped. (The paper's Table 1 has
no entries with fewer than 4 crossings and the trivial 0_1 is not in Fig. 16.)

Method (nothing typed from memory). Page 3 is rendered at 500 dpi; each diagram is cut out, closed over 1 px render
seams, skeletonised, split at the two trivalent vertices into strand pieces; the free ends of the pieces are the
ends of under-crossing gaps and are paired by distance (< 90 px), alignment of the tangents (cos > 0.6) and an
over-strand between them. Edges are followed from vertex to vertex through the pairs. The 2D curves are resampled
and smoothed lightly (ends fixed), scaled to a diagram width of about 4. Each crossing is found from the traced 2D
curves; the under-strand is the one that passes through a gap. Lift: `z = +0.22` on the over-strand and `-0.22` on
the under-strand at each crossing, cosine ramps between consecutive crossings, `z = 0` at the vertices.
Crossing separation 0.44; the smallest 3D clearance of the graph is 0.30 (6_1), the tube radius is 0.08
(`min(0.3 clearance, 0.08)`), so a tube diameter is 0.16 and the over/under offset is 2.75 diameters.

Checks (all pass; `ishii check`, `tests/ishii.rs`):
1. the number of gap pairs found by the tracer equals the crossing number in the entry's name (4, 5 or 6), with no
   unpaired gap end and no undecided crossing (this is an independent check: the number is not an input);
2. the overlay of the traced curves on the diagram was inspected for all 21 (gap segments drawn thick, continuity of
   the strands through each gap; `runs_tmp/ishii/out/*_ov.png`, not committed);
3. `check_projection`: the crossings of the stored 3D polylines in the xy plane (recomputed in Rust from the
   decimated points) are as many as the entry's crossing number, and the over/under sequence along each edge
   (with the index of the crossed edge) equals the one recorded from the diagram. A mirror image (z -> -z) fails it
   (test `mirror_image_fails_the_diagram_check`; 6_12 and 6_13 have the same projection and differ in signs).

Limits. The lifted graph is the diagram's spatial graph (planar diagram plus crossing information); the paper
shows irreducible handlebody-knots up to equivalence, our graphs are those drawings. Vertex positions, edge shapes
are as drawn. I did not compute an invariant to confirm that two entries are inequivalent. The classification of a
vertex as the junction rests on skeleton junctions; each diagram yields exactly two.

Edges: `e<k> v<a>-><b>`: from vertex a to b; `O<j>` / `U<j>`: the edge is over / under at a crossing with edge j (j = k
for a self-crossing), in order along the edge. v0, v1 are the two vertices in the order of the tracer.

| entry | crossings | graph | crossings along each edge |
|---|---|---|---|
| 4_1 | 4 | handcuff graph | e0 v1->v1: O2 U2<br>e1 v0->v0: O2 U2<br>e2 v1->v0: O1 U1 O0 U0 |
| 5_1 | 5 | handcuff graph | e0 v1->v0: O2 U2 O2<br>e1 v1->v1: U2 O2<br>e2 v0->v0: U0 O0 U0 O1 U1 |
| 5_2 | 5 | theta curve | e0 v1->v0: U1 O1<br>e1 v1->v0: O1 U0 O2 U2 O0 U1<br>e2 v0->v1: U1 O1 |
| 5_3 | 5 | theta curve | e0 v1->v0: O1 U2<br>e1 v0->v1: O2 U0 O2 U2<br>e2 v0->v1: O1 U1 O0 U1 |
| 5_4 | 5 | handcuff graph | e0 v1->v1: O1 U1<br>e1 v0->v0: U0 O0<br>e2 v1->v0: O2 U2 O2 U2 O2 U2 |
| 6_1 | 6 | handcuff graph | e0 v1->v1: O1 U1<br>e1 v1->v0: O2 U2 O2 U2 O0 U0<br>e2 v0->v0: U1 O1 U1 O1 |
| 6_2 | 6 | handcuff graph | e0 v0->v1: O2 U2<br>e1 v0->v0: O2 U2 O2 U2<br>e2 v1->v1: O1 U1 O0 U0 O1 U1 |
| 6_3 | 6 | handcuff graph | e0 v0->v0: O1 U1 O1 U2 O1<br>e1 v1->v1: U0 O0 U0 O2 U0<br>e2 v1->v0: O0 U1 |
| 6_4 | 6 | handcuff graph | e0 v0->v1: O2 U1<br>e1 v0->v0: U2 O0 U2 O2 U2<br>e2 v1->v1: O1 U0 O1 U1 O1 |
| 6_5 | 6 | theta curve | e0 v0->v1: U1 O1 U1 U2 O2 U2<br>e1 v0->v1: O0 U0 O0<br>e2 v1->v0: O0 U0 O0 |
| 6_6 | 6 | theta curve | e0 v1->v0: O2 U2 O2 U1 O1 U1<br>e1 v1->v0: O0 U0 O0<br>e2 v0->v1: U0 O0 U0 |
| 6_7 | 6 | handcuff graph | e0 v1->v1: O2 U2 O1 U2 O2<br>e1 v0->v1: O2 U0<br>e2 v0->v0: U0 O0 U1 O0 U0 |
| 6_8 | 6 | theta curve | e0 v0->v1: O1 U2 O1<br>e1 v0->v1: O1 U0 O2 U0 O2 U1<br>e2 v1->v0: U1 O0 U1 |
| 6_9 | 6 | theta curve | e0 v1->v0: O0 U2 O2 U0<br>e1 v0->v1: O2 U2<br>e2 v1->v0: O2 U1 O0 U0 O1 U2 |
| 6_10 | 6 | handcuff graph | e0 v1->v0: U0 O0 O1 U1 O2 U2 U0 O0<br>e1 v1->v1: U0 O0<br>e2 v0->v0: O0 U0 |
| 6_11 | 6 | handcuff graph | e0 v0->v1: U0 O0 U2 O2 U0 O0<br>e1 v1->v1: O2 U2<br>e2 v0->v0: U1 O1 U0 O0 |
| 6_12 | 6 | theta curve | e0 v1->v0: O1 U1<br>e1 v1->v0: U1 O1 U0 O2 U2 O0 U1 O1<br>e2 v0->v1: O1 U1 |
| 6_13 | 6 | theta curve | e0 v0->v1: O1 U1<br>e1 v1->v0: U1 O1 O0 U2 O2 U0 U1 O1<br>e2 v1->v0: O1 U1 |
| 6_14 | 6 | theta curve | e0 v0->v1: O2 U2 O2<br>e1 v1->v0: O1 U1 O1 U1 O1 U1<br>e2 v1->v0: U0 O0 U0 |
| 6_15 | 6 | theta curve | e0 v1->v0: O2 U2 O2<br>e1 v1->v0: O1 U1 O1 U1 O1 U1<br>e2 v0->v1: U0 O0 U0 |
| 6_16 | 6 | handcuff graph | e0 v0->v1: O0 U0 O0 U0 O0 U0 O0 U0<br>e1 v0->v0: O2 U2<br>e2 v1->v1: O1 U1 |

Theta curves: 5_2, 5_3, 6_5, 6_6, 6_8, 6_9, 6_12..6_15. Handcuff graphs: 4_1, 5_1, 5_4, 6_1..6_4, 6_7, 6_10, 6_11, 6_16.

## Meshing

`graph_tube` on the three edges. A tube of radius 0.08 around these graphs (length about 15 to 30) is **not
resolved by the marching-tetrahedra grid at a few thousand faces**: 4000 faces gives 18 to 43 components, 8000 faces
1 to 14 components, 10 000 one entry out of 21, 12 000 thirteen of 21, **14 000 all 21** (genus 2, one component, no
self-intersection by `check::self_intersections` and `ccd::intersects`). All runs below start at 14 000 faces
(13 800 to 14 400 as meshed). A thicker tube (radius up to 0.12) was not tried: the clearance is 0.3 to 0.44, so the
blended tubes would nearly touch.

## Flow (first run; superseded by the 6000-step rerun at the end of this file)

`ishii run NAME --faces 14000 --threads 3` (21 jobs in parallel, `runs_tmp/ishii/full3/`): barycenter + area constraints,
`FlowOptions::default()` otherwise (Hs descent, Barnes-Hut, CCD on, remeshing on), `StopRule::Off` in the library and
the stop rule of NOTES-FIDELITY in the example: mean energy over 200 steps differs from the mean over the 200 before
by < 1e-5 relative, cap 3000 steps; no plane attractor. Times are wall minutes on a machine at load 60 to 160 (not
clean). Final inertia moments: eigenvalues (ascending) of the area-weighted second-moment tensor about the area
centroid, with the fixed `eig3`. Column 2: H = handcuff graph, T = theta curve.

| entry | g | faces | genus | steps | status | min | final energy | moment eigenvalues | self-intersects |
|---|---|---|---|---|---|---|---|---|---|
| 4_1 | H | 13876->15016 | 2->2 | 497 | Converged | 38 | 2524.7 | 0.1001 0.1603 0.2711 | false |
| 5_1 | H | 13880->14878 | 2->2 | 3000 | MaxIterations | 127 | 3084.4 | 0.1166 0.1564 0.3148 | false |
| 5_2 | T | 14048->15408 | 2->2 | 1287 | Converged | 82 | 3123.8 | 0.1007 0.1978 0.2317 | false |
| 5_3 | T | 14036->14868 | 2->2 | 1458 | Converged | 89 | 3091.1 | 0.1778 0.1813 0.1825 | false |
| 5_4 | H | 13984->15162 | 2->2 | 1212 | Converged | 80 | 3430.2 | 0.1269 0.1450 0.3872 | false |
| 6_1 | H | 13944->15100 | 2->2 | 3000 | MaxIterations | 132 | 3853.5 | 0.1107 0.1434 0.4069 | false |
| 6_2 | H | 14052->15244 | 2->2 | 2551 | Converged | 122 | 4003.6 | 0.1175 0.1453 0.2980 | false |
| 6_3 | H | 13948->14634 | 2->2 | 1962 | Converged | 110 | 3979.3 | 0.1407 0.1691 0.2564 | false |
| 6_4 | H | 13976->15052 | 2->2 | 2370 | Converged | 117 | 3376.3 | 0.1296 0.2198 0.2352 | false |
| 6_5 | T | 13940->14976 | 2->2 | 2183 | Converged | 116 | 3887.8 | 0.1158 0.1681 0.3222 | false |
| 6_6 | T | 13772->14704 | 2->2 | 1723 | Converged | 104 | 4290.6 | 0.0901 0.1845 0.3207 | false |
| 6_7 | H | 14076->14664 | 2->2 | 2455 | Converged | 120 | 4027.7 | 0.1449 0.1980 0.1996 | false |
| 6_8 | T | 13992->14868 | 2->2 | 3000 | MaxIterations | 128 | 3749.9 | 0.1310 0.1943 0.2766 | false |
| 6_9 | T | 14128->15100 | 2->2 | 1371 | Converged | 91 | 3668.5 | 0.1465 0.1838 0.2938 | false |
| 6_10 | H | 13920->15238 | 2->2 | 3000 | MaxIterations | 130 | 3477.3 | 0.1028 0.2130 0.3579 | false |
| 6_11 | H | 14028->14666 | 2->2 | 2207 | Converged | 117 | 4069.9 | 0.1235 0.1544 0.2947 | false |
| 6_12 | T | 14040->15122 | 2->2 | 2136 | Converged | 115 | 3693.9 | 0.1315 0.2099 0.3042 | false |
| 6_13 | T | 13872->15294 | 2->2 | 1398 | Converged | 91 | 3360.1 | 0.1272 0.2151 0.3026 | false |
| 6_14 | T | 14388->15876 | 2->2 | 1732 | Converged | 107 | 4245.8 | 0.0675 0.1861 0.4832 | false |
| 6_15 | T | 14132->15266 | 2->2 | 1602 | Converged | 99 | 4233.9 | 0.1132 0.1758 0.3925 | false |
| 6_16 | H | 13804->15020 | 2->2 | 1102 | Converged | 81 | 4193.1 | 0.1480 0.1581 0.4095 | false |

Reading.
* All 21 flows ran without a step error or `Stuck`; every final surface has genus 2, one component, and `ccd::intersects`
  is false at the end (CCD was on throughout, so the surfaces stay embedded: the handlebody-knot type is kept).
* 17 converged by the windowed rule (497 to 2551 steps). **5_1, 6_1, 6_8, 6_10 hit the 3000-step cap** (energy still
  drifting down slowly: they are not converged; the numbers are the state at step 3000). Reproducer:
  `ishii run 5_1 --faces 14000 --threads 3 --max-steps 3000`.
* The energies are large (2.5e3 to 4.3e3) because the area is the start tube's area (fixed by the constraint) and
  the energy is not scale invariant: they are comparable between entries only. The tube is not resolved at a few
  thousand faces (see Meshing), so every flow has about 15 000 faces, not a few thousand.
* The ordering by energy (4_1 lowest, 2525; 6_6 4291 and 6_14 4246 among the highest) agrees with nothing in the
  paper (Fig. 16 gives no numbers); I did not compare shapes with the figure's renderings.
* Symmetry: only 5_3 reaches an (almost) triple moment (0.178 0.181 0.183), 6_7 a double one (0.198 0.200); the others
  are anisotropic (largest/smallest 2 to 7). These are indicators, not a classification.
* The 4_1 run stops at 497 steps after the energy plateaus at 2524.7 from step 150 on (2525.1 at 125): a very flat
  start, but it is a converged state under the stated rule, not verified as a global minimizer (one start per entry).
* OBJ start and final meshes: `runs_tmp/ishii/obj/<name>_{start,final}.obj` (not committed).


## Rerun to 6000 steps and determinism check (2026-11-02, `runs_tmp/ishii2/`)

`ishii run NAME --faces 14000 --threads 3 --max-steps 6000 --out runs_tmp/ishii2/NAME`: the previous rerun
(`runs_tmp/ishii/`) with the cap raised from 3000 to 6000. `--out` now also writes `NAME/steps.csv` (one line per step, `step::CSV_HEADER`:
energy, faces, tau, tau_ccd, halvings, remesh counts, ...). All 21 jobs ran in one queue with the 12 genus runs (NOTES-GENUS.md), at
most 63 threads at once; load average 10.8 at start, 2.4 at the end (about 50 during the runs). Times are wall minutes at that load.

| entry | faces | steps | status | min | final energy | moment eigenvalues | self-intersects |
|---|---|---|---|---|---|---|---|
| 4_1 | 13876->15238 | 766 | Converged | 34 | 2526.8 | 0.1003 0.1601 0.2712 | false |
| 5_1 | 13880->14918 | 3358 | Converged | 123 | 3078.8 | 0.1104 0.1613 0.3162 | false |
| 5_2 | 14048->15422 | 3356 | Converged | 122 | 3109.5 | 0.1273 0.1715 0.2169 | false |
| 5_3 | 14036->14854 | 3844 | Converged | 129 | 3076.8 | 0.1774 0.1811 0.1823 | false |
| 5_4 | 13984->15160 | 2210 | Converged | 97 | 3424.4 | 0.1182 0.1467 0.4017 | false |
| 6_1 | 13944->15080 | 3342 | Converged | 129 | 3847.6 | 0.1131 0.1367 0.4108 | false |
| 6_2 | 14052->15350 | 2568 | Converged | 114 | 4005.6 | 0.1180 0.1448 0.2986 | false |
| 6_3 | 13948->14672 | 1653 | Converged | 83 | 3989.2 | 0.1407 0.1690 0.2559 | false |
| 6_4 | 13976->14966 | 1162 | Converged | 58 | 3387.9 | 0.1287 0.2188 0.2374 | false |
| 6_5 | 13940->15010 | 2221 | Converged | 97 | 3878.4 | 0.1145 0.1686 0.3249 | false |
| 6_6 | 13772->14790 | 1872 | Converged | 95 | 4298.7 | 0.0891 0.1854 0.3212 | false |
| 6_7 | 14076->14640 | 1043 | Converged | 55 | 4046.8 | 0.1452 0.1978 0.2001 | false |
| 6_8 | 13992->14856 | 1019 | Converged | 51 | 3843.3 | 0.1006 0.2096 0.2991 | false |
| 6_9 | 14128->14980 | 1191 | Converged | 62 | 3658.9 | 0.1471 0.1843 0.2933 | false |
| 6_10 | 13920->15218 | 3113 | Converged | 124 | 3473.9 | 0.1111 0.2131 0.3542 | false |
| 6_11 | 14028->14650 | 1422 | Converged | 68 | 4080.4 | 0.1237 0.1530 0.2972 | false |
| 6_12 | 14040->15114 | 2491 | Converged | 94 | 3684.9 | 0.1321 0.2093 0.3032 | false |
| 6_13 | 13872->15184 | 2492 | Converged | 90 | 3355.6 | 0.1266 0.2160 0.2999 | false |
| 6_14 | 14388->15768 | 1908 | Converged | 80 | 4242.2 | 0.0682 0.1857 0.4849 | false |
| 6_15 | 14132->15270 | 1088 | Converged | 54 | 4228.2 | 0.1132 0.1743 0.3965 | false |
| 6_16 | 13804->14858 | 1896 | Converged | 75 | 4192.9 | 0.1507 0.1566 0.4076 | false |

* **21 of 21 converged** (windowed rule, 1e-5 over 200 steps). The five that capped at 3000 (5_1, 5_2, 5_3, 6_1, 6_10)
  converged at 3358, 3356, 3844, 3342 and 3113 steps; nothing reaches the new cap. Final surfaces: genus 2, one component,
  `ccd::intersects` false for all.
* Moment spreads stay anisotropic except 5_3 (0.1774 0.1811 0.1823, triple within 3 %) and 6_7 (double, 0.1978 0.2001).
* Wall time 34 to 129 minutes per run.

### Determinism

The previous rerun printed energy and face count every 25 steps (`PROGRESS` lines). For each entry the new run's `steps.csv`
(`energy_after`, `faces`) was compared at every multiple of 25 that both runs have (up to 3000 for the capped ones, up to the
converged step otherwise; 25 to 120 points per entry, 1670 in total; energies compared after rounding to the 7 significant
digits of the log):

**All 21 entries are identical at every common checkpoint** (energy and face count), and every entry that converged before
also stops at the same step with the same status (4_1 766, 5_4 2210, 6_2 2568, 6_3 1653, ... 6_16 1896). So the long flow is
reproducible run to run for this binary, 3 threads, and different machine load (about 50 here against 2 to 150 before); no divergence was
found, hence nothing to locate with a two-run diff. What the check cannot say: (a) only every 25th step was compared, so a
difference that heals within 25 steps would be missed (it is unlikely to heal exactly in a chaotic flow); (b) the table of the "Flow" section (17 converged, 4_1 at 497
steps, 6_8 capped) is the first run, whose logs are lost; it differs from the second run (4_1 at 766, 6_8 at 1019) and from this
one, which agree with each other. Two runs of the same build agree to the 7 printed digits at every checkpoint, so the first-versus-second
difference is not load noise within a build; a different build or start mesh in the first run is the remaining explanation, and
it cannot be checked now. Not investigated further.
