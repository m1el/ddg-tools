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

## Flow

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
