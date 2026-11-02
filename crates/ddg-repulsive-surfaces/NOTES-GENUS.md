# NOTES-GENUS: Fig. 15 run to convergence (g = 0 to 11)

Command: `ishii genus G --threads 3 --max-steps 8000 --out runs_tmp/genus2/gG` (`examples/ishii.rs`).
Same start and options as WP15's `apps genusNN`: `shapes::genus_surface(G, 3000 + 600 G)` (linear arrangement), area and barycenter
constraints, `FlowOptions::default()` (Hs descent, Barnes-Hut, CCD, remeshing), exponents (6, 12), no plane. The stop rule is the one of
`ishii run`: mean energy over 200 steps differs from the mean over the previous 200 by less than 1e-5 relative; cap 8000.
OBJ frames every 50 steps (`frame_NNNNN.obj`), `final.obj`, `steps.csv` per run in `runs_tmp/genus2/gG/`, log `gG.log`.
12 jobs ran in one queue with the 21 Fig. 16 jobs, 3 threads each, at most 63 threads at once (load average about 50).
Energy E is the flow's own value at the last step (Barnes-Hut, not the exact all-pairs value of NOTES-FIDELITY). Moments: eigenvalues of
the area-weighted second-moment tensor about the area centroid (fixed `eig3`), ascending.

| g | faces | steps | status | min | E | moments | ratios | embedded (`ccd::intersects`) |
|---|---|---|---|---|---|---|---|---|
| 0 | 2976 -> 3264 | 521 | Converged | 3 | 2.468 | 0.3682 0.3699 0.3707 | 1.005 1.002 | yes |
| 1 | 3676 -> 4008 | 660 | Converged | 5 | 31.85 | 0.1402 0.5538 0.5659 | 3.95 1.02 | yes |
| 2 | 4196 -> 4538 | 1013 | Converged | 9 | 38.33 | 0.5134 0.5671 0.9504 | 1.10 1.68 | yes |
| 3 | 4924 -> 5118 | 1167 | Converged | 12 | 40.09 | 0.8774 0.9102 0.9526 | 1.04 1.05 | yes |
| 4 | 5388 -> 5706 | 1136 | Converged | 14 | 46.78 | 1.0589 1.1258 1.2074 | 1.06 1.07 | yes |
| 5 | 6028 -> 6306 | 1364 | Converged | 19 | 48.65 | 1.3222 1.3242 1.3407 | 1.00 1.01 | yes |
| 6 | 6968 -> 7234 | 2500 | Converged | 41 | 55.98 | 1.4432 1.4883 1.7061 | 1.03 1.15 | yes |
| 7 | 7184 -> 7822 | 1509 | Converged | 27 | 60.89 | 1.4589 1.4836 2.3644 | 1.02 1.59 | yes |
| 8 | 6428 -> 6988 | 1805 | Converged | 31 | 71.23 | 1.4860 1.5500 2.9983 | 1.04 1.93 | yes |
| 9 | 9260 -> 9988 | 1493 | Converged | 36 | 68.95 | 1.6113 1.6392 3.4539 | 1.02 2.11 | yes |
| 10 | 8240 -> 9228 | 1512 | Converged | 35 | 80.05 | 1.6333 1.7969 3.8807 | 1.10 2.16 | yes |
| 11 | 10136 -> 11058 | 1711 | Converged | 48 | 81.90 | 1.7327 2.0335 4.1232 | 1.17 2.03 | yes |

Reading.
* All 12 converged by the windowed rule (521 to 2500 steps); none reached the 8000 cap. All final genera equal g, all surfaces embedded. The
  rule can stop on a long plateau (g = 6: 56.0 at 2000 steps, 55.98 at 2500), see NOTES-FIDELITY on slow saddles; it was not checked
  whether lower states exist (one start per genus, as in the paper's setup).
* g = 3 reaches the near-tetrahedral shape of NOTES-FIDELITY (moments within 9 %, E 40.09 against 39.9 there; it left the linear
  arrangement between steps 500 and 750: F 52.1 at 500, 40.3 at 750), in 1167 steps. g = 2 repeats NOTES-FIDELITY (1013 steps, moments 0.513 0.567 0.950).
* Symmetry (moments, multiplicity at 3 %): g = 0 triple (round); g = 1 double (torus of revolution, 0.554 0.566); g = 5 triple (1 %,
  the cube); g = 3 nearly triple (9 % spread, tetrahedron-like), g = 4 not (14 %); g = 6 double within 3 %; g = 7 double
  (1.459 1.484, long axis 2.36), g = 9 double (1.61 1.64, long axis 3.45), i.e. prism-like stacks; g = 8, 10, 11 are elongated with
  no multiplicity within 3 %. The paper's symmetric shapes (prism for g = 2, cube for g = 5, ...) appear for g = 0, 1, 5 and approximately
  for 3; g = 2 is again a bent saddle.
* WP15's 800-step shapes for g = 7 to 11 were unfinished: the energies fell (flow energy now against exact energy then, so only
  indicative: g = 7 65.6 -> 60.9, g = 11 90.8 -> 81.9) and g = 7 now has a double moment; g = 8 to 11 are still elongated (one moment
  1.9 to 2.2 times the next), as final states of this flow, not as an unfinished run.
* Files: frames and finals are in `runs_tmp/genus2/` (not committed); the sheet is `showcase/media/repulsive_surfaces_genus.png`.
