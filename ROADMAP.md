# Roadmap and known gaps

State of 2026-10-14, branch `intrinsic-triangulations`. Numbers and evidence are in the documents this
file points to; this is the list of what is scheduled, what was decided, and what is known to be
open. Add to it; move an entry to the crate's README or SPEC when it is closed.

## Paused and finished (2026-10-17)

Work started on the user's request ("2. implement mapping to the sphere. 4. implement
applications listed in the paper"), checkpointed as unmerged branches of this repository, each
with a notes file saying what is done, what remains and how to resume:

* Done and merged 2026-10-17: BFF mapping to the sphere (`wip-bff-sphere`, finished in
  `44412af`; see "Landed" and `crates/ddg-bff/README.md`).
* Done (2026-10-17): `wip-apps-5123` (6e59c99), applications 5.1–5.3. Merged with the next one
  into `crates/ddg-intrinsic-apps`; see "Landed".
* Done (2026-10-17): `wip-apps-amr-vector` (f730a89), applications 5.4–5.5 and
  `ddg_intrinsic::refine_with`; ddg-py part compiled, registered and tested. See "Landed".
* Done: the `check-formulas.sh --build` failure both app branches reported ("unknown
  Ddg.movedLength" in `Intrinsic/Move.lean`) was the stale Lean cache copied into their clones;
  a build with fresh project outputs passes.
`master` includes all of this work (moved up 2026-10-17); the agents' worktrees and branches,
including three misplaced in the authors' demo repository, were removed the same day.

## Scheduled

Nothing (item 1, the feedback loop of the cones solver, landed on 2026-10-26; see below).

## Landed since the first version of this file

* **Feedback loop of the cones solver** (2026-10-26, branch `cones-feedback`; formerly Scheduled
  1; `crates/ddg-optimal-cones/README.md`, "The feedback loop" and "Published UVs with the
  feedback loop"; `crates/ddg-uv/README.md`, "The cone feedback loop"): `Feedback`,
  `scale_mismatch`, `solve_with_feedback` in the solver crate, `UnwrapOptions::cone_feedback` in
  the pipeline, `cone_feedback=` in Python, "Cone Feedback" in Blender. The solver's linearized
  `u` is compared with the log area factor the flattening realizes (first on the surface without
  cones, since a cone at a needle tip hides it), and the distortion weight goes to zero where
  they differ by more than 20 log units. Judged by the 23-mesh `D*` comparison: `camelshotnoise`
  C1 goes from `+982 %` to `−20.2 %` (`u = 0`) and from `+57 944 %` to `−27.1 %` (optimized),
  exactly the hand-set tips; every other C1 row identical, the 11 reproducing meshes untouched.
  On by default (identical maps on 42 meshes without needles, catches needles the spike check
  misses; +29 % median time on meshes with boundary, 0 on closed ones). Lean:
  `cone-feedback-mismatch` proved, `cone-feedback-rule` cited as a design choice. Open: closed
  surfaces get no cone-free probe (no flattening without cones), so a needle on a closed piece
  whose tip received a cone is found only by the spike check; every extra round runs the
  reference method (zero weights disable the factor-once method).

* **Mapping to the sphere** (2026-10-17, `crates/ddg-bff`, README "Mapping to the sphere"):
  `map_to_sphere` maps a closed connected genus-0 mesh conformally to the unit sphere (port of
  `BFF::mapToSphere`: punctured disk, radius search, stereographic projection, Möbius centering
  after Baden, Crane, Kazhdan 2018). Default pole: the most regular vertex star (the reference's
  vertex 0 gives flipped triangles on a noisy ellipsoid). Differential test against the
  reference's own program; Lean `DdgProofs/Bff/Sphere.lean`; Python `ddg_tools.map_to_sphere`;
  Blender operator "Conformal Map to Sphere". Open: injectivity is reported, not guaranteed
  (brain: 38 flips, as the reference); the 1e-4 difference of the reference's in-place puncture
  is unexplained.
* **Applications of intrinsic triangulations** (2026-10-17, new crate `crates/ddg-intrinsic-apps`,
  README and SPEC N30–N49; Python `geodesic_distance`, `steiner_tree`,
  `smoothest_direction_field`, `adaptive_solve`): 5.1 Steiner trees (1–6 % shorter than the
  extrinsic MST), 5.2 FEM condition numbers (4-1 subdivision orders of magnitude worse than iDR
  and iODT on bad inputs), 5.3 the heat method against an exact polyhedral distance written for
  the crate, 5.4 adaptive refinement (residual estimator, bulk marking, `refine_with`), 5.5
  connection Laplacians, harmonic fields and the vector maximum principle, smoothest direction
  fields. Open:
  * the extrinsic ODT of Fig. 17 (5.2 compares against the 4-1 subdivision only) and Fig. 16
    (`L²` error of a Poisson problem) are not reproduced;
  * cut graphs for genus above zero (5.1 cuts a disk only on genus-zero surfaces);
  * Blender: nothing uses the applications; the 5.1 cut as seams needs its intrinsic edges
    inserted into the mesh first;
  * the vector-field study finds Fig. 21's large energy drop (input to iDT to iDR) only on the
    shot-noise camel; the paper's mesh is not in the data set (crate README, 5.5);
  * timings: done on 2026-10-19 (quiet machine, single thread; README 5.3 to 5.5). The AMR
    problem rows, the disk maps and the vector-field times were re-measured; the parallel
    other `apps_study` tables were not touched.

* **Disconnected meshes in the cones solver** (2026-10-17, `crates/ddg-optimal-cones`, README
  "Disconnected meshes"): every entry point accepts several components and returns the exact
  solution of the canonical problem on the whole mesh (masses divided by the total area, one
  constant of `u` and one Gauss–Bonnet sum per closed component), solved per component and
  certified as a whole, with a report per component. Lean: `component-separation`,
  `component-lambda-rescaling` (`DdgProofs/Components.lean`).
* **`reoptimize_cones`** (2026-10-17, same crate, README "Re-optimizing given cones"; Python
  `ddg_tools.reoptimize_cones`): the least-distortion angles for given cone positions, by a direct
  solve of the linear optimality system or, up to 32 free sites, the reduced system of the
  factor-once method; the generic solver when a bound is active. On the 23 published cone files
  it never gives more distortion than the published angles.

* **The optimal Delaunay variant** of the intrinsic triangulations paper (2026-10-14,
  `crates/ddg-intrinsic`, SPEC N19 to N25, README "The optimal Delaunay triangulation"):
  `move_vertex` (Alg. 10: an inserted vertex inside its flat star, an exact test of the star
  instead of the paper's height bound, refused beyond it; boundary vertices along the boundary)
  and `optimal_delaunay` (target-length splits, Chen and Holst's update with their line search,
  flips; options for the target, step, stopping, boundary and centres; a report with the
  iterations and the angle and area distribution). Lean: four new formula ids proved or cited,
  an oracle task for the move and its replay. Measured (12 meshes): at the vertex count of
  refinement the areas are far more even (coefficient of variation 0.31 to 0.58 against 0.31 to
  2.9 at 25°), the smallest angles are those of the input (an input vertex cannot move), the
  cotan Laplacian is mostly better conditioned than after refinement; it costs 3 to 40 times
  the time of refinement. In the cut flow (`flow_on_bad_meshes.rs`, TREAT=odt, odt15) it is an explicit
  option of `ddg_uv::cutting::IntrinsicOptions`; the default stays Delaunay flips (see the
  README for the comparison).

* **The way back of a segment along an input edge** (`crates/ddg-intrinsic`, SPEC N29): it failed
  on 2 of 900 cuts of the remeshing study (bob's control mesh, start plane = mirror plane;
  `remesh_study repro bob_control ka4 dirichlet 1 1`, and `ka2 dirichlet 0.5 1`), because the
  ends on the input edge were placed in one of its two faces by rounding. Point queries now return
  a place within `1e-12` of an input edge as a point of that edge, and a segment along a side
  that lies along an input edge is not traced; both reproducers come back, and
  `tests/way_back_along_input_edges.rs` covers the class.
* Flipped UV triangles after the Blender Apply: the library's triangulation of the split polygons
  is stored and used by the unwrap (`blender/README.md`).
* **Flipped UV triangles with iODT at a large cut offset** (2026-10-17, `crates/ddg-uv`, README
  "A large cut offset, and the optimal Delaunay triangulation"): the Blender spot atlas (8
  clusters, weight 2, 60 steps, iODT, Cut Offset 0.2) had 4 of 7354 flipped UV triangles; the
  stage is the flattening of the chart (the triangles have the orientation of their parent in
  space; the same with cones off), the cause the order of two rules of the way back: the rule of
  `min_offset` moved every point first, and the merge of bends measured its displacement at the
  moved places, so bends the rule pushed off a side stayed (iODT brings about 7 times as many
  bends as flips) and became spikes with slivers and shrunk parts next to them. The merge now
  measures at the found places, and removes the excursions across a side whose points the rule
  moved. `flip_study`, 44 runs per offset (spot, blub, bob, camelbad atlases and circles; flips,
  refined, both with iODT): 123 flipped triangles at offset 0.2 before, 3 now; 55 and 57 at
  1e-3. The same mechanism made the two refined cases at offset 0.2 of the Blender test (one
  flipped triangle each, a bend "0.103 of an edge" from its chord: that was the rule's motion);
  they are 0 now, and the Blender tests assert zero flips on all five variants again (`-odt`
  included) and run `the_cut_offset_moves_points_away_from_vertices` with iODT.
* **Intrinsic triangulations** (`crates/ddg-intrinsic`, contract `SPEC.md` there): signposts,
  flips, Delaunay flipping, insertion, refinement with boundary, overlay, simplicial export,
  the way back for a cut; Lean oracle and proofs for its twelve formula ids. The cut flow runs on
  the Delaunay flips by default in `ddg_uv::cutting`, the bindings and the Blender operators;
  refinement and the mesh's own triangulation are options. Measured answer on the angle bound
  (`crates/ddg-intrinsic/examples/flow_on_bad_meshes.rs`, logs in `runs_tmp/intrinsic-study/`,
  not in the repository): flipping alone brings the systematically obtuse meshes back to
  regular-grid quality in the evaluation, the redistancing and the exact laws; refinement adds
  30 % to 140 % vertices and helps only where the input has few vertices for thin geometry; 30°
  buys nothing clear over 25°.
* Step and stop rule of the cut flow after the first version: SPEC D26 (the local limit with the
  length step on the mesh), D27 (a rising energy is its own stop reason and sets the run back),
  D29 (non-finite numbers are errors), D30 (faces with a moved crossing leave the mean slope),
  D31 (convergence by motion needs a small predicted motion).

## Open after that work (2026-10-14)

* The merge of bends of a cut that comes back from an intrinsic triangulation has a default
  tolerance of 0.1 of an edge and moves the cut by up to that (measured on the cut the flow
  found; the rule of `min_offset` moves points by its own amount on top); smaller values let
  flipped UV slivers return (spot atlas, weight 2, offset 0.2, iODT: 0.05 leaves none, 0.02 leaves 12 and 7).
* The rule of `min_offset` moves every point of the cut on its own. At a large offset (0.2) it
  still makes, after the fix of 2026-10-17: a spike of an excursion whose bend lies further than
  the merge tolerance from its chord (blub atlas, flips with iODT: 1 flipped triangle), and of a
  triple point close to a side with two of its pieces to that side (spot atlas at ripple
  wavelength 1.0, flips with iODT: 1; the spot atlas after 20 steps at weights 2 and 10,
  refined without iODT: 1 and 2; a triple point cannot be merged, it would have to move into
  the neighbouring triangle with new crossings); and it refuses (`InvalidDescription`, "edge points ... are
  both at t = 0.8") where two crossings of one edge close to the same vertex end at the same
  parameter: 10 of the 44 `flip_study` runs at 0.2 (blub, camelbad and its atlas, on every
  triangulation), with or without iODT. A rule that keeps the order of the points on an edge would
  remove these refusals (not the spikes). Not seen in the Blender tests (spot, offset 0.2).
* `flip_study` at the default offset: the atlases of blub and camelbad at ripple wavelength 0.5
  (scaled Dirichlet) flip 3 to 15 UV triangles on every triangulation, iODT or not (strongly
  shrunk charts, fans around a vertex); not studied further. The bindings do not report
  `DescribeReport::n_removed_excursions` yet (`crates/ddg-py/src/cutting.rs`, next to
  `n_merged_nodes`).
* Refused, not handled, by the insertion of `ddg_uv::cutting`: a region that fits inside one
  input triangle, a cut exactly through a vertex, a piece along a diagonal; in the add-on, a
  loop of junctions inside one polygon that touches no edge point.
* The simplicial export loses the angle bound where it must split edges (thin tubes: restored
  by `ExportOptions::refine`; spikes: not), `ddg-intrinsic` SPEC N15. Termination of the
  refinement is argued for closed surfaces and only observed with boundary and spikes (N1, N4).
* From the study: long energy runs on flipped sheared grids still do not settle; the
  balanced-latitude bias on latitude-longitude spheres is unchanged; one outlier of the speed on
  1:4 cells refined to 30° is unexplained; the energy of the same start differs between the
  input mesh and the flipped one on bad meshes (on a flat mesh this is the cut, not the
  evaluation: the zero set of the interpolant differs; the distortion term on curved meshes was
  not compared).
* After D30 the 600-step run on camelshotnoise (stop rule off) ended at 3475 where it had
  reached 2967. Not a regression: that run is chaotic. With length weights perturbed by 1e-5 to
  1e-3 of their value, the code before D30 ends anywhere between 2967 and 3990 and the current
  code between 3214 and 3620 (six runs each, logs in `runs_tmp/camelshotnoise-d30/`); with the
  default stop rule both return the same state (step 16, energy 4524.7). Some README and SPEC §6
  tables carry numbers from before D30 and say so.
* The Blender session cache (a session kept between Start, Flow and Apply) is not tested across
  file loads or object deletion; `blender/tests/bench_overhead.py` still does not finish.
* The comparison with the authors' program was repeated on 2026-10-17 with the code of D32
  (below, and `crates/ddg-variational-cuts/README.md`).

* From the optimal Delaunay work (2026-10-14): a move beyond the star of a vertex (the paper's
  general route: remove the vertex, retriangulate, insert again) is not implemented (SPEC N20);
  iODT cannot improve angles at input vertices, so on meshes whose bad angles sit there
  (latitude-longitude poles, cells 1:16, spikes) it leaves them; its default target length (the
  mean input edge) adds 1 to 2 times the input's vertex count on curved meshes and lowers the
  smallest angle of an already regular mesh (grid 45° to 39°, icosphere 54° to 29°); with it the
  flow stalled once at its start energy (latitude-longitude 41x82, `α_L = 3`), which flips and
  refinement do not; the bindings and the
  Blender add-on do not expose the option yet; the accuracy of signposts after many operations
  with thin intermediate triangles (SPEC N25) is bounded by measurement, not by an argument.

## Decided (2026-10-14)

* Two regions: the default stages are the ones that meet the exact laws (SPEC D12, D14–D16,
  D22); the paper's stages stay selectable (`FlowOptions::paper`).
* More than two regions: the recommended stages stay the paper's (`FlowOptions::for_regions`),
  although the evidence became mixed after D22 (SPEC §8.5, "After D22").
* No remeshing except inserting the cut (paper Sec. 6.1 is not implemented). Measured since
  (`crates/ddg-intrinsic/README.md`, "Remeshing for the cut flow"; cut SPEC D32, intrinsic SPEC
  N28; `crates/ddg-uv/examples/remesh_study.rs`, 1 260 runs on 10 meshes): uniform,
  curvature-adaptive (the paper's rule, at matched vertex counts) and iODT refinement of the
  intrinsic triangulation, with the curvature at the input vertices or spread from the input
  faces, against the flips. Refinement does not make the flow converge (1.2 to 3 times the steps,
  never by motion), costs 11 to 70 times per step, and is even or worse on spot, bob, blub and the
  camels; only spread-curvature refinement of meshes with a few hundred vertices is better
  (−2 to −3 % median energy on a fine yardstick). The decision stands; no option was added.

## Known gaps: the cut flow (`crates/ddg-variational-cuts`)

Method
* The robust step rule (SPEC D22) rests on runs, not on a stability analysis; the quantile level
  0.9 was read from tables. The energy is less monotone with it (spot, Hencky, `α_L = 2`: rises
  in 23 of 150 steps, up to 8.7e-3). With the default stages the rate after the local limit is
  not always a descent direction (SPEC D25a, machine-checked counterexample).
* On irregular meshes the cut never comes to rest: it keeps moving by 2–4 % of an edge per step
  while the energy stands still; runs end by the energy criterion (SPEC §4, "Stopping").
* Dirichlet energy with too small a length weight does not settle (spot, `α_L = 1`: a period-two
  cycle; README, Fig. 5 with 0.18: two pieces). Meshes with spikes of curvature stay noisy
  (camelshotnoise: energy rises in 256 of 600 steps).
* Networks: the length of a relaxing honeycomb is not monotone and junctions stay a few degrees
  off 120° at the mesh widths tested, with either stage set; cell areas drift by 1–3 % (SPEC
  §8.5, §8.6). Which stage set is better for atlases is open (D20).
* Obtuse meshes: with every triangle obtuse (sheared grids from 135°) redistancing and the
  evaluation degrade by large factors; the default flow still meets the circle law to 1–2.5 %
  but the cut is rough and runs do not settle. Flipping to a Delaunay triangulation removes the
  measured damage but not the unsettled runs; why 117° is the noisiest case, and the bias on
  latitude–longitude spheres, are unexplained (README "Mesh quality"; logs of the study were not
  kept in the repository).
* The smoothing length of the flux recovery (`c = 0.5`, SPEC D12) comes from exact solutions on
  spheres and a flat annulus; where the normal derivative varies over about five edges it biases
  it by up to 15 %.
* Symmetry holds only to the discretization error without symmetry classes when quad diagonals
  tie (SPEC §3.11); a symmetry that permutes regions is not supported.
* A pointwise energy can be given in Rust only (not from Python).

Verification
* The Lean oracle checks one step. Its own energy-decrease check for several regions covers
  three regions on one mesh, and its finite-difference check for several regions has no triple
  point and a 10 % tolerance. Its sort for the quantile reverses equal values (SPEC D25e).
* Cited, not machine-checked: the Yamabe equation, the shape derivatives, Céa's method, the
  level-set evolution, the ripple linearization, the first variation of a network at a junction,
  the stiffness matrix of the weighted Laplacian as an integral, convergence of fast marching
  (`lean/proofs/CITED.md`). Hypotheses of proved statements that are not derived from the mesh:
  angle sums of 2π inside a flat face. Nothing is proved about rounding, cut extraction, the
  combinatorics of the cut mesh, or the marching order; Rust and Lean are linked by the formula
  tags and by reading.

Comparison with the authors' program (black box, `runs_tmp/blackbox/`, not in the repository)
* Repeated on 2026-10-17 (code of D32; caller's triangles and intrinsic Delaunay; Hencky and
  Dirichlet; the reference also run to rest; its UV maps measured); tables in
  `crates/ddg-variational-cuts/README.md`, "Comparison with the authors' program". Measured by
  the flattening of each cut, we are below the reference after 300 steps everywhere (4–14 %) and
  below the reference at rest on sphere, bob and blub (2–3 %), but 2.2 % above it on spot.
  Started from the reference's cut at rest, our flow ends 0.4–5 % below it (spot 3–4.5 %; spot,
  bob, blub; both length weights, stage sets and triangulations) and 5–8 % below our own spot
  runs: the gap was a local minimum. Open: our start from 8 clusters of normals stops in worse basins than the
  reference's long runs; the energy is not monotone along a run (single steps rise, up to 4 %).
* The reference's Hencky value is about 10³ times `∫u²` on the unit-area surface (930–1110
  against its own UV map; identical under scaling, unchanged under 1:4 subdivision). The factor
  itself is not explained from outside; the length-weight correspondence still uses 900 (with
  10³ the weights would be 10 % smaller). The old spread "440 to 1010" was the face-label
  transfer losing small charts, which also makes the evaluation of its cuts in our energy favour
  us; the flattening yardstick does not transfer anything.
* The reference's Dirichlet runs do not come to rest with its default fixed step (the objective
  oscillates by a factor 2 over 8 000 steps); its line search and other step sizes were not
  tried, nor its local-scale, symmetry, area and visibility options.
* Pictures: `runs_tmp/blackbox/pics_2026-10-17/contact_sheet.png` (all shapes, both energies).

## Known gaps: pipeline, bindings, Blender

* The Blender add-on keeps the library's triangulation of the polygons an applied cut split in a
  mesh property; a mesh whose vertices or faces were added or deleted since falls back to
  Blender's triangulation (with a warning), which can flip a few thin UV triangles again
  (`blender/README.md`, "stored triangulation").
* The geodesic-circle start is two-region only and does not reach across a seam; symmetry
  classes skip vertices that a seam splits (`crates/ddg-uv/README.md`).
* Polygons with many collinear boundary vertices are rejected ("a diagonal of a face coincides
  with another edge"); the insertion could in principle produce such a diagonal itself (no test
  has).
* Atlas quality was not evaluated in the pipeline tests: they run tens of steps, and with cones
  off the distortion of a chart can be high (bunny 4.4 against 0.29 with optimal cones).
* Blender: meshes with shape keys are refused; an orthographic camera is treated as a point; the
  facing term of the camera penalty uses object-space normals (non-uniform scale untested);
  penalty, visibility and symmetry options were tested on primitives only; the advanced options
  sit behind a toggle because `layout.panel` crashed Blender 5.0.1 in the redo popup; every
  operator call builds a new session; several regions were not benchmarked in Blender;
  `blender/tests/bench_overhead.py` does not finish within ten minutes.
* The penalty strength of the Flow operator is multiplied by the length weight (a choice made
  here, stated in the tooltip).

* `ddg-uv` solves every piece (component after cutting the seams) with its masses normalized by
  the piece's own area and the user's `λ`, which is the whole cut mesh's canonical problem with
  `λ · A_c / A` per piece: every chart gets the cones it would get alone. The cones solver now
  solves disconnected meshes with the whole-area normalization. Decided 2026-10-17 by the user
  ("uv pipeline is doing good work normalizing each mesh"): the pipeline keeps the per-piece
  normalization (`crates/ddg-optimal-cones/README.md`, "Disconnected meshes").
* Blender has no operator for re-optimizing given cones. The Unwrap operator's "pinned" vertex
  group (cones at no cost) plus a "forbidden" group holding every other vertex already asks the
  pipeline for that problem (with the pipeline's angle bounds: positive only by default, capped
  below `2π`), and then cuts and flattens with the result; a "cones only at the pinned vertices"
  option would build the forbidden set itself. Not added: the add-on's job is the UV map, and
  angles alone have no place in it.

## Developability of Triangle Meshes (`crates/ddg-developability`, branch `developability`, 2026-10-19)

Landed: the six energies with gradients, L-BFGS with the Lewis–Overton line search, gradient
descent, 4-1 subdivision, small-angle remeshing, the flow driver with exact face provenance, cut
and flatten (Sec. 4.4) with the angle-defect terminals (SPEC D30), the admissibility floor (D31,
D32); Lean: 43 `dev-` ids (30 proved, 13 cited), the Float oracle (two 400k-case fuzz runs, 0
disagreements); Python `make_developable`, `developability`, `developable_seams`,
`developable_unwrap`; Blender "Make Developable", "Show Developability", "Developable Seams + UV".
Clean room: the authors' GPL code is quarantined (mode 000) and was never read. Contract and
studies: `crates/ddg-developability/SPEC.md` (§12 holds the studies), README, `blender/README.md`
("Measured").

Open:
* **Global UV overlap.** Done (2026-10-19): `ddg-uv` splits overlapping charts automatically
  (`UnwrapOptions::split_overlaps`, on in "Developable Seams + UV", off in Conformal Unwrap).
  At the default flow the overlap falls from 1-40 % to at most 0.33 %, 0 mirrored triangles; some
  raw dense input meshes keep 1.4-6.5 % (budget of 32 cuts, undone cuts). Open: try a second pair of
  layers when a cut is undone; rerun only the changed charts instead of the whole pipeline
  (3 to 10 times the unsplit time). See `crates/ddg-uv/README.md`.
* **Black-box comparison with the authors' program**: done 2026-10-19
  (`crates/ddg-developability/COMPARISON.md`; libigl supplied by the user, a GL header shim the
  user approved; no reading of their code). Same mesh, same energy: their values equal ours to
  all printed digits once boundary vertices are excluded (they omit them; we count them by
  default, D8, the user's choice); their L∞ covariance and L∞ combinatorial are ½ × ours; their
  L2 combinatorial is a different normalization (per star, arcs ≥ 2 faces). Our flows are lower
  in 25 of 30 intrinsic/E_λ runs at 200 steps, 72/72 runs valid against 64/72 of theirs; our
  cut + flatten distorts far less and flips nothing. The slow L∞ covariance
  (`CovarianceMax`) flow (spot control, no remesh: 1.70 after 200 iterations against their 0.35)
  is explained (2026-10-20, two-team clean room, COMPARISON.md §3.3): their line search starts at
  2 × the previous step (our new opt-in `FirstTrial::GrowPrevious`, D37: 1.10) and has no step
  cap (with a cap of 2 ℓ̄: 0.062 at their level of distortion). D37 stays opt-in (the user,
  2026-10-26: the more reliable rule is the default). Open: the `develop` example's `RESULT energy=` shows the `--energy` default when `--rounds` sets
  another energy.
* **Timings** measured on a quiet machine, single thread, one job at a time (energy, eig and flow
  benches with 4/16/64 threads on 2026-10-19; the Blender end-to-end table of `blender/README.md`
  (all 57 rows, including the long flows N x 5000, 0 sub, 3000 and the converged control meshes) and
  the 132 jobs of SPEC 12.5 on 2026-10-26, pinned to one core, load 1.1 to 2.8). SPEC 12.1 to 12.4
  hold only wall-clock totals of parallel studies and were not re-run. Open: the thread scaling of
  the energy is only 2 to 10 times at 64 threads (SPEC 12.7), not investigated.
* The dense meshes do not converge at the default budget (every round ends by
  `max_iterations`); only the control meshes approach zero energy.
* Done (SPEC D34): a converged flow ends `Stationary` (was `line_search_failed`).
* Long flows on control meshes leave faces below the small angle that remeshing cannot repair
  (`bad_faces`): halved by D35 (62 -> 35 over the 45 robustness flows), not gone; none on coarse
  uniform inputs.
* **Quality against the paper's figures**: investigated, `crates/ddg-developability/QUALITY.md`.
  The cause is the input (a subdivision cage instead of a coarse uniform mesh: thin parts fold flat
  into double sheets) and `E_λ` in the coarse phase (it flattens coarse meshes); the paper's recipe
  (uniform ~1k triangles, `E_P` then `E_λ`, two subdivisions) gives its look (spot: Hausdorff 0.045
  of the diagonal, volume 0.96). The default flow is now that recipe (D36). Done since: the Blender
  operator and Python have a coarse energy (`E_P` by default), Blender's own Remesh (voxel) +
  Decimate to ~1k triangles gives inputs as good as the study's (`blender/README.md`), and the
  showcase runs the paper's recipe. The "short
  creases" of the in-place flows are mostly the hatching of seams that cross triangles in the
  showcase's edge shading.

## Repulsive Curves (`crates/ddg-repulsive-curves`, branch `repulsive-curves`, 2026-10-31)

Landed:
* the edge-pair tangent-point energy (exact and Barnes–Hut) with its differential;
* the fractional Sobolev metric, and the L², H¹, H² metrics for comparison;
* every constraint and penalty of Sec. 8.1, length schedules and subdivision;
* the collision-guarded and normalized Armijo line searches with CCD;
* L-BFGS and NCG for the soft problems;
* the block cluster tree and multigrid on curve networks, and `Backend::Auto`;
* the near-pair quadrature (not in the paper, on by default; SPEC §3.4);
* a reader for the authors' scenes, and the application scenes in `showcase/`;
* Lean: 73 `rc-` ids (55 proved, 18 cited), the Float oracle and golden fixtures;
* Python `RepulsiveSession` and `repulsive_flow`; Blender "Repulsive Flow (Live)" and "(Steps)".

Fig. 20 is beaten at every size (README, SPEC §17k and §17l). The Sec. 7 comparison on all 228
knots is in SPEC §17s (results of record): 0 stuck runs for Hs, soft Hs NCG and soft Hs L-BFGS,
where the paper has 6.5 %, 3.3 % and 14.8 %. The L², H¹ and H² failure rates are within 4 points
of the paper's.

Open:
* **Oracle fuzz**: `fuzz_against_the_oracle` (ignored) reports 47 disagreements in the
  chord-iteration count of one step (0 here, 1 in the oracle). They predate the near-pair work
  and are not explained (SPEC §17o).
* `Session::step` ignores `FlowOptions::time_limit`; only `Session::run` honours it.
* The parameter sweeps of SPEC §14.3 (quadrature, projection tolerance, Barnes–Hut, BCT and
  multigrid parameters) were not run; the defaults rest on the measurements of §17k–§17l.
* Multigrid:
  * the block cluster tree product is slow on the 4,296-vertex vascular scene, and its
    multigrid diverges there (SPEC §17g);
  * the solve and the hierarchy build stop scaling past 8–16 threads (capped at 16,
    `MgOptions::max_threads`);
  * knot 96 is 2–5× above the paper's multigrid curve at one thread (§17k).
* Evaluation leftovers (§17s):
  * one soft Hs L-BFGS run (trefoil0016) ends 2.5 % above its target;
  * one soft H¹ NCG run (knot0001) ends 1.82× its first-evaluation energy.
* The authors' `TangentInterpolation` scene collapses an end edge after about 60 steps (the
  showcase uses its own scene; SPEC §17j).
* `showcase/media/` is not committed.

## Integer Coordinates and FlipOut (`crates/ddg-integer-coords`, `crates/ddg-flip-geodesics`, branch `integer-coordinates`, 2026-10-31)

Contracts: `crates/ddg-integer-coords/SPEC.md` (decisions I1–I28) and
`crates/ddg-flip-geodesics/SPEC.md` (F1–F30).

Landed:
* the shared connectivity `ddg_intrinsic::Complex` and the trait `Triangulation`, implemented by
  both representations; Delaunay flipping and refinement are written once over the trait (IT N50–N52);
* integer coordinates:
  * flips, edge and face splits, vertex removal and moves; geometric tracing; the common
    subdivision with `P₀`, `P₁` and point queries;
  * refinement (the paper's rules) and `L²` transfer;
  * `certify()`, and exact differential tests against the signposts;
* FlipOut:
  * paths, loops and networks, waypoints and edge sets, extraction and the `T⁰` certificate;
  * on both representations, bitwise identical;
  * the shortening loop is within 0.97–1.42× of geometry-central's, with the same lengths on all
    4498 compared pairs;
* I27 resolved (aa096cb): a vertex with a straight star (made by an edge split) where removal by flips is
  stuck, and an inserted boundary vertex, are removed by merging the two halves
  (`Complex::merge_edge`, the inverse of `split_edge`, IT N53; integers by `ic-merge-edge-normal`,
  cited until tier C); both representations alike. Split vertices stuck: 860 of 2588 right after
  the split → 0; after Delaunay flipping 1898 → 74, all at stars with an edge from the vertex to
  itself;
* Lean:
  * the integer formulas (tier A), the one-face model (tier B) and Thm. D.1 (B′);
  * the integer oracle with 23 golden fixtures.

Open:
* FlipOut's loop is still 1.0–1.4× gc: F2's corner sums need the corner angles at every vertex a
  replacement puts on the path, where gc reads signposts (FlipOut SPEC §9.1).
* I28: the integer edge split refuses a point within rounding of a crossing, so refinement skips
  that face where the signposts split and snap (1 of 256 differential runs). Undecided whether to
  snap to the crossing. The export and the optimal Delaunay iteration split such an edge at
  `½ ± 2⁻¹⁶` instead (I29; torus 5×3).
* I29 **done** (2026-11-01, bb5f92c, d3c3999, 0e6c3a6, 7bf6014, e1ebbde, and the docs commit after them):
  the export, the way back of segments and the optimal Delaunay iteration on integer coordinates
  (written once over `Exportable` / `OdtMoves`, signposts bitwise unchanged); `ddg-uv`, Python
  `CutSession` and the Blender cut default to integer coordinates (`representation`). The cut is
  the signposts' to rounding on every measured case; the optimal Delaunay setup costs 1.2–1.3× on
  integers (3.4–4.2× before moves inside one input face skipped the arcs; IC SPEC §4.10, §4.11).
* FlipOut R5: refinement still refuses constrained edges, so refinement that keeps FlipOut paths as
  edges (F9) cannot run.
* Fig. 14: the trends are reproduced, but the factors are not (2.7× and 4.1× against the paper's
  11.5× and 28.1×), because the paper's random-mesh generator is not published (IC SPEC §10.2).
* Signpost refinement is about 1 % slower (instruction count) since it was made generic (IT N51).
* Not yet started:
  * integer coordinates: WP9 (Lean tier C, the whole-triangulation model), WP10 (the Thingi10k
    evaluation and the reference demo build), WP11 (consumers by default, Python, Blender "show
    refinement"; tangent frames I26 to decide);
  * FlipOut: F4 (Bézier, iCDT), F5/F6 (Lean), F7 (evaluation: Thingi10k, MPZ, Princeton; the
    reference driver exists in `runs_tmp/fo-reference-driver/`), F8 (Python, Blender), F9.
* Pending formula ids (untagged in the code until their proofs land): `ic-roundabout-update`,
  `ic-identify-edge`, `ic-trace-terminates`, `ic-strip-intersection`,
  `ic-extract-intrinsic-edge`, `ic-split-crossed-edge`, and the 13 `fo-` ids.
* `ic-vertex-removal-delta` (removal by flips in a Δ-complex): open in the paper and here; refused
  at run time with `RemovalStuck`.

## Housekeeping

* `ddg_core::distance`: the DOI and issue of Bornemann, Rasch 2006 were not confirmed
  independently (authors, title, journal, volume, pages, year were).
* Nothing is pushed beyond `origin/master` at 29ff52d (developability). `master` includes the
  repulsive-curves work since 2026-10-31.
