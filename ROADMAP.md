# Roadmap and known gaps

State of 2026-10-14, branch `variational-cuts`. Numbers and evidence are in the documents this
file points to; this is the list of what is scheduled, what was decided, and what is known to be
open. Add to it; move an entry to the crate's README or SPEC when it is closed.

## Scheduled

1. **Intrinsic triangulations**: "Navigating Intrinsic Triangulations", Sharp, Soliman, Crane,
   ACM TOG 2019 (`../NavigatingIntrinsicTriangulations/paper.pdf`). The next paper. Its first use
   here: an intrinsic Delaunay triangulation for the stages of the cut flow that suffer on
   obtuse meshes (redistancing first, then the evaluation). Evidence:
   `crates/ddg-variational-cuts/README.md`, "Mesh quality", and `examples/bad_mesh_study.rs`.
   Decided for it (2026-10-14): the authors' MIT demo may be used as a reference; scope is the
   signpost data structure, flips and Delaunay flipping, vertex insertion and Delaunay refinement
   (with boundaries, which the paper does not cover), the overlay and point queries; the optimal
   Delaunay variant comes later; the whole cut flow is to run on the intrinsic triangulation, the
   cut is traced back, refinement that adds vertices is allowed; no visible effect in Blender; a
   detected wrong result (a failed trace) is an error, never repaired silently. Branches: `master`
   moves up to `variational-cuts` when item 1 lands, the work is on `intrinsic-triangulations`.

## Decided (2026-10-14)

* Two regions: the default stages are the ones that meet the exact laws (SPEC D12, D14–D16,
  D22); the paper's stages stay selectable (`FlowOptions::paper`).
* More than two regions: the recommended stages stay the paper's (`FlowOptions::for_regions`),
  although the evidence became mixed after D22 (SPEC §8.5, "After D22").
* No remeshing except inserting the cut (paper Sec. 6.1 is not implemented).

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
* Done before D22; not repeated since. The reference's Hencky value is about 900 times ours on
  the same partition (440 to 1010 across shapes), unexplained; the length-weight correspondence
  was calibrated through it on one shape. The reference runs were not converged (300 fixed
  steps). Only the Hencky energy was compared, only its cut in our energy, and not its
  flattening.
* Pictures were looked at for spot only (two regions; eight clusters; the comparison), all
  before D22.

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

## Housekeeping

* `ddg_core::distance`: the DOI and issue of Bornemann, Rasch 2006 were not confirmed
  independently (authors, title, journal, volume, pages, year were).
* Nothing is pushed; `master` points at `optimal-cones`.
