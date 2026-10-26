# ddg-tools

Rust implementations of discrete differential geometry papers, with Blender add-ons built on them.
The papers (PDFs) and third-party reference code live next to this repo in `/var/data/discrete-geometry/`.
Treat everything outside this repo as read-only reference.

First paper: "Optimal Cone Singularities for Conformal Flattening", Soliman, Slepčev, Crane, ACM TOG 2018
(`../OptimalCones/paper.pdf`, thesis `../OptimalCones/YousufSolimanMSThesis.pdf`).
Its cones are turned into UV maps with "Boundary First Flattening", Sawhney & Crane, ACM TOG 2017
(reference C++, MIT: `../ddg-course/boundary-first-flattening`).

Second paper: "Variational Surface Cutting", Sharp, Crane, ACM TOG 2018 (`../VariationalCuts/paper.pdf`,
author's derivation notes `../VariationalCuts/cea_tutorial.pdf`, `cea_derivations.pdf`).
**Clean room:** it is implemented from the text only. The authors' code is GPL; never open, search,
download or paraphrase it (it is not in this tree; do not look for it, and do not fetch the
`variational-surface-cutting` repository from the web). Papers that the text cites may be downloaded and cited.
The contract for this paper is `crates/ddg-variational-cuts/SPEC.md`.

Third paper: "Navigating Intrinsic Triangulations", Sharp, Soliman, Crane, ACM TOG 2019
(`../NavigatingIntrinsicTriangulations/paper.pdf`; reference: the authors' demo in the same directory, MIT,
which may be read and ported with its licence notice). Contract: `crates/ddg-intrinsic/SPEC.md`.
`ROADMAP.md` lists what is scheduled and the known gaps.

Fourth paper: "Developability of Triangle Meshes", Stein, Grinspun, Crane, ACM TOG 2018
(`../DiscreteDevelopable/paper.pdf`). **Clean room**, as for the second paper: the authors' code (GPL,
with Kopp's LGPL 3x3 eigensolver) is quarantined outside this tree; never open, search, download or
paraphrase it, and do not copy Kopp's published C code (his paper may be read and cited). Contract:
`crates/ddg-developability/SPEC.md`; formula ids start with `dev-`.

## Layout

```
crates/ddg-core            mesh topology + validation, geometry, cotan operators,
                           sparse linear algebra wrappers, OBJ I/O, test-mesh generators
crates/ddg-optimal-cones   the optimal-cones solver (the paper)
crates/ddg-bff             boundary first flattening + cutting
crates/ddg-variational-cuts  the cut flow of the second paper; `SPEC.md` is its contract
crates/ddg-intrinsic       intrinsic triangulations (signposts, flips, Delaunay, refinement, overlay);
                           `SPEC.md` is its contract; the authors' MIT demo may be used as a reference
crates/ddg-intrinsic-apps  applications of intrinsic triangulations (paper Sec. 5): Steiner trees, FEM,
                           heat-method distance, adaptive refinement, tangent fields; `SPEC.md` N30-N49
crates/ddg-developability  the developability flow of the fourth paper (energies, L-BFGS, subdivision,
                           remeshing, cut and flatten); `SPEC.md` is its contract
crates/ddg-uv              end-to-end pipeline: polygons -> charts -> cones -> cut -> flatten -> per-corner UVs;
                           `cutting`: the cut flow on polygon meshes with seams, insertion of the cut, unwrap
crates/ddg-py              PyO3 bindings (Python module `ddg_tools`)
blender/                   Blender extension: mesh -> UV map; cut operators (Start, Flow, Apply); map to sphere;
                           Make Developable, Show Developability, Developable Seams + UV
lean/oracle/               Lake package, no Mathlib: scalar-generic formula definitions and a Float
                           executable used as a differential-testing oracle
lean/proofs/               Lake package on Mathlib: machine-checked derivations of the formulas the Rust code uses
scripts/                   helper scripts
data/                      gitignored; `data/optimal-cones/optimal-cones-data/data/*.obj` + `*-cones.txt`;
                           `data/models/` for meshes unpacked from `/var/data/discrete-geometry/models/*.zip`
                           (Keenan Crane's model repository: spot, bob, blub, ...)
runs_tmp/                  gitignored; logs of long-running commands
```

## Ground rules

- Rust stable (see `rust-toolchain.toml`), edition 2024, MIT license. `f64` everywhere; stored indices are `u32`.
- Correctness first, then performance. Performance floor for the cones solver: at least match the paper
  (100–150k triangles in 20–25 s on a 2018 laptop); we should be far below that. The cut flow has no
  published timings: be performance aware (reuse factorizations, no quadratic loops), measure later.
- APIs are explicit: plain option structs, no hidden global state, no silent normalization or clamping.
  Every solve returns a report (status, iteration counts, residuals).
- Bad input never panics. Return error enums whose messages say what is wrong and where
  (e.g. "edge (12, 57) is shared by 3 faces: 4, 9, 31").
- Non-manifold meshes are rejected. Polygons are triangulated for the algorithm only; the caller's
  polygons are never changed, and results are reported per original polygon corner.
- Credit every paper: each crate's README and crate-level rustdoc cite the paper it implements;
  ported code keeps the original license notice.
- Python goes through `uv` (`uv run`, `uvx`), never `pip`.
- This machine has 64 cores. Independent jobs (one run per mesh, parameter sweeps, studies, fixture
  generation) run in parallel across them (`xargs -P`, a thread pool), one log per job; never one
  after the other. Exception: timings and benches that report single-thread numbers run alone.
- Long-running commands write to a log under `runs_tmp/` (`cmd > runs_tmp/x.log 2>&1`), then read the log.
  Never pipe a slow command into `tail`/`head`. Wait on a PID with `kill -0 $PID`, not `pgrep -f`.

## Several agents work in this checkout at the same time

- Stay inside the paths you were assigned. Do not edit another crate to make yours compile; report it instead.
- Use your own build directory so builds do not block each other:
  `export CARGO_TARGET_DIR=$PWD/target-<task>` and build/test with `-p <crate>`.
- Commit often, in small commits, on the current branch (whatever is checked out; `variational-cuts`
  for the second paper), and only your own paths:
  `git add <paths> && git commit -m "<msg>" -- <paths>`.
  Never `git add -A`, `git commit -a`, stash, reset, rebase, clean, or switch branches.
  If `.git/index.lock` exists, wait a second and retry.
- End every commit message with the `Co-Authored-By` line given in your own session instructions, so that it
  names the model that actually made the commit.
- Shared files (`Cargo.toml`, `Cargo.lock`, `README.md`, `CLAUDE.md`, `.gitignore`): re-read right before
  editing and keep edits minimal. Commit `Cargo.lock` together with the change that caused it.

## Math conventions (the contract between all crates, the oracle and the proofs)

The conventions on `L`, the masses, `u` and the angle defects hold for every paper. The canonical
problem below is the first paper's; the second paper's contract is `crates/ddg-variational-cuts/SPEC.md`
(same `L`, same sign of `u`; its formula ids start with `cut-`).

Mesh: manifold oriented triangle mesh, vertices `V`, boundary vertices `B` (may be empty), interior `I = V \ B`.

- Corner angles come from edge lengths. Angle defect `Ω_i = 2π − Σ(corner angles at i)` for `i ∈ I`;
  boundary geodesic curvature `k_i = π − Σ(corner angles at i)` for `i ∈ B`.
  Discrete Gauss–Bonnet: `Σ_I Ω + Σ_B k = 2πχ`.
- `L` is the cotan Laplacian, symmetric **positive semidefinite**:
  `L_ij = −½ Σ cot(angle opposite edge ij)` over the one or two triangles on edge `ij`, `L_ii = −Σ_{j≠i} L_ij`.
- `m_i = ⅓ Σ(area of triangles at i)` (lumped mass). **Area normalization:** masses are divided by the
  total area so `Σ m_i = 1`. `L` and `Ω` are scale invariant. `λ` is dimensionless under this
  normalization (paper Sec. 6.2: default ½, useful range 0.1–1).
- `u` is the **true** log conformal scale factor: flattened edge length `≈ e^{(u_i+u_j)/2} ℓ_ij`.
  The linearized discrete Yamabe equation is `(L u)_i = µ_i − Ω_i` for `i ∈ I`, with `u_i = 0` on `B`.
  The paper prints `L u = Ω − µ` (Eq. 17), which with its own PSD `L` is `−u`; the cones are identical.
  BFF uses our sign (`A u = −(K − C)`).
- A **cone angle** is the target curvature `µ_i = 2π − (angle sum around i in the flattened cone metric)`.
  Positive means a cone tip. This matches the dataset's `*-cones.txt`. On a closed mesh `Σ µ = 2πχ`.

**Canonical discrete problem.** Find `µ ∈ R^I` (and the `u` it determines) minimizing

```
½ Σ_i wE_i m_i u_i²  +  λ Σ_i wR_i |µ_i|
subject to   (L u)_i = µ_i − Ω_i  (i ∈ I),   u_i = 0  (i ∈ B),   lo_i ≤ µ_i ≤ hi_i
```

- `wE_i ≥ 0` weights distortion (default 1; 0 excludes the vertex from the distortion measure).
- `wR_i ∈ [0, +∞]` weights cone cost (default 1; 0 is a free cone site; `+∞` forbids a cone).
- `lo_i ≤ 0 ≤ hi_i` bound the cone angle (default unbounded).
- Prescribed cones fix `µ_i` to a given angle at no cost.
- Closed mesh (`B = ∅`): `u` is defined up to an additive constant, which is optimized too.

This problem statement is the ground truth. The paper's algorithm (regularized stage, then semismooth
Newton, Fig. 12) is one way to solve it; any correct solver returns the same `u` and `µ`
(unique when all `wE_i > 0`). Where the paper is silent or inconsistent (closed meshes, angle bounds,
signs, the mass factor in the regularized stage), derive the answer from this problem and say so in the docs.

## Formula traceability

Every Rust function that implements a mathematical formula carries a doc line `[formula: <kebab-id>]`.
The matching Lean theorem in `lean/proofs` carries the same tag in its docstring. Results that are
cited rather than machine-checked (the paper's smooth theory) are listed in `lean/proofs/CITED.md`
with the same id, an exact reference, and the reason. `scripts/check-formulas.sh` cross-checks the tags.

## Verification layers

1. Optimality certificate: every solve checks its own optimality conditions and duality gap.
2. Generic convex solver on the canonical problem for small meshes (independent algorithm).
3. Lean oracle (`lean/oracle`): independent Float implementation for differential fuzzing, plus committed golden fixtures.
4. Lean proofs (`lean/proofs`): machine-checked derivations of the formulas.

For the cut flow there is no unique minimizer to certify. Its layers are: analytic solutions (sphere
cap, flat disk), finite differences of the discrete energy against the predicted gradient, energy
decrease and invariances, the Lean oracle for one flow step, and the Lean proofs (SPEC.md §6, §10).
