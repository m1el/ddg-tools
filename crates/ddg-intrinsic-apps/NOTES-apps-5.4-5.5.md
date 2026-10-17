# Checkpoint: paper Sec. 5.4 (AMR) and 5.5 (tangent fields)

Worktree: the clone `ddg-tools/` inside
`.../navigating-intrinsic-triangulations-demo/.claude/worktrees/agent-aed156ad98b2d7138/`
(the isolation worktree was of the demo repo, so ddg-tools was cloned into it), branch
`apps-amr-vector`, based on 75b6e2c. `data` there is a symlink to the main checkout's data/.

## Done (committed)
* `ddg-intrinsic`: `refine_with(options, &dyn Fn(&tri, face) -> bool)` (extra per-face
  criterion); `refine` = `refine_with` with a false criterion (behaviour unchanged).
* `ddg-intrinsic-apps` (`amr.rs`, `vector.rs`, `vector/{complex,connection,knoeppel,eigen,linsolve}.rs`):
  - vector: transport from signposts (boundary frames unnormalized, s=1), FD connection
    Laplacian (L_ij = −w r_ji), Knöppel FE connection Laplacian (own derivation, exact
    monomial integration, unwrapped face holonomy), block inverse iteration smoothest field,
    singularities (index sum = nχ), harmonic field with boundary values (tangent/normal/given,
    interpolated along boundary), vector maximum principle check, 3D output at input vertices.
  - amr: Green's function, heat kernel (one backward Euler step), harmonic disk map; residual
    estimator (h_T = diameter, not A² as the supplement prints), bulk marking (supplement sums or
    Dörfler squares), `refine_marked`, `adaptive_solve`, `l2_error` (point queries),
    `input_vertex_error`, `map_quality`.
  - tests: tests/vector.rs (7), tests/amr.rs (7), unit tests (6): all pass
    (`cargo test --release -p ddg-intrinsic-apps`).
  - examples: amr_study (AMR vs uniform at equal L² error, reference = adaptive + uniform),
    vector_study (smoothest | harmonic).
* Lean: `lean/proofs/DdgProofs/Intrinsic/AppsAmrVector.lean` (holonomy around a face, FD edge
  energy, singular values from conformal parts); compiles alone (`lake env lean`). CITED.md
  entries app-*; DdgProofs.lean import, Audit.lean lines. `scripts/check-formulas.sh` passes.

## WIP (committed, not built)
* ddg-py: `src/intrinsic_vector.rs` (pure Rust, with unit tests), `src/python/intrinsic_vector.rs`
  (glue), Cargo.toml deps, `pub mod intrinsic_vector` in lib.rs. NOT yet: `mod intrinsic_vector;`
  and `intrinsic_vector::register(m)?;` in `src/python.rs`, stubs in `ddg_tools.pyi`, a python
  test file and its entry in `scripts/py-test.sh`. Not compiled yet.

## Measured so far (runs_tmp/s3, v1, vh1; parallel runs, timings not alone)
* Heat kernel, spot control mesh (188 v): AMR 1740 v / 0.039 s reaches L² error 7.1e-3; uniform
  needs ~153k v / 2.2 s (88× vertices, 55× time). bob control: 50–60× vertices.
* Green's function, spot: 2–10× fewer vertices, time 0.5–7× (AMR wins below ~7e-3 error).
* Disk map: input of jittered sheared grid has 1 flipped triangle, iDT and AMR none; AMR
  1.3–1.6× fewer vertices, not faster.
* Smoothest fields (λ·area, input → iDT → iDR): modest except camelshotnoise n=1:
  60.9 → 58.6 → 34.6. Round sphere: λ ≈ n as in the smooth case.
* Max principle (harmonic field, FD): sheared grids with shear ≥4: 10–17 interior vertices
  outside the cone on the input, 0 on iDT and iDR.

## Next steps
1. Finish ddg-py registration in python.rs, stubs, `python-tests/test_intrinsic_vector.py`,
   add it to the loop in `scripts/py-test.sh`; run `scripts/py-test.sh`.
2. Lean build: `scripts/check-formulas.sh --build` failed in DdgProofs/Intrinsic/Move.lean
   (unknown `Ddg.movedLength`): the copied `lean/oracle/.lake` is stale; rebuild the oracle
   (`cd lean/oracle && lake build`) then rerun. My file itself compiles.
3. Docs: crate README.md sections 5.4/5.5 (tables above, rerun timings alone), SPEC.md N40–N49
   (decisions listed in the module docs: N40 frames/storage, N41 boundary s=1, N42 FD sign,
   N43 FE model + unwrapped Ω, N44 eigen solver, N45 max principle / boundary values, N46 AMR
   problems, N47 estimator h_T and point sources, N48 marking + refine_with hook, N49
   measurement protocol), ROADMAP.md.
4. `cargo test --workspace --release`, clippy on touched crates.

## Commands
```
S=/tmp/claude-1000/-var-data-discrete-geometry/938f789e-bf77-4d0d-ba4f-bd1ba09d3694/scratchpad/amrvec
$S/run.sh t cargo test --release -p ddg-intrinsic-apps
$S/run.sh b cargo build --release -p ddg-intrinsic-apps --examples
$S/study.sh s4 $S/jobs3.txt          # AMR studies in parallel (logs runs_tmp/s4/)
$S/par.sh v2 $S/jobs_v1.txt          # smoothest-field study
target-apps2/release/examples/vector_study harmonic sheared:12:4:0.45:2
```
