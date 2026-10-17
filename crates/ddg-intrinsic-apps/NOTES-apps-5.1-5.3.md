# Notes: applications 5.1–5.3 (checkpoint)

Working copy: a clone of ddg-tools inside the agent worktree (the isolation worktree was made of
the demo repo, not ddg-tools):
`/var/data/discrete-geometry/NavigatingIntrinsicTriangulations/navigating-intrinsic-triangulations-demo/.claude/worktrees/agent-a825df4ddbb500d57/ddg-tools`,
branch `apps-5123`, based on `75b6e2c`. Merge with `git fetch <that path> apps-5123`.

## Done

- Crate `crates/ddg-intrinsic-apps` (lib.rs: one `pub mod` line per module): `exact.rs` (exact
  polyhedral distance, own window propagation with Xin–Wang vertex filter), `steiner.rs` (5.1),
  `fem.rs` (5.2: stiffness, Galerkin/lumped mass, 4-1 subdivision, Lanczos condition number,
  `PinnedLaplacian`), `heat.rs` (5.3: `HeatSolver`, `geodesic_distance`, `value_at`).
- Tests (all pass, release): exact 6, heat 6, fem 6, steiner 5. Bench `benches/heat.rs` builds
  (not yet run alone). Example `examples/apps_study.rs` (steiner|fem|heat per mesh); logs of the
  runs in `runs_tmp/study-*.log` of the clone; tables are in README.md.
- ddg-py: `geodesic_distance`, `steiner_tree` (`src/intrinsic_apps.rs`,
  `src/python/intrinsic_apps.rs`), stubs, `python-tests/test_intrinsic_apps.py`, py-test.sh loop
  entry. `scripts/py-test.sh`: PASS x3 (44, 34, 5 tests).
- Lean: `lean/proofs/DdgProofs/IntrinsicApps/HeatFem.lean` (app-heat-gradient,
  app-heat-divergence, app-fem-galerkin-mass, app-fem-midpoint-subdivision), checked alone with
  `lake env lean` (OK); import in DdgProofs.lean, Audit entries; CITED.md entry
  app-steiner-fermat-point. `scripts/check-formulas.sh`: OK.
- Cited papers downloaded to `../NavigatingIntrinsicTriangulations/cited/applications-5.1-5.3/`.
- README.md first full draft (with the three tables).

## Remaining

1. `check-formulas.sh --build` FAILED in `DdgProofs/Intrinsic/Move.lean` (unknown
   `Ddg.movedLength`): the copied `.lake` of the oracle (lean/oracle/.lake, copied from the main
   checkout) is probably stale; not caused by HeatFem.lean. Next: `cd lean/oracle && lake build`,
   then rerun `scripts/check-formulas.sh --build`. `.lake` dirs in the clone are copies (8 GB),
   untracked.
2. README corrections: in 5.2 text, iODT is lower than iDR in 11 of 14 matched comparisons
   (iDR lower: bob and bump at two levels, sheared torus at one level); the sheared square at one
   level has iODT at 865 vertices (bisection missed 1089). In 5.3 "well-shaped meshes ... 1 to
   2 %" should read 1 to 2.3 %; the m=10/m=0.1 remark was measured at the far corner only.
3. Write `crates/ddg-intrinsic-apps/SPEC.md` (sections 5.1–5.3, decisions N30–N39):
   N30 Fermat barycentrics `a_k csc(θ_k+π/3)` (paper prints trilinears); N31 SLL ratio priority,
   MST-adjacent candidates, boundary loops pre-joined, no reflip, genus>0 not a disk; N32 refused
   insertion drops triangle and reruns; N33 cut_to_disk definition; N34 FEM matrices and
   condition number (Lanczos, pinned shift-invert, three masses); N35 heat details (lumped M,
   t=m h² of the triangulation, δ rhs, zero X on flat faces counted, pinned Poisson instead of
   gc's 1e-6 shift, gc's shift rule, BCs, Dirichlet refuses boundary sources); N36 value_at /
   inserted vertex places; N37 exact ground truth (own algorithm, flat vertices as pseudo-sources,
   thin-window drop 1e-12, validation); N38 metrics, iDR r/2 (angle bound alone does not reduce
   h), boundary-sliver finding on the sheared square; N39 Python scope, Blender nothing (natural
   use: 5.1 cut as seams needs insertion).
4. `crates/ddg-intrinsic/SPEC.md` §2: one line "the applications now live in ddg-intrinsic-apps".
   `ROADMAP.md`: landed + open (extrinsic ODT for Fig. 17, genus>0 cut graphs, seams in Blender).
   ddg-py README: short section on the two functions.
5. Run alone: `cargo bench -p ddg-intrinsic-apps --bench heat`; `cargo test --workspace
   --release` (log); `cargo clippy -p ddg-intrinsic-apps -p ddg-py --all-targets`.

## Commands (from the clone; the sandbox guard needs scripts, not `cd … && cargo`)

```text
export CARGO_TARGET_DIR=$PWD/target-apps1
cargo test --release -p ddg-intrinsic-apps > runs_tmp/t-apps.log 2>&1
cargo build --release -p ddg-intrinsic-apps --examples
target-apps1/release/examples/apps_study heat data/models/spot/spot/spot_triangulated.obj
CARGO_TARGET_DIR=$PWD/target-apps1-py scripts/py-test.sh
scripts/check-formulas.sh [--build]
```
