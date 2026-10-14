# Working notes: vertex repositioning and iODT (checkpoint; fold into SPEC at the end)

Branch: worktree of `intrinsic-triangulations` (fast-forwarded to a0773fd). Build dir `target-iodt`.
`data/` holds symlinks to the main checkout's data. Logs: `runs_tmp/odt/` (study),
`runs_tmp/flow/` (cut flow study, `run_one.sh SUB TREAT`, jobs.txt).

## Done (committed)
* move_vertex, odt_target, star_min_shape, optimal_delaunay, shape_statistics; tests
  `tests/reposition.rs`, `tests/odt.rs`; bench `benches/retriangulation.rs`; example
  `examples/odt_study.rs`; SPEC N19-N23 (+ N24 draft); ddg-uv `IntrinsicOptions::optimal_delaunay`
  with state v3; flow study TREAT=odt, odt15.
* Lean proofs + oracle `move` task: delegated to a subagent (lean/ only).

## Findings to write up
* N24/N25: thin intermediate triangles (shape down to the N5 limit 1e-6) during iODT moves made
  signposts inaccurate (eps/shape ~ 1e-10); vertices hugging input edges then failed N8 traces on
  blub. Fix: shape floor of moves in iODT (`min_shape` 1e-2 or the star's current smallest);
  then snap 1e-12 suffices (blub mismatch 1.3e-11). Snap bound by shortest edge at the vertex.
* Study (first run, before the shape floor): iODT at matched counts: area cv and max/min much
  better, min angle worse (input small angles stay), 99% largest angle lower, max angle higher
  on camels (162-164 vs 130); kappa ~ iDT, below iDR; time 30-600x iDR.

## To do
- [ ] rerun study (running), rerun flow study odt/odt15 after final code, speed alone
- [ ] SPEC N24/N25 final, README section + table, ROADMAP, delete this file
- [ ] oracle replay in tests/oracle.rs (after the subagent reports the format)
- [ ] workspace tests, clippy, check-formulas, lake build
- [ ] remove examples/odt_debug.rs (scratch, not committed)
