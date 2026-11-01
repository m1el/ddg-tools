# WP10 notes (Lean proofs and citations for Repulsive Surfaces)

Contract: `crates/ddg-repulsive-surfaces/SPEC.md` §12, §13.6. Namespace `DdgProofs.RS`; vectors,
`vnorm` and the motion calculus are reused from `DdgProofs.RC` (selective `open RC (...)`).

## Done

- `Basic.lean`: rs-exponents, rs-face-quantities (corner/opp/nvec/area/normal/bary/dualArea).
- `Geometry.lean`: rs-face-area-gradient, rs-face-normal-variation (along any motion).
- `Kernel.lean`: rs-kernel, rs-kernel-gradient (α > 1, valid at ν = 0 via `hasDerivAt_abs_rpow`).
- `Energy.lean`: rs-energy, rs-energy-scaling, rs-energy-gradient, rs-gradient-identities.

## Remaining

BarnesHut, Metric, Precond, Saddle, Constraints, Penalties, Step, Remesh, Ccd, Optimize; CITED.md
section; Audit.lean entries; check-formulas.

## How to resume

Worktree setup: symlink `lean/proofs/.lake/packages/*` to the main checkout's packages.
`cd lean/proofs && lake build DdgProofs.RepulsiveSurfaces.<Prev>` then
`lake env lean DdgProofs/RepulsiveSurfaces/<File>.lean > ../../runs_tmp/x.log 2>&1` per file.
