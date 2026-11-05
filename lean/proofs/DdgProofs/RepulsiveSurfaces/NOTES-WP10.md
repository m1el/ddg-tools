# WP10 notes (Lean proofs and citations for Repulsive Surfaces)

Contract: `crates/ddg-repulsive-surfaces/SPEC.md` §12, §13.6. Namespace `DDG.RS`; vectors,
`vnorm` and the motion calculus are reused from `DDG.RC` (selective `open RC (...)`).

## Done (all SPEC §12 ids)

Proved (35): every "proved" id of §12 plus `rs-fractional-preconditioner` (tagged in §6.3, missing
from §12):

- `Basic.lean`: rs-exponents, rs-face-quantities.
- `Geometry.lean`: rs-face-area-gradient, rs-face-normal-variation (along any motion).
- `Kernel.lean`: rs-kernel, rs-kernel-gradient (α > 1, valid at ν = 0 via `hasDerivAt_abs_rpow`).
- `Energy.lean`: rs-energy, rs-energy-scaling, rs-energy-gradient, rs-gradient-identities.
- `BarnesHut.lean`: rs-bh-energy, rs-bh-aggregate-adjoint, rs-box-distance.
- `Metric.lean`: rs-df-gradient (+ unit-edge counterexample), rs-sobolev-forms,
  rs-low-order-symmetrisation, rs-kernel-matrix-form, rs-hierarchical-null-space,
  rs-projector-average.
- `Precond.lean`: rs-preconditioner-order, rs-fractional-preconditioner. After D33 the id names the
  default augmented form `P = L′⁻¹ [[F, Cᵀ],[C, 0]] L′⁻¹` (`aug_precond_props`: equals
  `L′⁻¹ + L′⁻¹ blockdiag(F−K, 0) L′⁻¹`, symmetric, nonsingular iff the middle factor is,
  `C (P [r; b])_x = b` for all `r, b`, no extra hypothesis). The WP10 theorem for the printed form
  (`frac_precond_props`: symmetric, PSD, `C z_x = 0`) is kept, untagged (reference mode).
- `Saddle.lean`: rs-schur-direction, rs-schur-correction, rs-nested-saddle.
- `Constraints.lean`: rs-area-gradient, rs-volume-gradient (factor 3 on closed fans), rs-barycenter.
- `Penalties.lean`: rs-penalties, rs-willmore (structure; cotan/mass derivatives by FD).
- `Step.lean`: rs-armijo-reference. `Remesh.lean`: rs-circumcenter, rs-smoothing-weights,
  rs-delaunay-test. `Ccd.lean`: rs-ccd-advancement. `Optimize.lean`: rs-lbfgs-two-loop.

Cited (8, `CITED.md` § Repulsive Surfaces): rs-exponent-range, rs-differential-order,
rs-fractional-operator, rs-preconditioner-approx, rs-aqp, rs-bcqn, rs-energy-consistency,
rs-chen-holst. Nothing was moved from "proved" to "cited".

Audit: the tagged theorems are in `Audit.lean` (`#print axioms RS.…`); `AuditAll.lean` covers
every declaration automatically. `scripts/check-formulas.sh --audit-all` passes.

## SPEC errors found

- §6.2 / `rs-schur-correction`: `h = A′⁻¹C_NᵀS⁻¹Φ` with `S = −C_N A′⁻¹ C_Nᵀ` has `C_N h = −Φ`, so
  it must be **added** (paper Eq. 21), not subtracted (`schur_correction`, `schur_correction_sign`).
- §12 lacks the id `rs-fractional-preconditioner` used in §6.3 (proved here).

## How to resume

Worktree setup: symlink `lean/proofs/.lake/packages/*` to the main checkout's packages.
`cd lean/proofs && lake build DdgProofs.RepulsiveSurfaces` (log to `runs_tmp/`), or per file
`lake env lean DdgProofs/RepulsiveSurfaces/<File>.lean > ../../runs_tmp/x.log 2>&1` after building
its imports.
