# Developability oracle: progress notes

Done
- `Vec` (V3, cyclic Jacobi), `Mesh` (validation incl. duplicate faces, stars, face data), `Energy` (six
  energies + subgradients), `Step` (first L-BFGS step with weak Wolfe, 4-1 subdivision), `Format`
  (`dev`), `Gen` (`dev-gen`), `Selftest` (`dev-selftest`), `Fixtures` (`fixtures/dev`, 48 cases).
- Independent of `Formulas/Developability.lean` and of the Rust crate (written from SPEC 1-5 only).
- Rust: `crates/ddg-developability/tests/oracle.rs` (fixtures always, fuzz `#[ignore]`).

Decisions
- Margins: E_P / pi_max / lambda_max / intrinsic max = second best minus best value (pi_max and
  intrinsic max also the gap to the runner-up pair); E_lambda / intrinsic = (l2 - l1)/tr. +inf for
  boundary stars with one face. Boundary stars with at most two faces contribute exactly nothing (SPEC D33, margin +inf).
  (Before D33 the oracle reported the real eigenvalue gap there and the Rust test skipped them.)
- The generator makes hinge stars with at least one face on each side (a side without ring points
  makes a degenerate face through the apex).
- Step test skips: ties anywhere, E0 < 1e-8, mean_edge_length * |g|_inf < 1e-6, or a Wolfe decision
  closer than 1e-7 (relative to the first order change); trial counts may differ by one for capped
  steps (cap knife edge, last bit of t_max).

- SPEC 4.9 (floor): step input line `floor` (default 1e-4). Start check first (`near-degenerate-face`,
  relative slack 1e-6 on sin), trials with a face below the floor are +inf (counted `inadmissible`, no
  energy evaluation), blocked vertices (face min angle < 2 floor at the START point) have zero gradient;
  the blocked set stays fixed during the line search (recomputing it at the trial point disagreed with
  Rust on the dev-step-blocked fixture: the trial count differed).

Findings so far
- Debugging a 30 % gradient disagreement on near-hinge stars (E ~ 1e-18) showed that the adjoint formula
  of SPEC 2/3.4 evaluated in float64 is limited by the arbitrary eigenvector of 2-face boundary stars
  with nearly coplanar faces (spurious y = 2 theta (x.N) x ~ normal angle); both implementations have it.
