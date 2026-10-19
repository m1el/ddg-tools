# Developability oracle: progress notes

- Done: Vec (V3, Jacobi), Mesh (validation, stars), Energy (six energies + subgradients), Step
  (first L-BFGS step with weak Wolfe, 4-1 subdivision), Format (`dev` command), Gen (`dev-gen`),
  Selftest (`dev-selftest`: closed forms of SPEC 10.3 to 1e-15, brute force over all subsets and a
  dense sphere search, finite differences, invariances, subdivision, line search).
- Independent of `Formulas/Developability.lean` and of the Rust crate (written from SPEC 1-5 only).
- Margins: E_P / pi_max / lambda_max / intrinsic max = second best minus best candidate value
  (pi_max and intrinsic max: also the gap to the runner-up pair); E_lambda / intrinsic =
  (l2 - l1) / tr; +inf for boundary stars with <= 2 faces.
- Todo: fixtures, Rust test, fuzz run.
