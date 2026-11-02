# WP7 notes (penalties, obstacles, implicit surfaces, Willmore, Coulomb)

Build/test: `CARGO_TARGET_DIR=$PWD/target-wp7 cargo test --release -p ddg-repulsive-surfaces --test penalties`
(0.3 s release). clippy `--all-targets` clean; `scripts/check-formulas.sh` OK.

## Done
- `src/implicit.rs`: re-exports the curves crate's `ImplicitSurface`, `Sphere`, `Torus`, `Plane`, `Union`;
  adds `Cylinder` (infinite) and `FramedTorus` (centre, axis, R, r); exact signed distances, zero gradient
  at singular points.
- `src/penalties/`: `Penalty` enum (every variant carries `weight`), `Quadrature`, `Penalty::validate`,
  `evaluate_one` (weighted value + weighted gradient), `evaluate` (the hook), `simple.rs` (area, volume,
  powers, squared error, boundary length/curvature), `willmore.rs` (`willmore_and_gradient`,
  `willmore_matrix`), `potentials.rs` (implicit obstacle/attractor, `coulomb_energy`).
- `src/obstacle.rs`: `StaticSet` (`from_obstacle`, `from_points`), `CrossEvaluator` (build / refit /
  energy / energy_and_differential / metric_diagonal), `obstacle_diagonal`, `ObstacleKind`.
- `tests/penalties.rs` (16 tests).

## Decisions
- Hook signature changed (energy/mod.rs, one call): `penalties::evaluate(pool, surface, positions,
  penalties, &exponents, eval, &mut differential)`. Tangent-point obstacles need the energy's exponents and
  theta (theta = 0 for `EnergyEval::Exact`, leaf size 8).
- `Penalty` no longer derives `PartialEq` (it holds `Arc<dyn ImplicitSurface>`); `Debug` prints name + weight.
- A penalty whose weighted value is not finite adds no gradient (implicit obstacle with d <= 0 gives +inf).
- Weights must be finite and >= 0; targets are explicit numbers (the session resolves "initial value").
- Quadrature `DualArea` includes the variation of `a_v` in the gradient.
- Mesh obstacle TP: lists built once (`CrossEvaluator::build`), frozen under `refit`; two-sided = mesh leaves vs
  static nodes + static leaves vs mesh face-tree nodes. Point cloud one-sided. Coulomb obstacle:
  `sum_S sum_T a_S a_T |X_S - X_T|^-p` between barycenters/centroids (not the curves crate's
  point-to-triangle distance), same BH lists.
- Willmore: `P = sum_{i interior} |(Wf)_i|^2 / m_i`; about 16 pi on the unit sphere (= 4 int H^2; the
  reference's value is 1/4 of it). `willmore_matrix(.., interior_only)`: true = Hessian-like term of the
  penalty (interior rows), false = full `W M^-1 W` (H^2 metric).
- Formula tags used: rs-penalties, rs-willmore, rs-bh-energy, rs-bh-aggregate-adjoint, rs-volume-gradient.
  Proposed new ids (not added; Lean side is not mine): rs-mesh-obstacle, rs-obstacle-diagonal,
  rs-coulomb-energy.

## Remaining
Nothing planned. The session (WP5) should use `CrossEvaluator` directly in the line search (frozen lists);
`evaluate` rebuilds trees per call.
