# Developability proofs — progress notes

Contract: `crates/ddg-developability/SPEC.md` §9 (43 `dev-` ids: 30 proved, 13 cited).
Oracle definitions: `lean/oracle/DdgOracle/Formulas/Developability.lean` (namespace `Ddg.Dev`).

Workflow: `cd lean/proofs && lake env lean DdgProofs/Developability/X.lean` checks one file
(needs the imports built: `lake build DdgProofs.Developability.Vec` etc.).

| file | ids | status |
|---|---|---|
| Vec.lean | (infrastructure: coords, `HasDerivAt3`, `atan2`) | done |
| Basic.lean | face-frame, corner-angle, area/normal/angle-gradient, backprop-face, vertex-normal | todo |
| Invariance.lean | gradient-invariance | todo |
| Combinatorial.lean | star-bipartitions, combinatorial-width/-gradient/-max | todo |
| Covariance.lean | covariance-matrix, eigenvalue-supergradient, covariance-gradient, angle-weight-retessellation | todo |
| MaxCovariance.lean | spherical-circumcenter, max-covariance-gradient | todo |
| Intrinsic.lean | log-map, log-map-gradient, tangent-eigen-derivative | todo |
| Width.lean | intrinsic-width, width-altitude-gradient, min-width-gradient | todo |
| Eig.lean | jacobi-rotation | todo |
| Optimizer.lean | weak-wolfe-curvature, lewis-overton-line-search | todo |
| Subdivision.lean | midpoint-subdivision | todo |
| Remesh.lean | small-angle-classification | todo |
| RegularCone.lean | regular-cone-values | todo (may drop: SPEC allows) |
