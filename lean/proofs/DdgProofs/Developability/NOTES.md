# Developability proofs — progress notes

Contract: `crates/ddg-developability/SPEC.md` §9 (43 `dev-` ids: 30 proved, 13 cited).
Oracle definitions: `lean/oracle/DdgOracle/Formulas/Developability.lean` (namespace `Ddg.Dev`).
Statements of every id: `lean/proofs/README.md`, section "Developability of Triangle Meshes";
cited ids: `lean/proofs/CITED.md`, same heading.

Workflow: `cd lean/proofs && lake env lean DdgProofs/Developability/X.lean` checks one file
(the imports must be built). `scripts/check-formulas.sh --audit-all` passes (all 30 proved,
13 cited, 61 audit entries, only the standard axioms).

Status: **done.** All files compile without warnings, are imported from `DdgProofs.lean`, and
their main theorems are in `Audit.lean`.

| file | ids |
|---|---|
| Vec.lean | infrastructure: coordinates, `HasDerivAt3`, `atan2` (cos/sin, range, scaling, derivative) |
| Basic.lean | face-frame, corner-angle, area/normal/angle-gradient, backprop-face, vertex-normal |
| lean/ddg/DDG/Geometry/RigidInvariance.lean | gradient-invariance |
| Combinatorial.lean | star-bipartitions, combinatorial-width/-gradient/-max |
| Covariance.lean | covariance-matrix, eigenvalue-supergradient, covariance-gradient, angle-weight-retessellation |
| MaxCovariance.lean | spherical-circumcenter, max-covariance-gradient |
| Intrinsic.lean | log-map, log-map-gradient, tangent-eigen-derivative |
| Width.lean | intrinsic-width, width-altitude-gradient, min-width-gradient |
| Eig.lean | jacobi-rotation |
| Optimizer.lean | weak-wolfe-curvature, lewis-overton-line-search |
| Subdivision.lean | midpoint-subdivision |
| Remesh.lean | small-angle-classification |
| RegularCone.lean | regular-cone-values (E_λ and intrinsic only; E_P, λ_max, intrinsic max not proved) |

Open (possible later work): closed forms of SPEC §10.3 for E_P, λ_max, intrinsic max; the
"edge-connected = cyclic interval" equivalence (taken as the definition in Combinatorial.lean);
Rust `eig.rs` does not yet carry `[formula: dev-jacobi-rotation]`.
