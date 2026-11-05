import DDG.Geometry
import DDG.Optimization
import DDG.Mesh
import DDG.Laplace
import DDG.Intrinsic
import DDG.VectorFields
import DDG.FEM
import DDG.NormalCoords
import DDG.Delaunay
import DDG.Conformal
import DDG.OptimalCones
import DDG.Repulsive
import DDG.SDF

/-!
# DDG: a formally verified library of discrete differential geometry

The umbrella library: it imports the fourteen part roots, each of which is a Lake library of its own
(`DDGGeometry`, `DDGOptimization`, `DDGMesh`, `DDGLaplace`, `DDGIntrinsic`, `DDGVectorFields`,
`DDGFEM`, `DDGNormalCoords`, `DDGDelaunayPlanar`, `DDGDelaunay`, `DDGConformal`, `DDGOptimalCones`,
`DDGRepulsive`, `DDGSDF`) that can be built without the others it does not depend on. The part
graph is in `parts.toml` and is checked by `scripts/check-lean-dag.py`.

Direct dependencies of the parts (Mathlib below all of them): `Mesh` on `Geometry`; `Laplace` on
`Mesh`; `Intrinsic`, `Repulsive`, `SDF` on `Geometry`; `VectorFields` on `Mesh`; `FEM` on `Laplace`;
`NormalCoords` on `Intrinsic`; `Delaunay.Planar` on `Geometry`, `FEM`; `Delaunay` on `Delaunay.Planar`,
`Intrinsic`, `FEM`, `Laplace`, `Mesh`; `Conformal` on `FEM`, `Laplace`, `Mesh`; `OptimalCones` on
`Optimization`, `Conformal`. `Optimization` has no dependencies besides Mathlib.

Conventions: all proofs are complete and the library adds no axioms: every declaration depends on at most
`propext`, `Classical.choice` and `Quot.sound`; `autoImplicit = false`; definitions are native over `ℝ`/`ℤ`
(the library does not mention the Float oracle). Theorems carry docstring tags (`formula:` followed by a stable id) that
link them to formulas of the papers and of the `ddg-tools` Rust code. See `README.md` for the part
table, references and build instructions.
-/
