import DDG.Repulsive.Curves.BarnesHutAdjoint
import DDG.Repulsive.Curves.Constraints
import DDG.Repulsive.Curves.DualLengths
import DDG.Repulsive.Curves.Energy
import DDG.Repulsive.Curves.Gram
import DDG.Repulsive.Curves.Hierarchy
import DDG.Repulsive.Curves.Kernel
import DDG.Repulsive.Curves.KernelGradient
import DDG.Repulsive.Curves.Multigrid
import DDG.Repulsive.Curves.Penalties
import DDG.Repulsive.Curves.Saddle
import DDG.Repulsive.Surfaces.BarnesHut
import DDG.Repulsive.Surfaces.Basic
import DDG.Repulsive.Surfaces.Constraints
import DDG.Repulsive.Surfaces.Energy
import DDG.Repulsive.Surfaces.Geometry
import DDG.Repulsive.Surfaces.Kernel
import DDG.Repulsive.Surfaces.Metric
import DDG.Repulsive.Surfaces.Penalties
import DDG.Repulsive.Surfaces.Precond
import DDG.Repulsive.Surfaces.Saddle

/-!
# Repulsive Curves (Yu, Schumacher, Crane, ACM TOG 40(2), 2021)

Machine-checked derivations of the formulas of `crates/ddg-repulsive-curves` (contract
`crates/ddg-repulsive-curves/SPEC.md`, formula ids `rc-…`, §12). The smooth theory and the
convergence and error bounds are cited in `lean/proofs/CITED.md`, section "Repulsive Curves".
All theorems are in `namespace DDG.RC`.
-/

/-!
# Repulsive Surfaces (Yu, Brakensiek, Schumacher, Crane, ACM TOG 40(6), 2021)

Machine-checked derivations of the formulas of `crates/ddg-repulsive-surfaces` (contract
`crates/ddg-repulsive-surfaces/SPEC.md`, formula ids `rs-…`, §12). The smooth theory, the
approximation properties and the convergence of the optimizers are cited in
`lean/proofs/CITED.md`, section "Repulsive Surfaces". All theorems are in
`namespace DDG.RS`; vectors, norms and the curve-motion calculus come from the curves
package (`DDG.RC`).
-/
