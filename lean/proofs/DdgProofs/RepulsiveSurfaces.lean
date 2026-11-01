import DdgProofs.RepulsiveSurfaces.Basic
import DdgProofs.RepulsiveSurfaces.Geometry
import DdgProofs.RepulsiveSurfaces.Kernel
import DdgProofs.RepulsiveSurfaces.Energy
import DdgProofs.RepulsiveSurfaces.Metric
import DdgProofs.RepulsiveSurfaces.Saddle
import DdgProofs.RepulsiveSurfaces.Constraints
import DdgProofs.RepulsiveSurfaces.BarnesHut
import DdgProofs.RepulsiveSurfaces.Precond
import DdgProofs.RepulsiveSurfaces.Step
import DdgProofs.RepulsiveSurfaces.Remesh
import DdgProofs.RepulsiveSurfaces.Ccd
import DdgProofs.RepulsiveSurfaces.Optimize

/-!
# Repulsive Surfaces (Yu, Brakensiek, Schumacher, Crane, ACM TOG 40(6), 2021)

Machine-checked derivations of the formulas of `crates/ddg-repulsive-surfaces` (contract
`crates/ddg-repulsive-surfaces/SPEC.md`, formula ids `rs-…`, §12). The smooth theory, the
approximation properties and the convergence of the optimizers are cited in
`lean/proofs/CITED.md`, section "Repulsive Surfaces". All theorems are in
`namespace DdgProofs.RS`; vectors, norms and the curve-motion calculus come from the curves
package (`DdgProofs.RC`).
-/
