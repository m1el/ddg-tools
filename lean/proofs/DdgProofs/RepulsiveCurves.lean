import DdgProofs.RepulsiveCurves.Basic
import DdgProofs.RepulsiveCurves.Kernel
import DdgProofs.RepulsiveCurves.KernelGradient
import DdgProofs.RepulsiveCurves.Energy
import DdgProofs.RepulsiveCurves.Gram
import DdgProofs.RepulsiveCurves.Saddle
import DdgProofs.RepulsiveCurves.Constraints
import DdgProofs.RepulsiveCurves.Penalties
import DdgProofs.RepulsiveCurves.Schedule
import DdgProofs.RepulsiveCurves.Ccd
import DdgProofs.RepulsiveCurves.Hierarchy
import DdgProofs.RepulsiveCurves.Multigrid
import DdgProofs.RepulsiveCurves.BarnesHutAdjoint

/-!
# Repulsive Curves (Yu, Schumacher, Crane, ACM TOG 40(2), 2021)

Machine-checked derivations of the formulas of `crates/ddg-repulsive-curves` (contract
`crates/ddg-repulsive-curves/SPEC.md`, formula ids `rc-…`, §12). The smooth theory and the
convergence and error bounds are cited in `lean/proofs/CITED.md`, section "Repulsive Curves".
All theorems are in `namespace DdgProofs.RC`.
-/
