import DDG.SDF.Basic
import DDG.SDF.Blend
import DDG.SDF.Circle2d
import DDG.SDF.Features
import DDG.SDF.MongeFit
import DDG.SDF.Operations
import DDG.SDF.Torus

/-!
# Points as Tori (Feng, Gkioulekas, Crane, ACM TOG 45(4), Article 53, 2026)

Machine-checked derivations of the formulas of `crates/ddg-points-as-tori` (contract
`crates/ddg-points-as-tori/SPEC.md`, formula ids `pat-…`, §13). The smooth theory (Laplace's
method, the signed Hopf–Cole transformation), the network and the baselines' approximation
properties are cited in `lean/proofs/CITED.md`, section "Points as Tori". All theorems are in
`namespace DDG.PAT`; vectors are `RC.V3`.
-/
