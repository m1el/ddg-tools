import DdgProofs.PointsAsTori.Basic
import DdgProofs.PointsAsTori.Features
import DdgProofs.PointsAsTori.Fitting
import DdgProofs.PointsAsTori.Torus
import DdgProofs.PointsAsTori.Blend
import DdgProofs.PointsAsTori.Apps

/-!
# Points as Tori (Feng, Gkioulekas, Crane, ACM TOG 45(4), Article 53, 2026)

Machine-checked derivations of the formulas of `crates/ddg-points-as-tori` (contract
`crates/ddg-points-as-tori/SPEC.md`, formula ids `pat-…`, §13). The smooth theory (Laplace's
method, the signed Hopf–Cole transformation), the network and the baselines' approximation
properties are cited in `lean/proofs/CITED.md`, section "Points as Tori". All theorems are in
`namespace DdgProofs.PAT`; vectors are `RC.V3`.
-/
