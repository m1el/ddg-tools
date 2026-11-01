import DdgProofs.IntrinsicSimplification.Flatten
import DdgProofs.IntrinsicSimplification.Monotone
import DdgProofs.IntrinsicSimplification.Curvature
import DdgProofs.IntrinsicSimplification.Interpolation
import DdgProofs.IntrinsicSimplification.Removal
import DdgProofs.IntrinsicSimplification.Wasserstein
import DdgProofs.IntrinsicSimplification.Transport
import DdgProofs.IntrinsicSimplification.Operators
import DdgProofs.IntrinsicSimplification.Oracle

/-!
# Surface Simplification using Intrinsic Error Metrics (Liu, Gillespie, Chislett, Sharp, Jacobson, Crane, ACM TOG 42(4), 2023)

Tiers A, A/B and B′ of `crates/ddg-intrinsic-simplification/SPEC.md` §13 (formula ids `iem-…`).
Cited results: `CITED.md`, section "Intrinsic Simplification". All theorems are in
`namespace DdgProofs.IntrinsicSimplification`.
-/
