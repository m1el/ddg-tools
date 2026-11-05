import DdgProofs.Oracle.Scalar

/-!
# Oracle bridges: the intrinsic-triangulation formulas at `ℝ`

The oracle's (`lean/oracle`) scalar-generic definitions of the Delaunay test, instantiated at
`ℝ`, written in terms of the native quantities of `DdgProofs/Triangle.lean`: the excess
`θ_k + θ_l − π` (`delaunayExcess_eq`, also in the evaluation order of the Rust code,
`delaunayExcess_eq'`), the test with a tolerance (`delaunayHolds_iff`), the cotan weight as the
sum of two cotangents (`cotanWeight_eq`, `cotanWeightFromLengths_eq`), and the two opposite angles
of a diamond (`oppositeAngles_eq`).

App layer: these bridges stay with the oracle. (Library step S2b: moved here from
`DdgProofs/Intrinsic/Delaunay.lean`; step S3 adds the bridges of the other intrinsic formulas.)
-/

noncomputable section

namespace DDG

namespace Intrinsic

open scoped DDG.Oracle
open Real

theorem delaunayExcess_eq (θk θl : ℝ) : Ddg.delaunayExcess θk θl = θk + θl - π := by
  simp only [Ddg.delaunayExcess, Oracle.pi_eq]

theorem delaunayHolds_iff (tol θk θl : ℝ) :
    Ddg.delaunayHolds tol θk θl = true ↔ Ddg.delaunayExcess θk θl ≤ tol := by
  simp only [Ddg.delaunayHolds, decide_eq_true_eq]

theorem cotanWeight_eq (x y : ℝ) : Ddg.cotanWeight x y = x + y := rfl

theorem cotanWeightFromLengths_eq (a b c d e : ℝ) :
    Ddg.cotanWeightFromLengths a b c d e = cornerCot a b c + cornerCot a d e := by
  simp only [Ddg.cotanWeightFromLengths, Ddg.cotanWeight, Oracle.cornerCot_eq]

theorem oppositeAngles_eq {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    Ddg.oppositeAngles a b c d e = (cornerAngle a b c, cornerAngle a d e) := by
  simp only [Ddg.oppositeAngles, Oracle.cornerAngle_eq hk, Oracle.cornerAngle_eq hl]

/-- The order in which the Rust code evaluates the excess, `(θ_k − π) + θ_l`, is the same
number. -/
theorem delaunayExcess_eq' (θk θl : ℝ) : Ddg.delaunayExcess θk θl = (θk - π) + θl := by
  rw [delaunayExcess_eq]
  ring

end Intrinsic

end DDG
