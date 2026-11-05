import DdgProofs.Oracle.Geometry
import DDG.Intrinsic.Barycentric
import DDG.Intrinsic.Signpost
import DDG.Intrinsic.Trace
import DdgProofs.Delaunay.Local
import DdgProofs.Delaunay.ODT

/-!
# Oracle bridges: the intrinsic-triangulation formulas at `ℝ`

The oracle's (`lean/oracle`, `DdgOracle/Formulas/Intrinsic.lean`) scalar-generic definitions of
the intrinsic-triangulation formulas, instantiated at `ℝ`, are the native definitions of the
library (`DdgProofs/Intrinsic/*`, `DdgProofs/Delaunay/Local.lean`, `DdgProofs/Delaunay/ODT.lean`,
namespace `DDG.Intrinsic`). One-way bridges `Ddg.foo x = DDG.Intrinsic.foo x`: the theorems of the
library apply to the oracle's numbers after rewriting with them. Only the oracle's corner angle
needs a hypothesis (a nondegenerate triangle, `oppositeAngles_eq`), because its half-angle form
with `atan` agrees with the law of cosines only there.

App layer. (Library step S2b: the Delaunay bridges moved here from
`DdgProofs/Intrinsic/Delaunay.lean`; step S3: all bridges, in namespace `DDG.Oracle`.)
-/

noncomputable section

namespace DDG

namespace Oracle

open Real

/-! ## Unfolding and flipping a diamond -/

theorem diamondX_eq (a b c : ℝ) : Ddg.diamondX a b c = Intrinsic.diamondX a b c := by
  simp only [Ddg.diamondX, Intrinsic.diamondX, two_eq]

theorem diamondH_eq (a b c : ℝ) : Ddg.diamondH a b c = Intrinsic.diamondH a b c := by
  simp only [Ddg.diamondH, Intrinsic.diamondH, two_eq, triArea_eq]

theorem flippedLength_eq (xk hk xl hl : ℝ) :
    Ddg.flippedLength xk hk xl hl = Intrinsic.flippedLength xk hk xl hl := by
  simp only [Ddg.flippedLength, Intrinsic.flippedLength, sqrt_eq]

theorem flippedLengthFromLengths_eq (a b c d e : ℝ) :
    Ddg.flippedLengthFromLengths a b c d e = Intrinsic.flippedLengthFromLengths a b c d e := by
  simp only [Ddg.flippedLengthFromLengths, Intrinsic.flippedLengthFromLengths, flippedLength_eq,
    diamondX_eq, diamondH_eq]

theorem flipArea1_eq (xk hk xl hl : ℝ) :
    Ddg.flipArea1 xk hk xl hl = Intrinsic.flipArea1 xk hk xl hl := rfl

theorem flipArea2_eq (a xk hk xl hl : ℝ) :
    Ddg.flipArea2 a xk hk xl hl = Intrinsic.flipArea2 a xk hk xl hl := rfl

theorem flipConvex_eq (a xk hk xl hl : ℝ) :
    Ddg.flipConvex a xk hk xl hl = Intrinsic.flipConvex a xk hk xl hl := by
  rw [Bool.eq_iff_iff, Intrinsic.flipConvex_iff]
  simp only [Ddg.flipConvex, Bool.and_eq_true, decide_eq_true_eq, zero_eq, flipArea1_eq,
    flipArea2_eq]

theorem unfoldApex_eq (ax ay dx dy lab lbc lca : ℝ) :
    Ddg.unfoldApex ax ay dx dy lab lbc lca = Intrinsic.unfoldApex ax ay dx dy lab lbc lca := by
  simp only [Ddg.unfoldApex, Intrinsic.unfoldApex, two_eq, triArea_eq]

/-! ## The Delaunay condition and the cotan weight -/

theorem delaunayExcess_eq (θk θl : ℝ) :
    Ddg.delaunayExcess θk θl = Intrinsic.delaunayExcess θk θl := by
  simp only [Ddg.delaunayExcess, Intrinsic.delaunayExcess, pi_eq]

theorem delaunayHolds_eq (tol θk θl : ℝ) :
    Ddg.delaunayHolds tol θk θl = Intrinsic.delaunayHolds tol θk θl := by
  rw [Bool.eq_iff_iff, Intrinsic.delaunayHolds_iff]
  simp only [Ddg.delaunayHolds, decide_eq_true_eq, delaunayExcess_eq]

theorem cotanWeight_eq (x y : ℝ) : Ddg.cotanWeight x y = Intrinsic.cotanWeight x y := rfl

theorem cotanWeightFromLengths_eq (a b c d e : ℝ) :
    Ddg.cotanWeightFromLengths a b c d e = Intrinsic.cotanWeightFromLengths a b c d e := by
  simp only [Ddg.cotanWeightFromLengths, Intrinsic.cotanWeightFromLengths, cotanWeight_eq,
    cornerCot_eq]

theorem oppositeAngles_eq {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    Ddg.oppositeAngles a b c d e = Intrinsic.oppositeAngles a b c d e := by
  simp only [Ddg.oppositeAngles, Intrinsic.oppositeAngles, cornerAngle_eq hk, cornerAngle_eq hl]

/-! ## Angle normalization and the signpost update -/

theorem angleScale_eq (b : Bool) (Θ : ℝ) : Ddg.angleScale b Θ = Intrinsic.angleScale b Θ := by
  cases b <;> simp only [Ddg.angleScale, Intrinsic.angleScale, two_eq, pi_eq] <;> rfl

theorem reduceFullTurn_eq (x : ℝ) : Ddg.reduceFullTurn x = Intrinsic.reduceFullTurn x := by
  simp only [Ddg.reduceFullTurn, Intrinsic.reduceFullTurn, two_eq, pi_eq]

theorem signpostUpdate_eq (b : Bool) (φ θ Θ : ℝ) :
    Ddg.signpostUpdate b φ θ Θ = Intrinsic.signpostUpdate b φ θ Θ := by
  cases b <;> simp only [Ddg.signpostUpdate, Intrinsic.signpostUpdate, angleScale_eq,
    reduceFullTurn_eq] <;> rfl

/-! ## Tracing -/

theorem axisCrossingT_eq (yT yH : ℝ) : Ddg.axisCrossingT yT yH = Intrinsic.axisCrossingT yT yH :=
  rfl

theorem axisCrossingX_eq (xT xH t : ℝ) :
    Ddg.axisCrossingX xT xH t = Intrinsic.axisCrossingX xT xH t := rfl

theorem apexAboveAxis_eq (y : ℝ) : Ddg.apexAboveAxis y = Intrinsic.apexAboveAxis y := by
  rw [Bool.eq_iff_iff, Intrinsic.apexAboveAxis_iff]
  simp only [Ddg.apexAboveAxis, decide_eq_true_eq, zero_eq]

/-! ## Barycentric coordinates, distances, the circumcentre -/

theorem signedArea2_eq (ax ay bx bY cx cy : ℝ) :
    Ddg.signedArea2 ax ay bx bY cx cy = Intrinsic.signedArea2 ax ay bx bY cx cy := rfl

theorem barycentric_eq (ax ay bx bY cx cy px py : ℝ) :
    Ddg.barycentric ax ay bx bY cx cy px py = Intrinsic.barycentric ax ay bx bY cx cy px py := rfl

theorem pointFromDistances_eq (lab lbc lca ra rb rc : ℝ) :
    Ddg.pointFromDistances lab lbc lca ra rb rc =
      Intrinsic.pointFromDistances lab lbc lca ra rb rc := by
  simp only [Ddg.pointFromDistances, Intrinsic.pointFromDistances, two_eq, diamondX_eq,
    diamondH_eq]

theorem planeDistance_eq (px py qx qy : ℝ) :
    Ddg.planeDistance px py qx qy = Intrinsic.planeDistance px py qx qy := by
  simp only [Ddg.planeDistance, Intrinsic.planeDistance, sqrt_eq]

theorem cornerDistances_eq (l01 l12 l20 b0 b1 b2 : ℝ) :
    Ddg.cornerDistances l01 l12 l20 b0 b1 b2 = Intrinsic.cornerDistances l01 l12 l20 b0 b1 b2 := by
  simp only [Ddg.cornerDistances, Intrinsic.cornerDistances, planeDistance_eq, zero_eq,
    diamondX_eq, diamondH_eq]

theorem edgeSplitLengths_eq (lij bi bj : ℝ) :
    Ddg.edgeSplitLengths lij bi bj = Intrinsic.edgeSplitLengths lij bi bj := rfl

theorem circumWeight_eq (ak ak1 ak2 : ℝ) :
    Ddg.circumWeight ak ak1 ak2 = Intrinsic.circumWeight ak ak1 ak2 := rfl

theorem circumcenterBary_eq (a0 a1 a2 : ℝ) :
    Ddg.circumcenterBary a0 a1 a2 = Intrinsic.circumcenterBary a0 a1 a2 := rfl

/-! ## Moving a vertex, the optimal Delaunay triangulation -/

theorem movedLength_eq (l r c s : ℝ) : Ddg.movedLength l r c s = Intrinsic.movedLength l r c s := by
  simp only [Ddg.movedLength, Intrinsic.movedLength, sqrt_eq]

theorem movedDoubleArea_eq (A2 : ℝ) (v d : ℝ × ℝ) :
    Ddg.movedDoubleArea A2 v d = Intrinsic.movedDoubleArea A2 v d := rfl

theorem weightedCircumcenter_eq (lab lac lbc A : ℝ) :
    Ddg.weightedCircumcenter lab lac lbc A = Intrinsic.weightedCircumcenter lab lac lbc A := by
  have e8 : (@OfNat.ofNat ℝ 8 Ddg.instOfNatScalar) = (8 : ℝ) := by rw [ofNat_eq]
  simp only [Ddg.weightedCircumcenter, Intrinsic.weightedCircumcenter, two_eq, e8]

theorem weightedBarycenter_eq (lab lac lbc A : ℝ) :
    Ddg.weightedBarycenter lab lac lbc A = Intrinsic.weightedBarycenter lab lac lbc A := by
  simp only [Ddg.weightedBarycenter, Intrinsic.weightedBarycenter, two_eq, three_eq, diamondX_eq]

theorem odtDisplacement_eq (w : ℝ × ℝ) (total step : ℝ) :
    Ddg.odtDisplacement w total step = Intrinsic.odtDisplacement w total step := rfl

end Oracle

end DDG
