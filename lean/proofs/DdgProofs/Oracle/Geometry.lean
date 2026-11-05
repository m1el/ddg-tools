import DdgProofs.Oracle.Scalar
import DDG.Geometry.Barycentric

/-!
# Oracle bridges: the level crossing and the barycentric formulas at `ℝ`

The oracle's (`lean/oracle`, `DdgOracle/Formulas/Cuts.lean`) scalar-generic definitions of the
crossing parameter, the barycentric lengths, the barycentric inner product, the cotangent of a
sub-triangle and its area, instantiated at `ℝ`, are the native definitions of
`DDG/Geometry/Barycentric.lean` (library). One-way bridges `Ddg.foo x = DDG.foo x`: the
theorems of the library apply to the oracle's numbers after rewriting with them.

App layer (library step S3).
-/

noncomputable section

namespace DDG

namespace Oracle

theorem crossingParam_eq (φa φb : ℝ) : Ddg.crossingParam φa φb = levelCrossingParam φa φb := rfl

theorem baryLength2_eq (lab lbc lca da db dc : ℝ) :
    Ddg.baryLength2 lab lbc lca da db dc = baryLength2 lab lbc lca da db dc := rfl

theorem baryLength_eq (lab lbc lca da db dc : ℝ) :
    Ddg.baryLength lab lbc lca da db dc = baryLength lab lbc lca da db dc := by
  simp only [Ddg.baryLength, baryLength, sqrt_eq, zero_eq, baryLength2_eq]

theorem det3_eq (x0 x1 x2 y0 y1 y2 z0 z1 z2 : ℝ) :
    Ddg.det3 x0 x1 x2 y0 y1 y2 z0 z1 z2 = det3 x0 x1 x2 y0 y1 y2 z0 z1 z2 := rfl

theorem subTriangleArea_eq (det A : ℝ) : Ddg.subTriangleArea det A = subTriangleArea det A := by
  simp only [Ddg.subTriangleArea, subTriangleArea, abs_eq]

theorem baryInner_eq (lab lbc lca da db dc ea eb ec : ℝ) :
    Ddg.baryInner lab lbc lca da db dc ea eb ec = baryInner lab lbc lca da db dc ea eb ec := by
  simp only [Ddg.baryInner, baryInner, two_eq]

theorem cotFromInner_eq (x A : ℝ) : Ddg.cotFromInner x A = cotFromInner x A := by
  simp only [Ddg.cotFromInner, cotFromInner, two_eq]

end Oracle

end DDG
