import DdgProofs.Cuts.Crossing
import DDG.Geometry.HopfLax

/-!
# Cut flow: the oracle's Hopf–Lax update is the library's

`crates/ddg-variational-cuts/SPEC.md` §3.12, D10 (Bornemann, Rasch, 2006). The mathematics of the
local update on a triangle, `min_{λ ∈ [0,1]} λ T_j + (1 − λ) T_k + ‖x_i − (λ x_j + (1 − λ) x_k)‖`
(`cut-redistance`: the closed form is the minimum, it is bounded by the edge terms, and it is at
least `min(T_j, T_k)`), is in the library, `DDG/Geometry/HopfLax.lean` (`DDG.hopfLaxTriangle`,
`hopfLaxTriangle_isLeast`, `hopfLaxTriangle_isLeast_intrinsic`, ...). (Library step S6: moved
there.)

What stays here is the one-way bridge from the oracle's scalar-generic definitions
(`lean/oracle/DdgOracle/Formulas/Cuts.lean`, instantiated at `ℝ`) to the library's, and the
oracle forms of the main statements, so that the oracle's `Ddg.hopfLaxTriangle` (the function the
cut flow's redistancing calls) is shown to be the minimum.

What is **not** covered: the marching order and the single pass (a definition, SPEC §3.12), and
that the scheme converges to the geodesic distance (`cut-fast-marching-convergence`, cited).
-/

noncomputable section

namespace DDG

namespace Cuts

open scoped DDG.Oracle

/-- The oracle's edge term at `ℝ` is the library's. -/
theorem hopfLaxEdge_eq (Tj lij : ℝ) : Ddg.hopfLaxEdge Tj lij = hopfLaxEdge Tj lij := rfl

/-- **The oracle's `hopfLaxTriangle` at `ℝ` is the library's.** -/
theorem hopfLaxTriangle_eq (lij lik ljk Tj Tk : ℝ) :
    Ddg.hopfLaxTriangle lij lik ljk Tj Tk = hopfLaxTriangle lij lik ljk Tj Tk := by
  unfold Ddg.hopfLaxTriangle hopfLaxTriangle hlValue hlQ hlLam
  simp only [Oracle.sqrt_eq, Oracle.abs_eq, Oracle.zero_eq, Oracle.one_eq, Oracle.two_eq]

/-- **The oracle's closed form is the minimum** over the edge, for points of an inner product
space. -/
theorem Ddg_hopfLaxTriangle_isLeast {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (xi xj xk : E) (Tj Tk : ℝ) :
    IsLeast (hopfLaxObjective xi xj xk Tj Tk '' Set.Icc 0 1)
      (Ddg.hopfLaxTriangle (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk) := by
  rw [hopfLaxTriangle_eq]
  exact hopfLaxTriangle_isLeast xi xj xk Tj Tk

/-- The same from the side lengths alone, with the oracle's `baryLength`. -/
theorem Ddg_hopfLaxTriangle_isLeast_intrinsic {lij lik ljk : ℝ} (h : IsTriangle lij lik ljk)
    (Tj Tk : ℝ) :
    IsLeast ((fun lam : ℝ => lam * Tj + (1 - lam) * Tk +
        Ddg.baryLength lij ljk lik 1 (-lam) (-(1 - lam))) '' Set.Icc 0 1)
      (Ddg.hopfLaxTriangle lij lik ljk Tj Tk) := by
  rw [hopfLaxTriangle_eq]
  simpa only [Oracle.baryLength_eq] using hopfLaxTriangle_isLeast_intrinsic h Tj Tk

end Cuts

end DDG
