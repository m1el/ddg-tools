import Mathlib

/-!
# Weak Wolfe steps satisfy the BFGS curvature condition

`crates/ddg-developability/SPEC.md` §4.2 (L-BFGS: a pair is stored only if `sᵀy > 0`) and §4.3
(Lewis–Overton Algorithm 2.6; A. S. Lewis, M. L. Overton, *Nonsmooth optimization via
quasi-Newton methods*, Math. Program. A 141, 2013).

* `Dev.weakWolfe_curvature` — a step satisfying the weak Wolfe condition `C(t)` along a descent
  direction gives `sᵀy > t (c2 − 1) gᵀd > 0` (so the BFGS update keeps `H` positive definite).

The statement is on real vectors and uses no oracle definition. (It was in
`DdgProofs.Developability.Optimizer`, where its hypothesis was the oracle's
`curvatureHolds c2 (gᵀd) (g₊ᵀd) = true`; that form is now `Dev.weakWolfe_curvature_of_curvatureHolds`
there.)
-/

noncomputable section

namespace DDG

namespace Dev

/-- **Weak Wolfe implies the curvature condition** (SPEC §4.2): for a step `t > 0` along a
descent direction `d` (`gᵀd < 0`) with `c2 < 1`, if the weak Wolfe condition
`C(t): c2 gᵀd < g₊ᵀd` holds, then the
BFGS pair `s = t d`, `y = g₊ − g` has `sᵀy > t (c2 − 1) gᵀd > 0`.

[formula: dev-weak-wolfe-curvature] -/
theorem weakWolfe_curvature {n : ℕ} (g gp d : Fin n → ℝ) {t c2 : ℝ} (ht : 0 < t)
    (hdesc : g ⬝ᵥ d < 0) (hc2 : c2 < 1)
    (hC : c2 * (g ⬝ᵥ d) < gp ⬝ᵥ d) :
    t * (c2 - 1) * (g ⬝ᵥ d) < (t • d) ⬝ᵥ (gp - g) ∧ 0 < (t • d) ⬝ᵥ (gp - g) := by
  have e : (t • d) ⬝ᵥ (gp - g) = t * (gp ⬝ᵥ d - g ⬝ᵥ d) := by
    rw [smul_dotProduct, dotProduct_sub, dotProduct_comm d gp, dotProduct_comm d g, smul_eq_mul]
  rw [e]
  have h1 : (c2 - 1) * (g ⬝ᵥ d) < gp ⬝ᵥ d - g ⬝ᵥ d := by linarith
  have h2 : 0 < (c2 - 1) * (g ⬝ᵥ d) := mul_pos_of_neg_of_neg (by linarith) hdesc
  constructor
  · have := mul_lt_mul_of_pos_left h1 ht
    linarith
  · exact mul_pos ht (by linarith)

end Dev

end DDG
