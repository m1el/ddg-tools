import DDG.Optimization.LBFGS
import DDG.Optimization.PiecewiseAffine
import DDG.Optimization.SoftThreshold
import DDG.Optimization.Wasserstein
import DDG.Optimization.Wolfe

/-!
# DDG.Optimization: generic convex and numerical optimization facts

**Scope.** Statements about optimization with no geometry in them; this part can be used on its
own. Proofs are over `ℝ` in exact arithmetic.

**Main results.**
* `SoftThreshold`: the soft threshold and the clamp, the fixed-point characterization
  `eq_softThreshold_add_iff`, the subdifferential of `t·|x|` on an interval
  (`isSubgradientWithinAt_abs_iff`).
* `PiecewiseAffine`: a semismooth Newton step that keeps its pieces solves a piecewise-affine
  equation exactly (`newton_stageII_exact`, `newton_stageI_exact`).
* `Wolfe`: a weak-Wolfe step satisfies the BFGS curvature condition `sᵀy > 0`
  (`weakWolfe_curvature`).
* `LBFGS`: the L-BFGS two-loop recursion is the BFGS inverse update (`twoLoop_eq`), which is
  positive definite (`bfgsInv_posDef`).
* `Wasserstein`: the memoryless cost of a vertex removal is a 1-Wasserstein distance
  (`memorylessCost_eq`, `wasserstein_eq`).

**Sources.** M. Soliman, D. Slepčev, K. Crane, *Optimal Cone Singularities for Conformal
Flattening*, ACM TOG 2018 (Eqs. 18–19, 21, Props. B.2–B.3); A. S. Lewis, M. L. Overton,
*Nonsmooth optimization via quasi-Newton methods*, Math. Program. 2013; C. Yu, K. Brakensiek, H. Schumacher, K. Crane, *Repulsive
Surfaces*, ACM TOG 2021, Sec. 8.2; H.-T. D. Liu, M. Gillespie, B. Chislett, N. Sharp, A. Jacobson,
K. Crane, *Surface Simplification using Intrinsic Error Metrics*, ACM TOG 2023, Sec. 5.1;
G. Peyré, M. Cuturi, *Computational Optimal Transport*, 2019.

**Depends on:** Mathlib only.
-/
