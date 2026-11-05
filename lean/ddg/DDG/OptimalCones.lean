import DDG.OptimalCones.BoundaryProblem
import DDG.OptimalCones.Certificate
import DDG.OptimalCones.Closed
import DDG.OptimalCones.Components
import DDG.OptimalCones.Duality
import DDG.OptimalCones.MeshProblem
import DDG.OptimalCones.Newton
import DDG.OptimalCones.Optimality
import DDG.OptimalCones.Problem

/-!
# DDG.OptimalCones: the convex cone-placement problem

**Scope.** The canonical discrete problem of optimal cone singularities for conformal flattening
(a convex problem with a cotan-Laplacian state equation), its optimality theory and the facts that
let a solver certify a solution. Exact real arithmetic; solver control flow is in the app layer.

**Main results.**
* `Problem`: `Site`, `ConeProblem`, objective, feasibility, the optimality system `KKT`.
* `Optimality`: the optimality system is sufficient (`KKT.objective_add_le`) and, for invertible
  `L`, necessary (`IsMinimizer.existsUnique_kkt`); uniqueness (`IsMinimizer.unique`); the
  least-distortion guarantee of Sec. 4.5.
* `Duality`: the discrete dual, Fenchel–Young, weak and strong duality, the duality-gap
  certificate (`duality_gap_certificate`). `Certificate`: the certificate for inexact data
  (`certificate_enclosure`, `dual_lower_bound_of_clip`), with counterexamples.
* `Newton`: the residual maps, semismooth Newton on fixed active sets and its stopping rule
  (`paper_stageI_stopping_rule`, `paper_stageII_stopping_rule`), the regularized problem.
* `Closed`: closed meshes, the optimal additive constant, area normalization of `λ`.
  `Components`: the problem separates over connected components.
* `MeshProblem`, `BoundaryProblem`: the canonical problem of a mesh, with boundary unknowns.

**Sources.** M. Soliman, D. Slepčev, K. Crane, *Optimal Cone Singularities for Conformal
Flattening*, ACM TOG 37(4), 2018 (Problem 10 / A.1, Eq. 17–24, Sec. 4.5, 5.3, App. B, C).

**Depends on:** `DDG.Optimization`, `DDG.Conformal`.
-/
