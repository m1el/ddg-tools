import DDG.FEM.BetaExp
import DDG.FEM.Elements
import DDG.FEM.Flux
import DDG.FEM.Heat
import DDG.FEM.Tiling
import DDG.FEM.Transfer
import DDG.FEM.TriangleIntegral

/-!
# DDG.FEM: integrals over triangles and finite element formulas

**Scope.** Lebesgue integrals over a triangle of the plane and the finite element quantities they
give: the P1 Galerkin mass and stiffness matrices, tilings, element formulas of the heat method,
the connection-Laplacian exponential integrals, side fluxes, and the transfer of attributes.

**Main results.**
* `TriangleIntegral`: the standard triangle, the Beta integral (`integral_stdTri`), the affine
  change of variables (`integral_triangle`), area, barycentric coordinates and the area-coordinate
  rule (`integral_bary_monomial`).
* `Elements`: Galerkin mass, cotan and weighted stiffness and their sums over a mesh, the
  conformal energy of a linear map, Chen–Holst Lemma 4.1.
* `Tiling`: additivity over a tiling by triangles with disjoint interiors (`IsTiling`,
  `integral_biUnion_triangle`), Lemma 4.1 over a tiling, the transfer mass as an integral.
* `Flux`: the flux of a linear function through a side. `Heat`: gradient and divergence of a
  triangle (`divergence_corner`).
* `BetaExp`: the series and recurrence for `∫₀¹ s^p (1 − s)^m e^{iΩs} ds` (`betaExp`,
  `integral_bary_slice_exp`).
* `Transfer`: the normal equations of an L² attribute transfer (`transfer_minimizer_iff`).

**Sources.** L. Chen, M. Holst, *Efficient mesh optimization schemes based on Optimal Delaunay
Triangulations*, CMAME 2011 (Lemma 4.1); N. Sharp, Y. Soliman, K. Crane, *Navigating Intrinsic
Triangulations*, ACM TOG 2019, Sec. 5.2–5.3; K. Crane, C. Weischedel, M. Wardetzky, *Geodesics in
Heat*, ACM TOG 2013; M. Gillespie, N. Sharp, K. Crane, *Integer Coordinates for Intrinsic Geometry
Processing*, ACM TOG 2021, Sec. 4.3 (transfer); F. Knöppel et al., *Globally Optimal Direction
Fields*, ACM TOG 2013 (`BetaExp`).

**Depends on:** `DDG.Laplace`.
-/
