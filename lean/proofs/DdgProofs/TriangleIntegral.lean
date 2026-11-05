import DDG.FEM.TriangleIntegral
import DDG.FEM.Elements
import DDG.FEM.Tiling
import DDG.FEM.Flux
import DDG.FEM.BetaExp

/-!
# Integrals over a triangle

`Basic`: the standard triangle, Fubini, the Beta integral, the affine change of variables, the
area, the convex hull, barycentric coordinates, the slice formula and the area-coordinate rule.
`Fem`: the finite element integrals built on them (Galerkin mass, cotan stiffness, weighted
stiffness, the exponential slice of the connection Laplacian, the conformal energy of a linear
map, Chen–Holst Lemma 4.1, the element norms of the AMR estimator) and their sums over a mesh.
`Tiling`: additivity of the integral over triangles with disjoint interiors, Lemma 4.1 over a
tiling and the descent of the interpolation error under the ODT update, the transfer mass matrix
as an integral. `Flux`: the flux of a linear function through a side (AMR jump residuals).
`BetaExp`: the series and the recurrence by which `beta_exp` (`knoeppel.rs`) evaluates
`∫₀¹ s^p (1 − s)^m e^{iΩs} ds`.
-/
