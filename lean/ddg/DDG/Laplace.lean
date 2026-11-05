import DDG.Laplace.Cotan
import DDG.Laplace.Dirichlet
import DDG.Laplace.Kernel
import DDG.Laplace.Poisson
import DDG.Laplace.Weighted

/-!
# DDG.Laplace: the cotan Laplacian and linear boundary problems

**Scope.** The cotan formula and the cotan Laplacian of a triangle mesh, lumped mass, the weighted
variant, and the linear algebra of Dirichlet, Poisson, Neumann and Dirichlet-to-Neumann problems.

**Main results.**
* `Cotan`: the cotan formula for an affine function (`cotan_formula`), the Laplacian
  `TriMesh.laplacian` (symmetric `laplacian_isSymm`, positive semidefinite
  `TriMesh.laplacian_posSemidef`, kernel the constants `TriMesh.laplacian_mulVec_eq_zero_iff`,
  interior block positive definite `TriMesh.laplacian_interior_posDef`), lumped mass (`sum_mass`).
* `Dirichlet`: the Dirichlet problem and the harmonic extension for a square matrix
  (`existsUnique_isDirichletSolution`, `existsUnique_harmonicExtension`, energy minimality).
* `Poisson`: Poisson and Neumann problems of a mesh (`existsUnique_poisson`,
  `neumannData_spec`, `existsUnique_neumann`); Dirichlet-to-Neumann.
* `Kernel`: a linear system with zero column sums has one redundant equation
  (`mulVec_eq_of_sum_eq`).
* `Weighted`: the weighted cotan Laplacian `TriMesh.weightedLaplacian`.

**Sources.** M. Sawhney, K. Crane, *Boundary First Flattening*, ACM TOG 2017 (Poisson,
Dirichlet-to-Neumann and Neumann problems); N. Sharp, K. Crane, *Variational Surface Cutting*,
ACM TOG 2018 (the weighted Laplacian).

**Depends on:** `DDG.Mesh`.
-/
