import DDG.Conformal.BoundaryCurve
import DDG.Conformal.CETM.Curvature
import DDG.Conformal.CETM.Flatten
import DDG.Conformal.CETM.Interpolation
import DDG.Conformal.CETM.Monotone
import DDG.Conformal.Cherrier
import DDG.Conformal.ConformalEnergy
import DDG.Conformal.Cut
import DDG.Conformal.Distortion
import DDG.Conformal.Sphere
import DDG.Conformal.Yamabe

/-!
# DDG.Conformal: discrete conformal geometry

**Scope.** Discrete conformal changes of metric and the flattening of surfaces: the linearized
Yamabe equation, Boundary First Flattening (BFF), distortion measures, and the CETM vertex
flattening used for intrinsic simplification.

**Main results.**
* `Yamabe`: under `ℓ̃_ij = e^{(u_i + u_j)/2} ℓ_ij` the corner angles change to first order by the
  cotan weights (`hasDerivAt_cornerAngle_conformal`, `hasDerivAt_angleDefect_conformalScale`).
* `Cherrier`: the exterior-angle sum (`gauss_bonnet_conformal`), Cherrier's formula
  (`hasDerivAt_geodesicCurvature_conformalScale`) and the Hilbert transform from Neumann data.
* `ConformalEnergy`: the conformal energy of a flattening over an oriented mesh as a shoelace sum.
* `BoundaryCurve`: dual and target boundary lengths, closing the exterior angles, the boundary
  tangents and curve. `Cut`: cut meshes, angle sums at cones and the Euler characteristic after
  cutting (`cut_euler_characteristic`). `Sphere`: stereographic projection and the inversions
  behind the Möbius centering.
* `Distortion`: singular values and the quasi-conformal error of a triangle map (`qcError`).
* `CETM/*`: flattening one vertex, monotonicity of its angle sum (`angleSum_strictAnti_of_valid`),
  conservation of curvature (`flatten_curvature_conserved`), the redistribution weights
  (Eq. 9) and projective interpolation (Eq. 10).

**Sources.** M. Sawhney, K. Crane, *Boundary First Flattening*, ACM TOG 2017 (Sec. 4.3, 4.5, 6.5,
7.5, Eq. 15–20); B. Springborn, P. Schröder, U. Pinkall, *Conformal equivalence of triangle
meshes*, ACM TOG 2008 (Sec. 3, 4.1, App. B, D); M. Baden, K. Crane, M. Kazhdan, *Möbius
Registration*, CGF 37(5), SGP 2018 (Alg. 1); H.-T. D. Liu, M. Gillespie, B. Chislett, N. Sharp,
A. Jacobson, K. Crane, *Surface Simplification using Intrinsic Error Metrics*, ACM TOG 42(4), 2023
(Sec. 4.1, 5.4, 7.1).

**Depends on:** `DDG.FEM`, `DDG.Laplace`, `DDG.Mesh`.
-/
