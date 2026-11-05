import DDG.SDF.Basic
import DDG.SDF.Blend
import DDG.SDF.Circle2d
import DDG.SDF.Features
import DDG.SDF.MongeFit
import DDG.SDF.Operations
import DDG.SDF.Torus

/-!
# DDG.SDF: signed distance functions from points (Points as Tori)

**Scope.** Machine-checked derivations of the formulas of Points as Tori: signed distance to a
torus, fitting a torus to a point and its normal from a local quadratic patch, blending the torus
SDFs of neighbouring points, and the operations on the resulting field. The smooth theory
(Laplace's method, the signed Hopf–Cole transformation) and the network are cited
(`lean/proofs/CITED.md`). All theorems are in namespace `DDG.PAT`; vectors are `DDG.V3`.

**Main results.**
* `Torus`: the exact torus SDF `torusSDF`, its gradient and stable form (`torus_sdf_gradient`,
  `torus_sdf_stable`), and the contact of the fitted torus with the patch
  (`torus_contact_value`, `torus_contact_normal`).
* `MongeFit`: the fundamental forms of a Monge patch, mean and Gauss curvature, principal
  curvatures and directions, torus radii and sign (`fundamental_forms`, `mean_gauss`,
  `principal_pair`, `torus_radii`). `Features`: neighbourhood scale, local frame and features.
* `Blend`: the softmin blend, shift invariance, convexity and exact gradient (`blend_shift`,
  `blend_gradient_exact`, `blend_laplacian`).
* `Operations`: orientation metric, marching-cubes interpolation (`mc_interpolation`), Booleans
  with Lipschitz bounds, revolve and extrude, sphere tracing (`sphere_trace_step_safe`).
* `Circle2d`: the 2D variant with circles (`circle_fit`, `circle_sdf`).

**Sources.** L. Feng, I. Gkioulekas, K. Crane, *Points as Tori: Fast Pointwise Signed Distance for
Point Clouds*, ACM TOG 45(4), Article 53, 2026 (Eq. 23–26, Sec. 4.1, App. C, D); the contract is
`crates/ddg-points-as-tori/SPEC.md`. The Rust crate was developed with the authors' MIT-licensed
implementation `points-as-tori` as reference (see its README).

**Depends on:** `DDG.Geometry`.
-/
