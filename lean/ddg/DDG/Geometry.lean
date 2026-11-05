import DDG.Geometry.Advancement
import DDG.Geometry.AdvancementPairs
import DDG.Geometry.Barycentric
import DDG.Geometry.ClosestPoints
import DDG.Geometry.HopfLax
import DDG.Geometry.Plane
import DDG.Geometry.PlaneLemmas
import DDG.Geometry.RigidInvariance
import DDG.Geometry.SegmentCrossing
import DDG.Geometry.SignedAngle
import DDG.Geometry.SimilarityFit
import DDG.Geometry.Stable
import DDG.Geometry.Triangle
import DDG.Geometry.Vec3

/-!
# DDG.Geometry: triangles, planes, closest points

**Scope.** Euclidean geometry of triangles given by their edge lengths, and the plane and space
primitives the other parts build on; no meshes. The library is intrinsic: a triangle quantity is a
function of its three side lengths (Sharp, Soliman, Crane 2019, Sec. 5.1.1).

**Main results.**
* `Triangle`: `cornerAngle`, `cornerCot`, `triArea` and their relations (law of cosines
  `cos_cornerAngle`, Heron `sixteen_mul_triArea_sq`, angle sum
  `cornerAngle_add_cornerAngle_add_cornerAngle`, scaling), and the link to Mathlib's Euclidean geometry.
* `Stable`: the same quantities in the arrangements a numerical implementation evaluates
  (Kahan's parenthesization with sorted sides, `atan2` forms: `kahanCornerCot_eq_cornerCot`).
* `Plane`, `PlaneLemmas`, `SignedAngle`, `SegmentCrossing`, `SimilarityFit`: plane primitives
  (`dist2`, signed double area `sda`, the reference layout of a triangle), signed corner angles,
  proper crossing of two segments, the least-squares similarity fit of two point sets.
* `Barycentric`: barycentric combinations, level crossings (`crossingParam`), lengths and areas
  inside a face (`baryLength`).
* `HopfLax`: the Hopf–Lax (local eikonal) update on a triangle in closed form is the minimum over
  the edge (`hopfLaxTriangle_isLeast`; `cut-redistance`).
* `Vec3`, `RigidInvariance`: vectors of `ℝ³` and derivatives along motions; invariance under rigid
  motions and the identities `Σ g = 0`, `Σ x × g = 0` of a gradient (`gradient_invariance`).
* `ClosestPoints`, `Advancement`, `AdvancementPairs`: closest points of segments and triangles
  (`segSegDist_attained`, `closestPtTri_isClosest`) and conservative advancement
  (`primDist_advance`, `advance_safe`).

**Sources.** K. Crane, *Discrete Differential Geometry: An Applied Introduction*; N. Sharp,
Y. Soliman, K. Crane, *Navigating Intrinsic Triangulations*, ACM TOG 2019; C. Ericson,
*Real-Time Collision Detection*, 2004, Sec. 5.1.2, 5.1.5, 5.1.9; M. Li, D. M. Kaufman, C. Jiang,
*Codimensional Incremental Potential Contact*, ACM TOG 2021 (Additive CCD); M. Sawhney, K. Crane,
*Boundary First Flattening*, ACM TOG 2017 (plane measures); N. Sharp, K. Crane, *Variational Surface
Cutting*, ACM TOG 2018 (barycentric lengths); C. Yu, H. Schumacher, K. Crane, *Repulsive Curves*,
ACM TOG 2021 and C. Yu, K. Brakensiek, H. Schumacher, K. Crane, *Repulsive Surfaces*, ACM TOG 2021
(vectors, advancement); O. Stein, E. Grinspun, K. Crane, *Developability of Triangle Meshes*,
ACM TOG 2018 (rigid invariance).

**Depends on:** Mathlib only.
-/
