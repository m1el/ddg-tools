import DDG.Delaunay.BoundaryFan
import DDG.Delaunay.FlipQueue
import DDG.Delaunay.FromTriMesh
import DDG.Delaunay.Local
import DDG.Delaunay.Metric
import DDG.Delaunay.ODT
import DDG.Delaunay.Planar
import DDG.Delaunay.PlanarBridge
import DDG.Delaunay.Potential
import DDG.Delaunay.Rippa
import DDG.Delaunay.Surface.Angle
import DDG.Delaunay.Surface.Count
import DDG.Delaunay.Surface.Fan
import DDG.Delaunay.Surface.Finite
import DDG.Delaunay.Surface.FlipEdge
import DDG.Delaunay.Surface.FlipFace
import DDG.Delaunay.Surface.Layout
import DDG.Delaunay.Surface.Overlay
import DDG.Delaunay.Surface.Plane
import DDG.Delaunay.Surface.Quad
import DDG.Delaunay.Surface.Relayout
import DDG.Delaunay.Surface.Segment
import DDG.Delaunay.Termination
import DDG.Delaunay.Planar

/-!
# DDG.Delaunay: intrinsic Delaunay triangulations

**Scope.** Delaunay triangulations of piecewise flat surfaces, from one diamond up to the flip
algorithm on a whole surface, and optimal Delaunay triangulations. The planar theory is the
separate part `DDG.Delaunay.Planar`, to which `PlanarBridge` connects.

**Main results (layers 1 to 4).**
* `Local` (layer 1, one diamond from five lengths): `θ_k + θ_l ≤ π` iff `cot θ_k + cot θ_l ≥ 0` iff
  the other apex is not inside the circumcircle (`delaunayHolds_iff`,
  `delaunay_iff_cotanWeight_nonneg`), and the flip of a non-Delaunay edge (`flipLength`).
  `Rippa`: Rippa's theorem (`rippa_identity`, `rippa_le`). `Potential`: the lifting energy and
  Musin's harmonic index decrease strictly under a Delaunay flip.
* `Metric`, `FromTriMesh`, `BoundaryFan` (layer 2, on gluing data of `DDG.Mesh`): edge-length
  metrics `IsMetric`, cone angles, Gauss–Bonnet (`IsMetric.gauss_bonnet`), the flip keeps a
  metric and every cone angle, a `TriMesh` as gluing data, at most one boundary halfedge per vertex.
* `Termination` (layer 3): the flip algorithm (`flip_algorithm`, `exists_delaunay_of_acc`);
  `FlipQueue`: the tolerant flip loop with a queue (`loop_result`).
* `Surface/*` (layer 4): the overlay of a triangulation on the original one (`Overlay`,
  `Overlay.flip`), a straight segment of bounded length crosses boundedly many edges, hence
  finitely many triangulations with bounded edges (`MetricGluing.boundedLengthFinite`) and
  termination of the flip algorithm unconditionally (`loop_acc'`). The plan and route are in
  `DDG/Delaunay/LAYER4-PLAN.md`.
* `PlanarBridge`: the one-diamond criterion is the planar in-circle test.
* `ODT`: optimal Delaunay triangulations, the vertex update and the interpolation error
  (`odtDisplacement`, `odtEnergy_step_le`, `IsTiling.integral_abs_odt_error`).

**Sources.** C. Indermitte, Th. M. Liebling, M. Troyanov, H. Clémençon, *Voronoi diagrams on
piecewise flat surfaces and an application to biological growth*, TCS 263, 2001; A. I. Bobenko,
B. Springborn, *A discrete Laplace–Beltrami operator for simplicial surfaces*, DCG 38(4), 2007;
F. Fisher, B. Springborn, A. I. Bobenko, P. Schröder, *An algorithm for the construction of
intrinsic Delaunay triangulations*, Computing 81, 2007; N. Sharp, Y. Soliman, K. Crane,
*Navigating Intrinsic Triangulations*, ACM TOG 38(4), 2019; S. Rippa, 1990; I. Musin, 1997;
L. Chen, J. Xu, *Optimal Delaunay triangulations*, 2004; L. Chen, M. Holst, *Efficient mesh
optimization schemes based on Optimal Delaunay Triangulations*, CMAME 200, 2011.

**Depends on:** `DDG.Delaunay.Planar`, `DDG.Intrinsic`, `DDG.FEM`, `DDG.Laplace`, `DDG.Mesh`.
-/
