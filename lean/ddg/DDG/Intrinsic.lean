import DDG.Intrinsic.Barycentric
import DDG.Intrinsic.Diamond
import DDG.Intrinsic.FlatVertexRemoval
import DDG.Intrinsic.FlipAngles
import DDG.Intrinsic.Mollification
import DDG.Intrinsic.Move
import DDG.Intrinsic.Signpost
import DDG.Intrinsic.Strip
import DDG.Intrinsic.Trace
import DDG.Intrinsic.VertexRemoval

/-!
# DDG.Intrinsic: the local theory of intrinsic triangulations

**Scope.** Closed forms of the single operations on an intrinsic triangulation (a triangulation
given by edge lengths and signposts), each a statement about a diamond, a star or one tracing step;
the data structure and global algorithms built on them are not modelled here.

**Main results.**
* `Diamond`: unfolding a triangle and a diamond, uniqueness up to a rigid motion
  (`diamond_unique`), the flipped edge length (`flippedLengthFromLengths`), strict convexity.
* `FlipAngles`: what a flip does to the corner angles (`flip_angle_sums`).
* `Signpost`: angle normalization and the signpost update around interior and boundary vertices
  and under a flip (`signpostUpdate`, `flip_signpost_at_k`).
* `Trace`: one step of a tracing query (`axisCrossingT`, `exit_cases`).
* `Barycentric`: barycentric coordinates from lengths, distances to the corners, edge split,
  circumcentre (`circumcenterBary`).
* `Move`: moving a vertex, the weighted circumcentre (`weightedCircumcenter_equidistant`).
* `Strip`, `VertexRemoval`, `FlatVertexRemoval`: the intersection of an unfolded strip, removal of
  an inserted or a flat vertex by flips (`exists_flippable`, `flatStar_exists_flippable`).
* `Mollification`: intrinsic mollification of edge lengths (`isTriangle_mollified`).

**Sources.** N. Sharp, Y. Soliman, K. Crane, *Navigating Intrinsic Triangulations*, ACM TOG 2019
(Sec. 2.3, 3.2, 3.3, 4, Alg. 1–10); M. Gillespie, N. Sharp, K. Crane, *Integer Coordinates for
Intrinsic Geometry Processing*, ACM TOG 2021 (Sec. 3.5, App. D, Thm. D.1); N. Sharp, K. Crane,
*You Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges*, ACM TOG 2020 (FlipOut);
H.-T. D. Liu et al., *Surface Simplification using Intrinsic Error Metrics*, ACM TOG 2023
(Sec. 4.2); L. Chen, M. Holst, CMAME 2011 (vertex move); M. Sawhney, K. Crane, *Boundary First
Flattening*, ACM TOG 2017 (mollification).

**Depends on:** `DDG.Geometry`.
-/
