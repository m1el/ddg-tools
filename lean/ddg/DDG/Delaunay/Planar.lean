import DDG.Delaunay.Planar.Basic
import DDG.Delaunay.Planar.Encroachment
import DDG.Delaunay.Planar.Existence
import DDG.Delaunay.Planar.Flip
import DDG.Delaunay.Planar.FlipCount
import DDG.Delaunay.Planar.Insertion
import DDG.Delaunay.Planar.Lemma
import DDG.Delaunay.Planar.Termination
import DDG.Delaunay.Planar.Uniqueness

/-!
# DDG.Delaunay.Planar: Delaunay triangulations of planar point sets

**Scope.** The classical planar theory, in coordinates (points of `ℝ × ℝ`, orientation by the signed
double area `sda`, in-circle determinant `inCircle`); no mesh data structure. This part is
independent of the surface theory and can be used alone.

**Main results.**
* `Basic`: triangles, circumcircles and the lifting map in coordinates; a non-Delaunay diamond is
  strictly convex.
* `Flip`: triangulations of a finite point set (`IsTriangulation`), the Lawson flip
  (`DelaunayFlip`, `isTriangulation_of_delaunayFlip`).
* `Termination`: Lawson's algorithm terminates, by the lifting potential `liftPotential`.
* `Lemma`: **Delaunay's lemma** (`delaunay_lemma`: locally Delaunay implies empty circumcircles);
  from every triangulation flips reach a Delaunay one (`exists_flips_to_delaunay`).
* `Existence`: every non-collinear finite point set has a triangulation (`exists_triangulation`),
  hence a Delaunay triangulation (`exists_delaunay`); `Insertion`: inserting a point outside the
  hull (`isTriangulation_insert`).
* `Uniqueness`: unique when no four points are cocircular (`delaunay_unique`,
  `existsUnique_delaunay`); any two locally Delaunay triangulations coincide
  (`eq_of_isLocallyDelaunay`).
* `FlipCount`: an edge flipped away never reappears (`flipped_edge_not_mem`), so a flip sequence
  has at most `C(n, 2)` flips (`flip_count_le`).
* `Encroachment`: Ruppert's encroachment and the circumradius.

**Sources.** B. Delaunay, 1934; C. L. Lawson, 1977; H. Edelsbrunner, *Geometry and Topology for Mesh
Generation*, Cambridge University Press, 2001; M. de Berg, O. Cheong, M. van Kreveld, M. Overmars,
*Computational Geometry: Algorithms and Applications*, Ch. 9; J. Ruppert, *A Delaunay refinement
algorithm for quality 2-dimensional mesh generation*, 1995.

**Depends on:** `DDG.Geometry` and, for the triangle as a set and its measure only, `DDG.FEM`.
-/
