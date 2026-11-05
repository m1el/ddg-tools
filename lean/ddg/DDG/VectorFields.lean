import DDG.VectorFields.Holonomy
import DDG.VectorFields.PoincareHopf

/-!
# DDG.VectorFields: holonomy and the discrete Poincaré–Hopf theorem

**Scope.** Tangent vector fields on a triangulated surface in the signpost frames of an intrinsic
triangulation: the closed forms of parallel transport and the connection Laplacian, the index of a
face, and the index sum.

**Main results.**
* `Holonomy`: transport around a face (`transport_around_face`), the entries of the connection
  Laplacian (`edge_entries_hermitian`), the singular values of a conformal map from its parts
  (`singular_values_from_parts`).
* `PoincareHopf`: the integer index of a face (`faceIndex_isInt`) and the Poincaré–Hopf sum
  on a closed oriented mesh (`sum_faceIndex`, `sum_round_faceIndexFixed`): `index_sum = nχ`.

**Sources.** N. Sharp, Y. Soliman, K. Crane, *Navigating Intrinsic Triangulations*, ACM TOG 2019,
Sec. 5.4–5.5 (the vector heat method, VHM Eq. 10); F. Knöppel, K. Crane, U. Pinkall, P. Schröder,
*Globally Optimal Direction Fields*, ACM TOG 2013 (Eq. 19, App. B).

**Depends on:** `DDG.Mesh`.
-/
