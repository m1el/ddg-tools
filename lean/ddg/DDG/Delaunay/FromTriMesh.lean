import DDG.Delaunay.Metric
import DDG.Mesh.Orientation

/-!
# A simplicial `TriMesh` as gluing data

The simplicial meshes of `DDG/Mesh/GaussBonnet.lean` (`TriMesh`: faces as vertex triples) with
a consistent orientation (`TriMesh.Oriented`, `DDG/Mesh/Orientation.lean`) are gluing data in the
sense of `DDG/Mesh/Gluing.lean`: the halfedges are the same pairs `(f, k)`, and the twin
is `TriMesh.twin` (the reversed halfedge, or the halfedge itself on the boundary). This is Rust
`Complex::from_tri_mesh` (`crates/ddg-intrinsic/src/complex.rs`): same faces, halfedges, edges and
vertices, with the same ids.

* `toGluing_isLabeling`: the tails of the halfedges (Rust `faces`) are a vertex labelling.
* `numEdges_toGluing`: the edges of the gluing data are the edges of the mesh (unordered vertex
  pairs of `GaussBonnet.lean`); no hypothesis beyond the orientation is needed.
* `eulerChar_toGluing`: **the same Euler characteristic**, `χ(gluing) = χ(TriMesh)`, provided that
  the corners at every vertex form **one fan** (`FansConnected`) and every vertex lies in a face
  (`NoIsolatedVertex`) — Rust `Complex::check` tests both. These are genuine hypotheses:
  `TriMesh.Manifold` (every edge in at most two faces, every boundary vertex on two boundary
  edges) and `Oriented` do not exclude a *pinched* vertex, e.g. two tetrahedra sharing one vertex
  (`|V| = 7`, `|E| = 12`, `|F| = 8`, `χ = 3`), where the gluing has two vertices at the pinch
  (`χ = 4`, two spheres). The edge count needs no such hypothesis.
* `toGluing_closed`: a closed valid oriented mesh has no boundary halfedges, and
  `toGluing_isBoundary_iff`: in general the boundary halfedges are the boundary edges of
  `GaussBonnet.lean`.
* `isMetric_toGluing`, `coneAngle_toGluing`: with symmetric lengths the mesh is a metric, and the
  cone angle of the class of a vertex is `TriMesh.angleSum` of the vertex; so the Gauss–Bonnet
  theorem of the gluing data (`IsMetric.gauss_bonnet_closed`) and `TriMesh.gauss_bonnet_closed`
  are the same statement for such meshes.
-/

namespace DDG

open Real

namespace TriMesh

open Delaunay

variable {V F : Type*} (M : TriMesh V F) (hO : M.Oriented)

/-- The gluing data of an oriented mesh: the twin of `TriMesh.twin`.

[formula: it-complex-from-mesh] -/
noncomputable def toGluing : Gluing F := ⟨M.twin, M.twin_twin hO⟩

@[simp] theorem toGluing_twin (h : Halfedge F) : (M.toGluing hO).twin h = M.twin h := rfl

/-- The tail of a halfedge (Rust `faces[f][k]`). -/
def tail (h : Halfedge F) : V := (M.halfedge h).1

theorem tail_next (h : Halfedge F) : M.tail h.next = (M.halfedge h).2 := rfl

theorem hasTwin_of_twin_ne {h : Halfedge F} (hb : M.twin h ≠ h) : M.HasTwin h := by
  by_contra ht
  exact hb (M.twin_of_not_hasTwin ht)

/-- The tails are a vertex labelling of the gluing data. -/
theorem toGluing_isLabeling : (M.toGluing hO).IsLabeling M.tail := by
  rintro x y ⟨hb, rfl⟩
  rw [tail_next, toGluing_twin, M.halfedge_twin (M.hasTwin_of_twin_ne hb)]
  rfl

/-- The corners at every vertex form one fan: equal tails are the same vertex class (Rust
`Complex::check`, "the corners of every vertex form one fan"). -/
def FansConnected : Prop := (M.toGluing hO).Faithful M.tail

/-- Every vertex is the tail of a halfedge (lies in a face). -/
def NoIsolatedVertex : Prop := Function.Surjective M.tail

/-- The vertex count of the gluing data is `|V|`. -/
theorem numVertices_toGluing (hfan : M.FansConnected hO) (hiso : M.NoIsolatedVertex) :
    (M.toGluing hO).numVertices = Nat.card V :=
  (M.toGluing_isLabeling hO).numVertices_eq hfan hiso

variable [DecidableEq V] [Fintype F]

theorem mem_edges_halfedge (h : Halfedge F) :
    s((M.halfedge h).1, (M.halfedge h).2) ∈ M.edges := by
  simp only [edges, Finset.mem_biUnion, Finset.mem_univ, true_and]
  exact ⟨h.1, M.halfedge_mem_faceEdges h⟩

/-- The unordered edge of a halfedge. -/
def edgeOfHalfedge (h : Halfedge F) : M.edges :=
  ⟨s((M.halfedge h).1, (M.halfedge h).2), M.mem_edges_halfedge h⟩

theorem edgeOfHalfedge_eq_iff (x y : Halfedge F) :
    M.edgeOfHalfedge x = M.edgeOfHalfedge y ↔ (M.toGluing hO).edgeSetoid x y := by
  rw [Gluing.edgeSetoid_r, toGluing_twin, edgeOfHalfedge, edgeOfHalfedge, Subtype.mk.injEq,
    Sym2.eq_iff]
  constructor
  · rintro (⟨h1, h2⟩ | ⟨h1, h2⟩)
    · exact Or.inl (hO (Prod.ext h1 h2)).symm
    · right
      have ht : M.HasTwin x := ⟨y, Prod.ext h2.symm h1.symm⟩
      apply hO
      rw [M.halfedge_twin ht]
      exact Prod.ext h2.symm h1.symm
  · rintro (rfl | rfl)
    · exact Or.inl ⟨rfl, rfl⟩
    · by_cases ht : M.HasTwin x
      · rw [M.halfedge_twin ht]
        exact Or.inr ⟨rfl, rfl⟩
      · rw [M.twin_of_not_hasTwin ht]
        exact Or.inl ⟨rfl, rfl⟩

/-- **The edges of the gluing data are the edges of the mesh.** -/
noncomputable def edgeEquivToGluing : (M.toGluing hO).Edge ≃ M.edges :=
  Equiv.ofBijective
    (Quotient.lift M.edgeOfHalfedge fun x y h => (M.edgeOfHalfedge_eq_iff hO x y).2 h)
    ⟨by
      intro p q
      induction p using Quotient.inductionOn
      induction q using Quotient.inductionOn
      intro h
      exact Quotient.sound ((M.edgeOfHalfedge_eq_iff hO _ _).1 h), by
      rintro ⟨e, he⟩
      simp only [edges, Finset.mem_biUnion, Finset.mem_univ, true_and] at he
      obtain ⟨f, hf⟩ := he
      induction e using Sym2.inductionOn with
      | hf x y =>
        obtain ⟨k, hk | hk⟩ := (M.mem_faceEdges_iff f x y).1 hf
        · refine ⟨(M.toGluing hO).edgeOf (f, k), ?_⟩
          simp only [Gluing.edgeOf, Quotient.lift_mk, edgeOfHalfedge, hk]
        · refine ⟨(M.toGluing hO).edgeOf (f, k), ?_⟩
          simp only [Gluing.edgeOf, Quotient.lift_mk, edgeOfHalfedge, hk, Subtype.mk.injEq]
          exact Sym2.eq_swap⟩

theorem numEdges_toGluing : (M.toGluing hO).numEdges = M.edges.card := by
  rw [Gluing.numEdges, Nat.card_congr (M.edgeEquivToGluing hO), Nat.card_eq_fintype_card,
    Fintype.card_coe]

/-- **An oriented mesh whose vertex stars are single fans has the Euler characteristic of its
gluing data.**

[formula: it-complex-from-mesh] [formula: euler-characteristic] -/
theorem eulerChar_toGluing [Fintype V] (hfan : M.FansConnected hO)
    (hiso : M.NoIsolatedVertex) : (M.toGluing hO).eulerChar = M.eulerChar := by
  rw [Gluing.eulerChar, eulerChar, M.numVertices_toGluing hO hfan hiso, M.numEdges_toGluing hO,
    Nat.card_eq_fintype_card]

/-- The boundary halfedges of the gluing data are the boundary edges of the mesh. -/
theorem toGluing_isBoundary_iff (hM : M.Valid) (h : Halfedge F) :
    (M.toGluing hO).twin h = h ↔ s((M.halfedge h).1, (M.halfedge h).2) ∈ M.boundaryEdges := by
  rw [← M.not_hasTwin_iff hM hO h, toGluing_twin]
  constructor
  · intro e ht
    exact M.twin_ne hM ht e
  · exact M.twin_of_not_hasTwin

/-- A closed valid oriented mesh gives closed gluing data. -/
theorem toGluing_closed [Fintype V] (hM : M.Valid)
    (hclosed : M.boundaryVertices = ∅) : (M.toGluing hO).Closed := fun h =>
  M.twin_ne hM (M.hasTwin_of_closed hM hO hclosed h)

/-! ## Lengths and angles -/

omit [DecidableEq V] [Fintype F] in
/-- The edge lengths of the mesh on its halfedges. -/
def halfedgeLength (h : Halfedge F) : ℝ := M.len (M.halfedge h).1 (M.halfedge h).2

omit [DecidableEq V] [Fintype F] in
/-- With symmetric lengths a valid oriented mesh is a metric on its gluing data. -/
theorem isMetric_toGluing (hM : M.Valid) (hsymm : ∀ x y, M.len x y = M.len y x) :
    (M.toGluing hO).IsMetric M.halfedgeLength := by
  refine ⟨fun h => ?_, fun f => ?_⟩
  · rw [toGluing_twin]
    by_cases ht : M.HasTwin h
    · rw [halfedgeLength, halfedgeLength, M.halfedge_twin ht]
      exact hsymm _ _
    · rw [M.twin_of_not_hasTwin ht]
  · exact (hM.tri f).rotate.rotate

omit [DecidableEq V] [Fintype F] in
theorem cornerAngleAt_toGluing (f : F) :
    Gluing.cornerAngleAt M.halfedgeLength (f, 0) = M.angleA f ∧
      Gluing.cornerAngleAt M.halfedgeLength (f, 1) = M.angleB f ∧
      Gluing.cornerAngleAt M.halfedgeLength (f, 2) = M.angleC f := by
  refine ⟨?_, ?_, ?_⟩ <;>
  · simp only [Gluing.cornerAngleAt, halfedgeLength, Halfedge.next, Halfedge.prev, angleA,
      angleB, angleC, la, lb, lc]
    exact Intrinsic.cornerAngle_swap _ _ _

omit [Fintype F] in
theorem angleAt_eq_sum_corners (hM : M.Valid) (f : F) (v : V) :
    ((if M.tail (f, 0) = v then M.angleA f else 0) + (if M.tail (f, 1) = v then M.angleB f else 0) +
      (if M.tail (f, 2) = v then M.angleC f else 0)) = M.angleAt f v := by
  have hab := hM.ab f
  have hbc := hM.bc f
  have hca := hM.ca f
  have hba := (hM.ab f).symm
  have hcb := (hM.bc f).symm
  have hac := (hM.ca f).symm
  simp only [tail, halfedge_zero, halfedge_one, halfedge_two, angleAt, eq_comm (a := v)]
  by_cases h1 : M.a f = v <;> by_cases h2 : M.b f = v <;> by_cases h3 : M.c f = v <;> simp_all

/-- **The cone angle of the class of a vertex is its angle sum** `TriMesh.angleSum`, on a valid
oriented mesh with connected vertex fans. So `Σ_v (2π − Θ_v)` of the gluing data is
`Σ_v angleDefect v` of the mesh, and the two Gauss–Bonnet theorems agree. -/
theorem coneAngle_toGluing (hM : M.Valid) (hfan : M.FansConnected hO)
    (h : Halfedge F) :
    (M.toGluing hO).coneAngle M.halfedgeLength ((M.toGluing hO).vertexOf h) =
      M.angleSum (M.tail h) := by
  have hl := M.toGluing_isLabeling hO
  have key : ∀ x, (M.toGluing hO).vertexOf x = (M.toGluing hO).vertexOf h ↔ M.tail x = M.tail h :=
    fun x => ⟨hl.eq_of_vertexOf, hfan x h⟩
  rw [Gluing.coneAngle, Finset.filter_congr (fun x _ => key x), Finset.sum_filter, angleSum,
    ← Finset.univ_product_univ, Finset.sum_product]
  refine Finset.sum_congr rfl fun f _ => ?_
  obtain ⟨c0, c1, c2⟩ := M.cornerAngleAt_toGluing f
  rw [Fin.sum_univ_three, c0, c1, c2]
  exact M.angleAt_eq_sum_corners hM f (M.tail h)

end TriMesh

/-- The consistently oriented tetrahedron of `Orientation.lean` as gluing data is closed. -/
example : (exampleTetrahedronOriented.toGluing (by intro h h' e; revert h h'; decide)).Closed :=
  exampleTetrahedronOriented.toGluing_closed _
    ⟨by decide, by decide, by decide, fun f => by
      simp only [TriMesh.la, TriMesh.lb, TriMesh.lc, exampleTetrahedronOriented]
      exact ⟨by norm_num, by norm_num, by norm_num⟩⟩ (by decide)

end DDG
