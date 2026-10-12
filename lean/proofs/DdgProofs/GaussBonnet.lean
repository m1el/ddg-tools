import DdgProofs.Cotan

/-!
# Discrete curvature and discrete Gauss–Bonnet

For a triangle mesh with intrinsic edge lengths (`TriMesh`):

* `angleDefect v = 2π − Σ(corner angles at v)` (paper Eq. 14), for interior vertices;
* `geodesicCurvature v = π − Σ(corner angles at v)` (paper Eq. 15), for boundary vertices;
* `gauss_bonnet`: `Σ_{interior} Ω + Σ_{boundary} k = 2π χ`, `χ = |V| − |E| + |F|`.

The combinatorial hypotheses are explicit (`TriMesh.Manifold`):

1. every edge lies in at most two faces;
2. every boundary vertex lies on exactly two boundary edges;

where an edge is an unordered pair of vertices of a face, a *boundary edge* is an edge lying in
exactly one face, and a *boundary vertex* is an end point of a boundary edge. Together with
`TriMesh.Valid` (three distinct vertices and the strict triangle inequalities in every face) these
give the two counting identities `2|E| = 3|F| + |E_∂|` and `|E_∂| = |B|` that the proof needs.
Orientability is not needed.
-/

noncomputable section

namespace DdgProofs

open Real Matrix

namespace TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype V] [Fintype F] (M : TriMesh V F)

/-! ## Angles and curvature -/

/-- Interior angle of the face `f` at the vertex `v`; `0` if `v` is not a vertex of `f`. -/
def angleAt (f : F) (v : V) : ℝ :=
  if v = M.a f then M.angleA f
  else if v = M.b f then M.angleB f
  else if v = M.c f then M.angleC f
  else 0

/-- Sum of the corner angles at the vertex `v`. -/
def angleSum (v : V) : ℝ := ∑ f, M.angleAt f v

/-- Angle defect `Ω_v = 2π − Σ(corner angles at v)` (paper Eq. 14). -/
def angleDefect (v : V) : ℝ := 2 * π - M.angleSum v

/-- Geodesic curvature `k_v = π − Σ(corner angles at v)` of a boundary vertex (paper Eq. 15). -/
def geodesicCurvature (v : V) : ℝ := π - M.angleSum v

/-- The corner angles of one face add up to `π`, so all corner angles add up to `π·|F|`. -/
theorem sum_angleSum (hM : M.Valid) : ∑ v, M.angleSum v = π * (Fintype.card F : ℝ) := by
  unfold angleSum
  rw [Finset.sum_comm]
  have hface : ∀ f, ∑ v, M.angleAt f v = π := by
    intro f
    have hab := hM.ab f
    have hbc := hM.bc f
    have hca := hM.ca f
    have h3 : ∀ v, M.angleAt f v = (if v = M.a f then M.angleA f else 0) +
        (if v = M.b f then M.angleB f else 0) + (if v = M.c f then M.angleC f else 0) := by
      intro v
      unfold angleAt
      by_cases h1 : v = M.a f <;> by_cases h2 : v = M.b f <;> by_cases h3 : v = M.c f <;>
        simp_all
    simp only [h3, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
    exact cornerAngle_add_cornerAngle_add_cornerAngle (hM.tri f)
  simp only [hface, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  ring

omit [Fintype V] in
/-- The angle sums, hence the angle defects, are invariant under uniform scaling.

[formula: scale-invariance] -/
theorem angleDefect_scale {s : ℝ} (hs : s ≠ 0) (v : V) :
    (M.scale s).angleDefect v = M.angleDefect v := by
  simp only [angleDefect, angleSum, angleAt, M.angleA_scale hs, M.angleB_scale hs,
    M.angleC_scale hs]
  rfl

/-! ## Combinatorics -/

/-- The three edges of a face, as unordered pairs of vertices. -/
def faceEdges (f : F) : Finset (Sym2 V) :=
  {s(M.a f, M.b f), s(M.b f, M.c f), s(M.c f, M.a f)}

/-- The set of edges of the mesh. -/
def edges : Finset (Sym2 V) := Finset.univ.biUnion M.faceEdges

/-- Number of faces containing the edge `e`. -/
def edgeDegree (e : Sym2 V) : ℕ := (Finset.univ.filter fun f => e ∈ M.faceEdges f).card

/-- Boundary edges: the edges that lie in exactly one face. -/
def boundaryEdges : Finset (Sym2 V) := M.edges.filter fun e => M.edgeDegree e = 1

/-- Boundary vertices: the end points of boundary edges. All other vertices are interior. -/
def boundaryVertices : Finset V := Finset.univ.filter fun v => ∃ e ∈ M.boundaryEdges, v ∈ e

/-- Euler characteristic `χ = |V| − |E| + |F|`. -/
def eulerChar : ℤ := (Fintype.card V : ℤ) - (M.edges.card : ℤ) + (Fintype.card F : ℤ)

/-- The combinatorial manifold conditions used by discrete Gauss–Bonnet: every edge lies in at
most two faces, and every boundary vertex lies on exactly two boundary edges. -/
structure Manifold : Prop where
  edge_le_two : ∀ e ∈ M.edges, M.edgeDegree e ≤ 2
  boundary_two : ∀ v ∈ M.boundaryVertices, (M.boundaryEdges.filter fun e => v ∈ e).card = 2

omit [Fintype V] [Fintype F] in
theorem card_faceEdges (hM : M.Valid) (f : F) : (M.faceEdges f).card = 3 := by
  unfold faceEdges
  rw [Finset.card_eq_three]
  refine ⟨_, _, _, ?_, ?_, ?_, rfl⟩
  · intro h
    rw [Sym2.eq_iff] at h
    rcases h with ⟨h1, _⟩ | ⟨h1, _⟩
    · exact hM.ab f h1
    · exact hM.ca f h1.symm
  · intro h
    rw [Sym2.eq_iff] at h
    rcases h with ⟨h1, _⟩ | ⟨_, h2⟩
    · exact hM.ca f h1.symm
    · exact hM.bc f h2
  · intro h
    rw [Sym2.eq_iff] at h
    rcases h with ⟨h1, _⟩ | ⟨h1, _⟩
    · exact hM.bc f h1
    · exact hM.ab f h1.symm

omit [Fintype V] in
theorem exists_ne_of_mem_edges (hM : M.Valid) {e : Sym2 V} (he : e ∈ M.edges) :
    ∃ x y, x ≠ y ∧ e = s(x, y) := by
  simp only [edges, Finset.mem_biUnion, Finset.mem_univ, true_and, faceEdges,
    Finset.mem_insert, Finset.mem_singleton] at he
  obtain ⟨f, rfl | rfl | rfl⟩ := he
  · exact ⟨_, _, hM.ab f, rfl⟩
  · exact ⟨_, _, hM.bc f, rfl⟩
  · exact ⟨_, _, hM.ca f, rfl⟩

omit [Fintype V] in
theorem edgeDegree_pos {e : Sym2 V} (he : e ∈ M.edges) : 1 ≤ M.edgeDegree e := by
  simp only [edges, Finset.mem_biUnion, Finset.mem_univ, true_and] at he
  obtain ⟨f, hf⟩ := he
  exact Finset.card_pos.mpr ⟨f, by simp [hf]⟩

omit [Fintype V] in
/-- Counting incidences of edges and faces: `Σ_e (number of faces at e) = 3|F|`. -/
theorem sum_edgeDegree (hM : M.Valid) : ∑ e ∈ M.edges, M.edgeDegree e = 3 * Fintype.card F := by
  have h := Finset.sum_card_bipartiteAbove_eq_sum_card_bipartiteBelow (s := M.edges)
    (t := (Finset.univ : Finset F)) (fun e f => e ∈ M.faceEdges f)
  have h2 : ∀ f, M.edges.bipartiteBelow (fun e f => e ∈ M.faceEdges f) f = M.faceEdges f := by
    intro f
    ext e
    simp only [Finset.mem_bipartiteBelow, edges, Finset.mem_biUnion, Finset.mem_univ, true_and]
    exact ⟨fun h => h.2, fun h => ⟨⟨f, h⟩, h⟩⟩
  simp only [h2, M.card_faceEdges hM, Finset.sum_const, Finset.card_univ, smul_eq_mul] at h
  rw [mul_comm]
  exact h

/-- First counting identity of a manifold triangle mesh: `2|E| = 3|F| + |E_∂|`. -/
theorem two_mul_card_edges (hM : M.Valid) (hman : M.Manifold) :
    2 * M.edges.card = 3 * Fintype.card F + M.boundaryEdges.card := by
  have h1 := M.sum_edgeDegree hM
  have h2 : ∀ e ∈ M.edges, M.edgeDegree e + (if M.edgeDegree e = 1 then 1 else 0) = 2 := by
    intro e he
    have := M.edgeDegree_pos he
    have := hman.edge_le_two e he
    split_ifs with h <;> omega
  have h3 : ∑ e ∈ M.edges, (M.edgeDegree e + (if M.edgeDegree e = 1 then 1 else 0)) =
      2 * M.edges.card := by
    rw [Finset.sum_congr rfl h2, Finset.sum_const, smul_eq_mul, mul_comm]
  rw [Finset.sum_add_distrib, h1, ← Finset.card_filter] at h3
  unfold boundaryEdges
  omega

/-- Second counting identity: the boundary has as many edges as vertices, `|E_∂| = |B|`. -/
theorem card_boundaryEdges_eq (hM : M.Valid) (hman : M.Manifold) :
    M.boundaryEdges.card = M.boundaryVertices.card := by
  have h := Finset.sum_card_bipartiteAbove_eq_sum_card_bipartiteBelow (s := M.boundaryEdges)
    (t := M.boundaryVertices) (fun e v => v ∈ e)
  have hL : ∀ e ∈ M.boundaryEdges,
      (M.boundaryVertices.bipartiteAbove (fun e v => v ∈ e) e).card = 2 := by
    intro e he
    obtain ⟨x, y, hxy, rfl⟩ := M.exists_ne_of_mem_edges hM (Finset.mem_filter.mp he).1
    have hx : x ∈ M.boundaryVertices := by
      simp only [boundaryVertices, Finset.mem_filter, Finset.mem_univ, true_and]
      exact ⟨_, he, Sym2.mem_mk_left x y⟩
    have hy : y ∈ M.boundaryVertices := by
      simp only [boundaryVertices, Finset.mem_filter, Finset.mem_univ, true_and]
      exact ⟨_, he, Sym2.mem_mk_right x y⟩
    have hset : M.boundaryVertices.bipartiteAbove (fun e v => v ∈ e) s(x, y) = {x, y} := by
      ext v
      simp only [Finset.mem_bipartiteAbove, Sym2.mem_iff, Finset.mem_insert,
        Finset.mem_singleton]
      constructor
      · exact fun h => h.2
      · rintro (rfl | rfl)
        · exact ⟨hx, Or.inl rfl⟩
        · exact ⟨hy, Or.inr rfl⟩
    rw [hset, Finset.card_pair hxy]
  have hR : ∀ v ∈ M.boundaryVertices,
      (M.boundaryEdges.bipartiteBelow (fun e v => v ∈ e) v).card = 2 :=
    fun v hv => hman.boundary_two v hv
  rw [Finset.sum_congr rfl hL, Finset.sum_congr rfl hR] at h
  simp only [Finset.sum_const, smul_eq_mul] at h
  omega

/-! ## Gauss–Bonnet -/

/-- **Discrete Gauss–Bonnet.** For a valid manifold triangle mesh,
`Σ_{v interior} Ω_v + Σ_{v boundary} k_v = 2π·χ`, with `χ = |V| − |E| + |F|`.

The angle defect `Ω_v = 2π − Σθ`, the boundary curvature `k_v = π − Σθ` and the Euler
characteristic are definitions (`angleDefect`, `geodesicCurvature`, `eulerChar`); this theorem is
the statement that ties the three together, so it also carries their ids.

[formula: gauss-bonnet] [formula: angle-defect] [formula: boundary-geodesic-curvature]
[formula: euler-characteristic] -/
theorem gauss_bonnet (hM : M.Valid) (hman : M.Manifold) :
    ∑ v ∈ M.boundaryVerticesᶜ, M.angleDefect v + ∑ v ∈ M.boundaryVertices, M.geodesicCurvature v =
      2 * π * (M.eulerChar : ℝ) := by
  have hsum := M.sum_angleSum hM
  have hsplit : ∑ v ∈ M.boundaryVerticesᶜ, M.angleSum v + ∑ v ∈ M.boundaryVertices, M.angleSum v =
      ∑ v, M.angleSum v := by
    rw [add_comm]
    exact Finset.sum_add_sum_compl _ _
  have hE := M.two_mul_card_edges hM hman
  have hB := M.card_boundaryEdges_eq hM hman
  have hE' : (2 : ℝ) * M.edges.card = 3 * Fintype.card F + M.boundaryVertices.card := by
    rw [← hB]
    exact_mod_cast hE
  have hcard : (M.boundaryVerticesᶜ.card : ℝ) =
      (Fintype.card V : ℝ) - (M.boundaryVertices.card : ℝ) := by
    rw [Finset.card_compl, Nat.cast_sub (Finset.card_le_univ _)]
  simp only [angleDefect, geodesicCurvature, Finset.sum_sub_distrib, Finset.sum_const,
    nsmul_eq_mul, eulerChar]
  push_cast
  rw [hcard]
  linear_combination (-1 : ℝ) * hsplit + (-1 : ℝ) * hsum + π * hE'

/-- **Gauss–Bonnet on a closed mesh**: with no boundary, `Σ_v Ω_v = 2π·χ`.

[formula: gauss-bonnet] -/
theorem gauss_bonnet_closed (hM : M.Valid) (hman : M.Manifold)
    (hclosed : M.boundaryVertices = ∅) : ∑ v, M.angleDefect v = 2 * π * (M.eulerChar : ℝ) := by
  have h := M.gauss_bonnet hM hman
  rw [hclosed] at h
  simpa using h

/-- **Total cone angle on a closed mesh.** If `L u = μ − Ω` on a closed valid manifold mesh
(`L` the cotan Laplacian, `Ω` the angle defects) then `Σ_v μ_v = 2π·χ`.

[formula: closed-mesh-cone-sum] -/
theorem sum_cone_eq_two_pi_eulerChar (hM : M.Valid) (hman : M.Manifold)
    (hclosed : M.boundaryVertices = ∅) {u μ : V → ℝ}
    (h : M.laplacian *ᵥ u = μ - M.angleDefect) : ∑ v, μ v = 2 * π * (M.eulerChar : ℝ) := by
  have h1 : (fun _ => (1 : ℝ)) ⬝ᵥ (M.laplacian *ᵥ u) = 0 := by
    rw [Matrix.dotProduct_mulVec, M.laplacian_vecMul_one]
    simp
  rw [h] at h1
  simp only [dotProduct, Pi.sub_apply, one_mul, Finset.sum_sub_distrib] at h1
  rw [← M.gauss_bonnet_closed hM hman hclosed]
  linarith

end TriMesh

/-! ## The hypotheses are satisfiable -/

/-- A single equilateral triangle: a mesh with boundary. -/
def exampleTriangle : TriMesh (Fin 3) (Fin 1) :=
  ⟨fun _ => 0, fun _ => 1, fun _ => 2, fun _ _ => 1⟩

/-- A regular tetrahedron: a closed mesh. -/
def exampleTetrahedron : TriMesh (Fin 4) (Fin 4) :=
  ⟨![0, 0, 0, 1], ![1, 2, 3, 2], ![2, 3, 1, 3], fun _ _ => 1⟩

/-- The single triangle satisfies `Valid` and `Manifold`, and all its vertices are boundary
vertices. -/
example : exampleTriangle.Valid ∧ exampleTriangle.Manifold ∧
    exampleTriangle.boundaryVertices = Finset.univ := by
  refine ⟨⟨by decide, by decide, by decide, fun f => ?_⟩, ⟨by decide, by decide⟩, by decide⟩
  simp only [TriMesh.la, TriMesh.lb, TriMesh.lc, exampleTriangle]
  exact ⟨by norm_num, by norm_num, by norm_num⟩

/-- The tetrahedron satisfies `Valid` and `Manifold`, and has no boundary. -/
example : exampleTetrahedron.Valid ∧ exampleTetrahedron.Manifold ∧
    exampleTetrahedron.boundaryVertices = ∅ := by
  refine ⟨⟨by decide, by decide, by decide, fun f => ?_⟩, ⟨by decide, by decide⟩, by decide⟩
  simp only [TriMesh.la, TriMesh.lb, TriMesh.lc, exampleTetrahedron]
  exact ⟨by norm_num, by norm_num, by norm_num⟩

end DdgProofs
