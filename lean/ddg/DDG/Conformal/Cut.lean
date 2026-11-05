import DDG.Conformal.Cherrier
import DDG.Laplace.Cotan

/-!
# Cutting a mesh: bookkeeping, and the angle sums at cones

`crates/ddg-bff/src/cut.rs` cuts a mesh along edges: the faces stay, every vertex on the cut is
replaced by one copy per sector, and each face corner is re-indexed to its copy. A per-vertex
array is lifted through the map `ρ : (vertices of the cut mesh) → (original vertices)`
(`CutMesh::vertex_to_original`), and the edge lengths stay those of the original mesh.

`IsCutOf M' M ρ` is exactly that description: same faces, corners compatible with `ρ`, lengths
pulled back. Under it:

* per-face quantities (side lengths, cotangents, angles, areas) are equal;
* `sum_fiber_angleSum`, `sum_fiber_laplacian`: the angle sums, and `L'(u∘ρ)`, of the copies of a
  vertex add up to the angle sum, resp. `(L u)`, of the original vertex;
* **`bff-cone-angle-sum`** (BFF Proposition 2), `sum_fiber_interior_angle`: with the exterior
  angles `k̃_p = k'_p + (L'(u∘ρ))_p` of the cut surface, `Σ_{p ↦ v} (π − k̃_p) = 2π − (Ω_v + (L u)_v)`;
  if `u` solves `(L u)_v = μ_v − Ω_v` this is `2π − μ_v`: the interior angles of the layout at the
  copies of a cone add up to its cone angle sum;
* `curvature_copy_eq`: at a vertex that was not split the first-order curvature
  `Ω + L u` is unchanged, so the cut surface is flat away from the cut when `(L u)_v = −Ω_v` there,
  and Gauss–Bonnet on the cut surface (`TriMesh.sum_exteriorAngle` applied to `M'`) gives
  `Σ k̃ = 2πχ(M')`;
* **`cut-euler-characteristic`**, `cut_euler_characteristic`:
  `χ(M') = χ(M) + Σ_v (s_v − 1) − Σ_e (c_e − 1)`, with `s_v` the number of copies of the vertex `v`
  and `c_e` the number of copies of the edge `e`; `card_edge_copies_le`: `c_e` is at most the number
  of faces at `e`, so on a manifold mesh `c_e − 1 ∈ {0, 1}` and the last sum is the number of cut
  edges.

Not proved here: that a chart with `χ = 1` and one boundary loop is a disk (classification of
surfaces, cited), and that the planner's cuts produce such charts (cited; checked at run time).
-/

noncomputable section

namespace DDG

open Matrix

namespace TriMesh

variable {V V' F : Type*}

/-- `M'` is `M` cut open: the same faces, with corners re-indexed to copies of the vertices
(`ρ` maps a copy to its original), and the edge lengths of `M`. -/
structure IsCutOf (M' : TriMesh V' F) (M : TriMesh V F) (ρ : V' → V) : Prop where
  a : ∀ f, ρ (M'.a f) = M.a f
  b : ∀ f, ρ (M'.b f) = M.b f
  c : ∀ f, ρ (M'.c f) = M.c f
  len : ∀ x y, M'.len x y = M.len (ρ x) (ρ y)

namespace IsCutOf

variable {M' : TriMesh V' F} {M : TriMesh V F} {ρ : V' → V}

theorem la (h : IsCutOf M' M ρ) (f : F) : M'.la f = M.la f := by
  simp only [TriMesh.la, h.len, h.b, h.c]

theorem lb (h : IsCutOf M' M ρ) (f : F) : M'.lb f = M.lb f := by
  simp only [TriMesh.lb, h.len, h.c, h.a]

theorem lc (h : IsCutOf M' M ρ) (f : F) : M'.lc f = M.lc f := by
  simp only [TriMesh.lc, h.len, h.a, h.b]

theorem cotA (h : IsCutOf M' M ρ) (f : F) : M'.cotA f = M.cotA f := by
  simp only [TriMesh.cotA, h.la, h.lb, h.lc]

theorem cotB (h : IsCutOf M' M ρ) (f : F) : M'.cotB f = M.cotB f := by
  simp only [TriMesh.cotB, h.la, h.lb, h.lc]

theorem cotC (h : IsCutOf M' M ρ) (f : F) : M'.cotC f = M.cotC f := by
  simp only [TriMesh.cotC, h.la, h.lb, h.lc]

theorem angleA (h : IsCutOf M' M ρ) (f : F) : M'.angleA f = M.angleA f := by
  simp only [TriMesh.angleA, h.la, h.lb, h.lc]

theorem angleB (h : IsCutOf M' M ρ) (f : F) : M'.angleB f = M.angleB f := by
  simp only [TriMesh.angleB, h.la, h.lb, h.lc]

theorem angleC (h : IsCutOf M' M ρ) (f : F) : M'.angleC f = M.angleC f := by
  simp only [TriMesh.angleC, h.la, h.lb, h.lc]

/-- Cutting does not change the area of a face (so not the total area either). -/
theorem area (h : IsCutOf M' M ρ) (f : F) : M'.area f = M.area f := by
  simp only [TriMesh.area, h.la, h.lb, h.lc]

/-- A cut of a valid mesh is valid. -/
theorem valid (h : IsCutOf M' M ρ) (hM : M.Valid) : M'.Valid where
  ab f hab := hM.ab f (by rw [← h.a, ← h.b, hab])
  bc f hbc := hM.bc f (by rw [← h.b, ← h.c, hbc])
  ca f hca := hM.ca f (by rw [← h.c, ← h.a, hca])
  tri f := by rw [h.la, h.lb, h.lc]; exact hM.tri f

end IsCutOf

variable [DecidableEq V] [DecidableEq V'] [Fintype V] [Fintype V'] [Fintype F]
variable {M' : TriMesh V' F} {M : TriMesh V F} {ρ : V' → V}

omit [Fintype V] in
/-- Summing an indicator of one copy over the copies of `v`. -/
private theorem sum_fiber_ite (ρ : V' → V) (v : V) (x : V') (r : ℝ) :
    ∑ p ∈ Finset.univ.filter (fun p => ρ p = v), (if p = x then r else 0) =
      if ρ x = v then r else 0 := by
  rw [Finset.sum_ite_eq']
  simp

omit [Fintype V] [Fintype F] in
/-- With three distinct corners, the angle of a face at a vertex is a sum of three indicator
terms. -/
theorem angleAt_eq_indicator (M : TriMesh V F) (hM : M.Valid) (f : F) (v : V) :
    M.angleAt f v = (if v = M.a f then M.angleA f else 0) + (if v = M.b f then M.angleB f else 0) +
      (if v = M.c f then M.angleC f else 0) := by
  have hab := hM.ab f
  have hbc := hM.bc f
  have hca := hM.ca f
  unfold angleAt
  by_cases h1 : v = M.a f <;> by_cases h2 : v = M.b f <;> by_cases h3 : v = M.c f <;> simp_all

omit [Fintype V] [Fintype F] in
/-- The angles of a face at the copies of `v` add up to its angle at `v`. -/
theorem sum_fiber_angleAt (h : IsCutOf M' M ρ) (hM : M.Valid) (f : F) (v : V) :
    ∑ p ∈ Finset.univ.filter (fun p => ρ p = v), M'.angleAt f p = M.angleAt f v := by
  simp only [M'.angleAt_eq_indicator (h.valid hM), M.angleAt_eq_indicator hM,
    Finset.sum_add_distrib, sum_fiber_ite, h.a, h.b, h.c, h.angleA, h.angleB, h.angleC]
  simp only [eq_comm]

omit [Fintype V] in
/-- **The angle sums of the copies of a vertex add up to its angle sum.** -/
theorem sum_fiber_angleSum (h : IsCutOf M' M ρ) (hM : M.Valid) (v : V) :
    ∑ p ∈ Finset.univ.filter (fun p => ρ p = v), M'.angleSum p = M.angleSum v := by
  simp only [angleSum]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun f _ => sum_fiber_angleAt h hM f v

/-- Row `i` of the cotan Laplacian as a sum over faces. -/
theorem laplacian_mulVec_eq_sum_faces (M : TriMesh V F) (u : V → ℝ) (i : V) :
    (M.laplacian *ᵥ u) i = 1 / 2 * ∑ f, ∑ j, M.cotOpp f i j * (u i - u j) := by
  rw [laplacian_mulVec_apply, Finset.sum_comm, Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  simp only [offDiag]
  rw [← Finset.sum_mul]
  ring

/-- **`L'(u∘ρ)` summed over the copies of a vertex is `(L u)` at the vertex**: the cotan
Laplacian is a sum of per-face contributions at the corners, and the corners of a vertex are
distributed over its copies. -/
theorem sum_fiber_laplacian (h : IsCutOf M' M ρ) (hM : M.Valid) (u : V → ℝ) (v : V) :
    ∑ p ∈ Finset.univ.filter (fun p => ρ p = v), (M'.laplacian *ᵥ (u ∘ ρ)) p =
      (M.laplacian *ᵥ u) v := by
  have hM' := h.valid hM
  simp only [laplacian_mulVec_eq_sum_faces, ← Finset.mul_sum]
  congr 1
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun f _ => ?_
  have e1 : ∀ p : V', ∑ q, M'.cotOpp f p q * ((u ∘ ρ) p - (u ∘ ρ) q) =
      (if p = M'.a f then M.cotC f * (u (M.a f) - u (M.b f)) +
        M.cotB f * (u (M.a f) - u (M.c f)) else 0) +
      (if p = M'.b f then M.cotC f * (u (M.b f) - u (M.a f)) +
        M.cotA f * (u (M.b f) - u (M.c f)) else 0) +
      (if p = M'.c f then M.cotA f * (u (M.c f) - u (M.b f)) +
        M.cotB f * (u (M.c f) - u (M.a f)) else 0) := by
    intro p
    rw [M'.sum_cotOpp_mul hM' f p (fun q => (u ∘ ρ) p - (u ∘ ρ) q)]
    simp only [Function.comp, h.a, h.b, h.c, h.cotA, h.cotB, h.cotC]
    rcases eq_or_ne p (M'.a f) with h1 | h1
    · subst h1
      have hb : M'.a f ≠ M'.b f := hM'.ab f
      have hc : M'.a f ≠ M'.c f := (hM'.ca f).symm
      simp only [hb, hc, h.a, ite_true, ite_false]
      ring
    · rcases eq_or_ne p (M'.b f) with h2 | h2
      · subst h2
        have hc : M'.b f ≠ M'.c f := hM'.bc f
        simp only [h1, hc, h.b, ite_true, ite_false]
        ring
      · rcases eq_or_ne p (M'.c f) with h3 | h3
        · subst h3
          simp only [h1, h2, h.c, ite_true, ite_false]
          ring
        · simp only [h1, h2, h3, ite_false]
          ring
  have e2 : ∑ j, M.cotOpp f v j * (u v - u j) =
      (if v = M.a f then M.cotC f * (u (M.a f) - u (M.b f)) +
        M.cotB f * (u (M.a f) - u (M.c f)) else 0) +
      (if v = M.b f then M.cotC f * (u (M.b f) - u (M.a f)) +
        M.cotA f * (u (M.b f) - u (M.c f)) else 0) +
      (if v = M.c f then M.cotA f * (u (M.c f) - u (M.b f)) +
        M.cotB f * (u (M.c f) - u (M.a f)) else 0) := by
    rw [M.sum_cotOpp_mul hM f v (fun j => u v - u j)]
    rcases eq_or_ne v (M.a f) with h1 | h1
    · subst h1
      have hb : M.a f ≠ M.b f := hM.ab f
      have hc : M.a f ≠ M.c f := (hM.ca f).symm
      simp only [hb, hc, ite_true, ite_false]
      ring
    · rcases eq_or_ne v (M.b f) with h2 | h2
      · subst h2
        have hc : M.b f ≠ M.c f := hM.bc f
        simp only [h1, hc, ite_true, ite_false]
        ring
      · rcases eq_or_ne v (M.c f) with h3 | h3
        · subst h3
          simp only [h1, h2, ite_true, ite_false]
          ring
        · simp only [h1, h2, h3, ite_false]
          ring
  simp only [e1, e2, Finset.sum_add_distrib, sum_fiber_ite, h.a, h.b, h.c]
  simp only [eq_comm]

/-- **The angle sums at a cone** (BFF Proposition 2). Let `M'` be a cut of the valid mesh `M`
and `u` scale factors on `M`, lifted to the cut mesh as `u ∘ ρ`. With the exterior angles of the
cut surface `k̃_p = k'_p + (L'(u∘ρ))_p` (`k'_p = π − Σθ` over the corners of the copy `p`), the
interior angles `π − k̃_p` at the copies of a vertex `v` add up to

`Σ_{p : ρ p = v} (π − k̃_p) = 2π − (Ω_v + (L u)_v)`.

[formula: bff-cone-angle-sum] -/
theorem sum_fiber_interior_angle (h : IsCutOf M' M ρ) (hM : M.Valid) (u : V → ℝ) (v : V) :
    ∑ p ∈ Finset.univ.filter (fun p => ρ p = v),
        (Real.pi - (M'.geodesicCurvature p + (M'.laplacian *ᵥ (u ∘ ρ)) p)) =
      2 * Real.pi - (M.angleDefect v + (M.laplacian *ᵥ u) v) := by
  have h1 := sum_fiber_angleSum h hM v
  have h2 := sum_fiber_laplacian h hM u v
  have e : ∀ p : V', Real.pi - (M'.geodesicCurvature p + (M'.laplacian *ᵥ (u ∘ ρ)) p) =
      M'.angleSum p - (M'.laplacian *ᵥ (u ∘ ρ)) p := by
    intro p; simp only [geodesicCurvature]; ring
  simp only [e, Finset.sum_sub_distrib, h1, h2, angleDefect]
  ring

/-- **… equal the cone's angle sum `2π − μ_v`** when `u` solves the scale-factor equation
`(L u)_v = μ_v − Ω_v` at `v` (a vertex without cone has `μ_v = 0`: its copies' interior angles
add up to `2π`).

[formula: bff-cone-angle-sum] -/
theorem sum_fiber_interior_angle_cone (h : IsCutOf M' M ρ) (hM : M.Valid) {u μ : V → ℝ} (v : V)
    (hu : (M.laplacian *ᵥ u) v = μ v - M.angleDefect v) :
    ∑ p ∈ Finset.univ.filter (fun p => ρ p = v),
        (Real.pi - (M'.geodesicCurvature p + (M'.laplacian *ᵥ (u ∘ ρ)) p)) =
      2 * Real.pi - μ v := by
  rw [sum_fiber_interior_angle h hM u v, hu]; ring

/-- **Away from the cut nothing changes**: at a copy `p` that is the only copy of its vertex,
the first-order curvature of the cut surface equals that of the original surface,
`Ω'_p + (L'(u∘ρ))_p = Ω_v + (L u)_v` (`v = ρ p`). So if `(L u)_v = −Ω_v` there, the cut surface is
flat at `p` to first order, which is the hypothesis of `TriMesh.sum_exteriorAngle` for `M'`. -/
theorem curvature_copy_eq (h : IsCutOf M' M ρ) (hM : M.Valid) (u : V → ℝ) (p : V')
    (hp : Finset.univ.filter (fun q => ρ q = ρ p) = {p}) :
    M'.angleDefect p + (M'.laplacian *ᵥ (u ∘ ρ)) p =
      M.angleDefect (ρ p) + (M.laplacian *ᵥ u) (ρ p) := by
  have h1 := sum_fiber_angleSum h hM (ρ p)
  have h2 := sum_fiber_laplacian h hM u (ρ p)
  rw [hp, Finset.sum_singleton] at h1 h2
  simp only [angleDefect]
  linarith

/-! ## The Euler characteristic after cutting -/

omit [DecidableEq V'] [Fintype V] [Fintype V'] [Fintype F] in
/-- The three edges of a face of the cut mesh project to the three edges of the face. -/
theorem image_faceEdges [DecidableEq V'] (h : IsCutOf M' M ρ) (f : F) :
    (M'.faceEdges f).image (Sym2.map ρ) = M.faceEdges f := by
  simp [faceEdges, Finset.image_insert, h.a, h.b, h.c]

omit [Fintype V] [Fintype V'] in
/-- The edges of the cut mesh project onto the edges of the mesh. -/
theorem image_edges (h : IsCutOf M' M ρ) : M'.edges.image (Sym2.map ρ) = M.edges := by
  unfold edges
  rw [Finset.biUnion_image]
  exact Finset.biUnion_congr rfl fun f _ => image_faceEdges h f

omit [DecidableEq V] [DecidableEq V'] [Fintype F] in
/-- The number of vertices added by cutting: `|V'| − |V| = Σ_v (s_v − 1)`, `s_v` the number of
copies of `v` (every vertex has at least one copy). -/
theorem card_vertices_cut [DecidableEq V] (hsurj : Function.Surjective ρ) :
    (Fintype.card V' : ℤ) - Fintype.card V =
      ∑ v : V, (((Finset.univ.filter (fun p => ρ p = v)).card : ℤ) - 1) := by
  have h1 : (Finset.univ : Finset V').card =
      ∑ v ∈ (Finset.univ : Finset V').image ρ, (Finset.univ.filter (fun p => ρ p = v)).card :=
    Finset.card_eq_sum_card_image ρ Finset.univ
  have h2 : (Finset.univ : Finset V').image ρ = Finset.univ :=
    Finset.image_univ_of_surjective hsurj
  rw [h2] at h1
  rw [Finset.sum_sub_distrib]
  simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_one]
  have h3 : (Fintype.card V' : ℤ) =
      ∑ v : V, ((Finset.univ.filter (fun p => ρ p = v)).card : ℤ) := by
    rw [← Finset.card_univ, h1]; push_cast; rfl
  rw [h3]

omit [Fintype V] [Fintype V'] in
/-- The number of edges added by cutting: `|E'| − |E| = Σ_e (c_e − 1)`, `c_e` the number of copies
of the edge `e` in the cut mesh. -/
theorem card_edges_cut (h : IsCutOf M' M ρ) :
    (M'.edges.card : ℤ) - M.edges.card =
      ∑ e ∈ M.edges, (((M'.edges.filter (fun e' => Sym2.map ρ e' = e)).card : ℤ) - 1) := by
  have h1 : M'.edges.card =
      ∑ e ∈ M'.edges.image (Sym2.map ρ), (M'.edges.filter (fun e' => Sym2.map ρ e' = e)).card :=
    Finset.card_eq_sum_card_image (Sym2.map ρ) M'.edges
  rw [image_edges h] at h1
  rw [Finset.sum_sub_distrib]
  simp only [Finset.sum_const, nsmul_eq_mul, mul_one]
  rw [h1]; push_cast; rfl

/-- **The Euler characteristic after cutting**: `χ(M') = χ(M) + a − c` with
`a = Σ_v (s_v − 1)` the number of added vertex copies and `c = Σ_e (c_e − 1)` the number of added
edges (`cut_euler_characteristic` in `cut.rs`: faces are unchanged, every cut edge becomes two
edges — see `card_edge_copies_le` for "two").

[formula: cut-euler-characteristic] -/
theorem cut_euler_characteristic (h : IsCutOf M' M ρ) (hsurj : Function.Surjective ρ) :
    M'.eulerChar = M.eulerChar +
      ∑ v : V, (((Finset.univ.filter (fun p => ρ p = v)).card : ℤ) - 1) -
      ∑ e ∈ M.edges, (((M'.edges.filter (fun e' => Sym2.map ρ e' = e)).card : ℤ) - 1) := by
  rw [← card_vertices_cut hsurj, ← card_edges_cut h]
  simp only [eulerChar]
  ring

omit [Fintype V] [Fintype V'] in
/-- **An edge has at most as many copies as it has faces.** For a cut of a valid mesh, the
number of edges of the cut mesh that project to `e` is at most the number of faces of `M` at `e`
(each copy lies in a face at `e`, and a face contains only one copy). On a manifold mesh that
number is `1` for a boundary edge and at most `2` for an interior edge, so an edge is either
uncut (one copy) or cut (two copies). -/
theorem card_edge_copies_le (h : IsCutOf M' M ρ) (hM : M.Valid) (e : Sym2 V) :
    (M'.edges.filter (fun e' => Sym2.map ρ e' = e)).card ≤ M.edgeDegree e := by
  classical
  have hM' := h.valid hM
  rcases isEmpty_or_nonempty F with hF | hF
  · have : M'.edges = ∅ := by simp [edges, Finset.univ_eq_empty]
    simp [this]
  -- choose, for every copy, a face of the cut mesh that contains it
  have hex : ∀ e' ∈ M'.edges.filter (fun e' => Sym2.map ρ e' = e), ∃ f, e' ∈ M'.faceEdges f := by
    intro e' he'
    have := (Finset.mem_filter.mp he').1
    simp only [edges, Finset.mem_biUnion, Finset.mem_univ, true_and] at this
    exact this
  choose! face hface using hex
  unfold edgeDegree
  refine Finset.card_le_card_of_injOn face (fun e' he' => ?_) (fun e₁ h₁ e₂ h₂ heq => ?_)
  · have h1 := hface e' he'
    have h2 := (Finset.mem_filter.mp he').2
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq]
    rw [← image_faceEdges h, ← h2]
    exact Finset.mem_image_of_mem _ h1
  · -- two copies in the same face that project to the same edge are equal
    have m1 := hface e₁ h₁
    have m2 := hface e₂ h₂
    rw [heq] at m1
    have p1 := (Finset.mem_filter.mp (Finset.mem_coe.mp h₁)).2
    have p2 := (Finset.mem_filter.mp (Finset.mem_coe.mp h₂)).2
    have hinj : Set.InjOn (Sym2.map ρ) (M'.faceEdges (face e₂) : Set (Sym2 V')) := by
      apply Finset.injOn_of_card_image_eq
      rw [image_faceEdges h, card_faceEdges M hM, card_faceEdges M' hM']
    exact hinj m1 m2 (by rw [p1, p2])

end TriMesh

end DDG
