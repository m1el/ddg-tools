import DDG.Mesh.GluingFlip
import DDG.Delaunay.Local

/-!
# Edge-length metrics on gluing data, cone angles, and the flip

A **piecewise flat metric** on a triangulated surface given by gluing data
(`DDG/Mesh/Gluing.lean`) is a length for every edge such that every face satisfies the
strict triangle inequalities (Bobenko–Springborn 2007, Sec. 2; Fisher–Springborn–Bobenko–Schröder
2007, Sec. 2; Sharp–Soliman–Crane 2019, Sec. 2.1): `IsMetric ℓ`, with `ℓ : Halfedge F → ℝ` equal on
twins. Every face is then a Euclidean triangle and the surface is flat away from the vertices.

* `cornerAngleAt ℓ h`: the angle of the corner `h` (at the tail of `h`, between the sides `h` and
  `prev h`, opposite `next h`), by the law of cosines (`cornerAngle` of `Triangle.lean`; Rust
  `corner_angle`, the per-corner angles of `IntrinsicTriangulation`).
* `coneAngle ℓ v`: the sum of the corner angles at the vertex `v` (Rust `angle_sums`, the
  "`Θ`" of the signpost paper Sec. 3.1). The angle defect is `2π − Θ`.
* `sum_coneAngle`: the corner angles add up to `π·|F|`; `gauss_bonnet_closed`: **discrete
  Gauss–Bonnet** `Σ_v (2π − Θ_v) = 2π χ` on a closed surface, for any gluing data (self-loops,
  multiple edges, faces glued to themselves). It ties to `TriMesh.gauss_bonnet_closed` of
  `GaussBonnet.lean` through `DDG/Delaunay/FromTriMesh.lean`.

**The flip with lengths** (`flipLen`): the new edge `kl` gets `flipLength a b c d e` of
`DDG/Delaunay/Local.lean` (layer 1; equal to the oracle's and the Rust
`flippedLengthFromLengths`, `flipLength_eq_flippedLengthFromLengths`), with the diamond lengths
`a = ℓ(h₀)`, `b = ℓ(prev h₀)`, `c = ℓ(next h₀)`, `d = ℓ(next (twin h₀))`, `e = ℓ(prev (twin h₀))`
(`diamondLengths` reads them off); every other edge keeps its length (Rust `edge_lengths[e] =
new_length`, all other entries unchanged; edge ids are kept, halfedge ids move by `σ`).

* `IsMetric.flip`: if the diamond is strictly convex (`IsConvexDiamond`, layer 1), the result is
  again a metric; in particular after the flip of a non-Delaunay edge
  (`IsMetric.flip_of_not_isLocallyDelaunay`), whose new edge is locally Delaunay
  (`flip_edge_isLocallyDelaunay`).
* `coneAngle_flip`: **a flip keeps the cone angle of every vertex** of the surface, through the
  vertex bijection of `Flip.lean`. The one-diamond angle additivity is layer 1's
  `flip_corner_angles`; this file supplies the bookkeeping over the whole surface (which corners
  are at which vertex, also when `i, j, k, l` coincide).
-/

namespace DDG.Delaunay

open Real

namespace Gluing

variable {F : Type*} (G : Gluing F)

/-- `ℓ` is a **piecewise flat metric**: twins have the same length, and the three sides of every
face satisfy the strict triangle inequalities. -/
structure IsMetric (ℓ : Halfedge F → ℝ) : Prop where
  twin : ∀ h, ℓ (G.twin h) = ℓ h
  tri : ∀ f, IsTriangle (ℓ (f, 0)) (ℓ (f, 1)) (ℓ (f, 2))

omit G in
/-- The triangle inequalities of a face read from any of its corners. -/
theorem isTriangle_face_iff (ℓ : Halfedge F → ℝ) (h : Halfedge F) :
    IsTriangle (ℓ (h.1, 0)) (ℓ (h.1, 1)) (ℓ (h.1, 2)) ↔ IsTriangle (ℓ h) (ℓ h.next) (ℓ h.prev) := by
  obtain ⟨f, k⟩ := h
  fin_cases k
  · exact Iff.rfl
  · exact ⟨fun t => t.rotate, fun t => t.rotate.rotate⟩
  · exact ⟨fun t => t.rotate.rotate, fun t => t.rotate⟩

variable {G} in
theorem IsMetric.tri_at {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (h : Halfedge F) :
    IsTriangle (ℓ h) (ℓ h.next) (ℓ h.prev) :=
  (isTriangle_face_iff ℓ h).1 (hm.tri h.1)

omit G in
/-- The angle of the corner `h` of a face: at the tail of `h`, between the sides `h` and
`prev h`, opposite the side `next h`. -/
noncomputable def cornerAngleAt (ℓ : Halfedge F → ℝ) (h : Halfedge F) : ℝ :=
  cornerAngle (ℓ h.next) (ℓ h) (ℓ h.prev)

/-- The **cone angle** `Θ_v`: the sum of the corner angles at the vertex `v`. -/
noncomputable def coneAngle [Fintype F] (ℓ : Halfedge F → ℝ) (v : G.Vertex) : ℝ :=
  ∑ h ∈ Finset.univ.filter (fun h => G.vertexOf h = v), cornerAngleAt ℓ h

variable {G} in
/-- The three corner angles of a face add up to `π`. -/
theorem IsMetric.sum_cornerAngleAt_face {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (f : F) :
    ∑ k : Fin 3, cornerAngleAt ℓ (f, k) = π := by
  rw [Fin.sum_univ_three]
  have t := hm.tri f
  have := cornerAngle_add_cornerAngle_add_cornerAngle t
  simp only [cornerAngleAt, Halfedge.next, Halfedge.prev]
  have e0 : ((0 : Fin 3) + 1) = 1 := rfl
  have e1 : ((1 : Fin 3) + 1) = 2 := rfl
  have e2 : ((2 : Fin 3) + 1) = 0 := rfl
  have e3 : ((0 : Fin 3) + 2) = 2 := rfl
  have e4 : ((1 : Fin 3) + 2) = 0 := rfl
  have e5 : ((2 : Fin 3) + 2) = 1 := rfl
  rw [e0, e1, e2, e3, e4, e5, Intrinsic.cornerAngle_swap (ℓ (f, 1)),
    Intrinsic.cornerAngle_swap (ℓ (f, 2)), Intrinsic.cornerAngle_swap (ℓ (f, 0))]
  linarith

variable {G} in
/-- **The corner angles add up to `π·|F|`**, so the cone angles do. -/
theorem IsMetric.sum_coneAngle [Fintype F] {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) :
    ∑ v, G.coneAngle ℓ v = π * Fintype.card F := by
  unfold coneAngle
  rw [Finset.sum_fiberwise Finset.univ G.vertexOf (cornerAngleAt ℓ), ← Finset.univ_product_univ,
    Finset.sum_product]
  simp only [hm.sum_cornerAngleAt_face, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  ring

variable {G} in
/-- **Discrete Gauss–Bonnet** on a closed surface given by gluing data:
`Σ_v (2π − Θ_v) = 2π χ`, `χ = V − E + F`.

[formula: gauss-bonnet] [formula: it-complex-euler] -/
theorem IsMetric.gauss_bonnet_closed [Fintype F] {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ)
    (hc : G.Closed) : ∑ v, (2 * π - G.coneAngle ℓ v) = 2 * π * G.eulerChar := by
  have hE := G.two_mul_numEdges_of_closed hc
  have hE' : (2 : ℝ) * G.numEdges = 3 * Fintype.card F := by exact_mod_cast hE
  rw [Finset.sum_sub_distrib, hm.sum_coneAngle, Finset.sum_const, Finset.card_univ, eulerChar,
    numVertices, Nat.card_eq_fintype_card, nsmul_eq_mul]
  push_cast
  linear_combination π * hE'

/-! ## The flip with lengths -/

/-- The five lengths `(a, b, c, d, e)` of the diamond of `h₀` in the notation of layer 1
(`DDG/Delaunay/Local.lean`): `a = ℓ_ij`, `b = ℓ_ik`, `c = ℓ_jk`, `d = ℓ_il`, `e = ℓ_jl`. -/
def diamondLengths (ℓ : Halfedge F → ℝ) (h₀ : Halfedge F) : ℝ × ℝ × ℝ × ℝ × ℝ :=
  (ℓ h₀, ℓ h₀.prev, ℓ h₀.next, ℓ (G.twin h₀).next, ℓ (G.twin h₀).prev)

/-- The length of the flipped edge: `flipLength a b c d e` of the diamond. -/
noncomputable def newLength (ℓ : Halfedge F → ℝ) (h₀ : Halfedge F) : ℝ :=
  flipLength (ℓ h₀) (ℓ h₀.prev) (ℓ h₀.next) (ℓ (G.twin h₀).next) (ℓ (G.twin h₀).prev)

/-- The diamond of `h₀` is strictly convex (layer 1's `IsConvexDiamond` of its five lengths). -/
def ConvexAt (ℓ : Halfedge F → ℝ) (h₀ : Halfedge F) : Prop :=
  IsConvexDiamond (ℓ h₀) (ℓ h₀.prev) (ℓ h₀.next) (ℓ (G.twin h₀).next) (ℓ (G.twin h₀).prev)

/-- The edge of `h₀` is locally Delaunay (layer 1's `IsLocallyDelaunay` of its five lengths). -/
def LocallyDelaunayAt (ℓ : Halfedge F → ℝ) (h₀ : Halfedge F) : Prop :=
  IsLocallyDelaunay (ℓ h₀) (ℓ h₀.prev) (ℓ h₀.next) (ℓ (G.twin h₀).next) (ℓ (G.twin h₀).prev)

variable [DecidableEq F]

/-- The lengths after the flip (Rust `IntrinsicTriangulation::flip_edge`): the new edge gets
`newLength`, every surviving halfedge keeps its length under its new id. -/
noncomputable def flipLen (ℓ : Halfedge F → ℝ) (h₀ : Halfedge F) : Halfedge F → ℝ :=
  fun x => if x = h₀.next ∨ x = (G.twin h₀).next then G.newLength ℓ h₀
    else ℓ ((G.flipPerm h₀).symm x)

variable {G} {h₀ : Halfedge F}

section FlipLen

variable (hf : G.Flippable h₀) (ℓ : Halfedge F → ℝ)
include hf

/-- Every halfedge keeps its length under its new id; the flipped edge gets the new length. -/
theorem flipLen_flipPerm (u : Halfedge F) :
    G.flipLen ℓ h₀ (G.flipPerm h₀ u) =
      if u = h₀ ∨ u = G.twin h₀ then G.newLength ℓ h₀ else ℓ u := by
  have e1 : G.flipPerm h₀ u = h₀.next ↔ u = h₀ := by
    rw [← flipPerm_self hf, (G.flipPerm h₀).injective.eq_iff]
  have e2 : G.flipPerm h₀ u = (G.twin h₀).next ↔ u = G.twin h₀ := by
    rw [← flipPerm_twin hf, (G.flipPerm h₀).injective.eq_iff]
  simp only [flipLen, e1, e2, Equiv.symm_apply_apply]

omit hf in
theorem flipLen_next : G.flipLen ℓ h₀ h₀.next = G.newLength ℓ h₀ := by
  simp [flipLen]

omit hf in
theorem flipLen_twin_next : G.flipLen ℓ h₀ (G.twin h₀).next = G.newLength ℓ h₀ := by
  simp [flipLen]

theorem flipLen_self : G.flipLen ℓ h₀ h₀ = ℓ (G.twin h₀).next := by
  have := flipLen_flipPerm hf ℓ (G.twin h₀).next
  rw [flipPerm_twin_next hf] at this
  rw [this, ite_eq_right_of_eq_false]
  simp only [eq_iff_iff, iff_false, not_or]
  exact ⟨hf.ne' (G.twin h₀).next h₀, Halfedge.next_ne _⟩

theorem flipLen_twin : G.flipLen ℓ h₀ (G.twin h₀) = ℓ h₀.next := by
  have := flipLen_flipPerm hf ℓ h₀.next
  rw [flipPerm_next hf] at this
  rw [this, ite_eq_right_of_eq_false]
  simp only [eq_iff_iff, iff_false, not_or]
  exact ⟨Halfedge.next_ne _, hf.ne h₀.next (G.twin h₀)⟩

theorem flipLen_of_ne {x : Halfedge F} (x0 : x ≠ h₀) (x1 : x ≠ h₀.next) (x2 : x ≠ G.twin h₀)
    (x3 : x ≠ (G.twin h₀).next) : G.flipLen ℓ h₀ x = ℓ x := by
  conv_lhs => rw [← flipPerm_of_ne x0 x1 x2 x3]
  rw [flipLen_flipPerm hf, ite_eq_right_of_eq_false]
  simp only [eq_iff_iff, iff_false, not_or]
  exact ⟨x0, x2⟩

theorem flipLen_prev : G.flipLen ℓ h₀ h₀.prev = ℓ h₀.prev :=
  flipLen_of_ne hf ℓ (Halfedge.prev_ne _) (Halfedge.prev_ne_next _) (hf.ne h₀.prev (G.twin h₀))
    (hf.ne h₀.prev (G.twin h₀).next)

theorem flipLen_twin_prev : G.flipLen ℓ h₀ (G.twin h₀).prev = ℓ (G.twin h₀).prev :=
  flipLen_of_ne hf ℓ (hf.ne' (G.twin h₀).prev h₀) (hf.ne' (G.twin h₀).prev h₀.next)
    (Halfedge.prev_ne _) (Halfedge.prev_ne_next _)

/-- Outside the two faces of the diamond nothing changes. -/
theorem flipLen_of_face {x : Halfedge F} (x0 : x.1 ≠ h₀.1) (x1 : x.1 ≠ (G.twin h₀).1) :
    G.flipLen ℓ h₀ x = ℓ x :=
  flipLen_of_ne hf ℓ (fun e => x0 (by rw [e])) (fun e => x0 (by rw [e]; rfl))
    (fun e => x1 (by rw [e])) (fun e => x1 (by rw [e]; rfl))

end FlipLen

omit [DecidableEq F] in
/-- The two triangles of the diamond: `ijk` with sides `(a, b, c)`, `jil` with `(a, d, e)`. -/
theorem IsMetric.diamond_isTriangle {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (h₀ : Halfedge F) :
    IsTriangle (ℓ h₀) (ℓ h₀.prev) (ℓ h₀.next) ∧
      IsTriangle (ℓ h₀) (ℓ (G.twin h₀).next) (ℓ (G.twin h₀).prev) := by
  refine ⟨(hm.tri_at h₀).swap, ?_⟩
  have := hm.tri_at (G.twin h₀)
  rwa [hm.twin] at this

/-- **The flip of a strictly convex diamond gives a metric again**: the new faces `ilk` and `jkl`
satisfy the strict triangle inequalities (layer 1's `flipLength_isTriangle`), and the lengths
are equal on the new twins.

[formula: it-flip-length] -/
theorem IsMetric.flip {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (hf : G.Flippable h₀)
    (hconv : G.ConvexAt ℓ h₀) : (G.flip h₀).IsMetric (G.flipLen ℓ h₀) := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  obtain ⟨t1, t2⟩ := flipLength_isTriangle hk hl hconv
  refine ⟨fun x => ?_, fun g => ?_⟩
  · obtain ⟨u, rfl⟩ := (G.flipPerm h₀).surjective x
    rw [flip_twin_apply, flipLen_flipPerm hf, flipLen_flipPerm hf, hm.twin]
    congr 1
    apply propext
    constructor
    · rintro (e | e)
      · exact Or.inr (by rw [← e, G.twin_twin])
      · exact Or.inl (G.twin_injective e)
    · rintro (e | e)
      · exact Or.inr (by rw [e])
      · exact Or.inl (by rw [e, G.twin_twin])
  · by_cases g0 : g = h₀.1
    · subst g0
      refine (isTriangle_face_iff _ h₀).2 ?_
      rw [flipLen_self hf, flipLen_next, flipLen_prev hf]
      exact t1.rotate.rotate
    by_cases g1 : g = (G.twin h₀).1
    · subst g1
      refine (isTriangle_face_iff _ (G.twin h₀)).2 ?_
      rw [flipLen_twin hf, flipLen_twin_next, flipLen_twin_prev hf]
      exact t2.rotate.rotate
    rw [flipLen_of_face hf ℓ g0 g1, flipLen_of_face hf ℓ g0 g1, flipLen_of_face hf ℓ g0 g1]
    exact hm.tri g

/-- **Every non-Delaunay edge can be flipped**, and the result is a metric (B–S Prop. 11, through
layer 1's `isConvexDiamond_of_not_isLocallyDelaunay`). -/
theorem IsMetric.flip_of_not_locallyDelaunay {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ)
    (hf : G.Flippable h₀) (hnd : ¬G.LocallyDelaunayAt ℓ h₀) :
    (G.flip h₀).IsMetric (G.flipLen ℓ h₀) := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  exact hm.flip hf (isConvexDiamond_of_not_isLocallyDelaunay hk hl hnd)

omit [DecidableEq F] in
/-- The local Delaunay condition is the angle condition at the two corners opposite the edge:
`θ_k + θ_l ≤ π` with `θ_k` the corner `prev h₀`, `θ_l` the corner `prev (twin h₀)`. -/
theorem IsMetric.locallyDelaunayAt_iff {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (h₀ : Halfedge F) :
    G.LocallyDelaunayAt ℓ h₀ ↔
      cornerAngleAt ℓ h₀.prev + cornerAngleAt ℓ (G.twin h₀).prev ≤ π := by
  simp only [LocallyDelaunayAt, IsLocallyDelaunay, cornerAngleAt, Halfedge.next_prev,
    Halfedge.prev_prev, hm.twin]
  rw [Intrinsic.cornerAngle_swap (ℓ h₀) (ℓ (G.twin h₀).prev)]

/-- **After the flip of a non-Delaunay edge the new edge is locally Delaunay**, strictly
(layer 1's `flip_isLocallyDelaunay`): the two corners opposite the new edge `next h₀ = l → k`,
at `i` in `ilk` and at `j` in `jkl`, sum to less than `π`. -/
theorem IsMetric.flip_edge_isLocallyDelaunay {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ)
    (hf : G.Flippable h₀) (hnd : ¬G.LocallyDelaunayAt ℓ h₀) :
    cornerAngleAt (G.flipLen ℓ h₀) h₀.next.prev +
        cornerAngleAt (G.flipLen ℓ h₀) ((G.flip h₀).twin h₀.next).prev < π := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  have := flip_isLocallyDelaunay hk hl hnd
  rw [flip_twin_next hf]
  simp only [cornerAngleAt, Halfedge.prev_next]
  rw [flipLen_self hf, flipLen_next, flipLen_prev hf, flipLen_twin hf, flipLen_twin_next,
    flipLen_twin_prev hf,
    Intrinsic.cornerAngle_swap (G.newLength ℓ h₀) (ℓ h₀.prev) (ℓ (G.twin h₀).next)]
  unfold newLength
  linarith

/-! ## Cone angles are preserved -/

omit [DecidableEq F] in
/-- A sum over all halfedges of a function that vanishes outside the two faces of the diamond
is the sum of its six values there. -/
theorem Flippable.sum_eq_six [Fintype F] [DecidableEq F] (hf : G.Flippable h₀)
    (g : Halfedge F → ℝ) (hg : ∀ y : Halfedge F, y.1 ≠ h₀.1 → y.1 ≠ (G.twin h₀).1 → g y = 0) :
    ∑ y, g y = g h₀ + g h₀.next + g h₀.prev + g (G.twin h₀) + g (G.twin h₀).next +
      g (G.twin h₀).prev := by
  set D : Finset (Halfedge F) := {h₀, h₀.next, h₀.prev, G.twin h₀, (G.twin h₀).next,
    (G.twin h₀).prev}
  have hsub := Finset.sum_subset (Finset.subset_univ D) (f := g) (by
    intro y _ hy
    apply hg
    · intro e
      rcases Halfedge.eq_or_eq_next_or_eq_prev e with rfl | rfl | rfl <;> simp [D] at hy
    · intro e
      rcases Halfedge.eq_or_eq_next_or_eq_prev e with rfl | rfl | rfl <;> simp [D] at hy)
  rw [← hsub]
  have c1 := hf.ne h₀ (G.twin h₀); have c2 := hf.ne h₀ (G.twin h₀).next
  have c3 := hf.ne h₀ (G.twin h₀).prev; have c4 := hf.ne h₀.next (G.twin h₀)
  have c5 := hf.ne h₀.next (G.twin h₀).next; have c6 := hf.ne h₀.next (G.twin h₀).prev
  have c7 := hf.ne h₀.prev (G.twin h₀); have c8 := hf.ne h₀.prev (G.twin h₀).next
  have c9 := hf.ne h₀.prev (G.twin h₀).prev
  rw [Finset.sum_insert (by simp [c1, c2, c3]), Finset.sum_insert (by simp [c4, c5, c6]),
    Finset.sum_insert (by simp [c7, c8, c9]), Finset.sum_insert (by simp),
    Finset.sum_insert (by simp), Finset.sum_singleton]
  ring

/-- **A flip keeps the cone angle of every vertex**: `Θ'(φ v) = Θ(v)` for the vertex bijection
`φ = flipVertexEquiv`, for every strictly convex diamond. The four vertices of the diamond may
coincide in any way. The corner bookkeeping is here; the angle additivity in one diamond is layer
1's `flip_corner_angles` (`DDG/Delaunay/Local.lean`).

[formula: it-angle-normalization] -/
theorem IsMetric.coneAngle_flip [Fintype F] {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ)
    (hf : G.Flippable h₀) (hconv : G.ConvexAt ℓ h₀) (v : G.Vertex) :
    (G.flip h₀).coneAngle (G.flipLen ℓ h₀) (flipVertexEquiv hf v) = G.coneAngle ℓ v := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  obtain ⟨A1, A2, A3, A4⟩ := flip_corner_angles hk hl hconv
  unfold coneAngle
  rw [Finset.sum_filter, Finset.sum_filter, ← sub_eq_zero, ← Finset.sum_sub_distrib]
  have hφ : ∀ x, (G.flip h₀).vertexOf (G.flipCorner h₀ x) = flipVertexEquiv hf v ↔
      G.vertexOf x = v := fun x => by
    rw [← flipVertexEquiv_vertexOf hf x]; exact (flipVertexEquiv hf).injective.eq_iff
  rw [hf.sum_eq_six]
  · -- the six corners of the diamond
    obtain ⟨d1, d2, d3, d4, d5, d6⟩ := flip_diamond_vertices hf
    have hi : (G.flip h₀).vertexOf h₀ = flipVertexEquiv hf v ↔ G.vertexOf h₀ = v := by
      rw [← d1]; exact (flipVertexEquiv hf).injective.eq_iff
    have hj : (G.flip h₀).vertexOf (G.twin h₀) = flipVertexEquiv hf v ↔
        G.vertexOf (G.twin h₀) = v := by
      rw [← d2]; exact (flipVertexEquiv hf).injective.eq_iff
    have hk' : (G.flip h₀).vertexOf h₀.prev = flipVertexEquiv hf v ↔
        G.vertexOf h₀.prev = v := by
      rw [← d3]; exact (flipVertexEquiv hf).injective.eq_iff
    have hl' : (G.flip h₀).vertexOf (G.twin h₀).prev = flipVertexEquiv hf v ↔
        G.vertexOf (G.twin h₀).prev = v := by
      rw [← d4]; exact (flipVertexEquiv hf).injective.eq_iff
    have hl2 : (G.flip h₀).vertexOf h₀.next = flipVertexEquiv hf v ↔
        G.vertexOf (G.twin h₀).prev = v := by
      rw [d5]; exact (flipVertexEquiv hf).injective.eq_iff
    have hk2 : (G.flip h₀).vertexOf (G.twin h₀).next = flipVertexEquiv hf v ↔
        G.vertexOf h₀.prev = v := by
      rw [(G.flip h₀).vertexOf_eq_of_glued (flip_glued_k hf)]; exact hk'
    have oi : G.vertexOf (G.twin h₀).next = G.vertexOf h₀ :=
      (G.vertexOf_eq_of_glued (glued_i hf)).symm
    have oj : G.vertexOf h₀.next = G.vertexOf (G.twin h₀) :=
      (G.vertexOf_eq_of_glued (glued_j hf)).symm
    simp only [hi, hj, hk', hl', hl2, hk2, oi, oj]
    simp only [cornerAngleAt, Halfedge.next_next, Halfedge.prev_next, Halfedge.next_prev,
      Halfedge.prev_prev]
    rw [flipLen_self hf, flipLen_next, flipLen_prev hf, flipLen_twin hf, flipLen_twin_next,
      flipLen_twin_prev hf, hm.twin]
    have s1 := Intrinsic.cornerAngle_swap (G.newLength ℓ h₀) (ℓ h₀.prev) (ℓ (G.twin h₀).next)
    have s2 := Intrinsic.cornerAngle_swap (ℓ h₀.prev) (ℓ (G.twin h₀).next) (G.newLength ℓ h₀)
    have s3 := Intrinsic.cornerAngle_swap (ℓ (G.twin h₀).next) (G.newLength ℓ h₀) (ℓ h₀.prev)
    have s4 := Intrinsic.cornerAngle_swap (ℓ (G.twin h₀).next) (ℓ (G.twin h₀).prev) (ℓ h₀)
    have s5 := Intrinsic.cornerAngle_swap (ℓ (G.twin h₀).prev) (ℓ h₀) (ℓ (G.twin h₀).next)
    have s6 := Intrinsic.cornerAngle_swap (ℓ h₀) (ℓ (G.twin h₀).next) (ℓ (G.twin h₀).prev)
    simp only [newLength] at *
    split_ifs <;> linarith
  · -- outside the diamond nothing changes
    intro y y0 y1
    have yn : y.next.1 ≠ h₀.1 := y0
    have yn' : y.next.1 ≠ (G.twin h₀).1 := y1
    have yp : y.prev.1 ≠ h₀.1 := y0
    have yp' : y.prev.1 ≠ (G.twin h₀).1 := y1
    have hc : G.flipCorner h₀ y = y :=
      flipCorner_of_ne (fun e => y1 (by rw [e]; rfl)) (fun e => y0 (by rw [e]; rfl))
    have := hφ y
    rw [hc] at this
    simp only [this, cornerAngleAt, flipLen_of_face hf ℓ y0 y1, flipLen_of_face hf ℓ yn yn',
      flipLen_of_face hf ℓ yp yp', sub_self]

end Gluing

end DDG.Delaunay
