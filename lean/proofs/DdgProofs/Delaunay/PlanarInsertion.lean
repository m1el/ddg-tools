import DdgProofs.Delaunay.PlanarLemma

/-!
# Planar Delaunay triangulations, part 5: inserting a point outside the convex hull

Let `T` be a triangulation of the finite point set `S` (`IsTriangulation`,
`DdgProofs/Delaunay/PlanarFlip.lean`) and let `p` lie outside the convex hull of `S`. Coning `p` to
every boundary edge of `T` that is visible from `p` gives a triangulation of `S ∪ {p}`
(`isTriangulation_insert`). This is the insertion step of the incremental (sweep) construction of a
triangulation; `DdgProofs/Delaunay/PlanarExistence.lean` uses it to prove that every finite point
set that is not collinear has a triangulation.

## Visible edges

The edge `ab` of a triangle `{a, b, c} ∈ T` is **visible from `p`** (`Visible S T p a b c`) when
`p` is strictly on the other side of the line `ab` than `c`, and all of `S` (hence all of `conv S`,
`Visible.halfplane`) is on `c`'s closed side: `ab` lies on a supporting line of `conv S` that
separates `conv S` from `p`. The new triangles are the **cones** `{a, b, p}` over the visible edges
(`cones S T p`), and the new triangulation is `T ∪ cones S T p`.

## The proof

* **Disjoint interiors.** A cone lies strictly beyond its supporting line, `conv S` before it, so a
  cone and an old triangle have disjoint interiors (`disjoint_cone_old`). For two cones over `ab` and
  `a'b'`: a common interior point `z` lies on the segment from `p` to a point `q` of the open edge
  `ab` and on the segment from `p` to a point `q'` of the open edge `a'b'`; comparing the two
  supporting lines along the ray from `p` through `z` shows `q = q'` (each of `q, q'` is in `conv S`,
  so neither lies strictly between `p` and the other), and then adjacency
  (`eq_insert_of_mem_edge`, `DdgProofs/Delaunay/PlanarLemma.lean`) and the supporting line force
  `{a, b} = {a', b'}` (`disjoint_cone_cone`).
* **The vertex condition.** A point of `S` in a cone lies on the supporting line, hence on the edge
  `ab`, hence is `a` or `b` (`vertex_cone`).
* **A boundary edge is a supporting line** (`exists_across`): if a point of `conv S` lies strictly
  beyond the edge `ab` of `{a, b, c} ∈ T`, then `T` has a triangle `{a, b, d}` across `ab` (the
  segments from a relative interior point of `ab` to that point are covered by the finitely many
  closed triangles; adjacency identifies the one that contains the relative interior point).
* **The cover** (`cover_insert`). A point `x` of `conv(S ∪ {p}) \ conv S` lies on a segment from
  some `y ∈ conv S` to `p`. Walking from `y` to `x`, the walk leaves `conv S` for the last time
  through an edge `ab` of a triangle `{a, b, c}` (`exit_edge`); when the line `xp` contains no point
  of `S` (a *generic* `x`), the exit point is a relative interior point of `ab`, no triangle of `T` lies
  across `ab` (the walk would continue in it), so `ab` is visible from `p` and `x` lies in the cone
  `{a, b, p}` (`exists_cone_of_generic`). The generic points are dense (finitely many lines through
  `p` have empty interior), and the union of the triangles is closed, which covers the remaining
  points.

## References

* M. de Berg, O. Cheong, M. van Kreveld, M. Overmars, *Computational Geometry: Algorithms and
  Applications*, 3rd ed., Springer, 2008, Ch. 9.1 and Ch. 1 (incremental construction; the points
  outside the current hull see the hull edges they lie beyond).
* H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge University Press, 2001,
  Sec. 1.1 (incremental construction of a triangulation by sweeping the points in sorted order and
  coning each new point to the visible edges of the hull).

The bibliographic data were not checked against the works here; the statements are proved from
the definitions of `DdgProofs/Delaunay/PlanarFlip.lean`, not taken from the works.
-/

noncomputable section

namespace DDG

namespace Delaunay

namespace Planar

open MeasureTheory Set Filter Topology

/-! ## Coordinates along an edge -/

/-- `sda` on the line `ab` in two triangles `a b c` and `a b p` (first barycentric coordinate). -/
theorem sda_line_left (a b c p s : ℝ × ℝ) :
    sda s b c * sda a b p - sda s b p * sda a b c = sda p b c * sda a b s := by
  unfold sda; ring

/-- `sda` on the line `ab` in two triangles `a b c` and `a b p` (second barycentric coordinate). -/
theorem sda_line_right (a b c p s : ℝ × ℝ) :
    sda a s c * sda a b p - sda a s p * sda a b c = sda a p c * sda a b s := by
  unfold sda; ring

/-- For `s` on the line `ab`, the first barycentric coordinate of `s` has the same sign in the
triangles `a b c` and `a b d`. -/
theorem edge_transfer₁ {a b s : ℝ × ℝ} (hs : sda a b s = 0) (c d : ℝ × ℝ) :
    (sda a b d * sda s b d) * (sda a b c * sda a b c) =
      (sda a b d * sda a b d) * (sda a b c * sda s b c) := by
  have I := sda_line_left a b c d s
  linear_combination (-(sda a b d * sda a b c)) * I +
    (-(sda a b d * sda a b c * sda d b c)) * hs

/-- For `s` on the line `ab`, the second barycentric coordinate of `s` has the same sign in the
triangles `a b c` and `a b d`. -/
theorem edge_transfer₂ {a b s : ℝ × ℝ} (hs : sda a b s = 0) (c d : ℝ × ℝ) :
    (sda a b d * sda a s d) * (sda a b c * sda a b c) =
      (sda a b d * sda a b d) * (sda a b c * sda a s c) := by
  have I := sda_line_right a b c d s
  linear_combination (-(sda a b d * sda a b c)) * I +
    (-(sda a b d * sda a b c * sda a d c)) * hs

theorem pos_of_edge_transfer {X Y K M : ℝ} (h : X * (K * K) = (M * M) * Y) (hM : M ≠ 0)
    (hY : 0 < Y) : 0 < X := by
  have := mul_pos (mul_self_pos.2 hM) hY
  rw [← h] at this
  exact pos_of_mul_pos_left this (mul_self_nonneg _)

theorem nonneg_of_edge_transfer {X Y K M : ℝ} (h : X * (K * K) = (M * M) * Y) (hK : K ≠ 0)
    (hY : 0 ≤ Y) : 0 ≤ X := by
  have := mul_nonneg (mul_self_nonneg M) hY
  rw [← h] at this
  exact (mul_nonneg_iff_of_pos_right (mul_self_pos.2 hK)).1 this

/-- Three reals: `0 ≤ D Z`, `0 ≤ P Z` and `D P < 0` force `Z = 0`. -/
theorem eq_zero_of_signs {D P Z : ℝ} (h1 : 0 ≤ D * Z) (h2 : 0 ≤ P * Z) (h : D * P < 0) :
    Z = 0 := by
  by_contra hZ
  have hZ2 : 0 < Z ^ 2 := by positivity
  nlinarith [mul_nonneg h1 h2]

theorem div_pos_of_mul_pos' {P X : ℝ} (hP : P ≠ 0) (h : 0 < P * X) : 0 < X / P := by
  have : X / P = (P * X) / (P * P) := (mul_div_mul_left X P hP).symm
  rw [this]
  exact div_pos h (mul_self_pos.2 hP)

/-- A relative interior point of the edge `ab` of a triangle of `T` is not a point of `S`. -/
theorem not_mem_of_mem_edge {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {a b c q : ℝ × ℝ} (hu : {a, b, c} ∈ T)
    (hq1 : 0 < sda a b c * sda q b c) (hq2 : 0 < sda a b c * sda a q c) (hq3 : sda a b q = 0) :
    q ∉ S := by
  have hD := hT.sda_ne_zero hu
  intro hqS
  have hqu : q ∈ hull {a, b, c} := by
    rw [hull_triple]; exact (mem_hull3_iff hD).2 ⟨hq1.le, hq2.le, by rw [hq3, mul_zero]⟩
  have := hT.vertex _ hu q hqS hqu
  simp only [Finset.mem_insert, Finset.mem_singleton] at this
  rcases this with rfl | rfl | rfl
  · rw [sda_self₁₂, mul_zero] at hq2; exact lt_irrefl _ hq2
  · rw [sda_self₁₂, mul_zero] at hq1; exact lt_irrefl _ hq1
  · exact hD hq3

/-! ## A boundary edge lies on a supporting line -/

/-- **An edge with points of `conv S` beyond it has a neighbour across it.** Let `q` be a relative
interior point of the edge `ab` of `{a, b, c} ∈ T`, and let `x ∈ conv S` be strictly on the other
side of `ab` than `c`. Then `T` contains a triangle `{a, b, d}` with `d` strictly across `ab`. -/
theorem exists_across {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {a b c q x : ℝ × ℝ} (hu : {a, b, c} ∈ T)
    (hq1 : 0 < sda a b c * sda q b c) (hq2 : 0 < sda a b c * sda a q c) (hq3 : sda a b q = 0)
    (hx : x ∈ convexHull ℝ (S : Set (ℝ × ℝ))) (hxs : sda a b c * sda a b x < 0) :
    ∃ d, {a, b, d} ∈ T ∧ sda a b c * sda a b d < 0 := by
  have hD := hT.sda_ne_zero hu
  have hqS := not_mem_of_mem_edge hT hu hq1 hq2 hq3
  have hqK : q ∈ convexHull ℝ (S : Set (ℝ × ℝ)) := hT.hull_subset hu (by
    rw [hull_triple]; exact (mem_hull3_iff hD).2 ⟨hq1.le, hq2.le, by rw [hq3, mul_zero]⟩)
  obtain ⟨φ, hφ⟩ : ∃ φ : ℝ → ℝ × ℝ, ∀ ε, φ ε = (1 - ε) • q + ε • x := ⟨_, fun _ => rfl⟩
  have hφc : Continuous φ := by
    have : φ = fun ε => (1 - ε) • q + ε • x := funext hφ
    rw [this]; fun_prop
  have hφK : ∀ ε ∈ Ioc (0 : ℝ) 1, φ ε ∈ convexHull ℝ (S : Set (ℝ × ℝ)) := fun ε hε => by
    rw [hφ]; exact convex_convexHull ℝ _ hqK hx (by linarith [hε.2]) hε.1.le (by ring)
  have hφs : ∀ ε, sda a b (φ ε) = ε * sda a b x := fun ε => by
    rw [hφ, sda_combo₃ _ _ _ _ (by ring), hq3]; ring
  have hcov : Ioc (0 : ℝ) 1 ⊆ ⋃ u ∈ (T : Set (Finset (ℝ × ℝ))), φ ⁻¹' hull u ∩ Ioc 0 1 := by
    intro ε hε
    have := hφK ε hε
    rw [← hT.cover] at this
    obtain ⟨u, hu, hm⟩ := Set.mem_iUnion₂.1 this
    exact Set.mem_biUnion (Finset.mem_coe.2 hu) ⟨hm, hε⟩
  have hcl : (0 : ℝ) ∈ closure (⋃ u ∈ (T : Set (Finset (ℝ × ℝ))), φ ⁻¹' hull u ∩ Ioc 0 1) :=
    closure_mono hcov (by rw [closure_Ioc zero_ne_one]; exact ⟨le_rfl, zero_le_one⟩)
  rw [T.finite_toSet.closure_biUnion] at hcl
  obtain ⟨u, hu', hcu⟩ := Set.mem_iUnion₂.1 hcl
  have hne : (φ ⁻¹' hull u ∩ Ioc 0 1).Nonempty := by
    by_contra hemp
    rw [Set.not_nonempty_iff_eq_empty] at hemp
    rw [hemp, closure_empty] at hcu
    exact hcu
  obtain ⟨ε, hm, hε⟩ := hne
  have hcA : closure (φ ⁻¹' hull u ∩ Ioc 0 1) ⊆ φ ⁻¹' hull u :=
    closure_minimal Set.inter_subset_left
      (((Set.toFinite _).isCompact_convexHull (𝕜 := ℝ)).isClosed.preimage hφc)
  have hqu : q ∈ hull u := by
    have := hcA hcu
    simp only [Set.mem_preimage, hφ, zero_smul, add_zero, sub_zero, one_smul] at this
    exact this
  have hne' : u ≠ {a, b, c} := by
    rintro rfl
    simp only [Set.mem_preimage] at hm
    rw [hull_triple] at hm
    have h3 := ((mem_hull3_iff hD).1 hm).2.2
    rw [hφs] at h3
    have : sda a b c * (ε * sda a b x) = ε * (sda a b c * sda a b x) := by ring
    rw [this] at h3
    nlinarith [hε.1]
  obtain ⟨d, rfl, hd⟩ := eq_insert_of_mem_edge hT hu hq1 hq2 hq3 hqS (Finset.mem_coe.1 hu') hne' hqu
  exact ⟨d, Finset.mem_coe.1 hu', hd⟩

/-! ## Visible edges and cones -/

/-- The edge `ab` of the triangle `{a, b, c} ∈ T` is **visible from `p`**: `p` is strictly across
the line `ab` from `c`, and every point of `S` is on `c`'s closed side. -/
def Visible (S : Finset (ℝ × ℝ)) (T : Finset (Finset (ℝ × ℝ))) (p a b c : ℝ × ℝ) : Prop :=
  {a, b, c} ∈ T ∧ sda a b c * sda a b p < 0 ∧ ∀ s ∈ S, 0 ≤ sda a b c * sda a b s

/-- A closed half-plane `{z | 0 ≤ k · sda a b z}` is convex. -/
theorem convex_halfplane (a b : ℝ × ℝ) (k : ℝ) : Convex ℝ {z : ℝ × ℝ | 0 ≤ k * sda a b z} := by
  intro y hy w hw s t hs ht hst
  simp only [Set.mem_ofPred_eq] at hy hw ⊢
  rw [sda_combo₃ _ _ _ _ hst]
  nlinarith [mul_nonneg hs hy, mul_nonneg ht hw]

/-- If all of `S` is on one closed side of the line `ab`, so is `conv S`. -/
theorem halfplane_of_forall {S : Finset (ℝ × ℝ)} {a b : ℝ × ℝ} {k : ℝ}
    (h : ∀ s ∈ S, 0 ≤ k * sda a b s) {z : ℝ × ℝ} (hz : z ∈ convexHull ℝ (S : Set (ℝ × ℝ))) :
    0 ≤ k * sda a b z :=
  convexHull_min (fun s hs => h s hs) (convex_halfplane a b k) hz

theorem Visible.halfplane {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))} {p a b c : ℝ × ℝ}
    (h : Visible S T p a b c) {z : ℝ × ℝ} (hz : z ∈ convexHull ℝ (S : Set (ℝ × ℝ))) :
    0 ≤ sda a b c * sda a b z :=
  halfplane_of_forall h.2.2 hz

theorem Visible.sda_ne_zero {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))} {p a b c : ℝ × ℝ}
    (h : Visible S T p a b c) : sda a b p ≠ 0 := by
  intro h0; have := h.2.1; rw [h0, mul_zero] at this; exact lt_irrefl _ this

open Classical in
/-- The **cones** `{a, b, p}` over the edges `ab` of `T` that are visible from `p`. -/
def cones (S : Finset (ℝ × ℝ)) (T : Finset (Finset (ℝ × ℝ))) (p : ℝ × ℝ) :
    Finset (Finset (ℝ × ℝ)) :=
  ((S ×ˢ S ×ˢ S).filter fun q => Visible S T p q.1 q.2.1 q.2.2).image fun q => {q.1, q.2.1, p}

theorem mem_cones {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))} (hT : IsTriangulation S T)
    {p : ℝ × ℝ} {r : Finset (ℝ × ℝ)} :
    r ∈ cones S T p ↔ ∃ a b c, Visible S T p a b c ∧ r = {a, b, p} := by
  classical
  unfold cones
  simp only [Finset.mem_image, Finset.mem_filter, Finset.mem_product]
  constructor
  · rintro ⟨⟨a, b, c⟩, ⟨-, hv⟩, rfl⟩
    exact ⟨a, b, c, hv, rfl⟩
  · rintro ⟨a, b, c, hv, rfl⟩
    exact ⟨(a, b, c), ⟨⟨hT.mem_S hv.1 (by simp), hT.mem_S hv.1 (by simp),
      hT.mem_S hv.1 (by simp)⟩, hv⟩, rfl⟩

/-! ## Disjoint interiors -/

/-- **A cone and an old triangle have disjoint interiors**: the cone is strictly beyond the
supporting line, the triangle in `conv S` before it. -/
theorem disjoint_cone_old {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {p a b c : ℝ × ℝ} (hv : Visible S T p a b c) {u : Finset (ℝ × ℝ)}
    (hu : u ∈ T) : Disjoint (interior (hull {a, b, p})) (interior (hull u)) := by
  rw [Set.disjoint_left]
  intro z hz1 hz2
  rw [hull_triple] at hz1
  have h1 := (pos_of_mem_interior_hull3 hv.sda_ne_zero hz1).2.2
  have h2 := hv.halfplane (hT.hull_subset hu (interior_subset hz2))
  have h0 := eq_zero_of_signs h2 h1.le hv.2.1
  rw [h0, mul_zero] at h1
  exact lt_irrefl _ h1

/-- An interior point `z` of a nondegenerate cone `{a, b, p}` lies on the open segment from `p` to
a point of the open segment `ab`. -/
theorem interior_cone {a b p z : ℝ × ℝ} (hP : sda a b p ≠ 0)
    (hz : z ∈ interior (convexHull ℝ {a, b, p})) :
    ∃ g l : ℝ, 0 < g ∧ g < 1 ∧ 0 < l ∧ l < 1 ∧ z = g • p + (1 - g) • ((1 - l) • a + l • b) := by
  obtain ⟨h1, h2, h3⟩ := pos_of_mem_interior_hull3 hP hz
  have hα := div_pos_of_mul_pos' hP h1
  have hβ := div_pos_of_mul_pos' hP h2
  have hγ := div_pos_of_mul_pos' hP h3
  have hsum : sda z b p / sda a b p + sda a z p / sda a b p + sda a b z / sda a b p = 1 := by
    rw [← add_div, ← add_div, ← sda_four, div_self hP]
  have hz' := sda_cramer hP z
  set α := sda z b p / sda a b p
  set β := sda a z p / sda a b p
  set γ := sda a b z / sda a b p
  have hab : α + β ≠ 0 := by linarith
  refine ⟨γ, β / (α + β), hγ, by linarith, div_pos hβ (by linarith),
    (div_lt_one (by linarith)).2 (by linarith), ?_⟩
  have h1γ : 1 - γ = α + β := by linarith
  rw [h1γ]
  conv_lhs => rw [hz']
  ext
  · simp only [Prod.fst_add, Prod.smul_fst, smul_eq_mul]
    field_simp
    ring
  · simp only [Prod.snd_add, Prod.smul_snd, smul_eq_mul]
    field_simp
    ring

/-- **Two different cones have disjoint interiors.** -/
theorem disjoint_cone_cone {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {p a b c a' b' c' : ℝ × ℝ} (hv : Visible S T p a b c)
    (hv' : Visible S T p a' b' c') (hne : ({a, b, p} : Finset (ℝ × ℝ)) ≠ {a', b', p}) :
    Disjoint (interior (hull {a, b, p})) (interior (hull {a', b', p})) := by
  rw [Set.disjoint_left]
  intro z hz hz'
  rw [hull_triple] at hz hz'
  have hD := hT.sda_ne_zero hv.1
  obtain ⟨g, l, hg0, hg1, hl0, hl1, hzq⟩ := interior_cone hv.sda_ne_zero hz
  obtain ⟨g', l', hg0', hg1', hl0', hl1', hzq'⟩ := interior_cone hv'.sda_ne_zero hz'
  set q := (1 - l) • a + l • b with hq
  set q' := (1 - l') • a' + l' • b' with hq'
  have fq : sda a b q = 0 := by
    rw [hq, sda_combo₃ _ _ _ _ (by ring), sda_self₁₃, sda_self₂₃]; ring
  have fq' : sda a' b' q' = 0 := by
    rw [hq', sda_combo₃ _ _ _ _ (by ring), sda_self₁₃, sda_self₂₃]; ring
  have hqu : q ∈ convexHull ℝ {a, b, c} :=
    convex_convexHull ℝ _ (subset_convexHull ℝ _ (by simp)) (subset_convexHull ℝ _ (by simp))
      (by linarith) hl0.le (by ring)
  have hqu' : q' ∈ convexHull ℝ {a', b', c'} :=
    convex_convexHull ℝ _ (subset_convexHull ℝ _ (by simp)) (subset_convexHull ℝ _ (by simp))
      (by linarith) hl0'.le (by ring)
  have hqK : q ∈ convexHull ℝ (S : Set (ℝ × ℝ)) := hT.hull_subset hv.1 (by rwa [hull_triple])
  have hqK' : q' ∈ convexHull ℝ (S : Set (ℝ × ℝ)) := hT.hull_subset hv'.1 (by rwa [hull_triple])
  -- the two supporting lines along the ray from `p` through `z`
  have e1 : sda a' b' z = g * sda a' b' p + (1 - g) * sda a' b' q := by
    rw [hzq, sda_combo₃ _ _ _ _ (by ring)]
  have e1' : sda a' b' z = g' * sda a' b' p := by
    rw [hzq', sda_combo₃ _ _ _ _ (by ring), fq']; ring
  have e2 : sda a b z = g' * sda a b p + (1 - g') * sda a b q' := by
    rw [hzq', sda_combo₃ _ _ _ _ (by ring)]
  have e2' : sda a b z = g * sda a b p := by
    rw [hzq, sda_combo₃ _ _ _ _ (by ring), fq]; ring
  have s1 := hv'.halfplane hqK
  have s2 := hv.halfplane hqK'
  have k1 : (g' - g) * (sda a' b' c' * sda a' b' p) = (1 - g) * (sda a' b' c' * sda a' b' q) := by
    linear_combination sda a' b' c' * (e1'.symm.trans e1)
  have k2 : (g - g') * (sda a b c * sda a b p) = (1 - g') * (sda a b c * sda a b q') := by
    linear_combination sda a b c * (e2'.symm.trans e2)
  have i1 : g' ≤ g := by nlinarith [mul_nonneg (sub_nonneg.2 hg1.le) s1, hv'.2.1]
  have i2 : g ≤ g' := by nlinarith [mul_nonneg (sub_nonneg.2 hg1'.le) s2, hv.2.1]
  have hgg : g = g' := le_antisymm i2 i1
  have hqq : q' = q := by
    have h := hzq'.symm.trans hzq
    rw [hgg] at h
    have h' := add_left_cancel h
    exact smul_right_injective (ℝ × ℝ) (sub_ne_zero.2 (ne_of_gt hg1')) h'
  rw [hqq] at fq' hqu'
  -- `q` is a relative interior point of the edge `ab` of `{a, b, c}`
  have eq1 : sda q b c = (1 - l) * sda a b c := by
    rw [hq, sda_combo₁ _ _ _ _ (by ring), sda_self₁₂]; ring
  have eq2 : sda a q c = l * sda a b c := by
    rw [hq, sda_combo₂ _ _ _ _ (by ring), sda_self₁₂]; ring
  have hDD := mul_self_pos.2 hD
  have hq1 : 0 < sda a b c * sda q b c := by
    rw [eq1, mul_left_comm]; exact mul_pos (by linarith) hDD
  have hq2 : 0 < sda a b c * sda a q c := by
    rw [eq2, mul_left_comm]; exact mul_pos hl0 hDD
  have hqS := not_mem_of_mem_edge hT hv.1 hq1 hq2 fq
  by_cases ht : ({a', b', c'} : Finset (ℝ × ℝ)) = {a, b, c}
  · have ha' : a' ∈ ({a, b, c} : Finset (ℝ × ℝ)) := ht ▸ by simp
    have hb' : b' ∈ ({a, b, c} : Finset (ℝ × ℝ)) := ht ▸ by simp
    have hab' : a' ≠ b' := sda_ne_zero_left (hT.sda_ne_zero hv'.1)
    have hqa : sda a q c ≠ 0 := fun h0 => by rw [h0, mul_zero] at hq2; exact lt_irrefl _ hq2
    have hqb : sda q b c ≠ 0 := fun h0 => by rw [h0, mul_zero] at hq1; exact lt_irrefl _ hq1
    simp only [Finset.mem_insert, Finset.mem_singleton] at ha' hb'
    rcases ha' with rfl | rfl | rfl <;> rcases hb' with rfl | rfl | rfl
    · exact hab' rfl
    · exact hne rfl
    · exact hqa (by unfold sda at fq' ⊢; linear_combination -fq')
    · exact hne (Finset.insert_comm _ _ _)
    · exact hab' rfl
    · exact hqb (by unfold sda at fq' ⊢; linear_combination fq')
    · exact hqa (by unfold sda at fq' ⊢; linear_combination fq')
    · exact hqb (by unfold sda at fq' ⊢; linear_combination -fq')
    · exact hab' rfl
  · obtain ⟨d, hd, hds⟩ := eq_insert_of_mem_edge hT hv.1 hq1 hq2 fq hqS hv'.1 ht
      (by rw [hull_triple]; exact hqu')
    have hdS : d ∈ S := hT.mem_S hv'.1 (by rw [hd]; simp)
    have := hv.2.2 d hdS
    linarith

/-! ## The vertex condition -/

/-- **A point of `S` in a cone is a corner of it** (it lies on the supporting line, hence on the
edge, hence is a corner of the old triangle there). -/
theorem vertex_cone {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))} (hT : IsTriangulation S T)
    {p a b c : ℝ × ℝ} (hv : Visible S T p a b c) {s : ℝ × ℝ} (hs : s ∈ S)
    (hsr : s ∈ convexHull ℝ {a, b, p}) : s ∈ ({a, b, p} : Finset (ℝ × ℝ)) := by
  have hP := hv.sda_ne_zero
  have hD := hT.sda_ne_zero hv.1
  obtain ⟨g1, g2, g3⟩ := (mem_hull3_iff hP).1 hsr
  have hZ := eq_zero_of_signs (hv.2.2 s hs) g3 hv.2.1
  have hsu : s ∈ hull {a, b, c} := by
    rw [hull_triple]
    refine (mem_hull3_iff hD).2 ⟨?_, ?_, by rw [hZ, mul_zero]⟩
    · exact nonneg_of_edge_transfer (edge_transfer₁ hZ p c) hP g1
    · exact nonneg_of_edge_transfer (edge_transfer₂ hZ p c) hP g2
  have := hT.vertex _ hv.1 s hs hsu
  simp only [Finset.mem_insert, Finset.mem_singleton] at this ⊢
  rcases this with h | h | h
  · exact Or.inl h
  · exact Or.inr (Or.inl h)
  · exact absurd (h ▸ hZ) hD

/-! ## The cover -/

/-- One step of the cover argument, with the walk `γ τ = (1 − τ) y + τ x` leaving `conv S` for the
last time at `q = γ τ'` through the edge `ab` of `{a, b, c} ∈ T`, beyond which `x` lies strictly.
Here `x = α p + (1 − α) y` with `y ∈ conv S`, and `q ∉ S`. Then `ab` is visible from `p` and `x`
lies in the cone `{a, b, p}`. -/
theorem cone_step {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))} (hT : IsTriangulation S T)
    {p x y q a b c : ℝ × ℝ} {α τ' : ℝ} {γ : ℝ → ℝ × ℝ} (hγ : ∀ τ, γ τ = (1 - τ) • y + τ • x)
    (hτ0 : 0 ≤ τ') (hτ1 : τ' < 1) (hq : γ τ' = q)
    (hout : ∀ τ ∈ Ioc τ' 1, γ τ ∉ convexHull ℝ (S : Set (ℝ × ℝ)))
    (hy : y ∈ convexHull ℝ (S : Set (ℝ × ℝ))) (hα0 : 0 ≤ α) (hα1 : α ≤ 1)
    (hx : x = α • p + (1 - α) • y) (hqS : q ∉ S) (hu : {a, b, c} ∈ T)
    (hqu : q ∈ convexHull ℝ {a, b, c}) (hq3 : sda a b q = 0)
    (hxs : sda a b c * sda a b x < 0) :
    ∃ a b c, Visible S T p a b c ∧ x ∈ convexHull ℝ {a, b, p} := by
  have hD := hT.sda_ne_zero hu
  obtain ⟨hy1, hy2, -⟩ := (mem_hull3_iff hD).1 hqu
  have hy1' : 0 < sda a b c * sda q b c := lt_of_le_of_ne hy1 fun h0 => hqS (by
    rw [eq_of_sda_eq_zero hD ((mul_eq_zero.1 h0.symm).resolve_left hD) hq3]
    exact hT.mem_S hu (by simp))
  have hy2' : 0 < sda a b c * sda a q c := lt_of_le_of_ne hy2 fun h0 => hqS (by
    rw [eq_of_sda_eq_zero' hD ((mul_eq_zero.1 h0.symm).resolve_left hD) hq3]
    exact hT.mem_S hu (by simp))
  -- the walk beyond `q`
  have hγτ : ∀ τ, (1 - τ') * sda a b (γ τ) = (τ - τ') * sda a b x := by
    intro τ
    have e1 : sda a b (γ τ) = (1 - τ) * sda a b y + τ * sda a b x := by
      rw [hγ, sda_combo₃ _ _ _ _ (by ring)]
    have e2 : sda a b q = (1 - τ') * sda a b y + τ' * sda a b x := by
      rw [← hq, hγ, sda_combo₃ _ _ _ _ (by ring)]
    rw [e1]
    rw [hq3] at e2
    linear_combination (-(1 - τ)) * e2
  have hγc : Continuous γ := by
    have : γ = fun τ => (1 - τ) • y + τ • x := funext hγ
    rw [this]; fun_prop
  -- no triangle across `ab`, so all of `S` is on `c`'s side
  have hall : ∀ s ∈ S, 0 ≤ sda a b c * sda a b s := by
    intro s hs
    by_contra hneg
    push Not at hneg
    obtain ⟨d, hd, hdD⟩ := exists_across hT hu hy1' hy2' hq3 (subset_convexHull ℝ _ hs) hneg
    have hDd := hT.sda_ne_zero hd
    have hO : IsOpen {z : ℝ × ℝ | 0 < sda a b d * sda z b d ∧ 0 < sda a b d * sda a z d} :=
      (isOpen_lt continuous_const (continuous_const.mul (continuous_sda₁ b d))).inter
        (isOpen_lt continuous_const (continuous_const.mul (continuous_sda₂ a d)))
    have hqO : γ τ' ∈ {z : ℝ × ℝ | 0 < sda a b d * sda z b d ∧ 0 < sda a b d * sda a z d} := by
      rw [hq]
      exact ⟨pos_of_edge_transfer (edge_transfer₁ hq3 c d) hDd hy1',
        pos_of_edge_transfer (edge_transfer₂ hq3 c d) hDd hy2'⟩
    have ev := ((hγc.tendsto τ').eventually (hO.mem_nhds hqO)).filter_mono (nhdsWithin_le_nhds (s := Ioi τ'))
    obtain ⟨τ, hτO, hτI⟩ := (ev.and (Ioo_mem_nhdsGT hτ1)).exists
    have hdx : 0 < sda a b d * sda a b x := by
      have h := mul_pos_of_neg_of_neg hdD hxs
      have : 0 < (sda a b c * sda a b c) * (sda a b d * sda a b x) := by linarith
      exact pos_of_mul_pos_right this (mul_self_nonneg _)
    have h3 : 0 < sda a b d * sda a b (γ τ) := by
      have h := hγτ τ
      have : (1 - τ') * (sda a b d * sda a b (γ τ)) = (τ - τ') * (sda a b d * sda a b x) := by
        linear_combination sda a b d * h
      have hpos : 0 < (1 - τ') * (sda a b d * sda a b (γ τ)) := by
        rw [this]; exact mul_pos (by linarith [hτI.1]) hdx
      exact pos_of_mul_pos_right hpos (by linarith)
    have hmem : γ τ ∈ convexHull ℝ {a, b, d} := (mem_hull3_iff hDd).2 ⟨hτO.1.le, hτO.2.le, h3.le⟩
    exact hout τ ⟨hτI.1, hτI.2.le⟩ (hT.hull_subset hd (by rw [hull_triple]; exact hmem))
  -- `p` is strictly beyond `ab`
  have hpside : sda a b c * sda a b p < 0 := by
    have hyD := halfplane_of_forall hall hy
    have e : sda a b x = α * sda a b p + (1 - α) * sda a b y := by
      rw [hx, sda_combo₃ _ _ _ _ (by ring)]
    rw [e] at hxs
    by_contra hge
    push Not at hge
    nlinarith [mul_nonneg hα0 hge, mul_nonneg (sub_nonneg.2 hα1) hyD]
  refine ⟨a, b, c, ⟨hu, hpside, hall⟩, ?_⟩
  -- `x` is on the segment from `p` to `q`
  have hP : sda a b p ≠ 0 := by intro h0; rw [h0, mul_zero] at hpside; exact lt_irrefl _ hpside
  have hqP : q ∈ convexHull ℝ {a, b, p} :=
    (mem_hull3_iff hP).2 ⟨nonneg_of_edge_transfer (edge_transfer₁ hq3 c p) hD hy1,
      nonneg_of_edge_transfer (edge_transfer₂ hq3 c p) hD hy2, by rw [hq3, mul_zero]⟩
  have hpP : p ∈ convexHull ℝ {a, b, p} := subset_convexHull ℝ _ (by simp)
  have key : (1 - α * τ') • x = (α * (1 - τ')) • p + (1 - α) • q := by
    rw [← hq, hγ]
    have h1 := congrArg Prod.fst hx
    have h2 := congrArg Prod.snd hx
    simp only [Prod.fst_add, Prod.smul_fst, smul_eq_mul, Prod.snd_add, Prod.smul_snd] at h1 h2
    ext
    · simp only [Prod.fst_add, Prod.smul_fst, smul_eq_mul]
      linear_combination (1 - τ') * h1
    · simp only [Prod.snd_add, Prod.smul_snd, smul_eq_mul]
      linear_combination (1 - τ') * h2
  have hne : 1 - α * τ' ≠ 0 := by nlinarith
  have hpos : 0 < 1 - α * τ' := by nlinarith
  have hxe : x = ((1 - α * τ')⁻¹ * (α * (1 - τ'))) • p + ((1 - α * τ')⁻¹ * (1 - α)) • q := by
    have : x = (1 - α * τ')⁻¹ • ((1 - α * τ') • x) := by
      rw [smul_smul, inv_mul_cancel₀ hne, one_smul]
    rw [this, key, smul_add, smul_smul, smul_smul]
  rw [hxe]
  refine convex_convexHull ℝ _ hpP hqP ?_ ?_ ?_
  · exact mul_nonneg (inv_nonneg.2 hpos.le) (mul_nonneg hα0 (by linarith))
  · exact mul_nonneg (inv_nonneg.2 hpos.le) (by linarith)
  · rw [← mul_add, inv_mul_eq_one₀ hne]; ring

/-- **A generic point of `conv(S ∪ {p}) \ conv S` is in a cone.** If `p ∉ conv S` and the line
through `x` and `p` contains no point of `S`, every `x ∈ conv(S ∪ {p})` outside `conv S` lies in a
cone `{a, b, p}` over an edge `ab` visible from `p`. -/
theorem exists_cone_of_generic {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (hSne : S.Nonempty) {p x : ℝ × ℝ}
    (hx : x ∈ convexHull ℝ (insert p (S : Set (ℝ × ℝ))))
    (hxS : x ∉ convexHull ℝ (S : Set (ℝ × ℝ))) (hgen : ∀ s ∈ S, sda x p s ≠ 0) :
    ∃ a b c, Visible S T p a b c ∧ x ∈ convexHull ℝ {a, b, p} := by
  rw [convexHull_insert (Finset.coe_nonempty.2 hSne)] at hx
  obtain ⟨p', hp', y, hy, hxy⟩ := mem_convexJoin.1 hx
  rw [Set.mem_singleton_iff] at hp'
  subst hp'
  obtain ⟨α, β, hα, hβ, hαβ, hxe⟩ := hxy
  have hβ' : β = 1 - α := by linarith
  subst hβ'
  have hxe' : x = α • p' + (1 - α) • y := hxe.symm
  have hα1 : α ≤ 1 := by linarith
  -- the walk from `y` to `x`
  obtain ⟨γ, hγ⟩ : ∃ γ : ℝ → ℝ × ℝ, ∀ τ, γ τ = (1 - τ) • y + τ • x := ⟨_, fun _ => rfl⟩
  have hγc : Continuous γ := by
    have : γ = fun τ => (1 - τ) • y + τ • x := funext hγ
    rw [this]; fun_prop
  have hKc : IsClosed (convexHull ℝ (S : Set (ℝ × ℝ))) :=
    (S.finite_toSet.isCompact_convexHull (𝕜 := ℝ)).isClosed
  let J : Set ℝ := Icc 0 1 ∩ γ ⁻¹' convexHull ℝ (S : Set (ℝ × ℝ))
  have hJc : IsClosed J := isClosed_Icc.inter (hKc.preimage hγc)
  have hJcpt : IsCompact J := isCompact_Icc.of_isClosed_subset hJc Set.inter_subset_left
  have h0J : (0 : ℝ) ∈ J := ⟨⟨le_rfl, zero_le_one⟩, by
    show γ 0 ∈ convexHull ℝ (S : Set (ℝ × ℝ)); rw [hγ]; simpa using hy⟩
  have hτJ : sSup J ∈ J := hJcpt.sSup_mem ⟨0, h0J⟩
  have hmax : ∀ τ ∈ J, τ ≤ sSup J := fun τ hτ => le_csSup hJcpt.bddAbove hτ
  set τ' := sSup J
  obtain ⟨⟨hτ0, hτ1⟩, hτK⟩ := hτJ
  simp only [Set.mem_preimage] at hτK
  have hγ1 : γ 1 = x := by rw [hγ]; simp
  have hτ1' : τ' < 1 := lt_of_le_of_ne hτ1 fun h => hxS (by rw [← hγ1, ← h]; exact hτK)
  have hout : ∀ τ ∈ Ioc τ' 1, γ τ ∉ convexHull ℝ (S : Set (ℝ × ℝ)) := fun τ hτ hm =>
    not_le.2 hτ.1 (hmax τ ⟨⟨by linarith [hτ.1], hτ.2⟩, hm⟩)
  rw [← hT.cover] at hτK
  obtain ⟨t, ht, hqt⟩ := Set.mem_iUnion₂.1 hτK
  obtain ⟨a, b, c, hD, rfl⟩ := hT.tri t ht
  rw [hull_triple] at hqt
  -- the exit point is not a point of `S`
  have hqS : γ τ' ∉ S := by
    intro hm
    refine hgen _ hm ?_
    rw [hγ, sda_combo₃ _ _ _ _ (by ring), sda_self₁₃, hxe', sda_combo₁ _ _ _ _ (by ring),
      sda_self₁₂, sda_self₁₃]
    ring
  have hγ' : ∀ τ, γ τ = y + τ • (x - y) := fun τ => by
    rw [hγ, smul_sub, sub_smul, one_smul]; abel
  have hout' : ∀ τ ∈ Ioc τ' 1, γ τ ∉ convexHull ℝ {a, b, c} := fun τ hτ hm =>
    hout τ hτ (hT.hull_subset ht (by rw [hull_triple]; exact hm))
  rcases exit_edge hD hγ' hτ1' hqt hout' with ⟨h3, hs⟩ | ⟨h1, hs⟩ | ⟨h2, hs⟩
  · exact cone_step hT hγ hτ0 hτ1' rfl hout hy hα hα1 hxe' hqS ht hqt h3 hs
  · refine cone_step (a := b) (b := c) (c := a) hT hγ hτ0 hτ1' rfl hout hy hα hα1 hxe' hqS
      (by convert ht using 1; ext; simp; tauto) (by rw [← hull_rot a b c]; exact hqt)
      (by rw [← h1]; unfold sda; ring) ?_
    have : sda b c a * sda b c x = sda a b c * sda x b c := by unfold sda; ring
    rwa [this]
  · refine cone_step (a := c) (b := a) (c := b) hT hγ hτ0 hτ1' rfl hout hy hα hα1 hxe' hqS
      (by convert ht using 1; ext; simp; tauto)
      (by rw [← hull_rot b c a, ← hull_rot a b c]; exact hqt) (by rw [← h2]; unfold sda; ring) ?_
    have : sda c a b * sda c a x = sda a b c * sda a x c := by unfold sda; ring
    rwa [this]

/-- **The cover.** The old triangles and the cones cover `conv(S ∪ {p})`. -/
theorem cover_insert {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (hSne : S.Nonempty) {p : ℝ × ℝ}
    (hp : p ∉ convexHull ℝ (S : Set (ℝ × ℝ))) :
    ⋃ t ∈ T ∪ cones S T p, hull t = convexHull ℝ ((insert p S : Finset (ℝ × ℝ)) : Set (ℝ × ℝ)) := by
  classical
  have hpS : p ∉ S := fun h => hp (subset_convexHull ℝ _ h)
  have hsub : convexHull ℝ (S : Set (ℝ × ℝ)) ⊆
      convexHull ℝ ((insert p S : Finset (ℝ × ℝ)) : Set (ℝ × ℝ)) :=
    convexHull_mono (Finset.coe_subset.2 (Finset.subset_insert _ _))
  apply Set.Subset.antisymm
  · refine Set.iUnion₂_subset fun t ht => ?_
    rcases Finset.mem_union.1 ht with ht | ht
    · exact (hT.hull_subset ht).trans hsub
    · obtain ⟨a, b, c, hv, rfl⟩ := (mem_cones hT).1 ht
      refine convexHull_mono ?_
      intro z hz
      simp only [Finset.coe_insert, Finset.coe_singleton, Set.mem_insert_iff,
        Set.mem_singleton_iff] at hz ⊢
      rcases hz with rfl | rfl | rfl
      · exact Or.inr (hT.mem_S hv.1 (by simp))
      · exact Or.inr (hT.mem_S hv.1 (by simp))
      · exact Or.inl rfl
  · set K := ⋃ t ∈ T ∪ cones S T p, hull t
    set C := convexHull ℝ ((insert p S : Finset (ℝ × ℝ)) : Set (ℝ × ℝ))
    have hKc : IsClosed K := isClosed_biUnion_finset fun t _ =>
      ((Set.toFinite _).isCompact_convexHull (𝕜 := ℝ)).isClosed
    -- generic points are covered
    have hgen : ∀ x ∈ C, (∀ s ∈ S, sda x p s ≠ 0) → x ∈ K := by
      intro x hx hg
      by_cases hxS : x ∈ convexHull ℝ (S : Set (ℝ × ℝ))
      · rw [← hT.cover] at hxS
        obtain ⟨t, ht, hm⟩ := Set.mem_iUnion₂.1 hxS
        exact Set.mem_biUnion (Finset.mem_union_left _ ht) hm
      · have hx' : x ∈ convexHull ℝ (insert p (S : Set (ℝ × ℝ))) := by
          rw [← Finset.coe_insert]; exact hx
        obtain ⟨a, b, c, hv, hm⟩ := exists_cone_of_generic hT hSne hx' hxS hg
        exact Set.mem_biUnion (Finset.mem_union_right _ ((mem_cones hT).2 ⟨a, b, c, hv, rfl⟩))
          (by rw [hull_triple]; exact hm)
    -- generic points are dense
    have hG : Dense {z : ℝ × ℝ | ∀ s ∈ S, sda z p s ≠ 0} := by
      have h := interior_biUnion_eq_empty S (fun s => {z : ℝ × ℝ | sda z p s = 0})
        (fun s _ => isClosed_eq (continuous_sda₁ p s) continuous_const)
        (fun s hs => interior_line_eq_empty fun h => hpS (h ▸ hs))
      have hd := interior_eq_empty_iff_dense_compl.1 h
      convert hd using 1
      ext z
      simp
    -- `C` has interior points
    have hCi : (interior C).Nonempty := by
      obtain ⟨s, hs⟩ := hSne
      have : s ∈ convexHull ℝ (S : Set (ℝ × ℝ)) := subset_convexHull ℝ _ hs
      rw [← hT.cover] at this
      obtain ⟨t, ht, -⟩ := Set.mem_iUnion₂.1 this
      obtain ⟨a, b, c, hD, rfl⟩ := hT.tri t ht
      obtain ⟨z, hz⟩ := interior_hull3_nonempty hD
      refine ⟨z, interior_mono ((hT.hull_subset ht).trans hsub) ?_⟩
      rw [hull_triple]; exact hz
    intro x hx
    have h1 : x ∈ closure (interior C) := by
      rw [(convex_convexHull ℝ _).closure_interior_eq_closure_of_nonempty_interior hCi]
      exact subset_closure hx
    have h2 : closure (interior C) ⊆ closure (interior C ∩ {z | ∀ s ∈ S, sda z p s ≠ 0}) :=
      closure_minimal (hG.open_subset_closure_inter isOpen_interior) isClosed_closure
    have h3 : closure (interior C ∩ {z | ∀ s ∈ S, sda z p s ≠ 0}) ⊆ K :=
      closure_minimal (fun z hz => hgen z (interior_subset hz.1) hz.2) hKc
    exact h3 (h2 h1)

/-! ## The insertion theorem -/

/-- **Inserting a point outside the convex hull.** If `T` is a triangulation of the nonempty finite
point set `S` and `p ∉ conv S`, then `T` together with the cones from `p` over the edges of `T`
visible from `p` is a triangulation of `S ∪ {p}`.

[formula: planar-triangulation-insertion] -/
theorem isTriangulation_insert {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (hSne : S.Nonempty) {p : ℝ × ℝ}
    (hp : p ∉ convexHull ℝ (S : Set (ℝ × ℝ))) :
    IsTriangulation (insert p S) (T ∪ cones S T p) := by
  classical
  have hpS : p ∉ S := fun h => hp (subset_convexHull ℝ _ h)
  refine ⟨?_, ?_, cover_insert hT hSne hp, ?_, ?_⟩
  · intro t ht
    rcases Finset.mem_union.1 ht with ht | ht
    · exact hT.tri t ht
    · obtain ⟨a, b, c, hv, rfl⟩ := (mem_cones hT).1 ht
      exact ⟨a, b, p, hv.sda_ne_zero, rfl⟩
  · intro t ht
    rcases Finset.mem_union.1 ht with ht | ht
    · exact (hT.subset t ht).trans (Finset.subset_insert _ _)
    · obtain ⟨a, b, c, hv, rfl⟩ := (mem_cones hT).1 ht
      intro z hz
      simp only [Finset.mem_insert, Finset.mem_singleton] at hz ⊢
      rcases hz with rfl | rfl | rfl
      · exact Or.inr (hT.mem_S hv.1 (by simp))
      · exact Or.inr (hT.mem_S hv.1 (by simp))
      · exact Or.inl rfl
  · intro t ht u hu htu
    rw [Finset.coe_union] at ht hu
    rcases ht with ht | ht <;> rcases hu with hu | hu
    · exact hT.disjoint ht hu htu
    · obtain ⟨a, b, c, hv, rfl⟩ := (mem_cones hT).1 hu
      exact (disjoint_cone_old hT hv ht).symm
    · obtain ⟨a, b, c, hv, rfl⟩ := (mem_cones hT).1 ht
      exact disjoint_cone_old hT hv hu
    · obtain ⟨a, b, c, hv, rfl⟩ := (mem_cones hT).1 ht
      obtain ⟨a', b', c', hv', rfl⟩ := (mem_cones hT).1 hu
      exact disjoint_cone_cone hT hv hv' htu
  · intro t ht s hs hst
    rcases Finset.mem_union.1 ht with htT | htC
    · rcases Finset.mem_insert.1 hs with rfl | hs
      · exact absurd (hT.hull_subset htT hst) hp
      · exact hT.vertex t htT s hs hst
    · obtain ⟨a, b, c, hv, rfl⟩ := (mem_cones hT).1 htC
      rcases Finset.mem_insert.1 hs with rfl | hs
      · simp
      · rw [hull_triple] at hst
        exact vertex_cone hT hv hs hst

end Planar

end Delaunay

end DDG
