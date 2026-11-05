import DdgProofs.Delaunay.PlanarBasic

/-!
# Planar Delaunay triangulations, part 2: triangulations of a point set and the Lawson flip

## The model

A **triangle** is a three-element set `t = {a, b, c}` of points of the plane `ℝ × ℝ` with
`sda a b c ≠ 0` (`IsTri`; the order of the corners is irrelevant, `IsTri.sda_ne_zero`). A
**triangulation of a finite point set `S`** (`IsTriangulation S T`) is a finite set `T` of
triangles with corners in `S` such that

* `cover`: the triangles cover the convex hull of `S`, `⋃_{t ∈ T} conv t = conv S`;
* `disjoint`: two different triangles have disjoint interiors;
* `vertex`: a point of `S` that lies in a triangle is one of its corners.

The last condition is the usual convention that the vertex set of the triangulation is all of `S`
and that no point of `S` lies on an edge or inside a triangle except as a corner: every point of
`S` lies in `conv S`, hence in some triangle, hence is a corner (de Berg et al., Ch. 9: a
triangulation of `P` is a maximal planar subdivision with vertex set `P`). Without it the Delaunay
property would be false for triangulations that ignore some of the points, and the adjacency of
triangles across an edge (`DdgProofs/Delaunay/PlanarLemma.lean`) would fail for points on edges.
It is preserved by the flip (`lawson_flip`). The interiors are those of `ℝ × ℝ` with its product
topology (the same as the Euclidean one); the boundaries of triangles are null sets for the
Lebesgue measure `volume`, which is how the flip keeps interiors disjoint
(`eq_empty_of_subset_union_hull`). This is the model of `DdgProofs/TriangleIntegral/Tiling.lean`
(`IsTiling`: union and disjoint interiors) with corners in `S` and the vertex condition.

Two triangles `{i, j, k}` and `{i, j, l}` of `T` with `k ≠ l` share the **interior edge** `ij`;
their third corners lie strictly on opposite sides of `ij` (`IsTriangulation.opposite`).

## Main statements

* `flip_hull_union` — for a strictly convex quadrilateral (the segments `ij` and `kl` cross,
  `sda i j k · sda i j l < 0` and `sda k l i · sda k l j < 0`) the triangles on the two diagonals
  cover the same set: `conv{i,k,l} ∪ conv{j,k,l} = conv{i,j,k} ∪ conv{i,j,l}`.
* **`lawson_flip`** — **the Lawson flip.** Replacing `{i,j,k}`, `{i,j,l}` by `{i,k,l}`, `{j,k,l}`
  in a triangulation of `S` (`flip T i j k l`) gives again a triangulation of `S` when `ijkl` is
  strictly convex: same union, disjoint interiors (also with all the other triangles, which are
  unchanged), the vertex condition.
* `isTriangulation_of_delaunayFlip` — a **Delaunay flip** (`DelaunayFlip T T'`: flip an interior
  edge `ij` with `l` strictly inside the circumcircle of `i j k`) always applies: such a diamond is
  strictly convex (`crossing_of_inCircumcircle`, the planar form of layer 1's
  `isConvexDiamond_of_not_isLocallyDelaunay`), so the result is a triangulation of `S`.

`IsLocallyDelaunay T` (no interior edge has its fourth corner strictly inside a circumcircle) and
`IsDelaunay S T` (no point of `S` strictly inside the circumcircle of any triangle) are defined
here; termination of Delaunay flipping is in `DdgProofs/Delaunay/PlanarTermination.lean`, and
Delaunay's lemma (locally Delaunay ⟹ Delaunay) in `DdgProofs/Delaunay/PlanarLemma.lean`.

## References

* C. L. Lawson, *Software for C¹ surface interpolation*, in: Mathematical Software III, Academic
  Press, 1977, 161–194 (the flip of the diagonal of a strictly convex quadrilateral).
* M. de Berg, O. Cheong, M. van Kreveld, M. Overmars, *Computational Geometry: Algorithms and
  Applications*, 3rd ed., Springer, 2008, Ch. 9 (triangulations of point sets, edge flips,
  legal edges).
* H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge University Press, 2001,
  Ch. 1.
* J.-F. Dufourd, Y. Bertot, *Formal study of plane Delaunay triangulation*, ITP 2010, LNCS 6172,
  211–226 (a Coq formalization of the planar flip algorithm on hypermaps; related work, a
  combinatorial model rather than the point-set model here).

The bibliographic data were not checked against the works here.
-/

noncomputable section

namespace DdgProofs

namespace Delaunay

namespace Planar

open MeasureTheory Set

/-! ## Triangles of a point set -/

/-- The (closed) triangle spanned by a finite set of points: its convex hull. -/
abbrev hull (t : Finset (ℝ × ℝ)) : Set (ℝ × ℝ) := convexHull ℝ (t : Set (ℝ × ℝ))

theorem hull_triple (a b c : ℝ × ℝ) : hull {a, b, c} = convexHull ℝ {a, b, c} := by
  simp [hull]

/-- `t` is a nondegenerate triangle: `t = {a, b, c}` with `sda a b c ≠ 0`. -/
def IsTri (t : Finset (ℝ × ℝ)) : Prop := ∃ a b c : ℝ × ℝ, sda a b c ≠ 0 ∧ t = {a, b, c}

/-- The orientation of a triangle is nonzero in **every** order of its corners. -/
theorem IsTri.sda_ne_zero {t : Finset (ℝ × ℝ)} (ht : IsTri t) {i j k : ℝ × ℝ}
    (hijk : t = {i, j, k}) : sda i j k ≠ 0 := by
  obtain ⟨a, b, c, h, rfl⟩ := ht
  have hab := sda_ne_zero_left h
  have hbc := sda_ne_zero_right h
  have hac := sda_ne_zero_outer h
  have mem : ∀ x, x ∈ ({a, b, c} : Finset (ℝ × ℝ)) ↔ x ∈ ({i, j, k} : Finset (ℝ × ℝ)) := by
    intro x; rw [hijk]
  have hi := (mem i).2 (by simp)
  have hj := (mem j).2 (by simp)
  have hk := (mem k).2 (by simp)
  have ha := (mem a).1 (by simp)
  have hb := (mem b).1 (by simp)
  have hc := (mem c).1 (by simp)
  simp only [Finset.mem_insert, Finset.mem_singleton] at hi hj hk ha hb hc
  rcases hi with rfl | rfl | rfl <;> rcases hj with rfl | rfl | rfl <;>
    rcases hk with rfl | rfl | rfl <;>
    first
    | (intro h0; apply h; unfold sda at h0 ⊢; linarith)
    | (exfalso; aesop)

/-- The three corners of a triangle are distinct. -/
theorem IsTri.ne {t : Finset (ℝ × ℝ)} (ht : IsTri t) {i j k : ℝ × ℝ} (hijk : t = {i, j, k}) :
    i ≠ j ∧ j ≠ k ∧ i ≠ k := by
  have h := ht.sda_ne_zero hijk
  exact ⟨sda_ne_zero_left h, sda_ne_zero_right h, sda_ne_zero_outer h⟩

/-- **A triangulation of the finite point set `S`** (see the module docstring). -/
structure IsTriangulation (S : Finset (ℝ × ℝ)) (T : Finset (Finset (ℝ × ℝ))) : Prop where
  /-- every member is a nondegenerate triangle -/
  tri : ∀ t ∈ T, IsTri t
  /-- with corners in `S` -/
  subset : ∀ t ∈ T, t ⊆ S
  /-- the triangles cover the convex hull of `S` -/
  cover : ⋃ t ∈ T, hull t = convexHull ℝ (S : Set (ℝ × ℝ))
  /-- different triangles have disjoint interiors -/
  disjoint : (T : Set (Finset (ℝ × ℝ))).Pairwise fun t u =>
    Disjoint (interior (hull t)) (interior (hull u))
  /-- a point of `S` in a triangle is a corner of it -/
  vertex : ∀ t ∈ T, ∀ s ∈ S, s ∈ hull t → s ∈ t

namespace IsTriangulation

variable {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}

theorem sda_ne_zero (hT : IsTriangulation S T) {i j k : ℝ × ℝ} (h : {i, j, k} ∈ T) :
    sda i j k ≠ 0 :=
  (hT.tri _ h).sda_ne_zero rfl

theorem mem_S (hT : IsTriangulation S T) {t : Finset (ℝ × ℝ)} (ht : t ∈ T) {x : ℝ × ℝ}
    (hx : x ∈ t) : x ∈ S :=
  hT.subset t ht hx

theorem hull_subset (hT : IsTriangulation S T) {t : Finset (ℝ × ℝ)} (ht : t ∈ T) :
    hull t ⊆ convexHull ℝ (S : Set (ℝ × ℝ)) := by
  rw [← hT.cover]
  exact Set.subset_biUnion_of_mem (u := fun t => hull t) ht

/-- Two triangles that share the corners `i, j` and differ in the third corner are different. -/
theorem triple_ne {i j k l : ℝ × ℝ} (hkl : k ≠ l) (hl : sda i j l ≠ 0) :
    ({i, j, k} : Finset (ℝ × ℝ)) ≠ {i, j, l} := by
  intro h
  have : l ∈ ({i, j, k} : Finset (ℝ × ℝ)) := by rw [h]; simp
  simp only [Finset.mem_insert, Finset.mem_singleton] at this
  rcases this with rfl | rfl | rfl
  · exact hl (sda_self₁₃ _ _)
  · exact hl (sda_self₂₃ _ _)
  · exact hkl rfl

/-- **The two triangles on an interior edge lie on opposite sides of it.** -/
theorem opposite (hT : IsTriangulation S T) {i j k l : ℝ × ℝ} (h1 : {i, j, k} ∈ T)
    (h2 : {i, j, l} ∈ T) (hkl : k ≠ l) : sda i j k * sda i j l < 0 := by
  have hD1 := hT.sda_ne_zero h1
  have hD2 := hT.sda_ne_zero h2
  by_contra hcon
  push Not at hcon
  have hpos : 0 < sda i j k * sda i j l := lt_of_le_of_ne hcon (mul_ne_zero hD1 hD2).symm
  have hdisj : Disjoint (interior (hull {i, j, k})) (interior (hull {i, j, l})) :=
    hT.disjoint h1 h2 (triple_ne hkl hD2)
  rw [hull_triple, hull_triple] at hdisj
  -- a point near the midpoint of `ij`, on the common side
  set m : ℝ × ℝ := (1 / 2 : ℝ) • i + (1 / 2 : ℝ) • j
  have em1 : sda m j k = sda i j k / 2 := by
    simp only [m]; unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst,
      Prod.smul_snd, smul_eq_mul]; ring
  have em2 : sda i m k = sda i j k / 2 := by
    simp only [m]; unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst,
      Prod.smul_snd, smul_eq_mul]; ring
  have em3 : sda m j l = sda i j l / 2 := by
    simp only [m]; unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst,
      Prod.smul_snd, smul_eq_mul]; ring
  have em4 : sda i m l = sda i j l / 2 := by
    simp only [m]; unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst,
      Prod.smul_snd, smul_eq_mul]; ring
  have hU : IsOpen {y : ℝ × ℝ | 0 < sda i j k * sda y j k ∧ 0 < sda i j k * sda i y k ∧
      0 < sda i j l * sda y j l ∧ 0 < sda i j l * sda i y l} :=
    (isOpen_lt continuous_const (continuous_const.mul (continuous_sda₁ j k))).inter
      ((isOpen_lt continuous_const (continuous_const.mul (continuous_sda₂ i k))).inter
      ((isOpen_lt continuous_const (continuous_const.mul (continuous_sda₁ j l))).inter
      (isOpen_lt continuous_const (continuous_const.mul (continuous_sda₂ i l)))))
  have hm : m ∈ {y : ℝ × ℝ | 0 < sda i j k * sda y j k ∧ 0 < sda i j k * sda i y k ∧
      0 < sda i j l * sda y j l ∧ 0 < sda i j l * sda i y l} := by
    have p1 := mul_self_pos.2 hD1
    have p2 := mul_self_pos.2 hD2
    simp only [Set.mem_ofPred_eq, em1, em2, em3, em4]
    refine ⟨?_, ?_, ?_, ?_⟩ <;> nlinarith
  obtain ⟨s, hs, hspos⟩ := ((exists_pos_near (hU.mem_nhds hm) (k - m)).and
    self_mem_nhdsWithin).exists
  have hspos' : (0 : ℝ) < s := hspos
  have e3 : sda i j (m + s • (k - m)) = s * sda i j k := by
    have h0 : sda i j m = 0 := by
      simp only [m]; unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst,
        Prod.smul_snd, smul_eq_mul]; ring
    have : m + s • (k - m) = (1 - s) • m + s • k := by
      rw [smul_sub, sub_smul, one_smul]; abel
    rw [this, sda_combo₃ _ _ _ _ (by ring), h0]; ring
  have i1 : m + s • (k - m) ∈ interior (convexHull ℝ {i, j, k}) :=
    mem_interior_hull3 hD1 hs.1 hs.2.1 (by rw [e3]; nlinarith [mul_self_pos.2 hD1])
  have i2 : m + s • (k - m) ∈ interior (convexHull ℝ {i, j, l}) :=
    mem_interior_hull3 hD2 hs.2.2.1 hs.2.2.2 (by rw [e3]; nlinarith)
  exact Set.disjoint_left.1 hdisj i1 i2

end IsTriangulation

/-! ## The two diagonals of a convex quadrilateral -/

/-- One half of `flip_hull_union`, with `sda i j k > 0`. -/
theorem hull_subset_flip_pos {i j k l : ℝ × ℝ} (hk : sda i j k * sda i j l < 0)
    (hc : sda k l i * sda k l j < 0) (hD : 0 < sda i j k) :
    convexHull ℝ {i, j, k} ⊆ convexHull ℝ {i, k, l} ∪ convexHull ℝ {j, k, l} := by
  obtain ⟨hl, hki, hkj⟩ := quad_signs hk hc hD
  intro x hx
  obtain ⟨x1, x2, x3⟩ := (mem_hull3_iff hD.ne').1 hx
  have y1 : 0 ≤ sda x j k := (mul_nonneg_iff_of_pos_left hD).1 x1
  have y2 : 0 ≤ sda i x k := (mul_nonneg_iff_of_pos_left hD).1 x2
  have y3 : 0 ≤ sda i j x := (mul_nonneg_iff_of_pos_left hD).1 x3
  have r1 : sda i k l = sda k l i := by unfold sda; ring
  have r2 : sda j k l = sda k l j := by unfold sda; ring
  have r3 : sda x k l = sda k l x := by unfold sda; ring
  rcases le_total (sda k l x) 0 with hx0 | hx0
  · left
    have hD1 : sda i k l ≠ 0 := by rw [r1]; exact hki.ne
    refine (mem_hull3_iff hD1).2 ⟨?_, ?_, ?_⟩
    · rw [r1, r3]; exact mul_nonneg_of_nonpos_of_nonpos hki.le hx0
    · have e : sda i j k * sda i x l = sda i x k * sda i j l + sda i j x * sda i k l := by
        unfold sda; ring
      have hil : sda i k l < 0 := by rw [r1]; exact hki
      have : sda i x l ≤ 0 := by
        by_contra hpos; push Not at hpos
        nlinarith [mul_nonpos_of_nonneg_of_nonpos y2 hl.le, mul_nonpos_of_nonneg_of_nonpos y3 hil.le,
          mul_pos hD hpos]
      exact mul_nonneg_of_nonpos_of_nonpos hil.le this
    · have e : sda i k x = -sda i x k := by unfold sda; ring
      rw [e, r1]; nlinarith
  · right
    have hD2 : sda j k l ≠ 0 := by rw [r2]; exact hkj.ne'
    refine (mem_hull3_iff hD2).2 ⟨?_, ?_, ?_⟩
    · rw [r2, r3]; exact mul_nonneg hkj.le hx0
    · have e : sda i j k * sda j x l = sda x j k * sda j i l + sda i j x * sda j k l := by
        unfold sda; ring
      have e' : sda j i l = -sda i j l := by unfold sda; ring
      have hjl : 0 < sda j k l := by rw [r2]; exact hkj
      have : 0 ≤ sda j x l := by
        by_contra hneg; push Not at hneg
        nlinarith [mul_nonneg y1 (neg_nonneg.2 hl.le), mul_nonneg y3 hjl.le, mul_pos hD (neg_pos.2 hneg)]
      exact mul_nonneg hjl.le this
    · have e : sda j k x = sda x j k := by unfold sda; ring
      rw [e, r2]; exact mul_nonneg hkj.le y1

theorem hull_subset_flip {i j k l : ℝ × ℝ} (hk : sda i j k * sda i j l < 0)
    (hc : sda k l i * sda k l j < 0) :
    convexHull ℝ {i, j, k} ⊆ convexHull ℝ {i, k, l} ∪ convexHull ℝ {j, k, l} := by
  have hD : sda i j k ≠ 0 := by rintro h0; simp [h0] at hk
  rcases hD.lt_or_gt with hneg | hpos
  · have hk' : sda j i k * sda j i l < 0 := by
      have : sda j i k * sda j i l = sda i j k * sda i j l := by unfold sda; ring
      linarith
    have hc' : sda k l j * sda k l i < 0 := by linarith [mul_comm (sda k l i) (sda k l j)]
    have hpos' : 0 < sda j i k := by rw [sda_swap_left]; linarith
    have := hull_subset_flip_pos hk' hc' hpos'
    have e : ({j, i, k} : Set (ℝ × ℝ)) = {i, j, k} := Set.insert_comm _ _ _
    rw [e] at this
    exact this.trans (Set.union_comm _ _).le
  · exact hull_subset_flip_pos hk hc hpos

theorem hull_rot (a b c : ℝ × ℝ) : convexHull ℝ {a, b, c} = convexHull ℝ {b, c, a} := by
  congr 1; ext; simp; tauto

theorem hull_swap (a b c : ℝ × ℝ) : convexHull ℝ {a, b, c} = convexHull ℝ {b, a, c} := by
  congr 1; ext; simp; tauto

theorem hull_swap₂₃ (a b c : ℝ × ℝ) : convexHull ℝ {a, b, c} = convexHull ℝ {a, c, b} := by
  congr 1; ext; simp; tauto

/-- **The two triangulations of a strictly convex quadrilateral cover the same set.** -/
theorem flip_hull_union {i j k l : ℝ × ℝ} (hk : sda i j k * sda i j l < 0)
    (hc : sda k l i * sda k l j < 0) :
    convexHull ℝ {i, k, l} ∪ convexHull ℝ {j, k, l} = convexHull ℝ {i, j, k} ∪ convexHull ℝ {i, j, l} := by
  have hk2 : sda i j l * sda i j k < 0 := by linarith [mul_comm (sda i j l) (sda i j k)]
  have hc2 : sda l k i * sda l k j < 0 := by
    have : sda l k i * sda l k j = sda k l i * sda k l j := by unfold sda; ring
    linarith
  have hc3 : sda k l j * sda k l i < 0 := by linarith [mul_comm (sda k l i) (sda k l j)]
  have hk3 : sda j i k * sda j i l < 0 := by
    have : sda j i k * sda j i l = sda i j k * sda i j l := by unfold sda; ring
    linarith
  apply Set.Subset.antisymm
  · apply Set.union_subset
    · have := hull_subset_flip hc hk
      rwa [hull_rot k l i, hull_rot l i k, hull_rot k i j, hull_rot l i j] at this
    · have := hull_subset_flip hc3 hk3
      rwa [hull_rot k l j, hull_rot l j k, hull_rot k j i, hull_swap j i k, hull_rot l j i,
        hull_swap j i l] at this
  · apply Set.union_subset
    · exact hull_subset_flip hk hc
    · have := hull_subset_flip hk2 hc2
      rwa [hull_swap₂₃ i l k, hull_swap₂₃ j l k] at this

/-! ## The flip -/

/-- **The flip of the interior edge `ij`** with opposite corners `k, l`: the triangles `{i,j,k}`,
`{i,j,l}` are replaced by `{i,k,l}`, `{j,k,l}`. -/
def flip (T : Finset (Finset (ℝ × ℝ))) (i j k l : ℝ × ℝ) : Finset (Finset (ℝ × ℝ)) :=
  insert {i, k, l} (insert {j, k, l} ((T.erase {i, j, k}).erase {i, j, l}))

theorem mem_flip {T : Finset (Finset (ℝ × ℝ))} {i j k l : ℝ × ℝ} {t : Finset (ℝ × ℝ)} :
    t ∈ flip T i j k l ↔ t = {i, k, l} ∨ t = {j, k, l} ∨ (t ≠ {i, j, l} ∧ t ≠ {i, j, k} ∧ t ∈ T) := by
  simp [flip, Finset.mem_erase]

/-- An open subset of the two old triangles that misses the interiors of both is empty; so a new
triangle (inside the two old ones) has interior disjoint from every other triangle. -/
theorem disjoint_interior_of_subset {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {t₁ t₂ u : Finset (ℝ × ℝ)} (h₁ : t₁ ∈ T) (h₂ : t₂ ∈ T) (hu : u ∈ T)
    (hu₁ : u ≠ t₁) (hu₂ : u ≠ t₂) {K : Set (ℝ × ℝ)} (hK : K ⊆ hull t₁ ∪ hull t₂) :
    Disjoint (interior K) (interior (hull u)) := by
  rw [Set.disjoint_iff_inter_eq_empty]
  refine eq_empty_of_subset_union_hull (convex_convexHull ℝ _) (convex_convexHull ℝ _)
    (isOpen_interior.inter isOpen_interior) (fun x hx => hK (interior_subset hx.1)) ?_ ?_
  · exact Set.disjoint_of_subset_left Set.inter_subset_right (hT.disjoint hu h₁ hu₁)
  · exact Set.disjoint_of_subset_left Set.inter_subset_right (hT.disjoint hu h₂ hu₂)

/-- **The Lawson flip.** In a triangulation of `S`, let `{i,j,k}` and `{i,j,l}` (`k ≠ l`) share the
edge `ij`, and let the quadrilateral be strictly convex (`k, l` are on opposite sides of `ij`
automatically, `IsTriangulation.opposite`; here also `i, j` on opposite sides of `kl`). Then
`flip T i j k l` is again a triangulation of `S`.

[formula: planar-lawson-flip] -/
theorem lawson_flip {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))} (hT : IsTriangulation S T)
    {i j k l : ℝ × ℝ} (h1 : {i, j, k} ∈ T) (h2 : {i, j, l} ∈ T) (hkl : k ≠ l)
    (hc : sda k l i * sda k l j < 0) : IsTriangulation S (flip T i j k l) := by
  have hk := hT.opposite h1 h2 hkl
  have hD2 := hT.sda_ne_zero h2
  have h12 := IsTriangulation.triple_ne hkl hD2
  have hni : sda k l i ≠ 0 := by rintro h0; simp [h0] at hc
  have hnj : sda k l j ≠ 0 := by rintro h0; simp [h0] at hc
  have r1 : sda i k l = sda k l i := by unfold sda; ring
  have r2 : sda j k l = sda k l j := by unfold sda; ring
  have hU := flip_hull_union hk hc
  have hn1 : hull {i, k, l} ⊆ hull {i, j, k} ∪ hull {i, j, l} := by
    rw [hull_triple, hull_triple, hull_triple, ← hU]; exact Set.subset_union_left
  have hn2 : hull {j, k, l} ⊆ hull {i, j, k} ∪ hull {i, j, l} := by
    rw [hull_triple, hull_triple, hull_triple, ← hU]; exact Set.subset_union_right
  have hS : ∀ x ∈ ({i, j, k} : Finset (ℝ × ℝ)) ∪ {i, j, l}, x ∈ S := by
    intro x hx
    rcases Finset.mem_union.1 hx with hx | hx
    · exact hT.mem_S h1 hx
    · exact hT.mem_S h2 hx
  have hT' : T = insert {i, j, k} (insert {i, j, l} ((T.erase {i, j, k}).erase {i, j, l})) := by
    rw [Finset.insert_erase (Finset.mem_erase.2 ⟨h12.symm, h2⟩), Finset.insert_erase h1]
  -- the interiors of the two new triangles are on opposite sides of `kl`
  have hnew : Disjoint (interior (hull {i, k, l})) (interior (hull {j, k, l})) := by
    rw [Set.disjoint_left]
    intro x hx1 hx2
    have e1 : hull {i, k, l} = convexHull ℝ {k, l, i} := by
      rw [hull_triple]; congr 1; ext; simp; tauto
    have e2 : hull {j, k, l} = convexHull ℝ {k, l, j} := by
      rw [hull_triple]; congr 1; ext; simp; tauto
    rw [e1] at hx1
    rw [e2] at hx2
    have p1 := (pos_of_mem_interior_hull3 hni hx1).2.2
    have p2 := (pos_of_mem_interior_hull3 hnj hx2).2.2
    have e : (sda k l i * sda k l x) * (sda k l j * sda k l x) =
        (sda k l i * sda k l j) * (sda k l x) ^ 2 := by ring
    have := mul_pos p1 p2
    nlinarith [mul_nonpos_of_nonpos_of_nonneg hc.le (sq_nonneg (sda k l x))]
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro t ht
    rcases mem_flip.1 ht with rfl | rfl | ⟨-, -, ht⟩
    · exact ⟨i, k, l, by rw [r1]; exact hni, rfl⟩
    · exact ⟨j, k, l, by rw [r2]; exact hnj, rfl⟩
    · exact hT.tri t ht
  · intro t ht
    rcases mem_flip.1 ht with rfl | rfl | ⟨-, -, ht⟩
    · intro x hx
      apply hS
      simp only [Finset.mem_insert, Finset.mem_singleton] at hx
      rcases hx with rfl | rfl | rfl <;> simp
    · intro x hx
      apply hS
      simp only [Finset.mem_insert, Finset.mem_singleton] at hx
      rcases hx with rfl | rfl | rfl <;> simp
    · exact hT.subset t ht
  · have hcov : ⋃ t ∈ T, hull t =
        hull {i, j, k} ∪ hull {i, j, l} ∪ ⋃ t ∈ (T.erase {i, j, k}).erase {i, j, l}, hull t := by
      conv_lhs => rw [hT']
      simp only [Finset.set_biUnion_insert, Set.union_assoc]
    rw [← hT.cover, hcov, flip]
    simp only [Finset.set_biUnion_insert]
    rw [← Set.union_assoc, hull_triple, hull_triple, hull_triple, hull_triple, hU]
  · have hR : ∀ u ∈ (T.erase {i, j, k}).erase {i, j, l}, u ∈ T ∧ u ≠ {i, j, k} ∧ u ≠ {i, j, l} := by
      intro u hu
      simp only [Finset.mem_erase] at hu
      exact ⟨hu.2.2, hu.2.1, hu.1⟩
    simp only [flip, Finset.coe_insert]
    rw [Set.pairwise_insert, Set.pairwise_insert]
    refine ⟨⟨?_, ?_⟩, ?_⟩
    · exact hT.disjoint.mono (fun u hu => (hR u hu).1)
    · intro u hu _
      obtain ⟨huT, hu1, hu2⟩ := hR u hu
      have := disjoint_interior_of_subset hT h1 h2 huT hu1 hu2 hn2
      exact ⟨this, this.symm⟩
    · intro u hu hne
      rcases hu with rfl | hu
      · exact ⟨hnew, hnew.symm⟩
      · obtain ⟨huT, hu1, hu2⟩ := hR u hu
        have := disjoint_interior_of_subset hT h1 h2 huT hu1 hu2 hn1
        exact ⟨this, this.symm⟩
  · intro t ht s hs hst
    have hnot : ∀ {a b : ℝ × ℝ}, sda k l a * sda k l b < 0 → a ∉ convexHull ℝ {k, l, b} := by
      intro a b hab ha
      have hb : sda k l b ≠ 0 := by rintro h0; simp [h0] at hab
      have := ((mem_hull3_iff hb).1 ha).2.2
      nlinarith [mul_comm (sda k l a) (sda k l b)]
    rcases mem_flip.1 ht with rfl | rfl | ⟨-, -, ht⟩
    · have hs' := hn1 hst
      have hsj : s ≠ j := by
        rintro rfl
        refine hnot (by linarith [mul_comm (sda k l i) (sda k l s)]) (?_ : s ∈ convexHull ℝ {k, l, i})
        rw [hull_triple] at hst
        have e : ({i, k, l} : Set (ℝ × ℝ)) = {k, l, i} := by ext; simp; tauto
        rwa [e] at hst
      rcases hs' with hs' | hs'
      · have := hT.vertex _ h1 s hs hs'
        simp only [Finset.mem_insert, Finset.mem_singleton] at this ⊢
        tauto
      · have := hT.vertex _ h2 s hs hs'
        simp only [Finset.mem_insert, Finset.mem_singleton] at this ⊢
        tauto
    · have hs' := hn2 hst
      have hsi : s ≠ i := by
        rintro rfl
        refine hnot hc (?_ : s ∈ convexHull ℝ {k, l, j})
        rw [hull_triple] at hst
        have e : ({j, k, l} : Set (ℝ × ℝ)) = {k, l, j} := by ext; simp; tauto
        rwa [e] at hst
      rcases hs' with hs' | hs'
      · have := hT.vertex _ h1 s hs hs'
        simp only [Finset.mem_insert, Finset.mem_singleton] at this ⊢
        tauto
      · have := hT.vertex _ h2 s hs hs'
        simp only [Finset.mem_insert, Finset.mem_singleton] at this ⊢
        tauto
    · exact hT.vertex t ht s hs hst

/-! ## Delaunay flips and the Delaunay properties -/

/-- **Locally Delaunay**: at every interior edge `ij` of `T` (triangles `{i,j,k}`, `{i,j,l}`,
`k ≠ l`) the fourth corner `l` is not strictly inside the circumcircle of `i, j, k`. -/
def IsLocallyDelaunay (T : Finset (Finset (ℝ × ℝ))) : Prop :=
  ∀ i j k l : ℝ × ℝ, {i, j, k} ∈ T → {i, j, l} ∈ T → k ≠ l → ¬ InCircumcircle i j k l

/-- **Delaunay** (empty circumcircles): no point of `S` lies strictly inside the circumcircle of a
triangle of `T`. -/
def IsDelaunay (S : Finset (ℝ × ℝ)) (T : Finset (Finset (ℝ × ℝ))) : Prop :=
  ∀ a b c : ℝ × ℝ, {a, b, c} ∈ T → ∀ p ∈ S, ¬ InCircumcircle a b c p

/-- A Delaunay triangulation is locally Delaunay. -/
theorem IsDelaunay.isLocallyDelaunay {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (h : IsDelaunay S T) : IsLocallyDelaunay T :=
  fun i j k l h1 h2 _ => h i j k h1 l (hT.mem_S h2 (by simp))

/-- **A Delaunay flip**: `T'` is obtained from `T` by flipping an interior edge `ij` whose fourth
corner `l` is strictly inside the circumcircle of `i, j, k`. -/
def DelaunayFlip (T T' : Finset (Finset (ℝ × ℝ))) : Prop :=
  ∃ i j k l : ℝ × ℝ, {i, j, k} ∈ T ∧ {i, j, l} ∈ T ∧ k ≠ l ∧ InCircumcircle i j k l ∧
    T' = flip T i j k l

/-- **A non-locally-Delaunay edge can always be flipped**: the result of a Delaunay flip of a
triangulation of `S` is a triangulation of `S`. -/
theorem isTriangulation_of_delaunayFlip {S : Finset (ℝ × ℝ)} {T T' : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (h : DelaunayFlip T T') : IsTriangulation S T' := by
  obtain ⟨i, j, k, l, h1, h2, hkl, hin, rfl⟩ := h
  exact lawson_flip hT h1 h2 hkl (crossing_of_inCircumcircle (hT.opposite h1 h2 hkl) hin)

/-- A triangulation admits no Delaunay flip exactly when it is locally Delaunay. -/
theorem isLocallyDelaunay_iff_forall_not_delaunayFlip {T : Finset (Finset (ℝ × ℝ))} :
    IsLocallyDelaunay T ↔ ∀ T', ¬ DelaunayFlip T T' := by
  constructor
  · rintro h T' ⟨i, j, k, l, h1, h2, hkl, hin, -⟩
    exact h i j k l h1 h2 hkl hin
  · intro h i j k l h1 h2 hkl hin
    exact h _ ⟨i, j, k, l, h1, h2, hkl, hin, rfl⟩

end Planar

end Delaunay

end DdgProofs
