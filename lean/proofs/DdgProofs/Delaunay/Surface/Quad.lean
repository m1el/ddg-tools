import DdgProofs.Delaunay.Surface.Segment

/-!
# Delaunay layer 4, planar toolkit: a strictly convex quadrilateral and its two diagonals

The diamond of a flip: points `P, Q` (the ends of the old diagonal) and `K, L` (the apexes) with
`K` to the left of `P → Q`, `L` to the right, and the quadrilateral strictly convex:
`sda P L K > 0` and `sda Q K L > 0` (the new triangles `P L K`, `Q K L` are positively
oriented; `Intrinsic.flipArea_realized`). Everything is an evaluation of `sda`, affine in each
argument (`sda_combo₃`, `sda_four`), at the corners.

* `cross_mem_segment` (**the old diagonal is crossed inside**): a point of the line `PQ` on a
  segment from the triangle `PQK` to the triangle `QPL` lies on the closed segment `PQ`;
* `hull_new_subset`, `int_new_subset` (**the new triangles lie in the old diamond**):
  `hull PLK ⊆ hull PQK ∪ hull QPL` and `int PLK ⊆ int PQK ∪ open PQ ∪ int QPL`;
* `openSegment_old_subset`: the old open diagonal lies in `int PLK ∪ open LK ∪ int QKL`;
* `symm`: the configuration `(L, K, P, Q)` is again a strictly convex diamond (the new diagonal
  `L → K` with `P` on its left and `Q` on its right), which gives the symmetric statements;
* `Tri.hull_rot`, `Tri.int_rot`: hull and interior do not depend on the numbering of the corners;
* `Tri.exists_ball_int_of_side`: near a point of the open side `m`, every point strictly on the
  inner side of that side is interior;
* `exists_sda_zero_openSegment`: a segment whose ends are on opposite sides of a line meets it
  in an interior point; `exists_mem_sda_ne` : a nonempty open set is not contained in a line.
-/

noncomputable section

namespace DdgProofs.Delaunay.Surface

open DdgProofs.Intrinsic (dist2_comm)
open DdgProofs.Delaunay.Planar
open Set Filter Topology

/-- The triangle with the corners `A, B, C`. -/
def tri3 (A B C : Pt) : Tri := ![A, B, C]

@[simp] theorem tri3_zero (A B C : Pt) : tri3 A B C 0 = A := rfl
@[simp] theorem tri3_one (A B C : Pt) : tri3 A B C 1 = B := rfl
@[simp] theorem tri3_two (A B C : Pt) : tri3 A B C 2 = C := rfl

theorem tri3_pos {A B C : Pt} (h : 0 < sda A B C) : (tri3 A B C).Pos := h

/-! ### Signs at the corners of `tri3` -/

theorem mem_hull_tri3 {A B C : Pt} (h : 0 < sda A B C) {x : Pt} :
    x ∈ (tri3 A B C).hull ↔ 0 ≤ sda B C x ∧ 0 ≤ sda C A x ∧ 0 ≤ sda A B x := by
  rw [Tri.mem_hull_iff (tri3_pos h), Fin.forall_fin_succ, Fin.forall_fin_two]
  exact Iff.rfl

theorem mem_int_tri3 {A B C : Pt} (h : 0 < sda A B C) {x : Pt} :
    x ∈ (tri3 A B C).int ↔ 0 < sda B C x ∧ 0 < sda C A x ∧ 0 < sda A B x := by
  rw [Tri.mem_int_iff (tri3_pos h), Fin.forall_fin_succ, Fin.forall_fin_two]
  exact Iff.rfl

theorem mem_openSegment_tri3 {A B C : Pt} (h : 0 < sda A B C) {x : Pt} :
    x ∈ openSegment ℝ A B ↔ sda A B x = 0 ∧ 0 < sda B C x ∧ 0 < sda C A x :=
  Tri.mem_side_iff (T := tri3 A B C) (k := 0) (tri3_pos h)

/-- `sda U V ·` at `x`, times the area, is the combination of its corner values with the
barycentric numerators of `x`. -/
theorem sda_mul_eq_bary (A B C U V x : Pt) :
    sda A B C * sda U V x =
      sda B C x * sda U V A + sda C A x * sda U V B + sda A B x * sda U V C := by
  unfold sda; ring

/-- An affine `sda U V ·` that is `≥ 0` at the corners is `≥ 0` on the closed triangle. -/
theorem sda_nonneg_of_hull {A B C U V x : Pt} (h : 0 < sda A B C) (hx : x ∈ (tri3 A B C).hull)
    (ha : 0 ≤ sda U V A) (hb : 0 ≤ sda U V B) (hc : 0 ≤ sda U V C) : 0 ≤ sda U V x := by
  rw [mem_hull_tri3 h] at hx
  have e := sda_mul_eq_bary A B C U V x
  have : 0 ≤ sda A B C * sda U V x := by
    rw [e]
    exact add_nonneg (add_nonneg (mul_nonneg hx.1 ha) (mul_nonneg hx.2.1 hb))
      (mul_nonneg hx.2.2 hc)
  exact (mul_nonneg_iff_of_pos_left h).1 this

/-- An affine `sda U V ·` that is `≥ 0` at the corners and `> 0` at one of them is `> 0` on the
open triangle. -/
theorem sda_pos_of_int {A B C U V x : Pt} (h : 0 < sda A B C) (hx : x ∈ (tri3 A B C).int)
    (ha : 0 ≤ sda U V A) (hb : 0 ≤ sda U V B) (hc : 0 ≤ sda U V C)
    (hs : 0 < sda U V A ∨ 0 < sda U V B ∨ 0 < sda U V C) : 0 < sda U V x := by
  rw [mem_int_tri3 h] at hx
  have e := sda_mul_eq_bary A B C U V x
  have h1 := mul_nonneg hx.1.le ha
  have h2 := mul_nonneg hx.2.1.le hb
  have h3 := mul_nonneg hx.2.2.le hc
  have : 0 < sda A B C * sda U V x := by
    rw [e]
    rcases hs with hs | hs | hs
    · have := mul_pos hx.1 hs; linarith
    · have := mul_pos hx.2.1 hs; linarith
    · have := mul_pos hx.2.2 hs; linarith
  exact (mul_pos_iff_of_pos_left h).1 this

/-- On an open segment, an affine `sda U V ·` that is `≥ 0` at the ends and `> 0` at one of
them is `> 0`. -/
theorem sda_pos_of_openSegment {A B U V y : Pt} (hy : y ∈ openSegment ℝ A B)
    (ha : 0 ≤ sda U V A) (hb : 0 ≤ sda U V B) (hs : 0 < sda U V A ∨ 0 < sda U V B) :
    0 < sda U V y := by
  obtain ⟨a, b, ha0, hb0, hab, rfl⟩ := hy
  rw [sda_combo₃ _ _ _ _ hab]
  rcases hs with hs | hs
  · nlinarith [mul_nonneg hb0.le hb]
  · nlinarith [mul_nonneg ha0.le ha]

/-- A point of the closed triangle on the line of the side `A B` is on the closed side. -/
theorem mem_segment_of_sda {A B C z : Pt} (h : 0 < sda A B C) (h1 : 0 ≤ sda B C z)
    (h2 : 0 ≤ sda C A z) (h0 : sda A B z = 0) : z ∈ segment ℝ A B := by
  have hc := sda_cramer h.ne' z
  have e1 : sda z B C = sda B C z := by unfold sda; ring
  have e2 : sda A z C = sda C A z := by unfold sda; ring
  have hs := sda_four A B C z
  rw [e1, e2, h0, zero_div, zero_smul, add_zero] at hc
  refine ⟨sda B C z / sda A B C, sda C A z / sda A B C, div_nonneg h1 h.le,
    div_nonneg h2 h.le, ?_, hc.symm⟩
  rw [← add_div, div_eq_one_iff_eq h.ne']
  rw [e1, e2, h0] at hs
  linarith

namespace Tri

/-- The hull of a triangle read from the corner `a`. -/
theorem hull_rot (T : Tri) (a : Fin 3) : T.hull = (tri3 (T a) (T (a + 1)) (T (a + 2))).hull := by
  fin_cases a
  · rfl
  · show convexHull ℝ {T 0, T 1, T 2} = convexHull ℝ {T 1, T 2, T 0}
    congr 1
    ext x
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    tauto
  · show convexHull ℝ {T 0, T 1, T 2} = convexHull ℝ {T 2, T 0, T 1}
    congr 1
    ext x
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    tauto

theorem int_rot (T : Tri) (a : Fin 3) : T.int = (tri3 (T a) (T (a + 1)) (T (a + 2))).int := by
  rw [int, int, hull_rot T a]

theorem side_rot (T : Tri) (a : Fin 3) : T.side a = openSegment ℝ (T a) (T (a + 1)) := rfl

/-- Near a point of the open side `m`, a point strictly inside the line of that side is
interior. -/
theorem exists_ball_int_of_side {T : Tri} (hT : T.Pos) {m : Fin 3} {y : Pt} (hy : y ∈ T.side m) :
    ∃ ε > 0, ∀ z, dist z y < ε → 0 < T.bary (m + 2) z → z ∈ T.int := by
  obtain ⟨_, h1, h2⟩ := (mem_side_iff hT).1 hy
  have hU : IsOpen {z : Pt | 0 < T.bary m z ∧ 0 < T.bary (m + 1) z} :=
    (isOpen_lt continuous_const (continuous_bary _)).inter
      (isOpen_lt continuous_const (continuous_bary _))
  obtain ⟨ε, hε, hball⟩ := Metric.isOpen_iff.1 hU y ⟨h1, h2⟩
  refine ⟨ε, hε, fun z hz h3 => ?_⟩
  have hz' : z ∈ Metric.ball y ε := hz
  obtain ⟨k1, k2⟩ := hball hz'
  rw [mem_int_iff hT]
  intro j
  rcases (by omega : j = m ∨ j = m + 1 ∨ j = m + 2) with rfl | rfl | rfl
  · exact k1
  · exact k2
  · exact h3

end Tri

/-- A segment whose ends are on opposite sides of the line `PQ` meets the line at an interior
point of the segment. -/
theorem exists_sda_zero_openSegment {P Q w y : Pt} (hw : sda P Q w < 0) (hy : 0 < sda P Q y) :
    ∃ z ∈ openSegment ℝ w y, sda P Q z = 0 := by
  have hd : sda P Q w - sda P Q y < 0 := by linarith
  have hθ0 : 0 < sda P Q w / (sda P Q w - sda P Q y) := div_pos_of_neg_of_neg hw hd
  have hθ1 : sda P Q w / (sda P Q w - sda P Q y) < 1 := by
    rw [div_lt_one_of_neg hd]; linarith
  refine ⟨(1 - sda P Q w / (sda P Q w - sda P Q y)) • w +
    (sda P Q w / (sda P Q w - sda P Q y)) • y,
    ⟨_, _, by linarith, hθ0, by ring, rfl⟩, ?_⟩
  rw [sda_combo₃ _ _ _ _ (by ring)]
  have hd' : sda P Q w - sda P Q y ≠ 0 := hd.ne
  field_simp
  ring

/-- A nonempty open set has a point off the line `PQ`. -/
theorem exists_mem_sda_ne {P Q : Pt} (hPQ : P ≠ Q) {U : Set Pt} (hU : IsOpen U) {z : Pt}
    (hz : z ∈ U) : ∃ z' ∈ U, sda P Q z' ≠ 0 := by
  by_cases h0 : sda P Q z = 0
  · have hc : 0 < sqDist2 P Q := sqDist2_pos hPQ
    obtain ⟨s, hsU, hs⟩ :=
      ((exists_pos_near (hU.mem_nhds hz) ((-(Q.2 - P.2), Q.1 - P.1) : Pt)).and
        (eventually_mem_nhdsWithin (a := (0 : ℝ)) (s := Ioi 0))).exists
    refine ⟨_, hsU, ?_⟩
    have e : sda P Q (z + s • ((-(Q.2 - P.2), Q.1 - P.1) : Pt)) =
        sda P Q z + s * sqDist2 P Q := by
      simp only [sda, sqDist2, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
        smul_eq_mul]
      ring
    rw [e, h0, zero_add]
    exact (mul_pos (mem_Ioi.1 hs) hc).ne'
  · exact ⟨z, hz, h0⟩

/-! ## The convex diamond -/

/-- A strictly convex diamond with the diagonal `P → Q`: `K` on the left, `L` on the right,
and the other diagonal `L K` inside (the triangles `P L K` and `Q K L` positively oriented). -/
structure ConvexQuad (P Q K L : Pt) : Prop where
  hK : 0 < sda P Q K
  hL : 0 < sda Q P L
  hc1 : 0 < sda P L K
  hc2 : 0 < sda Q K L

namespace ConvexQuad

variable {P Q K L : Pt} (h : ConvexQuad P Q K L)
include h

theorem hL' : sda P Q L < 0 := by
  have := h.hL; rw [sda_swap_left] at this; linarith

/-- The new diagonal `L → K` with `P` on the left and `Q` on the right is again a strictly convex
diamond. -/
theorem symm : ConvexQuad L K P Q := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [sda_rotate]; exact h.hc1
  · rw [sda_rotate]; exact h.hc2
  · rw [sda_rotate, sda_rotate]; exact h.hL
  · rw [sda_rotate, sda_rotate]; exact h.hK

theorem ne_PQ : P ≠ Q := sda_ne_zero_left h.hK.ne'
theorem ne_LK : L ≠ K := sda_ne_zero_left h.symm.hK.ne'

/-! ### Corner values -/

theorem pKPQ : 0 < sda K P Q := by
  have e : sda K P Q = sda P Q K := by unfold sda; ring
  rw [e]; exact h.hK

theorem pQKP : 0 < sda Q K P := by
  have e : sda Q K P = sda P Q K := by unfold sda; ring
  rw [e]; exact h.hK

theorem pKPL : 0 < sda K P L := by
  have e : sda K P L = sda P L K := by unfold sda; ring
  rw [e]; exact h.hc1

theorem pLKP : 0 < sda L K P := by
  have e : sda L K P = sda P L K := by unfold sda; ring
  rw [e]; exact h.hc1

theorem pLQP : 0 < sda L Q P := by
  have e : sda L Q P = sda Q P L := by unfold sda; ring
  rw [e]; exact h.hL

theorem pPLQ : 0 < sda P L Q := by
  have e : sda P L Q = sda Q P L := by unfold sda; ring
  rw [e]; exact h.hL

theorem pLQK : 0 < sda L Q K := by
  have e : sda L Q K = sda Q K L := by unfold sda; ring
  rw [e]; exact h.hc2

theorem pKLQ : 0 < sda K L Q := by
  have e : sda K L Q = sda Q K L := by unfold sda; ring
  rw [e]; exact h.hc2

/-- A point strictly inside the four outer sides is in the open old diamond. -/
theorem old_of_four {y : Pt} (h1 : 0 < sda K P y) (h2 : 0 < sda P L y) (h3 : 0 < sda L Q y)
    (h4 : 0 < sda Q K y) :
    y ∈ (tri3 P Q K).int ∪ openSegment ℝ P Q ∪ (tri3 Q P L).int := by
  rcases lt_trichotomy (sda P Q y) 0 with hn | h0 | hp
  · refine Or.inr ((mem_int_tri3 h.hL).2 ⟨h2, h3, ?_⟩)
    rw [sda_swap_left]; linarith
  · exact Or.inl (Or.inr ((mem_openSegment_tri3 h.hK).2 ⟨h0, h4, h1⟩))
  · exact Or.inl (Or.inl ((mem_int_tri3 h.hK).2 ⟨h4, h1, hp⟩))

/-- A point strictly inside the four outer sides is in the open new diamond. -/
theorem new_of_four {y : Pt} (h1 : 0 < sda K P y) (h2 : 0 < sda P L y) (h3 : 0 < sda L Q y)
    (h4 : 0 < sda Q K y) :
    y ∈ (tri3 P L K).int ∪ openSegment ℝ L K ∪ (tri3 Q K L).int := by
  rcases lt_trichotomy (sda L K y) 0 with hn | h0 | hp
  · refine Or.inr ((mem_int_tri3 h.hc2).2 ⟨?_, h3, h4⟩)
    rw [sda_swap_left]; linarith
  · exact Or.inl (Or.inr ((mem_openSegment_tri3 h.pLKP).2 ⟨h0, h1, h2⟩))
  · exact Or.inl (Or.inl ((mem_int_tri3 h.hc1).2 ⟨hp, h1, h2⟩))

/-- **The old diagonal is crossed inside**: a point of the line `PQ` on the closed segment
between a point of `PQK` and a point of `QPL` lies on the closed segment `PQ`. -/
theorem cross_mem_segment {y w z : Pt} (hy : y ∈ (tri3 P Q K).hull) (hw : w ∈ (tri3 Q P L).hull)
    (hz : z ∈ segment ℝ y w) (h0 : sda P Q z = 0) : z ∈ segment ℝ P Q := by
  have a := (mem_hull_tri3 h.hK).1 hy
  have hKP : 0 ≤ sda K P w :=
    sda_nonneg_of_hull h.hL hw h.pKPQ.le (sda_self₂₃ K P).ge h.pKPL.le
  have hQK : 0 ≤ sda Q K w :=
    sda_nonneg_of_hull h.hL hw (sda_self₁₃ Q K).ge h.pQKP.le h.hc2.le
  obtain ⟨s, t, hs, ht, hst, rfl⟩ := hz
  refine mem_segment_of_sda h.hK ?_ ?_ h0
  · rw [sda_combo₃ _ _ _ _ hst]
    exact add_nonneg (mul_nonneg hs a.1) (mul_nonneg ht hQK)
  · rw [sda_combo₃ _ _ _ _ hst]
    exact add_nonneg (mul_nonneg hs a.2.1) (mul_nonneg ht hKP)

/-- The new triangle `P L K` lies in the old diamond. -/
theorem hull_new_subset : (tri3 P L K).hull ⊆ (tri3 P Q K).hull ∪ (tri3 Q P L).hull := by
  intro y hy
  have a := (mem_hull_tri3 h.hc1).1 hy
  have hQK : 0 ≤ sda Q K y :=
    sda_nonneg_of_hull h.hc1 hy h.pQKP.le h.hc2.le (sda_self₂₃ Q K).ge
  have hLQ : 0 ≤ sda L Q y :=
    sda_nonneg_of_hull h.hc1 hy h.pLQP.le (sda_self₁₃ L Q).ge h.pLQK.le
  by_cases hp : 0 ≤ sda P Q y
  · exact Or.inl ((mem_hull_tri3 h.hK).2 ⟨hQK, a.2.1, hp⟩)
  · refine Or.inr ((mem_hull_tri3 h.hL).2 ⟨a.2.2, hLQ, ?_⟩)
    rw [sda_swap_left]; linarith

/-- The open new triangle `P L K` lies in the open old diamond. -/
theorem int_new_subset :
    (tri3 P L K).int ⊆ (tri3 P Q K).int ∪ openSegment ℝ P Q ∪ (tri3 Q P L).int := by
  intro y hy
  have a := (mem_int_tri3 h.hc1).1 hy
  have hQK : 0 < sda Q K y :=
    sda_pos_of_int h.hc1 hy h.pQKP.le h.hc2.le (sda_self₂₃ Q K).ge (Or.inl h.pQKP)
  have hLQ : 0 < sda L Q y :=
    sda_pos_of_int h.hc1 hy h.pLQP.le (sda_self₁₃ L Q).ge h.pLQK.le (Or.inl h.pLQP)
  exact h.old_of_four a.2.1 a.2.2 hLQ hQK

/-- The old open diagonal lies in the open new diamond. -/
theorem openSegment_old_subset :
    openSegment ℝ P Q ⊆ (tri3 P L K).int ∪ openSegment ℝ L K ∪ (tri3 Q K L).int := by
  intro y hy
  have a := (mem_openSegment_tri3 h.hK).1 hy
  have hy' : y ∈ openSegment ℝ Q P := by rw [openSegment_symm]; exact hy
  have b := (mem_openSegment_tri3 h.hL).1 hy'
  exact h.new_of_four a.2.2 b.2.1 b.2.2 a.2.1

/-- The other new triangle `Q K L` lies in the old diamond (by symmetry). -/
theorem hull_new_subset' : (tri3 Q K L).hull ⊆ (tri3 P Q K).hull ∪ (tri3 Q P L).hull := by
  intro y hy
  have a := (mem_hull_tri3 h.hc2).1 hy
  have hKP : 0 ≤ sda K P y :=
    sda_nonneg_of_hull h.hc2 hy h.pKPQ.le (sda_self₁₃ K P).ge h.pKPL.le
  have hPL : 0 ≤ sda P L y :=
    sda_nonneg_of_hull h.hc2 hy h.pPLQ.le h.hc1.le (sda_self₂₃ P L).ge
  by_cases hp : 0 ≤ sda P Q y
  · exact Or.inl ((mem_hull_tri3 h.hK).2 ⟨a.2.2, hKP, hp⟩)
  · refine Or.inr ((mem_hull_tri3 h.hL).2 ⟨hPL, a.2.1, ?_⟩)
    rw [sda_swap_left]; linarith

theorem int_new_subset' :
    (tri3 Q K L).int ⊆ (tri3 P Q K).int ∪ openSegment ℝ P Q ∪ (tri3 Q P L).int := by
  intro y hy
  have a := (mem_int_tri3 h.hc2).1 hy
  have hKP : 0 < sda K P y :=
    sda_pos_of_int h.hc2 hy h.pKPQ.le (sda_self₁₃ K P).ge h.pKPL.le (Or.inl h.pKPQ)
  have hPL : 0 < sda P L y :=
    sda_pos_of_int h.hc2 hy h.pPLQ.le h.hc1.le (sda_self₂₃ P L).ge (Or.inl h.pPLQ)
  exact h.old_of_four hKP hPL a.2.1 a.2.2

/-- The open old diamond lies in the closed... : every point of the open old diamond lies in
one of the open new triangles or on the new open diagonal. -/
theorem int_old_subset :
    (tri3 P Q K).int ∪ openSegment ℝ P Q ∪ (tri3 Q P L).int ⊆
      (tri3 P L K).int ∪ openSegment ℝ L K ∪ (tri3 Q K L).int := by
  rintro y ((hy | hy) | hy)
  · have a := (mem_int_tri3 h.hK).1 hy
    have hPL : 0 < sda P L y :=
      sda_pos_of_int h.hK hy (sda_self₁₃ P L).ge h.pPLQ.le h.hc1.le (Or.inr (Or.inl h.pPLQ))
    have hLQ : 0 < sda L Q y :=
      sda_pos_of_int h.hK hy h.pLQP.le (sda_self₂₃ L Q).ge h.pLQK.le (Or.inl h.pLQP)
    exact h.new_of_four a.2.1 hPL hLQ a.1
  · exact h.openSegment_old_subset hy
  · have a := (mem_int_tri3 h.hL).1 hy
    have hKP : 0 < sda K P y :=
      sda_pos_of_int h.hL hy h.pKPQ.le (sda_self₂₃ K P).ge h.pKPL.le (Or.inl h.pKPQ)
    have hQK : 0 < sda Q K y :=
      sda_pos_of_int h.hL hy (sda_self₁₃ Q K).ge h.pQKP.le h.hc2.le (Or.inr (Or.inl h.pQKP))
    exact h.new_of_four hKP a.1 a.2.1 hQK

/-- Points of the open side `P L` are strictly on `L`'s side of `PQ`. -/
theorem sda_neg_of_mem_PL {y : Pt} (hy : y ∈ openSegment ℝ P L) : sda P Q y < 0 := by
  have := sda_pos_of_openSegment (U := Q) (V := P) hy (sda_self₂₃ Q P).ge h.hL.le (Or.inr h.hL)
  rw [sda_swap_left] at this; linarith

/-- Points of the open side `Q K` are strictly on `K`'s side of `PQ`. -/
theorem sda_pos_of_mem_QK {y : Pt} (hy : y ∈ openSegment ℝ Q K) : 0 < sda P Q y :=
  sda_pos_of_openSegment hy (sda_self₂₃ P Q).ge h.hK.le (Or.inr h.hK)

/-- Points of the open side `K P` are strictly on `K`'s side of `PQ`. -/
theorem sda_pos_of_mem_KP {y : Pt} (hy : y ∈ openSegment ℝ K P) : 0 < sda P Q y :=
  sda_pos_of_openSegment hy h.hK.le (sda_self₁₃ P Q).ge (Or.inl h.hK)

/-- Points of the open side `L Q` are strictly on `L`'s side of `PQ`. -/
theorem sda_neg_of_mem_LQ {y : Pt} (hy : y ∈ openSegment ℝ L Q) : sda P Q y < 0 := by
  have := sda_pos_of_openSegment (U := Q) (V := P) hy h.hL.le (sda_self₁₃ Q P).ge (Or.inl h.hL)
  rw [sda_swap_left] at this; linarith

/-- The open new diagonal is in the open old diamond. -/
theorem openSegment_new_subset :
    openSegment ℝ L K ⊆ (tri3 P Q K).int ∪ openSegment ℝ P Q ∪ (tri3 Q P L).int := by
  intro y hy
  have a := (mem_openSegment_tri3 h.pLKP).1 hy
  have hy' : y ∈ openSegment ℝ K L := by rw [openSegment_symm]; exact hy
  have b := (mem_openSegment_tri3 h.pKLQ).1 hy'
  exact h.old_of_four a.2.1 a.2.2 b.2.1 b.2.2

end ConvexQuad

end DdgProofs.Delaunay.Surface
