import DDG.Geometry.Plane

/-!
# Integer coordinates: chords of a convex polygon, cevians (tier B)

`crates/ddg-integer-coords/SPEC.md` §9.1–9.2. The face model (`FaceModel.lean`) declares two arcs
disjoint when their ends do not interleave on the boundary cycle of the triangle. This file shows
that this is the geometric statement for straight arcs. Plane points are `ℝ × ℝ`, orientation is
the workspace's `sda` (twice the signed area, `signed-area-2d`), and `lerp a b s = a + s (b − a)`.

* **`ic-chord-interleave`** —
  * four points in strictly convex position in the cyclic order `a, b, c, d` (`Convex4`: the four
    orientations of consecutive triples are positive): the chords with interleaved ends, `ac` and
    `bd`, meet in exactly one point, which is interior to both (`Convex4.diagonals_cross`); the
    chords with non-interleaved ends, `ab, cd` and `ad, bc`, are disjoint as closed segments
    (`Convex4.sides_disjoint_ab_cd`, `Convex4.sides_disjoint_ad_bc`);
  * points on the boundary of a positively oriented triangle, given as `bpt V k u` (the point at
    parameter `u ∈ [0, 1)` of side `k`, from corner `k` to corner `k + 1`, i.e. the face model's
    convention with the corner `k` at `u = 0`): the orientation of three of them is the determinant
    of their barycentric rows times the triangle's (`sda_bpt`); three of them in strict boundary
    order (`BLt`, lexicographic in `(k, u)`) have a nonnegative orientation, which is positive unless
    the three lie on one closed side (`sda_bpt_pos`); hence four points in strict boundary order of
    which no three lie on one closed side are in strictly convex position (`convex4_bpt`), and the
    chord statement applies to them (`bpt_chords`): two straight arcs of a triangle with four
    distinct ends cross if and only if their ends interleave on the boundary.
* **`ic-cevians-cross`** — in a nondegenerate triangle the segments from two different corners to
  points of their open opposite sides meet at a point interior to both (`cevians_cross`), and a
  segment cutting corner `0` (from the open side `0` to the open side `2`) meets every segment from
  corner `0` to a point of the open opposite side (`corner_arc_crosses_cevian`): the geometric
  reason why at most one corner of a face emanates and why an emanating corner is cut by no arc
  (the combinatorial halves are `FaceModel.Valid.emanating_unique`,
  `FaceModel.Valid.cut_eq_zero_of_emanating`).
-/

namespace DDG

namespace NormalCoords

/-! ## Segments -/

/-- The point at parameter `s` on the segment from `a` to `b`. -/
def lerp (a b : ℝ × ℝ) (s : ℝ) : ℝ × ℝ := (a.1 + s * (b.1 - a.1), a.2 + s * (b.2 - a.2))

/-- `p` lies on the closed segment `ab`. -/
def OnSeg (a b p : ℝ × ℝ) : Prop := ∃ s : ℝ, 0 ≤ s ∧ s ≤ 1 ∧ p = lerp a b s

/-- `sda` is affine in a point of a segment. -/
theorem sda_lerp (c d a b : ℝ × ℝ) (s : ℝ) :
    sda c d (lerp a b s) = (1 - s) * sda c d a + s * sda c d b := by
  unfold sda lerp; ring

theorem sda_lerp_self (b d : ℝ × ℝ) (t : ℝ) : sda b d (lerp b d t) = 0 := by
  unfold sda lerp; ring

theorem sda_cycle (a b c : ℝ × ℝ) : sda b c a = sda a b c := by unfold sda; ring

theorem sda_swap (a b c : ℝ × ℝ) : sda a c b = -sda a b c := by unfold sda; ring

/-! ## Four points in strictly convex position -/

/-- `a, b, c, d` are in strictly convex position, counter-clockwise in this cyclic order. -/
def Convex4 (a b c d : ℝ × ℝ) : Prop :=
  0 < sda a b c ∧ 0 < sda b c d ∧ 0 < sda c d a ∧ 0 < sda d a b

/-- Points of a closed segment lie on the positive side of a line that has both ends there. -/
theorem sda_pos_of_onSeg {c d a b p : ℝ × ℝ} (ha : 0 < sda c d a) (hb : 0 < sda c d b)
    (hp : OnSeg a b p) : 0 < sda c d p := by
  obtain ⟨s, h0, h1, rfl⟩ := hp
  rw [sda_lerp]
  rcases eq_or_lt_of_le h0 with rfl | h0'
  · simpa using ha
  · nlinarith

/-- **Non-interleaved chords are disjoint** (`ab` and `cd`). [formula: ic-chord-interleave] -/
theorem Convex4.sides_disjoint_ab_cd {a b c d : ℝ × ℝ} (h : Convex4 a b c d) :
    ¬ ∃ p, OnSeg a b p ∧ OnSeg c d p := by
  rintro ⟨p, hab, ⟨t, -, -, rfl⟩⟩
  have hpos := sda_pos_of_onSeg (c := c) (d := d) h.2.2.1
    (by rw [sda_cycle]; exact h.2.1) hab
  rw [sda_lerp_self] at hpos
  exact lt_irrefl _ hpos

/-- **Non-interleaved chords are disjoint** (`ad` and `bc`). [formula: ic-chord-interleave] -/
theorem Convex4.sides_disjoint_ad_bc {a b c d : ℝ × ℝ} (h : Convex4 a b c d) :
    ¬ ∃ p, OnSeg a d p ∧ OnSeg b c p := by
  rintro ⟨p, had, ⟨t, -, -, rfl⟩⟩
  have hpos := sda_pos_of_onSeg (c := b) (d := c) (by rw [sda_cycle]; exact h.1) h.2.1 had
  rw [sda_lerp_self] at hpos
  exact lt_irrefl _ hpos

/-- The parameter of the crossing on the chord `ac`. -/
theorem lerp_cross_eq {a b c d : ℝ × ℝ} (hQ : sda d a b + sda b c d ≠ 0) :
    lerp a c (sda d a b / (sda d a b + sda b c d)) =
      lerp b d (sda a b c / (sda a b c + sda c d a)) := by
  have e : sda a b c + sda c d a = sda d a b + sda b c d := by unfold sda; ring
  rw [e]
  unfold lerp
  ext <;> simp only <;> field_simp <;> unfold sda <;> ring

/-- **Interleaved chords cross**: the diagonals `ac` and `bd` of a strictly convex quadrilateral
meet in exactly one point, at parameters strictly between `0` and `1` on both.
[formula: ic-chord-interleave] -/
theorem Convex4.diagonals_cross {a b c d : ℝ × ℝ} (h : Convex4 a b c d) :
    (∃ s t : ℝ, 0 < s ∧ s < 1 ∧ 0 < t ∧ t < 1 ∧ lerp a c s = lerp b d t) ∧
      ∀ p q, OnSeg a c p → OnSeg b d p → OnSeg a c q → OnSeg b d q → p = q := by
  obtain ⟨h1, h2, h3, h4⟩ := h
  refine ⟨⟨_, _, div_pos h4 (by linarith), (div_lt_one (by linarith)).2 (by linarith),
    div_pos h1 (by linarith), (div_lt_one (by linarith)).2 (by linarith),
    lerp_cross_eq (by linarith)⟩, ?_⟩
  -- a point of `ac` on the line `bd` has a determined parameter
  have key : ∀ p, OnSeg a c p → OnSeg b d p → p = lerp a c (sda d a b / (sda d a b + sda b c d)) := by
    rintro p ⟨s, -, -, rfl⟩ ⟨t, -, -, ht⟩
    have e := congrArg (sda b d) ht
    rw [sda_lerp_self, sda_lerp] at e
    have ea : sda b d a = sda d a b := by unfold sda; ring
    have ec : sda b d c = -sda b c d := by unfold sda; ring
    rw [ea, ec] at e
    have hs : s = sda d a b / (sda d a b + sda b c d) := by
      rw [eq_div_iff (by linarith)]; linarith
    rw [hs]
  intro p q hp1 hp2 hq1 hq2
  rw [key p hp1 hp2, key q hq1 hq2]

/-! ## Points on the boundary of a triangle -/

/-- The point at parameter `u` of side `k` of the triangle `V` (side `k` runs from corner `k` to
corner `k + 1`; `u = 0` is corner `k`). -/
def bpt (V : Fin 3 → ℝ × ℝ) (k : Fin 3) (u : ℝ) : ℝ × ℝ := lerp (V k) (V (k + 1)) u

/-- The barycentric row of `bpt V k u`. -/
def brow (k : Fin 3) (u : ℝ) : ℝ × ℝ × ℝ :=
  if k = 0 then (1 - u, u, 0) else if k = 1 then (0, 1 - u, u) else (u, 0, 1 - u)

/-- Determinant of three barycentric rows. -/
def det3 (p q r : ℝ × ℝ × ℝ) : ℝ :=
  p.1 * (q.2.1 * r.2.2 - q.2.2 * r.2.1) - p.2.1 * (q.1 * r.2.2 - q.2.2 * r.1) +
    p.2.2 * (q.1 * r.2.1 - q.2.1 * r.1)

/-- The point with barycentric row `w`. -/
def baryPt (V : Fin 3 → ℝ × ℝ) (w : ℝ × ℝ × ℝ) : ℝ × ℝ :=
  (w.1 * (V 0).1 + w.2.1 * (V 1).1 + w.2.2 * (V 2).1,
    w.1 * (V 0).2 + w.2.1 * (V 1).2 + w.2.2 * (V 2).2)

theorem bpt_eq_baryPt (V : Fin 3 → ℝ × ℝ) (k : Fin 3) (u : ℝ) :
    bpt V k u = baryPt V (brow k u) := by
  fin_cases k <;> simp [bpt, brow, baryPt, lerp] <;> constructor <;> ring_nf

/-- Orientation of three points given by barycentric rows that sum to one. -/
theorem sda_baryPt (V : Fin 3 → ℝ × ℝ) (p q r : ℝ × ℝ × ℝ) (hp : p.1 + p.2.1 + p.2.2 = 1)
    (hq : q.1 + q.2.1 + q.2.2 = 1) (hr : r.1 + r.2.1 + r.2.2 = 1) :
    sda (baryPt V p) (baryPt V q) (baryPt V r) = det3 p q r * sda (V 0) (V 1) (V 2) := by
  obtain ⟨p0, p1, p2⟩ := p
  obtain ⟨q0, q1, q2⟩ := q
  obtain ⟨r0, r1, r2⟩ := r
  simp only at hp hq hr
  have ep : p0 = 1 - p1 - p2 := by linarith
  have eq : q0 = 1 - q1 - q2 := by linarith
  have er : r0 = 1 - r1 - r2 := by linarith
  subst ep eq er
  unfold sda baryPt det3; ring

theorem brow_sum (k : Fin 3) (u : ℝ) : (brow k u).1 + (brow k u).2.1 + (brow k u).2.2 = 1 := by
  fin_cases k <;> simp [brow]

/-- **Orientation of boundary points** is the determinant of their barycentric rows times the
orientation of the triangle. [formula: ic-chord-interleave] -/
theorem sda_bpt (V : Fin 3 → ℝ × ℝ) (k1 k2 k3 : Fin 3) (u1 u2 u3 : ℝ) :
    sda (bpt V k1 u1) (bpt V k2 u2) (bpt V k3 u3) =
      det3 (brow k1 u1) (brow k2 u2) (brow k3 u3) * sda (V 0) (V 1) (V 2) := by
  rw [bpt_eq_baryPt, bpt_eq_baryPt, bpt_eq_baryPt]
  exact sda_baryPt V _ _ _ (brow_sum _ _) (brow_sum _ _) (brow_sum _ _)

/-- Strict boundary order: lexicographic in `(side, parameter)`. -/
def BLt (k1 : Fin 3) (u1 : ℝ) (k2 : Fin 3) (u2 : ℝ) : Prop := k1 < k2 ∨ (k1 = k2 ∧ u1 < u2)

/-- Three boundary points (in boundary order) lie on one closed side: all on side `j`, or two on
side `j` and the third the corner `j + 1` (for `j = 2` that corner, `V 0`, comes first). -/
def OnOneSide (k1 : Fin 3) (u1 : ℝ) (k2 : Fin 3) (_u2 : ℝ) (k3 : Fin 3) (u3 : ℝ) : Prop :=
  (k1 = k2 ∧ k2 = k3) ∨ (k1 = k2 ∧ k3 = k1 + 1 ∧ k1 ≠ 2 ∧ u3 = 0) ∨
    (k1 = 0 ∧ u1 = 0 ∧ k2 = 2 ∧ k3 = 2)

/-- **Three boundary points in boundary order are positively oriented unless they lie on one
closed side** (parameters in `[0, 1)`). [formula: ic-chord-interleave] -/
theorem det3_brow_pos {k1 k2 k3 : Fin 3} {u1 u2 u3 : ℝ} (h1 : 0 ≤ u1) (h1' : u1 < 1)
    (h2 : 0 ≤ u2) (h2' : u2 < 1) (h3 : 0 ≤ u3) (h3' : u3 < 1) (o12 : BLt k1 u1 k2 u2)
    (o23 : BLt k2 u2 k3 u3) (hn : ¬ OnOneSide k1 u1 k2 u2 k3 u3) :
    0 < det3 (brow k1 u1) (brow k2 u2) (brow k3 u3) := by
  unfold BLt at o12 o23
  unfold OnOneSide at hn
  fin_cases k1 <;> fin_cases k2 <;> fin_cases k3 <;>
    simp [brow, det3] at o12 o23 hn ⊢ <;>
    first
    | nlinarith [mul_pos (sub_pos.2 h1') (sub_pos.2 h2'), mul_nonneg h1 h2,
        mul_pos (mul_pos (sub_pos.2 h1') (sub_pos.2 h2')) (sub_pos.2 h3'),
        mul_nonneg (mul_nonneg h1 h2) h3]
    | (have h3p := lt_of_le_of_ne h3 (Ne.symm hn); nlinarith [mul_pos h3p (sub_pos.2 o12)])
    | (have h1p := lt_of_le_of_ne h1 (Ne.symm hn); nlinarith [mul_pos h1p (sub_pos.2 o23)])

theorem BLt.trans {k1 k2 k3 : Fin 3} {u1 u2 u3 : ℝ} (h12 : BLt k1 u1 k2 u2)
    (h23 : BLt k2 u2 k3 u3) : BLt k1 u1 k3 u3 := by
  unfold BLt at *
  rcases h12 with h | ⟨rfl, h⟩ <;> rcases h23 with h' | ⟨rfl, h'⟩
  · exact Or.inl (h.trans h')
  · exact Or.inl h
  · exact Or.inl h'
  · exact Or.inr ⟨rfl, h.trans h'⟩

/-- **Four boundary points in boundary order, no three on one closed side, are in strictly convex
position** (for a positively oriented triangle; parameters in `[0, 1)`).
[formula: ic-chord-interleave] -/
theorem convex4_bpt (V : Fin 3 → ℝ × ℝ) (hV : 0 < sda (V 0) (V 1) (V 2)) {k1 k2 k3 k4 : Fin 3}
    {u1 u2 u3 u4 : ℝ} (h1 : 0 ≤ u1 ∧ u1 < 1) (h2 : 0 ≤ u2 ∧ u2 < 1) (h3 : 0 ≤ u3 ∧ u3 < 1)
    (h4 : 0 ≤ u4 ∧ u4 < 1) (o12 : BLt k1 u1 k2 u2) (o23 : BLt k2 u2 k3 u3) (o34 : BLt k3 u3 k4 u4)
    (n123 : ¬ OnOneSide k1 u1 k2 u2 k3 u3) (n124 : ¬ OnOneSide k1 u1 k2 u2 k4 u4)
    (n134 : ¬ OnOneSide k1 u1 k3 u3 k4 u4) (n234 : ¬ OnOneSide k2 u2 k3 u3 k4 u4) :
    Convex4 (bpt V k1 u1) (bpt V k2 u2) (bpt V k3 u3) (bpt V k4 u4) := by
  have pos : ∀ {a b c : Fin 3} {x y z : ℝ}, 0 ≤ x ∧ x < 1 → 0 ≤ y ∧ y < 1 → 0 ≤ z ∧ z < 1 →
      BLt a x b y → BLt b y c z → ¬ OnOneSide a x b y c z → 0 < sda (bpt V a x) (bpt V b y) (bpt V c z) :=
    fun hx hy hz o o' n => by
      rw [sda_bpt]; exact mul_pos (det3_brow_pos hx.1 hx.2 hy.1 hy.2 hz.1 hz.2 o o' n) hV
  refine ⟨pos h1 h2 h3 o12 o23 n123, pos h2 h3 h4 o23 o34 n234, ?_, ?_⟩
  · rw [sda_cycle]; exact pos h1 h3 h4 (o12.trans o23) o34 n134
  · rw [← sda_cycle]; exact pos h1 h2 h4 o12 (o23.trans o34) n124

/-- **Straight arcs of a triangle cross exactly when their ends interleave**: for four boundary
points in boundary order, no three on one closed side, the chords `p₁p₃` and `p₂p₄` meet in exactly
one point, interior to both, while `p₁p₂, p₃p₄` and `p₁p₄, p₂p₃` are disjoint.
[formula: ic-chord-interleave] -/
theorem bpt_chords (V : Fin 3 → ℝ × ℝ) (hV : 0 < sda (V 0) (V 1) (V 2)) {k1 k2 k3 k4 : Fin 3}
    {u1 u2 u3 u4 : ℝ} (h1 : 0 ≤ u1 ∧ u1 < 1) (h2 : 0 ≤ u2 ∧ u2 < 1) (h3 : 0 ≤ u3 ∧ u3 < 1)
    (h4 : 0 ≤ u4 ∧ u4 < 1) (o12 : BLt k1 u1 k2 u2) (o23 : BLt k2 u2 k3 u3) (o34 : BLt k3 u3 k4 u4)
    (n123 : ¬ OnOneSide k1 u1 k2 u2 k3 u3) (n124 : ¬ OnOneSide k1 u1 k2 u2 k4 u4)
    (n134 : ¬ OnOneSide k1 u1 k3 u3 k4 u4) (n234 : ¬ OnOneSide k2 u2 k3 u3 k4 u4) :
    let p1 := bpt V k1 u1; let p2 := bpt V k2 u2; let p3 := bpt V k3 u3; let p4 := bpt V k4 u4
    ((∃ s t : ℝ, 0 < s ∧ s < 1 ∧ 0 < t ∧ t < 1 ∧ lerp p1 p3 s = lerp p2 p4 t) ∧
      ∀ p q, OnSeg p1 p3 p → OnSeg p2 p4 p → OnSeg p1 p3 q → OnSeg p2 p4 q → p = q) ∧
    (¬ ∃ p, OnSeg p1 p2 p ∧ OnSeg p3 p4 p) ∧ (¬ ∃ p, OnSeg p1 p4 p ∧ OnSeg p2 p3 p) := by
  intro p1 p2 p3 p4
  have hc := convex4_bpt V hV h1 h2 h3 h4 o12 o23 o34 n123 n124 n134 n234
  exact ⟨hc.diagonals_cross, hc.sides_disjoint_ab_cd, hc.sides_disjoint_ad_bc⟩

theorem bpt_zero (V : Fin 3 → ℝ × ℝ) (k : Fin 3) : bpt V k 0 = V k := by
  simp [bpt, lerp]

/-- **Two cevians cross**: in a positively oriented triangle, the segment from corner `0` to the
point `X` at parameter `u ∈ (0, 1)` of the opposite side `1` (from corner `1` to corner `2`) and the
segment from corner `1` to the point `Y` at parameter `w ∈ (0, 1)` of its opposite side `2` meet at
one point, interior to both. [formula: ic-cevians-cross] -/
theorem cevians_cross (V : Fin 3 → ℝ × ℝ) (hV : 0 < sda (V 0) (V 1) (V 2)) {u w : ℝ}
    (hu : 0 < u ∧ u < 1) (hw : 0 < w ∧ w < 1) :
    ∃ s t : ℝ, 0 < s ∧ s < 1 ∧ 0 < t ∧ t < 1 ∧
      lerp (V 0) (bpt V 1 u) s = lerp (V 1) (bpt V 2 w) t := by
  have := (bpt_chords V hV (k1 := 0) (k2 := 1) (k3 := 1) (k4 := 2) (u1 := 0) (u2 := 0)
    (u3 := u) (u4 := w) ⟨le_rfl, one_pos⟩ ⟨le_rfl, one_pos⟩ ⟨hu.1.le, hu.2⟩ ⟨hw.1.le, hw.2⟩
    (Or.inl (by decide)) (Or.inr ⟨rfl, hu.1⟩) (Or.inl (by decide))
    (by simp [OnOneSide]) (by simp [OnOneSide]) (by simp [OnOneSide])
    (by simp [OnOneSide]; exact hw.1.ne')).1.1
  simpa only [bpt_zero] using this

/-- **An arc that cuts a corner crosses every arc from that corner to the opposite side**: in a
positively oriented triangle, the segment from the point `P` at parameter `u ∈ (0, 1)` of side `0`
to the point `Q` at parameter `w ∈ (0, 1)` of side `2` (cutting corner `0`) and the segment from
corner `0` to the point `X` at parameter `v ∈ (0, 1)` of side `1` meet at one point, interior to
both. [formula: ic-cevians-cross] -/
theorem corner_arc_crosses_cevian (V : Fin 3 → ℝ × ℝ) (hV : 0 < sda (V 0) (V 1) (V 2))
    {u v w : ℝ} (hu : 0 < u ∧ u < 1) (hv : 0 < v ∧ v < 1) (hw : 0 < w ∧ w < 1) :
    ∃ s t : ℝ, 0 < s ∧ s < 1 ∧ 0 < t ∧ t < 1 ∧
      lerp (V 0) (bpt V 1 v) s = lerp (bpt V 0 u) (bpt V 2 w) t := by
  have := (bpt_chords V hV (k1 := 0) (k2 := 0) (k3 := 1) (k4 := 2) (u1 := 0) (u2 := u)
    (u3 := v) (u4 := w) ⟨le_rfl, one_pos⟩ ⟨hu.1.le, hu.2⟩ ⟨hv.1.le, hv.2⟩ ⟨hw.1.le, hw.2⟩
    (Or.inr ⟨rfl, hu.1⟩) (Or.inl (by decide)) (Or.inl (by decide))
    (by simp [OnOneSide]; exact hv.1.ne') (by simp [OnOneSide]) (by simp [OnOneSide])
    (by simp [OnOneSide])).1.1
  simpa only [bpt_zero] using this

end NormalCoords

end DDG
