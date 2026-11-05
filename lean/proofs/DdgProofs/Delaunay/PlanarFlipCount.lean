import DdgProofs.Delaunay.PlanarUniqueness

/-!
# Planar Delaunay triangulations, part 6: at most `C(n, 2)` Lawson flips

**The flip count.** Every sequence of Delaunay flips starting from a triangulation of a finite
point set `S` with `n` points has at most `n(n − 1)/2 = C(n, 2)` flips (`flip_count_le`), because
**an edge that has been flipped away never reappears** (`flipped_edge_not_mem`), and every flip
removes an edge `{i, j}` with `i, j ∈ S`. This is the polynomial bound of Lawson's algorithm
(Edelsbrunner, *Geometry and Topology for Mesh Generation*, Sec. 1.3; de Berg et al., Ch. 9,
notes and comments, after Lawson 1977 and Sibson 1978); `flip_sequence_length_le` of
`DdgProofs/Delaunay/PlanarTermination.lean` bounded the length only by the number of sets of
candidate triangles.

## The proof: the lifted surface goes down

Lift every triangle `t` of `T` to its lifted plane `liftPlane t` (`PlanarUniqueness.lean`; the
plane through the lifted corners `(v, |v|²)`). Instead of the piecewise linear lifted function on
the hull (which would need the agreement of neighbouring planes on shared edges everywhere), the
proof uses the relation **`LiftBelow T T'`**: for every point `x` and every triangle `t ∈ T`
containing `x` there is a triangle `t' ∈ T'` containing `x` whose lifted plane at `x` is not
higher, `liftPlane t' x ≤ liftPlane t x`. It is reflexive and transitive, and

* **a Delaunay flip lowers the lifted surface** (`liftBelow_flip`): the quadrilateral `i k j l`
  is covered by both pairs of triangles (`flip_hull_union`); a new triangle, say `i k l`, and an
  old one, say `i j k`, share the corners `i, k`, so their planes differ by an affine function
  that vanishes on the line `ik` and is positive at `l` (`l` is inside the circumcircle of
  `i j k`): on `i k l` the new plane is below (`liftPlane_le_of_mem`, from
  `liftPlane_sub_mul`). Along any sequence of Delaunay flips, `LiftBelow` holds between the first
  and the last triangulation (`liftBelow_of_reach`).
* **just after the flip of `ij` the surface is strictly below the lifted segment `îĵ`** at the
  crossing point `y = (1 − σ) i + σ j` of the diagonals (`exists_crossing`, `0 < σ < 1`): the new
  triangle `i k l` contains `y` and its plane at `y` is strictly below the plane of `i j k`, whose
  value at `y` is `(1 − σ)|i|² + σ|j|²` (`liftPlane_segment`).
* **while `ij` is an edge, the surface at `y` is the lifted segment**
  (`liftPlane_eq_of_mem_edge`): `y` is a relative interior point of the edge and not a point of
  `S`, so a triangle containing `y` is one of the two triangles on the edge (the adjacency lemma
  `eq_insert_of_mem_edge` of `PlanarLemma.lean`), and both planes take the value
  `(1 − σ)|i|² + σ|j|²` at `y`.

So if `ij` were an edge again later, some triangle containing `y` would have its plane at `y`
strictly below the lifted segment and equal to it — impossible. The flipped edges of a sequence
are therefore distinct two-element subsets of `S`.

## References

* C. L. Lawson, *Software for C¹ surface interpolation*, in: Mathematical Software III, Academic
  Press, 1977, 161–194.
* R. Sibson, *Locally equiangular triangulations*, The Computer Journal 21(3), 1978, 243–245.
* H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge University Press, 2001,
  Sec. 1.3 (the flip algorithm takes at most `C(n, 2)` flips: the lifted surface moves down, an
  edge once flipped away never comes back).
* M. de Berg, O. Cheong, M. van Kreveld, M. Overmars, *Computational Geometry: Algorithms and
  Applications*, 3rd ed., Springer, 2008, Ch. 9 (legal triangulations; notes and comments).

The bibliographic data were not checked against the works here.
-/

noncomputable section

namespace DDG

namespace Delaunay

namespace Planar

open Set Relation

/-! ## Two lifted planes on a common edge -/

/-- The difference of the lifted planes of `p q s` and `p q r` (common edge `pq`) is an affine
function vanishing on the line `pq`, with the value `circleLift p q s r` at `r`:
`sda p q r · (L_pqs(x) − L_pqr(x)) = sda p q x · circleLift p q s r`. -/
theorem liftPlane_sub_mul {p q r s : ℝ × ℝ} (hr : sda p q r ≠ 0) (hs : sda p q s ≠ 0)
    (x : ℝ × ℝ) :
    sda p q r * (liftPlane p q s x - liftPlane p q r x) = sda p q x * circleLift p q s r := by
  have e1 := liftPlane_bary p q s p q r x
  have e2 := liftPlane_bary p q r p q r x
  rw [liftPlane_self₁ hs, liftPlane_self₂ hs] at e1
  rw [liftPlane_self₁ hr, liftPlane_self₂ hr, liftPlane_self₃ hr] at e2
  rw [circleLift_eq_liftPlane_sub hs]
  linear_combination e1 - e2

/-- **On the new triangle the new plane is below the old one.** If `r` is strictly inside the
circumcircle of `p q s`, then on the triangle `p q r` the lifted plane of `p q r` is not above
that of `p q s`. -/
theorem liftPlane_le_of_mem {p q r s x : ℝ × ℝ} (hr : sda p q r ≠ 0) (hs : sda p q s ≠ 0)
    (hin : InCircumcircle p q s r) (hx : x ∈ convexHull ℝ {p, q, r}) :
    liftPlane p q r x ≤ liftPlane p q s x := by
  have e := liftPlane_sub_mul hr hs x
  have hc := (circleLift_pos_iff hs r).2 hin
  have h3 := ((mem_hull3_iff hr).1 hx).2.2
  have hDD : 0 < sda p q r * sda p q r := mul_self_pos.2 hr
  have h : (sda p q r * sda p q r) * (liftPlane p q s x - liftPlane p q r x) =
      (sda p q r * sda p q x) * circleLift p q s r := by rw [mul_assoc, e]; ring
  have : 0 ≤ (sda p q r * sda p q r) * (liftPlane p q s x - liftPlane p q r x) := by
    rw [h]; exact mul_nonneg h3 hc.le
  have := (mul_nonneg_iff_of_pos_left hDD).1 this
  linarith

/-- Strictly below off the common edge. -/
theorem liftPlane_lt_of_pos {p q r s x : ℝ × ℝ} (hr : sda p q r ≠ 0) (hs : sda p q s ≠ 0)
    (hin : InCircumcircle p q s r) (hx : 0 < sda p q r * sda p q x) :
    liftPlane p q r x < liftPlane p q s x := by
  have e := liftPlane_sub_mul hr hs x
  have hc := (circleLift_pos_iff hs r).2 hin
  have h : (sda p q r * sda p q r) * (liftPlane p q s x - liftPlane p q r x) =
      (sda p q r * sda p q x) * circleLift p q s r := by rw [mul_assoc, e]; ring
  have : 0 < (sda p q r * sda p q r) * (liftPlane p q s x - liftPlane p q r x) := by
    rw [h]; exact mul_pos hx hc
  have := pos_of_mul_pos_right this (mul_self_nonneg _)
  linarith

/-- **The lifted plane on an edge** is the lifted segment: at `(1 − σ) i + σ j` the plane of a
triangle `i j c` takes the value `(1 − σ)|i|² + σ|j|²`. -/
theorem liftPlane_segment {i j c : ℝ × ℝ} (h : sda i j c ≠ 0) (σ : ℝ) :
    liftPlane i j c ((1 - σ) • i + σ • j) = (1 - σ) * nsq i + σ * nsq j := by
  have h1 : (1 - σ) + σ = 1 := by ring
  unfold liftPlane
  rw [sda_combo₁ _ _ _ _ h1, sda_combo₂ _ _ _ _ h1, sda_combo₃ _ _ _ _ h1]
  simp only [sda_self₁₂, sda_self₁₃, sda_self₂₃]
  field_simp
  ring

/-! ## The quadrilateral of a Delaunay flip -/

theorem inCircle_swap₃₄ (a b c d : ℝ × ℝ) : inCircle a b d c = -inCircle a b c d := by
  unfold inCircle sda sqDist2; ring

/-- For two triangles on opposite sides of `ij`, `l` is inside the circumcircle of `i j k` iff `k`
is inside that of `i j l` (one direction). -/
theorem inCircumcircle_opposite {i j k l : ℝ × ℝ} (hk : sda i j k * sda i j l < 0)
    (hin : InCircumcircle i j k l) : InCircumcircle i j l k := by
  unfold InCircumcircle at hin ⊢
  rw [inCircle_swap₃₄]
  have h1 : 0 < (sda i j k * sda i j k) * (sda i j l * -inCircle i j k l) := by
    have e : (sda i j k * sda i j k) * (sda i j l * -inCircle i j k l) =
        -(sda i j k * sda i j l) * (sda i j k * inCircle i j k l) := by ring
    rw [e]; exact mul_pos (neg_pos.2 hk) hin
  exact pos_of_mul_pos_right h1 (mul_self_nonneg _)

/-- Signs in a strictly convex quadrilateral `i k j l` (diagonals `ij` and `kl`): `j` and `l` lie
on the same side of the side `ik` (`sda i k l · sda i k j > 0`), and `j` and `k` on the same side
of the side `il` (`sda i k l · sda i j l > 0`), written with `sda i k l = sda k l i`. -/
theorem cross_signs {i j k l : ℝ × ℝ} (hk : sda i j k * sda i j l < 0)
    (hc : sda k l i * sda k l j < 0) :
    0 < -(sda k l i * sda i j k) ∧ 0 < sda k l i * sda i j l := by
  have hD : sda i j k ≠ 0 := by rintro h0; simp [h0] at hk
  rcases hD.lt_or_gt with hneg | hpos
  · have hk' : sda j i k * sda j i l < 0 := by
      have : sda j i k * sda j i l = sda i j k * sda i j l := by unfold sda; ring
      linarith
    have hc' : sda k l j * sda k l i < 0 := by linarith [mul_comm (sda k l i) (sda k l j)]
    have hpos' : 0 < sda j i k := by rw [sda_swap_left]; linarith
    obtain ⟨h1, -, h3⟩ := quad_signs hk' hc' hpos'
    rw [sda_swap_left] at h1
    exact ⟨by nlinarith, by nlinarith⟩
  · obtain ⟨h1, h2, -⟩ := quad_signs hk hc hpos
    exact ⟨by nlinarith, by nlinarith⟩

/-- **The crossing point of the diagonals**: `y = (1 − σ) i + σ j` with `0 < σ < 1` on the line
`kl`. -/
theorem exists_crossing {i j k l : ℝ × ℝ} (hc : sda k l i * sda k l j < 0) :
    ∃ σ : ℝ, 0 < σ ∧ σ < 1 ∧ sda k l ((1 - σ) • i + σ • j) = 0 := by
  set A := sda k l i with hA
  set B := sda k l j with hB
  have hAB : A - B ≠ 0 := by
    intro h
    have : A = B := by linarith
    rw [this] at hc
    nlinarith [mul_self_nonneg B]
  refine ⟨A / (A - B), ?_, ?_, ?_⟩
  · rcases lt_trichotomy A 0 with h | h | h
    · have : 0 < B := by nlinarith
      exact div_pos_of_neg_of_neg h (by linarith)
    · rw [h] at hc; simp at hc
    · have : B < 0 := by nlinarith
      exact div_pos h (by linarith)
  · rcases lt_trichotomy A 0 with h | h | h
    · have : 0 < B := by nlinarith
      exact (div_lt_one_of_neg (by linarith)).2 (by linarith)
    · rw [h] at hc; simp at hc
    · have : B < 0 := by nlinarith
      exact (div_lt_one (by linarith)).2 (by linarith)
  · rw [sda_combo₃ _ _ _ _ (by ring : (1 - A / (A - B)) + A / (A - B) = 1), ← hA, ← hB]
    field_simp
    ring

/-! ## The lifted surface goes down -/

/-- **The lifted surface of `T'` is nowhere above that of `T`**: every triangle `t ∈ T` containing
`x` can be matched with a triangle `t' ∈ T'` containing `x` whose lifted plane at `x` is not
higher. -/
def LiftBelow (T T' : Finset (Finset (ℝ × ℝ))) : Prop :=
  ∀ x a b c, {a, b, c} ∈ T → x ∈ convexHull ℝ {a, b, c} →
    ∃ a' b' c', {a', b', c'} ∈ T' ∧ x ∈ convexHull ℝ {a', b', c'} ∧
      liftPlane a' b' c' x ≤ liftPlane a b c x

theorem LiftBelow.refl (T : Finset (Finset (ℝ × ℝ))) : LiftBelow T T :=
  fun _ a b c h hx => ⟨a, b, c, h, hx, le_rfl⟩

theorem LiftBelow.trans {T T' T'' : Finset (Finset (ℝ × ℝ))} (h1 : LiftBelow T T')
    (h2 : LiftBelow T' T'') : LiftBelow T T'' := by
  intro x a b c h hx
  obtain ⟨a', b', c', h', hx', hle⟩ := h1 x a b c h hx
  obtain ⟨a'', b'', c'', h'', hx'', hle'⟩ := h2 x a' b' c' h' hx'
  exact ⟨a'', b'', c'', h'', hx'', hle'.trans hle⟩

theorem convexHull_mem_of_eq {a b c a' b' c' x : ℝ × ℝ}
    (he : ({a, b, c} : Finset (ℝ × ℝ)) = {a', b', c'}) (hx : x ∈ convexHull ℝ {a, b, c}) :
    x ∈ convexHull ℝ {a', b', c'} := by
  rw [← hull_triple, ← he, hull_triple]; exact hx

/-- **A Delaunay flip lowers the lifted surface.**

[formula: planar-flip-count] -/
theorem liftBelow_flip {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {i j k l : ℝ × ℝ} (h1 : {i, j, k} ∈ T) (h2 : {i, j, l} ∈ T)
    (hkl : k ≠ l) (hin : InCircumcircle i j k l) : LiftBelow T (flip T i j k l) := by
  have hk := hT.opposite h1 h2 hkl
  have hc := crossing_of_inCircumcircle hk hin
  have hin' := inCircumcircle_opposite hk hin
  have hD1 := hT.sda_ne_zero h1
  have hD2 := hT.sda_ne_zero h2
  have hni : sda i k l ≠ 0 := by
    rw [show sda i k l = sda k l i by unfold sda; ring]; rintro h0; simp [h0] at hc
  have hnj : sda j k l ≠ 0 := by
    rw [show sda j k l = sda k l j by unfold sda; ring]; rintro h0; simp [h0] at hc
  have hU := flip_hull_union hk hc
  have m1 : ({i, k, l} : Finset (ℝ × ℝ)) ∈ flip T i j k l := mem_flip.2 (Or.inl rfl)
  have m2 : ({j, k, l} : Finset (ℝ × ℝ)) ∈ flip T i j k l := mem_flip.2 (Or.inr (Or.inl rfl))
  have hD1' : sda i k j ≠ 0 := by
    rw [show sda i k j = -sda i j k by unfold sda; ring]; exact neg_ne_zero.2 hD1
  have hD1r : sda j k i ≠ 0 := by rw [sda_rotate]; exact hD1
  have hD2' : sda i l j ≠ 0 := by
    rw [show sda i l j = -sda i j l by unfold sda; ring]; exact neg_ne_zero.2 hD2
  have hD2r : sda j l i ≠ 0 := by rw [sda_rotate]; exact hD2
  have hni' : sda i l k ≠ 0 := by
    rw [show sda i l k = -sda i k l by unfold sda; ring]; exact neg_ne_zero.2 hni
  have hnj' : sda j l k ≠ 0 := by
    rw [show sda j l k = -sda j k l by unfold sda; ring]; exact neg_ne_zero.2 hnj
  -- the four comparisons of a new and an old plane
  have c11 : ∀ x, x ∈ convexHull ℝ {i, k, l} → liftPlane i k l x ≤ liftPlane i j k x := by
    intro x hx
    have := liftPlane_le_of_mem hni hD1' (inCircumcircle_swap₂₃.2 hin) hx
    rwa [liftPlane_eq_of_eq hD1 (by ext; simp; tauto : ({i, j, k} : Finset (ℝ × ℝ)) = {i, k, j})]
      at this
  have c21 : ∀ x, x ∈ convexHull ℝ {j, k, l} → liftPlane j k l x ≤ liftPlane i j k x := by
    intro x hx
    have := liftPlane_le_of_mem hnj hD1r (inCircumcircle_rotate.2 hin) hx
    rwa [liftPlane_eq_of_eq hD1 (by ext; simp; tauto : ({i, j, k} : Finset (ℝ × ℝ)) = {j, k, i})]
      at this
  have c12 : ∀ x, x ∈ convexHull ℝ {i, k, l} → liftPlane i k l x ≤ liftPlane i j l x := by
    intro x hx
    have hx' : x ∈ convexHull ℝ {i, l, k} := by rw [hull_swap₂₃]; exact hx
    have := liftPlane_le_of_mem hni' hD2' (inCircumcircle_swap₂₃.2 hin') hx'
    rwa [liftPlane_eq_of_eq hD2 (by ext; simp; tauto : ({i, j, l} : Finset (ℝ × ℝ)) = {i, l, j}),
      liftPlane_eq_of_eq hni (by ext; simp; tauto : ({i, k, l} : Finset (ℝ × ℝ)) = {i, l, k})]
      at this
  have c22 : ∀ x, x ∈ convexHull ℝ {j, k, l} → liftPlane j k l x ≤ liftPlane i j l x := by
    intro x hx
    have hx' : x ∈ convexHull ℝ {j, l, k} := by rw [hull_swap₂₃]; exact hx
    have := liftPlane_le_of_mem hnj' hD2r (inCircumcircle_rotate.2 hin') hx'
    rwa [liftPlane_eq_of_eq hD2 (by ext; simp; tauto : ({i, j, l} : Finset (ℝ × ℝ)) = {j, l, i}),
      liftPlane_eq_of_eq hnj (by ext; simp; tauto : ({j, k, l} : Finset (ℝ × ℝ)) = {j, l, k})]
      at this
  intro x a b c habc hx
  by_cases hm : ({a, b, c} : Finset (ℝ × ℝ)) ∈ flip T i j k l
  · exact ⟨a, b, c, hm, hx, le_rfl⟩
  have hD := hT.sda_ne_zero habc
  have hold : ({a, b, c} : Finset (ℝ × ℝ)) = {i, j, k} ∨
      ({a, b, c} : Finset (ℝ × ℝ)) = {i, j, l} := by
    by_contra hno
    push Not at hno
    exact hm (mem_flip.2 (Or.inr (Or.inr ⟨hno.2, hno.1, habc⟩)))
  rcases hold with he | he
  · have hx' := convexHull_mem_of_eq he hx
    rw [← liftPlane_eq_of_eq hD he]
    have : x ∈ convexHull ℝ {i, k, l} ∪ convexHull ℝ {j, k, l} := by
      rw [hU]; exact Or.inl hx'
    rcases this with hx1 | hx2
    · exact ⟨i, k, l, m1, hx1, c11 x hx1⟩
    · exact ⟨j, k, l, m2, hx2, c21 x hx2⟩
  · have hx' := convexHull_mem_of_eq he hx
    rw [← liftPlane_eq_of_eq hD he]
    have : x ∈ convexHull ℝ {i, k, l} ∪ convexHull ℝ {j, k, l} := by
      rw [hU]; exact Or.inr hx'
    rcases this with hx1 | hx2
    · exact ⟨i, k, l, m1, hx1, c12 x hx1⟩
    · exact ⟨j, k, l, m2, hx2, c22 x hx2⟩

theorem liftBelow_delaunayFlip {S : Finset (ℝ × ℝ)} {T T' : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (h : DelaunayFlip T T') : LiftBelow T T' := by
  obtain ⟨i, j, k, l, h1, h2, hkl, hin, rfl⟩ := h
  exact liftBelow_flip hT h1 h2 hkl hin

/-- Along any sequence of Delaunay flips the lifted surface goes down. -/
theorem liftBelow_of_reach {S : Finset (ℝ × ℝ)} {T T' : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (h : ReflTransGen DelaunayFlip T T') :
    IsTriangulation S T' ∧ LiftBelow T T' := by
  induction h with
  | refl => exact ⟨hT, LiftBelow.refl _⟩
  | tail _ hstep ih =>
    obtain ⟨hU, hB⟩ := ih
    exact ⟨isTriangulation_of_delaunayFlip hU hstep, hB.trans (liftBelow_delaunayFlip hU hstep)⟩

/-! ## While `ij` is an edge, the lifted surface contains the lifted segment -/

/-- **On an edge the lifted surface is the lifted segment.** If `{i, j, c}` is a triangle of a
triangulation of `S` and `y = (1 − σ) i + σ j` (`0 < σ < 1`) is not a point of `S`, then every
triangle of the triangulation that contains `y` has the value `(1 − σ)|i|² + σ|j|²` at `y`. -/
theorem liftPlane_eq_of_mem_edge {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {i j c : ℝ × ℝ} (hu : {i, j, c} ∈ T) {σ : ℝ} (hσ0 : 0 < σ)
    (hσ1 : σ < 1) (hyS : (1 - σ) • i + σ • j ∉ S) {a b d : ℝ × ℝ} (h : {a, b, d} ∈ T)
    (hy : (1 - σ) • i + σ • j ∈ convexHull ℝ {a, b, d}) :
    liftPlane a b d ((1 - σ) • i + σ • j) = (1 - σ) * nsq i + σ * nsq j := by
  have hD := hT.sda_ne_zero hu
  have h1 : (1 - σ) + σ = 1 := by ring
  by_cases he : ({a, b, d} : Finset (ℝ × ℝ)) = {i, j, c}
  · rw [liftPlane_eq_of_eq hD he.symm, liftPlane_segment hD]
  · have hy1 : 0 < sda i j c * sda ((1 - σ) • i + σ • j) j c := by
      rw [sda_combo₁ _ _ _ _ h1, sda_self₁₂]
      have := mul_pos (sub_pos.2 hσ1) (mul_self_pos.2 hD)
      nlinarith
    have hy2 : 0 < sda i j c * sda i ((1 - σ) • i + σ • j) c := by
      rw [sda_combo₂ _ _ _ _ h1, sda_self₁₂]
      have := mul_pos hσ0 (mul_self_pos.2 hD)
      nlinarith
    have hy3 : sda i j ((1 - σ) • i + σ • j) = 0 := by
      rw [sda_combo₃ _ _ _ _ h1, sda_self₁₃, sda_self₂₃]; ring
    obtain ⟨d', hd', hs⟩ := eq_insert_of_mem_edge hT hu hy1 hy2 hy3 hyS h he
      (by rw [hull_triple]; exact hy)
    have hD' : sda i j d' ≠ 0 := by rintro h0; simp [h0] at hs
    rw [liftPlane_eq_of_eq hD' hd'.symm, liftPlane_segment hD']

/-! ## A flipped edge never reappears -/

/-- After the Delaunay flip of `ij`, no triangulation below the new lifted surface has `ij` as an
edge. -/
theorem not_mem_of_flipped {S : Finset (ℝ × ℝ)} {T T₂ : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {i j k l : ℝ × ℝ} (h1 : {i, j, k} ∈ T) (h2 : {i, j, l} ∈ T)
    (hkl : k ≠ l) (hin : InCircumcircle i j k l) (hT₂ : IsTriangulation S T₂)
    (hB : LiftBelow (flip T i j k l) T₂) (c : ℝ × ℝ) : {i, j, c} ∉ T₂ := by
  intro hc2
  have hk := hT.opposite h1 h2 hkl
  have hc := crossing_of_inCircumcircle hk hin
  obtain ⟨σ, hσ0, hσ1, hy⟩ := exists_crossing hc
  have h1' : (1 - σ) + σ = 1 := by ring
  set y := (1 - σ) • i + σ • j with hydef
  have hD1 := hT.sda_ne_zero h1
  have hA : sda k l i ≠ 0 := by rintro h0; simp [h0] at hc
  have hB' : sda k l j ≠ 0 := by rintro h0; simp [h0] at hc
  have r1 : sda i k l = sda k l i := by unfold sda; ring
  have r2 : sda i k j = -sda i j k := by unfold sda; ring
  have hn1 : sda i k l ≠ 0 := by rw [r1]; exact hA
  obtain ⟨s1, s2⟩ := cross_signs hk hc
  -- coordinates of the crossing point
  have eik : sda i k y = σ * sda i k j := by
    rw [hydef, sda_combo₃ _ _ _ _ h1', sda_self₁₃]; ring
  have eil : sda i y l = σ * sda i j l := by
    rw [hydef, sda_combo₂ _ _ _ _ h1', sda_self₁₂]; ring
  have ekl : sda y k l = 0 := by rw [← sda_rotate]; exact hy
  have eyj : sda y j k = (1 - σ) * sda i j k := by
    rw [hydef, sda_combo₁ _ _ _ _ h1', sda_self₁₂]; ring
  have eiy : sda i y k = σ * sda i j k := by
    rw [hydef, sda_combo₂ _ _ _ _ h1', sda_self₁₂]; ring
  have eij : sda i j y = 0 := by
    rw [hydef, sda_combo₃ _ _ _ _ h1', sda_self₁₃, sda_self₂₃]; ring
  have hstrict : 0 < sda i k l * sda i k y := by
    rw [eik, r1, r2]
    have := mul_pos hσ0 s1
    linarith
  have hy_mem : y ∈ convexHull ℝ {i, k, l} := by
    refine (mem_hull3_iff hn1).2 ⟨by rw [ekl, mul_zero], ?_, hstrict.le⟩
    rw [eil, r1]
    have := mul_pos hσ0 s2
    nlinarith
  -- the crossing point is not a point of `S`
  have hyS : y ∉ S := by
    intro hyS
    have hyh : y ∈ hull {i, j, k} := by
      rw [hull_triple]
      refine (mem_hull3_iff hD1).2 ⟨?_, ?_, by rw [eij, mul_zero]⟩
      · rw [eyj]
        have := mul_pos (sub_pos.2 hσ1) (mul_self_pos.2 hD1)
        nlinarith
      · rw [eiy]
        have := mul_pos hσ0 (mul_self_pos.2 hD1)
        nlinarith
    have := hT.vertex _ h1 y hyS hyh
    simp only [Finset.mem_insert, Finset.mem_singleton] at this
    rcases this with h | h | h
    · rw [h] at hy; exact hA hy
    · rw [h] at hy; exact hB' hy
    · rw [h] at eij; exact hD1 eij
  -- just after the flip the surface at `y` is strictly below the lifted segment
  have hlt : liftPlane i k l y < (1 - σ) * nsq i + σ * nsq j := by
    have hD1' : sda i k j ≠ 0 := by rw [r2]; exact neg_ne_zero.2 hD1
    have := liftPlane_lt_of_pos hn1 hD1' (inCircumcircle_swap₂₃.2 hin) hstrict
    rwa [liftPlane_eq_of_eq hD1 (by ext; simp; tauto : ({i, j, k} : Finset (ℝ × ℝ)) = {i, k, j}),
      hydef, liftPlane_segment hD1] at this
  obtain ⟨a', b', c', hm, hym, hle⟩ := hB y i k l (mem_flip.2 (Or.inl rfl)) hy_mem
  have hval := liftPlane_eq_of_mem_edge hT₂ hc2 hσ0 hσ1 hyS hm hym
  linarith

/-- **An edge flipped away by a Delaunay flip never reappears** in any later triangulation of a
sequence of Delaunay flips.

[formula: planar-flip-count] -/
theorem flipped_edge_not_mem {S : Finset (ℝ × ℝ)} {T T₂ : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {i j k l : ℝ × ℝ} (h1 : {i, j, k} ∈ T) (h2 : {i, j, l} ∈ T)
    (hkl : k ≠ l) (hin : InCircumcircle i j k l)
    (h : ReflTransGen DelaunayFlip (flip T i j k l) T₂) (c : ℝ × ℝ) : {i, j, c} ∉ T₂ := by
  obtain ⟨hT₂, hB⟩ := liftBelow_of_reach
    (isTriangulation_of_delaunayFlip hT ⟨i, j, k, l, h1, h2, hkl, hin, rfl⟩) h
  exact not_mem_of_flipped hT h1 h2 hkl hin hT₂ hB c

/-! ## The count -/

/-- **At most `C(n, 2)` Delaunay flips.** A sequence of Delaunay flips starting from a
triangulation of `S` has at most `|S|(|S| − 1)/2` flips: the flipped edges are distinct
(`flipped_edge_not_mem`) two-element subsets of `S`.

[formula: planar-flip-count] -/
theorem flip_count_le {S : Finset (ℝ × ℝ)} {f : ℕ → Finset (Finset (ℝ × ℝ))} {n : ℕ}
    (h0 : IsTriangulation S (f 0)) (hstep : ∀ m < n, DelaunayFlip (f m) (f (m + 1))) :
    n ≤ S.card.choose 2 := by
  classical
  have htri : ∀ m ≤ n, IsTriangulation S (f m) := by
    intro m
    induction m with
    | zero => intro _; exact h0
    | succ m ih =>
      intro hm
      exact isTriangulation_of_delaunayFlip (ih (by omega)) (hstep m (by omega))
  have hreach : ∀ m m', m ≤ m' → m' ≤ n → ReflTransGen DelaunayFlip (f m) (f m') := by
    intro m m' hle
    induction m', hle using Nat.le_induction with
    | base => intro _; exact ReflTransGen.refl
    | succ m' _ ih => intro hm'; exact (ih (by omega)).tail (hstep m' (by omega))
  have hex : ∀ m : Fin n, ∃ e ∈ S.powersetCard 2, ∃ i j k l : ℝ × ℝ, e = {i, j} ∧
      {i, j, k} ∈ f m ∧ {i, j, l} ∈ f m ∧ k ≠ l ∧ InCircumcircle i j k l ∧
      f ((m : ℕ) + 1) = flip (f m) i j k l := by
    intro m
    obtain ⟨i, j, k, l, h1, h2, hkl, hin, hf⟩ := hstep m m.2
    have hT := htri m m.2.le
    have hij := sda_ne_zero_left (hT.sda_ne_zero h1)
    refine ⟨{i, j}, Finset.mem_powersetCard.2 ⟨?_, Finset.card_pair hij⟩, i, j, k, l, rfl, h1, h2,
      hkl, hin, hf⟩
    intro v hv
    simp only [Finset.mem_insert, Finset.mem_singleton] at hv
    rcases hv with rfl | rfl
    · exact hT.mem_S h1 (by simp)
    · exact hT.mem_S h1 (by simp)
  choose e he hdata using hex
  have key : ∀ m m' : Fin n, (m : ℕ) < m' → e m ≠ e m' := by
    intro m m' hlt heq
    obtain ⟨i, j, k, l, hei, h1, h2, hkl, hin, hf⟩ := hdata m
    obtain ⟨i', j', k', l', hei', h1', -⟩ := hdata m'
    have hr : ReflTransGen DelaunayFlip (flip (f m) i j k l) (f m') := by
      rw [← hf]; exact hreach (m + 1) m' (by omega) m'.2.le
    refine flipped_edge_not_mem (htri m m.2.le) h1 h2 hkl hin hr k' ?_
    have hpair : ({i', j'} : Finset (ℝ × ℝ)) = {i, j} := by rw [← hei', ← heq, hei]
    have hset : ({i, j, k'} : Finset (ℝ × ℝ)) = {i', j', k'} := by
      ext v
      have := congrArg (v ∈ ·) hpair
      simp only [Finset.mem_insert, Finset.mem_singleton, eq_iff_iff] at this ⊢
      tauto
    rw [hset]; exact h1'
  have hinj : Function.Injective e := by
    intro m m' heq
    rcases lt_trichotomy (m : ℕ) m' with h | h | h
    · exact absurd heq (key m m' h)
    · exact Fin.ext h
    · exact absurd heq.symm (key m' m h)
  have := Finset.card_le_card_of_injOn e (s := Finset.univ) (t := S.powersetCard 2)
    (fun m _ => he m) hinj.injOn
  simpa [Finset.card_powersetCard] using this

end Planar

end Delaunay

end DDG
