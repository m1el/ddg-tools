import DdgProofs.Delaunay.PlanarLemma

/-!
# Planar Delaunay triangulations, part 5: uniqueness up to cocircular points

**Uniqueness.** If no four points of the finite point set `S` are cocircular, then `S` has at
most one Delaunay triangulation (`delaunay_unique`); with existence
(`DdgProofs/Delaunay/PlanarLemma.lean`), exactly one as soon as `S` has some triangulation
(`existsUnique_delaunay`), and every sequence of Delaunay flips that cannot be continued ends at
it, whatever the order of the flips (`eq_of_isLocallyDelaunay`). Without the genericity
hypothesis: **two Delaunay triangulations of `S` differ only in triangles whose circumcircle
passes through a fourth point of `S`** (`exists_cocircular_of_not_mem`).

## Genericity

`NoFourCocircular S`: for every nondegenerate triangle `a b c` with corners in `S` and every other
point `p ∈ S`, `inCircle a b c p ≠ 0` — `p` is not on the circumcircle of `a b c`
(`inCircumcircle_iff_lt`, `power_eq_of_equidistant`: `inCircle a b c p = 0` iff `p` is on the
circle). This is the usual "no four points on a circle" (de Berg et al., Sec. 9.2), required only
for triples that span a triangle: collinear triples are not excluded.

## The proof: the lifted planes

The **lifted plane** of a triangle `a b c` is the affine function `liftPlane a b c` with the values
`|v|²` at the three corners (`liftPlane_self₁`, …); it does not depend on the order of the corners
(`liftPlane_eq_of_eq`), and `circleLift a b c x = liftPlane a b c x − |x|²`
(`circleLift_eq_liftPlane_sub`): `x` is inside the circumcircle iff the lifted plane is above the
paraboloid at `x`. Affine functions satisfy `sda p q r · g(x) = sda x q r · g(p) + sda p x r · g(q)
+ sda p q x · g(r)` (`liftPlane_bary`: Cramer's rule).

Let `t ∈ T` and `t' ∈ T'` be triangles of two Delaunay triangulations, `x` an interior point of
`t` that lies in `t'` (`T'` covers the hull). The affine function `M = liftPlane t − liftPlane t'`
is `≥ 0` at the corners of `t` (they are not strictly inside the circumcircle of `t'`) and `≤ 0`
at the corners of `t'`. Writing `M(x)` in the barycentric coordinates of `t` (all positive) and
of `t'` (all nonnegative) gives `M(x) ≥ 0 ≥ M(x)`; so `M` vanishes at the three corners of `t`,
hence everywhere, and every corner of `t'` lies on the circumcircle of `t`
(`inCircle_eq_zero_of_mem`). If `t ≠ t'`, one of these corners is not a corner of `t`. This is the
lifting proof (Edelsbrunner, Sec. 1.3; de Berg et al., Thm. 9.6–9.7: the Delaunay triangulation is
the projection of the lower convex hull of the lifted points, unique when no four points are
cocircular), with the hull replaced by the comparison of two lifted planes at one point.

## References

* B. Delaunay, *Sur la sphère vide*, Bull. Acad. Sci. URSS (1934), 793–800.
* M. de Berg, O. Cheong, M. van Kreveld, M. Overmars, *Computational Geometry: Algorithms and
  Applications*, 3rd ed., Springer, 2008, Sec. 9.2, Theorems 9.6 and 9.7.
* H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge University Press, 2001,
  Sec. 1.3 (the lifting map; uniqueness in general position).

The bibliographic data were not checked against the works here.
-/

noncomputable section

namespace DdgProofs

namespace Delaunay

namespace Planar

open Set

/-! ## The lifted plane of a triangle -/

/-- **The lifted plane** of the triangle `a b c`: the affine function of `x` with the values `|v|²`
at the corners, `Σ_v λ_v(x) |v|²` with the barycentric coordinates `λ`. -/
def liftPlane (a b c x : ℝ × ℝ) : ℝ :=
  (nsq a * sda x b c + nsq b * sda a x c + nsq c * sda a b x) / sda a b c

/-- `circleLift` (positive exactly inside the circumcircle) is the height of the lifted plane above
the paraboloid. -/
theorem circleLift_eq_liftPlane_sub {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) (x : ℝ × ℝ) :
    circleLift a b c x = liftPlane a b c x - nsq x := by
  unfold circleLift liftPlane
  rw [inCircle_eq_lift]
  field_simp

/-- **The lifted plane is affine** (Cramer's rule): for every triangle `p q r`,
`sda p q r · L(x) = sda x q r · L(p) + sda p x r · L(q) + sda p q x · L(r)`. -/
theorem liftPlane_bary (a b c p q r x : ℝ × ℝ) :
    sda p q r * liftPlane a b c x = sda x q r * liftPlane a b c p + sda p x r * liftPlane a b c q +
      sda p q x * liftPlane a b c r := by
  have key : sda p q r * (nsq a * sda x b c + nsq b * sda a x c + nsq c * sda a b x) =
      sda x q r * (nsq a * sda p b c + nsq b * sda a p c + nsq c * sda a b p) +
      sda p x r * (nsq a * sda q b c + nsq b * sda a q c + nsq c * sda a b q) +
      sda p q x * (nsq a * sda r b c + nsq b * sda a r c + nsq c * sda a b r) := by
    unfold sda; ring
  unfold liftPlane
  simp only [div_eq_mul_inv]
  linear_combination (sda a b c)⁻¹ * key

theorem liftPlane_self₁ {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : liftPlane a b c a = nsq a := by
  unfold liftPlane; rw [sda_self₁₂, sda_self₁₃]; field_simp; ring

theorem liftPlane_self₂ {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : liftPlane a b c b = nsq b := by
  unfold liftPlane; rw [sda_self₁₂, sda_self₂₃]; field_simp; ring

theorem liftPlane_self₃ {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : liftPlane a b c c = nsq c := by
  unfold liftPlane; rw [sda_self₁₃, sda_self₂₃]; field_simp; ring

/-- Two affine functions, given as lifted planes, that agree at the corners of a nondegenerate
triangle agree everywhere. -/
theorem liftPlane_eq_of_agree {a b c a' b' c' p q r : ℝ × ℝ} (h : sda p q r ≠ 0)
    (hp : liftPlane a b c p = liftPlane a' b' c' p) (hq : liftPlane a b c q = liftPlane a' b' c' q)
    (hr : liftPlane a b c r = liftPlane a' b' c' r) (x : ℝ × ℝ) :
    liftPlane a b c x = liftPlane a' b' c' x := by
  have e1 := liftPlane_bary a b c p q r x
  have e2 := liftPlane_bary a' b' c' p q r x
  rw [hp, hq, hr] at e1
  exact mul_left_cancel₀ h (e1.trans e2.symm)

/-- **The lifted plane depends only on the set of corners.** -/
theorem liftPlane_eq_of_eq {a b c a' b' c' : ℝ × ℝ} (h : sda a b c ≠ 0)
    (he : ({a, b, c} : Finset (ℝ × ℝ)) = {a', b', c'}) (x : ℝ × ℝ) :
    liftPlane a' b' c' x = liftPlane a b c x := by
  have h' : sda a' b' c' ≠ 0 := IsTri.sda_ne_zero ⟨a, b, c, h, rfl⟩ he
  have val : ∀ v ∈ ({a, b, c} : Finset (ℝ × ℝ)), liftPlane a' b' c' v = liftPlane a b c v := by
    intro v hv
    have hv' := hv
    rw [he] at hv'
    have e : liftPlane a' b' c' v = nsq v := by
      simp only [Finset.mem_insert, Finset.mem_singleton] at hv'
      rcases hv' with rfl | rfl | rfl
      · exact liftPlane_self₁ h'
      · exact liftPlane_self₂ h'
      · exact liftPlane_self₃ h'
    rw [e]
    simp only [Finset.mem_insert, Finset.mem_singleton] at hv
    rcases hv with rfl | rfl | rfl
    · exact (liftPlane_self₁ h).symm
    · exact (liftPlane_self₂ h).symm
    · exact (liftPlane_self₃ h).symm
  exact liftPlane_eq_of_agree h (val a (by simp)) (val b (by simp)) (val c (by simp)) x

/-- Not strictly inside the circumcircle: the lifted plane is not above the paraboloid. -/
theorem liftPlane_le_of_not_inCircumcircle {a b c p : ℝ × ℝ} (h : sda a b c ≠ 0)
    (hp : ¬ InCircumcircle a b c p) : liftPlane a b c p ≤ nsq p := by
  rw [← circleLift_pos_iff h, not_lt, circleLift_eq_liftPlane_sub h] at hp
  linarith

/-! ## Two Delaunay triangles that overlap have the same circumcircle -/

/-- **Two empty-circle triangles that overlap are cocircular.** Let `a b c` and `a' b' c'` be
nondegenerate triangles whose corners are not strictly inside the other's circumcircle, and let an
interior point `x` of `a b c` lie in `a' b' c'`. Then the corners of `a' b' c'` lie on the
circumcircle of `a b c`. -/
theorem inCircle_eq_zero_of_mem {a b c a' b' c' x : ℝ × ℝ} (hD : sda a b c ≠ 0)
    (hD' : sda a' b' c' ≠ 0)
    (ha : ¬ InCircumcircle a' b' c' a) (hb : ¬ InCircumcircle a' b' c' b)
    (hc : ¬ InCircumcircle a' b' c' c)
    (ha' : ¬ InCircumcircle a b c a') (hb' : ¬ InCircumcircle a b c b')
    (hc' : ¬ InCircumcircle a b c c')
    (hx : x ∈ interior (convexHull ℝ {a, b, c})) (hx' : x ∈ convexHull ℝ {a', b', c'}) :
    inCircle a b c a' = 0 ∧ inCircle a b c b' = 0 ∧ inCircle a b c c' = 0 := by
  -- `M = liftPlane (a b c) − liftPlane (a' b' c')`
  set M : ℝ × ℝ → ℝ := fun y => liftPlane a b c y - liftPlane a' b' c' y with hM
  have bary : ∀ p q r y : ℝ × ℝ, sda p q r * M y = sda y q r * M p + sda p y r * M q +
      sda p q y * M r := by
    intro p q r y
    simp only [hM]
    linear_combination liftPlane_bary a b c p q r y - liftPlane_bary a' b' c' p q r y
  have Ma : 0 ≤ M a := by
    simp only [hM, liftPlane_self₁ hD]; linarith [liftPlane_le_of_not_inCircumcircle hD' ha]
  have Mb : 0 ≤ M b := by
    simp only [hM, liftPlane_self₂ hD]; linarith [liftPlane_le_of_not_inCircumcircle hD' hb]
  have Mc : 0 ≤ M c := by
    simp only [hM, liftPlane_self₃ hD]; linarith [liftPlane_le_of_not_inCircumcircle hD' hc]
  have Ma' : M a' ≤ 0 := by
    simp only [hM, liftPlane_self₁ hD']; linarith [liftPlane_le_of_not_inCircumcircle hD ha']
  have Mb' : M b' ≤ 0 := by
    simp only [hM, liftPlane_self₂ hD']; linarith [liftPlane_le_of_not_inCircumcircle hD hb']
  have Mc' : M c' ≤ 0 := by
    simp only [hM, liftPlane_self₃ hD']; linarith [liftPlane_le_of_not_inCircumcircle hD hc']
  obtain ⟨p1, p2, p3⟩ := pos_of_mem_interior_hull3 hD hx
  obtain ⟨q1, q2, q3⟩ := (mem_hull3_iff hD').1 hx'
  have hDD : 0 < sda a b c * sda a b c := mul_self_pos.2 hD
  have hDD' : 0 < sda a' b' c' * sda a' b' c' := mul_self_pos.2 hD'
  -- `M x ≥ 0` from the side of `a b c`, `M x ≤ 0` from the side of `a' b' c'`
  have e1 : (sda a b c * sda a b c) * M x = (sda a b c * sda x b c) * M a +
      (sda a b c * sda a x c) * M b + (sda a b c * sda a b x) * M c := by
    rw [mul_assoc, bary a b c x]; ring
  have e2 : (sda a' b' c' * sda a' b' c') * M x = (sda a' b' c' * sda x b' c') * M a' +
      (sda a' b' c' * sda a' x c') * M b' + (sda a' b' c' * sda a' b' x) * M c' := by
    rw [mul_assoc, bary a' b' c' x]; ring
  have t1 := mul_nonneg p1.le Ma
  have t2 := mul_nonneg p2.le Mb
  have t3 := mul_nonneg p3.le Mc
  have s1 := mul_nonpos_of_nonneg_of_nonpos q1 Ma'
  have s2 := mul_nonpos_of_nonneg_of_nonpos q2 Mb'
  have s3 := mul_nonpos_of_nonneg_of_nonpos q3 Mc'
  have hMx1 : 0 ≤ M x := by
    by_contra hneg; push Not at hneg; nlinarith [mul_neg_of_pos_of_neg hDD hneg]
  have hMx2 : M x ≤ 0 := by
    by_contra hpos; push Not at hpos; nlinarith [mul_pos hDD' hpos]
  have hMx : M x = 0 := le_antisymm hMx2 hMx1
  rw [hMx, mul_zero] at e1
  -- so `M` vanishes at `a, b, c`, hence everywhere
  have za : M a = 0 := by
    by_contra hne
    have := mul_pos p1 (lt_of_le_of_ne Ma (Ne.symm hne)); linarith
  have zb : M b = 0 := by
    by_contra hne
    have := mul_pos p2 (lt_of_le_of_ne Mb (Ne.symm hne)); linarith
  have zc : M c = 0 := by
    by_contra hne
    have := mul_pos p3 (lt_of_le_of_ne Mc (Ne.symm hne)); linarith
  have zero : ∀ y, M y = 0 := by
    intro y
    have := bary a b c y
    rw [za, zb, zc] at this
    simpa [hD] using this
  have cl : ∀ v, liftPlane a' b' c' v = nsq v → inCircle a b c v = 0 := by
    intro v hv
    have h0 := zero v
    simp only [hM, hv] at h0
    have : circleLift a b c v = 0 := by rw [circleLift_eq_liftPlane_sub hD]; linarith
    unfold circleLift at this
    simpa [hD] using this
  exact ⟨cl a' (liftPlane_self₁ hD'), cl b' (liftPlane_self₂ hD'), cl c' (liftPlane_self₃ hD')⟩

/-! ## Uniqueness -/

/-- **Two Delaunay triangulations differ only at cocircular points.** If `{a, b, c}` is a triangle
of a Delaunay triangulation `T` of `S` but not of another Delaunay triangulation `T'` of `S`, then
some point of `S` other than `a, b, c` lies on the circumcircle of `a b c`.

[formula: planar-delaunay-uniqueness] -/
theorem exists_cocircular_of_not_mem {S : Finset (ℝ × ℝ)} {T T' : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (hD : IsDelaunay S T) (hT' : IsTriangulation S T')
    (hD' : IsDelaunay S T') {a b c : ℝ × ℝ} (ht : {a, b, c} ∈ T) (hn : {a, b, c} ∉ T') :
    ∃ p ∈ S, p ∉ ({a, b, c} : Finset (ℝ × ℝ)) ∧ inCircle a b c p = 0 := by
  have hDabc := hT.sda_ne_zero ht
  obtain ⟨x, hx⟩ := interior_hull3_nonempty hDabc
  have hxS : x ∈ convexHull ℝ (S : Set (ℝ × ℝ)) :=
    hT.hull_subset ht (by rw [hull_triple]; exact interior_subset hx)
  rw [← hT'.cover] at hxS
  obtain ⟨t', ht', hxt'⟩ := Set.mem_iUnion₂.1 hxS
  obtain ⟨a', b', c', hD'', rfl⟩ := hT'.tri t' ht'
  rw [hull_triple] at hxt'
  have hSa := hT.mem_S ht (x := a) (by simp)
  have hSb := hT.mem_S ht (x := b) (by simp)
  have hSc := hT.mem_S ht (x := c) (by simp)
  have hSa' := hT'.mem_S ht' (x := a') (by simp)
  have hSb' := hT'.mem_S ht' (x := b') (by simp)
  have hSc' := hT'.mem_S ht' (x := c') (by simp)
  obtain ⟨za, zb, zc⟩ := inCircle_eq_zero_of_mem hDabc hD'' (hD' _ _ _ ht' a hSa)
    (hD' _ _ _ ht' b hSb) (hD' _ _ _ ht' c hSc) (hD _ _ _ ht a' hSa') (hD _ _ _ ht b' hSb')
    (hD _ _ _ ht c' hSc') hx hxt'
  -- some corner of `a' b' c'` is not a corner of `a b c`
  by_contra hno
  push Not at hno
  have hsub : ({a', b', c'} : Finset (ℝ × ℝ)) ⊆ {a, b, c} := by
    intro v hv
    simp only [Finset.mem_insert, Finset.mem_singleton] at hv
    rcases hv with rfl | rfl | rfl
    · by_contra hv; exact hno _ hSa' hv za
    · by_contra hv; exact hno _ hSb' hv zb
    · by_contra hv; exact hno _ hSc' hv zc
  have heq : ({a', b', c'} : Finset (ℝ × ℝ)) = {a, b, c} :=
    Finset.eq_of_subset_of_card_le hsub (by
      rw [(hT.tri _ ht).card, (hT'.tri _ ht').card])
  exact hn (heq ▸ ht')

/-- **General position**: no point of `S` lies on the circumcircle of a nondegenerate triangle
spanned by three other points of `S`. -/
def NoFourCocircular (S : Finset (ℝ × ℝ)) : Prop :=
  ∀ a ∈ S, ∀ b ∈ S, ∀ c ∈ S, ∀ p ∈ S, sda a b c ≠ 0 → p ∉ ({a, b, c} : Finset (ℝ × ℝ)) →
    inCircle a b c p ≠ 0

/-- In general position, every triangle of a Delaunay triangulation is a triangle of every other
Delaunay triangulation. -/
theorem subset_of_isDelaunay {S : Finset (ℝ × ℝ)} (hg : NoFourCocircular S)
    {T T' : Finset (Finset (ℝ × ℝ))} (hT : IsTriangulation S T) (hD : IsDelaunay S T)
    (hT' : IsTriangulation S T') (hD' : IsDelaunay S T') : T ⊆ T' := by
  intro t ht
  by_contra hn
  obtain ⟨a, b, c, h, rfl⟩ := hT.tri t ht
  obtain ⟨p, hp, hpn, h0⟩ := exists_cocircular_of_not_mem hT hD hT' hD' ht hn
  exact hg a (hT.mem_S ht (by simp)) b (hT.mem_S ht (by simp)) c (hT.mem_S ht (by simp)) p hp h
    hpn h0

/-- **Uniqueness of the Delaunay triangulation** (Delaunay 1934; de Berg et al., Thm. 9.7): if no
four points of `S` are cocircular, any two Delaunay triangulations of `S` are equal.

[formula: planar-delaunay-uniqueness] -/
theorem delaunay_unique {S : Finset (ℝ × ℝ)} (hg : NoFourCocircular S)
    {T T' : Finset (Finset (ℝ × ℝ))} (hT : IsTriangulation S T) (hD : IsDelaunay S T)
    (hT' : IsTriangulation S T') (hD' : IsDelaunay S T') : T = T' :=
  Finset.Subset.antisymm (subset_of_isDelaunay hg hT hD hT' hD')
    (subset_of_isDelaunay hg hT' hD' hT hD)

/-- In general position, any two locally Delaunay triangulations of `S` are equal: Lawson's
algorithm ends at the same triangulation whatever the order of the flips. -/
theorem eq_of_isLocallyDelaunay {S : Finset (ℝ × ℝ)} (hg : NoFourCocircular S)
    {T T' : Finset (Finset (ℝ × ℝ))} (hT : IsTriangulation S T) (hL : IsLocallyDelaunay T)
    (hT' : IsTriangulation S T') (hL' : IsLocallyDelaunay T') : T = T' :=
  delaunay_unique hg hT (delaunay_lemma hT hL) hT' (delaunay_lemma hT' hL')

/-- **Existence and uniqueness**: a finite point set in general position that has a triangulation
has exactly one Delaunay triangulation. -/
theorem existsUnique_delaunay {S : Finset (ℝ × ℝ)} (hg : NoFourCocircular S)
    (h : ∃ T, IsTriangulation S T) : ∃! T, IsTriangulation S T ∧ IsDelaunay S T := by
  obtain ⟨T, hT, hD⟩ := exists_delaunay_of_triangulation h
  exact ⟨T, ⟨hT, hD⟩, fun T' ⟨hT', hD'⟩ => delaunay_unique hg hT' hD' hT hD⟩

end Planar

end Delaunay

end DdgProofs
