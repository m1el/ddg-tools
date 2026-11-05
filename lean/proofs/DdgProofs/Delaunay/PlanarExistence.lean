import DdgProofs.Delaunay.PlanarInsertion

/-!
# Planar Delaunay triangulations, part 6: existence of a triangulation

**Every finite point set `S ⊂ ℝ²` that is not collinear has a triangulation**
(`exists_triangulation`), and hence, by Lawson's algorithm and Delaunay's lemma
(`exists_delaunay_of_triangulation`, `DdgProofs/Delaunay/PlanarLemma.lean`), a **Delaunay
triangulation** (`exists_delaunay`). Non-collinearity (`NonCollinear S`: three points of `S` with
`sda a b c ≠ 0`) is necessary for a nonempty `S`: the triangles of a triangulation are nondegenerate
and cover `conv S ⊇ S` (`nonCollinear_of_isTriangulation`); `exists_triangulation_iff` records both
directions.

## The route: removing extreme points

This is the incremental construction ("sweep") of de Berg et al., Ch. 9.1 and Edelsbrunner, Sec. 1.1,
run backwards as an induction on `|S|`:

* `|S| = 3`: the three points are a nondegenerate triangle (`isTriangulation_triangle`).
* `|S| ≥ 4`: the lexicographically largest point `e⁺` of `S` (largest first coordinate, then largest
  second coordinate) and the lexicographically smallest `e⁻` are **extreme**: `e ∉ conv(S \ {e})`
  (`exists_extreme`; the set of points lexicographically below `e` is convex). At least one of
  `S \ {e⁺}`, `S \ {e⁻}` is still not collinear (`exists_extreme_nonCollinear`): if `S \ {e⁺}` lies
  on a line `L`, then `e⁺ ∉ L`, and `S \ {e⁻}` contains `e⁺` and at least two points of `L`. By
  induction `S \ {e}` has a triangulation, and inserting `e`, which lies outside its convex hull,
  by coning it to the visible boundary edges (`isTriangulation_insert`,
  `DdgProofs/Delaunay/PlanarInsertion.lean`) triangulates `S`.

The choice of the removed point avoids two cases of the general incremental algorithm: only points
outside the current hull are inserted (no splitting of triangles or edges), and the induction never
meets a collinear set (no fan over collinear points).

## References

* M. de Berg, O. Cheong, M. van Kreveld, M. Overmars, *Computational Geometry: Algorithms and
  Applications*, 3rd ed., Springer, 2008, Ch. 9.1 (every point set has a triangulation; its number
  of triangles) and Ch. 1 (the incremental convex hull in lexicographic order).
* H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge University Press, 2001,
  Sec. 1.1 (triangulations of point sets exist: sort the points and cone each new point to the
  visible part of the hull) and Sec. 1.2 (flipping to the Delaunay triangulation).
* C. L. Lawson, *Software for C¹ surface interpolation*, in: Mathematical Software III, Academic
  Press, 1977, 161–194 (a first triangulation, then flips to the Delaunay one).

The bibliographic data were not checked against the works here; the statements are proved from
the definitions of `DdgProofs/Delaunay/PlanarFlip.lean`, not taken from the works.
-/

noncomputable section

namespace DdgProofs

namespace Delaunay

namespace Planar

open Set

/-- `S` is **not collinear**: it contains three points spanning a nondegenerate triangle. -/
def NonCollinear (S : Finset (ℝ × ℝ)) : Prop :=
  ∃ a ∈ S, ∃ b ∈ S, ∃ c ∈ S, sda a b c ≠ 0

/-! ## Extreme points -/

/-- **A lexicographically extreme point.** For `σ ≠ 0` (in practice `σ = ±1`), the point `e` of `S`
that maximizes `σ e.1`, and among those `σ e.2`, is not in the convex hull of the other points. -/
theorem exists_extreme (S : Finset (ℝ × ℝ)) (hS : S.Nonempty) {σ : ℝ} (hσ : σ ≠ 0) :
    ∃ e ∈ S, (∀ s ∈ S, σ * s.1 ≤ σ * e.1) ∧
      e ∉ convexHull ℝ ((S.erase e : Finset (ℝ × ℝ)) : Set (ℝ × ℝ)) := by
  classical
  obtain ⟨m, hm, hmax⟩ := S.exists_max_image (fun z => σ * z.1) hS
  have hne : (S.filter fun z => σ * z.1 = σ * m.1).Nonempty := ⟨m, Finset.mem_filter.2 ⟨hm, rfl⟩⟩
  obtain ⟨e, he, hemax⟩ := (S.filter fun z => σ * z.1 = σ * m.1).exists_max_image
    (fun z => σ * z.2) hne
  rw [Finset.mem_filter] at he
  refine ⟨e, he.1, fun s hs => he.2 ▸ hmax s hs, ?_⟩
  set M := σ * m.1
  let H : Set (ℝ × ℝ) := {z | σ * z.1 < M ∨ (σ * z.1 = M ∧ σ * z.2 < σ * e.2)}
  have hsub : ((S.erase e : Finset (ℝ × ℝ)) : Set (ℝ × ℝ)) ⊆ H := by
    intro s hs
    rw [Finset.mem_coe, Finset.mem_erase] at hs
    rcases (hmax s hs.2).lt_or_eq with h | h
    · exact Or.inl h
    · refine Or.inr ⟨h, lt_of_le_of_ne (hemax s (Finset.mem_filter.2 ⟨hs.2, h⟩)) fun h2 => hs.1 ?_⟩
      ext
      · exact mul_left_cancel₀ hσ (h.trans he.2.symm)
      · exact mul_left_cancel₀ hσ h2
  have hHc : Convex ℝ H := by
    intro z hz w hw s t hs ht hst
    rcases hs.eq_or_lt with rfl | hs'
    · have : t = 1 := by linarith
      subst this
      simpa using hw
    rcases ht.eq_or_lt with rfl | ht'
    · have : s = 1 := by linarith
      subst this
      simpa using hz
    have hz1 : σ * z.1 ≤ M := by rcases hz with h | h; exact h.le; exact h.1.le
    have hw1 : σ * w.1 ≤ M := by rcases hw with h | h; exact h.le; exact h.1.le
    show σ * (s • z + t • w).1 < M ∨
      (σ * (s • z + t • w).1 = M ∧ σ * (s • z + t • w).2 < σ * e.2)
    simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    have hM : s * M + t * M = M := by rw [← add_mul, hst, one_mul]
    have hE : s * (σ * e.2) + t * (σ * e.2) = σ * e.2 := by rw [← add_mul, hst, one_mul]
    rcases hz with hz | ⟨hz1', hz2⟩
    · left; nlinarith [mul_lt_mul_of_pos_left hz hs', mul_le_mul_of_nonneg_left hw1 ht, hM]
    rcases hw with hw | ⟨hw1', hw2⟩
    · left; nlinarith [mul_lt_mul_of_pos_left hw ht', mul_le_mul_of_nonneg_left hz1 hs, hM]
    right
    constructor
    · have : σ * (s * z.1 + t * w.1) = s * (σ * z.1) + t * (σ * w.1) := by ring
      rw [this, hz1', hw1']; linear_combination M * hst
    · nlinarith [mul_lt_mul_of_pos_left hz2 hs', mul_lt_mul_of_pos_left hw2 ht', hE]
  have heH : e ∉ H := by
    rintro (h | ⟨-, h⟩)
    · rw [he.2] at h; exact lt_irrefl _ h
    · exact lt_irrefl _ h
  exact fun he' => heH (convexHull_min hsub hHc he')

/-! ## Removing an extreme point keeps a non-collinear set non-collinear -/

/-- Points `w ≠ v` on the line `uv` (`sda u v w = 0`): if `p, v, u` is a nondegenerate triangle,
so is `p, v, w`. -/
theorem sda_ne_zero_of_collinear {p u v w : ℝ × ℝ} (h : sda p v u ≠ 0) (hcol : sda u v w = 0)
    (hw : w ≠ v) : sda p v w ≠ 0 := by
  intro h0
  have i1 : sda p v u * (w.1 - v.1) - sda p v w * (u.1 - v.1) = (v.1 - p.1) * sda u v w := by
    unfold sda; ring
  have i2 : sda p v u * (w.2 - v.2) - sda p v w * (u.2 - v.2) = (v.2 - p.2) * sda u v w := by
    unfold sda; ring
  rw [h0, hcol] at i1 i2
  have j1 : sda p v u * (w.1 - v.1) = 0 := by linarith
  have j2 : sda p v u * (w.2 - v.2) = 0 := by linarith
  apply hw
  ext
  · linarith [(mul_eq_zero.1 j1).resolve_left h]
  · linarith [(mul_eq_zero.1 j2).resolve_left h]

/-- **Some extreme point can be removed.** A non-collinear set of at least four points has a point
`e` outside the convex hull of the others such that the others are still not collinear. -/
theorem exists_extreme_nonCollinear {S : Finset (ℝ × ℝ)} (hS : NonCollinear S)
    (hcard : 4 ≤ S.card) :
    ∃ e ∈ S, e ∉ convexHull ℝ ((S.erase e : Finset (ℝ × ℝ)) : Set (ℝ × ℝ)) ∧
      NonCollinear (S.erase e) := by
  classical
  have hSne : S.Nonempty := Finset.card_pos.1 (by omega)
  obtain ⟨p, hpS, hpmax, hpext⟩ := exists_extreme S hSne one_ne_zero
  by_cases h : NonCollinear (S.erase p)
  · exact ⟨p, hpS, hpext, h⟩
  obtain ⟨e, heS, hemax, heext⟩ := exists_extreme S hSne (neg_ne_zero.2 one_ne_zero)
  refine ⟨e, heS, heext, ?_⟩
  -- `e ≠ p`: otherwise all of `S` is on a vertical line
  have hep : e ≠ p := by
    rintro rfl
    obtain ⟨a, ha, b, hb, c, hc, habc⟩ := hS
    have k : ∀ s ∈ S, s.1 = e.1 := fun s hs => by
      have h1 := hpmax s hs
      have h2 := hemax s hs
      linarith
    apply habc
    unfold sda
    rw [k a ha, k b hb, k c hc]
    ring
  unfold NonCollinear at h
  push Not at h
  -- a nondegenerate triangle of `S` has the corner `p`
  obtain ⟨u, hu, v, hv, huv⟩ : ∃ u ∈ S.erase p, ∃ v ∈ S.erase p, sda p u v ≠ 0 := by
    obtain ⟨a, ha, b, hb, c, hc, habc⟩ := hS
    by_cases hap : a = p
    · subst hap
      exact ⟨b, Finset.mem_erase.2 ⟨fun h => habc (by rw [h, sda_self₁₂]), hb⟩,
        c, Finset.mem_erase.2 ⟨fun h => habc (by rw [h, sda_self₁₃]), hc⟩, habc⟩
    by_cases hbp : b = p
    · subst hbp
      refine ⟨c, Finset.mem_erase.2 ⟨fun h => habc (by rw [h, sda_self₂₃]), hc⟩,
        a, Finset.mem_erase.2 ⟨hap, ha⟩, ?_⟩
      rwa [sda_rotate]
    by_cases hcp : c = p
    · subst hcp
      refine ⟨a, Finset.mem_erase.2 ⟨hap, ha⟩, b, Finset.mem_erase.2 ⟨hbp, hb⟩, ?_⟩
      rwa [← sda_rotate]
    exact absurd (h a (Finset.mem_erase.2 ⟨hap, ha⟩) b (Finset.mem_erase.2 ⟨hbp, hb⟩) c
      (Finset.mem_erase.2 ⟨hcp, hc⟩)) habc
  obtain ⟨hup, huS⟩ := Finset.mem_erase.1 hu
  obtain ⟨hvp, hvS⟩ := Finset.mem_erase.1 hv
  have hu_ne_v : u ≠ v := sda_ne_zero_right huv
  have hpe : p ∈ S.erase e := Finset.mem_erase.2 ⟨hep.symm, hpS⟩
  -- a fourth point
  obtain ⟨w, hw⟩ : (((S.erase p).erase u).erase v).Nonempty := by
    apply Finset.card_pos.1
    rw [Finset.card_erase_of_mem (Finset.mem_erase.2 ⟨hu_ne_v.symm, hv⟩),
      Finset.card_erase_of_mem hu, Finset.card_erase_of_mem hpS]
    omega
  obtain ⟨hwv, hw⟩ := Finset.mem_erase.1 hw
  obtain ⟨hwu, hw⟩ := Finset.mem_erase.1 hw
  obtain ⟨hwp, hwS⟩ := Finset.mem_erase.1 hw
  by_cases hue : u = e
  · subst hue
    have h' : sda p v w ≠ 0 := sda_ne_zero_of_collinear
      (by rw [show sda p v u = -sda p u v by unfold sda; ring]; exact neg_ne_zero.2 huv)
      (h u hu v hv w (Finset.mem_erase.2 ⟨hwp, hwS⟩)) hwv
    exact ⟨p, hpe, v, Finset.mem_erase.2 ⟨hu_ne_v.symm, hvS⟩, w,
      Finset.mem_erase.2 ⟨hwu, hwS⟩, h'⟩
  by_cases hve : v = e
  · subst hve
    have h' : sda p u w ≠ 0 := sda_ne_zero_of_collinear huv
      (h v hv u hu w (Finset.mem_erase.2 ⟨hwp, hwS⟩)) hwu
    exact ⟨p, hpe, u, Finset.mem_erase.2 ⟨hu_ne_v, huS⟩, w,
      Finset.mem_erase.2 ⟨hwv, hwS⟩, h'⟩
  exact ⟨p, hpe, u, Finset.mem_erase.2 ⟨hue, huS⟩, v, Finset.mem_erase.2 ⟨hve, hvS⟩, huv⟩

/-! ## Existence -/

theorem NonCollinear.nonempty {S : Finset (ℝ × ℝ)} (hS : NonCollinear S) : S.Nonempty := by
  obtain ⟨a, ha, -⟩ := hS
  exact ⟨a, ha⟩

theorem exists_triangulation_of_card (n : ℕ) :
    ∀ S : Finset (ℝ × ℝ), S.card = n → NonCollinear S → ∃ T, IsTriangulation S T := by
  classical
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro S hn hS
    by_cases hle : S.card ≤ 3
    · obtain ⟨a, ha, b, hb, c, hc, habc⟩ := hS
      have hab := sda_ne_zero_left habc
      have hbc := sda_ne_zero_right habc
      have hac := sda_ne_zero_outer habc
      have h3 : ({a, b, c} : Finset (ℝ × ℝ)) = S := by
        refine Finset.eq_of_subset_of_card_le ?_ ?_
        · intro x hx
          simp only [Finset.mem_insert, Finset.mem_singleton] at hx
          rcases hx with rfl | rfl | rfl <;> assumption
        · rw [Finset.card_insert_of_notMem (by simp [hab, hac]), Finset.card_pair hbc]
          exact hle
      rw [← h3]
      exact ⟨_, isTriangulation_triangle habc⟩
    · push Not at hle
      obtain ⟨e, he, hext, hnc⟩ := exists_extreme_nonCollinear hS hle
      obtain ⟨T, hT⟩ := ih (S.erase e).card (by rw [Finset.card_erase_of_mem he]; omega) _ rfl hnc
      have := isTriangulation_insert hT hnc.nonempty hext
      rw [Finset.insert_erase he] at this
      exact ⟨_, this⟩

/-- **Every finite point set that is not collinear has a triangulation** (triangles with corners in
`S`, covering `conv S`, with disjoint interiors, and with every point of `S` that lies in a triangle a
corner of it; `IsTriangulation`).

[formula: planar-triangulation-existence] -/
theorem exists_triangulation (S : Finset (ℝ × ℝ))
    (hS : ∃ a ∈ S, ∃ b ∈ S, ∃ c ∈ S, sda a b c ≠ 0) : ∃ T, IsTriangulation S T :=
  exists_triangulation_of_card _ S rfl hS

/-- **Every finite point set that is not collinear has a Delaunay triangulation**: a triangulation
with no point of `S` strictly inside the circumcircle of any triangle. (A first triangulation,
`exists_triangulation`, then Lawson's flips, `exists_flips_to_delaunay`.)

[formula: planar-delaunay-existence] -/
theorem exists_delaunay (S : Finset (ℝ × ℝ)) (hS : ∃ a ∈ S, ∃ b ∈ S, ∃ c ∈ S, sda a b c ≠ 0) :
    ∃ T, IsTriangulation S T ∧ IsDelaunay S T :=
  exists_delaunay_of_triangulation (exists_triangulation S hS)

/-- The hypothesis is necessary: a nonempty point set with a triangulation is not collinear. -/
theorem nonCollinear_of_isTriangulation {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (hS : S.Nonempty) : NonCollinear S := by
  obtain ⟨s, hs⟩ := hS
  have : s ∈ convexHull ℝ (S : Set (ℝ × ℝ)) := subset_convexHull ℝ _ hs
  rw [← hT.cover] at this
  obtain ⟨t, ht, -⟩ := Set.mem_iUnion₂.1 this
  obtain ⟨a, b, c, hD, rfl⟩ := hT.tri t ht
  exact ⟨a, hT.mem_S ht (by simp), b, hT.mem_S ht (by simp), c, hT.mem_S ht (by simp), hD⟩

/-- The empty point set is triangulated by no triangles. -/
theorem isTriangulation_empty : IsTriangulation ∅ ∅ := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> simp

/-- **A finite point set has a triangulation iff it is empty or not collinear.** -/
theorem exists_triangulation_iff (S : Finset (ℝ × ℝ)) :
    (∃ T, IsTriangulation S T) ↔ S = ∅ ∨ ∃ a ∈ S, ∃ b ∈ S, ∃ c ∈ S, sda a b c ≠ 0 := by
  constructor
  · rintro ⟨T, hT⟩
    rcases S.eq_empty_or_nonempty with h | h
    · exact Or.inl h
    · exact Or.inr (nonCollinear_of_isTriangulation hT h)
  · rintro (rfl | h)
    · exact ⟨∅, isTriangulation_empty⟩
    · exact exists_triangulation S h

end Planar

end Delaunay

end DdgProofs
