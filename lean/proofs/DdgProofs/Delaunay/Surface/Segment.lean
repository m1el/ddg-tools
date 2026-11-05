import DdgProofs.Delaunay.Surface.Plane

/-!
# Delaunay layer 4, planar toolkit: straight segments through triangles

The lemmas the walk of `Overlay.lean` needs about the ray `z t = p + t • d` and a positively
oriented triangle `T` (`Plane.lean`):

* the parameters at which the ray lies in `T` form a closed interval (`exists_exit`); the exit
  point is not interior (`exit_not_mem_int`); between a point of `T` and an interior point the
  ray is interior (`mem_int_of_between`);
* a ray that passes through the interior and leaves through the open side `k` crosses its line
  transversally, towards the outside (`exit_side_cross_neg`, `not_mem_hull_after_side`), and a
  ray that crosses the line of the open side `k` towards the inside enters the interior
  (`mem_int_after_side`); it cannot enter and leave through the same side
  (`side_ne_of_int_between`);
* the triangle unfolded across side `k` (`T'` with `T' k' = T (k+1)`, `T' (k'+1) = T k`) shares
  the open side (`side_across`), has the opposite numerator (`bary_across`), and lies on the other
  side (`not_mem_hull_of_mem_int_across`);
* a ray from a corner that never enters the interior leaves at a corner
  (`exit_corner_of_not_int`);
* finitely many closed convex sets covering a segment: one of them contains an initial piece
  (`exists_initial`); finitely many closed sets: those containing points near `z` contain `z`
  (`exists_ball_subset`); a direction into the wedge of two transversal lines
  (`exists_wedge_dir`).
-/

noncomputable section

namespace DDG.Delaunay.Surface

open Set Filter Topology

/-- The ray `p + t • d`. -/
def ray (p d : Pt) (t : ℝ) : Pt := p + t • d

theorem ray_zero (p d : Pt) : ray p d 0 = p := by simp [ray]

theorem ray_one (p d : Pt) : ray p d 1 = p + d := by simp [ray]

theorem continuous_ray (p d : Pt) : Continuous (ray p d) := by unfold ray; fun_prop

theorem ray_combo (p d : Pt) {t₀ t₁ a b : ℝ} (h : a + b = 1) :
    ray p d (a * t₀ + b * t₁) = a • ray p d t₀ + b • ray p d t₁ := by
  ext
  · simp only [ray, Prod.fst_add, Prod.smul_fst, smul_eq_mul]
    linear_combination (-p.1) * h
  · simp only [ray, Prod.snd_add, Prod.smul_snd, smul_eq_mul]
    linear_combination (-p.2) * h

/-- A convex set containing the ray at `s` and at `u` contains it on `[s, u]`. -/
private theorem ray_mem_of_le {C : Set Pt} (hC : Convex ℝ C) {p d : Pt} {s u t : ℝ}
    (hs : ray p d s ∈ C) (hu : ray p d u ∈ C) (h1 : s ≤ t) (h2 : t ≤ u) : ray p d t ∈ C := by
  rcases h1.eq_or_lt with rfl | h1'
  · exact hs
  have hsu : 0 < u - s := by linarith
  have hab : (u - t) / (u - s) + (t - s) / (u - s) = 1 := by
    rw [← add_div, div_eq_one_iff_eq hsu.ne']; ring
  have key : t = (u - t) / (u - s) * s + (t - s) / (u - s) * u := by
    rw [div_mul_eq_mul_div, div_mul_eq_mul_div, ← add_div, eq_div_iff hsu.ne']; ring
  rw [key, ray_combo p d hab]
  exact hC hs hu (div_nonneg (by linarith) hsu.le) (div_nonneg (by linarith) hsu.le) hab

theorem ray_injective {p d : Pt} (hd : d ≠ 0) : Function.Injective (ray p d) := by
  intro s t h
  simp only [ray, add_right_inj] at h
  exact smul_left_injective ℝ hd h

namespace Tri

variable {T T' : Tri} {p d : Pt}

/-- The numerator of corner `k + 2` is affine along a ray, with slope `cross2 (T (k+1) − T k) d`
(the side `k`, which it vanishes on). -/
theorem bary_ray (k : Fin 3) (t : ℝ) :
    T.bary (k + 2) (ray p d t) = T.bary (k + 2) p + t * cross2 (T (k + 1) - T k) d := by
  rw [ray, bary_affine, show k + 2 + 2 = k + 1 by omega, show k + 2 + 1 = k by omega]

/-- Some side of `T` is crossed outwards by a nonzero direction. -/
theorem exists_cross_neg (hT : T.Pos) (hd : d ≠ 0) : ∃ k, cross2 (T (k + 1) - T k) d < 0 := by
  by_contra hne
  push Not at hne
  have h0 := hne 0
  have h1 := hne 1
  have h2 := hne 2
  simp only [Fin.isValue, Fin.reduceAdd] at h0 h1 h2
  have hsum : cross2 (T 1 - T 0) d + cross2 (T 2 - T 1) d + cross2 (T 0 - T 2) d = 0 := by
    simp only [cross2, Prod.fst_sub, Prod.snd_sub]; ring
  have e0 : cross2 (T 1 - T 0) d = 0 := by linarith
  have e2 : cross2 (T 0 - T 2) d = 0 := by linarith
  have hT' : 0 < sda (T 0) (T 1) (T 2) := hT
  simp only [cross2, Prod.fst_sub, Prod.snd_sub] at e0 e2
  simp only [sda] at hT'
  apply hd
  have y1 : d.1 * (((T 1).1 - (T 0).1) * ((T 2).2 - (T 0).2) -
      ((T 1).2 - (T 0).2) * ((T 2).1 - (T 0).1)) = 0 := by
    linear_combination ((T 2).1 - (T 0).1) * e0 + ((T 1).1 - (T 0).1) * e2
  have y2 : d.2 * (((T 1).1 - (T 0).1) * ((T 2).2 - (T 0).2) -
      ((T 1).2 - (T 0).2) * ((T 2).1 - (T 0).1)) = 0 := by
    linear_combination ((T 2).2 - (T 0).2) * e0 + ((T 1).2 - (T 0).2) * e2
  ext
  · exact (mul_eq_zero.1 y1).resolve_right hT'.ne'
  · exact (mul_eq_zero.1 y2).resolve_right hT'.ne'

/-- The parameters at which the ray lies in the closed triangle, from `t₀` on, form a closed
interval `[t₀, t₁]`. -/
theorem exists_exit (hT : T.Pos) (hd : d ≠ 0) {t₀ : ℝ} (h₀ : ray p d t₀ ∈ T.hull) :
    ∃ t₁, t₀ ≤ t₁ ∧ ∀ t, (t₀ ≤ t ∧ ray p d t ∈ T.hull) ↔ t ∈ Icc t₀ t₁ := by
  obtain ⟨k, hk⟩ := exists_cross_neg hT hd
  have hne : ({t | t₀ ≤ t ∧ ray p d t ∈ T.hull} : Set ℝ).Nonempty := ⟨t₀, le_rfl, h₀⟩
  have hbdd : BddAbove {t | t₀ ≤ t ∧ ray p d t ∈ T.hull} := by
    refine ⟨T.bary (k + 2) p / -cross2 (T (k + 1) - T k) d, ?_⟩
    rintro t ⟨-, ht⟩
    have h := (mem_hull_iff hT).1 ht (k + 2)
    rw [bary_ray] at h
    rw [le_div_iff₀ (by linarith)]
    nlinarith
  have hcl : IsClosed {t | t₀ ≤ t ∧ ray p d t ∈ T.hull} :=
    isClosed_Ici.inter (isClosed_hull.preimage (continuous_ray p d))
  refine ⟨sSup {t | t₀ ≤ t ∧ ray p d t ∈ T.hull}, le_csSup hbdd ⟨le_rfl, h₀⟩,
    fun t => ⟨fun ht => ⟨ht.1, le_csSup hbdd ht⟩, fun ht => ⟨ht.1, ?_⟩⟩⟩
  exact ray_mem_of_le convex_hull h₀ (hcl.csSup_mem hne hbdd).2 ht.1 ht.2

/-- Between a point of the closed triangle and an interior point the ray is interior. -/
theorem mem_int_of_between (hT : T.Pos) {t₀ t₁ t : ℝ} (h₀ : ray p d t₀ ∈ T.hull)
    (h₁ : ray p d t₁ ∈ T.int) (ht : t ∈ Ioo (min t₀ t₁) (max t₀ t₁)) : ray p d t ∈ T.int := by
  have hcl : ray p d t₀ ∈ closure T.hull := subset_closure h₀
  apply convex_hull.openSegment_interior_closure_subset_interior h₁ hcl
  rw [openSegment_eq_image]
  rcases lt_trichotomy t₀ t₁ with hlt | heq | hgt
  · rw [min_eq_left hlt.le, max_eq_right hlt.le] at ht
    have hpos : 0 < t₁ - t₀ := by linarith
    have hθ : (t₁ - t) / (t₁ - t₀) * (t₁ - t₀) = t₁ - t := div_mul_cancel₀ _ hpos.ne'
    refine ⟨(t₁ - t) / (t₁ - t₀), ⟨div_pos (by linarith [ht.2]) hpos,
      (div_lt_one hpos).2 (by linarith [ht.1])⟩, ?_⟩
    show (1 - (t₁ - t) / (t₁ - t₀)) • ray p d t₁ + ((t₁ - t) / (t₁ - t₀)) • ray p d t₀ =
      ray p d t
    rw [← ray_combo p d (by ring)]
    congr 1
    linear_combination (-1 : ℝ) * hθ
  · subst heq
    rw [min_self, max_self] at ht
    exact absurd (ht.1.trans ht.2) (lt_irrefl _)
  · rw [min_eq_right hgt.le, max_eq_left hgt.le] at ht
    have hpos : 0 < t₀ - t₁ := by linarith
    have hθ : (t - t₁) / (t₀ - t₁) * (t₀ - t₁) = t - t₁ := div_mul_cancel₀ _ hpos.ne'
    refine ⟨(t - t₁) / (t₀ - t₁), ⟨div_pos (by linarith [ht.1]) hpos,
      (div_lt_one hpos).2 (by linarith [ht.2])⟩, ?_⟩
    show (1 - (t - t₁) / (t₀ - t₁)) • ray p d t₁ + ((t - t₁) / (t₀ - t₁)) • ray p d t₀ =
      ray p d t
    rw [← ray_combo p d (by ring)]
    congr 1
    linear_combination hθ

/-- The exit point is not interior, for `t₀ ≤ t₁`. -/
private theorem exit_not_mem_int_of_le {t₀ t₁ : ℝ} (hle : t₀ ≤ t₁)
    (hS : ∀ t, (t₀ ≤ t ∧ ray p d t ∈ T.hull) ↔ t ∈ Icc t₀ t₁) : ray p d t₁ ∉ T.int := by
  intro h
  have hm : ray p d ⁻¹' T.int ∈ 𝓝 t₁ :=
    (continuous_ray p d).continuousAt.preimage_mem_nhds (isOpen_int.mem_nhds h)
  obtain ⟨ε, hε, hball⟩ := Metric.mem_nhds_iff.1 hm
  have hdist : dist (t₁ + ε / 2) t₁ < ε := by
    rw [Real.dist_eq, abs_of_pos (by linarith)]; linarith
  have hi : ray p d (t₁ + ε / 2) ∈ T.int := hball hdist
  have := ((hS (t₁ + ε / 2)).1 ⟨by linarith, int_subset_hull hi⟩).2
  linarith

/-- The exit point is not interior. -/
theorem exit_not_mem_int (hT : T.Pos) {t₀ t₁ : ℝ} (hle : t₀ ≤ t₁)
    (hS : ∀ t, (t₀ ≤ t ∧ ray p d t ∈ T.hull) ↔ t ∈ Icc t₀ t₁) : ray p d t₁ ∉ T.int :=
  exit_not_mem_int_of_le hle hS

/-- A ray that is interior on `(t₀, t₁)` and on the open side `k` at `t₁` crosses the line of
that side outwards. -/
theorem exit_side_cross_neg (hT : T.Pos) {t₀ t₁ : ℝ} {k : Fin 3} (hlt : t₀ < t₁)
    (hint : ∀ t ∈ Ioo t₀ t₁, ray p d t ∈ T.int) (hside : ray p d t₁ ∈ T.side k) :
    cross2 (T (k + 1) - T k) d < 0 := by
  have e1 := ((mem_side_iff hT).1 hside).1
  have hm : (t₀ + t₁) / 2 ∈ Ioo t₀ t₁ := ⟨by linarith, by linarith⟩
  have e2 := (mem_int_iff hT).1 (hint _ hm) (k + 2)
  rw [bary_ray] at e1 e2
  have h3 : (t₁ - (t₀ + t₁) / 2) * cross2 (T (k + 1) - T k) d < 0 := by linarith
  by_contra hc
  push Not at hc
  exact absurd h3 (not_lt.2 (mul_nonneg (by linarith) hc))

/-- After crossing the line of side `k` outwards the ray is outside the closed triangle. -/
theorem not_mem_hull_after_side (hT : T.Pos) {t₁ : ℝ} {k : Fin 3} (hside : ray p d t₁ ∈ T.side k)
    (hc : cross2 (T (k + 1) - T k) d < 0) {t : ℝ} (ht : t₁ < t) : ray p d t ∉ T.hull := by
  intro h
  have e1 := ((mem_side_iff hT).1 hside).1
  have e2 := (mem_hull_iff hT).1 h (k + 2)
  rw [bary_ray] at e1 e2
  have h3 : (t - t₁) * cross2 (T (k + 1) - T k) d < 0 := mul_neg_of_pos_of_neg (by linarith) hc
  linarith

/-- A ray that is on the open side `k` at `t₀` and interior on `(t₀, t₁)` crosses the line of
that side inwards. -/
theorem entry_side_cross_pos (hT : T.Pos) {t₀ t₁ : ℝ} {k : Fin 3} (hlt : t₀ < t₁)
    (hint : ∀ t ∈ Ioo t₀ t₁, ray p d t ∈ T.int) (hside : ray p d t₀ ∈ T.side k) :
    0 < cross2 (T (k + 1) - T k) d := by
  have e1 := ((mem_side_iff hT).1 hside).1
  have hm : (t₀ + t₁) / 2 ∈ Ioo t₀ t₁ := ⟨by linarith, by linarith⟩
  have e2 := (mem_int_iff hT).1 (hint _ hm) (k + 2)
  rw [bary_ray] at e1 e2
  have h3 : 0 < ((t₀ + t₁) / 2 - t₀) * cross2 (T (k + 1) - T k) d := by linarith
  by_contra hc
  push Not at hc
  exact absurd h3 (not_lt.2 (mul_nonpos_of_nonneg_of_nonpos (by linarith) hc))

/-- Before crossing the line of side `k` inwards the ray is outside the closed triangle. -/
theorem not_mem_hull_before_side (hT : T.Pos) {t₀ : ℝ} {k : Fin 3} (hside : ray p d t₀ ∈ T.side k)
    (hc : 0 < cross2 (T (k + 1) - T k) d) {t : ℝ} (ht : t < t₀) : ray p d t ∉ T.hull := by
  intro h
  have e1 := ((mem_side_iff hT).1 hside).1
  have e2 := (mem_hull_iff hT).1 h (k + 2)
  rw [bary_ray] at e1 e2
  have h3 : (t - t₀) * cross2 (T (k + 1) - T k) d < 0 := mul_neg_of_neg_of_pos (by linarith) hc
  linarith

/-- Two numerators positive at `t₁` stay positive near `t₁`. -/
private theorem exists_near_pos (j₁ j₂ : Fin 3) {t₁ : ℝ} (h₁ : 0 < T.bary j₁ (ray p d t₁))
    (h₂ : 0 < T.bary j₂ (ray p d t₁)) :
    ∃ ε > 0, ∀ t, |t - t₁| < ε → 0 < T.bary j₁ (ray p d t) ∧ 0 < T.bary j₂ (ray p d t) := by
  have hc : ∀ j, Continuous fun t => T.bary j (ray p d t) :=
    fun j => (continuous_bary (T := T) j).comp (continuous_ray p d)
  have hm : {t | 0 < T.bary j₁ (ray p d t)} ∩ {t | 0 < T.bary j₂ (ray p d t)} ∈ 𝓝 t₁ :=
    Filter.inter_mem ((isOpen_lt continuous_const (hc j₁)).mem_nhds h₁)
      ((isOpen_lt continuous_const (hc j₂)).mem_nhds h₂)
  obtain ⟨ε, hε, hball⟩ := Metric.mem_nhds_iff.1 hm
  exact ⟨ε, hε, fun t ht => hball (by rw [Metric.mem_ball, Real.dist_eq]; exact ht)⟩

/-- A ray through the open side `k` at `t₁` that crosses its line inwards is interior just
after `t₁`. -/
theorem mem_int_after_side (hT : T.Pos) {t₁ : ℝ} {k : Fin 3} (hside : ray p d t₁ ∈ T.side k)
    (hc : 0 < cross2 (T (k + 1) - T k) d) :
    ∃ ε > 0, ∀ t ∈ Ioo t₁ (t₁ + ε), ray p d t ∈ T.int := by
  obtain ⟨e0, e1, e2⟩ := (mem_side_iff hT).1 hside
  obtain ⟨ε, hε, hnear⟩ := exists_near_pos k (k + 1) e1 e2
  refine ⟨ε, hε, fun t ht => ?_⟩
  have hn := hnear t (by rw [abs_of_pos (by linarith [ht.1])]; linarith [ht.2])
  rw [mem_int_iff hT]
  intro j
  rcases (by omega : j = k ∨ j = k + 1 ∨ j = k + 2) with rfl | rfl | rfl
  · exact hn.1
  · exact hn.2
  · rw [bary_ray] at e0 ⊢
    have h3 : 0 < (t - t₁) * cross2 (T (k + 1) - T k) d := mul_pos (by linarith [ht.1]) hc
    linarith

/-- A ray cannot be on the same open side at two parameters and interior in between. -/
theorem side_ne_of_int_between (hT : T.Pos) {t₀ t₁ t : ℝ} {k : Fin 3} (h₀ : ray p d t₀ ∈ T.side k)
    (h₁ : ray p d t₁ ∈ T.side k) (hlt : t₀ < t₁) (ht : t ∈ Ioo t₀ t₁)
    (hint : ray p d t ∈ T.int) : False := by
  have e0 := ((mem_side_iff hT).1 h₀).1
  have e1 := ((mem_side_iff hT).1 h₁).1
  have e := (mem_int_iff hT).1 hint (k + 2)
  rw [bary_ray] at e0 e1 e
  have hc : cross2 (T (k + 1) - T k) d = 0 := by
    have h : (t₁ - t₀) * cross2 (T (k + 1) - T k) d = 0 := by linarith
    exact (mul_eq_zero.1 h).resolve_left (by linarith)
  rw [hc, mul_zero, add_zero] at e0 e
  linarith

/-- A ray from the corner `k` of `T` (`p = T k`) that is never interior on `(0, t₁)`, with
`[0, t₁]` its parameters inside `T` and `0 < t₁`, leaves `T` at a corner. -/
theorem exit_corner_of_not_int (hT : T.Pos) (hd : d ≠ 0) {k : Fin 3} {t₁ : ℝ} (hp : p = T k)
    (hS : ∀ t, (0 ≤ t ∧ ray p d t ∈ T.hull) ↔ t ∈ Icc 0 t₁) (ht₁ : 0 < t₁)
    (hni : ∀ t ∈ Ioo 0 t₁, ray p d t ∉ T.int) : ∃ j, ray p d t₁ = T j := by
  have hh : ray p d t₁ ∈ T.hull := ((hS t₁).2 ⟨ht₁.le, le_rfl⟩).2
  have hi : ray p d t₁ ∉ T.int := exit_not_mem_int_of_le ht₁.le hS
  by_contra hc
  push Not at hc
  obtain ⟨j, hj⟩ := exists_side_of_mem_hull hT hh hi hc
  obtain ⟨e0, e1, e2⟩ := (mem_side_iff hT).1 hj
  obtain ⟨ε, hε, hnear⟩ := exists_near_pos j (j + 1) e1 e2
  rw [bary_ray] at e0
  rcases eq_or_ne (cross2 (T (j + 1) - T j) d) 0 with hz | hz
  · -- the ray runs along the line of side `j`: it is still in `T` just after `t₁`
    have hn := hnear (t₁ + ε / 2) (by rw [abs_of_pos (by linarith)]; linarith)
    have hm : ray p d (t₁ + ε / 2) ∈ T.hull := by
      rw [mem_hull_iff hT]
      intro i
      rcases (by omega : i = j ∨ i = j + 1 ∨ i = j + 2) with rfl | rfl | rfl
      · exact hn.1.le
      · exact hn.2.le
      · rw [bary_ray, hz]; rw [hz] at e0; linarith
    have := ((hS (t₁ + ε / 2)).1 ⟨by linarith, hm⟩).2
    linarith
  · -- the ray crosses the line of side `j` outwards: it is interior just before `t₁`
    have h0 := (mem_hull_iff hT).1 ((hS 0).2 ⟨le_rfl, ht₁.le⟩).2 (j + 2)
    rw [bary_ray, zero_mul, add_zero] at h0
    have hneg : cross2 (T (j + 1) - T j) d < 0 := by
      rcases hz.lt_or_gt with h | h
      · exact h
      · nlinarith
    have hmin : 0 < min ε t₁ := lt_min hε ht₁
    have hmin1 : min ε t₁ ≤ ε := min_le_left _ _
    have hmin2 : min ε t₁ ≤ t₁ := min_le_right _ _
    have hsI : t₁ - min ε t₁ / 2 ∈ Ioo 0 t₁ := ⟨by linarith, by linarith⟩
    have hn := hnear (t₁ - min ε t₁ / 2) (by rw [abs_of_neg (by linarith)]; linarith)
    apply hni _ hsI
    rw [mem_int_iff hT]
    intro i
    rcases (by omega : i = j ∨ i = j + 1 ∨ i = j + 2) with rfl | rfl | rfl
    · exact hn.1
    · exact hn.2
    · rw [bary_ray]
      have h3 : 0 < (t₁ - min ε t₁ / 2 - t₁) * cross2 (T (j + 1) - T j) d :=
        mul_pos_of_neg_of_neg (by linarith) hneg
      linarith

/-! ### The triangle unfolded across a side -/

section Across

variable {k k' : Fin 3} (h0 : T' k' = T (k + 1)) (h1 : T' (k' + 1) = T k)
include h0 h1

theorem side_across : T'.side k' = T.side k := by
  simp only [side, h0, h1]
  exact openSegment_symm ℝ _ _

theorem bary_across (x : Pt) : T'.bary (k' + 2) x = -T.bary (k + 2) x := by
  simp only [bary, show k' + 2 + 1 = k' by omega, show k' + 2 + 2 = k' + 1 by omega,
    show k + 2 + 1 = k by omega, show k + 2 + 2 = k + 1 by omega, h0, h1]
  exact sda_swap_left _ _ _

theorem cross_across (w : Pt) :
    cross2 (T' (k' + 1) - T' k') w = -cross2 (T (k + 1) - T k) w := by
  rw [h0, h1]
  simp only [cross2, Prod.fst_sub, Prod.snd_sub]
  ring

omit h0 h1 in
/-- A point interior to the unfolded triangle is outside the closed original one. -/
theorem not_mem_hull_of_mem_int_across (hT : T.Pos) (hT' : T'.Pos) (h0 : T' k' = T (k + 1))
    (h1 : T' (k' + 1) = T k) {x : Pt} (hx : x ∈ T'.int) : x ∉ T.hull := by
  intro hx'
  rw [mem_int_iff hT'] at hx
  rw [mem_hull_iff hT] at hx'
  have := hx (k' + 2)
  rw [bary_across h0 h1] at this
  linarith [hx' (k + 2)]

end Across

end Tri

/-! ### Finite covers -/

/-- Finitely many closed convex sets covering the segment `ray p d '' [0, 1]`: one of them
contains an initial piece `ray p d '' [0, δ]`. -/
theorem exists_initial {ι : Type*} (S : Finset ι) (C : ι → Set Pt)
    (hcl : ∀ i ∈ S, IsClosed (C i)) (hcv : ∀ i ∈ S, Convex ℝ (C i)) {p d : Pt}
    (hcov : ∀ t ∈ Icc (0 : ℝ) 1, ∃ i ∈ S, ray p d t ∈ C i) :
    ∃ i ∈ S, ∃ δ > 0, ∀ t ∈ Icc 0 δ, ray p d t ∈ C i := by
  by_contra H
  push Not at H
  have hev : ∀ i ∈ S, ∀ᶠ t in 𝓝[>] (0 : ℝ), ray p d t ∉ C i := by
    intro i hi
    by_cases h0 : ray p d 0 ∈ C i
    · refine eventually_nhdsWithin_of_forall fun t (ht : 0 < t) hrt => ?_
      obtain ⟨s, hs, hns⟩ := H i hi t ht
      exact hns (ray_mem_of_le (hcv i hi) h0 hrt hs.1 hs.2)
    · exact eventually_nhdsWithin_of_eventually_nhds
        (((hcl i hi).isOpen_compl.preimage (continuous_ray p d)).mem_nhds h0)
  have hall := (Filter.eventually_all_finset S).2 hev
  have hlt : ∀ᶠ t in 𝓝[>] (0 : ℝ), t < 1 :=
    eventually_nhdsWithin_of_eventually_nhds (eventually_lt_nhds zero_lt_one)
  have hpos : ∀ᶠ t in 𝓝[>] (0 : ℝ), t ∈ Ioi 0 := self_mem_nhdsWithin
  obtain ⟨t, ht, ht1, ht0⟩ := (hall.and (hlt.and hpos)).exists
  obtain ⟨i, hi, hti⟩ := hcov t ⟨le_of_lt ht0, ht1.le⟩
  exact ht i hi hti

/-- Finitely many closed sets: every point near enough to `z` lies only in sets that contain
`z`. -/
theorem exists_ball_subset {ι : Type*} (S : Finset ι) (C : ι → Set Pt)
    (hcl : ∀ i ∈ S, IsClosed (C i)) (z : Pt) :
    ∃ ε > 0, ∀ y, dist y z < ε → ∀ i ∈ S, y ∈ C i → z ∈ C i := by
  have hev : ∀ i ∈ S, ∀ᶠ y in 𝓝 z, y ∈ C i → z ∈ C i := by
    intro i hi
    by_cases hz : z ∈ C i
    · exact Filter.Eventually.of_forall fun _ _ => hz
    · exact Filter.mem_of_superset ((hcl i hi).isOpen_compl.mem_nhds hz)
        fun y hy hyC => absurd hyC hy
  obtain ⟨ε, hε, h⟩ := Metric.eventually_nhds_iff.1 ((Filter.eventually_all_finset S).2 hev)
  exact ⟨ε, hε, fun y hy => h hy⟩

/-- A direction strictly on the right of `a` and strictly on the left of `b`, for two
non-parallel vectors. -/
theorem exists_wedge_dir {a b : Pt} (h : cross2 a b ≠ 0) :
    ∃ w : Pt, cross2 a w < 0 ∧ 0 < cross2 b w := by
  refine ⟨(-(1 / cross2 a b)) • (a + b), ?_, ?_⟩
  · have e : cross2 a ((-(1 / cross2 a b)) • (a + b)) = -(1 / cross2 a b * cross2 a b) := by
      simp only [cross2, Prod.smul_fst, Prod.smul_snd, Prod.fst_add, Prod.snd_add, smul_eq_mul]
      ring
    rw [e, one_div, inv_mul_cancel₀ h]
    norm_num
  · have e : cross2 b ((-(1 / cross2 a b)) • (a + b)) = 1 / cross2 a b * cross2 a b := by
      simp only [cross2, Prod.smul_fst, Prod.smul_snd, Prod.fst_add, Prod.snd_add, smul_eq_mul]
      ring
    rw [e, one_div, inv_mul_cancel₀ h]
    norm_num

end DDG.Delaunay.Surface
