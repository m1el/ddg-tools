import DdgProofs.Delaunay.Surface.FlipFace
import DdgProofs.Delaunay.Surface.Relayout

/-!
# Delaunay layer 4: the flip of an overlaid triangulation, edge level and assembly

Continues `FlipFace.lean`. The halfedges of the flipped triangulation `t'` are those of `t`
renamed by the 4-cycle `σ = flipPerm h₀` (`DdgProofs/Delaunay/Flip.lean`): `σ h₀ = h₀.next` is
the new edge `l → k`, and every **outer** side of the diamond (`h₀.next = j→k`, `h₀.prev = k→i`,
`(twin h₀).next = i→l`, `(twin h₀).prev = l→j`) keeps its two end points under the new layouts
(`outer_lay`); the cells of the new face at a point of an outer side are the cells of the old
face of that side (`outer_germAt`). Far halfedges (in faces other than `f₀, f₁`) are untouched.
So the matching motion `layΔ` of every edge of `t'` other than the new one is the old one
(`layΔ_flip_eq`), the edge germs are the old germs, and the new edge has the identity motion
with the germ of the whole diamond (`germ_D`). The symmetric statements for `f₁` come from
`FlipOK.symm` and `flip_twin_eq` (the flip of the twin halfedge is the same flip).

`FlipOK.flipOverlay` assembles the overlay of `t'`; `Overlay.flip` (for any overlay, after
relaying `f₁` into the plane of `f₀`, `Relayout.lean`) is the theorem **a convex flip keeps an
overlay**.
-/

noncomputable section

namespace DdgProofs.Delaunay.Surface

open DdgProofs.Intrinsic (dist2_comm)
open DdgProofs.Delaunay.Planar
open Set
open scoped Classical Topology

variable {F : Type*}

/-! ## The flip of the twin halfedge is the same flip -/

namespace Gluing

variable [DecidableEq F] {G : Gluing F} {h₀ : Halfedge F}

omit [DecidableEq F] in
theorem flippable_twin (hf : G.Flippable h₀) : G.Flippable (G.twin h₀) := by
  unfold Gluing.Flippable; rw [G.twin_twin]; exact Ne.symm hf

theorem flipPerm_twin_eq (hf : G.Flippable h₀) : G.flipPerm (G.twin h₀) = G.flipPerm h₀ := by
  have hf' := flippable_twin hf
  have e0 := Gluing.flipPerm_twin hf'
  have e1 := Gluing.flipPerm_twin_next hf'
  have e2 := Gluing.flipPerm_self hf'
  have e3 := Gluing.flipPerm_next hf'
  rw [G.twin_twin] at e0 e1 e3
  refine Equiv.ext fun x => ?_
  by_cases x0 : x = h₀
  · rw [x0, e0, Gluing.flipPerm_self hf]
  by_cases x1 : x = h₀.next
  · rw [x1, e1, Gluing.flipPerm_next hf]
  by_cases x2 : x = G.twin h₀
  · rw [x2, e2, Gluing.flipPerm_twin hf]
  by_cases x3 : x = (G.twin h₀).next
  · rw [x3, e3, Gluing.flipPerm_twin_next hf]
  rw [Gluing.flipPerm_of_ne x0 x1 x2 x3,
    Gluing.flipPerm_of_ne (h₀ := G.twin h₀) x2 x3 (by rwa [G.twin_twin])
      (by rwa [G.twin_twin])]

theorem flip_twin_eq (hf : G.Flippable h₀) : G.flip (G.twin h₀) = G.flip h₀ := by
  have h := flipPerm_twin_eq hf
  have ext : ∀ G₁ G₂ : Gluing F, G₁.twin = G₂.twin → G₁ = G₂ := by
    rintro ⟨a, _⟩ ⟨b, _⟩ e
    cases e
    rfl
  apply ext
  funext x
  show G.flipPerm (G.twin h₀) (G.twin ((G.flipPerm (G.twin h₀)).symm x)) =
    G.flipPerm h₀ (G.twin ((G.flipPerm h₀).symm x))
  rw [h]

theorem newLength_twin {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (hf : G.Flippable h₀) :
    G.newLength ℓ (G.twin h₀) = G.newLength ℓ h₀ := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  unfold Gluing.newLength
  rw [G.twin_twin, hm.twin h₀, flipLength_eq_flippedLengthFromLengths hl.swap hk.swap,
    flipLength_eq_flippedLengthFromLengths hk hl]
  exact (DdgProofs.Intrinsic.flippedLengthFromLengths_comm _ _ _ _ _).trans
    (DdgProofs.Intrinsic.flippedLengthFromLengths_mirror hk.pos_a.ne' _ _ _ _)

theorem flipLen_twin_eq {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (hf : G.Flippable h₀) :
    G.flipLen ℓ (G.twin h₀) = G.flipLen ℓ h₀ := by
  funext x
  simp only [Gluing.flipLen, flipPerm_twin_eq hf, newLength_twin hm hf, G.twin_twin]
  by_cases h1 : x = h₀.next ∨ x = (G.twin h₀).next
  · simp only [h1, h1.symm, ite_true]
  · have h2 : ¬(x = (G.twin h₀).next ∨ x = h₀.next) := fun h => h1 h.symm
    simp only [h1, h2, ite_false]

end Gluing

namespace Overlay

variable {s t : MetricGluing F} {O : Overlay s t} {h₀ : Halfedge F}

/-! ### Outer and far halfedges -/

/-- `m` is an outer side of the diamond of `h₀`. -/
def Outer (t : MetricGluing F) (h₀ m : Halfedge F) : Prop :=
  m = h₀.next ∨ m = h₀.prev ∨ m = (t.G.twin h₀).next ∨ m = (t.G.twin h₀).prev

/-- `m` is far from the diamond of `h₀`: in a face other than `f₀`, `f₁`. -/
def Far (t : MetricGluing F) (h₀ m : Halfedge F) : Prop := m.1 ≠ h₀.1 ∧ m.1 ≠ (t.G.twin h₀).1

theorem far_flipPerm {m : Halfedge F} (hm : Far t h₀ m) :
    t.G.flipPerm h₀ m = m := by
  obtain ⟨h0, h1⟩ := hm
  apply Gluing.flipPerm_of_ne
  · exact fun e => h0 (by rw [e])
  · exact fun e => h0 (by rw [e]; rfl)
  · exact fun e => h1 (by rw [e])
  · exact fun e => h1 (by rw [e]; rfl)

theorem far_flipLay {m : Halfedge F} (hm : Far t h₀ m) :
    O.flipLay h₀ m.1 = O.lay m.1 := flipLay_other hm.1 hm.2

theorem far_flipCells {m : Halfedge F} (hm : Far t h₀ m) :
    O.flipCells h₀ m.1 = O.cells m.1 := flipCells_other hm.1 hm.2

/-- Points of the open segment `(w, y)` come arbitrarily close to `y`. -/
theorem exists_openSegment_near (w y : Pt) {ε : ℝ} (hε : 0 < ε) :
    ∃ z ∈ openSegment ℝ w y, dist z y < ε := by
  have h1 : ∀ᶠ r in 𝓝[>] (0 : ℝ), dist (y + r • (w - y)) y < ε :=
    exists_pos_near (Metric.ball_mem_nhds y hε) (w - y)
  obtain ⟨r, hr, hr1⟩ := (h1.and (Ioo_mem_nhdsGT (by norm_num : (0 : ℝ) < 1))).exists
  refine ⟨y + r • (w - y), ⟨r, 1 - r, hr1.1, by linarith [hr1.2], by ring, ?_⟩, hr⟩
  ext <;> simp <;> ring

/-- A cell containing a point `y` of the open side `m` of a positively oriented triangle `T`,
with an interior point strictly on `T`'s side of that line, overlaps the open triangle. -/
theorem overlap_of_side {c : Cell F} {T : Tri} (hT : T.Pos) {m : Fin 3} {y : Pt}
    (hy : y ∈ T.side m) (hyc : y ∈ (c.tri s).hull) {w : Pt} (hwc : w ∈ (c.tri s).int)
    (hw : 0 < T.bary (m + 2) w) : ((c.tri s).int ∩ T.int).Nonempty := by
  obtain ⟨ε, hε, hball⟩ := Tri.exists_ball_int_of_side hT hy
  obtain ⟨z, hz, hzd⟩ := exists_openSegment_near w y hε
  have hzc : z ∈ (c.tri s).int :=
    Convex.openSegment_interior_closure_subset_interior Tri.convex_hull hwc (subset_closure hyc) hz
  have hy0 : T.bary (m + 2) y = 0 := ((Tri.mem_side_iff hT).1 hy).1
  obtain ⟨a, b, ha, hb, hab, rfl⟩ := hz
  refine ⟨_, hzc, hball _ hzd ?_⟩
  rw [Tri.bary_combo _ _ _ hab, hy0, mul_zero, add_zero]
  exact mul_pos ha hw

/-- The germ of a union with a moved set of cells. -/
theorem germAt_union_image (A B : Finset (Cell F)) (Δ : Rigid) (y : Pt) :
    germAt s (A ∪ B.image (Cell.mv Δ)) y =
      germAt s A y ∪ (germAt s B (Δ.inv y)).image (Cell.mv Δ) := by
  have h1 : germAt s (A ∪ B.image (Cell.mv Δ)) y =
      germAt s A y ∪ germAt s (B.image (Cell.mv Δ)) y := by
    unfold germAt; exact Finset.filter_union _ _ _
  have h2 := germAt_image_mv s B Δ (Δ.inv y)
  rw [Rigid.apply_inv_apply] at h2
  rw [h1, h2]

namespace FlipOK

variable (hF : FlipOK O h₀)
include hF

/-! ### Germs on the sides of the new face `f₀` -/

/-! ### Outer sides keep their end points and their germs -/

/-- On the open side `P L` of the new face `f₀` (the old side `(twin h₀).next` of `f₁`), the
cells of the new face are the old cells of `f₁`. -/
theorem germAt_flipCells_f₀_of_mem_PL {y : Pt} (hy : y ∈ openSegment ℝ (O.dP h₀) (O.dL h₀)) :
    germAt s (O.flipCells h₀ h₀.1) y = germAt s (O.cells (t.G.twin h₀).1) y := by
  have hq := hF.quad
  ext c
  rw [flipCells_f₀, mem_germAt, mem_germAt, Finset.mem_filter]
  constructor
  · rintro ⟨⟨hcD, -⟩, hyc⟩
    refine ⟨?_, hyc⟩
    rcases mem_Dcells.1 hcD with hc0 | hc1
    · have hs : (tri3 (O.dQ h₀) (O.dP h₀) (O.dL h₀)).side 1 = openSegment ℝ (O.dP h₀) (O.dL h₀) :=
        rfl
      have hyL : y ∈ (tri3 (O.dQ h₀) (O.dP h₀) (O.dL h₀)).hull :=
        Tri.side_subset_hull (tri3_pos hq.hL) 1 (hs ▸ hy)
      exact hF.mem_cells_of_sda_neg hcD hyc hyL (hq.sda_neg_of_mem_PL hy)
    · exact hc1
  · rintro ⟨hc1, hyc⟩
    refine ⟨⟨mem_Dcells.2 (Or.inr hc1), ?_⟩, hyc⟩
    obtain ⟨w, hwc, hwf⟩ := O.nonjunk _ c hc1
    rw [hF.lay_f₁_eq, rot_int] at hwf
    have hw := ((mem_int_tri3 hq.hL).1 hwf).1
    exact overlap_of_side (T := tri3 (O.dP h₀) (O.dL h₀) (O.dK h₀)) (m := 0)
      (tri3_pos hq.hc1) hy hyc hwc hw

/-- On the open side `K P` of the new face `f₀` (the old side `h₀.prev` of `f₀`), the cells of
the new face are the old cells of `f₀`. -/
theorem germAt_flipCells_f₀_of_mem_KP {y : Pt} (hy : y ∈ openSegment ℝ (O.dK h₀) (O.dP h₀)) :
    germAt s (O.flipCells h₀ h₀.1) y = germAt s (O.cells h₀.1) y := by
  have hq := hF.quad
  ext c
  rw [flipCells_f₀, mem_germAt, mem_germAt, Finset.mem_filter]
  constructor
  · rintro ⟨⟨hcD, -⟩, hyc⟩
    refine ⟨?_, hyc⟩
    rcases mem_Dcells.1 hcD with hc0 | hc1
    · exact hc0
    · have hs : (tri3 (O.dP h₀) (O.dQ h₀) (O.dK h₀)).side 2 = openSegment ℝ (O.dK h₀) (O.dP h₀) :=
        rfl
      have hyK : y ∈ (tri3 (O.dP h₀) (O.dQ h₀) (O.dK h₀)).hull :=
        Tri.side_subset_hull (tri3_pos hq.hK) 2 (hs ▸ hy)
      exact hF.mem_cells_of_sda_pos hcD hyc hyK (hq.sda_pos_of_mem_KP hy)
  · rintro ⟨hc0, hyc⟩
    refine ⟨⟨mem_Dcells.2 (Or.inl hc0), ?_⟩, hyc⟩
    obtain ⟨w, hwc, hwf⟩ := O.nonjunk _ c hc0
    rw [hF.lay_f₀_eq, rot_int] at hwf
    have hw := ((mem_int_tri3 hq.hK).1 hwf).2.1
    exact overlap_of_side (T := tri3 (O.dP h₀) (O.dL h₀) (O.dK h₀)) (m := 2)
      (tri3_pos hq.hc1) hy hyc hwc hw

/-- On the new open diagonal `L K`, the cells of the two new faces together are the cells of
the diamond. -/
theorem germAt_flipCells_union_of_mem_LK {y : Pt} (hy : y ∈ openSegment ℝ (O.dL h₀) (O.dK h₀)) :
    germAt s (O.flipCells h₀ h₀.1 ∪ O.flipCells h₀ (t.G.twin h₀).1) y =
      germAt s (O.Dcells h₀) y := by
  have hq := hF.quad
  ext c
  rw [mem_germAt, mem_germAt, Finset.mem_union, flipCells_f₀, flipCells_f₁ hF.flippable,
    Finset.mem_filter, Finset.mem_filter]
  constructor
  · rintro ⟨(⟨hcD, -⟩ | ⟨hcD, -⟩), hyc⟩ <;> exact ⟨hcD, hyc⟩
  · rintro ⟨hcD, hyc⟩
    refine ⟨?_, hyc⟩
    have hy1 : y ∈ (tri3 (O.dP h₀) (O.dL h₀) (O.dK h₀)).side 1 := hy
    have hy2 : y ∈ (tri3 (O.dQ h₀) (O.dK h₀) (O.dL h₀)).side 1 := by
      show y ∈ openSegment ℝ (O.dK h₀) (O.dL h₀)
      rw [openSegment_symm]; exact hy
    obtain ⟨ε₁, hε₁, hb₁⟩ := Tri.exists_ball_int_of_side (tri3_pos hq.hc1) hy1
    obtain ⟨ε₂, hε₂, hb₂⟩ := Tri.exists_ball_int_of_side (tri3_pos hq.hc2) hy2
    have hε : 0 < min ε₁ ε₂ := lt_min hε₁ hε₂
    obtain ⟨b, hbc, hbd⟩ :=
      Metric.mem_closure_iff.1 (Tri.hull_subset_closure_int (Cell.tri_pos c) hyc) _ hε
    have hU : IsOpen ((c.tri s).int ∩ Metric.ball y (min ε₁ ε₂)) :=
      Tri.isOpen_int.inter Metric.isOpen_ball
    obtain ⟨z, ⟨hzc, hzb⟩, hz⟩ := exists_mem_sda_ne hq.ne_LK hU
      ⟨hbc, by rw [Metric.mem_ball, dist_comm]; exact hbd⟩
    have hzd : dist z y < min ε₁ ε₂ := hzb
    rcases lt_or_gt_of_ne hz with hn | hp
    · right
      refine ⟨hcD, z, hzc, hb₂ z (lt_of_lt_of_le hzd (min_le_right _ _)) ?_⟩
      show 0 < sda (O.dK h₀) (O.dL h₀) z
      rw [sda_swap_left]; linarith
    · left
      exact ⟨hcD, z, hzc, hb₁ z (lt_of_lt_of_le hzd (min_le_left _ _)) hp⟩

/-- An outer side keeps its two end points: the side `σ m` of the new layout is the segment of
the side `m` of the old one. -/
theorem outer_lay {m : Halfedge F} (hm : Outer t h₀ m) :
    O.flipLay h₀ (t.G.flipPerm h₀ m).1 (t.G.flipPerm h₀ m).2 = O.lay m.1 m.2 ∧
      O.flipLay h₀ (t.G.flipPerm h₀ m).1 ((t.G.flipPerm h₀ m).2 + 1) = O.lay m.1 (m.2 + 1) := by
  have hf := hF.flippable
  rcases hm with rfl | rfl | rfl | rfl
  · rw [Gluing.flipPerm_next hf, flipLay_f₁ hf, rot_self, rot_add_one]
    exact ⟨rfl, congrArg (O.lay h₀.1) (show h₀.2 + 2 = h₀.2 + 1 + 1 by omega)⟩
  · rw [Gluing.flipPerm_of_ne (Halfedge.prev_ne _) (Halfedge.prev_ne_next _)
      (hf.ne h₀.prev (t.G.twin h₀)) (hf.ne h₀.prev (t.G.twin h₀).next)]
    show O.flipLay h₀ h₀.1 (h₀.2 + 2) = O.lay h₀.1 (h₀.2 + 2) ∧
      O.flipLay h₀ h₀.1 (h₀.2 + 2 + 1) = O.lay h₀.1 (h₀.2 + 2 + 1)
    rw [flipLay_f₀, rot_add_two, show h₀.2 + 2 + 1 = h₀.2 by omega, rot_self]
    exact ⟨rfl, rfl⟩
  · rw [Gluing.flipPerm_twin_next hf, flipLay_f₀, rot_self, rot_add_one]
    exact ⟨hF.lay_twin_head.symm, congrArg (O.lay (t.G.twin h₀).1)
      (show (t.G.twin h₀).2 + 2 = (t.G.twin h₀).2 + 1 + 1 by omega)⟩
  · rw [Gluing.flipPerm_of_ne (hf.ne' (t.G.twin h₀).prev h₀) (hf.ne' (t.G.twin h₀).prev h₀.next)
      (Halfedge.prev_ne _) (Halfedge.prev_ne_next _)]
    show O.flipLay h₀ (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2) =
        O.lay (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2) ∧
      O.flipLay h₀ (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2 + 1) =
        O.lay (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2 + 1)
    rw [flipLay_f₁ hf, rot_add_two, show (t.G.twin h₀).2 + 2 + 1 = (t.G.twin h₀).2 by omega,
      rot_self]
    exact ⟨rfl, hF.lay_twin_tail.symm⟩

/-- The cells of the new face at a point of an outer side are the cells of the old face. -/
theorem outer_germAt {m : Halfedge F} (hm : Outer t h₀ m) {y : Pt}
    (hy : y ∈ (O.lay m.1).side m.2) :
    germAt s (O.flipCells h₀ (t.G.flipPerm h₀ m).1) y = germAt s (O.cells m.1) y := by
  have hf := hF.flippable
  rcases hm with rfl | rfl | rfl | rfl
  · rw [Gluing.flipPerm_next hf]
    have hy' : y ∈ openSegment ℝ (O.dP (t.G.twin h₀)) (O.dL (t.G.twin h₀)) := by
      rw [hF.dP_twin, hF.dL_twin]
      have : y ∈ openSegment ℝ (O.lay h₀.1 (h₀.2 + 1)) (O.lay h₀.1 (h₀.2 + 1 + 1)) := hy
      rwa [show h₀.2 + 1 + 1 = h₀.2 + 2 by omega] at this
    have := hF.symm.germAt_flipCells_f₀_of_mem_PL hy'
    rw [hF.flipCells_twin, t.G.twin_twin] at this
    exact this
  · rw [Gluing.flipPerm_of_ne (Halfedge.prev_ne _) (Halfedge.prev_ne_next _)
      (hf.ne h₀.prev (t.G.twin h₀)) (hf.ne h₀.prev (t.G.twin h₀).next)]
    have hy' : y ∈ openSegment ℝ (O.dK h₀) (O.dP h₀) := by
      have : y ∈ openSegment ℝ (O.lay h₀.1 (h₀.2 + 2)) (O.lay h₀.1 (h₀.2 + 2 + 1)) := hy
      rwa [show h₀.2 + 2 + 1 = h₀.2 by omega] at this
    exact hF.germAt_flipCells_f₀_of_mem_KP hy'
  · rw [Gluing.flipPerm_twin_next hf]
    have hy' : y ∈ openSegment ℝ (O.dP h₀) (O.dL h₀) := by
      have : y ∈ openSegment ℝ (O.lay (t.G.twin h₀).1 ((t.G.twin h₀).2 + 1))
          (O.lay (t.G.twin h₀).1 ((t.G.twin h₀).2 + 1 + 1)) := hy
      rwa [hF.lay_twin_head, show (t.G.twin h₀).2 + 1 + 1 = (t.G.twin h₀).2 + 2 by omega] at this
    exact hF.germAt_flipCells_f₀_of_mem_PL hy'
  · rw [Gluing.flipPerm_of_ne (hf.ne' (t.G.twin h₀).prev h₀) (hf.ne' (t.G.twin h₀).prev h₀.next)
      (Halfedge.prev_ne _) (Halfedge.prev_ne_next _)]
    have hy' : y ∈ openSegment ℝ (O.dK (t.G.twin h₀)) (O.dP (t.G.twin h₀)) := by
      rw [hF.dK_twin, hF.dP_twin]
      have : y ∈ openSegment ℝ (O.lay (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2))
          (O.lay (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2 + 1)) := hy
      rwa [show (t.G.twin h₀).2 + 2 + 1 = (t.G.twin h₀).2 by omega, hF.lay_twin_tail] at this
    have := hF.symm.germAt_flipCells_f₀_of_mem_KP hy'
    rw [hF.flipCells_twin] at this
    exact this

omit hF in
/-- Every halfedge is `h₀`, `twin h₀`, outer or far. -/
theorem outer_or_far (m : Halfedge F) :
    m = h₀ ∨ m = t.G.twin h₀ ∨ Outer t h₀ m ∨ Far t h₀ m := by
  by_cases h0 : m.1 = h₀.1
  · rcases Halfedge.eq_or_eq_next_or_eq_prev h0 with rfl | rfl | rfl
    · exact Or.inl rfl
    · exact Or.inr (Or.inr (Or.inl (Or.inl rfl)))
    · exact Or.inr (Or.inr (Or.inl (Or.inr (Or.inl rfl))))
  by_cases h1 : m.1 = (t.G.twin h₀).1
  · rcases Halfedge.eq_or_eq_next_or_eq_prev h1 with rfl | rfl | rfl
    · exact Or.inr (Or.inl rfl)
    · exact Or.inr (Or.inr (Or.inl (Or.inr (Or.inr (Or.inl rfl)))))
    · exact Or.inr (Or.inr (Or.inl (Or.inr (Or.inr (Or.inr rfl)))))
  exact Or.inr (Or.inr (Or.inr ⟨h0, h1⟩))

/-- The twin of an outer side is outer or far (never on the flipped edge). -/
theorem twin_outer {m : Halfedge F} (hm : Outer t h₀ m) :
    Outer t h₀ (t.G.twin m) ∨ Far t h₀ (t.G.twin m) := by
  have hf := hF.flippable
  have h1 : t.G.twin m ≠ h₀ := by
    intro e
    rw [t.G.twin_eq_iff] at e
    rcases hm with h | h | h | h <;> rw [h] at e
    · exact (hf.ne h₀.next (t.G.twin h₀)) e
    · exact (hf.ne h₀.prev (t.G.twin h₀)) e
    · exact Halfedge.next_ne _ e
    · exact Halfedge.prev_ne _ e
  have h2 : t.G.twin m ≠ t.G.twin h₀ := by
    intro e
    have e' := t.G.twin_injective e
    rcases hm with h | h | h | h <;> rw [h] at e'
    · exact Halfedge.next_ne _ e'
    · exact Halfedge.prev_ne _ e'
    · exact (hf.ne' (t.G.twin h₀).next h₀) e'
    · exact (hf.ne' (t.G.twin h₀).prev h₀) e'
  rcases outer_or_far (t := t) (h₀ := h₀) (t.G.twin m) with e | e | h | h
  · exact absurd e h1
  · exact absurd e h2
  · exact Or.inl h
  · exact Or.inr h

theorem twin_far {m : Halfedge F} (hm : Far t h₀ m) :
    Outer t h₀ (t.G.twin m) ∨ Far t h₀ (t.G.twin m) := by
  have h1 : t.G.twin m ≠ h₀ := fun e => hm.2 (by rw [t.G.twin_eq_iff.1 e])
  have h2 : t.G.twin m ≠ t.G.twin h₀ := fun e => hm.1 (by rw [t.G.twin_injective e])
  rcases outer_or_far (t := t) (h₀ := h₀) (t.G.twin m) with e | e | h | h
  · exact absurd e h1
  · exact absurd e h2
  · exact Or.inl h
  · exact Or.inr h

theorem twin_outer_or_far {m : Halfedge F} (hm : Outer t h₀ m ∨ Far t h₀ m) :
    Outer t h₀ (t.G.twin m) ∨ Far t h₀ (t.G.twin m) :=
  hm.elim hF.twin_outer hF.twin_far

theorem lay_of_outer_or_far {m : Halfedge F} (hm : Outer t h₀ m ∨ Far t h₀ m) :
    O.flipLay h₀ (t.G.flipPerm h₀ m).1 (t.G.flipPerm h₀ m).2 = O.lay m.1 m.2 ∧
      O.flipLay h₀ (t.G.flipPerm h₀ m).1 ((t.G.flipPerm h₀ m).2 + 1) = O.lay m.1 (m.2 + 1) := by
  rcases hm with hm | hm
  · exact hF.outer_lay hm
  · rw [far_flipPerm hm, far_flipLay hm]; exact ⟨rfl, rfl⟩

theorem side_of_outer_or_far {m : Halfedge F} (hm : Outer t h₀ m ∨ Far t h₀ m) :
    (O.flipLay h₀ (t.G.flipPerm h₀ m).1).side (t.G.flipPerm h₀ m).2 = (O.lay m.1).side m.2 := by
  obtain ⟨e1, e2⟩ := hF.lay_of_outer_or_far hm
  unfold Tri.side
  rw [e1, e2]

theorem germAt_of_outer_or_far {m : Halfedge F} (hm : Outer t h₀ m ∨ Far t h₀ m) {y : Pt}
    (hy : y ∈ (O.lay m.1).side m.2) :
    germAt s (O.flipCells h₀ (t.G.flipPerm h₀ m).1) y = germAt s (O.cells m.1) y := by
  rcases hm with hm | hm
  · exact hF.outer_germAt hm hy
  · rw [far_flipPerm hm, far_flipCells hm]

/-- The length of the new edge is the distance of the laid out `K` and `L`. -/
theorem newLength_eq : t.G.newLength t.ℓ h₀ = dist2 (O.dK h₀) (O.dL h₀) := by
  have hq := hF.quad
  obtain ⟨hk, hl⟩ := t.isMetric.diamond_isTriangle h₀
  have e1 : t.ℓ h₀ = dist2 (O.dP h₀) (O.dQ h₀) := (O.lay_dist h₀.1 h₀.2).symm
  have e2 : t.ℓ h₀.prev = dist2 (O.dP h₀) (O.dK h₀) := by
    have := O.lay_dist h₀.1 (h₀.2 + 2)
    rw [show h₀.2 + 2 + 1 = h₀.2 by omega, dist2_comm] at this
    exact this.symm
  have e3 : t.ℓ h₀.next = dist2 (O.dQ h₀) (O.dK h₀) := by
    have := O.lay_dist h₀.1 (h₀.2 + 1)
    rw [show h₀.2 + 1 + 1 = h₀.2 + 2 by omega] at this
    exact this.symm
  have e4 : t.ℓ (t.G.twin h₀).next = dist2 (O.dP h₀) (O.dL h₀) := by
    have := O.lay_dist (t.G.twin h₀).1 ((t.G.twin h₀).2 + 1)
    rw [show (t.G.twin h₀).2 + 1 + 1 = (t.G.twin h₀).2 + 2 by omega, hF.lay_twin_head] at this
    exact this.symm
  have e5 : t.ℓ (t.G.twin h₀).prev = dist2 (O.dQ h₀) (O.dL h₀) := by
    have := O.lay_dist (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2)
    rw [show (t.G.twin h₀).2 + 2 + 1 = (t.G.twin h₀).2 by omega, hF.lay_twin_tail,
      dist2_comm] at this
    exact this.symm
  unfold Gluing.newLength
  rw [flipLength_eq_flippedLengthFromLengths hk hl, e1, e2, e3, e4, e5]
  exact DdgProofs.Intrinsic.flippedLengthFromLengths_realized hq.hK hq.hL'

/-! ### The flipped state -/

variable {t' : MetricGluing F} (hG : t'.G = t.G.flip h₀) (hℓ : t'.ℓ = t.G.flipLen t.ℓ h₀)
include hG

theorem twin_flip (m : Halfedge F) :
    t'.G.twin (t.G.flipPerm h₀ m) = t.G.flipPerm h₀ (t.G.twin m) := by
  rw [hG]; exact Gluing.flip_twin_apply _ _ _

/-- The matching motion of a side of `t'` that is not the new edge is the old one. -/
theorem layΔ_flip_eq {m : Halfedge F}
    (hm : Outer t h₀ m ∨ Far t h₀ m) :
    layΔ t' (O.flipLay h₀) (t.G.flipPerm h₀ m) = layΔ t O.lay m := by
  obtain ⟨a1, a2⟩ := hF.lay_of_outer_or_far (hF.twin_outer_or_far hm)
  obtain ⟨b1, b2⟩ := hF.lay_of_outer_or_far hm
  unfold layΔ
  rw [twin_flip hF hG m, a1, a2, b1, b2]

theorem twin_flip_next : t'.G.twin h₀.next = (t.G.twin h₀).next := by
  rw [hG]; exact Gluing.flip_twin_next hF.flippable

/-- The matching motion of the new edge is the identity. -/
theorem layΔ_flip_next : layΔ t' (O.flipLay h₀) h₀.next = Rigid.id := by
  unfold layΔ
  rw [twin_flip_next hF hG]
  show rigidOfSegment (O.flipLay h₀ (t.G.twin h₀).1 ((t.G.twin h₀).2 + 1))
      (O.flipLay h₀ (t.G.twin h₀).1 ((t.G.twin h₀).2 + 1 + 1))
      (O.flipLay h₀ h₀.1 (h₀.2 + 1 + 1)) (O.flipLay h₀ h₀.1 (h₀.2 + 1)) = Rigid.id
  rw [flipLay_f₀, flipLay_f₁ hF.flippable,
    show (t.G.twin h₀).2 + 1 + 1 = (t.G.twin h₀).2 + 2 by omega,
    show h₀.2 + 1 + 1 = h₀.2 + 2 by omega]
  simp only [rot_add_one, rot_add_two]
  exact rigidOfSegment_self hF.quad.ne_LK.symm

/-- The germ along an old side of `t'` (not the new edge) is the old germ. -/
theorem germ_edge_flip_of_outer_or_far {m : Halfedge F}
    (hm : Outer t h₀ m ∨ Far t h₀ m)
    (hb : t'.G.twin (t.G.flipPerm h₀ m) ≠ t.G.flipPerm h₀ m) {y : Pt}
    (hy : y ∈ (O.flipLay h₀ (t.G.flipPerm h₀ m).1).side (t.G.flipPerm h₀ m).2) :
    ConsistentGerm s (germAt s (O.flipCells h₀ (t.G.flipPerm h₀ m).1 ∪
      (O.flipCells h₀ (t'.G.twin (t.G.flipPerm h₀ m)).1).image
        (Cell.mv (layΔ t' (O.flipLay h₀) (t.G.flipPerm h₀ m)))) y) y := by
  have hb' : t.G.twin m ≠ m := by
    rw [twin_flip hF hG m] at hb
    exact fun e => hb (by rw [e])
  have hy' : y ∈ (O.lay m.1).side m.2 := by rw [← hF.side_of_outer_or_far hm]; exact hy
  have hm' := hF.twin_outer_or_far hm
  have hy2 : (layΔ t O.lay m).inv y ∈ (O.lay (t.G.twin m).1).side (t.G.twin m).2 := by
    rw [← Tri.mem_map_side_iff (layΔ t O.lay m) (O.lay_pos (t.G.twin m).1),
      Rigid.apply_inv_apply, O.side_layΔ m]
    exact hy'
  rw [twin_flip hF hG m, layΔ_flip_eq hF hG hm, germAt_union_image,
    hF.germAt_of_outer_or_far hm hy', hF.germAt_of_outer_or_far hm' hy2,
    ← germAt_union_image]
  exact O.germ_edge m hb' y hy'

/-- The germ along the new edge. -/
theorem germ_edge_flip_next {y : Pt} (hy : y ∈ (O.flipLay h₀ h₀.1).side (h₀.2 + 1)) :
    ConsistentGerm s (germAt s (O.flipCells h₀ h₀.1 ∪
      (O.flipCells h₀ (t'.G.twin h₀.next).1).image
        (Cell.mv (layΔ t' (O.flipLay h₀) h₀.next))) y) y := by
  have hmv : (Cell.mv Rigid.id : Cell F → Cell F) = id := funext Cell.mv_id
  rw [layΔ_flip_next hF hG, twin_flip_next hF hG, hmv, Finset.image_id, Halfedge.next_fst]
  have hy' : y ∈ openSegment ℝ (O.dL h₀) (O.dK h₀) := by
    have : y ∈ openSegment ℝ (O.flipLay h₀ h₀.1 (h₀.2 + 1)) (O.flipLay h₀ h₀.1 (h₀.2 + 1 + 1)) :=
      hy
    rwa [flipLay_f₀, show h₀.2 + 1 + 1 = h₀.2 + 2 by omega, rot_add_one, rot_add_two] at this
  rw [hF.germAt_flipCells_union_of_mem_LK hy']
  exact hF.germ_D (hF.quad.openSegment_new_subset hy')

theorem germ_bdry_flip_of_outer_or_far {m : Halfedge F}
    (hm : Outer t h₀ m ∨ Far t h₀ m)
    (hb : t'.G.twin (t.G.flipPerm h₀ m) = t.G.flipPerm h₀ m) {y : Pt}
    (hy : y ∈ (O.flipLay h₀ (t.G.flipPerm h₀ m).1).side (t.G.flipPerm h₀ m).2) :
    BoundaryGerm s (germAt s (O.flipCells h₀ (t.G.flipPerm h₀ m).1) y) y := by
  have hb' : t.G.twin m = m := by
    rw [twin_flip hF hG m] at hb
    exact (t.G.flipPerm h₀).injective hb
  have hy' : y ∈ (O.lay m.1).side m.2 := by rw [← hF.side_of_outer_or_far hm]; exact hy
  rw [hF.germAt_of_outer_or_far hm hy']
  exact O.germ_bdry m hb' y hy'

include hℓ in
/-- The new layouts have the new lengths. -/
theorem flipLay_dist (f : F) (k : Fin 3) :
    dist2 (O.flipLay h₀ f k) (O.flipLay h₀ f (k + 1)) = t'.ℓ (f, k) := by
  have hf := hF.flippable
  rw [hℓ]
  rcases eq_or_ne f h₀.1 with rfl | h0
  · rw [flipLay_f₀]
    rcases (by omega : k = h₀.2 ∨ k = h₀.2 + 1 ∨ k = h₀.2 + 2) with rfl | rfl | rfl
    · rw [rot_self, rot_add_one, show ((h₀.1, h₀.2) : Halfedge F) = h₀ from rfl,
        Gluing.flipLen_self hf]
      have := O.lay_dist (t.G.twin h₀).1 ((t.G.twin h₀).2 + 1)
      rw [show (t.G.twin h₀).2 + 1 + 1 = (t.G.twin h₀).2 + 2 by omega, hF.lay_twin_head] at this
      exact this
    · rw [show h₀.2 + 1 + 1 = h₀.2 + 2 by omega, rot_add_one, rot_add_two,
        show ((h₀.1, h₀.2 + 1) : Halfedge F) = h₀.next from rfl, Gluing.flipLen_next,
        hF.newLength_eq, dist2_comm]
    · rw [show h₀.2 + 2 + 1 = h₀.2 by omega, rot_add_two, rot_self,
        show ((h₀.1, h₀.2 + 2) : Halfedge F) = h₀.prev from rfl, Gluing.flipLen_prev hf]
      have := O.lay_dist h₀.1 (h₀.2 + 2)
      rw [show h₀.2 + 2 + 1 = h₀.2 by omega] at this
      exact this
  rcases eq_or_ne f (t.G.twin h₀).1 with rfl | h1
  · rw [flipLay_f₁ hf]
    rcases (by omega : k = (t.G.twin h₀).2 ∨ k = (t.G.twin h₀).2 + 1 ∨ k = (t.G.twin h₀).2 + 2)
      with rfl | rfl | rfl
    · rw [rot_self, rot_add_one,
        show (((t.G.twin h₀).1, (t.G.twin h₀).2) : Halfedge F) = t.G.twin h₀ from rfl,
        Gluing.flipLen_twin hf]
      have := O.lay_dist h₀.1 (h₀.2 + 1)
      rw [show h₀.2 + 1 + 1 = h₀.2 + 2 by omega] at this
      exact this
    · rw [show (t.G.twin h₀).2 + 1 + 1 = (t.G.twin h₀).2 + 2 by omega, rot_add_one, rot_add_two,
        show (((t.G.twin h₀).1, (t.G.twin h₀).2 + 1) : Halfedge F) = (t.G.twin h₀).next from rfl,
        Gluing.flipLen_twin_next, hF.newLength_eq]
    · rw [show (t.G.twin h₀).2 + 2 + 1 = (t.G.twin h₀).2 by omega, rot_add_two, rot_self,
        show (((t.G.twin h₀).1, (t.G.twin h₀).2 + 2) : Halfedge F) = (t.G.twin h₀).prev from rfl,
        Gluing.flipLen_twin_prev hf]
      have := O.lay_dist (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2)
      rw [show (t.G.twin h₀).2 + 2 + 1 = (t.G.twin h₀).2 by omega, hF.lay_twin_tail] at this
      exact this
  rw [flipLay_other h0 h1, Gluing.flipLen_of_face hf t.ℓ (x := (f, k)) h0 h1]
  exact O.lay_dist f k

omit hG in
theorem flipLay_pos (f : F) : (O.flipLay h₀ f).Pos := by
  by_cases h0 : f = h₀.1
  · subst h0; rw [flipLay_f₀]; exact rot_pos hF.quad.hc1
  by_cases h1 : f = (t.G.twin h₀).1
  · subst h1; rw [flipLay_f₁ hF.flippable]; exact rot_pos hF.quad.hc2
  rw [flipLay_other h0 h1]; exact O.lay_pos f

include hℓ in
/-- **The flipped overlay.** -/
def flipOverlay : Overlay s t' where
  lay := O.flipLay h₀
  lay_dist := flipLay_dist hF hG hℓ
  lay_pos := flipLay_pos hF
  cells := O.flipCells h₀
  cover := by
    intro f
    by_cases h0 : f = h₀.1
    · subst h0; exact hF.cover_f₀
    by_cases h1 : f = (t.G.twin h₀).1
    · subst h1
      have := hF.symm.cover_f₀
      rwa [hF.flipLay_twin, hF.flipCells_twin] at this
    rw [flipLay_other h0 h1, flipCells_other h0 h1]; exact O.cover f
  nonjunk := by
    intro f
    by_cases h0 : f = h₀.1
    · subst h0; exact hF.nonjunk_f₀
    by_cases h1 : f = (t.G.twin h₀).1
    · subst h1
      have := hF.symm.nonjunk_f₀
      rwa [hF.flipLay_twin, hF.flipCells_twin] at this
    rw [flipLay_other h0 h1, flipCells_other h0 h1]; exact O.nonjunk f
  germ_int := by
    intro f
    by_cases h0 : f = h₀.1
    · subst h0; exact hF.germ_int_f₀
    by_cases h1 : f = (t.G.twin h₀).1
    · subst h1
      have := hF.symm.germ_int_f₀
      rwa [hF.flipLay_twin, hF.flipCells_twin] at this
    rw [flipLay_other h0 h1, flipCells_other h0 h1]; exact O.germ_int f
  germ_edge := by
    intro x hb y hy
    obtain ⟨m, rfl⟩ := (t.G.flipPerm h₀).surjective x
    rcases outer_or_far (t := t) (h₀ := h₀) m with rfl | rfl | hm | hm
    · rw [Gluing.flipPerm_self hF.flippable] at hy ⊢
      exact hF.germ_edge_flip_next hG hy
    · rw [Gluing.flipPerm_twin hF.flippable] at hy ⊢
      have hG' : t'.G = t.G.flip (t.G.twin h₀) := by
        rw [hG, Surface.Gluing.flip_twin_eq hF.flippable]
      have hy' : y ∈ (O.flipLay (t.G.twin h₀) (t.G.twin h₀).1).side ((t.G.twin h₀).2 + 1) := by
        rw [hF.flipLay_twin]; exact hy
      have := hF.symm.germ_edge_flip_next hG' hy'
      rwa [hF.flipLay_twin, hF.flipCells_twin] at this
    · exact hF.germ_edge_flip_of_outer_or_far hG (Or.inl hm) hb hy
    · exact hF.germ_edge_flip_of_outer_or_far hG (Or.inr hm) hb hy
  germ_bdry := by
    intro x hb y hy
    obtain ⟨m, rfl⟩ := (t.G.flipPerm h₀).surjective x
    rcases outer_or_far (t := t) (h₀ := h₀) m with e | e | hm | hm
    · exfalso
      rw [e, twin_flip hF hG h₀, Gluing.flipPerm_twin hF.flippable,
        Gluing.flipPerm_self hF.flippable] at hb
      exact (hF.flippable.ne' (t.G.twin h₀).next h₀.next) hb
    · exfalso
      rw [e, twin_flip hF hG (t.G.twin h₀), t.G.twin_twin, Gluing.flipPerm_self hF.flippable,
        Gluing.flipPerm_twin hF.flippable] at hb
      exact (hF.flippable.ne h₀.next (t.G.twin h₀).next) hb
    · exact hF.germ_bdry_flip_of_outer_or_far hG (Or.inl hm) hb hy
    · exact hF.germ_bdry_flip_of_outer_or_far hG (Or.inr hm) hb hy
  corner := by
    intro f
    by_cases h0 : f = h₀.1
    · subst h0; exact hF.corner_f₀
    by_cases h1 : f = (t.G.twin h₀).1
    · subst h1
      have := hF.symm.corner_f₀
      rwa [hF.flipLay_twin, hF.flipCells_twin] at this
    rw [flipLay_other h0 h1, flipCells_other h0 h1]; exact O.corner f

end FlipOK

/-- **A convex flip keeps an overlay**: if `t` is overlaid on `s` and `t'` is `t` with the edge
`h₀` flipped (a strictly convex diamond), then `t'` is overlaid on `s`. -/
theorem flip (O : Overlay s t) {t' : MetricGluing F} (h : MetricGluing.ConvexFlip t t') :
    Nonempty (Overlay s t') := by
  obtain ⟨h₀, hf, hconv, hG, hℓ⟩ := h
  let O₁ := O.relayout (t.G.twin h₀).1 (layΔ t O.lay h₀)
  have hF : FlipOK O₁ h₀ := ⟨hf, hconv, O.relayout_layΔ_eq_id hf⟩
  exact ⟨hF.flipOverlay hG hℓ⟩

end Overlay

end DdgProofs.Delaunay.Surface
