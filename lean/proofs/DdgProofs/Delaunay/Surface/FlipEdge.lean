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
open scoped Classical

variable {F : Type*}

/-! ## The flip of the twin halfedge is the same flip -/

namespace Gluing

variable [DecidableEq F] {G : Gluing F} {h₀ : Halfedge F}

theorem flipPerm_twin_eq (hf : G.Flippable h₀) : G.flipPerm (G.twin h₀) = G.flipPerm h₀ := by
  sorry

theorem flip_twin_eq (hf : G.Flippable h₀) : G.flip (G.twin h₀) = G.flip h₀ := by
  sorry

theorem newLength_twin {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (hf : G.Flippable h₀) :
    G.newLength ℓ (G.twin h₀) = G.newLength ℓ h₀ := by
  sorry

theorem flipLen_twin_eq {ℓ : Halfedge F → ℝ} (hm : G.IsMetric ℓ) (hf : G.Flippable h₀) :
    G.flipLen ℓ (G.twin h₀) = G.flipLen ℓ h₀ := by
  sorry

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

namespace FlipOK

variable (hF : FlipOK O h₀)
include hF

/-! ### Germs on the sides of the new face `f₀` -/

/-! ### Outer sides keep their end points and their germs -/

/-- On the open side `P L` of the new face `f₀` (the old side `(twin h₀).next` of `f₁`), the
cells of the new face are the old cells of `f₁`. -/
theorem germAt_flipCells_f₀_of_mem_PL {y : Pt} (hy : y ∈ openSegment ℝ (O.dP h₀) (O.dL h₀)) :
    germAt s (O.flipCells h₀ h₀.1) y = germAt s (O.cells (t.G.twin h₀).1) y := by
  sorry

/-- On the open side `K P` of the new face `f₀` (the old side `h₀.prev` of `f₀`), the cells of
the new face are the old cells of `f₀`. -/
theorem germAt_flipCells_f₀_of_mem_KP {y : Pt} (hy : y ∈ openSegment ℝ (O.dK h₀) (O.dP h₀)) :
    germAt s (O.flipCells h₀ h₀.1) y = germAt s (O.cells h₀.1) y := by
  sorry

/-- On the new open diagonal `L K`, the cells of the two new faces together are the cells of
the diamond. -/
theorem germAt_flipCells_union_of_mem_LK {y : Pt} (hy : y ∈ openSegment ℝ (O.dL h₀) (O.dK h₀)) :
    germAt s (O.flipCells h₀ h₀.1 ∪ O.flipCells h₀ (t.G.twin h₀).1) y =
      germAt s (O.Dcells h₀) y := by
  sorry

/-- An outer side keeps its two end points: the side `σ m` of the new layout is the segment of
the side `m` of the old one. -/
theorem outer_lay {m : Halfedge F} (hm : Outer t h₀ m) :
    O.flipLay h₀ (t.G.flipPerm h₀ m).1 (t.G.flipPerm h₀ m).2 = O.lay m.1 m.2 ∧
      O.flipLay h₀ (t.G.flipPerm h₀ m).1 ((t.G.flipPerm h₀ m).2 + 1) = O.lay m.1 (m.2 + 1) := by
  sorry

/-- The cells of the new face at a point of an outer side are the cells of the old face. -/
theorem outer_germAt {m : Halfedge F} (hm : Outer t h₀ m) {y : Pt}
    (hy : y ∈ (O.lay m.1).side m.2) :
    germAt s (O.flipCells h₀ (t.G.flipPerm h₀ m).1) y = germAt s (O.cells m.1) y := by
  sorry

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
  sorry

theorem twin_far {m : Halfedge F} (hm : Far t h₀ m) :
    Outer t h₀ (t.G.twin m) ∨ Far t h₀ (t.G.twin m) := by
  sorry

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
  sorry

/-- The matching motion of the new edge is the identity. -/
theorem layΔ_flip_next : layΔ t' (O.flipLay h₀) h₀.next = Rigid.id := by
  sorry

/-- The germ along an old side of `t'` (not the new edge) is the old germ. -/
theorem germ_edge_flip_of_outer_or_far {m : Halfedge F}
    (hm : Outer t h₀ m ∨ Far t h₀ m)
    (hb : t'.G.twin (t.G.flipPerm h₀ m) ≠ t.G.flipPerm h₀ m) {y : Pt}
    (hy : y ∈ (O.flipLay h₀ (t.G.flipPerm h₀ m).1).side (t.G.flipPerm h₀ m).2) :
    ConsistentGerm s (germAt s (O.flipCells h₀ (t.G.flipPerm h₀ m).1 ∪
      (O.flipCells h₀ (t'.G.twin (t.G.flipPerm h₀ m)).1).image
        (Cell.mv (layΔ t' (O.flipLay h₀) (t.G.flipPerm h₀ m)))) y) y := by
  sorry

/-- The germ along the new edge. -/
theorem germ_edge_flip_next {y : Pt} (hy : y ∈ (O.flipLay h₀ h₀.1).side (h₀.2 + 1)) :
    ConsistentGerm s (germAt s (O.flipCells h₀ h₀.1 ∪
      (O.flipCells h₀ (t'.G.twin h₀.next).1).image
        (Cell.mv (layΔ t' (O.flipLay h₀) h₀.next))) y) y := by
  sorry

theorem germ_bdry_flip_of_outer_or_far {m : Halfedge F}
    (hm : Outer t h₀ m ∨ Far t h₀ m)
    (hb : t'.G.twin (t.G.flipPerm h₀ m) = t.G.flipPerm h₀ m) {y : Pt}
    (hy : y ∈ (O.flipLay h₀ (t.G.flipPerm h₀ m).1).side (t.G.flipPerm h₀ m).2) :
    BoundaryGerm s (germAt s (O.flipCells h₀ (t.G.flipPerm h₀ m).1) y) y := by
  sorry

include hℓ in
/-- The new layouts have the new lengths. -/
theorem flipLay_dist (f : F) (k : Fin 3) :
    dist2 (O.flipLay h₀ f k) (O.flipLay h₀ f (k + 1)) = t'.ℓ (f, k) := by
  sorry

omit hG in
theorem flipLay_pos (f : F) : (O.flipLay h₀ f).Pos := by
  sorry

include hℓ in
/-- **The flipped overlay.** -/
def flipOverlay : Overlay s t' where
  lay := O.flipLay h₀
  lay_dist := flipLay_dist hF hG hℓ
  lay_pos := flipLay_pos hF
  cells := O.flipCells h₀
  cover := by sorry
  nonjunk := by sorry
  germ_int := by sorry
  germ_edge := by sorry
  germ_bdry := by sorry
  corner := by sorry

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
