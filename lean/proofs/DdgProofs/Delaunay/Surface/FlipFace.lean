import DdgProofs.Delaunay.Surface.Overlay
import DdgProofs.Delaunay.Surface.Quad

/-!
# Delaunay layer 4: the flip of an overlaid triangulation, face level

The flip of the edge `h₀ = (f₀, a)` of `t` (`DdgProofs/Delaunay/Flip.lean`, `Metric.lean`) on an
overlay `O : Overlay s t` whose two faces `f₀ = (i, j, k)` and `f₁ = (j, i, l)` are laid out
**in one plane** (`FlipOK.coplanar`: the motion `layΔ` matching the two copies of the edge `ij`
is the identity, so `lay f₁` is the triangle `j, i, l` next to `lay f₀ = i, j, k`). With the
diamond strictly convex (`FlipOK.convex`), the four points `P = i, Q = j, K = k, L = l` form a
`ConvexQuad` (`Quad.lean`).

The flipped overlay: the faces `f₀, f₁` get the layouts `(i, l, k)` and `(j, k, l)` (`flipLay`,
with the corner numbering of the Rust flip: `i→l` keeps the slot `a`, `j→k` the slot `b`), and
the cells of the diamond `Dcells = cells f₀ ∪ cells f₁` are distributed to the two new faces
by overlap (`flipCells`). Nothing else changes.

This file proves the face-level invariants of the new face `f₀`: the cells of the diamond that
reach a point of `f₀`'s side of the old diagonal are cells of `f₀` (`mem_cells_of_sda_pos`, the
key lemma, by crossing the old diagonal at an interior point of the cell and the edge germ
there), the germs at points of the open old diamond are consistent (`germ_D`), and `cover`,
`nonjunk`, `corner`, `germ_int` for the new face. The symmetric statements for `f₁` follow from
`FlipOK.symm` (the flip of the twin halfedge is the same flip) in `FlipEdge.lean`.
-/

noncomputable section

namespace DDG.Delaunay.Surface

open DDG.Intrinsic (dist2_comm)
open DDG.Delaunay.Planar
open Set
open scoped Classical

variable {F : Type*}

/-- The triangle with corner `a ↦ A`, `a + 1 ↦ B`, `a + 2 ↦ C`. -/
def rot (a : Fin 3) (A B C : Pt) : Tri := fun m => if m = a then A else if m = a + 1 then B else C

@[simp] theorem rot_self (a : Fin 3) (A B C : Pt) : rot a A B C a = A := by simp [rot]
@[simp] theorem rot_add_one (a : Fin 3) (A B C : Pt) : rot a A B C (a + 1) = B := by
  simp [rot, show a + 1 ≠ a by omega]
@[simp] theorem rot_add_two (a : Fin 3) (A B C : Pt) : rot a A B C (a + 2) = C := by
  simp [rot, show a + 2 ≠ a by omega, show a + 2 ≠ a + 1 by omega]

theorem rot_eq_self (T : Tri) (a : Fin 3) : rot a (T a) (T (a + 1)) (T (a + 2)) = T := by
  funext m
  rcases (by omega : m = a ∨ m = a + 1 ∨ m = a + 2) with rfl | rfl | rfl <;> simp

theorem rot_hull (a : Fin 3) (A B C : Pt) : (rot a A B C).hull = (tri3 A B C).hull := by
  rw [Tri.hull_rot (rot a A B C) a]; simp

theorem rot_int (a : Fin 3) (A B C : Pt) : (rot a A B C).int = (tri3 A B C).int := by
  rw [Tri.int, Tri.int, rot_hull]

theorem rot_pos {a : Fin 3} {A B C : Pt} (h : 0 < sda A B C) : (rot a A B C).Pos := by
  unfold Tri.Pos
  fin_cases a <;> simp [rot] <;>
    first | exact h | (rw [sda_rotate]; exact h) | (rw [← sda_rotate]; exact h)

theorem sda_of_rot_pos {a : Fin 3} {A B C : Pt} (h : (rot a A B C).Pos) : 0 < sda A B C := by
  unfold Tri.Pos at h
  fin_cases a <;> simp [rot] at h <;> linarith [sda_rotate A B C, sda_rotate C A B]

theorem cornerAngle_swap' (a b c : ℝ) : cornerAngle a b c = cornerAngle a c b := by
  unfold cornerAngle
  congr 1
  ring

/-- A consistent germ containing a cell with the point interior is that single cell. -/
theorem ConsistentGerm.eq_singleton_of_mem_int {s : MetricGluing F} {Γ : Finset (Cell F)}
    {y : Pt} (h : ConsistentGerm s Γ y) {c : Cell F} (hc : c ∈ Γ) (hy : y ∈ (c.tri s).int) :
    Γ = {c} := by
  rcases h with ⟨c₁, rfl, -⟩ | ⟨c₁, k₁, -, rfl, hy₁⟩
  · rw [Finset.mem_singleton] at hc; rw [hc]
  · exfalso
    rw [Finset.mem_insert, Finset.mem_singleton] at hc
    rcases hc with rfl | rfl
    · exact Tri.not_mem_int_of_mem_side (Cell.tri_pos _) hy₁ hy
    · rw [← Cell.across_side] at hy₁
      exact Tri.not_mem_int_of_mem_side (Cell.tri_pos _) hy₁ hy

/-- A closed positively oriented triangle containing a point of an open triangle overlaps it. -/
theorem inter_int_nonempty_of_mem_hull {T U : Tri} (hT : T.Pos) {y : Pt} (hy : y ∈ T.hull)
    (hU : y ∈ U.int) : (T.int ∩ U.int).Nonempty := by
  obtain ⟨x, hxU, hxT⟩ :=
    mem_closure_iff.1 (Tri.hull_subset_closure_int hT hy) U.int Tri.isOpen_int hU
  exact ⟨x, hxT, hxU⟩

namespace Overlay

variable {s t : MetricGluing F} (O : Overlay s t) (h₀ : Halfedge F)

/-- The corners of the diamond: `P = i`, `Q = j`, `K = k` from `f₀`, `L = l` from `f₁`. -/
abbrev dP : Pt := O.lay h₀.1 h₀.2
abbrev dQ : Pt := O.lay h₀.1 (h₀.2 + 1)
abbrev dK : Pt := O.lay h₀.1 (h₀.2 + 2)
abbrev dL : Pt := O.lay (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2)

/-- The hypotheses of the flip on an overlay: the edge is flippable with a strictly convex
diamond, and the two faces are laid out in one plane. -/
structure FlipOK : Prop where
  flippable : t.G.Flippable h₀
  convex : t.G.ConvexAt t.ℓ h₀
  coplanar : layΔ t O.lay h₀ = Rigid.id

/-- The cells of the diamond. -/
def Dcells : Finset (Cell F) := O.cells h₀.1 ∪ O.cells (t.G.twin h₀).1

/-- The layouts after the flip. -/
def flipLay : F → Tri := fun f =>
  if f = h₀.1 then rot h₀.2 (O.dP h₀) (O.dL h₀) (O.dK h₀)
  else if f = (t.G.twin h₀).1 then rot (t.G.twin h₀).2 (O.dQ h₀) (O.dK h₀) (O.dL h₀)
  else O.lay f

/-- The cells after the flip. -/
def flipCells : F → Finset (Cell F) := fun f =>
  if f = h₀.1 then (O.Dcells h₀).filter fun c => ((c.tri s).int ∩ (tri3 (O.dP h₀) (O.dL h₀) (O.dK h₀)).int).Nonempty
  else if f = (t.G.twin h₀).1 then
    (O.Dcells h₀).filter fun c => ((c.tri s).int ∩ (tri3 (O.dQ h₀) (O.dK h₀) (O.dL h₀)).int).Nonempty
  else O.cells f

variable {O h₀}

theorem flipLay_f₀ : O.flipLay h₀ h₀.1 = rot h₀.2 (O.dP h₀) (O.dL h₀) (O.dK h₀) := by
  simp [flipLay]

theorem flipLay_f₁ (hf : t.G.Flippable h₀) :
    O.flipLay h₀ (t.G.twin h₀).1 = rot (t.G.twin h₀).2 (O.dQ h₀) (O.dK h₀) (O.dL h₀) := by
  simp [flipLay, show (t.G.twin h₀).1 ≠ h₀.1 from hf]

theorem flipLay_other {f : F} (h0 : f ≠ h₀.1) (h1 : f ≠ (t.G.twin h₀).1) :
    O.flipLay h₀ f = O.lay f := by
  simp [flipLay, h0, h1]

theorem flipCells_f₀ : O.flipCells h₀ h₀.1 = (O.Dcells h₀).filter
    fun c => ((c.tri s).int ∩ (tri3 (O.dP h₀) (O.dL h₀) (O.dK h₀)).int).Nonempty := by
  simp [flipCells]

theorem flipCells_f₁ (hf : t.G.Flippable h₀) : O.flipCells h₀ (t.G.twin h₀).1 = (O.Dcells h₀).filter
    fun c => ((c.tri s).int ∩ (tri3 (O.dQ h₀) (O.dK h₀) (O.dL h₀)).int).Nonempty := by
  simp [flipCells, show (t.G.twin h₀).1 ≠ h₀.1 from hf]

theorem flipCells_other {f : F} (h0 : f ≠ h₀.1) (h1 : f ≠ (t.G.twin h₀).1) :
    O.flipCells h₀ f = O.cells f := by
  simp [flipCells, h0, h1]

theorem mem_Dcells {c : Cell F} :
    c ∈ O.Dcells h₀ ↔ c ∈ O.cells h₀.1 ∨ c ∈ O.cells (t.G.twin h₀).1 := by
  simp [Dcells]

namespace FlipOK

variable (hF : FlipOK O h₀)
include hF

theorem lay_twin_tail : O.lay (t.G.twin h₀).1 (t.G.twin h₀).2 = O.dQ h₀ := by
  have := O.layΔ_tail h₀
  rwa [hF.coplanar, Rigid.id_apply] at this

theorem lay_twin_head : O.lay (t.G.twin h₀).1 ((t.G.twin h₀).2 + 1) = O.dP h₀ := by
  have := O.layΔ_head h₀
  rwa [hF.coplanar, Rigid.id_apply] at this

/-- The old layout of `f₀` is `i, j, k` from the slot `a`. -/
theorem lay_f₀_eq : O.lay h₀.1 = rot h₀.2 (O.dP h₀) (O.dQ h₀) (O.dK h₀) :=
  (rot_eq_self _ _).symm

/-- The old layout of `f₁` is `j, i, l` from the slot `b`. -/
theorem lay_f₁_eq : O.lay (t.G.twin h₀).1 = rot (t.G.twin h₀).2 (O.dQ h₀) (O.dP h₀) (O.dL h₀) := by
  rw [← hF.lay_twin_tail, ← hF.lay_twin_head]
  exact (rot_eq_self _ _).symm

/-- The four points form a strictly convex diamond. -/
theorem quad : ConvexQuad (O.dP h₀) (O.dQ h₀) (O.dK h₀) (O.dL h₀) := by
  have hK : 0 < sda (O.dP h₀) (O.dQ h₀) (O.dK h₀) := by
    have := O.lay_pos h₀.1
    rw [hF.lay_f₀_eq] at this
    exact sda_of_rot_pos this
  have hL : 0 < sda (O.dQ h₀) (O.dP h₀) (O.dL h₀) := by
    have := O.lay_pos (t.G.twin h₀).1
    rw [hF.lay_f₁_eq] at this
    exact sda_of_rot_pos this
  have hL' : sda (O.dP h₀) (O.dQ h₀) (O.dL h₀) < 0 := by
    have := hL
    rw [sda_swap_left] at this; linarith
  have hconv := hF.convex
  unfold Gluing.ConvexAt at hconv
  obtain ⟨hk, hl⟩ := t.isMetric.diamond_isTriangle h₀
  rw [isConvexDiamond_iff_diamondConvex hk hl] at hconv
  unfold DDG.Intrinsic.DiamondConvex at hconv
  rw [DDG.Intrinsic.flipConvex_iff] at hconv
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
  rw [e1, e2, e3, e4, e5] at hconv
  obtain ⟨a1, a2⟩ := DDG.Intrinsic.flipArea_realized hK hL'
  rw [a1, a2] at hconv
  exact ⟨hK, hL, hconv.1, hconv.2⟩

/-- The flip of the twin halfedge is the same flip: the hypotheses hold for it too. -/
theorem symm : FlipOK O (t.G.twin h₀) := by
  refine ⟨?_, ?_, ?_⟩
  · show (t.G.twin (t.G.twin h₀)).1 ≠ (t.G.twin h₀).1
    rw [t.G.twin_twin]
    exact Ne.symm hF.flippable
  · have hconv := hF.convex
    unfold Gluing.ConvexAt IsConvexDiamond at hconv ⊢
    rw [t.G.twin_twin, t.isMetric.twin]
    obtain ⟨h1, h2⟩ := hconv
    constructor
    · rw [cornerAngle_swap' (t.ℓ (t.G.twin h₀).next), cornerAngle_swap' (t.ℓ h₀.prev)]
      linarith
    · rw [cornerAngle_swap' (t.ℓ (t.G.twin h₀).prev), cornerAngle_swap' (t.ℓ h₀.next)]
      linarith
  · refine (O.layΔ_unique (t.G.twin h₀) (g := Rigid.id) ?_ ?_).symm
    · rw [Rigid.id_apply, t.G.twin_twin]; exact hF.lay_twin_head.symm
    · rw [Rigid.id_apply, t.G.twin_twin]; exact hF.lay_twin_tail.symm

theorem f₁_ne : (t.G.twin h₀).1 ≠ h₀.1 := hF.flippable

/-! ### Transport of the layouts across the symmetry -/

theorem dP_twin : O.dP (t.G.twin h₀) = O.dQ h₀ := hF.lay_twin_tail
theorem dQ_twin : O.dQ (t.G.twin h₀) = O.dP h₀ := hF.lay_twin_head
theorem dK_twin : O.dK (t.G.twin h₀) = O.dL h₀ := rfl
theorem dL_twin : O.dL (t.G.twin h₀) = O.dK h₀ := by simp [dL, t.G.twin_twin]

theorem flipLay_twin : O.flipLay (t.G.twin h₀) = O.flipLay h₀ := by
  funext f
  unfold flipLay
  rw [hF.dP_twin, hF.dQ_twin, hF.dK_twin, hF.dL_twin, t.G.twin_twin]
  by_cases h1 : f = (t.G.twin h₀).1
  · subst h1; simp [hF.f₁_ne]
  · by_cases h2 : f = h₀.1
    · subst h2; simp [h1]
    · simp [h1, h2]

theorem Dcells_twin : O.Dcells (t.G.twin h₀) = O.Dcells h₀ := by
  simp [Dcells, t.G.twin_twin, Finset.union_comm]

theorem flipCells_twin : O.flipCells (t.G.twin h₀) = O.flipCells h₀ := by
  funext f
  unfold flipCells
  rw [hF.Dcells_twin, hF.dP_twin, hF.dQ_twin, hF.dK_twin, hF.dL_twin, t.G.twin_twin]
  by_cases h1 : f = (t.G.twin h₀).1
  · subst h1; simp [hF.f₁_ne]
  · by_cases h2 : f = h₀.1
    · subst h2; simp [h1]
    · simp [h1, h2]

/-- The germ at a point of the open old diagonal: the edge germ with the coplanar layouts. -/
theorem germ_edge_D {y : Pt} (hy : y ∈ openSegment ℝ (O.dP h₀) (O.dQ h₀)) :
    ConsistentGerm s (germAt s (O.Dcells h₀) y) y := by
  have hΓ := O.germ_edge h₀ hF.flippable.twin_ne y hy
  have hmv : (Cell.mv Rigid.id : Cell F → Cell F) = id := funext Cell.mv_id
  rw [hF.coplanar, hmv, Finset.image_id] at hΓ
  exact hΓ

/-! ### Crossing the old diagonal -/

/-- **The key lemma**: a cell of the diamond that contains a point of the closed triangle
`PQK` strictly on `K`'s side of the diagonal is a cell of `f₀`. -/
theorem mem_cells_of_sda_pos {c : Cell F} (hc : c ∈ O.Dcells h₀) {y : Pt}
    (hy : y ∈ (c.tri s).hull) (hyK : y ∈ (tri3 (O.dP h₀) (O.dQ h₀) (O.dK h₀)).hull)
    (hpos : 0 < sda (O.dP h₀) (O.dQ h₀) y) : c ∈ O.cells h₀.1 := by
  rcases mem_Dcells.1 hc with hc0 | hc1
  · exact hc0
  have hq := hF.quad
  obtain ⟨w, hwc, hwf⟩ := O.nonjunk _ c hc1
  rw [hF.lay_f₁_eq, rot_int] at hwf
  have hw : sda (O.dP h₀) (O.dQ h₀) w < 0 := by
    have := ((mem_int_tri3 hq.hL).1 hwf).2.2
    rw [sda_swap_left] at this; linarith
  obtain ⟨z, hzs, hz0⟩ := exists_sda_zero_openSegment hw hpos
  have hzc : z ∈ (c.tri s).int :=
    Convex.openSegment_interior_closure_subset_interior Tri.convex_hull hwc (subset_closure hy) hzs
  have hzseg : z ∈ segment ℝ (O.dP h₀) (O.dQ h₀) :=
    hq.cross_mem_segment hyK (Tri.int_subset_hull hwf)
      (by rw [segment_symm]; exact openSegment_subset_segment ℝ w y hzs) hz0
  have hnc : ∀ k : Fin 3, O.lay (t.G.twin h₀).1 k ≠ z := by
    intro k hk
    have hkc : O.lay (t.G.twin h₀).1 k ∈ (c.tri s).hull := by
      rw [hk]; exact Tri.int_subset_hull hzc
    obtain ⟨j, hj⟩ := O.corner _ k c hc1 hkc
    apply Tri.corner_not_mem_int (Cell.tri_pos c) j
    rw [hj, hk]; exact hzc
  have hzo : z ∈ openSegment ℝ (O.dP h₀) (O.dQ h₀) :=
    mem_openSegment_of_ne_left_right (fun e => hnc _ (hF.lay_twin_head.trans e))
      (fun e => hnc _ (hF.lay_twin_tail.trans e)) hzseg
  have hΓ := hF.germ_edge_D hzo
  have hcΓ : c ∈ germAt s (O.Dcells h₀) z := mem_germAt.2 ⟨hc, Tri.int_subset_hull hzc⟩
  have hsing := hΓ.eq_singleton_of_mem_int hcΓ hzc
  have hzhull : z ∈ (O.lay h₀.1).hull := Tri.side_subset_hull (O.lay_pos _) h₀.2 hzo
  obtain ⟨c', hc', hzc'⟩ := Set.mem_iUnion₂.1 (O.cover _ hzhull)
  have hc'Γ : c' ∈ germAt s (O.Dcells h₀) z := mem_germAt.2 ⟨mem_Dcells.2 (Or.inl hc'), hzc'⟩
  rw [hsing, Finset.mem_singleton] at hc'Γ
  rw [← hc'Γ]; exact hc'

/-- The symmetric statement: a cell of the diamond containing a point of `QPL` strictly on `L`'s
side is a cell of `f₁`. -/
theorem mem_cells_of_sda_neg {c : Cell F} (hc : c ∈ O.Dcells h₀) {y : Pt}
    (hy : y ∈ (c.tri s).hull) (hyL : y ∈ (tri3 (O.dQ h₀) (O.dP h₀) (O.dL h₀)).hull)
    (hneg : sda (O.dP h₀) (O.dQ h₀) y < 0) : c ∈ O.cells (t.G.twin h₀).1 := by
  have := hF.symm.mem_cells_of_sda_pos (c := c) (y := y)
  rw [hF.Dcells_twin, dP_twin hF, dQ_twin hF, dK_twin hF] at this
  apply this hc hy hyL
  rw [sda_swap_left]; linarith

/-- The germ of the diamond's cells at any point of the open old diamond is consistent. -/
theorem germ_D {y : Pt}
    (hy : y ∈ (tri3 (O.dP h₀) (O.dQ h₀) (O.dK h₀)).int ∪ openSegment ℝ (O.dP h₀) (O.dQ h₀) ∪
      (tri3 (O.dQ h₀) (O.dP h₀) (O.dL h₀)).int) :
    ConsistentGerm s (germAt s (O.Dcells h₀) y) y := by
  have hq := hF.quad
  rcases hy with (hy | hy) | hy
  · have e : germAt s (O.Dcells h₀) y = germAt s (O.cells h₀.1) y := by
      ext c
      simp only [mem_germAt, mem_Dcells]
      constructor
      · rintro ⟨hc | hc, hyc⟩
        · exact ⟨hc, hyc⟩
        · exact ⟨hF.mem_cells_of_sda_pos (mem_Dcells.2 (Or.inr hc)) hyc (Tri.int_subset_hull hy)
            ((mem_int_tri3 hq.hK).1 hy).2.2, hyc⟩
      · rintro ⟨hc, hyc⟩; exact ⟨Or.inl hc, hyc⟩
    rw [e]
    apply O.germ_int h₀.1 y
    rw [hF.lay_f₀_eq, rot_int]; exact hy
  · exact hF.germ_edge_D hy
  · have hneg : sda (O.dP h₀) (O.dQ h₀) y < 0 := by
      have := ((mem_int_tri3 hq.hL).1 hy).2.2
      rw [sda_swap_left] at this; linarith
    have e : germAt s (O.Dcells h₀) y = germAt s (O.cells (t.G.twin h₀).1) y := by
      ext c
      simp only [mem_germAt, mem_Dcells]
      constructor
      · rintro ⟨hc | hc, hyc⟩
        · exact ⟨hF.mem_cells_of_sda_neg (mem_Dcells.2 (Or.inl hc)) hyc (Tri.int_subset_hull hy)
            hneg, hyc⟩
        · exact ⟨hc, hyc⟩
      · rintro ⟨hc, hyc⟩; exact ⟨Or.inr hc, hyc⟩
    rw [e]
    apply O.germ_int (t.G.twin h₀).1 y
    rw [hF.lay_f₁_eq, rot_int]; exact hy

/-! ### The new face `f₀` -/

theorem flipLay_f₀_pos : (O.flipLay h₀ h₀.1).Pos := by
  rw [flipLay_f₀]; exact rot_pos hF.quad.hc1

theorem flipLay_f₀_hull : (O.flipLay h₀ h₀.1).hull = (tri3 (O.dP h₀) (O.dL h₀) (O.dK h₀)).hull := by
  rw [flipLay_f₀, rot_hull]

theorem flipLay_f₀_int : (O.flipLay h₀ h₀.1).int = (tri3 (O.dP h₀) (O.dL h₀) (O.dK h₀)).int := by
  rw [flipLay_f₀, rot_int]

/-- Every cell of the diamond containing a point of the open new face is a cell of the new
face. -/
theorem germAt_flipCells_f₀_of_mem_int {y : Pt} (hy : y ∈ (tri3 (O.dP h₀) (O.dL h₀) (O.dK h₀)).int) :
    germAt s (O.flipCells h₀ h₀.1) y = germAt s (O.Dcells h₀) y := by
  ext c
  rw [flipCells_f₀, mem_germAt, mem_germAt, Finset.mem_filter]
  constructor
  · rintro ⟨⟨hc, -⟩, hyc⟩; exact ⟨hc, hyc⟩
  · rintro ⟨hc, hyc⟩
    exact ⟨⟨hc, inter_int_nonempty_of_mem_hull (Cell.tri_pos c) hyc hy⟩, hyc⟩

theorem cover_f₀ : (O.flipLay h₀ h₀.1).hull ⊆ ⋃ c ∈ O.flipCells h₀ h₀.1, (c.tri s).hull := by
  intro y hy
  rw [hF.flipLay_f₀_hull] at hy
  have hq := hF.quad
  obtain ⟨ε, hε, hball⟩ := exists_ball_germ (s := s) (O.Dcells h₀) y
  obtain ⟨y', hy'int, hdist⟩ :=
    Metric.mem_closure_iff.1 (Tri.hull_subset_closure_int (tri3_pos hq.hc1) hy) ε hε
  have hy'h := Tri.int_subset_hull hy'int
  obtain ⟨c, hcD, hy'c⟩ : ∃ c ∈ O.Dcells h₀, y' ∈ (c.tri s).hull := by
    rcases hq.hull_new_subset hy'h with h | h
    · have h' : y' ∈ (O.lay h₀.1).hull := by rw [hF.lay_f₀_eq, rot_hull]; exact h
      obtain ⟨c, hc, hc'⟩ := Set.mem_iUnion₂.1 (O.cover _ h')
      exact ⟨c, mem_Dcells.2 (Or.inl hc), hc'⟩
    · have h' : y' ∈ (O.lay (t.G.twin h₀).1).hull := by rw [hF.lay_f₁_eq, rot_hull]; exact h
      obtain ⟨c, hc, hc'⟩ := Set.mem_iUnion₂.1 (O.cover _ h')
      exact ⟨c, mem_Dcells.2 (Or.inr hc), hc'⟩
  have hcy : c ∈ germAt s (O.Dcells h₀) y :=
    hball y' (by rw [dist_comm]; exact hdist) (mem_germAt.2 ⟨hcD, hy'c⟩)
  have hcF : c ∈ O.flipCells h₀ h₀.1 := by
    rw [flipCells_f₀, Finset.mem_filter]
    exact ⟨hcD, inter_int_nonempty_of_mem_hull (Cell.tri_pos c) hy'c hy'int⟩
  exact Set.mem_iUnion₂.2 ⟨c, hcF, (mem_germAt.1 hcy).2⟩

theorem nonjunk_f₀ : ∀ c ∈ O.flipCells h₀ h₀.1, ((c.tri s).int ∩ (O.flipLay h₀ h₀.1).int).Nonempty := by
  intro c hc
  rw [flipCells_f₀, Finset.mem_filter] at hc
  rw [flipLay_f₀_int hF]
  exact hc.2

theorem germ_int_f₀ : ∀ y ∈ (O.flipLay h₀ h₀.1).int,
    ConsistentGerm s (germAt s (O.flipCells h₀ h₀.1) y) y := by
  intro y hy
  rw [flipLay_f₀_int hF] at hy
  rw [germAt_flipCells_f₀_of_mem_int hF hy]
  exact hF.germ_D (hF.quad.int_new_subset hy)

theorem corner_f₀ : ∀ k : Fin 3, ∀ c ∈ O.flipCells h₀ h₀.1,
    O.flipLay h₀ h₀.1 k ∈ (c.tri s).hull → ∃ j, c.tri s j = O.flipLay h₀ h₀.1 k := by
  intro k c hc hk
  have hq := hF.quad
  rw [flipCells_f₀, Finset.mem_filter] at hc
  have hcD := hc.1
  rw [flipLay_f₀] at hk ⊢
  rcases (by omega : k = h₀.2 ∨ k = h₀.2 + 1 ∨ k = h₀.2 + 2) with rfl | rfl | rfl
  · rw [rot_self] at hk ⊢
    rcases mem_Dcells.1 hcD with h0 | h1
    · exact O.corner h₀.1 h₀.2 c h0 hk
    · rw [← hF.lay_twin_head] at hk ⊢
      exact O.corner _ _ c h1 hk
  · rw [rot_add_one] at hk ⊢
    have h1 : c ∈ O.cells (t.G.twin h₀).1 :=
      hF.mem_cells_of_sda_neg hcD hk
        (Tri.corner_mem_hull (T := tri3 (O.dQ h₀) (O.dP h₀) (O.dL h₀)) 2) hq.hL'
    exact O.corner (t.G.twin h₀).1 ((t.G.twin h₀).2 + 2) c h1 hk
  · rw [rot_add_two] at hk ⊢
    have h0 : c ∈ O.cells h₀.1 :=
      hF.mem_cells_of_sda_pos hcD hk
        (Tri.corner_mem_hull (T := tri3 (O.dP h₀) (O.dQ h₀) (O.dK h₀)) 2) hq.hK
    exact O.corner h₀.1 (h₀.2 + 2) c h0 hk

end FlipOK

end Overlay

end DDG.Delaunay.Surface
