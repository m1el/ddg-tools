import DdgProofs.Delaunay.Surface.Angle
import DdgProofs.Delaunay.Termination
import Mathlib.Data.Set.Finite.List

/-!
# Delaunay layer 4: layouts of the faces of a metric triangulation, cells, germs and words

For a metric triangulation `s : MetricGluing F` (`DdgProofs/Delaunay/Termination.lean`):

* `sLayout s σ` is the **standard layout** of the face `σ` in the plane: corner `0` at the
  origin, corner `1` on the positive `x`-axis, corner `2` above (the `apexK` of
  `DDG/Intrinsic/Diamond.lean`); it has the lengths of `s` and is positively oriented.
* `unfold s e` for a halfedge `e = (σ, k)` glued to `e' = (σ', k')` is the rigid motion that
  places the standard layout of `σ'` next to that of `σ` across the side `k`: `e'` runs against
  `e`, so the tail of `e'` goes to the head of `e` and the head of `e'` to the tail of `e`.
* A **cell** `(σ, g)` is a face of `s` placed in the plane by a rigid motion `g`; `Cell.tri` is
  the placed triangle, `Cell.across k` the cell unfolded across its side `k`, `Cell.mv g` the cell
  moved by `g`.
* A **germ** at a point `y` is a finite set of cells; it is *consistent* when it is one cell with
  `y` interior, or two cells, one the unfolding of the other across an interior edge of `s`, with
  `y` on that open side (`ConsistentGerm`); a *boundary germ* is one cell with `y` on an open
  side that is a boundary edge of `s` (`BoundaryGerm`). Both are carried along by rigid motions.
* A **word** is a face, a corner of it, a list of exit sides, and a corner of the final face;
  `Word.devLength` is the distance between the first corner and the last one when the faces
  along the word are unfolded in the plane one after the other (`chainEnd`). Words with at most
  `N` exit sides are finitely many (`Word.finite_length_le`), so their lengths form a finite set.
* `angleAt_tri`: the angle of a placed face at its corner `j` is the corner angle
  `cornerAngleAt` of `s` (`DdgProofs/Delaunay/Metric.lean`).
-/

noncomputable section

namespace DDG.Delaunay.Surface

open DDG.Intrinsic (apexK diamond_isometric diamond_sda dist2_comm)
open DDG.Delaunay.Planar
open Set
open scoped Classical

variable {F : Type*}

/-! ## Standard layouts -/

/-- The standard layout of the face `f` for the lengths `ℓ`: `(0, 0)`, `(ℓ (f, 0), 0)`, and the
apex above the edge. -/
def sLayout (ℓ : Halfedge F → ℝ) (f : F) : Tri :=
  ![(0, 0), (ℓ (f, 0), 0), apexK (ℓ (f, 0)) (ℓ (f, 2)) (ℓ (f, 1))]

section Layout

variable {G : Gluing F} {ℓ : Halfedge F → ℝ}

theorem sLayout_dist (hm : G.IsMetric ℓ) (f : F) (k : Fin 3) :
    dist2 (sLayout ℓ f k) (sLayout ℓ f (k + 1)) = ℓ (f, k) := by
  obtain ⟨-, -, e1, e2, e3, -, -⟩ := diamond_isometric (hm.tri f).swap (hm.tri f).swap
  fin_cases k
  · exact e1
  · exact e3
  · simpa [sLayout, dist2_comm] using e2

theorem sLayout_pos (hm : G.IsMetric ℓ) (f : F) : (sLayout ℓ f).Pos := by
  obtain ⟨h1, -, -, -, -, -, -⟩ := diamond_isometric (hm.tri f).swap (hm.tri f).swap
  have ha := (hm.tri f).pos_a
  show 0 < sda (0, 0) (ℓ (f, 0), 0) (apexK (ℓ (f, 0)) (ℓ (f, 2)) (ℓ (f, 1)))
  rw [apexK, (diamond_sda _ _ _ 0 0).1]
  exact mul_pos ha h1

theorem sLayout_ne (hm : G.IsMetric ℓ) (f : F) {j k : Fin 3} (h : j ≠ k) :
    sLayout ℓ f j ≠ sLayout ℓ f k := by
  intro e
  have hpos := sLayout_pos hm f
  rw [Tri.Pos] at hpos
  fin_cases j <;> fin_cases k <;> simp at h e <;>
    · rw [e] at hpos
      simp [sda_self₁₂, sda_self₂₃, sda_self₁₃] at hpos

end Layout

/-! ## Unfolding across an edge -/

/-- The rigid motion placing the standard layout of the face across the halfedge `e = (σ, k)`
next to the standard layout of `σ`: the twin `e' = (σ', k')` runs against `e`. -/
def unfold (s : MetricGluing F) (e : Halfedge F) : Rigid :=
  rigidOfSegment (sLayout s.ℓ (s.G.twin e).1 (s.G.twin e).2)
    (sLayout s.ℓ (s.G.twin e).1 ((s.G.twin e).2 + 1))
    (sLayout s.ℓ e.1 (e.2 + 1)) (sLayout s.ℓ e.1 e.2)

section Unfold

variable (s : MetricGluing F) (e : Halfedge F)

theorem unfold_dist_eq :
    dist2 (sLayout s.ℓ (s.G.twin e).1 (s.G.twin e).2)
        (sLayout s.ℓ (s.G.twin e).1 ((s.G.twin e).2 + 1)) =
      dist2 (sLayout s.ℓ e.1 (e.2 + 1)) (sLayout s.ℓ e.1 e.2) := by
  rw [sLayout_dist s.isMetric, dist2_comm, sLayout_dist s.isMetric, s.isMetric.twin]

theorem unfold_tail :
    unfold s e (sLayout s.ℓ (s.G.twin e).1 (s.G.twin e).2) = sLayout s.ℓ e.1 (e.2 + 1) :=
  rigidOfSegment_fst (unfold_dist_eq s e) (sLayout_ne s.isMetric _ (by omega))

theorem unfold_head :
    unfold s e (sLayout s.ℓ (s.G.twin e).1 ((s.G.twin e).2 + 1)) = sLayout s.ℓ e.1 e.2 :=
  rigidOfSegment_snd (unfold_dist_eq s e) (sLayout_ne s.isMetric _ (by omega))

end Unfold

/-! ## Cells -/

/-- A placed face: a face of `s` and the rigid motion applied to its standard layout. -/
abbrev Cell (F : Type*) := F × Rigid

namespace Cell

variable (s : MetricGluing F)

/-- The placed triangle. -/
def tri (c : Cell F) : Tri := (sLayout s.ℓ c.1).map c.2

/-- The cell unfolded across its side `k`. -/
def across (c : Cell F) (k : Fin 3) : Cell F :=
  ((s.G.twin (c.1, k)).1, c.2.comp (unfold s (c.1, k)))

/-- The cell moved by `g`. -/
def mv (g : Rigid) (c : Cell F) : Cell F := (c.1, g.comp c.2)

variable {s}

theorem tri_pos (c : Cell F) : (c.tri s).Pos := Tri.map_pos _ (sLayout_pos s.isMetric _)

theorem tri_apply (c : Cell F) (k : Fin 3) : c.tri s k = c.2 (sLayout s.ℓ c.1 k) := rfl

theorem tri_dist (c : Cell F) (k : Fin 3) : dist2 (c.tri s k) (c.tri s (k + 1)) = s.ℓ (c.1, k) := by
  rw [tri_apply, tri_apply, Rigid.dist2_apply, sLayout_dist s.isMetric]

theorem tri_mv (g : Rigid) (c : Cell F) : (c.mv g).tri s = (c.tri s).map g := by
  simp [tri, mv, Tri.map_map]

@[simp] theorem mv_fst (g : Rigid) (c : Cell F) : (c.mv g).1 = c.1 := rfl

theorem across_fst (c : Cell F) (k : Fin 3) : (c.across s k).1 = (s.G.twin (c.1, k)).1 := rfl

theorem mv_across (g : Rigid) (c : Cell F) (k : Fin 3) :
    (c.across s k).mv g = (c.mv g).across s k := by
  simp [across, mv, Rigid.comp_assoc]

theorem mv_mv (g h : Rigid) (c : Cell F) : (c.mv h).mv g = c.mv (g.comp h) := by
  simp [mv, Rigid.comp_assoc]

@[simp] theorem mv_id (c : Cell F) : c.mv Rigid.id = c := by
  obtain ⟨f, g⟩ := c
  simp [mv, Rigid.id_comp]

theorem mv_injective (g : Rigid) : Function.Injective (mv g : Cell F → Cell F) := by
  intro c c' h
  have := congrArg (mv g.inv) h
  rw [mv_mv, mv_mv] at this
  rwa [Rigid.inv_comp, mv_id, mv_id] at this

/-- The unfolded cell shares the side: its corner `k'` (the tail of the twin) is the head of
side `k`, its corner `k' + 1` the tail. -/
theorem across_tail (c : Cell F) (k : Fin 3) :
    (c.across s k).tri s (s.G.twin (c.1, k)).2 = c.tri s (k + 1) := by
  simp only [tri, across, Tri.map, Rigid.comp_apply]
  rw [unfold_tail]

theorem across_head (c : Cell F) (k : Fin 3) :
    (c.across s k).tri s ((s.G.twin (c.1, k)).2 + 1) = c.tri s k := by
  simp only [tri, across, Tri.map, Rigid.comp_apply]
  rw [unfold_head]

theorem across_side (c : Cell F) (k : Fin 3) :
    ((c.across s k).tri s).side (s.G.twin (c.1, k)).2 = (c.tri s).side k :=
  Tri.side_across (across_tail c k) (across_head c k)

/-- Unfolding back and forth: `(c.across k).across k'` is `c` again. -/
theorem across_across (c : Cell F) (k : Fin 3) (hb : s.G.twin (c.1, k) ≠ (c.1, k)) :
    (c.across s k).across s (s.G.twin (c.1, k)).2 = c := by
  obtain ⟨f, g⟩ := c
  simp only [across, s.G.twin_twin]
  congr 1
  -- the composite motion fixes the two corners of the side
  rw [Rigid.comp_assoc]
  have hA : unfold s (s.G.twin (f, k)) (sLayout s.ℓ f (k + 1)) =
      sLayout s.ℓ (s.G.twin (f, k)).1 (s.G.twin (f, k)).2 := by
    have := unfold_head s (s.G.twin (f, k))
    rwa [s.G.twin_twin] at this
  have hB : unfold s (s.G.twin (f, k)) (sLayout s.ℓ f k) =
      sLayout s.ℓ (s.G.twin (f, k)).1 ((s.G.twin (f, k)).2 + 1) := by
    have := unfold_tail s (s.G.twin (f, k))
    rwa [s.G.twin_twin] at this
  have h2 : (unfold s (f, k)).comp (unfold s (s.G.twin (f, k))) = Rigid.id := by
    apply Rigid.ext_of_two (A := sLayout s.ℓ f (k + 1)) (B := sLayout s.ℓ f k)
      (sLayout_ne s.isMetric f (by omega))
    · rw [Rigid.comp_apply, Rigid.id_apply, hA, unfold_tail]
    · rw [Rigid.comp_apply, Rigid.id_apply, hB, unfold_head]
  rw [h2, Rigid.comp_id]

end Cell

/-! ## Germs -/

variable (s : MetricGluing F)

/-- A consistent germ at `y`: one cell with `y` interior, or a cell and its unfolding across an
interior edge of `s` with `y` on that open side. -/
def ConsistentGerm (Γ : Finset (Cell F)) (y : Pt) : Prop :=
  (∃ c : Cell F, Γ = {c} ∧ y ∈ (c.tri s).int) ∨
    (∃ (c : Cell F) (k : Fin 3), s.G.twin (c.1, k) ≠ (c.1, k) ∧ Γ = {c, c.across s k} ∧
      y ∈ (c.tri s).side k)

/-- A boundary germ at `y`: one cell with `y` on an open side that is a boundary edge of `s`. -/
def BoundaryGerm (Γ : Finset (Cell F)) (y : Pt) : Prop :=
  ∃ (c : Cell F) (k : Fin 3), s.G.twin (c.1, k) = (c.1, k) ∧ Γ = {c} ∧ y ∈ (c.tri s).side k

variable {s}

theorem Cell.mem_int_mv (g : Rigid) (c : Cell F) (y : Pt) :
    g y ∈ ((c.mv g).tri s).int ↔ y ∈ (c.tri s).int := by
  rw [tri_mv]; exact Tri.mem_map_int_iff g (tri_pos c)

theorem Cell.mem_side_mv (g : Rigid) (c : Cell F) (k : Fin 3) (y : Pt) :
    g y ∈ ((c.mv g).tri s).side k ↔ y ∈ (c.tri s).side k := by
  rw [tri_mv]; exact Tri.mem_map_side_iff g (tri_pos c)

theorem Cell.mem_hull_mv (g : Rigid) (c : Cell F) (y : Pt) :
    g y ∈ ((c.mv g).tri s).hull ↔ y ∈ (c.tri s).hull := by
  rw [tri_mv]; exact Tri.mem_map_hull_iff g (tri_pos c)

/-- Consistent germs are carried along by rigid motions. -/
theorem ConsistentGerm.mv {Γ : Finset (Cell F)} {y : Pt} (h : ConsistentGerm s Γ y)
    (g : Rigid) : ConsistentGerm s (Γ.image (Cell.mv g)) (g y) := by
  rcases h with ⟨c, rfl, hy⟩ | ⟨c, k, hb, rfl, hy⟩
  · exact Or.inl ⟨c.mv g, by simp, (Cell.mem_int_mv g c y).2 hy⟩
  · refine Or.inr ⟨c.mv g, k, hb, ?_, (Cell.mem_side_mv g c k y).2 hy⟩
    simp [Finset.image_insert, Cell.mv_across]

theorem BoundaryGerm.mv {Γ : Finset (Cell F)} {y : Pt} (h : BoundaryGerm s Γ y) (g : Rigid) :
    BoundaryGerm s (Γ.image (Cell.mv g)) (g y) := by
  obtain ⟨c, k, hb, rfl, hy⟩ := h
  exact ⟨c.mv g, k, hb, by simp, (Cell.mem_side_mv g c k y).2 hy⟩

/-- In a consistent germ every cell contains the point, and the point is neither a corner of a
cell nor ... : the point lies in the closed triangle of each cell. -/
theorem ConsistentGerm.mem_hull {Γ : Finset (Cell F)} {y : Pt} (h : ConsistentGerm s Γ y)
    {c : Cell F} (hc : c ∈ Γ) : y ∈ (c.tri s).hull := by
  rcases h with ⟨c', rfl, hy⟩ | ⟨c', k, -, rfl, hy⟩
  · rw [Finset.mem_singleton] at hc; subst hc; exact Tri.int_subset_hull hy
  · rw [Finset.mem_insert, Finset.mem_singleton] at hc
    rcases hc with rfl | rfl
    · exact Tri.side_subset_hull (Cell.tri_pos _) _ hy
    · rw [← Cell.across_side] at hy
      exact Tri.side_subset_hull (Cell.tri_pos _) _ hy

/-- The point of a consistent germ is not a corner of any of its cells. -/
theorem ConsistentGerm.ne_corner {Γ : Finset (Cell F)} {y : Pt} (h : ConsistentGerm s Γ y)
    {c : Cell F} (hc : c ∈ Γ) (j : Fin 3) : y ≠ c.tri s j := by
  rcases h with ⟨c', rfl, hy⟩ | ⟨c', k, -, rfl, hy⟩
  · rw [Finset.mem_singleton] at hc; subst hc
    rintro rfl; exact Tri.corner_not_mem_int (Cell.tri_pos _) j hy
  · rw [Finset.mem_insert, Finset.mem_singleton] at hc
    rcases hc with rfl | rfl
    · exact Tri.ne_corner_of_mem_side (Cell.tri_pos _) hy j
    · rw [← Cell.across_side] at hy
      exact Tri.ne_corner_of_mem_side (Cell.tri_pos _) hy j

theorem BoundaryGerm.ne_corner {Γ : Finset (Cell F)} {y : Pt} (h : BoundaryGerm s Γ y)
    {c : Cell F} (hc : c ∈ Γ) (j : Fin 3) : y ≠ c.tri s j := by
  obtain ⟨c', k, -, rfl, hy⟩ := h
  rw [Finset.mem_singleton] at hc; subst hc
  exact Tri.ne_corner_of_mem_side (Cell.tri_pos _) hy j

/-- A cell of a consistent germ with the point on its open side `k` has that side interior
(`twin ≠ self`) and its unfolding across `k` in the germ. -/
theorem ConsistentGerm.across_mem {Γ : Finset (Cell F)} {y : Pt} (h : ConsistentGerm s Γ y)
    {c : Cell F} (hc : c ∈ Γ) {k : Fin 3} (hy : y ∈ (c.tri s).side k) :
    s.G.twin (c.1, k) ≠ (c.1, k) ∧ c.across s k ∈ Γ := by
  rcases h with ⟨c', rfl, hy'⟩ | ⟨c', k', hb, rfl, hy'⟩
  · rw [Finset.mem_singleton] at hc; subst hc
    exact absurd hy' (Tri.not_mem_int_of_mem_side (Cell.tri_pos _) hy)
  · rw [Finset.mem_insert, Finset.mem_singleton] at hc
    rcases hc with rfl | rfl
    · -- `y` is on the open sides `k` and `k'` of `c`: they are the same side
      have hk : k = k' := Tri.side_eq_of_mem (Cell.tri_pos _) hy hy'
      subst hk
      exact ⟨hb, by simp⟩
    · -- `y` on side `k` of `c'.across k'`, which is the side `k'` of `c'` as a set
      have hb' : s.G.twin ((c'.across s k').1, (s.G.twin (c'.1, k')).2) ≠
          ((c'.across s k').1, (s.G.twin (c'.1, k')).2) := by
        rw [Cell.across_fst, Prod.mk.eta, s.G.twin_twin]
        exact fun e => hb e.symm
      have hk : k = (s.G.twin (c'.1, k')).2 := by
        rw [← Cell.across_side] at hy'
        exact Tri.side_eq_of_mem (Cell.tri_pos _) hy hy'
      subst hk
      refine ⟨hb', ?_⟩
      rw [Cell.across_across c' k' hb]
      simp

/-! ## Words and their developed length -/

/-- The end of a chain of unfoldings: starting from the cell `c`, unfold across the listed
sides in order. -/
def chainEnd (s : MetricGluing F) : Cell F → List (Fin 3) → Cell F
  | c, [] => c
  | c, k :: ks => chainEnd s (c.across s k) ks

theorem chainEnd_mv (g : Rigid) (c : Cell F) (ks : List (Fin 3)) :
    chainEnd s (c.mv g) ks = (chainEnd s c ks).mv g := by
  induction ks generalizing c with
  | nil => rfl
  | cons k ks ih => simp only [chainEnd]; rw [← Cell.mv_across, ih]

/-- A word: the first face and corner, the exit sides, and the last corner. -/
structure Word (F : Type*) where
  f₀ : F
  c₀ : Fin 3
  es : List (Fin 3)
  c₁ : Fin 3

namespace Word

/-- The number of sides crossed. -/
def numCross (w : Word F) : ℕ := w.es.length

/-- The developed length of a word: the distance between the first corner and the last corner,
all faces unfolded in the plane from the standard layout of the first face. -/
def devLength (s : MetricGluing F) (w : Word F) : ℝ :=
  dist2 (Cell.tri s (w.f₀, Rigid.id) w.c₀) (Cell.tri s (chainEnd s (w.f₀, Rigid.id) w.es) w.c₁)

/-- The developed length is the distance between the corners of the first and the last cell of
any chain with these exit sides. -/
theorem devLength_eq (s : MetricGluing F) (w : Word F) (g : Rigid) :
    w.devLength s =
      dist2 (Cell.tri s (w.f₀, g) w.c₀) (Cell.tri s (chainEnd s (w.f₀, g) w.es) w.c₁) := by
  have e : ((w.f₀, g) : Cell F) = Cell.mv g (w.f₀, Rigid.id) := by
    simp [Cell.mv, Rigid.comp_id]
  rw [e, chainEnd_mv, Cell.tri_mv, Cell.tri_mv, devLength]
  simp only [Tri.map]
  rw [Rigid.dist2_apply]

/-- The words with at most `N` crossings are finitely many. -/
theorem finite_length_le [Finite F] (N : ℕ) : {w : Word F | w.numCross ≤ N}.Finite := by
  let e : Word F → F × Fin 3 × List (Fin 3) × Fin 3 := fun w => (w.f₀, w.c₀, w.es, w.c₁)
  have he : Function.Injective e := by
    rintro ⟨a, b, c, d⟩ ⟨a', b', c', d'⟩ h
    simp only [e, Prod.mk.injEq] at h
    obtain ⟨rfl, rfl, rfl, rfl⟩ := h
    rfl
  have : {w : Word F | w.numCross ≤ N} ⊆ e ⁻¹'
      (Set.univ ×ˢ Set.univ ×ˢ {l : List (Fin 3) | l.length ≤ N} ×ˢ Set.univ) := by
    intro w hw
    simp only [Set.mem_preimage, Set.mem_prod, Set.mem_univ, true_and, Set.mem_ofPred_eq, e,
      and_true]
    exact hw
  refine Set.Finite.subset ?_ this
  exact Set.Finite.preimage he.injOn (Set.Finite.prod Set.finite_univ (Set.Finite.prod
    Set.finite_univ (Set.Finite.prod (List.finite_length_le (Fin 3) N) Set.finite_univ)))

end Word

/-! ## Angles of placed faces -/

/-- The angle of a placed face at its corner `j` is the corner angle of `s` there. -/
theorem angleAt_tri (c : Cell F) (j : Fin 3) :
    angleAt (c.tri s j) (c.tri s (j + 1)) (c.tri s (j + 2)) = Gluing.cornerAngleAt s.ℓ (c.1, j) := by
  unfold angleAt Gluing.cornerAngleAt
  have h1 := Cell.tri_dist (s := s) c (j + 1)
  rw [show j + 1 + 1 = j + 2 by omega] at h1
  have h2 := Cell.tri_dist (s := s) c (j + 2)
  rw [show j + 2 + 1 = j by omega] at h2
  rw [h1, Cell.tri_dist, dist2_comm (c.tri s j), h2]
  rfl

end DDG.Delaunay.Surface
