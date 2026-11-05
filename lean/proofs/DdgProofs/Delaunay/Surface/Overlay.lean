import DdgProofs.Delaunay.Surface.Layout

/-!
# Delaunay layer 4: the overlay of a triangulation on the original one, and the walk

An **overlay** `Overlay s t` of a metric triangulation `t` on `s` (both on the face type `F`)
witnesses that `t` is a geodesic triangulation of the piecewise flat surface of `s` with the same
vertices, without ever building that surface:

* every face `τ` of `t` is laid out in the plane (`lay τ`, a positively oriented triangle with
  the lengths of `t`), and carries a finite set of **cells** (`cells τ`): faces of `s` placed in
  the same plane (`Layout.lean`);
* the cells cover the face (`cover`); every cell overlaps the open face (`nonjunk`);
* at every point `y` of the open face, the cells containing `y` form a consistent germ
  (`germ_int`): one cell with `y` interior, or a cell and its unfolding across an interior edge of
  `s` with `y` on that side;
* at every point `y` of an open side `x` of `τ` the same holds for the cells of `τ` together with
  the cells of the face across `x`, placed next to `τ` by the rigid motion `Δ x` that matches the
  two copies of the side (`germ_edge`); on a boundary side the cells of `τ` containing `y` form
  a boundary germ (`germ_bdry`);
* a corner of `τ` that lies in a cell is a corner of that cell (`corner`).

So `s`'s edges run through `t`'s faces as straight segments glued flatly across the edges of
`t`, and `t`'s vertices are vertices of `s` — the overlay of two geodesic triangulations of one
surface. Everything is local and planar; "phantom" placements (developing around a cone point)
are excluded because every condition is triggered at a point.

**The walk** (`Overlay.exists_chain`): along every halfedge `x = p → q` of `t` the ray from `p`
to `q` passes through a chain of cells of `x`'s face (`Chain`): it starts at a corner of the first
cell, is interior to each cell between two crossings, crosses from a cell into its unfolding
across an interior edge of `s` at a point of that open side, and ends at a corner of the last
cell. Consequently `t.ℓ x` is the developed length of the word of the chain
(`Chain.devLength_eq`, `Overlay.exists_word`). The bound on the number of crossings by the
length is in `Fan.lean`.
-/

noncomputable section

namespace DdgProofs.Delaunay.Surface

open DdgProofs.Intrinsic (dist2_comm)
open DdgProofs.Delaunay.Planar
open Set
open scoped Classical

variable {F : Type*}

/-- The cells of `S` whose closed triangle contains `y`. -/
def germAt (s : MetricGluing F) (S : Finset (Cell F)) (y : Pt) : Finset (Cell F) :=
  S.filter fun c => y ∈ (c.tri s).hull

theorem mem_germAt {s : MetricGluing F} {S : Finset (Cell F)} {y : Pt} {c : Cell F} :
    c ∈ germAt s S y ↔ c ∈ S ∧ y ∈ (c.tri s).hull := by
  simp [germAt]

/-- The triangle of the face `f` of `t` laid out by `lay`, with the neighbouring face across the
halfedge `x` placed next to it: the rigid motion taking the twin side onto the side `x`
(reversed). -/
def layΔ (t : MetricGluing F) (lay : F → Tri) (x : Halfedge F) : Rigid :=
  rigidOfSegment (lay (t.G.twin x).1 (t.G.twin x).2) (lay (t.G.twin x).1 ((t.G.twin x).2 + 1))
    (lay x.1 (x.2 + 1)) (lay x.1 x.2)

/-- An overlay of `t` on `s`. -/
structure Overlay (s t : MetricGluing F) where
  /-- The layout of each face of `t`. -/
  lay : F → Tri
  lay_dist : ∀ (f : F) (k : Fin 3), dist2 (lay f k) (lay f (k + 1)) = t.ℓ (f, k)
  lay_pos : ∀ f, (lay f).Pos
  /-- The cells of each face of `t`: faces of `s` placed in the plane of the face. -/
  cells : F → Finset (Cell F)
  cover : ∀ f, (lay f).hull ⊆ ⋃ c ∈ cells f, (c.tri s).hull
  nonjunk : ∀ f, ∀ c ∈ cells f, ((c.tri s).int ∩ (lay f).int).Nonempty
  germ_int : ∀ f, ∀ y ∈ (lay f).int, ConsistentGerm s (germAt s (cells f) y) y
  germ_edge : ∀ x : Halfedge F, t.G.twin x ≠ x → ∀ y ∈ (lay x.1).side x.2,
    ConsistentGerm s
      (germAt s (cells x.1 ∪ (cells (t.G.twin x).1).image (Cell.mv (layΔ t lay x))) y) y
  germ_bdry : ∀ x : Halfedge F, t.G.twin x = x → ∀ y ∈ (lay x.1).side x.2,
    BoundaryGerm s (germAt s (cells x.1) y) y
  corner : ∀ (f : F) (k : Fin 3), ∀ c ∈ cells f, lay f k ∈ (c.tri s).hull →
    ∃ j, c.tri s j = lay f k

namespace Overlay

variable {s t : MetricGluing F} (O : Overlay s t)

theorem lay_ne (f : F) {j k : Fin 3} (h : j ≠ k) : O.lay f j ≠ O.lay f k := by
  intro e
  have hpos := O.lay_pos f
  rw [Tri.Pos] at hpos
  fin_cases j <;> fin_cases k <;> simp at h e <;>
    · rw [e] at hpos
      simp [sda_self₁₂, sda_self₂₃, sda_self₁₃] at hpos

theorem layΔ_dist (x : Halfedge F) :
    dist2 (O.lay (t.G.twin x).1 (t.G.twin x).2) (O.lay (t.G.twin x).1 ((t.G.twin x).2 + 1)) =
      dist2 (O.lay x.1 (x.2 + 1)) (O.lay x.1 x.2) := by
  rw [O.lay_dist, dist2_comm, O.lay_dist, t.isMetric.twin]

/-- `Δ x` takes the tail of the twin to the head of `x`. -/
theorem layΔ_tail (x : Halfedge F) :
    layΔ t O.lay x (O.lay (t.G.twin x).1 (t.G.twin x).2) = O.lay x.1 (x.2 + 1) :=
  rigidOfSegment_fst (O.layΔ_dist x) (O.lay_ne _ (by omega))

/-- `Δ x` takes the head of the twin to the tail of `x`. -/
theorem layΔ_head (x : Halfedge F) :
    layΔ t O.lay x (O.lay (t.G.twin x).1 ((t.G.twin x).2 + 1)) = O.lay x.1 x.2 :=
  rigidOfSegment_snd (O.layΔ_dist x) (O.lay_ne _ (by omega))

/-- `Δ x` is the only rigid motion that does so. -/
theorem layΔ_unique (x : Halfedge F) {g : Rigid}
    (h0 : g (O.lay (t.G.twin x).1 (t.G.twin x).2) = O.lay x.1 (x.2 + 1))
    (h1 : g (O.lay (t.G.twin x).1 ((t.G.twin x).2 + 1)) = O.lay x.1 x.2) :
    g = layΔ t O.lay x :=
  rigidOfSegment_unique (O.layΔ_dist x) (O.lay_ne _ (by omega)) h0 h1

/-- The placed neighbour shares the open side. -/
theorem side_layΔ (x : Halfedge F) :
    ((O.lay (t.G.twin x).1).map (layΔ t O.lay x)).side (t.G.twin x).2 = (O.lay x.1).side x.2 :=
  Tri.side_across (T := O.lay x.1) (k := x.2) (by simp [Tri.map, O.layΔ_tail])
    (by simp [Tri.map, O.layΔ_head])

/-! ### Local determination of germs -/

/-- Every point near `z` lies only in cells of `S` that contain `z`. -/
theorem exists_ball_germ (S : Finset (Cell F)) (z : Pt) :
    ∃ ε > 0, ∀ y, dist y z < ε → germAt s S y ⊆ germAt s S z := by
  obtain ⟨ε, hε, h⟩ := exists_ball_subset S (fun c => (c.tri s).hull) (fun _ _ => Tri.isClosed_hull) z
  refine ⟨ε, hε, fun y hy c hc => ?_⟩
  rw [mem_germAt] at hc ⊢
  exact ⟨hc.1, h y hy c hc.1 hc.2⟩

/-- A corner of a face never lies in the open interior nor on an open side of a cell of its
germ... : the point of a consistent germ is not a corner of its cells (restated). -/
theorem corner_ne_of_germ_int (f : F) {y : Pt} (hy : y ∈ (O.lay f).int) {c : Cell F}
    (hc : c ∈ O.cells f) (hyc : y ∈ (c.tri s).hull) (j : Fin 3) : y ≠ c.tri s j :=
  (O.germ_int f y hy).ne_corner (mem_germAt.2 ⟨hc, hyc⟩) j

/-! ## Chains -/

/-- A chain of cells along the ray `ray p d` from the parameter `t₀` on: either the current cell
`c` contains the rest of the ray up to the parameter `1`, where the ray ends at a corner of `c`
(`last`), or the ray is interior to `c` on `(t₀, t₁)`, leaves `c` through its open side `k` at
`t₁ < 1`, and continues in the cell unfolded across `k` (`step`). -/
inductive Chain (s : MetricGluing F) (p d : Pt) : Cell F → ℝ → Type _
  | last (c : Cell F) (t₀ : ℝ) (j : Fin 3) (hq : ray p d 1 = c.tri s j)
      (hS : ∀ t ∈ Icc t₀ 1, ray p d t ∈ (c.tri s).hull) : Chain s p d c t₀
  | step (c : Cell F) (t₀ t₁ : ℝ) (k : Fin 3) (hlt : t₀ < t₁) (h₁ : t₁ < 1)
      (hint : ∀ t ∈ Ioo t₀ t₁, ray p d t ∈ (c.tri s).int) (hside : ray p d t₁ ∈ (c.tri s).side k)
      (rest : Chain s p d (c.across s k) t₁) : Chain s p d c t₀

namespace Chain

variable {p d : Pt}

/-- The exit sides, in order. -/
def exits : {c : Cell F} → {t₀ : ℝ} → Chain s p d c t₀ → List (Fin 3)
  | _, _, last _ _ _ _ _ => []
  | _, _, step _ _ _ k _ _ _ _ rest => k :: rest.exits

/-- The last cell. -/
def endCell : {c : Cell F} → {t₀ : ℝ} → Chain s p d c t₀ → Cell F
  | c, _, last _ _ _ _ _ => c
  | _, _, step _ _ _ _ _ _ _ _ rest => rest.endCell

/-- The corner of the last cell at the end of the ray. -/
def endCorner : {c : Cell F} → {t₀ : ℝ} → Chain s p d c t₀ → Fin 3
  | _, _, last _ _ j _ _ => j
  | _, _, step _ _ _ _ _ _ _ _ rest => rest.endCorner

/-- The number of crossings. -/
def numCross {c : Cell F} {t₀ : ℝ} (C : Chain s p d c t₀) : ℕ := C.exits.length

theorem endCell_eq {c : Cell F} {t₀ : ℝ} (C : Chain s p d c t₀) :
    C.endCell = chainEnd s c C.exits := by
  induction C with
  | last => rfl
  | step c t₀ t₁ k _ _ _ _ rest ih => simp only [endCell, exits, chainEnd]; exact ih

theorem ray_one_eq {c : Cell F} {t₀ : ℝ} (C : Chain s p d c t₀) :
    ray p d 1 = C.endCell.tri s C.endCorner := by
  induction C with
  | last c t₀ j hq hS => exact hq
  | step _ _ _ _ _ _ _ _ rest ih => exact ih

/-- The word of a chain that starts at the corner `j₀` of its first cell. -/
def word {c : Cell F} {t₀ : ℝ} (C : Chain s p d c t₀) (j₀ : Fin 3) : Word F :=
  ⟨c.1, j₀, C.exits, C.endCorner⟩

/-- **The length of the segment is the developed length of the word of the chain**: for a chain
from the parameter `0` whose start is the corner `j₀` of the first cell. -/
theorem devLength_eq {c : Cell F} (C : Chain s p d c 0) {j₀ : Fin 3} (hp : p = c.tri s j₀) :
    (C.word j₀).devLength s = dist2 p (ray p d 1) := by
  subst hp
  rw [Word.devLength_eq s _ c.2, C.ray_one_eq, C.endCell_eq]
  rfl

end Chain

/-! ## The walk -/

section Walk

variable {f : F} {k : Fin 3}

/-- The data of the walk along the halfedge `x = (f, k)` of `t`: the ray from the corner `k` of
the face to the corner `k + 1`. -/
def rayDir (f : F) (k : Fin 3) : Pt := O.lay f (k + 1) - O.lay f k

theorem rayDir_ne_zero : O.rayDir f k ≠ 0 := by
  intro h
  apply O.lay_ne f (j := k + 1) (k := k) (by omega)
  rwa [rayDir, sub_eq_zero] at h

theorem ray_rayDir_one : ray (O.lay f k) (O.rayDir f k) 1 = O.lay f (k + 1) := by
  simp [ray, rayDir]

/-- The points of the ray strictly between `0` and `1` lie on the open side `k`. -/
theorem ray_mem_side {u : ℝ} (hu : u ∈ Ioo (0 : ℝ) 1) :
    ray (O.lay f k) (O.rayDir f k) u ∈ (O.lay f).side k := by
  rw [Tri.side, openSegment_eq_image]
  refine ⟨u, hu, ?_⟩
  simp only [ray, rayDir]
  ext <;> simp <;> ring

/-- The points of the ray in `[0, 1]` lie in the closed face. -/
theorem ray_mem_hull {u : ℝ} (hu : u ∈ Icc (0 : ℝ) 1) :
    ray (O.lay f k) (O.rayDir f k) u ∈ (O.lay f).hull := by
  rcases hu.1.lt_or_eq with h0 | rfl
  · rcases hu.2.lt_or_eq with h1 | rfl
    · exact Tri.side_subset_hull (O.lay_pos f) k (O.ray_mem_side ⟨h0, h1⟩)
    · rw [ray_rayDir_one]; exact Tri.corner_mem_hull _
  · rw [ray_zero]; exact Tri.corner_mem_hull _

/-- **The step of the walk**: the ray at the parameter `t₁ ∈ (0, 1)` lies on the open side `k'`
of the cell `c ∈ cells f`, having been interior to `c` just before; then the side is an
interior edge of `s` and the unfolded cell `c.across k'` is a cell of `f`. -/
theorem across_mem_cells {c : Cell F} (hc : c ∈ O.cells f) {t₀ t₁ : ℝ} (hlt : t₀ < t₁)
    (ht₁ : t₁ ∈ Ioo (0 : ℝ) 1)
    (hint : ∀ u ∈ Ioo t₀ t₁, ray (O.lay f k) (O.rayDir f k) u ∈ (c.tri s).int) {k' : Fin 3}
    (hside : ray (O.lay f k) (O.rayDir f k) t₁ ∈ (c.tri s).side k') :
    s.G.twin (c.1, k') ≠ (c.1, k') ∧ c.across s k' ∈ O.cells f := by
  sorry

/-- **The walk**: along the halfedge `(f, k)` of `t` there is a chain of cells of `f` from a
cell having the corner `k` of the face as a corner. -/
theorem exists_chain (f : F) (k : Fin 3) :
    ∃ (c : Cell F) (j₀ : Fin 3) (_ : c ∈ O.cells f),
      O.lay f k = c.tri s j₀ ∧ Nonempty (Chain s (O.lay f k) (O.rayDir f k) c 0) := by
  sorry

/-- **Every edge length of `t` is the developed length of a word of `s`.** -/
theorem exists_word (f : F) (k : Fin 3) :
    ∃ (c : Cell F) (j₀ : Fin 3) (C : Chain s (O.lay f k) (O.rayDir f k) c 0),
      c ∈ O.cells f ∧ (C.word j₀).devLength s = t.ℓ (f, k) := by
  obtain ⟨c, j₀, hc, hp, ⟨C⟩⟩ := O.exists_chain f k
  refine ⟨c, j₀, C, hc, ?_⟩
  rw [C.devLength_eq hp, ray_rayDir_one, O.lay_dist]

end Walk

end Overlay

end DdgProofs.Delaunay.Surface
