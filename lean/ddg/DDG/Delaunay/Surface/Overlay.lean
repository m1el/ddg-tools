import DDG.Delaunay.Surface.Layout

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

namespace DDG.Delaunay.Surface

open DDG.Delaunay.Planar
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

/-- The cells of a consistent germ at `y` that contains `c` with `y` on its open side `k` are `c`
and `c.across k`. -/
private theorem germ_mem_pair {Γ : Finset (Cell F)} {y : Pt} (h : ConsistentGerm s Γ y)
    {c : Cell F} (hc : c ∈ Γ) {k : Fin 3} (hy : y ∈ (c.tri s).side k) {c' : Cell F}
    (hc' : c' ∈ Γ) : c' = c ∨ c' = c.across s k := by
  rcases h with ⟨c₁, rfl, -⟩ | ⟨c₁, k₁, hb, rfl, hy₁⟩
  · rw [Finset.mem_singleton] at hc hc'
    exact Or.inl (hc'.trans hc.symm)
  · rw [Finset.mem_insert, Finset.mem_singleton] at hc hc'
    rcases hc with rfl | rfl
    · have hk : k = k₁ := Tri.side_eq_of_mem (Cell.tri_pos _) hy hy₁
      subst hk
      exact hc'
    · have hk : k = (s.G.twin (c₁.1, k₁)).2 := by
        rw [← Cell.across_side] at hy₁
        exact Tri.side_eq_of_mem (Cell.tri_pos _) hy hy₁
      subst hk
      rw [Cell.across_across c₁ k₁ hb]
      exact hc'.symm

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
  set p := O.lay f k with hp
  set d := O.rayDir f k with hd
  have hz : ray p d t₁ ∈ (O.lay f).side k := O.ray_mem_side ht₁
  have hcross := Tri.exit_side_cross_neg (Cell.tri_pos c) hlt hint hside
  obtain ⟨ε₀, hε₀, hball⟩ := exists_ball_germ (s := s) (O.cells f) (ray p d t₁)
  obtain ⟨δ, hδ, hcont⟩ := Metric.continuous_iff.1 (continuous_ray p d) t₁ ε₀ hε₀
  have hm : 0 < min δ (1 - t₁) := lt_min hδ (by linarith [ht₁.2])
  have hm1 : min δ (1 - t₁) ≤ δ := min_le_left _ _
  have hm2 : min δ (1 - t₁) ≤ 1 - t₁ := min_le_right _ _
  have hu1 : t₁ < t₁ + min δ (1 - t₁) / 2 := by linarith
  have hu2 : t₁ + min δ (1 - t₁) / 2 < 1 := by linarith
  have hdist : dist (ray p d (t₁ + min δ (1 - t₁) / 2)) (ray p d t₁) < ε₀ :=
    hcont _ (by rw [Real.dist_eq, abs_of_pos (by linarith)]; linarith)
  obtain ⟨c', hc', hyc'⟩ : ∃ c' ∈ O.cells f,
      ray p d (t₁ + min δ (1 - t₁) / 2) ∈ (c'.tri s).hull := by
    obtain ⟨i, hi, h⟩ :=
      Set.mem_iUnion₂.1 (O.cover f (O.ray_mem_hull ⟨by linarith [ht₁.1], hu2.le⟩))
    exact ⟨i, hi, h⟩
  have hc'z : c' ∈ germAt s (O.cells f) (ray p d t₁) :=
    hball _ hdist (mem_germAt.2 ⟨hc', hyc'⟩)
  have hne : c' ≠ c := by
    rintro rfl
    exact Tri.not_mem_hull_after_side (Cell.tri_pos c') hside hcross hu1 hyc'
  have hzc : ray p d t₁ ∈ (c.tri s).hull := Tri.side_subset_hull (Cell.tri_pos c) k' hside
  by_cases htw : t.G.twin (f, k) = (f, k)
  · obtain ⟨c₁, k₁, -, hΓ, -⟩ : BoundaryGerm s (germAt s (O.cells f) (ray p d t₁)) (ray p d t₁) :=
      O.germ_bdry (f, k) htw _ hz
    have hcz : c ∈ germAt s (O.cells f) (ray p d t₁) := mem_germAt.2 ⟨hc, hzc⟩
    rw [hΓ, Finset.mem_singleton] at hcz hc'z
    exact absurd (hc'z.trans hcz.symm) hne
  · have hΓ : ConsistentGerm s (germAt s (O.cells f ∪ (O.cells (t.G.twin (f, k)).1).image
        (Cell.mv (layΔ t O.lay (f, k)))) (ray p d t₁)) (ray p d t₁) :=
      O.germ_edge (f, k) htw _ hz
    have hcΓ := mem_germAt.2 ⟨Finset.mem_union_left
      ((O.cells (t.G.twin (f, k)).1).image (Cell.mv (layΔ t O.lay (f, k)))) hc, hzc⟩
    have hc'Γ := mem_germAt.2 ⟨Finset.mem_union_left
      ((O.cells (t.G.twin (f, k)).1).image (Cell.mv (layΔ t O.lay (f, k)))) hc',
      (mem_germAt.1 hc'z).2⟩
    refine ⟨(hΓ.across_mem hcΓ hside).1, ?_⟩
    rcases germ_mem_pair hΓ hcΓ hside hc'Γ with h | h
    · exact absurd h hne
    · rw [← h]; exact hc'

/-- A point of the ray strictly between `0` and `1` lying in a cell of `f` is not a corner of
that cell. -/
private theorem ray_ne_corner {c : Cell F} (hc : c ∈ O.cells f) {u : ℝ} (hu : u ∈ Ioo (0 : ℝ) 1)
    (hh : ray (O.lay f k) (O.rayDir f k) u ∈ (c.tri s).hull) (j : Fin 3) :
    ray (O.lay f k) (O.rayDir f k) u ≠ c.tri s j := by
  have hz := O.ray_mem_side (f := f) (k := k) hu
  by_cases htw : t.G.twin (f, k) = (f, k)
  · have hΓ : BoundaryGerm s (germAt s (O.cells f) (ray (O.lay f k) (O.rayDir f k) u))
        (ray (O.lay f k) (O.rayDir f k) u) := O.germ_bdry (f, k) htw _ hz
    exact hΓ.ne_corner (mem_germAt.2 ⟨hc, hh⟩) j
  · have hΓ : ConsistentGerm s (germAt s (O.cells f ∪ (O.cells (t.G.twin (f, k)).1).image
        (Cell.mv (layΔ t O.lay (f, k)))) (ray (O.lay f k) (O.rayDir f k) u))
        (ray (O.lay f k) (O.rayDir f k) u) := O.germ_edge (f, k) htw _ hz
    exact hΓ.ne_corner (mem_germAt.2 ⟨Finset.mem_union_left _ hc, hh⟩) j

/-- The recursion of the walk, by strong induction on the number of cells the ray still meets
after `t₀`: from a cell `c` containing the ray on `[0, δ]` (`t₀ = 0`), or containing the ray at
`t₀` and interior just after it, there is a chain. -/
private theorem chain_aux (n : ℕ) : ∀ (c : Cell F) (t₀ : ℝ),
    ((O.cells f).filter fun c' =>
      ∃ u, t₀ < u ∧ ray (O.lay f k) (O.rayDir f k) u ∈ (c'.tri s).hull).card = n →
    c ∈ O.cells f → 0 ≤ t₀ → t₀ < 1 →
    ((t₀ = 0 ∧ ∃ δ > 0, ∀ u ∈ Icc 0 δ, ray (O.lay f k) (O.rayDir f k) u ∈ (c.tri s).hull) ∨
      (ray (O.lay f k) (O.rayDir f k) t₀ ∈ (c.tri s).hull ∧
        ∃ ε > 0, ∀ u ∈ Ioo t₀ (t₀ + ε), ray (O.lay f k) (O.rayDir f k) u ∈ (c.tri s).int)) →
    Nonempty (Chain s (O.lay f k) (O.rayDir f k) c t₀) := by
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro c t₀ hn hc h0 h1 H
  set p := O.lay f k with hp
  set d := O.rayDir f k with hd
  have hT := Cell.tri_pos (s := s) c
  have hr₀ : ray p d t₀ ∈ (c.tri s).hull := by
    rcases H with ⟨rfl, δ, hδ, hH⟩ | ⟨h, -⟩
    · exact hH 0 ⟨le_rfl, hδ.le⟩
    · exact h
  obtain ⟨t₁, hle, hS⟩ := Tri.exists_exit hT O.rayDir_ne_zero hr₀
  have hlt : t₀ < t₁ := by
    rcases H with ⟨rfl, δ, hδ, hH⟩ | ⟨-, ε, hε, hH⟩
    · have := ((hS δ).1 ⟨hδ.le, hH δ ⟨hδ.le, le_rfl⟩⟩).2
      linarith
    · have := ((hS (t₀ + ε / 2)).1
        ⟨by linarith, Tri.int_subset_hull (hH _ ⟨by linarith, by linarith⟩)⟩).2
      linarith
  have hz : ray p d t₁ ∈ (c.tri s).hull := ((hS t₁).2 ⟨hle, le_rfl⟩).2
  rcases le_or_gt 1 t₁ with h1t | h1t
  · -- the ray ends in `c`
    have hq : ray p d 1 ∈ (c.tri s).hull := ((hS 1).2 ⟨h1.le, h1t⟩).2
    have hq' : O.lay f (k + 1) ∈ (c.tri s).hull := by
      rw [← O.ray_rayDir_one (f := f) (k := k)]; exact hq
    obtain ⟨j, hj⟩ := O.corner f (k + 1) c hc hq'
    exact ⟨Chain.last c t₀ j (by rw [O.ray_rayDir_one, hj])
      fun u hu => ((hS u).2 ⟨hu.1, hu.2.trans h1t⟩).2⟩
  · have ht₁ : t₁ ∈ Ioo (0 : ℝ) 1 := ⟨by linarith, h1t⟩
    have hnc := O.ray_ne_corner hc ht₁ hz
    have hzi : ray p d t₁ ∉ (c.tri s).int := Tri.exit_not_mem_int hT hle hS
    by_cases hex : ∃ u ∈ Ioo t₀ t₁, ray p d u ∈ (c.tri s).int
    · obtain ⟨u, hu, hui⟩ := hex
      have hint : ∀ v ∈ Ioo t₀ t₁, ray p d v ∈ (c.tri s).int := by
        intro v hv
        rcases lt_trichotomy v u with hvu | rfl | hvu
        · exact Tri.mem_int_of_between hT hr₀ hui
            ⟨by rw [min_eq_left hu.1.le]; exact hv.1, by rw [max_eq_right hu.1.le]; exact hvu⟩
        · exact hui
        · exact Tri.mem_int_of_between hT hz hui
            ⟨by rw [min_eq_right hu.2.le]; exact hvu, by rw [max_eq_left hu.2.le]; exact hv.2⟩
      obtain ⟨k', hk'⟩ := Tri.exists_side_of_mem_hull hT hz hzi hnc
      obtain ⟨-, hc'⟩ := O.across_mem_cells hc hlt ht₁ hint hk'
      have hcross := Tri.exit_side_cross_neg hT hlt hint hk'
      have hside' : ray p d t₁ ∈ ((c.across s k').tri s).side (s.G.twin (c.1, k')).2 := by
        rw [Cell.across_side]; exact hk'
      have hcross' : 0 < cross2 ((c.across s k').tri s ((s.G.twin (c.1, k')).2 + 1) -
          (c.across s k').tri s (s.G.twin (c.1, k')).2) d := by
        rw [Tri.cross_across (Cell.across_tail c k') (Cell.across_head c k')]
        linarith
      obtain ⟨ε, hε, hH⟩ := Tri.mem_int_after_side (Cell.tri_pos _) hside' hcross'
      have hsub : ((O.cells f).filter fun c' => ∃ u, t₁ < u ∧ ray p d u ∈ (c'.tri s).hull) ⊆
          ((O.cells f).filter fun c' => ∃ u, t₀ < u ∧ ray p d u ∈ (c'.tri s).hull) := by
        intro c'' hm
        obtain ⟨hm1, u, hu, hu'⟩ := Finset.mem_filter.1 hm
        exact Finset.mem_filter.2 ⟨hm1, u, by linarith, hu'⟩
      have hlt' : ((O.cells f).filter fun c' =>
          ∃ u, t₁ < u ∧ ray p d u ∈ (c'.tri s).hull).card < n := by
        rw [← hn]
        apply Finset.card_lt_card
        refine (Finset.ssubset_iff_of_subset hsub).2
          ⟨c, Finset.mem_filter.2 ⟨hc, t₁, hlt, hz⟩, fun hm => ?_⟩
        obtain ⟨-, u, hu, hu'⟩ := Finset.mem_filter.1 hm
        have := ((hS u).1 ⟨by linarith, hu'⟩).2
        linarith
      obtain ⟨C⟩ := ih _ hlt' (c.across s k') t₁ rfl hc' ht₁.1.le h1t
        (Or.inr ⟨Tri.side_subset_hull (Cell.tri_pos _) _ hside', ε, hε, hH⟩)
      exact ⟨Chain.step c t₀ t₁ k' hlt h1t hint hk' C⟩
    · push Not at hex
      exfalso
      rcases H with ⟨rfl, δ, hδ, hH⟩ | ⟨-, ε, hε, hH⟩
      · have hp0 : O.lay f k ∈ (c.tri s).hull := by
          have := hr₀
          rwa [ray_zero] at this
        obtain ⟨j₀, hj₀⟩ := O.corner f k c hc hp0
        obtain ⟨j, hj⟩ := Tri.exit_corner_of_not_int hT O.rayDir_ne_zero hj₀.symm hS hlt hex
        exact hnc j hj
      · have hm : 0 < min ε (t₁ - t₀) := lt_min hε (by linarith)
        have hm1 : min ε (t₁ - t₀) ≤ ε := min_le_left _ _
        have hm2 : min ε (t₁ - t₀) ≤ t₁ - t₀ := min_le_right _ _
        exact hex (t₀ + min ε (t₁ - t₀) / 2) ⟨by linarith, by linarith⟩
          (hH _ ⟨by linarith, by linarith⟩)

/-- **The walk**: along the halfedge `(f, k)` of `t` there is a chain of cells of `f` from a
cell having the corner `k` of the face as a corner. -/
theorem exists_chain (f : F) (k : Fin 3) :
    ∃ (c : Cell F) (j₀ : Fin 3) (_ : c ∈ O.cells f),
      O.lay f k = c.tri s j₀ ∧ Nonempty (Chain s (O.lay f k) (O.rayDir f k) c 0) := by
  obtain ⟨c, hc, δ, hδ, hH⟩ := exists_initial (O.cells f) (fun c => (c.tri s).hull)
    (fun _ _ => Tri.isClosed_hull) (fun _ _ => Tri.convex_hull)
    (p := O.lay f k) (d := O.rayDir f k)
    (fun u hu => by
      obtain ⟨i, hi, h⟩ := Set.mem_iUnion₂.1 (O.cover f (O.ray_mem_hull (f := f) (k := k) hu))
      exact ⟨i, hi, h⟩)
  have h0 : O.lay f k ∈ (c.tri s).hull := by
    have := hH 0 ⟨le_rfl, hδ.le⟩
    rwa [ray_zero] at this
  obtain ⟨j₀, hj₀⟩ := O.corner f k c hc h0
  exact ⟨c, j₀, hc, hj₀.symm,
    O.chain_aux _ c 0 rfl hc le_rfl one_pos (Or.inl ⟨rfl, δ, hδ, hH⟩)⟩

/-- **Every edge length of `t` is the developed length of a word of `s`.** -/
theorem exists_word (f : F) (k : Fin 3) :
    ∃ (c : Cell F) (j₀ : Fin 3) (C : Chain s (O.lay f k) (O.rayDir f k) c 0),
      c ∈ O.cells f ∧ (C.word j₀).devLength s = t.ℓ (f, k) := by
  obtain ⟨c, j₀, hc, hp, ⟨C⟩⟩ := O.exists_chain f k
  refine ⟨c, j₀, C, hc, ?_⟩
  rw [C.devLength_eq hp, ray_rayDir_one, O.lay_dist]

end Walk

end Overlay

end DDG.Delaunay.Surface
