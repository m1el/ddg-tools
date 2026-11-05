import DdgProofs.Delaunay.Surface.Overlay

/-!
# Delaunay layer 4: moving the layout of one face

An overlay lays out every face of `t` in its own plane; the plane of a face can be changed by a
rigid motion `R` without changing anything else: the face's layout and its cells move by `R`
(`Overlay.relayout`). The motions `layΔ` matching neighbouring faces change accordingly
(`relayout_layΔ_*`). This is used to lay out the two faces of a flip in one plane
(`FlipFace.lean`): relaying `f₁` by `layΔ t lay h₀` makes `layΔ` of `h₀` the identity
(`relayout_layΔ_eq_id`).
-/

noncomputable section

namespace DdgProofs.Delaunay.Surface

open DdgProofs.Delaunay.Planar
open Set
open scoped Classical

variable {F : Type*}

/-- Moving a finite set of cells commutes with taking the germ at a point. -/
theorem germAt_image_mv (s : MetricGluing F) (S : Finset (Cell F)) (R : Rigid) (y : Pt) :
    germAt s (S.image (Cell.mv R)) (R y) = (germAt s S y).image (Cell.mv R) := by
  ext c
  simp only [mem_germAt, Finset.mem_image]
  constructor
  · rintro ⟨⟨c₀, hc₀, rfl⟩, h⟩
    exact ⟨c₀, ⟨hc₀, (Cell.mem_hull_mv R c₀ y).1 h⟩, rfl⟩
  · rintro ⟨c₀, ⟨hc₀, h⟩, rfl⟩
    exact ⟨⟨c₀, hc₀, rfl⟩, (Cell.mem_hull_mv R c₀ y).2 h⟩

private theorem rigid_inv_id : Rigid.id.inv = Rigid.id := by
  have := Rigid.inv_comp Rigid.id
  rwa [Rigid.comp_id] at this

namespace Overlay

variable {s t : MetricGluing F} (O : Overlay s t)

/-! ### Moving every face by its own rigid motion -/

/-- The matching motions after laying out every face `f` by `Q f`. -/
theorem layΔ_map (Q : F → Rigid) (x : Halfedge F) :
    layΔ t (fun f => (O.lay f).map (Q f)) x =
      (Q x.1).comp ((layΔ t O.lay x).comp (Q (t.G.twin x).1).inv) := by
  symm
  unfold layΔ
  apply rigidOfSegment_unique
  · simp only [Tri.map, Rigid.dist2_apply]
    exact O.layΔ_dist x
  · simp only [Tri.map, ne_eq, Rigid.apply_inj]
    exact O.lay_ne _ (by omega)
  · simp only [Tri.map, Rigid.comp_apply, Rigid.inv_apply_apply]
    rw [← layΔ, O.layΔ_tail]
  · simp only [Tri.map, Rigid.comp_apply, Rigid.inv_apply_apply]
    rw [← layΔ, O.layΔ_head]

/-- The overlay with every face `f` laid out by `Q f`. -/
def remap (Q : F → Rigid) : Overlay s t where
  lay f := (O.lay f).map (Q f)
  lay_dist := by
    intro f k
    simp only [Tri.map, Rigid.dist2_apply]
    exact O.lay_dist f k
  lay_pos f := Tri.map_pos _ (O.lay_pos f)
  cells f := (O.cells f).image (Cell.mv (Q f))
  cover := by
    intro f y hy
    rw [Tri.map_hull _ (O.lay_pos f)] at hy
    obtain ⟨y₀, hy₀, rfl⟩ := hy
    obtain ⟨c, hc, hyc⟩ := Set.mem_iUnion₂.1 (O.cover f hy₀)
    exact Set.mem_iUnion₂.2 ⟨c.mv (Q f), Finset.mem_image_of_mem _ hc,
      (Cell.mem_hull_mv _ c y₀).2 hyc⟩
  nonjunk := by
    intro f c hc
    obtain ⟨c₀, hc₀, rfl⟩ := Finset.mem_image.1 hc
    obtain ⟨z, hz1, hz2⟩ := O.nonjunk f c₀ hc₀
    exact ⟨Q f z, (Cell.mem_int_mv _ _ _).2 hz1, (Tri.mem_map_int_iff _ (O.lay_pos f)).2 hz2⟩
  germ_int := by
    intro f y hy
    rw [Tri.map_int _ (O.lay_pos f)] at hy
    obtain ⟨y₀, hy₀, rfl⟩ := hy
    rw [germAt_image_mv]
    exact (O.germ_int f y₀ hy₀).mv _
  germ_edge := by
    intro x hx y hy
    rw [Tri.map_side _ (O.lay_pos _)] at hy
    obtain ⟨y₀, hy₀, rfl⟩ := hy
    rw [O.layΔ_map Q x]
    have hU : (O.cells x.1).image (Cell.mv (Q x.1)) ∪
        ((O.cells (t.G.twin x).1).image (Cell.mv (Q (t.G.twin x).1))).image
          (Cell.mv ((Q x.1).comp ((layΔ t O.lay x).comp (Q (t.G.twin x).1).inv))) =
        (O.cells x.1 ∪ (O.cells (t.G.twin x).1).image (Cell.mv (layΔ t O.lay x))).image
          (Cell.mv (Q x.1)) := by
      rw [Finset.image_union, Finset.image_image, Finset.image_image]
      congr 2
      funext c
      simp only [Function.comp, Cell.mv_mv, Rigid.comp_assoc, Rigid.inv_comp, Rigid.comp_id]
    rw [hU, germAt_image_mv]
    exact (O.germ_edge x hx y₀ hy₀).mv _
  germ_bdry := by
    intro x hx y hy
    rw [Tri.map_side _ (O.lay_pos _)] at hy
    obtain ⟨y₀, hy₀, rfl⟩ := hy
    rw [germAt_image_mv]
    exact (O.germ_bdry x hx y₀ hy₀).mv _
  corner := by
    intro f k c hc h
    obtain ⟨c₀, hc₀, rfl⟩ := Finset.mem_image.1 hc
    have h' : O.lay f k ∈ (c₀.tri s).hull := (Cell.mem_hull_mv (Q f) c₀ _).1 h
    obtain ⟨j, hj⟩ := O.corner f k c₀ hc₀ h'
    refine ⟨j, ?_⟩
    rw [Cell.tri_mv]
    show Q f (c₀.tri s j) = Q f (O.lay f k)
    rw [hj]

theorem remap_lay (Q : F → Rigid) : (O.remap Q).lay = fun f => (O.lay f).map (Q f) := rfl

theorem remap_cells (Q : F → Rigid) :
    (O.remap Q).cells = fun f => (O.cells f).image (Cell.mv (Q f)) := rfl

variable (g : F) (R : Rigid)

/-- The motion of the face `f` when only `g` is moved by `R`. -/
def rf (f : F) : Rigid := if f = g then R else Rigid.id

theorem ite_lay (f : F) :
    (if f = g then (O.lay g).map R else O.lay f) = (O.lay f).map (rf g R f) := by
  unfold rf
  split_ifs with h
  · subst h; rfl
  · simp

theorem ite_cells (f : F) :
    (if f = g then (O.cells g).image (Cell.mv R) else O.cells f) =
      (O.cells f).image (Cell.mv (rf g R f)) := by
  unfold rf
  split_ifs with h
  · subst h; rfl
  · ext c
    simp only [Finset.mem_image, Cell.mv_id, exists_eq_right]

/-- The overlay with the face `g` laid out by `R` instead. -/
def relayout : Overlay s t where
  lay f := if f = g then (O.lay g).map R else O.lay f
  lay_dist := by
    have h := (O.remap (rf g R)).lay_dist
    simp only [remap_lay] at h
    simpa only [ite_lay, ite_cells] using h
  lay_pos := by
    have h := (O.remap (rf g R)).lay_pos
    simp only [remap_lay] at h
    simpa only [ite_lay, ite_cells] using h
  cells f := if f = g then (O.cells g).image (Cell.mv R) else O.cells f
  cover := by
    have h := (O.remap (rf g R)).cover
    simp only [remap_lay, remap_cells] at h
    simpa only [ite_lay, ite_cells] using h
  nonjunk := by
    have h := (O.remap (rf g R)).nonjunk
    simp only [remap_lay, remap_cells] at h
    simpa only [ite_lay, ite_cells] using h
  germ_int := by
    have h := (O.remap (rf g R)).germ_int
    simp only [remap_lay, remap_cells] at h
    simpa only [ite_lay, ite_cells] using h
  germ_edge := by
    have h := (O.remap (rf g R)).germ_edge
    simp only [remap_lay, remap_cells] at h
    simpa only [ite_lay, ite_cells] using h
  germ_bdry := by
    have h := (O.remap (rf g R)).germ_bdry
    simp only [remap_lay, remap_cells] at h
    simpa only [ite_lay, ite_cells] using h
  corner := by
    have h := (O.remap (rf g R)).corner
    simp only [remap_lay, remap_cells] at h
    simpa only [ite_lay, ite_cells] using h

theorem relayout_lay_self : (O.relayout g R).lay g = (O.lay g).map R := by
  simp [relayout]

theorem relayout_lay_other {f : F} (h : f ≠ g) : (O.relayout g R).lay f = O.lay f := by
  simp [relayout, h]

theorem relayout_cells_self : (O.relayout g R).cells g = (O.cells g).image (Cell.mv R) := by
  simp [relayout]

theorem relayout_cells_other {f : F} (h : f ≠ g) : (O.relayout g R).cells f = O.cells f := by
  simp [relayout, h]

/-- The matching motion of a halfedge of another face whose twin lies in `g`. -/
theorem relayout_layΔ_twin {x : Halfedge F} (hx : x.1 ≠ g) (hx' : (t.G.twin x).1 = g) :
    layΔ t (O.relayout g R).lay x = (layΔ t O.lay x).comp R.inv := by
  have hL : (O.relayout g R).lay = fun f => (O.lay f).map (rf g R f) :=
    funext fun f => O.ite_lay g R f
  rw [hL, layΔ_map]
  simp only [rf, hx, hx', ↓reduceIte, Rigid.id_comp]

/-- The matching motion of a halfedge of `g` whose twin lies in another face. -/
theorem relayout_layΔ_self {x : Halfedge F} (hx : x.1 = g) (hx' : (t.G.twin x).1 ≠ g) :
    layΔ t (O.relayout g R).lay x = R.comp (layΔ t O.lay x) := by
  have hL : (O.relayout g R).lay = fun f => (O.lay f).map (rf g R f) :=
    funext fun f => O.ite_lay g R f
  rw [hL, layΔ_map]
  simp only [rf, hx, hx', ↓reduceIte, rigid_inv_id, Rigid.comp_id]

/-- Relaying the face across `h₀` by the matching motion of `h₀` makes that motion the
identity. -/
theorem relayout_layΔ_eq_id {h₀ : Halfedge F} (hf : t.G.Flippable h₀) :
    layΔ t (O.relayout (t.G.twin h₀).1 (layΔ t O.lay h₀)).lay h₀ = Rigid.id := by
  rw [relayout_layΔ_twin O _ _ (Ne.symm hf) rfl, Rigid.comp_inv]

end Overlay

end DdgProofs.Delaunay.Surface
