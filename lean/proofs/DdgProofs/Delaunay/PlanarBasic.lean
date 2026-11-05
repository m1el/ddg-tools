import DdgProofs.TriangleIntegral.Basic

/-!
# Planar Delaunay triangulations, part 1: triangles, circles and the lifting map in coordinates

The coordinate geometry behind the planar theory of `DdgProofs/Delaunay/Planar*.lean`. Points are
pairs `ℝ × ℝ` (the plane model of `DdgProofs/TriangleIntegral/Tiling.lean`); orientation is the
signed double area `sda a b c = (b − a) × (c − a)` of `DdgProofs/Bff/Quality.lean`, squared
lengths are `sqDist2`. Nothing here uses the sup metric of the type `ℝ × ℝ`.

## Main statements

* **Triangles.** For `sda a b c ≠ 0` the convex hull of `a, b, c` is the set of points that lie, for
  every side, on the closed side of the opposite corner (`mem_hull3_iff`; the barycentric coordinates
  are `sda x b c / sda a b c`, … and `x = Σ λ_v v`, `sda_cramer`). A point with all three of them
  positive is interior (`mem_interior_hull3`), and an interior point has all three positive
  (`pos_of_mem_interior_hull3`).
* **Interiors and null sets.** An open set covered by two triangles and disjoint from both of
  their interiors is empty (`eq_empty_of_subset_union_hull`): the boundary of a convex set is a
  null set (`Convex.addHaar_frontier`), a nonempty open set is not.
* **The in-circle determinant.** `inCircle a b c d = Σ_cyc |d a|² sda(d, b, c)` is the
  `4 × 4` lifting determinant (`inCircle_eq_lift`); it is alternating in `a, b, c`.
  `InCircumcircle a b c d :⟺ sda(a, b, c) · inCircle a b c d > 0` (independent of the order of
  `a, b, c`, `inCircumcircle_swap₁₂`, `inCircumcircle_rotate`). For **every** point `o` equidistant
  from `a, b, c` the power of `d` is `|o d|² − |o a|² = −inCircle a b c d / sda a b c`
  (`power_eq_of_equidistant`); the circumcenter exists (`circumcenter`, `circumcenter_equidistant`),
  so `InCircumcircle a b c d` says exactly that `d` lies in the **open circumdisk**
  (`inCircumcircle_iff_lt`). In terms of the lifting `x ↦ (x, |x|²)`: `inCircle` is `sda(a,b,c)`
  times the height of `d̂` below the plane through `â, b̂, ĉ` (`inCircle_eq_bary`).
* **A non-Delaunay diamond is convex.** If `c` and `d` lie on opposite sides of `ab` and `d` is
  inside the circumcircle of `a b c`, then `a` and `b` lie on opposite sides of `cd`
  (`crossing_of_inCircumcircle`) — the planar form of `isConvexDiamond_of_not_isLocallyDelaunay`
  in `DdgProofs/Delaunay/Local.lean`; it rests on the identity `inCircle_mul_sda_eq`.
* **Monotonicity of the power across an edge** (`circleLift_le_of_not_inCircumcircle`): if the
  edge `ab` of `a b c` and `a b d` is locally Delaunay, then every point `p` strictly on `d`'s side
  of `ab` has a smaller power with respect to the circle of `a b d` than to that of `a b c`. This is
  the step of the walk in the proof of Delaunay's lemma; it rests on the three-term identity
  `inCircle_pluecker`.
* **The lifted volume of a diamond.** For a strictly convex quadrilateral,
  `Σ area(t) · Σ_{v ∈ t} |v|²` decreases by `sda(a,b,c) · inCircle a b c d / (2 |sda(a,b,c)|)` when
  the diagonal `ab` is replaced by `cd` (`liftQuad_flip`); this is positive exactly when `d` is
  inside the circumcircle of `a b c`.

## References

* B. Delaunay, *Sur la sphère vide*, Bull. Acad. Sci. URSS (1934), 793–800 (the empty-sphere
  criterion).
* C. L. Lawson, *Software for C¹ surface interpolation*, in: Mathematical Software III, Academic
  Press, 1977, 161–194 (the flip algorithm; the in-circle test of a convex quadrilateral).
* H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge University Press, 2001,
  Ch. 1–2 (the lifting map; the in-circle determinant; the power of a point).
* M. de Berg, O. Cheong, M. van Kreveld, M. Overmars, *Computational Geometry: Algorithms and
  Applications*, 3rd ed., Springer, 2008, Ch. 9 (legal edges, Lemma 9.4 / Theorem 9.6).

The bibliographic data were not checked against the works here; the statements are proved from
the definitions in this file, not taken from the works.
-/

noncomputable section

namespace DDG

namespace Delaunay

namespace Planar

open MeasureTheory Set Filter Topology

/-! ## Signed area -/

theorem sda_self₁₂ (a c : ℝ × ℝ) : sda a a c = 0 := by unfold sda; ring

theorem sda_self₂₃ (a b : ℝ × ℝ) : sda a b b = 0 := by unfold sda; ring

theorem sda_self₁₃ (a b : ℝ × ℝ) : sda a b a = 0 := by unfold sda; ring

/-- The four-point identity `sda a b c = sda x b c + sda a x c + sda a b x`. -/
theorem sda_four (a b c x : ℝ × ℝ) : sda a b c = sda x b c + sda a x c + sda a b x := by
  unfold sda; ring

/-- `sda a b c · x = sda x b c · a + sda a x c · b + sda a b x · c` (Cramer's rule), first
coordinate. -/
theorem sda_cramer_fst (a b c x : ℝ × ℝ) :
    sda a b c * x.1 = sda x b c * a.1 + sda a x c * b.1 + sda a b x * c.1 := by unfold sda; ring

theorem sda_cramer_snd (a b c x : ℝ × ℝ) :
    sda a b c * x.2 = sda x b c * a.2 + sda a x c * b.2 + sda a b x * c.2 := by unfold sda; ring

/-- **Cramer's rule.** `x = Σ λ_v v` with the barycentric coordinates `λ_a = sda x b c / D`, … -/
theorem sda_cramer {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) (x : ℝ × ℝ) :
    x = (sda x b c / sda a b c) • a + (sda a x c / sda a b c) • b + (sda a b x / sda a b c) • c := by
  ext
  · simp only [Prod.fst_add, Prod.smul_fst, smul_eq_mul]
    field_simp
    linear_combination sda_cramer_fst a b c x
  · simp only [Prod.snd_add, Prod.smul_snd, smul_eq_mul]
    field_simp
    linear_combination sda_cramer_snd a b c x

/-- `sda` is affine in its last argument. -/
theorem sda_combo₃ (a b x y : ℝ × ℝ) {s t : ℝ} (h : s + t = 1) :
    sda a b (s • x + t • y) = s * sda a b x + t * sda a b y := by
  unfold sda
  simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  have : s = 1 - t := by linarith
  subst this
  ring

/-- `sda` is affine in its first argument. -/
theorem sda_combo₁ (b c x y : ℝ × ℝ) {s t : ℝ} (h : s + t = 1) :
    sda (s • x + t • y) b c = s * sda x b c + t * sda y b c := by
  unfold sda
  simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  have : s = 1 - t := by linarith
  subst this
  ring

/-- `sda` is affine in its middle argument. -/
theorem sda_combo₂ (a c x y : ℝ × ℝ) {s t : ℝ} (h : s + t = 1) :
    sda a (s • x + t • y) c = s * sda a x c + t * sda a y c := by
  unfold sda
  simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  have : s = 1 - t := by linarith
  subst this
  ring

/-- Translating the last argument: `sda a b (y + w) = sda a b y + (b − a) × w`. -/
theorem sda_add₃ (a b y w : ℝ × ℝ) :
    sda a b (y + w) = sda a b y + ((b.1 - a.1) * w.2 - (b.2 - a.2) * w.1) := by
  unfold sda; simp only [Prod.fst_add, Prod.snd_add]; ring

theorem sda_ne_zero_left {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : a ≠ b := by
  rintro rfl; exact h (sda_self₁₂ _ _)

theorem sda_ne_zero_right {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : b ≠ c := by
  rintro rfl; exact h (sda_self₂₃ _ _)

theorem sda_ne_zero_outer {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : a ≠ c := by
  rintro rfl; exact h (sda_self₁₃ _ _)

theorem continuous_sda₁ (b c : ℝ × ℝ) : Continuous fun x => sda x b c := by
  unfold sda; fun_prop

theorem continuous_sda₂ (a c : ℝ × ℝ) : Continuous fun x => sda a x c := by
  unfold sda; fun_prop

theorem continuous_sda₃ (a b : ℝ × ℝ) : Continuous fun x => sda a b x := by
  unfold sda; fun_prop

/-- If `x` lies on the lines `ab` and `bc` of a nondegenerate triangle, it is `b`. -/
theorem eq_of_sda_eq_zero {a b c x : ℝ × ℝ} (h : sda a b c ≠ 0) (h1 : sda x b c = 0)
    (h3 : sda a b x = 0) : x = b := by
  have h2 : sda a x c = sda a b c := by have := sda_four a b c x; linarith
  rw [sda_cramer h x, h1, h2, h3, div_self h]
  simp

/-! ## The convex hull of three points -/

/-- `x` lies on the closed side of `bc` that contains `a`, etc.: the hull of three points with
`sda a b c ≠ 0`, written with signs (no division). -/
theorem mem_hull3_iff {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) {x : ℝ × ℝ} :
    x ∈ convexHull ℝ {a, b, c} ↔
      0 ≤ sda a b c * sda x b c ∧ 0 ≤ sda a b c * sda a x c ∧ 0 ≤ sda a b c * sda a b x := by
  constructor
  · intro hx
    refine convexHull_min (t := {y | 0 ≤ sda a b c * sda y b c ∧ 0 ≤ sda a b c * sda a y c ∧
      0 ≤ sda a b c * sda a b y}) ?_ ?_ hx
    · intro y hy
      simp only [mem_insert_iff, mem_singleton_iff] at hy
      rcases hy with rfl | rfl | rfl
      all_goals refine ⟨?_, ?_, ?_⟩ <;> simp [sda_self₁₂, sda_self₁₃, sda_self₂₃, mul_self_nonneg]
    · intro y hy z hz s t hs ht hst
      simp only [Set.mem_ofPred_eq] at hy hz ⊢
      rw [sda_combo₁ _ _ _ _ hst, sda_combo₂ _ _ _ _ hst, sda_combo₃ _ _ _ _ hst]
      refine ⟨?_, ?_, ?_⟩ <;> nlinarith [mul_nonneg hs hy.1, mul_nonneg ht hz.1,
        mul_nonneg hs hy.2.1, mul_nonneg ht hz.2.1, mul_nonneg hs hy.2.2, mul_nonneg ht hz.2.2]
  · rintro ⟨h1, h2, h3⟩
    have hD2 : 0 < sda a b c * sda a b c := mul_self_pos.2 h
    have key : ∀ y : ℝ, 0 ≤ sda a b c * y → 0 ≤ y / sda a b c := fun y hy => by
      have : y / sda a b c = (sda a b c * y) / (sda a b c * sda a b c) := by
        field_simp
      rw [this]
      exact div_nonneg hy hD2.le
    rw [sda_cramer h x]
    have hs : sda x b c / sda a b c + sda a x c / sda a b c + sda a b x / sda a b c = 1 := by
      rw [← add_div, ← add_div, ← sda_four, div_self h]
    have hc := convex_convexHull ℝ ({a, b, c} : Set (ℝ × ℝ))
    have ha : a ∈ convexHull ℝ ({a, b, c} : Set (ℝ × ℝ)) := subset_convexHull ℝ _ (by simp)
    have hb : b ∈ convexHull ℝ ({a, b, c} : Set (ℝ × ℝ)) := subset_convexHull ℝ _ (by simp)
    have hc' : c ∈ convexHull ℝ ({a, b, c} : Set (ℝ × ℝ)) := subset_convexHull ℝ _ (by simp)
    have := hc.sum_mem (t := Finset.univ) (w := ![sda x b c / sda a b c, sda a x c / sda a b c,
      sda a b x / sda a b c]) (z := ![a, b, c]) (by
        intro i _
        fin_cases i
        · exact key _ h1
        · exact key _ h2
        · exact key _ h3) (by simp [Fin.sum_univ_three, hs]) (by
        intro i _
        fin_cases i <;> simp [ha, hb, hc'])
    simpa [Fin.sum_univ_three, add_assoc] using this

/-- A point with three positive barycentric coordinates is interior. -/
theorem mem_interior_hull3 {a b c x : ℝ × ℝ} (h : sda a b c ≠ 0) (h1 : 0 < sda a b c * sda x b c)
    (h2 : 0 < sda a b c * sda a x c) (h3 : 0 < sda a b c * sda a b x) :
    x ∈ interior (convexHull ℝ {a, b, c}) := by
  have hU : IsOpen {y : ℝ × ℝ | 0 < sda a b c * sda y b c ∧ 0 < sda a b c * sda a y c ∧
      0 < sda a b c * sda a b y} := by
    exact (isOpen_lt continuous_const (continuous_const.mul (continuous_sda₁ b c))).inter
      ((isOpen_lt continuous_const (continuous_const.mul (continuous_sda₂ a c))).inter
      (isOpen_lt continuous_const (continuous_const.mul (continuous_sda₃ a b))))
  refine interior_maximal (fun y hy => ?_) hU ⟨h1, h2, h3⟩
  exact (mem_hull3_iff h).2 ⟨hy.1.le, hy.2.1.le, hy.2.2.le⟩

/-- Along `s ↦ x + s • w`, every property that holds near `x` holds for some small `s > 0`. -/
theorem exists_pos_near {P : ℝ × ℝ → Prop} {x : ℝ × ℝ} (hP : ∀ᶠ y in 𝓝 x, P y) (w : ℝ × ℝ) :
    ∀ᶠ s in 𝓝[>] (0 : ℝ), P (x + s • w) := by
  have hc : Tendsto (fun s : ℝ => x + s • w) (𝓝[>] (0 : ℝ)) (𝓝 x) := by
    have hcont : Continuous fun s : ℝ => x + s • w := by fun_prop
    have h0 := hcont.tendsto (0 : ℝ)
    simp only [zero_smul, add_zero] at h0
    exact h0.mono_left nhdsWithin_le_nhds
  exact hc.eventually hP

/-- On an interior point the third barycentric coordinate is positive. -/
theorem pos_of_mem_interior_hull3₃ {a b c x : ℝ × ℝ} (h : sda a b c ≠ 0)
    (hx : x ∈ interior (convexHull ℝ {a, b, c})) : 0 < sda a b c * sda a b x := by
  have hx' := (mem_hull3_iff h).1 (interior_subset hx)
  rcases hx'.2.2.lt_or_eq with hlt | heq
  · exact hlt
  · exfalso
    have hz : sda a b x = 0 := by
      rcases mul_eq_zero.1 heq.symm with h0 | h0
      · exact absurd h0 h
      · exact h0
    obtain ⟨s, hs, hspos⟩ := ((exists_pos_near (isOpen_interior.mem_nhds hx) (x - c)).and
      self_mem_nhdsWithin).exists
    have hmem := (mem_hull3_iff h).1 (interior_subset hs)
    have : sda a b (x + s • (x - c)) = -(s * sda a b c) := by
      have e : x + s • (x - c) = (1 + s) • x + (-s) • c := by
        rw [smul_sub, add_smul, one_smul, neg_smul, sub_eq_add_neg, add_assoc]
      rw [e, sda_combo₃ _ _ _ _ (by ring), hz]
      ring
    rw [this] at hmem
    have hD2 : 0 < sda a b c * sda a b c := mul_self_pos.2 h
    have hs' : (0 : ℝ) < s := hspos
    nlinarith [hmem.2.2]

/-- **On an interior point all three barycentric coordinates are positive.** -/
theorem pos_of_mem_interior_hull3 {a b c x : ℝ × ℝ} (h : sda a b c ≠ 0)
    (hx : x ∈ interior (convexHull ℝ {a, b, c})) :
    0 < sda a b c * sda x b c ∧ 0 < sda a b c * sda a x c ∧ 0 < sda a b c * sda a b x := by
  refine ⟨?_, ?_, pos_of_mem_interior_hull3₃ h hx⟩
  · have e : ({a, b, c} : Set (ℝ × ℝ)) = {b, c, a} := by
      rw [insert_comm, pair_comm, insert_comm]
    have h' : sda b c a ≠ 0 := by
      have : sda b c a = sda a b c := by unfold sda; ring
      rwa [this]
    have := pos_of_mem_interior_hull3₃ h' (by rwa [← e])
    have e1 : sda b c a = sda a b c := by unfold sda; ring
    have e2 : sda b c x = sda x b c := by unfold sda; ring
    rwa [e1, e2] at this
  · have e : ({a, b, c} : Set (ℝ × ℝ)) = {c, a, b} := by
      rw [pair_comm, insert_comm]
    have h' : sda c a b ≠ 0 := by
      have : sda c a b = sda a b c := by unfold sda; ring
      rwa [this]
    have := pos_of_mem_interior_hull3₃ h' (by rwa [← e])
    have e1 : sda c a b = sda a b c := by unfold sda; ring
    have e2 : sda c a x = sda a x c := by unfold sda; ring
    rwa [e1, e2] at this

theorem isClosed_hull3 (a b c : ℝ × ℝ) : IsClosed (convexHull ℝ ({a, b, c} : Set (ℝ × ℝ))) :=
  ((Set.toFinite _).isCompact_convexHull (𝕜 := ℝ)).isClosed

/-- The centroid of a nondegenerate triangle is interior; the interior is nonempty. -/
theorem interior_hull3_nonempty {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) :
    (interior (convexHull ℝ ({a, b, c} : Set (ℝ × ℝ)))).Nonempty := by
  refine ⟨(1 / 3 : ℝ) • a + (1 / 3 : ℝ) • b + (1 / 3 : ℝ) • c, mem_interior_hull3 h ?_ ?_ ?_⟩ <;>
  · have hD2 : 0 < sda a b c * sda a b c := mul_self_pos.2 h
    have e1 : sda ((1 / 3 : ℝ) • a + (1 / 3 : ℝ) • b + (1 / 3 : ℝ) • c) b c = sda a b c / 3 := by
      unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
      ring
    have e2 : sda a ((1 / 3 : ℝ) • a + (1 / 3 : ℝ) • b + (1 / 3 : ℝ) • c) c = sda a b c / 3 := by
      unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
      ring
    have e3 : sda a b ((1 / 3 : ℝ) • a + (1 / 3 : ℝ) • b + (1 / 3 : ℝ) • c) = sda a b c / 3 := by
      unfold sda; simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
      ring
    first | rw [e1] | rw [e2] | rw [e3]
    nlinarith

/-- A boundary point of a nondegenerate triangle is a limit of interior points. -/
theorem hull3_subset_closure_interior {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) :
    convexHull ℝ ({a, b, c} : Set (ℝ × ℝ)) ⊆ closure (interior (convexHull ℝ {a, b, c})) := by
  rw [(convex_convexHull ℝ _).closure_interior_eq_closure_of_nonempty_interior
    (interior_hull3_nonempty h)]
  exact subset_closure

/-! ## Interiors and null sets -/

/-- **An open set covered by two convex closed sets and disjoint from both interiors is empty.**
(The boundaries are null sets, `Convex.addHaar_frontier`, and a nonempty open set is not.) -/
theorem eq_empty_of_subset_union_hull {K₁ K₂ U : Set (ℝ × ℝ)} (hK₁ : Convex ℝ K₁)
    (hK₂ : Convex ℝ K₂) (hU : IsOpen U) (hsub : U ⊆ K₁ ∪ K₂) (h₁ : Disjoint U (interior K₁))
    (h₂ : Disjoint U (interior K₂)) : U = ∅ := by
  have hf : U ⊆ frontier K₁ ∪ frontier K₂ := by
    intro x hx
    rcases hsub hx with hx1 | hx2
    · exact Or.inl ⟨subset_closure hx1, fun hi => Set.disjoint_left.1 h₁ hx hi⟩
    · exact Or.inr ⟨subset_closure hx2, fun hi => Set.disjoint_left.1 h₂ hx hi⟩
  have h0 : volume U = 0 := measure_mono_null hf
    (measure_union_null (hK₁.addHaar_frontier volume) (hK₂.addHaar_frontier volume))
  exact (hU.measure_eq_zero_iff volume).1 h0

/-! ## The in-circle determinant -/

/-- The **in-circle determinant** `det [a − d, |a − d|²; b − d, |b − d|²; c − d, |c − d|²]`. -/
def inCircle (a b c d : ℝ × ℝ) : ℝ :=
  sqDist2 d a * sda d b c + sqDist2 d b * sda d c a + sqDist2 d c * sda d a b

/-- `d` lies strictly inside the circle through `a, b, c` (for `sda a b c ≠ 0`; the sign of
`sda a b c` makes the test independent of the orientation). -/
def InCircumcircle (a b c d : ℝ × ℝ) : Prop := 0 < sda a b c * inCircle a b c d

/-- The squared norm `|x|²`, the height of the lifting map `x ↦ (x, |x|²)`. -/
def nsq (x : ℝ × ℝ) : ℝ := x.1 ^ 2 + x.2 ^ 2

/-- **`inCircle` is the lifting determinant** `det [x_v, y_v, |v|², 1]_{v = a,b,c,d}`. -/
theorem inCircle_eq_lift (a b c d : ℝ × ℝ) :
    inCircle a b c d = nsq a * sda d b c + nsq b * sda a d c + nsq c * sda a b d - nsq d * sda a b c := by
  unfold inCircle nsq sqDist2 sda; ring

theorem inCircle_rotate (a b c d : ℝ × ℝ) : inCircle b c a d = inCircle a b c d := by
  unfold inCircle; ring

theorem inCircle_swap₁₂ (a b c d : ℝ × ℝ) : inCircle b a c d = -inCircle a b c d := by
  unfold inCircle sda sqDist2; ring

theorem inCircle_self₁ (a b c : ℝ × ℝ) : inCircle a b c a = 0 := by
  unfold inCircle sda sqDist2; ring

theorem inCircle_self₂ (a b c : ℝ × ℝ) : inCircle a b c b = 0 := by
  unfold inCircle sda sqDist2; ring

theorem inCircle_self₃ (a b c : ℝ × ℝ) : inCircle a b c c = 0 := by
  unfold inCircle sda sqDist2; ring

theorem inCircumcircle_rotate {a b c d : ℝ × ℝ} : InCircumcircle b c a d ↔ InCircumcircle a b c d := by
  unfold InCircumcircle; rw [sda_rotate, inCircle_rotate]

theorem inCircumcircle_swap₁₂ {a b c d : ℝ × ℝ} : InCircumcircle b a c d ↔ InCircumcircle a b c d := by
  unfold InCircumcircle; rw [sda_swap_left, inCircle_swap₁₂, neg_mul_neg]

theorem inCircumcircle_swap₂₃ {a b c d : ℝ × ℝ} : InCircumcircle a c b d ↔ InCircumcircle a b c d := by
  rw [← inCircumcircle_rotate, inCircumcircle_swap₁₂, inCircumcircle_rotate, inCircumcircle_rotate]

/-- A corner is never strictly inside its own circumcircle. -/
theorem not_inCircumcircle_self₁ (a b c : ℝ × ℝ) : ¬ InCircumcircle a b c a := by
  simp [InCircumcircle, inCircle_self₁]

theorem not_inCircumcircle_self₂ (a b c : ℝ × ℝ) : ¬ InCircumcircle a b c b := by
  simp [InCircumcircle, inCircle_self₂]

theorem not_inCircumcircle_self₃ (a b c : ℝ × ℝ) : ¬ InCircumcircle a b c c := by
  simp [InCircumcircle, inCircle_self₃]

/-- **The lifting picture.** With the barycentric coordinates `λ = (sda d b c, sda a d c,
sda a b d) / sda a b c` of `d`, `inCircle a b c d / sda a b c = Σ λ_v |v|² − |d|²`: the height of
the plane through the lifted corners above the lifted `d`. -/
theorem inCircle_eq_bary {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) (d : ℝ × ℝ) :
    inCircle a b c d / sda a b c =
      sda d b c / sda a b c * nsq a + sda a d c / sda a b c * nsq b + sda a b d / sda a b c * nsq c
        - nsq d := by
  rw [inCircle_eq_lift]
  field_simp

/-! ## Circumcircle and power of a point -/

/-- **The power of a point.** For every `o` equidistant from `a, b, c`,
`|o d|² − |o a|² = −inCircle a b c d / sda a b c`. -/
theorem power_eq_of_equidistant {a b c o : ℝ × ℝ} (h : sda a b c ≠ 0)
    (hb : sqDist2 o a = sqDist2 o b) (hc : sqDist2 o a = sqDist2 o c) (d : ℝ × ℝ) :
    sqDist2 o d - sqDist2 o a = -inCircle a b c d / sda a b c := by
  rw [eq_div_iff h]
  have : (sqDist2 o d - sqDist2 o a) * sda a b c + inCircle a b c d =
      sda a c d * (sqDist2 o a - sqDist2 o b) + sda a d b * (sqDist2 o a - sqDist2 o c) := by
    unfold sqDist2 inCircle sda sqDist2; ring
  linear_combination this + sda a c d * hb + sda a d b * hc

/-- The circumcenter `a + (|u|² v^⊥ − |v|² u^⊥) / (2 u × v)`, `u = b − a`, `v = c − a`. -/
def circumcenter (a b c : ℝ × ℝ) : ℝ × ℝ :=
  (a.1 + ((c.2 - a.2) * sqDist2 a b - (b.2 - a.2) * sqDist2 a c) / (2 * sda a b c),
   a.2 + ((b.1 - a.1) * sqDist2 a c - (c.1 - a.1) * sqDist2 a b) / (2 * sda a b c))

/-- The circumcenter is equidistant from the three corners. -/
theorem circumcenter_equidistant {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) :
    sqDist2 (circumcenter a b c) a = sqDist2 (circumcenter a b c) b ∧
      sqDist2 (circumcenter a b c) a = sqDist2 (circumcenter a b c) c := by
  have h2 : 2 * sda a b c ≠ 0 := mul_ne_zero two_ne_zero h
  constructor <;>
  · unfold circumcenter sqDist2
    simp only
    field_simp
    unfold sda
    ring

/-- **`InCircumcircle a b c d` says that `d` lies in the open circumdisk**: closer to the
circumcenter than the corners are. -/
theorem inCircumcircle_iff_lt {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) (d : ℝ × ℝ) :
    InCircumcircle a b c d ↔ sqDist2 (circumcenter a b c) d < sqDist2 (circumcenter a b c) a := by
  have e := power_eq_of_equidistant h (circumcenter_equidistant h).1
    (circumcenter_equidistant h).2 d
  have hD2 : 0 < sda a b c * sda a b c := mul_self_pos.2 h
  have hs : sqDist2 (circumcenter a b c) d < sqDist2 (circumcenter a b c) a ↔
      sqDist2 (circumcenter a b c) d - sqDist2 (circumcenter a b c) a < 0 := sub_neg.symm
  rw [InCircumcircle, hs, e, neg_div, neg_lt_zero]
  constructor
  · intro hp
    have : inCircle a b c d / sda a b c = (sda a b c * inCircle a b c d) / (sda a b c * sda a b c) := by
      field_simp
    rw [this]; exact div_pos hp hD2
  · intro hp
    have : sda a b c * inCircle a b c d = (inCircle a b c d / sda a b c) * (sda a b c * sda a b c) := by
      field_simp
    rw [this]; exact mul_pos hp hD2

/-! ## A non-Delaunay diamond is strictly convex -/

/-- The identity behind `crossing_of_inCircumcircle`: with `α = sda d b c`, `γ = sda a b d`,
`sda a b c · inCircle a b c d = −|α (a − b) − γ (b − c)|² + sda a b c (α |ab|² + γ |bc|²)`. -/
theorem inCircle_mul_sda_eq (a b c d : ℝ × ℝ) :
    sda a b c * inCircle a b c d =
      -((sda d b c * (a.1 - b.1) - sda a b d * (b.1 - c.1)) ^ 2 +
        (sda d b c * (a.2 - b.2) - sda a b d * (b.2 - c.2)) ^ 2) +
      sda a b c * (sda d b c * sqDist2 a b + sda a b d * sqDist2 b c) := by
  unfold inCircle sda sqDist2; ring

/-- `d` inside the circumcircle of `a b c` and across `ab` from `c`: then `d` is on `a`'s side of
`bc` (strictly). -/
theorem sda_pos_of_inCircumcircle {a b c d : ℝ × ℝ} (hk : sda a b c * sda a b d < 0)
    (hin : InCircumcircle a b c d) : 0 < sda a b c * sda d b c := by
  by_contra hle
  push Not at hle
  have e := inCircle_mul_sda_eq a b c d
  unfold InCircumcircle at hin
  have h1 := sqDist2_nonneg a b
  have h2 := sqDist2_nonneg b c
  nlinarith [sq_nonneg (sda d b c * (a.1 - b.1) - sda a b d * (b.1 - c.1)),
    sq_nonneg (sda d b c * (a.2 - b.2) - sda a b d * (b.2 - c.2)), mul_nonpos_of_nonneg_of_nonpos h1 hle,
    mul_nonneg h2 (neg_nonneg.2 hk.le)]

/-- **A diamond whose edge is not locally Delaunay is strictly convex.** If `c, d` lie on opposite
sides of `ab` and `d` is strictly inside the circumcircle of `a b c`, then `a, b` lie on opposite
sides of `cd`: the segments `ab` and `cd` cross, and the edge can be flipped. -/
theorem crossing_of_inCircumcircle {a b c d : ℝ × ℝ} (hk : sda a b c * sda a b d < 0)
    (hin : InCircumcircle a b c d) : sda c d a * sda c d b < 0 := by
  have hα := sda_pos_of_inCircumcircle hk hin
  have hk' : sda b a c * sda b a d < 0 := by
    have : sda b a c * sda b a d = sda a b c * sda a b d := by unfold sda; ring
    linarith
  have hβ := sda_pos_of_inCircumcircle hk' (inCircumcircle_swap₁₂.2 hin)
  have e1 : sda c d a = -sda a d c := by unfold sda; ring
  have e2 : sda c d b = sda d b c := by unfold sda; ring
  have e3 : sda b a c * sda d a c = sda a b c * sda a d c := by unfold sda; ring
  rw [e3] at hβ
  rw [e1, e2]
  have hD : sda a b c ≠ 0 := by rintro h0; simp [h0] at hk
  have hD2 : 0 < sda a b c * sda a b c := mul_self_pos.2 hD
  nlinarith [mul_pos hα hβ]

/-! ## The power of a point across a locally Delaunay edge -/

/-- The three-term (Plücker-type) identity
`inCircle a b c p · sda a b d − inCircle a b d p · sda a b c = inCircle a b c d · sda a b p`. -/
theorem inCircle_pluecker (a b c d p : ℝ × ℝ) :
    inCircle a b c p * sda a b d - inCircle a b d p * sda a b c = inCircle a b c d * sda a b p := by
  unfold inCircle sda sqDist2; ring

/-- `inCircle a b c p / sda a b c`, the negative power of `p` with respect to the circumcircle
of `a, b, c` (`power_eq_of_equidistant`); positive exactly inside. -/
def circleLift (a b c p : ℝ × ℝ) : ℝ := inCircle a b c p / sda a b c

theorem circleLift_rotate (a b c p : ℝ × ℝ) : circleLift b c a p = circleLift a b c p := by
  unfold circleLift; rw [inCircle_rotate, sda_rotate]

theorem circleLift_swap₁₂ (a b c p : ℝ × ℝ) : circleLift b a c p = circleLift a b c p := by
  unfold circleLift; rw [inCircle_swap₁₂, sda_swap_left, neg_div_neg_eq]

theorem circleLift_pos_iff {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) (p : ℝ × ℝ) :
    0 < circleLift a b c p ↔ InCircumcircle a b c p := by
  have hD2 : 0 < sda a b c * sda a b c := mul_self_pos.2 h
  unfold circleLift InCircumcircle
  have : inCircle a b c p / sda a b c = (sda a b c * inCircle a b c p) / (sda a b c * sda a b c) := by
    field_simp
  rw [this, div_pos_iff_of_pos_right hD2]

/-- **Monotonicity across a locally Delaunay edge.** Let `a b c` and `a b d` lie on opposite sides
of `ab`, with `d` not strictly inside the circumcircle of `a b c`. Then for every `p` strictly on
`d`'s side of `ab`, the circle of `a b d` has a larger `circleLift` at `p` (a smaller power) than
the circle of `a b c`. In the lifting picture: the two lifted planes agree over the line `ab`, and
beyond it the plane of `a b d` is the higher one. -/
theorem circleLift_le_of_not_inCircumcircle {a b c d p : ℝ × ℝ} (hk : sda a b c * sda a b d < 0)
    (hloc : ¬ InCircumcircle a b c d) (hp : 0 < sda a b d * sda a b p) :
    circleLift a b c p ≤ circleLift a b d p := by
  unfold InCircumcircle at hloc
  push Not at hloc
  have hD1 : sda a b c ≠ 0 := by rintro h0; simp [h0] at hk
  have hD2 : sda a b d ≠ 0 := by rintro h0; simp [h0] at hk
  have e := inCircle_pluecker a b c d p
  have hprod : 0 ≤ inCircle a b c d * sda a b p := by
    have h1 : (sda a b c * inCircle a b c d) * (sda a b d * sda a b p) ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg hloc hp.le
    have h2 : (sda a b c * sda a b d) * (inCircle a b c d * sda a b p) ≤ 0 := by linarith
    by_contra hneg
    push Not at hneg
    nlinarith [mul_pos_of_neg_of_neg hk hneg]
  unfold circleLift
  rw [← sub_nonneg]
  have : inCircle a b d p / sda a b d - inCircle a b c p / sda a b c =
      -(inCircle a b c d * sda a b p) / (sda a b c * sda a b d) := by
    field_simp
    linear_combination -e
  rw [this]
  exact div_nonneg_of_nonpos (neg_nonpos.2 hprod) hk.le

/-! ## The lifted volume of a convex quadrilateral -/

/-- `area · Σ |v|²` for the triangle `a b c`: three times the volume under the lifted triangle
(the integral of the linear interpolant of `|x|²`). -/
def liftTri (a b c : ℝ × ℝ) : ℝ := |sda a b c| / 2 * (nsq a + nsq b + nsq c)

/-- Two reals of opposite signs whose difference `y − x` is positive: `x < 0 < y`. -/
theorem neg_pos_of_mul_neg_of_sub_pos {x y : ℝ} (h : x * y < 0) (hd : 0 < y - x) : x < 0 ∧ 0 < y := by
  rcases lt_trichotomy x 0 with hx | hx | hx
  · rcases lt_trichotomy y 0 with hy | hy | hy
    · nlinarith
    · subst hy; simp at h
    · exact ⟨hx, hy⟩
  · subst hx; simp at h
  · nlinarith

/-- The orientations of a strictly convex quadrilateral `a, d, b, c` (the segments `ab` and `cd`
cross): if `sda a b c > 0` then `sda a b d < 0`, `sda c d a < 0 < sda c d b`. -/
theorem quad_signs {a b c d : ℝ × ℝ} (hk : sda a b c * sda a b d < 0)
    (hc : sda c d a * sda c d b < 0) (hD : 0 < sda a b c) :
    sda a b d < 0 ∧ sda c d a < 0 ∧ 0 < sda c d b := by
  have h1 : sda a b d < 0 := by
    by_contra h; push Not at h; nlinarith [mul_nonneg hD.le h]
  have e : sda c d b - sda c d a = sda a b c - sda a b d := by unfold sda; ring
  exact ⟨h1, neg_pos_of_mul_neg_of_sub_pos hc (by linarith)⟩

/-- **The lifted volume of a convex quadrilateral under the two diagonals.** For a strictly convex
quadrilateral with diagonals `ab` and `cd`, replacing the triangles `a b c`, `a b d` by `a c d`,
`b c d` changes `Σ area · Σ|v|²` by `−sda(a,b,c) · inCircle a b c d / (2 |sda(a,b,c)|)`. This
is positive (a strict decrease) exactly when `d` is strictly inside the circumcircle of `a b c`. -/
theorem liftQuad_flip {a b c d : ℝ × ℝ} (hk : sda a b c * sda a b d < 0)
    (hc : sda c d a * sda c d b < 0) :
    liftTri a b c + liftTri a b d - (liftTri a c d + liftTri b c d) =
      sda a b c * inCircle a b c d / (2 * |sda a b c|) := by
  have hD : sda a b c ≠ 0 := by rintro h0; simp [h0] at hk
  have r1 : sda a c d = sda c d a := by unfold sda; ring
  have r2 : sda b c d = sda c d b := by unfold sda; ring
  unfold liftTri
  rw [r1, r2]
  rcases hD.lt_or_gt with hneg | hpos
  · have hk' : sda b a c * sda b a d < 0 := by
      have : sda b a c * sda b a d = sda a b c * sda a b d := by unfold sda; ring
      linarith
    have hc' : sda c d b * sda c d a < 0 := by linarith [mul_comm (sda c d a) (sda c d b)]
    have hpos' : 0 < sda b a c := by rw [sda_swap_left]; linarith
    obtain ⟨h1, h2, h3⟩ := quad_signs hk' hc' hpos'
    rw [sda_swap_left] at h1
    rw [abs_of_neg hneg, abs_of_pos (by linarith : 0 < sda a b d), abs_of_neg h2, abs_of_pos h3]
    field_simp
    unfold inCircle nsq sda sqDist2
    ring
  · obtain ⟨h1, h2, h3⟩ := quad_signs hk hc hpos
    rw [abs_of_pos hpos, abs_of_neg h1, abs_of_neg h2, abs_of_pos h3]
    field_simp
    unfold inCircle nsq sda sqDist2
    ring

end Planar

end Delaunay

end DDG
