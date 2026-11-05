import DDG.Delaunay.Surface.Overlay
import DDG.Delaunay.Surface.Count

/-!
# Delaunay layer 4: a straight segment of bounded length crosses boundedly many edges

For a chain of cells along a straight segment (`Chain`, `Overlay.lean`) the number of
crossings is bounded by a function `N s L` of the length `L` of the segment and of three
constants of `s` only: the smallest corner angle `αMin`, the smallest sine of a corner angle
`sinMin`, and the shortest edge `ℓMin`.

Every crossing but the first is a **chord** through a cell from one open side to another; the
two sides meet at a corner of the cell, the chord's **pivot** (`pivotIdx`). The angle at the pivot
between the two end points of the chord is the corner angle of `s` there (`angleAt_pivot`),
at least `αMin`. Consecutive chords with the same pivot (as a point of the plane) are a fan
around it, and a straight line through a fan subtends less than `π` (`sum_angleAt_ray_lt_pi`,
`Angle.lean`): a run of chords with one pivot has fewer than `π / αMin` members. When the pivot
changes between two consecutive chords, the two pivots are the two ends of the crossed edge, so
the two chords together are at least `sinMin · ℓMin` long (`dist2_ge_mul_sin`). The counting
lemma `fan_count` (`Count.lean`) turns this into `numCross ≤ 1 + M (1 + 2L / c₀)` with
`M = ⌈π/αMin⌉`, `c₀ = sinMin · ℓMin` (`Chain.numCross_le`).
-/

noncomputable section

namespace DDG.Delaunay.Surface

open DDG.Delaunay.Planar
open Set Real
open scoped Classical

variable {F : Type*}

/-! ## The constants of `s` -/

section Constants

variable [Fintype F] [Nonempty F] (s : MetricGluing F)

/-- The smallest corner angle of `s`. -/
def αMin : ℝ := Finset.univ.inf' Finset.univ_nonempty fun h : Halfedge F => Gluing.cornerAngleAt s.ℓ h

/-- The smallest sine of a corner angle of `s`. -/
def sinMin : ℝ := Finset.univ.inf' Finset.univ_nonempty fun h : Halfedge F =>
  sin (Gluing.cornerAngleAt s.ℓ h)

/-- The shortest edge of `s`. -/
def ℓMin : ℝ := Finset.univ.inf' Finset.univ_nonempty fun h : Halfedge F => s.ℓ h

theorem αMin_le (h : Halfedge F) : αMin s ≤ Gluing.cornerAngleAt s.ℓ h :=
  Finset.inf'_le _ (Finset.mem_univ h)

theorem sinMin_le (h : Halfedge F) : sinMin s ≤ sin (Gluing.cornerAngleAt s.ℓ h) :=
  Finset.inf'_le _ (Finset.mem_univ h)

theorem ℓMin_le (h : Halfedge F) : ℓMin s ≤ s.ℓ h := Finset.inf'_le _ (Finset.mem_univ h)

theorem αMin_pos : 0 < αMin s := by
  obtain ⟨h, -, e⟩ := Finset.exists_mem_eq_inf' (Finset.univ_nonempty (α := Halfedge F))
    fun h : Halfedge F => Gluing.cornerAngleAt s.ℓ h
  rw [αMin, e]
  exact s.isMetric.cornerAngleAt_pos h

theorem sinMin_pos : 0 < sinMin s := by
  obtain ⟨h, -, e⟩ := Finset.exists_mem_eq_inf' (Finset.univ_nonempty (α := Halfedge F))
    fun h : Halfedge F => sin (Gluing.cornerAngleAt s.ℓ h)
  rw [sinMin, e]
  exact sin_pos_of_pos_of_lt_pi (s.isMetric.cornerAngleAt_pos h)
    (cornerAngle_lt_pi (s.isMetric.tri_at h).rotate.swap)

theorem ℓMin_pos : 0 < ℓMin s := by
  obtain ⟨h, -, e⟩ := Finset.exists_mem_eq_inf' (Finset.univ_nonempty (α := Halfedge F))
    fun h : Halfedge F => s.ℓ h
  rw [ℓMin, e]
  exact (s.isMetric.tri_at h).pos_a

/-- The maximal number of chords around one pivot. -/
def blockBound : ℕ := ⌈π / αMin s⌉₊

/-- The minimal length of two consecutive chords with different pivots. -/
def switchCost : ℝ := sinMin s * ℓMin s

theorem switchCost_pos : 0 < switchCost s := mul_pos (sinMin_pos s) (ℓMin_pos s)

/-- **The crossing bound**: a segment of length `L` crosses at most `N s L` edges of `s`. -/
def crossBound (L : ℝ) : ℕ := 1 + blockBound s * (1 + ⌈2 * L / switchCost s⌉₊)

end Constants

/-! ## Chords of a chain -/

/-- The common corner of the sides `kin` and `kout` (`kin ≠ kout`) of a triangle. -/
def pivotIdx (kin kout : Fin 3) : Fin 3 := if kout = kin + 1 then kin + 1 else kin

theorem pivotIdx_mem {kin kout : Fin 3} (h : kin ≠ kout) :
    pivotIdx kin kout = kout ∨ pivotIdx kin kout = kout + 1 := by
  unfold pivotIdx
  split_ifs with h1
  · exact Or.inl h1.symm
  · right; omega

/-- The pivot is also a corner of the entry side. -/
theorem pivotIdx_mem' (kin kout : Fin 3) :
    pivotIdx kin kout = kin ∨ pivotIdx kin kout = kin + 1 := by
  unfold pivotIdx
  split_ifs
  · exact Or.inr rfl
  · exact Or.inl rfl

/-- A point of an open segment `(A, B)` is `A + a • (B − A)` with `0 < a < 1`. -/
theorem exists_eq_add_smul_of_mem_openSegment {A B z : Pt} (hz : z ∈ openSegment ℝ A B) :
    ∃ a : ℝ, 0 < a ∧ a < 1 ∧ z = A + a • (B - A) := by
  rw [openSegment_eq_image] at hz
  obtain ⟨θ, ⟨h0, h1⟩, rfl⟩ := hz
  refine ⟨θ, h0, h1, ?_⟩
  ext <;> simp <;> ring

/-- The angle at the pivot between a point of the entry side and a point of the exit side is
the corner angle of `s` at the pivot. -/
theorem angleAt_pivot {s : MetricGluing F} (c : Cell F) {kin kout : Fin 3} (h : kin ≠ kout)
    {z₀ z₁ : Pt} (h₀ : z₀ ∈ (c.tri s).side kin) (h₁ : z₁ ∈ (c.tri s).side kout) :
    angleAt (c.tri s (pivotIdx kin kout)) z₀ z₁ =
      Gluing.cornerAngleAt s.ℓ (c.1, pivotIdx kin kout) := by
  unfold pivotIdx
  split_ifs with hk
  · subst hk
    rw [Tri.side_symm] at h₀
    obtain ⟨a, ha, -, rfl⟩ := exists_eq_add_smul_of_mem_openSegment h₀
    obtain ⟨b, hb, -, rfl⟩ := exists_eq_add_smul_of_mem_openSegment h₁
    rw [angleAt_smul _ _ _ ha hb, angleAt_comm]
    have := angleAt_tri (s := s) c (kin + 1)
    rwa [show kin + 1 + 2 = kin by omega] at this
  · have hk2 : kout = kin + 2 := by omega
    subst hk2
    rw [Tri.side_symm, show kin + 2 + 1 = kin by omega] at h₁
    obtain ⟨a, ha, -, rfl⟩ := exists_eq_add_smul_of_mem_openSegment h₀
    obtain ⟨b, hb, -, rfl⟩ := exists_eq_add_smul_of_mem_openSegment h₁
    rw [angleAt_smul _ _ _ ha hb]
    exact angleAt_tri (s := s) c kin

/-- Scaling a point along its ray from `v` scales its distance from `v`. -/
theorem dist2_add_smul (v P : Pt) {a : ℝ} (ha : 0 ≤ a) :
    dist2 v (v + a • (P - v)) = a * dist2 v P := by
  unfold dist2
  rw [show (v + a • (P - v)).1 - v.1 = a * (P.1 - v.1) by simp,
    show (v + a • (P - v)).2 - v.2 = a * (P.2 - v.2) by simp,
    show (a * (P.1 - v.1)) ^ 2 + (a * (P.2 - v.2)) ^ 2 =
      a ^ 2 * ((P.1 - v.1) ^ 2 + (P.2 - v.2) ^ 2) by ring,
    Real.sqrt_mul (sq_nonneg a), Real.sqrt_sq ha]

/-- The distances from the ends of a segment to a point of it add up. -/
theorem dist2_add_of_mem_openSegment {A B z : Pt} (hz : z ∈ openSegment ℝ A B) :
    dist2 A z + dist2 z B = dist2 A B := by
  obtain ⟨a, ha, ha1, rfl⟩ := exists_eq_add_smul_of_mem_openSegment hz
  have e : A + a • (B - A) = B + (1 - a) • (A - B) := by ext <;> simp <;> ring
  rw [dist2_add_smul _ _ ha.le, e, dist2_comm (B + (1 - a) • (A - B)) B,
    dist2_add_smul _ _ (by linarith), dist2_comm B A]
  ring

theorem dist2_ray_ray (p d : Pt) (a b : ℝ) :
    dist2 (ray p d a) (ray p d b) = |b - a| * dist2 p (p + d) := by
  unfold dist2 ray
  rw [show (p + b • d).1 - (p + a • d).1 = (b - a) * d.1 by simp; ring,
    show (p + b • d).2 - (p + a • d).2 = (b - a) * d.2 by simp; ring,
    show ((b - a) * d.1) ^ 2 + ((b - a) * d.2) ^ 2 = (b - a) ^ 2 * (d.1 ^ 2 + d.2 ^ 2) by ring,
    Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq_eq_abs]
  congr 2; simp

/-- A point on the line of a ray is a point of the ray. -/
theorem exists_eq_ray_of_cross_eq_zero {p d v : Pt} (hd : d ≠ 0) (h : cross2 d (v - p) = 0) :
    ∃ u, v = ray p d u := by
  simp only [cross2, Prod.fst_sub, Prod.snd_sub] at h
  by_cases h1 : d.1 = 0
  · have h2 : d.2 ≠ 0 := fun h2 => hd (Prod.ext h1 h2)
    refine ⟨(v.2 - p.2) / d.2, ?_⟩
    ext
    · simp only [ray, Prod.fst_add, Prod.smul_fst, smul_eq_mul, h1, mul_zero] at h ⊢
      have : v.1 - p.1 = 0 := by
        have := mul_eq_zero.1 (show d.2 * (v.1 - p.1) = 0 by linarith)
        exact this.resolve_left h2
      linarith
    · simp only [ray, Prod.snd_add, Prod.smul_snd, smul_eq_mul]
      field_simp
      ring
  · refine ⟨(v.1 - p.1) / d.1, ?_⟩
    ext
    · simp only [ray, Prod.fst_add, Prod.smul_fst, smul_eq_mul]
      field_simp
      ring
    · simp only [ray, Prod.snd_add, Prod.smul_snd, smul_eq_mul]
      field_simp
      linarith

/-- The distances from the two ends of a segment, in either order, to a point of it add up. -/
theorem dist2_add_of_pair {A B z P₀ P₁ : Pt} (hz : z ∈ openSegment ℝ A B) (h₀ : P₀ = A ∨ P₀ = B)
    (h₁ : P₁ = B ∨ P₁ = A) (hne : P₀ ≠ P₁) : dist2 P₀ z + dist2 P₁ z = dist2 A B := by
  have hs := dist2_add_of_mem_openSegment hz
  rcases h₀ with h₀ | h₀ <;> rcases h₁ with h₁ | h₁
  · rw [h₀, h₁, dist2_comm B z]; exact hs
  · exact absurd (h₀.trans h₁.symm) hne
  · exact absurd (h₀.trans h₁.symm) hne
  · rw [h₀, h₁, add_comm, dist2_comm B z]; exact hs

/-- A chord of the cell `c` along the ray: entered at `t₀` through the open side `kin`, left at
`t₁` through the open side `kout`, interior in between. -/
structure IsChord (s : MetricGluing F) (p d : Pt) (c : Cell F) (t₀ t₁ : ℝ) (kin kout : Fin 3) :
    Prop where
  lt : t₀ < t₁
  int : ∀ u ∈ Ioo t₀ t₁, ray p d u ∈ (c.tri s).int
  entry : ray p d t₀ ∈ (c.tri s).side kin
  exit : ray p d t₁ ∈ (c.tri s).side kout

namespace IsChord

variable {s : MetricGluing F} {p d : Pt} {c : Cell F} {t₀ t₁ : ℝ} {kin kout : Fin 3}

/-- A chord enters and leaves through different sides. -/
theorem ne (h : IsChord s p d c t₀ t₁ kin kout) : kin ≠ kout := by
  rintro rfl
  have hm : (t₀ + t₁) / 2 ∈ Ioo t₀ t₁ := ⟨by linarith [h.lt], by linarith [h.lt]⟩
  exact Tri.side_ne_of_int_between c.tri_pos h.entry h.exit h.lt hm (h.int _ hm)

/-- No corner of the cell is on the line of a chord. -/
theorem cross_ne (h : IsChord s p d c t₀ t₁ kin kout) (hd : d ≠ 0) (j : Fin 3) :
    cross2 d (c.tri s j - p) ≠ 0 := by
  intro hc
  obtain ⟨u, hu⟩ := exists_eq_ray_of_cross_eq_zero hd hc
  have hT := c.tri_pos (s := s)
  have hmem := Tri.corner_mem_hull (T := c.tri s) j
  rw [hu] at hmem
  rcases lt_trichotomy u t₀ with hlt | rfl | hgt
  · exact Tri.not_mem_hull_before_side hT h.entry
      (Tri.entry_side_cross_pos hT h.lt h.int h.entry) hlt hmem
  · exact Tri.ne_corner_of_mem_side hT h.entry j hu.symm
  · rcases lt_trichotomy u t₁ with hlt' | rfl | hgt'
    · exact Tri.corner_not_mem_int hT j (hu ▸ h.int u ⟨hgt, hlt'⟩)
    · exact Tri.ne_corner_of_mem_side hT h.exit j hu.symm
    · exact Tri.not_mem_hull_after_side hT h.exit
        (Tri.exit_side_cross_neg hT h.lt h.int h.exit) hgt' hmem

/-- The angle of a chord at its pivot is the corner angle of `s` there. -/
theorem angle (h : IsChord s p d c t₀ t₁ kin kout) :
    angleAt (c.tri s (pivotIdx kin kout)) (ray p d t₀) (ray p d t₁) =
      Gluing.cornerAngleAt s.ℓ (c.1, pivotIdx kin kout) :=
  angleAt_pivot c h.ne h.entry h.exit

/-- The chord bound towards the exit point. -/
theorem exit_le [Fintype F] [Nonempty F] (h : IsChord s p d c t₀ t₁ kin kout) (hd : d ≠ 0) :
    sinMin s * dist2 (c.tri s (pivotIdx kin kout)) (ray p d t₁) ≤
      dist2 (ray p d t₀) (ray p d t₁) := by
  have hs : sda (c.tri s (pivotIdx kin kout)) (ray p d t₀) (ray p d t₁) ≠ 0 := by
    rw [sda_ray]; exact mul_ne_zero (sub_ne_zero.2 h.lt.ne') (h.cross_ne hd _)
  have h1 := dist2_ge_mul_sin hs
  rw [h.angle] at h1
  have h2 := sinMin_le s (c.1, pivotIdx kin kout)
  nlinarith [DDG.dist2_nonneg (c.tri s (pivotIdx kin kout)) (ray p d t₁)]

/-- The chord bound towards the entry point. -/
theorem entry_le [Fintype F] [Nonempty F] (h : IsChord s p d c t₀ t₁ kin kout) (hd : d ≠ 0) :
    sinMin s * dist2 (c.tri s (pivotIdx kin kout)) (ray p d t₀) ≤
      dist2 (ray p d t₀) (ray p d t₁) := by
  have hs : sda (c.tri s (pivotIdx kin kout)) (ray p d t₁) (ray p d t₀) ≠ 0 := by
    rw [sda_ray]; exact mul_ne_zero (sub_ne_zero.2 h.lt.ne) (h.cross_ne hd _)
  have h1 := dist2_ge_mul_sin hs
  rw [angleAt_comm, h.angle, dist2_comm (ray p d t₁)] at h1
  have h2 := sinMin_le s (c.1, pivotIdx kin kout)
  nlinarith [DDG.dist2_nonneg (c.tri s (pivotIdx kin kout)) (ray p d t₀)]

end IsChord

/-- The records of the steps of a chain: the cell, its entry and exit parameters, its exit
side. -/
structure Step (F : Type*) where
  c : Cell F
  t₀ : ℝ
  t₁ : ℝ
  k : Fin 3

namespace Overlay.Chain

variable {s : MetricGluing F} {p d : Pt}

/-- The steps of a chain. -/
def steps : {c : Cell F} → {t₀ : ℝ} → Chain s p d c t₀ → List (Step F)
  | _, _, last _ _ _ _ _ => []
  | _, _, step c t₀ t₁ k _ _ _ _ rest => ⟨c, t₀, t₁, k⟩ :: rest.steps

theorem steps_length {c : Cell F} {t₀ : ℝ} (C : Chain s p d c t₀) :
    C.steps.length = C.numCross := by
  induction C with
  | last => rfl
  | step _ _ _ _ _ _ _ _ rest ih => simp only [steps, numCross, exits, List.length_cons]; rw [ih]; rfl

/-- The facts recorded about every step. -/
theorem steps_forall {c : Cell F} {t₀ : ℝ} (C : Chain s p d c t₀) :
    ∀ st ∈ C.steps, st.t₀ < st.t₁ ∧ st.t₁ < 1 ∧
      (∀ u ∈ Ioo st.t₀ st.t₁, ray p d u ∈ (st.c.tri s).int) ∧
      ray p d st.t₁ ∈ (st.c.tri s).side st.k := by
  induction C with
  | last => simp [steps]
  | step c t₀ t₁ k hlt h₁ hint hside rest ih =>
    intro st hst
    simp only [steps, List.mem_cons] at hst
    rcases hst with rfl | hst
    · exact ⟨hlt, h₁, hint, hside⟩
    · exact ih st hst

/-- Consecutive steps are linked: the next cell is the unfolding across the exit side, entered
at the exit parameter. -/
theorem steps_isChain {c : Cell F} {t₀ : ℝ} (C : Chain s p d c t₀) :
    List.IsChain (fun a b : Step F => b.c = a.c.across s a.k ∧ b.t₀ = a.t₁) C.steps := by
  induction C with
  | last => simp [steps]
  | step c t₀ t₁ k hlt h₁ hint hside rest ih =>
    simp only [steps]
    cases rest with
    | last c' t₀' j hq hS => simp [steps]
    | step c' t₀' t₁' k' _ _ _ _ rest' =>
      simp only [steps, List.isChain_cons_cons] at ih ⊢
      exact ⟨by simp, ih⟩

theorem steps_head {c : Cell F} {t₀ : ℝ} (C : Chain s p d c t₀) :
    ∀ st ∈ C.steps.head?, st.c = c ∧ st.t₀ = t₀ := by
  cases C with
  | last => simp [steps]
  | step => simp [steps]

/-- **The crossing bound for a chain**: along a segment of length `L = dist2 p (ray p d 1)` from
the parameter `0`, with `d ≠ 0`. -/
theorem numCross_le [Fintype F] [Nonempty F] {c : Cell F} (C : Chain s p d c 0) (hd : d ≠ 0) :
    C.numCross ≤ crossBound s (dist2 p (ray p d 1)) := by
  have hF := steps_forall C
  have hL := steps_isChain C
  have hH := steps_head C
  rw [← steps_length]
  obtain ⟨S, hS⟩ : ∃ S, S = C.steps := ⟨_, rfl⟩
  rw [← hS] at hF hL hH ⊢
  by_cases hN : S.length = 0
  · rw [hN]; exact Nat.zero_le _
  obtain ⟨n, hn⟩ := Nat.exists_eq_succ_of_ne_zero hN
  obtain ⟨st, hst⟩ : ∃ st : ℕ → Step F, ∀ j (hj : j < S.length), st j = S[j] :=
    ⟨fun j => S.getD j S[0], fun j hj => List.getD_eq_getElem _ _ hj⟩
  have hF' : ∀ j, j < n + 1 → (st j).t₀ < (st j).t₁ ∧ (st j).t₁ < 1 ∧
      (∀ u ∈ Ioo (st j).t₀ (st j).t₁, ray p d u ∈ ((st j).c.tri s).int) ∧
      ray p d (st j).t₁ ∈ ((st j).c.tri s).side (st j).k := by
    intro j hj
    rw [hst j (by omega)]
    exact hF _ (List.getElem_mem _)
  have hlink : ∀ j, j < n → (st (j + 1)).c = (st j).c.across s (st j).k ∧
      (st (j + 1)).t₀ = (st j).t₁ := by
    intro j hj
    rw [hst j (by omega), hst (j + 1) (by omega)]
    exact hL.getElem j (by omega)
  have hhead : (st 0).t₀ = 0 := by
    rw [hst 0 (by omega)]
    refine (hH S[0] ?_).2
    simp [List.head?_eq_getElem?]
  obtain ⟨kinN, hkinN⟩ : ∃ f : ℕ → Fin 3, ∀ j, f j = (s.G.twin ((st j).c.1, (st j).k)).2 :=
    ⟨_, fun _ => rfl⟩
  have hchord : ∀ j, j < n →
      IsChord s p d (st (j + 1)).c (st (j + 1)).t₀ (st (j + 1)).t₁ (kinN j) (st (j + 1)).k := by
    intro j hj
    obtain ⟨h1, -, h3, h4⟩ := hF' (j + 1) (by omega)
    obtain ⟨e1, e2⟩ := hlink j hj
    refine ⟨h1, h3, ?_, h4⟩
    rw [e2, e1, hkinN, Cell.across_side]
    exact (hF' j (by omega)).2.2.2
  obtain ⟨vN, hvN⟩ : ∃ f : ℕ → Pt, ∀ j,
      f j = (st (j + 1)).c.tri s (pivotIdx (kinN j) (st (j + 1)).k) := ⟨_, fun _ => rfl⟩
  obtain ⟨lamN, hlamN⟩ : ∃ f : ℕ → ℝ, ∀ j,
      f j = dist2 (ray p d (st (j + 1)).t₀) (ray p d (st (j + 1)).t₁) := ⟨_, fun _ => rfl⟩
  have hlam : ∀ i : Fin n, 0 ≤ lamN i.1 := fun i => by
    rw [hlamN]; exact DDG.dist2_nonneg _ _
  -- the fan bound
  have hrun : ∀ i j : Fin n, i ≤ j → (∀ m, i ≤ m → m ≤ j → vN m.1 = vN i.1) →
      j.1 + 1 - i.1 ≤ blockBound s := by
    intro i j hij hv
    have hij' : i.1 ≤ j.1 := hij
    have hjn : j.1 < n := j.2
    obtain ⟨m, hm⟩ : ∃ m, m = j.1 + 1 - i.1 := ⟨_, rfl⟩
    rw [← hm]
    let t : Fin (m + 1) → ℝ := fun q => (st (i.1 + q.1)).t₁
    have hmono : StrictMono t := by
      refine Fin.strictMono_iff_lt_succ.2 fun q => ?_
      show (st (i.1 + q.1)).t₁ < (st (i.1 + (q.1 + 1))).t₁
      have hq := q.2
      have e := (hlink (i.1 + q.1) (by omega)).2
      have l := (hF' (i.1 + q.1 + 1) (by omega)).1
      rw [show i.1 + (q.1 + 1) = i.1 + q.1 + 1 by ring]
      linarith
    have hvp : cross2 d (vN i.1 - p) ≠ 0 := by
      rw [hvN]; exact (hchord i.1 i.2).cross_ne hd _
    have hlt := sum_angleAt_ray_lt_pi hvp m t hmono
    have hge : ∀ q : Fin m,
        αMin s ≤ angleAt (vN i.1) (ray p d (t q.castSucc)) (ray p d (t q.succ)) := by
      intro q
      have hq := q.2
      have hmq : i.1 + q.1 < n := by omega
      have hvq : vN (i.1 + q.1) = vN i.1 :=
        hv ⟨i.1 + q.1, hmq⟩ (by simp only [Fin.le_def]; omega) (by simp only [Fin.le_def]; omega)
      rw [← hvq]
      show αMin s ≤ angleAt (vN (i.1 + q.1)) (ray p d (st (i.1 + q.1)).t₁)
        (ray p d (st (i.1 + (q.1 + 1))).t₁)
      rw [hvN, ← (hlink (i.1 + q.1) hmq).2, show i.1 + (q.1 + 1) = i.1 + q.1 + 1 by ring,
        (hchord _ hmq).angle]
      exact αMin_le s _
    have hsum : (m : ℝ) * αMin s ≤
        ∑ q : Fin m, angleAt (vN i.1) (ray p d (t q.castSucc)) (ray p d (t q.succ)) := by
      have := Finset.sum_le_sum (fun q (_ : q ∈ Finset.univ) => hge q)
      simpa using this
    have hlt2 : (m : ℝ) < π / αMin s := by
      rw [lt_div_iff₀ (αMin_pos s)]; linarith
    have hc := Nat.le_ceil (π / αMin s)
    have : (m : ℝ) < blockBound s := by unfold blockBound; linarith
    exact_mod_cast this.le
  -- the switch cost
  have hsw : ∀ i : Fin n, ∀ h : i.1 + 1 < n, vN i.1 ≠ vN (i.1 + 1) →
      switchCost s ≤ lamN i.1 + lamN (i.1 + 1) := by
    intro i hi hne
    have ch0 := hchord i.1 i.2
    have ch1 := hchord (i.1 + 1) hi
    obtain ⟨e1, e2⟩ := hlink (i.1 + 1) hi
    have hz : ray p d (st (i.1 + 1)).t₁ ∈ openSegment ℝ ((st (i.1 + 1)).c.tri s (st (i.1 + 1)).k)
        ((st (i.1 + 1)).c.tri s ((st (i.1 + 1)).k + 1)) := ch0.exit
    have hAB := Cell.tri_dist (s := s) (st (i.1 + 1)).c (st (i.1 + 1)).k
    have hℓ := ℓMin_le s ((st (i.1 + 1)).c.1, (st (i.1 + 1)).k)
    have l0 := ch0.exit_le hd
    have l1 := ch1.entry_le hd
    rw [e2] at l1
    rw [hvN, hvN] at hne
    rw [hlamN, hlamN]
    have hv0 : (st (i.1 + 1)).c.tri s (pivotIdx (kinN i.1) (st (i.1 + 1)).k) =
          (st (i.1 + 1)).c.tri s (st (i.1 + 1)).k ∨
        (st (i.1 + 1)).c.tri s (pivotIdx (kinN i.1) (st (i.1 + 1)).k) =
          (st (i.1 + 1)).c.tri s ((st (i.1 + 1)).k + 1) := by
      rcases pivotIdx_mem ch0.ne with h | h
      · left; rw [h]
      · right; rw [h]
    have hv1 : (st (i.1 + 1 + 1)).c.tri s (pivotIdx (kinN (i.1 + 1)) (st (i.1 + 1 + 1)).k) =
          (st (i.1 + 1)).c.tri s ((st (i.1 + 1)).k + 1) ∨
        (st (i.1 + 1 + 1)).c.tri s (pivotIdx (kinN (i.1 + 1)) (st (i.1 + 1 + 1)).k) =
          (st (i.1 + 1)).c.tri s (st (i.1 + 1)).k := by
      rcases pivotIdx_mem' (kinN (i.1 + 1)) (st (i.1 + 1 + 1)).k with h | h
      · left; rw [h, e1, hkinN]; exact Cell.across_tail _ _
      · right; rw [h, e1, hkinN]; exact Cell.across_head _ _
    have hsum := dist2_add_of_pair hz hv0 hv1 hne
    rw [hAB] at hsum
    rw [e2]
    have hs := sinMin_pos s
    have h3 := mul_le_mul_of_nonneg_left hℓ hs.le
    have h4 := congrArg (sinMin s * ·) hsum
    simp only [mul_add] at h4
    unfold switchCost
    linarith
  have hfc := fan_count (fun i : Fin n => vN i.1) (fun i : Fin n => lamN i.1) (switchCost_pos s)
    hlam hrun hsw
  -- the total length of the chords
  have hD := DDG.dist2_nonneg p (p + d)
  have hΛ : ∑ i : Fin n, lamN i.1 ≤ dist2 p (ray p d 1) := by
    rw [Fin.sum_univ_eq_sum_range (fun j => lamN j) n, ray_one]
    have e : ∀ j ∈ Finset.range n,
        lamN j = ((st (j + 1)).t₁ - (st j).t₁) * dist2 p (p + d) := by
      intro j hj
      rw [Finset.mem_range] at hj
      rw [hlamN, dist2_ray_ray, abs_of_pos (sub_pos.2 (hchord j hj).lt), (hlink j hj).2]
    rw [Finset.sum_congr rfl e, ← Finset.sum_mul,
      Finset.sum_range_sub (fun j => (st j).t₁) n]
    have h1 := (hF' n (by omega)).2.1
    have h2 := (hF' 0 (by omega)).1
    rw [hhead] at h2
    nlinarith
  have hΛ0 : 0 ≤ ∑ i : Fin n, lamN i.1 := Finset.sum_nonneg fun i _ => hlam i
  have hc := switchCost_pos s
  have h1 : 2 * (∑ i : Fin n, lamN i.1) / switchCost s ≤
      (⌈2 * dist2 p (ray p d 1) / switchCost s⌉₊ : ℝ) :=
    (div_le_div_of_nonneg_right (by linarith) hc.le).trans (Nat.le_ceil _)
  have h2 : (n : ℝ) ≤
      ((blockBound s * (1 + ⌈2 * dist2 p (ray p d 1) / switchCost s⌉₊) : ℕ) : ℝ) := by
    push_cast
    refine hfc.trans ?_
    have : (0 : ℝ) ≤ blockBound s := Nat.cast_nonneg _
    nlinarith
  have h3 : n ≤ blockBound s * (1 + ⌈2 * dist2 p (ray p d 1) / switchCost s⌉₊) := by
    exact_mod_cast h2
  unfold crossBound
  omega

end Overlay.Chain

/-- **Every edge of `t` is the developed length of a word of `s` with at most
`crossBound s (t.ℓ x)` crossings.** -/
theorem Overlay.exists_word_le [Fintype F] [Nonempty F] {s t : MetricGluing F} (O : Overlay s t)
    (x : Halfedge F) :
    ∃ w : Word F, w.numCross ≤ crossBound s (t.ℓ x) ∧ w.devLength s = t.ℓ x := by
  obtain ⟨c, j₀, C, -, hw⟩ := O.exists_word x.1 x.2
  refine ⟨C.word j₀, ?_, hw⟩
  have := C.numCross_le O.rayDir_ne_zero
  rw [O.ray_rayDir_one, O.lay_dist] at this
  exact this

end DDG.Delaunay.Surface
