import Mathlib
import DDG.Geometry.Advancement

/-!
# Closest points of segments and triangles (Ericson 2004, §5.1.2, §5.1.5, §5.1.9)

The distances that `crates/ddg-repulsive-curves/src/ccd.rs` computes for the conservative
advancement, and the closest point of a triangle of `obstacle.rs` and of
`crates/ddg-points-as-tori/src/truth.rs`, are the true infimum distances (exact real arithmetic,
any real inner product space; the cross products of the Rust code are rewritten by Lagrange's
identity in `ℝ³`, see the end of the file).

* `psDist p a b` (`point_segment_distance`): `‖(p − a) − t (b − a)‖` with `t` the clamped
  projection parameter; `psDist_le`: it is at most the distance to every point of the segment,
  and it is attained (`psParam_mem`).
* `segSegDist` (the minimum of the four endpoint distances and of the interior pair of the lines
  when it lies in `[0,1]²`), `segSegDist_le`, `segSegDist_attained`: the distance of two segments.
  `segSegCode` is the fixed Rust routine (Ericson's re-projection added): `segSegCode_eq`.
* `segSegThr ε` is the former Rust routine (interior pair only when `ae − b² > ε·ae`): it is never
  below the distance, and `segSegThr_sq_le`: `segSegThr² ≤ d² + (ε/4)·‖q₁ − p₁‖²` (survey
  finding 1; the constant is attained by two segments crossing at their midpoints).
* `ptTriDist` (`point_triangle_distance`: the three edge distances and, when the projection onto
  the plane lies in the triangle, the plane distance): `ptTriDist_le`, `ptTriDist_attained`; it
  holds for degenerate (collinear) triangles too.
* `closestPtTri` (Ericson's Voronoi-region routine of `truth.rs`, with its degenerate fallback):
  `closestPtTri_isClosest`; `closestPtTriRC` (the copy in `obstacle.rs`, no guards) agrees with it
  for non-degenerate triangles (`closestPtTriRC_eq`).
* `primDist_segSeg`, `primDist_ptTri`: these are `primDist` of `Advancement.lean` for the corner positions
  at time `τ`; `segSeg_advance_safe`, `ptTri_advance_safe`: the advancement of `ccd.rs` with the
  computed distances keeps the primitives apart (with corner speeds measured relative to any
  common velocity, `primDist_sub_common`).
-/

noncomputable section

namespace DDG

namespace RC

open Set

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

local notation "⟪" x ", " y "⟫" => @inner ℝ _ _ x y

/-- `x.clamp(0.0, 1.0)` of the Rust code. -/
def clamp01 (x : ℝ) : ℝ := max 0 (min x 1)

theorem clamp01_mem (x : ℝ) : clamp01 x ∈ Icc (0 : ℝ) 1 :=
  ⟨le_max_left _ _, max_le zero_le_one (min_le_right _ _)⟩

theorem clamp01_of_mem {x : ℝ} (hx : x ∈ Icc (0 : ℝ) 1) : clamp01 x = x := by
  unfold clamp01; rw [min_eq_left hx.2, max_eq_right hx.1]

/-- The clamped value is the point of `[0, 1]` nearest to `x`. -/
theorem clamp01_sq_le (x : ℝ) {s : ℝ} (hs : s ∈ Icc (0 : ℝ) 1) :
    (clamp01 x - x) ^ 2 ≤ (s - x) ^ 2 := by
  obtain ⟨h0, h1⟩ := hs
  unfold clamp01
  rcases le_total x 0 with hx | hx
  · rw [min_eq_left (hx.trans zero_le_one), max_eq_left hx]; nlinarith
  rcases le_total x 1 with hx1 | hx1
  · rw [min_eq_left hx1, max_eq_right hx]; nlinarith
  · rw [min_eq_right hx1, max_eq_right zero_le_one]; nlinarith

/-! ## Point–segment -/

/-- The point `a + s (b − a)` of the segment `a b`. -/
def segPt (a b : E) (s : ℝ) : E := a + s • (b - a)

/-- The parameter of `point_segment_distance`: `clamp(⟪p − a, b − a⟫/e)` if `e = ‖b − a‖² > 0`,
else `0`. -/
def psParam (p a b : E) : ℝ :=
  if 0 < ⟪b - a, b - a⟫ then clamp01 (⟪p - a, b - a⟫ / ⟪b - a, b - a⟫) else 0

/-- `point_segment_distance(p, a, b)`. -/
def psDist (p a b : E) : ℝ := ‖(p - a) - psParam p a b • (b - a)‖

theorem psParam_mem (p a b : E) : psParam p a b ∈ Icc (0 : ℝ) 1 := by
  unfold psParam; split_ifs
  · exact clamp01_mem _
  · exact ⟨le_rfl, zero_le_one⟩

theorem psDist_eq (p a b : E) : psDist p a b = ‖p - segPt a b (psParam p a b)‖ := by
  unfold psDist segPt; congr 1; abel

/-- `‖x − s y‖² = ‖y‖² (s − s₀)² + const` with `s₀ = ⟪x, y⟫/‖y‖²`. -/
theorem norm_sub_smul_sq (x y : E) (s : ℝ) :
    ‖x - s • y‖ ^ 2 = ‖x‖ ^ 2 - 2 * s * ⟪x, y⟫ + s ^ 2 * ⟪y, y⟫ := by
  rw [norm_sub_sq_real, inner_smul_right, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs,
    real_inner_self_eq_norm_sq]
  ring

theorem norm_le_of_sq_le {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) (h : x ^ 2 ≤ y ^ 2) : x ≤ y :=
  (sq_le_sq₀ hx hy).1 h

/-- **Point–segment distance**: `psDist p a b` is at most `‖p − x‖` for every point `x` of the
segment.

[formula: rc-closest-points] -/
theorem psDist_le (p a b : E) {s : ℝ} (hs : s ∈ Icc (0 : ℝ) 1) :
    psDist p a b ≤ ‖p - segPt a b s‖ := by
  have hx : p - segPt a b s = (p - a) - s • (b - a) := by unfold segPt; abel
  rw [hx]
  unfold psDist psParam
  split_ifs with he
  · refine norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) ?_
    rw [norm_sub_smul_sq, norm_sub_smul_sq]
    set e := ⟪b - a, b - a⟫
    set k := ⟪p - a, b - a⟫
    have h := clamp01_sq_le (k / e) hs
    have e1 : ∀ r : ℝ, ‖p - a‖ ^ 2 - 2 * r * k + r ^ 2 * e =
        ‖p - a‖ ^ 2 + e * ((r - k / e) ^ 2 - (k / e) ^ 2) := fun r => by
      field_simp; ring
    rw [e1, e1]
    nlinarith [mul_le_mul_of_nonneg_left h he.le]
  · have he0 : ⟪b - a, b - a⟫ = 0 := le_antisymm (not_lt.1 he) real_inner_self_nonneg
    have hba : b - a = 0 := by rwa [inner_self_eq_zero (𝕜 := ℝ)] at he0
    simp [hba]

/-! ## The exit lemma: a path that leaves a region crosses its boundary -/

/-- A continuous function that is `≥ 0` at `0` and `< 0` at `L ≥ 0` has a zero in `[0, L]`. -/
theorem exists_zero_of_nonneg_of_neg {g : ℝ → ℝ} (hg : Continuous g) {L : ℝ} (hL : 0 ≤ L)
    (h0 : 0 ≤ g 0) (hL' : g L < 0) : ∃ l ∈ Icc 0 L, g l = 0 := by
  have := intermediate_value_Icc' hL hg.continuousOn ⟨hL'.le, h0⟩
  obtain ⟨l, hl, hgl⟩ := this
  exact ⟨l, hl, hgl⟩

/-- The unit square `[0, 1]²` of the two segment parameters. -/
abbrev InSq (s t : ℝ) : Prop := s ∈ Icc (0 : ℝ) 1 ∧ t ∈ Icc (0 : ℝ) 1

/-- A point of the square on its boundary. -/
def OnSqBdry (s t : ℝ) : Prop := s = 0 ∨ s = 1 ∨ t = 0 ∨ t = 1

/-- A segment of parameters from a point of the square to a point outside it (or on its
boundary) meets the boundary. -/
theorem exists_sq_bdry {s₀ t₀ u v L : ℝ} (hL : 0 ≤ L) (h₀ : InSq s₀ t₀)
    (hend : ¬ InSq (s₀ + L * u) (t₀ + L * v) ∨
      (InSq (s₀ + L * u) (t₀ + L * v) ∧ OnSqBdry (s₀ + L * u) (t₀ + L * v))) :
    ∃ l ∈ Icc 0 L, InSq (s₀ + l * u) (t₀ + l * v) ∧ OnSqBdry (s₀ + l * u) (t₀ + l * v) := by
  rcases hend with hout | hin
  · let g : ℝ → ℝ := fun l =>
      min (min (s₀ + l * u) (1 - (s₀ + l * u))) (min (t₀ + l * v) (1 - (t₀ + l * v)))
    have hg : Continuous g := by fun_prop
    have hg0 : 0 ≤ g 0 := by
      obtain ⟨⟨a1, a2⟩, ⟨b1, b2⟩⟩ := h₀
      simp only [g, zero_mul, add_zero]
      exact le_min (le_min a1 (by linarith)) (le_min b1 (by linarith))
    have hgL : g L < 0 := by
      simp only [g]
      by_contra hc
      push Not at hc
      apply hout
      have h1 := (le_min_iff.1 hc).1
      have h2 := (le_min_iff.1 hc).2
      exact ⟨⟨(le_min_iff.1 h1).1, by linarith [(le_min_iff.1 h1).2]⟩,
        ⟨(le_min_iff.1 h2).1, by linarith [(le_min_iff.1 h2).2]⟩⟩
    obtain ⟨l, hl, hgl⟩ := exists_zero_of_nonneg_of_neg hg hL hg0 hgL
    refine ⟨l, hl, ?_, ?_⟩
    · have hc : 0 ≤ g l := hgl.ge
      simp only [g] at hc
      have h1 := (le_min_iff.1 hc).1
      have h2 := (le_min_iff.1 hc).2
      exact ⟨⟨(le_min_iff.1 h1).1, by linarith [(le_min_iff.1 h1).2]⟩,
        ⟨(le_min_iff.1 h2).1, by linarith [(le_min_iff.1 h2).2]⟩⟩
    · simp only [g] at hgl
      unfold OnSqBdry
      rcases min_choice (min (s₀ + l * u) (1 - (s₀ + l * u)))
          (min (t₀ + l * v) (1 - (t₀ + l * v))) with h | h <;> rw [h] at hgl
      · rcases min_choice (s₀ + l * u) (1 - (s₀ + l * u)) with h' | h' <;> rw [h'] at hgl
        · exact Or.inl hgl
        · exact Or.inr (Or.inl (by linarith))
      · rcases min_choice (t₀ + l * v) (1 - (t₀ + l * v)) with h' | h' <;> rw [h'] at hgl
        · exact Or.inr (Or.inr (Or.inl hgl))
        · exact Or.inr (Or.inr (Or.inr (by linarith)))
  · exact ⟨L, ⟨hL, le_rfl⟩, hin⟩

/-! ## Segment–segment -/

section SegSeg

variable (p1 q1 p2 q2 : E)

/-- `X(s, t) = (p₁ + s d₁) − (p₂ + t d₂)`, `d₁ = q₁ − p₁`, `d₂ = q₂ − p₂`. -/
def ssVec (s t : ℝ) : E := segPt p1 q1 s - segPt p2 q2 t

/-- The scalars of the Rust code: `a = d₁·d₁`, `e = d₂·d₂`, `b = d₁·d₂`, `c = d₁·r`,
`f = d₂·r` with `r = p₁ − p₂`, `denom = ae − b²`. -/
def ssA : ℝ := ⟪q1 - p1, q1 - p1⟫
def ssE : ℝ := ⟪q2 - p2, q2 - p2⟫
def ssB : ℝ := ⟪q1 - p1, q2 - p2⟫
def ssC : ℝ := ⟪q1 - p1, p1 - p2⟫
def ssF : ℝ := ⟪q2 - p2, p1 - p2⟫
def ssDen : ℝ := ssA p1 q1 * ssE p2 q2 - ssB p1 q1 p2 q2 ^ 2

/-- The closest points of the two lines: `s* = (bf − ce)/denom`, `t* = (af − bc)/denom`. -/
def ssS : ℝ := (ssB p1 q1 p2 q2 * ssF p1 p2 q2 - ssC p1 q1 p2 * ssE p2 q2) / ssDen p1 q1 p2 q2
def ssT : ℝ := (ssA p1 q1 * ssF p1 p2 q2 - ssB p1 q1 p2 q2 * ssC p1 q1 p2) / ssDen p1 q1 p2 q2

/-- The minimum of the four endpoint–segment distances. -/
def ssBest : ℝ :=
  min (min (min (psDist p1 p2 q2) (psDist q1 p2 q2)) (psDist p2 p1 q1)) (psDist q2 p1 q1)

/-- The segment–segment distance as the minimum over the candidate pairs: the four endpoint
distances and, when `a, e, denom > 0` and the closest points of the lines lie in both
segments, the interior pair. -/
def segSegDist : ℝ :=
  if 0 < ssA p1 q1 ∧ 0 < ssE p2 q2 ∧ 0 < ssDen p1 q1 p2 q2 ∧
      InSq (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2) then
    min (ssBest p1 q1 p2 q2) ‖ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2)‖
  else ssBest p1 q1 p2 q2

/-- **`segment_segment_distance` of `ccd.rs`** (exact arithmetic): the endpoint distances and the
pair `s₀ = clamp(s*)`, `t₁ = clamp((b s₀ + f)/e)`, `s₂ = clamp((b t₁ − c)/a)`, when
`a, e, |d₁ × d₂|² > 0` (`|d₁ × d₂|² = denom` and the cross-product form of `s*` is `ssS`:
`cross_formulas_segSeg`).

[formula: rc-closest-points] -/
def segSegCode : ℝ :=
  if 0 < ssA p1 q1 ∧ 0 < ssE p2 q2 ∧ 0 < ssDen p1 q1 p2 q2 then
    min (ssBest p1 q1 p2 q2)
      ‖ssVec p1 q1 p2 q2
        (clamp01 ((ssB p1 q1 p2 q2 * clamp01 ((ssB p1 q1 p2 q2 * clamp01 (ssS p1 q1 p2 q2) +
            ssF p1 p2 q2) / ssE p2 q2) - ssC p1 q1 p2) / ssA p1 q1))
        (clamp01 ((ssB p1 q1 p2 q2 * clamp01 (ssS p1 q1 p2 q2) + ssF p1 p2 q2) / ssE p2 q2))‖
  else ssBest p1 q1 p2 q2

/-- **The former `segment_segment_distance`** (before the fix of survey finding 1): the interior
pair only when `denom > ε·a·e` (`ε = 1e-14`). -/
def segSegThr (ε : ℝ) : ℝ :=
  if 0 < ssA p1 q1 ∧ 0 < ssE p2 q2 ∧ ε * ssA p1 q1 * ssE p2 q2 < ssDen p1 q1 p2 q2 ∧
      InSq (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2) then
    min (ssBest p1 q1 p2 q2) ‖ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2)‖
  else ssBest p1 q1 p2 q2

theorem ssVec_eq (s t : ℝ) :
    ssVec p1 q1 p2 q2 s t = (p1 - p2) + s • (q1 - p1) - t • (q2 - p2) := by
  unfold ssVec segPt; abel

variable {p1 q1 p2 q2}

/-- The endpoint candidates bound every boundary pair of the square. -/
theorem ssBest_le_of_bdry {s t : ℝ} (h : InSq s t) (hb : OnSqBdry s t) :
    ssBest p1 q1 p2 q2 ≤ ‖ssVec p1 q1 p2 q2 s t‖ := by
  unfold ssBest
  rcases hb with rfl | rfl | rfl | rfl
  · have e : ssVec p1 q1 p2 q2 0 t = p1 - segPt p2 q2 t := by unfold ssVec segPt; simp
    rw [e]
    exact (min_le_left _ _).trans ((min_le_left _ _).trans
      ((min_le_left _ _).trans (psDist_le p1 p2 q2 h.2)))
  · have e : ssVec p1 q1 p2 q2 1 t = q1 - segPt p2 q2 t := by unfold ssVec segPt; simp
    rw [e]
    exact (min_le_left _ _).trans ((min_le_left _ _).trans
      ((min_le_right _ _).trans (psDist_le q1 p2 q2 h.2)))
  · have e : ‖ssVec p1 q1 p2 q2 s 0‖ = ‖p2 - segPt p1 q1 s‖ := by
      rw [← norm_neg]; unfold ssVec segPt; congr 1; simp
    rw [e]
    exact (min_le_left _ _).trans ((min_le_right _ _).trans (psDist_le p2 p1 q1 h.1))
  · have e : ‖ssVec p1 q1 p2 q2 s 1‖ = ‖q2 - segPt p1 q1 s‖ := by
      rw [← norm_neg]; unfold ssVec segPt; congr 1; simp
    rw [e]
    exact (min_le_right _ _).trans (psDist_le q2 p1 q1 h.1)

/-- Stationarity of the closest points of two lines, for abstract `r, d₁, d₂`. -/
theorem star_inner_aux (r d1 d2 : E)
    (hD : ⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2 ≠ 0) :
    ⟪r + ((⟪d1, d2⟫ * ⟪d2, r⟫ - ⟪d1, r⟫ * ⟪d2, d2⟫) / (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2)) • d1 -
        ((⟪d1, d1⟫ * ⟪d2, r⟫ - ⟪d1, d2⟫ * ⟪d1, r⟫) / (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2)) • d2,
      d1⟫ = 0 ∧
    ⟪r + ((⟪d1, d2⟫ * ⟪d2, r⟫ - ⟪d1, r⟫ * ⟪d2, d2⟫) / (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2)) • d1 -
        ((⟪d1, d1⟫ * ⟪d2, r⟫ - ⟪d1, d2⟫ * ⟪d1, r⟫) / (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2)) • d2,
      d2⟫ = 0 := by
  simp only [inner_sub_left, inner_add_left, inner_smul_left, RCLike.conj_to_real]
  rw [real_inner_comm d1 r, real_inner_comm d2 r, real_inner_comm d1 d2]
  generalize ⟪d1, d1⟫ = A at hD ⊢
  generalize ⟪d2, d2⟫ = C at hD ⊢
  generalize ⟪d1, d2⟫ = B at hD ⊢
  generalize ⟪d1, r⟫ = P
  generalize ⟪d2, r⟫ = Q
  have hS : (B * Q - P * C) / (A * C - B ^ 2) * (A * C - B ^ 2) = B * Q - P * C :=
    div_mul_cancel₀ _ hD
  have hT : (A * Q - B * P) / (A * C - B ^ 2) * (A * C - B ^ 2) = A * Q - B * P :=
    div_mul_cancel₀ _ hD
  constructor <;> refine mul_left_cancel₀ hD ?_ <;> rw [mul_zero]
  · linear_combination A * hS - B * hT
  · linear_combination B * hS - C * hT

/-- The closest points of the lines are stationary: `X(s*, t*) ⊥ d₁, d₂`. -/
theorem ssVec_star_inner (hD : ssDen p1 q1 p2 q2 ≠ 0) :
    ⟪ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2), q1 - p1⟫ = 0 ∧
      ⟪ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2), q2 - p2⟫ = 0 := by
  rw [ssVec_eq]
  exact star_inner_aux (p1 - p2) (q1 - p1) (q2 - p2) hD

/-- **The quadratic splits at the closest points of the lines**:
`‖X(s, t)‖² = ‖X(s*, t*)‖² + ‖(s − s*) d₁ − (t − t*) d₂‖²`. -/
theorem ssVec_sq_split (hD : ssDen p1 q1 p2 q2 ≠ 0) (s t : ℝ) :
    ‖ssVec p1 q1 p2 q2 s t‖ ^ 2 =
      ‖ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2)‖ ^ 2 +
        ‖(s - ssS p1 q1 p2 q2) • (q1 - p1) - (t - ssT p1 q1 p2 q2) • (q2 - p2)‖ ^ 2 := by
  obtain ⟨h1, h2⟩ := ssVec_star_inner hD
  have e : ssVec p1 q1 p2 q2 s t = ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2) +
      ((s - ssS p1 q1 p2 q2) • (q1 - p1) - (t - ssT p1 q1 p2 q2) • (q2 - p2)) := by
    rw [ssVec_eq, ssVec_eq]; module
  rw [e, norm_add_sq_real, inner_sub_right, inner_smul_right, inner_smul_right, h1, h2]
  ring

theorem ssA_nonneg : 0 ≤ ssA p1 q1 := real_inner_self_nonneg
theorem ssE_nonneg : 0 ≤ ssE p2 q2 := real_inner_self_nonneg

/-- `‖x‖² = ⟪x, x⟫`, expanded bilinearly: `‖a d₂ − b d₁‖² = a·(ae − b²)`. -/
theorem norm_sq_smul_sub_aux (d1 d2 : E) :
    ‖⟪d1, d1⟫ • d2 - ⟪d1, d2⟫ • d1‖ ^ 2 = ⟪d1, d1⟫ * (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2) := by
  rw [← real_inner_self_eq_norm_sq]
  simp only [inner_sub_left, inner_sub_right, inner_smul_left, inner_smul_right,
    RCLike.conj_to_real]
  rw [real_inner_comm d2 d1]; ring

/-- `‖a d₂ − b d₁‖² = a·denom`. -/
theorem norm_sq_ssA_smul_sub :
    ‖ssA p1 q1 • (q2 - p2) - ssB p1 q1 p2 q2 • (q1 - p1)‖ ^ 2 =
      ssA p1 q1 * ssDen p1 q1 p2 q2 :=
  norm_sq_smul_sub_aux (q1 - p1) (q2 - p2)

/-- `denom ≥ 0` (Cauchy–Schwarz). -/
theorem ssDen_nonneg : 0 ≤ ssDen p1 q1 p2 q2 := by
  unfold ssDen ssA ssE ssB
  nlinarith [real_inner_mul_inner_self_le (q1 - p1) (q2 - p2)]

/-- **Segment–segment distance, lower bound**: `segSegDist ≤ ‖X(s, t)‖` for all
`(s, t) ∈ [0, 1]²`.

[formula: rc-closest-points] -/
theorem segSegDist_le {s t : ℝ} (h : InSq s t) :
    segSegDist p1 q1 p2 q2 ≤ ‖ssVec p1 q1 p2 q2 s t‖ := by
  unfold segSegDist
  split_ifs with hc
  · obtain ⟨-, -, hD, -⟩ := hc
    refine (min_le_right _ _).trans (norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) ?_)
    rw [ssVec_sq_split hD.ne' s t]
    nlinarith [sq_nonneg ‖(s - ssS p1 q1 p2 q2) • (q1 - p1) - (t - ssT p1 q1 p2 q2) • (q2 - p2)‖]
  · by_cases ha : 0 < ssA p1 q1
    · by_cases he : 0 < ssE p2 q2
      · by_cases hD : 0 < ssDen p1 q1 p2 q2
        · -- the closest points of the lines are outside the square: walk towards them
          have hout : ¬ InSq (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2) := fun hi => hc ⟨ha, he, hD, hi⟩
          obtain ⟨l, hl, hin, hbd⟩ := exists_sq_bdry (s₀ := s) (t₀ := t)
            (u := ssS p1 q1 p2 q2 - s) (v := ssT p1 q1 p2 q2 - t) zero_le_one h
            (Or.inl (by simpa using hout))
          refine (ssBest_le_of_bdry hin hbd).trans
            (norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) ?_)
          rw [ssVec_sq_split hD.ne', ssVec_sq_split hD.ne' s t]
          have e : (s + l * (ssS p1 q1 p2 q2 - s) - ssS p1 q1 p2 q2) • (q1 - p1) -
              (t + l * (ssT p1 q1 p2 q2 - t) - ssT p1 q1 p2 q2) • (q2 - p2) =
              (1 - l) • ((s - ssS p1 q1 p2 q2) • (q1 - p1) -
                (t - ssT p1 q1 p2 q2) • (q2 - p2)) := by module
          rw [e, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs]
          have : (1 - l) ^ 2 ≤ 1 := by nlinarith [hl.1, hl.2]
          nlinarith [sq_nonneg ‖(s - ssS p1 q1 p2 q2) • (q1 - p1) -
            (t - ssT p1 q1 p2 q2) • (q2 - p2)‖]
        · -- parallel lines: slide along the direction (b, a), which does not change X
          have hD0 : ssDen p1 q1 p2 q2 = 0 := le_antisymm (not_lt.1 hD) ssDen_nonneg
          have hpar : ssA p1 q1 • (q2 - p2) - ssB p1 q1 p2 q2 • (q1 - p1) = 0 := by
            have := norm_sq_ssA_smul_sub (p1 := p1) (q1 := q1) (p2 := p2) (q2 := q2)
            rw [hD0, mul_zero] at this
            exact norm_eq_zero.1 (pow_eq_zero_iff two_ne_zero |>.1 this)
          have hconst : ∀ l : ℝ, ssVec p1 q1 p2 q2 (s + l * ssB p1 q1 p2 q2)
              (t + l * ssA p1 q1) = ssVec p1 q1 p2 q2 s t := fun l => by
            rw [ssVec_eq, ssVec_eq]
            have : (s + l * ssB p1 q1 p2 q2) • (q1 - p1) - (t + l * ssA p1 q1) • (q2 - p2) =
                s • (q1 - p1) - t • (q2 - p2) -
                  l • (ssA p1 q1 • (q2 - p2) - ssB p1 q1 p2 q2 • (q1 - p1)) := by module
            rw [add_sub_assoc, this, hpar, smul_zero, sub_zero, add_sub_assoc]
          have hL : 0 ≤ 2 / ssA p1 q1 := by positivity
          obtain ⟨l, -, hin, hbd⟩ := exists_sq_bdry (s₀ := s) (t₀ := t) (u := ssB p1 q1 p2 q2)
            (v := ssA p1 q1) hL h (Or.inl (fun hi => by
              have := hi.2.2
              rw [div_mul_cancel₀ _ ha.ne'] at this
              linarith [h.2.1]))
          rw [← hconst l]
          exact ssBest_le_of_bdry hin hbd
      · have he0 : ssE p2 q2 = 0 := le_antisymm (not_lt.1 he) ssE_nonneg
        have hd : q2 - p2 = 0 := by
          have h' : ⟪q2 - p2, q2 - p2⟫ = 0 := he0
          rwa [inner_self_eq_zero (𝕜 := ℝ)] at h'
        have e : ssVec p1 q1 p2 q2 s t = ssVec p1 q1 p2 q2 s 0 := by
          rw [ssVec_eq, ssVec_eq, hd]; simp
        rw [e]
        exact ssBest_le_of_bdry ⟨h.1, le_rfl, zero_le_one⟩ (Or.inr (Or.inr (Or.inl rfl)))
    · have ha0 : ssA p1 q1 = 0 := le_antisymm (not_lt.1 ha) ssA_nonneg
      have hd : q1 - p1 = 0 := by
        have h' : ⟪q1 - p1, q1 - p1⟫ = 0 := ha0
        rwa [inner_self_eq_zero (𝕜 := ℝ)] at h'
      have e : ssVec p1 q1 p2 q2 s t = ssVec p1 q1 p2 q2 0 t := by
        rw [ssVec_eq, ssVec_eq, hd]; simp
      rw [e]
      exact ssBest_le_of_bdry ⟨⟨le_rfl, zero_le_one⟩, h.2⟩ (Or.inl rfl)

/-- The endpoint candidate is the distance of a pair of the square. -/
theorem ssBest_attained :
    ∃ s t, InSq s t ∧ ssBest p1 q1 p2 q2 = ‖ssVec p1 q1 p2 q2 s t‖ := by
  have h1 : psDist p1 p2 q2 = ‖ssVec p1 q1 p2 q2 0 (psParam p1 p2 q2)‖ := by
    rw [psDist_eq]; unfold ssVec segPt; simp
  have h2 : psDist q1 p2 q2 = ‖ssVec p1 q1 p2 q2 1 (psParam q1 p2 q2)‖ := by
    rw [psDist_eq]; unfold ssVec segPt; simp
  have h3 : psDist p2 p1 q1 = ‖ssVec p1 q1 p2 q2 (psParam p2 p1 q1) 0‖ := by
    rw [psDist_eq, ← norm_neg]; unfold ssVec segPt; congr 1; simp
  have h4 : psDist q2 p1 q1 = ‖ssVec p1 q1 p2 q2 (psParam q2 p1 q1) 1‖ := by
    rw [psDist_eq, ← norm_neg]; unfold ssVec segPt; congr 1; simp
  have i0 : (0 : ℝ) ∈ Icc (0 : ℝ) 1 := ⟨le_rfl, zero_le_one⟩
  have i1 : (1 : ℝ) ∈ Icc (0 : ℝ) 1 := ⟨zero_le_one, le_rfl⟩
  unfold ssBest
  rcases min_choice (min (min (psDist p1 p2 q2) (psDist q1 p2 q2)) (psDist p2 p1 q1))
      (psDist q2 p1 q1) with h | h <;> rw [h]
  · rcases min_choice (min (psDist p1 p2 q2) (psDist q1 p2 q2)) (psDist p2 p1 q1) with h' | h' <;>
      rw [h']
    · rcases min_choice (psDist p1 p2 q2) (psDist q1 p2 q2) with h'' | h'' <;> rw [h'']
      · exact ⟨_, _, ⟨i0, psParam_mem _ _ _⟩, h1⟩
      · exact ⟨_, _, ⟨i1, psParam_mem _ _ _⟩, h2⟩
    · exact ⟨_, _, ⟨psParam_mem _ _ _, i0⟩, h3⟩
  · exact ⟨_, _, ⟨psParam_mem _ _ _, i1⟩, h4⟩

/-- **Segment–segment distance, attained**: `segSegDist` is the distance of a pair of points of
the two segments; with `segSegDist_le` it is the minimum distance.

[formula: rc-closest-points] -/
theorem segSegDist_attained :
    ∃ s t, InSq s t ∧ segSegDist p1 q1 p2 q2 = ‖ssVec p1 q1 p2 q2 s t‖ := by
  unfold segSegDist
  split_ifs with hc
  · rcases min_choice (ssBest p1 q1 p2 q2)
        ‖ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2)‖ with h | h <;> rw [h]
    · exact ssBest_attained
    · exact ⟨_, _, hc.2.2.2, rfl⟩
  · exact ssBest_attained

theorem ssBest_ge : segSegDist p1 q1 p2 q2 ≤ ssBest p1 q1 p2 q2 := by
  obtain ⟨s, t, h, e⟩ := ssBest_attained (p1 := p1) (q1 := q1) (p2 := p2) (q2 := q2)
  rw [e]; exact segSegDist_le h

theorem valley_t_aux (A E B P Q : ℝ) (he : 0 < E) (hD : 0 < A * E - B ^ 2) :
    (B * ((B * Q - P * E) / (A * E - B ^ 2)) + Q) / E = (A * Q - B * P) / (A * E - B ^ 2) := by
  generalize hDd : A * E - B ^ 2 = D at *
  have hD' := hD.ne'
  field_simp
  rw [← hDd]; ring

theorem valley_s_aux (A E B P Q : ℝ) (ha : 0 < A) (hD : 0 < A * E - B ^ 2) :
    (B * ((A * Q - B * P) / (A * E - B ^ 2)) - P) / A = (B * Q - P * E) / (A * E - B ^ 2) := by
  generalize hDd : A * E - B ^ 2 = D at *
  have hD' := hD.ne'
  field_simp
  rw [← hDd]; ring

/-- **The fixed Rust routine is exact**: `segSegCode = segSegDist` (exact arithmetic). The
re-projected pair is a pair of the square (so never below the distance), and it is the interior
pair `(s*, t*)` whenever that lies in the square.

[formula: rc-closest-points] -/
theorem segSegCode_eq : segSegCode p1 q1 p2 q2 = segSegDist p1 q1 p2 q2 := by
  apply le_antisymm
  · unfold segSegCode
    split_ifs with hc
    · by_cases hi : InSq (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2)
      · obtain ⟨ha, he, hD⟩ := hc
        have hS := hi.1
        have hT := hi.2
        have e1 : (ssB p1 q1 p2 q2 * ssS p1 q1 p2 q2 + ssF p1 p2 q2) / ssE p2 q2 =
            ssT p1 q1 p2 q2 := valley_t_aux _ _ _ _ _ he hD
        have e2 : (ssB p1 q1 p2 q2 * ssT p1 q1 p2 q2 - ssC p1 q1 p2) / ssA p1 q1 =
            ssS p1 q1 p2 q2 := valley_s_aux _ _ _ _ _ ha hD
        rw [clamp01_of_mem hS, e1, clamp01_of_mem hT, e2, clamp01_of_mem hS]
        unfold segSegDist; rw [ite_eq_left ⟨ha, he, hD, hi⟩]
      · unfold segSegDist; rw [ite_eq_right (fun h => hi h.2.2.2)]; exact min_le_left _ _
    · unfold segSegDist; rw [ite_eq_right (fun h => hc ⟨h.1, h.2.1, h.2.2.1⟩)]
  · unfold segSegCode
    split_ifs
    · exact le_min ssBest_ge (segSegDist_le ⟨clamp01_mem _, clamp01_mem _⟩)
    · exact ssBest_ge

/-! ### The former routine (survey finding 1) -/

/-- The former routine is never below the distance (every candidate is a pair of the
square). -/
theorem segSegDist_le_thr {ε : ℝ} (hε : 0 ≤ ε) :
    segSegDist p1 q1 p2 q2 ≤ segSegThr p1 q1 p2 q2 ε := by
  unfold segSegThr
  split_ifs with hc
  · obtain ⟨ha, he, hD, hi⟩ := hc
    have : 0 < ssDen p1 q1 p2 q2 :=
      lt_of_le_of_lt (mul_nonneg (mul_nonneg hε ssA_nonneg) ssE_nonneg) hD
    unfold segSegDist; rw [ite_eq_left ⟨ha, he, this, hi⟩]
  · exact ssBest_ge

theorem valley_aux (d1 d2 : E) (l : ℝ) (he : 0 < ⟪d2, d2⟫) :
    ‖l • (d1 - (⟪d1, d2⟫ / ⟪d2, d2⟫) • d2)‖ ^ 2 =
      l ^ 2 * ((⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2) / ⟪d2, d2⟫) := by
  rw [norm_smul, mul_pow, Real.norm_eq_abs, sq_abs, ← real_inner_self_eq_norm_sq]
  simp only [inner_sub_left, inner_sub_right, inner_smul_left, inner_smul_right,
    RCLike.conj_to_real]
  rw [real_inner_comm d1 d2]
  generalize ⟪d2, d2⟫ = C at *
  generalize ⟪d1, d2⟫ = B
  generalize ⟪d1, d1⟫ = A
  have := he.ne'
  field_simp
  ring

/-- Along the valley `t = t* + (b/e)(s − s*)` the squared distance grows by `(denom/e)(s − s*)²`. -/
theorem ssVec_sq_valley (hD : ssDen p1 q1 p2 q2 ≠ 0) (he : 0 < ssE p2 q2) (l : ℝ) :
    ‖ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2 + l * 1)
        (ssT p1 q1 p2 q2 + l * (ssB p1 q1 p2 q2 / ssE p2 q2))‖ ^ 2 =
      ‖ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2)‖ ^ 2 +
        l ^ 2 * (ssDen p1 q1 p2 q2 / ssE p2 q2) := by
  rw [ssVec_sq_split hD]
  congr 1
  have e : (ssS p1 q1 p2 q2 + l * 1 - ssS p1 q1 p2 q2) • (q1 - p1) -
      (ssT p1 q1 p2 q2 + l * (ssB p1 q1 p2 q2 / ssE p2 q2) - ssT p1 q1 p2 q2) • (q2 - p2) =
      l • ((q1 - p1) - (ssB p1 q1 p2 q2 / ssE p2 q2) • (q2 - p2)) := by module
  rw [e]
  exact valley_aux (q1 - p1) (q2 - p2) l he

/-- **The bound of survey finding 1**: for `ε ≥ 0` the former routine satisfies
`segSegThr² ≤ d² + (ε/4)·‖q₁ − p₁‖²`, `d = segSegDist` the distance. (Proof: when the interior
pair is skipped, `0 < denom ≤ ε a e`; walking from `(s*, t*)` along the valley to the nearer of
`s = 0`, `s = 1` reaches the boundary with `|s − s*| ≤ ½`, where the squared distance is at most
`d² + denom/(4e)`. Two segments crossing at their midpoints with `denom = ε a e` attain it.) -/
theorem segSegThr_sq_le {ε : ℝ} (hε : 0 ≤ ε) :
    segSegThr p1 q1 p2 q2 ε ^ 2 ≤ segSegDist p1 q1 p2 q2 ^ 2 + ε / 4 * ssA p1 q1 := by
  have hA := ssA_nonneg (p1 := p1) (q1 := q1)
  by_cases hc : 0 < ssA p1 q1 ∧ 0 < ssE p2 q2 ∧ 0 < ssDen p1 q1 p2 q2 ∧
      InSq (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2)
  · obtain ⟨ha, he, hD, hi⟩ := hc
    -- the distance is the interior pair
    have hdist : segSegDist p1 q1 p2 q2 =
        ‖ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2) (ssT p1 q1 p2 q2)‖ := by
      apply le_antisymm
      · unfold segSegDist; rw [ite_eq_left ⟨ha, he, hD, hi⟩]; exact min_le_right _ _
      · obtain ⟨s, t, h, e⟩ := segSegDist_attained (p1 := p1) (q1 := q1) (p2 := p2) (q2 := q2)
        rw [e]
        refine norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) ?_
        rw [ssVec_sq_split hD.ne' s t]
        nlinarith [sq_nonneg ‖(s - ssS p1 q1 p2 q2) • (q1 - p1) -
          (t - ssT p1 q1 p2 q2) • (q2 - p2)‖]
    by_cases hthr : ε * ssA p1 q1 * ssE p2 q2 < ssDen p1 q1 p2 q2
    · have : segSegThr p1 q1 p2 q2 ε = segSegDist p1 q1 p2 q2 := by
        unfold segSegThr segSegDist; rw [ite_eq_left ⟨ha, he, hthr, hi⟩, ite_eq_left ⟨ha, he, hD, hi⟩]
      rw [this]; nlinarith
    · have hthr' : ssDen p1 q1 p2 q2 ≤ ε * ssA p1 q1 * ssE p2 q2 := not_lt.1 hthr
      have hval : segSegThr p1 q1 p2 q2 ε = ssBest p1 q1 p2 q2 := by
        unfold segSegThr; rw [ite_eq_right (fun h => hthr h.2.2.1)]
      rw [hval, hdist]
      -- walk along the valley towards the nearer end of the first segment
      have key : ∃ l : ℝ, l ^ 2 ≤ 1 / 4 ∧ ssBest p1 q1 p2 q2 ≤
          ‖ssVec p1 q1 p2 q2 (ssS p1 q1 p2 q2 + l * 1)
            (ssT p1 q1 p2 q2 + l * (ssB p1 q1 p2 q2 / ssE p2 q2))‖ := by
        rcases le_total (ssS p1 q1 p2 q2) (1 / 2) with hs | hs
        · -- towards s = 0: parameter `-λ`, `λ ∈ [0, s*]`
          obtain ⟨l, hl, hin, hbd⟩ := exists_sq_bdry (s₀ := ssS p1 q1 p2 q2)
            (t₀ := ssT p1 q1 p2 q2) (u := -1) (v := -(ssB p1 q1 p2 q2 / ssE p2 q2))
            (L := ssS p1 q1 p2 q2) hi.1.1 hi (by
              by_cases h : InSq (ssS p1 q1 p2 q2 + ssS p1 q1 p2 q2 * -1)
                  (ssT p1 q1 p2 q2 + ssS p1 q1 p2 q2 * -(ssB p1 q1 p2 q2 / ssE p2 q2))
              · exact Or.inr ⟨h, Or.inl (by ring)⟩
              · exact Or.inl h)
          refine ⟨-l, by nlinarith [hl.1, hl.2], ?_⟩
          have e1 : ssS p1 q1 p2 q2 + -l * 1 = ssS p1 q1 p2 q2 + l * -1 := by ring
          have e2 : ssT p1 q1 p2 q2 + -l * (ssB p1 q1 p2 q2 / ssE p2 q2) =
              ssT p1 q1 p2 q2 + l * -(ssB p1 q1 p2 q2 / ssE p2 q2) := by ring
          rw [e1, e2]; exact ssBest_le_of_bdry hin hbd
        · obtain ⟨l, hl, hin, hbd⟩ := exists_sq_bdry (s₀ := ssS p1 q1 p2 q2)
            (t₀ := ssT p1 q1 p2 q2) (u := 1) (v := ssB p1 q1 p2 q2 / ssE p2 q2)
            (L := 1 - ssS p1 q1 p2 q2) (by linarith [hi.1.2]) hi (by
              by_cases h : InSq (ssS p1 q1 p2 q2 + (1 - ssS p1 q1 p2 q2) * 1)
                  (ssT p1 q1 p2 q2 + (1 - ssS p1 q1 p2 q2) * (ssB p1 q1 p2 q2 / ssE p2 q2))
              · exact Or.inr ⟨h, Or.inr (Or.inl (by ring))⟩
              · exact Or.inl h)
          exact ⟨l, by nlinarith [hl.1, hl.2], ssBest_le_of_bdry hin hbd⟩
      obtain ⟨l, hl, hle⟩ := key
      have h0 : (0 : ℝ) ≤ ssBest p1 q1 p2 q2 := by unfold ssBest psDist; positivity
      have hb := pow_le_pow_left₀ h0 hle 2
      rw [ssVec_sq_valley hD.ne' he] at hb
      have hq : ssDen p1 q1 p2 q2 / ssE p2 q2 ≤ ε * ssA p1 q1 := by
        rw [div_le_iff₀ he]; linarith
      have hq0 : 0 ≤ ssDen p1 q1 p2 q2 / ssE p2 q2 := div_nonneg ssDen_nonneg he.le
      have := mul_le_mul hl hq hq0 (by norm_num : (0 : ℝ) ≤ 1 / 4)
      linarith
  · have : segSegThr p1 q1 p2 q2 ε = segSegDist p1 q1 p2 q2 := by
      apply le_antisymm
      · unfold segSegThr segSegDist
        rw [ite_eq_right hc]
        split_ifs with h
        · exact (hc ⟨h.1, h.2.1, lt_of_le_of_lt (mul_nonneg (mul_nonneg hε ssA_nonneg)
            ssE_nonneg) h.2.2.1, h.2.2.2⟩).elim
        · exact le_rfl
      · exact segSegDist_le_thr hε
    rw [this]; nlinarith

end SegSeg

/-! ## Point–triangle -/

/-- The point `a + v (b − a) + w (c − a)` of the plane of `a b c`. -/
def triPt (a b c : E) (v w : ℝ) : E := a + v • (b - a) + w • (c - a)

/-- The triangle in the coordinates `(v, w)`. -/
abbrev InTri (v w : ℝ) : Prop := 0 ≤ v ∧ 0 ≤ w ∧ v + w ≤ 1

/-- A point of the triangle on one of its edges. -/
def OnTriBdry (v w : ℝ) : Prop := v = 0 ∨ w = 0 ∨ v + w = 1

/-- A segment of coordinates from a point of the triangle to a point outside meets an edge. -/
theorem exists_tri_bdry {v₀ w₀ u z L : ℝ} (hL : 0 ≤ L) (h₀ : InTri v₀ w₀)
    (hout : ¬ InTri (v₀ + L * u) (w₀ + L * z)) :
    ∃ l ∈ Icc 0 L, InTri (v₀ + l * u) (w₀ + l * z) ∧ OnTriBdry (v₀ + l * u) (w₀ + l * z) := by
  let g : ℝ → ℝ := fun l => min (min (v₀ + l * u) (w₀ + l * z)) (1 - (v₀ + l * u) - (w₀ + l * z))
  have hg : Continuous g := by fun_prop
  have hg0 : 0 ≤ g 0 := by
    obtain ⟨a1, a2, a3⟩ := h₀
    simp only [g, zero_mul, add_zero]
    exact le_min (le_min a1 a2) (by linarith)
  have hgL : g L < 0 := by
    simp only [g]
    by_contra hc
    push Not at hc
    apply hout
    have h1 := (le_min_iff.1 hc).1
    have h2 := (le_min_iff.1 hc).2
    exact ⟨(le_min_iff.1 h1).1, (le_min_iff.1 h1).2, by linarith⟩
  obtain ⟨l, hl, hgl⟩ := exists_zero_of_nonneg_of_neg hg hL hg0 hgL
  refine ⟨l, hl, ?_, ?_⟩
  · have hc : 0 ≤ g l := hgl.ge
    simp only [g] at hc
    have h1 := (le_min_iff.1 hc).1
    have h2 := (le_min_iff.1 hc).2
    exact ⟨(le_min_iff.1 h1).1, (le_min_iff.1 h1).2, by linarith⟩
  · simp only [g] at hgl
    unfold OnTriBdry
    rcases min_choice (min (v₀ + l * u) (w₀ + l * z)) (1 - (v₀ + l * u) - (w₀ + l * z)) with h | h <;>
      rw [h] at hgl
    · rcases min_choice (v₀ + l * u) (w₀ + l * z) with h' | h' <;> rw [h'] at hgl
      · exact Or.inl hgl
      · exact Or.inr (Or.inl hgl)
    · exact Or.inr (Or.inr (by linarith))

section PtTri

variable (p a b c : E)

/-- The scalars: `A = ab·ab`, `B = ab·ac`, `C = ac·ac`, `P = ap·ab`, `Q = ap·ac`,
`N = AC − B² = |ab × ac|²` (Lagrange). -/
def tA : ℝ := ⟪b - a, b - a⟫
def tB : ℝ := ⟪b - a, c - a⟫
def tC : ℝ := ⟪c - a, c - a⟫
def tP : ℝ := ⟪p - a, b - a⟫
def tQ : ℝ := ⟪p - a, c - a⟫
def tN : ℝ := tA a b * tC a c - tB a b c ^ 2

/-- The coordinates of the projection of `p` onto the plane. -/
def tV : ℝ := (tP p a b * tC a c - tQ p a c * tB a b c) / tN a b c
def tW : ℝ := (tQ p a c * tA a b - tP p a b * tB a b c) / tN a b c

/-- The three edge distances. -/
def ptEdges : ℝ := min (min (psDist p a b) (psDist p b c)) (psDist p c a)

/-- **`point_triangle_distance` of `ccd.rs`** (exact arithmetic): the edge distances, and the
distance to the projection onto the plane when the triangle is non-degenerate and the
projection lies in it (in `ℝ³` that is `|ap·n|/|n|`, `n = ab × ac`: `cross_formulas_ptTri`).

[formula: rc-closest-points] -/
def ptTriDist : ℝ :=
  if tN a b c ≤ 0 then ptEdges p a b c
  else if InTri (tV p a b c) (tW p a b c) then
    min (ptEdges p a b c) ‖p - triPt a b c (tV p a b c) (tW p a b c)‖
  else ptEdges p a b c

theorem triPt_sub (v w : ℝ) : p - triPt a b c v w = (p - a) - v • (b - a) - w • (c - a) := by
  unfold triPt; abel

variable {p a b c}

theorem ptEdges_le_of_bdry {v w : ℝ} (h : InTri v w) (hb : OnTriBdry v w) :
    ptEdges p a b c ≤ ‖p - triPt a b c v w‖ := by
  obtain ⟨h1, h2, h3⟩ := h
  unfold ptEdges
  rcases hb with rfl | rfl | hvw
  · have e : triPt a b c 0 w = segPt c a (1 - w) := by unfold triPt segPt; module
    rw [e]
    exact (min_le_right _ _).trans (psDist_le p c a ⟨by linarith, by linarith⟩)
  · have e : triPt a b c v 0 = segPt a b v := by unfold triPt segPt; module
    rw [e]
    exact (min_le_left _ _).trans ((min_le_left _ _).trans
      (psDist_le p a b ⟨h1, by linarith⟩))
  · have e : triPt a b c v w = segPt b c w := by
      unfold triPt segPt
      have : v = 1 - w := by linarith
      subst this; module
    rw [e]
    exact (min_le_left _ _).trans ((min_le_right _ _).trans
      (psDist_le p b c ⟨h2, by linarith⟩))

theorem tri_star_aux (x d1 d2 : E) (hD : ⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2 ≠ 0) :
    ⟪x - ((⟪x, d1⟫ * ⟪d2, d2⟫ - ⟪x, d2⟫ * ⟪d1, d2⟫) / (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2)) • d1 -
        ((⟪x, d2⟫ * ⟪d1, d1⟫ - ⟪x, d1⟫ * ⟪d1, d2⟫) / (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2)) • d2,
      d1⟫ = 0 ∧
    ⟪x - ((⟪x, d1⟫ * ⟪d2, d2⟫ - ⟪x, d2⟫ * ⟪d1, d2⟫) / (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2)) • d1 -
        ((⟪x, d2⟫ * ⟪d1, d1⟫ - ⟪x, d1⟫ * ⟪d1, d2⟫) / (⟪d1, d1⟫ * ⟪d2, d2⟫ - ⟪d1, d2⟫ ^ 2)) • d2,
      d2⟫ = 0 := by
  simp only [inner_sub_left, inner_smul_left, RCLike.conj_to_real]
  rw [real_inner_comm d1 d2]
  generalize ⟪d1, d1⟫ = A at hD ⊢
  generalize ⟪d2, d2⟫ = C at hD ⊢
  generalize ⟪d1, d2⟫ = B at hD ⊢
  generalize ⟪x, d1⟫ = P
  generalize ⟪x, d2⟫ = Q
  have hS : (P * C - Q * B) / (A * C - B ^ 2) * (A * C - B ^ 2) = P * C - Q * B :=
    div_mul_cancel₀ _ hD
  have hT : (Q * A - P * B) / (A * C - B ^ 2) * (A * C - B ^ 2) = Q * A - P * B :=
    div_mul_cancel₀ _ hD
  constructor <;> refine mul_left_cancel₀ hD ?_ <;> rw [mul_zero]
  · linear_combination -A * hS - B * hT
  · linear_combination -B * hS - C * hT

/-- The projection onto the plane is stationary: `p − π ⊥ ab, ac`. -/
theorem triPt_star_inner (hN : tN a b c ≠ 0) :
    ⟪p - triPt a b c (tV p a b c) (tW p a b c), b - a⟫ = 0 ∧
      ⟪p - triPt a b c (tV p a b c) (tW p a b c), c - a⟫ = 0 := by
  rw [triPt_sub]
  exact tri_star_aux (p - a) (b - a) (c - a) hN

/-- `‖p − x‖² = ‖p − π‖² + ‖(v − v*) ab + (w − w*) ac‖²` for every point `x` of the plane. -/
theorem triPt_sq_split (hN : tN a b c ≠ 0) (v w : ℝ) :
    ‖p - triPt a b c v w‖ ^ 2 = ‖p - triPt a b c (tV p a b c) (tW p a b c)‖ ^ 2 +
      ‖(v - tV p a b c) • (b - a) + (w - tW p a b c) • (c - a)‖ ^ 2 := by
  obtain ⟨h1, h2⟩ := triPt_star_inner (p := p) hN
  have e : p - triPt a b c v w = (p - triPt a b c (tV p a b c) (tW p a b c)) -
      ((v - tV p a b c) • (b - a) + (w - tW p a b c) • (c - a)) := by
    rw [triPt_sub, triPt_sub]; module
  rw [e, norm_sub_sq_real, inner_add_right, inner_smul_right, inner_smul_right, h1, h2]
  ring

theorem tA_nonneg : 0 ≤ tA a b := real_inner_self_nonneg

theorem tN_nonneg : 0 ≤ tN a b c := by
  unfold tN tA tB tC
  nlinarith [real_inner_mul_inner_self_le (b - a) (c - a)]

/-- **Point–triangle distance, lower bound**: `ptTriDist ≤ ‖p − x‖` for every point `x` of the
(closed) triangle, degenerate triangles included.

[formula: rc-closest-points] -/
theorem ptTriDist_le {v w : ℝ} (h : InTri v w) : ptTriDist p a b c ≤ ‖p - triPt a b c v w‖ := by
  by_cases hN : tN a b c ≤ 0
  · -- collinear corners: slide along an affine dependence to an edge
    unfold ptTriDist; rw [ite_eq_left hN]
    have hN0 : tN a b c = 0 := le_antisymm hN tN_nonneg
    by_cases hA : 0 < tA a b
    · have hpar : tA a b • (c - a) - tB a b c • (b - a) = 0 := by
        have := norm_sq_smul_sub_aux (b - a) (c - a)
        have h' : ‖tA a b • (c - a) - tB a b c • (b - a)‖ ^ 2 = tA a b * tN a b c := this
        rw [hN0, mul_zero] at h'
        exact norm_eq_zero.1 (pow_eq_zero_iff two_ne_zero |>.1 h')
      have hconst : ∀ l : ℝ, triPt a b c (v + l * tB a b c) (w + l * -tA a b) =
          triPt a b c v w := fun l => by
        unfold triPt
        have : (v + l * tB a b c) • (b - a) + (w + l * -tA a b) • (c - a) =
            v • (b - a) + w • (c - a) - l • (tA a b • (c - a) - tB a b c • (b - a)) := by module
        rw [add_assoc, this, hpar, smul_zero, sub_zero, ← add_assoc]
      obtain ⟨l, -, hin, hbd⟩ := exists_tri_bdry (v₀ := v) (w₀ := w) (u := tB a b c)
        (z := -tA a b) (L := 2 / tA a b) (by positivity) h (fun hi => by
          have := hi.2.1
          rw [mul_neg, div_mul_cancel₀ _ hA.ne'] at this
          linarith [h.2.2, h.1])
      rw [← hconst l]
      exact ptEdges_le_of_bdry hin hbd
    · have hA0 : tA a b = 0 := le_antisymm (not_lt.1 hA) tA_nonneg
      have hd : b - a = 0 := by
        have h' : ⟪b - a, b - a⟫ = 0 := hA0
        rwa [inner_self_eq_zero (𝕜 := ℝ)] at h'
      have hconst : ∀ l : ℝ, triPt a b c (v + l * 1) (w + l * 0) = triPt a b c v w := fun l => by
        unfold triPt; rw [hd]; simp
      obtain ⟨l, -, hin, hbd⟩ := exists_tri_bdry (v₀ := v) (w₀ := w) (u := 1) (z := 0) (L := 2)
        (by norm_num) h (fun hi => by have := hi.2.2; linarith [h.1, h.2.1])
      rw [← hconst l]
      exact ptEdges_le_of_bdry hin hbd
  · have hN' : tN a b c ≠ 0 := by intro h0; exact hN h0.le
    unfold ptTriDist; rw [ite_eq_right hN]
    split_ifs with hi
    · refine (min_le_right _ _).trans (norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) ?_)
      rw [triPt_sq_split hN' v w]
      nlinarith [sq_nonneg ‖(v - tV p a b c) • (b - a) + (w - tW p a b c) • (c - a)‖]
    · obtain ⟨l, hl, hin, hbd⟩ := exists_tri_bdry (v₀ := v) (w₀ := w) (u := tV p a b c - v)
        (z := tW p a b c - w) zero_le_one h (by simpa using hi)
      refine (ptEdges_le_of_bdry hin hbd).trans
        (norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) ?_)
      rw [triPt_sq_split hN', triPt_sq_split hN' v w]
      have e : (v + l * (tV p a b c - v) - tV p a b c) • (b - a) +
          (w + l * (tW p a b c - w) - tW p a b c) • (c - a) =
          (1 - l) • ((v - tV p a b c) • (b - a) + (w - tW p a b c) • (c - a)) := by module
      rw [e, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs]
      have : (1 - l) ^ 2 ≤ 1 := by nlinarith [hl.1, hl.2]
      nlinarith [sq_nonneg ‖(v - tV p a b c) • (b - a) + (w - tW p a b c) • (c - a)‖]

theorem ptEdges_attained : ∃ v w, InTri v w ∧ ptEdges p a b c = ‖p - triPt a b c v w‖ := by
  obtain ⟨m1, m2⟩ := psParam_mem p a b
  obtain ⟨n1, n2⟩ := psParam_mem p b c
  obtain ⟨o1, o2⟩ := psParam_mem p c a
  unfold ptEdges
  rcases min_choice (min (psDist p a b) (psDist p b c)) (psDist p c a) with h | h <;> rw [h]
  · rcases min_choice (psDist p a b) (psDist p b c) with h' | h' <;> rw [h', psDist_eq]
    · exact ⟨psParam p a b, 0, ⟨m1, le_rfl, by linarith⟩, by unfold triPt segPt; simp⟩
    · refine ⟨1 - psParam p b c, psParam p b c, ⟨by linarith, n1, by linarith⟩, ?_⟩
      unfold triPt segPt; congr 1; module
  · rw [psDist_eq]
    refine ⟨0, 1 - psParam p c a, ⟨le_rfl, by linarith, by linarith⟩, ?_⟩
    unfold triPt segPt; congr 1; module

/-- **Point–triangle distance, attained**: `ptTriDist` is the distance to a point of the
triangle; with `ptTriDist_le` it is the minimum distance.

[formula: rc-closest-points] -/
theorem ptTriDist_attained : ∃ v w, InTri v w ∧ ptTriDist p a b c = ‖p - triPt a b c v w‖ := by
  unfold ptTriDist
  split_ifs with hN hi
  · exact ptEdges_attained
  · rcases min_choice (ptEdges p a b c) ‖p - triPt a b c (tV p a b c) (tW p a b c)‖ with h | h <;>
      rw [h]
    · exact ptEdges_attained
    · exact ⟨_, _, hi, rfl⟩
  · exact ptEdges_attained

end PtTri

/-! ## The Voronoi-region routine (Ericson §5.1.5): `truth.rs`, `obstacle.rs` -/

section Voronoi

variable (p a b c : E)

/-- The dot products of the routine: `d₁ = ab·ap`, `d₂ = ac·ap`, `d₃ = ab·bp`, `d₄ = ac·bp`,
`d₅ = ab·cp`, `d₆ = ac·cp`, and `v_c = d₁d₄ − d₃d₂`, `v_b = d₅d₂ − d₁d₆`, `v_a = d₃d₆ − d₅d₄`. -/
def vd1 : ℝ := ⟪b - a, p - a⟫
def vd2 : ℝ := ⟪c - a, p - a⟫
def vd3 : ℝ := ⟪b - a, p - b⟫
def vd4 : ℝ := ⟪c - a, p - b⟫
def vd5 : ℝ := ⟪b - a, p - c⟫
def vd6 : ℝ := ⟪c - a, p - c⟫
def vVc : ℝ := vd1 p a b * vd4 p a b c - vd3 p a b * vd2 p a c
def vVb : ℝ := vd5 p a b c * vd2 p a c - vd1 p a b * vd6 p a c
def vVa : ℝ := vd3 p a b * vd6 p a c - vd5 p a b c * vd4 p a b c

/-- `closest_point_on_segment` of `truth.rs`. -/
def cpSeg (p a b : E) : E :=
  if ⟪b - a, b - a⟫ = 0 then a else a + clamp01 (⟪p - a, b - a⟫ / ⟪b - a, b - a⟫) • (b - a)

/-- `x`, replaced by `y` if `y` is strictly nearer to `p` (the update step of the loop in
`truth.rs`). -/
def nearer (p x y : E) : E := if ‖p - y‖ ^ 2 < ‖p - x‖ ^ 2 then y else x

/-- The degenerate fallback of `truth.rs`: the nearest of the three edge points (the first one
on ties). -/
def nearestEdgePt : E := nearer p (nearer p (cpSeg p a b) (cpSeg p b c)) (cpSeg p c a)

/-- **`closest_point_on_triangle` of `crates/ddg-points-as-tori/src/truth.rs`** (exact
arithmetic; `den.is_finite()` is always true for reals): the Voronoi regions of `a`, `b`, the
edge `ab`, `c`, the edge `ac`, the edge `bc`, then the projection onto the plane, then the
nearest edge point for a degenerate triangle.

[formula: pat-truth-closest-point] -/
def closestPtTri : E :=
  if vd1 p a b ≤ 0 ∧ vd2 p a c ≤ 0 then a
  else if 0 ≤ vd3 p a b ∧ vd4 p a b c ≤ vd3 p a b then b
  else if vVc p a b c ≤ 0 ∧ 0 ≤ vd1 p a b ∧ vd3 p a b ≤ 0 ∧ 0 < vd1 p a b - vd3 p a b then
    a + (vd1 p a b / (vd1 p a b - vd3 p a b)) • (b - a)
  else if 0 ≤ vd6 p a c ∧ vd5 p a b c ≤ vd6 p a c then c
  else if vVb p a b c ≤ 0 ∧ 0 ≤ vd2 p a c ∧ vd6 p a c ≤ 0 ∧ 0 < vd2 p a c - vd6 p a c then
    a + (vd2 p a c / (vd2 p a c - vd6 p a c)) • (c - a)
  else if vVa p a b c ≤ 0 ∧ 0 ≤ vd4 p a b c - vd3 p a b ∧ 0 ≤ vd5 p a b c - vd6 p a c ∧
      0 < (vd4 p a b c - vd3 p a b) + (vd5 p a b c - vd6 p a c) then
    b + ((vd4 p a b c - vd3 p a b) / ((vd4 p a b c - vd3 p a b) + (vd5 p a b c - vd6 p a c))) •
      (c - b)
  else if 0 < vVa p a b c + vVb p a b c + vVc p a b c then
    a + (vVb p a b c / (vVa p a b c + vVb p a b c + vVc p a b c)) • (b - a) +
      (vVc p a b c / (vVa p a b c + vVb p a b c + vVc p a b c)) • (c - a)
  else nearestEdgePt p a b c

/-- **`closest_point_on_triangle` of `crates/ddg-repulsive-curves/src/obstacle.rs`** (no
degenerate guards; the obstacle rejects degenerate triangles).

[formula: rc-closest-points] -/
def closestPtTriRC : E :=
  if vd1 p a b ≤ 0 ∧ vd2 p a c ≤ 0 then a
  else if 0 ≤ vd3 p a b ∧ vd4 p a b c ≤ vd3 p a b then b
  else if vVc p a b c ≤ 0 ∧ 0 ≤ vd1 p a b ∧ vd3 p a b ≤ 0 then
    a + (vd1 p a b / (vd1 p a b - vd3 p a b)) • (b - a)
  else if 0 ≤ vd6 p a c ∧ vd5 p a b c ≤ vd6 p a c then c
  else if vVb p a b c ≤ 0 ∧ 0 ≤ vd2 p a c ∧ vd6 p a c ≤ 0 then
    a + (vd2 p a c / (vd2 p a c - vd6 p a c)) • (c - a)
  else if vVa p a b c ≤ 0 ∧ 0 ≤ vd4 p a b c - vd3 p a b ∧ 0 ≤ vd5 p a b c - vd6 p a c then
    b + ((vd4 p a b c - vd3 p a b) / ((vd4 p a b c - vd3 p a b) + (vd5 p a b c - vd6 p a c))) •
      (c - b)
  else
    a + (vVb p a b c * (1 / (vVa p a b c + vVb p a b c + vVc p a b c))) • (b - a) +
      (vVc p a b c * (1 / (vVa p a b c + vVb p a b c + vVc p a b c))) • (c - a)

variable {p a b c}

/-- The routine's dot products in terms of `P = ap·ab`, `Q = ap·ac`, `A`, `B`, `C`. -/
theorem vd_eqs :
    vd1 p a b = tP p a b ∧ vd2 p a c = tQ p a c ∧ vd3 p a b = tP p a b - tA a b ∧
      vd4 p a b c = tQ p a c - tB a b c ∧ vd5 p a b c = tP p a b - tB a b c ∧
      vd6 p a c = tQ p a c - tC a c := by
  have e1 : p - b = (p - a) - (b - a) := by abel
  have e2 : p - c = (p - a) - (c - a) := by abel
  unfold vd1 vd2 vd3 vd4 vd5 vd6 tP tQ tA tB tC
  rw [e1, e2]
  generalize p - a = x
  generalize b - a = u
  generalize c - a = z
  simp only [inner_sub_right]
  rw [real_inner_comm x u, real_inner_comm x z, real_inner_comm u z]
  exact ⟨rfl, rfl, rfl, rfl, rfl, rfl⟩

/-- **Variational inequality, sufficient**: a point `q` of the triangle with
`⟪p − q, x − q⟫ ≤ 0` at the three corners is the closest point. -/
theorem vi_suff {q : E} {v₀ w₀ : ℝ} (hq : q = triPt a b c v₀ w₀)
    (ha : ⟪p - q, a - q⟫ ≤ 0) (hb : ⟪p - q, b - q⟫ ≤ 0) (hc : ⟪p - q, c - q⟫ ≤ 0)
    {v w : ℝ} (h : InTri v w) : ‖p - q‖ ≤ ‖p - triPt a b c v w‖ := by
  have e : triPt a b c v w - q = (1 - v - w) • (a - q) + v • (b - q) + w • (c - q) := by
    rw [hq]; unfold triPt; module
  have hk : ⟪p - q, triPt a b c v w - q⟫ ≤ 0 := by
    rw [e, inner_add_right, inner_add_right, inner_smul_right, inner_smul_right,
      inner_smul_right]
    obtain ⟨h1, h2, h3⟩ := h
    nlinarith [mul_nonneg (by linarith : (0 : ℝ) ≤ 1 - v - w) (neg_nonneg.2 ha),
      mul_nonneg h1 (neg_nonneg.2 hb), mul_nonneg h2 (neg_nonneg.2 hc)]
  refine norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) ?_
  have e2 : p - triPt a b c v w = (p - q) - (triPt a b c v w - q) := by abel
  rw [e2, norm_sub_sq_real (p - q) (triPt a b c v w - q)]
  nlinarith [sq_nonneg ‖triPt a b c v w - q‖]

/-- **Variational inequality, necessary**: at a closest point `q` of the triangle,
`⟪p − q, x − q⟫ ≤ 0` for every point `x` of the triangle. -/
theorem vi_nec {v₀ w₀ : ℝ} (h₀ : InTri v₀ w₀)
    (hmin : ∀ v w, InTri v w → ‖p - triPt a b c v₀ w₀‖ ≤ ‖p - triPt a b c v w‖)
    {v w : ℝ} (h : InTri v w) :
    ⟪p - triPt a b c v₀ w₀, triPt a b c v w - triPt a b c v₀ w₀⟫ ≤ 0 := by
  set q := triPt a b c v₀ w₀
  set x := triPt a b c v w
  set k := ⟪p - q, x - q⟫
  set m := ‖x - q‖ ^ 2
  have key : ∀ l : ℝ, 0 ≤ l → l ≤ 1 → 2 * l * k ≤ l ^ 2 * m := by
    intro l hl0 hl1
    have hin : InTri (v₀ + l * (v - v₀)) (w₀ + l * (w - w₀)) := by
      obtain ⟨a1, a2, a3⟩ := h₀; obtain ⟨b1, b2, b3⟩ := h
      refine ⟨by nlinarith, by nlinarith, by nlinarith⟩
    have hpt : triPt a b c (v₀ + l * (v - v₀)) (w₀ + l * (w - w₀)) = q + l • (x - q) := by
      simp only [q, x]; unfold triPt; module
    have := hmin _ _ hin
    rw [hpt] at this
    have h2 := pow_le_pow_left₀ (norm_nonneg _) this 2
    have e : p - (q + l • (x - q)) = (p - q) - l • (x - q) := by abel
    rw [e, norm_sub_smul_sq] at h2
    simp only [m, ← real_inner_self_eq_norm_sq] at *
    linarith
  by_contra hk
  push Not at hk
  have hm : 0 ≤ m := sq_nonneg _
  have hl0 : 0 ≤ k / (k + m) := div_nonneg hk.le (by linarith)
  have hl1 : k / (k + m) ≤ 1 := (div_le_one (by linarith)).2 (by linarith)
  have := key _ hl0 hl1
  have hkm : 0 < k + m := by linarith
  rw [div_pow] at this
  have h' : 2 * k * k * (k + m) ≤ k ^ 2 * m := by
    have := mul_le_mul_of_nonneg_right this (sq_nonneg (k + m))
    field_simp at this
    nlinarith [this]
  nlinarith [mul_pos hk hk, mul_pos (mul_pos hk hk) hkm]

/-- The coverage of the plane by the six Voronoi regions, in real arithmetic: `(v₀, w₀)` are the
coordinates of a point of the triangle satisfying the variational inequality at the corners
(`Ia, Ib, Ic ≤ 0`, written with `Y₁ = ⟪p − q, ab⟫`, `Y₂ = ⟪p − q, ac⟫`); if none of the six
region tests holds and `AC − B² > 0`, then `v_a, v_b, v_c ≥ 0`. -/
theorem voronoi_cover_real {A B C P Q v₀ w₀ : ℝ} (hA : 0 ≤ A) (hC : 0 ≤ C)
    (hN : 0 < A * C - B ^ 2) (h₀ : InTri v₀ w₀)
    (hIa : -(v₀ * (P - v₀ * A - w₀ * B) + w₀ * (Q - v₀ * B - w₀ * C)) ≤ 0)
    (hIb : (1 - v₀) * (P - v₀ * A - w₀ * B) - w₀ * (Q - v₀ * B - w₀ * C) ≤ 0)
    (hIc : -v₀ * (P - v₀ * A - w₀ * B) + (1 - w₀) * (Q - v₀ * B - w₀ * C) ≤ 0)
    (t1 : ¬(P ≤ 0 ∧ Q ≤ 0)) (t2 : ¬(0 ≤ P - A ∧ Q - B ≤ P - A))
    (t3 : ¬(A * Q - B * P ≤ 0 ∧ 0 ≤ P ∧ P - A ≤ 0 ∧ 0 < P - (P - A)))
    (t4 : ¬(0 ≤ Q - C ∧ P - B ≤ Q - C))
    (t5 : ¬(C * P - B * Q ≤ 0 ∧ 0 ≤ Q ∧ Q - C ≤ 0 ∧ 0 < Q - (Q - C)))
    (t6 : ¬(A * C - B ^ 2 - (C - B) * P - (A - B) * Q ≤ 0 ∧ 0 ≤ (Q - B) - (P - A) ∧
      0 ≤ (P - B) - (Q - C) ∧ 0 < ((Q - B) - (P - A)) + ((P - B) - (Q - C)))) :
    0 ≤ A * C - B ^ 2 - (C - B) * P - (A - B) * Q ∧ 0 ≤ C * P - B * Q ∧ 0 ≤ A * Q - B * P := by
  obtain ⟨hv, hw, hvw⟩ := h₀
  have hApos : 0 < A := by
    rcases hA.eq_or_lt with h | h
    · subst h; nlinarith [sq_nonneg B]
    · exact h
  have hCpos : 0 < C := by
    rcases hC.eq_or_lt with h | h
    · subst h; nlinarith [sq_nonneg B]
    · exact h
  have hK : 0 < A + C - 2 * B := by nlinarith [sq_nonneg (A - C), sq_nonneg (A + C - 2 * B)]
  set Y1 := P - v₀ * A - w₀ * B with hY1
  set Y2 := Q - v₀ * B - w₀ * C with hY2
  -- the weighted sum of the three inequalities is `0`
  have hsum : (1 - v₀ - w₀) * (-(v₀ * Y1 + w₀ * Y2)) + v₀ * ((1 - v₀) * Y1 - w₀ * Y2) +
      w₀ * (-v₀ * Y1 + (1 - w₀) * Y2) = 0 := by ring
  have hα : 0 ≤ 1 - v₀ - w₀ := by linarith
  -- each weighted term is `≤ 0`, so each is `0`
  have ta : (1 - v₀ - w₀) * (-(v₀ * Y1 + w₀ * Y2)) ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hα hIa
  have tb : v₀ * ((1 - v₀) * Y1 - w₀ * Y2) ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hv hIb
  have tc : w₀ * (-v₀ * Y1 + (1 - w₀) * Y2) ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hw hIc
  have za : (1 - v₀ - w₀) * (-(v₀ * Y1 + w₀ * Y2)) = 0 := by linarith
  have zb : v₀ * ((1 - v₀) * Y1 - w₀ * Y2) = 0 := by linarith
  have zc : w₀ * (-v₀ * Y1 + (1 - w₀) * Y2) = 0 := by linarith
  rcases hw.eq_or_lt with hw0 | hw0
  · -- on the edge `ab` (w₀ = 0)
    subst hw0
    simp only [zero_mul, add_zero, sub_zero] at *
    rcases hv.eq_or_lt with hv0 | hv0
    · subst hv0; simp only [zero_mul, sub_zero, one_mul] at *
      exact absurd ⟨by linarith, by linarith⟩ t1
    rcases hα.eq_or_lt with hα0 | hα0
    · have : v₀ = 1 := by linarith
      subst this; simp only [one_mul, sub_self, zero_mul, neg_mul] at *
      exact absurd ⟨by linarith, by linarith⟩ t2
    · have hY1z : Y1 = 0 := by
        have : (1 - v₀) * (-(v₀ * Y1)) = 0 := za
        rcases mul_eq_zero.1 this with h | h
        · linarith
        · rcases mul_eq_zero.1 (neg_eq_zero.1 h) with h' | h'
          · linarith
          · exact h'
      have hP : P = v₀ * A := by linarith
      refine absurd ⟨?_, by nlinarith, by nlinarith, by linarith⟩ t3
      have : Y2 ≤ 0 := by nlinarith
      nlinarith
  rcases hv.eq_or_lt with hv0 | hv0
  · -- on the edge `ac` (v₀ = 0, w₀ > 0)
    subst hv0
    simp only [zero_mul, add_zero, sub_zero, zero_add, one_mul, neg_mul] at *
    rcases hα.eq_or_lt with hα0 | hα0
    · have : w₀ = 1 := by linarith
      subst this; simp only [one_mul, sub_self, zero_mul, add_zero] at *
      exact absurd ⟨by linarith, by linarith⟩ t4
    · have hY2z : Y2 = 0 := by
        have h' : (1 - w₀) * (-(w₀ * Y2)) = 0 := by linarith
        rcases mul_eq_zero.1 h' with h | h
        · linarith
        · rcases mul_eq_zero.1 (neg_eq_zero.1 h) with h'' | h''
          · linarith
          · exact h''
      have hQ : Q = w₀ * C := by linarith
      refine absurd ⟨?_, by nlinarith, by nlinarith, by linarith⟩ t5
      have : Y1 ≤ 0 := by nlinarith
      nlinarith
  rcases hα.eq_or_lt with hα0 | hα0
  · -- on the edge `bc` (v₀ + w₀ = 1)
    have hv1 : v₀ = 1 - w₀ := by linarith
    have hb0 : (1 - v₀) * Y1 - w₀ * Y2 = 0 := by
      rcases mul_eq_zero.1 zb with h | h
      · linarith
      · exact h
    have hc0 : -v₀ * Y1 + (1 - w₀) * Y2 = 0 := by
      rcases mul_eq_zero.1 zc with h | h
      · linarith
      · exact h
    have hcon : Y2 - Y1 = 0 := by rw [hv1] at hb0 hc0; nlinarith
    refine absurd ⟨?_, ?_, ?_, ?_⟩ t6
    · have hva : A * C - B ^ 2 - (C - B) * P - (A - B) * Q =
          (A + C - 2 * B) * (-(v₀ * Y1 + w₀ * Y2)) +
            (A * w₀ - A - 2 * B * w₀ + B + C * w₀) * (Y2 - Y1) := by
        rw [hY1, hY2, hv1]; ring
      rw [hva, hcon, mul_zero, add_zero]
      exact mul_nonpos_of_nonneg_of_nonpos hK.le hIa
    · have : (Q - B) - (P - A) = w₀ * (A + C - 2 * B) + (Y2 - Y1) := by
        rw [hY1, hY2, hv1]; ring
      rw [this, hcon]; nlinarith
    · have : (P - B) - (Q - C) = v₀ * (A + C - 2 * B) - (Y2 - Y1) := by
        rw [hY1, hY2, hv1]; ring
      rw [this, hcon]; nlinarith
    · linarith
  · -- inside: `p − q ⊥ ab, ac`, so `(v_a, v_b, v_c) = (AC − B²)(1 − v₀ − w₀, v₀, w₀)`
    have ha0 : -(v₀ * Y1 + w₀ * Y2) = 0 := by
      rcases mul_eq_zero.1 za with h | h
      · linarith
      · exact h
    have hb0 : (1 - v₀) * Y1 - w₀ * Y2 = 0 := by
      rcases mul_eq_zero.1 zb with h | h
      · linarith
      · exact h
    have hc0 : -v₀ * Y1 + (1 - w₀) * Y2 = 0 := by
      rcases mul_eq_zero.1 zc with h | h
      · linarith
      · exact h
    have hY1z : Y1 = 0 := by linarith
    have hY2z : Y2 = 0 := by linarith
    have hP : P = v₀ * A + w₀ * B := by linarith
    have hQ : Q = v₀ * B + w₀ * C := by linarith
    rw [hP, hQ]
    refine ⟨?_, ?_, ?_⟩
    · have : A * C - B ^ 2 - (C - B) * (v₀ * A + w₀ * B) - (A - B) * (v₀ * B + w₀ * C) =
          (A * C - B ^ 2) * (1 - v₀ - w₀) := by ring
      rw [this]; positivity
    · have : C * (v₀ * A + w₀ * B) - B * (v₀ * B + w₀ * C) = (A * C - B ^ 2) * v₀ := by ring
      rw [this]; positivity
    · have : A * (v₀ * B + w₀ * C) - B * (v₀ * A + w₀ * B) = (A * C - B ^ 2) * w₀ := by ring
      rw [this]; positivity

/-- The variational-inequality terms of a point `q = triPt v₀ w₀` in scalars. -/
theorem vi_scalars (v₀ w₀ : ℝ) :
    ⟪p - triPt a b c v₀ w₀, a - triPt a b c v₀ w₀⟫ =
        -(v₀ * (tP p a b - v₀ * tA a b - w₀ * tB a b c) +
          w₀ * (tQ p a c - v₀ * tB a b c - w₀ * tC a c)) ∧
      ⟪p - triPt a b c v₀ w₀, b - triPt a b c v₀ w₀⟫ =
        (1 - v₀) * (tP p a b - v₀ * tA a b - w₀ * tB a b c) -
          w₀ * (tQ p a c - v₀ * tB a b c - w₀ * tC a c) ∧
      ⟪p - triPt a b c v₀ w₀, c - triPt a b c v₀ w₀⟫ =
        -v₀ * (tP p a b - v₀ * tA a b - w₀ * tB a b c) +
          (1 - w₀) * (tQ p a c - v₀ * tB a b c - w₀ * tC a c) := by
  have ea : a - triPt a b c v₀ w₀ = -(v₀ • (b - a) + w₀ • (c - a)) := by unfold triPt; abel
  have eb : b - triPt a b c v₀ w₀ = (1 - v₀) • (b - a) - w₀ • (c - a) := by unfold triPt; module
  have ec : c - triPt a b c v₀ w₀ = -v₀ • (b - a) + (1 - w₀) • (c - a) := by unfold triPt; module
  rw [ea, eb, ec, triPt_sub]
  unfold tP tQ tA tB tC
  generalize p - a = x
  generalize b - a = u
  generalize c - a = z
  simp only [inner_neg_right, inner_add_right, inner_sub_right, inner_smul_right,
    inner_sub_left, inner_smul_left, RCLike.conj_to_real]
  rw [real_inner_comm z u]
  refine ⟨by ring, by ring, by ring⟩

/-- A point of the triangle with the three scalar variational inequalities is the closest
point. -/
theorem closest_of_scalars {A B C P Q v₀ w₀ : ℝ} (eA : tA a b = A) (eB : tB a b c = B)
    (eC : tC a c = C) (eP : tP p a b = P) (eQ : tQ p a c = Q) (h₀ : InTri v₀ w₀)
    (hIa : -(v₀ * (P - v₀ * A - w₀ * B) + w₀ * (Q - v₀ * B - w₀ * C)) ≤ 0)
    (hIb : (1 - v₀) * (P - v₀ * A - w₀ * B) - w₀ * (Q - v₀ * B - w₀ * C) ≤ 0)
    (hIc : -v₀ * (P - v₀ * A - w₀ * B) + (1 - w₀) * (Q - v₀ * B - w₀ * C) ≤ 0) :
    (∃ v w, InTri v w ∧ triPt a b c v₀ w₀ = triPt a b c v w) ∧
      ∀ v w, InTri v w → ‖p - triPt a b c v₀ w₀‖ ≤ ‖p - triPt a b c v w‖ := by
  obtain ⟨ea, eb, ec⟩ := vi_scalars (p := p) (a := a) (b := b) (c := c) v₀ w₀
  rw [eA, eB, eC, eP, eQ] at ea eb ec
  exact ⟨⟨v₀, w₀, h₀, rfl⟩, fun v w h => vi_suff rfl (ea ▸ hIa) (eb ▸ hIb) (ec ▸ hIc) h⟩

omit [InnerProductSpace ℝ E] in
theorem nearer_le_left (p x y : E) : ‖p - nearer p x y‖ ≤ ‖p - x‖ := by
  unfold nearer; split_ifs with h
  · exact norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) h.le
  · exact le_rfl

omit [InnerProductSpace ℝ E] in
theorem nearer_le_right (p x y : E) : ‖p - nearer p x y‖ ≤ ‖p - y‖ := by
  unfold nearer; split_ifs with h
  · exact le_rfl
  · exact norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) (not_lt.1 h)

omit [InnerProductSpace ℝ E] in
theorem nearer_cases (p x y : E) : nearer p x y = x ∨ nearer p x y = y := by
  unfold nearer; split_ifs <;> simp

theorem cpSeg_eq (p a b : E) : cpSeg p a b = segPt a b (psParam p a b) := by
  unfold cpSeg psParam segPt
  by_cases h : ⟪b - a, b - a⟫ = 0
  · rw [ite_eq_left h, ite_eq_right (by rw [h]; exact lt_irrefl 0)]; simp
  · rw [ite_eq_right h, ite_eq_left (lt_of_le_of_ne real_inner_self_nonneg (Ne.symm h))]

theorem vden_eq : vVa p a b c + vVb p a b c + vVc p a b c = tN a b c := by
  obtain ⟨e1, e2, e3, e4, e5, e6⟩ := vd_eqs (p := p) (a := a) (b := b) (c := c)
  unfold vVa vVb vVc tN; rw [e1, e2, e3, e4, e5, e6]; ring

theorem proj_Y1 (A B C P Q : ℝ) (hN : 0 < A * C - B ^ 2) :
    P - (C * P - B * Q) / (A * C - B ^ 2) * A - (A * Q - B * P) / (A * C - B ^ 2) * B = 0 := by
  generalize hDd : A * C - B ^ 2 = D at *
  have := hN.ne'
  field_simp
  rw [← hDd]; ring

theorem proj_Y2 (A B C P Q : ℝ) (hN : 0 < A * C - B ^ 2) :
    Q - (C * P - B * Q) / (A * C - B ^ 2) * B - (A * Q - B * P) / (A * C - B ^ 2) * C = 0 := by
  generalize hDd : A * C - B ^ 2 = D at *
  have := hN.ne'
  field_simp
  rw [← hDd]; ring

/-- **The Voronoi-region routine of `truth.rs` returns the closest point of the triangle**
(exact arithmetic, any triangle, degenerate ones included): the result is a point of the
triangle, and no point of the triangle is closer to `p`.

[formula: pat-truth-closest-point] -/
theorem closestPtTri_isClosest :
    (∃ v w, InTri v w ∧ closestPtTri p a b c = triPt a b c v w) ∧
      ∀ v w, InTri v w → ‖p - closestPtTri p a b c‖ ≤ ‖p - triPt a b c v w‖ := by
  obtain ⟨e1, e2, e3, e4, e5, e6⟩ := vd_eqs (p := p) (a := a) (b := b) (c := c)
  have hA0 := tA_nonneg (a := a) (b := b)
  have hC0 : 0 ≤ tC a c := real_inner_self_nonneg
  unfold closestPtTri
  unfold vVc vVb vVa
  rw [e1, e2, e3, e4, e5, e6]
  set A := tA a b
  set B := tB a b c
  set C := tC a c
  set P := tP p a b
  set Q := tQ p a c
  have ta : triPt a b c 0 0 = a := by unfold triPt; simp
  have tb : triPt a b c 1 0 = b := by unfold triPt; simp
  have tc : triPt a b c 0 1 = c := by unfold triPt; simp
  split_ifs with t1 t2 t3 t4 t5 t6 t7
  · -- the corner `a`
    have h := closest_of_scalars (p := p) rfl rfl rfl rfl rfl (v₀ := 0) (w₀ := 0)
      ⟨le_rfl, le_rfl, by norm_num⟩ (by simp) (by simpa using t1.1) (by simpa using t1.2)
    rwa [ta] at h
  · -- the corner `b`
    have h := closest_of_scalars (p := p) rfl rfl rfl rfl rfl (v₀ := 1) (w₀ := 0)
      ⟨zero_le_one, le_rfl, by norm_num⟩ (by simp; linarith [t2.1]) (by simp)
      (by simp; linarith [t2.2])
    rwa [tb] at h
  · -- the edge `ab`
    obtain ⟨h1, h2, h3, h4⟩ := t3
    have hA : 0 < A := by linarith
    have ht : P - (P - A) = A := by ring
    rw [ht]
    have e : triPt a b c (P / A) 0 = a + (P / A) • (b - a) := by unfold triPt; simp
    rw [← e]
    have hY1 : P - P / A * A - 0 * B = 0 := by field_simp; ring
    refine closest_of_scalars (p := p) rfl rfl rfl rfl rfl ⟨div_nonneg h2 hA.le, le_rfl, by
      rw [add_zero, div_le_one hA]; linarith⟩ (by rw [hY1]; simp) (by rw [hY1]; simp) ?_
    rw [hY1]
    have : Q - P / A * B - 0 * C = (A * Q - B * P) / A := by field_simp; ring
    rw [this]; simp only [mul_zero, zero_add, sub_zero, one_mul]
    exact div_nonpos_of_nonpos_of_nonneg (by linarith) hA.le
  · -- the corner `c`
    have h := closest_of_scalars (p := p) rfl rfl rfl rfl rfl (v₀ := 0) (w₀ := 1)
      ⟨le_rfl, zero_le_one, by norm_num⟩ (by simp; linarith [t4.1])
      (by simp; linarith [t4.2]) (by simp)
    rwa [tc] at h
  · -- the edge `ac`
    obtain ⟨h1, h2, h3, h4⟩ := t5
    have hC : 0 < C := by linarith
    have ht : Q - (Q - C) = C := by ring
    rw [ht]
    have e : triPt a b c 0 (Q / C) = a + (Q / C) • (c - a) := by unfold triPt; simp
    rw [← e]
    have hY2 : Q - 0 * B - Q / C * C = 0 := by field_simp; ring
    refine closest_of_scalars (p := p) rfl rfl rfl rfl rfl ⟨le_rfl, div_nonneg h2 hC.le, by
      rw [zero_add, div_le_one hC]; linarith⟩ (by rw [hY2]; simp) ?_ (by rw [hY2]; simp)
    rw [hY2]
    have : P - 0 * A - Q / C * B = (C * P - B * Q) / C := by field_simp; ring
    rw [this]; simp only [mul_zero, sub_zero, one_mul]
    exact div_nonpos_of_nonpos_of_nonneg (by linarith) hC.le
  · -- the edge `bc`
    obtain ⟨h1, h2, h3, h4⟩ := t6
    set K := (Q - B - (P - A)) + (P - B - (Q - C)) with hK
    set t := (Q - B - (P - A)) / K
    have ht0 : 0 ≤ t := div_nonneg h2 h4.le
    have ht1 : t ≤ 1 := (div_le_one h4).2 (by linarith)
    have e : triPt a b c (1 - t) t = b + t • (c - b) := by unfold triPt; module
    rw [← e]
    have hKe : K = A + C - 2 * B := by rw [hK]; ring
    have htK : t * K = Q - B - (P - A) := div_mul_cancel₀ _ h4.ne'
    have hcon : (Q - (1 - t) * B - t * C) - (P - (1 - t) * A - t * B) = 0 := by
      rw [hKe] at htK; linear_combination -htK
    have hva : A * C - B ^ 2 - (C - B) * P - (A - B) * Q =
        (A + C - 2 * B) * (-((1 - t) * (P - (1 - t) * A - t * B) +
          t * (Q - (1 - t) * B - t * C))) +
        (A * t - A - 2 * B * t + B + C * t) *
          ((Q - (1 - t) * B - t * C) - (P - (1 - t) * A - t * B)) := by ring
    have h1' : (P - A) * (Q - C) - (P - B) * (Q - B) = A * C - B ^ 2 - (C - B) * P - (A - B) * Q :=
      by ring
    rw [h1', hva, hcon, mul_zero, add_zero] at h1
    have hK' : 0 < A + C - 2 * B := hKe ▸ h4
    refine closest_of_scalars (p := p) rfl rfl rfl rfl rfl ⟨by linarith, ht0, by linarith⟩
      ?_ ?_ ?_
    · by_contra hc
      push Not at hc
      nlinarith [mul_pos hK' hc]
    · have : (Q - (1 - t) * B - t * C) = (P - (1 - t) * A - t * B) := by linarith
      rw [this, show 1 - (1 - t) = t by ring, sub_self]
    · have : (Q - (1 - t) * B - t * C) = (P - (1 - t) * A - t * B) := by linarith
      rw [this]; exact le_of_eq (by ring)
  · -- inside: the projection onto the plane
    have hN : 0 < A * C - B ^ 2 := by
      have := vden_eq (p := p) (a := a) (b := b) (c := c)
      unfold vVa vVb vVc at this; rw [e1, e2, e3, e4, e5, e6] at this
      have h' : A * C - B ^ 2 = tN a b c := rfl
      rw [h', ← this]
      linarith [t7]
    have hden : (P - A) * (Q - C) - (P - B) * (Q - B) + ((P - B) * Q - P * (Q - C)) +
        (P * (Q - B) - (P - A) * Q) = A * C - B ^ 2 := by ring
    rw [hden]
    -- a closest point exists (`ptTriDist`), and satisfies the variational inequality
    obtain ⟨v₀, w₀, h₀, hd⟩ := ptTriDist_attained (p := p) (a := a) (b := b) (c := c)
    have hmin : ∀ v w, InTri v w → ‖p - triPt a b c v₀ w₀‖ ≤ ‖p - triPt a b c v w‖ :=
      fun v w h => hd ▸ ptTriDist_le h
    have via := vi_nec h₀ hmin (v := 0) (w := 0) ⟨le_rfl, le_rfl, by norm_num⟩
    have vib := vi_nec h₀ hmin (v := 1) (w := 0) ⟨zero_le_one, le_rfl, by norm_num⟩
    have vic := vi_nec h₀ hmin (v := 0) (w := 1) ⟨le_rfl, zero_le_one, by norm_num⟩
    rw [ta] at via; rw [tb] at vib; rw [tc] at vic
    obtain ⟨ea, eb, ec⟩ := vi_scalars (p := p) (a := a) (b := b) (c := c) v₀ w₀
    rw [ea] at via; rw [eb] at vib; rw [ec] at vic
    obtain ⟨ca, cb, cc⟩ := voronoi_cover_real hA0 hC0 hN h₀ via vib vic t1 t2
      (by intro h; exact t3 ⟨by linarith [h.1], h.2.1, h.2.2.1, h.2.2.2⟩) t4
      (by intro h; exact t5 ⟨by linarith [h.1], h.2.1, h.2.2.1, h.2.2.2⟩)
      (by intro h; exact t6 ⟨by linarith [h.1], h.2.1, h.2.2.1, h.2.2.2⟩)
    have e : triPt a b c ((C * P - B * Q) / (A * C - B ^ 2)) ((A * Q - B * P) / (A * C - B ^ 2)) =
        a + (((P - B) * Q - P * (Q - C)) / (A * C - B ^ 2)) • (b - a) +
        ((P * (Q - B) - (P - A) * Q) / (A * C - B ^ 2)) • (c - a) := by
      unfold triPt; congr 3 <;> ring
    rw [← e]
    have hY1 := proj_Y1 A B C P Q hN
    have hY2 := proj_Y2 A B C P Q hN
    refine closest_of_scalars (p := p) rfl rfl rfl rfl rfl
      ⟨div_nonneg cb hN.le, div_nonneg cc hN.le, ?_⟩
      (by rw [hY1, hY2]; simp) (by rw [hY1, hY2]; simp) (by rw [hY1, hY2]; simp)
    rw [← add_div, div_le_one hN]; linarith
  · -- degenerate: the nearest edge point
    have hN : tN a b c ≤ 0 := by
      have := vden_eq (p := p) (a := a) (b := b) (c := c)
      unfold vVa vVb vVc at this; rw [e1, e2, e3, e4, e5, e6] at this
      rw [← this]; linarith [not_lt.1 t7]
    have hdist : ptTriDist p a b c = ptEdges p a b c := by unfold ptTriDist; rw [ite_eq_left hN]
    -- the three candidates are points of the triangle at the edge distances
    have q1 : ∃ v w, InTri v w ∧ cpSeg p a b = triPt a b c v w ∧
        ‖p - cpSeg p a b‖ = psDist p a b := by
      obtain ⟨m1, m2⟩ := psParam_mem p a b
      exact ⟨psParam p a b, 0, ⟨m1, le_rfl, by linarith⟩, by
        rw [cpSeg_eq]; unfold triPt segPt; simp, by rw [cpSeg_eq, psDist_eq]⟩
    have q2 : ∃ v w, InTri v w ∧ cpSeg p b c = triPt a b c v w ∧
        ‖p - cpSeg p b c‖ = psDist p b c := by
      obtain ⟨m1, m2⟩ := psParam_mem p b c
      exact ⟨1 - psParam p b c, psParam p b c, ⟨by linarith, m1, by linarith⟩, by
        rw [cpSeg_eq]; unfold triPt segPt; module, by rw [cpSeg_eq, psDist_eq]⟩
    have q3 : ∃ v w, InTri v w ∧ cpSeg p c a = triPt a b c v w ∧
        ‖p - cpSeg p c a‖ = psDist p c a := by
      obtain ⟨m1, m2⟩ := psParam_mem p c a
      exact ⟨0, 1 - psParam p c a, ⟨le_rfl, by linarith, by linarith⟩, by
        rw [cpSeg_eq]; unfold triPt segPt; module, by rw [cpSeg_eq, psDist_eq]⟩
    have hle : ‖p - nearestEdgePt p a b c‖ ≤ ptEdges p a b c := by
      obtain ⟨-, -, -, -, d1⟩ := q1
      obtain ⟨-, -, -, -, d2⟩ := q2
      obtain ⟨-, -, -, -, d3⟩ := q3
      unfold nearestEdgePt ptEdges
      rw [← d1, ← d2, ← d3]
      refine le_min (le_min ?_ ?_) (nearer_le_right _ _ _)
      · exact (nearer_le_left _ _ _).trans (nearer_le_left _ _ _)
      · exact (nearer_le_left _ _ _).trans (nearer_le_right _ _ _)
    refine ⟨?_, fun v w h => hle.trans (hdist ▸ ptTriDist_le h)⟩
    obtain ⟨v1, w1, i1, r1, -⟩ := q1
    obtain ⟨v2, w2, i2, r2, -⟩ := q2
    obtain ⟨v3, w3, i3, r3, -⟩ := q3
    unfold nearestEdgePt
    rcases nearer_cases p (nearer p (cpSeg p a b) (cpSeg p b c)) (cpSeg p c a) with h | h <;>
      rw [h]
    · rcases nearer_cases p (cpSeg p a b) (cpSeg p b c) with h' | h' <;> rw [h']
      · exact ⟨v1, w1, i1, r1⟩
      · exact ⟨v2, w2, i2, r2⟩
    · exact ⟨v3, w3, i3, r3⟩

/-- **The copy in `obstacle.rs` agrees with it for non-degenerate triangles** (`AC − B² > 0`;
`Obstacle` rejects degenerate triangles), hence also returns the closest point.

[formula: rc-closest-points] -/
theorem closestPtTriRC_eq (hN : 0 < tN a b c) : closestPtTriRC p a b c = closestPtTri p a b c := by
  obtain ⟨e1, e2, e3, e4, e5, e6⟩ := vd_eqs (p := p) (a := a) (b := b) (c := c)
  have hden := vden_eq (p := p) (a := a) (b := b) (c := c)
  have hA0 := tA_nonneg (a := a) (b := b)
  have hC0 : 0 ≤ tC a c := real_inner_self_nonneg
  have hN' : 0 < tA a b * tC a c - tB a b c ^ 2 := hN
  have hA : 0 < tA a b := by
    rcases hA0.eq_or_lt with h | h
    · rw [← h] at hN'; nlinarith [sq_nonneg (tB a b c)]
    · exact h
  have hC : 0 < tC a c := by
    rcases hC0.eq_or_lt with h | h
    · rw [← h] at hN'; nlinarith [sq_nonneg (tB a b c)]
    · exact h
  have g3 : 0 < vd1 p a b - vd3 p a b := by rw [e1, e3]; linarith
  have g5 : 0 < vd2 p a c - vd6 p a c := by rw [e2, e6]; linarith
  have g6 : 0 < (vd4 p a b c - vd3 p a b) + (vd5 p a b c - vd6 p a c) := by
    rw [e3, e4, e5, e6]
    nlinarith [sq_nonneg (tA a b - tC a c), sq_nonneg (tA a b + tC a c - 2 * tB a b c)]
  have g7 : 0 < vVa p a b c + vVb p a b c + vVc p a b c := hden ▸ hN
  unfold closestPtTri closestPtTriRC
  simp only [g3, g5, g6, g7, and_true, ite_true, mul_one_div]

end Voronoi

/-! ## The distances are `primDist` of `Advancement.lean`; safety of the advancement -/

section Prim

theorem weights2_point (P : Fin 2 → E)
    (w : Weights (Fin 2)) : ∑ i, w.1 i • P i = segPt (P 0) (P 1) (w.1 1) := by
  have h := w.2.2
  rw [Fin.sum_univ_two] at h
  rw [Fin.sum_univ_two]; unfold segPt
  have : w.1 0 = 1 - w.1 1 := by linarith
  rw [this]; module

theorem weights2_mem (w : Weights (Fin 2)) : w.1 1 ∈ Icc (0 : ℝ) 1 := by
  have h := w.2.2
  rw [Fin.sum_univ_two] at h
  exact ⟨w.2.1 1, by linarith [w.2.1 0]⟩

/-- The weights `(1 − s, s)`. -/
def weights2 {s : ℝ} (hs : s ∈ Icc (0 : ℝ) 1) : Weights (Fin 2) :=
  ⟨![1 - s, s], fun i => by fin_cases i <;> simp [hs.1, hs.2], by simp [Fin.sum_univ_two]⟩

theorem weights3_point (P : Fin 3 → E)
    (w : Weights (Fin 3)) : ∑ i, w.1 i • P i = triPt (P 0) (P 1) (P 2) (w.1 1) (w.1 2) := by
  have h := w.2.2
  rw [Fin.sum_univ_three] at h
  rw [Fin.sum_univ_three]; unfold triPt
  have : w.1 0 = 1 - w.1 1 - w.1 2 := by linarith
  rw [this]; module

theorem weights3_mem (w : Weights (Fin 3)) : InTri (w.1 1) (w.1 2) := by
  have h := w.2.2
  rw [Fin.sum_univ_three] at h
  exact ⟨w.2.1 1, w.2.1 2, by linarith [w.2.1 0]⟩

/-- The weights `(1 − v − w, v, w)`. -/
def weights3 {v w : ℝ} (h : InTri v w) : Weights (Fin 3) :=
  ⟨![1 - v - w, v, w], fun i => by
    obtain ⟨h1, h2, h3⟩ := h
    fin_cases i <;> simp <;> linarith, by simp [Fin.sum_univ_three]; ring⟩

theorem weights1_point (P : Fin 1 → E)
    (w : Weights (Fin 1)) : ∑ i, w.1 i • P i = P 0 := by
  have h := w.2.2
  rw [Fin.sum_univ_one] at h ⊢
  rw [h, one_smul]

/-- The corner positions at time `τ`: `p_i − τ a_i`. -/
def posAt {ι F : Type*} [AddCommGroup F] [Module ℝ F] (p a : ι → F) (τ : ℝ) : ι → F :=
  fun i => p i - τ • a i

theorem primPoint_eq_sum {ι F : Type*} [Fintype ι] [NormedAddCommGroup F] [NormedSpace ℝ F]
    (p a : ι → F) (w : Weights ι) (τ : ℝ) : primPoint p a w τ = ∑ i, w.1 i • posAt p a τ i := rfl

theorem primDist_bdd {ι κ F : Type*} [Fintype ι] [Fintype κ] [NormedAddCommGroup F]
    [NormedSpace ℝ F] (p a : ι → F) (q b : κ → F) (τ : ℝ) :
    BddBelow (Set.range fun x : Weights ι × Weights κ =>
      ‖primPoint p a x.1 τ - primPoint q b x.2 τ‖) :=
  ⟨0, by rintro _ ⟨x, rfl⟩; exact norm_nonneg _⟩

/-- **The computed segment–segment distance is `primDist`** (exact arithmetic): for two
linearly moving segments, `segment_segment_distance` of the corner positions at time `τ` is the
infimum over the barycentric weights.

[formula: rc-closest-points] -/
theorem primDist_segSeg (p a q b : Fin 2 → E) (τ : ℝ) :
    primDist p a q b τ = segSegCode (posAt p a τ 0) (posAt p a τ 1) (posAt q b τ 0)
      (posAt q b τ 1) := by
  rw [segSegCode_eq]
  apply le_antisymm
  · obtain ⟨s, t, h, e⟩ := segSegDist_attained (p1 := posAt p a τ 0) (q1 := posAt p a τ 1)
      (p2 := posAt q b τ 0) (q2 := posAt q b τ 1)
    rw [e]
    refine (ciInf_le (primDist_bdd p a q b τ) (weights2 h.1, weights2 h.2)).trans (le_of_eq ?_)
    simp only [primPoint_eq_sum, weights2_point, ssVec]
    rfl
  · refine le_ciInf fun x => ?_
    simp only [primPoint_eq_sum, weights2_point]
    exact segSegDist_le ⟨weights2_mem x.1, weights2_mem x.2⟩

/-- **The computed point–triangle distance is `primDist`** (a moving point against a moving
triangle).

[formula: rc-closest-points] -/
theorem primDist_ptTri (p a : Fin 1 → E) (q b : Fin 3 → E) (τ : ℝ) :
    primDist p a q b τ = ptTriDist (posAt p a τ 0) (posAt q b τ 0) (posAt q b τ 1)
      (posAt q b τ 2) := by
  apply le_antisymm
  · obtain ⟨v, w, h, e⟩ := ptTriDist_attained (p := posAt p a τ 0) (a := posAt q b τ 0)
      (b := posAt q b τ 1) (c := posAt q b τ 2)
    rw [e]
    have w1 : Weights (Fin 1) := ⟨fun _ => 1, fun _ => zero_le_one, by simp⟩
    refine (ciInf_le (primDist_bdd p a q b τ) (w1, weights3 h)).trans (le_of_eq ?_)
    simp only [primPoint_eq_sum, weights1_point, weights3_point]
    rfl
  · refine le_ciInf fun x => ?_
    simp only [primPoint_eq_sum, weights1_point, weights3_point]
    exact ptTriDist_le (weights3_mem x.2)

theorem primPoint_sub_common {ν : Type*} [Fintype ν] (r c : ν → E) (m : E) (w : Weights ν)
    (τ : ℝ) : primPoint r (fun i => c i - m) w τ = primPoint r c w τ + τ • m := by
  unfold primPoint
  have : ∑ i, w.1 i • (r i - τ • (c i - m)) = ∑ i, w.1 i • (r i - τ • c i) +
      (∑ i, w.1 i) • (τ • m) := by
    rw [Finset.sum_smul, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    module
  rw [this, w.2.2, one_smul]

/-- `primDist` does not change when the same velocity `m` is subtracted from every corner
(`ccd.rs` measures the corner speeds relative to the mean velocity of the pair). -/
theorem primDist_sub_common {ι κ : Type*} [Fintype ι] [Fintype κ] (p a : ι → E) (q b : κ → E)
    (m : E) (τ : ℝ) :
    primDist p (fun i => a i - m) q (fun j => b j - m) τ = primDist p a q b τ := by
  unfold primDist
  congr 1; ext x
  rw [primPoint_sub_common, primPoint_sub_common]; congr 1; abel

/-- **CCD safety for two moving segments, end to end** (exact arithmetic): corner velocities
`−a_i`, `−b_j` (positions `p − τ a`), measured relative to a common velocity `m` with speeds
`≤ l_A`, `≤ l_B`; if the computed distance at `τ_k` exceeds the target `g > 0`, then on
`[τ_k, τ_k + (d_k − g)/(l_A + l_B)]` the computed distance stays `≥ g` and the segments are
disjoint. This is the advance of `ccd.rs` (and of the repulsive-surfaces CCD for edge pairs).

[formula: rc-closest-points] -/
theorem segSeg_advance_safe (p a q b : Fin 2 → E) (m : E) {lA lB : ℝ}
    (ha : ∀ i, ‖a i - m‖ ≤ lA) (hb : ∀ j, ‖b j - m‖ ≤ lB) (hl : 0 < lA + lB) {τk g : ℝ}
    (hg0 : 0 < g)
    (hg : g < segSegCode (posAt p a τk 0) (posAt p a τk 1) (posAt q b τk 0) (posAt q b τk 1)) :
    ∀ τ, τk ≤ τ → τ ≤ τk + (segSegCode (posAt p a τk 0) (posAt p a τk 1) (posAt q b τk 0)
        (posAt q b τk 1) - g) / (lA + lB) →
      g ≤ segSegCode (posAt p a τ 0) (posAt p a τ 1) (posAt q b τ 0) (posAt q b τ 1) ∧
      ∀ s t, InSq s t → segPt (posAt p a τ 0) (posAt p a τ 1) s ≠
        segPt (posAt q b τ 0) (posAt q b τ 1) t := by
  intro τ h1 h2
  have hf : ∀ σ, primDist p (fun i => a i - m) q (fun j => b j - m) σ =
      segSegCode (posAt p a σ 0) (posAt p a σ 1) (posAt q b σ 0) (posAt q b σ 1) := fun σ => by
    rw [primDist_sub_common, primDist_segSeg]
  have hlip : ∀ σ, τk ≤ σ → primDist p (fun i => a i - m) q (fun j => b j - m) τk -
      (lA + lB) * (σ - τk) ≤ primDist p (fun i => a i - m) q (fun j => b j - m) σ := by
    intro σ hσ
    have := primDist_lipschitz p (fun i => a i - m) q (fun j => b j - m) ha hb τk σ
    rwa [abs_of_nonneg (by linarith), mul_comm] at this
  have hsafe := (advance_safe hl hlip (by rw [hf]; exact hg)).2 τ h1 (by rw [hf]; exact h2)
  rw [hf] at hsafe
  refine ⟨hsafe.1, fun s t hst he => ?_⟩
  have := (segSegCode_eq (p1 := posAt p a τ 0) (q1 := posAt p a τ 1) (p2 := posAt q b τ 0)
    (q2 := posAt q b τ 1)) ▸ segSegDist_le hst
  unfold ssVec at this
  rw [he, sub_self, norm_zero] at this
  linarith [hsafe.1]

/-- **CCD safety for a moving point and a moving triangle, end to end** (exact arithmetic), as
`segSeg_advance_safe` (for `ccd.rs` the triangle is a rigid obstacle, `b = m`).

[formula: rc-closest-points] -/
theorem ptTri_advance_safe (p a : Fin 1 → E) (q b : Fin 3 → E) (m : E) {lA lB : ℝ}
    (ha : ∀ i, ‖a i - m‖ ≤ lA) (hb : ∀ j, ‖b j - m‖ ≤ lB) (hl : 0 < lA + lB) {τk g : ℝ}
    (hg0 : 0 < g)
    (hg : g < ptTriDist (posAt p a τk 0) (posAt q b τk 0) (posAt q b τk 1) (posAt q b τk 2)) :
    ∀ τ, τk ≤ τ → τ ≤ τk + (ptTriDist (posAt p a τk 0) (posAt q b τk 0) (posAt q b τk 1)
        (posAt q b τk 2) - g) / (lA + lB) →
      g ≤ ptTriDist (posAt p a τ 0) (posAt q b τ 0) (posAt q b τ 1) (posAt q b τ 2) ∧
      ∀ v w, InTri v w → posAt p a τ 0 ≠
        triPt (posAt q b τ 0) (posAt q b τ 1) (posAt q b τ 2) v w := by
  intro τ h1 h2
  have hf : ∀ σ, primDist p (fun i => a i - m) q (fun j => b j - m) σ =
      ptTriDist (posAt p a σ 0) (posAt q b σ 0) (posAt q b σ 1) (posAt q b σ 2) := fun σ => by
    rw [primDist_sub_common, primDist_ptTri]
  have hlip : ∀ σ, τk ≤ σ → primDist p (fun i => a i - m) q (fun j => b j - m) τk -
      (lA + lB) * (σ - τk) ≤ primDist p (fun i => a i - m) q (fun j => b j - m) σ := by
    intro σ hσ
    have := primDist_lipschitz p (fun i => a i - m) q (fun j => b j - m) ha hb τk σ
    rwa [abs_of_nonneg (by linarith), mul_comm] at this
  have hsafe := (advance_safe hl hlip (by rw [hf]; exact hg)).2 τ h1 (by rw [hf]; exact h2)
  rw [hf] at hsafe
  refine ⟨hsafe.1, fun v w hvw he => ?_⟩
  have := ptTriDist_le (p := posAt p a τ 0) (a := posAt q b τ 0) (b := posAt q b τ 1)
    (c := posAt q b τ 2) hvw
  rw [← he, sub_self, norm_zero] at this
  linarith [hsafe.1]

end Prim

/-! ## The cross-product forms of the Rust code (`ℝ³`, Lagrange's identity) -/

section Cross

/-- `cross(a, b)` of the Rust code. -/
def cross3 (u v : EuclideanSpace ℝ (Fin 3)) : EuclideanSpace ℝ (Fin 3) :=
  !₂[u 1 * v 2 - u 2 * v 1, u 2 * v 0 - u 0 * v 2, u 0 * v 1 - u 1 * v 0]

theorem inner3 (u v : EuclideanSpace ℝ (Fin 3)) : ⟪u, v⟫ = u 0 * v 0 + u 1 * v 1 + u 2 * v 2 := by
  simp [PiLp.inner_apply, Fin.sum_univ_three, mul_comm]

/-- `N x − (PC − QB) u − (QA − PB) z = (x·(u × z)) (u × z)` in `ℝ³` (`N = AC − B²`): the
component of `x` orthogonal to the plane of `u, z`. -/
theorem plane_dist_key (x u z : EuclideanSpace ℝ (Fin 3)) :
    (⟪u, u⟫ * ⟪z, z⟫ - ⟪u, z⟫ ^ 2) • x - (⟪x, u⟫ * ⟪z, z⟫ - ⟪x, z⟫ * ⟪u, z⟫) • u -
      (⟪x, z⟫ * ⟪u, u⟫ - ⟪x, u⟫ * ⟪u, z⟫) • z = ⟪x, cross3 u z⟫ • cross3 u z := by
  simp only [inner3]
  ext i
  fin_cases i <;> simp [cross3] <;> ring

/-- Lagrange's identity `(u × v)·(w × x) = (u·w)(v·x) − (u·x)(v·w)`. -/
theorem inner_cross3 (u v w x : EuclideanSpace ℝ (Fin 3)) :
    ⟪cross3 u v, cross3 w x⟫ = ⟪u, w⟫ * ⟪v, x⟫ - ⟪u, x⟫ * ⟪v, w⟫ := by
  simp only [inner3, cross3]; simp; ring

/-- **The cross-product scalars of `segment_segment_distance`**: `|d₁ × d₂|² = ae − b²` and
`((p₂ − p₁) × d₂)·(d₁ × d₂) / |d₁ × d₂|² = (bf − ce)/(ae − b²) = s*`.

[formula: rc-closest-points] -/
theorem cross_formulas_segSeg (p1 q1 p2 q2 : EuclideanSpace ℝ (Fin 3)) :
    ⟪cross3 (q1 - p1) (q2 - p2), cross3 (q1 - p1) (q2 - p2)⟫ = ssDen p1 q1 p2 q2 ∧
      ⟪cross3 (p2 - p1) (q2 - p2), cross3 (q1 - p1) (q2 - p2)⟫ /
          ⟪cross3 (q1 - p1) (q2 - p2), cross3 (q1 - p1) (q2 - p2)⟫ = ssS p1 q1 p2 q2 := by
  have e : p2 - p1 = -(p1 - p2) := by abel
  refine ⟨?_, ?_⟩
  · rw [inner_cross3]; unfold ssDen ssA ssE ssB; rw [real_inner_comm (q1 - p1) (q2 - p2)]; ring
  · rw [inner_cross3, inner_cross3, e, inner_neg_left, inner_neg_left]
    unfold ssS ssDen ssA ssE ssB ssC ssF
    rw [real_inner_comm (q1 - p1) (q2 - p2), real_inner_comm (p1 - p2) (q1 - p1),
      real_inner_comm (p1 - p2) (q2 - p2)]
    congr 1 <;> ring

/-- **The cross-product formulas of `point_triangle_distance`**: with `n = ab × ac`,
`|n|² = AC − B²`, `(ap × ac)·n/|n|² = v*`, `(ab × ap)·n/|n|² = w*`, and for `|n|² > 0` the plane
distance `|ap·n|/|n|` is `‖p − π‖`.

[formula: rc-closest-points] -/
theorem cross_formulas_ptTri (p a b c : EuclideanSpace ℝ (Fin 3)) :
    ⟪cross3 (b - a) (c - a), cross3 (b - a) (c - a)⟫ = tN a b c ∧
      ⟪cross3 (p - a) (c - a), cross3 (b - a) (c - a)⟫ /
          ⟪cross3 (b - a) (c - a), cross3 (b - a) (c - a)⟫ = tV p a b c ∧
      ⟪cross3 (b - a) (p - a), cross3 (b - a) (c - a)⟫ /
          ⟪cross3 (b - a) (c - a), cross3 (b - a) (c - a)⟫ = tW p a b c ∧
      (0 < tN a b c → |⟪p - a, cross3 (b - a) (c - a)⟫| /
          √⟪cross3 (b - a) (c - a), cross3 (b - a) (c - a)⟫ =
        ‖p - triPt a b c (tV p a b c) (tW p a b c)‖) := by
  have hN : ⟪cross3 (b - a) (c - a), cross3 (b - a) (c - a)⟫ = tN a b c := by
    rw [inner_cross3]; unfold tN tA tB tC; rw [real_inner_comm (b - a) (c - a)]; ring
  refine ⟨hN, ?_, ?_, fun hpos => ?_⟩
  · rw [inner_cross3, hN]; unfold tV tP tQ tB tC
    rw [real_inner_comm (b - a) (c - a)]
  · rw [inner_cross3, hN]; unfold tW tP tQ tA tB
    congr 1; ring
  · -- `N (p − π) = (ap·n) n`, a polynomial identity in the coordinates
    set n := cross3 (b - a) (c - a)
    have hv : tV p a b c * tN a b c = tP p a b * tC a c - tQ p a c * tB a b c :=
      div_mul_cancel₀ _ hpos.ne'
    have hw : tW p a b c * tN a b c = tQ p a c * tA a b - tP p a b * tB a b c :=
      div_mul_cancel₀ _ hpos.ne'
    have key : tN a b c • (p - triPt a b c (tV p a b c) (tW p a b c)) = ⟪p - a, n⟫ • n := by
      rw [triPt_sub]
      have : tN a b c • ((p - a) - tV p a b c • (b - a) - tW p a b c • (c - a)) =
          tN a b c • (p - a) - (tV p a b c * tN a b c) • (b - a) -
            (tW p a b c * tN a b c) • (c - a) := by module
      rw [this, hv, hw]
      unfold tN tA tB tC tP tQ
      exact plane_dist_key (p - a) (b - a) (c - a)
    have hn : ‖n‖ = √(tN a b c) := by rw [← hN, norm_eq_sqrt_real_inner]
    have := congrArg norm key
    rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos hpos, hn] at this
    have hs : 0 < √(tN a b c) := Real.sqrt_pos.2 hpos
    have hsq : √(tN a b c) * √(tN a b c) = tN a b c := Real.mul_self_sqrt hpos.le
    have h3 : √(tN a b c) * (√(tN a b c) * ‖p - triPt a b c (tV p a b c) (tW p a b c)‖) =
        √(tN a b c) * |⟪p - a, n⟫| := by rw [← mul_assoc, hsq, this, mul_comm]
    have h4 := mul_left_cancel₀ hs.ne' h3
    rw [hN, div_eq_iff hs.ne']; linarith [h4]

end Cross

end RC

end DDG
