import Mathlib
import DdgProofs.RepulsiveCurves.Ccd

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
* `primDist_segSeg`, `primDist_ptTri`: these are `primDist` of `Ccd.lean` for the corner positions
  at time `τ`; `segSeg_advance_safe`, `ptTri_advance_safe`: the advancement of `ccd.rs` with the
  computed distances keeps the primitives apart (with corner speeds measured relative to any
  common velocity, `primDist_sub_common`).
-/

noncomputable section

namespace DdgProofs

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
        unfold segSegDist; rw [if_pos ⟨ha, he, hD, hi⟩]
      · unfold segSegDist; rw [if_neg (fun h => hi h.2.2.2)]; exact min_le_left _ _
    · unfold segSegDist; rw [if_neg (fun h => hc ⟨h.1, h.2.1, h.2.2.1⟩)]
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
    unfold segSegDist; rw [if_pos ⟨ha, he, this, hi⟩]
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
      · unfold segSegDist; rw [if_pos ⟨ha, he, hD, hi⟩]; exact min_le_right _ _
      · obtain ⟨s, t, h, e⟩ := segSegDist_attained (p1 := p1) (q1 := q1) (p2 := p2) (q2 := q2)
        rw [e]
        refine norm_le_of_sq_le (norm_nonneg _) (norm_nonneg _) ?_
        rw [ssVec_sq_split hD.ne' s t]
        nlinarith [sq_nonneg ‖(s - ssS p1 q1 p2 q2) • (q1 - p1) -
          (t - ssT p1 q1 p2 q2) • (q2 - p2)‖]
    by_cases hthr : ε * ssA p1 q1 * ssE p2 q2 < ssDen p1 q1 p2 q2
    · have : segSegThr p1 q1 p2 q2 ε = segSegDist p1 q1 p2 q2 := by
        unfold segSegThr segSegDist; rw [if_pos ⟨ha, he, hthr, hi⟩, if_pos ⟨ha, he, hD, hi⟩]
      rw [this]; nlinarith
    · have hthr' : ssDen p1 q1 p2 q2 ≤ ε * ssA p1 q1 * ssE p2 q2 := not_lt.1 hthr
      have hval : segSegThr p1 q1 p2 q2 ε = ssBest p1 q1 p2 q2 := by
        unfold segSegThr; rw [if_neg (fun h => hthr h.2.2.1)]
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
        rw [if_neg hc]
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
    unfold ptTriDist; rw [if_pos hN]
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
    unfold ptTriDist; rw [if_neg hN]
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

end RC

end DdgProofs
