import DdgOracle.Formulas
import DdgProofs.Optimization.PiecewiseAffine
import DdgProofs.SoftThreshold
import DdgProofs.Triangle

/-!
# The oracle's scalar-generic formulas at the real numbers

`lean/oracle` writes every formula once, against the interface `Ddg.Scalar`, and runs it at
`Float`. This file supplies the instance for `ℝ` and identifies the oracle's scalar and
per-vertex definitions (namespace `Ddg`) with the math-level definitions of this package
(namespace `DDG`), so that the theorems can be read as statements about the code that runs.

The identifications hold for exact real arithmetic; nothing is claimed about rounding.
-/

noncomputable section

namespace DDG

namespace Oracle

/-- The real numbers as a `Ddg.Scalar`: the obvious instance (`sqrt = Real.sqrt`,
`atan = Real.arctan`, `pi = Real.pi`, order and lattice operations those of `ℝ`).

The instance is `scoped`: it is active inside `namespace DDG.Oracle` and after
`open DDG.Oracle`, and nowhere else, so importing this package does not change how
`+`, `/`, `<`, `max` on `ℝ` are elaborated in a downstream file. -/
scoped instance instScalarReal : Ddg.Scalar ℝ where
  ofNat := fun n => (n : ℝ)
  sqrt := Real.sqrt
  atan := Real.arctan
  pi := Real.pi
  decLt := fun a b => Classical.dec (a < b)
  decLe := fun a b => Classical.dec (a ≤ b)

/-! ## Numerals and elementary functions -/

theorem ofNat_eq (n : ℕ) : (@OfNat.ofNat ℝ n Ddg.instOfNatScalar) = (n : ℝ) := rfl

@[scoped simp] theorem zero_eq : (@OfNat.ofNat ℝ 0 Ddg.instOfNatScalar) = (0 : ℝ) := by
  rw [ofNat_eq]; norm_num
@[scoped simp] theorem one_eq : (@OfNat.ofNat ℝ 1 Ddg.instOfNatScalar) = (1 : ℝ) := by
  rw [ofNat_eq]; norm_num
@[scoped simp] theorem two_eq : (@OfNat.ofNat ℝ 2 Ddg.instOfNatScalar) = (2 : ℝ) := by
  rw [ofNat_eq]
@[scoped simp] theorem three_eq : (@OfNat.ofNat ℝ 3 Ddg.instOfNatScalar) = (3 : ℝ) := by
  rw [ofNat_eq]
@[scoped simp] theorem four_eq : (@OfNat.ofNat ℝ 4 Ddg.instOfNatScalar) = (4 : ℝ) := by
  rw [ofNat_eq]

@[scoped simp] theorem sqrt_eq (x : ℝ) : Ddg.Scalar.sqrt x = √x := rfl
@[scoped simp] theorem atan_eq (x : ℝ) : Ddg.Scalar.atan x = Real.arctan x := rfl
@[scoped simp] theorem pi_eq : (Ddg.Scalar.pi : ℝ) = Real.pi := rfl
@[scoped simp] theorem abs_eq (x : ℝ) : Ddg.Scalar.abs x = |x| := (abs_eq_max_neg (a := x)).symm

/-! ## Soft threshold, clamp, bounded threshold -/

/-- The oracle's soft threshold is `DDG.softThreshold`. -/
theorem softThreshold_eq (t x : ℝ) : Ddg.softThreshold t x = softThreshold t x := by
  simp [Ddg.softThreshold, softThreshold]

/-- The oracle's `clamp lo hi x = min hi (max lo x)` is the clamp to `[lo, hi]` (for `lo ≤ hi`). -/
theorem clamp_eq {lo hi : ℝ} (h : lo ≤ hi) (x : ℝ) :
    Ddg.clamp lo hi x = clampE (lo : EReal) (hi : EReal) x := by
  unfold Ddg.clamp clampE
  rcases le_total x lo with h1 | h1
  · rw [max_eq_left h1, min_eq_right h,
      min_eq_right (by exact_mod_cast le_trans h1 h : ((x : ℝ) : EReal) ≤ hi),
      max_eq_left (by exact_mod_cast h1 : ((x : ℝ) : EReal) ≤ lo), EReal.toReal_coe]
  · rw [max_eq_right h1]
    rcases le_total x hi with h2 | h2
    · rw [min_eq_right h2, min_eq_right (by exact_mod_cast h2 : ((x : ℝ) : EReal) ≤ hi),
        max_eq_right (by exact_mod_cast h1 : (lo : EReal) ≤ ((x : ℝ) : EReal)), EReal.toReal_coe]
    · rw [min_eq_left h2, min_eq_left (by exact_mod_cast h2 : (hi : EReal) ≤ ((x : ℝ) : EReal)),
        max_eq_right (by exact_mod_cast h : (lo : EReal) ≤ hi), EReal.toReal_coe]

/-- The oracle's `InSubdiff t lo hi μ φ` is the subgradient inclusion
`φ ∈ ∂(t|·| + ι_[lo,hi])(μ)` (for `t ≥ 0`, `lo ≤ 0 ≤ hi`). -/
theorem inSubdiff_iff {t lo hi : ℝ} (ht : 0 ≤ t) (hlo : lo ≤ 0) (hhi : 0 ≤ hi) (μ φ : ℝ) :
    Ddg.InSubdiff t lo hi μ φ ↔
      IsSubgradientWithinAt (fun x => t * |x|) (boundSet (lo : EReal) (hi : EReal)) φ μ := by
  rw [isSubgradientWithinAt_abs_iff ht (by exact_mod_cast hlo) (by exact_mod_cast hhi)]
  unfold Ddg.InSubdiff boundSet
  simp only [zero_eq, Set.mem_ofPred_eq, EReal.coe_le_coe_iff, EReal.coe_lt_coe_iff, and_assoc]
  have e1 : (if μ < 0 then -t else t) = (if 0 ≤ μ then t else -t) := by
    split_ifs with h1 h2 h2 <;> first | rfl | (exfalso; linarith)
  rw [e1]

/-- **Fixed-point form on the oracle's definitions** (paper Eq. 19 with bounds): for `t ≥ 0`,
`lo ≤ 0 ≤ hi`, `InSubdiff t lo hi μ φ  ⟺  μ = boundedThreshold t lo hi (φ + μ)`.

[formula: cone-fixed-point-bounded] -/
theorem inSubdiff_iff_eq_boundedThreshold {t lo hi : ℝ} (ht : 0 ≤ t) (hlo : lo ≤ 0)
    (hhi : 0 ≤ hi) (μ φ : ℝ) :
    Ddg.InSubdiff t lo hi μ φ ↔ μ = Ddg.boundedThreshold t lo hi (φ + μ) := by
  rw [inSubdiff_iff ht hlo hhi,
    isSubgradientWithinAt_abs_iff_eq_clampE ht (by exact_mod_cast hlo) (by exact_mod_cast hhi)]
  unfold Ddg.boundedThreshold
  rw [clamp_eq (le_trans hlo hhi), softThreshold_eq]

/-- **Regularized cone update on the oracle's definitions**: for `t ≥ 0`, `k > 0`, `lo ≤ 0 ≤ hi`,
`InSubdiff t lo hi μ (φ − μ/k)  ⟺  μ = regularizedMeasure k t lo hi φ`.
With `k = m/γ` this is the Stage I formula `μ = clamp((m/γ)·P_t(φ))` including the mass factor.

[formula: regularized-cone-update] -/
theorem inSubdiff_sub_div_iff_eq_regularizedMeasure {t k lo hi : ℝ} (ht : 0 ≤ t) (hk : 0 < k)
    (hlo : lo ≤ 0) (hhi : 0 ≤ hi) (μ φ : ℝ) :
    Ddg.InSubdiff t lo hi μ (φ - μ / k) ↔ μ = Ddg.regularizedMeasure k t lo hi φ := by
  rw [inSubdiff_iff ht hlo hhi]
  have h := isSubgradientWithinAt_abs_sub_iff_eq_clampE ht (inv_pos.mpr hk)
    (lo := (lo : EReal)) (hi := (hi : EReal)) (by exact_mod_cast hlo) (by exact_mod_cast hhi) φ μ
  have e1 : φ - k⁻¹ * μ = φ - μ / k := by field_simp
  have e2 : softThreshold t φ / k⁻¹ = k * softThreshold t φ := by field_simp
  rw [e1, e2] at h
  rw [h]
  unfold Ddg.regularizedMeasure
  rw [clamp_eq (le_trans hlo hhi), softThreshold_eq]

/-- **The oracle's regularized fixed point with `μ` kept as an unknown** (not in the paper):
for `t ≥ 0`, `m > 0`, `γ ≥ 0`, `lo ≤ 0 ≤ hi` and `c = m/(m + γ)`,
`μ = regularizedMeasure c t lo hi (φ + μ)  ⟺  InSubdiff t lo hi μ (φ − (γ/m)·μ)`,
i.e. it is equivalent to the optimality condition `φ − (γ/m)μ ∈ ∂g(μ)` of the regularized
problem; for `γ = 0` it is the exact cone condition.

[formula: regularized-cone-update] -/
theorem eq_regularizedMeasure_add_iff {t m γ lo hi : ℝ} (ht : 0 ≤ t) (hm : 0 < m) (hγ : 0 ≤ γ)
    (hlo : lo ≤ 0) (hhi : 0 ≤ hi) (μ φ : ℝ) :
    μ = Ddg.regularizedMeasure (m / (m + γ)) t lo hi (φ + μ) ↔
      Ddg.InSubdiff t lo hi μ (φ - γ / m * μ) := by
  have hc : 0 < m / (m + γ) := div_pos hm (by linarith)
  rw [← inSubdiff_sub_div_iff_eq_regularizedMeasure ht hc hlo hhi μ (φ + μ)]
  have e : φ + μ - μ / (m / (m + γ)) = φ - γ / m * μ := by
    have : m + γ ≠ 0 := by linarith
    field_simp
    ring
  rw [e]

/-! ## The oracle's active sets (`coneState`) -/

/-- **The oracle's pieces are correct**: for `t ≥ 0` and `lo ≤ hi`, the map
`x ↦ clamp lo hi (k·P_t(x))` equals, at every `x`, the value of the piece `coneState k t lo hi x`
at `x`.

[formula: soft-threshold-active-sets] -/
theorem regularizedMeasure_eq_value {t lo hi : ℝ} (ht : 0 ≤ t) (hlh : lo ≤ hi) (k x : ℝ) :
    Ddg.regularizedMeasure k t lo hi x =
      (Ddg.coneState k t lo hi x).value k t lo hi x := by
  unfold Ddg.regularizedMeasure Ddg.coneState Ddg.clamp
  rw [softThreshold_eq]
  split_ifs with h1 h2 h3 h4
  · -- upper
    simp only [Ddg.ConeState.value]
    exact min_eq_left (le_trans h1 (le_max_right _ _))
  · -- lower
    simp only [Ddg.ConeState.value]
    rw [max_eq_left h2, min_eq_right hlh]
  · -- pos
    simp only [Ddg.ConeState.value]
    have hs : softThreshold t x = x - t := by
      rcases h3 with h3 | ⟨h3, h3'⟩
      · exact softThreshold_of_gt ht h3
      · simp only [zero_eq] at h3 h3'
        have : t = 0 := le_antisymm h3 ht
        subst this
        rw [softThreshold_zero]; ring
    rw [max_eq_right (not_le.mp h2).le, min_eq_right (not_le.mp h1).le, hs]
  · -- neg
    simp only [Ddg.ConeState.value]
    have hs : softThreshold t x = x + t := softThreshold_of_lt ht h4
    rw [max_eq_right (not_le.mp h2).le, min_eq_right (not_le.mp h1).le, hs]
  · -- zero
    simp only [Ddg.ConeState.value]
    unfold Ddg.ActivePos at h3
    unfold Ddg.ActiveNeg at h4
    push Not at h3 h4
    simp only [zero_eq] at h3 ⊢
    have hx : x ≤ t := h3.1
    have hs : softThreshold t x = 0 := softThreshold_of_abs_le (abs_le.mpr ⟨h4, hx⟩)
    rw [hs, mul_zero] at h1 h2 ⊢
    rw [max_eq_right (not_le.mp h2).le, min_eq_right (not_le.mp h1).le]

/-- On a fixed piece the value is affine in `x` with slope `k·indicator`. -/
theorem value_affine (k t lo hi x x' : ℝ) (σ : Ddg.ConeState) :
    σ.value k t lo hi x' = σ.value k t lo hi x + k * σ.indicator * (x' - x) := by
  cases σ <;> simp [Ddg.ConeState.value, Ddg.ConeState.indicator] <;> ring

/-- **The oracle's map `x ↦ clamp lo hi (k·P_t(x))` is piecewise affine with respect to
`coneState`**, with generalized derivative `k·indicator`: the fact behind the stopping rule
"active sets unchanged ⟹ exact solution".

[formula: soft-threshold-active-sets] -/
theorem regularizedMeasure_isPiecewiseAffine {t lo hi : ℝ} (ht : 0 ≤ t) (hlh : lo ≤ hi)
    (k : ℝ) :
    IsPiecewiseAffine (Ddg.regularizedMeasure k t lo hi) (Ddg.coneState k t lo hi)
      (fun x => k * (Ddg.coneState k t lo hi x).indicator) := by
  intro x x' h
  rw [regularizedMeasure_eq_value ht hlh, regularizedMeasure_eq_value ht hlh, h]
  exact value_affine k t lo hi x x' _

/-- When `lo < hi`, a piece and its canonical name have the same value everywhere. -/
theorem value_canon_of_lt {t lo hi : ℝ} (ht : 0 ≤ t) (h : lo < hi) (k x : ℝ)
    (σ : Ddg.ConeState) : (σ.canon t lo hi).value k t lo hi x = σ.value k t lo hi x := by
  unfold Ddg.ConeState.canon
  rw [ite_eq_left h]
  cases σ with
  | zero => rfl
  | pos => rfl
  | neg =>
    dsimp only
    simp only [zero_eq]
    split_ifs with h2
    · rfl
    · have ht0 : t = 0 := le_antisymm (not_lt.mp h2) ht
      simp [Ddg.ConeState.value, ht0]
  | lower =>
    dsimp only
    simp only [zero_eq]
    split_ifs with h2
    · rfl
    · push Not at h2
      have h0 : lo = 0 := le_antisymm h2.1 h2.2
      simp [Ddg.ConeState.value, h0]
  | upper =>
    dsimp only
    simp only [zero_eq]
    split_ifs with h2
    · rfl
    · push Not at h2
      have h0 : hi = 0 := le_antisymm h2.1 h2.2
      simp [Ddg.ConeState.value, h0]

/-- The canonical name of the piece a point lies on has the same value at that point. -/
theorem value_canon_coneState {t lo hi : ℝ} (ht : 0 ≤ t) (hlh : lo ≤ hi) (k x : ℝ) :
    ((Ddg.coneState k t lo hi x).canon t lo hi).value k t lo hi x =
      (Ddg.coneState k t lo hi x).value k t lo hi x := by
  by_cases h1 : lo < hi
  · exact value_canon_of_lt ht h1 k x _
  · have heq : lo = hi := le_antisymm hlh (not_lt.mp h1)
    unfold Ddg.ConeState.canon
    rw [ite_eq_right h1]
    unfold Ddg.coneState
    by_cases h2 : hi ≤ k * Ddg.softThreshold t x
    · rw [ite_eq_left h2]
    · have h3 : k * Ddg.softThreshold t x ≤ lo := by rw [heq]; exact le_of_not_ge h2
      rw [ite_eq_right h2, ite_eq_left h3]
      simp [Ddg.ConeState.value, heq]

/-- When `lo < hi`, a piece and its canonical name have the same slope indicator. -/
theorem indicator_canon_of_lt {t lo hi : ℝ} (h : lo < hi) (σ : Ddg.ConeState) :
    ((σ.canon t lo hi).indicator : ℝ) = σ.indicator := by
  unfold Ddg.ConeState.canon
  rw [ite_eq_left h]
  cases σ with
  | zero => rfl
  | pos => rfl
  | neg =>
    dsimp only
    simp only [zero_eq]
    split_ifs <;> rfl
  | lower =>
    dsimp only
    simp only [zero_eq]
    split_ifs <;> rfl
  | upper =>
    dsimp only
    simp only [zero_eq]
    split_ifs <;> rfl

/-- At a vertex with `lo = hi` (the oracle's encoding of a forbidden or prescribed vertex) the
map `x ↦ clamp lo hi (k·P_t(x))` is the constant `hi`, and the piece of every point is one of the
two clamped ones, with slope indicator `0`. -/
theorem regularizedMeasure_of_eq {t lo hi : ℝ} (h : lo = hi) (k x : ℝ) :
    Ddg.regularizedMeasure k t lo hi x = hi ∧
      ((Ddg.coneState k t lo hi x).indicator : ℝ) = 0 := by
  subst h
  constructor
  · unfold Ddg.regularizedMeasure Ddg.clamp
    exact min_eq_left (le_max_left _ _)
  · unfold Ddg.coneState
    by_cases h2 : lo ≤ k * Ddg.softThreshold t x
    · rw [ite_eq_left h2]; simp [Ddg.ConeState.indicator]
    · rw [ite_eq_right h2, ite_eq_left (le_of_not_ge h2)]; simp [Ddg.ConeState.indicator]

/-- **The oracle's map `x ↦ clamp lo hi (k·P_t(x))` is affine between two points whose pieces
have the same canonical name** (`ConeState.canon`), with slope `k·indicator` of the piece of the
first point. Canonical names identify pieces that impose the same equation: the two clamped
pieces of a vertex with `lo = hi`, a clamped piece whose bound is `0` and the piece `zero`, and
the pieces `pos`, `neg` of a free site (`t = 0`). Equality of the raw pieces is a special case.

[formula: soft-threshold-active-sets] -/
theorem regularizedMeasure_affine_of_canon_eq {t lo hi : ℝ} (ht : 0 ≤ t) (hlh : lo ≤ hi)
    (k x x' : ℝ)
    (h : (Ddg.coneState k t lo hi x').canon t lo hi = (Ddg.coneState k t lo hi x).canon t lo hi) :
    Ddg.regularizedMeasure k t lo hi x' =
      Ddg.regularizedMeasure k t lo hi x +
        k * (Ddg.coneState k t lo hi x).indicator * (x' - x) := by
  by_cases hlt : lo < hi
  · rw [regularizedMeasure_eq_value ht hlh k x', ← value_canon_coneState ht hlh k x', h,
      value_canon_of_lt ht hlt k x' _, value_affine k t lo hi x x' _,
      ← regularizedMeasure_eq_value ht hlh k x]
  · have heq : lo = hi := le_antisymm hlh (not_lt.mp hlt)
    rw [(regularizedMeasure_of_eq heq k x').1, (regularizedMeasure_of_eq heq k x).1,
      (regularizedMeasure_of_eq heq k x).2]
    ring

/-- `boundedThreshold` is `regularizedMeasure` with `k = 1`. -/
theorem boundedThreshold_eq_regularizedMeasure (t lo hi x : ℝ) :
    Ddg.boundedThreshold t lo hi x = Ddg.regularizedMeasure 1 t lo hi x := by
  simp [Ddg.boundedThreshold, Ddg.regularizedMeasure]

/-! ## Triangle formulas -/

/-- The oracle's Heron area is `DDG.triArea`. -/
theorem triArea_eq (a b c : ℝ) : Ddg.triArea a b c = triArea a b c := by
  simp [Ddg.triArea, triArea]

/-- The oracle's cotangent is `DDG.cornerCot`. -/
theorem cornerCot_eq (a b c : ℝ) : Ddg.cornerCot a b c = cornerCot a b c := by
  unfold Ddg.cornerCot cornerCot
  rw [triArea_eq]
  simp only [four_eq]
  ring_nf

/-- **The oracle's corner angle** (half-angle form with `atan`) is the angle of the law of
cosines, for every nondegenerate triangle.

[formula: corner-angle] -/
theorem cornerAngle_eq {a b c : ℝ} (h : IsTriangle a b c) :
    Ddg.cornerAngle a b c = cornerAngle a b c := by
  rw [cornerAngle_eq_two_arctan h]
  simp [Ddg.cornerAngle]

/-- The oracle's `dist` on `Point ℝ` is the Euclidean distance in `ℝ³`.

[formula: edge-length-euclidean] -/
theorem dist_eq (p q : Ddg.Point ℝ) :
    Ddg.dist p q = dist (!₂[p.x, p.y, p.z] : EuclideanSpace ℝ (Fin 3)) !₂[q.x, q.y, q.z] := by
  rw [dist_eq_sqrt_sum_sq]
  simp [Ddg.dist, sq]

end Oracle

end DDG
