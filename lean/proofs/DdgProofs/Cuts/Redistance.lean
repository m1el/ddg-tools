import DdgProofs.Cuts.Crossing

/-!
# Cut flow: the local Hopf–Lax update of the redistancing

`crates/ddg-variational-cuts/SPEC.md` §3.12, D10 (Bornemann, Rasch, *Finite-element
discretization of static Hamilton–Jacobi equations based on a local variational principle*,
Comput. Vis. Sci. 9, 2006). For a triangle `i j k` with accepted values `T_j`, `T_k` the tentative
value at `i` is

`min_{λ ∈ [0,1]}  λ T_j + (1 − λ) T_k + ‖x_i − (λ x_j + (1 − λ) x_k)‖`.

**`cut-redistance`**: the oracle's closed form `Ddg.hopfLaxTriangle`, evaluated at the three
side lengths, **is** this minimum:

* `hopfLaxTriangle_isLeast`: for points `x_i, x_j, x_k` of any real inner product space (the
  plane of the triangle; no nondegeneracy assumption, `x_j = x_k` included),
  `hopfLaxTriangle (dist x_i x_j) (dist x_i x_k) (dist x_j x_k) T_j T_k` is the least value of the
  objective on `[0, 1]`: it is attained (`hopfLaxTriangle_mem`) and it is a lower bound
  (`hopfLaxTriangle_le`);
* `hopfLaxTriangle_isLeast_intrinsic`: the same from the side lengths alone, for lengths that
  satisfy the strict triangle inequalities; the distance from `x_i` to the point of the edge is
  the oracle's `baryLength` (`cut-barycentric-length`);
* `hopfLaxTriangle_le_edge`: it never exceeds the two edge terms `hopfLaxEdge T_j ℓ_ij`,
  `hopfLaxEdge T_k ℓ_ik`; `min_le_hopfLaxTriangle`: it is at least `min(T_j, T_k)`.

The proof: the objective is convex in `λ` (`hopfLaxObjective_convex`); if `|T_j − T_k| < ℓ_jk` its
minimum over the whole line is at the stationary point of the closed form
(`hl_stationary_le`, a Cauchy–Schwarz inequality), and the minimum over `[0, 1]` of a convex
function is at the clamped unconstrained minimizer (`clamp_min_of_convex`); otherwise the
objective is monotone (`hopfLaxObjective_endpoint_le`) and the minimum is at an end point.

What is **not** covered: the marching order and the single pass (a definition, SPEC §3.12), and
that the scheme converges to the geodesic distance (`cut-fast-marching-convergence`, cited).
-/

noncomputable section

namespace DdgProofs

namespace Cuts

open scoped DdgProofs.Oracle
open Real

/-! ## The oracle's closed form at `ℝ` -/

/-- The oracle's edge term at `ℝ`: `T_j + ℓ_ij`. -/
theorem hopfLaxEdge_eq (Tj lij : ℝ) : Ddg.hopfLaxEdge Tj lij = Tj + lij := rfl

/-- The radicand of the closed form: `|d|² + 2λ d·e + λ²|e|²` with `d = x_k − x_i`,
`e = x_j − x_k`, written with the side lengths exactly as the oracle does. -/
def hlQ (lij lik ljk lam : ℝ) : ℝ :=
  lik * lik + 2 * lam * ((lij * lij - lik * lik - ljk * ljk) / 2) + lam * lam * (ljk * ljk)

/-- The oracle's local function `value`. -/
def hlValue (lij lik ljk Tj Tk lam : ℝ) : ℝ :=
  Tk + lam * (Tj - Tk) + √(max 0 (hlQ lij lik ljk lam))

/-- The oracle's stationary point `λ* = λ₀ − ΔT h / (a √(a² − ΔT²))`. -/
def hlLam (lij lik ljk Tj Tk : ℝ) : ℝ :=
  -((lij * lij - lik * lik - ljk * ljk) / 2) / (ljk * ljk) -
    (Tj - Tk) * √(max 0 (lik * lik - (lij * lij - lik * lik - ljk * ljk) / 2 *
      ((lij * lij - lik * lik - ljk * ljk) / 2) / (ljk * ljk))) /
      (ljk * √(ljk * ljk - (Tj - Tk) * (Tj - Tk)))

/-- **The oracle's `hopfLaxTriangle` at `ℝ`**, with its `let`s named. -/
theorem hopfLaxTriangle_eq (lij lik ljk Tj Tk : ℝ) :
    Ddg.hopfLaxTriangle lij lik ljk Tj Tk =
      if |Tj - Tk| < ljk then hlValue lij lik ljk Tj Tk (max 0 (min 1 (hlLam lij lik ljk Tj Tk)))
      else min (hlValue lij lik ljk Tj Tk 0) (hlValue lij lik ljk Tj Tk 1) := by
  unfold Ddg.hopfLaxTriangle hlValue hlQ hlLam
  simp only [Oracle.sqrt_eq, Oracle.abs_eq, Oracle.zero_eq, Oracle.one_eq, Oracle.two_eq]

/-! ## One-dimensional facts -/

/-- A convex function does not exceed its value at `y` between a global minimizer and `y`. -/
theorem convex_between {F : ℝ → ℝ}
    (hconv : ∀ x y θ : ℝ, 0 ≤ θ → θ ≤ 1 → F ((1 - θ) * x + θ * y) ≤ (1 - θ) * F x + θ * F y)
    {m y θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ 1) (hmy : F m ≤ F y) :
    F ((1 - θ) * m + θ * y) ≤ F y := by
  have hc := hconv m y θ h0 h1
  have hp := mul_nonneg (sub_nonneg.mpr h1) (sub_nonneg.mpr hmy)
  linarith

/-- **The minimum of a convex function over `[0, 1]` is at the clamped unconstrained
minimizer.** -/
theorem clamp_min_of_convex {F : ℝ → ℝ}
    (hconv : ∀ x y θ : ℝ, 0 ≤ θ → θ ≤ 1 → F ((1 - θ) * x + θ * y) ≤ (1 - θ) * F x + θ * F y)
    {m : ℝ} (hm : ∀ x, F m ≤ F x) {lam : ℝ} (h0 : 0 ≤ lam) (h1 : lam ≤ 1) :
    F (max 0 (min 1 m)) ≤ F lam := by
  by_cases hm0 : m ≤ 0
  · rw [min_eq_right (by linarith : m ≤ 1), max_eq_left hm0]
    rcases eq_or_lt_of_le hm0 with heq | hlt
    · rw [← heq]
      exact hm lam
    · have hd : 0 < lam - m := by linarith
      have hθ0 : 0 ≤ -m / (lam - m) := div_nonneg (by linarith) hd.le
      have hθ1 : -m / (lam - m) ≤ 1 := by
        rw [div_le_one hd]
        linarith
      have h := convex_between hconv hθ0 hθ1 (hm lam)
      have e : (1 - -m / (lam - m)) * m + -m / (lam - m) * lam = 0 := by
        have := hd.ne'
        field_simp
        ring
      rwa [e] at h
  · push Not at hm0
    by_cases hm1 : 1 ≤ m
    · rw [min_eq_left hm1, max_eq_right zero_le_one]
      rcases eq_or_lt_of_le hm1 with heq | hlt
      · rw [heq]
        exact hm lam
      · have hd : 0 < m - lam := by linarith
        have hθ0 : 0 ≤ (m - 1) / (m - lam) := div_nonneg (by linarith) hd.le
        have hθ1 : (m - 1) / (m - lam) ≤ 1 := by
          rw [div_le_one hd]
          linarith
        have h := convex_between hconv hθ0 hθ1 (hm lam)
        have e : (1 - (m - 1) / (m - lam)) * m + (m - 1) / (m - lam) * lam = 1 := by
          have := hd.ne'
          field_simp
          ring
        rwa [e] at h
    · push Not at hm1
      rw [min_eq_right hm1.le, max_eq_right hm0.le]
      exact hm lam

/-- **The stationary point of the closed form is a global minimizer** of
`λ ↦ λ ΔT + √(|d|² + 2λ d·e + λ² a²)` when `|ΔT| < a` and `h² = |d|² − (d·e)²/a² ≥ 0`.
With `s = a(λ − λ₀)` this is the Cauchy–Schwarz inequality
`(ΔT/a) s + √(s² + h²) ≥ h √(1 − ΔT²/a²)`, with equality at the stationary point. -/
theorem hl_stationary_le {dd de a dT : ℝ} (ha : 0 < a) (hh : 0 ≤ dd - de * de / (a * a))
    (hdT : |dT| < a) (lam : ℝ) :
    (-de / (a * a) - dT * √(dd - de * de / (a * a)) / (a * √(a * a - dT * dT))) * dT +
        √(dd + 2 * (-de / (a * a) - dT * √(dd - de * de / (a * a)) /
            (a * √(a * a - dT * dT))) * de +
          (-de / (a * a) - dT * √(dd - de * de / (a * a)) / (a * √(a * a - dT * dT))) *
            (-de / (a * a) - dT * √(dd - de * de / (a * a)) / (a * √(a * a - dT * dT))) *
            (a * a)) ≤
      lam * dT + √(dd + 2 * lam * de + lam * lam * (a * a)) := by
  obtain ⟨hl, hu⟩ := abs_lt.mp hdT
  have hane : a ≠ 0 := ha.ne'
  have hpos : 0 < a * a - dT * dT := by nlinarith
  obtain ⟨R, hR⟩ : ∃ R, R = √(a * a - dT * dT) := ⟨_, rfl⟩
  obtain ⟨h, hh'⟩ : ∃ h, h = √(dd - de * de / (a * a)) := ⟨_, rfl⟩
  rw [← hR, ← hh']
  have hRpos : 0 < R := by rw [hR]; exact Real.sqrt_pos.mpr hpos
  have hRne : R ≠ 0 := hRpos.ne'
  have hR2 : R ^ 2 = a ^ 2 - dT ^ 2 := by
    rw [hR, Real.sq_sqrt hpos.le]
    ring
  have hh0 : 0 ≤ h := by rw [hh']; exact Real.sqrt_nonneg _
  have hh2 : h ^ 2 = dd - de * de / (a * a) := by rw [hh', Real.sq_sqrt hh]
  obtain ⟨w, hw⟩ : ∃ w, w = h / (a * R) := ⟨_, rfl⟩
  have hw0 : 0 ≤ w := by rw [hw]; positivity
  have hhw : h = w * a * R := by rw [hw]; field_simp
  obtain ⟨lam0, hlam0⟩ : ∃ lam0, lam0 = -de / (a * a) := ⟨_, rfl⟩
  have hlamS : -de / (a * a) - dT * h / (a * R) = lam0 - dT * w := by
    rw [hlam0, hw]
    ring
  rw [hlamS]
  have hQ : ∀ x : ℝ, dd + 2 * x * de + x * x * (a * a) = a ^ 2 * (x - lam0) ^ 2 + h ^ 2 := by
    intro x
    rw [hh2, hlam0]
    field_simp
    ring
  rw [hQ, hQ]
  have hS : a ^ 2 * (lam0 - dT * w - lam0) ^ 2 + h ^ 2 = (w * a ^ 2) ^ 2 := by
    rw [hhw]
    linear_combination (w ^ 2 * a ^ 2) * hR2
  rw [hS, Real.sqrt_sq (by positivity)]
  have key : (w * R ^ 2 - (lam - lam0) * dT) ^ 2 ≤ a ^ 2 * (lam - lam0) ^ 2 + h ^ 2 := by
    have e : a ^ 2 * (lam - lam0) ^ 2 + h ^ 2 - (w * R ^ 2 - (lam - lam0) * dT) ^ 2 =
        R ^ 2 * ((lam - lam0) + w * dT) ^ 2 := by
      rw [hhw]
      linear_combination (-((lam - lam0) ^ 2 + w ^ 2 * R ^ 2)) * hR2
    have : 0 ≤ R ^ 2 * ((lam - lam0) + w * dT) ^ 2 := by positivity
    linarith
  have h1 := (le_abs_self _).trans (Real.abs_le_sqrt key)
  have e2 : (lam0 - dT * w) * dT + w * a ^ 2 = lam * dT + (w * R ^ 2 - (lam - lam0) * dT) := by
    rw [hR2]
    ring
  linarith

/-! ## The objective for points of an inner product space -/

section Geometry

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The Hopf–Lax objective on the edge `jk`:
`λ T_j + (1 − λ) T_k + ‖x_i − (λ x_j + (1 − λ) x_k)‖`. -/
def hopfLaxObjective (xi xj xk : E) (Tj Tk lam : ℝ) : ℝ :=
  lam * Tj + (1 - lam) * Tk + ‖xi - (lam • xj + (1 - lam) • xk)‖

/-- The radicand of the closed form is the squared distance from `x_i` to the point
`λ x_j + (1 − λ) x_k` of the line `jk`. -/
theorem hlQ_eq_norm_sq (xi xj xk : E) (lam : ℝ) :
    hlQ (dist xi xj) (dist xi xk) (dist xj xk) lam = ‖xi - (lam • xj + (1 - lam) • xk)‖ ^ 2 := by
  have e1 : xi - xj = -((xk - xi) + (xj - xk)) := by abel
  have e2 : xi - xk = -(xk - xi) := by abel
  have e3 : xi - (lam • xj + (1 - lam) • xk) = -((xk - xi) + lam • (xj - xk)) := by module
  simp only [dist_eq_norm]
  rw [e1, e2, e3, norm_neg, norm_neg, norm_neg]
  generalize xk - xi = d
  generalize xj - xk = e
  have h1 := norm_add_sq_real d e
  have h2 := norm_add_sq_real d (lam • e)
  rw [real_inner_smul_right, norm_smul, Real.norm_eq_abs, mul_pow, sq_abs] at h2
  unfold hlQ
  rw [h2]
  linear_combination lam * h1

theorem hlQ_nonneg (xi xj xk : E) (lam : ℝ) :
    0 ≤ hlQ (dist xi xj) (dist xi xk) (dist xj xk) lam := by
  rw [hlQ_eq_norm_sq]
  positivity

/-- The oracle's `value` is the objective. -/
theorem hlValue_eq_objective (xi xj xk : E) (Tj Tk lam : ℝ) :
    hlValue (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk lam =
      hopfLaxObjective xi xj xk Tj Tk lam := by
  unfold hlValue hopfLaxObjective
  rw [hlQ_eq_norm_sq, max_eq_right (sq_nonneg _), Real.sqrt_sq (norm_nonneg _)]
  ring

/-- **The objective is convex in `λ`.** -/
theorem hopfLaxObjective_convex (xi xj xk : E) (Tj Tk : ℝ) (x y : ℝ) {θ : ℝ} (h0 : 0 ≤ θ)
    (h1 : θ ≤ 1) :
    hopfLaxObjective xi xj xk Tj Tk ((1 - θ) * x + θ * y) ≤
      (1 - θ) * hopfLaxObjective xi xj xk Tj Tk x + θ * hopfLaxObjective xi xj xk Tj Tk y := by
  unfold hopfLaxObjective
  have e : xi - (((1 - θ) * x + θ * y) • xj + (1 - ((1 - θ) * x + θ * y)) • xk) =
      (1 - θ) • (xi - (x • xj + (1 - x) • xk)) + θ • (xi - (y • xj + (1 - y) • xk)) := by
    module
  have hn : ‖xi - (((1 - θ) * x + θ * y) • xj + (1 - ((1 - θ) * x + θ * y)) • xk)‖ ≤
      (1 - θ) * ‖xi - (x • xj + (1 - x) • xk)‖ + θ * ‖xi - (y • xj + (1 - y) • xk)‖ := by
    rw [e]
    refine (norm_add_le _ _).trans_eq ?_
    rw [norm_smul, norm_smul, Real.norm_of_nonneg h0,
      Real.norm_of_nonneg (by linarith : (0 : ℝ) ≤ 1 - θ)]
  linarith

/-- The distance to a point of the line `jk` is Lipschitz in the parameter, with constant
`ℓ_jk`. -/
theorem norm_sub_point_le (xi xj xk : E) (x y : ℝ) :
    ‖xi - (x • xj + (1 - x) • xk)‖ ≤
      ‖xi - (y • xj + (1 - y) • xk)‖ + |x - y| * dist xj xk := by
  have e : xi - (x • xj + (1 - x) • xk) =
      (xi - (y • xj + (1 - y) • xk)) + (y - x) • (xj - xk) := by module
  rw [e, dist_eq_norm]
  calc ‖(xi - (y • xj + (1 - y) • xk)) + (y - x) • (xj - xk)‖
      ≤ ‖xi - (y • xj + (1 - y) • xk)‖ + ‖(y - x) • (xj - xk)‖ := norm_add_le _ _
    _ = ‖xi - (y • xj + (1 - y) • xk)‖ + |x - y| * ‖xj - xk‖ := by
        rw [norm_smul, Real.norm_eq_abs, abs_sub_comm]

/-- **If `|T_j − T_k| ≥ ℓ_jk` the objective is monotone**, so its minimum over `[0, 1]` is at an
end point. -/
theorem hopfLaxObjective_endpoint_le (xi xj xk : E) (Tj Tk : ℝ)
    (h : dist xj xk ≤ |Tj - Tk|) {lam : ℝ} (h0 : 0 ≤ lam) (h1 : lam ≤ 1) :
    min (hopfLaxObjective xi xj xk Tj Tk 0) (hopfLaxObjective xi xj xk Tj Tk 1) ≤
      hopfLaxObjective xi xj xk Tj Tk lam := by
  rcases le_abs'.mp h with h' | h'
  · refine (min_le_right _ _).trans ?_
    have hn := norm_sub_point_le xi xj xk 1 lam
    rw [abs_of_nonneg (by linarith : (0 : ℝ) ≤ 1 - lam)] at hn
    have hp := mul_nonneg (sub_nonneg.mpr h1)
      (show 0 ≤ -(Tj - Tk) - dist xj xk by linarith)
    unfold hopfLaxObjective
    linarith
  · refine (min_le_left _ _).trans ?_
    have hn := norm_sub_point_le xi xj xk 0 lam
    rw [abs_of_nonpos (by linarith : (0 : ℝ) - lam ≤ 0)] at hn
    have hp := mul_nonneg h0 (show 0 ≤ (Tj - Tk) - dist xj xk by linarith)
    unfold hopfLaxObjective
    linarith

/-- **If `|T_j − T_k| < ℓ_jk` the oracle's stationary point minimizes the objective over the
whole line.** -/
theorem hlValue_stationary_le (xi xj xk : E) (Tj Tk : ℝ) (h : |Tj - Tk| < dist xj xk)
    (lam : ℝ) :
    hlValue (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk
        (hlLam (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk) ≤
      hlValue (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk lam := by
  have ha : 0 < dist xj xk := lt_of_le_of_lt (abs_nonneg _) h
  have hane : dist xj xk ≠ 0 := ha.ne'
  have hq0 := hlQ_nonneg xi xj xk
  have hh : 0 ≤ dist xi xk * dist xi xk -
      (dist xi xj * dist xi xj - dist xi xk * dist xi xk - dist xj xk * dist xj xk) / 2 *
        ((dist xi xj * dist xi xj - dist xi xk * dist xi xk - dist xj xk * dist xj xk) / 2) /
        (dist xj xk * dist xj xk) := by
    have h0 := hq0 (-((dist xi xj * dist xi xj - dist xi xk * dist xi xk -
      dist xj xk * dist xj xk) / 2) / (dist xj xk * dist xj xk))
    have e : hlQ (dist xi xj) (dist xi xk) (dist xj xk)
        (-((dist xi xj * dist xi xj - dist xi xk * dist xi xk -
          dist xj xk * dist xj xk) / 2) / (dist xj xk * dist xj xk)) =
        dist xi xk * dist xi xk -
          (dist xi xj * dist xi xj - dist xi xk * dist xi xk - dist xj xk * dist xj xk) / 2 *
            ((dist xi xj * dist xi xj - dist xi xk * dist xi xk -
              dist xj xk * dist xj xk) / 2) / (dist xj xk * dist xj xk) := by
      unfold hlQ
      field_simp
      ring
    rwa [e] at h0
  have key := hl_stationary_le ha hh h lam
  unfold hlValue hlLam
  rw [max_eq_right hh, max_eq_right (hq0 _), max_eq_right (hq0 _)]
  unfold hlQ
  linarith

/-- **Lower bound**: the closed form does not exceed the objective at any `λ ∈ [0, 1]`.

[formula: cut-redistance] -/
theorem hopfLaxTriangle_le (xi xj xk : E) (Tj Tk : ℝ) {lam : ℝ} (h0 : 0 ≤ lam) (h1 : lam ≤ 1) :
    Ddg.hopfLaxTriangle (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk ≤
      hopfLaxObjective xi xj xk Tj Tk lam := by
  rw [hopfLaxTriangle_eq]
  by_cases h : |Tj - Tk| < dist xj xk
  · rw [ite_eq_left h, hlValue_eq_objective]
    refine clamp_min_of_convex (F := hopfLaxObjective xi xj xk Tj Tk)
      (fun x y θ hθ0 hθ1 => hopfLaxObjective_convex xi xj xk Tj Tk x y hθ0 hθ1)
      (fun x => ?_) h0 h1
    rw [← hlValue_eq_objective, ← hlValue_eq_objective]
    exact hlValue_stationary_le xi xj xk Tj Tk h x
  · rw [ite_eq_right h, hlValue_eq_objective, hlValue_eq_objective]
    exact hopfLaxObjective_endpoint_le xi xj xk Tj Tk (not_lt.mp h) h0 h1

/-- **Attainment**: the closed form is the objective at some `λ ∈ [0, 1]`.

[formula: cut-redistance] -/
theorem hopfLaxTriangle_mem (xi xj xk : E) (Tj Tk : ℝ) :
    ∃ lam : ℝ, 0 ≤ lam ∧ lam ≤ 1 ∧
      Ddg.hopfLaxTriangle (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk =
        hopfLaxObjective xi xj xk Tj Tk lam := by
  rw [hopfLaxTriangle_eq]
  by_cases h : |Tj - Tk| < dist xj xk
  · rw [ite_eq_left h, hlValue_eq_objective]
    exact ⟨_, le_max_left _ _, max_le zero_le_one (min_le_left _ _), rfl⟩
  · rw [ite_eq_right h, hlValue_eq_objective, hlValue_eq_objective]
    rcases min_choice (hopfLaxObjective xi xj xk Tj Tk 0) (hopfLaxObjective xi xj xk Tj Tk 1)
      with h' | h'
    · exact ⟨0, le_rfl, zero_le_one, h'⟩
    · exact ⟨1, zero_le_one, le_rfl, h'⟩

/-- **The Hopf–Lax update in closed form** (SPEC §3.12). For points `x_i, x_j, x_k` of a real
inner product space, the oracle's `hopfLaxTriangle` of the three distances and the two values
`T_j`, `T_k` is the minimum over `λ ∈ [0, 1]` of
`λ T_j + (1 − λ) T_k + ‖x_i − (λ x_j + (1 − λ) x_k)‖`.

[formula: cut-redistance] -/
theorem hopfLaxTriangle_isLeast (xi xj xk : E) (Tj Tk : ℝ) :
    IsLeast (hopfLaxObjective xi xj xk Tj Tk '' Set.Icc 0 1)
      (Ddg.hopfLaxTriangle (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk) := by
  constructor
  · obtain ⟨lam, h0, h1, h⟩ := hopfLaxTriangle_mem xi xj xk Tj Tk
    exact ⟨lam, ⟨h0, h1⟩, h.symm⟩
  · rintro _ ⟨lam, ⟨h0, h1⟩, rfl⟩
    exact hopfLaxTriangle_le xi xj xk Tj Tk h0 h1

/-- At `λ = 1` the objective is the edge term of `j`, at `λ = 0` that of `k`. -/
theorem hopfLaxObjective_one (xi xj xk : E) (Tj Tk : ℝ) :
    hopfLaxObjective xi xj xk Tj Tk 1 = Ddg.hopfLaxEdge Tj (dist xi xj) := by
  unfold hopfLaxObjective
  rw [hopfLaxEdge_eq, dist_eq_norm]
  simp

theorem hopfLaxObjective_zero (xi xj xk : E) (Tj Tk : ℝ) :
    hopfLaxObjective xi xj xk Tj Tk 0 = Ddg.hopfLaxEdge Tk (dist xi xk) := by
  unfold hopfLaxObjective
  rw [hopfLaxEdge_eq, dist_eq_norm]
  simp

/-- **The triangle update never exceeds either edge term.**

[formula: cut-redistance] -/
theorem hopfLaxTriangle_le_edge (xi xj xk : E) (Tj Tk : ℝ) :
    Ddg.hopfLaxTriangle (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk ≤
        Ddg.hopfLaxEdge Tj (dist xi xj) ∧
      Ddg.hopfLaxTriangle (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk ≤
        Ddg.hopfLaxEdge Tk (dist xi xk) := by
  rw [← hopfLaxObjective_one xi xj xk Tj Tk, ← hopfLaxObjective_zero xi xj xk Tj Tk]
  exact ⟨hopfLaxTriangle_le xi xj xk Tj Tk zero_le_one le_rfl,
    hopfLaxTriangle_le xi xj xk Tj Tk le_rfl zero_le_one⟩

/-- **Causality**: the update is at least the smaller of the two values it is computed from. -/
theorem min_le_hopfLaxTriangle (xi xj xk : E) (Tj Tk : ℝ) :
    min Tj Tk ≤ Ddg.hopfLaxTriangle (dist xi xj) (dist xi xk) (dist xj xk) Tj Tk := by
  obtain ⟨lam, h0, h1, h⟩ := hopfLaxTriangle_mem xi xj xk Tj Tk
  rw [h]
  unfold hopfLaxObjective
  have hn := norm_nonneg (xi - (lam • xj + (1 - lam) • xk))
  have hj : min Tj Tk ≤ Tj := min_le_left _ _
  have hk : min Tj Tk ≤ Tk := min_le_right _ _
  nlinarith [mul_nonneg h0 (sub_nonneg.mpr hj), mul_nonneg (sub_nonneg.mpr h1) (sub_nonneg.mpr hk)]

end Geometry

/-! ## From the side lengths alone -/

/-- Three lengths that satisfy the strict triangle inequalities are the side lengths of a
triangle in the plane (`ℂ`). -/
theorem exists_triangle_points {lij lik ljk : ℝ} (h : IsTriangle lij lik ljk) :
    ∃ xi xj xk : ℂ, dist xi xj = lij ∧ dist xi xk = lik ∧ dist xj xk = ljk := by
  have hij := h.pos_a
  have hik := h.pos_b
  have hjk := h.pos_c
  obtain ⟨x, hx⟩ : ∃ x : ℝ, x = (lik ^ 2 + ljk ^ 2 - lij ^ 2) / (2 * ljk) := ⟨_, rfl⟩
  have hx' : 2 * x * ljk = lik ^ 2 + ljk ^ 2 - lij ^ 2 := by
    rw [hx]
    field_simp
  have hH := h.heron_pos
  have hy0 : 0 ≤ lik ^ 2 - x ^ 2 := by
    have e : 4 * ljk ^ 2 * (lik ^ 2 - x ^ 2) =
        (lij + lik + ljk) * (lik - lij + ljk) * (lij - lik + ljk) * (lij - ljk + lik) := by
      linear_combination (-(2 * x * ljk + (lik ^ 2 + ljk ^ 2 - lij ^ 2))) * hx'
    have h1 : 0 ≤ 4 * ljk ^ 2 * (lik ^ 2 - x ^ 2) := by
      rw [e]
      exact hH.le
    exact nonneg_of_mul_nonneg_right h1 (by positivity)
  refine ⟨(x : ℂ) + (√(lik ^ 2 - x ^ 2) : ℝ) * Complex.I, (ljk : ℂ), 0, ?_, ?_, ?_⟩
  · have e : (x : ℂ) + (√(lik ^ 2 - x ^ 2) : ℝ) * Complex.I - (ljk : ℂ) =
        ((x - ljk : ℝ) : ℂ) + (√(lik ^ 2 - x ^ 2) : ℝ) * Complex.I := by
      push_cast
      ring
    rw [Complex.dist_eq, e, Complex.norm_add_mul_I, Real.sq_sqrt hy0,
      show (x - ljk) ^ 2 + (lik ^ 2 - x ^ 2) = lij ^ 2 by linear_combination (-1 : ℝ) * hx',
      Real.sqrt_sq hij.le]
  · rw [Complex.dist_eq, sub_zero, Complex.norm_add_mul_I, Real.sq_sqrt hy0,
      show x ^ 2 + (lik ^ 2 - x ^ 2) = lik ^ 2 by ring, Real.sqrt_sq hik.le]
  · rw [Complex.dist_eq, sub_zero, Complex.norm_real, Real.norm_of_nonneg hjk.le]

/-- **The Hopf–Lax update from the side lengths.** If `ℓ_ij, ℓ_ik, ℓ_jk` satisfy the strict
triangle inequalities, the oracle's `hopfLaxTriangle` is the minimum over `λ ∈ [0, 1]` of
`λ T_j + (1 − λ) T_k + (distance from the vertex i to the point λ j + (1 − λ) k)`, the distance
being the oracle's `baryLength` in the face `i j k` (barycentric difference `(1, −λ, −(1 − λ))`).

[formula: cut-redistance] -/
theorem hopfLaxTriangle_isLeast_intrinsic {lij lik ljk : ℝ} (h : IsTriangle lij lik ljk)
    (Tj Tk : ℝ) :
    IsLeast ((fun lam : ℝ => lam * Tj + (1 - lam) * Tk +
        Ddg.baryLength lij ljk lik 1 (-lam) (-(1 - lam))) '' Set.Icc 0 1)
      (Ddg.hopfLaxTriangle lij lik ljk Tj Tk) := by
  obtain ⟨xi, xj, xk, hij, hik, hjk⟩ := exists_triangle_points h
  have key := hopfLaxTriangle_isLeast xi xj xk Tj Tk
  rw [hij, hik, hjk] at key
  have e : (fun lam : ℝ => lam * Tj + (1 - lam) * Tk +
      Ddg.baryLength lij ljk lik 1 (-lam) (-(1 - lam))) = hopfLaxObjective xi xj xk Tj Tk := by
    funext lam
    unfold hopfLaxObjective
    have hb := baryLength_eq_norm xi xj xk (da := 1) (db := -lam) (dc := -(1 - lam)) (by ring)
    rw [← dist_eq_norm, ← dist_eq_norm, ← dist_eq_norm, hij, hjk, dist_comm xk xi, hik] at hb
    have e2 : (1 : ℝ) • xi + (-lam) • xj + (-(1 - lam)) • xk =
        xi - (lam • xj + (1 - lam) • xk) := by module
    rw [hb, e2]
  rw [e]
  exact key

end Cuts

end DdgProofs
