import Mathlib

/-!
# Applications of intrinsic triangulations (Sec. 5.2, 5.3): one triangle's formulas

`crates/ddg-intrinsic-apps` (`SPEC.md`, decisions N30–N39); the paper Sharp, Soliman, Crane 2019,
Sec. 5.2 (finite elements) and 5.3 (the heat method of Crane, Weischedel, Wardetzky 2013).

A triangle is given by the coordinates of its corners `p₀, p₁, p₂` in the plane of its layout;
`D = (p₁ − p₀) × (p₂ − p₀)` is twice its signed area (positive counter-clockwise).

* **`app-heat-gradient`** — `Σ_k u_k J(p_{k+2} − p_{k+1}) = D ∇u` for the linear function
  `u(x, y) = a x + b y + c` sampled at the corners (`gradSum_linear`), `J` the quarter turn. So
  the direction of `gradSum` is that of `∇u` for `D > 0`, and the unit field `−∇u/|∇u|` of the heat
  method's step II is `−gradSum/|gradSum|` (`gradSum_normalized`).
* **`app-heat-divergence`** — the integrated divergence of a constant field `X` at corner `k`,
  `½ (cot θ_{k+2} (p_{k+1} − p_k)·X + cot θ_{k+1} (p_{k+2} − p_k)·X)` (Crane et al. 2013, Eq. 8,
  the cotangents as dot over cross product), is `−½ J(p_{k+2} − p_{k+1})·X`, which is
  `−A ∇φ_k · X` (`divergence_corner`); over the three corners it sums to zero
  (`divergence_sum_zero`), so the Poisson problem of step III is solvable on every component.
* **`app-fem-galerkin-mass`** — the Galerkin mass of one triangle, `A/6` on the diagonal and
  `A/12` off it: its rows sum to the lumped mass `A/3` (`galerkin_row_sum`), and its quadratic
  form is `(A/12)((Σu)² + Σu²)` (`galerkin_form_eq`), so it lies between `A/12 Σu²` and
  `A/3 Σu²` (`galerkin_form_bounds`): positive definite, and within a factor four of the lumped
  mass.
* **`app-fem-midpoint-subdivision`** — the regular 4-1 subdivision from lengths alone: the half
  of a side has half its length and the segment between the midpoints of two sides has half the
  length of the third side (`midpoint_half`, `midline`), in any real normed space.
-/

noncomputable section

namespace DDG

namespace IntrinsicApps

/-- `Σ_k u_k J(p_{k+2} − p_{k+1})` with `J (x, y) = (−y, x)`: twice the area times the gradient of
the linear interpolant of `u₀, u₁, u₂` (the factor is left out; it does not change the
direction). As in `heat.rs` (`unit_descent`).

[formula: app-heat-gradient] -/
def gradSum (x₀ y₀ x₁ y₁ x₂ y₂ u₀ u₁ u₂ : ℝ) : ℝ × ℝ :=
  (u₀ * -(y₂ - y₁) + u₁ * -(y₀ - y₂) + u₂ * -(y₁ - y₀),
    u₀ * (x₂ - x₁) + u₁ * (x₀ - x₂) + u₂ * (x₁ - x₀))

/-- Twice the signed area of the triangle. -/
def twiceArea (x₀ y₀ x₁ y₁ x₂ y₂ : ℝ) : ℝ :=
  (x₁ - x₀) * (y₂ - y₀) - (y₁ - y₀) * (x₂ - x₀)

/-- **The gradient formula is exact on linear functions**: sampling `u = a x + b y + c` at the
corners gives `gradSum = D · (a, b)`.

[formula: app-heat-gradient] -/
theorem gradSum_linear (x₀ y₀ x₁ y₁ x₂ y₂ a b c : ℝ) :
    gradSum x₀ y₀ x₁ y₁ x₂ y₂ (a * x₀ + b * y₀ + c) (a * x₁ + b * y₁ + c)
        (a * x₂ + b * y₂ + c) =
      (twiceArea x₀ y₀ x₁ y₁ x₂ y₂ * a, twiceArea x₀ y₀ x₁ y₁ x₂ y₂ * b) := by
  unfold gradSum twiceArea
  ext <;> ring

/-- **The normalized field does not depend on the factor**: for a counter-clockwise triangle
(`D > 0`) the unit vector of `gradSum` is that of the gradient `(a, b)`.

[formula: app-heat-gradient] -/
theorem gradSum_normalized (x₀ y₀ x₁ y₁ x₂ y₂ a b c : ℝ)
    (hD : 0 < twiceArea x₀ y₀ x₁ y₁ x₂ y₂) :
    let g := gradSum x₀ y₀ x₁ y₁ x₂ y₂ (a * x₀ + b * y₀ + c) (a * x₁ + b * y₁ + c)
        (a * x₂ + b * y₂ + c)
    (‖(g.1, g.2)‖ : ℝ) = twiceArea x₀ y₀ x₁ y₁ x₂ y₂ * ‖((a, b) : ℝ × ℝ)‖ := by
  intro g
  have h := gradSum_linear x₀ y₀ x₁ y₁ x₂ y₂ a b c
  simp only [g, h]
  rw [show ((twiceArea x₀ y₀ x₁ y₁ x₂ y₂ * a, twiceArea x₀ y₀ x₁ y₁ x₂ y₂ * b) : ℝ × ℝ) =
      twiceArea x₀ y₀ x₁ y₁ x₂ y₂ • ((a, b) : ℝ × ℝ) by ext <;> simp [smul_eq_mul]]
  rw [norm_smul, Real.norm_of_nonneg hD.le]

/-- The integrated divergence of the constant field `(X, Y)` at the corner `p₀` of the triangle
`p₀ p₁ p₂`, the cotangents written as dot product over twice the signed area:
`½ (cot θ₂ (p₁ − p₀)·X + cot θ₁ (p₂ − p₀)·X)` with `cot θ₂ = (p₀ − p₂)·(p₁ − p₂) / D`,
`cot θ₁ = (p₀ − p₁)·(p₂ − p₁) / D`.

[formula: app-heat-divergence] -/
def divergenceCorner (x₀ y₀ x₁ y₁ x₂ y₂ X Y : ℝ) : ℝ :=
  let D := twiceArea x₀ y₀ x₁ y₁ x₂ y₂
  let cot₂ := ((x₀ - x₂) * (x₁ - x₂) + (y₀ - y₂) * (y₁ - y₂)) / D
  let cot₁ := ((x₀ - x₁) * (x₂ - x₁) + (y₀ - y₁) * (y₂ - y₁)) / D
  (cot₂ * ((x₁ - x₀) * X + (y₁ - y₀) * Y) + cot₁ * ((x₂ - x₀) * X + (y₂ - y₀) * Y)) / 2

/-- **Crane et al.'s divergence is the weak divergence** `−A ∇φ₀ · X = −½ J(p₂ − p₁) · X`.

[formula: app-heat-divergence] -/
theorem divergence_corner (x₀ y₀ x₁ y₁ x₂ y₂ X Y : ℝ) (hD : twiceArea x₀ y₀ x₁ y₁ x₂ y₂ ≠ 0) :
    divergenceCorner x₀ y₀ x₁ y₁ x₂ y₂ X Y = -((-(y₂ - y₁)) * X + (x₂ - x₁) * Y) / 2 := by
  unfold divergenceCorner
  simp only
  unfold twiceArea at hD ⊢
  field_simp
  ring

/-- **The divergence sums to zero over the corners** (the three cyclic orders of the corners),
so `L φ = −div X` is solvable.

[formula: app-heat-divergence] -/
theorem divergence_sum_zero (x₀ y₀ x₁ y₁ x₂ y₂ X Y : ℝ)
    (hD : twiceArea x₀ y₀ x₁ y₁ x₂ y₂ ≠ 0) :
    divergenceCorner x₀ y₀ x₁ y₁ x₂ y₂ X Y + divergenceCorner x₁ y₁ x₂ y₂ x₀ y₀ X Y +
        divergenceCorner x₂ y₂ x₀ y₀ x₁ y₁ X Y = 0 := by
  have h₁ : twiceArea x₁ y₁ x₂ y₂ x₀ y₀ ≠ 0 := by
    unfold twiceArea at hD ⊢; intro h; apply hD; linarith
  have h₂ : twiceArea x₂ y₂ x₀ y₀ x₁ y₁ ≠ 0 := by
    unfold twiceArea at hD ⊢; intro h; apply hD; linarith
  rw [divergence_corner _ _ _ _ _ _ _ _ hD, divergence_corner _ _ _ _ _ _ _ _ h₁,
    divergence_corner _ _ _ _ _ _ _ _ h₂]
  ring

/-- The quadratic form of the Galerkin mass matrix of one triangle of area `A`: `A/6` for a
corner with itself, `A/12` for two different corners.

[formula: app-fem-galerkin-mass] -/
def galerkinForm (A u₀ u₁ u₂ : ℝ) : ℝ :=
  A / 6 * (u₀ ^ 2 + u₁ ^ 2 + u₂ ^ 2) + A / 12 * (2 * u₀ * u₁ + 2 * u₁ * u₂ + 2 * u₂ * u₀)

/-- **Every row sums to the lumped mass** `A/3`.

[formula: app-fem-galerkin-mass] -/
theorem galerkin_row_sum (A : ℝ) : A / 6 + A / 12 + A / 12 = A / 3 := by ring

/-- [formula: app-fem-galerkin-mass] -/
theorem galerkin_form_eq (A u₀ u₁ u₂ : ℝ) :
    galerkinForm A u₀ u₁ u₂ = A / 12 * ((u₀ + u₁ + u₂) ^ 2 + (u₀ ^ 2 + u₁ ^ 2 + u₂ ^ 2)) := by
  unfold galerkinForm; ring

/-- **Between a quarter of the lumped mass and the lumped mass** (`A > 0`): positive definite,
and the eigenvalues of `L x = λ M x` with the Galerkin mass lie between those with the lumped
mass and four times them (the bound used in `tests/fem.rs`).

[formula: app-fem-galerkin-mass] -/
theorem galerkin_form_bounds (A u₀ u₁ u₂ : ℝ) (hA : 0 < A) :
    A / 12 * (u₀ ^ 2 + u₁ ^ 2 + u₂ ^ 2) ≤ galerkinForm A u₀ u₁ u₂ ∧
      galerkinForm A u₀ u₁ u₂ ≤ A / 3 * (u₀ ^ 2 + u₁ ^ 2 + u₂ ^ 2) := by
  rw [galerkin_form_eq]
  constructor
  · nlinarith [sq_nonneg (u₀ + u₁ + u₂)]
  · nlinarith [sq_nonneg (u₀ - u₁), sq_nonneg (u₁ - u₂), sq_nonneg (u₂ - u₀)]

/-- **The half of a side** has half its length.

[formula: app-fem-midpoint-subdivision] -/
theorem midpoint_half {E : Type*} [SeminormedAddCommGroup E] [NormedSpace ℝ E] (a b : E) :
    ‖a - (1 / 2 : ℝ) • (a + b)‖ = ‖a - b‖ / 2 := by
  have : a - (1 / 2 : ℝ) • (a + b) = (1 / 2 : ℝ) • (a - b) := by module
  rw [this, norm_smul]
  norm_num
  ring

/-- **The midline**: the segment between the midpoints of the sides `ab` and `ac` has half the
length of `bc`.

[formula: app-fem-midpoint-subdivision] -/
theorem midline {E : Type*} [SeminormedAddCommGroup E] [NormedSpace ℝ E] (a b c : E) :
    ‖(1 / 2 : ℝ) • (a + b) - (1 / 2 : ℝ) • (a + c)‖ = ‖b - c‖ / 2 := by
  rw [← smul_sub, show a + b - (a + c) = b - c by abel, norm_smul]
  norm_num
  ring

end IntrinsicApps

end DDG
