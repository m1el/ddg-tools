import DdgProofs.TriangleIntegral.Basic
import DdgProofs.Cotan
import DdgProofs.Cuts.Weighted
import DdgProofs.IntrinsicApps.HeatFem
import DdgProofs.Intrinsic.Move
import DdgProofs.Laplace.Cotan

/-!
# Finite element integrals over a triangle

The integrals behind the finite element matrices of the library, as Lebesgue integrals over a
triangle `T` of the plane (`DdgProofs/TriangleIntegral/Basic.lean`), with `λᵢ = bary p i` the hat
functions of the corners and `A = area p` the Heron area of the Euclidean side lengths.

* **`app-fem-galerkin-mass`** — `∫_T λᵢ λⱼ = A/6` (`i = j`) and `A/12` (`i ≠ j`)
  (`integral_bary_mul_bary`); `∫_T λᵢ = A/3`, the lumped mass, is the row sum
  (`sum_integral_bary_mul_bary`); `∫_T (Σ uᵢ λᵢ)² = galerkinForm A u₀ u₁ u₂`
  (`integral_interp_sq`).
* **`cotan-fem-stiffness`** — the gradient of `λᵢ` is the constant `J(p_{i+2} − p_{i+1})/D`
  (`grad_bary`); `∫_T ∇f·∇g` of two linear interpolants is the cotan form
  `½ Σ cot θ_k (f_i − f_j)(g_i − g_j)` of the side lengths (`integral_grad_dot_grad`), so
  `∫_T ∇λᵢ·∇λⱼ = −½ cot θ_k` for `i ≠ j` (`integral_grad_bary_dot`) and
  `∫_T |∇f|² = cotEnergy` (`integral_grad_sq`); summed over the faces of a mesh, each laid out
  isometrically, it is `uᵀ L u` (`TriMesh.sum_integral_grad_sq`).
* **`cut-weighted-fem-stiffness`** — `∫_T a ∇f·∇g = ā · (cotan form)` for linear `a`, `ā` the mean
  of its corner values (`integral_weighted_grad_dot_grad`); over a mesh, `uᵀ L_ā u`
  (`TriMesh.sum_integral_weighted_grad_sq`).
* **`fem-area-coordinate-rule`** — the area-coordinate rule and the slice formula of
  `Basic.lean` with the area: `integral_bary_monomial_area`, and the slice formula for
  `G(s) = e^{iΩs}` (`integral_bary_slice_exp`), which is what `monomial_integral` of
  `crates/ddg-intrinsic-apps/src/vector/knoeppel.rs` evaluates.
* **`bff-conformal-energy-triangle`** — on a counter-clockwise triangle,
  `2∫_T |∂̄f|² = ½(∫|∇a|² + ∫|∇b|²) − ½ sda(f(p₀), f(p₁), f(p₂))` for `f = a + ib` linear
  (`integral_dbarSq`), zero iff `f` satisfies the Cauchy–Riemann equations
  (`integral_dbarSq_eq_zero_iff`).
* **`it-odt-interpolation-error`** — Chen–Holst Lemma 4.1 on one triangle, for `u = |x|²`:
  `∫_T |u_I − u| = (A/3) Σ|pᵢ|² − ∫_T |x|² = (A/12)(ℓ₀₁² + ℓ₁₂² + ℓ₂₀²)`
  (`integral_abs_odt_error`, `integral_odt_error_eq_sides`), and the star energy `odtEnergy` of
  `Intrinsic/Move.lean` is `Σ_k ∫_{T_k} u_I` (`odtEnergy_eq_sum_integral`).
* **`app-amr-element-norms`** — for linear `R`: `‖R‖²_T = (A/6)(Σ rᵢ² + Σ_{i<j} rᵢ rⱼ)` and
  `‖R − R̄‖²_T = ‖R‖²_T − A R̄²`, `R̄` the corner mean (`integral_amr_sq`, `integral_amr_osc`), as
  `amr::estimate` evaluates them.
-/

noncomputable section

namespace DDG

namespace TriangleIntegral

open MeasureTheory Set Real
open scoped Nat

/-! ## Helpers -/

section Helpers

variable (p : Fin 3 → ℝ × ℝ)

theorem continuous_triMap : Continuous (triMap p) := by
  unfold triMap; fun_prop

theorem isCompact_triangle : IsCompact (triangle p) :=
  isCompact_stdTri.image (continuous_triMap p)

theorem integrableOn_triangle {E : Type*} [NormedAddCommGroup E] {f : ℝ × ℝ → E}
    (hf : Continuous f) : IntegrableOn f (triangle p) :=
  hf.continuousOn.integrableOn_compact (isCompact_triangle p)

theorem continuous_bary (i : Fin 3) : Continuous (bary p i) := by
  have h0 : bary p 0 = fun x => sda x (p 1) (p 2) / D p := funext (bary_zero p)
  have h1 : bary p 1 = fun x => sda (p 0) x (p 2) / D p := funext (bary_one p)
  have h2 : bary p 2 = fun x => sda (p 0) (p 1) x / D p := funext (bary_two p)
  fin_cases i
  · show Continuous (bary p 0); rw [h0]; unfold sda; fun_prop
  · show Continuous (bary p 1); rw [h1]; unfold sda; fun_prop
  · show Continuous (bary p 2); rw [h2]; unfold sda; fun_prop

/-- The Heron area of the Euclidean side lengths `ℓ₁₂, ℓ₂₀, ℓ₀₁`. -/
def area : ℝ := triArea (dist2 (p 1) (p 2)) (dist2 (p 2) (p 0)) (dist2 (p 0) (p 1))

theorem abs_D_eq_two_mul_area : |D p| = 2 * area p := abs_sda_eq_two_mul_triArea _ _ _

/-- The area is the Lebesgue measure of the triangle. -/
theorem volume_real_triangle_eq_area {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) :
    volume.real (triangle p) = area p :=
  volume_real_triangle_eq_triArea hD

/-- The linear interpolant `Σ uᵢ λᵢ` of the corner values `u`. -/
def interp (u : Fin 3 → ℝ) (x : ℝ × ℝ) : ℝ := ∑ i, u i * bary p i x

theorem continuous_interp (u : Fin 3 → ℝ) : Continuous (interp p u) := by
  unfold interp
  exact continuous_finsetSum _ fun i _ => continuous_const.mul (continuous_bary p i)

/-- The cotangent of the corner `0` from the side lengths is `(p₁ − p₀)·(p₂ − p₀)/|D|`. -/
theorem cornerCot_eq_dot_div (p : Fin 3 → ℝ × ℝ) :
    cornerCot (dist2 (p 1) (p 2)) (dist2 (p 2) (p 0)) (dist2 (p 0) (p 1)) =
      dot2 (p 1 - p 0) (p 2 - p 0) / |D p| := by
  unfold cornerCot
  rw [dist2_sq, dist2_sq, dist2_sq, show 4 * triArea (dist2 (p 1) (p 2)) (dist2 (p 2) (p 0))
    (dist2 (p 0) (p 1)) = 2 * |D p| by rw [abs_D_eq_two_mul_area, area]; ring]
  rw [show sqDist2 (p 2) (p 0) + sqDist2 (p 0) (p 1) - sqDist2 (p 1) (p 2) =
      2 * dot2 (p 1 - p 0) (p 2 - p 0) by
    simp only [sqDist2, dot2, Prod.fst_sub, Prod.snd_sub]; ring]
  rw [mul_div_mul_left _ _ two_ne_zero]

end Helpers

variable {p : Fin 3 → ℝ × ℝ}

/-- An integrand that is a monomial in the barycentric coordinates. -/
theorem integral_eq_of_monomial (hD : D p ≠ 0) {f : ℝ × ℝ → ℝ} (a b c : ℕ)
    (h : ∀ x, f x = bary p 0 x ^ a * bary p 1 x ^ b * bary p 2 x ^ c) :
    ∫ x in triangle p, f x = 2 * area p * ((a ! * b ! * c ! : ℕ) / ((a + b + c + 2)! : ℕ) : ℝ) := by
  rw [show f = _ from funext h, integral_bary_monomial hD, abs_D_eq_two_mul_area]

/-! ## `app-fem-galerkin-mass`: the mass matrix -/

/-- **The lumped mass is the integral of the hat function**: `∫_T λᵢ = A/3`.

[formula: app-fem-galerkin-mass] -/
theorem integral_bary (hD : D p ≠ 0) (i : Fin 3) : ∫ x in triangle p, bary p i x = area p / 3 := by
  fin_cases i
  · rw [integral_eq_of_monomial hD 1 0 0 (fun x => by simp)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 0 1 0 (fun x => by simp)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 0 0 1 (fun x => by simp)]; norm_num [Nat.factorial]; ring

/-- **The Galerkin mass matrix is the matrix of the integrals `∫_T λᵢ λⱼ`**: `A/6` on the
diagonal and `A/12` off it, as `fem::galerkin_mass_matrix` assembles it.

[formula: app-fem-galerkin-mass] -/
theorem integral_bary_mul_bary (hD : D p ≠ 0) (i j : Fin 3) :
    ∫ x in triangle p, bary p i x * bary p j x = if i = j then area p / 6 else area p / 12 := by
  fin_cases i <;> fin_cases j <;> simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk, Fin.isValue]
  · rw [integral_eq_of_monomial hD 2 0 0 (fun x => by ring)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 1 1 0 (fun x => by ring)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 1 0 1 (fun x => by ring)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 1 1 0 (fun x => by ring)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 0 2 0 (fun x => by ring)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 0 1 1 (fun x => by ring)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 1 0 1 (fun x => by ring)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 0 1 1 (fun x => by ring)]; norm_num [Nat.factorial]; ring
  · rw [integral_eq_of_monomial hD 0 0 2 (fun x => by ring)]; norm_num [Nat.factorial]; ring

/-- **The rows of the Galerkin mass sum to the lumped mass**: `Σⱼ ∫_T λᵢ λⱼ = ∫_T λᵢ`.

[formula: app-fem-galerkin-mass] -/
theorem sum_integral_bary_mul_bary (hD : D p ≠ 0) (i : Fin 3) :
    ∑ j, ∫ x in triangle p, bary p i x * bary p j x = ∫ x in triangle p, bary p i x := by
  rw [integral_bary hD, Fin.sum_univ_three, integral_bary_mul_bary hD, integral_bary_mul_bary hD,
    integral_bary_mul_bary hD]
  fin_cases i <;> simp <;> ring

/-- **The Galerkin quadratic form is the `L²` norm of the interpolant**:
`∫_T (Σ uᵢ λᵢ)² = (A/6) Σ uᵢ² + (A/12) Σ_{i≠j} uᵢ uⱼ`.

[formula: app-fem-galerkin-mass] -/
theorem integral_interp_sq (hD : D p ≠ 0) (u : Fin 3 → ℝ) :
    ∫ x in triangle p, interp p u x ^ 2 = IntrinsicApps.galerkinForm (area p) (u 0) (u 1) (u 2) := by
  have e : ∀ x, interp p u x ^ 2 = ∑ i, ∑ j, u i * u j * (bary p i x * bary p j x) := by
    intro x
    rw [interp, sq, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring
  simp_rw [e]
  have hij : ∀ i j, Integrable (fun x => u i * u j * (bary p i x * bary p j x))
      (volume.restrict (triangle p)) := fun i j =>
    integrableOn_triangle p (continuous_const.mul ((continuous_bary p i).mul (continuous_bary p j)))
  rw [integral_finsetSum _ fun i _ => integrable_finsetSum _ fun j _ => hij i j]
  simp_rw [integral_finsetSum _ fun j _ => hij _ j, integral_const_mul, integral_bary_mul_bary hD]
  simp only [Fin.sum_univ_three, Fin.isValue, IntrinsicApps.galerkinForm]
  simp
  ring

/-- **The integral of a linear interpolant is the area times the mean of its corner values.** -/
theorem integral_interp (hD : D p ≠ 0) (u : Fin 3 → ℝ) :
    ∫ x in triangle p, interp p u x = area p * ((u 0 + u 1 + u 2) / 3) := by
  unfold interp
  rw [integral_finsetSum]
  · simp_rw [integral_const_mul, integral_bary hD, Fin.sum_univ_three]
    ring
  · exact fun i _ => integrableOn_triangle p (continuous_const.mul (continuous_bary p i))

/-! ## `cotan-fem-stiffness`: gradients and the stiffness matrix -/

/-- The gradient `(∂₁f, ∂₂f)` of a function on the plane, read off its Fréchet derivative. -/
def grad (f : ℝ × ℝ → ℝ) (x : ℝ × ℝ) : ℝ × ℝ := (fderiv ℝ f x (1, 0), fderiv ℝ f x (0, 1))

/-- The side opposite the corner `i`, `p_{i+2} − p_{i+1}`. -/
def side (p : Fin 3 → ℝ × ℝ) (i : Fin 3) : ℝ × ℝ := p (i + 2) - p (i + 1)

/-- The length of the side opposite the corner `i`, `|p_{i+1} p_{i+2}|`. -/
def sideLen (p : Fin 3 → ℝ × ℝ) (i : Fin 3) : ℝ := dist2 (p (i + 1)) (p (i + 2))

/-- The cotangent of the corner `k` from the three side lengths (`cornerCot`, the corner opposite
its first argument). -/
def cotCorner (p : Fin 3 → ℝ × ℝ) (k : Fin 3) : ℝ :=
  cornerCot (sideLen p k) (sideLen p (k + 1)) (sideLen p (k + 2))

section SideLemmas

variable (p)

theorem side_zero : side p 0 = p 2 - p 1 := rfl
theorem side_one : side p 1 = p 0 - p 2 := rfl
theorem side_two : side p 2 = p 1 - p 0 := rfl
theorem sideLen_zero : sideLen p 0 = dist2 (p 1) (p 2) := rfl
theorem sideLen_one : sideLen p 1 = dist2 (p 2) (p 0) := rfl
theorem sideLen_two : sideLen p 2 = dist2 (p 0) (p 1) := rfl

theorem cotCorner_zero : cotCorner p 0 = dot2 (p 1 - p 0) (p 2 - p 0) / |D p| :=
  cornerCot_eq_dot_div p

theorem cotCorner_one : cotCorner p 1 = dot2 (p 2 - p 1) (p 0 - p 1) / |D p| := by
  have h := cornerCot_eq_dot_div (rot p)
  simp only [D_rot, rot_zero, rot_one, rot_two] at h
  exact h

theorem cotCorner_two : cotCorner p 2 = dot2 (p 0 - p 2) (p 1 - p 2) / |D p| := by
  have h := cornerCot_eq_dot_div (rot (rot p))
  simp only [D_rot, rot_zero, rot_one, rot_two] at h
  exact h

end SideLemmas

/-- The gradient of the hat function `λᵢ`: `J(p_{i+2} − p_{i+1})/D` with `J(x, y) = (−y, x)`
(`∇b_m = J d[m]/(2A)` in `knoeppel.rs`). -/
def baryGrad (p : Fin 3 → ℝ × ℝ) (i : Fin 3) : ℝ × ℝ := (-(side p i).2 / D p, (side p i).1 / D p)

/-- The hat function is affine with linear part `baryGrad`. -/
theorem bary_affine (p : Fin 3 → ℝ × ℝ) (i : Fin 3) (x : ℝ × ℝ) :
    bary p i x = bary p i 0 + (baryGrad p i).1 * x.1 + (baryGrad p i).2 * x.2 := by
  fin_cases i <;>
    simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk, Fin.isValue, bary_zero, bary_one,
      bary_two, baryGrad, side_zero, side_one, side_two] <;>
    simp only [sda, div_eq_mul_inv, Prod.fst_zero, Prod.snd_zero, Prod.fst_sub, Prod.snd_sub] <;>
    ring

theorem hasFDerivAt_bary (p : Fin 3 → ℝ × ℝ) (i : Fin 3) (x : ℝ × ℝ) :
    HasFDerivAt (bary p i) (Intrinsic.toDual (baryGrad p i)) x := by
  have h : bary p i = fun y => bary p i 0 + Intrinsic.toDual (baryGrad p i) y := by
    funext y; rw [Intrinsic.toDual_apply, bary_affine p i y]; ring
  rw [h]
  exact (Intrinsic.toDual _).hasFDerivAt.const_add _

/-- **The gradient of the hat function is the constant `baryGrad`.**

[formula: cotan-fem-stiffness] -/
theorem grad_bary (p : Fin 3 → ℝ × ℝ) (i : Fin 3) (x : ℝ × ℝ) : grad (bary p i) x = baryGrad p i := by
  simp [grad, (hasFDerivAt_bary p i x).fderiv, Intrinsic.toDual_apply]

/-- The gradient of the interpolant `Σ uᵢ λᵢ`. -/
def gradInterp (p : Fin 3 → ℝ × ℝ) (u : Fin 3 → ℝ) : ℝ × ℝ := ∑ i, u i • baryGrad p i

theorem hasFDerivAt_interp (p : Fin 3 → ℝ × ℝ) (u : Fin 3 → ℝ) (x : ℝ × ℝ) :
    HasFDerivAt (interp p u) (Intrinsic.toDual (gradInterp p u)) x := by
  have h : interp p u = fun y => interp p u 0 + Intrinsic.toDual (gradInterp p u) y := by
    funext y
    rw [Intrinsic.toDual_apply]
    simp only [interp, gradInterp, Fin.sum_univ_three, Prod.fst_add, Prod.snd_add, Prod.smul_fst,
      Prod.smul_snd, smul_eq_mul]
    rw [bary_affine p 0 y, bary_affine p 1 y, bary_affine p 2 y]
    ring
  rw [h]
  exact (Intrinsic.toDual _).hasFDerivAt.const_add _

/-- **The gradient of a linear interpolant is the constant `Σ uᵢ ∇λᵢ`.** -/
theorem grad_interp (p : Fin 3 → ℝ × ℝ) (u : Fin 3 → ℝ) (x : ℝ × ℝ) :
    grad (interp p u) x = gradInterp p u := by
  simp [grad, (hasFDerivAt_interp p u x).fderiv, Intrinsic.toDual_apply]

/-- The cotan bilinear form of a triangle with sides `a b c` (`a` opposite the corner `0`):
`½ (cot θ₀ (u₁ − u₂)(w₁ − w₂) + cot θ₁ (u₂ − u₀)(w₂ − w₀) + cot θ₂ (u₀ − u₁)(w₀ − w₁))`. -/
def cotForm (a b c : ℝ) (u w : Fin 3 → ℝ) : ℝ :=
  (cornerCot a b c * (u 1 - u 2) * (w 1 - w 2) + cornerCot b c a * (u 2 - u 0) * (w 2 - w 0) +
    cornerCot c a b * (u 0 - u 1) * (w 0 - w 1)) / 2

/-- On the diagonal the cotan form is the cotan energy. -/
theorem cotForm_self (a b c : ℝ) (u : Fin 3 → ℝ) :
    cotForm a b c u u = cotEnergy a b c (u 0) (u 1) (u 2) := by
  unfold cotForm cotEnergy; ring

/-- The constant integrand: area times the dot product of the gradients is the cotan form. -/
theorem area_mul_dot_gradInterp (hD : D p ≠ 0) (u w : Fin 3 → ℝ) :
    area p * dot2 (gradInterp p u) (gradInterp p w) =
      cotForm (sideLen p 0) (sideLen p 1) (sideLen p 2) u w := by
  have hA : area p = |D p| / 2 := by rw [abs_D_eq_two_mul_area]; ring
  have c0 := cotCorner_zero p
  have c1 := cotCorner_one p
  have c2 := cotCorner_two p
  simp only [cotCorner, Fin.isValue, Fin.reduceAdd, zero_add] at c0 c1 c2
  rw [cotForm, c0, c1, c2, hA]
  simp only [gradInterp, baryGrad, Fin.sum_univ_three, side_zero, side_one, side_two, dot2,
    Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul, Prod.fst_sub,
    Prod.snd_sub]
  set d := D p with hd
  have hd' : -d ≠ 0 := neg_ne_zero.mpr hD
  rcases abs_choice d with h | h <;> rw [h] <;> field_simp <;> ring

/-- **The stiffness matrix is the matrix of the integrals `∫_T ∇f·∇g`**: for the linear
interpolants of `u` and `w`, `∫_T ∇f·∇g = ½ Σ_k cot θ_k (u_i − u_j)(w_i − w_j)`, the cotangents
from the three side lengths.

[formula: cotan-fem-stiffness] -/
theorem integral_grad_dot_grad (hD : D p ≠ 0) (u w : Fin 3 → ℝ) :
    ∫ x in triangle p, dot2 (grad (interp p u) x) (grad (interp p w) x) =
      cotForm (sideLen p 0) (sideLen p 1) (sideLen p 2) u w := by
  simp_rw [grad_interp]
  rw [integral_triangle_const hD, smul_eq_mul, ← area_mul_dot_gradInterp hD,
    abs_D_eq_two_mul_area]
  ring

/-- **The Dirichlet energy of a linear interpolant is the cotan energy**
`∫_T |∇f|² = ½ Σ cot θ_k (u_i − u_j)²` (`[formula: cotan-formula]` as a Lebesgue integral).

[formula: cotan-fem-stiffness] -/
theorem integral_grad_sq (hD : D p ≠ 0) (u : Fin 3 → ℝ) :
    ∫ x in triangle p, dot2 (grad (interp p u) x) (grad (interp p u) x) =
      cotEnergy (sideLen p 0) (sideLen p 1) (sideLen p 2) (u 0) (u 1) (u 2) := by
  rw [integral_grad_dot_grad hD, cotForm_self]

theorem interp_single (p : Fin 3 → ℝ × ℝ) (i : Fin 3) : interp p (Pi.single i 1) = bary p i := by
  funext x
  simp [interp, Pi.single_apply]

/-- **The off-diagonal stiffness entries**: `∫_T ∇λᵢ·∇λ_{i+1} = −½ cot θ_{i+2}`, the corner
opposite the edge `{i, i+1}`.

[formula: cotan-fem-stiffness] -/
theorem integral_grad_bary_dot_succ (hD : D p ≠ 0) (i : Fin 3) :
    ∫ x in triangle p, dot2 (grad (bary p i) x) (grad (bary p (i + 1)) x) = -cotCorner p (i + 2) / 2 := by
  rw [← interp_single, ← interp_single, integral_grad_dot_grad hD]
  fin_cases i <;> simp [cotForm, cotCorner]

/-- **The diagonal stiffness entries**: `∫_T |∇λᵢ|² = ½ (cot θ_{i+1} + cot θ_{i+2})`.

[formula: cotan-fem-stiffness] -/
theorem integral_grad_bary_sq (hD : D p ≠ 0) (i : Fin 3) :
    ∫ x in triangle p, dot2 (grad (bary p i) x) (grad (bary p i) x) =
      (cotCorner p (i + 1) + cotCorner p (i + 2)) / 2 := by
  rw [← interp_single, integral_grad_dot_grad hD]
  fin_cases i
  · simp [cotForm, cotCorner]
  · simp [cotForm, cotCorner]; ring
  · simp [cotForm, cotCorner]

/-! ## `cut-weighted-fem-stiffness`: a linear coefficient -/

/-- **The weighted stiffness of one triangle**: for linear `a`, `f`, `g`,
`∫_T a ∇f·∇g = ā · (cotan form)`, `ā = (a₀ + a₁ + a₂)/3` the mean of the corner values of `a`
(the gradient is constant, and a linear function integrates to its mean times the area).

[formula: cut-weighted-fem-stiffness] -/
theorem integral_weighted_grad_dot_grad (hD : D p ≠ 0) (a u w : Fin 3 → ℝ) :
    ∫ x in triangle p, interp p a x * dot2 (grad (interp p u) x) (grad (interp p w) x) =
      (a 0 + a 1 + a 2) / 3 * cotForm (sideLen p 0) (sideLen p 1) (sideLen p 2) u w := by
  simp_rw [grad_interp]
  rw [integral_mul_const, integral_interp hD, ← area_mul_dot_gradInterp hD]
  ring

/-! ## `fem-area-coordinate-rule`: the rule with the area, and the exponential slice -/

/-- **The area-coordinate rule** `∫_T λ₀ᵃ λ₁ᵇ λ₂ᶜ = 2A a! b! c!/(a + b + c + 2)!`, `A` the Heron
area of the side lengths; for `a + b + c = 2` it is the Galerkin mass.

[formula: fem-area-coordinate-rule] -/
theorem integral_bary_monomial_heron (hD : D p ≠ 0) (a b c : ℕ) :
    ∫ x in triangle p, bary p 0 x ^ a * bary p 1 x ^ b * bary p 2 x ^ c =
      2 * area p * ((a ! * b ! * c ! : ℕ) / ((a + b + c + 2)! : ℕ) : ℝ) :=
  integral_eq_of_monomial hD a b c fun _ => rfl

/-- **The exponential slice formula** that `monomial_integral` (`knoeppel.rs`) evaluates:
`∫_T λ₀^{e₀} λ₁^{e₁} λ₂^{e₂} e^{iΩλ_c} = 2A q! r!/(q + r + 1)! · ∫₀¹ s^{e_c} (1 − s)^{q+r+1} e^{iΩs} ds`,
`q = e_{c+1}`, `r = e_{c+2}`; the one-dimensional integral is `beta_exp(e_c, q + r + 1, Ω)`.

[formula: fem-area-coordinate-rule] -/
theorem integral_bary_slice_exp (hD : D p ≠ 0) (e : Fin 3 → ℕ) (Ω : ℝ) (c : Fin 3) :
    ∫ x in triangle p, (bary p 0 x ^ e 0 * bary p 1 x ^ e 1 * bary p 2 x ^ e 2) •
        Complex.exp (((Ω * bary p c x : ℝ) : ℂ) * Complex.I) =
      (2 * area p * (((e (c + 1))! * (e (c + 2))! : ℕ) / ((e (c + 1) + e (c + 2) + 1)! : ℕ) : ℝ)) •
        ∫ s in (0 : ℝ)..1, (s ^ e c * (1 - s) ^ (e (c + 1) + e (c + 2) + 1)) •
          Complex.exp (((Ω * s : ℝ) : ℂ) * Complex.I) := by
  have h := integral_bary_slice hD e (fun s => Complex.exp (((Ω * s : ℝ) : ℂ) * Complex.I))
    (by fun_prop) c
  rw [h, abs_D_eq_two_mul_area]

/-! ## `bff-conformal-energy-triangle`: the conformal energy of a linear map -/

/-- `|∂̄f|²` for `f = a + ib`: `∂̄f = ½((a_x − b_y) + i(a_y + b_x))`. -/
def dbarSq (a b : ℝ × ℝ → ℝ) (x : ℝ × ℝ) : ℝ :=
  ((grad a x).1 - (grad b x).2) ^ 2 / 4 + ((grad a x).2 + (grad b x).1) ^ 2 / 4

/-- **The conformal energy of a linear map on a counter-clockwise triangle** (Hutchinson; Lévy et
al.; BFF Eq. 15 per triangle): for `f = a + ib` linear with corner values `(aᵢ, bᵢ)`,
`2 ∫_T |∂̄f|² = ½ (∫_T |∇a|² + ∫_T |∇b|²) − ½ sda(f(p₀), f(p₁), f(p₂))`, the Dirichlet energies
being the cotan energies of `a` and `b` and the last term the signed area of the image triangle.

[formula: bff-conformal-energy-triangle] -/
theorem integral_dbarSq (hD : 0 < D p) (ua ub : Fin 3 → ℝ) :
    2 * ∫ x in triangle p, dbarSq (interp p ua) (interp p ub) x =
      (cotEnergy (sideLen p 0) (sideLen p 1) (sideLen p 2) (ua 0) (ua 1) (ua 2) +
          cotEnergy (sideLen p 0) (sideLen p 1) (sideLen p 2) (ub 0) (ub 1) (ub 2)) / 2 -
        sda (ua 0, ub 0) (ua 1, ub 1) (ua 2, ub 2) / 2 := by
  simp only [dbarSq, grad_interp]
  rw [integral_triangle_const hD.ne', smul_eq_mul, ← cotForm_self, ← cotForm_self,
    ← area_mul_dot_gradInterp hD.ne', ← area_mul_dot_gradInterp hD.ne', abs_of_pos hD,
    show area p = D p / 2 by rw [← abs_of_pos hD, abs_D_eq_two_mul_area]; ring]
  simp only [gradInterp, baryGrad, Fin.sum_univ_three, side_zero, side_one, side_two, dot2,
    Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul, Prod.fst_sub,
    Prod.snd_sub]
  set d := D p with hd
  have hd0 : d ≠ 0 := hD.ne'
  field_simp
  rw [hd]
  simp only [D, sda]
  ring

/-- **The conformal energy vanishes exactly for the Cauchy–Riemann equations** `a_x = b_y`,
`a_y = −b_x` (an orientation-preserving similarity, or a constant).

[formula: bff-conformal-energy-triangle] -/
theorem integral_dbarSq_eq_zero_iff (hD : D p ≠ 0) (ua ub : Fin 3 → ℝ) :
    ∫ x in triangle p, dbarSq (interp p ua) (interp p ub) x = 0 ↔
      (gradInterp p ua).1 = (gradInterp p ub).2 ∧ (gradInterp p ua).2 = -(gradInterp p ub).1 := by
  simp only [dbarSq, grad_interp]
  rw [integral_triangle_const hD, smul_eq_mul]
  have hA : 0 < |D p| / 2 := by positivity
  constructor
  · intro h
    have h0 := (mul_eq_zero.1 h).resolve_left hA.ne'
    have h1 : ((gradInterp p ua).1 - (gradInterp p ub).2) ^ 2 = 0 := by
      nlinarith [sq_nonneg ((gradInterp p ua).1 - (gradInterp p ub).2),
        sq_nonneg ((gradInterp p ua).2 + (gradInterp p ub).1)]
    have h2 : ((gradInterp p ua).2 + (gradInterp p ub).1) ^ 2 = 0 := by
      nlinarith [sq_nonneg ((gradInterp p ua).1 - (gradInterp p ub).2),
        sq_nonneg ((gradInterp p ua).2 + (gradInterp p ub).1)]
    constructor
    · linarith [pow_eq_zero_iff (n := 2) (two_ne_zero) |>.1 h1]
    · linarith [pow_eq_zero_iff (n := 2) (two_ne_zero) |>.1 h2]
  · rintro ⟨h1, h2⟩
    rw [h1, h2]
    ring

/-! ## `it-odt-interpolation-error`: Chen–Holst Lemma 4.1 on a triangle -/

open Intrinsic in
/-- The linear interpolant of `u(x) = |x|²`. -/
def odtInterp (p : Fin 3 → ℝ × ℝ) : ℝ × ℝ → ℝ := interp p fun i => nsq (p i)

open Intrinsic in
/-- `u_I − u = Σ_{i<j} λᵢ λⱼ |pᵢ − pⱼ|²` for `u = |x|²` (everywhere in the plane). -/
theorem odtInterp_sub (hD : D p ≠ 0) (x : ℝ × ℝ) :
    odtInterp p x - nsq x = bary p 0 x * bary p 1 x * sqDist2 (p 0) (p 1) +
      bary p 1 x * bary p 2 x * sqDist2 (p 1) (p 2) + bary p 2 x * bary p 0 x * sqDist2 (p 2) (p 0) := by
  have hsum := sum_bary x hD
  have hx := sum_bary_smul x hD
  have hx1 := congrArg Prod.fst hx
  have hx2 := congrArg Prod.snd hx
  simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul] at hx1 hx2
  simp only [odtInterp, interp, Fin.sum_univ_three, nsq, sqDist2]
  set l0 := bary p 0 x
  set l1 := bary p 1 x
  set l2 := bary p 2 x
  rw [← hx1, ← hx2]
  linear_combination (-(l0 * ((p 0).1 ^ 2 + (p 0).2 ^ 2) + l1 * ((p 1).1 ^ 2 + (p 1).2 ^ 2) +
    l2 * ((p 2).1 ^ 2 + (p 2).2 ^ 2)) +
    (l0 * (p 0).1 + l1 * (p 1).1 + l2 * (p 2).1) * (0 : ℝ)) * hsum

theorem bary_nonneg_of_mem (hD : D p ≠ 0) {x : ℝ × ℝ} (hx : x ∈ triangle p) (i : Fin 3) :
    0 ≤ bary p i x := by
  obtain ⟨q, ⟨h1, h2, h3⟩, rfl⟩ := hx
  fin_cases i
  · show 0 ≤ bary p 0 _; rw [bary_triMap_zero hD]; linarith
  · show 0 ≤ bary p 1 _; rw [bary_triMap_one hD]; exact h1
  · show 0 ≤ bary p 2 _; rw [bary_triMap_two hD]; exact h2

open Intrinsic in
/-- **The interpolant of the convex `|x|²` lies above it on the triangle.** -/
theorem nsq_le_odtInterp (hD : D p ≠ 0) {x : ℝ × ℝ} (hx : x ∈ triangle p) : nsq x ≤ odtInterp p x := by
  have h := odtInterp_sub hD x
  have := bary_nonneg_of_mem hD hx
  have h0 := this 0
  have h1 := this 1
  have h2 := this 2
  have := sqDist2_nonneg (p 0) (p 1)
  have := sqDist2_nonneg (p 1) (p 2)
  have := sqDist2_nonneg (p 2) (p 0)
  have : 0 ≤ odtInterp p x - nsq x := by rw [h]; positivity
  linarith

theorem measurableSet_triangle (p : Fin 3 → ℝ × ℝ) : MeasurableSet (triangle p) :=
  (isCompact_triangle p).isClosed.measurableSet

theorem continuous_odtInterp (p : Fin 3 → ℝ × ℝ) : Continuous (odtInterp p) :=
  continuous_interp p _

theorem continuous_nsq : Continuous Intrinsic.nsq := by
  unfold Intrinsic.nsq; fun_prop

open Intrinsic in
/-- **The interpolation error of `|x|²` on a triangle** in closed form:
`∫_T |u_I − u| = (A/12)(ℓ₁₂² + ℓ₂₀² + ℓ₀₁²)`.

[formula: it-odt-interpolation-error] -/
theorem integral_odt_error_eq_sides (hD : D p ≠ 0) :
    ∫ x in triangle p, |odtInterp p x - nsq x| =
      area p / 12 * (sideLen p 0 ^ 2 + sideLen p 1 ^ 2 + sideLen p 2 ^ 2) := by
  rw [setIntegral_congr_fun (measurableSet_triangle p) (fun x hx =>
    abs_of_nonneg (sub_nonneg.2 (nsq_le_odtInterp hD hx)))]
  simp_rw [odtInterp_sub hD]
  have hi : ∀ (i j : Fin 3) (c : ℝ), Integrable (fun x => bary p i x * bary p j x * c)
      (volume.restrict (triangle p)) := fun i j c =>
    integrableOn_triangle p (((continuous_bary p i).mul (continuous_bary p j)).mul continuous_const)
  rw [integral_add ?_ (hi 2 0 _), integral_add (hi 0 1 _) (hi 1 2 _),
    integral_mul_const, integral_mul_const, integral_mul_const, integral_bary_mul_bary hD,
    integral_bary_mul_bary hD, integral_bary_mul_bary hD, sideLen_zero, sideLen_one, sideLen_two,
    dist2_sq, dist2_sq, dist2_sq]
  · simp only [Fin.isValue, Fin.reduceEq, ite_false]
    rw [show sqDist2 (p 2) (p 0) = sqDist2 (p 0) (p 2) by unfold sqDist2; ring]
    ring
  · exact (hi 0 1 _).add (hi 1 2 _)

open Intrinsic in
/-- **Chen–Holst Lemma 4.1 on one triangle** (`d = 2`, uniform density): for `u(x) = |x|²` and its
linear interpolant `u_I`, `∫_T |u_I − u| = (1/3) Σᵢ |pᵢ|² |T| − ∫_T |x|²`.

[formula: it-odt-interpolation-error] -/
theorem integral_abs_odt_error (hD : D p ≠ 0) :
    ∫ x in triangle p, |odtInterp p x - nsq x| =
      (nsq (p 0) + nsq (p 1) + nsq (p 2)) * area p / 3 - ∫ x in triangle p, nsq x := by
  rw [setIntegral_congr_fun (measurableSet_triangle p) (fun x hx =>
    abs_of_nonneg (sub_nonneg.2 (nsq_le_odtInterp hD hx))),
    integral_sub (integrableOn_triangle p (continuous_odtInterp p))
      (integrableOn_triangle p continuous_nsq), odtInterp, integral_interp hD]
  ring

open Intrinsic in
/-- **The star energy of `Intrinsic/Move.lean` is the integral of the interpolant**: with
counter-clockwise triangles `T_k = (x, p_k, p_{k+1})`, `odtEnergy p n x = Σ_k ∫_{T_k} u_I`; with
`integral_abs_odt_error` this is Lemma 4.1 for the star, `Σ_k ∫_{T_k} |u_I − u| + Σ_k ∫_{T_k} |y|²`.

[formula: it-odt-interpolation-error] -/
theorem odtEnergy_eq_sum_integral (q : ℕ → ℝ × ℝ) (n : ℕ) (x : ℝ × ℝ)
    (hpos : ∀ k < n, 0 < sda x (q k) (q (k + 1))) :
    odtEnergy q n x = ∑ k ∈ Finset.range n, ∫ y in triangle ![x, q k, q (k + 1)],
      odtInterp ![x, q k, q (k + 1)] y := by
  unfold odtEnergy
  refine Finset.sum_congr rfl fun k hk => ?_
  have h := hpos k (Finset.mem_range.1 hk)
  have hD : D ![x, q k, q (k + 1)] = sda x (q k) (q (k + 1)) := rfl
  rw [odtInterp, integral_interp (by rw [hD]; exact h.ne'),
    show area ![x, q k, q (k + 1)] = |D ![x, q k, q (k + 1)]| / 2 by
      rw [abs_D_eq_two_mul_area]; ring, hD, abs_of_pos h]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons]
  ring

open Intrinsic in
/-- **Lemma 4.1 for a star** (as a sum over its triangles): `S(x) = Σ_k ∫_{T_k} |u_I − u| +
Σ_k ∫_{T_k} |y|²`, so on a fixed tiling of the star region `S` and the interpolation error differ by
a constant. -/
theorem odtEnergy_eq_sum_error_add (q : ℕ → ℝ × ℝ) (n : ℕ) (x : ℝ × ℝ)
    (hpos : ∀ k < n, 0 < sda x (q k) (q (k + 1))) :
    odtEnergy q n x = ∑ k ∈ Finset.range n, ((∫ y in triangle ![x, q k, q (k + 1)],
      |odtInterp ![x, q k, q (k + 1)] y - nsq y|) + ∫ y in triangle ![x, q k, q (k + 1)], nsq y) := by
  rw [odtEnergy_eq_sum_integral q n x hpos]
  refine Finset.sum_congr rfl fun k hk => ?_
  have h := hpos k (Finset.mem_range.1 hk)
  have hD : D ![x, q k, q (k + 1)] ≠ 0 := h.ne'
  rw [integral_abs_odt_error hD, odtInterp, integral_interp hD]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons]
  ring

/-! ## `app-amr-element-norms`: the element norms of the residual estimator -/

/-- **The `L²` norm of a linear residual**, as `amr::estimate` evaluates it:
`‖R‖²_T = (A/6)(r₀² + r₁² + r₂² + r₀r₁ + r₁r₂ + r₂r₀)`.

[formula: app-amr-element-norms] -/
theorem integral_amr_sq (hD : D p ≠ 0) (r : Fin 3 → ℝ) :
    ∫ x in triangle p, interp p r x ^ 2 =
      area p / 6 * (r 0 * r 0 + r 1 * r 1 + r 2 * r 2 + r 0 * r 1 + r 1 * r 2 + r 2 * r 0) := by
  rw [integral_interp_sq hD, IntrinsicApps.galerkinForm]
  ring

/-- **The oscillation**: `‖R − R̄‖²_T = ‖R‖²_T − A R̄²` with `R̄ = (r₀ + r₁ + r₂)/3`, which is the
mean of `R` over `T` (`integral_interp`).

[formula: app-amr-element-norms] -/
theorem integral_amr_osc (hD : D p ≠ 0) (r : Fin 3 → ℝ) :
    ∫ x in triangle p, (interp p r x - (r 0 + r 1 + r 2) / 3) ^ 2 =
      (∫ x in triangle p, interp p r x ^ 2) - area p * ((r 0 + r 1 + r 2) / 3) ^ 2 := by
  set m := (r 0 + r 1 + r 2) / 3
  have e : ∀ x, (interp p r x - m) ^ 2 = interp p r x ^ 2 - (2 * m) * interp p r x + m ^ 2 :=
    fun x => by ring
  simp_rw [e]
  have i1 : Integrable (fun x => interp p r x ^ 2) (volume.restrict (triangle p)) :=
    integrableOn_triangle p ((continuous_interp p r).pow 2)
  have i2 : Integrable (fun x => 2 * m * interp p r x) (volume.restrict (triangle p)) :=
    integrableOn_triangle p (continuous_const.mul (continuous_interp p r))
  rw [integral_add ?_ (integrableOn_triangle p continuous_const), integral_sub i1 i2,
    integral_const_mul, integral_interp hD, integral_triangle_const hD, smul_eq_mul,
    abs_D_eq_two_mul_area]
  · ring
  · exact i1.sub i2

end TriangleIntegral

/-! ## Meshes: the integrals of all faces, each laid out isometrically -/

namespace TriMesh

open TriangleIntegral Matrix

variable {V F : Type*} (M : TriMesh V F)

/-- An isometric layout of every face of `M` in the plane: `P f` has the corners `a, b, c` of `f`
in this order, nonzero signed area, and the side lengths of `M` (`|p₁p₂| = ℓ_bc`, `|p₂p₀| = ℓ_ca`,
`|p₀p₁| = ℓ_ab`). -/
structure IsLayout (P : F → Fin 3 → ℝ × ℝ) : Prop where
  nondeg : ∀ f, D (P f) ≠ 0
  la : ∀ f, sideLen (P f) 0 = M.la f
  lb : ∀ f, sideLen (P f) 1 = M.lb f
  lc : ∀ f, sideLen (P f) 2 = M.lc f

/-- The values of `u` at the corners `a, b, c` of `f`. -/
def cornerVals (u : V → ℝ) (f : F) : Fin 3 → ℝ := ![u (M.a f), u (M.b f), u (M.c f)]

/-- **Every valid mesh has an isometric layout** (the reference layout of `layout_isometric` for
every face), so the statements below are about the edge lengths alone. -/
theorem exists_isLayout (hM : M.Valid) : ∃ P, M.IsLayout P := by
  refine ⟨fun f => ![(0, 0), (M.lc f, 0), (layoutX (M.lc f) (M.la f) (M.lb f),
    layoutY (M.lc f) (triArea (M.lc f) (M.la f) (M.lb f)))], ⟨fun f => ?_, fun f => ?_,
    fun f => ?_, fun f => ?_⟩⟩ <;>
  have ht : IsTriangle (M.lc f) (M.la f) (M.lb f) := (hM.tri f).rotate.rotate <;>
  obtain ⟨_, h01, h02, h12⟩ := layout_isometric ht rfl
  · show sda _ _ _ ≠ 0
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons]
    rw [sda_layout ht.pos_a.ne']
    exact (mul_pos two_pos (triArea_pos ht)).ne'
  · rw [sideLen_zero]
    simpa only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons] using h12
  · rw [sideLen_one, Intrinsic.dist2_comm]
    simpa only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons] using h02
  · rw [sideLen_two]
    simpa only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons] using h01

variable [DecidableEq V] [Fintype V] [Fintype F]

/-- **The cotan Laplacian is the stiffness matrix**: for every isometric layout of the faces,
`Σ_f ∫_{T_f} |∇u_f|² = uᵀ L u`, `u_f` the linear interpolant of `u` on the face `f`.

[formula: cotan-fem-stiffness] -/
theorem sum_integral_grad_sq (hM : M.Valid) {P : F → Fin 3 → ℝ × ℝ} (hP : M.IsLayout P)
    (u : V → ℝ) :
    ∑ f, ∫ x in triangle (P f), dot2 (grad (interp (P f) (M.cornerVals u f)) x)
        (grad (interp (P f) (M.cornerVals u f)) x) = u ⬝ᵥ (M.laplacian *ᵥ u) := by
  rw [M.dotProduct_laplacian_mulVec_eq_sum_faceEnergy hM]
  refine Finset.sum_congr rfl fun f _ => ?_
  rw [integral_grad_sq (hP.nondeg f), hP.la, hP.lb, hP.lc]
  rfl

omit [DecidableEq V] in
/-- Polarization for a symmetric matrix. -/
theorem dotProduct_mulVec_polarize {A : Matrix V V ℝ} (hA : A.IsSymm) (u w : V → ℝ) :
    u ⬝ᵥ (A *ᵥ w) = ((u + w) ⬝ᵥ (A *ᵥ (u + w)) - (u - w) ⬝ᵥ (A *ᵥ (u - w))) / 4 := by
  have hs : w ⬝ᵥ (A *ᵥ u) = u ⬝ᵥ (A *ᵥ w) := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.vecMul_transpose, hA.eq, dotProduct_comm]
  simp only [Matrix.mulVec_add, Matrix.mulVec_sub, dotProduct_add, add_dotProduct,
    dotProduct_sub, sub_dotProduct, hs]
  ring

/-- **The weighted cotan matrix is the stiffness matrix of a linear coefficient**: for every
isometric layout, `Σ_f ∫_{T_f} a_f ∇u_f·∇w_f = uᵀ L_ā w` with `ā` the face means of `a`
(`TriMesh.weightedLaplacian (M.faceMean a)`, SPEC §3.6, D11).

[formula: cut-weighted-fem-stiffness] -/
theorem sum_integral_weighted_grad_dot_grad (hM : M.Valid) {P : F → Fin 3 → ℝ × ℝ}
    (hP : M.IsLayout P) (a u w : V → ℝ) :
    ∑ f, ∫ x in triangle (P f), interp (P f) (M.cornerVals a f) x *
        dot2 (grad (interp (P f) (M.cornerVals u f)) x) (grad (interp (P f) (M.cornerVals w f)) x) =
      u ⬝ᵥ (M.weightedLaplacian (M.faceMean a) *ᵥ w) := by
  rw [dotProduct_mulVec_polarize (M.weightedLaplacian_isSymm _),
    M.dotProduct_weightedLaplacian_mulVec hM, M.dotProduct_weightedLaplacian_mulVec hM,
    ← Finset.sum_sub_distrib, Finset.sum_div]
  refine Finset.sum_congr rfl fun f _ => ?_
  rw [integral_weighted_grad_dot_grad (hP.nondeg f), hP.la, hP.lb, hP.lc]
  simp only [faceMean, faceEnergy, cotEnergy, cotForm, cornerVals, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Pi.add_apply,
    Pi.sub_apply]
  ring

/-- The quadratic form: `Σ_f ∫_{T_f} a_f |∇u_f|² = uᵀ L_ā u`.

[formula: cut-weighted-fem-stiffness] -/
theorem sum_integral_weighted_grad_sq (hM : M.Valid) {P : F → Fin 3 → ℝ × ℝ}
    (hP : M.IsLayout P) (a u : V → ℝ) :
    ∑ f, ∫ x in triangle (P f), interp (P f) (M.cornerVals a f) x *
        dot2 (grad (interp (P f) (M.cornerVals u f)) x) (grad (interp (P f) (M.cornerVals u f)) x) =
      u ⬝ᵥ (M.weightedLaplacian (M.faceMean a) *ᵥ u) :=
  M.sum_integral_weighted_grad_dot_grad hM hP a u u

end TriMesh

end DDG
