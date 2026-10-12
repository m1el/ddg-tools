import DdgProofs.Optimality

/-!
# Closed meshes, and the area normalization of `λ`

On a closed mesh the cotan Laplacian has zero row and column sums (`L·1 = 0`, `1ᵀL = 0`).

* `Feasible.sum_cone_eq_sum_defect`: the state equation forces `Σ μ = Σ Ω` (which is `2πχ` by
  discrete Gauss–Bonnet, see `DdgProofs.GaussBonnet`);
* `distortion_add_const`, `distortion_add_const_eq`: the distortion as a function of the
  additive constant of `u`, and the optimal constant;
* `KKT.weighted_mean_eq_zero`, `IsMinimizer.weighted_mean_eq_zero`: at the optimum the
  `wE·m`-weighted mean of `u` vanishes.

Scaling:

* `isMinimizer_scaleMass_iff`: multiplying all masses by `c > 0` is the same as dividing `λ`
  by `c`; `normalizeMass_scaleMass`: the mass-normalized problem does not see the factor at all.
-/

noncomputable section

namespace Ddg

open Matrix

namespace ConeProblem

variable {ι : Type*} [Fintype ι] {P : ConeProblem ι} {u φ μ : ι → ℝ}

/-! ## Closed meshes -/

/-- **Total cone angle on a closed mesh.** If the columns of `L` sum to zero (`1ᵀL = 0`) and
`L u = μ − Ω`, then `Σ_i μ_i = Σ_i Ω_i`.

[formula: closed-mesh-cone-sum] -/
theorem sum_cone_eq_sum_defect (hcol : (fun _ => (1 : ℝ)) ᵥ* P.L = 0)
    (h : P.L *ᵥ u = μ - P.Ω) : ∑ i, μ i = ∑ i, P.Ω i := by
  have h1 : (fun _ => (1 : ℝ)) ⬝ᵥ (P.L *ᵥ u) = 0 := by
    rw [Matrix.dotProduct_mulVec, hcol]; simp
  rw [h] at h1
  simp only [dotProduct, Pi.sub_apply, one_mul, Finset.sum_sub_distrib] at h1
  linarith

/-- Every feasible pair on a closed mesh has `Σ μ = Σ Ω`.

[formula: closed-mesh-cone-sum] -/
theorem Feasible.sum_cone_eq_sum_defect (hcol : (fun _ => (1 : ℝ)) ᵥ* P.L = 0)
    (h : P.Feasible u μ) : ∑ i, μ i = ∑ i, P.Ω i :=
  ConeProblem.sum_cone_eq_sum_defect hcol h.1

/-- **The distortion as a function of the additive constant**:
`E(u + c) = E(u) + c·Σ wE m u + (c²/2)·Σ wE m`.

[formula: closed-mesh-optimal-constant] -/
theorem distortion_add_const (P : ConeProblem ι) (u : ι → ℝ) (c : ℝ) :
    P.distortion (fun i => u i + c) =
      P.distortion u + c * ∑ i, P.wE i * P.m i * u i + c ^ 2 / 2 * ∑ i, P.wE i * P.m i := by
  unfold distortion
  simp only [Finset.mul_sum]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- **The optimal additive constant.** If `Σ wE m > 0`, then with
`c₀ = −(Σ wE m u)/(Σ wE m)` (minus the weighted mean of `u`),
`E(u + c) = E(u + c₀) + ½·(Σ wE m)·(c − c₀)²` for every `c`. So `c₀` is the unique minimizer
of the distortion over the additive constant, and `u + c₀` has zero weighted mean.

[formula: closed-mesh-optimal-constant] -/
theorem distortion_add_const_eq (P : ConeProblem ι) (hW : 0 < ∑ i, P.wE i * P.m i) (u : ι → ℝ)
    (c : ℝ) :
    P.distortion (fun i => u i + c) =
      P.distortion (fun i => u i + -(∑ i, P.wE i * P.m i * u i) / ∑ i, P.wE i * P.m i) +
        1 / 2 * (∑ i, P.wE i * P.m i) * (c - -(∑ i, P.wE i * P.m i * u i) / ∑ i, P.wE i * P.m i) ^ 2
    := by
  rw [distortion_add_const, distortion_add_const]
  have hne := hW.ne'
  generalize (∑ i, P.wE i * P.m i * u i) = a
  generalize (∑ i, P.wE i * P.m i) = W at hW hne
  field_simp
  ring

/-- The weighted mean of `u + c₀` is zero for the optimal constant `c₀`. -/
theorem weighted_mean_add_optimal_const (P : ConeProblem ι) (hW : 0 < ∑ i, P.wE i * P.m i)
    (u : ι → ℝ) :
    ∑ i, P.wE i * P.m i * (u i + -(∑ i, P.wE i * P.m i * u i) / ∑ i, P.wE i * P.m i) = 0 := by
  simp only [mul_add, Finset.sum_add_distrib, ← Finset.sum_mul]
  have hne := hW.ne'
  field_simp
  ring

/-- **At a solution of the optimality system on a closed mesh the weighted mean of `u`
vanishes**: if `L·1 = 0` then `Σ_i wE_i m_i u_i = 0`. (Sum the adjoint equation over all
vertices.)

[formula: closed-mesh-optimal-constant] -/
theorem KKT.weighted_mean_eq_zero (hrow : P.L *ᵥ (fun _ => (1 : ℝ)) = 0) (h : P.KKT u φ μ) :
    ∑ i, P.wE i * P.m i * u i = 0 := by
  have h1 : (P.Lᵀ *ᵥ φ) ⬝ᵥ (fun _ => (1 : ℝ)) = 0 := by
    rw [Matrix.mulVec_transpose, ← Matrix.dotProduct_mulVec, hrow]; simp
  rw [h.adjoint] at h1
  simp only [dotProduct, mul_one, Finset.sum_neg_distrib] at h1
  linarith

/-- **At any minimizer on a closed mesh the weighted mean of `u` vanishes**: if `L·1 = 0` and
`(u, μ)` is a global minimizer then `Σ_i wE_i m_i u_i = 0`. (No optimality system assumed: this is
optimality with respect to the additive constant alone.)

[formula: closed-mesh-optimal-constant] -/
theorem IsMinimizer.weighted_mean_eq_zero (hrow : P.L *ᵥ (fun _ => (1 : ℝ)) = 0)
    (h : P.IsMinimizer u μ) : ∑ i, P.wE i * P.m i * u i = 0 := by
  set a := ∑ i, P.wE i * P.m i * u i with ha
  set W := ∑ i, P.wE i * P.m i with hW
  have key : ∀ θ : ℝ, 0 < θ → θ ≤ 1 → 0 ≤ θ * (-(a ^ 2)) + θ ^ 2 * (a ^ 2 / 2 * W) := by
    intro θ _ _
    have h1 := h.2 _ μ (Feasible.add_const hrow h.1 (-(θ * a)))
    unfold objective at h1
    rw [distortion_add_const] at h1
    have : -(θ * a) * a + (-(θ * a)) ^ 2 / 2 * W = θ * (-(a ^ 2)) + θ ^ 2 * (a ^ 2 / 2 * W) := by
      ring
    linarith
  have h2 := nonneg_of_forall_small key
  have h3 : a ^ 2 = 0 := le_antisymm (by linarith) (sq_nonneg a)
  exact pow_eq_zero_iff (two_ne_zero) |>.mp h3

/-! ## Scaling of the masses and the normalization of `λ` -/

/-- The same problem with all masses multiplied by `c`. Under uniform scaling of the mesh by a
factor `s`, `L` and `Ω` do not change and the masses are multiplied by `c = s²`
(`DdgProofs.Triangle`, `DdgProofs.Cotan`). -/
def scaleMass (P : ConeProblem ι) (c : ℝ) : ConeProblem ι :=
  { P with m := fun i => c * P.m i }

/-- The same problem with a different `λ`. -/
def withLam (P : ConeProblem ι) (lam : ℝ) : ConeProblem ι := { P with lam := lam }

/-- The same problem with masses divided by the total mass (= total area), so `Σ m = 1`. -/
def normalizeMass (P : ConeProblem ι) : ConeProblem ι :=
  { P with m := fun i => P.m i / ∑ j, P.m j }

theorem objective_scaleMass (P : ConeProblem ι) {c : ℝ} (hc : c ≠ 0) (u μ : ι → ℝ) :
    (P.scaleMass c).objective u μ = c * (P.withLam (P.lam / c)).objective u μ := by
  simp only [objective, distortion, coneCost, scaleMass, withLam, Finset.mul_sum]
  have : ∑ i, 1 / 2 * (P.wE i * (c * P.m i) * u i ^ 2) =
      ∑ i, c * (1 / 2 * (P.wE i * P.m i * u i ^ 2)) :=
    Finset.sum_congr rfl fun i _ => by ring
  rw [this, mul_add, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum]
  field_simp

/-- **Why `λ` must be normalized by the area.** Multiplying all masses by `c > 0` (scaling the
mesh by `√c`) with `λ` kept fixed gives the same minimizers as keeping the masses and replacing
`λ` by `λ/c`. So an un-normalized `λ` is not scale invariant: it scales like an area.

[formula: lambda-area-normalization] -/
theorem isMinimizer_scaleMass_iff (P : ConeProblem ι) {c : ℝ} (hc : 0 < c) (u μ : ι → ℝ) :
    (P.scaleMass c).IsMinimizer u μ ↔ (P.withLam (P.lam / c)).IsMinimizer u μ := by
  have hobj := P.objective_scaleMass hc.ne'
  constructor
  · rintro ⟨hf, hmin⟩
    refine ⟨hf, fun v ν hv => ?_⟩
    have := hmin v ν hv
    rw [hobj, hobj] at this
    exact le_of_mul_le_mul_left this hc
  · rintro ⟨hf, hmin⟩
    refine ⟨hf, fun v ν hv => ?_⟩
    rw [hobj, hobj]
    exact mul_le_mul_of_nonneg_left (hmin v ν hv) hc.le

/-- **The area-normalized problem is scale invariant**: normalizing the masses by their sum
erases a common factor `c ≠ 0`, so with `Σ m = 1` the same `λ` gives the same problem on a
uniformly scaled mesh.

[formula: lambda-area-normalization] -/
theorem normalizeMass_scaleMass (P : ConeProblem ι) {c : ℝ} (hc : c ≠ 0) :
    (P.scaleMass c).normalizeMass = P.normalizeMass := by
  simp only [normalizeMass, scaleMass, ← Finset.mul_sum]
  congr 1
  funext i
  rw [mul_div_mul_left _ _ hc]

/-- After normalization the masses sum to one (when the total mass is nonzero). -/
theorem sum_normalizeMass (P : ConeProblem ι) (h : ∑ j, P.m j ≠ 0) :
    ∑ i, P.normalizeMass.m i = 1 := by
  simp only [normalizeMass, ← Finset.sum_div]
  exact div_self h

end ConeProblem

end Ddg
