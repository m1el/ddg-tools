import DdgProofs.PointsAsTori.Torus

/-!
# Points as Tori: the blend (Eq. 25–26; App. D; SPEC §7.4)

For a finite neighbor set `ι`, data points `p i`, per-torus SDFs `g i` and parameters `λ, σ`:

```
ρ_i = |x − p_i|,   w_i = exp(−λ(ρ_i − σ)),   φ = Σ w_i g_i / Σ w_i                  (Eq. 25)
```

* **Shift** (`blend_shift`): `φ` does not depend on `σ` (Eq. 25 = the self-normalized sum of
  Eq. 1 with the same `λ`), also with area weights.
* **Weight range** (`blend_weight_range`): with `σ = ½ max ρ_j` and `λ = C/σ` every exponent is in
  `[−C, C]`; the farthest point has exponent `−C`. (The SPEC's remark "the max weight is
  `e^{λσ} ≥ 1` for the nearest point" holds only when `ρ_min ≤ σ`; in general the weights are
  `≥ e^{−C}`, so `W ≥ |𝒩| e^{−C} > 0`, `blend_weight_sum_pos`.)
* **Convexity** (`blend_convex`).
* **Gradient**: a general softmax-average rule (`blend_hasFDerivAt`): for weights `exp(e_i)`
  with `∇e_i = E_i`, `∇φ = Σ w_i (∇g_i + (g_i − φ) E_i)/W`. With `λ, σ` frozen this is App. D
  (`blend_gradient`, `∇φ = A + φB`); with `σ = ½ρ_j`, `λ = C/σ` for a fixed argmax `j` it is the
  frozen gradient plus `(∇λ) Σ w(σ − ρ)(g − φ)/W` (`blend_gradient_exact`).
-/

noncomputable section

namespace DdgProofs

namespace PAT

open Matrix Real Finset
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff vnorm_smul
  dot_self_nonneg)

variable {ι : Type*} [Fintype ι]

/-- The blend of Eq. 25 with shift `σ` and (optional) area weights `a`. -/
def blend (a : ι → ℝ) (p : ι → V3) (g : ι → V3 → ℝ) (lam σ : ℝ) (x : V3) : ℝ :=
  (∑ i, a i * Real.exp (-lam * (vnorm (x - p i) - σ)) * g i x) /
    ∑ i, a i * Real.exp (-lam * (vnorm (x - p i) - σ))

/-- **The shift cancels** (Eq. 25 = Eq. 1 discretized, the self-normalized sum with weights
`a_i e^{−λρ_i}`): `φ_σ = φ_0` for every `σ`.

[formula: pat-blend-shift] -/
theorem blend_shift (a : ι → ℝ) (p : ι → V3) (g : ι → V3 → ℝ) (lam σ : ℝ) (x : V3) :
    blend a p g lam σ x = blend a p g lam 0 x := by
  unfold blend
  have h : ∀ i, Real.exp (-lam * (vnorm (x - p i) - σ)) =
      Real.exp (lam * σ) * Real.exp (-lam * (vnorm (x - p i) - 0)) := by
    intro i; rw [← Real.exp_add]; congr 1; ring
  simp only [h]
  have e1 : ∑ i, a i * (Real.exp (lam * σ) * Real.exp (-lam * (vnorm (x - p i) - 0))) * g i x =
      Real.exp (lam * σ) * ∑ i, a i * Real.exp (-lam * (vnorm (x - p i) - 0)) * g i x := by
    rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun i _ => by ring
  have e2 : ∑ i, a i * (Real.exp (lam * σ) * Real.exp (-lam * (vnorm (x - p i) - 0))) =
      Real.exp (lam * σ) * ∑ i, a i * Real.exp (-lam * (vnorm (x - p i) - 0)) := by
    rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun i _ => by ring
  rw [e1, e2, mul_div_mul_left _ _ (Real.exp_pos _).ne']

/-- **Weight range** (Eq. 26, the per-query rule): with `σ = ½ max ρ_j > 0` (here: `ρ_i ≤ 2σ`
for all `i`) and `λ = C/σ`, every exponent `−λ(ρ_i − σ)` is in `[−C, C]`, so
`w_i ∈ [e^{−C}, e^{C}]`; a point attaining the maximum (`ρ_j = 2σ`) has exponent exactly `−C`.

[formula: pat-blend-weight-range] -/
theorem blend_weight_range {C σ ρ : ℝ} (hC : 0 ≤ C) (hσ : 0 < σ) (h0 : 0 ≤ ρ) (h2 : ρ ≤ 2 * σ) :
    -C ≤ -(C / σ) * (ρ - σ) ∧ -(C / σ) * (ρ - σ) ≤ C ∧
      Real.exp (-C) ≤ Real.exp (-(C / σ) * (ρ - σ)) ∧
      Real.exp (-(C / σ) * (ρ - σ)) ≤ Real.exp C ∧
      (ρ = 2 * σ → -(C / σ) * (ρ - σ) = -C) := by
  have e : -(C / σ) * (ρ - σ) = C * (1 - ρ / σ) := by field_simp; ring
  have h1 : ρ / σ ≤ 2 := by rw [div_le_iff₀ hσ]; linarith
  have h1' : 0 ≤ ρ / σ := div_nonneg h0 hσ.le
  have hlo : -C ≤ -(C / σ) * (ρ - σ) := by rw [e]; nlinarith
  have hhi : -(C / σ) * (ρ - σ) ≤ C := by rw [e]; nlinarith
  refine ⟨hlo, hhi, Real.exp_le_exp.mpr hlo, Real.exp_le_exp.mpr hhi, fun h => ?_⟩
  rw [e, h]; field_simp; ring

/-- With positive weights bounded below by `e^{−C}`, the denominator is at least
`|ι| e^{−C} > 0` (no underflow to `0`).

[formula: pat-blend-weight-range] -/
theorem blend_weight_sum_pos [Nonempty ι] {C : ℝ} (w : ι → ℝ) (hw : ∀ i, Real.exp (-C) ≤ w i) :
    (Fintype.card ι : ℝ) * Real.exp (-C) ≤ ∑ i, w i ∧ 0 < ∑ i, w i := by
  have h := Finset.card_nsmul_le_sum Finset.univ w (Real.exp (-C)) (fun i _ => hw i)
  simp only [Finset.card_univ, nsmul_eq_mul] at h
  refine ⟨h, lt_of_lt_of_le ?_ h⟩
  have : (0 : ℝ) < Fintype.card ι := by exact_mod_cast Fintype.card_pos
  positivity

/-- **The blend is a convex combination**: with positive weights, `min g ≤ φ ≤ max g`.

[formula: pat-blend-convex] -/
theorem blend_convex [Nonempty ι] (w : ι → ℝ) (hw : ∀ i, 0 < w i) (g : ι → ℝ) {lo hi : ℝ}
    (hlo : ∀ i, lo ≤ g i) (hhi : ∀ i, g i ≤ hi) :
    lo ≤ (∑ i, w i * g i) / ∑ i, w i ∧ (∑ i, w i * g i) / ∑ i, w i ≤ hi := by
  have hW : 0 < ∑ i, w i := Finset.sum_pos (fun i _ => hw i) Finset.univ_nonempty
  constructor
  · rw [le_div_iff₀ hW, Finset.mul_sum]
    exact Finset.sum_le_sum fun i _ => by nlinarith [hw i, hlo i]
  · rw [div_le_iff₀ hW, Finset.mul_sum]
    exact Finset.sum_le_sum fun i _ => by nlinarith [hw i, hhi i]

/-! ## Gradients -/

/-- Quotient rule for gradients. -/
theorem hasFDerivAt_div_grad {N W : V3 → ℝ} {dN dW : V3} {x : V3}
    (hN : HasFDerivAt N (dualOf dN) x) (hW : HasFDerivAt W (dualOf dW) x) (h0 : W x ≠ 0) :
    HasFDerivAt (fun y => N y / W y) (dualOf ((W x)⁻¹ • dN - (N x / W x ^ 2) • dW)) x := by
  have := hN.mul (hasFDerivAt_inv_comp hW h0)
  refine (this.congr_fderiv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
  · ext h
    simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply, dualOf_apply,
      smul_eq_mul, sub_dotProduct, smul_dotProduct]
    field_simp
    ring
  · simp [div_eq_mul_inv]

theorem hasFDerivAt_sum_grad {f : ι → V3 → ℝ} {d : ι → V3} {x : V3}
    (h : ∀ i, HasFDerivAt (f i) (dualOf (d i)) x) :
    HasFDerivAt (fun y => ∑ i, f i y) (dualOf (∑ i, d i)) x := by
  have := HasFDerivAt.sum (u := Finset.univ) (fun i _ => h i)
  refine (this.congr_fderiv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
  · ext v; simp [sum_dotProduct]
  · simp

/-- **Gradient of a softmax average** (the core of App. D): for weights `w_i = exp(e_i)` with
`∇e_i = E_i` and functions `g_i` with `∇g_i = G_i`,
`∇(Σ w g/Σ w) = Σ w_i (G_i + (g_i − φ) E_i) / Σ w`. -/
theorem blend_hasFDerivAt [Nonempty ι] {e g : ι → V3 → ℝ} {E G : ι → V3} {x : V3}
    (he : ∀ i, HasFDerivAt (e i) (dualOf (E i)) x) (hg : ∀ i, HasFDerivAt (g i) (dualOf (G i)) x) :
    HasFDerivAt (fun y => (∑ i, Real.exp (e i y) * g i y) / ∑ i, Real.exp (e i y))
      (dualOf ((∑ i, Real.exp (e i x))⁻¹ • ∑ i, Real.exp (e i x) •
        (G i + (g i x - (∑ j, Real.exp (e j x) * g j x) / ∑ j, Real.exp (e j x)) • E i))) x := by
  have hw : ∀ i, HasFDerivAt (fun y => Real.exp (e i y)) (dualOf (Real.exp (e i x) • E i)) x := by
    intro i
    refine ((he i).exp).congr_fderiv ?_
    ext v; simp [smul_dotProduct]
  have hwg : ∀ i, HasFDerivAt (fun y => Real.exp (e i y) * g i y)
      (dualOf (Real.exp (e i x) • G i + g i x • (Real.exp (e i x) • E i))) x := by
    intro i
    refine ((hw i).mul (hg i)).congr_fderiv ?_
    ext v; simp [add_dotProduct, smul_dotProduct]; try ring
  have hW0 : (∑ i, Real.exp (e i x)) ≠ 0 :=
    (Finset.sum_pos (fun i _ => Real.exp_pos _) Finset.univ_nonempty).ne'
  refine (hasFDerivAt_div_grad (hasFDerivAt_sum_grad hwg) (hasFDerivAt_sum_grad hw) hW0).congr_fderiv ?_
  congr 1
  set W := ∑ i, Real.exp (e i x)
  set N := ∑ i, Real.exp (e i x) * g i x
  ext k
  simp only [Pi.sub_apply, Pi.smul_apply, Pi.add_apply, Finset.sum_apply, smul_eq_mul,
    Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  field_simp
  ring

/-- The unit vector `r̂_i = (x − p_i)/|x − p_i|`. -/
def rhat (p x : V3) : V3 := (vnorm (x - p))⁻¹ • (x - p)

theorem hasFDerivAt_dist {p x : V3} (h : x ≠ p) :
    HasFDerivAt (fun y => vnorm (y - p)) (dualOf (rhat p x)) x := by
  have := hasFDerivAt_vnorm_comp ((hasFDerivAt_id x).sub_const p) (sub_ne_zero.mpr h)
  refine this.congr_fderiv ?_
  ext v; simp [rhat, smul_dotProduct]

/-- **The blend gradient with `λ`, `σ` frozen** (App. D, the reference's `∇φ`): away from the
data points, `∇φ = A + φB` with `A = Σ w(∇g − λ r̂ g)/W` and `B = λ Σ w r̂/W`.

[formula: pat-blend-gradient] -/
theorem blend_gradient [Nonempty ι] (p : ι → V3) {g : ι → V3 → ℝ} {G : ι → V3} (lam σ : ℝ)
    {x : V3} (hx : ∀ i, x ≠ p i) (hg : ∀ i, HasFDerivAt (g i) (dualOf (G i)) x) :
    HasFDerivAt (blend (fun _ => 1) p g lam σ)
      (dualOf ((∑ i, Real.exp (-lam * (vnorm (x - p i) - σ)))⁻¹ •
          ∑ i, Real.exp (-lam * (vnorm (x - p i) - σ)) • (G i - (lam * g i x) • rhat (p i) x) +
        blend (fun _ => 1) p g lam σ x •
          ((lam * (∑ i, Real.exp (-lam * (vnorm (x - p i) - σ)))⁻¹) •
            ∑ i, Real.exp (-lam * (vnorm (x - p i) - σ)) • rhat (p i) x))) x := by
  have he : ∀ i, HasFDerivAt (fun y => -lam * (vnorm (y - p i) - σ))
      (dualOf ((-lam) • rhat (p i) x)) x := by
    intro i
    refine (((hasFDerivAt_dist (hx i)).sub_const σ).const_mul (-lam)).congr_fderiv ?_
    ext v; simp [smul_dotProduct]
  have h := blend_hasFDerivAt he hg
  refine (h.congr_fderiv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
  · congr 1
    simp only [blend, one_mul]
    ext k
    simp only [Pi.sub_apply, Pi.smul_apply, Pi.add_apply, Finset.sum_apply, smul_eq_mul,
      Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    ring
  · simp [blend]

/-- **The exact blend gradient** (ours, SPEC §7.4, `GradientMode::Exact`): where the farthest
neighbor `j` is fixed (`σ(x) = ½|x − p_j|`, `λ(x) = C/σ(x)`), the gradient of
`x ↦ φ_{λ(x), σ(x)}(x)` is the frozen gradient plus `(∂λ/∂x) Σ w(σ − ρ)(g − φ)/W` with
`∂λ/∂x = −(C/σ²) ½ r̂_j`; the contribution of `∇σ` cancels.

[formula: pat-blend-gradient-exact] -/
theorem blend_gradient_exact [Nonempty ι] (p : ι → V3) {g : ι → V3 → ℝ} {G : ι → V3} (C : ℝ)
    (j : ι) {x : V3} (hx : ∀ i, x ≠ p i) (hg : ∀ i, HasFDerivAt (g i) (dualOf (G i)) x) :
    let σ := vnorm (x - p j) / 2
    let lam := C / σ
    let w := fun i => Real.exp (-lam * (vnorm (x - p i) - σ))
    let W := ∑ i, w i
    let φ := blend (fun _ => 1) p g lam σ x
    HasFDerivAt (fun y => blend (fun _ => 1) p g (C / (vnorm (y - p j) / 2)) (vnorm (y - p j) / 2) y)
      (dualOf (W⁻¹ • ∑ i, w i • (G i - (lam * g i x) • rhat (p i) x) +
        φ • ((lam * W⁻¹) • ∑ i, w i • rhat (p i) x) +
        (W⁻¹ * ∑ i, w i * (σ - vnorm (x - p i)) * (g i x - φ)) •
          ((-(C / σ ^ 2) * (1 / 2)) • rhat (p j) x))) x := by
  intro σ lam w W φ
  have hσpos : 0 < σ := by
    have := vnorm_pos_iff.mpr (sub_ne_zero.mpr (hx j)); simp only [σ]; positivity
  have hS : HasFDerivAt (fun y => vnorm (y - p j) / 2) (dualOf ((1 / 2 : ℝ) • rhat (p j) x)) x := by
    refine (((hasFDerivAt_dist (hx j)).const_mul (1 / 2 : ℝ)).congr_fderiv ?_).congr_of_eventuallyEq
      (Filter.Eventually.of_forall fun y => by ring)
    ext v; simp [smul_dotProduct]
  have hL : HasFDerivAt (fun y => C / (vnorm (y - p j) / 2))
      (dualOf ((-(C / σ ^ 2) * (1 / 2)) • rhat (p j) x)) x := by
    have := (hasFDerivAt_inv_comp hS hσpos.ne').const_mul C
    refine (this.congr_fderiv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
    · ext v; simp [smul_dotProduct, σ]; ring
    · simp [div_eq_mul_inv]
  set dL := (-(C / σ ^ 2) * (1 / 2)) • rhat (p j) x
  set dS := (1 / 2 : ℝ) • rhat (p j) x
  have he : ∀ i, HasFDerivAt (fun y => -(C / (vnorm (y - p j) / 2)) *
      (vnorm (y - p i) - vnorm (y - p j) / 2))
      (dualOf (-(vnorm (x - p i) - σ) • dL - lam • (rhat (p i) x - dS))) x := by
    intro i
    refine (hL.neg.mul ((hasFDerivAt_dist (hx i)).sub hS)).congr_fderiv ?_
    ext v
    simp only [_root_.add_apply, _root_.smul_apply, _root_.neg_apply, _root_.sub_apply,
      dualOf_apply, smul_eq_mul, sub_dotProduct, smul_dotProduct, neg_dotProduct, Pi.neg_apply,
      Pi.sub_apply, lam, σ]
    ring
  have h := blend_hasFDerivAt he hg
  refine (h.congr_fderiv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
  · congr 1
    have hW0 : W ≠ 0 := (Finset.sum_pos (fun i _ => Real.exp_pos _) Finset.univ_nonempty).ne'
    have hφ : φ = (∑ i, w i * g i x) / W := by simp only [φ, blend, one_mul]; rfl
    -- `Σ w (g − φ) = 0`
    have hzero : ∑ i, w i * (g i x - φ) = 0 := by
      rw [hφ]; simp only [mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul]
      field_simp; ring
    change (∑ i, w i)⁻¹ • ∑ i, w i • (G i + (g i x - (∑ j, w j * g j x) / ∑ j, w j) •
      (-(vnorm (x - p i) - σ) • dL - lam • (rhat (p i) x - dS))) = _
    rw [← hφ]
    ext k
    simp only [Pi.sub_apply, Pi.smul_apply, Pi.add_apply, Pi.neg_apply, Finset.sum_apply,
      smul_eq_mul]
    have hsplit : ∑ i, w i * (G i k + (g i x - φ) *
        (-(vnorm (x - p i) - σ) * dL k - lam * (rhat (p i) x k - dS k))) =
        ∑ i, w i * (G i k - lam * g i x * rhat (p i) x k) +
        (φ * lam) * ∑ i, w i * rhat (p i) x k +
        dL k * ∑ i, w i * (σ - vnorm (x - p i)) * (g i x - φ) +
        (lam * dS k) * ∑ i, w i * (g i x - φ) := by
      simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun i _ => by ring
    rw [hsplit, hzero]
    ring
  · simp [blend]

end PAT

end DdgProofs
