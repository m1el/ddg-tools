import DDG.SDF.Torus

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

namespace DDG

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

/-! ## The Laplacian (App. D, with our `Δg`) -/

/-- `F` has divergence `d` at `x`: `F` is differentiable at `x` and the trace of its derivative
is `d`. -/
def HasDiv (F : V3 → V3) (d : ℝ) (x : V3) : Prop :=
  ∃ D : V3 →L[ℝ] V3, HasFDerivAt F D x ∧ trace3 D = d

theorem trace3_sum (D : ι → V3 →L[ℝ] V3) : trace3 (∑ i, D i) = ∑ i, trace3 (D i) := by
  simp only [trace3, ContinuousLinearMap.sum_apply, Finset.sum_apply]
  exact Finset.sum_comm

theorem HasDiv.add {F G : V3 → V3} {d1 d2 : ℝ} {x : V3} (h1 : HasDiv F d1 x) (h2 : HasDiv G d2 x) :
    HasDiv (fun y => F y + G y) (d1 + d2) x := by
  obtain ⟨D1, hD1, t1⟩ := h1
  obtain ⟨D2, hD2, t2⟩ := h2
  exact ⟨D1 + D2, hD1.add hD2, by rw [trace3_add, t1, t2]⟩

theorem HasDiv.sum {F : ι → V3 → V3} {d : ι → ℝ} {x : V3} (h : ∀ i, HasDiv (F i) (d i) x) :
    HasDiv (fun y => ∑ i, F i y) (∑ i, d i) x := by
  choose D hD ht using h
  refine ⟨∑ i, D i, ?_, by rw [trace3_sum]; exact Finset.sum_congr rfl fun i _ => ht i⟩
  have := HasFDerivAt.sum (u := Finset.univ) (fun i _ => hD i)
  refine this.congr_of_eventuallyEq (Filter.Eventually.of_forall fun y => ?_)
  simp

/-- Product rule: `div(c F) = c div F + ∇c·F`. -/
theorem HasDiv.smul {c : V3 → ℝ} {gc : V3} {F : V3 → V3} {d : ℝ} {x : V3}
    (hc : HasFDerivAt c (dualOf gc) x) (h : HasDiv F d x) :
    HasDiv (fun y => c y • F y) (c x * d + gc ⬝ᵥ F x) x := by
  obtain ⟨D, hD, t⟩ := h
  refine ⟨c x • D + (dualOf gc).smulRight (F x), hc.smul hD, ?_⟩
  rw [trace3_add, trace3_smul, trace3_smulRight, t, dualOf_apply]

theorem HasDiv.const_smul {F : V3 → V3} {d : ℝ} {x : V3} (k : ℝ) (h : HasDiv F d x) :
    HasDiv (fun y => k • F y) (k * d) x := by
  have := HasDiv.smul (hasFDerivAt_const k x |>.congr_fderiv (by ext; simp) :
    HasFDerivAt (fun _ => k) (dualOf 0) x) h
  simpa using this

/-- `div r̂ = 2/ρ` in `ℝ³` (`(d − 1)/ρ` with `d = 3`). -/
theorem hasDiv_rhat {p x : V3} (hx : x ≠ p) : HasDiv (rhat p) (2 / vnorm (x - p)) x := by
  have hρ : 0 < vnorm (x - p) := vnorm_pos_iff.mpr (sub_ne_zero.mpr hx)
  have hid : HasDiv (fun y => y - p) 3 x :=
    ⟨ContinuousLinearMap.id ℝ V3, (hasFDerivAt_id x).sub_const p, trace3_id⟩
  have hinv := hasFDerivAt_inv_comp (hasFDerivAt_dist hx) hρ.ne'
  have hinv' : HasFDerivAt (fun y => (vnorm (y - p))⁻¹)
      (dualOf (-(vnorm (x - p) ^ 2)⁻¹ • rhat p x)) x := by
    refine hinv.congr_fderiv ?_; ext v; simp [smul_dotProduct]
  have := HasDiv.smul hinv' hid
  have e : 2 / vnorm (x - p) = (vnorm (x - p))⁻¹ * 3 +
      (-(vnorm (x - p) ^ 2)⁻¹ • rhat p x) ⬝ᵥ (x - p) := by
    simp only [rhat, smul_dotProduct, smul_eq_mul]
    rw [← vnorm_mul_self]
    field_simp
    ring
  rw [e]
  exact this

/-- The frozen blend gradient field `F = A + φB` (`blend_gradient`) as a function of `x`. -/
def blendGradField (p : ι → V3) (g : ι → V3 → ℝ) (Gf : ι → V3 → V3) (lam σ : ℝ) (y : V3) : V3 :=
  (∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)))⁻¹ •
      ∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)) • (Gf i y - (lam * g i y) • rhat (p i) y) +
    blend (fun _ => 1) p g lam σ y •
      ((lam * (∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)))⁻¹) •
        ∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)) • rhat (p i) y)

/-- **The Laplacian of the blend** (App. D, with frozen `λ, σ` and our `Δg_i`): if every
`∇g_i = Gf_i` has divergence `Δg_i` at `x` (off the data points), then the gradient field
`F = A + φB` of `blend_gradient` has divergence
`Σ w(Δg − 2λ r̂·∇g + λ g(λ − 2/ρ))/W + A·B + F·B + φ(λ Σ w(2/ρ − λ)/W + B·B)`.

[formula: pat-blend-laplacian] -/
theorem blend_laplacian [Nonempty ι] (p : ι → V3) {g : ι → V3 → ℝ} {Gf : ι → V3 → V3}
    {Lg : ι → ℝ} (lam σ : ℝ) {x : V3} (hx : ∀ i, x ≠ p i)
    (hg : ∀ i, HasFDerivAt (g i) (dualOf (Gf i x)) x) (hG : ∀ i, HasDiv (Gf i) (Lg i) x) :
    let w := fun i => Real.exp (-lam * (vnorm (x - p i) - σ))
    let W := ∑ i, w i
    let A := W⁻¹ • ∑ i, w i • (Gf i x - (lam * g i x) • rhat (p i) x)
    let B := (lam * W⁻¹) • ∑ i, w i • rhat (p i) x
    let φ := blend (fun _ => 1) p g lam σ x
    HasDiv (blendGradField p g Gf lam σ)
      ((∑ i, w i * (Lg i - 2 * lam * (rhat (p i) x ⬝ᵥ Gf i x) +
          g i x * lam * (lam - 2 / vnorm (x - p i)))) / W + A ⬝ᵥ B +
        blendGradField p g Gf lam σ x ⬝ᵥ B +
        φ * (lam * (∑ i, w i * (2 / vnorm (x - p i) - lam)) / W + B ⬝ᵥ B)) x := by
  intro w W A B φ
  have hρ : ∀ i, 0 < vnorm (x - p i) := fun i => vnorm_pos_iff.mpr (sub_ne_zero.mpr (hx i))
  have hr1 : ∀ i, rhat (p i) x ⬝ᵥ rhat (p i) x = 1 := by
    intro i; rw [rhat]; exact RC.unit_dot_self (sub_ne_zero.mpr (hx i))
  -- weights
  have hw : ∀ i, HasFDerivAt (fun y => Real.exp (-lam * (vnorm (y - p i) - σ)))
      (dualOf ((-lam * w i) • rhat (p i) x)) x := by
    intro i
    refine ((((hasFDerivAt_dist (hx i)).sub_const σ).const_mul (-lam)).exp).congr_fderiv ?_
    ext v; simp [smul_dotProduct, w]; ring
  have hW : HasFDerivAt (fun y => ∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)))
      (dualOf (∑ i, (-lam * w i) • rhat (p i) x)) x := hasFDerivAt_sum_grad hw
  have hW0 : W ≠ 0 := (Finset.sum_pos (fun i _ => Real.exp_pos _) Finset.univ_nonempty).ne'
  have hWi : HasFDerivAt (fun y => (∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)))⁻¹)
      (dualOf ((-(W ^ 2)⁻¹) • ∑ i, (-lam * w i) • rhat (p i) x)) x := by
    refine (hasFDerivAt_inv_comp hW hW0).congr_fderiv ?_; ext v; simp [smul_dotProduct, W, w]
  have hgi : ∀ i, HasFDerivAt (fun y => lam * g i y) (dualOf (lam • Gf i x)) x := by
    intro i; refine ((hg i).const_mul lam).congr_fderiv ?_; ext v; simp [smul_dotProduct]
  -- the divergence of the two sums
  have hSA : HasDiv (fun y => ∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)) •
      (Gf i y - (lam * g i y) • rhat (p i) y))
      (∑ i, (w i * (Lg i + (-(lam * g i x * (2 / vnorm (x - p i))) +
        -((lam • Gf i x) ⬝ᵥ rhat (p i) x))) +
        ((-lam * w i) • rhat (p i) x) ⬝ᵥ (Gf i x - (lam * g i x) • rhat (p i) x))) x := by
    refine HasDiv.sum fun i => HasDiv.smul (hw i) ?_
    have h2 := HasDiv.smul (hgi i) (hasDiv_rhat (hx i))
    have h3 := (hG i).add (h2.const_smul (-1))
    convert h3 using 1
    · funext y; simp [sub_eq_add_neg]
    · ring
  have hSB : HasDiv (fun y => ∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)) • rhat (p i) y)
      (∑ i, (w i * (2 / vnorm (x - p i)) + ((-lam * w i) • rhat (p i) x) ⬝ᵥ rhat (p i) x)) x :=
    HasDiv.sum fun i => HasDiv.smul (hw i) (hasDiv_rhat (hx i))
  have hφ := blend_gradient p (g := g) (G := fun i => Gf i x) lam σ hx hg
  have h1 := HasDiv.smul hWi hSA
  have hlW : HasFDerivAt (fun y => lam * (∑ i, Real.exp (-lam * (vnorm (y - p i) - σ)))⁻¹)
      (dualOf (lam • ((-(W ^ 2)⁻¹) • ∑ i, (-lam * w i) • rhat (p i) x))) x := by
    refine (hWi.const_mul lam).congr_fderiv ?_; ext v; simp [smul_dotProduct]
  have h2 := HasDiv.smul hφ (HasDiv.smul hlW hSB)
  have h := h1.add h2
  convert h using 1
  · funext y; rfl
  -- the algebra of App. D
  set SA := ∑ i, w i • (Gf i x - (lam * g i x) • rhat (p i) x) with hSAdef
  set SB := ∑ i, w i • rhat (p i) x with hSBdef
  have eA : ∑ i, (w i * (Lg i + (-(lam * g i x * (2 / vnorm (x - p i))) +
        -((lam • Gf i x) ⬝ᵥ rhat (p i) x))) +
        ((-lam * w i) • rhat (p i) x) ⬝ᵥ (Gf i x - (lam * g i x) • rhat (p i) x)) =
      ∑ i, w i * (Lg i - 2 * lam * (rhat (p i) x ⬝ᵥ Gf i x) +
        g i x * lam * (lam - 2 / vnorm (x - p i))) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [smul_dotProduct, dotProduct_sub, dotProduct_smul, hr1, smul_eq_mul,
      dotProduct_comm (Gf i x) (rhat (p i) x)]
    ring
  have eB : ∑ i, (w i * (2 / vnorm (x - p i)) + ((-lam * w i) • rhat (p i) x) ⬝ᵥ rhat (p i) x) =
      ∑ i, w i * (2 / vnorm (x - p i) - lam) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [smul_dotProduct, hr1, smul_eq_mul]; ring
  have eN : ∑ i, (-lam * w i) • rhat (p i) x = (-lam) • SB := by
    rw [hSBdef, Finset.smul_sum]
    exact Finset.sum_congr rfl fun i _ => by rw [smul_smul]
  change _ = W⁻¹ * _ + _ + (_ + _)
  rw [eA, eN]
  have hW' : (∑ i, Real.exp (-lam * (vnorm (x - p i) - σ))) = W := rfl
  have hφ' : blend (fun _ => 1) p g lam σ x = φ := rfl
  have hF : blendGradField p g Gf lam σ x = W⁻¹ • SA + φ • ((lam * W⁻¹) • SB) := rfl
  have hA : A = W⁻¹ • SA := rfl
  have hB : B = (lam * W⁻¹) • SB := rfl
  rw [hW', hφ', hF, eB, hA, hB]
  simp only [smul_dotProduct, dotProduct_smul, add_dotProduct, dotProduct_add, smul_eq_mul,
    dotProduct_comm SB SA]
  field_simp
  ring

end PAT

end DDG
