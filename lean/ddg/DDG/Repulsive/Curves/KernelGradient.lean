import DDG.Repulsive.Curves.Kernel

/-!
# Repulsive Curves: the gradient of the tangent-point kernel (SPEC §3.2, D24)

With `d = p − q`, `n = P_T d`, `ν = |T × d|` (`= |n|` for a unit `T`), `n̂ = n/ν`, `ρ = |d|`:

```
∇_p k = α ν^{α−1} n̂ / ρ^β − β k d / ρ²,      ∇_q k = −∇_p k
∇_T k = −α ν^{α−1} (T·d) n̂ / ρ^β               (the part along T^⊥)
```

Here `ν^{α−1} n̂` is written `ν^{α−2} n` (`kgradP`, `kgradT`), which is the same for `ν > 0` and is `0`
at `ν = 0` (Lean's `0 · x = 0`), the value SPEC D24 prescribes.

* `hasDerivAt_kernel`: along any differentiable curves `p(s), q(s), T(s)` with `T` a unit vector and
  `T·T' = 0` at `s = t` (unit-length variations of `T`, e.g. the tangent of a moving edge,
  `tangent_dot_deriv`) and `ν > 0` there,
  `d/ds k(p, q, T) = ∇_p k · (p' − q') + ∇_T k · T'`.
* `tendsto_kgrad`: for `α > 1` both gradients, as functions of `(p, q, T)` with `|T| = 1`, tend to
  `0` at a point with `ν = 0`, `p ≠ q`, and are `0` there: `ν^{α−1} n̂ := 0` is the continuous
  extension (SPEC §3.2; the reference instead cuts at `ν < 1e-10`, D24).
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Real Filter Topology

/-- `∇_p k = α ν^{α−2} P_T d / ρ^β − β k d / ρ²` (`= α ν^{α−1} n̂ / ρ^β − β k d / ρ²` for `ν > 0`).

[formula: rc-kernel-gradient] -/
def kgradP (α β : ℝ) (p q T : V3) : V3 :=
  (α * vnorm (crossProduct T (p - q)) ^ (α - 2) / vnorm (p - q) ^ β) • proj T (p - q)
    - (β * kernel α β p q T / vnorm (p - q) ^ 2) • (p - q)

/-- `∇_T k = −α ν^{α−2} (T·d) P_T d / ρ^β` (`= −α ν^{α−1} (T·d) n̂ / ρ^β` for `ν > 0`).

[formula: rc-kernel-gradient] -/
def kgradT (α β : ℝ) (p q T : V3) : V3 :=
  -(α * vnorm (crossProduct T (p - q)) ^ (α - 2) * (T ⬝ᵥ (p - q)) / vnorm (p - q) ^ β) •
    proj T (p - q)

/-- For `ν > 0` the definitions are the SPEC's formulas with `n̂ = n/ν`, `n = P_T d`.

[formula: rc-kernel-gradient] -/
theorem kgrad_eq_spec {α β : ℝ} {p q T : V3} (hν : 0 < vnorm (crossProduct T (p - q))) :
    kgradP α β p q T =
        (α * vnorm (crossProduct T (p - q)) ^ (α - 1) / vnorm (p - q) ^ β) •
            ((vnorm (crossProduct T (p - q)))⁻¹ • proj T (p - q))
          - (β * kernel α β p q T / vnorm (p - q) ^ 2) • (p - q) ∧
      kgradT α β p q T =
        -(α * vnorm (crossProduct T (p - q)) ^ (α - 1) * (T ⬝ᵥ (p - q)) / vnorm (p - q) ^ β) •
          ((vnorm (crossProduct T (p - q)))⁻¹ • proj T (p - q)) := by
  have h1 : vnorm (crossProduct T (p - q)) ^ (α - 1) =
      vnorm (crossProduct T (p - q)) ^ (α - 2) * vnorm (crossProduct T (p - q)) := by
    rw [show α - 1 = (α - 2) + 1 by ring, Real.rpow_add_one hν.ne']
  constructor
  · rw [kgradP, smul_smul, h1]
    congr 2
    field_simp
  · rw [kgradT, smul_smul, h1]
    congr 1
    field_simp

/-- **Derivative of the kernel along curves** (SPEC §3.2): for curves `p, q, T` with `|T(t)| = 1`,
`T(t)·T'(t) = 0` and `ν = |T × (p − q)| > 0` at `t`,
`d/ds k(p, q, T) = ∇_p k · p' + ∇_q k · q' + ∇_T k · T'` with `∇_q k = −∇_p k`.

[formula: rc-kernel-gradient] -/
theorem hasDerivAt_kernel {α β : ℝ} {p q T : ℝ → V3} {p' q' T' : V3} {t : ℝ}
    (hp : HasDerivAt p p' t) (hq : HasDerivAt q q' t) (hT : HasDerivAt T T' t)
    (hunit : T t ⬝ᵥ T t = 1) (horth : T t ⬝ᵥ T' = 0)
    (hν : crossProduct (T t) (p t - q t) ≠ 0) :
    HasDerivAt (fun s => kernel α β (p s) (q s) (T s))
      (kgradP α β (p t) (q t) (T t) ⬝ᵥ (p' - q') + kgradT α β (p t) (q t) (T t) ⬝ᵥ T') t := by
  have hd : HasDerivAt (fun s => p s - q s) (p' - q') t := hp.sub hq
  have hc := hasDerivAt_cross hT hd
  have hd0 : p t - q t ≠ 0 := by
    intro h; apply hν; rw [h]; simp
  have hνpos : 0 < vnorm (crossProduct (T t) (p t - q t)) := vnorm_pos_iff.mpr hν
  have hρpos : 0 < vnorm (p t - q t) := vnorm_pos_iff.mpr hd0
  have hn := hasDerivAt_vnorm hc hν
  have hr := hasDerivAt_vnorm hd hd0
  have hna := hn.rpow_const (p := α) (Or.inl hνpos.ne')
  have hrb := hr.rpow_const (p := β) (Or.inl hρpos.ne')
  have hk := hna.div hrb (Real.rpow_pos_of_pos hρpos β).ne'
  refine hk.congr_deriv ?_
  set ν := vnorm (crossProduct (T t) (p t - q t)) with hν_def
  set ρ := vnorm (p t - q t) with hρ_def
  set d := p t - q t with hd_def
  have e1 : crossProduct (T t) d ⬝ᵥ (crossProduct T' d + crossProduct (T t) (p' - q')) =
      d ⬝ᵥ (p' - q') - (T t ⬝ᵥ d) * (T t ⬝ᵥ (p' - q')) - (T t ⬝ᵥ d) * (d ⬝ᵥ T') := by
    rw [dotProduct_add, cross_dot_cross, cross_dot_cross, hunit, horth, dotProduct_comm d (T t)]
    ring
  have e2 : kgradP α β (p t) (q t) (T t) ⬝ᵥ (p' - q') =
      α * ν ^ (α - 2) / ρ ^ β * (d ⬝ᵥ (p' - q') - (T t ⬝ᵥ d) * (T t ⬝ᵥ (p' - q')))
        - β * (ν ^ α / ρ ^ β) / ρ ^ 2 * (d ⬝ᵥ (p' - q')) := by
    rw [kgradP, sub_dotProduct, smul_dotProduct, smul_dotProduct, proj_dot, smul_eq_mul,
      smul_eq_mul]
    rfl
  have e3 : kgradT α β (p t) (q t) (T t) ⬝ᵥ T' =
      -(α * ν ^ (α - 2) * (T t ⬝ᵥ d) / ρ ^ β) * (d ⬝ᵥ T') := by
    simp only [kgradT, smul_dotProduct, proj_dot, horth, smul_eq_mul]
    ring
  have r1 : ν ^ (α - 1) = ν ^ (α - 2) * ν := by
    rw [show α - 1 = (α - 2) + 1 by ring, Real.rpow_add_one hνpos.ne']
  have r2 : ν ^ α = ν ^ (α - 2) * ν ^ 2 := by
    rw [show α = (α - 2) + 2 by ring, Real.rpow_add hνpos, Real.rpow_two]
    ring_nf
  have r3 : ρ ^ (β - 1) = ρ ^ β / ρ := Real.rpow_sub_one hρpos.ne' β
  have hB : 0 < ρ ^ β := Real.rpow_pos_of_pos hρpos β
  rw [e2, e3, e1, r1, r2, r3]
  field_simp
  ring

/-! ## Continuity at `ν = 0` -/

theorem continuous_vnorm : Continuous vnorm := by
  unfold vnorm; fun_prop

/-- The sup norm of `Fin 3 → ℝ` is at most the Euclidean norm. -/
theorem norm_le_vnorm (n : V3) : ‖n‖ ≤ vnorm n := by
  refine (pi_norm_le_iff_of_nonneg (vnorm_nonneg n)).mpr fun i => ?_
  rw [Real.norm_eq_abs]
  apply Real.abs_le_sqrt
  simp only [dotProduct, Fin.sum_univ_three]
  fin_cases i <;> simp <;> nlinarith [mul_self_nonneg (n 0), mul_self_nonneg (n 1),
    mul_self_nonneg (n 2)]

/-- `|n|^{α−2} n → 0` as `n → 0`, for `α > 1`. -/
theorem tendsto_rpow_smul {X : Type*} [TopologicalSpace X] {n : X → V3} {x₀ : X} {α : ℝ}
    (hα : 1 < α) (hn : ContinuousAt n x₀) (h0 : n x₀ = 0) :
    Tendsto (fun x => vnorm (n x) ^ (α - 2) • n x) (𝓝 x₀) (𝓝 0) := by
  refine squeeze_zero_norm' (a := fun x => vnorm (n x) ^ (α - 1))
    (Eventually.of_forall fun x => ?_) ?_
  · rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (vnorm_nonneg _) _)]
    rcases (vnorm_nonneg (n x)).eq_or_lt with h | h
    · have : n x = 0 := vnorm_eq_zero_iff.mp h.symm
      simp only [this, norm_zero, mul_zero]
      exact Real.rpow_nonneg (vnorm_nonneg _) _
    · calc vnorm (n x) ^ (α - 2) * ‖n x‖ ≤ vnorm (n x) ^ (α - 2) * vnorm (n x) :=
            mul_le_mul_of_nonneg_left (norm_le_vnorm _) (Real.rpow_nonneg (vnorm_nonneg _) _)
        _ = vnorm (n x) ^ (α - 1) := by
            rw [← Real.rpow_add_one h.ne']; ring_nf
  · have hc : Tendsto (fun x => vnorm (n x)) (𝓝 x₀) (𝓝 0) := by
      have := (continuous_vnorm.continuousAt (x := n x₀)).tendsto.comp hn.tendsto
      rwa [h0, vnorm_zero] at this
    have := (Real.continuousAt_rpow_const 0 (α - 1) (Or.inr (by linarith))).tendsto.comp hc
    rwa [Real.zero_rpow (by linarith)] at this

/-- Unit tangents with two points: the domain of the kernel gradient. -/
abbrev UnitTriple := {x : V3 × V3 × V3 // x.2.2 ⬝ᵥ x.2.2 = 1}

/-- **Continuity at `ν = 0`** (SPEC §3.2, D24): for `α > 1`, at `(p₀, q₀, T₀)` with `|T₀| = 1`,
`p₀ ≠ q₀` and `T₀ × (p₀ − q₀) = 0`, both `∇_p k` and `∇_T k` are `0` and tend to `0` as `(p, q, T)`
tends to `(p₀, q₀, T₀)` with `|T| = 1`.

[formula: rc-kernel-gradient] -/
theorem tendsto_kgrad {α : ℝ} (β : ℝ) (hα : 1 < α) {x₀ : UnitTriple}
    (hρ : x₀.1.1 ≠ x₀.1.2.1) (hν : crossProduct x₀.1.2.2 (x₀.1.1 - x₀.1.2.1) = 0) :
    Tendsto (fun x : UnitTriple => kgradP α β x.1.1 x.1.2.1 x.1.2.2) (𝓝 x₀) (𝓝 0) ∧
      kgradP α β x₀.1.1 x₀.1.2.1 x₀.1.2.2 = 0 ∧
      Tendsto (fun x : UnitTriple => kgradT α β x.1.1 x.1.2.1 x.1.2.2) (𝓝 x₀) (𝓝 0) ∧
      kgradT α β x₀.1.1 x₀.1.2.1 x₀.1.2.2 = 0 := by
  -- the pieces as functions on the domain
  set d : UnitTriple → V3 := fun x => x.1.1 - x.1.2.1 with hd
  set n : UnitTriple → V3 := fun x => proj x.1.2.2 (d x) with hn
  have hcont_d : Continuous d := by rw [hd]; fun_prop
  have hcont_n : Continuous n := by rw [hn, hd]; unfold proj; fun_prop
  have hν_eq : ∀ x : UnitTriple, vnorm (crossProduct x.1.2.2 (d x)) = vnorm (n x) :=
    fun x => (cross_vnorm_sq_unit x.2 (d x)).2
  have hn0 : n x₀ = 0 := by
    rw [← vnorm_eq_zero_iff, ← hν_eq]; exact vnorm_eq_zero_iff.mpr hν
  have hρ0 : 0 < vnorm (d x₀) := vnorm_pos_iff.mpr (sub_ne_zero.mpr hρ)
  have hρc : ContinuousAt (fun x => vnorm (d x)) x₀ := (continuous_vnorm.comp hcont_d).continuousAt
  have hρβ : ∀ γ : ℝ, Tendsto (fun x => vnorm (d x) ^ γ) (𝓝 x₀) (𝓝 (vnorm (d x₀) ^ γ)) :=
    fun γ => (Real.continuousAt_rpow_const _ γ (Or.inl hρ0.ne')).tendsto.comp hρc.tendsto
  have hsing := tendsto_rpow_smul hα hcont_n.continuousAt hn0
  -- the kernel tends to 0
  have hk : Tendsto (fun x : UnitTriple => kernel α β x.1.1 x.1.2.1 x.1.2.2) (𝓝 x₀) (𝓝 0) := by
    have hν' : Tendsto (fun x => vnorm (n x)) (𝓝 x₀) (𝓝 0) := by
      have := (continuous_vnorm.comp hcont_n).continuousAt (x := x₀)
      rw [ContinuousAt, Function.comp_apply, hn0, vnorm_zero] at this
      exact this
    have hνα := (Real.continuousAt_rpow_const 0 α (Or.inr (by linarith))).tendsto.comp hν'
    rw [Real.zero_rpow (by linarith)] at hνα
    have := hνα.div (hρβ β) (Real.rpow_pos_of_pos hρ0 β).ne'
    rw [zero_div] at this
    refine this.congr fun x => ?_
    simp only [kernel, Function.comp_apply, Pi.div_apply, ← hν_eq, hd]
  have hk0 : kernel α β x₀.1.1 x₀.1.2.1 x₀.1.2.2 = 0 := by
    rw [kernel, hν, vnorm_zero, Real.zero_rpow (by linarith), zero_div]
  have hTd : Continuous fun x : UnitTriple => x.1.2.2 ⬝ᵥ d x := by rw [hd]; fun_prop
  refine ⟨?_, ?_, ?_, ?_⟩
  · have hA : Tendsto (fun x => α / vnorm (d x) ^ β) (𝓝 x₀) (𝓝 (α / vnorm (d x₀) ^ β)) :=
      tendsto_const_nhds.div (hρβ β) (Real.rpow_pos_of_pos hρ0 β).ne'
    have hB : Tendsto (fun x : UnitTriple => β * kernel α β x.1.1 x.1.2.1 x.1.2.2 / vnorm (d x) ^ 2)
        (𝓝 x₀) (𝓝 0) := by
      have h2 : Tendsto (fun x => vnorm (d x) ^ 2) (𝓝 x₀) (𝓝 (vnorm (d x₀) ^ 2)) :=
        hρc.tendsto.pow 2
      have := (hk.const_mul β).div h2 (by positivity)
      rwa [mul_zero, zero_div] at this
    have := (hA.smul hsing).sub (hB.smul hcont_d.continuousAt.tendsto)
    rw [smul_zero, zero_smul, sub_zero] at this
    refine this.congr fun x => ?_
    simp only [kgradP, smul_smul, hn, hd, (cross_vnorm_sq_unit x.2 (x.1.1 - x.1.2.1)).2]
    congr 2
    ring
  · simp only [kgradP, hk0, hν, vnorm_zero]
    have : proj x₀.1.2.2 (x₀.1.1 - x₀.1.2.1) = 0 := hn0
    simp [this]
  · have hA : Tendsto (fun x : UnitTriple => -(α * (x.1.2.2 ⬝ᵥ d x) / vnorm (d x) ^ β)) (𝓝 x₀)
        (𝓝 (-(α * (x₀.1.2.2 ⬝ᵥ d x₀) / vnorm (d x₀) ^ β))) :=
      ((tendsto_const_nhds.mul hTd.continuousAt.tendsto).div (hρβ β)
        (Real.rpow_pos_of_pos hρ0 β).ne').neg
    have := hA.smul hsing
    rw [smul_zero] at this
    refine this.congr fun x => ?_
    simp only [kgradT, smul_smul, hn, hd, (cross_vnorm_sq_unit x.2 (x.1.1 - x.1.2.1)).2]
    congr 1
    ring
  · have : proj x₀.1.2.2 (x₀.1.1 - x₀.1.2.1) = 0 := hn0
    simp [kgradT, this]

end RC

end DDG
