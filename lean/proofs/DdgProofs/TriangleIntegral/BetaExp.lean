import DdgProofs.TriangleIntegral.Fem

/-!
# The one-dimensional integral of the finite element connection Laplacian (`app-vector-beta-exp`)

`monomial_integral` of `crates/ddg-intrinsic-apps/src/vector/knoeppel.rs` reduces the integrals of
the local matrices to `B(p, m, Ω) = ∫₀¹ s^p (1 − s)^m e^{iΩs} ds` (`integral_bary_slice_exp`), which
`beta_exp` evaluates in one of two ways. Both are exact identities, proved here:

* **`|Ω| < 3`, the power series** (`hasSum_betaExp`): `B = Σ_k t_k` with
  `t_k = (iΩ)^k/k! · (p + k)! m!/(p + k + m + 1)!` (the exponential series integrated term by
  term against the Beta integral), `t_0 = p! m!/(p + m + 1)!` and
  `t_{k+1} = t_k · iΩ (p + k + 1)/((k + 1)(p + k + m + 2))` (`betaExpTerm_succ`), the running
  product of the code; for `k + 1 ≥ |Ω|` the terms do not grow, `|t_{k+1}| ≤ (|Ω|/(k + 1)) |t_k|`
  (`norm_betaExpTerm_succ_le`), and the tail after the term `K` is at most
  `|t_K| ρ/(1 − ρ)`, `ρ = |Ω|/(K + 1) < 1` (`norm_tail_betaExp_le`). The code stops after the
  first term with `|t_k| < 10⁻¹⁸ |Σ|`; this bounds the truncation error by
  `10⁻¹⁸ |Σ| ρ/(1 − ρ)` whenever `ρ < 1` there (not a statement about rounding).
* **`|Ω| ≥ 3`, the recurrence** (any `Ω ≠ 0`): `K_k = ∫₀¹ s^k e^{iΩs} ds` has
  `K_0 = (e^{iΩ} − 1)/(iΩ)` and `K_k = (e^{iΩ} − k K_{k−1})/(iΩ)` (`expMoment_zero`,
  `expMoment_succ`, integration by parts), and `B = Σ_{l ≤ m} (−1)^l C(m, l) K_{p+l}`
  (`betaExp_eq_sum_expMoment`), with `C(m, l + 1) = C(m, l)(m − l)/(l + 1)` the code's running
  binomial (`choose_succ_eq`).

With it, `monomial_integral` is `2A q! r!/(q + r + 1)! · B(e_c, q + r + 1, Ω)` exactly
(`integral_bary_monomial_cis`).

What these identities do not cover: rounding (the doc comment of `beta_exp` on the growth of the
terms and the amplification of the recurrence is about rounding), and that the integrands of
SPEC N43 are those of the gauge model (`app-vector-fem-connection-laplacian`, cited).
-/

noncomputable section

namespace DdgProofs

namespace TriangleIntegral

open intervalIntegral
open scoped Nat

/-- `e^{iΩs}`. -/
def cis (Ω s : ℝ) : ℂ := Complex.exp (((Ω * s : ℝ) : ℂ) * Complex.I)

/-- `B(p, m, Ω) = ∫₀¹ s^p (1 − s)^m e^{iΩs} ds`, the integral of `integral_bary_slice_exp`. -/
def betaExp (p m : ℕ) (Ω : ℝ) : ℂ := ∫ s in (0 : ℝ)..1, (s ^ p * (1 - s) ^ m) • cis Ω s

/-- `K_k = ∫₀¹ s^k e^{iΩs} ds`. -/
def expMoment (k : ℕ) (Ω : ℝ) : ℂ := ∫ s in (0 : ℝ)..1, (s ^ k) • cis Ω s

@[fun_prop]
theorem continuous_cis (Ω : ℝ) : Continuous (cis Ω) := by
  unfold cis; fun_prop

theorem hasDerivAt_cis (Ω s : ℝ) : HasDerivAt (cis Ω) (cis Ω s * (Ω * Complex.I)) s := by
  have h : HasDerivAt (fun s : ℝ => ((Ω * s : ℝ) : ℂ) * Complex.I) ((Ω : ℂ) * Complex.I) s := by
    have := (((hasDerivAt_id s).const_mul Ω).ofReal_comp).mul_const Complex.I
    simpa using this
  exact h.cexp

theorem cis_zero (Ω : ℝ) : cis Ω 0 = 1 := by simp [cis]

theorem cis_one (Ω : ℝ) : cis Ω 1 = Complex.exp (Ω * Complex.I) := by simp [cis]

/-- **`monomial_integral` with `beta_exp`**: `∫_T λ₀^{e₀} λ₁^{e₁} λ₂^{e₂} e^{iΩλ_c} =
2A q! r!/(q + r + 1)! · B(e_c, q + r + 1, Ω)`, `q = e_{c+1}`, `r = e_{c+2}`.

[formula: app-vector-beta-exp] -/
theorem integral_bary_monomial_cis {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (e : Fin 3 → ℕ) (Ω : ℝ)
    (c : Fin 3) :
    ∫ x in triangle p, (bary p 0 x ^ e 0 * bary p 1 x ^ e 1 * bary p 2 x ^ e 2) • cis Ω (bary p c x) =
      (2 * area p * (((e (c + 1))! * (e (c + 2))! : ℕ) / ((e (c + 1) + e (c + 2) + 1)! : ℕ) : ℝ)) •
        betaExp (e c) (e (c + 1) + e (c + 2) + 1) Ω :=
  integral_bary_slice_exp hD e Ω c

/-! ## The recurrence -/

/-- **The first moment** `K_0 = (e^{iΩ} − 1)/(iΩ)`. 

[formula: app-vector-beta-exp] -/
theorem expMoment_zero {Ω : ℝ} (hΩ : Ω ≠ 0) :
    expMoment 0 Ω = (Complex.exp (Ω * Complex.I) - 1) / (Ω * Complex.I) := by
  have hI : (Ω : ℂ) * Complex.I ≠ 0 := mul_ne_zero (Complex.ofReal_ne_zero.2 hΩ) Complex.I_ne_zero
  have h := integral_eq_sub_of_hasDerivAt (a := 0) (b := 1) (fun s _ => hasDerivAt_cis Ω s)
    ((by fun_prop : Continuous fun s => cis Ω s * (Ω * Complex.I)).intervalIntegrable _ _)
  rw [intervalIntegral.integral_mul_const, cis_one, cis_zero] at h
  rw [eq_div_iff hI, ← h]
  simp [expMoment]

/-- **The recurrence** `K_{k+1} = (e^{iΩ} − (k + 1) K_k)/(iΩ)` (integration by parts). 

[formula: app-vector-beta-exp] -/
theorem expMoment_succ {Ω : ℝ} (hΩ : Ω ≠ 0) (k : ℕ) :
    expMoment (k + 1) Ω =
      (Complex.exp (Ω * Complex.I) - ((k + 1 : ℕ) : ℂ) * expMoment k Ω) / (Ω * Complex.I) := by
  have hI : (Ω : ℂ) * Complex.I ≠ 0 := mul_ne_zero (Complex.ofReal_ne_zero.2 hΩ) Complex.I_ne_zero
  set G : ℝ → ℂ := fun s => ((s ^ (k + 1) : ℝ) : ℂ) * cis Ω s
  have hG : ∀ s, HasDerivAt G (((((k + 1 : ℕ) : ℝ) * s ^ k : ℝ) : ℂ) * cis Ω s +
      ((s ^ (k + 1) : ℝ) : ℂ) * (cis Ω s * (Ω * Complex.I))) s := by
    intro s
    have h1 : HasDerivAt (fun s : ℝ => ((s ^ (k + 1) : ℝ) : ℂ)) (((((k + 1 : ℕ) : ℝ) * s ^ k : ℝ)) : ℂ) s := by
      have := (hasDerivAt_pow (k + 1) s).ofReal_comp
      simpa using this
    exact h1.mul (hasDerivAt_cis Ω s)
  have h := integral_eq_sub_of_hasDerivAt (a := 0) (b := 1) (fun s _ => hG s)
    ((by fun_prop : Continuous fun s : ℝ => (((((k + 1 : ℕ) : ℝ) * s ^ k : ℝ)) : ℂ) * cis Ω s +
      ((s ^ (k + 1) : ℝ) : ℂ) * (cis Ω s * (Ω * Complex.I))).intervalIntegrable _ _)
  have i1 : IntervalIntegrable (fun s : ℝ => (((((k + 1 : ℕ) : ℝ) * s ^ k : ℝ)) : ℂ) * cis Ω s)
      MeasureTheory.volume 0 1 := (by fun_prop : Continuous fun s : ℝ =>
        (((((k + 1 : ℕ) : ℝ) * s ^ k : ℝ)) : ℂ) * cis Ω s).intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun s : ℝ => ((s ^ (k + 1) : ℝ) : ℂ) * (cis Ω s * (Ω * Complex.I)))
      MeasureTheory.volume 0 1 := (by fun_prop : Continuous fun s : ℝ =>
        ((s ^ (k + 1) : ℝ) : ℂ) * (cis Ω s * (Ω * Complex.I))).intervalIntegrable _ _
  rw [intervalIntegral.integral_add i1 i2] at h
  have e1 : ∫ s in (0 : ℝ)..1, (((((k + 1 : ℕ) : ℝ) * s ^ k : ℝ)) : ℂ) * cis Ω s =
      ((k + 1 : ℕ) : ℂ) * expMoment k Ω := by
    rw [expMoment, ← intervalIntegral.integral_const_mul]
    congr 1; funext s
    rw [Complex.real_smul]; push_cast; ring
  have e2 : ∫ s in (0 : ℝ)..1, ((s ^ (k + 1) : ℝ) : ℂ) * (cis Ω s * (Ω * Complex.I)) =
      (Ω * Complex.I) * expMoment (k + 1) Ω := by
    rw [expMoment, ← intervalIntegral.integral_const_mul]
    congr 1; funext s
    rw [Complex.real_smul]; ring
  rw [e1, e2] at h
  simp only [G, one_pow, Complex.ofReal_one, one_mul, cis_one, ne_eq, Nat.add_eq_zero_iff,
    one_ne_zero, and_false, not_false_eq_true, zero_pow, Complex.ofReal_zero, zero_mul,
    sub_zero] at h
  rw [eq_div_iff hI]
  linear_combination h

/-- **The binomial expansion** of `(1 − s)^m`: `B(p, m, Ω) = Σ_{l ≤ m} (−1)^l C(m, l) K_{p+l}`. 

[formula: app-vector-beta-exp] -/
theorem betaExp_eq_sum_expMoment (p m : ℕ) (Ω : ℝ) :
    betaExp p m Ω = ∑ l ∈ Finset.range (m + 1), ((-1) ^ l * (m.choose l : ℝ)) • expMoment (p + l) Ω := by
  have e : ∀ s : ℝ, (s ^ p * (1 - s) ^ m) • cis Ω s =
      ∑ l ∈ Finset.range (m + 1), ((-1) ^ l * (m.choose l : ℝ)) • ((s ^ (p + l)) • cis Ω s) := by
    intro s
    rw [show (1 - s) = (-s + 1) by ring, add_pow, Finset.mul_sum, Finset.sum_smul]
    refine Finset.sum_congr rfl fun l _ => ?_
    rw [smul_smul, neg_pow, one_pow, mul_one, pow_add]
    congr 1
    ring
  simp_rw [betaExp, e]
  rw [intervalIntegral.integral_finsetSum]
  · refine Finset.sum_congr rfl fun l _ => ?_
    rw [intervalIntegral.integral_smul, expMoment]
  · intro l _
    exact (by fun_prop : Continuous fun s : ℝ =>
      ((-1) ^ l * (m.choose l : ℝ)) • ((s ^ (p + l)) • cis Ω s)).intervalIntegrable _ _

/-- The running binomial of `beta_exp`: `C(m, l + 1) = C(m, l) (m − l)/(l + 1)` for `l ≤ m`. 

[formula: app-vector-beta-exp] -/
theorem choose_succ_eq (m l : ℕ) (hl : l ≤ m) :
    (m.choose (l + 1) : ℝ) = (m.choose l : ℝ) * ((m : ℝ) - l) / (l + 1) := by
  have h := Nat.choose_succ_right_eq m l
  rw [eq_div_iff (by positivity)]
  have : ((m - l : ℕ) : ℝ) = (m : ℝ) - l := by rw [Nat.cast_sub hl]
  rw [← this]
  exact_mod_cast h

/-! ## The power series -/

/-- The term `t_k = (iΩ)^k/k! · (p + k)! m!/(p + k + m + 1)!` of the series. -/
def betaExpTerm (p m : ℕ) (Ω : ℝ) (k : ℕ) : ℂ :=
  ((Ω : ℂ) * Complex.I) ^ k / (k ! : ℂ) * (((p + k)! * m ! : ℕ) / ((p + k + m + 1)! : ℕ) : ℝ)

theorem betaExpTerm_zero (p m : ℕ) (Ω : ℝ) :
    betaExpTerm p m Ω 0 = ((p ! * m ! : ℕ) / ((p + m + 1)! : ℕ) : ℝ) := by
  simp [betaExpTerm]

/-- **The running product**: `t_{k+1} = t_k · iΩ (p + k + 1)/((k + 1)(p + k + m + 2))`. 

[formula: app-vector-beta-exp] -/
theorem betaExpTerm_succ (p m : ℕ) (Ω : ℝ) (k : ℕ) :
    betaExpTerm p m Ω (k + 1) = betaExpTerm p m Ω k *
      ((Ω : ℂ) * Complex.I * (((p + k + 1 : ℕ) : ℝ) / (((k + 1 : ℕ) : ℝ) * ((p + k + m + 2 : ℕ) : ℝ)) : ℝ)) := by
  simp only [betaExpTerm]
  rw [show p + (k + 1) = p + k + 1 by ring, show p + k + 1 + m + 1 = (p + k + m + 1) + 1 by ring,
    Nat.factorial_succ (p + k), Nat.factorial_succ (p + k + m + 1), Nat.factorial_succ k]
  push_cast
  have h1 : (k ! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos k).ne'
  have h2 : ((p + k + m + 1)! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos _).ne'
  have h3 : ((k : ℂ) + 1) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero k
  have h4 : ((p : ℂ) + k + m + 2) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero (p + k + m + 1)
  have h5 : ((p : ℂ) + k + m + 1 + 1) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero (p + k + m + 1)
  field_simp
  ring

/-- **The series of `beta_exp`**: `B(p, m, Ω) = Σ_k t_k` (the exponential series integrated
term by term; each term is a Beta integral). 

[formula: app-vector-beta-exp] -/
theorem hasSum_betaExp (p m : ℕ) (Ω : ℝ) : HasSum (betaExpTerm p m Ω) (betaExp p m Ω) := by
  have hterm : ∀ k : ℕ, ∫ s in (0 : ℝ)..1, (s ^ p * (1 - s) ^ m) •
      ((((Ω * s : ℝ) : ℂ) * Complex.I) ^ k / (k ! : ℂ)) = betaExpTerm p m Ω k := by
    intro k
    have e : ∀ s : ℝ, (s ^ p * (1 - s) ^ m) • ((((Ω * s : ℝ) : ℂ) * Complex.I) ^ k / (k ! : ℂ)) =
        ((s ^ (p + k) * (1 - s) ^ m : ℝ) : ℂ) * (((Ω : ℂ) * Complex.I) ^ k / (k ! : ℂ)) := by
      intro s
      rw [Complex.real_smul]
      push_cast
      ring
    simp_rw [e]
    rw [intervalIntegral.integral_mul_const, intervalIntegral.integral_ofReal,
      integral_pow_mul_sub_pow 1 (p + k) m, one_pow, one_mul, betaExpTerm,
      show p + k + m + 1 = p + k + m + 1 from rfl]
    ring
  have hlim : ∀ s : ℝ, HasSum (fun k : ℕ => (s ^ p * (1 - s) ^ m) •
      ((((Ω * s : ℝ) : ℂ) * Complex.I) ^ k / (k ! : ℂ))) ((s ^ p * (1 - s) ^ m) • cis Ω s) := by
    intro s
    have := NormedSpace.expSeries_div_hasSum_exp (𝔸 := ℂ) (((Ω * s : ℝ) : ℂ) * Complex.I)
    rw [← Complex.exp_eq_exp_ℂ] at this
    exact this.const_smul _
  have h := intervalIntegral.hasSum_integral_of_dominated_convergence (μ := MeasureTheory.volume)
    (a := 0) (b := 1)
    (F := fun k s => (s ^ p * (1 - s) ^ m) • ((((Ω * s : ℝ) : ℂ) * Complex.I) ^ k / (k ! : ℂ)))
    (f := fun s => (s ^ p * (1 - s) ^ m) • cis Ω s)
    (fun k _ => |Ω| ^ k / k !)
    (fun k => (by fun_prop : Continuous fun s : ℝ =>
      (s ^ p * (1 - s) ^ m) • ((((Ω * s : ℝ) : ℂ) * Complex.I) ^ k / (k ! : ℂ))).aestronglyMeasurable)
    ?_ (MeasureTheory.ae_of_all _ fun _ _ => Real.summable_pow_div_factorial |Ω|)
    intervalIntegrable_const (MeasureTheory.ae_of_all _ fun s _ => hlim s)
  · simp only [hterm] at h
    exact h
  · intro k
    refine MeasureTheory.ae_of_all _ fun s hs => ?_
    rw [Set.uIoc_of_le zero_le_one] at hs
    have hs0 : 0 ≤ s := hs.1.le
    have hs1 : s ≤ 1 := hs.2
    have h1s : 0 ≤ 1 - s := by linarith
    have hw : ‖s ^ p * (1 - s) ^ m‖ ≤ 1 := by
      rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
      calc s ^ p * (1 - s) ^ m ≤ 1 * 1 :=
            mul_le_mul (pow_le_one₀ hs0 hs1) (pow_le_one₀ h1s (by linarith)) (by positivity)
              zero_le_one
        _ = 1 := one_mul 1
    have hz : ‖(((Ω * s : ℝ) : ℂ) * Complex.I) ^ k / (k ! : ℂ)‖ ≤ |Ω| ^ k / k ! := by
      rw [norm_div, norm_pow, norm_mul, Complex.norm_I, mul_one, Complex.norm_real,
        Complex.norm_natCast, Real.norm_eq_abs, abs_mul, abs_of_nonneg hs0]
      gcongr
      exact mul_le_of_le_one_right (abs_nonneg _) hs1
    rw [norm_smul]
    calc ‖s ^ p * (1 - s) ^ m‖ * ‖(((Ω * s : ℝ) : ℂ) * Complex.I) ^ k / (k ! : ℂ)‖
        ≤ 1 * (|Ω| ^ k / k !) := by gcongr
      _ = |Ω| ^ k / k ! := one_mul _

/-- **The terms decrease once `k + 1 ≥ |Ω|`**: `|t_{k+1}| ≤ (|Ω|/(k + 1)) |t_k|`. -/
theorem norm_betaExpTerm_succ_le (p m : ℕ) (Ω : ℝ) (k : ℕ) :
    ‖betaExpTerm p m Ω (k + 1)‖ ≤ |Ω| / (k + 1) * ‖betaExpTerm p m Ω k‖ := by
  have hk : (0 : ℝ) < k + 1 := by positivity
  have hr0 : (0 : ℝ) ≤ ((p + k + 1 : ℕ) : ℝ) / (((k + 1 : ℕ) : ℝ) * ((p + k + m + 2 : ℕ) : ℝ)) := by
    positivity
  have hr : ((p + k + 1 : ℕ) : ℝ) / (((k + 1 : ℕ) : ℝ) * ((p + k + m + 2 : ℕ) : ℝ)) ≤ 1 / (k + 1) := by
    push_cast
    rw [div_le_div_iff₀ (by positivity) hk]
    nlinarith
  have hn := norm_nonneg (betaExpTerm p m Ω k)
  have h2 := mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hr (abs_nonneg Ω)) hn
  rw [betaExpTerm_succ]
  simp only [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs]
  rw [abs_of_nonneg hr0]
  have e : |Ω| / (k + 1) * ‖betaExpTerm p m Ω k‖ = ‖betaExpTerm p m Ω k‖ * (|Ω| * (1 / (k + 1))) := by
    ring
  rw [e]
  linarith

/-- **The truncation error**: after the term `K`, with `ρ = |Ω|/(K + 1) < 1`,
`|B − Σ_{k ≤ K} t_k| ≤ |t_K| ρ/(1 − ρ)`. 

[formula: app-vector-beta-exp] -/
theorem norm_tail_betaExp_le (p m : ℕ) (Ω : ℝ) (K : ℕ) (hρ : |Ω| / (K + 1) < 1) :
    ‖betaExp p m Ω - ∑ k ∈ Finset.range (K + 1), betaExpTerm p m Ω k‖ ≤
      ‖betaExpTerm p m Ω K‖ * ((|Ω| / (K + 1)) / (1 - |Ω| / (K + 1))) := by
  set ρ := |Ω| / (K + 1) with hρdef
  have hρ0 : 0 ≤ ρ := by positivity
  -- the terms after `K` are dominated by a geometric sequence
  have hdom : ∀ j : ℕ, ‖betaExpTerm p m Ω (K + 1 + j)‖ ≤ ‖betaExpTerm p m Ω K‖ * ρ ^ (j + 1) := by
    intro j
    induction j with
    | zero =>
      simpa [mul_comm] using norm_betaExpTerm_succ_le p m Ω K
    | succ j ih =>
      have h := norm_betaExpTerm_succ_le p m Ω (K + 1 + j)
      have hle : |Ω| / ((↑(K + 1 + j) : ℝ) + 1) ≤ ρ := by
        rw [hρdef]
        apply div_le_div_of_nonneg_left (abs_nonneg _) (by positivity)
        push_cast; linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)]
      calc ‖betaExpTerm p m Ω (K + 1 + (j + 1))‖ = ‖betaExpTerm p m Ω (K + 1 + j + 1)‖ := rfl
        _ ≤ |Ω| / ((↑(K + 1 + j) : ℝ) + 1) * ‖betaExpTerm p m Ω (K + 1 + j)‖ := h
        _ ≤ ρ * (‖betaExpTerm p m Ω K‖ * ρ ^ (j + 1)) := by gcongr
        _ = ‖betaExpTerm p m Ω K‖ * ρ ^ (j + 1 + 1) := by ring
  have hs := hasSum_betaExp p m Ω
  have htail := (hasSum_nat_add_iff' (K + 1)).2 hs
  have hgeom : HasSum (fun j : ℕ => ‖betaExpTerm p m Ω K‖ * ρ ^ (j + 1))
      (‖betaExpTerm p m Ω K‖ * (ρ / (1 - ρ))) := by
    have := (hasSum_geometric_of_lt_one hρ0 hρ).mul_left (‖betaExpTerm p m Ω K‖ * ρ)
    convert this using 1
    · funext j; ring
    · field_simp
  have hsum : Summable fun j => ‖betaExpTerm p m Ω (j + (K + 1))‖ :=
    Summable.of_nonneg_of_le (fun _ => norm_nonneg _)
      (fun j => by rw [add_comm j]; exact hdom j) hgeom.summable
  rw [← htail.tsum_eq]
  calc ‖∑' j, betaExpTerm p m Ω (j + (K + 1))‖ ≤ ∑' j, ‖betaExpTerm p m Ω (j + (K + 1))‖ :=
        norm_tsum_le_tsum_norm hsum
    _ ≤ ∑' j : ℕ, ‖betaExpTerm p m Ω K‖ * ρ ^ (j + 1) :=
        hsum.tsum_le_tsum (fun j => by rw [add_comm j]; exact hdom j) hgeom.summable
    _ = ‖betaExpTerm p m Ω K‖ * (ρ / (1 - ρ)) := hgeom.tsum_eq

end TriangleIntegral

end DdgProofs
