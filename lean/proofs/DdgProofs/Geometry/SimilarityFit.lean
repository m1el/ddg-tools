import DdgProofs.Bff.Quality

/-!
# Least-squares similarity fit of two plane point sets

Transcription of `fit_similarity` in `crates/ddg-bff/src/quality.rs`.

## `similarity-fit` — `fit_similarity`

* `simFit_residual_eq_add` — `Σ|a'p_i + b' − q_i|² = Σ|a p_i + b − q_i|² +
  Σ|(a' − a)p_i + (b' − b)|²`;
* `simFit_residual_le` — the fitted `(a, b)` minimizes the sum of squared residuals;
* `simFit_residual_eq_iff` — and is the only minimizer;
* `toC_fitMean`, `fitNorm_eq`, `fitRe_eq`, `fitIm_eq`, `toC_fitCoeff`, `toC_fitTranslation` —
  the real sums of the Rust code are the real and imaginary parts of the complex expressions;
* `fit_similarity_optimal`, `fit_similarity_unique` — hence the evaluated coefficient and
  translation are the unique least-squares optimum;
* `fitNorm_eq_zero_iff` — `norm = 0` iff all `p_i` coincide.


Nothing is claimed about floating-point rounding.
-/

noncomputable section

namespace DDG

open Real

/-! ## Similarity fit -/

section SimilarityFit

open ComplexConjugate

variable {ι : Type*} [Fintype ι]

/-- Centroid `p̄ = (Σ p_i)/n` of finitely many complex numbers. -/
def centroid (p : ι → ℂ) : ℂ := (∑ i, p i) / (Fintype.card ι : ℂ)

/-- The fitted coefficient `a = Σ (q_i − q̄)·conj(p_i − p̄) / Σ |p_i − p̄|²`. -/
def simFitA (p q : ι → ℂ) : ℂ :=
  (∑ i, (q i - centroid q) * conj (p i - centroid p)) /
    ((∑ i, ‖p i - centroid p‖ ^ 2 : ℝ) : ℂ)

/-- The fitted translation `b = q̄ − a·p̄`. -/
def simFitB (p q : ι → ℂ) : ℂ := centroid q - simFitA p q * centroid p

theorem sum_sub_centroid [Nonempty ι] (p : ι → ℂ) : ∑ i, (p i - centroid p) = 0 := by
  have hn : (Fintype.card ι : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul, centroid,
    mul_div_cancel₀ _ hn, sub_self]

/-- The residual of the fit is `a (p_i − p̄) − (q_i − q̄)`. -/
theorem simFit_residual (p q : ι → ℂ) (i : ι) :
    simFitA p q * p i + simFitB p q - q i =
      simFitA p q * (p i - centroid p) - (q i - centroid q) := by
  unfold simFitB
  ring

/-- The residuals of the fit sum to zero. -/
theorem simFit_sum_residual [Nonempty ι] (p q : ι → ℂ) :
    ∑ i, (simFitA p q * p i + simFitB p q - q i) = 0 := by
  simp only [simFit_residual]
  rw [Finset.sum_sub_distrib, ← Finset.mul_sum, sum_sub_centroid, sum_sub_centroid, mul_zero,
    sub_zero]

/-- The residuals of the fit are orthogonal to the points: `Σ r_i·conj(p_i) = 0` (the normal
equation for `a`). -/
theorem simFit_sum_residual_mul_conj [Nonempty ι] (p q : ι → ℂ)
    (hp : ∑ i, ‖p i - centroid p‖ ^ 2 ≠ 0) :
    ∑ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i) = 0 := by
  have hN : ((∑ i, ‖p i - centroid p‖ ^ 2 : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hp
  have hcentred : ∑ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i - centroid p) = 0 := by
    have hterm : ∀ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i - centroid p) =
        simFitA p q * ((p i - centroid p) * conj (p i - centroid p)) -
          (q i - centroid q) * conj (p i - centroid p) := fun i => by
      rw [simFit_residual]
      ring
    have hsq : ∑ i, (p i - centroid p) * conj (p i - centroid p) =
        ((∑ i, ‖p i - centroid p‖ ^ 2 : ℝ) : ℂ) := by
      rw [Complex.ofReal_sum]
      exact Finset.sum_congr rfl fun i _ => by rw [Complex.mul_conj', Complex.ofReal_pow]
    simp only [hterm]
    rw [Finset.sum_sub_distrib, ← Finset.mul_sum, hsq, simFitA, div_mul_cancel₀ _ hN, sub_self]
  have hsplit : ∑ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i) =
      ∑ i, (simFitA p q * p i + simFitB p q - q i) * conj (p i - centroid p) +
        (∑ i, (simFitA p q * p i + simFitB p q - q i)) * conj (centroid p) := by
    rw [Finset.sum_mul, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [map_sub]
    ring
  rw [hsplit, hcentred, simFit_sum_residual, zero_mul, add_zero]

/-- **Pythagoras for the similarity fit.** For any similarity `z ↦ a'z + b'`, the sum of squared
residuals exceeds that of the fitted `(a, b)` by `Σ |(a' − a) p_i + (b' − b)|²`. Hypotheses: at
least one point, and not all `p_i` equal (`Σ |p_i − p̄|² ≠ 0`).

[formula: similarity-fit] -/
theorem simFit_residual_eq_add [Nonempty ι] (p q : ι → ℂ)
    (hp : ∑ i, ‖p i - centroid p‖ ^ 2 ≠ 0) (a' b' : ℂ) :
    ∑ i, ‖a' * p i + b' - q i‖ ^ 2 =
      ∑ i, ‖simFitA p q * p i + simFitB p q - q i‖ ^ 2 +
        ∑ i, ‖(a' - simFitA p q) * p i + (b' - simFitB p q)‖ ^ 2 := by
  have hr0 := simFit_sum_residual p q
  have hr1 := simFit_sum_residual_mul_conj p q hp
  generalize simFitA p q = a at *
  generalize simFitB p q = b at *
  have hsum : ∑ i, (a * p i + b - q i) * conj ((a' - a) * p i + (b' - b)) =
      conj (a' - a) * ∑ i, (a * p i + b - q i) * conj (p i) +
        conj (b' - b) * ∑ i, (a * p i + b - q i) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [map_add, map_mul]
    ring
  have hcross : ∑ i, ((a * p i + b - q i) * conj ((a' - a) * p i + (b' - b))).re = 0 := by
    rw [← Complex.re_sum, hsum, hr0, hr1, mul_zero, mul_zero, add_zero, Complex.zero_re]
  have hterm : ∀ i, ‖a' * p i + b' - q i‖ ^ 2 =
      ‖a * p i + b - q i‖ ^ 2 + ‖(a' - a) * p i + (b' - b)‖ ^ 2 +
        2 * ((a * p i + b - q i) * conj ((a' - a) * p i + (b' - b))).re := by
    intro i
    rw [Complex.sq_norm, Complex.sq_norm, Complex.sq_norm, ← Complex.normSq_add]
    congr 1
    ring
  simp only [hterm]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.mul_sum, hcross, mul_zero,
    add_zero]

/-- **The similarity fit is the least-squares optimum**: no similarity `z ↦ a'z + b'` has a
smaller sum of squared residuals `Σ |a' p_i + b' − q_i|²` than the fitted
`a = Σ (q_i − q̄)·conj(p_i − p̄) / Σ |p_i − p̄|²`, `b = q̄ − a·p̄`.

[formula: similarity-fit] -/
theorem simFit_residual_le [Nonempty ι] (p q : ι → ℂ)
    (hp : ∑ i, ‖p i - centroid p‖ ^ 2 ≠ 0) (a' b' : ℂ) :
    ∑ i, ‖simFitA p q * p i + simFitB p q - q i‖ ^ 2 ≤ ∑ i, ‖a' * p i + b' - q i‖ ^ 2 := by
  rw [simFit_residual_eq_add p q hp a' b']
  exact le_add_of_nonneg_right (Finset.sum_nonneg fun i _ => sq_nonneg _)

/-- **The optimum is unique**: a similarity attains the minimal sum of squared residuals only
if it is the fitted one.

[formula: similarity-fit] -/
theorem simFit_residual_eq_iff [Nonempty ι] (p q : ι → ℂ)
    (hp : ∑ i, ‖p i - centroid p‖ ^ 2 ≠ 0) (a' b' : ℂ) :
    ∑ i, ‖a' * p i + b' - q i‖ ^ 2 = ∑ i, ‖simFitA p q * p i + simFitB p q - q i‖ ^ 2 ↔
      a' = simFitA p q ∧ b' = simFitB p q := by
  constructor
  · intro h
    rw [simFit_residual_eq_add p q hp a' b'] at h
    have h0 : ∑ i, ‖(a' - simFitA p q) * p i + (b' - simFitB p q)‖ ^ 2 = 0 := by linarith
    rw [Finset.sum_eq_zero_iff_of_nonneg fun i _ => sq_nonneg _] at h0
    have hi : ∀ i, (a' - simFitA p q) * p i + (b' - simFitB p q) = 0 := fun i => by
      simpa using h0 i (Finset.mem_univ i)
    have hn : (Fintype.card ι : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
    have hs : ∑ i, ((a' - simFitA p q) * p i + (b' - simFitB p q)) = 0 :=
      Finset.sum_eq_zero fun i _ => hi i
    rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_const, Finset.card_univ,
      nsmul_eq_mul] at hs
    have hc : (a' - simFitA p q) * centroid p + (b' - simFitB p q) = 0 := by
      unfold centroid
      field_simp
      linear_combination hs
    have ha : a' - simFitA p q = 0 := by
      by_contra hne
      apply hp
      refine Finset.sum_eq_zero fun i _ => ?_
      have : (a' - simFitA p q) * (p i - centroid p) = 0 := by linear_combination hi i - hc
      rw [(mul_eq_zero.mp this).resolve_left hne, norm_zero]
      norm_num
    obtain ⟨i⟩ := ‹Nonempty ι›
    have hb := hi i
    rw [ha, zero_mul, zero_add] at hb
    exact ⟨sub_eq_zero.mp ha, sub_eq_zero.mp hb⟩
  · rintro ⟨rfl, rfl⟩
    rfl

/-! ### The real sums of `fit_similarity` -/

/-- `mean` in `fit_similarity`: `(Σ p_i[0] / n, Σ p_i[1] / n)`. -/
def fitMean (p : ι → ℝ × ℝ) : ℝ × ℝ :=
  ((∑ i, (p i).1) / Fintype.card ι, (∑ i, (p i).2) / Fintype.card ι)

/-- `re = Σ (qx·px + qy·py)` on centred coordinates. -/
def fitRe (p q : ι → ℝ × ℝ) : ℝ :=
  ∑ i, (((q i).1 - (fitMean q).1) * ((p i).1 - (fitMean p).1) +
    ((q i).2 - (fitMean q).2) * ((p i).2 - (fitMean p).2))

/-- `im = Σ (qy·px − qx·py)` on centred coordinates. -/
def fitIm (p q : ι → ℝ × ℝ) : ℝ :=
  ∑ i, (((q i).2 - (fitMean q).2) * ((p i).1 - (fitMean p).1) -
    ((q i).1 - (fitMean q).1) * ((p i).2 - (fitMean p).2))

/-- `norm = Σ (px·px + py·py)` on centred coordinates. -/
def fitNorm (p : ι → ℝ × ℝ) : ℝ :=
  ∑ i, (((p i).1 - (fitMean p).1) * ((p i).1 - (fitMean p).1) +
    ((p i).2 - (fitMean p).2) * ((p i).2 - (fitMean p).2))

/-- `(ar, ai) = (re / norm, im / norm)`. -/
def fitCoeff (p q : ι → ℝ × ℝ) : ℝ × ℝ := (fitRe p q / fitNorm p, fitIm p q / fitNorm p)

/-- `translation = (qm[0] − (ar·pm[0] − ai·pm[1]), qm[1] − (ai·pm[0] + ar·pm[1]))`. -/
def fitTranslation (p q : ι → ℝ × ℝ) : ℝ × ℝ :=
  ((fitMean q).1 - ((fitCoeff p q).1 * (fitMean p).1 - (fitCoeff p q).2 * (fitMean p).2),
    (fitMean q).2 - ((fitCoeff p q).2 * (fitMean p).1 + (fitCoeff p q).1 * (fitMean p).2))

/-- The residual length `d = hypot(dx, dy)` of one pair of points under the similarity with
coefficient `ar + i·ai` and translation `tr`: `dx = ar·px − ai·py + tr[0] − qx`,
`dy = ai·px + ar·py + tr[1] − qy`. -/
def fitDist (ar ai : ℝ) (tr p q : ℝ × ℝ) : ℝ :=
  √((ar * p.1 - ai * p.2 + tr.1 - q.1) ^ 2 + (ai * p.1 + ar * p.2 + tr.2 - q.2) ^ 2)

/-- `sum = Σ d·d`, the sum of squared residual lengths (`rms = √(sum / n)`). -/
def fitSqError (ar ai : ℝ) (tr : ℝ × ℝ) (p q : ι → ℝ × ℝ) : ℝ :=
  ∑ i, fitDist ar ai tr (p i) (q i) * fitDist ar ai tr (p i) (q i)

/-- The mean is the centroid.

[formula: similarity-fit] -/
theorem toC_fitMean (p : ι → ℝ × ℝ) : toC (fitMean p) = centroid fun i => toC (p i) := by
  apply Complex.ext <;> simp [toC, fitMean, centroid, Complex.div_natCast_re,
    Complex.div_natCast_im]

/-- `norm = Σ |p_i − p̄|²`.

[formula: similarity-fit] -/
theorem fitNorm_eq (p : ι → ℝ × ℝ) :
    fitNorm p = ∑ i, ‖toC (p i) - centroid fun i => toC (p i)‖ ^ 2 := by
  rw [← toC_fitMean]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Complex.sq_norm, Complex.normSq_apply]
  simp [toC]

/-- `re` is the real part of `Σ (q_i − q̄)·conj(p_i − p̄)`.

[formula: similarity-fit] -/
theorem fitRe_eq (p q : ι → ℝ × ℝ) :
    fitRe p q = (∑ i, (toC (q i) - centroid fun i => toC (q i)) *
      conj (toC (p i) - centroid fun i => toC (p i))).re := by
  rw [← toC_fitMean, ← toC_fitMean, Complex.re_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp [toC]
  ring

/-- `im` is the imaginary part of `Σ (q_i − q̄)·conj(p_i − p̄)`.

[formula: similarity-fit] -/
theorem fitIm_eq (p q : ι → ℝ × ℝ) :
    fitIm p q = (∑ i, (toC (q i) - centroid fun i => toC (q i)) *
      conj (toC (p i) - centroid fun i => toC (p i))).im := by
  rw [← toC_fitMean, ← toC_fitMean, Complex.im_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp [toC]
  ring

/-- `ar + i·ai` is the complex coefficient `a` of the fit.

[formula: similarity-fit] -/
theorem toC_fitCoeff (p q : ι → ℝ × ℝ) :
    toC (fitCoeff p q) = simFitA (fun i => toC (p i)) (fun i => toC (q i)) := by
  unfold simFitA
  rw [← fitNorm_eq]
  apply Complex.ext
  · rw [Complex.div_ofReal_re, ← fitRe_eq]
    rfl
  · rw [Complex.div_ofReal_im, ← fitIm_eq]
    rfl

/-- The translation formula is `b = q̄ − a·p̄` in real and imaginary parts.

[formula: similarity-fit] -/
theorem toC_fitTranslation (p q : ι → ℝ × ℝ) :
    toC (fitTranslation p q) = simFitB (fun i => toC (p i)) (fun i => toC (q i)) := by
  unfold simFitB
  rw [← toC_fitCoeff, ← toC_fitMean, ← toC_fitMean]
  apply Complex.ext <;> simp [toC, fitTranslation]
  ring

/-- The residual length is `|a·p + b − q|` for `a = ar + i·ai`, `b = tr`. -/
theorem fitDist_eq_norm (ar ai : ℝ) (tr p q : ℝ × ℝ) :
    fitDist ar ai tr p q = ‖(⟨ar, ai⟩ : ℂ) * toC p + toC tr - toC q‖ := by
  rw [Complex.norm_eq_sqrt_sq_add_sq]
  unfold fitDist
  congr 1
  simp [toC]
  ring

theorem fitSqError_eq (ar ai : ℝ) (tr : ℝ × ℝ) (p q : ι → ℝ × ℝ) :
    fitSqError ar ai tr p q = ∑ i, ‖(⟨ar, ai⟩ : ℂ) * toC (p i) + toC tr - toC (q i)‖ ^ 2 := by
  unfold fitSqError
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [fitDist_eq_norm, sq]

theorem fitNorm_nonneg (p : ι → ℝ × ℝ) : 0 ≤ fitNorm p :=
  Finset.sum_nonneg fun _ _ => add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)

/-- `norm = 0` (the case in which `fit_similarity` returns `None`) exactly when all points `p_i`
coincide.

[formula: similarity-fit] -/
theorem fitNorm_eq_zero_iff [Nonempty ι] (p : ι → ℝ × ℝ) :
    fitNorm p = 0 ↔ ∀ i j, p i = p j := by
  have hmean : fitNorm p = 0 ↔ ∀ i, p i = fitMean p := by
    unfold fitNorm
    rw [Finset.sum_eq_zero_iff_of_nonneg fun i _ =>
      add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)]
    constructor
    · intro h i
      have hi := h i (Finset.mem_univ i)
      have h1 : (p i).1 - (fitMean p).1 = 0 := by
        nlinarith [mul_self_nonneg ((p i).1 - (fitMean p).1),
          mul_self_nonneg ((p i).2 - (fitMean p).2)]
      have h2 : (p i).2 - (fitMean p).2 = 0 := by
        nlinarith [mul_self_nonneg ((p i).1 - (fitMean p).1),
          mul_self_nonneg ((p i).2 - (fitMean p).2)]
      exact Prod.ext (sub_eq_zero.mp h1) (sub_eq_zero.mp h2)
    · intro h i _
      rw [← h i]
      ring
  rw [hmean]
  constructor
  · intro h i j
    rw [h i, h j]
  · intro h i
    have hn : (Fintype.card ι : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
    unfold fitMean
    rw [Finset.sum_congr rfl fun j _ => congrArg Prod.fst (h j i),
      Finset.sum_congr rfl fun j _ => congrArg Prod.snd (h j i)]
    simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    exact Prod.ext (by field_simp) (by field_simp)

/-- **`fit_similarity` returns the least-squares similarity.** With at least one point and not
all `p_i` equal (`norm ≠ 0`; otherwise the function returns `None`), the evaluated coefficient
`(ar, ai)` and translation give a sum of squared residuals `sum = Σ d·d` that is at most that
of any other coefficient `(ar', ai')` and translation `tr'`. (Hence also the smallest `rms`.)

[formula: similarity-fit] -/
theorem fit_similarity_optimal [Nonempty ι] (p q : ι → ℝ × ℝ) (hp : fitNorm p ≠ 0)
    (ar' ai' : ℝ) (tr' : ℝ × ℝ) :
    fitSqError (fitCoeff p q).1 (fitCoeff p q).2 (fitTranslation p q) p q ≤
      fitSqError ar' ai' tr' p q := by
  rw [fitSqError_eq, fitSqError_eq]
  rw [fitNorm_eq] at hp
  have h := simFit_residual_le (fun i => toC (p i)) (fun i => toC (q i)) hp ⟨ar', ai'⟩ (toC tr')
  rw [← toC_fitCoeff, ← toC_fitTranslation] at h
  exact h

/-- The same with uniqueness: another coefficient and translation give the same sum of squared
residuals only if they are the evaluated ones.

[formula: similarity-fit] -/
theorem fit_similarity_unique [Nonempty ι] (p q : ι → ℝ × ℝ) (hp : fitNorm p ≠ 0)
    (ar' ai' : ℝ) (tr' : ℝ × ℝ)
    (h : fitSqError ar' ai' tr' p q =
      fitSqError (fitCoeff p q).1 (fitCoeff p q).2 (fitTranslation p q) p q) :
    (ar', ai') = fitCoeff p q ∧ tr' = fitTranslation p q := by
  rw [fitSqError_eq, fitSqError_eq] at h
  rw [fitNorm_eq] at hp
  have h' := (simFit_residual_eq_iff (fun i => toC (p i)) (fun i => toC (q i)) hp ⟨ar', ai'⟩
    (toC tr')).mp
  rw [← toC_fitCoeff, ← toC_fitTranslation] at h'
  obtain ⟨ha, hb⟩ := h' h
  exact ⟨toC_injective ha, toC_injective hb⟩

end SimilarityFit

end DDG
