import DdgProofs.IntrinsicSimplification.Monotone
import DdgProofs.IntrinsicSimplification.Curvature
import DdgProofs.IntrinsicSimplification.Interpolation
import DdgProofs.IntrinsicSimplification.Transport
import DdgProofs.Oracle.Scalar

/-!
# Intrinsic simplification: the oracle's definitions are the ones of the theorems

`lean/oracle/DdgOracle/Formulas/IntrinsicSimplification.lean` writes the formulas the `iem`
executable runs once, against `Ddg.Scalar`. At `ℝ` (the instance of `DdgProofs/Oracle/Scalar.lean`)
they are the definitions the theorems of this directory are about:

* `iemExtAngle_eq` — the oracle's extended angle is `cornerAngle` (Mathlib's clamped `arccos`) for
  positive sides; so `cornerAngle_sum_pos` and the monotonicity results hold for it;
* `iemNewtonStep_eq` — Eq. 4 (`newtonStep`), given `Θ_i(u)` and the denominator `H(u)`;
* `iemWeight_eq` — Eq. 9 with the flat rule (`weights`);
* `iemMassUpdate_eq`, `iemVectorUpdate_eq` — Eq. 5 and Eq. 8 (`massUpdate`, `vectorUpdate`, pairs
  read as complex numbers);
* `iemProjInterp_eq` — Eq. 10 (`projInterp`).
-/

noncomputable section

namespace DdgProofs

namespace IntrinsicSimplification

open Real
open scoped DdgProofs.Oracle

/-- [formula: iem-flatten-monotone] -/
theorem iemExtAngle_eq {a b c : ℝ} (ha : 0 < a) (hb : 0 < b) (hc : 0 < c) :
    Ddg.iemExtAngle a b c = cornerAngle a b c := by
  unfold Ddg.iemExtAngle
  split_ifs with h1 h2 h3
  · exact (cornerAngle_eq_pi_of_ge hb hc h1).symm
  · rw [Oracle.zero_eq]; exact (cornerAngle_eq_zero_of_ge hb hc ha.le (Or.inl h2)).symm
  · rw [Oracle.zero_eq]; exact (cornerAngle_eq_zero_of_ge hb hc ha.le (Or.inr h3)).symm
  · push Not at h1 h2 h3
    exact Oracle.cornerAngle_eq ⟨h1, by linarith, by linarith⟩

/-- [formula: iem-flatten-newton] -/
theorem iemNewtonStep_eq {ι : Type*} [Fintype ι] (S : ι → StarFace) (target u : ℝ) :
    Ddg.iemNewtonStep u target (angleSum S u) (hessian S u) = newtonStep S target u := rfl

/-- [formula: iem-redistribution-weights] -/
theorem iemWeight_eq {ι : Type*} [Fintype ι] (dK : ι → ℝ) (Ki tol : ℝ) (j : ι) :
    Ddg.iemWeight (dK j) (∑ l, |dK l|) Ki tol (Fintype.card ι : ℝ) = weights dK Ki tol j := by
  unfold Ddg.iemWeight weights
  simp only [Oracle.abs_eq, Oracle.one_eq]

/-- [formula: iem-memoryless-cost] -/
theorem iemMassUpdate_eq (mj a mi : ℝ) : Ddg.iemMassUpdate mj a mi = massUpdate mj a mi := rfl

/-- [formula: iem-center-of-mass-update] -/
theorem iemVectorUpdate_eq (mj a mi : ℝ) (R ti e tj : ℝ × ℝ) :
    (⟨(Ddg.iemVectorUpdate mj a mi R ti e tj).1, (Ddg.iemVectorUpdate mj a mi R ti e tj).2⟩ : ℂ) =
      vectorUpdate mj a mi ⟨R.1, R.2⟩ ⟨ti.1, ti.2⟩ ⟨e.1, e.2⟩ ⟨tj.1, tj.2⟩ := by
  by_cases h : 0 < massUpdate mj a mi
  · have h' : massUpdate mj a mi ≠ 0 := h.ne'
    simp only [Ddg.iemVectorUpdate, vectorUpdate, Oracle.zero_eq, iemMassUpdate_eq, h, ite_true]
    apply Complex.ext <;>
      simp [Ddg.cmul, Complex.div_re, Complex.div_im, Complex.normSq_ofReal] <;>
      field_simp
  · simp only [Ddg.iemVectorUpdate, vectorUpdate, Oracle.zero_eq, iemMassUpdate_eq, h, ite_false]

/-- [formula: iem-projective-interpolation] -/
theorem iemProjInterp_eq (w0 w1 w2 b0 b1 b2 : ℝ) :
    Ddg.iemProjInterp w0 w1 w2 b0 b1 b2 =
      (projInterp ![w0, w1, w2] ![b0, b1, b2] 0, projInterp ![w0, w1, w2] ![b0, b1, b2] 1,
        projInterp ![w0, w1, w2] ![b0, b1, b2] 2) := by
  simp [Ddg.iemProjInterp, projInterp, Fin.sum_univ_three]

end IntrinsicSimplification

end DdgProofs
