import DdgProofs.Developability.Vec

/-!
# Developability: one Jacobi rotation (Rutishauser's formulas)

`crates/ddg-developability/SPEC.md` §3.10 (`eig::sym_eig3`, `eig::sym_eig2`: cyclic Jacobi with
Rutishauser's rotation). The oracle's `Ddg.Dev.jacobiT`, `jacobiC`, `jacobiStep` at `ℝ`.

A symmetric 3x3 matrix is written in the basis order `(p, q, r)` of the rotation plane `(p, q)`
and the third index `r` (`PlaneEntries`); since `(p, q, r)` is a permutation of `(0, 1, 2)` the
statements below cover the three planes `(0,1)`, `(0,2)`, `(1,2)` the code visits. The rotation
is `J` with `J_pp = c`, `J_pq = s`, `J_qp = −s`, `J_qq = c`, `J_rr = 1` (the code updates its
eigenvector columns by `V ← V J`).

* `jacobiT_quadratic` — `t = sgn θ/(|θ| + √(θ² + 1))` solves `t² + 2θt − 1 = 0`;
  `abs_jacobiT_le_one` — `|t| ≤ 1` (it is the root of smaller modulus);
* `rotated_eq_conj` — the explicit entries `rotated c s m` are those of `Jᵀ M J`;
* `jacobiStep_eq_rotated` — for `a_pq ≠ 0` the oracle's update (the code's) is `Jᵀ M J` with
  `c = 1/√(t² + 1)`, `s = t c`; in particular it zeroes `a_pq`;
* `jacobiStep_trace`, `jacobiStep_frobenius` — trace and Frobenius norm are kept;
* `jacobi2_eigen` — the 2x2 case (`sym_eig2`): one rotation diagonalizes, the eigenpairs are
  `(a_pp − t a_pq, (c, −s))` and `(a_qq + t a_pq, (s, c))`, an orthonormal pair.

Not covered: the overflow guard `t = 1/(2θ)` for `|θ| > 1e150`, the power-of-two pre-scaling,
the threshold that sets tiny `a_pq` to zero, convergence of the sweeps.
-/

noncomputable section

namespace DdgProofs

namespace Dev

open scoped DdgProofs.Oracle
open Real Matrix

/-! ## The rotation parameters -/

/-- The oracle's `θ = (a_qq − a_pp)/(2 a_pq)`. -/
def jacobiTheta (app aqq apq : ℝ) : ℝ := (aqq - app) / (2 * apq)

theorem jacobiT_eq (app aqq apq : ℝ) :
    Ddg.Dev.jacobiT app aqq apq =
      (if 0 ≤ jacobiTheta app aqq apq then 1 else -1) /
        (|jacobiTheta app aqq apq| + √(jacobiTheta app aqq apq * jacobiTheta app aqq apq + 1)) := by
  simp [Ddg.Dev.jacobiT, Ddg.Dev.sgn, jacobiTheta]

theorem jacobiC_eq (t : ℝ) : Ddg.Dev.jacobiC t = 1 / √(t * t + 1) := by
  simp [Ddg.Dev.jacobiC]

/-- `sgn θ/(|θ| + √(θ² + 1))` is a root of `t² + 2θt − 1`. -/
theorem sgn_div_quadratic (θ : ℝ) :
    let t := (if 0 ≤ θ then (1 : ℝ) else -1) / (|θ| + √(θ * θ + 1))
    t * t + 2 * θ * t - 1 = 0 := by
  intro t
  have hpos0 : 0 < θ * θ + 1 := by nlinarith [mul_self_nonneg θ]
  have hr : √(θ * θ + 1) * √(θ * θ + 1) = θ * θ + 1 := mul_self_sqrt hpos0.le
  have hr1 : 0 < √(θ * θ + 1) := sqrt_pos.mpr hpos0
  have ht : t = (if 0 ≤ θ then (1 : ℝ) else -1) / (|θ| + √(θ * θ + 1)) := rfl
  clear_value t
  generalize √(θ * θ + 1) = r at hr hr1 ht
  by_cases h : 0 ≤ θ
  · rw [ite_eq_left h, abs_of_nonneg h] at ht
    have hd : 0 < θ + r := by linarith
    rw [ht]
    field_simp
    linear_combination (-1 : ℝ) * hr
  · have hlt : θ < 0 := lt_of_not_ge h
    rw [ite_eq_right h, abs_of_neg hlt] at ht
    have hd : 0 < -θ + r := by linarith
    rw [ht]
    field_simp
    linear_combination (-1 : ℝ) * hr

/-- **Rutishauser's `t` solves `t² + 2θt − 1 = 0`**, `θ = (a_qq − a_pp)/(2a_pq)`. -/
theorem jacobiT_quadratic (app aqq apq : ℝ) :
    let t := Ddg.Dev.jacobiT app aqq apq
    t * t + 2 * jacobiTheta app aqq apq * t - 1 = 0 := by
  rw [jacobiT_eq]
  exact sgn_div_quadratic _

/-- `|t| ≤ 1`: Rutishauser's `t` is the root of smaller modulus (the product of the roots is `−1`),
i.e. the rotation angle is at most `π/4`. -/
theorem abs_jacobiT_le_one (app aqq apq : ℝ) : |Ddg.Dev.jacobiT app aqq apq| ≤ 1 := by
  rw [jacobiT_eq]
  set θ := jacobiTheta app aqq apq
  have h1 : 1 ≤ √(θ * θ + 1) := by
    have := Real.sqrt_le_sqrt (show (1 : ℝ) ≤ θ * θ + 1 by nlinarith [mul_self_nonneg θ])
    rwa [Real.sqrt_one] at this
  have hd : 1 ≤ |θ| + √(θ * θ + 1) := by linarith [abs_nonneg θ]
  have hpos : 0 < |θ| + √(θ * θ + 1) := by linarith
  rw [abs_div, abs_of_pos hpos, div_le_one hpos]
  split_ifs <;> simp <;> linarith

/-- `c² (t² + 1) = 1` for `c = 1/√(t² + 1)`. -/
theorem jacobiC_sq (t : ℝ) : Ddg.Dev.jacobiC t ^ 2 * (t * t + 1) = 1 := by
  have h0 : 0 < t * t + 1 := by nlinarith [mul_self_nonneg t]
  rw [jacobiC_eq, div_pow, sq_sqrt h0.le]
  field_simp

/-! ## The rotated matrix -/

/-- The symmetric matrix of `m` in the basis order `(p, q, r)`. -/
def toMatrix (m : Ddg.Dev.PlaneEntries ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  !![m.app, m.apq, m.arp; m.apq, m.aqq, m.arq; m.arp, m.arq, m.arr]

/-- The Jacobi rotation `J` in the basis order `(p, q, r)`. -/
def rotMatrix (c s : ℝ) : Matrix (Fin 3) (Fin 3) ℝ := !![c, s, 0; -s, c, 0; 0, 0, 1]

/-- The entries of `Jᵀ M J`, written out. -/
def rotated (c s : ℝ) (m : Ddg.Dev.PlaneEntries ℝ) : Ddg.Dev.PlaneEntries ℝ :=
  ⟨c ^ 2 * m.app - 2 * c * s * m.apq + s ^ 2 * m.aqq,
   s ^ 2 * m.app + 2 * s * c * m.apq + c ^ 2 * m.aqq, m.arr,
   c * s * (m.app - m.aqq) + (c ^ 2 - s ^ 2) * m.apq, c * m.arp - s * m.arq,
   s * m.arp + c * m.arq⟩

theorem planeEntries_ext {m m' : Ddg.Dev.PlaneEntries ℝ} (h1 : m.app = m'.app)
    (h2 : m.aqq = m'.aqq) (h3 : m.arr = m'.arr) (h4 : m.apq = m'.apq) (h5 : m.arp = m'.arp)
    (h6 : m.arq = m'.arq) : m = m' := by
  cases m; cases m'; simp_all

/-- `rotated c s m` is `Jᵀ M J`. -/
theorem rotated_eq_conj (c s : ℝ) (m : Ddg.Dev.PlaneEntries ℝ) :
    toMatrix (rotated c s m) = (rotMatrix c s)ᵀ * toMatrix m * rotMatrix c s := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [toMatrix, rotMatrix, rotated, Matrix.mul_apply, Fin.sum_univ_three] <;> ring

/-- `J` is a rotation for `c² + s² = 1`: `Jᵀ J = 1`. -/
theorem rotMatrix_orthogonal {c s : ℝ} (h : c ^ 2 + s ^ 2 = 1) :
    (rotMatrix c s)ᵀ * rotMatrix c s = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rotMatrix, Matrix.mul_apply, Fin.sum_univ_three] <;>
    first | ring1 | linear_combination h

/-- **One Jacobi rotation** (SPEC §3.10, `eig::sym_eig3`): for `a_pq ≠ 0`, the oracle's
`jacobiStep` (the code's update `a_pp ← a_pp − t a_pq`, `a_qq ← a_qq + t a_pq`, `a_pq ← 0`,
`a_rp ← c a_rp − s a_rq`, `a_rq ← s a_rp + c a_rq` with Rutishauser's `t`, `c = 1/√(t² + 1)`,
`s = t c`) equals `Jᵀ M J` for the rotation `J` (`rotated`, `rotated_eq_conj`); so it is an
orthogonal similarity (same eigenvalues), it zeroes `a_pq`, and it keeps the trace and the
Frobenius norm (`jacobiStep_trace`, `jacobiStep_frobenius`).

[formula: dev-jacobi-rotation] -/
theorem jacobiStep_eq_rotated (m : Ddg.Dev.PlaneEntries ℝ) (h : m.apq ≠ 0) :
    let t := Ddg.Dev.jacobiT m.app m.aqq m.apq
    let c := Ddg.Dev.jacobiC t
    Ddg.Dev.jacobiStep m = rotated c (t * c) m ∧ (Ddg.Dev.jacobiStep m).apq = 0 := by
  intro t c
  have hc : c ^ 2 * (t * t + 1) = 1 := jacobiC_sq t
  have hq : t * t + 2 * jacobiTheta m.app m.aqq m.apq * t - 1 = 0 := jacobiT_quadratic _ _ _
  have hθ : 2 * jacobiTheta m.app m.aqq m.apq * m.apq = m.aqq - m.app := by
    unfold jacobiTheta; field_simp
  set θ := jacobiTheta m.app m.aqq m.apq
  have hstep : Ddg.Dev.jacobiStep m =
      ⟨m.app - t * m.apq, m.aqq + t * m.apq, m.arr, 0, c * m.arp - t * c * m.arq,
        t * c * m.arp + c * m.arq⟩ := by
    simp [Ddg.Dev.jacobiStep, t, c]
  refine ⟨?_, by rw [hstep]⟩
  rw [hstep]
  refine planeEntries_ext ?_ ?_ ?_ ?_ ?_ ?_ <;> dsimp only [rotated] <;>
  first
  | ring1
  | linear_combination (-(m.app - t * m.apq)) * hc - c ^ 2 * t * m.apq * hq + c ^ 2 * t ^ 2 * hθ
  | linear_combination (m.app - t * m.apq) * hc + c ^ 2 * t * m.apq * hq - c ^ 2 * t ^ 2 * hθ
  | linear_combination (-(m.aqq + t * m.apq)) * hc + c ^ 2 * t * m.apq * hq - c ^ 2 * t ^ 2 * hθ
  | linear_combination (m.aqq + t * m.apq) * hc - c ^ 2 * t * m.apq * hq + c ^ 2 * t ^ 2 * hθ
  | linear_combination c ^ 2 * m.apq * hq - c ^ 2 * t * hθ

/-- The trace is kept by a Jacobi rotation. -/
theorem jacobiStep_trace (m : Ddg.Dev.PlaneEntries ℝ) :
    (Ddg.Dev.jacobiStep m).app + (Ddg.Dev.jacobiStep m).aqq + (Ddg.Dev.jacobiStep m).arr =
      m.app + m.aqq + m.arr := by
  simp [Ddg.Dev.jacobiStep]

/-- Squared Frobenius norm of a symmetric matrix by its entries. -/
def frobSq (m : Ddg.Dev.PlaneEntries ℝ) : ℝ :=
  m.app ^ 2 + m.aqq ^ 2 + m.arr ^ 2 + 2 * m.apq ^ 2 + 2 * m.arp ^ 2 + 2 * m.arq ^ 2

/-- The Frobenius norm is kept by a Jacobi rotation (`a_pq ≠ 0`). -/
theorem jacobiStep_frobenius (m : Ddg.Dev.PlaneEntries ℝ) (h : m.apq ≠ 0) :
    frobSq (Ddg.Dev.jacobiStep m) = frobSq m := by
  set t := Ddg.Dev.jacobiT m.app m.aqq m.apq
  set c := Ddg.Dev.jacobiC t
  have hc : c ^ 2 * (t * t + 1) = 1 := jacobiC_sq t
  have hq : t * t + 2 * jacobiTheta m.app m.aqq m.apq * t - 1 = 0 := jacobiT_quadratic _ _ _
  have hθ : 2 * jacobiTheta m.app m.aqq m.apq * m.apq = m.aqq - m.app := by
    unfold jacobiTheta; field_simp
  set θ := jacobiTheta m.app m.aqq m.apq
  have hstep : Ddg.Dev.jacobiStep m =
      ⟨m.app - t * m.apq, m.aqq + t * m.apq, m.arr, 0, c * m.arp - t * c * m.arq,
        t * c * m.arp + c * m.arq⟩ := by
    simp [Ddg.Dev.jacobiStep, t, c]
  rw [hstep]
  simp only [frobSq]
  linear_combination 2 * (m.arp ^ 2 + m.arq ^ 2) * hc + 2 * m.apq ^ 2 * hq
    - 2 * t * m.apq * hθ

/-! ## The 2x2 case (`eig::sym_eig2`) -/

/-- **One rotation diagonalizes a symmetric 2x2 matrix** `[[p, b], [b, q]]`, `b ≠ 0`: with
Rutishauser's `t`, `c = 1/√(t² + 1)`, `s = t c`, the vectors `(c, −s)` and `(s, c)` are
orthonormal eigenvectors with the eigenvalues `p − t b` and `q + t b`. -/
theorem jacobi2_eigen (p q b : ℝ) (hb : b ≠ 0) :
    let t := Ddg.Dev.jacobiT p q b
    let c := Ddg.Dev.jacobiC t
    let s := t * c
    (p * c + b * (-s) = (p - t * b) * c ∧ b * c + q * (-s) = (p - t * b) * (-s)) ∧
    (p * s + b * c = (q + t * b) * s ∧ b * s + q * c = (q + t * b) * c) ∧
    c * c + s * s = 1 ∧ c * s + (-s) * c = 0 := by
  intro t c s
  have hc : c ^ 2 * (t * t + 1) = 1 := jacobiC_sq t
  have hq : t * t + 2 * jacobiTheta p q b * t - 1 = 0 := jacobiT_quadratic _ _ _
  have hθ : 2 * jacobiTheta p q b * b = q - p := by
    unfold jacobiTheta; field_simp
  set θ := jacobiTheta p q b
  refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩, ?_, ?_⟩ <;> simp only [s] <;>
  first
  | ring1
  | linear_combination hc
  | linear_combination (-(c * b)) * hq + c * t * hθ

end Dev

end DdgProofs
