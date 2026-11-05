import DdgProofs.IntrinsicSimplification.Removal
import DdgProofs.IntrinsicSimplification.Wasserstein

/-!
# Intrinsic simplification: error vectors, transport, vector interpolation

`crates/ddg-intrinsic-simplification/SPEC.md` §4.3, §7.3, §13; the paper Sec. 3.2.2, Sec. 5.1.1,
5.2, 5.3 (Eq. 7, 8), Sec. 7.3.1.

The plane is `ℂ`. A tangent frame at a vertex of a flat region is a unit complex number `ρ`
(the direction of the vertex's reference direction): a vector with coordinates `v` in the frame
is the plane vector `ρ v`.

* **`iem-center-of-mass-update`** — Eq. 7 in the plane: if `t_i`, `t_j` point from `x_i`, `x_j` to
  the centres of mass `c_i`, `c_j` of the masses `m_i`, `m_j`, then `t̃_j` of Eq. 7 points from `x_j`
  to the centre of mass of `α m_i` at `c_i` and `m_j` at `c_j` (`eq7_centerOfMass`), which is the
  centre of mass of the combined distributions (`centerOfMass_combine`). Eq. 8, written in the
  frames (`R_ij = ρ_j⁻¹ ρ_i`, `e_ji = ρ_j⁻¹(x_i − x_j)`), is Eq. 7 (`eq8_eq_eq7`); with SPEC §4.3's
  `vectorUpdate` (`t̃_j = t_j` unless `m̃_j > 0`) the same holds whenever `m̃_j > 0`.
* **`iem-transport-rotation`** — in a flat region, if the edge `ij` leaves `i` at the angle `φ_ij`
  of the frame of `i` and arrives at `j` from the angle `φ_ji` of the frame of `j` (it leaves `j`
  in the opposite direction), then `ρ_j⁻¹ ρ_i = e^{i(φ_ji + π − φ_ij)} = R_ij`, so `R_ij` maps the
  coordinates of a vector in the frame of `i` to its coordinates in the frame of `j`
  (`transport_rotation`). On a curved surface this is the definition of the discrete connection
  (`app-vector-transport`, cited); here it is a theorem about the flattened star.
* **`iem-vector-interpolation`** — Sec. 7.3.1: laying out the face `ijk` with `i` at the origin and
  `ij` along the positive axis, the oriented edges `ij`, `jk`, `ki` point at the angles
  `β_ij = 0`, `β_jk = π − θ_j`, `β_ki = (π − θ_j) + (π − θ_k)` (`beta_angles`, from the lengths), and
  `z_i = e^{iβ_ij}/w_ij` turns the unit vector `w_ij` into `e^{iβ_ij}` and is a rotation
  (`z_rotation`).
-/

noncomputable section

namespace DDG

namespace IntrinsicSimplification

open Real Finset Complex

/-! ## `iem-center-of-mass-update` -/

/-- Eq. 7: `t̃_j = (α m_i (t_i + e_ji) + m_j t_j)/(α m_i + m_j)`. -/
def eq7 (a mi mj : ℝ) (ti tj e : ℂ) : ℂ :=
  ((a * mi : ℝ) * (ti + e) + (mj : ℂ) * tj) / ((a * mi + mj : ℝ) : ℂ)

/-- **Eq. 7 moves the error vector to the new centre of mass.** With `c_i = x_i + t_i`,
`c_j = x_j + t_j`, `e_ji = x_i − x_j` and `α m_i + m_j ≠ 0`:
`x_j + t̃_j = (α m_i c_i + m_j c_j)/(α m_i + m_j)`.

[formula: iem-center-of-mass-update] -/
theorem eq7_centerOfMass {a mi mj : ℝ} (h : a * mi + mj ≠ 0) (xi xj ti tj : ℂ) :
    xj + eq7 a mi mj ti tj (xi - xj) =
      ((a * mi : ℝ) * (xi + ti) + (mj : ℂ) * (xj + tj)) / ((a * mi + mj : ℝ) : ℂ) := by
  have h' : ((a * mi + mj : ℝ) : ℂ) ≠ 0 := by exact_mod_cast h
  unfold eq7
  field_simp
  push_cast
  ring

variable {Y : Type*} [Fintype Y]

/-- The centre of mass of the masses `w` at the points `p`. -/
def centerOfMass (w : Y → ℝ) (p : Y → ℂ) : ℂ := (∑ y, (w y : ℂ) * p y) / ((∑ y, w y : ℝ) : ℂ)

/-- **Combining two distributions**: the centre of mass of `α w₁ + w₂` is the mass-weighted mean of
the two centres (`M₁ = Σ w₁`, `M₂ = Σ w₂`, both non-zero; the formula needs `α M₁ + M₂ ≠ 0` to
mean anything). So Eq. 7 keeps
`x_j + t_j` at the centre of mass of all the fine mass merged into `j` (in the plane).

[formula: iem-center-of-mass-update] -/
theorem centerOfMass_combine {w₁ w₂ : Y → ℝ} (p : Y → ℂ) (a : ℝ) (h₁ : ∑ y, w₁ y ≠ 0)
    (h₂ : ∑ y, w₂ y ≠ 0) :
    centerOfMass (fun y => a * w₁ y + w₂ y) p =
      ((a * ∑ y, w₁ y : ℝ) * centerOfMass w₁ p + ((∑ y, w₂ y : ℝ) : ℂ) * centerOfMass w₂ p) /
        ((a * ∑ y, w₁ y + ∑ y, w₂ y : ℝ) : ℂ) := by
  have h₁' : ((∑ y, w₁ y : ℝ) : ℂ) ≠ 0 := by exact_mod_cast h₁
  have h₂' : ((∑ y, w₂ y : ℝ) : ℂ) ≠ 0 := by exact_mod_cast h₂
  have hs : ∑ y, (a * w₁ y + w₂ y) = a * ∑ y, w₁ y + ∑ y, w₂ y := by
    rw [Finset.sum_add_distrib, Finset.mul_sum]
  unfold centerOfMass
  rw [hs]
  congr 1
  have e1 : ((a * ∑ y, w₁ y : ℝ) : ℂ) * ((∑ y, (w₁ y : ℂ) * p y) / ((∑ y, w₁ y : ℝ) : ℂ)) =
      (a : ℂ) * ∑ y, (w₁ y : ℂ) * p y := by
    push_cast at h₁' ⊢; field_simp
  have e2 : ((∑ y, w₂ y : ℝ) : ℂ) * ((∑ y, (w₂ y : ℂ) * p y) / ((∑ y, w₂ y : ℝ) : ℂ)) =
      ∑ y, (w₂ y : ℂ) * p y := by
    field_simp
  rw [e1, e2, Finset.mul_sum, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl; intro y _; push_cast; ring

/-- **Eq. 8 is Eq. 7 written in the frames.** With frames `ρ_i`, `ρ_j` (non-zero), the transport
`R_ij = ρ_j⁻¹ ρ_i`, the edge vector `e_ji = ρ_j⁻¹ (x_i − x_j)` in the frame of `j`, and the vectors
`t_i`, `t_j` in their frames: the plane vector of Eq. 8's `t̃_j` is Eq. 7's in the plane. Hence it
points from `x_j` to the combined centre of mass (`eq7_centerOfMass`).

[formula: iem-center-of-mass-update] -/
theorem eq8_eq_eq7 {a mi mj : ℝ} {ρi ρj : ℂ} (hj : ρj ≠ 0) (xi xj ti tj : ℂ) :
    ρj * eq7 a mi mj (ρj⁻¹ * ρi * ti) tj (ρj⁻¹ * (xi - xj)) =
      eq7 a mi mj (ρi * ti) (ρj * tj) (xi - xj) := by
  unfold eq7
  field_simp

/-- SPEC §4.3's update (`vectorUpdate`, Eq. 8 with `R` and `e` given) is Eq. 7's formula when
`m̃_j > 0`. -/
theorem vectorUpdate_eq_eq7 {mj a mi : ℝ} (h : 0 < massUpdate mj a mi) (R ti e tj : ℂ) :
    vectorUpdate mj a mi R ti e tj = eq7 a mi mj (R * ti) tj e := by
  unfold vectorUpdate eq7
  simp only [h, ite_true]
  unfold massUpdate
  push_cast
  ring_nf

/-! ## `iem-transport-rotation` -/

/-- **`R_ij = e^{i(φ_ji + π − φ_ij)}` changes frames along the edge.** In a flat region let the
frames of `i`, `j` be `ρ_i`, `ρ_j` and let the edge `ij` have the plane direction `δ ≠ 0`: it leaves
`i` at the angle `φ_ij` of the frame of `i` (`ρ_i e^{iφ_ij} = δ`) and leaves `j` towards `i` at
the angle `φ_ji` of the frame of `j` (`ρ_j e^{iφ_ji} = −δ`). Then `ρ_j⁻¹ ρ_i = R_ij`, and a plane
vector with coordinates `v` in the frame of `i` has the coordinates `R_ij v` in the frame of `j`.

[formula: iem-transport-rotation] -/
theorem transport_rotation {ρi ρj δ : ℂ} {φij φji : ℝ} (hδ : δ ≠ 0)
    (hi : ρi * exp (φij * I) = δ) (hj : ρj * exp (φji * I) = -δ) (v : ℂ) :
    ρj⁻¹ * ρi = exp ((φji + π - φij : ℝ) * I) ∧
      ρi * v = ρj * (exp ((φji + π - φij : ℝ) * I) * v) := by
  have hρj : ρj ≠ 0 := by
    rintro rfl
    have : -δ = 0 := by rw [← hj]; simp
    exact hδ (neg_eq_zero.mp this)
  have hρi : ρi = δ * exp (-(φij * I)) := by
    rw [← hi, mul_assoc, ← Complex.exp_add]; simp
  have hρj' : ρj = -δ * exp (-(φji * I)) := by
    rw [← hj, mul_assoc, ← Complex.exp_add]; simp
  have key : ρj⁻¹ * ρi = exp ((φji + π - φij : ℝ) * I) := by
    rw [hρj', hρi]
    have hpi : exp ((π : ℂ) * I) = -1 := Complex.exp_pi_mul_I
    push_cast
    rw [show ((φji : ℂ) + π - φij) * I = φji * I + π * I + -(φij * I) by ring, Complex.exp_add,
      Complex.exp_add, hpi]
    have he : exp (-(φji * I)) ≠ 0 := Complex.exp_ne_zero _
    rw [Complex.exp_neg (φji * I)]
    field_simp
  refine ⟨key, ?_⟩
  rw [← key, ← mul_assoc, ← mul_assoc, mul_inv_cancel₀ hρj, one_mul]

/-! ## `iem-vector-interpolation` -/

/-- **The angles `β` of Sec. 7.3.1.** The face `ijk` (sides `a = ℓ_jk`, `b = ℓ_ki`, `c = ℓ_ij`, a
strict triangle) laid out with `i = 0`, `j = (c, 0)`, `k = b (cos θ_i, sin θ_i)`: the oriented
edges `jk` and `ki` are `a` and `b` times the unit vectors at the angles `β_jk = π − θ_j` and
`β_ki = (π − θ_j) + (π − θ_k)`; the edge `ij` is at `β_ij = 0`.

[formula: iem-vector-interpolation] -/
theorem beta_angles {a b c : ℝ} (h : IsTriangle a b c) :
    let J := polar c 0
    let K := polar b (cornerAngle a b c)
    (K.1 - J.1, K.2 - J.2) = polar a (π - cornerAngle b c a) ∧
      ((0 : ℝ) - K.1, (0 : ℝ) - K.2) =
        polar b ((π - cornerAngle b c a) + (π - cornerAngle c a b)) := by
  intro J K
  have ha := h.pos_a
  have hb := h.pos_b
  have hc := h.pos_c
  have hsum := cornerAngle_add_cornerAngle_add_cornerAngle h
  have cI := cos_cornerAngle h
  have cJ := cos_cornerAngle h.rotate
  have sI := sin_cornerAngle h
  have sJ := sin_cornerAngle h.rotate
  rw [triArea_rotate] at sJ
  constructor
  · simp only [J, K, polar, Real.cos_zero, Real.sin_zero, mul_one, mul_zero, sub_zero,
      Real.cos_pi_sub, Real.sin_pi_sub, Prod.mk.injEq]
    constructor
    · rw [cI, cJ]; field_simp; ring
    · rw [sI, sJ]; field_simp
  · have hβ : (π - cornerAngle b c a) + (π - cornerAngle c a b) = cornerAngle a b c + π := by
      linarith
    simp only [K, polar, hβ, Real.cos_add_pi, Real.sin_add_pi, Prod.mk.injEq]
    constructor <;> ring

/-- **The rotations `z`.** For a unit vector `w` (the direction of the edge `ij` at `i` in the
unnormalized frame) `z = e^{iβ}/w` is a rotation (`‖z‖ = 1`) taking `w` to `e^{iβ}`, the direction
of the edge in the face's frame.

[formula: iem-vector-interpolation] -/
theorem z_rotation {w : ℂ} (hw : ‖w‖ = 1) (β : ℝ) :
    ‖exp (β * I) / w‖ = 1 ∧ exp (β * I) / w * w = exp (β * I) := by
  have hw0 : w ≠ 0 := by rintro rfl; simp at hw
  refine ⟨?_, div_mul_cancel₀ _ hw0⟩
  rw [norm_div, hw, Complex.norm_exp_ofReal_mul_I, div_one]

end IntrinsicSimplification

end DDG
