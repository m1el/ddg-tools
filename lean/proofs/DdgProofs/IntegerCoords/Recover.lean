import DdgProofs.Geometry.Barycentric
import Mathlib.Algebra.BigOperators.Fin

/-!
# Integer coordinates: positions of an inserted vertex (tier A)

`crates/ddg-integer-coords/SPEC.md` §4.4 (face split, App. A), §9.2; the paper Sec. 3.3, Eq. 5,
App. A (Eq. 15–16).

* **`ic-recover-barycentric`** — App. A. A region `R` of a face `ijk` of `T¹` lies inside one face
  `abc` of `T⁰`; both are flat, so the barycentric coordinates of one point with respect to the two
  faces are related by an affine change: `u = A v` with `A` the `3 × 3` matrix whose column `i` holds
  the coordinates of corner `i` of `ijk` with respect to `abc` (`combo_change`: the plane point with
  coordinates `v` on `ijk` is the plane point with coordinates `A v` on `abc`, for any corner positions).
  If the corners `m = 0 … ρ−1` of `R` are known exactly on both sides (`u^(m) = A v^(m)`), then for
  **any** solution `ξ` of Eq. 15 (`Σ ξ_m v^(m) = v`, not only the minimum-norm one), Eq. 16
  (`u := Σ ξ_m u^(m)`) gives `A v`, the exact coordinates of the inserted point on `abc`
  (`recover_eq`), and `Σ ξ_m = 1` when the rows and `v` sum to one (`recover_sum_xi`). Nothing
  depends on `ρ` (3 to 6 in the paper) or on the solution chosen; how well the floating-point
  solve does is not claimed.
* **`ic-tangent-length`** — Eq. 5: a tangent vector with barycentric components `w`
  (`w_i + w_j + w_k = 0`) has `‖w‖² = −ℓ²_ij w_i w_j − ℓ²_jk w_j w_k − ℓ²_ki w_k w_i`
  (`tangent_length_sq`). This is the same statement as `cut-barycentric-length`
  (`DDG.Cuts.baryLength2_eq_norm_sq`), restated in the paper's form; the new edge lengths of a
  face split, `ℓ_im = ‖v − e_i‖`, are its value at `w = v − e_i` (`split_length`), which is what
  `it-point-distances` (`DDG.Intrinsic.cornerDistances_eq_baryLength`) proves for the
  intrinsic crate's `cornerDistances` on the unfolded triangle; here it is stated for any three
  points of a real inner product space.
-/

namespace DDG

namespace IntegerCoords

open Finset

/-! ## Recovering the position on `T⁰` (App. A) -/

/-- The plane point with barycentric coordinates `w` with respect to the corners `P`. -/
def combo (P : Fin 3 → ℝ × ℝ) (w : Fin 3 → ℝ) : ℝ × ℝ :=
  (∑ a, w a * (P a).1, ∑ a, w a * (P a).2)

/-- The coordinates on `abc` of the point with coordinates `v` on `ijk`, for the change matrix
`A` (column `i`: the coordinates of corner `i` of `ijk` on `abc`). -/
def changeCoords (A : Fin 3 → Fin 3 → ℝ) (v : Fin 3 → ℝ) : Fin 3 → ℝ := fun a => ∑ i, A a i * v i

/-- **The change of barycentric coordinates is linear.** If corner `i` of the triangle of `T¹` is
the point with coordinates `A · i` on the triangle `W` of `T⁰`, the point with coordinates `v` on the
triangle of `T¹` is the point with coordinates `A v` on `W`. [formula: ic-recover-barycentric] -/
theorem combo_change (W : Fin 3 → ℝ × ℝ) (A : Fin 3 → Fin 3 → ℝ) (v : Fin 3 → ℝ) :
    combo (fun i => combo W (fun a => A a i)) v = combo W (changeCoords A v) := by
  simp only [combo, changeCoords, Fin.sum_univ_three]
  ext <;> simp only <;> ring

/-- **Eq. 15–16 reproduce the point exactly, for any solution `ξ`.** The corners `m` of the region
have coordinates `v m` on the face of `T¹` and `u m = A (v m)` on the face of `T⁰`; if `ξ` solves
Eq. 15, `Σ_m ξ_m v^(m) = v`, then Eq. 16, `Σ_m ξ_m u^(m)`, is `A v`.
[formula: ic-recover-barycentric] -/
theorem recover_eq {ρ : ℕ} (A : Fin 3 → Fin 3 → ℝ) (vs us : Fin ρ → Fin 3 → ℝ)
    (hus : ∀ m, us m = changeCoords A (vs m)) (ξ : Fin ρ → ℝ) (v : Fin 3 → ℝ)
    (h15 : ∀ i, ∑ m, ξ m * vs m i = v i) :
    (fun a => ∑ m, ξ m * us m a) = changeCoords A v := by
  funext a
  simp only [hus, changeCoords, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [← h15 i, Finset.mul_sum]
  exact Finset.sum_congr rfl (fun m _ => by ring)

/-- With rows and `v` summing to one, any solution of Eq. 15 sums to one (so Eq. 16 is an affine
combination). [formula: ic-recover-barycentric] -/
theorem recover_sum_xi {ρ : ℕ} (vs : Fin ρ → Fin 3 → ℝ) (hrows : ∀ m, ∑ i, vs m i = 1)
    (ξ : Fin ρ → ℝ) (v : Fin 3 → ℝ) (hv : ∑ i, v i = 1) (h15 : ∀ i, ∑ m, ξ m * vs m i = v i) :
    ∑ m, ξ m = 1 := by
  rw [← hv]
  simp only [← h15]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  rw [← Finset.mul_sum, hrows m, mul_one]

/-- The recovered coordinates sum to one when the columns of `A` do (the corners of `ijk` are
points of the plane of `abc`) and `v` sums to one. [formula: ic-recover-barycentric] -/
theorem changeCoords_sum (A : Fin 3 → Fin 3 → ℝ) (hA : ∀ i, ∑ a, A a i = 1) (v : Fin 3 → ℝ)
    (hv : ∑ i, v i = 1) : ∑ a, changeCoords A v a = 1 := by
  simp only [changeCoords]
  rw [Finset.sum_comm]
  simp only [← Finset.sum_mul, hA, one_mul, hv]

/-! ## Lengths in barycentric coordinates (Eq. 5) -/

section Lengths

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- **Eq. 5**: a vector with barycentric components `w_i + w_j + w_k = 0` has squared length
`−ℓ²_ij w_i w_j − ℓ²_jk w_j w_k − ℓ²_ki w_k w_i`. The restatement of `cut-barycentric-length`.
[formula: ic-tangent-length] -/
theorem tangent_length_sq (Pi Pj Pk : E) {wi wj wk : ℝ} (h : wi + wj + wk = 0) :
    ‖wi • Pi + wj • Pj + wk • Pk‖ ^ 2 =
      -(‖Pi - Pj‖ ^ 2 * wi * wj) - ‖Pj - Pk‖ ^ 2 * wj * wk - ‖Pk - Pi‖ ^ 2 * wk * wi := by
  rw [← Cuts.baryLength2_eq_norm_sq Pi Pj Pk h]
  unfold baryLength2
  ring

/-- **The lengths of a face split** (paper Sec. 3.3): the new vertex `m` with coordinates `v`
(`Σ v = 1`) is at distance `√(Eq. 5 at w = v − e_i)` from corner `i`, i.e.
`ℓ_im = √(−ℓ²_ij (v_i − 1) v_j − ℓ²_jk v_j v_k − ℓ²_ki v_k (v_i − 1))`.
[formula: ic-tangent-length] -/
theorem split_length (Pi Pj Pk : E) {vi vj vk : ℝ} (h : vi + vj + vk = 1) :
    dist (vi • Pi + vj • Pj + vk • Pk) Pi =
      Real.sqrt (-(‖Pi - Pj‖ ^ 2 * (vi - 1) * vj) - ‖Pj - Pk‖ ^ 2 * vj * vk -
        ‖Pk - Pi‖ ^ 2 * vk * (vi - 1)) := by
  have e : vi • Pi + vj • Pj + vk • Pk - Pi = (vi - 1) • Pi + vj • Pj + vk • Pk := by
    simp only [sub_smul, one_smul]; abel
  rw [dist_eq_norm, e, ← tangent_length_sq Pi Pj Pk (by linarith), Real.sqrt_sq (norm_nonneg _)]

end Lengths

end IntegerCoords

end DDG
