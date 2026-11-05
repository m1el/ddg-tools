import DdgProofs.Triangle

/-!
# Barycentric combinations: squared lengths and Heron's polynomial

* **`cut-barycentric-length`** — for a finite family of points `P_i` of a real inner product
  space and coefficients `d_i` with `Σ d_i = 0`,
  `‖Σ d_i P_i‖² = −½ Σ_{i,j} d_i d_j ‖P_i − P_j‖²` (`norm_sq_sum_smul_of_sum_eq_zero`): the length
  of a barycentric difference is a function of the pairwise distances.
* Heron's formula in the squared side lengths, `triArea a b c = √(heronPoly a² b² c²)/4`
  (`heronPoly`, `triArea_eq_sqrt_heronPoly`).

The oracle-level formulas built on these (`Ddg.baryLength`, `Ddg.subTriangleArea`, `Ddg.baryInner`)
are still in `DdgProofs/Cuts/Crossing.lean`; step S3 gives them native definitions here.
(Library step S2b: moved here from `DdgProofs/Cuts/Crossing.lean`; the names keep the
`DDG.Cuts` namespace until the renames of step S4.)
-/

noncomputable section

namespace DDG

namespace Cuts

open Real

section Barycentric

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- **Squared length of a combination with coefficient sum zero**, for any finite family of
points of a real inner product space:
`‖Σ d_i P_i‖² = −½ Σ_{i,j} d_i d_j ‖P_i − P_j‖² = −Σ_{i<j} d_i d_j ‖P_i − P_j‖²`.

[formula: cut-barycentric-length] -/
theorem norm_sq_sum_smul_of_sum_eq_zero {ι : Type*} [Fintype ι] (P : ι → E) (d : ι → ℝ)
    (hd : ∑ i, d i = 0) :
    ‖∑ i, d i • P i‖ ^ 2 = -(1 / 2) * ∑ i, ∑ j, d i * d j * ‖P i - P j‖ ^ 2 := by
  have h1 : ∀ i j, d i * d j * ‖P i - P j‖ ^ 2 =
      d j * (d i * ‖P i‖ ^ 2) - 2 * (d i * d j * inner ℝ (P i) (P j)) +
        d i * (d j * ‖P j‖ ^ 2) := by
    intro i j
    rw [norm_sub_sq_real]
    ring
  have h2 : ‖∑ i, d i • P i‖ ^ 2 = ∑ i, ∑ j, d i * d j * inner ℝ (P i) (P j) := by
    rw [← real_inner_self_eq_norm_sq, sum_inner]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [inner_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [real_inner_smul_left, real_inner_smul_right]
    ring
  have e1 : ∑ i, ∑ j, d j * (d i * ‖P i‖ ^ 2) = 0 := by
    refine Finset.sum_eq_zero fun i _ => ?_
    rw [← Finset.sum_mul, hd, zero_mul]
  have e3 : ∑ i, ∑ j, d i * (d j * ‖P j‖ ^ 2) = 0 := by
    rw [Finset.sum_comm]
    exact e1
  have h3 : ∑ i, ∑ j, d i * d j * ‖P i - P j‖ ^ 2 =
      ∑ i, ∑ j, d j * (d i * ‖P i‖ ^ 2) - 2 * ∑ i, ∑ j, d i * d j * inner ℝ (P i) (P j) +
        ∑ i, ∑ j, d i * (d j * ‖P j‖ ^ 2) := by
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun j _ => h1 i j
  rw [h2, h3, e1, e3]
  ring

end Barycentric

/-- Sixteen times the squared area of a triangle as a function of its squared side lengths:
`2xy + 2yz + 2zx − x² − y² − z²` (Heron). -/
def heronPoly (x y z : ℝ) : ℝ := 2 * x * y + 2 * y * z + 2 * z * x - x ^ 2 - y ^ 2 - z ^ 2

/-- Heron's formula in terms of the squared sides. -/
theorem triArea_eq_sqrt_heronPoly (a b c : ℝ) :
    triArea a b c = √(heronPoly (a ^ 2) (b ^ 2) (c ^ 2)) / 4 := by
  unfold triArea heronPoly
  congr 2
  ring

end Cuts

end DDG
