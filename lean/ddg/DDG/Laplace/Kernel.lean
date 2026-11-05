import Mathlib

/-!
# A linear system with zero column sums: one equation is implied by the others

For a matrix `A` over a finite index type with `1ᵀA = 0` (for instance the cotan Laplacian,
`TriMesh.laplacian_vecMul_one`):

* `Solver.mulVec_eq_of_sum_eq`: if `Σ b = 0` and `(A u)_i = b_i` for every `i ≠ i₀`, then also
  `(A u)_{i₀} = b_{i₀}` (the dropped equation of a pinned solve on a closed mesh);
* `Solver.sum_eq_zero_of_mulVec_eq`: conversely, `A u = b` forces `Σ b = 0`.

(Moved out of `DdgProofs.OptimalCones.Solver`; the declarations keep their names.)
-/

noncomputable section

namespace DDG

open Matrix

namespace Solver

variable {ι : Type*} [Fintype ι]

/-- **An equation of a system with zero column sums is implied by the others.** If `1ᵀA = 0`,
`Σ b = 0`, and `(A u)_i = b_i` for every `i ≠ i₀`, then also `(A u)_{i₀} = b_{i₀}`.

This is why, on a closed mesh, the code may drop one equation: the state equation of the vertex
whose adjoint value it pins (`newton.rs`, when no vertex is active, after the cone angles were made
to sum to `Σ Ω`), and the first row of the Laplacian in `evaluate.rs` (`LaplaceSolver`).

[formula: gauss-bonnet-adjoint-shift] -/
theorem mulVec_eq_of_sum_eq {A : Matrix ι ι ℝ} (hcol : (fun _ => (1 : ℝ)) ᵥ* A = 0)
    {u b : ι → ℝ} (hb : ∑ i, b i = 0) (i₀ : ι) (h : ∀ i, i ≠ i₀ → (A *ᵥ u) i = b i) :
    (A *ᵥ u) i₀ = b i₀ := by
  classical
  have h1 : ∑ i, (A *ᵥ u) i = 0 := by
    have : (fun _ => (1 : ℝ)) ⬝ᵥ (A *ᵥ u) = 0 := by
      rw [Matrix.dotProduct_mulVec, hcol]; simp
    simpa [dotProduct] using this
  have h2 : ∑ i, ((A *ᵥ u) i - b i) = (A *ᵥ u) i₀ - b i₀ := by
    rw [Finset.sum_eq_single i₀]
    · intro i _ hi; rw [h i hi]; ring
    · intro hi; exact absurd (Finset.mem_univ i₀) hi
  rw [Finset.sum_sub_distrib, h1, hb] at h2
  linarith

/-- Conversely, a system with zero column sums is solvable only if its right-hand side sums to
zero. -/
theorem sum_eq_zero_of_mulVec_eq {A : Matrix ι ι ℝ} (hcol : (fun _ => (1 : ℝ)) ᵥ* A = 0)
    {u b : ι → ℝ} (h : A *ᵥ u = b) : ∑ i, b i = 0 := by
  have : (fun _ => (1 : ℝ)) ⬝ᵥ (A *ᵥ u) = 0 := by
    rw [Matrix.dotProduct_mulVec, hcol]; simp
  rw [h] at this
  simpa [dotProduct] using this

end Solver

end DDG
