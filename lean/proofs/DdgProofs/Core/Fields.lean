import DDG.Mesh.GaussBonnet
import DDG.Laplace.Cotan

/-!
# Curvature fields and Laplacian entries in the form the Rust code evaluates them

`crates/ddg-core` stores the curvatures as per-vertex fields that vanish on the other class of
vertices and assembles the cotan Laplacian entrywise. This file proves that those arrangements
are the quantities of `DDG.GaussBonnet` and `DDG.Cotan`:

* `TriMesh.gauss_bonnet_fields` — `angle_defects`, `boundary_geodesic_curvatures`;
* `TriMesh.laplacian_apply_of_ne`, `TriMesh.laplacian_apply_self` — off-diagonal and diagonal
  entries of the cotan Laplacian.

(The Kahan and `atan2` forms are in `DdgProofs.Evaluated`.)

Nothing is claimed about floating-point rounding.
-/

noncomputable section

namespace DDG

open Real

/-! ## Curvature fields -/

namespace TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype V] [Fintype F] (M : TriMesh V F)

/-- The angle defect as a per-vertex field (`angle_defects` in `crates/ddg-core`):
`Ω_v = 2π − Σθ` at interior vertices and `0` at boundary vertices. -/
def angleDefectField (v : V) : ℝ := if v ∈ M.boundaryVertices then 0 else M.angleDefect v

/-- The boundary curvature as a per-vertex field (`boundary_geodesic_curvatures`):
`k_v = π − Σθ` at boundary vertices and `0` at interior vertices. -/
def geodesicCurvatureField (v : V) : ℝ :=
  if v ∈ M.boundaryVertices then M.geodesicCurvature v else 0

/-- **Discrete Gauss–Bonnet for the two per-vertex fields**: the sum of all entries of the
angle-defect field plus the sum of all entries of the boundary-curvature field is `2π·χ`.

[formula: gauss-bonnet] [formula: angle-defect] [formula: boundary-geodesic-curvature] -/
theorem gauss_bonnet_fields (hM : M.Valid) (hman : M.Manifold) :
    ∑ v, M.angleDefectField v + ∑ v, M.geodesicCurvatureField v =
      2 * π * (M.eulerChar : ℝ) := by
  rw [← M.gauss_bonnet hM hman]
  congr 1
  · unfold angleDefectField
    rw [← Finset.sum_filter_add_sum_filter_not Finset.univ (fun v => v ∈ M.boundaryVertices)]
    have h1 : ∑ v ∈ Finset.univ.filter (fun v => v ∈ M.boundaryVertices),
        (if v ∈ M.boundaryVertices then 0 else M.angleDefect v) = 0 :=
      Finset.sum_eq_zero fun v hv => by simp [(Finset.mem_filter.mp hv).2]
    have h2 : Finset.univ.filter (fun v => v ∉ M.boundaryVertices) = M.boundaryVerticesᶜ := by
      ext v; simp
    rw [h1, zero_add, h2]
    exact Finset.sum_congr rfl fun v hv => by simp [Finset.mem_compl.mp hv]
  · unfold geodesicCurvatureField
    rw [← Finset.sum_filter_add_sum_filter_not Finset.univ (fun v => v ∈ M.boundaryVertices)]
    have h1 : ∑ v ∈ Finset.univ.filter (fun v => v ∉ M.boundaryVertices),
        (if v ∈ M.boundaryVertices then M.geodesicCurvature v else 0) = 0 :=
      Finset.sum_eq_zero fun v hv => by simp [(Finset.mem_filter.mp hv).2]
    have h2 : Finset.univ.filter (fun v => v ∈ M.boundaryVertices) = M.boundaryVertices := by
      ext v; simp
    rw [h1, add_zero, h2]
    exact Finset.sum_congr rfl fun v hv => by simp [hv]

/-! ## Entries of the cotan Laplacian -/

/-- **Off-diagonal entries**: for `i ≠ j`, `L_ij = −½ Σ_f (cot of the corner of f opposite the
edge ij)`; a face that does not contain the edge contributes `0`, so on a manifold mesh the sum
has one or two nonzero terms.

[formula: cotan-laplacian] -/
theorem laplacian_apply_of_ne {i j : V} (h : i ≠ j) :
    M.laplacian i j = -(1 / 2) * ∑ f, M.cotOpp f i j := by
  simp only [laplacian, offDiag, Matrix.of_apply, h, ite_false]
  ring

/-- **Diagonal entries**: `L_ii = −Σ_{j ≠ i} L_ij`.

[formula: cotan-laplacian] -/
theorem laplacian_apply_self (i : V) :
    M.laplacian i i = -∑ j ∈ Finset.univ.erase i, M.laplacian i j := by
  have h1 : ∑ j ∈ Finset.univ.erase i, M.laplacian i j =
      ∑ j ∈ Finset.univ.erase i, (if j = i then 0 else M.offDiag i j) := by
    refine Finset.sum_congr rfl fun j hj => ?_
    have hne : j ≠ i := Finset.ne_of_mem_erase hj
    simp [laplacian, hne, Ne.symm hne]
  rw [h1, Finset.sum_erase _ (by simp)]
  simp [laplacian]

end TriMesh

end DDG
