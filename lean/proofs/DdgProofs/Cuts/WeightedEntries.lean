import DDG.Laplace.Weighted
import DdgProofs.Oracle.Mesh

/-!
# Cut flow: the oracle's weighted-Laplacian entries

`crates/ddg-variational-cuts/SPEC.md` §3.6, §3.10, D11. The oracle's `Ddg.triangleMean`,
`Ddg.weightedOffDiag` and `Ddg.weightedDiag` are the face mean and the entries of the weighted
cotan Laplacian `TriMesh.weightedLaplacian` of `DDG.Laplace.Weighted`:

* `Cuts.triangleMean_eq`, `Cuts.triangleMean_nonneg`, `Cuts.triangleMean_const`;
* `Cuts.weightedOffDiag_eq`, `Cuts.weightedDiag_eq`, `Cuts.weighted_quadratic_nonneg`: the
  oracle's entries are the entries of this matrix.

(Split off `DDG.Laplace.Weighted`; the declarations keep their names.)
-/

noncomputable section

namespace DDG

open Matrix

/-! ## The oracle's entries -/

namespace Cuts

open scoped DDG.Oracle

variable {faces : List Ddg.Face} {len : ℕ → ℕ → ℝ} {nv : ℕ}

/-- The oracle's `triangleMean` at `ℝ`. -/
theorem triangleMean_eq (x y z : ℝ) : Ddg.triangleMean x y z = (x + y + z) / 3 := by
  simp only [Ddg.triangleMean, Oracle.three_eq]

/-- The mean of nonnegative weights is nonnegative. -/
theorem triangleMean_nonneg {x y z : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) (hz : 0 ≤ z) :
    0 ≤ Ddg.triangleMean x y z := by
  rw [triangleMean_eq]
  positivity

/-- The mean of a constant is that constant. -/
theorem triangleMean_const (c : ℝ) : Ddg.triangleMean c c c = c := by
  rw [triangleMean_eq]
  ring

/-- **The oracle's off-diagonal entry is the entry of `TriMesh.weightedLaplacian`** of the mesh
it describes, with the face weights `w̄` read along the list of faces. -/
theorem weightedOffDiag_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv)
    (wbar : Ddg.Face → ℝ) (i j : Fin nv) :
    Ddg.weightedOffDiag faces len wbar i j =
      (Oracle.toTriMesh faces len nv h).weightedOffDiag (fun k => wbar (faces.get k)) i j := by
  unfold Ddg.weightedOffDiag TriMesh.weightedOffDiag
  rw [Oracle.sumList_map_eq]
  simp only [Oracle.cotOpposite_eq h, Oracle.two_eq]

/-- The oracle's entry of the weighted matrix: `weightedDiag` on the diagonal, `weightedOffDiag`
off it (the executable assembles the matrix from these two). -/
def weightedEntry (faces : List Ddg.Face) (len : ℕ → ℕ → ℝ) (nv : ℕ) (wbar : Ddg.Face → ℝ)
    (i j : ℕ) : ℝ :=
  if i = j then Ddg.weightedDiag nv (Ddg.weightedOffDiag faces len wbar) i
  else Ddg.weightedOffDiag faces len wbar i j

/-- **The oracle's weighted matrix is `TriMesh.weightedLaplacian`.** -/
theorem weightedEntry_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv)
    (wbar : Ddg.Face → ℝ) (i j : Fin nv) :
    weightedEntry faces len nv wbar i j =
      (Oracle.toTriMesh faces len nv h).weightedLaplacian (fun k => wbar (faces.get k)) i j := by
  unfold weightedEntry Ddg.weightedDiag Ddg.laplacianDiag TriMesh.weightedLaplacian lapOfOffDiag
  rw [Oracle.sumTo_eq]
  simp only [Matrix.of_apply, Fin.val_inj, Oracle.zero_eq]
  by_cases hij : i = j
  · subst hij
    simp only [ite_true]
    congr 1
    refine Finset.sum_congr rfl fun k _ => ?_
    by_cases hk : k = i
    · simp [hk]
    · simp only [hk, ite_false]
      exact weightedOffDiag_eq h wbar i k
  · simp only [hij, ite_false]
    exact weightedOffDiag_eq h wbar i j

/-- **The oracle's weighted matrix is positive semidefinite for nonnegative face weights**:
`Σ_i u_i Σ_j (L_w)_ij u_j ≥ 0` on every valid mesh.

[formula: cut-weighted-laplacian] -/
theorem weighted_quadratic_nonneg (hv : Oracle.MeshValid faces len nv) {wbar : Ddg.Face → ℝ}
    (hw : ∀ f ∈ faces, 0 ≤ wbar f) (u : ℕ → ℝ) :
    0 ≤ ∑ i : Fin nv, u i * ∑ j : Fin nv, weightedEntry faces len nv wbar i j * u j := by
  have hpsd := ((Oracle.toTriMesh faces len nv hv.lt).weightedLaplacian_posSemidef
    (Oracle.toTriMesh_valid hv) (w := fun k => wbar (faces.get k))
    (fun k => hw _ (List.get_mem faces k))).dotProduct_mulVec_nonneg (Oracle.res nv u)
  rw [star_trivial] at hpsd
  refine hpsd.trans_eq ?_
  simp only [dotProduct, Matrix.mulVec, Oracle.res]
  refine Finset.sum_congr rfl fun i _ => ?_
  congr 1
  exact Finset.sum_congr rfl fun j _ => by rw [weightedEntry_eq hv.lt wbar i j]

/-- **The rows of the oracle's weighted matrix sum to zero.**

[formula: cut-weighted-laplacian] -/
theorem weighted_row_sum (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv)
    (wbar : Ddg.Face → ℝ) (i : Fin nv) :
    ∑ j : Fin nv, weightedEntry faces len nv wbar i j = 0 := by
  have := congrFun ((Oracle.toTriMesh faces len nv h).weightedLaplacian_mulVec_one
    (fun k => wbar (faces.get k))) i
  simp only [Matrix.mulVec, dotProduct, mul_one, Pi.zero_apply] at this
  rw [← this]
  exact Finset.sum_congr rfl fun j _ => weightedEntry_eq h wbar i j

/-- **The oracle's weighted matrix is symmetric.**

[formula: cut-weighted-laplacian] -/
theorem weightedEntry_symm (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv)
    (wbar : Ddg.Face → ℝ) (i j : Fin nv) :
    weightedEntry faces len nv wbar i j = weightedEntry faces len nv wbar j i := by
  rw [weightedEntry_eq h, weightedEntry_eq h]
  exact congrFun (congrFun ((Oracle.toTriMesh faces len nv h).weightedLaplacian_isSymm
    (fun k => wbar (faces.get k))).eq.symm i) j

/-- **With the weight `1` the oracle's weighted off-diagonal entry is `laplacianOffDiag`.** -/
theorem weightedOffDiag_one (i j : ℕ) :
    Ddg.weightedOffDiag faces len (fun _ => (1 : ℝ)) i j = Ddg.laplacianOffDiag faces len i j := by
  unfold Ddg.weightedOffDiag Ddg.laplacianOffDiag
  simp only [one_mul]

end Cuts

end DDG
