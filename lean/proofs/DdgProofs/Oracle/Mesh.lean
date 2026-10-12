import DdgProofs.Oracle.Problem
import DdgProofs.GaussBonnet

/-!
# The mesh formulas of the oracle: cotan Laplacian, curvature, lumped mass

The oracle (`lean/oracle`, `DdgOracle/Formulas/Mesh.lean`) gives a mesh as a list of faces with
vertex indices in `Nat`, an edge-length function `len : Nat → Nat → α` and a vertex count `nv`.
This file identifies its `laplacian`, `angleDefect`, `geodesicCurvature` and `lumpedMass` at the
real numbers with the corresponding quantities of the `TriMesh` it describes
(`Oracle.toTriMesh`), and restates the theorems of `DdgProofs.Cotan` and `DdgProofs.GaussBonnet`
about the oracle's definitions.
-/

noncomputable section

namespace DdgProofs

namespace Oracle

open Matrix

/-- The oracle's `sumList` of a mapped list is the sum over the positions of the list. -/
theorem sumList_map_eq {β : Type*} (l : List β) (g : β → ℝ) :
    Ddg.sumList (l.map g) = ∑ k : Fin l.length, g (l.get k) := by
  induction l with
  | nil => simp [Ddg.sumList]
  | cons x xs ih =>
    simp only [List.map_cons, Ddg.sumList, List.length_cons]
    rw [Fin.sum_univ_succ, ih]
    simp

/-- Validity of an oracle mesh: vertex indices in range, three distinct vertices per face, and
strict triangle inequalities for the three edge lengths of every face. -/
structure MeshValid (faces : List Ddg.Face) (len : ℕ → ℕ → ℝ) (nv : ℕ) : Prop where
  lt : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv
  ne : ∀ f ∈ faces, f.a ≠ f.b ∧ f.b ≠ f.c ∧ f.c ≠ f.a
  tri : ∀ f ∈ faces, IsTriangle (len f.b f.c) (len f.c f.a) (len f.a f.b)

/-- The oracle's mesh data as a `TriMesh` on the vertices `Fin nv` and the faces
`Fin faces.length`. -/
def toTriMesh (faces : List Ddg.Face) (len : ℕ → ℕ → ℝ) (nv : ℕ)
    (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) : TriMesh (Fin nv) (Fin faces.length) where
  a := fun k => ⟨(faces.get k).a, (h _ (List.get_mem _ _)).1⟩
  b := fun k => ⟨(faces.get k).b, (h _ (List.get_mem _ _)).2.1⟩
  c := fun k => ⟨(faces.get k).c, (h _ (List.get_mem _ _)).2.2⟩
  len := fun i j => len i j

variable {faces : List Ddg.Face} {len : ℕ → ℕ → ℝ} {nv : ℕ}

theorem toTriMesh_valid (hv : MeshValid faces len nv) : (toTriMesh faces len nv hv.lt).Valid where
  ab := fun _ h => (hv.ne _ (List.get_mem _ _)).1 (congrArg Fin.val h)
  bc := fun _ h => (hv.ne _ (List.get_mem _ _)).2.1 (congrArg Fin.val h)
  ca := fun _ h => (hv.ne _ (List.get_mem _ _)).2.2 (congrArg Fin.val h)
  tri := fun _ => hv.tri _ (List.get_mem _ _)

/-! ## The cotan Laplacian -/

theorem cotOpposite_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv)
    (k : Fin faces.length) (i j : Fin nv) :
    (faces.get k).cotOpposite len i j = (toTriMesh faces len nv h).cotOpp k i j := by
  have hA : (faces.get k).EdgeAB i j ↔
      ((i = (toTriMesh faces len nv h).a k ∧ j = (toTriMesh faces len nv h).b k) ∨
        (i = (toTriMesh faces len nv h).b k ∧ j = (toTriMesh faces len nv h).a k)) := by
    simp [Ddg.Face.EdgeAB, toTriMesh, Fin.ext_iff]
  have hB : (faces.get k).EdgeBC i j ↔
      ((i = (toTriMesh faces len nv h).b k ∧ j = (toTriMesh faces len nv h).c k) ∨
        (i = (toTriMesh faces len nv h).c k ∧ j = (toTriMesh faces len nv h).b k)) := by
    simp [Ddg.Face.EdgeBC, toTriMesh, Fin.ext_iff]
  have hC : (faces.get k).EdgeCA i j ↔
      ((i = (toTriMesh faces len nv h).c k ∧ j = (toTriMesh faces len nv h).a k) ∨
        (i = (toTriMesh faces len nv h).a k ∧ j = (toTriMesh faces len nv h).c k)) := by
    simp [Ddg.Face.EdgeCA, toTriMesh, Fin.ext_iff]
  unfold Ddg.Face.cotOpposite TriMesh.cotOpp
  simp only [hA, hB, hC, Ddg.Face.cotA, Ddg.Face.cotB, Ddg.Face.cotC, cornerCot_eq, zero_eq,
    TriMesh.cotA, TriMesh.cotB, TriMesh.cotC, TriMesh.la, TriMesh.lb, TriMesh.lc]
  rfl

theorem laplacianOffDiag_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) (i j : Fin nv) :
    Ddg.laplacianOffDiag faces len i j = (toTriMesh faces len nv h).offDiag i j := by
  unfold Ddg.laplacianOffDiag TriMesh.offDiag
  rw [sumList_map_eq]
  simp only [cotOpposite_eq h, two_eq]

/-- **The oracle's cotan Laplacian is `TriMesh.laplacian`** of the mesh it describes. -/
theorem laplacian_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) (i j : Fin nv) :
    Ddg.laplacian faces len nv i j = (toTriMesh faces len nv h).laplacian i j := by
  unfold Ddg.laplacian Ddg.laplacianDiag TriMesh.laplacian
  rw [sumTo_eq]
  simp only [Matrix.of_apply, Fin.val_inj, zero_eq]
  by_cases hij : i = j
  · subst hij
    simp only [ite_true]
    congr 1
    refine Finset.sum_congr rfl fun k _ => ?_
    by_cases hk : k = i
    · simp [hk]
    · simp only [hk, ite_false]
      exact laplacianOffDiag_eq h i k
  · simp only [hij, ite_false]
    exact laplacianOffDiag_eq h i j

/-- **The oracle's cotan Laplacian is symmetric.**

[formula: cotan-laplacian] -/
theorem laplacian_symm (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) {i j : ℕ}
    (hi : i < nv) (hj : j < nv) :
    Ddg.laplacian faces len nv i j = Ddg.laplacian faces len nv j i := by
  have e1 := laplacian_eq (len := len) h ⟨i, hi⟩ ⟨j, hj⟩
  have e2 := laplacian_eq (len := len) h ⟨j, hj⟩ ⟨i, hi⟩
  simp only at e1 e2
  rw [e1, e2]
  exact congrFun (congrFun (toTriMesh faces len nv h).laplacian_isSymm.eq.symm ⟨i, hi⟩) ⟨j, hj⟩

/-- **The rows of the oracle's cotan Laplacian sum to zero.**

[formula: cotan-laplacian] -/
theorem laplacian_row_sum (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) {i : ℕ}
    (hi : i < nv) : Ddg.sumTo nv (fun j => Ddg.laplacian faces len nv i j) = 0 := by
  rw [sumTo_eq]
  have := congrFun (toTriMesh faces len nv h).laplacian_mulVec_one ⟨i, hi⟩
  simp only [Matrix.mulVec, dotProduct, mul_one, Pi.zero_apply] at this
  rw [← this]
  exact Finset.sum_congr rfl fun j _ => laplacian_eq h ⟨i, hi⟩ j

/-- The oracle's quadratic form `Σ_i u_i Σ_j L_ij u_j` is `uᵀ L u` of the `TriMesh`. -/
theorem quadratic_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) (u : ℕ → ℝ) :
    Ddg.sumTo nv (fun i => u i * Ddg.sumTo nv (fun j => Ddg.laplacian faces len nv i j * u j)) =
      res nv u ⬝ᵥ ((toTriMesh faces len nv h).laplacian *ᵥ res nv u) := by
  rw [sumTo_eq]
  simp only [dotProduct, Matrix.mulVec, res]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [sumTo_eq]
  congr 1
  exact Finset.sum_congr rfl fun j _ => by rw [laplacian_eq h i j]

/-- **The oracle's cotan Laplacian is positive semidefinite** on every valid mesh (no Delaunay
assumption): `Σ_i u_i Σ_j L_ij u_j ≥ 0`.

[formula: cotan-laplacian-psd] -/
theorem laplacian_quadratic_nonneg (hv : MeshValid faces len nv) (u : ℕ → ℝ) :
    0 ≤ Ddg.sumTo nv
      (fun i => u i * Ddg.sumTo nv (fun j => Ddg.laplacian faces len nv i j * u j)) := by
  rw [quadratic_eq hv.lt]
  have := ((toTriMesh faces len nv hv.lt).laplacian_posSemidef
    (toTriMesh_valid hv)).dotProduct_mulVec_nonneg (res nv u)
  rwa [star_trivial] at this

/-! ## Curvature -/

theorem angleAt_eq (hv : MeshValid faces len nv) (k : Fin faces.length) (v : Fin nv) :
    (faces.get k).angleAt len v = (toTriMesh faces len nv hv.lt).angleAt k v := by
  have ht := hv.tri _ (List.get_mem faces k)
  unfold Ddg.Face.angleAt TriMesh.angleAt
  simp only [Ddg.Face.angleA, Ddg.Face.angleB, Ddg.Face.angleC, cornerAngle_eq ht,
    cornerAngle_eq ht.rotate, cornerAngle_eq ht.rotate.rotate, zero_eq]
  simp only [toTriMesh, Fin.ext_iff]
  rfl

theorem angleSum_eq (hv : MeshValid faces len nv) (v : Fin nv) :
    Ddg.angleSum faces len v = (toTriMesh faces len nv hv.lt).angleSum v := by
  unfold Ddg.angleSum TriMesh.angleSum
  rw [sumList_map_eq]
  exact Finset.sum_congr rfl fun k _ => angleAt_eq hv k v

/-- The oracle's angle defect `2π − Σθ` is `TriMesh.angleDefect`.

[formula: angle-defect] -/
theorem angleDefect_eq (hv : MeshValid faces len nv) (v : Fin nv) :
    Ddg.angleDefect faces len v = (toTriMesh faces len nv hv.lt).angleDefect v := by
  unfold Ddg.angleDefect TriMesh.angleDefect
  rw [angleSum_eq hv]
  simp

/-- The oracle's boundary curvature `π − Σθ` is `TriMesh.geodesicCurvature`.

[formula: boundary-geodesic-curvature] -/
theorem geodesicCurvature_eq (hv : MeshValid faces len nv) (v : Fin nv) :
    Ddg.geodesicCurvature faces len v = (toTriMesh faces len nv hv.lt).geodesicCurvature v := by
  unfold Ddg.geodesicCurvature TriMesh.geodesicCurvature
  rw [angleSum_eq hv]
  simp

/-- **Discrete Gauss–Bonnet for the oracle's curvatures.** For a valid mesh whose `TriMesh` is a
combinatorial manifold, the oracle's angle defects of the interior vertices and geodesic
curvatures of the boundary vertices add up to `2π·χ`.

[formula: gauss-bonnet] -/
theorem gauss_bonnet (hv : MeshValid faces len nv)
    (hman : (toTriMesh faces len nv hv.lt).Manifold) :
    ∑ v ∈ (toTriMesh faces len nv hv.lt).boundaryVerticesᶜ, Ddg.angleDefect faces len v +
      ∑ v ∈ (toTriMesh faces len nv hv.lt).boundaryVertices, Ddg.geodesicCurvature faces len v =
      2 * Real.pi * ((toTriMesh faces len nv hv.lt).eulerChar : ℝ) := by
  rw [← (toTriMesh faces len nv hv.lt).gauss_bonnet (toTriMesh_valid hv) hman]
  congr 1
  · exact Finset.sum_congr rfl fun v _ => angleDefect_eq hv v
  · exact Finset.sum_congr rfl fun v _ => geodesicCurvature_eq hv v

/-! ## Lumped mass -/

theorem area_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) (k : Fin faces.length) :
    (faces.get k).area len = (toTriMesh faces len nv h).area k := by
  unfold Ddg.Face.area TriMesh.area
  rw [triArea_eq]
  exact (triArea_rotate (len (faces.get k).c (faces.get k).a)
    (len (faces.get k).a (faces.get k).b) (len (faces.get k).b (faces.get k).c)).trans
    (triArea_rotate _ _ _)

theorem areaAt_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) (v : Fin nv) :
    Ddg.areaAt faces len v = (toTriMesh faces len nv h).vertexArea v := by
  unfold Ddg.areaAt TriMesh.vertexArea
  rw [sumList_map_eq]
  refine Finset.sum_congr rfl fun k _ => ?_
  unfold Ddg.Face.areaAt
  simp only [Ddg.Face.Has, area_eq h, zero_eq]
  simp only [toTriMesh, Fin.ext_iff]

theorem totalArea_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) :
    Ddg.totalArea faces len = (toTriMesh faces len nv h).totalArea := by
  unfold Ddg.totalArea TriMesh.totalArea
  rw [sumList_map_eq]
  exact Finset.sum_congr rfl fun k _ => area_eq h k

/-- The oracle's `lumpedMass` is the lumped mass of the `TriMesh` divided by its total area. -/
theorem lumpedMass_eq (h : ∀ f ∈ faces, f.a < nv ∧ f.b < nv ∧ f.c < nv) (v : Fin nv) :
    Ddg.lumpedMass faces len v =
      (toTriMesh faces len nv h).mass v / (toTriMesh faces len nv h).totalArea := by
  unfold Ddg.lumpedMass Ddg.massOf TriMesh.mass
  rw [areaAt_eq h, totalArea_eq h]

/-- **The oracle's normalized lumped masses sum to one** on a valid mesh with at least one face.

[formula: lumped-mass] -/
theorem sum_lumpedMass (hv : MeshValid faces len nv) (hne : faces ≠ []) :
    Ddg.sumTo nv (fun v => Ddg.lumpedMass faces len v) = 1 := by
  rw [sumTo_eq]
  have hM := toTriMesh_valid hv
  have : Nonempty (Fin faces.length) := ⟨⟨0, List.length_pos_of_ne_nil hne⟩⟩
  have hA : (toTriMesh faces len nv hv.lt).totalArea ≠ 0 :=
    ((toTriMesh faces len nv hv.lt).totalArea_pos hM).ne'
  rw [← (toTriMesh faces len nv hv.lt).sum_mass_div_totalArea hM hA]
  exact Finset.sum_congr rfl fun v _ => lumpedMass_eq hv.lt v

end Oracle

end DdgProofs
