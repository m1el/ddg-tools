import DdgProofs.Oracle.Mesh

/-!
# Cut flow: the weighted cotan Laplacian

`crates/ddg-variational-cuts/SPEC.md` §3.6, §3.10, D11. For a per-face weight `w̄_f`,

`L_w = Σ_f w̄_f · (cotan Laplacian of the single face f)`,

i.e. `(L_w)_ij = −½ Σ_f w̄_f cot(corner of f opposite ij)` for `i ≠ j` and zero row sums. With
`w̄_f = (w_i + w_j + w_k)/3` the mean of a per-vertex weight (`Ddg.triangleMean`) this is the
matrix `L_a` of the scaled Dirichlet energy and the matrix `L_w` of the length step.

**`cut-weighted-laplacian`**:

* `TriMesh.weightedLaplacian_isSymm`, `weightedLaplacian_mulVec_one`: symmetric, zero row sums;
* `TriMesh.dotProduct_weightedLaplacian_mulVec`: `xᵀ L_w x = Σ_f w̄_f · (cotan energy of x on f)`,
  the cotan energy being the Dirichlet energy of the linear interpolant (`cotan_formula`);
* `TriMesh.weightedLaplacian_posSemidef`: positive semidefinite for `w̄ ≥ 0` on every valid mesh
  (also with obtuse triangles); `weightedLaplacian_faceMean_posSemidef` for vertex weights `≥ 0`;
* `TriMesh.weightedLaplacian_const`: `w̄ ≡ c` gives `c·L`; `weightedLaplacian_faceMean_const`:
  a constant vertex weight `c` gives `c·L`;
* `edgeRule_const`, `edgeRuleEnergy_neg_example`: the paper's edge rule (cotangent of edge `ij`
  times `(a_i + a_j)/2`, Sec. 5.3) also gives `c·L` for a constant weight, but its quadratic
  form is negative for a nonnegative weight on an obtuse triangle (D11);
* `Cuts.weightedOffDiag_eq`, `Cuts.weightedDiag_eq`, `Cuts.weighted_quadratic_nonneg`: the
  oracle's entries `Ddg.weightedOffDiag`, `Ddg.weightedDiag` are the entries of this matrix.

What is not formalized: that `ā_t · (cotan energy)` is the integral `∫_t a |∇f|²` for piecewise
linear `a` and `f` (the mean value of a linear function over a triangle; the integral over a
triangle is not formalized anywhere in this package, see `README.md`).
-/

noncomputable section

namespace DDG

open Matrix

/-! ## A matrix with given off-diagonal entries and zero row sums -/

section OffDiag

variable {V : Type*} [DecidableEq V] [Fintype V]

/-- The matrix with off-diagonal entries `o i j` and zero row sums. -/
def lapOfOffDiag (o : V → V → ℝ) : Matrix V V ℝ :=
  Matrix.of fun i j => if i = j then -∑ k, (if k = i then 0 else o i k) else o i j

theorem lapOfOffDiag_isSymm {o : V → V → ℝ} (ho : ∀ i j, o i j = o j i) :
    (lapOfOffDiag o).IsSymm := by
  ext i j
  simp only [lapOfOffDiag, Matrix.transpose_apply, Matrix.of_apply]
  by_cases h : i = j
  · subst h; simp
  · rw [ite_eq_right h, ite_eq_right (Ne.symm h)]
    exact ho j i

/-- Row `i` applied to `u`: `Σ_j (−o_ij)(u_i − u_j)`. -/
theorem lapOfOffDiag_mulVec_apply (o : V → V → ℝ) (u : V → ℝ) (i : V) :
    (lapOfOffDiag o *ᵥ u) i = ∑ j, -o i j * (u i - u j) := by
  set p : V → ℝ := fun j => if j = i then 0 else o i j with hp
  have h1 : (lapOfOffDiag o *ᵥ u) i = ∑ j, ((if i = j then -∑ k, p k else 0) + p j) * u j := by
    simp only [Matrix.mulVec, dotProduct, lapOfOffDiag, Matrix.of_apply]
    refine Finset.sum_congr rfl fun j _ => ?_
    by_cases h : i = j
    · subst h; simp [hp]
    · simp [hp, h, Ne.symm h]
  have h2 : ∑ j, -o i j * (u i - u j) = ∑ j, -p j * (u i - u j) := by
    refine Finset.sum_congr rfl fun j _ => ?_
    by_cases h : j = i
    · subst h; simp
    · simp [hp, h]
  rw [h1, h2]
  simp only [add_mul, Finset.sum_add_distrib, ite_mul, zero_mul, Finset.sum_ite_eq,
    Finset.mem_univ, ite_true, neg_mul, mul_sub, Finset.sum_neg_distrib, Finset.sum_sub_distrib,
    ← Finset.sum_mul]
  ring

theorem lapOfOffDiag_mulVec_one (o : V → V → ℝ) :
    lapOfOffDiag o *ᵥ (fun _ => (1 : ℝ)) = 0 := by
  ext i
  rw [lapOfOffDiag_mulVec_apply]
  simp

/-- The quadratic form as a sum over ordered pairs. -/
theorem dotProduct_lapOfOffDiag_mulVec {o : V → V → ℝ} (ho : ∀ i j, o i j = o j i)
    (u : V → ℝ) :
    u ⬝ᵥ (lapOfOffDiag o *ᵥ u) = 1 / 2 * ∑ i, ∑ j, -o i j * (u i - u j) ^ 2 := by
  have h1 : u ⬝ᵥ (lapOfOffDiag o *ᵥ u) = ∑ i, ∑ j, -o i j * (u i * (u i - u j)) := by
    simp only [dotProduct, lapOfOffDiag_mulVec_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring
  have h2 : ∑ i, ∑ j, -o i j * (u i * (u i - u j)) =
      ∑ i, ∑ j, -o i j * (u j * (u j - u i)) := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by rw [ho j i]
  have h3 : ∑ i, ∑ j, -o i j * (u i - u j) ^ 2 =
      ∑ i, ∑ j, -o i j * (u i * (u i - u j)) + ∑ i, ∑ j, -o i j * (u j * (u j - u i)) := by
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [h1, h3, ← h2]
  ring

/-- Scaling the off-diagonal entries scales the matrix. -/
theorem lapOfOffDiag_smul (c : ℝ) (o : V → V → ℝ) :
    lapOfOffDiag (fun i j => c * o i j) = c • lapOfOffDiag o := by
  ext i j
  simp only [lapOfOffDiag, Matrix.of_apply, Matrix.smul_apply, smul_eq_mul]
  by_cases h : i = j
  · simp only [h, ite_true, mul_neg, Finset.mul_sum, mul_ite, mul_zero]
  · simp only [h, ite_false]

end OffDiag

namespace TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype F] (M : TriMesh V F)

/-- Off-diagonal entry of the weighted cotan Laplacian:
`−½ Σ_f w̄_f · cot(corner of f opposite ij)`. -/
def weightedOffDiag (w : F → ℝ) (i j : V) : ℝ := -(∑ f, w f * M.cotOpp f i j) / 2

theorem weightedOffDiag_symm (w : F → ℝ) (i j : V) :
    M.weightedOffDiag w i j = M.weightedOffDiag w j i := by
  unfold weightedOffDiag
  simp only [M.cotOpp_symm _ i j]

theorem weightedOffDiag_const (c : ℝ) (i j : V) :
    M.weightedOffDiag (fun _ => c) i j = c * M.offDiag i j := by
  unfold weightedOffDiag offDiag
  rw [← Finset.mul_sum]
  ring

/-- Mean of a per-vertex weight over the three vertices of a face (`Ddg.triangleMean`). -/
def faceMean (a : V → ℝ) (f : F) : ℝ := (a (M.a f) + a (M.b f) + a (M.c f)) / 3

variable [Fintype V]

/-- The weighted cotan Laplacian `L_w = Σ_f w̄_f · (cotan Laplacian of the face f)`:
off-diagonal entries `weightedOffDiag`, zero row sums. -/
def weightedLaplacian (w : F → ℝ) : Matrix V V ℝ := lapOfOffDiag (M.weightedOffDiag w)

/-- The cotan Laplacian is the case of the off-diagonal entries `offDiag`. -/
theorem laplacian_eq_lapOfOffDiag : M.laplacian = lapOfOffDiag M.offDiag := rfl

/-- **The weighted Laplacian is symmetric.**

[formula: cut-weighted-laplacian] -/
theorem weightedLaplacian_isSymm (w : F → ℝ) : (M.weightedLaplacian w).IsSymm :=
  lapOfOffDiag_isSymm (M.weightedOffDiag_symm w)

/-- **Zero row sums**: `L_w · 1 = 0`.

[formula: cut-weighted-laplacian] -/
theorem weightedLaplacian_mulVec_one (w : F → ℝ) :
    M.weightedLaplacian w *ᵥ (fun _ => (1 : ℝ)) = 0 :=
  lapOfOffDiag_mulVec_one _

/-- Zero column sums: `1ᵀ L_w = 0`. -/
theorem weightedLaplacian_vecMul_one (w : F → ℝ) :
    (fun _ => (1 : ℝ)) ᵥ* M.weightedLaplacian w = 0 := by
  rw [← Matrix.mulVec_transpose, (M.weightedLaplacian_isSymm w).eq, weightedLaplacian_mulVec_one]

/-- **The quadratic form is the weighted sum of the face energies**:
`xᵀ L_w x = Σ_f w̄_f · ½ (cot_a (x_b − x_c)² + cot_b (x_c − x_a)² + cot_c (x_a − x_b)²)`.
By `cotan_formula` the face energy is the Dirichlet energy `A_f ‖∇x‖²` of the linear
interpolant on that face.

[formula: cut-weighted-laplacian] -/
theorem dotProduct_weightedLaplacian_mulVec (hM : M.Valid) (w : F → ℝ) (u : V → ℝ) :
    u ⬝ᵥ (M.weightedLaplacian w *ᵥ u) = ∑ f, w f * M.faceEnergy f u := by
  rw [weightedLaplacian, dotProduct_lapOfOffDiag_mulVec (M.weightedOffDiag_symm w)]
  have e : ∀ i j, -M.weightedOffDiag w i j * (u i - u j) ^ 2 =
      ∑ f, 1 / 2 * (w f * (M.cotOpp f i j * (u i - u j) ^ 2)) := by
    intro i j
    unfold weightedOffDiag
    rw [neg_div, neg_neg, Finset.sum_div, Finset.sum_mul]
    exact Finset.sum_congr rfl fun f _ => by ring
  have h1 : ∑ i, ∑ j, -M.weightedOffDiag w i j * (u i - u j) ^ 2 =
      ∑ f, ∑ i, ∑ j, 1 / 2 * (w f * (M.cotOpp f i j * (u i - u j) ^ 2)) := by
    simp only [e]
    calc ∑ i, ∑ j, ∑ f, 1 / 2 * (w f * (M.cotOpp f i j * (u i - u j) ^ 2))
        = ∑ i, ∑ f, ∑ j, 1 / 2 * (w f * (M.cotOpp f i j * (u i - u j) ^ 2)) :=
          Finset.sum_congr rfl fun i _ => Finset.sum_comm
      _ = ∑ f, ∑ i, ∑ j, 1 / 2 * (w f * (M.cotOpp f i j * (u i - u j) ^ 2)) := Finset.sum_comm
  rw [h1, Finset.mul_sum]
  refine Finset.sum_congr rfl fun f _ => ?_
  have h2 : ∑ i, ∑ j, 1 / 2 * (w f * (M.cotOpp f i j * (u i - u j) ^ 2)) =
      1 / 2 * (w f * ∑ i, ∑ j, M.cotOpp f i j * (u i - u j) ^ 2) := by
    simp only [Finset.mul_sum]
  rw [h2, M.sum_sum_cotOpp_mul hM f (fun i j => (u i - u j) ^ 2)]
  simp only [faceEnergy, cotEnergy, cotA, cotB, cotC]
  ring

/-- **The weighted Laplacian is positive semidefinite for nonnegative face weights**, on every
valid mesh (no Delaunay or acuteness assumption).

[formula: cut-weighted-laplacian] -/
theorem weightedLaplacian_posSemidef (hM : M.Valid) {w : F → ℝ} (hw : ∀ f, 0 ≤ w f) :
    (M.weightedLaplacian w).PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ?_ fun u => ?_
  · exact M.weightedLaplacian_isSymm w
  · rw [star_trivial, M.dotProduct_weightedLaplacian_mulVec hM]
    exact Finset.sum_nonneg fun f _ => mul_nonneg (hw f) (cotEnergy_nonneg (hM.tri f) _ _ _)

/-- **A constant face weight `c` gives `c·L`.**

[formula: cut-weighted-laplacian] -/
theorem weightedLaplacian_const (c : ℝ) :
    M.weightedLaplacian (fun _ => c) = c • M.laplacian := by
  rw [weightedLaplacian, laplacian_eq_lapOfOffDiag, ← lapOfOffDiag_smul]
  congr 1
  funext i j
  exact M.weightedOffDiag_const c i j

/-- With the weight `1` it is the cotan Laplacian. -/
theorem weightedLaplacian_one : M.weightedLaplacian (fun _ => 1) = M.laplacian := by
  rw [weightedLaplacian_const, one_smul]

/-- **Vertex weights.** With `w̄_f` the mean of a nonnegative per-vertex weight over the three
vertices of `f` (SPEC §3.6, §3.10), the matrix is positive semidefinite.

[formula: cut-weighted-laplacian] -/
theorem weightedLaplacian_faceMean_posSemidef (hM : M.Valid) {a : V → ℝ} (ha : ∀ v, 0 ≤ a v) :
    (M.weightedLaplacian (M.faceMean a)).PosSemidef :=
  M.weightedLaplacian_posSemidef hM fun f => by
    unfold faceMean
    have := ha (M.a f)
    have := ha (M.b f)
    have := ha (M.c f)
    positivity

/-- **A constant vertex weight `c` gives `c·L`.**

[formula: cut-weighted-laplacian] -/
theorem weightedLaplacian_faceMean_const (c : ℝ) :
    M.weightedLaplacian (M.faceMean fun _ => c) = c • M.laplacian := by
  have : M.faceMean (fun _ => c) = fun _ => c := by
    funext f
    unfold faceMean
    ring
  rw [this, weightedLaplacian_const]

/-- The paper's edge rule (Sec. 5.3): the cotangent weight of the edge `ij` is multiplied by
`(a_i + a_j)/2`. -/
def edgeRuleLaplacian (a : V → ℝ) : Matrix V V ℝ :=
  lapOfOffDiag fun i j => (a i + a j) / 2 * M.offDiag i j

/-- For a constant weight the paper's edge rule and the face-mean rule agree: both give `c·L`
(SPEC D11). -/
theorem edgeRule_const (c : ℝ) : M.edgeRuleLaplacian (fun _ => c) = c • M.laplacian := by
  rw [edgeRuleLaplacian, laplacian_eq_lapOfOffDiag, ← lapOfOffDiag_smul]
  congr 1
  funext i j
  ring

end TriMesh

/-! ## The edge rule is not positive semidefinite on an obtuse triangle -/

/-- The quadratic form of the paper's edge rule on one triangle with side lengths `la lb lc`
(opposite the vertices `a b c`), vertex weights `wa wb wc` and vertex values `x y z`:
each term of the cotan energy is multiplied by the mean of the weight over its edge. -/
def edgeRuleEnergy (la lb lc wa wb wc x y z : ℝ) : ℝ :=
  1 / 2 * ((wb + wc) / 2 * cornerCot la lb lc * (y - z) ^ 2 +
    (wc + wa) / 2 * cornerCot lb lc la * (z - x) ^ 2 +
    (wa + wb) / 2 * cornerCot lc la lb * (x - y) ^ 2)

/-- The side lengths `4, 13, 15` form a triangle of area `24` (obtuse at the corner opposite the
side `15`). -/
theorem triArea_4_13_15 : triArea 4 13 15 = 24 := by
  unfold triArea
  rw [show ((4 : ℝ) + 13 + 15) * (13 - 4 + 15) * (4 - 13 + 15) * (4 - 15 + 13) = 96 ^ 2 by
    norm_num, Real.sqrt_sq (by norm_num)]
  norm_num

/-- **The edge rule is indefinite on an obtuse triangle** (SPEC D11): on the triangle with sides
`4, 13, 15`, the nonnegative vertex weight `(1, 1, 0)` (zero at the obtuse corner) and the values
`(1, −1, −1)` give a negative edge-rule energy, while the face-mean rule gives
`w̄ · cotEnergy ≥ 0` for every nonnegative weight (`weightedLaplacian_posSemidef`). -/
theorem edgeRuleEnergy_neg_example : edgeRuleEnergy 4 13 15 1 1 0 1 (-1) (-1) < 0 := by
  have h1 : triArea 13 15 4 = 24 := by rw [triArea_rotate, triArea_4_13_15]
  have h2 : triArea 15 4 13 = 24 := by rw [triArea_rotate, h1]
  unfold edgeRuleEnergy cornerCot
  rw [triArea_4_13_15, h1, h2]
  norm_num

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
