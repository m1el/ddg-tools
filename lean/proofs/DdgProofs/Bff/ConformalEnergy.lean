import DDG.Mesh.Orientation
import DdgProofs.TriangleIntegral.Fem
import DdgProofs.Bff.Curve
import DDG.Laplace.Cotan

/-!
# The conformal energy of a flattening over an oriented mesh (BFF Eqs. 15–16)

`bff-conformal-energy`: for an oriented valid mesh whose boundary halfedges form the loop
`(next, emb)` (`TriMesh.IsBoundaryLoop`, the surface on the left of every boundary edge, as
`hilbert_transform` of `crates/ddg-bff/src/disk.rs` walks it), and any map `f = a + ib` given by
its vertex values:

* **shoelace over the mesh** (`TriMesh.sum_sda_eq_two_mul_loopSignedArea`): the signed image
  areas of the faces add up to the signed area of the image of the boundary polygon,
  `Σ_f ½ sda(f(a_f), f(b_f), f(c_f)) = 𝒜(a, b)`; the contributions of the interior edges cancel
  because the two faces at such an edge traverse it in opposite directions
  (`TriMesh.sum_halfedge_antisymm`);
* **the energy** (`TriMesh.two_mul_sum_integral_dbarSq`): with every face laid out isometrically
  and counter-clockwise and `f` linear on each face,
  `2 Σ_f ∫_{T_f} |∂̄f|² = ½(aᵀLa + bᵀLb) − 𝒜(a, b) = loopConformalEnergy`, the energy whose
  minimizers over `b` are the solutions of the Hilbert-transform system
  (`hilbertTransform_minimizes_conformalEnergy_iff`, `[formula: bff-hilbert-transform]`);
* `E_C(a, b) ≥ 0`, and `E_C(a, b) = 0` iff `f` satisfies the Cauchy–Riemann equations on every
  face (`TriMesh.loopConformalEnergy_eq_zero_iff`).

The per-face identity is `[formula: bff-conformal-energy-triangle]` (`TriangleIntegral/Fem.lean`).
-/

noncomputable section

namespace DDG

namespace TriMesh

open TriangleIntegral MeasureTheory

variable {V F β : Type*} (M : TriMesh V F)

/-- The boundary of `M` is the loop `(next, emb)`, traversed in the direction of the faces: `emb`
is injective, every loop edge `emb j → emb (next j)` is a halfedge without a twin, and every
halfedge without a twin is a loop edge. (`next` may have several cycles: several boundary
loops.) -/
structure IsBoundaryLoop (next : Equiv.Perm β) (emb : β → V) : Prop where
  injective : Function.Injective emb
  edge : ∀ j, ∃ h, ¬M.HasTwin h ∧ M.halfedge h = (emb j, emb (next j))
  cover : ∀ h, ¬M.HasTwin h → ∃ j, M.halfedge h = (emb j, emb (next j))

/-- The shoelace form of the signed double area of a triangle. -/
theorem sda_eq_shoelace (p q r : ℝ × ℝ) :
    sda p q r = (p.1 * q.2 - q.1 * p.2) + (q.1 * r.2 - r.1 * q.2) + (r.1 * p.2 - p.1 * r.2) := by
  unfold sda
  ring

/-- **The signed image areas of the faces add up to the signed area of the boundary polygon**
(the shoelace formula over an oriented mesh): for every `a, b`,
`Σ_f sda((a, b)(a_f), (a, b)(b_f), (a, b)(c_f)) = 2·𝒜(a, b)`, `𝒜 = loopSignedArea`.

[formula: bff-conformal-energy] -/
theorem sum_sda_eq_two_mul_loopSignedArea [Fintype F] [Fintype β] (hO : M.Oriented)
    {next : Equiv.Perm β} {emb : β → V} (hL : M.IsBoundaryLoop next emb) (a b : V → ℝ) :
    ∑ f, sda (a (M.a f), b (M.a f)) (a (M.b f), b (M.b f)) (a (M.c f), b (M.c f)) =
      2 * loopSignedArea next (a ∘ emb) (b ∘ emb) := by
  classical
  set w : V → V → ℝ := fun x y => a x * b y - a y * b x with hw
  have hanti : ∀ x y, w y x = -w x y := fun x y => by simp only [hw]; ring
  have hface : ∀ f, sda (a (M.a f), b (M.a f)) (a (M.b f), b (M.b f)) (a (M.c f), b (M.c f)) =
      ∑ k : Fin 3, w (M.halfedge (f, k)).1 (M.halfedge (f, k)).2 := by
    intro f
    rw [Fin.sum_univ_three, halfedge_zero, halfedge_one, halfedge_two, sda_eq_shoelace]
  rw [Finset.sum_congr rfl fun f _ => hface f, ← Fintype.sum_prod_type',
    M.sum_halfedge_antisymm hO w hanti]
  have hbij : ∑ j, w (emb j) (emb (next j)) =
      ∑ h ∈ Finset.univ.filter (fun h => ¬M.HasTwin h), w (M.halfedge h).1 (M.halfedge h).2 := by
    refine Finset.sum_nbij (fun j => (hL.edge j).choose) ?_ ?_ ?_ ?_
    · intro j _
      simpa using (hL.edge j).choose_spec.1
    · intro j _ j' _ e
      have h1 := (hL.edge j).choose_spec.2
      have h2 := (hL.edge j').choose_spec.2
      simp only at e
      rw [e, h2] at h1
      exact hL.injective (congrArg Prod.fst h1).symm
    · intro h hh
      simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hh
      obtain ⟨j, hj⟩ := hL.cover h hh
      refine ⟨j, by simp, hO ?_⟩
      rw [(hL.edge j).choose_spec.2, hj]
    · intro j _
      rw [(hL.edge j).choose_spec.2]
  rw [← hbij, loopSignedArea, Finset.mul_sum, ← Finset.mul_sum]
  simp only [hw, Function.comp]
  ring

/-- **Every valid mesh has a counter-clockwise isometric layout** (the reference layout of
`exists_isLayout`, which has `D = 2·area > 0`). -/
theorem exists_isLayout_pos (hM : M.Valid) : ∃ P, M.IsLayout P ∧ ∀ f, 0 < D (P f) := by
  refine ⟨fun f => ![(0, 0), (M.lc f, 0), (layoutX (M.lc f) (M.la f) (M.lb f),
    layoutY (M.lc f) (triArea (M.lc f) (M.la f) (M.lb f)))], ⟨fun f => ?_, fun f => ?_,
    fun f => ?_, fun f => ?_⟩, fun f => ?_⟩ <;>
  have ht : IsTriangle (M.lc f) (M.la f) (M.lb f) := (hM.tri f).rotate.rotate <;>
  obtain ⟨_, h01, h02, h12⟩ := layout_isometric ht rfl
  · show sda _ _ _ ≠ 0
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons]
    rw [sda_layout ht.pos_a.ne']
    exact (mul_pos two_pos (triArea_pos ht)).ne'
  · rw [sideLen_zero]
    simpa only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons] using h12
  · rw [sideLen_one, Intrinsic.dist2_comm]
    simpa only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons] using h02
  · rw [sideLen_two]
    simpa only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons] using h01
  · show 0 < sda _ _ _
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons]
    rw [sda_layout ht.pos_a.ne']
    exact mul_pos two_pos (triArea_pos ht)

variable [DecidableEq V] [Fintype V] [Fintype F] [Fintype β]

/-- **The conformal energy of BFF is twice the integral of `|∂̄f|²`** (BFF Eq. 15; Hutchinson;
Lévy et al.): on a valid oriented mesh with boundary loop `(next, emb)`, every face laid out
isometrically and counter-clockwise, and `f = a + ib` linear on every face,
`2 Σ_f ∫_{T_f} |∂̄f|² = ½(aᵀLa + bᵀLb) − 𝒜(a, b)`, the energy `E_C` whose minimizers over `b`
`hilbert_transform` computes (`hilbertTransform_cotan`).

[formula: bff-conformal-energy] -/
theorem two_mul_sum_integral_dbarSq (hM : M.Valid) (hO : M.Oriented) {next : Equiv.Perm β}
    {emb : β → V} (hL : M.IsBoundaryLoop next emb) {P : F → Fin 3 → ℝ × ℝ} (hP : M.IsLayout P)
    (hpos : ∀ f, 0 < D (P f)) (a b : V → ℝ) :
    2 * ∑ f, ∫ x in triangle (P f),
        dbarSq (interp (P f) (M.cornerVals a f)) (interp (P f) (M.cornerVals b f)) x =
      loopConformalEnergy next emb M.laplacian a b := by
  rw [Finset.mul_sum]
  have hf : ∀ f, 2 * ∫ x in triangle (P f),
        dbarSq (interp (P f) (M.cornerVals a f)) (interp (P f) (M.cornerVals b f)) x =
      (M.faceEnergy f a + M.faceEnergy f b) / 2 -
        sda (a (M.a f), b (M.a f)) (a (M.b f), b (M.b f)) (a (M.c f), b (M.c f)) / 2 := by
    intro f
    rw [integral_dbarSq (hpos f), hP.la, hP.lb, hP.lc]
    rfl
  rw [Finset.sum_congr rfl fun f _ => hf f, Finset.sum_sub_distrib, ← Finset.sum_div,
    ← Finset.sum_div, Finset.sum_add_distrib, ← M.dotProduct_laplacian_mulVec_eq_sum_faceEnergy hM,
    ← M.dotProduct_laplacian_mulVec_eq_sum_faceEnergy hM,
    M.sum_sda_eq_two_mul_loopSignedArea hO hL a b, loopConformalEnergy]
  ring

/-- **`E_C ≥ 0`, with equality exactly for maps that satisfy the Cauchy–Riemann equations on
every face** (discrete conformal maps, BFF Sec. 4.3.3).

[formula: bff-conformal-energy] -/
theorem loopConformalEnergy_eq_zero_iff (hM : M.Valid) (hO : M.Oriented) {next : Equiv.Perm β}
    {emb : β → V} (hL : M.IsBoundaryLoop next emb) {P : F → Fin 3 → ℝ × ℝ} (hP : M.IsLayout P)
    (hpos : ∀ f, 0 < D (P f)) (a b : V → ℝ) :
    0 ≤ loopConformalEnergy next emb M.laplacian a b ∧
      (loopConformalEnergy next emb M.laplacian a b = 0 ↔
        ∀ f, (gradInterp (P f) (M.cornerVals a f)).1 = (gradInterp (P f) (M.cornerVals b f)).2 ∧
          (gradInterp (P f) (M.cornerVals a f)).2 = -(gradInterp (P f) (M.cornerVals b f)).1) := by
  rw [← M.two_mul_sum_integral_dbarSq hM hO hL hP hpos a b]
  have hnn : ∀ f, 0 ≤ ∫ x in triangle (P f),
      dbarSq (interp (P f) (M.cornerVals a f)) (interp (P f) (M.cornerVals b f)) x := fun f =>
    setIntegral_nonneg (measurableSet_triangle _) fun x _ => by unfold dbarSq; positivity
  refine ⟨mul_nonneg two_pos.le (Finset.sum_nonneg fun f _ => hnn f), ?_⟩
  rw [mul_eq_zero, or_iff_right two_ne_zero,
    Finset.sum_eq_zero_iff_of_nonneg fun f _ => hnn f]
  simp only [Finset.mem_univ, true_implies]
  exact forall_congr' fun f => integral_dbarSq_eq_zero_iff (hpos f).ne' _ _

end TriMesh

end DDG
