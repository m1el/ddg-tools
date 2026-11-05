import DdgProofs.Triangle

/-!
# The cotan formula, the cotan Laplacian and the lumped mass

* `cotan_formula`: on a Euclidean triangle the Dirichlet energy `Area·|∇f|²` of an affine function
  is `½ Σ cot(θ_k)(f_i − f_j)²`, the sum over the three edges `ij`, `θ_k` the angle opposite;
* `cotEnergy_nonneg`, `cotEnergy_eq_zero_iff`: this quadratic form is positive semidefinite with
  kernel the constants — for every triangle, also an obtuse one (negative cotangents);
* for the assembled matrix `L` of a triangle mesh (`TriMesh.laplacian`, entries as in `CLAUDE.md`):
  symmetric, zero row sums, `uᵀLu = Σ_faces (face energy)`, positive semidefinite (no Delaunay
  assumption), kernel equal to the constants on a connected mesh, interior block positive
  definite when the boundary is nonempty;
* `TriMesh.sum_mass`: the lumped masses sum to the total area;
* invariance under uniform scaling.
-/

noncomputable section

namespace DDG

open Matrix

/-! ## One triangle -/

/-- The cotan energy of three vertex values `x y z` on a triangle with side lengths `a b c`,
`a` opposite the vertex carrying `x`, `b` opposite `y`, `c` opposite `z`:
`½ (cot A·(y − z)² + cot B·(z − x)² + cot C·(x − y)²)`. -/
def cotEnergy (a b c x y z : ℝ) : ℝ :=
  1 / 2 * (cornerCot a b c * (y - z) ^ 2 + cornerCot b c a * (z - x) ^ 2 +
    cornerCot c a b * (x - y) ^ 2)

section OneTriangle

variable {a b c : ℝ}

private theorem cot_mul_area (h : IsTriangle a b c) :
    cornerCot a b c * (4 * triArea a b c) = b ^ 2 + c ^ 2 - a ^ 2 ∧
    cornerCot b c a * (4 * triArea a b c) = c ^ 2 + a ^ 2 - b ^ 2 ∧
    cornerCot c a b * (4 * triArea a b c) = a ^ 2 + b ^ 2 - c ^ 2 := by
  have hT := (triArea_pos h).ne'
  refine ⟨?_, ?_, ?_⟩
  · unfold cornerCot; field_simp
  · unfold cornerCot; rw [triArea_rotate]; field_simp
  · unfold cornerCot; rw [triArea_rotate, triArea_rotate]; field_simp

/-- Completing the square in the cotan energy. -/
theorem cotEnergy_identity (h : IsTriangle a b c) (x y z : ℝ) :
    16 * b ^ 2 * triArea a b c * cotEnergy a b c x y z =
      (2 * b ^ 2 * (y - z) + (a ^ 2 + b ^ 2 - c ^ 2) * (z - x)) ^ 2 +
        16 * triArea a b c ^ 2 * (z - x) ^ 2 := by
  have hT := sixteen_mul_triArea_sq h
  obtain ⟨cA, cB, cC⟩ := cot_mul_area h
  unfold cotEnergy
  linear_combination (2 * b ^ 2 * (y - z) ^ 2) * cA + (2 * b ^ 2 * (z - x) ^ 2) * cB +
    (2 * b ^ 2 * (x - y) ^ 2) * cC - (z - x) ^ 2 * hT

/-- **The cotan energy of a triangle is nonnegative**, whatever the shape of the triangle
(individual cotangents may be negative).

[formula: cotan-laplacian-psd] -/
theorem cotEnergy_nonneg (h : IsTriangle a b c) (x y z : ℝ) : 0 ≤ cotEnergy a b c x y z := by
  have hid := cotEnergy_identity h x y z
  have hpos : 0 < 16 * b ^ 2 * triArea a b c := by
    have := h.pos_b
    have := triArea_pos h
    positivity
  by_contra hneg
  push Not at hneg
  have h1 : 16 * b ^ 2 * triArea a b c * cotEnergy a b c x y z < 0 :=
    mul_neg_of_pos_of_neg hpos hneg
  nlinarith [sq_nonneg (2 * b ^ 2 * (y - z) + (a ^ 2 + b ^ 2 - c ^ 2) * (z - x)),
    sq_nonneg (triArea a b c * (z - x))]

/-- The cotan energy of a triangle vanishes only on constants. -/
theorem cotEnergy_eq_zero_iff (h : IsTriangle a b c) (x y z : ℝ) :
    cotEnergy a b c x y z = 0 ↔ x = y ∧ y = z := by
  constructor
  · intro h0
    have hid := cotEnergy_identity h x y z
    rw [h0, mul_zero] at hid
    have hb := h.pos_b
    have hT := triArea_pos h
    have h1 : 0 ≤ (2 * b ^ 2 * (y - z) + (a ^ 2 + b ^ 2 - c ^ 2) * (z - x)) ^ 2 := sq_nonneg _
    have h2 : 0 ≤ 16 * triArea a b c ^ 2 * (z - x) ^ 2 := by positivity
    have h3 : 16 * triArea a b c ^ 2 * (z - x) ^ 2 = 0 := by linarith
    have h4 : (2 * b ^ 2 * (y - z) + (a ^ 2 + b ^ 2 - c ^ 2) * (z - x)) ^ 2 = 0 := by linarith
    have hzx : z - x = 0 := by
      have : (z - x) ^ 2 = 0 := by
        rcases mul_eq_zero.mp h3 with h5 | h5
        · exfalso; have : 0 < 16 * triArea a b c ^ 2 := by positivity
          linarith
        · exact h5
      exact pow_eq_zero_iff (two_ne_zero) |>.mp this
    have h5 := pow_eq_zero_iff (two_ne_zero) |>.mp h4
    rw [hzx, mul_zero, add_zero] at h5
    have hyz : y - z = 0 := by
      rcases mul_eq_zero.mp h5 with h6 | h6
      · exfalso; have : 0 < 2 * b ^ 2 := by positivity
        linarith
      · exact h6
    constructor <;> linarith
  · rintro ⟨rfl, rfl⟩
    simp [cotEnergy]

end OneTriangle

/-- **The cotan formula.** Take a nondegenerate triangle with vertices `p₁`, `p₂ = p₁ + e₁`,
`p₃ = p₁ + e₂` in a Euclidean space and an affine function `f` on its plane with gradient
`g = α e₁ + β e₂` and `f(p₁) = 0`, so `f(p₂) = ⟪g, e₁⟫`, `f(p₃) = ⟪g, e₂⟫`. Then

`Area · ‖g‖² = ½ (cot θ₁ (f₂ − f₃)² + cot θ₂ (f₃ − f₁)² + cot θ₃ (f₁ − f₂)²)`,

i.e. the Dirichlet energy `∫_T |∇f|² dA` (the integrand is constant) equals the cotan energy, all
quantities on the right computed from the three edge lengths.

[formula: cotan-formula] -/
theorem cotan_formula {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] (e₁ e₂ : E)
    (α β : ℝ) (h : IsTriangle ‖e₂ - e₁‖ ‖e₂‖ ‖e₁‖) :
    triArea ‖e₂ - e₁‖ ‖e₂‖ ‖e₁‖ * ‖α • e₁ + β • e₂‖ ^ 2 =
      cotEnergy ‖e₂ - e₁‖ ‖e₂‖ ‖e₁‖ 0 (inner ℝ (α • e₁ + β • e₂) e₁)
        (inner ℝ (α • e₁ + β • e₂) e₂) := by
  have hT := sixteen_mul_triArea_sq h
  have hTpos := triArea_pos h
  obtain ⟨cA, cB, cC⟩ := cot_mul_area h
  have hc2 : inner ℝ e₁ e₁ = ‖e₁‖ ^ 2 := real_inner_self_eq_norm_sq e₁
  have hb2 : inner ℝ e₂ e₂ = ‖e₂‖ ^ 2 := real_inner_self_eq_norm_sq e₂
  have hsymm : inner ℝ e₂ e₁ = inner ℝ e₁ e₂ := real_inner_comm _ _
  have ha2 : ‖e₂ - e₁‖ ^ 2 = ‖e₂‖ ^ 2 - 2 * inner ℝ e₁ e₂ + ‖e₁‖ ^ 2 := by
    rw [norm_sub_sq_real, hsymm]
  have hg : ‖α • e₁ + β • e₂‖ ^ 2 =
      α ^ 2 * ‖e₁‖ ^ 2 + 2 * α * β * inner ℝ e₁ e₂ + β ^ 2 * ‖e₂‖ ^ 2 := by
    rw [← real_inner_self_eq_norm_sq]
    simp only [inner_add_left, inner_add_right, inner_smul_left, inner_smul_right, hc2, hb2,
      hsymm, RCLike.conj_to_real]
    ring
  have hf2 : inner ℝ (α • e₁ + β • e₂) e₁ = α * ‖e₁‖ ^ 2 + β * inner ℝ e₁ e₂ := by
    simp only [inner_add_left, inner_smul_left, hc2, hsymm, RCLike.conj_to_real]
  have hf3 : inner ℝ (α • e₁ + β • e₂) e₂ = α * inner ℝ e₁ e₂ + β * ‖e₂‖ ^ 2 := by
    simp only [inner_add_left, inner_smul_left, hb2, RCLike.conj_to_real]
  rw [hg, hf2, hf3]
  generalize inner ℝ e₁ e₂ = F at *
  generalize ‖e₂ - e₁‖ = a at *
  generalize ‖e₂‖ = b at *
  generalize ‖e₁‖ = c at *
  generalize triArea a b c = T at *
  have h4 : a ^ 4 = (b ^ 2 - 2 * F + c ^ 2) ^ 2 := by rw [← ha2]; ring
  have hT' : 4 * T ^ 2 = b ^ 2 * c ^ 2 - F ^ 2 := by
    linear_combination (1 / 4) * hT + (1 / 2) * (b ^ 2 + c ^ 2) * ha2 - (1 / 4) * h4
  have cA' : cornerCot a b c * (4 * T) = 2 * F := by linear_combination cA - ha2
  have cB' : cornerCot b c a * (4 * T) = 2 * c ^ 2 - 2 * F := by linear_combination cB + ha2
  have cC' : cornerCot c a b * (4 * T) = 2 * b ^ 2 - 2 * F := by linear_combination cC + ha2
  have key : 8 * T * cotEnergy a b c 0 (α * c ^ 2 + β * F) (α * F + β * b ^ 2) =
      8 * T * (T * (α ^ 2 * c ^ 2 + 2 * α * β * F + β ^ 2 * b ^ 2)) := by
    unfold cotEnergy
    linear_combination ((α * c ^ 2 + β * F) - (α * F + β * b ^ 2)) ^ 2 * cA' +
      (α * F + β * b ^ 2) ^ 2 * cB' + (α * c ^ 2 + β * F) ^ 2 * cC' -
      2 * (α ^ 2 * c ^ 2 + 2 * α * β * F + β ^ 2 * b ^ 2) * hT'
  exact (mul_left_cancel₀ (by positivity) key).symm

/-! ## The assembled cotan Laplacian -/

/-- A triangle mesh given intrinsically: every face `f` has three vertices `a f`, `b f`, `c f`,
and `len i j` is the length of the edge between vertices `i` and `j` (only evaluated on pairs of
vertices of a common face). -/
structure TriMesh (V F : Type*) where
  a : F → V
  b : F → V
  c : F → V
  len : V → V → ℝ

namespace TriMesh

variable {V F : Type*} (M : TriMesh V F)

/-- Length of the side of `f` opposite its corner `a`. -/
def la (f : F) : ℝ := M.len (M.b f) (M.c f)
/-- Length of the side of `f` opposite its corner `b`. -/
def lb (f : F) : ℝ := M.len (M.c f) (M.a f)
/-- Length of the side of `f` opposite its corner `c`. -/
def lc (f : F) : ℝ := M.len (M.a f) (M.b f)

/-- Cotangent of the corner of `f` at `a`. -/
def cotA (f : F) : ℝ := cornerCot (M.la f) (M.lb f) (M.lc f)
/-- Cotangent of the corner of `f` at `b`. -/
def cotB (f : F) : ℝ := cornerCot (M.lb f) (M.lc f) (M.la f)
/-- Cotangent of the corner of `f` at `c`. -/
def cotC (f : F) : ℝ := cornerCot (M.lc f) (M.la f) (M.lb f)

/-- Interior angle of `f` at `a`. -/
def angleA (f : F) : ℝ := cornerAngle (M.la f) (M.lb f) (M.lc f)
/-- Interior angle of `f` at `b`. -/
def angleB (f : F) : ℝ := cornerAngle (M.lb f) (M.lc f) (M.la f)
/-- Interior angle of `f` at `c`. -/
def angleC (f : F) : ℝ := cornerAngle (M.lc f) (M.la f) (M.lb f)

/-- Area of `f`. -/
def area (f : F) : ℝ := triArea (M.la f) (M.lb f) (M.lc f)

/-- Hypotheses on the mesh used below: every face has three distinct vertices and its three edge
lengths satisfy the strict triangle inequalities. (No manifoldness, no consistency of the lengths
across faces and no Delaunay condition are needed for the results of this file.) -/
structure Valid : Prop where
  ab : ∀ f, M.a f ≠ M.b f
  bc : ∀ f, M.b f ≠ M.c f
  ca : ∀ f, M.c f ≠ M.a f
  tri : ∀ f, IsTriangle (M.la f) (M.lb f) (M.lc f)

/-- The cotan energy of `u` on the face `f`:
`½ (cot_a (u_b − u_c)² + cot_b (u_c − u_a)² + cot_c (u_a − u_b)²)`. -/
def faceEnergy (f : F) (u : V → ℝ) : ℝ :=
  cotEnergy (M.la f) (M.lb f) (M.lc f) (u (M.a f)) (u (M.b f)) (u (M.c f))

variable [DecidableEq V]

/-- Cotangent of the corner of `f` opposite the edge `{i, j}`; `0` if `{i, j}` is not an edge
of `f`. -/
def cotOpp (f : F) (i j : V) : ℝ :=
  if (i = M.a f ∧ j = M.b f) ∨ (i = M.b f ∧ j = M.a f) then M.cotC f
  else if (i = M.b f ∧ j = M.c f) ∨ (i = M.c f ∧ j = M.b f) then M.cotA f
  else if (i = M.c f ∧ j = M.a f) ∨ (i = M.a f ∧ j = M.c f) then M.cotB f
  else 0

theorem cotOpp_symm (f : F) (i j : V) : M.cotOpp f i j = M.cotOpp f j i := by
  unfold cotOpp
  have h1 : ((i = M.a f ∧ j = M.b f) ∨ (i = M.b f ∧ j = M.a f)) ↔
      ((j = M.a f ∧ i = M.b f) ∨ (j = M.b f ∧ i = M.a f)) := by tauto
  have h2 : ((i = M.b f ∧ j = M.c f) ∨ (i = M.c f ∧ j = M.b f)) ↔
      ((j = M.b f ∧ i = M.c f) ∨ (j = M.c f ∧ i = M.b f)) := by tauto
  have h3 : ((i = M.c f ∧ j = M.a f) ∨ (i = M.a f ∧ j = M.c f)) ↔
      ((j = M.c f ∧ i = M.a f) ∨ (j = M.a f ∧ i = M.c f)) := by tauto
  simp only [h1, h2, h3]

variable [Fintype F]

/-- Off-diagonal entry of the cotan Laplacian (`CLAUDE.md`):
`L_ij = −½ Σ cot(angle opposite edge ij)` over the faces containing the edge. -/
def offDiag (i j : V) : ℝ := -(∑ f, M.cotOpp f i j) / 2

theorem offDiag_symm (i j : V) : M.offDiag i j = M.offDiag j i := by
  unfold offDiag
  simp only [M.cotOpp_symm _ i j]

variable [Fintype V]

/-- The cotan Laplacian (`CLAUDE.md`): `L_ij = −½ Σ cot(angle opposite edge ij)` for `i ≠ j`,
`L_ii = −Σ_{j ≠ i} L_ij`. -/
def laplacian : Matrix V V ℝ :=
  Matrix.of fun i j => if i = j then -∑ k, (if k = i then 0 else M.offDiag i k) else M.offDiag i j

/-- **The cotan Laplacian is symmetric.**

[formula: cotan-laplacian] -/
theorem laplacian_isSymm : M.laplacian.IsSymm := by
  ext i j
  simp only [laplacian, Matrix.transpose_apply, Matrix.of_apply]
  by_cases h : i = j
  · subst h; simp
  · rw [ite_eq_right h, ite_eq_right (Ne.symm h)]
    exact M.offDiag_symm j i

/-- Row `i` of `L` applied to `u`: `(L u)_i = Σ_j w_ij (u_i − u_j)` with the edge weights
`w_ij = −L_ij = ½ Σ cot`. -/
theorem laplacian_mulVec_apply (u : V → ℝ) (i : V) :
    (M.laplacian *ᵥ u) i = ∑ j, -M.offDiag i j * (u i - u j) := by
  set o : V → ℝ := fun j => if j = i then 0 else M.offDiag i j with ho
  have h1 : (M.laplacian *ᵥ u) i = ∑ j, ((if i = j then -∑ k, o k else 0) + o j) * u j := by
    simp only [Matrix.mulVec, dotProduct, laplacian, Matrix.of_apply]
    refine Finset.sum_congr rfl fun j _ => ?_
    by_cases h : i = j
    · subst h; simp [ho]
    · simp [ho, h, Ne.symm h]
  have h2 : ∑ j, -M.offDiag i j * (u i - u j) = ∑ j, -o j * (u i - u j) := by
    refine Finset.sum_congr rfl fun j _ => ?_
    by_cases h : j = i
    · subst h; simp
    · simp [ho, h]
  rw [h1, h2]
  simp only [add_mul, Finset.sum_add_distrib, ite_mul, zero_mul, Finset.sum_ite_eq,
    Finset.mem_univ, ite_true, neg_mul, mul_sub, Finset.sum_neg_distrib, Finset.sum_sub_distrib,
    ← Finset.sum_mul]
  ring

/-- **Zero row sums**: `L·1 = 0`.

[formula: cotan-laplacian] -/
theorem laplacian_mulVec_one : M.laplacian *ᵥ (fun _ => (1 : ℝ)) = 0 := by
  ext i
  rw [laplacian_mulVec_apply]
  simp

/-- **Zero column sums**: `1ᵀL = 0`. -/
theorem laplacian_vecMul_one : (fun _ => (1 : ℝ)) ᵥ* M.laplacian = 0 := by
  rw [← Matrix.mulVec_transpose, M.laplacian_isSymm.eq, laplacian_mulVec_one]

/-- The quadratic form as a sum over ordered pairs of vertices. -/
theorem dotProduct_laplacian_mulVec (u : V → ℝ) :
    u ⬝ᵥ (M.laplacian *ᵥ u) = 1 / 2 * ∑ i, ∑ j, -M.offDiag i j * (u i - u j) ^ 2 := by
  have h1 : u ⬝ᵥ (M.laplacian *ᵥ u) = ∑ i, ∑ j, -M.offDiag i j * (u i * (u i - u j)) := by
    simp only [dotProduct, laplacian_mulVec_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring
  have h2 : ∑ i, ∑ j, -M.offDiag i j * (u i * (u i - u j)) =
      ∑ i, ∑ j, -M.offDiag i j * (u j * (u j - u i)) := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by
      rw [M.offDiag_symm j i]
  have h3 : ∑ i, ∑ j, -M.offDiag i j * (u i - u j) ^ 2 =
      ∑ i, ∑ j, -M.offDiag i j * (u i * (u i - u j)) +
        ∑ i, ∑ j, -M.offDiag i j * (u j * (u j - u i)) := by
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [h1, h3, ← h2]
  ring

omit [Fintype F] [Fintype V] in
/-- With three distinct vertices, `cotOpp` is a sum of six indicator terms. -/
theorem cotOpp_eq (hM : M.Valid) (f : F) (i j : V) :
    M.cotOpp f i j =
      M.cotC f * ((if i = M.a f then 1 else 0) * (if j = M.b f then 1 else 0) +
        (if i = M.b f then 1 else 0) * (if j = M.a f then 1 else 0)) +
      M.cotA f * ((if i = M.b f then 1 else 0) * (if j = M.c f then 1 else 0) +
        (if i = M.c f then 1 else 0) * (if j = M.b f then 1 else 0)) +
      M.cotB f * ((if i = M.c f then 1 else 0) * (if j = M.a f then 1 else 0) +
        (if i = M.a f then 1 else 0) * (if j = M.c f then 1 else 0)) := by
  have hab := hM.ab f
  have hbc := hM.bc f
  have hca := hM.ca f
  unfold cotOpp
  by_cases hia : i = M.a f <;> by_cases hib : i = M.b f <;> by_cases hic : i = M.c f <;>
    by_cases hja : j = M.a f <;> by_cases hjb : j = M.b f <;> by_cases hjc : j = M.c f <;>
    simp_all

private theorem sum_sum_ind (x y : V) (g : V → V → ℝ) :
    ∑ i, ∑ j, ((if i = x then (1 : ℝ) else 0) * (if j = y then 1 else 0)) * g i j = g x y := by
  simp

omit [Fintype F] in
/-- Summing `cotOpp f i j · g i j` over all ordered pairs `(i, j)` picks out the six ordered
edges of the face. -/
theorem sum_sum_cotOpp_mul (hM : M.Valid) (f : F) (g : V → V → ℝ) :
    ∑ i, ∑ j, M.cotOpp f i j * g i j =
      M.cotC f * (g (M.a f) (M.b f) + g (M.b f) (M.a f)) +
      M.cotA f * (g (M.b f) (M.c f) + g (M.c f) (M.b f)) +
      M.cotB f * (g (M.c f) (M.a f) + g (M.a f) (M.c f)) := by
  have key : ∀ i j, M.cotOpp f i j * g i j =
      M.cotC f * (((if i = M.a f then (1 : ℝ) else 0) * (if j = M.b f then 1 else 0)) * g i j) +
      M.cotC f * (((if i = M.b f then (1 : ℝ) else 0) * (if j = M.a f then 1 else 0)) * g i j) +
      M.cotA f * (((if i = M.b f then (1 : ℝ) else 0) * (if j = M.c f then 1 else 0)) * g i j) +
      M.cotA f * (((if i = M.c f then (1 : ℝ) else 0) * (if j = M.b f then 1 else 0)) * g i j) +
      M.cotB f * (((if i = M.c f then (1 : ℝ) else 0) * (if j = M.a f then 1 else 0)) * g i j) +
      M.cotB f * (((if i = M.a f then (1 : ℝ) else 0) * (if j = M.c f then 1 else 0)) * g i j) := by
    intro i j
    rw [M.cotOpp_eq hM]
    ring
  simp only [key, Finset.sum_add_distrib, ← Finset.mul_sum, sum_sum_ind]
  ring

/-- **The quadratic form of the cotan Laplacian is the sum of the face energies**:
`uᵀ L u = Σ_f ½ (cot_a (u_b − u_c)² + cot_b (u_c − u_a)² + cot_c (u_a − u_b)²)`.
By `cotan_formula` each summand is the Dirichlet energy of the piecewise linear interpolant on
that face.

[formula: cotan-laplacian] -/
theorem dotProduct_laplacian_mulVec_eq_sum_faceEnergy (hM : M.Valid) (u : V → ℝ) :
    u ⬝ᵥ (M.laplacian *ᵥ u) = ∑ f, M.faceEnergy f u := by
  rw [dotProduct_laplacian_mulVec]
  have e : ∀ i j, -M.offDiag i j * (u i - u j) ^ 2 =
      ∑ f, 1 / 2 * (M.cotOpp f i j * (u i - u j) ^ 2) := by
    intro i j
    unfold offDiag
    rw [neg_div, neg_neg, Finset.sum_div, Finset.sum_mul]
    exact Finset.sum_congr rfl fun f _ => by ring
  have h1 : ∑ i, ∑ j, -M.offDiag i j * (u i - u j) ^ 2 =
      ∑ f, ∑ i, ∑ j, 1 / 2 * (M.cotOpp f i j * (u i - u j) ^ 2) := by
    simp only [e]
    calc ∑ i, ∑ j, ∑ f, 1 / 2 * (M.cotOpp f i j * (u i - u j) ^ 2)
        = ∑ i, ∑ f, ∑ j, 1 / 2 * (M.cotOpp f i j * (u i - u j) ^ 2) :=
          Finset.sum_congr rfl fun i _ => Finset.sum_comm
      _ = ∑ f, ∑ i, ∑ j, 1 / 2 * (M.cotOpp f i j * (u i - u j) ^ 2) := Finset.sum_comm
  rw [h1, Finset.mul_sum]
  refine Finset.sum_congr rfl fun f _ => ?_
  have h2 : ∑ i, ∑ j, 1 / 2 * (M.cotOpp f i j * (u i - u j) ^ 2) =
      1 / 2 * ∑ i, ∑ j, M.cotOpp f i j * (u i - u j) ^ 2 := by
    simp only [Finset.mul_sum]
  rw [h2, M.sum_sum_cotOpp_mul hM f (fun i j => (u i - u j) ^ 2)]
  simp only [faceEnergy, cotEnergy, cotA, cotB, cotC]
  ring

/-- **The cotan Laplacian is positive semidefinite**, for every valid mesh — no Delaunay or
acuteness assumption.

[formula: cotan-laplacian-psd] -/
theorem laplacian_posSemidef (hM : M.Valid) : M.laplacian.PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ?_ fun u => ?_
  · exact M.laplacian_isSymm
  · rw [star_trivial, M.dotProduct_laplacian_mulVec_eq_sum_faceEnergy hM]
    exact Finset.sum_nonneg fun f _ => cotEnergy_nonneg (hM.tri f) _ _ _

/-- Two vertices are adjacent if they belong to a common face. -/
def Adj (i j : V) : Prop :=
  ∃ f, (i = M.a f ∨ i = M.b f ∨ i = M.c f) ∧ (j = M.a f ∨ j = M.b f ∨ j = M.c f)

/-- The mesh is connected: any two vertices are joined by a chain of faces. -/
def Connected : Prop := ∀ i j : V, Relation.ReflTransGen M.Adj i j

/-- **Kernel of the cotan Laplacian on a connected mesh**: `L u = 0` iff `u` is constant.

[formula: cotan-laplacian-kernel] -/
theorem laplacian_mulVec_eq_zero_iff (hM : M.Valid) (hconn : M.Connected) (u : V → ℝ) :
    M.laplacian *ᵥ u = 0 ↔ ∀ i j, u i = u j := by
  constructor
  · intro h0
    have hq : ∑ f, M.faceEnergy f u = 0 := by
      rw [← M.dotProduct_laplacian_mulVec_eq_sum_faceEnergy hM, h0]; simp
    have hf : ∀ f, M.faceEnergy f u = 0 := fun f =>
      (Finset.sum_eq_zero_iff_of_nonneg
        (fun f _ => cotEnergy_nonneg (hM.tri f) _ _ _)).mp hq f (Finset.mem_univ f)
    have hadj : ∀ i j, M.Adj i j → u i = u j := by
      rintro i j ⟨f, hi, hj⟩
      obtain ⟨h1, h2⟩ := (cotEnergy_eq_zero_iff (hM.tri f) _ _ _).mp (hf f)
      rcases hi with rfl | rfl | rfl <;> rcases hj with rfl | rfl | rfl <;> simp_all
    intro i j
    induction hconn i j with
    | refl => rfl
    | tail _ hbc ih => exact ih.trans (hadj _ _ hbc)
  · intro hconst
    ext i
    rw [laplacian_mulVec_apply]
    simp only [Pi.zero_apply]
    exact Finset.sum_eq_zero fun j _ => by rw [hconst i j]; ring

/-- **Dirichlet reduction.** If `u` vanishes outside the set `p` of interior vertices, then for an
interior vertex `i` the full product `(L u)_i` is the product with the interior block
`L_II = L.submatrix`. -/
theorem laplacian_mulVec_of_zero_off (p : V → Prop) [DecidablePred p] (u : V → ℝ)
    (hu : ∀ v, ¬p v → u v = 0) (i : {v // p v}) :
    (M.laplacian *ᵥ u) i = ((M.laplacian.submatrix Subtype.val Subtype.val) *ᵥ
      (fun j : {v // p v} => u j)) i := by
  simp only [Matrix.mulVec, dotProduct, Matrix.submatrix_apply]
  rw [← Finset.sum_subtype (Finset.univ.filter p) (by simp)
    (fun j => M.laplacian i j * u j)]
  rw [Finset.sum_filter]
  refine Finset.sum_congr rfl fun j _ => ?_
  by_cases h : p j
  · simp [h]
  · simp [h, hu j h]

/-- **The interior block is positive definite** on a connected mesh with at least one boundary
(non-interior) vertex: the hypothesis of the Dirichlet case of the optimality theorems.

[formula: cotan-laplacian-interior-posdef] -/
theorem laplacian_interior_posDef (hM : M.Valid) (hconn : M.Connected) (p : V → Prop)
    [DecidablePred p] (hB : ∃ v, ¬p v) :
    (M.laplacian.submatrix (Subtype.val : {v // p v} → V) Subtype.val).PosDef := by
  refine Matrix.PosDef.of_dotProduct_mulVec_pos ?_ fun x hx => ?_
  · exact (M.laplacian_isSymm.submatrix _ : _)
  · rw [star_trivial]
    set u : V → ℝ := fun v => if h : p v then x ⟨v, h⟩ else 0 with hu
    have hux : (fun j : {v // p v} => u j) = x := by
      funext j; simp [hu, j.2]
    have hoff : ∀ v, ¬p v → u v = 0 := fun v hv => by simp [hu, hv]
    have hq : x ⬝ᵥ ((M.laplacian.submatrix Subtype.val Subtype.val) *ᵥ x) =
        u ⬝ᵥ (M.laplacian *ᵥ u) := by
      have e1 : u ⬝ᵥ (M.laplacian *ᵥ u) =
          ∑ i : {v // p v}, u i * (M.laplacian *ᵥ u) i := by
        simp only [dotProduct]
        rw [← Finset.sum_subtype (Finset.univ.filter p) (by simp)
          (fun i => u i * (M.laplacian *ᵥ u) i), Finset.sum_filter]
        refine Finset.sum_congr rfl fun i _ => ?_
        by_cases h : p i
        · simp [h]
        · simp [h, hoff i h]
      rw [e1]
      simp only [dotProduct]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [M.laplacian_mulVec_of_zero_off p u hoff i, hux]
      congr 1
      exact (congrFun hux i).symm
    rw [hq]
    have hnn := (M.laplacian_posSemidef hM).dotProduct_mulVec_nonneg u
    rw [star_trivial] at hnn
    rcases hnn.lt_or_eq with hpos | hzero
    · exact hpos
    · exfalso
      have hLu : M.laplacian *ᵥ u = 0 := by
        have := (M.laplacian_posSemidef hM).dotProduct_mulVec_zero_iff (x := u)
        rw [star_trivial] at this
        exact this.mp hzero.symm
      have hconst := (M.laplacian_mulVec_eq_zero_iff hM hconn u).mp hLu
      obtain ⟨b, hb⟩ := hB
      apply hx
      funext j
      have := hconst j b
      rw [hoff b hb] at this
      rw [← hux]
      exact this

/-! ## Lumped mass -/

/-- Sum of the areas of the faces containing the vertex `v`. -/
def vertexArea (v : V) : ℝ :=
  ∑ f, if v = M.a f ∨ v = M.b f ∨ v = M.c f then M.area f else 0

/-- The lumped mass `m_v = ⅓ Σ (area of the triangles at v)` (paper Eq. 16). -/
def mass (v : V) : ℝ := M.vertexArea v / 3

/-- Total area of the mesh. -/
def totalArea : ℝ := ∑ f, M.area f

/-- **The lumped masses sum to the total area.**

[formula: lumped-mass] -/
theorem sum_mass (hM : M.Valid) : ∑ v, M.mass v = M.totalArea := by
  unfold mass vertexArea totalArea
  rw [← Finset.sum_div, Finset.sum_comm]
  rw [Finset.sum_div]
  refine Finset.sum_congr rfl fun f _ => ?_
  have hab := hM.ab f
  have hbc := hM.bc f
  have hca := hM.ca f
  have h3 : ∀ v : V, (if v = M.a f ∨ v = M.b f ∨ v = M.c f then M.area f else 0) =
      (if v = M.a f then M.area f else 0) + (if v = M.b f then M.area f else 0) +
        (if v = M.c f then M.area f else 0) := by
    intro v
    by_cases h1 : v = M.a f <;> by_cases h2 : v = M.b f <;> by_cases h3 : v = M.c f <;> simp_all
  simp only [h3, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  ring

/-- **Area normalization**: dividing the lumped masses by the total area makes them sum to `1`.

[formula: lumped-mass] -/
theorem sum_mass_div_totalArea (hM : M.Valid) (hA : M.totalArea ≠ 0) :
    ∑ v, M.mass v / M.totalArea = 1 := by
  rw [← Finset.sum_div, M.sum_mass hM, div_self hA]

omit [DecidableEq V] [Fintype V] in
/-- The total area of a valid mesh with at least one face is positive. -/
theorem totalArea_pos (hM : M.Valid) [Nonempty F] : 0 < M.totalArea :=
  Finset.sum_pos (fun f _ => triArea_pos (hM.tri f)) Finset.univ_nonempty

/-! ## Uniform scaling -/

omit [DecidableEq V] [Fintype F] [Fintype V] in
/-- The mesh with all edge lengths multiplied by `s`. -/
def scale (s : ℝ) : TriMesh V F := { M with len := fun i j => s * M.len i j }

/-- **The cotan Laplacian is invariant under uniform scaling.**

[formula: scale-invariance] -/
theorem laplacian_scale {s : ℝ} (hs : s ≠ 0) : (M.scale s).laplacian = M.laplacian := by
  have hcot : ∀ f i j, (M.scale s).cotOpp f i j = M.cotOpp f i j := by
    intro f i j
    simp only [cotOpp, cotA, cotB, cotC, la, lb, lc, scale, cornerCot_smul hs]
  ext i j
  simp only [laplacian, offDiag, hcot, Matrix.of_apply]

omit [Fintype V] in
/-- **The lumped masses scale with the square of the scaling factor.**

[formula: scale-invariance] -/
theorem mass_scale (s : ℝ) (v : V) : (M.scale s).mass v = s ^ 2 * M.mass v := by
  simp only [mass, vertexArea, area, la, lb, lc, scale, triArea_smul]
  rw [← mul_div_assoc, Finset.mul_sum]
  congr 1
  refine Finset.sum_congr rfl fun f _ => ?_
  split_ifs <;> simp

omit [DecidableEq V] [Fintype F] [Fintype V] in
/-- The corner angles are invariant under uniform scaling. -/
theorem angleA_scale {s : ℝ} (hs : s ≠ 0) (f : F) : (M.scale s).angleA f = M.angleA f := by
  simp only [angleA, la, lb, lc, scale, cornerAngle_smul hs]

omit [DecidableEq V] [Fintype F] [Fintype V] in
theorem angleB_scale {s : ℝ} (hs : s ≠ 0) (f : F) : (M.scale s).angleB f = M.angleB f := by
  simp only [angleB, la, lb, lc, scale, cornerAngle_smul hs]

omit [DecidableEq V] [Fintype F] [Fintype V] in
theorem angleC_scale {s : ℝ} (hs : s ≠ 0) (f : F) : (M.scale s).angleC f = M.angleC f := by
  simp only [angleC, la, lb, lc, scale, cornerAngle_smul hs]

end TriMesh

end DDG
