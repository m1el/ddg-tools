import DdgProofs.Intrinsic.Barycentric

/-!
# Intrinsic simplification: anisotropy, Galerkin coarsening, distortion, quadrics

`crates/ddg-intrinsic-simplification/SPEC.md` §8, §9.2, §10.3, §10.4, §13; the paper Sec. 8.3,
8.4, 8.6.2, 8.7; Khodakovsky, Litke, Schröder 2003, Sec. 2; Garland, Heckbert 1997, Sec. 5.

* **`iem-anisotropic-scale`** — the factor `(1 − τ) + (τ/2)((u_i·ê_ij)² + (u_j·ê_ji)²)` lies in
  `[1 − τ, 1]` for `τ ∈ [0, 1]` and unit vectors (`anisoFactor_mem`, any real inner product
  space), and is positive for `τ < 1`.
* **`iem-galerkin-spd`** — `Pᵀ A P` is positive definite for positive definite `A` and `P` of full
  column rank (`galerkin_posDef`); a prolongation in which every coarse vertex has the row of a
  fine vertex equal to its unit vector has full column rank (`injective_of_unit_rows`).
* **`iem-distortion`** — the affine map between two plane triangles (`p` counter-clockwise, with
  the edge vectors `e₁ = p₁ − p₀`, `e₂ = p₂ − p₀`, and any `q`): its linear part `J` (`affJ`) maps
  `e₁ ↦ q₁ − q₀`, `e₂ ↦ q₂ − q₀`; `det J = [q]/[p]` (signed areas), so `|det J|` is the area ratio
  `A(ℓ̂)/A(ℓ)` of the Heron areas; `‖J‖²_F = Σ_k cot θ_k ℓ̂_k² / (2A)` with the corner cotangents of
  the input lengths (`affJ_frobenius`): both from the lengths only. The two singular values
  `σ₁ ≥ σ₂ ≥ 0` have `σ₁² + σ₂² = ‖J‖²_F`, `σ₁σ₂ = |det J|` and are the square roots of
  `(s ± √(s² − 4d²))/2` (`singularValues`): these are the eigenvalues `γ` of the metric tensor
  `g = JᵀJ` (their sum is `tr g`, their product `det g`). Khodakovsky et al. define the
  anisotropic and area distortion with `γ_max/γ_min = (σ₁/σ₂)²` and `γ_max γ_min = (σ₁σ₂)²`;
  SPEC §10.3 uses `σ₁/σ₂` and `log(σ₁σ₂)` (the same orderings, different scales).
* **`iem-qem-quadric`** — the quadric `(n nᵀ, d n, d²)` of the plane `n·x + d = 0` (`‖n‖ = 1`)
  evaluates to `(n·v + d)²`, the squared distance to the plane: attained at the foot
  `v − (n·v + d) n` and a lower bound for every point of the plane (`planeQuadric_eval`,
  `planeQuadric_dist`); quadrics add (`evalQ_add`); a solution of the 3×3 system `A v = −b` of a
  quadric with symmetric positive semidefinite `A` minimizes it, uniquely if `A` is positive
  definite (`evalQ_optimal`).
-/

noncomputable section

namespace DdgProofs

namespace IntrinsicSimplification

open Real Matrix

/-! ## `iem-anisotropic-scale` -/

/-- The length factor of Sec. 8.6.2, with `x = u_i·ê_ij`, `y = u_j·ê_ji`. -/
def anisoFactor (τ x y : ℝ) : ℝ := (1 - τ) + τ / 2 * (x ^ 2 + y ^ 2)

/-- **The factor lies in `[1 − τ, 1]`** for `τ ∈ [0, 1]` and `|x|, |y| ≤ 1`.

[formula: iem-anisotropic-scale] -/
theorem anisoFactor_bounds {τ x y : ℝ} (h0 : 0 ≤ τ) (_h1 : τ ≤ 1) (hx : |x| ≤ 1) (hy : |y| ≤ 1) :
    1 - τ ≤ anisoFactor τ x y ∧ anisoFactor τ x y ≤ 1 := by
  have hx2 : x ^ 2 ≤ 1 := by rw [← sq_abs]; nlinarith [abs_nonneg x]
  have hy2 : y ^ 2 ≤ 1 := by rw [← sq_abs]; nlinarith [abs_nonneg y]
  unfold anisoFactor
  constructor <;> nlinarith [sq_nonneg x, sq_nonneg y]

/-- **For unit vectors** `u_i, ê_ij, u_j, ê_ji` of a real inner product space the factor lies in
`[1 − τ, 1]`, and it is positive when `τ < 1`.

[formula: iem-anisotropic-scale] -/
theorem anisoFactor_mem {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] {τ : ℝ}
    (h0 : 0 ≤ τ) (h1 : τ ≤ 1) {ui eij uj eji : E} (hui : ‖ui‖ = 1) (heij : ‖eij‖ = 1)
    (huj : ‖uj‖ = 1) (heji : ‖eji‖ = 1) :
    1 - τ ≤ anisoFactor τ (inner ℝ ui eij) (inner ℝ uj eji) ∧
      anisoFactor τ (inner ℝ ui eij) (inner ℝ uj eji) ≤ 1 ∧
      (τ < 1 → 0 < anisoFactor τ (inner ℝ ui eij) (inner ℝ uj eji)) := by
  have hx : |inner ℝ ui eij| ≤ 1 := by
    have := abs_real_inner_le_norm ui eij; rw [hui, heij, one_mul] at this; exact this
  have hy : |inner ℝ uj eji| ≤ 1 := by
    have := abs_real_inner_le_norm uj eji; rw [huj, heji, one_mul] at this; exact this
  obtain ⟨b1, b2⟩ := anisoFactor_bounds h0 h1 hx hy
  exact ⟨b1, b2, fun h => by linarith⟩

/-! ## `iem-galerkin-spd` -/

/-- **`Pᵀ A P` is SPD** for SPD `A` and `P` of full column rank.

[formula: iem-galerkin-spd] -/
theorem galerkin_posDef {n m : Type*} [Fintype n] [Fintype m] {A : Matrix n n ℝ} (hA : A.PosDef)
    {P : Matrix n m ℝ} (hP : Function.Injective P.mulVec) : (Pᵀ * A * P).PosDef := by
  have := hA.conjTranspose_mul_mul_same hP
  rwa [conjTranspose_eq_transpose_of_trivial] at this

/-- **The rows of the coarse vertices give full column rank**: if every column `c` (coarse vertex)
has a row (the fine vertex that is `c`) equal to the unit vector `e_c`, then `P` is injective.

[formula: iem-galerkin-spd] -/
theorem injective_of_unit_rows {n m : Type*} [Fintype m] [DecidableEq m] {P : Matrix n m ℝ}
    (h : ∀ c, ∃ r, ∀ c', P r c' = if c' = c then 1 else 0) : Function.Injective P.mulVec := by
  intro x y hxy
  funext c
  obtain ⟨r, hr⟩ := h c
  have := congrFun hxy r
  simp only [mulVec, dotProduct, hr, ite_mul, one_mul, zero_mul, Finset.sum_ite_eq',
    Finset.mem_univ, ite_true] at this
  exact this

/-! ## `iem-distortion` -/

/-- The linear part of the affine map taking the triangle `p₀p₁p₂` to `q₀q₁q₂`: `J = F E⁻¹` with
the edge matrices `E = [p₁ − p₀, p₂ − p₀]`, `F = [q₁ − q₀, q₂ − q₀]`. Entries `J₁₁ J₁₂ J₂₁ J₂₂`. -/
def affJ (p₀ p₁ p₂ q₀ q₁ q₂ : ℝ × ℝ) : ℝ × ℝ × ℝ × ℝ :=
  let e1 := p₁ - p₀
  let e2 := p₂ - p₀
  let f1 := q₁ - q₀
  let f2 := q₂ - q₀
  let D := sda p₀ p₁ p₂
  ((f1.1 * e2.2 - f2.1 * e1.2) / D, (f2.1 * e1.1 - f1.1 * e2.1) / D,
    (f1.2 * e2.2 - f2.2 * e1.2) / D, (f2.2 * e1.1 - f1.2 * e2.1) / D)

/-- Apply a 2×2 matrix given by its entries. -/
def apply2 (J : ℝ × ℝ × ℝ × ℝ) (v : ℝ × ℝ) : ℝ × ℝ :=
  (J.1 * v.1 + J.2.1 * v.2, J.2.2.1 * v.1 + J.2.2.2 * v.2)

/-- Determinant of a 2×2 matrix given by its entries. -/
def det2 (J : ℝ × ℝ × ℝ × ℝ) : ℝ := J.1 * J.2.2.2 - J.2.1 * J.2.2.1

/-- Squared Frobenius norm. -/
def frob2 (J : ℝ × ℝ × ℝ × ℝ) : ℝ := J.1 ^ 2 + J.2.1 ^ 2 + J.2.2.1 ^ 2 + J.2.2.2 ^ 2

variable {p₀ p₁ p₂ : ℝ × ℝ} (q₀ q₁ q₂ : ℝ × ℝ)

/-- `J` maps the edges of the first triangle to those of the second.

[formula: iem-distortion] -/
theorem affJ_maps (h : sda p₀ p₁ p₂ ≠ 0) :
    apply2 (affJ p₀ p₁ p₂ q₀ q₁ q₂) (p₁ - p₀) = q₁ - q₀ ∧
      apply2 (affJ p₀ p₁ p₂ q₀ q₁ q₂) (p₂ - p₀) = q₂ - q₀ := by
  obtain ⟨x0, y0⟩ := p₀
  obtain ⟨x1, y1⟩ := p₁
  obtain ⟨x2, y2⟩ := p₂
  unfold apply2 affJ sda at *
  simp only [Prod.mk_sub_mk, Prod.ext_iff] at h ⊢
  generalize hD : (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0) = D at h ⊢
  refine ⟨⟨?_, ?_⟩, ?_, ?_⟩ <;> (field_simp; subst hD; ring)

/-- **`det J` is the ratio of the signed areas**; with `p` counter-clockwise, `|det J|` is the
ratio of the Heron areas of the two triangles' lengths.

[formula: iem-distortion] -/
theorem affJ_det (h : 0 < sda p₀ p₁ p₂) :
    det2 (affJ p₀ p₁ p₂ q₀ q₁ q₂) = sda q₀ q₁ q₂ / sda p₀ p₁ p₂ ∧
      |det2 (affJ p₀ p₁ p₂ q₀ q₁ q₂)| =
        triArea (dist2 q₁ q₂) (dist2 q₂ q₀) (dist2 q₀ q₁) /
          triArea (dist2 p₁ p₂) (dist2 p₂ p₀) (dist2 p₀ p₁) := by
  have hd : det2 (affJ p₀ p₁ p₂ q₀ q₁ q₂) = sda q₀ q₁ q₂ / sda p₀ p₁ p₂ := by
    obtain ⟨x0, y0⟩ := p₀
    obtain ⟨x1, y1⟩ := p₁
    obtain ⟨x2, y2⟩ := p₂
    obtain ⟨a0, b0⟩ := q₀
    obtain ⟨a1, b1⟩ := q₁
    obtain ⟨a2, b2⟩ := q₂
    unfold det2 affJ
    unfold sda at *
    simp only [Prod.mk_sub_mk] at h ⊢
    generalize hD : (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0) = D at h ⊢
    have hD0 : D ≠ 0 := h.ne'
    field_simp
    subst hD
    ring
  refine ⟨hd, ?_⟩
  rw [hd, abs_div, abs_of_pos h]
  have hq := abs_sda_eq_two_mul_triArea q₀ q₁ q₂
  have hp := abs_sda_eq_two_mul_triArea p₀ p₁ p₂
  rw [abs_of_pos h] at hp
  rw [hq, hp]
  have : 0 < triArea (dist2 p₁ p₂) (dist2 p₂ p₀) (dist2 p₀ p₁) := by linarith
  field_simp

/-- **`‖J‖²_F` from the lengths**: `(cot θ₀ ℓ̂₁₂² + cot θ₁ ℓ̂₂₀² + cot θ₂ ℓ̂₀₁²)/(2A)`, with the
corner cotangents and the Heron area `A` of the first triangle's lengths and the lengths `ℓ̂` of
the second.

[formula: iem-distortion] -/
theorem affJ_frobenius (h : 0 < sda p₀ p₁ p₂) :
    let l12 := dist2 p₁ p₂
    let l20 := dist2 p₂ p₀
    let l01 := dist2 p₀ p₁
    frob2 (affJ p₀ p₁ p₂ q₀ q₁ q₂) =
      (cornerCot l12 l20 l01 * dist2 q₁ q₂ ^ 2 + cornerCot l20 l01 l12 * dist2 q₂ q₀ ^ 2 +
        cornerCot l01 l12 l20 * dist2 q₀ q₁ ^ 2) / (2 * triArea l12 l20 l01) := by
  intro l12 l20 l01
  have hp := abs_sda_eq_two_mul_triArea p₀ p₁ p₂
  rw [abs_of_pos h] at hp
  have hA : triArea l12 l20 l01 = sda p₀ p₁ p₂ / 2 := by simp only [l12, l20, l01]; linarith
  have r1 : triArea l20 l01 l12 = triArea l12 l20 l01 := triArea_rotate l12 l20 l01
  have r2 : triArea l01 l12 l20 = triArea l12 l20 l01 := by
    rw [triArea_rotate l20 l01 l12, r1]
  unfold cornerCot
  rw [r1, r2, hA]
  simp only [l12, l20, l01, dist2_sq, sqDist2]
  have hs : sda p₀ p₁ p₂ ≠ 0 := h.ne'
  unfold frob2 affJ
  simp only [Prod.fst_sub, Prod.snd_sub]
  unfold sda at *
  field_simp
  ring

/-- **The singular values.** For any 2×2 matrix with `s = ‖J‖²_F`, `d = det J`: `s² − 4d² ≥ 0`, and
`γ± = (s ± √(s² − 4d²))/2` satisfy `γ₊ + γ₋ = s`, `γ₊ γ₋ = d²`, `0 ≤ γ₋ ≤ γ₊` — the eigenvalues of
`JᵀJ` (trace `s`, determinant `d²`); `σ = √γ` has `σ₁² + σ₂² = s` and `σ₁σ₂ = |d|`.

[formula: iem-distortion] -/
theorem singularValues (J : ℝ × ℝ × ℝ × ℝ) :
    let s := frob2 J
    let d := det2 J
    let γp := (s + √(s ^ 2 - 4 * d ^ 2)) / 2
    let γm := (s - √(s ^ 2 - 4 * d ^ 2)) / 2
    0 ≤ s ^ 2 - 4 * d ^ 2 ∧ γp + γm = s ∧ γp * γm = d ^ 2 ∧ 0 ≤ γm ∧ γm ≤ γp ∧
      √γp ^ 2 + √γm ^ 2 = s ∧ √γp * √γm = |d| := by
  intro s d γp γm
  have hdisc : 0 ≤ s ^ 2 - 4 * d ^ 2 := by
    have : s ^ 2 - 4 * d ^ 2 = ((J.1 - J.2.2.2) ^ 2 + (J.2.1 + J.2.2.1) ^ 2) *
        ((J.1 + J.2.2.2) ^ 2 + (J.2.1 - J.2.2.1) ^ 2) := by
      simp only [s, d, frob2, det2]; ring
    rw [this]; positivity
  have hsq := sq_sqrt hdisc
  have hs0 : 0 ≤ s := by simp only [s, frob2]; positivity
  have hroot : √(s ^ 2 - 4 * d ^ 2) ≤ s := by
    rw [sqrt_le_left hs0]; nlinarith [sq_nonneg d]
  have hsum : γp + γm = s := by simp only [γp, γm]; ring
  have hprod : γp * γm = d ^ 2 := by
    simp only [γp, γm]; nlinarith [hsq]
  have hm : 0 ≤ γm := by simp only [γm]; linarith
  have hmp : γm ≤ γp := by simp only [γm, γp]; linarith [sqrt_nonneg (s ^ 2 - 4 * d ^ 2)]
  have hp : 0 ≤ γp := le_trans hm hmp
  refine ⟨hdisc, hsum, hprod, hm, hmp, ?_, ?_⟩
  · rw [sq_sqrt hp, sq_sqrt hm, hsum]
  · rw [← sqrt_mul hp, hprod, sqrt_sq_eq_abs]

/-! ## `iem-qem-quadric` -/

/-- A quadric `v ↦ vᵀAv + 2bᵀv + c` on `ℝ³`. -/
structure Quadric where
  /-- the matrix -/
  A : Matrix (Fin 3) (Fin 3) ℝ
  /-- the linear part -/
  b : Fin 3 → ℝ
  /-- the constant -/
  c : ℝ

namespace Quadric

/-- Evaluate a quadric. -/
def eval (Q : Quadric) (v : Fin 3 → ℝ) : ℝ := v ⬝ᵥ (Q.A *ᵥ v) + 2 * (Q.b ⬝ᵥ v) + Q.c

/-- The sum of two quadrics. -/
def add (Q R : Quadric) : Quadric := ⟨Q.A + R.A, Q.b + R.b, Q.c + R.c⟩

/-- Quadrics add: the quadric of a vertex is the sum over its planes.

[formula: iem-qem-quadric] -/
theorem evalQ_add (Q R : Quadric) (v : Fin 3 → ℝ) : (Q.add R).eval v = Q.eval v + R.eval v := by
  simp only [eval, add, add_mulVec, dotProduct_add, add_dotProduct]
  ring

/-- The fundamental quadric of the plane `n·x + d = 0` (Garland–Heckbert's `K_p = ppᵀ`). -/
def plane (n : Fin 3 → ℝ) (d : ℝ) : Quadric := ⟨vecMulVec n n, d • n, d ^ 2⟩

/-- **The plane quadric evaluates to `(n·v + d)²`.**

[formula: iem-qem-quadric] -/
theorem planeQuadric_eval (n : Fin 3 → ℝ) (d : ℝ) (v : Fin 3 → ℝ) :
    (plane n d).eval v = (n ⬝ᵥ v + d) ^ 2 := by
  simp only [eval, plane, vecMulVec_mulVec, dotProduct_smul, smul_dotProduct, smul_eq_mul,
    op_smul_eq_mul]
  rw [dotProduct_comm v n]
  ring

/-- **It is the squared distance to the plane** (`‖n‖ = 1`): the foot `v − (n·v + d) n` lies on the
plane at the squared distance `(n·v + d)²`, and every point of the plane is at least that far.

[formula: iem-qem-quadric] -/
theorem planeQuadric_dist {n : Fin 3 → ℝ} (hn : n ⬝ᵥ n = 1) (d : ℝ) (v : Fin 3 → ℝ) :
    let foot := v - (n ⬝ᵥ v + d) • n
    n ⬝ᵥ foot + d = 0 ∧ (v - foot) ⬝ᵥ (v - foot) = (plane n d).eval v ∧
      ∀ x, n ⬝ᵥ x + d = 0 → (plane n d).eval v ≤ (v - x) ⬝ᵥ (v - x) := by
  intro foot
  rw [planeQuadric_eval]
  refine ⟨?_, ?_, fun x hx => ?_⟩
  · simp only [foot, dotProduct_sub, dotProduct_smul, smul_eq_mul, hn]; ring
  · simp only [foot, sub_sub_cancel, smul_dotProduct, dotProduct_smul, smul_eq_mul, hn]; ring
  · -- `n·(v − x) = n·v + d`, and Cauchy–Schwarz with `‖n‖ = 1`
    have hnv : n ⬝ᵥ (v - x) = n ⬝ᵥ v + d := by rw [dotProduct_sub]; linarith
    have cs : (n ⬝ᵥ (v - x)) ^ 2 ≤ (n ⬝ᵥ n) * ((v - x) ⬝ᵥ (v - x)) := by
      set w := v - x
      have key : 0 ≤ (w - (n ⬝ᵥ w) • n) ⬝ᵥ (w - (n ⬝ᵥ w) • n) := by
        simp only [dotProduct, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
        exact Finset.sum_nonneg (fun i _ => mul_self_nonneg _)
      simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, smul_eq_mul,
        hn] at key
      rw [dotProduct_comm w n] at key
      rw [hn]
      nlinarith
    rw [hn, one_mul, hnv] at cs
    exact cs

/-- **The optimal point solves the 3×3 system.** For a quadric with symmetric `A`, if `A v* = −b`
then `Q(v) = Q(v*) + (v − v*)ᵀ A (v − v*)`; hence `v*` minimizes `Q` when `A` is positive
semidefinite, and is the unique minimizer when `A` is positive definite.

[formula: iem-qem-quadric] -/
theorem evalQ_optimal (Q : Quadric) (hsym : Q.Aᵀ = Q.A) {vs : Fin 3 → ℝ}
    (hv : Q.A *ᵥ vs = -Q.b) (v : Fin 3 → ℝ) :
    Q.eval v = Q.eval vs + (v - vs) ⬝ᵥ (Q.A *ᵥ (v - vs)) ∧
      (Q.A.PosSemidef → Q.eval vs ≤ Q.eval v) ∧ (Q.A.PosDef → v ≠ vs → Q.eval vs < Q.eval v) := by
  have hb : Q.b = -(Q.A *ᵥ vs) := by rw [hv, neg_neg]
  have hsymm : ∀ x y : Fin 3 → ℝ, x ⬝ᵥ (Q.A *ᵥ y) = y ⬝ᵥ (Q.A *ᵥ x) := fun x y => by
    rw [dotProduct_mulVec, ← mulVec_transpose, hsym, dotProduct_comm]
  have hexp : Q.eval v = Q.eval vs + (v - vs) ⬝ᵥ (Q.A *ᵥ (v - vs)) := by
    simp only [eval, hb, mulVec_sub, dotProduct_sub, sub_dotProduct, neg_dotProduct]
    rw [dotProduct_comm (Q.A *ᵥ vs) v, dotProduct_comm (Q.A *ᵥ vs) vs, hsymm v vs]
    ring
  refine ⟨hexp, fun hA => ?_, fun hA hne => ?_⟩
  · have := hA.dotProduct_mulVec_nonneg (v - vs)
    simp only [star_trivial] at this
    linarith
  · have := hA.dotProduct_mulVec_pos (sub_ne_zero.mpr hne)
    simp only [star_trivial] at this
    linarith

end Quadric

end IntrinsicSimplification

end DdgProofs
