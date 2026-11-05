import DDG.Geometry.Vec3

/-!
# Repulsive Surfaces: exponents, face quantities, dual areas

"Repulsive Surfaces", Yu, Brakensiek, Schumacher, Crane, ACM TOG 40(6), Article 268, 2021
("the paper"); contract `crates/ddg-repulsive-surfaces/SPEC.md` ("SPEC"), §2.1, §2.2.

Vectors are `V3 = Fin 3 → ℝ` with Mathlib's `dotProduct` (`⬝ᵥ`) and `crossProduct` (`⨯₃`); `vnorm`
is the Euclidean norm (all from `DDG.RC`, the curves package). A triangle mesh has vertices
`V`, faces `F` and corners `fv : F → Fin 3 → V` (`S = (i, j, k) = (fv S 0, fv S 1, fv S 2)`, the
caller's order; nothing here needs the mesh to be manifold or oriented). Positions are
`f : V → V3`. Per face `S`:

```
n_S = (f_j − f_i) × (f_k − f_i),   a_S = ½|n_S|,   N_S = n_S/|n_S|,   X_S = (f_i + f_j + f_k)/3
e_a = f_{a+2} − f_{a+1}            (the edge opposite corner a: e_jk, e_ki, e_ij)
a_v = ⅓ Σ_{S ∋ v} a_S
```

* `rs-exponents` (SPEC §2.2, D3): `s = (β − 2)/α`, `s = 2 − 2/p` for `α = p`, `β = 2p`;
  `2σ + 2 = 2s`; `2(2 − s) + 2 = 6 − 2s`; `s ∈ (1, 2) ⇔ α + 2 < β < 2α + 2` (for `α > 0`).
* `rs-face-quantities`: `a_S = ½|n_S|`, `|N_S| = 1`, `n_S = 2 a_S N_S`, `N_S·n_S = 2a_S`,
  `N_S ⊥ e_a`, `Σ_a e_a = 0`, `n_S` is the same from every corner, `Σ_v a_v = Σ_S a_S`.
-/

noncomputable section

namespace DDG

namespace RS

open Matrix Real Finset
open RC (V3 vnorm vnorm_nonneg vnorm_sq vnorm_mul_self vnorm_eq_zero_iff vnorm_pos_iff
  dot_self_nonneg unit_dot_self)

/-! ## Exponents -/

/-- The order parameter `s = (β − 2)/α` of the energy (SPEC §1, §2.2).

[formula: rs-exponents] -/
def sOrder (α β : ℝ) : ℝ := (β - 2) / α

/-- **The exponent bookkeeping of SPEC §2.2** (D3): for `α = p`, `β = 2p`, `s = 2 − 2/p`; with
`σ = s − 1`, the high-order kernel exponent `2σ + 2` is `2s`; the preconditioner's kernel exponent
`2(2 − s) + 2` is `6 − 2s`.

[formula: rs-exponents] -/
theorem exponents_bookkeeping (α β : ℝ) {p : ℝ} (hp : p ≠ 0) :
    sOrder p (2 * p) = 2 - 2 / p ∧
      2 * (sOrder α β - 1) + 2 = 2 * sOrder α β ∧
      2 * (2 - sOrder α β) + 2 = 6 - 2 * sOrder α β := by
  refine ⟨?_, by ring, by ring⟩
  unfold sOrder
  field_simp

/-- **The admissible range** (SPEC §2.2): for `α > 0`, `s = (β − 2)/α ∈ (1, 2)` iff
`α + 2 < β < 2α + 2`; then `σ = s − 1 ∈ (0, 1)`.

[formula: rs-exponents] -/
theorem exponents_range {α : ℝ} (hα : 0 < α) (β : ℝ) :
    (1 < sOrder α β ∧ sOrder α β < 2 ↔ α + 2 < β ∧ β < 2 * α + 2) ∧
      (1 < sOrder α β ∧ sOrder α β < 2 → 0 < sOrder α β - 1 ∧ sOrder α β - 1 < 1) := by
  unfold sOrder
  refine ⟨?_, fun h => ⟨by linarith [h.1], by linarith [h.2]⟩⟩
  rw [one_lt_div hα, div_lt_iff₀ hα]
  constructor <;> rintro ⟨h1, h2⟩ <;> constructor <;> linarith

/-- The paper's default `p = 6` (`α = 6`, `β = 12`): `s = 5/3`, critical exponent `β − α = 6 > 4`.

[formula: rs-exponents] -/
theorem exponents_default : sOrder 6 12 = 5 / 3 ∧ (12 : ℝ) - 6 > 4 := by
  unfold sOrder; norm_num

/-! ## Face quantities -/

section faces

variable {V F : Type*}

/-- The corner `a` of face `S`: `f_{fv S a}`. -/
def corner (f : V → V3) (fv : F → Fin 3 → V) (S : F) (a : Fin 3) : V3 := f (fv S a)

/-- The edge opposite corner `a`: `e_a = f_{a+2} − f_{a+1}` (`e_jk`, `e_ki`, `e_ij`; not
normalized, SPEC D2). -/
def opp (f : V → V3) (fv : F → Fin 3 → V) (S : F) (a : Fin 3) : V3 :=
  corner f fv S (a + 2) - corner f fv S (a + 1)

/-- The area vector `n_S = (f_j − f_i) × (f_k − f_i)`. -/
def nvec (f : V → V3) (fv : F → Fin 3 → V) (S : F) : V3 :=
  (corner f fv S 1 - corner f fv S 0) ⨯₃ (corner f fv S 2 - corner f fv S 0)

/-- The face area `a_S = ½|n_S|`.

[formula: rs-face-quantities] -/
def area (f : V → V3) (fv : F → Fin 3 → V) (S : F) : ℝ := vnorm (nvec f fv S) / 2

/-- The unit normal `N_S = n_S/|n_S|` (`0` for a degenerate face).

[formula: rs-face-quantities] -/
def normal (f : V → V3) (fv : F → Fin 3 → V) (S : F) : V3 :=
  (vnorm (nvec f fv S))⁻¹ • nvec f fv S

/-- The barycenter `X_S = (f_i + f_j + f_k)/3`.

[formula: rs-face-quantities] -/
def bary (f : V → V3) (fv : F → Fin 3 → V) (S : F) : V3 :=
  (1 / 3 : ℝ) • (corner f fv S 0 + corner f fv S 1 + corner f fv S 2)

/-- A face is nondegenerate when `n_S ≠ 0` (positive area). -/
def Nondeg (f : V → V3) (fv : F → Fin 3 → V) (S : F) : Prop := nvec f fv S ≠ 0

/-- Componentwise expansion of vector identities in `ℝ³`. -/
theorem v3_ext {u w : V3} (h0 : u 0 = w 0) (h1 : u 1 = w 1) (h2 : u 2 = w 2) : u = w := by
  ext i; fin_cases i <;> assumption

@[simp] theorem fin3_add_two_zero : (0 : Fin 3) + 2 = 2 := rfl
@[simp] theorem fin3_add_one_zero : (0 : Fin 3) + 1 = 1 := rfl
@[simp] theorem fin3_add_two_one : (1 : Fin 3) + 2 = 0 := rfl
@[simp] theorem fin3_add_one_one : (1 : Fin 3) + 1 = 2 := rfl
@[simp] theorem fin3_add_two_two : (2 : Fin 3) + 2 = 1 := rfl
@[simp] theorem fin3_add_one_two : (2 : Fin 3) + 1 = 0 := rfl

/-- `n_S` as the sum of the cross products of consecutive corners. -/
theorem nvec_eq_sum (f : V → V3) (fv : F → Fin 3 → V) (S : F) :
    nvec f fv S = corner f fv S 0 ⨯₃ corner f fv S 1 + corner f fv S 1 ⨯₃ corner f fv S 2
      + corner f fv S 2 ⨯₃ corner f fv S 0 := by
  apply v3_ext <;> simp [nvec, cross_apply] <;> ring

/-- `n_S` is the same from every corner: `n_S = e_{a+2} × (−e_{a+1})`, i.e.
`(f_{a+1} − f_a) × (f_{a+2} − f_a)`. -/
theorem nvec_from_corner (f : V → V3) (fv : F → Fin 3 → V) (S : F) (a : Fin 3) :
    nvec f fv S = (corner f fv S (a + 1) - corner f fv S a) ⨯₃
      (corner f fv S (a + 2) - corner f fv S a) := by
  fin_cases a <;> apply v3_ext <;> simp [nvec, cross_apply] <;> ring

/-- The three opposite edges sum to zero. -/
theorem sum_opp (f : V → V3) (fv : F → Fin 3 → V) (S : F) : ∑ a, opp f fv S a = 0 := by
  simp only [Fin.sum_univ_three, opp]
  simp only [fin3_add_two_zero, fin3_add_one_zero, fin3_add_two_one, fin3_add_one_one,
    fin3_add_two_two, fin3_add_one_two]
  abel

/-- `n_S ⊥ e_a`. -/
theorem nvec_dot_opp (f : V → V3) (fv : F → Fin 3 → V) (S : F) (a : Fin 3) :
    nvec f fv S ⬝ᵥ opp f fv S a = 0 := by
  fin_cases a <;> simp [nvec, opp, cross_apply, dotProduct, Fin.sum_univ_three] <;> ring

theorem area_nonneg (f : V → V3) (fv : F → Fin 3 → V) (S : F) : 0 ≤ area f fv S :=
  div_nonneg (vnorm_nonneg _) (by norm_num)

theorem area_pos {f : V → V3} {fv : F → Fin 3 → V} {S : F} (h : Nondeg f fv S) :
    0 < area f fv S :=
  div_pos (vnorm_pos_iff.mpr h) (by norm_num)

/-- **The face quantities of SPEC §2.1**: for a nondegenerate face, `|N_S| = 1`, `n_S = 2a_S N_S`,
`N_S · n_S = 2a_S`, `N_S ⊥ e_a` for every corner `a`, and `Σ_a e_a = 0`.

[formula: rs-face-quantities] -/
theorem face_quantities {f : V → V3} {fv : F → Fin 3 → V} {S : F} (h : Nondeg f fv S) :
    normal f fv S ⬝ᵥ normal f fv S = 1 ∧
      nvec f fv S = (2 * area f fv S) • normal f fv S ∧
      normal f fv S ⬝ᵥ nvec f fv S = 2 * area f fv S ∧
      (∀ a, normal f fv S ⬝ᵥ opp f fv S a = 0) ∧
      ∑ a, opp f fv S a = 0 := by
  have hpos := vnorm_pos_iff.mpr h
  refine ⟨unit_dot_self h, ?_, ?_, fun a => ?_, sum_opp f fv S⟩
  · rw [normal, smul_smul, area]
    field_simp
    simp
  · rw [normal, smul_dotProduct, ← vnorm_mul_self, area, smul_eq_mul]
    field_simp
  · rw [normal, smul_dotProduct, nvec_dot_opp, smul_zero]

/-! ## Dual (barycentric) areas -/

variable [Fintype V] [Fintype F] [DecidableEq V]

/-- The barycentric dual area `a_v = ⅓ Σ_{S ∋ v} a_S` (a face is counted once per corner at `v`).

[formula: rs-face-quantities] -/
def dualArea (f : V → V3) (fv : F → Fin 3 → V) (v : V) : ℝ :=
  (1 / 3) * ∑ S, ∑ a : Fin 3, if fv S a = v then area f fv S else 0

omit [Fintype F] in
/-- `Σ_v Σ_a [fv S a = v] c = 3c`: every face has three corners. -/
theorem sum_corners_const (fv : F → Fin 3 → V) (S : F) (c : ℝ) :
    ∑ v : V, ∑ a : Fin 3, (if fv S a = v then c else 0) = 3 * c := by
  rw [Finset.sum_comm]
  simp only [Finset.sum_ite_eq, Finset.mem_univ, ite_true, Fin.sum_univ_three]
  ring

/-- **The dual areas sum to the total area**: `Σ_v a_v = Σ_S a_S`.

[formula: rs-face-quantities] -/
theorem sum_dualArea (f : V → V3) (fv : F → Fin 3 → V) :
    ∑ v, dualArea f fv v = ∑ S, area f fv S := by
  simp only [dualArea, ← Finset.mul_sum]
  rw [Finset.sum_comm]
  simp only [sum_corners_const, ← Finset.mul_sum]
  ring

omit [Fintype V] in
theorem dualArea_nonneg (f : V → V3) (fv : F → Fin 3 → V) (v : V) : 0 ≤ dualArea f fv v := by
  unfold dualArea
  refine mul_nonneg (by norm_num) (Finset.sum_nonneg fun S _ => Finset.sum_nonneg fun a _ => ?_)
  split_ifs
  · exact area_nonneg f fv S
  · exact le_rfl

end faces

end RS

end DDG
