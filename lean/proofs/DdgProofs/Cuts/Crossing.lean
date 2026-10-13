import DdgProofs.Oracle.Scalar
import DdgProofs.Cotan

/-!
# Cut flow: crossings, lengths and areas inside a face

*Variational Surface Cutting* (Sharp, Crane, ACM TOG 2018); the contract is
`crates/ddg-variational-cuts/SPEC.md`, §3.2 and §3.3. All statements are about the oracle's
scalar-generic definitions (`lean/oracle/DdgOracle/Formulas/Cuts.lean`) instantiated at `ℝ`.

* **`cut-edge-crossing`** (`Ddg.crossingParam`, `Ddg.clampCrossing`, `Ddg.crossingMoved`):
  `t = φ_a/(φ_a − φ_b)` is the only zero of the linear interpolant `(1 − t)φ_a + tφ_b`
  (`interp_crossingParam`, `interp_eq_zero_iff`); `0 ≤ t < 1` if `φ_a ≤ 0 < φ_b` and `0 < t ≤ 1`
  if `φ_b ≤ 0 < φ_a` (`crossingParam_mem_Ico`, `crossingParam_mem_Ioc`); the clamp lies in
  `[δ, 1 − δ]` for `δ ≤ ½`, is the identity exactly on that interval, moves a parameter of
  `[0, 1]` by at most `δ`, and `crossingMoved` says exactly when it moved (`clampCrossing_mem`,
  `clampCrossing_eq_self_iff`, `abs_clampCrossing_sub_le`, `crossingMoved_eq_true_iff`).
* **`cut-barycentric-length`** (`Ddg.baryLength2`, `Ddg.baryLength`): for points `P_a, P_b, P_c`
  of a real inner product space and a barycentric difference `d` with `Σ d = 0`,
  `‖d_a P_a + d_b P_b + d_c P_c‖² = −(d_a d_b ℓ_ab² + d_b d_c ℓ_bc² + d_c d_a ℓ_ca²)`
  (`baryLength2_eq_norm_sq`, `baryLength_eq_dist`; any number of points:
  `norm_sq_sum_smul_of_sum_eq_zero`). For side lengths that satisfy the triangle inequalities the
  expression is nonnegative, so the floor at `0` in `baryLength` is inactive in exact arithmetic
  (`baryLength2_nonneg`).
* **`cut-subtriangle-area`** (`Ddg.det3`, `Ddg.subTriangleArea`): the Heron area of the triangle
  whose three sides are computed by `baryLength` from three barycentric rows is
  `|det(rows)| · A_f` (`subTriangleArea_eq_triArea`; for points: `triArea_bary_points`). The
  three sub-triangles of a cut face have nonnegative determinants that sum to `1`, for either
  diagonal of the quad (`cutFace_det_sum_pc`, `cutFace_det_sum_bq`).
-/

noncomputable section

namespace DdgProofs

namespace Cuts

open scoped DdgProofs.Oracle
open Real

/-! ## Crossing of an edge by the zero level set -/

/-- The oracle's crossing parameter at `ℝ`. -/
theorem crossingParam_eq (φa φb : ℝ) : Ddg.crossingParam φa φb = φa / (φa - φb) := rfl

/-- **The crossing parameter is a zero of the linear interpolant**: with
`t = φ_a/(φ_a − φ_b)`, `(1 − t)·φ_a + t·φ_b = 0`. The cut vertex is the point
`(1 − t)·a + t·b` of the edge.

[formula: cut-edge-crossing] -/
theorem interp_crossingParam {φa φb : ℝ} (h : φa ≠ φb) :
    (1 - Ddg.crossingParam φa φb) * φa + Ddg.crossingParam φa φb * φb = 0 := by
  have h' : φa - φb ≠ 0 := sub_ne_zero.mpr h
  rw [crossingParam_eq]
  field_simp
  ring

/-- **… and the only one.**

[formula: cut-edge-crossing] -/
theorem interp_eq_zero_iff {φa φb : ℝ} (h : φa ≠ φb) (t : ℝ) :
    (1 - t) * φa + t * φb = 0 ↔ t = Ddg.crossingParam φa φb := by
  have h' : φa - φb ≠ 0 := sub_ne_zero.mpr h
  rw [crossingParam_eq, eq_div_iff h']
  constructor <;> intro h1 <;> linarith

/-- **Range of the parameter, `a` in region 0**: if `φ_a ≤ 0 < φ_b` then `0 ≤ t < 1`.

[formula: cut-edge-crossing] -/
theorem crossingParam_mem_Ico {φa φb : ℝ} (ha : φa ≤ 0) (hb : 0 < φb) :
    0 ≤ Ddg.crossingParam φa φb ∧ Ddg.crossingParam φa φb < 1 := by
  have h : φa - φb < 0 := by linarith
  rw [crossingParam_eq]
  refine ⟨div_nonneg_of_nonpos ha h.le, ?_⟩
  rw [div_lt_one_of_neg h]
  linarith

/-- **Range of the parameter, `b` in region 0**: if `φ_b ≤ 0 < φ_a` then `0 < t ≤ 1`.

[formula: cut-edge-crossing] -/
theorem crossingParam_mem_Ioc {φa φb : ℝ} (ha : 0 < φa) (hb : φb ≤ 0) :
    0 < Ddg.crossingParam φa φb ∧ Ddg.crossingParam φa φb ≤ 1 := by
  have h : 0 < φa - φb := by linarith
  rw [crossingParam_eq]
  refine ⟨div_pos ha h, ?_⟩
  rw [div_le_one h]
  linarith

/-- On a crossed edge (`region(a) ≠ region(b)`, the region being `0` where `φ ≤ 0`) the two end
values differ and the parameter lies in `[0, 1]`; it is `0` exactly when `φ_a = 0` and `1`
exactly when `φ_b = 0`. -/
theorem crossingParam_mem_Icc {φa φb : ℝ} (h : (φa ≤ 0 ∧ 0 < φb) ∨ (φb ≤ 0 ∧ 0 < φa)) :
    φa ≠ φb ∧ 0 ≤ Ddg.crossingParam φa φb ∧ Ddg.crossingParam φa φb ≤ 1 ∧
      (Ddg.crossingParam φa φb = 0 ↔ φa = 0) ∧ (Ddg.crossingParam φa φb = 1 ↔ φb = 0) := by
  have hne : φa ≠ φb := by
    rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;> intro h3 <;> linarith
  have h' : φa - φb ≠ 0 := sub_ne_zero.mpr hne
  refine ⟨hne, ?_, ?_, ?_, ?_⟩
  · rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact (crossingParam_mem_Ico h1 h2).1
    · exact (crossingParam_mem_Ioc h2 h1).1.le
  · rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact (crossingParam_mem_Ico h1 h2).2.le
    · exact (crossingParam_mem_Ioc h2 h1).2
  · rw [crossingParam_eq, div_eq_zero_iff]
    exact ⟨fun h1 => h1.resolve_right h', fun h1 => Or.inl h1⟩
  · rw [crossingParam_eq, div_eq_one_iff_eq h']
    constructor <;> intro h1 <;> linarith

/-- The oracle's clamp at `ℝ`. -/
theorem clampCrossing_eq (δ t : ℝ) : Ddg.clampCrossing δ t = max δ (min (1 - δ) t) := by
  simp only [Ddg.clampCrossing, Oracle.one_eq]

/-- **The clamped parameter lies in `[δ, 1 − δ]`** (for `δ ≤ ½`, so that the interval is not
empty).

[formula: cut-edge-crossing] -/
theorem clampCrossing_mem {δ : ℝ} (hδ : δ ≤ 1 / 2) (t : ℝ) :
    δ ≤ Ddg.clampCrossing δ t ∧ Ddg.clampCrossing δ t ≤ 1 - δ := by
  rw [clampCrossing_eq]
  exact ⟨le_max_left _ _, max_le (by linarith) (min_le_left _ _)⟩

/-- For `0 < δ ≤ ½` the clamped parameter is strictly inside the edge. -/
theorem clampCrossing_mem_Ioo {δ : ℝ} (hδ0 : 0 < δ) (hδ : δ ≤ 1 / 2) (t : ℝ) :
    0 < Ddg.clampCrossing δ t ∧ Ddg.clampCrossing δ t < 1 := by
  obtain ⟨h1, h2⟩ := clampCrossing_mem hδ t
  constructor <;> linarith

/-- **The clamp is the identity exactly on `[δ, 1 − δ]`.**

[formula: cut-edge-crossing] -/
theorem clampCrossing_eq_self_iff {δ : ℝ} (hδ : δ ≤ 1 / 2) (t : ℝ) :
    Ddg.clampCrossing δ t = t ↔ δ ≤ t ∧ t ≤ 1 - δ := by
  constructor
  · intro h
    have := clampCrossing_mem hδ t
    rw [h] at this
    exact this
  · rintro ⟨h1, h2⟩
    rw [clampCrossing_eq, min_eq_right h2, max_eq_right h1]

/-- **`crossingMoved` reports exactly the crossings that the clamp moved.**

[formula: cut-edge-crossing] -/
theorem crossingMoved_eq_true_iff {δ : ℝ} (hδ : δ ≤ 1 / 2) (t : ℝ) :
    Ddg.crossingMoved δ t = true ↔ Ddg.clampCrossing δ t ≠ t := by
  rw [Ne, clampCrossing_eq_self_iff hδ]
  simp only [Ddg.crossingMoved, Bool.or_eq_true, decide_eq_true_eq, Oracle.one_eq]
  constructor
  · rintro (h | h) ⟨h1, h2⟩ <;> linarith
  · intro h
    by_contra hcon
    push Not at hcon
    exact h ⟨hcon.1, hcon.2⟩

/-- **The clamp moves a parameter of `[0, 1]` by at most `δ`.**

[formula: cut-edge-crossing] -/
theorem abs_clampCrossing_sub_le {δ t : ℝ} (hδ0 : 0 ≤ δ) (hδ : δ ≤ 1 / 2) (ht0 : 0 ≤ t)
    (ht1 : t ≤ 1) : |Ddg.clampCrossing δ t - t| ≤ δ := by
  rw [clampCrossing_eq, abs_le]
  rcases le_total (1 - δ) t with h | h
  · rw [min_eq_left h, max_eq_right (by linarith)]
    constructor <;> linarith
  · rw [min_eq_right h]
    rcases le_total δ t with h' | h'
    · rw [max_eq_right h']
      constructor <;> linarith
    · rw [max_eq_left h']
      constructor <;> linarith

/-! ## Lengths inside a face from barycentric coordinates -/

section Barycentric

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The oracle's `baryLength2` at `ℝ`. -/
theorem baryLength2_eq (lab lbc lca da db dc : ℝ) :
    Ddg.baryLength2 lab lbc lca da db dc =
      -(da * db * lab ^ 2 + db * dc * lbc ^ 2 + dc * da * lca ^ 2) := by
  unfold Ddg.baryLength2
  ring

/-- The oracle's `baryLength` at `ℝ`. -/
theorem baryLength_eq (lab lbc lca da db dc : ℝ) :
    Ddg.baryLength lab lbc lca da db dc = √(max 0 (Ddg.baryLength2 lab lbc lca da db dc)) := by
  simp only [Ddg.baryLength, Oracle.sqrt_eq, Oracle.zero_eq]

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

/-- **Barycentric length formula** (SPEC §3.3). For three points `P_a, P_b, P_c` of a real inner
product space and a barycentric difference `(d_a, d_b, d_c)` with `d_a + d_b + d_c = 0`,
`‖d_a P_a + d_b P_b + d_c P_c‖² = −(d_a d_b ℓ_ab² + d_b d_c ℓ_bc² + d_c d_a ℓ_ca²)`,
the right-hand side being the oracle's `baryLength2` of the three side lengths.

[formula: cut-barycentric-length] -/
theorem baryLength2_eq_norm_sq (Pa Pb Pc : E) {da db dc : ℝ} (h : da + db + dc = 0) :
    Ddg.baryLength2 ‖Pa - Pb‖ ‖Pb - Pc‖ ‖Pc - Pa‖ da db dc =
      ‖da • Pa + db • Pb + dc • Pc‖ ^ 2 := by
  have hdc : dc = -da - db := by linarith
  subst hdc
  have e : da • Pa + db • Pb + (-da - db) • Pc = da • (Pa - Pc) + db • (Pb - Pc) := by
    simp only [smul_sub, sub_smul, neg_smul]
    abel
  have hab : Pa - Pb = (Pa - Pc) - (Pb - Pc) := by abel
  have hca : ‖Pc - Pa‖ = ‖Pa - Pc‖ := norm_sub_rev _ _
  rw [e, baryLength2_eq, hab, hca]
  generalize Pa - Pc = x
  generalize Pb - Pc = y
  rw [norm_sub_sq_real, norm_add_sq_real, norm_smul, norm_smul, real_inner_smul_left,
    real_inner_smul_right]
  simp only [Real.norm_eq_abs, mul_pow, sq_abs]
  ring

/-- The same without the square: the oracle's `baryLength` is the norm. -/
theorem baryLength_eq_norm (Pa Pb Pc : E) {da db dc : ℝ} (h : da + db + dc = 0) :
    Ddg.baryLength ‖Pa - Pb‖ ‖Pb - Pc‖ ‖Pc - Pa‖ da db dc = ‖da • Pa + db • Pb + dc • Pc‖ := by
  rw [baryLength_eq, baryLength2_eq_norm_sq Pa Pb Pc h, max_eq_right (sq_nonneg _),
    Real.sqrt_sq (norm_nonneg _)]

/-- **Distance of two points given by barycentric coordinates.** For `x = Σ x_i P_i` and
`y = Σ y_i P_i` with `Σ x_i = Σ y_i = 1`, the oracle's `baryLength` of the three side lengths and
the coordinate difference `x − y` is the distance of the two points.

[formula: cut-barycentric-length] -/
theorem baryLength_eq_dist (Pa Pb Pc : E) {xa xb xc ya yb yc : ℝ} (hx : xa + xb + xc = 1)
    (hy : ya + yb + yc = 1) :
    Ddg.baryLength (dist Pa Pb) (dist Pb Pc) (dist Pc Pa) (xa - ya) (xb - yb) (xc - yc) =
      dist (xa • Pa + xb • Pb + xc • Pc) (ya • Pa + yb • Pb + yc • Pc) := by
  simp only [dist_eq_norm]
  rw [baryLength_eq_norm Pa Pb Pc (by linarith : (xa - ya) + (xb - yb) + (xc - yc) = 0)]
  congr 1
  simp only [sub_smul]
  abel

end Barycentric

/-- **The intrinsic formula is nonnegative.** If the three side lengths satisfy the strict
triangle inequalities and `d_a + d_b + d_c = 0`, then `baryLength2 ≥ 0`: the floor at `0` in the
oracle's `baryLength` only guards against rounding. -/
theorem baryLength2_nonneg {lab lbc lca : ℝ} (h : IsTriangle lab lbc lca) {da db dc : ℝ}
    (hd : da + db + dc = 0) : 0 ≤ Ddg.baryLength2 lab lbc lca da db dc := by
  have hdc : dc = -da - db := by linarith
  subst hdc
  rw [baryLength2_eq]
  have hH := h.heron_pos
  have hc := h.pos_c
  have key : 4 * lca ^ 2 *
      -(da * db * lab ^ 2 + db * (-da - db) * lbc ^ 2 + (-da - db) * da * lca ^ 2) =
      (2 * lca ^ 2 * da + (lbc ^ 2 + lca ^ 2 - lab ^ 2) * db) ^ 2 +
        (lab + lbc + lca) * (lbc - lab + lca) * (lab - lbc + lca) * (lab - lca + lbc) *
          db ^ 2 := by
    ring
  have h1 : 0 ≤ 4 * lca ^ 2 *
      -(da * db * lab ^ 2 + db * (-da - db) * lbc ^ 2 + (-da - db) * da * lca ^ 2) := by
    rw [key]
    positivity
  have h2 : 0 < 4 * lca ^ 2 := by positivity
  exact nonneg_of_mul_nonneg_right h1 h2

/-! ## Areas of sub-triangles -/

/-- Sixteen times the squared area of a triangle as a function of its squared side lengths:
`2xy + 2yz + 2zx − x² − y² − z²` (Heron). -/
def heronPoly (x y z : ℝ) : ℝ := 2 * x * y + 2 * y * z + 2 * z * x - x ^ 2 - y ^ 2 - z ^ 2

/-- Heron's formula in terms of the squared sides. -/
theorem triArea_eq_sqrt_heronPoly (a b c : ℝ) :
    triArea a b c = √(heronPoly (a ^ 2) (b ^ 2) (c ^ 2)) / 4 := by
  unfold triArea heronPoly
  congr 2
  ring

/-- The oracle's `det3` at `ℝ` is the determinant of the `3 × 3` matrix of the rows. -/
theorem det3_eq (x0 x1 x2 y0 y1 y2 z0 z1 z2 : ℝ) :
    Ddg.det3 x0 x1 x2 y0 y1 y2 z0 z1 z2 =
      Matrix.det !![x0, x1, x2; y0, y1, y2; z0, z1, z2] := by
  rw [Matrix.det_fin_three]
  simp only [Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val, Ddg.det3]
  ring

/-- The oracle's `subTriangleArea` at `ℝ`. -/
theorem subTriangleArea_eq (det A : ℝ) : Ddg.subTriangleArea det A = |det| * A := by
  simp only [Ddg.subTriangleArea, Oracle.abs_eq]

/-- **Heron's polynomial of a sub-triangle.** For three barycentric rows (each with sum `1`),
sixteen times the squared area computed from the three squared side lengths `baryLength2` is
`det(rows)²` times sixteen times the squared area of the parent. A polynomial identity: no
hypothesis on the side lengths. -/
theorem heronPoly_baryLength2 (lab lbc lca : ℝ) {x0 x1 x2 y0 y1 y2 z0 z1 z2 : ℝ}
    (hx : x0 + x1 + x2 = 1) (hy : y0 + y1 + y2 = 1) (hz : z0 + z1 + z2 = 1) :
    heronPoly (Ddg.baryLength2 lab lbc lca (x0 - y0) (x1 - y1) (x2 - y2))
        (Ddg.baryLength2 lab lbc lca (y0 - z0) (y1 - z1) (y2 - z2))
        (Ddg.baryLength2 lab lbc lca (z0 - x0) (z1 - x1) (z2 - x2)) =
      Ddg.det3 x0 x1 x2 y0 y1 y2 z0 z1 z2 ^ 2 * heronPoly (lab ^ 2) (lbc ^ 2) (lca ^ 2) := by
  have hx2 : x2 = 1 - x0 - x1 := by linarith
  have hy2 : y2 = 1 - y0 - y1 := by linarith
  have hz2 : z2 = 1 - z0 - z1 := by linarith
  subst hx2 hy2 hz2
  simp only [baryLength2_eq]
  unfold Ddg.det3 heronPoly
  ring

/-- **Area of a sub-triangle** (SPEC §3.3). Let the parent face have side lengths
`ℓ_ab, ℓ_bc, ℓ_ca` (strict triangle inequalities) and area `A_f`, and let three points of the
face be given by barycentric rows `x, y, z` (each with sum `1`). The Heron area of the triangle
whose side lengths are the oracle's `baryLength` of the three coordinate differences equals the
oracle's `subTriangleArea (det3 rows) A_f = |det(rows)| · A_f`.

[formula: cut-subtriangle-area] -/
theorem subTriangleArea_eq_triArea {lab lbc lca : ℝ} (h : IsTriangle lab lbc lca)
    {x0 x1 x2 y0 y1 y2 z0 z1 z2 : ℝ}
    (hx : x0 + x1 + x2 = 1) (hy : y0 + y1 + y2 = 1) (hz : z0 + z1 + z2 = 1) :
    Ddg.subTriangleArea (Ddg.det3 x0 x1 x2 y0 y1 y2 z0 z1 z2) (triArea lab lbc lca) =
      triArea (Ddg.baryLength lab lbc lca (x0 - y0) (x1 - y1) (x2 - y2))
        (Ddg.baryLength lab lbc lca (y0 - z0) (y1 - z1) (y2 - z2))
        (Ddg.baryLength lab lbc lca (z0 - x0) (z1 - x1) (z2 - x2)) := by
  have hq1 := baryLength2_nonneg h (by linarith : (x0 - y0) + (x1 - y1) + (x2 - y2) = 0)
  have hq2 := baryLength2_nonneg h (by linarith : (y0 - z0) + (y1 - z1) + (y2 - z2) = 0)
  have hq3 := baryLength2_nonneg h (by linarith : (z0 - x0) + (z1 - x1) + (z2 - x2) = 0)
  have hp : ∀ da db dc : ℝ, 0 ≤ Ddg.baryLength2 lab lbc lca da db dc →
      Ddg.baryLength lab lbc lca da db dc ^ 2 = Ddg.baryLength2 lab lbc lca da db dc := by
    intro da db dc hq
    rw [baryLength_eq, max_eq_right hq, Real.sq_sqrt hq]
  have hT := triArea_pos h
  have h16 : heronPoly (lab ^ 2) (lbc ^ 2) (lca ^ 2) = (4 * triArea lab lbc lca) ^ 2 := by
    have := sixteen_mul_triArea_sq h
    unfold heronPoly
    linear_combination (-1 : ℝ) * this
  rw [triArea_eq_sqrt_heronPoly (Ddg.baryLength lab lbc lca (x0 - y0) (x1 - y1) (x2 - y2)),
    hp _ _ _ hq1, hp _ _ _ hq2, hp _ _ _ hq3, heronPoly_baryLength2 lab lbc lca hx hy hz, h16,
    ← mul_pow, Real.sqrt_sq_eq_abs, abs_mul, abs_of_pos (by positivity : 0 < 4 * triArea lab lbc lca),
    subTriangleArea_eq]
  ring

/-- **Area of a triangle spanned by three points given in barycentric coordinates.** For a
nondegenerate triangle `P_a P_b P_c` in a real inner product space and three points
`X = Σ x_i P_i`, `Y = Σ y_i P_i`, `Z = Σ z_i P_i` (coefficient sums `1`), the Heron area of
`XYZ` is `|det(x; y; z)|` times the Heron area of `P_a P_b P_c`.

[formula: cut-subtriangle-area] -/
theorem triArea_bary_points {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (Pa Pb Pc : E) (h : IsTriangle (dist Pa Pb) (dist Pb Pc) (dist Pc Pa))
    {x0 x1 x2 y0 y1 y2 z0 z1 z2 : ℝ}
    (hx : x0 + x1 + x2 = 1) (hy : y0 + y1 + y2 = 1) (hz : z0 + z1 + z2 = 1) :
    triArea (dist (x0 • Pa + x1 • Pb + x2 • Pc) (y0 • Pa + y1 • Pb + y2 • Pc))
        (dist (y0 • Pa + y1 • Pb + y2 • Pc) (z0 • Pa + z1 • Pb + z2 • Pc))
        (dist (z0 • Pa + z1 • Pb + z2 • Pc) (x0 • Pa + x1 • Pb + x2 • Pc)) =
      |Ddg.det3 x0 x1 x2 y0 y1 y2 z0 z1 z2| *
        triArea (dist Pa Pb) (dist Pb Pc) (dist Pc Pa) := by
  rw [← baryLength_eq_dist Pa Pb Pc hx hy, ← baryLength_eq_dist Pa Pb Pc hy hz,
    ← baryLength_eq_dist Pa Pb Pc hz hx, ← subTriangleArea_eq_triArea h hx hy hz,
    subTriangleArea_eq]

/-- **The sub-triangles of a cut face fill it**, quad split by the diagonal `p–c`. Face `abc`
with cut vertices `p = (1 − t)a + t b`, `q = (1 − s)a + s c`: the determinants of the barycentric
rows of `(a, p, q)`, `(p, b, c)`, `(p, c, q)` are `t s`, `1 − t`, `t (1 − s)`; they are
nonnegative for `t, s ∈ [0, 1]` (orientation is inherited) and sum to `1`, so the areas
`|det| · A_f` sum to `A_f`.

[formula: cut-subtriangle-area] -/
theorem cutFace_det_sum_pc (t s : ℝ) :
    Ddg.det3 1 0 0 (1 - t) t 0 (1 - s) 0 s = t * s ∧
    Ddg.det3 (1 - t) t 0 0 1 0 0 0 1 = 1 - t ∧
    Ddg.det3 (1 - t) t 0 0 0 1 (1 - s) 0 s = t * (1 - s) ∧
    t * s + (1 - t) + t * (1 - s) = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> (try unfold Ddg.det3) <;> ring

/-- The same for the diagonal `b–q`: `(a, p, q)`, `(p, b, q)`, `(b, c, q)` have determinants
`t s`, `(1 − t) s`, `1 − s`, which sum to `1`.

[formula: cut-subtriangle-area] -/
theorem cutFace_det_sum_bq (t s : ℝ) :
    Ddg.det3 1 0 0 (1 - t) t 0 (1 - s) 0 s = t * s ∧
    Ddg.det3 (1 - t) t 0 0 1 0 (1 - s) 0 s = (1 - t) * s ∧
    Ddg.det3 0 1 0 0 0 1 (1 - s) 0 s = 1 - s ∧
    t * s + (1 - t) * s + (1 - s) = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> (try unfold Ddg.det3) <;> ring

end Cuts

end DdgProofs
