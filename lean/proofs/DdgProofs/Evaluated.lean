import DdgProofs.GaussBonnet

/-!
# The formulas in the form the Rust code evaluates them

`crates/ddg-core` evaluates the triangle formulas in numerically stable arrangements (Kahan's
parenthesization with sorted sides, `atan2` of two square roots) and stores the curvatures as
per-vertex fields that vanish on the other class of vertices. This file proves that, in exact
arithmetic, those arrangements are the quantities of `DdgProofs.Triangle` and
`DdgProofs.GaussBonnet`:

* `kahanArea_of_sorted_rearrangement` — `triangle_area`: sort the sides, test
  `z − (x − y) > 0`, evaluate `¼√((x+(y+z))(z−(x−y))(z+(x−y))(x+(y−z)))`;
* `kahanCornerAngle_eq_cornerAngle` — `corner_angle`: `2·atan2(√num, √den)`;
* `kahanCornerCot_eq_cornerCot` — `corner_cotan`: `(c² + (b − a)(b + a))/(4A)`;
* `TriMesh.gauss_bonnet_fields` — `angle_defects`, `boundary_geodesic_curvatures`.

Nothing is claimed about floating-point rounding.
-/

noncomputable section

namespace DdgProofs

open Real

/-! ## `atan2` -/

/-- The two-argument arctangent: `atan2 y x` is the argument of the complex number `x + i·y`
(the angle of the point `(x, y)`, in `(−π, π]`). -/
def atan2 (y x : ℝ) : ℝ := Complex.arg ⟨x, y⟩

/-- For `x > 0`, `atan2 y x = arctan (y / x)`. -/
theorem atan2_eq_arctan {x : ℝ} (hx : 0 < x) (y : ℝ) : atan2 y x = arctan (y / x) := by
  unfold atan2
  rw [Complex.arg_of_re_nonneg (x := ⟨x, y⟩) hx.le, arctan_eq_arcsin,
    Complex.norm_eq_sqrt_sq_add_sq]
  congr 1
  have h1 : 1 + (y / x) ^ 2 = (x ^ 2 + y ^ 2) / x ^ 2 := by field_simp
  rw [h1, sqrt_div (by positivity), sqrt_sq hx.le]
  have h2 : 0 < √(x ^ 2 + y ^ 2) := sqrt_pos.mpr (by positivity)
  field_simp

/-! ## Rearrangements of the three sides -/

/-- `(x, y, z)` is a rearrangement of `(a, b, c)`. -/
def IsRearrangement (x y z a b c : ℝ) : Prop :=
  (x = a ∧ y = b ∧ z = c) ∨ (x = a ∧ y = c ∧ z = b) ∨ (x = b ∧ y = a ∧ z = c) ∨
    (x = b ∧ y = c ∧ z = a) ∨ (x = c ∧ y = a ∧ z = b) ∨ (x = c ∧ y = b ∧ z = a)

variable {a b c x y z : ℝ}

/-- The area does not depend on the order of the sides. -/
theorem triArea_of_rearrangement (h : IsRearrangement x y z a b c) :
    triArea x y z = triArea a b c := by
  rcases h with ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl⟩ |
    ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl⟩ <;> first | rfl | (unfold triArea; congr 2; ring)

/-- Being the side lengths of a triangle does not depend on the order of the sides. -/
theorem isTriangle_of_rearrangement (h : IsRearrangement x y z a b c) :
    IsTriangle x y z ↔ IsTriangle a b c := by
  rcases h with ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl⟩ |
    ⟨rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl⟩ <;>
    exact ⟨fun ⟨h1, h2, h3⟩ => ⟨by linarith, by linarith, by linarith⟩,
      fun ⟨h1, h2, h3⟩ => ⟨by linarith, by linarith, by linarith⟩⟩

/-- **The validity test of `sorted_valid`.** For sides sorted `x ≥ y ≥ z`, the three strict
triangle inequalities hold iff `z − (x − y) > 0`. (In particular the test implies `z > 0`.) -/
theorem isTriangle_iff_slack_pos (hxy : y ≤ x) (hyz : z ≤ y) :
    IsTriangle x y z ↔ 0 < z - (x - y) := by
  constructor
  · rintro ⟨h1, _, _⟩; linarith
  · intro h
    exact ⟨by linarith, by linarith, by linarith⟩

/-! ## Area -/

/-- Heron's formula in Kahan's parenthesization, as `triangle_area` evaluates it for sides
sorted `x ≥ y ≥ z`: `¼ √((x + (y + z))·(z − (x − y))·(z + (x − y))·(x + (y − z)))`. -/
def kahanArea (x y z : ℝ) : ℝ :=
  1 / 4 * √((x + (y + z)) * (z - (x - y)) * (z + (x - y)) * (x + (y - z)))

/-- Kahan's arrangement is Heron's formula (for any three numbers, sorted or not). -/
theorem kahanArea_eq_triArea (x y z : ℝ) : kahanArea x y z = triArea x y z := by
  unfold kahanArea triArea
  have h : (x + (y + z)) * (z - (x - y)) * (z + (x - y)) * (x + (y - z)) =
      (x + y + z) * (y - x + z) * (x - y + z) * (x - z + y) := by ring
  rw [h]
  ring

/-- **`triangle_area` computes the area.** Let `(x, y, z)` be the sides `(a, b, c)` sorted in
descending order. Then the test `z − (x − y) > 0` holds exactly when `a, b, c` are the side
lengths of a nondegenerate triangle, and Kahan's expression in `x, y, z` equals the Heron area of
`a, b, c` (which is `½·b·c·sin A`, see `sin_cornerAngle`, `half_mul_sin_angle_eq_triArea`).

[formula: triangle-area] -/
theorem kahanArea_of_sorted_rearrangement (h : IsRearrangement x y z a b c) (hxy : y ≤ x)
    (hyz : z ≤ y) :
    (IsTriangle a b c ↔ 0 < z - (x - y)) ∧ kahanArea x y z = triArea a b c :=
  ⟨(isTriangle_of_rearrangement h).symm.trans (isTriangle_iff_slack_pos hxy hyz),
    (kahanArea_eq_triArea x y z).trans (triArea_of_rearrangement h)⟩

/-- `16 A² = (a + b + c)(−a + b + c)(a − b + c)(a + b − c)`.

[formula: triangle-area] -/
theorem sixteen_mul_triArea_sq_eq (h : IsTriangle a b c) :
    16 * triArea a b c ^ 2 = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) := by
  have := sq_four_mul_triArea h
  nlinarith

/-! ## Corner angle -/

/-- **Half-angle form with `atan2`**: for a nondegenerate triangle the angle opposite `a` is
`2·atan2(√((a−b+c)(a+b−c)), √((a+b+c)(−a+b+c)))`.

[formula: corner-angle] -/
theorem two_mul_atan2_sqrt_eq_cornerAngle (h : IsTriangle a b c) :
    2 * atan2 (√((a - b + c) * (a - c + b))) (√((a + b + c) * (b - a + c))) =
      cornerAngle a b c := by
  have f1 : 0 < a + b + c := by linarith [h.pos_a, h.pos_b, h.pos_c]
  have f2 : 0 < b - a + c := by linarith [h.lt_a]
  have f3 : 0 < a - b + c := by linarith [h.lt_b]
  have f4 : 0 < a - c + b := by linarith [h.lt_c]
  have hden : 0 < √((a + b + c) * (b - a + c)) := sqrt_pos.mpr (by positivity)
  rw [cornerAngle_eq_two_arctan h, atan2_eq_arctan hden,
    sqrt_div (by positivity : 0 ≤ (a - b + c) * (a - c + b))]

/-- The corner angle as `corner_angle` in `crates/ddg-core` evaluates it (Kahan's arrangement):
with `c = opposite` and the other two sides ordered `a ≥ b`,
`mu = if b ≥ c then c − (a − b) else b − (a − c)`, `num = ((a − b) + c)·mu`,
`den = (a + (b + c))·((a − c) + b)`, result `2·atan2(√num, √den)`. -/
def kahanCornerAngle (opposite side1 side2 : ℝ) : ℝ :=
  2 * atan2
    (√(((max side1 side2 - min side1 side2) + opposite) *
      (if opposite ≤ min side1 side2 then opposite - (max side1 side2 - min side1 side2)
        else min side1 side2 - (max side1 side2 - opposite))))
    (√((max side1 side2 + (min side1 side2 + opposite)) *
      ((max side1 side2 - opposite) + min side1 side2)))

/-- **`corner_angle` computes the law-of-cosines angle**: for a nondegenerate triangle, Kahan's
arrangement with `atan2` equals `arccos((b² + c² − a²)/(2bc))`, `a` the opposite side.

[formula: corner-angle] -/
theorem kahanCornerAngle_eq_cornerAngle (h : IsTriangle a b c) :
    kahanCornerAngle a b c = cornerAngle a b c := by
  rw [← two_mul_atan2_sqrt_eq_cornerAngle h]
  unfold kahanCornerAngle
  have hmu : (if a ≤ min b c then a - (max b c - min b c) else min b c - (max b c - a)) =
      a - max b c + min b c := by
    split_ifs <;> ring
  rw [hmu]
  rcases le_total b c with hbc | hbc
  · rw [max_eq_right hbc, min_eq_left hbc]
    congr 3 <;> ring
  · rw [max_eq_left hbc, min_eq_right hbc]
    congr 3 <;> ring

/-! ## Cotangent -/

/-- The cotangent as `corner_cotan` evaluates it: with the other two sides ordered `b ≥ c`,
`(c·c + (b − a)·(b + a)) / (4·A)`, `a = opposite`. -/
def kahanCornerCot (opposite side1 side2 : ℝ) : ℝ :=
  (min side1 side2 * min side1 side2 +
      (max side1 side2 - opposite) * (max side1 side2 + opposite)) /
    (4 * triArea opposite side1 side2)

/-- **`corner_cotan` computes `(b² + c² − a²)/(4A)`**, which is the cotangent of the angle
opposite `a` (`cornerCot_eq_cot`).

[formula: corner-cot] -/
theorem kahanCornerCot_eq_cornerCot (a b c : ℝ) : kahanCornerCot a b c = cornerCot a b c := by
  unfold kahanCornerCot cornerCot
  rcases le_total b c with hbc | hbc
  · rw [max_eq_right hbc, min_eq_left hbc]
    congr 1
    ring
  · rw [max_eq_left hbc, min_eq_right hbc]
    congr 1
    ring

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

end DdgProofs
