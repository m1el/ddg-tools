import DdgProofs.Developability.Vec
import DdgProofs.Intrinsic.Diamond

/-!
# Developability: classification of bad triangles for the small-angle remeshing (SPEC §5.3)

`crates/ddg-developability/SPEC.md` §5.3, decision D21 (paper Sec. 4.3, Fig. 18). A triangle is
*bad* if a corner angle is below `θ_s`, `0 < θ_s < π/3`. Two small angles make a *cap*, whose third
angle exceeds `π − 2θ_s` and which is flipped across its longest edge; one small angle makes a
*needle*, collapsed along its shortest edge.

* `smallAngleCount_le_two`, `third_angle_of_two_small`: for any three angles with sum `π` and
  `θ_s ≤ π/3`;
* `cornerAngleAt_eq_cornerAngle`: the corner angle `atan2(|a × b|, a·b)` of the code is the angle
  of the law of cosines of the three edge lengths (`DdgProofs.cornerAngle`), so the three corner
  angles of a nondegenerate face add up to `π` (`cornerAngleAt_add`);
* `cornerAngle_lt_iff`, `cornerAngleAt_lt_iff`: the larger of two corner angles is opposite the
  longer edge;
* `smallAngle_classification`: the statement used by `remesh::classify`.
-/

noncomputable section

namespace DdgProofs

namespace Dev

open scoped DdgProofs.Oracle
open Real Ddg.Dev

/-! ## Angles only -/

/-- With `θ_s ≤ π/3`, three angles with sum `π` are not all below `θ_s`. -/
theorem smallAngleCount_le_two {θs a0 a1 a2 : ℝ} (hθ : θs ≤ π / 3) (hsum : a0 + a1 + a2 = π) :
    smallAngleCount θs a0 a1 a2 ≤ 2 := by
  unfold smallAngleCount
  split_ifs with h0 h1 h2 <;> first | omega | (exfalso; linarith)

/-- Two angles below `θ_s ≤ π/3` (sum `π`): the third exceeds `π − 2θ_s ≥ θ_s` and is the
largest. -/
theorem third_angle_of_two_small {θs a0 a1 a2 : ℝ} (hθ : θs ≤ π / 3) (hsum : a0 + a1 + a2 = π)
    (h0 : a0 < θs) (h1 : a1 < θs) : π - 2 * θs < a2 ∧ a0 < a2 ∧ a1 < a2 := by
  refine ⟨by linarith, by linarith, by linarith⟩

/-! ## Larger angle opposite the longer side -/

theorem cornerAngle_comm (a b c : ℝ) : cornerAngle a b c = cornerAngle a c b := by
  unfold cornerAngle; congr 1; ring_nf

/-- In a triangle with sides `a b c`, the angle opposite `a` is smaller than the angle opposite
`c` iff `a < c`. -/
theorem cornerAngle_lt_iff {a b c : ℝ} (h : IsTriangle a b c) :
    cornerAngle a b c < cornerAngle c a b ↔ a < c := by
  have ha := h.pos_a
  have hb := h.pos_b
  have hc := h.pos_c
  have b1 := cos_arg_bounds h
  have b2 := cos_arg_bounds h.rotate.rotate
  unfold cornerAngle
  rw [strictAntiOn_arccos.lt_iff_gt ⟨b1.1.le, b1.2.le⟩ ⟨b2.1.le, b2.2.le⟩]
  rw [div_lt_div_iff₀ (by positivity) (by positivity)]
  have hfac : 0 < (a + c) ^ 2 - b ^ 2 := by
    have := h.lt_a; have := h.lt_b; have := h.lt_c; nlinarith
  have key : (b ^ 2 + c ^ 2 - a ^ 2) * (2 * a * b) - (a ^ 2 + b ^ 2 - c ^ 2) * (2 * b * c) =
      2 * b * ((c - a) * ((a + c) ^ 2 - b ^ 2)) := by ring
  constructor
  · intro hlt
    by_contra hge
    replace hge := not_lt.mp hge
    have : (c - a) * ((a + c) ^ 2 - b ^ 2) ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg (by linarith) hfac.le
    nlinarith
  · intro hlt
    have : 0 < (c - a) * ((a + c) ^ 2 - b ^ 2) := mul_pos (by linarith) hfac
    nlinarith

/-- The same with the angle opposite `a` against the one opposite `b`. -/
theorem cornerAngle_lt_iff' {a b c : ℝ} (h : IsTriangle a b c) :
    cornerAngle a b c < cornerAngle b c a ↔ a < b := by
  have h' : IsTriangle a c b := ⟨by linarith [h.lt_a], by linarith [h.lt_c], by linarith [h.lt_b]⟩
  rw [cornerAngle_comm a b c, cornerAngle_comm b c a]
  exact cornerAngle_lt_iff h'

/-! ## Corner angles of positions -/

theorem norm_sub_comm (p q : V3) : norm (sub p q) = norm (sub q p) := by
  rw [norm_def, norm_def]; congr 1; simp only [norm2_def, sub_x, sub_y, sub_z]; ring

theorem faceNormal_rotate (p q r : V3) : faceNormal q r p = faceNormal p q r := by
  ext <;> simp [faceNormal] <;> ring

/-- The corner angle of positions is the law-of-cosines angle of the edge lengths:
`cornerAngleAt p q r = cornerAngle |r − q| |p − r| |q − p|` (the lengths of the edges opposite
`p, q, r`, in this cyclic order), for a nondegenerate face. -/
theorem cornerAngleAt_eq_cornerAngle (p q r : V3) (hn : faceNormal p q r ≠ zero) :
    cornerAngleAt p q r = cornerAngle (norm (sub r q)) (norm (sub p r)) (norm (sub q p)) := by
  rw [norm_sub_comm p r]
  set a := sub q p with ha_def
  set b := sub r p with hb_def
  have hn' : cross a b ≠ zero := hn
  have ha : a ≠ zero := by
    intro h; apply hn'; rw [h]; ext <;> simp
  have hb : b ≠ zero := by
    intro h; apply hn'; rw [h]; ext <;> simp
  have hapos := (norm_pos_iff' a).mpr ha
  have hbpos := (norm_pos_iff' b).mpr hb
  have hL := norm2_cross_add_sq_dot a b
  have hy := norm_sq (cross a b)
  have hxy : dot a b ^ 2 + norm (cross a b) ^ 2 = (norm a * norm b) ^ 2 := by
    rw [hy, mul_pow, norm_sq, norm_sq]; linarith
  have hpos : 0 < dot a b ^ 2 + norm (cross a b) ^ 2 := by
    rw [hxy]; positivity
  have hcos := cos_atan2 (norm_nonneg' (cross a b)) hpos
  rw [hxy, sqrt_sq (by positivity)] at hcos
  have hrange := atan2_mem_Icc (norm_nonneg' (cross a b)) (dot a b)
  have hc2 : norm (sub r q) ^ 2 = norm b ^ 2 + norm a ^ 2 - 2 * dot a b := by
    rw [norm_sq, norm_sq, norm_sq]
    simp only [ha_def, hb_def, norm2_def, dot_def, sub_x, sub_y, sub_z]; ring
  have hA : cornerAngleAt p q r = Ddg.Dev.atan2 (norm (cross a b)) (dot a b) := by
    simp only [cornerAngleAt, ha_def, hb_def]
  rw [hA, ← arccos_cos hrange.1 hrange.2, hcos]
  unfold cornerAngle
  rw [hc2]
  congr 1
  field_simp
  ring

/-- The edge lengths `|r − q|, |p − r|, |q − p|` of a nondegenerate face satisfy the strict
triangle inequalities. -/
theorem isTriangle_of_faceNormal (p q r : V3) (hn : faceNormal p q r ≠ zero) :
    IsTriangle (norm (sub r q)) (norm (sub p r)) (norm (sub q p)) := by
  apply Intrinsic.isTriangle_of_heron_pos (norm_nonneg' _) (norm_nonneg' _) (norm_nonneg' _)
  set x := norm (sub r q)
  set y := norm (sub p r)
  set z := norm (sub q p)
  have hP : (x + y + z) * (y - x + z) * (x - y + z) * (x - z + y) =
      2 * x ^ 2 * y ^ 2 + 2 * y ^ 2 * z ^ 2 + 2 * z ^ 2 * x ^ 2 - (x ^ 2) ^ 2 - (y ^ 2) ^ 2 -
        (z ^ 2) ^ 2 := by ring
  rw [hP, norm_sq, norm_sq, norm_sq]
  have h4 : 0 < norm2 (faceNormal p q r) := (norm2_pos_iff _).mpr hn
  have : 2 * norm2 (sub r q) * norm2 (sub p r) + 2 * norm2 (sub p r) * norm2 (sub q p) +
      2 * norm2 (sub q p) * norm2 (sub r q) - norm2 (sub r q) ^ 2 - norm2 (sub p r) ^ 2 -
      norm2 (sub q p) ^ 2 = 4 * norm2 (faceNormal p q r) := by
    simp only [faceNormal, norm2_def, sub_x, sub_y, sub_z, cross_x, cross_y, cross_z]; ring
  rw [this]; linarith

/-- The three corner angles of a nondegenerate face add up to `π`. -/
theorem cornerAngleAt_add (fi fj fk : V3) (hn : faceNormal fi fj fk ≠ zero) :
    cornerAngleAt fi fj fk + cornerAngleAt fj fk fi + cornerAngleAt fk fi fj = π := by
  have hn2 : faceNormal fj fk fi ≠ zero := by rw [faceNormal_rotate]; exact hn
  have hn3 : faceNormal fk fi fj ≠ zero := by rw [faceNormal_rotate]; exact hn2
  rw [cornerAngleAt_eq_cornerAngle _ _ _ hn, cornerAngleAt_eq_cornerAngle _ _ _ hn2,
    cornerAngleAt_eq_cornerAngle _ _ _ hn3]
  exact cornerAngle_add_cornerAngle_add_cornerAngle (isTriangle_of_faceNormal fi fj fk hn)

/-- **Larger corner angle opposite the longer edge**: in a nondegenerate face `(f_i, f_j, f_k)`
with `ℓ_i = |f_k − f_j|`, `ℓ_j = |f_i − f_k|`, `ℓ_k = |f_j − f_i|` (each opposite its corner):
`θ_i < θ_k ↔ ℓ_i < ℓ_k` and `θ_i < θ_j ↔ ℓ_i < ℓ_j`. -/
theorem cornerAngleAt_lt_iff (fi fj fk : V3) (hn : faceNormal fi fj fk ≠ zero) :
    (cornerAngleAt fi fj fk < cornerAngleAt fk fi fj ↔ norm (sub fk fj) < norm (sub fj fi)) ∧
      (cornerAngleAt fi fj fk < cornerAngleAt fj fk fi ↔ norm (sub fk fj) < norm (sub fi fk)) := by
  have hn2 : faceNormal fj fk fi ≠ zero := by rw [faceNormal_rotate]; exact hn
  have hn3 : faceNormal fk fi fj ≠ zero := by rw [faceNormal_rotate]; exact hn2
  have hT := isTriangle_of_faceNormal fi fj fk hn
  rw [cornerAngleAt_eq_cornerAngle _ _ _ hn, cornerAngleAt_eq_cornerAngle _ _ _ hn2,
    cornerAngleAt_eq_cornerAngle _ _ _ hn3]
  exact ⟨cornerAngle_lt_iff hT, cornerAngle_lt_iff' hT⟩

/-- **Classification of a bad triangle** (SPEC §5.3): for a nondegenerate face with corner angles
`θ_i, θ_j, θ_k` (`cornerAngleAt`) and `θ_s ≤ π/3`, at most two angles are below `θ_s`; if two are,
`θ_i, θ_j < θ_s` say, then `θ_k > π − 2θ_s` and the edge `f_i f_j` opposite `θ_k` is strictly the
longest edge (so the flip of a cap is the flip of its longest edge). The smaller of two angles is
opposite the shorter edge (`cornerAngleAt_lt_iff`), so a needle's smallest angle is opposite its
shortest edge.

[formula: dev-small-angle-classification] -/
theorem smallAngle_classification (θs : ℝ) (fi fj fk : V3) (hn : faceNormal fi fj fk ≠ zero)
    (hθ : θs ≤ π / 3) :
    smallAngleCount θs (cornerAngleAt fi fj fk) (cornerAngleAt fj fk fi)
        (cornerAngleAt fk fi fj) ≤ 2 ∧
      (cornerAngleAt fi fj fk < θs → cornerAngleAt fj fk fi < θs →
        π - 2 * θs < cornerAngleAt fk fi fj ∧ norm (sub fk fj) < norm (sub fj fi) ∧
          norm (sub fi fk) < norm (sub fj fi)) := by
  have hsum := cornerAngleAt_add fi fj fk hn
  refine ⟨smallAngleCount_le_two hθ hsum, fun h0 h1 => ?_⟩
  obtain ⟨h2, h02, h12⟩ := third_angle_of_two_small hθ hsum h0 h1
  have hn2 : faceNormal fj fk fi ≠ zero := by rw [faceNormal_rotate]; exact hn
  refine ⟨h2, (cornerAngleAt_lt_iff fi fj fk hn).1.mp h02, ?_⟩
  -- in the rotated face `(f_j, f_k, f_i)`: `θ_j < θ_k ↔ ℓ_j < ℓ_k`
  have := (cornerAngleAt_lt_iff fj fk fi hn2).2
  exact this.mp h12
