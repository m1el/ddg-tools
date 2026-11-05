import DdgProofs.Oracle.Scalar

/-!
# Developability: vectors of `ℝ³`, derivatives along curves, `atan2`

Infrastructure for `DdgProofs/Developability/*` ("Developability of Triangle Meshes", Stein,
Grinspun, Crane 2018; contract `crates/ddg-developability/SPEC.md`).

The oracle's formulas (`lean/oracle/DdgOracle/Formulas/Developability.lean`, namespace `Ddg.Dev`)
work on `Ddg.Point α`. Here they are taken at `ℝ` (the scoped instance of
`DdgProofs/Oracle/Scalar.lean`), and:

* every vector operation is unfolded to coordinates by `simp` lemmas, so that identities between
  vectors are proved by `ext` and `ring`;
* a **derivative along a curve** of a vector-valued function is a derivative of each coordinate
  (`HasDerivAt3`). All derivative statements of this development are of this kind: the curves are
  arbitrary differentiable curves `t ↦ f(t)` of positions (or of normals, angles), so the
  statements contain the directional derivatives (`f(t) = f + tδ`) and the chain rule along any
  differentiable motion;
* `Ddg.Dev.atan2`, for `y > 0`, is `π/2 − arctan(x/y)`; its cosine and sine are `x/r`, `y/r`; its
  derivative along a curve is `(x y' − y x')/(x² + y²)`.
-/

noncomputable section

namespace DDG

namespace Dev

open scoped DDG.Oracle
open Real Filter Topology

/-- Vectors of `ℝ³` as the oracle stores them. -/
abbrev V3 := Ddg.Point ℝ

@[ext] theorem V3.ext {p q : V3} (hx : p.x = q.x) (hy : p.y = q.y) (hz : p.z = q.z) : p = q := by
  cases p; cases q; simp_all

/-! ## Coordinates of the vector operations -/

section coords

variable (p q : V3) (s : ℝ)

@[simp] theorem add_x : (Ddg.Dev.add p q).x = p.x + q.x := rfl
@[simp] theorem add_y : (Ddg.Dev.add p q).y = p.y + q.y := rfl
@[simp] theorem add_z : (Ddg.Dev.add p q).z = p.z + q.z := rfl
@[simp] theorem sub_x : (Ddg.Dev.sub p q).x = p.x - q.x := rfl
@[simp] theorem sub_y : (Ddg.Dev.sub p q).y = p.y - q.y := rfl
@[simp] theorem sub_z : (Ddg.Dev.sub p q).z = p.z - q.z := rfl
@[simp] theorem neg_x : (Ddg.Dev.neg p).x = -p.x := rfl
@[simp] theorem neg_y : (Ddg.Dev.neg p).y = -p.y := rfl
@[simp] theorem neg_z : (Ddg.Dev.neg p).z = -p.z := rfl
@[simp] theorem smul_x : (Ddg.Dev.smul s p).x = s * p.x := rfl
@[simp] theorem smul_y : (Ddg.Dev.smul s p).y = s * p.y := rfl
@[simp] theorem smul_z : (Ddg.Dev.smul s p).z = s * p.z := rfl
@[simp] theorem zero_x : (Ddg.Dev.zero : V3).x = 0 := by simp [Ddg.Dev.zero]
@[simp] theorem zero_y : (Ddg.Dev.zero : V3).y = 0 := by simp [Ddg.Dev.zero]
@[simp] theorem zero_z : (Ddg.Dev.zero : V3).z = 0 := by simp [Ddg.Dev.zero]
@[simp] theorem cross_x : (Ddg.Dev.cross p q).x = p.y * q.z - p.z * q.y := rfl
@[simp] theorem cross_y : (Ddg.Dev.cross p q).y = p.z * q.x - p.x * q.z := rfl
@[simp] theorem cross_z : (Ddg.Dev.cross p q).z = p.x * q.y - p.y * q.x := rfl

theorem dot_def : Ddg.Dev.dot p q = p.x * q.x + p.y * q.y + p.z * q.z := rfl
theorem norm2_def : Ddg.Dev.norm2 p = p.x * p.x + p.y * p.y + p.z * p.z := rfl
theorem norm_def : Ddg.Dev.norm p = √(Ddg.Dev.norm2 p) := rfl
theorem det3_def (a b c : V3) : Ddg.Dev.det3 a b c = Ddg.Dev.dot (Ddg.Dev.cross a b) c := rfl
theorem normalize_def : Ddg.Dev.normalize p = Ddg.Dev.smul (1 / Ddg.Dev.norm p) p := by
  simp [Ddg.Dev.normalize]

@[simp] theorem scalar_ofNat (n : ℕ) : (Ddg.Scalar.ofNat n : ℝ) = n := rfl

end coords

/-- Unfolds `dot`, `norm2` to coordinates. -/
theorem dot_eq (p q : V3) : Ddg.Dev.dot p q = p.x * q.x + p.y * q.y + p.z * q.z := rfl

/-! ## Algebra -/

section algebra

variable (a b c d : V3)

open Ddg.Dev in
theorem norm2_nonneg : 0 ≤ norm2 a := by
  rw [norm2_def]; nlinarith [mul_self_nonneg a.x, mul_self_nonneg a.y, mul_self_nonneg a.z]

open Ddg.Dev in
theorem norm_nonneg' : 0 ≤ norm a := Real.sqrt_nonneg _

open Ddg.Dev in
theorem norm_sq : norm a ^ 2 = norm2 a := by
  rw [norm_def, sq_sqrt (norm2_nonneg a)]

open Ddg.Dev in
theorem norm_mul_self : norm a * norm a = norm2 a := by
  rw [← sq, norm_sq]

open Ddg.Dev in
theorem norm2_eq_zero_iff : norm2 a = 0 ↔ a = zero := by
  constructor
  · intro h
    rw [norm2_def] at h
    have hx : a.x = 0 := by nlinarith [mul_self_nonneg a.x, mul_self_nonneg a.y, mul_self_nonneg a.z]
    have hy : a.y = 0 := by nlinarith [mul_self_nonneg a.x, mul_self_nonneg a.y, mul_self_nonneg a.z]
    have hz : a.z = 0 := by nlinarith [mul_self_nonneg a.x, mul_self_nonneg a.y, mul_self_nonneg a.z]
    ext <;> simp [hx, hy, hz]
  · rintro rfl; simp [norm2_def]

open Ddg.Dev in
theorem norm2_pos_iff : 0 < norm2 a ↔ a ≠ zero := by
  rw [lt_iff_le_and_ne, ne_comm, Ne, norm2_eq_zero_iff]
  exact ⟨fun h => h.2, fun h => ⟨norm2_nonneg a, h⟩⟩

open Ddg.Dev in
theorem norm_pos_iff' : 0 < norm a ↔ a ≠ zero := by
  rw [norm_def, Real.sqrt_pos, norm2_pos_iff]

open Ddg.Dev in
/-- Lagrange's identity `|a × b|² + (a·b)² = |a|²|b|²`. -/
theorem norm2_cross_add_sq_dot : norm2 (cross a b) + dot a b ^ 2 = norm2 a * norm2 b := by
  simp only [norm2_def, dot_def, cross_x, cross_y, cross_z]; ring

open Ddg.Dev in
theorem dot_comm : dot a b = dot b a := by simp only [dot_def]; ring

open Ddg.Dev in
theorem cross_anticomm : cross a b = neg (cross b a) := by ext <;> simp <;> ring

open Ddg.Dev in
theorem dot_cross_self_left : dot (cross a b) a = 0 := by simp only [dot_def, cross_x, cross_y, cross_z]; ring

open Ddg.Dev in
theorem dot_cross_self_right : dot (cross a b) b = 0 := by simp only [dot_def, cross_x, cross_y, cross_z]; ring

open Ddg.Dev in
/-- `a × (b × c) = (a·c) b − (a·b) c`. -/
theorem cross_cross : cross a (cross b c) = sub (smul (dot a c) b) (smul (dot a b) c) := by
  ext <;> simp [dot_def] <;> ring

open Ddg.Dev in
/-- Scalar triple product, cyclic. -/
theorem dot_cross_cyclic : dot (cross a b) c = dot (cross b c) a := by
  simp only [dot_def, cross_x, cross_y, cross_z]; ring

open Ddg.Dev in
theorem norm_smul_of_nonneg {s : ℝ} (hs : 0 ≤ s) : norm (smul s a) = s * norm a := by
  rw [norm_def, norm_def]
  have : norm2 (smul s a) = s ^ 2 * norm2 a := by simp only [norm2_def, smul_x, smul_y, smul_z]; ring
  rw [this, sqrt_mul (sq_nonneg s), sqrt_sq hs]

open Ddg.Dev in
/-- `normalize a` is a unit vector for `a ≠ 0`. -/
theorem norm2_normalize (h : a ≠ zero) : norm2 (normalize a) = 1 := by
  have hp := (norm_pos_iff' a).mpr h
  have h2 := norm_mul_self a
  rw [normalize_def]
  have : norm2 (smul (1 / norm a) a) = norm2 a / (norm a * norm a) := by
    simp only [norm2_def, smul_x, smul_y, smul_z]; field_simp
  rw [this, h2, div_self ((norm2_pos_iff a).mpr h).ne']

open Ddg.Dev in
theorem norm_normalize (h : a ≠ zero) : norm (normalize a) = 1 := by
  rw [norm_def, norm2_normalize a h, sqrt_one]

open Ddg.Dev in
theorem normalize_coords :
    (normalize a) = smul (1 / norm a) a := normalize_def a

end algebra

/-! ## Derivatives along curves -/

/-- `F : ℝ → ℝ³` has the derivative `F'` at `t`: each coordinate has. -/
structure HasDerivAt3 (F : ℝ → V3) (F' : V3) (t : ℝ) : Prop where
  x : HasDerivAt (fun s => (F s).x) F'.x t
  y : HasDerivAt (fun s => (F s).y) F'.y t
  z : HasDerivAt (fun s => (F s).z) F'.z t

namespace HasDerivAt3

open Ddg.Dev

variable {F G : ℝ → V3} {F' G' : V3} {t : ℝ}

theorem const (p : V3) (t : ℝ) : HasDerivAt3 (fun _ => p) zero t :=
  ⟨by convert hasDerivAt_const t p.x using 1; simp,
    by convert hasDerivAt_const t p.y using 1; simp,
    by convert hasDerivAt_const t p.z using 1; simp⟩

/-- The line `s ↦ p + s δ`. -/
theorem line (p δ : V3) (t : ℝ) : HasDerivAt3 (fun s => add p (smul s δ)) δ t := by
  refine ⟨?_, ?_, ?_⟩ <;> simp only [add_x, add_y, add_z, smul_x, smul_y, smul_z]
  · simpa using ((hasDerivAt_id t).mul_const δ.x).const_add p.x
  · simpa using ((hasDerivAt_id t).mul_const δ.y).const_add p.y
  · simpa using ((hasDerivAt_id t).mul_const δ.z).const_add p.z

theorem add (hF : HasDerivAt3 F F' t) (hG : HasDerivAt3 G G' t) :
    HasDerivAt3 (fun s => Ddg.Dev.add (F s) (G s)) (Ddg.Dev.add F' G') t :=
  ⟨hF.x.add hG.x, hF.y.add hG.y, hF.z.add hG.z⟩

theorem sub (hF : HasDerivAt3 F F' t) (hG : HasDerivAt3 G G' t) :
    HasDerivAt3 (fun s => Ddg.Dev.sub (F s) (G s)) (Ddg.Dev.sub F' G') t :=
  ⟨hF.x.sub hG.x, hF.y.sub hG.y, hF.z.sub hG.z⟩

theorem neg (hF : HasDerivAt3 F F' t) :
    HasDerivAt3 (fun s => Ddg.Dev.neg (F s)) (Ddg.Dev.neg F') t :=
  ⟨hF.x.neg, hF.y.neg, hF.z.neg⟩

theorem smul {c : ℝ → ℝ} {c' : ℝ} (hc : HasDerivAt c c' t) (hF : HasDerivAt3 F F' t) :
    HasDerivAt3 (fun s => Ddg.Dev.smul (c s) (F s))
      (Ddg.Dev.add (Ddg.Dev.smul c' (F t)) (Ddg.Dev.smul (c t) F')) t :=
  ⟨hc.mul hF.x, hc.mul hF.y, hc.mul hF.z⟩

theorem smul_const (c : ℝ) (hF : HasDerivAt3 F F' t) :
    HasDerivAt3 (fun s => Ddg.Dev.smul c (F s)) (Ddg.Dev.smul c F') t :=
⟨hF.x.const_mul c, hF.y.const_mul c, hF.z.const_mul c⟩

theorem cross (hF : HasDerivAt3 F F' t) (hG : HasDerivAt3 G G' t) :
    HasDerivAt3 (fun s => Ddg.Dev.cross (F s) (G s))
      (Ddg.Dev.add (Ddg.Dev.cross F' (G t)) (Ddg.Dev.cross (F t) G')) t := by
  refine ⟨?_, ?_, ?_⟩
  · convert (hF.y.mul hG.z).sub (hF.z.mul hG.y) using 1 <;> first | rfl | (simp; ring)
  · convert (hF.z.mul hG.x).sub (hF.x.mul hG.z) using 1 <;> first | rfl | (simp; ring)
  · convert (hF.x.mul hG.y).sub (hF.y.mul hG.x) using 1 <;> first | rfl | (simp; ring)

theorem dot (hF : HasDerivAt3 F F' t) (hG : HasDerivAt3 G G' t) :
    HasDerivAt (fun s => Ddg.Dev.dot (F s) (G s))
      (Ddg.Dev.dot F' (G t) + Ddg.Dev.dot (F t) G') t := by
  have := ((hF.x.mul hG.x).add (hF.y.mul hG.y)).add (hF.z.mul hG.z)
  convert this using 1 <;> first | rfl | (simp only [dot_def]; ring)

theorem dot_const (w : V3) (hF : HasDerivAt3 F F' t) :
    HasDerivAt (fun s => Ddg.Dev.dot w (F s)) (Ddg.Dev.dot w F') t := by
  have := (const w t).dot hF
  convert this using 1
  simp [dot_def]

theorem norm2 (hF : HasDerivAt3 F F' t) :
    HasDerivAt (fun s => Ddg.Dev.norm2 (F s)) (2 * Ddg.Dev.dot (F t) F') t := by
  have := hF.dot hF
  convert this using 1 <;> first | rfl | (rw [dot_comm F' (F t)]; ring)

theorem norm (hF : HasDerivAt3 F F' t) (h0 : F t ≠ zero) :
    HasDerivAt (fun s => Ddg.Dev.norm (F s)) (Ddg.Dev.dot (F t) F' / Ddg.Dev.norm (F t)) t := by
  have h2 := hF.norm2
  have hpos := (norm2_pos_iff (F t)).mpr h0
  have := h2.sqrt hpos.ne'
  convert this using 1 <;> first | rfl | (rw [norm_def]; field_simp)

/-- Derivative of `F/|F|`: `(F' − (N·F') N)/|F|`, `N = F/|F|`. -/
theorem normalize (hF : HasDerivAt3 F F' t) (h0 : F t ≠ zero) :
    HasDerivAt3 (fun s => Ddg.Dev.normalize (F s))
      (Ddg.Dev.smul (1 / Ddg.Dev.norm (F t))
        (Ddg.Dev.sub F' (Ddg.Dev.smul (Ddg.Dev.dot (Ddg.Dev.normalize (F t)) F')
          (Ddg.Dev.normalize (F t))))) t := by
  have hn := hF.norm h0
  have hpos := (norm_pos_iff' (F t)).mpr h0
  have hinv : HasDerivAt (fun s => 1 / Ddg.Dev.norm (F s))
      (-(Ddg.Dev.dot (F t) F' / Ddg.Dev.norm (F t)) / Ddg.Dev.norm (F t) ^ 2) t := by
    have := (hasDerivAt_const t (1 : ℝ)).div hn hpos.ne'
    convert this using 1; ring
  have := HasDerivAt3.smul hinv hF
  simp only [normalize_def]
  convert this using 1
  ext <;> simp [dot_def] <;> field_simp <;> ring

theorem of_eq {F' F'' : V3} (h : HasDerivAt3 F F' t) (e : F' = F'') : HasDerivAt3 F F'' t :=
  e ▸ h

end HasDerivAt3

/-! ## `atan2` -/

section atan2

open Ddg.Dev

theorem atan2_def (y x : ℝ) : Ddg.Dev.atan2 y x =
    if 0 < x then arctan (y / x)
    else if x < 0 then (if y < 0 then arctan (y / x) - π else arctan (y / x) + π)
    else if 0 < y then π / 2 else if y < 0 then -(π / 2) else 0 := by
  simp [Ddg.Dev.atan2]

/-- For `y > 0`: `atan2(y, x) = π/2 − arctan(x/y)`. -/
theorem atan2_of_pos {y : ℝ} (hy : 0 < y) (x : ℝ) : Ddg.Dev.atan2 y x = π / 2 - arctan (x / y) := by
  rw [atan2_def]
  rcases lt_trichotomy x 0 with hx | rfl | hx
  · have hxy : y / x < 0 := div_neg_of_pos_of_neg hy hx
    have := arctan_inv_of_neg hxy
    rw [inv_div] at this
    rw [ite_eq_right (not_lt.mpr hx.le), ite_eq_left hx, ite_eq_right (not_lt.mpr hy.le), this]; ring
  · simp [hy]
  · have hxy : 0 < y / x := div_pos hy hx
    have := arctan_inv_of_pos hxy
    rw [inv_div] at this
    rw [ite_eq_left hx, this]; ring

theorem atan2_zero_of_pos {x : ℝ} (hx : 0 < x) : Ddg.Dev.atan2 0 x = 0 := by
  rw [atan2_def, ite_eq_left hx]; simp

theorem atan2_zero_of_neg {x : ℝ} (hx : x < 0) : Ddg.Dev.atan2 0 x = π := by
  rw [atan2_def, ite_eq_right (not_lt.mpr hx.le), ite_eq_left hx]; simp

/-- `cos (atan2 y x) = x/√(x² + y²)` for `y ≥ 0`, `(x, y) ≠ 0`. -/
theorem cos_atan2 {y x : ℝ} (hy : 0 ≤ y) (h : 0 < x ^ 2 + y ^ 2) :
    cos (Ddg.Dev.atan2 y x) = x / √(x ^ 2 + y ^ 2) := by
  rcases hy.lt_or_eq with hy | rfl
  · rw [atan2_of_pos hy, cos_pi_div_two_sub, sin_arctan]
    have : 1 + (x / y) ^ 2 = (x ^ 2 + y ^ 2) / y ^ 2 := by field_simp; ring
    rw [this, sqrt_div' _ (sq_nonneg y), sqrt_sq hy.le]
    field_simp
  · rcases lt_trichotomy x 0 with hx | rfl | hx
    · rw [atan2_zero_of_neg hx, cos_pi]
      have : √(x ^ 2 + 0 ^ 2) = -x := by
        rw [show x ^ 2 + 0 ^ 2 = (-x) ^ 2 by ring, sqrt_sq (by linarith)]
      rw [this, div_neg, div_self hx.ne]
    · simp at h
    · rw [atan2_zero_of_pos hx, cos_zero]
      have : √(x ^ 2 + 0 ^ 2) = x := by
        rw [show x ^ 2 + 0 ^ 2 = x ^ 2 by ring, sqrt_sq hx.le]
      rw [this, div_self hx.ne']

/-- `sin (atan2 y x) = y/√(x² + y²)` for `y ≥ 0`, `(x, y) ≠ 0`. -/
theorem sin_atan2 {y x : ℝ} (hy : 0 ≤ y) (h : 0 < x ^ 2 + y ^ 2) :
    sin (Ddg.Dev.atan2 y x) = y / √(x ^ 2 + y ^ 2) := by
  rcases hy.lt_or_eq with hy | rfl
  · rw [atan2_of_pos hy, sin_pi_div_two_sub, cos_arctan]
    have : 1 + (x / y) ^ 2 = (x ^ 2 + y ^ 2) / y ^ 2 := by field_simp; ring
    rw [this, sqrt_div' _ (sq_nonneg y), sqrt_sq hy.le]
    have : 0 < √(x ^ 2 + y ^ 2) := sqrt_pos.mpr h
    field_simp
  · rcases lt_trichotomy x 0 with hx | rfl | hx
    · rw [atan2_zero_of_neg hx, sin_pi]; simp
    · simp at h
    · rw [atan2_zero_of_pos hx, sin_zero]; simp

/-- `atan2 y x ∈ [0, π]` for `y ≥ 0`. -/
theorem atan2_mem_Icc {y : ℝ} (hy : 0 ≤ y) (x : ℝ) :
    0 ≤ Ddg.Dev.atan2 y x ∧ Ddg.Dev.atan2 y x ≤ π := by
  rcases hy.lt_or_eq with hy | rfl
  · rw [atan2_of_pos hy]
    have h1 := arctan_lt_pi_div_two (x / y)
    have h2 := neg_pi_div_two_lt_arctan (x / y)
    constructor <;> linarith
  · rw [atan2_def]
    split_ifs <;> simp_all [pi_pos.le]

/-- `atan2 (c y) (c x) = atan2 y x` for `c > 0`. -/
theorem atan2_smul {c : ℝ} (hc : 0 < c) (y x : ℝ) :
    Ddg.Dev.atan2 (c * y) (c * x) = Ddg.Dev.atan2 y x := by
  have e : c * y / (c * x) = y / x := by
    rcases eq_or_ne x 0 with rfl | hx
    · simp
    · field_simp
  rw [atan2_def, atan2_def, e]
  simp only [mul_pos_iff_of_pos_left hc, mul_neg_iff, hc, true_and, not_lt_of_gt hc, false_and,
    or_false]

/-- **Derivative of `atan2`** along a curve with `y > 0`: `(x y' − y x')/(x² + y²)`. -/
theorem HasDerivAt.atan2 {x y : ℝ → ℝ} {x' y' t : ℝ} (hx : HasDerivAt x x' t)
    (hy : HasDerivAt y y' t) (hpos : 0 < y t) :
    HasDerivAt (fun s => Ddg.Dev.atan2 (y s) (x s))
      ((x t * y' - y t * x') / (x t ^ 2 + y t ^ 2)) t := by
  have hev : ∀ᶠ s in 𝓝 t, 0 < y s := hy.continuousAt.eventually (lt_mem_nhds hpos)
  have heq : (fun s => Ddg.Dev.atan2 (y s) (x s)) =ᶠ[𝓝 t]
      (fun s => π / 2 - arctan (x s / y s)) := by
    filter_upwards [hev] with s hs
    rw [atan2_of_pos hs]
  have hd := ((hx.div hy hpos.ne').arctan).const_sub (π / 2)
  simp only [Pi.div_apply] at hd
  refine (hd.congr_of_eventuallyEq heq).congr_deriv ?_
  have h1 : y t ^ 2 ≠ 0 := by positivity
  have h2 : x t ^ 2 + y t ^ 2 ≠ 0 := by positivity
  field_simp
  ring

end atan2

end Dev

end DDG
