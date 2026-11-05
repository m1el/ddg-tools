import DDG.Geometry.Plane

/-!
# Plane lemmas: distances and signed areas of points of `ℝ × ℝ`

Small facts about the plane primitives `dist2`, `sda` and `triArea` of
`DDG/Geometry/Plane.lean` and `DDG/Geometry/Triangle.lean`, used by the intrinsic
triangulations, the Delaunay theory and the finite element integrals:

* `dist2` is symmetric and nonnegative, and its square written out (`dist2_comm`,
  `dist2_nonneg`, `dist2_sq'`);
* a plane triangle of nonzero signed area has distinct corners and its sides satisfy the strict
  triangle inequalities (`dist2_pos_of_sda_ne`, `isTriangle_of_sda_ne`), and its Heron area is
  half the absolute signed double area (`two_mul_triArea_dist2`);
* the quarter turn `rot`, the squared norm `nsq` and the linear form `toDual G` of a vector
  (`toDual_apply`).

(Library step S2b: moved here from `DDG/Intrinsic/Diamond.lean` and, `rot`, `nsq`,
`toDual`, from `DDG/Intrinsic/Move.lean`; the names keep their
`DDG.Intrinsic` namespace until the renames of step S4.)
-/

noncomputable section

namespace DDG

namespace Intrinsic

open Real

theorem dist2_comm (p q : ℝ × ℝ) : dist2 p q = dist2 q p := by
  unfold dist2
  congr 1
  ring

theorem dist2_nonneg (p q : ℝ × ℝ) : 0 ≤ dist2 p q := sqrt_nonneg _

/-- `dist2 p q ^ 2` written out. -/
theorem dist2_sq' (p q : ℝ × ℝ) : dist2 p q ^ 2 = (q.1 - p.1) ^ 2 + (q.2 - p.2) ^ 2 :=
  dist2_sq p q

/-- A triangle of nonzero signed area has two different first corners. -/
theorem dist2_pos_of_sda_ne {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : 0 < dist2 a b := by
  rcases (dist2_nonneg a b).lt_or_eq with hpos | hzero
  · exact hpos
  · exfalso
    have h2 := dist2_sq' a b
    rw [← hzero] at h2
    have h1 : (b.1 - a.1) ^ 2 = 0 := by nlinarith [sq_nonneg (b.1 - a.1), sq_nonneg (b.2 - a.2)]
    have h3 : (b.2 - a.2) ^ 2 = 0 := by nlinarith [sq_nonneg (b.1 - a.1), sq_nonneg (b.2 - a.2)]
    have e1 : b.1 - a.1 = 0 := by simpa using h1
    have e2 : b.2 - a.2 = 0 := by simpa using h3
    apply h
    unfold sda
    rw [e1, e2]
    ring

/-- **The sides of a plane triangle of nonzero area satisfy the strict triangle
inequalities.** -/
theorem isTriangle_of_sda_ne {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) :
    IsTriangle (dist2 b c) (dist2 c a) (dist2 a b) := by
  apply isTriangle_of_heron_pos (dist2_nonneg _ _) (dist2_nonneg _ _) (dist2_nonneg _ _)
  have key : (dist2 b c + dist2 c a + dist2 a b) * (dist2 c a - dist2 b c + dist2 a b) *
      (dist2 b c - dist2 c a + dist2 a b) * (dist2 b c - dist2 a b + dist2 c a) =
      (2 * sda a b c) ^ 2 := by
    calc _ = 2 * dist2 b c ^ 2 * dist2 c a ^ 2 + 2 * dist2 c a ^ 2 * dist2 a b ^ 2 +
          2 * dist2 a b ^ 2 * dist2 b c ^ 2 - (dist2 b c ^ 2) ^ 2 - (dist2 c a ^ 2) ^ 2 -
          (dist2 a b ^ 2) ^ 2 := by ring
      _ = _ := by
        rw [dist2_sq, dist2_sq, dist2_sq]
        linear_combination (-1 : ℝ) * sda_sq_mul_four a b c
  rw [key]
  positivity

/-- The Heron area of the sides of a plane triangle, in the order of the oracle's diamond
(`a = |ij|`, `b = |ik|`, `c = |jk|`), is half the absolute signed double area. -/
theorem two_mul_triArea_dist2 (i j k : ℝ × ℝ) :
    2 * triArea (dist2 i j) (dist2 i k) (dist2 j k) = |sda i j k| := by
  rw [abs_sda_eq_two_mul_triArea i j k, dist2_comm k i]
  unfold triArea
  congr 3
  ring

/-! ## The quarter turn, the squared norm, a vector as a linear form -/

/-- The quarter turn `rot u = (−u₂, u₁)`. -/
def rot (u : ℝ × ℝ) : ℝ × ℝ := (-u.2, u.1)

/-- The squared norm `|u|² = u₁² + u₂²`. -/
def nsq (u : ℝ × ℝ) : ℝ := u.1 ^ 2 + u.2 ^ 2

/-- The linear form `v ↦ G₁ v₁ + G₂ v₂` (the gradient `G` as a derivative). -/
def toDual (G : ℝ × ℝ) : (ℝ × ℝ) →L[ℝ] ℝ :=
  G.1 • ContinuousLinearMap.fst ℝ ℝ ℝ + G.2 • ContinuousLinearMap.snd ℝ ℝ ℝ

@[simp] theorem toDual_apply (G v : ℝ × ℝ) : toDual G v = G.1 * v.1 + G.2 * v.2 := by
  simp [toDual]

end Intrinsic

end DDG
