import DdgProofs.Delaunay.Local
import DdgProofs.Delaunay.PlanarBasic

/-!
# Planar Delaunay triangulations: the bridge to the one-diamond theory (layer 1)

The planar theory (`DdgProofs/Delaunay/Planar*.lean`) is written in coordinates with the in-circle
determinant; layer 1 (`DdgProofs/Delaunay/Local.lean`) is written intrinsically, with the five
edge lengths of a diamond. For a diamond of the plane they are the same criterion:

* `inCircle_eq_dot` — `inCircle i j k l = ⟨i − k, j − k⟩ sda(i,j,l) − ⟨i − l, j − l⟩ sda(i,j,k)`;
* `cornerCot_dist2` — the cotangent of the corner at `k` from the three Euclidean lengths is
  `⟨i − k, j − k⟩ / |sda(i,j,k)|`;
* **`isLocallyDelaunay_dist2_iff`** — for `k, l` on opposite sides of `ij`, layer 1's
  `IsLocallyDelaunay |ij| |ik| |jk| |il| |jl|` (`θ_k + θ_l ≤ π`) holds iff `l` is not strictly
  inside the circumcircle of `i, j, k` (`¬ InCircumcircle i j k l`), via
  `sda(i,j,k) · inCircle i j k l = −|sda(i,j,k)|² |sda(i,j,l)| (cot θ_k + cot θ_l)`.

So the planar `IsLocallyDelaunay T` of `DdgProofs/Delaunay/PlanarFlip.lean` is layer 1's local
criterion at every interior edge of `T`, and Delaunay's lemma (`delaunay_lemma`) and the
termination of flipping (`exists_flips_to_locallyDelaunay`) are statements about it.
-/

noncomputable section

namespace DDG

namespace Delaunay

namespace Planar

open Real

/-- The in-circle determinant through two dot products. -/
theorem inCircle_eq_dot (i j k l : ℝ × ℝ) :
    inCircle i j k l = ((i.1 - k.1) * (j.1 - k.1) + (i.2 - k.2) * (j.2 - k.2)) * sda i j l -
      ((i.1 - l.1) * (j.1 - l.1) + (i.2 - l.2) * (j.2 - l.2)) * sda i j k := by
  unfold inCircle sda sqDist2; ring

/-- The Euclidean side lengths `|ij|, |ik|, |jk|` of a nondegenerate plane triangle form a
triangle in the sense of `DdgProofs/Triangle.lean`. -/
theorem isTriangle_dist2 {i j k : ℝ × ℝ} (h : sda i j k ≠ 0) :
    IsTriangle (dist2 i j) (dist2 i k) (dist2 j k) := by
  have h' : sda k i j ≠ 0 := by
    have : sda k i j = sda i j k := by unfold sda; ring
    rwa [this]
  have t := Intrinsic.isTriangle_of_sda_ne h'
  rw [Intrinsic.dist2_comm k i] at t
  exact ⟨by linarith [t.lt_a], by linarith [t.lt_c], by linarith [t.lt_b]⟩

/-- **The cotangent of a corner in coordinates**: `cot θ_k = ⟨i − k, j − k⟩ / |sda(i,j,k)|`. -/
theorem cornerCot_dist2 {i j k : ℝ × ℝ} (h : sda i j k ≠ 0) :
    cornerCot (dist2 i j) (dist2 i k) (dist2 j k) =
      ((i.1 - k.1) * (j.1 - k.1) + (i.2 - k.2) * (j.2 - k.2)) / |sda i j k| := by
  have hA : 4 * triArea (dist2 i j) (dist2 i k) (dist2 j k) = 2 * |sda i j k| := by
    have e := abs_sda_eq_two_mul_triArea k i j
    rw [Intrinsic.dist2_comm k i, show sda k i j = sda i j k by unfold sda; ring,
      triArea_swap] at e
    linarith
  have hα : 0 < |sda i j k| := abs_pos.2 h
  unfold cornerCot
  rw [hA, dist2_sq, dist2_sq, dist2_sq]
  unfold sqDist2
  field_simp
  ring

/-- **Layer 1's local Delaunay criterion is the planar one.** For `k` and `l` on opposite sides of
`ij`, the diamond with lengths `|ij|, |ik|, |jk|, |il|, |jl|` is locally Delaunay
(`θ_k + θ_l ≤ π`) iff `l` is not strictly inside the circumcircle of `i, j, k`.

[formula: planar-local-delaunay-bridge] -/
theorem isLocallyDelaunay_dist2_iff {i j k l : ℝ × ℝ} (hk : sda i j k * sda i j l < 0) :
    _root_.DDG.Delaunay.IsLocallyDelaunay (dist2 i j) (dist2 i k) (dist2 j k) (dist2 i l)
      (dist2 j l) ↔ ¬ InCircumcircle i j k l := by
  have hD : sda i j k ≠ 0 := by rintro h0; simp [h0] at hk
  have hD' : sda i j l ≠ 0 := by rintro h0; simp [h0] at hk
  rw [isLocallyDelaunay_iff_cot_nonneg (isTriangle_dist2 hD) (isTriangle_dist2 hD'),
    cornerCot_dist2 hD, cornerCot_dist2 hD', InCircumcircle, inCircle_eq_dot, not_lt]
  set dk := (i.1 - k.1) * (j.1 - k.1) + (i.2 - k.2) * (j.2 - k.2)
  set dl := (i.1 - l.1) * (j.1 - l.1) + (i.2 - l.2) * (j.2 - l.2)
  have hα : 0 < |sda i j k| := abs_pos.2 hD
  have hβ : 0 < |sda i j l| := abs_pos.2 hD'
  have h1 : |sda i j k| * |sda i j l| = -(sda i j k * sda i j l) := by
    rw [← abs_mul, abs_of_neg hk]
  have h2 : |sda i j k| * |sda i j k| = sda i j k * sda i j k := abs_mul_abs_self _
  have key : sda i j k * (dk * sda i j l - dl * sda i j k) =
      -|sda i j k| * (dk * |sda i j l| + dl * |sda i j k|) := by
    linear_combination dk * h1 + dl * h2
  have hsum : dk / |sda i j k| + dl / |sda i j l| =
      (dk * |sda i j l| + dl * |sda i j k|) / (|sda i j k| * |sda i j l|) := by
    field_simp
  rw [key, hsum, le_div_iff₀ (mul_pos hα hβ), zero_mul]
  constructor <;> intro h <;> nlinarith

end Planar

end Delaunay

end DDG
