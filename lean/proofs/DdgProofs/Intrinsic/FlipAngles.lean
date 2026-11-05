import DdgProofs.Intrinsic.Diamond

/-!
# Intrinsic triangulations: what a flip does to the corner angles

`crates/ddg-intrinsic/SPEC.md` §4 (Delaunay flipping), §5 item 2; the paper Sec. 4.1.

Notation of `DdgProofs/Intrinsic/Diamond.lean`: an edge `ij` of length `a` with the triangles
`ijk` (`b = ℓ_ik`, `c = ℓ_jk`) and `jil` (`d = ℓ_il`, `e = ℓ_jl`); the angles opposite the edge
are `θ_k = cornerAngle a b c` and `θ_l = cornerAngle a d e`. After a flip the edge `kl` of length
`f` has the diamond `(f, c, e, b, d)` (apexes `j` and `i`).

* the angles of the quadrilateral at the two ends of the edge sum to `2π − (θ_k + θ_l)`
  (`quad_angle_sum`); the diamond of the flipped edge is strictly convex too (`flip_convex`);
* **`it-angle-normalization`** — a flip keeps the angle sum at each of its four vertices: the
  new corner at `i` is the sum of the two old ones, the old corner at `k` is the sum of the two
  new ones, and likewise at `j` and `l` (`flip_angle_i`, `flip_angle_j`, `flip_angle_sums`).

The Delaunay criterion of one diamond (the cotan weight, the excess, the flip of a violating edge)
is in `DdgProofs/Delaunay/Local.lean`.
-/

noncomputable section

namespace DDG

namespace Intrinsic

open Real

/-- The angles of the quadrilateral at the two ends of the edge sum to `2π − (θ_k + θ_l)`. -/
theorem quad_angle_sum {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e) :
    (cornerAngle c a b + cornerAngle e a d) + (cornerAngle b c a + cornerAngle d e a) =
      2 * π - (cornerAngle a b c + cornerAngle a d e) := by
  have h1 := cornerAngle_add_cornerAngle_add_cornerAngle hk
  have h2 := cornerAngle_add_cornerAngle_add_cornerAngle hl
  linarith

/-- The new corner at `i` (in the triangle `ilk`) is the sum of the two old corners at `i`,
when that sum is at most `π`. -/
theorem flip_angle_i {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle c a b + cornerAngle e a d ≤ π) :
    cornerAngle (flippedLengthFromLengths a b c d e) b d =
      cornerAngle c a b + cornerAngle e a d := by
  have hb := hk.pos_b
  have hd := hl.pos_b
  have a1 := cornerAngle_pos hk.rotate.rotate
  have a3 := cornerAngle_pos hl.rotate.rotate
  have hsq := flippedLengthFromLengths_sq hk hl
  have e1 : (b ^ 2 + d ^ 2 - flippedLengthFromLengths a b c d e ^ 2) / (2 * b * d) =
      cos (cornerAngle c a b + cornerAngle e a d) := by
    rw [hsq]
    field_simp
    ring
  rw [show cornerAngle (flippedLengthFromLengths a b c d e) b d =
      arccos ((b ^ 2 + d ^ 2 - flippedLengthFromLengths a b c d e ^ 2) / (2 * b * d))
      from rfl, e1, Real.arccos_cos (by linarith) h]

/-- The new corner at `j` (in the triangle `jkl`) is the sum of the two old corners at `j`. -/
theorem flip_angle_j {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle b c a + cornerAngle d e a ≤ π) :
    cornerAngle (flippedLengthFromLengths a b c d e) c e =
      cornerAngle b c a + cornerAngle d e a := by
  have hc := hk.pos_c
  have he := hl.pos_c
  have b1 := cornerAngle_pos hk.rotate
  have b3 := cornerAngle_pos hl.rotate
  have hsq := flippedLengthFromLengths_sq' hk hl
  have e1 : (c ^ 2 + e ^ 2 - flippedLengthFromLengths a b c d e ^ 2) / (2 * c * e) =
      cos (cornerAngle b c a + cornerAngle d e a) := by
    rw [hsq]
    field_simp
    ring
  rw [show cornerAngle (flippedLengthFromLengths a b c d e) c e =
      arccos ((c ^ 2 + e ^ 2 - flippedLengthFromLengths a b c d e ^ 2) / (2 * c * e))
      from rfl, e1, Real.arccos_cos (by linarith) h]

/-- **The diamond of the flipped edge is strictly convex too** (it is the same quadrilateral):
with `f` the flipped length, `DiamondConvex f c e b d`. -/
theorem flip_convex {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    DiamondConvex (flippedLengthFromLengths a b c d e) c e b d := by
  have ha := hk.pos_a
  obtain ⟨hk0, hl0, -, e2, e3, e4, e5⟩ := diamond_isometric hk hl
  obtain ⟨c1, c2⟩ := diamondConvex_sda hconv
  have hJ : 0 < sda (apexK a b c) (apexL a d e) (a, 0) := by
    rw [sda_rotate]
    exact c2
  have hI : sda (apexK a b c) (apexL a d e) (0, 0) < 0 := by
    rw [sda_rotate, sda_swap_right]
    linarith
  obtain ⟨r1, r2⟩ := flipArea_realized hJ hI
  rw [dist2_apex, dist2_comm (apexK a b c) (a, 0), dist2_comm (apexL a d e) (a, 0),
    dist2_comm (apexK a b c) (0, 0), dist2_comm (apexL a d e) (0, 0), e2, e3, e4, e5] at r1 r2
  obtain ⟨s1, s2⟩ := diamond_sda a (diamondX a b c) (diamondH a b c)
    (diamondX a d e) (diamondH a d e)
  have t1 : sda (apexK a b c) (0, 0) (a, 0) = a * diamondH a b c := by
    rw [sda_rotate, sda_rotate]
    exact s1
  have t2 : sda (apexL a d e) (a, 0) (0, 0) = a * diamondH a d e := by
    rw [sda_rotate, sda_swap_right]
    unfold apexL
    rw [s2]
    ring
  unfold DiamondConvex
  rw [flipConvex_iff, r1, r2, t1, t2]
  exact ⟨mul_pos ha hk0, mul_pos ha hl0⟩

/-- **A flip keeps the angle sum at each of the four vertices of the diamond.** For a strictly
convex diamond with flipped length `f`: the corner of the new triangle `ilk` at `i` is the sum of
the two old corners at `i`; the corner of `jkl` at `j` is the sum of the two old corners at `j`;
the old corner at `k` is the sum of the corners of `jkl` and `ilk` at `k`; the old corner at `l`
is the sum of their corners at `l`. So the totals `Θ` that the angle normalization uses do not
change (SPEC §3, §5 item 2).

[formula: it-angle-normalization] -/
theorem flip_angle_sums {a b c d e : ℝ} (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : DiamondConvex a b c d e) :
    cornerAngle (flippedLengthFromLengths a b c d e) b d =
        cornerAngle c a b + cornerAngle e a d ∧
      cornerAngle (flippedLengthFromLengths a b c d e) c e =
        cornerAngle b c a + cornerAngle d e a ∧
      cornerAngle e (flippedLengthFromLengths a b c d e) c +
          cornerAngle d (flippedLengthFromLengths a b c d e) b = cornerAngle a b c ∧
      cornerAngle c e (flippedLengthFromLengths a b c d e) +
          cornerAngle b d (flippedLengthFromLengths a b c d e) = cornerAngle a d e := by
  obtain ⟨h1, h2⟩ := (diamondConvex_iff_angles hk hl).mp hconv
  obtain ⟨t1, t2⟩ := flip_isTriangle hk hl hconv
  have hk' := t2.swap
  have hconv' := flip_convex hk hl hconv
  obtain ⟨g1, g2⟩ := (diamondConvex_iff_angles hk' t1).mp hconv'
  have hback := flip_back hk hl hconv
  have k1 := flip_angle_i hk' t1 g1.le
  have k2 := flip_angle_j hk' t1 g2.le
  rw [hback] at k1 k2
  refine ⟨flip_angle_i hk hl h1.le, flip_angle_j hk hl h2.le, ?_, ?_⟩
  · rw [← k1, cornerAngle_swap]
  · rw [← k2, cornerAngle_swap]

end Intrinsic

end DDG
