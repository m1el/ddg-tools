import DdgProofs.Delaunay.Rippa
import DDG.Laplace.Cotan

/-!
# Delaunay triangulations, layer 1: potentials that decrease under a Delaunay flip

Two real-valued functions of the edge lengths of a triangle, summed over the triangles of a
triangulation, decrease strictly with every flip of an edge that is not locally Delaunay and do
not change when an edge with `θ_k + θ_l = π` is flipped. Only the two triangles of the diamond
change, so each statement is about one diamond `(a, b, c, d, e)` and its flip `(f, c, e, b, d)`
(notation of `DdgProofs/Delaunay/Local.lean`).

## The lifting energy `A (a² + b² + c²)`

`liftEnergy a b c = A (a² + b² + c²)`, `A = triArea a b c`. On a planar triangle, `liftEnergy/12`
is the integral of `u_I − u` for `u(x) = |x|²` and `u_I` its linear interpolant (Chen–Holst 2011,
Lemma 4.1; proved in `DdgProofs/Delaunay/ODT.lean`, `it-odt-interpolation-error`): the
volume between the paraboloid and the triangle lifted to it. For a planar point set, the sum over
a triangulation is the volume under the lifted triangulation, minimal for the lower convex hull,
i.e. the Delaunay triangulation (Chen–Xu 2004; Musin 1997). It is intrinsic: it depends on the
lengths only, so it is defined for a triangulation of any piecewise flat surface.

* `liftEnergy_flip` — **the exact change under a flip**: for a strictly convex diamond,
  `ΔP = P_before − P_after = −8 A_ijk A_jil (cot θ_k + cot θ_l)`.
* `liftEnergy_flip_lt` (strict decrease for a non-Delaunay edge), `liftEnergy_flip_eq_of_tie`,
  `liftEnergy_flip_lt_iff` (for a strictly convex diamond the flip lowers it exactly when the edge
  is not locally Delaunay).

This is the potential to use for termination: the identity has no new-triangle denominators,
and the change is the product of the two old areas and the cotan weight.

## Musin's harmonic index `(a² + b² + c²)/A`

`harmonicIndex a b c = (a² + b² + c²)/A = 4 (cot α + cot β + cot γ)`: the trace of the cotan
Laplacian of a triangulation is `Σ_T harmonicIndex T / 4` (each corner cotangent is counted at
the two ends of its opposite edge, halved). B–S (Prop. 12 with Rippa's theorem) use it for the
termination of the flip algorithm on piecewise flat surfaces. Here it comes from Rippa's identity
summed over the four hat functions of the diamond (`harmonicIndex_flip`):
`ΔH = 4 κ (ℓ(e_i)² + ℓ(e_j)² + ℓ(e_k)² + ℓ(e_l)²)`, `κ = −(cot θ_k + cot θ_l)/(2 A_ilk A_jkl)`;
`harmonicIndex_flip_lt`, `harmonicIndex_flip_eq_of_tie`.

## References

* O. R. Musin, *Properties of the Delaunay triangulation*, Proc. 13th Annual Symposium on
  Computational Geometry (SoCG '97), ACM, 1997 (the harmonic index and other functionals
  minimized by the Delaunay triangulation; not read here, cited as B–S cite it).
* A. I. Bobenko, B. A. Springborn, *A discrete Laplace–Beltrami operator for simplicial
  surfaces*, Discrete Comput. Geom. 38(4), 740–756, 2007, Sec. 2 (the harmonic index on
  piecewise flat surfaces and termination of the flip algorithm, Prop. 12, Thm. 15 of arXiv v3).
* L. Chen, J. Xu, *Optimal Delaunay triangulations*, J. Comput. Math. 22(2), 299–308, 2004;
  L. Chen, M. Holst, *Efficient mesh optimization schemes based on Optimal Delaunay
  Triangulations*, CMAME 200(9–12), 967–984, 2011, Lemma 4.1 (the interpolation error of `|x|²`).
* C. Indermitte, Th. M. Liebling, M. Troyanov, H. Clémençon, Theoret. Comput. Sci. 263, 2001,
  Sec. 4, Thm. 1 (termination by a third potential, the sum of the areas of the circumcircles;
  not formalized here).
-/

noncomputable section

namespace DDG

namespace Delaunay

open Real

variable {a b c d e : ℝ}

/-! ## The lifting energy -/

/-- The lifting energy of a triangle with sides `a, b, c`: `A (a² + b² + c²)` (twelve times the
interpolation error of `|x|²` over the triangle). -/
def liftEnergy (a b c : ℝ) : ℝ := triArea a b c * (a ^ 2 + b ^ 2 + c ^ 2)

theorem liftEnergy_pos (h : IsTriangle a b c) : 0 < liftEnergy a b c := by
  have := triArea_pos h
  have := h.pos_a
  unfold liftEnergy
  positivity

/-- **The change of the lifting energy under a flip.** For a strictly convex diamond with
flipped length `f`:
`(P_ijk + P_jil) − (P_ilk + P_jkl) = −8 A_ijk A_jil (cot θ_k + cot θ_l)`.

[formula: it-delaunay-potential] -/
theorem liftEnergy_flip (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    liftEnergy a b c + liftEnergy a d e -
        (liftEnergy (flipLength a b c d e) b d + liftEnergy (flipLength a b c d e) e c) =
      -8 * triArea a b c * triArea a d e * (cornerCot a b c + cornerCot a d e) := by
  obtain ⟨X, H, Y, G, -, -, -, -, eb, ec, ed, ee, ef, ak, al, an1, an2⟩ :=
    layout_coords hk hl hconv
  have hk4 : 4 * triArea a b c * cornerCot a b c = b ^ 2 + c ^ 2 - a ^ 2 := by
    have := (triArea_pos hk).ne'
    unfold cornerCot; field_simp
  have hl4 : 4 * triArea a d e * cornerCot a d e = d ^ 2 + e ^ 2 - a ^ 2 := by
    have := (triArea_pos hl).ne'
    unfold cornerCot; field_simp
  have hr : -8 * triArea a b c * triArea a d e * (cornerCot a b c + cornerCot a d e) =
      -2 * triArea a d e * (b ^ 2 + c ^ 2 - a ^ 2) - 2 * triArea a b c * (d ^ 2 + e ^ 2 - a ^ 2) := by
    linear_combination (-2 * triArea a d e) * hk4 + (-2 * triArea a b c) * hl4
  rw [hr]
  unfold liftEnergy
  generalize flipLength a b c d e = f at *
  rw [eb, ec, ed, ee, ef, show triArea a b c = a * H / 2 by linarith,
    show triArea a d e = a * G / 2 by linarith, show triArea f b d = (X * G + Y * H) / 2 by linarith,
    show triArea f e c = ((a - X) * G + (a - Y) * H) / 2 by linarith]
  ring

/-- **The lifting energy decreases strictly with the flip of a non-Delaunay edge.**

[formula: it-delaunay-potential] -/
theorem liftEnergy_flip_lt (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : ¬ IsLocallyDelaunay a b c d e) :
    liftEnergy (flipLength a b c d e) b d + liftEnergy (flipLength a b c d e) e c <
      liftEnergy a b c + liftEnergy a d e := by
  have hconv := isConvexDiamond_of_not_isLocallyDelaunay hk hl h
  have hid := liftEnergy_flip hk hl hconv
  have hc := (pi_lt_iff_cot_neg hk hl).mp (not_isLocallyDelaunay_iff.mp h)
  have hA := triArea_pos hk
  have hB := triArea_pos hl
  have : 0 < -8 * triArea a b c * triArea a d e * (cornerCot a b c + cornerCot a d e) := by
    have : 0 < triArea a b c * triArea a d e * -(cornerCot a b c + cornerCot a d e) :=
      mul_pos (mul_pos hA hB) (by linarith)
    linarith
  linarith

/-- At a tie the flip does not change the lifting energy. -/
theorem liftEnergy_flip_eq_of_tie (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle a b c + cornerAngle a d e = π) :
    liftEnergy (flipLength a b c d e) b d + liftEnergy (flipLength a b c d e) e c =
      liftEnergy a b c + liftEnergy a d e := by
  have hid := liftEnergy_flip hk hl (isConvexDiamond_of_pi_le hk hl h.ge)
  rw [(eq_pi_iff_cot_eq_zero hk hl).mp h, mul_zero] at hid
  linarith

/-- For a strictly convex diamond the flip lowers the lifting energy exactly when the edge is
not locally Delaunay (and raises it when the edge is strictly Delaunay). -/
theorem liftEnergy_flip_lt_iff (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    liftEnergy (flipLength a b c d e) b d + liftEnergy (flipLength a b c d e) e c <
        liftEnergy a b c + liftEnergy a d e ↔ ¬ IsLocallyDelaunay a b c d e := by
  refine ⟨fun hlt hD => ?_, liftEnergy_flip_lt hk hl⟩
  have hid := liftEnergy_flip hk hl hconv
  have hc := (isLocallyDelaunay_iff_cot_nonneg hk hl).mp hD
  have hA := triArea_pos hk
  have hB := triArea_pos hl
  have : 0 ≤ triArea a b c * triArea a d e * (cornerCot a b c + cornerCot a d e) :=
    mul_nonneg (mul_pos hA hB).le hc
  linarith

/-! ## Musin's harmonic index -/

/-- Musin's harmonic index of a triangle: `(a² + b² + c²)/A`. -/
def harmonicIndex (a b c : ℝ) : ℝ := (a ^ 2 + b ^ 2 + c ^ 2) / triArea a b c

/-- The harmonic index is four times the sum of the cotangents of the three corners. -/
theorem harmonicIndex_eq_cot (h : IsTriangle a b c) :
    harmonicIndex a b c = 4 * (cornerCot a b c + cornerCot b c a + cornerCot c a b) := by
  have hT := (triArea_pos h).ne'
  unfold harmonicIndex cornerCot
  rw [triArea_rotate, triArea_rotate b c a]
  field_simp
  ring

/-- The cotan energies of the four hat functions of a diamond add up to the sums of the corner
cotangents of its two triangles. -/
theorem sum_diamondEnergy_hat (a b c d e : ℝ) :
    diamondEnergy a b c d e 1 0 0 0 + diamondEnergy a b c d e 0 1 0 0 +
        diamondEnergy a b c d e 0 0 1 0 + diamondEnergy a b c d e 0 0 0 1 =
      (cornerCot a b c + cornerCot b c a + cornerCot c a b) +
        (cornerCot a d e + cornerCot d e a + cornerCot e a d) := by
  unfold diamondEnergy cotEnergy
  ring

/-- **The change of the harmonic index under a flip**, from Rippa's identity at the four hat
functions: with `κ = −(cot θ_k + cot θ_l)/(2 A_ilk A_jkl)`,
`ΔH = 4 κ (A_jkl² + A_ilk² + A_jil² + A_ijk²)`.

[formula: it-delaunay-potential] -/
theorem harmonicIndex_flip (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    harmonicIndex a b c + harmonicIndex a d e -
        (harmonicIndex (flipLength a b c d e) b d + harmonicIndex (flipLength a b c d e) c e) =
      4 * (-(cornerCot a b c + cornerCot a d e) /
          (2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e)) *
        (triArea (flipLength a b c d e) c e ^ 2 + triArea (flipLength a b c d e) b d ^ 2 +
          triArea a d e ^ 2 + triArea a b c ^ 2) := by
  obtain ⟨t1, t2⟩ := flipLength_isTriangle hk hl hconv
  have r1 := rippa_identity hk hl hconv 1 0 0 0
  have r2 := rippa_identity hk hl hconv 0 1 0 0
  have r3 := rippa_identity hk hl hconv 0 0 1 0
  have r4 := rippa_identity hk hl hconv 0 0 0 1
  have s1 := sum_diamondEnergy_hat a b c d e
  have s2 := sum_diamondEnergy_hat (flipLength a b c d e) c e b d
  rw [harmonicIndex_eq_cot hk, harmonicIndex_eq_cot hl, harmonicIndex_eq_cot t1,
    harmonicIndex_eq_cot t2.swap]
  unfold rippaForm at r1 r2 r3 r4
  generalize -(cornerCot a b c + cornerCot a d e) /
      (2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e) = κ at *
  linear_combination 4 * r1 + 4 * r2 + 4 * r3 + 4 * r4 - 4 * s1 + 4 * s2

/-- **The harmonic index decreases strictly with the flip of a non-Delaunay edge** (B–S,
Prop. 12).

[formula: it-delaunay-potential] -/
theorem harmonicIndex_flip_lt (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : ¬ IsLocallyDelaunay a b c d e) :
    harmonicIndex (flipLength a b c d e) b d + harmonicIndex (flipLength a b c d e) c e <
      harmonicIndex a b c + harmonicIndex a d e := by
  have hconv := isConvexDiamond_of_not_isLocallyDelaunay hk hl h
  obtain ⟨t1, t2⟩ := flipLength_isTriangle hk hl hconv
  have hid := harmonicIndex_flip hk hl hconv
  have hc := (pi_lt_iff_cot_neg hk hl).mp (not_isLocallyDelaunay_iff.mp h)
  have hA := triArea_pos hk
  have hD : 0 < 2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e := by
    have := triArea_pos t1
    have := triArea_pos t2.swap
    positivity
  have hκ : 0 < -(cornerCot a b c + cornerCot a d e) /
      (2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e) :=
    div_pos (by linarith) hD
  have hS : 0 < triArea (flipLength a b c d e) c e ^ 2 + triArea (flipLength a b c d e) b d ^ 2 +
      triArea a d e ^ 2 + triArea a b c ^ 2 := by positivity
  have := mul_pos (mul_pos (by norm_num : (0 : ℝ) < 4) hκ) hS
  linarith

/-- At a tie the flip does not change the harmonic index. -/
theorem harmonicIndex_flip_eq_of_tie (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle a b c + cornerAngle a d e = π) :
    harmonicIndex (flipLength a b c d e) b d + harmonicIndex (flipLength a b c d e) c e =
      harmonicIndex a b c + harmonicIndex a d e := by
  have hid := harmonicIndex_flip hk hl (isConvexDiamond_of_pi_le hk hl h.ge)
  rw [(eq_pi_iff_cot_eq_zero hk hl).mp h, neg_zero, zero_div, mul_zero, zero_mul] at hid
  linarith

end Delaunay

end DDG
