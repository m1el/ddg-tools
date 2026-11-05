import DdgProofs.Delaunay.Local
import DdgProofs.Cotan
import DdgProofs.Laplace.Cotan

/-!
# Delaunay triangulations, layer 1: Rippa's theorem for one diamond

For vertex values `u_i, u_j, u_k, u_l` on the diamond `(a, b, c, d, e)` of
`DdgProofs/Delaunay/Local.lean`, the Dirichlet energy of the piecewise linear interpolant on the
two triangles `ijk`, `jil` is the cotan energy (`cotEnergy`, `DdgProofs/Cotan.lean`, which
`cotan_formula` identifies with `∫ |∇u|²` over each triangle)
`diamondEnergy a b c d e u_i u_j u_k u_l`. After the flip it is the energy of the new diamond
`(f, c, e, b, d)` with the values read in its order `(u_k, u_l, u_j, u_i)`.

## Main statements

* `rippa_identity` — **the difference is an exact rank-one form**: for a strictly convex diamond,
  `E_before − E_after = −(cot θ_k + cot θ_l) / (2 A_ilk A_jkl) · ℓ(u)²` with the linear form
  `ℓ(u) = A_jil u_k + A_ijk u_l − A_jkl u_i − A_ilk u_j` (the affine dependence of the four
  points of the layout; it vanishes on constants, `rippaForm_const`, and on all affine functions
  of the layout).
* `rippa_le` — Rippa's theorem for one diamond: if the old edge is not locally Delaunay
  (`θ_k + θ_l > π`), or at a tie, the flip does not increase the energy, for every `u`;
  `rippa_lt` — strictly, unless `ℓ(u) = 0`; `rippa_eq_of_tie` — at a tie the two energies agree
  for every `u`; `forall_rippa_le_iff` — conversely, for a strictly convex diamond the flip
  lowers or keeps the energy for every `u` exactly when `θ_k + θ_l ≥ π`.

The global statement (the Delaunay triangulation minimizes the Dirichlet energy over all
triangulations with the same vertices) needs the flip algorithm on a whole triangulation and
stays cited (`CITED.md`, `it-rippa-theorem`).

## References

* S. Rippa, *Minimal roughness property of the Delaunay triangulation*, Comput. Aided Geom.
  Design 7(6), 489–497, 1990 (the plane).
* A. I. Bobenko, B. A. Springborn, *A discrete Laplace–Beltrami operator for simplicial
  surfaces*, Discrete Comput. Geom. 38(4), 740–756, 2007, Sec. 2 ("Rippa's Theorem", between
  Definition 13 and Proposition 14 of arXiv v3): the same statement for piecewise flat surfaces,
  whose proof is local to one diamond.
* L. Chen, J. Xu, *Optimal Delaunay triangulations*, J. Comput. Math. 22(2), 299–308, 2004
  (the rank-one structure of the difference of two interpolants on a quadrilateral).

The identity is checked here as a polynomial identity in the coordinates of the planar layout
(`layout_coords`); the "rank ≤ 1" argument of the literature (both energies agree on affine
functions and the difference vanishes on them) is not needed.
-/

noncomputable section

namespace DDG

namespace Delaunay

open Real

variable {a b c d e : ℝ}

/-- The Dirichlet (cotan) energy of the values `u_i, u_j, u_k, u_l` on the two triangles `ijk`
(sides `a = ℓ_ij`, `b = ℓ_ik`, `c = ℓ_jk`) and `jil` (`d = ℓ_il`, `e = ℓ_jl`) of a diamond. -/
def diamondEnergy (a b c d e ui uj uk ul : ℝ) : ℝ :=
  cotEnergy a b c uk uj ui + cotEnergy a d e ul uj ui

/-- The linear form of Rippa's identity: `ℓ(u) = A_jil u_k + A_ijk u_l − A_jkl u_i − A_ilk u_j`,
with `A_jkl = triArea f c e`, `A_ilk = triArea f b d` the areas of the new triangles. -/
def rippaForm (a b c d e ui uj uk ul : ℝ) : ℝ :=
  triArea a d e * uk + triArea a b c * ul - triArea (flipLength a b c d e) c e * ui -
    triArea (flipLength a b c d e) b d * uj

/-- The energy after the flip, written out: the triangles `ilk` and `jkl`. -/
theorem diamondEnergy_flip_eq (f ui uj uk ul : ℝ) :
    diamondEnergy f c e b d uk ul uj ui = cotEnergy f b d ui ul uk + cotEnergy f e c uj uk ul := by
  unfold diamondEnergy cotEnergy cornerCot
  rw [triArea_swap e c f, triArea_swap f e c, triArea_swap c f e]
  ring

private theorem area_rot {x y z A : ℝ} (h : 2 * triArea x y z = A) :
    triArea x y z = A / 2 ∧ triArea y z x = A / 2 ∧ triArea z x y = A / 2 := by
  refine ⟨by linarith, ?_, ?_⟩
  · rw [triArea_rotate]; linarith
  · rw [triArea_rotate, triArea_rotate]; linarith

/-- **Rippa's identity for one diamond.** For a strictly convex diamond with flipped length
`f = flipLength a b c d e` and any vertex values,
`E_before(u) − E_after(u) = −(cot θ_k + cot θ_l) / (2 A_ilk A_jkl) · ℓ(u)²`.

[formula: it-rippa-diamond] -/
theorem rippa_identity (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) (ui uj uk ul : ℝ) :
    diamondEnergy a b c d e ui uj uk ul -
        diamondEnergy (flipLength a b c d e) c e b d uk ul uj ui =
      -(cornerCot a b c + cornerCot a d e) /
          (2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e) *
        rippaForm a b c d e ui uj uk ul ^ 2 := by
  obtain ⟨X, H, Y, G, hH, hG, h1, h2, eb, ec, ed, ee, ef, ak, al, an1, an2⟩ :=
    layout_coords hk hl hconv
  have ha := hk.pos_a
  rw [triArea_swap] at an2
  obtain ⟨k1, k2, k3⟩ := area_rot ak
  obtain ⟨l1, l2, l3⟩ := area_rot al
  obtain ⟨m1, m2, m3⟩ := area_rot an1
  obtain ⟨n1, n2, n3⟩ := area_rot an2
  unfold diamondEnergy rippaForm cotEnergy cornerCot
  generalize flipLength a b c d e = f at *
  rw [k1, k2, k3, l1, l2, l3, m1, m2, m3, n1, n2, n3, eb, ec, ed, ee, ef]
  have hH' := hH.ne'
  have hG' := hG.ne'
  have ha' := ha.ne'
  generalize hP : X * G + Y * H = P at *
  generalize hQ : (a - X) * G + (a - Y) * H = Q at *
  have h1' := h1.ne'
  have h2' := h2.ne'
  field_simp
  subst hP hQ
  ring

/-- The form vanishes on constants (it vanishes on every affine function of the layout). -/
theorem rippaForm_const (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) (t : ℝ) : rippaForm a b c d e t t t t = 0 := by
  have h := flip_area hk hl hconv
  rw [triArea_swap (flipLength a b c d e) c e] at h
  unfold rippaForm
  linear_combination (-t) * h

/-- The factor of Rippa's identity has the sign of `θ_k + θ_l − π`. -/
private theorem rippa_factor_sign (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    0 < 2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e := by
  obtain ⟨t1, t2⟩ := flipLength_isTriangle hk hl hconv
  have := triArea_pos t1
  have := triArea_pos t2.swap
  positivity

/-- **Rippa's theorem for one diamond**: if `θ_k + θ_l ≥ π` (the old edge is not locally
Delaunay, or a tie), the flip does not increase the Dirichlet energy, for every `u`.

[formula: it-rippa-diamond] -/
theorem rippa_le (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : π ≤ cornerAngle a b c + cornerAngle a d e) (ui uj uk ul : ℝ) :
    diamondEnergy (flipLength a b c d e) c e b d uk ul uj ui ≤
      diamondEnergy a b c d e ui uj uk ul := by
  have hconv := isConvexDiamond_of_pi_le hk hl h
  have hid := rippa_identity hk hl hconv ui uj uk ul
  have hD := rippa_factor_sign hk hl hconv
  have hc : cornerCot a b c + cornerCot a d e ≤ 0 := by
    rcases h.lt_or_eq with h' | h'
    · exact ((pi_lt_iff_cot_neg hk hl).mp h').le
    · exact ((eq_pi_iff_cot_eq_zero hk hl).mp h'.symm).le
  have : 0 ≤ -(cornerCot a b c + cornerCot a d e) /
      (2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e) *
        rippaForm a b c d e ui uj uk ul ^ 2 :=
    mul_nonneg (div_nonneg (by linarith) hD.le) (sq_nonneg _)
  linarith

/-- **Strict decrease**: if the old edge is not locally Delaunay, the flip lowers the energy of
every `u` with `ℓ(u) ≠ 0` (for instance every hat function of `i, j, k, l`). -/
theorem rippa_lt (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : ¬ IsLocallyDelaunay a b c d e) {ui uj uk ul : ℝ}
    (hu : rippaForm a b c d e ui uj uk ul ≠ 0) :
    diamondEnergy (flipLength a b c d e) c e b d uk ul uj ui <
      diamondEnergy a b c d e ui uj uk ul := by
  have hconv := isConvexDiamond_of_not_isLocallyDelaunay hk hl h
  have hid := rippa_identity hk hl hconv ui uj uk ul
  have hD := rippa_factor_sign hk hl hconv
  have hc := (pi_lt_iff_cot_neg hk hl).mp (not_isLocallyDelaunay_iff.mp h)
  have : 0 < -(cornerCot a b c + cornerCot a d e) /
      (2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e) *
        rippaForm a b c d e ui uj uk ul ^ 2 :=
    mul_pos (div_pos (by linarith) hD) (by positivity)
  linarith

/-- **At a tie the two triangulations of the diamond have the same energy** for every `u`
(B–S: the cotan Laplacian does not depend on how the Delaunay tessellation is triangulated). -/
theorem rippa_eq_of_tie (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (h : cornerAngle a b c + cornerAngle a d e = π) (ui uj uk ul : ℝ) :
    diamondEnergy (flipLength a b c d e) c e b d uk ul uj ui =
      diamondEnergy a b c d e ui uj uk ul := by
  have hconv := isConvexDiamond_of_pi_le hk hl h.ge
  have hid := rippa_identity hk hl hconv ui uj uk ul
  rw [(eq_pi_iff_cot_eq_zero hk hl).mp h, neg_zero, zero_div, zero_mul] at hid
  linarith

/-- **The converse**: for a strictly convex diamond, the flip lowers or keeps the energy of every
`u` exactly when `θ_k + θ_l ≥ π`. (If the edge is strictly Delaunay, the hat function of `k`
has a larger energy after the flip.)

[formula: it-rippa-diamond] -/
theorem forall_rippa_le_iff (hk : IsTriangle a b c) (hl : IsTriangle a d e)
    (hconv : IsConvexDiamond a b c d e) :
    (∀ ui uj uk ul : ℝ, diamondEnergy (flipLength a b c d e) c e b d uk ul uj ui ≤
        diamondEnergy a b c d e ui uj uk ul) ↔
      π ≤ cornerAngle a b c + cornerAngle a d e := by
  refine ⟨fun hall => ?_, fun h => rippa_le hk hl h⟩
  by_contra hlt
  push Not at hlt
  have hc : 0 < cornerCot a b c + cornerCot a d e := by
    by_contra hc
    push Not at hc
    rcases hc.lt_or_eq with hc | hc
    · have := (pi_lt_iff_cot_neg hk hl).mpr hc; linarith
    · have := (eq_pi_iff_cot_eq_zero hk hl).mpr hc; linarith
  have hid := rippa_identity hk hl hconv 0 0 1 0
  have hD := rippa_factor_sign hk hl hconv
  have hT := triArea_pos hl
  have hform : rippaForm a b c d e 0 0 1 0 = triArea a d e := by unfold rippaForm; ring
  rw [hform] at hid
  have : -(cornerCot a b c + cornerCot a d e) /
      (2 * triArea (flipLength a b c d e) b d * triArea (flipLength a b c d e) c e) *
        triArea a d e ^ 2 < 0 :=
    mul_neg_of_neg_of_pos (div_neg_of_neg_of_pos (by linarith) hD) (by positivity)
  have := hall 0 0 1 0
  linarith

end Delaunay

end DDG
