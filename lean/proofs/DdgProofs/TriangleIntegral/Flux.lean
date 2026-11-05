import DdgProofs.TriangleIntegral.Fem

/-!
# The flux of a linear function through a side (`app-amr-edge-flux`)

The jump residuals of `amr::estimate` (`crates/ddg-intrinsic-apps/src/amr.rs`) are sums of the
fluxes `side_flux`: for the side `k` of a face (from the corner `k` to the corner `k + 1`, the
opposite corner `o = k + 2`), `cot θ_{k+1} (u_k − u_o) + cot θ_k (u_{k+1} − u_o)`, the cotangents
from the face lengths `l_c = |p_c p_{c+1}|` as `corner_cotan(l_{c+1}, l_c, l_{c+2})`.

* `sideNormal p k` — `ℓ ν`, the outward normal of the side `k` scaled by its length: orthogonal to
  the side (`dot_sideNormal_side`), of length `ℓ` (`dot_sideNormal_self`), pointing away from the
  opposite corner (`dot_sideNormal_opposite`), for either orientation of the triangle;
* `dot_grad_sideNormal` — **the flux** `ℓ ∇u·ν = cot θ_{k+1} (u_k − u_o) + cot θ_k (u_{k+1} − u_o)`
  for the linear interpolant `u` of the corner values (the gradient from `fderiv`), and
  `side_flux_eq` the same with the cotangents in the argument order of `side_flux`;
* `sideNormal_twin`, `side_flux_add_twin` — the two faces of an interior side, consistently
  oriented, have opposite normals, so the sum of their fluxes (`flux[h] + flux[t]` in `estimate`)
  is `ℓ (∇u_T − ∇u_{T'})·ν_T`, the length times the jump of the normal derivative.
-/

noncomputable section

namespace DdgProofs

namespace TriangleIntegral

/-- The outward normal of the side `k` (from `p_k` to `p_{k+1}`) times its length:
`sign(D) · (e₂, −e₁)` with `e = p_{k+1} − p_k`. -/
def sideNormal (p : Fin 3 → ℝ × ℝ) (k : Fin 3) : ℝ × ℝ :=
  (|D p| / D p) • ((p (k + 1) - p k).2, -(p (k + 1) - p k).1)

/-- The face lengths in the order of the Rust code: `l_c = |p_c p_{c+1}|`, the side of the
halfedge `3f + c`. -/
def faceLen (p : Fin 3 → ℝ × ℝ) (c : Fin 3) : ℝ := dist2 (p c) (p (c + 1))

theorem sign_sq {d : ℝ} (hd : d ≠ 0) : (|d| / d) ^ 2 = 1 := by
  rw [div_pow, sq_abs, div_self (pow_ne_zero 2 hd)]

theorem dot_sideNormal_side (p : Fin 3 → ℝ × ℝ) (k : Fin 3) :
    dot2 (sideNormal p k) (p (k + 1) - p k) = 0 := by
  simp only [sideNormal, dot2, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  ring

theorem dot_sideNormal_self {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (k : Fin 3) :
    dot2 (sideNormal p k) (sideNormal p k) = sqDist2 (p k) (p (k + 1)) := by
  have h := sign_sq hD
  simp only [sideNormal, dot2, sqDist2, Prod.smul_fst, Prod.smul_snd, smul_eq_mul, Prod.fst_sub,
    Prod.snd_sub]
  linear_combination ((p (k + 1)).1 - (p k).1) ^ 2 * h + ((p (k + 1)).2 - (p k).2) ^ 2 * h

/-- The normal points away from the opposite corner: `ℓν · (p_{k+2} − p_k) = −|D| < 0`. -/
theorem dot_sideNormal_opposite {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (k : Fin 3) :
    dot2 (sideNormal p k) (p (k + 2) - p k) = -|D p| := by
  have e : ∀ k : Fin 3, sda (p k) (p (k + 1)) (p (k + 2)) = D p := by
    intro k
    fin_cases k
    · rfl
    · exact sda_rotate _ _ _
    · simp only [Fin.reduceFinMk, Fin.isValue, Fin.reduceAdd, D]
      rw [← sda_rotate, sda_rotate]
  have h : dot2 (sideNormal p k) (p (k + 2) - p k) =
      -(|D p| / D p) * sda (p k) (p (k + 1)) (p (k + 2)) := by
    simp only [sideNormal, dot2, sda, Prod.smul_fst, Prod.smul_snd, smul_eq_mul, Prod.fst_sub,
      Prod.snd_sub]
    ring
  rw [h, e k]
  field_simp

/-- **The flux of a linear function through a side**: for the linear interpolant `u` of the
corner values, `ℓ ∇u · ν = cot θ_{k+1} (u_k − u_{k+2}) + cot θ_k (u_{k+1} − u_{k+2})`, the
cotangents of the corners at the ends of the side, from the side lengths.

[formula: app-amr-edge-flux] -/
theorem dot_grad_sideNormal {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (u : Fin 3 → ℝ) (k : Fin 3)
    (x : ℝ × ℝ) :
    dot2 (grad (interp p u) x) (sideNormal p k) =
      cotCorner p (k + 1) * (u k - u (k + 2)) + cotCorner p k * (u (k + 1) - u (k + 2)) := by
  rw [grad_interp]
  have c0 := cotCorner_zero p
  have c1 := cotCorner_one p
  have c2 := cotCorner_two p
  fin_cases k <;>
    simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk, Fin.isValue, Fin.reduceAdd, zero_add,
      c0, c1, c2, gradInterp, baryGrad, sideNormal, side_zero, side_one, side_two, dot2,
      Fin.sum_univ_three, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul,
      Prod.fst_sub, Prod.snd_sub] <;>
    · set d := D p with hd
      have hD' : -d ≠ 0 := neg_ne_zero.mpr hD
      rcases abs_choice d with h | h <;> rw [h] <;> field_simp <;> ring

theorem cornerCot_swap (a b c : ℝ) : cornerCot a c b = cornerCot a b c := by
  simp only [cornerCot, triArea_swap]
  ring

/-- The cotangent as `side_flux` evaluates it, `corner_cotan(l_{c+1}, l_c, l_{c+2})`, is the
cotangent of the corner `c`. -/
theorem cornerCot_faceLen (p : Fin 3 → ℝ × ℝ) (c : Fin 3) :
    cornerCot (faceLen p (c + 1)) (faceLen p c) (faceLen p (c + 2)) = cotCorner p c := by
  rw [cotCorner, ← cornerCot_swap]
  fin_cases c <;> rfl

/-- **`side_flux` is the flux**: with the cotangents from the face lengths in the argument order
of `crates/ddg-intrinsic-apps/src/amr.rs`,
`cot_{k+1} (u_k − u_o) + cot_k (u_{k+1} − u_o) = ℓ ∇u · ν` for the side `k` (from the corner `k`
to `k + 1`, `o = k + 2`) and the outward normal `ν`.

[formula: app-amr-edge-flux] -/
theorem side_flux_eq {p : Fin 3 → ℝ × ℝ} (hD : D p ≠ 0) (u : Fin 3 → ℝ) (k : Fin 3) (x : ℝ × ℝ) :
    cornerCot (faceLen p (k + 1 + 1)) (faceLen p (k + 1)) (faceLen p (k + 1 + 2)) *
        (u k - u (k + 2)) +
      cornerCot (faceLen p (k + 1)) (faceLen p k) (faceLen p (k + 2)) * (u (k + 1) - u (k + 2)) =
      dot2 (grad (interp p u) x) (sideNormal p k) := by
  rw [cornerCot_faceLen, cornerCot_faceLen, dot_grad_sideNormal hD]

/-- **Twin sides have opposite normals**: if the side `k'` of `p'` is the side `k` of `p`
traversed backwards and the two triangles have the same orientation, `ℓν' = −ℓν`. -/
theorem sideNormal_twin {p p' : Fin 3 → ℝ × ℝ} {k k' : Fin 3} (h0 : p' k' = p (k + 1))
    (h1 : p' (k' + 1) = p k) (hs : 0 < D p * D p') : sideNormal p' k' = -sideNormal p k := by
  have hsign : |D p'| / D p' = |D p| / D p := by
    rcases lt_or_gt_of_ne (show D p ≠ 0 by rintro h; simp [h] at hs) with h | h
    · have h' : D p' < 0 := by nlinarith
      rw [abs_of_neg h, abs_of_neg h', neg_div_self h.ne, neg_div_self h'.ne]
    · have h' : 0 < D p' := by nlinarith
      rw [abs_of_pos h, abs_of_pos h', div_self h.ne', div_self h'.ne']
  simp only [sideNormal, h0, h1, hsign]
  ext <;> simp only [Prod.smul_fst, Prod.smul_snd, Prod.fst_neg, Prod.snd_neg, Prod.fst_sub,
    Prod.snd_sub, smul_eq_mul] <;> ring

/-- **The jump through an interior side**: the fluxes of the two faces add up to
`ℓ (∇u_T − ∇u_{T'}) · ν_T` (`flux[h] + flux[t]` in `amr::estimate`).

[formula: app-amr-edge-flux] -/
theorem side_flux_add_twin {p p' : Fin 3 → ℝ × ℝ} {k k' : Fin 3} (h0 : p' k' = p (k + 1))
    (h1 : p' (k' + 1) = p k) (hs : 0 < D p * D p') (u u' : Fin 3 → ℝ) (x : ℝ × ℝ) :
    dot2 (grad (interp p u) x) (sideNormal p k) + dot2 (grad (interp p' u') x) (sideNormal p' k') =
      dot2 (grad (interp p u) x - grad (interp p' u') x) (sideNormal p k) := by
  rw [sideNormal_twin h0 h1 hs]
  simp only [dot2, Prod.fst_sub, Prod.snd_sub, Prod.fst_neg, Prod.snd_neg]
  ring

end TriangleIntegral

end DdgProofs
