import DdgProofs.Developability.Basic

/-!
# Developability: the max covariance energy `λ_max` (Eq. 6, App. B.4)

`crates/ddg-developability/SPEC.md` §3.5, decisions D5, D11. The energy is
`λ_max = min_{|u|=1} max_σ ⟨u, N_σ⟩²`; its candidates are spherical circumcenters of triples of
sites `±N_σ`. That the minimum is attained at a candidate (the spherical Voronoi argument of the
SPEC) is cited (`dev-max-covariance-voronoi`, `CITED.md`). Proved here, about the oracle's
`Ddg.Dev.siteNormal`, `circumcenter`, `circumG`, `maxCovSensitivities`:

* the circumcenter formula: `v = n/|n|`, `n = (b − a) × (c − a)`, is equidistant from the three
  sites, `⟨v,a⟩ = ⟨v,b⟩ = ⟨v,c⟩ = det(a,b,c)/|n|`, and for unit sites the only unit vectors
  equidistant from them are `±v` (`dev-spherical-circumcenter`);
* the sensitivities of `e = g²`, `g = ⟨v, a⟩`, along arbitrary motions of the three sites
  (`dev-max-covariance-gradient`); the site sign of D5 is the chain rule for `q = s_q N_σ`
  (`hasDerivAt_circumG_sq_signed`).
-/

noncomputable section

namespace DDG

namespace Dev

open scoped DDG.Oracle
open Real Ddg.Dev

theorem siteNormal_eq_faceNormal (a b c : V3) : siteNormal a b c = faceNormal a b c := rfl

theorem circumcenter_eq (a b c : V3) : circumcenter a b c = faceUnitNormal a b c := rfl

/-- **Spherical circumcenter** (SPEC §3.5, D5): for sites `a, b, c` with
`n = (b − a) × (c − a) ≠ 0` and `v = n/|n|`,
`⟨v, a⟩ = ⟨v, b⟩ = ⟨v, c⟩ = det(a, b, c)/|n| = g`, and `|n| = 2 A_abc` (twice the Euclidean area of
the triangle `abc`). Hence `v` (and `−v`) is equidistant on the sphere from three unit sites; the
paper's "spherical centroid" is this circumcenter.

[formula: dev-spherical-circumcenter] -/
theorem circumcenter_equidistant (a b c : V3) (h : siteNormal a b c ≠ zero) :
    dot (circumcenter a b c) a = circumG a b c ∧ dot (circumcenter a b c) b = circumG a b c ∧
    dot (circumcenter a b c) c = circumG a b c ∧ norm2 (circumcenter a b c) = 1 ∧
    norm (siteNormal a b c) = 2 * faceArea a b c := by
  have hp := (norm_pos_iff' _).mpr h
  have key : ∀ q : V3, dot (siteNormal a b c) q = det3 a b c →
      dot (circumcenter a b c) q = circumG a b c := by
    intro q hq
    unfold circumcenter circumG
    rw [normalize_def, dot_smul_l, hq]
    ring
  refine ⟨key a ?_, key b ?_, key c ?_, norm2_normalize _ h, ?_⟩
  · simp only [siteNormal, det3_def, dot_def, cross_x, cross_y, cross_z, sub_x, sub_y, sub_z]; ring
  · simp only [siteNormal, det3_def, dot_def, cross_x, cross_y, cross_z, sub_x, sub_y, sub_z]; ring
  · simp only [siteNormal, det3_def, dot_def, cross_x, cross_y, cross_z, sub_x, sub_y, sub_z]; ring
  · unfold faceArea; rw [siteNormal_eq_faceNormal]; simp only [Oracle.two_eq]; ring

/-- If `u ⊥ p` and `u ⊥ q`, then `|p × q|² u = ((p × q)·u) (p × q)`. -/
theorem smul_par_of_perp (p q u : V3) (hp : dot u p = 0) (hq : dot u q = 0) :
    smul (norm2 (cross p q)) u = smul (dot (cross p q) u) (cross p q) := by
  rw [dot_def] at hp hq
  ext <;> simp only [norm2_def, dot_def, smul_x, smul_y, smul_z, cross_x, cross_y, cross_z]
  · linear_combination ((p.z * q.x - p.x * q.z) * p.z - (p.x * q.y - p.y * q.x) * p.y) * hq - ((p.z * q.x - p.x * q.z) * q.z - (p.x * q.y - p.y * q.x) * q.y) * hp
  · linear_combination ((p.x * q.y - p.y * q.x) * p.x - (p.y * q.z - p.z * q.y) * p.z) * hq - ((p.x * q.y - p.y * q.x) * q.x - (p.y * q.z - p.z * q.y) * q.z) * hp
  · linear_combination ((p.y * q.z - p.z * q.y) * p.y - (p.z * q.x - p.x * q.z) * p.x) * hq - ((p.y * q.z - p.z * q.y) * q.y - (p.z * q.x - p.x * q.z) * q.x) * hp

/-- **The circumcenter is the only equidistant point** (up to sign): a unit vector `u` with
`⟨u, a⟩ = ⟨u, b⟩ = ⟨u, c⟩` is `v` or `−v` (for `n ≠ 0`). Together with
`circumcenter_equidistant`: the Voronoi vertices of the sites `{±N_σ}` are among the circumcenters
`±v` of triples of sites (the candidates of SPEC §3.5).

[formula: dev-spherical-circumcenter] -/
theorem eq_circumcenter_of_equidistant (a b c u : V3) (h : siteNormal a b c ≠ zero)
    (hu : norm2 u = 1) (hab : dot u a = dot u b) (hac : dot u a = dot u c) :
    u = circumcenter a b c ∨ u = neg (circumcenter a b c) := by
  have hup : dot u (sub b a) = 0 := by rw [dot_sub_r]; linarith
  have huq : dot u (sub c a) = 0 := by rw [dot_sub_r]; linarith
  have hpar := smul_par_of_perp _ _ u hup huq
  rw [← siteNormal] at hpar
  unfold circumcenter
  generalize siteNormal a b c = n at h hpar ⊢
  have hp := (norm_pos_iff' n).mpr h
  have hP := norm_sq n
  -- (n·u)² = |n|²
  have hsq : dot n u ^ 2 = norm n ^ 2 := by
    have := congrArg (fun w => dot w u) hpar
    simp only [dot_smul_l] at this
    rw [show dot u u = norm2 u from rfl, hu] at this
    rw [hP, sq]; linarith [this]
  have hP0 : norm2 n ≠ 0 := by rw [← hP]; positivity
  have hu' : u = smul (dot n u / norm2 n) n := by
    ext
    · have := congrArg Ddg.Point.x hpar; simp only [smul_x] at this ⊢; field_simp; linarith
    · have := congrArg Ddg.Point.y hpar; simp only [smul_y] at this ⊢; field_simp; linarith
    · have := congrArg Ddg.Point.z hpar; simp only [smul_z] at this ⊢; field_simp; linarith
  have hr := hp.ne'
  rcases sq_eq_sq_iff_eq_or_eq_neg.mp hsq with h1 | h1
  · left
    rw [hu', h1, ← hP]
    rw [normalize_def]
    ext <;> simp only [smul_x, smul_y, smul_z] <;> field_simp
  · right
    rw [hu', h1, ← hP]
    rw [normalize_def]
    ext <;> simp only [smul_x, smul_y, smul_z, neg_x, neg_y, neg_z] <;> field_simp

/-! ## Sensitivities (SPEC §3.5, D5) -/

section gradient

variable {a b c : ℝ → V3} {da db dc : V3} {t : ℝ}

/-- **Gradient of `λ_max` at a generic minimizer** (SPEC §3.5, corrected App. B.4, D5): with
`v = n/|n|` the circumcenter of the active sites `a, b, c`, `g = ⟨v, a⟩` (`= det(a,b,c)/|n|`) and
`e_v = g²`, along any motion of the three sites
`d(g²) = Σ_q (∂e_v/∂q) · dq`, `∂e_v/∂q = 2g [[q = a] v + ((e_q × v)·a/|n|) v]`,
`e_a = c − b`, `e_b = a − c`, `e_c = b − a` (`maxCovSensitivities`). No unit length of the sites is
needed.

[formula: dev-max-covariance-gradient] -/
theorem hasDerivAt_circumG_sq (ha : HasDerivAt3 a da t) (hb : HasDerivAt3 b db t)
    (hc : HasDerivAt3 c dc t) (hn : siteNormal (a t) (b t) (c t) ≠ zero) :
    HasDerivAt (fun s => dot (circumcenter (a s) (b s) (c s)) (a s) ^ 2)
      (dot (maxCovSensitivities (a t) (b t) (c t)).1 da +
        dot (maxCovSensitivities (a t) (b t) (c t)).2.1 db +
        dot (maxCovSensitivities (a t) (b t) (c t)).2.2 dc) t := by
  have dv := hasDerivAt3_faceUnitNormal ha hb hc hn
  have dg := dv.dot ha
  have h2 := dg.pow 2
  refine h2.congr_deriv ?_
  unfold maxCovSensitivities circumSensitivity normalDerivative
  rw [faceUnitNormal_cycle (a t) (b t) (c t), faceUnitNormal_cycle (b t) (c t) (a t),
    faceUnitNormal_cycle (a t) (b t) (c t), faceArea_cycle (a t) (b t) (c t),
    faceArea_cycle (b t) (c t) (a t), faceArea_cycle (a t) (b t) (c t)]
  have hA : 2 * faceArea (a t) (b t) (c t) = norm (siteNormal (a t) (b t) (c t)) := by
    unfold faceArea; simp only [Oracle.two_eq]; rw [siteNormal_eq_faceNormal]; ring
  rw [hA, ← circumcenter_eq]
  simp only [Bool.false_eq_true, ↓reduceIte]
  generalize circumcenter (a t) (b t) (c t) = v
  generalize norm (siteNormal (a t) (b t) (c t)) = r
  simp only [dot_add_l, dot_sub_l, dot_smul_l, dot_cross_l, Oracle.two_eq, Nat.cast_ofNat]
  simp only [dot_def, zero_x, zero_y, zero_z]
  ring

/-- **The site sign** (D5): a site `q = s_q N_σ` (`s_q = ±1`, the antipode for `s_q = −1`) passes
`s_q ∂e_v/∂q` to the normal `N_σ`. Stated for three sites `s_a N_a, s_b N_b, s_c N_c` moving with
their normals; for candidates of the second kind (`(N_1, −N_1, N_2)`) the motions of `a` and `b`
are those of one face and the two contributions add up.

[formula: dev-max-covariance-gradient] -/
theorem hasDerivAt_circumG_sq_signed {Na Nb Nc : ℝ → V3} {dNa dNb dNc : V3} (sa sb sc : ℝ)
    (ha : HasDerivAt3 Na dNa t) (hb : HasDerivAt3 Nb dNb t) (hc : HasDerivAt3 Nc dNc t)
    (hn : siteNormal (smul sa (Na t)) (smul sb (Nb t)) (smul sc (Nc t)) ≠ zero) :
    HasDerivAt (fun s => dot (circumcenter (smul sa (Na s)) (smul sb (Nb s)) (smul sc (Nc s)))
        (smul sa (Na s)) ^ 2)
      (dot (smul sa (maxCovSensitivities (smul sa (Na t)) (smul sb (Nb t)) (smul sc (Nc t))).1)
          dNa +
        dot (smul sb (maxCovSensitivities (smul sa (Na t)) (smul sb (Nb t))
          (smul sc (Nc t))).2.1) dNb +
        dot (smul sc (maxCovSensitivities (smul sa (Na t)) (smul sb (Nb t))
          (smul sc (Nc t))).2.2) dNc) t := by
  have h := hasDerivAt_circumG_sq (ha.smul_const sa) (hb.smul_const sb) (hc.smul_const sc) hn
  refine h.congr_deriv ?_
  simp only [dot_smul_l, dot_smul_r]

end gradient

end Dev

end DDG
