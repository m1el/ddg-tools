import DdgProofs.Developability.Vec

/-!
# Developability: the per-face quantities and their derivatives

`crates/ddg-developability/SPEC.md` §1–§2; paper ("Developability of Triangle Meshes", Stein,
Grinspun, Crane 2018) Sec. 3.1, App. B.1 (Eqs. 7–9), App. B.5 (vertex normal).

Statements are about the oracle's definitions (`Ddg.Dev.faceNormal`, `faceArea`,
`faceUnitNormal`, `cornerAngleAt`, `areaGradient`, `normalDerivative`, `angleGradientI/J/K`,
`backpropW`, `backpropFace`, `vectorArea`, `ringVectorArea`) at `ℝ`.

Derivatives are taken **along arbitrary differentiable motions** `t ↦ (f_i(t), f_j(t), f_k(t))` of
the three corners (`HasDerivAt3`, coordinatewise), with velocities `δ_i, δ_j, δ_k`. The derivative
is then `Σ_p g_p · δ_p`, i.e. `g_p` is the gradient with respect to `f_p` (moving one corner only
is the case of the other two velocities `0`; the directional derivative is the case of lines).
This is the chain-rule form in which the Rust code uses the gradients.

The corrections of the SPEC are part of the statements: the normal derivative has the factor
`1/(2A)` (D2, the paper prints `1/A`), the angle gradient divides by the *squared* edge length
(D3, the paper prints `|f_i − f_j|`).
-/

noncomputable section

namespace DDG

namespace Dev

open scoped DDG.Oracle
open Real Ddg.Dev

/-! ## Helpers -/

section helpers

variable (a b c : V3)

theorem faceNormal_def : faceNormal a b c = cross (sub b a) (sub c a) := rfl

theorem faceUnitNormal_def : faceUnitNormal a b c = normalize (faceNormal a b c) := rfl

/-- The face normal is invariant under a cyclic relabeling of the corners. -/
theorem faceNormal_cycle : faceNormal b c a = faceNormal a b c := by
  ext <;> simp [faceNormal_def] <;> ring

theorem faceUnitNormal_cycle : faceUnitNormal b c a = faceUnitNormal a b c := by
  rw [faceUnitNormal_def, faceUnitNormal_def, faceNormal_cycle]

theorem faceArea_cycle : faceArea b c a = faceArea a b c := by
  unfold faceArea; rw [faceNormal_cycle]

/-- `n · (f_j − f_i) = 0`, `n · (f_k − f_i) = 0`, `n · (f_k − f_j) = 0`. -/
theorem dot_faceNormal_edges :
    dot (faceNormal a b c) (sub b a) = 0 ∧ dot (faceNormal a b c) (sub c a) = 0 ∧
      dot (faceNormal a b c) (sub c b) = 0 ∧ dot (faceNormal a b c) (sub a c) = 0 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> simp only [faceNormal_def, dot_def, cross_x, cross_y, cross_z,
    sub_x, sub_y, sub_z] <;> ring

/-- `y = |n|` written as `N · n`: used to keep a single factor `1/|n|` in rational identities. -/
theorem norm_eq_dot_normalize (h : a ≠ zero) : norm a = dot (normalize a) a := by
  have hp := (norm_pos_iff' a).mpr h
  have h2 := norm_mul_self a
  rw [normalize_def]
  have : dot (smul (1 / norm a) a) a = norm2 a / norm a := by
    simp only [dot_def, norm2_def, smul_x, smul_y, smul_z]; field_simp
  rw [this, ← h2]; field_simp

/-- A nondegenerate face has distinct corners. -/
theorem ne_of_faceNormal_ne (h : faceNormal a b c ≠ zero) :
    sub b a ≠ zero ∧ sub c a ≠ zero ∧ sub c b ≠ zero := by
  have key : ∀ p q : V3, sub p q = zero → p = q := by
    intro p q h0
    have hx := congrArg Ddg.Point.x h0
    have hy := congrArg Ddg.Point.y h0
    have hz := congrArg Ddg.Point.z h0
    simp at hx hy hz
    ext <;> linarith
  refine ⟨fun h0 => h ?_, fun h0 => h ?_, fun h0 => h ?_⟩ <;> rw [key _ _ h0] <;>
    ext <;> simp [faceNormal_def] <;> ring

/-- `N · (e × δ) = (N × e) · δ` and the projector identity behind Eq. 8:
for a unit `N` and `e ⊥ N`, `(e × δ) − (N · (e × δ)) N = (N · δ) (e × N)`. -/
theorem proj_cross {N e δ : V3} (hN : norm2 N = 1) (he : dot N e = 0) :
    sub (cross e δ) (smul (dot N (cross e δ)) N) = smul (dot N δ) (cross e N) := by
  rw [norm2_def] at hN
  rw [dot_def] at he
  ext <;> simp only [sub_x, sub_y, sub_z, smul_x, smul_y, smul_z, cross_x, cross_y, cross_z,
    dot_def]
  · linear_combination (-(e.y * δ.z - e.z * δ.y)) * hN + (N.y * δ.z - N.z * δ.y) * he
  · linear_combination (-(e.z * δ.x - e.x * δ.z)) * hN + (N.z * δ.x - N.x * δ.z) * he
  · linear_combination (-(e.x * δ.y - e.y * δ.x)) * hN + (N.x * δ.y - N.y * δ.x) * he

/-! Linearity of `dot` and `cross`, for rewriting vector expressions into sums of atoms. -/

theorem dot_add_l (u v w : V3) : dot (add u v) w = dot u w + dot v w := by
  simp only [dot_def, add_x, add_y, add_z]; ring
theorem dot_add_r (u v w : V3) : dot w (add u v) = dot w u + dot w v := by
  simp only [dot_def, add_x, add_y, add_z]; ring
theorem dot_sub_l (u v w : V3) : dot (sub u v) w = dot u w - dot v w := by
  simp only [dot_def, sub_x, sub_y, sub_z]; ring
theorem dot_sub_r (u v w : V3) : dot w (sub u v) = dot w u - dot w v := by
  simp only [dot_def, sub_x, sub_y, sub_z]; ring
theorem dot_smul_l (c : ℝ) (u w : V3) : dot (smul c u) w = c * dot u w := by
  simp only [dot_def, smul_x, smul_y, smul_z]; ring
theorem dot_smul_r (c : ℝ) (u w : V3) : dot w (smul c u) = c * dot w u := by
  simp only [dot_def, smul_x, smul_y, smul_z]; ring
theorem dot_neg_l (u w : V3) : dot (neg u) w = -dot u w := by
  simp only [dot_def, neg_x, neg_y, neg_z]; ring
theorem dot_neg_r (u w : V3) : dot w (neg u) = -dot w u := by
  simp only [dot_def, neg_x, neg_y, neg_z]; ring
/-- `(w × e) · δ = w · (e × δ)`. -/
theorem dot_cross_l (w e δ : V3) : dot (cross w e) δ = dot w (cross e δ) := by
  simp only [dot_def, cross_x, cross_y, cross_z]; ring
theorem cross_smul_l (c : ℝ) (u v : V3) : cross (smul c u) v = smul c (cross u v) := by
  ext <;> simp <;> ring
theorem cross_add_r (u v w : V3) : cross w (add u v) = add (cross w u) (cross w v) := by
  ext <;> simp <;> ring
theorem cross_add_l (u v w : V3) : cross (add u v) w = add (cross u w) (cross v w) := by
  ext <;> simp <;> ring

/-- The polynomial identity behind Eq. 9 (corrected, D3): with `n = a × b`, `x = a·b`,
`dn = da × b + a × db`, `x' = da·b + a·db`:
`x (n·dn) − |n|² x' = |b|² (n × (−a))·da + |a|² (n × b)·db`. -/
theorem angle_poly (a b da db : V3) :
    dot a b * dot (cross a b) (add (cross da b) (cross a db)) -
        norm2 (cross a b) * (dot da b + dot a db) =
      norm2 b * dot (cross (cross a b) (neg a)) da + norm2 a * dot (cross (cross a b) b) db := by
  simp only [dot_def, norm2_def, add_x, add_y, add_z, cross_x, cross_y, cross_z, neg_x, neg_y,
    neg_z]
  ring

end helpers

/-! ## Face frame -/

/-- **Face frame** (SPEC §1): `n = (f_j − f_i) × (f_k − f_i) = f_i × f_j + f_j × f_k + f_k × f_i`
(cyclic form), `N = n/|n|` is a unit vector orthogonal to the three edges, and the area
`A = |n|/2` equals Heron's formula `triArea` of the three edge lengths (the area the rest of the
workspace computes from lengths).

[formula: dev-face-frame] -/
theorem face_frame (fi fj fk : V3) (h : faceNormal fi fj fk ≠ zero) :
    faceNormal fi fj fk = add (add (cross fi fj) (cross fj fk)) (cross fk fi) ∧
    norm2 (faceUnitNormal fi fj fk) = 1 ∧
    dot (faceUnitNormal fi fj fk) (sub fj fi) = 0 ∧
    dot (faceUnitNormal fi fj fk) (sub fk fi) = 0 ∧
    dot (faceUnitNormal fi fj fk) (sub fk fj) = 0 ∧
    faceArea fi fj fk = Ddg.triArea (Ddg.dist fj fk) (Ddg.dist fk fi) (Ddg.dist fi fj) := by
  have hN := norm2_normalize _ h
  have hd := dot_faceNormal_edges fi fj fk
  have hdot : ∀ e : V3, dot (faceNormal fi fj fk) e = 0 → dot (faceUnitNormal fi fj fk) e = 0 := by
    intro e he
    rw [faceUnitNormal_def, normalize_def]
    rw [dot_def] at he ⊢
    simp only [smul_x, smul_y, smul_z]
    linear_combination (1 / norm (faceNormal fi fj fk)) * he
  refine ⟨?_, hN, hdot _ hd.1, hdot _ hd.2.1, hdot _ hd.2.2.1, ?_⟩
  · ext <;> simp [faceNormal_def] <;> ring
  · -- Heron
    set a := Ddg.dist fj fk
    set b := Ddg.dist fk fi
    set c := Ddg.dist fi fj
    have sq_dist : ∀ p q : V3, Ddg.dist p q ^ 2 =
        (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z) := by
      intro p q
      show (√_) ^ 2 = _
      rw [sq_sqrt]
      nlinarith [mul_self_nonneg (p.x - q.x), mul_self_nonneg (p.y - q.y),
        mul_self_nonneg (p.z - q.z)]
    have ha := sq_dist fj fk
    have hb := sq_dist fk fi
    have hc := sq_dist fi fj
    have hprod : (a + b + c) * (b - a + c) * (a - b + c) * (a - c + b) =
        4 * norm2 (faceNormal fi fj fk) := by
      have e1 : (a + b + c) * (b - a + c) * (a - b + c) * (a - c + b) =
          2 * a ^ 2 * b ^ 2 + 2 * b ^ 2 * c ^ 2 + 2 * c ^ 2 * a ^ 2 - (a ^ 2) ^ 2 - (b ^ 2) ^ 2 -
            (c ^ 2) ^ 2 := by ring
      rw [e1, ha, hb, hc]
      simp only [norm2_def, faceNormal_def, cross_x, cross_y, cross_z, sub_x, sub_y, sub_z]
      ring
    unfold Ddg.triArea faceArea
    simp only [Oracle.sqrt_eq, Oracle.four_eq, Oracle.two_eq]
    rw [hprod, sqrt_mul (by norm_num : (0 : ℝ) ≤ 4), show √(4 : ℝ) = 2 by
      rw [show (4 : ℝ) = 2 ^ 2 by norm_num, sqrt_sq (by norm_num)], norm_def]
    ring

/-! ## Corner angle -/

/-- **Corner angle** (SPEC §1, D13): for `a = f_j − f_i ≠ 0`, `b = f_k − f_i ≠ 0`, the angle
`θ = atan2(|a × b|, a · b)` lies in `[0, π]` and has `cos θ = a·b/(|a||b|)`,
`sin θ = |a × b|/(|a||b|)`; hence it is the Euclidean angle `arccos(a·b/(|a||b|))` between `a`
and `b` (`cornerAngleAt_eq_arccos`).

[formula: dev-corner-angle] -/
theorem cornerAngleAt_spec (fi fj fk : V3) (ha : sub fj fi ≠ zero) (hb : sub fk fi ≠ zero) :
    0 ≤ cornerAngleAt fi fj fk ∧ cornerAngleAt fi fj fk ≤ π ∧
    cos (cornerAngleAt fi fj fk) =
      dot (sub fj fi) (sub fk fi) / (norm (sub fj fi) * norm (sub fk fi)) ∧
    sin (cornerAngleAt fi fj fk) =
      norm (cross (sub fj fi) (sub fk fi)) / (norm (sub fj fi) * norm (sub fk fi)) := by
  set a := sub fj fi
  set b := sub fk fi
  have hL := norm2_cross_add_sq_dot a b
  have hy : 0 ≤ norm (cross a b) := norm_nonneg' _
  have hpa := (norm2_pos_iff a).mpr ha
  have hpb := (norm2_pos_iff b).mpr hb
  have hr : dot a b ^ 2 + norm (cross a b) ^ 2 = norm2 a * norm2 b := by
    rw [norm_sq]; linarith
  have hpos : 0 < dot a b ^ 2 + norm (cross a b) ^ 2 := by rw [hr]; positivity
  have hsq : √(dot a b ^ 2 + norm (cross a b) ^ 2) = norm a * norm b := by
    rw [hr, sqrt_mul (norm2_nonneg a), norm_def, norm_def]
  have hI := atan2_mem_Icc hy (dot a b)
  refine ⟨hI.1, hI.2, ?_, ?_⟩
  · show cos (Ddg.Dev.atan2 _ _) = _
    rw [cos_atan2 hy hpos, hsq]
  · show sin (Ddg.Dev.atan2 _ _) = _
    rw [sin_atan2 hy hpos, hsq]

/-- The corner angle is `arccos(a·b/(|a||b|))`, the Euclidean angle between the edge vectors.

[formula: dev-corner-angle] -/
theorem cornerAngleAt_eq_arccos (fi fj fk : V3) (ha : sub fj fi ≠ zero) (hb : sub fk fi ≠ zero) :
    cornerAngleAt fi fj fk =
      arccos (dot (sub fj fi) (sub fk fi) / (norm (sub fj fi) * norm (sub fk fi))) := by
  obtain ⟨h0, hπ, hc, -⟩ := cornerAngleAt_spec fi fj fk ha hb
  rw [← hc, arccos_cos h0 hπ]

/-! ## Derivatives (Eqs. 7–9, D2, D3) -/

section derivatives

variable {fi fj fk : ℝ → V3} {δi δj δk : V3} {t : ℝ}

/-- `dn = e_i × δ_i + e_j × δ_j + e_k × δ_k` (`n` is bilinear), `e_i = f_k − f_j`,
`e_j = f_i − f_k`, `e_k = f_j − f_i`. -/
theorem hasDerivAt3_faceNormal (hi : HasDerivAt3 fi δi t) (hj : HasDerivAt3 fj δj t)
    (hk : HasDerivAt3 fk δk t) :
    HasDerivAt3 (fun s => faceNormal (fi s) (fj s) (fk s))
      (add (add (cross (sub (fk t) (fj t)) δi) (cross (sub (fi t) (fk t)) δj))
        (cross (sub (fj t) (fi t)) δk)) t := by
  have := (hj.sub hi).cross (hk.sub hi)
  refine this.of_eq ?_
  ext <;> simp <;> ring

/-- **Area gradient** (Eq. 7): along any motion of the corners,
`dA = Σ_p (½ N × e_p) · δ_p`, i.e. `∇_{f_p} A = ½ N × e_p` for each corner `p` (`areaGradient`
taken in the cyclic order starting at `p`).

[formula: dev-area-gradient] -/
theorem hasDerivAt_faceArea (hi : HasDerivAt3 fi δi t) (hj : HasDerivAt3 fj δj t)
    (hk : HasDerivAt3 fk δk t) (hn : faceNormal (fi t) (fj t) (fk t) ≠ zero) :
    HasDerivAt (fun s => faceArea (fi s) (fj s) (fk s))
      (dot (areaGradient (fi t) (fj t) (fk t)) δi + dot (areaGradient (fj t) (fk t) (fi t)) δj +
        dot (areaGradient (fk t) (fi t) (fj t)) δk) t := by
  have h1 := ((hasDerivAt3_faceNormal hi hj hk).norm hn).div_const 2
  refine h1.congr_deriv ?_
  have hp := (norm_pos_iff' (faceNormal (fi t) (fj t) (fk t))).mpr hn
  unfold areaGradient
  rw [faceUnitNormal_cycle (fi t) (fj t) (fk t), faceUnitNormal_cycle (fj t) (fk t) (fi t),
    faceUnitNormal_cycle (fi t) (fj t) (fk t), faceUnitNormal_def, normalize_def]
  simp only [dot_def, add_x, add_y, add_z, smul_x, smul_y, smul_z, cross_x, cross_y, cross_z,
    sub_x, sub_y, sub_z, Oracle.one_eq, Oracle.two_eq]
  field_simp
  ring

/-- **Normal derivative** (Eq. 8 with the factor `1/(2A)`, D2): along any motion of the corners,
`dN = Σ_p (N · δ_p)(e_p × N)/(2A)`; for one corner, `∂N/∂f_p = (e_p × N) Nᵀ/(2A)`.

[formula: dev-normal-gradient] -/
theorem hasDerivAt3_faceUnitNormal (hi : HasDerivAt3 fi δi t) (hj : HasDerivAt3 fj δj t)
    (hk : HasDerivAt3 fk δk t) (hn : faceNormal (fi t) (fj t) (fk t) ≠ zero) :
    HasDerivAt3 (fun s => faceUnitNormal (fi s) (fj s) (fk s))
      (add (add (normalDerivative (fi t) (fj t) (fk t) δi)
        (normalDerivative (fj t) (fk t) (fi t) δj)) (normalDerivative (fk t) (fi t) (fj t) δk))
      t := by
  have h1 := (hasDerivAt3_faceNormal hi hj hk).normalize hn
  refine h1.of_eq ?_
  unfold normalDerivative
  rw [faceUnitNormal_cycle (fi t) (fj t) (fk t), faceUnitNormal_cycle (fj t) (fk t) (fi t),
    faceUnitNormal_cycle (fi t) (fj t) (fk t), faceArea_cycle (fi t) (fj t) (fk t),
    faceArea_cycle (fj t) (fk t) (fi t), faceArea_cycle (fi t) (fj t) (fk t)]
  have hA : 2 * faceArea (fi t) (fj t) (fk t) = norm (faceNormal (fi t) (fj t) (fk t)) := by
    unfold faceArea; simp only [Oracle.two_eq]; ring
  rw [hA, faceUnitNormal_def]
  set n := faceNormal (fi t) (fj t) (fk t)
  set N := normalize n with hNdef
  have hN : norm2 N = 1 := norm2_normalize n hn
  have hp := (norm_pos_iff' n).mpr hn
  -- N ⊥ each edge
  have hd := dot_faceNormal_edges (fi t) (fj t) (fk t)
  have hperp : ∀ e : V3, dot n e = 0 → dot N e = 0 := by
    intro e he
    rw [hNdef, normalize_def]
    rw [dot_def] at he ⊢
    simp only [smul_x, smul_y, smul_z]
    linear_combination (1 / norm n) * he
  have pi := proj_cross (δ := δi) hN (hperp _ hd.2.2.1)
  have pj := proj_cross (δ := δj) hN (hperp _ hd.2.2.2)
  have pk := proj_cross (δ := δk) hN (hperp _ hd.1)
  change _ = add (add (smul (dot N δi / norm n) (cross (sub (fk t) (fj t)) N))
    (smul (dot N δj / norm n) (cross (sub (fi t) (fk t)) N)))
    (smul (dot N δk / norm n) (cross (sub (fj t) (fi t)) N))
  have ex := congrArg Ddg.Point.x pi; have ey := congrArg Ddg.Point.y pi
  have ez := congrArg Ddg.Point.z pi
  have fx := congrArg Ddg.Point.x pj; have fy := congrArg Ddg.Point.y pj
  have fz := congrArg Ddg.Point.z pj
  have gx := congrArg Ddg.Point.x pk; have gy := congrArg Ddg.Point.y pk
  have gz := congrArg Ddg.Point.z pk
  simp only [sub_x, sub_y, sub_z, smul_x, smul_y, smul_z, cross_x, cross_y, cross_z,
    dot_def] at ex ey ez fx fy fz gx gy gz
  ext <;> simp only [sub_x, sub_y, sub_z, smul_x, smul_y, smul_z, cross_x, cross_y, cross_z,
    add_x, add_y, add_z, dot_def] <;> field_simp
  · linear_combination ex + fx + gx
  · linear_combination ey + fy + gy
  · linear_combination ez + fz + gz

/-- **Angle gradient** (Eq. 9 with the squared lengths, D3): along any motion of the corners,
`dθ_i = ∇_{f_i}θ_i · δ_i + ∇_{f_j}θ_i · δ_j + ∇_{f_k}θ_i · δ_k` with
`∇_{f_j}θ_i = N × (f_i − f_j)/|f_i − f_j|²`, `∇_{f_k}θ_i = N × (f_k − f_i)/|f_k − f_i|²`,
`∇_{f_i}θ_i = −(∇_{f_j}θ_i + ∇_{f_k}θ_i)`.

[formula: dev-angle-gradient] -/
theorem hasDerivAt_cornerAngleAt (hi : HasDerivAt3 fi δi t) (hj : HasDerivAt3 fj δj t)
    (hk : HasDerivAt3 fk δk t) (hn : faceNormal (fi t) (fj t) (fk t) ≠ zero) :
    HasDerivAt (fun s => cornerAngleAt (fi s) (fj s) (fk s))
      (dot (angleGradientI (fi t) (fj t) (fk t)) δi + dot (angleGradientJ (fi t) (fj t) (fk t)) δj +
        dot (angleGradientK (fi t) (fj t) (fk t)) δk) t := by
  have ha := hj.sub hi
  have hb := hk.sub hi
  have hy := (ha.cross hb).norm hn
  have hx := ha.dot hb
  have hpos := (norm_pos_iff' _).mpr hn
  have h1 := HasDerivAt.atan2 hx hy hpos
  refine h1.congr_deriv ?_
  obtain ⟨hab, hcb, -⟩ := ne_of_faceNormal_ne _ _ _ hn
  have hpa := (norm2_pos_iff _).mpr hab
  have hpc := (norm2_pos_iff _).mpr hcb
  set a := sub (fj t) (fi t) with ha_def
  set b := sub (fk t) (fi t) with hb_def
  have hn' : faceNormal (fi t) (fj t) (fk t) = cross a b := rfl
  have hfi : sub (fi t) (fj t) = neg a := by
    ext <;> simp only [ha_def, sub_x, sub_y, sub_z, neg_x, neg_y, neg_z] <;> ring
  have hsa : norm2 (neg a) = norm2 a := by simp [norm2_def]
  -- the target, linearized
  have hT : dot (angleGradientI (fi t) (fj t) (fk t)) δi +
      dot (angleGradientJ (fi t) (fj t) (fk t)) δj + dot (angleGradientK (fi t) (fj t) (fk t)) δk =
      1 / norm (cross a b) * (1 / norm2 a * dot (cross (cross a b) (neg a)) (sub δj δi) +
        1 / norm2 b * dot (cross (cross a b) b) (sub δk δi)) := by
    unfold angleGradientI angleGradientJ angleGradientK
    rw [hfi, faceUnitNormal_def, hn', normalize_def, hsa, ← hb_def]
    simp only [dot_add_l, dot_neg_l, dot_smul_l, cross_smul_l, dot_sub_r, Oracle.one_eq]
    ring
  rw [hT]
  rw [hn'] at hpos
  have hL := norm2_cross_add_sq_dot a b
  have hr2 := norm_sq (cross a b)
  have key := angle_poly a b (sub δj δi) (sub δk δi)
  generalize norm (cross a b) = r at hpos hr2 ⊢
  generalize dot (cross a b) (add (cross (sub δj δi) b) (cross a (sub δk δi))) = X at key ⊢
  generalize dot (cross (cross a b) (neg a)) (sub δj δi) = U at key ⊢
  generalize dot (cross (cross a b) b) (sub δk δi) = W at key ⊢
  generalize dot (sub δj δi) b + dot a (sub δk δi) = Y at key ⊢
  generalize norm2 (cross a b) = P at key hL hr2
  generalize dot a b = x at key hL ⊢
  generalize norm2 a = A at key hL hpa ⊢
  generalize norm2 b = B at key hL hpc ⊢
  have hden : x ^ 2 + r ^ 2 = A * B := by rw [hr2]; linarith
  rw [hden]
  have e1 : (x * (X / r) - r * Y) / (A * B) = (x * X - r ^ 2 * Y) / (r * (A * B)) := by
    field_simp
  rw [e1, hr2, key]
  field_simp

/-- **Adjoint (back-propagation) form** (SPEC §2): for sensitivities `y` of the unit normal,
`ζ` of `n`, and `s_i, s_j, s_k` of the three corner angles, the function
`F = y · N + ζ · n + s_i θ_i + s_j θ_j + s_k θ_k` of the face has, along any motion of its
corners, the derivative `g_i · δ_i + g_j · δ_j + g_k · δ_k` with `(g_i, g_j, g_k) =
backpropFace f_i f_j f_k y ζ s_i s_j s_k`, i.e. `g_p = w × e_p + (angle terms)`,
`w = (y − (y·N) N)/|n| + ζ`. (For the vertex normals, `ζ = ½ Σ_u z_u`.)

[formula: dev-backprop-face] -/
theorem hasDerivAt_backpropFace (hi : HasDerivAt3 fi δi t) (hj : HasDerivAt3 fj δj t)
    (hk : HasDerivAt3 fk δk t) (hn : faceNormal (fi t) (fj t) (fk t) ≠ zero)
    (y zeta : V3) (si sj sk : ℝ) :
    HasDerivAt (fun s => dot y (faceUnitNormal (fi s) (fj s) (fk s)) +
        dot zeta (faceNormal (fi s) (fj s) (fk s)) +
        si * cornerAngleAt (fi s) (fj s) (fk s) + sj * cornerAngleAt (fj s) (fk s) (fi s) +
        sk * cornerAngleAt (fk s) (fi s) (fj s))
      (dot (backpropFace (fi t) (fj t) (fk t) y zeta si sj sk).1 δi +
        dot (backpropFace (fi t) (fj t) (fk t) y zeta si sj sk).2.1 δj +
        dot (backpropFace (fi t) (fj t) (fk t) y zeta si sj sk).2.2 δk) t := by
  have hn2 : faceNormal (fj t) (fk t) (fi t) ≠ zero := by rwa [faceNormal_cycle]
  have hn3 : faceNormal (fk t) (fi t) (fj t) ≠ zero := by
    rwa [faceNormal_cycle, faceNormal_cycle]
  have dn := hasDerivAt3_faceNormal hi hj hk
  have dN := dn.normalize hn
  have dYN := dN.dot_const y
  have dZn := dn.dot_const zeta
  have ai := (hasDerivAt_cornerAngleAt hi hj hk hn).const_mul si
  have aj := (hasDerivAt_cornerAngleAt hj hk hi hn2).const_mul sj
  have ak := (hasDerivAt_cornerAngleAt hk hi hj hn3).const_mul sk
  have hsum := (((dYN.add dZn).add ai).add aj).add ak
  refine hsum.congr_deriv ?_
  unfold backpropFace backpropW
  simp only [faceUnitNormal_def]
  generalize normalize (faceNormal (fi t) (fj t) (fk t)) = N
  generalize norm (faceNormal (fi t) (fj t) (fk t)) = r
  generalize angleGradientI (fi t) (fj t) (fk t) = g1
  generalize angleGradientJ (fi t) (fj t) (fk t) = g2
  generalize angleGradientK (fi t) (fj t) (fk t) = g3
  generalize angleGradientI (fj t) (fk t) (fi t) = g4
  generalize angleGradientJ (fj t) (fk t) (fi t) = g5
  generalize angleGradientK (fj t) (fk t) (fi t) = g6
  generalize angleGradientI (fk t) (fi t) (fj t) = g7
  generalize angleGradientJ (fk t) (fi t) (fj t) = g8
  generalize angleGradientK (fk t) (fi t) (fj t) = g9
  simp only [dot_add_l, dot_add_r, dot_sub_l, dot_sub_r, dot_smul_l, dot_smul_r,
    dot_cross_l, Oracle.one_eq]
  ring

end derivatives

/-! ## Vertex normal -/

section vertexNormal

/-- `sumTo` is the sum over `Finset.range`. -/
theorem sumTo_eq_sum_range (n : ℕ) (f : ℕ → ℝ) : Ddg.sumTo n f = ∑ i ∈ Finset.range n, f i := by
  induction n with
  | zero => simp [Ddg.sumTo]
  | succ n ih => rw [Finset.sum_range_succ, ← ih]; rfl

/-- A cyclic shift does not change a sum over `0 … k − 1`. -/
theorem sum_range_succ_mod (k : ℕ) (g : ℕ → ℝ) :
    ∑ m ∈ Finset.range k, g ((m + 1) % k) = ∑ m ∈ Finset.range k, g m := by
  rcases k with _ | k
  · simp
  · rw [Finset.sum_range_succ, Finset.sum_range_succ']
    have h1 : (k + 1) % (k + 1) = 0 := Nat.mod_self _
    rw [h1]
    congr 1
    apply Finset.sum_congr rfl
    intro m hm
    rw [Finset.mem_range] at hm
    rw [Nat.mod_eq_of_lt (by omega)]

/-- `Σ_m (c m + (g m − g (m + 1 mod k))) = Σ_m c m`. -/
theorem sum_telescope (k : ℕ) (c g : ℕ → ℝ) :
    ∑ m ∈ Finset.range k, (c m + (g m - g ((m + 1) % k))) = ∑ m ∈ Finset.range k, c m := by
  rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, sum_range_succ_mod, sub_self, add_zero]

/-- **Vertex normal** (SPEC §1, App. B.5): the vector area `S_v = ½ Σ_m n_m` of the star of an
interior vertex `v` with the cyclic ring `w_0 … w_{k−1}` (faces `(v, w_m, w_{m+1 mod k})`) is
`½ Σ_m w_m × w_{m+1}`: it does not depend on the position of `v` (the terms in `f_v` telescope
around the cycle).

[formula: dev-vertex-normal] -/
theorem ringVectorArea_eq (v : V3) (k : ℕ) (ring : ℕ → V3) :
    ringVectorArea v k ring =
      smul (1 / 2) (vsumTo k (fun m => cross (ring m) (ring ((m + 1) % k)))) := by
  unfold ringVectorArea vectorArea
  simp only [Oracle.one_eq, Oracle.two_eq]
  refine congrArg (Ddg.Dev.smul (α := ℝ) _) ?_
  ext
  · show Ddg.sumTo k _ = Ddg.sumTo k _
    rw [sumTo_eq_sum_range, sumTo_eq_sum_range]
    refine Eq.trans ?_ (sum_telescope k (fun m => (cross (ring m) (ring ((m + 1) % k))).x)
      (fun m => v.y * (ring m).z - v.z * (ring m).y))
    refine Finset.sum_congr rfl (fun m _ => ?_)
    simp [faceNormal_def]; ring
  · show Ddg.sumTo k _ = Ddg.sumTo k _
    rw [sumTo_eq_sum_range, sumTo_eq_sum_range]
    refine Eq.trans ?_ (sum_telescope k (fun m => (cross (ring m) (ring ((m + 1) % k))).y)
      (fun m => v.z * (ring m).x - v.x * (ring m).z))
    refine Finset.sum_congr rfl (fun m _ => ?_)
    simp [faceNormal_def]; ring
  · show Ddg.sumTo k _ = Ddg.sumTo k _
    rw [sumTo_eq_sum_range, sumTo_eq_sum_range]
    refine Eq.trans ?_ (sum_telescope k (fun m => (cross (ring m) (ring ((m + 1) % k))).z)
      (fun m => v.x * (ring m).y - v.y * (ring m).x))
    refine Finset.sum_congr rfl (fun m _ => ?_)
    simp [faceNormal_def]; ring

/-- Unit vertex normal: `N_v = S_v/|S_v|` is a unit vector when `S_v ≠ 0`.

[formula: dev-vertex-normal] -/
theorem vertexNormal_unit (n : ℕ) (face : ℕ → V3 × V3 × V3) (h : vectorArea n face ≠ zero) :
    norm2 (vertexNormal n face) = 1 :=
  norm2_normalize _ h

end vertexNormal

end Dev

end DDG
