import DdgProofs.Delaunay.Surface.Segment
import DdgProofs.Intrinsic.Delaunay

/-!
# Delaunay layer 4, planar toolkit: angles at a point

The angle `angleAt v P Q` between the rays from `v` to `P` and to `Q`, defined from the three
lengths by the law of cosines (`cornerAngle` of `DdgProofs/Triangle.lean`: the angle opposite
`|PQ|` in the triangle with sides `|PQ|, |vP|, |vQ|`). What the crossing bound of
`Overlay.lean` needs:

* the angle does not change when `P`, `Q` move along their rays (`angleAt_smul`);
* **additivity along a line**: for `Q` strictly between `P` and `R`, and `v` off the line,
  `∠PvR = ∠PvQ + ∠QvR` (`angleAt_add`); hence the angles of consecutive points of a ray add up to
  the angle of the end points (`sum_angleAt_ray`), which is `< π` (`angleAt_lt_pi`): **a straight
  line through a fan of wedges at `v` subtends less than `π`**;
* the chord bound `|PQ| ≥ |vQ| sin ∠PvQ` (`dist2_ge_mul_sin`).

The proofs are algebraic: with `a = P − v`, `b = Q − v`, `c = R − v`, `cos ∠PvQ = a·b/(|a||b|)`
and `sin ∠PvQ = |a × b|/(|a||b|)`; additivity is the identities
`(a·b)(b·c) − (a×b)(b×c) = |b|²(a·c)` and, for `b = (1−μ)a + μc`,
`|a×b|(b·c) + (a·b)|b×c| = |a×c| |b|²`.
-/

noncomputable section

namespace DDG.Delaunay.Surface

open Real Set DDG.Intrinsic

/-- The angle at `v` between `P` and `Q`, from the three lengths. -/
def angleAt (v P Q : Pt) : ℝ := cornerAngle (dist2 P Q) (dist2 v P) (dist2 v Q)

theorem angleAt_comm (v P Q : Pt) : angleAt v P Q = angleAt v Q P := by
  unfold angleAt
  rw [dist2_comm, Intrinsic.cornerAngle_swap]

theorem angleAt_nonneg (v P Q : Pt) : 0 ≤ angleAt v P Q := arccos_nonneg _

theorem angleAt_le_pi (v P Q : Pt) : angleAt v P Q ≤ π := arccos_le_pi _

/-- The triangle `v P Q` is nondegenerate when `sda v P Q ≠ 0`. -/
theorem isTriangle_angleAt {v P Q : Pt} (h : sda v P Q ≠ 0) :
    IsTriangle (dist2 P Q) (dist2 v P) (dist2 v Q) := by
  have := (isTriangle_of_sda_ne h).swap
  rwa [dist2_comm Q v] at this

theorem angleAt_pos {v P Q : Pt} (h : sda v P Q ≠ 0) : 0 < angleAt v P Q :=
  cornerAngle_pos (isTriangle_angleAt h)

theorem angleAt_lt_pi {v P Q : Pt} (h : sda v P Q ≠ 0) : angleAt v P Q < π :=
  cornerAngle_lt_pi (isTriangle_angleAt h)

/-- `cos ∠PvQ = (P − v)·(Q − v) / (|vP| |vQ|)` for a nondegenerate triangle. -/
theorem cos_angleAt {v P Q : Pt} (h : sda v P Q ≠ 0) :
    cos (angleAt v P Q) = dot2 (P - v) (Q - v) / (dist2 v P * dist2 v Q) := by
  have hc : cos (angleAt v P Q) = _ := cos_cornerAngle (isTriangle_angleAt h)
  have key : dist2 v P ^ 2 + dist2 v Q ^ 2 - dist2 P Q ^ 2 = 2 * dot2 (P - v) (Q - v) := by
    rw [dist2_sq', dist2_sq', dist2_sq']
    simp only [dot2, Prod.fst_sub, Prod.snd_sub]
    ring
  rw [hc, key, mul_assoc, mul_div_mul_left _ _ two_ne_zero]

/-- `sin ∠PvQ = |sda v P Q| / (|vP| |vQ|)` for a nondegenerate triangle. -/
theorem sin_angleAt {v P Q : Pt} (h : sda v P Q ≠ 0) :
    sin (angleAt v P Q) = |sda v P Q| / (dist2 v P * dist2 v Q) := by
  have hs : sin (angleAt v P Q) = _ := sin_cornerAngle (isTriangle_angleAt h)
  rw [hs, ← triArea_rotate, two_mul_triArea_dist2]

private theorem dist2_smul (v P : Pt) {a : ℝ} (ha : 0 < a) :
    dist2 v (v + a • (P - v)) = a * dist2 v P := by
  unfold dist2
  rw [show (v + a • (P - v)).1 - v.1 = a * (P.1 - v.1) by simp,
    show (v + a • (P - v)).2 - v.2 = a * (P.2 - v.2) by simp,
    show (a * (P.1 - v.1)) ^ 2 + (a * (P.2 - v.2)) ^ 2 =
      a ^ 2 * ((P.1 - v.1) ^ 2 + (P.2 - v.2) ^ 2) by ring,
    Real.sqrt_mul (sq_nonneg a), Real.sqrt_sq ha.le]

/-- The angle does not change when the points move along their rays from `v`. -/
theorem angleAt_smul (v P Q : Pt) {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    angleAt v (v + a • (P - v)) (v + b • (Q - v)) = angleAt v P Q := by
  unfold angleAt cornerAngle
  rw [dist2_smul v P ha, dist2_smul v Q hb]
  congr 1
  have e1 := dist2_sq' v P
  have e2 := dist2_sq' v Q
  have e3 := dist2_sq' P Q
  have e4 := dist2_sq' (v + a • (P - v)) (v + b • (Q - v))
  simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, Prod.fst_sub,
    Prod.snd_sub, smul_eq_mul] at e4
  have hnum : (a * dist2 v P) ^ 2 + (b * dist2 v Q) ^ 2 -
      dist2 (v + a • (P - v)) (v + b • (Q - v)) ^ 2 =
      (a * b) * (dist2 v P ^ 2 + dist2 v Q ^ 2 - dist2 P Q ^ 2) := by
    linear_combination (a ^ 2 - a * b) * e1 + (b ^ 2 - a * b) * e2 - e4 + a * b * e3
  rw [hnum, show 2 * (a * dist2 v P) * (b * dist2 v Q) = (a * b) * (2 * dist2 v P * dist2 v Q) by
    ring, mul_div_mul_left _ _ (mul_pos ha hb).ne']

/-- **Additivity along a line**: for `Q` strictly between `P` and `R` and `v` off the line,
`∠PvR = ∠PvQ + ∠QvR`. -/
theorem angleAt_add {v P Q R : Pt} (h : sda v P R ≠ 0) (hQ : Q ∈ openSegment ℝ P R) :
    angleAt v P R = angleAt v P Q + angleAt v Q R := by
  rw [openSegment_eq_image] at hQ
  obtain ⟨μ, ⟨hμ0, hμ1⟩, rfl⟩ := hQ
  set Q := (1 - μ) • P + μ • R with hQdef
  have hμ1' : 0 < 1 - μ := sub_pos.2 hμ1
  have hPQ : sda v P Q = μ * sda v P R := by
    simp only [hQdef, sda, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    ring
  have hQR : sda v Q R = (1 - μ) * sda v P R := by
    simp only [hQdef, sda, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    ring
  have h1 : sda v P Q ≠ 0 := by rw [hPQ]; exact mul_ne_zero hμ0.ne' h
  have h2 : sda v Q R ≠ 0 := by rw [hQR]; exact mul_ne_zero hμ1'.ne' h
  have dP : 0 < dist2 v P := dist2_pos_of_sda_ne h
  have dR : 0 < dist2 v R := dist2_pos_of_sda_ne
    (show sda v R P ≠ 0 by rw [sda_swap_right]; exact neg_ne_zero.2 h)
  have dQ : 0 < dist2 v Q := dist2_pos_of_sda_ne h2
  have eQ : dist2 v Q ^ 2 = (Q - v).1 ^ 2 + (Q - v).2 ^ 2 := dist2_sq' v Q
  have aPQ : |sda v P Q| = μ * |sda v P R| := by rw [hPQ, abs_mul, abs_of_pos hμ0]
  have aQR : |sda v Q R| = (1 - μ) * |sda v P R| := by rw [hQR, abs_mul, abs_of_pos hμ1']
  have aa' : |sda v P R| ^ 2 = sda v P R * sda v P R := by rw [sq_abs]; ring
  have I1 : dot2 (P - v) (Q - v) * dot2 (Q - v) (R - v) - μ * (1 - μ) * (sda v P R * sda v P R) =
      ((Q - v).1 ^ 2 + (Q - v).2 ^ 2) * dot2 (P - v) (R - v) := by
    simp only [hQdef, dot2, sda, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
      smul_eq_mul, Prod.fst_sub, Prod.snd_sub]
    ring
  have I2 : μ * dot2 (Q - v) (R - v) + (1 - μ) * dot2 (P - v) (Q - v) =
      (Q - v).1 ^ 2 + (Q - v).2 ^ 2 := by
    simp only [hQdef, dot2, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
      smul_eq_mul, Prod.fst_sub, Prod.snd_sub]
    ring
  have hcos : cos (angleAt v P Q + angleAt v Q R) = cos (angleAt v P R) := by
    rw [cos_add, cos_angleAt h1, cos_angleAt h2, sin_angleAt h1, sin_angleAt h2, cos_angleAt h,
      aPQ, aQR]
    field_simp
    linear_combination I1 - dot2 (P - v) (R - v) * eQ - μ * (1 - μ) * aa'
  have hsin : sin (angleAt v P Q + angleAt v Q R) = sin (angleAt v P R) := by
    rw [sin_add, cos_angleAt h1, cos_angleAt h2, sin_angleAt h1, sin_angleAt h2, sin_angleAt h,
      aPQ, aQR]
    field_simp
    linear_combination I2 - eQ
  have hpos : 0 < sin (angleAt v P R) := sin_pos_of_pos_of_lt_pi (angleAt_pos h) (angleAt_lt_pi h)
  have hlt : angleAt v P Q + angleAt v Q R < π := by
    by_contra hge
    push Not at hge
    have h3 : angleAt v P Q + angleAt v Q R - π ≤ π := by
      linarith [angleAt_lt_pi h1, angleAt_lt_pi h2]
    have h4 := sin_nonneg_of_nonneg_of_le_pi (sub_nonneg.2 hge) h3
    rw [sin_sub_pi] at h4
    linarith
  exact injOn_cos ⟨angleAt_nonneg v P R, angleAt_le_pi v P R⟩
    ⟨add_nonneg (angleAt_nonneg v P Q) (angleAt_nonneg v Q R), hlt.le⟩ hcos.symm

/-- The chord bound: `|PQ| ≥ |vQ| sin ∠PvQ`. -/
theorem dist2_ge_mul_sin {v P Q : Pt} (h : sda v P Q ≠ 0) :
    dist2 v Q * sin (angleAt v P Q) ≤ dist2 P Q := by
  have hP : 0 < dist2 v P := dist2_pos_of_sda_ne h
  have hQ : 0 < dist2 v Q := dist2_pos_of_sda_ne
    (show sda v Q P ≠ 0 by rw [sda_swap_right]; exact neg_ne_zero.2 h)
  have hc : cos (angleAt v P Q) = _ := cos_cornerAngle (isTriangle_angleAt h)
  have hc' : 2 * dist2 v P * dist2 v Q * cos (angleAt v P Q) =
      dist2 v P ^ 2 + dist2 v Q ^ 2 - dist2 P Q ^ 2 := by
    rw [hc]; field_simp
  have hsin : 0 ≤ sin (angleAt v P Q) :=
    sin_nonneg_of_nonneg_of_le_pi (angleAt_nonneg v P Q) (angleAt_le_pi v P Q)
  have hs2 := sin_sq (angleAt v P Q)
  have hsq : (dist2 v Q * sin (angleAt v P Q)) ^ 2 ≤ dist2 P Q ^ 2 := by
    nlinarith [sq_nonneg (dist2 v P - dist2 v Q * cos (angleAt v P Q))]
  have h0 : 0 ≤ dist2 P Q := dist2_nonneg _ _
  have h1 : 0 ≤ dist2 v Q * sin (angleAt v P Q) := mul_nonneg hQ.le hsin
  nlinarith

/-- `sda` of `v` with two points of a ray: `(t₁ − t₀) · cross2 d (v − p)`, up to sign:
`sda v (ray p d t₀) (ray p d t₁) = (t₁ − t₀) * cross2 d (v - p)`. -/
theorem sda_ray (v p d : Pt) (t₀ t₁ : ℝ) :
    sda v (ray p d t₀) (ray p d t₁) = (t₁ - t₀) * cross2 d (v - p) := by
  simp only [sda, ray, cross2, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    smul_eq_mul, Prod.fst_sub, Prod.snd_sub]
  ring

private theorem angleAt_self {v x : Pt} (hx : 0 < dist2 v x) : angleAt v x x = 0 := by
  unfold angleAt cornerAngle
  have h0 : dist2 x x = 0 := by simp [dist2]
  rw [h0, show (dist2 v x ^ 2 + dist2 v x ^ 2 - 0 ^ 2) / (2 * dist2 v x * dist2 v x) = 1 by
    rw [div_eq_one_iff_eq (by positivity)]; ring, arccos_one]

private theorem dist2_ray_pos {v p d : Pt} (hv : cross2 d (v - p) ≠ 0) (t : ℝ) :
    0 < dist2 v (ray p d t) := by
  apply dist2_pos_of_sda_ne (c := ray p d (t + 1))
  rw [sda_ray]
  simpa using hv

private theorem ray_mem_openSegment (p d : Pt) {t₀ t t₁ : ℝ} (h0 : t₀ < t) (h1 : t < t₁) :
    ray p d t ∈ openSegment ℝ (ray p d t₀) (ray p d t₁) := by
  rw [openSegment_eq_image]
  have hne : t₁ - t₀ ≠ 0 := by linarith
  refine ⟨(t - t₀) / (t₁ - t₀), ⟨div_pos (by linarith) (by linarith),
    (div_lt_one (by linarith)).2 (by linarith)⟩, ?_⟩
  ext <;> simp only [ray, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul] <;>
    field_simp <;> ring

/-- **The angles of consecutive points of a ray add up to the angle of its end points**: for
`t 0 < t 1 < ⋯ < t m` and `v` off the line, `Σ_i ∠(z(t i), v, z(t (i+1))) = ∠(z(t 0), v, z(t m))`.
-/
theorem sum_angleAt_ray {v p d : Pt} (hv : cross2 d (v - p) ≠ 0) (m : ℕ) (t : Fin (m + 1) → ℝ)
    (hmono : StrictMono t) :
    ∑ i : Fin m, angleAt v (ray p d (t i.castSucc)) (ray p d (t i.succ)) =
      angleAt v (ray p d (t 0)) (ray p d (t (Fin.last m))) := by
  induction m with
  | zero =>
    simp only [Finset.univ_eq_empty, Finset.sum_empty]
    exact (angleAt_self (dist2_ray_pos hv _)).symm
  | succ m ih =>
    rw [Fin.sum_univ_castSucc]
    have ih' := ih (t ∘ Fin.castSucc) (hmono.comp (Fin.strictMono_castSucc))
    simp only [Function.comp_apply, Fin.castSucc_zero] at ih'
    simp only [Fin.succ_castSucc, ih', Fin.succ_last]
    rcases Nat.eq_zero_or_pos m with rfl | hm
    · have : (Fin.last 0).castSucc = (0 : Fin 2) := rfl
      rw [this, angleAt_self (dist2_ray_pos hv _), zero_add]
    · have ha : t 0 < t (Fin.last m).castSucc := hmono (by simp [Fin.lt_def]; omega)
      have hb : t (Fin.last m).castSucc < t (Fin.last (m + 1)) := hmono (by simp [Fin.lt_def])
      symm
      apply angleAt_add
      · rw [sda_ray]
        exact mul_ne_zero (sub_ne_zero.2 (ha.trans hb).ne') hv
      · exact ray_mem_openSegment p d ha hb

/-- Hence a straight line through a fan of `m` wedges at `v` subtends less than `π`. -/
theorem sum_angleAt_ray_lt_pi {v p d : Pt} (hv : cross2 d (v - p) ≠ 0) (m : ℕ)
    (t : Fin (m + 1) → ℝ) (hmono : StrictMono t) :
    ∑ i : Fin m, angleAt v (ray p d (t i.castSucc)) (ray p d (t i.succ)) < π := by
  rcases Nat.eq_zero_or_pos m with rfl | hm
  · simp [pi_pos]
  rw [sum_angleAt_ray hv m t hmono]
  apply angleAt_lt_pi
  rw [sda_ray]
  refine mul_ne_zero ?_ hv
  have h : t 0 < t (Fin.last m) := hmono (by simp [Fin.lt_def, hm])
  exact sub_ne_zero.2 h.ne'

end DDG.Delaunay.Surface
