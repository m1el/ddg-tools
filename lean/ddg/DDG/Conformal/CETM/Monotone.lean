import DDG.Conformal.CETM.Flatten

/-!
# Intrinsic simplification: the angle sum of the vertex being flattened is monotone

`crates/ddg-intrinsic-simplification/SPEC.md` §3.4, S9, §13 (`iem-flatten-monotone`); the paper
Sec. 4.1; Springborn, Schröder, Pinkall 2008, Sec. 3.1 (extended angles, which make the angle
sum a continuous function of the scale factors on all of `ℝ`).

The model is the one of `Flatten.lean`. Every face with positive sides has one of four shapes
(`shape`), with `s = e^u`:

* no corner at `i`: its share of `Θ_i` is `0`;
* one corner at `i`: `arccos(P − Q/s)` with `Q > 0`, `P ≥ 1` (`OneCorner`): non-increasing,
  `→ π` as `u → −∞` (eventually equal: the face degenerates with the angle `π` at `i`), `→ 0`
  as `u → +∞`;
* two corners at `i` (a self-edge `ii`): `π − arccos(P − Q s)` (`TwoCorners`), the angle sum `π`
  minus the angle at the other corner: non-increasing, `→ π` at `−∞`, `→ 0` at `+∞` (eventually);
* a self-face: `π`, constant.

The argument of `arccos` lies strictly between `−1` and `1` exactly where the face is a strict
triangle; there the share is strictly decreasing. Summing over the faces gives the statements of
SPEC §3.4 (`iem-flatten-monotone`).
-/

noncomputable section

namespace DDG

namespace IntrinsicSimplification

open Real Filter Topology

/-- `arccos` is strictly decreasing between two arguments of which the smaller is below `1` and
the larger above `−1` (the clamp of `Real.arccos` makes it constant outside `[−1, 1]`). -/
theorem arccos_lt_arccos_of_lt {x y : ℝ} (hxy : x < y) (hx : x < 1) (hy : -1 < y) :
    arccos y < arccos x := by
  have ex : arccos x = arccos (max x (-1)) := by
    rcases le_total x (-1) with h | h
    · rw [max_eq_right h, arccos_eq_pi.mpr h, arccos_neg_one]
    · rw [max_eq_left h]
  have ey : arccos y = arccos (min y 1) := by
    rcases le_total y 1 with h | h
    · rw [min_eq_left h]
    · rw [min_eq_right h, arccos_eq_zero.mpr h, arccos_one]
  rw [ex, ey]
  apply arccos_lt_arccos (le_max_right _ _)
  · rcases le_total x (-1) with h | h
    · rw [max_eq_right h]; exact lt_min hy (by norm_num)
    · rw [max_eq_left h]; exact lt_min hxy hx
  · exact min_le_right _ _

theorem one_le_ratio {b c : ℝ} (hb : 0 < b) (hc : 0 < c) : 1 ≤ (b ^ 2 + c ^ 2) / (2 * b * c) := by
  rw [le_div_iff₀ (by positivity)]
  nlinarith [sq_nonneg (b - c)]

theorem exp_half_sq (u : ℝ) : exp (u / 2) ^ 2 = exp u := by
  rw [sq, ← exp_add]; ring_nf

/-- The cosine argument of the angle at a corner at `i` whose two sides are spokes. -/
theorem ratio_one {a b c : ℝ} (hb : 0 < b) (hc : 0 < c) (u : ℝ) :
    ((b * exp (u / 2)) ^ 2 + (c * exp (u / 2)) ^ 2 - a ^ 2) /
        (2 * (b * exp (u / 2)) * (c * exp (u / 2))) =
      (b ^ 2 + c ^ 2) / (2 * b * c) - a ^ 2 / (2 * b * c) * exp (-u) := by
  have e := exp_half_sq u
  have hp := exp_pos (u / 2)
  rw [exp_neg, ← e]
  field_simp

/-- The cosine argument of the angle at the corner opposite a self-edge `ii`. -/
theorem ratio_two {a b c : ℝ} (ha : 0 < a) (hb : 0 < b) (u : ℝ) :
    ((a * exp (u / 2)) ^ 2 + (b * exp (u / 2)) ^ 2 - (c * exp u) ^ 2) /
        (2 * (a * exp (u / 2)) * (b * exp (u / 2))) =
      (a ^ 2 + b ^ 2) / (2 * a * b) - c ^ 2 / (2 * a * b) * exp u := by
  have e := exp_half_sq u
  have hp := exp_pos (u / 2)
  rw [← e]
  field_simp

namespace StarFace

variable (f : StarFace)

/-- A face with one corner at `i`: its share of `Θ_i` is `arccos(P − Q e^{−u})`. -/
structure OneCorner (P Q : ℝ) : Prop where
  hQ : 0 < Q
  hP : 1 ≤ P
  eq : ∀ u, f.angleAtI u = arccos (P - Q * exp (-u))
  valid : ∀ u, f.Valid u → -1 < P - Q * exp (-u) ∧ P - Q * exp (-u) < 1

/-- A face with two corners at `i`: its share of `Θ_i` is `π − arccos(P − Q e^{u})`. -/
structure TwoCorners (P Q : ℝ) : Prop where
  hQ : 0 < Q
  hP : 1 ≤ P
  eq : ∀ u, f.angleAtI u = π - arccos (P - Q * exp u)
  valid : ∀ u, f.Valid u → -1 < P - Q * exp u ∧ P - Q * exp u < 1

/-- The face with its corners relabelled `1, 2, 0 ↦ 0, 1, 2`. -/
def rot : StarFace := ⟨f.l1, f.l2, f.l0, f.c1, f.c2, f.c0⟩

theorem rot_angleAtI (u : ℝ) : f.rot.angleAtI u = f.angleAtI u := by
  simp only [angleAtI, rot, θ0, θ1, θ2, L0, L1, L2]
  ring

theorem rot_valid {u : ℝ} (h : f.Valid u) : f.rot.Valid u := by
  have := h.rotate
  simpa [Valid, rot, L0, L1, L2] using this

theorem rot_pos (h : f.Pos) : f.rot.Pos := ⟨h.2.1, h.2.2, h.1⟩

theorem oneCorner_of_rot {P Q : ℝ} (h : f.rot.OneCorner P Q) : f.OneCorner P Q :=
  ⟨h.hQ, h.hP, fun u => by rw [← f.rot_angleAtI]; exact h.eq u,
    fun u hv => h.valid u (f.rot_valid hv)⟩

theorem twoCorners_of_rot {P Q : ℝ} (h : f.rot.TwoCorners P Q) : f.TwoCorners P Q :=
  ⟨h.hQ, h.hP, fun u => by rw [← f.rot_angleAtI]; exact h.eq u,
    fun u hv => h.valid u (f.rot_valid hv)⟩

theorem exp_self_edge (u : ℝ) : exp (u * (1 + 1) / 2) = exp u := by congr 1; ring

/-- Corner 0 is the only corner at `i`. -/
theorem oneCorner_base (hp : f.Pos) (h0 : f.c0 = true) (h1 : f.c1 = false) (h2 : f.c2 = false) :
    f.OneCorner ((f.l1 ^ 2 + f.l2 ^ 2) / (2 * f.l1 * f.l2)) (f.l0 ^ 2 / (2 * f.l1 * f.l2)) := by
  obtain ⟨p0, p1, p2⟩ := hp
  have hL0 : ∀ u, f.L0 u = f.l0 := fun u => by simp [L0, h1, h2]
  have hL1 : ∀ u, f.L1 u = f.l1 * exp (u / 2) := fun u => by simp [L1, h2, h0]
  have hL2 : ∀ u, f.L2 u = f.l2 * exp (u / 2) := fun u => by simp [L2, h0, h1]
  have hA : ∀ u, f.angleAtI u = f.θ0 u := fun u => by simp [angleAtI, h0, h1, h2]
  refine ⟨by positivity, one_le_ratio p1 p2, fun u => ?_, fun u hv => ?_⟩
  · rw [hA, θ0, hL0, hL1, hL2, cornerAngle, ratio_one p1 p2]
  · have := cos_arg_bounds hv
    rwa [hL0, hL1, hL2, ratio_one p1 p2] at this

/-- Corners 0 and 1 are at `i`; side 2 is the self-edge. -/
theorem twoCorners_base (hp : f.Pos) (h0 : f.c0 = true) (h1 : f.c1 = true) (h2 : f.c2 = false) :
    f.TwoCorners ((f.l0 ^ 2 + f.l1 ^ 2) / (2 * f.l0 * f.l1)) (f.l2 ^ 2 / (2 * f.l0 * f.l1)) := by
  obtain ⟨p0, p1, p2⟩ := hp
  have hL0 : ∀ u, f.L0 u = f.l0 * exp (u / 2) := fun u => by simp [L0, h1, h2]
  have hL1 : ∀ u, f.L1 u = f.l1 * exp (u / 2) := fun u => by simp [L1, h2, h0]
  have hL2 : ∀ u, f.L2 u = f.l2 * exp u := fun u => by simp [L2, h0, h1, exp_self_edge]
  have hA : ∀ u, f.angleAtI u = π - f.θ2 u := fun u => by
    have hs := cornerAngle_sum_pos (a := f.L0 u) (b := f.L1 u) (c := f.L2 u)
      (by rw [hL0]; positivity) (by rw [hL1]; positivity) (by rw [hL2]; positivity)
    simp only [angleAtI, h0, h1, h2, ind_true, ind_false, θ0, θ1, θ2]
    linarith
  refine ⟨by positivity, one_le_ratio p0 p1, fun u => ?_, fun u hv => ?_⟩
  · rw [hA, θ2, hL0, hL1, hL2, cornerAngle, ratio_two p0 p1]
  · have := cos_arg_bounds hv.rotate.rotate
    rwa [hL0, hL1, hL2, ratio_two p0 p1] at this

/-- **The shapes.** A face of the star with positive sides that is not a self-face has one or
two corners at `i` and the closed form `OneCorner` or `TwoCorners`.

[formula: iem-flatten-monotone] -/
theorem shape (hp : f.Pos) (hs : f.InStar) (hn : ¬ f.SelfFace) :
    (∃ P Q, f.OneCorner P Q) ∨ (∃ P Q, f.TwoCorners P Q) := by
  obtain ⟨l0, l1, l2, c0, c1, c2⟩ := f
  unfold InStar at hs
  unfold SelfFace at hn
  simp only at hs hn
  cases c0 <;> cases c1 <;> cases c2 <;> simp at hs hn
  · -- corner 2 only
    left
    exact ⟨_, _, oneCorner_of_rot _ (oneCorner_of_rot _
      (oneCorner_base _ (rot_pos _ (rot_pos _ hp)) rfl rfl rfl))⟩
  · -- corner 1 only
    left
    exact ⟨_, _, oneCorner_of_rot _ (oneCorner_base _ (rot_pos _ hp) rfl rfl rfl)⟩
  · -- corners 1 and 2
    right
    exact ⟨_, _, twoCorners_of_rot _ (twoCorners_base _ (rot_pos _ hp) rfl rfl rfl)⟩
  · -- corner 0 only
    left
    exact ⟨_, _, oneCorner_base _ hp rfl rfl rfl⟩
  · -- corners 2 and 0
    right
    exact ⟨_, _, twoCorners_of_rot _ (twoCorners_of_rot _
      (twoCorners_base _ (rot_pos _ (rot_pos _ hp)) rfl rfl rfl))⟩
  · -- corners 0 and 1
    right
    exact ⟨_, _, twoCorners_base _ hp rfl rfl rfl⟩

theorem angleAtI_self_face (hp : f.Pos) (h : f.SelfFace) (u : ℝ) : f.angleAtI u = π := by
  obtain ⟨h0, h1, h2⟩ := h
  have e0 : f.L0 u = f.l0 * exp u := by simp [L0, h1, h2, exp_self_edge]
  have e1 : f.L1 u = f.l1 * exp u := by simp [L1, h2, h0, exp_self_edge]
  have e2 : f.L2 u = f.l2 * exp u := by simp [L2, h0, h1, exp_self_edge]
  have := cornerAngle_sum_pos (a := f.L0 u) (b := f.L1 u) (c := f.L2 u)
    (by rw [e0]; exact mul_pos hp.1 (exp_pos u)) (by rw [e1]; exact mul_pos hp.2.1 (exp_pos u))
    (by rw [e2]; exact mul_pos hp.2.2 (exp_pos u))
  simpa [angleAtI, h0, h1, h2, θ0, θ1, θ2] using this

theorem angleAtI_not_in_star (h : ¬ f.InStar) (u : ℝ) : f.angleAtI u = 0 := by
  unfold InStar at h
  push Not at h
  simp [angleAtI, h.1, h.2.1, h.2.2]

/-! ### One face -/

theorem continuous_angleAtI (hp : f.Pos) : Continuous f.angleAtI := by
  by_cases hs : f.InStar
  · by_cases hn : f.SelfFace
    · rw [funext (f.angleAtI_self_face hp hn)]; exact continuous_const
    · rcases f.shape hp hs hn with ⟨P, Q, h⟩ | ⟨P, Q, h⟩
      · rw [funext h.eq]; exact continuous_arccos.comp (by fun_prop)
      · rw [funext h.eq]; exact continuous_const.sub (continuous_arccos.comp (by fun_prop))
  · rw [funext (f.angleAtI_not_in_star hs)]; exact continuous_const

theorem antitone_angleAtI (hp : f.Pos) : Antitone f.angleAtI := by
  by_cases hs : f.InStar
  · by_cases hn : f.SelfFace
    · rw [funext (f.angleAtI_self_face hp hn)]; exact antitone_const
    · rcases f.shape hp hs hn with ⟨P, Q, h⟩ | ⟨P, Q, h⟩
      · intro u v huv
        rw [h.eq, h.eq]
        apply arccos_le_arccos
        have := exp_le_exp.mpr (neg_le_neg huv)
        nlinarith [h.hQ]
      · intro u v huv
        rw [h.eq, h.eq]
        have := exp_le_exp.mpr huv
        have : arccos (P - Q * exp u) ≤ arccos (P - Q * exp v) :=
          arccos_le_arccos (by nlinarith [h.hQ])
        linarith
  · rw [funext (f.angleAtI_not_in_star hs)]; exact antitone_const

/-- **Strict decrease of one face's share** between `u₁ < u₂` when the face is a strict triangle
at one of them. -/
theorem angleAtI_lt (hp : f.Pos) (hs : f.InStar) (hn : ¬ f.SelfFace) {u₁ u₂ : ℝ}
    (h12 : u₁ < u₂) (hv : f.Valid u₁ ∨ f.Valid u₂) : f.angleAtI u₂ < f.angleAtI u₁ := by
  rcases f.shape hp hs hn with ⟨P, Q, h⟩ | ⟨P, Q, h⟩
  · rw [h.eq, h.eq]
    have he := exp_lt_exp.mpr (neg_lt_neg h12)
    have hlt : P - Q * exp (-u₁) < P - Q * exp (-u₂) := by nlinarith [h.hQ]
    rcases hv with hv | hv
    · have := h.valid u₁ hv
      exact arccos_lt_arccos_of_lt hlt this.2 (by linarith [this.1])
    · have := h.valid u₂ hv
      exact arccos_lt_arccos_of_lt hlt (by linarith [this.2]) this.1
  · rw [h.eq, h.eq]
    have he := exp_lt_exp.mpr h12
    have hlt : P - Q * exp u₂ < P - Q * exp u₁ := by nlinarith [h.hQ]
    rcases hv with hv | hv
    · have := h.valid u₁ hv
      have := arccos_lt_arccos_of_lt hlt (by linarith [this.2]) this.1
      linarith
    · have := h.valid u₂ hv
      have := arccos_lt_arccos_of_lt hlt this.2 (by linarith [this.1])
      linarith

theorem tendsto_angleAtI_atBot (hp : f.Pos) (hs : f.InStar) :
    Tendsto f.angleAtI atBot (𝓝 π) := by
  by_cases hn : f.SelfFace
  · rw [funext (f.angleAtI_self_face hp hn)]; exact tendsto_const_nhds
  rcases f.shape hp hs hn with ⟨P, Q, h⟩ | ⟨P, Q, h⟩
  · have h1 : Tendsto (fun u : ℝ => exp (-u)) atBot atTop :=
      tendsto_exp_atTop.comp tendsto_neg_atBot_atTop
    have h2 : Tendsto (fun u : ℝ => P - Q * exp (-u)) atBot atBot := by
      have := (h1.const_mul_atTop h.hQ)
      have := tendsto_neg_atTop_atBot.comp this
      exact (tendsto_atBot_add_const_left _ P this).congr (fun u => by simp [sub_eq_add_neg])
    have h3 : ∀ᶠ u in atBot, f.angleAtI u = π := by
      filter_upwards [h2.eventually (eventually_le_atBot (-1))] with u hu
      rw [h.eq, arccos_eq_pi.mpr hu]
    exact tendsto_const_nhds.congr' (EventuallyEq.symm h3)
  · have h1 : Tendsto (fun u : ℝ => P - Q * exp u) atBot (𝓝 (P - Q * 0)) :=
      (tendsto_exp_atBot.const_mul Q).const_sub P
    have h2 := (continuous_arccos.tendsto _).comp h1
    rw [mul_zero, sub_zero, arccos_eq_zero.mpr h.hP] at h2
    have h3 := h2.const_sub π
    rw [sub_zero] at h3
    exact h3.congr (fun u => (h.eq u).symm)

theorem tendsto_angleAtI_atTop (hp : f.Pos) (hn : ¬ f.SelfFace) :
    Tendsto f.angleAtI atTop (𝓝 0) := by
  by_cases hs : f.InStar
  swap
  · rw [funext (f.angleAtI_not_in_star hs)]; exact tendsto_const_nhds
  rcases f.shape hp hs hn with ⟨P, Q, h⟩ | ⟨P, Q, h⟩
  · have h1 : Tendsto (fun u : ℝ => P - Q * exp (-u)) atTop (𝓝 (P - Q * 0)) :=
      (tendsto_exp_neg_atTop_nhds_zero.const_mul Q).const_sub P
    have h2 := (continuous_arccos.tendsto _).comp h1
    rw [mul_zero, sub_zero, arccos_eq_zero.mpr h.hP] at h2
    exact h2.congr (fun u => (h.eq u).symm)
  · have h2 : Tendsto (fun u : ℝ => P - Q * exp u) atTop atBot := by
      have := tendsto_neg_atTop_atBot.comp (tendsto_exp_atTop.const_mul_atTop h.hQ)
      exact (tendsto_atBot_add_const_left _ P this).congr (fun u => by simp [sub_eq_add_neg])
    have h3 : ∀ᶠ u in atTop, f.angleAtI u = 0 := by
      filter_upwards [h2.eventually (eventually_le_atBot (-1))] with u hu
      rw [h.eq, arccos_eq_pi.mpr hu, sub_self]
    exact tendsto_const_nhds.congr' (EventuallyEq.symm h3)

end StarFace

/-! ### The star -/

variable {ι : Type*} [Fintype ι] (S : ι → StarFace)

/-- **`Θ_i` is continuous** in the flattening scale `u` (extended angles).

[formula: iem-flatten-monotone] -/
theorem angleSum_continuous (hp : ∀ f, (S f).Pos) : Continuous (angleSum S) := by
  unfold angleSum
  exact continuous_finsetSum _ (fun f _ => (S f).continuous_angleAtI (hp f))

/-- **`Θ_i` is non-increasing** on all of `ℝ`.

[formula: iem-flatten-monotone] -/
theorem angleSum_antitone (hp : ∀ f, (S f).Pos) : Antitone (angleSum S) := by
  intro u v huv
  unfold angleSum
  exact Finset.sum_le_sum (fun f _ => (S f).antitone_angleAtI (hp f) huv)

/-- **Strictly decreasing where some face is valid**: for `u₁ < u₂`, if some face of the star that
is not a self-face is a strict triangle at `u₁` or at `u₂`, then `Θ_i(u₂) < Θ_i(u₁)`.

[formula: iem-flatten-monotone] -/
theorem angleSum_strictAnti_of_valid (hp : ∀ f, (S f).Pos) {u₁ u₂ : ℝ} (h12 : u₁ < u₂) (g : ι)
    (hs : (S g).InStar) (hn : ¬ (S g).SelfFace) (hv : (S g).Valid u₁ ∨ (S g).Valid u₂) :
    angleSum S u₂ < angleSum S u₁ := by
  unfold angleSum
  exact Finset.sum_lt_sum (fun f _ => (S f).antitone_angleAtI (hp f) h12.le)
    ⟨g, Finset.mem_univ _, (S g).angleAtI_lt (hp g) hs hn h12 hv⟩

/-- **`Θ_i → deg·π` as `u → −∞`**, `deg` the number of faces of the star (each face tends to the
degenerate face with the angle `π` at `i`).

[formula: iem-flatten-monotone] -/
theorem tendsto_angleSum_atBot (hp : ∀ f, (S f).Pos) (hs : ∀ f, (S f).InStar) :
    Tendsto (angleSum S) atBot (𝓝 (Fintype.card ι * π)) := by
  have := tendsto_finsetSum (Finset.univ : Finset ι)
    (fun f _ => (S f).tendsto_angleAtI_atBot (hp f) (hs f))
  simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at this
  exact this.congr (fun u => by simp [angleSum])

/-- **`Θ_i → 0` as `u → +∞`** when the star has no self-face.

[formula: iem-flatten-monotone] -/
theorem tendsto_angleSum_atTop (hp : ∀ f, (S f).Pos) (hn : ∀ f, ¬ (S f).SelfFace) :
    Tendsto (angleSum S) atTop (𝓝 0) := by
  have := tendsto_finsetSum (Finset.univ : Finset ι)
    (fun f _ => (S f).tendsto_angleAtI_atTop (hp f) (hn f))
  simp only [Finset.sum_const_zero] at this
  exact this.congr (fun u => by simp [angleSum])

/-- **The bracket of S9 contains a root**: if `Θ_i(u_lo) ≥ Θ̃ ≥ Θ_i(u_hi)` with `u_lo ≤ u_hi`,
some `u ∈ [u_lo, u_hi]` has `Θ_i(u) = Θ̃`.

[formula: iem-flatten-monotone] -/
theorem exists_root_of_bracket (hp : ∀ f, (S f).Pos) {lo hi target : ℝ} (hlh : lo ≤ hi)
    (hlo : target ≤ angleSum S lo) (hhi : angleSum S hi ≤ target) :
    ∃ u ∈ Set.Icc lo hi, angleSum S u = target :=
  intermediate_value_Icc' hlh (angleSum_continuous S hp).continuousOn ⟨hhi, hlo⟩

/-- **A root exists** for every target `0 < Θ̃ < deg·π` (interior `Θ̃ = 2π` needs `deg ≥ 3`,
boundary `Θ̃ = π` needs `deg ≥ 2`), when the star has no self-face.

[formula: iem-flatten-monotone] -/
theorem exists_root (hp : ∀ f, (S f).Pos) (hs : ∀ f, (S f).InStar) (hn : ∀ f, ¬ (S f).SelfFace)
    {target : ℝ} (h0 : 0 < target) (h1 : target < Fintype.card ι * π) :
    ∃ u, angleSum S u = target := by
  obtain ⟨lo, hlo⟩ := ((tendsto_angleSum_atBot S hp hs).eventually (lt_mem_nhds h1)).exists
  obtain ⟨hi, hhi⟩ := ((tendsto_angleSum_atTop S hp hn).eventually (gt_mem_nhds h0)).exists
  have hlh : lo ≤ hi := by
    by_contra hc
    push Not at hc
    have := angleSum_antitone S hp hc.le
    linarith
  obtain ⟨u, -, hu⟩ := exists_root_of_bracket S hp hlh hlo.le hhi.le
  exact ⟨u, hu⟩

/-- **The root is unique** when some face of the star (not a self-face) is a strict triangle at it.

[formula: iem-flatten-monotone] -/
theorem root_unique_of_valid (hp : ∀ f, (S f).Pos) {u₁ u₂ : ℝ}
    (h : angleSum S u₁ = angleSum S u₂) (g : ι) (hs : (S g).InStar) (hn : ¬ (S g).SelfFace)
    (hv : (S g).Valid u₁) : u₁ = u₂ := by
  rcases lt_trichotomy u₁ u₂ with h12 | h12 | h12
  · have := angleSum_strictAnti_of_valid S hp h12 g hs hn (Or.inl hv); linarith
  · exact h12
  · have := angleSum_strictAnti_of_valid S hp h12 g hs hn (Or.inr hv); linarith

end IntrinsicSimplification

end DDG
