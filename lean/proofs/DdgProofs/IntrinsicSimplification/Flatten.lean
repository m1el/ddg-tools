import DdgProofs.Yamabe

/-!
# Intrinsic simplification: flattening one vertex (paper Sec. 4.1, Eq. 3, 4)

`crates/ddg-intrinsic-simplification/SPEC.md` §3.4, S9, §13; the paper Sec. 4.1; Springborn,
Schröder, Pinkall 2008 ("CETM"), Sec. 3.1 (extended angles), Eq. 9, 10.

**The star.** The vertex `i` being flattened has the star faces `S f` (`f : ι`, any finite index
type; in a Δ-complex a face may have one, two or three corners at `i`). A face is a `StarFace`:
its three sides `l0 l1 l2` (side `k` is opposite corner `k`, i.e. it joins the corners `k + 1`
and `k + 2`) and three flags `c0 c1 c2` saying which corners are `i`. Flattening scales every
length by Eq. 3 with `u_j = 0` for `j ≠ i`: side `k` gets `e^{u m_k / 2}`, `m_k` the number of its
two ends that are `i` (`0` for an outer edge, `1` for a spoke, `2` for a self-edge `ii`). The
angle sum of `i` is `Θ_i(u) = Σ_f Σ_{corners k at i} θ_k(u)` (`angleSum`).

**Angles.** `cornerAngle a b c = arccos((b² + c² − a²)/(2bc))` with Mathlib's `Real.arccos`, which
is constant `π` below `−1` and `0` above `1`. For positive lengths this *is* the extended angle of
CETM Sec. 3.1: a face violating the triangle inequality gets the angles `π, 0, 0`
(`cornerAngle_eq_pi_of_ge`, `cornerAngle_eq_zero_of_ge`), and the three angles of any face with
positive sides sum to `π` (`cornerAngle_sum_pos`).

* **`iem-flatten-newton`** — `hasDerivAt_angleSum`: where every face is a strict triangle,
  `dΘ_i/du = −½ Σ_f hessTerm f`, with `hessTerm` the sum, over the corners of `f` at `i`, of the
  cotangents of the angles opposite its sides **that are not self-edges**. For a face with one
  corner at `i` that is the paper's `cot θ^k_ij + cot θ^j_ki` (`hessTerm_one_corner`); so Eq. 4 is
  Newton's step for `Θ_i(u) = Θ̃_i` (`newtonStep_eq`, `newtonStep_eq_paper`). For a face with two
  corners at `i` (a self-edge `ii` in the star) the per-corner sum of Eq. 4 (SPEC §3.4 "each face
  once per corner at `i`") is **not** the derivative: the cotangent opposite the self-edge must be
  dropped (`hessTerm_two_corners`; it only changes the speed of Newton's method, not the root).
  The derivative is negative (`hessTerm_pos`) when no face is a self-face.
* **`iem-flatten-monotone`** — `angleSum_continuous`, `angleSum_antitone`: `Θ_i` is continuous and
  non-increasing on all of `ℝ`; `tendsto_angleSum_atBot`: `→ (number of faces)·π` as `u → −∞`;
  `tendsto_angleSum_atTop`: `→ 0` as `u → +∞` (no self-face); `angleSum_strictAnti_of_valid`:
  strictly decreasing between two values of `u` at one of which some face (not a self-face) is a
  strict triangle; `exists_root_of_bracket` (the solver's bracket of S9 contains a root),
  `exists_root` (a root exists for every target in `(0, deg·π)`), `root_unique_of_valid`.
-/

noncomputable section

namespace DdgProofs

namespace IntrinsicSimplification

open Real Filter Topology

/-! ## Extended angles -/

/-- A side at least the sum of the other two: the angle opposite it is `π`. -/
theorem cornerAngle_eq_pi_of_ge {a b c : ℝ} (hb : 0 < b) (hc : 0 < c) (h : b + c ≤ a) :
    cornerAngle a b c = π := by
  unfold cornerAngle
  rw [arccos_eq_pi, div_le_iff₀ (by positivity)]
  nlinarith

/-- A side at most the difference of the other two: the angle opposite it is `0`. -/
theorem cornerAngle_eq_zero_of_ge {a b c : ℝ} (hb : 0 < b) (hc : 0 < c) (ha : 0 ≤ a)
    (h : a + c ≤ b ∨ a + b ≤ c) : cornerAngle a b c = 0 := by
  unfold cornerAngle
  rw [arccos_eq_zero, le_div_iff₀ (by positivity)]
  rcases h with h | h <;> nlinarith

/-- **Extended angle sum.** For any positive sides (a strict triangle or not) the three angles sum
to `π`; a degenerate face has the angles `π, 0, 0` (CETM Sec. 3.1).

[formula: iem-flatten-monotone] -/
theorem cornerAngle_sum_pos {a b c : ℝ} (ha : 0 < a) (hb : 0 < b) (hc : 0 < c) :
    cornerAngle a b c + cornerAngle b c a + cornerAngle c a b = π := by
  by_cases h : IsTriangle a b c
  · exact cornerAngle_add_cornerAngle_add_cornerAngle h
  have h' : b + c ≤ a ∨ c + a ≤ b ∨ a + b ≤ c := by
    by_contra hn
    push Not at hn
    exact h ⟨hn.1, hn.2.1, hn.2.2⟩
  rcases h' with h1 | h1 | h1
  · rw [cornerAngle_eq_pi_of_ge hb hc h1, cornerAngle_eq_zero_of_ge hc ha hb.le (Or.inr (by linarith)),
      cornerAngle_eq_zero_of_ge ha hb hc.le (Or.inl (by linarith))]
    ring
  · rw [cornerAngle_eq_zero_of_ge hb hc ha.le (Or.inl (by linarith)), cornerAngle_eq_pi_of_ge hc ha h1,
      cornerAngle_eq_zero_of_ge ha hb hc.le (Or.inr (by linarith))]
    ring
  · rw [cornerAngle_eq_zero_of_ge hb hc ha.le (Or.inr (by linarith)),
      cornerAngle_eq_zero_of_ge hc ha hb.le (Or.inl (by linarith)), cornerAngle_eq_pi_of_ge ha hb h1]
    ring

/-! ## The star of the vertex being flattened -/

/-- `1` if the flag is set, else `0`. -/
def ind (b : Bool) : ℝ := if b then 1 else 0

@[simp] theorem ind_true : ind true = 1 := rfl
@[simp] theorem ind_false : ind false = 0 := rfl

/-- A face of the star of the vertex `i`: the sides `l0 l1 l2` (side `k` opposite corner `k`) and
the flags `c0 c1 c2` (corner `k` is `i`). -/
structure StarFace where
  /-- side opposite corner 0 (joins corners 1 and 2) -/
  l0 : ℝ
  /-- side opposite corner 1 (joins corners 2 and 0) -/
  l1 : ℝ
  /-- side opposite corner 2 (joins corners 0 and 1) -/
  l2 : ℝ
  /-- corner 0 is `i` -/
  c0 : Bool
  /-- corner 1 is `i` -/
  c1 : Bool
  /-- corner 2 is `i` -/
  c2 : Bool

namespace StarFace

variable (f : StarFace)

/-- The sides after the flattening scale `u` (Eq. 3 with `u_j = 0` for `j ≠ i`): side `k` times
`e^{u m_k/2}`, `m_k` the number of its ends at `i`. -/
def L0 (u : ℝ) : ℝ := f.l0 * exp (u * (ind f.c1 + ind f.c2) / 2)
/-- see `L0` -/
def L1 (u : ℝ) : ℝ := f.l1 * exp (u * (ind f.c2 + ind f.c0) / 2)
/-- see `L0` -/
def L2 (u : ℝ) : ℝ := f.l2 * exp (u * (ind f.c0 + ind f.c1) / 2)

/-- The angles at the three corners after scaling (extended angles). -/
def θ0 (u : ℝ) : ℝ := cornerAngle (f.L0 u) (f.L1 u) (f.L2 u)
/-- see `θ0` -/
def θ1 (u : ℝ) : ℝ := cornerAngle (f.L1 u) (f.L2 u) (f.L0 u)
/-- see `θ0` -/
def θ2 (u : ℝ) : ℝ := cornerAngle (f.L2 u) (f.L0 u) (f.L1 u)

/-- Cotangents of the three angles after scaling. -/
def cot0 (u : ℝ) : ℝ := cornerCot (f.L0 u) (f.L1 u) (f.L2 u)
/-- see `cot0` -/
def cot1 (u : ℝ) : ℝ := cornerCot (f.L1 u) (f.L2 u) (f.L0 u)
/-- see `cot0` -/
def cot2 (u : ℝ) : ℝ := cornerCot (f.L2 u) (f.L0 u) (f.L1 u)

/-- The face's share of the angle sum of `i`: the angles at its corners that are `i`. -/
def angleAtI (u : ℝ) : ℝ := ind f.c0 * f.θ0 u + ind f.c1 * f.θ1 u + ind f.c2 * f.θ2 u

/-- The face's term of the denominator of Eq. 4: for each corner at `i`, the cotangents of the
angles opposite its two sides, leaving out a side that is a self-edge (both ends at `i`). Corner 0
has the sides 2 (to corner 1, opposite corner 2) and 1 (to corner 2, opposite corner 1). -/
def hessTerm (u : ℝ) : ℝ :=
  ind f.c0 * ((if f.c1 then 0 else f.cot2 u) + (if f.c2 then 0 else f.cot1 u)) +
  ind f.c1 * ((if f.c2 then 0 else f.cot0 u) + (if f.c0 then 0 else f.cot2 u)) +
  ind f.c2 * ((if f.c0 then 0 else f.cot1 u) + (if f.c1 then 0 else f.cot0 u))

/-- The scaled face is a strict triangle. -/
def Valid (u : ℝ) : Prop := IsTriangle (f.L0 u) (f.L1 u) (f.L2 u)

/-- Positive sides. -/
def Pos : Prop := 0 < f.l0 ∧ 0 < f.l1 ∧ 0 < f.l2

/-- A self-face: all three corners are `i`. -/
def SelfFace : Prop := f.c0 = true ∧ f.c1 = true ∧ f.c2 = true

/-- The face belongs to the star: some corner is `i`. -/
def InStar : Prop := f.c0 = true ∨ f.c1 = true ∨ f.c2 = true

/-! ### `iem-flatten-newton`: the derivative of the angle sum -/

theorem L0_add (u t : ℝ) :
    f.L0 (u + t) = f.L0 u * exp (t * ((ind f.c1 + ind f.c2) / 2)) := by
  unfold L0; rw [mul_assoc, ← exp_add]; congr 2; ring
theorem L1_add (u t : ℝ) :
    f.L1 (u + t) = f.L1 u * exp (t * ((ind f.c2 + ind f.c0) / 2)) := by
  unfold L1; rw [mul_assoc, ← exp_add]; congr 2; ring
theorem L2_add (u t : ℝ) :
    f.L2 (u + t) = f.L2 u * exp (t * ((ind f.c0 + ind f.c1) / 2)) := by
  unfold L2; rw [mul_assoc, ← exp_add]; congr 2; ring

/-- From the derivative at `0` of `t ↦ g (u + t)` to the derivative of `g` at `u`. -/
theorem hasDerivAt_of_shift {g : ℝ → ℝ} {d u : ℝ} (h : HasDerivAt (fun t => g (u + t)) d 0) :
    HasDerivAt g d u := by
  have h' : HasDerivAt (fun t => g (u + t)) d (u - u) := by rwa [sub_self]
  have := h'.comp_sub_const u u
  simpa using this

theorem hasDerivAt_θ0 {u : ℝ} (h : f.Valid u) :
    HasDerivAt f.θ0 (1 / 2 * f.cot2 u * (ind f.c1 - ind f.c0) +
      1 / 2 * f.cot1 u * (ind f.c2 - ind f.c0)) u := by
  apply hasDerivAt_of_shift
  have := hasDerivAt_cornerAngle_conformal h (ind f.c0) (ind f.c1) (ind f.c2)
  simp only [θ0, cot1, cot2, L0_add, L1_add, L2_add]
  exact this

theorem hasDerivAt_θ1 {u : ℝ} (h : f.Valid u) :
    HasDerivAt f.θ1 (1 / 2 * f.cot0 u * (ind f.c2 - ind f.c1) +
      1 / 2 * f.cot2 u * (ind f.c0 - ind f.c1)) u := by
  apply hasDerivAt_of_shift
  have := hasDerivAt_cornerAngle_conformal h.rotate (ind f.c1) (ind f.c2) (ind f.c0)
  simp only [θ1, cot0, cot2, L0_add, L1_add, L2_add]
  exact this

theorem hasDerivAt_θ2 {u : ℝ} (h : f.Valid u) :
    HasDerivAt f.θ2 (1 / 2 * f.cot1 u * (ind f.c0 - ind f.c2) +
      1 / 2 * f.cot0 u * (ind f.c1 - ind f.c2)) u := by
  apply hasDerivAt_of_shift
  have := hasDerivAt_cornerAngle_conformal h.rotate.rotate (ind f.c2) (ind f.c0) (ind f.c1)
  simp only [θ2, cot0, cot1, L0_add, L1_add, L2_add]
  exact this

/-- **The derivative of one face's share of `Θ_i`** is `−½ hessTerm`. -/
theorem hasDerivAt_angleAtI {u : ℝ} (h : f.Valid u) :
    HasDerivAt f.angleAtI (-(1 / 2) * f.hessTerm u) u := by
  have h0 := (f.hasDerivAt_θ0 h).const_mul (ind f.c0)
  have h1 := (f.hasDerivAt_θ1 h).const_mul (ind f.c1)
  have h2 := (f.hasDerivAt_θ2 h).const_mul (ind f.c2)
  have := (h0.add h1).add h2
  unfold angleAtI
  convert this using 1
  unfold hessTerm
  cases f.c0 <;> cases f.c1 <;> cases f.c2 <;> simp <;> ring

/-- For a face with exactly one corner at `i` (corner 0, the simplicial case) the term is the
paper's `cot θ^k_ij + cot θ^j_ki`: the cotangents of the two other corners.

[formula: iem-flatten-newton] -/
theorem hessTerm_one_corner (h : f.c0 = true ∧ f.c1 = false ∧ f.c2 = false) (u : ℝ) :
    f.hessTerm u = f.cot1 u + f.cot2 u := by
  obtain ⟨h0, h1, h2⟩ := h
  simp [hessTerm, h0, h1, h2]; ring

/-- For a face with two corners at `i` (corners 0 and 1; side 2 is a self-edge `ii`) the term is
`cot θ_0 + cot θ_1`, the cotangents of the two corners **at `i`**. The per-corner sum of Eq. 4,
`(cot θ_1 + cot θ_2) + (cot θ_0 + cot θ_2)`, has the extra `2 cot θ_2`.

[formula: iem-flatten-newton] -/
theorem hessTerm_two_corners (h : f.c0 = true ∧ f.c1 = true ∧ f.c2 = false) (u : ℝ) :
    f.hessTerm u = f.cot0 u + f.cot1 u := by
  obtain ⟨h0, h1, h2⟩ := h
  simp [hessTerm, h0, h1, h2]; ring

/-- A self-face contributes nothing (its angles at `i` always sum to `π`). -/
theorem hessTerm_self_face (h : f.SelfFace) (u : ℝ) : f.hessTerm u = 0 := by
  obtain ⟨h0, h1, h2⟩ := h
  simp [hessTerm, h0, h1, h2]

/-- Two cotangents of a strict triangle have a positive sum: `cot A + cot B = 2c²/(4·Area)`. -/
theorem cornerCot_add_pos {a b c : ℝ} (h : IsTriangle a b c) :
    0 < cornerCot a b c + cornerCot b c a := by
  have hT := triArea_pos h
  have hc := h.pos_c
  have e : triArea b c a = triArea a b c := triArea_rotate a b c
  unfold cornerCot
  rw [e, ← add_div]
  apply div_pos _ (by positivity)
  nlinarith

/-- **The Newton denominator is positive**: a strict triangle in the star that is not a self-face
has a positive term (so `Θ_i` is strictly decreasing where all faces are valid).

[formula: iem-flatten-newton] -/
theorem hessTerm_pos {u : ℝ} (h : f.Valid u) (hs : f.InStar) (hn : ¬ f.SelfFace) :
    0 < f.hessTerm u := by
  have p01 := cornerCot_add_pos h
  have p12 := cornerCot_add_pos h.rotate
  have p20 := cornerCot_add_pos h.rotate.rotate
  unfold hessTerm cot0 cot1 cot2
  unfold InStar at hs
  unfold SelfFace at hn
  revert hs hn
  cases f.c0 <;> cases f.c1 <;> cases f.c2 <;> simp <;> linarith

end StarFace

/-! ### The star -/

variable {ι : Type*} [Fintype ι]

/-- The angle sum `Θ_i(u)` of the vertex being flattened. -/
def angleSum (S : ι → StarFace) (u : ℝ) : ℝ := ∑ f, (S f).angleAtI u

/-- The denominator of Eq. 4 (twice the second derivative of CETM's energy). -/
def hessian (S : ι → StarFace) (u : ℝ) : ℝ := 1 / 2 * ∑ f, (S f).hessTerm u

/-- **`iem-flatten-newton`.** Where every face of the star is a strict triangle,
`dΘ_i/du = −½ Σ_f hessTerm f` (Springborn et al. 2008 Eq. 9, 10 restricted to one vertex).

[formula: iem-flatten-newton] -/
theorem hasDerivAt_angleSum (S : ι → StarFace) {u : ℝ} (h : ∀ f, (S f).Valid u) :
    HasDerivAt (angleSum S) (-hessian S u) u := by
  have := HasDerivAt.fun_sum (u := Finset.univ) (fun f _ => (S f).hasDerivAt_angleAtI (h f))
  unfold hessian
  rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
  convert this using 1
  · funext x; simp [angleSum]
  · congr 1; funext f; ring

/-- Eq. 4: `u ← u − (Θ̃ − Θ(u)) / (½ Σ hessTerm)`. -/
def newtonStep (S : ι → StarFace) (target u : ℝ) : ℝ :=
  u - (target - angleSum S u) / hessian S u

/-- **Eq. 4 is Newton's step** for `Θ_i(u) = Θ̃_i`: `u − g(u)/g'(u)` with `g = Θ_i − Θ̃_i`.

[formula: iem-flatten-newton] -/
theorem newtonStep_eq (S : ι → StarFace) (target : ℝ) {u : ℝ} (h : ∀ f, (S f).Valid u) :
    HasDerivAt (fun v => angleSum S v - target) (-hessian S u) u ∧
      newtonStep S target u = u - (angleSum S u - target) / (-hessian S u) := by
  refine ⟨(hasDerivAt_angleSum S h).sub_const target, ?_⟩
  unfold newtonStep
  rw [div_neg]
  ring

/-- **The paper's Eq. 4** for a simplicial star (every face has the single corner 0 at `i`):
the denominator is `½ Σ (cot θ^k_ij + cot θ^j_ki)`, positive where the faces are valid.

[formula: iem-flatten-newton] -/
theorem newtonStep_eq_paper (S : ι → StarFace) (target u : ℝ)
    (h1 : ∀ f, (S f).c0 = true ∧ (S f).c1 = false ∧ (S f).c2 = false) :
    newtonStep S target u =
      u - (target - angleSum S u) / (1 / 2 * ∑ f, ((S f).cot1 u + (S f).cot2 u)) := by
  unfold newtonStep hessian
  simp only [fun f => (S f).hessTerm_one_corner (h1 f) u]

/-- The denominator of Eq. 4 is positive when every face is valid, belongs to the star and none is
a self-face, and the star is not empty.

[formula: iem-flatten-newton] -/
theorem hessian_pos (S : ι → StarFace) [Nonempty ι] {u : ℝ} (h : ∀ f, (S f).Valid u)
    (hs : ∀ f, (S f).InStar) (hn : ∀ f, ¬ (S f).SelfFace) : 0 < hessian S u := by
  unfold hessian
  have := Finset.sum_pos (s := Finset.univ) (fun f _ => (S f).hessTerm_pos (h f) (hs f) (hn f))
    Finset.univ_nonempty
  positivity

end IntrinsicSimplification

end DdgProofs
