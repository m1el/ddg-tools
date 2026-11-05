import DdgProofs.Intrinsic.Diamond
import DDG.Geometry.Barycentric

/-!
# Intrinsic triangulations: one step of a tracing query

`crates/ddg-intrinsic/SPEC.md` §4 (Tracing), N7; the paper Sec. 3.2.2.

The triangles along a path are unfolded into one plane in which the path is the positive
x-axis. A triangle is entered through an edge with one end above the axis (`up`, `y > 0`) and
one below (`down`, `y < 0`); its third vertex `apex` is unfolded to the left of `up → down`
(`unfoldApex`, `DdgProofs/Intrinsic/Diamond.lean`), and the sign of its `y` decides through which
of the two other edges the path leaves.

* **`it-ray-edge-intersection`** — `axisCrossingT y_tail y_head` is the parameter
  at which the line of the edge meets the axis, and the only one; the point of the edge at that
  parameter is `(axisCrossingX …, 0)` (`axisCrossing_on_axis`, `axisCrossing_unique`). If the
  ends are not strictly on the same side (`y_tail ≤ 0 ≤ y_head` or the reverse, not both zero)
  the parameter lies in `[0, 1]`, so the point is on the edge; it is `0` or `1` exactly when
  that end lies on the axis (`axisCrossingT_mem_Icc`, `axisCrossingT_mem_Ioo`). The parameter
  is the level-crossing parameter of `DDG/Geometry/Barycentric.lean`
  (`axisCrossingT_eq_crossingParam`, `cut-edge-crossing`). The exit point lies ahead of the entry point on the path
  (`exit_ahead_above`, `exit_ahead_below`).
* **`it-unfold-apex`** (`apexAboveAxis`) — a closed segment meets the axis exactly
  when its ends are not strictly on the same side (`meetsAxis_iff`). Hence, for an entry edge
  with `up` above and `down` below: the edge `down–apex` meets the axis iff `y_apex ≥ 0`, the
  edge `apex–up` iff `y_apex ≤ 0`; one of them always does, and both do exactly when the apex
  lies on the axis, the corner case (`exit_cases`).

Definitions: the closed forms of the Rust code on `ℝ` (`axisCrossingT`, `axisCrossingX`,
`apexAboveAxis`) and `MeetsAxis` (a closed segment has a point on the x-axis).
-/

noncomputable section

namespace DDG

namespace Intrinsic

open Real

/-! ## The crossing of an edge with the axis -/

/-- Parameter at which the segment from `(x_t, y_t)` to `(x_h, y_h)` crosses the x-axis (the ray
of a tracing step): `t = y_tail/(y_tail − y_head)`. -/
def axisCrossingT (yTail yHead : ℝ) : ℝ := yTail / (yTail - yHead)

/-- Abscissa of that crossing: `x = x_tail + t (x_head − x_tail)`. -/
def axisCrossingX (xTail xHead t : ℝ) : ℝ := xTail + t * (xHead - xTail)

/-- The path leaves through the edge from the new apex to the vertex on the other side of the
axis: `true` iff the apex has `y > 0`. -/
def apexAboveAxis (yApex : ℝ) : Bool := decide (0 < yApex)

theorem axisCrossingT_eq (yT yH : ℝ) : axisCrossingT yT yH = yT / (yT - yH) := rfl

theorem axisCrossingX_eq (xT xH t : ℝ) : axisCrossingX xT xH t = xT + t * (xH - xT) := rfl

/-- The parameter of a tracing step is the level-crossing parameter of the cut flow
(`cut-edge-crossing`), with the heights over the path in place of the level-set values. -/
theorem axisCrossingT_eq_crossingParam (yT yH : ℝ) :
    axisCrossingT yT yH = levelCrossingParam yT yH := rfl

/-- **The stated parameters are the intersection of the axis with the line of the edge**: with
`t = axisCrossingT y_tail y_head` (`y_tail ≠ y_head`), the point `tail + t (head − tail)` is
`(axisCrossingX x_tail x_head t, 0)`.

[formula: it-ray-edge-intersection] -/
theorem axisCrossing_on_axis {tail head : ℝ × ℝ} (h : tail.2 ≠ head.2) :
    (tail.1 + axisCrossingT tail.2 head.2 * (head.1 - tail.1),
        tail.2 + axisCrossingT tail.2 head.2 * (head.2 - tail.2)) =
      (axisCrossingX tail.1 head.1 (axisCrossingT tail.2 head.2), 0) := by
  have h' : tail.2 - head.2 ≠ 0 := sub_ne_zero.mpr h
  rw [axisCrossingX_eq, axisCrossingT_eq]
  ext
  · rfl
  · dsimp only
    field_simp
    ring

/-- **… and the only one**: the line of the edge meets the axis at no other parameter.

[formula: it-ray-edge-intersection] -/
theorem axisCrossing_unique {tail head : ℝ × ℝ} (h : tail.2 ≠ head.2) (s : ℝ) :
    tail.2 + s * (head.2 - tail.2) = 0 ↔ s = axisCrossingT tail.2 head.2 := by
  have h' : tail.2 - head.2 ≠ 0 := sub_ne_zero.mpr h
  rw [axisCrossingT_eq, eq_div_iff h']
  constructor <;> intro h1 <;> linarith

/-- **The point is on the edge**: if the tail is not above and the head not below the axis
(not both on it), the parameter lies in `[0, 1]`; it is `0` exactly when the tail is on the axis
and `1` exactly when the head is.

[formula: it-ray-edge-intersection] -/
theorem axisCrossingT_mem_Icc {yT yH : ℝ} (hT : yT ≤ 0) (hH : 0 ≤ yH) (hne : yT < yH) :
    0 ≤ axisCrossingT yT yH ∧ axisCrossingT yT yH ≤ 1 ∧
      (axisCrossingT yT yH = 0 ↔ yT = 0) ∧ (axisCrossingT yT yH = 1 ↔ yH = 0) := by
  have hd : yT - yH < 0 := by linarith
  rw [axisCrossingT_eq]
  refine ⟨div_nonneg_of_nonpos hT hd.le, ?_, ?_, ?_⟩
  · rw [div_le_one_of_neg hd]
    linarith
  · rw [div_eq_zero_iff]
    exact ⟨fun h => h.resolve_right hd.ne, Or.inl⟩
  · rw [div_eq_one_iff_eq hd.ne]
    constructor <;> intro h <;> linarith

/-- The same with the tail above and the head below. -/
theorem axisCrossingT_mem_Icc' {yT yH : ℝ} (hT : 0 ≤ yT) (hH : yH ≤ 0) (hne : yH < yT) :
    0 ≤ axisCrossingT yT yH ∧ axisCrossingT yT yH ≤ 1 ∧
      (axisCrossingT yT yH = 0 ↔ yT = 0) ∧ (axisCrossingT yT yH = 1 ↔ yH = 0) := by
  have hd : 0 < yT - yH := by linarith
  rw [axisCrossingT_eq]
  refine ⟨div_nonneg hT hd.le, ?_, ?_, ?_⟩
  · rw [div_le_one hd]
    linarith
  · rw [div_eq_zero_iff]
    exact ⟨fun h => h.resolve_right hd.ne', Or.inl⟩
  · rw [div_eq_one_iff_eq hd.ne']
    constructor <;> intro h <;> linarith

/-- With the ends strictly on opposite sides the crossing is interior to the edge. -/
theorem axisCrossingT_mem_Ioo {yT yH : ℝ} (hT : yT < 0) (hH : 0 < yH) :
    0 < axisCrossingT yT yH ∧ axisCrossingT yT yH < 1 := by
  obtain ⟨h0, h1, e0, e1⟩ := axisCrossingT_mem_Icc hT.le hH.le (by linarith)
  refine ⟨lt_of_le_of_ne h0 fun h => ?_, lt_of_le_of_ne h1 fun h => ?_⟩
  · exact hT.ne (e0.mp h.symm)
  · exact hH.ne' (e1.mp h)

/-! ## Which edge the path leaves through -/

theorem apexAboveAxis_iff (y : ℝ) : apexAboveAxis y = true ↔ 0 < y := by
  simp only [apexAboveAxis, decide_eq_true_eq]

/-- The closed segment from `p` to `q` has a point on the x-axis. -/
def MeetsAxis (p q : ℝ × ℝ) : Prop := ∃ t : ℝ, 0 ≤ t ∧ t ≤ 1 ∧ p.2 + t * (q.2 - p.2) = 0

/-- A segment meets the axis exactly when its ends are not strictly on the same side. -/
theorem meetsAxis_iff (p q : ℝ × ℝ) :
    MeetsAxis p q ↔ (p.2 ≤ 0 ∧ 0 ≤ q.2) ∨ (q.2 ≤ 0 ∧ 0 ≤ p.2) := by
  constructor
  · rintro ⟨t, h0, h1, h⟩
    by_contra hcon
    push Not at hcon
    obtain ⟨c1, c2⟩ := hcon
    rcases lt_trichotomy p.2 0 with hp | hp | hp
    · have hq := c1 hp.le
      have a1 : t * q.2 ≤ 0 := mul_nonpos_of_nonneg_of_nonpos h0 hq.le
      have a2 : (1 - t) * p.2 ≤ 0 := mul_nonpos_of_nonneg_of_nonpos (sub_nonneg.mpr h1) hp.le
      have a3 : t * q.2 = 0 := by linarith
      have ht : t = 0 := (mul_eq_zero.mp a3).resolve_right hq.ne
      rw [ht] at h
      linarith
    · have hq := c1 hp.le
      have := c2 hq.le
      linarith
    · have hq : 0 < q.2 := by
        by_contra hq'
        linarith [c2 (not_lt.mp hq')]
      have a1 : 0 ≤ t * q.2 := mul_nonneg h0 hq.le
      have a2 : 0 ≤ (1 - t) * p.2 := mul_nonneg (sub_nonneg.mpr h1) hp.le
      have a3 : t * q.2 = 0 := by linarith
      have ht : t = 0 := (mul_eq_zero.mp a3).resolve_right hq.ne'
      rw [ht] at h
      linarith
  · rintro (⟨hp, hq⟩ | ⟨hq, hp⟩)
    · rcases eq_or_lt_of_le (hp.trans hq) with he | hlt
      · exact ⟨0, le_rfl, zero_le_one, by rw [zero_mul, add_zero]; linarith⟩
      · obtain ⟨h0, h1, -, -⟩ := axisCrossingT_mem_Icc hp hq hlt
        refine ⟨axisCrossingT p.2 q.2, h0, h1, ?_⟩
        exact (axisCrossing_unique hlt.ne _).mpr rfl
    · rcases eq_or_lt_of_le (hq.trans hp) with he | hlt
      · exact ⟨0, le_rfl, zero_le_one, by rw [zero_mul, add_zero]; linarith⟩
      · obtain ⟨h0, h1, -, -⟩ := axisCrossingT_mem_Icc' hp hq hlt
        refine ⟨axisCrossingT p.2 q.2, h0, h1, ?_⟩
        exact (axisCrossing_unique hlt.ne' _).mpr rfl

/-- **The conditions that pick the exit edge are exclusive except at a corner.** For an entry
edge with `up` above and `down` below the axis and the third vertex `apex`: the edge
`down–apex` meets the axis iff `y_apex ≥ 0` (`apexAboveAxis`, or the apex on the
axis); the edge `apex–up` meets it iff `y_apex ≤ 0`; one of them always does; both do exactly
when the apex lies on the axis.

[formula: it-unfold-apex] -/
theorem exit_cases {up down : ℝ × ℝ} (hu : 0 < up.2) (hd : down.2 < 0) (apex : ℝ × ℝ) :
    (MeetsAxis down apex ↔ 0 ≤ apex.2) ∧ (MeetsAxis apex up ↔ apex.2 ≤ 0) ∧
      (MeetsAxis down apex ∨ MeetsAxis apex up) ∧
      ((MeetsAxis down apex ∧ MeetsAxis apex up) ↔ apex.2 = 0) ∧
      (apexAboveAxis apex.2 = true → MeetsAxis down apex ∧ ¬ MeetsAxis apex up) := by
  have e1 : MeetsAxis down apex ↔ 0 ≤ apex.2 := by
    rw [meetsAxis_iff]
    constructor
    · rintro (⟨-, h⟩ | ⟨-, h⟩)
      · exact h
      · linarith
    · exact fun h => Or.inl ⟨hd.le, h⟩
  have e2 : MeetsAxis apex up ↔ apex.2 ≤ 0 := by
    rw [meetsAxis_iff]
    constructor
    · rintro (⟨h, -⟩ | ⟨h, -⟩)
      · exact h
      · linarith
    · exact fun h => Or.inl ⟨h, hu.le⟩
  refine ⟨e1, e2, ?_, ?_, ?_⟩
  · rw [e1, e2]
    exact le_total 0 apex.2
  · rw [e1, e2]
    exact ⟨fun h => le_antisymm h.2 h.1, fun h => ⟨h.ge, h.le⟩⟩
  · intro h
    rw [apexAboveAxis_iff] at h
    rw [e1, e2]
    exact ⟨h.le, not_le.mpr h⟩

/-! ## The path moves forward -/

/-- **Leaving through `down–apex`** (the apex above the axis, to the left of `up → down`): the
exit abscissa exceeds the entry abscissa by `t · [up, down, apex] / (y_up − y_down)`, with `t`
the parameter of the exit crossing and `[·]` the signed double area. So the exit point lies
strictly ahead.

[formula: it-ray-edge-intersection] -/
theorem exit_ahead_above {up down apex : ℝ × ℝ} (hu : 0 < up.2) (hd : down.2 < 0)
    (ha : 0 < apex.2) (hleft : 0 < sda up down apex) :
    axisCrossingX down.1 apex.1 (axisCrossingT down.2 apex.2) -
        axisCrossingX up.1 down.1 (axisCrossingT up.2 down.2) =
      axisCrossingT down.2 apex.2 * sda up down apex / (up.2 - down.2) ∧
    axisCrossingX up.1 down.1 (axisCrossingT up.2 down.2) <
      axisCrossingX down.1 apex.1 (axisCrossingT down.2 apex.2) := by
  have h1 : up.2 - down.2 ≠ 0 := by linarith
  have h2 : down.2 - apex.2 ≠ 0 := by linarith
  have e : axisCrossingX down.1 apex.1 (axisCrossingT down.2 apex.2) -
        axisCrossingX up.1 down.1 (axisCrossingT up.2 down.2) =
      axisCrossingT down.2 apex.2 * sda up down apex / (up.2 - down.2) := by
    simp only [axisCrossingX_eq, axisCrossingT_eq]
    unfold sda
    field_simp
    ring
  refine ⟨e, ?_⟩
  have ht := (axisCrossingT_mem_Ioo hd ha).1
  have : 0 < axisCrossingT down.2 apex.2 * sda up down apex / (up.2 - down.2) :=
    div_pos (mul_pos ht hleft) (by linarith)
  linarith

/-- **Leaving through `apex–up`** (the apex below the axis): the exit abscissa exceeds the
entry abscissa by `(1 − t) · [up, down, apex] / (y_up − y_down)`.

[formula: it-ray-edge-intersection] -/
theorem exit_ahead_below {up down apex : ℝ × ℝ} (hu : 0 < up.2) (hd : down.2 < 0)
    (ha : apex.2 < 0) (hleft : 0 < sda up down apex) :
    axisCrossingX apex.1 up.1 (axisCrossingT apex.2 up.2) -
        axisCrossingX up.1 down.1 (axisCrossingT up.2 down.2) =
      (1 - axisCrossingT apex.2 up.2) * sda up down apex / (up.2 - down.2) ∧
    axisCrossingX up.1 down.1 (axisCrossingT up.2 down.2) <
      axisCrossingX apex.1 up.1 (axisCrossingT apex.2 up.2) := by
  have h1 : up.2 - down.2 ≠ 0 := by linarith
  have h2 : apex.2 - up.2 ≠ 0 := by linarith
  have e : axisCrossingX apex.1 up.1 (axisCrossingT apex.2 up.2) -
        axisCrossingX up.1 down.1 (axisCrossingT up.2 down.2) =
      (1 - axisCrossingT apex.2 up.2) * sda up down apex / (up.2 - down.2) := by
    simp only [axisCrossingX_eq, axisCrossingT_eq]
    unfold sda
    field_simp
    ring
  refine ⟨e, ?_⟩
  have ht := (axisCrossingT_mem_Ioo ha hu).2
  have : 0 < (1 - axisCrossingT apex.2 up.2) * sda up down apex / (up.2 - down.2) :=
    div_pos (mul_pos (by linarith) hleft) (by linarith)
  linarith

end Intrinsic

end DDG
