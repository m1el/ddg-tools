import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring

/-!
# Integer coordinates: intersecting the unfolded strip (`ic-strip-intersection`, tier A)

`crates/ddg-integer-coords/SPEC.md` §3.1 (Alg. 2 lines 6–17, the geometry of `ExtractCurve` and
`ExtractIntrinsicEdge`), §3.2 (I5); the Rust `Chord::frame` and `Chord::hit`
(`crates/ddg-integer-coords/src/geometry.rs`).

The faces crossed by a curve are unfolded into the plane; the curve (or the intrinsic edge) is the
straight segment from `p` to `q = p + d`, `len = ‖d‖ > 0`. Each interior edge of the strip is
crossed from the segment's right to its left: its tail `x_t` lies below the segment's line and its
head `x_h` above. `Chord::frame` writes a point in the frame in which the segment is the positive
x-axis from the origin: `frame x = (((x − p)·d)/len, (d × (x − p))/len)`. With `a = frame x_t`,
`b = frame x_h`, the edge parameter is `t = a_y / (a_y − b_y)` (`axis_crossing`,
`it-ray-edge-intersection`) and the segment parameter `s = (a_x + t (b_x − a_x)) / len`.

* `hit_point`: the point `x_t + t (x_h − x_t)` of the edge **is** the point `p + s d` of the
  segment: the two computed parameters name one point of the plane, where the lines meet.
* `hit_unique`: it is the only point where the edge's line meets the segment's line, and
  `0 < t < 1` (the point is inside the edge) when the tail is strictly below and the head strictly
  above (`hit_mem`); so the segment's line meets each interior edge of the strip once.
* `hit_reverse`: followed from the other end (the segment from `q` to `p`, crossing the twin
  halfedge `x_h → x_t`, as the backward trace of Alg. 2 and the reversed edge see it), the
  parameters are `1 − t` and `1 − s`: both directions give the same crossing.

That `0 ≤ s ≤ 1` (the crossing lies on the segment, not only on its line) holds when the strip is the
face sequence of the geodesic with exact lengths; with rounding the code clamps and reports
(I5). This is not part of the statement.
-/

noncomputable section

namespace DDG.Intrinsic

namespace Strip

/-- `Chord::frame`: the point `x` in the frame of the segment from `p` along `d` (`len = ‖d‖`). -/
def frame (p d : ℝ × ℝ) (len : ℝ) (x : ℝ × ℝ) : ℝ × ℝ :=
  (((x.1 - p.1) * d.1 + (x.2 - p.2) * d.2) / len, (d.1 * (x.2 - p.2) - d.2 * (x.1 - p.1)) / len)

/-- The edge parameter of `Chord::hit` (`axis_crossing`). -/
def hitT (a b : ℝ × ℝ) : ℝ := a.2 / (a.2 - b.2)

/-- The segment parameter of `Chord::hit`. -/
def hitS (len : ℝ) (a b : ℝ × ℝ) : ℝ := (a.1 + hitT a b * (b.1 - a.1)) / len

variable {p d xt xh : ℝ × ℝ} {len : ℝ}

/-- **The two parameters name one point**: the point of the edge at `t` is the point of the
segment at `s`. [formula: ic-strip-intersection] -/
theorem hit_point (hlen : 0 < len) (hd : len ^ 2 = d.1 ^ 2 + d.2 ^ 2)
    (hne : (frame p d len xt).2 ≠ (frame p d len xh).2) :
    let t := hitT (frame p d len xt) (frame p d len xh)
    let s := hitS len (frame p d len xt) (frame p d len xh)
    xt.1 + t * (xh.1 - xt.1) = p.1 + s * d.1 ∧ xt.2 + t * (xh.2 - xt.2) = p.2 + s * d.2 := by
  intro t s
  have hl : len ≠ 0 := hlen.ne'
  set a := frame p d len xt with ha
  set b := frame p d len xh with hb
  have hden : a.2 - b.2 ≠ 0 := sub_ne_zero.mpr hne
  have ha1 : a.1 * len = (xt.1 - p.1) * d.1 + (xt.2 - p.2) * d.2 := by
    simp only [ha, frame]; field_simp
  have hb1 : b.1 * len = (xh.1 - p.1) * d.1 + (xh.2 - p.2) * d.2 := by
    simp only [hb, frame]; field_simp
  have ha2 : a.2 * len = d.1 * (xt.2 - p.2) - d.2 * (xt.1 - p.1) := by
    simp only [ha, frame]; field_simp
  have hb2 : b.2 * len = d.1 * (xh.2 - p.2) - d.2 * (xh.1 - p.1) := by
    simp only [hb, frame]; field_simp
  have ht : t * (a.2 - b.2) = a.2 := by simp only [t, hitT, ← ha, ← hb]; field_simp
  have hs : s * len = a.1 + t * (b.1 - a.1) := by simp only [s, t, hitS, ← ha, ← hb]; field_simp
  clear_value t s a b
  have z : a.2 + t * (b.2 - a.2) = 0 := by linear_combination (-1 : ℝ) * ht
  set X1 := xt.1 + t * (xh.1 - xt.1) with hX1
  set X2 := xt.2 + t * (xh.2 - xt.2) with hX2
  have k1 : s * len * len = (X1 - p.1) * d.1 + (X2 - p.2) * d.2 := by
    rw [hs, hX1, hX2]; linear_combination (1 - t) * ha1 + t * hb1
  have k2 : d.1 * (X2 - p.2) - d.2 * (X1 - p.1) = 0 := by
    rw [hX1, hX2]; linear_combination len * z - (1 - t) * ha2 - t * hb2
  have hlen2 : len * len ≠ 0 := mul_ne_zero hl hl
  constructor
  · have e : (X1 - p.1 - s * d.1) * (len * len) = 0 := by
      linear_combination (X1 - p.1) * hd - d.1 * k1 - d.2 * k2
    rcases mul_eq_zero.mp e with e | e
    · linarith
    · exact absurd e hlen2
  · have e : (X2 - p.2 - s * d.2) * (len * len) = 0 := by
      linear_combination (X2 - p.2) * hd - d.2 * k1 + d.1 * k2
    rcases mul_eq_zero.mp e with e | e
    · linarith
    · exact absurd e hlen2

/-- **The only meeting point**: a point of the edge's line at parameter `u` is on the segment's
line (frame height `0`) only for `u = t`. [formula: ic-strip-intersection] -/
theorem hit_unique (hne : (frame p d len xt).2 ≠ (frame p d len xh).2) (u : ℝ) :
    (frame p d len xt).2 + u * ((frame p d len xh).2 - (frame p d len xt).2) = 0 ↔
      u = hitT (frame p d len xt) (frame p d len xh) := by
  have hden : (frame p d len xt).2 - (frame p d len xh).2 ≠ 0 := sub_ne_zero.mpr hne
  simp only [hitT]
  rw [eq_div_iff hden]
  constructor <;> intro h <;> linarith

/-- **Inside the edge**: with the tail strictly below the segment's line and the head strictly
above, `0 < t < 1`. [formula: ic-strip-intersection] -/
theorem hit_mem {a b : ℝ × ℝ} (ha : a.2 < 0) (hb : 0 < b.2) : 0 < hitT a b ∧ hitT a b < 1 := by
  have hd : a.2 - b.2 < 0 := by linarith
  simp only [hitT]
  constructor
  · exact div_pos_of_neg_of_neg ha hd
  · rw [div_lt_one_of_neg hd]; linarith

/-- **Both directions give the same crossing**: the segment from `q = p + d` back to `p` (direction
`−d`), crossing the twin halfedge `x_h → x_t`, has the parameters `1 − t` (from the twin's tail)
and `1 − s` (from `q`). [formula: ic-strip-intersection] -/
theorem hit_reverse (hlen : 0 < len) (hd : len ^ 2 = d.1 ^ 2 + d.2 ^ 2)
    (hne : (frame p d len xt).2 ≠ (frame p d len xh).2) :
    let q : ℝ × ℝ := (p.1 + d.1, p.2 + d.2)
    let d' : ℝ × ℝ := (-d.1, -d.2)
    hitT (frame q d' len xh) (frame q d' len xt) = 1 - hitT (frame p d len xt) (frame p d len xh) ∧
      hitS len (frame q d' len xh) (frame q d' len xt) =
        1 - hitS len (frame p d len xt) (frame p d len xh) := by
  intro q d'
  have hl : len ≠ 0 := hlen.ne'
  have hden : (frame p d len xt).2 - (frame p d len xh).2 ≠ 0 := sub_ne_zero.mpr hne
  have h1 : (frame q d' len xh).2 = -(frame p d len xh).2 := by
    simp only [frame, q, d']; field_simp; ring
  have h2 : (frame q d' len xt).2 = -(frame p d len xt).2 := by
    simp only [frame, q, d']; field_simp; ring
  have h3 : (frame q d' len xh).1 = len - (frame p d len xh).1 := by
    simp only [frame, q, d']; field_simp; rw [hd]; ring
  have h4 : (frame q d' len xt).1 = len - (frame p d len xt).1 := by
    simp only [frame, q, d']; field_simp; rw [hd]; ring
  have ht : hitT (frame q d' len xh) (frame q d' len xt) =
      1 - hitT (frame p d len xt) (frame p d len xh) := by
    simp only [hitT, h1, h2]
    have : -(frame p d len xh).2 - -(frame p d len xt).2 ≠ 0 := by
      intro h; apply hden; linarith
    field_simp
    ring
  refine ⟨ht, ?_⟩
  simp only [hitS, ht, h3, h4]
  simp only [hitT]
  field_simp
  ring

end Strip

end DDG.Intrinsic
