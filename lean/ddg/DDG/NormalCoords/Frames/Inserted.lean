import DDG.NormalCoords.Frames.Angle

/-!
# Input-aligned frames at an inserted vertex (`ic-input-aligned-frame`, part 3 of 4)

`crates/ddg-integer-coords/SPEC.md` §4.12 (I26); the Rust `segment_direction` (cases
`ChordEndPlace::InFace` and `OnEdge`) and `along_direction` (case `a ∈ V*`),
`crates/ddg-integer-coords/src/frames.rs`; the frames of IT-SPEC N7, N10 (a point of face `g`:
angle zero along halfedge `3g`; a point of an edge: zero along the edge from its lower-numbered
vertex; a boundary point: zero along the boundary halfedge that leaves it, `π` along the arriving
one). An inserted vertex is flat (`Θ = 2π` inside, `π` on the boundary), so its scale is `1`
(`Intrinsic.angleScale_flat`): the stored direction is the plain angle in the layout.

`IsDir e d ψ` says that `ψ` is the direction of `d` in the frame whose zero direction is `e`:
`‖e‖ d = ‖d‖ e^{iψ} e`. It determines `ψ` in `[0, 2π)` (`isDir_unique`), and
`wrap (angleFrom e d)` is that `ψ` (`isDir_wrap`).

* `frame_in_face`: inside input face `g` the Rust `wrap(angle_from(c₁ − c₀, d))` is the direction
  of `d` in the frame of `g` (zero along side `0`, the halfedge `3g`).
* `frame_on_edge`: on an interior input edge `E`, with `x` the side of the face the segment
  leaves into (`0 ≤ x × d`: then `α = angle_from(x, d) ∈ [0, π]`) and `e` the direction of `E` from
  `edge_vertices(E)[0]` to `[1]`: `x` is a positive multiple of `e` (`β = 0`) or a negative one
  (`β = π`), and `wrap(β + α)` is the direction of `d` in the frame of `E`.
* `frame_on_boundary_edge`: on a boundary input edge, `x` the input boundary halfedge (the
  surface on its left, so `0 ≤ x × d`): `clamp(α, 0, π) = α ∈ [0, π]`, the direction of `d` in the
  boundary frame.
* `frame_along_edge` (`n = −1`): along `E` from parameter `t` to `t' ≠ t`, the Rust `0` (for
  `t' > t`) or `π` is the direction of the segment in the frame of `E`; `frame_along_boundary`:
  along the boundary halfedge that leaves the vertex, `0`.
-/

noncomputable section

namespace DDG.NormalCoords

namespace Frames

open Complex

/-- `ψ` is the direction of `d` in the frame whose zero direction is `e`. -/
def IsDir (e d : ℝ × ℝ) (ψ : ℝ) : Prop :=
  (‖toC e‖ : ℂ) * toC d = ‖toC d‖ * cexp (ψ * I) * toC e

lemma exp_sub_two_pi (ψ : ℝ) : cexp (((ψ - 2 * Real.pi : ℝ) : ℂ) * I) = cexp ((ψ : ℂ) * I) := by
  push_cast
  rw [sub_mul, Complex.exp_sub, Complex.exp_two_pi_mul_I, div_one]

/-- **`wrap (angleFrom e d)` is a direction of `d`.** -/
theorem isDir_wrap (e d : ℝ × ℝ) : IsDir e d (wrap (angleFrom e d)) := by
  unfold IsDir
  rw [wrap_exp]
  exact angleFrom_spec e d

/-- **The direction in `[0, 2π)` is unique** (`e, d ≠ 0`): it is `wrap (angleFrom e d)`. -/
theorem isDir_unique {e d : ℝ × ℝ} (he : e ≠ 0) (hd : d ≠ 0) {ψ : ℝ} (h0 : 0 ≤ ψ)
    (h1 : ψ < 2 * Real.pi) (h : IsDir e d ψ) : ψ = wrap (angleFrom e d) := by
  unfold IsDir at h
  rcases le_or_gt ψ Real.pi with hp | hp
  · rw [angleFrom_unique he hd ⟨by linarith [Real.pi_pos], hp⟩ h, wrap_eq_self h0 h1]
  · have h' : (‖toC e‖ : ℂ) * toC d =
        ‖toC d‖ * cexp (((ψ - 2 * Real.pi : ℝ) : ℂ) * I) * toC e := by
      rw [exp_sub_two_pi]; exact h
    rw [angleFrom_unique he hd ⟨by linarith, by linarith⟩ h', wrap, toIcoMod_sub]
    exact ((toIcoMod_eq_self Real.two_pi_pos).mpr ⟨h0, by linarith⟩).symm

/-- **Inside input face `g`**: the Rust `wrap(angle_from(c₁ − c₀, d))` is the direction of `d` in
the frame of the face (zero along halfedge `3g`, unscaled), the only one in `[0, 2π)`.
[formula: ic-input-aligned-frame] -/
theorem frame_in_face {c₀ c₁ d : ℝ × ℝ} (hx : c₁ - c₀ ≠ 0) (hd : d ≠ 0) :
    let φ := wrap (angleFrom (c₁ - c₀) d)
    0 ≤ φ ∧ φ < 2 * Real.pi ∧ IsDir (c₁ - c₀) d φ ∧
      ∀ ψ, 0 ≤ ψ → ψ < 2 * Real.pi → IsDir (c₁ - c₀) d ψ → ψ = φ :=
  ⟨(wrap_mem _).1, (wrap_mem _).2, isDir_wrap _ _, fun _ h0 h1 h => isDir_unique hx hd h0 h1 h⟩

lemma toC_smul (c : ℝ) (e : ℝ × ℝ) : toC (c • e) = (c : ℂ) * toC e := by
  apply Complex.ext <;> simp [toC]

/-- A direction relative to `x = c e` is one relative to `e`, turned by `0` (`c > 0`) or `π`
(`c < 0`). -/
lemma isDir_of_smul {e d : ℝ × ℝ} {c : ℝ} (hc : c ≠ 0) {α β : ℝ}
    (hβ : (0 < c ∧ β = 0) ∨ (c < 0 ∧ β = Real.pi)) (h : IsDir (c • e) d α) :
    IsDir e d (β + α) := by
  unfold IsDir at h ⊢
  rw [toC_smul, norm_mul, Complex.norm_real, Real.norm_eq_abs] at h
  have hsign : (c : ℂ) = ((|c| : ℝ) : ℂ) * cexp (β * I) := by
    rcases hβ with ⟨hc, rfl⟩ | ⟨hc, rfl⟩
    · rw [abs_of_pos hc]; simp
    · rw [abs_of_neg hc, Complex.exp_pi_mul_I]; push_cast; ring
  have hc' : ((|c| : ℝ) : ℂ) ≠ 0 := by exact_mod_cast (abs_pos.mpr hc).ne'
  apply mul_left_cancel₀ hc'
  push_cast at h ⊢
  rw [add_mul, Complex.exp_add]
  linear_combination h + (‖toC d‖ : ℂ) * cexp (α * I) * toC e * hsign

/-- **On an interior input edge `E`** (Rust `OnEdge`, interior): the segment leaves into the face
of the side `x = c e` of `E` (`0 ≤ x × d`), `e` the direction of `E` from `edge_vertices(E)[0]`,
`c > 0` when `x` runs from that vertex (`β = 0`), `c < 0` otherwise (`β = π`). Then
`α = angle_from(x, d) ∈ [0, π]` and `wrap(β + α)` is the direction of `d` in the frame of `E`
(zero along `e`), the only one in `[0, 2π)`. [formula: ic-input-aligned-frame] -/
theorem frame_on_edge {e d : ℝ × ℝ} (he : e ≠ 0) (hd : d ≠ 0) {c β : ℝ}
    (hβ : (0 < c ∧ β = 0) ∨ (c < 0 ∧ β = Real.pi)) (hin : 0 ≤ cross (c • e) d) :
    let φ := wrap (β + angleFrom (c • e) d)
    (0 ≤ angleFrom (c • e) d ∧ angleFrom (c • e) d ≤ Real.pi) ∧ 0 ≤ φ ∧ φ < 2 * Real.pi ∧
      IsDir e d φ ∧ ∀ ψ, 0 ≤ ψ → ψ < 2 * Real.pi → IsDir e d ψ → ψ = φ := by
  intro φ
  have hc : c ≠ 0 := by rcases hβ with ⟨h, -⟩ | ⟨h, -⟩ <;> [exact h.ne'; exact h.ne]
  have h1 : IsDir e d (β + angleFrom (c • e) d) := isDir_of_smul hc hβ (angleFrom_spec _ _)
  have hφ : IsDir e d φ := by
    unfold IsDir at h1 ⊢; rw [wrap_exp]; exact h1
  have hm := wrap_mem (β + angleFrom (c • e) d)
  refine ⟨angleFrom_mem_Icc hin, hm.1, hm.2, hφ, fun ψ h0 h2 hψ => ?_⟩
  rw [isDir_unique he hd h0 h2 hψ]
  exact (isDir_unique he hd hm.1 hm.2 hφ).symm

/-- **On a boundary input edge** (Rust `OnEdge`, boundary): `x` is the input boundary halfedge,
the surface on its left, so the segment has `0 ≤ x × d`; the Rust `clamp(angle_from(x, d), 0, π)`
does nothing and is the direction of `d` in the boundary frame (zero along `x`), in `[0, π]`.
[formula: ic-input-aligned-frame] -/
theorem frame_on_boundary_edge {x d : ℝ × ℝ} (hin : 0 ≤ cross x d) :
    max 0 (min (angleFrom x d) Real.pi) = angleFrom x d ∧
      0 ≤ angleFrom x d ∧ angleFrom x d ≤ Real.pi ∧ IsDir x d (angleFrom x d) := by
  obtain ⟨h0, h1⟩ := angleFrom_mem_Icc hin
  exact ⟨by rw [min_eq_left h1, max_eq_right h0], h0, h1, angleFrom_spec x d⟩

/-- **Along input edge `E`** (`n = −1`, Rust `along_direction` at an inserted vertex): from the
point at parameter `t` of `E` (`A + t e`, `e = B − A` from `edge_vertices(E)[0] = A`) to the
point at `t' ≠ t`, the Rust direction `0` if `t' > t`, `π` if `t' < t`, is the direction of the
segment in the frame of `E`. [formula: ic-input-aligned-frame] -/
theorem frame_along_edge {e : ℝ × ℝ} {t t' : ℝ} (hne : t' ≠ t) :
    IsDir e ((t' - t) • e) (if t' > t then 0 else Real.pi) := by
  unfold IsDir
  rw [toC_smul, norm_mul, Complex.norm_real, Real.norm_eq_abs]
  split_ifs with h
  · rw [abs_of_pos (sub_pos.mpr h)]; simp; ring
  · rw [abs_of_neg (sub_neg.mpr (lt_of_le_of_ne (not_lt.mp h) hne)), Complex.exp_pi_mul_I]
    push_cast; ring

/-- **Along the boundary halfedge that leaves the vertex** (the Rust `0`): the segment is a
positive multiple of the input boundary halfedge `x`. [formula: ic-input-aligned-frame] -/
theorem frame_along_boundary {x : ℝ × ℝ} {c : ℝ} (hc : 0 < c) : IsDir x (c • x) 0 := by
  unfold IsDir
  rw [toC_smul, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hc]
  simp; ring

end Frames

end DDG.NormalCoords
