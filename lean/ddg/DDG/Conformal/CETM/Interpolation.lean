import DDG.Conformal.CETM.Flatten

/-!
# Intrinsic simplification: projective interpolation of tracked points (Eq. 10)

`crates/ddg-intrinsic-simplification/SPEC.md` §3.4 (Commit), §7.1, §13; the paper Sec. 7.1,
"Vertex flattening", Eq. 10; Springborn, Schröder, Pinkall 2008, Sec. 3.4 and App. D.

A point of a face of the star with barycentric coordinates `b` (corners `0, 1, 2`) is moved by
Eq. 10 to `b̃_k = w_k b_k / Σ_m w_m b_m` with `w_k = e^{u_i}` for every corner equal to `i` and
`w_k = 1` otherwise (`eq10Weights`; in a Δ-complex the factor goes on every corner equal to `i`).
Everything below holds for **any positive per-corner weights** (`projInterp`), so it holds for
Eq. 10 as printed and for its inverse:

* the map sends the closed simplex onto itself bijectively, with inverse the map with the
  weights `1/w` (`projInterp_mapsTo`, `projInterp_inv`, `projInterp_bijOn`);
* it fixes the corners and keeps every side (`projInterp_single`, `projInterp_eq_zero_iff`);
* on a side its value depends only on the weights and coordinates of the two ends
  (`projInterp_on_side`), and the weights are a function of the vertex, so the two faces on an
  edge move the points of the edge in the same way: the map is continuous across edges (C⁰).

**A reading of the paper.** Springborn et al. 2008 (Sec. 3.4, App. D) interpolate the *image* with
the homogeneous weights `e^{−u_k}`: the point with coordinates `b` in the old face goes to the point
with coordinates `∝ e^{−u_k} b_k` in the new face, and only that map sends the circumcircle of the
old face onto the circumcircle of the new one (their Proposition). Eq. 10 (and the reference
code, `flatten_interior_vertex_and_cost.cpp:350`, `BC *= exp(u)`) uses `e^{+u_i}`: it is the
inverse of the CETM map. Both are continuous bijections that fix the corners (the theorems here),
so the tracked map stays bijective either way; only the CETM sign is "the projective interpolation
scheme of Springborn et al." (a float check: `lean/oracle`, `iem-selftest`).

[formula: iem-projective-interpolation]
-/

noncomputable section

namespace DDG

namespace IntrinsicSimplification

open Real Finset

/-- The closed triangle in barycentric coordinates. -/
def Simplex : Set (Fin 3 → ℝ) := {b | (∀ k, 0 ≤ b k) ∧ ∑ k, b k = 1}

/-- The projective re-weighting of barycentric coordinates by positive corner weights `w`. -/
def projInterp (w b : Fin 3 → ℝ) : Fin 3 → ℝ := fun k => w k * b k / ∑ m, w m * b m

/-- The weights of Eq. 10: `e^{u}` at the corners that are `i`, `1` elsewhere. -/
def eq10Weights (atI : Fin 3 → Bool) (u : ℝ) : Fin 3 → ℝ :=
  fun k => if atI k then exp u else 1

theorem eq10Weights_pos (atI : Fin 3 → Bool) (u : ℝ) (k : Fin 3) : 0 < eq10Weights atI u k := by
  unfold eq10Weights; split_ifs <;> positivity

/-- Eq. 10 as printed: `(e^{u} b_i, b_j, b_k)/(e^{u} b_i + b_j + b_k)` for a face with the
single corner `0` at `i`.

[formula: iem-projective-interpolation] -/
theorem projInterp_eq10 (u : ℝ) (b : Fin 3 → ℝ) :
    projInterp (eq10Weights ![true, false, false] u) b =
      ![exp u * b 0 / (exp u * b 0 + b 1 + b 2), b 1 / (exp u * b 0 + b 1 + b 2),
        b 2 / (exp u * b 0 + b 1 + b 2)] := by
  funext k
  fin_cases k <;> simp [projInterp, eq10Weights, Fin.sum_univ_three]

variable {w : Fin 3 → ℝ}

theorem sum_pos_of_simplex (hw : ∀ k, 0 < w k) {b : Fin 3 → ℝ} (hb : b ∈ Simplex) :
    0 < ∑ m, w m * b m := by
  obtain ⟨h0, h1⟩ := hb
  have : ∃ k, 0 < b k := by
    by_contra hc
    push Not at hc
    have : ∑ k, b k ≤ 0 := Finset.sum_nonpos (fun k _ => hc k)
    linarith
  obtain ⟨k, hk⟩ := this
  calc 0 < w k * b k := mul_pos (hw k) hk
    _ ≤ ∑ m, w m * b m :=
      Finset.single_le_sum (f := fun m => w m * b m) (fun m _ => mul_nonneg (hw m).le (h0 m))
        (Finset.mem_univ k)

/-- **The simplex is mapped into itself.**

[formula: iem-projective-interpolation] -/
theorem projInterp_mapsTo (hw : ∀ k, 0 < w k) :
    Set.MapsTo (projInterp w) (Simplex) (Simplex) := by
  intro b hb
  have hs := sum_pos_of_simplex hw hb
  refine ⟨fun k => div_nonneg (mul_nonneg (hw k).le (hb.1 k)) hs.le, ?_⟩
  simp only [projInterp]
  rw [← Finset.sum_div, div_self hs.ne']

/-- **The weights `1/w` undo the map** (on the simplex).

[formula: iem-projective-interpolation] -/
theorem projInterp_inv (hw : ∀ k, 0 < w k) {b : Fin 3 → ℝ} (hb : b ∈ Simplex) :
    projInterp (fun k => (w k)⁻¹) (projInterp w b) = b := by
  have hs := sum_pos_of_simplex hw hb
  funext k
  have hk : ∀ m, (w m)⁻¹ * (w m * b m / ∑ n, w n * b n) = b m / ∑ n, w n * b n := fun m => by
    field_simp [(hw m).ne']
  simp only [projInterp, hk]
  rw [← Finset.sum_div, hb.2, one_div]
  field_simp

/-- **Eq. 10 is a bijection of each face onto itself.**

[formula: iem-projective-interpolation] -/
theorem projInterp_bijOn (hw : ∀ k, 0 < w k) :
    Set.BijOn (projInterp w) (Simplex) (Simplex) := by
  have hw' : ∀ k, 0 < (fun k => (w k)⁻¹) k := fun k => inv_pos.mpr (hw k)
  refine Set.InvOn.bijOn ⟨fun b hb => ?_, fun b hb => ?_⟩ (projInterp_mapsTo hw)
    (projInterp_mapsTo hw')
  · exact projInterp_inv hw hb
  · have := projInterp_inv hw' hb
    simpa using this

/-- **The corners are fixed.**

[formula: iem-projective-interpolation] -/
theorem projInterp_single (hw : ∀ k, 0 < w k) (k : Fin 3) :
    projInterp w (Pi.single k 1) = Pi.single k 1 := by
  funext m
  have hs : ∑ n, w n * (Pi.single k (1 : ℝ) : Fin 3 → ℝ) n = w k := by
    simp [Pi.single_apply]
  simp only [projInterp, hs]
  by_cases h : m = k
  · subst h; simp [(hw m).ne']
  · simp [h]

/-- **Every side is kept**: a coordinate is zero after the map exactly when it was zero.

[formula: iem-projective-interpolation] -/
theorem projInterp_eq_zero_iff (hw : ∀ k, 0 < w k) {b : Fin 3 → ℝ}
    (hb : b ∈ Simplex) (k : Fin 3) : projInterp w b k = 0 ↔ b k = 0 := by
  have hs := sum_pos_of_simplex hw hb
  simp only [projInterp, div_eq_zero_iff, mul_eq_zero, (hw k).ne', false_or, hs.ne', or_false]

/-- **Continuity across edges (C⁰).** On the side opposite corner `k` (`b_k = 0`) the value at the
two other corners `p q` is `w_p b_p / (w_p b_p + w_q b_q)` and `w_q b_q / (…)`: it depends only on
the weights of the two ends of the side, which are those of the two vertices (`e^{u}` at `i`, `1`
elsewhere), so the two faces on an edge agree on it.

[formula: iem-projective-interpolation] -/
theorem projInterp_on_side {b : Fin 3 → ℝ} {k p q : Fin 3} (hpq : p ≠ q) (hkp : k ≠ p)
    (hkq : k ≠ q) (hk : b k = 0) :
    projInterp w b p = w p * b p / (w p * b p + w q * b q) ∧
      projInterp w b q = w q * b q / (w p * b p + w q * b q) := by
  have hs : ∑ m, w m * b m = w p * b p + w q * b q := by
    have huniv : (Finset.univ : Finset (Fin 3)) = {k, p, q} := by
      ext m; fin_cases k <;> fin_cases p <;> fin_cases q <;> fin_cases m <;> simp_all
    rw [huniv, Finset.sum_insert (by simp [hkp, hkq]), Finset.sum_pair hpq, hk, mul_zero, zero_add]
  simp only [projInterp, hs, and_self]

end IntrinsicSimplification

end DDG
