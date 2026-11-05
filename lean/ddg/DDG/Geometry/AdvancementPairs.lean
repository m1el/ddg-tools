import DDG.Geometry.Advancement

/-!
# Repulsive Surfaces: conservative advancement with moving triangles (SPEC §7.3, D20)

`rs-ccd-advancement` extends the curves package's `rc-ccd-advancement` (`RC.primDist_lipschitz`,
`RC.advance_safe`, which are stated for primitives with any finite set of corners) to the pairs of
the surface CCD: a vertex against a triangle (`Fin 1` against `Fin 3` corners) and an edge against
an edge (`Fin 2` against `Fin 2`), in any real normed space (for `ℝ³`,
`EuclideanSpace ℝ (Fin 3)`). The corners move linearly, `p_i(τ) = p_i − τ u_i`; with the bound
`l = max_i |u_i| + max_j |u'_j|` on the relative motion (SPEC §7.3) and a gap `g ≥ 0` below the
current distance `d(τ_k)`, the distance stays `≥ g` on `[τ_k, τ_k + (d(τ_k) − g)/l]`: the advanced
time is a lower bound of the first time the pair comes within `g` (in particular, of the first
contact). The same holds for an obstacle triangle (corner speeds `0`).
-/

noncomputable section

namespace DDG

namespace RS

open RC (Weights primDist primDist_lipschitz advance_safe)

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- **Conservative advancement for two moving primitives with a gap** (SPEC §7.3): with corner
speed bounds `l_A`, `l_B` (`l = l_A + l_B > 0`) and `g < d(τ_k)`, the distance is at least `g` on
`[τ_k, τ_k + (d(τ_k) − g)/l]`, and the advanced time is later than `τ_k`.

[formula: rs-ccd-advancement] -/
theorem ccd_advance {ι κ : Type*} [Fintype ι] [Fintype κ] [DecidableEq ι] [DecidableEq κ]
    [Nonempty ι] [Nonempty κ] (p a : ι → E) (q b : κ → E) {lA lB : ℝ}
    (ha : ∀ i, ‖a i‖ ≤ lA) (hb : ∀ j, ‖b j‖ ≤ lB) (hl : 0 < lA + lB) {τk g : ℝ}
    (hg : g < primDist p a q b τk) :
    τk < τk + (primDist p a q b τk - g) / (lA + lB) ∧
      ∀ τ, τk ≤ τ → τ ≤ τk + (primDist p a q b τk - g) / (lA + lB) → g ≤ primDist p a q b τ := by
  have h := advance_safe (f := primDist p a q b) hl (fun τ hτ => ?_) hg
  · exact ⟨h.1, fun τ h1 h2 => (h.2 τ h1 h2).1⟩
  · have := primDist_lipschitz p a q b ha hb τk τ
    rwa [abs_of_nonneg (by linarith), mul_comm] at this

/-- **Vertex against a moving triangle** (SPEC §7.3): the bound uses `max |u_vertex| +
max |u_triangle corner|`.

[formula: rs-ccd-advancement] -/
theorem ccd_vertex_triangle (p a : Fin 1 → E) (q b : Fin 3 → E) {lA lB : ℝ}
    (ha : ∀ i, ‖a i‖ ≤ lA) (hb : ∀ j, ‖b j‖ ≤ lB) (hl : 0 < lA + lB) {τk g : ℝ}
    (hg : g < primDist p a q b τk) :
    ∀ τ, τk ≤ τ → τ ≤ τk + (primDist p a q b τk - g) / (lA + lB) → g ≤ primDist p a q b τ :=
  (ccd_advance p a q b ha hb hl hg).2

/-- **Edge against edge** (SPEC §7.3).

[formula: rs-ccd-advancement] -/
theorem ccd_edge_edge (p a : Fin 2 → E) (q b : Fin 2 → E) {lA lB : ℝ}
    (ha : ∀ i, ‖a i‖ ≤ lA) (hb : ∀ j, ‖b j‖ ≤ lB) (hl : 0 < lA + lB) {τk g : ℝ}
    (hg : g < primDist p a q b τk) :
    ∀ τ, τk ≤ τ → τ ≤ τk + (primDist p a q b τk - g) / (lA + lB) → g ≤ primDist p a q b τ :=
  (ccd_advance p a q b ha hb hl hg).2

end RS

end DDG
