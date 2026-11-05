import DDG.NormalCoords.Complex.QuadSide

/-! # `ic-flip-normal` (C part), outer side `ki`: see `Complex/Quad.lean`. -/

namespace DDG.IntegerCoords

namespace Glue

variable {nA nB : Fin 3 → ℕ} {m : ℕ}

set_option maxHeartbeats 8000000 in
/-- **The re-cut arcs are the same curves, side `ki`**: the curve entering the quadrilateral at
crossing `p` of `ki` leaves it at the same place through the old faces and through the new ones. -/
theorem Hyp.outer_ki (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) {p : ℕ} (hp : p < nA 2) :
    oldExit nA nB .ki p = newExit nA nB m .ki p := by
  glue_side

end Glue

end DDG.IntegerCoords
