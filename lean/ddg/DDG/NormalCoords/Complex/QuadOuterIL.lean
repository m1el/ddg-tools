import DDG.NormalCoords.Complex.QuadSide

/-! # `ic-flip-normal` (C part), outer side `il`: see `Complex/Quad.lean`. -/

namespace DDG.IntegerCoords

namespace Glue

variable {nA nB : Fin 3 → ℕ} {m : ℕ}

set_option maxHeartbeats 8000000 in
/-- **The re-cut arcs are the same curves, side `il`**: the curve entering the quadrilateral at
crossing `p` of `il` leaves it at the same place through the old faces and through the new ones. -/
theorem Hyp.outer_il (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) {p : ℕ} (hp : p < nB 1) :
    oldExit nA nB .il p = newExit nA nB m .il p := by
  glue_side

end Glue

end DDG.IntegerCoords
