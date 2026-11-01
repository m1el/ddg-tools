import DdgProofs.IntegerCoords.Complex.QuadSide

/-! # `ic-flip-normal` (C part), outer side `lj`: see `Complex/Quad.lean`. -/

namespace DdgProofs.IntegerCoords

namespace Glue

variable {nA nB : Fin 3 → ℕ} {m : ℕ}

set_option maxHeartbeats 8000000 in
/-- **The re-cut arcs are the same curves, side `lj`**: the curve entering the quadrilateral at
crossing `p` of `lj` leaves it at the same place through the old faces and through the new ones. -/
theorem Hyp.outer_lj (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) {p : ℕ} (hp : p < nB 2) :
    oldExit nA nB .lj p = newExit nA nB m .lj p := by
  glue_side

end Glue

end DdgProofs.IntegerCoords
