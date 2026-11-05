import DdgProofs.IntegerCoords.Complex.QuadSide

/-! # `ic-flip-normal` (C part), outer side `jk`: see `Complex/Quad.lean`. -/

namespace DDG.IntegerCoords

namespace Glue

variable {nA nB : Fin 3 → ℕ} {m : ℕ}

set_option maxHeartbeats 8000000 in
/-- **The re-cut arcs are the same curves, side `jk`**: the curve entering the quadrilateral at
crossing `p` of `jk` leaves it at the same place through the old faces and through the new ones. -/
theorem Hyp.outer_jk (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) {p : ℕ} (hp : p < nA 1) :
    oldExit nA nB .jk p = newExit nA nB m .jk p := by
  glue_side

end Glue

end DDG.IntegerCoords
