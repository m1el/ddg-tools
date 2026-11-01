import DdgProofs.IntegerCoords.Complex.QuadCount
import DdgProofs.IntegerCoords.Complex.QuadOuterJK
import DdgProofs.IntegerCoords.Complex.QuadOuterKI
import DdgProofs.IntegerCoords.Complex.QuadOuterIL
import DdgProofs.IntegerCoords.Complex.QuadOuterLJ

/-! # `ic-flip-normal` (the C part): the four outer sides together (`Complex/Quad.lean`). -/

namespace DdgProofs.IntegerCoords

namespace Glue

variable {nA nB : Fin 3 → ℕ} {m : ℕ}

/-- **The re-cut arcs are the same curves** (`ic-flip-normal`, the C part): on every outer side, the
curve entering the quadrilateral at any crossing leaves it at the same boundary point through the
old faces (across `ij`) and through the new faces (across `kl`). [formula: ic-flip-normal] -/
theorem Hyp.outer (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) (s : QS) {p : ℕ}
    (hp : p < qcount nA nB s) : oldExit nA nB s p = newExit nA nB m s p := by
  cases s
  · exact h.outer_jk hD hp
  · exact h.outer_ki hD hp
  · exact h.outer_il hD hp
  · exact h.outer_lj hD hp

end Glue

end DdgProofs.IntegerCoords
