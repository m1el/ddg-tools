import DdgProofs.IntegerCoords.Complex.Quad

/-! # The new face `F₁` of the flip is valid (`Complex/Quad.lean`); a file of its own so that the
two case analyses build in parallel. -/

namespace DdgProofs.IntegerCoords

open Ddg.IntegerCoords (faceValid)

namespace Glue

variable {nA nB : Fin 3 → ℕ} {m : ℕ}

set_option maxHeartbeats 4000000 in
/-- **The new faces are valid**: the corner counts of the re-cut are admissible for the new side
counts (when there is at most one arc from `k` to `l`, i.e. `n_kl ≥ −1`). -/
theorem Hyp.admissible1 (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) : Admissible (n1 nA nB m) (c1 nA nB) := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  refine ⟨fun k => ?_, fun k => ?_, fun k => ?_⟩ <;>
    rcases fin3_cases k with rfl | rfl | rfl <;>
    simp only [em, n1, c1, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] <;>
    glue_facts h <;> (clear dD; omega)

lemma Hyp.cc1 (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) : cc (n1 nA nB m) = c1 nA nB := by
  have ha := h.admissible1 hD
  have hv : faceValid (n1 nA nB m 0) (n1 nA nB m 1) (n1 nA nB m 2) = true :=
    (FaceModel.exists_valid_iff _).1 ⟨_, rfl, canon_valid ha⟩
  funext k
  have e1 := ha.c_eq k
  have e2 := (admissible_cc hv).c_eq k
  exact_mod_cast e2.trans e1.symm

lemma Hyp.model1 (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) :
    (model (n1 nA nB m)).partner = canon (n1 nA nB m) (c1 nA nB) := by
  simp only [model, h.cc1 hD]

end Glue

end DdgProofs.IntegerCoords
