import DDG.Geometry.Vec3

/-!
# Repulsive Curves: dual lengths

"Repulsive Curves", Yu, Schumacher, Crane, ACM TOG 40(2), 2021; contract
`crates/ddg-repulsive-curves/SPEC.md` §2.1 (D23).

* `RC.dualLen` — the dual length (lumped mass) of a vertex, `m_v = ½ Σ_{I ∋ v} ℓ_I`;
* `RC.sum_dualLen` — the dual lengths sum to the total length (`rc-dual-length`);
* `RC.dualLen_nonneg`.

(Split off `DDG.Geometry.Vec3`; the declarations keep their names.)
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Real Filter Topology

section edges

variable {V E : Type*}

/-! ## Dual lengths -/

variable [Fintype V] [Fintype E] [DecidableEq V]

/-- The dual length (lumped mass) of a vertex, `m_v = ½ Σ_{I ∋ v} ℓ_I` (SPEC §2.1, D23); an edge
is counted once per endpoint at `v`. -/
def dualLen (γ : V → V3) (ep : E → Fin 2 → V) (v : V) : ℝ :=
  (1 / 2) * ∑ I, ∑ a : Fin 2, if ep I a = v then len γ ep I else 0

/-- **The dual lengths sum to the total length**: `Σ_v m_v = Σ_I ℓ_I`.

[formula: rc-dual-length] -/
theorem sum_dualLen (γ : V → V3) (ep : E → Fin 2 → V) :
    ∑ v, dualLen γ ep v = ∑ I, len γ ep I := by
  have h : ∀ I : E, ∑ v : V, ∑ a : Fin 2, (if ep I a = v then len γ ep I else 0) =
      2 * len γ ep I := by
    intro I
    rw [Finset.sum_comm]
    simp only [Finset.sum_ite_eq, Finset.mem_univ, ite_true, Fin.sum_univ_two]
    ring
  calc ∑ v, dualLen γ ep v
      = (1 / 2) * ∑ I, ∑ v : V, ∑ a : Fin 2, (if ep I a = v then len γ ep I else 0) := by
        simp only [dualLen, ← Finset.mul_sum]
        rw [Finset.sum_comm]
    _ = ∑ I, len γ ep I := by
        simp only [h, Finset.mul_sum]
        refine Finset.sum_congr rfl fun I _ => by ring

omit [Fintype V] in
theorem dualLen_nonneg (γ : V → V3) (ep : E → Fin 2 → V) (v : V) : 0 ≤ dualLen γ ep v := by
  unfold dualLen
  refine mul_nonneg (by norm_num) (Finset.sum_nonneg fun I _ => Finset.sum_nonneg fun a _ => ?_)
  split_ifs
  · exact len_nonneg γ ep I
  · exact le_rfl

end edges

end RC

end DDG
